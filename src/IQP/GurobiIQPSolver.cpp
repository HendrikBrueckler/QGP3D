#ifdef QGP3D_WITH_GUROBI
#ifndef QGP3D_WITHOUT_IQP

#include "QGP3D/IQP/GurobiIQPSolver.hpp"

namespace qgp3d
{
namespace impl
{

GurobiIQPSolver::GurobiIQPSolver(
    TetMeshProps& meshProps, double scaling, double varLowerBound, double maxSeconds, double individualArcFactor)
    : TetMeshNavigator(meshProps), TetMeshManipulator(meshProps), MCMeshNavigator(meshProps),
      MCMeshManipulator(meshProps), BaseIQPSolver(scaling, varLowerBound, maxSeconds, individualArcFactor),
      _env(true), _model(
        (_env.set(GRB_IntParam_LogToConsole, true),_env.start(), _env))
{
}

void GurobiIQPSolver::setupDefaultOptions()
{
#ifdef NDEBUG
    _model.set(GRB_IntParam_LogToConsole, true);
#else
    _model.set(GRB_IntParam_LogToConsole, true);
#endif
    _model.set(GRB_IntParam_Threads, 1);
    _model.set(GRB_DoubleParam_TimeLimit, _maxSeconds);
    _model.set(GRB_IntParam_MIPFocus, 3);
    _model.set(GRB_IntParam_Cuts, 2);
}

void GurobiIQPSolver::setupVariables()
{
    auto& mc = mcMeshProps().mesh();

    _arc2var.clear();
    // One integer length variable per arc
    for (EH arc : mc.edges())
    {
        double length = mcMeshProps().isAllocated<ARC_INT_LENGTH>() ? mcMeshProps().get<ARC_INT_LENGTH>(arc)
                                                                    : mcMeshProps().get<ARC_DBL_LENGTH>(arc);
        DLOG(INFO) << "Arc " << arc << " continuous length: " << length;

        // Configure gurobi var
        GRBVar x = _model.addVar(
            _varLowerBound, GRB_INFINITY, 0., GRB_INTEGER, std::string("Arc ") + std::to_string(arc.idx()));
        x.set(GRB_DoubleAttr_Start, length);
        _arc2var[arc] = {x};
    }
}

void GurobiIQPSolver::setupObjective()
{
    auto& mc = mcMeshProps().mesh();

    // Objective Function minimizes deviation from scaled seamless param
    GRBQuadExpr objective = 0.;
    if (_individualArcFactor < 1.0)
    {
        for (CH block : mc.cells())
        {
            // For each axis U/V/W
            for (UVWDir axis : {UVWDir::NEG_V_NEG_W, UVWDir::NEG_U_NEG_W, UVWDir::NEG_U_NEG_V})
            {
                GRBLinExpr varSum = 0.;
                double lenSum = 0;
                // WARNING this way of using range based for loop (chaining 2 function calls)
                // depends on the first calls not returning rvalues!
                for (EH arc : mcMeshProps().ref<BLOCK_EDGE_ARCS>(block).at(axis))
                {
                    varSum += _arc2var.at(arc);
                    lenSum += mcMeshProps().get<ARC_DBL_LENGTH>(arc);
                }
                objective += (varSum - _scaling * lenSum) * (varSum - _scaling * lenSum) * (1.0 - _individualArcFactor);
            }
        }
    }
    for (EH arc : mc.edges())
    {
        GRBVar var = _arc2var.at(arc);
        double len = mcMeshProps().get<ARC_DBL_LENGTH>(arc);
        objective += (var - _scaling * len) * (var - _scaling * len) * _individualArcFactor;
    }

    _model.setObjective(objective, GRB_MINIMIZE);
}

void GurobiIQPSolver::setupConstraints()
{
    auto& mc = mcMeshProps().mesh();

    vector<CriticalLink> criticalLinks;
    map<EH, int> a2criticalLinkIdx;
    map<VH, vector<int>> n2criticalLinksOut;
    map<VH, vector<int>> n2criticalLinksIn;
    getCriticalLinks(criticalLinks, a2criticalLinkIdx, n2criticalLinksOut, n2criticalLinksIn, true);

    // Each patches opposite arc lengths must match
    for (FH patch : mc.faces())
    {
        HFH hp = mc.halfface_handle(patch, 0);
        if (mc.is_boundary(hp))
            hp = mc.opposite_halfface_handle(hp);

        auto dir2orderedHas = halfpatchHalfarcsByDir(hp);
        assert(dir2orderedHas.size() == 4);
        UVWDir dirHpNormal = halfpatchNormalDir(hp);

        for (UVWDir dirSide : decompose(~(dirHpNormal | -dirHpNormal), {UVWDir::POS_U, UVWDir::POS_V, UVWDir::POS_W}))
        {
            UVWDir dirSideOpp = -dirSide;

            GRBLinExpr side = 0;
            for (HEH ha : dir2orderedHas[dirSide])
                side += _arc2var.at(mc.edge_handle(ha));

            GRBLinExpr sideOpp = 0;
            for (HEH ha : dir2orderedHas[-dirSide])
                sideOpp += _arc2var.at(mc.edge_handle(ha));

            std::string constraintName = "Patch" + std::to_string(patch.idx()) + "dir"
                                         + std::to_string(static_cast<uint8_t>(dirSide | dirSideOpp));
            _model.addConstr(side, GRB_EQUAL, sideOpp, constraintName);
        }
    }

    if (_varLowerBound < 0.0)
    {
        // Enforce no block with negative extension along any axis
        for (CH b : mc.cells())
        {
            for (UVWDir dir : {UVWDir::POS_V_POS_W, UVWDir::POS_U_POS_W, UVWDir::POS_U_POS_V})
            {
                GRBLinExpr sumExpr = 0;
                auto& arcs = mcMeshProps().ref<BLOCK_EDGE_ARCS>(b).at(dir);
                for (EH a : arcs)
                    sumExpr += _arc2var.at(a);
                std::string constraintName
                    = "Block" + std::to_string(b.idx()) + "dir" + std::to_string(static_cast<uint8_t>(~(dir | -dir)));
                _model.addConstr(sumExpr, GRB_GREATER_EQUAL, 0.0, constraintName);
            }
        }
    }

    // Enforce critical links of length >= 1
    for (auto& criticalLink : criticalLinks)
    {
        if (criticalLink.pathHas.empty())
            continue;
        GRBLinExpr sumExpr = 0;
        for (HEH ha : criticalLink.pathHas)
            sumExpr += _arc2var.at(mc.edge_handle(ha));

        std::string constraintName = "Separation" + std::to_string(_nSeparationConstraints);
        if (criticalLink.nFrom == criticalLink.nTo)
            _model.addConstr(sumExpr, GRB_GREATER_EQUAL, 3.0, constraintName);
        else
            _model.addConstr(sumExpr, GRB_GREATER_EQUAL, 1.0, constraintName);
        _nSeparationConstraints++;
    }

    // Additional code to avoid selfadjacency
    {
        for (CH b : mc.cells())
        {
            auto& dir2ps = mcMeshProps().ref<BLOCK_FACE_PATCHES>(b);
            for (auto& kv : dir2ps)
            {
                if (isNeg(kv.first))
                    continue;

                UVWDir dir = kv.first;
                auto& ps = kv.second;
                auto& psOpp = dir2ps.at(-dir);
                set<CH> bs, bsOpp;
                for (FH p : ps)
                {
                    HFH hp = mc.halfface_handle(p, 0);
                    if (mc.incident_cell(hp) == b)
                        hp = mc.opposite_halfface_handle(hp);
                    CH bNext = mc.incident_cell(hp);
                    if (!bNext.is_valid())
                        continue;
                    bs.insert(bNext);
                }
                for (FH p : psOpp)
                {
                    HFH hp = mc.halfface_handle(p, 0);
                    if (mc.incident_cell(hp) == b)
                        hp = mc.opposite_halfface_handle(hp);
                    CH bNext = mc.incident_cell(hp);
                    if (!bNext.is_valid())
                        continue;
                    bsOpp.insert(bNext);
                }
                bool mustBeNonZero = bs.count(b) != 0 || bsOpp.count(b) != 0 || containsSomeOf(bs, bsOpp);
                if (mustBeNonZero)
                {
                    // Constrain non-zero length along dir
                    GRBLinExpr sumExpr = 0;
                    auto& dir2as = mcMeshProps().ref<BLOCK_EDGE_ARCS>(b);
                    for (EH a : dir2as.at(decompose(~(dir | -dir), DIM_2_DIRS)[0]))
                        sumExpr += _arc2var.at(a);
                    std::string constraintName = "Selfadjacency" + std::to_string(_nSeparationConstraints);
                    _model.addConstr(sumExpr, GRB_GREATER_EQUAL, 1.0, constraintName);
                    _nSeparationConstraints++;

                    continue;
                }
            }
        }
    }
}

void GurobiIQPSolver::enableQuickSolve()
{
    _model.set(GRB_DoubleParam_MIPGap, 0.99);
}

void GurobiIQPSolver::enableExactSolve()
{
    _model.set(GRB_DoubleParam_MIPGap, 1e-4);
}

void GurobiIQPSolver::finalizeSetup()
{
    // Nothing necessary here
}

BaseIQPSolver::RetCode GurobiIQPSolver::solve()
{
    auto& mc = mcMeshProps().mesh();
    try
    {
        _model.optimize();
        auto obj = _model.getObjective();
        DLOG(INFO) << "Gurobi solved IQP with final objective value of " << obj.getValue();

        int status = _model.get(GRB_IntAttr_Status);
        if (status != 2 && status != 9 && status != 10)
        {
            LOG(ERROR) << "Bad status return by GUROBI solver";
            return BaseIQPSolver::SOLVER_ERROR;
        }

        for (EH arc : mc.edges())
        {
            auto& grbvar = _arc2var.at(arc);
            mcMeshProps().set<ARC_INT_LENGTH>(arc, (int)std::round(grbvar.get(GRB_DoubleAttr_X)));
        }
    }
    catch (GRBException e)
    {
        LOG(ERROR) << "Gurobi exception, errcode: " << e.getErrorCode();
        LOG(ERROR) << "Gurobi error message: " << e.getMessage();
        return SOLVER_ERROR;
    }

    return BaseIQPSolver::SUCCESS;
}

double GurobiIQPSolver::objectiveValue() const
{
    return _model.getObjective().getValue();
}

void GurobiIQPSolver::useCurrentAsWarmStart()
{
    for (auto& kv : _arc2var)
        kv.second.set(GRB_DoubleAttr_Start, mcMeshProps().get<ARC_INT_LENGTH>(kv.first));
}

void GurobiIQPSolver::addConstraints(const vector<vector<pair<int, EH>>>& nonZeroSum)
{
    for (auto& aColl : nonZeroSum)
    {
        vector<GRBVar> vars;
        vector<double> coeff;
        for (auto sign2a : aColl)
        {
            vars.push_back(_arc2var.at(sign2a.second));
            coeff.push_back(sign2a.first);
        }
        GRBLinExpr sumExpr = 0;
        sumExpr.addTerms(&coeff[0], &vars[0], vars.size());

        std::string constraintName = "Separation" + std::to_string(_nSeparationConstraints);
        _model.addConstr(sumExpr, GRB_GREATER_EQUAL, 1.0, constraintName);
        _nSeparationConstraints++;
    }
}

} // namespace impl
} // namespace qgp3d

#endif
#endif
