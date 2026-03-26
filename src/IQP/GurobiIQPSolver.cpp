#ifdef QGP3D_WITH_GUROBI
#ifndef QGP3D_WITHOUT_IQP

#include "QGP3D/IQP/GurobiIQPSolver.hpp"
#include "QGP3D/ObjectiveBuilder.hpp"

namespace qgp3d
{
namespace impl
{

GurobiIQPSolver::GurobiIQPSolver(
    TetMeshProps& meshProps, double varLowerBound, double maxSeconds)
    : TetMeshNavigator(meshProps), TetMeshManipulator(meshProps), MCMeshNavigator(meshProps),
      MCMeshManipulator(meshProps), BaseIQPSolver(varLowerBound, maxSeconds), _env(true),
      _model((_env.set(GRB_IntParam_LogToConsole, true), _env.start(), _env))
{
}

void GurobiIQPSolver::setupDefaultOptions()
{
#ifdef NDEBUG
    _model.set(GRB_IntParam_LogToConsole, false);
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

void GurobiIQPSolver::setupObjective(QuadraticObjective& obj)
{
    auto& mc = mcMeshProps().mesh();

    std::vector<EH> idx2arc;
    for (EH a : mc.edges())
        idx2arc.push_back(a);

    // Objective expression
    GRBQuadExpr objective = obj.f0;

    for (int i = 0; i < (int)idx2arc.size(); i++)
    {
        EH a = idx2arc[i];
        objective += obj.grad0(i) * _arc2var.at(a);
    }

    for (int k = 0; k < obj.hess.outerSize(); ++k) {
      for (Eigen::SparseMatrix<double>::InnerIterator it(obj.hess, k); it; ++it) {
        int i = it.row();
        int j = it.col();
        EH ai = idx2arc[i];
        if (i == j) {
            auto& var = _arc2var.at(ai);
            objective += (0.5 * it.value()) * var * var;
        } else if (j > i) {
            EH aj = idx2arc[j];
            objective += it.value() * _arc2var.at(ai) * _arc2var.at(aj);
        }
      }
    }

    _model.setObjective(objective, GRB_MINIMIZE);
}

void GurobiIQPSolver::setupConstraints()
{
    auto& mc = mcMeshProps().mesh();
    const bool SINGULARITIES_FIXED = !meshProps().isAllocated<ALGO_VARIANT>() || (meshProps().get<ALGO_VARIANT>() % 2);

    vector<bool> isCriticalNode, isCriticalArc, isCriticalPatch;
    vector<CriticalEntity> criticalEntities;
    map<EH, int> a2criticalLinkIdx;
    map<FH, int> p2criticalRegionIdx;
    map<VH, vector<int>> n2criticalLinksOut;
    map<VH, vector<int>> n2criticalLinksIn;
    getCriticalEntities(isCriticalNode,
                        isCriticalArc,
                        isCriticalPatch,
                        criticalEntities,
                        a2criticalLinkIdx,
                        p2criticalRegionIdx,
                        n2criticalLinksOut,
                        n2criticalLinksIn,
                        true,
                        SINGULARITIES_FIXED,
                        true);

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
    for (auto& criticalEntity : criticalEntities)
    {
        if (criticalEntity.pathHas.empty())
            continue;
        GRBLinExpr sumExpr = 0;
        for (HEH ha : criticalEntity.pathHas)
            sumExpr += _arc2var.at(mc.edge_handle(ha));

        std::string constraintName = "Separation" + std::to_string(_nSeparationConstraints);
        _model.addConstr(sumExpr, GRB_GREATER_EQUAL, 1.0, constraintName);
        _nSeparationConstraints++;
    }

    for (auto& criticalEntity : criticalEntities)
    {
        if (criticalEntity.regionPs.empty())
            continue;
        set<EH> as;
        for (FH p : criticalEntity.regionPs)
            for (EH a : mc.face_edges(p))
                as.insert(a);
        GRBLinExpr sumExpr = 0;
        for (EH a : as)
            sumExpr += _arc2var.at(a);

        std::string constraintName = "Separation" + std::to_string(_nSeparationConstraints);
        _model.addConstr(sumExpr, GRB_GREATER_EQUAL, 1.0, constraintName);
        _nSeparationConstraints++;
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
        double rhs = 1.0;
        for (auto sign2a : aColl)
        {
            if (sign2a.second.is_valid())
            {
                vars.push_back(_arc2var.at(sign2a.second));
                coeff.push_back(sign2a.first);
            }
            else
                rhs = sign2a.first;
        }
        GRBLinExpr sumExpr = 0;
        sumExpr.addTerms(&coeff[0], &vars[0], vars.size());

        std::string constraintName = "Separation" + std::to_string(_nSeparationConstraints);
        _model.addConstr(sumExpr, GRB_GREATER_EQUAL, rhs, constraintName);
        _nSeparationConstraints++;
    }
}

} // namespace impl
} // namespace qgp3d

#endif
#endif
