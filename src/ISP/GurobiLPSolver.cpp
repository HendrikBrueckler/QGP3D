#ifdef QGP3D_WITH_GUROBI

#include "QGP3D/ISP/GurobiLPSolver.hpp"

namespace qgp3d
{
namespace impl
{

GurobiLPSolver::GurobiLPSolver(const TetMeshProps& meshProps, double scaling, const ISPQuantizer::Decomposition& decomp)
    : TetMeshNavigator(meshProps), MCMeshNavigator(meshProps), BaseLPSolver(meshProps, scaling, decomp), _lpenv(true)
{
    _lpenv.set(GRB_IntParam_LogToConsole, false);
    _lpenv.set(GRB_IntParam_Threads, 1);
    _lpenv.start();
}

void GurobiLPSolver::setupLPBase()
{
    auto& mcMesh = mcMeshProps().mesh();
    _subproblem2bundle2vars.resize(_decomp.subproblem2bundles.size());
    for (int subproblem = 0; subproblem < (int)_decomp.subproblem2bundles.size(); subproblem++)
    {
        DLOG(INFO) << "Creating model for subproblem " << subproblem;
        std::unique_ptr<GRBModel> mptr = std::make_unique<GRBModel>(_lpenv);
        GRBModel& model = *mptr;
        model.set(GRB_IntParam_Method, 0); // Primal simplex

        // Variables
        for (int bundle : _decomp.subproblem2bundles[subproblem])
        {
            // Configure gurobi var
            GRBVar x
                = model.addVar(0.0, GRB_INFINITY, 0., GRB_CONTINUOUS, std::string("Bundle ") + std::to_string(bundle));
            x.set(GRB_DoubleAttr_Start, 0.0);
            GRBVar y = model.addVar(
                0.0, GRB_INFINITY, 0., GRB_CONTINUOUS, std::string("BundleMinus ") + std::to_string(bundle));
            y.set(GRB_DoubleAttr_Start, 0.0);
            _subproblem2bundle2vars[subproblem][bundle] = {x, y};
        }

        // Equality Constraints
        for (auto& kv : _decomp.subproblem2patches[subproblem])
        {
            FH p = kv.first;
            HFH hp0 = mcMesh.halfface_handle(p, 0);
            UVWDir dir = decompose(kv.second, DIM_1_DIRS)[0];

            auto& side2has = _decomp.hp2hasByDir.at(hp0);
            if (side2has.at(dir).size() == 1 && side2has.at(-dir).size() == 1
                && _decomp.arc2bundle.at(mcMesh.edge_handle(side2has.at(dir).front()))
                       == _decomp.arc2bundle.at(mcMesh.edge_handle(side2has.at(-dir).front())))
                continue;

            GRBLinExpr sum = 0;
            for (UVWDir side : {dir, -dir})
                for (HEH ha : side2has.at(side))
                // if (asVisited.count(mcMesh.edge_handle(ha)) != 0)
                {
                    auto& varPair = _subproblem2bundle2vars[subproblem].at(_decomp.arc2bundle.at(mcMesh.edge_handle(ha)));
                    sum += (side == dir ? 1 : -1) * (varPair.first - varPair.second);
                }

            std::string constraintName
                = "Patch" + std::to_string(p.idx()) + "dir" + std::to_string(static_cast<uint8_t>(dir | -dir));
            model.addConstr(sum, GRB_EQUAL, 0, constraintName);
        }

        _subproblem2model.emplace_back(std::move(mptr));
    }
}

void GurobiLPSolver::setupDynamicObjective(int subproblem)
{
    auto& model = *_subproblem2model[subproblem];
    GRBLinExpr objective = 0.;
    for (int bundle : _decomp.subproblem2bundles.at(subproblem))
    {
        auto& varPair = _subproblem2bundle2vars[subproblem].at(bundle);
        double weightAdd = weight(bundle, true);
        double weightSub = weight(bundle, false);
        objective += _decomp.bundle2arcs[bundle].size() * (weightAdd * varPair.first + weightSub * varPair.second);
    }
    model.setObjective(objective, GRB_MINIMIZE);
}

void GurobiLPSolver::setupPumpConstraints(int subproblem, int bundle, bool inflate)
{
    auto& model = *_subproblem2model[subproblem];
    auto& varPair = _subproblem2bundle2vars[subproblem].at(bundle);

    _temporaryConstraints.push_back(
        model.addConstr(inflate ? varPair.first : varPair.second, GRB_GREATER_EQUAL, 1.0, "CycleGenerator"));
    _temporaryConstraints.push_back(
        model.addConstr(inflate ? varPair.second : varPair.first, GRB_EQUAL, 0.0, "CycleGenerator2"));
}

void GurobiLPSolver::setupInflationOnlyConstraints(int subproblem)
{
    auto& model = *_subproblem2model[subproblem];
    for (int bundle : _decomp.subproblem2bundles.at(subproblem))
    {
        auto& varPair = _subproblem2bundle2vars[subproblem].at(bundle);
        _temporaryConstraints.push_back(model.addConstr(varPair.second, GRB_EQUAL, 0.0, "OnlyAdd"));
    }
}

bool GurobiLPSolver::setupSepNonnegConstraints(int subproblem)
{
    auto& model = *_subproblem2model[subproblem];
    bool allFulfilled = true;
    for (auto& aColl : _dynamicConstraints)
    {
        map<int, int> bundle2sign;
        for (auto sign2a : aColl)
        {
            if (sign2a.second.is_valid()
                && _decomp.subproblem2bundles.at(subproblem).count(_decomp.arc2bundle.at(sign2a.second)) != 0)
            {
                if (bundle2sign.count(_decomp.arc2bundle.at(sign2a.second)) == 0)
                    bundle2sign[_decomp.arc2bundle.at(sign2a.second)] = sign2a.first;
                else
                    bundle2sign[_decomp.arc2bundle.at(sign2a.second)] += sign2a.first;
            }
        }
        bool any = false;
        for (auto& kv : bundle2sign)
            if (kv.second != 0)
                any = true;
        if (!any)
            continue;
        int sum = 0;
        for (auto& sign2a : aColl)
            if (sign2a.second.is_valid())
                sum += sign2a.first * mcMeshProps().get<ARC_INT_LENGTH>(sign2a.second);
            else
                sum += 1 - sign2a.first;
        bool fulfilled = sum > 0;

        vector<GRBVar> vars;
        vector<double> coeff;
        double consts = 0;
        for (auto sign2a : aColl)
        {
            if (!sign2a.second.is_valid())
            {
                consts += 1 - sign2a.first;
                continue;
            }
            if (_decomp.subproblem2bundles.at(subproblem).count(_decomp.arc2bundle.at(sign2a.second)) == 0)
                continue;
            auto& varPair = _subproblem2bundle2vars[subproblem].at(_decomp.arc2bundle.at(sign2a.second));
            vars.push_back(varPair.first);
            coeff.push_back(sign2a.first);
            vars.push_back(varPair.second);
            coeff.push_back(sign2a.first * -1);
            consts += mcMeshProps().get<ARC_INT_LENGTH>(sign2a.second) * sign2a.first;
        }
        GRBLinExpr sumExpr = 0;
        sumExpr.addTerms(&coeff[0], &vars[0], vars.size());
        sumExpr += consts;

        std::string constraintName = "Separation";
        if (!fulfilled)
            allFulfilled = false;

        _temporaryConstraints.push_back(model.addConstr(sumExpr, GRB_GREATER_EQUAL, 1.0, constraintName));
    }
    return allFulfilled;
}

bool GurobiLPSolver::solve(int subproblem)
{
    auto& model = *_subproblem2model[subproblem];
    model.optimize();
    int status = model.get(GRB_IntAttr_Status);
    return status == 2 || status == 9;
}

double GurobiLPSolver::solution(int subproblem, int bundle)
{
    auto& varPair = _subproblem2bundle2vars[subproblem].at(bundle);
    double plus = varPair.first.get(GRB_DoubleAttr_X);
    double minus = varPair.second.get(GRB_DoubleAttr_X);
    return plus - minus;
}

void GurobiLPSolver::reconstrainToNextInt(int subproblem, int bundle)
{
    auto& model = *_subproblem2model[subproblem];
    auto& varPair = _subproblem2bundle2vars[subproblem].at(bundle);
    double plus = varPair.first.get(GRB_DoubleAttr_X);
    double minus = varPair.second.get(GRB_DoubleAttr_X);
    double val = plus - minus;
    _temporaryConstraints.push_back(model.addConstr(
        val > 0 ? varPair.first : varPair.second, GRB_EQUAL, std::ceil(std::abs(val)), "CycleGenerator"));
    _temporaryConstraints.push_back(
        model.addConstr(val > 0 ? varPair.second : varPair.first, GRB_EQUAL, 0.0, "CycleGenerator2"));
}

void GurobiLPSolver::removeTemporaryConstraints(int subproblem)
{
    auto& model = *_subproblem2model[subproblem];
    for (auto& constr : _temporaryConstraints)
        model.remove(constr);
    _temporaryConstraints.clear();
}

} // namespace impl
} // namespace qgp3d

#endif
