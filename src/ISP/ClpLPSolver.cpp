#ifndef QGP3D_WITH_GUROBI

#include "QGP3D/ISP/ClpLPSolver.hpp"

namespace qgp3d
{
namespace impl
{

ClpLPSolver::ClpLPSolver(const TetMeshProps& meshProps, double scaling, const ISPQuantizer::Decomposition& decomp)
    : TetMeshNavigator(meshProps), MCMeshNavigator(meshProps), BaseLPSolver(meshProps, scaling, decomp)
{
}

void ClpLPSolver::setupLPBase()
{
    auto& mcMesh = mcMeshProps().mesh();
    _subproblem2bundle2vars.resize(_decomp.subproblem2bundles.size());
    for (int subproblem = 0; subproblem < (int)_decomp.subproblem2bundles.size(); subproblem++)
    {
        DLOG(INFO) << "Creating model for subproblem " << subproblem;
        std::unique_ptr<ClpSimplex> mptr = std::make_unique<ClpSimplex>();
        ClpSimplex& model = *mptr;
#ifdef NDEBUG
        model.setLogLevel(0);
#endif

        DLOG(INFO) << "Creating " << _decomp.subproblem2bundles[subproblem].size() * 2 << " variables for subproblem " << subproblem;
        model.resize(0, 2 * _decomp.subproblem2bundles[subproblem].size());
        int currentIndex = 0;

        for (int j : _decomp.subproblem2bundles[subproblem])
        {
            _subproblem2bundle2vars[subproblem][j] = {currentIndex, currentIndex + 1};
            // Add variable
            model.setColumnLower(currentIndex, 0.0);
            model.setColumnUpper(currentIndex, COIN_DBL_MAX);
            currentIndex++;
            // Subtract variable
            model.setColumnLower(currentIndex, 0.0);
            model.setColumnUpper(currentIndex, COIN_DBL_MAX);
            currentIndex++;
        }

        // Equality Constraints
        for (auto& kv : _decomp.subproblem2patches[subproblem])
        {
            FH p = kv.first;
            HFH hp0 = mcMesh.halfface_handle(p, 0);
            UVWDir dir = decompose(kv.second, DIM_1_DIRS)[0];

            auto& side2has = _decomp.hp2hasByDir[hp0];
            if (side2has[dir].size() == 1 && side2has[-dir].size() == 1
                && _decomp.arc2bundle[mcMesh.edge_handle(side2has[dir].front())]
                       == _decomp.arc2bundle[mcMesh.edge_handle(side2has[-dir].front())])
                continue;

            vector<int> vars;
            vector<double> coeffs;
            for (UVWDir side : {dir, -dir})
                for (HEH ha : side2has[side])
                    // if (asVisited.count(mcMesh.edge_handle(ha)) != 0)
                    {
                        auto& varPair = _subproblem2bundle2vars[subproblem].at(_decomp.arc2bundle[mcMesh.edge_handle(ha)]);
                        vars.push_back(varPair.first);
                        coeffs.push_back(side == dir ? 1.0 : -1.0);
                        vars.push_back(varPair.second);
                        coeffs.push_back(side == dir ? -1.0 : 1.0);
                    }

            model.addRow(vars.size(), vars.data(), coeffs.data(), 0.0, 0.0);
        }
        _subproblem2numberOfRowsBase.push_back(mptr->numberRows());
        _subproblem2model.emplace_back(std::move(mptr));
    }
}

void ClpLPSolver::setupDynamicObjective(int subproblem)
{
    auto& model = *_subproblem2model[subproblem];
    for (int bundle : _decomp.subproblem2bundles.at(subproblem))
    {
        auto& varPair = _subproblem2bundle2vars[subproblem].at(bundle);
        double weightAdd = weight(bundle, true);
        double weightSub = weight(bundle, false);
        model.setObjectiveCoefficient(varPair.first, weightAdd * _decomp.bundle2arcs[bundle].size());
        model.setObjectiveCoefficient(varPair.second, weightSub * _decomp.bundle2arcs[bundle].size());
    }
}

void ClpLPSolver::setupPumpConstraints(int subproblem, int bundle, bool inflate)
{
    auto& model = *_subproblem2model[subproblem];
    auto& varPair = _subproblem2bundle2vars[subproblem].at(bundle);

    double coeff = 1.0;
    model.addRow(1, inflate ? &varPair.first : &varPair.second, &coeff, 1.0, COIN_DBL_MAX);
    model.addRow(1, inflate ? &varPair.second : &varPair.first, &coeff, 0.0, 0.0);
}

void ClpLPSolver::setupInflationOnlyConstraints(int subproblem)
{
    auto& model = *_subproblem2model[subproblem];
    for (int bundle : _decomp.subproblem2bundles.at(subproblem))
    {
        auto& varPair = _subproblem2bundle2vars[subproblem].at(bundle);
        double coeff = 1.0;
        model.addRow(1, &varPair.second, &coeff, 0.0, 0.0);
    }
}

bool ClpLPSolver::setupSepNonnegConstraints(int subproblem)
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

        vector<int> vars;
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
        if (!fulfilled)
            allFulfilled = false;

        model.addRow(vars.size(), &vars[0], &coeff[0], 1.0 - consts, COIN_DBL_MAX);
    }
    return allFulfilled;
}

bool ClpLPSolver::solve(int subproblem)
{
    auto& model = *_subproblem2model[subproblem];
    model.primal();
    int status = model.status();
    return status == 0;
}

double ClpLPSolver::solution(int subproblem, int bundle)
{
    auto& model = *_subproblem2model[subproblem];
    double* sol = model.primalColumnSolution();
    auto& varPair = _subproblem2bundle2vars[subproblem].at(bundle);
    double plus = sol[varPair.first];
    double minus = sol[varPair.second];
    return plus - minus;
}

void ClpLPSolver::reconstrainToNextInt(int subproblem, int bundle)
{
    auto& model = *_subproblem2model[subproblem];
    double* sol = model.primalColumnSolution();
    auto& varPair = _subproblem2bundle2vars[subproblem].at(bundle);
    double plus = sol[varPair.first];
    double minus = sol[varPair.second];
    double val = plus - minus;
    double coeff = 1.0;
    model.addRow(
        1, val > 0 ? &varPair.first : &varPair.second, &coeff, std::ceil(std::abs(val)), std::ceil(std::abs(val)));
    model.addRow(1, val > 0 ? &varPair.second : &varPair.first, &coeff, 0.0, 0.0);
}

void ClpLPSolver::removeTemporaryConstraints(int subproblem)
{
    auto& model = *_subproblem2model[subproblem];

    vector<int> rows(model.numberRows() - _subproblem2numberOfRowsBase[subproblem]);
    std::iota(rows.begin(), rows.end(), _subproblem2numberOfRowsBase[subproblem]);
    model.deleteRows(rows.size(), rows.data());
}

} // namespace impl
} // namespace qgp3d

#endif
