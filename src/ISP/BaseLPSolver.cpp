#include "QGP3D/ISP/BaseLPSolver.hpp"

namespace qgp3d
{

double BaseLPSolver::weight(int bundle, bool inflating)
{
    auto& mcMesh = mcMeshProps().mesh();

    EH a = *_decomp.bundle2arcs[bundle].begin();
    double xopt = mcMeshProps().get<ARC_DBL_LENGTH>(a) * _scaling;
    int xcurr = mcMeshProps().get<ARC_INT_LENGTH>(a);
    double d = (inflating ? 1 : -1) * (xopt - xcurr);
    if (1 <= d)
        return 1.0 / (d + 1);
    else if (0 <= d && d <= 1)
        return mcMesh.n_logical_edges() / (d + 1);
    else
        return mcMesh.n_logical_edges() * mcMesh.n_logical_edges() * (1 - d);
}

map<int, double> BaseLPSolver::integerSheet(int bundleToChange, bool add, bool fixing, bool useGlobalProblem)
{
    map<int, double> sheet;
    int subproblem = -1;

    if (useGlobalProblem)
        subproblem = _decomp.subproblem2bundles.size() - 1;
    else
        subproblem = _decomp.bundle2subproblem.at(bundleToChange);

    bool allFulfilled = setupSepNonnegConstraints(subproblem);
    if (fixing && allFulfilled)
    {
        removeTemporaryConstraints(subproblem);
        return sheet;
    }

    if (!fixing)
        setupPumpConstraints(subproblem, bundleToChange, add);
    else if (add)
        setupInflationOnlyConstraints(subproblem);

    setupDynamicObjective(subproblem);

    int nonZero = 0;

    {
        DLOG(INFO) << "Optimizing LP model";
        bool success = solve(subproblem);
        if (!success)
        {
            removeTemporaryConstraints(subproblem);
            return sheet;
        }
    }

    int lcd = 1;
    // Enforce integerness
    set<int> factors;
    bool solved = false;
    while (!solved)
    {
        sheet.clear();
        solved = true;
        lcd = 1;
        int maxBundle = -1;
        double minDelta = DBL_MAX;
        for (int bundle : _decomp.subproblem2bundles.at(subproblem))
        {
            double val = solution(subproblem, bundle);
            if (std::abs(val) > 1e-6)
            {
                sheet[bundle] = val;
                nonZero++;
                if (std::abs(std::round(val) - val) > 1e-6)
                {
                    solved = false;
                    double delta = std::ceil(std::abs(val)) - std::abs(val);
                    if (delta < minDelta)
                    {
                        minDelta = delta;
                        maxBundle = bundle;
                    }

                    int factor = -1;
                    for (int j = 2; j < 1024; j++)
                        if (std::abs(std::round(j * val) - j * val) < 1e-6)
                        {
                            factor = j;
                            break;
                        }
                    assert(factor >= 2);
                    lcd = std::lcm(lcd, factor);
                }
            }
        }

        if (!solved)
        {
            reconstrainToNextInt(subproblem, maxBundle);
            DLOG(INFO) << "Reoptimizing LP model";
            bool success = solve(subproblem);
            if (!success)
            {
                // instead scale by lcd
                for (auto& kv: sheet)
                    kv.second *= lcd;
                break;
            }
        }
    }

    removeTemporaryConstraints(subproblem);

    // At this point sheet is pure integer, now we can find optimal integer step size which again yields integer
    int factor = optimalFactor(sheet);
    if (factor == 0)
        return {};

    for (auto& kv: sheet)
        kv.second *= factor;
    return sheet;
}

int BaseLPSolver::optimalFactor(const map<int, double>& sheet) const
{
    double dqDotQL = 0.0;
    double dqSq = 0.0;
    for (auto& kv : sheet)
    {
        int j = kv.first;
        EH aj = *_decomp.bundle2arcs[j].begin();
        double d = kv.second;

        double xopt2 = mcMeshProps().get<ARC_DBL_LENGTH>(aj) * _scaling;
        int xcurr2 = mcMeshProps().get<ARC_INT_LENGTH>(aj);

        dqDotQL += d * (xcurr2 - xopt2)* _decomp.bundle2arcs[j].size();
        dqSq += d*d* _decomp.bundle2arcs[j].size();
    }
    int optDelta = std::round(-dqDotQL / dqSq);

    // Clamp to feasible region
    int maxFactor = INT_MAX;
    int minFactor = INT_MIN;
    for (auto& coll : _dynamicConstraints)
    {
        bool affected = false;
        int sum = 0;
        int rhs = 1;
        for (auto& sign2a : coll)
        {
            if (sign2a.second.is_valid())
            {
                if (sheet.count(_decomp.arc2bundle.at(sign2a.second)) != 0)
                {
                    affected = true;
                    sum += sign2a.first * (int)std::round(sheet.at(_decomp.arc2bundle.at(sign2a.second)));
                }
                rhs -= sign2a.first * mcMeshProps().get<ARC_INT_LENGTH>(sign2a.second);
            }
            else
                rhs -= (1 - sign2a.first);
        }
        if (affected)
        {
            double factor = (double)rhs / sum;
            if (sum < 0)
                maxFactor = std::min(maxFactor, (int)std::floor(factor));
            else if (sum > 0)
                minFactor = std::max(minFactor, (int)std::ceil(factor));
            else if (sum == 0)
            {
                if (rhs > 0)
                {
                    minFactor = INT_MAX;
                    maxFactor = INT_MIN;
                }
            }
            if (maxFactor < minFactor) // Infeasible
            {
                return 0;
            }
        }
    }

    return std::clamp(optDelta, minFactor, maxFactor);
}

} // namespace qgp3d
