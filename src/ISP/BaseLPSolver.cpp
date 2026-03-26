#include "QGP3D/ISP/BaseLPSolver.hpp"

namespace qgp3d
{

double BaseLPSolver::weight(int bundle, bool inflating)
{
    auto& mcMesh = mcMeshProps().mesh();

    map<int, double> partialSheet = {{bundle, inflating ? 1 : -1}};

    double weight = 0.0;
    if (_alternativeWeights)
    {
        double d1 = optimalDoubleFactor(partialSheet);
        EH a = *_decomp.bundle2arcs[bundle].begin();
        double xopt = _currentContinuousQuadraticOpt(_decomp.a2idx.at(a));
        int xcurr = mcMeshProps().get<ARC_INT_LENGTH>(a);
        double d2 = (inflating ? 1 : -1) * (xopt - xcurr);
        double d = std::max(d1, d2);
        if (1 <= d)
            weight = 1.0 / (d + 1);
        else if (0 <= d)
            weight = mcMesh.n_logical_edges() / (d + 1);
        else
            weight = mcMesh.n_logical_edges() * mcMesh.n_logical_edges() * (1-d);
    }
    else
    {
        EH a = *_decomp.bundle2arcs[bundle].begin();
        double xopt = _currentContinuousQuadraticOpt(_decomp.a2idx.at(a));
        int xcurr = mcMeshProps().get<ARC_INT_LENGTH>(a);
        double d = (inflating ? 1 : -1) * (xopt - xcurr);
        if (1 <= d)
            weight = 1.0 / (d + 1);
        else if (0 <= d)
            weight = mcMesh.n_logical_edges() / (d + 1);
        else
            weight = mcMesh.n_logical_edges() * mcMesh.n_logical_edges() * (1 - d);
    }

    return weight * _decomp.bundle2arcs[bundle].size();
}

map<int, double>
BaseLPSolver::integerSheet(int bundleToChange, bool add, bool fixing, bool useGlobalProblem, double dampingFactor)
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
    // Enforce integrality
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
                for (auto& kv : sheet)
                    kv.second *= lcd;
                break;
            }
        }
    }

    removeTemporaryConstraints(subproblem);

    // At this point sheet is pure integer, now we can find optimal integer step size which again yields integer
    int factor = optimalCleanedFactor(sheet, dampingFactor);
    if (factor == 0)
        return {};

    for (auto& kv : sheet)
        kv.second *= factor;
    return sheet;
}

map<int, double> BaseLPSolver::integerSheet(int bundleI, int bundleJ, double dampingFactor)
{
    map<int, double> sheet;
    int subproblem = _decomp.subproblem2bundles.size() - 1;

    setupSepNonnegConstraints(subproblem);

    setupPumpConstraints(subproblem, bundleI, true);
    setupPumpConstraints(subproblem, bundleJ, false);

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
                for (auto& kv : sheet)
                    kv.second *= lcd;
                break;
            }
        }
    }

    removeTemporaryConstraints(subproblem);

    // At this point sheet is pure integer, now we can find optimal integer step size which again yields integer
    int factor = optimalCleanedFactor(sheet, dampingFactor);
    if (factor == 0)
        return {};

    for (auto& kv : sheet)
        kv.second *= factor;
    return sheet;
}

void BaseLPSolver::useAlternativeWeights(bool altWeights)
{
    _alternativeWeights = altWeights;
}

double BaseLPSolver::optimalDoubleFactor(const map<int, double>& sheet) const
{
    Eigen::SparseVector<double> deltaQ(mcMeshProps().mesh().n_logical_edges());
    for (auto& kv : sheet)
        for (EH a : _decomp.bundle2arcs[kv.first])
            deltaQ.insert(_decomp.a2idx.at(a)) = kv.second;

    return -deltaQ.dot(_currentGrad) / deltaQ.dot(_currentHess * deltaQ);
}

int BaseLPSolver::optimalCleanedFactor(const map<int, double>& sheet, double dampingFactor) const
{
    double dblFactor = optimalDoubleFactor(sheet);
    dblFactor = std::max(std::min(0.5, dblFactor), dblFactor * dampingFactor);
    int optFactor = std::round(dblFactor);

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

    return std::clamp(optFactor, minFactor, maxFactor);
}

} // namespace qgp3d
