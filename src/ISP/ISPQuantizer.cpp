#include "QGP3D/ISP/ISPQuantizer.hpp"

#include "QGP3D/ObjectiveBuilder.hpp"

#ifdef QGP3D_WITH_GUROBI
#include "QGP3D/ISP/GurobiLPSolver.hpp"
#else
#include "QGP3D/ISP/ClpLPSolver.hpp"
#endif

#include <Eigen/Sparse>

namespace qgp3d
{

ISPQuantizer::ISPQuantizer(TetMeshProps& meshProps_, StructurePreserver& sep, ObjectiveFunction& obj)
    : TetMeshNavigator(meshProps_), TetMeshManipulator(meshProps_), MCMeshNavigator(meshProps_),
      MCMeshManipulator(meshProps_), _sep(sep), _obj(obj)
{
    decomposeIntoSubproblems();
}

namespace
{

struct BundlePrio
{
    BundlePrio(int _i, double _prio, bool _add, int _state) : i(_i), prio(_prio), add(_add), state(_state)
    {
    }

    int i;
    double prio;
    bool add;
    int state;
};

struct LeastWeightCompare
{
    bool operator()(const BundlePrio& b1, const BundlePrio& b2) const
    {
        // When true is returned, b2 gets popped before b1
        return b1.prio < b2.prio;
    }
};
using BundleQueue = std::priority_queue<BundlePrio, std::deque<BundlePrio>, LeastWeightCompare>;

struct BiBundlePrio
{
    // Assume adding i and subtracting j
    BiBundlePrio(int _i, int _j, double _prio, int _state) : i(_i), j(_j), prio(_prio), state(_state)
    {
    }

    int i;
    int j;
    double prio;
    int state;
};

struct DoubleLeastWeightCompare
{
    bool operator()(const BiBundlePrio& b1, const BiBundlePrio& b2) const
    {
        // When true is returned, b2 gets popped before b1
        return b1.prio < b2.prio;
    }
};
using BiBundleQueue = std::priority_queue<BiBundlePrio, std::deque<BiBundlePrio>, DoubleLeastWeightCompare>;

} // namespace

ISPQuantizer::RetCode ISPQuantizer::quantize(double varLowerBound)
{
    auto& mcMesh = mcMeshProps().mesh();

    if (!mcMeshProps().isAllocated<ARC_INT_LENGTH>())
        mcMeshProps().allocate<ARC_INT_LENGTH>(0);

    // Create and setup LP solver
#ifdef QGP3D_WITH_GUROBI
    impl::GurobiLPSolver sheetFinder(meshProps(), _decomp, _currentGrad, _currentHess, _currentContinuousQuadraticOpt);
#else
    impl::ClpLPSolver sheetFinder(meshProps(), _decomp, _currentGrad, _currentHess, _currentContinuousQuadraticOpt);
#endif
    sheetFinder.setupLPBase();

    // Main algo
    vector<vector<pair<int, EH>>> dynamicConstraints;
    vector<vector<pair<int, EH>>> simpleDynamicConstraints;

    cacheGradHess(_obj);

    double currentObj = greedyDescent(sheetFinder, false, false);

    int iter = 0;

    bool validatedSolution = false;
    bool recheckAfterThoroughSolve = false;
    bool uncheckedChanges = true;
    while (!validatedSolution || recheckAfterThoroughSolve)
    {
        vector<vector<pair<int, EH>>> newConstraints;
        _sep.violatedSimpleConstraints(varLowerBound, newConstraints);
        simpleDynamicConstraints.insert(simpleDynamicConstraints.end(), newConstraints.begin(), newConstraints.end());

        if (newConstraints.empty())
        {
            if (dynamicConstraints.empty() && !_sep.previousStructuralConstraints().empty())
                newConstraints = _sep.previousStructuralConstraints();
            else
            {
                if (uncheckedChanges)
                {
                    _sep.violatedStructuralConstraints(newConstraints);
                    uncheckedChanges = false;
                }
            }
        }

        validatedSolution = newConstraints.empty();
        DLOG(INFO) << "Iter " << iter << " solution valid? " << newConstraints.empty();

        if (!validatedSolution)
        {
            dynamicConstraints.insert(dynamicConstraints.end(), newConstraints.begin(), newConstraints.end());
            sheetFinder.setDynamicConstraints(dynamicConstraints);
            currentObj = makeFeasible(sheetFinder);

            if (numViolatedConstraints(sheetFinder.dynamicConstraints()) > 0)
            {
                LOG(WARNING) << "Switching to failsafe separation, due to no feasible integer solution found with "
                                "standard separation";
                vector<vector<pair<int, EH>>> failsafeDynamicConstraints = _sep.previousFailsafeStructuralConstraints();
                failsafeDynamicConstraints.insert(
                    failsafeDynamicConstraints.end(), simpleDynamicConstraints.begin(), simpleDynamicConstraints.end());
                sheetFinder.setDynamicConstraints(failsafeDynamicConstraints);
                currentObj = makeFeasible(sheetFinder);

                if (numViolatedConstraints(sheetFinder.dynamicConstraints()) > 0)
                    throw std::logic_error("Should be feasible now");

                sheetFinder.setDynamicConstraints(dynamicConstraints);
            }

            currentObj = greedyDescent(sheetFinder, true, false);

            uncheckedChanges = true;
            iter++;
        }

        if (validatedSolution && !recheckAfterThoroughSolve)
        {
            double objNew = greedyDescent(sheetFinder, true, true);
            if (objNew != currentObj)
                uncheckedChanges = true;
            currentObj = objNew;
            recheckAfterThoroughSolve = true;
        }
        else
            recheckAfterThoroughSolve = false;
    }
    DLOG(INFO) << "Needed " << iter << " sheet pump iterations to separate critical entities";
    DLOG(INFO) << "Final obj after separating critical arcs: " << currentObj;

    int minArcLength = 10000;
    for (EH a : mcMesh.edges())
        minArcLength = std::min(minArcLength, mcMeshProps().get<ARC_INT_LENGTH>(a));
    LOG(INFO) << "Sheet-pump-quantized MC, nHexes: " << _sep.numHexesInQuantization() << ", iterations: " << iter + 1
              << ", minArcLength " << minArcLength << ", objective: " << currentObj;

    return SUCCESS;
}

void ISPQuantizer::decomposeIntoSubproblems()
{
    auto& mcMesh = mcMeshProps().mesh();

    // Mini precompute
    _decomp.a2idx.clear();
    _decomp.idx2a.clear();
    for (EH a : mcMesh.edges())
    {
        _decomp.a2idx[a] = _decomp.a2idx.size();
        _decomp.idx2a.push_back(a);
    }

    _decomp.hp2hasByDir.clear();
    for (HFH hp : mcMesh.halffaces())
        _decomp.hp2hasByDir[hp] = halfpatchHalfarcsByDir(hp);

    // Group arcs into arc bundles constrained to be of equal length
    set<EH> bundledArcs;
    for (auto aStart : mcMesh.edges())
    {
        if (bundledArcs.count(aStart) != 0)
            continue;
        set<pair<FH, UVWDir>> psVisited;
        set<EH> asVisited({aStart});
        list<EH> aQ({aStart});
        while (!aQ.empty())
        {
            EH aCurr = aQ.front();
            aQ.pop_front();

            for (FH p : mcMesh.edge_faces(aCurr))
            {
                HFH hp0 = mcMesh.halfface_handle(p, 0);
                auto& side2has = _decomp.hp2hasByDir[hp0];

                UVWDir dirOpp = -findMatching(side2has,
                                              [&](const pair<const UVWDir, vector<HEH>>& kv)
                                              { return containsSomeOf(kv.second, mcMesh.edge_halfedges(aCurr)); })
                                     .first;

                if (side2has.at(dirOpp).size() != 1 || side2has.at(-dirOpp).size() != 1)
                    continue;
                if (psVisited.count({p, dirOpp | -dirOpp}) != 0)
                    continue;

                psVisited.insert({p, dirOpp | -dirOpp});

                for (HEH haNext : side2has.at(dirOpp))
                {
                    EH aNext = mcMesh.edge_handle(haNext);
                    if (asVisited.count(aNext) != 0)
                        continue;
                    aQ.push_back(aNext);
                    asVisited.insert(aNext);
                }
            }
        }
        for (auto a : asVisited)
            _decomp.arc2bundle[a] = _decomp.bundle2arcs.size();
        _decomp.bundle2arcs.emplace_back(asVisited);
        bundledArcs.insert(asVisited.begin(), asVisited.end());
    }

    _decomp.bundle2subproblem.resize(_decomp.bundle2arcs.size());

    // Group bundles into independent subproblems
    map<EH, int> a2multiBundle;
    set<int> multiBundled;
    int nMultiBundles = 0;
    for (int i = 0; i < (int)_decomp.bundle2arcs.size(); i++)
    {
        if (multiBundled.count(i) != 0)
            continue;
        DLOG(INFO) << "Gathering subproblem for bundle " << i;
        auto& bundleArcs = _decomp.bundle2arcs[i];
        EH aStart = *bundleArcs.begin();

        // One integer length variable per bundle
        set<pair<FH, UVWDir>> psVisited;
        set<EH> asVisited({aStart});
        list<EH> aQ({aStart});
        while (!aQ.empty())
        {
            EH aCurr = aQ.front();
            aQ.pop_front();

            for (FH p : mcMesh.edge_faces(aCurr))
            {
                HFH hp0 = mcMesh.halfface_handle(p, 0);
                auto& side2has = _decomp.hp2hasByDir[hp0];

                UVWDir dirOpp = -findMatching(side2has,
                                              [&](const pair<const UVWDir, vector<HEH>>& kv)
                                              { return containsSomeOf(kv.second, mcMesh.edge_halfedges(aCurr)); })
                                     .first;
                if (psVisited.count({p, dirOpp | -dirOpp}) != 0)
                    continue;

                psVisited.insert({p, dirOpp | -dirOpp});

                for (HEH haNext : side2has.at(dirOpp))
                {
                    EH aNext = mcMesh.edge_handle(haNext);
                    if (asVisited.count(aNext) != 0)
                        continue;
                    aQ.push_back(aNext);
                    asVisited.insert(aNext);
                }
                for (HEH haNext : side2has.at(-dirOpp))
                {
                    EH aNext = mcMesh.edge_handle(haNext);
                    if (asVisited.count(aNext) != 0)
                        continue;
                    aQ.push_back(aNext);
                    asVisited.insert(aNext);
                }
            }
        }
        set<int> localBundles;
        for (EH a : asVisited)
        {
            localBundles.insert(_decomp.arc2bundle.at(a));
            a2multiBundle[a] = nMultiBundles;
        }
        multiBundled.insert(localBundles.begin(), localBundles.end());
        nMultiBundles++;

        for (int j : localBundles)
        {
            _decomp.bundle2subproblem[j] = _decomp.subproblem2bundles.size();
        }
        _decomp.subproblem2patches.push_back(psVisited);
        _decomp.subproblem2bundles.push_back(localBundles);
    }

    // Make one globally connected problem
    set<pair<FH, UVWDir>> patches;
    set<int> bundles;
    for (auto& psVisited : _decomp.subproblem2patches)
        patches.insert(psVisited.begin(), psVisited.end());
    for (auto& bs : _decomp.subproblem2bundles)
        bundles.insert(bs.begin(), bs.end());
    _decomp.subproblem2patches.push_back(patches);
    _decomp.subproblem2bundles.push_back(bundles);

    DLOG(INFO) << "SUBPROBLEMS: " << nMultiBundles << " + 1 global problem instance";
}

double ISPQuantizer::greedyDescent(BaseLPSolver& sheetFinder, bool useGlobalProblem, bool thoroughSolve)
{
    double currentObj = obj();

    auto getPriority = [&, this](int i)
    {
        double d1 = 0.0;
        {
            auto a = *_decomp.bundle2arcs[i].begin();
            double xopt = _currentContinuousQuadraticOpt[_decomp.a2idx.at(a)];
            int xcurr = mcMeshProps().get<ARC_INT_LENGTH>(a);
            d1 = xopt - xcurr;
        }
        return d1;
    };

    int skipped = 0;

    bool redoWithBiBundleSheets = true;
    bool improving = true;
    while (improving)
    {
        while (improving)
        {
            improving = false;

            vector<int> bundle2lastState(_decomp.bundle2arcs.size(), 0);
            vector<vector<int>> w2bundle2sheetState(2, vector<int>(_decomp.bundle2arcs.size(), -1));
            vector<vector<map<int, double>>> w2bundle2sheet(2, vector<map<int, double>>(_decomp.bundle2arcs.size()));
            vector<vector<map<int, double>>> w2bundle2sheetOpp(2, vector<map<int, double>>(_decomp.bundle2arcs.size()));

            BundleQueue bQ;
            for (int i = 0; i < (int)_decomp.bundle2arcs.size(); i++)
            {
                double prio = getPriority(i);
                bQ.push(BundlePrio(i, std::abs(prio), !std::signbit(prio), 0));
            }

            while (!bQ.empty() && (thoroughSolve || (double)skipped < std::max(0.5 * _decomp.bundle2arcs.size(), 30.0)))
            {
                BundlePrio top = bQ.top();
                bQ.pop();
                int i = top.i;
                if (bundle2lastState[i] != top.state)
                    continue;

                // Sheet testing
                map<int, double> improvingSheet;
                for (bool alternativeWeights : {false, true})
                {
                    // Lazy sheet update
                    if (w2bundle2sheetState[alternativeWeights][i] != bundle2lastState[i])
                    {
                        sheetFinder.useAlternativeWeights(alternativeWeights);
                        w2bundle2sheet[alternativeWeights][i]
                            = sheetFinder.integerSheet(i, top.add, false, useGlobalProblem);
                        if (thoroughSolve)
                            w2bundle2sheetOpp[alternativeWeights][i]
                                = sheetFinder.integerSheet(i, !top.add, false, useGlobalProblem);
                        w2bundle2sheetState[alternativeWeights][i] = bundle2lastState[i];
                    }
                    for (bool add : {top.add, !top.add})
                    {
                        if (!thoroughSolve && add != top.add)
                            continue;
                        auto sheet = (add == top.add) ? w2bundle2sheet[alternativeWeights][i]
                                                      : w2bundle2sheetOpp[alternativeWeights][i];

                        if (sheet.empty())
                            continue;

                        double deltaObj = objDelta(sheet);
                        if (deltaObj < -1e-6)
                        {
                            DLOG(INFO) << "Found a " << (alternativeWeights ? "alternative " : "")
                                      << "sheet thats worth " << (add ? "adding." : "subtracting.") << " currentObj "
                                      << currentObj << " newObj " << currentObj + deltaObj;
                            // Execute update
                            currentObj = currentObj + deltaObj;
                            improvingSheet = sheet;
                            break;
                        }
                    }
                    if (!improvingSheet.empty())
                        break;
                }

                // Applying sheet
                if (!improvingSheet.empty())
                {
                    improving = true;
                    redoWithBiBundleSheets = true;

                    applyChange(improvingSheet);

                    // This invalidates queue prio but timestamps below guarantee old entries are rejected
                    cacheGradHess(_obj);

                    std::set<int> affBundles;
                    for (bool alternativeWeights : {false, true})
                    {
                        for (int i2 : affectedBundles(
                                 improvingSheet, w2bundle2sheet[alternativeWeights], sheetFinder.dynamicConstraints()))
                            affBundles.insert(i2);
                        if (thoroughSolve)
                            for (int i2 : affectedBundles(improvingSheet,
                                                          w2bundle2sheetOpp[alternativeWeights],
                                                          sheetFinder.dynamicConstraints()))
                                affBundles.insert(i2);
                    }
                    for (int i2 : affBundles)
                    {
                        double prio = getPriority(i2);
                        bQ.push(BundlePrio(i2, std::abs(prio), !std::signbit(prio), ++bundle2lastState[i2]));
                    }
                    skipped = 0;
                }
                else
                    skipped++;
            }
        }

        if (thoroughSolve && redoWithBiBundleSheets)
        {
            redoWithBiBundleSheets = false;
            DLOG(INFO) << "Bi-Bundle Sheeting";
            // Assume inflating i, deflating j
            auto getBiPriority
                = [&, this](int i, int j) { return sheetFinder.optimalDoubleFactor({{i, 1}, {j, -1}}); };

            vector<map<pairTT<int>, map<int, double>>> w2biBundle2sheet(2);
            map<pairTT<int>, int> biBundle2lastState;
            map<pairTT<int>, int> biBundle2sheetState;

            BiBundleQueue bQ;
            set<pairTT<int>> biBundles;
            for (int i = 0; i < _currentHess.outerSize(); ++i)
                for (Eigen::SparseMatrix<double>::InnerIterator it(_currentHess, i); it; ++it)
                {
                    if (it.col() > it.row())
                    {
                        int bundle1 = _decomp.arc2bundle[_decomp.idx2a[it.row()]];
                        int bundle2 = _decomp.arc2bundle[_decomp.idx2a[it.col()]];
                        if (bundle1 != bundle2)
                            biBundles.insert({std::min(bundle1, bundle2), std::max(bundle1, bundle2)});
                    }
                }

            for (auto& [i, j] : biBundles)
            {
                biBundle2lastState[{i, j}] = 0;
                biBundle2lastState[{j, i}] = 0;
                biBundle2sheetState[{i, j}] = -1;
                biBundle2sheetState[{j, i}] = -1;
                bQ.push(BiBundlePrio(i, j, getBiPriority(i, j), 0));
                bQ.push(BiBundlePrio(j, i, getBiPriority(j, i), 0));
            }

            while (!bQ.empty())
            {
                BiBundlePrio top = bQ.top();
                bQ.pop();
                int i = top.i;
                int j = top.j;
                if (biBundle2lastState[{i, j}] != top.state)
                    continue;

                if (biBundle2sheetState[{i, j}] != biBundle2lastState[{i, j}])
                {
                    for (bool alternativeWeights : {false, true})
                    {
                        sheetFinder.useAlternativeWeights(alternativeWeights);
                        w2biBundle2sheet[alternativeWeights][{i, j}] = sheetFinder.integerSheet(i, j);
                    }
                    biBundle2sheetState[{i, j}] = biBundle2lastState[{i, j}];
                }

                map<int, double> improvingSheet;

                for (bool alternativeWeights : {false, true})
                {
                    auto sheet = w2biBundle2sheet[alternativeWeights][{i, j}];

                    if (sheet.empty())
                        continue;

                    double deltaObj = objDelta(sheet);
                    if (deltaObj < -1e-6)
                    {
                        DLOG(INFO) << "Found a " << (alternativeWeights ? "alternative " : "")
                                  << "bi-bundle sheet thats worth " << "adding."
                                  << " currentObj " << currentObj << " newObj " << currentObj + deltaObj;
                        // Execute update
                        currentObj = currentObj + deltaObj;
                        improvingSheet = sheet;
                        break;
                    }
                }
                // Applying sheet
                if (!improvingSheet.empty())
                {
                    improving = true;

                    applyChange(improvingSheet);

                    // This invalidates queue prio but timestamps below guarantee old entries are rejected
                    cacheGradHess(_obj);

                    std::set<int> affBundles;
                    for (bool alternativeWeights : {false, true})
                    {
                        auto set1 = affectedBundles(
                            improvingSheet, w2biBundle2sheet[alternativeWeights], sheetFinder.dynamicConstraints());
                        for (auto i2 : set1)
                            affBundles.insert(i2);
                    }
                    for (auto& [nextI, nextJ] : biBundles)
                    {
                        if (affBundles.count(nextI) || affBundles.count(nextJ))
                        {
                            bQ.push(BiBundlePrio(nextI,
                                                     nextJ,
                                                     getBiPriority(nextI, nextJ),
                                                     ++biBundle2lastState[{nextI, nextJ}]));
                            bQ.push(BiBundlePrio(nextJ,
                                                     nextI,
                                                     getBiPriority(nextJ, nextI),
                                                     ++biBundle2lastState[{nextJ, nextI}]));
                        }
                    }
                    skipped = 0;
                }
                else
                    skipped++;
            }
        }
    }

    return currentObj;
}

double ISPQuantizer::makeFeasible(BaseLPSolver& sheetFinder)
{
    double currentObj = obj();
    sheetFinder.useAlternativeWeights(false);

    {
        int nUnfulfilled = numViolatedConstraints(sheetFinder.dynamicConstraints());
        DLOG(INFO) << nUnfulfilled << " separation constraints unfulfilled, inflating separation sheets";
    }

    vector<map<int, double>> bundle2sheet(_decomp.bundle2arcs.size());
    vector<int> bundle2lastState(_decomp.bundle2arcs.size(), 0);
    vector<int> bundle2sheetState(_decomp.bundle2arcs.size(), -1);

    bool change = true;
    while (change)
    {
        change = false;
        set<int> handledBundles;
        int bestBundle = -1;
        double bestDeltaObj = DBL_MAX;
        double bestRelDeltaObj = DBL_MAX;

        map<int, double> bestSheet = sheetFinder.integerSheet(0, false, true, true);
        if (bestSheet.empty())
            bestSheet = sheetFinder.integerSheet(0, true, true, true);
        if (!bestSheet.empty())
        {
            bestBundle = 0;
            bundle2sheet[0] = bestSheet;
            double deltaObj = objDelta(bestSheet);
            bestDeltaObj = deltaObj;
            bestRelDeltaObj = deltaObj / numViolatedConstraints(sheetFinder.dynamicConstraints());
        }
        else
        {
            // For each subproblem individually check if constraints violated and fixable
            // Avoid global problem instance here, which is last in _decomp.subproblem2bundles
            for (int subproblem = 0; subproblem < (int)_decomp.subproblem2bundles.size() - 1; subproblem++)
            {
                int bundle = *_decomp.subproblem2bundles[subproblem].begin();
                if (bundle2sheetState[bundle] != bundle2lastState[bundle])
                {
                    bundle2sheet[bundle] = sheetFinder.integerSheet(bundle, false, true, false);
                    bundle2sheetState[bundle] = bundle2lastState[bundle];
                }
                // Determine sheet:
                auto& sheet = bundle2sheet[bundle];

                int feasibilityImprovement = 0;
                for (bool tryAddOnly : {false, true})
                {
                    if (tryAddOnly && !sheet.empty())
                    {
                        int maxDenominator = 0;
                        for (auto& kv : sheet)
                            maxDenominator = std::max(std::abs((int)std::round(kv.second)), maxDenominator);
                        if (maxDenominator > 1)
                            sheet = sheetFinder.integerSheet(bundle, true, true, false);
                        else
                            continue;
                    }
                    if (sheet.empty())
                        continue;

                    // Does sheet fulfill all its associated separation constraints?

                    int nUnfulfilledPre = numViolatedConstraints(sheetFinder.dynamicConstraints());
                    applyChange(sheet);
                    int nUnfulfilled = numViolatedConstraints(sheetFinder.dynamicConstraints());
                    // revert
                    applyChange(sheet, false);

                    if (nUnfulfilled < nUnfulfilledPre)
                    {
                        feasibilityImprovement = nUnfulfilledPre - nUnfulfilled;
                        break;
                    }
                }
                if (feasibilityImprovement <= 0)
                    continue;

                double deltaObj = objDelta(sheet);
                if (deltaObj / feasibilityImprovement < bestRelDeltaObj)
                {
                    bestBundle = bundle;
                    bestDeltaObj = deltaObj;
                    bestRelDeltaObj = deltaObj / feasibilityImprovement;
                }
            }
        }

        if (bestBundle >= 0)
        {
            auto& sheet = bundle2sheet[bestBundle];
            applyChange(sheet);
            cacheGradHess(_obj);

            int nUnfulfilled = numViolatedConstraints(sheetFinder.dynamicConstraints());
            DLOG(INFO) << "Found a separating sheet, " << nUnfulfilled
                      << " separation constraints still unfulfilled. currentObj " << currentObj << " newObj "
                      << currentObj + bestDeltaObj;
            currentObj += bestDeltaObj;
            if (nUnfulfilled == 0)
                break;

            for (int i2 : affectedBundles(sheet, bundle2sheet, sheetFinder.dynamicConstraints()))
                bundle2lastState[i2]++;

            change = true;
        }
    }
    return currentObj;
}

int ISPQuantizer::numViolatedConstraints(const vector<vector<pair<int, EH>>>& constraints) const
{
    int nUnfulfilled = 0;
    for (auto& coll : constraints)
    {
        int sum = 0;
        for (auto& sign2a : coll)
            if (sign2a.second.is_valid())
                sum += sign2a.first * mcMeshProps().get<ARC_INT_LENGTH>(sign2a.second);
            else
                sum += (1 - sign2a.first);
        if (sum <= 0)
            nUnfulfilled++;
    }
    return nUnfulfilled;
}

double ISPQuantizer::obj() const
{
    Eigen::VectorXd x(mcMeshProps().mesh().n_logical_edges());
    int i = 0;
    for (EH a : mcMeshProps().mesh().edges())
        x(i++) = mcMeshProps().get<ARC_INT_LENGTH>(a);

    return _obj.function_value(x);
}

double ISPQuantizer::objDelta(const map<int, double>& sheet)
{
    applyChange(sheet);
    double newObj = obj();
    applyChange(sheet, false);
    return newObj - obj();
}

set<int> ISPQuantizer::affectedBundles(const map<int, double>& sheet,
                                       const vector<map<int, double>>& bundle2sheet,
                                       const vector<vector<pair<int, EH>>>& constraints) const
{
    set<int> bundleSet;
    for (auto& kv : sheet)
        bundleSet.insert(kv.first);
    for (auto& coll : constraints)
    {
        if (containsMatching(
                coll,
                [&bundleSet, this](const pair<const int, EH>& sign2a)
                { return sign2a.second.is_valid() && bundleSet.count(_decomp.arc2bundle.at(sign2a.second)) != 0; }))
            for (auto& sign2a : coll)
                if (sign2a.second.is_valid())
                    bundleSet.insert(_decomp.arc2bundle.at(sign2a.second));
    }
    set<int> bundleSet2;
    for (int i2 = 0; i2 < (int)bundle2sheet.size(); i2++)
    {
        if (containsMatching(bundle2sheet[i2],
                             [&bundleSet](const pair<const int, double>& kv)
                             { return bundleSet.count(kv.first) != 0; }))
            bundleSet2.insert(i2);
    }
    bundleSet2.insert(bundleSet.begin(), bundleSet.end());
    return bundleSet2;
}

set<int> ISPQuantizer::affectedBundles(const map<int, double>& sheet,
                                       const map<pairTT<int>, map<int, double>>& biBundle2sheet,
                                       const vector<vector<pair<int, EH>>>& constraints) const
{
    set<int> bundleSet;
    for (auto& kv : sheet)
        bundleSet.insert(kv.first);
    for (auto& coll : constraints)
    {
        if (containsMatching(
                coll,
                [&bundleSet, this](const pair<const int, EH>& sign2a)
                { return sign2a.second.is_valid() && bundleSet.count(_decomp.arc2bundle.at(sign2a.second)) != 0; }))
            for (auto& sign2a : coll)
                if (sign2a.second.is_valid())
                    bundleSet.insert(_decomp.arc2bundle.at(sign2a.second));
    }
    set<int> bundleSet2;
    for (auto& [biBundle, sheet2] : biBundle2sheet)
    {
        if (containsMatching(
                sheet2, [&bundleSet](const pair<const int, double>& kv) { return bundleSet.count(kv.first) != 0; }))
        {
            bundleSet2.insert(biBundle.first);
            bundleSet2.insert(biBundle.second);
        }
    }
    bundleSet2.insert(bundleSet.begin(), bundleSet.end());
    return bundleSet2;
}

void ISPQuantizer::cacheGradHess(ObjectiveFunction& obj)
{
    Eigen::VectorXd x(mcMeshProps().mesh().n_logical_edges());
    int i = 0;
    for (EH a : mcMeshProps().mesh().edges())
        x(i++) = mcMeshProps().get<ARC_INT_LENGTH>(a);
    obj.gradient(x, _currentGrad);
    obj.hessian(x, _currentHess);
    if (!obj.is_hessian_const() || _currentContinuousQuadraticOpt.rows() == 0)
    {
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> ldlt(_currentHess);
        if (ldlt.info() != Eigen::Success)
        {
            LOG(ERROR) << "Ill conditioned objective hessian";
            _currentContinuousQuadraticOpt = Eigen::VectorXd(_decomp.idx2a.size());
            for (i = 0; i < (int)_decomp.idx2a.size(); i++)
                _currentContinuousQuadraticOpt(i) = mcMeshProps().get<ARC_DBL_LENGTH>(_decomp.idx2a[i]);
        }
        else
            _currentContinuousQuadraticOpt = x + ldlt.solve(-_currentGrad);
    }
}

void ISPQuantizer::applyChange(const map<int, double>& sheet, bool add)
{
    for (auto& kv : sheet)
        for (EH a : _decomp.bundle2arcs[kv.first])
            mcMeshProps().ref<ARC_INT_LENGTH>(a) += (add ? 1 : -1) * std::round(kv.second);
}

} // namespace qgp3d
