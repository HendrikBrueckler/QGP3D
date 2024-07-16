#include "QGP3D/ISP/ISPQuantizer.hpp"

#ifdef QGP3D_WITH_GUROBI
#include "QGP3D/ISP/GurobiLPSolver.hpp"
#else
#include "QGP3D/ISP/ClpLPSolver.hpp"
#endif

namespace qgp3d
{

ISPQuantizer::ISPQuantizer(TetMeshProps& meshProps, SeparationChecker& sep)
    : TetMeshNavigator(meshProps), TetMeshManipulator(meshProps), MCMeshNavigator(meshProps),
      MCMeshManipulator(meshProps), _sep(sep)
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

struct GreatestWeightCompare
{
    bool operator()(const BundlePrio& b1, const BundlePrio& b2) const
    {
        // When true is returned, b1 gets lower priority than b2
        return b1.prio > b2.prio;
    }
};
using BundleQueue = std::priority_queue<BundlePrio, std::deque<BundlePrio>, GreatestWeightCompare>;

} // namespace

ISPQuantizer::RetCode ISPQuantizer::quantize(double scaling, double varLowerBound)
{
    auto& mcMesh = mcMeshProps().mesh();

    // Allocate properties
    bool wasAllocated = mcMeshProps().isAllocated<ARC_DBL_LENGTH>();
    if (!wasAllocated)
    {
        mcMeshProps().allocate<ARC_DBL_LENGTH>(0.0);
        for (EH arc : mcMesh.edges())
        {
            // Determine current arc length
            double length = 0.0;
            for (HEH he : mcMeshProps().ref<ARC_MESH_HALFEDGES>(arc))
                length += edgeLengthUVW<CHART>(meshProps().mesh().edge_handle(he));
            mcMeshProps().set<ARC_DBL_LENGTH>(arc, length);
        }
    }

    if (!mcMeshProps().isAllocated<ARC_INT_LENGTH>())
        mcMeshProps().allocate<ARC_INT_LENGTH>(0);

        // Create and setup LP solver
#ifdef QGP3D_WITH_GUROBI
    impl::GurobiLPSolver sheetFinder(meshProps(), scaling, _decomp);
#else
    impl::ClpLPSolver sheetFinder(meshProps(), scaling, _decomp);
#endif
    sheetFinder.setupLPBase();

    // Compute critical link structure
    vector<CriticalLink> criticalLinks;
    map<EH, int> a2criticalLinkIdx;
    map<VH, vector<int>> n2criticalLinksOut;
    map<VH, vector<int>> n2criticalLinksIn;
    getCriticalLinks(criticalLinks, a2criticalLinkIdx, n2criticalLinksOut, n2criticalLinksIn, true);

    vector<bool> isCriticalArc(mcMesh.n_edges(), false);
    vector<bool> isCriticalNode(mcMesh.n_vertices(), false);
    vector<bool> isCriticalPatch(mcMesh.n_faces(), false);
    for (auto& kv : a2criticalLinkIdx)
        isCriticalArc[kv.first.idx()] = true;
    for (VH n : mcMesh.vertices())
    {
        auto type = mcMeshProps().nodeType(n);
        if (type.first == SingularNodeType::SINGULAR || type.second == FeatureNodeType::FEATURE
            || type.second == FeatureNodeType::SEMI_FEATURE_SINGULAR_BRANCH)
            isCriticalNode[n.idx()] = true;
    }
    for (FH p : mcMesh.faces())
        isCriticalPatch[p.idx()] = mcMesh.is_boundary(p)
                                   || (mcMeshProps().isAllocated<IS_FEATURE_F>() && mcMeshProps().get<IS_FEATURE_F>(p));

    // Main algo
    vector<vector<pair<int, EH>>> dynamicConstraints;
    vector<vector<pair<int, EH>>> simpleDynamicConstraints;

    double currentObj = greedyDescent(sheetFinder, scaling, false);

    DLOG(INFO) << "Obj before separation checking: " << currentObj;

    int iter = 0;

    bool validatedSolution = false;
    bool useGlobalProblem = false;
    while (!validatedSolution || useGlobalProblem)
    {
        if (useGlobalProblem)
            currentObj = greedyDescent(sheetFinder, scaling, true);

        auto constraints = violatedSimpleConstraints(varLowerBound, criticalLinks);
        simpleDynamicConstraints.insert(simpleDynamicConstraints.end(), constraints.begin(), constraints.end());

        if (constraints.empty())
        {
            if (dynamicConstraints.empty() && !_sep.previousSeparationViolatingPaths().empty())
                constraints = _sep.previousSeparationViolatingPaths();
            else
            {
                _sep.findSeparationViolatingPaths(
                    criticalLinks, isCriticalArc, isCriticalNode, isCriticalPatch, constraints);
                DLOG(INFO) << "Found unseparated features? " << !constraints.empty();
            }
        }

        validatedSolution = constraints.empty();
        DLOG(INFO) << "Iter " << iter << " solution valid? " << constraints.empty();

        if (!validatedSolution)
        {
            useGlobalProblem = false;
            dynamicConstraints.insert(dynamicConstraints.end(), constraints.begin(), constraints.end());
            sheetFinder.setDynamicConstraints(dynamicConstraints);
            currentObj = makeFeasible(sheetFinder, scaling);

            if (numViolatedConstraints(sheetFinder.dynamicConstraints()) > 0)
            {
                DLOG(WARNING) << "Switching to failsafe separation, due to no feasible integer solution found with "
                                 "standard separation";
                vector<vector<pair<int, EH>>> failsafeDynamicConstraints = _sep.failsafeSeparationViolatingPaths();
                failsafeDynamicConstraints.insert(
                    failsafeDynamicConstraints.end(), simpleDynamicConstraints.begin(), simpleDynamicConstraints.end());
                sheetFinder.setDynamicConstraints(failsafeDynamicConstraints);
                currentObj = makeFeasible(sheetFinder, scaling);

                sheetFinder.setDynamicConstraints(dynamicConstraints);
            }

            currentObj = greedyDescent(sheetFinder, scaling, false);

            iter++;
        }
        else if (!dynamicConstraints.empty())
            useGlobalProblem = !useGlobalProblem;
    }
    DLOG(INFO) << "Needed " << iter << " sheet pump iterations to separate critical entities";
    DLOG(INFO) << "Final obj after separating critical arcs: " << currentObj;

    int minArcLength = 10000;
    for (EH a : mcMesh.edges())
        minArcLength = std::min(minArcLength, mcMeshProps().get<ARC_INT_LENGTH>(a));
    LOG(INFO) << "Sheet-pump-quantized MC, nHexes: " << _sep.numHexesInQuantization() << ", objective: " << currentObj
              << ", iterations: " << iter + 1 << ", minArcLength " << minArcLength << std::endl;

    return SUCCESS;
}

void ISPQuantizer::decomposeIntoSubproblems()
{
    auto& mcMesh = mcMeshProps().mesh();
    // Mini precompute
    _decomp.hp2hasByDir.clear();
    for (auto hp : mcMesh.halffaces())
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

double ISPQuantizer::greedyDescent(BaseLPSolver& sheetFinder, double scaling, bool useGlobalProblem)
{
    double currentObj = obj(scaling);

    auto getPriority = [&, this](int i)
    {
        auto a = *_decomp.bundle2arcs[i].begin();
        double xopt = mcMeshProps().get<ARC_DBL_LENGTH>(a) * scaling;
        int xcurr = mcMeshProps().get<ARC_INT_LENGTH>(a);
        return xopt - xcurr;
    };

    int skipped = 0;

    bool improvement = true;
    while (improvement)
    {
        improvement = false;

        vector<map<int, double>> bundle2sheet(_decomp.bundle2arcs.size());
        vector<int> bundle2lastState(_decomp.bundle2arcs.size(), 0);
        vector<int> bundle2sheetState(_decomp.bundle2arcs.size(), -1);

        BundleQueue bQ;
        for (int i = 0; i < (int)_decomp.bundle2arcs.size(); i++)
        {
            double prio = getPriority(i);
            bQ.push(BundlePrio(i, std::abs(prio), !std::signbit(prio), 0));
        }

        while (!bQ.empty() && (double)skipped < std::max(0.5 * _decomp.bundle2arcs.size(), 30.0))
        {
            BundlePrio top = bQ.top();
            bQ.pop();
            int i = top.i;
            if (bundle2lastState[i] != top.state)
                continue;

            if (bundle2sheetState[i] != bundle2lastState[i])
            {
                bundle2sheet[i] = sheetFinder.integerSheet(i, top.add, false, useGlobalProblem);
                bundle2sheetState[i] = bundle2lastState[i];
            }

            auto sheet = bundle2sheet[i];

            if (sheet.empty())
                continue;

            double deltaObj = 0.0;
            for (auto& kv : sheet)
            {
                int j = kv.first;
                EH aj = *_decomp.bundle2arcs[j].begin();
                double d = kv.second;

                double xopt = mcMeshProps().get<ARC_DBL_LENGTH>(aj) * scaling;
                int xcurr = mcMeshProps().get<ARC_INT_LENGTH>(aj);

                deltaObj += d * (d + 2 * (xcurr - xopt)) * _decomp.bundle2arcs[j].size();
            }
            if (deltaObj < -1e-6)
            {
                DLOG(INFO) << "Found a sheet thats worth " << (top.add ? "adding." : "subtracting.") << " currentObj "
                           << currentObj << " newObj " << currentObj + deltaObj;
                // Execute update
                currentObj = currentObj + deltaObj;
                improvement = true;

                for (auto& kv : sheet)
                    for (EH a : _decomp.bundle2arcs[kv.first])
                        mcMeshProps().ref<ARC_INT_LENGTH>(a) += std::round(kv.second);

                for (int i2 : affectedBundles(sheet, bundle2sheet, sheetFinder.dynamicConstraints()))
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

    return currentObj;
}

double ISPQuantizer::makeFeasible(BaseLPSolver& sheetFinder, double scaling)
{
    double currentObj = obj(scaling);

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
            double deltaObj = 0.0;
            for (auto& kv : bestSheet)
            {
                int j = kv.first;
                EH aj = *_decomp.bundle2arcs[j].begin();
                double d = kv.second;

                double xopt = mcMeshProps().get<ARC_DBL_LENGTH>(aj) * scaling;
                int xcurr = mcMeshProps().get<ARC_INT_LENGTH>(aj);

                deltaObj += d * (d + 2 * (xcurr - xopt)) * _decomp.bundle2arcs[j].size();
            }
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
                    for (auto& kv : sheet)
                        for (EH a : _decomp.bundle2arcs[kv.first])
                            mcMeshProps().ref<ARC_INT_LENGTH>(a) += kv.second;
                    int nUnfulfilled = numViolatedConstraints(sheetFinder.dynamicConstraints());

                    // revert test
                    for (auto& kv : sheet)
                        for (EH a : _decomp.bundle2arcs[kv.first])
                            mcMeshProps().ref<ARC_INT_LENGTH>(a) -= kv.second;

                    if (nUnfulfilled < nUnfulfilledPre)
                    {
                        feasibilityImprovement = nUnfulfilledPre - nUnfulfilled;
                        break;
                    }
                }
                if (feasibilityImprovement <= 0)
                    continue;

                double deltaObj = 0.0;
                for (auto& kv : sheet)
                {
                    int j = kv.first;
                    EH aj = *_decomp.bundle2arcs[j].begin();
                    double d = kv.second;

                    double xopt = mcMeshProps().get<ARC_DBL_LENGTH>(aj) * scaling;
                    int xcurr = mcMeshProps().get<ARC_INT_LENGTH>(aj);

                    deltaObj += d * (d + 2 * (xcurr - xopt)) * _decomp.bundle2arcs[j].size();
                }
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
            for (auto& kv : bundle2sheet[bestBundle])
                for (EH a : _decomp.bundle2arcs[kv.first])
                    mcMeshProps().ref<ARC_INT_LENGTH>(a) += std::round(kv.second);

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

vector<vector<pair<int, EH>>> ISPQuantizer::violatedSimpleConstraints(double varLowerBound,
                                                                      vector<CriticalLink>& criticalLinks)
{
    auto& mcMesh = mcMeshProps().mesh();

    vector<vector<pair<int, EH>>> newConstraints;
    bool foundVarBelowBounds = containsMatching(mcMesh.edges(),
                                                [this, varLowerBound](const EH& a)
                                                { return mcMeshProps().get<ARC_INT_LENGTH>(a) < varLowerBound; });
    if (foundVarBelowBounds)
    {
        for (EH a : mcMesh.edges())
        {
            newConstraints.emplace_back();
            newConstraints.back().push_back({1, a});
            newConstraints.back().push_back({varLowerBound, EH()});
        }
    }

    bool collapsedLink
        = containsMatching(criticalLinks,
                           [this](const CriticalLink& link)
                           {
                               if (link.pathHas.empty())
                                   return false;
                               int sum = 0;
                               for (HEH ha : link.pathHas)
                                   sum += mcMeshProps().get<ARC_INT_LENGTH>(mcMeshProps().mesh().edge_handle(ha));
                               return sum <= 0;
                           });
    if (collapsedLink)
    {
        for (auto& criticalLink : criticalLinks)
        {
            if (criticalLink.pathHas.empty())
                continue;
            newConstraints.emplace_back();
            for (HEH ha : criticalLink.pathHas)
                newConstraints.back().push_back({1, mcMesh.edge_handle(ha)});
        }
    }
    bool negBlock
        = containsMatching(mcMesh.cells(),
                           [this](const CH& b)
                           {
                               for (UVWDir dir : {UVWDir::POS_V_POS_W, UVWDir::POS_U_POS_W, UVWDir::POS_U_POS_V})
                               {
                                   int length = 0;
                                   auto& arcs = mcMeshProps().ref<BLOCK_EDGE_ARCS>(b).at(dir);
                                   for (EH a : arcs)
                                       length += mcMeshProps().get<ARC_INT_LENGTH>(a);
                                   if (length < 0)
                                       return true;
                               }
                               return false;
                           });
    if (negBlock)
    {
        for (CH b : mcMesh.cells())
        {
            for (UVWDir dir : {UVWDir::POS_V_POS_W, UVWDir::POS_U_POS_W, UVWDir::POS_U_POS_V})
            {
                newConstraints.emplace_back();
                auto& nonZeroSum = newConstraints.back();
                auto& arcs = mcMeshProps().ref<BLOCK_EDGE_ARCS>(b).at(dir);
                for (EH a : arcs)
                    nonZeroSum.push_back({1, a});
                // Abuse this as marker that constraint is >= 0
                nonZeroSum.push_back({0, EH()});
            }
        }
    }
    DLOG(INFO) << "Found var below bounds? " << foundVarBelowBounds;
    DLOG(INFO) << "Found critical link? " << collapsedLink;
    DLOG(INFO) << "Found negative block? " << negBlock;

    return newConstraints;
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

double ISPQuantizer::obj(double scaling) const
{
    double obj = 0.0;
    for (EH a : mcMeshProps().mesh().edges())
    {
        double xopt = mcMeshProps().get<ARC_DBL_LENGTH>(a) * scaling;
        int xcurr = mcMeshProps().get<ARC_INT_LENGTH>(a);
        obj += (xopt - xcurr) * (xopt - xcurr);
    }
    return obj;
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
        if (containsMatching(coll,
                             [&bundleSet, this](const pair<const int, EH>& sign2a) {
                                 return sign2a.second.is_valid()
                                        && bundleSet.count(_decomp.arc2bundle.at(sign2a.second)) != 0;
                             }))
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

} // namespace qgp3d
