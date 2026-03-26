#include "QGP3D/StructurePreserver.hpp"

#include "QGP3D/ConstraintExtractor.hpp"

namespace qgp3d
{

StructurePreserver::StructurePreserver(TetMeshProps& meshProps, int maxSingIdx)
    : TetMeshNavigator(meshProps), MCMeshNavigator(meshProps), _maxSingIdx(maxSingIdx)
{
    fixedSingularities = !meshProps.isAllocated<ALGO_VARIANT>() || (meshProps.get<ALGO_VARIANT>() % 2);
    paddedSingularities = !meshProps.isAllocated<ALGO_VARIANT>() || (meshProps.get<ALGO_VARIANT>() / 4) % 2;

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
                        fixedSingularities,
                        true);

    isCriticalArcFromRegion = isCriticalArc;
    if (paddedSingularities)
        for (EH a : mcMeshProps().mesh().edges())
            if (mcMeshProps().get<IS_SINGULAR>(a))
                isCriticalArcFromRegion[a.idx()] = true;

    vector<bool> trash1_, trash2_;
    map<FH, int> trash3_;
    map<VH, vector<int>> trash4_, trash5_;
    getCriticalEntities(trash1_,
                        isSingularArc,
                        trash2_,
                        singularEntities,
                        a2singularLinkIdx,
                        trash3_,
                        trash4_,
                        trash5_,
                        false,
                        true,
                        false);
}

long StructurePreserver::numHexesInQuantization() const
{
    const MCMesh& mc = mcMeshProps().mesh();
    long nHexes = 0;
    for (CH b : mc.cells())
    {
        long nBlockHexes = 1;
        for (UVWDir dir : {UVWDir::NEG_U_NEG_V, UVWDir::NEG_U_NEG_W, UVWDir::NEG_V_NEG_W})
        {
            long arcLen = 0;
            // WARNING this way of using range based for loop (chaining 2 function calls)
            // depends on the first calls not returning rvalues!
            for (EH a : mcMeshProps().ref<BLOCK_EDGE_ARCS>(b).at(dir))
                arcLen += mcMeshProps().get<ARC_INT_LENGTH>(a);
            nBlockHexes *= arcLen;
        }
        nHexes += nBlockHexes;
    }
    return nHexes;
}

StructurePreserver::RetCode
StructurePreserver::violatedSimpleConstraints(double varLowerBound, vector<vector<pair<int, EH>>>& nonZeroSumArcs) const
{
    auto& mcMesh = mcMeshProps().mesh();

    bool foundVarBelowBounds = varLowerBound > INT_MIN
                               && containsMatching(mcMesh.edges(),
                                                   [this, varLowerBound](const EH& a)
                                                   { return mcMeshProps().get<ARC_INT_LENGTH>(a) < varLowerBound; });
    bool collapsedLink
        = containsMatching(criticalEntities,
                           [this](const CriticalEntity& link)
                           {
                               if (link.pathHas.empty())
                                   return false;
                               int sum = 0;
                               for (HEH ha : link.pathHas)
                                   sum += mcMeshProps().get<ARC_INT_LENGTH>(mcMeshProps().mesh().edge_handle(ha));
                               return sum <= 0;
                           });

    bool collapsedRegion = containsMatching(criticalEntities,
                                            [&, this](const CriticalEntity& region)
                                            {
                                                if (region.regionPs.empty())
                                                    return false;
                                                for (FH p : region.regionPs)
                                                    for (EH a : mcMesh.face_edges(p))
                                                        if (!isZeroArc(a))
                                                            return false;
                                                return true;
                                            });

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

    if (foundVarBelowBounds || collapsedLink || collapsedRegion || negBlock)
    {
        if (varLowerBound > INT_MIN)
            for (EH a : mcMesh.edges())
            {
                nonZeroSumArcs.emplace_back();
                nonZeroSumArcs.back().push_back({1, a});
                nonZeroSumArcs.back().push_back({varLowerBound, EH()});
            }
        for (auto& criticalLink : criticalEntities)
        {
            if (criticalLink.pathHas.empty())
                continue;
            nonZeroSumArcs.emplace_back();
            for (HEH ha : criticalLink.pathHas)
                nonZeroSumArcs.back().push_back({1, mcMesh.edge_handle(ha)});
        }
        for (auto& region : criticalEntities)
        {
            if (region.regionPs.empty())
                continue;
            nonZeroSumArcs.emplace_back();
            set<EH> as;
            for (FH p : region.regionPs)
                for (EH a : mcMesh.face_edges(p))
                    as.insert(a);
            for (EH a : as)
                nonZeroSumArcs.back().push_back({1, a});
        }
        for (CH b : mcMesh.cells())
        {
            for (UVWDir dir : {UVWDir::POS_V_POS_W, UVWDir::POS_U_POS_W, UVWDir::POS_U_POS_V})
            {
                nonZeroSumArcs.emplace_back();
                auto& nonZeroSum = nonZeroSumArcs.back();
                auto& arcs = mcMeshProps().ref<BLOCK_EDGE_ARCS>(b).at(dir);
                for (EH a : arcs)
                    nonZeroSum.push_back({1, a});
                // Abuse this as marker that constraint is >= 0
                nonZeroSum.push_back({0, EH()});
            }
        }
    }

    DLOG(INFO) << "Found var below bounds? " << foundVarBelowBounds;
    DLOG(INFO) << "Found collapsed critical link? " << collapsedLink;
    DLOG(INFO) << "Found collapsed critical region? " << collapsedRegion;
    DLOG(INFO) << "Found negative block? " << negBlock;

    return SUCCESS;
}

StructurePreserver::RetCode
StructurePreserver::violatedStructuralConstraints(vector<vector<pair<int, EH>>>& nonZeroSumArcs)
{
    auto& mcMesh = mcMeshProps().mesh();

    nonZeroSumArcs.clear();
    vector<vector<pair<int, EH>>> failsafeNonZeroSumArcs;

    if (fixedSingularities)
    {
        for (auto& criticalEntity : criticalEntities)
        {
            vector<vector<pair<int, EH>>> test;

            vector<set<HEH>> sector2has;
            set<EH> allowedPreBranchArcs;

            WeaklyMonotonousPath pStart;
            pStart.branchedOff = false;
            pStart.length = 0;
            pStart.path = {};
            pStart.delta = Vec3i(0, 0, 0);
            pStart.walkedDirs = UVWDir::NONE;
            pStart.transCurrent = Transition();

            if (criticalEntity.dim == 0)
            {
                // Each halfedge leaving critical node is a single possible direction
                for (HEH ha : mcMesh.outgoing_halfedges(criticalEntity.nFrom))
                    sector2has.push_back({ha});

                pStart.n = criticalEntity.nFrom;
                pStart.bRefCurrent = *mcMesh.vc_iter(criticalEntity.nFrom);
                pStart.dirs1 = UVWDir::NONE; // No directions are within critical entity for nodes
                pStart.monotonousDirs = ~pStart.dirs1;

                pStart.deltaMin = Vec3i(0, 0, 0);
                pStart.deltaMax = Vec3i(0, 0, 0);
            }
            else if (criticalEntity.dim == 1)
            {
                HEH criticalStartHa = criticalEntity.pathHas.front();
                allowedPreBranchArcs.insert(mcMesh.edge_handle(criticalStartHa));

                // Categorize all vertical has by direction
                for (HFH criticalStartHp : mcMesh.halfedge_halffaces(criticalStartHa))
                {
                    sector2has.emplace_back();

                    auto& has = sector2has.back();

                    CH bRef = mcMesh.incident_cell(criticalStartHp);
                    if (!bRef.is_valid())
                        bRef = mcMesh.incident_cell(mcMesh.opposite_halfface_handle(criticalStartHp));

                    UVWDir dirStart = halfarcDirInBlock(criticalStartHa, bRef);

                    HEH haOrth = mcMesh.prev_halfedge_in_halfface(criticalStartHa, criticalStartHp);
                    while (halfarcDirInBlock(haOrth, bRef) == dirStart)
                    {
                        allowedPreBranchArcs.insert(mcMesh.edge_handle(haOrth));
                        haOrth = mcMesh.prev_halfedge_in_halfface(haOrth, criticalStartHp);
                    }
                    haOrth = mcMesh.opposite_halfedge_handle(haOrth);
                    has.insert(haOrth);

                    HEH haCurr = criticalStartHa;
                    HFH hpCurr = criticalStartHp;

                    do
                    {
                        bRef = mcMesh.incident_cell(hpCurr);
                        if (!bRef.is_valid())
                            bRef = mcMesh.incident_cell(mcMesh.opposite_halfface_handle(hpCurr));

                        dirStart = halfarcDirInBlock(haCurr, bRef);

                        while (halfarcDirInBlock(mcMesh.next_halfedge_in_halfface(haCurr, hpCurr), bRef) == dirStart)
                        {
                            haCurr = mcMesh.next_halfedge_in_halfface(haCurr, hpCurr);
                            allowedPreBranchArcs.insert(mcMesh.edge_handle(haCurr));
                        }
                        haOrth = mcMesh.next_halfedge_in_halfface(haCurr, hpCurr);
                        if (has.count(haOrth))
                            break;
                        has.insert(haOrth);

                        HEH haOrthOpp = mcMesh.opposite_halfedge_handle(haOrth);
                        hpCurr
                            = findMatching(mcMesh.halfedge_halffaces(haOrthOpp),
                                           [&](const HFH& hpNext)
                                           {
                                               return mcMesh.face_handle(hpNext) != mcMesh.face_handle(hpCurr)
                                                      && contains(criticalEntity.pathHas,
                                                                  mcMesh.next_halfedge_in_halfface(haOrthOpp, hpNext));
                                           });
                        if (hpCurr.is_valid())
                            haCurr = mcMesh.next_halfedge_in_halfface(haOrthOpp, hpCurr);
                    } while (hpCurr.is_valid() && haCurr != criticalStartHa);
                }

                pStart.n = criticalEntity.nFrom;
                pStart.bRefCurrent = *mcMesh.hec_iter(criticalStartHa);
                UVWDir dirEntity = halfarcDirInBlock(criticalStartHa, pStart.bRefCurrent);
                pStart.dirs1 = dirEntity | -dirEntity;
                pStart.monotonousDirs = ~pStart.dirs1;

                int length = 0;
                for (HEH ha : criticalEntity.pathHas)
                    length += (_considerDoubleIntLength ? 2 : 1)
                              * mcMeshProps().get<ARC_INT_LENGTH>(mcMeshProps().mesh().edge_handle(ha));
                pStart.deltaMin = (isNeg(dirEntity) ? length * toVec(dirEntity) : Vec3i(0, 0, 0));
                pStart.deltaMax = (isNeg(dirEntity) ? Vec3i(0, 0, 0) : length * toVec(dirEntity));
                if (criticalEntity.nFrom == criticalEntity.nTo && dirEntity != UVWDir::NONE)
                {
                    pStart.deltaMin[toCoord(dirEntity)] = INT_MIN;
                    pStart.deltaMax[toCoord(dirEntity)] = INT_MAX;
                }
            }
            else if (criticalEntity.dim == 2)
            {
                continue;
            }

            for (int sector = 0; sector < (int)sector2has.size(); sector++)
            {
                traceExhaustPaths(
                    pStart, sector2has[sector], allowedPreBranchArcs, nonZeroSumArcs, failsafeNonZeroSumArcs, true);
            }
        }
    }
    else
    {
        precompute();

        // Cyclic singular links may not collapse
        {
            vector<vector<pair<int, EH>>> forcedNonZeroSum;
            for (auto& singularLink : singularEntities)
                if (singularLink.dim == 1 && singularLink.nFrom == singularLink.nTo)
                {
                    int length = 0;
                    for (HEH ha : singularLink.pathHas)
                        length += mcMeshProps().get<ARC_INT_LENGTH>(mcMesh.edge_handle(ha));
                    if (length == 0)
                    {
                        failsafeNonZeroSumArcs.emplace_back();
                        auto& constraint = failsafeNonZeroSumArcs.back();
                        for (HEH ha : singularLink.pathHas)
                            constraint.push_back({1, mcMesh.edge_handle(ha)});
                        constraint.push_back({1, EH()});
                        nonZeroSumArcs.emplace_back();
                        nonZeroSumArcs.back() = constraint;
                        DLOG(INFO) << "Pumping cyclic singular link " << singularLink.id;
                    }
                }
        }

        // Different critical entities may not merge
        {
            set<VH> nsSkippable;
            set<EH> asSkippable;

            // For each critical link s1 find paths connecting s1 to other critical links s2 or surface patches p2
            // Then check for overlaps between these, accumulating the quantized edge lengths along the path as deltas.
            for (auto& criticalEntity : criticalEntities)
            {
                if (criticalEntity.dim == 0)
                {
                    IntegerGridGraphNode c0(criticalEntity.nFrom);
                    CH bRef = *mcMesh.vc_iter(c0.n);

                    set<IntegerGridGraphNode> problematicNodes;
                    for (auto& crit2 : criticalEntities)
                    {
                        if (crit2.dim == 0 && crit2.nFrom != c0.n)
                            problematicNodes.emplace(crit2.nFrom);
                        else if (crit2.dim == 1 && crit2.nFrom != c0.n && crit2.nTo != c0.n)
                            for (HEH ha : crit2.pathHas)
                            {
                                problematicNodes.emplace(mcMesh.from_vertex_handle(ha));
                                problematicNodes.emplace(mcMesh.to_vertex_handle(ha));
                                problematicNodes.emplace(mcMesh.edge_handle(ha));
                            }
                        else if (crit2.dim == 2
                                 && !containsMatching(crit2.regionPs,
                                                      [&, this](const FH& p)
                                                      { return contains(mcMesh.face_vertices(p), c0.n); }))
                            for (FH p : crit2.regionPs)
                            {
                                problematicNodes.emplace(p);
                                for (HEH ha : mcMesh.halfface_halfedges(mcMesh.halfface_handle(p, 0)))
                                {
                                    problematicNodes.emplace(mcMesh.from_vertex_handle(ha));
                                    problematicNodes.emplace(mcMesh.to_vertex_handle(ha));
                                    problematicNodes.emplace(mcMesh.edge_handle(ha));
                                }
                            }
                    }
                    set<IntegerGridGraphNode> potentiallyUnproblematicNodes;
                    potentiallyUnproblematicNodes.insert(c0);
                    for (auto& crit2 : criticalEntities)
                    {
                        if (crit2.dim == 1 && (crit2.nFrom == c0.n || crit2.nTo == c0.n))
                            for (HEH ha : crit2.pathHas)
                            {
                                if (!problematicNodes.count(mcMesh.from_vertex_handle(ha)))
                                    potentiallyUnproblematicNodes.emplace(mcMesh.from_vertex_handle(ha));
                                if (!problematicNodes.count(mcMesh.to_vertex_handle(ha)))
                                    potentiallyUnproblematicNodes.emplace(mcMesh.to_vertex_handle(ha));
                                if (!problematicNodes.count(mcMesh.edge_handle(ha)))
                                    potentiallyUnproblematicNodes.emplace(mcMesh.edge_handle(ha));
                            }
                        else if (crit2.dim == 2
                                 && containsMatching(crit2.regionPs,
                                                     [&, this](const FH& p)
                                                     { return contains(mcMesh.face_vertices(p), c0.n); }))
                            for (FH p : crit2.regionPs)
                            {
                                if (!problematicNodes.count(p))
                                    potentiallyUnproblematicNodes.emplace(p);
                                for (HEH ha : mcMesh.halfface_halfedges(mcMesh.halfface_handle(p, 0)))
                                {
                                    if (!problematicNodes.count(mcMesh.from_vertex_handle(ha)))
                                        potentiallyUnproblematicNodes.emplace(mcMesh.from_vertex_handle(ha));
                                    if (!problematicNodes.count(mcMesh.to_vertex_handle(ha)))
                                        potentiallyUnproblematicNodes.emplace(mcMesh.to_vertex_handle(ha));
                                    if (!problematicNodes.count(mcMesh.edge_handle(ha)))
                                        potentiallyUnproblematicNodes.emplace(mcMesh.edge_handle(ha));
                                }
                            }
                    }

                    DLOG(INFO) << "Tracing from critical node " << criticalEntity.id << " (" << criticalEntity.nFrom
                               << " / " << mcMeshProps().get<NODE_MESH_VERTEX>(criticalEntity.nFrom) << ")";
                    set<IntegerGridGraphNode> unproblematicNodes;
                    traceUnproblematicPaths(c0,
                                            bRef,
                                            _b2n2igm[bRef.idx()].at(c0.n),
                                            _b2n2igm[bRef.idx()].at(c0.n),
                                            potentiallyUnproblematicNodes,
                                            unproblematicNodes);

                    if (traceExhaustPaths(c0,
                                          bRef,
                                          _b2n2igm[bRef.idx()].at(c0.n),
                                          _b2n2igm[bRef.idx()].at(c0.n),
                                          unproblematicNodes,
                                          nonZeroSumArcs,
                                          failsafeNonZeroSumArcs))
                        DLOG(INFO) << "Separating feature node " << criticalEntity.id;
                    nsSkippable.insert(c0.n);
                }
                else if (criticalEntity.dim == 1)
                {
                    set<IntegerGridGraphNode> problematicNodes;
                    for (auto& crit2 : criticalEntities)
                    {
                        if (crit2.dim == 0 && crit2.nFrom != criticalEntity.nFrom && crit2.nFrom != criticalEntity.nTo)
                            problematicNodes.emplace(crit2.nFrom);
                        else if (crit2.dim == 1 && crit2.nFrom != criticalEntity.nFrom
                                 && crit2.nTo != criticalEntity.nFrom && crit2.nFrom != criticalEntity.nTo
                                 && crit2.nTo != criticalEntity.nTo)
                            for (HEH ha : crit2.pathHas)
                            {
                                problematicNodes.emplace(mcMesh.from_vertex_handle(ha));
                                problematicNodes.emplace(mcMesh.to_vertex_handle(ha));
                                problematicNodes.emplace(mcMesh.edge_handle(ha));
                            }
                        else if (crit2.dim == 2
                                 && !containsMatching(crit2.regionPs,
                                                      [&, this](const FH& p)
                                                      {
                                                          return containsSomeOf(
                                                              mcMesh.face_vertices(p),
                                                              std::set<VH>({criticalEntity.nFrom, criticalEntity.nTo}));
                                                      }))
                            for (FH p : crit2.regionPs)
                            {
                                problematicNodes.emplace(p);
                                for (HEH ha : mcMesh.halfface_halfedges(mcMesh.halfface_handle(p, 0)))
                                {
                                    problematicNodes.emplace(mcMesh.from_vertex_handle(ha));
                                    problematicNodes.emplace(mcMesh.to_vertex_handle(ha));
                                    problematicNodes.emplace(mcMesh.edge_handle(ha));
                                }
                            }
                    }
                    set<IntegerGridGraphNode> potentiallyUnproblematicNodes;
                    for (auto& crit2 : criticalEntities)
                    {
                        if (crit2.dim == 1
                            && (crit2.nFrom == criticalEntity.nFrom || crit2.nTo == criticalEntity.nFrom
                                || crit2.nFrom == criticalEntity.nTo || crit2.nTo == criticalEntity.nTo))
                            for (HEH ha : crit2.pathHas)
                            {
                                if (!problematicNodes.count(mcMesh.from_vertex_handle(ha)))
                                    potentiallyUnproblematicNodes.emplace(mcMesh.from_vertex_handle(ha));
                                if (!problematicNodes.count(mcMesh.to_vertex_handle(ha)))
                                    potentiallyUnproblematicNodes.emplace(mcMesh.to_vertex_handle(ha));
                                if (!problematicNodes.count(mcMesh.edge_handle(ha)))
                                    potentiallyUnproblematicNodes.emplace(mcMesh.edge_handle(ha));
                            }
                        else if (crit2.dim == 2
                                 && containsMatching(crit2.regionPs,
                                                     [&, this](const FH& p)
                                                     {
                                                         return containsSomeOf(
                                                             mcMesh.face_vertices(p),
                                                             std::set<VH>({criticalEntity.nFrom, criticalEntity.nTo}));
                                                     }))
                            for (FH p : crit2.regionPs)
                            {
                                if (!problematicNodes.count(p))
                                    potentiallyUnproblematicNodes.emplace(p);
                                for (HEH ha : mcMesh.halfface_halfedges(mcMesh.halfface_handle(p, 0)))
                                {
                                    if (!problematicNodes.count(mcMesh.from_vertex_handle(ha)))
                                        potentiallyUnproblematicNodes.emplace(mcMesh.from_vertex_handle(ha));
                                    if (!problematicNodes.count(mcMesh.to_vertex_handle(ha)))
                                        potentiallyUnproblematicNodes.emplace(mcMesh.to_vertex_handle(ha));
                                    if (!problematicNodes.count(mcMesh.edge_handle(ha)))
                                        potentiallyUnproblematicNodes.emplace(mcMesh.edge_handle(ha));
                                }
                            }
                    }

                    for (HEH ha : criticalEntity.pathHas)
                    {
                        EH a = mcMesh.edge_handle(ha);
                        int l = (_considerDoubleIntLength ? 2 : 1) * mcMeshProps().get<ARC_INT_LENGTH>(a);
                        CH bRef = *mcMesh.hec_iter(ha);

                        auto& n2IGM = _b2n2igm[bRef.idx()];

                        for (VH n : mcMesh.halfedge_vertices(ha))
                        {
                            if (!nsSkippable.count(n))
                            {
                                DLOG(INFO) << "Tracing from critical link " << criticalEntity.id << " arc " << a
                                           << " startnode " << n;

                                set<IntegerGridGraphNode> unproblematicNodes;
                                traceUnproblematicPaths({n},
                                                        bRef,
                                                        n2IGM.at(n),
                                                        n2IGM.at(n),
                                                        potentiallyUnproblematicNodes,
                                                        unproblematicNodes);

                                if (traceExhaustPaths({n},
                                                      bRef,
                                                      n2IGM.at(n),
                                                      n2IGM.at(n),
                                                      unproblematicNodes,
                                                      nonZeroSumArcs,
                                                      failsafeNonZeroSumArcs))
                                    DLOG(INFO) << "Separating feature curve " << criticalEntity.id;
                                nsSkippable.insert(n);
                                for (HEH haSkip : criticalEntity.pathHas)
                                    for (VH nSkip : mcMesh.halfedge_vertices(haSkip))
                                    {
                                        if (unproblematicNodes.count(nSkip))
                                            nsSkippable.insert(nSkip);
                                    }
                            }
                        }

                        if (l >= 2 && !asSkippable.count(a))
                        {
                            Vec3i igmMin = _b2a2igmMin[bRef.idx()].at(a);
                            Vec3i igmMax = _b2a2igmMax[bRef.idx()].at(a);
                            Vec3i delta = (igmMax - igmMin) / l;
                            igmMin += delta;
                            igmMax -= delta;

                            set<IntegerGridGraphNode> unproblematicNodes;
                            traceUnproblematicPaths(
                                {a}, bRef, igmMin, igmMax, potentiallyUnproblematicNodes, unproblematicNodes);

                            DLOG(INFO) << "Tracing from critical link " << criticalEntity.id << " arc " << a;
                            if (traceExhaustPaths({a},
                                                  bRef,
                                                  igmMin,
                                                  igmMax,
                                                  unproblematicNodes,
                                                  nonZeroSumArcs,
                                                  failsafeNonZeroSumArcs))
                                DLOG(INFO) << "Separating feature curve " << criticalEntity.id;
                            asSkippable.insert(a);
                        }
                    }
                }
                else // criticalEntity.dim == 2
                {
                    struct TempSwapper
                    {
                        vector<bool>&_a, _b;
                        TempSwapper(vector<bool>& a, vector<bool>& b) : _a(a), _b(b)
                        {
                            std::swap(_a, _b);
                        }
                        ~TempSwapper()
                        {
                            std::swap(_a, _b);
                        }
                    };
                    TempSwapper swp(isCriticalArc, isCriticalArcFromRegion);

                    set<VH> nsRegion;
                    set<EH> asRegion;
                    for (FH p : criticalEntity.regionPs)
                    {
                        for (VH n : mcMesh.face_vertices(p))
                            nsRegion.insert(n);
                        for (EH a : mcMesh.face_edges(p))
                            asRegion.insert(a);
                    }
                    set<IntegerGridGraphNode> problematicNodes;
                    for (auto& crit2 : criticalEntities)
                    {
                        if (crit2.dim == 0 && !nsRegion.count(crit2.nFrom))
                            problematicNodes.emplace(crit2.nFrom);
                        else if (crit2.dim == 1 && !nsRegion.count(crit2.nFrom) && !nsRegion.count(crit2.nTo))
                            for (HEH ha : crit2.pathHas)
                            {
                                problematicNodes.emplace(mcMesh.from_vertex_handle(ha));
                                problematicNodes.emplace(mcMesh.to_vertex_handle(ha));
                                problematicNodes.emplace(mcMesh.edge_handle(ha));
                            }
                        else if (crit2.dim == 2
                                 && !containsMatching(crit2.regionPs,
                                                      [&, this](const FH& p)
                                                      { return containsSomeOf(mcMesh.face_vertices(p), nsRegion); }))
                            for (FH p : crit2.regionPs)
                            {
                                problematicNodes.emplace(p);
                                for (HEH ha : mcMesh.halfface_halfedges(mcMesh.halfface_handle(p, 0)))
                                {
                                    problematicNodes.emplace(mcMesh.from_vertex_handle(ha));
                                    problematicNodes.emplace(mcMesh.to_vertex_handle(ha));
                                    problematicNodes.emplace(mcMesh.edge_handle(ha));
                                }
                            }
                    }
                    if (paddedSingularities)
                        for (auto& crit2 : singularEntities)
                        {
                            if (crit2.dim == 1 && !nsRegion.count(crit2.nFrom) && !nsRegion.count(crit2.nTo))
                                for (HEH ha : crit2.pathHas)
                                {
                                    problematicNodes.emplace(mcMesh.from_vertex_handle(ha));
                                    problematicNodes.emplace(mcMesh.to_vertex_handle(ha));
                                    problematicNodes.emplace(mcMesh.edge_handle(ha));
                                }
                        }

                    set<IntegerGridGraphNode> potentiallyUnproblematicNodes;
                    for (auto& crit2 : criticalEntities)
                    {
                        if (crit2.dim == 1 && (nsRegion.count(crit2.nFrom) || nsRegion.count(crit2.nTo)))
                            for (HEH ha : crit2.pathHas)
                            {
                                if (!problematicNodes.count(mcMesh.from_vertex_handle(ha)))
                                    potentiallyUnproblematicNodes.emplace(mcMesh.from_vertex_handle(ha));
                                if (!problematicNodes.count(mcMesh.to_vertex_handle(ha)))
                                    potentiallyUnproblematicNodes.emplace(mcMesh.to_vertex_handle(ha));
                                if (!problematicNodes.count(mcMesh.edge_handle(ha)))
                                    potentiallyUnproblematicNodes.emplace(mcMesh.edge_handle(ha));
                            }
                        else if (crit2.dim == 2
                                 && containsMatching(crit2.regionPs,
                                                     [&, this](const FH& p)
                                                     { return containsSomeOf(mcMesh.face_vertices(p), nsRegion); }))
                            for (FH p : crit2.regionPs)
                            {
                                if (!problematicNodes.count(p))
                                    potentiallyUnproblematicNodes.emplace(p);
                                for (HEH ha : mcMesh.halfface_halfedges(mcMesh.halfface_handle(p, 0)))
                                {
                                    if (!problematicNodes.count(mcMesh.from_vertex_handle(ha)))
                                        potentiallyUnproblematicNodes.emplace(mcMesh.from_vertex_handle(ha));
                                    if (!problematicNodes.count(mcMesh.to_vertex_handle(ha)))
                                        potentiallyUnproblematicNodes.emplace(mcMesh.to_vertex_handle(ha));
                                    if (!problematicNodes.count(mcMesh.edge_handle(ha)))
                                        potentiallyUnproblematicNodes.emplace(mcMesh.edge_handle(ha));
                                }
                            }
                    }
                    if (paddedSingularities)
                        for (auto& crit2 : singularEntities)
                            if (crit2.dim == 1 && (nsRegion.count(crit2.nFrom) || nsRegion.count(crit2.nTo)))
                                for (HEH ha : crit2.pathHas)
                                {
                                    if (!problematicNodes.count(mcMesh.from_vertex_handle(ha)))
                                        potentiallyUnproblematicNodes.emplace(mcMesh.from_vertex_handle(ha));
                                    if (!problematicNodes.count(mcMesh.to_vertex_handle(ha)))
                                        potentiallyUnproblematicNodes.emplace(mcMesh.to_vertex_handle(ha));
                                    if (!problematicNodes.count(mcMesh.edge_handle(ha)))
                                        potentiallyUnproblematicNodes.emplace(mcMesh.edge_handle(ha));
                                }

                    for (FH p : criticalEntity.regionPs)
                    {
                        HFH hp = mcMesh.halfface_handle(p, 0);
                        if (mcMesh.is_boundary(hp))
                            hp = mcMesh.opposite_halfface_handle(hp);

                        CH bRef = mcMesh.incident_cell(hp);
                        auto& n2IGM = _b2n2igm[bRef.idx()];

                        for (EH a : mcMesh.face_edges(p))
                        {
                            int l = (_considerDoubleIntLength ? 2 : 1) * mcMeshProps().get<ARC_INT_LENGTH>(a);
                            for (VH n : mcMesh.edge_vertices(a))
                            {
                                if (nsSkippable.count(n))
                                    continue;

                                set<IntegerGridGraphNode> unproblematicNodes;
                                traceUnproblematicPaths({n},
                                                        bRef,
                                                        n2IGM.at(n),
                                                        n2IGM.at(n),
                                                        potentiallyUnproblematicNodes,
                                                        unproblematicNodes);

                                DLOG(INFO) << "Tracing from critical region " << criticalEntity.id << " node " << n;
                                if (traceExhaustPaths({n},
                                                      bRef,
                                                      n2IGM.at(n),
                                                      n2IGM.at(n),
                                                      unproblematicNodes,
                                                      nonZeroSumArcs,
                                                      failsafeNonZeroSumArcs))
                                    DLOG(INFO) << "Separating feature region " << criticalEntity.id;

                                nsSkippable.insert(n);
                                for (VH nSkip : nsRegion)
                                {
                                    if (unproblematicNodes.count(nSkip))
                                        nsSkippable.insert(nSkip);
                                }
                                if (!failsafeNonZeroSumArcs.empty())
                                    for (EH aSkip : asRegion)
                                    {
                                        if (unproblematicNodes.count(aSkip))
                                            asSkippable.insert(aSkip);
                                    }
                            }
                            if (!asSkippable.count(a) && l >= 2)
                            {
                                Vec3i igmMin = _b2a2igmMin[bRef.idx()].at(a);
                                Vec3i igmMax = _b2a2igmMax[bRef.idx()].at(a);
                                Vec3i delta = (igmMax - igmMin) / l;
                                igmMin += delta;
                                igmMax -= delta;
                                set<IntegerGridGraphNode> unproblematicNodes;
                                traceUnproblematicPaths(
                                    {a}, bRef, igmMin, igmMax, potentiallyUnproblematicNodes, unproblematicNodes);

                                DLOG(INFO) << "Tracing from critical region " << criticalEntity.id << " arc " << a;
                                if (traceExhaustPaths({a},
                                                      bRef,
                                                      igmMin,
                                                      igmMax,
                                                      unproblematicNodes,
                                                      nonZeroSumArcs,
                                                      failsafeNonZeroSumArcs))
                                    DLOG(INFO) << "Separating feature region " << criticalEntity.id;
                                asSkippable.insert(a);
                                for (VH nSkip : nsRegion)
                                {
                                    if (unproblematicNodes.count(nSkip))
                                        nsSkippable.insert(nSkip);
                                }
                                if (!failsafeNonZeroSumArcs.empty())
                                    for (EH aSkip : asRegion)
                                    {
                                        if (unproblematicNodes.count(aSkip))
                                            asSkippable.insert(aSkip);
                                    }
                            }
                        }

                        Vec3i igmMin = _b2p2igmMin[bRef.idx()].at(p);
                        Vec3i igmMax = _b2p2igmMax[bRef.idx()].at(p);
                        Vec3i delta = toVec(toDir(igmMax - igmMin));
                        if (dim(toDir(delta)) < 2)
                            continue;
                        igmMin += delta;
                        igmMax -= delta;
                        if (igmMin[0] <= igmMax[0] && igmMin[1] <= igmMax[0] && igmMin[2] <= igmMax[2])
                        {
                            set<IntegerGridGraphNode> unproblematicNodes;
                            unproblematicNodes.emplace(p);
                            DLOG(INFO) << "Tracing from critical region " << criticalEntity.id << " patch " << p;
                            std::swap(isCriticalArc, isCriticalArcFromRegion);
                            if (traceExhaustPaths({p},
                                                  bRef,
                                                  igmMin,
                                                  igmMax,
                                                  unproblematicNodes,
                                                  nonZeroSumArcs,
                                                  failsafeNonZeroSumArcs))
                                DLOG(INFO) << "Separating feature region " << criticalEntity.id;
                            std::swap(isCriticalArc, isCriticalArcFromRegion);
                        }
                    }
                }
            }
        }

        // Incontractible loops may not be zero-length
        {
            vector<vector<FH>> cutSurfaceComponents;
            vector<int> p2componentIdx;
            ConstraintExtractor(meshProps()).getCutSurfaces(cutSurfaceComponents, p2componentIdx);

            set<VH> nsCritical;
            for (auto& component : cutSurfaceComponents)
            {
                for (FH p : component)
                    for (VH n : mcMesh.face_vertices(p))
                        nsCritical.insert(n);
            }

            set<VH> nsVisited;
            for (VH nCritical : nsCritical)
            {
                if (nsVisited.count(nCritical))
                    continue;
                IntegerGridGraphNode c0(nCritical);
                CH bRef = *mcMesh.vc_iter(nCritical);
                DLOG(INFO) << "Tracing from cut surface node " << nCritical;

                if (findIncontractibleZeroLoops(c0,
                                                bRef,
                                                _b2n2igm[bRef.idx()].at(c0.n),
                                                _b2n2igm[bRef.idx()].at(c0.n),
                                                nsVisited,
                                                nonZeroSumArcs,
                                                failsafeNonZeroSumArcs))
                    DLOG(INFO) << "Pumping incontractible loop(s)";
            }
        }

        // singularities may not cluster to form out of bounds singularity indices
        if (failsafeNonZeroSumArcs.empty())
            findSingularityIndexViolationPaths(nonZeroSumArcs, failsafeNonZeroSumArcs);

        // Closed regions may not completely collapse (should already be caught by all the above)
        if (failsafeNonZeroSumArcs.empty())
            for (auto& region : criticalEntities)
            {
                if (region.regionPs.empty()
                    || containsMatching(region.regionPs, [&, this](const FH& p) { return !isZeroPatch(p); }))
                    continue;
                set<EH> as;
                for (FH p : region.regionPs)
                    if (containsMatching(mcMesh.face_edges(p), [&, this](const EH& a) { return !isZeroArc(a); }))
                    {
                        for (auto& kv : halfpatchHalfarcsByDir(mcMesh.halfface_handle(p, 0)))
                            if (!containsMatching(
                                    kv.second, [&, this](const HEH& ha) { return !isZeroArc(mcMesh.edge_handle(ha)); }))
                            {
                                for (HEH ha : kv.second)
                                {
                                    as.insert(mcMesh.edge_handle(ha));
                                }
                                break;
                            }
                    }
                failsafeNonZeroSumArcs.emplace_back();
                auto& constraint = failsafeNonZeroSumArcs.back();
                for (EH a : as)
                    constraint.push_back({1, a});
                nonZeroSumArcs.emplace_back();
                nonZeroSumArcs.back() = constraint;
                DLOG(INFO) << "Pumping region " << region.id;
            }
    }

    _failsafeStructuralConstraints.insert(
        _failsafeStructuralConstraints.end(), failsafeNonZeroSumArcs.begin(), failsafeNonZeroSumArcs.end());
    _allStructuralConstraints.insert(_allStructuralConstraints.end(), nonZeroSumArcs.begin(), nonZeroSumArcs.end());

    return SUCCESS;
}

void StructurePreserver::traceExhaustPaths(const WeaklyMonotonousPath& pStart,
                                           const set<HEH>& allowedBranchOffs,
                                           const set<EH>& allowedPreBranchArcs,
                                           vector<vector<pair<int, EH>>>& nonZeroSumArcs,
                                           vector<vector<pair<int, EH>>>& failsafeNonZeroSumArcs,
                                           bool returnOnViolation) const
{
    using PathQueue = std::
        priority_queue<WeaklyMonotonousPath, std::deque<WeaklyMonotonousPath>, GreaterMonotonousPathLengthCompare>;

    const MCMesh& mcMesh = mcMeshProps().mesh();

    vector<bool> nsVisited(mcMesh.n_vertices(), false);
    vector<bool> nsInitialized(mcMesh.n_vertices(), false);
    PathQueue pathQ;
    pathQ.push(pStart);
    while (!pathQ.empty())
    {
        auto pathCurrent = pathQ.top();
        pathQ.pop();

        if ((pathCurrent.branchedOff && nsVisited[pathCurrent.n.idx()])
            || (!pathCurrent.branchedOff && nsInitialized[pathCurrent.n.idx()]))
            continue;
        if (pathCurrent.branchedOff)
            nsVisited[pathCurrent.n.idx()] = true;
        else
            nsInitialized[pathCurrent.n.idx()] = true;

        map<CH, Transition> b2trans;
        map<HEH, vector<CH>> ha2bRef;
        determineNextHalfarcs(pathCurrent, b2trans, ha2bRef);

        // Check for separation violations by current path
        if ((pathCurrent.walkedDirs & pathCurrent.monotonousDirs) != UVWDir::NONE)
        {
            UVWDir violationDirs = UVWDir::NONE;
            // Check whether startinterval-criticalnode pairs overlap
            if (isCriticalNode[pathCurrent.n.idx()]
                && bboxOverlap(pathCurrent.deltaMin, pathCurrent.deltaMax, pathCurrent.delta, pathCurrent.delta))
                violationDirs = pathCurrent.walkedDirs & ~pathCurrent.dirs1;

            if (violationDirs == UVWDir::NONE)
            {
                // Check whether startinterval-criticallink pairs overlap
                for (const auto& kv : ha2bRef)
                {
                    HEH ha2 = kv.first;
                    EH a2 = mcMesh.edge_handle(ha2);
                    if (isCriticalArc[a2.idx()])
                    {
                        CH bRef2 = kv.second.front();
                        Transition trans2 = b2trans.at(bRef2);

                        UVWDir dir2 = trans2.invert().rotate(halfarcDirInBlock(ha2, bRef2));

                        UVWDir possibleViolationDir = pathCurrent.walkedDirs & ~(dir2 | -dir2) & ~pathCurrent.dirs1;
                        if ((possibleViolationDir & pathCurrent.monotonousDirs) == UVWDir::NONE)
                            continue;

                        // Measure just overlap with the arc, not the link
                        int length = (_considerDoubleIntLength ? 2 : 1) * mcMeshProps().get<ARC_INT_LENGTH>(a2);

                        if (length < 0)
                        {
                            length = -length;
                            dir2 = -dir2;
                        }
                        Vec3i deltaMin2 = pathCurrent.delta + (isNeg(dir2) ? length * toVec(dir2) : Vec3i(0, 0, 0));
                        Vec3i deltaMax2 = pathCurrent.delta + (isNeg(dir2) ? Vec3i(0, 0, 0) : length * toVec(dir2));

                        // Check for overlap
                        if (bboxOverlap(pathCurrent.deltaMin, pathCurrent.deltaMax, deltaMin2, deltaMax2))
                        {
                            assert(possibleViolationDir != UVWDir::NONE);
                            violationDirs = possibleViolationDir;
                            break;
                        }
                    }
                }
            }

            if (violationDirs == UVWDir::NONE)
            {
                // Check for overlaps between startinterval-surfacepatch pairs
                for (HFH hp2 : mcMesh.vertex_halffaces(pathCurrent.n))
                {
                    if (isCriticalPatch[mcMesh.face_handle(hp2).idx()])
                    {
                        HFH hp2opp = mcMesh.opposite_halfface_handle(hp2);
                        CH bRef2 = mcMesh.incident_cell(hp2opp);

                        if (b2trans.find(bRef2) == b2trans.end())
                            continue;
                        Transition trans2 = b2trans.at(bRef2);

                        UVWDir hpDirs = UVWDir::NONE;
                        for (HEH ha : mcMesh.halfface_halfedges(hp2opp))
                            hpDirs = hpDirs | halfarcDirInBlock(ha, bRef2);

                        UVWDir hpDirsLocal = trans2.invert().rotate(hpDirs);
                        assert(dim(hpDirsLocal) == 2);

                        UVWDir possibleViolationDir = pathCurrent.walkedDirs & ~hpDirsLocal & ~pathCurrent.dirs1;
                        if ((possibleViolationDir & pathCurrent.monotonousDirs) == UVWDir::NONE)
                            continue;

                        if (checkPatchOverlap(pathCurrent, hp2opp, trans2))
                        {
                            assert(possibleViolationDir != UVWDir::NONE);
                            violationDirs = possibleViolationDir;
                            break;
                        }
                    }
                }
            }

            if (violationDirs != UVWDir::NONE)
            {
                vector<EH> posSignArcs;
                vector<EH> negSignArcs;
                vector<EH> monotonousDirArcs;
                for (UVWDir dim1dirPos : {UVWDir::POS_U, UVWDir::POS_V, UVWDir::POS_W})
                {
                    UVWDir dim1dirAny = dim1dirPos | -dim1dirPos;
                    UVWDir intersection = violationDirs & dim1dirAny;
                    if (intersection != UVWDir::NONE)
                    {
                        if (intersection == dim1dirAny)
                        {
                            const auto& posArcs = pathCurrent.dir2walkedArcs.at(dim1dirPos);
                            const auto& negArcs = pathCurrent.dir2walkedArcs.at(-dim1dirPos);
                            assert(posArcs.size() > 0);
                            assert(negArcs.size() > 0);
                            double lengthPos = 0;
                            double lengthNeg = 0;
                            for (EH a : posArcs)
                                lengthPos += mcMeshProps().get<ARC_DBL_LENGTH>(a);
                            for (EH a : negArcs)
                                lengthNeg += mcMeshProps().get<ARC_DBL_LENGTH>(a);
                            if (lengthPos > lengthNeg)
                            {
                                posSignArcs.insert(posSignArcs.end(), posArcs.begin(), posArcs.end());
                                negSignArcs.insert(negSignArcs.end(), negArcs.begin(), negArcs.end());
                            }
                            else
                            {
                                posSignArcs.insert(posSignArcs.end(), negArcs.begin(), negArcs.end());
                                negSignArcs.insert(negSignArcs.end(), posArcs.begin(), posArcs.end());
                            }
                        }
                        else
                        {
                            auto& arcs = pathCurrent.dir2walkedArcs.at(intersection);
                            posSignArcs.insert(posSignArcs.end(), arcs.begin(), arcs.end());
                            monotonousDirArcs.insert(monotonousDirArcs.end(), arcs.begin(), arcs.end());
                        }
                    }
                }
                assert(posSignArcs.size() > 0);
                if (posSignArcs.empty())
                    continue;
                nonZeroSumArcs.emplace_back();
                auto& nonZeroSum = nonZeroSumArcs.back();
                for (EH a : posSignArcs)
                    nonZeroSum.push_back({1, a});
                for (EH a : negSignArcs)
                    nonZeroSum.push_back({-1, a});

                failsafeNonZeroSumArcs.emplace_back();
                auto& failsafeNonZeroSum = failsafeNonZeroSumArcs.back();
                for (EH a : monotonousDirArcs)
                    failsafeNonZeroSum.push_back({1, a});

                assert(pathCurrent.branchedOff);
                if (returnOnViolation)
                    return;
            }
        }

        for (const auto& kv : ha2bRef)
        {
            HEH ha = kv.first;
            VH nTo = mcMesh.to_vertex_handle(ha);
            if (pathCurrent.branchedOff && nsVisited[nTo.idx()])
                continue;

            for (CH bRef : kv.second)
            {
                auto& trans = b2trans.at(bRef);
                WeaklyMonotonousPath nextP = pathCurrent;
                nextP.bRefCurrent = bRef;
                nextP.transCurrent = trans;

                if (!checkP0Containment(nextP))
                    continue;

                nextP.n = nTo;

                UVWDir walkedDir = trans.invert().rotate(halfarcDirInBlock(ha, nextP.bRefCurrent));

                // Record branch-off from link, allow only limited set of edges
                if (!nextP.branchedOff)
                {
                    if (allowedPreBranchArcs.count(mcMesh.edge_handle(ha)))
                    {
                        if ((walkedDir & nextP.dirs1) == UVWDir::NONE)
                            continue; // Must be some unnecessary transition, there must be another one
                    }
                    if ((walkedDir & nextP.dirs1) == UVWDir::NONE)
                    {
                        if (allowedBranchOffs.count(ha))
                        {
                            nextP.branchedOff = true;
                            nextP.monotonousDirs = nextP.monotonousDirs & walkedDir;
                        }
                        else
                            // Dont allow branch-off at another arc than those allowed
                            continue;
                    }
                }
                if (nextP.branchedOff)
                    nextP.path.emplace_back(ha);

                nextP.delta += (_considerDoubleIntLength ? 2 : 1)
                               * mcMeshProps().get<ARC_INT_LENGTH>(mcMesh.edge_handle(ha)) * toVec(walkedDir);

                nextP.walkedDirs = nextP.walkedDirs | walkedDir;
                nextP.monotonousDirs = nextP.monotonousDirs & ~(-walkedDir);

                // Do not walk unnecessary circles
                if (nextP.monotonousDirs == UVWDir::NONE)
                    continue;

                nextP.dir2walkedArcs[walkedDir].emplace_back(mcMesh.edge_handle(ha));

                // Walk along the link for free, but other arcs accumulate distance
                // if (nextP.branchedOff)
                if ((walkedDir & nextP.dirs1) == UVWDir::NONE)
                    nextP.length += mcMeshProps().get<ARC_DBL_LENGTH>(mcMesh.edge_handle(ha));

                pathQ.push(nextP);
                break; // Only push the edge once (still need to check each bRef)
            }
        }
    }
}

void StructurePreserver::determineNextHalfarcs(const WeaklyMonotonousPath& pathCurrent,
                                               map<CH, Transition>& b2trans,
                                               map<HEH, vector<CH>>& ha2bRef) const
{
    auto& mcMesh = mcMeshProps().mesh();

    b2trans = map<CH, Transition>({{pathCurrent.bRefCurrent, {pathCurrent.transCurrent}}});

    // Floodfill blocks around n, storing current transition for each expanded block
    list<pair<CH, Transition>> bQ({{pathCurrent.bRefCurrent, pathCurrent.transCurrent}});

    while (!bQ.empty())
    {
        auto b2t = bQ.front();
        bQ.pop_front();

        for (HFH hp : mcMesh.cell_halffaces(b2t.first))
        {
            HFH hpOpp = mcMesh.opposite_halfface_handle(hp);
            CH bNext = mcMesh.incident_cell(hpOpp);
            if (!bNext.is_valid() || b2trans.find(bNext) != b2trans.end())
                continue;
            if (!contains(mcMesh.halfface_vertices(hp), pathCurrent.n))
                continue;

            // Check if overlapping with current path
            if (!checkPatchOverlap(pathCurrent, hp, b2t.second))
                continue;

            Transition trans = b2t.second.chain(mcMeshProps().hpTransition<PATCH_TRANSITION>(hp));

            bool exists = b2trans.find(bNext) != b2trans.end();
            if (exists)
                continue;
            b2trans[bNext] = trans;

            bQ.push_back({bNext, trans});
        }
    }
    ha2bRef.clear();
    for (HEH ha : mcMesh.outgoing_halfedges(pathCurrent.n))
    {
        if (!pathCurrent.path.empty() && pathCurrent.path.back() == mcMesh.opposite_halfedge_handle(ha))
            continue;
        for (CH b : mcMesh.halfedge_cells(ha))
        {
            auto it = b2trans.find(b);
            if (it != b2trans.end())
            {
                ha2bRef[ha].emplace_back(b);
            }
        }
    }
}

void StructurePreserver::traceUnproblematicPaths(const IntegerGridGraphNode& c0,
                                                 const CH& bRef0,
                                                 const Vec3i& i0Min,
                                                 const Vec3i& i0Max,
                                                 const set<IntegerGridGraphNode>& potentiallyUnproblematicNodes,
                                                 set<IntegerGridGraphNode>& unproblematicNodes)
{
    auto& mcMesh = mcMeshProps().mesh();
    using PathQueue = std::
        priority_queue<IntegerGridGraphPath, std::deque<IntegerGridGraphPath>, GreaterIntegerGridGraphPathCompare>;

    set<IntegerGridGraphNode> csVisited;

    PathQueue pQ;
    IntegerGridGraphPath pStart = initializePathStart(c0, bRef0, i0Min, i0Max, false);
    pQ.push(pStart);

    while (!pQ.empty())
    {
        IntegerGridGraphPath pIn = pQ.top();
        pQ.pop();

        if (csVisited.count(pIn.c))
            continue;
        csVisited.insert(pIn.c);
        unproblematicNodes.insert(pIn.c);

        map<CH, vector<Transition>> bRef2trans;
        fetchTransitionsAroundEntity(pIn.c, pIn.bRef, bRef2trans);

        VH nRef = pIn.c.nodes(mcMesh).front();
        for (const auto& [bRef, transitions] : bRef2trans)
        {
            auto& trans = transitions.front();
            // for (const Transition& trans : transitions)
            {
                Vec3i shift = _b2n2igm[bRef.idx()].at(nRef) - trans.rotate(_b2n2igm[pIn.bRef.idx()].at(nRef));

                DLOG(INFO) << "Shift determined: " << shift;

                Vec3i iLocal = trans.rotate(pIn.i) + shift;
                Vec3i iMinLocal = trans.rotate(pIn.iMin) + shift;
                Vec3i iMaxLocal = trans.rotate(pIn.iMax) + shift;
                Vec3Q uLocal = trans.apply(pIn.u);
                Vec3Q u0Local = trans.apply(pIn.u0);
                if (pIn.c.n.is_valid())
                    if (uLocal != _b2n2uvw[bRef.idx()].at(pIn.c.n))
                        throw std::logic_error("What?");
                Transition transTotal = pIn.trans.chain(trans);

                map<IntegerGridGraphNode, pairTT<Vec3i>> csOut2interval;
                auto& p2igmMin = _b2p2igmMin[bRef.idx()];
                auto& p2igmMax = _b2p2igmMax[bRef.idx()];
                auto& a2igmMin = _b2a2igmMin[bRef.idx()];
                auto& a2igmMax = _b2a2igmMax[bRef.idx()];
                auto& n2IGM = _b2n2igm[bRef.idx()];
                for (HFH hp : mcMesh.cell_halffaces(bRef))
                {
                    FH p = mcMesh.face_handle(hp);
                    Vec3i pMin = p2igmMin.at(p);
                    Vec3i pMax = p2igmMax.at(p);
                    if (!bboxOverlap(iMinLocal, iMaxLocal, pMin, pMax))
                        continue;
                    UVWDir patchPlusDir = ~halfpatchNormalDir(hp) & UVWDir::POS_U_POS_V_POS_W;
                    pairTT<Vec3i> patchInnerOverlap = {pMin + toVec(patchPlusDir), pMax - toVec(patchPlusDir)};
                    if (bboxOverlap(iMinLocal, iMaxLocal, patchInnerOverlap.first, patchInnerOverlap.second))
                    {
                        for (int i = 0; i < 3; i++)
                        {
                            patchInnerOverlap.first[i] = std::max(patchInnerOverlap.first[i], iMinLocal[i]);
                            patchInnerOverlap.second[i] = std::min(patchInnerOverlap.second[i], iMaxLocal[i]);
                        }
                        csOut2interval[p] = patchInnerOverlap;
                    }
                    else
                        csOut2interval[p] = {Vec3i(1, 1, 1), Vec3i(-1, -1, -1)};

                    for (HEH ha : mcMesh.halfface_halfedges(hp))
                    {
                        EH a = mcMesh.edge_handle(ha);
                        if (csOut2interval.count(a))
                            continue;
                        Vec3i aMin = a2igmMin.at(a);
                        Vec3i aMax = a2igmMax.at(a);
                        if (!bboxOverlap(iMinLocal, iMaxLocal, aMin, aMax))
                            continue;
                        UVWDir arcPlusDir = halfarcDirInBlock(ha, bRef);
                        arcPlusDir = (arcPlusDir | -arcPlusDir) & UVWDir::POS_U_POS_V_POS_W;
                        pairTT<Vec3i> arcInnerOverlap = {aMin + toVec(arcPlusDir), aMax - toVec(arcPlusDir)};
                        if (bboxOverlap(iMinLocal, iMaxLocal, arcInnerOverlap.first, arcInnerOverlap.second))
                        {
                            for (int i = 0; i < 3; i++)
                            {
                                arcInnerOverlap.first[i] = std::max(arcInnerOverlap.first[i], iMinLocal[i]);
                                arcInnerOverlap.second[i] = std::min(arcInnerOverlap.second[i], iMaxLocal[i]);
                            }
                            csOut2interval[a] = arcInnerOverlap;
                        }
                        else
                            csOut2interval[a] = {Vec3i(1, 1, 1), Vec3i(-1, -1, -1)};

                        for (VH n : mcMesh.halfedge_vertices(ha))
                        {
                            Vec3i igm = n2IGM.at(n);
                            if (bboxOverlap(iMinLocal, iMaxLocal, igm, igm))
                                csOut2interval[n] = {igm, igm};
                        }
                    }
                }

                for (const auto& kv : csOut2interval)
                {
                    auto& cOut = kv.first;
                    auto& range = kv.second;
                    if (range.first[0] > range.second[0])
                        continue;
                    if (csVisited.count(cOut))
                        continue;
                    if (!potentiallyUnproblematicNodes.count(cOut))
                        continue;

                    if (cOut.n.is_valid())
                    {
                        if (pIn.c.n.is_valid())
                        {
                            if (cOut.n != pIn.c.n
                                && !containsMatching(mcMesh.outgoing_halfedges(pIn.c.n),
                                                     [&, this](const HEH& ha)
                                                     {
                                                         return mcMesh.to_vertex_handle(ha) == cOut.n
                                                                && potentiallyUnproblematicNodes.count(
                                                                    mcMesh.edge_handle(ha));
                                                     })
                                && !containsMatching(mcMesh.vertex_faces(pIn.c.n),
                                                     [&, this](const FH& p)
                                                     {
                                                         return potentiallyUnproblematicNodes.count(p)
                                                                && contains(mcMesh.face_vertices(p), cOut.n);
                                                     }))
                                continue;
                        }
                        else if (pIn.c.a.is_valid())
                        {
                            if (!containsMatching(mcMesh.edge_faces(pIn.c.a),
                                                  [&, this](const FH& p)
                                                  {
                                                      return potentiallyUnproblematicNodes.count(p)
                                                             && contains(mcMesh.face_vertices(p), cOut.n);
                                                  }))
                                continue;
                        }
                        else if (pIn.c.p.is_valid())
                        {
                            DLOG(WARNING) << "I think this should not happen";
                            continue;
                        }
                    }
                    else if (cOut.a.is_valid())
                    {
                        if (pIn.c.n.is_valid())
                        {
                            if (!containsMatching(mcMesh.vertex_faces(pIn.c.n),
                                                  [&, this](const FH& p)
                                                  {
                                                      return potentiallyUnproblematicNodes.count(p)
                                                             && contains(mcMesh.face_edges(p), cOut.a);
                                                  }))
                                continue;
                        }
                        else if (pIn.c.a.is_valid() && pIn.c.a != cOut.a)
                        {
                            if (!containsMatching(mcMesh.edge_faces(pIn.c.a),
                                                  [&, this](const FH& p)
                                                  {
                                                      return potentiallyUnproblematicNodes.count(p)
                                                             && contains(mcMesh.face_edges(p), cOut.a);
                                                  }))
                                continue;
                        }
                        else if (pIn.c.p.is_valid())
                        {
                            DLOG(WARNING) << "I think this should not happen";
                            continue;
                        }
                    }
                    else
                    {
                        continue;
                    }

                    auto [nOut, minSqDist] = findClosestNode(cOut, bRef, uLocal);
                    Vec3Q uOut = _b2n2uvw[bRef.idx()].at(nOut);
                    Vec3Q delta = uOut - uLocal;
                    double lengthOut = pIn.length + Vec3Q2d(delta).norm();

                    IntegerGridGraphPath pOut;
                    pOut.i = iLocal;
                    pOut.iMin = range.first;
                    pOut.iMax = range.second;
                    pOut.u = uOut;
                    pOut.u0 = u0Local;
                    pOut.c = cOut;
                    pOut.bRef = bRef;
                    pOut.trans = transTotal;
                    pOut.hops = pIn.hops + 1;
                    pOut.length = lengthOut;
                    pQ.push(pOut);
                }
            }
        }
    }
}

bool StructurePreserver::traceExhaustPaths(const IntegerGridGraphNode& c0,
                                           const CH& bRef0,
                                           const Vec3i& i0Min,
                                           const Vec3i& i0Max,
                                           const set<IntegerGridGraphNode>& unproblematicNodes,
                                           vector<vector<pair<int, EH>>>& nonZeroSumArcs,
                                           vector<vector<pair<int, EH>>>& failsafeNonZeroSumArcs) const
{
    bool sep = false;
    auto& mcMesh = mcMeshProps().mesh();
    using PathQueue = std::
        priority_queue<IntegerGridGraphPath, std::deque<IntegerGridGraphPath>, GreaterIntegerGridGraphPathCompare>;

    set<IntegerGridGraphNode> csVisited; // Prevent visiting, if shorter path exists)

    set<pair<IntegerGridGraphNode, CH>> cornersVisited; // only necessary if considering cut surface cycles

    PathQueue pQ;
    IntegerGridGraphPath pStart = initializePathStart(c0, bRef0, i0Min, i0Max, true);
    pQ.push(pStart);

    while (!pQ.empty())
    {
        IntegerGridGraphPath pIn = pQ.top();
        pQ.pop();

        if (csVisited.count(pIn.c))
            continue;
        csVisited.insert(pIn.c);

        if (pIn.c.n.is_valid() && isCriticalNode[pIn.c.n.idx()])
            if (isStructureAlteringZeroPath(pIn, unproblematicNodes, nonZeroSumArcs, failsafeNonZeroSumArcs))
            {
                sep = true;
                continue;
            }

        map<CH, vector<Transition>> bRef2trans;
        fetchTransitionsAroundEntity(pIn.c, pIn.bRef, bRef2trans);

        VH nRef = pIn.c.nodes(mcMesh).front();
        for (const auto& [bRef, transitions] : bRef2trans)
        {
            auto& trans = transitions.front();
            // for (const Transition& trans : transitions)
            {
                Vec3i shift = _b2n2igm[bRef.idx()].at(nRef) - trans.rotate(_b2n2igm[pIn.bRef.idx()].at(nRef));
                Vec3i iLocal = trans.rotate(pIn.i) + shift;
                Vec3i iMinLocal = trans.rotate(pIn.iMin) + shift;
                Vec3i iMaxLocal = trans.rotate(pIn.iMax) + shift;
                Vec3Q uLocal = trans.apply(pIn.u);
                Vec3Q u0Local = trans.apply(pIn.u0);
                Transition transTotal = pIn.trans.chain(trans);

                map<IntegerGridGraphNode, pairTT<Vec3i>> csOut2interval;
                auto& p2igmMin = _b2p2igmMin[bRef.idx()];
                auto& p2igmMax = _b2p2igmMax[bRef.idx()];
                auto& a2igmMin = _b2a2igmMin[bRef.idx()];
                auto& a2igmMax = _b2a2igmMax[bRef.idx()];
                auto& n2IGM = _b2n2igm[bRef.idx()];
                for (HFH hp : mcMesh.cell_halffaces(bRef))
                {
                    FH p = mcMesh.face_handle(hp);
                    Vec3i pMin = p2igmMin.at(p);
                    Vec3i pMax = p2igmMax.at(p);
                    if (!bboxOverlap(iMinLocal, iMaxLocal, pMin, pMax))
                        continue;
                    UVWDir patchPlusDir = ~halfpatchNormalDir(hp) & UVWDir::POS_U_POS_V_POS_W;
                    pairTT<Vec3i> patchInnerOverlap = {pMin + toVec(patchPlusDir), pMax - toVec(patchPlusDir)};
                    if (bboxOverlap(iMinLocal, iMaxLocal, patchInnerOverlap.first, patchInnerOverlap.second))
                    {
                        for (int i = 0; i < 3; i++)
                        {
                            patchInnerOverlap.first[i] = std::max(patchInnerOverlap.first[i], iMinLocal[i]);
                            patchInnerOverlap.second[i] = std::min(patchInnerOverlap.second[i], iMaxLocal[i]);
                        }
                        csOut2interval[p] = patchInnerOverlap;
                    }
                    else
                        csOut2interval[p] = {Vec3i(1, 1, 1), Vec3i(-1, -1, -1)};

                    for (HEH ha : mcMesh.halfface_halfedges(hp))
                    {
                        EH a = mcMesh.edge_handle(ha);
                        if (csOut2interval.count(a))
                            continue;
                        Vec3i aMin = a2igmMin.at(a);
                        Vec3i aMax = a2igmMax.at(a);
                        if (!bboxOverlap(iMinLocal, iMaxLocal, aMin, aMax))
                            continue;
                        UVWDir arcPlusDir = halfarcDirInBlock(ha, bRef);
                        arcPlusDir = (arcPlusDir | -arcPlusDir) & UVWDir::POS_U_POS_V_POS_W;
                        pairTT<Vec3i> arcInnerOverlap = {aMin + toVec(arcPlusDir), aMax - toVec(arcPlusDir)};
                        if (bboxOverlap(iMinLocal, iMaxLocal, arcInnerOverlap.first, arcInnerOverlap.second))
                        {
                            for (int i = 0; i < 3; i++)
                            {
                                arcInnerOverlap.first[i] = std::max(arcInnerOverlap.first[i], iMinLocal[i]);
                                arcInnerOverlap.second[i] = std::min(arcInnerOverlap.second[i], iMaxLocal[i]);
                            }
                            csOut2interval[a] = arcInnerOverlap;
                        }
                        else
                            csOut2interval[a] = {Vec3i(1, 1, 1), Vec3i(-1, -1, -1)};

                        for (VH n : mcMesh.halfedge_vertices(ha))
                        {
                            Vec3i igm = n2IGM.at(n);
                            if (bboxOverlap(iMinLocal, iMaxLocal, igm, igm))
                                csOut2interval[n] = {igm, igm};
                        }
                    }
                }

                for (const auto& kv : csOut2interval)
                {
                    auto& cOut = kv.first;
                    auto& range = kv.second;
                    if (range.first[0] > range.second[0])
                        continue;
                    if (csVisited.count(cOut))
                        continue;

                    auto [nOut, minSqDist] = findClosestNode(cOut, bRef, uLocal);
                    Vec3Q uOut = _b2n2uvw[bRef.idx()].at(nOut);
                    Vec3Q delta = uOut - uLocal;
                    double lengthOut = pIn.length + Vec3Q2d(delta).norm();
                    Precursor precursorOut = {};
                    precursorOut.cPrev = pIn.c;
                    precursorOut.transPrev = pIn.trans;
                    precursorOut.bRefPrev = pIn.bRef;
                    precursorOut.deltaToPrev = -delta;
                    precursorOut.transToPrev = trans.invert();

                    IntegerGridGraphPath pOut;
                    pOut.i = iLocal;
                    pOut.iMin = range.first;
                    pOut.iMax = range.second;
                    pOut.u = uOut;
                    pOut.u0 = u0Local;
                    pOut.c = cOut;
                    pOut.bRef = bRef;
                    pOut.trans = transTotal;
                    pOut.precursors = pIn.precursors;
                    pOut.precursors.push_back(precursorOut);
                    pOut.hops = pIn.hops + 1;
                    pOut.length = lengthOut;
                    pQ.push(pOut);
                }
            }
        }
    }

    return sep;
}

bool StructurePreserver::isStructureAlteringZeroPath(const IntegerGridGraphPath& pCurrent,
                                                     const set<IntegerGridGraphNode>& unproblematicNodes,
                                                     vector<vector<pair<int, EH>>>& nonZeroSumArcs,
                                                     vector<vector<pair<int, EH>>>& failsafeNonZeroSumArcs) const
{
    auto& mcMesh = mcMeshProps().mesh();
    bool sepViolated = false;

    if (pCurrent.c.n.is_valid() && !unproblematicNodes.count(pCurrent.c.n))
    {
        if (isCriticalNode[pCurrent.c.n.idx()])
            sepViolated = true;
        if (!sepViolated)
            for (HEH ha : mcMesh.outgoing_halfedges(pCurrent.c.n))
            {
                EH a = mcMesh.edge_handle(ha);
                if (isCriticalArc[a.idx()] && !unproblematicNodes.count(a))
                {
                    sepViolated = true;
                    break;
                }
            }
        if (!sepViolated)
            for (FH p : mcMesh.vertex_faces(pCurrent.c.n))
                if (isCriticalPatch[p.idx()] && !unproblematicNodes.count(p))
                {
                    sepViolated = true;
                    break;
                }
    }
    else if (pCurrent.c.a.is_valid() && !unproblematicNodes.count(pCurrent.c.a))
    {
        if (isCriticalArc[pCurrent.c.a.idx()] && !unproblematicNodes.count(pCurrent.c.a))
            sepViolated = true;
        if (!sepViolated)
            for (FH p : mcMesh.edge_faces(pCurrent.c.a))
                if (isCriticalPatch[p.idx()] && !unproblematicNodes.count(p))
                {
                    sepViolated = true;
                    break;
                }
    }
    else if (pCurrent.c.p.is_valid())
    {
        if (isCriticalPatch[pCurrent.c.p.idx()] && !unproblematicNodes.count(pCurrent.c.p))
            sepViolated = true;
    }
    if (!sepViolated)
        return false;

    vector<vector<pair<HEH, UVWDir>>> pathHaDirs;
    reconstructArcPaths(pCurrent, pathHaDirs, true);

    set<EH> as;
    for (int i = 0; i < (int)pathHaDirs.size(); i++)
        for (auto& [ha, dir] : pathHaDirs[i])
            as.insert(mcMesh.edge_handle(ha));

    // Dirty workaround for something that went wrong earlier
    if (as.empty())
        return false;

    failsafeNonZeroSumArcs.emplace_back();
    auto& failsafeConstraint = failsafeNonZeroSumArcs.back();
    for (EH a : as)
        failsafeConstraint.push_back({1, a});
    nonZeroSumArcs.emplace_back();
    auto& normalConstraint = nonZeroSumArcs.back();
    normalConstraint = failsafeConstraint;

    return true;
}

void StructurePreserver::reconstructArcPaths(const IntegerGridGraphPath& pCurrent,
                                             vector<vector<pair<HEH, UVWDir>>>& pathHaDirs,
                                             bool excludeNonZeroArcs) const
{
    auto& mcMesh = mcMeshProps().mesh();

    // First, reconstruct the full sequence, starting from current c, ending at path start
    vector<IntegerGridGraphNode> cs; // sequence of implicit graph nodes
    vector<Vec3Q> uvws;              // sequence of UVWs in current (!) coordinate system
    vector<CH> bsRef;                // reference block sequence
    vector<Transition> transitions;  // sequence of transitions for (i.e. "into") each block in sequence (from current)
    vector<VH> ns;                   // sequence of nodes, closest in UVW coordinates to each c in sequence

    // Fill the vectors in reverse (starting from cCurrent, and then traverse them in reverse afterwards
    {
        IntegerGridGraphNode c = pCurrent.c;
        Vec3Q uvw = pCurrent.u;
        CH bRef = pCurrent.bRef;
        Transition trans{};
        auto [n, minSqDist] = findClosestNode(c, bRef, uvw);
        assert(minSqDist == 0.0);

        cs.push_back(c);
        uvws.push_back(uvw);
        bsRef.push_back(bRef);
        transitions.push_back(trans);
        ns.push_back(n);

        for (int precursorIndex = pCurrent.precursors.size() - 1; precursorIndex > 0; precursorIndex--)
        {
            const Precursor& pre = pCurrent.precursors[precursorIndex];

            c = pre.cPrev;
            uvw = uvw + trans.invert().rotate(pre.deltaToPrev);
            bRef = pre.bRefPrev;
            trans = trans.chain(pre.transToPrev);
            std::tie(n, minSqDist) = findClosestNode(c, bRef, trans.apply(uvw));
            assert(minSqDist == 0.0);

            cs.push_back(c);
            uvws.push_back(uvw);
            bsRef.push_back(bRef);
            transitions.push_back(trans);
            ns.push_back(n);
        }
    }

    // Then, from the reconstructed sequence, construct an arc-based shortest path (blockwise, then glue
    // together)
    pathHaDirs = vector<vector<pair<HEH, UVWDir>>>(ns.size() - 1);

    struct Path
    {
        VH nCurr = VH();
        HEH precursor = HEH();
        double length = 0.0;
    };
    auto pathCmp = [](const Path& p1, const Path& p2) { return p1.length > p2.length; };
    for (int i = 0; i < (int)ns.size() - 1; i++)
    {
        VH nStart = ns[i];
        VH nEnd = ns[i + 1];
        if (nStart == nEnd)
            continue;

        CH bRef = bsRef[i];
        assert(contains(mcMesh.cell_vertices(bRef), nStart) && contains(mcMesh.cell_vertices(bRef), nEnd));

        UVWDir deltaDir = toDir(_b2n2uvw[bRef.idx()].at(nEnd) - _b2n2uvw[bRef.idx()].at(nStart));
        deltaDir = deltaDir | -deltaDir;
        assert(dim(deltaDir) != 0);

        // Now find shortest path between nStart and nEnd on boundary of bRef
        using ArcPathQueue = std::priority_queue<Path, std::deque<Path>, decltype(pathCmp)>;

        ArcPathQueue apQ(pathCmp);

        map<VH, HEH> n2precursor;
        Path p0;
        p0.nCurr = nEnd;
        p0.precursor = HEH();
        p0.length = 0.0;
        apQ.push(p0);

        while (!apQ.empty())
        {
            Path pIn = apQ.top();
            apQ.pop();

            if (n2precursor.count(pIn.nCurr))
                continue;

            n2precursor[pIn.nCurr] = pIn.precursor;
            if (pIn.nCurr == nStart)
                break;

            for (HEH haOut : mcMesh.outgoing_halfedges(pIn.nCurr))
            {
                VH nTo = mcMesh.to_vertex_handle(haOut);
                if (n2precursor.count(nTo) || !contains(mcMesh.halfedge_cells(haOut), bRef))
                    continue;
                Path pOut;
                pOut.nCurr = nTo;
                pOut.precursor = mcMesh.opposite_halfedge_handle(haOut);
                pOut.length = pIn.length + mcMeshProps().get<ARC_DBL_LENGTH>(mcMesh.edge_handle(haOut));
                apQ.push(pOut);
            }
        }

        assert(n2precursor.count(nStart));
        HEH precursor = n2precursor.at(nStart);
        while (precursor.is_valid())
        {
            if ((halfarcDirInBlock(precursor, bRef) & deltaDir) != UVWDir::NONE
                && (!excludeNonZeroArcs || mcMeshProps().get<ARC_INT_LENGTH>(mcMesh.edge_handle(precursor)) == 0))
                pathHaDirs[i].push_back(
                    {precursor, -transitions[i].invert().rotate(halfarcDirInBlock(precursor, bRef))});
            precursor = n2precursor.at(mcMesh.to_vertex_handle(precursor));
        }
    }

    assert(pCurrent.u - pCurrent.u0 == uvws[0] - uvws.back());
}

void StructurePreserver::precompute()
{
    auto& mcMesh = mcMeshProps().mesh();

    _b2n2igm = vector<map<VH, Vec3i>>(mcMesh.n_cells());
    _b2n2uvw = vector<map<VH, Vec3Q>>(mcMesh.n_cells());
    _b2a2igmMin = vector<map<EH, Vec3i>>(mcMesh.n_cells());
    _b2a2igmMax = vector<map<EH, Vec3i>>(mcMesh.n_cells());
    _b2p2igmMin = vector<map<FH, Vec3i>>(mcMesh.n_cells());
    _b2p2igmMax = vector<map<FH, Vec3i>>(mcMesh.n_cells());

    for (CH b : mcMesh.cells())
    {
        for (VH n : mcMesh.cell_vertices(b))
        {
            _b2n2uvw[b.idx()][n] = nodeUVWinBlock(n, b);
            _bn2transitionsAroundNode[{b, n}] = determineTransitionsAroundNode(n, b, Transition());
        }

        VH n0 = *mcMesh.cv_iter(b);
        list<pair<VH, Vec3i>> nQ;
        nQ.emplace_back(n0, Vec3i(0, 0, 0));
        _b2n2igm[b.idx()][n0] = Vec3i(0, 0, 0);
        while (!nQ.empty())
        {
            auto [n, i] = nQ.front();
            nQ.pop_front();

            for (HEH haOut : mcMesh.outgoing_halfedges(n))
            {
                VH nTo = mcMesh.to_vertex_handle(haOut);
                if (_b2n2igm[b.idx()].count(nTo) || !contains(mcMesh.halfedge_cells(haOut), b))
                    continue;
                Vec3i igm = i
                            + (_considerDoubleIntLength ? 2 : 1)
                                  * mcMeshProps().get<ARC_INT_LENGTH>(mcMesh.edge_handle(haOut))
                                  * toVec(halfarcDirInBlock(haOut, b));
                nQ.push_back({nTo, igm});
                _b2n2igm[b.idx()][nTo] = igm;
            }
        }

        for (EH a : mcMesh.cell_edges(b))
        {
            Vec3i igmMin(INT_MAX, INT_MAX, INT_MAX);
            Vec3i igmMax(INT_MIN, INT_MIN, INT_MIN);
            for (VH n : mcMesh.edge_vertices(a))
            {
                Vec3i& igm = _b2n2igm[b.idx()][n];
                for (int coord = 0; coord < 3; coord++)
                {
                    igmMin[coord] = std::min(igmMin[coord], igm[coord]);
                    igmMax[coord] = std::max(igmMax[coord], igm[coord]);
                }
            }
            _b2a2igmMin[b.idx()][a] = igmMin;
            _b2a2igmMax[b.idx()][a] = igmMax;
            _ba2transitionsAroundArc[{b, a}] = determineTransitionsAroundArc(a, b, Transition());
        }

        for (FH p : mcMesh.cell_faces(b))
        {
            Vec3i igmMin(INT_MAX, INT_MAX, INT_MAX);
            Vec3i igmMax(INT_MIN, INT_MIN, INT_MIN);
            for (VH n : mcMesh.face_vertices(p))
            {
                Vec3i& igm = _b2n2igm[b.idx()][n];
                for (int coord = 0; coord < 3; coord++)
                {
                    igmMin[coord] = std::min(igmMin[coord], igm[coord]);
                    igmMax[coord] = std::max(igmMax[coord], igm[coord]);
                }
            }
            _b2p2igmMin[b.idx()][p] = igmMin;
            _b2p2igmMax[b.idx()][p] = igmMax;
        }
    }
}

void StructurePreserver::findSingularityIndexViolationPaths(vector<vector<pair<int, EH>>>& nonZeroSumArcs,
                                                            vector<vector<pair<int, EH>>>& failsafeNonZeroSumArcs) const
{
    const MCMesh& mcMesh = mcMeshProps().mesh();

    set<int> visitedLinks;
    for (auto& singularLink : singularEntities)
    {
        if (singularLink.dim != 1 || visitedLinks.count(singularLink.id))
            continue;

        bool violation = false;
        // Now, sample the link in integer steps
        for (HEH ha : singularLink.pathHas)
        {
            EH a = mcMesh.edge_handle(ha);
            int l = (_considerDoubleIntLength ? 2 : 1) * mcMeshProps().get<ARC_INT_LENGTH>(a);
            if (l == 0)
                continue;
            CH bRef = *mcMesh.hec_iter(ha);

            auto& n2IGM = _b2n2igm[bRef.idx()];
            Vec3i igmStart = n2IGM.at(mcMesh.from_vertex_handle(ha));
            UVWDir dirHa = halfarcDirInBlock(ha, bRef);
            for (int i = 1; i < l; i++)
            {
                if (i % 2 == 0)
                    continue;
                Vec3i i0 = igmStart + toVec(dirHa) * i;
                DLOG(INFO) << "Index-tracing from singular link " << singularLink.id << " arc " << a << " between "
                           << mcMeshProps().get<NODE_MESH_VERTEX>(mcMesh.from_vertex_handle(ha)) << " - "
                           << mcMeshProps().get<NODE_MESH_VERTEX>(mcMesh.to_vertex_handle(ha)) << " interval " << i
                           << " / " << l;
                violation = violation
                            || traceSingularityIndexPaths({a},
                                                          bRef,
                                                          i0,
                                                          dirHa | -dirHa,
                                                          singularLink.id,
                                                          nonZeroSumArcs,
                                                          failsafeNonZeroSumArcs,
                                                          visitedLinks);
                if (violation)
                    break;
            }
            if (violation)
                break;
        }
    }
}

bool StructurePreserver::traceSingularityIndexPaths(const IntegerGridGraphNode& c0,
                                                    const CH& bRef0,
                                                    const Vec3i& i0,
                                                    const UVWDir dirsParallel0,
                                                    int startID,
                                                    vector<vector<pair<int, EH>>>& nonZeroSumArcs,
                                                    vector<vector<pair<int, EH>>>& failsafeNonZeroSumArcs,
                                                    set<int>& visitedLinks) const
{
    auto& mcMesh = mcMeshProps().mesh();
    using PathQueue = std::
        priority_queue<IntegerGridGraphPath, std::deque<IntegerGridGraphPath>, GreaterIntegerGridGraphPathCompare>;

    VH n0;
    double minSqDist = DBL_MAX;
    for (VH n : c0.nodes(mcMesh))
    {
        double sqDist = (i0 - _b2n2igm[bRef0.idx()].at(n)).sqrnorm();
        if (sqDist < minSqDist)
        {
            minSqDist = sqDist;
            n0 = n;
        }
    }
    Vec3Q u0 = _b2n2uvw[bRef0.idx()].at(n0);

    map<IntegerGridGraphNode, double> c2minDist;      // Prevent visiting, if shorter path exists
    map<IntegerGridGraphNode, Precursor> c2precursor; // Prevent visiting, if path different from shortest path
    map<IntegerGridGraphNode, set<Vec3i>>
        c2clearedTransitions; // Prevent visiting, if transitions along path identical to previous ones

    PathQueue pQ;
    IntegerGridGraphPath pStart = {};
    pStart.i = i0;
    pStart.u = u0;
    pStart.u0 = u0;
    pStart.c = c0;
    pStart.bRef = bRef0;
    pStart.trans = Transition();
    pStart.dirs = UVWDir::NONE;
    pStart.dirsParallel = dirsParallel0;
    pStart.hops = 0;
    pStart.length = 0.0;
    Precursor emptyPre = {};
    pStart.precursors = {emptyPre};
    pStart.criticalStartID = startID;
    pStart.partLength = 0.0;
    pStart.dirsTotal = UVWDir::NONE;
    pQ.push(pStart);

    map<pairTT<int>, IntegerGridGraphPath> linkPair2connectingPath;
    set<int> linksCluster;
    linksCluster.insert(startID);

    map<pairTT<EH>, IntegerGridGraphPath> aPair2connectingPath;
    set<EH> asCluster;
    asCluster.insert(c0.a);

    while (!pQ.empty())
    {
        IntegerGridGraphPath pIn = pQ.top();
        pQ.pop();
        assert(!pIn.c.n.is_valid());

        if (!c2precursor.count(pIn.c))
        {
            c2precursor[pIn.c] = pIn.precursors.back();
            c2minDist[pIn.c] = pIn.length;
        }
        else if (c2minDist.at(pIn.c) < pIn.length || c2precursor.at(pIn.c).cPrev != pIn.precursors.back().cPrev
                 || c2clearedTransitions.at(pIn.c).count(pIn.trans.rotation))
            continue;
        c2clearedTransitions[pIn.c].insert(pIn.trans.rotation);

        {
            int criticalLinkID = -1;
            EH criticalArc;
            UVWDir violationDirs = UVWDir::NONE;
            UVWDir monotonousDirsIn = pIn.dirs & ~(-pIn.dirs);
            if (pIn.c.n.is_valid())
            {
                for (HEH ha : mcMesh.outgoing_halfedges(pIn.c.n))
                    if (isSingularArc[mcMesh.edge_handle(ha).idx()]
                        // && !linksCluster.count(a2singularLinkIdx.at(mcMesh.edge_handle(ha)))
                        && !asCluster.count(mcMesh.edge_handle(ha))
                        && mcMeshProps().get<ARC_INT_LENGTH>(mcMesh.edge_handle(ha)) != 0)
                    {
                        CH bRef = *mcMesh.hec_iter(ha);
                        for (auto& trans : _bn2transitionsAroundNode.at({pIn.bRef, pIn.c.n}).at(bRef))
                        {
                            UVWDir dirA = trans.invert().rotate(halfarcDirInBlock(ha, bRef));
                            dirA = dirA | -dirA;
                            if (dirA == pIn.dirsParallel && ((~pIn.dirsParallel & monotonousDirsIn) != UVWDir::NONE))
                            {
                                DLOG(INFO) << "Found unseparated critical arc" << mcMesh.edge_handle(ha) << " at node "
                                           << pIn.c.n << " / " << mcMeshProps().get<NODE_MESH_VERTEX>(pIn.c.n);
                                violationDirs = (~dirA & ~pIn.dirsParallel & pIn.dirs);
                                criticalLinkID = a2singularLinkIdx.at(mcMesh.edge_handle(ha));
                                criticalArc = mcMesh.edge_handle(ha);
                                break;
                            }
                        }
                        if (violationDirs != UVWDir::NONE)
                            break;
                    }
            }
            else if (pIn.c.a.is_valid())
            {
                if (isSingularArc[pIn.c.a.idx()]
                    // && !linksCluster.count(a2singularLinkIdx.at(pIn.c.a)
                    && !asCluster.count(pIn.c.a) && mcMeshProps().get<ARC_INT_LENGTH>(pIn.c.a) != 0)
                {
                    UVWDir dirA = halfarcDirInBlock(mcMesh.halfedge_handle(pIn.c.a, 0), pIn.bRef);
                    dirA = dirA | -dirA;
                    if (dirA == pIn.dirsParallel && ((~pIn.dirsParallel & monotonousDirsIn) != UVWDir::NONE))
                    {
                        DLOG(INFO) << "Found unseparated critical arc" << pIn.c.a;
                        violationDirs = (~dirA & ~pIn.dirsParallel & pIn.dirs);
                        criticalLinkID = a2singularLinkIdx.at(pIn.c.a);
                        criticalArc = pIn.c.a;
                    }
                }
            }
            if (violationDirs != UVWDir::NONE
                // && !linksCluster.count(criticalLinkID)
                && !asCluster.count(criticalArc))
            {
                linksCluster.insert(criticalLinkID);
                asCluster.insert(criticalArc);
                linkPair2connectingPath[{pIn.criticalStartID, criticalLinkID}] = pIn;
                aPair2connectingPath[{c0.a, criticalArc}] = pIn;

                pIn.u0 = pIn.u;
                pIn.dirs = UVWDir::NONE;
                pIn.precursors = {emptyPre};
                pIn.criticalStartID = criticalLinkID;
                pIn.partLength = 0.0;
            }
        }

        map<CH, vector<Transition>> bRef2trans;
        fetchTransitionsAroundEntity(pIn.c, pIn.bRef, bRef2trans);

        VH nRef = pIn.c.nodes(mcMesh).front();
        for (const auto& [bRef, transitions] : bRef2trans)
        {
            for (const Transition& trans : transitions)
            {
                if (trans.rotation[0] == -1 || trans.rotation[1] == -2 || trans.rotation[2] == -3)
                    continue;
                Vec3i shift = _b2n2igm[bRef.idx()].at(nRef) - trans.rotate(_b2n2igm[pIn.bRef.idx()].at(nRef));

                DLOG(INFO) << "Shift determined: " << shift;

                Vec3i iLocal = trans.rotate(pIn.i) + shift;
                Vec3Q uLocal = trans.apply(pIn.u);
                Vec3Q u0Local = trans.apply(pIn.u0);
                if (pIn.c.n.is_valid())
                    if (uLocal != _b2n2uvw[bRef.idx()].at(pIn.c.n))
                        throw std::logic_error("What?");
                Transition transTotal = pIn.trans.chain(trans);
                UVWDir dirsLocal = trans.rotate(pIn.dirs);
                UVWDir dirsTotalLocal = trans.rotate(pIn.dirsTotal);
                UVWDir dirsParallelLocal = trans.rotate(pIn.dirsParallel);

                set<IntegerGridGraphNode> csOut;
                auto& p2igmMin = _b2p2igmMin[bRef.idx()];
                auto& p2igmMax = _b2p2igmMax[bRef.idx()];
                auto& a2igmMin = _b2a2igmMin[bRef.idx()];
                auto& a2igmMax = _b2a2igmMax[bRef.idx()];
                auto& n2IGM = _b2n2igm[bRef.idx()];
                for (FH p : mcMesh.cell_faces(bRef))
                {
                    if (!bboxOverlap(iLocal, iLocal, p2igmMin.at(p), p2igmMax.at(p)))
                        continue;
                    DLOG(INFO) << "Overlap " << iLocal << " with patch " << p << " [" << p2igmMin.at(p) << ", "
                               << p2igmMax.at(p) << "]";

                    bool arcOrNodeAdded = false;
                    for (EH a : mcMesh.face_edges(p))
                    {
                        if (!bboxOverlap(iLocal, iLocal, a2igmMin.at(a), a2igmMax.at(a)))
                            continue;
                        DLOG(INFO) << "Overlap " << iLocal << " with arc " << a << " [" << a2igmMin.at(a) << ", "
                                   << a2igmMax.at(a) << "]";
                        arcOrNodeAdded = true;
                        bool nodeAdded = false;
                        for (VH n : mcMesh.edge_vertices(a))
                        {
                            if (iLocal != n2IGM.at(n))
                                continue;
                            DLOG(INFO) << "Overlap " << iLocal << " with node " << n;
                            nodeAdded = true;
                            csOut.insert({n});
                        }
                        if (!nodeAdded)
                            csOut.insert({a});
                    }
                    if (!arcOrNodeAdded)
                        csOut.insert({p});
                }

                for (const IntegerGridGraphNode& cOut : csOut)
                {
                    if (cOut == pIn.c)
                        continue;
                    if ((c2precursor.count(cOut) && c2precursor.at(cOut).cPrev != pIn.c)
                        || (c2clearedTransitions.count(cOut)
                            && c2clearedTransitions.at(cOut).count(transTotal.rotation)))
                        continue;
                    auto [nOut, minSqDist2] = findClosestNode(cOut, bRef, uLocal);
                    Vec3Q uOut = _b2n2uvw[bRef.idx()].at(nOut);
                    Vec3Q delta = uOut - uLocal;
                    UVWDir dirsOut = dirsLocal | toDir(delta);
                    UVWDir dirsTotalOut = dirsTotalLocal | toDir(delta);
                    double deltaLength = Vec3Q2d(delta).norm();
                    double lengthOut = pIn.length + deltaLength;
                    bool monoDirExists = (dirsTotalOut == UVWDir::NONE) || (dirsTotalOut != -dirsTotalOut);
                    if (!monoDirExists || (c2minDist.count(cOut) && lengthOut > c2minDist.at(cOut)))
                        continue;
                    Precursor precursorOut = {};
                    precursorOut.cPrev = pIn.c;
                    precursorOut.transPrev = pIn.trans;
                    precursorOut.bRefPrev = pIn.bRef;
                    precursorOut.deltaToPrev = -delta;
                    precursorOut.transToPrev = trans.invert();

                    IntegerGridGraphPath pOut;
                    pOut.i = iLocal;
                    pOut.u = uOut;
                    pOut.u0 = u0Local;
                    pOut.c = cOut;
                    pOut.bRef = bRef;
                    pOut.trans = transTotal;
                    pOut.dirs = dirsOut;
                    pOut.dirsParallel = dirsParallelLocal;
                    pOut.precursors = pIn.precursors;
                    pOut.precursors.push_back(precursorOut);
                    pOut.hops = pIn.hops + 1;
                    pOut.length = lengthOut;
                    pOut.criticalStartID = pIn.criticalStartID;
                    pStart.partLength = pIn.partLength + deltaLength;
                    pQ.push(pOut);
                }
                if (pIn.dirs == UVWDir::NONE)
                {
                    // No need to enumerate all possible transitions if no direction has been walked yet
                    break;
                }
            }
        }
    }

    int clusterIdxSum = 0;
    for (int id : linksCluster)
        visitedLinks.insert(id);
    for (EH a : asCluster)
        clusterIdxSum += mcMeshProps().get<IS_SINGULAR>(a);

    if (clusterIdxSum > _maxSingIdx || clusterIdxSum < _minSingIdx)
    {
        bool above = clusterIdxSum > _maxSingIdx;

        DLOG(INFO) << "Starting from ID " << startID << " found singularity cluster of size " << asCluster.size()
                   << " with idx sum: " << clusterIdxSum;

        map<int, set<int>> link2neighborLinks;
        for (auto& kv : linkPair2connectingPath)
        {
            auto& [link1, link2] = kv.first;
            link2neighborLinks[link1].insert(link2);
            link2neighborLinks[link2].insert(link1);
        }

        map<EH, set<EH>> a2neighborAs;
        for (auto& kv : aPair2connectingPath)
        {
            auto& [a1, a2] = kv.first;
            a2neighborAs[a1].insert(a2);
            a2neighborAs[a2].insert(a1);
        }

        set<pairTT<EH>> clusterSplittingPaths;
        for (auto& kv : aPair2connectingPath)
        {
            auto& [a1, a2] = kv.first;
            // graphsearch build 2 clusters

            int idxSumFront = 0;
            int idxSumBack = 0;

            set<EH> frontAs;
            set<EH> backAs;
            frontAs.insert(a1);
            backAs.insert(a2);
            for (bool front : {true, false})
            {
                int idxSum = front ? idxSumFront : idxSumBack;
                set<EH>& asVisited = front ? frontAs : backAs;
                list<EH> aQ({*asVisited.begin()});
                while (!aQ.empty())
                {
                    EH a = aQ.front();
                    DLOG(INFO) << "Link " << a << " contributes " << mcMeshProps().get<IS_SINGULAR>(a) << " to "
                               << (front ? "front" : "back");
                    idxSum += mcMeshProps().get<IS_SINGULAR>(a);
                    aQ.pop_front();
                    for (EH aNext : a2neighborAs.at(a))
                        if (!asVisited.count(aNext))
                        {
                            aQ.push_back(aNext);
                            asVisited.insert(aNext);
                        }
                }
            }

            if ((above && std::max(idxSumFront, idxSumBack) < clusterIdxSum)
                || (!above && std::min(idxSumFront, idxSumBack) > clusterIdxSum))
                clusterSplittingPaths.insert(kv.first);
        }
        DLOG(INFO) << "linkPair2connectingPath contains";
        for (auto& kv : linkPair2connectingPath)
            DLOG(INFO) << kv.first.first << ", " << kv.first.second;
        DLOG(INFO) << "Clustersplittingpath contains";
        for (auto path : clusterSplittingPaths)
            DLOG(INFO) << path.first << ", " << path.second;

        if (clusterSplittingPaths.empty())
            return false;

        map<EH, int> a2sign;
        map<EH, int> a2signFailsafe;

        for (pairTT<EH> aPair : clusterSplittingPaths)
        {
            auto& pCurrent = aPair2connectingPath.at(aPair);
            assert(pCurrent.c.a.is_valid());

#ifndef NDEBUG
            UVWDir violationDirs = UVWDir::NONE;
#endif
            UVWDir nonParallelDirs = UVWDir::NONE;
            UVWDir monotonousDirsIn = pCurrent.dirs & ~(-pCurrent.dirs);

            assert((monotonousDirsIn & ~(pCurrent.dirsParallel)) != UVWDir::NONE);
            if (isSingularArc[pCurrent.c.a.idx()])
            {
                UVWDir dirA = halfarcDirInBlock(mcMesh.halfedge_handle(pCurrent.c.a, 0), pCurrent.bRef);
                dirA = dirA | -dirA;
                if ((~dirA & ~pCurrent.dirsParallel & monotonousDirsIn) != UVWDir::NONE)
                {
#ifndef NDEBUG
                    violationDirs = (~dirA & ~pCurrent.dirsParallel & pCurrent.dirs);
#endif
                    nonParallelDirs = (~dirA & ~pCurrent.dirsParallel);
                }
            }
            assert(violationDirs != UVWDir::NONE);

            vector<vector<pair<HEH, UVWDir>>> pathHaDirs;
            reconstructArcPaths(pCurrent, pathHaDirs, false);

            UVWDir walkedDirs = UVWDir::NONE;
            for (int i = 0; i < (int)pathHaDirs.size(); i++)
                for (auto& [ha, dir] : pathHaDirs[i])
                    walkedDirs = walkedDirs | dir;
            UVWDir monotonousDirs = walkedDirs & ~(-walkedDirs);
            // assuming there must be at least one monotonous direction on an arcpath around the boundary of the blocks
            if (monotonousDirs == UVWDir::NONE)
            {
                DLOG(WARNING) << "Could not find arc path that is monotonous, skipping...";
                continue;
            }

            Vec3Q deltaUVW = pCurrent.u - pCurrent.u0;
            UVWDir signGivingDir = toDir(deltaUVW);
            // Now craft constraint and only use edges along nonParallelDirs
            for (int i = 0; i < (int)pathHaDirs.size(); i++)
                for (auto& [ha, dir] : pathHaDirs[i])
                {
                    // 1st filter: is non-parallel dir?
                    if ((dir & nonParallelDirs) == UVWDir::NONE)
                        continue;
                    if ((-dir & monotonousDirs) != UVWDir::NONE)
                        throw std::logic_error("Misconception about directions...");
                    if ((dir & monotonousDirs) != UVWDir::NONE)
                        a2signFailsafe[mcMesh.edge_handle(ha)] += 1;

                    if ((dir & signGivingDir) != UVWDir::NONE)
                        a2sign[mcMesh.edge_handle(ha)] += 1;
                    else if ((-dir & signGivingDir) != UVWDir::NONE)
                        a2sign[mcMesh.edge_handle(ha)] -= 1;
                    else // Signs dont matter, POS gets +, NEG gets -
                        a2sign[mcMesh.edge_handle(ha)] += (isNeg(dir) ? -1 : 1);
                }
        }

        nonZeroSumArcs.emplace_back();
        failsafeNonZeroSumArcs.emplace_back();
        vector<pair<int, EH>>& nonZeroSumConstraint = nonZeroSumArcs.back();
        vector<pair<int, EH>>& failsafeNonZeroSumConstraint = failsafeNonZeroSumArcs.back();

        for (auto& [aa, sign] : a2signFailsafe)
        {
            if (sign == 0)
                DLOG(WARNING) << "Tangled arcpath on block boundary, should still be okay";
            else
                failsafeNonZeroSumConstraint.push_back({sign, aa});
        }
        nonZeroSumConstraint = failsafeNonZeroSumConstraint;
        return true;
    }

    return false;
}

bool StructurePreserver::findIncontractibleZeroLoops(const IntegerGridGraphNode& c0,
                                                     const CH& bRef0,
                                                     const Vec3i& i0Min,
                                                     const Vec3i& i0Max,
                                                     set<VH>& nsVisited,
                                                     vector<vector<pair<int, EH>>>& nonZeroSumArcs,
                                                     vector<vector<pair<int, EH>>>& failsafeNonZeroSumArcs) const
{
    auto& mcMesh = mcMeshProps().mesh();
    using PathQueue = std::
        priority_queue<IntegerGridGraphPath, std::deque<IntegerGridGraphPath>, GreaterIntegerGridGraphPathCompare>;

    struct Point
    {
        IntegerGridGraphNode c;
        CH bRef;
        Vec3i igm;

        bool operator==(const Point& other) const
        {
            return other.c == c && other.bRef == bRef && other.igm == igm;
        }
    };
    struct GreaterPointCompare
    {
        bool operator()(const Point& p1, const Point& p2) const
        {
            Vec3i v1 = Vec3i(p1.c.n.idx(), p1.c.a.idx(), p1.c.p.idx());
            Vec3i v2 = Vec3i(p2.c.n.idx(), p2.c.a.idx(), p2.c.p.idx());
            if (v1 != v2)
                return v1 < v2;
            if (p1.bRef.idx() != p2.bRef.idx())
                return p1.bRef.idx() < p2.bRef.idx();
            return p1.igm < p2.igm;
        }
    };
    struct Link
    {
        IntegerGridGraphNode from;
        IntegerGridGraphNode to;
        IntegerGridGraphNode via;
        CH bRef;
        Vec3i igm;

        bool operator==(const Link& other) const
        {
            return other.from == from && other.to == to && other.via == via && other.bRef == bRef && other.igm == igm;
        }
    };
    struct GreaterLinkCompare
    {
        bool operator()(const Link& l1, const Link& l2) const
        {
            Vec3i v1 = Vec3i(l1.from.n.idx(), l1.from.a.idx(), l1.from.p.idx());
            Vec3i v2 = Vec3i(l2.from.n.idx(), l2.from.a.idx(), l2.from.p.idx());
            if (v1 != v2)
                return v1 < v2;
            v1 = Vec3i(l1.to.n.idx(), l1.to.a.idx(), l1.to.p.idx());
            v2 = Vec3i(l2.to.n.idx(), l2.to.a.idx(), l2.to.p.idx());
            if (v1 != v2)
                return v1 < v2;
            v1 = Vec3i(l1.via.n.idx(), l1.via.a.idx(), l1.via.p.idx());
            v2 = Vec3i(l2.via.n.idx(), l2.via.a.idx(), l2.via.p.idx());
            if (v1 != v2)
                return v1 < v2;
            if (l1.bRef.idx() != l2.bRef.idx())
                return l1.bRef.idx() < l2.bRef.idx();
            return l1.igm < l2.igm;
        }
    };
    set<Link, GreaterLinkCompare> linksVisited;
    set<Link, GreaterLinkCompare> linksLoop;

    map<Point, IntegerGridGraphPath, GreaterPointCompare> point2path;

    int nTraversed = 0;
    int nLoop = 0;
    {
        PathQueue pQ;
        IntegerGridGraphPath pStart = initializePathStart(c0, bRef0, i0Min, i0Max, true);
        pQ.push(pStart);

        while (!pQ.empty())
        {
            IntegerGridGraphPath pIn = pQ.top();
            pQ.pop();

            // Tree/Cotree decomp
            {
                map<CH, vector<Transition>> bRef2trans;
                IntegerGridGraphNode cStart = pIn.precursors.back().cPrev;
                IntegerGridGraphNode via;

                list<Link> linksToAdd;
                if (cStart.is_valid() && pIn.c != cStart)
                {
                    for (VH n : cStart.nodes(mcMesh))
                        if (!contains(mcMesh.cell_vertices(pIn.bRef), n))
                            throw std::logic_error("Bad Path");
                    set<VH> nodes;
                    for (VH n : pIn.c.nodes(mcMesh))
                        nodes.insert(n);
                    for (VH n : cStart.nodes(mcMesh))
                        nodes.insert(n);
                    if (nodes.size() == 2)
                    {
                        if (!cStart.n.is_valid() || !pIn.c.n.is_valid())
                            throw std::logic_error("Unexpected node count");
                        for (EH a : mcMesh.cell_edges(pIn.bRef))
                            if (contains(mcMesh.edge_vertices(a), *nodes.begin())
                                && contains(mcMesh.edge_vertices(a), *(++nodes.begin())))
                            {
                                via.a = a;
                                bRef2trans = _ba2transitionsAroundArc.at({pIn.bRef, a});
                                break;
                            }
                    }
                    if (bRef2trans.empty())
                    {
                        for (FH p : mcMesh.cell_faces(pIn.bRef))
                        {
                            bool containsAll = true;
                            for (VH n : nodes)
                                if (!contains(mcMesh.face_vertices(p), n))
                                {
                                    containsAll = false;
                                    break;
                                }
                            if (containsAll)
                            {
                                via.p = p;
                                bRef2trans[pIn.bRef] = {Transition()};
                                HFH hp = mcMesh.halfface_handle(p, 0);
                                if (mcMesh.incident_cell(hp) == pIn.bRef)
                                    hp = mcMesh.opposite_halfface_handle(hp);
                                if (!mcMesh.is_boundary(hp))
                                    bRef2trans[mcMesh.incident_cell(hp)]
                                        = {mcMeshProps().hpTransition<PATCH_TRANSITION>(
                                            mcMesh.opposite_halfface_handle(hp))};
                                break;
                            }
                        }
                    }
                    if (bRef2trans.empty())
                        bRef2trans[pIn.bRef] = {Transition()};

                    Link link1 = {cStart, pIn.c, via, pIn.bRef, pIn.i};
                    Link link2 = {pIn.c, cStart, via, pIn.bRef, pIn.i};
                    if (!linksVisited.count(link1))
                    {
                        if (!bRef2trans.empty())
                            for (const auto& [bRef, transitions] : bRef2trans)
                            {
                                link1.bRef = link2.bRef = bRef;
                                if (bRef == pIn.bRef)
                                {
                                    link1.igm = link2.igm = pIn.i;
                                }
                                else
                                {
                                    VH nRef = via.nodes(mcMesh).front();
                                    auto& trans = transitions.front();
                                    Vec3i shift = _b2n2igm[bRef.idx()].at(nRef)
                                                  - trans.rotate(_b2n2igm[pIn.bRef.idx()].at(nRef));
                                    Vec3i iLocal = trans.rotate(pIn.i) + shift;
                                    link1.igm = link2.igm = iLocal;
                                }
                                assert(!linksVisited.count(link1) && !linksVisited.count(link2));
                                linksToAdd.push_back(link1);
                                linksToAdd.push_back(link2);
                            }
                        if (!point2path.count({pIn.c, pIn.bRef, pIn.i}))
                            nTraversed++;
                        else
                        {
                            if (!linksLoop.count(link1))
                                nLoop++;
                        }
                    }
                }
                if (point2path.count({pIn.c, pIn.bRef, pIn.i}))
                {
#ifndef NDEBUG
                    for (auto& l : linksToAdd)
                        assert(!linksVisited.count(l));
#endif
                    linksLoop.insert(linksToAdd.begin(), linksToAdd.end());
                }
                else
                {
#ifndef NDEBUG
                    for (auto& l : linksToAdd)
                        assert(!linksLoop.count(l));
#endif
                    linksVisited.insert(linksToAdd.begin(), linksToAdd.end());
                }
            }

            // Keep this as is
            if (point2path.count({pIn.c, pIn.bRef, pIn.i}))
                continue;

            if (pIn.c.n.is_valid())
                nsVisited.insert(pIn.c.n);

            map<CH, vector<Transition>> bRef2trans;
            fetchTransitionsAroundEntity(pIn.c, pIn.bRef, bRef2trans);
            VH nRef = pIn.c.nodes(mcMesh).front();

            for (const auto& [bRef, transitions] : bRef2trans)
            {
                auto& trans = transitions.front();
                // for (const Transition& trans : transitions)
                {
                    Vec3i shift = _b2n2igm[bRef.idx()].at(nRef) - trans.rotate(_b2n2igm[pIn.bRef.idx()].at(nRef));

                    Vec3i iLocal = trans.rotate(pIn.i) + shift;
                    Vec3i iMinLocal = trans.rotate(pIn.iMin) + shift;
                    Vec3i iMaxLocal = trans.rotate(pIn.iMax) + shift;
                    Vec3Q uLocal = trans.apply(pIn.u);
                    Vec3Q u0Local = trans.apply(pIn.u0);
                    if (pIn.c.n.is_valid())
                        if (uLocal != _b2n2uvw[bRef.idx()].at(pIn.c.n))
                            throw std::logic_error("What?");
                    Transition transTotal = pIn.trans.chain(trans);

                    point2path[{pIn.c, bRef, iLocal}] = pIn;

                    map<IntegerGridGraphNode, pairTT<Vec3i>> csOut2interval;
                    auto& p2igmMin = _b2p2igmMin[bRef.idx()];
                    auto& p2igmMax = _b2p2igmMax[bRef.idx()];
                    auto& a2igmMin = _b2a2igmMin[bRef.idx()];
                    auto& a2igmMax = _b2a2igmMax[bRef.idx()];
                    auto& n2IGM = _b2n2igm[bRef.idx()];
                    for (HFH hp : mcMesh.cell_halffaces(bRef))
                    {
                        FH p = mcMesh.face_handle(hp);
                        Vec3i pMin = p2igmMin.at(p);
                        Vec3i pMax = p2igmMax.at(p);
                        if (!bboxOverlap(iMinLocal, iMaxLocal, pMin, pMax))
                            continue;
                        UVWDir patchPlusDir = ~halfpatchNormalDir(hp) & UVWDir::POS_U_POS_V_POS_W;
                        pairTT<Vec3i> patchInnerOverlap = {pMin + toVec(patchPlusDir), pMax - toVec(patchPlusDir)};
                        if (bboxOverlap(iMinLocal, iMaxLocal, patchInnerOverlap.first, patchInnerOverlap.second))
                        {
                            for (int i = 0; i < 3; i++)
                            {
                                patchInnerOverlap.first[i] = std::max(patchInnerOverlap.first[i], iMinLocal[i]);
                                patchInnerOverlap.second[i] = std::min(patchInnerOverlap.second[i], iMaxLocal[i]);
                            }
                            csOut2interval[p] = patchInnerOverlap;
                        }
                        else
                            csOut2interval[p] = {Vec3i(1, 1, 1), Vec3i(-1, -1, -1)};

                        for (HEH ha : mcMesh.halfface_halfedges(hp))
                        {
                            EH a = mcMesh.edge_handle(ha);
                            if (csOut2interval.count(a))
                                continue;
                            Vec3i aMin = a2igmMin.at(a);
                            Vec3i aMax = a2igmMax.at(a);
                            if (!bboxOverlap(iMinLocal, iMaxLocal, aMin, aMax))
                                continue;
                            UVWDir arcPlusDir = halfarcDirInBlock(ha, bRef);
                            arcPlusDir = (arcPlusDir | -arcPlusDir) & UVWDir::POS_U_POS_V_POS_W;
                            pairTT<Vec3i> arcInnerOverlap = {aMin + toVec(arcPlusDir), aMax - toVec(arcPlusDir)};
                            if (bboxOverlap(iMinLocal, iMaxLocal, arcInnerOverlap.first, arcInnerOverlap.second))
                            {
                                for (int i = 0; i < 3; i++)
                                {
                                    arcInnerOverlap.first[i] = std::max(arcInnerOverlap.first[i], iMinLocal[i]);
                                    arcInnerOverlap.second[i] = std::min(arcInnerOverlap.second[i], iMaxLocal[i]);
                                }
                                csOut2interval[a] = arcInnerOverlap;
                            }
                            else
                                csOut2interval[a] = {Vec3i(1, 1, 1), Vec3i(-1, -1, -1)};

                            for (VH n : mcMesh.halfedge_vertices(ha))
                            {
                                Vec3i igm = n2IGM.at(n);
                                if (bboxOverlap(iMinLocal, iMaxLocal, igm, igm))
                                    csOut2interval[n] = {igm, igm};
                            }
                        }
                    }

                    for (const auto& kv : csOut2interval)
                    {
                        auto& cOut = kv.first;
                        auto& range = kv.second;
                        if (range.first[0] > range.second[0])
                            continue;

                        auto [nOut, minSqDist] = findClosestNode(cOut, bRef, uLocal);
                        Vec3Q uOut = _b2n2uvw[bRef.idx()].at(nOut);
                        Vec3Q delta = uOut - uLocal;
                        double lengthOut = pIn.length + Vec3Q2d(delta).norm();
                        Precursor precursorOut = {};
                        precursorOut.cPrev = pIn.c;
                        precursorOut.transPrev = pIn.trans;
                        precursorOut.bRefPrev = pIn.bRef;
                        precursorOut.deltaToPrev = -delta;
                        precursorOut.transToPrev = trans.invert();

                        IntegerGridGraphPath pOut;
                        pOut.i = iLocal;
                        pOut.iMin = range.first;
                        pOut.iMax = range.second;
                        pOut.u = uOut;
                        pOut.u0 = u0Local;
                        pOut.c = cOut;
                        pOut.bRef = bRef;
                        pOut.trans = transTotal;
                        pOut.precursors = pIn.precursors;
                        pOut.precursors.push_back(precursorOut);
                        pOut.hops = pIn.hops + 1;
                        pOut.length = lengthOut;
                        pQ.push(pOut);
                    }
                }
            }
        }
    }

    for (const auto& link : linksVisited)
        if (linksLoop.count(link))
            throw std::logic_error("Link both visited and a loop?");

    DLOG(INFO) << "So far so good, traversed links: " << nTraversed << "(" << linksVisited.size()
               << "), looplinks: " << nLoop << "(" << linksLoop.size() << ")";

    map<Point, list<Link>, GreaterPointCompare> point2linksVisitedOut;
    map<Point, list<Link>, GreaterPointCompare> point2linksLoopOut;
    for (const auto& link : linksVisited)
    {
        Point pFrom = {link.from, link.bRef, link.igm};
        point2linksVisitedOut[pFrom].push_back(link);
    }
    for (const auto& link : linksLoop)
    {
        Point pFrom = {link.from, link.bRef, link.igm};
        point2linksLoopOut[pFrom].push_back(link);
    }

    map<CH, list<list<Point>>> bRef2clusters;

    set<Point, GreaterPointCompare> pointsVisited;
    for (auto& [p0, linksOut] : point2linksVisitedOut)
    {
        if (pointsVisited.count(p0))
            continue;
        list<Point> pointQ;
        pointQ.push_back(p0);
        bRef2clusters[p0.bRef].emplace_back();
        bRef2clusters[p0.bRef].back().push_back(p0);

        while (!pointQ.empty())
        {
            Point p = pointQ.front();
            pointQ.pop_front();
            for (const auto& link : point2linksVisitedOut[p])
            {
                if (link.bRef != p0.bRef || link.igm != p0.igm)
                    continue;
                Point pTo = {link.to, link.bRef, link.igm};
                if (!pointsVisited.count(pTo))
                {
                    pointsVisited.insert(pTo);
                    bRef2clusters[p0.bRef].back().push_back(pTo);
                    pointQ.push_back(pTo);
                }
            }
        }
    }

    list<Link> loopLinksToProcess(linksLoop.begin(), linksLoop.end());
    while (!loopLinksToProcess.empty())
    {
        // Check if stays within one cluster
        Link linkLoop = loopLinksToProcess.front();
        loopLinksToProcess.pop_front();
        if (!linksLoop.count(linkLoop))
            continue;

        bool sameCluster = false;
        {
            Point pFrom = {linkLoop.from, linkLoop.bRef, linkLoop.igm};
            Point pTo = {linkLoop.to, linkLoop.bRef, linkLoop.igm};
            for (const auto& cluster : bRef2clusters[linkLoop.bRef])
                if (std::find(cluster.begin(), cluster.end(), pFrom) != cluster.end()
                    && std::find(cluster.begin(), cluster.end(), pTo) != cluster.end())
                {
                    sameCluster = true;
                    break;
                }
        }
        if (sameCluster)
        {
            list<Link> linkClones;
            {
                map<CH, vector<Transition>> bRef2trans;
                if (linkLoop.via.a.is_valid())
                    bRef2trans = _ba2transitionsAroundArc.at({linkLoop.bRef, linkLoop.via.a});
                else if (linkLoop.via.p.is_valid())
                {
                    bRef2trans[linkLoop.bRef] = {Transition()};
                    HFH hp = mcMesh.halfface_handle(linkLoop.via.p, 0);
                    if (mcMesh.incident_cell(hp) == linkLoop.bRef)
                        hp = mcMesh.opposite_halfface_handle(hp);
                    if (!mcMesh.is_boundary(hp))
                        bRef2trans[mcMesh.incident_cell(hp)]
                            = {mcMeshProps().hpTransition<PATCH_TRANSITION>(mcMesh.opposite_halfface_handle(hp))};
                }
                else
                    bRef2trans[linkLoop.bRef] = {Transition()};

                Link linkClone1 = {linkLoop.from, linkLoop.to, linkLoop.via, linkLoop.bRef, linkLoop.igm};
                Link linkClone2 = {linkLoop.to, linkLoop.from, linkLoop.via, linkLoop.bRef, linkLoop.igm};

                for (const auto& [bRef, transitions] : bRef2trans)
                {
                    linkClone1.bRef = linkClone2.bRef = bRef;
                    if (bRef == linkLoop.bRef)
                    {
                        linkClone1.igm = linkClone2.igm = linkLoop.igm;
                    }
                    else
                    {
                        VH nRef = linkLoop.via.nodes(mcMesh).front();
                        auto& trans = transitions.front();
                        Vec3i shift
                            = _b2n2igm[bRef.idx()].at(nRef) - trans.rotate(_b2n2igm[linkLoop.bRef.idx()].at(nRef));
                        Vec3i iLocal = trans.rotate(linkLoop.igm) + shift;
                        linkClone1.igm = linkClone2.igm = iLocal;
                    }
                    linkClones.push_back(linkClone1);
                    linkClones.push_back(linkClone2);
                }
            }
            for (const auto& l : linkClones)
            {
                assert(!linksVisited.count(l) && linksLoop.count(l));
                linksLoop.erase(l);
                linksVisited.insert(l);
                Point pFrom = {l.from, l.bRef, l.igm};
                Point pTo = {l.to, l.bRef, l.igm};
                point2linksVisitedOut[pFrom].remove(l);
                point2linksLoopOut[pFrom].remove(l);

                auto& clusters = bRef2clusters[l.bRef];
                auto it1 = clusters.end();
                auto it2 = clusters.end();
                for (auto it = clusters.begin(); it != clusters.end(); it++)
                {
                    const auto& cluster = *it;
                    if (std::find(cluster.begin(), cluster.end(), pFrom) != cluster.end())
                        it1 = it;
                    if (std::find(cluster.begin(), cluster.end(), pTo) != cluster.end())
                        it2 = it;
                    if (it1 != clusters.end() && it2 != clusters.end())
                        break;
                }
                if (it1 == clusters.end() && it2 == clusters.end())
                {
                    clusters.emplace_back();
                    clusters.back().push_back(pFrom);
                    clusters.back().push_back(pTo);

                    for (auto p : {pFrom, pTo})
                        for (const auto& lNext : point2linksLoopOut[p])
                            if (lNext.bRef == l.bRef && lNext.igm == l.igm)
                                loopLinksToProcess.push_back(lNext);
                }
                else if (it1 == clusters.end())
                {
                    // Add pFrom to it2
                    auto& cluster2 = *it2;
                    cluster2.push_back(pFrom);
                    for (const auto& p : cluster2)
                        for (const auto& lNext : point2linksLoopOut[p])
                            if (lNext.bRef == l.bRef && lNext.igm == l.igm)
                                loopLinksToProcess.push_back(lNext);
                }
                else if (it2 == clusters.end())
                {
                    // Add pTo to it1
                    auto& cluster1 = *it1;
                    cluster1.push_back(pTo);
                    for (const auto& p : cluster1)
                        for (const auto& lNext : point2linksLoopOut[p])
                            if (lNext.bRef == l.bRef && lNext.igm == l.igm)
                                loopLinksToProcess.push_back(lNext);
                }
                else if (it1 != it2)
                {
                    // Merge clusters
                    auto& cluster1 = *it1;
                    auto& cluster2 = *it2;
                    for (const auto& p : cluster2)
                        cluster1.push_back(p);
                    clusters.erase(it2);
                    for (const auto& p : cluster1)
                        for (const auto& lNext : point2linksLoopOut[p])
                            if (lNext.bRef == l.bRef && lNext.igm == l.igm)
                                loopLinksToProcess.push_back(lNext);
                }
            }
        }
        else
            continue;
    }

#ifndef NDEBUG
    for (const auto& link : linksVisited)
        assert(!linksLoop.count(link));
#endif
    DLOG(INFO) << "Even better, remaining cotree links " << linksLoop.size() << ", tree links " << linksVisited.size();
    list<Link> remainingLinks(linksLoop.begin(), linksLoop.end());
    while (!remainingLinks.empty())
    {
        // Check if stays within one cluster
        Link linkLoop = remainingLinks.front();
        remainingLinks.pop_front();

        bool sameCluster = false;
        {
            Point pFrom = {linkLoop.from, linkLoop.bRef, linkLoop.igm};
            Point pTo = {linkLoop.to, linkLoop.bRef, linkLoop.igm};
            for (const auto& cluster : bRef2clusters[linkLoop.bRef])
                if (std::find(cluster.begin(), cluster.end(), pFrom) != cluster.end()
                    && std::find(cluster.begin(), cluster.end(), pTo) != cluster.end())
                {
                    sameCluster = true;
                    break;
                }
        }

        if (sameCluster)
            throw std::logic_error("Link should have been processed earlier");
    }

    list<set<EH>> constraintArcSets;
    for (auto& linkLoop : linksLoop)
    {
        Point pFrom = {linkLoop.from, linkLoop.bRef, linkLoop.igm};
        Point pTo = {linkLoop.to, linkLoop.bRef, linkLoop.igm};

        IntegerGridGraphPath& pFromPath = point2path.at(pFrom);
        IntegerGridGraphPath& pToPath = point2path.at(pTo);

        IntegerGridGraphPath pLink;
        // Link start
        {
            pLink.i = linkLoop.igm;
            pLink.iMin = linkLoop.igm;
            pLink.iMax = linkLoop.igm;
            pLink.u0 = _b2n2uvw[linkLoop.bRef.idx()].at(linkLoop.from.nodes(mcMesh).front());
            pLink.bRef = linkLoop.bRef;
            pLink.trans = Transition();
            Precursor pre = {};
            pLink.precursors = {pre};
        }
        // Link end
        {
            auto [nOut, minSqDist] = findClosestNode(linkLoop.to, pLink.bRef, pLink.u0);
            Vec3Q uOut = _b2n2uvw[pLink.bRef.idx()].at(nOut);
            Vec3Q delta = uOut - pLink.u0;
            double lengthOut = Vec3Q2d(delta).norm();
            Precursor precursorOut = {};
            precursorOut.cPrev = linkLoop.from;
            precursorOut.transPrev = Transition();
            precursorOut.bRefPrev = pLink.bRef;
            precursorOut.deltaToPrev = -delta;
            precursorOut.transToPrev = Transition();

            pLink.u = uOut;
            pLink.c = linkLoop.to;
            pLink.precursors.push_back(precursorOut);
            pLink.hops = 1;
            pLink.length = lengthOut;
        }

        set<EH> as;
        set<EH> asRemove;
        for (IntegerGridGraphPath* pPtr : {&pFromPath, &pToPath, &pLink})
        {
            vector<vector<pair<HEH, UVWDir>>> pathHaDirs;
            reconstructArcPaths(*pPtr, pathHaDirs, true);

            for (int i = 0; i < (int)pathHaDirs.size(); i++)
                for (auto& [ha, dir] : pathHaDirs[i])
                {
                    if (as.count(mcMesh.edge_handle(ha)))
                        asRemove.insert(mcMesh.edge_handle(ha));
                    else
                        as.insert(mcMesh.edge_handle(ha));
                }
        }
        // Dirty workaround for something that went wrong earlier
        if (as.empty())
            continue;

        std::set<EH> asFinal;
        std::set_difference(
            as.begin(), as.end(), asRemove.begin(), asRemove.end(), std::inserter(asFinal, asFinal.end()));

        if (asFinal.empty())
            continue;

        for (auto it = constraintArcSets.begin(); it != constraintArcSets.end();)
        {
            set<EH> intersection;
            std::set_intersection(it->begin(),
                                  it->end(),
                                  asFinal.begin(),
                                  asFinal.end(),
                                  std::inserter(intersection, intersection.begin()));
            if (intersection.size() == asFinal.size())
            {
                it = constraintArcSets.erase(it);
            }
            else if (intersection.size() == it->size())
            {
                asFinal.clear();
                break;
            }
            else
                ++it;
        }
        if (!asFinal.empty())
            constraintArcSets.push_back(asFinal);
    }
    DLOG(INFO) << "Final constraint arc sets: " << constraintArcSets.size();
    for (const auto& as : constraintArcSets)
    {
        failsafeNonZeroSumArcs.emplace_back();
        auto& failsafeConstraint = failsafeNonZeroSumArcs.back();
        for (EH a : as)
            failsafeConstraint.push_back({1, a});
        nonZeroSumArcs.emplace_back();
        auto& normalConstraint = nonZeroSumArcs.back();
        normalConstraint = failsafeConstraint;
    }

    return !linksLoop.empty();
}

void StructurePreserver::fetchTransitionsAroundEntity(const IntegerGridGraphNode& c,
                                                      const CH& bRef,
                                                      map<CH, vector<Transition>>& bRef2trans) const
{
    auto& mcMesh = mcMeshProps().mesh();
    if (c.n.is_valid())
        bRef2trans = _bn2transitionsAroundNode.at({bRef, c.n});
    else if (c.a.is_valid())
        bRef2trans = _ba2transitionsAroundArc.at({bRef, c.a});
    else // c.p.is_valid()
    {
        bRef2trans[bRef] = {Transition()};
        HFH hp = mcMesh.halfface_handle(c.p, 0);
        if (mcMesh.incident_cell(hp) == bRef)
            hp = mcMesh.opposite_halfface_handle(hp);
        if (!mcMesh.is_boundary(hp))
            bRef2trans[mcMesh.incident_cell(hp)]
                = {mcMeshProps().hpTransition<PATCH_TRANSITION>(mcMesh.opposite_halfface_handle(hp))};
    }
}

StructurePreserver::IntegerGridGraphPath StructurePreserver::initializePathStart(const IntegerGridGraphNode& c0,
                                                                                 const CH& bRef0,
                                                                                 const Vec3i& i0Min,
                                                                                 const Vec3i& i0Max,
                                                                                 bool includePrecursors) const
{
    auto& mcMesh = mcMeshProps().mesh();
    IntegerGridGraphPath pStart = {};
    pStart.i = i0Min;
    pStart.iMin = i0Min;
    pStart.iMax = i0Max;
    pStart.u = _b2n2uvw[bRef0.idx()].at(c0.nodes(mcMesh).front());
    pStart.u0 = pStart.u;
    pStart.c = c0;
    pStart.bRef = bRef0;
    pStart.trans = Transition();
    pStart.hops = 0;
    pStart.length = 0.0;
    if (includePrecursors)
    {
        Precursor pre = {};
        pStart.precursors = {pre};
    }
    return pStart;
}

std::pair<VH, double>
StructurePreserver::findClosestNode(const IntegerGridGraphNode& c, const CH& bRef, const Vec3Q& uvw) const
{
    auto& mcMesh = mcMeshProps().mesh();
    double minSqDist = DBL_MAX;
    VH closestNode;
    for (VH n : c.nodes(mcMesh))
    {
        Vec3Q uvwC = _b2n2uvw[bRef.idx()].at(n);
        double sqDist = Vec3Q2d(uvw - uvwC).sqrnorm();
        if (sqDist < minSqDist)
        {
            minSqDist = sqDist;
            closestNode = n;
        }
    }
    return {closestNode, minSqDist};
}

bool StructurePreserver::checkArcOverlap(const WeaklyMonotonousPath& pathCurrent,
                                         const HEH& ha,
                                         const CH& bRef,
                                         const Transition& trans) const
{
    auto& mcMesh = mcMeshProps().mesh();

    UVWDir dir2 = trans.invert().rotate(halfarcDirInBlock(ha, bRef));
    int arcLen = (_considerDoubleIntLength ? 2 : 1) * mcMeshProps().get<ARC_INT_LENGTH>(mcMesh.edge_handle(ha));

    Vec3i deltaMin2 = isNeg(dir2) ? pathCurrent.delta + arcLen * toVec(dir2) : pathCurrent.delta;
    Vec3i deltaMax2 = isNeg(dir2) ? pathCurrent.delta : pathCurrent.delta + arcLen * toVec(dir2);

    // Check for overlap
    return bboxOverlap(pathCurrent.deltaMin, pathCurrent.deltaMax, deltaMin2, deltaMax2);
}

bool StructurePreserver::checkPatchOverlap(const WeaklyMonotonousPath& pathCurrent,
                                           const HFH& hp,
                                           const Transition& trans) const
{
    auto& mcMesh = mcMeshProps().mesh();

    CH bRef2 = mcMesh.incident_cell(hp);

    Vec3i deltaMin2 = pathCurrent.delta;
    Vec3i deltaMax2 = pathCurrent.delta;

    HEH ha1 = findMatching(mcMesh.halfface_halfedges(hp),
                           [&](const HEH& ha) { return mcMesh.from_vertex_handle(ha) == pathCurrent.n; });
    assert(ha1.is_valid());
    HEH haCurr = ha1;
    Vec3i deltaCurr = pathCurrent.delta;

    Transition transInv = trans.invert();
    do
    {
        UVWDir dir = transInv.rotate(halfarcDirInBlock(haCurr, bRef2));
        EH aCurr = mcMesh.edge_handle(haCurr);
        int lengthHa = (_considerDoubleIntLength ? 2 : 1) * mcMeshProps().get<ARC_INT_LENGTH>(aCurr);
        deltaCurr += lengthHa * toVec(dir);
        for (int coord = 0; coord < 3; coord++)
        {
            deltaMax2[coord] = std::max(deltaMax2[coord], deltaCurr[coord]);
            deltaMin2[coord] = std::min(deltaMin2[coord], deltaCurr[coord]);
        }
        haCurr = mcMesh.next_halfedge_in_halfface(haCurr, hp);
    } while (haCurr != ha1);

    assert(deltaCurr == pathCurrent.delta);
    assert(deltaMin2[0] <= deltaMax2[0] && deltaMin2[1] <= deltaMax2[1] && deltaMin2[2] <= deltaMax2[2]);

    return bboxOverlap(pathCurrent.deltaMin, pathCurrent.deltaMax, deltaMin2, deltaMax2);
}

bool StructurePreserver::checkP0Containment(const WeaklyMonotonousPath& pathCurrent) const
{
    auto& mcMesh = mcMeshProps().mesh();

    CH bRef = pathCurrent.bRefCurrent;

    // Find the closest corner of bRef

    // Graph search
    list<pair<VH, Vec3i>> nQ;
    nQ.push_back({pathCurrent.n, Vec3i(0, 0, 0)});
    map<VH, Vec3i> n2displacement({{pathCurrent.n, Vec3i(0, 0, 0)}});
    while (!nQ.empty())
    {
        VH n = nQ.front().first;
        Vec3i displacement = nQ.front().second;
        nQ.pop_front();

        for (HEH ha : mcMesh.outgoing_halfedges(n))
        {
            bool inBRef = false;
            for (CH b : mcMesh.halfedge_cells(ha))
                if (b == bRef)
                {
                    inBRef = true;
                    break;
                }
            if (inBRef)
            {
                VH nTo = mcMesh.to_vertex_handle(ha);
                if (n2displacement.find(nTo) == n2displacement.end())
                {
                    UVWDir dirHa = pathCurrent.transCurrent.invert().rotate(halfarcDirInBlock(ha, bRef));
                    double length = (_considerDoubleIntLength ? 2 : 1)
                                    * mcMeshProps().get<ARC_INT_LENGTH>(mcMesh.edge_handle(ha));
                    Vec3i newDisplacement = displacement + toVec(dirHa) * length;
                    n2displacement[nTo] = newDisplacement;
                    nQ.push_back({nTo, newDisplacement});
                }
            }
        }
    }

    auto& dir2n = mcMeshProps().ref<BLOCK_CORNER_NODES>(bRef);
    VH nMin = dir2n.at(UVWDir::NEG_U_NEG_V_NEG_W);
    VH nMax = dir2n.at(UVWDir::POS_U_POS_V_POS_W);

    return bboxOverlap(pathCurrent.deltaMin,
                       pathCurrent.deltaMax,
                       pathCurrent.delta + n2displacement.at(nMin),
                       pathCurrent.delta + n2displacement.at(nMax));
}

bool StructurePreserver::bboxOverlap(const Vec3i& min1, const Vec3i& max1, const Vec3i& min2, const Vec3i& max2)
{
    int coordTouchingOrOverlaps = 0;
    for (int i = 0; i < 3; i++)
    {
        if (max1[i] < min1[i] || max2[i] < min2[i])
            return false;
        // min2[i] is in [min1[i], max1[i]]
        // OR max2[i] is in [min1[i], max1[i]]
        // OR min1[i] is in [min2[i], max2[i]]
        // OR max1[i] is in [min2[i], max2[i]]
        coordTouchingOrOverlaps
            += (min1[i] <= min2[i] && min2[i] <= max1[i]) || (min1[i] <= max2[i] && max2[i] <= max1[i])
               || (min2[i] <= min1[i] && min1[i] <= max2[i]) || (min2[i] <= max1[i] && max1[i] <= max2[i]);

        if (coordTouchingOrOverlaps == i)
            return false;
    }

    return true;
}

void StructurePreserver::debugViewConstraintPath(const vector<pair<int, EH>>& constraintPath) const
{
#ifdef MC3D_WITH_VIEWER
    auto& mcMesh = mcMeshProps().mesh();
    set<EH> asPos;
    set<EH> asNeg;
    for (auto& kv : constraintPath)
        if (kv.first > 0)
            asPos.insert(kv.second);
        else
            asNeg.insert(kv.second);
    auto& tetMesh = meshProps().mesh();
    set<HFH> restHfs;
    for (FH p : mcMesh.faces())
        for (HFH hf : mcMeshProps().ref<PATCH_MESH_HALFFACES>(p))
            for (HFH hf2 : tetMesh.face_halffaces(tetMesh.face_handle(hf)))
                restHfs.insert(hf2);
    set<EH> restEs, posEs, negEs;
    for (EH a : mcMesh.edges())
        for (HEH he : mcMeshProps().ref<ARC_MESH_HALFEDGES>(a))
        {
            if (asPos.count(a))
                posEs.insert(tetMesh.edge_handle(he));
            else if (asNeg.count(a))
                negEs.insert(tetMesh.edge_handle(he));
            else
                restEs.insert(tetMesh.edge_handle(he));
        }
    set<VH> restVs;
    for (VH n : mcMesh.vertices())
        restVs.insert(mcMeshProps().ref<NODE_MESH_VERTEX>(n));
    debugView({}, {}, restHfs, {}, negEs, posEs, restEs, {}, {}, restVs);
#endif
}

} // namespace qgp3d
