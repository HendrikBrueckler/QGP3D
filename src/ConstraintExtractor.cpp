#include "QGP3D/ConstraintExtractor.hpp"

#include "QGP3D/MCQuantizer.hpp"

#include "OpenVolumeMesh/FileManager/FileManager.hh"

namespace qgp3d
{

ConstraintExtractor::ConstraintExtractor(const TetMeshProps& meshProps)
    : TetMeshNavigator(meshProps), MCMeshNavigator(meshProps)
{
}

vector<vector<HEH>> ConstraintExtractor::getCriticalSkeletonArcs()
{
    auto& mcMesh = mcMeshProps().mesh();

    assignNodeTypes();

    // Mark cut patches.
    markCutPatches();

    vector<vector<HEH>> nodeTreePaths;
    // Mark one node per cut patch
    _isSkeletonArc = vector<bool>(mcMesh.n_edges(), false);

    // Initialize MSP with all arcs
    {
        vector<bool> nVisited(mcMesh.n_edges(), false);

        nVisited[mcMesh.v_iter()->idx()] = true;
        list<VH> nQ({*mcMesh.v_iter()});
        while (!nQ.empty())
        {
            VH n = nQ.front();
            nQ.pop_front();

            for (HEH ha : mcMesh.outgoing_halfedges(n))
            {
                VH nTo = mcMesh.to_vertex_handle(ha);
                if (!nVisited[nTo.idx()])
                {
                    _isSkeletonArc[mcMesh.edge_handle(ha).idx()] = true;
                    nQ.push_back(nTo);
                    nVisited[nTo.idx()] = true;
                }
            }
        }
    }

    // Shrink from all nodes with only one arc, if node is not singular
    list<EH> aQ;
    for (EH a : mcMesh.edges())
        for (VH n : mcMesh.edge_vertices(a))
        {
            if (countIncidentSkeletonArcs(n) == 1 && _constraintNodeType[n.idx()] == NONE)
            {
                _isSkeletonArc[a.idx()] = false;
                aQ.push_back(a);
            }
        }

    // shrink while possible
    while (!aQ.empty())
    {
        // get next element
        EH a = aQ.front();
        aQ.pop_front();
        for (VH n : mcMesh.edge_vertices(a))
            if (countIncidentSkeletonArcs(n) == 1 && _constraintNodeType[n.idx()] == NONE)
                for (EH aNext : mcMesh.vertex_edges(n))
                    if (_isSkeletonArc[aNext.idx()])
                    {
                        _isSkeletonArc[aNext.idx()] = false;
                        aQ.push_back(aNext);
                    }
    }

    // Extract paths from MSP
    {
        // TODO handle case where no leaf node exists
        VH leaf = findMatching(mcMesh.vertices(), [this](const VH& n) { return countIncidentSkeletonArcs(n) == 1; });
        if (leaf.is_valid())
        {
            // exhaustive graph search from leaf
            vector<bool> nVisited(mcMesh.n_vertices(), false);
            nVisited[leaf.idx()] = true;
            list<VH> nQ({leaf});
            list<vector<HEH>> pathQ({vector<HEH>()});
            while (!nQ.empty())
            {
                VH n = nQ.front();
                auto path = pathQ.front();
                nQ.pop_front();
                pathQ.pop_front();

                if (_constraintNodeType[n.idx()] != NONE && n != leaf)
                {
                    nodeTreePaths.emplace_back(path);
                    // Only reset the path if node is natively singular, if it is only artificial we do not constrain
                    // its position in all coordinates, so we can not cut the incoming and outgoing paths in two
                    if (_constraintNodeType[n.idx()] == NATIVE)
                        path = {};
                }

                for (HEH ha : mcMesh.outgoing_halfedges(n))
                {
                    if (!_isSkeletonArc[mcMesh.edge_handle(ha).idx()])
                        continue;
                    VH nTo = mcMesh.to_vertex_handle(ha);
                    if (!nVisited[nTo.idx()])
                    {
                        nVisited[nTo.idx()] = true;
                        nQ.push_back(nTo);
                        pathQ.emplace_back(path);
                        pathQ.back().push_back(ha);
                    }
                }
            }
        }
    }

    LOG(INFO) << nodeTreePaths.size() << " paths in spanning tree of critical nodes";

    // Get cycles constraining cut graph
    vector<vector<HEH>> cycles = getCutSurfaceCycles();
    LOG(INFO) << cycles.size() << " cycles through " << _cutSurfaces.size() << " cut surfaces determined";
    for (auto& cycle : cycles)
        nodeTreePaths.emplace_back(cycle);

    return nodeTreePaths;
}

vector<ConstraintExtractor::TetPathConstraint>
ConstraintExtractor::getTetPathConstraints(const vector<vector<HEH>>& haSequences)
{
    assert(meshProps().isAllocated<CHART_ORIG>());
    assert(meshProps().isAllocated<TRANSITION_ORIG>());
    assert(meshProps().isAllocated<IS_ORIGINAL_V>());

    auto& mcMesh = mcMeshProps().mesh();
    auto& tetMesh = meshProps().mesh();
    vector<TetPathConstraint> tetPathConstraints;

    map<CH, CH> tet2parent;
    assert(meshProps().isAllocated<CHILD_CELLS>());
    const auto& childTetProp = meshProps().prop<CHILD_CELLS>();
    for (auto kv : childTetProp)
        for (CH child : kv.second)
            tet2parent[child] = kv.first;

    map<EH, EH> e2parent;
    assert(meshProps().isAllocated<CHILD_EDGES>());
    const auto& childEdgeProp = meshProps().prop<CHILD_EDGES>();
    for (auto kv : childEdgeProp)
        for (EH child : kv.second)
            e2parent[child] = kv.first;

    map<FH, FH> f2parent;
    assert(meshProps().isAllocated<CHILD_FACES>());
    const auto& childFaceProp = meshProps().prop<CHILD_FACES>();
    for (auto kv : childFaceProp)
        for (FH child : kv.second)
            f2parent[child] = kv.first;

    for (auto& haSequence : haSequences)
    {
        list<HEH> pathHas(haSequence.begin(), haSequence.end());

        tetPathConstraints.emplace_back();
        auto& path = tetPathConstraints.back();

        auto pathHes = list<HEH>();
        for (HEH ha : pathHas)
        {
            auto hesHa = mcMeshProps().haHalfedges(ha);
            pathHes.insert(pathHes.end(), hesHa.begin(), hesHa.end());
        }

        VH nStart = mcMesh.from_vertex_handle(haSequence.front());
        VH nEnd = mcMesh.to_vertex_handle(haSequence.back());

        ConstraintNodeType startNodeType = _constraintNodeType[nStart.idx()];
        ConstraintNodeType endNodeType = _constraintNodeType[nEnd.idx()];

        assert(startNodeType != NONE && endNodeType != NONE);

        list<CH> pathTets;
        for (HEH he : pathHes)
        {
            // First tet may be any tet
            if (pathTets.empty())
            {
                pathTets.emplace_back(*tetMesh.hec_iter(he));
                continue;
            }

            VH vConn = tetMesh.from_vertex_handle(he);

            auto connectingTets = connectingTetPath(pathTets.back(), he, vConn);
            pathTets.insert(pathTets.end(), connectingTets.begin(), connectingTets.end());
        }
        assert(!pathTets.empty());

        path.vFrom = tetMesh.from_vertex_handle(pathHes.front());
        path.vTo = tetMesh.to_vertex_handle(pathHes.back());

        // If either of from/to is not original (which can happen, if from/to is an artificial singular node
        // on a cyclic singularity) snap it to one of original edges/faces corners
        // this can introduce new tets at the start/end of pathtets, that are not necessarily incident to first/last
        // path halfedge!
        determineEquivalentEndpoints(path, pathTets, e2parent, f2parent);

        // Now get the original tet sequence
        for (CH tet : pathTets)
        {
            CH tetOrig = tet;
            auto it = tet2parent.find(tetOrig);
            while (it != tet2parent.end())
            {
                tetOrig = it->second;
                it = tet2parent.find(tetOrig);
            }
            if (path.pathOrigTets.empty() || path.pathOrigTets.back() != tetOrig)
                path.pathOrigTets.emplace_back(tetOrig);
        }
        assert(!path.pathOrigTets.empty());

        vector<CH> pathCheck;
        path.offset = Vec3i(0, 0, 0);
        Vec3i ignore(0, 0, 0);
        {
            // Walk on original path, but use coordinate system given by (possibly) shifted path
            Transition transOrigCurr;

            // Walk on original path, as only that gives us the correct arc start/end distances
            auto itHe = pathHes.begin();
            auto itTet = pathTets.begin();

            // IF FROM IS ARTIFICIAL SINGULAR CYCLE SPLIT, DO NOT CONSTRAIN DISPLACEMENT ALONG THE CYCLE
            // (ONLY IF WE CONSTRAIN THE CYCLE ITSELF)
            if (startNodeType == ARTIFICIAL_ON_LINK && nStart != nEnd)
                ignore[getCycleCoord(*itTet, tetMesh.from_vertex_handle(*itHe), transOrigCurr)] = true;
            else if (startNodeType == ARTIFICIAL_ON_BOUNDARY && nStart != nEnd)
            {
                auto boundaryCoords = getBoundaryCoords(*itTet, tetMesh.from_vertex_handle(*itHe), transOrigCurr);
                ignore[boundaryCoords.first] = true;
                ignore[boundaryCoords.second] = true;
            }

            // Walk to first he (necessary because we might have prepended some connector tets that are not incident on
            // first he)
            walkToMatchingTet(tetMesh.edge_handle(*itHe), itTet, transOrigCurr);

            for (HEH ha : pathHas)
            {
                {
                    // We are at the first halfedge of the current halfarc
                    // We can use this to get the arc direction under the orig param
                    CH tet = *itTet;
                    HEH he = *itHe;
                    auto& chartOrig = meshProps().ref<CHART_ORIG>(tet);
                    auto uvwOrig1 = chartOrig.at(meshProps().mesh().from_vertex_handle(he));
                    auto uvwOrig2 = chartOrig.at(meshProps().mesh().to_vertex_handle(he));
                    UVWDir dir = transOrigCurr.invert().rotate(toDir(uvwOrig2 - uvwOrig1));

                    if (dim(dir) != 1)
                        throw std::logic_error("invalid arc dir");

                    // DLOG(INFO) << "Path " << path.vFrom << "->" << path.vTo << ": accumulating offset "
                    //            << Vec3d(toVec(dir)) * mcMeshProps().get<ARC_INT_LENGTH>(mcMesh.edge_handle(ha))
                    //            << " of ha " << ha;
                    path.offset += toVec(dir) * mcMeshProps().get<ARC_INT_LENGTH>(mcMesh.edge_handle(ha));
                }

                // Go along pathHes until vTo is reached;
                VH nTo = mcMesh.to_vertex_handle(ha);
                VH vTo = mcMeshProps().get<NODE_MESH_VERTEX>(nTo);
                while (tetMesh.to_vertex_handle(*itHe) != vTo)
                {
                    itHe++;
                    if (itHe == pathHes.end())
                        throw std::logic_error("Unconnected edge path");
                    // Go along pathTets until itHe is reached
                    walkToMatchingTet(tetMesh.edge_handle(*itHe), itTet, transOrigCurr);
                }

                itHe++;
                if (itHe != pathHes.end())
                    walkToMatchingTet(tetMesh.edge_handle(*itHe), itTet, transOrigCurr);

                // IF END IS ARTIFICIAL SINGULAR CYCLE SPLIT, DO NOT CONSTRAIN DISPLACEMENT ALONG THE CYCLE
                // (ONLY IF WE CONSTRAIN THE CYCLE ITSELF)
                else if (endNodeType == ARTIFICIAL_ON_LINK && nStart != nEnd)
                {
                    auto itHe2 = itHe;
                    itHe2--;

                    ignore[getCycleCoord(*itTet, tetMesh.to_vertex_handle(*itHe2), transOrigCurr)] = true;
                }
                else if (startNodeType == ARTIFICIAL_ON_BOUNDARY && nStart != nEnd)
                {
                    auto boundaryCoords = getBoundaryCoords(*itTet, tetMesh.from_vertex_handle(*itHe), transOrigCurr);
                    ignore[boundaryCoords.first] = true;
                    ignore[boundaryCoords.second] = true;
                }
            }
            // This assertion does not hold true, because we might have appended connector tets.
            // However, they are irrelevant for offset as all path arc directions have already been determined.
            // assert(*itTet == pathTets.back());
            assert(itHe == pathHes.end());
        }
        for (int i = 0; i < 3; i++)
            if (ignore[i] == 1)
                path.offset[i] = INT_MAX;
    }

    return tetPathConstraints;
}

list<CH> ConstraintExtractor::connectingTetPath(const CH& tetStart, const HEH& heTarget, const VH& vConn) const
{
    assert(heTarget.is_valid());
    auto& tetMesh = meshProps().mesh();
    set<CH> tetsVisited({tetStart});
    list<list<CH>> tetQ({{tetStart}});
    bool hasHe = false;
    while (!tetQ.empty())
    {
        auto tets = tetQ.front();
        tetQ.pop_front();

        hasHe = contains(tetMesh.cell_edges(tets.back()), tetMesh.edge_handle(heTarget));
        if (hasHe)
        {
            tets.pop_front(); // First cell is tetStart
            return tets;
        }

        for (CH tetNext : tetMesh.cell_cells(tets.back()))
        {
            if (tetsVisited.find(tetNext) != tetsVisited.end() || !contains(tetMesh.cell_vertices(tetNext), vConn))
                continue;
            auto tetsPlus = tets;
            tetsPlus.emplace_back(tetNext);
            tetsVisited.insert(tetNext);
            tetQ.emplace_back(tetsPlus);
        }
    }
    assert(hasHe);
    return {};
}

list<CH> ConstraintExtractor::connectingTetPath(const CH& tetStart, const HFH& hfTarget, const VH& vConn) const
{
    assert(hfTarget.is_valid());
    auto& tetMesh = meshProps().mesh();
    set<CH> tetsVisited({tetStart});
    list<list<CH>> tetQ({{tetStart}});
    bool hasHf = false;
    while (!tetQ.empty())
    {
        auto tets = tetQ.front();
        tetQ.pop_front();

        hasHf = contains(tetMesh.cell_halffaces(tets.back()), hfTarget);
        if (hasHf)
        {
            tets.pop_front(); // First cell is tetStart
            return tets;
        }

        for (CH tetNext : tetMesh.cell_cells(tets.back()))
        {
            if (tetsVisited.find(tetNext) != tetsVisited.end() || !contains(tetMesh.cell_vertices(tetNext), vConn))
                continue;
            auto tetsPlus = tets;
            tetsPlus.emplace_back(tetNext);
            tetsVisited.insert(tetNext);
            tetQ.emplace_back(tetsPlus);
        }
    }
    assert(hasHf);
    return {};
}

Transition ConstraintExtractor::transitionAlongPath(const list<CH>& path) const
{
    if (path.size() < 2)
        return Transition();

    auto& tetMesh = meshProps().mesh();
    Transition transTotal;
    for (auto itTet = path.begin(); itTet != --path.end();)
    {
        CH tetPre = *itTet;
        itTet++;
        CH tetPost = *itTet;
        // Accumulate transition
        HFH hfBetween = findMatching(
            tetMesh.cell_halffaces(tetPre),
            [&](const HFH& hf) { return tetMesh.incident_cell(tetMesh.opposite_halfface_handle(hf)) == tetPost; });
        assert(hfBetween.is_valid());
        Transition trans = meshProps().hfTransition<TRANSITION_ORIG>(hfBetween);
        transTotal = transTotal.chain(trans);
    }
    return transTotal;
}

int ConstraintExtractor::countCutPatches(const EH& a) const
{
    auto& mcMesh = mcMeshProps().mesh();
    int n = 0;
    for (FH p : mcMesh.edge_faces(a))
        if (_isCutPatch[p.idx()])
            n++;
    return n;
}

bool ConstraintExtractor::isInActualCut(const VH& n) const
{
    auto& mcMesh = mcMeshProps().mesh();
    for (FH p : mcMesh.vertex_faces(n))
        if (_isCutPatch[p.idx()] && !mcMesh.is_boundary(p))
            return true;
    return false;
}

int ConstraintExtractor::countIncidentSkeletonArcs(const VH& n) const
{
    auto& mcMesh = mcMeshProps().mesh();
    int count = 0;
    for (EH a : mcMesh.vertex_edges(n))
        if (_isSkeletonArc[a.idx()])
            count++;
    return count;
}

void ConstraintExtractor::getCutSurfaces(vector<vector<FH>>& cutSurfaces, vector<int>& p2cutSurface)
{
    assignNodeTypes();

    markCutPatches();
    cutSurfaces = _cutSurfaces;
    p2cutSurface = _cutSurfaceID;
}

void ConstraintExtractor::assignNodeTypes()
{
    auto& mcMesh = mcMeshProps().mesh();

    vector<CriticalLink> criticalLinks;
    map<EH, int> a2criticalLinkIdx;
    map<VH, vector<int>> n2criticalLinksOut;
    map<VH, vector<int>> n2criticalLinksIn;
    getCriticalLinks(criticalLinks, a2criticalLinkIdx, n2criticalLinksOut, n2criticalLinksIn, true);

    _arcIsCritical = vector<bool>(mcMesh.n_edges(), false);
    for (auto& kv : a2criticalLinkIdx)
        _arcIsCritical[kv.first.idx()] = true;

    int nCircular = 0;
    for (auto path : criticalLinks)
        if (path.cyclic)
            nCircular++;
    LOG(INFO) << nCircular << " cyclic links out of " << criticalLinks.size();

    // Mark singular nodes (including those inserted into circular arcs)
    _constraintNodeType = vector<ConstraintNodeType>(mcMesh.n_vertices(), NONE);
    for (auto* coll : {&n2criticalLinksIn, &n2criticalLinksOut})
        for (auto& kv : *coll)
        {
            auto type = mcMeshProps().nodeType(kv.first);
            if (type.first == SingularNodeType::SINGULAR || type.second == FeatureNodeType::FEATURE
                || type.second == FeatureNodeType::SEMI_FEATURE_SINGULAR_BRANCH)
                _constraintNodeType[kv.first.idx()] = NATIVE;
            else
                _constraintNodeType[kv.first.idx()] = ARTIFICIAL_ON_LINK;
        }

    vector<BoundaryRegion> boundaryRegions;
    map<HFH, int> hp2boundaryID;
    getBoundaryRegions(boundaryRegions, hp2boundaryID);

    int nBorderless = 0;
    for (auto& region : boundaryRegions)
        if (region.boundaryHas.empty())
            nBorderless++;
    LOG(INFO) << nBorderless << " borderless boundary regions out of " << boundaryRegions.size();

    // Check for boundary regions with no critical nodes or critical arcs in them
    for (auto& region : boundaryRegions)
    {
        if (containsMatching(region.ns,
                             [&](const VH& n)
                             { return n2criticalLinksIn.count(n) != 0 || n2criticalLinksOut.count(n) != 0; }))
            continue;
        assert(!region.ns.empty());
        _constraintNodeType[region.ns.begin()->idx()] = ARTIFICIAL_ON_BOUNDARY;
    }

    int nNative = 0;
    int nArtificialOnLink = 0;
    int nArtificialOnBoundary = 0;
    for (VH n : mcMesh.vertices())
        if (_constraintNodeType[n.idx()] == NATIVE)
            nNative++;
        else if (_constraintNodeType[n.idx()] == ARTIFICIAL_ON_BOUNDARY)
            nArtificialOnBoundary++;
        else if (_constraintNodeType[n.idx()] == ARTIFICIAL_ON_LINK)
            nArtificialOnLink++;

    LOG(INFO) << nNative << " native singular nodes, " << nArtificialOnLink
              << " singular nodes inserted on cyclic singular links, " << nArtificialOnBoundary
              << " singular nodes inserted on borderless boundary region ("
              << nNative + nArtificialOnBoundary + nArtificialOnLink << " total)";
}

bool ConstraintExtractor::buildBlockSpanningTree(VH& nRoot,
                                                 CH& bRoot,
                                                 vector<CH>& cellTreePrecursor,
                                                 int cutSurfaceToAvoid)
{
    auto& mcMesh = mcMeshProps().mesh();

    cellTreePrecursor = vector<CH>(mcMesh.n_cells(), CH(-1));

    // 1. start with trivial cut (all faces are cut)
    if (cutSurfaceToAvoid < 0)
        _isCutPatch = vector<bool>(mcMesh.n_faces(), true);

    // 2. reduce cut-faces by dual spanning tree

    // Find strictly singular node
    nRoot = findMatching(mcMesh.vertices(),
                         [&](const VH& n)
                         {
                             return _constraintNodeType[n.idx()] == NATIVE
                                    && (cutSurfaceToAvoid < 0
                                        || !containsMatching(mcMesh.vertex_faces(n),
                                                             [&, this](const FH& p)
                                                             { return _cutSurfaceID[p.idx()] == cutSurfaceToAvoid; }));
                         });

    // Else, find artificial singular node
    if (!nRoot.is_valid())
        nRoot = findMatching(mcMesh.vertices(),
                             [&](const VH& n)
                             {
                                 return (_constraintNodeType[n.idx()] == ARTIFICIAL_ON_LINK
                                         || _constraintNodeType[n.idx()] == ARTIFICIAL_ON_BOUNDARY)
                                        && (cutSurfaceToAvoid < 0
                                            || !containsMatching(mcMesh.vertex_faces(n),
                                                                 [&, this](const FH& p) {
                                                                     return _cutSurfaceID[p.idx()] == cutSurfaceToAvoid;
                                                                 }));
                             });

    // Else fail
    if (!nRoot.is_valid())
        return false;

    bRoot = *mcMesh.vc_iter(nRoot);
    // add seed
    vector<bool> bVisited(mcMesh.n_cells(), false);
    bVisited[bRoot.idx()] = true;
    list<CH> bQ({bRoot});

    // grow tree until all cells in connected component visited
    while (!bQ.empty())
    {
        // get opposite of next halfface
        CH b = bQ.front();
        bQ.pop_front();

        // if b == bRoot, handle first halffaces that contain nRoot, so nRoot does not land on cut surface
        list<HFH> hps;
        for (HFH hp : mcMesh.cell_halffaces(b))
        {
            if (b == bRoot)
            {
                if (contains(mcMesh.halfface_vertices(hp), nRoot))
                    hps.push_front(hp);
                else
                    hps.push_back(hp);
            }
            else
                hps.push_back(hp);
        }

        for (HFH hp : hps)
        {
            hp = mcMesh.opposite_halfface_handle(hp);
            // only process if not on boundary
            if (!mcMesh.is_boundary(hp) && (cutSurfaceToAvoid < 0 || !_isCutPatch[mcMesh.face_handle(hp).idx()]))
            {
                // check if cell already visited
                CH bNext = mcMesh.incident_cell(hp);
                if (!bVisited[bNext.idx()])
                {
                    if (cutSurfaceToAvoid < 0)
                        _isCutPatch[mcMesh.face_handle(hp).idx()] = false; // remove face from cutgraph
                    // process cell
                    bVisited[bNext.idx()] = true;
                    bQ.push_back(bNext);
                    _cellTreePrecursor[bNext.idx()] = b;
                }
            }
        }
    }
    return true;
}

void ConstraintExtractor::markCutPatches()
{
    auto& mcMesh = mcMeshProps().mesh();

    // From AlgoHex' generate_cut_surface()
    if (!buildBlockSpanningTree(_nRoot, _bRoot, _cellTreePrecursor, -1))
        throw std::logic_error("No leaf node in spanning tree");

    {
        // 3. further shrink cut surface
        list<FH> pQ;
        // initialize queue with edges (1) not on boundary, and (2) adjacent to one cut face
        for (FH p : mcMesh.faces())
            if (containsMatching(mcMesh.face_edges(p), [&, this](const EH& a) { return countCutPatches(a) == 1; }))
            {
                pQ.push_back(p);
                _isCutPatch[p.idx()] = false;
            }

        // shrink while possible
        while (!pQ.empty())
        {
            // get next element
            FH p = pQ.front();
            pQ.pop_front();

            for (EH a : mcMesh.face_edges(p))
                if (countCutPatches(a) == 1)
                    for (FH pNext : mcMesh.edge_faces(a))
                        if (_isCutPatch[pNext.idx()])
                        {
                            _isCutPatch[pNext.idx()] = false;
                            pQ.push_back(pNext);
                        }
        }
    }

    for (EH a : mcMesh.edges())
        if (countCutPatches(a) == 1)
            throw std::logic_error("ERROR: there are shrinking candidates left but queue is empty!");

    // Remove boundary
    for (FH p : mcMesh.faces())
        if (mcMesh.is_boundary(p))
        {
            assert(_isCutPatch[p.idx()]);
            _isCutPatch[p.idx()] = false;
        }

    // 4. decompose cut surface into manifold pieces
    _cutSurfaces.clear();
    _cutSurfaceID = vector<int>(mcMesh.n_faces(), -1);
    vector<bool> pVisited(mcMesh.n_faces(), false);
    for (FH p : mcMesh.faces())
        if (_isCutPatch[p.idx()] && !pVisited[p.idx()])
        {
            // grow a new manifold patch
            _cutSurfaces.emplace_back();
            _cutSurfaces.back().push_back(p);
            _cutSurfaceID[p.idx()] = _cutSurfaces.size() - 1;
            pVisited[p.idx()] = true;

            list<FH> pQ({p});
            while (!pQ.empty())
            {
                FH pCurr = pQ.front();
                pQ.pop_front();

                for (EH a : mcMesh.face_edges(pCurr))
                {
                    if (!mcMesh.is_boundary(a) && countCutPatches(a) == 2)
                    {
                        for (FH pNext : mcMesh.edge_faces(a))
                        {
                            if (_isCutPatch[pNext.idx()] && !pVisited[pNext.idx()])
                            {
                                _cutSurfaces.back().push_back(pNext);
                                _cutSurfaceID[pNext.idx()] = _cutSurfaces.size() - 1;
                                pVisited[pNext.idx()] = true;
                                pQ.push_back(pNext);
                            }
                        }
                    }
                }
            }
        }
}

vector<HEH> ConstraintExtractor::cycleThroughCutSurface(
    const VH& nStart, const FH& pStart, const VH& nRoot, const CH& bRoot, const vector<CH>& cellTreePrecursor) const
{
    auto& mcMesh = mcMeshProps().mesh();

    vector<HEH> cycle;
    for (bool forward : {false, true})
    {
        CH b = mcMesh.incident_cell(mcMesh.halfface_handle(pStart, forward));
        VH n = nStart;
        list<HEH> halfpath;
        // Find a path connecting nStart to nRoot and going through all the cells in the sequence
        while (b != bRoot)
        {
            CH bNext = cellTreePrecursor[b.idx()];

            FH pInterface = mcMesh.face_handle(
                findMatching(mcMesh.cell_halffaces(b),
                             [&, this](const HFH& hp)
                             {
                                 return !_isCutPatch[mcMesh.face_handle(hp).idx()]
                                        && mcMesh.incident_cell(mcMesh.opposite_halfface_handle(hp)) == bNext;
                             }));
            set<VH> nsTarget;
            for (VH nTarget : mcMesh.face_vertices(pInterface))
                nsTarget.insert(nTarget);

            auto connector = pathThroughBlock(n, nsTarget, b);
            halfpath.insert(halfpath.end(), connector.begin(), connector.end());

            n = connector.empty() ? n : mcMesh.to_vertex_handle(connector.back());
            if (nsTarget.count(n) == 0)
                throw std::logic_error("Found vertex not in nsTarget");
            b = bNext;
        }
        auto connector = pathThroughBlock(n, {nRoot}, bRoot);
        halfpath.insert(halfpath.end(), connector.begin(), connector.end());

        if (forward)
            cycle.insert(cycle.end(), halfpath.begin(), halfpath.end());
        else
        {
            halfpath.reverse();
            for (HEH ha : halfpath)
                cycle.emplace_back(mcMesh.opposite_halfedge_handle(ha));
        }
    }

    if (mcMesh.from_vertex_handle(cycle.front()) != mcMesh.to_vertex_handle(cycle.back()))
        throw std::logic_error("Unclosed cycle");
    return cycle;
}

vector<vector<HEH>> ConstraintExtractor::getCutSurfaceCycles()
{
    auto& mcMesh = mcMeshProps().mesh();

    if (_cellTreePrecursor[_bRoot.idx()].is_valid())
        throw std::logic_error("Root has precursor");

    vector<vector<HEH>> cycles;
    for (auto& cutSurface : _cutSurfaces)
    {
        cycles.emplace_back();
        auto& cycle = cycles.back();
        // Get any face and trace path back to root in both directions
        if (cutSurface.empty())
            throw std::logic_error("Empty cut surface");
        FH pStart = cutSurface.front();
        VH nStart = *mcMesh.fv_iter(pStart);

        bool nRootOnCutSurface = containsSomeOf(cutSurface, mcMesh.vertex_faces(_nRoot));
        if (!nRootOnCutSurface)
            cycle = cycleThroughCutSurface(nStart, pStart, _nRoot, _bRoot, _cellTreePrecursor);
        else
        {
            // Find ANOTHER singular node, suitable as root node
            vector<CH> cellTreePrecursorAlt(mcMesh.n_cells(), CH(-1));
            VH nRootAlt;
            CH bRootAlt;
            if (!buildBlockSpanningTree(
                    nRootAlt, bRootAlt, cellTreePrecursorAlt, _cutSurfaceID[cutSurface.front().idx()]))
                throw std::logic_error("No leaf node in spanning tree");

            cycle = cycleThroughCutSurface(nStart, pStart, nRootAlt, bRootAlt, cellTreePrecursorAlt);
        }
    }
    return cycles;
}

list<HEH> ConstraintExtractor::pathThroughBlock(const VH& nStart, const set<VH>& nsTarget, const CH& b) const
{
    auto& mcMesh = mcMeshProps().mesh();

    vector<HEH> inHa(mcMesh.n_vertices(), HEH(-1));
    vector<bool> nVisited(mcMesh.n_vertices(), false);
    list<VH> nQ({nStart});
    nVisited[nStart.idx()] = true;
    while (!nQ.empty())
    {
        VH n = nQ.front();
        nQ.pop_front();

        if (nsTarget.count(n) != 0)
        {
            list<HEH> path;
            VH nCurr = n;
            while (nCurr != nStart)
            {
                HEH ha = inHa[nCurr.idx()];
                path.push_front(ha);
                nCurr = mcMesh.from_vertex_handle(ha);
            }
            return path;
        }

        for (HEH ha : mcMesh.outgoing_halfedges(n))
        {
            if (!contains(mcMesh.halfedge_cells(ha), b))
                continue;
            VH nTo = mcMesh.to_vertex_handle(ha);
            if (!nVisited[nTo.idx()])
            {
                inHa[nTo.idx()] = ha;
                nVisited[nTo.idx()] = true;
                nQ.push_back(nTo);
            }
        }
    }
    throw std::logic_error("No path found");
}

void ConstraintExtractor::determineEquivalentEndpoints(TetPathConstraint& path,
                                                       list<CH>& pathTets,
                                                       const map<EH, EH>& e2parent,
                                                       const map<FH, FH>& f2parent) const
{
    auto& tetMesh = meshProps().mesh();

    for (bool from : {true, false})
    {
        auto& v = from ? path.vFrom : path.vTo;
        assert(meshProps().get<MC_NODE>(v).is_valid());
        if (meshProps().get<IS_ORIGINAL_V>(v))
            continue;

        ConstraintNodeType type = _constraintNodeType[meshProps().get<MC_NODE>(v).idx()];

        assert(type != NATIVE);
        assert(type != NONE);
        // If on singular edge: shift to original edge endpoint
        if (type == ARTIFICIAL_ON_LINK)
        {
            // Find original edge
            EH eSingular
                = findMatching(tetMesh.vertex_edges(v), [&](const EH& e) { return meshProps().get<IS_SINGULAR>(e); });

            // pave connection to pathTets start/end tets
            auto connector = connectingTetPath(
                from ? pathTets.front() : pathTets.back(), tetMesh.halfedge_handle(eSingular, 0), v);
            if (from)
                pathTets.insert(pathTets.begin(), connector.rbegin(), connector.rend());
            else
                pathTets.insert(pathTets.end(), connector.begin(), connector.end());

            // Restore original
            auto it = e2parent.find(eSingular);
            while (it != e2parent.end())
            {
                eSingular = it->second;
                it = e2parent.find(eSingular);
            }
            assert(eSingular.is_valid());

            // Snap to any vertex of original
            assert(meshProps().get<IS_ORIGINAL_V>(tetMesh.edge_vertices(eSingular)[0]));
            v = tetMesh.edge_vertices(eSingular)[0];
        }
        else if (type == ARTIFICIAL_ON_BOUNDARY)
        {
            // Find original face
            HFH hfBoundary
                = findMatching(tetMesh.vertex_halffaces(v), [&](const HFH& hf) { return tetMesh.is_boundary(hf); });

            // pave connection to pathTets start/end tets
            auto connector = connectingTetPath(
                from ? pathTets.front() : pathTets.back(), tetMesh.opposite_halfface_handle(hfBoundary), v);
            if (from)
                pathTets.insert(pathTets.begin(), connector.rbegin(), connector.rend());
            else
                pathTets.insert(pathTets.end(), connector.begin(), connector.end());

            // Restore original
            FH fBoundary = tetMesh.face_handle(hfBoundary);
            auto it = f2parent.find(fBoundary);
            while (it != f2parent.end())
            {
                fBoundary = it->second;
                it = f2parent.find(fBoundary);
            }
            assert(fBoundary.is_valid());

            // Snap to any vertex of original
            assert(meshProps().get<IS_ORIGINAL_V>(*tetMesh.fv_iter(fBoundary)));
            v = *tetMesh.fv_iter(fBoundary);
        }
    }
}

void ConstraintExtractor::walkToMatchingTet(const EH& e, list<CH>::iterator& itTet, Transition& transOrigCurr) const
{
    auto& tetMesh = meshProps().mesh();
    bool hasHe = contains(tetMesh.cell_edges(*itTet), e);
    list<CH> connectingTets({*itTet});
    while (!hasHe)
    {
        itTet++;
        connectingTets.push_back({*itTet});
        CH tetPost = *itTet;
        hasHe = contains(tetMesh.cell_edges(tetPost), e);
    }
    transOrigCurr = transOrigCurr.chain(transitionAlongPath(connectingTets));
}

std::pair<int, int>
ConstraintExtractor::getBoundaryCoords(const CH& tetStart, const VH& vConn, const Transition& transOrigCurr) const
{
    auto& tetMesh = meshProps().mesh();

    // Get direction of boundary normal
    HFH boundaryHf;
    vector<bool> tetVisited(tetMesh.n_cells(), false);
    forVertexNeighbourTetsInBlock(vConn,
                                  tetStart,
                                  [&boundaryHf, vConn, this](const CH& tet)
                                  {
                                      boundaryHf = findMatching(
                                          meshProps().mesh().cell_halffaces(tet),
                                          [&](const HFH& hf)
                                          {
                                              return contains(meshProps().mesh().halfface_vertices(hf), vConn)
                                                     && meshProps().mesh().is_boundary(
                                                         meshProps().mesh().opposite_halfface_handle(hf));
                                          });
                                      return boundaryHf.is_valid();
                                  });
    if (!boundaryHf.is_valid())
        throw std::logic_error("Artificial singular node that is not incident on a cyclic singularity");
    auto connectingTets = connectingTetPath(tetStart, boundaryHf, vConn);
    connectingTets.push_front(tetStart);
    int coord
        = toCoord(transOrigCurr.chain(transitionAlongPath(connectingTets)).invert().rotate(normalDirUVW(boundaryHf)));
    return {coord == 0 ? 1 : 0, coord == 2 ? 1 : 2};
}

int ConstraintExtractor::getCycleCoord(const CH& tetStart, const VH& vConn, const Transition& transOrigCurr) const
{
    auto& tetMesh = meshProps().mesh();

    // Get direction of cycle
    HEH cycleHe;
    CH cycleTet;
    vector<bool> tetVisited(tetMesh.n_cells(), false);
    forVertexNeighbourTetsInBlock(vConn,
                                  tetStart,
                                  [&cycleHe, &cycleTet, vConn, &tetMesh, this](const CH& tet)
                                  {
                                      for (HEH he : tetMesh.cell_halfedges(tet))
                                      {
                                          EH a = meshProps().get<MC_ARC>(tetMesh.edge_handle(he));
                                          if (tetMesh.from_vertex_handle(he) == vConn && a.is_valid()
                                              && _arcIsCritical[a.idx()])
                                          {
                                              cycleHe = he;
                                              cycleTet = tet;
                                              return true;
                                          }
                                      }
                                      return false;
                                  });
    if (!cycleHe.is_valid())
        throw std::logic_error("Artificial singular node that is not incident on a cyclic singularity");
    auto connectingTets = connectingTetPath(tetStart, cycleHe, vConn);
    connectingTets.push_front(tetStart);
    auto& chartOrig = meshProps().ref<CHART_ORIG>(cycleTet);
    Vec3Q uvwOrig1 = chartOrig.at(meshProps().mesh().from_vertex_handle(cycleHe));
    Vec3Q uvwOrig2 = chartOrig.at(meshProps().mesh().to_vertex_handle(cycleHe));
    return toCoord(
        transOrigCurr.chain(transitionAlongPath(connectingTets)).invert().rotate(toDir(uvwOrig1 - uvwOrig2)));
}

}; // namespace mc3d
