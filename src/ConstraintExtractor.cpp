#include "QGP3D/ConstraintExtractor.hpp"

#include "QGP3D/MCQuantizer.hpp"

#include "OpenVolumeMesh/FileManager/FileManager.hh"

namespace qgp3d
{

ConstraintExtractor::ConstraintExtractor(const TetMeshProps& meshProps)
    : TetMeshNavigator(meshProps), MCMeshNavigator(meshProps)
{
}

vector<vector<OVM::HalfEdgeHandle>> ConstraintExtractor::getCriticalSkeletonArcs()
{
    auto& mcMesh = _mcMeshPropsC.mesh;

    assignNodeTypes();

    // Mark cut patches.
    markCutPatches();

    vector<vector<OVM::HalfEdgeHandle>> nodeTreePaths;
    // Mark one node per cut patch
    _isSkeletonArc = vector<bool>(mcMesh.n_edges(), false);

    // Initialize MSP with all arcs
    {
        vector<bool> nVisited(mcMesh.n_edges(), false);

        nVisited[mcMesh.v_iter()->idx()] = true;
        list<OVM::VertexHandle> nQ({*mcMesh.v_iter()});
        while (!nQ.empty())
        {
            auto n = nQ.front();
            nQ.pop_front();

            for (auto ha : mcMesh.outgoing_halfedges(n))
            {
                auto nTo = mcMesh.to_vertex_handle(ha);
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
    list<OVM::EdgeHandle> aQ;
    for (auto a : mcMesh.edges())
        for (auto n : mcMesh.edge_vertices(a))
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
        auto a = aQ.front();
        aQ.pop_front();
        for (auto n : mcMesh.edge_vertices(a))
            if (countIncidentSkeletonArcs(n) == 1 && _constraintNodeType[n.idx()] == NONE)
                for (auto aNext : mcMesh.vertex_edges(n))
                    if (_isSkeletonArc[aNext.idx()])
                    {
                        _isSkeletonArc[aNext.idx()] = false;
                        aQ.push_back(aNext);
                    }
    }

    // Extract paths from MSP
    {
        OVM::VertexHandle leaf;
        // Start from a leaf:
        for (auto n : mcMesh.vertices())
        {
            if (countIncidentSkeletonArcs(n) != 1)
                continue;
            if (_constraintNodeType[n.idx()] == NONE)
                throw std::logic_error("Invalid MSP leaf");
            leaf = n;
            break;
        }
        // TODO handle case where no leaf node exists
        if (leaf.is_valid())
        {
            // exhaustive graph search from leaf
            vector<bool> nVisited(mcMesh.n_vertices(), false);
            nVisited[leaf.idx()] = true;
            list<OVM::VertexHandle> nQ({leaf});
            list<vector<OVM::HalfEdgeHandle>> pathQ({vector<OVM::HalfEdgeHandle>()});
            while (!nQ.empty())
            {
                auto n = nQ.front();
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

                for (auto ha : mcMesh.outgoing_halfedges(n))
                {
                    if (!_isSkeletonArc[mcMesh.edge_handle(ha).idx()])
                        continue;
                    auto nTo = mcMesh.to_vertex_handle(ha);
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
    vector<vector<OVM::HalfEdgeHandle>> cycles = getCutSurfaceCycles();
    LOG(INFO) << cycles.size() << " cycles through " << _cutSurfaces.size() << " cut surfaces determined";
    for (auto& cycle : cycles)
        nodeTreePaths.emplace_back(cycle);

    return nodeTreePaths;
}

vector<ConstraintExtractor::TetPathConstraint>
ConstraintExtractor::getTetPathConstraints(const vector<vector<OVM::HalfEdgeHandle>>& haSequences)
{
    assert(_meshPropsC.isAllocated<CHART_ORIG>());
    assert(_meshPropsC.isAllocated<TRANSITION_ORIG>());
    assert(_meshPropsC.isAllocated<IS_ORIGINAL_VTX>());

    auto& mcMesh = _mcMeshPropsC.mesh;
    auto& tetMesh = _meshPropsC.mesh;
    vector<TetPathConstraint> tetPathConstraints;

    map<OVM::CellHandle, OVM::CellHandle> tet2parent;
    assert(_meshPropsC.isAllocated<CHILD_CELLS>());
    const auto& childTetProp = _meshPropsC.prop<CHILD_CELLS>();
    for (auto kv : childTetProp)
        for (auto child : kv.second)
            tet2parent[child] = kv.first;

    map<OVM::EdgeHandle, OVM::EdgeHandle> e2parent;
    assert(_meshPropsC.isAllocated<CHILD_EDGES>());
    const auto& childEdgeProp = _meshPropsC.prop<CHILD_EDGES>();
    for (auto kv : childEdgeProp)
        for (auto child : kv.second)
            e2parent[child] = kv.first;

    map<OVM::FaceHandle, OVM::FaceHandle> f2parent;
    assert(_meshPropsC.isAllocated<CHILD_FACES>());
    const auto& childFaceProp = _meshPropsC.prop<CHILD_FACES>();
    for (auto kv : childFaceProp)
        for (auto child : kv.second)
            f2parent[child] = kv.first;

    for (auto& haSequence : haSequences)
    {
        list<OVM::HalfEdgeHandle> pathHas(haSequence.begin(), haSequence.end());

        tetPathConstraints.emplace_back();
        auto& path = tetPathConstraints.back();

        auto pathHes = list<OVM::HalfEdgeHandle>();
        for (auto ha : pathHas)
        {
            auto haHes = _mcMeshPropsC.haHalfedges(ha);
            pathHes.insert(pathHes.end(), haHes.begin(), haHes.end());
        }

        auto nStart = mcMesh.from_vertex_handle(haSequence.front());
        auto nEnd = mcMesh.to_vertex_handle(haSequence.back());

        ConstraintNodeType startNodeType = _constraintNodeType[nStart.idx()];
        ConstraintNodeType endNodeType = _constraintNodeType[nEnd.idx()];

        assert(startNodeType != NONE && endNodeType != NONE);

        list<OVM::CellHandle> pathTets;
        for (auto he : pathHes)
        {
            // First tet may be any tet
            if (pathTets.empty())
            {
                pathTets.emplace_back(*tetMesh.hec_iter(he));
                continue;
            }

            auto vConn = tetMesh.from_vertex_handle(he);

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
        determineEquivalentEndpoints(path, pathTets, e2parent, f2parent, tet2parent);

        // Now get the original tet sequence
        for (auto tet : pathTets)
        {
            auto tetOrig = tet;
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

        vector<OVM::CellHandle> pathCheck;
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

            for (auto ha : pathHas)
            {
                {
                    // We are at the first halfedge of the current halfarc
                    // We can use this to get the arc direction under the orig param
                    auto tet = *itTet;
                    auto he = *itHe;
                    auto& chartOrig = _meshPropsC.ref<CHART_ORIG>(tet);
                    auto uvwOrig1 = chartOrig.at(_meshPropsC.mesh.from_vertex_handle(he));
                    auto uvwOrig2 = chartOrig.at(_meshPropsC.mesh.to_vertex_handle(he));
                    auto dir = transOrigCurr.invert().rotate(toDir(uvwOrig2 - uvwOrig1));

                    if (dim(dir) != 1)
                        throw std::logic_error("invalid arc dir");

                    // DLOG(INFO) << "Path " << path.vFrom << "->" << path.vTo << ": accumulating offset "
                    //            << Vec3d(toVec(dir)) * _mcMeshPropsC.get<ARC_INT_LENGTH>(mcMesh.edge_handle(ha))
                    //            << " of ha " << ha;
                    path.offset += toVec(dir) * _mcMeshPropsC.get<ARC_INT_LENGTH>(mcMesh.edge_handle(ha));
                }

                // Go along pathHes until vTo is reached;
                auto nTo = mcMesh.to_vertex_handle(ha);
                auto vTo = _mcMeshPropsC.get<NODE_MESH_VERTEX>(nTo);
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

                // if (he == pathHes.back())
                // {
                //     DLOG(INFO) << "FINAL ROTATION " << transOrigCurr.invert().rotation;
                // }
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

list<OVM::CellHandle> ConstraintExtractor::connectingTetPath(const OVM::CellHandle& tetStart,
                                                             const OVM::HalfEdgeHandle& heTarget,
                                                             const OVM::VertexHandle& vConn) const
{
    assert(heTarget.is_valid());
    auto& tetMesh = _meshPropsC.mesh;
    set<OVM::CellHandle> tetsVisited({tetStart});
    list<list<OVM::CellHandle>> tetQ({{tetStart}});
    bool hasHe = false;
    while (!tetQ.empty())
    {
        auto tets = tetQ.front();
        tetQ.pop_front();

        for (auto e : tetMesh.cell_edges(tets.back()))
            if (e == tetMesh.edge_handle(heTarget))
            {
                hasHe = true;
                break;
            }
        if (hasHe)
        {
            tets.pop_front(); // First cell is tetStart
            return tets;
        }

        for (auto tetNext : tetMesh.cell_cells(tets.back()))
        {
            if (tetsVisited.find(tetNext) != tetsVisited.end())
                continue;
            bool hasVConn = false;
            for (auto v : tetMesh.cell_vertices(tetNext))
                if (v == vConn)
                {
                    hasVConn = true;
                    break;
                }
            if (!hasVConn)
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

list<OVM::CellHandle> ConstraintExtractor::connectingTetPath(const OVM::CellHandle& tetStart,
                                                             const OVM::HalfFaceHandle& hfTarget,
                                                             const OVM::VertexHandle& vConn) const
{
    assert(hfTarget.is_valid());
    auto& tetMesh = _meshPropsC.mesh;
    set<OVM::CellHandle> tetsVisited({tetStart});
    list<list<OVM::CellHandle>> tetQ({{tetStart}});
    bool hasHf = false;
    while (!tetQ.empty())
    {
        auto tets = tetQ.front();
        tetQ.pop_front();

        for (auto hf : tetMesh.cell_halffaces(tets.back()))
            if (hf == hfTarget)
            {
                hasHf = true;
                break;
            }
        if (hasHf)
        {
            tets.pop_front(); // First cell is tetStart
            return tets;
        }

        for (auto tetNext : tetMesh.cell_cells(tets.back()))
        {
            if (tetsVisited.find(tetNext) != tetsVisited.end())
                continue;
            bool hasVConn = false;
            for (auto v : tetMesh.cell_vertices(tetNext))
                if (v == vConn)
                {
                    hasVConn = true;
                    break;
                }
            if (!hasVConn)
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

Transition ConstraintExtractor::transitionAlongPath(const list<OVM::CellHandle>& path) const
{
    if (path.size() < 2)
        return Transition();

    auto& tetMesh = _meshPropsC.mesh;
    Transition transTotal;
    for (auto itTet = path.begin(); itTet != --path.end();)
    {
        auto tetPre = *itTet;
        itTet++;
        auto tetPost = *itTet;
        // Accumulate transition
        OVM::HalfFaceHandle hfBetween;
        for (auto hf : tetMesh.cell_halffaces(tetPre))
            if (tetMesh.incident_cell(tetMesh.opposite_halfface_handle(hf)) == tetPost)
            {
                hfBetween = hf;
                break;
            }
        assert(hfBetween.is_valid());
        auto trans = _meshPropsC.hfTransitionOrig(hfBetween);
        transTotal = transTotal.chain(trans);
    }
    return transTotal;
}

int ConstraintExtractor::countCutPatches(const OVM::EdgeHandle& a) const
{
    auto& mcMesh = _mcMeshPropsC.mesh;
    int n = 0;
    for (auto p : mcMesh.edge_faces(a))
        if (_isCutPatch[p.idx()])
            n++;
    return n;
}

bool ConstraintExtractor::isOnActualCut(const OVM::VertexHandle& n) const
{
    auto& mcMesh = _mcMeshPropsC.mesh;
    for (auto p : mcMesh.vertex_faces(n))
        if (_isCutPatch[p.idx()] && !mcMesh.is_boundary(p))
            return true;
    return false;
}

int ConstraintExtractor::countIncidentSkeletonArcs(const OVM::VertexHandle& n) const
{
    auto& mcMesh = _mcMeshPropsC.mesh;
    int count = 0;
    for (auto a : mcMesh.vertex_edges(n))
        if (_isSkeletonArc[a.idx()])
            count++;
    return count;
}

void ConstraintExtractor::getCutSurfaces(vector<vector<OVM::FaceHandle>>& cutSurfaces, vector<int>& p2cutSurface)
{
    assignNodeTypes();

    markCutPatches();
    cutSurfaces = _cutSurfaces;
    p2cutSurface = _cutSurfaceID;
}

void ConstraintExtractor::assignNodeTypes()
{
    auto& mcMesh = _mcMeshPropsC.mesh;

    vector<MCQuantizer::CriticalLink> criticalLinks;
    map<OVM::EdgeHandle, int> a2criticalLinkIdx;
    map<OVM::VertexHandle, vector<int>> n2criticalLinksOut;
    map<OVM::VertexHandle, vector<int>> n2criticalLinksIn;
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
            auto type = nodeType(kv.first);
            if (type.first == SingularNodeType::SINGULAR || type.second == FeatureNodeType::FEATURE
                || (type.first == SingularNodeType::SEMI_SINGULAR && type.second == FeatureNodeType::SEMI_FEATURE))
                _constraintNodeType[kv.first.idx()] = NATIVE;
            else
                _constraintNodeType[kv.first.idx()] = ARTIFICIAL_ON_LINK;
        }

    vector<BoundaryRegion> boundaryRegions;
    map<OVM::HalfFaceHandle, int> hp2boundaryID;
    getBoundaryRegions(boundaryRegions, hp2boundaryID);

    int nBorderless = 0;
    for (auto& region : boundaryRegions)
        if (region.boundaryHas.empty())
            nBorderless++;
    LOG(INFO) << nBorderless << " borderless boundary regions out of " << boundaryRegions.size();

    // Check for boundary regions with no critical nodes or critical arcs in them
    for (auto& region : boundaryRegions)
    {
        bool hasCriticalNode = false;
        for (auto n : region.ns)
            if (n2criticalLinksIn.count(n) != 0 || n2criticalLinksOut.count(n) != 0)
            {
                hasCriticalNode = true;
                break;
            }
        if (hasCriticalNode)
            continue;
        assert(!region.ns.empty());
        _constraintNodeType[region.ns.begin()->idx()] = ARTIFICIAL_ON_BOUNDARY;
    }

    int nNative = 0;
    int nArtificialOnLink = 0;
    int nArtificialOnBoundary = 0;
    for (auto n : mcMesh.vertices())
        if (_constraintNodeType[n.idx()] == NATIVE)
            nNative++;
        else if (_constraintNodeType[n.idx()] == ARTIFICIAL_ON_BOUNDARY)
            nArtificialOnBoundary++;
        else if (_constraintNodeType[n.idx()] == ARTIFICIAL_ON_LINK)
            nArtificialOnLink++;

    int nBoth = 0;
    for (auto n : mcMesh.vertices())
        if (nodeType(n).first == SingularNodeType::SINGULAR && nodeType(n).second == FeatureNodeType::FEATURE)
            nBoth++;

    LOG(INFO) << "Nodes that are feature and singular: " << nBoth;
    LOG(INFO) << nNative << " native singular nodes, " << nArtificialOnLink
              << " singular nodes inserted on cyclic singular links, " << nArtificialOnBoundary
              << " singular nodes inserted on borderless boundary region ("
              << nNative + nArtificialOnBoundary + nArtificialOnLink << " total)";
}

void ConstraintExtractor::markCutPatches()
{
    auto& mcMesh = _mcMeshPropsC.mesh;
    _cellTreePrecursor = vector<OVM::CellHandle>(mcMesh.n_cells(), OVM::CellHandle(-1));

    // From AlgoHex' generate_cut_surface()

    // 1. start with trivial cut (all faces are cut)
    _isCutPatch = vector<bool>(mcMesh.n_faces(), true);

    {
        // 2. reduce cut-faces by dual spanning tree
        vector<bool> bVisited(mcMesh.n_cells(), false);

        // Find strictly singular node
        for (auto n : mcMesh.vertices())
            if (_constraintNodeType[n.idx()] == NATIVE)
            {
                _bRoot = *mcMesh.vc_iter(n);
                _nRoot = n;
                break;
            }

        // Else, find artificial singular node
        if (!_nRoot.is_valid())
        {
            for (auto n : mcMesh.vertices())
                if (_constraintNodeType[n.idx()] == ARTIFICIAL_ON_LINK
                    || _constraintNodeType[n.idx()] == ARTIFICIAL_ON_BOUNDARY)
                {
                    _bRoot = *mcMesh.vc_iter(n);
                    _nRoot = n;
                    break;
                }
        }

        // Else fail
        assert(_nRoot.is_valid());
        assert(_bRoot.is_valid());
        // add seed
        list<OVM::CellHandle> bQ({_bRoot});
        bVisited[_bRoot.idx()] = true;

        // grow tree until all cells in connected component visited
        while (!bQ.empty())
        {
            // get opposite of next halfface
            auto b = bQ.front();
            bQ.pop_front();

            // if b == bRoot, handle first halffaces that contain nRoot, so nRoot does not land on cut surface
            list<OVM::HalfFaceHandle> hps;
            for (auto hp : mcMesh.cell_halffaces(b))
            {
                if (b == _bRoot)
                {
                    bool hasNRoot = false;
                    for (auto n : mcMesh.halfface_vertices(hp))
                        if (n == _nRoot)
                        {
                            hasNRoot = true;
                            break;
                        }
                    if (hasNRoot)
                        hps.push_front(hp);
                    else
                        hps.push_back(hp);
                }
                else
                    hps.push_back(hp);
            }

            for (auto hp : hps)
            {
                hp = mcMesh.opposite_halfface_handle(hp);
                // only process if not on boundary
                if (!mcMesh.is_boundary(hp))
                {
                    // check if cell already visited
                    OVM::CellHandle bNext = mcMesh.incident_cell(hp);
                    if (!bVisited[bNext.idx()])
                    {
                        // remove face from cutgraph
                        _isCutPatch[mcMesh.face_handle(hp).idx()] = false;
                        // process cell
                        bVisited[bNext.idx()] = true;
                        bQ.push_back(bNext);
                        _cellTreePrecursor[bNext.idx()] = b;
                    }
                }
            }
        }

#ifndef NDEBUG
        // assert all cells have been visited (so all have a precursor except _bRoot)
        for (auto b : mcMesh.cells())
            assert(bVisited[b.idx()]);
#endif
    }
    {
        // 3. further shrink cut surface
        list<OVM::FaceHandle> pQ;
        // initialize queue with edges (1) not on boundary, and (2) adjacent to one cut face
        for (auto p : mcMesh.faces())
            for (auto a : mcMesh.face_edges(p))
                if (countCutPatches(a) == 1)
                {
                    pQ.push_back(p);
                    _isCutPatch[p.idx()] = false;
                    break;
                }

        // shrink while possible
        while (!pQ.empty())
        {
            // get next element
            auto p = pQ.front();
            pQ.pop_front();

            for (auto a : mcMesh.face_edges(p))
                if (countCutPatches(a) == 1)
                    for (auto pNext : mcMesh.edge_faces(a))
                        if (_isCutPatch[pNext.idx()])
                        {
                            _isCutPatch[pNext.idx()] = false;
                            pQ.push_back(pNext);
                        }
        }
    }

    for (auto a : mcMesh.edges())
        if (countCutPatches(a) == 1)
            throw std::logic_error("ERROR: there are shrinking candidates left but queue is empty!");

    // Remove boundary
    for (auto p : mcMesh.faces())
        if (mcMesh.is_boundary(p))
        {
            assert(_isCutPatch[p.idx()]);
            _isCutPatch[p.idx()] = false;
        }

    // 4. decompose cut surface into manifold pieces
    _cutSurfaces.clear();
    _cutSurfaceID = vector<int>(mcMesh.n_faces(), -1);
    vector<bool> pVisited(mcMesh.n_faces(), false);
    for (auto p : mcMesh.faces())
        if (_isCutPatch[p.idx()] && !pVisited[p.idx()])
        {
            // grow a new manifold patch
            _cutSurfaces.emplace_back();
            _cutSurfaces.back().push_back(p);
            _cutSurfaceID[p.idx()] = _cutSurfaces.size() - 1;
            pVisited[p.idx()] = true;

            list<OVM::FaceHandle> pQ({p});
            vector<OVM::FaceHandle> fhs;
            while (!pQ.empty())
            {
                auto pCurr = pQ.front();
                pQ.pop_front();

                for (auto a : mcMesh.face_edges(pCurr))
                {
                    if (!mcMesh.is_boundary(a) && countCutPatches(a) == 2)
                    {
                        for (auto pNext : mcMesh.edge_faces(a))
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

vector<vector<OVM::HalfEdgeHandle>> ConstraintExtractor::getCutSurfaceCycles()
{
    auto& mcMesh = _mcMeshPropsC.mesh;

    if (_cellTreePrecursor[_bRoot.idx()].is_valid())
        throw std::logic_error("Root has precursor");

    vector<vector<OVM::HalfEdgeHandle>> cycles;
    for (auto& cutSurface : _cutSurfaces)
    {
        cycles.emplace_back();
        auto& cycle = cycles.back();
        // Get any face and trace path back to root in both directions
        if (cutSurface.empty())
            throw std::logic_error("Empty cut surface");
        auto pStart = cutSurface.front();
        auto nStart = *mcMesh.fv_iter(pStart);

        bool nRootOnCutSurface = false;
        for (auto p : cutSurface)
        {
            for (auto n : mcMesh.face_vertices(p))
                if (n == _nRoot)
                {
                    nRootOnCutSurface = true;
                    break;
                }
            if (nRootOnCutSurface)
                break;
        }
        if (!nRootOnCutSurface)
        {
            for (bool forward : {false, true})
            {
                auto b = mcMesh.incident_cell(mcMesh.halfface_handle(pStart, forward));
                auto n = nStart;
                list<OVM::HalfEdgeHandle> halfpath;
                // Find a path connecting nStart to nRoot and going through all the cells in the sequence
                while (b != _bRoot)
                {
                    auto bNext = _cellTreePrecursor[b.idx()];

                    OVM::FaceHandle pInterface;
                    for (auto hp : mcMesh.cell_halffaces(b))
                        if (!_isCutPatch[mcMesh.face_handle(hp).idx()]
                            && mcMesh.incident_cell(mcMesh.opposite_halfface_handle(hp)) == bNext)
                        {
                            pInterface = mcMesh.face_handle(hp);
                            break;
                        }
                    set<OVM::VertexHandle> nsTarget;
                    for (auto nTarget : mcMesh.face_vertices(pInterface))
                        nsTarget.insert(nTarget);

                    auto connector = pathThroughBlock(n, nsTarget, b);
                    halfpath.insert(halfpath.end(), connector.begin(), connector.end());

                    n = connector.empty() ? n : mcMesh.to_vertex_handle(connector.back());
                    if (nsTarget.count(n) == 0)
                        throw std::logic_error("Found vertex not in nsTarget");
                    b = bNext;
                }
                auto connector = pathThroughBlock(n, {_nRoot}, _bRoot);
                halfpath.insert(halfpath.end(), connector.begin(), connector.end());

                if (forward)
                {
                    cycle.insert(cycle.end(), halfpath.begin(), halfpath.end());
                }
                else
                {
                    halfpath.reverse();
                    for (auto ha : halfpath)
                        cycle.emplace_back(mcMesh.opposite_halfedge_handle(ha));
                }
            }
            if (mcMesh.from_vertex_handle(cycle.front()) != mcMesh.to_vertex_handle(cycle.back()))
                throw std::logic_error("Unclosed cycle");
        }
        else
        {
            // Find ANOTHER singular node, suitable as root node
            OVM::VertexHandle nRootAlt;
            OVM::CellHandle bRootAlt;
            vector<OVM::CellHandle> cellTreePrecursorAlt(mcMesh.n_cells(), OVM::CellHandle(-1));

            // Find strictly singular node
            for (auto n : mcMesh.vertices())
                if (_constraintNodeType[n.idx()] == NATIVE)
                {
                    bool onCutSurface = false;
                    for (auto p : mcMesh.vertex_faces(n))
                        if (_cutSurfaceID[p.idx()] == _cutSurfaceID[cutSurface.front().idx()])
                        {
                            onCutSurface = true;
                            break;
                        }
                    if (onCutSurface)
                        continue;
                    nRootAlt = n;
                    bRootAlt = *mcMesh.vc_iter(n);
                    break;
                }

            // Else, find artificial singular node
            if (!nRootAlt.is_valid())
            {
                for (auto n : mcMesh.vertices())
                    if (_constraintNodeType[n.idx()] == ARTIFICIAL_ON_LINK
                        || _constraintNodeType[n.idx()] == ARTIFICIAL_ON_BOUNDARY)
                    {
                        bool onCutSurface = false;
                        for (auto p : mcMesh.vertex_faces(n))
                            if (_cutSurfaceID[p.idx()] == _cutSurfaceID[cutSurface.front().idx()])
                            {
                                onCutSurface = true;
                                break;
                            }
                        if (onCutSurface)
                            continue;
                        nRootAlt = n;
                        bRootAlt = *mcMesh.vc_iter(n);
                        break;
                    }
            }

            // Else fail
            assert(nRootAlt.is_valid());
            assert(bRootAlt.is_valid());

            // Build new cell precursor tree, but do not cross any established cut surface

            // add seed
            vector<bool> bVisited(mcMesh.n_cells(), false);
            list<OVM::CellHandle> bQ({bRootAlt});
            bVisited[bRootAlt.idx()] = true;

            // grow tree until all cells in connected component visited
            while (!bQ.empty())
            {
                // get opposite of next halfface
                auto b = bQ.front();
                bQ.pop_front();

                // if b == bRoot, handle first halffaces that contain nRoot, so nRoot does not land on cut surface
                list<OVM::HalfFaceHandle> hps;
                for (auto hp : mcMesh.cell_halffaces(b))
                {
                    if (b == bRootAlt)
                    {
                        bool hasNRoot = false;
                        for (auto n : mcMesh.halfface_vertices(hp))
                            if (n == nRootAlt)
                            {
                                hasNRoot = true;
                                break;
                            }
                        if (hasNRoot)
                            hps.push_front(hp);
                        else
                            hps.push_back(hp);
                    }
                    else
                        hps.push_back(hp);
                }

                for (auto hp : hps)
                {
                    hp = mcMesh.opposite_halfface_handle(hp);
                    auto p = mcMesh.face_handle(hp);
                    // only process if not on boundary and not on cut surface
                    if (!mcMesh.is_boundary(hp) && !_isCutPatch[p.idx()])
                    {
                        // check if cell already visited
                        OVM::CellHandle bNext = mcMesh.incident_cell(hp);
                        if (!bVisited[bNext.idx()])
                        {
                            // process cell
                            bVisited[bNext.idx()] = true;
                            bQ.push_back(bNext);
                            cellTreePrecursorAlt[bNext.idx()] = b;
                        }
                    }
                }
            }
#ifndef NDEBUG
            // assert all cells have been visited (so all have a precursor except bRootAlt)
            for (auto b : mcMesh.cells())
                assert(bVisited[b.idx()]);
#endif

            for (bool forward : {false, true})
            {
                auto b = mcMesh.incident_cell(mcMesh.halfface_handle(pStart, forward));
                auto n = nStart;
                list<OVM::HalfEdgeHandle> halfpath;
                // Find a path connecting nStart to nRoot and going through all the cells in the sequence
                while (b != bRootAlt)
                {
                    auto bNext = cellTreePrecursorAlt[b.idx()];

                    OVM::FaceHandle pInterface;
                    for (auto hp : mcMesh.cell_halffaces(b))
                        if (!_isCutPatch[mcMesh.face_handle(hp).idx()]
                            && mcMesh.incident_cell(mcMesh.opposite_halfface_handle(hp)) == bNext)
                        {
                            pInterface = mcMesh.face_handle(hp);
                            break;
                        }
                    set<OVM::VertexHandle> nsTarget;
                    for (auto nTarget : mcMesh.face_vertices(pInterface))
                        nsTarget.insert(nTarget);

                    auto connector = pathThroughBlock(n, nsTarget, b);
                    halfpath.insert(halfpath.end(), connector.begin(), connector.end());

                    n = connector.empty() ? n : mcMesh.to_vertex_handle(connector.back());
                    if (nsTarget.count(n) == 0)
                        throw std::logic_error("Found vertex not in nsTarget");
                    b = bNext;
                }
                auto connector = pathThroughBlock(n, {nRootAlt}, bRootAlt);
                halfpath.insert(halfpath.end(), connector.begin(), connector.end());

                if (forward)
                {
                    cycle.insert(cycle.end(), halfpath.begin(), halfpath.end());
                }
                else
                {
                    halfpath.reverse();
                    for (auto ha : halfpath)
                        cycle.emplace_back(mcMesh.opposite_halfedge_handle(ha));
                }
            }
            if (mcMesh.from_vertex_handle(cycle.front()) != mcMesh.to_vertex_handle(cycle.back()))
                throw std::logic_error("Unclosed cycle");
        }
    }
    return cycles;
}

list<OVM::HalfEdgeHandle> ConstraintExtractor::pathThroughBlock(const OVM::VertexHandle& nStart,
                                                                const set<OVM::VertexHandle>& nsTarget,
                                                                const OVM::CellHandle& b)
{
    auto& mcMesh = _mcMeshPropsC.mesh;

    bool incidentOnB = false;
    for (auto b2 : mcMesh.vertex_cells(nStart))
        if (b2 == b)
            incidentOnB = true;
    if (!incidentOnB)
        throw std::logic_error("NStart Not incident on b");
    for (auto n : nsTarget)
    {
        incidentOnB = false;
        for (auto b2 : mcMesh.vertex_cells(n))
            if (b2 == b)
                incidentOnB = true;
        if (!incidentOnB)
            throw std::logic_error("Nstarget Not incident on b");
    }

    vector<OVM::HalfEdgeHandle> inHa(mcMesh.n_vertices(), OVM::HalfEdgeHandle(-1));
    vector<bool> nVisited(mcMesh.n_vertices(), false);
    list<OVM::VertexHandle> nQ({nStart});
    nVisited[nStart.idx()] = true;
    while (!nQ.empty())
    {
        auto n = nQ.front();
        nQ.pop_front();

        if (nsTarget.count(n) != 0)
        {
            list<OVM::HalfEdgeHandle> path;
            auto nCurr = n;
            while (nCurr != nStart)
            {
                auto ha = inHa[nCurr.idx()];
                path.push_front(ha);
                nCurr = mcMesh.from_vertex_handle(ha);
            }
            return path;
        }

        for (auto ha : mcMesh.outgoing_halfedges(n))
        {
            bool allowed = false;
            for (auto bOther : mcMesh.halfedge_cells(ha))
                if (bOther == b)
                {
                    allowed = true;
                    break;
                }
            if (!allowed)
                continue;
            auto nTo = mcMesh.to_vertex_handle(ha);
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
                                                       list<OVM::CellHandle>& pathTets,
                                                       const map<OVM::EdgeHandle, OVM::EdgeHandle>& e2parent,
                                                       const map<OVM::FaceHandle, OVM::FaceHandle>& f2parent,
                                                       const map<OVM::CellHandle, OVM::CellHandle>& tet2parent) const
{
    auto& tetMesh = _meshPropsC.mesh;

    (void)tet2parent;

    for (bool from : {true, false})
    {
        auto& v = from ? path.vFrom : path.vTo;
        assert(_meshPropsC.get<MC_NODE>(v).is_valid());
        if (_meshPropsC.get<IS_ORIGINAL_VTX>(v))
            continue;

        ConstraintNodeType type = _constraintNodeType[_meshPropsC.get<MC_NODE>(v).idx()];

        assert(type != NATIVE);
        assert(type != NONE);
        // If on singular edge: shift to original edge endpoint
        if (type == ARTIFICIAL_ON_LINK)
        {
            // Find original edge
            OVM::EdgeHandle eSingular;
            for (auto e : tetMesh.vertex_edges(v))
            {
                if (_meshPropsC.get<MC_ARC>(e).is_valid() && _mcMeshPropsC.get<IS_SINGULAR>(_meshPropsC.get<MC_ARC>(e)))
                {
                    eSingular = e;
                    break;
                }
            }

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
            assert(_meshPropsC.get<IS_ORIGINAL_VTX>(tetMesh.edge_vertices(eSingular)[0]));
            v = tetMesh.edge_vertices(eSingular)[0];
        }
        else if (type == ARTIFICIAL_ON_BOUNDARY)
        {
            // Find original face
            OVM::HalfFaceHandle hfBoundary;
            for (auto hf : tetMesh.vertex_halffaces(v))
            {
                if (tetMesh.is_boundary(hf))
                {
                    hfBoundary = hf;
                    break;
                }
            }

            // pave connection to pathTets start/end tets
            auto connector = connectingTetPath(
                from ? pathTets.front() : pathTets.back(), tetMesh.opposite_halfface_handle(hfBoundary), v);
            if (from)
                pathTets.insert(pathTets.begin(), connector.rbegin(), connector.rend());
            else
                pathTets.insert(pathTets.end(), connector.begin(), connector.end());

            // Restore original
            auto fBoundary = tetMesh.face_handle(hfBoundary);
            auto it = f2parent.find(fBoundary);
            while (it != f2parent.end())
            {
                fBoundary = it->second;
                it = f2parent.find(fBoundary);
            }
            assert(fBoundary.is_valid());

            // Snap to any vertex of original
            assert(_meshPropsC.get<IS_ORIGINAL_VTX>(*tetMesh.fv_iter(fBoundary)));
            v = *tetMesh.fv_iter(fBoundary);
        }
    }
}

void ConstraintExtractor::walkToMatchingTet(const OVM::EdgeHandle& e,
                                            list<OVM::CellHandle>::iterator& itTet,
                                            Transition& transOrigCurr) const
{
    auto& tetMesh = _meshPropsC.mesh;
    bool hasHe = false;
    for (auto e2 : tetMesh.cell_edges(*itTet))
        if (e2 == e)
        {
            hasHe = true;
            break;
        }
    list<OVM::CellHandle> connectingTets({*itTet});
    while (!hasHe)
    {
        itTet++;
        connectingTets.push_back({*itTet});
        auto tetPost = *itTet;
        for (auto e2 : tetMesh.cell_edges(tetPost))
            if (e2 == e)
            {
                hasHe = true;
                break;
            }
    }
    transOrigCurr = transOrigCurr.chain(transitionAlongPath(connectingTets));
}

std::pair<int, int> ConstraintExtractor::getBoundaryCoords(const OVM::CellHandle& tetStart,
                                                           const OVM::VertexHandle& vConn,
                                                           const Transition& transOrigCurr) const
{
    auto& tetMesh = _meshPropsC.mesh;

    // Get direction of boundary normal
    OVM::HalfFaceHandle boundaryHf;
    vector<bool> tetVisited(tetMesh.n_cells(), false);
    forVertexNeighbourTetsInBlock(vConn,
                                  tetStart,
                                  [&boundaryHf, vConn, &tetMesh](const OVM::CellHandle& tet)
                                  {
                                      for (auto hf : tetMesh.cell_halffaces(tet))
                                      {
                                          bool hasVConn = false;
                                          for (auto v : tetMesh.halfface_vertices(hf))
                                              if (v == vConn)
                                              {
                                                  hasVConn = true;
                                                  break;
                                              }
                                          if (hasVConn && tetMesh.is_boundary(tetMesh.opposite_halfface_handle(hf)))
                                          {
                                              boundaryHf = hf;
                                              return true;
                                          }
                                      }
                                      return false;
                                  });
    if (!boundaryHf.is_valid())
        throw std::logic_error("Artificial singular node that is not incident on a cyclic singularity");
    auto connectingTets = connectingTetPath(tetStart, boundaryHf, vConn);
    connectingTets.push_front(tetStart);
    auto coord = toCoord(transOrigCurr.chain(transitionAlongPath(connectingTets))
                             .invert()
                             .rotate(axisAlignedHalfFaceNormal(boundaryHf)));
    return {coord == 0 ? 1 : 0, coord == 2 ? 1 : 2};
}

int ConstraintExtractor::getCycleCoord(const OVM::CellHandle& tetStart,
                                       const OVM::VertexHandle& vConn,
                                       const Transition& transOrigCurr) const
{
    auto& tetMesh = _meshPropsC.mesh;

    // Get direction of cycle
    OVM::HalfEdgeHandle cycleHe;
    OVM::CellHandle cycleTet;
    vector<bool> tetVisited(tetMesh.n_cells(), false);
    forVertexNeighbourTetsInBlock(vConn,
                                  tetStart,
                                  [&cycleHe, &cycleTet, vConn, &tetMesh, this](const OVM::CellHandle& tet)
                                  {
                                      for (auto he : tetMesh.cell_halfedges(tet))
                                      {
                                          auto a = _meshPropsC.get<MC_ARC>(tetMesh.edge_handle(he));
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
    auto& chartOrig = _meshPropsC.ref<CHART_ORIG>(cycleTet);
    auto uvwOrig1 = chartOrig.at(_meshPropsC.mesh.from_vertex_handle(cycleHe));
    auto uvwOrig2 = chartOrig.at(_meshPropsC.mesh.to_vertex_handle(cycleHe));
    return toCoord(
        transOrigCurr.chain(transitionAlongPath(connectingTets)).invert().rotate(toDir(uvwOrig1 - uvwOrig2)));
}

}; // namespace qgp3d
