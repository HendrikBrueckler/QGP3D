#include "QGP3D/SeparationChecker.hpp"

namespace qgp3d
{

bool SeparationChecker::GreaterPathLengthCompare::operator()(const WeaklyMonotonousPath& p1,
                                                             const WeaklyMonotonousPath& p2) const
{
    return p1.length > p2.length || (p1.length == p2.length && p1.path.size() > p2.path.size());
}

SeparationChecker::SeparationChecker(TetMeshProps& meshProps) : TetMeshNavigator(meshProps), MCMeshNavigator(meshProps)
{
}

int SeparationChecker::numHexesInQuantization() const
{
    const MCMesh& mc = mcMeshProps().mesh();
    int nHexes = 0;
    for (CH b : mc.cells())
    {
        int nBlockHexes = 1;
        for (UVWDir dir : {UVWDir::NEG_U_NEG_V, UVWDir::NEG_U_NEG_W, UVWDir::NEG_V_NEG_W})
        {
            int arcLen = 0;
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

void
SeparationChecker::findSeparationViolatingPaths(vector<CriticalLink>& criticalLinks,
                                                const vector<bool>& arcIsCritical,
                                                const vector<bool>& nodeIsCritical,
                                                const vector<bool>& patchIsCritical,
                                                vector<vector<pair<int, EH>>>& nonZeroSumArcs)
{
    nonZeroSumArcs.clear();
    vector<vector<pair<int, EH>>> failsafeNonZeroSumArcs;
    // For each critical link s1 find paths connecting s1 to other critical links s2 or surface patches p2
    // Then check for overlaps between these, accumulating the quantized edge lengths along the path as deltas.
    for (auto& criticalLink : criticalLinks)
    {
        // This has been shifted from calling functions to here. Makes no sense to have caller do this
        criticalLink.length = 0;
        for (HEH ha : criticalLink.pathHas)
            criticalLink.length += mcMeshProps().get<ARC_INT_LENGTH>(mcMeshProps().mesh().edge_handle(ha));

        traceExhaustPaths(criticalLink, arcIsCritical, nodeIsCritical, patchIsCritical, nonZeroSumArcs, failsafeNonZeroSumArcs);
    }

    _failsafeSeparatingPaths.insert(_failsafeSeparatingPaths.end(), failsafeNonZeroSumArcs.begin(), failsafeNonZeroSumArcs.end());
    _allSeparatingPaths.insert(_allSeparatingPaths.end(), nonZeroSumArcs.begin(), nonZeroSumArcs.end());
}

void SeparationChecker::traceExhaustPaths(const CriticalLink& criticalLink1,
                                                                const vector<bool>& arcIsCritical,
                                                                const vector<bool>& nodeIsCritical,
                                                                const vector<bool>& patchIsCritical,
                                                                vector<vector<pair<int, EH>>>& nonZeroSumArcs,
                                                                vector<vector<pair<int, EH>>>& failsafeNonZeroSumArcs) const
{
    const MCMesh& mcMesh = mcMeshProps().mesh();

    HEH criticalStartHa = criticalLink1.pathHas.empty() ? HEH() : criticalLink1.pathHas.front();

    vector<vector<HEH>> dir2has;
    map<HEH, int> haOrth2dir;
    if (criticalStartHa.is_valid())
    {
        // Categorize all vertical has by direction
        for (HFH criticalStartHp : mcMesh.halfedge_halffaces(criticalStartHa))
        {
            dir2has.emplace_back();

            auto& has = dir2has.back();

            CH bRef = mcMesh.incident_cell(criticalStartHp);
            if (!bRef.is_valid())
                bRef = mcMesh.incident_cell(mcMesh.opposite_halfface_handle(criticalStartHp));

            UVWDir dirStart = halfarcDirInBlock(criticalStartHa, bRef);

            HEH haOrth = mcMesh.prev_halfedge_in_halfface(criticalStartHa, criticalStartHp);
            while (halfarcDirInBlock(haOrth, bRef) == dirStart)
                haOrth = mcMesh.prev_halfedge_in_halfface(haOrth, criticalStartHp);
            haOrth = mcMesh.opposite_halfedge_handle(haOrth);
            has.emplace_back(haOrth);

            HEH haCurr = criticalStartHa;
            HFH hpCurr = criticalStartHp;

            do
            {
                bRef = mcMesh.incident_cell(hpCurr);
                if (!bRef.is_valid())
                    bRef = mcMesh.incident_cell(mcMesh.opposite_halfface_handle(hpCurr));

                dirStart = halfarcDirInBlock(haCurr, bRef);

                while (halfarcDirInBlock(mcMesh.next_halfedge_in_halfface(haCurr, hpCurr), bRef) == dirStart)
                    haCurr = mcMesh.next_halfedge_in_halfface(haCurr, hpCurr);
                haOrth = mcMesh.next_halfedge_in_halfface(haCurr, hpCurr);
                if (haOrth == has.front())
                    break;
                has.emplace_back(haOrth);

                HEH haOrthOpp = mcMesh.opposite_halfedge_handle(haOrth);
                hpCurr = findMatching(mcMesh.halfedge_halffaces(haOrthOpp),
                                      [&](const HFH& hpNext)
                                      {
                                          return mcMesh.face_handle(hpNext) != mcMesh.face_handle(hpCurr)
                                                 && std::find(criticalLink1.pathHas.begin(),
                                                              criticalLink1.pathHas.end(),
                                                              mcMesh.next_halfedge_in_halfface(haOrthOpp, hpNext))
                                                        != criticalLink1.pathHas.end();
                                      });
                if (hpCurr.is_valid())
                    haCurr = mcMesh.next_halfedge_in_halfface(haOrthOpp, hpCurr);
            } while (hpCurr.is_valid() && haCurr != criticalStartHa);
        }
        for (int i = 0; i < (int)dir2has.size(); i++)
            for (HEH ha : dir2has[i])
                haOrth2dir[ha] = i;
    }
    else
    {
        for (HEH ha : mcMesh.outgoing_halfedges(criticalLink1.nFrom))
            dir2has.push_back({ha});
        for (int i = 0; i < (int)dir2has.size(); i++)
            for (HEH ha : dir2has[i])
                haOrth2dir[ha] = i;
    }

    VH criticalStartN = criticalLink1.nFrom;
    CH bRefStart
        = criticalStartHa.is_valid() ? *mcMesh.hec_iter(criticalStartHa) : *mcMesh.vc_iter(criticalLink1.nFrom);
    UVWDir dirStartHa = criticalStartHa.is_valid() ? halfarcDirInBlock(criticalStartHa, bRefStart) : UVWDir::NONE;

    WeaklyMonotonousPath pStart;
    pStart.branchedOff = false;
    pStart.length = 0;
    pStart.n = criticalStartN;
    pStart.dirs1 = dirStartHa | -dirStartHa;
    pStart.path = {};
    pStart.monotonousDirs = ~pStart.dirs1; // Do not search along the link but orthogonally
    pStart.walkedDirs = UVWDir::NONE;
    pStart.deltaMin = (isNeg(dirStartHa) ? criticalLink1.length * toVec(dirStartHa) : Vec3i(0, 0, 0));
    pStart.deltaMax = (isNeg(dirStartHa) ? Vec3i(0, 0, 0) : criticalLink1.length * toVec(dirStartHa));
    if (criticalLink1.cyclic && dirStartHa != UVWDir::NONE)
    {
        pStart.deltaMin[toCoord(dirStartHa)] = INT_MIN;
        pStart.deltaMax[toCoord(dirStartHa)] = INT_MAX;
    }
    assert(pStart.deltaMin[0] <= pStart.deltaMax[0] && pStart.deltaMin[1] <= pStart.deltaMax[1]
           && pStart.deltaMin[2] <= pStart.deltaMax[2]);
    pStart.delta = Vec3i(0, 0, 0);
    pStart.bRefCurrent = bRefStart;
    pStart.transCurrent = Transition();

    using PathQueue
        = std::priority_queue<WeaklyMonotonousPath, std::deque<WeaklyMonotonousPath>, GreaterPathLengthCompare>;

    for (int dirIdx = 0; dirIdx < (int)dir2has.size(); dirIdx++)
    {
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
            determineNextHalfedges(pathCurrent, b2trans, ha2bRef);

            // Check for separation violations by current path
            if ((pathCurrent.walkedDirs & pathCurrent.monotonousDirs) != UVWDir::NONE)
            {
                UVWDir violationDir = UVWDir::NONE;
                // Check whether startinterval-criticalnode pairs overlap
                if (nodeIsCritical[pathCurrent.n.idx()]
                    && bboxOverlap(pathCurrent.deltaMin, pathCurrent.deltaMax, pathCurrent.delta, pathCurrent.delta))
                    violationDir = pathCurrent.walkedDirs & ~pathCurrent.dirs1;

                if (violationDir == UVWDir::NONE)
                {
                    // Check whether startinterval-criticallink pairs overlap
                    for (const auto& kv : ha2bRef)
                    {
                        HEH ha2 = kv.first;
                        EH a2 = mcMesh.edge_handle(ha2);
                        if (arcIsCritical[a2.idx()])
                        {
                            CH bRef2 = kv.second.front();
                            Transition trans2 = b2trans.at(bRef2);

                            UVWDir dir2 = trans2.invert().rotate(halfarcDirInBlock(ha2, bRef2));

                            UVWDir possibleViolationDir = pathCurrent.walkedDirs & ~(dir2 | -dir2) & ~pathCurrent.dirs1;
                            if ((possibleViolationDir & pathCurrent.monotonousDirs) == UVWDir::NONE)
                                continue;

                            // Measure just overlap with the arc, not the link
                            int length = mcMeshProps().get<ARC_INT_LENGTH>(a2);

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
                                violationDir = possibleViolationDir;
                                break;
                            }
                        }
                    }
                }

                if (violationDir == UVWDir::NONE)
                {
                    // Check for overlaps between startinterval-surfacepatch pairs
                    for (HFH hp2 : mcMesh.vertex_halffaces(pathCurrent.n))
                    {
                        if (patchIsCritical[mcMesh.face_handle(hp2).idx()])
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
                                violationDir = possibleViolationDir;
                                break;
                            }
                        }
                    }
                }

                if (violationDir != UVWDir::NONE)
                {
                    vector<EH> posSignArcs;
                    vector<EH> negSignArcs;
                    vector<EH> monotonousDirArcs;
                    for (UVWDir dim1dirPos : {UVWDir::POS_U, UVWDir::POS_V, UVWDir::POS_W})
                    {
                        UVWDir dim1dirAny = dim1dirPos | -dim1dirPos;
                        UVWDir intersection = violationDir & dim1dirAny;
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
                    break;
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
                        auto itDir = haOrth2dir.find(ha);
                        if (itDir != haOrth2dir.end())
                        {
                            if (itDir->second != dirIdx)
                                continue;
                            nextP.branchedOff = true;

                            nextP.monotonousDirs = nextP.monotonousDirs & walkedDir;
                        }
                        else if ((walkedDir & nextP.dirs1) == UVWDir::NONE)
                        {
                            // Dont allow branch-off at another arc than those allowed
                            continue;
                        }
                    }
                    if (nextP.branchedOff)
                        nextP.path.emplace_back(ha);

                    nextP.delta += mcMeshProps().get<ARC_INT_LENGTH>(mcMesh.edge_handle(ha)) * toVec(walkedDir);

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
}

bool SeparationChecker::bboxOverlap(const Vec3i& min1, const Vec3i& max1, const Vec3i& min2, const Vec3i& max2)
{
    int coordTouchingOrOverlaps = 0;
    for (int i = 0; i < 3; i++)
    {
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

void SeparationChecker::determineNextHalfedges(const WeaklyMonotonousPath& pathCurrent,
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

bool SeparationChecker::checkArcOverlap(const WeaklyMonotonousPath& pathCurrent,
                                        const HEH& ha,
                                        const CH& bRef,
                                        const Transition& trans) const
{
    auto& mcMesh = mcMeshProps().mesh();

    UVWDir dir2 = trans.invert().rotate(halfarcDirInBlock(ha, bRef));
    int arcLen = mcMeshProps().get<ARC_INT_LENGTH>(mcMesh.edge_handle(ha));

    Vec3i deltaMin2 = isNeg(dir2) ? pathCurrent.delta + arcLen * toVec(dir2) : pathCurrent.delta;
    Vec3i deltaMax2 = isNeg(dir2) ? pathCurrent.delta : pathCurrent.delta + arcLen * toVec(dir2);

    // Check for overlap
    return bboxOverlap(pathCurrent.deltaMin, pathCurrent.deltaMax, deltaMin2, deltaMax2);
}

bool SeparationChecker::checkPatchOverlap(const WeaklyMonotonousPath& pathCurrent,
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
        int lengthHa = mcMeshProps().get<ARC_INT_LENGTH>(aCurr);
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

bool SeparationChecker::checkP0Containment(const WeaklyMonotonousPath& pathCurrent) const
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
                    double length = mcMeshProps().get<ARC_INT_LENGTH>(mcMesh.edge_handle(ha));
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

} // namespace qgp3d
