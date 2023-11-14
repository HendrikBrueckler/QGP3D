#include "QGP3D/MCQuantizer.hpp"

#include <gurobi_c++.h>
#include <queue>

namespace qgp3d
{

bool MCQuantizer::GreaterPathLengthCompare::operator()(const WeaklyMonotonousPath& p1,
                                                       const WeaklyMonotonousPath& p2) const
{
    return p1.length > p2.length || (p1.length == p2.length && p1.path.size() > p2.path.size());
}

MCQuantizer::MCQuantizer(TetMeshProps& meshProps)
    : TetMeshNavigator(meshProps), TetMeshManipulator(meshProps), MCMeshNavigator(meshProps),
      MCMeshManipulator(meshProps)
{
}

MCQuantizer::RetCode MCQuantizer::quantizeArcLengths(double scaling, bool allowZero, bool allowNegative)
{
    const MCMesh& mc = mcMeshProps().mesh();
    const TetMesh& mesh = meshProps().mesh();

    try
    {
        GRBEnv env = GRBEnv(true);
#ifdef NDEBUG
        env.set(GRB_IntParam_LogToConsole, false);
#endif
        env.set(GRB_DoubleParam_TimeLimit, 5 * 60);
        env.start();

        GRBModel model = GRBModel(env);

        bool wasAllocated = mcMeshProps().isAllocated<ARC_DBL_LENGTH>();
        if (!wasAllocated)
            mcMeshProps().allocate<ARC_DBL_LENGTH>(0.0);
        // One integer length variable per arc
        std::map<EH, GRBVar> arc2var;
        for (EH arc : mc.edges())
        {
            if (!wasAllocated)
            {
                // Determine current arc length
                double length = 0.0;
                for (HEH he : mcMeshProps().ref<ARC_MESH_HALFEDGES>(arc))
                    length += edgeLengthUVW<CHART>(mesh.edge_handle(he));
                mcMeshProps().set<ARC_DBL_LENGTH>(arc, length);
            }
            double length = mcMeshProps().get<ARC_DBL_LENGTH>(arc);

            // Configure gurobi var
            GRBVar x = model.addVar(allowZero ? (allowNegative ? -GRB_INFINITY : 0.0) : 0.9,
                                    GRB_INFINITY,
                                    0.,
                                    GRB_INTEGER,
                                    std::string("Arc ") + std::to_string(arc.idx()));
            x.set(GRB_DoubleAttr_Start, length);
            arc2var[arc] = {x};
        }

        GRBQuadExpr objective = 0.;
        // Objective Function minimizes deviation from scaled seamless param
        for (CH block : mc.cells())
        {
            // For each axis U/V/W
            for (UVWDir axis : {UVWDir::NEG_V_NEG_W, UVWDir::NEG_U_NEG_W, UVWDir::NEG_U_NEG_V})
            {
                GRBLinExpr varSum = 0.;
                double lenSum = 0;
                // WARNING this way of using range based for loop (chaining 2 function calls)
                // depends on the first calls not returning rvalues!
                for (EH arc : mcMeshProps().ref<BLOCK_EDGE_ARCS>(block).at(axis))
                {
                    varSum += arc2var.at(arc);
                    lenSum += mcMeshProps().get<ARC_DBL_LENGTH>(arc);
                }
                objective += (varSum - scaling * lenSum) * (varSum - scaling * lenSum);
            }
        }

        model.setObjective(objective, GRB_MINIMIZE);

        // Each patches opposite arc lengths must match
        for (FH patch : mc.faces())
        {
            HFH hp = mc.halfface_handle(patch, 0);
            if (mc.is_boundary(hp))
                hp = mc.opposite_halfface_handle(hp);

            auto dir2orderedHas = halfpatchHalfarcsByDir(hp);
            assert(dir2orderedHas.size() == 4);
            UVWDir dirHpNormal = halfpatchNormalDir(hp);

            for (UVWDir dirSide :
                 decompose(~(dirHpNormal | -dirHpNormal), {UVWDir::POS_U, UVWDir::POS_V, UVWDir::POS_W}))
            {
                UVWDir dirSideOpp = -dirSide;

                GRBLinExpr side = 0;
                for (HEH ha : dir2orderedHas.at(dirSide))
                    side += arc2var.at(mc.edge_handle(ha));

                GRBLinExpr sideOpp = 0;
                for (HEH ha : dir2orderedHas.at(-dirSide))
                    sideOpp += arc2var.at(mc.edge_handle(ha));

                std::string constraintName = "Patch" + std::to_string(patch.idx()) + "dir"
                                             + std::to_string(static_cast<uint8_t>(dirSide | dirSideOpp));
                model.addConstr(side, GRB_EQUAL, sideOpp, constraintName);
            }
        }

        if (allowNegative)
        {
            // Enforce no block with negative extension along any axis
            for (CH b : mc.cells())
            {
                for (UVWDir dir : {UVWDir::POS_V_POS_W, UVWDir::POS_U_POS_W, UVWDir::POS_U_POS_V})
                {
                    GRBLinExpr sumExpr = 0;
                    auto& arcs = mcMeshProps().ref<BLOCK_EDGE_ARCS>(b).at(dir);
                    for (EH a : arcs)
                        sumExpr += arc2var.at(a);
                    std::string constraintName = "Block" + std::to_string(b.idx()) + "dir"
                                                 + std::to_string(static_cast<uint8_t>(~(dir | -dir)));
                    model.addConstr(sumExpr, GRB_GREATER_EQUAL, 0.0, constraintName);
                }
            }
        }

        // Enforce critical links of length >= 1
        vector<CriticalLink> criticalLinks;
        map<EH, int> a2criticalLinkIdx;
        map<VH, vector<int>> n2criticalLinksOut;
        map<VH, vector<int>> n2criticalLinksIn;

        int nSeparationConstraints = 0;
        getCriticalLinks(criticalLinks, a2criticalLinkIdx, n2criticalLinksOut, n2criticalLinksIn, true);

        for (auto& criticalLink : criticalLinks)
        {
            if (criticalLink.pathHas.empty())
                continue;
            GRBLinExpr sumExpr = 0;
            for (HEH ha : criticalLink.pathHas)
                sumExpr += arc2var.at(mc.edge_handle(ha));

            std::string constraintName = "Separation" + std::to_string(nSeparationConstraints);
            if (criticalLink.nFrom == criticalLink.nTo)
                model.addConstr(sumExpr, GRB_GREATER_EQUAL, 3.0, constraintName);
            else
                model.addConstr(sumExpr, GRB_GREATER_EQUAL, 1.0, constraintName);
            nSeparationConstraints++;
        }

        // Additional code to avoid selfadjacency
        {
            for (CH b : mc.cells())
            {
                auto& dir2ps = mcMeshProps().ref<BLOCK_FACE_PATCHES>(b);
                for (auto& kv : dir2ps)
                {
                    if (isNeg(kv.first))
                        continue;

                    UVWDir dir = kv.first;
                    auto& ps = kv.second;
                    auto& psOpp = dir2ps.at(-dir);
                    set<CH> bs, bsOpp;
                    for (FH p : ps)
                    {
                        HFH hp = mc.halfface_handle(p, 0);
                        if (mc.incident_cell(hp) == b)
                            hp = mc.opposite_halfface_handle(hp);
                        CH bNext = mc.incident_cell(hp);
                        if (!bNext.is_valid())
                            continue;
                        bs.insert(bNext);
                    }
                    for (FH p : psOpp)
                    {
                        HFH hp = mc.halfface_handle(p, 0);
                        if (mc.incident_cell(hp) == b)
                            hp = mc.opposite_halfface_handle(hp);
                        CH bNext = mc.incident_cell(hp);
                        if (!bNext.is_valid())
                            continue;
                        bsOpp.insert(bNext);
                    }
                    bool mustBeNonZero = bs.count(b) != 0 || bsOpp.count(b) != 0
                                         || containsSomeOf(bs, bsOpp);
                    if (mustBeNonZero)
                    {
                        // Constrain non-zero length along dir
                        GRBLinExpr sumExpr = 0;
                        auto& dir2as = mcMeshProps().ref<BLOCK_EDGE_ARCS>(b);
                        for (EH a : dir2as.at(decompose(~(dir | -dir), DIM_2_DIRS)[0]))
                            sumExpr += arc2var.at(a);
                        std::string constraintName = "Selfadjacency" + std::to_string(nSeparationConstraints);
                        model.addConstr(sumExpr, GRB_GREATER_EQUAL, 1.0, constraintName);
                        nSeparationConstraints++;

                        continue;
                    }
                }
            }
        }

        int iter = 0;
        // Iteratively enforce violated separation constraints
        bool validSolution = false;
        while (!validSolution)
        {
            DLOG(INFO) << "Gurobi solving with " << nSeparationConstraints << " nonzero sum constraints";
            model.optimize();
            auto obj = model.getObjective();
            DLOG(INFO) << "Solved with final objective value of " << obj.getValue();
            iter++;

            int status = model.get(GRB_IntAttr_Status);
            if (status != 2 && status != 9)
            {
                LOG(ERROR) << "Bad status return by GUROBI solver";
                return SOLVER_ERROR;
            }
            mcMeshProps().allocate<ARC_INT_LENGTH>();
            for (EH arc : mc.edges())
            {
                auto& grbvar = arc2var.at(arc);
                mcMeshProps().set<ARC_INT_LENGTH>(arc, (int)std::round(grbvar.get(GRB_DoubleAttr_X)));
            }

            int minArcLength = 10000;
            for (EH a : mc.edges())
                minArcLength = std::min(minArcLength, mcMeshProps().get<ARC_INT_LENGTH>(a));
            DLOG(INFO) << "Min quantized arc length " << minArcLength;

            for (auto& criticalLink : criticalLinks)
            {
                criticalLink.length = 0;
                for (HEH ha : criticalLink.pathHas)
                    criticalLink.length += mcMeshProps().get<ARC_INT_LENGTH>(mcMeshProps().mesh().edge_handle(ha));
            }

            if (!allowZero)
                validSolution = true;
            else
            {
                vector<bool> isCriticalArc(mc.n_edges(), false);
                vector<bool> isCriticalNode(mc.n_vertices(), false);
                vector<bool> isCriticalPatch(mc.n_faces(), false);
                for (auto& kv : a2criticalLinkIdx)
                    isCriticalArc[kv.first.idx()] = true;
                for (VH n : mc.vertices())
                {
                    auto type = mcMeshProps().nodeType(n);
                    if (type.first == SingularNodeType::SINGULAR || type.second == FeatureNodeType::FEATURE
                        || type.second == FeatureNodeType::SEMI_FEATURE_SINGULAR_BRANCH)
                        isCriticalNode[n.idx()] = true;
                }
                for (FH p : mc.faces())
                    isCriticalPatch[p.idx()]
                        = mc.is_boundary(p)
                          || (mcMeshProps().isAllocated<IS_FEATURE_F>() && mcMeshProps().get<IS_FEATURE_F>(p));

                vector<vector<pair<int, EH>>> forcedNonZeroSum;
                findSeparationViolatingPaths(
                    criticalLinks, isCriticalArc, isCriticalNode, isCriticalPatch, forcedNonZeroSum);
                validSolution = forcedNonZeroSum.size() == 0;
                for (auto& aColl : forcedNonZeroSum)
                {
                    assert(!aColl.empty());
                    vector<GRBVar> vars;
                    vector<double> coeff;
                    for (auto sign2a : aColl)
                    {
                        vars.push_back(arc2var.at(sign2a.second));
                        coeff.push_back(sign2a.first);
                    }
                    GRBLinExpr sumExpr = 0;
                    sumExpr.addTerms(&coeff[0], &vars[0], vars.size());

                    std::string constraintName = "Separation" + std::to_string(nSeparationConstraints);
                    model.addConstr(sumExpr, GRB_GREATER_EQUAL, 0.9, constraintName);
                    nSeparationConstraints++;
                }
                if (!validSolution)
                {
                    int nHexes = numHexesInQuantization();
                    DLOG(INFO) << "Invalid intermediate solution has " << nHexes << " hexes";
                }
            }
        }
        int minArcLength = 10000;
        for (EH a : mc.edges())
            minArcLength = std::min(minArcLength, mcMeshProps().get<ARC_INT_LENGTH>(a));
        int nHexes = numHexesInQuantization();
        LOG(INFO) << "Quantized MC, nHexes: " << nHexes << ", objective: " << model.getObjective().getValue()
                  << ", iterations: " << iter << ", minArcLength " << minArcLength;
    }
    catch (GRBException e)
    {
        LOG(ERROR) << "Gurobi exception, errcode: " << e.getErrorCode();
        LOG(ERROR) << "Gurobi error message: " << e.getMessage();
        return SOLVER_ERROR;
    }

    return SUCCESS;
}

int MCQuantizer::numHexesInQuantization() const
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

MCQuantizer::RetCode MCQuantizer::findSeparationViolatingPaths(const vector<CriticalLink>& criticalLinks,
                                                               const vector<bool>& arcIsCritical,
                                                               const vector<bool>& nodeIsCritical,
                                                               const vector<bool>& patchIsCritical,
                                                               vector<vector<pair<int, EH>>>& nonZeroSumArcs) const
{
    nonZeroSumArcs.clear();

    // For each critical link s1 find paths connecting s1 to other critical links s2 or surface patches p2
    // Then check for overlaps between these, accumulating the quantized edge lengths along the path as deltas.
    for (auto& criticalLink : criticalLinks)
    {
        auto ret = traceExhaustPaths(criticalLink, arcIsCritical, nodeIsCritical, patchIsCritical, nonZeroSumArcs);
        if (ret != SUCCESS)
            return ret;
    }

    return SUCCESS;
}

MCQuantizer::RetCode MCQuantizer::traceExhaustPaths(const CriticalLink& criticalLink1,
                                                    const vector<bool>& arcIsCritical,
                                                    const vector<bool>& nodeIsCritical,
                                                    const vector<bool>& patchIsCritical,
                                                    vector<vector<pair<int, EH>>>& nonZeroSumArcs) const
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
    return SUCCESS;
}

bool MCQuantizer::bboxOverlap(const Vec3i& min1, const Vec3i& max1, const Vec3i& min2, const Vec3i& max2)
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

MCQuantizer::RetCode MCQuantizer::determineNextHalfedges(const WeaklyMonotonousPath& pathCurrent,
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
    return SUCCESS;
}

bool MCQuantizer::checkArcOverlap(const WeaklyMonotonousPath& pathCurrent,
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

bool MCQuantizer::checkPatchOverlap(const WeaklyMonotonousPath& pathCurrent,
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

bool MCQuantizer::checkP0Containment(const WeaklyMonotonousPath& pathCurrent) const
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

} // namespace mc3d
