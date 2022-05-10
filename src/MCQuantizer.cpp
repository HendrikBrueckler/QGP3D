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
    const MCMesh& mc = _mcMeshProps.mesh;
    const TetMesh& mesh = _meshProps.mesh;

    try
    {
        GRBEnv env = GRBEnv(true);
#ifdef NDEBUG
        env.set(GRB_IntParam_LogToConsole, false);
#endif
        env.set(GRB_DoubleParam_TimeLimit, 5 * 60);
        env.start();

        GRBModel model = GRBModel(env);

        _mcMeshProps.allocate<ARC_DBL_LENGTH>(0.0);
        // One integer length variable per arc
        std::map<OVM::EdgeHandle, GRBVar> arc2var;
        for (auto arc : mc.edges())
        {
            // Determine current arc length
            double length = 0.0;
            for (auto he : _mcMeshProps.ref<ARC_MESH_HALFEDGES>(arc))
                length += edgeLengthUVW(mesh.edge_handle(he));
            _mcMeshProps.set<ARC_DBL_LENGTH>(arc, length);

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
        for (auto block : mc.cells())
        {
            // For each axis U/V/W
            for (auto axis : {UVWDir::NEG_V_NEG_W, UVWDir::NEG_U_NEG_W, UVWDir::NEG_U_NEG_V})
            {
                GRBLinExpr varSum = 0.;
                double lenSum = 0;
                for (auto arc : _mcMeshProps.ref<BLOCK_EDGE_ARCS>(block).at(axis))
                {
                    varSum += arc2var.at(arc);
                    lenSum += _mcMeshProps.get<ARC_DBL_LENGTH>(arc);
                }
                objective += (varSum - scaling * lenSum) * (varSum - scaling * lenSum);
            }
        }

        model.setObjective(objective, GRB_MINIMIZE);

        // Each patches opposite arc lengths must match
        for (auto patch : mc.faces())
        {
            auto halfpatch = mc.halfface_handle(patch, 0);
            if (mc.is_boundary(halfpatch))
                halfpatch = mc.opposite_halfface_handle(halfpatch);

            auto dir2orderedHas = halfpatchHalfarcsByDir(halfpatch);

            assert(dir2orderedHas.size() == 4);
            map<UVWDir, bool> checked;
            for (const auto& kv : dir2orderedHas)
            {
                auto dir = kv.first;
                auto halfarcs = kv.second;
                assert(!halfarcs.empty());
                assert(dir2orderedHas.find(-dir) != dir2orderedHas.end());
                if (checked[dir])
                    continue;

                std::vector<GRBVar> vars1, vars2;
                std::vector<double> coeff1, coeff2;
                for (auto halfarc : halfarcs)
                {
                    auto arc = mc.edge_handle(halfarc);
                    vars1.push_back(arc2var.at(arc));
                    coeff1.push_back(1.);
                }
                checked[dir] = true;

                for (auto halfarc : dir2orderedHas.at(-dir))
                {
                    auto arc = mc.edge_handle(halfarc);
                    vars2.push_back(arc2var.at(arc));
                    coeff2.push_back(1.);
                }
                checked[-dir] = true;

                GRBLinExpr edge1 = 0;
                edge1.addTerms(&coeff1[0], &vars1[0], vars1.size());

                GRBLinExpr edge2 = 0;
                edge2.addTerms(&coeff2[0], &vars2[0], vars2.size());

                std::string constraintName
                    = "Patch" + std::to_string(patch.idx()) + "dir" + std::to_string(static_cast<uint8_t>(dir));
                model.addConstr(edge1, GRB_EQUAL, edge2, constraintName);
            }
            assert(checked.size() == 4);
        }

        if (allowNegative)
        {
            // Enforce no block with negative extension along any axis
            for (auto b : mc.cells())
            {
                int i = 0;
                for (auto dir : {UVWDir::POS_V_POS_W, UVWDir::POS_U_POS_W, UVWDir::POS_U_POS_V})
                {
                    auto& arcs = _mcMeshPropsC.ref<BLOCK_EDGE_ARCS>(b).at(dir);
                    std::vector<GRBVar> vars;
                    std::vector<double> coeff;
                    for (auto a : arcs)
                    {
                        vars.emplace_back(arc2var.at(a));
                        coeff.emplace_back(1.0);
                    }
                    GRBLinExpr sumExpr = 0;
                    sumExpr.addTerms(&coeff[0], &vars[0], vars.size());

                    std::string constraintName = "BlockLength" + std::to_string(b.idx()) + "dir" + std::to_string(i);
                    model.addConstr(sumExpr, GRB_GREATER_EQUAL, 0.0, constraintName);
                    i++;
                }
            }
        }

        // Enforce singular links of length >= 1
        vector<SingularLink> singularLinks;
        map<OVM::EdgeHandle, int> aSing2singularLinkIdx;
        map<OVM::VertexHandle, vector<int>> n2singularLinksOut;
        map<OVM::VertexHandle, vector<int>> n2singularLinksIn;

        int nSeparationConstraints = 0;
        getSingularLinks(singularLinks, aSing2singularLinkIdx, n2singularLinksOut, n2singularLinksIn);

        for (auto& singPath : singularLinks)
        {
            std::vector<GRBVar> vars;
            std::vector<double> coeff;
            for (auto ha : singPath.pathHas)
            {
                auto arc = mc.edge_handle(ha);
                vars.push_back(arc2var.at(arc));
                coeff.push_back(1.);
            }
            GRBLinExpr sumExpr = 0;
            sumExpr.addTerms(&coeff[0], &vars[0], vars.size());

            std::string constraintName = "Separation" + std::to_string(nSeparationConstraints);
            if (singPath.cyclic)
                model.addConstr(sumExpr, GRB_GREATER_EQUAL, 2.9, constraintName);
            else
                model.addConstr(sumExpr, GRB_GREATER_EQUAL, 0.9, constraintName);
            nSeparationConstraints++;
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
            _mcMeshProps.allocate<ARC_INT_LENGTH>();
            for (auto arc : mc.edges())
            {
                auto& grbvar = arc2var.at(arc);
                _mcMeshProps.set<ARC_INT_LENGTH>(arc, (int)std::round(grbvar.get(GRB_DoubleAttr_X)));
            }

            int minArcLength = 10000;
            for (auto a : mc.edges())
                minArcLength = std::min(minArcLength, _mcMeshPropsC.get<ARC_INT_LENGTH>(a));
            DLOG(INFO) << "Min quantized arc length " << minArcLength;

            for (auto& singPath : singularLinks)
            {
                singPath.length = 0;
                for (auto ha : singPath.pathHas)
                    singPath.length += _mcMeshPropsC.get<ARC_INT_LENGTH>(_mcMeshProps.mesh.edge_handle(ha));
            }

            if (!allowZero)
                validSolution = true;
            else
            {
                vector<vector<std::pair<int, OVM::EdgeHandle>>> forcedNonZeroSum;
                findSeparationViolatingPaths(singularLinks, forcedNonZeroSum);
                validSolution = forcedNonZeroSum.size() == 0;
                for (auto& aColl : forcedNonZeroSum)
                {
                    assert(!aColl.empty());
                    std::vector<GRBVar> vars;
                    std::vector<double> coeff;
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
        for (auto a : mc.edges())
            minArcLength = std::min(minArcLength, _mcMeshPropsC.get<ARC_INT_LENGTH>(a));
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
    const MCMesh& mc = _mcMeshProps.mesh;
    int nHexes = 0;
    for (auto b : mc.cells())
    {
        int nBlockHexes = 1;
        for (auto dir : {UVWDir::NEG_U_NEG_V, UVWDir::NEG_U_NEG_W, UVWDir::NEG_V_NEG_W})
        {
            int arcLen = 0;
            for (auto a : _mcMeshPropsC.ref<BLOCK_EDGE_ARCS>(b).at(dir))
                arcLen += _mcMeshPropsC.get<ARC_INT_LENGTH>(a);
            nBlockHexes *= arcLen;
        }
        nHexes += nBlockHexes;
    }
    return nHexes;
}

MCQuantizer::RetCode
MCQuantizer::findSeparationViolatingPaths(const vector<SingularLink>& singularLinks,
                                          vector<vector<std::pair<int, OVM::EdgeHandle>>>& nonZeroSumArcs) const
{
    nonZeroSumArcs.clear();

    // For each singular link s1 find paths connecting s1 to other singular links s2 or surface patches p2
    // Then check for overlaps between these, accumulating the quantized edge lengths along the path as deltas.
    for (auto& singPath : singularLinks)
    {
        auto ret = traceExhaustPaths(singPath, nonZeroSumArcs);
        if (ret != SUCCESS)
            return ret;
    }

    return SUCCESS;
}

MCQuantizer::RetCode
MCQuantizer::traceExhaustPaths(const SingularLink& singPath1,
                               vector<vector<std::pair<int, OVM::EdgeHandle>>>& nonZeroSumArcs) const
{
    const MCMesh& mcMesh = _mcMeshProps.mesh;

    auto singStartHa = singPath1.pathHas.front();

    vector<vector<OVM::HalfEdgeHandle>> dir2has;
    // Categorize all vertical halfarcs by direction
    for (auto singStartHp : mcMesh.halfedge_halffaces(singStartHa))
    {
        dir2has.emplace_back();

        auto& has = dir2has.back();

        auto bRef = mcMesh.incident_cell(singStartHp);
        if (!bRef.is_valid())
            bRef = mcMesh.incident_cell(mcMesh.opposite_halfface_handle(singStartHp));

        auto dirStart = halfarcDirInBlock(singStartHa, bRef);

        auto haOrth = mcMesh.prev_halfedge_in_halfface(singStartHa, singStartHp);
        while (halfarcDirInBlock(haOrth, bRef) == dirStart)
            haOrth = mcMesh.prev_halfedge_in_halfface(haOrth, singStartHp);
        haOrth = mcMesh.opposite_halfedge_handle(haOrth);
        has.emplace_back(haOrth);

        auto haCurr = singStartHa;
        auto hpCurr = singStartHp;

        bool foundNext = true;
        do
        {
            bRef = mcMesh.incident_cell(hpCurr);
            if (!bRef.is_valid())
                bRef = mcMesh.incident_cell(mcMesh.opposite_halfface_handle(hpCurr));

            dirStart = halfarcDirInBlock(haCurr, bRef);

            haOrth = mcMesh.next_halfedge_in_halfface(haCurr, hpCurr);
            while (halfarcDirInBlock(haOrth, bRef) == dirStart)
                haOrth = mcMesh.next_halfedge_in_halfface(haOrth, hpCurr);
            has.emplace_back(haOrth);

            foundNext = false;
            auto haOrthOpp = mcMesh.opposite_halfedge_handle(haOrth);
            for (auto hpNext : mcMesh.halfedge_halffaces(haOrthOpp))
            {
                if (mcMesh.opposite_halfface_handle(hpNext) == hpCurr)
                    continue;
                auto haNext = mcMesh.next_halfedge_in_halfface(haOrthOpp, hpNext);
                if (_mcMeshPropsC.get<ARC_IS_SINGULAR>(mcMesh.edge_handle(haNext))
                    && nodeType(mcMesh.from_vertex_handle(haNext)) == NodeType::SEMI_SINGULAR
                    && mcMesh.from_vertex_handle(haNext) == mcMesh.to_vertex_handle(haCurr))
                {
                    haCurr = haNext;
                    hpCurr = hpNext;
                    foundNext = true;
                    break;
                }
            }
        } while (foundNext && haCurr != singStartHa);
    }

    map<OVM::HalfEdgeHandle, int> haOrth2dir;
    for (int i = 0; i < (int)dir2has.size(); i++)
        for (auto ha : dir2has[i])
            haOrth2dir[ha] = i;

    auto singStartN = singPath1.nFrom;
    auto bRefStart = *mcMesh.hec_iter(singStartHa);
    auto dirStartHa = halfarcDirInBlock(singStartHa, bRefStart);

    WeaklyMonotonousPath pStart;
    pStart.branchedOff = false;
    pStart.length = 0;
    pStart.singPathIdx1 = singPath1.id;
    pStart.n = singStartN;
    pStart.dirs1 = dirStartHa | -dirStartHa;
    pStart.path = {};
    pStart.monotonousDirs = ~pStart.dirs1; // Do not search along the singularity but orthogonally
    pStart.walkedDirs = UVWDir::NONE;
    pStart.deltaMin = (isNeg(dirStartHa) ? singPath1.length * toVec(dirStartHa) : Vec3i(0, 0, 0));
    pStart.deltaMax = (isNeg(dirStartHa) ? Vec3i(0, 0, 0) : singPath1.length * toVec(dirStartHa));
    if (singPath1.cyclic)
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
        // vector<UVWDir> nsVisitedByDirs(mcMesh.n_vertices(), UVWDir::NONE);
        vector<bool> nsVisited(mcMesh.n_vertices(), false);
        vector<bool> nsInitialized(mcMesh.n_vertices(), false);
        PathQueue pathQ;
        pathQ.push(pStart);
        while (!pathQ.empty())
        {
            auto currentP = pathQ.top();
            pathQ.pop();

            // if ((~nsVisitedByDirs[currentP.n.idx()] & currentP.monotonousDirs) == UVWDir::NONE)
            //     continue;
            // nsVisitedByDirs[currentP.n.idx()] = nsVisitedByDirs[currentP.n.idx()] | currentP.monotonousDirs;

            if ((currentP.branchedOff && nsVisited[currentP.n.idx()])
                || (!currentP.branchedOff && nsInitialized[currentP.n.idx()]))
                continue;
            if (currentP.branchedOff)
                nsVisited[currentP.n.idx()] = true;
            else
                nsInitialized[currentP.n.idx()] = true;

            map<OVM::CellHandle, Transition> b2trans;
            map<OVM::HalfEdgeHandle, vector<OVM::CellHandle>> ha2bRef;
            determineNextHalfedges(currentP, b2trans, ha2bRef);

            // Check for separation violations by current path
            if ((currentP.walkedDirs & currentP.monotonousDirs) != UVWDir::NONE)
            {
                UVWDir violationDir = UVWDir::NONE;
                // Check whether startinterval-singularnode pairs overlap
                if (nodeType(currentP.n) == NodeType::SINGULAR)
                {
                    int nDimFullOverlap = 0;
                    if (bboxOverlap(
                            currentP.deltaMin, currentP.deltaMax, currentP.delta, currentP.delta, nDimFullOverlap))
                    {
                        violationDir = currentP.walkedDirs & ~currentP.dirs1;
                        assert(violationDir != UVWDir::NONE);
                    }
                }

                if (violationDir == UVWDir::NONE)
                {
                    // Check whether startinterval-singularlink pairs overlap
                    for (const auto& kv : ha2bRef)
                    {
                        auto ha2 = kv.first;
                        auto a2 = mcMesh.edge_handle(ha2);
                        if (_mcMeshPropsC.get<ARC_IS_SINGULAR>(a2))
                        {
                            auto bRef2 = kv.second.front();
                            auto trans2 = b2trans.at(bRef2);

                            auto dir2 = trans2.invert().rotate(halfarcDirInBlock(ha2, bRef2));

                            auto possibleViolationDir = currentP.walkedDirs & ~(dir2 | -dir2) & ~currentP.dirs1;
                            if ((possibleViolationDir & currentP.monotonousDirs) == UVWDir::NONE)
                                continue;

                            // Measure just overlap with the arc, not the singular link
                            int length = _mcMeshPropsC.get<ARC_INT_LENGTH>(a2);

                            if (length < 0)
                            {
                                length = -length;
                                dir2 = -dir2;
                            }
                            Vec3i deltaMin2 = currentP.delta + (isNeg(dir2) ? length * toVec(dir2) : Vec3i(0, 0, 0));
                            Vec3i deltaMax2 = currentP.delta + (isNeg(dir2) ? Vec3i(0, 0, 0) : length * toVec(dir2));

                            // Check for overlap
                            int nDimFullOverlap = 0;
                            if (bboxOverlap(
                                    currentP.deltaMin, currentP.deltaMax, deltaMin2, deltaMax2, nDimFullOverlap))
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
                    if (mcMesh.is_boundary(currentP.n))
                    {
                        // Check for overlaps between startinterval-surfacepatch pairs
                        for (auto hp2 : mcMesh.vertex_halffaces(currentP.n))
                        {
                            if (mcMesh.is_boundary(hp2))
                            {
                                auto hp2opp = mcMesh.opposite_halfface_handle(hp2);
                                auto bRef2 = mcMesh.incident_cell(hp2opp);

                                if (b2trans.find(bRef2) == b2trans.end())
                                    continue;
                                auto trans2 = b2trans.at(bRef2);

                                auto hpDirs = UVWDir::NONE;
                                for (auto ha : mcMesh.halfface_halfedges(hp2opp))
                                    hpDirs = hpDirs | halfarcDirInBlock(ha, bRef2);

                                auto hpDirsLocal = trans2.invert().rotate(hpDirs);
                                assert(dim(hpDirsLocal) == 2);

                                auto possibleViolationDir = currentP.walkedDirs & ~hpDirsLocal & ~currentP.dirs1;
                                if ((possibleViolationDir & currentP.monotonousDirs) == UVWDir::NONE)
                                    continue;

                                if (checkPatchOverlap(currentP, hp2opp, trans2))
                                {
                                    assert(possibleViolationDir != UVWDir::NONE);
                                    violationDir = possibleViolationDir;
                                    break;
                                }
                            }
                        }
                    }
                }

                if (violationDir != UVWDir::NONE)
                {
                    vector<OVM::EdgeHandle> posSignArcs;
                    vector<OVM::EdgeHandle> negSignArcs;
                    for (auto dim1dirPos : {UVWDir::POS_U, UVWDir::POS_V, UVWDir::POS_W})
                    {
                        auto dim1dirAny = dim1dirPos | -dim1dirPos;
                        auto intersection = violationDir & dim1dirAny;
                        if (intersection != UVWDir::NONE)
                        {
                            if (intersection == dim1dirAny)
                            {
                                const auto& posArcs = currentP.dir2walkedArcs.at(dim1dirPos);
                                const auto& negArcs = currentP.dir2walkedArcs.at(-dim1dirPos);
                                assert(posArcs.size() > 0);
                                assert(negArcs.size() > 0);
                                double lengthPos = 0;
                                double lengthNeg = 0;
                                for (auto a : posArcs)
                                    lengthPos += _mcMeshPropsC.get<ARC_DBL_LENGTH>(a);
                                for (auto a : negArcs)
                                    lengthNeg += _mcMeshPropsC.get<ARC_DBL_LENGTH>(a);
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
                                auto& arcs = currentP.dir2walkedArcs.at(intersection);
                                posSignArcs.insert(posSignArcs.end(), arcs.begin(), arcs.end());
                            }
                        }
                    }
                    assert(posSignArcs.size() > 0);
                    if (posSignArcs.empty())
                        continue;
                    nonZeroSumArcs.emplace_back();
                    auto& nonZeroSum = nonZeroSumArcs.back();
                    for (auto a : posSignArcs)
                        nonZeroSum.push_back({1, a});
                    for (auto a : negSignArcs)
                        nonZeroSum.push_back({-1, a});

                    assert(currentP.branchedOff);

                    // return SUCCESS;
                    // LOG(INFO) << "S " << startSingPath << " added constraint in dir " << dirIdx;
                    break;
                }
            }

            for (const auto& kv : ha2bRef)
            {
                auto ha = kv.first;
                auto nTo = mcMesh.to_vertex_handle(ha);
                if (currentP.branchedOff && nsVisited[nTo.idx()])
                    continue;
                // if ((~nsVisitedByDirs[nTo.idx()] & currentP.monotonousDirs) == UVWDir::NONE)
                //     continue;

                for (auto bRef : kv.second)
                {
                    auto& trans = b2trans.at(bRef);
                    WeaklyMonotonousPath nextP = currentP;
                    nextP.bRefCurrent = bRef;
                    nextP.transCurrent = trans;

                    if (!checkP0Containment(currentP))
                        continue;

                    nextP.n = nTo;
                    nextP.path.emplace_back(ha);

                    auto walkedDir = trans.invert().rotate(halfarcDirInBlock(ha, nextP.bRefCurrent));

                    // Record branch-off from singularity, allow only limited set of edges
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

                    nextP.delta += _mcMeshPropsC.get<ARC_INT_LENGTH>(mcMesh.edge_handle(ha)) * toVec(walkedDir);

                    nextP.walkedDirs = nextP.walkedDirs | walkedDir;
                    nextP.monotonousDirs = nextP.monotonousDirs & ~(-walkedDir);

                    // Do not walk unnecessary circles
                    if (nextP.monotonousDirs == UVWDir::NONE)
                        continue;

                    nextP.dir2walkedArcs[walkedDir].emplace_back(mcMesh.edge_handle(ha));

                    // Walk along the singularity for free, but other arcs accumulate distance
                    // if (nextP.branchedOff)
                    if ((walkedDir & nextP.dirs1) == UVWDir::NONE)
                        nextP.length += _mcMeshPropsC.get<ARC_DBL_LENGTH>(mcMesh.edge_handle(ha));

                    pathQ.push(nextP);
                    break; // Only push the edge once (still need to check each bRef)
                }
            }
        }
    }
    return SUCCESS;
}

bool MCQuantizer::bboxOverlap(
    const Vec3i& min1, const Vec3i& max1, const Vec3i& min2, const Vec3i& max2, int& nDimFullOverlap)
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

        // min2[i] is in (min1[i], max1[i])
        // OR max2[i] is in (min1[i], max1[i])
        // OR min1[i] is in (min2[i], max2[i])
        // OR max1[i] is in (min2[i], max2[i])
        nDimFullOverlap += (min1[i] < min2[i] && min2[i] < max1[i]) || (min1[i] < max2[i] && max2[i] < max1[i])
                           || (min2[i] < min1[i] && min1[i] < max2[i]) || (min2[i] < max1[i] && max1[i] < max2[i]);
    }

    return true;
}

MCQuantizer::RetCode
MCQuantizer::determineNextHalfedges(const WeaklyMonotonousPath& currentP,
                                    map<OVM::CellHandle, Transition>& b2trans,
                                    map<OVM::HalfEdgeHandle, vector<OVM::CellHandle>>& ha2bRef) const
{
    auto& mcMesh = _mcMeshPropsC.mesh;

    b2trans = map<OVM::CellHandle, Transition>({{currentP.bRefCurrent, {currentP.transCurrent}}});

    // Floodfill blocks around n, storing current transition for each expanded block
    list<std::pair<OVM::CellHandle, Transition>> bQ({{currentP.bRefCurrent, currentP.transCurrent}});

    while (!bQ.empty())
    {
        auto b2t = bQ.front();
        bQ.pop_front();

        for (auto hp : mcMesh.cell_halffaces(b2t.first))
        {
            auto hpOpp = mcMesh.opposite_halfface_handle(hp);
            auto bNext = mcMesh.incident_cell(hpOpp);
            if (!bNext.is_valid() || b2trans.find(bNext) != b2trans.end())
                continue;
            bool hasN = false;
            for (auto n2 : mcMesh.halfface_vertices(hp))
                if (n2 == currentP.n)
                {
                    hasN = true;
                    break;
                }
            if (!hasN)
                continue;

            // Check if overlapping with current path
            if (!checkPatchOverlap(currentP, hp, b2t.second))
                continue;

            auto trans = b2t.second.chain(_mcMeshPropsC.hpTransition(hp));

            bool exists = b2trans.find(bNext) != b2trans.end();
            if (exists)
                continue;
            b2trans[bNext] = trans;

            bQ.push_back({bNext, trans});
        }
    }
    ha2bRef.clear();
    for (auto ha : mcMesh.outgoing_halfedges(currentP.n))
    {
        if (!currentP.path.empty() && currentP.path.back() == mcMesh.opposite_halfedge_handle(ha))
            continue;
        for (auto b : mcMesh.halfedge_cells(ha))
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

bool MCQuantizer::checkArcOverlap(const WeaklyMonotonousPath& currentP,
                                  const OVM::HalfEdgeHandle& ha,
                                  const OVM::CellHandle& bRef,
                                  const Transition& trans) const
{
    auto& mcMesh = _mcMeshPropsC.mesh;

    auto dir2 = trans.invert().rotate(halfarcDirInBlock(ha, bRef));
    int arcLen = _mcMeshPropsC.get<ARC_INT_LENGTH>(mcMesh.edge_handle(ha));

    Vec3i deltaMin2 = isNeg(dir2) ? currentP.delta + arcLen * toVec(dir2) : currentP.delta;
    Vec3i deltaMax2 = isNeg(dir2) ? currentP.delta : currentP.delta + arcLen * toVec(dir2);

    // Check for overlap
    int nDimFullOverlap = 0;
    return bboxOverlap(currentP.deltaMin, currentP.deltaMax, deltaMin2, deltaMax2, nDimFullOverlap);
}

bool MCQuantizer::checkPatchOverlap(const WeaklyMonotonousPath& currentP,
                                    const OVM::HalfFaceHandle& hp,
                                    const Transition& trans) const
{
    auto& mcMesh = _mcMeshPropsC.mesh;

    auto bRef2 = mcMesh.incident_cell(hp);

    Vec3i deltaMin2 = currentP.delta;
    Vec3i deltaMax2 = currentP.delta;

    OVM::HalfEdgeHandle ha1;
    for (auto haOther : mcMesh.halfface_halfedges(hp))
        if (mcMesh.from_vertex_handle(haOther) == currentP.n)
        {
            ha1 = haOther;
            break;
        }
    assert(ha1.is_valid());
    auto haCurr = ha1;
    auto deltaCurr = currentP.delta;

    auto transInv = trans.invert();
    do
    {
        auto dir = transInv.rotate(halfarcDirInBlock(haCurr, bRef2));
        auto aCurr = mcMesh.edge_handle(haCurr);
        auto haLength = _mcMeshPropsC.get<ARC_INT_LENGTH>(aCurr);
        deltaCurr += haLength * toVec(dir);
        for (int coord = 0; coord < 3; coord++)
        {
            deltaMax2[coord] = std::max(deltaMax2[coord], deltaCurr[coord]);
            deltaMin2[coord] = std::min(deltaMin2[coord], deltaCurr[coord]);
        }
        haCurr = mcMesh.next_halfedge_in_halfface(haCurr, hp);
    } while (haCurr != ha1);

    assert(deltaCurr == currentP.delta);
    assert(deltaMin2[0] <= deltaMax2[0] && deltaMin2[1] <= deltaMax2[1] && deltaMin2[2] <= deltaMax2[2]);

    int nDimFullOverlap = 0;
    return bboxOverlap(currentP.deltaMin, currentP.deltaMax, deltaMin2, deltaMax2, nDimFullOverlap);
}

bool MCQuantizer::checkP0Containment(const WeaklyMonotonousPath& currentP) const
{
    auto& mcMesh = _mcMeshPropsC.mesh;

    auto bRef = currentP.bRefCurrent;

    // Find the closest corner of bRef

    // Graph search
    list<std::pair<OVM::VertexHandle, Vec3i>> nQ;
    nQ.push_back({currentP.n, Vec3i(0, 0, 0)});
    map<OVM::VertexHandle, Vec3i> n2displacement({{currentP.n, Vec3i(0, 0, 0)}});
    while (!nQ.empty())
    {
        auto n = nQ.front().first;
        auto displacement = nQ.front().second;
        nQ.pop_front();

        for (auto ha : mcMesh.outgoing_halfedges(n))
        {
            bool inBRef = false;
            for (auto b : mcMesh.halfedge_cells(ha))
                if (b == bRef)
                {
                    inBRef = true;
                    break;
                }
            if (inBRef)
            {
                auto nTo = mcMesh.to_vertex_handle(ha);
                if (n2displacement.find(nTo) == n2displacement.end())
                {
                    auto dirHa = halfarcDirInBlock(ha, bRef);
                    auto length = _mcMeshPropsC.get<ARC_INT_LENGTH>(mcMesh.edge_handle(ha));
                    auto newDisplacement = displacement + toVec(dirHa) * length;
                    n2displacement[nTo] = newDisplacement;
                    nQ.push_back({nTo, newDisplacement});
                }
            }
        }
    }

    auto& dir2n = _mcMeshPropsC.ref<BLOCK_CORNER_NODES>(bRef);
    auto nMin = dir2n.at(UVWDir::NEG_U_NEG_V_NEG_W);
    auto nMax = dir2n.at(UVWDir::POS_U_POS_V_POS_W);

    int dimFullOverlap = 0;
    return bboxOverlap(currentP.deltaMin,
                       currentP.deltaMax,
                       currentP.delta + n2displacement.at(nMin),
                       currentP.delta + n2displacement.at(nMax),
                       dimFullOverlap);
}

} // namespace qgp3d
