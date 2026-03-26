#include "QGP3D/ObjectiveBuilder.hpp"

namespace qgp3d
{

ObjectiveBuilder::ObjectiveBuilder(TetMeshProps& meshProps) : TetMeshNavigator(meshProps), MCMeshNavigator(meshProps)
{
}

QuadraticObjective ObjectiveBuilder::arcLengthDeviationObjective(double scaling) const
{
    auto& mcMesh = mcMeshProps().mesh();
    double f0 = 0.0;
    Eigen::VectorXd grad0(mcMesh.n_logical_edges());
    Eigen::SparseMatrix<double> hess(mcMesh.n_logical_edges(), mcMesh.n_logical_edges());
    vector<Eigen::Triplet<double>> triplets;
    int i = 0;
    for (EH a : mcMesh.edges())
    {
        double l = mcMeshProps().get<ARC_DBL_LENGTH>(a);
        f0 += scaling * scaling * l * l;
        grad0(i) = -2 * l * scaling;
        triplets.push_back({i, i, 2});
        i++;
    }
    hess.setFromTriplets(triplets.begin(), triplets.end());

    return QuadraticObjective(f0, grad0, hess);
}

QuadraticObjective ObjectiveBuilder::simplifiedDistortionObjective(double scaling,
                                                                   double individualArcRegularization) const
{
    auto& mcMesh = mcMeshProps().mesh();
    double f0 = 0.0;
    Eigen::VectorXd grad0(mcMesh.n_logical_edges());
    Eigen::SparseMatrix<double> hess(mcMesh.n_logical_edges(), mcMesh.n_logical_edges());

    {
        set<VH> seeds;
        getCriticalSeeds(seeds);
        LOG(INFO) << "Critical seeds: " << seeds.size();
        map<VH, HEH> n2precursors;
        growVoronoiCellsMC(seeds, n2precursors);
        map<pairTT<VH>, vector<set<EH>>> seedPairs2interfaces;
        groupInterfacesMC(seeds, n2precursors, seedPairs2interfaces);
        vector<vector<HEH>> paths;
        vector<double> pathWeights; // Normalized (sum over all paths = 1), will denormalize if split into components
        getInterfaceTraversingPathsMC(n2precursors, seedPairs2interfaces, paths, pathWeights);

        map<EH, int> a2idx;
        for (EH a : mcMesh.edges())
            a2idx[a] = a2idx.size();

        vector<map<EH, int>> pathComponents;
        vector<double> pathComponentWeights;
        // Get block path
        for (int i = 0; i < (int)paths.size(); i++)
        {
            Transition trans;
            auto& path = paths[i];
            double weight = pathWeights[i];
            VH nCurrent = mcMesh.from_vertex_handle(path.front());
            CH bCurrent = *mcMesh.hec_iter(path.front());

            vector<map<EH, int>> components(3);
            for (HEH ha : path)
            {
                // DETERMINE CURRENT COORDINATE WRT INITIAL BLOCK, I.E. TRANSFORMED BY INVERSE ACCUMULATED
                // TRANSITION
                auto b2trans = determineTransitionsAroundNode(nCurrent, bCurrent, trans);

                nCurrent = mcMesh.to_vertex_handle(ha);
                bCurrent = *mcMesh.hec_iter(ha);
                trans = b2trans.at(bCurrent).front();

                UVWDir dir = toDir(trans.invert().rotate(toVec(halfarcDirInBlock(ha, bCurrent))));
                int coord = toCoord(dir);
                int sign = isNeg(dir) ? -1 : 1;
                components[coord][mcMesh.edge_handle(ha)] = sign;
            }
            for (int coord = 0; coord < 3; coord++)
                if (components[coord].size() > 0)
                {
                    pathComponents.push_back(components[coord]);
                    pathComponentWeights.push_back(weight);
                }
        }

        Eigen::Map<Eigen::VectorXd> massDiag(pathComponentWeights.data(), pathComponentWeights.size());
        Eigen::DiagonalMatrix<double, Eigen::Dynamic> mass(massDiag);
        Eigen::SparseMatrix<double> S(mass.cols(), a2idx.size());
        vector<Eigen::Triplet<double>> tripletsS;
        for (int i = 0; i < (int)pathComponents.size(); i++)
            for (auto& kv : pathComponents[i])
                tripletsS.emplace_back(i, a2idx.at(kv.first), kv.second);
        S.setFromTriplets(tripletsS.begin(), tripletsS.end());

        Eigen::VectorXd l(a2idx.size());
        for (EH a : mcMesh.edges())
            l(a2idx.at(a)) = mcMeshProps().get<ARC_DBL_LENGTH>(a) * scaling;

        hess = 2 * S.transpose() * mass * S;
        grad0 = -(hess * l);
        f0 = 0.5 * l.dot(hess * l);
    }

    // This is already quasi-normalized (somewhere between average and 3*average path length deviation)
    double weightPathDeviation = (1.0 - individualArcRegularization);

    // Divide by nArcs to get the "average" quadratic arc length deviation energy
    double weightArcDeviation = individualArcRegularization / mcMesh.n_logical_edges();
    auto regObj = arcLengthDeviationObjective(scaling);

    auto obj = QuadraticObjective(weightPathDeviation * f0 + weightArcDeviation * regObj.f0,
                                  weightPathDeviation * grad0 + weightArcDeviation * regObj.grad0,
                                  weightPathDeviation * hess + weightArcDeviation * regObj.hess);

    return obj;
}

ObjectiveBuilder::RetCode ObjectiveBuilder::getCriticalSeeds(set<VH>& seeds) const
{
    seeds.clear();
    auto& mcMesh = mcMeshProps().mesh();

    // check if all critical entities have enough seeds, if not insert additional ones
    vector<bool> isCriticalNode, isCriticalArc, isCriticalPatch;
    vector<CriticalEntity> criticalEntities;
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
                        false,
                        true);

    for (VH n : mcMesh.vertices())
        if (isCriticalNode[n.idx()])
            seeds.insert(n);

    for (auto crit : criticalEntities)
    {
        if (crit.dim == 1)
        {
            if (crit.pathHas.empty())
                throw std::logic_error("Empty path");
            set<VH> nsSeed;
            for (HEH ha : crit.pathHas)
                for (VH n : mcMesh.halfedge_vertices(ha))
                    if (seeds.count(n))
                        nsSeed.insert(n);
            if (nsSeed.size() < 2)
            {
                if (!crit.nFrom.is_valid() || !crit.nTo.is_valid())
                    throw std::logic_error("Invalid start or end");
                nsSeed.insert(crit.nFrom);
                nsSeed.insert(crit.nTo);
                seeds.insert(crit.nFrom);
                seeds.insert(crit.nTo);
                if (nsSeed.size() == 1)
                {
                    // Find maximally distant node around cycle
                    double targetLength = 0.0;
                    for (auto ha : crit.pathHas)
                        targetLength += 0.5 * mcMeshProps().get<ARC_DBL_LENGTH>(mcMesh.edge_handle(ha));
                    double currentLength = 0.0;
                    for (auto ha : crit.pathHas)
                    {
                        double nextLength = currentLength + mcMeshProps().get<ARC_DBL_LENGTH>(mcMesh.edge_handle(ha));
                        if (currentLength <= targetLength && nextLength > targetLength)
                        {
                            if (nextLength - targetLength < targetLength - currentLength
                                && !nsSeed.count(mcMesh.to_vertex_handle(ha)))
                            {
                                nsSeed.insert(mcMesh.to_vertex_handle(ha));
                                seeds.insert(mcMesh.to_vertex_handle(ha));
                            }
                            else
                            {
                                nsSeed.insert(mcMesh.from_vertex_handle(ha));
                                seeds.insert(mcMesh.from_vertex_handle(ha));
                            }
                            break;
                        }
                    }
                }
            }
        }
        else if (crit.dim == 2)
        {
            if (crit.regionPs.empty())
                throw std::logic_error("Empty region");
            set<VH> nsSeed;
            set<VH> nsAll;
            set<EH> asAll;
            for (FH p : crit.regionPs)
            {
                for (VH n : mcMesh.face_vertices(p))
                {
                    nsAll.insert(n);
                    if (seeds.count(n))
                        nsSeed.insert(n);
                }
                for (EH a : mcMesh.face_edges(p))
                    asAll.insert(a);
            }
            if (nsSeed.size() < 4)
            {
                if (nsSeed.size() == 0)
                {
                    for (VH n : nsAll)
                    {
                        if (mcMeshProps().nodeType(n).first == SingularNodeType::SINGULAR)
                        {
                            nsSeed.insert(n);
                            seeds.insert(n);
                            break;
                        }
                    }
                }

                while (nsSeed.size() < 4)
                {
                    // Dijkstra from all current nodes
                    struct Path
                    {
                        VH nCurr = VH();
                        double length = 0.0;
                    };
                    auto pathCmp = [](const Path& p1, const Path& p2) { return p1.length > p2.length; };

                    using PathQueue = std::priority_queue<Path, std::deque<Path>, decltype(pathCmp)>;

                    map<VH, double> n2dist;
                    PathQueue pQ(pathCmp);
                    for (auto n : nsSeed)
                    {
                        Path p0;
                        p0.nCurr = n;
                        p0.length = 0.0;
                        pQ.push(p0);
                    }
                    while (!pQ.empty())
                    {
                        Path pIn = pQ.top();
                        pQ.pop();

                        if (n2dist.count(pIn.nCurr))
                            continue;

                        n2dist[pIn.nCurr] = pIn.length;

                        for (HEH haOut : mcMesh.outgoing_halfedges(pIn.nCurr))
                        {
                            VH nTo = mcMesh.to_vertex_handle(haOut);
                            if (n2dist.count(nTo) || nsSeed.count(nTo) || !asAll.count(mcMesh.edge_handle(haOut)))
                                continue;
                            Path pOut;
                            pOut.nCurr = nTo;
                            pOut.length = pIn.length + mcMeshProps().get<ARC_DBL_LENGTH>(mcMesh.edge_handle(haOut));
                            pQ.push(pOut);
                        }
                    }

                    // Find furthest
                    VH nFurthest;
                    double furthest = 0.0;
                    VH nFurthestSingular;
                    double furthestSingular = 0.0;
                    for (auto& kv : n2dist)
                    {
                        if (kv.second > furthest)
                        {
                            furthest = kv.second;
                            nFurthest = kv.first;
                        }
                        if (mcMeshProps().nodeType(kv.first).first == SingularNodeType::SINGULAR
                            && kv.second > furthestSingular)
                        {
                            furthestSingular = kv.second;
                            nFurthestSingular = kv.first;
                        }
                    }
                    if (furthestSingular > 0.0)
                    {
                        if (nsSeed.count(nFurthestSingular))
                            throw std::logic_error("Error, impossible");
                        nsSeed.insert(nFurthestSingular);
                        seeds.insert(nFurthestSingular);
                    }
                    else if (furthest > 0.0)
                    {
                        if (nsSeed.count(nFurthest))
                            throw std::logic_error("Error, also impossible");
                        nsSeed.insert(nFurthest);
                        seeds.insert(nFurthest);
                    }
                    else
                    {
                        VH nAny;
                        for (auto n : nsAll)
                            if (!nsSeed.count(n))
                            {
                                nAny = n;
                                break;
                            }
                        if (!nAny.is_valid())
                            LOG(WARNING) << "NOT ENOUGH NODES IN FEATURE REGION";
                        nsSeed.insert(nAny);
                        seeds.insert(nAny);
                    }
                    for (auto n : nsSeed)
                        if (!nsAll.count(n))
                            throw std::logic_error("Node reached that should be unreachable");
                }
            }
        }
    }

    return SUCCESS;
}

namespace
{

struct Path
{
    VH nSource;
    HEH haCurrent;
    CH bCurrent;

    double pathLength;
    Vec3d pathDelta;
};

struct PathCompare
{
    bool operator()(const Path& a, const Path& b) const
    {
        return a.pathLength > b.pathLength;
    }
};

} // namespace

ObjectiveBuilder::RetCode ObjectiveBuilder::growVoronoiCellsMC(const set<VH>& seeds, map<VH, HEH>& n2precursors) const
{
    auto& mcMesh = mcMeshProps().mesh();
    n2precursors.clear();

    using PathQueue = std::priority_queue<Path, std::deque<Path>, PathCompare>;
    PathQueue pQ;

    for (VH seed : seeds)
    {
        Path p = {};
        p.nSource = seed;
        p.haCurrent = HEH();
        p.bCurrent = *mcMesh.vc_iter(seed);
        p.pathLength = 0.0;
        p.pathDelta = Vec3d();
        pQ.push(p);
    }

    set<VH> nsVisited;
    while (!pQ.empty())
    {
        Path p = pQ.top();
        pQ.pop();

        VH nCurrent = p.haCurrent.is_valid() ? mcMesh.to_vertex_handle(p.haCurrent) : p.nSource;

        if (n2precursors.count(nCurrent))
            continue;
        n2precursors[nCurrent] = p.haCurrent;

        auto b2trans = determineTransitionsAroundNode(nCurrent, p.bCurrent, Transition());

        for (HEH haOut : mcMesh.outgoing_halfedges(nCurrent))
        {
            CH bRef = *mcMesh.hec_iter(haOut);
            Vec3d pathDelta = Vec3Q2d(b2trans.at(bRef).front().apply(Vec3Q(p.pathDelta)));
            VH nTo = mcMesh.to_vertex_handle(haOut);
            if (n2precursors.count(nTo))
                continue;
            Path pNew = {};
            pNew.nSource = pNew.nSource;
            pNew.haCurrent = haOut;
            pNew.bCurrent = bRef;
            pNew.pathLength = p.pathLength + mcMeshProps().get<ARC_DBL_LENGTH>(mcMesh.edge_handle(haOut));
            pNew.pathDelta = pathDelta
                             + Vec3d(toVec(halfarcDirInBlock(haOut, bRef)))
                                   * mcMeshProps().get<ARC_DBL_LENGTH>(mcMesh.edge_handle(haOut));
            pQ.push(pNew);
        }
    }

    return SUCCESS;
}

ObjectiveBuilder::RetCode
ObjectiveBuilder::groupInterfacesMC(const set<VH>& seeds,
                                    const map<VH, HEH>& n2precursors,
                                    map<pairTT<VH>, vector<set<EH>>>& seedPairs2interfaces) const
{
    auto& mcMesh = mcMeshProps().mesh();
    seedPairs2interfaces.clear();
    if (seeds.empty())
        return SUCCESS;

    map<VH, VH> n2seed;
    for (VH n : mcMesh.vertices())
    {
        if (n2seed.count(n))
            continue;
        vector<VH> ns = {{n}};
        HEH precursor = n2precursors.at(ns.back());
        while (precursor.is_valid())
        {
            ns.push_back(mcMesh.from_vertex_handle(precursor));
            precursor = n2precursors.at(ns.back());
        }

        assert(seeds.count(ns.back()));

        for (VH nChain : ns)
            n2seed[nChain] = ns.back();
    }
    map<EH, pairTT<VH>> a2seeds;
    for (EH a : mcMesh.edges())
    {
        auto ns = mcMesh.edge_vertices(a);
        a2seeds[a] = {n2seed[ns[0]], n2seed[ns[1]]};
    }

    map<EH, int> a2component;
    int currentComponentIdx = 0;
    for (EH aStart : mcMesh.edges())
    {
        if (a2component.count(aStart))
            continue;
        HEH haStart = mcMesh.halfedge_handle(aStart, 0);
        auto origins = a2seeds.at(aStart);
        if (origins.first.idx() == origins.second.idx())
            continue;
        if (origins.first.idx() > origins.second.idx())
        {
            std::swap(origins.first, origins.second);
            haStart = mcMesh.opposite_halfedge_handle(haStart);
        }
        auto& currentInterface = seedPairs2interfaces[origins].emplace_back();

        list<HEH> haQ = {{haStart}};
        currentInterface.insert(aStart);
        a2component[aStart] = currentComponentIdx;

        while (!haQ.empty())
        {
            HEH haCurrent = haQ.front();
            haQ.pop_front();

            for (HFH hp : mcMesh.halfedge_halffaces(haCurrent))
            {
                vector<HEH> haCycle;
                HEH haNext = haCurrent;
                do
                {
                    haCycle.push_back(haNext);
                    haNext = mcMesh.next_halfedge_in_halfface(haNext, hp);
                } while (haNext != haCurrent);

                bool foundNext = false;
                for (int i = 1; i < (int)haCycle.size(); i++)
                {
                    haNext = haCycle[i];
                    auto originsNext = a2seeds.at(mcMesh.edge_handle(haNext));
                    if (originsNext.first == originsNext.second)
                        continue;
                    if (haNext.idx() % 2 != 0)
                        std::swap(originsNext.first, originsNext.second);
                    if (originsNext.second == origins.first)
                    {
                        foundNext = true;
                        if (!a2component.count(mcMesh.edge_handle(haNext)))
                        {
                            haQ.push_back(mcMesh.opposite_halfedge_handle(haNext));
                            currentInterface.insert(mcMesh.edge_handle(haNext));
                            a2component[mcMesh.edge_handle(haNext)] = currentComponentIdx;
                        }
                    }
                    break;
                }
                if (!foundNext)
                {
                    for (int i = haCycle.size() - 1; i > 0; i--)
                    {
                        haNext = mcMesh.opposite_halfedge_handle(haCycle[i]);
                        auto originsNext = a2seeds.at(mcMesh.edge_handle(haNext));
                        if (originsNext.first == originsNext.second)
                            continue;
                        if (haNext.idx() % 2 != 0)
                            std::swap(originsNext.first, originsNext.second);
                        if (originsNext.second == origins.second)
                        {
                            foundNext = true;
                            if (!a2component.count(mcMesh.edge_handle(haNext)))
                            {
                                haQ.push_back(haNext);
                                currentInterface.insert(mcMesh.edge_handle(haNext));
                                a2component[mcMesh.edge_handle(haNext)] = currentComponentIdx;
                            }
                        }
                        break;
                    }
                }
            }
        }
        currentComponentIdx++;
    }

    return SUCCESS;
}

ObjectiveBuilder::RetCode
ObjectiveBuilder::getInterfaceTraversingPathsMC(const map<VH, HEH>& n2precursors,
                                                const map<pairTT<VH>, vector<set<EH>>>& seedPairs2interfaces,
                                                vector<vector<HEH>>& paths,
                                                vector<double>& pathWeights) const
{
    auto& mcMesh = mcMeshProps().mesh();
    paths.clear();
    pathWeights.clear();

    auto computeBlockCenter = [&mcMesh, this](const CH& b) -> Vec3d
    {
        VH n1 = mcMeshProps().get<BLOCK_CORNER_NODES>(b).at(UVWDir::NEG_U | UVWDir::NEG_V | UVWDir::NEG_W);
        VH n2 = mcMeshProps().get<BLOCK_CORNER_NODES>(b).at(UVWDir::POS_U | UVWDir::POS_V | UVWDir::POS_W);

        return 0.5 * (Vec3Q2d(nodeUVWinBlock(n1, b)) + Vec3Q2d(nodeUVWinBlock(n2, b)));
    };

    auto computePatchCenter = [&mcMesh, this](const HFH& hp) -> Vec3d
    {
        auto hasByDir = halfpatchHalfarcsByDir(hp);
        VH n1 = mcMesh.from_vertex_handle(hasByDir.begin()->second.front());
        VH n2 = mcMesh.from_vertex_handle(hasByDir.at(-hasByDir.begin()->first).front());

        CH b = mcMesh.incident_cell(hp);
        return 0.5 * (Vec3Q2d(nodeUVWinBlock(n1, b)) + Vec3Q2d(nodeUVWinBlock(n2, b)));
    };

    auto computeArcCenter = [&mcMesh, this](const EH& a, const CH& b) -> Vec3d
    {
        VH n1 = mcMesh.from_vertex_handle(mcMesh.halfedge_handle(a, 0));
        VH n2 = mcMesh.from_vertex_handle(mcMesh.halfedge_handle(a, 1));
        return 0.5 * (Vec3Q2d(nodeUVWinBlock(n1, b)) + Vec3Q2d(nodeUVWinBlock(n2, b)));
    };

    map<CH, Vec3d> b2center;
    map<HFH, Vec3d> hp2center;
    map<pair<CH, EH>, Vec3d> a2center;
    for (auto& kv : seedPairs2interfaces)
        for (auto& interface : kv.second)
        {
            double totalWeight = 0.0;
            for (EH a : interface)
            {
                for (CH b : mcMesh.edge_cells(a))
                {
                    auto bCenterIt = b2center.find(b);
                    if (bCenterIt == b2center.end())
                        bCenterIt = b2center.insert({b, computeBlockCenter(b)}).first;
                    HEH ha0 = mcMesh.halfedge_handle(a, 0);
                    HFH hp1;
                    HFH hp2;
                    for (HFH hp : mcMesh.halfedge_halffaces(ha0))
                        if (mcMesh.incident_cell(hp) == b)
                        {
                            hp1 = hp;
                            hp2 = safeAdjacentHalffaceInBlock(hp, ha0);
                            break;
                        }
                    auto hpCenterIt1 = hp2center.find(hp1);
                    if (hpCenterIt1 == hp2center.end())
                    {
                        Vec3d center = computePatchCenter(hp1);
                        hpCenterIt1 = hp2center.insert({hp1, center}).first;
                        hp2center.insert(
                            {mcMesh.opposite_halfface_handle(hp1),
                             Vec3Q2d(mcMeshProps().hpTransition<PATCH_TRANSITION>(hp1).apply(Vec3Q(center)))});
                    }
                    auto hpCenterIt2 = hp2center.find(hp2);
                    if (hpCenterIt2 == hp2center.end())
                    {
                        Vec3d center = computePatchCenter(hp2);
                        hpCenterIt2 = hp2center.insert({hp2, center}).first;
                        hp2center.insert(
                            {mcMesh.opposite_halfface_handle(hp2),
                             Vec3Q2d(mcMeshProps().hpTransition<PATCH_TRANSITION>(hp2).apply(Vec3Q(center)))});
                    }

                    auto aCenterIt = a2center.find({b, a});
                    if (aCenterIt == a2center.end())
                        aCenterIt = a2center.insert({std::make_pair(b, a), computeArcCenter(a, b)}).first;

                    // Compute area: 0.5 * ((p1-a) x (b-a)).norm() + 0.5 * ((p2 - a) x (b-a)).norm()
                    Vec3d ba = bCenterIt->second - aCenterIt->second;
                    Vec3d hp1a = hpCenterIt1->second - aCenterIt->second;
                    Vec3d hp2a = hpCenterIt2->second - aCenterIt->second;
                    double area = 0.5 * ((hp1a % ba).norm() + (hp2a % ba).norm());
                    totalWeight += area;
                }
            }

            list<HEH> pathMin;
            double lenMin = DBL_MAX;
            for (EH aStart : interface)
            {
                auto nsStart = mcMesh.edge_vertices(aStart);

                list<HEH> path({mcMesh.halfedge_handle(aStart, 0)});

                HEH precursorFront = n2precursors.at(nsStart[0]);
                while (precursorFront.is_valid())
                {
                    path.push_front(precursorFront);
                    precursorFront = n2precursors.at(mcMesh.from_vertex_handle(precursorFront));
                }
                HEH precursorBack = n2precursors.at(nsStart[1]);
                while (precursorBack.is_valid())
                {
                    path.push_back(mcMesh.opposite_halfedge_handle(precursorBack));
                    precursorBack = n2precursors.at(mcMesh.from_vertex_handle(precursorBack));
                }
                double len = 0.0;
                for (HEH ha : path)
                    len += mcMeshProps().get<ARC_DBL_LENGTH>(mcMesh.edge_handle(ha));
                if (len < lenMin)
                {
                    lenMin = len;
                    pathMin = path;
                }
            }

            paths.emplace_back(pathMin.begin(), pathMin.end());
            pathWeights.push_back(totalWeight / lenMin);
        }

    // Normalize weights
    double sumOfWeights = 0.0;
    for (double weight : pathWeights)
        sumOfWeights += weight;
    for (double& weight : pathWeights)
        weight /= sumOfWeights;

    return SUCCESS;
}

void ObjectiveBuilder::debugViewVoronoiCells(const set<VH>& seeds, const map<VH, HEH>& n2precursors) const
{
#ifdef MC3D_WITH_VIEWER
    auto& mcMesh = mcMeshProps().mesh();
    map<VH, VH> n2seed;
    for (VH n : mcMesh.vertices())
    {
        if (n2seed.count(n))
            continue;
        vector<VH> ns = {{n}};
        HEH precursor = n2precursors.at(ns.back());
        while (precursor.is_valid())
        {
            ns.push_back(mcMesh.from_vertex_handle(precursor));
            precursor = n2precursors.at(ns.back());
        }

        assert(seeds.count(ns.back()));

        for (VH nChain : ns)
            n2seed[nChain] = ns.back();
    }
    map<EH, pairTT<VH>> a2seeds;
    for (EH a : mcMesh.edges())
    {
        auto ns = mcMesh.edge_vertices(a);
        a2seeds[a] = {n2seed[ns[0]], n2seed[ns[1]]};
    }
    auto& tetMesh = meshProps().mesh();
    set<HFH> patchHfs;
    for (FH p : mcMesh.faces())
        for (HFH hf : mcMeshProps().ref<PATCH_MESH_HALFFACES>(p))
            for (HFH hf2 : tetMesh.face_halffaces(tetMesh.face_handle(hf)))
                patchHfs.insert(hf2);
    set<EH> homoEs, heteroEs;
    for (EH a : mcMesh.edges())
    {
        auto& aSeeds = a2seeds.at(a);
        if (aSeeds.first != aSeeds.second)
            for (HEH he : mcMeshProps().ref<ARC_MESH_HALFEDGES>(a))
                heteroEs.insert(tetMesh.edge_handle(he));
        else
            for (HEH he : mcMeshProps().ref<ARC_MESH_HALFEDGES>(a))
                homoEs.insert(tetMesh.edge_handle(he));
    }
    set<VH> nodeVs, seedVs;
    for (VH n : mcMesh.vertices())
    {
        if (seeds.count(n))
            seedVs.insert(mcMeshProps().ref<NODE_MESH_VERTEX>(n));
        else
            nodeVs.insert(mcMeshProps().ref<NODE_MESH_VERTEX>(n));
    }
    debugView({}, patchHfs, {}, {}, heteroEs, {}, homoEs, {}, seedVs, nodeVs);
#endif
}

void ObjectiveBuilder::debugViewPaths(const map<pairTT<VH>, vector<set<EH>>>& seedPairs2interfaces,
                                      const vector<vector<HEH>>& paths) const
{
#ifdef MC3D_WITH_VIEWER
    auto& mcMesh = mcMeshProps().mesh();
    auto& tetMesh = meshProps().mesh();
    set<HFH> patchHfs;
    for (FH p : mcMesh.faces())
        for (HFH hf : mcMeshProps().ref<PATCH_MESH_HALFFACES>(p))
            for (HFH hf2 : tetMesh.face_halffaces(tetMesh.face_handle(hf)))
                patchHfs.insert(hf2);
    set<EH> pathEs;
    for (EH a : mcMesh.edges())
    {
        for (auto& path : paths)
        {
            bool found = false;
            for (auto ha : path)
                if (mcMesh.edge_handle(ha) == a)
                {
                    found = true;
                    break;
                }
            if (found)
            {
                for (HEH he : mcMeshProps().ref<ARC_MESH_HALFEDGES>(a))
                    pathEs.insert(tetMesh.edge_handle(he));
                break;
            }
        }
    }
    set<VH> nodeVs, seedVs;
    for (VH n : mcMesh.vertices())
    {
        bool found = false;
        for (auto& kv : seedPairs2interfaces)
            if (kv.first.first == n || kv.first.second == n)
                found = true;

        if (found)
            seedVs.insert(mcMeshProps().ref<NODE_MESH_VERTEX>(n));
        else
            nodeVs.insert(mcMeshProps().ref<NODE_MESH_VERTEX>(n));
    }
    debugView({}, patchHfs, {}, {}, pathEs, {}, {}, {}, seedVs, nodeVs);
#endif
}

} // namespace qgp3d
