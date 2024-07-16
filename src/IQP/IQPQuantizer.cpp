#ifndef QGP3D_WITHOUT_IQP

#include "QGP3D/IQP/IQPQuantizer.hpp"

#define INDIVIDUAL_ARC_FACTOR 1.0

#ifdef QGP3D_WITH_GUROBI
#include "QGP3D/IQP/GurobiIQPSolver.hpp"
#include <gurobi_c++.h>

#else // QGP3D with CLP or BONMIN

#include "QGP3D/IQP/BonminIQPSolver.hpp"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseQR>

#include "BonBonminSetup.hpp"
#include "BonCbc.hpp"
#include "BonIpoptSolver.hpp"
#include "BonOsiTMINLPInterface.hpp"
#include "BonTMINLP.hpp"
#include "BonEcpCuts.hpp"
#include "BonOACutGenerator2.hpp"
#include "BonOaNlpOptim.hpp"
#include "BonBonminSetup.hpp"
#endif

namespace qgp3d{

IQPQuantizer::IQPQuantizer(TetMeshProps& meshProps, SeparationChecker& sep)
    : TetMeshNavigator(meshProps), TetMeshManipulator(meshProps), MCMeshNavigator(meshProps),
      MCMeshManipulator(meshProps), _sep(sep)
{
}

IQPQuantizer::RetCode IQPQuantizer::quantize(double scaling, double varLowerBound, int maxSecondsIQP)
{
    const MCMesh& mc = mcMeshProps().mesh();

    map<EH, int> previousSolution;

    bool wasAllocated = mcMeshProps().isAllocated<ARC_DBL_LENGTH>();
    if (!wasAllocated)
    {
        mcMeshProps().allocate<ARC_DBL_LENGTH>(0.0);
        for (EH arc : mc.edges())
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
    else
        for (EH a: mcMeshProps().mesh().edges())
            previousSolution[a] = mcMeshProps().get<ARC_INT_LENGTH>(a);

    vector<CriticalLink> criticalLinks;
    map<EH, int> a2criticalLinkIdx;
    map<VH, vector<int>> n2criticalLinksOut;
    map<VH, vector<int>> n2criticalLinksIn;
    getCriticalLinks(criticalLinks, a2criticalLinkIdx, n2criticalLinksOut, n2criticalLinksIn, true);

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



#ifdef QGP3D_WITH_GUROBI
    impl::GurobiIQPSolver iqp(meshProps(), scaling, std::max(-GRB_INFINITY, varLowerBound), maxSecondsIQP, INDIVIDUAL_ARC_FACTOR);
#else
    impl::BonminIQPSolver iqp(meshProps(), scaling, varLowerBound, maxSecondsIQP, INDIVIDUAL_ARC_FACTOR);
#endif
    iqp.setupDefaultOptions();
    iqp.setupVariables();
    iqp.setupObjective();
    iqp.setupConstraints();
    iqp.finalizeSetup();

    int iter = 0;
    double obj = 0.0;
    // Iteratively enforce violated separation constraints
    bool validSolution = false;

    bool nextSolveExact = true;
    if (varLowerBound <= 0.0)
    {
        // If separation nontrivial and not predetermined, dont waste too much time on first few solves
        if (_sep.previousSeparationViolatingPaths().empty())
            nextSolveExact = false;
        else
            iqp.addConstraints(_sep.previousSeparationViolatingPaths());
    }

    while (!validSolution || nextSolveExact)
    {
        iter++;
        if (nextSolveExact)
            iqp.enableExactSolve();
        else
            iqp.enableQuickSolve();
        auto ret = iqp.solve();

        if (ret != BaseIQPSolver::SUCCESS)
        {
            LOG(WARNING) << "IQP solver could not find a valid solution in given time, reverting to previous solution (if existing)";
            for (auto& kv: previousSolution)
                mcMeshProps().set<ARC_INT_LENGTH>(kv.first, kv.second);
            return SOLVER_ERROR;
        }

        obj = iqp.objectiveValue();

        int minArcLength = 10000;
        for (EH a : mc.edges())
            minArcLength = std::min(minArcLength, mcMeshProps().get<ARC_INT_LENGTH>(a));
        DLOG(INFO) << "Min quantized arc length " << minArcLength;

        if (varLowerBound > 0.0)
        {
            validSolution = true;
            nextSolveExact = false;
        }
        else
        {
            vector<vector<pair<int, EH>>> forcedNonZeroSum;
            _sep.findSeparationViolatingPaths(
                criticalLinks, isCriticalArc, isCriticalNode, isCriticalPatch, forcedNonZeroSum);
            validSolution = forcedNonZeroSum.size() == 0;

            if (!validSolution)
            {
                iqp.addConstraints(forcedNonZeroSum);
                iqp.useCurrentAsWarmStart();
                nextSolveExact = false;
            }
            else
            {
                for (EH a: mcMeshProps().mesh().edges())
                    previousSolution[a] = mcMeshProps().get<ARC_INT_LENGTH>(a);
                nextSolveExact = !nextSolveExact;
            }
        }
    }

    int minArcLength = 10000;
    for (EH a : mc.edges())
        minArcLength = std::min(minArcLength, mcMeshProps().get<ARC_INT_LENGTH>(a));
    int nHexes = _sep.numHexesInQuantization();

    LOG(INFO) << "Quantized MC, nHexes: " << nHexes << ", iterations: " << iter
                << ", minArcLength: " << minArcLength << ", objective value: " << obj << std::endl;

    return SUCCESS;
}

} // namespace qgp3d

#endif
