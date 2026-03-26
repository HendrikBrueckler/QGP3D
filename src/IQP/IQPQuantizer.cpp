#ifndef QGP3D_WITHOUT_IQP

#include "QGP3D/IQP/IQPQuantizer.hpp"
#include "QGP3D/ObjectiveBuilder.hpp"
#include "QGP3D/ObjectiveFunction.hpp"

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
#include "BonEcpCuts.hpp"
#include "BonIpoptSolver.hpp"
#include "BonOACutGenerator2.hpp"
#include "BonOaNlpOptim.hpp"
#include "BonOsiTMINLPInterface.hpp"
#include "BonTMINLP.hpp"
#endif

namespace qgp3d
{

IQPQuantizer::IQPQuantizer(TetMeshProps& meshProps_, StructurePreserver& sep, QuadraticObjective& obj)
    : TetMeshNavigator(meshProps_), TetMeshManipulator(meshProps_), MCMeshNavigator(meshProps_),
      MCMeshManipulator(meshProps_), _sep(sep), _obj(obj)
{
}

IQPQuantizer::RetCode IQPQuantizer::quantize(double varLowerBound, int maxSecondsIQP)
{
    const MCMesh& mc = mcMeshProps().mesh();

    map<EH, int> previousValidSolution;
    if (!mcMeshProps().isAllocated<ARC_INT_LENGTH>())
        mcMeshProps().allocate<ARC_INT_LENGTH>(0);
    else
        for (EH a : mcMeshProps().mesh().edges())
            previousValidSolution[a] = mcMeshProps().get<ARC_INT_LENGTH>(a);

    map<EH, int> previousSolution;
    for (EH a : mcMeshProps().mesh().edges())
        previousSolution[a] = mcMeshProps().get<ARC_INT_LENGTH>(a);

#ifdef QGP3D_WITH_GUROBI
    impl::GurobiIQPSolver iqp(meshProps(), std::max(-GRB_INFINITY, varLowerBound), maxSecondsIQP);
#else
    impl::BonminIQPSolver iqp(meshProps(), varLowerBound, maxSecondsIQP);
#endif
    iqp.setupDefaultOptions();
    iqp.setupVariables();
    iqp.setupObjective(_obj);
    iqp.setupConstraints();
    iqp.finalizeSetup();

    int iter = 0;
    // Iteratively enforce violated separation constraints
    bool validSolution = false;

    bool nextSolveExact = true;
    if (varLowerBound <= 0.0)
    {
        // If separation nontrivial and not predetermined, dont waste too much time on first few solves
        if (_sep.previousStructuralConstraints().empty())
            nextSolveExact = false;
        else
            iqp.addConstraints(_sep.previousStructuralConstraints());
    }

    bool uncheckedChanges = true;
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
            LOG(WARNING) << "IQP solver could not find a valid solution in given time, reverting to previous solution "
                            "(if existing)";
            for (auto& kv : previousValidSolution)
                mcMeshProps().set<ARC_INT_LENGTH>(kv.first, kv.second);
            return SOLVER_ERROR;
        }

        for (auto& kv : previousSolution)
            if (mcMeshProps().get<ARC_INT_LENGTH>(kv.first) != kv.second)
                uncheckedChanges = true;
        for (EH a : mcMeshProps().mesh().edges())
            previousSolution[a] = mcMeshProps().get<ARC_INT_LENGTH>(a);

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
            if (uncheckedChanges)
            {
                _sep.violatedStructuralConstraints(forcedNonZeroSum);
                uncheckedChanges = false;
            }
            validSolution = forcedNonZeroSum.size() == 0;

            if (!validSolution)
            {
                iqp.addConstraints(forcedNonZeroSum);
                iqp.useCurrentAsWarmStart();
                nextSolveExact = false;
            }
            else
            {
                previousValidSolution = previousSolution;
                nextSolveExact = !nextSolveExact;
                if (nextSolveExact)
                    iqp.useCurrentAsWarmStart();
            }
        }
    }

    return SUCCESS;
}

} // namespace qgp3d

#endif
