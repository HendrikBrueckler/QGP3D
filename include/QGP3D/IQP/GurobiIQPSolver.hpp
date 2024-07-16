#ifdef QGP3D_WITH_GUROBI
#ifndef QGP3D_WITHOUT_IQP

#ifndef QGP3D_GUROBIIQPSOLVER_HPP
#define QGP3D_GUROBIIQPSOLVER_HPP

#include <gurobi_c++.h>

#include "QGP3D/IQP/BaseIQPSolver.hpp"
#include <MC3D/Mesh/MCMeshManipulator.hpp>

namespace qgp3d
{
namespace impl
{

/**
 * @brief IQP quantization solving based on GUROBI
 */
class GurobiIQPSolver : public BaseIQPSolver, public virtual MCMeshManipulator
{
  public:
    /**
     * @brief Create new IQP solver instance based on Bonmin
     *
     * @param meshProps IN: mesh whose MC to quantize
     * @param scaling IN: scale target lengths by this factor for quantization
     * @param varLowerBound IN: lower bound for arc lengths
     * @param maxSeconds IN: time limit for solver in seconds
     * @param individualArcFactor IN: objective = this * <arc-length-deviation> + (1-this) * <block-length-deviation>
     */
    GurobiIQPSolver(TetMeshProps& meshProps,
                    double scaling,
                    double varLowerBound,
                    double maxSeconds = 180,
                    double individualArcFactor = 1.0);

    /// @brief See \ref BaseIQPSolver for usage

    ///@{
    virtual void setupDefaultOptions();
    virtual void setupVariables();
    virtual void setupObjective();
    virtual void setupConstraints();
    virtual void finalizeSetup();
    virtual void enableQuickSolve();
    virtual void enableExactSolve();
    virtual BaseIQPSolver::RetCode solve();
    virtual double objectiveValue() const;
    virtual void useCurrentAsWarmStart();
    virtual void addConstraints(const vector<vector<pair<int, EH>>>& nonZeroSum);
    ///@}

    virtual ~GurobiIQPSolver()
    {
    }

  private:
    GRBEnv _env;                     // Gurobi environment
    GRBModel _model;                 // Gurobi problem instance
    map<EH, GRBVar> _arc2var;        // Matching of arcs to IQP variables
    int _nSeparationConstraints = 0; // Total number of separation constraints added to the IQP
};

} // namespace impl
} // namespace qgp3d

#endif
#endif
#endif
