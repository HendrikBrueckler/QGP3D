#ifndef QGP3D_WITHOUT_IQP

#ifndef QGP3D_BASEIQPSOLVER_HPP
#define QGP3D_BASEIQPSOLVER_HPP

#include <MC3D/Types.hpp>

namespace qgp3d
{
using namespace mc3d;

/**
 * @brief Base interface for IQP quantization solvers
 */
class BaseIQPSolver
{
  public:
    enum RetCode
    {
        SUCCESS = 0,
        SOLVER_ERROR = 20, // Unspecific solver error. Returned when no solution found within time limit.
    };

    /**
     * @brief Initialize common parameters in base class
     *
     * @param scaling IN: scale target lengths by this factor for quantization
     * @param varLowerBound IN: lower bound for arc lengths
     * @param maxSeconds IN: time limit for solver in seconds
     * @param individualArcFactor IN: objective = this * <arc-length-deviation> + (1-this) * <block-length-deviation>
     */
    BaseIQPSolver(double scaling = 1.0,
                  double varLowerBound = 0.0,
                  double maxSeconds = 300,
                  double individualArcFactor = 1.0)
        : _scaling(scaling), _varLowerBound(varLowerBound), _maxSeconds(maxSeconds),
          _individualArcFactor(individualArcFactor)
    {
    }

    virtual ~BaseIQPSolver() = default;

    /**
     * @brief Configure the solver here
     */
    virtual void setupDefaultOptions() = 0;

    /**
     * @brief Create variables, bounds and initial values
     */
    virtual void setupVariables() = 0;

    /**
     * @brief Configure objective function
     */
    virtual void setupObjective() = 0;

    /**
     * @brief Configure constraints
     */
    virtual void setupConstraints() = 0;

    /**
     * @brief Apply setup and prepare for solve
     */
    virtual void finalizeSetup() = 0;

    /**
     * @brief Switch to some form of quick feasible but suboptimal solving
     */
    virtual void enableQuickSolve() = 0;

    /**
     * @brief Switch to exact solving with optimality gap <1e-4
     */
    virtual void enableExactSolve() = 0;

    /**
     * @brief Solve the problem as configured
     *
     * @return RetCode SUCCESS or SOLVER_ERROR
     */
    virtual RetCode solve() = 0;

    /**
     * @brief Objective value of last solve
     *
     * @return double objective value
     */
    virtual double objectiveValue() const = 0;

    /**
     * @brief Use the previous solution as warm start for next solve
     */
    virtual void useCurrentAsWarmStart() = 0;

    /**
     * @brief Add additional constraints (after setup was finalized)
     * @param nonZeroSum IN: constraints as encoded in \ref SeparationChecker
     */
    virtual void addConstraints(const vector<vector<pair<int, EH>>>& nonZeroSum) = 0;

  protected:
    double _scaling;             // scale target lengths by this factor for quantization
    double _varLowerBound;       // lower bound for arc lengths
    double _maxSeconds;          // time limit for solver in seconds
    double _individualArcFactor; // objective = this * <arc-length-deviation> + (1-this) * <block-length-deviation>
};

} // namespace qgp3d
#endif

#endif
