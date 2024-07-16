#ifndef QGP3D_BASELPSOLVER_HPP
#define QGP3D_BASELPSOLVER_HPP

#include <MC3D/Mesh/MCMeshNavigator.hpp>

#include "QGP3D/ISP/ISPQuantizer.hpp"

namespace qgp3d
{
using namespace mc3d;

/**
 * @brief // Class managing LP solvers used in integer-sheet-pump quantization
 */
class BaseLPSolver : public virtual MCMeshNavigator
{
  public:
    virtual ~BaseLPSolver()
    {
    }

    /**
     * @brief Creates the base setup for LP solvers
     *
     * @param meshProps IN: mesh being quantized
     * @param scaling IN: scaling factor for target lengths
     * @param decomp IN: decomposition into subproblems
     */
    BaseLPSolver(const TetMeshProps& meshProps, double scaling, const ISPQuantizer::Decomposition& decomp)
        : TetMeshNavigator(meshProps), MCMeshNavigator(meshProps), _scaling(scaling), _decomp(decomp)
    {
    }

    /**
     * @brief Setup the static parts of the LP solver, excluding dynamic parts. Should be called before any method
     *        below.
     */
    virtual void setupLPBase() = 0;

    /**
     * @brief Get mixed-sign sheet passing through variable bundle \p i , respecting
     *        separation constraints. Scaling factor is already chosen optimally.
     *        If no solution could be found, returned sheet is empty.
     *
     * @param bundleToChange IN: variable (bundle) sheet should pass through
     * @param add IN: whether variable \p i should be inflated (else deflated). If \p fixing , means inflation-only
     * @param fixing IN: whether inequality constraints are violated by current solution and need fixing
     * @param useGlobalProblem IN: whether to compute solution on global problem instead of subproblem
     * @return map<int, double> mapping from bundles to integer deltas
     */
    map<int, double> integerSheet(int bundleToChange, bool add, bool fixing, bool useGlobalProblem);

    /**
     * @brief Set dynamic inequality constraints (encoded non-intuitively)
     *
     * @param newConstraints IN: new dynamic constraints
     */
    void setDynamicConstraints(const vector<vector<pair<int, EH>>>& newConstraints)
    {
        _dynamicConstraints = newConstraints;
    }

    /**
     * @brief Get dynamic inequality constraints (encoded non-intuitively)
     *
     * @return const vector<vector<pair<int, EH>>>& dynamic constraints previously added
     */
    const vector<vector<pair<int, EH>>>& dynamicConstraints() const
    {
        return _dynamicConstraints;
    }

  protected:
    /**
     * @brief Optimal factor by which to scale a given sheet to minimize deviation objective
     *
     * @param sheet IN: integer sheet
     * @return int optimal integer scaling factor (within feasible bounds, 0 if always infeasible)
     */
    int optimalFactor(const map<int, double>& sheet) const;

    /**
     * @brief Unfavorability weights for LP variables
     *
     * @param bundle IN: bundle variable index
     * @param inflating IN: whether inflating or deflating weight should be returned
     * @return double unfavorability weight of \p bundle
     */
    double weight(int bundle, bool inflating);

    /// The following are called from integerSheet() and need to be implemented in child classes

    //@{
    /**
     * @brief Setup the objective function, that takes the current quantization into account.
     *        Called before each call to solve()
     *
     * @param subproblem IN: subproblem to setup
     */
    virtual void setupDynamicObjective(int subproblem) = 0;

    /**
     * @brief Setup temporary pump constraints that enforce inflation or deflation of given bundle
     *
     * @param subproblem IN: subproblem to setup
     * @param bundle IN: bundle to constrain
     * @param inflate IN: whether to force inflation (else deflation)
     */
    virtual void setupPumpConstraints(int subproblem, int bundle, bool inflate) = 0;

    /**
     * @brief Setup temporary constraints that disable deflation.
     *
     * @param subproblem IN: subproblem to setup
     */
    virtual void setupInflationOnlyConstraints(int subproblem) = 0;

    /**
     * @brief Setup temporary inequality constraints enforcing nonnegativity of variables/blocks and separation.
     *
     * @param subproblem IN: subproblem to setup
     */
    virtual bool setupSepNonnegConstraints(int subproblem) = 0;

    /**
     * @brief Solve the previously setup problem
     *
     * @param subproblem IN: subproblem to solve
     * @return true if feasible and successful
     * @return false else
     */
    virtual bool solve(int subproblem) = 0;

    /**
     * @brief Last solution value for given bundle obtained via solve()
     *
     * @param subproblem IN: subproblem solved before
     * @param bundle IN: bundle to query
     * @return double LP solution for \p bundle
     */
    virtual double solution(int subproblem, int bundle) = 0;

    /**
     * @brief Add temporary constraint that, based on the previous solution, fixes \p bundle to the next integer value
     *        further away from 0.
     *
     * @param subproblem IN: subproblem solved before
     * @param bundle IN: bundle to constrain
     */
    virtual void reconstrainToNextInt(int subproblem, int bundle) = 0;

    /**
     * @brief Remove all temporary constraints from a problem
     *
     * @param subproblem IN: subproblem to remove constraints from
     */
    virtual void removeTemporaryConstraints(int subproblem) = 0;
    //@}

    double _scaling;                                   // scaling factor for target lengths
    const ISPQuantizer::Decomposition& _decomp;        // decomposition into subproblems passed from outside
    vector<vector<pair<int, EH>>> _dynamicConstraints; // inequality constraints dynamically and incrementally added
};

} // namespace qgp3d

#endif
