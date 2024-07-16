#ifdef QGP3D_WITH_GUROBI

#ifndef QGP3D_GUROBILPSOLVER_HPP
#define QGP3D_GUROBILPSOLVER_HPP

#include <MC3D/Mesh/MCMeshNavigator.hpp>

#include <gurobi_c++.h>

#include "QGP3D/ISP/BaseLPSolver.hpp"

namespace qgp3d
{
namespace impl
{
using namespace mc3d;

/**
 * @brief Class managing LP solving via GUROBI for integer-sheet-pump quantization
 */
class GurobiLPSolver : public BaseLPSolver
{
  public:
    virtual ~GurobiLPSolver()
    {
    }

    /**
     * @brief Create a new CLP LP solver instance
     *
     * @param meshProps IN: mesh being quantized
     * @param scaling IN: scaling factor for target lengths
     * @param decomp IN: decomposition into subproblems
     */
    GurobiLPSolver(const TetMeshProps& meshProps, double scaling, const ISPQuantizer::Decomposition& decomp);

    /// @brief See \ref BaseIQPSolver for usage

    ///@{
    virtual void setupLPBase();

  protected:
    virtual void setupDynamicObjective(int subproblem);
    virtual void setupPumpConstraints(int subproblem, int bundle, bool inflate);
    virtual void setupInflationOnlyConstraints(int subproblem);
    virtual bool setupSepNonnegConstraints(int subproblem);
    virtual bool solve(int subproblem);
    virtual double solution(int subproblem, int bundle);
    virtual void reconstrainToNextInt(int subproblem, int bundle);
    virtual void removeTemporaryConstraints(int subproblem);
    ///@}

    GRBEnv _lpenv;                                            // Gurobi environment
    vector<std::unique_ptr<GRBModel>> _subproblem2model;      // Different LP instance per subproblem
    vector<map<int, pairTT<GRBVar>>> _subproblem2bundle2vars; // Different set of variables per subproblem

    vector<GRBConstr> _temporaryConstraints; // Temporary constraint storage, to remove them later
};
} // namespace impl

} // namespace qgp3d

#endif
#endif
