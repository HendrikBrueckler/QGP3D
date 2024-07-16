#ifndef QGP3D_WITH_GUROBI

#ifndef QGP3D_CLPLPSOLVER_HPP
#define QGP3D_CLPLPSOLVER_HPP

#include "QGP3D/ISP/BaseLPSolver.hpp"
#include <MC3D/Mesh/MCMeshNavigator.hpp>

#include "ClpSimplex.hpp"

namespace qgp3d
{
namespace impl
{
using namespace mc3d;

/**
 * @brief Class managing LP solving via CLP for integer-sheet-pump quantization
 */
class ClpLPSolver : public BaseLPSolver
{
  public:
    virtual ~ClpLPSolver()
    {
    }

    /**
     * @brief Create a new CLP LP solver instance
     *
     * @param meshProps IN: mesh being quantized
     * @param scaling IN: scaling factor for target lengths
     * @param decomp IN: decomposition into subproblems
     */
    ClpLPSolver(const TetMeshProps& meshProps, double scaling, const ISPQuantizer::Decomposition& decomp);

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

    vector<std::unique_ptr<ClpSimplex>> _subproblem2model; // Different LP instance per subproblem
    vector<map<int, pairTT<int>>> _subproblem2bundle2vars; // Different set of variables per subproblem

    vector<int> _subproblem2numberOfRowsBase; // To remove temporary constraints, cut constraint matrix down to this
};

} // namespace impl
} // namespace qgp3d

#endif
#endif
