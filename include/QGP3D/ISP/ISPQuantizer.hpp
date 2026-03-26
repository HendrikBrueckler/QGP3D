#ifndef QGP3D_ISPQUANTIZER_HPP
#define QGP3D_ISPQUANTIZER_HPP

#include <MC3D/Mesh/MCMeshManipulator.hpp>

#include "QGP3D/ObjectiveFunction.hpp"
#include "QGP3D/StructurePreserver.hpp"

namespace qgp3d
{
using namespace mc3d;
class BaseLPSolver;

/**
 * @brief Class to execute quantization using the greedy [I]nteger[S]heet[P]ump algorithm
 */
class ISPQuantizer : public virtual MCMeshManipulator
{
  public:
    enum RetCode
    {
        SUCCESS = 0,
        SOLVER_ERROR = 20,
    };

    /**
     * @brief Struct to store information about (mostly) independent subproblems
     */
    struct Decomposition
    {
        map<EH, int> a2idx; // Arc IDs are not contiguous after MC modification, this maps arcs to contiguous indices
        vector<EH> idx2a;   // Arc IDs are not contiguous after MC modification, this maps indices back to arcs
        map<HFH, map<UVWDir, vector<HEH>>> hp2hasByDir;   // To speedup halfpatchHalfarcsByDir query
        map<EH, int> arc2bundle;                          // Arcs grouped into bundles, constrained to equal length
        vector<set<EH>> bundle2arcs;                      // Inverse map from bundle to arcs
        vector<int> bundle2subproblem;                    // Assignment of bundles to subproblems
        vector<set<int>> subproblem2bundles;              // Map of subproblems to contained bundles
        vector<set<pair<FH, UVWDir>>> subproblem2patches; // Map of subproblems to contained patches
    };

    /**
     * @brief Create an instance handling quantization via integer sheet pump
     *
     * @param meshProps IN: mesh whose MC to quantize
     * @param sep IN: StructurePreserver to use for lazy separation. May come with predetermined separation constraints
     * @param obj IN: objective function
     */
    ISPQuantizer(TetMeshProps& meshProps, StructurePreserver& sep, ObjectiveFunction& obj);

    /**
     * @brief Compute quantization
     *
     * @param varLowerBound IN: lower bound for arc lengths
     * @return RetCode SUCCESS or error code
     */
    RetCode quantize(double varLowerBound = 0.0);

  protected:
    /**
     * @brief Decompose MC domain into problems (mostly) independent, as specified above, stored internally
     */
    void decomposeIntoSubproblems();

    /**
     * @brief Perform a greedy descent using sheet inflation/deflation operators obtained from LP solves
     *
     * @param sheetFinder IN: instance to handle LP solving (polymorphic interface)
     * @param scaling IN: scale target lengths by this factor for quantization
     * @param useGlobalProblem IN: whether to compute solution on global problem instead of subproblem
     * @param individualArcFactor IN: importance of individual arc lengths vs critical path lengths
     * @return double objective value after descent
     */
    double greedyDescent(BaseLPSolver& sheetFinder, bool useGlobalProblem, bool thoroughSolve = false);

    /**
     * @brief Use sheet inflation/deflation operators obtained from LP solves to recover quantization feasibility
     *        (after additional constraints were added)
     *
     * @param sheetFinder IN: instance to handle LP solving (polymorphic interface)
     * @param scaling IN: scale target lengths by this factor for quantization
     * @param individualArcFactor IN: importance of individual arc lengths vs critical path lengths
     * @return double objective value after descent
     */
    double makeFeasible(BaseLPSolver& sheetFinder);

    /**
     * @brief Number of constraints in \p constraints violated by the current quantization
     *
     * @param constraints IN: constraints for which to query violation
     * @return int number of violated constraints
     */
    int numViolatedConstraints(const vector<vector<pair<int, EH>>>& constraints) const;

    /**
     * @brief Objective value of current quantization
     *
     * @param scaling IN: scaling factor for target lengths
     * @param individualArcFactor IN: importance of individual arc lengths vs critical path lengths
     * @return double objective value
     */
    double obj() const;

    /**
     * @brief Delta of objective value of current quantization when applying sheet
     *
     * @param scaling IN: scaling factor for target lengths
     * @param individualArcFactor IN: importance of individual arc lengths vs critical path lengths
     * @param sheet IN: sheet to apply
     * @return double objective value
     */
    double objDelta(const map<int, double>& sheet);

    /**
     * @brief Set of other bundles affected directly/semi-directly by quantization change via sheet.
     *
     * @param sheet IN: sheet previously applied
     * @param bundle2sheet IN: map of bundles to most recently determined sheets in/deflating them
     * @param constraints IN: set of constraints connecting different bundles
     * @return set<int> bundles affected by sheet
     */
    set<int> affectedBundles(const map<int, double>& sheet,
                             const vector<map<int, double>>& bundle2sheet,
                             const vector<vector<pair<int, EH>>>& constraints) const;

    /**
     * @brief Set of other bundles affected directly/semi-directly by quantization change via sheet.
     *
     * @param sheet IN: sheet previously applied
     * @param biBundle2sheet IN: map of bi-bundles to most recently determined sheets in/deflating them
     * @param constraints IN: set of constraints connecting different bundles
     * @return set<int> bundles affected by sheet
     */
    set<int> affectedBundles(const map<int, double>& sheet,
                             const map<pairTT<int>, map<int, double>>& biBundle2sheet,
                             const vector<vector<pair<int, EH>>>& constraints) const;

    /**
     * @brief Cache gradient and hessian of current quantization
     *
     * @param obj IN: Objective function
     */
    void cacheGradHess(ObjectiveFunction& obj);

    /**
     * @brief Apply sheet to underlying quantization
     *
     * @param sheet IN: integer sheet
     * @param add IN: whether to add or subtract
     */
    void applyChange(const map<int, double>& sheet, bool add = true);

    StructurePreserver& _sep;                       // Structure preserver (managed externally)
    ObjectiveFunction& _obj;                        // Objective function (managed externally)
    Decomposition _decomp;                          // Decomposition of the MC domain into quantization subproblems
    Eigen::VectorXd _currentGrad;                   // Objective gradient at current quantization
    Eigen::SparseMatrix<double> _currentHess;       // Objective hessian at current quantization
    Eigen::VectorXd _currentContinuousQuadraticOpt; // Optimum of quadratic objective approximation (trivial if
                                                    // target arc lengths given)
};

} // namespace qgp3d
#endif
