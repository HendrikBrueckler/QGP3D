#ifndef QGP3D_ISPQUANTIZER_HPP
#define QGP3D_ISPQUANTIZER_HPP

#include <MC3D/Mesh/MCMeshManipulator.hpp>

#include "QGP3D/SeparationChecker.hpp"

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
     * @param sep IN: separationchecker to use for lazy separation. May come with predetermined separation constraints
     */
    ISPQuantizer(TetMeshProps& meshProps, SeparationChecker& sep);

    /**
     * @brief Compute quantization
     *
     * @param scaling IN: scale target lengths by this factor for quantization
     * @param varLowerBound IN: lower bound for arc lengths
     * @return RetCode SUCCESS or error code
     */
    RetCode quantize(double scaling = 1.0, double varLowerBound = 0.0);

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
     * @return double objective value after descent
     */
    double greedyDescent(BaseLPSolver& sheetFinder, double scaling, bool useGlobalProblem);

    /**
     * @brief Use sheet inflation/deflation operators obtained from LP solves to recover quantization feasibility
     *        (after additional constraints were added)
     *
     * @param sheetFinder IN: instance to handle LP solving (polymorphic interface)
     * @param scaling IN: scale target lengths by this factor for quantization
     * @return double objective value after descent
     */
    double makeFeasible(BaseLPSolver& sheetFinder, double scaling);

    /**
     * @brief Get all simple constraints (non-separation) violated by the current quantization
     *
     * @param varLowerBound IN: lower bound of variables
     * @param criticalLinks IN: critical links previously determined
     * @return vector<vector<pair<int, EH>>> violated constraints (encoded non-intuitively)
     */
    vector<vector<pair<int, EH>>> violatedSimpleConstraints(double varLowerBound, vector<CriticalLink>& criticalLinks);

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
     * @return double objective value
     */
    double obj(double scaling) const;

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

    SeparationChecker& _sep; // Separation checker given from outside
    Decomposition _decomp;   // Decomposition of the MC domain into quantization subproblems
};

} // namespace qgp3d
#endif
