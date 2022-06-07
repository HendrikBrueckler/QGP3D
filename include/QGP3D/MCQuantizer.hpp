#ifndef QGP3D_MCQUANTIZER_HPP
#define QGP3D_MCQUANTIZER_HPP

#include <MC3D/Mesh/MCMeshManipulator.hpp>

namespace qgp3d
{
using namespace mc3d;

/**
 * @brief Class to manage the quantization of the MC and its arc lengths.
 *        First step towards an IGM.
 *
 */
class MCQuantizer : public virtual MCMeshManipulator
{
  public:
    enum RetCode
    {
        SUCCESS = 0,
        SOLVER_ERROR = 20,
    };

    /**
     * @brief Create an instance that manages the quantization of the MC associated with \p meshProps
     *
     * @param meshProps IN/OUT: mesh whose MC should be quantized
     */
    MCQuantizer(TetMeshProps& meshProps);

    /**
     * @brief Actually assign quantized consistent arc lengths to each arc in the MC.
     *        \p scaling determines the objective function, which tries to minimize
     *        |origLength * scaling - quantizedLength|^2
     *
     * @param scaling IN: relative factor by which to approximately scale the original arc lengths
     * @param allowZero IN: whether to allow 0-length arcs
     * @return RetCode SUCCESS or SOLVER_ERROR
     */
    RetCode quantizeArcLengths(double scaling, bool allowZero = false, bool allowNegative = false);

    /**
     * @brief Count the number of hex cells in a mesh extractable from a quantized parametrization
     *        conforming to the current quantization
     *
     * @return int number of hex cells
     */
    int numHexesInQuantization() const;

    /**
     * @brief Struct to store information about path expanded during MC search
     */
    struct WeaklyMonotonousPath
    {
        OVM::VertexHandle n; // Current node
        int singPathIdx1;    // Index of source singular link
        UVWDir dirs1;        // Direction of extension of source

        map<UVWDir, vector<OVM::EdgeHandle>> dir2walkedArcs;
        vector<OVM::HalfEdgeHandle> path;
        // For efficiency dont store separate paths for each direction. Instead store the possible
        // monotonous directions that can still be realized given the current path
        UVWDir monotonousDirs;
        UVWDir walkedDirs; // Opposite dirs of these are never in monotonousDirs

        bool branchedOff; // Whether the path has left the starting singular link

        Vec3i deltaMin;
        Vec3i deltaMax;
        // If this is within [deltaMin, deltaMax] upon connecting 2 non-regular vertices, its a separation violation
        // If this goes outside of [deltaMin, deltaMax] for each coordinate anywhere on the path, we can stop tracing
        // this path, as there is no weakly monotonous connection to an endpoint within [deltaMin, deltaMax] anymore
        Vec3i delta;

        OVM::CellHandle bRefCurrent;
        Transition transCurrent;

        Q length = 0;
    };

    /**
     * @brief Compare paths based on their length
     */
    struct GreaterPathLengthCompare
    {
        bool operator()(const WeaklyMonotonousPath& p1, const WeaklyMonotonousPath& p2) const;
    };

    /**
     * @brief This will find 0-length paths connecting any singularlink-singularlink or singularlink-boundary pairs
     *        and store the MC arcs with associated +/- sign info into \p nonZeroSumArcs . The sum of these (including
     *        sign) may not be 0. This way, separation can be enforced. It is not guaranteed, that this method will find
     *        all violations. You have to iteratively call it and optimize using GUROBI until \p nonZeroSumArcs is
     *        returned empty.
     *
     * @param singularLinks IN: a collection of singular links that can be precomputed
     * @param nonZeroSumArcs OUT: MC arcs with associated +/- sign info into \p nonZeroSumArcs . The sum of these arcs'
     *                            lengths (including sign) may not be 0.
     * @return RetCode SUCCESS or error code
     */
    RetCode findSeparationViolatingPaths(const vector<SingularLink>& singularLinks,
                                         vector<vector<std::pair<int, OVM::EdgeHandle>>>& nonZeroSumArcs) const;

  protected:

    /**
     * @brief Given an initialization \p pStart , trace all possible paths monotonous in at least
     *        one coordinate until a separation violation is encountered and registered in \p nonZeroSumArcs
     *
     * @param singPath1 IN: singular link to use as start for MC graph search
     * @param nonZeroSumArcs IN/OUT: new separation violation is registered here
     * @return RetCode SUCCESS or error code
     */
    RetCode traceExhaustPaths(const SingularLink& singPath1,
                              vector<vector<std::pair<int, OVM::EdgeHandle>>>& nonZeroSumArcs) const;

    /**
     * @brief During monotonous path expansion, this is used to determine the set of valid
     *        outgoing halfedges and their reference blocks/transitions, so that the path of unfolded
     *        blocks always parametrically overlaps (i.e. 0-length distance away) with the source MC element.
     *
     * @param currentP IN: current path state
     * @param b2trans OUT: map of transition to each block reachable within P0 of \p currentP
     * @param ha2bRef OUT: mapping of outgoing next halfedges to their reference blocks/transitions
     * @return RetCode SUCCESS or error code
     */
    RetCode determineNextHalfedges(const WeaklyMonotonousPath& currentP,
                                   map<OVM::CellHandle, Transition>& b2trans,
                                   map<OVM::HalfEdgeHandle, vector<OVM::CellHandle>>& ha2bRef) const;

    /**
     * @brief Check if source of \p currentP and \p ha (must be connected to endpoint of
     *        \p currentP ) overlap parametrically.
     *
     * @param currentP IN: path currently expanded
     * @param ha IN: halfarc connecting onward from currentP.n
     * @param bRef IN: reference block of \p ha
     * @param trans IN: reference transition of \p bRef
     * @return true if overlapping parametrically
     * @return false else
     */
    bool checkArcOverlap(const WeaklyMonotonousPath& currentP,
                         const OVM::HalfEdgeHandle& ha,
                         const OVM::CellHandle& bRef,
                         const Transition& trans) const;

    /**
     * @brief Check if source of \p currentP and \p hp (must be connected to endpoint of
     *        \p currentP ) overlap parametrically.
     *
     * @param currentP IN: path currently expanded
     * @param hp IN: halfpatch connecting onward from currentP.n
     * @param trans IN: reference transition of cell incident on opposite patch of \p hp
     * @return true if overlapping parametrically
     * @return false else
     */
    bool checkPatchOverlap(const WeaklyMonotonousPath& currentP,
                           const OVM::HalfFaceHandle& hp,
                           const Transition& trans) const;

    /**
     * @brief Check if the current path tip is still contained in P0 region.
     *
     * @param currentP IN: current path
     * @return true if contained in P0 region
     * @return false else
     */
    bool checkP0Containment(const WeaklyMonotonousPath& currentP) const;

    /**
     * @brief Check whether bounding boxes overlap, and if so, in how many dimensions they are
     *        actually overlapping (and not just touching).
     *
     * @param min1 IN: Min point of bbox1
     * @param max1 IN: Max point of bbox1
     * @param min2 IN: Min point of bbox2
     * @param max2 IN: Max point of bbox2
     * @param nDimFullOverlap OUT: in how many dimensions they are actually overlapping (and not just touching)
     * @return true if overlapping
     * @return false else
     */
    static bool
    bboxOverlap(const Vec3i& min1, const Vec3i& max1, const Vec3i& min2, const Vec3i& max2, int& nDimFullOverlap);
};

} // namespace qgp3d

#endif
