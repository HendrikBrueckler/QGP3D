#ifndef QGP3D_SEPARATIONCHECKER_HPP
#define QGP3D_SEPARATIONCHECKER_HPP

#include <MC3D/Mesh/MCMeshNavigator.hpp>

namespace qgp3d
{
using namespace mc3d;

/**
 * @brief Class handling separation in a quantization
 */
class SeparationChecker : public virtual MCMeshNavigator
{
  public:
    /**
     * @brief Create new instance checking the separation of singularities, boundaries and features in a quantization
     *
     * @param meshProps IN: mesh with quantized MC
     */
    SeparationChecker(TetMeshProps& meshProps);

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
        VH n;         // Current node
        UVWDir dirs1; // Direction of extension of source

        map<UVWDir, vector<EH>> dir2walkedArcs;
        vector<HEH> path;
        // For efficiency dont store separate paths for each direction. Instead store the possible
        // monotonous directions that can still be realized given the current path
        UVWDir monotonousDirs;
        UVWDir walkedDirs; // Opposite dirs of these are never in monotonousDirs

        bool branchedOff; // Whether the path has left the starting critical link

        Vec3i deltaMin;
        Vec3i deltaMax;
        // If this is within [deltaMin, deltaMax] upon connecting 2 non-regular vertices, its a separation violation
        // If this goes outside of [deltaMin, deltaMax] for each coordinate anywhere on the path, we can stop tracing
        // this path, as there is no weakly monotonous connection to an endpoint within [deltaMin, deltaMax] anymore
        Vec3i delta;

        CH bRefCurrent;
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
     * @brief This will find 0-length paths connecting any criticallink-criticallink or criticallink-boundary pairs
     *        and store the MC arcs with associated +/- sign info into \p nonZeroSumArcs . The sum of these (including
     *        sign) may not be 0. This way, separation can be enforced. It is not guaranteed, that this method will find
     *        all violations. You have to iteratively call it and optimize using GUROBI until \p nonZeroSumArcs is
     *        returned empty.
     *
     * @param criticalLinks IN: a collection of critical links that can be precomputed
     * @param nonZeroSumArcs OUT: MC arcs with associated +/- sign info into \p nonZeroSumArcs . The sum of these arcs'
     *                            lengths (including sign) may not be 0.
     * @return RetCode SUCCESS or error code
     */
    void findSeparationViolatingPaths(vector<CriticalLink>& criticalLinks,
                                      const vector<bool>& arcIsCritical,
                                      const vector<bool>& nodeIsCritical,
                                      const vector<bool>& patchIsCritical,
                                      vector<vector<pair<int, EH>>>& nonZeroSumArcs);

    /**
     * @brief Return all paths previously found by this separationchecker via findSeparationViolatingPaths
     *
     * @return const vector<vector<pair<int, EH>>>& previous separation violating paths
     */
    const vector<vector<pair<int, EH>>>& previousSeparationViolatingPaths()
    {
        return _allSeparatingPaths;
    }

    /**
     * @brief Return positive-coefficient-only parts of paths previously found by this separationchecker via
     *        findSeparationViolatingPaths
     *
     * @return const vector<vector<pair<int, EH>>>& failsafe separation violating paths
     */
    const vector<vector<pair<int, EH>>>& failsafeSeparationViolatingPaths()
    {
        return _failsafeSeparatingPaths;
    }

  protected:
    /**
     * @brief Given an initialization \p pStart , trace all possible paths monotonous in at least
     *        one coordinate until a separation violation is encountered and registered in \p nonZeroSumArcs
     *
     * @param criticalPath1 IN: critical link to use as start for MC graph search
     * @param nonZeroSumArcs IN/OUT: new separation violation is registered here
     * @param failsafeNonZeroSumArcs IN/OUT: positive-coefficient-only version of \p nonZeroSumArcs
     * @return RetCode SUCCESS or error code
     */
    void traceExhaustPaths(const CriticalLink& criticalPath1,
                           const vector<bool>& arcIsCritical,
                           const vector<bool>& nodeIsCritical,
                           const vector<bool>& patchIsCritical,
                           vector<vector<std::pair<int, EH>>>& nonZeroSumArcs,
                           vector<vector<std::pair<int, EH>>>& failsafeNonZeroSumArcs) const;

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
    void determineNextHalfedges(const WeaklyMonotonousPath& currentP,
                                map<CH, Transition>& b2trans,
                                map<HEH, vector<CH>>& ha2bRef) const;

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
    bool
    checkArcOverlap(const WeaklyMonotonousPath& currentP, const HEH& ha, const CH& bRef, const Transition& trans) const;

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
    bool checkPatchOverlap(const WeaklyMonotonousPath& currentP, const HFH& hp, const Transition& trans) const;

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
     * @return true if overlapping
     * @return false else
     */
    static bool bboxOverlap(const Vec3i& min1, const Vec3i& max1, const Vec3i& min2, const Vec3i& max2);

    vector<vector<std::pair<int, EH>>> _allSeparatingPaths;      // Full set of separation violating paths
    vector<vector<std::pair<int, EH>>> _failsafeSeparatingPaths; // Full set of separation violating paths (monotonous)
};

} // namespace qgp3d

#endif
