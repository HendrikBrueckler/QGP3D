#ifndef QGP3D_STRUCTUREPRESERVER_HPP
#define QGP3D_STRUCTUREPRESERVER_HPP

#include <MC3D/Mesh/MCMeshNavigator.hpp>

namespace qgp3d
{
using namespace mc3d;

/**
 * @brief Class handling separation and structure preservation for quantization
 */
class StructurePreserver : public virtual MCMeshNavigator
{
  protected:
    // Forward declarations for easier overview
    struct WeaklyMonotonousPath;               // Tracks separable paths in fixed-singularity setting
    struct GreaterMonotonousPathLengthCompare; // For priority queue ordering

    struct IntegerGridGraphPath;               // Tracks separable paths in flexible-singularity setting
    struct IntegerGridGraphNode;               // Nodes of the integer grid graph
    struct Precursor;                          // Precursor backtracking information within integer-grid graph
    struct GreaterIntegerGridGraphPathCompare; // For priority queue ordering

  public:
    enum RetCode
    {
        SUCCESS = 0,
    };

    /**
     * @brief Create new instance checking the separation of singularities, boundaries and features in a quantization
     *
     * @param meshProps IN: mesh with quantized MC
     * @param maxSingIdx IN: maximum singularity index when allowing flexible singularities
     */
    StructurePreserver(TetMeshProps& meshProps, int maxSingIdx = 2);

    /**
     * @brief Count the number of hex cells in a mesh extractable from a quantized parametrization
     *        conforming to the current quantization
     *
     * @return int number of hex cells
     */
    long numHexesInQuantization() const;

    /**
     * @brief Get all simple constraints violated by the current quantization
     *
     * @param varLowerBound IN: lower bound of variables
     * @param nonZeroSumArcs OUT: MC arcs with associated +/- sign info into \p nonZeroSumArcs . The sum of these arcs'
     *                            lengths (including sign) may not be 0.
     */
    RetCode violatedSimpleConstraints(double varLowerBound, vector<vector<pair<int, EH>>>& nonZeroSumArcs) const;

    /**
     * @brief This will find any 0-length paths that imply a change of structure under the quantization metric,
     *        and store the path's arcs with associated +/- sign info into \p nonZeroSumArcs.
     *        The sum of these (including sign) must be >0. It is not guaranteed that this method will find all
     *        violations at once, but it will always some if existing.
     *
     * @param nonZeroSumArcs OUT: MC arcs with associated +/- sign info into \p nonZeroSumArcs . The sum of these arcs'
     *                            lengths (including sign) may not be 0.
     * @return RetCode SUCCESS or error code
     */
    RetCode violatedStructuralConstraints(vector<vector<pair<int, EH>>>& nonZeroSumArcs);

    /**
     * @brief Return all paths previously found by this StructurePreserver via violatedStructuralConstraints
     *
     * @return const vector<vector<pair<int, EH>>>& previous separation violating paths
     */
    const vector<vector<pair<int, EH>>>& previousStructuralConstraints()
    {
        return _allStructuralConstraints;
    }

    /**
     * @brief Return positive-coefficient-only parts of paths previously found by this StructurePreserver via
     *        violatedStructuralConstraints
     *
     * @return const vector<vector<pair<int, EH>>>& failsafe separation violating paths
     */
    const vector<vector<pair<int, EH>>>& previousFailsafeStructuralConstraints()
    {
        return _failsafeStructuralConstraints;
    }

  protected:
    /**
     * @brief Given an initialization \p pStart , trace all possible paths monotonous in at least
     *        one coordinate until a separation violation is encountered and registered in \p nonZeroSumArcs
     *
     * @param pStart IN: initial path and information about startpoint (start critical entity)
     * @param allowedBranchOffs IN: all arcs allowed for leaving the initial critical entity
     * @param allowedPreBranchArcs IN: all arcs allowed before leaving the initial critical entity
     * @param isCriticalNode IN: markings of critical nodes
     * @param isCriticalArc IN: markings of critical arcs
     * @param isCriticalPatch IN: markings of critical patches
     * @param nonZeroSumArcs IN/OUT: new structural constraints are registered here
     * @param failsafeNonZeroSumArcs IN/OUT: positive-coefficient-only version of \p nonZeroSumArcs
     * @param returnOnViolation IN: whether to return on first separation violation encountered
     * @return RetCode SUCCESS or error code
     */
    void traceExhaustPaths(const WeaklyMonotonousPath& pStart,
                           const set<HEH>& allowedBranchOffs,
                           const set<EH>& allowedPreBranchArcs,
                           vector<vector<pair<int, EH>>>& nonZeroSumArcs,
                           vector<vector<pair<int, EH>>>& failsafeNonZeroSumArcs,
                           bool returnOnViolation = true) const;

    /**
     * @brief During monotonous path expansion, this is used to determine the set of valid
     *        outgoing halfearcs and their reference blocks/transitions, so that the path of unfolded
     *        blocks always parametrically overlaps (i.e. 0-length distance away) with the source MC element.
     *
     * @param currentP IN: current path state
     * @param b2trans OUT: map of transition to each block reachable within P0 of \p currentP
     * @param ha2bRef OUT: mapping of outgoing next halfarcs to their reference blocks/transitions
     * @return RetCode SUCCESS or error code
     */
    void determineNextHalfarcs(const WeaklyMonotonousPath& currentP,
                               map<CH, Transition>& b2trans,
                               map<HEH, vector<CH>>& ha2bRef) const;

    /**
     * @brief Given an initial integer-grid graph node \p c0 , trace all zero-paths from \p c0 staying within the union
     *        of starting critical entity and any other entity sharing a boundary with it.
     *
     * @param c0 IN: start graphnode
     * @param bRef0 IN: initial reference block
     * @param i0Min IN: lower bound of 1D/2D interval of covered integer grid nodes
     * @param i0Max IN: upper bound of 1D/2D interval of covered integer grid nodes
     * @param potentiallyUnproblematicNodes IN: all integer grid graph nodes that might be unproblematic
     * @param unproblematicNodes OUT: actually unproblematic nodes
     */
    void traceUnproblematicPaths(const IntegerGridGraphNode& c0,
                                 const CH& bRef0,
                                 const Vec3i& i0Min,
                                 const Vec3i& i0Max,
                                 const set<IntegerGridGraphNode>& potentiallyUnproblematicNodes,
                                 set<IntegerGridGraphNode>& unproblematicNodes);

    /**
     * @brief Given an initial integer-grid graph node \p c0 , trace all possible paths until a structural
     *        violation is encountered and registered in \p nonZeroSumArcs or the search space is exhausted.
     *
     * @param c0 IN: start graphnode
     * @param bRef0 IN: initial reference block
     * @param i0Min IN: lower bound of 1D/2D interval of covered integer grid nodes
     * @param i0Max IN: upper bound of 1D/2D interval of covered integer grid nodes
     * @param unproblematicNodes IN: actually unproblematic nodes
     * @param nonZeroSumArcs IN/OUT: new structural constraints are registered here
     * @param failsafeNonZeroSumArcs IN/OUT: positive-coefficient-only version of \p nonZeroSumArcs
     * @return true if violation found
     * @return false otherwise
     */
    bool traceExhaustPaths(const IntegerGridGraphNode& c0,
                           const CH& bRef0,
                           const Vec3i& i0Min,
                           const Vec3i& i0Max,
                           const set<IntegerGridGraphNode>& unproblematicNodes,
                           vector<vector<pair<int, EH>>>& nonZeroSumArcs,
                           vector<vector<pair<int, EH>>>& failsafeNonZeroSumArcs) const;

    /**
     * @brief Check if zero path alters the core domain structure (features or topology).
     *
     * @param pCurrent IN: current zero path
     * @param unproblematicNodes IN: unproblematic nodes
     * @param nonZeroSumArcs IN/OUT: new structural constraint is registered here
     * @param failsafeNonZeroSumArcs IN/OUT: positive-coefficient-only version of \p nonZeroSumArcs
     * @return true if path is structure altering
     * @return false otherwise
     */
    bool isStructureAlteringZeroPath(const IntegerGridGraphPath& pCurrent,
                                     const set<IntegerGridGraphNode>& unproblematicNodes,
                                     vector<vector<pair<int, EH>>>& nonZeroSumArcs,
                                     vector<vector<pair<int, EH>>>& failsafeNonZeroSumArcs) const;

    /**
     * @brief Reconstructs an arc path from an integer-grid graph path.
     *
     * @param pCurrent IN: input integer grid graph path
     * @param pathHaDirs OUT: per-block half-arc sequences and directions
     * @param excludeNonZeroArcs IN: if true, arcs with non-zero integer length are not inserted into arc path
     */
    void reconstructArcPaths(const IntegerGridGraphPath& pCurrent,
                             vector<vector<pair<HEH, UVWDir>>>& pathHaDirs,
                             bool excludeNonZeroArcs) const;

    /**
     * @brief Precompute internal integer-grid mapping and transition tables used by path tracing.
     */
    void precompute();

    /**
     * @brief Find paths that lead to singularity clusters with out-of-bounds total index and add constraints
     *        that enforce splitting cluster into less offending subclusters.
     *
     * @param nonZeroSumArcs IN/OUT: new structural constraints are registered here
     * @param failsafeNonZeroSumArcs IN/OUT: positive-coefficient-only version of \p nonZeroSumArcs
     */
    void findSingularityIndexViolationPaths(vector<vector<pair<int, EH>>>& nonZeroSumArcs,
                                            vector<vector<pair<int, EH>>>& failsafeNonZeroSumArcs) const;

    /**
     * @brief Explore integer-grid graph paths for detecting singularity index violations.
     *
     * @param c0 IN: starting graph node
     * @param bRef0 IN: starting block reference
     * @param i0 IN: starting integer-grid position
     * @param dirsParallel0 IN: initial allowed parallel directions
     * @param startID IN: critical link ID from which path starts
     * @param nonZeroSumArcs IN/OUT: new structural constraints are registered here
     * @param failsafeNonZeroSumArcs IN/OUT: positive-coefficient-only version of \p nonZeroSumArcs
     * @param visitedLinks IN/OUT: track links already visited
     * @return true if a violating cluster is found
     * @return false otherwise
     */
    bool traceSingularityIndexPaths(const IntegerGridGraphNode& c0,
                                    const CH& bRef0,
                                    const Vec3i& i0,
                                    const UVWDir dirsParallel0,
                                    int startID,
                                    vector<vector<pair<int, EH>>>& nonZeroSumArcs,
                                    vector<vector<pair<int, EH>>>& failsafeNonZeroSumArcs,
                                    set<int>& visitedLinks) const;

    /**
     * @brief Detect (topologically) incontractible zero-loops passing through c0. Detected by checking
     *        simple-connectedness of the expanded portion of the integer-grid graph.
     *
     * @param c0 IN: start graph node
     * @param bRef0 IN: reference block for start node
     * @param i0Min IN: lower bound of 1D/2D interval of covered integer grid nodes
     * @param i0Max IN: upper bound of 1D/2D interval of covered integer grid nodes
     * @param nsVisited IN/OUT: node set visited by search
     * @param nonZeroSumArcs IN/OUT: new structural constraints are registered here
     * @param failsafeNonZeroSumArcs IN/OUT: positive-coefficient-only version of \p nonZeroSumArcs
     * @return true if incontractible loop found
     * @return false otherwise
     */
    bool findIncontractibleZeroLoops(const IntegerGridGraphNode& c0,
                                     const CH& bRef0,
                                     const Vec3i& i0Min,
                                     const Vec3i& i0Max,
                                     set<VH>& nsVisited,
                                     vector<vector<pair<int, EH>>>& nonZeroSumArcs,
                                     vector<vector<pair<int, EH>>>& failsafeNonZeroSumArcs) const;

    /**
     * @brief Fetch transitions around an entity (node, arc, or patch) for path expansion.
     *
     * @param c IN: the integer grid graph node
     * @param bRef IN: reference block
     * @param bRef2trans OUT: map of blocks to their transition vectors
     */
    void fetchTransitionsAroundEntity(const IntegerGridGraphNode& c,
                                      const CH& bRef,
                                      map<CH, vector<Transition>>& bRef2trans) const;

    /**
     * @brief Initialize a starting path for integer grid graph traversal.
     *
     * @param c0 IN: start graph node
     * @param bRef0 IN: initial reference block
     * @param i0Min IN: lower bound of interval
     * @param i0Max IN: upper bound of interval
     * @param includePrecursors IN: whether to include empty precursor
     * @return IntegerGridGraphPath initialized path
     */
    IntegerGridGraphPath initializePathStart(const IntegerGridGraphNode& c0,
                                             const CH& bRef0,
                                             const Vec3i& i0Min,
                                             const Vec3i& i0Max,
                                             bool includePrecursors = false) const;

    /**
     * @brief Find the node in the given entity that is closest to the given UVW position.
     *
     * @param c IN: the integer grid graph node
     * @param bRef IN: reference block
     * @param uvw IN: the UVW position
     * @return pair<VH, double> the closest node and its squared distance
     */
    std::pair<VH, double> findClosestNode(const IntegerGridGraphNode& c, const CH& bRef, const Vec3Q& uvw) const;

    /**
     * @brief Check if source of \p currentP and \p ha (must be connected to endpoint of
     *        \p currentP ) overlap parametrically.
     *
     * @param currentP IN: path currently expanded
     * @param ha IN: halfarc connecting onward from currentP.n
     * @param bRef IN: reference block of \p ha
     * @param trans IN: reference transition of \p bRef
     * @return true if overlapping parametrically
     * @return false otherwise
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
     * @return false otherwise
     */
    bool checkPatchOverlap(const WeaklyMonotonousPath& currentP, const HFH& hp, const Transition& trans) const;

    /**
     * @brief Check if the current path tip is still contained in P0 region.
     *
     * @param currentP IN: current path
     * @return true if contained in P0 region
     * @return false otherwise
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
     * @return false otherwise
     */
    static bool bboxOverlap(const Vec3i& min1, const Vec3i& max1, const Vec3i& min2, const Vec3i& max2);

    /**
     * @brief Visualize a constraint path in an external viewer (only if compiled with viewer support).
     *
     * @param constraintPath IN: arc-list constraint path to be shown
     */
    void debugViewConstraintPath(const vector<pair<int, EH>>& constraintPath) const;

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
    struct GreaterMonotonousPathLengthCompare
    {
        bool operator()(const WeaklyMonotonousPath& p1, const WeaklyMonotonousPath& p2) const
        {
            return p1.length > p2.length || (p1.length == p2.length && p1.path.size() > p2.path.size());
        }
    };

    struct IntegerGridGraphNode
    {
        IntegerGridGraphNode()
        {
        }
        IntegerGridGraphNode(VH n_) : n(n_)
        {
        }
        IntegerGridGraphNode(EH a_) : a(a_)
        {
        }
        IntegerGridGraphNode(FH p_) : p(p_)
        {
        }

        VH n = {};
        EH a = {};
        FH p = {};

        bool operator==(const IntegerGridGraphNode& other) const
        {
            return other.n == n && other.a == a && other.p == p;
        }
        bool operator!=(const IntegerGridGraphNode& other) const
        {
            return !(*this == other);
        }

        bool is_valid() const
        {
            return n.is_valid() || a.is_valid() || p.is_valid();
        }

        vector<VH> nodes(const MCMesh& mcMesh) const
        {
            if (n.is_valid())
                return {n};
            if (a.is_valid())
            {
                auto ns = mcMesh.edge_vertices(a);
                return {ns[0], ns[1]};
            }
            if (p.is_valid())
                return mcMesh.get_halfface_vertices(mcMesh.halfface_handle(p, 0));
            return {};
        }

        bool operator<(const IntegerGridGraphNode& other) const
        {
            return Vec3i(n.idx(), a.idx(), p.idx()) < Vec3i(other.n.idx(), other.a.idx(), other.p.idx());
        }
    };

    struct Precursor
    {
        StructurePreserver::IntegerGridGraphNode cPrev = {};
        Transition transPrev = {};
        CH bRefPrev = {};

        Vec3Q deltaToPrev = {};
        Transition transToPrev = {};
    };

    struct IntegerGridGraphPath
    {
        Vec3i i = Vec3i(0, 0, 0);
        Vec3i iMin = Vec3i(0, 0, 0);
        Vec3i iMax = Vec3i(0, 0, 0);
        Vec3Q u = Vec3Q(0, 0, 0);
        Vec3Q u0 = Vec3Q(0, 0, 0);
        StructurePreserver::IntegerGridGraphNode c = {};
        CH bRef = {};
        Transition trans = {};
        UVWDir dirs = UVWDir::NONE;
        UVWDir dirsParallel = UVWDir::NONE;
        int hops = 0;
        double length = 0.0;
        vector<Precursor> precursors;

        int criticalStartID = -1;
        double partLength = 0.0;
        UVWDir dirsTotal = UVWDir::NONE;
    };

    struct GreaterIntegerGridGraphPathCompare
    {
        bool operator()(const IntegerGridGraphPath& p1, const IntegerGridGraphPath& p2) const
        {
            return p1.length > p2.length
                   || (p1.length == p2.length
                       && (p1.hops > p2.hops
                           || (p1.hops == p2.hops
                               && (decompose(p1.dirs, DIM_1_DIRS).size() > decompose(p2.dirs, DIM_1_DIRS).size()
                                   || (decompose(p1.dirs, DIM_1_DIRS).size() == decompose(p2.dirs, DIM_1_DIRS).size()
                                       && (decompose(p1.dirs & -p1.dirs, DIM_1_DIRS).size()
                                           > decompose(p2.dirs & -p2.dirs, DIM_1_DIRS).size()))))));
            // Last 2 lines probably overengineered, but oh well
        }
    };

    // Criticality bookkeeping:
    vector<CriticalEntity> criticalEntities; // All critical entities (nodes/arcs/patches)
    vector<bool> isCriticalArc;              // Flag per MC arc if it is critical for structure-preservation
    vector<bool> isCriticalArcFromRegion;    // Flag per MC arc if it is critical for potential boundary-padding
    vector<bool> isCriticalNode;             // Flag per MC vertex if it is critical for structure-preservation
    vector<bool> isCriticalPatch;            // Flag per MC patch if it is critical for structure-preservation

    // Singularity cluster bookkeeping:
    vector<CriticalEntity> singularEntities; // All singular links
    map<EH, int> a2singularLinkIdx;          // Map arc to singular link index
    vector<bool> isSingularArc;              // Flag per MC arc if it is singular

    // Singularity mode settings:
    bool fixedSingularities = true;  // Whether singularities are fixed (separated and incollapsible)
    bool paddedSingularities = true; // Whether singularities are padded (separated at least from boundary)

    // Block / entity mapping caches for efficient graph traversal:
    vector<map<VH, Vec3i>> _b2n2igm;    // Block -> (vertex -> integer grid position) mapping
    vector<map<VH, Vec3Q>> _b2n2uvw;    // Block -> (vertex -> UVW coordinate) mapping
    vector<map<EH, Vec3i>> _b2a2igmMin; // Block -> (arc -> integer grid min bound) mapping
    vector<map<EH, Vec3i>> _b2a2igmMax; // Block -> (arc -> integer grid max bound) mapping
    vector<map<FH, Vec3i>> _b2p2igmMin; // Block -> (patch -> integer grid min bound) mapping
    vector<map<FH, Vec3i>> _b2p2igmMax; // Block -> (patch -> integer grid max bound) mapping

    // Transition caches for neighbor lookup:
    map<pair<CH, VH>, map<CH, vector<Transition>>> _bn2transitionsAroundNode; // Node-centric transitions by block
    map<pair<CH, EH>, map<CH, vector<Transition>>> _ba2transitionsAroundArc;  // Arc-centric transitions by block

    // Result collection for structural constraints paths:
    vector<vector<std::pair<int, EH>>> _allStructuralConstraints;      // All explicitly added structural constraints
    vector<vector<std::pair<int, EH>>> _failsafeStructuralConstraints; // Positive-sign-only structural constraints

    // Singularity index bounds
    const int _minSingIdx = -1; // Fixed minimum singularity index
    const int _maxSingIdx;      // Configurable maximum singularity index

    // Path-length precision control:
    const bool _considerDoubleIntLength = true; // Whether to assume integer arc lengths twice as long (for edge cases)
};

} // namespace qgp3d

#endif
