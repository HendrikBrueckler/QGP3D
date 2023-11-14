#ifndef QGP3D_CONSTRAINTEXTRACTOR_HPP
#define QGP3D_CONSTRAINTEXTRACTOR_HPP

#include "MC3D/Mesh/MCMeshManipulator.hpp"

#include "QGP3D/PathConstraint.hpp"
namespace qgp3d
{
using namespace mc3d;

/**
 * @brief Class that manages the extraction of the MC based arc-tree of singular links and boundary regions
 */
class ConstraintExtractor : public MCMeshNavigator
{
  public:
    struct TetPathConstraint
    {
        VH vFrom;                // Path start vertex
        VH vTo;                  // Path end vertex
        vector<CH> pathOrigTets; // Halffaces traversed by path between start and end
        Vec3i offset;            // Integer offset between start and end in coordinate system of tetFrom
    };

    /**
     * @brief Creates an instance that manages NodeTree extraction from \p meshProps
     *
     * @param meshProps IN: mesh whose singular-node NodeTree to extract
     */
    ConstraintExtractor(const TetMeshProps& meshProps);

    /**
     * @brief Get the segments of the arc-sceleton connecting all singular links of the MC.
     *        Segments are maximum sequences of (half)arcs connecting 2 (possibly pseudo-)singular
     *        nodes. Pseudosingular nodes are virtual (fake) singular nodes inserted on cyclic singular
     *        links.
     *
     * @return vector<vector<HEH>> singularity arc skeleton segments
     */
    vector<vector<HEH>> getCriticalSkeletonArcs();

    /**
     * @brief From a set of connected arc segments together with an underlying quantization, generate
     *        3D constraints, that fix the difference in parametrization of the two end nodes when accumulated
     *        along the arc segment.
     *
     * @param paths IN: connected arc segments
     * @return vector<PathConstraint> constraints
     */
    vector<TetPathConstraint> getTetPathConstraints(const vector<vector<HEH>>& paths);

    /**
     * @brief Get the cut surfaces, that cut the MC into a topological ball
     *
     * @param cutSurfaces OUT: manifold patch collections that together form the cut surface
     * @param p2cutSurface OUT: tagging of patches (patch idx => cut surface idx)
     */
    void getCutSurfaces(vector<vector<FH>>& cutSurfaces, vector<int>& p2cutSurface);

  protected:
    /**
     * @brief Determine which nodes are critical and which are pseudocritical (inserted on critical elements
     *        that have no actual critical node)
     */
    void assignNodeTypes();

    /**
     * @brief Mark patches of the cut graph (cutting mesh to topological ball)
     */
    void markCutPatches();

    /**
     * @brief Count cut patches incident on arc
     *
     * @param a IN: arc
     * @return int number of cut surface patches incident on \p a
     */
    int countCutPatches(const EH& a) const;

    /**
     * @brief Whether node is on an actual cut of the cut graph
     *
     * @param n IN: node
     * @return true if on cut
     * @return false else
     */
    bool isInActualCut(const VH& n) const;

    /**
     * @brief Count number of skeleton arcs incident on node
     *
     * @param n IN: node
     * @return int number of skeleton arcs incident on \p n
     */
    int countIncidentSkeletonArcs(const VH& n) const;

    /**
     * @brief Get a path from start node to any target node through block b.
     *        \p n and at least one element of \p nsTarget must be in \p b .
     *
     * @param n IN: start node
     * @param nsTarget IN: target nodes
     * @param b IN: block
     * @return list<HEH> halfarcs connecting \p n to any of \p nsTarget
     */
    list<HEH> pathThroughBlock(const VH& n, const set<VH>& nsTarget, const CH& b) const;

    /**
     * @brief Get one cyclic halfarc path through each manifold piece of the cut surface.
     *
     * @return vector<vector<HEH>> cyclic halfarc paths
     */
    vector<vector<HEH>> getCutSurfaceCycles();

    /**
     * @brief Path of tets connecting the current tet to any tet (including that tet) incident on a target halfedge
     *
     * @param tetStart IN: start tet
     * @param heTarget IN: target halfedge
     * @param vConn IN: pivot vertex
     * @return list<CH> tet path
     */
    list<CH> connectingTetPath(const CH& tetStart, const HEH& heTarget, const VH& vConn) const;

    /**
     * @brief Path of tets connecting the current tet to any tet (including that tet) incident on a target halfface
     *
     * @param tetStart IN: start tet
     * @param hfTarget IN: target halfface
     * @param vConn IN: pivot vertex
     * @return list<CH>  tet path
     */
    list<CH> connectingTetPath(const CH& tetStart, const HFH& hfTarget, const VH& vConn) const;

    /**
     * @brief Get the accumulated transition along a path of tets
     *
     * @param path IN: tet path
     * @return Transition accumulated transition along path
     */
    Transition transitionAlongPath(const list<CH>& path) const;

    /**
     * @brief Snap the endpoints of the given path (which might be non-original")
     *        to one of their enclosing original vertices (2 on edge, 3 on face).
     *        This assumes that both endpoints are either
     *
     * @param path IN: path with set endpoints, OUT: path with endpoints snapped to original vertices
     * @param pathTets IN: tet sequence between previous path endpoints, OUT: tet sequence between new endpoints
     * @param e2parent IN: edge split hierarchy
     * @param f2parent IN: face split hierarchy
     */
    void determineEquivalentEndpoints(TetPathConstraint& path,
                                      list<CH>& pathTets,
                                      const map<EH, EH>& e2parent,
                                      const map<FH, FH>& f2parent) const;

    /**
     * @brief Increment tet iterator until current element is incident on e
     *
     * @param e IN: edge
     * @param itTet IN: iterator in tet sequence, OUT: (possibly incremented) iterator
     * @param transOrigCurr IN: current transition (in orig charts), OUT: updated transition
     */
    void walkToMatchingTet(const EH& e, list<CH>::iterator& itTet, Transition& transOrigCurr) const;

    /**
     * @brief Get the coordinate index [0 u, 1 v, 2 w] of the direction in which a cyclic singularity
     *        located at \p vConn extends.
     *
     * @param tetStart IN: some tet incident on \p vConn
     * @param vConn IN: vertex on cyclic singularity
     * @param transOrigCurr IN: current transition inside coord system of \p tetStart
     * @return int coordinate index
     */
    int getCycleCoord(const CH& tetStart, const VH& vConn, const Transition& transOrigCurr) const;

    /**
     * @brief Get the coordinate indices [0 u, 1 v, 2 w] of the directions in which a borderless
     *        boundary region located at \p vConn extends.
     *
     * @param tetStart IN: some tet incident on \p vConn
     * @param vConn IN: vertex on borderless boundary region
     * @param transOrigCurr IN: current transition inside coord system of \p tetStart
     * @return std::pair<int, int> coordinate indices
     */
    std::pair<int, int> getBoundaryCoords(const CH& tetStart, const VH& vConn, const Transition& transOrigCurr) const;

    /**
     * @brief Build a dual spanning tree. Fails in case no leaf node exists.
     *
     * @param nRoot OUT: tree root
     * @param bRoot OUT: tree root block
     * @param cellTreePrecursor OUT: spanning tree encoded by precursor for each block
     * @param cutSurfaceToAvoid IN: wether a cut surface has been precomputed and the root node should avoid that surface
     * @return true if valid result
     * @return false if failure
     */
    bool buildBlockSpanningTree(VH& nRoot, CH& bRoot, vector<CH>& cellTreePrecursor, int cutSurfaceToAvoid);

    /**
     * @brief Find a cycle of halfarcs passing through a given node \p nStart incident on a cut surface patch \p pStart .
     *
     * @param nStart IN: node on cut surface
     * @param pStart IN: patch of cut surface
     * @param nRoot IN: root node to use
     * @param bRoot IN: root block to use
     * @param cellTreePrecursor IN: spanning tree encoded by precursor for each block
     * @return vector<HEH> ordered halfarc cycle
     */
    vector<HEH> cycleThroughCutSurface(const VH& nStart, const FH& pStart, const VH& nRoot, const CH& bRoot, const vector<CH>& cellTreePrecursor) const;

    /**
     * @brief Types of nodes marked as "singular"
     */
    enum ConstraintNodeType
    {
        NONE,
        NATIVE,
        ARTIFICIAL_ON_LINK,
        ARTIFICIAL_ON_BOUNDARY
    };

    VH _nRoot;                                      // Root node of dual MC spanning tree
    CH _bRoot;                                      // Root block of dual MC spanning tree
    vector<CH> _cellTreePrecursor;                  // Precursor in dual spanning tree
    vector<bool> _isSkeletonArc;                    // whether a given arc is part of the constraint skeleton
    vector<bool> _isCutPatch;                       // whether a given patch is part of the cut surface
    vector<int> _cutSurfaceID;                      // which manifold cut surface patch a given patch is part of
    vector<vector<FH>> _cutSurfaces;                // manifold collections of cut surface patches
    vector<ConstraintNodeType> _constraintNodeType; // tags for different types of constraint tree nodes
    vector<bool> _arcIsCritical;                    // tags for critical arcs
};

} // namespace mc3d

#endif
