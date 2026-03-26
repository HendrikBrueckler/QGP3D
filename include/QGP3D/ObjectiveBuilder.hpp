#ifndef QGP3D_OBJECTIVEBUILDER_HPP
#define QGP3D_OBJECTIVEBUILDER_HPP

#include "QGP3D/ObjectiveFunction.hpp"
#include <MC3D/Mesh/MCMeshNavigator.hpp>

namespace qgp3d
{
using namespace mc3d;

/**
 * @brief Class helping with the formulation of a quantization objective function
 */
class ObjectiveBuilder : public virtual MCMeshNavigator
{
  public:
    enum RetCode
    {
        SUCCESS = 0,
    };

    /**
     * @brief Create new instance checking the separation of singularities, boundaries and features in a quantization
     *
     * @param meshProps IN: mesh with quantized MC
     */
    ObjectiveBuilder(TetMeshProps& meshProps);

    /**
     * @brief Get a quantization objective function simply measuring the sum over squared arc length deviations.
     *
     * @param scaling IN: Assume parametrization to be scaled by this factor
     * @return QuadraticObjective quantization objective
     */
    QuadraticObjective arcLengthDeviationObjective(double scaling = 1.0) const;

    /**
     * @brief Get a quantization objective function measuring the sum over squared feature-feature-path length
     * deviations. These paths are coarse approximations of the edges in a Delaunay tetrahedrization of features.
     *
     * @param scaling IN: Assume parametrization to be scaled by this factor
     * @param individualArcRegularization IN: obj = (1-this) * pathDeviationSum + this * arcDeviationSum
     * @return QuadraticObjective quantization objective
     */
    QuadraticObjective simplifiedDistortionObjective(double scaling = 1.0,
                                                     double individualArcRegularization = 0.01) const;

    /**
     * @brief Get the nodes that serve as seeds for constructing the (parametrically) Voronoi diagram of critical
     * (feature) entities for \ref simplifiedDistortionObjective.
     *
     * @param seeds OUT: seed nodes
     * @return RetCode SUCCESS
     */
    RetCode getCriticalSeeds(set<VH>& seeds) const;

    /**
     * @brief From given \p seeds grow cells of Voronoi diagram on node-arc graph of MC mesh (essentially Dijkstra).
     *
     * @param seeds IN: seed nodes
     * @param n2precursors OUT: the precursor from which each node was conquered
     * @return RetCode SUCCESS
     */
    RetCode growVoronoiCellsMC(const set<VH>& seeds, map<VH, HEH>& n2precursors) const;

    /**
     * @brief Given the paths from which each node was conquered, and the corresponding seed node from which it was
     * conquered, group the interface arcs (connecting nodes conquered from different seeds) by local connectivity.
     *
     * @param seeds IN: seed nodes
     * @param n2precursors IN: the precursor from which each node was conquered
     * @param seedPairs2interfaces OUT: for each pair of adjacent seeds, the spatially distinct interfaces via which
     * they are connected
     * @return RetCode SUCCESS
     */
    RetCode groupInterfacesMC(const set<VH>& seeds,
                              const map<VH, HEH>& n2precursors,
                              map<pairTT<VH>, vector<set<EH>>>& seedPairs2interfaces) const;

    /**
     * @brief For each interface, get a representative path, representing the Delaunay tetrahedrization edge.
     *
     * @param n2precursors IN: the precursor from which each node was conquered
     * @param seedPairs2interfaces in: for each pair of adjacent seeds, the spatially distinct interfaces via which they
     * are connected
     * @param paths OUT: paths
     * @param pathWeights OUT: weights per path (volume integral factor for distortion objective sum)
     * @return RetCode SUCCESS
     */
    RetCode getInterfaceTraversingPathsMC(const map<VH, HEH>& n2precursors,
                                          const map<pairTT<VH>, vector<set<EH>>>& seedPairs2interfaces,
                                          vector<vector<HEH>>& paths,
                                          vector<double>& pathWeights) const;

    /**
     * @brief Debug visualization of seed nodes as spheres and interface/non-interface arcs as polylines (cylinders).
     *
     * @param seeds IN: seed nodes
     * @param n2precursors IN: the precursor from which each node was conquered
     */
    void debugViewVoronoiCells(const set<VH>& seeds, const map<VH, HEH>& n2precursors) const;

    /**
     * @brief Debug visualization of seed nodes as spheres and interface/non-interface arcs as polylines (cylinders).
     *
     * @param seeds IN: seed nodes
     * @param n2precursors IN: the precursor from which each node was conquered
     */
    void debugViewPaths(const map<pairTT<VH>, vector<set<EH>>>& seedPairs2interfaces,
                        const vector<vector<HEH>>& paths) const;
};

} // namespace qgp3d

#endif
