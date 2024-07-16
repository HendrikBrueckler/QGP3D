#ifndef QGP3D_QUANTIZER_HPP
#define QGP3D_QUANTIZER_HPP

#include <MC3D/Mesh/TetMeshProps.hpp>

#include "QGP3D/PathConstraint.hpp"

namespace qgp3d
{
using namespace mc3d;

/**
 * @brief Class managing the generation of quantization constraints (used to transform a
 *        seamless parametrization into an integer-grid map).
 *
 */
class Quantizer
{
  public:
    enum RetCode
    {
        SUCCESS = 0,
        MISSING_CHART = 4,          // Tet without parametrization
        INVALID_CHART = 5,          // Error in parametrization
        INVALID_SINGULARITY = 8,    // Singular curve not axis-aligned
        MC_CONSTRUCTION_ERROR = 21, // Motorcycle Complex construction failed
        MC_QUANTIZATION_ERROR = 22, // Motorcycle Complex quantization failed
    };

    /**
     * @brief Create an instance managing the generation of quantization constraints (used to transform a
     *        seamless parametrization into an integer-grid map).
     *
     * @param tetMesh IN: tet mesh
     */
    Quantizer(const TetMesh& tetMesh);

    /**
     * @brief Set the seamless parametrization for a given tet corner
     *
     * @param tet IN: tet
     * @param corner IN: vertex incident on \p tet
     * @param uvw IN: parametric coordinates of \p corner in chart of \p tet
     */
    void setParam(const CH& tet, const VH& corner, const Vec3d& uvw);

    /**
     * @brief Mark feature faces (no feature is default).
     *
     * @param f IN: face
     * @param isFeature IN: whether \p f is a feature
     */
    void setFeature(const FH& f, bool isFeature);

    /**
     * @brief Mark feature edges (no feature is default).
     *
     * @param e IN: edge
     * @param isFeature IN: whether \p e is a feature
     */
    void setFeature(const EH& e, bool isFeature);

    /**
     * @brief Mark feature vertices (no feature is default).
     *
     * @param v IN: vertex
     * @param isFeature IN: whether \p v is a feature
     */
    void setFeature(const VH& v, bool isFeature);

    /**
     * @brief Generate constraints for the quantization of the input seamless parametrization.
     *
     * @param scaleFactor IN: factor by which to scale the input seamless parametrization
     *                        (0.0 = coarsest possible IGM)
     * @param constraints OUT: complete set of integer spacing constraints fixing all integer DOFs
     * @param nHexes OUT: number of hexahedra in hex mesh implied by quantization
     * @param maxSecondsIQP IN: optional time limit for IQP solvers. 0 disables IQP solvers and uses greedy solve only.
     * @return RetCode SUCCESS or error code
     */
    RetCode quantize(double scaleFactor, vector<PathConstraint>& constraints, int& nHexes, int maxSecondsIQP = 300);

  private:
    const TetMesh& _tetMesh; // ref to input mesh
    TetMesh _meshCopy;       // internal copy
    MCMesh _mcMesh;          // Motorcycle Complex mesh
    TetMeshProps _meshProps; // property manager for quantization
};

} // namespace qgp3d

#endif
