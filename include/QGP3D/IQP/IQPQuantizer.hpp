#ifndef QGP3D_WITHOUT_IQP

#ifndef QGP3D_IQPQUANTIZER_HPP
#define QGP3D_IQPQUANTIZER_HPP

#include <MC3D/Mesh/MCMeshManipulator.hpp>

#include "QGP3D/ObjectiveFunction.hpp"
#include "QGP3D/StructurePreserver.hpp"

namespace qgp3d
{
using namespace mc3d;

/**
 * @brief Class to execute quantization using an IQP solver
 */
class IQPQuantizer : public virtual MCMeshManipulator
{
  public:
    enum RetCode
    {
        SUCCESS = 0,
        SOLVER_ERROR = 20,
    };

    /**
     * @brief Create an instance handling quantization via IQP solvers
     *
     * @param meshProps IN: mesh whose MC to quantize
     * @param sep IN: StructurePreserver to use for lazy separation. May come with predetermined separation constraints
     */
    IQPQuantizer(TetMeshProps& meshProps, StructurePreserver& sep, QuadraticObjective& obj);

    /**
     * @brief Compute quantization
     *
     * @param scaling IN: scale target lengths by this factor for quantization
     * @param varLowerBound IN: lower bound for arc lengths
     * @param maxSeconds IN: time limit in seconds
     * @return RetCode SUCCESS or error code
     */
    RetCode quantize(double varLowerBound = 0.0, int maxSeconds = 300.0);

  protected:
    StructurePreserver& _sep;  // Separation checker managed externally
    QuadraticObjective& _obj; // The objective function to minimize managed externally
};

} // namespace qgp3d
#endif

#endif
