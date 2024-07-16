#ifndef QGP3D_WITHOUT_IQP

#ifndef QGP3D_IQPQUANTIZER_HPP
#define QGP3D_IQPQUANTIZER_HPP

#include <MC3D/Mesh/MCMeshManipulator.hpp>

#include "QGP3D/SeparationChecker.hpp"

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
     * @param sep IN: separationchecker to use for lazy separation. May come with predetermined separation constraints
     */
    IQPQuantizer(TetMeshProps& meshProps, SeparationChecker& sep);

    /**
     * @brief Compute quantization
     *
     * @param scaling IN: scale target lengths by this factor for quantization
     * @param varLowerBound IN: lower bound for arc lengths
     * @param maxSeconds IN: time limit in seconds
     * @return RetCode SUCCESS or error code
     */
    RetCode quantize(double scaling = 1.0, double varLowerBound = 0.0, int maxSeconds = 300.0);

  protected:
    SeparationChecker& _sep; // Separation checker given from outside
};

} // namespace qgp3d
#endif

#endif
