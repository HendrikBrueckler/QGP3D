#ifndef QGP3D_PATHCONSTRAINT_HPP
#define QGP3D_PATHCONSTRAINT_HPP

namespace qgp3d
{
using namespace mc3d;

/**
 * @brief Path along which a fixed differential constraint (difference in param between start/end accumulated along
 *        path in coordinate system of tetFrom) is set up.
 */
struct PathConstraint
{
    VH vFrom; // Path start vertex
    VH vTo;   // Path end vertex

    CH tetFrom; // Path start reference tet
    CH tetTo;   // Path end reference tet

    vector<OVM::HalfFaceHandle> pathHalffaces; // Halffaces traversed by path between start and end

    Vec3i offset; // Integer offset between start and end in coordinate system of tetFrom
};
} // namespace qgp3d

#endif
