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
    OVM::VertexHandle vFrom; // Path start vertex
    OVM::VertexHandle vTo;   // Path end vertex

    OVM::CellHandle tetFrom; // Path start reference tet
    OVM::CellHandle tetTo;   // Path end reference tet

    vector<OVM::HalfFaceHandle> pathHalffaces; // Halffaces traversed by path between start and end

    Vec3i offset; // Integer offset between start and end in coordinate system of tetFrom
};
} // namespace qgp3d

#endif
