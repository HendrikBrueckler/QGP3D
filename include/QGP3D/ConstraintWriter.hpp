#ifndef QGP3D_CONSTRAINTWRITER_HPP
#define QGP3D_CONSTRAINTWRITER_HPP

#include <MC3D/Mesh/TetMeshNavigator.hpp>

#include <fstream>
#include <map>
#include <string>

namespace qgp3d
{
using namespace mc3d;

class ConstraintWriter : public TetMeshNavigator
{
  public:
    enum RetCode
    {
        SUCCESS = 0,
        FILE_INACCESSIBLE = 1, // Could not access file
    };

    /**
     * @brief Create a writer, that reads a mesh with parametrization \p meshProps
     *  and writes it in .hexex format to \p meshProps.
     *
     * @param meshProps IN: mesh to write
     * @param fileName IN: file to write to
     */
    ConstraintWriter(const TetMeshProps& meshProps, const std::string& fileName);

    /**
     * @brief Write constraints/paths to specified file
     *
     * @return RetCode
     */
    RetCode writeTetPathConstraints();

  private:
    const std::string _fileName;
    std::ofstream _os;
    bool _exact;

    map<int, int> _vtx2idx; // Internal map of vtx mesh idx to vtx output idx

    /**
     * @brief  Check if file is writeable
     *
     * @return RetCode SUCCESS or FILE_INACCESSIBLE
     */
    RetCode checkFile() const;
};

} // namespace qgp3d

#endif
