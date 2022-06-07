#include "QGP3D/ConstraintWriter.hpp"

#include "QGP3D/ConstraintExtractor.hpp"

#include <iomanip>

namespace qgp3d
{

ConstraintWriter::ConstraintWriter(const TetMeshProps& meshProps, const std::string& fileName)
    : TetMeshNavigator(meshProps), _fileName(fileName), _os()
{
}

ConstraintWriter::RetCode ConstraintWriter::writeTetPathConstraints()
{
    _os = std::ofstream(_fileName);
    auto ret = checkFile();
    if (ret != RetCode::SUCCESS)
        return ret;
    ConstraintExtractor extr(_meshPropsC);

    auto haSequences = extr.getSingularitySkeletonArcs();
    auto tetPathConstraints = extr.getTetPathConstraints(haSequences);

    LOG(INFO) << "Writing " << tetPathConstraints.size() << " constraints to file " << _fileName;

    auto& mcMeshProps = *_meshPropsC.get<MC_MESH_PROPS>();

    auto& mc = mcMeshProps.mesh;
    int nHexes = 0;
    for (auto b : mc.cells())
    {
        int nBlockHexes = 1;
        for (auto dir : {UVWDir::NEG_U_NEG_V, UVWDir::NEG_U_NEG_W, UVWDir::NEG_V_NEG_W})
        {
            int arcLen = 0;
            for (auto a : mcMeshProps.ref<BLOCK_EDGE_ARCS>(b).at(dir))
                arcLen += mcMeshProps.get<ARC_INT_LENGTH>(a);
            nBlockHexes *= arcLen;
        }
        nHexes += nBlockHexes;
    }

    _os << nHexes << std::endl;
    _os << tetPathConstraints.size() << std::endl;
    for (auto& path : tetPathConstraints)
    {
        _os << path.vFrom << " " << path.vTo << std::endl;
        _os << path.offset << std::endl;
        _os << path.pathOrigTets.size() << " ";
        for (auto it = path.pathOrigTets.begin(); it != path.pathOrigTets.end() - 1; it++)
            _os << it->idx() << " ";
        _os << path.pathOrigTets.back() << std::endl;
    }

    return SUCCESS;
}

ConstraintWriter::RetCode ConstraintWriter::checkFile() const
{
    if (!_os.good())
    {
        LOG(ERROR) << "Could not write to file " << _fileName;
        return FILE_INACCESSIBLE;
    }
    return SUCCESS;
}

} // namespace qgp3d
