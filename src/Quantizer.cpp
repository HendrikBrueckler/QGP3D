#include "QGP3D/Quantizer.hpp"

#include <MC3D/Interface/MCGenerator.hpp>

#include "QGP3D/ConstraintExtractor.hpp"
#include "QGP3D/MCQuantizer.hpp"

#include "TrulySeamless3D/trulyseamless.h"

namespace qgp3d
{

Quantizer::Quantizer(const TetMesh& tetMesh) : _tetMesh(tetMesh), _meshCopy(_tetMesh), _mcMesh(), _meshProps(_meshCopy, _mcMesh)
{
    _meshProps.allocate<CHART>();
}

void Quantizer::setParam(const OVM::CellHandle& tet, const OVM::VertexHandle& corner, const Vec3d& uvw)
{
    assert(set<OVM::VertexHandle>(_meshCopy.tet_vertices(tet).first, _meshCopy.tet_vertices(tet).second).count(corner) != 0);
    _meshProps.ref<CHART>(tet)[corner] = Vec3Q(uvw);
}

void Quantizer::setFeature(const OVM::FaceHandle& f, bool isFeature)
{
    if (!_meshProps.isAllocated<IS_FEATURE_F>())
        _meshProps.allocate<IS_FEATURE_F>(false);
    _meshProps.set<IS_FEATURE_F>(f, isFeature);
}

void Quantizer::setFeature(const OVM::EdgeHandle& e, bool isFeature)
{
    if (!_meshProps.isAllocated<IS_FEATURE_E>())
        _meshProps.allocate<IS_FEATURE_E>(false);
    _meshProps.set<IS_FEATURE_E>(e, isFeature);
}

void Quantizer::setFeature(const OVM::VertexHandle& v, bool isFeature)
{
    if (!_meshProps.isAllocated<IS_FEATURE_V>())
        _meshProps.allocate<IS_FEATURE_V>(false);
    _meshProps.set<IS_FEATURE_V>(v, isFeature);
}

Quantizer::RetCode Quantizer::quantize(double scaleFactor, vector<PathConstraint>& constraints, int& nHexes)
{
    MCGenerator mcGen(_meshProps);
    bool invalidCharts = false;
    for (auto tet : _meshCopy.cells())
    {
        if (mcGen.volumeUVW(tet) == 0)
        {
            LOG(ERROR) << "Degenerate parametrization in tet " << tet;
            invalidCharts = true;
        }
        else if (mcGen.volumeUVW(tet) < 0)
        {
            LOG(ERROR) << "Flipped parametrization in tet " << tet;
            invalidCharts = true;
        }
    }

    if (invalidCharts)
        return INVALID_CHART;

    try
    {
        TS3D::TrulySeamless3D sanitizer(_meshCopy);
        for (auto tet : _meshCopy.cells())
            for (auto v : _meshCopy.tet_vertices(tet))
                sanitizer.parameter(tet, v) = Vec3Q2d(_meshProps.ref<CHART>(tet).at(v));
        if (!sanitizer.init() || !sanitizer.sanitize(0.0, true))
        {
            LOG(ERROR) << "Sanitization failed";
            return INVALID_CHART;
        }
        for (auto tet : _meshCopy.cells())
            for (auto v : _meshCopy.tet_vertices(tet))
                _meshProps.ref<CHART>(tet).at(v) = sanitizer.parameter(tet, v);
    }
    catch (std::runtime_error& e)
    {
        LOG(ERROR) << "Sanitization failed: " << e.what();
        return INVALID_CHART;
    }
    LOG(INFO) << "Sanitization successful";

    auto retGen = mcGen.traceMC(true, true, false, true);

    switch (retGen)
    {
    case MCGenerator::MISSING_CHART:
        return MISSING_CHART;
    case MCGenerator::INVALID_SINGULARITY:
        return INVALID_SINGULARITY;
    case MCGenerator::SPAWNING_FAILED:
    case MCGenerator::TRACING_FAILED:
    case MCGenerator::SPLITTING_FAILED:
    case MCGenerator::BUILDING_MC_FAILED:
        return MC_CONSTRUCTION_ERROR;
    case MCGenerator::SUCCESS:
        break; // noop
    }

    auto retQuant = MCQuantizer(_meshProps).quantizeArcLengths(scaleFactor, true, true);

    switch (retQuant)
    {
    case MCQuantizer::SOLVER_ERROR:
        return MC_QUANTIZATION_ERROR;
    case MCQuantizer::SUCCESS:
        break; // noop
    }

    ConstraintExtractor extr(_meshProps);

    auto haSequences = extr.getCriticalSkeletonArcs();
    auto tetPathConstraints = extr.getTetPathConstraints(haSequences);

    constraints.clear();
    for (auto& tetPath : tetPathConstraints)
    {
        constraints.emplace_back();
        constraints.back().vFrom = tetPath.vFrom;
        constraints.back().vTo = tetPath.vTo;
        constraints.back().tetFrom = tetPath.pathOrigTets.front();
        constraints.back().tetTo = tetPath.pathOrigTets.back();
        constraints.back().offset = tetPath.offset;

        auto& tets = tetPath.pathOrigTets;
        for (int i = 0; i + 1 < (int)tets.size(); i++)
        {
            auto tet1 = tets[i];
            auto tet2 = tets[i + 1];
            for (auto hf : _tetMesh.cell_halffaces(tet1))
                if (_tetMesh.incident_cell(_tetMesh.opposite_halfface_handle(hf)) == tet2)
                    constraints.back().pathHalffaces.push_back(hf);
        }
        assert(constraints.back().pathHalffaces.size() + 1 == tetPath.pathOrigTets.size());
    }

    auto& mcMeshProps = *_meshProps.get<MC_MESH_PROPS>();
    nHexes = 0;
    for (auto b : mcMeshProps.mesh.cells())
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

    return SUCCESS;
}

} // namespace qgp3d
