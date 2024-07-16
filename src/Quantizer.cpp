#include "QGP3D/Quantizer.hpp"

#include "QGP3D/ConstraintExtractor.hpp"
#include "QGP3D/IQP/IQPQuantizer.hpp"
#include "QGP3D/ISP/ISPQuantizer.hpp"

#include <MC3D/Interface/MCGenerator.hpp>

#include <TS3D/trulyseamless.h>

namespace qgp3d
{

Quantizer::Quantizer(const TetMesh& tetMesh)
    : _tetMesh(tetMesh), _meshCopy(_tetMesh), _mcMesh(), _meshProps(_meshCopy, _mcMesh)
{
    _meshProps.allocate<CHART>();
}

void Quantizer::setParam(const OVM::CellHandle& tet, const VH& corner, const Vec3d& uvw)
{
    assert(set<VH>(_meshCopy.tet_vertices(tet).first, _meshCopy.tet_vertices(tet).second).count(corner) != 0);
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

void Quantizer::setFeature(const VH& v, bool isFeature)
{
    if (!_meshProps.isAllocated<IS_FEATURE_V>())
        _meshProps.allocate<IS_FEATURE_V>(false);
    _meshProps.set<IS_FEATURE_V>(v, isFeature);
}

Quantizer::RetCode
Quantizer::quantize(double scaleFactor, vector<PathConstraint>& constraints, int& nHexes, int maxSecondsIQP)
{
#ifdef QGP3D_WITHOUT_IQP
    (void)maxSecondsIQP;
#endif
    MCGenerator mcGen(_meshProps);
    bool invalidCharts = false;
    for (auto tet : _meshCopy.cells())
    {
        Q vol = mcGen.rationalVolumeUVW(tet);
        if (vol == 0)
        {
            LOG(ERROR) << "Degenerate parametrization in tet " << tet;
            invalidCharts = true;
        }
        else if (vol < 0)
        {
            LOG(ERROR) << "Flipped parametrization in tet " << tet;
            invalidCharts = true;
        }
    }

    if (invalidCharts)
        return INVALID_CHART;

    if (_meshProps.isAllocated<IS_FEATURE_F>())
    {
        if (!_meshProps.isAllocated<IS_FEATURE_E>())
            _meshProps.allocate<IS_FEATURE_E>();
        // Make features consistent:
        // Mark each edge that has != 0 or != 2 feature patches as feature
        for (auto e : _meshCopy.edges())
        {
            if (_meshProps.get<IS_FEATURE_E>(e))
                continue;

            int nFeatureFaces = 0;
            for (auto f : _meshCopy.edge_faces(e))
                if (_meshProps.get<IS_FEATURE_F>(f))
                    nFeatureFaces++;
            if (nFeatureFaces != 2 && nFeatureFaces != 0)
                _meshProps.set<IS_FEATURE_E>(e, true);
        }
    }

    if (_meshProps.isAllocated<IS_FEATURE_E>())
    {
        if (!_meshProps.isAllocated<IS_FEATURE_V>())
            _meshProps.allocate<IS_FEATURE_V>();
        // Mark each vertex that has != 0 or != 2 feature edges as feature
        for (auto v : _meshCopy.vertices())
        {
            if (_meshProps.get<IS_FEATURE_V>(v))
                continue;

            int nFeatureEdges = 0;
            for (auto e : _meshCopy.vertex_edges(v))
                if (_meshProps.get<IS_FEATURE_E>(e))
                    nFeatureEdges++;
            if (nFeatureEdges != 2 && nFeatureEdges != 0)
                _meshProps.set<IS_FEATURE_V>(v, true);
        }
    }

    try
    {
        TS3D::TrulySeamless3D sanitizer(_meshCopy);
        for (auto tet : _meshCopy.cells())
            for (auto v : _meshCopy.tet_vertices(tet))
                sanitizer.setParam(tet, v, Vec3Q2d(_meshProps.ref<CHART>(tet).at(v)));
        if (_meshProps.isAllocated<IS_FEATURE_E>())
            for (auto e : _meshCopy.edges())
                if (_meshProps.get<IS_FEATURE_E>(e))
                    sanitizer.setFeature(e);
        if (_meshProps.isAllocated<IS_FEATURE_V>())
            for (auto v : _meshCopy.vertices())
                if (_meshProps.get<IS_FEATURE_V>(v))
                    sanitizer.setFeature(v);
        if (_meshProps.isAllocated<IS_FEATURE_F>())
            for (auto f : _meshCopy.faces())
                if (_meshProps.get<IS_FEATURE_F>(f))
                    sanitizer.setFeature(f);
        if (!sanitizer.init() || !sanitizer.sanitize(0.0, true))
        {
            LOG(ERROR) << "Sanitization failed";
            return INVALID_CHART;
        }
        for (auto tet : _meshCopy.cells())
            for (auto v : _meshCopy.tet_vertices(tet))
                _meshProps.ref<CHART>(tet).at(v) = sanitizer.getParam(tet, v);
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

    retGen = mcGen.reduceMC(true, true, true);
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

    SeparationChecker sep(_meshProps);
    ISPQuantizer isp(_meshProps, sep);
    isp.quantize(scaleFactor, -DBL_MAX);

#ifndef QGP3D_WITHOUT_IQP
    IQPQuantizer iqp(_meshProps, sep);
    IQPQuantizer::RetCode retIQP = IQPQuantizer::SUCCESS;
    if (maxSecondsIQP > 0)
        retIQP = iqp.quantize(scaleFactor, -DBL_MAX, maxSecondsIQP);

    (void)retIQP;
#endif

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
    for (auto b : mcMeshProps.mesh().cells())
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
