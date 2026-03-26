#ifndef QGP3D_WITH_GUROBI
#ifndef QGP3D_WITHOUT_IQP

#include "QGP3D/IQP/BonminIQPSolver.hpp"

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseQR>

#include "BonCbc.hpp"
#include "BonEcpCuts.hpp"
#include "BonIpoptSolver.hpp"
#include "BonOACutGenerator2.hpp"
#include "BonOaNlpOptim.hpp"
#include "BonOsiTMINLPInterface.hpp"

namespace qgp3d
{
namespace impl
{
BonminIQPSolver::BonminIQPSolver(
    TetMeshProps& meshProps, double varLowerBound, double maxSeconds)
    : TetMeshNavigator(meshProps), TetMeshManipulator(meshProps), MCMeshNavigator(meshProps),
      MCMeshManipulator(meshProps), BaseIQPSolver(varLowerBound, maxSeconds),
      _instance(new BonminIQP(meshProps, varLowerBound))
{
}

void BonminIQPSolver::setupDefaultOptions()
{
    for (auto* setup: {&_setupExact, &_setupQuick})
    {
        setup->initializeOptionsAndJournalist();
        setup->options()->SetNumericValue("bonmin.time_limit", _maxSeconds); // changes bonmin's time limit

        // Problem properties
        setup->options()->SetBoolValue("hessian_constant", true);
        setup->options()->SetBoolValue("jac_c_constant", true);
        setup->options()->SetBoolValue("jac_d_constant", true);
        setup->options()->SetBoolValue("expect_infeasible_problem", false);
        setup->options()->SetStringValue("expect_infeasible_problem", "no");

        // Finetuning
        setup->options()->SetStringValue("algorithm", "B-BB");
        setup->options()->SetStringValue("mu_oracle", "loqo");
        setup->options()->SetStringValue("nlp_failure_behavior", "fathom"); // default
        setup->options()->SetStringValue("node_comparison", "best-bound");  // default
        setup->options()->SetStringValue("variable_selection",
                                        "qp-strong-branching");   // default: strong-branching
        setup->options()->SetStringValue("nlp_solver", "Ipopt");   // default
        setup->options()->SetStringValue("warm_start", "optimum"); // default: none

        // Tolerances
        setup->options()->SetNumericValue("integer_tolerance", 1e-3);
    }
    _setupExact.options()->SetNumericValue("allowable_fraction_gap",
                                    1e-4); // obj considered optimal if gap small
    _setupQuick.options()->SetNumericValue("allowable_fraction_gap",
                                    0.99); // obj considered optimal if gap small
    _setupQuick.options()->SetIntegerValue("solution_limit", 1);
}

void BonminIQPSolver::setupVariables()
{
    // One integer length variable per arc
    _instance->_arc2idx.clear();
    for (EH arc : mcMeshProps().mesh().edges())
    {
        int idx = _instance->_arc2idx.size();
        _instance->_arc2idx[arc] = idx;
    }
}

void BonminIQPSolver::setupObjective(QuadraticObjective& obj)
{
    auto& mc = mcMeshProps().mesh();

    _instance->_objectiveBase = obj.f0;
    _instance->_objectiveLin = obj.grad0;
    _instance->_hessianSparse = obj.hessian;
}

void BonminIQPSolver::setupConstraints()
{
    auto& mc = mcMeshProps().mesh();
    const bool SINGULARITIES_FIXED = !meshProps().isAllocated<ALGO_VARIANT>() || (meshProps().get<ALGO_VARIANT>() % 2);

    vector<bool> isCriticalNode, isCriticalArc, isCriticalPatch;
    vector<CriticalEntity> criticalEntities;
    map<EH, int> a2criticalLinkIdx;
    map<FH, int> p2criticalRegionIdx;
    map<VH, vector<int>> n2criticalLinksOut;
    map<VH, vector<int>> n2criticalLinksIn;
    getCriticalEntities(isCriticalNode,
                        isCriticalArc,
                        isCriticalPatch,
                        criticalEntities,
                        a2criticalLinkIdx,
                        p2criticalRegionIdx,
                        n2criticalLinksOut,
                        n2criticalLinksIn,
                        true,
                        SINGULARITIES_FIXED,
                        true);

    _instance->_constraints.clear();

    {
        Eigen::MatrixXd equalityMatrix(2 * mc.n_logical_faces(), _instance->_arc2idx.size());
        equalityMatrix.setZero();
        int globalIdx = 0;
        // Each patches opposite arc lengths must match
        for (FH patch : mc.faces())
        {
            HFH hp = mc.halfface_handle(patch, 0);
            if (mc.is_boundary(hp))
                hp = mc.opposite_halfface_handle(hp);

            auto dir2orderedHas = halfpatchHalfarcsByDir(hp);
            assert(dir2orderedHas.size() == 4);
            UVWDir dirHpNormal = halfpatchNormalDir(hp);

            for (UVWDir dirSide :
                 decompose(~(dirHpNormal | -dirHpNormal), {UVWDir::POS_U, UVWDir::POS_V, UVWDir::POS_W}))
            {
                set<EH> esA;
                for (HEH ha : dir2orderedHas[dirSide])
                    esA.insert(mc.edge_handle(ha));
                set<EH> esB;
                for (HEH ha : dir2orderedHas[-dirSide])
                    esB.insert(mc.edge_handle(ha));
                if (esA == esB)
                    continue;

                for (HEH ha : dir2orderedHas[dirSide])
                    equalityMatrix(globalIdx, _instance->_arc2idx.at(mc.edge_handle(ha))) = 1;

                for (HEH ha : dir2orderedHas[-dirSide])
                    equalityMatrix(globalIdx, _instance->_arc2idx.at(mc.edge_handle(ha))) = -1;

                globalIdx++;
            }
        }
        Eigen::FullPivLU<Eigen::MatrixXd> decomp(equalityMatrix.transpose());
        Eigen::MatrixXd image = decomp.image(equalityMatrix.transpose()).transpose();

        for (int i = 0; i < image.rows(); i++)
        {
            vector<int> vars;
            vector<double> coeffs;
            for (int j = 0; j < image.cols(); j++)
                if (std::abs(image(i, j)) > 1e-6)
                {
                    vars.push_back(j);
                    coeffs.push_back(std::round(image(i, j)));
                }
            _instance->_constraints.push_back({vars, coeffs, 0.0, true});
        }
    }

    if (_varLowerBound < 0.0)
    {
        // Enforce no block with negative extension along any axis
        for (CH b : mc.cells())
            for (UVWDir dir : {UVWDir::POS_V_POS_W, UVWDir::POS_U_POS_W, UVWDir::POS_U_POS_V})
            {
                vector<int> vars;
                vector<double> coeffs;
                auto& arcs = mcMeshProps().ref<BLOCK_EDGE_ARCS>(b).at(dir);
                for (EH a : arcs)
                {
                    vars.push_back(_instance->_arc2idx.at(a));
                    coeffs.push_back(1.0);
                }
                _instance->_constraints.push_back({vars, coeffs, 0.0, false});
            }
    }

    // Enforce critical links of length >= 1
    for (auto& criticalLink : criticalEntities)
    {
        if (criticalLink.pathHas.empty())
            continue;

        vector<int> vars;
        vector<double> coeffs;
        for (HEH ha : criticalLink.pathHas)
        {
            vars.push_back(_instance->_arc2idx.at(mc.edge_handle(ha)));
            coeffs.push_back(1.0);
        }

        _instance->_constraints.push_back({vars, coeffs, 1.0, false});
    }
}

void BonminIQPSolver::finalizeSetup()
{
    _setupExact.initialize(_instance);
    _setupQuick.initialize(_instance);
#ifdef NDEBUG
    _setupExact.continuousSolver()->messageHandler()->setLogLevel(0);
    _setupExact.nonlinearSolver()->messageHandler()->setLogLevel(0);
    _setupQuick.continuousSolver()->messageHandler()->setLogLevel(0);
    _setupQuick.nonlinearSolver()->messageHandler()->setLogLevel(0);
#endif
}

void BonminIQPSolver::enableQuickSolve()
{
    _quickSolve = true;
}

void BonminIQPSolver::enableExactSolve()
{
    _quickSolve = false;
}


void BonminIQPSolver::useCurrentAsWarmStart()
{
    // No extra effort needed
}

BaseIQPSolver::RetCode BonminIQPSolver::solve()
{
    try
    {
        Bonmin::Bab bb;
        bb(_quickSolve ? _setupQuick : _setupExact);
        std::cout << std::flush;

        _lastObjective = bb.bestObj();
        DLOG(INFO) << "Bonmin solved IQP with final objective value of " << _lastObjective;

        int status = bb.mipStatus();
        if (status != Bonmin::Bab::FeasibleOptimal && status != Bonmin::Bab::Feasible)
        {
            LOG(ERROR) << "Bad status return by Bonmin solver";
            return BaseIQPSolver::SOLVER_ERROR;
        }

        auto solution = bb.bestSolution();
        for (EH arc : mcMeshProps().mesh().edges())
            mcMeshProps().set<ARC_INT_LENGTH>(arc, (int)std::round(solution[_instance->_arc2idx.at(arc)]));
    }
    catch (...)
    {
        LOG(ERROR) << "Bonmin exception";
        return SOLVER_ERROR;
    }

    return BaseIQPSolver::SUCCESS;
}

double BonminIQPSolver::objectiveValue() const
{
    return _lastObjective;
}

void BonminIQPSolver::addConstraints(const vector<vector<pair<int, EH>>>& nonZeroSum)
{
    for (auto& aColl : nonZeroSum)
    {
        vector<int> vars;
        vector<double> coeffs;
        double rhs = 1.0;
        for (auto sign2a : aColl)
        {
            if (sign2a.second.is_valid())
            {
                vars.push_back(_instance->_arc2idx.at(sign2a.second));
                coeffs.push_back(sign2a.first);
            }
            else
                rhs = sign2a.first;
        }
        _instance->_constraints.push_back({vars, coeffs, rhs, false});
    }
    _setupQuick.nonlinearSolver()->setModel(_instance);
    _setupExact.nonlinearSolver()->setModel(_instance);
}

BonminIQPSolver::BonminIQP::BonminIQP(const TetMeshProps& meshProps, double varLowerBound)
    : TetMeshNavigator(meshProps), MCMeshNavigator(meshProps), Bonmin::TMINLP(),
      _varLowerBound(varLowerBound)
{
}

BonminIQPSolver::BonminIQP::BonminIQP(const BonminIQPSolver::BonminIQP& other)
    : TetMeshNavigator(other.meshProps()), MCMeshNavigator(other.meshProps()), Bonmin::TMINLP(),
      _varLowerBound(other._varLowerBound), _arc2idx(other._arc2idx),
      _objectiveBase(other._objectiveBase), _objectiveLin(other._objectiveLin), _hessianSparse(other._hessianSparse),
      _constraints(other._constraints)
{
}

bool BonminIQPSolver::BonminIQP::get_variables_types(Ipopt::Index n, Bonmin::TMINLP::VariableType* var_types)
{
    for (int i = 0; i < n; i++)
        var_types[i] = INTEGER;
    return true;
}

bool BonminIQPSolver::BonminIQP::get_variables_linearity(Ipopt::Index n, Ipopt::TNLP::LinearityType* var_types)
{
    for (int i = 0; i < n; i++)
        var_types[i] = Ipopt::TNLP::NON_LINEAR;
    return true;
}

bool BonminIQPSolver::BonminIQP::get_constraints_linearity(Ipopt::Index m, Ipopt::TNLP::LinearityType* const_types)
{
    for (int i = 0; i < m; i++)
        const_types[i] = Ipopt::TNLP::LINEAR;
    return true;
}

bool BonminIQPSolver::BonminIQP::get_nlp_info(Ipopt::Index& n,
                                       Ipopt::Index& m,
                                       Ipopt::Index& nnz_jac_g,
                                       Ipopt::Index& nnz_h_lag,
                                       Ipopt::TNLP::IndexStyleEnum& index_style)
{
    n = _arc2idx.size();

    m = _constraints.size();

    nnz_h_lag = 0;
    for (int i = 0; i < _hessianSparse.outerSize(); ++i)
        for (Eigen::SparseMatrix<double>::InnerIterator it(_hessianSparse, i); it; ++it)
            if (it.row() >= it.col())
                ++nnz_h_lag;

    int nNonZero = 0;
    for (auto& constraint : _constraints)
        nNonZero += constraint.vars.size();

    nnz_jac_g = nNonZero;

    index_style = Ipopt::TNLP::C_STYLE;
    return true;
}

bool BonminIQPSolver::BonminIQP::get_bounds_info(
    Ipopt::Index n, Ipopt::Number* x_l, Ipopt::Number* x_u, Ipopt::Index m, Ipopt::Number* g_l, Ipopt::Number* g_u)
{
    (void)m;
    for (int i = 0; i < n; i++)
    {
        x_l[i] = _varLowerBound;
        x_u[i] = DBL_MAX;
    }
    for (int i = 0; i < (int)_constraints.size(); i++)
    {
        if (_constraints[i].equality)
            g_l[i] = g_u[i] = _constraints[i].rhs;
        else
        {
            g_l[i] = _constraints[i].rhs;
            g_u[i] = DBL_MAX;
        }
    }
    return true;
}

bool BonminIQPSolver::BonminIQP::get_starting_point(Ipopt::Index n,
                                             bool init_x,
                                             Ipopt::Number* x,
                                             bool init_z,
                                             Ipopt::Number* z_L,
                                             Ipopt::Number* z_U,
                                             Ipopt::Index m,
                                             bool init_lambda,
                                             Ipopt::Number* lambda)
{
    (void)n;
    (void)init_x;
    (void)init_z;
    (void)z_L;
    (void)z_U;
    (void)m;
    (void)init_lambda;
    (void)lambda;
    if (mcMeshProps().isAllocated<ARC_INT_LENGTH>())
        for (EH arc : mcMeshProps().mesh().edges())
            x[_arc2idx.at(arc)] = mcMeshProps().get<ARC_INT_LENGTH>(arc);
    else
    {
        Eigen::VectorXd x = Eigen::VectorXd::Zero(mcMeshProps().mesh().n_logical_edges());

        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> ldlt(_hessianSparse);
        if (ldlt.info() != Eigen::Success)
        {
            LOG(ERROR) << "Ill conditioned objective hessian for some reason";
            throw std::logic_error("");
        }

        Eigen::VectorXd sol = ldlt.solve(-_objectiveLin);

        int i  = 0;
        for (EH a: mcMeshProps().mesh().edges())
            x[i++] = sol(i);
    }

    return true;
}

bool BonminIQPSolver::BonminIQP::eval_f(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Number& obj_value)
{
    (void)new_x;
    Eigen::Map<const Eigen::VectorXd> xd(x, n);
    obj_value = 0.5 * xd.dot(_hessianSparse * xd) + _objectiveLin.dot(xd) + _objectiveBase;

    return true;
}

bool BonminIQPSolver::BonminIQP::eval_grad_f(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Number* grad_f)
{
    (void)new_x;
    Eigen::Map<const Eigen::VectorXd> xd(x, n);
    Eigen::Map<Eigen::VectorXd> g(grad_f, n);
    g = _hessianSparse * xd + _objectiveLin;
    return true;
}

bool BonminIQPSolver::BonminIQP::eval_g(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Index m, Ipopt::Number* g)
{
    (void)n;
    (void)new_x;
    (void)m;
    for (int i = 0; i < (int)_constraints.size(); i++)
    {
        g[i] = 0.0;
        for (int j = 0; j < (int)_constraints[i].vars.size(); j++)
            g[i] += x[_constraints[i].vars[j]] * _constraints[i].coeffs[j];
    }

    return true;
}

bool BonminIQPSolver::BonminIQP::eval_jac_g(Ipopt::Index n,
                                     const Ipopt::Number* x,
                                     bool new_x,
                                     Ipopt::Index m,
                                     Ipopt::Index nele_jac,
                                     Ipopt::Index* iRow,
                                     Ipopt::Index* jCol,
                                     Ipopt::Number* values)
{
    (void)n;
    (void)x;
    (void)new_x;
    (void)m;
    (void)nele_jac;
    if (values == nullptr)
    {
        int gi = 0;
        for (int i = 0; i < (int)_constraints.size(); i++)
            for (int j = 0; j < (int)_constraints[i].vars.size(); j++)
            {
                iRow[gi] = i;
                jCol[gi] = _constraints[i].vars[j];
                gi++;
            }
    }
    else
    {
        int gi = 0;
        for (int i = 0; i < (int)_constraints.size(); i++)
            for (int j = 0; j < (int)_constraints[i].vars.size(); j++)
            {
                values[gi] = _constraints[i].coeffs[j];
                gi++;
            }
    }
    return true;
}

bool BonminIQPSolver::BonminIQP::eval_h(Ipopt::Index n,
                                 const Ipopt::Number* x,
                                 bool new_x,
                                 Ipopt::Number obj_factor,
                                 Ipopt::Index m,
                                 const Ipopt::Number* lambda,
                                 bool new_lambda,
                                 Ipopt::Index nele_hess,
                                 Ipopt::Index* iRow,
                                 Ipopt::Index* jCol,
                                 Ipopt::Number* values)
{
    (void)n;
    (void)x;
    (void)new_x;
    (void)m;
    (void)lambda;
    (void)new_lambda;
    (void)nele_hess;
    if (values == NULL)
    {
        int gi = 0;
        for (int i = 0; i < _hessianSparse.outerSize(); ++i)
            for (Eigen::SparseMatrix<double>::InnerIterator it(_hessianSparse, i); it; ++it)
                if (it.row() >= it.col())
                {
                    iRow[gi] = it.row();
                    jCol[gi] = it.col();
                    ++gi;
                }
    }
    else
    {
        int gi = 0;
        for (int i = 0; i < _hessianSparse.outerSize(); ++i)
            for (Eigen::SparseMatrix<double>::InnerIterator it(_hessianSparse, i); it; ++it)
                if (it.row() >= it.col())
                {
                    values[gi] = obj_factor * it.value();
                    ++gi;
                }
    }
    return true;
}

void BonminIQPSolver::BonminIQP::finalize_solution(Bonmin::TMINLP::SolverReturn status,
                                            Ipopt::Index n,
                                            const Ipopt::Number* x,
                                            Ipopt::Number obj_value)
{
    (void)n;
    (void)x;
    (void)status;
    (void)obj_value;
}

} // namespace impl
} // namespace qgp3d

#endif
#endif
