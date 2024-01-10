#include "QGP3D/MCQuantizer.hpp"

#include <queue>

#define INDIVIDUAL_ARC_FACTOR 0.01
#define TIME_LIMIT 600

#ifdef QGP3D_WITH_GUROBI

#include <gurobi_c++.h>

#else // QGP3D with bonmin#include "BonBonminSetup.hpp"

#include "BonCbc.hpp"
#include "BonIpoptSolver.hpp"
#include "BonOsiTMINLPInterface.hpp"
#include "BonTMINLP.hpp"

#include "BonEcpCuts.hpp"
#include "BonOACutGenerator2.hpp"
#include "BonOaNlpOptim.hpp"

#include "BonBonminSetup.hpp"

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace impl
{
using namespace Ipopt;
using namespace Bonmin;
using namespace mc3d;
using std::vector;

class QuantizationProblem : public TMINLP
{
  private:
    struct Constraint
    {
        Constraint(const vector<int>& vars_, const vector<double>& coeffs_, double rhs_, bool equality_)
            : vars(vars_), coeffs(coeffs_), rhs(rhs_), equality(equality_)
        {
        }

        vector<int> vars;
        vector<double> coeffs;
        double rhs;
        bool equality;
    };

    bool _printSol;
    double _scaling;
    bool _allowZeros;
    bool _allowNeg;
    int n_nonzeros_jac;
    int n_nonzeros_hess;
    TetMeshProps& _meshProps;

    vector<vector<double>> _hessianDense;
    Eigen::VectorXd _objectiveLin;
    Eigen::SparseMatrix<double> _hessianSparse;
    double _objectiveBase;
    vector<Constraint> _constraints;

  public:
    map<EH, int> _arc2idx;

    void addConstraint(const vector<int>& vars, const vector<double>& coeffs, double rhs, bool equality)
    {
        _constraints.push_back({vars, coeffs, rhs, equality});
    }

    QuantizationProblem(TetMeshProps& meshProps, double scaling, bool allowZeros, bool allowNeg)
        : _printSol(false), _scaling(scaling), _allowZeros(allowZeros), _allowNeg(allowNeg), _meshProps(meshProps)
    {
        auto& mcMeshProps = *_meshProps.get<MC_MESH_PROPS>();
        // Construct all objective and constraint info here
        auto& mc = mcMeshProps.mesh();

        // 0. VARIABLES

        // One integer length variable per arc
        _arc2idx.clear();
        for (EH arc : mc.edges())
        {
            int idx = _arc2idx.size();
            _arc2idx[arc] = idx;
        }

        // 1. OBJECTIVE
        _objectiveBase = 0.0;
        _objectiveLin = Eigen::VectorXd(_arc2idx.size());
        _objectiveLin.setZero();
        _hessianDense = vector<vector<double>>(_arc2idx.size(), vector<double>(_arc2idx.size(), 0.0));
        for (CH block : mc.cells())
        {
            for (UVWDir axis : {UVWDir::NEG_V_NEG_W, UVWDir::NEG_U_NEG_W, UVWDir::NEG_U_NEG_V})
            {
                vector<int> indices;
                double lenSum = 0.0;
                for (EH arc : mcMeshProps.ref<BLOCK_EDGE_ARCS>(block).at(axis))
                {
                    indices.push_back(_arc2idx.at(arc));
                    lenSum += mcMeshProps.get<ARC_DBL_LENGTH>(arc);
                }
                for (int i = 0; i < indices.size(); i++)
                {
                    _objectiveLin[indices[i]] += -2 * lenSum * _scaling;
                    _hessianDense[indices[i]][indices[i]] += 2;
                    for (int j = i + 1; j < indices.size(); j++)
                    {
                        _hessianDense[indices[i]][indices[j]] += 2;
                        _hessianDense[indices[j]][indices[i]] += 2;
                    }
                }
                _objectiveBase += _scaling * lenSum * _scaling * lenSum;
            }
        }
        for (EH arc : mc.edges())
        {
            int index = _arc2idx.at(arc);
            double lenSum = mcMeshProps.get<ARC_DBL_LENGTH>(arc);
            _objectiveLin[index] += -2 * lenSum * _scaling * INDIVIDUAL_ARC_FACTOR;
            _hessianDense[index][index] += 2 * INDIVIDUAL_ARC_FACTOR;
            _objectiveBase += _scaling * lenSum * _scaling * lenSum * INDIVIDUAL_ARC_FACTOR;
        }
        vector<Eigen::Triplet<double>> triplets;
        for (int i = 0; i < _arc2idx.size(); i++)
            for (int j = 0; j < _arc2idx.size(); j++)
                if (_hessianDense[i][j] > 0)
                    triplets.push_back({i, j, _hessianDense[i][j]});

        _hessianSparse.resize(_arc2idx.size(), _arc2idx.size());
        _hessianSparse.setFromTriplets(triplets.begin(), triplets.end());

        // 2. CONSTRAINTS

        _constraints.clear();

        {
            Eigen::MatrixXd equalityMatrix(2 * mc.n_logical_faces(), _arc2idx.size());
            equalityMatrix.setZero();
            int globalIdx = 0;
            // Each patches opposite arc lengths must match
            for (FH patch : mc.faces())
            {
                HFH hp = mc.halfface_handle(patch, 0);
                if (mc.is_boundary(hp))
                    hp = mc.opposite_halfface_handle(hp);

                auto dir2orderedHas = MCMeshNavigator(_meshProps).halfpatchHalfarcsByDir(hp);
                assert(dir2orderedHas.size() == 4);
                UVWDir dirHpNormal = MCMeshNavigator(_meshProps).halfpatchNormalDir(hp);

                for (UVWDir dirSide :
                     decompose(~(dirHpNormal | -dirHpNormal), {UVWDir::POS_U, UVWDir::POS_V, UVWDir::POS_W}))
                {
                    UVWDir dirSideOpp = -dirSide;
                    set<EH> esA;
                    for (HEH ha : dir2orderedHas[dirSide])
                        esA.insert(mc.edge_handle(ha));
                    set<EH> esB;
                    for (HEH ha : dir2orderedHas[-dirSide])
                        esB.insert(mc.edge_handle(ha));
                    if (esA == esB)
                        continue;

                    for (HEH ha : dir2orderedHas[dirSide])
                        equalityMatrix(globalIdx, _arc2idx.at(mc.edge_handle(ha))) = 1;

                    for (HEH ha : dir2orderedHas[-dirSide])
                        equalityMatrix(globalIdx, _arc2idx.at(mc.edge_handle(ha))) = -1;

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
                {
                    if (std::abs(image(i, j)) > 1e-6)
                    {
                        vars.push_back(j);
                        coeffs.push_back(std::round(image(i, j)));
                    }
                }
                _constraints.push_back({vars, coeffs, 0.0, true});
            }
        }

        if (_allowNeg)
        {
            // Enforce no block with negative extension along any axis
            for (CH b : mc.cells())
            {
                for (UVWDir dir : {UVWDir::POS_V_POS_W, UVWDir::POS_U_POS_W, UVWDir::POS_U_POS_V})
                {
                    vector<int> vars;
                    vector<double> coeffs;
                    auto& arcs = mcMeshProps.ref<BLOCK_EDGE_ARCS>(b).at(dir);
                    for (EH a : arcs)
                    {
                        vars.push_back(_arc2idx.at(a));
                        coeffs.push_back(1.0);
                    }
                    _constraints.push_back({vars, coeffs, 0.0, false});
                }
            }
        }

        // Enforce critical links of length >= 1
        vector<MCMeshNavigator::CriticalLink> criticalLinks;
        map<EH, int> a2criticalLinkIdx;
        map<VH, vector<int>> n2criticalLinksOut;
        map<VH, vector<int>> n2criticalLinksIn;

        int nSeparationConstraints = 0;
        MCMeshNavigator(_meshProps)
            .getCriticalLinks(criticalLinks, a2criticalLinkIdx, n2criticalLinksOut, n2criticalLinksIn, true);

        for (auto& criticalLink : criticalLinks)
        {
            if (criticalLink.pathHas.empty())
                continue;

            vector<int> vars;
            vector<double> coeffs;
            for (HEH ha : criticalLink.pathHas)
            {
                vars.push_back(_arc2idx.at(mc.edge_handle(ha)));
                coeffs.push_back(1.0);
            }

            if (criticalLink.nFrom == criticalLink.nTo)
                _constraints.push_back({vars, coeffs, 3.0, false});
            else
                _constraints.push_back({vars, coeffs, 1.0, false});
            nSeparationConstraints++;
        }

        // Additional code to avoid selfadjacency
        {
            for (CH b : mc.cells())
            {
                auto& dir2ps = mcMeshProps.ref<BLOCK_FACE_PATCHES>(b);
                for (auto& kv : dir2ps)
                {
                    if (isNeg(kv.first))
                        continue;

                    UVWDir dir = kv.first;
                    auto& ps = kv.second;
                    auto& psOpp = dir2ps.at(-dir);
                    set<CH> bs, bsOpp;
                    for (FH p : ps)
                    {
                        HFH hp = mc.halfface_handle(p, 0);
                        if (mc.incident_cell(hp) == b)
                            hp = mc.opposite_halfface_handle(hp);
                        CH bNext = mc.incident_cell(hp);
                        if (!bNext.is_valid())
                            continue;
                        bs.insert(bNext);
                    }
                    for (FH p : psOpp)
                    {
                        HFH hp = mc.halfface_handle(p, 0);
                        if (mc.incident_cell(hp) == b)
                            hp = mc.opposite_halfface_handle(hp);
                        CH bNext = mc.incident_cell(hp);
                        if (!bNext.is_valid())
                            continue;
                        bsOpp.insert(bNext);
                    }
                    bool mustBeNonZero = bs.count(b) != 0 || bsOpp.count(b) != 0 || containsSomeOf(bs, bsOpp);
                    if (mustBeNonZero)
                    {
                        // Constrain non-zero length along dir
                        vector<int> vars;
                        vector<double> coeffs;
                        auto& dir2as = mcMeshProps.ref<BLOCK_EDGE_ARCS>(b);
                        for (EH a : dir2as.at(decompose(~(dir | -dir), DIM_2_DIRS)[0]))
                        {
                            vars.push_back(_arc2idx.at(a));
                            coeffs.push_back(1.0);
                        }
                        _constraints.push_back({vars, coeffs, 1.0, false});
                        nSeparationConstraints++;

                        continue;
                    }
                }
            }
        }
    }

    void addConstraint(const Constraint& constraint)
    {
        _constraints.push_back(constraint);
    }

    /// virtual destructor.
    virtual ~QuantizationProblem()
    {
    }

    /** Copy constructor.*/
    QuantizationProblem(const QuantizationProblem& other)
        : _printSol(other._printSol), _scaling(other._scaling), _allowZeros(other._allowZeros),
          _allowNeg(other._allowNeg), n_nonzeros_jac(other.n_nonzeros_jac), n_nonzeros_hess(other.n_nonzeros_hess),
          _meshProps(other._meshProps), _hessianDense(other._hessianDense), _objectiveLin(other._objectiveLin),
          _hessianSparse(other._hessianSparse), _objectiveBase(other._objectiveBase), _constraints(other._constraints),
          _arc2idx(other._arc2idx)
    {
    }

    /** Assignment operator. no data = nothing to assign*/
    // QuantizationProblem& operator=(const QuantizationProblem&) {}

    /** \name Overloaded functions specific to a TMINLP.*/
    //@{
    /** Pass the type of the variables (INTEGER, BINARY, CONTINUOUS) to the optimizer.
       \param n size of var_types (has to be equal to the number of variables in the problem)
    \param var_types types of the variables (has to be filled by function).
    */
    virtual bool get_variables_types(Ipopt::Index n, VariableType* var_types)
    {
        for (int i = 0; i < n; i++)
            var_types[i] = INTEGER;
        return true;
    }

    /** Pass info about linear and nonlinear variables.*/
    virtual bool get_variables_linearity(Ipopt::Index n, Ipopt::TNLP::LinearityType* var_types)
    {
        for (int i = 0; i < n; i++)
            var_types[i] = Ipopt::TNLP::NON_LINEAR;
        return true;
    }

    /** Pass the type of the constraints (LINEAR, NON_LINEAR) to the optimizer.
    \param m size of const_types (has to be equal to the number of constraints in the problem)
    \param const_types types of the constraints (has to be filled by function).
    */
    virtual bool get_constraints_linearity(Ipopt::Index m, Ipopt::TNLP::LinearityType* const_types)
    {
        for (int i = 0; i < m; i++)
            const_types[i] = Ipopt::TNLP::LINEAR;
        return true;
    }
    //@}

    /** \name Overloaded functions defining a TNLP.
     * This group of function implement the various elements needed to define and solve a TNLP.
     * They are the same as those in a standard Ipopt NLP problem*/
    //@{
    /** Method to pass the main dimensions of the problem to Ipopt.
          \param n number of variables in problem.
          \param m number of constraints.
          \param nnz_jac_g number of nonzeroes in Jacobian of constraints system.
          \param nnz_h_lag number of nonzeroes in Hessian of the Lagrangean.
          \param index_style indicate wether arrays are numbered from 0 (C-style) or
          from 1 (Fortran).
          \return true in case of success.*/
    virtual bool get_nlp_info(Ipopt::Index& n,
                              Ipopt::Index& m,
                              Ipopt::Index& nnz_jac_g,
                              Ipopt::Index& nnz_h_lag,
                              TNLP::IndexStyleEnum& index_style)
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

        index_style = TNLP::C_STYLE;
        return true;
    }

    /** Method to pass the bounds on variables and constraints to Ipopt.
         \param n size of x_l and x_u (has to be equal to the number of variables in the problem)
         \param x_l lower bounds on variables (function should fill it).
         \param x_u upper bounds on the variables (function should fill it).
         \param m size of g_l and g_u (has to be equal to the number of constraints in the problem).
         \param g_l lower bounds of the constraints (function should fill it).
         \param g_u upper bounds of the constraints (function should fill it).
    \return true in case of success.*/
    virtual bool get_bounds_info(Ipopt::Index n, Number* x_l, Number* x_u, Ipopt::Index m, Number* g_l, Number* g_u)
    {
        for (int i = 0; i < n; i++)
        {
            x_l[i] = (_allowNeg ? -DBL_MAX : (_allowZeros ? 0 : 1));
            x_u[i] = DBL_MAX;
        }
        for (int i = 0; i < _constraints.size(); i++)
        {
            if (_constraints[i].equality)
            {
                g_l[i] = g_u[i] = _constraints[i].rhs;
            }
            else
            {
                g_l[i] = _constraints[i].rhs;
                g_u[i] = DBL_MAX;
            }
        }
        return true;
    }

    /** Method to to pass the starting point for optimization to Ipopt.
      \param init_x do we initialize primals?
      \param x pass starting primal points (function should fill it if init_x is 1).
      \param m size of lambda (has to be equal to the number of constraints in the problem).
      \param init_lambda do we initialize duals of constraints?
      \param lambda lower bounds of the constraints (function should fill it).
      \return true in case of success.*/
    virtual bool get_starting_point(Ipopt::Index n,
                                    bool init_x,
                                    Number* x,
                                    bool init_z,
                                    Number* z_L,
                                    Number* z_U,
                                    Ipopt::Index m,
                                    bool init_lambda,
                                    Number* lambda)
    {
        auto& mcMeshProps = *_meshProps.get<MC_MESH_PROPS>();
        if (mcMeshProps.isAllocated<ARC_INT_LENGTH>())
            for (EH arc : mcMeshProps.mesh().edges())
                x[_arc2idx.at(arc)] = mcMeshProps.get<ARC_INT_LENGTH>(arc);
        else
            for (EH arc : mcMeshProps.mesh().edges())
                x[_arc2idx.at(arc)] = _scaling * mcMeshProps.get<ARC_DBL_LENGTH>(arc);

        return true;
    }

    /** Method which compute the value of the objective function at point x.
      \param n size of array x (has to be the number of variables in the problem).
      \param x point where to evaluate.
      \param new_x Is this the first time we evaluate functions at this point?
      (in the present context we don't care).
      \param obj_value value of objective in x (has to be computed by the function).
      \return true in case of success.*/
    virtual bool eval_f(Ipopt::Index n, const Number* x, bool new_x, Number& obj_value)
    {
        Eigen::Map<const Eigen::VectorXd> xd(x, n);
        obj_value = 0.5 * xd.dot(_hessianSparse * xd) + _objectiveLin.dot(xd) + _objectiveBase;

        // std::cout << "Eval f input " << xd.transpose() << std::endl;
        // std::cout << "Eval f objectivelin " << _objectiveLin.transpose() << std::endl;
        // std::cout << "Eval f objectivequad " << Eigen::MatrixXd(_hessianSparse) << std::endl;
        // std::cout << "Eval f returned " << obj_value << std::endl;
        // std::cout << "quad term: " << 0.5 * xd.dot(_hessianSparse * xd) << std::endl;
        // std::cout << "lin term: " << _objectiveLin.dot(xd) << std::endl;
        // std::cout << "base term: " << _objectiveBase << std::endl;
        return true;
    }

    /** Method which compute the gradient of the objective at a point x.
      \param n size of array x (has to be the number of variables in the problem).
      \param x point where to evaluate.
      \param new_x Is this the first time we evaluate functions at this point?
      (in the present context we don't care).
      \param grad_f gradient of objective taken in x (function has to fill it).
      \return true in case of success.*/
    virtual bool eval_grad_f(Ipopt::Index n, const Number* x, bool new_x, Number* grad_f)
    {
        Eigen::Map<const Eigen::VectorXd> xd(x, n);
        Eigen::Map<Eigen::VectorXd> g(grad_f, n);
        g = _hessianSparse * xd + _objectiveLin;
        return true;
    }

    /** Method which compute the value of the functions defining the constraints at a point
      x.
      \param n size of array x (has to be the number of variables in the problem).
      \param x point where to evaluate.
      \param new_x Is this the first time we evaluate functions at this point?
      (in the present context we don't care).
      \param m size of array g (has to be equal to the number of constraints in the problem)
      \param g values of the constraints (function has to fill it).
      \return true in case of success.*/
    virtual bool eval_g(Ipopt::Index n, const Number* x, bool new_x, Ipopt::Index m, Number* g)
    {
        for (int i = 0; i < _constraints.size(); i++)
        {
            g[i] = 0.0;
            for (int j = 0; j < _constraints[i].vars.size(); j++)
                g[i] += x[_constraints[i].vars[j]] * _constraints[i].coeffs[j];
        }

        return true;
    }

    /** Method to compute the Jacobian of the functions defining the constraints.
      If the parameter values==NULL fill the arrays iCol and jRow which store the position of
      the non-zero element of the Jacobian.
      If the paramenter values!=NULL fill values with the non-zero elements of the Jacobian.
      \param n size of array x (has to be the number of variables in the problem).
      \param x point where to evaluate.
      \param new_x Is this the first time we evaluate functions at this point?
      (in the present context we don't care).
      \return true in case of success.*/
    virtual bool eval_jac_g(Ipopt::Index n,
                            const Number* x,
                            bool new_x,
                            Ipopt::Index m,
                            Ipopt::Index nele_jac,
                            Ipopt::Index* iRow,
                            Ipopt::Index* jCol,
                            Number* values)
    {
        if (values == nullptr)
        {
            int gi = 0;
            for (int i = 0; i < _constraints.size(); i++)
            {
                for (int j = 0; j < _constraints[i].vars.size(); j++)
                {
                    iRow[gi] = i;
                    jCol[gi] = _constraints[i].vars[j];
                    gi++;
                }
            }
        }
        else
        {
            int gi = 0;
            for (int i = 0; i < _constraints.size(); i++)
            {
                for (int j = 0; j < _constraints[i].vars.size(); j++)
                {
                    values[gi] = _constraints[i].coeffs[j];
                    gi++;
                }
            }
        }
        return true;
    }

    /** Method to compute the Jacobian of the functions defining the constraints.
      If the parameter values==NULL fill the arrays iCol and jRow which store the position of
      the non-zero element of the Jacobian.
      If the paramenter values!=NULL fill values with the non-zero elements of the Jacobian.
      \param n size of array x (has to be the number of variables in the problem).
      \param x point where to evaluate.
      \param new_x Is this the first time we evaluate functions at this point?
      (in the present context we don't care).
      \param m size of array g (has to be equal to the number of constraints in the problem)
      \param grad_f values of the constraints (function has to fill it).
      \return true in case of success.*/
    virtual bool eval_h(Ipopt::Index n,
                        const Number* x,
                        bool new_x,
                        Number obj_factor,
                        Ipopt::Index m,
                        const Number* lambda,
                        bool new_lambda,
                        Ipopt::Index nele_hess,
                        Ipopt::Index* iRow,
                        Ipopt::Index* jCol,
                        Number* values)
    {
        if (values == NULL)
        {
            int gi = 0;
            for (int i = 0; i < _hessianSparse.outerSize(); ++i)
                for (Eigen::SparseMatrix<double>::InnerIterator it(_hessianSparse, i); it; ++it)
                {
                    if (it.row() >= it.col())
                    {
                        iRow[gi] = it.row();
                        jCol[gi] = it.col();
                        ++gi;
                    }
                }
        }
        else
        {
            int gi = 0;
            for (int i = 0; i < _hessianSparse.outerSize(); ++i)
                for (Eigen::SparseMatrix<double>::InnerIterator it(_hessianSparse, i); it; ++it)
                {
                    if (it.row() >= it.col())
                    {
                        values[gi] = obj_factor * it.value();
                        ++gi;
                    }
                }
        }
        return true;
    }

    /** Method called by Ipopt at the end of optimization.*/
    virtual void finalize_solution(TMINLP::SolverReturn status, Ipopt::Index n, const Number* x, Number obj_value)
    {
        // Maybe print some interesting stuff here?
        DLOG(INFO) << "Problem status: " << status;
        DLOG(INFO) << "Objective value: " << obj_value;
    }

    //@}

    virtual const SosInfo* sosConstraints() const
    {
        return NULL;
    }
    virtual const BranchingInfo* branchingInfo() const
    {
        return NULL;
    }

    void printSolutionAtEndOfAlgorithm()
    {
        _printSol = true;
    }
};
} // namespace impl
#endif

namespace qgp3d
{

bool MCQuantizer::GreaterPathLengthCompare::operator()(const WeaklyMonotonousPath& p1,
                                                       const WeaklyMonotonousPath& p2) const
{
    return p1.length > p2.length || (p1.length == p2.length && p1.path.size() > p2.path.size());
}

MCQuantizer::MCQuantizer(TetMeshProps& meshProps)
    : TetMeshNavigator(meshProps), TetMeshManipulator(meshProps), MCMeshNavigator(meshProps),
      MCMeshManipulator(meshProps)
{
}

MCQuantizer::RetCode MCQuantizer::quantizeArcLengths(double scaling, bool allowZero, bool allowNegative)
{
    const MCMesh& mc = mcMeshProps().mesh();
    const TetMesh& mesh = meshProps().mesh();

    bool wasAllocated = mcMeshProps().isAllocated<ARC_DBL_LENGTH>();
    if (!wasAllocated)
        mcMeshProps().allocate<ARC_DBL_LENGTH>(0.0);
    if (!wasAllocated)
    {
        for (EH arc : mc.edges())
        {
            // Determine current arc length
            double length = 0.0;
            for (HEH he : mcMeshProps().ref<ARC_MESH_HALFEDGES>(arc))
                length += edgeLengthUVW<CHART>(mesh.edge_handle(he));
            mcMeshProps().set<ARC_DBL_LENGTH>(arc, length);
        }
    }

    vector<CriticalLink> criticalLinks;
    map<EH, int> a2criticalLinkIdx;
    map<VH, vector<int>> n2criticalLinksOut;
    map<VH, vector<int>> n2criticalLinksIn;
    getCriticalLinks(criticalLinks, a2criticalLinkIdx, n2criticalLinksOut, n2criticalLinksIn, true);

    try
    {
#ifdef QGP3D_WITH_GUROBI
        GRBEnv env = GRBEnv(true);
#ifdef NDEBUG
        env.set(GRB_IntParam_LogToConsole, false);
#endif
        env.set(GRB_DoubleParam_TimeLimit, TIME_LIMIT);
        env.set(GRB_IntParam_Threads, 1);
        env.start();

        GRBModel model = GRBModel(env);

        // One integer length variable per arc
        std::map<EH, GRBVar> arc2var;
        for (EH arc : mc.edges())
        {
            double length = mcMeshProps().get<ARC_DBL_LENGTH>(arc);
            DLOG(INFO) << "Arc " << arc << " continuous length: " << length;

            // Configure gurobi var
            GRBVar x = model.addVar(allowZero ? (allowNegative ? -GRB_INFINITY : 0.0) : 0.9,
                                    GRB_INFINITY,
                                    0.,
                                    GRB_INTEGER,
                                    std::string("Arc ") + std::to_string(arc.idx()));
            x.set(GRB_DoubleAttr_Start, length);
            arc2var[arc] = {x};
        }

        // Objective Function minimizes deviation from scaled seamless param
        GRBQuadExpr objective = 0.;
        for (CH block : mc.cells())
        {
            // For each axis U/V/W
            for (UVWDir axis : {UVWDir::NEG_V_NEG_W, UVWDir::NEG_U_NEG_W, UVWDir::NEG_U_NEG_V})
            {
                GRBLinExpr varSum = 0.;
                double lenSum = 0;
                // WARNING this way of using range based for loop (chaining 2 function calls)
                // depends on the first calls not returning rvalues!
                for (EH arc : mcMeshProps().ref<BLOCK_EDGE_ARCS>(block).at(axis))
                {
                    varSum += arc2var.at(arc);
                    lenSum += mcMeshProps().get<ARC_DBL_LENGTH>(arc);
                }
                objective += (varSum - scaling * lenSum) * (varSum - scaling * lenSum);
            }
        }
        for (EH arc : mc.edges())
        {
            GRBVar var = arc2var.at(arc);
            double len = mcMeshProps().get<ARC_DBL_LENGTH>(arc);
            objective += (var - scaling * len) * (var - scaling * len) * INDIVIDUAL_ARC_FACTOR;
        }

        model.setObjective(objective, GRB_MINIMIZE);

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
                UVWDir dirSideOpp = -dirSide;

                GRBLinExpr side = 0;
                for (HEH ha : dir2orderedHas[dirSide])
                    side += arc2var.at(mc.edge_handle(ha));

                GRBLinExpr sideOpp = 0;
                for (HEH ha : dir2orderedHas[-dirSide])
                    sideOpp += arc2var.at(mc.edge_handle(ha));

                std::string constraintName = "Patch" + std::to_string(patch.idx()) + "dir"
                                             + std::to_string(static_cast<uint8_t>(dirSide | dirSideOpp));
                model.addConstr(side, GRB_EQUAL, sideOpp, constraintName);
            }
        }

        if (allowNegative)
        {
            // Enforce no block with negative extension along any axis
            for (CH b : mc.cells())
            {
                for (UVWDir dir : {UVWDir::POS_V_POS_W, UVWDir::POS_U_POS_W, UVWDir::POS_U_POS_V})
                {
                    GRBLinExpr sumExpr = 0;
                    auto& arcs = mcMeshProps().ref<BLOCK_EDGE_ARCS>(b).at(dir);
                    for (EH a : arcs)
                        sumExpr += arc2var.at(a);
                    std::string constraintName = "Block" + std::to_string(b.idx()) + "dir"
                                                 + std::to_string(static_cast<uint8_t>(~(dir | -dir)));
                    model.addConstr(sumExpr, GRB_GREATER_EQUAL, 0.0, constraintName);
                }
            }
        }

        // Enforce critical links of length >= 1

        int nSeparationConstraints = 0;

        for (auto& criticalLink : criticalLinks)
        {
            if (criticalLink.pathHas.empty())
                continue;
            GRBLinExpr sumExpr = 0;
            for (HEH ha : criticalLink.pathHas)
                sumExpr += arc2var.at(mc.edge_handle(ha));

            std::string constraintName = "Separation" + std::to_string(nSeparationConstraints);
            if (criticalLink.nFrom == criticalLink.nTo)
                model.addConstr(sumExpr, GRB_GREATER_EQUAL, 3.0, constraintName);
            else
                model.addConstr(sumExpr, GRB_GREATER_EQUAL, 1.0, constraintName);
            nSeparationConstraints++;
        }

        // Additional code to avoid selfadjacency
        {
            for (CH b : mc.cells())
            {
                auto& dir2ps = mcMeshProps().ref<BLOCK_FACE_PATCHES>(b);
                for (auto& kv : dir2ps)
                {
                    if (isNeg(kv.first))
                        continue;

                    UVWDir dir = kv.first;
                    auto& ps = kv.second;
                    auto& psOpp = dir2ps.at(-dir);
                    set<CH> bs, bsOpp;
                    for (FH p : ps)
                    {
                        HFH hp = mc.halfface_handle(p, 0);
                        if (mc.incident_cell(hp) == b)
                            hp = mc.opposite_halfface_handle(hp);
                        CH bNext = mc.incident_cell(hp);
                        if (!bNext.is_valid())
                            continue;
                        bs.insert(bNext);
                    }
                    for (FH p : psOpp)
                    {
                        HFH hp = mc.halfface_handle(p, 0);
                        if (mc.incident_cell(hp) == b)
                            hp = mc.opposite_halfface_handle(hp);
                        CH bNext = mc.incident_cell(hp);
                        if (!bNext.is_valid())
                            continue;
                        bsOpp.insert(bNext);
                    }
                    bool mustBeNonZero = bs.count(b) != 0 || bsOpp.count(b) != 0 || containsSomeOf(bs, bsOpp);
                    if (mustBeNonZero)
                    {
                        // Constrain non-zero length along dir
                        GRBLinExpr sumExpr = 0;
                        auto& dir2as = mcMeshProps().ref<BLOCK_EDGE_ARCS>(b);
                        for (EH a : dir2as.at(decompose(~(dir | -dir), DIM_2_DIRS)[0]))
                            sumExpr += arc2var.at(a);
                        std::string constraintName = "Selfadjacency" + std::to_string(nSeparationConstraints);
                        model.addConstr(sumExpr, GRB_GREATER_EQUAL, 1.0, constraintName);
                        nSeparationConstraints++;

                        continue;
                    }
                }
            }
        }
#else
        Ipopt::SmartPtr<impl::QuantizationProblem> iqp
            = new impl::QuantizationProblem(meshProps(), scaling, allowZero, allowNegative);

        Bonmin::BonminSetup bonmin;
        bonmin.initializeOptionsAndJournalist();
        bonmin.options()->SetNumericValue("bonmin.time_limit", TIME_LIMIT); // changes bonmin's time limit

        // Problem properties
        bonmin.options()->SetBoolValue("hessian_constant", true);
        bonmin.options()->SetBoolValue("jac_c_constant", true);
        bonmin.options()->SetBoolValue("jac_d_constant", true);
        bonmin.options()->SetBoolValue("expect_infeasible_problem", false);
        bonmin.options()->SetStringValue("expect_infeasible_problem", "no");

        // Finetuning
        bonmin.options()->SetStringValue("algorithm", "B-BB");
        bonmin.options()->SetStringValue("mu_oracle", "loqo");
        bonmin.options()->SetStringValue("nlp_failure_behavior", "fathom"); // default
        bonmin.options()->SetStringValue("node_comparison", "best-bound");  // default
        bonmin.options()->SetStringValue("variable_selection",
                                         "qp-strong-branching");   // default: strong-branching
        bonmin.options()->SetStringValue("nlp_solver", "Ipopt");   // default
        bonmin.options()->SetStringValue("warm_start", "optimum"); // default: none

        // Tolerances
        // bonmin.options()->SetNumericValue("allowable_fraction_gap", 1e-3); // obj considered optimal if gap small
        bonmin.options()->SetNumericValue("integer_tolerance", 1e-3);

#ifdef NDEBUG
        bonmin.options()->SetIntegerValue("lp_log_level", 0);
        bonmin.options()->SetIntegerValue("nlp_log_level", 0);
        bonmin.options()->SetIntegerValue("bb_log_level", 0);
        bonmin.options()->SetIntegerValue("print_level", 0);
#endif
        bonmin.initialize(iqp);
#ifdef NDEBUG
        bonmin.continuousSolver()->messageHandler()->setLogLevel(0);
        bonmin.nonlinearSolver()->messageHandler()->setLogLevel(0);
#endif
#endif

        int iter = 0;
        double bestObj = 0.0;
        // Iteratively enforce violated separation constraints
        bool validSolution = false;
        while (!validSolution)
        {
#ifdef QGP3D_WITH_GUROBI
            DLOG(INFO) << "Gurobi solving with " << nSeparationConstraints << " nonzero sum constraints";
            model.optimize();
            auto obj = model.getObjective();
            DLOG(INFO) << "Solved with final objective value of " << obj.getValue();
            iter++;

            int status = model.get(GRB_IntAttr_Status);
            if (status != 2 && status != 9)
            {
                LOG(ERROR) << "Bad status return by GUROBI solver";
                return SOLVER_ERROR;
            }
            bestObj = model.getObjective().getValue();
            mcMeshProps().allocate<ARC_INT_LENGTH>();
            for (EH arc : mc.edges())
            {
                auto& grbvar = arc2var.at(arc);
                mcMeshProps().set<ARC_INT_LENGTH>(arc, (int)std::round(grbvar.get(GRB_DoubleAttr_X)));
            }
#else

            Bonmin::Bab bb;
            bb(bonmin);

            bestObj = bb.bestObj();
            DLOG(INFO) << "Solved with final objective value of " << bestObj;
            iter++;

            int status = bb.mipStatus();
            if (status != Bonmin::Bab::FeasibleOptimal && status != Bonmin::Bab::Feasible)
            {
                LOG(ERROR) << "Bad status return by Bonmin solver";
                return SOLVER_ERROR;
            }
            auto solution = bb.bestSolution();
            mcMeshProps().allocate<ARC_INT_LENGTH>();
            for (EH arc : mc.edges())
                mcMeshProps().set<ARC_INT_LENGTH>(arc, (int)std::round(solution[iqp->_arc2idx.at(arc)]));
#endif

            int minArcLength = 10000;
            for (EH a : mc.edges())
                minArcLength = std::min(minArcLength, mcMeshProps().get<ARC_INT_LENGTH>(a));
            DLOG(INFO) << "Min quantized arc length " << minArcLength;

            for (auto& criticalLink : criticalLinks)
            {
                criticalLink.length = 0;
                for (HEH ha : criticalLink.pathHas)
                    criticalLink.length += mcMeshProps().get<ARC_INT_LENGTH>(mcMeshProps().mesh().edge_handle(ha));
            }

            if (!allowZero)
                validSolution = true;
            else
            {
                vector<bool> isCriticalArc(mc.n_edges(), false);
                vector<bool> isCriticalNode(mc.n_vertices(), false);
                vector<bool> isCriticalPatch(mc.n_faces(), false);
                for (auto& kv : a2criticalLinkIdx)
                    isCriticalArc[kv.first.idx()] = true;
                for (VH n : mc.vertices())
                {
                    auto type = mcMeshProps().nodeType(n);
                    if (type.first == SingularNodeType::SINGULAR || type.second == FeatureNodeType::FEATURE
                        || type.second == FeatureNodeType::SEMI_FEATURE_SINGULAR_BRANCH)
                        isCriticalNode[n.idx()] = true;
                }
                for (FH p : mc.faces())
                    isCriticalPatch[p.idx()]
                        = mc.is_boundary(p)
                          || (mcMeshProps().isAllocated<IS_FEATURE_F>() && mcMeshProps().get<IS_FEATURE_F>(p));

                vector<vector<pair<int, EH>>> forcedNonZeroSum;
                findSeparationViolatingPaths(
                    criticalLinks, isCriticalArc, isCriticalNode, isCriticalPatch, forcedNonZeroSum);
                validSolution = forcedNonZeroSum.size() == 0;
                for (auto& aColl : forcedNonZeroSum)
                {
                    assert(!aColl.empty());
#ifdef QGP3D_WITH_GUROBI
                    vector<GRBVar> vars;
                    vector<double> coeff;
                    for (auto sign2a : aColl)
                    {
                        vars.push_back(arc2var.at(sign2a.second));
                        coeff.push_back(sign2a.first);
                    }
                    GRBLinExpr sumExpr = 0;
                    sumExpr.addTerms(&coeff[0], &vars[0], vars.size());

                    std::string constraintName = "Separation" + std::to_string(nSeparationConstraints);
                    model.addConstr(sumExpr, GRB_GREATER_EQUAL, 0.9, constraintName);
                    nSeparationConstraints++;
#else
                    vector<int> vars;
                    vector<double> coeffs;
                    for (auto sign2a : aColl)
                    {
                        vars.push_back(iqp->_arc2idx.at(sign2a.second));
                        coeffs.push_back(sign2a.first);
                    }
                    iqp->addConstraint(vars, coeffs, 1.0, false);
                    bonmin.nonlinearSolver()->setModel(iqp);
#endif
                }
                if (!validSolution)
                {
                    int nHexes = numHexesInQuantization();
                    DLOG(INFO) << "Invalid intermediate solution has " << nHexes << " hexes";
                }
            }
        }
        int minArcLength = 10000;
        for (EH a : mc.edges())
            minArcLength = std::min(minArcLength, mcMeshProps().get<ARC_INT_LENGTH>(a));
        int nHexes = numHexesInQuantization();

        LOG(INFO) << "Quantized MC, nHexes: " << nHexes << ", objective: " << bestObj << ", iterations: " << iter
                  << ", minArcLength " << minArcLength;
#ifdef QGP3D_WITH_GUROBI
    }
    catch (GRBException e)
    {
        LOG(ERROR) << "Gurobi exception, errcode: " << e.getErrorCode();
        LOG(ERROR) << "Gurobi error message: " << e.getMessage();
        return SOLVER_ERROR;
    }
#else
    }
    catch (...)
    {
        LOG(ERROR) << "Bonmin exception";
        return SOLVER_ERROR;
    }
#endif

    return SUCCESS;
}

int MCQuantizer::numHexesInQuantization() const
{
    const MCMesh& mc = mcMeshProps().mesh();
    int nHexes = 0;
    for (CH b : mc.cells())
    {
        int nBlockHexes = 1;
        for (UVWDir dir : {UVWDir::NEG_U_NEG_V, UVWDir::NEG_U_NEG_W, UVWDir::NEG_V_NEG_W})
        {
            int arcLen = 0;
            // WARNING this way of using range based for loop (chaining 2 function calls)
            // depends on the first calls not returning rvalues!
            for (EH a : mcMeshProps().ref<BLOCK_EDGE_ARCS>(b).at(dir))
                arcLen += mcMeshProps().get<ARC_INT_LENGTH>(a);
            nBlockHexes *= arcLen;
        }
        nHexes += nBlockHexes;
    }
    return nHexes;
}

MCQuantizer::RetCode MCQuantizer::findSeparationViolatingPaths(const vector<CriticalLink>& criticalLinks,
                                                               const vector<bool>& arcIsCritical,
                                                               const vector<bool>& nodeIsCritical,
                                                               const vector<bool>& patchIsCritical,
                                                               vector<vector<pair<int, EH>>>& nonZeroSumArcs) const
{
    nonZeroSumArcs.clear();

    // For each critical link s1 find paths connecting s1 to other critical links s2 or surface patches p2
    // Then check for overlaps between these, accumulating the quantized edge lengths along the path as deltas.
    for (auto& criticalLink : criticalLinks)
    {
        auto ret = traceExhaustPaths(criticalLink, arcIsCritical, nodeIsCritical, patchIsCritical, nonZeroSumArcs);
        if (ret != SUCCESS)
            return ret;
    }

    return SUCCESS;
}

MCQuantizer::RetCode MCQuantizer::traceExhaustPaths(const CriticalLink& criticalLink1,
                                                    const vector<bool>& arcIsCritical,
                                                    const vector<bool>& nodeIsCritical,
                                                    const vector<bool>& patchIsCritical,
                                                    vector<vector<pair<int, EH>>>& nonZeroSumArcs) const
{
    const MCMesh& mcMesh = mcMeshProps().mesh();

    HEH criticalStartHa = criticalLink1.pathHas.empty() ? HEH() : criticalLink1.pathHas.front();

    vector<vector<HEH>> dir2has;
    map<HEH, int> haOrth2dir;
    if (criticalStartHa.is_valid())
    {
        // Categorize all vertical has by direction
        for (HFH criticalStartHp : mcMesh.halfedge_halffaces(criticalStartHa))
        {
            dir2has.emplace_back();

            auto& has = dir2has.back();

            CH bRef = mcMesh.incident_cell(criticalStartHp);
            if (!bRef.is_valid())
                bRef = mcMesh.incident_cell(mcMesh.opposite_halfface_handle(criticalStartHp));

            UVWDir dirStart = halfarcDirInBlock(criticalStartHa, bRef);

            HEH haOrth = mcMesh.prev_halfedge_in_halfface(criticalStartHa, criticalStartHp);
            while (halfarcDirInBlock(haOrth, bRef) == dirStart)
                haOrth = mcMesh.prev_halfedge_in_halfface(haOrth, criticalStartHp);
            haOrth = mcMesh.opposite_halfedge_handle(haOrth);
            has.emplace_back(haOrth);

            HEH haCurr = criticalStartHa;
            HFH hpCurr = criticalStartHp;

            do
            {
                bRef = mcMesh.incident_cell(hpCurr);
                if (!bRef.is_valid())
                    bRef = mcMesh.incident_cell(mcMesh.opposite_halfface_handle(hpCurr));

                dirStart = halfarcDirInBlock(haCurr, bRef);

                while (halfarcDirInBlock(mcMesh.next_halfedge_in_halfface(haCurr, hpCurr), bRef) == dirStart)
                    haCurr = mcMesh.next_halfedge_in_halfface(haCurr, hpCurr);
                haOrth = mcMesh.next_halfedge_in_halfface(haCurr, hpCurr);
                if (haOrth == has.front())
                    break;
                has.emplace_back(haOrth);

                HEH haOrthOpp = mcMesh.opposite_halfedge_handle(haOrth);
                hpCurr = findMatching(mcMesh.halfedge_halffaces(haOrthOpp),
                                      [&](const HFH& hpNext)
                                      {
                                          return mcMesh.face_handle(hpNext) != mcMesh.face_handle(hpCurr)
                                                 && std::find(criticalLink1.pathHas.begin(),
                                                              criticalLink1.pathHas.end(),
                                                              mcMesh.next_halfedge_in_halfface(haOrthOpp, hpNext))
                                                        != criticalLink1.pathHas.end();
                                      });
                if (hpCurr.is_valid())
                    haCurr = mcMesh.next_halfedge_in_halfface(haOrthOpp, hpCurr);
            } while (hpCurr.is_valid() && haCurr != criticalStartHa);
        }
        for (int i = 0; i < (int)dir2has.size(); i++)
            for (HEH ha : dir2has[i])
                haOrth2dir[ha] = i;
    }
    else
    {
        for (HEH ha : mcMesh.outgoing_halfedges(criticalLink1.nFrom))
            dir2has.push_back({ha});
        for (int i = 0; i < (int)dir2has.size(); i++)
            for (HEH ha : dir2has[i])
                haOrth2dir[ha] = i;
    }

    VH criticalStartN = criticalLink1.nFrom;
    CH bRefStart
        = criticalStartHa.is_valid() ? *mcMesh.hec_iter(criticalStartHa) : *mcMesh.vc_iter(criticalLink1.nFrom);
    UVWDir dirStartHa = criticalStartHa.is_valid() ? halfarcDirInBlock(criticalStartHa, bRefStart) : UVWDir::NONE;

    WeaklyMonotonousPath pStart;
    pStart.branchedOff = false;
    pStart.length = 0;
    pStart.n = criticalStartN;
    pStart.dirs1 = dirStartHa | -dirStartHa;
    pStart.path = {};
    pStart.monotonousDirs = ~pStart.dirs1; // Do not search along the link but orthogonally
    pStart.walkedDirs = UVWDir::NONE;
    pStart.deltaMin = (isNeg(dirStartHa) ? criticalLink1.length * toVec(dirStartHa) : Vec3i(0, 0, 0));
    pStart.deltaMax = (isNeg(dirStartHa) ? Vec3i(0, 0, 0) : criticalLink1.length * toVec(dirStartHa));
    if (criticalLink1.cyclic && dirStartHa != UVWDir::NONE)
    {
        pStart.deltaMin[toCoord(dirStartHa)] = INT_MIN;
        pStart.deltaMax[toCoord(dirStartHa)] = INT_MAX;
    }
    assert(pStart.deltaMin[0] <= pStart.deltaMax[0] && pStart.deltaMin[1] <= pStart.deltaMax[1]
           && pStart.deltaMin[2] <= pStart.deltaMax[2]);
    pStart.delta = Vec3i(0, 0, 0);
    pStart.bRefCurrent = bRefStart;
    pStart.transCurrent = Transition();

    using PathQueue
        = std::priority_queue<WeaklyMonotonousPath, std::deque<WeaklyMonotonousPath>, GreaterPathLengthCompare>;

    for (int dirIdx = 0; dirIdx < (int)dir2has.size(); dirIdx++)
    {
        vector<bool> nsVisited(mcMesh.n_vertices(), false);
        vector<bool> nsInitialized(mcMesh.n_vertices(), false);
        PathQueue pathQ;
        pathQ.push(pStart);
        while (!pathQ.empty())
        {
            auto pathCurrent = pathQ.top();
            pathQ.pop();

            if ((pathCurrent.branchedOff && nsVisited[pathCurrent.n.idx()])
                || (!pathCurrent.branchedOff && nsInitialized[pathCurrent.n.idx()]))
                continue;
            if (pathCurrent.branchedOff)
                nsVisited[pathCurrent.n.idx()] = true;
            else
                nsInitialized[pathCurrent.n.idx()] = true;

            map<CH, Transition> b2trans;
            map<HEH, vector<CH>> ha2bRef;
            determineNextHalfedges(pathCurrent, b2trans, ha2bRef);

            // Check for separation violations by current path
            if ((pathCurrent.walkedDirs & pathCurrent.monotonousDirs) != UVWDir::NONE)
            {
                UVWDir violationDir = UVWDir::NONE;
                // Check whether startinterval-criticalnode pairs overlap
                if (nodeIsCritical[pathCurrent.n.idx()]
                    && bboxOverlap(pathCurrent.deltaMin, pathCurrent.deltaMax, pathCurrent.delta, pathCurrent.delta))
                    violationDir = pathCurrent.walkedDirs & ~pathCurrent.dirs1;

                if (violationDir == UVWDir::NONE)
                {
                    // Check whether startinterval-criticallink pairs overlap
                    for (const auto& kv : ha2bRef)
                    {
                        HEH ha2 = kv.first;
                        EH a2 = mcMesh.edge_handle(ha2);
                        if (arcIsCritical[a2.idx()])
                        {
                            CH bRef2 = kv.second.front();
                            Transition trans2 = b2trans.at(bRef2);

                            UVWDir dir2 = trans2.invert().rotate(halfarcDirInBlock(ha2, bRef2));

                            UVWDir possibleViolationDir = pathCurrent.walkedDirs & ~(dir2 | -dir2) & ~pathCurrent.dirs1;
                            if ((possibleViolationDir & pathCurrent.monotonousDirs) == UVWDir::NONE)
                                continue;

                            // Measure just overlap with the arc, not the link
                            int length = mcMeshProps().get<ARC_INT_LENGTH>(a2);

                            if (length < 0)
                            {
                                length = -length;
                                dir2 = -dir2;
                            }
                            Vec3i deltaMin2 = pathCurrent.delta + (isNeg(dir2) ? length * toVec(dir2) : Vec3i(0, 0, 0));
                            Vec3i deltaMax2 = pathCurrent.delta + (isNeg(dir2) ? Vec3i(0, 0, 0) : length * toVec(dir2));

                            // Check for overlap
                            if (bboxOverlap(pathCurrent.deltaMin, pathCurrent.deltaMax, deltaMin2, deltaMax2))
                            {
                                assert(possibleViolationDir != UVWDir::NONE);
                                violationDir = possibleViolationDir;
                                break;
                            }
                        }
                    }
                }

                if (violationDir == UVWDir::NONE)
                {
                    // Check for overlaps between startinterval-surfacepatch pairs
                    for (HFH hp2 : mcMesh.vertex_halffaces(pathCurrent.n))
                    {
                        if (patchIsCritical[mcMesh.face_handle(hp2).idx()])
                        {
                            HFH hp2opp = mcMesh.opposite_halfface_handle(hp2);
                            CH bRef2 = mcMesh.incident_cell(hp2opp);

                            if (b2trans.find(bRef2) == b2trans.end())
                                continue;
                            Transition trans2 = b2trans.at(bRef2);

                            UVWDir hpDirs = UVWDir::NONE;
                            for (HEH ha : mcMesh.halfface_halfedges(hp2opp))
                                hpDirs = hpDirs | halfarcDirInBlock(ha, bRef2);

                            UVWDir hpDirsLocal = trans2.invert().rotate(hpDirs);
                            assert(dim(hpDirsLocal) == 2);

                            UVWDir possibleViolationDir = pathCurrent.walkedDirs & ~hpDirsLocal & ~pathCurrent.dirs1;
                            if ((possibleViolationDir & pathCurrent.monotonousDirs) == UVWDir::NONE)
                                continue;

                            if (checkPatchOverlap(pathCurrent, hp2opp, trans2))
                            {
                                assert(possibleViolationDir != UVWDir::NONE);
                                violationDir = possibleViolationDir;
                                break;
                            }
                        }
                    }
                }

                if (violationDir != UVWDir::NONE)
                {
                    vector<EH> posSignArcs;
                    vector<EH> negSignArcs;
                    for (UVWDir dim1dirPos : {UVWDir::POS_U, UVWDir::POS_V, UVWDir::POS_W})
                    {
                        UVWDir dim1dirAny = dim1dirPos | -dim1dirPos;
                        UVWDir intersection = violationDir & dim1dirAny;
                        if (intersection != UVWDir::NONE)
                        {
                            if (intersection == dim1dirAny)
                            {
                                const auto& posArcs = pathCurrent.dir2walkedArcs.at(dim1dirPos);
                                const auto& negArcs = pathCurrent.dir2walkedArcs.at(-dim1dirPos);
                                assert(posArcs.size() > 0);
                                assert(negArcs.size() > 0);
                                double lengthPos = 0;
                                double lengthNeg = 0;
                                for (EH a : posArcs)
                                    lengthPos += mcMeshProps().get<ARC_DBL_LENGTH>(a);
                                for (EH a : negArcs)
                                    lengthNeg += mcMeshProps().get<ARC_DBL_LENGTH>(a);
                                if (lengthPos > lengthNeg)
                                {
                                    posSignArcs.insert(posSignArcs.end(), posArcs.begin(), posArcs.end());
                                    negSignArcs.insert(negSignArcs.end(), negArcs.begin(), negArcs.end());
                                }
                                else
                                {
                                    posSignArcs.insert(posSignArcs.end(), negArcs.begin(), negArcs.end());
                                    negSignArcs.insert(negSignArcs.end(), posArcs.begin(), posArcs.end());
                                }
                            }
                            else
                            {
                                auto& arcs = pathCurrent.dir2walkedArcs.at(intersection);
                                posSignArcs.insert(posSignArcs.end(), arcs.begin(), arcs.end());
                            }
                        }
                    }
                    assert(posSignArcs.size() > 0);
                    if (posSignArcs.empty())
                        continue;
                    nonZeroSumArcs.emplace_back();
                    auto& nonZeroSum = nonZeroSumArcs.back();
                    for (EH a : posSignArcs)
                        nonZeroSum.push_back({1, a});
                    for (EH a : negSignArcs)
                        nonZeroSum.push_back({-1, a});

                    assert(pathCurrent.branchedOff);
                    break;
                }
            }

            for (const auto& kv : ha2bRef)
            {
                HEH ha = kv.first;
                VH nTo = mcMesh.to_vertex_handle(ha);
                if (pathCurrent.branchedOff && nsVisited[nTo.idx()])
                    continue;

                for (CH bRef : kv.second)
                {
                    auto& trans = b2trans.at(bRef);
                    WeaklyMonotonousPath nextP = pathCurrent;
                    nextP.bRefCurrent = bRef;
                    nextP.transCurrent = trans;

                    if (!checkP0Containment(nextP))
                        continue;

                    nextP.n = nTo;

                    UVWDir walkedDir = trans.invert().rotate(halfarcDirInBlock(ha, nextP.bRefCurrent));

                    // Record branch-off from link, allow only limited set of edges
                    if (!nextP.branchedOff)
                    {
                        auto itDir = haOrth2dir.find(ha);
                        if (itDir != haOrth2dir.end())
                        {
                            if (itDir->second != dirIdx)
                                continue;
                            nextP.branchedOff = true;

                            nextP.monotonousDirs = nextP.monotonousDirs & walkedDir;
                        }
                        else if ((walkedDir & nextP.dirs1) == UVWDir::NONE)
                        {
                            // Dont allow branch-off at another arc than those allowed
                            continue;
                        }
                    }
                    if (nextP.branchedOff)
                        nextP.path.emplace_back(ha);

                    nextP.delta += mcMeshProps().get<ARC_INT_LENGTH>(mcMesh.edge_handle(ha)) * toVec(walkedDir);

                    nextP.walkedDirs = nextP.walkedDirs | walkedDir;
                    nextP.monotonousDirs = nextP.monotonousDirs & ~(-walkedDir);

                    // Do not walk unnecessary circles
                    if (nextP.monotonousDirs == UVWDir::NONE)
                        continue;

                    nextP.dir2walkedArcs[walkedDir].emplace_back(mcMesh.edge_handle(ha));

                    // Walk along the link for free, but other arcs accumulate distance
                    // if (nextP.branchedOff)
                    if ((walkedDir & nextP.dirs1) == UVWDir::NONE)
                        nextP.length += mcMeshProps().get<ARC_DBL_LENGTH>(mcMesh.edge_handle(ha));

                    pathQ.push(nextP);
                    break; // Only push the edge once (still need to check each bRef)
                }
            }
        }
    }
    return SUCCESS;
}

bool MCQuantizer::bboxOverlap(const Vec3i& min1, const Vec3i& max1, const Vec3i& min2, const Vec3i& max2)
{
    int coordTouchingOrOverlaps = 0;
    for (int i = 0; i < 3; i++)
    {
        // min2[i] is in [min1[i], max1[i]]
        // OR max2[i] is in [min1[i], max1[i]]
        // OR min1[i] is in [min2[i], max2[i]]
        // OR max1[i] is in [min2[i], max2[i]]
        coordTouchingOrOverlaps
            += (min1[i] <= min2[i] && min2[i] <= max1[i]) || (min1[i] <= max2[i] && max2[i] <= max1[i])
               || (min2[i] <= min1[i] && min1[i] <= max2[i]) || (min2[i] <= max1[i] && max1[i] <= max2[i]);

        if (coordTouchingOrOverlaps == i)
            return false;
    }

    return true;
}

MCQuantizer::RetCode MCQuantizer::determineNextHalfedges(const WeaklyMonotonousPath& pathCurrent,
                                                         map<CH, Transition>& b2trans,
                                                         map<HEH, vector<CH>>& ha2bRef) const
{
    auto& mcMesh = mcMeshProps().mesh();

    b2trans = map<CH, Transition>({{pathCurrent.bRefCurrent, {pathCurrent.transCurrent}}});

    // Floodfill blocks around n, storing current transition for each expanded block
    list<pair<CH, Transition>> bQ({{pathCurrent.bRefCurrent, pathCurrent.transCurrent}});

    while (!bQ.empty())
    {
        auto b2t = bQ.front();
        bQ.pop_front();

        for (HFH hp : mcMesh.cell_halffaces(b2t.first))
        {
            HFH hpOpp = mcMesh.opposite_halfface_handle(hp);
            CH bNext = mcMesh.incident_cell(hpOpp);
            if (!bNext.is_valid() || b2trans.find(bNext) != b2trans.end())
                continue;
            if (!contains(mcMesh.halfface_vertices(hp), pathCurrent.n))
                continue;

            // Check if overlapping with current path
            if (!checkPatchOverlap(pathCurrent, hp, b2t.second))
                continue;

            Transition trans = b2t.second.chain(mcMeshProps().hpTransition<PATCH_TRANSITION>(hp));

            bool exists = b2trans.find(bNext) != b2trans.end();
            if (exists)
                continue;
            b2trans[bNext] = trans;

            bQ.push_back({bNext, trans});
        }
    }
    ha2bRef.clear();
    for (HEH ha : mcMesh.outgoing_halfedges(pathCurrent.n))
    {
        if (!pathCurrent.path.empty() && pathCurrent.path.back() == mcMesh.opposite_halfedge_handle(ha))
            continue;
        for (CH b : mcMesh.halfedge_cells(ha))
        {
            auto it = b2trans.find(b);
            if (it != b2trans.end())
            {
                ha2bRef[ha].emplace_back(b);
            }
        }
    }
    return SUCCESS;
}

bool MCQuantizer::checkArcOverlap(const WeaklyMonotonousPath& pathCurrent,
                                  const HEH& ha,
                                  const CH& bRef,
                                  const Transition& trans) const
{
    auto& mcMesh = mcMeshProps().mesh();

    UVWDir dir2 = trans.invert().rotate(halfarcDirInBlock(ha, bRef));
    int arcLen = mcMeshProps().get<ARC_INT_LENGTH>(mcMesh.edge_handle(ha));

    Vec3i deltaMin2 = isNeg(dir2) ? pathCurrent.delta + arcLen * toVec(dir2) : pathCurrent.delta;
    Vec3i deltaMax2 = isNeg(dir2) ? pathCurrent.delta : pathCurrent.delta + arcLen * toVec(dir2);

    // Check for overlap
    return bboxOverlap(pathCurrent.deltaMin, pathCurrent.deltaMax, deltaMin2, deltaMax2);
}

bool MCQuantizer::checkPatchOverlap(const WeaklyMonotonousPath& pathCurrent,
                                    const HFH& hp,
                                    const Transition& trans) const
{
    auto& mcMesh = mcMeshProps().mesh();

    CH bRef2 = mcMesh.incident_cell(hp);

    Vec3i deltaMin2 = pathCurrent.delta;
    Vec3i deltaMax2 = pathCurrent.delta;

    HEH ha1 = findMatching(mcMesh.halfface_halfedges(hp),
                           [&](const HEH& ha) { return mcMesh.from_vertex_handle(ha) == pathCurrent.n; });
    assert(ha1.is_valid());
    HEH haCurr = ha1;
    Vec3i deltaCurr = pathCurrent.delta;

    Transition transInv = trans.invert();
    do
    {
        UVWDir dir = transInv.rotate(halfarcDirInBlock(haCurr, bRef2));
        EH aCurr = mcMesh.edge_handle(haCurr);
        int lengthHa = mcMeshProps().get<ARC_INT_LENGTH>(aCurr);
        deltaCurr += lengthHa * toVec(dir);
        for (int coord = 0; coord < 3; coord++)
        {
            deltaMax2[coord] = std::max(deltaMax2[coord], deltaCurr[coord]);
            deltaMin2[coord] = std::min(deltaMin2[coord], deltaCurr[coord]);
        }
        haCurr = mcMesh.next_halfedge_in_halfface(haCurr, hp);
    } while (haCurr != ha1);

    assert(deltaCurr == pathCurrent.delta);
    assert(deltaMin2[0] <= deltaMax2[0] && deltaMin2[1] <= deltaMax2[1] && deltaMin2[2] <= deltaMax2[2]);

    return bboxOverlap(pathCurrent.deltaMin, pathCurrent.deltaMax, deltaMin2, deltaMax2);
}

bool MCQuantizer::checkP0Containment(const WeaklyMonotonousPath& pathCurrent) const
{
    auto& mcMesh = mcMeshProps().mesh();

    CH bRef = pathCurrent.bRefCurrent;

    // Find the closest corner of bRef

    // Graph search
    list<pair<VH, Vec3i>> nQ;
    nQ.push_back({pathCurrent.n, Vec3i(0, 0, 0)});
    map<VH, Vec3i> n2displacement({{pathCurrent.n, Vec3i(0, 0, 0)}});
    while (!nQ.empty())
    {
        VH n = nQ.front().first;
        Vec3i displacement = nQ.front().second;
        nQ.pop_front();

        for (HEH ha : mcMesh.outgoing_halfedges(n))
        {
            bool inBRef = false;
            for (CH b : mcMesh.halfedge_cells(ha))
                if (b == bRef)
                {
                    inBRef = true;
                    break;
                }
            if (inBRef)
            {
                VH nTo = mcMesh.to_vertex_handle(ha);
                if (n2displacement.find(nTo) == n2displacement.end())
                {
                    UVWDir dirHa = pathCurrent.transCurrent.invert().rotate(halfarcDirInBlock(ha, bRef));
                    double length = mcMeshProps().get<ARC_INT_LENGTH>(mcMesh.edge_handle(ha));
                    Vec3i newDisplacement = displacement + toVec(dirHa) * length;
                    n2displacement[nTo] = newDisplacement;
                    nQ.push_back({nTo, newDisplacement});
                }
            }
        }
    }

    auto& dir2n = mcMeshProps().ref<BLOCK_CORNER_NODES>(bRef);
    VH nMin = dir2n.at(UVWDir::NEG_U_NEG_V_NEG_W);
    VH nMax = dir2n.at(UVWDir::POS_U_POS_V_POS_W);

    return bboxOverlap(pathCurrent.deltaMin,
                       pathCurrent.deltaMax,
                       pathCurrent.delta + n2displacement.at(nMin),
                       pathCurrent.delta + n2displacement.at(nMax));
}

} // namespace qgp3d
