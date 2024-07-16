#ifndef QGP3D_WITH_GUROBI
#ifndef QGP3D_WITHOUT_IQP

#ifndef QGP3D_BONMINIQPSOLVER_HPP
#define QGP3D_BONMINIQPSOLVER_HPP

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseQR>

#include "BonBonminSetup.hpp"
#include "BonTMINLP.hpp"

#include "QGP3D/IQP/BaseIQPSolver.hpp"

namespace qgp3d
{
namespace impl
{

/**
 * @brief IQP quantization solving based on open-source IQP solver BONMIN
 */
class BonminIQPSolver : public BaseIQPSolver, public virtual MCMeshManipulator
{
  public:
    /**
     * @brief Create new IQP solver instance based on Bonmin
     *
     * @param meshProps IN: mesh whose MC to quantize
     * @param scaling IN: scale target lengths by this factor for quantization
     * @param varLowerBound IN: lower bound for arc lengths
     * @param maxSeconds IN: time limit for solver in seconds
     * @param individualArcFactor IN: objective = this * <arc-length-deviation> + (1-this) * <block-length-deviation>
     */
    BonminIQPSolver(TetMeshProps& meshProps,
                    double scaling,
                    double varLowerBound,
                    double maxSeconds = 180,
                    double individualArcFactor = 1.0);

    /// @brief See \ref BaseIQPSolver for usage

    ///@{
    virtual void setupDefaultOptions();
    virtual void setupVariables();
    virtual void setupObjective();
    virtual void setupConstraints();
    virtual void finalizeSetup();
    virtual void enableQuickSolve();
    virtual void enableExactSolve();
    virtual BaseIQPSolver::RetCode solve();
    virtual double objectiveValue() const;
    virtual void addConstraints(const vector<vector<pair<int, EH>>>& nonZeroSum);
    virtual void useCurrentAsWarmStart();
    ///@}

    virtual ~BonminIQPSolver()
    {
    }

  private:
    /**
     * @brief Convenient encoding of constraints
     */
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

    /**
     * @brief Needed for interfacing with Bonmin
     */
    class BonminIQP : public Bonmin::TMINLP, public virtual MCMeshNavigator
    {
      public:
        double _scaling;
        double _varLowerBound;

        map<EH, int> _arc2idx;

        double _objectiveBase;
        Eigen::VectorXd _objectiveLin;
        Eigen::SparseMatrix<double> _hessianSparse;
        vector<Constraint> _constraints;

        BonminIQP(const TetMeshProps& meshProps, double scaling, double varLowerBound);

        /// virtual destructor.
        virtual ~BonminIQP()
        {
        }

        /** Copy constructor.*/
        BonminIQP(const BonminIQP& other);

        /** \name Overloaded functions specific to a TMINLP.*/
        //@{
        /** Pass the type of the variables (INTEGER, BINARY, CONTINUOUS) to the optimizer.
           \param n size of var_types (has to be equal to the number of variables in the problem)
        \param var_types types of the variables (has to be filled by function).
        */
        virtual bool get_variables_types(Ipopt::Index n, Bonmin::TMINLP::VariableType* var_types);

        /** Pass info about linear and nonlinear variables.*/
        virtual bool get_variables_linearity(Ipopt::Index n, Ipopt::TNLP::LinearityType* var_types);

        /** Pass the type of the constraints (LINEAR, NON_LINEAR) to the optimizer.
        \param m size of const_types (has to be equal to the number of constraints in the problem)
        \param const_types types of the constraints (has to be filled by function).
        */
        virtual bool get_constraints_linearity(Ipopt::Index m, Ipopt::TNLP::LinearityType* const_types);
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
                                  Ipopt::TNLP::IndexStyleEnum& index_style);

        /** Method to pass the bounds on variables and constraints to Ipopt.
             \param n size of x_l and x_u (has to be equal to the number of variables in the problem)
             \param x_l lower bounds on variables (function should fill it).
             \param x_u upper bounds on the variables (function should fill it).
             \param m size of g_l and g_u (has to be equal to the number of constraints in the problem).
             \param g_l lower bounds of the constraints (function should fill it).
             \param g_u upper bounds of the constraints (function should fill it).
        \return true in case of success.*/
        virtual bool get_bounds_info(Ipopt::Index n,
                                     Ipopt::Number* x_l,
                                     Ipopt::Number* x_u,
                                     Ipopt::Index m,
                                     Ipopt::Number* g_l,
                                     Ipopt::Number* g_u);

        /** Method to to pass the starting point for optimization to Ipopt.
          \param init_x do we initialize primals?
          \param x pass starting primal points (function should fill it if init_x is 1).
          \param m size of lambda (has to be equal to the number of constraints in the problem).
          \param init_lambda do we initialize duals of constraints?
          \param lambda lower bounds of the constraints (function should fill it).
          \return true in case of success.*/
        virtual bool get_starting_point(Ipopt::Index n,
                                        bool init_x,
                                        Ipopt::Number* x,
                                        bool init_z,
                                        Ipopt::Number* z_L,
                                        Ipopt::Number* z_U,
                                        Ipopt::Index m,
                                        bool init_lambda,
                                        Ipopt::Number* lambda);

        /** Method which compute the value of the objective function at point x.
          \param n size of array x (has to be the number of variables in the problem).
          \param x point where to evaluate.
          \param new_x Is this the first time we evaluate functions at this point?
          (in the present context we don't care).
          \param obj_value value of objective in x (has to be computed by the function).
          \return true in case of success.*/
        virtual bool eval_f(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Number& obj_value);

        /** Method which compute the gradient of the objective at a point x.
          \param n size of array x (has to be the number of variables in the problem).
          \param x point where to evaluate.
          \param new_x Is this the first time we evaluate functions at this point?
          (in the present context we don't care).
          \param grad_f gradient of objective taken in x (function has to fill it).
          \return true in case of success.*/
        virtual bool eval_grad_f(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Number* grad_f);

        /** Method which compute the value of the functions defining the constraints at a point
          x.
          \param n size of array x (has to be the number of variables in the problem).
          \param x point where to evaluate.
          \param new_x Is this the first time we evaluate functions at this point?
          (in the present context we don't care).
          \param m size of array g (has to be equal to the number of constraints in the problem)
          \param g values of the constraints (function has to fill it).
          \return true in case of success.*/
        virtual bool eval_g(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Index m, Ipopt::Number* g);

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
                                const Ipopt::Number* x,
                                bool new_x,
                                Ipopt::Index m,
                                Ipopt::Index nele_jac,
                                Ipopt::Index* iRow,
                                Ipopt::Index* jCol,
                                Ipopt::Number* values);

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
                            const Ipopt::Number* x,
                            bool new_x,
                            Ipopt::Number obj_factor,
                            Ipopt::Index m,
                            const Ipopt::Number* lambda,
                            bool new_lambda,
                            Ipopt::Index nele_hess,
                            Ipopt::Index* iRow,
                            Ipopt::Index* jCol,
                            Ipopt::Number* values);

        /** Method called by Ipopt at the end of optimization.*/
        virtual void finalize_solution(Bonmin::TMINLP::SolverReturn status,
                                       Ipopt::Index n,
                                       const Ipopt::Number* x,
                                       Ipopt::Number obj_value);

        //@}

        virtual const Bonmin::TMINLP::SosInfo* sosConstraints() const
        {
            return NULL;
        }
        virtual const Bonmin::TMINLP::BranchingInfo* branchingInfo() const
        {
            return NULL;
        }
    };

    Ipopt::SmartPtr<BonminIQP> _instance; // The Bonmin problem instance
    Bonmin::BonminSetup _setupQuick;      // used for solving quick but suboptimal
    Bonmin::BonminSetup _setupExact;      // used for solving exact with tight optimality gap
    double _lastObjective = DBL_MAX;      // objective of last solution
    bool _quickSolve = false;             // which setup to use for solving
};

} // namespace impl
} // namespace qgp3d

#endif
#endif
#endif
