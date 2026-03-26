#ifndef QGP3D_OBJECTIVEFUNCTION_HPP
#define QGP3D_OBJECTIVEFUNCTION_HPP

#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace qgp3d
{
/***
 * Abstract objective function interface, exposing evaluation, gradient and hessian.
 */
class ObjectiveFunction
{
  public:
    ObjectiveFunction()
    {
    }

    virtual double function_value(const Eigen::VectorXd& x) = 0;
    virtual void gradient(const Eigen::VectorXd& x, Eigen::VectorXd& grad) = 0;
    virtual void hessian(const Eigen::VectorXd& x, Eigen::SparseMatrix<double>& hess) = 0;
    virtual bool is_hessian_const() = 0;

    ~ObjectiveFunction() {};
};

/**
 * @brief Quadratic objective function of arbitrary input dimension.
 */
class QuadraticObjective : public ObjectiveFunction
{
  public:
    QuadraticObjective(double f0_, const Eigen::VectorXd& grad0_, const Eigen::SparseMatrix<double>& hess_)
        : f0(f0_), grad0(grad0_), hess(hess_)
    {
    }

    inline virtual double function_value(const Eigen::VectorXd& x) override
    {
        return f0 + grad0.dot(x) + 0.5 * x.dot(hess * x);
    }

    inline virtual void gradient(const Eigen::VectorXd& x, Eigen::VectorXd& grad) override
    {
        grad = grad0 + hess * x;
    }

    inline virtual void hessian(const Eigen::VectorXd& x, Eigen::SparseMatrix<double>& hess_) override
    {
        (void)x;
        hess_ = hess;
    }

    inline bool is_hessian_const()
    {
        return true;
    }

    const double f0;
    const Eigen::VectorXd grad0;
    const Eigen::SparseMatrix<double> hess;
};

} // namespace qgp3d
#endif
