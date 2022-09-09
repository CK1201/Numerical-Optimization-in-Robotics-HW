#include "socp/lbfgs.hpp"
#include <Eigen/Eigen>
#include <iostream>
#include <iomanip>

class SOCPExample
{
public:
    int run()
    {
        N_ = 7;
        m_ = 7;
        rho_ = 1;
        gamma_ = 1;
        beta_ = 1e5;
        int ret, times = 0;
        double finalCost, e_cons = 1e-4, e_prec = 1e-4, res_cons = 1e5, res_prec = 1e5;
        Eigen::VectorXd x(N_), x_opt(N_);
        x_opt << -0.127286, -0.506097, -1.01317, -1.77744, -3.06097, -5.66462, -13.7682;

        /* Set matrix */
        f_ = Eigen::VectorXd::Zero(N_);
        A_ = Eigen::MatrixXd::Zero(m_, N_);
        b_ = Eigen::VectorXd::Zero(m_);
        c_ = Eigen::VectorXd::Zero(N_);
        d_ = 1.0;
        A_ba_ = Eigen::MatrixXd::Zero(m_ + 1, N_);
        b_ba_ = Eigen::VectorXd::Zero(m_ + 1);
        mu_ = Eigen::VectorXd::Zero(m_ + 1);
        for (int i = 0; i < N_; ++i){
            x(i) = 0;
            f_(i) = i + 1;
            A_(i, i) = 7 - i;
            b_(i) = i * 2 + 1;
            A_ba_(i + 1, i) = 7 - i;
            b_ba_(i + 1) = i * 2 + 1;
        }
        c_(0) = 1;
        A_ba_(0, 0) = 1;
        b_ba_(0) = 1;

        /* Set the minimization parameters */
        lbfgs::lbfgs_parameter_t params;
        params.mem_size = 16;
        params.past = 0;
        params.g_epsilon = 1e-8;
        params.min_step = 1e-32;
        params.delta = 1e-4;

        while(res_cons > e_cons || res_prec > e_prec){
            ++times;
            /* Start minimization */
            ret = lbfgs::lbfgs_optimize(x,
                                            finalCost,
                                            costFunction,
                                            nullptr,
                                            monitorProgress,
                                            this,
                                            params);
            Eigen::VectorXd v = mu_ / rho_ - A_ba_ * x - b_ba_;
            Eigen::VectorXd g = f_ - rho_ * A_ba_.transpose() * proj2SOC(v);
            res_prec = g.lpNorm<Eigen::Infinity>();
            
            Eigen::VectorXd proj = mu_ / rho_ - proj2SOC(v);
            res_cons = proj.lpNorm<Eigen::Infinity>();

            // params.g_epsilon = std::max(params.g_epsilon * 0.1, 1e-8);

            /* Update */
            mu_ = proj2SOC(mu_ - rho_ * (A_ba_ * x + b_ba_));
            rho_ = std::min(beta_, (1 + gamma_) * rho_);

            std::cout << "================================" << std::endl
                    << "iter times: " << times << std::endl
                    // << "res_cons: " << res_cons << std::endl
                    // << "res_prec: " << g.transpose() << std::endl
                    << "L-BFGS Optimization Returned: " << ret << std::endl
                    << "\033[32m" << "      Minimized Cost: " << f_.dot(x) << "\033[0m" << std::endl
                    << "given Minimized Cost: " << -156.9589 << std::endl
                    << "\033[32m" << "      Optimal Variables: " << x.transpose() << "\033[0m" << std::endl
                    << "given Optimal Variables: " << x_opt.transpose() << std::endl << std::endl;
        }
        return ret;
    }

private:
    int N_, m_;
    double d_, rho_, gamma_, beta_;
    Eigen::MatrixXd A_, A_ba_;
    Eigen::VectorXd f_, b_, c_, b_ba_, mu_;

    Eigen::VectorXd proj2SOC(Eigen::VectorXd v){
        double v0 = v(0);
        int N = v.size();
        Eigen::VectorXd v1 = v.segment(1, N - 1), proj;
        double v1norm = v1.norm();
        if(v0 <= -v1norm){
            proj = Eigen::VectorXd::Zero(N);
        }else if (v0 >= v1norm){
            proj = v;
        }else{
            double coeff = (v0 + v1norm) / 2 / v1norm;
            proj = v * coeff;
            proj(0) = v1norm * coeff;
        }
        return proj;
    }

    Eigen::MatrixXd grad_proj2SOC(Eigen::VectorXd v){
        double v0 = v(0);
        int N = v.size();
        Eigen::VectorXd v1 = v.segment(1, N - 1);
        Eigen::MatrixXd grad, v1_mat = Eigen::MatrixXd::Zero(N - 1, 1);
        for(int i = 0; i < N -1; ++i){
            v1_mat(i, 0) = v1(i);
        }
        
        double v1norm = v1.norm();
        if(v0 <= -v1norm){
            grad = Eigen::MatrixXd::Zero(N, N);
        }else if (v0 >= v1norm){
            grad = Eigen::MatrixXd::Identity(N, N);
        }else{
            grad = Eigen::MatrixXd::Zero(N, N);
            grad(0, 0) = 0.5;
            grad.block(1, 0, N - 1, 1) = v1_mat / 2 / v1norm;
            grad.block(0, 1, 1, N - 1) = v1_mat.transpose() / 2 / v1norm;
            grad.block(1, 1, N - 1, N - 1) = Eigen::MatrixXd::Identity(N - 1, N - 1) * (v0 + v1norm) / 2 / v1norm - v1_mat * v1_mat.transpose() * v0 / 2 / v1norm / v1norm / v1norm;
        }
        return grad;
    }

    Eigen::MatrixXd Hessian(Eigen::VectorXd x){
        return rho_ * A_ba_.transpose() * grad_proj2SOC(mu_ - rho_ * (A_ba_ * x + b_ba_)) * A_ba_;
        // return hess;
    }

    static double costFunction(void *instance,
                               const Eigen::VectorXd &x,
                               Eigen::VectorXd &g)
    {
        SOCPExample &obj = *(SOCPExample *)instance;
        // const int n = x.size();
        auto v = obj.mu_ / obj.rho_ - obj.A_ba_ * x - obj.b_ba_;
        auto proj = obj.proj2SOC(v);
        double fx = obj.f_.dot(x) + obj.rho_ / 2 * proj.squaredNorm();
        g = obj.f_ - obj.rho_ * obj.A_ba_.transpose() * proj; //  * obj.grad_proj2SOC(v)
        return fx;
    }

    static int monitorProgress(void *instance,
                               const Eigen::VectorXd &x,
                               const Eigen::VectorXd &g,
                               const double fx,
                               const double step,
                               const int k,
                               const int ls)
    {
        // std::cout << std::setprecision(4)
        //           << "================================" << std::endl
        //           << "Iteration: " << k << std::endl
        //           << "Function Value: " << fx << std::endl
        //           << "Gradient Inf Norm: " << g.cwiseAbs().maxCoeff() << std::endl
        //           << "Variables: " << std::endl
        //           << x.transpose() << std::endl;
        return 0;
    }
};

int main(int argc, char **argv)
{
    SOCPExample example;
    return example.run();
}
