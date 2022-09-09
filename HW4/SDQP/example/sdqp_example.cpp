#include <iostream>

#include "sdqp/sdqp.hpp"

using namespace std;
using namespace Eigen;

int main(int argc, char **argv)
{
    const int m = 5, d = 3;
    Eigen::Matrix<double, d, d> Q, Identity3, Q_prox, QT;
    Eigen::Matrix<double, d, 1> c, c_prox;
    Eigen::Matrix<double, d, 1> x, x_ba, x_opt;        // decision variables
    Eigen::Matrix<double, -1, d> A(m, d); // constraint matrix
    Eigen::VectorXd b(m);                 // constraint bound
    x_ba << -10.0, 4.0, 10.0;
    x_opt << -103.0 / 97.0, -93.0 / 97.0, 95.0 / 97.0;
    x = x_ba;
    Identity3 << 1.0, 0.0, 0.0, 
                 0.0, 1.0, 0.0, 
                 0.0, 0.0, 1.0;
    double rho = 1;

    Q << 8.0, -6.0, 2.0, 
         -6.0, 6.0, -3.0, 
         2.0, -3.0, 2.0;
    // QT << 2.0, 0.0, 2.0,
    //       -1.0, 1.0, -2.0,
    //       0.0, -1.0, 1.0;
    // Q = (Q * QT) / 2;
    // auto eigenval = ((Q + QT) / 2).eigenvalues().real();
    // double min_eigenval = 1e5;
    // for(int i = 0; i < d; ++i){
    //     min_eigenval = std::min(min_eigenval, eigenval(i));
    // }
    // min_eigenval = std::min(0.0, min_eigenval);
    

    c << 1.0, 3.0, -2.0;

    A << 0.0, -1.0, -2.0,
        -1.0, 1.0, -3.0,
        1.0, -2.0, 0.0,
        -1.0, -2.0, -1.0,
        3.0, 5.0, 1.0;
    b << -1.0, 2.0, 7.0, 2.0, -1.0;

    double tol = 1;
    int times = 0;
    while (tol > 1e-5){
        ++times;
        Q_prox = Q + Identity3 / rho;
        c_prox = c - x_ba / rho;
        double minobj = sdqp::sdqp<d>(Q_prox, c_prox, A, b, x);
        tol = (x - x_ba).norm() / std::max(1.0, x.norm());

        std::cout << "======================" << std::endl;
        std::cout << "iter times: " << times << std::endl;
        std::cout << "tol: " << tol << std::endl;
        std::cout << "Q_prox eigenvalues: " << Q_prox.eigenvalues().real().transpose() << std::endl;
        std::cout << "\033[32m"
              << "      optimal sol: " << x.transpose() << "\033[0m" << std::endl;
        std::cout << "given optimal sol: " << x_opt.transpose() << std::endl;
        std::cout << "\033[32m"
              << "      optimal obj: " << 0.5 * (Q * x).dot(x) + c.dot(x) << "\033[0m" << std::endl;
        std::cout << "given optimal obj: " << -295.0 / 97.0 << std::endl;
        std::cout << "cons precision: " << (A * x - b).transpose() << std::endl << std::endl;
        rho = std::min(rho * 10, 1e6);
        x_ba = x;
    }

    

    return 0;
}
