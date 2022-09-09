#include <iostream>

#include "sdqp/sdqp.hpp"

using namespace std;
using namespace Eigen;

int main(int argc, char **argv)
{
    int m = 10;
    Eigen::Matrix<double, 2, 2> Q;
    Eigen::Matrix<double, 2, 1> c;
    Eigen::Matrix<double, 2, 1> x, z;        // decision variables
    Eigen::Matrix<double, -1, 2> A(m, 2); // constraint matrix
    Eigen::VectorXd b = -Eigen::VectorXd::Ones(m);                 // constraint bound

    Q << 1.0, 0.0, 0.0, 1.0;
    c << 0.0, 0.0;

    A << 2.83404400940515, 2.83838902880659,
        3.44064898688432, 3.37043900079352,
        2.00022874963469, 2.40890449946304,
        2.60466514526368, 3.75623487278189,
        2.29351178163423, 2.05477518639585,
        2.18467718953760, 3.34093502035680,
        2.37252042275534, 2.83460960473425,
        2.69112145408610, 3.11737965689150,
        2.79353494846134, 2.28077387719047,
        3.07763346800671, 2.39620297816976;
    A = -A;
    // std::cout << A << std::endl;
    // std::cout << b << std::endl;

    double minobj = sdqp::sdqp<2>(Q, c, A, b, z);

    x = z / z.dot(z);

    std::cout << "optimal sol: " << x.transpose() << std::endl;
    std::cout << "min collision distance: " << x.norm() << std::endl;
    std::cout << "cons precision: " << (A * z - b).maxCoeff() << std::endl;

    

    return 0;
}
