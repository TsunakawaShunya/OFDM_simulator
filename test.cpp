#include <Eigen/Core>
#include <Eigen/Eigen>
#include <iostream>

int main() {
    Eigen::MatrixXd mat(3, 3);
    mat << 1, 1, 1,
           1, 1, 1,
           1, 1, 1;
    
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(mat);
    Eigen::MatrixXd Lambda = es.eigenvalues().asDiagonal();

    // 1e-10以下の成分を0にする
    for (int i = 0; i < Lambda.rows(); ++i) {
        if (std::abs(Lambda(i, i)) < 1e-10) {
            Lambda(i, i) = 0;
        }
    }

    std::cout << "Lambda = " << Lambda << std::endl;
    Eigen::MatrixXd U = es.eigenvectors();
    std::cout << "U = " << U << std::endl;

    std::cout << "mat =  " << U * Lambda * U.inverse() << std::endl;
    
    std::cout << "A_q = " << U * Lambda.cwiseSqrt() << std::endl;
    return 0;
}