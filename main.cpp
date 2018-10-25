#include <iostream>

#include <eigen3/Eigen/Dense>
//#include <Eigen/QR>


using namespace Eigen;
using namespace std;

/* @brief QR decomposition from Cholesky decomposition
 * @param[in] A An $m \times n$ matrix
 * @param[out] R The upper triangular matrix from the QR decomposition of $A$
 * @param[out] Q The orthogonal matrix from the QR decomposition of $A$
 */
void CholeskyQR(const Eigen::MatrixXd &A, Eigen::MatrixXd &R, Eigen::MatrixXd &Q) {

    MatrixXd A2 = A.transpose() * A;

    MatrixXd L = A2.llt().matrixL();
    R = L.transpose();
    Q = A * (L.inverse()).transpose();

}

/* @brief Direct QR decomposition
 * @param[in] A An $m \times n$ matrix
 * @param[out] R The upper triangular matrix from the QR decomposition of $A$
 * @param[out] Q The orthogonal matrix from the QR decomposition of $A$
 */
void DirectQR(const Eigen::MatrixXd &A, Eigen::MatrixXd &R, Eigen::MatrixXd &Q) {

//    ColPivHouseholderQR<MatrixXd> qr(A);
//    MatrixXd thinQ(MatrixXd::Identity(A.rows(), A.cols()));
//    Q = qr.householderQ().setLength(qr.nonzeroPivots()) * thinQ;




    HouseholderQR<MatrixXd> qr(A);
    MatrixXd thinQ(MatrixXd::Identity(A.rows(),A.cols()));
    Q = qr.householderQ() * thinQ;

   MatrixXd R2 = qr.matrixQR().triangularView<Upper>();
   long size = Q.cols();
    R = R2.block(0,0,size,size);

}

int main() {
    int m = 3;
    int n = 2;

    Eigen::MatrixXd A(m, n);
    double epsilon = 1.e-8;
//    A << 3, 5, 1, 9, 7, 1;
    A << 1, 1, epsilon, 0, 0, epsilon;
    std::cout << "A =" << std::endl << A << std::endl;

    Eigen::MatrixXd R, Q;

    CholeskyQR(A, R, Q);
    std::cout << "CholeskyQR: ===========" << std::endl;
    std::cout << "R =" << std::endl << R << std::endl;
    std::cout << "Q =" << std::endl << Q << std::endl;

    DirectQR(A, R, Q);
    std::cout << "DirectQR: =============" << std::endl;
    std::cout << "R =" << std::endl << R << std::endl;
    std::cout << "Q =" << std::endl << Q << std::endl;

    return 0;
}
