#define _USE_MATH_DEFINES
#include <iostream>
#include <iomanip>
#include <fstream>
using namespace std;

#include <vector>
#include <array>
#include <algorithm>

#include <Eigen/Dense>
#include <Eigen/SVD>
#include <Eigen/QR>
#include <complex>
#include <cmath>
#include <bit>
#include <chrono>
#include <limits>



///// Transformation matrix from complex spinor to real spinor

Eigen::MatrixXcd makeCpxSpinorToRealSpinor() {
    using namespace std::complex_literals;

    // the 2x2 block
    Eigen::Matrix2cd block;
    block << 1.0, -1.0i,
            1.0, 1.0i;

    // the 24x24 zero matrix
    Eigen::MatrixXcd mat = Eigen::MatrixXcd::Zero(24, 24);

    // Place the 2x2 block diagonally
    for (int k = 0; k < 12; k++) {
        mat.block<2,2>(2*k, 2*k) = block;
    }

    // Normalize
    mat *= 1.0 / std::sqrt(2.0);

    return mat;
}

Eigen::MatrixXcd CpxSpinorToRealSpinor = makeCpxSpinorToRealSpinor();



///// the Clifford algebra action on spinors

std::vector<std::vector<std::array<double,2>>> PsiActionList(24, std::vector<std::array<double,2>>(4096)); // loaded later in main()

Eigen::VectorXcd PsiAction(int i, const Eigen::VectorXcd& spinor) {
    int dim = 1 << 12; // 2^12
    Eigen::VectorXcd result = Eigen::VectorXcd::Zero(dim);
    for (int d = 0; d < dim; d++) {
        int index = static_cast<int>(PsiActionList[i][d][0]);
        double coeff = PsiActionList[i][d][1];
        result[d] = coeff * spinor[index];
    }
    return result;
}

Eigen::VectorXcd PsiLinCombiAction(const Eigen::VectorXcd& coeffList,
                                   const Eigen::VectorXcd& spinor) {
    int dim = spinor.size();
    Eigen::VectorXcd result = Eigen::VectorXcd::Zero(dim);
    int terms = coeffList.size();
    for (int i = 0; i < terms; i++) {
        result += coeffList[i] * PsiAction(i, spinor);
    }
    return result;
}

Eigen::VectorXcd SO24LiftAction(const std::vector<double>& angles,
                                const Eigen::MatrixXd& eigenmat,
                                const Eigen::VectorXcd& spinor) {
    Eigen::MatrixXcd eigenmatCpx = CpxSpinorToRealSpinor * eigenmat;
    Eigen::VectorXcd temp = spinor;

    for (int i = 12; i >= 1; i--) {
        Eigen::VectorXcd cosTemp = std::cos(angles[i-1] * M_PI) * temp;

        Eigen::VectorXcd coeffList1 = eigenmatCpx.col(2*i-1);
        Eigen::VectorXcd coeffList2 = eigenmatCpx.col(2*i-2);

        temp = PsiLinCombiAction(coeffList1, temp);
        temp = PsiLinCombiAction(coeffList2, temp);

        temp = std::sin(angles[i-1] * M_PI) * temp;
        temp = cosTemp + temp;
    }

    return temp;
}



///// Retain only spinors with positive chirality (against expectations, a bit slower than above)

std::vector<std::vector<std::array<double,2>>> PsiPsiActionList(276, std::vector<std::array<double,2>>(2048));

Eigen::VectorXcd PsiPsiAction(int x, const Eigen::VectorXcd& spinorChi) {
    int dim = 1 << 11; // 2^11
    Eigen::VectorXcd result = Eigen::VectorXcd::Zero(dim);
    for (int d = 0; d < dim; d++) {
        int index = static_cast<int>(PsiPsiActionList[x][d][0]);
        double coeff = PsiPsiActionList[x][d][1];
        result[d] = coeff * spinorChi[index];
    }
    return result;
}

Eigen::VectorXcd PsiPsiLinCombiAction(const Eigen::VectorXcd& PsiPsiCoeffList,
                                   const Eigen::VectorXcd& spinorChi) {
    int dim = spinorChi.size();
    Eigen::VectorXcd result = Eigen::VectorXcd::Zero(dim);
    int terms = PsiPsiCoeffList.size();
    for (int i = 0; i < terms-1; i++) {
        result += PsiPsiCoeffList[i] * PsiPsiAction(i, spinorChi);
    }
    result += PsiPsiCoeffList[terms-1] * spinorChi; // constant term
    return result;
}

// Multiplication of two linear combinations of Psi's
Eigen::VectorXcd MultiOfPsiLinCombi(const Eigen::VectorXcd& coeffList1,
                                    const Eigen::VectorXcd& coeffList2) {
    int n = 24;
    Eigen::VectorXcd result((n*(n-1))/2 + 1); // number of PsiPsi's + constant term
    int idx = 0;

    // coefficient of Psi_i Psi_j
    for (int i = 0; i < n-1; ++i) { // i = 0..22
        for (int j = i+1; j < n; ++j) { // j = i+1..23
            result(idx) = coeffList1(i) * coeffList2(j) - coeffList1(j) * coeffList2(i);
            idx++;
        }
    }

    // constant term: -2 * Sum[coeffList1[[2 i]] coeffList2[[2 i -1]], {i, 1, 12}] in Mathematica (1-indexed)
    std::complex<double> constant = 0.0;
    for (int i = 0; i < 12; ++i) {
        constant += coeffList1(2*i + 1) * coeffList2(2*i);
    }
    constant *= -2.0;
    result(idx) = constant;

    return result;
}

Eigen::VectorXcd SO24LiftActionChi(const std::vector<double>& angles,
                                   const Eigen::MatrixXd& eigenmat,
                                   const Eigen::VectorXcd& spinorChi) {
    Eigen::MatrixXcd eigenmatCpx = CpxSpinorToRealSpinor * eigenmat;
    Eigen::VectorXcd temp = spinorChi;

    for (int i = 12; i >= 1; --i) {
        Eigen::VectorXcd cosTemp = std::cos(angles[i-1] * M_PI) * temp;

        Eigen::VectorXcd coeffList1 = eigenmatCpx.col(2*i - 2);
        Eigen::VectorXcd coeffList2 = eigenmatCpx.col(2*i - 1);

        Eigen::VectorXcd multiCoeffs = MultiOfPsiLinCombi(coeffList1, coeffList2);

        temp = PsiPsiLinCombiAction(multiCoeffs, temp);

        temp = std::sin(angles[i-1] * M_PI) * temp;
        temp = cosTemp + temp;
    }

    return temp;
}

// Make positive-chiral spinor into full spinor
Eigen::VectorXcd ChiSpinToFullSpin(const Eigen::VectorXcd& ChiSpin) {
    int dim = 1 << 12; // 2^12
    Eigen::VectorXcd result = Eigen::VectorXcd::Zero(dim);
    for (int d = 0; d < dim; d++) {
        if (__builtin_popcount(d) % 2 == 0) {
            result(d) = ChiSpin(d >> 1);
        }
    }
    return result;
}



///// Utilities /////

// matrix power
template <typename MatrixType>
MatrixType matrixPower(const MatrixType& A, unsigned int n) {
    if (n == 0) {
        return MatrixType::Identity(A.rows(), A.cols());
    }
    if (n == 1) {
        return A;
    }
    MatrixType half = matrixPower(A, n / 2);
    if (n % 2 == 0) {
        return half * half;
    } else {
        return half * half * A;
    }
}





///// main /////

int main() {

    ////////////////////////
    ///// Initial settings
    ////////////////////////

    std::ifstream fin("cppPsiActionList.txt");
    for (int i = 0; i < 24; i++) {
        for (int d = 0; d < 4096; d++) {
            fin >> PsiActionList[i][d][0] >> PsiActionList[i][d][1];
        }
    }
    fin.close();

    fin.open("cppPsiPsiActionList.txt");
    for (int i = 0; i < 276; i++) {
        for (int d = 0; d < 2048; d++) {
            fin >> PsiPsiActionList[i][d][0] >> PsiPsiActionList[i][d][1];
        }
    }
    fin.close();


    
    ////////////////////////
    ///// Find the supercurrent as the common kernel of A-I and B-I
    ////////////////////////

    ///// Load files of the generators A and B

    std::vector<double> anglesA(12);
    fin.open("cppanglesA.txt");
    for (int i=0; i<12; i++) {
        fin >> anglesA.at(i);
    }
    fin.close();

    std::vector<double> anglesB(12);
    fin.open("cppanglesB.txt");
    for (int i=0; i<12; i++) {
        fin >> anglesB.at(i);
    }
    fin.close();

    Eigen::MatrixXd eigenmatA(24, 24);
    fin.open("cppeigenmatA.txt");
    for (int i = 0; i < 24; i++)
        for (int j = 0; j < 24; j++)
            fin >> eigenmatA(i, j);
    fin.close();

    Eigen::MatrixXd eigenmatB(24, 24);
    fin.open("cppeigenmatB.txt");
    for (int i = 0; i < 24; i++)
        for (int j = 0; j < 24; j++)
            fin >> eigenmatB(i, j);
    fin.close();

    ///// Initialize variables

    int n = 2048;
    std::ofstream fout;
    
    ////////////////////////
    // Step 1: Generate representation matrices of A and B on positive-chiral spinors
    ////////////////////////

    ///// Generate representation matrices of A and B on positive-chiral spinors

    // Generate representation matrix of A on positive-chiral spinors
    Eigen::MatrixXcd A(n, n);
    for (int i = 0; i < n; i++) {
        Eigen::VectorXcd basisChi = Eigen::VectorXcd::Zero(n);
        basisChi[i] = std::complex<double>(1.0, 0.0);
        Eigen::VectorXcd v = SO24LiftActionChi(anglesA, eigenmatA, basisChi);
        A.col(i) = v; // record as a column vector
        cout << i << " "; // log
    }


    // Generate representation matrix of B on positive-chiral spinors
    Eigen::MatrixXcd B(n, n);
    for (int i = 0; i < n; i++) {
        Eigen::VectorXcd basisChi = Eigen::VectorXcd::Zero(n);
        basisChi[i] = std::complex<double>(1.0, 0.0);
        Eigen::VectorXcd v = SO24LiftActionChi(anglesB, eigenmatB, basisChi);
        B.col(i) = v; // record as a column vector
        cout << i << " "; // log
    }

    

    ////////////////////////
    // Step 2: Compute the common kernel of A-I and B-I
    ////////////////////////

    

    ///// Find a basis QA of Ker(A-I) and a basis QB of Ker(B-I)

    // A-I and B-I
    Eigen::MatrixXcd AI = A - Eigen::MatrixXcd::Identity(n,n);
    Eigen::MatrixXcd BI = B - Eigen::MatrixXcd::Identity(n,n);

    // We use LU decomposition to construct bases of kernels as follows.
    Eigen::FullPivLU<Eigen::MatrixXcd> luAI(AI);
    luAI.setThreshold(1e-12);
    Eigen::MatrixXcd QA = luAI.kernel(); // basis of kernel of A-I
    
    Eigen::FullPivLU<Eigen::MatrixXcd> luBI(BI);
    luBI.setThreshold(1e-12);
    Eigen::MatrixXcd QB = luBI.kernel(); // basis of kernel of B-I

    cout << "rank Ker(A-I) = " << QA.cols() << ", rank Ker(B-I) = " << QB.cols() << endl;

    // orthonormalize the bases of kernels
    auto orthonormalize = [&](const Eigen::MatrixXcd& X){
        Eigen::HouseholderQR<Eigen::MatrixXcd> qr(X);
        Eigen::MatrixXcd Q = Eigen::MatrixXcd::Identity(X.rows(), X.cols());
        Q = qr.householderQ() * Q;
        return Q;
    };
    QA = orthonormalize(QA);
    QB = orthonormalize(QB);

    // Sanity checks
    std::cout << "||QA^*QA - I|| = " << (QA.adjoint()*QA - Eigen::MatrixXcd::Identity(QA.cols(), QA.cols())).norm() << "\n";
    std::cout << "||QB^*QB - I|| = " << (QB.adjoint()*QB - Eigen::MatrixXcd::Identity(QB.cols(), QB.cols())).norm() << "\n";
    fout.open("test_precisionQAQB.txt");
    fout.setf(std::ios::scientific);
    fout << std::setprecision(std::numeric_limits<double>::max_digits10);
    for (int j=0; j<QA.cols(); ++j) {
        fout << "A-res col " << j << " = " << (A*QA.col(j) - QA.col(j)).norm() << "\n";
    }
    for (int j=0; j<QB.cols(); ++j) {
        fout << "B-res col " << j << " = " << (B*QB.col(j) - QB.col(j)).norm() << "\n";
    }
    fout.close();

    ///// Solve QA * x = QB * y, by finding kernel of C = [QA, -QB]

    const int p = QA.cols();
    const int q = QB.cols();
    Eigen::MatrixXcd C(n, p + q); // C = [QA, -QB] (size: n x (p+q))
    C.leftCols(p) = QA;
    C.rightCols(q) = -QB;

    // SVD of C = [QA, -QB] = U S V*
    // Just in the case of p+q > n, it is safe to require FullV.
    Eigen::JacobiSVD<Eigen::MatrixXcd> svdC(C, Eigen::ComputeThinU | Eigen::ComputeFullV);

    // Sanity check
    const Eigen::VectorXd sC = svdC.singularValues();
    fout.open("test_SingularValuesQA-QB.txt");
    fout.setf(std::ios::scientific);
    fout << std::setprecision(std::numeric_limits<double>::max_digits10);
    fout << "size: " << sC.size() << endl;
    for (int i=0; i<sC.size(); i++) fout << sC[i] << "\n"; // The last singular value must be (numerically) zero.
    fout.close();

    // z = [x^T, y^T]^T
    const Eigen::VectorXcd z = svdC.matrixV().col(svdC.matrixV().cols() - 1);
    Eigen::VectorXcd x = z.topRows(p);
    Eigen::VectorXcd y = z.bottomRows(q);

    ///// The supercurrent is u = QA * x = QB * y

    Eigen::VectorXcd uA = QA * x;
    Eigen::VectorXcd uB = QB * y;
    Eigen::VectorXcd u = uA +uB; // take average to stabilize u numerically

    // Normalize u
    if (uA.norm() > 0) uA /= uA.norm();
    if (uB.norm() > 0) uB /= uB.norm();
    if (u.norm() > 0) u /= u.norm();

    // Sanity checks; everything must be numerically zero
    cout << "||C z|| = " << (C * z).norm() << "\n";
    cout << "||QA x - QB y|| = " << (uA - uB).norm() << "\n";
    cout << "||A uA - uA|| = " << (A*uA - uA).norm() << "\n";
    cout << "||A u - u|| = " << (A*u - uA).norm() << "\n";
    cout << "||B uB - uB|| = " << (B*uB - uB).norm() << "\n";
    cout << "||B u - u|| = " << (B*u - u).norm() << "\n";

    // Save the supercurrent 

    fout.open("datacpp_supercurrent_ChiSpin.txt");
    fout.setf(std::ios::scientific);
    fout << std::setprecision(std::numeric_limits<double>::max_digits10);
    for (int i = 0; i < u.size(); i++) { fout << u(i).real() << " " << u(i).imag() << endl; }
    fout.close();

    

    ////////////////////////
    // Step 3: Sanity checks
    ////////////////////////
	
    // Make positive-chiral spinor into full spinor
    sc = ChiSpinToFullSpin(sc);

    ///// Sanity checks of the supercurrent
    
    cout << "Sanity Checks of the supercurrent" << endl;

    // invariant under the generators A and B
    Eigen::VectorXcd A_sc = SO24LiftAction(anglesA, eigenmatA, sc);
    cout << (sc - A_sc).norm() << endl; // must be 0
    cout << (sc + A_sc).norm() << endl; // must be 2
    Eigen::VectorXcd B_sc = SO24LiftAction(anglesB, eigenmatB, sc);
    cout << (sc - B_sc).norm() << endl; // must be 0
    cout << (sc + B_sc).norm() << endl; // must be 2


    

}