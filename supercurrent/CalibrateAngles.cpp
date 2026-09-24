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

///// Change here /////

string PairName = "co0";
string InputFileName = PairName + "_cpp.txt";
string OutputFileName = PairName + "_dataIsCalibrationNeeded.m";
string OutputListName = PairName + "IsCalibrationNeeded"; // In Mathematica, the underscore '_' cannot be used in variable names.
string LogFileName = PairName + "_logActionOnSupercurrent.txt";



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


    
    ////////////////////////
    ///// Calibrate the angles for each representative of SL(2,Z)-orbit
    ////////////////////////

    ///// Initialize variables

    int n = 2048;
    std::ofstream fout, flog;


    ///// Load the supercurrent

    Eigen::VectorXcd sc(n);
    fin.open("datacpp_supercurrent_ChiSpin.txt");
    for (int i = 0; i < n; i++) {
        double re, im;
        fin >> re >> im;
        sc(i) = complex<double>(re, im);
    }
    fin.close();

    // Make positive-chiral spinor into full spinor
    sc = ChiSpinToFullSpin(sc);

    cout << "The supercurrent is prepared." << endl;



    ///// Calibrate the angles of commuting pairs

    double eps = 1e-11;

    int num_pairs;
    vector<double> angles_first(12), angles_second(12);
    Eigen::MatrixXd eigenmat(24, 24);

    fin.open(InputFileName);
    fout.open(OutputFileName);
    fout << OutputListName << " = {";
    flog.open(LogFileName);
    
    fin >> num_pairs;

    for (int r = 0; r < num_pairs; r++) {
        cout << r << " "; // log
        if (r > 0) {
            fout << "," << "\n";
        }

        // Load data
        for (int i = 0; i < 12; i++)
            fin >> angles_first.at(i);
        for (int i = 0; i < 12; i++)
            fin >> angles_second.at(i);
        for (int i = 0; i < 24; i++)
            for (int j = 0; j < 24; j++)
                fin >> eigenmat(i, j);

        // Calculate whether the supercurrent is preserved or reversed
        Eigen::VectorXcd first_sc = SO24LiftAction(angles_first, eigenmat, sc);
        // One of the followings must be 0, and the other 2.
        double first_preserved = (sc - first_sc).norm();
        double first_reversed = (sc + first_sc).norm();
        Eigen::VectorXcd second_sc = SO24LiftAction(angles_second, eigenmat, sc);
        // One of the followings must be 0, and the other 2.
        double second_preserved = (sc - second_sc).norm();
        double second_reversed = (sc + second_sc).norm();

        flog << r << ": " << "\n";
        flog << first_preserved << " " << first_reversed << " " << second_preserved << " " << second_reversed << "\n";

        // Judge whether the angles need to be calibrated or not
        fout << "{";
        // Deal with the first set of angles
        if (first_preserved < eps && abs(2 - first_reversed) < eps) {
            fout << 0;
        }
        else if (abs(2 - first_preserved) < eps && first_reversed < eps)
        {
            fout << 1;
        }
        else {
            cout << "Something is wrong at the first set of angles of " << r << "-th pair" << "\n";
            cout << first_preserved << " " << abs(2 - first_reversed) << " " << abs(2 - first_preserved) << " " << first_reversed << "\n";
            break;
        }
        fout << ", ";
        // Deal with the second set of angles
        if (second_preserved < eps && abs(2 - second_reversed) < eps) {
            fout << 0;
        }
        else if (abs(2 - second_preserved) < eps && second_reversed < eps)
        {
            fout << 1;
        }
        else {
            cout << "Something is wrong at the second set of angles of " << r << "-th pair" << "\n";
            cout << second_preserved << " " << abs(2 - second_reversed) << " " << abs(2 - second_preserved) << " " << second_reversed << "\n";
            break;
        }
        fout << "}";
    }

    flog.close();
    fout << "}";
    fout.close();
    fin.close();
}