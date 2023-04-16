#include <iostream>
#include <fstream>
#include "Matrix.h"
#include "LUDecomposition.h"

using namespace std;

int main()
{
    cout << "Hello world!" << endl;
    unsigned n = 8;
    std::ifstream file("Matrix3D.txt");
    Matrix m = Matrix(file, n);
    cout<<m<<endl;
//    LUDecomposition lu(m);
//    cout<<lu.getL()<<endl;
//    cout<<lu.getU()<<endl;
//    cout<<lu.getLU()<<endl;
//    cout<<lu.getL()*lu.getU()<<endl;
    cout<<Matrix::eye(n)<<endl;
    cout<<Matrix::ones(n)<<endl;
    cout<<Matrix::eye(n)+Matrix::ones(n)<<endl;
    cout<<Matrix::eye(n)-Matrix::ones(n)<<endl;
    cout<<Matrix::eye(n)+3.14<<endl;
    cout<<Matrix::eye(n)-3.14<<endl;
    return 0;
}
