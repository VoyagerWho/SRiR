#include <iostream>
#include <fstream>
#include <mpi.h>
#include "src/Matrix.h"
#include "src/LUDecomposition.h"

using namespace std;

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    int myid, numProcs;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
		unsigned n = 8;
		if(myid==0)
		{
			std::ifstream file("Matrix3D.txt");
			Matrix m = Matrix(file, n);
			cout<<"A:\n"<<m<<endl;
			LUDecomposition lu(m[0], n, myid, numProcs);
			cout<<"L:\n"<<lu.getL()<<endl;
			cout<<"U:\n"<<lu.getU()<<endl;
			cout<<"[LU]:\n"<<lu.getLU()<<endl;
			Matrix mlu = lu.getL()*lu.getU();
			cout<<"LU:\n"<<mlu<<endl;
			cout<<"MSE: "<<mse(m, mlu)<<endl;
		}
		else
		{
			LUDecomposition lu(nullptr, n, myid, numProcs);
		}
		MPI_Finalize();
    return 0;
}
