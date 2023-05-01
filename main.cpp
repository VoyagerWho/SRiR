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
	unsigned n = 1600;
	if (myid == 0)
	{
		std::ifstream file("Matrix.txt");
		Matrix m = Matrix(file, n);
		// cout<<"A:\n"<<m<<endl;
		double startTime = MPI_Wtime();
		LUDecomposition lu(m[0], n, myid, numProcs);
		// LUDecomposition lu(m);
		double endTime = MPI_Wtime();
		// cout<<"L:\n"<<lu.getL()<<endl;
		// cout<<"U:\n"<<lu.getU()<<endl;
		// cout<<"[LU]:\n"<<lu.getLU()<<endl;
		Matrix mlu = lu.getL() * lu.getU();
		// cout<<"LU:\n"<<mlu<<endl;
		cout << "MSE: " << mse(m, mlu) << endl;
		cout << "deltaTime: " << endTime - startTime << endl;
	}
	else
	{
		LUDecomposition lu(nullptr, n, myid, numProcs);
	}
	MPI_Finalize();
	return 0;
}