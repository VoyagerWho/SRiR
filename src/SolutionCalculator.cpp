#include <math.h>
#include "SolutionCalculator.h"

SolutionCalculator::SolutionCalculator(const LUDecomposition& _lu, std::ifstream* _bVecFile, 
                                        const int _myid, const int _numProcs, const unsigned _n):
                                        L(_lu.getL()), U(_lu.getU()), myid(_myid), 
                                        numProcs(_numProcs), n(_n)
{
    assignColumnsToProcess();

    x = new double[n]();
    y = new double[n]();
    b = new double[n]();
    if(myid == 0)
    {
        double *ptr = b;
        for (unsigned i = 0; i < n; ++i)
        {
            *_bVecFile >> *ptr++;
        }
    }

    MPI_Bcast(b, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
}

SolutionCalculator::~SolutionCalculator()
{
    delete[] y;
    delete[] x;
    delete[] b;
}

void SolutionCalculator::printSolutionVector() const
{
    for(int i = 0; i < (int)n; i++)
    {
        printf("% 8.3f\n ", x[i]);
    }
}

void SolutionCalculator::run()
{
    calculateYVector();
    calculateXVector();
    MPI_Barrier(MPI_COMM_WORLD);
}

double SolutionCalculator::mse(std::ifstream* targetVec)
{
    double* target = new double[n]();
    double *ptr = target;
    for (unsigned i = 0; i < n; ++i)
    {
        *targetVec >> *ptr++;
    }

	double error = 0;
	for (int i = 0; i < (int)n; i++) {
		error += pow((target[i] - x[i]), 2);
	}
    delete[] target;
	return error / n;
}

//-----------------------------------
// Private methods definitions
//-----------------------------------

void SolutionCalculator::calculateYVector()
{
    double* yPartial = new double[n]();
    bool ownCurrentRow;
    double partialRowSum;
    double globalSum;
    for(int row = 0; row < (int)n; row++)
    {
        if (myid == 0 && !(row % 100))
            std::cout << "y: " << row << "\n";
        partialRowSum = 0;
        globalSum = 0;
        ownCurrentRow = false;

        for(unsigned collumnOwned = 0; collumnOwned < maxNumberOfCollumnsOwned; collumnOwned++)
        {
            if((mycols[collumnOwned] < row) && mycols[collumnOwned] != -1)
            {
                partialRowSum += L[row][mycols[collumnOwned]] * yPartial[mycols[collumnOwned]];
            }
            if(mycols[collumnOwned] == row)
                ownCurrentRow = true;
        }
        
        MPI_Allreduce( &partialRowSum, &globalSum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
        if(ownCurrentRow && abs(L[row][row]) > 0.000001)
        {
            yPartial[row] = (b[row] - globalSum) / L[row][row];
        }
    }

    MPI_Allreduce(yPartial, y, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    delete[] yPartial;
}

void SolutionCalculator::calculateXVector()
{
    double* xPartial = new double[n]();
    bool ownCurrentRow;
    double partialRowSum;
    double globalSum;

    for(int row = n-1; row >= 0; row--)
    {
        if (myid == 0 && !(row % 100))
            std::cout << "x: " << row << "\n";
        partialRowSum = 0;
        globalSum = 0;
        ownCurrentRow = false;

        for(unsigned collumnOwned = 0; collumnOwned < maxNumberOfCollumnsOwned; collumnOwned++)
        {
            if((mycols[collumnOwned] > row) && mycols[collumnOwned] != -1)
            {
                partialRowSum += U[row][mycols[collumnOwned]] * xPartial[mycols[collumnOwned]];
            }
            if(mycols[collumnOwned] == row)
                ownCurrentRow = true;
        }
        
        MPI_Allreduce( &partialRowSum, &globalSum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
        if(ownCurrentRow && abs(L[row][row]) > 0.000001)
        {
            xPartial[row] = (y[row] - globalSum) / U[row][row];
        }
    }

    MPI_Allreduce(xPartial, x, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    delete[] xPartial;
}

//Cyclic assignement of columns
void SolutionCalculator::assignColumnsToProcess()
{
    maxNumberOfCollumnsOwned = ceil((double)n/numProcs);
    mycols = new int[maxNumberOfCollumnsOwned]();
    for(unsigned i = 0; i < maxNumberOfCollumnsOwned; i++)
    {
        mycols[i] = -1;
    }

    int numberOfCollumnsAssigned = 0;
    for(unsigned col = 0; col < n; col++)
    {
        if(col == (unsigned)(numberOfCollumnsAssigned*numProcs + myid))
        {
            mycols[numberOfCollumnsAssigned] = col;
            numberOfCollumnsAssigned++;
        }
    }
}