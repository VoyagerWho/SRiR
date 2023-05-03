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

double* SolutionCalculator::getSolutionVector() const
{
    return x;
}

void SolutionCalculator::run()
{
    calculateYVector();
    // calculateXVector();
}

//-----------------------------------
// Private methods definitions
//-----------------------------------

void SolutionCalculator::calculateYVector()
{
    bool ownCurrentRow;
    double partialRowSum;
    double globalSum;
    for(int row = 0; row < (int)n; row++)
    {
        std::cout << "[" << myid << "] begin row " << row << std::endl;
        partialRowSum = 0;
        globalSum = 0;
        ownCurrentRow = false;

        for(unsigned collumnOwned = 0; collumnOwned < maxNumberOfCollumnsOwned; collumnOwned++)
        {
            // std::cout << "[" << myid << "] searching for collumn " << collumnOwned << std::endl;
            if((mycols[collumnOwned] < row) && mycols[collumnOwned] != -1)
            {
                // std::cout << "[" << myid << "] own collumn " << mycols[collumnOwned] << std::endl;
                partialRowSum += L[row][mycols[collumnOwned]] * y[mycols[collumnOwned]];
            }
            if(mycols[collumnOwned] == row)
                ownCurrentRow = true;
        }
        
        std::cout << "[" << myid << "] executing reduce now " << std::endl;
        MPI_Allreduce( &partialRowSum, &globalSum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
        std::cout << "[" << myid << "] reduce complete, globalSum: " << globalSum << std::endl;
        if(ownCurrentRow && L[row][row] != 0.0)
        {
            std::cout << "[" << myid << "] calculating y of row " << row << std::endl;
            y[row] = (b[row] - globalSum) / L[row][row];
        }
    }

}

void SolutionCalculator::calculateXVector()
{

    std::cout << "[" << myid << "] beggin x calc" << std::endl;
    for(int col = n-1; col >= 0; col--)
    {   
        if(U[col][col] != 0.0)
            x[col] = y[col] / U[col][col];
        for(int row = col-1; row >= 0; row--)
        {
            y[row] -= U[row][col]*x[col];
        }
    }
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
            std::cout << "[" << myid << "] I own col " << col << std::endl;
            mycols[numberOfCollumnsAssigned] = col;
            numberOfCollumnsAssigned++;
        }
    }
}