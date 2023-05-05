#ifndef SOLUTIONCALCULATOR_H
#define SOLUTIONCALCULATOR_H

#include "LUDecomposition.h"
#include <fstream>

class SolutionCalculator{
public:
    SolutionCalculator(const LUDecomposition& _lu, std::ifstream* _bVecFile, const int _myid, const int _numProcs, const unsigned _n);
    ~SolutionCalculator();
    void run();
    void printSolutionVector() const;

    double mse(std::ifstream* targetVec);

private:
    void calculateYVector();
    void calculateXVector();
    void assignColumnsToProcess();

    Matrix L;
    Matrix U;
    double* y;
    double* yAll;
    double* x;
    double* b;

    const int myid;
    const int numProcs;
    const unsigned n;
    unsigned maxNumberOfCollumnsOwned;
    int* mycols;

};

#endif // SOLUTIONCALCULATOR_H