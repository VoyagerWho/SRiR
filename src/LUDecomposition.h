#ifndef LUDECOMPOSITION_H
#define LUDECOMPOSITION_H
#include "Matrix.h"
//#include <mpi.h>

class LUDecomposition
{
  public:
    LUDecomposition(const Matrix& A);
    LUDecomposition(const double* orgA, const unsigned matSize, const int myid, const int numOfProcs);
    ~LUDecomposition();
    Matrix getL() const;
    Matrix getU() const;
    Matrix getLU() const;

  protected:
  private:
    void transpose(const double* matrix, double* transposed, unsigned matSize) const;
    void getColumns(const double *matrix, double *col, unsigned matSize, int processNum, unsigned columns) const;
    void getRows(const double *matrix, double *row, unsigned matSize, int processNum, unsigned rows) const;
    Matrix u;
    Matrix l;
};

#endif // LUDECOMPOSITION_H
