#ifndef LUDECOMPOSITION_H
#define LUDECOMPOSITION_H
#include "Matrix.h"
#include <mpi.h>

class LUDecomposition
{
public:
  LUDecomposition(const Matrix &A);
  LUDecomposition(const double *orgA, const unsigned matSize, const int myid, const int numOfProcs);
  ~LUDecomposition();
  Matrix getL() const;
  Matrix getU() const;
  Matrix getLU() const;

protected:
private:
  Matrix u;
  Matrix l;
};

#endif // LUDECOMPOSITION_H
