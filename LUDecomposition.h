#ifndef LUDECOMPOSITION_H
#define LUDECOMPOSITION_H
#include "Matrix.h"

class LUDecomposition
{
  public:
    LUDecomposition(const Matrix& A);
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
