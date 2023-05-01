#include "LUDecomposition.h"

LUDecomposition::LUDecomposition(const Matrix &A)
    : u(A), l(Matrix(A.size()))
{
  unsigned n = u.size();
  for (unsigned i = 0; i < n; ++i)
  {
    l[i][i] = 1.0;
  }
  for (unsigned k = 0; k < n - 1; ++k)
  {
    for (unsigned i = k + 1; i < n; ++i)
      l[i][k] = u[i][k] / u[k][k];
    for (unsigned j = k + 1; j < n; ++j)
    {
      for (unsigned i = k + 1; i < n; ++i)
        u[i][j] = u[i][j] - l[i][k] * u[k][j];
    }
  }
  for (unsigned i = 0; i < n - 1; ++i)
  {
    for (unsigned j = i + 1; j < n; ++j)
    {
      u[j][i] = 0;
    }
  }
}

LUDecomposition::LUDecomposition(const double* orgA, const unsigned matSize, const int myid, const int numOfProcs)
    : u(Matrix(matSize)), l(Matrix(matSize))
{
  const unsigned rows = matSize / numOfProcs;

  double *A = new double[matSize*matSize]();
  double *Ut = nullptr;
  double *locL = new double[rows * matSize]();
  double *locU = new double[rows * matSize]();
  if (myid == 0)
  {
    Ut = new double[matSize * matSize]();
    double* ptr = A;
    for(const double* i=orgA;unsigned(i-orgA)<matSize*matSize;++i)
    {
      *ptr++ = *i;
    }
  }
  MPI_Bcast(A, matSize*matSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  u = Matrix(A, matSize);

  for (unsigned k = 0; k < matSize; ++k)
  {
    for (unsigned i = myid*rows; i < (myid+1)*rows; ++i)
    {
      if (i < k)
      {
        l[i][k] = 0.0;
      }
      else if (i == k)
      {
        l[i][k] = 1.0;
      }
      else
      {
        l[i][k] = u[i][k] / u[k][k] ;
      }

    }

    getRows(l[0], locL, matSize, myid, rows);
    MPI_Gather(locL, matSize*rows, MPI_DOUBLE, l[0], matSize*rows, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(l[0], matSize*matSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);


    for (unsigned i = myid*rows; i < (myid+1)*rows; ++i)
    {
      if (i < k)
      {
        u[k][i] = 0.0;
      }
      else
      {
        for (unsigned j = k+1; j < matSize; ++j)
          u[j][i] = u[j][i] - l[j][k] * u[k][i];
      }

    }
    getColumns(u[0], locU, matSize, myid, rows);
    MPI_Gather(locU, matSize*rows, MPI_DOUBLE, Ut, matSize*rows, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if (myid == 0)
      transpose(Ut, u[0], matSize);
    MPI_Bcast(u[0], matSize*matSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  }
  delete[] A;
  if (myid == 0)
      delete[] Ut;
  delete[] locL;
  delete[] locU;
}

LUDecomposition::~LUDecomposition()
{
  // dtor
}

void LUDecomposition::transpose(const double *matrix, double *transposed, unsigned matSize) const
{
  for (unsigned i = 0; i < matSize; ++i)
  {
    for (unsigned j = 0; j < matSize; ++j)
    {
      transposed[i * matSize + j] = matrix[i * matSize + j];
    }
  }
}

void LUDecomposition::getColumns(const double *matrix, double *col, unsigned matSize, int processNum, unsigned columns) const
{
  for(unsigned i=0;i<columns;++i)
  {
    for(unsigned j=0;j<matSize;++j)
    {
      col[i*matSize+j]=matrix[j*matSize + processNum*columns + i];
    }
  }
}

void LUDecomposition::getRows(const double *matrix, double *row, unsigned matSize, int processNum, unsigned rows) const
{
  for(unsigned i=0;i<rows;++i)
  {
    for(unsigned j=0;j<matSize;++j)
    {
      row[i*matSize+j]=matrix[(processNum*rows+i)*matSize + j];
    }
  }
}

Matrix LUDecomposition::getL() const
{
  return l;
}

Matrix LUDecomposition::getU() const
{
  return u;
}

Matrix LUDecomposition::getLU() const
{
  Matrix lu(u);
  for (unsigned i = 0; i < lu.size() - 1; ++i)
  {
    for (unsigned j = i + 1; j < lu.size(); ++j)
    {
      lu[j][i] = l[j][i];
    }
  }
  return lu;
}
