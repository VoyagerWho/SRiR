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

LUDecomposition::LUDecomposition(const double *orgA, const unsigned matSize, const int myid, const int numOfProcs)
    : u(Matrix(matSize)), l(Matrix(matSize))
{
  const unsigned rows = matSize / numOfProcs;
  double *A = new double[matSize * matSize]();
  double *row = new double[numOfProcs * rows]();
  double *myRow = new double[rows]();
  if (myid == 0)
  {
    double *ptr = A;
    for (const double *i = orgA; unsigned(i - orgA) < matSize * matSize; ++i)
    {
      *ptr++ = *i;
    }
  }
  MPI_Bcast(A, matSize * matSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  for (unsigned i = 0; i < matSize; ++i)
  {
    if (myid == 0 && !(i % 100))
      std::cout << i << "\n";
    for (unsigned j = myid * rows; j < (myid + 1) * rows; ++j)
    {
      if (j < i)
      {
        l[j][i] = 0;
      }
      else
      {
        l[j][i] = A[j * matSize + i];
        for (unsigned k = 0; k < i; ++k)
        {
          l[j][i] = l[j][i] - l[j][k] * u[k][i];
        }
      }
      myRow[j - myid * rows] = l[j][i];
    }
    MPI_Gather(myRow, rows, MPI_DOUBLE, row, rows, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(row, matSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    for (unsigned j = 0; j < matSize; ++j)
      l[j][i] = row[j];
    for (unsigned j = myid * rows; j < (myid + 1) * rows; ++j)
    {
      if (j < i)
      {
        u[i][j] = 0;
      }
      else if (j == i)
      {
        u[i][j] = 1;
      }
      else
      {
        u[i][j] = A[i * matSize + j] / l[i][i];
        for (unsigned k = 0; k < i; ++k)
        {
          double multi = l[i][k] * u[k][j];
          u[i][j] = u[i][j] - (multi / l[i][i]);
        }
      }
      myRow[j - myid * rows] = u[i][j];
    }
    MPI_Gather(myRow, rows, MPI_DOUBLE, row, rows, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(row, matSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    for (unsigned j = 0; j < matSize; ++j)
      u[i][j] = row[j];
  }

  delete[] A;
  delete[] row;
  delete[] myRow;
}

LUDecomposition::~LUDecomposition()
{
  // dtor
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
