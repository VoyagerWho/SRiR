#include "Matrix.h"

Matrix::Matrix(unsigned newSize)
    : matSize(newSize)
{
  matrix = new double[matSize * matSize]();
}
Matrix::Matrix(const Matrix &other)
    : matSize(other.matSize)
{
  matrix = new double[matSize * matSize]();
  double *ptr = matrix;
  for (double *i = other.matrix; unsigned(i - other.matrix) < matSize * matSize; ++i)
  {
    *ptr++ = *i;
  }
}

Matrix::Matrix(Matrix &&other)
    : matSize(std::move(other.matSize)), matrix(std::move(other.matrix))
{
}

Matrix::Matrix(const double *rawArray, unsigned newSize)
    : matSize(newSize)
{
  matrix = new double[matSize * matSize]();
  if (rawArray)
  {
    double *ptr = matrix;
    for (const double *i = rawArray; unsigned(i - rawArray) < matSize * matSize; ++i)
    {
      *ptr++ = *i;
    }
  }
}
Matrix::Matrix(std::ifstream &file, unsigned newSize)
    : matSize(newSize)
{
  matrix = new double[matSize * matSize]();
  double *ptr = matrix;
  for (unsigned i = 0; i < matSize * matSize; ++i)
  {
    file >> *ptr++;
  }
}

Matrix::~Matrix()
{
  if (matrix)
    delete[] matrix;
}
unsigned Matrix::size() const
{
  return matSize;
}

const double *Matrix::operator[](unsigned pos) const
{
  return matrix + (pos * matSize);
}

double *Matrix::operator[](unsigned pos)
{
  return matrix + (pos * matSize);
}

Matrix Matrix::operator*(const Matrix &other)
{
  Matrix result(matSize);
  for (unsigned i = 0; i < matSize; ++i)
  {
    for (unsigned j = 0; j < matSize; ++j)
    {
      for (unsigned k = 0; k < matSize; ++k)
      {
        result[i][j] += (*this)[i][k] * other[k][j];
      }
    }
  }
  return result;
}

Matrix Matrix::operator+(const Matrix &other)
{
  Matrix result(matSize);
  for (unsigned i = 0; i < matSize; ++i)
  {
    for (unsigned j = 0; j < matSize; ++j)
    {
      result[i][j] = (*this)[i][j] + other[i][j];
    }
  }
  return result;
}

Matrix Matrix::operator-(const Matrix &other)
{
  Matrix result(matSize);
  for (unsigned i = 0; i < matSize; ++i)
  {
    for (unsigned j = 0; j < matSize; ++j)
    {
      result[i][j] = (*this)[i][j] - other[i][j];
    }
  }
  return result;
}

Matrix Matrix::operator+(double value)
{
  Matrix result(matSize);
  for (unsigned i = 0; i < matSize; ++i)
  {
    for (unsigned j = 0; j < matSize; ++j)
    {
      result[i][j] = (*this)[i][j] + value;
    }
  }
  return result;
}

Matrix Matrix::operator-(double value)
{
  Matrix result(matSize);
  for (unsigned i = 0; i < matSize; ++i)
  {
    for (unsigned j = 0; j < matSize; ++j)
    {
      result[i][j] = (*this)[i][j] - value;
    }
  }
  return result;
}

const Matrix &Matrix::operator=(const Matrix &other)
{
  matSize = other.matSize;
  double *ptr = matrix;
  for (double *i = other.matrix; unsigned(i - other.matrix) < matSize * matSize; ++i)
  {
    *ptr++ = *i;
  }
  return *this;
}

Matrix Matrix::eye(unsigned newSize)
{
  Matrix mat(newSize);
  for (unsigned i = 0; i < newSize; ++i)
  {
    mat[i][i] = 1.0;
  }
  return mat;
}

Matrix Matrix::zeros(unsigned newSize)
{
  return Matrix(newSize);
}

Matrix Matrix::ones(unsigned newSize)
{
  Matrix mat(newSize);
  double *ptr = mat.matrix;
  for (unsigned i = 0; i < mat.matSize * mat.matSize; ++i)
  {
    *ptr++ = 1;
  }
  return mat;
}

std::ostream &operator<<(std::ostream &o, const Matrix &matrix)
{
  for (unsigned i = 0; i < matrix.size(); ++i)
  {
    o << "| ";
    for (unsigned j = 0; j < matrix.size(); ++j)
    {
      printf("% 8.3f\t| ", matrix[i][j]);
    }
    o << '\n';
  }
  return o;
}

double mse(const Matrix &a, const Matrix &b)
{
  const unsigned n = a.size();
  const double *aptr = a[0];
  double sum = 0.0;
  for (const double *bptr = b[0]; unsigned(bptr - b[0]) < n * n; ++bptr)
  {
    sum += (*aptr - *bptr) * (*aptr - *bptr);
    ++aptr;
  }
  return sum / (n * n);
}
