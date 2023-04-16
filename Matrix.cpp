#include "Matrix.h"

Matrix::Matrix(unsigned newSize)
:matSize(newSize)
{
  matrix = matrix_type(matSize);
  for(unsigned i=0;i<matSize;++i)
  {
    matrix[i] = row_type(matSize);
    for(unsigned j=0;j<matSize;++j)
    {
      matrix[i][j] = 0.0;
    }
  }
}
Matrix::Matrix(const Matrix& other)
:matSize(other.matSize)
{
  matrix = matrix_type(matSize);
  for(unsigned i=0;i<matSize;++i)
  {
    matrix[i] = row_type(matSize);
    for(unsigned j=0;j<matSize;++j)
    {
      matrix[i][j] = other.matrix[i][j];
    }
  }
}

Matrix::Matrix(Matrix&& other)
:matSize(other.matSize), matrix(std::move(other.matrix))
{

}

Matrix::Matrix(const matrix_type& array2d)
:matSize(array2d.size())
{
  matrix = matrix_type(matSize);
  for(const auto& vec:array2d)
  {
    matrix.push_back(row_type(matSize));
    for(const auto& val:vec)
    {
      matrix.back().push_back(val);
    }
  }
}
Matrix::Matrix(std::ifstream& file, unsigned newSize)
:matSize(newSize)
{
  matrix = matrix_type(matSize);
  for(unsigned i=0;i<matSize;++i)
  {
    matrix[i] = row_type(matSize);
    for(unsigned j=0;j<matSize;++j)
    {
      file >> matrix[i][j];
    }
  }
}

Matrix::~Matrix()
{

}
unsigned Matrix::size() const
{
  return matSize;
}

const row_type Matrix::operator[](unsigned pos) const
{
  return matrix[pos];
}

row_type& Matrix::operator[](unsigned pos)
{
  return matrix[pos];
}

Matrix Matrix::operator*(const Matrix& other)
{
  Matrix result(matSize);
  for(unsigned i=0;i<matSize;++i)
  {
    for(unsigned j=0;j<matSize;++j)
    {
      for(unsigned k=0;k<matSize;++k)
      {
        result[i][j]+=matrix[i][k]*other[k][j];
      }
    }
  }
  return result;
}

Matrix Matrix::operator+(const Matrix& other)
{
  Matrix result(matSize);
  for(unsigned i=0;i<matSize;++i)
  {
    for(unsigned j=0;j<matSize;++j)
    {
      result[i][j]=matrix[i][j]+other[i][j];
    }
  }
  return result;
}

Matrix Matrix::operator-(const Matrix& other)
{
  Matrix result(matSize);
  for(unsigned i=0;i<matSize;++i)
  {
    for(unsigned j=0;j<matSize;++j)
    {
      result[i][j]=matrix[i][j]-other[i][j];
    }
  }
  return result;
}

Matrix Matrix::operator+(double value)
{
  Matrix result(matSize);
  for(unsigned i=0;i<matSize;++i)
  {
    for(unsigned j=0;j<matSize;++j)
    {
      result[i][j]=matrix[i][j]+value;
    }
  }
  return result;
}

Matrix Matrix::operator-(double value)
{
  Matrix result(matSize);
  for(unsigned i=0;i<matSize;++i)
  {
    for(unsigned j=0;j<matSize;++j)
    {
      result[i][j]=matrix[i][j]-value;
    }
  }
  return result;
}

const Matrix& Matrix::operator=(const Matrix& other)
{
  matSize = other.matSize;
  matrix = other.matrix;
  return *this;
}

Matrix Matrix::eye(unsigned newSize)
{
  Matrix mat(newSize);
  for(unsigned i=0; i<newSize; ++i)
  {
      mat[i][i]=1.0;
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
  for(auto& vec:mat.matrix)
    for(auto& val:vec)
      val = 1.0;
  return mat;
}

std::ostream& operator<<(std::ostream& o, const Matrix& matrix)
{
  for(unsigned i=0;i<matrix.size();++i)
  {
    o<<"| ";
    for(unsigned j=0;j<matrix.size();++j)
    {
      printf("% 8.3f\t| ", matrix[i][j]);
    }
    o<<'\n';
  }
  return o;
}
