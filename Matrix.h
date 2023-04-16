#ifndef MATRIX_H
#define MATRIX_H
#include <iostream>
#include <fstream>
#include <functional>
#include <vector>

using row_type = std::vector<double>;
using matrix_type = std::vector<row_type>;

class Matrix
{
  public:
    Matrix(unsigned newSize);
    Matrix(const Matrix& other);
    Matrix(Matrix&& other);
    Matrix(const matrix_type& array2d);
    Matrix(std::ifstream& file, unsigned newSize);
    ~Matrix();
    unsigned size() const;

    const row_type operator[](unsigned pos) const;
    row_type& operator[](unsigned pos);
    Matrix operator*(const Matrix& other);
    Matrix operator+(const Matrix& other);
    Matrix operator+(double);
    Matrix operator-(const Matrix& other);
    Matrix operator-(double);
    const Matrix& operator=(const Matrix& other);

    static Matrix eye(unsigned newSize);
    static Matrix zeros(unsigned newSize);
    static Matrix ones(unsigned newSize);
  protected:
  private:
    unsigned matSize;
    matrix_type matrix; //!< Member variable "matrix"
};

std::ostream& operator<<(std::ostream& o, const Matrix& matrix);

#endif // MATRIX_H
