#ifndef MATRIX_H
#define MATRIX_H
#include <iostream>
#include <fstream>
#include <functional>

class Matrix
{
public:
  Matrix(unsigned newSize);
  Matrix(const Matrix &other);
  Matrix(Matrix &&other);
  Matrix(const double *rawArray, unsigned newSize);
  Matrix(std::ifstream &file, unsigned newSize);
  ~Matrix();
  unsigned size() const;

  const double *operator[](unsigned pos) const;
  double *operator[](unsigned pos);
  Matrix operator*(const Matrix &other);
  Matrix operator+(const Matrix &other);
  Matrix operator+(double);
  Matrix operator-(const Matrix &other);
  Matrix operator-(double);
  const Matrix &operator=(const Matrix &other);

  static Matrix eye(unsigned newSize);
  static Matrix zeros(unsigned newSize);
  static Matrix ones(unsigned newSize);

protected:
private:
  unsigned matSize;
  double *matrix = nullptr; //!< Member variable "matrix"
};

std::ostream &operator<<(std::ostream &o, const Matrix &matrix);

double mse(const Matrix &a, const Matrix &b);

#endif // MATRIX_H
