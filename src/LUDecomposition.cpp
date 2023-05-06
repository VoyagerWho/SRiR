#include "LUDecomposition.h"

LUDecomposition::LUDecomposition(const double *orgA, const unsigned matSize, const int myid, const int numOfProcs)
    : u(Matrix(matSize)), l(Matrix(matSize))
{
  // wyznaczenie maksymalnej liczby wierszy/kolumn przydzielonych do jednego węzła
  const unsigned rows = matSize / numOfProcs + 1;
  // alokacja pamięci podręcznej
  double *A = new double[matSize * matSize]();
  double *row = new double[numOfProcs * rows]();
  double *myRow = new double[rows]();
  // kopiowanie macierzy A i rozesłanie jej do wszystkich węzłów
  if (myid == 0)
  {
    double *ptr = A;
    for (const double *i = orgA; unsigned(i - orgA) < matSize * matSize; ++i)
    {
      *ptr++ = *i;
    }
  }
  MPI_Bcast(A, matSize * matSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  // iteracja po zależnym wymiarze macierzy
  // wiersze dla macierzy L
  // kolumny dla macierzy U
  for (unsigned i = 0; i < matSize; ++i)
  {
    // wiadomość kontrolna informująca o stanie działania
    if (myid == 0 && !(i % 100))
      std::cout << i << "\n";

    // iteracja po niezależnych kolumnach macierzy L
    for (unsigned j = myid * rows; (j < (myid + 1) * rows) && (j < matSize); ++j)
    {
      if (j < i)
      {
        // zerowanie elementów ponad diagonalą
        l[j][i] = 0;
      }
      else
      {
        // wyznaczenie współczynników pod diagonalą
        l[j][i] = A[j * matSize + i];
        for (unsigned k = 0; k < i; ++k)
        {
          l[j][i] = l[j][i] - l[j][k] * u[k][i];
        }
      }
      // zapis zmodyfikowanych wartości przez węzeł
      myRow[j - myid * rows] = l[j][i];
    }
    // zebranie zmodyfikowanej kolumny macierzy L
    MPI_Gather(myRow, rows, MPI_DOUBLE, row, rows, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    // aktualizacja wartości kolumny macierzy L
    MPI_Bcast(row, matSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    for (unsigned j = 0; j < matSize; ++j)
      l[j][i] = row[j];

    // iteracja po niezależnych wierszach macierzy U
    for (unsigned j = myid * rows; (j < (myid + 1) * rows) && (j < matSize); ++j)
    {
      if (j < i)
      {
        // zerowanie wartości pod diagonalą
        u[i][j] = 0;
      }
      else if (j == i)
      {
        // ustawienie jedynki na diagonali
        u[i][j] = 1;
      }
      else
      {
        // wyznaczenie współczynników nad diagonalą
        u[i][j] = A[i * matSize + j] / l[i][i];
        for (unsigned k = 0; k < i; ++k)
        {
          double multi = l[i][k] * u[k][j];
          u[i][j] = u[i][j] - (multi / l[i][i]);
        }
      }
      // zapis zmodyfikowanych wartości przez węzeł
      myRow[j - myid * rows] = u[i][j];
    }
    // zebranie zmodyfikowanego wiersza macierzy U
    MPI_Gather(myRow, rows, MPI_DOUBLE, row, rows, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    // aktualizacja wartości wiersza macierzy U
    MPI_Bcast(row, matSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    for (unsigned j = 0; j < matSize; ++j)
      u[i][j] = row[j];
  }
  // dealokacja pamięci podręcznej
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
    for (unsigned j = i; j < lu.size(); ++j)
    {
      lu[j][i] = l[j][i];
    }
  }
  return lu;
}
