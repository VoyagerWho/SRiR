#ifndef LUDECOMPOSITION_H
#define LUDECOMPOSITION_H
#include "Matrix.h"
#include <upcxx/upcxx.hpp>

/**
 * @brief Klasa reprezentująca algorytm rozkładu LU dla macierzy kwadratowych
 */
class LUDecomposition
{

public:
  /**
   * @brief Konstruktor wykonujący rozkład LU macierzy poprzez klaster MPI
   *
   * @param org Stały wskaźnik do rozkładanej macierzy
   * @param matSize Rozmiar rozkładanej macierzy
   * @param myid Identyfikator węzła wewnątrz komunikatora MPI_COMM_WORLD
   * @param numOfProcs Łączna liczba węzłów wewnątrz komunikatora MPI_COMM_WORLD
   */
  LUDecomposition(const double *orgA, const unsigned matSize, const int myid, const int numOfProcs);

  /**
   * @brief Destruktor
   */
  ~LUDecomposition();

  /**
   * @brief Funkcja zwracająca kopię macierzy L
   *
   * @return Kopia macierzy L
   */
  Matrix getL() const;

  /**
   * @brief Funkcja zwracająca kopię macierzy U
   *
   * @return Kopia macierzy U
   */
  Matrix getU() const;

  /**
   * @brief Funkcja zwracająca złożenie macierzy L i U w jedną
   *
   * @return Złożenie L i U
   */
  Matrix getLU() const;

protected:
private:
  Matrix u_done; //!< Lokalna macierz U
  Matrix l_done; //!< Lokalna macierz L
};

#endif // LUDECOMPOSITION_H
