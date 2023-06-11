#ifndef SOLUTIONCALCULATOR_H
#define SOLUTIONCALCULATOR_H

#include "LUDecomposition.h"
#include <fstream>
#include <upcxx/upcxx.hpp>

/**
 * @brief Klasa rozwiązująca układ równań liniowych Ax = b z wykorzystaniem rozkładu LU
 */
class SolutionCalculator
{
public:
    /**
     * @brief Konstruktor przygotowujący węzły klastera MPI
     *
     * @param _lu Stała referencja do rozłożonej macierzy
     * @param _bVecFile Wskaźnik do pliku z wektorem b
     * @param _myid Identyfikator węzła wewnątrz komunikatora MPI_COMM_WORLD
     * @param _numOfProcs Łączna liczba węzłów wewnątrz komunikatora MPI_COMM_WORLD
     * @param _n Wymiarowość problemu
     */
    SolutionCalculator(const LUDecomposition &_lu, std::ifstream *_bVecFile, const int _myid, const int _numProcs, const unsigned _n);

    /**
     * @brief Destruktor zwalniający pamięć
     */
    ~SolutionCalculator();

    /**
     * @brief Procedura przeprowadzająca kolejne etapy rozwiązania układu równań
     *
     */
    void run();

    /**
     * @brief Procedura wyświetlająca wektor x z podstawowym formatowaniem
     *
     */
    void printSolutionVector() const;

    /**
     * @brief Funkcja wyznaczająca błąd średniokwadratowy od oczekiwanego wyniku
     *
     * @param targetVec Wskaźnik do pliku z oczekiwanym wektorem x
     * @return Wartość błędu średniokwadratowego
     */
    double mse(std::ifstream *targetVec);

private:
    /**
     * @brief Procedura rozwiązująca układ Ly = b
     */
    void calculateYVector();

    /**
     * @brief Procedura rozwiązująca układ Ux = y
     */
    void calculateXVector();

    /**
     * @brief Procedura rozdzielająca kolumny rozwiązania do procesu
     */
    void assignColumnsToProcess();

    Matrix L;                          //!< Lokalna macierz L
    Matrix U;                          //!< Lokalna macierz U
    double *y;                         //!< Lokalna przestrzeń na wektor y
    double *x;                         //!< Lokalna przestrzeń na wektor x
    double *b;                         //!< Lokalna przestrzeń na wektor b
    const int myid;                    //!< Identyfikator węzła wewnątrz komunikatora MPI_COMM_WORLD
    const int numProcs;                //!< Łączna liczba węzłów wewnątrz komunikatora MPI_COMM_WORLD
    const unsigned n;                  //!< Wymiarowość problemu (Rozmiar macierzy)
    unsigned maxNumberOfCollumnsOwned; //!< Maksymalna liczba kolumn przydzielonych węzłowi
    int *mycols;                       //!< Lista przydzielonych kolumn
};

#endif // SOLUTIONCALCULATOR_H