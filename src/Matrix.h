#ifndef MATRIX_H
#define MATRIX_H
#include <iostream>
#include <fstream>
#include <functional>

/**
 * @brief Klasa reprezentująca macierz kwadratową nxn
 */
class Matrix
{

public:
  /**
   * @brief Konstruktor macierzy o zadanym rozmiarze
   * Macierz wypełniona zerami
   * @param newSize Rozmiar nowej macierzy
   */
  Matrix(unsigned newSize);
  /**
   * @brief Konstruktor kopiujący
   *
   * @param other Macierz do skopiowania.
   */
  Matrix(const Matrix &other);

  /**
   * @brief Konstruktor przenoszący
   *
   * @param other Macierz do przeniesienia
   */
  Matrix(Matrix &&other);

  /**
   * @brief Konstruktor macierzy o zadanym rozmiarze i wartościach początkowych
   *
   * @param rawArray Wartości macierzy zapisane wierszami
   * @param newSize Rozmiar macierzy
   */
  Matrix(const double *rawArray, unsigned newSize);

  /**
   * @brief  Konstruktor macierzy o zadanym rozmiarze i wartościach początkowych z pliku tekstowego
   *
   * @param file Plik tekstowy z wartościami macierzy zapisanymi wierszami
   * @param newSize Rozmiar macierzy
   */
  Matrix(std::ifstream &file, unsigned newSize);

  /**
   * @brief Dstruktor ze zwolnieniem pamięci
   */
  ~Matrix();

  /**
   * @brief Funkcja do pobrania rozmiaru macierzy
   *
   * @return Rozmiar macierzy
   */
  unsigned size() const;

  /**
   * @brief Operator dostępu do wierszy macierzy (do odczytu)
   *
   * @param pos Pozycja wiersza
   * @return Stały wskaznić na wiersz
   */
  const double *operator[](unsigned pos) const;

  /**
   * @brief Operator dostępu do wierszy macierzy (do zapisu)
   * Pozwala na stosowanie notacji mat[][]
   * @param pos Pozycja wiersza
   * @return Wskaznić na wiersz
   */
  double *operator[](unsigned pos);

  /**
   * @brief Operator mnożenia macierzy
   * Wolny i nie wydajny - używać rozważnie
   * @param other Domnażana macierz
   * @return Macierz wynikowa
   */
  Matrix operator*(const Matrix &other);

  /**
   * @brief Operator sumy macierzy po elementach
   *
   * @param other Macierz do dodania
   * @return Macierz wynikowa
   */
  Matrix operator+(const Matrix &other);

  /**
   * @brief Operator dodania biasu do wszystkich elementów macierzy
   *
   * @param value Wartość biasu
   * @return Macierz wynikowa
   */
  Matrix operator+(double value);

  /**
   * @brief Operator różnicy macierzy po elementach
   *
   * @param other Odejmowana macierz
   * @return Macierz wynikowa
   */
  Matrix operator-(const Matrix &other);

  /**
   * @brief Operator odejmowania biasu od wszystkich elementów macierzy
   *
   * @param value Wartość biasu
   * @return Macierz wynikowa
   */
  Matrix operator-(double value);

  /**
   * @brief Kopiujący operator przypisania wartości
   *
   * @param other Kopiowana macierz
   * @return Referencja to bierzącej macierzy
   */
  const Matrix &operator=(const Matrix &other);

  /**
   * @brief Funkcja statyczna do tworzenia macierzy jednostkowej
   *
   * @param newSize Rozmiar nowej macierzy
   * @return Nowa macierz
   */
  static Matrix eye(unsigned newSize);

  /**
   * @brief Funkcja statyczna do tworzenia macierzy zer
   *
   * @param newSize Rozmiar nowej macierzy
   * @return Nowa macierz
   */
  static Matrix zeros(unsigned newSize);

  /**
   * @brief Funkcja statyczna do tworzenia macierzy jedynek
   *
   * @param newSize Rozmiar nowej macierzy
   * @return Nowa macierz
   */
  static Matrix ones(unsigned newSize);

protected:
private:
  unsigned matSize;         //!< Rozmiar macierzy
  double *matrix = nullptr; //!< Adres zaalokowanej pamięci
};

/**
 * @brief Wypisanie macierzy z prostym formatowaniem
 *
 * @param o Referencja do strumienia wyjściowego
 * @param matrix Referencja do wyświetlanej macierzy
 * @return Referencje do otrzymanego strumienia wyjściowego
 */
std::ostream &operator<<(std::ostream &o, const Matrix &matrix);

/**
 * @brief Wyznaczenie błędu średniokwadratowego między dwoma macierzami
 *
 * @param a Pierwsza macierz
 * @param b Druga macierz
 * @return Błąd średniokwadratowy operacji a-b;
 */
double mse(const Matrix &a, const Matrix &b);

#endif // MATRIX_H
