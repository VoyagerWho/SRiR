#include "LUDecomposition.h"

LUDecomposition::LUDecomposition(const double *orgA, const unsigned matSize, const int myid, const int numOfProcs)
    : u(Matrix(matSize)), l(Matrix(matSize))
{

  // stworzenie macierzy wejściowej i przypisanie jej wartości
  upcxx::dist_object<upcxx::global_ptr<double>> A_g(upcxx::new_array<double>(matSize*matSize));
  double *A = A_g->local();
  upcxx::global_ptr<double> l0 = nullptr;
  upcxx::global_ptr<double> u0 = nullptr;
  if (upcxx::rank_me() == 0) {
    l0 = upcxx::new_array<double>(matSize*matSize);
    u0 = upcxx::new_array<double>(matSize*matSize);
    for(unsigned i = 0; i < matSize; i++)
    {
      for(unsigned j = 0; j < matSize; j++)
      {
        A[(i * matSize) + j] = orgA[(i * matSize) + j];
      }
    }
  }
  upcxx::barrier();	


  //rozesłanie macierzy wejściowej do wszystkich węzłów
  upcxx::future<> fut_A = upcxx::broadcast(A, matSize * matSize, 0);
  fut_A.wait();
  upcxx::barrier();	

  //wskaźniki na obiektu l i u węzła 0 do zapisywania wyliczeń
  l0 = upcxx::broadcast(l0, 0).wait();
  u0 = upcxx::broadcast(u0, 0).wait();

  // wyznaczenie maksymalnej liczby wierszy/kolumn przydzielonych do jednego węzła
  const unsigned rows = matSize / numOfProcs + 1;

  //zmienne pomocnicze
  double lii = 0;
  unsigned j = 0;
  std::vector<double> buffer;
  buffer.reserve(matSize);
  std::vector<upcxx::future<>> fut_l0;
  std::vector<upcxx::future<>> fut_u0;

  // iteracja po zależnym wymiarze macierzy
  // wiersze dla macierzy L
  // kolumny dla macierzy U
  for (unsigned i = 0; i < matSize; ++i)
  {
    // wiadomość kontrolna informująca o stanie działania
    if (myid == 0 && !(i % 100))
      std::cout << i << "\n";

    // iteracja po niezależnych kolumnach macierzy L
    for (j = myid * rows; (j < (myid + 1) * rows) && (j < matSize); ++j)
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
      //asynchroniczne wysłanie wyniku węzłowi 0
      fut_l0.push_back(upcxx::rput(l[j][i], l0 + j*matSize + i));
    }
    //upewnienie się, że wszystkie wyniki dotarły i zostały zapisane
    for(auto future: fut_l0)
      future.wait();
    fut_l0.clear();

    //zebranie przed chwilą wyliczonej kolumny macierzy l przez węzeł 0 o rozesłanie jej wszystkim węzłom
    upcxx::barrier();	
    if (upcxx::rank_me() == 0) {
      double* myl0 = l0.local();
      for (size_t j = i; j < matSize; j++)
      {
        buffer[j] = myl0[j * matSize + i];
      }
    }
    upcxx::future<> fut_l = upcxx::broadcast( buffer.data(), matSize, 0); 
    fut_l.wait();
    for (size_t j = i; j < matSize; j++){
      l[j][i] = buffer[j];
    }
    upcxx::barrier();

    // iteracja po niezależnych wierszach macierzy U
    for (j = myid * rows; (j < (myid + 1) * rows) && (j < matSize); ++j)
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
        lii = l[i][i];
        if(std::abs(lii) > 1e-14){
          u[i][j] = A[i * matSize + j] / lii;
          for (unsigned k = 0; k < i; ++k)
          {
            double multi = l[i][k] * u[k][j];
            u[i][j] = u[i][j] - (multi / lii);
          }
        }
      }
      //asynchroniczne wysłanie wyniku węzłowi 0
      fut_u0.push_back(upcxx::rput(u[i][j], u0 + i*matSize + j));
    }
    //upewnienie się, że wszystkie wyniki dotarły i zostały zapisane
    for(auto future: fut_u0)
      future.wait();
    fut_u0.clear();

    //zebranie przed chwilą wyliczonego wiersza macierzy u przez węzeł 0 o rozesłanie jej wszystkim węzłom
    upcxx::barrier();	
    if (upcxx::rank_me() == 0) {
      double* myu0 = u0.local();
      for (size_t j = i; j < matSize; j++)
      {
        buffer[j] = myu0[i * matSize + j];
      }
    }
    upcxx::future<> fut_u = upcxx::broadcast( buffer.data(), matSize, 0); 
    fut_u.wait();
    for (size_t j = i; j < matSize; j++){
      u[i][j] = buffer[j];
    }
    upcxx::barrier();	
  }
  
  // dealokacja pamięci
  upcxx::delete_(*A_g);
  if(upcxx::rank_me() == 0)
  {
    upcxx::delete_(l0);
    upcxx::delete_(u0);
  }
  upcxx::barrier();
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
