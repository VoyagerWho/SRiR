#include "LUDecomposition.h"

LUDecomposition::LUDecomposition(const double *orgA, const unsigned matSize, const int myid, const int numOfProcs)
    : u_done(Matrix(matSize)), l_done(Matrix(matSize))
{

  // stworzenie macierzy wejściowej i przypisanie jej wartości
  upcxx::dist_object<upcxx::global_ptr<double>> A_g(upcxx::new_array<double>(matSize*matSize));
  double *A = A_g->local();
  if (upcxx::rank_me() == 0) {
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

  //stworzenie rozłożonych po węzłach obiektów przechowujących lokalne macierze l i u
  upcxx::dist_object<upcxx::global_ptr<double>> l_g(upcxx::new_array<double>(matSize*matSize));
  upcxx::dist_object<upcxx::global_ptr<double>> u_g(upcxx::new_array<double>(matSize*matSize));
  double *l = l_g->local();
  double *u = u_g->local();
  //wskaźniki na obiektu l i u węzła 0 do zapisywania wyliczeń
  upcxx::global_ptr<double> l0 = l_g.fetch(0).wait();
  upcxx::global_ptr<double> u0 = u_g.fetch(0).wait();

  // wyznaczenie maksymalnej liczby wierszy/kolumn przydzielonych do jednego węzła
  const unsigned rows = matSize / numOfProcs + 1;

  //zmienne pomocnicze
  double lii = 0;
  unsigned j = 0;
  std::vector<double> buffer;
  buffer.reserve(matSize);

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
        l[j * matSize + i] = 0; 
      }
      else
      {
        // wyznaczenie współczynników pod diagonalą
        l[j * matSize + i] = A[j * matSize + i];
        for (unsigned k = 0; k < i; ++k)
        {
          l[j * matSize + i] = l[j * matSize + i] - l[j * matSize + k] * u[k * matSize + i];
        }
      }
      //wysłanie wyniku węzłowi 0
      upcxx::rput(l[j * matSize + i], l0 + j*matSize + i).wait();
    }

    //zebranie przed chwilą wyliczonej kolumny macierzy l przez węzeł 0 o rozesłanie jej wszystkim węzłom
    upcxx::barrier();	
    if (upcxx::rank_me() == 0) {
      for (size_t j = i; j < matSize; j++)
      {
        buffer[j] = l[j * matSize + i];
      }
    }
    upcxx::future<> fut_l = upcxx::broadcast( buffer.data(), matSize, 0); 
    fut_l.wait();
    if(upcxx::rank_me() != 0)
    {
      for (size_t j = i; j < matSize; j++){
        l[j * matSize + i] = buffer[j];
      }
    }
    upcxx::barrier();

    // iteracja po niezależnych wierszach macierzy U
    for (j = myid * rows; (j < (myid + 1) * rows) && (j < matSize); ++j)
    {
      if (j < i)
      {
        // zerowanie wartości pod diagonalą
        u[i * matSize + j] = 0;
      }
      else if (j == i)
      {
        // ustawienie jedynki na diagonali
        u[i * matSize + j] = 1;
      }
      else
      {
        // wyznaczenie współczynników nad diagonalą
        lii = l[i*matSize + i];
        if(std::abs(lii) > 1e-32){
          u[i * matSize + j] = A[i * matSize + j] / lii;
          for (unsigned k = 0; k < i; ++k)
          {
            double multi = l[i * matSize + k] * u[k * matSize + j];
            u[i * matSize + j] = u[i * matSize + j] - (multi / lii);
          }
        }
      }
      //wysłanie wyniku węzłowi 0
      upcxx::rput( u[i * matSize + j], u0 + i*matSize + j).wait();
    }
    
    //zebranie przed chwilą wyliczonego wiersza macierzy u przez węzeł 0 o rozesłanie jej wszystkim węzłom
    upcxx::barrier();	
    if (upcxx::rank_me() == 0) {
      for (size_t j = i; j < matSize; j++)
      {
        buffer[j] = u[i * matSize + j];
      }
    }
    upcxx::future<> fut_u = upcxx::broadcast( buffer.data(), matSize, 0); 
    fut_u.wait();
    if(upcxx::rank_me() != 0)
    {
      for (size_t j = i; j < matSize; j++){
        u[i * matSize + j] = buffer[j];
      }
    }
    upcxx::barrier();	
  }

  //przepisanie przez węzel 0 wyliczonych macierzy do lokalnej pamięci
  if (upcxx::rank_me() == 0) {
    for(unsigned i = 0; i < matSize; i++)
    {
      for(unsigned j = 0; j < matSize; j++)
      {
        l_done[i][j] = l[(i * matSize) + j];
        u_done[i][j] = u[(i * matSize) + j];
      }
    }
  }
  // dealokacja pamięci
  upcxx::delete_(*A_g);
  upcxx::delete_(*l_g);
  upcxx::delete_(*u_g);
  upcxx::barrier();
}

LUDecomposition::~LUDecomposition()
{
  // dtor
}

Matrix LUDecomposition::getL() const
{
  return l_done;
}

Matrix LUDecomposition::getU() const
{
  return u_done;
}

Matrix LUDecomposition::getLU() const
{
  Matrix lu(u_done);
  for (unsigned i = 0; i < lu.size() - 1; ++i)
  {
    for (unsigned j = i; j < lu.size(); ++j)
    {
      lu[j][i] = l_done[j][i];
    }
  }
  return lu;
}
