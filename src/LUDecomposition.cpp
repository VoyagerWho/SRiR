#include "LUDecomposition.h"

LUDecomposition::LUDecomposition(const double *orgA, const unsigned matSize, const int myid, const int numOfProcs)
    : u_done(Matrix(matSize)), l_done(Matrix(matSize))
{
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
  upcxx::future<> fut_A = upcxx::broadcast(A, matSize * matSize, 0);
  fut_A.wait();
  upcxx::barrier();	

  upcxx::dist_object<upcxx::global_ptr<double>> l_g(upcxx::new_array<double>(matSize*matSize));
  upcxx::dist_object<upcxx::global_ptr<double>> u_g(upcxx::new_array<double>(matSize*matSize));
  double *l = l_g->local();
  double *u = u_g->local();
  upcxx::global_ptr<double> l0 = l_g.fetch(0).wait();
  upcxx::global_ptr<double> u0 = u_g.fetch(0).wait();

  // wyznaczenie maksymalnej liczby wierszy/kolumn przydzielonych do jednego węzła
  const unsigned rows = matSize / numOfProcs ;
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
          // std::cout << "[ " << myid << "] j: " << j << " k: " << k << " u: " << l[j*matSize + k] << std::endl;
        }
      }
      upcxx::rput(l[j * matSize + i], l0 + j*matSize + i).wait();
    }

    //launch a bulk broadcast of element data from rank 0
    upcxx::barrier();	
    if (upcxx::rank_me() == 0) {
      for (size_t j = i; j < matSize; j++)
      {
        // std::cout << "[" << myid << "] " << l[j * matSize + i] << std::endl;
        buffer[j] = l[j * matSize + i];
      }
    }
    // launch a bulk broadcast of element data from rank 0
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
          // std::cout << "[ " << myid << "] i: " << i << " j: " << j << " l: " << l[i*matSize + i] << std::endl;
          for (unsigned k = 0; k < i; ++k)
          {
            double multi = l[i * matSize + k] * u[k * matSize + j];
            u[i * matSize + j] = u[i * matSize + j] - (multi / lii);
            // std::cout << "[ " << myid << "] k: " << k << " j: " << j << " l: " << u[k*matSize + j] << std::endl;
          }
        }
      }
      // std::cout << "[" << myid << "] " << u[i * matSize + j] << std::endl;
      upcxx::rput( u[i * matSize + j], u0 + i*matSize + j).wait();
    }
    
    upcxx::barrier();	
    if (upcxx::rank_me() == 0) {
      for (size_t j = i; j < matSize; j++)
      {
        // std::cout << "[" << myid << "] " << u[i * matSize + j] << std::endl;
        buffer[j] = u[i * matSize + j];
      }
    }
    // launch a bulk broadcast of element data from rank 0
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

  // dealokacja pamięci
  if (upcxx::rank_me() == 0) {
    for(unsigned i = 0; i < matSize; i++)
    {
      for(unsigned j = 0; j < matSize; j++)
      {
        // std::cout << l[(i*matSize + j)] << " ";
        l_done[i][j] = l[(i * matSize) + j];
        u_done[i][j] = u[(i * matSize) + j];
      }
    }
  }
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
