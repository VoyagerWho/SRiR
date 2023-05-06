#include <math.h>
#include "SolutionCalculator.h"

SolutionCalculator::SolutionCalculator(const LUDecomposition& _lu, std::ifstream* _bVecFile, 
                                        const int _myid, const int _numProcs, const unsigned _n):
                                        L(_lu.getL()), U(_lu.getU()), myid(_myid), 
                                        numProcs(_numProcs), n(_n)
{
    assignColumnsToProcess();

    // alokacja pamięci dla wektorów x, y i b
    x = new double[n]();
    y = new double[n]();
    b = new double[n]();

    // proces master kopiuje wektor b do zaalokowanej pamięci
    if(myid == 0)
    {
        double *ptr = b;
        for (unsigned i = 0; i < n; ++i)
        {
            *_bVecFile >> *ptr++;
        }
    }

    // rozesłanie wektora b reszcie procesów
    MPI_Bcast(b, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

SolutionCalculator::~SolutionCalculator()
{
    // zwolnienie zaalokowanej pamięci
    delete[] y;
    delete[] x;
    delete[] b;
}

void SolutionCalculator::printSolutionVector() const
{
    // pętla wypisująca każdy element wektora rozwiązań x
    for(int i = 0; i < (int)n; i++)
    {
        printf("% 8.3f\n ", x[i]);
    }
}

void SolutionCalculator::run()
{
    // wyliczenie kolejno wektora y oraz x
    calculateYVector();
    calculateXVector();
}

double SolutionCalculator::mse(std::ifstream* targetVec)
{
    // alokacja pamięci na wektor do którego porównywany będzie wynik programu
    double* target = new double[n]();
    double *ptr = target;
    // przapisanie wartości z pliku do wektora porównawczego
    for (unsigned i = 0; i < n; ++i)
    {
        *targetVec >> *ptr++;
    }

    // wyliczenie błędy średniokwadratowego
	double error = 0;
	for (int i = 0; i < (int)n; i++) {
		error += pow((target[i] - x[i]), 2);
	}

    delete[] target;
	return error / n;
}

//-----------------------------------
// Private methods definitions
//-----------------------------------

void SolutionCalculator::calculateYVector()
{
    // inicjalizacja zmiennych wykorzystywanych do policzenia wektora y
    double* yPartial = new double[n]();
    bool ownCurrentRow;
    double partialRowSum;
    double globalSum;
    // pętla po każdym wierszu
    for(int row = 0; row < (int)n; row++)
    {
        //wypisanie co setny wiersz do konsoli przez proces master aktualnego stanu wyliczania y
        if (myid == 0 && !(row % 100))
            std::cout << "y: " << row << "\n";
        // wyzerowanie zmiennych,
        partialRowSum = 0;
        globalSum = 0;
        ownCurrentRow = false;

        // przejście po kolumnach zapisanych w mycols
        for(unsigned collumnOwned = 0; collumnOwned < maxNumberOfCollumnsOwned; collumnOwned++)
        {
            // jeżeli indeks kolumny jest mniejszy niż indeks kolumny z aktualnie wylicznym y, proces wylicza swoją część rozwiązania
            if((mycols[collumnOwned] < row) && mycols[collumnOwned] != -1)
            {
                partialRowSum += L[row][mycols[collumnOwned]] * yPartial[mycols[collumnOwned]];
            }
            // jeżeli posiadana kolumna jest tą która powinna wyliczyć aktualny element y, zmiana stanu flagi
            if(mycols[collumnOwned] == row)
                ownCurrentRow = true;
        }
        
        // zebranie od procesów częściowych sum wiersza i rozesłanie sumy całkowitej do każdego z procesów
        MPI_Allreduce( &partialRowSum, &globalSum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
        // proces który posiada kolumne z aktualnie wyliczanym elementem y, wylicza go
        if(ownCurrentRow && abs(L[row][row]) > 0.000001)
        {
            yPartial[row] = (b[row] - globalSum) / L[row][row];
        }
    }

    // zebranie elemetów wektora z procesów w cały wektor y i rozesłanie go procesom
    MPI_Allreduce(yPartial, y, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    delete[] yPartial;
}

void SolutionCalculator::calculateXVector()
{
    // inicjalizacja zmiennych wykorzystywanych do policzenia wektora x
    double* xPartial = new double[n]();
    bool ownCurrentRow;
    double partialRowSum;
    double globalSum;

    // pętla po każdym wierszu
    for(int row = n-1; row >= 0; row--)
    {
        //wypisanie co setny wiersz do konsoli przez proces master aktualnego stanu wyliczania x
        if (myid == 0 && !(row % 100))
            std::cout << "x: " << row << "\n";
        // wyzerowanie zmiennych,
        partialRowSum = 0;
        globalSum = 0;
        ownCurrentRow = false;

        // przejście po kolumnach zapisanych w mycols
        for(unsigned collumnOwned = 0; collumnOwned < maxNumberOfCollumnsOwned; collumnOwned++)
        {
            // jeżeli indeks kolumny jest większy niż indeks kolumny z aktualnie wylicznym x, proces wylicza swoją część rozwiązania
            if((mycols[collumnOwned] > row) && mycols[collumnOwned] != -1)
            {
                partialRowSum += U[row][mycols[collumnOwned]] * xPartial[mycols[collumnOwned]];
            }
            // jeżeli posiadana kolumna jest tą która powinna wyliczyć aktualny element x, zmiana stanu flagi
            if(mycols[collumnOwned] == row)
                ownCurrentRow = true;
        }
        
        // zebranie od procesów częściowych sum wiersza i rozesłanie sumy całkowitej do każdego z procesów
        MPI_Allreduce( &partialRowSum, &globalSum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
        // proces który posiada kolumne z aktualnie wyliczanym elementem x, wylicza go
        if(ownCurrentRow && abs(L[row][row]) > 0.000001)
        {
            xPartial[row] = (y[row] - globalSum) / U[row][row];
        }
    }

    // zebranie elemetów wektora z procesów w cały wektor x i wysłanie go procesowi master
    MPI_Reduce(xPartial, x, n, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    delete[] xPartial;
}

//Cyclic assignement of columns
void SolutionCalculator::assignColumnsToProcess()
{
    // wyliczenie maksymalnej ilość column jaką może posiadać pojedyńczy proces
    maxNumberOfCollumnsOwned = ceil((double)n/numProcs);
    // alokacja pamieci na tablice z indeksami kolumn którymi zajmie sie proces
    mycols = new int[maxNumberOfCollumnsOwned]();
    // inicjalizacja z wartością -1 na wypadek przypisania procesowi mniejszej
    // ilości kolumn niż jest na to miejsce w tablicy
    for(unsigned i = 0; i < maxNumberOfCollumnsOwned; i++)
    {
        mycols[i] = -1;
    }

    // cykliczne przypisanie kolumn procesom w celu bardziej równomiernego
    // rozłożenia ilości obliczeń
    int numberOfCollumnsAssigned = 0;
    for(unsigned col = 0; col < n; col++)
    {
        if(col == (unsigned)(numberOfCollumnsAssigned*numProcs + myid))
        {
            mycols[numberOfCollumnsAssigned] = col;
            numberOfCollumnsAssigned++;
        }
    }
}