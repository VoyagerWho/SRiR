#include <iostream>
#include <fstream>
#include <getopt.h>
#include <string.h>
#include "src/Matrix.h"
#include "src/LUDecomposition.h"
#include "src/SolutionCalculator.h"

using namespace std;

/**
 * @brief Struktura przechowująca argumenty wywołania
 */
struct Arguments
{
	int matrixSize = -1;
	int compareX = 0;
	char fileNameScheme[256] = {0};
	char matrixFileName[256] = {'M', 'a', 't', 'r', 'i', 'x'};
	char bFileName[256] = {'b'};
	char xFileName[256] = {'x'};
};

/**
 * @brief Funkcja wyznaczająca wartości argumentów wywołania
 *
 * @param argc Liczba argumentów
 * @param argv Lista argumentów
 * @param args Lista przetworzonych argumentów
 * @param myid Identyfikator procesu MPI
 */
void handleArguments(int argc, char **argv, Arguments &args, int myid)
{
	int opt;
	// lista argumentów podawanych przy wywołaniu programu
	struct option longopts[] =
		{
			{"matrixSize", required_argument, &args.matrixSize, 's'},
			{"fileNameScheme", required_argument, NULL, 'm'},
			{"compareX", no_argument, &args.compareX, 1},
			{0}};
	// pętla przypisująca argumenty do właściwych pól struktury Arguments
	while (true)
	{
		opt = getopt_long(argc, argv, "xs:m:", longopts, 0);
		if (opt == -1)
			break;
		switch (opt)
		{
		case 's':
			args.matrixSize = atoi(optarg);
			if (myid == 0)
				std::cout << "matrix size: " << args.matrixSize << std::endl;
			break;
		case 'm':
			strncpy(args.fileNameScheme, optarg ? optarg : "Test.txt", sizeof(args.fileNameScheme));
			args.fileNameScheme[sizeof(args.fileNameScheme) - 1] = '\0';
			strcat(args.matrixFileName, args.fileNameScheme);
			strcat(args.bFileName, args.fileNameScheme);
			strcat(args.xFileName, args.fileNameScheme);
			if (myid == 0)
			{
				std::cout << "file name scheme: " << args.fileNameScheme << std::endl;
				std::cout << "matrix file name: " << args.matrixFileName << std::endl;
				std::cout << "b file name: " << args.bFileName << std::endl;
			}
			break;
		case 'x':
			args.compareX = true;
			break;
		case '?':
			if (myid == 0)
				std::cout << "argument error!" << std::endl;
			exit(1);
		default:
			break;
		}
	}
	// w przypadku nie podania wymaganego argumentu, zamknięcie programu
	if (myid == 0)
	{
		if (args.matrixSize == -1)
		{
			std::cout << "matrix size not specified" << std::endl;
			exit(1);
		}
		if (args.fileNameScheme[0] == false)
		{
			std::cout << "file name scheme not specified" << std::endl;
			exit(1);
		}
		if (args.compareX)
		{
			std::cout << "mse x comparison active" << std::endl;
		}
	}
}

int main(int argc, char **argv)
{
	Arguments args;

	MPI_Init(&argc, &argv);
	int myid, numProcs;
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

	handleArguments(argc, argv, args, myid);

	if (myid == 0)
	{
		// czynności węzła master
		// wczytanie macierzy A
		std::ifstream file(args.matrixFileName);
		Matrix m = Matrix(file, args.matrixSize);

		// pomiar czasu wykonania rozkładu LU
		double startTime = MPI_Wtime();
		LUDecomposition lu(m[0], args.matrixSize, myid, numProcs);
		double endTime = MPI_Wtime();
		// sprzwdzenie poprawności rozkładu
		// A = L*U
		Matrix mlu = lu.getL() * lu.getU();
		cout << "LU MSE: " << mse(m, mlu) << endl;
		cout << "LU deltaTime: " << endTime - startTime << endl;

		// wczytanie wektora b
		std::ifstream bVec(args.bFileName);
		startTime = MPI_Wtime();
		// rozwiązanie układu równań Ax = b
		SolutionCalculator solutionCalculator(lu, &bVec, myid, numProcs, args.matrixSize);
		solutionCalculator.run();
		endTime = MPI_Wtime();
		// solutionCalculator.printSolutionVector();
		cout << "Solution deltaTime: " << endTime - startTime << endl;

		if (args.compareX)
		{
			// porównanie wyniku z wartością oczekiwaną
			std::ifstream targetVecFile(args.xFileName);
			cout << "Soulution MSE: " << solutionCalculator.mse(&targetVecFile) << endl;
		}
	}
	else
	{
		// czynności wezłów worker
		// rozłożenie macierzy
		LUDecomposition lu(nullptr, args.matrixSize, myid, numProcs);
		// rozwiązanie układu równań
		SolutionCalculator solutionCalculator(lu, nullptr, myid, numProcs, args.matrixSize);
		solutionCalculator.run();
	}
	MPI_Finalize();
	return 0;
}