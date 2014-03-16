#include "stdafx.h"

#ifndef _DEBUG
#include <mpi.h>
#else
#include <omp.h>
#endif

#ifndef GROUP2
#include "LEG1sol.h"
#else
#include "LLSsol.h"
#endif

#include "common.h"
#include <iostream>
#include <fstream>
#include <iomanip>


double WANTED = 1;
double MAXDELTA = 1e-5;
double MAXRESID = 1e+6;

double const EPSMACH = 2.2e-16;
const int CONJ = 0;

int const MAXITER = 10000;
int const MAXMATRIXINT = 10;

std::ofstream logger;
char fileName[50];


void setStaticFdelta() {
	for (int j = 0; j < N; j++)
		for (int l = 0; l < N; l++)
			for (int s = 0; s < N; s++) {
				fdelta[s][j][j][l][l][s] += 1.0 / 3;
				fdelta[j][j][l][l][s][s] -= 1.0 / 3;
			}
}

void prepareBCmul() {
	for (int t = 0; t < T; t++)
		for (int k = 0; k < N; k++)
			for (int l = 0; l < N; l++)
				for (int r = 0; r < N; r++)
					for (int s = 0; s < N; s++)
						bcmul[t][k][l][r][s] = (*b[t])[k][l] * (*c[t])[r][s];
}

void setRandomInitialBC() {
	int n = MAXMATRIXINT;
	if (CONJ) {
		for (int k = 0; k < N; k++) {
			for (int i = 0; i < N - 1; i++)
			for (int j = 0; j < i + 2; j++){
				(*a[k])[i][j] = Complex(std::rand() % (2 * n) - n,
					(std::rand() % (2 * n) - n));
				(*a[k])[2 * i%N][2 * j%N] = std::conj((*a[k])[i][j]);
				(*b[k])[i][j] = Complex(std::rand() % (2 * n) - n,
					(std::rand() % (2 * n) - n));
				(*b[k])[2 * i%N][2 * j%N] = std::conj((*b[k])[i][j]);
				(*c[k])[i][j] = Complex(std::rand() % (2 * n) - n,
					(std::rand() % (2 * n) - n));
				(*c[k])[2 * i%N][2 * j%N] = std::conj((*c[k])[i][j]);
			}
			(*a[k])[0][0] = std::rand() % (2 * n) - n;
			(*b[k])[0][0] = std::rand() % (2 * n) - n;
			(*c[k])[0][0] = std::rand() % (2 * n) - n;
		}
	}
	else {
		int n = MAXMATRIXINT;
		for (int k = 0; k < N; k++)
		for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++){
			(*a[k])[i][j] = Complex(std::rand() % (2 * n) - n,
				(std::rand() % (2 * n) - n));
			(*b[k])[i][j] = Complex(std::rand() % (2 * n) - n,
				(std::rand() % (2 * n) - n));
			(*c[k])[i][j] = Complex(std::rand() % (2 * n) - n,
				(std::rand() % (2 * n) - n));
		}
	}
}

double getResidual() {
	double resid = 0;
	for (int i = 0; i < N; i++)
	for (int j = 0; j < N; j++){
		BruteForceKLRSstart{
		Complex sum = 0;
		for (int t = 0; t < T; t++){
			sum += (*a[t])[i][j] * bcmul[t][k][l][r][s];
#ifdef GROUP2
			if (!i ^ !j ^ !k ^ !l ^ !r ^ !s)
				sum -= (*a[t])[ii][jj] * bcmul[t][kk][ll][rr][ss];
			else
#endif
				sum += (*a[t])[ii][jj] * bcmul[t][kk][ll][rr][ss];
		}
		sum -= fdelta[i][j][k][l][r][s];
		double temp = std::norm(sum);
		resid += temp;
	}BruteForceKLRSend
	}
	return resid / 2;
}

void swapABC() {
	ComplexArray* temp;
	for (int t = 0; t < T; t++) {
		temp = a[t];
		a[t] = b[t];
		b[t] = c[t];
		c[t] = temp;
	}
}

double findMaxInMatrices(ComplexArray* m[T]){
	double max = 0;
	for (int t = 0; t < T; t++)
	for (int i = 0; i < N; i++)
	for (int j = 0; j < N; j++) {
		if (abs((*m[t])[i][j].imag()) > max)
			max = abs((*m[t])[i][j].imag());
		if (abs((*m[t])[i][j].real()) > max)
			max = abs((*m[t])[i][j].real());
	}
	return max;
}

void mulMatrix(ComplexArray* m[T], double alfa){
	for (int t = 0; t < T; t++)
	for (int i = 0; i < N; i++)
	for (int j = 0; j < N; j++) 
		(*m[t])[i][j] *= alfa;
}

void normalizeABC(){
	double aMax = findMaxInMatrices(a);
	double bMax = findMaxInMatrices(b);
	double cMax = findMaxInMatrices(c);
	mulMatrix(a, std::pow(bMax*cMax / aMax / aMax, (double)1 / 3));
	mulMatrix(b, std::pow(aMax*cMax / bMax / bMax, (double)1 / 3));
	mulMatrix(c, std::pow(bMax*aMax / cMax / cMax, (double)1 / 3));
}

void printMatrix(char* desc, int z, int n, int m, Complex* a) {
	logger << desc << z << std::endl;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++)
			logger << a[i*m + j] << "  ";
		logger << std::endl;
	}
}

void printMatrixesABC(){
	for (int t = 0; t < T; t++) {
		printMatrix("a", t, N, N, (Complex*)a[t]);
		printMatrix("b", t, N, N, (Complex*)b[t]);
		printMatrix("c", t, N, N, (Complex*)c[t]);
	}
}

void readMatrix(std::ifstream& input, int n, int m, Complex* a) {
	for (int i = 0; i < n; i++)
	for (int j = 0; j < m; j++)
		input >> a[i*m + j];
}

void readMatrixesABC(char* fileName){
	std::ifstream inputFile;
	inputFile.open(fileName, std::ifstream::in);
	for (int t = 0; t < T; t++) {
		readMatrix(inputFile, N, N, (Complex*)a[t]);
		readMatrix(inputFile, N, N, (Complex*)b[t]);
		readMatrix(inputFile, N, N, (Complex*)c[t]);
	}
}

/*--------------------------------------------------------------------------------*/

/*returns 0 if OK. 1,2 if there is a problem with solution and makes resudual = WANTED + 1.*/
int nextApprox(double& delta, double& residual){
	prepareBCmul();
	double newResid;

#ifdef GROUP2
	if (prepareLLSSolPart()) {
		residual = WANTED + 1;
		return 1;
	}
	newResid = getResidual();
	if (newResid < 0) {
		residual = WANTED + 2;
		return 2;
	}
	delta = residual - newResid;
	if (delta < -MAXDELTA) {
		residual = WANTED + 3;
		return 3;
	}
	residual = newResid;
#else
	if (prepareLEsol()) {
		residual = WANTED + 1;
		return 1;
	}
	newResid = getResidual();
	if (newResid < 0) {
		residual = WANTED + 2;
		return 2;
	}
	delta = residual - newResid;
	if (delta < -EPSMACH) {
		residual = WANTED + 3;
		return 3;
	}
	residual = newResid;
#endif

#ifdef _DEBUG
	logger.open(fileName, std::ios::app);
	logger << "resid  = " << residual << std::endl;
	printMatrixesABC();
	logger.close();

	if (prepareLLSSolPart()) {
		residual = WANTED + 1;
		return 1;
	}
	newResid = getResidual();

	logger.open(fileName, std::ios::app);
	logger << "resid  = " << newResid << std::endl;
	printMatrixesABC();
	logger.close();
#endif
	swapABC();
	return 0;
}

void printApprox(double residual) {
	if (residual <= WANTED) {
		if (residual > EPSMACH)
			WANTED = residual + MAXDELTA / 100;
		normalizeABC();
		prepareBCmul();
		double newResidual = getResidual();
		logger.open(fileName, std::ios::app);
		logger << "resid  = " << residual << std::endl;
		printMatrixesABC();
		logger.close();
	}
}

void findSolutions() {
	double residual, delta;
	long attempt = 0;

	while (true) {
		setRandomInitialBC();

		delta = MAXDELTA + 1;
		residual = MAXRESID + 1;

		long cycle = 0;
		while (abs(delta) > MAXDELTA) {
#ifdef _DEBUG 
			double t = omp_get_wtime();
#endif
			if (nextApprox(delta, residual))
				break;
			cycle++;
			if (cycle % 5000 == 0)
				std::cout << "cycle = " << cycle << std::endl;

#ifdef _DEBUG
			std::cout << omp_get_wtime() - t << std::endl;
			t = omp_get_wtime();
#endif
		}

		printApprox(residual);
		attempt++;
		if (attempt % 100 == 0){
			std::cout << "attempt = " << attempt << std::endl;
			std::cout << "resid = " << residual << std::endl;
		}
	}
}

void morePrecise() {
	double residual = MAXRESID, delta = 1;

	readMatrixesABC("input.txt");
	normalizeABC();

	for (int trial = 0; trial < MAXITER; trial++) {
		std::cout << trial << std::endl;
		if (nextApprox(delta, residual)){
			std::cout << delta << std::endl;
			logger.open(fileName, std::ios::app);
			logger << "resid  = " << residual << " delta  = " << delta << std::endl;
			logger << "trial  = " << trial << std::endl;
			printMatrixesABC();
			logger.close();
			break;
		}
	}
	for (int trial = 0; trial < 10; trial++) {
		nextApprox(delta, residual);
		logger.open(fileName, std::ios::app);
		logger << "resid  = " << residual << " delta  = " << delta << std::endl;
		logger << "trial  = " << trial << std::endl;
		printMatrixesABC();
		logger.close();
	}
}

int main(int argc, char* argv[]) {
	int k = 0;
#ifndef _DEBUG
	MPI_Init(&argc, &argv);                       /* starts MPI */
	MPI_Comm_rank(MPI_COMM_WORLD, &k); /* get current process id */
#endif

	srand(static_cast<unsigned int>(time(NULL)) + k*k * 1000);

#ifdef PRECISE
	sprintf_s(fileName, "loggerPrecise_%dproc.txt", k);
#else
	sprintf_s(fileName, "loggerFind_maxint%d_%dproc.txt", MAXMATRIXINT, k);
#endif

	logger.open(fileName);
	logger << std::fixed;
	logger.precision(15);
	logger.close();

	setStaticFdelta();

	for (int i = 0; i < 3; i++) {
		a[i] = new ComplexArray[1];
		b[i] = new ComplexArray[1];
		c[i] = new ComplexArray[1];
	}

#ifndef GROUP2
	setStaticBruteForce();
	setStaticNumerationLSVariables();
#endif


#ifdef PRECISE
	morePrecise();
#else
	findSolutions();
#endif


#ifndef _DEBUG
	MPI_Finalize();
#endif
	return(0);
}


