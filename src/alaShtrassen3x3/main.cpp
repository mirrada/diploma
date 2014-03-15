#include "stdafx.h"

#ifndef _DEBUG
#include "mpi.h"
#else
#include "omp.h"
#endif

#include <stdio.h>
#include <string.h>
#include <complex>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <map>
#include <array>
#include <iomanip>
#include "mkl.h"

#define BruteForceKLRSstart for (int k = 0; k < N; k++)\
for (int l = 0; l < N; l++)\
for (int r = 0; r < N; r++){\
	int ii = (2 * i) % N; \
	int jj = (2 * j) % N; \
	int kk = (2 * k) % N; \
	int ll = (2 * l) % N; \
	int rr = (2 * r) % N; \
	int s = (i + k + r + 2 * l + 2 * j) % N; \
	int ss = (2 * s) % N;
#define BruteForceKLRSend   }

double WANTED = 1;
int MAXINSOLUTION = 1000000;
double MAXDELTA = 0.00001;

double const EPSMACH = 2.2e-16;
int const N = 3;
int const T = 3;
const int CONJ = 1;

int const MAXITER = 1000000;
int const MAXMATRIXINT = 100;

typedef std::complex<double> Complex;
typedef Complex ComplexArray[N][N];

ComplexArray* a[T];
ComplexArray* b[T];
ComplexArray* c[T];

Complex bcmul[T][N][N][N][N]; /* t, k, l, r, s pre-calculation b_tkl*c_trs */

std::ofstream logger;
char fileName[50];


double fdelta[N][N][N][N][N][N] = { 0 };

void setStaticFdelta() {
	for (int j = 0; j < N; j++)
	for (int l = 0; l < N; l++)
	for (int s = 0; s < N; s++) {
		fdelta[s][j][j][l][l][s] += 1.0 / 3;
		fdelta[j][j][l][l][s][s] -= 1.0 / 3;
	}
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

void prepareBCmul() {
	for (int t = 0; t < T; t++)
	for (int k = 0; k < N; k++)
	for (int l = 0; l < N; l++)
	for (int r = 0; r < N; r++)
	for (int s = 0; s < N; s++)
		bcmul[t][k][l][r][s] = (*b[t])[k][l] * (*c[t])[r][s];
}

double getResidual() {
	double resid = 0;
	for (int i = 0; i < N; i++)
	for (int j = 0; j < N; j++){
		BruteForceKLRSstart{
		Complex sum = 0;
		for (int t = 0; t < T; t++){
			sum += (*a[t])[i][j] * bcmul[t][k][l][r][s];
#ifdef GROUP1
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

/*------------------------------for LE solution------------------------------*/

int const SIZEB = 2; /*columns in right hand matrix b. first for system a_t00, second for a_t11+a_t22*/
int const LSNUM = 3; /*number of systems for i!=j*/
int const ROWSBIGLS = 6; /*size of systems for i!=j*/
Complex LS0011[T][T]; /*LESystem for i=0,j=0 and at00+at11*/
Complex b0011[T][SIZEB];  /*right hand vecors  for i=0,j=0 (the first column) and a00+a11 (the second column)*/

Complex LS22[T][T]; /*LESystem for at00-at11*/
Complex b22[T]; /*right hand vecors at00-at11*/

Complex LSIJ[LSNUM][ROWSBIGLS][ROWSBIGLS];  /*LESystem for i!=j*/
Complex bIJ[LSNUM][ROWSBIGLS];   /*right hand vecor i!=j*/

struct IndexesKLRS {
	int k, l, r, s;
	int kk, ll, rr, ss; // 0 -> 0, 1->2, 2->3:  k->kk
	IndexesKLRS(int a = 0, int b = 0, int c = 0, int d = 0) : k(a), l(b), r(c), s(d) {
		kk = (N - k) % N;
		ll = (N - l) % N;
		rr = (N - r) % N;
		ss = (N - s) % N;
	}
};

struct IndexesIJT {
	int i, j, t;
	IndexesIJT(int a = 0, int b = 0, int c = 0) : i(a), j(b), t(c) {}
};

class NumerationLSVariables {
	int toLS[N][N][N];
	IndexesIJT toA[ROWSBIGLS];
public:
	NumerationLSVariables(bool(*condition)(int, int)) {
		int m = 0;
		for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++) {
			if (condition(i, j)) {
				for (int t = 0; t < T; t++) {
					toLS[i][j][t] = m * T + t;
					toA[toLS[i][j][t]] = IndexesIJT(i, j, t);
				}
				m++;
			}
		}
	}
	int getOrder(int i, int j, int t) {
		return toLS[i][j][t];
	}
	IndexesIJT getIndexes(int m) {
		return toA[m];
	}
};

IndexesKLRS  bruteForce[N][N*N*N];

void setStaticBruteForce() {
	for (int j_i = 0; j_i < N; j_i++) {// j-i = k-l+r-s, list of suitable klrs
		int m = 0;
		for (int k = 0; k < N; k++)
		for (int l = 0; l < N; l++)
		for (int r = 0; r < N; r++) {
			int s = (k + 2 * l + r + 2 * j_i) % N;
			bruteForce[j_i][m] = IndexesKLRS(k, l, r, s);
			m++;
		}
	}
}

/*make LESystem for a_tii */
void fillLSII() {
	for (int t0 = 0; t0 < T; t0++) {
		for (int t = 0; t < T; t++) {
			Complex sum = 0;
			Complex sumdiff = 0;
			for (IndexesKLRS d : bruteForce[0]) {
				sum += bcmul[t][d.k][d.l][d.r][d.s] *
					std::conj(bcmul[t0][d.k][d.l][d.r][d.s] + bcmul[t0][d.kk][d.ll][d.rr][d.ss]);
				sumdiff += bcmul[t][d.k][d.l][d.r][d.s] *
					std::conj(bcmul[t0][d.k][d.l][d.r][d.s] - bcmul[t0][d.kk][d.ll][d.rr][d.ss]);
			}
			LS0011[t0][t] = sum;
			LS22[t0][t] = sumdiff;
		}
		b0011[t0][0] = b0011[t0][1] = b22[t0] = 0;
		for (IndexesKLRS d : bruteForce[0]) {
			b0011[t0][0] += fdelta[0][0][d.k][d.l][d.r][d.s] * std::conj(bcmul[t0][d.k][d.l][d.r][d.s]);
			b0011[t0][1] += fdelta[1][1][d.k][d.l][d.r][d.s] *
				std::conj(bcmul[t0][d.k][d.l][d.r][d.s] + bcmul[t0][d.kk][d.ll][d.rr][d.ss]);
			b22[t0] += fdelta[1][1][d.k][d.l][d.r][d.s] *
				std::conj(bcmul[t0][d.k][d.l][d.r][d.s] - bcmul[t0][d.kk][d.ll][d.rr][d.ss]);
		}
	}
}

/*make LESystem for a_tij : i!=j */
void fillLSIJ(int i, int j, NumerationLSVariables num) {
	int ii = 2 * i % N;
	int jj = 2 * j % N;
	for (int t0 = 0; t0 < T; t0++) {
		int row = num.getOrder(i, j, t0);
		int row2 = num.getOrder(ii, jj, t0);
		for (int t = 0; t < T; t++) {
			Complex sum = 0;
			Complex sumsum = 0;
			Complex sum2 = 0;
			Complex sumsum2 = 0;
			for (IndexesKLRS d : bruteForce[(j + 2 * i) % N]) {
				sum += bcmul[t][d.k][d.l][d.r][d.s] * std::conj(bcmul[t0][d.k][d.l][d.r][d.s]);
				sumsum += bcmul[t][d.kk][d.ll][d.rr][d.ss] * std::conj(bcmul[t0][d.k][d.l][d.r][d.s]);
				sum2 += bcmul[t][d.k][d.l][d.r][d.s] * std::conj(bcmul[t0][d.kk][d.ll][d.rr][d.ss]);
				sumsum2 += bcmul[t][d.kk][d.ll][d.rr][d.ss] * std::conj(bcmul[t0][d.kk][d.ll][d.rr][d.ss]);
			}
			LSIJ[i][row][num.getOrder(i, j, t)] = sum;
			LSIJ[i][row][num.getOrder(ii, jj, t)] = sumsum;
			LSIJ[i][row2][num.getOrder(i, j, t)] = sum2;
			LSIJ[i][row2][num.getOrder(ii, jj, t)] = sumsum2;
		}
		bIJ[i][row] = 0;
		bIJ[i][row2] = 0;
		for (IndexesKLRS d : bruteForce[(j + 2 * i) % N]) {
			bIJ[i][row] += fdelta[i][j][d.k][d.l][d.r][d.s] * std::conj(bcmul[t0][d.k][d.l][d.r][d.s]);
			bIJ[i][row2] += fdelta[i][j][d.k][d.l][d.r][d.s] * std::conj(bcmul[t0][d.kk][d.ll][d.rr][d.ss]);
		}
	}
}

void getAijFromSolution(NumerationLSVariables num[LSNUM]) {
	for (int t = 0; t < T; t++) {
		(*a[t])[0][0] = b0011[t][0];
		(*a[t])[1][1] = (b0011[t][1] + b22[t]) / 2.0;
		(*a[t])[2][2] = (b0011[t][1] - b22[t]) / 2.0;
	}
	for (int lsnum = 0; lsnum < LSNUM; lsnum++)
	for (int m = 0; m < ROWSBIGLS; m++) {
		auto d = num[lsnum].getIndexes(m);
		(*a[d.t])[d.i][d.j] = bIJ[lsnum][m];
	}
}

int prepareLEsol(NumerationLSVariables numIJ[3]) {
	lapack_int info[LSNUM + 2];
	lapack_int ipiv[ROWSBIGLS];

	fillLSII();
	fillLSIJ(0, 1, numIJ[0]);
	fillLSIJ(1, 0, numIJ[1]);
	fillLSIJ(2, 1, numIJ[2]);

	info[0] = LAPACKE_zgesv(LAPACK_ROW_MAJOR, T, SIZEB, (lapack_complex_double*)LS0011, T, ipiv, (lapack_complex_double*)b0011, SIZEB);
	info[1] = LAPACKE_zgesv(LAPACK_ROW_MAJOR, T, 1, (lapack_complex_double*)LS22, T, ipiv, (lapack_complex_double*)b22, 1);
	for (int i = 0; i < LSNUM; i++)
		info[i + 2] = LAPACKE_zgesv(LAPACK_ROW_MAJOR, ROWSBIGLS, 1, (lapack_complex_double*)(LSIJ[i]), ROWSBIGLS, ipiv, (lapack_complex_double*)(bIJ[i]), 1);

	for (int i = 0; i < LSNUM + 2; i++)
	if (info[i])
		return -1;

	getAijFromSolution(numIJ);

	return 0;
}

/*---------------------------------for LLS solution------------------------------*/
int const NROWpart = 27;
int const NROW = 5 * NROWpart;

std::array<std::array<Complex, N*N*T>, NROW> LLST;
std::array<Complex, NROW> bLLS;

std::array<std::array<Complex, 6>, NROWpart> LLST6;
std::array<std::array<Complex, 3>, NROWpart> LLST3;
std::array<Complex, NROWpart> bLLSpart;

/*int col = i*N*N + j*N;  int colcol = ii*N*N + jj*N;*/
template<unsigned int ROW, unsigned int COL>
void fillLLSsub(int i, int j, int rowStart, int col, int colcol, std::array<std::array<Complex, COL>, ROW> &arr, std::array<Complex, ROW> &b) {
	BruteForceKLRSstart{
	int row = rowStart + k*N*N + l*N + r;
	for (int t = 0; t < T; t++) {
		double one = 1;
		arr[row][col + t] += bcmul[t][k][l][r][s];
#ifdef GROUP1
		if (!i ^ !j ^ !k ^ !l ^ !r ^ !s)
			one = -1;
#endif
		arr[row][colcol + t] += one * bcmul[t][kk][ll][rr][ss];
		b[row] = fdelta[i][j][k][l][r][s];
	}
}BruteForceKLRSend
}

template<unsigned int ROW, unsigned int COL>
int prepareLLSSolOnePart(int i, int j, int colcol, std::array<std::array<Complex, COL>, ROW> &arr){
	memset(arr.data(), 0, sizeof(arr));

	fillLLSsub(i, j, 0, 0, colcol, arr, bLLSpart);
	int info = LAPACKE_zgels(LAPACK_ROW_MAJOR, 'N', ROW, COL, 1,
		(lapack_complex_double*)arr.data(), COL,
		(lapack_complex_double*)bLLSpart.data(), 1);

	if (info)
		return info;
	for (unsigned int k = 0; k < COL; k++)
	if (abs(bLLSpart[k]) > MAXINSOLUTION)
		bLLSpart[k] = 0;
	for (int t = 0; t < T; t++) {
		(*a[t])[i][j] = bLLSpart[t];
		(*a[t])[2 * i%N][2 * j%N] = bLLSpart[t + colcol];
	}
	return 0;
}

int prepareLLSSolPart(){
	if (std::rand() % 2)
	if (prepareLLSSolOnePart(0, 0, 0, LLST3))
		return -1;
	if (std::rand() % 2)
	if (prepareLLSSolOnePart(0, 1, T, LLST6))
		return -1;
	if (std::rand() % 2)
	if (prepareLLSSolOnePart(1, 0, T, LLST6))
		return -1;
	if (std::rand() % 2)
	if (prepareLLSSolOnePart(1, 1, T, LLST6))
		return -1;
	if (std::rand() % 2)
	if (prepareLLSSolOnePart(2, 1, T, LLST6))
		return -1;
}

void fillLLS(){
	memset(LLST.data(), 0, sizeof(LLST));
	fillLLSsub(0, 0, 0, 0, 0, LLST, bLLS);
	fillLLSsub(0, 1, 1 * N*N*N, 3, 6, LLST, bLLS);
	fillLLSsub(1, 0, 2 * N*N*N, 9, 18, LLST, bLLS);
	fillLLSsub(1, 1, 3 * N*N*N, 12, 24, LLST, bLLS);
	fillLLSsub(2, 1, 4 * N*N*N, 21, 15, LLST, bLLS);
}

void getAijFromLLS(){
	for (int i = 0; i < N; i++)
	for (int j = 0; j < N; j++)
	for (int t = 0; t < T; t++)
		(*a[t])[i][j] = bLLS[i*N*N + j*N + t];
}

int prepareLLSSol(){
	fillLLS();

	int info = LAPACKE_zgels(LAPACK_ROW_MAJOR, 'N', NROW, N*N*T, 1, (lapack_complex_double*)LLST.data(), N*N*T,
		(lapack_complex_double*)bLLS.data(), 1);

	if (info)
		return info;

	for (int i = 0; i < N*N*T; i++)
	if (abs(bLLS[i]) > MAXINSOLUTION) 
		bLLS[i] = 0;
	
	getAijFromLLS();
	return 0;
}

/*--------------------------------------------------------------------------------*/

/*returns 0 if OK. 1,2 if there is a problem with solution and makes resudual = WANTED + 1.*/
int nextApprox(double& delta, double& residual, NumerationLSVariables numIJ[3]){
	prepareBCmul();
	double newResid;

#ifdef GROUP1
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
	if (prepareLEsol(numIJ)) {
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
		logger.open(fileName, std::ios::app);
		logger << "resid  = " << residual << std::endl;
		printMatrixesABC();
		logger.close();
	}
}

void findSolutions(NumerationLSVariables numIJ[3]) {
	double residual, delta;
	long attempt = 0;

	while (true) {
		setRandomInitialBC();

		delta = MAXDELTA + 1;
		residual = MAXINSOLUTION + 1;

		long cycle = 0;
		while (abs(delta) > MAXDELTA) {
#ifdef _DEBUG 
			double t = omp_get_wtime();
#endif
			if (nextApprox(delta, residual, numIJ))
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

void morePrecise(NumerationLSVariables numIJ[3]) {
	double residual = MAXINSOLUTION, delta = 1;

	readMatrixesABC("input.txt");

	for (int trial = 0; trial < MAXITER; trial++) {
		std::cout << trial << std::endl;
		if (nextApprox(delta, residual, numIJ)){
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
		nextApprox(delta, residual, numIJ);
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
	setStaticBruteForce();

	auto condIJ01 = [](int i, int j){return (i == 0 && j != 0); };
	auto condIJ10 = [](int i, int j){return (i != 0 && j == 0); };
	auto condIJ21 = [](int i, int j){return (i != j && i * j != 0); };
	NumerationLSVariables numIJ[3] = { NumerationLSVariables(condIJ01),
		NumerationLSVariables(condIJ10),
		NumerationLSVariables(condIJ21) };

	for (int i = 0; i < 3; i++) {
		a[i] = new ComplexArray[1];
		b[i] = new ComplexArray[1];
		c[i] = new ComplexArray[1];
	}

#ifdef PRECISE
	morePrecise(numIJ);
#else
	findSolutions(numIJ);
#endif


#ifndef _DEBUG
	MPI_Finalize();
#endif
	return(0);
}


