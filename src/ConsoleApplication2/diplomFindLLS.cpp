#include "stdafx.h"
#include <stdio.h>
#include <string.h>
#include <complex>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <map>
#include <iomanip>
#include "lapacke.h"
#include "mpi.h"


double WANTED = 2;
int MAXINSOLUTION = 1000000;
double MAXDELTA = 0.00001;
bool COMPLEX = 1;
int MAXMATRIXINT = 3;



double const EPSMACH = 2.2e-16;
int const MAXITER = 1000;
int const N = 3;
int const T = 3;
bool const random = 1;

typedef std::complex<double> Complex;
typedef Complex ComplexArray[N][N];
ComplexArray* a[T];
ComplexArray* b[T];
ComplexArray* c[T];

Complex bcmul[T][N][N][N][N]; /* t, k, l, r, s pre-calculation b_tkl*c_trs */

std::ofstream logger;
char fileName[40];


double fdelta[N][N][N][N][N][N] = { 0 };

void setStaticFdelta() {
//#pragma omp parallel for
	for (int j = 0; j < N; j++)
	for (int l = 0; l < N; l++)
	for (int s = 0; s < N; s++) {
		fdelta[s][j][j][l][l][s] += 1.0 / 3;
		fdelta[j][j][l][l][s][s] -= 1.0 / 3;
	}
}

void setRandomInitialBC() {
	int n = MAXMATRIXINT;
	for (int i = 0; i < 3; i++)
	for (int j = 0; j < 3; j++)
	for (int k = 0; k < 3; k++) {
		(*b[i])[j][k] = Complex(std::rand() % (2 * n) - n,
			COMPLEX * (std::rand() % (2 * n) - n));
		(*c[i])[j][k] = Complex(std::rand() % (2 * n) - n,
			COMPLEX * (std::rand() % (2 * n) - n));
	}
}

void prepareBCmul() {
//#pragma omp parallel for 
	for (int t = 0; t < T; t++)
	for (int k = 0; k < N; k++)
	for (int l = 0; l < N; l++)
	for (int r = 0; r < N; r++)
	for (int s = 0; s < N; s++)
		bcmul[t][k][l][r][s] = (*b[t])[k][l] * (*c[t])[r][s];
}

double getResidual() {
	double resid = 0;
//#pragma omp parallel for reduction(+ : resid) 
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++){
			for (int k = 0; k < N; k++){
				for (int l = 0; l < N; l++){
					for (int r = 0; r < N; r++){
						int ii = (N - i) % N;
						int jj = (N - j) % N;
						int kk = (N - k) % N;
						int ll = (N - l) % N;
						int rr = (N - r) % N;
						int s = (i + k + r + 2 * l + 2 * j) % N;
						int ss = (N - s) % N;
						Complex sum = 0;
						for (int t = 0; t < T; t++){
							sum += (*a[t])[i][j] * bcmul[t][k][l][r][s];
							sum += (*a[t])[ii][jj] * bcmul[t][kk][ll][rr][ss];
						}
						sum -= fdelta[i][j][k][l][r][s];
						resid += std::norm(sum);
					}
				}
			}
		}
	}
	return resid / 2;
}

void swapABC() {
	Complex(*temp)[N][N];
	for (int t = 0; t < T; t++) {
		temp = a[t];
		a[t] = b[t];
		b[t] = c[t];
		c[t] = temp;
	}
}

void print_matrix(char* desc, int z, int n, int m, Complex* a) {
	lapack_int i, j;
	logger << desc << z << std::endl;
	logger << "{ ";
	for (i = 0; i < n; i++) {
		for (j = 0; j < m; j++)
			logger << "Complex" << a[i*m + j] << ", ";
		logger << std::endl << "  ";
	}
	logger << "}" << std::endl;
}

void printMatrixesABC(){
	for (int t = 0; t < T; t++) {
		print_matrix("a", t, N, N, (Complex*)a[t]);
		print_matrix("b", t, N, N, (Complex*)b[t]);
		print_matrix("c", t, N, N, (Complex*)c[t]);
	}
}

/*------------------------------for LE solution---------------------------------------------------------------------------------------*/


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

	return 0;
}

/*---------------------------------for LLS solution------------------------------------------------------------------------------------*/
int const NROW = 135;
Complex  LLST[N*N*T][NROW] = { 0 };
Complex bLLS[NROW] = { 0 };

void fillLLSsub(int i, int j, int rowStart) {
	int jj = 2 * j % N;
	int ii = 2 * i % N;
	int col = i*N*N + j*N;
	int colcol = ii*N*N + jj*N;
	for (int k = 0; k < N; k++)
	for (int l = 0; l < N; l++)
	for (int r = 0; r < N; r++) {
		int s = (i + k + r + 2 * (l + j)) % N;
		int kk = (2 * k) % N;
		int ll = (2 * l) % N;
		int rr = (2 * r) % N;
		int ss = (2 * s) % N;
		int row = rowStart + k*N*N + l*N + r;
		for (int t = 0; t < T; t++) {
			LLST[col + t][row] = bcmul[t][k][l][r][s];
			LLST[colcol + t][row] = bcmul[t][kk][ll][rr][ss];
			bLLS[row] = fdelta[i][j][k][l][r][s];
		}
	}
}

void fillLLS(){

	memset(LLST, 0, sizeof(LLST));
	//int i = 0; int j = 0;
//#pragma omp parallel for 
	for (int k = 0; k < N; k++)
	for (int l = 0; l < N; l++)
	for (int r = 0; r < N; r++) {
		int s = (k + r + 2 * l) % N;
		int kk = (2 * k) % N;
		int ll = (2 * l) % N;
		int rr = (2 * r) % N;
		int ss = (2 * s) % N;
		int row = k*N*N + l*N + r;
		for (int t = 0; t < T; t++) {
			LLST[t][row] = bcmul[t][k][l][r][s] + bcmul[t][kk][ll][rr][ss];
			bLLS[row] = fdelta[0][0][k][l][r][s];
		}
	}

//#pragma omp parallel sections
	{
//#pragma omp section
		{
			fillLLSsub(0, 1, N*N*N);
		}
//#pragma omp section
		{
			fillLLSsub(1, 0, 2 * N*N*N);
		}
//#pragma omp section
		{
			fillLLSsub(1, 1, 3 * N*N*N);
		}
//#pragma omp section
		{
			fillLLSsub(2, 1, 4 * N*N*N);
		}
	}
}

void getAijFromLLS(){
//#pragma omp parallel for 
	for (int i = 0; i < N; i++)
	for (int j = 0; j < N; j++)
	for (int t = 0; t < T; t++)
		(*a[t])[i][j] = bLLS[i*N*N + j*N + t];
}

double prepareLLSSol(){

	fillLLS();
	//transposeLLS();

	int n4 = NROW;
	int n2 = N*N*T;
	int nrhs = 1;
	lapack_complex_double* work;
	lapack_complex_double wkopt;
	int lwork = -1;
	lapack_int info;

	LAPACK_zgels("No transpose", &n4, &n2, &nrhs, (lapack_complex_double*)LLST, &n4,
		(lapack_complex_double*)bLLS, &n4, &wkopt, &lwork, &info);
	lwork = (int)wkopt.real;
	work = new lapack_complex_double[lwork];
	/* Solve the equations LLS*X = bLLS */
	LAPACK_zgels("No transpose", &n4, &n2, &nrhs, (lapack_complex_double*)LLST, &n4,
		(lapack_complex_double*)bLLS, &n4, work, &lwork, &info);
	delete[] work;

	if (info)
		return -1;

	for (int i = 0; i < n2; i++)
	if (abs(bLLS[i]) > MAXINSOLUTION) 
		return -1;

	double norm = 0.0;
//#pragma omp parallel for reduction(+ : norm)
	for (int i = n2; i < n4; i++)
		norm += std::norm(bLLS[i]);
	return norm;
}

/*---------------------------------------------------------------------------------------------------------------------------------------*/


/*returns 0 if OK. 1,2 if there is a problem with solution and makes resudual = WANTED + 1.
*/
int nextApprox(double& delta, double& residual, NumerationLSVariables numIJ[3]){
	prepareBCmul();
	double newResid;

	/*prepareLLSSol();
	getAijFromLLS();

	newResid = getResidual();
	if (newResid < 0) {
		residual = WANTED + 1;
		return 1;
	}
	delta = residual - newResid;
	if (delta < -EPSMACH) {
		residual = WANTED + 1;
		return 2;
	}
	residual = newResid;*/

	prepareLEsol(numIJ);
	getAijFromSolution(numIJ);

	newResid = getResidual();
	if (newResid < 0) {
		residual = WANTED + 1;
		return 1;
	}
	delta = residual - newResid;
	if (delta < -EPSMACH) {
		residual = WANTED + 1;
		return 2;
	}
	residual = newResid;
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

void byRandom(NumerationLSVariables numIJ[3]) {
	double residual, delta;
	long attempt = 0;

	while (true) {
		setRandomInitialBC();

		delta = MAXDELTA + 1;
		residual = MAXINSOLUTION + 1;

		long cycle = 0;
		while (abs(delta) > MAXDELTA) {
			double t = omp_get_wtime();
			if (nextApprox(delta, residual, numIJ))
				break;
			cycle++;
			if (cycle % 5000 == 0)
				std::cout << "cycle = " << cycle << std::endl;

			std::cout << omp_get_wtime() - t << std::endl;
			t = omp_get_wtime();
		}

		printApprox(residual);
		attempt++;
		if (attempt % 100 == 0){
			std::cout << "attempt = " << attempt << std::endl;
			std::cout << "resid = " << residual << std::endl;
		}
	}
}

void hardCode(NumerationLSVariables numIJ[3]) {
	double residual, delta;
	for (int trial = 0; trial < MAXITER; trial++) {
		delta = 1;
		residual = 10;
		prepareBCmul();
		nextApprox(delta, residual, numIJ);
		printApprox(residual);
	}
}

int main(int argc, char* argv[]) {
	int k = 0;

	MPI_Init(&argc, &argv);                       /* starts MPI */
	MPI_Comm_rank(MPI_COMM_WORLD, &k); /* get current process id */

	for (int i = 0; i < 3; i++) {
		b[i] = new ComplexArray[1];
		c[i] = new ComplexArray[1];
		a[i] = new ComplexArray[1];
	}

	srand(time(NULL) + k*k*1000);

	sprintf_s(fileName, "logger%f_%d_%d_%dproc.txt", WANTED, MAXMATRIXINT, COMPLEX, k);

	logger.open(fileName);
	logger << std::fixed;
	logger.precision(5);
	logger.close();

	setStaticFdelta();
	setStaticBruteForce();

	auto condIJ01 = [](int i, int j){return (i == 0 && j != 0); };
	auto condIJ10 = [](int i, int j){return (i != 0 && j == 0); };
	auto condIJ21 = [](int i, int j){return (i != j && i * j != 0); };
	NumerationLSVariables numIJ[3] = { NumerationLSVariables(condIJ01),
		NumerationLSVariables(condIJ10),
		NumerationLSVariables(condIJ21) };

	if (random)
		logger.precision(5);
	else
		logger.precision(15);
	logger.close();


	if (random) {
		byRandom(numIJ);
	}
	else {
		hardCode(numIJ);
	}

	for (int i = 0; i < 3; i++) {
		delete[] a[i];
		delete[] b[i];
		delete[] c[i];
	}

	MPI_Finalize();
	return(0);
}


