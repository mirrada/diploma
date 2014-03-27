#include "stdafx.h"

#ifndef _DEBUG
#include <mpi.h>
#else
#include <omp.h>
#endif

#ifndef GROUP2
#include "LEG1sol.h"
#endif

#include "LLSsol.h"
#include "common.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>


ComplexArray* a[T];
ComplexArray* b[T];
ComplexArray* c[T];

Complex bcmul[T][N][N][N][N] = { 0 };		/* t, k, l, r, s pre-calculation b_tkl*c_trs */
double fdelta[N][N][N][N][N][N] = { 0 };
double penalty[T][N][N][N][N][N][N] = { 0 };

double WANTED[8] = { 1, 1, 1, 1, 1, 1, 1, 1 };

double const MAXDELTA = 1e-5;
double const MAXRESID = 1e+6;
double const PENALTYINF = 1e+20;
double const PENALTYNULL = 1e-20;

int const CONJ = 0;
int const MAXITER = 100000;
int const MAXMATRIXINT = 10;

std::ofstream logger;

FileNameString fileNames[8];



void setStaticFdelta() {
	for (int j = 0; j < N; j++)
		for (int l = 0; l < N; l++)
			for (int s = 0; s < N; s++) {
				fdelta[s][j][j][l][l][s] += 1.0 / 3;
				fdelta[j][j][l][l][s][s] -= 1.0 / 3;
			}
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

void swapPenalty(){
	double temp[T][N][N][N][N][N][N] = { 0 };
	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
			for (int k = 0; k < N; k++)
				for (int l = 0; l < N; l++)
					for (int r = 0; r < N; r++){
						int s = (i + k + r + 2 * l + 2 * j) % N;
						for (int t = 0; t < T; t++)
							temp[t][k][l][r][s][i][j] = penalty[t][i][j][k][l][r][s];
					}
	memcpy(penalty, temp, sizeof(temp));
}

double getResidual(bool withPenalty) {
	double resid = 0;

	for (int i = 0; i < 2; i++)
		for (int j = 0; j < std::min(N, i + 2); j++)
			for (int k = 0; k < std::min(N, i + j + 2); k++)
				for (int l = 0; l < std::min(N, i + j + k + 2); l++)
					for (int r = 0; r < std::min(N, i + j + k + l + 2); r++){
						int ii = (2 * i) % N;
						int jj = (2 * j) % N;
						int kk = (2 * k) % N;
						int ll = (2 * l) % N;
						int rr = (2 * r) % N;
						int s = (i + k + r + 2 * l + 2 * j) % N;
						int ss = (2 * s) % N;
						Complex sum = 0;
						for (int t = 0; t < T; t++){
							sum += (*a[t])[i][j] * bcmul[t][k][l][r][s];

#ifdef GROUP2
							if (!i ^ !j ^ !k ^ !l ^ !r ^ !s)
								sum -= (*a[t])[ii][jj] * bcmul[t][kk][ll][rr][ss];
							else
#endif
								sum += (*a[t])[ii][jj] * bcmul[t][kk][ll][rr][ss];
							if (withPenalty) {
								resid += std::norm((*a[t])[i][j] * bcmul[t][k][l][r][s] * penalty[t][i][j][k][l][r][s]);
								if (i + j + k + l + r + s != 0) {
									resid += std::norm((*a[t])[ii][jj] * bcmul[t][kk][ll][rr][ss] * penalty[t][ii][jj][kk][ll][rr][ss]);
								}
							}
						}
						sum -= fdelta[i][j][k][l][r][s];
						double temp = std::norm(sum);
						/*if (fdelta[i][j][k][l][r][s] > EPSMACH || fdelta[i][j][k][l][r][s] < -EPSMACH) {
							resid += temp;
							}*/
						resid += temp;
					}
	return resid;
}

void printForTrefilov() {
	logger.open(fileNames[0], std::ios::app);
	for (int i = 0; i < 2; i++)
		for (int j = 0; j < std::min(N, i + 2); j++)
			for (int k = 0; k < std::min(N, i + j + 2); k++)
				for (int l = 0; l < std::min(N, i + j + k + 2); l++)
					for (int r = 0; r < std::min(N, i + j + k + l + 2); r++){
						int ii = (2 * i) % N;
						int jj = (2 * j) % N;
						int kk = (2 * k) % N;
						int ll = (2 * l) % N;
						int rr = (2 * r) % N;
						int s = (i + k + r + 2 * l + 2 * j) % N;
						int ss = (2 * s) % N;
						logger << i << j << k << l << r << s << std::endl;
						for (int t = 0; t < T; t++){
							logger << " + " << (*a[t])[i][j] * bcmul[t][k][l][r][s];
#ifdef GROUP2
							if (!i ^ !j ^ !k ^ !l ^ !r ^ !s)
								logger << "-" << (*a[t])[ii][jj] * bcmul[t][kk][ll][rr][ss];
							else
#endif
								logger << "+" << (*a[t])[ii][jj] * bcmul[t][kk][ll][rr][ss];
						}
						logger << '=' << fdelta[i][j][k][l][r][s] << std::endl;
					}
	logger.close();
}

Complex getTrace(ComplexArray& a) {
	return a[0][0] + a[1][1] + a[2][2];
}

Complex getMulTraces(int t) {
	return getTrace(*a[t])*getTrace(*b[t])*getTrace(*c[t]);
}

/*--------------------------------------------------------------------------------*/

void printApprox(char fileName[50], double residual, double &wanted) {
	if (residual <= wanted) {
		if (residual > EPSMACH) {
			wanted = residual + MAXDELTA / 100;
			WANTED[2] = 1;
			WANTED[3] = 1;
			WANTED[6] = 0.28;
		}
		normalizeABC();
		logger.open(fileName, std::ios::app);
		logger << "resid  = " << residual << std::endl;
		logger << "traces = " << getMulTraces(0) << "; " << getMulTraces(1) << "; " << getMulTraces(2) << std::endl;
		printMatrixesABC();
		logger.close();
	}
}

double thirdMax(int t, int i, int j, int k, int l, int r, int s, double(&maxArr)[3]){
	double mul = std::norm((*a[t])[i][j] * bcmul[t][k][l][r][s]);
	if (mul > maxArr[0]){ maxArr[2] = maxArr[1]; maxArr[1] = maxArr[0]; maxArr[0] = mul; }
	else
		if (mul > maxArr[1]) { maxArr[2] = maxArr[1]; maxArr[1] = mul; }
		else
			if (mul > maxArr[2]) { maxArr[2] = mul; }
	return mul;
};

bool setPenalty(int t, int i, int j, int k, int l, int r, int s, double resid, double penaltyVal, double(&maxArr)[3]){
	double mul = std::norm((*a[t])[i][j] * bcmul[t][k][l][r][s]);
	double p = pow(maxArr[2] / mul, 2);
	if (mul < PENALTYNULL) {
		double na = std::norm((*a[t])[i][j]);
		double nb = std::norm((*b[t])[k][l]);
		double nc = std::norm((*c[t])[r][s]);
		if (na * 100 < nb && na * 100 < nc) {
			(*a[t])[i][j] = 0;
			penalty[t][i][j][k][l][r][s] = PENALTYINF*penaltyVal;
			return true;
		}
		if (nb * 100 < na && nb * 100 < nc) {
			(*b[t])[k][l] = 0;
			penalty[t][i][j][k][l][r][s] = PENALTYINF*penaltyVal;
			return true;
		}
		if (nc * 100 < nb && nc * 100 < na) {
			(*c[t])[r][s] = 0;
			penalty[t][i][j][k][l][r][s] = PENALTYINF*penaltyVal;
			return true;
		}
		penalty[t][i][j][k][l][r][s] = 10 * penaltyVal;
	}
	else {
		if (mul > resid)
			penalty[t][i][j][k][l][r][s] = 0;
		else
			penalty[t][i][j][k][l][r][s] = std::min(p * penaltyVal, penaltyVal * 10);
	}
	return false;
};

int setPenaltyFunction(double resid, double penaltyVal){
	int maxNotNull;
	bool not000000 = true;
	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
			for (int k = 0; k < N; k++)
				for (int l = 0; l < N; l++)
					for (int r = 0; r < N; r++){
						int ii = (2 * i) % N;
						int jj = (2 * j) % N;
						int kk = (2 * k) % N;
						int ll = (2 * l) % N;
						int rr = (2 * r) % N;
						int s = (i + k + r + 2 * l + 2 * j) % N;
						int ss = (2 * s) % N;
						double maxArr[3] = { 0 };
						for (int t = 0; t < T; t++){
							thirdMax(t, i, j, k, l, r, s, maxArr);
							thirdMax(t, ii, jj, kk, ll, rr, ss, maxArr);
						}
						int not0count = 6;
						for (int t = 0; t < T; t++){
							not0count -= setPenalty(t, i, j, k, l, r, s, resid, penaltyVal, maxArr);
							not0count -= setPenalty(t, ii, jj, kk, ll, rr, ss, resid, penaltyVal, maxArr);
						}
						not000000 = ((bool)not0count) && not000000;
						maxNotNull = std::max(not0count, maxNotNull);
					}
	if (not000000)
		return maxNotNull;
	else
		return 7;
}


/*returns 0 if OK. 1,2 if there is a problem with solution and makes resudual = WANTED + 1.*/
int nextApprox(double& delta, double& residual){
	prepareBCmul();

	double theirResid = prepareLLSSolPart();
	double newResid = getResidual(true);
	double diff = theirResid - newResid;

	if (newResid < 0) {
		residual = MAXRESID + 1;
		return 1;
	}
	if (newResid < 0) {
		residual = MAXRESID + 2;
		return 2;
	}
	delta = residual - newResid;
	if (delta < -EPSMACH) {
		residual = MAXRESID + 3;
		return 3;
	}
	residual = newResid;

#ifdef _DEBUG
	prepareBCmul();
	newResid = getResidual(false);
#endif
	swapABC();
	return 0;
}

void getSolutionWithNulls(double resid){
	double residSystem = 0.2;
	double penaltyVal = 0.0001;
	double delta = MAXDELTA + 1;
	double residual;
	int countNot0 = 6;
	while (residSystem < 0.4 && residSystem > 0.1 && penaltyVal < 1000000 && countNot0 > 2 && countNot0 < 7) {
		long cycle = 0;
		delta = MAXDELTA + 1;
		residual = MAXRESID;
		while (abs(delta) > MAXDELTA) {
			if (nextApprox(delta, residual))
				break;
			swapPenalty();
			cycle++;
			if (cycle % 500 == 0)
				std::cout << "cycle = " << cycle << " penalty = " << penaltyVal << std::endl;
		}
		prepareBCmul();
		residSystem = getResidual(false);
		countNot0 = setPenaltyFunction(residSystem, penaltyVal);
		printApprox(fileNames[countNot0], residSystem, WANTED[countNot0]);
		penaltyVal *= 1.01;
	}
}

void findSolutions() {
	double residual, delta;
	long attempt = 0;

	while (true) {
		setRandomInitialBC();
		memset(penalty, 0, sizeof(penalty));

		delta = MAXDELTA + 1;
		residual = MAXRESID + 1;

		long cycle = 0;
		while (abs(delta) > MAXDELTA) {
#ifdef _DEBUG 
			//double t = omp_get_wtime();
#endif
			if (nextApprox(delta, residual))
				break;
			cycle++;
			if (cycle % 5000 == 0)
				std::cout << "cycle = " << cycle << std::endl;

#ifdef _DEBUG
			//std::cout << omp_get_wtime() - t << std::endl;
			//t = omp_get_wtime();
#endif
		}
		prepareBCmul();
		double residSystem = getResidual(false);
		printApprox(fileNames[6], residSystem, WANTED[6]);
		if (residSystem < 1)
			getSolutionWithNulls(residual);
		attempt++;
		if (attempt % 10 == 0){
			std::cout << "attempt = " << attempt << std::endl;
			std::cout << "resid = " << residSystem << std::endl;
		}
	}
}

void morePrecise() {
	double residual = MAXRESID, delta = 1;

	readMatrixesABC("input.txt");
	normalizeABC();

	prepareBCmul();
	printForTrefilov();

	for (int trial = 0; trial < MAXITER; trial++) {
		std::cout << trial << std::endl;
		if (nextApprox(delta, residual)){
			std::cout << delta << std::endl;
			logger.open(fileNames[0], std::ios::app);
			logger << "resid  = " << residual << " delta  = " << delta << std::endl;
			logger << "trial  = " << trial << std::endl;
			printMatrixesABC();
			logger.close();
			break;
		}
	}
	for (int trial = 0; trial < 10; trial++) {
		nextApprox(delta, residual);
		logger.open(fileNames[0], std::ios::app);
		logger << "resid  = " << residual << " delta  = " << delta << std::endl;
		logger << "trial  = " << trial << std::endl;
		printMatrixesABC();
		logger.close();
	}
}

int main(int argc, char* argv[]) {
	int k = 0;
#ifndef _DEBUG
#ifndef PRECISE
	MPI_Init(&argc, &argv);                       /* starts MPI */
	MPI_Comm_rank(MPI_COMM_WORLD, &k); /* get current process id */
#endif
#endif

	srand(static_cast<unsigned int>(time(NULL)) + k*k * 1000);

#ifdef PRECISE
	sprintf_s(fileNames[0], "loggerPrecise_%dproc.txt", k);
	logger.open(fileNames[0]);
	logger << std::fixed;
	logger.precision(15);
	logger.close();
#else
	for (int i = 0; i < 8; i++) {
		sprintf_s(fileNames[i], "loggerFind%d_maxint%d_%dproc.txt", i, MAXMATRIXINT, k);
		logger.open(fileNames[i]);
		logger << std::fixed;
		logger.precision(15);
		logger.close();
	}
#endif

	setStaticFdelta();

	for (int i = 0; i < T; i++) {
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
#ifndef PRECISE
	MPI_Finalize();
#endif
#endif
	return(0);
}


