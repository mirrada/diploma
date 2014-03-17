#include "stdafx.h"
#include "LEG1sol.h"



Complex LS0011[T][T]; /*LESystem for i=0,j=0 and at00+at11*/
Complex b0011[T][SIZEB];  /*right hand vecors  for i=0,j=0 (the first column) and a00+a11 (the second column)*/

Complex LS22[T][T]; /*LESystem for at00-at11*/
Complex b22[T]; /*right hand vecors at00-at11*/

Complex LSIJ[LSNUM][ROWSBIGLS][ROWSBIGLS];  /*LESystem for i!=j*/
Complex bIJ[LSNUM][ROWSBIGLS];   /*right hand vecor i!=j*/

IndexesKLRS  bruteForce[N][N*N*N];

NumerationLSVariables numIJ[3];

NumerationLSVariables::NumerationLSVariables(bool(*condition)(int, int)) {
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

void setStaticNumerationLSVariables() {
	auto condIJ01 = [](int i, int j){return (i == 0 && j != 0); };
	auto condIJ10 = [](int i, int j){return (i != 0 && j == 0); };
	auto condIJ21 = [](int i, int j){return (i != j && i * j != 0); };
	numIJ[0] = NumerationLSVariables(condIJ01);
	numIJ[1] = NumerationLSVariables(condIJ10);
	numIJ[2] = NumerationLSVariables(condIJ21);
}

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

void getAijFromSolution() {
	for (int t = 0; t < T; t++) {
		(*a[t])[0][0] = b0011[t][0];
		(*a[t])[1][1] = (b0011[t][1] + b22[t]) / 2.0;
		(*a[t])[2][2] = (b0011[t][1] - b22[t]) / 2.0;
	}
	for (int lsnum = 0; lsnum < LSNUM; lsnum++)
		for (int m = 0; m < ROWSBIGLS; m++) {
			auto d = numIJ[lsnum].getIndexes(m);
			(*a[d.t])[d.i][d.j] = bIJ[lsnum][m];
		}
}

int prepareLEsol() {
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

	getAijFromSolution();

	return 0;
}