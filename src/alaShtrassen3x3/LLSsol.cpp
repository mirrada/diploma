#include "stdafx.h"
#include "LLSsol.h"


int MAXINSOLUTION = 1000000;

int const NROWpart6 = 27;
int const NROWpart3 = 14;
int const addition13_3 = 0;
int const addition13_6 = 0;


Complex LLST6[NROWpart6 + 6 * NROWpart6 + addition13_6][6];
Complex bLLS6[NROWpart6 + 6 * NROWpart6 + addition13_6];

Complex LLST3[NROWpart3 + 6 * NROWpart3 + addition13_3][3];
Complex bLLS3[NROWpart3 + 6 * NROWpart3 + addition13_3];



/*int col = i*N*N + j*N;  int shift = ii*N*N + jj*N;*/
template<unsigned int ROWS, unsigned int COLS>
void fillLLSsub(int i, int j, int notIJ00shift, Complex(&arr)[ROWS][COLS], Complex(&b)[ROWS]) {
	int rowInit = 0;
	int rows = ROWS;
	int cols = COLS;
	//int addition13 = 0;
	for (int k = 0; k < (notIJ00shift ? N : 2); k++) {
		for (int l = 0; l < (notIJ00shift ? N : k + 2); l++)
			for (int r = 0; r < (notIJ00shift ? N : N - !(k + l)); r++){
				int s = (i + k + r + 2 * l + 2 * j) % N;
				int ii = (2 * i) % N;
				int jj = (2 * j) % N;
				int kk = (2 * k) % N;
				int ll = (2 * l) % N;
				int rr = (2 * r) % N;
				int ss = (2 * s) % N;
				int row = (rowInit)*(6 + 1);
				rowInit++;
				for (int t = 0; t < T; t++) {
					arr[row][t] += bcmul[t][k][l][r][s];
#ifdef GROUP2
					if (!i ^ !j ^ !k ^ !l ^ !r ^ !s)
						arr[row][notIJ00shift + t] -= bcmul[t][kk][ll][rr][ss];
					else
#endif
						arr[row][notIJ00shift + t] += bcmul[t][kk][ll][rr][ss];
					b[row] = fdelta[i][j][k][l][r][s];

					arr[row + t + 1][t] = penalty[t][i][j][k][l][r][s] * bcmul[t][k][l][r][s];
					b[row + t + 1] = 0;
					if (i + j + k + l + r + s != 0) {
						arr[row + t + 1 + T][notIJ00shift + t] = penalty[t][ii][jj][kk][ll][rr][ss] * bcmul[t][kk][ll][rr][ss];
						b[row + t + 1 + T] = 0;
					}
				}/*
				if (fdelta[i][j][k][l][r][s] > EPSMACH || fdelta[i][j][k][l][r][s] < -EPSMACH) {
					addition13++;
					memcpy(arr[rows - addition13], arr[row], sizeof(arr[row]));
					b[rows - addition13] = b[row];
				}*/
			}
	}
}

template<unsigned int ROWS, unsigned int COLS>
double prepareLLSSolOnePart(int i, int j, int shift, Complex(&arr)[ROWS][COLS], Complex(&b)[ROWS]){
	memset(arr, 0, sizeof(arr));	
	memset(b, 0, sizeof(b));

	fillLLSsub(i, j, shift, arr, b);

	int info = LAPACKE_zgels(LAPACK_ROW_MAJOR, 'N', ROWS, COLS, 1,
		(lapack_complex_double*)arr, COLS,
		(lapack_complex_double*)b, 1);

	if (info)
		return -1;
	for (unsigned int k = 0; k < COLS; k++)
		if (abs(b[k]) > MAXINSOLUTION)
			b[k] = 0;

	double resid = 0;
	for (unsigned int k = COLS; k < ROWS; k++)
		resid += std::norm(b[k]);

	for (int t = 0; t < T; t++) {
		(*a[t])[i][j] = b[t];
		(*a[t])[2 * i%N][2 * j%N] = b[t + shift];
	}

	return resid;
}

double prepareLLSSolPart(){
	double r00 = prepareLLSSolOnePart(0, 0, 0, LLST3, bLLS3);
	double r01 = prepareLLSSolOnePart(0, 1, T, LLST6, bLLS6);
	double r10 = prepareLLSSolOnePart(1, 0, T, LLST6, bLLS6);
	double r11 = prepareLLSSolOnePart(1, 1, T, LLST6, bLLS6);
	double r21 = prepareLLSSolOnePart(2, 1, T, LLST6, bLLS6);
	if (r00 < 0 || r01 < 0 || r10 < 0 || r11 < 0 || r21 < 0)
		return -1;
	else
		return (r00)+(r01)+(r10)+(r11)+(r21);
}
