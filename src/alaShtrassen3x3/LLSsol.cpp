#include "stdafx.h"
#include "LLSsol.h"


int MAXINSOLUTION = 1000000;

int const NROWpart = 27;
int const NROW = 5 * NROWpart;

Complex LLST[NROW][N*N*T];
Complex bLLS[NROW];

Complex LLST6[NROWpart][6];
Complex LLST3[NROWpart][3];
Complex bLLSpart[NROWpart];

/*int col = i*N*N + j*N;  int colcol = ii*N*N + jj*N;*/
template<unsigned int ROW, unsigned int COL>
void fillLLSsub(int i, int j, int rowStart, int col, int colcol, Complex(&arr)[ROW][COL], Complex(&b)[ROW]) {
	BruteForceKLRSstart{
		int row = rowStart + k*N*N + l*N + r;
		for (int t = 0; t < T; t++) {
			double one = 1;
			arr[row][col + t] += bcmul[t][k][l][r][s];
#ifdef GROUP2
			if (!i ^ !j ^ !k ^ !l ^ !r ^ !s)
				one = -1;
#endif
			arr[row][colcol + t] += one * bcmul[t][kk][ll][rr][ss];
			b[row] = fdelta[i][j][k][l][r][s];
		}
	}BruteForceKLRSend
}

template<unsigned int ROW, unsigned int COL>
int prepareLLSSolOnePart(int i, int j, int colcol, Complex(&arr)[ROW][COL]){
	memset(arr, 0, sizeof(arr));

	fillLLSsub(i, j, 0, 0, colcol, arr, bLLSpart);
	int info = LAPACKE_zgels(LAPACK_ROW_MAJOR, 'N', ROW, COL, 1,
		(lapack_complex_double*)arr, COL,
		(lapack_complex_double*)bLLSpart, 1);

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
	return prepareLLSSolOnePart(0, 0, 0, LLST3)||
	(prepareLLSSolOnePart(0, 1, T, LLST6))||
	(prepareLLSSolOnePart(1, 0, T, LLST6))||
	(prepareLLSSolOnePart(1, 1, T, LLST6))||
	prepareLLSSolOnePart(2, 1, T, LLST6);
}

void fillLLS(){
	memset(LLST, 0, sizeof(LLST));
	fillLLSsub(0, 0, 0, 0, 0, LLST, bLLS);
	fillLLSsub(0, 1, 1 * N*N*N, 3, 6, LLST, bLLS);
	fillLLSsub(1, 0, 2 * N*N*N, 9, 18, LLST, bLLS);
	fillLLSsub(1, 1, 3 * N*N*N, 12, 24, LLST, bLLS);
	fillLLSsub(2, 1, 4 * N*N*N, 21, 15, LLST, bLLS);
}

int prepareLLSSol(){
	fillLLS();

	int info = LAPACKE_zgels(LAPACK_ROW_MAJOR, 'N', NROW, N*N*T, 1, (lapack_complex_double*)LLST, N*N*T,
		(lapack_complex_double*)bLLS, 1);

	if (info)
		return info;

	for (int i = 0; i < N*N*T; i++)
		if (abs(bLLS[i]) > MAXINSOLUTION)
			bLLS[i] = 0;

	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
			for (int t = 0; t < T; t++)
				(*a[t])[i][j] = bLLS[i*N*N + j*N + t];
	return 0;
}