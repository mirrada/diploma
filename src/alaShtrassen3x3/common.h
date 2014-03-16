#pragma once

int const N = 3;
int const T = 3;

typedef std::complex<double> Complex;
typedef Complex ComplexArray[N][N];

static ComplexArray* a[T];
static ComplexArray* b[T];
static ComplexArray* c[T];

static double fdelta[N][N][N][N][N][N] = { 0 };
static Complex bcmul[T][N][N][N][N] = { 0 }; /* t, k, l, r, s pre-calculation b_tkl*c_trs */


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