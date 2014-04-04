#pragma once

int const N = 3;
int const T = 3;
double const EPSMACH = 2.2e-15;

typedef std::complex<double> Complex;
typedef Complex ComplexArray[N][N];
typedef char FileNameString[50];

extern ComplexArray* a[T];
extern ComplexArray* b[T];
extern ComplexArray* c[T];

extern double fdelta[N][N][N][N][N][N];
extern Complex bcmul[T][N][N][N][N]; /* t, k, l, r, s pre-calculation b_tkl*c_trs */
extern double penalty[T][N][N][N][N][N][N];

struct IndexesIJKLRS {
	int i, j, k, l, r, s;
	int ii, jj, kk, ll, rr, ss; // 0 -> 0, 1->2, 2->3:  k->kk
	IndexesIJKLRS(int e = 0, int f = 0, int a = 0, int b = 0, int c = 0, int d = 0) : i(e), j(f), k(a), l(b), r(c), s(d) {
		ii = (2 * i) % N;
		jj = (2 * j) % N;
		kk = (2 * k) % N;
		ll = (2 * l) % N;
		rr = (2 * r) % N;
		ss = (2 * s) % N;
	}
};

int const BRUTEFORCESIZE = 122;
extern IndexesIJKLRS bruteForce[BRUTEFORCESIZE];