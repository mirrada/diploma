#pragma once

int const N = 3;
int const T = 3;

typedef std::complex<double> Complex;
typedef Complex ComplexArray[N][N];
typedef char FileNameString[50];

extern ComplexArray* a[T];
extern ComplexArray* b[T];
extern ComplexArray* c[T];

extern double fdelta[N][N][N][N][N][N];
extern Complex bcmul[T][N][N][N][N]; /* t, k, l, r, s pre-calculation b_tkl*c_trs */
extern double penalty[T][N][N][N][N][N][N];