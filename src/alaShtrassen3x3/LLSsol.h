#pragma once
#include "common.h"
#include <mkl.h>

extern ComplexArray* a[T];
extern ComplexArray* b[T];
extern ComplexArray* c[T];

extern double fdelta[N][N][N][N][N][N];
extern Complex bcmul[T][N][N][N][N]; /* t, k, l, r, s pre-calculation b_tkl*c_trs */

extern double penalty[T][N][N][N][N][N][N];

double prepareLLSSolPart(); 
