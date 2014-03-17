#pragma once
#include "common.h"
#include <mkl.h>

extern ComplexArray* a[T];
extern ComplexArray* b[T];
extern ComplexArray* c[T];

extern double fdelta[N][N][N][N][N][N];
extern Complex bcmul[T][N][N][N][N]; /* t, k, l, r, s pre-calculation b_tkl*c_trs */

int const SIZEB = 2; /*columns in right hand matrix b. first for system a_t00, second for a_t11+a_t22*/
int const LSNUM = 3; /*number of systems for i!=j*/
int const ROWSBIGLS = 6; /*size of systems for i!=j*/

struct IndexesKLRS {
	int k, l, r, s;
	int kk, ll, rr, ss; // 0 -> 0, 1->2, 2->3:  k->kk
	IndexesKLRS(int a = 0, int b = 0, int c = 0, int d = 0) : k(a), l(b), r(c), s(d) {
		kk = (2* k) % N;
		ll = (2* l) % N;
		rr = (2* r) % N;
		ss = (2* s) % N;
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
	NumerationLSVariables(){};
	NumerationLSVariables(bool(*condition)(int i, int j));
	int getOrder(int i, int j, int t) { return toLS[i][j][t]; }
	IndexesIJT getIndexes(int m) { return toA[m]; }
};

void setStaticNumerationLSVariables();

void setStaticBruteForce();

int prepareLEsol();