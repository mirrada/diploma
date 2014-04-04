#include "stdafx.h"
#include "fmSolver.h"


const double PI = 3.141592653589793;


Complex temp;
IndexesIJKLRS tempIJKLRS;
int tempT;
Complex abcArr[T][N][N][N][N][N][N];
double maxInRow[BRUTEFORCESIZE];
double MYINF = 1e+10;
double GREATERTHAN0 = 1e-5;


lprec* lp;
int const COLUMNS = T*N*N * 3;

std::ofstream logger1;

FileNameString fileName = "test.txt";


void prepareABC(){
	int i, j, k, l, r, s;
	for (int n = 0; n < BRUTEFORCESIZE; n++) {
		i = bruteForce[n].i;
		j = bruteForce[n].j;
		k = bruteForce[n].k;
		l = bruteForce[n].l;
		r = bruteForce[n].r;
		s = bruteForce[n].s;
		for (int t = 0; t < T; t++)
			abcArr[t][i][j][k][l][r][s] = (*a[t])[i][j] * bcmul[t][k][l][r][s];
		i = bruteForce[n].ii;
		j = bruteForce[n].jj;
		k = bruteForce[n].kk;
		l = bruteForce[n].ll;
		r = bruteForce[n].rr;
		s = bruteForce[n].ss;
		for (int t = 0; t < T; t++)
			abcArr[t][i][j][k][l][r][s] = (*a[t])[i][j] * bcmul[t][k][l][r][s];
	}
}

void prepareMaxInRow(){
	int i, j, k, l, r, s, ii, jj, kk, ll, rr, ss;
	for (int n = 0; n < BRUTEFORCESIZE; n++) {
		i = bruteForce[n].i;
		j = bruteForce[n].j;
		k = bruteForce[n].k;
		l = bruteForce[n].l;
		r = bruteForce[n].r;
		s = bruteForce[n].s;
		ii = bruteForce[n].ii;
		jj = bruteForce[n].jj;
		kk = bruteForce[n].kk;
		ll = bruteForce[n].ll;
		rr = bruteForce[n].rr;
		ss = bruteForce[n].ss;
		maxInRow[n] = abs(abcArr[0][i][j][k][l][r][s]);
		for (int t = 0; t < T; t++) {
			maxInRow[n] = max(maxInRow[n], abs(abcArr[t][i][j][k][l][r][s]));
			maxInRow[n] = max(maxInRow[n], abs(abcArr[t][ii][jj][kk][ll][rr][ss]));
		}
	}
}

void makeRow2(int t, IndexesIJKLRS index, int t1, IndexesIJKLRS index1, bool greater) {
	double sparserow[6] = { 1, 1, 1, -1, -1, -1 }; /* must be the number of non-zero values */
	int colno[6];
	int type;
	double rhs;

	colno[0] = 1 + t*N*N + index.i*N + index.j;
	colno[1] = 1 + T*N*N + t*N*N + index.k*N + index.l;
	colno[2] = 1 + 2 * T*N*N + t*N*N + index.r*N + index.s;
	colno[3] = 1 + t1*N*N + index1.i*N + index1.j;
	colno[4] = 1 + T*N*N + t1*N*N + index1.k*N + index1.l;
	colno[5] = 1 + 2 * T*N*N + t1*N*N + index1.r*N + index1.s;

	if (greater) {
		type = GE;
		rhs = GREATERTHAN0;
	}
	else {
		type = EQ;
		rhs = 0;
	}

	add_constraintex(lp, 6, sparserow, colno, type, rhs);
}

void makeRow1(int t, IndexesIJKLRS index, bool greater) {
	double sparserow[3] = { 1, 1, 1 }; /* must be the number of non-zero values */
	int colno[3];
	int type;
	double rhs;

	colno[0] = 1 + t*N*N + index.i*N + index.j;
	colno[1] = 1 + T*N*N + t*N*N + index.k*N + index.l;
	colno[2] = 1 + 2 * T*N*N + t*N*N + index.r*N + index.s;

	if (greater) {
		type = GE;
		rhs = GREATERTHAN0;
	}
	else {
		type = EQ;
		rhs = 0;
	}

	add_constraintex(lp, 3, sparserow, colno, type, rhs);
}

void findFuncApproxSolutionRec(int n);

void recOneRow(int n, int m, int cond, bool has0) {
	int i, j, k, l, r, s, t;
	Complex abc;
	if (m < 3) {
		i = bruteForce[n].i;
		j = bruteForce[n].j;
		k = bruteForce[n].k;
		l = bruteForce[n].l;
		r = bruteForce[n].r;
		s = bruteForce[n].s;
		t = m;
		abc = abcArr[t][i][j][k][l][r][s];
	}
	else {
		i = bruteForce[n].ii;
		j = bruteForce[n].jj;
		k = bruteForce[n].kk;
		l = bruteForce[n].ll;
		r = bruteForce[n].rr;
		s = bruteForce[n].ss;
		t = m % 3;
#ifdef GROUP2
		if (!i ^ !j ^ !k ^ !l ^ !r ^ !s)
			abc = -abcArr[t][i][j][k][l][r][s];
		else
#endif
			abc = abcArr[t][i][j][k][l][r][s];
	}

	if (m >= 6) {
		if (abs(fdelta[i][j][k][l][r][s]) < EPSMACH) { // sum=0.000
			if (cond != 1)
				findFuncApproxSolutionRec(n + 1);
		}
		else { // sum=0.333
			if (has0 && cond != 1)
				findFuncApproxSolutionRec(n + 1);
		}
		return;
	}


	if (abs(abc) <= EPSMACH) { // = 0
		recOneRow(n, m + 1, cond, has0);
	}
	else {
		if (abs(abc) < 1 + EPSMACH) { // |x| < 1
			makeRow1(t, IndexesIJKLRS(i, j, k, l, r, s), true);
			recOneRow(n, m + 1, cond, has0);
			//set_add_rowmode(lp, false);
			del_constraint(lp, get_Nrows(lp));
			//set_add_rowmode(lp, true);
		}

		if (abs(fdelta[i][j][k][l][r][s]) < EPSMACH)  {// sum = 0.0000
			if (abs(abc) > 1 - EPSMACH)  // |x| > 1
				if (cond == 0) {
					temp = abc;
					tempT = t;
					tempIJKLRS = IndexesIJKLRS(i, j, k, l, r, s);
					recOneRow(n, m + 1, 1, has0);
				}
		}
		else { // sum = +-0.333333
			if (!has0)
				if (maxInRow[n] / abs(abc) < 10) {
					has0 = true;
					makeRow1(t, IndexesIJKLRS(i, j, k, l, r, s), false);
					recOneRow(n, m + 1, cond, true);
					//set_add_rowmode(lp, false);
					del_constraint(lp, get_Nrows(lp));
					//set_add_rowmode(lp, true);
				}
			if (abs(abc) > 5e-1)
				if (cond == 0) {
					temp = abc;
					tempT = t;
					tempIJKLRS = IndexesIJKLRS(i, j, k, l, r, s);
					//recOneRow(n, m + 1, 1, true);
					recOneRow(n, m + 1, 1, false);
				}
		}

	}

	if (cond == 1) {
		double diff = abs(temp) / abs(abc);
		if (diff < 10 && diff > 0.10)
			if (temp.imag()*abc.imag() < EPSMACH &&
				temp.real()*abc.real() < EPSMACH) { // from different sectors -> <-
				makeRow2(t, IndexesIJKLRS(i, j, k, l, r, s), tempT, tempIJKLRS, false);
				recOneRow(n, m + 1, 2, has0);
				//set_add_rowmode(lp, false);
				del_constraint(lp, get_Nrows(lp));
				//set_add_rowmode(lp, true);
			}
	}

}

void printMatrix(char* desc, int z, int n, int m, int* a) {
	logger1 << desc << z << std::endl;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++)
			logger1 << a[i*m + j] << "  ";
		logger1 << std::endl;
	}
}

void printOneRow(int t, int i, int j, int k, int l, int r, int s, double v[COLUMNS]){

	if ((abs((*a[t])[i][j]) < EPSMACH && abs((*a[t])[i][j]) > -EPSMACH) ||
		(abs((*b[t])[k][l]) < EPSMACH && abs((*b[t])[k][l]) > -EPSMACH) ||
		(abs((*c[t])[r][s]) < EPSMACH && abs((*c[t])[r][s]) > -EPSMACH)) {
		logger1 << "0";
	}
	else {
		int ia = t*N*N + i*N + j;
		int ib = T*N*N + t*N*N + k*N + l;
		int ic = 2 * T*N*N + t*N*N + r*N + s;
		logger1 << "x^(" << v[ia] << "+" << v[ib] << "+" << v[ic] << "=" << v[ia] + v[ib] + v[ic] << ")";
	}
}

void printRationalSolution(double v[COLUMNS]) {
	int i, j, k, l, r, s, ii, jj, kk, ll, rr, ss;
	logger1.open(fileName, std::ios::app);
	logger1 << std::fixed;
	logger1.precision(10);
	for (int n = 0; n < BRUTEFORCESIZE; n++) {
		i = bruteForce[n].i;
		j = bruteForce[n].j;
		k = bruteForce[n].k;
		l = bruteForce[n].l;
		r = bruteForce[n].r;
		s = bruteForce[n].s;
		ii = bruteForce[n].ii;
		jj = bruteForce[n].jj;
		kk = bruteForce[n].kk;
		ll = bruteForce[n].ll;
		rr = bruteForce[n].rr;
		ss = bruteForce[n].ss;
		logger1 << i << j << k << l << r << s << std::endl;
		for (int t = 0; t < T; t++){
			logger1 << " + ";
			printOneRow(t, i, j, k, l, r, s, v);

#ifdef GROUP2
			if (!i ^ !j ^ !k ^ !l ^ !r ^ !s)
				logger1 << " - ";
			else
#endif
				logger1 << " + ";
			printOneRow(t, ii, jj, kk, ll, rr, ss, v);
		}
		logger1 << " = " << fdelta[i][j][k][l][r][s] << std::endl;
	}
	logger1.close();
}


int count = 0;
void findFuncApproxSolutionRec(int n) {
	if (n < BRUTEFORCESIZE) {
		recOneRow(n, 0, 0, false);
		//std::cout << "n " << n << std::endl;
	}
	else {
		count++;

		std::cout << "count " << count << std::endl;

		//set_add_rowmode(lp, false);
		//print_lp(lp);

		int info = solve(lp);
		if (info != 2) {
			print_solution(lp, 81);
			double vars[COLUMNS];
			get_variables(lp, vars);
			printRationalSolution(vars);
		}
		//set_add_rowmode(lp, true);
	}
}


void managetIJKLRS(int t, int i, int j, int k, int l, int r, int s){
	int ia = 1 + t*N*N + i*N + j;
	int ib = 1 + T*N*N + t*N*N + k*N + l;
	int ic = 1 + 2 * T*N*N + t*N*N + r*N + s;
	if (abs((*a[t])[i][j]) < EPSMACH && abs((*a[t])[i][j]) > -EPSMACH) {
		set_obj(lp, ia, 0);
		set_bounds(lp, ia, 0, 0);
	}
	else
		set_bounds(lp, ia, -MYINF, MYINF);
	if (abs((*b[t])[k][l]) < EPSMACH && abs((*b[t])[k][l]) > -EPSMACH){
		set_obj(lp, ib, 0);
		set_bounds(lp, ib, 0, 0);
	}
	else
		set_bounds(lp, ib, -MYINF, MYINF);
	if (abs((*c[t])[r][s]) < EPSMACH && abs((*c[t])[r][s]) > -EPSMACH){
		set_obj(lp, ic, 0);
		set_bounds(lp, ic, 0, 0);
	}
	else
		set_bounds(lp, ic, -MYINF, MYINF);
	std::string name;
	name = std::to_string(t) + "a" + std::to_string(i) + std::to_string(j);
	set_col_name(lp, ia, const_cast<char*>(name.c_str()));
	name = std::to_string(t) + "b" + std::to_string(k) + std::to_string(l);
	set_col_name(lp, ib, const_cast<char*>(name.c_str()));
	name = std::to_string(t) + "c" + std::to_string(r) + std::to_string(s);
	set_col_name(lp, ic, const_cast<char*>(name.c_str()));
}

void setBoundedNamesObj() {
	double row[COLUMNS + 1];
	int i, j, k, l, r, s, ii, jj, kk, ll, rr, ss;
	for (int i = 1; i < COLUMNS + 1; i++) {
		row[i] = 1.0;
	}
	set_obj_fn(lp, row);
	for (int n = 0; n < BRUTEFORCESIZE; n++) {
		i = bruteForce[n].i;
		j = bruteForce[n].j;
		k = bruteForce[n].k;
		l = bruteForce[n].l;
		r = bruteForce[n].r;
		s = bruteForce[n].s;
		ii = bruteForce[n].ii;
		jj = bruteForce[n].jj;
		kk = bruteForce[n].kk;
		ll = bruteForce[n].ll;
		rr = bruteForce[n].rr;
		ss = bruteForce[n].ss;
		for (int t = 0; t < T; t++){
			managetIJKLRS(t, i, j, k, l, r, s);
			managetIJKLRS(t, ii, jj, kk, ll, rr, ss);
		}
	}
}


void findFuncApproxSolution() {
	count = 0;
	lp = make_lp(0, COLUMNS);
	set_minim(lp);
	set_outputfile(lp, fileName);
	setBoundedNamesObj();
	//set_add_rowmode(lp, true);
	//print_lp(lp);
	prepareABC();
	prepareMaxInRow();
	findFuncApproxSolutionRec(0);
	delete_lp(lp);
}