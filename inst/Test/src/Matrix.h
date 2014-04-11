#ifndef MATRIX_H
#define MATRIX_H

using namespace std;

#include <cmath>
#include <vector>


double inverse(vector<vector<double> > A, vector<vector<double> > &Ainv);
void matrixMult(const vector<vector<double> > &A,const vector<vector<double> > &B,
		vector<vector<double> > &C);
void matrixMult(const vector<vector<double> > &A,const vector<double> &B,
		vector<double> &C);
void outerProduct(const vector<vector<double> > &A,vector<vector<double> > &AAT);
double quadratic(const vector<vector<double> > &A,const vector<double> &x);



#endif
