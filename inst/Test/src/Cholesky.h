#ifndef CHOLESKY_H
#define CHOLESKY_H

#include "Matrix.h"

				// A class that holds L - the cholesky
				// decomposition of a square, positive
                                // definite matrix
class Cholesky
{
  const int N;			        // Dimensions of the problem (N*N)
  vector<vector<double> > L;		// L*L lower diagonal matrix L
  
public:
  Cholesky(const vector<vector<double> > &A,int &err);		// Decompose Matrix A, of N*N
  ~Cholesky(void);
   
  const vector<vector<double> > &q_L(void) const  {return L;};
  int q_nrows(void) const               { return L.size();};
  int q_ncols(void) const               { return L[0].size();};
        
};




class Cholesky_solve
{
  vector<double> result;
  
 public:
  Cholesky_solve(const Cholesky &chol,const vector<double> &b);
  const vector<double> &solution(void) const    {return result;};
};



#endif
