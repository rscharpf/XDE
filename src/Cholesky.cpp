#include <iostream>
#include <stdio.h>
#include <stdlib.h>

#include "Cholesky.h"


//
// FUNCTION: Cholesky::Cholesky
//
// PURPOSE: Constructor
//
Cholesky::Cholesky(const vector<vector<double> > &A,int &err) : N(A.size())
{
  L.resize(N);
  int i;
  for (i = 0; i < N; i++)
    L[i].resize(N);

  err = 0;
  if(N != A[0].size())
    {
      cout << "Cholesky: Matrix must be square !" << "\n";
      exit(-1);
    }
  
  int j;
  for (i = 0; i < N; i++)
    for (j = 0; j < N; j++)
      L[i][j] = 0.0;
  
  int k;
  double sum;
  
  for (i = 0; i < N; i++)
    {
      for (j = i; j < N; j++)
	{
	  for (sum = A[i][j],k=i-1; k >= 0; k--)
	    sum -= L[i][k] * L[j][k];
	  
	  if (i == j && sum <= 0.0)
	    {
	      err = 1;
	      //	      fprintf(stderr,"Cholesky: Matrix is not positive definite !\n");
	    }
	  
	  if (i == j)
	    L[j][i] = sqrt(sum);
	  else
	    L[j][i] = sum / L[i][i];
	}
    }
}



//
// FUNCTION: Cholesky::~Cholesky
//
// PURPOSE: Destructor
//
Cholesky::~Cholesky(void)
{
  return;
}







//
// FUNCTION: Cholesky_solve::Cholesky_solve
//
// PURPOSE: Constructor
//
Cholesky_solve::Cholesky_solve(const Cholesky &chol,const vector<double> &b) :
  result(b.size())
{
  const vector<vector<double> > L(chol.q_L());
  //  assert(L.size() == L[0].size());
  //  assert(L.size() == b.size());
  
  int n = b.size();
  
  int i,k;
  double sum;
  for (i = 0; i < n; i++)
    {
      for (sum = b[i], k = i - 1; k >= 0; k--)
	sum -= L[i][k] * result[k];
      result[i] = sum / L[i][i];
    }
  

  for (i = n - 1; i >= 0; i--)
    {
      for (sum = result[i], k = i + 1; k < n; k++)
	sum -= L[k][i] * result[k];
      result[i] = sum / L[i][i];
    }
  
  return;
}


  

