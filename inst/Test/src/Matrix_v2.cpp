#include "Matrix_v2.h"

double inverse(vector<vector<double> > A, vector<vector<double> > &Ainv) 
// 
// Computes the inverse of A, returns the determinant
//
{
  int n = A.size();
  int i;
  //  for (i = 0; i < n; i++)
  //    assert(n == A[i].size());
  
  Ainv.resize(n);
  for (i = 0; i < n; i++)
    Ainv[i].resize(n);

  vector<double> scale(n);
  vector<vector<double> > b;
  b.resize(n);
  for (i = 0; i < n; i++)
    b[i].resize(n);
  vector<int> index(n+1);

  // Initialize b to identy matrix
  int j;
  for (i = 0; i < n; i++)
    {
      for (j = 0; j < n; j++)
	b[i][j] = 0.0;
      b[i][i] = 1.0;
    }


  // Set scale factor, scale[i] = max( |a[i][j]| ), for each row
  for (i = 0; i < n; i++) 
    {
      index[i] = i;			  
      double scalemax = 0.;
      for(j = 0; j < n; j++) 
	scalemax = (scalemax > fabs(A[i][j])) ? scalemax : fabs(A[i][j]);
      scale[i] = scalemax;
    }

  // Loop over rows k = 0, ..., (n-2)
  int signDet = 1;
  int k;
  for (k = 0; k < n - 1; k++) 
    {
      //* Select pivot row from max( |a[j][k]/s[j]| )
      double ratiomax = 0.0;
      int jPivot = k;
      for (i = k; i < n; i++) 
	{
	  double ratio = fabs(A[index[i]][k])/scale[index[i]];
	  if(ratio > ratiomax) {
	    jPivot=i;
	    ratiomax = ratio;
	  }
	}

      //* Perform pivoting using row index list
      int indexJ = index[k];
      if (jPivot != k) {	            // Pivot
	indexJ = index[jPivot];
	index[jPivot] = index[k];   // Swap index jPivot and k
	index[k] = indexJ;
	signDet *= -1;		    // Flip sign of determinant
      }

      //* Perform forward elimination
      for (i = k+1; i < n; i++) 
	{
	  double coeff = A[index[i]][k]/A[indexJ][k];
	  for(j = k+1; j < n; j++)
	    A[index[i]][j] -= coeff*A[indexJ][j];
	  A[index[i]][k] = coeff;
	  for (j = 0; j < n; j++) 
	    b[index[i]][j] -= A[index[i]][k]*b[indexJ][j];
	}
    }
  
  //* Compute determinant as product of diagonal elements
  double determ = signDet;	   // Sign of determinant
  for (i = 0; i < n; i++)
    determ *= A[index[i]][i];
  
  //* Perform backsubstitution
  for (k=0; k < n; k++) 
    {
      Ainv[n-1][k] = b[index[n-1]][k]/A[index[n-1]][n-1];
      for (i = n-2; i >= 0; i--) 
	{
	  double sum = b[index[i]][k];
	  for (j=i+1; j < n; j++)
	    sum -= A[index[i]][j]*Ainv[j][k];
	  Ainv[i][k] = sum/A[index[i]][i];
	}
    }

  return determ;
}




double inverseLnDeterminant(vector<vector<double> > A, vector<vector<double> > &Ainv) 
// 
// Computes the inverse of A, returns ln(|determinant|)
//
{
  int n = A.size();
  int i;
  //  for (i = 0; i < n; i++)
  //    assert(n == A[i].size());
  
  Ainv.resize(n);
  for (i = 0; i < n; i++)
    Ainv[i].resize(n);

  vector<double> scale(n);
  vector<vector<double> > b;
  b.resize(n);
  for (i = 0; i < n; i++)
    b[i].resize(n);
  vector<int> index(n+1);

  // Initialize b to identy matrix
  int j;
  for (i = 0; i < n; i++)
    {
      for (j = 0; j < n; j++)
	b[i][j] = 0.0;
      b[i][i] = 1.0;
    }


  // Set scale factor, scale[i] = max( |a[i][j]| ), for each row
  for (i = 0; i < n; i++) 
    {
      index[i] = i;			  
      double scalemax = 0.;
      for(j = 0; j < n; j++) 
	scalemax = (scalemax > fabs(A[i][j])) ? scalemax : fabs(A[i][j]);
      scale[i] = scalemax;
    }

  // Loop over rows k = 0, ..., (n-2)
  int signDet = 1;
  int k;
  for (k = 0; k < n - 1; k++) 
    {
      //* Select pivot row from max( |a[j][k]/s[j]| )
      double ratiomax = 0.0;
      int jPivot = k;
      for (i = k; i < n; i++) 
	{
	  double ratio = fabs(A[index[i]][k])/scale[index[i]];
	  if(ratio > ratiomax) {
	    jPivot=i;
	    ratiomax = ratio;
	  }
	}

      //* Perform pivoting using row index list
      int indexJ = index[k];
      if (jPivot != k) {	            // Pivot
	indexJ = index[jPivot];
	index[jPivot] = index[k];   // Swap index jPivot and k
	index[k] = indexJ;
	signDet *= -1;		    // Flip sign of determinant
      }

      //* Perform forward elimination
      for (i = k+1; i < n; i++) 
	{
	  double coeff = A[index[i]][k]/A[indexJ][k];
	  for(j = k+1; j < n; j++)
	    A[index[i]][j] -= coeff*A[indexJ][j];
	  A[index[i]][k] = coeff;
	  for (j = 0; j < n; j++) 
	    b[index[i]][j] -= A[index[i]][k]*b[indexJ][j];
	}
    }
  
  //* Compute determinant as product of diagonal elements
  double determ = 0.0;
  for (i = 0; i < n; i++)
    determ += log(abs(A[index[i]][i]));
  
  //* Perform backsubstitution
  for (k=0; k < n; k++) 
    {
      Ainv[n-1][k] = b[index[n-1]][k]/A[index[n-1]][n-1];
      for (i = n-2; i >= 0; i--) 
	{
	  double sum = b[index[i]][k];
	  for (j=i+1; j < n; j++)
	    sum -= A[index[i]][j]*Ainv[j][k];
	  Ainv[i][k] = sum/A[index[i]][i];
	}
    }

  return determ;
}




void matrixMult(const vector<vector<double> > &A,const vector<vector<double> > &B,
		vector<vector<double> > &C)
{
  int m = A.size();
  int n = A[0].size();
  int k = B[0].size();

  int i,j,r;
  //  for (i = 0; i < m; i++)
  //    assert(A[i].size() == n);
  //  assert(B.size() == n);
  //  for (i = 0; i < n; i++)
  //    assert(B[i].size() == k);

  C.resize(m);
  for (i = 0; i < m; i++)
    C[i].resize(k);

  for (i = 0; i < m; i++)
    for (j = 0; j < k; j++)
      {
	C[i][j] = 0.0;
	for (r = 0; r < n; r++)
	  C[i][j] += A[i][r] * B[r][j];
      }

  return;
}


void matrixMult(const vector<vector<double> > &A,const vector<double> &B,
		vector<double> &C)
{
  int m = A.size();
  int n = A[0].size();
  int k = 1;

  int i,j,r;
  //  for (i = 0; i < m; i++)
  //    assert(A[i].size() == n);
  //  assert(B.size() == n);

  C.resize(m);

  for (i = 0; i < m; i++)
    {
      C[i] = 0.0;
      for (r = 0; r < n; r++)
	C[i] += A[i][r] * B[r];
    }

  return;
}




void outerProduct(const vector<vector<double> > &A,vector<vector<double> > &AAT)
{
  int m = A.size();
  int n = A[0].size();
  int i;
  //  for (i = 0; i < m; i++)
  //    assert(A[i].size() == n);

  AAT.resize(m);
  for (i = 0; i < m; i++)
    AAT[i].resize(m);

  int j,k;
  for (i = 0; i < m; i++)
    for (j = 0; j < m; j++)
      {
	AAT[i][j] = 0.0;
	for (k = 0; k < n; k++)
	  AAT[i][j] += A[i][k] * A[j][k];
      }

  return;
}




double quadratic(const vector<vector<double> > &A,const vector<double> &x)
{
  int n = A.size();
  int i;
  //  for (i = 0; i < n; i++) 
  //    assert(A[i].size() == n);
  //  assert(x.size() == n);

  double sum = 0.0;
  int j;
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      sum += A[i][j] * x[i] * x[j];

  return sum;
}



void quadratic(const vector<vector<double> > &A,const vector<vector<double> > &B,vector<vector<double> > &C) {
  // computes C = A^T * B * A

  int n = B.size();
  int m = A[0].size();

  C.resize(m);
  int i;
  for (i = 0; i < m; i++)
    C[i].resize(m);

  int j;
  for (i = 0; i < m; i++)
    for (j = 0; j < m; j++) {
      C[i][j] = 0.0;
      int k,l;
      for (k = 0; k < n; k++)
	for (l = 0; l < n; l++)
	  C[i][j] += A[k][i] * B[k][l] * A[l][j];
    }
  
  return;
}




void quadratic2(const vector<vector<double> > &A,const vector<vector<double> > &B,vector<vector<double> > &C) {
  // computes C = A * B * A^T

  int n = B.size();
  int m = A.size();

  C.resize(m);
  int i;
  for (i = 0; i < m; i++)
    C[i].resize(m);

  int j;
  for (i = 0; i < m; i++)
    for (j = 0; j < m; j++) {
      C[i][j] = 0.0;
      int k,l;
      for (k = 0; k < n; k++)
	for (l = 0; l < n; l++)
	  C[i][j] += A[i][k] * B[k][l] * A[j][l];
    }
  
  return;
}




