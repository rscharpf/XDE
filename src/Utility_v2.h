#ifndef INDEX_V2_H
#define INDEX_V2_H

#include "Matrix.h"

inline int qg2index(int q,int g,int Q,int G) {
  int index = g * Q + q;
  
  return index;
}


inline int qs2index(int q,int s,int Q,const int *S) {
  int index = 0;
  int p;
  for (p = 0; p < q; p++)
    index += S[p];
  
  index += s;
  
  return index;
}


inline int qgs2index(int q,int g,int s,int Q,int G,const int *S) {
  int index = 0;

  

  return index;
}


inline int qq2index(int q1,int q2,int Q) {
  int index = q1 * Q + q2;

  return index;
}


inline void makeSigma(std::vector<std::vector<double> > &Sigma,const int Q,
		      const double gamma2,const double *tau2,
		      const double *a,const double *sigma2g,
		      const double *r) {
  Sigma.resize(Q);
  int q;
  for (q = 0; q < Q; q++) {
    Sigma[q].resize(Q);
    Sigma[q][q] = gamma2 * tau2[q];
    Sigma[q][q] *= exp(a[q] * log(sigma2g[q]));
  }

  int q1,q2;
  for (q1 = 0; q1 < Q; q1++)
    for (q2 = q1 + 1; q2 < Q; q2++) {
      int k = qq2index(q1,q2,Q);
      Sigma[q1][q2] = r[k] * sqrt(Sigma[q1][q1] * Sigma[q2][q2]);

      Sigma[q2][q1] = Sigma[q1][q2];
    }

  return;
}



inline void makeSigma(std::vector<std::vector<double> > &Sigma,
		      const vector<int> &on,const int Q,
		      const double gamma2,const double *tau2,
		      const double *a,const double *sigma2g,
		      const double *r) {
  int dim = 0;
  int q; 
  for (q = 0; q < Q; q++)
    dim += on[q];

  Sigma.resize(dim);
  int k = 0;
  for (q = 0; q < Q; q++) {
    if (on[q] == 1) {
      Sigma[k].resize(dim);
      Sigma[k][k] = gamma2 * tau2[q];
      Sigma[k][k] *= exp(a[q] * log(sigma2g[q]));
      k++;
    }
  }

  int q1,q2;
  int k1 = 0;
  for (q1 = 0; q1 < Q; q1++) {
    if (on[q1] == 1) {
      int k2 = 0;
      for (q2 = 0; q2 < Q; q2++) {
	if (on[q2] == 1) {
	  if (k1 != k2) {
	    int k = qq2index(q1,q2,Q);
	    Sigma[k1][k2] = r[k] * (sqrt(Sigma[k1][k1] * Sigma[k2][k2]));

	    k2++;
	  }
	}
      }
      
      k1++;
    }
  }

  return;
}



inline double nuGibbs(double *nu,int Q,int G,const int *S,double gamma2,
		      const double *tau2Rho,const double *a,const double *rho,
		      const double *sigma2,const double *phi,
		      const int *psi,const double *x,
		      const int *delta,const double *Delta,Random &ran) {
  double pot = 0.0;
  
  int g;
  for (g = 0; g < G; g++) {
    //
    // compute prior covariance matrix
    //
    
    std::vector<std::vector<double> > var;
    makeSigma(var,Q,gamma2,tau2Rho,a,sigma2 + g * Q,rho);
    
    //
    // define prior mean
    //
    
    std::vector<double> Mean(Q);
    int q;
    for (q = 0; q < Q; q++) Mean[q] = 0.0;
	  
    //
    // compute extra linear and quadratic terms
    //

    std::vector<double> lin;
    std::vector<double> quad;
    lin.resize(Q);
    quad.resize(Q);
    int s;
    for (q = 0; q < Q; q++) {
      int kqg = qg2index(q,g,Q,G);
      double var0 = sigma2[kqg] * phi[kqg];
      double var1 = sigma2[kqg] / phi[kqg];
      for (s = 0; s < S[q]; s++) {
	int kqs = qs2index(q,s,Q,S);
	double variance = psi[kqs] == 0 ? var0 : var1;
	quad[q] += 1.0 / variance;
	int kqgs = qgs2index(q,g,s,Q,G,S);
	lin[q] += (x[kqgs] - delta[kqg] * (2.0 * psi[kqs] - 1.0) * 
		   Delta[kqg]) / variance;
      }
    }
	  
    //
    // Update parameters based on available observations 
    //
    
    std::vector<std::vector<double> > varInv;
    double detPrior = inverse(var,varInv);
	  
    for (q = 0; q < Q; q++) {
      varInv[q][q] += quad[q];
      Mean[q] += lin[q];
    }
	  
    double detPosterior = 1.0 /inverse(varInv,var);
    std::vector<double> mean;
    matrixMult(var,Mean,mean);
	  
    //
    // Draw new values
    //
    
    std::vector<double> vv(Q);
    vv = ran.MultiGaussian(var,mean);
    for (q = 0; q < Q; q++) {
      int k = qg2index(q,g,Q,G);
      nu[k] = vv[q];
    }

    pot += 0.5 * log(detPrior) - 0.5 * log(detPosterior);
    pot += - 0.5 * quadratic(varInv,mean);
  }

  return pot;
}
  


inline double DeltaGibbs(double *Delta,int Q,int G,const int *S,double c2,
			 const double *tau2R,const double *b,const double *r,
			 const double *sigma2,const double *phi,
			 const int *psi,const double *x,
			 const int *delta,const double *nu,Random &ran) {
  double pot = 0.0;
  
  int g;
  for (g = 0; g < G; g++) {
    //
    // compute prior covariance matrix
    //

    int dim = 0;
    vector<int> on(Q,0);
    int q;
    for (q = 0; q < Q; q++) {
      int index = qg2index(q,g,Q,G);
      if (delta[index] == 1) {
	on[q] = 1;
	dim++;
      }
    }

    if (dim > 0) {
      std::vector<std::vector<double> > var;
      makeSigma(var,on,Q,c2,tau2R,b,sigma2 + g * Q,r);
      
      //
      // define prior mean
      //
      
      std::vector<double> Mean(dim);
      int k;
      for (k = 0; k < dim; k++) Mean[k] = 0.0;
      
      //
      // compute extra linear and quadratic terms
      //
      
      vector<double> mean(dim,0.0);
      
      std::vector<double> lin;
      std::vector<double> quad;
      lin.resize(dim);
      quad.resize(dim);
      int s;
      k = 0;
      for (q = 0; q < Q; q++)  if (on[q] == 1) {
	  int index = qg2index(q,g,Q,G);
	  double var0 = sigma2[index] * phi[index];
	  double var1 = sigma2[index] / phi[index];
	  int s;
	  for (s = 0; s < S[q]; s++)
	    {
	      double variance = psi[index] == 0 ? var0 : var1;
	      quad[k] += 1.0 / variance;
	      int xIndex = qgs2index(q,g,s,Q,G,S);
	      lin[k] += (2.0 * psi[index] - 1.0) * 
		(x[xIndex] - nu[index]) / variance;
	    }
	  k++;
	}
      
      //
      // Update parameters based on available observations 
      //
      
      std::vector<std::vector<double> > varInv;
      double detPrior = inverse(var,varInv);
      for (k = 0; k < dim; k++) {
	Mean[k] += lin[k];
	varInv[k][k] += quad[k];
      }
      double detPosterior = 1.0 / inverse(varInv,var);
      matrixMult(var,Mean,mean);
      
      //
      // Draw new values
      //
      
      std::vector<double> vv(dim);
      vv = ran.MultiGaussian(var,mean);
      k = 0;
      for (q = 0; q < Q; q++) {
	if (on[q] == 1) {
	  int index = qg2index(q,g,Q,G);
	  Delta[k] = vv[q];
	  k++;
	}
      }
      
      pot += 0.5 * log(detPrior) - 0.5 * log(detPosterior);
      pot += - 0.5 * quadratic(varInv,mean);
    }
  }

  return pot;
}
  


#endif
