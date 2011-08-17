#include <iostream>
#include <fstream>
#include <stdio.h>

#include "Random.h"
#include "Matrix.h"
#include "Utility_v2.h"
#include "Potential_v2.h"


int qg2index(int q,int g,int Q,int G) {
  int index = g * Q + q;
  
  return index;
}


int sq2index(int s,int q,const int *S,int Q) {
  int index = 0;
  int p;
  for (p = 0; p < q; p++)
    index += S[p];
  
  index += s;
  
  return index;
}


int sqg2index(int s,int q,int g,const int *S,int Q,int G) {
  int sumS = 0;
  int p;
  for (p = 0; p < Q; p++)
    sumS += S[p];
  int index = g * sumS;

  for (p = 0; p < q; p++)
    index += S[p];

  index += s;
  
  return index;
}


int qg2indexNew(int q,int g,int Q,int G) {
  int index = q * G + g;
  
  return index;
}


int sq2indexNew(int s,int q,const int *S,int Q) {
  int index = 0;
  int qq;
  for (qq = 0; qq < q; qq++)
    index += S[qq];
  
  index += s;
  
  return index;
}


int sqg2indexNew(int s,int q,int g,const int *S,int Q,int G) {
  int index = 0;
  int qq;
  for (qq = 0; qq < q; qq++)
    index += G * S[qq];
  
  index += S[q] * g + s;
  
  return index;
}


int qq2index(int q1,int q2,int Q) {
  if (q1 > q2) {
    int temp = q1;
    q1 = q2;
    q2 = temp;
  }

  int index = 0;
  int p1;
  for (p1 = 0; p1 < q1; p1++)
    index += Q - p1 - 1;

  index += q2 - q1 - 1;
  
  return index;
}


void makeSigma(int g,int G,std::vector<std::vector<double> > &Sigma,const int Q,
	       const double gamma2,const double *tau2,
	       const double *a,const double *sigma2,
	       const double *r) {
  Sigma.resize(Q);
  int q;
  for (q = 0; q < Q; q++) {
    int kqg = qg2index(q,g,Q,G);
    Sigma[q].resize(Q);
    Sigma[q][q] = gamma2 * tau2[q];
    Sigma[q][q] *= exp(a[q] * log(sigma2[kqg]));
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



void makeSigma(int g,int G,std::vector<std::vector<double> > &Sigma,
	       const std::vector<int> &on,const int Q,
	       const double gamma2,const double *tau2,
	       const double *a,const double *sigma2,
	       const double *r) {
  int dim = 0;
  int q; 
  for (q = 0; q < Q; q++)
    dim += on[q];

  Sigma.resize(dim);
  int k = 0;
  for (q = 0; q < Q; q++) {
    if (on[q] == 1) {
      int kqg = qg2index(q,g,Q,G);
      Sigma[k].resize(dim);
      Sigma[k][k] = gamma2 * tau2[q];
      Sigma[k][k] *= exp(a[q] * log(sigma2[kqg]));
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
	  if (q2 > q1) {
	    int k = qq2index(q1,q2,Q);
	    Sigma[k1][k2] = r[k] * sqrt(Sigma[k1][k1] * Sigma[k2][k2]);
	    Sigma[k2][k1] = Sigma[k1][k2];
	  }
	  k2++;
	}
      }
      
      k1++;
    }
  }

  return;
}



double nuGibbs(double *nu,int Q,int G,const int *S,double gamma2,
	       const double *tau2Rho,const double *a,const double *rho,
	       const double *sigma2,const double *phi,
	       const int *psi,const double *x,
	       const int *delta,const double *Delta,Random &ran,int draw) {
  double pot = 0.0;
  
  int g;
  for (g = 0; g < G; g++) {
    //
    // compute prior covariance matrix
    //
    
    std::vector<std::vector<double> > var;
    makeSigma(g,G,var,Q,gamma2,tau2Rho,a,sigma2,rho);
    
    //
    // define prior mean
    //
    
    std::vector<double> Mean(Q,0.0);
	  
    //
    // compute extra linear and quadratic terms
    //

    std::vector<double> lin(Q,0.0);
    std::vector<double> quad(Q,0.0);
    int s;
    int q;
    for (q = 0; q < Q; q++) {
      int kqg = qg2index(q,g,Q,G);
      double var0 = sigma2[kqg] * phi[kqg];
      double var1 = sigma2[kqg] / phi[kqg];
      for (s = 0; s < S[q]; s++) {
	int ksq = sq2index(s,q,S,Q);
	double variance = psi[ksq] == 0 ? var0 : var1;
	quad[q] += 1.0 / variance;
	int ksqg = sqg2index(s,q,g,S,Q,G);
	lin[q] += (x[ksqg] - delta[kqg] * (2.0 * psi[ksq] - 1.0) * 
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
    std::vector<double> mean(Q,0.0);
    matrixMult(var,Mean,mean);
	  
    //
    // Draw new values
    //
    
    std::vector<double> vv(Q,0.0);

    if (draw == 1)
      vv = ran.MultiGaussian(var,mean);
    else {
      for (q = 0; q < Q; q++) {
	int kqg = qg2index(q,g,Q,G);
	vv[q] = nu[kqg];
      }
    }

    pot += ran.PotentialMultiGaussian(var,mean,vv);
    
    if (draw == 1) {
      for (q = 0; q < Q; q++) {
	int kqg = qg2index(q,g,Q,G);
	nu[kqg] = vv[q];
      }
    }
  }

  return pot;
}
  


double DeltaGibbs(int g,double *Delta,int Q,int G,const int *S,double c2,
		  const double *tau2R,const double *b,const double *r,
		  const double *sigma2,const double *phi,
		  const int *psi,const double *x,
		  const int *delta,const double *nu,Random &ran,int draw) {
  double pot = 0.0;
  
  //
  // compute prior covariance matrix
  //
  
  int dim = 0;
  std::vector<int> on(Q,0);
  int q;
  for (q = 0; q < Q; q++) {
    int kqg = qg2index(q,g,Q,G);
    if (delta[kqg] == 1) {
      on[q] = 1;
      dim++;
    }
  }
  
  if (dim > 0) {
    std::vector<std::vector<double> > var;
    makeSigma(g,G,var,on,Q,c2,tau2R,b,sigma2,r);
    
    //
    // define prior mean
    //
    
    std::vector<double> Mean(dim,0.0);
    std::vector<double> meanPrior(Mean);
    
    //
    // compute extra linear and quadratic terms
    //
    
    std::vector<double> mean(dim,0.0);
    
    std::vector<double> lin(dim,0.0);
    std::vector<double> quad(dim,0.0);
    int s;
    int k = 0;
    for (q = 0; q < Q; q++) {
      if (on[q] == 1) {
	int kqg = qg2index(q,g,Q,G);
	double var0 = sigma2[kqg] * phi[kqg];
	double var1 = sigma2[kqg] / phi[kqg];
	int s;
	for (s = 0; s < S[q]; s++)
	  {
	    int ksq = sq2index(s,q,S,Q);
	    double variance = psi[ksq] == 0 ? var0 : var1;
	    quad[k] += 1.0 / variance;
	    int xIndex = sqg2index(s,q,g,S,Q,G);
	    lin[k] += (2.0 * psi[ksq] - 1.0) * (x[xIndex] - nu[kqg]) / variance;
	  }
	k++;
      }
    }
    
    //
    // Update parameters based on available observations 
    //
    
    std::vector<std::vector<double> > varInv;
    double detPrior = inverse(var,varInv);
    std::vector<std::vector<double> > covInvPrior(varInv);
    for (k = 0; k < dim; k++) {
      Mean[k] += lin[k];
      varInv[k][k] += quad[k];
    }
    double detPosterior = 1.0 / inverse(varInv,var);
    matrixMult(var,Mean,mean);
    
    //
    // Draw new values
    //
    
    std::vector<double> vv(dim,0.0);
    if (draw == 1) 
      vv = ran.MultiGaussian(var,mean);
    else {
      k = 0;
      for (q = 0; q < Q; q++) {
	if (on[q] == 1) {
	  int kqg = qg2index(q,g,Q,G);
	  vv[k] = Delta[kqg];
	  k++;
	}
      }
    }

    pot += ran.PotentialMultiGaussian(var,mean,vv);

    if (draw == 1) {
      k = 0;
      for (q = 0; q < Q; q++) {
	if (on[q] == 1) {
	  int kqg = qg2index(q,g,Q,G);
	  Delta[kqg] = vv[k];
	  k++;
	}
      }
    }

    /*
    double potTest = 0.0;
    potTest -= ran.PotentialMultiGaussian(var,mean,vv);
    potTest += ran.PotentialMultiGaussian(var,mean,vvOld);

    potTest += potentialDDeltag(g,Q,G,on,Delta,c2,b,r,tau2R,sigma2);
    for (q = 0; q < Q; q++)
      potTest += potentialXqg(q,g,Q,G,S,x,psi,nu,delta,Delta,sigma2,phi);

    for (q = 0; q < Q; q++) {
      int k = qg2index(q,g,Q,G);
      Delta[k] = vvOld[q];
    }

    potTest -= potentialDDeltag(g,Q,G,on,Delta,c2,b,r,tau2R,sigma2);
    for (q = 0; q < Q; q++)
      potTest -= potentialXqg(q,g,Q,G,S,x,psi,nu,delta,Delta,sigma2,phi);

    cout << "nuDDelta: " << potTest << endl;
    */
  }
  
  return pot;
}
  


double DeltaGibbs(double *Delta,int Q,int G,const int *S,double c2,
		  const double *tau2R,const double *b,const double *r,
		  const double *sigma2,const double *phi,
		  const int *psi,const double *x,
		  const int *delta,const double *nu,Random &ran,int draw) {
  double pot = 0.0;
  
  int g;
  for (g = 0; g < G; g++)
    pot += DeltaGibbs(g,Delta,Q,G,S,c2,tau2R,b,r,sigma2,phi,psi,x,
		      delta,nu,ran,draw);
  
  return pot;
}
  
