#ifndef POTENTIAL_V2_H
#define POTENTIAL_V2_H

#include <vector>

using namespace std;

#include "Random.h"
#include "Matrix.h"
#include "Utility_v2.h"


inline double potentialXqg(int q,
			   int g,
			   int Q,
			   int G,
			   const int *S,
			   const double *x,
			   const int *psi,
			   const double *nu,
			   const int *delta,
			   const double *Delta,
			   const double *sigma2,
			   const double *phi) {
  double pot = 0.0;

  unsigned int seed = 1;
  Random ran(seed);

  int index = qg2index(q,g,Q,G);
  double var0 = sigma2[index] * phi[index];
  double var1 = sigma2[index] / phi[index];
  double mm = nu[index];

  if (delta[index] != 0) {
    int s;
    for (s = 0; s < psi[q]; s++) {
      double mean;
      double var;
      if (psi[index] == 0) {
	mean = mm - Delta[index];
	var = var0;
      }
      else {
	mean = mm + Delta[index];
	var = var1;
      }

      int xIndex = sqg2index(s,q,g,S,Q,G);
      pot += ran.PotentialGaussian(var,mean,x[xIndex]);
    }
  }
  else {
    int s;
    for (s = 0; s < psi[q]; s++) {
      double var = psi[index] == 0 ? var0 : var1;
      int xIndex = sqg2index(s,q,g,S,Q,G);
      pot += ran.PotentialGaussian(var,mm,x[xIndex]);
    }
  }

  return pot;
}







inline double potentialX(int Q,
			 int G,
			 const int *S,
			 const double *x,
			 const int *psi,
			 const double *nu,
			 const int *delta,
			 const double *Delta,
			 const double *sigma2,
			 const double *phi) {
  double pot = 0.0;
  
  int q,g;
  for (q = 0; q < Q; q++)
    for (g = 0; g < G; g++)
      pot += potentialXqg(q,g,Q,G,S,x,psi,nu,delta,Delta,sigma2,phi);

  return pot;
}




inline double potentialNug(int Q,
			   const double *nug,
			   double gamma2,
			   const double *a,
			   const double *rho,
			   const double *tau2Rho,
			   const double *sigma2g) {
  double pot = 0.0;
  
  unsigned int seed = 1;
  Random ran(seed);
  
  std::vector<std::vector<double> > Sigma;
  makeSigma(Sigma,Q,gamma2,tau2Rho,a,sigma2g,rho);

  std::vector<std::vector<double> > SSigma;
  double determinant = inverse(Sigma,SSigma);

  std::vector<double> value(Q,0.0);
  int q;
  for (q = 0; q < Q; q++)
    value[q] = nug[q];

  pot += ran.PotentialMultiGaussian(SSigma,determinant,value);
  
  return pot;
}





inline double potentialNu(int Q,
			  int G,
			  const double *nu,
			  double gamma2,
			  const double *a,
			  const double *rho,
			  const double *tau2Rho,
			  const double *sigma2) {
  double pot = 0.0;

  int g;
  for (g = 0; g < G; g++)
    pot += potentialNug(Q,nu + g * Q,gamma2,a,rho,tau2Rho,sigma2 + g * Q);
    
  return pot;
}




inline double potentialDDeltag(int Q,
			       const std::vector<int> &on,
			       const double *Deltag,
			       double c2,
			       const double *b,
			       const double *r,
			       const double *tau2R,
			       const double *sigma2g) {
  double pot = 0.0;

  int dim = 0;
  int q;
  for (q = 0; q < Q; q++) dim += on[q];
  
  unsigned int seed = 1;
  Random ran(seed);
  
  std::vector<std::vector<double> > Sigma;
  makeSigma(Sigma,on,Q,c2,tau2R,b,sigma2g,r);

  std::vector<std::vector<double> > SSigma;
  double determinant = inverse(Sigma,SSigma);

  std::vector<double> value(dim,0.0);
  int k = 0;
  for (q = 0; q < Q; q++) {
    if (on[q] == 1) {
      value[k] = Deltag[q];
      k++;
    }
  }

  pot += ran.PotentialMultiGaussian(SSigma,determinant,value);
  
  return pot;
}
			      



inline double potentialDDelta(int Q,
			      int G,
			      const int *delta,
			      const double *Delta,
			      double c2,
			      const double *b,
			      const double *r,
			      const double *tau2R,
			      const double *sigma2) {
  double pot = 0.0;

  int g;
  for (g = 0; g < G; g++) {
    int nOn = 0;
    std::vector<int> on(Q,0);
    int q;
    for (q = 0; q < Q; q++) {
      int index = qg2index(q,g,Q,G);
      on[q] = delta[index];
      nOn += on[q];
    }

    if (nOn > 0)
      pot += potentialDDeltag(Q,on,Delta + g * Q,c2,b,r,tau2R,sigma2 + g * Q);
  }
  
  return pot;
}



inline double potentialA(int Q,
			 const double *a,
			 double pA0,
			 double pA1,
			 double alphaA,
			 double betaA) {
  double pot = 0.0;

  unsigned int seed = 1;
  Random ran(seed);
  
  int q;
  for (q = 0; q < Q; q++) {
    if (a[q] == 0.0)
      pot += - log(pA0);
    else if (a[q] == 1.0)
      pot += - log(pA1);
    else {
      pot += - log(1.0 - pA0 - pA1);
      pot += ran.PotentialBeta(alphaA,betaA,a[q]);
    }
  }
  
  return pot;
}

		     

inline double potentialB(int Q,
			 const double *b,
			 double pB0,
			 double pB1,
			 double alphaB,
			 double betaB) {
  double pot = 0.0;

  unsigned int seed = 1;
  Random ran(seed);
  
  int q;
  for (q = 0; q < Q; q++) {
    if (b[q] == 0.0)
      pot += - log(pB0);
    else if (b[q] == 1.0)
      pot += - log(pB1);
    else {
      pot += - log(1.0 - pB0 - pB1);
      pot += ran.PotentialBeta(alphaB,betaB,b[q]);
    }
  }
  
  return pot;
}


inline double potentialR(int Q,
			 const double *r,
			 double nuR) {
  unsigned int seed = 1;
  Random ran(seed);

  vector<vector<double> > rMatrix;
  rMatrix.resize(Q);
  int qq = 0;
  int p,q;
  for (p = 0; p < Q; p++) {
    rMatrix[p].resize(Q);
    for (q = 0; q < Q; q++) {
      rMatrix[p][q] = r[qq];
      qq++;
    }
  }

  double pot = ran.PotentialCorrelationStandardInverseWishart(nuR,rMatrix);
  
  return pot;
}


inline double potentialRho(int Q,
			   const double *rho,
			   double nuRho) {
  unsigned int seed = 1;
  Random ran(seed);
  
  vector<vector<double> > rhoMatrix;
  rhoMatrix.resize(Q);
  int qq = 0;
  int p,q;
  for (p = 0; p < Q; p++) {
    rhoMatrix[p].resize(Q);
    for (q = 0; q < Q; q++) {
      rhoMatrix[p][q] = rho[qq];
      qq++;
    }
  }
  
  double pot = ran.PotentialCorrelationStandardInverseWishart(nuRho,rhoMatrix);
  
  return pot;
}



inline double potentialSigma2qg(int q,
				int g,
				int Q,
				int G,
				const double *sigma2,
				const double *l,
				const double *t) {
  unsigned int seed = 1;
  Random ran(seed);
  
  double param2 = l[q] / t[q];
  double param1 = l[q] * param2;

  int kqg = qg2index(q,g,Q,G);
  double pot = ran.PotentialGamma(param1,param2,sigma2[kqg]);
  
  return pot;
}



inline double potentialPhiqg(int q,
			     int g,
			     int Q,
			     int G,
			     const double *phi,
			     const double *lambda,
			     const double *theta) {
  unsigned int seed = 1;
  Random ran(seed);
  
  double param2 = lambda[q] / theta[q];
  double param1 = lambda[q] * param2;

  int kqg = qg2index(q,g,Q,G);
  double pot = ran.PotentialGamma(param1,param2,phi[kqg]);
  
  return pot;
}





inline double potentialC2(void) {
  return 0.0;
}


inline double potentialGamma2(void) {
  return 0.0;
}



inline double potentialTau2Rho(void) {
  return 0.0;
}


inline double potentialTau2R(void) {
  return 0.0;
}


inline double potentialT(void) {
  return 0.0;
}


inline double potentialL(void) {
  return 0.0;
}


inline double potentialTheta(void) {
  return 0.0;
}


inline double potentialLambda(void) {
  return 0.0;
}







#endif
