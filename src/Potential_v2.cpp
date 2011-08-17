#include <vector>

#include "Random.h"
#include "Matrix.h"
#include "Potential_v2.h"
#include "Utility_v2.h"


double potentialXqg(int q,
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

  int kqg = qg2index(q,g,Q,G);
  double var0 = sigma2[kqg] * phi[kqg];
  double var1 = sigma2[kqg] / phi[kqg];
  double mm = nu[kqg];

  if (delta[kqg] != 0) {
    int s;
    for (s = 0; s < S[q]; s++) {
      double mean;
      double var;
      
      int ksq = sq2index(s,q,S,Q);
      if (psi[ksq] == 0) {
	mean = mm - Delta[kqg];
	var = var0;
      }
      else {
	mean = mm + Delta[kqg];
	var = var1;
      }

      int xIndex = sqg2index(s,q,g,S,Q,G);
      pot += ran.PotentialGaussian(var,mean,x[xIndex]);
    }
  }
  else {
    int s;
    for (s = 0; s < S[q]; s++) {
      int ksq = sq2index(s,q,S,Q);
      double var = psi[ksq] == 0 ? var0 : var1;
      int xIndex = sqg2index(s,q,g,S,Q,G);
      pot += ran.PotentialGaussian(var,mm,x[xIndex]);
    }
  }

  return pot;
}







double potentialXg(int g, 
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
  
  int q;
  for (q = 0; q < Q; q++)
    pot += potentialXqg(q,g,Q,G,S,x,psi,nu,delta,Delta,sigma2,phi);

  return pot;
}




double potentialX(int Q,
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
  
  int g;
  for (g = 0; g < G; g++)
    pot += potentialXg(g,Q,G,S,x,psi,nu,delta,Delta,sigma2,phi);

  return pot;
}




double potentialNug(int g,
		    int Q,
		    int G,
		    const double *nu,
		    double gamma2,
		    const double *a,
		    const double *rho,
		    const double *tau2Rho,
		    const double *sigma2) {
  double pot = 0.0;
  
  unsigned int seed = 1;
  Random ran(seed);
  
  std::vector<std::vector<double> > Sigma;
  makeSigma(g,G,Sigma,Q,gamma2,tau2Rho,a,sigma2,rho);

  std::vector<std::vector<double> > SSigma;
  double determinant = inverse(Sigma,SSigma);

  std::vector<double> value(Q,0.0);
  int q;
  for (q = 0; q < Q; q++) {
    int kqg = qg2index(q,g,Q,G);
    value[q] = nu[kqg];
  }

  pot += ran.PotentialMultiGaussian(SSigma,determinant,value);
  
  return pot;
}





double potentialNu(int Q,
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
    pot += potentialNug(g,Q,G,nu,gamma2,a,rho,tau2Rho,sigma2);
    
  return pot;
}




double potentialDDeltag(int g,
			int Q,
			int G,
			const std::vector<int> &on,
			const double *Delta,
			double c2,
			const double *b,
			const double *r,
			const double *tau2R,
			const double *sigma2) {
  double pot = 0.0;

  int dim = 0;
  int q;
  for (q = 0; q < Q; q++) dim += on[q];
  
  unsigned int seed = 1;
  Random ran(seed);
  
  std::vector<std::vector<double> > Sigma;
  makeSigma(g,G,Sigma,on,Q,c2,tau2R,b,sigma2,r);

  std::vector<std::vector<double> > SSigma;
  double determinant = inverse(Sigma,SSigma);

  std::vector<double> value(dim,0.0);
  int k = 0;
  for (q = 0; q < Q; q++) {
    if (on[q] == 1) {
      int kqg = qg2index(q,g,Q,G);
      value[k] = Delta[kqg];
      k++;
    }
  }

  pot += ran.PotentialMultiGaussian(SSigma,determinant,value);
  
  return pot;
}
			      



double potentialDDeltag(int g,
			int Q,
			int G,
			const int *delta,
			const double *Delta,
			double c2,
			const double *b,
			const double *r,
			const double *tau2R,
			const double *sigma2) {
  
  std::vector<int> on(Q,0);
  int q;
  for (q = 0; q < Q; q++) {
    int kqg = qg2index(q,g,Q,G);
    on[q] = delta[kqg];
  }


  return potentialDDeltag(g,Q,G,on,Delta,c2,b,r,tau2R,sigma2);
}

double potentialDDelta(int Q,
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
      int kqg = qg2index(q,g,Q,G);
      on[q] = delta[kqg];
      nOn += on[q];
    }

    if (nOn > 0)
      pot += potentialDDeltag(g,Q,G,on,Delta,c2,b,r,tau2R,sigma2);
  }
  
  return pot;
}



double potentialA(int Q,
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

		     

double potentialB(int Q,
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


double potentialR(int Q,
		  const double *r,
		  double nuR) {
  unsigned int seed = 1;
  Random ran(seed);

  vector<vector<double> > rMatrix;
  rMatrix.resize(Q);
  int p,q;
  for (p = 0; p < Q; p++)
    rMatrix[p].resize(Q);

  for (p = 0; p < Q; p++) {
    rMatrix[p][p] = 1.0;
    for (q = p + 1; q < Q; q++) {
      int kqq = qq2index(p,q,Q);
      rMatrix[p][q] = r[kqq];
      rMatrix[q][p] = r[kqq];
    }
  }
  

  double pot = ran.PotentialCorrelationStandardInverseWishart(nuR,rMatrix);
  
  return pot;
}


double potentialRho(int Q,
		    const double *rho,
		    double nuRho) {
  unsigned int seed = 1;
  Random ran(seed);
  
  vector<vector<double> > rhoMatrix;
  rhoMatrix.resize(Q);
  int p,q;
  for (p = 0; p < Q; p++)
    rhoMatrix[p].resize(Q);
  
  for (p = 0; p < Q; p++) {
    rhoMatrix[p][p] = 1.0;
    for (q = p + 1; q < Q; q++) {
      int kqq = qq2index(p,q,Q);
      rhoMatrix[p][q] = rho[kqq];
      rhoMatrix[q][p] = rho[kqq];
    }
  }
  
  
  double pot = ran.PotentialCorrelationStandardInverseWishart(nuRho,rhoMatrix);
  
  return pot;
}



double potentialSigma2qg(int q,
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


double potentialSigma2(int Q,
		       int G,
		       const double *sigma2,
		       const double *l,
		       const double *t) {
  double pot = 0.0;

  int q,g;
  for (q = 0; q < Q; q++)
    for (g = 0; g < G; g++)
      pot += potentialSigma2qg(q,g,Q,G,sigma2,l,t);

  return pot;
}



double potentialPhiqg(int q,
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



double potentialPhi(int Q,
		    int G,
		    const double *phi,
		    const double *lambda,
		    const double *theta) {
  double pot = 0.0;
  
  int q,g;
  for (q = 0; q < Q; q++)
    for (g = 0; g < G; g++)
      pot += potentialPhiqg(q,g,Q,G,phi,lambda,theta);

  return pot;
}



double potentialDeltag(int g,
		       int Q,
		       int G,
		       const int *delta,
		       const double *xi) {
  double pot = 0.0;

  int q;
  for (q = 0; q < Q; q++) {
    int kqg = qg2index(q,g,Q,G);
    if (delta[kqg] == 1)
      pot += - log(xi[q]);
    else
      pot += - log(1.0 - xi[q]);
  }

  return pot;
}


double potentialDeltag_onedelta(int g,
				int Q,
				int G,
				const int *delta,
				const double *xi) {
  double pot = 0.0;
  
  int kqg = qg2index(0,g,Q,G);
  if (delta[kqg] == 1)
    pot += - log(xi[0]);
  else
    pot += - log(1.0 - xi[0]);
  
  return pot;
}


double potentialDelta(int Q,
		      int G,
		      const int *delta,
		      const double *xi) {
  double pot = 0.0;
  
  int g;
  for (g = 0; g < G; g++)
    pot += potentialDeltag(g,Q,G,delta,xi);

  return pot;
}
    



double potentialDelta_onedelta(int Q,
			       int G,
			       const int *delta,
			       const double *xi) {
  double pot = 0.0;
  
  int g;
  for (g = 0; g < G; g++)
    pot += potentialDeltag_onedelta(g,Q,G,delta,xi);

  return pot;
}
    


double potentialXi(int Q,
		   const double *xi,
		   double alphaXi,
		   double betaXi) {
  unsigned int seed = 1;
  Random ran(seed);

  double pot = 0.0;

  int q;
  for (q = 0; q < Q; q++)
    pot += ran.PotentialBeta(alphaXi,betaXi,xi[q]);

  return pot;
}



double potentialXi_onedelta(const double *xi,
			    double alphaXi,
			    double betaXi) {
  unsigned int seed = 1;
  Random ran(seed);

  double pot = ran.PotentialBeta(alphaXi,betaXi,xi[0]);

  return pot;
}




double potentialC2(void) {
  return 0.0;
}


double potentialGamma2(void) {
  return 0.0;
}



double potentialTau2Rho(void) {
  return 0.0;
}


double potentialTau2R(void) {
  return 0.0;
}


double potentialT(void) {
  return 0.0;
}


double potentialL(void) {
  return 0.0;
}


double potentialTheta(void) {
  return 0.0;
}


double potentialLambda(void) {
  return 0.0;
}

