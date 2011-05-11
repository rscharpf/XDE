#include <string>
#include <sstream>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <limits.h>

#include "Random.h"
#include "Utility_v2.h"
#include "Potential_v2.h"
#include "Update_v2.h"


int main(void) {
  unsigned int seed = 1784378;
  Random ran(seed);

  // initialise fixed hyper-parameters

  int G = 1000;
  int Q = 4;
  int *S = (int *) calloc(Q,sizeof(int));
  S[0] = 40;
  S[1] = 40;
  S[2] = 40;
  S[3] = 40;
  int sumS = S[0] + S[1] + S[2] + S[3];

  double alphaA = 1.0;
  double betaA = 1.0;
  double pA0 = 0.1;
  double pA1 = 0.1;

  double alphaB = 1.0;
  double betaB = 1.0;
  double pB0 = 0.1;
  double pB1 = 0.1;

  double nuR = 1.0 + Q;
  double nuRho = 1.0 + Q;
  
  double alphaXi = 1.0;
  double betaXi = 1.0;

  double c2Max = 1.0;

  // initialise clicinal variables


  int *psi = (int *) calloc(Q * sumS,sizeof(int));
  int k;
  for (k = 0; k < Q * sumS; k++)
    psi[k] = (ran.Unif01() <= 0.5);

  // Initialise parameters to be simulated

  double gamma2 = 0.1 * 0.1;
  double c2 = 0.1 * 0.1;
  double *tau2Rho = (double *) calloc(Q,sizeof(double));
  int q;
  for (q = 0; q < Q; q++)
    tau2Rho[q] = 0.25;
  double *tau2R = (double *) calloc(Q,sizeof(double));
  for (q = 0; q < Q; q++)
    tau2R[q] = 0.25;
  double *a = (double *) calloc(Q,sizeof(double));
  for (q = 0; q < Q; q++)
    a[q] = 1.0;
  double *b = (double *) calloc(Q,sizeof(double));
  for (q = 0; q < Q; q++)
    b[q] = 0.0;

  double *l = (double *) calloc(Q,sizeof(double));
  double *t = (double *) calloc(Q,sizeof(double));
  for (q = 0; q < Q; q++) {
    l[q] = 0.5 * 0.5;
    t[q] = 0.1 * 0.1;
  }

  double *sigma2 = (double *) calloc(Q * G,sizeof(double));
  int g;
  for (q = 0; q < Q; q++)
    for (g = 0; g < G; g++) {
      int k = qg2index(q,g,Q,G);
      double param2 = l[q] / t[q];
      double param1 = l[q] * param2;
      
      sigma2[k] = ran.Gamma(param1,param2);
    }


  double *lambda = (double *) calloc(Q,sizeof(double));
  double *theta = (double *) calloc(Q,sizeof(double));
  for (q = 0; q < Q; q++) {
    lambda[q] = 0.25 * 0.25;
    theta[q] = 0.05 * 0.05;
  }

  double *phi = (double *) calloc(Q * G,sizeof(double));
  for (q = 0; q < Q; q++)
    for (g = 0; g < G; g++) {
      int k = qg2index(q,g,Q,G);
      double param2 = lambda[q] / theta[q];
      double param1 = lambda[q] * param2;
      
      phi[k] = ran.Gamma(param1,param2);
    }

  double *rho = (double *) calloc(Q*Q,sizeof(double));
  int q1,q2;
  for (q1 = 0; q1 < Q; q1++) {
    int k = qq2index(q1,q1,Q);
    rho[k] = 1.0;
    int q2;
    for (q2 = q1 + 1; q2 < Q; q2++) {
      int k = qq2index(q1,q2,Q);
      rho[k] = 0.5;
      k = qq2index(q2,q1,Q);
      rho[k] = 0.5;
    }
  }    

  double *r = (double *) calloc(Q*Q,sizeof(double));
  for (q1 = 0; q1 < Q; q1++) {
    int k = qq2index(q1,q1,Q);
    r[k] = 1.0;
    int q2;
    for (q2 = q1 + 1; q2 < Q; q2++) {
      int k = qq2index(q1,q2,Q);
      r[k] = 0.5;
      k = qq2index(q2,q1,Q);
      r[k] = 0.5;
    }
  }    

  double *nu = (double *) calloc(Q * G,sizeof(double));
  for (g = 0; g < G; g++) {
    std::vector<std::vector<double> > Sigma;
    makeSigma(Sigma,Q,gamma2,tau2Rho,a,sigma2 + Q * g,rho);
    std::vector<double> zero(Q,0.0);
    std::vector<double> rr(ran.MultiGaussian(Sigma,zero));
    
    int q;
    for (q = 0; q < Q; q++) {
      int k = qg2index(q,g,Q,G);
      nu[k] = rr[q];
    }
  }
  
  
  int *delta = (int *) calloc(Q * G,sizeof(int));
  int oneDelta = 1;
  double xi = 0.2;
  for (g = 0; g < G; g++) {
    int on = (ran.Unif01() <= xi);
    
    int q;
    for (q = 0; q < Q; q++) {
      int k = qg2index(q,g,Q,G);
      delta[k] = on;
    }
  }


  double *Delta = (double *) calloc(Q * G,sizeof(double));
  for (g = 0; g < G; g++) {
    std::vector<std::vector<double> > Sigma;
    makeSigma(Sigma,Q,c2,tau2R,b,sigma2 + Q * g,r);
    std::vector<double> zero(Q,0.0);
    std::vector<double> rr(ran.MultiGaussian(Sigma,zero));
    
    int q;
    for (q = 0; q < Q; q++) {
      int k = qg2index(q,g,Q,G);
      if (delta[k] == 1)
	Delta[k] = rr[q];
      else
	Delta[k] = 0.0;
    }
  }
   
  // run Metropolis-Hastings updates

  int nIt = 10000;
  for (k = 0; k < nIt; k++) {
    int nTry = Q;
    int nAccept = 0;

    double epsilonA = 0.05;
    updateA(&seed,epsilonA,nTry,&nAccept,a,Q,G,pA0,pA1,alphaA,betaA,
	    nu,gamma2,rho,tau2Rho,sigma2);
    cout << nTry << " " << nAccept << endl;
  }
  
  return 0;
}
