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
  unsigned int seed = 17843783;
  Random ran(seed);

  // initialise fixed hyper-parameters

  int G = 1000;
  int Q = 4;
  int *S = (int *) calloc(Q,sizeof(int));
  S[0] = 200;
  S[1] = 200;
  S[2] = 200;
  S[3] = 200;
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


  int *psi = (int *) calloc(sumS,sizeof(int));
  int k;
  for (k = 0; k < sumS; k++)
    psi[k] = (ran.Unif01() <= 0.5);

  // Initialise parameters to be simulated

  double gamma2 = 0.1 * 0.1;
  double c2 = 0.1 * 0.1;
  double *tau2Rho = (double *) calloc(Q,sizeof(double));
  int q;
  for (q = 0; q < Q; q++)
    tau2Rho[q] = 1.0;
  double *tau2R = (double *) calloc(Q,sizeof(double));
  for (q = 0; q < Q; q++)
    tau2R[q] = 1.0;
  double *a = (double *) calloc(Q,sizeof(double));
  for (q = 0; q < Q; q++) {
    double u = ran.Unif01();
    if (u < pA0)
      a[q] = 0.0;
    else if (u < pA0 + pA1)
      a[q] = 1.0;
    else
      a[q] = ran.Unif01();

    a[q] = 0.5;
  }
  double *b = (double *) calloc(Q,sizeof(double));
  for (q = 0; q < Q; q++) {
    double u = ran.Unif01();
    if (u < pB0)
      b[q] = 0.0;
    else if (u < pB0 + pB1)
      b[q] = 1.0;
    else
      b[q] = ran.Unif01();

    b[q] = 0.5;
  }

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
  int oneDelta = 0;
  double *xi = (double *) calloc(Q,sizeof(double));
  for (q = 0; q < Q; q++)
    xi[q] = 0.5;
  for (g = 0; g < G; g++) {
    int q;
    int on = (ran.Unif01() <= xi[0]);
    for (q = 0; q < Q; q++) {
      int k = qg2index(q,g,Q,G);
      delta[k] = on;
    }
  }


  double *Delta = (double *) calloc(Q * G,sizeof(double));
  for (g = 0; g < G; g++) {
    int nOn = 0;
    vector<int> on(Q,0);
    int q;
    for (q = 0; q < Q; q++) {
      int k = qg2index(q,g,Q,G);
      if (delta[k] == 1) {
	nOn++;
	on[q] = 1;
      }
    }
    
    if (nOn > 0) {
      std::vector<std::vector<double> > Sigma;
      makeSigma(Sigma,on,Q,c2,tau2R,b,sigma2 + Q * g,r);
      std::vector<double> zero(Q,0.0);
      std::vector<double> rr(ran.MultiGaussian(Sigma,zero));
      
      int k = 0;
      int q;
      for (q = 0; q < Q; q++) {
	if (on[q] == 1) {
	  int kqg = qg2index(q,g,Q,G);
	  Delta[kqg] = rr[k];
	  k++;
	}
      }
    }
  }


  double *x = (double *) calloc(G * sumS,sizeof(double));
  for (g = 0; g < G; g++)
    for (q = 0; q < Q; q++) {
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

	  int ksqg = sqg2index(s,q,g,S,Q,G);
	  x[ksqg] = mean + sqrt(var) * ran.Norm01();
	}
      }
      else {
	int s;
	double mean = mm;
	for (s = 0; s < S[q]; s++) {
	  int kqs = sq2index(s,q,S,Q);
	  double var = psi[kqs] == 0 ? var0 : var1;
	  
	  int ksqg = sqg2index(s,q,g,S,Q,G);
	  x[ksqg] = mean + sqrt(var) * ran.Norm01();
	}
      }
    }

  // run Metropolis-Hastings updates

  for (q = 0; q < Q; q++) 
    xi[q] = 0.01;

  int nIt = 10000;
  for (k = 0; k < nIt; k++) {
    int nTry = Q * 10;
    int nAccept = 0;
    double epsilonA = 0.05;


    updateA(&seed,nTry,&nAccept,epsilonA,a,Q,G,nu,gamma2,rho,sigma2,tau2Rho,
	    pA0,pA1,alphaA,betaA);
    cout << "updateA: " << nTry << " " << nAccept << endl;
    cout << "a: " << a[0] << " " << a[1] << " " << a[2] << " " << a[3] << endl;


    nTry = Q * 10;
    nAccept = 0;
    double epsilonB = 0.05;
    updateB(&seed,nTry,&nAccept,epsilonB,b,Q,G,delta,Delta,c2,r,sigma2,tau2R,
	    pB0,pB1,alphaB,betaB);
    cout << "updateB: " << nTry << " " << nAccept << endl;
    cout << "b: " << b[0] << " " << b[1] << " " << b[2] << " " << b[3] << endl;
    

    nTry = 1 * 10;
    nAccept = 0;
    double epsilonTau2RhoNu = 0.02;
    updateTau2RhoNu(&seed,nTry,&nAccept,epsilonTau2RhoNu,tau2Rho,
		    nu,Q,G,S,x,psi,delta,Delta,gamma2,rho,sigma2,
		    phi,a);
    cout << "updateTau2RhoNu: " << nTry << " " << nAccept << endl;

    
    nTry = 1 * 10;
    nAccept = 0;
    double epsilonTau2RDDelta = 0.02;
    updateTau2RDDelta(&seed,nTry,&nAccept,epsilonTau2RDDelta,tau2R,
		      Delta,Q,G,S,x,psi,nu,delta,c2,r,sigma2,
		      phi,b);
    cout << "updateTau2RDDelta: " << nTry << " " << nAccept << endl;
    

    nTry = 1;
    nAccept = 0;
    updateNu(&seed,&nAccept,nu,Q,G,S,x,psi,delta,Delta,gamma2,rho,sigma2,
	     phi,tau2Rho,a);
    cout << "updateNu: " << nTry << " " << nAccept << endl;


    nTry = 1;
    nAccept = 0;
    updateDDelta(&seed,&nAccept,Delta,Q,G,S,x,psi,nu,delta,c2,r,sigma2,
		phi,tau2R,b);
    cout << "updateDelta: " << nTry << " " << nAccept << endl;

  
    nTry = 5;
    nAccept = 0;
    updateC2(&seed,nTry,&nAccept,&c2,Q,G,delta,Delta,r,sigma2,tau2R,b,c2Max);
    cout << "updateC2: " << nTry << " " << nAccept << endl;
    cout << "c2: " << c2 << endl;


    nTry = 1;
    nAccept = 0;
    updateGamma2(&seed,&nAccept,&gamma2,Q,G,nu,rho,sigma2,tau2Rho,a);
    cout << "updateC2: " << nTry << " " << nAccept << endl;
    cout << "gamma2: " << gamma2 << endl;


    nTry = 5;
    nAccept = 0;
    double epsilonRC2 = 0.01;
    updateRC2(&seed,nTry,&nAccept,epsilonRC2,r,&c2,Q,G,delta,Delta,sigma2,
	      tau2R,b,nuR,c2Max);
    cout << "updateRC2: " << nTry << " " << nAccept << endl;
    cout << "c2: " << c2 << endl;


    nTry = 5;
    nAccept = 0;
    double epsilonRhoGamma2 = 0.01;
    updateRhoGamma2(&seed,nTry,&nAccept,epsilonRhoGamma2,rho,&gamma2,
		    Q,G,nu,sigma2,tau2Rho,a,nuRho);
    cout << "updateRhoGamma2: " << nTry << " " << nAccept << endl;
    cout << "gamma2: " << gamma2 << endl;


    nTry = 5 * Q * G;
    nAccept = 0;
    double epsilonSigma2 = 0.5;
    updateSigma2(&seed,nTry,&nAccept,epsilonSigma2,sigma2,Q,G,S,x,psi,nu,
		 delta,Delta,c2,gamma2,r,rho,phi,t,l,tau2R,tau2Rho,a,b);
    cout << "updateSigma2: " << nTry << " " << nAccept << endl;


    nTry = 5 * Q * G;
    nAccept = 0;
    double epsilonPhi = 0.5;
    updatePhi(&seed,nTry,&nAccept,epsilonPhi,phi,Q,G,S,x,psi,nu,
	      delta,Delta,sigma2,theta,lambda);
    cout << "updatePhi: " << nTry << " " << nAccept << endl;


    nTry = 5 * Q;
    nAccept = 0;
    double epsilonTheta = 0.1;
    updateTheta(&seed,nTry,&nAccept,epsilonTheta,theta,Q,G,phi,lambda);
    cout << "updateTheta: " << nTry << " " << nAccept << endl;


    nTry = 5 * Q;
    nAccept = 0;
    double epsilonLambda = 0.1;
    updateLambda(&seed,nTry,&nAccept,epsilonLambda,lambda,Q,G,phi,theta);
    cout << "updateLambda: " << nTry << " " << nAccept << endl;


    nTry = 5 * Q;
    nAccept = 0;
    double epsilonT = 0.1;
    updateT(&seed,nTry,&nAccept,epsilonT,t,Q,G,sigma2,l);
    cout << "updateT: " << nTry << " " << nAccept << endl;


    nTry = 5 * Q;
    nAccept = 0;
    double epsilonL = 0.1;
    updateL(&seed,nTry,&nAccept,epsilonT,l,Q,G,sigma2,t);
    cout << "updateL: " << nTry << " " << nAccept << endl;


    nTry = Q;
    nAccept = 0;
    updateXi(&seed,&nAccept,xi,Q,G,delta,alphaXi,betaXi);
    cout << "updateXi: " << nTry << " " << nAccept << endl;
    cout << "xi: " << xi[0] << " " << xi[1] << " " << 
      xi[2] << " " << xi[3] << endl;


  }
  
  return 0;
}
