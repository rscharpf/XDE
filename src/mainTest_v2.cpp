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
#include "Potential_v2.h"
#include "Utility_v2.h"
#include "Update_v2.h"


int main(void) {
  unsigned int seed = 1744334;
  Random ran(seed);

  // initialise fixed hyper-parameters

  int G = 1000;
  int Q = 3;
  int *S = (int *) calloc(Q,sizeof(int));
  S[0] = 100;
  S[1] = 50;
  S[2] = 75;
  int sumS = 0;
  int q;
  for (q = 0; q < Q; q++)
    sumS += S[q];

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

  double c2Max = 10.0;

  // initialise clicinal variables


  int *psi = (int *) calloc(sumS,sizeof(int));
  int s;
  for (q = 0; q < Q; q++)
    for (s = 0; s < S[q]; s++) {
      int ksq = sq2index(s,q,S,Q);
      psi[ksq] = (ran.Unif01() <= 0.5);
    }

  // Initialise parameters to be simulated

  double gamma2 = 0.5 * 0.5;
  double c2 = 0.5 * 0.5;
  double *tau2Rho = (double *) calloc(Q,sizeof(double));
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
  }
  a[0] = 0.0;
  a[1] = 0.5;
  a[2] = 1.0;
  double *b = (double *) calloc(Q,sizeof(double));
  for (q = 0; q < Q; q++) {
    double u = ran.Unif01();
    if (u < pB0)
      b[q] = 0.0;
    else if (u < pB0 + pB1)
      b[q] = 1.0;
    else
      b[q] = ran.Unif01();
  }
  b[0] = 1.0;
  b[1] = 0.0;
  b[2] = 0.5;

  double *l = (double *) calloc(Q,sizeof(double));
  double *t = (double *) calloc(Q,sizeof(double));
  for (q = 0; q < Q; q++) {
    l[q] = 1.0;
    t[q] = 0.5 * 0.5;
  }

  ofstream sigma2file("sigma2.txt");
  double *sigma2 = (double *) calloc(Q * G,sizeof(double));
  int g;
  for (q = 0; q < Q; q++)
    for (g = 0; g < G; g++) {
      int k = qg2index(q,g,Q,G);
      double param2 = l[q] / t[q];
      double param1 = l[q] * param2;
    
      sigma2[k] = ran.Gamma(param1,param2);
      sigma2file << q << " " << g << " " << sigma2[k] << endl;
    }
  sigma2file.close();


  double *lambda = (double *) calloc(Q,sizeof(double));
  double *theta = (double *) calloc(Q,sizeof(double));
  for (q = 0; q < Q; q++) {
    lambda[q] = 0.9;
    theta[q] = 0.1 * 0.1;
  }

  ofstream phifile("phi.txt");
  double *phi = (double *) calloc(Q * G,sizeof(double));
  for (q = 0; q < Q; q++)
    for (g = 0; g < G; g++) {
      int k = qg2index(q,g,Q,G);
      double param2 = lambda[q] / theta[q];
      double param1 = lambda[q] * param2;
      
      phi[k] = ran.Gamma(param1,param2);
      phifile << q << " " << g << " " << phi[k] << " " << 
	sigma2[k] * phi[k] << " " << sigma2[k] / phi[k] << endl;
    }
  phifile.close();

  double *rho = (double *) calloc(Q * (Q - 1) / 2,sizeof(double));
  int q1,q2;
  for (q1 = 0; q1 < Q; q1++) {
    int q2;
    for (q2 = q1 + 1; q2 < Q; q2++) {
      int k = qq2index(q1,q2,Q);
      rho[k] = 0.5;
    }
  }    

  double *r = (double *) calloc(Q * (Q - 1) / 2,sizeof(double));
  for (q1 = 0; q1 < Q; q1++) {
    int q2;
    for (q2 = q1 + 1; q2 < Q; q2++) {
      int k = qq2index(q1,q2,Q);
      r[k] = 0.5;
    }
  }    

  double *nu = (double *) calloc(Q * G,sizeof(double));
  for (g = 0; g < G; g++) {
    std::vector<std::vector<double> > Sigma;
    makeSigma(g,G,Sigma,Q,gamma2,tau2Rho,a,sigma2,rho);
    std::vector<double> zero(Q,0.0);
    std::vector<double> rr(ran.MultiGaussian(Sigma,zero));
    
    int q;
    for (q = 0; q < Q; q++) {
      int kqg = qg2index(q,g,Q,G);
      nu[kqg] = rr[q];
    }
  }
  
  
  int *delta = (int *) calloc(Q * G,sizeof(int));
  int *deltaTrue = (int *) calloc(Q * G,sizeof(int));
  int oneDelta = 1;
  double *xi = (double *) calloc(Q,sizeof(double));
  for (q = 0; q < Q; q++)
    xi[q] = 0.5;
  for (g = 0; g < G; g++) {
    int q;
    int on = (ran.Unif01() <= xi[0]);
    for (q = 0; q < Q; q++) {
      int kqg = qg2index(q,g,Q,G);
      delta[kqg] = on;
      deltaTrue[kqg] = on;
    }
  }


  double *Delta = (double *) calloc(Q * G,sizeof(double));
  for (g = 0; g < G; g++) {
    int nOn = 0;
    vector<int> on(Q,0);
    int q;
    for (q = 0; q < Q; q++) {
      int kqg = qg2index(q,g,Q,G);
      if (delta[kqg] == 1) {
	nOn++;
	on[q] = 1;
      }
    }
    
    if (nOn > 0) {
      std::vector<std::vector<double> > Sigma;
      makeSigma(g,G,Sigma,on,Q,c2,tau2R,b,sigma2,r);
      std::vector<double> zero(nOn,0.0);
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

  ofstream xfile("x.txt");
  double *x = (double *) calloc(G * sumS,sizeof(double));
  for (q = 0; q < Q; q++) {
    int tt = 5;
    for (g = 0; g < G; g++) {
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
	  xfile << q << " " << g << " " << s << " " << x[ksqg] << endl;
	}
      }
      else {
	int s;
	double mean = mm;
	for (s = 0; s < S[q]; s++) {
	  int ksq = sq2index(s,q,S,Q);
	  double var = psi[ksq] == 0 ? var0 : var1;
	  
	  int ksqg = sqg2index(s,q,g,S,Q,G);
	  x[ksqg] = mean + sqrt(var) * ran.Norm01();
	  xfile << q << " " << g << " " << s << " " << x[ksqg] << endl;
	}
      }
    }
    xfile.close();
  }

  // run Metropolis-Hastings updates


  ofstream outt("deltaTest_randomstart.txt");
  for (q = 0; q < Q; q++)
    for (g = 67; g < 67 + 8; g++)
      outt << delta[qg2index(q,g,Q,G)] << " ";
  outt << endl;
  
  // random start
  
  seed = 78434821;

  // Initialise parameters to be simulated
  
  gamma2 = 0.1 * 0.1;
  c2 = 0.1 * 0.1;
  
  
  for (q = 0; q < Q; q++)
    tau2Rho[q] = 1.0;
  int ss;
  for (ss = 0; ss < 25; ss++) {
    int q1 = (int) (ran.Unif01() * Q);
    int q2 = (int) (ran.Unif01() * Q);
    if (q1 != q2) {
      double u = 0.5 + ran.Unif01();
      tau2Rho[q1] *= u;
      tau2Rho[q2] /= u;
    }
  }
  

  for (q = 0; q < Q; q++)
    tau2R[q] = 1.0;
  for (s = 0; s < 25; s++) {
    int q1 = (int) (ran.Unif01() * Q);
    int q2 = (int) (ran.Unif01() * Q);
    if (q1 != q2) {
      double u = 0.5 + ran.Unif01();
      tau2R[q1] *= u;
      tau2R[q2] /= u;
    }
  }


  
  for (q = 0; q < Q; q++) {
    double u = ran.Unif01();
    if (u < pA0)
      a[q] = 0.0;
    else if (u < pA0 + pA1)
      a[q] = 1.0;
    else
      a[q] = ran.Unif01();
  }

  for (q = 0; q < Q; q++) {
    double u = ran.Unif01();
    if (u < pB0)
      b[q] = 0.0;
    else if (u < pB0 + pB1)
      b[q] = 1.0;
    else
      b[q] = ran.Unif01();
  }
  
  
  for (q = 0; q < Q; q++) {
    l[q] = 0.02;
    t[q] = 0.001 * 0.001;
  }
  
  
  for (q = 0; q < Q; q++)
    for (g = 0; g < G; g++) {
      int k = qg2index(q,g,Q,G);
      double param2 = l[q] / t[q];
      double param1 = l[q] * param2;
    
      sigma2[k] = ran.Gamma(param1,param2);
    }
  
  
  for (q = 0; q < Q; q++) {
    lambda[q] = 5.0;
    theta[q] = 1.0 * 1.0;
  }
  
  
  for (q = 0; q < Q; q++)
    for (g = 0; g < G; g++) {
      int k = qg2index(q,g,Q,G);
      double param2 = lambda[q] / theta[q];
      double param1 = lambda[q] * param2;
      
      phi[k] = ran.Gamma(param1,param2);
    }
  
  
  for (q1 = 0; q1 < Q; q1++) {
    int q2;
    for (q2 = q1 + 1; q2 < Q; q2++) {
      int k = qq2index(q1,q2,Q);
      rho[k] = 0.05;
    }
  }    
  
  for (q1 = 0; q1 < Q; q1++) {
    int q2;
    for (q2 = q1 + 1; q2 < Q; q2++) {
      int k = qq2index(q1,q2,Q);
      r[k] = 0.01;
    }
  }    
  
  
  for (g = 0; g < G; g++) {
    std::vector<std::vector<double> > Sigma;
    makeSigma(g,G,Sigma,Q,gamma2,tau2Rho,a,sigma2,rho);
    std::vector<double> zero(Q,0.0);
    std::vector<double> rr(ran.MultiGaussian(Sigma,zero));
    
    int q;
    for (q = 0; q < Q; q++) {
      int k = qg2index(q,g,Q,G);
      nu[k] = rr[q];
    }
  }
  
  
  for (q = 0; q < Q; q++)
    xi[q] = 0.1;
  for (g = 0; g < G; g++) {
    int q;
    int on = (ran.Unif01() <= xi[0]);
    for (q = 0; q < Q; q++) {
      int k = qg2index(q,g,Q,G);
      delta[k] = on;
    }
  }


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
      makeSigma(g,G,Sigma,on,Q,c2,tau2R,b,sigma2,r);
      std::vector<double> zero(Q,0.0);
      std::vector<double> rr(ran.MultiGaussian(Sigma,zero));
      
      int k = 0;
      int q;
      for (q = 0; q < Q; q++) {
	if (on[q] == 1) {
	  int kqg = qg2index(q,g,Q,G);
	  Delta[kqg] = rr[k];
	  k++;	}
      }
    }
  }
  
  
  // end random start
  
  int nCorrect = 0;
  int tt;
  for (tt = 0; tt < Q * G; tt++)
    nCorrect += (delta[tt] == deltaTrue[tt]);
  
    double pot = potentialX(Q,G,S,x,psi,nu,delta,Delta,sigma2,phi) +
      potentialNu(Q,G,nu,gamma2,a,rho,tau2Rho,sigma2) + 
      potentialDDelta(Q,G,delta,Delta,c2,b,r,tau2R,sigma2) + 
      potentialA(Q,a,pA0,pA1,alphaA,betaA) + 
      potentialB(Q,b,pB0,pB1,alphaB,betaB) + 
      potentialR(Q,r,nuR) + 
      potentialRho(Q,rho,nuRho) + 
      potentialSigma2(Q,G,sigma2,l,t) +
      potentialPhi(Q,G,phi,lambda,theta);


    if (oneDelta == 1) 
      pot += potentialDelta_onedelta(Q,G,delta,xi) +
	potentialXi_onedelta(xi,alphaXi,betaXi);
    else
      pot += potentialDelta(Q,G,delta,xi) +
	potentialXi(Q,xi,alphaXi,betaXi);

    cout << endl << "pot: " << pot << " " <<
      potentialX(Q,G,S,x,psi,nu,delta,Delta,sigma2,phi)	<< " " << 
      ((double) nCorrect) / ((double) (Q * G)) << endl << endl;

    cout << "potentialX: " << 
      potentialX(Q,G,S,x,psi,nu,delta,Delta,sigma2,phi) << endl;
    cout << "potentialNu: " << 
      potentialNu(Q,G,nu,gamma2,a,rho,tau2Rho,sigma2) << endl;
    cout << "potentialDDelta: " << 
      potentialDDelta(Q,G,delta,Delta,c2,b,r,tau2R,sigma2) << endl;
    cout << "potentialA: " << potentialA(Q,a,pA0,pA1,alphaA,betaA) << endl;
    cout << "potentialB: " << potentialB(Q,b,pB0,pB1,alphaB,betaB) << endl;
    cout << "potentialR: " << potentialR(Q,r,nuR) << endl;
    cout << "potentialRho: " << potentialRho(Q,rho,nuRho) << endl;
    cout << "potentialSigma2 " << potentialSigma2(Q,G,sigma2,l,t) << endl;
    cout << "potentialPhi: " << potentialPhi(Q,G,phi,lambda,theta) << endl;
    cout << "potentialDelta: ";
    if (oneDelta == 1) 
      cout << potentialDelta_onedelta(Q,G,delta,xi) << endl;
    else
      cout << potentialDelta(Q,G,delta,xi) << endl;
    cout << "potentialXi: ";
    if (oneDelta == 1) 
      cout << potentialXi_onedelta(xi,alphaXi,betaXi) << endl;
    else
      cout << potentialXi(Q,xi,alphaXi,betaXi) << endl;
    cout << endl;



  int nIt = 1000000;
  int k;
  for (k = 0; k < nIt; k++) {
    int nTry = Q;
    int nAccept = 0;
    double epsilonANu = 0.1;

    updateANu(&seed,nTry,&nAccept,epsilonANu,a,nu,Q,G,S,x,psi,delta,Delta,
	      gamma2,rho,sigma2,phi,tau2Rho,pA0,pA1,alphaA,betaA);
    cout << "updateANu: " << nTry << " " << nAccept << endl;
    cout << "a: ";
    for (q = 0; q < Q; q++)
      cout << a[q] << " ";
    cout << endl;

  
    nTry = Q;
    nAccept = 0;
    double epsilonBDDelta = 0.1;
    updateBDDelta(&seed,nTry,&nAccept,epsilonBDDelta,b,Delta,Q,G,S,x,psi,
		  nu,delta,c2,r,sigma2,phi,tau2R,pB0,pB1,alphaB,betaB);
    cout << "updateBDDelta: " << nTry << " " << nAccept << endl;
    cout << "b: ";
    for (q = 0; q < Q; q++)
      cout << b[q] << " ";
    cout << endl;


    nTry = Q;
    nAccept = 0;
    double epsilonTau2RhoNu = 0.02;
    updateTau2RhoNu(&seed,nTry,&nAccept,epsilonTau2RhoNu,tau2Rho,
		    nu,Q,G,S,x,psi,delta,Delta,gamma2,rho,sigma2,
		    phi,a);
    cout << "updateTau2RhoNu: " << nTry << " " << nAccept << endl;
    cout << "tau2Rho: ";
    for (q = 0; q < Q; q++)
      cout << tau2Rho[q] << " ";
    cout << endl;

    
    nTry = Q;
    nAccept = 0;
    double epsilonTau2RDDelta = 0.02;
    updateTau2RDDelta(&seed,nTry,&nAccept,epsilonTau2RDDelta,tau2R,
		      Delta,Q,G,S,x,psi,nu,delta,c2,r,sigma2,
		      phi,b);
    cout << "updateTau2RDDelta: " << nTry << " " << nAccept << endl;
    cout << "tau2R: ";
    for (q = 0; q < Q; q++)
      cout << tau2R[q] << " ";
    cout << endl;


    nTry = Q * G;
    nAccept = 0;
    updateNu(&seed,&nAccept,nu,Q,G,S,x,psi,delta,Delta,gamma2,rho,sigma2,
	     phi,tau2Rho,a);
    cout << "updateNu: " << nTry << " " << nAccept << endl;
    
    
    nTry = Q * G;
    nAccept = 0;
    updateDDelta(&seed,&nAccept,Delta,Q,G,S,x,psi,nu,delta,c2,r,sigma2,
		phi,tau2R,b);
    cout << "updateDDelta: " << nTry << " " << nAccept << endl;


    nTry = 1;
    nAccept = 0;
    updateC2(&seed,nTry,&nAccept,&c2,Q,G,delta,Delta,r,sigma2,tau2R,b,c2Max);
    cout << "updateC2: " << nTry << " " << nAccept << endl;
    cout << "c2: " << c2 << endl;

  
    nTry = 4;
    nAccept = 0;
    double epsilonC2DDelta = 0.2;
    updateC2DDelta(&seed,nTry,&nAccept,epsilonC2DDelta,&c2,Delta,Q,G,S,x,psi,
		   nu,delta,r,sigma2,phi,tau2R,b,c2Max);
    cout << "updateC2DDelta: " << nTry << " " << nAccept << endl;
    cout << "c2: " << c2 << endl;
  
    
    nTry = 1;
    nAccept = 0;
    updateGamma2(&seed,&nAccept,&gamma2,Q,G,nu,rho,sigma2,tau2Rho,a);
    cout << "updateGamma2: " << nTry << " " << nAccept << endl;
    cout << "gamma2: " << gamma2 << endl;


    nTry = 4;
    nAccept = 0;
    double epsilonGamma2Nu = 0.2;
    updateGamma2Nu(&seed,nTry,&nAccept,epsilonGamma2Nu,&gamma2,nu,Q,G,S,x,psi,
		   delta,Delta,rho,sigma2,phi,tau2Rho,a);
    cout << "updateGamma2Nu: " << nTry << " " << nAccept << endl;
    cout << "gamma2: " << gamma2 << endl;

    
    
    nTry = Q * (Q - 1) / 2;
    nAccept = 0;
    double epsilonRC2 = 0.2;
    updateRC2(&seed,nTry,&nAccept,epsilonRC2,r,&c2,Q,G,delta,Delta,sigma2,
	      tau2R,b,nuR,c2Max);
    cout << "updateRC2: " << nTry << " " << nAccept << endl;
    cout << "c2: " << c2 << endl;
    cout << "r: ";
    int q1;
    int q2;
    for (q1 = 0; q1 < Q; q1++)
      for (q2 = q1 + 1; q2 < Q; q2++) {
	int kqq = qq2index(q1,q2,Q);
	cout << r[kqq] << " ";
      }
    cout << endl;
    
    
    nTry = Q * (Q - 1) / 2;
    nAccept = 0;
    double epsilonRDDelta = 0.2;
    updateRDDelta(&seed,nTry,&nAccept,epsilonRDDelta,r,Delta,Q,G,S,x,psi,
		  nu,delta,c2,sigma2,phi,tau2R,b,nuR);
    cout << "updateRDDelta: " << nTry << " " << nAccept << endl;
    cout << "r: ";
    for (q1 = 0; q1 < Q; q1++)
      for (q2 = q1 + 1; q2 < Q; q2++) {
	int kqq = qq2index(q1,q2,Q);
	cout << r[kqq] << " ";
      }
    cout << endl;
    

    nTry = Q * (Q - 1) / 2;
    nAccept = 0;
    double epsilonRhoGamma2 = 0.2;
    updateRhoGamma2(&seed,nTry,&nAccept,epsilonRhoGamma2,rho,&gamma2,
		    Q,G,nu,sigma2,tau2Rho,a,nuRho);
    cout << "updateRhoGamma2: " << nTry << " " << nAccept << endl;
    cout << "gamma2: " << gamma2 << endl;
    cout << "rho: ";
    for (q1 = 0; q1 < Q; q1++)
      for (q2 = q1 + 1; q2 < Q; q2++) {
	int kqq = qq2index(q1,q2,Q);
	cout << rho[kqq] << " ";
      }
    cout << endl;


    nTry = Q * (Q - 1) / 2;
    nAccept = 0;
    double epsilonRhoNu = 0.2;
    updateRhoNu(&seed,nTry,&nAccept,epsilonRhoNu,rho,nu,
		Q,G,S,x,psi,delta,Delta,gamma2,sigma2,phi,tau2Rho,a,nuRho);
    cout << "updateRhoNu: " << nTry << " " << nAccept << endl;
    cout << "rho: ";
    for (q1 = 0; q1 < Q; q1++)
      for (q2 = q1 + 1; q2 < Q; q2++) {
	int kqq = qq2index(q1,q2,Q);
	cout << rho[kqq] << " ";
      }
    cout << endl;

    
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
    

    nTry = 250 * Q;
    nAccept = 0;
    double epsilonTheta = 0.4;
    updateTheta(&seed,nTry,&nAccept,epsilonTheta,theta,Q,G,phi,lambda);
    cout << "updateTheta: " << nTry << " " << nAccept << endl;
    cout << "theta: ";
    for (q = 0; q < Q; q++)
      cout << theta[q] << " ";
    cout << endl;


    
    nTry = 250 * Q;
    nAccept = 0;
    double epsilonLambda = 0.01;
    updateLambda(&seed,nTry,&nAccept,epsilonLambda,lambda,Q,G,phi,theta);
    cout << "updateLambda: " << nTry << " " << nAccept << endl;
    cout << "lambda: ";
    for (q = 0; q < Q; q++)
      cout << lambda[q] << " ";
    cout << endl;


    nTry = 10 * Q;
    nAccept = 0;
    double epsilonLambdaPhi = 0.025;
    updateLambdaPhi(&seed,nTry,&nAccept,epsilonLambdaPhi,lambda,phi,Q,G,S,
		    x,psi,nu,delta,Delta,sigma2,theta);
    cout << "updateLambdaPhi: " << nTry << " " << nAccept << endl;
    cout << "lambda: ";
    for (q = 0; q < Q; q++)
      cout << lambda[q] << " ";
    cout << endl;


    
    nTry = 10 * Q;
    nAccept = 0;
    double epsilonThetaPhi = 0.05;
    updateThetaPhi(&seed,nTry,&nAccept,epsilonThetaPhi,theta,phi,Q,G,S,
		    x,psi,nu,delta,Delta,sigma2,lambda);
    cout << "updateThetaPhi: " << nTry << " " << nAccept << endl;
    cout << "theta: ";
    for (q = 0; q < Q; q++)
      cout << theta[q] << " ";
    cout << endl;
    


   
    nTry = 250 * Q;
    nAccept = 0;
    double epsilonT = 0.4;
    updateT(&seed,nTry,&nAccept,epsilonT,t,Q,G,sigma2,l);
    cout << "updateT: " << nTry << " " << nAccept << endl;
    cout << "t: ";
    for (q = 0; q < Q; q++)
      cout << t[q] << " ";
    cout << endl;
    

    nTry = 250 * Q;
    nAccept = 0;
    double epsilonL = 0.05;
    updateL(&seed,nTry,&nAccept,epsilonL,l,Q,G,sigma2,t);
    cout << "updateL: " << nTry << " " << nAccept << endl;
    cout << "l: ";
    for (q = 0; q < Q; q++)
      cout << l[q] << " ";
    cout << endl;


    
    nTry = 10 * Q;
    nAccept = 0;
    double epsilonLSigma2 = 0.025;
    updateLSigma2(&seed,nTry,&nAccept,epsilonLSigma2,l,sigma2,Q,G,S,x,psi,nu,
		  delta,Delta,c2,gamma2,r,rho,phi,t,tau2R,tau2Rho,a,b);
    cout << "updateLSigma2: " << nTry << " " << nAccept << endl;
    cout << "l: ";
    for (q = 0; q < Q; q++)
      cout << l[q] << " ";
    cout << endl;
        
    
    nTry = 10 * Q;
    nAccept = 0;
    double epsilonTSigma2 = 0.05;
    updateTSigma2(&seed,nTry,&nAccept,epsilonTSigma2,t,sigma2,Q,G,S,x,psi,nu,
		  delta,Delta,c2,gamma2,r,rho,phi,l,tau2R,tau2Rho,a,b);
    cout << "updateTSigma2: " << nTry << " " << nAccept << endl;
    cout << "t: ";
    for (q = 0; q < Q; q++)
      cout << t[q] << " ";
    cout << endl;
    
    
    

    nTry = Q;
    nAccept = 0;
    if (oneDelta == 1)
      updateXi_onedelta(&seed,&nAccept,xi,Q,G,delta,alphaXi,betaXi);
    else
      updateXi(&seed,&nAccept,xi,Q,G,delta,alphaXi,betaXi);
    cout << "updateXi: " << nTry << " " << nAccept << endl;
    cout << "xi: ";
    for (q = 0; q < Q; q++)
      cout << xi[q] << " ";
    cout << endl;


    nTry = G;
    nAccept = 0;
    if (oneDelta == 1)
      updateDeltaDDelta_onedelta(&seed,nTry,&nAccept,delta,Delta,Q,G,S,x,
				 psi,nu,c2,r,sigma2,phi,tau2R,xi,b);
    else
      updateDeltaDDelta(&seed,nTry,&nAccept,delta,Delta,Q,G,S,x,
			psi,nu,c2,r,sigma2,phi,tau2R,xi,b);
    cout << "updateDeltaDDelta: " << nTry << " " << nAccept << endl;
    vector<int> nOn(Q,0);
    for (q = 0; q < Q; q++)
      for (g = 0; g < G; g++) {
	int kqg = qg2index(q,g,Q,G);
	nOn[q] += delta[kqg];
      }
    cout << "nOn: ";
    for (q = 0; q < Q; q++)
      cout << ((double) nOn[q]) / ((double) G) << " ";
    cout << endl;


    int nCorrect = 0;
    int tt;
    for (tt = 0; tt < Q * G; tt++)
      nCorrect += (delta[tt] == deltaTrue[tt]);

    double pot = potentialX(Q,G,S,x,psi,nu,delta,Delta,sigma2,phi) +
      potentialNu(Q,G,nu,gamma2,a,rho,tau2Rho,sigma2) + 
      potentialDDelta(Q,G,delta,Delta,c2,b,r,tau2R,sigma2) + 
      potentialA(Q,a,pA0,pA1,alphaA,betaA) + 
      potentialB(Q,b,pB0,pB1,alphaB,betaB) + 
      potentialR(Q,r,nuR) + 
      potentialRho(Q,rho,nuRho) + 
      potentialSigma2(Q,G,sigma2,l,t) +
      potentialPhi(Q,G,phi,lambda,theta);


    if (oneDelta == 1) 
      pot += potentialDelta_onedelta(Q,G,delta,xi) +
	potentialXi_onedelta(xi,alphaXi,betaXi);
    else
      pot += potentialDelta(Q,G,delta,xi) +
	potentialXi(Q,xi,alphaXi,betaXi);

    cout << endl << "pot: " << pot << " " <<
      potentialX(Q,G,S,x,psi,nu,delta,Delta,sigma2,phi)	<< " " << 
      ((double) nCorrect) / ((double) (Q * G)) << endl << endl;

    cout << "potentialX: " << 
      potentialX(Q,G,S,x,psi,nu,delta,Delta,sigma2,phi) << endl;
    cout << "potentialNu: " << 
      potentialNu(Q,G,nu,gamma2,a,rho,tau2Rho,sigma2) << endl;
    cout << "potentialDDelta: " << 
      potentialDDelta(Q,G,delta,Delta,c2,b,r,tau2R,sigma2) << endl;
    cout << "potentialA: " << potentialA(Q,a,pA0,pA1,alphaA,betaA) << endl;
    cout << "potentialB: " << potentialB(Q,b,pB0,pB1,alphaB,betaB) << endl;
    cout << "potentialR: " << potentialR(Q,r,nuR) << endl;
    cout << "potentialRho: " << potentialRho(Q,rho,nuRho) << endl;
    cout << "potentialSigma2 " << potentialSigma2(Q,G,sigma2,l,t) << endl;
    cout << "potentialPhi: " << potentialPhi(Q,G,phi,lambda,theta) << endl;
    cout << "potentialDelta: ";
    if (oneDelta == 1) 
      cout << potentialDelta_onedelta(Q,G,delta,xi) << endl;
    else
      cout << potentialDelta(Q,G,delta,xi) << endl;
    cout << "potentialXi: ";
    if (oneDelta == 1) 
      cout << potentialXi_onedelta(xi,alphaXi,betaXi) << endl;
    else
      cout << potentialXi(Q,xi,alphaXi,betaXi) << endl;
    cout << endl;
    
    
    for (q = 0; q < Q; q++)
      for (g = 67; g < 67 + 8; g++)
	outt << delta[qg2index(q,g,Q,G)] << " ";
    outt << endl;
    outt.flush();
    
    
    
    
  }
  
  return 0;
}
