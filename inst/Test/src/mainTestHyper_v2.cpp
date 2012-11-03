#include <string>
#include <sstream>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <limits.h>

#include "Random_v2.h"
#include "Potential_v2.h"
#include "Utility_v2.h"
#include "Update_v2.h"

extern "C" {

  // int main(void) {
  void initializeParams(int *nIt, int *seedR, int *G, int *Q, int *S,
			int *psi,
			double *alphaA,
			double *betaA,
			double *pA0,
			double *pA1,
			double *alphaB,
			double *betaB,
			double *pB0,
			double *pB1,
			double *nuR,
			double *nuRho,
			double *alphaXi,
			double *betaXi,
			double *c2Max,
			double *sigma2,
			double *tau2Rho,
			double *gamma2,
			double *tau2R,
			int *simulateExpression, // whether to simulate data from the prior
			//double *x, // expression data
			int *aOut,
			int *sigma2Out,
			int *simulateSigma2,
			int *oldCliqueInput,
			int *oldComponentsInput,
			int *nClique,
			int *nOldClique,
			int *nTotalInCliqueInput){//22
    unsigned int seed;
    seed = *seedR;
    Random ran(seed);
    int sumS = 0;
    int q;
    for (q = 0; q < *Q; q++)
      sumS += S[q];

    vector<int> oldClique;
    vector<vector<int> > oldComponents;
    vector<int> clique(*G,0);
    vector<int> nNewInClique(0,0);

    // goal: to pass cliques from R to C, we need to transform the input parameters
    //       oldCliqueInput and oldComponentsInput
    //       This transformation requires variables nClique, oldCliqueInput, nOldClique
    // what is nClique?
    // transformGraph(nClique, oldCliqueInput, nOldClique, oldComponentsInput, oldClique, oldComponents);

    // This section defines oldClique, oldComponents, and nTotalInClique.
    // The remaining vars are temporary
    clique[0] = 0;
    int g;
    for (g = 1; g < *G; g++) {
      if (ran.Unif01() < 0.95)
	clique[g] = clique[g - 1];
      else
	clique[g] = clique[g - 1] + 1;
    }
    nNewInClique.resize(clique[*G - 1] + 1);
    int i;
    for (i = 0; i < nNewInClique.size(); i++)
      nNewInClique[i] = 0;
    for (g = 0; g < *G; g++)
      nNewInClique[clique[g]] += 1;

    oldClique.resize(nNewInClique.size());
    oldClique[0] = -1;
    oldComponents.resize(nNewInClique.size());
    oldComponents[0].resize(0);
    vector<int> nTotalInClique(nNewInClique);
    int c;
    for (c = 1; c < nNewInClique.size(); c++) {
      oldClique[c] = (int) (ran.Unif01() * ((double) (c - 1)));
      double u = ran.Unif01();
      int nOld = (int) (u * u * ((double) (nTotalInClique[oldClique[c]] - 1)));
      nTotalInClique[c] += nOld;
      oldComponents[c].resize(nOld);
      int nSampled = 0;
      while (nSampled < nOld) {
	oldComponents[c][nSampled] = (int) (ran.Unif01() * nTotalInClique[oldClique[c]]);
	int existAlready = 0;
	int i;
	for (i = 0; i < nSampled; i++) {
	  if (oldComponents[c][i] == oldComponents[c][nSampled])
	    existAlready = 1;
	}
	if (existAlready == 0)
	  nSampled++;
      }
    }
    for (c = 0; c < nTotalInClique.size(); c++)
      cout << "nTotalInClique: " << c << " " << nTotalInClique[c] << endl;
    // initialise clinical variables

    int s;

    double *a = (double *) calloc(*Q,sizeof(double));
    for (q = 0; q < *Q; q++) {
      double u = ran.Unif01();
      if (u < *pA0)
	a[q] = 0.0;
      else if (u < *pA0 + *pA1)
	a[q] = 1.0;
      else
	a[q] = ran.Unif01();
    }
    a[0] = 0.0;
    a[1] = 0.5;
    a[2] = 1.0;
    double *b = (double *) calloc(*Q,sizeof(double));
    for (q = 0; q < *Q; q++) {
      double u = ran.Unif01();
      if (u < *pB0)
	b[q] = 0.0;
      else if (u < *pB0 + *pB1)
	b[q] = 1.0;
      else
	b[q] = ran.Unif01();
    }
    b[0] = 1.0;
    b[1] = 0.0;
    b[2] = 0.5;

    double *l = (double *) calloc(*Q,sizeof(double));
    double *t = (double *) calloc(*Q,sizeof(double));
    for (q = 0; q < *Q; q++) {
      l[q] = 1.0;
      t[q] = 0.5 * 0.5;
    }

    // RS : simulating sigma2 from the prior.
    //      - perhaps add an option for simulating sigma2 from prior
    //         or using the value of sigma2 passed in from R
    if (*simulateSigma2 == 1){
      for (q = 0; q < *Q; q++){
	for (g = 0; g < *G; g++) {
	  int k = qg2index(q,g,*Q,*G);
	  double param2 = l[q] / t[q];
	  double param1 = l[q] * param2;
	}
      }
    }

    double *lambda = (double *) calloc(*Q,sizeof(double));
    double *theta = (double *) calloc(*Q,sizeof(double));
    for (q = 0; q < *Q; q++) {
      lambda[q] = 0.9;
      theta[q] = 0.1 * 0.1;
    }

    ofstream phifile("phi.txt");
    double *phi = (double *) calloc(*Q * *G,sizeof(double));
    for (q = 0; q < *Q; q++)
      for (g = 0; g < *G; g++) {
	int k = qg2index(q,g,*Q,*G);
	double param2 = lambda[q] / theta[q];
	double param1 = lambda[q] * param2;
	phi[k] = ran.Gamma(param1,param2);

	phifile << q << " " << g << " " << phi[k] << " " << sigma2[k] * phi[k] << " " << sigma2[k] / phi[k] << endl;
      }
    phifile.close();

    double *rho = (double *) calloc(*Q * (*Q - 1) / 2,sizeof(double));
    int q1,q2;
    for (q1 = 0; q1 < *Q; q1++) {
      int q2;
      for (q2 = q1 + 1; q2 < *Q; q2++) {
	int k = qq2index(q1,q2,*Q);
	rho[k] = 0.5;
      }
    }

    double *r = (double *) calloc(*Q * (*Q - 1) / 2,sizeof(double));
    for (q1 = 0; q1 < *Q; q1++) {
      int q2;
      for (q2 = q1 + 1; q2 < *Q; q2++) {
	int k = qq2index(q1,q2,*Q);
	r[k] = 0.5;
      }
    }

    double *nu = (double *) calloc(*Q * *G,sizeof(double));
    for (g = 0; g < *G; g++) {
      std::vector<std::vector<double> > Sigma;
      makeSigma(g,*G,Sigma,*Q,*gamma2,tau2Rho,a,sigma2,rho);
      std::vector<double> zero(*Q,0.0);
      std::vector<double> rr(ran.MultiGaussian(Sigma,zero));

      int q;
      for (q = 0; q < *Q; q++) {
	int kqg = qg2index(q,g,*Q,*G);
	nu[kqg] = rr[q];
      }
    }


    int *delta = (int *) calloc(*Q * *G,sizeof(int));

    // simulate delta from an MRF

    vector<vector<int> > neighbour;
    neighbour.resize(*G);

    ifstream in("R.txt");
    for (g = 0; g < *G; g++) {
      neighbour[g].resize(0);
      int gg;
      for (gg = 0; gg < *G; gg++) {
	int vv;
	in >> vv;
	if (vv == 1 && g != gg) {
	  neighbour[g].push_back(gg);
	}
      }
    }

    vector<double> potOn(*G,0.0);
    vector<double> potOff(*G,0.0);

    double alpha = -0.2; alpha = 0.0;
    double beta = 3.0; beta = 0.42;
    double betag = 0.25; betag = 2.0;

    double *xi = (double *) calloc(*Q,sizeof(double));
    for (q = 0; q < *Q; q++) {
      xi[q] = 0.5;
    }

    int oneDelta = 1; oneDelta = 1;
    int *dd = NULL;
    if (oneDelta == 1) {
      dd = (int *) calloc(*G,sizeof(int));
      for (g = 0; g < *G; g++)
	dd[g] = (ran.Unif01() <= xi[0]);
    }
    else {
      dd = (int *) calloc(*Q * *G,sizeof(int));
      for (g = 0; g < *Q * *G; g++)
	dd[g] = (ran.Unif01() <= xi[0]);
    }

    int *deltaTrue = (int *) calloc(*Q * *G,sizeof(int));
    for (g = 0; g < *G; g++) {
      int q;
      for (q = 0; q < *Q; q++) {
	int kqg = qg2index(q,g,*Q,*G);
	if (oneDelta == 1) {
	  delta[kqg] = dd[g];
	  deltaTrue[kqg] = dd[g];
	}
	else {
	  delta[kqg] = dd[kqg];
	  deltaTrue[kqg] = dd[kqg];
	}
      }
    }

    vector<vector<vector<double> > > D;
    D.resize(oldComponents.size());
    for (c = 0; c < D.size(); c++) {
      D[c].resize(nTotalInClique[c]);
      int g1;
      for (g1 = 0; g1 < D[c].size(); g1++) {
	D[c][g1].resize(nTotalInClique[c]);
	int g2;
	for (g2 = 0; g2 < D[c][g1].size(); g2++)
	  D[c][g1][g2] = 2.0 * (g1 == g2);
      }
    }

    double df = 1.0;

    vector<vector<vector<double> > > Omega(ran.HyperInverseWishart(df,D,oldClique,oldComponents));

    vector<vector<double> > zero;
    zero.resize(*G);
    for (g = 0; g < *G; g++) {
      zero[g].resize(*Q);
      for (q = 0; q < *Q; q++)
	zero[g][q] = 0.0;
    }
    vector<vector<double> > R;
    R.resize(*Q);
    for (q = 0; q < *Q; q++) {
      R[q].resize(*Q);
    }
    for (q = 0; q < *Q; q++) {
      R[q][q] = tau2R[q];
      int p;
      for (p = q + 1; p < *Q; p++) {
	R[q][p] = sqrt(tau2R[p] * tau2R[q]) * r[qq2index(p,q,*Q)];
	R[p][q] = R[q][p];
      }
    }

    vector<vector<double> > DeltaStar(ran.MatrixVariateNormal(zero,R,Omega,oldClique,oldComponents));

    double *Delta = (double *) calloc(*Q * *G,sizeof(double));
    for (g = 0; g < *G; g++) {
      for (q = 0; q < *Q; q++)
	Delta[qg2index(q,g,*Q,*G)] = DeltaStar[g][q] * exp(0.5 * b[q] * log(sigma2[qg2index(q,g,*Q,*G)]));
    }

    // I believe this is simulating the expression data
    //if(*simulateExpression == 1){
      ofstream xfile("x.txt");
      double *x = (double *) calloc(*G * sumS,sizeof(double));
      for (q = 0; q < *Q; q++) {
	int tt = 5;
	for (g = 0; g < *G; g++) {
	  int kqg = qg2index(q,g,*Q,*G);

	  double var0 = sigma2[kqg] * phi[kqg];
	  double var1 = sigma2[kqg] / phi[kqg];
	  double mm = nu[kqg];

	  if (delta[kqg] != 0) {
	    int s;
	    for (s = 0; s < S[q]; s++) {
	      double mean;
	      double var;
	      int ksq = sq2index(s,q,S,*Q);
	      if (psi[ksq] == 0) {
		mean = mm - Delta[kqg];
		var = var0;
	      }
	      else {
		mean = mm + Delta[kqg];
		var = var1;
	      }

	      int ksqg = sqg2index(s,q,g,S,*Q,*G);
	      x[ksqg] = mean + sqrt(var) * ran.Norm01();
	      xfile << q << " " << g << " " << s << " " << x[ksqg] << endl;
	    }
	  }
	  else {
	    int s;
	    double mean = mm;
	    for (s = 0; s < S[q]; s++) {
	      int ksq = sq2index(s,q,S,*Q);
	      double var = psi[ksq] == 0 ? var0 : var1;

	      int ksqg = sqg2index(s,q,g,S,*Q,*G);
	      x[ksqg] = mean + sqrt(var) * ran.Norm01();
	      xfile << q << " " << g << " " << s << " " << x[ksqg] << endl;
	    }
	  }
	}
	xfile.close();
      }
      //}

    // random start
    // seed = 784348215;
    // Initialise parameters to be simulated
    //gamma2 = 0.1 * 0.1;
    for (q = 0; q < *Q; q++)
      tau2Rho[q] = 1.0;
    int ss;
    for (ss = 0; ss < 25; ss++) {
      int q1 = (int) (ran.Unif01() * *Q);
      int q2 = (int) (ran.Unif01() * *Q);
      if (q1 != q2) {
	double u = 0.5 + ran.Unif01();
	tau2Rho[q1] *= u;
	tau2Rho[q2] /= u;
      }
    }


    for (q = 0; q < *Q; q++)
      tau2R[q] = 1.0;
    for (s = 0; s < 25; s++) {
      int q1 = (int) (ran.Unif01() * *Q);
      int q2 = (int) (ran.Unif01() * *Q);
      if (q1 != q2) {
	double u = 0.5 + ran.Unif01();
	tau2R[q1] *= u;
	tau2R[q2] /= u;
      }
    }



  for (q = 0; q < *Q; q++) {
    double u = ran.Unif01();
    if (u < *pA0)
      a[q] = 0.0;
    else if (u < *pA0 + *pA1)
      a[q] = 1.0;
    else
      a[q] = ran.Unif01();
  }

  for (q = 0; q < *Q; q++) {
    double u = ran.Unif01();
    if (u < *pB0)
      b[q] = 0.0;
    else if (u < *pB0 + *pB1)
      b[q] = 1.0;
    else
      b[q] = ran.Unif01();
  }


  for (q = 0; q < *Q; q++) {
    l[q] = 0.02;
    t[q] = 0.001 * 0.001;
  }

  if(*simulateSigma2 == 1){
    for (q = 0; q < *Q; q++)
      for (g = 0; g < *G; g++) {
	int k = qg2index(q,g,*Q,*G);
	double param2 = l[q] / t[q];
	double param1 = l[q] * param2;
	sigma2[k] = ran.Gamma(param1,param2);
      }
  }
  // else assume supplied sigma2 is a valid variance

  for (q = 0; q < *Q; q++) {
    lambda[q] = 5.0;
    theta[q] = 1.0 * 1.0;
  }


  for (q = 0; q < *Q; q++)
    for (g = 0; g < *G; g++) {
      int k = qg2index(q,g,*Q,*G);
      double param2 = lambda[q] / theta[q];
      double param1 = lambda[q] * param2;

      phi[k] = ran.Gamma(param1,param2);
    }


  for (q1 = 0; q1 < *Q; q1++) {
    int q2;
    for (q2 = q1 + 1; q2 < *Q; q2++) {
      int k = qq2index(q1,q2,*Q);
      rho[k] = 0.05;
    }
  }

  for (q1 = 0; q1 < *Q; q1++) {
    int q2;
    for (q2 = q1 + 1; q2 < *Q; q2++) {
      int k = qq2index(q1,q2,*Q);
      r[k] = 0.01;
    }
  }


  for (g = 0; g < *G; g++) {
    std::vector<std::vector<double> > Sigma;
    makeSigma(g,*G,Sigma,*Q,*gamma2,tau2Rho,a,sigma2,rho);
    std::vector<double> zero(*Q,0.0);
    std::vector<double> rr(ran.MultiGaussian(Sigma,zero));

    int q;
    for (q = 0; q < *Q; q++) {
      int k = qg2index(q,g,*Q,*G);
      nu[k] = rr[q];
    }
  }

  for (q = 0; q < *Q; q++)
    xi[q] = 0.25;
  for (g = 0; g < *G; g++) {
    int q;
    int on = (ran.Unif01() <= xi[0]);
    for (q = 0; q < *Q; q++) {
      int k = qg2index(q,g,*Q,*G);
      delta[k] = on;
    }
  }

  alpha = 0.0;
  beta = 0.01;
  betag = 0.01;

  // end random start
  ofstream aFile("a.txt");
  ofstream sigma2File("sigma2.txt");
  // run Metropolis-Hastings updates


  int nCorrect = 0;
  int tt;
  for (tt = 0; tt < *Q * *G; tt++)
    nCorrect += (delta[tt] == deltaTrue[tt]);
  cout << "fraction correct: " << ((double) nCorrect) / ((double) (*Q * *G)) << endl;

  double pot = potentialX(*Q,*G,S,x,psi,nu,delta,Delta,sigma2,phi) +
    potentialNu(*Q,*G,nu,*gamma2,a,rho,tau2Rho,sigma2) +
    potentialDDeltaStar_HyperInverseWishart(Delta,b,sigma2,tau2R,r,*Q,*G,Omega,oldClique,oldComponents) +
    potentialOmega_HyperInverseWishart(Omega,D,df,oldClique,oldComponents) +
    potentialA(*Q,a,*pA0,*pA1,*alphaA,*betaA) +
    potentialB(*Q,b,*pB0,*pB1,*alphaB,*betaB) +
    potentialR(*Q,r,*nuR) +
    potentialRho(*Q,rho,*nuRho) +
    potentialSigma2(*Q,*G,sigma2,l,t) +
    potentialPhi(*Q,*G,phi,lambda,theta) +
    potentialDelta(*Q,*G,delta,xi) +
    potentialXi(*Q,xi,*alphaXi,*betaXi);


  //int nIt = 1000000;
    int k;
    for (k = 0; k < *nIt; k++) {

      if (k == 0) {
	cout << endl << "pot: " << pot << " " <<
	  potentialX(*Q,*G,S,x,psi,nu,delta,Delta,sigma2,phi)	<< " " <<
	  ((double) nCorrect) / ((double) (*Q * *G)) << endl << endl;

	cout << "potentialX: " <<
	  potentialX(*Q,*G,S,x,psi,nu,delta,Delta,sigma2,phi) << endl;
	cout << "potentialNu: " <<
	  potentialNu(*Q,*G,nu,*gamma2,a,rho,tau2Rho,sigma2) << endl;
	cout << "potentialDDeltaStar: " <<
	  potentialDDeltaStar_HyperInverseWishart(Delta,b,sigma2,tau2R,r,*Q,*G,Omega,oldClique,oldComponents) << endl;
	cout << "potentialOmega: " <<
	  potentialOmega_HyperInverseWishart(Omega,D,df,oldClique,oldComponents) << endl;
	cout << "potentialA: " << potentialA(*Q,a,*pA0,*pA1,*alphaA,*betaA) << endl;
	cout << "potentialB: " << potentialB(*Q,b,*pB0,*pB1,*alphaB,*betaB) << endl;
	cout << "potentialR: " << potentialR(*Q,r,*nuR) << endl;
	cout << "potentialRho: " << potentialRho(*Q,rho,*nuRho) << endl;
	cout << "potentialSigma2 " << potentialSigma2(*Q,*G,sigma2,l,t) << endl;
	cout << "potentialPhi: " << potentialPhi(*Q,*G,phi,lambda,theta) << endl;
	cout << "potentialDelta: " << potentialDelta(*Q,*G,delta,xi) << endl;
	cout << "potentialXi: " << potentialXi(*Q,xi,*alphaXi,*betaXi) << endl;
	cout << endl;
      }

      int nTry = *Q;
      int nAccept = 0;
      double epsilonANu = 0.01;

      // We have R updates for many of these

      updateANu(&seed,nTry,&nAccept,epsilonANu,a,nu,*Q,*G,S,x,psi,delta,Delta,
		*gamma2,rho,sigma2,phi,tau2Rho,*pA0,*pA1,*alphaA,*betaA);
      cout << "updateANu: " << nTry << " " << nAccept << endl;
      cout << "a: ";
      if(*aOut == 1){
	for (q = 0; q < *Q; q++)
	  aFile << a[q] << " ";
	aFile << endl;
	aFile.flush();
	aFile.close();
      }
      else{
	for (q = 0; q < *Q; q++)
	  cout << a[q] << " ";
	cout << endl;
      }


      if (oneDelta == 1) {
	nTry = *G;
	nAccept = 0;
	updateDelta_HyperInverseWishart_onedelta(&seed,nTry,&nAccept,delta,*Q,*G,S,x,psi,
						 nu,Delta,r,sigma2,phi,xi,b);
	cout << "updateDelta_onedelta: " << nTry << " " << nAccept << endl;


	nTry = *Q;
	nAccept = 0;
	updateXi_onedelta(&seed,&nAccept,xi,*Q,*G,delta,*alphaXi,*betaXi);
	cout << "updateXi_onedelta: " << nTry << " " << nAccept << endl;
	cout << "xi: ";
	for (q = 0; q < *Q; q++)
	  cout << xi[q] << " ";
	cout << endl;

      }
      else {
	nTry = *Q * *G;
	nAccept = 0;
	updateDelta_HyperInverseWishart(&seed,nTry,&nAccept,delta,*Q,*G,S,x,psi,
					nu,Delta,r,sigma2,phi,xi,b);
	cout << "updateDelta: " << nTry << " " << nAccept << endl;


	nTry = *Q;
	nAccept = 0;
	updateXi(&seed,&nAccept,xi,*Q,*G,delta,*alphaXi,*betaXi);
	cout << "updateXi: " << nTry << " " << nAccept << endl;
	cout << "xi: ";
	for (q = 0; q < *Q; q++)
	  cout << xi[q] << " ";
	cout << endl;
      }



      nTry = *Q;
      nAccept = 0;
      double epsilonBDDelta = 0.1;
      updateBDDeltaStar_HyperInverseWishart(&seed,nTry,&nAccept,epsilonBDDelta,b,Delta,*Q,*G,S,x,psi,
					    nu,delta,r,sigma2,phi,tau2R,*pB0,*pB1,*alphaB,*betaB,Omega,
					    oldClique,oldComponents);
      cout << "updateBDDelta_HyperInverseWishart: " << nTry << " " << nAccept << endl;
      cout << "b: ";
      for (q = 0; q < *Q; q++)
	cout << b[q] << " ";
      cout << endl;


      nTry = *Q * (*Q - 1) / 2;
      nAccept = 0;
      double epsilonRDDelta = 0.1;
      updateRDDeltaStar_HyperInverseWishart(&seed,nTry,&nAccept,epsilonRDDelta,r,Delta,*Q,*G,S,x,psi,
					    nu,delta,sigma2,phi,tau2R,b,*nuR,Omega,
					    oldClique,oldComponents);
      cout << "updateRDDelta_HyperInverseWishart: " << nTry << " " << nAccept << endl;
      cout << "r: ";
      for (q = 0; q < *Q; q++)
	cout << r[q] << " ";
      cout << endl;



      nTry = *Q;
      nAccept = 0;
      double epsilonTau2RhoNu = 0.02;
      updateTau2RhoNu(&seed,nTry,&nAccept,epsilonTau2RhoNu,tau2Rho,
		      nu,*Q,*G,S,x,psi,delta,Delta,*gamma2,rho,sigma2,
		      phi,a);
      cout << "updateTau2RhoNu: " << nTry << " " << nAccept << endl;
      cout << "tau2Rho: ";
      for (q = 0; q < *Q; q++)
	cout << tau2Rho[q] << " ";
      cout << endl;


      nTry = *Q;
      nAccept = 0;
      double epsilonTau2RDDelta = 0.02;
      updateTau2RDDeltaStar_HyperInverseWishart(&seed,nTry,&nAccept,epsilonTau2RDDelta,tau2R,
						Delta,*Q,*G,S,x,psi,nu,delta,r,sigma2,phi,b,Omega,
						oldClique,oldComponents);
      cout << "updateTau2RDDelta_HyperInverseWishart: " << nTry << " " << nAccept << endl;
      cout << "tau2R: ";
      for (q = 0; q < *Q; q++)
	cout << tau2R[q] << " ";
      cout << endl;


      nTry = *Q * *G;
      nAccept = 0;
      updateNu(&seed,&nAccept,nu,*Q,*G,S,x,psi,delta,Delta,*gamma2,rho,sigma2,
	       phi,tau2Rho,a);
      cout << "updateNu: " << nTry << " " << nAccept << endl;


      nTry = 1;
      nAccept = 0;
      updateOmega_HyperInverseWishart(&seed,&nAccept,Omega,*Q,*G,Delta,r,sigma2,tau2R,b,df,D,
				      oldClique,oldComponents);
      cout << "updateOmega_HyperInverseWishart: " << nTry << " " << nAccept << endl;


      nTry = 1;
      nAccept = 0;
      updateDDeltaStar_HyperInverseWishart(&seed,&nAccept,Delta,*Q,*G,S,x,psi,nu,delta,r,sigma2,phi,tau2R,b,
					   Omega,oldClique,oldComponents);
      cout << "updateDeltaStar_HyperInverseWishart: " << nTry << " " << nAccept << endl;



      nTry = 1;
      nAccept = 0;
      updateGamma2(&seed,&nAccept,&(*gamma2),*Q,*G,nu,rho,sigma2,tau2Rho,a);
      cout << "updateGamma2: " << nTry << " " << nAccept << endl;
      cout << "gamma2: " << gamma2 << endl;


      nTry = 4;
      nAccept = 0;
      double epsilonGamma2Nu = 0.2;
      updateGamma2Nu(&seed,nTry,&nAccept,epsilonGamma2Nu,&(*gamma2),nu,*Q,*G,S,x,psi,
		     delta,Delta,rho,sigma2,phi,tau2Rho,a);
      cout << "updateGamma2Nu: " << nTry << " " << nAccept << endl;
      cout << "gamma2: " << gamma2 << endl;


      nTry = *Q * (*Q - 1) / 2;
      nAccept = 0;
      double epsilonRhoGamma2 = 0.2;
      updateRhoGamma2(&seed,nTry,&nAccept,epsilonRhoGamma2,rho,&(*gamma2),
		      *Q,*G,nu,sigma2,tau2Rho,a,*nuRho);
      cout << "updateRhoGamma2: " << nTry << " " << nAccept << endl;
      cout << "gamma2: " << gamma2 << endl;
      cout << "rho: ";
      for (q1 = 0; q1 < *Q; q1++)
	for (q2 = q1 + 1; q2 < *Q; q2++) {
	  int kqq = qq2index(q1,q2,*Q);
	  cout << rho[kqq] << " ";
	}
      cout << endl;


      nTry = *Q * (*Q - 1) / 2;
      nAccept = 0;
      double epsilonRhoNu = 0.2;
      updateRhoNu(&seed,nTry,&nAccept,epsilonRhoNu,rho,nu,
		  *Q,*G,S,x,psi,delta,Delta,*gamma2,sigma2,phi,tau2Rho,a,*nuRho);
      cout << "updateRhoNu: " << nTry << " " << nAccept << endl;
      cout << "rho: ";
      for (q1 = 0; q1 < *Q; q1++)
	for (q2 = q1 + 1; q2 < *Q; q2++) {
	  int kqq = qq2index(q1,q2,*Q);
	  cout << rho[kqq] << " ";
	}
      cout << endl;


      nTry = *G;
      nAccept = 0;
      double epsilonSigma2 = 0.5;
      updateSigma2_HyperInverseWishart(&seed,nTry,&nAccept,epsilonSigma2,sigma2,*Q,*G,S,x,psi,nu,
				       delta,Delta,*gamma2,r,rho,phi,t,l,tau2R,tau2Rho,a,b,
				       Omega,oldClique,oldComponents);
      cout << "updateSigma2_HyperInverseWishart: " << nTry << " " << nAccept << endl;


      nTry = 5 * *Q * *G;
      nAccept = 0;
      double epsilonPhi = 0.5;
      updatePhi(&seed,nTry,&nAccept,epsilonPhi,phi,*Q,*G,S,x,psi,nu,
		delta,Delta,sigma2,theta,lambda);
      cout << "updatePhi: " << nTry << " " << nAccept << endl;


      nTry = 250 * *Q;
      nAccept = 0;
      double epsilonTheta = 0.4;
      updateTheta(&seed,nTry,&nAccept,epsilonTheta,theta,*Q,*G,phi,lambda);
      cout << "updateTheta: " << nTry << " " << nAccept << endl;
      cout << "theta: ";
      for (q = 0; q < *Q; q++)
	cout << theta[q] << " ";
      cout << endl;



      nTry = 250 * *Q;
      nAccept = 0;
      double epsilonLambda = 0.01;
      updateLambda(&seed,nTry,&nAccept,epsilonLambda,lambda,*Q,*G,phi,theta);
      cout << "updateLambda: " << nTry << " " << nAccept << endl;
      cout << "lambda: ";
      for (q = 0; q < *Q; q++)
	cout << lambda[q] << " ";
      cout << endl;


      nTry = 10 * *Q;
      nAccept = 0;
      double epsilonLambdaPhi = 0.025;
      updateLambdaPhi(&seed,nTry,&nAccept,epsilonLambdaPhi,lambda,phi,*Q,*G,S,
		    x,psi,nu,delta,Delta,sigma2,theta);
      cout << "updateLambdaPhi: " << nTry << " " << nAccept << endl;
      cout << "lambda: ";
      for (q = 0; q < *Q; q++)
	cout << lambda[q] << " ";
      cout << endl;



      nTry = 10 * *Q;
      nAccept = 0;
      double epsilonThetaPhi = 0.05;
      updateThetaPhi(&seed,nTry,&nAccept,epsilonThetaPhi,theta,phi,*Q,*G,S,
		     x,psi,nu,delta,Delta,sigma2,lambda);
      cout << "updateThetaPhi: " << nTry << " " << nAccept << endl;
      cout << "theta: ";
      for (q = 0; q < *Q; q++)
	cout << theta[q] << " ";
      cout << endl;


      nTry = 250 * *Q;
      nAccept = 0;
      double epsilonT = 0.4;
      updateT(&seed,nTry,&nAccept,epsilonT,t,*Q,*G,sigma2,l);
      cout << "updateT: " << nTry << " " << nAccept << endl;
      cout << "t: ";
      for (q = 0; q < *Q; q++)
	cout << t[q] << " ";
      cout << endl;


      nTry = 250 * *Q;
      nAccept = 0;
      double epsilonL = 0.05;
      updateL(&seed,nTry,&nAccept,epsilonL,l,*Q,*G,sigma2,t);
      cout << "updateL: " << nTry << " " << nAccept << endl;
      cout << "l: ";
      for (q = 0; q < *Q; q++)
	cout << l[q] << " ";
      cout << endl;



      nTry = 10 * *Q;
      nAccept = 0;
      double epsilonLSigma2 = 0.025;
      updateLSigma2_HyperInverseWishart(&seed,nTry,&nAccept,epsilonLSigma2,l,sigma2,*Q,*G,S,x,psi,nu,
					delta,Delta,*gamma2,r,rho,phi,t,tau2R,tau2Rho,a,b,
					Omega,oldClique,oldComponents);
      cout << "updateLSigma2_HyperInverseWishart: " << nTry << " " << nAccept << endl;
      cout << "l: ";
      for (q = 0; q < *Q; q++)
	cout << l[q] << " ";
      cout << endl;


      nTry = 10 * *Q;
      nAccept = 0;
      double epsilonTSigma2 = 0.05;
      updateTSigma2_HyperInverseWishart(&seed,nTry,&nAccept,epsilonTSigma2,t,sigma2,*Q,*G,S,x,psi,nu,
					delta,Delta,*gamma2,r,rho,phi,l,tau2R,tau2Rho,a,b,
					Omega,oldClique,oldComponents);
      cout << "updateTSigma2_HyperInverseWishart: " << nTry << " " << nAccept << endl;
      cout << "t: ";
      for (q = 0; q < *Q; q++)
	cout << t[q] << " ";
      cout << endl;


      int nCorrect = 0;
      int tt;
      for (tt = 0; tt < *Q * *G; tt++)
	nCorrect += (delta[tt] == deltaTrue[tt]);
      cout << "fraction correct: " << ((double) nCorrect) / ((double) (*Q * *G)) << endl;

      double pot = potentialX(*Q,*G,S,x,psi,nu,delta,Delta,sigma2,phi) +
	potentialNu(*Q,*G,nu,*gamma2,a,rho,tau2Rho,sigma2) +
	potentialDDeltaStar_HyperInverseWishart(Delta,b,sigma2,tau2R,r,*Q,*G,Omega,oldClique,oldComponents) +
	potentialA(*Q,a,*pA0,*pA1,*alphaA,*betaA) +
	potentialB(*Q,b,*pB0,*pB1,*alphaB,*betaB) +
	potentialR(*Q,r,*nuR) +
	potentialRho(*Q,rho,*nuRho) +
	potentialSigma2(*Q,*G,sigma2,l,t) +
	potentialPhi(*Q,*G,phi,lambda,theta);


      if (oneDelta == 1)
	pot += potentialDelta_onedelta(*Q,*G,delta,xi) +
	  potentialXi_onedelta(xi,*alphaXi,*betaXi);
      else
	pot += potentialDelta(*Q,*G,delta,xi) +
	  potentialXi(*Q,xi,*alphaXi,*betaXi);

      cout << endl << "pot: " << pot << " " <<
	potentialX(*Q,*G,S,x,psi,nu,delta,Delta,sigma2,phi)	<< " " <<
	((double) nCorrect) / ((double) (*Q * *G)) << endl << endl;

      cout << "potentialX: " <<
	potentialX(*Q,*G,S,x,psi,nu,delta,Delta,sigma2,phi) << endl;
      cout << "potentialNu: " <<
	potentialNu(*Q,*G,nu,*gamma2,a,rho,tau2Rho,sigma2) << endl;
      cout << "potentialDDeltaStar: " <<
	potentialDDeltaStar_HyperInverseWishart(Delta,b,sigma2,tau2R,r,*Q,*G,Omega,oldClique,oldComponents) << endl;
      cout << "potentialOmega: " <<
	potentialOmega_HyperInverseWishart(Omega,D,df,oldClique,oldComponents) << endl;
      cout << "potentialA: " << potentialA(*Q,a,*pA0,*pA1,*alphaA,*betaA) << endl;
      cout << "potentialB: " << potentialB(*Q,b,*pB0,*pB1,*alphaB,*betaB) << endl;
      cout << "potentialR: " << potentialR(*Q,r,*nuR) << endl;
      cout << "potentialRho: " << potentialRho(*Q,rho,*nuRho) << endl;
      cout << "potentialSigma2 " << potentialSigma2(*Q,*G,sigma2,l,t) << endl;
      cout << "potentialPhi: " << potentialPhi(*Q,*G,phi,lambda,theta) << endl;
      cout << "potentialDelta: " << potentialDelta(*Q,*G,delta,xi) << endl;
      cout << "potentialXi: " << potentialXi(*Q,xi,*alphaXi,*betaXi) << endl;
      cout << endl;

      if(*sigma2Out == 1){
	for (q = 0; q < *Q; q++)
	  for (g = 0; g < *G; g++) {
	    int k = qg2index(q,g,*Q,*G);
	    sigma2File << sigma2[k] << " ";
	  }
	sigma2File << endl;
	sigma2File.flush();
      }
      else {
	for (q = 0; q < *Q; q++)
	  for (g = 0; g < *G; g++) {
	    int k = qg2index(q,g,*Q,*G);
	    sigma2File << sigma2[k] << " ";
	  }
      }
    } // MH-iteration
    aFile.close();
    sigma2File.close();

    //return 0;
  }
}
