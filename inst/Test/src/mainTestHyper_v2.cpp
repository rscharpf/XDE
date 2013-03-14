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
  void initializeParams(int *nIt,
			int *seedR,
			int *G,
			int *Q,
			int *S,
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
			double *aIn,
			double *bIn,
			double *lIn,
			double *tIn,
			double *lambdaIn,
			double *thetaIn,
			double *xiIn,
			int *simulateExpression, // whether to simulate data from the prior
			double *x,
			int *writeout,
			int *simulateSigma2,
			int *oldCliqueInput,
			int *oldComponentsInput,
			int *nClique,
			int *nOldClique,
			int *nTotalInCliqueInput,
			int *clique,
			int *nNewInCliqueInput){//22
    unsigned int seed;
    seed = *seedR;
    Random ran(seed);
    int sumS = 0;
    int q;
    for (q = 0; q < *Q; q++)
      sumS += S[q];

    vector<int> oldClique;
    vector<vector<int> > oldComponents;
    vector<int> nNewInClique(0,0);
    int g;
    nNewInClique.resize(clique[*G - 1] + 1);
    int i;
    for (i = 0; i < nNewInClique.size(); i++)
      nNewInClique[i] = nNewInCliqueInput[i];
    oldClique.resize(nNewInClique.size());
    for(i = 0; i < nNewInClique.size(); i++)
      oldClique[i] = oldCliqueInput[i];
    oldComponents.resize(nNewInClique.size());
    oldComponents[0].resize(0);
    vector<int> nTotalInClique(nNewInClique);
    int c;
    for (c = 1; c < nNewInClique.size(); c++) {
      int nOld=nTotalInCliqueInput[c] - nNewInClique[c];
      nTotalInClique[c] = nTotalInCliqueInput[c];
      oldComponents[c].resize(nOld);
      int nSampled=0;
      int nr = 0;
      while(nSampled < nOld){
	oldComponents[c][nSampled]=oldComponentsInput[nr];
	nSampled++;
	nr++;
      }
    }
    for (c = 0; c < nTotalInClique.size(); c++)
      cout << "nTotalInClique: " << c << " " << nTotalInClique[c] << endl;
    // initialise clinical variables

    int s;

    double *a = (double *) calloc(*Q,sizeof(double));
    double *b = (double *) calloc(*Q,sizeof(double));
    double *l = (double *) calloc(*Q,sizeof(double));
    double *t = (double *) calloc(*Q,sizeof(double));
    double *lambda = (double *) calloc(*Q,sizeof(double));
    double *theta = (double *) calloc(*Q,sizeof(double));
    double *xi = (double *) calloc(*Q, sizeof(double));

    for (q = 0; q < *Q; q++) {
      a[q] = aIn[q];
      b[q] = bIn[q];
      l[q] = lIn[q];
      t[q] = tIn[q];
      lambda[q] = lambdaIn[q];
      theta[q] = thetaIn[q];
      xi[q]=xiIn[q];
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

    // Simulating expression data
    ofstream xfile("x.txt");
    //double *x = (double *) calloc(*G * sumS,sizeof(double));
    if(*simulateExpression == 1){
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
  ofstream bFile("b.txt");
  ofstream lFile("l.txt");
  ofstream lambdaFile("lambda.txt");
  ofstream rFile("r.txt");
  ofstream rhoFile("rho.txt");
  ofstream sigma2File("sigma2.txt");
  ofstream tFile("t.txt");
  ofstream tau2RhoFile("tau2Rho.txt");
  ofstream tau2RFile("tau2R.txt");
  ofstream thetaFile("theta.txt");
  ofstream xiFile("xi.txt");

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

      updateANu(&seed,nTry,&nAccept,epsilonANu,a,nu,*Q,*G,S,x,psi,delta,Delta,
		*gamma2,rho,sigma2,phi,tau2Rho,*pA0,*pA1,*alphaA,*betaA);

      if (oneDelta == 1) {
	nTry = *G;
	nAccept = 0;
	updateDelta_HyperInverseWishart_onedelta(&seed,nTry,&nAccept,delta,*Q,*G,S,x,psi,
						 nu,Delta,r,sigma2,phi,xi,b);
	//cout << "updateDelta_onedelta: " << nTry << " " << nAccept << endl;

	nTry = *Q;
	nAccept = 0;
	updateXi_onedelta(&seed,&nAccept,xi,*Q,*G,delta,*alphaXi,*betaXi);

      }
      else {
	nTry = *Q * *G;
	nAccept = 0;
	updateDelta_HyperInverseWishart(&seed,nTry,&nAccept,delta,*Q,*G,S,x,psi,
					nu,Delta,r,sigma2,phi,xi,b);
	//cout << "updateDelta: " << nTry << " " << nAccept << endl;


	nTry = *Q;
	nAccept = 0;
	updateXi(&seed,&nAccept,xi,*Q,*G,delta,*alphaXi,*betaXi);
      }

      nTry = *Q;
      nAccept = 0;
      double epsilonBDDelta = 0.1;
      updateBDDeltaStar_HyperInverseWishart(&seed,nTry,&nAccept,epsilonBDDelta,b,Delta,*Q,*G,S,x,psi,
					    nu,delta,r,sigma2,phi,tau2R,*pB0,*pB1,*alphaB,*betaB,Omega,
					    oldClique,oldComponents);

      nTry = *Q * (*Q - 1) / 2;
      nAccept = 0;
      double epsilonRDDelta = 0.1;
      updateRDDeltaStar_HyperInverseWishart(&seed,nTry,&nAccept,epsilonRDDelta,r,Delta,*Q,*G,S,x,psi,
					    nu,delta,sigma2,phi,tau2R,b,*nuR,Omega,
					    oldClique,oldComponents);

      nTry = *Q;
      nAccept = 0;
      double epsilonTau2RhoNu = 0.02;
      updateTau2RhoNu(&seed,nTry,&nAccept,epsilonTau2RhoNu,tau2Rho,
		      nu,*Q,*G,S,x,psi,delta,Delta,*gamma2,rho,sigma2,
		      phi,a);

      nTry = *Q;
      nAccept = 0;
      double epsilonTau2RDDelta = 0.02;
      updateTau2RDDeltaStar_HyperInverseWishart(&seed,nTry,&nAccept,epsilonTau2RDDelta,tau2R,
						Delta,*Q,*G,S,x,psi,nu,delta,r,sigma2,phi,b,Omega,
						oldClique,oldComponents);

      nTry = *Q * *G;
      nAccept = 0;
      updateNu(&seed,&nAccept,nu,*Q,*G,S,x,psi,delta,Delta,*gamma2,rho,sigma2,
	       phi,tau2Rho,a);
      //cout << "updateNu: " << nTry << " " << nAccept << endl;

      nTry = 1;
      nAccept = 0;
      updateOmega_HyperInverseWishart(&seed,&nAccept,Omega,*Q,*G,Delta,r,sigma2,tau2R,b,df,D,
				      oldClique,oldComponents);
      cout << "updateOmega_HyperInverseWishart: " << nTry << " " << nAccept << endl;


      nTry = 1;
      nAccept = 0;
      updateDDeltaStar_HyperInverseWishart(&seed,&nAccept,Delta,*Q,*G,S,x,psi,nu,delta,r,sigma2,phi,tau2R,b,
					   Omega,oldClique,oldComponents);
      //cout << "updateDeltaStar_HyperInverseWishart: " << nTry << " " << nAccept << endl;

      nTry = 1;
      nAccept = 0;
      updateGamma2(&seed,&nAccept,&(*gamma2),*Q,*G,nu,rho,sigma2,tau2Rho,a);
      /*
      cout << "updateGamma2: " << nTry << " " << nAccept << endl;
      cout << "gamma2: " << gamma2 << endl;
      */


      nTry = 4;
      nAccept = 0;
      double epsilonGamma2Nu = 0.2;
      updateGamma2Nu(&seed,nTry,&nAccept,epsilonGamma2Nu,&(*gamma2),nu,*Q,*G,S,x,psi,
		     delta,Delta,rho,sigma2,phi,tau2Rho,a);
      /*
      cout << "updateGamma2Nu: " << nTry << " " << nAccept << endl;
      cout << "gamma2: " << gamma2 << endl;
      */


      nTry = *Q * (*Q - 1) / 2;
      nAccept = 0;
      double epsilonRhoGamma2 = 0.2;
      updateRhoGamma2(&seed,nTry,&nAccept,epsilonRhoGamma2,rho,&(*gamma2),
		      *Q,*G,nu,sigma2,tau2Rho,a,*nuRho);

      nTry = *Q * (*Q - 1) / 2;
      nAccept = 0;
      double epsilonRhoNu = 0.2;
      updateRhoNu(&seed,nTry,&nAccept,epsilonRhoNu,rho,nu,
		  *Q,*G,S,x,psi,delta,Delta,*gamma2,sigma2,phi,tau2Rho,a,*nuRho);

      nTry = *G;
      nAccept = 0;
      double epsilonSigma2 = 0.5;
      updateSigma2_HyperInverseWishart(&seed,nTry,&nAccept,epsilonSigma2,sigma2,*Q,*G,S,x,psi,nu,
				       delta,Delta,*gamma2,r,rho,phi,t,l,tau2R,tau2Rho,a,b,
				       Omega,oldClique,oldComponents);
      //cout << "updateSigma2_HyperInverseWishart: " << nTry << " " << nAccept << endl;


      nTry = 5 * *Q * *G;
      nAccept = 0;
      double epsilonPhi = 0.5;
      updatePhi(&seed,nTry,&nAccept,epsilonPhi,phi,*Q,*G,S,x,psi,nu,
		delta,Delta,sigma2,theta,lambda);
      //cout << "updatePhi: " << nTry << " " << nAccept << endl;


      nTry = 250 * *Q;
      nAccept = 0;
      double epsilonTheta = 0.4;
      updateTheta(&seed,nTry,&nAccept,epsilonTheta,theta,*Q,*G,phi,lambda);

      nTry = 250 * *Q;
      nAccept = 0;
      double epsilonLambda = 0.01;
      updateLambda(&seed,nTry,&nAccept,epsilonLambda,lambda,*Q,*G,phi,theta);

      nTry = 10 * *Q;
      nAccept = 0;
      double epsilonLambdaPhi = 0.025;
      updateLambdaPhi(&seed,nTry,&nAccept,epsilonLambdaPhi,lambda,phi,*Q,*G,S,
		    x,psi,nu,delta,Delta,sigma2,theta);

      nTry = 10 * *Q;
      nAccept = 0;
      double epsilonThetaPhi = 0.05;
      updateThetaPhi(&seed,nTry,&nAccept,epsilonThetaPhi,theta,phi,*Q,*G,S,
		     x,psi,nu,delta,Delta,sigma2,lambda);

      nTry = 250 * *Q;
      nAccept = 0;
      double epsilonT = 0.4;
      updateT(&seed,nTry,&nAccept,epsilonT,t,*Q,*G,sigma2,l);

      nTry = 250 * *Q;
      nAccept = 0;
      double epsilonL = 0.05;
      updateL(&seed,nTry,&nAccept,epsilonL,l,*Q,*G,sigma2,t);

      nTry = 10 * *Q;
      nAccept = 0;
      double epsilonLSigma2 = 0.025;
      updateLSigma2_HyperInverseWishart(&seed,nTry,&nAccept,epsilonLSigma2,l,sigma2,*Q,*G,S,x,psi,nu,
					delta,Delta,*gamma2,r,rho,phi,t,tau2R,tau2Rho,a,b,
					Omega,oldClique,oldComponents);

      nTry = 10 * *Q;
      nAccept = 0;
      double epsilonTSigma2 = 0.05;
      updateTSigma2_HyperInverseWishart(&seed,nTry,&nAccept,epsilonTSigma2,t,sigma2,*Q,*G,S,x,psi,nu,
					delta,Delta,*gamma2,r,rho,phi,l,tau2R,tau2Rho,a,b,
					Omega,oldClique,oldComponents);

      int nCorrect = 0;
      int tt;
      for (tt = 0; tt < *Q * *G; tt++)
	nCorrect += (delta[tt] == deltaTrue[tt]);
      //cout << "fraction correct: " << ((double) nCorrect) / ((double) (*Q * *G)) << endl;

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

      if(writeout[1] == 1){
	for (q = 0; q < *Q; q++){
	  aFile << a[q] << " ";
	  bFile << b[q] << " ";
	  lFile << l[q] << " ";
	  lambdaFile << lambda[q] << " ";
	  rFile << r[q] << " ";
	  rhoFile << rho[q] << " ";
	  tFile << t[q] << " ";
	  tau2RFile << tau2R[q] << " ";
	  tau2RhoFile << tau2Rho[q] << " ";
	  thetaFile << theta[q] << " ";
	  xiFile << xi[q] << " ";
	}
	aFile << endl;
	aFile.flush();
	bFile << endl;
	bFile.flush();
	lFile << endl;
	lFile.flush();
	lambdaFile << endl;
	lambdaFile.flush();
	rFile << endl;
	rFile.flush();
	tFile << endl;
	tFile.flush();
	rhoFile << endl;
	rhoFile.flush();
	tau2RFile << endl;
	tau2RFile.flush();
	tau2RhoFile << endl;
	tau2RhoFile.flush();
	thetaFile << endl;
	thetaFile.flush();
	xiFile << endl;
	xiFile.flush();
	for (g = 0; g < *G; g++) {
	  int k = qg2index(q,g,*Q,*G);
	  sigma2File << sigma2[k] << " ";
	}
	sigma2File << endl;
	sigma2File.flush();
      }
    } // MH-iteration
    aFile.close();
    bFile.close();
    lFile.close();
    lambdaFile.close();
    rFile.close();
    tFile.close();
    rhoFile.close();
    sigma2File.close();
    tau2RhoFile.close();
    tau2RFile.close();
    thetaFile.close();
    xiFile.close();
  }
}
