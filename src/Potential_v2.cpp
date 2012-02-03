#include <vector>

#include "Random_v2.h"
#include "Matrix_v2.h"
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
			      



double potentialDDeltaStar_HyperInverseWishart(const double *Delta,
					       const double *b,
					       const double *sigma2,
					       const double *tau2R,
					       const double *r,
					       int Q,int G,
					       const vector<vector<vector<double> > > &Omega,
					       const vector<int> &oldClique,
					       const vector<vector<int> > &oldComponents) {
  vector<vector<double> > zero;
  zero.resize(G);
  int q,g;
  for (g = 0; g < G; g++) {
    zero[g].resize(Q);
    for (q = 0; q < Q; q++)
      zero[g][q] = 0.0;
  }
  
  vector<vector<double> > R;
  R.resize(Q);
  for (q = 0; q < Q; q++) {
    R[q].resize(Q);
  }
  for (q = 0; q < Q; q++) {
    R[q][q] = tau2R[q];
    int p;
    for (p = q + 1; p < Q; p++) {
      R[q][p] = sqrt(tau2R[p] * tau2R[q]) * r[qq2index(p,q,Q)];
      R[p][q] = R[q][p];
    }
  }
  
  vector<vector<double> > DeltaStar;
  DeltaStar.resize(G);
  for (g = 0; g < G; g++) {
    DeltaStar[g].resize(Q);
    for (q = 0; q < Q; q++) {
      DeltaStar[g][q] = Delta[qg2index(q,g,Q,G)] / exp(0.5 * b[q] * log(sigma2[qg2index(q,g,Q,G)]));
    }
  }

  Random ran(1);
  double pot = ran.PotentialMatrixVariateNormal(zero,R,Omega,oldClique,oldComponents,DeltaStar);
  
  return pot;
}
  



double potentialDDeltaStar_HyperInverseWishart(int gene,
					       const double *Delta,
					       const double *b,
					       const double *sigma2,
					       const double *tau2R,
					       const double *r,
					       int Q,int G,
					       const vector<vector<vector<double> > > &Omega,
					       const vector<int> &oldClique,
					       const vector<vector<int> > &oldComponents) {
  vector<vector<double> > R;
  R.resize(Q);
  int q,g;
  for (q = 0; q < Q; q++) {
    R[q].resize(Q);
  }
  for (q = 0; q < Q; q++) {
    R[q][q] = tau2R[q];
    int p;
    for (p = q + 1; p < Q; p++) {
      R[q][p] = sqrt(tau2R[p] * tau2R[q]) * r[qq2index(p,q,Q)];
      R[p][q] = R[q][p];
    }
  }
  
  vector<vector<double> > DeltaStar;
  DeltaStar.resize(G);
  for (g = 0; g < G; g++) {
    DeltaStar[g].resize(Q);
    for (q = 0; q < Q; q++) {
      DeltaStar[g][q] = Delta[qg2index(q,g,Q,G)] / exp(0.5 * b[q] * log(sigma2[qg2index(q,g,Q,G)]));
    }
  }

  Random ran(1);
  double pot = 0.0;
  
  vector<vector<double> > UU(DeltaStar);

  // allocate space and initialise temporal storage of U in blocks

  int k,i,j;
  vector<vector<vector<double> > > UBlocks;
  UBlocks.resize(Omega.size());
  for (k = 0; k < UBlocks.size(); k++) {
    UBlocks[k].resize(Omega[k].size());
    for (i = 0; i < UBlocks[k].size(); i++)
      UBlocks[k][i].resize(R.size());
  }

  vector<vector<int> > theGene;
  theGene.resize(UBlocks.size());
  for (k = 0; k < theGene.size(); k++) 
    theGene[k].resize(UBlocks[k].size());  


  int first = 0;
  for (i = 0; i < Omega[0].size(); i++)
    for (j = 0; j < UU[first].size(); j++) {
      UBlocks[0][i][j] = UU[first + i][j];
      theGene[0][i] = (first + i == gene);
    }
  first += Omega[0].size();

  for (k = 1; k < Omega.size(); k++) {
    for (i = 0; i < oldComponents[k].size(); i++)
      for (j = 0; j < UU[first].size(); j++) {
	UBlocks[k][i][j] = UBlocks[oldClique[k]][oldComponents[k][i]][j];
	theGene[k][i] = theGene[oldClique[k]][oldComponents[k][i]];
      }
    
    for (i = 0; i < Omega[k].size() - oldComponents[k].size(); i++)
      for (j = 0; j < UU[first].size(); j++) {
	UBlocks[k][i + oldComponents[k].size()][j] = UU[first + i][j];
	theGene[k][i + oldComponents[k].size()] = (first + i == gene);
      }
    first += Omega[k].size() - oldComponents[k].size();
  }
  
  // add potential for each clique
  
  for (k = 0; k < Omega.size(); k++) {
    int include = 0;
    for (i = 0; i < theGene[k].size(); i++)
      if (theGene[k][i] == 1) include = 1;
    
    if (include == 1)
      pot += ran.PotentialMatrixVariateNormal(Omega[k],R,UBlocks[k]);
  }
  
  // subtract potential for each separator

  for (k = 1; k < Omega.size(); k++) 
    if (oldComponents[k].size() > 0) {
      int include = 0;

      vector<vector<double> > OmegaSub;
      vector<vector<double> > USub;
      OmegaSub.resize(oldComponents[k].size());
      USub.resize(oldComponents[k].size());
      for (i = 0; i < OmegaSub.size(); i++) {
	OmegaSub[i].resize(oldComponents[k].size());
	for (j = 0; j < OmegaSub[i].size(); j++)
	  OmegaSub[i][j] = Omega[k][i][j];
      }
      for (i = 0; i < USub.size(); i++) {
	USub[i].resize(UBlocks[k][i].size());
	for (j = 0; j < USub[i].size(); j++)
	  USub[i][j] = UBlocks[k][i][j];

	if (theGene[k][i] == 1) include == 1;
      }
      
      if (include == 1)
	pot -= ran.PotentialMatrixVariateNormal(OmegaSub,R,USub);
    }
  

  return pot;
}
  



double potentialOmega_HyperInverseWishart(const vector<vector<vector<double> > > &Omega,
					  const vector<vector<vector<double> > > &D,
					  double df,
					  const vector<int> &oldClique,
					  const vector<vector<int> > &oldComponents) {
  Random ran(1);
  double pot = ran.PotentialHyperInverseWishart(df,D,oldClique,oldComponents,Omega);

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



double potentialDelta_MRF1_onedelta(int Q,
				    int G,
				    const int *delta,
				    const vector<vector<int> > &neighbour,
				    double eta0,
				    double omega0,
				    double kappa) {
  int *dd = (int *) calloc(G,sizeof(int));
  int g;
  for (g = 0; g < G; g++) {
    int kqg = qg2index(0,g,Q,G);
    dd[g] = delta[kqg];
  }
  vector<double> potMarg(G,0);

  unsigned int dummy = 1;
  double pot = perfectMRF1_onedelta(dd,G,neighbour,potMarg,potMarg,
				    eta0,omega0,kappa,&dummy,0);
  free(dd);

  return pot;
}



double potentialDelta_MRF2_onedelta(int Q,
				    int G,
				    const int *delta,
				    const vector<vector<int> > &neighbour,
				    double alpha,
				    double beta) {
  int *dd = (int *) calloc(G,sizeof(int));
  int g;
  for (g = 0; g < G; g++) {
    int kqg = qg2index(0,g,Q,G);
    dd[g] = delta[kqg];
  }
  vector<double> potMarg(G,0);

  unsigned int dummy = 1;
  double pot = perfectMRF2_onedelta(dd,G,neighbour,potMarg,potMarg,
				    alpha,beta,&dummy,0);
  free(dd);

  return pot;
}



double potentialDelta_MRF2(int Q,
			   int G,
			   const int *delta,
			   const vector<vector<int> > &neighbour,
			   double alpha,
			   double beta,
			   double betag) {
  vector<double> potMarg(Q * G,0);

  int *dd = (int *) calloc(Q * G,sizeof(int));
  int k;
  for (k = 0; k < Q * G; k++)
    dd[k] = delta[k];

  unsigned int dummy = 1;
  double pot = perfectMRF2(dd,Q,G,neighbour,potMarg,potMarg,
			   alpha,beta,betag,&dummy,0);

  free(dd);

  return pot;
}



double potentialEta0(double eta0,double alphaEta,double betaEta) {
  double pot = 0.0;

  unsigned int seed = 1;
  Random ran(seed);

  pot += ran.PotentialBeta(alphaEta,betaEta,eta0);

  return pot;
}




double potentialOmega0(double omega0,double pOmega0,double lambdaOmega) {
  double pot = 0.0;
  
  if (omega0 = 0.0)
    pot += - log(pOmega0);
  else {
    pot += - log(1.0 - pOmega0);
    //    pot += - log(lambdaOmega) + lambdaOmega * omega0;
  }

  return pot;
}




double potentialKappa(double kappa,double lambdaKappa) {
  double pot = 0.0;
  
  pot += - log(lambdaKappa) + lambdaKappa * kappa;

  pot = 0.0;
    
  return pot;
}




double potentialAlpha(void) {
  return 0.0;
}


double potentialBeta(void) {
  return 0.0;
}



double potentialBetag(void) {
  return 0.0;
}
