#include <R.h>
#include "Update_v2.h"
#include "Utility_v2.h"

extern "C" {
  
  // updates that can be used for all model variants

  // nu|... ~ N(0, Sigma)
  // a is the power conjugate
  void updateANu(int *seed,
  	 const int *nTry,
		 int *nAccept,
		 const double *epsilon,
		 double *a,
		 double *nu,
		 const int *Q,
		 const int *G,
		 const int *S,
		 const double *x,
		 const int *psi,
		 const int *delta,
		 const double *Delta,
		 const double *gamma2,
		 const double *rho,
		 const double *sigma2,
		 const double *phi,
		 const double *tau2Rho,
		 const double *pA0,
		 const double *pA1,
		 const double *alphaA,
		 const double *betaA) {
    unsigned int seedU = (unsigned int) *seed;

    updateANu(&seedU,*nTry,nAccept,*epsilon,a,nu,*Q,*G,S,x,psi,delta,Delta,
	      *gamma2,rho,sigma2,phi,tau2Rho,*pA0,*pA1,*alphaA,*betaA);

    *seed = (int) seedU;

    return;
  }


  void updateTau2RhoNu(int *seed,
		       const int *nTry,
		       int *nAccept,
		       const double *epsilon,
		       double *tau2Rho,
		       double *nu,
		       const int *Q,
		       const int *G,
		       const int *S,
		       const double *x,
		       const int *psi,
		       const int *delta,
		       const double *Delta,
		       const double *gamma2,
		       const double *rho,
		       const double *sigma2,
		       const double *phi,
		       const double *a) {
    unsigned int seedU = (unsigned int) *seed;

    updateTau2RhoNu(&seedU,*nTry,nAccept,*epsilon,tau2Rho,nu,*Q,*G,S,x,
		    psi,delta,Delta,*gamma2,rho,sigma2,phi,a);

    *seed = (int) seedU;

    return;
  }


  void updateNu(int *seed,
		int *nAccept,
		double *nu,
		const int *Q,
		const int *G,
		const int *S,
		const double *x,
		const int *psi,
		const int *delta,
		const double *Delta,
		const double *gamma2,
		const double *rho,
		const double *sigma2,
		const double *phi,
		const double *tau2Rho,
		const double *a) {
    unsigned int seedU = (unsigned int) *seed;

    updateNu(&seedU,nAccept,nu,*Q,*G,S,x,psi,delta,Delta,*gamma2,
	     rho,sigma2,phi,tau2Rho,a);

    *seed = (int) seedU;

    return;
  }


  void updateGamma2(int *seed,
		    int *nAccept,
		    double *gamma2,
		    const int *Q,
		    const int *G,
		    const double *nu,
		    const double *rho,
		    const double *sigma2,
		    const double *tau2Rho,
		    const double *a) {
    unsigned int seedU = (unsigned int) *seed;

    updateGamma2(&seedU,nAccept,gamma2,*Q,*G,nu,rho,sigma2,tau2Rho,a);

    *seed = (int) seedU;

    return;
  }


  void updateRhoNu(int *seed,
		   const int *nTry,
		   int *nAccept,
		   const double *epsilon,
		   double *rho,
		   double *nu,
		   const int *Q,
		   const int *G,
		   const int *S,
		   const double *x,
		   const int *psi,
		   const int *delta,
		   const double *Delta,
		   const double *gamma2,
		   const double *sigma2,
		   const double *phi,
		   const double *tau2Rho,
		   const double *a,
		   const double *nuRho) {
    unsigned int seedU = (unsigned int) *seed;

    updateRhoNu(&seedU,*nTry,nAccept,*epsilon,rho,nu,*Q,*G,S,x,psi,delta,
		Delta,*gamma2,sigma2,phi,tau2Rho,a,*nuRho);

    *seed = (int) seedU;

    return;
  }


  void updateRhoGamma2(int *seed,
		       const int *nTry,
		       int *nAccept,
		       const double *epsilon,
		       double *rho,
		       double *gamma2,
		       const int *Q,
		       const int *G,
		       const double *nu,
		       const double *sigma2,
		       const double *tau2Rho,
		       const double *a,
		       const double *nuRho) {
    unsigned int seedU = (unsigned int) *seed;

    updateRhoGamma2(&seedU,*nTry,nAccept,*epsilon,rho,gamma2,*Q,*G,nu,
	      sigma2,tau2Rho,a,*nuRho);

    *seed = (int) seedU;

    return;
  }


  void updatePhi(int *seed,
		 const int *nTry,
		 int *nAccept,
		 const double *epsilon,
		 double *phi,
		 const int *Q,
		 const int *G,
		 const int *S,
		 const double *x,
		 const int *psi,
		 const double *nu,
		 const int *delta,
		 const double *Delta,
		 const double *sigma2,
		 const double *theta,
		 const double *lambda) {
    unsigned int seedU = (unsigned int) *seed;

    updatePhi(&seedU,*nTry,nAccept,*epsilon,phi,*Q,*G,S,x,psi,nu,
		 delta,Delta,sigma2,theta,lambda);

    *seed = (int) seedU;

    return;
  }


  void updateTheta(int *seed,
		   const int *nTry,
		   int *nAccept,
		   const double *epsilon,
		   double *theta,
		   const int *Q,
		   const int *G,
		   const double *phi,
		   const double *lambda) {
    unsigned int seedU = (unsigned int) *seed;

    updateTheta(&seedU,*nTry,nAccept,*epsilon,theta,*Q,*G,phi,lambda);

    *seed = (int) seedU;

    return;
  }


  void updateLambda(int *seed,
		    const int *nTry,
		    int *nAccept,
		    const double *epsilon,
		    double *lambda,
		    const int *Q,
		    const int *G,
		    const double *phi,
		    const double *theta) {
    unsigned int seedU = (unsigned int) *seed;

    updateLambda(&seedU,*nTry,nAccept,*epsilon,lambda,*Q,*G,phi,theta);

    *seed = (int) seedU;

    return;
  }


  void updateT(int *seed,
	       const int *nTry,
	       int *nAccept,
	       const double *epsilon,
	       double *t,
	       const int *Q,
	       const int *G,
	       const double *sigma2,
	       const double *l) {
    unsigned int seedU = (unsigned int) *seed;

    updateT(&seedU,*nTry,nAccept,*epsilon,t,*Q,*G,sigma2,l);

    *seed = (int) seedU;

    return;
  }


  void updateL(int *seed,
	       const int *nTry,
	       int *nAccept,
	       const double *epsilon,
	       double *l,
	       const int *Q,
	       const int *G,
	       const double *sigma2,
	       const double *t) {
    unsigned int seedU = (unsigned int) *seed;

    updateL(&seedU,*nTry,nAccept,*epsilon,l,*Q,*G,sigma2,t); 

    *seed = (int) seedU;

    return;
  }


  void updateLambdaPhi(unsigned int *seed,
		       const int *nTry,
		       int *nAccept,
		       const double *epsilon,
		       double *lambda,
		       double *phi,
		       const int *Q,
		       const int *G,
		       const int *S,
		       const double *x,
		       const int *psi,
		       const double *nu,
		       const int *delta,
		       const double *Delta,
		       const double *sigma2,
		       const double *theta) {
    unsigned int seedU = (unsigned int) *seed;

    updateLambdaPhi(&seedU,*nTry,nAccept,*epsilon,lambda,phi,*Q,*G,S,x,psi,
		    nu,delta,Delta,sigma2,theta);

    *seed = (int) seedU;

    return;
  }


  void updateThetaPhi(unsigned int *seed,
		      const int *nTry,
		      int *nAccept,
		      const double *epsilon,
		      double *theta,
		      double *phi,
		      const int *Q,
		      const int *G,
		      const int *S,
		      const double *x,
		      const int *psi,
		      const double *nu,
		      const int *delta,
		      const double *Delta,
		      const double *sigma2,
		      const double *lambda) {
    unsigned int seedU = (unsigned int) *seed;

    updateThetaPhi(&seedU,*nTry,nAccept,*epsilon,theta,phi,*Q,*G,S,x,psi,
		   nu,delta,Delta,sigma2,lambda);

    *seed = (int) seedU;

    return;
  }


  // finished: updates that can be used by all model variants


  // updates only for model A



  void updateXi_MA(int *seed,
		   int *nAccept,
		   double *xi,
		   const int *Q,
		   const int *G,
		   const int *delta,
		   const double *alphaXi,
		   const double *betaXi) {
    unsigned int seedU = (unsigned int) *seed;

    updateXi(&seedU,nAccept,xi,*Q,*G,delta,*alphaXi,*betaXi);

    *seed = (int) seedU;

    return;
  }


  // finished: updates that can be used only for model A


  // updates only for model B


  void updateXi_MB(int *seed,
		   int *nAccept,
		   double *xi,
		   const int *Q,
		   const int *G,
		   const int *delta,
		   const double *alphaXi,
		   const double *betaXi) {
    unsigned int seedU = (unsigned int) *seed;

    updateXi_onedelta(&seedU,nAccept,xi,*Q,*G,delta,*alphaXi,*betaXi);

    *seed = (int) seedU;

    return;
  }


  // finished: updates for model B
  
 
  // updates only for model C


  void updateAlpha_MC(unsigned int *seed,
		      const int *nTry,
		      int *nAccept,
		      const double *epsilon,
		      double *alpha,
		      const int *Q,
		      const int *G,
		      const int *delta,
		      const int *nNeighbour,
		      const int *neighbour,
		      const double *beta,
		      const double *betag) {
    unsigned int seedU = (unsigned int) *seed;
    
    vector<vector<int> > nn;
    nn.resize(*G);
    int g;
    for (g = 0; g < *G; g++) nn[g].resize(0);

    int i;
    for (i = 0; i < *nNeighbour; i++) {
      int g1 = neighbour[2 * i];
      int g2 = neighbour[2 * i + 1];
      nn[g1].push_back(g2);
      nn[g2].push_back(g1);
    }

    double betaCopy = *beta;
    double betagCopy = *betag;
    updateAlphaBetaBetag_MRF2(&seedU,*nTry,nAccept,*epsilon,0.0,0.0,alpha,
			      &betaCopy,&betagCopy,*Q,*G,delta,nn);

    *seed = (int) seedU;

    return;
  }


  void updateBeta_MC(unsigned int *seed,
		     const int *nTry,
		     int *nAccept,
		     const double *epsilon,
		     double *beta,
		     const int *Q,
		     const int *G,
		     const int *delta,
		     const int *nNeighbour,
		     const int *neighbour,
		     const double *alpha,
		     const double *betag) {
    unsigned int seedU = (unsigned int) *seed;

    vector<vector<int> > nn;
    nn.resize(*G);
    int g;
    for (g = 0; g < *G; g++) nn[g].resize(0);

    int i;
    for (i = 0; i < *nNeighbour; i++) {
      int g1 = neighbour[2 * i];
      int g2 = neighbour[2 * i + 1];
      nn[g1].push_back(g2);
      nn[g2].push_back(g1);
    }

    double alphaCopy = *alpha;
    double betagCopy = *betag;
    updateAlphaBetaBetag_MRF2(&seedU,*nTry,nAccept,0.0,*epsilon,0.0,&alphaCopy,
			      beta,&betagCopy,*Q,*G,delta,nn);

    *seed = (int) seedU;

    return;
  }


  void updateBetag_MC(unsigned int *seed,
		      const int *nTry,
		      int *nAccept,
		      const double *epsilon,
		      double *betag,
		      const int *Q,
		      const int *G,
		      const int *delta,
		      const int *nNeighbour,
		      const int *neighbour,
		      const double *alpha,
		      const double *beta) {
    unsigned int seedU = (unsigned int) *seed;

    vector<vector<int> > nn;
    nn.resize(*G);
    int g;
    for (g = 0; g < *G; g++) nn[g].resize(0);

    int i;
    for (i = 0; i < *nNeighbour; i++) {
      int g1 = neighbour[2 * i];
      int g2 = neighbour[2 * i + 1];
      nn[g1].push_back(g2);
      nn[g2].push_back(g1);
    }

    double alphaCopy = *alpha;
    double betaCopy = *beta;
    updateAlphaBetaBetag_MRF2(&seedU,*nTry,nAccept,0.0,0.0,*epsilon,&alphaCopy,
			      &betaCopy,betag,*Q,*G,delta,nn);

    *seed = (int) seedU;

    return;
  }


  // finished: updates for model C


  // updates only for model D


  void updateAlpha_MD(unsigned int *seed,
		      const int *nTry,
		      int *nAccept,
		      const double *epsilon,
		      double *alpha,
		      const int *Q,
		      const int *G,
		      const int *delta,
		      const int *nNeighbour,
		      const int *neighbour,
		      const double *beta) {
    unsigned int seedU = (unsigned int) *seed;

    vector<vector<int> > nn;
    nn.resize(*G);
    int g;
    for (g = 0; g < *G; g++) nn[g].resize(0);

    int i;
    for (i = 0; i < *nNeighbour; i++) {
      int g1 = neighbour[2 * i];
      int g2 = neighbour[2 * i + 1];
      nn[g1].push_back(g2);
      nn[g2].push_back(g1);
    }

    double betaCopy = *beta;
    updateAlphaBeta_MRF2_onedelta(&seedU,*nTry,nAccept,*epsilon,0.0,alpha,
				  &betaCopy,*Q,*G,delta,nn);

    *seed = (int) seedU;

    return;
  }


  void updateBeta_MD(unsigned int *seed,
		     const int *nTry,
		     int *nAccept,
		     const double *epsilon,
		     double *beta,
		     const int *Q,
		     const int *G,
		     const int *delta,
		     const int *nNeighbour,
		     const int *neighbour,
		     const double *alpha) {
    unsigned int seedU = (unsigned int) *seed;
    
    vector<vector<int> > nn;
    nn.resize(*G);
    int g;
    for (g = 0; g < *G; g++) nn[g].resize(0);
    
    int i;
    for (i = 0; i < *nNeighbour; i++) {
      int g1 = neighbour[2 * i];
      int g2 = neighbour[2 * i + 1];
      nn[g1].push_back(g2);
      nn[g2].push_back(g1);
    }
    
    double alphaCopy = *alpha;
    updateAlphaBeta_MRF2_onedelta(&seedU,*nTry,nAccept,0.0,*epsilon,&alphaCopy,
				  beta,*Q,*G,delta,nn);
    
    *seed = (int) seedU;
    
    return;
  }
  
  
  // finished: updates for model D



  // updates only for model I


  void updateBDDelta_MI(int *seed,
			const int *nTry,
			int *nAccept,
			const double *epsilon,
			double *b,
			double *Delta,
			const int *Q,
			const int *G,
			const int *S,
			const double *x,
			const int *psi,
			const double *nu,
			const int *delta,
			const double *c2,
			const double *r,
			const double *sigma2,
			const double *phi,
			const double *tau2R,
			const double *pB0,
			const double *pB1,
			const double *alphaB,
			const double *betaB) {
    unsigned int seedU = (unsigned int) *seed;
    
    updateBDDelta(&seedU,*nTry,nAccept,*epsilon,b,Delta,*Q,*G,S,x,psi,nu,
		  delta,*c2,r,sigma2,phi,tau2R,*pB0,*pB1,*alphaB,*betaB);

    *seed = (int) seedU;

    return;
  }

  

  void updateTau2RDDelta_MI(int *seed,
			    const int *nTry,
			    int *nAccept,
			    const double *epsilon,
			    double *tau2R,
			    double *Delta,
			    const int *Q,
			    const int *G,
			    const int *S,
			    const double *x,
			    const int *psi,
			    const double *nu,
			    const int *delta,
			    const double *c2,
			    const double *r,
			    const double *sigma2,
			    const double *phi,
			    const double *b) {
    unsigned int seedU = (unsigned int) *seed;

    updateTau2RDDelta(&seedU,*nTry,nAccept,*epsilon,tau2R,Delta,*Q,*G,S,x,
		      psi,nu,delta,*c2,r,sigma2,phi,b);

    *seed = (int) seedU;

    return;
  }


  void updateDDelta_MI(int *seed,
		       int *nAccept,
		       double *Delta,
		       const int *Q,
		       const int *G,
		       const int *S,
		       const double *x,
		       const int *psi,
		       const double *nu,
		       const int *delta,
		       const double *c2,
		       const double *r,
		       const double *sigma2,
		       const double *phi,
		       const double *tau2R,
		       const double *b) {
    unsigned int seedU = (unsigned int) *seed;

    updateDDelta(&seedU,nAccept,Delta,*Q,*G,S,x,psi,nu,delta,*c2,
		 r,sigma2,phi,tau2R,b);

    *seed = (int) seedU;

    return;
  }


  void updateC2_MI(int *seed,
		   const int *nTry,
		   int *nAccept,
		   double *c2,
		   const int *Q,
		   const int *G,
		   const int *delta,
		   const double *Delta,
		   const double *r,
		   const double *sigma2,
		   const double *tau2R,
		   const double *b,
		   const double *c2Max) {
    unsigned int seedU = (unsigned int) *seed;

    updateC2(&seedU,*nTry,nAccept,c2,*Q,*G,delta,Delta,r,
	     sigma2,tau2R,b,*c2Max);

    *seed = (int) seedU;

    return;
  }


  void updateC2DDelta_MI(int *seed,
			 const int *nTry,
			 int *nAccept,
			 const double *epsilon,
			 double *c2,
			 double *Delta,
			 const int *Q,
			 const int *G,
			 const int *S,
			 const double *x,
			 const int *psi,
			 const double *nu,
			 const int *delta,
			 const double *r,
			 const double *sigma2,
			 const double *phi,
			 const double *tau2R,
			 const double *b,
			 const double *c2Max) {
    unsigned int seedU = (unsigned int) *seed;

    updateC2DDelta(&seedU,*nTry,nAccept,*epsilon,c2,Delta,*Q,*G,S,x,psi,
		   nu,delta,r,sigma2,phi,tau2R,b,*c2Max);

    *seed = (int) seedU;

    return;
  }


  void updateGamma2Nu_MI(int *seed,
			 const int *nTry,
			 int *nAccept,
			 const double *epsilon,
			 double *gamma2,
			 double *nu,
			 const int *Q,
			 const int *G,
			 const int *S,
			 const double *x,
			 const int *psi,
			 const int *delta,
			 const double *Delta,
			 const double *rho,
			 const double *sigma2,
			 const double *phi,
			 const double *tau2Rho,
			 const double *a) {
    unsigned int seedU = (unsigned int) *seed;

    updateGamma2Nu(&seedU,*nTry,nAccept,*epsilon,gamma2,nu,*Q,*G,S,x,psi,
		   delta,Delta,rho,sigma2,phi,tau2Rho,a);

    *seed = (int) seedU;

    return;
  }


  void updateRDDelta_MI(int *seed,
			const int *nTry,
			int *nAccept,
			const double *epsilon,
			double *r,
			double *Delta,
			const int *Q,
			const int *G,
			const int *S,
			const double *x,
			const int *psi,
			const double *nu,
			const int *delta,
			const double *c2,
			const double *sigma2,
			const double *phi,
			const double *tau2R,
			const double *b,
			const double *nuR) {
    unsigned int seedU = (unsigned int) *seed;
    
    updateRDDelta(&seedU,*nTry,nAccept,*epsilon,r,Delta,*Q,*G,S,x,psi,
		  nu,delta,*c2,sigma2,phi,tau2R,b,*nuR);

    *seed = (int) seedU;

    return;
  }


  void updateRC2_MI(int *seed,
		    const int *nTry,
		    int *nAccept,
		    const double *epsilon,
		    double *r,
		    double *c2,
		    const int *Q,
		    const int *G,
		    const int *delta,
		    const double *Delta,
		    const double *sigma2,
		    const double *tau2R,
		    const double *b,
		    const double *nuR,
		    const double *c2Max) {
    unsigned int seedU = (unsigned int) *seed;

    updateRC2(&seedU,*nTry,nAccept,*epsilon,r,c2,*Q,*G,delta,Delta,
		 sigma2,tau2R,b,*nuR,*c2Max);

    *seed = (int) seedU;

    return;
  }


  void updateSigma2_MI(int *seed,
		       const int *nTry,
		       int *nAccept,
		       const double *epsilon,
		       double *sigma2,
		       const int *Q,
		       const int *G,
		       const int *S,
		       const double *x,
		       const int *psi,
		       const double *nu,
		       const int *delta,
		       const double *Delta,
		       const double *c2,
		       const double *gamma2,
		       const double *r,
		       const double *rho,
		       const double *phi,
		       const double *t,
		       const double *l,
		       const double *tau2R,
		       const double *tau2Rho,
		       const double *a,
		       const double *b) {
    unsigned int seedU = (unsigned int) *seed;

    updateSigma2(&seedU,*nTry,nAccept,*epsilon,sigma2,*Q,*G,S,x,psi,nu,
		 delta,Delta,*c2,*gamma2,r,rho,phi,t,l,tau2R,tau2Rho,a,b);

    *seed = (int) seedU;

    return;
  }


  void updateLSigma2_MI(unsigned int *seed,
			const int *nTry,
			int *nAccept,
			const double *epsilon,
			double *l,
			double *sigma2,
			const int *Q,
			const int *G,
			const int *S,
			const double *x,
			const int *psi,
			const double *nu,
			const int *delta,
			const double *Delta,
			const double *c2,
			const double *gamma2,
			const double *r,
			const double *rho,
			const double *phi,
			const double *t,
			const double *tau2R,
			const double *tau2Rho,
			const double *a,
			const double *b) {
    unsigned int seedU = (unsigned int) *seed;
    
    updateLSigma2(&seedU,*nTry,nAccept,*epsilon,l,sigma2,*Q,*G,S,x,psi,nu,
		  delta,Delta,*c2,*gamma2,r,rho,phi,t,tau2R,tau2Rho,a,b);

    *seed = (int) seedU;

    return;
  }



  void updateTSigma2_MI(unsigned int *seed,
			const int *nTry,
			int *nAccept,
			const double *epsilon,
			double *t,
			double *sigma2,
			const int *Q,
			const int *G,
			const int *S,
			const double *x,
			const int *psi,
			const double *nu,
			const int *delta,
			const double *Delta,
			const double *c2,
			const double *gamma2,
			const double *r,
			const double *rho,
			const double *phi,
			const double *l,
			const double *tau2R,
			const double *tau2Rho,
			const double *a,
			const double *b) {
    unsigned int seedU = (unsigned int) *seed;

    updateTSigma2(&seedU,*nTry,nAccept,*epsilon,t,sigma2,*Q,*G,S,x,psi,nu,
		  delta,Delta,*c2,*gamma2,r,rho,phi,l,tau2R,tau2Rho,a,b); // should it be updateTSigma2?

    *seed = (int) seedU;

    return;
  }


  // finished: updates for model I




  // updates for model II


  void updateSigma2_MII(int *seed,
			const int *nTry,
			int *nAccept,
			const double *epsilon,
			double *sigma2,
			const int *Q,
			const int *G,
			const int *S,
			const double *x,
			const int *psi,
			const double *nu,
			const int *delta,
			const double *Delta,
			const double *gamma2,
			const double *r,
			const double *rho,
			const double *phi,
			const double *t,
			const double *l,
			const double *tau2R,
			const double *tau2Rho,
			const double *a,
			const double *b,
			const double *Omega,
			const int *nClique,
			const int *oldClique,
			const int *nOldComponents,
			const int *nNewComponents,
			const int *oldComponents) {
    unsigned int seedU = (unsigned int) *seed;

    
    vector<vector<vector<double> > > OmegaTransformed;
    vector<int> oldCliqueTransformed;
    vector<vector<int> > oldComponentsTransformed;
    transformGraph(nClique,oldClique,nOldComponents,oldComponents,oldCliqueTransformed,oldComponentsTransformed);
    transformOmega(nClique,nOldComponents,nNewComponents,Omega,OmegaTransformed);
    
    updateSigma2_HyperInverseWishart(&seedU,*nTry,nAccept,*epsilon,sigma2,*Q,*G,S,x,psi,nu,
				     delta,Delta,*gamma2,r,rho,phi,t,l,tau2R,tau2Rho,a,b,
				     OmegaTransformed,oldCliqueTransformed,oldComponentsTransformed);

    *seed = (int) seedU;

    return;
  }

  
  void updateLSigma2_MII(unsigned int *seed,
			 const int *nTry,
			 int *nAccept,
			 const double *epsilon,
			 double *l,
			 double *sigma2,
			 const int *Q,
			 const int *G,
			 const int *S,
			 const double *x,
			 const int *psi,
			 const double *nu,
			 const int *delta,
			 const double *Delta,
			 const double *gamma2,
			 const double *r,
			 const double *rho,
			 const double *phi,
			 const double *t,
			 const double *tau2R,
			 const double *tau2Rho,
			 const double *a,
			 const double *b,
			 const double *Omega,
			 const int *nClique,
			 const int *oldClique,
			 const int *nOldComponents,
			 const int *nNewComponents,
			 const int *oldComponents) {
    unsigned int seedU = (unsigned int) *seed;

    vector<vector<vector<double> > > OmegaTransformed;
    vector<int> oldCliqueTransformed;
    vector<vector<int> > oldComponentsTransformed;
    transformGraph(nClique,oldClique,nOldComponents,oldComponents,oldCliqueTransformed,oldComponentsTransformed);
    transformOmega(nClique,nOldComponents,nNewComponents,Omega,OmegaTransformed);
    
    updateLSigma2_HyperInverseWishart(&seedU,*nTry,nAccept,*epsilon,l,sigma2,*Q,*G,S,x,psi,nu,
				      delta,Delta,*gamma2,r,rho,phi,t,tau2R,tau2Rho,a,b,
				      OmegaTransformed,oldCliqueTransformed,oldComponentsTransformed);
    
    *seed = (int) seedU;
    
    return;
  }


  void updateTSigma2_MII(unsigned int *seed,
			 const int *nTry,
			 int *nAccept,
			 const double *epsilon,
			 double *t,
			 double *sigma2,
			 const int *Q,
			 const int *G,
			 const int *S,
			 const double *x,
			 const int *psi,
			 const double *nu,
			 const int *delta,
			 const double *Delta,
			 const double *gamma2,
			 const double *r,
			 const double *rho,
			 const double *phi,
			 const double *l,
			 const double *tau2R,
			 const double *tau2Rho,
			 const double *a,
			 const double *b,
			 const double *Omega,
			 const int *nClique,
			 const int *oldClique,
			 const int *nOldComponents,
			 const int *nNewComponents,
			 const int *oldComponents) {
    unsigned int seedU = (unsigned int) *seed;
    
    vector<vector<vector<double> > > OmegaTransformed;
    vector<int> oldCliqueTransformed;
    vector<vector<int> > oldComponentsTransformed;
    transformGraph(nClique,oldClique,nOldComponents,oldComponents,oldCliqueTransformed,oldComponentsTransformed);
    transformOmega(nClique,nOldComponents,nNewComponents,Omega,OmegaTransformed);
    
    updateTSigma2_HyperInverseWishart(&seedU,*nTry,nAccept,*epsilon,t,sigma2,*Q,*G,S,x,psi,nu,
				      delta,Delta,*gamma2,r,rho,phi,l,tau2R,tau2Rho,a,b,
				      OmegaTransformed,oldCliqueTransformed,oldComponentsTransformed);
    
    *seed = (int) seedU;
    
    return;
  }



  void updateOmega_MII(unsigned int *seed,
		       int *nAccept,
		       double *Omega,
		       const int *Q,
		       const int *G,
		       const double *Delta,
		       const double *r,
		       const double *sigma2,
		       const double *tau2R,
		       const double *b,
		       const double *zeta,
		       const double *D,
		       const int *nClique,
		       const int *oldClique,
		       const int *nOldComponents,
		       const int *nNewComponents,
		       const int *oldComponents) {
    unsigned int seedU = (unsigned int) *seed;

		       
    vector<vector<vector<double> > > OmegaTransformed;
    vector<vector<vector<double> > > DTransformed;
    vector<int> oldCliqueTransformed;
    vector<vector<int> > oldComponentsTransformed;
    transformGraph(nClique,oldClique,nOldComponents,oldComponents,oldCliqueTransformed,oldComponentsTransformed);
    transformOmega(nClique,nOldComponents,nNewComponents,D,DTransformed);
    transformOmega(nClique,nOldComponents,nNewComponents,Omega,OmegaTransformed);

    updateOmega_HyperInverseWishart(&seedU,nAccept,OmegaTransformed,*Q,*G,Delta,r,sigma2,
				    tau2R,b,*zeta,DTransformed,oldCliqueTransformed,
				    oldComponentsTransformed);

    inverseTransformOmega(OmegaTransformed,Omega);
		       
    *seed = (int) seedU;
  }


  void updateDeltaStar_MII(unsigned int *seed,
			   int *nAccept,
			   double *Delta,
			   const int *Q,
			   const int *G,
			   const int *S,
			   const double *x,
			   const int *psi,
			   const double *nu,
			   const int *delta,
			   const double *r,
			   const double *sigma2,
			   const double *phi,
			   const double *tau2R,
			   const double *b,
			   const double *Omega,
			   const int *nClique,
			   const int *oldClique,
			   const int *nOldComponents,
			   const int *nNewComponents,
			   const int *oldComponents) {
    unsigned int seedU = (unsigned int) *seed;
    
    
    vector<vector<vector<double> > > OmegaTransformed;
    vector<int> oldCliqueTransformed;
    vector<vector<int> > oldComponentsTransformed;
    transformGraph(nClique,oldClique,nOldComponents,oldComponents,oldCliqueTransformed,oldComponentsTransformed);
    transformOmega(nClique,nOldComponents,nNewComponents,Omega,OmegaTransformed);
    
    updateDDeltaStar_HyperInverseWishart(&seedU,nAccept,Delta,*Q,*G,S,x,psi,nu,delta,r,sigma2,
					 phi,tau2R,b,OmegaTransformed,oldCliqueTransformed,
					 oldComponentsTransformed);
    
    *seed = (int) seedU;
  }
  

  void updateTau2RDDeltaStar_MII(unsigned int *seed,
				 int *nTry,
				 int *nAccept,
				 double *epsilon,
				 double *tau2R,
				 double *Delta,
				 const int *Q,
				 const int *G,
				 const int *S,
				 const double *x,
				 const int *psi,
				 const double *nu,
				 const int *delta,
				 const double *r,
				 const double *sigma2,
				 const double *phi,
				 const double *b,
				 const double *Omega,
				 const int *nClique,
				 const int *oldClique,
				 const int *nOldComponents,
				 const int *nNewComponents,
				 const int *oldComponents) {
    unsigned int seedU = (unsigned int) *seed;
    
    
    vector<vector<vector<double> > > OmegaTransformed;
    vector<int> oldCliqueTransformed;
    vector<vector<int> > oldComponentsTransformed;
    transformGraph(nClique,oldClique,nOldComponents,oldComponents,oldCliqueTransformed,oldComponentsTransformed);
    transformOmega(nClique,nOldComponents,nNewComponents,Omega,OmegaTransformed);
    
    updateTau2RDDeltaStar_HyperInverseWishart(&seedU,*nTry,nAccept,*epsilon,tau2R,Delta,*Q,
					      *G,S,x,psi,nu,delta,r,sigma2,phi,b,
					      OmegaTransformed,oldCliqueTransformed,
					      oldComponentsTransformed);
    
    *seed = (int) seedU;
  }


  void updateBDDeltaStar_MII(unsigned int *seed,
			     int *nTry,
			     int *nAccept,
			     double *epsilon,
			     double *b,
			     double *Delta,
			     const int *Q,
			     const int *G,
			     const int *S,
			     const double *x,
			     const int *psi,
			     const double *nu,
			     const int *delta,
			     const double *r,
			     const double *sigma2,
			     const double *phi,
			     const double *tau2R,
			     const double *pB0,
			     const double *pB1,
			     const double *alphaB,
			     const double *betaB,
			     const double *Omega,
			     const int *nClique,
			     const int *oldClique,
			     const int *nOldComponents,
			     const int *nNewComponents,
			     const int *oldComponents) {
    unsigned int seedU = (unsigned int) *seed;
    
    
    vector<vector<vector<double> > > OmegaTransformed;
    vector<int> oldCliqueTransformed;
    vector<vector<int> > oldComponentsTransformed;
    transformGraph(nClique,oldClique,nOldComponents,oldComponents,oldCliqueTransformed,oldComponentsTransformed);
    transformOmega(nClique,nOldComponents,nNewComponents,Omega,OmegaTransformed);
    
    updateBDDeltaStar_HyperInverseWishart(&seedU,*nTry,nAccept,*epsilon,b,Delta,*Q,
					  *G,S,x,psi,nu,delta,r,sigma2,phi,tau2R,
					  *pB0,*pB1,*alphaB,*betaB,
					  OmegaTransformed,oldCliqueTransformed,
					  oldComponentsTransformed);
    
    *seed = (int) seedU;
  }


  void updateRDDeltaStar_MII(unsigned int *seed,
			     int *nTry,
			     int *nAccept,
			     double *epsilon,
			     double *r,
			     double *Delta,
			     const int *Q,
			     const int *G,
			     const int *S,
			     const double *x,
			     const int *psi,
			     const double *nu,
			     const int *delta,
			     const double *sigma2,
			     const double *phi,
			     const double *tau2R,
			     const double *b,
			     const double *nuR,
			     const double *Omega,
			     const int *nClique,
			     const int *oldClique,
			     const int *nOldComponents,
			     const int *nNewComponents,
			     const int *oldComponents) {
    unsigned int seedU = (unsigned int) *seed;
    
    
    vector<vector<vector<double> > > OmegaTransformed;
    vector<int> oldCliqueTransformed;
    vector<vector<int> > oldComponentsTransformed;
    transformGraph(nClique,oldClique,nOldComponents,oldComponents,oldCliqueTransformed,oldComponentsTransformed);
    transformOmega(nClique,nOldComponents,nNewComponents,Omega,OmegaTransformed);
    
    updateRDDeltaStar_HyperInverseWishart(&seedU,*nTry,nAccept,*epsilon,r,Delta,*Q,
					  *G,S,x,psi,nu,delta,sigma2,phi,tau2R,b,*nuR,
					  OmegaTransformed,oldCliqueTransformed,
					  oldComponentsTransformed);
    
    *seed = (int) seedU;
  }


  // finished: updates for model II



  // updates only for model AI


  void updateDeltaDDelta_MAI(int *seed,
			     const int *nTry,
			     int *nAccept,
			     int *delta,
			     double *Delta,
			     const int *Q,
			     const int *G,
			     const int *S,
			     const double *x,
			     const int *psi,
			     const double *nu,
			     const double *c2,
			     const double *r,
			     const double *sigma2,
			     const double *phi,
			     const double *tau2R,
			     const double *xi,
			     const double *b) {
    unsigned int seedU = (unsigned int) *seed;

    updateDeltaDDelta(&seedU,*nTry,nAccept,delta,Delta,*Q,*G,S,x,psi,
		      nu,*c2,r,sigma2,phi,tau2R,xi,b);

    *seed = (int) seedU;

    return;
  }

  // finished: updates for model AI
  


  // updates only for model BI


  void updateDeltaDDelta_MBI(int *seed,
			     const int *nTry,
			     int *nAccept,
			     int *delta,
			     double *Delta,
			     const int *Q,
			     const int *G,
			     const int *S,
			     const double *x,
			     const int *psi,
			     const double *nu,
			     const double *c2,
			     const double *r,
			     const double *sigma2,
			     const double *phi,
			     const double *tau2R,
			     const double *xi,
			     const double *b) {
    unsigned int seedU = (unsigned int) *seed;
    
    updateDeltaDDelta_onedelta(&seedU,*nTry,nAccept,delta,Delta,*Q,*G,S,x,psi,
			       nu,*c2,r,sigma2,phi,tau2R,xi,b);

    *seed = (int) seedU;

    return;
  }


  // finished: updates for model BI


  // updates only for model CI


  void updateDeltaDDelta_MCI(unsigned int *seed,
			     const int *nTry,
			     int *nAccept,
			     int *delta,
			     double *Delta,
			     const int *Q,
			     const int *G,
			     const int *S,
			     const double *x,
			     const int *psi,
			     const double *nu,
			     const double *c2,
			     const double *r,
			     const double *sigma2,
			     const double *phi,
			     const double *tau2R,
			     const double *b,
			     const int *nNeighbour,
			     const int *neighbour,
			     const double *alpha,
			     const double *beta,
			     const double *betag) {
    unsigned int seedU = (unsigned int) *seed;

    vector<vector<int> > nn;
    nn.resize(*G);
    int g;
    for (g = 0; g < *G; g++) nn[g].resize(0);

    int i;
    for (i = 0; i < *nNeighbour; i++) {
      int g1 = neighbour[2 * i];
      int g2 = neighbour[2 * i + 1];
      nn[g1].push_back(g2);
      nn[g2].push_back(g1);
    }

    updateDeltaDDelta_MRF2(&seedU,*nTry,nAccept,delta,Delta,
			   *Q,*G,S,x,psi,nu,*c2,r,sigma2,phi,
			   tau2R,b,nn,*alpha,*beta,*betag);

    *seed = (int) seedU;

    return;
  }


  // finished: update for model CI


  // updates only for model DI


  void updateDeltaDDelta_MDI(unsigned int *seed,
			     const int *nTry,
			     int *nAccept,
			     int *delta,
			     double *Delta,
			     const int *Q,
			     const int *G,
			     const int *S,
			     const double *x,
			     const int *psi,
			     const double *nu,
			     const double *c2,
			     const double *r,
			     const double *sigma2,
			     const double *phi,
			     const double *tau2R,
			     const double *b,
			     const int *nNeighbour,
			     const int *neighbour,
			     const double *alpha,
			     const double *beta) {
    unsigned int seedU = (unsigned int) *seed;
    
    vector<vector<int> > nn;
    nn.resize(*G);
    int g;
    for (g = 0; g < *G; g++) nn[g].resize(0);
    
    int i;
    for (i = 0; i < *nNeighbour; i++) {
      int g1 = neighbour[2 * i];
      int g2 = neighbour[2 * i + 1];
      nn[g1].push_back(g2);
      nn[g2].push_back(g1);
    }
    
    updateDeltaDDelta_MRF2_onedelta(&seedU,*nTry,nAccept,delta,Delta,
				    *Q,*G,S,x,psi,nu,*c2,r,sigma2,phi,
				    tau2R,b,nn,*alpha,*beta);
    
    *seed = (int) seedU;
    
    return;
  }

  // finished: update for model DI


  // updates only for model AII


  void updateDelta_MAII(unsigned int *seed,
			const int *nTry,
			int *nAccept,
			int *delta,
			const int *Q,
			const int *G,
			const int *S,
			const double *x,
			const int *psi,
			const double *nu,
			const double *Delta,
			const double *r,
			const double *sigma2,
			const double *phi,
			const double *xi,
			const double *b) {
    unsigned int seedU = (unsigned int) *seed;
    
    updateDelta_HyperInverseWishart(&seedU,*nTry,nAccept,delta,*Q,*G,S,x,psi,
				    nu,Delta,r,sigma2,phi,xi,b);
    
    *seed = (int) seedU;
    
    return;
  }
  
  // finished: update for model AII  
  



  // updates only for model BII

  void updateDelta_MBII(unsigned int *seed,
			const int *nTry,
			int *nAccept,
			int *delta,
			const int *Q,
			const int *G,
			const int *S,
			const double *x,
			const int *psi,
			const double *nu,
			const double *Delta,
			const double *r,
			const double *sigma2,
			const double *phi,
			const double *xi,
			const double *b) {
    unsigned int seedU = (unsigned int) *seed;
    
    updateDelta_HyperInverseWishart_onedelta(&seedU,*nTry,nAccept,delta,*Q,*G,S,x,psi,
					     nu,Delta,r,sigma2,phi,xi,b);
    
    *seed = (int) seedU;
    
    return;
  }
  
  // finished: update for model BII  
  


  // updates only for model CII

  void updateDelta_MCII(unsigned int *seed,
			const int *nTry,
			int *nAccept,
			int *delta,
			const int *Q,
			const int *G,
			const int *S,
			const double *x,
			const int *psi,
			const double *nu,
			const double *Delta,
			const double *r,
			const double *sigma2,
			const double *phi,
			const double *b,
			const int *nNeighbour,
			const int *neighbour,
			const double *alpha,
			const double *beta,
			const double *betag) {
    unsigned int seedU = (unsigned int) *seed;
    
    vector<vector<int> > nn;
    nn.resize(*G);
    int g;
    for (g = 0; g < *G; g++) nn[g].resize(0);
    
    int i;
    for (i = 0; i < *nNeighbour; i++) {
      int g1 = neighbour[2 * i];
      int g2 = neighbour[2 * i + 1];
      nn[g1].push_back(g2);
      nn[g2].push_back(g1);
    }
    
    updateDelta_HyperInverseWishart_MRF2(&seedU,*nTry,nAccept,delta,*Q,*G,S,x,psi,
					 nu,Delta,r,sigma2,phi,b,nn,*alpha,
					 *beta,*betag);
    
    *seed = (int) seedU;
    
    return;
  }
  
  // finished: update for model CII  
  


  // updates only for model DII

  void updateDelta_MDII(unsigned int *seed,
			const int *nTry,
			int *nAccept,
			int *delta,
			const int *Q,
			const int *G,
			const int *S,
			const double *x,
			const int *psi,
			const double *nu,
			const double *Delta,
			const double *r,
			const double *sigma2,
			const double *phi,
			const double *b,
			const int *nNeighbour,
			const int *neighbour,
			const double *alpha,
			const double *beta) {
    unsigned int seedU = (unsigned int) *seed;
    
    vector<vector<int> > nn;
    nn.resize(*G);
    int g;
    for (g = 0; g < *G; g++) nn[g].resize(0);
    
    int i;
    for (i = 0; i < *nNeighbour; i++) {
      int g1 = neighbour[2 * i];
      int g2 = neighbour[2 * i + 1];
      nn[g1].push_back(g2);
      nn[g2].push_back(g1);
    }
    
    updateDelta_HyperInverseWishart_MRF2_onedelta(&seedU,*nTry,nAccept,delta,*Q,*G,S,x,psi,
						  nu,Delta,r,sigma2,phi,b,nn,*alpha,*beta);
    
    *seed = (int) seedU;
    
    return;
  }
  
  // finished: update for model DII  



} // extern "C"
