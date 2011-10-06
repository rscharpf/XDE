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
  // diagonal elements
  for (q = 0; q < Q; q++) {
    int kqg = qg2index(q,g,Q,G);
    Sigma[q].resize(Q);
    Sigma[q][q] = gamma2 * tau2[q];
    Sigma[q][q] *= exp(a[q] * log(sigma2[kqg]));
  }

  int q1,q2;
  // off-diagonal elements
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



void updateMRF1perfect_onedelta(int g,vector<int> &valueLower,
				vector<int> &valueUpper,
				const vector<double> &potOn,
				const vector<double> &potOff,
				const vector<vector<int> > &neighbour,
				double eta0,double omega0,double kappa,
				Random &ran) {
  double potLower = potOff[g] - potOn[g];
  double potUpper = potOff[g] - potOn[g];

  // add potential for clique centered at gene g

  int n = neighbour[g].size();
  double omega;
  if (n > 0)
    omega = omega0 * ((double) n) / (kappa + ((double) n));
  else
    omega = 0.0;

  double meanLower = 0.0;
  double meanUpper = 0.0;
  int gg;
  for (gg = 0; gg < neighbour[g].size(); gg++) {
    meanLower += (double) valueLower[neighbour[g][gg]];
    meanUpper += (double) valueUpper[neighbour[g][gg]];
  }

  if (neighbour[g].size() > 0) {
    meanLower /= (double) neighbour[g].size();
    meanUpper /= (double) neighbour[g].size();

    meanLower = (1.0 - omega) * eta0 + omega * meanLower;
    meanUpper = (1.0 - omega) * eta0 + omega * meanUpper;
  }
  else {
    meanLower = eta0;
    meanUpper = eta0;
  }

  potUpper += - log(1.0 - meanLower) + log(meanLower);
  potLower += - log(1.0 - meanUpper) + log(meanUpper);

  // add potential for cliques centered at each gene connected to g

  for (gg = 0; gg < neighbour[g].size(); gg++) {
    int gene = neighbour[g][gg];
    int n = neighbour[gene].size();
    double omega;
    if (n > 0)
      omega = omega0 * ((double) n) / (kappa + ((double) n));
    else
      omega = 0.0;

    double meanLower = 0.0;
    double meanUpper = 0.0;
    int ggg;
    for (ggg = 0; ggg < neighbour[gene].size(); ggg++) {
      if (neighbour[gene][ggg] != g) {
	meanLower += (double) valueLower[neighbour[gene][ggg]];
	meanUpper += (double) valueUpper[neighbour[gene][ggg]];
      }
    }

    meanLower /= (double) neighbour[gene].size();
    meanUpper /= (double) neighbour[gene].size();

    meanLower = (1.0 - omega) * eta0 + omega * meanLower;
    meanUpper = (1.0 - omega) * eta0 + omega * meanUpper;

    double extra = omega / ((double) neighbour[gene].size());

    if (valueLower[gene] == 0 && valueUpper[gene] == 0) {
      potUpper += - log(1.0 - meanUpper) + log(1.0 - meanUpper - extra);
      potLower += - log(1.0 - meanLower) + log(1.0 - meanLower - extra);
    }
    else if (valueLower[gene] == 1 && valueUpper[gene] == 1) {
      potUpper += - log(meanUpper) + log(meanUpper + extra);
      potLower += - log(meanLower) + log(meanLower + extra);
    }
    else {
      double p0Upper = - log(1.0 - meanUpper) + log(1.0 - meanUpper - extra);
      double p0Lower = - log(1.0 - meanLower) + log(1.0 - meanLower - extra);

      double p1Upper = - log(meanUpper) + log(meanUpper + extra);
      double p1Lower = - log(meanLower) + log(meanLower + extra);

      if (p0Lower < p1Lower)
	potLower += p1Lower;
      else
	potLower += p0Lower;

      if (p0Upper < p1Upper)
	potUpper += p0Upper;
      else
	potUpper += p1Upper;
    }
  }

  double probLower;
  if (potUpper > 0.0)
    probLower = 1.0 / (1.0 + exp(- potUpper));
  else
    probLower = exp(potUpper) / (1.0 + exp(potUpper));

  double probUpper;
  if (potLower > 0.0)
    probUpper = 1.0 / (1.0 + exp(- potLower));
  else
    probUpper = exp(potLower) / (1.0 + exp(potLower));

  double u = ran.Unif01();
  if (u < probLower)
    valueLower[g] = 1;
  else
    valueLower[g] = 0;

  if (u < probUpper)
    valueUpper[g] = 1;
  else
    valueUpper[g] = 0;

  return;
}




double perfectMRF1_onedelta(int *delta,int G,
			    const vector<vector<int> > &neighbour,
			    const vector<double> &potOn,
			    const vector<double> &potOff,
			    double eta0,double omega0,double kappa,
			    unsigned int *seed,int draw) {
  unsigned int finalStart = *seed;

  if (draw == 1) {
    vector<int> start(1,-1);
    vector<unsigned int> seeds(1,finalStart);

    unsigned int nextSeed;
    int finished = 0;
    while (finished == 0) {
      vector<int> valueLower(G,0);
      vector<int> valueUpper(G,1);

      int b;
      for (b = start.size() - 1; b >= 0; b--) {
	int first = start[b];
	int last;
	if (b > 0)
	  last = start[b-1];
	else
	  last = 0;

	Random ran(seeds[b]);
	int k;
	for (k = first; k < last; k++) {
	  int g;
	  for (g = 0; g < G; g++)
	    updateMRF1perfect_onedelta(g,valueLower,valueUpper,potOn,
				       potOff,neighbour,
				       eta0,omega0,kappa,ran);
	}

	unsigned int dummy = 1;
	if (b == start.size() - 1) nextSeed = ran.ChangeSeed(dummy);
      }

      int nUndef = 0;
      int g;
      for (g = 0; g < G; g++)
	nUndef += (valueLower[g] != valueUpper[g]);
      cout << "nUndef: " << nUndef << endl;

      if (nUndef == 0) {
	finished = 1;
	finalStart = nextSeed;
      }
      else {
	finished = 0;
	seeds.push_back(nextSeed);
	start.push_back(2*start[start.size() - 1]);
      }

      if (finished == 1) {
	for (g = 0; g < G; g++)
	  delta[g] = valueLower[g];
      }
    }

    *seed = nextSeed;
  }

  double pot = 0.0;
  int g;
  for (g = 0; g < G; g++) {
    if (delta[g] == 1)
      pot += potOn[g];
    else
      pot += potOff[g];

    int n = neighbour[g].size();
    double omega;
    if (n > 0)
      omega = omega0 * ((double) n) / (kappa + ((double) n));
    else
      omega = 0.0;
    int nOn = 0;
    int gg;
    for (gg = 0; gg < neighbour[g].size(); gg++)
      nOn += delta[neighbour[g][gg]];
    double fraction = ((double) nOn) / ((double) neighbour[g].size());
    double eta;
    if (omega > 0)
      eta = (1.0 - omega) * eta0 + omega * fraction;
    else
      eta = eta0;

    if (delta[g] == 1)
      pot += - log(eta);
    else
      pot += - log(1.0 - eta);
  }

  return pot;
}





void updateMRF2perfect_onedelta(int g,vector<int> &valueLower,
				vector<int> &valueUpper,
				const vector<double> &potOn,
				const vector<double> &potOff,
				const vector<vector<int> > &neighbour,
				double alpha,double beta,Random &ran) {
  double potLower = potOff[g] - potOn[g];
  double potUpper = potOff[g] - potOn[g];

  potLower += - alpha;
  potUpper += - alpha;

  int k;
  for (k = 0; k < neighbour[g].size(); k++) {
    int gg = neighbour[g][k];
    int ng = neighbour[g].size();
    int ngg = neighbour[gg].size();
    //    double w = beta * exp(- kappa * log((double) (ng * ngg)));
    double w = beta * (1.0 / ((double) ng) + 1.0 / ((double) ngg));

    if (valueLower[gg] == 0 && valueUpper[gg] == 0) {
      potLower += w;
      potUpper += w;
    }
    else if (valueLower[gg] == 1 && valueUpper[gg] == 1) {
      potLower += - w;
      potUpper += - w;
    }
    else {
      potLower += w;
      potUpper += - w;
    }
  }


  double probLower;
  double probUpper;
  if (potLower < 0.0)
    probLower = 1.0 / (1.0 + exp(potLower));
  else
    probLower = exp(- potLower) / (1.0 + exp(- potLower));

  if (potUpper < 0.0)
    probUpper = 1.0 / (1.0 + exp(potUpper));
  else
    probUpper = exp(- potUpper) / (1.0 + exp(- potUpper));

  double u = ran.Unif01();
  if (u < probLower)
    valueLower[g] = 1;
  else
    valueLower[g] = 0;

  if (u < probUpper)
    valueUpper[g] = 1;
  else
    valueUpper[g] = 0;

  return;
}




double perfectMRF2_onedelta(int *delta,int G,
			    const vector<vector<int> > &neighbour,
			    const vector<double> &potOn,
			    const vector<double> &potOff,
			    double alpha,double beta,
			    unsigned int *seed,int draw) {
  unsigned int finalStart = *seed;

  if (draw == 1) {
    vector<int> start(1,-1);
    vector<unsigned int> seeds(1,finalStart);

    unsigned int nextSeed;
    int finished = 0;
    while (finished == 0) {
      vector<int> valueLower(G,0);
      vector<int> valueUpper(G,1);

      int b;
      for (b = start.size() - 1; b >= 0; b--) {
	int first = start[b];
	int last;
	if (b > 0)
	  last = start[b-1];
	else
	  last = 0;

	Random ran(seeds[b]);
	int k;
	for (k = first; k < last; k++) {
	  int g;
	  for (g = 0; g < G; g++)
	    updateMRF2perfect_onedelta(g,valueLower,valueUpper,potOn,
				       potOff,neighbour,alpha,beta,ran);
	}

	unsigned int dummy = 1;
	if (b == start.size() - 1) nextSeed = ran.ChangeSeed(dummy);
      }

      int nUndef = 0;
      int g;
      for (g = 0; g < G; g++)
	nUndef += (valueLower[g] != valueUpper[g]);
      cout << "nUndef: " << nUndef << endl;

      if (nUndef == 0) {
	finished = 1;
	finalStart = nextSeed;
      }
      else {
	finished = 0;
	seeds.push_back(nextSeed);
	start.push_back(2*start[start.size() - 1]);
      }

      if (finished == 1) {
	for (g = 0; g < G; g++)
	  delta[g] = valueLower[g];
      }
    }

    *seed = nextSeed;
  }

  double pot = 0.0;
  int g;
  for (g = 0; g < G; g++) {
    if (delta[g] == 1)
      pot += - alpha + potOn[g];
    else
      pot += potOff[g];

    int k;
    for (k = 0; k < neighbour[g].size(); k++) {
      int gg = neighbour[g][k];
      if (delta[g] == delta[gg]) {
	int ng = neighbour[g].size();
	int ngg = neighbour[gg].size();
	//	double w = 0.5 * exp(- kappa * log((double) (ng * ngg)));
	double w = 1.0 / ((double) ng);

	pot += - beta * w;
      }
    }
  }

  return pot;
}




void updateMRF2perfect(int q,int g,int Q,int G,vector<int> &valueLower,
		       vector<int> &valueUpper,const vector<double> &potOn,
		       const vector<double> &potOff,
		       const vector<vector<int> > &neighbour,
		       double alpha,double beta,double betag,
		       Random &ran) {
  int kqg = qg2index(q,g,Q,G);
  double potLower = potOff[kqg] - potOn[kqg];
  double potUpper = potOff[kqg] - potOn[kqg];

  potLower += - alpha;
  potUpper += - alpha;

  int k;
  for (k = 0; k < neighbour[g].size(); k++) {
    int gg = neighbour[g][k];
    int ng = neighbour[g].size();
    int ngg = neighbour[gg].size();
    //    double w = beta * exp(- kappa * log((double) (ng * ngg)));
    double w = beta * (1.0 / ((double) ng) + 1.0 / ((double) ngg));


    int kqgg = qg2index(q,gg,Q,G);
    if (valueLower[kqgg] == 0 && valueUpper[kqgg] == 0) {
      potLower += w;
      potUpper += w;
    }
    else if (valueLower[kqgg] == 1 && valueUpper[kqgg] == 1) {
      potLower += - w;
      potUpper += - w;
    }
    else {
      potLower += w;
      potUpper += - w;
    }
  }

  int qq;
  for (qq = 0; qq < Q; qq++) {
    if (qq != q) {
      int kqqg = qg2index(qq,g,Q,G);

      if (valueLower[kqqg] == 0 && valueUpper[kqqg] == 0) {
	potLower += betag / ((double) (Q - 1));
	potUpper += betag / ((double) (Q - 1));
      }
      else if (valueLower[kqqg] == 1 && valueUpper[kqqg] == 1) {
	potLower += - betag / ((double) (Q - 1));
	potUpper += - betag / ((double) (Q - 1));
      }
      else {
	potLower += betag / ((double) (Q - 1));
	potUpper -= betag / ((double) (Q - 1));
      }
    }
  }

  double probLower;
  double probUpper;
  if (potLower < 0.0)
    probLower = 1.0 / (1.0 + exp(potLower));
  else
    probLower = exp(- potLower) / (1.0 + exp(- potLower));

  if (potUpper < 0.0)
    probUpper = 1.0 / (1.0 + exp(potUpper));
  else
    probUpper = exp(- potUpper) / (1.0 + exp(- potUpper));

  kqg = qg2index(q,g,Q,G);
  double u = ran.Unif01();
  if (u < probLower)
    valueLower[kqg] = 1;
  else
    valueLower[kqg] = 0;

  if (u < probUpper)
    valueUpper[kqg] = 1;
  else
    valueUpper[kqg] = 0;

  return;
}




double perfectMRF2(int *delta,int Q,int G,
		   const vector<vector<int> > &neighbour,
		   const vector<double> &potOn,
		   const vector<double> &potOff,
		   double alpha,double beta,
		   double betag,unsigned int *seed,int draw) {
  unsigned int finalStart = *seed;

  if (draw == 1) {
    vector<int> start(1,-1);
    vector<unsigned int> seeds(1,finalStart);

    unsigned int nextSeed;
    int finished = 0;
    while (finished == 0) {
      vector<int> valueLower(Q * G,0);
      vector<int> valueUpper(Q * G,1);

      int b;
      for (b = start.size() - 1; b >= 0; b--) {
	int first = start[b];
	int last;
	if (b > 0)
	  last = start[b-1];
	else
	  last = 0;

	Random ran(seeds[b]);
	int k;
	for (k = first; k < last; k++) {
	  int q,g;
	  for (q = 0; q < Q; q++)
	    for (g = 0; g < G; g++)
	      updateMRF2perfect(q,g,Q,G,valueLower,valueUpper,potOn,potOff,
				neighbour,alpha,beta,betag,ran);
	}

	unsigned int dummy = 1;
	if (b == start.size() - 1) nextSeed = ran.ChangeSeed(dummy);
      }

      int nUndef = 0;
      int q,g;
      for (q = 0; q < Q; q++)
	for (g = 0; g < G; g++) {
	  int kqg = qg2index(q,g,Q,G);
	  nUndef += (valueLower[kqg] != valueUpper[kqg]);
	}
      cout << "nUndef: " << nUndef << endl;

      if (nUndef == 0) {
	finished = 1;
	finalStart = nextSeed;
      }
      else {
	finished = 0;
	seeds.push_back(nextSeed);
	start.push_back(2*start[start.size() - 1]);
      }

      if (finished == 1) {
	for (q = 0; q < Q; q++)
	  for (g = 0; g < G; g++) {
	    int kqg = qg2index(q,g,Q,G);
	    delta[kqg] = valueLower[kqg];
	  }
      }
    }

    *seed = nextSeed;
  }

  double pot = 0.0;
  int q,g;
  for (q = 0; q < Q; q++)
    for (g = 0; g < G; g++) {
      int kqg = qg2index(q,g,Q,G);
      if (delta[kqg] == 1)
	pot += - alpha + potOn[kqg];
      else
	pot += potOff[kqg];

      int k;
      for (k = 0; k < neighbour[g].size(); k++) {
	int gg = neighbour[g][k];
	int kqgg = qg2index(q,gg,Q,G);
	if (delta[kqg] == delta[kqgg]) {
	  int ng = neighbour[g].size();
	  int ngg = neighbour[gg].size();
	  //	  double w = 0.5 * exp(- kappa * log((double) (ng * ngg)));
	  double w = 1.0 / ((double) ng);

	  pot += - beta * w;
	}
      }
    }

  int qq;
  for (q = 0; q < Q; q++)
    for (qq = q + 1; qq < Q; qq++)
      for (g = 0; g < G; g++) {
	int kqg = qg2index(q,g,Q,G);
	int kqqg = qg2index(qq,g,Q,G);
	if (delta[kqg] == delta[kqqg])
	  pot += - betag / ((double) (Q - 1));
      }

  return pot;
}



