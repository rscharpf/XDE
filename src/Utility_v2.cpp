#include <iostream>
#include <fstream>
#include <stdio.h>
#include <map>

#include "Random_v2.h"
#include "Matrix_v2.h"
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
      //      cout << "nUndef: " << nUndef << endl;

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
      //      cout << "nUndef: " << nUndef << endl;

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
      //      cout << "nUndef: " << nUndef << endl;

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



double OmegaGibbs(double df,const vector<vector<vector<double> > > &D,
		  const vector<int> &oldClique,
		  const vector<vector<int> > &oldComponents,
		  int Q,int G,const double *Delta,
		  const double *r,const double *sigma2,
		  const double *tau2R,const double *b,
		  vector<vector<vector<double> > > &Omega,
		  Random &ran,int draw) {
  double dfNew = df + ((double) Q);

  // construct R matrix

  vector<vector<double> > R;
  R.resize(Q);
  int p,q;
  for (p = 0; p < Q; p++) {
    R[p].resize(Q);
  }
  for (p = 0; p < Q; p++) {
    R[p][p] = tau2R[p];
    for (q = p + 1; q < Q; q++) {
      R[p][q] = sqrt(tau2R[p] * tau2R[q]) * r[qq2index(p,q,Q)];
      R[q][p] = R[p][q];
    }
  }
  vector<vector<double> > RInverse;

  inverse(R,RInverse);

  // compute updated D parameter

  vector<vector<vector<double> > > DNew(D);
  vector<vector<vector<double> > > DeltaStar;
  DeltaStar.resize(DNew.size());

  // update first component of D

  int cliqueSize = DNew[0].size();
  DeltaStar[0].resize(cliqueSize);
  int g;
  for (g = 0; g < cliqueSize; g++) {
    DeltaStar[0][g].resize(Q);
    for (q = 0; q < Q; q++)
      DeltaStar[0][g][q] = Delta[qg2index(q,g,Q,G)] / exp(0.5 * b[q] * log(sigma2[qg2index(q,g,Q,G)]));
  }

  vector<vector<double> > temp;
  quadratic2(DeltaStar[0],RInverse,temp);
  int g1,g2;
  for (g1 = 0; g1 < DNew[0].size(); g1++)
    for (g2 = 0; g2 < DNew[0][g1].size(); g2++)
      DNew[0][g1][g2] += temp[g1][g2];

  int first = cliqueSize;

  // update remaining components of D

  int k;
  for (k = 1; k < DNew.size(); k++) {
    int cliqueSize = DNew[k].size();
    DeltaStar[k].resize(cliqueSize);

    int g;
    for (g = 0; g < oldComponents[k].size(); g++) {
      DeltaStar[k][g].resize(Q);
      for (q = 0; q < Q; q++)
	DeltaStar[k][g][q] = DeltaStar[oldClique[k]][oldComponents[k][g]][q];
    }

    for (g = 0; g < cliqueSize - oldComponents[k].size(); g++) {
      DeltaStar[k][g + oldComponents[k].size()].resize(Q);
      int gg = g + oldComponents[k].size();
      for (q = 0; q < Q; q++) {
	int kqg = qg2index(q,g + first,Q,G);
	DeltaStar[k][gg][q] = Delta[kqg] / exp(0.5 * b[q] * log(sigma2[kqg]));
      }
    }

    vector<vector<double> > temp;
    quadratic2(DeltaStar[k],RInverse,temp);
    int g1,g2;
    for (g1 = 0; g1 < DNew[k].size(); g1++)
      for (g2 = 0; g2 < DNew[k][g1].size(); g2++)
	DNew[k][g1][g2] += temp[g1][g2];

    first += cliqueSize - oldComponents[k].size();
  }

  if (draw == 1)
    Omega = ran.HyperInverseWishart(dfNew,DNew,oldClique,oldComponents);

  double pot = ran.PotentialHyperInverseWishart(dfNew,DNew,oldClique,oldComponents,Omega);

  return pot;
}





double DeltaStarGibbs(const vector<int> &oldClique,
		      const vector<vector<int> > &oldComponents,
		      int Q, int G, const int *S,
		      double *Delta,const double *r,
		      const double *sigma2,const double *phi,
		      const double *tau2R,const double *b,
		      const double *nu,const int *delta,
		      const int *psi,const double *x,
		      const vector<vector<vector<double> > > &Omega,
		      Random &ran,int draw) {

  // compute prior presicion matrix

  //  cout << "start computing prior precision matrix" << endl;

  vector<vector<vector<double> > > OmegaInv(Omega);
  vector<double> OmegaDet(Omega.size(),0.0);
  int k;
  for (k = 0; k < OmegaInv.size(); k++)
    OmegaDet[k] = inverse(Omega[k],OmegaInv[k]);

  vector<vector<vector<double> > > OmegaSep;
  OmegaSep.resize(Omega.size());
  for (k = 0; k < OmegaSep.size(); k++) {
    OmegaSep[k].resize(oldComponents[k].size());
    int i;
    for (i = 0; i < oldComponents[k].size(); i++) {
      OmegaSep[k][i].resize(oldComponents[k].size());
      int j;
      for (j = 0; j < oldComponents[k].size(); j++) {
	OmegaSep[k][i][j] = Omega[oldClique[k]][oldComponents[k][i]][oldComponents[k][j]];
      }
    }
  }

  vector<vector<vector<double> > > OmegaSepInv(OmegaSep);
  vector<double> OmegaSepDet(Omega.size(),0.0);
  for (k = 0; k < OmegaSep.size(); k++) {
    if (OmegaSep[k].size() > 0)
      OmegaSepDet[k] = inverse(OmegaSep[k],OmegaSepInv[k]);
  }

  vector<map<int,double> > OmegaInvSparse;
  OmegaInvSparse.resize(G);
  int g;
  for (g = 0; g < G; g++)
    OmegaInvSparse[g].clear();

  vector<vector<int> > nr;
  nr.resize(Omega.size());
  g = 0;
  for (k = 0; k < Omega.size(); k++) {

    nr[k].resize(Omega[k].size());
    int gg;
    for (gg = 0; gg < oldComponents[k].size(); gg++)
      nr[k][gg] = nr[oldClique[k]][oldComponents[k][gg]];
    for (gg = oldComponents[k].size(); gg < Omega[k].size(); gg++) {
      nr[k][gg] = g;
      g++;
    }
  }

  // print to files
  /*
  for (k = 0; k < Omega.size(); k++) {
    char filename[120];
    sprintf(filename,"Omega-%d.txt",k);
    FILE *out = fopen(filename,"w");
    int i,j;
    for (i = 0; i < Omega[k].size(); i++) {
      for (j = 0; j < Omega[k][i].size(); j++)
	//fprintf(out,"%20.18e ",Omega[k][i][j]);
      //fprintf(out,"\n");
    }
    fclose(out);
  }

  for (k = 0; k < OmegaInv.size(); k++) {
    char filename[120];
    sprintf(filename,"OmegaInv-%d.txt",k);
    FILE *out = fopen(filename,"w");
    int i,j;
    for (i = 0; i < OmegaInv[k].size(); i++) {
      for (j = 0; j < OmegaInv[k][i].size(); j++)
	//fprintf(out,"%20.18e ",OmegaInv[k][i][j]);
      //fprintf(out,"\n");
    }
    fclose(out);
  }

  for (k = 1; k < OmegaSep.size(); k++) {
    char filename[120];
    sprintf(filename,"OmegaSep-%d.txt",k);
    FILE *out = fopen(filename,"w");
    int i,j;
    for (i = 0; i < OmegaSep[k].size(); i++) {
      for (j = 0; j < OmegaSep[k][i].size(); j++)
	//fprintf(out,"%20.18e ",OmegaSep[k][i][j]);
      //fprintf(out,"\n");
    }
    fclose(out);
  }

  for (k = 1; k < OmegaSepInv.size(); k++) {
    char filename[120];
    sprintf(filename,"OmegaSepInv-%d.txt",k);
    FILE *out = fopen(filename,"w");
    int i,j;
    for (i = 0; i < OmegaSepInv[k].size(); i++) {
      for (j = 0; j < OmegaSepInv[k][i].size(); j++)
	//fprintf(out,"%20.18e ",OmegaSepInv[k][i][j]);
      //fprintf(out,"\n");
    }
    fclose(out);
  }
  */
  // finished printing to files


  // establish and print large Omega
  /*
  vector<vector<double> > OmegaTotal;
  OmegaTotal.resize(G);
  for (g = 0; g < G; g++) {
    OmegaTotal[g].resize(G);
    int gg;
    for (gg = 0; gg < G; gg++)
      OmegaTotal[g][gg] = 0.0;
  }

  for (k = 0; k < OmegaInv.size(); k++) {
    int r;
    for (r = 0; r < nr[k].size(); r++) {
      int g = nr[k][r];
      int s;
      for (s = 0; s <= r; s++) {
	int gg = nr[k][s];
	double value = Omega[k][r][s];
	OmegaTotal[g][gg] = value;
	OmegaTotal[gg][g] = value;
      }
    }
  }

  {
    char filename[120];
    sprintf(filename,"OmegaTotal.txt");
    FILE *out = fopen(filename,"w");
    int i,j;
    for (i = 0; i < OmegaTotal.size(); i++) {
      for (j = 0; j < OmegaTotal[i].size(); j++)
	//fprintf(out,"%20.18e ",OmegaTotal[i][j]);
      //fprintf(out,"\n");
    }
    fclose(out);
  }
  */
  // finished printing files

  for (k = 0; k < OmegaInv.size(); k++) {
    int r;
    for (r = 0; r < nr[k].size(); r++) {
      int s;
      for (s = 0; s < nr[k].size(); s++) {
	int g = nr[k][r];
	int gg = nr[k][s];
	double value = OmegaInv[k][r][s];

	map<int,double>::iterator it;
	it = OmegaInvSparse[g].find(gg);
	if (it == OmegaInvSparse[g].end())
	  OmegaInvSparse[g].insert(pair<int,double>(gg,value));
	else {
	  OmegaInvSparse[g][gg] += value;
	}
      }
    }
  }

  for (k = 0; k < OmegaSepInv.size(); k++) {
    int r;
    for (r = 0; r < oldComponents[k].size(); r++) {
      int s;
      for (s = 0; s < oldComponents[k].size(); s++) {
	int g = nr[k][r];
	int gg = nr[k][s];
	double value = OmegaSepInv[k][r][s];

	OmegaInvSparse[g][gg] -= value;
      }
    }
  }

  // print OmegaInv to files
  /*
  {
    char filename[120];
    sprintf(filename,"OmegaTotalInv.txt");
    FILE *out = fopen(filename,"w");
    int i,j;
    for (i = 0; i < OmegaTotal.size(); i++) {
      for (j = 0; j < OmegaTotal[i].size(); j++) {
	double value = 0.0;
	if (OmegaInvSparse[i].find(j) != OmegaInvSparse[i].end())
	  value = OmegaInvSparse[i][j];
	//fprintf(out,"%20.18e ",value);
      }
      //fprintf(out,"\n");
    }
    fclose(out);
  }
  */
  // finished printing


  // establish covariance matrix R and its inverse

  vector<vector<double> > R;
  R.resize(Q);
  int p;
  for (p = 0; p < Q; p++) {
    R[p].resize(Q);
  }
  int q;
  for (p = 0; p < Q; p++) {
    R[p][p] = tau2R[p];
    for (q = p + 1; q < Q; q++) {
      R[p][q] = sqrt(tau2R[p] * tau2R[q]) * r[qq2index(p,q,Q)];
      R[q][p] = R[p][q];
    }
  }
  vector<vector<double> > RInverse;

  inverse(R,RInverse);

  // print RInverse to file
  /*
  {
    char filename[120];
    sprintf(filename,"RInverse.txt");
    FILE *out = fopen(filename,"w");
    int i,j;
    for (i = 0; i < RInverse.size(); i++) {
      for (j = 0; j < RInverse[i].size(); j++)
	//fprintf(out,"%20.18e ",RInverse[i][j]);
      //fprintf(out,"\n");
    }
    fclose(out);
  }
  */
  // finished printing RInverse


  // compute precision matrix of full conditional

  //  cout << "start computing precision matrix of full conditional" << endl;

  vector<map<int,double> > VinvSparse;
  VinvSparse.resize(G*Q);
  for (k = 0; k < G*Q; k++)
    VinvSparse[k].clear();

  for (g = 0; g < G; g++) {
    map<int,double>::iterator it;
    for (it = OmegaInvSparse[g].begin(); it != OmegaInvSparse[g].end(); it++) {
      int gg = it->first;
      double value = it->second;

      int qq;
      for (q = 0; q < Q; q++) {
	for (qq = 0; qq < Q; qq++) {
	  int index = g * Q + q;
      	  int indexindex = gg * Q + qq;
	  double cov = value * RInverse[q][qq];

	  VinvSparse[index].insert(pair<int,double>(indexindex,cov));
	}
      }
    }
  }

  // print precision matrix of full conditional
  /*
  {
    char filename[120];
    sprintf(filename,"VinvReduced.txt");
    FILE *out = fopen(filename,"w");
    int i,j;
    for (i = 0; i < VinvSparse.size(); i++) {
      for (j = 0; j < VinvSparse.size(); j++) {
	double value = 0.0;
	if (VinvSparse[i].find(j) != VinvSparse[i].end())
	  value = VinvSparse[i][j];
	//fprintf(out,"%20.18e ",value);
      }
      //fprintf(out,"\n");
    }
    fclose(out);
  }
  */
  // finished printing

  // add effect of data

  vector<double> l(VinvSparse.size(),0.0);
  vector<double> L(VinvSparse.size(),0.0);
  for (g = 0; g < G; g++)
    for (q = 0; q < Q; q++) {
      if (delta[qg2index(q,g,Q,G)] == 1) {
	int index = g * Q + q;

	double v0 = sigma2[qg2index(q,g,Q,G)] * phi[qg2index(q,g,Q,G)];
	double v1 = sigma2[qg2index(q,g,Q,G)] / phi[qg2index(q,g,Q,G)];
	double diag = 0.0;
	double ll = 0.0;
	int s;
	for (s = 0; s < S[q]; s++) {
	  double variance = psi[sq2index(s,q,S,Q)] == 0 ? v0 : v1;
	  diag += exp(b[q] * log(sigma2[qg2index(q,g,Q,G)])) / variance;

	  ll += (2.0 * psi[sq2index(s,q,S,Q)] - 1.0) * (x[sqg2index(s,q,g,S,Q,G)] - nu[qg2index(q,g,Q,G)]) / variance;
	}
	ll *= exp(0.5 * b[q] * log(sigma2[qg2index(q,g,Q,G)]));
	L[index] = diag;
	l[index] = ll;
      }
    }


  int index;
  for (index = 0; index < L.size(); index++)
    VinvSparse[index][index] += L[index];

  // print L and l
  /*
  {
    char filename[120];
    sprintf(filename,"l.txt");
    FILE *out = fopen(filename,"w");
    int i;
    for (i = 0; i < l.size(); i++) {
      //fprintf(out,"%20.18e\n",l[i]);
    }
    fclose(out);

    sprintf(filename,"L.txt");
    out = fopen(filename,"w");
    for (i = 0; i < L.size(); i++) {
      //fprintf(out,"%20.18e\n",L[i]);
    }
    fclose(out);

  }
  */
  // print precision matrix of full conditional
  /*
  {
    char filename[120];
    sprintf(filename,"Vinv.txt");
    FILE *out = fopen(filename,"w");
    int i,j;
    for (i = 0; i < VinvSparse.size(); i++) {
      for (j = 0; j < VinvSparse.size(); j++) {
	double value = 0.0;
	if (VinvSparse[i].find(j) != VinvSparse[i].end())
	  value = VinvSparse[i][j];
	//fprintf(out,"%20.18e ",value);
      }
      //fprintf(out,"\n");
    }
    fclose(out);
  }
  */
  // finished printing

  // establish a version of VinvSparse with the indices reversed

  vector<map<int,double> > VinvSparseReversed;
  VinvSparseReversed.resize(VinvSparse.size());
  for (k = 0; k < VinvSparseReversed.size(); k++)
    VinvSparseReversed[k].clear();
  for (k = 0; k < VinvSparseReversed.size(); k++) {
    map<int,double>::iterator it;
    for (it = VinvSparse[k].begin(); it != VinvSparse[k].end(); it++) {
      map<int,double>::iterator itextra = VinvSparse[k].end();

      int r = it->first;
      double value = it->second;

      int kk = VinvSparseReversed.size() - k - 1;
      int rr = VinvSparseReversed.size() - r - 1;
      VinvSparseReversed[rr].insert(pair<int,double>(kk,value));
    }
  }

  // print VinvSparseReversed
  /*
  {
    char filename[120];
    sprintf(filename,"VinvReversed.txt");
    FILE *out = fopen(filename,"w");
    int i,j;
    for (i = 0; i < VinvSparseReversed.size(); i++) {
      for (j = 0; j < VinvSparseReversed.size(); j++) {
	double value = 0.0;
	if (VinvSparseReversed[i].find(j) != VinvSparseReversed[i].end())
	  value = VinvSparseReversed[i][j];
	//fprintf(out,"%20.18e ",value);
      }
      //fprintf(out,"\n");
    }
    fclose(out);
  }
  */
  // finished printing

  // perform Cholesky decomposition for (the sparse matrix) VinvSparseReversed

  //  cout << "start Cholesky factorization" << endl;

  vector<map<int,double> > cholReversed;
  cholReversed.resize(Q * G);
  for (index = 0; index < cholReversed.size(); index++)
    cholReversed[index].clear();

  int N = cholReversed.size();
  int i;
  for (i = 0; i < N; i++) {
    map<int,double>::iterator it;
    for (it = VinvSparseReversed[i].find(i); it != VinvSparseReversed[i].end(); it++) {
      int j = it->first;
      double value = it->second;

      double sum = value;
      map<int,double>::iterator it2;
      for (it2 = cholReversed[i].begin(); it2 != cholReversed[i].end(); it2++) {
	if (it2->first < i) {
	  int k = it2->first;
	  double valuei = it2->second;
	  map<int,double>::iterator it3 = cholReversed[j].find(k);
	  if (it3 != cholReversed[j].end()) {
	    double valuej = it3->second;
	    sum -= valuei * valuej;
	  }
	}
      }

      //if (i == j && sum <= 0.0) {
	//fprintf(stderr,"DeltaStarGibbs: Matrix is not positive definite!\n");
      //exit(-1);
      //}

      if (i == j)
	cholReversed[j].insert(pair<int,double>(i,sqrt(sum)));
      else
	cholReversed[j].insert(pair<int,double>(i,sum / cholReversed[i][i]));
    }
  }
  /*
  // print (reversed) cholesky matrix

  {
    char filename[120];
    sprintf(filename,"cholReversed.txt");
    FILE *out = fopen(filename,"w");
    int i,j;
    for (i = 0; i < cholReversed.size(); i++) {
      for (j = 0; j < cholReversed.size(); j++) {
	double value = 0.0;
	if (cholReversed[i].find(j) != cholReversed[i].end())
	  value = cholReversed[i][j];
	//fprintf(out,"%20.18e ",value);
      }
      //fprintf(out,"\n");
    }
    fclose(out);
  }
  */
  // finished printing

  // establish er version of chol with the indices reversed back

  vector<map<int,double> > chol;
  chol.resize(cholReversed.size());
  for (k = 0; k < chol.size(); k++)
    chol[k].clear();
  for (k = 0; k < cholReversed.size(); k++) {
    map<int,double>::iterator it;
    for (it = cholReversed[k].begin(); it != cholReversed[k].end(); it++) {
      int r = it->first;
      double value = it->second;

      int kk = chol.size() - k - 1;
      int rr = chol.size() - r - 1;
      chol[rr].insert(pair<int,double>(kk,value));
    }
  }

  // print cholesky matrix
  /*
  {
    char filename[120];
    sprintf(filename,"chol.txt");
    FILE *out = fopen(filename,"w");
    int i,j;
    for (i = 0; i < chol.size(); i++) {
      for (j = 0; j < chol.size(); j++) {
	double value = 0.0;
	if (chol[i].find(j) != chol[i].end())
	  value = chol[i][j];
	//fprintf(out,"%20.18e ",value);
      }
      //fprintf(out,"\n");
    }
    fclose(out);
  }
  */
  // establish representation of the transpose of the cholesky matrix

  vector<map<int,double> > cholT;
  cholT.resize(chol.size());
  for (k = 0; k < cholT.size(); k++)
    cholT[k].clear();
  for (i = 0; i < chol.size(); i++) {
    map<int,double>::iterator it;
    for (it = chol[i].begin(); it != chol[i].end(); it++) {
      int j = it->first;
      double value = it->second;

      cholT[j].insert(pair<int,double>(i,value));
    }
  }

  // print cholT
  /*
  {
    char filename[120];
    sprintf(filename,"cholT.txt");
    FILE *out = fopen(filename,"w");
    int i,j;
    for (i = 0; i < cholT.size(); i++) {
      for (j = 0; j < cholT.size(); j++) {
	double value = 0.0;
	if (cholT[i].find(j) != cholT[i].end())
	  value = cholT[i][j];
	//fprintf(out,"%20.18e ",value);
      }
      //fprintf(out,"\n");
    }
    fclose(out);
  }
  */
  // finished printing

  //  cout << "start computing mean value" << endl;

  vector<double> u(l.size(),0.0);
  for (i = u.size() - 1; i >= 0; i--) {
    double diag = 0.0;
    double sum = 0.0;
    map<int,double>::iterator it;
    for (it = cholT[i].begin(); it != cholT[i].end(); it++) {
      int j= it->first;
      double value = it->second;

      if (i == j)
	diag = value;
      else
	sum += value * u[j];
    }
    u[i] = (l[i] - sum) / diag;
  }

  vector<double> mean(u.size(),0.0);
  for (i = 0; i < mean.size(); i++) {
    double diag = 0.0;
    double sum = 0.0;
    map<int,double>::iterator it;
    for (it = chol[i].begin(); it != chol[i].end(); it++) {
      int j = it->first;
      double value = it->second;

      if (j == i)
	diag = value;
      else
	sum += value * mean[j];
    }
    mean[i] = (u[i] - sum) / diag;
  }

  // print mean value
  /*
  {
    char filename[120];
    sprintf(filename,"mean.txt");
    FILE *out = fopen(filename,"w");
    int i;
    for (i = 0; i < mean.size(); i++) {
      //fprintf(out,"%20.18e\n",mean[i]);
    }
    fclose(out);
  }
  */
  // finished printing



  // generate a sample with zero mean, or compute the sample that should have been sampled

  //  cout << "start sampling" << endl;

  vector<double> sample(chol.size(),0.0);
  if (draw == 1) {
    vector<double> z(chol.size(),0.0);
    for (k = 0; k < z.size(); k++)
      z[k] = ran.Norm01();

    // print z
    /*
    {
      char filename[120];
      sprintf(filename,"z.txt");
      FILE *out = fopen(filename,"w");
      int i;
      for (i = 0; i < z.size(); i++) {
	//fprintf(out,"%20.18e\n",z[i]);
      }
      fclose(out);
    }
    */
    // finished printing



    for (i = 0; i < sample.size(); i++) {
      double diag = 0.0;
      double sum = 0.0;
      map<int,double>::iterator it;
      for (it = chol[i].begin(); it != chol[i].end(); it++) {
	int j = it->first;
	double value = it->second;

	if (j == i)
	  diag = value;
	else
	  sum += value * sample[j];
      }
      sample[i] = (z[i] - sum) / diag;
    }

    // print sample
    /*
  {
    char filename[120];
    sprintf(filename,"sample.txt");
    FILE *out = fopen(filename,"w");
    int i;
    for (i = 0; i < sample.size(); i++) {
      //fprintf(out,"%20.18e\n",sample[i]);
    }
    fclose(out);
  }
    */
  // finished printing

  }
  else {  // compute sample[i] necessary to generate the current values
    // compute DeltaStar

    vector<vector<double> > DeltaStar;
    DeltaStar.resize(G);
    for (g = 0; g < G; g++) {
      DeltaStar[g].resize(Q);
      for (q = 0; q < Q; q++) {
	DeltaStar[g][q] = Delta[qg2index(q,g,Q,G)] / exp(0.5 * b[q] * log(sigma2[qg2index(q,g,Q,G)]));
      }
    }

    // subtract mean value

    for (g = 0; g < G; g++)
      for (q = 0; q < Q; q++) {
	int index = g * Q + q;

	DeltaStar[g][q] -= mean[index];
      }

    // insert value in sample[]

    for (g = 0; g < G; g++)
      for (q = 0; q < Q; q++) {
	int index = g * Q + q;

	sample[index] = DeltaStar[g][q];
      }
  }

  // compute potential for sample

  double pot = 0.0;
  for (i = 0; i < sample.size(); i++) {
    double diag = 0.0;
    double sum = 0.0;
    map<int,double>::iterator it;
    for (it = chol[i].begin(); it != chol[i].end(); it++) {
      int j = it->first;
      double value = it->second;

      if (j == i)
	diag = value;
      else
	sum += value * sample[j];
    }
    double mean = - sum / diag;
    double variance = 1.0 / (diag * diag);
    pot += ran.PotentialGaussian(variance,mean,sample[i]);
  }

  if (draw == 1) { // add mean value and insert in data structure
    for (k = 0; k < Q * G; k++)
      sample[k] += mean[k];

    for (g = 0; g < G; g++)
      for (q = 0; q < Q; q++) {
	int index = g * Q + q;

	double DeltaStar = sample[index];

	Delta[qg2index(q,g,Q,G)] = DeltaStar * exp(0.5 * b[q] * log(sigma2[qg2index(q,g,Q,G)]));
      }
  }

  //  cout << "finished" << endl;

  return pot;
}




void transformGraph(const int *nClique,const int *oldClique,const int *nOldComponents,
		    const int *oldComponents,vector<int> &oldCliqueTransformed,
		    vector<vector<int> > &oldComponentsTransformed) {
  oldCliqueTransformed.resize(*nClique);
  oldComponentsTransformed.resize(*nClique);

  int nr = 0;
  int k;
  for (k = 0; k < *nClique; k++) {
    oldCliqueTransformed[k] = oldClique[k];

    oldComponentsTransformed[k].resize(nOldComponents[k]);
    int i;
    for (i = 0; i < nOldComponents[k]; i++) {
      oldComponentsTransformed[k][i] = oldComponents[nr];
      nr++;
    }
  }

  return;
}



void transformOmega(const int *nClique,const int *nOldComponents,
		    const int *nNewComponents,const double *Omega,
		    vector<vector<vector<double> > > &OmegaTransformed) {
  OmegaTransformed.resize(*nClique);

  int nr = 0;
  int k;
  for (k = 0; k < *nClique; k++) {
    int size = nOldComponents[k] + nNewComponents[k];
    OmegaTransformed[k].resize(size);
    int i,j;
    for (i = 0; i < size; i++)
      OmegaTransformed[k][i].resize(size);

    for (i = 0; i < size; i++)
      for (j = 0; j <= i; j++) {
	OmegaTransformed[k][i][j] = Omega[nr];
	OmegaTransformed[k][j][i] = Omega[nr];

	nr++;
      }
  }

  return;
}



void inverseTransformOmega(const vector<vector<vector<double> > > &OmegaTransformed,
			  double *Omega) {
  int nr = 0;
  int k;
  for (k = 0; k < OmegaTransformed.size(); k++) {
    int i,j;
    for (i = 0; i < OmegaTransformed[k].size(); i++)
      for (j = 0; j <= i; j++) {
	Omega[nr] = OmegaTransformed[k][i][j];

	nr++;
      }
  }

  return;
}
