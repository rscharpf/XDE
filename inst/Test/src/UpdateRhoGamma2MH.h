#ifndef UPDATERHOGAMMA2MH_H
#define UPDATERHOGAMMA2MH_H

#include "Update.h"
#include "Cholesky.h"


class UpdateRhoGamma2MH : public Update
{
 public:

  UpdateRhoGamma2MH(Structure *str,const Potential *model,double epsilon);
  ~UpdateRhoGamma2MH(void);

  int update(Random &ran);
  Update *copy(void) const;

 private:
  Structure *str;
  Potential *model;
};



inline UpdateRhoGamma2MH::UpdateRhoGamma2MH(Structure *str,const Potential *model,double epsilon) : Update(epsilon)
{
  this->str = str;
  this->model = model->copy();

  return;
}



inline UpdateRhoGamma2MH::~UpdateRhoGamma2MH(void)
{
  delete model;

  return;
}


inline Update *UpdateRhoGamma2MH::copy(void) const
{
  Update *u = new UpdateRhoGamma2MH(str,model,epsilon);

  return u;
} 



inline int UpdateRhoGamma2MH::update(Random &ran)
{
  int nAccept = 0;

  addTry();

  vector<vector<double> > oldRho;
  vector<vector<double> > newRho;
  oldRho.resize(str->Q);
  newRho.resize(str->Q);
  int p,q;
  for (p = 0; p < str->Q; p++)
    {
      oldRho[p].resize(str->Q);
      newRho[p].resize(str->Q);
      for (q = 0; q < str->Q; q++)
	{
	  oldRho[p][q] = str->rho[p][q];
	  newRho[p][q] = str->rho[p][q];
	}
    }

  double pot = 0.0;
  double u = epsilon * ran.Norm01();

  //
  // draw what element to change
  //

  vector<double> prob(str->Q);
  for (p = 0; p < str->Q; p++)
    prob[p] = 1.0 / ((double) str->Q);
  int pp = ran.Discrete(prob);
  prob.resize(str->Q - 1);
  for (p = 0; p < str->Q - 1; p++)
    prob[p] = 1.0 / ((double) (str->Q - 1));
  
  int qq = ran.Discrete(prob);
  qq += (qq >= pp);
  
  //
  // compute potential new correlation value
  //

  newRho[pp][qq] = oldRho[pp][qq] * exp(u) / 
    (1 - oldRho[pp][qq] + oldRho[pp][qq] * exp(u));
  newRho[qq][pp] = newRho[pp][qq];
  
  // check that potential new correlation matrix is positive definite
  //

  int err = 0;
  Cholesky chol(newRho,err);
  if (err != 0)
  return nAccept;

  //
  // compute Jacobian determinant
  //

  double y = oldRho[pp][qq];
  double xtilde = log(y) - log(1.0 - y) + u;
  double pot1;
  if (xtilde <= 0.0)
    pot1 = - log(1.0 + exp(xtilde));
  else
    pot1 = - xtilde - log(1.0 + exp(- xtilde));
  double potdytildedxtilde = - xtilde - 2.0 * pot1;
  
  double potdxtildedy = log(1.0 - y) + log(y);

  pot += potdytildedxtilde + potdxtildedy;
		 
  //
  // if any of the proposed new values are negative, reject the proposal
  //

  int isNeg = 0;
  for (p = 0; p < str->Q; p++)
    for (q = 0; q < str->Q; q++)
      isNeg += (newRho[p][q] < 0.0);

  if (isNeg > 0)
    return nAccept;


  //
  // propose new values for gamma2
  //
  
  double oldGamma2 = str->gamma2;
  double newGamma2;
  
  double s = -1.0;
  double lambda = 0.0;
  int g;
  for (g = 0; g < str->G; g++)
    {
      vector<vector<double> > varInv;
      varInv.resize(str->Q);
      int p,q;
      for (p = 0; p < str->Q; p++)
	varInv[p].resize(str->Q);
      for (p = 0; p < str->Q; p++)
	for (q = p; q < str->Q; q++)
	  {
	    varInv[p][q] = 1;
	    if (p != q) varInv[p][q] *= newRho[p][q];
	    varInv[p][q] *= sqrt(str->tau2Rho[p] * str->tau2Rho[q]);
	    varInv[p][q] *= exp(0.5 * (str->a[q] * log(str->sigma2[q][g]) + str->a[p] * log(str->sigma2[p][g])));
	    
	    varInv[q][p] = varInv[p][q];
	  }
      vector<vector<double> > var;
      inverse(varInv,var);
      
      vector<double> nu(str->Q);
      for (q = 0; q < str->Q; q++)
	nu[q] = str->nu[q][g];
      
      double value = quadratic(var,nu);
      lambda += 0.5 * value;
      s += 0.5 * str->Q;
    }
  
  newGamma2 = ran.InverseGamma(s,lambda);
  pot -= ran.PotentialInverseGamma(s,lambda,newGamma2);
  
  //
  // compute potential for inverse proposal for gamma2
  //
  
  s = -1.0;
  lambda = 0.0;
  for (g = 0; g < str->G; g++)
    {
      vector<vector<double> > varInv;
      varInv.resize(str->Q);
      int p,q;
      for (p = 0; p < str->Q; p++)
	varInv[p].resize(str->Q);
      for (p = 0; p < str->Q; p++)
	for (q = p; q < str->Q; q++)
	  {
	    varInv[p][q] = 1;
	    if (p != q) varInv[p][q] *= oldRho[p][q];
	    varInv[p][q] *= sqrt(str->tau2Rho[p] * str->tau2Rho[q]);
	    varInv[p][q] *= exp(0.5 * (str->a[q] * log(str->sigma2[q][g]) + str->a[p] * log(str->sigma2[p][g])));
	    
	    varInv[q][p] = varInv[p][q];
	  }
      vector<vector<double> > var;
      inverse(varInv,var);
      
      vector<double> nu(str->Q);
      for (q = 0; q < str->Q; q++)
	nu[q] = str->nu[q][g];
      
      double value = quadratic(var,nu);
      lambda += 0.5 * value;
      s += 0.5 * str->Q;
    }
  
  pot += ran.PotentialInverseGamma(s,lambda,oldGamma2);
  
  //
  // compute potential for new and old states
  //
  
  pot -= model->potential(ran);
  
  for (p = 0; p < str->Q; p++)
    for (q = 0; q < str->Q; q++)
      str->rho[p][q] = newRho[p][q];
  str->gamma2 = newGamma2;
  
  pot += model->potential(ran);
  
  for (p = 0; p < str->Q; p++)
    for (q = 0; q < str->Q; q++)
      str->rho[p][q] = oldRho[p][q];
  
  str->gamma2 = oldGamma2;
  
  if (ran.Unif01() <= exp(- pot))
    {
      nAccept++;
      addAccept();
      
      for (p = 0; p < str->Q; p++)
	for (q = 0; q < str->Q; q++)
	  str->rho[p][q] = newRho[p][q];
      str->gamma2 = newGamma2;
    }

  return nAccept;
}




#endif
