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
	oldRho[p][q] = str->rho[p][q];
    }

  //
  // propose new values for rho
  //

  vector<vector<double> > T(str->Q);
  for (p = 0; p < str->Q; p++)
    T[p].resize(str->Q);
  if (ran.Unif01() <= 0.5)
    T = ran.CorrelationStandardInverseWishart(str->Q,str->nuR);
  else
    {
      int m = str->Q;
      double deltaRho = -1.0 / (m - 1.0) + (1.0 + 1.0 / (m - 1.0)) * ran.Unif01();
      for (p = 0; p < str->Q; p++)
	T[p][p] = 1.0;
      for (p = 0; p < str->Q; p++)
	for (q = p+1; q < str->Q; q++)
	  {
	    T[p][q] = deltaRho;
	    T[q][p] = deltaRho;
	  }
    }

  double pot = 0.0;
  if (ran.Unif01() <= 0.5)
    {
      for (p = 0; p < str->Q; p++)
	for (q = 0; q < str->Q; q++)
	  {
	    if (p != q) 
	      newRho[p][q] = (1.0 - epsilon) * oldRho[p][q] + epsilon * T[p][q];
	    else
	      newRho[p][q] = 1.0;
	  }
      pot += - str->Q * (str->Q - 1.0) * log(1.0 - epsilon) / 2.0;
    }
  else
    {
      for (p = 0; p < str->Q; p++)
	for (q = 0; q < str->Q; q++)
	  {
	    if (p != q) 
	      newRho[p][q] = oldRho[p][q] / (1.0 - epsilon) - epsilon * T[p][q] / (1.0 - epsilon);
	    else
	      newRho[p][q] = 1.0;
	  }
      pot += str->Q * (str->Q - 1) * log(1.0 - epsilon) / 2.0;
    }


  //
  // if any of the proposed new values are negative, reject the proposal
  //

  int isNeg = 0;
  for (p = 0; p < str->Q; p++)
    for (q = 0; q < str->Q; q++)
      isNeg += (newRho[p][q] < 0.0);

  if (isNeg > 0)
      return nAccept;


  int err = 0;
  Cholesky chol(newRho,err);

  if (err == 0)  // if err == 1: proposal is not positive definite
    {
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
		varInv[p][q] *= sqrt(str->tau2[p] * str->tau2[q]);
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
		varInv[p][q] *= sqrt(str->tau2[p] * str->tau2[q]);
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
    }


  return nAccept;
}




#endif
