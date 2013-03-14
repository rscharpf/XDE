#ifndef UPDATERC2MH_H
#define UPDATERC2MH_H

#include "Update.h"
#include "Cholesky.h"

class UpdateRC2MH : public Update
{
 public:

  UpdateRC2MH(Structure *str,const Potential *model,double epsilon);
  ~UpdateRC2MH(void);

  int update(Random &ran);
  Update *copy(void) const;

 private:
  Structure *str;
  Potential *model;
};



inline UpdateRC2MH::UpdateRC2MH(Structure *str,const Potential *model,double epsilon) : Update(epsilon)
{
  this->str = str;
  this->model = model->copy();

  return;
}



inline UpdateRC2MH::~UpdateRC2MH(void)
{
  delete model;

  return;
}


inline Update *UpdateRC2MH::copy(void) const
{
  Update *u = new UpdateRC2MH(str,model,epsilon);

  return u;
} 



inline int UpdateRC2MH::update(Random &ran)
{
  int nAccept = 0;

  addTry();

  vector<vector<double> > oldR;
  vector<vector<double> > newR;
  oldR.resize(str->Q);
  newR.resize(str->Q);
  int p,q;
  for (p = 0; p < str->Q; p++)
    {
      oldR[p].resize(str->Q);
      newR[p].resize(str->Q);
      for (q = 0; q < str->Q; q++)
	{
	  oldR[p][q] = str->r[p][q];
	  newR[p][q] = str->r[p][q];
	}
    }

  //
  // propose new values for r
  //  

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

  newR[pp][qq] = oldR[pp][qq] * exp(u) / 
    (1 - oldR[pp][qq] + oldR[pp][qq] * exp(u));
  newR[qq][pp] = newR[pp][qq];
  
  //
  // check that potential new correlation matrix is positive definite
  //

  int err = 0;
  Cholesky chol(newR,err);
  if (err != 0)
    return nAccept;

  //
  // compute Jacobian determinant
  //

  double y = oldR[pp][qq];
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
      isNeg += (newR[p][q] < 0.0);

  if (isNeg > 0)
  return nAccept; 


  //
  // propose new value for c2
  //
  
  double oldC2 = str->c2;
  double newC2;
  
  double s = -1.0;
  double lambda = 0.0;
  
  int g;
  for (g = 0; g < str->G; g++)
    {
      vector<vector<double> > varInv;
      varInv.resize(str->Q);
      int p,q;
      for (p = 0; p < str->Q; p++) varInv[p].resize(str->Q);
      for (p = 0; p < str->Q; p++)
	for (q = p; q < str->Q; q++)
	  {
	    varInv[p][q] = 1;
	    if (p != q) varInv[p][q] *= newR[p][q];
	    varInv[p][q] *= sqrt(str->tau2R[p] * str->tau2R[q]);
	    varInv[p][q] *= exp(0.5 * (str->b[q] * log(str->sigma2[q][g]) + str->b[p] * log(str->sigma2[p][g])));
	    
	    varInv[q][p] = varInv[p][q];
	  }
      vector<vector<double> > var;
      inverse(varInv,var);
      
      vector<double> Delta(str->Q);
      for (q = 0; q < str->Q; q++)
	Delta[q] = str->Delta[q][g];
      
      double value = quadratic(var,Delta);
      lambda += 0.5 * value;
      s += 0.5 * str->Q;
    }
  
  if (s > 0.0)
    {
      newC2 = ran.InverseGamma(s,lambda);
      pot -= ran.PotentialInverseGamma(s,lambda,newC2);
    }
  else
    {
      newC2 = str->c2Max * ran.Unif01();
      pot -= - log(1.0 / str->c2Max);
    }
  
  //
  // compute potential for inverse proposal for c2
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
	    if (p != q) varInv[p][q] *= oldR[p][q];
	    varInv[p][q] *= sqrt(str->tau2R[p] * str->tau2R[q]);
	    varInv[p][q] *= exp(0.5 * (str->b[q] * log(str->sigma2[q][g]) + str->b[p] * log(str->sigma2[p][g])));
	    
	    varInv[q][p] = varInv[p][q];
	  }
      vector<vector<double> > var;
      inverse(varInv,var);
      
      vector<double> Delta(str->Q);
      for (q = 0; q < str->Q; q++)
	Delta[q] = str->Delta[q][g];
      
      double value = quadratic(var,Delta);
      lambda += 0.5 * value;
      s += 0.5 * str->Q;
    }
  
  if (s > 0.0)
    pot += ran.PotentialInverseGamma(s,lambda,oldC2);
  else
    pot += - log(1.0 / str->c2Max);
  
  //
  // compute potentials for new and old states
  //
  
  if (newC2 <= str->c2Max) // otherwise do not accept
    {
      pot -= model->potential(ran);
      
      for (p = 0; p < str->Q; p++)
	for (q = 0; q < str->Q; q++)
	  str->r[p][q] = newR[p][q];
      str->c2 = newC2;
      
      pot += model->potential(ran);
      
      for (p = 0; p < str->Q; p++)
	for (q = 0; q < str->Q; q++)
	  str->r[p][q] = oldR[p][q];
      str->c2 = oldC2;
      
      
      if (ran.Unif01() <= exp(- pot))
	{
	  nAccept++;
	  addAccept();
	  
	  for (p = 0; p < str->Q; p++)
	    for (q = 0; q < str->Q; q++)
	      str->r[p][q] = newR[p][q];
	  str->c2 = newC2;
	}
    }

  
  return nAccept;
}




#endif
