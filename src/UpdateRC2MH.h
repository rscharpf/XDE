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
	oldR[p][q] = str->r[p][q];
    }

  //
  // propose new values for r
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
	      newR[p][q] = (1.0 - epsilon) * oldR[p][q] + epsilon * T[p][q];
	    else
	      newR[p][q] = 1.0;
	  }
      pot += - str->Q * (str->Q - 1.0) * log(1.0 - epsilon) / 2.0;
    }
  else
    {
      for (p = 0; p < str->Q; p++)
	for (q = 0; q < str->Q; q++)
	  {
	    if (p != q) 
	      newR[p][q] = oldR[p][q] / (1.0 - epsilon) - epsilon * T[p][q] / (1.0 - epsilon);
	    else
	      newR[p][q] = 1.0;
	  }
      pot += str->Q * (str->Q - 1) * log(1.0 - epsilon) / 2.0;
    }


  //
  // if any of the proposed new values are negative, reject the proposal
  //

  int isNeg = 0;
  for (p = 0; p < str->Q; p++)
    for (q = 0; q < str->Q; q++)
      isNeg += (newR[p][q] < 0.0);

  if (isNeg > 0)
      return nAccept;


  int err = 0;
  Cholesky chol(newR,err);
  if (err == 0)  // if err == 1: proposal is not positive definite
    {
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
		varInv[p][q] *= sqrt(str->tau2[p] * str->tau2[q]);
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
		varInv[p][q] *= sqrt(str->tau2[p] * str->tau2[q]);
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
    }
  
  return nAccept;
}




#endif
