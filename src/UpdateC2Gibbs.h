#ifndef UPDATEC2GIBBS_H
#define UPDATEC2GIBBS_H

#include "Update.h"
#include "Matrix.h"


class UpdateC2Gibbs : public Update
{
 public:

  UpdateC2Gibbs(Structure *str,int check,const Potential *model);
  ~UpdateC2Gibbs(void);

  int update(Random &ran);
  Update *copy(void) const;

 private:
  Structure *str;
  int check;
  Potential *model;
};



inline UpdateC2Gibbs::UpdateC2Gibbs(Structure *str,int check,const Potential *model) : Update(0.0)
{
  this->str = str;
  this->check = check;
  this->model = model->copy();

  return;
}



inline UpdateC2Gibbs::~UpdateC2Gibbs(void)
{
  delete model;

  return;
}


inline Update *UpdateC2Gibbs::copy(void) const
{
  Update *u = new UpdateC2Gibbs(str,check,model);

  return u;
}



inline int UpdateC2Gibbs::update(Random &ran)
{
  int nAccept = 0;

  //
  // set prior parameters
  //

  double s = -1.0;
  double lambda = 0.0;

  //
  // update parameters based on available observations
  //

  int g;
  for (g = 0; g < str->G; g++)
    {
      if (str->delta[g] == 1)
	{
	  vector<vector<double> > varInv;
	  varInv.resize(str->Q);
	  int p,q;
	  for (p = 0; p < str->Q; p++) varInv[p].resize(str->Q);
	  for (p = 0; p < str->Q; p++)
	    for (q = p; q < str->Q; q++)
	      {
		varInv[p][q] = 1;
		if (p != q) varInv[p][q] *= str->r[p][q];
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
    }

  //
  // Draw new value
  //
  
  double newValue;
  if (s != -1.0 && s != 0.0) // there exist at least one delta == 1
    {
      int nTry = 0;
      do
	{
	  nTry++;
	  newValue = ran.InverseGamma(s,lambda);
	}
      while (newValue > str->c2Max && nTry < 100);
      if (nTry == 100) newValue = str->c2;   // propose an unchanged value!
    }
  else if (s == 0.0)
    {
      double fmax = exp(-1.0)/lambda;
      if (str->c2Max < lambda) fmax = exp(-lambda/str->c2Max)/str->c2Max;
      int nTry = 0;
      int accept = 0;
      do
	{
	  nTry++;
	  newValue = str->c2Max * ran.Unif01();
	  double alpha = (exp(-lambda/newValue)/newValue) / fmax;
	  accept = (ran.Unif01() <= alpha);
	  if (newValue > str->c2Max) accept = 0;
	}
      while (accept == 0 && nTry < 100);
      if (nTry == 100) newValue = str->c2;   // propose an unchanged value!
    }
  else
    newValue = str->c2Max * ran.Unif01();


  
  //
  // Check acceptance probability
  //
  
  if (check != 0)
    {
      double oldValue = str->c2;
      double pot = - model->potential(ran);
      pot -= ran.PotentialInverseGamma(s,lambda,newValue);
      
      str->c2 = newValue;

      pot += model->potential(ran);
      pot += ran.PotentialInverseGamma(s,lambda,oldValue);

      str->c2 = oldValue;
      
      if (pot >= 1.0e-6 || pot <= -1.0e-6)
	cout << "WARNING: Possible implementation error in UpdateC2Gibbs located. Check out!\n\n";
    }
  
  //
  // Set new value
  //
  
  str->c2 = newValue;

  //
  // draw new Delta values for genes with delta_g == 0
  //

  for (g = 0; g < str->G; g++)
    {
      if (str->delta[g] == 0)
	{
	  vector<vector<double> > var;
	  var.resize(str->Q);
	  int p,q;
	  for (p = 0; p < str->Q; p++) var[p].resize(str->Q);
	  for (p = 0; p < str->Q; p++)
	    {
	      var[p][p] = str->c2 * str->tau2[p];
	      var[p][p] *= exp(str->b[p] * log(str->sigma2[p][g]));
	    }
	  
	  for (p = 0; p < str->Q; p++)
	    for (q = p + 1; q < str->Q; q++)
	      {
		var[p][q] = str->c2;
		var[p][q] *= str->r[p][q];
		var[p][q] *= sqrt(str->tau2[p] * str->tau2[q]);
		var[p][q] *= exp(0.5 * (str->b[q] * log(str->sigma2[q][g]) + str->b[p] * log(str->sigma2[p][g])));
		
		var[q][p] = var[p][q];
	      }
	  
	  vector<double> mean(str->Q);
	  for (p = 0; p < str->Q; p++) mean[p] = 0.0;
	  
	  vector<double> newDDelta = ran.MultiGaussian(var,mean);
	  for (p = 0; p < str->Q; p++)
	    str->Delta[p][g] = newDDelta[p];
	}
    }
  
  
  addTry();
  addAccept();
  nAccept++;
  
  return nAccept;
}


#endif
