#ifndef UPDATETAU2RDDELTAMH_H
#define UPDATETAU2RDDELTAMH_H

#include "Update.h"


class UpdateTau2RDDeltaMH : public Update
{
 public:

  UpdateTau2RDDeltaMH(Structure *str,const Potential *model,double epsilon);
  ~UpdateTau2RDDeltaMH(void);

  int update(Random &ran);
  Update *copy(void) const;

 private:
  Structure *str;
  Potential *model;
};



inline UpdateTau2RDDeltaMH::UpdateTau2RDDeltaMH(Structure *str,const Potential *model,double epsilon) : Update(epsilon)
{
  this->str = str;
  this->model = model->copy();
  
  return;
}



inline UpdateTau2RDDeltaMH::~UpdateTau2RDDeltaMH(void)
{
  delete model;

  return;
}


inline Update *UpdateTau2RDDeltaMH::copy(void) const
{
  Update *u = new UpdateTau2RDDeltaMH(str,model,epsilon);

  return u;
}




inline int UpdateTau2RDDeltaMH::update(Random &ran)
{
  int nAccept = 0;

  if (str->Q > 1)
    {
      int q = (int) (str->Q * ran.Unif01());
      int p = (int) ((str->Q - 1) * ran.Unif01());
      if (p >= q) p++;
      
      double upper = 1.0 + epsilon;
      double lower = 1.0 / upper;
      
      double u = lower + (upper - lower) * ran.Unif01();
      vector<double> oldValues(str->Q);
      vector<double> newValues(str->Q);
      
      int i;
      for (i = 0; i < str->Q; i++)
	{
	  oldValues[i] = str->tau2R[i];
	  newValues[i] = str->tau2R[i];
	}
      
      newValues[q] *= u;
      newValues[p] /= u;
      
      double prod = 1.0;
      for (i = 0; i < str->Q; i++)
	prod *= newValues[i];
      
      prod = exp(log(prod) / str->Q);
      for (i = 0; i < str->Q; i++)
	newValues[i] /= prod;

      double pot = - log(1.0 / (u * u));

      //
      // propose new values for Delta from full conditionals
      //

      vector<vector<double> > newDDelta;
      vector<vector<double> > newDDeltaOldTau2;
      newDDelta.resize(str->Q);
      newDDeltaOldTau2.resize(str->Q);
      for (q = 0; q < str->Q; q++)
	{
	  newDDelta[q].resize(str->G);
	  newDDeltaOldTau2[q].resize(str->G);
	}
      
      int g;
      for (g = 0; g < str->G; g++)
	{
	  //
	  // compute prior covariance matrix
	  //

	  vector<vector<double> > var;
	  var.resize(str->Q);
	  int p,q;
	  for (p = 0; p < str->Q; p++)
	    var[p].resize(str->Q);
	  for (p = 0; p < str->Q; p++)
	    {
	      var[p][p] = str->c2 * newValues[p];
	      var[p][p] *= exp(str->b[p] * log(str->sigma2[p][g]));
	    }
	  
	  for (p = 0; p < str->Q; p++)
	    for (q = p + 1; q < str->Q; q++)
	      {
		var[p][q] = str->c2;
		var[p][q] *= str->r[p][q];
		var[p][q] *= sqrt(newValues[p] * newValues[q]);
		var[p][q] *= exp(0.5 * (str->b[q] * log(str->sigma2[q][g]) + str->b[p] * log(str->sigma2[p][g])));
		
		var[q][p] = var[p][q];
	      }
	  
	  //
	  // define prior mean
	  //
	  
	  vector<double> Mean(str->Q);
	  for (p = 0; p < str->Q; p++) Mean[p] = 0.0;
	  
	  //
	  // Update parameters based on available observations 
	  //
	  
	  vector<double> mean(str->Q);
	  for (p = 0; p < str->Q; p++) mean[p] = 0.0;
	  
	  int modify = 0;
	  for (q = 0; q < str->Q; q++)
	    {
	      if (str->delta[q][g] == 1)
		modify = 1;
	    }


	  //
	  // Compute extra linear and quadratic terms
	  //

	  vector<double> lin;
	  vector<double> quad;
	  lin.resize(str->Q);
	  quad.resize(str->Q);
	  for (q = 0; q < str->Q; q++)
	    {
	      lin[q] = 0.0;
	      quad[q] = 0.0;
	    }
	  if (modify != 0)
	    {
	      for (q = 0; q < str->Q; q++)
		{
		  if (str->delta[q][g] == 1)
		    {
		      double var0 = str->sigma2[q][g] * str->phi[q][g];
		      double var1 = str->sigma2[q][g] / str->phi[q][g];
		      int s;
		      for (s = 0; s < str->S[q]; s++)
			{
			  double variance = str->psi[q][s] == 0 ? var0 : var1;
			  quad[q] += 1.0 / variance;
			  lin[q] += (2.0 * str->psi[q][s] - 1.0) * 
			    (str->x[q][g][s] - str->nu[q][g]) / variance;
			}
		    }
		}
	    }

	  vector<vector<double> > varInv;
	  double detPrior = inverse(var,varInv);
	  double detPosterior = detPrior;
	  if (modify != 0)
	    {
	      detPrior = inverse(var,varInv);
	      
	      for (q = 0; q < str->Q; q++)
		{
		  Mean[q] += lin[q];
		  varInv[q][q] += quad[q];
		}

	      detPosterior = 1.0 / inverse(varInv,var);
	      matrixMult(var,Mean,mean);
	    }
	  
	  //
	  // Draw new values
	  //
	  
	  vector<double> vv(str->Q);
	  vv = ran.MultiGaussian(var,mean);
	  for (q = 0; q < str->Q; q++)
	    newDDelta[q][g] = vv[q];

	  pot += 0.5 * log(detPrior) - 0.5 * log(detPosterior);
	  pot += - 0.5 * quadratic(varInv,mean);

	  //
	  // add potential for reverse proposal
	  //

	  //
	  // compute prior covariance matrix
	  //

	  for (p = 0; p < str->Q; p++)
	    {
	      var[p][p] = str->c2 * oldValues[p];
	      var[p][p] *= exp(str->b[p] * log(str->sigma2[p][g]));
	    }
	  
	  for (p = 0; p < str->Q; p++)
	    for (q = p + 1; q < str->Q; q++)
	      {
		var[p][q] = str->c2;
		var[p][q] *= str->r[p][q];
		var[p][q] *= sqrt(oldValues[p] * oldValues[q]);
		var[p][q] *= exp(0.5 * (str->b[q] * log(str->sigma2[q][g]) + str->b[p] * log(str->sigma2[p][g])));
		
		var[q][p] = var[p][q];
	      }
	  
	  //
	  // define prior mean
	  //
	  
	  for (p = 0; p < str->Q; p++) Mean[p] = 0.0;
	  
	  //
	  // Update parameters based on available observations 
	  //
	  
	  for (p = 0; p < str->Q; p++) mean[p] = 0.0;
	  
	  detPrior = inverse(var,varInv);
	  detPosterior = detPrior;
	  if (modify != 0)
	    {
	      for (q = 0; q < str->Q; q++)
		{
		  Mean[q] += lin[q];
		  varInv[q][q] += quad[q];
		}
	      
	      detPosterior = 1.0 / inverse(varInv,var);
	      matrixMult(var,Mean,mean);
	    }
	  
	  vv = ran.MultiGaussian(var,mean);
	  for (q = 0; q < str->Q; q++)
	    newDDeltaOldTau2[q][g] = vv[q];
	  
	  pot -= 0.5 * log(detPrior) - 0.5 * log(detPosterior);
	  pot -= - 0.5 * quadratic(varInv,mean);
	}

      pot -= model->potential(ran);
      
      for (i = 0; i < str->Q; i++)
	str->tau2R[i] = newValues[i];
      
      pot += model->potential(ran);
      
      for (i = 0; i < str->Q; i++)
	str->tau2R[i] = oldValues[i];


      //      cout << "UpdateTau2RDDeltaMH: " << pot << "\n";

      addTry();
      if (ran.Unif01() <= exp(- pot))
	{
	  for (i = 0; i < str->Q; i++)
	    str->tau2R[i] = newValues[i];
	  for (q = 0; q < str->Q; q++)
	    for (g = 0; g < str->G; g++)
	      str->Delta[q][g] = newDDelta[q][g];
	  addAccept();
	  nAccept++;
	}
      else
	{
	  for (q = 0; q < str->Q; q++)
	    for (g = 0; g < str->G; g++)
	      str->Delta[q][g] = newDDeltaOldTau2[q][g];
	}
    }
  
  return nAccept;
}


#endif
