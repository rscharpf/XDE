#ifndef UPDATEDDELTAGIBBS_H
#define UPDATEDDELTAGIBBS_H

#include "Update.h"


class UpdateDDeltaGibbs : public Update
{
 public:

  UpdateDDeltaGibbs(Structure *str,int check,const Potential *model);
  ~UpdateDDeltaGibbs(void);

  int update(Random &ran);
  Update *copy(void) const;

 private:
  Structure *str;
  int check;
  Potential *model;
};



inline UpdateDDeltaGibbs::UpdateDDeltaGibbs(Structure *str,int check,const Potential *model) : Update(0.0)
{
  this->str = str;
  this->check = check;
  this->model = model->copy();

  return;
}



inline UpdateDDeltaGibbs::~UpdateDDeltaGibbs(void)
{
  delete model;

  return;
}


inline Update *UpdateDDeltaGibbs::copy(void) const
{
  Update *u = new UpdateDDeltaGibbs(str,check,model);

  return u;
}



inline int UpdateDDeltaGibbs::update(Random &ran)
{
  int nAccept = 0;

  int g;
  for (g = 0; g < str->G; g++)
    {
      addTry();

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

      if (modify != 0)
	{
	  vector<vector<double> > varInv;
	  inverse(var,varInv);
	  
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
		      varInv[q][q] += 1.0 / variance;
		      Mean[q] += (2.0 * str->psi[q][s] - 1.0) * (str->x[q][g][s] - str->nu[q][g]) / variance;
		    }
		}
	    }
	  
	  inverse(varInv,var);
	  matrixMult(var,Mean,mean);
	}

      //
      // Draw new values
      //

      vector<double> newValues(str->Q);
      newValues = ran.MultiGaussian(var,mean);

      //
      // Check acceptance probability
      //

      if (check != 0)
	{
	  vector<double> oldValues(str->Q);
	  double pot = - model->potential(ran);
	  pot -= ran.PotentialMultiGaussian(var,mean,newValues);

	  for (q = 0; q < str->Q; q++)
	    {
	      oldValues[q] = str->Delta[q][g];
	      str->Delta[q][g] = newValues[q];
	    }
	  
	  pot += model->potential(ran);
	  pot += ran.PotentialMultiGaussian(var,mean,oldValues);

	  for (q = 0; q < str->Q; q++)
	    str->Delta[q][g] = oldValues[q];

	  if (pot >= 1.0e-6 || pot <= -1.0e-6)
	    cout << "WARNING: Possible implementation error in UpdateDDeltaGibbs located. Check out!\n\n";
	}

      //
      // Set new values
      //

      for (q = 0; q < str->Q; q++)
	str->Delta[q][g] = newValues[q];

      addAccept();
      nAccept++;
    }

  return nAccept;
}


#endif
