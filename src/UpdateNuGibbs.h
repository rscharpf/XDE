#ifndef UPDATENUGIBBS_H
#define UPDATENUGIBBS_H

#include "Update.h"
#include "Matrix.h"


class UpdateNuGibbs : public Update
{
 public:

  UpdateNuGibbs(Structure *str,int check,const Potential *model);
  ~UpdateNuGibbs(void);

  int update(Random &ran);
  Update *copy(void) const;

 private:
  Structure *str;
  int check;
  Potential *model;
};



inline UpdateNuGibbs::UpdateNuGibbs(Structure *str,int check,const Potential *model) : Update(0.0)
{
  this->str = str;
  this->check = check;
  this->model = model->copy();

  return;
}



inline UpdateNuGibbs::~UpdateNuGibbs(void)
{
  delete model;

  return;
}


inline Update *UpdateNuGibbs::copy(void) const
{
  Update *u = new UpdateNuGibbs(str,check,model);

  return u;
}



inline int UpdateNuGibbs::update(Random &ran)
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
	  var[p][p] = str->gamma2 * str->tau2[p];
	  var[p][p] *= exp(str->a[p] * log(str->sigma2[p][g]));
	}

      for (p = 0; p < str->Q; p++)
	for (q = p + 1; q < str->Q; q++)
	  {
	    var[p][q] = str->gamma2;
	    var[p][q] *= str->rho[p][q];
	    var[p][q] *= sqrt(str->tau2[p] * str->tau2[q]);
	    var[p][q] *= exp(0.5 * (str->a[q] * log(str->sigma2[q][g]) + str->a[p] * log(str->sigma2[p][g])));

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

      vector<vector<double> > varInv;
      inverse(var,varInv);

      int s;
      for (q = 0; q < str->Q; q++)
	{
	  double var0 = str->sigma2[q][g] * str->phi[q][g];
	  double var1 = str->sigma2[q][g] / str->phi[q][g];
	  for (s = 0; s < str->S[q]; s++)
	    {
	      double variance = str->psi[q][s] == 0 ? var0 : var1;
	      varInv[q][q] += 1.0 / variance;
	      Mean[q] += (str->x[q][g][s] - str->delta[q][g] * (2.0 * str->psi[q][s] - 1.0) * 
			  str->Delta[q][g]) / variance;
	    }
	}
      
      inverse(varInv,var);
      vector<double> mean;
      matrixMult(var,Mean,mean);

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
	      oldValues[q] = str->nu[q][g];
	      str->nu[q][g] = newValues[q];
	    }

	  pot += model->potential(ran);
	  pot += ran.PotentialMultiGaussian(var,mean,oldValues);

	  for (q = 0; q < str->Q; q++)
	    str->nu[q][g] = oldValues[q];

	  if (pot >= 1.0e-6 || pot <= -1.0e-6)
	    cout << "WARNING: Possible implementation error in UpdateNuGibbs located. Check out!\n\n";
	}

      //
      // Set new values
      //

      for (q = 0; q < str->Q; q++)
	str->nu[q][g] = newValues[q];

      addAccept();
      nAccept++;
    }

  return nAccept;
}


#endif
