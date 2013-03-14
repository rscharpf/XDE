#ifndef UPDATEGAMMA2GIBBS_H
#define UPDATEGAMMA2GIBBS_H

#include "Update.h"
#include "Matrix.h"


class UpdateGamma2Gibbs : public Update
{
 public:

  UpdateGamma2Gibbs(Structure *str,int check,const Potential *model);
  ~UpdateGamma2Gibbs(void);

  int update(Random &ran);
  Update *copy(void) const;

 private:
  Structure *str;
  int check;
  Potential *model;
};



inline UpdateGamma2Gibbs::UpdateGamma2Gibbs(Structure *str,int check,const Potential *model) : Update(0.0)
{
  this->str = str;
  this->check = check;
  this->model = model->copy();

  return;
}



inline UpdateGamma2Gibbs::~UpdateGamma2Gibbs(void)
{
  delete model;

  return;
}


inline Update *UpdateGamma2Gibbs::copy(void) const
{
  Update *u = new UpdateGamma2Gibbs(str,check,model);

  return u;
}



inline int UpdateGamma2Gibbs::update(Random &ran)
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
      vector<vector<double> > varInv;
      varInv.resize(str->Q);
      int p,q;
      for (p = 0; p < str->Q; p++) varInv[p].resize(str->Q);
      for (p = 0; p < str->Q; p++)
	for (q = p; q < str->Q; q++)
	  {
	    varInv[p][q] = 1;
	    if (p != q) varInv[p][q] *= str->rho[p][q];
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

  //
  // Draw new value
  //

  double newValue = ran.InverseGamma(s,lambda);

  //
  // Check acceptance probability
  //

  if (check != 0)
    {
      double oldValue = str->gamma2;
      double pot = - model->potential(ran);
      pot -= ran.PotentialInverseGamma(s,lambda,newValue);

      str->gamma2 = newValue;

      pot += model->potential(ran);
      pot += ran.PotentialInverseGamma(s,lambda,oldValue);

      str->gamma2 = oldValue;

//      if (pot >= 1.0e-6 || pot <= -1.0e-6)
//	cout << "WARNING: Possible implementation error in UpdateGamma2Gibbs located. Check out!\n\n";
    }

  //
  // Set new value
  //

  str->gamma2 = newValue;

  addTry();
  addAccept();
  nAccept++;

  return nAccept;
}


#endif
