#ifndef UPDATEDELTAGIBBS_H
#define UPDATEDELTAGIBBS_H

#include "Update.h"
#include "PotentialDeltag.h"
#include "PotentialDDeltag.h"
#include "PotentialXqg.h"


class UpdateDeltaMH : public Update
{
 public:

  UpdateDeltaMH(Structure *str,int oneDelta);
  ~UpdateDeltaMH(void);

  int update(Random &ran);
  Update *copy(void) const;

 private:
  Structure *str;
  int oneDelta;

  vector<Potential *> model;
};



inline UpdateDeltaMH::UpdateDeltaMH(Structure *str, int oneDelta = 0) : Update(0.0)
{
  this->str = str;
  this->oneDelta = oneDelta;

  int g;
  for (g = 0; g < str->G; g++)
    {
      vector<Potential *> term;
      term.push_back(new PotentialDeltag(g,str,oneDelta));
      term.push_back(new PotentialDDeltag(g,str));
      int q;
      for (q = 0; q < str->Q; q++)
	term.push_back(new PotentialXqg(q,g,str));

      model.push_back(new PotentialSum(term));

      int i;
      for (i = 0; i < term.size(); i++)
	delete term[i];
    }

  return;
}



inline UpdateDeltaMH::~UpdateDeltaMH(void)
{
  int i;
  for (i = 0; i < model.size(); i++)
    delete model[i];

  return;
}


inline Update *UpdateDeltaMH::copy(void) const
{
  Update *u = new UpdateDeltaMH(str);

  return u;
}



inline int UpdateDeltaMH::update(Random &ran)
{
  int nAccept = 0;

  int g;
  for (g = 0; g < str->G; g++)
    {
      int nOn = 0;
      int q;
      for (q = 0; q < str->Q; q++)
	nOn += str->delta[q][g];

      if ((nOn == 0 || nOn == str->Q) && (oneDelta == 1))   // propose to change all delta's
	{
	  addTry();

	  int oldDelta = str->delta[0][g];
	  int newDelta = 1 - oldDelta;

	  //
	  // compute prior covariance matrix for Delta
	  //

	  vector<vector<double> > var;
	  var.resize(str->Q);
	  int p,q;
	  for (p = 0; p < str->Q; p++)
	    var[p].resize(str->Q);
	  for (p = 0; p < str->Q; p++)
	    {
	      var[p][p] = str->c2 * str->tau2R[p];
	      var[p][p] *= exp(str->b[p] * log(str->sigma2[p][g]));
	    }

	  for (p = 0; p < str->Q; p++)
	    for (q = p + 1; q < str->Q; q++)
	      {
		var[p][q] = str->c2;
		var[p][q] *= str->r[p][q];
		var[p][q] *= sqrt(str->tau2R[p] * str->tau2R[q]);
		var[p][q] *= exp(0.5 * (str->b[q] * log(str->sigma2[q][g]) + str->b[p] * log(str->sigma2[p][g])));
		var[q][p] = var[p][q];
	      }

	  //
	  // define prior mean
	  //

	  vector<double> mean(str->Q);
	  for (p = 0; p < str->Q; p++)
	    mean[p] = 0.0;

	  //
	  // Update parameters for delta[q][g] = 1
	  //

	  vector<vector<double> > var1;
	  inverse(var,var1);
	  vector<double> mean1(str->Q);
	  for (p = 0; p < str->Q; p++)
	    mean1[p] = mean[p];

	  int s;
	  for (q = 0; q < str->Q; q++)
	    {
	      double v0 = str->sigma2[q][g] * str->phi[q][g];
	      double v1 = str->sigma2[q][g] / str->phi[q][g];

	      for (s = 0; s < str->S[q]; s++)
		{
		  double variance = str->psi[q][s] == 0 ? v0 : v1;
		  var1[q][q] += 1.0 / variance;
		  mean1[q] += (2.0 * str->psi[q][s] - 1.0) * (str->x[q][g][s] - str->nu[q][g]) / variance;
		}
	    }

	  vector<vector<double> > var1Inv;
	  inverse(var1,var1Inv);
	  vector<double> mean1Mult;
	  matrixMult(var1Inv,mean1,mean1Mult);

	  //
	  // Draw new values for Delta[g]
	  //

	  vector<double> newDDelta(str->Q);
	  if (newDelta == 0)
	    newDDelta = ran.MultiGaussian(var,mean);
	  else
	    newDDelta = ran.MultiGaussian(var1Inv,mean1Mult);

	  vector<double> oldDDelta(str->Q);
	  for (q = 0; q < str->Q; q++)
	    oldDDelta[q] = str->Delta[q][g];

	  //
	  // Compute acceptance probability
	  //

	  double pot = 0.0;
	  if (newDelta == 0)
	    {
	      pot -= ran.PotentialMultiGaussian(var,mean,newDDelta);
	      pot += ran.PotentialMultiGaussian(var1Inv,mean1Mult,oldDDelta);
	    }
	  else
	    {
	      pot -= ran.PotentialMultiGaussian(var1Inv,mean1Mult,newDDelta);
	      pot += ran.PotentialMultiGaussian(var,mean,oldDDelta);
	    }

	  pot -= model[g]->potential(ran);

	  for (q = 0; q < str->Q; q++)
	    str->delta[q][g] = newDelta;
	  for (q = 0; q < str->Q; q++)
	    str->Delta[q][g] = newDDelta[q];

	  pot += model[g]->potential(ran);

	  for (q = 0; q < str->Q; q++)
	    str->delta[q][g] = oldDelta;
	  for (q = 0; q < str->Q; q++)
	    str->Delta[q][g] = oldDDelta[q];

	  if (ran.Unif01() <= exp(- pot))
	    {
	      addAccept();
	      nAccept++;
	      for (q = 0; q < str->Q; q++)
		str->delta[q][g] = newDelta;
	      for (q = 0; q < str->Q; q++)
		str->Delta[q][g] = newDDelta[q];
	    }
	}
      else if (oneDelta == 0)  // propose to change one delta at a time
	{
	  int q = (int) (ran.Unif01() * str->Q);
	  if (q >= str->Q) q = str->Q - 1;
	  addTry();

	  int oldDelta = str->delta[q][g];
	  int newDelta = 1 - oldDelta;

	  //
	  // compute prior covariance matrix for Delta
	  //

	  vector<vector<double> > var;
	  var.resize(str->Q);
	  int pp,qq;
	  for (pp = 0; pp < str->Q; pp++)
	    var[pp].resize(str->Q);
	  for (pp = 0; pp < str->Q; pp++)
	    {
	      var[pp][pp] = str->c2 * str->tau2R[pp];
	      var[pp][pp] *= exp(str->b[pp] * log(str->sigma2[pp][g]));
	    }

	  for (pp = 0; pp < str->Q; pp++)
	    for (qq = pp + 1; qq < str->Q; qq++)
	      {
		var[pp][qq] = str->c2;
		var[pp][qq] *= str->r[pp][qq];
		var[pp][qq] *= sqrt(str->tau2R[pp] * str->tau2R[qq]);
		var[pp][qq] *= exp(0.5 * (str->b[qq] * log(str->sigma2[qq][g]) + str->b[pp] * log(str->sigma2[pp][g])));
		var[qq][pp] = var[pp][qq];
	      }

	  //
	  // define prior mean
	  //

	  vector<double> mean(str->Q);
	  for (pp = 0; pp < str->Q; pp++)
	    mean[pp] = 0.0;

	  //
	  // Update parameters for delta[qq][g]=1 for qq != q
	  //

	  vector<vector<double> > var0;
	  inverse(var,var0);
	  vector<double> mean0(str->Q);
	  for (pp = 0; pp < str->Q; pp++)
	    mean0[pp] = mean[pp];

	  int s;
	  for (qq = 0; qq < str->Q; qq++)
	    {
	      if (qq != q && str->delta[qq][g] == 1)
		{
		  double v0 = str->sigma2[qq][g] * str->phi[qq][g];
		  double v1 = str->sigma2[qq][g] / str->phi[qq][g];

		  for (s = 0; s < str->S[qq]; s++)
		    {
		      double variance = str->psi[qq][s] == 0 ? v0 : v1;
		      var0[qq][qq] += 1.0 / variance;
		      mean0[qq] += (2.0 * str->psi[qq][s] - 1.0) * (str->x[qq][g][s] - str->nu[qq][g]) / variance;
		    }
		}
	    }

	  vector<vector<double> > var0Inv;
	  inverse(var0,var0Inv);
	  vector<double> mean0Mult;
	  matrixMult(var0Inv,mean0,mean0Mult);

	  //
	  // Update parameters for delta[q][g]=1
	  //

	  double v0 = str->sigma2[q][g] * str->phi[q][g];
	  double v1 = str->sigma2[q][g] / str->phi[q][g];

	  for (s = 0; s < str->S[q]; s++)
	    {
	      double variance = str->psi[q][s] == 0 ? v0 : v1;
	      var0[q][q] += 1.0 / variance;
	      mean0[q] += (2.0 * str->psi[q][s] - 1.0) * (str->x[q][g][s] - str->nu[q][g]) / variance;
	    }

	  vector<vector<double> > var1Inv;
	  inverse(var0,var1Inv);
	  vector<double> mean1Mult;
	  matrixMult(var1Inv,mean0,mean1Mult);

	  //
	  // Draw new values for Delta[g]

	  vector<double> newDDelta(str->Q);
	  if (newDelta == 0)
	    newDDelta = ran.MultiGaussian(var0Inv,mean0Mult);
	  else
	    newDDelta = ran.MultiGaussian(var1Inv,mean1Mult);

	  vector<double> oldDDelta(str->Q);
	  for (qq = 0; qq < str->Q; qq++)
	    oldDDelta[qq] = str->Delta[qq][g];

	  //
	  // compute acceptance probability
	  //

	  double pot = 0.0;
	  if (newDelta == 0)
	    {
	      pot -= ran.PotentialMultiGaussian(var0Inv,mean0Mult,newDDelta);
	      pot += ran.PotentialMultiGaussian(var1Inv,mean1Mult,oldDDelta);
	    }
	  else
	    {
	      pot -= ran.PotentialMultiGaussian(var1Inv,mean1Mult,newDDelta);
	      pot += ran.PotentialMultiGaussian(var0Inv,mean0Mult,oldDDelta);
	    }

	  pot -= model[g]->potential(ran);

	  str->delta[q][g] = newDelta;
	  for (qq = 0; qq < str->Q; qq++)
	    str->Delta[qq][g] = newDDelta[qq];

	  pot += model[g]->potential(ran);

	  str->delta[q][g] = oldDelta;
	  for (qq = 0; qq < str->Q; qq++)
	    str->Delta[qq][g] = oldDDelta[qq];

	  if (ran.Unif01() <= exp(- pot))
	    {
	      addAccept();
	      nAccept++;
	      str->delta[q][g] = newDelta;
	      for (qq = 0; qq < str->Q; qq++)
		str->Delta[qq][g] = newDDelta[qq];
	    }
	}
//      else  // this should never happen!
//	{
//	  cout << "A gene with different delta's are found even if oneDelta = 1 is set! Something is wrong.\n";
//	  exit(-1);
//	}
    }

  return nAccept;
}



#endif
