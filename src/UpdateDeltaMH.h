#ifndef UPDATEDELTAGIBBS_H
#define UPDATEDELTAGIBBS_H

#include "Update.h"
#include "PotentialDelta.h"
#include "PotentialDDeltag.h"
#include "PotentialXqg.h"


class UpdateDeltaMH : public Update
{
 public:

  UpdateDeltaMH(Structure *str);
  ~UpdateDeltaMH(void);

  int update(Random &ran);
  Update *copy(void) const;

 private:
  Structure *str;

  vector<Potential *> model;
};



inline UpdateDeltaMH::UpdateDeltaMH(Structure *str) : Update(0.0)
{
  this->str = str;

  int g;
  for (g = 0; g < str->G; g++)
    {
      vector<Potential *> term;
      term.push_back(new PotentialDelta(str));
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
      addTry();

      int oldDelta = str->delta[g];
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

      vector<double> mean(str->Q);
      for (p = 0; p < str->Q; p++) 
	mean[p] = 0.0;

      //
      // Update parameters for delta[g] = 1
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
      
      str->delta[g] = newDelta;
      for (q = 0; q < str->Q; q++)
	str->Delta[q][g] = newDDelta[q];

      pot += model[g]->potential(ran);

      str->delta[g] = oldDelta;
      for (q = 0; q < str->Q; q++)
	str->Delta[q][g] = oldDDelta[q];

      if (ran.Unif01() <= exp(- pot))
	{
	  addAccept();
	  nAccept++;
	  str->delta[g] = newDelta;
	  for (q = 0; q < str->Q; q++)
	    str->Delta[q][g] = newDDelta[q];
	}
    }
  
  return nAccept;
}


#endif
