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
      int nOn = 0;
      int q;
      for (q = 0; q < str->Q; q++)
	nOn += str->delta[q][g];

      if (nOn >= 1)
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
		varInv[p][q] *= sqrt(str->tau2R[p] * str->tau2R[q]);
		varInv[p][q] *= exp(0.5 * (str->b[q] * log(str->sigma2[q][g]) + str->b[p] * log(str->sigma2[p][g])));
		
		varInv[q][p] = varInv[p][q];
	      }

	  //
	  // pick out the elements of Delta and varInv that corresponds to delta_qg=1
	  //
	  
	  vector<double> Delta(nOn);
	  vector<vector<double> > varInvRed;
	  varInvRed.resize(nOn);
	  for (q = 0; q < nOn; q++)
	    varInvRed[q].resize(nOn);

	  int qq = 0;
	  for (q = 0; q < str->Q; q++)
	    {
	      if (str->delta[q][g] == 1)
		{
		  Delta[qq] = str->Delta[q][g];
		  qq++;
		}
	    }

	  qq = 0;
	  for (q = 0; q < str->Q; q++)
	    {
	      if (str->delta[q][g] == 1)
		{
		  int p,pp = 0;
		  for (p = 0; p < str->Q; p++)
		    {
		      if (str->delta[p][g] == 1)
			{
			  varInvRed[qq][pp] = varInv[q][p];
			  pp++;
			}
		    }
		  qq++;
		}
	    }


	  vector<vector<double> > var;
	  inverse(varInvRed,var);
	  
	  double value = quadratic(var,Delta);
	  lambda += 0.5 * value;
	  s += 0.5 * nOn;
	}
    }

  //
  // Draw new value
  //
  
  double newValue;
  if (s > 0.0) // there exist at least one delta == 1
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
      double fmax;
      if (lambda < str->c2Max)
	fmax = exp(-1.0)/lambda;
      else
	fmax = exp(-lambda/str->c2Max) / str->c2Max;
      int nTry = 0;
      int accept = 0;
      do
	{
	  nTry++;
	  newValue = str->c2Max * ran.Unif01();
	  double alpha = (exp(-lambda/newValue)/newValue) / fmax;
	  accept = (ran.Unif01() <= alpha);
	}
      while (accept == 0 && nTry < 100);
      if (nTry == 100) newValue = str->c2;   // propose an unchanged value!
    }
  else if (s == -0.5)
    {
      double fmax;
      if (lambda < 0.5 * str->c2Max)
	fmax = exp(-0.5)/sqrt(2*lambda);
      else
	fmax = exp(-lambda/str->c2Max) / sqrt(str->c2Max);
      int nTry = 0;
      int accept = 0;
      do 
	{
	  nTry++;
	  newValue = str->c2Max * ran.Unif01();
	  double alpha = (exp(-lambda/newValue)/sqrt(newValue)) / fmax;
	  accept = (ran.Unif01() < alpha);
	}
      while (accept == 0 && newValue < 100);
      if (nTry == 100) newValue = str->c2;   // propose an unchanged value!
    }
  else
    newValue = str->c2Max * ran.Unif01();

  //
  // Set new value
  //
  
  str->c2 = newValue;

  //
  // draw new Delta values for genes with delta_qg == 0
  //

  for (g = 0; g < str->G; g++)
    {
      int nOn = 0;
      int q;
      for (q = 0; q < str->Q; q++)
	nOn += str->delta[q][g];

      if (nOn < str->Q)
	{
	  vector<vector<double> > var;
	  var.resize(str->Q);
	  int p,q;
	  for (p = 0; p < str->Q; p++) var[p].resize(str->Q);
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
	  
	  vector<double> mean(str->Q);
	  for (p = 0; p < str->Q; p++) mean[p] = 0.0;

	  if (nOn == 0)  // draw all the Delta's
	    {
	      vector<double> newDDelta = ran.MultiGaussian(var,mean);
	      for (p = 0; p < str->Q; p++)
		str->Delta[p][g] = newDDelta[p];
	    }
	  else  // draw some Delta's conditional on the rest
	    {
	      vector<double> Delta1(nOn);
	      vector<double> Delta2(str->Q - nOn);
	      vector<double> mu1(nOn);
	      vector<double> mu2(str->Q - nOn);
	      vector<vector<double> > sigma11(nOn);
	      int q;
	      for (q = 0; q < nOn; q++)
		sigma11[q].resize(nOn);
	      vector<vector<double> > sigma12(nOn);
	      for (q = 0; q < nOn; q++)
		sigma12[q].resize(str->Q - nOn);
	      vector<vector<double> > sigma21(str->Q - nOn);
	      for (q = 0; q < str->Q - nOn; q++)
		sigma21[q].resize(nOn);
	      vector<vector<double> > sigma22(str->Q - nOn);
	      for (q = 0; q < str->Q - nOn; q++)
		sigma22[q].resize(str->Q - nOn);

	      int q1 = 0;
	      int q2 = 0;
	      for (q = 0; q < str->Q; q++)
		{
		  if (str->delta[q][g] == 1)
		    {
		      Delta1[q1] = str->Delta[q][g];
		      mu1[q1] = mean[q];
		      q1++;
		    }
		  else
		    {
		      Delta2[q2] = str->Delta[q][g];
		      mu2[q2] = mean[q];
		      q2++;
		    }
		}
	      
	      q1 = 0;
	      q2 = 0;
	      for (q = 0; q < str->Q; q++)
		{
		  if (str->delta[q][g] == 1)
		    {
		      int p1 = 0;
		      int p2 = 0;
		      int p;
		      for (p = 0; p < str->Q; p++)
			{
			  if (str->delta[p][g] == 1)
			    {
			      sigma11[q1][p1] = var[q][p];
			      p1++;
			    }
			  else
			    {
			      sigma12[q1][p2] = var[q][p];
			      p2++;
			    }
			}
		      q1++;
		    }
		  else
		    {
		      int p1 = 0;
		      int p2 = 0;
		      int p;
		      for (p = 0; p < str->Q; p++)
			{
			  if (str->delta[p][g] == 1)
			    {
			      sigma21[q2][p1] = var[q][p];
			      p1++;
			    }
			  else
			    {
			      sigma22[q2][p2] = var[q][p];
			      p2++;
			    }
			}
		      q2++;
		    }
		}
	      
	      vector<vector<double> > sigma11Inv;
	      inverse(sigma11,sigma11Inv);
	      vector<double> diff1(Delta1);
	      for (q = 0; q < nOn; q++)
		diff1[q] -= mu1[q];
	      vector<vector<double> > sigma21sigma11Inv;
	      matrixMult(sigma21,sigma11Inv,sigma21sigma11Inv);
	      vector<double> correct1;
	      matrixMult(sigma21sigma11Inv,diff1,correct1);
	      vector<double> muCond(mu2);
	      for (q = 0; q < str->Q - nOn; q++)
		muCond[q] += correct1[q];
	      vector<vector<double> > correctMatrix;
	      matrixMult(sigma21sigma11Inv,sigma12,correctMatrix);
	      vector<vector<double> > sigmaCond(str->Q - nOn);
	      for (q = 0; q < str->Q - nOn; q++)
		{
		  sigmaCond[q].resize(str->Q - nOn);
		  int p;
		  for (p = 0; p < str->Q - nOn; p++)
		    sigmaCond[q][p] = sigma22[q][p] - correctMatrix[q][p];
		}

	      vector<double> newDDelta = ran.MultiGaussian(sigmaCond,muCond);
	      int qq = 0;
	      for (q = 0; q < str->Q; q++)
		{
		  if (str->delta[q][g] == 0)
		    {
		      str->Delta[q][g] = newDDelta[qq];
		      qq++;
		    }
		}
	    }
	}
    }
  
  
  addTry();
  addAccept();
  nAccept++;
  
  return nAccept;
}


#endif
