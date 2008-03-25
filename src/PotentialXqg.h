#ifndef POTENTIALXQG_H
#define POTENTIALXQG_H

#include "Potential.h"
#include "Structure.h"
#include "Random.h"


class PotentialXqg : public Potential
{
 public:

  PotentialXqg(int q,int g,const Structure *str);
  virtual ~PotentialXqg(void);

  double potential(Random &ran) const;
  Potential *copy(void) const;

 private:
  int q,g;
  const Structure *str;
};



inline PotentialXqg::PotentialXqg(int q,int g,const Structure *str) : Potential()
{
  this->q = q;
  this->g = g;
  this->str = str;

  return;
}



inline PotentialXqg::~PotentialXqg(void)
{
  return;
}


inline Potential *PotentialXqg::copy(void) const
{
  Potential *pp = new PotentialXqg(q,g,str);

  return pp;
}


inline double PotentialXqg::potential(Random &ran) const
{
  double pot = 0.0;
  
  double var0 = str->sigma2[q][g] * str->phi[q][g];
  double var1 = str->sigma2[q][g] / str->phi[q][g];
  double mm = str->nu[q][g];

  if (str->delta[q][g] != 0)
    {
      int s;
      for (s = 0; s < str->S[q]; s++)
	{
	  double mean;
	  double var;
	  if (str->psi[q][s] == 0)
	    {
	      mean = mm - str->Delta[q][g];
	      var = var0;
	    }
	  else
	    {
	      mean = mm + str->Delta[q][g];
	      var = var1;
	    }

	  //	  double diff = str->x[q][g][s] - mean;
	  //	  pot += 0.5 * (diff * diff / var + log(2.0 * PI) + log(var));
	  pot += ran.PotentialGaussian(var,mean,str->x[q][g][s]);
	}
    }
  else
    {
      int s;
      double mean = mm;
      for (s = 0; s < str->S[q]; s++)
	{
	  double var = str->psi[q][s] == 0 ? var0 : var1;
	  
	  //	  double diff = str->x[q][g][s] - mean;
	  //	  pot += 0.5 * (diff * diff / var + log(2.0 * PI) + log(var));
	  pot += ran.PotentialGaussian(var,mean,str->x[q][g][s]);
	}
    }	  

  
  return pot;
}



#endif
