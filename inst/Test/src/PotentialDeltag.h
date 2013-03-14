#ifndef POTENTIALDELTAG_H
#define POTENTIALDELTAG_H

#include "Potential.h"
#include "Structure.h"
#include "Random.h"


class PotentialDeltag : public Potential
{
 public:

  PotentialDeltag(int g,const Structure *str,int oneDelta);
  virtual ~PotentialDeltag(void);

  double potential(Random &ran) const;
  Potential *copy(void) const;

 private:
  int g;
  const Structure *str;
  int oneDelta;
};



inline PotentialDeltag::PotentialDeltag(int g,const Structure *str,int oneDelta) : Potential()
{
  this->g = g;
  this->str = str;
  this->oneDelta = oneDelta;

  return;
}



inline PotentialDeltag::~PotentialDeltag(void)
{
  return;
}


inline Potential *PotentialDeltag::copy(void) const
{
  Potential *pp = new PotentialDeltag(g,str,oneDelta);

  return pp;
}


inline double PotentialDeltag::potential(Random &ran) const
{
  double pot = 0.0;
  
  if(oneDelta == 0)
    {
      int q;
      for (q = 0; q < str->Q; q++)
	{
	  if (str->delta[q][g] == 1)
	    pot += - log(str->xi[q]);
	  else
	    pot += - log(1.0 - str->xi[q]);
	}
    }
  else
    {
      if (str->delta[0][g] == 1)
	pot += - log(str->xi[0]);
      else
	pot += - log(1.0 - str->xi[0]);
    }
  
  return pot;
}



#endif
