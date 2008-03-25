#ifndef POTENTIALDELTA_H
#define POTENTIALDELTA_H

#include "Potential.h"
#include "Structure.h"
#include "Random.h"


class PotentialDelta : public Potential
{
 public:

  PotentialDelta(const Structure *str,int oneDelta);
  virtual ~PotentialDelta(void);

  double potential(Random &ran) const;
  Potential *copy(void) const;

 private:
  const Structure *str;
  int oneDelta;
};



inline PotentialDelta::PotentialDelta(const Structure *str,int oneDelta) : Potential()
{
  this->str = str;
  this->oneDelta = oneDelta;

  return;
}



inline PotentialDelta::~PotentialDelta(void)
{
  return;
}


inline Potential *PotentialDelta::copy(void) const
{
  Potential *pp = new PotentialDelta(str,oneDelta);

  return pp;
}


inline double PotentialDelta::potential(Random &ran) const
{
  double pot = 0.0;
  
  if(oneDelta == 0)
    {
      int q,g;
      for (q = 0; q < str->Q; q++)
	for (g = 0; g < str->G; g++)
	  {
	    if (str->delta[q][g] == 1)
	      pot += - log(str->xi[q]);
	    else
	      pot += - log(1.0 - str->xi[q]);
	  }
    }
  else
    {
      int g;
      for (g = 0; g < str->G; g++)
	{
	  if (str->delta[0][g] == 1)
	    pot += - log(str->xi[0]);
	  else
	    pot += - log(1.0 - str->xi[0]);
	}
    }
  
  return pot;
}



#endif
