#ifndef POTENTIALSIGMA2_H
#define POTENTIALSIGMA2_H

#include "Potential.h"
#include "PotentialSigma2qg.h"
#include "Structure.h"
#include "Random.h"


class PotentialSigma2 : public Potential
{
 public:

  PotentialSigma2(const Structure *str);
  virtual ~PotentialSigma2(void);

  double potential(Random &ran) const;
  Potential *copy(void) const;

 private:
  const Structure *str;
  vector<Potential *> model;
};



inline PotentialSigma2::PotentialSigma2(const Structure *str) : Potential()
{
  this->str = str;

  int q,g;
  for (q = 0; q < str->Q; q++)
    for (g = 0; g < str->G; g++)
      model.push_back(new PotentialSigma2qg(q,g,str));

  return;
}



inline PotentialSigma2::~PotentialSigma2(void)
{
  int i;
  for (i = 0; i < model.size(); i++)
    delete model[i];

  return;
}


inline Potential *PotentialSigma2::copy(void) const
{
  Potential *pp = new PotentialSigma2(str);

  return pp;
}


inline double PotentialSigma2::potential(Random &ran) const
{
  double pot = 0.0;

  int i;
  for (i = 0; i < model.size(); i++)
    pot += model[i]->potential(ran);

  return pot;
}



#endif
