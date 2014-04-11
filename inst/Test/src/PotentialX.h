#ifndef POTENTIALX_H
#define POTENTIALX_H


#include "Potential.h"
#include "PotentialXqg.h"
#include "Structure.h"
#include "Random.h"

class PotentialX : public Potential
{
 public:

  PotentialX(const Structure *str);
  virtual ~PotentialX(void);

  double potential(Random &ran) const;
  Potential *copy(void) const;

 private:
  const Structure *str;
  vector<Potential *> model;
};



inline PotentialX::PotentialX(const Structure *str) : Potential()
{
  this->str = str;

  int q,g;
  for (q = 0; q < str->Q; q++)
    for (g = 0; g < str->G; g++)
      model.push_back(new PotentialXqg(q,g,str));

  return;
}



inline PotentialX::~PotentialX(void)
{
  int i;
  for (i = 0; i < model.size(); i++)
    delete model[i];

  return;
}


inline Potential *PotentialX::copy(void) const
{
  Potential *pp = new PotentialX(str);

  return pp;
}


inline double PotentialX::potential(Random &ran) const
{
  double pot = 0.0;

  int i;
  for (i = 0; i < model.size(); i++)
    pot += model[i]->potential(ran);

  return pot;
} 



#endif
