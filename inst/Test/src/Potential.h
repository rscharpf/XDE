#ifndef POTENTIAL_H
#define POTENTIAL_H

class Random;


class Potential
{
 public:
  
  Potential(void);
  virtual ~Potential(void);

  virtual double potential(Random &ran) const = 0;
  virtual Potential *copy(void) const = 0;

 private:
};



inline Potential::Potential(void)
{
  return;
}



inline Potential::~Potential(void)
{
  return;
} 


#endif
