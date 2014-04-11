#ifndef UPDATE_H
#define UPDATE_H

#include "Structure.h"
#include "Potential.h"

class Update
{
 public:
  
  Update(double epsilon);
  virtual ~Update(void);

  virtual int update(Random &ran) = 0;
  virtual Update *copy(void) const = 0;
  virtual void setEpsilon(double epsilon) {this->epsilon = epsilon;};
  double getEpsilon(void) const {return epsilon;};

  void addTry(void) {nTry++;};
  void addAccept(void) {nAccept++;};
  virtual double acceptRate(void) {if (nTry == 0) return 0.0; else return ((double) nAccept) / ((double) nTry);};

 protected:
  double epsilon;

 private:
  int nTry;
  int nAccept;
};



inline Update::Update(double epsilon)
{
  this->epsilon = epsilon;

  nTry = 0;
  nAccept = 0;

  return;
}



inline Update::~Update(void)
{
  return;
}



#endif
