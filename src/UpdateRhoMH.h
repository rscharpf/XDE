#ifndef UPDATERHOMH_H
#define UPDATERHOMH_H

#include "Update.h"


class UpdateRhoMH : public Update
{
 public:

  UpdateRhoMH(Structure *str,const Potential *model,double epsilon);
  ~UpdateRhoMH(void);

  int update(Random &ran);
  Update *copy(void) const;

 private:
  Structure *str;
  Potential *model;
};



inline UpdateRhoMH::UpdateRhoMH(Structure *str,const Potential *model,double epsilon) : Update(epsilon)
{
  this->str = str;
  this->model = model->copy();

  return;
}



inline UpdateRhoMH::~UpdateRhoMH(void)
{
  delete model;

  return;
}


inline Update *UpdateRhoMH::copy(void) const
{
  Update *u = new UpdateRhoMH(str,model,epsilon);

  return u;
} 



inline int UpdateRhoMH::update(Random &ran)
{
  int nAccept = 0;

  addTry();

  vector<vector<double> > oldRho;
  vector<vector<double> > newRho;
  oldRho.resize(str->Q);
  newRho.resize(str->Q);
  int p,q;
  for (p = 0; p < str->Q; p++)
    {
      oldRho[p].resize(str->Q);
      newRho[p].resize(str->Q);
      for (q = 0; q < str->Q; q++)
	oldRho[p][q] = str->rho[p][q];
    }

  vector<vector<double> > T(str->Q);
  for (p = 0; p < str->Q; p++)
    T[p].resize(str->Q);
  T = ran.CorrelationStandardInverseWishart(str->Q,str->nuRho);

  double pot = 0.0;
  if (ran.Unif01() <= 0.5)
    {
      for (p = 0; p < str->Q; p++)
	for (q = 0; q < str->Q; q++)
	  {
	    if (p != q) 
	      newRho[p][q] = (1.0 - epsilon) * oldRho[p][q] + epsilon * T[p][q];
	    else
	      newRho[p][q] = 1.0;
	  }
      pot += - str->Q * (str->Q - 1.0) * log(1.0 - epsilon) / 2.0;
    }
  else
    {
      for (p = 0; p < str->Q; p++)
	for (q = 0; q < str->Q; q++)
	  {
	    if (p != q) 
	      newRho[p][q] = oldRho[p][q] / (1.0 - epsilon) - epsilon * T[p][q] / (1.0 - epsilon);
	    else
	      newRho[p][q] = 1.0;
	  }
      pot += str->Q * (str->Q - 1) * log(1.0 - epsilon) / 2.0;
    }

  pot -= model->potential(ran);

  for (p = 0; p < str->Q; p++)
    for (q = 0; q < str->Q; q++)
      str->rho[p][q] = newRho[p][q];
  
  pot += model->potential(ran);

  for (p = 0; p < str->Q; p++)
    for (q = 0; q < str->Q; q++)
      str->rho[p][q] = oldRho[p][q];

  if (ran.Unif01() <= exp(- pot))
    {
      nAccept++;
      addAccept();
      
      for (p = 0; p < str->Q; p++)
	for (q = 0; q < str->Q; q++)
	  str->rho[p][q] = newRho[p][q];
    }

  return nAccept;
}




#endif
