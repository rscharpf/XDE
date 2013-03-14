#ifndef UPDATETAU2RMH_H
#define UPDATETAU2RMH_H

#include "UpdateMultiplicativePositive.h"


class UpdateTau2RMH : public Update
{
 public:

  UpdateTau2RMH(Structure *str,const Potential *model,double epsilon);
  ~UpdateTau2RMH(void);

  int update(Random &ran);
  Update *copy(void) const;

 private:
  Structure *str;
  Potential *model;
};



inline UpdateTau2RMH::UpdateTau2RMH(Structure *str,const Potential *model,double epsilon) : Update(epsilon)
{
  this->str = str;
  this->model = model->copy();
  
  return;
}



inline UpdateTau2RMH::~UpdateTau2RMH(void)
{
  delete model;

  return;
}


inline Update *UpdateTau2RMH::copy(void) const
{
  Update *u = new UpdateTau2RMH(str,model,epsilon);

  return u;
}




inline int UpdateTau2RMH::update(Random &ran)
{
  int nAccept = 0;

  if (str->Q > 1)
    {
      int q = (int) (str->Q * ran.Unif01());
      int p = (int) ((str->Q - 1) * ran.Unif01());
      if (p >= q) p++;
      
      double upper = 1.0 + epsilon;
      double lower = 1.0 / upper;
      
      double u = lower + (upper - lower) * ran.Unif01();
      vector<double> oldValues(str->Q);
      vector<double> newValues(str->Q);
      
      int i;
      for (i = 0; i < str->Q; i++)
	{
	  oldValues[i] = str->tau2R[i];
	  newValues[i] = str->tau2R[i];
	}
      
      newValues[q] *= u;
      newValues[p] /= u;
      
      double prod = 1.0;
      for (i = 0; i < str->Q; i++)
	prod *= newValues[i];
      
      prod = exp(log(prod) / str->Q);
      for (i = 0; i < str->Q; i++)
	newValues[i] /= prod;
      
      
      
      double pot = - log(1.0 / (u * u));
      pot -= model->potential(ran);
      
      for (i = 0; i < str->Q; i++)
	str->tau2R[i] = newValues[i];
      
      pot += model->potential(ran);
      for (i = 0; i < str->Q; i++)
	str->tau2R[i] = oldValues[i];
      
      addTry();
      if (ran.Unif01() <= exp(- pot))
	{
	  for (i = 0; i < str->Q; i++)
	    str->tau2R[i] = newValues[i];
	  addAccept();
	  nAccept++;
	}
    }
  
  return nAccept;
}


#endif
