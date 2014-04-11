#ifndef UPDATERMH_H
#define UPDATERMH_H

#include "Update.h"
#include "LinAlg.h"
#include "Cholesky.h"

class UpdateRMH : public Update
{
 public:

  UpdateRMH(Structure *str,const Potential *model,double epsilon);
  ~UpdateRMH(void);

  int update(Random &ran);
  Update *copy(void) const;

 private:
  Structure *str;
  Potential *model;
};



inline UpdateRMH::UpdateRMH(Structure *str,const Potential *model,double epsilon) : Update(epsilon)
{
  this->str = str;
  this->model = model->copy();

  return;
}



inline UpdateRMH::~UpdateRMH(void)
{
  delete model;

  return;
}


inline Update *UpdateRMH::copy(void) const
{
  Update *u = new UpdateRMH(str,model,epsilon);

  return u;
} 



inline int UpdateRMH::update(Random &ran)
{
  int nAccept = 0;

  addTry();

  vector<vector<double> > oldR;
  vector<vector<double> > newR;
  oldR.resize(str->Q);
  newR.resize(str->Q);
  int p,q;
  for (p = 0; p < str->Q; p++)
    {
      oldR[p].resize(str->Q);
      newR[p].resize(str->Q);
      for (q = 0; q < str->Q; q++)
	oldR[p][q] = str->r[p][q];
    }
  
  vector<vector<double> > T(str->Q);
  for (p = 0; p < str->Q; p++)
    T[p].resize(str->Q);
  if (ran.Unif01() <= 0.5)
    T = ran.CorrelationStandardInverseWishart(str->Q,str->nuR);
  else
    {
      int m = str->Q;
      double deltaRho = -1.0 / (m - 1.0) + (1.0 + 1.0 / (m - 1.0)) * ran.Unif01();
      for (p = 0; p < str->Q; p++)
	T[p][p] = 1.0;
      for (p = 0; p < str->Q; p++)
	for (q = p+1; q < str->Q; q++)
	  {
	    T[p][q] = deltaRho;
	    T[q][p] = deltaRho;
	  }
    }
  
  
  double pot = 0.0;
  if (ran.Unif01() <= 0.5)
    {
      for (p = 0; p < str->Q; p++)
	for (q = 0; q < str->Q; q++)
	  {
	    if (p != q) 
	      newR[p][q] = (1.0 - epsilon) * oldR[p][q] + epsilon * T[p][q];
	    else
	      newR[p][q] = 1.0;
	  }
      pot += - str->Q * (str->Q - 1.0) * log(1.0 - epsilon) / 2.0;
    }
  else
    {
      for (p = 0; p < str->Q; p++)
	for (q = 0; q < str->Q; q++)
	  {
	    if (p != q) 
	      newR[p][q] = oldR[p][q] / (1.0 - epsilon) - epsilon * T[p][q] / (1.0 - epsilon);
	    else
	      newR[p][q] = 1.0;
	  }
      pot += str->Q * (str->Q - 1) * log(1.0 - epsilon) / 2.0;
    }
  
  int err = 0;
  Matrix newTry(str->Q,str->Q);
  for (p = 0; p < str->Q; p++)
    for (q = 0; q < str->Q; q++)
      newTry(p+1,q+1) = newR[p][q];
  Cholesky chol(newTry,err);
  if (err == 0)  // if err == 1: proposal is not positive definite
    {
      pot -= model->potential(ran);
      
      for (p = 0; p < str->Q; p++)
	for (q = 0; q < str->Q; q++)
	  str->r[p][q] = newR[p][q];
      
      pot += model->potential(ran);
      
      for (p = 0; p < str->Q; p++)
	for (q = 0; q < str->Q; q++)
	  str->r[p][q] = oldR[p][q];
      
      if (ran.Unif01() <= exp(- pot))
	{
	  nAccept++;
	  addAccept();
	  
	  for (p = 0; p < str->Q; p++)
	    for (q = 0; q < str->Q; q++)
	      str->r[p][q] = newR[p][q];
	}
    }

  return nAccept;
}




#endif
