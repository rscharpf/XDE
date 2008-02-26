#ifndef UPDATEXIGIBBS_H
#define UPDATEXIGIBBS_H

#include "Update.h"


class UpdateXiGibbs : public Update
{
 public:

  UpdateXiGibbs(Structure *str,int check,const Potential *model);
  ~UpdateXiGibbs(void);

  int update(Random &ran);
  Update *copy(void) const;

 private:
  Structure *str;
  int check;
  Potential *model;
};



inline UpdateXiGibbs::UpdateXiGibbs(Structure *str,int check,const Potential *model) : Update(0.0)
{
  this->str = str;
  this->check = check;
  this->model = model->copy();

  return;
}



inline UpdateXiGibbs::~UpdateXiGibbs(void)
{
  delete model;

  return;
}


inline Update *UpdateXiGibbs::copy(void) const
{
  Update *u = new UpdateXiGibbs(str,check,model);

  return u;
}



inline int UpdateXiGibbs::update(Random &ran)
{
  int nAccept = 0;

  //
  // set prior parameters
  //

  double alpha = str->alphaXi;
  double beta = str->betaXi;

  //
  // update parameters based on available observations
  //

  int g;
  for (g = 0; g < str->G; g++)
    {
      if (str->delta[g] == 1)
	alpha += 1.0;
      else
	beta += 1.0;
    }

  //
  // Draw new value
  //
  
  double newValue = ran.Beta(alpha,beta);

  //
  // Check acceptance probability
  //
  
  if (check != 0)
    {
      double oldValue = str->xi;
      double pot = - model->potential(ran);
      pot -= ran.PotentialBeta(alpha,beta,newValue);
      
      str->xi = newValue;

      pot += model->potential(ran);
      pot += ran.PotentialBeta(alpha,beta,oldValue);

      str->xi = oldValue;
      
      if (pot >= 1.0e-6 || pot <= -1.0e-6)
	cout << "WARNING: Possible implementation error in UpdateXiGibbs located. Check out!\n\n";
    }
  
  //
  // Set new value
  //
  
  str->xi = newValue;

  addTry();
  addAccept();
  nAccept++;
  
  return nAccept;
}


#endif
