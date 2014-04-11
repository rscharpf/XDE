#ifndef UPDATEMULTIPLICATIVEPOSITIVE_H
#define UPDATEMULTIPLICATIVEPOSITIVE_H

#include "Update.h"


class UpdateMultiplicativePositive : public Update
{
 public:

  UpdateMultiplicativePositive(const vector<Potential *> &model,const vector<double *> &variable,double epsilon);
  UpdateMultiplicativePositive(const Potential *model,const vector<double *> &variable,double epsilon);
  UpdateMultiplicativePositive(const Potential *model,double *variable,double epsilon);
  UpdateMultiplicativePositive(const vector<Potential *> &model,const vector<double *> &variable1,
			       const vector<double *> &variable2,double epsilon);
  UpdateMultiplicativePositive(const Potential *model,const vector<double *> &variable1,
			       const vector<double *> &variable2,double epsilon);
  UpdateMultiplicativePositive(const Potential *model,double *variable1,double *variable2,double epsilon);
  virtual ~UpdateMultiplicativePositive(void);

  int update(Random &ran);
  virtual Update *copy(void) const;

 protected:
  vector<Potential *> model;
  vector<double *> variable1;
  vector<double *> variable2;
};



inline UpdateMultiplicativePositive::UpdateMultiplicativePositive(const vector<Potential *> &model,
								  const vector<double *> &variable,
								  double epsilon) : Update(epsilon)
{
  if (model.size() != 1 && model.size() != variable.size())
    {
      cout << "ERROR: Internal error! Function \"";
      cout << "UpdateMultiplicativePositive::UpdateMultiplicativePositive\" is called with illegal values.\n";
      cout << "Aborting.\n";
      exit(-1);
    }


  this->model.resize(model.size());
  this->variable1.resize(variable.size());
  this->variable2.resize(variable.size());

  int i;
  for (i = 0; i < model.size(); i++)
    this->model[i] = model[i]->copy();

  for (i = 0; i < variable.size(); i++)
    {
      this->variable1[i] = variable[i];
      this->variable2[i] = NULL;
    }


  return;
}



inline UpdateMultiplicativePositive::UpdateMultiplicativePositive(const Potential *model,
								  const vector<double *> &variable,
								  double epsilon) : Update(epsilon)
{
  this->model.resize(1);
  this->model[0] = model->copy();

  this->variable1.resize(variable.size());
  this->variable2.resize(variable.size());

  int i;
  for (i = 0; i < variable.size(); i++)
    {
      this->variable1[i] = variable[i];
      this->variable2[i] = NULL;
    }


  return;
}



inline UpdateMultiplicativePositive::UpdateMultiplicativePositive(const Potential *model,
								  double *variable,
								  double epsilon) : Update(epsilon)
{
  this->model.resize(1);
  this->model[0] = model->copy();

  this->variable1.resize(1);
  this->variable2.resize(1);
  this->variable1[0] = variable;
  this->variable2[0] = NULL;

  return;
}



inline UpdateMultiplicativePositive::UpdateMultiplicativePositive(const vector<Potential *> &model,
								  const vector<double *> &variable1,
								  const vector<double *> &variable2,
								  double epsilon) : Update(epsilon)
{
  if (model.size() != 1 && (model.size() != variable1.size() || model.size() != variable2.size()))
    {
      cout << "ERROR: Internal error! Function \"";
      cout << "UpdateMultiplicativePositive::UpdateMultiplicativePositive\" is called with illegal values.\n";
      cout << "Aborting.\n";
      exit(-1);
    }
      
  
  this->model.resize(model.size());
  this->variable1.resize(variable1.size());
  this->variable2.resize(variable2.size());
  
  int i;
  for (i = 0; i < model.size(); i++)
    this->model[i] = model[i]->copy();
  
  for (i = 0; i < variable1.size(); i++)
    this->variable1[i] = variable1[i];

  for (i = 0; i < variable2.size(); i++)
    this->variable2[i] = variable2[i];

  return;
}



inline UpdateMultiplicativePositive::UpdateMultiplicativePositive(const Potential *model,
								  const vector<double *> &variable1,
								  const vector<double *> &variable2,
								  double epsilon) : Update(epsilon)
{
  if (variable1.size() != variable2.size())
    {
      cout << "ERROR: Internal error! Function \"";
      cout << "UpdateMultiplicativePositive::UpdateMultiplicativePositive\" is called with illegal values.\n";
      cout << "Aborting.\n";
      exit(-1);
    }


  this->model.resize(1);
  this->model[0] = model->copy();

  this->variable1.resize(variable1.size());
  this->variable2.resize(variable2.size());

  int i;
  for (i = 0; i < variable1.size(); i++)
    this->variable1[i] = variable1[i];
    
  for (i = 0; i < variable2.size(); i++)
    this->variable2[i] = variable2[i];
    
  return;
}



inline UpdateMultiplicativePositive::UpdateMultiplicativePositive(const Potential *model,
								  double *variable1,double *variable2,
								  double epsilon) : Update(epsilon)
{
  this->model.resize(1);
  this->model[0] = model->copy();

  this->variable1.resize(1);
  this->variable2.resize(1);
  this->variable1[0] = variable1;
  this->variable2[0] = variable2;

  return;
}



inline UpdateMultiplicativePositive::~UpdateMultiplicativePositive(void)
{
  return;
}




inline Update *UpdateMultiplicativePositive::copy(void) const
{
  Update *u = new UpdateMultiplicativePositive(model,variable1,variable2,epsilon);

  return u;
}







inline int UpdateMultiplicativePositive::update(Random &ran)
{
  int nAccept = 0;

  int i;
  for (i = 0; i < variable1.size(); i++)
    {
      addTry();
      double pot = 0.0;

      double upper = 1.0 + epsilon;
      double lower = 1.0 / upper;
      
      double valueOld1 = *(variable1[i]);
      double valueOld2 = 0.0;
      if (variable2[i] != NULL) valueOld2 = *(variable2[i]);
      double u = lower + (upper - lower) * ran.Unif01();
      double valueNew1 = valueOld1 * u;
      double valueNew2 = valueOld2 * u;

      if (variable2[i] == NULL)
	pot += - log(1.0 / u);
      else
	pot += - log(1.0);

      if (model.size() == 1)
	pot -= model[0]->potential(ran);
      else
	pot -= model[i]->potential(ran);

      *(variable1[i]) = valueNew1;
      if (variable2[i] != NULL) *(variable2[i]) = valueNew2;

      if (model.size() == 1)
	pot += model[0]->potential(ran);
      else
	pot += model[i]->potential(ran);

      *(variable1[i]) = valueOld1;
      if (variable2[i] != NULL) *(variable2[i]) = valueOld2;

      if (ran.Unif01() <= exp(- pot))
	{
	  *(variable1[i]) = valueNew1;
	  if (variable2[i] != NULL) *(variable2[i]) = valueNew2;
	  addAccept();
	  nAccept++;
	}
    }


  return nAccept;
}



#endif
