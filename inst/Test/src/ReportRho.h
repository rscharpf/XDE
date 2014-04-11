#ifndef REPORTRHO_H
#define REPORTRHO_H

#include "Report.h"


class ReportRho : public Report
{
 public:

  ReportRho(const string &filename);
  ReportRho(double *value);
  ~ReportRho(void);

  void report(const Structure *str);

 private:
  int writeToFile;

  double *value;
  int nr;
};



inline ReportRho::ReportRho(const string &filename) : Report(filename)
{
  writeToFile = 1;

  return;
}



inline ReportRho::ReportRho(double *value) : Report()
{
  writeToFile = 0;

  this->value = value;
  nr = 0;

  return;
}



inline ReportRho::~ReportRho(void)
{
  return;
}


inline void ReportRho::report(const Structure *str)
{
  if (writeToFile)
    {
      int i,j;
      for (i = 0; i < str->Q; i++)
	for (j = i+1; j < str->Q; j++)
	  out << str->rho[i][j] << " ";
      
      out << "\n";
      out.flush();
    }
  else
    {
      int i,j;
      for (i = 0; i < str->Q; i++)
	for (j = i+1; j < str->Q; j++)
	  {
	    value[nr] = str->rho[i][j];
	    nr++;
	  }
    }

  return;
}


#endif
