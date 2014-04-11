#ifndef REPORTTAU2Rho_H
#define REPORTTAU2Rho_H

#include "Report.h"


class ReportTau2Rho : public Report
{
 public:

  ReportTau2Rho(const string &filename);
  ReportTau2Rho(double *value);
  ~ReportTau2Rho(void);

  void report(const Structure *str);

 private:
  int writeToFile;

  double *value;
  int nr;
};



inline ReportTau2Rho::ReportTau2Rho(const string &filename) : Report(filename)
{
  writeToFile = 1;

  return;
}


inline ReportTau2Rho::ReportTau2Rho(double *value) : Report()
{
  writeToFile = 0;

  this->value = value;
  nr = 0;

  return;
}



inline ReportTau2Rho::~ReportTau2Rho(void)
{
  return;
}


inline void ReportTau2Rho::report(const Structure *str)
{
  if (writeToFile)
    {
      int q;
      for (q = 0; q < str->Q; q++)
	out << str->tau2Rho[q] << " ";
      
      out << "\n";
      out.flush();
    }
  else
    {
      int q;
      for (q = 0; q < str->Q; q++)
	{
	  value[nr] = str->tau2Rho[q];
	  nr++;
	}
    }

  return;
}


#endif
