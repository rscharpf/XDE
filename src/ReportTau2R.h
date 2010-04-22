#ifndef REPORTTAU2R_H
#define REPORTTAU2R_H

#include "Report.h"


class ReportTau2R : public Report
{
 public:

  ReportTau2R(const string &filename);
  ReportTau2R(double *value);
  ~ReportTau2R(void);

  void report(const Structure *str);

 private:
  int writeToFile;

  double *value;
  int nr;
};



inline ReportTau2R::ReportTau2R(const string &filename) : Report(filename)
{
  writeToFile = 1;

  return;
}


inline ReportTau2R::ReportTau2R(double *value) : Report()
{
  writeToFile = 0;

  this->value = value;
  nr = 0;

  return;
}



inline ReportTau2R::~ReportTau2R(void)
{
  return;
}


inline void ReportTau2R::report(const Structure *str)
{
  if (writeToFile)
    {
      int q;
      for (q = 0; q < str->Q; q++)
	out << str->tau2R[q] << " ";
      
      out << "\n";
      out.flush();
    }
  else
    {
      int q;
      for (q = 0; q < str->Q; q++)
	{
	  value[nr] = str->tau2R[q];
	  nr++;
	}
    }

  return;
}


#endif
