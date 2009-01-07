#ifndef REPORTTAU2_H
#define REPORTTAU2_H

#include "Report.h"


class ReportTau2 : public Report
{
 public:

  ReportTau2(const string &filename);
  ReportTau2(double *value);
  ~ReportTau2(void);

  void report(const Structure *str);

 private:
  int writeToFile;

  double *value;
  int nr;
};



inline ReportTau2::ReportTau2(const string &filename) : Report(filename)
{
  writeToFile = 1;

  return;
}


inline ReportTau2::ReportTau2(double *value) : Report()
{
  writeToFile = 0;

  this->value = value;
  nr = 0;

  return;
}



inline ReportTau2::~ReportTau2(void)
{
  return;
}


inline void ReportTau2::report(const Structure *str)
{
  if (writeToFile)
    {
      int q;
      for (q = 0; q < str->Q; q++)
	out << str->tau2[q] << " ";
      
      out << "\n";
      out.flush();
    }
  else
    {
      int q;
      for (q = 0; q < str->Q; q++)
	{
	  value[nr] = str->tau2[q];
	  nr++;
	}
    }

  return;
}


#endif
