#ifndef REPORTTHETA_H
#define REPORTTHETA_H

#include "Report.h"


class ReportTheta : public Report
{
 public:

  ReportTheta(const string &filename);
  ReportTheta(double *value);
  ~ReportTheta(void);

  void report(const Structure *str);

 private:
  int writeToFile;

  double *value;
  int nr;
};



inline ReportTheta::ReportTheta(const string &filename) : Report(filename)
{
  writeToFile = 1;

  return;
}



inline ReportTheta::ReportTheta(double *value) : Report()
{
  writeToFile = 0;

  this->value = value;
  nr = 0;

  return;
}


inline ReportTheta::~ReportTheta(void)
{
  return;
}


inline void ReportTheta::report(const Structure *str)
{
  if (writeToFile)
    {
      int q;
      for (q = 0; q < str->Q; q++)
	out << str->theta[q] << " ";
      
      out << "\n";
      out.flush();
    }
  else
    {
      int q;
      for (q = 0; q < str->Q; q++)
	{
	  value[nr] = str->theta[q];
	  nr++;
	}
    }

  return;
}


#endif
