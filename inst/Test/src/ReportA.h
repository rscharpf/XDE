#ifndef REPORTA_H
#define REPORTA_H

#include "Report.h"


class ReportA : public Report
{
 public:

  ReportA(const string &filename);
  ReportA(double *value);
  ~ReportA(void);

  void report(const Structure *str);

 private:
  int writeToFile;

  double *value;
  int nr;
};



inline ReportA::ReportA(const string &filename) : Report(filename)
{
  writeToFile = 1;

  return;
}


inline ReportA::ReportA(double *value) : Report()
{
  writeToFile = 0;

  this->value = value;
  nr = 0;

  return;
}



inline ReportA::~ReportA(void)
{
  return;
}


inline void ReportA::report(const Structure *str)
{
  if (writeToFile)
    {
      int q;
      for (q = 0; q < str->Q; q++)
	out << str->a[q] << " ";
      
      out << "\n";
      out.flush();
    }
  else
    {
      int q;
      for (q = 0; q < str->Q; q++)
	{
	  value[nr] = str->a[q];
	  nr++;
	}
    }

  return;
}


#endif
