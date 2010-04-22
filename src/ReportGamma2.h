#ifndef REPORTGAMMA2_H
#define REPORTGAMMA2_H

#include "Report.h"


class ReportGamma2 : public Report
{
 public:

  ReportGamma2(const string &filename);
  ReportGamma2(double *value);
  ~ReportGamma2(void);

  void report(const Structure *str);

 private:
  int writeToFile;

  double *value;
  int nr;
};



inline ReportGamma2::ReportGamma2(const string &filename) : Report(filename)
{
  writeToFile = 1;

  return;
}


inline ReportGamma2::ReportGamma2(double *value) : Report()
{
  writeToFile = 0;

  this->value = value;
  nr = 0;

  return;
}


inline ReportGamma2::~ReportGamma2(void)
{
  return;
}


inline void ReportGamma2::report(const Structure *str)
{
  if (writeToFile)
    {
      out << str->gamma2 << " ";
      
      out << "\n";
      out.flush();
    }
  else
    {
      value[nr] = str->gamma2;
      nr++;
    }

  return;
}


#endif
