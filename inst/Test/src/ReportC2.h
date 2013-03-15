#ifndef REPORTC2_H
#define REPORTC2_H

#include "Report.h"


class ReportC2 : public Report
{
 public:

  ReportC2(const string &filename);
  ReportC2(double *value);
  ~ReportC2(void);

  void report(const Structure *str);

 private:
  int writeToFile;

  double *value;
  int nr;
};



inline ReportC2::ReportC2(const string &filename) : Report(filename)
{
  writeToFile = 1;

  return;
}



inline ReportC2::ReportC2(double *value) : Report()
{
  writeToFile = 0;

  this->value = value;
  nr = 0;

  return;
}


inline ReportC2::~ReportC2(void)
{
  return;
}


inline void ReportC2::report(const Structure *str)
{
  if (writeToFile)
    {
      out << str->c2 << " ";
      
      out << "\n";
      out.flush();
    }
  else
    {
      value[nr] = str->c2;
      nr++;
    }

  return;
}


#endif
