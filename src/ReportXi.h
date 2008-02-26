#ifndef REPORTXI_H
#define REPORTXI_H

#include "Report.h"


class ReportXi : public Report
{
 public:

  ReportXi(const string &filename);
  ReportXi(double *value);
  ~ReportXi(void);

  void report(const Structure *str);

 private:
  int writeToFile;

  double *value;
  int nr;
};



inline ReportXi::ReportXi(const string &filename) : Report(filename)
{
  writeToFile = 1;

  return;
}



inline ReportXi::ReportXi(double *value) : Report()
{
  writeToFile = 0;

  this->value = value;
  nr = 0;

  return;
}


inline ReportXi::~ReportXi(void)
{
  return;
}


inline void ReportXi::report(const Structure *str)
{
  if (writeToFile)
    {
      out << str->xi << " ";
      
      out << "\n";
      out.flush();
    }
  else
    {
      value[nr] = str->xi;
      nr++;
    }

  return;
}


#endif
