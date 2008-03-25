#ifndef REPORTDIFFEXPRESSED_H
#define REPORTDIFFEXPRESSED_H

#include "Structure.h"
#include "Report.h"


class ReportDiffexpressed : public Report
{
 public:

  ReportDiffexpressed(const string &filename,Structure *str);
  ReportDiffexpressed(double *value,Structure *str);
  ~ReportDiffexpressed(void);

  void report(const Structure *str);
  void update(const Structure *str);

 private:
  string filename;
  int writeToFile;

  double *value;

  int nSample;
  vector<vector<int> > nOn;
};



inline ReportDiffexpressed::ReportDiffexpressed(const string &filename,Structure *str) : Report(filename)
{
  this->filename = filename;
  writeToFile = 1;
  file = 0;

  nSample = 0;
  nOn.resize(str->G);
  int g;
  for (g = 0; g < str->G; g++)
    nOn[g].resize(3);

  return;
}


inline ReportDiffexpressed::ReportDiffexpressed(double *value,Structure *str) : Report()
{
  writeToFile = 0;
  
  nSample = 0;
  this->value = value;
  
  nOn.resize(str->G);
  int g;
  for (g = 0; g < str->G; g++)
    nOn[g].resize(3);

  return;
}



inline ReportDiffexpressed::~ReportDiffexpressed(void)
{
  return;
}


inline void ReportDiffexpressed::update(const Structure *str)
{
  return;
}



inline void ReportDiffexpressed::report(const Structure *str)
{  
  /* what should report here now? */

  return;
}


#endif
