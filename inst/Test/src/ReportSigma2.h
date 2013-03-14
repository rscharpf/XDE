#ifndef REPORTSIGMA2_H
#define REPORTSIGMA2_H

#include "Report.h"


class ReportSigma2 : public Report
{
 public:

  ReportSigma2(const string &filename);
  ReportSigma2(double *value);
  ~ReportSigma2(void);

  void report(const Structure *str);

 private:
  int writeToFile;

  double *value;
  int nr;
};



inline ReportSigma2::ReportSigma2(const string &filename) : Report(filename)
{
  writeToFile = 1;

  return;
}



inline ReportSigma2::ReportSigma2(double *value) : Report()
{
  writeToFile = 0;

  this->value = value;
  nr = 0;

  return;
}



inline ReportSigma2::~ReportSigma2(void)
{
  return;
}


inline void ReportSigma2::report(const Structure *str)
{
  if (writeToFile)
    {
      int q,g;
      for (g = 0; g < str->G; g++)
	for (q = 0; q < str->Q; q++)
	  out << str->sigma2[q][g] << " ";
      
      out << "\n";
      out.flush();
    }
  else
    {
      int q,g;
      for (g = 0; g < str->G; g++)
	for (q = 0; q < str->Q; q++)
	  {
	    value[nr] = str->sigma2[q][g];
	    nr++;
	  }
    }
  
  return;
}


#endif
