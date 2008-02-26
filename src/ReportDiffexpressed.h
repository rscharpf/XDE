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
  nSample++;

  int g;
  for (g = 0; g < str->G; g++)
    {
      if (str->delta[g] == 1)
	{
	  nOn[g][0] += str->delta[g];
	  int concordant = 1;
	  int v = (str->Delta[0][g] >= 0.0);
	  int q;
	  for (q = 1; q < str->Q; q++)
	    {
	      if ((str->Delta[q][g] >= 0.0) != v)
		concordant = 0;
	    }
	  nOn[g][1] += concordant;
	  nOn[g][2] += 1 - concordant;
	}
    }

  return;
}



inline void ReportDiffexpressed::report(const Structure *str)
{  
  if (writeToFile)
    {
      out.close();
      out.open(filename.c_str());
      int g;
      for (g = 0; g < str->G; g++)
	out << ((double) nOn[g][0]) / ((double) nSample) << " ";
      for (g = 0; g < str->G; g++)
	out << ((double) nOn[g][1]) / ((double) nSample) << " ";
      for (g = 0; g < str->G; g++)
	out << ((double) nOn[g][2]) / ((double) nSample) << " ";
      
      out << "\n";
      out.flush();
    }
  else
    {
      int nr = 0;
      int g;
      for (g = 0; g < str->G; g++)
	{
	  value[nr] = ((double) nOn[g][0]) / ((double) nSample);
	  nr++;
	}
      for (g = 0; g < str->G; g++)
	{
	  value[nr] = ((double) nOn[g][1]) / ((double) nSample);
	  nr++;
	}
      for (g = 0; g < str->G; g++)
	{
	  value[nr] = ((double) nOn[g][2]) / ((double) nSample);
	  nr++;
	}
    }

  return;
}


#endif
