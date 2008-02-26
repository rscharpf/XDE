#include <string>
#include <sstream>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <iostream>
#include <fstream>
#include <stdio.h>

using namespace std;

#define CHECKGIBBS 0

#include "Random.h"
#include "Structure.h"

#include "Report.h"

#include "PotentialZero.h"
#include "PotentialX.h"
#include "PotentialNu.h"
#include "PotentialDDelta.h"
#include "PotentialA.h"
#include "PotentialB.h"
#include "PotentialR.h"
#include "PotentialRho.h"
#include "PotentialDelta.h"
#include "PotentialXi.h"
#include "PotentialSigma2.h"
#include "PotentialPhi.h"
#include "PotentialSum.h"

#include "UpdateNuGibbs.h"
#include "UpdateDDeltaGibbs.h"
#include "UpdateAMH.h"
#include "UpdateBMH.h"
#include "UpdateC2Gibbs.h"
#include "UpdateGamma2Gibbs.h"
#include "UpdateRC2MH.h"
#include "UpdateRhoGamma2MH.h"
#include "UpdateDeltaMH.h"
#include "UpdateXiGibbs.h"
#include "UpdateSigma2MH.h"
#include "UpdateTMH.h"
#include "UpdateLMH.h"
#include "UpdatePhiMH.h"
#include "UpdateThetaMH.h"
#include "UpdateLambdaMH.h"
#include "UpdateTau2MH.h"

#include "ReportDiffexpressed.h"

int main(int argc,char **argv)
{
  cout.setf(ios::scientific,ios::floatfield);

  if (argc != 8)
    {
      cout << "Usage: " << argv[0] << " <seed> <number of iterations> <data file> <parameter file> <initial values file> <update file> <output file>\n\n";
      exit(-1);
    }

  //
  // Initialise data and parameters
  //

  unsigned int seed;
  seed = atoi(argv[1]);
  Random ran(seed);


  int nIt;
  nIt = atoi(argv[2]);


  string filenameData(argv[3]);
  Structure str(filenameData,ran);

  //
  // Set non-default hyper-parameter values
  //

  string filenameParam(argv[4]);
  if (argv[4][0] != '=')
    str.setParameterValues(filenameParam);
  str.printParameterValues();
  
  //
  // Set non-default initial parameter values
  //

  string filenameInitial(argv[5]);
  if (argv[5][0] != '=')
    {
      cout << "Specified initial values:\n\n";
      str.setInitialValues(filenameInitial);
      cout << "\n\n";
    }


  //
  // Define potential functions for each variable
  //

  PotentialX potX(&str);
  PotentialNu potNu(&str);
  PotentialDDelta potDDelta(&str);
  PotentialA potA(&str);
  PotentialB potB(&str);
  PotentialZero potC2;
  PotentialZero potGamma2;
  PotentialR potR(&str);
  PotentialRho potRho(&str);
  PotentialDelta potDelta(&str);
  PotentialXi potXi(&str);
  PotentialSigma2 potSigma2(&str);
  PotentialZero potT;
  PotentialZero potL;
  PotentialPhi potPhi(&str);
  PotentialZero potTheta;
  PotentialZero potLambda;
  PotentialZero potTau2;

  //
  // Define the total potential function
  //

  vector<Potential *> term;
  term.push_back(&potX);
  term.push_back(&potNu);
  term.push_back(&potDDelta);
  term.push_back(&potA);
  term.push_back(&potB);
  term.push_back(&potC2);
  term.push_back(&potGamma2);
  term.push_back(&potR);
  term.push_back(&potRho);
  term.push_back(&potDelta);
  term.push_back(&potXi);
  term.push_back(&potSigma2);
  term.push_back(&potT);
  term.push_back(&potL);
  term.push_back(&potPhi);
  term.push_back(&potTheta);
  term.push_back(&potLambda);
  term.push_back(&potTau2);
  
  PotentialSum potTotal(term);
  
  vector<Potential *> potAll;
  potAll.push_back(&potTotal);
  int i;
  for (i = 0; i < term.size(); i++)
    potAll.push_back(term[i]);
  
  //
  // Define Metropolis-Hastings updates
  //
  
  vector<Update *> update;
  vector<int> nUpdate;
  Potential *potSum = NULL;
  
  //
  // Update for nu
  //
  
  update.push_back(new UpdateNuGibbs(&str,CHECKGIBBS,&potTotal));
  nUpdate.push_back(1);
  
  //
  // Update for Delta
  //
  
  update.push_back(new UpdateDDeltaGibbs(&str,CHECKGIBBS,&potTotal));
  nUpdate.push_back(1);
  
  //
  // Update for a
  //
  
  term.resize(0);
  term.push_back(&potA);
  term.push_back(&potNu);
  potSum = new PotentialSum(term);
  update.push_back(new UpdateAMH(&str,potSum,0.05));
  nUpdate.push_back(20);
  delete potSum;
  potSum = NULL;
  
  //
  // Update for b
  //
  
  term.resize(0);
  term.push_back(&potB);
  term.push_back(&potDDelta);
  potSum = new PotentialSum(term);
  update.push_back(new UpdateBMH(&str,potSum,0.05));
  nUpdate.push_back(20);
  delete potSum;
  potSum = NULL;
  
  //
  // Update for c2
  //
  
  update.push_back(new UpdateC2Gibbs(&str,CHECKGIBBS,&potTotal));
  nUpdate.push_back(1);
  
  //
  // Update for gamma2
  //
  
  update.push_back(new UpdateGamma2Gibbs(&str,CHECKGIBBS,&potTotal));
  nUpdate.push_back(1);
  
  //
  // Update for r
  //
  
  term.resize(0);
  term.push_back(&potR);
  term.push_back(&potC2);
  term.push_back(&potDDelta);
  potSum = new PotentialSum(term);
  update.push_back(new UpdateRC2MH(&str,potSum,0.01));
  nUpdate.push_back(20);
  delete potSum;
  potSum = NULL;
  
  //
  // Update for rho
  //
  
  term.resize(0);
  term.push_back(&potRho);
  term.push_back(&potGamma2);
  term.push_back(&potNu);
  potSum = new PotentialSum(term);
  update.push_back(new UpdateRhoGamma2MH(&str,potSum,0.01));
  nUpdate.push_back(20);
  delete potSum;
  potSum = NULL;
  
  //
  // Update for delta
  //
  
  update.push_back(new UpdateDeltaMH(&str));
  nUpdate.push_back(1);
  
  //
  // Update for xi
  //
  
  update.push_back(new UpdateXiGibbs(&str,CHECKGIBBS,&potTotal));
  nUpdate.push_back(1);
  
  //
  // Update for sigma2
  //
  
  update.push_back(new UpdateSigma2MH(&str,0.5));
  nUpdate.push_back(3);
  
  //
  // Update for t
  //
  
  update.push_back(new UpdateTMH(&str,&potT,0.1));
  nUpdate.push_back(10);
  
  //
  // Update for l
  //
  
  update.push_back(new UpdateLMH(&str,&potL,0.1));
  nUpdate.push_back(10);
  
  //
  // Update for phi
  //
  
  update.push_back(new UpdatePhiMH(&str,0.5));
  nUpdate.push_back(3);
  
  //
  // Update for theta
  //
  
  update.push_back(new UpdateThetaMH(&str,&potTheta,0.1));
  nUpdate.push_back(10);
  
  //
  // Update for lambda
  //
  
  update.push_back(new UpdateLambdaMH(&str,&potLambda,0.1));
  nUpdate.push_back(10);
  
  //
  // Update for tau2
  //
  
  term.resize(0);
  term.push_back(&potTau2);
  term.push_back(&potNu);
  term.push_back(&potDDelta);
  potSum = new PotentialSum(term);
  update.push_back(new UpdateTau2MH(&str,potSum,0.02));
  nUpdate.push_back(20);
  delete potSum;
  potSum = NULL;
  
  
  //
  // Set non-default number of updates of each type
  //
  
  string filenameUpdate(argv[6]);
  if (argv[6][0] != '=')
    str.setNumberOfUpdates(filenameUpdate,nUpdate,update);
  cout << "Number of updates in one iteration (epsilon):\n\n";
  cout << "nu:     " << nUpdate[0] << "\n";
  cout << "Delta:  " << nUpdate[1] << "\n";
  cout << "a:      " << nUpdate[2] << " (" << update[2]->getEpsilon() << ")\n";
  cout << "b:      " << nUpdate[3] << " (" << update[3]->getEpsilon() << ")\n";
  cout << "c2:     " << nUpdate[4] << "\n";
  cout << "gamma2: " << nUpdate[5] << "\n";
  cout << "r:      " << nUpdate[6] << " (" << update[6]->getEpsilon() << ")\n";
  cout << "rho:    " << nUpdate[7] << " (" << update[7]->getEpsilon() << ")\n";
  cout << "delta:  " << nUpdate[8] << "\n";
  cout << "xi:     " << nUpdate[9] << "\n";
  cout << "sigma2: " << nUpdate[10] << " (" << update[10]->getEpsilon() << ")\n";
  cout << "t:      " << nUpdate[11] << " (" << update[11]->getEpsilon() << ")\n";
  cout << "l:      " << nUpdate[12] << " (" << update[12]->getEpsilon() << ")\n";
  cout << "phi:    " << nUpdate[13] << " (" << update[13]->getEpsilon() << ")\n";
  cout << "theta:  " << nUpdate[14] << " (" << update[14]->getEpsilon() << ")\n";
  cout << "lambda: " << nUpdate[15] << " (" << update[15]->getEpsilon() << ")\n";
  cout << "tau2:   " << nUpdate[16] << " (" << update[16]->getEpsilon() << ")\n";
  cout << "\n\n";
  
  //
  // Define which output to generate
  //
  
  vector<Report *> report;
  int nBetweenReport = 1;
  string filenameOutput(argv[7]);
  ReportDiffexpressed *reportDiffexpressed = NULL;
  if (argv[7][0] != '=')
    reportDiffexpressed = str.setReports(filenameOutput,nBetweenReport,potAll,update,report);
  
  //
  // run Metropolis-Hastings algorithm
  //

  int nSinceReport = 0;
  for (i = 0; i < nIt; i++)
    {
      int k;
      for (k = 0; k < update.size(); k++)
	{
	  int t;
	  for (t = 0; t < nUpdate[k]; t++)
	    update[k]->update(ran);
	}

      if (reportDiffexpressed != NULL) reportDiffexpressed->update(&str);
      nSinceReport++;

      if (nSinceReport == nBetweenReport)
	{
	  nSinceReport = 0;

	  for (k = 0; k < report.size(); k++)
	    report[k]->report(&str);
	}
    }
  
  for (i = 0; i < report.size(); i++)
    delete report[i];

  for (i = 0; i < update.size(); i++)
    delete update[i];
  
  return 0;
}
