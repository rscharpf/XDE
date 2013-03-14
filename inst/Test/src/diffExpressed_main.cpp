//wrapper for etc.. R interface for the same?


#include <R.h>

#include <string>
#include <sstream>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <iostream>
#include <fstream>

using namespace std;

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
#include "UpdateTau2RDDeltaMH.h"
#include "UpdateTau2RhoNuMH.h"
#include "UpdateTau2RMH.h"
#include "UpdateTau2RhoMH.h"

#include "ReportDiffexpressed.h"

#define CHECKINPUT 0
#define CHECKGIBBS 0



extern "C" {

  void xdeLIN_main(int *Iseed,int *IshowIterations,int *InIt,int *IoneDelta,
		   int *IQ,int *IG,int *IS,
		   double *Ix,int *IPsi,int *IspecifiedInitialValues,
		   double *INu,double *IDDelta,double *IA,
		   double *IB,double *IC2,double *IGamma2,
		   double *IR,double *IRho,int *IDelta,
		   double *IXi,double *ISigma2,double *IT,
		   double *IL,double *IPhi,double *ITheta,
		   double *ILambda,double *ITau2R,double *ITau2Rho,double *Iparam,
		   int *InUpdate,double *Iepsilon,int *Ioutput,
		   int *IwriteToFile,char **Idirectory,
		   double *OvaluePotential,double *OvalueAcceptance,
		   double *OvalueNu,double *OvalueDDelta,
		   double *OvalueA,double *OvalueB,double *OvalueC2,
		   double *OvalueGamma2,double *OvalueR,
		   double *OvalueRho,int *OvalueDelta,double *OvalueXi,
		   double *OvalueSigma2,double *OvalueT,double *OvalueL,
		   double *OvaluePhi,double *OvalueTheta,
		   double *OvalueLambda,double *OvalueTau2R,double *OvalueTau2Rho,
		   double *OvalueProbDelta,double *OvalueDiffexpressed,
		   int *IwriteDiffexpressedTofile)
  {
    unsigned int seed = (unsigned int) *Iseed;
    if (CHECKINPUT)
      cout << "seed: " << seed << "\n";
    Random ran(seed);
    
    int nIt = *InIt;
    if (CHECKINPUT)
      cout << "nIt: " << nIt << "\n";
    
    int Q = *IQ;
    if (CHECKINPUT)
      cout << "Q: " << Q << "\n";
    
    int G = *IG;
    if (CHECKINPUT)
      cout << "G: " << G << "\n";

    int *S = IS;
    if (CHECKINPUT)
      {
	cout << "S: ";
	int q;
	for (q = 0; q < Q; q++)
	  cout << S[q] << " ";
	cout << "\n";
      }

    double *x = Ix;
    int *psi = IPsi;

    int oneDelta = *IoneDelta;

    //
    // Initialise data and parameters
    //

    Structure str(Q,G,S,x,psi,ran,CHECKINPUT,oneDelta);

    //
    // Set hyper-parameter values
    //

    double *param = Iparam;
    str.setParameterValues(param);
    if (CHECKINPUT) str.printParameterValues();

    //
    // Set initial parameter values
    //

    if (*IspecifiedInitialValues)
      str.setInitialValues(INu,IDDelta,IA,IB,IC2,IGamma2,IR,IRho,IDelta,
			   IXi,ISigma2,IT,IL,IPhi,ITheta,ILambda,ITau2R,ITau2Rho);

    //
    // Define potential functions for each variable
    //

#include "PotentialFunction.h"
    
    //
    // Set number of updates fo each type
    //

    str.setNumberOfUpdates(InUpdate,Iepsilon,nUpdate,update);
    if (CHECKINPUT)
      {
	cout << "Number of updates in one iteration (epsilon):\n\n";
	cout << "nu:     " << nUpdate[0] << " (" << update[0]->getEpsilon() << ")\n";
	cout << "Delta:  " << nUpdate[1] << " (" << update[1]->getEpsilon() << ")\n";
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
	cout << "tau2R:  " << nUpdate[16] << " (" << update[16]->getEpsilon() << ")\n";
	cout << "tau2Rho:" << nUpdate[17] << " (" << update[17]->getEpsilon() << ")\n";
	cout << "\n\n";
      }


    //
    // Define which output to generate
    //

    vector<Report *> report;
    int nBetweenReport = 1;

    ReportDiffexpressed *reportDiffexpressed = NULL;
    reportDiffexpressed = str.setReports(Ioutput,nBetweenReport,*IwriteToFile,Idirectory,
					 potAll,update,report,OvaluePotential,
					 OvalueAcceptance,OvalueNu,OvalueDDelta,
					 OvalueA,OvalueB,OvalueC2,OvalueGamma2,OvalueR,
					 OvalueRho,OvalueDelta,OvalueXi,OvalueSigma2,OvalueT,
					 OvalueL,OvaluePhi,OvalueTheta,OvalueLambda,OvalueTau2R,
					 OvalueTau2Rho,OvalueProbDelta,OvalueDiffexpressed,
					 IwriteDiffexpressedTofile,oneDelta);

    //
    // run Metropolis-Hastings algorithm
    //
    
    int nSinceReport = 0;
    for (i = 0; i < nIt; i++)
      {
	if (*IshowIterations) cout << i+1 << " ";
	cout.flush();
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
	    
	    int k;
	    for (k = 0; k < report.size(); k++)
	      report[k]->report(&str);
	  }
      }

    if (*IshowIterations) cout << "\n\n";
    
    for (i = 0; i < report.size(); i++)
      delete report[i];
    
    for (i = 0; i < update.size(); i++)
      delete update[i];


    str.setFinalValues(INu,IDDelta,IA,IB,IC2,IGamma2,IR,IRho,IDelta,
		       IXi,ISigma2,IT,IL,IPhi,ITheta,ILambda,ITau2R,ITau2Rho);
    

    seed = ran.ChangeSeed(seed);
    *Iseed = (int) seed;
    
    return;
  }
  
} // extern "C"
