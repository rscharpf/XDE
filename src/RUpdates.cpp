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

  void RUpdate(int *updateType,int *Iseed,int *InIt,int *IoneDelta,
	       int *IQ,int *IG,int *IS,
	       double *Ix,int *IPsi,
	       double *INu,double *IDDelta,double *IA,
	       double *IB,double *IC2,double *IGamma2,
	       double *IR,double *IRho,int *IDelta,
	       double *IXi,double *ISigma2,double *IT,
	       double *IL,double *IPhi,double *ITheta,
	       double *ILambda,double *ITau2R,double *ITau2Rho,double *Iparam,
	       double *Iepsilon)
  {
    unsigned int seed = (unsigned int) *Iseed;
    if (CHECKINPUT)
      cout << "seed: " << seed << "\n";
    Random ran(seed);
    
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
    // Initialise necessary variables and do the update
    //

    switch (*updateType) {
    case 0: // update for nu and tau2Rho
      Structure str(Q,G,S,x,psi,ran,CHECKINPUT,oneDelta);

      double *param = Iparam;
      str.setParameterValues(param);
      if (CHECKINPUT) str.printParameterValues();

      str.setInitialValues(INu,IDDelta,IA,IB,IC2,IGamma2,IR,IRho,IDelta,
			   IXi,ISigma2,IT,IL,IPhi,ITheta,ILambda,
			   ITau2R,ITau2Rho);
      
      PotentialZero potTau2R;
      
      UpdateTau2RhoNuMH update(&str,&potTau2Rho,0.02);
      update.setEpsilon(*Iepsilon);
      update.update(ran);

      str.setFinalValues(INu,IDDelta,IA,IB,IC2,IGamma2,IR,IRho,IDelta,
			 IXi,ISigma2,IT,IL,IPhi,ITheta,ILambda,ITau2R,ITau2Rho);
      break;
      
    case 1: // update for Delta and tau2R
      
      break;
    case 2: // and so on for each update type
    }

    
    seed = ran.ChangeSeed(seed);
    *Iseed = (int) seed;

    return;
  }
  
} // extern "C"
