#ifndef STRUCTURE_H
#define STRUCTURE_H

#include <string>
#include <sstream>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <iostream>
#include <fstream>

#include "Random.h"
class Potential;
class Update;
class Report;
class ReportDiffexpressed;

using namespace std;


class Structure
{
  friend class PotentialZero;
  friend class PotentialX;
  friend class PotentialNu;
  friend class PotentialDDelta;
  friend class PotentialA;
  friend class PotentialB;
  friend class PotentialR;
  friend class PotentialRho;
  friend class PotentialDelta;
  friend class PotentialDeltag;
  friend class PotentialXi;
  friend class PotentialSigma2;
  friend class PotentialPhi;

  friend class PotentialSigma2qg;
  friend class PotentialXqg;
  friend class PotentialNug;
  friend class PotentialDDeltag;
  friend class PotentialPhiqg;

  friend class UpdateAMH;
  friend class UpdateBMH;
  friend class UpdateC2Gibbs;
  friend class UpdateGamma2Gibbs;
  friend class UpdateRC2MH;
  friend class UpdateRhoGamma2MH;
  friend class UpdateDeltaMH;
  friend class UpdateXiGibbs;
  friend class UpdateSigma2MH;
  friend class UpdateTMH;
  friend class UpdateLMH;
  friend class UpdatePhiMH;
  friend class UpdateThetaMH;
  friend class UpdateLambdaMH;
  friend class UpdateTau2RDDeltaMH;
  friend class UpdateTau2RhoNuMH;
  friend class UpdateTau2RMH;
  friend class UpdateTau2RhoMH;

  friend class ReportNu;
  friend class ReportDDelta;
  friend class ReportA;
  friend class ReportB;
  friend class ReportC2;
  friend class ReportGamma2;
  friend class ReportR;
  friend class ReportRho;
  friend class ReportDelta;
  friend class ReportXi;
  friend class ReportSigma2;
  friend class ReportT;
  friend class ReportL;
  friend class ReportPhi;
  friend class ReportTheta;
  friend class ReportLambda;
  friend class ReportTau2R;
  friend class ReportTau2Rho;
  friend class ReportProbDelta;
  friend class ReportDiffexpressed;

 public:
  Structure(string &filename,Random &ran,int oneDelta);
  Structure(int P,int G,int *S,double *x,int *psi,Random &ran,int checkinput,
	    int oneDelta);
  virtual ~Structure(void);

  void setParameterValues(string &filename);
  void setParameterValues(double *param);

  void printParameterValues(void) const;

  void setInitialValues(string &filename);
  void setInitialValues(double *Nu,double *DDelta,double *A,
			double *B,double *C2,double *Gamma2,
			double *R,double *Rho,int *Delta,
			double *Xi,double *Sigma2,double *T,
			double *L,double *Phi,double *Theta,
			double *Lambda,double *Tau2R,double *Tau2Rho);
  void setFinalValues(double *Nu,double *DDelta,double *A,
		      double *B,double *C2,double *Gamma2,
		      double *R,double *Rho,int *Delta,
		      double *Xi,double *Sigma2,double *T,
		      double *L,double *Phi,double *Theta,
		      double *Lambda,double *Tau2R,double *Tau2Rho) const;


  void setNumberOfUpdates(string &filename,vector<int> &nUpdate,vector<Update *> &update) const;
  void setNumberOfUpdates(int *InUpdate,double *Iepsilon,
			  vector<int> &nUpdate,vector<Update *> &update);

  ReportDiffexpressed *setReports(string &filename,int &nBetweenReport,vector<Potential *> &potential,
		  vector<Update *> &update,vector<Report *> &report,int oneDelta);
  ReportDiffexpressed *setReports(int *output,int &nBetweenReport,int writeToFile,
		  char **directory,vector<Potential *> &potential,
		  vector<Update *> &update,vector<Report *> &report,
		  double *valuePotential,double *valueAcceptance,
		  double *valueNu,double *valueDDelta,double *valueA,
		  double *valueB,double *valueC2,double *valueGamma2,
		  double *valueR,double *valueRho,int *valueDelta,
		  double *valueXi,double *valueSigma2,double *valueT,
		  double *valueL,double *valuePhi,double *valueTheta,
		  double *valueLambda,double *valueTau2R,double *valueTau2Rho,
		  double *valueProbDelta,double *valueDiffexpressed,
		  int *writeDiffexpressedTofile,
		  int oneDelta);


 protected: 
  void readCorrelationValues(vector<vector<double> > &corr,const char *line,int bufSize,
			     const string varName,string &filename) const;
  void allocateSpace(void);
  void initialiseVariables(Random &ran,int oneDelta);

  //
  // fixed parameter and variables
  //

  int G;                // Number of genes
  int Q;                // Number of studies
  vector<int> S;        // S[q], Number of samples in study q
  
  vector<vector<vector<double> > > x;   // x[q][g][s], Expression value
  vector<vector<int> > psi;             // psi[q][s], clinical variable
  double alphaA,betaA,pA0,pA1;          // parameters in beta prior for a
  double alphaB,betaB,pB0,pB1;          // parameters in beta prior for b
  double nuR;                           // parameter in prior for r
  double nuRho;                         // parameter in prior for rho
  double alphaXi,betaXi;                // parameters in beta prior for xi
  double c2Max;                         // parameters in uniform prior for c2

  //
  // parameters to be simulated
  //

  vector<vector<double> > nu;           // nu[q][g]
  vector<vector<double> > Delta;        // Delta[q][g]
  vector<vector<int> > delta;           // delta[q][g]
  vector<double> a;                     // a[q]
  vector<double> b;                     // b[q]
  double c2; 
  double gamma2;
  vector<double> tau2R;                 // tau2R[q]
  vector<double> tau2Rho;               // tau2Rho[q]
  vector<vector<double> > r;            // r[q1][q2]
  vector<vector<double> > rho;          // rho[q1][q2]
  vector<double> xi;                    // xi[q]
  vector<vector<double> > sigma2;       // sigma2[q][g]
  vector<double> t;                     // t[q];
  vector<double> l;                     // l[q]
  vector<vector<double> > phi;          // phi[q][g]
  vector<double> theta;                 // theta[q]
  vector<double> lambda;                // lambda[q]
};






#endif
