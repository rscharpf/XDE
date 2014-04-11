//#include <strstream>
#include <cstdio>
#include <cstring>

#include "Structure.h"
#include "Cholesky.h"

#include "Potential.h"
#include "Update.h"

#include "ReportNu.h"
#include "ReportDDelta.h"
#include "ReportA.h"
#include "ReportB.h"
#include "ReportC2.h"
#include "ReportGamma2.h"
#include "ReportR.h"
#include "ReportRho.h"
#include "ReportDelta.h"
#include "ReportXi.h"
#include "ReportSigma2.h"
#include "ReportT.h"
#include "ReportL.h"
#include "ReportPhi.h"
#include "ReportTheta.h"
#include "ReportLambda.h"
#include "ReportTau2R.h"
#include "ReportTau2Rho.h"
#include "ReportPotential.h"
#include "ReportAcceptance.h"
#include "ReportProbDelta.h"
#include "ReportDiffexpressed.h"


Structure::Structure(string &filename,Random &ran,int oneDelta)
{
  ifstream in;
  in.open(filename.c_str());
//  if (in.fail())
//    {
//      //cout << "ERROR: Unable to open file " << filename.c_str() << ". Aborting.\n\n";
//      //exit(-1);
//    }

  Q = 0;
  in >> Q;
//  if (Q <= 0)
//    {
//      cout << "ERROR: Illegal number of studies, Q = " << Q << ", read in file " << filename << ". Aborting.\n\n";
//      exit(-1);
//    }
  S.resize(Q);

  vector<string> fileData(Q);
  vector<string> filePsi(Q);
  int q;
  for (q = 0; q < Q; q++)
    {
      in >> fileData[q];
//      if (in.fail())
//	{
//	  cout << "ERROR: Unable to read expression value file number " << q + 1 << " from file " << filename << ". Aborting.\n\n";
//	  exit(-1);
//	}

      in >> filePsi[q];
//      if (in.fail())
//	{
//	  cout << "ERROR: Unable to read clinical value file number " << q + 1 << " from file " << filename << ". Aborting.\n\n";
//	  exit(-1);
//	}
    }
  in.close();

  for (q = 0; q < Q; q++)
    {
      ifstream inData;
      inData.open(fileData[q].c_str());
//      if (inData.fail())
//	{
//	  cout << "ERROR: Error when reading " << filename.c_str() << ". Unable to open file " << fileData[q] << ". Aborting.\n\n";
//	  exit(-1);
//	}

      ifstream inPsi;
      inPsi.open(filePsi[q].c_str());
//      if (inPsi.fail())
//	{
//	  cout << "ERROR: Error when reading " << filename.c_str() << ". Unable to open file " << filePsi[q] << ". Aborting.\n\n";
//	  exit(-1);
//	}

      int GG = 0;
      inData >> GG;
//      if (GG <= 0)
//	{
//	  cout << "ERROR: Illegal number of genes, G = " << G << ", read in file " << fileData[q] << ". Aborting.\n\n";
//	  exit(-1);
//	}
      if (q == 0)
	G = GG;
      //      else
      //	{
//	  if (GG != G)
//	    {
//	      cout << "ERROR: Inconsistent number of genes in file " << fileData[0] << " and file " << fileData[q] << ". Aborting.\n\n";
//	      exit(-1);
//	    }
//	}


      S[q] = 0;
      inData >> S[q];
//      if (S[q] <= 0)
//	{
//	  cout << "ERROR: Illegal number of samples, S = " << S[q] << ", read from " << fileData[q] << ". Aborting.\n\n";
//	  exit(-1);
//	}
      int SS = 0;
      inPsi >> SS;
//      if (SS <= 0)
//	{
//	  cout << "ERROR: Illegal number of samples, S = " << SS << ", read from " << filePsi[q] << ". Aborting.\n\n";
//	  exit(-1);
//	}
//      if (SS != S[q])
//	{
//	  cout << "ERROR: Inconsistent number of samples read in " << fileData[q] << " and " << filePsi[q] << ". Aborting.\n\n";
//	  exit(-1);
//	}
      inData.close();
      inPsi.close();
    }

  //
  // Print data size
  //

//  cout << "Size of data set:\n\n";
//  cout << "P: " << Q << "\n";
//  cout << "G: " << G << "\n";
//  cout << "S: ";
//  for (q = 0; q < Q; q++)
//    cout << S[q] << " ";
//  cout << "\n\n\n";

  //
  // allocate neccesary space
  //

  allocateSpace();

  //
  // read expression values from files
  //

  for (q = 0; q < Q; q++)
    {
      ifstream in;
      in.open(fileData[q].c_str());
      int GG,SS;
      in >> GG;
      in >> SS;
      int g,s;
      for (g = 0; g < G; g++)
	for (s = 0; s < S[q]; s++)
	  in >> x[q][g][s];

      in.close();
    }

  //
  // read clinical variables from files
  //

  for (q = 0; q < Q; q++)
    {
      ifstream in;
      in.open(filePsi[q].c_str());
      int SS;
      in >> SS;
      int s;
      for (s = 0; s < S[q]; s++)
	{
	  in >> psi[q][s];
//	  if (psi[q][s] != 0 && psi[q][s] != 1)
//	    {
//	      cout << "ERROR: value different from 0 or 1 found in \"" << filePsi[q] << "\". Aborting.\n";
//	      exit(-1);
//	    }
	}

      in.close();
    }

  //
  // initialise remaining variables
  //

  initialiseVariables(ran,oneDelta);

  return;
}



Structure::Structure(int Q,int G,int *S,double *x,int *psi,Random &ran,int checkinput,int oneDelta)
{
  this->Q = Q;
  this->G = G;
  this->S.resize(this->Q);
  int q;
  for (q = 0; q < this->Q; q++)
    this->S[q] = S[q];

  //
  // allocate necessary space
  //

  allocateSpace();

  //
  // set expression values
  //

  int g,s;
  int nr = 0;
  for (q = 0; q < this->Q; q++)
    for (g = 0; g < this->G; g++)
      for (s = 0; s < this->S[q]; s++)
	{
	  this->x[q][g][s] = x[nr];
	  nr++;
	}

  //
  // set clinical variables
  //

  nr = 0;
  for (q = 0; q < this->Q; q++)
    for (s = 0; s < this->S[q]; s++)
      {
	this->psi[q][s] = psi[nr];
	nr++;
      }

  //  if (checkinput)
  //    {
      //cout << "Expression values:\n";
//      for (q = 0; q < this->Q; q++)
//	{
//	  cout << "first values of study " << q << "\n";
//	  cout << this->x[q][0][0] << " " << this->x[q][0][1] << "\n";
//	  cout << this->x[q][1][0] << " " << this->x[q][1][1] << "\n";
//	}
//      cout << "\n";

//      cout << "Clinical values:\n";
//      for (q = 0; q < this->Q; q++)
//	{
//	  cout << "study " << q << ": ";
//	  for (s = 0; s < this->S[q]; s++)
//	    cout << this->psi[q][s] << " ";
//	  cout << "\n";
//	}
//      cout << "\n";
//    }

  //
  // initialise remaining parameters
  //

  initialiseVariables(ran,oneDelta);

  return;
}








Structure::~Structure(void)
{
  return;
}



void Structure::allocateSpace(void)
{
  x.resize(Q);
  int q;
  for (q = 0; q < Q; q++)
    {
      x[q].resize(G);
      int g;
      for (g = 0; g < G; g++)
	x[q][g].resize(S[q]);
    }

  psi.resize(Q);
  for (q = 0; q < Q; q++)
    psi[q].resize(S[q]);

  nu.resize(Q);
  for (q = 0; q < Q; q++)
    nu[q].resize(G);

  Delta.resize(Q);
  for (q = 0; q < Q; q++)
    Delta[q].resize(G);

  delta.resize(Q);
  for (q = 0; q < Q; q++)
    delta[q].resize(G);
  xi.resize(Q);
  a.resize(Q);
  b.resize(Q);
  tau2R.resize(Q);
  tau2Rho.resize(Q);

  sigma2.resize(Q);
  for (q = 0; q < Q; q++)
    sigma2[q].resize(G);

  t.resize(Q);
  l.resize(Q);

  phi.resize(Q);
  for (q = 0; q < Q; q++)
    phi[q].resize(G);

  theta.resize(Q);
  lambda.resize(Q);


  r.resize(Q);
  rho.resize(Q);
  for (q = 0; q < Q; q++)
    {
      r[q].resize(Q);
      rho[q].resize(Q);
    }

  return;
}


void Structure::initialiseVariables(Random &ran,int oneDelta)
{
  //
  // initialise remaining fixed variables
  //

  alphaA = 1.0;
  betaA = 1.0;
  pA0 = 0.0;
  pA1 = 0.0;

  alphaB = 1.0;
  betaB = 1.0;
  pB0 = 0.0;
  pB1 = 0.0;

  nuR = Q + 1.0;
  nuRho = Q + 1.0;

  alphaXi = 1.0;
  betaXi = 1.0;

  c2Max = 1.0;

  //
  // initialise variable parameters
  //

  int isneg = 0;
  do
    {
      rho = ran.CorrelationStandardInverseWishart(Q,nuRho);
      r = ran.CorrelationStandardInverseWishart(Q,nuR);

      isneg = 0;
      int p,q;
      for (p = 0; p < Q; p++)
	for (q = 0; q < Q; q++)
	  {
	    isneg += (r[p][q] < 0.0);
	    isneg += (rho[p][q] < 0.0);
	  }
    }
  while (isneg > 0);


  gamma2 = 1.0;
  c2 = c2Max / 2.0;
  int q;
  for (q = 0; q < Q; q++)
    {
      tau2R[q] = 1.0;
      tau2Rho[q] = 1.0;
    }

  for (q = 0; q < Q; q++)
    {
      t[q] = 1.0;
      l[q] = 1.0;
    }

  for (q = 0; q < Q; q++)
    {
      int g;
      for (g = 0; g < G; g++)
	sigma2[q][g] = ran.InverseGamma(t[q],l[q]);
    }

  for (q = 0; q < Q; q++)
    {
      a[q] = ran.Beta(alphaA,betaA);
      b[q] = ran.Beta(alphaB,betaB);
    }

  int g;
  for (g = 0; g < G; g++)
    {
      vector<vector<double> > Sigma;
      Sigma.resize(Q);
      for (q = 0; q < Q; q++)
	Sigma[q].resize(Q);
      int p;
      for (p = 0; p < Q; p++)
	for (q = 0; q < Q; q++)
	  {
	    Sigma[p][q] = gamma2;
	    if (p != q) Sigma[p][q] *= rho[p][q];
	    Sigma[p][q] *= sqrt(tau2Rho[p] * tau2Rho[q]);
	    Sigma[p][q] *= exp(0.5 * (a[q] * log(sigma2[q][g]) + a[p] * log(sigma2[p][g])));
	  }
      vector<double> zero(Q,0);
      vector<double> value(Q,0);
      value = ran.MultiGaussian(Sigma,zero);
      for (q = 0; q < Q; q++)
	nu[q][g] = value[q];
    }


  for (g = 0; g < G; g++)
    {
      vector<vector<double> > R;
      R.resize(Q);
      for (q = 0; q < Q; q++)
	R[q].resize(Q);
      int p;
      for (p = 0; p < Q; p++)
	for (q = 0; q < Q; q++)
	  {
	    R[p][q] = c2;
	    if (p != q) R[p][q] *= r[p][q];
	    R[p][q] *= sqrt(tau2R[p] * tau2R[q]);
	    R[p][q] *= exp(0.5 * (b[q] * log(sigma2[q][g]) + b[p] * log(sigma2[p][g])));
	  }
      vector<double> zero(Q,0);
      vector<double> value(Q,0);
      value = ran.MultiGaussian(R,zero);
      for (q = 0; q < Q; q++)
	Delta[q][g] = value[q];
    }

  if (oneDelta == 0)
    {
      for (q = 0; q < Q; q++)
	xi[q] = ran.Beta(alphaXi,betaXi);
      for (q = 0; q < Q; q++)
	for (g = 0; g < G; g++)
	  delta[q][g] = (ran.Unif01() <= xi[q]);
    }
  else
    {
      xi[0] = ran.Beta(alphaXi,betaXi);
      for (q = 1; q < Q; q++)
	xi[q] = xi[0];
      for (g = 0; g < G; g++)
	{
	  int dd = (ran.Unif01() <= xi[0]);
	  for (q = 0; q < Q; q++)
	    delta[q][g] = dd;
	}
    }


  for (q = 0; q < Q; q++)
    {
      theta[q] = 1.0;
      lambda[q] = 1.0;
    }

  for (q = 0; q < Q; q++)
    for (g = 0; g < G; g++)
      phi[q][g] = ran.Gamma(theta[q],lambda[q]);

  return;
}




void Structure::setParameterValues(string &filename)
{
  ifstream in;
  in.open(filename.c_str());
//  if (in.fail())
//    {
//      cout << "ERROR: Unable to open file " << filename.c_str() << ". Aborting.\n\n";
//      exit(-1);
//    }

  //
  // read first line in file
  //

  int bufSize = 1000;
  char buf[bufSize];
  in.get(buf,bufSize,'\n');
  char c;
//  if (in.get(c) && c != '\n')
//    {
//      couta << "ERROR: Line 1 in file " << filename << " is too long. Maximum line length is ";
//      cout << bufSize << ". Aborting.\n\n";
//      exit(-1);
//    }

  char var1[bufSize];
  char var2[bufSize];
  char var3[bufSize];
  char var4[bufSize];
  int nRead = sscanf(buf,"%s %s %s %s",var1,var2,var3,var4);
  if (nRead != 4)
    {
      if (in.eof()) return;
//      cout << "ERROR: Error when reading line 1 in file " << filename << ". Aborting.\n\n";
//      exit(-1);
    }
  if (var1[0] != '=')
    {
      nRead = sscanf(var1,"%le",&alphaA);
//      if (nRead != 1)
//	{
//	  cout << "ERROR: Error when reading alphaA from file " << filename << ". Aborting.\n\n";
//	  exit(-1);
//	}
//      else if (alphaA <= 0.0)
//	{
//	  cout << "ERROR: Illegal value given for alphaA in file " << filename << ". Aborting.\n\n";
//	  exit(-1);
//	}
    }

  if (var2[0] != '=')
    {
      nRead = sscanf(var2,"%le",&betaA);
//      if (nRead != 1)
//	{
//	  cout << "ERROR: Error when reading betaA from file " << filename << ". Aborting.\n\n";
//	  exit(-1);
//	}
//      else if (betaA <= 0.0)
//	{
//	  cout << "ERROR: Illegal value given for betaA in file " << filename << ". Aborting.\n\n";
//	  exit(-1);
//	}
    }

  if (var3[0] != '=')
    {
      nRead = sscanf(var3,"%le",&pA0);
//      if (nRead != 1)
//	{
//	  cout << "ERROR: Error when reading pA0 from file " << filename << ". Aborting.\n\n";
//	  exit(-1);
//	}
//      else if (pA0 < 0.0 && pA0 > 1.0)
//	{
//	  cout << "ERROR: Illegal value given for pA0 in file " << filename << ". Aborting.\n\n";
//	  exit(-1);
//	}
    }

  if (var4[0] != '=')
    {
      nRead = sscanf(var4,"%le",&pA1);
//      if (nRead != 1)
//	{
//	  cout << "ERROR: Error when reading pA1 from file " << filename << ". Aborting.\n\n";
//	  exit(-1);
//	}
//      else if (pA1 < 0.0 && pA1 > 1.0)
//	{
//	  cout << "ERROR: Illegal value given for pA1 in file " << filename << ". Aborting.\n\n";
//	  exit(-1);
//	}
    }

//  if (pA0 + pA1 >= 1.0)
//    {
//      cout << "ERROR: Illegal values given for pA0 and pA1 in file " << filename << ". Aborting.\n\n";
//      exit(-1);
//    }


  //
  // read second line in file
  //


  in.get(buf,bufSize,'\n');
//  if (in.get(c) && c != '\n')
//    {
//      cout << "ERROR: Line 2 in file " << filename << " is too long. Maximum line length is ";
//      cout << bufSize << ". Aborting.\n\n";
//      exit(-1);
//    }

  nRead = sscanf(buf,"%s %s %s %s",var1,var2,var3,var4);
  if (nRead != 4)
    {
      if (in.eof()) return;
//      cout << "ERROR: Error when reading line 2 in file " << filename << ". Aborting.\n\n";
//      exit(-1);
    }
  if (var1[0] != '=')
    {
      nRead = sscanf(var1,"%le",&alphaB);
//      if (nRead != 1)
//	{
//	  cout << "ERROR: Error when reading alphaB from file " << filename << ". Aborting.\n\n";
//	  exit(-1);
//	}
//      else if (alphaB <= 0.0)
//	{
//	  cout << "ERROR: Illegal value given for alphaB in file " << filename << ". Aborting.\n\n";
//	  exit(-1);
//	}
    }

  if (var2[0] != '=')
    {
      nRead = sscanf(var2,"%le",&betaB);
//      if (nRead != 1)
//	{
//	  cout << "ERROR: Error when reading betaB from file " << filename << ". Aborting.\n\n";
//	  exit(-1);
//	}
//      else if (betaB <= 0.0)
//	{
//	  cout << "ERROR: Illegal value given for betaB in file " << filename << ". Aborting.\n\n";
//	  exit(-1);
//	}
    }

  if (var3[0] != '=')
    {
      nRead = sscanf(var3,"%le",&pB0);
//      if (nRead != 1)
//	{
//	  cout << "ERROR: Error when reading pB0 from file " << filename << ". Aborting.\n\n";
//	  exit(-1);
//	}
//      else if (pB0 < 0.0 && pB0 > 1.0)
//	{
//	  cout << "ERROR: Illegal value given for pB0 in file " << filename << ". Aborting.\n\n";
//	  exit(-1);
//	}
    }

  if (var4[0] != '=')
    {
      nRead = sscanf(var4,"%le",&pB1);
//      if (nRead != 1)
//	{
//	  cout << "ERROR: Error when reading pB1 from file " << filename << ". Aborting.\n\n";
//	  exit(-1);
//	}
//      else if (pB1 < 0.0 && pB1 > 1.0)
//	{
//	  cout << "ERROR: Illegal value given for pB1 in file " << filename << ". Aborting.\n\n";
//	  exit(-1);
//	}
    }

//  if (pB0 + pB1 >= 1.0)
//    {
//      cout << "ERROR: Illegal values given for pB0 and pB1 in file " << filename << ". Aborting.\n\n";
//      exit(-1);
//    }

  //
  // read third line in file
  //

  in.get(buf,bufSize,'\n');
//  if (in.get(c) && c != '\n')
//    {
//      cout << "ERROR: Line 3 in file " << filename << " is too long. Maximum line length is ";
//      cout << bufSize << ". Aborting.\n\n";
//      exit(-1);
//    }

  nRead = sscanf(buf,"%s",var1);
  if (nRead != 1)
    {
      if (in.eof()) return;
//      cout << "ERROR: Error when reading line 3 in file " << filename << ". Aborting.\n\n";
//      exit(-1);
    }
  if (var1[0] != '=')
    {
      nRead = sscanf(var1,"%le",&nuR);
//      if (nRead != 1)
//	{
//	  cout << "ERROR: Error when reading nuR from file " << filename << ". Aborting.\n\n";
//	  exit(-1);
//	}
//      else if (nuR <= 0.0)
//	{
//	  cout << "ERROR: Illegal value given for nuR in file " << filename << ". Aborting.\n\n";
//	  exit(-1);
//	}
    }



  //
  // read fourth line in file
  //

  in.get(buf,bufSize,'\n');
//  if (in.get(c) && c != '\n')
//    {
//      cout << "ERROR: Line 4 in file " << filename << " is too long. Maximum line length is ";
//      cout << bufSize << ". Aborting.\n\n";
//      exit(-1);
//    }

  nRead = sscanf(buf,"%s",var1);
  if (nRead != 1)
    {
      if (in.eof()) return;
//      cout << "ERROR: Error when reading line 4 in file " << filename << ". Aborting.\n\n";
//      exit(-1);
    }
  if (var1[0] != '=')
    {
      nRead = sscanf(var1,"%le",&nuRho);
//      if (nRead != 1)
//	{
//	  cout << "ERROR: Error when reading nuRho from file " << filename << ". Aborting.\n\n";
//	  exit(-1);
//	}
//      else if (nuR <= 0.0)
//	{
//	  cout << "ERROR: Illegal value given for nuRho in file " << filename << ". Aborting.\n\n";
//	  exit(-1);
//	}
    }


  //
  // read fifth line in file
  //


  in.get(buf,bufSize,'\n');
//  if (in.get(c) && c != '\n')
//    {
//      cout << "ERROR: Line 5 in file " << filename << " is too long. Maximum line length is ";
//      cout << bufSize << ". Aborting.\n\n";
//      exit(-1);
//    }

  nRead = sscanf(buf,"%s %s",var1,var2);
  if (nRead != 2)
    {
      if (in.eof()) return;
//      cout << "ERROR: Error when reading line 5 in file " << filename << ". Aborting.\n\n";
//      exit(-1);
    }
  if (var1[0] != '=')
    {
      nRead = sscanf(var1,"%le",&alphaXi);
//      if (nRead != 1)
//	{
//	  cout << "ERROR: Error when reading alphaXi from file " << filename << ". Aborting.\n\n";
//	  exit(-1);
//	}
//      else if (alphaXi <= 0.0)
//	{
//	  cout << "ERROR: Illegal value given for alphaXi in file " << filename << ". Aborting.\n\n";
//	  exit(-1);
//	}
    }

  if (var2[0] != '=')
    {
      nRead = sscanf(var2,"%le",&betaXi);
//      if (nRead != 1)
//	{
//	  cout << "ERROR: Error when reading betaXi from file " << filename << ". Aborting.\n\n";
//	  exit(-1);
//	}
//      else if (betaXi <= 0.0)
//	{
//	  cout << "ERROR: Illegal value given for betaXi in file " << filename << ". Aborting.\n\n";
//	  exit(-1);
//	}
    }

  //
  // read sixth line in file
  //

  in.get(buf,bufSize,'\n');
//  if (in.get(c) && c != '\n')
//    {
//      cout << "ERROR: Line 6 in file " << filename << " is too long. Maximum line length is ";
//      cout << bufSize << ". Aborting.\n\n";
//      exit(-1);
//    }

  nRead = sscanf(buf,"%s",var1);
  if (nRead != 1)
    {
      if (in.eof()) return;
//      cout << "ERROR: Error when reading line 6 in file " << filename << ". Aborting.\n\n";
//      exit(-1);
    }
  if (var1[0] != '=')
    {
      nRead = sscanf(var1,"%le",&c2Max);
//      if (nRead != 1)
//	{
//	  cout << "ERROR: Error when reading c2Max from file " << filename << ". Aborting.\n\n";
//	  exit(-1);
//	}
//      else if (c2Max <= 0.0)
//	{
//	  cout << "ERROR: Illegal value given for c2Max in file " << filename << ". Aborting.\n\n";
//	  exit(-1);
//	}
    }
  c2 = c2Max / 2.0;

  return;
}




void Structure::setParameterValues(double *param)
{
  int nr = 0;

  alphaA = param[nr];
  nr++;

  betaA = param[nr];
  nr++;

  pA0 = param[nr];
  nr++;

  pA1 = param[nr];
  nr++;

  alphaB = param[nr];
  nr++;

  betaB = param[nr];
  nr++;

  pB0 = param[nr];
  nr++;

  pB1 = param[nr];
  nr++;

  nuR = param[nr];
  nr++;

  nuRho = param[nr];
  nr++;

  alphaXi = param[nr];
  nr++;

  betaXi = param[nr];
  nr++;

  c2Max = param[nr];
  c2 = c2Max / 2.0;
  nr++;

  return;
}




//void Structure::printParameterValues(void) const
//{
//  cout << "Hyper-parameter values:\n\n";
//  cout << "alphaA =  " << alphaA << "\n";
//  cout << "betaA =   " << betaA << "\n";
//  cout << "pA0 =     " << pA0 << "\n";
//  cout << "pA1 =     " << pA1 << "\n";
//  cout << "alphaB =  " << alphaB << "\n";
//  cout << "betaB =   " << betaB << "\n";
//  cout << "pB0 =     " << pB0 << "\n";
//  cout << "pB1 =     " << pB1 << "\n";
//  cout << "nuR =     " << nuR << "\n";
//  cout << "nuRho =   " << nuRho << "\n";
//  cout << "alphaXi = " << alphaXi << "\n";
//  cout << "betaXi =  " << betaXi << "\n";
//  cout << "c2Max =   " << c2Max << "\n";
//  cout << "\n\n";
//
//  return;
//}




//void Structure::setInitialValues(string &filename)
//{
//  int tau2Specified = 0;
//
//  ifstream in;
//  in.open(filename.c_str());
//  if (in.fail())
//    {
//      cout << "ERROR: Unable to open file " << filename.c_str() << ". Aborting.\n\n";
//      exit(-1);
//    }
//
//  //
//  // read one and one line
//  //
//
//  int line;
//  for (line = 1; line <= 16 + Q; line++)
//    {
//      int bufSize = 1000;
//      char buf[bufSize];
//      in.get(buf,bufSize,'\n');
//      char c;
//      if (in.get(c) && c != '\n')
//	{
//	  cout << "ERROR: Line " << line << " in file " << filename << " is too long. Maximum line length is ";
//	  cout << bufSize << ". Aborting.\n\n";
//	  exit(-1);
//	}
//
//      char var[bufSize];
//      int nRead = sscanf(buf,"%s",var);
//      if (nRead != 1)
//	{
//	  if (in.eof())
//	    {
//	      if (line > 16 && tau2Specified)
//		{
//		  double prod = 1.0;
//		  int q;
//		  for (q = 0; q < Q; q++)
//		    prod *= tau2R[q];
//		  for (q = 0; q < Q; q++)
//		    tau2R[q] /= exp(log(prod) / Q);
//
//		  for (q = 0; q < Q; q++)
//		    cout << "tau2R[" << q+1 << "] = " << tau2R[q] << "\n";
//		}
//	      return;
//	    }
//	  cout << "ERROR: Error when reading line " << line << " in file " << filename << ". Aborting.\n\n";
//	  exit(-1);
//	}
//      string varName;
//      if (line == 1)
//	varName = "nu";
//      else if (line == 2)
//	varName = "Delta";
//      else if (line == 3)
//	varName = "a";
//      else if (line == 4)
//	varName = "b";
//      else if (line == 5)
//	varName = "c2";
//      else if (line == 6)
//	varName = "gamma2";
//      else if (line == 7)
//	varName = "r";
//      else if (line == 8)
//	varName = "rho";
//      else if (line == 9)
//	varName = "delta";
//      else if (line == 10)
//	varName = "xi";
//      else if (line == 11)
//	varName = "sigma2";
//      else if (line == 12)
//	varName = "t";
//      else if (line == 13)
//	varName = "l";
//      else if (line == 14)
//	varName = "phi";
//      else if (line == 15)
//	varName = "theta";
//      else if (line == 16)
//	varName = "lambda";
//      else
//	{
//	  char temp[120];
//	  sprintf(temp,"tau2R[%d]",line - 16);
//	  varName = temp;
//	}
//
//
//      if (var[0] != '=')
//	{
//	  double value = 0.0;
//	  nRead = sscanf(var,"%le",&value);
//	  if (nRead != 1)
//	    {
//	      cout << "ERROR: Error when reading " << varName << " from file " << filename << ". Aborting.\n\n";
//	      exit(-1);
//	    }
//	  if (line == 1)
//	    {
//	      int q,g;
//	      for (q = 0; q < Q; q++)
//		for (g = 0; g < G; g++)
//		  nu[q][g] = value;
//	    }
//	  else if (line == 2)
//	    {
//	      int q,g;
//	      for (q = 0; q < Q; q++)
//		for (g = 0; g < G; g++)
//		  Delta[q][g] = value;
//	    }
//	  else if (line == 3)
//	    {
//	      if (value < 0.0 || value > 1.0)
//		{
//		  cout << "ERROR: Illegal value given for " << varName << " in line " << line;
//		  cout << " in file " << filename << ". Aborting.";
//		  exit(-1);
//		}
//	      int q;
//	      for (q = 0; q < Q; q++) a[q] = value;
//	    }
//	  else if (line == 4)
//	    {
//	      if (value < 0.0 || value > 1.0)
//		{
//		  cout << "ERROR: Illegal value given for " << varName << " in line " << line;
//		  cout << " in file " << filename << ". Aborting.";
//		  exit(-1);
//		}
//	      int q;
//	      for (q = 0; q < Q; q++) b[q] = value;
//	    }
//	  else if (line == 5)
//	    {
//	      if (value < 0.0)
//		{
//		  cout << "ERROR: Illegal value given for " << varName << " in line " << line;
//		  cout << " in file " << filename << ". Aborting.";
//		  exit(-1);
//		}
//
//	      c2 = value;
//	    }
//	  else if (line == 6)
//	    {
//	      if (value < 0.0)
//		{
//		  cout << "ERROR: Illegal value given for " << varName << " in line " << line;
//		  cout << " in file " << filename << ". Aborting.";
//		  exit(-1);
//		}
//	      gamma2 = value;
//	    }
//	  else if (line == 7)
//	    {
//	      readCorrelationValues(r,buf,bufSize,varName,filename);
//
//	      cout << varName << " = ";
//	      int i,j;
//	      for (i = 0; i < r.size(); i++)
//		{
//		  if (i != 0) cout << "    ";
//		  cout << "[";
//		  for (j = 0; j < r.size(); j++)
//		    {
//		      cout << r[i][j];
//		      if (j < r.size() - 1)
//			cout << ",";
//		    }
//		  cout << "]\n";
//		}
//	    }
//	  else if (line == 8)
//	    {
//	      readCorrelationValues(rho,buf,bufSize,varName,filename);
//
//	      cout << varName << " = ";
//	      int i,j;
//	      for (i = 0; i < rho.size(); i++)
//		{
//		  if (i != 0) cout << "      ";
//		  cout << "[";
//		  for (j = 0; j < rho.size(); j++)
//		    {
//		      cout << rho[i][j];
//		      if (j < rho.size() - 1)
//			cout << ",";
//		    }
//		  cout << "]\n";
//		}
//	    }
//	  else if (line == 9)
//	    {
//	      if (value != 0.0 && value != 1.0)
//		{
//		  cout << "ERROR: Illegal value given for " << varName << " in line " << line;
//		  cout << " in file " << filename << ". Aborting.";
//		  exit(-1);
//		}
//	      int q,g;
//	      for (q = 0; q < Q; q++)
//		for (g = 0; g < G; g++)
//		  delta[q][g] = (int) value;
//
//	      cout << varName << " = " << (int) value << "\n";
//	    }
//	  else if (line == 10)
//	    {
//	      if (value <= 0.0 || value >= 1.0)
//		{
//		  cout << "ERROR: Illegal value given for " << varName << " in line " << line;
//		  cout << " in file " << filename << ". Aborting.";
//		  exit(-1);
//		}
//	      int q;
//	      for (q = 0; q < Q; q++)
//		xi[q] = value;
//	    }
//	  else if (line == 11)
//	    {
//	      if (value <= 0.0)
//		{
//		  cout << "ERROR: Illegal value given for " << varName << " in line " << line;
//		  cout << " in file " << filename << ". Aborting.";
//		  exit(-1);
//		}
//	      int q,g;
//	      for (q = 0; q < Q; q++)
//		for (g = 0; g < G; g++)
//		  sigma2[q][g] = value;
//	    }
//	  else if (line == 12)
//	    {
//	      if (value <= 0.0)
//		{
//		  cout << "ERROR: Illegal value given for " << varName << " in line " << line;
//		  cout << " in file " << filename << ". Aborting.";
//		  exit(-1);
//		}
//	      int q;
//	      for (q = 0; q < Q; q++) t[q] = value;
//	    }
//	  else if (line == 13)
//	    {
//	      if (value <= 0.0)
//		{
//		  cout << "ERROR: Illegal value given for " << varName << " in line " << line;
//		  cout << " in file " << filename << ". Aborting.";
//		  exit(-1);
//		}
//	      int q;
//	      for (q = 0; q < Q; q++) l[q] = value;
//	    }
//	  else if (line == 14)
//	    {
//	      if (value <= 0.0)
//		{
//		  cout << "ERROR: Illegal value given for " << varName << " in line " << line;
//		  cout << " in file " << filename << ". Aborting.";
//		  exit(-1);
//		}
//	      int q,g;
//	      for (q = 0; q < Q; q++)
//		for (g = 0; g < G; g++)
//		  phi[q][g] = value;
//	    }
//	  else if (line == 15)
//	    {
//	      if (value <= 0.0)
//		{
//		  cout << "ERROR: Illegal value given for " << varName << " in line " << line;
//		  cout << " in file " << filename << ". Aborting.";
//		  exit(-1);
//		}
//	      int q;
//	      for (q = 0; q < Q; q++) theta[q] = value;
//	    }
//	  else if (line == 16)
//	    {
//	      if (value <= 0.0)
//		{
//		  cout << "ERROR: Illegal value given for " << varName << " in line " << line;
//		  cout << " in file " << filename << ". Aborting.";
//		  exit(-1);
//		}
//	      int q;
//	      for (q = 0; q < Q; q++) lambda[q] = value;
//	    }
//	  else
//	    {
//	      tau2Specified = 1;
//	      if (value <= 0.0)
//		{
//		  cout << "ERROR: Illegal value given for " << varName << " in line " << line;
//		  cout << " in file " << filename << ". Aborting.";
//		  exit(-1);
//		}
//	      tau2R[line - 16] = value;
//	    }
//
//	  if (line <= 16 && line != 7 && line != 8 && line != 9)
//	    cout << varName << " = " << value << "\n";
//	}
//    }
//
//  if (tau2Specified)
//    {
//      double prod = 1.0;
//      int q;
//      for (q = 0; q < Q; q++)
//	prod *= tau2R[q];
//      for (q = 0; q < Q; q++)
//	tau2R[q] /= exp(log(prod) / Q);
//
//      for (q = 0; q < Q; q++)
//	cout << "tau2R[" << q+1 << "] = " << tau2R[q] << "\n";
//    }
//
//  int q;
//  for (q = 0; q < Q; q++)
//    tau2Rho[q] = tau2R[q];
//
//  return;
//}



void Structure::setInitialValues(double *Nu,double *DDelta,double *A,
				 double *B,double *C2,double *Gamma2,
				 double *R,double *Rho,int *Delta,
				 double *Xi,double *Sigma2,double *T,
				 double *L,double *Phi,double *Theta,
				 double *Lambda,double *Tau2R,double *Tau2Rho)
{
  int nr;

  int q,g;
  nr = 0;
  for (g = 0; g < G; g++)
    for (q = 0; q < Q; q++)
      {
	nu[q][g] = Nu[nr];
	nr++;
      }

  nr = 0;
  for (g = 0; g < G; g++)
    for (q = 0; q < Q; q++)
      {
	this->Delta[q][g] = DDelta[nr];
	nr++;
      }

  nr = 0;
  for (q = 0; q < Q; q++)
    {
      a[q] = A[nr];
      nr++;
    }

  nr = 0;
  for (q = 0; q < Q; q++)
    {
      b[q] = B[nr];
      nr++;
    }

  c2 = *C2;
  gamma2 = *Gamma2;

  nr = 0;
  int i,j;
  for (i = 0; i < Q; i++)
    for (j = i+1; j < Q; j++)
      {
	r[i][j] = R[nr];
	r[j][i] = R[nr];
	nr++;
      }


  nr = 0;
  for (i = 0; i < Q; i++)
    for (j = i+1; j < Q; j++)
      {
	rho[i][j] = Rho[nr];
	rho[j][i] = Rho[nr];
	nr++;
      }


  nr = 0;
  for (g = 0; g < G; g++)
    for (q = 0; q < Q; q++)
      {
	delta[q][g] = Delta[nr];
	nr++;
      }

  for (q = 0; q < Q; q++)
    xi[q] = Xi[q];

  nr = 0;
  for (g = 0; g < G; g++)
    for (q = 0; q < Q; q++)
      {
	sigma2[q][g] = Sigma2[nr];
	nr++;
      }


  nr = 0;
  for (q = 0; q < Q; q++)
    {
      t[q] = T[nr];
      nr++;
    }


  nr = 0;
  for (q = 0; q < Q; q++)
    {
      l[q] = L[nr];
      nr++;
    }


  nr = 0;
  for (g = 0; g < G; g++)
    for (q = 0; q < Q; q++)
      {
	phi[q][g] = Phi[nr];
	nr++;
      }


  nr = 0;
  for (q = 0; q < Q; q++)
    {
      theta[q] = Theta[nr];
      nr++;
    }


  nr = 0;
  for (q = 0; q < Q; q++)
    {
      lambda[q] = Lambda[nr];
      nr++;
    }


  nr = 0;
  for (q = 0; q < Q; q++)
    {
      tau2R[q] = Tau2R[nr];
      nr++;
    }

  nr = 0;
  for (q = 0; q < Q; q++)
    {
      tau2Rho[q] = Tau2Rho[nr];
      nr++;
    }

  return;
}




void Structure::setFinalValues(double *Nu,double *DDelta,double *A,
			       double *B,double *C2,double *Gamma2,
			       double *R,double *Rho,int *Delta,
			       double *Xi,double *Sigma2,double *T,
			       double *L,double *Phi,double *Theta,
			       double *Lambda,double *Tau2R,
			       double *Tau2Rho) const
{
  int nr;

  int q,g;
  nr = 0;
  for (g = 0; g < G; g++)
    for (q = 0; q < Q; q++)
      {
	Nu[nr] = nu[q][g];
	nr++;
      }

  nr = 0;
  for (g = 0; g < G; g++)
    for (q = 0; q < Q; q++)
      {
	DDelta[nr] = this->Delta[q][g];
	nr++;
      }

  nr = 0;
  for (q = 0; q < Q; q++)
    {
      A[nr] = a[q];
      nr++;
    }

  nr = 0;
  for (q = 0; q < Q; q++)
    {
      B[nr] = b[q];
      nr++;
    }

  *C2 = c2;
  *Gamma2 = gamma2;

  nr = 0;
  int i,j;
  for (i = 0; i < Q; i++)
    for (j = i+1; j < Q; j++)
      {
	R[nr] = r[i][j];
	nr++;
      }


  nr = 0;
  for (i = 0; i < Q; i++)
    for (j = i+1; j < Q; j++)
      {
	Rho[nr] = rho[i][j];
	nr++;
      }


  nr = 0;
  for (g = 0; g < G; g++)
    for (q = 0; q < Q; q++)
      {
	Delta[nr] = delta[q][g];
	nr++;
      }

  for (q = 0; q < Q; q++)
    Xi[q] = xi[q];

  nr = 0;
  for (g = 0; g < G; g++)
    for (q = 0; q < Q; q++)
      {
	Sigma2[nr] = sigma2[q][g];
	nr++;
      }


  nr = 0;
  for (q = 0; q < Q; q++)
    {
      T[nr] = t[q];
      nr++;
    }


  nr = 0;
  for (q = 0; q < Q; q++)
    {
      L[nr] = l[q];
      nr++;
    }


  nr = 0;
  for (g = 0; g < G; g++)
    for (q = 0; q < Q; q++)
      {
	Phi[nr] = phi[q][g];
	nr++;
      }


  nr = 0;
  for (q = 0; q < Q; q++)
    {
      Theta[nr] = theta[q];
      nr++;
    }


  nr = 0;
  for (q = 0; q < Q; q++)
    {
      Lambda[nr] = lambda[q];
      nr++;
    }


  nr = 0;
  for (q = 0; q < Q; q++)
    {
      Tau2R[nr] = tau2R[q];
      nr++;
    }

  nr = 0;
  for (q = 0; q < Q; q++)
    {
      Tau2Rho[nr] = tau2Rho[q];
      nr++;
    }

  return;
}




void Structure::setNumberOfUpdates(string &filename,vector<int> &nUpdate,vector<Update *> &update) const
{
  ifstream in;
  in.open(filename.c_str());
//  if (in.fail())
//    {
//      cout << "ERROR: Unable to open file " << filename.c_str() << ". Aborting.\n\n";
//      exit(-1);
//    }

  //
  // read one and one line
  //

  int line;
  for (line = 1; line <= 18; line++)
    {
      int bufSize = 1000;
      char buf[bufSize];
      in.get(buf,bufSize,'\n');
      char c;
//      if (in.get(c) && c != '\n')
//	{
//	  cout << "ERROR: Line " << line << " in file " << filename << " is too long. Maximum line length is ";
//	  cout << bufSize << ". Aborting.\n\n";
//	  exit(-1);
//	}

      char var1[bufSize];
      char var2[bufSize];
      int nRead = 0;
      if (line == 6 || line == 9 || line == 10)
	nRead = sscanf(buf,"%s",var1);
      else
	nRead = sscanf(buf,"%s %s",var1,var2);
      if (nRead < 1)
	{
	  if (in.eof()) return;
//	  cout << "ERROR: Error when reading line " << line << " in file " << filename << ". Aborting.\n\n";
//	  exit(-1);
	}
      string varName;
      if (line == 1)
	varName = "nu";
      else if (line == 2)
	varName = "Delta";
      else if (line == 3)
	varName = "a";
      else if (line == 4)
	varName = "b";
      else if (line == 5)
	varName = "c2";
      else if (line == 6)
	varName = "gamma2";
      else if (line == 7)
	varName = "r";
      else if (line == 8)
	varName = "rho";
      else if (line == 9)
	varName = "delta";
      else if (line == 10)
	varName = "xi";
      else if (line == 11)
	varName = "sigma2";
      else if (line == 12)
	varName = "t";
      else if (line == 13)
	varName = "l";
      else if (line == 14)
	varName = "phi";
      else if (line == 15)
	varName = "theta";
      else if (line == 16)
	varName = "lambda";
      else if (line == 17)
	varName = "tau2R";
      else
	varName = "tau2Rho";


      if (var1[0] != '=')
	{
	  int value = 0;
	  int nRead = sscanf(var1,"%d",&value);
//	  if (nRead != 1)
//	    {
//	      cout << "ERROR: Error when reading " << varName << " from file " << filename << ". Aborting.\n\n";
//	      exit(-1);
//	    }
//	  if (value < 0)
//	    {
//	      cout << "ERROR: Illegal value specified in line " << line << " in file " << filename << ". Aborting.\n\n";
//	      exit(-1);
//	    }
	  nUpdate[line - 1] = value;
	}

      if (line != 6 && line != 9 && line != 10)
	{
	  if (nRead < 2)
	    {
	      if (in.eof()) return;
//	      cout << "ERROR: Error when reading line " << line << " in file " << filename << ". Aborting.\n\n";
//	      exit(-1);
	    }

	  if (var2[0] != '=')
	    {
	      double value = 0.0;
	      int nRead = sscanf(var2,"%le",&value);
//	      if (nRead != 1)
//		{
//		  cout << "ERROR: Error when reading " << varName << " from file " << filename << ". Aborting.\n\n";
//		  exit(-1);
//		}
//	      if (value <= 0)
//		{
//		  cout << "ERROR: Illegal value specified in line " << line << " in file " << filename << ". Aborting.\n\n";
//		  exit(-1);
//		}
	      update[line - 1]->setEpsilon(value);
	    }
	}
    }

  return;
}




void Structure::setNumberOfUpdates(int *InUpdate,double *Iepsilon,
				   vector<int> &nUpdate,vector<Update *> &update)
{
  int l;
  for (l = 0; l < 18; l++)
    {
      nUpdate[l] = InUpdate[l];
      update[l]->setEpsilon(Iepsilon[l]);
    }

  return;
}



ReportDiffexpressed *Structure::setReports(string &filename,int &nBetweenReport,vector<Potential *> &potential,
			   vector<Update *> &update,vector<Report *> &report,int oneDelta)
{
  ReportDiffexpressed *reportDiffexpressed = NULL;

  ifstream in;
  in.open(filename.c_str());
//  if (in.fail())
//    {
//      cout << "ERROR: Unable to open file " << filename.c_str() << ". Aborting.\n\n";
//      exit(-1);
//    }

  //
  // read one and one line
  //

  int line;
  for (line = 1; line <= 22; line++)
    {
      int bufSize = 1000;
      char buf[bufSize];
      in.get(buf,bufSize,'\n');
      char c;
//      if (in.get(c) && c != '\n')
//	{
//	  cout << "ERROR: Line " << line << " in file " << filename << " is too long. Maximum line length is ";
//	  cout << bufSize << ". Aborting.\n\n";
//	  exit(-1);
//	}

      char var[bufSize];
      int nRead = sscanf(buf,"%s",var);
      if (nRead != 1)
	{
	  if (in.eof()) return NULL;
//	  cout << "ERROR: Error when reading line " << line << " in file " << filename << ". Aborting.\n\n";
//	  exit(-1);
	}
      string varName;
      if (line == 1)
	varName = "number of iterations between each output";
      else if (line == 2)
	varName = "potential values";
      else if (line == 3)
	varName = "acceptance rates";
      else if (line == 4)
	varName = "nu";
      else if (line == 5)
	varName = "Delta";
      else if (line == 6)
	varName = "a";
      else if (line == 7)
	varName = "b";
      else if (line == 8)
	varName = "c2";
      else if (line == 9)
	varName = "gamma2";
      else if (line == 10)
	varName = "r";
      else if (line == 11)
	varName = "rho";
      else if (line == 12)
	varName = "delta";
      else if (line == 13)
	varName = "xi";
      else if (line == 14)
	varName = "sigma2";
      else if (line == 15)
	varName = "t";
      else if (line == 16)
	varName = "l";
      else if (line == 17)
	varName = "phi";
      else if (line == 18)
	varName = "theta";
      else if (line == 19)
	varName = "lambda";
      else if (line == 20)
	varName = "tau2";
      else if (line == 21)
	varName = "probDelta";
      else
	varName = "diffexpressed";


      if (var[0] != '=')
	{
	  if (line == 1)
	    {
	      int value = 0;
	      nRead = sscanf(var,"%d",&value);
//	      if (nRead != 1)
//		{
//		  cout << "ERROR: Error when reading " << varName << " from file " << filename << ". Aborting.\n\n";
//		  exit(-1);
//		}
//	      if (value <= 0)
//		{
//		  cout << "ERROR: Illegal value specified in line " << line << " in file " << filename << ". Aborting.\n\n";
//		  exit(-1);
//		}
	      nBetweenReport = value;
	    }
	  else
	    {
	      char name[bufSize];
	      nRead = sscanf(var,"%s",name);
//	      if (nRead != 1)
//		{
//		  cout << "ERROR: Error when reading filename for " << varName << " from file " << filename << ". Aborting.\n\n";
//		  exit(-1);
//		}
	      string value(name);
	      if (line == 2)
		report.push_back(new ReportPotential(value,potential));
	      else if (line == 3)
		report.push_back(new ReportAcceptance(value,update));
	      else if (line == 4)
		report.push_back(new ReportNu(value));
	      else if (line == 5)
		report.push_back(new ReportDDelta(value));
	      else if (line == 6)
		report.push_back(new ReportA(value));
	      else if (line == 7)
		report.push_back(new ReportB(value));
	      else if (line == 8)
		report.push_back(new ReportC2(value));
	      else if (line == 9)
		report.push_back(new ReportGamma2(value));
	      else if (line == 10)
		report.push_back(new ReportR(value));
	      else if (line == 11)
		report.push_back(new ReportRho(value));
	      else if (line == 12)
		report.push_back(new ReportDelta(value));
	      else if (line == 13)
		report.push_back(new ReportXi(value));
	      else if (line == 14)
		report.push_back(new ReportSigma2(value));
	      else if (line == 15)
		report.push_back(new ReportT(value));
	      else if (line == 16)
		report.push_back(new ReportL(value));
	      else if (line == 17)
		report.push_back(new ReportPhi(value));
	      else if (line == 18)
		report.push_back(new ReportTheta(value));
	      else if (line == 19)
		report.push_back(new ReportLambda(value));
	      else if (line == 20)
		report.push_back(new ReportTau2R(value));
	      else if (line == 21)
		report.push_back(new ReportTau2Rho(value));
	      else if (line == 22)
		report.push_back(new ReportProbDelta(value,this,oneDelta));
	      else if (line == 23)
		{
		  reportDiffexpressed = new ReportDiffexpressed(value,this);
		  report.push_back(reportDiffexpressed);
		}
	    }
	}
    }

  return reportDiffexpressed;
}




ReportDiffexpressed *Structure::setReports(int *output,int &nBetweenReport,int writeToFile,
			   char **directory,vector<Potential *> &potential,
			   vector<Update *> &update,vector<Report *> &report,
			   double *valuePotential,double *valueAcceptance,
			   double *valueNu,double *valueDDelta,double *valueA,
			   double *valueB,double *valueC2,double *valueGamma2,
			   double *valueR,double *valueRho,int *valueDelta,
			   double *valueXi,double *valueSigma2,double *valueT,
			   double *valueL,double *valuePhi,double *valueTheta,
			   double *valueLambda,double *valueTau2R,
			   double *valueTau2Rho,double *valueProbDelta,
			   double *valueDiffexpressed,
			   int *writeDiffexpressedTofile,int oneDelta)
{
  int nr = 0;
  ReportDiffexpressed *reportDiffexpressed = NULL;

  nBetweenReport = output[nr];
  nr++;

  if (writeToFile == 0)
    {
      if (output[nr]) report.push_back(new ReportPotential(valuePotential,potential));
      nr++;

      if (output[nr]) report.push_back(new ReportAcceptance(valueAcceptance,update));
      nr++;

      if (output[nr]) report.push_back(new ReportNu(valueNu));
      nr++;

      if (output[nr]) report.push_back(new ReportDDelta(valueDDelta));
      nr++;

      if (output[nr]) report.push_back(new ReportA(valueA));
      nr++;

      if (output[nr]) report.push_back(new ReportB(valueB));
      nr++;

      if (output[nr]) report.push_back(new ReportC2(valueC2));
      nr++;

      if (output[nr]) report.push_back(new ReportGamma2(valueGamma2));
      nr++;

      if (output[nr]) report.push_back(new ReportR(valueR));
      nr++;

      if (output[nr]) report.push_back(new ReportRho(valueRho));
      nr++;

      if (output[nr]) report.push_back(new ReportDelta(valueDelta));
      nr++;

      if (output[nr]) report.push_back(new ReportXi(valueXi));
      nr++;

      if (output[nr]) report.push_back(new ReportSigma2(valueSigma2));
      nr++;

      if (output[nr]) report.push_back(new ReportT(valueT));
      nr++;

      if (output[nr]) report.push_back(new ReportL(valueL));
      nr++;

      if (output[nr]) report.push_back(new ReportPhi(valuePhi));
      nr++;

      if (output[nr]) report.push_back(new ReportTheta(valueTheta));
      nr++;

      if (output[nr]) report.push_back(new ReportLambda(valueLambda));
      nr++;

      if (output[nr]) report.push_back(new ReportTau2R(valueTau2R));
      nr++;

      if (output[nr]) report.push_back(new ReportTau2Rho(valueTau2Rho));
      nr++;

      if (output[nr]) report.push_back(new ReportProbDelta(valueProbDelta,this,oneDelta));
      nr++;

      if (output[nr])
	{
	  reportDiffexpressed = new ReportDiffexpressed(valueDiffexpressed,this);
	  report.push_back(reportDiffexpressed);
	}
      if (*writeDiffexpressedTofile == 1)
	{
	  string dir(*directory);
	  string filename;

	  filename.empty();
	  filename = "diffExpressed.log";
	  filename = dir + filename;
	  reportDiffexpressed = new ReportDiffexpressed(filename,this);
	  report.push_back(reportDiffexpressed);
	}
    }
  else
    {
      string dir(*directory);
      string filename;

      filename.empty();

      filename = "potential.log";
      filename = dir + filename;
      if (output[nr]) report.push_back(new ReportPotential(filename,potential));
      nr++;

      filename.empty();
      filename = "acceptance.log";
      filename = dir + filename;
      if (output[nr]) report.push_back(new ReportAcceptance(filename,update));
      nr++;

      filename.empty();
      filename = "nu.log";
      filename = dir + filename;
      if (output[nr]) report.push_back(new ReportNu(filename));
      nr++;

      filename.empty();
      filename = "DDelta.log";
      filename = dir + filename;
      if (output[nr]) report.push_back(new ReportDDelta(filename));
      nr++;

      filename.empty();
      filename = "a.log";
      filename = dir + filename;
      if (output[nr]) report.push_back(new ReportA(filename));
      nr++;

      filename.empty();
      filename = "b.log";
      filename = dir + filename;
      if (output[nr]) report.push_back(new ReportB(filename));
      nr++;

      filename.empty();
      filename = "c2.log";
      filename = dir + filename;
      if (output[nr]) report.push_back(new ReportC2(filename));
      nr++;

      filename.empty();
      filename = "gamma2.log";
      filename = dir + filename;
      if (output[nr]) report.push_back(new ReportGamma2(filename));
      nr++;

      filename.empty();
      filename = "r.log";
      filename = dir + filename;
      if (output[nr]) report.push_back(new ReportR(filename));
      nr++;

      filename.empty();
      filename = "rho.log";
      filename = dir + filename;
      if (output[nr]) report.push_back(new ReportRho(filename));
      nr++;

      filename.empty();
      filename = "delta.log";
      filename = dir + filename;
      if (output[nr]) report.push_back(new ReportDelta(filename));
      nr++;

      filename.empty();
      filename = "xi.log";
      filename = dir + filename;
      if (output[nr]) report.push_back(new ReportXi(filename));
      nr++;

      filename.empty();
      filename = "sigma2.log";
      filename = dir + filename;
      if (output[nr]) report.push_back(new ReportSigma2(filename));
      nr++;

      filename.empty();
      filename = "t.log";
      filename = dir + filename;
      if (output[nr]) report.push_back(new ReportT(filename));
      nr++;

      filename.empty();
      filename = "l.log";
      filename = dir + filename;
      if (output[nr]) report.push_back(new ReportL(filename));
      nr++;

      filename.empty();
      filename = "phi.log";
      filename = dir + filename;
      if (output[nr]) report.push_back(new ReportPhi(filename));
      nr++;

      filename.empty();
      filename = "theta.log";
      filename = dir + filename;
      if (output[nr]) report.push_back(new ReportTheta(filename));
      nr++;

      filename.empty();
      filename = "lambda.log";
      filename = dir + filename;
      if (output[nr]) report.push_back(new ReportLambda(filename));
      nr++;

      filename.empty();
      filename = "tau2R.log";
      filename = dir + filename;
      if (output[nr]) report.push_back(new ReportTau2R(filename));
      nr++;

      filename.empty();
      filename = "tau2Rho.log";
      filename = dir + filename;
      if (output[nr]) report.push_back(new ReportTau2Rho(filename));
      nr++;

      filename.empty();
      filename = "probDelta.log";
      filename = dir + filename;
      if (output[nr]) report.push_back(new ReportProbDelta(filename,this,oneDelta));
      nr++;

      filename.empty();
      filename = "diffExpressed.log";
      filename = dir + filename;
      if (output[nr])
	{
	  reportDiffexpressed = new ReportDiffexpressed(filename,this);
	  report.push_back(reportDiffexpressed);
	}
      nr++;
    }

  return reportDiffexpressed;
}




//void Structure::readCorrelationValues(vector<vector<double> > &corr,const char *line,int bufSize,
//				      const string varName,string &filename) const
//{
//  char text[bufSize];
//  strncpy(text,line,bufSize);
//  istrstream inputLine(text,bufSize);
//
//  int Q = corr.size();
//
//  double value;
//  inputLine >> value;
//  if (inputLine.fail())
//    {
//      cout << "ERROR: Error when reading " << varName << " from file " << filename << ". Aborting.\n\n";
//      exit(-1);
//    }
//
//  int i,j;
//  for (i = 0; i < Q; i++)
//    for (j = i+1; j < Q; j++)
//      {
//	corr[i][j] = value;
//	corr[j][i] = value;
//
//	//
//	// read next value
//	//
//
//	double vv;
//	inputLine >> vv;
//	if (!(inputLine.fail()))
//	  value = vv;
//      }
//
//  //
//  // check that the matrix is positive definite and all elements are non-negative
//  //
//
//  int isneg = 0;
//  for (i = 0; i < Q; i++)
//    for (j = 0; j < Q; j++)
//      isneg += (corr[i][j] < 0.0);
//  if (isneg > 0)
//    {
//      cout << "ERROR: Initial matrix specified for " << varName << " in file " << filename << " has negative values.\n";
//      cout << "Matrix specified:\n";
//      for (i = 0; i < Q; i++)
//	{
//	  cout << "[";
//	  for (j = 0; j < Q; j++)
//	    {
//	      cout << corr[i][j];
//	      if (j < Q - 1)
//		cout << ",";
//	    }
//	  cout << "]\n";
//	}
//      cout << "Aborting.\n\n";
//      exit(-1);
//    }
//
//  int err = 0;
//  vector<vector<double> > mat;
//  mat.resize(Q);
//  for (i = 0; i < Q; i++)
//    mat[i].resize(Q);
//  for (i = 0; i < Q; i++)
//    for (j = 0; j < Q; j++)
//      mat[i][j] = corr[i][j];
//  Cholesky chol(mat,err);
//  if (err != 0)
//    {
//      cout << "ERROR: Initial matrix specified for " << varName << " in file " << filename << " is not positive definite.\n";
//      cout << "Matrix specified:\n";
//      for (i = 0; i < Q; i++)
//	{
//	  cout << "[";
//	  for (j = 0; j < Q; j++)
//	    {
//	      cout << corr[i][j];
//	      if (j < Q - 1)
//		cout << ",";
//	    }
//	  cout << "]\n";
//	}
//      cout << "Aborting.\n\n";
//      exit(-1);
//    }
//
//
//  return;
//}
