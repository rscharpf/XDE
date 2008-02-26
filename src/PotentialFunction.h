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
     
     
     
