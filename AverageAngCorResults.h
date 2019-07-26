#include <TVector3.h>
#include <fstream>
#include <iostream>
#include <TLorentzVector.h>
#include <TRandom3.h>
#include <string>
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TGraph.h>

double ThetaAlphaStartAngle = 0.;

int NumberThetaAlphaPoints = 40;
int NumberPhiAlphaPoints = 360;
int NumberThetaDecayPoints = 180;
int NumberPhiDecayPoints = 180;
int NumberOfCHUCK3Angles = 100;
int NumberMonteCarloEvents = 100;

float DeltaThetaAlpha = 0.1; //Size of the ThetaAlpha steps

double Ex = 11.864; //MeV

double **ReadCrossSectionTable(char *InputFileName);
double CrossSectionCalculation(double **CrossSectionTable, double ThetaAlpha);
double CalculateThetaAlpha(int Iteration){return DeltaThetaAlpha * Iteration + ThetaAlphaStartAngle;}
TLorentzVector* CalculateEjectileRecoilVectors(double *Masses, double TBeam, double Ex, double ThetaAlphaCM, double PhiAlphaCM);
double*** ReadAngCorTable(char *InputFileName);
TH2F* MakeAngCorHistogram(double ***AngCorTable, double ThetaCM);
TLorentzVector* CalculateDecayResidualVectors(double *Masses, double Ex,TLorentzVector Recoil4Vector,double DecayTheta,double DecayPhi);
bool TestKinematics(double *Masses, double TBeam, double Ex, TLorentzVector* KinematicVectors, TLorentzVector* KinematicVectorsDecay);

TVector3 CollisionalCoMBoost(double *Masses, double TBeam);