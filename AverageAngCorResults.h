#include <TVector3.h>
#include <fstream>
#include <iostream>
#include <TLorentzVector.h>
#include <TRandom3.h>
#include <string>
#include <TH1.h>
#include <TFile.h>
#include <TTree.h>

int NumberThetaAlphaPoints = 40;
int NumberPhiAlphaPoints = 180;
int NumberThetaDecayPoints = 180;
int NumberOfCHUCK3Angles = 100;
int NumberMonteCarloEvents = 5000;

float DeltaThetaAlpha = 0.1; //Size of the ThetaAlpha steps

double Ex = 11.864; //MeV

double **ReadCrossSectionTable(char *InputFileName);
double CrossSectionCalculation(double **CrossSectionTable, double ThetaAlpha);
double CalculateThetaAlpha(int Iteration){return DeltaThetaAlpha * Iteration;}
TLorentzVector* CalculateEjectileRecoilVectors(double *Masses, double TBeam, double Ex, double ThetaAlphaCM, double PhiAlphaCM);
double*** ReadAngCorTable(char *InputFileName);
TH1F* MakeAngCorHistogram(double ***AngCorTable, double ThetaCM, double PhiAlphaCM);
TLorentzVector* CalculateDecayResidualVectors(double *Masses, double Ex,TLorentzVector Recoil4Vector,double DecayTheta,double DecayPhi);
bool TestKinematics(double *Masses, double TBeam, double Ex, TLorentzVector* KinematicVectors, TLorentzVector* KinematicVectorsDecay);