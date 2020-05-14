#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <TMath.h>
#include <TLorentzVector.h>
#include <TVector3.h>

double TBeam = 200, ExInitial = 11.728, ExFinal = 0; //MeV

double DoSimpleKinematicConversion(double *Masses, double TBeam, double ExInitial, double ExFinal, double ThetaDecayCM, double PhiDecayCM);

int main(int argc, char *argv[])
{
    if(argc!=2) 
    {
        fprintf(stderr, "usage: mcerr <filename for output>\n");
        exit(1);
    }
    
    std::cout << "Output file name: " << argv[1] << std::endl;
    
    double *Masses = new double[6];
    Masses[0] = 3728.400952;//4He
    Masses[1] = 22341.92265;//24Mg
    Masses[2] = Masses[0];//4He
    Masses[3] = Masses[1];//24Mg - recoil
    Masses[4] = Masses[0];//4He again but could be another decay particle
    Masses[5] = 18622.83825;//20Ne
    
    TFile *fout = new TFile(argv[1],"RECREATE");
    TTree *trout = new TTree("AngCorData","AngCorData");
    
    double ThetaDecayCM = 0, ThetaDecayLab = 0;
    double JacobeanFactor = 0;
    
    trout->Branch("ThetaDecayLab",&ThetaDecayLab);//
    trout->Branch("ThetaDecayCM",&ThetaDecayCM);//
    trout->Branch("ExInitial",&ExInitial);//
    trout->Branch("ExFinal",&ExFinal);
    trout->Branch("JacobeanFactor",&JacobeanFactor);
    
    
    for(ThetaDecayCM = 1;ThetaDecayCM<180;ThetaDecayCM+=1)
    {
        ThetaDecayLab = DoSimpleKinematicConversion(Masses, TBeam, ExInitial, ExFinal, ThetaDecayCM, 0);
        
        //these are the lab angles
        double ThetaDecayLower = DoSimpleKinematicConversion(Masses, TBeam, ExInitial, ExFinal, ThetaDecayCM - 0.5, 0);
        double ThetaDecayUpper = DoSimpleKinematicConversion(Masses, TBeam, ExInitial, ExFinal, ThetaDecayCM + 0.5, 0);
        
        JacobeanFactor = sin(ThetaDecayCM*TMath::Pi()/180.)/sin(ThetaDecayLab*TMath::Pi()/180.)/(ThetaDecayUpper - ThetaDecayLower);
        
        trout->Fill();
    }
    
    trout->Write();
    fout->Close();
    
    return 0;
}

double DoSimpleKinematicConversion(double *Masses, double TBeam, double ExInitial, double ExFinal, double ThetaDecayCM, double PhiDecayCM)
{
    double RecoilMass = Masses[3] + ExInitial; //Convert the recoil mass to the right value including the excitation-energy dependence
    
    double s = pow(Masses[0],2.) + pow(Masses[1],2.) + 2 * Masses[1] * (TBeam + Masses[0]);
    
    TLorentzVector f4MomentumCentreOfMass(0,0,0,sqrt(s));
    
    double ECM0 = (s + pow(Masses[0],2.) - pow(Masses[1],2.))/(2*sqrt(s));
    double ECM1 = (s + pow(Masses[1],2.) - pow(Masses[0],2.))/(2*sqrt(s));
    double ECM2 = (s + pow(Masses[2],2.) - pow(RecoilMass,2.))/(2*sqrt(s));
    double ECM3 = (s + pow(RecoilMass,2.) - pow(Masses[2],2.))/(2*sqrt(s));
    
    double PCM0 = sqrt(pow(ECM0,2.) - pow(Masses[0],2.));
    double PCM1 = sqrt(pow(ECM1,2.) - pow(Masses[1],2.));
    double PCM2 = sqrt(pow(ECM2,2.) - pow(Masses[2],2.));
    double PCM3 = sqrt(pow(ECM3,2.) - pow(RecoilMass,2.));
    
    TVector3 Lab3Momentum0(0,0,sqrt(pow(TBeam,2.) + 2 * TBeam * Masses[0]));
    TVector3 Lab3Momentum1(0,0,0);
    
    TLorentzVector Lab4Momentum0(Lab3Momentum0,Masses[0] + TBeam);
    TLorentzVector Lab4Momentum1(Lab3Momentum1,Masses[1]);
    
    double BetaCM = (Lab4Momentum0+Lab4Momentum1).Beta();
//     std::cout << "BetaCM = " << BetaCM << std::endl;
    
    TLorentzVector CoM4Momentum0 = Lab4Momentum0;
    CoM4Momentum0.Boost(TVector3(0,0,-BetaCM));
    
    TLorentzVector CoM4Momentum1 = Lab4Momentum1;
    CoM4Momentum1.Boost(TVector3(0,0,-BetaCM));
    
    TLorentzVector CoM4Momentum2 = TLorentzVector(PCM2,0,0,ECM2);
    
    TLorentzVector CoM4Momentum3 = CoM4Momentum0 + CoM4Momentum1 - CoM4Momentum2;
    
    double BetaRecoil = CoM4Momentum3.Beta();
//     std::cout << "BetaRecoil = " << BetaRecoil << std::endl;
    
    //Assume recoil is stationary in the centre-of-mass frame
//     std::cout << "Threshold for the decay: " << Masses[3] - Masses[4] - Masses[5] << std::endl;
    double TDecayParticle = Masses[5]/(Masses[4] + Masses[5]) * (ExInitial + (Masses[3] - Masses[4] - Masses[5]) - ExFinal);//kinetic energy of the decay particle
//     std::cout << "TDecayParticle: " << TDecayParticle << std::endl;
    
    
    double pDecayParticle = sqrt(TDecayParticle * (TDecayParticle + 2*Masses[4]));
    TVector3 DecayParticle3Momentum(pDecayParticle*sin(ThetaDecayCM*TMath::Pi()/180.)*cos(PhiDecayCM*TMath::Pi()/180.),
                           pDecayParticle*sin(ThetaDecayCM*TMath::Pi()/180.)*sin(PhiDecayCM*TMath::Pi()/180.),
                           pDecayParticle*cos(ThetaDecayCM*TMath::Pi()/180.));
    
    TLorentzVector DecayParticle4Momentum(DecayParticle3Momentum, TDecayParticle + Masses[4]);
 
    DecayParticle4Momentum.Boost(TVector3(0,0,-BetaRecoil));
    
    DecayParticle4Momentum.Boost(TVector3(0,0,BetaCM));
    
    return DecayParticle4Momentum.Theta()*180./TMath::Pi();
}
