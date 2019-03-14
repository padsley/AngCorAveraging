#include "AverageAngCorResults.h"

double *ReadCrossSectionTable(char InputFileName);

bool VerboseFlag = false;

//usage AverageAngCorResults CHUCK3CrossSectionFilePath PathToAngCorOutputFile

int main(int argc, char *argv[])
{
    TH1F *AngularCorrelationHistogram = new TH1F("AngularCorrelationHistogram","",181,0,180);
    
    if(argc!=4) 
    {
        fprintf(stderr, "usage: mcerr <filename for cross section> <filename for AngCor results> <filename for output>\n");
        exit(1);
    }
    std::cout << "Using cross section file: " << argv[1] << std::endl;
    std::cout << "Using AngCor input file: " << argv[2] << std::endl;
    std::cout << "Output file name: " << argv[3] << std::endl;
    
    double *Masses = new double[6];
    Masses[0] = 3728.400952;//4He
    Masses[1] = 22341.92265;//24Mg
    Masses[2] = Masses[0];//4He
    Masses[3] = Masses[1];//24Mg - recoil
    Masses[4] = Masses[0];//4He again but could be another decay particle
    Masses[5] = 18622.83825;//20Ne
    
    TRandom3 *randy = new TRandom3(0);
    
//     TLorentzVector *KinematicVectors, *KinematicVectorsDecay;
    
    //Read in CHUCK3 differential cross sections
    double **CrossSectionTable = ReadCrossSectionTable(argv[1]);
    
    //Read in the AngCor data table
    double ***AngCorTable = ReadAngCorTable(argv[2]);
    
    std::cout << "Read in all of the CHUCK3 and AngCor information" << std::endl;
    
    TFile *fout = new TFile(argv[3],"RECREATE");
    TTree *trout = new TTree("AngCorData","AngCorData");
    
    double ApertureX = 0, ApertureY = 0; //ApertureX is ThetaSCAT in K600 parlance (I think) but has been given a different name to avoid confusion - ApertureY should be PhiSCAT but this is annulled by the K600 focussing.
    
    double ThetaAlphaLab = 0, PhiAlphaLab = 0, ThetaDecayLab = 0;
    
    trout->Branch("ApertureX",&ApertureX);
    trout->Branch("ApertureY",&ApertureY);
    trout->Branch("ThetaAlphaLab",&ThetaAlphaLab);
    trout->Branch("PhiAlphaLab",&PhiAlphaLab);
    trout->Branch("ThetaDecayLab",&ThetaDecayLab);
    
    for(int i=0;i<NumberThetaAlphaPoints;i++)//Loop over the centre-of-mass theta angles
    {
        double ThetaAlphaCM = CalculateThetaAlpha(i);
        if(VerboseFlag)std::cout << "ThetaAlphaCM: " << ThetaAlphaCM << std::endl;
        
        double CrossSectionValue = CrossSectionCalculation(CrossSectionTable, ThetaAlphaCM);
        if(VerboseFlag)std::cout << "CrossSectionValue: " << CrossSectionValue << std::endl;
        for(int j=0;j<NumberPhiAlphaPoints;j++)//Loop over the centre-of-mass phi angles
        {
            double PhiAlphaCM = (double)j;
            
            //Get the Boost TVector3 for the recoil in the lab frame
            TLorentzVector *KinematicVectors = CalculateEjectileRecoilVectors(Masses, 200, Ex, ThetaAlphaCM, PhiAlphaCM);
            //KinematicVectors[0] is the kinematic vector for the ejectile
            //KinematicVectors[1] is the kinematic vector for the recoil
            ThetaAlphaLab = KinematicVectors[0].Theta()*180./TMath::Pi();
            PhiAlphaLab = KinematicVectors[0].Phi()*180./TMath::Pi();
            
            ApertureX = 180./TMath::Pi() * asin(sin(ThetaAlphaLab*TMath::Pi()/180.) * cos(PhiAlphaLab*TMath::Pi()/180.));
            ApertureY = 180./TMath::Pi() * asin(sin(ThetaAlphaLab*TMath::Pi()/180.) * sin(PhiAlphaLab*TMath::Pi()/180.));
            
            
            //Get the histogram for the 
            TH1F *hAngCor = MakeAngCorHistogram(AngCorTable, ThetaAlphaCM, PhiAlphaCM);
            
            if(KinematicVectors[0].Theta()<2.*TMath::Pi()/180.)
            {
                if(VerboseFlag)std::cout << KinematicVectors[0].Theta()*180./TMath::Pi() << "\t" << CrossSectionValue << std::endl;
                for(int n=0;n<NumberMonteCarloEvents;n++)
                {
                    //Get decay theta for the decay particle - this is according to the AngCor calculations table
                    double DecayTheta = hAngCor->GetRandom();
                    
                    //Get decay phi for the decay particle - in degrees
                    double DecayPhi = randy->Rndm() * 360;
                    
                    //The kinematic behaviour of the recoil is described by KinematicVectors[1]
                    double BetaRecoil = KinematicVectors[1].Beta();
                    TVector3 Lab3MomentumRecoil = KinematicVectors[1].Vect();
                    
                    //Calculate the decay particle and the residual particle's 4-momenta - these are in the LAB frame!!!!
                    if(Masses[3] != Masses[1]&& VerboseFlag)std::cout << "Masses[3] was permanently changed in the first kinematics loop: " << Masses[1] << "\t" << Masses[3] << "\t" << Ex << std::endl;
                    TLorentzVector *KinematicVectorsDecay = CalculateDecayResidualVectors(Masses,Ex,KinematicVectors[1],DecayTheta,DecayPhi);
                    
                    //Store the results - probably in a histogram
                    AngularCorrelationHistogram->Fill(KinematicVectorsDecay[0].Theta()*180./TMath::Pi(),CrossSectionValue*sin(ThetaAlphaCM*TMath::Pi()/180.));
                    
                    ThetaDecayLab = KinematicVectorsDecay[0].Theta()*180./TMath::Pi();
                    
                    delete [] KinematicVectorsDecay;
                    
                    trout->Fill();
                }
            }
            delete [] KinematicVectors;
            delete hAngCor;
            
        }
    }
    delete CrossSectionTable;//These aren't actually the right way of doing it. Baaaah.
    delete AngCorTable;
    
    
    AngularCorrelationHistogram->Write();
    fout->Close();
    
    return 0;
}

double **ReadCrossSectionTable(char *InputFileName)
{
    double **result = new double*[NumberOfCHUCK3Angles];
    for(int i=0;i<NumberOfCHUCK3Angles;i++)
    {
        result[i] = new double[2];
        result[i][0] = 0;
        result[i][1] = 0;
    }
    
    std::cout << "Using the cross section from file: " << InputFileName << std::endl;
    
    std::ifstream inputFile;
    inputFile.open(InputFileName);
    
    std::string str;
    
    bool ReadCrossSections = false;
    
    if(inputFile.is_open())
    {
        while(getline(inputFile,str))
        {
            if(!ReadCrossSections)
            {
                if(VerboseFlag)std::cout << str << std::endl;
                if(VerboseFlag)std::cout << "Substring: " << str.substr(0,11).data() << std::endl;
                if(strcmp(str.substr(0,11).c_str(),"0CHANNEL  2")==0)
                {
                    getline(inputFile,str); //Read in one extra line to skip the one with the "THETA" heading etc on
                    if(VerboseFlag)std::cout << "Substring2: " << str.substr(0,1) << std::endl;
                    if(strcmp(str.substr(0,1).c_str(),"+")==0)
                    {
                        std::cout << "Found CHUCK3 cross sections to load" << std::endl;
                        ReadCrossSections = true;
                    }
                }
            }
            else
            {
                double Theta = 0, CrossSection = 0;
                
                //                 getline(inputFile,str);
                if(VerboseFlag)printf("%s\n",str.c_str());
                if(strcmp(str.substr(0,1).c_str(),"*")==0)
                {
                    Theta = atof(str.substr(4,6).c_str());
                    if(VerboseFlag)std::cout << "Theta = " << Theta << std::endl;
                    CrossSection = atof(str.substr(14,10).c_str());
                    if(VerboseFlag)std::cout << "CrossSection = " << CrossSection << std::endl;
                    if(VerboseFlag)std::cout << "Index = " << (Theta/DeltaThetaAlpha) << "\t" << (int)round(Theta/DeltaThetaAlpha) << std::endl;
                    result[(int)round(Theta/DeltaThetaAlpha)][0] = Theta;
                    result[(int)round(Theta/DeltaThetaAlpha)][1] = CrossSection;
                    if(VerboseFlag)std::cout << "Done cross section for theta = " << Theta << std::endl;
                }
            }
        }
    }
    inputFile.close();
    std:: cout << "Read the cross section table" << std::endl;
    return result;
}

double*** ReadAngCorTable(char *InputFileName)
{
    std::cout << "Read AngCor Table" << std::endl;
    double*** result = new double**[NumberThetaAlphaPoints];
    for(int i=0;i<NumberThetaAlphaPoints;i++)
    {
        result[i] = new double*[NumberPhiAlphaPoints+1];
        for(int j=0;j<=NumberPhiAlphaPoints;j++)
        {
            result[i][j] = new double[NumberThetaDecayPoints+1];
            for(int k=0;k<=NumberThetaDecayPoints;k++)
            {
                result[i][j][k] = 0;
            }
        }
    }
    std::cout << "Created array for AngCor Table" << std::endl;
    std::ifstream inputFile;
    inputFile.open(InputFileName);
    
    if(inputFile.is_open())
    {
        double dummy0 = 0, dummy1 = 0, dummy2 = 0, dummy3 = 0;
        while(inputFile >> dummy0 >> dummy1 >> dummy2 >> dummy3)
        {
            if(VerboseFlag)std::cout << dummy0 << "\t" << dummy1 << "\t" << dummy2 << "\t" << dummy3 << std::endl;
            result[(int)(dummy0/DeltaThetaAlpha)][(int)dummy1][(int)dummy2] = dummy3;
        }
    }
    inputFile.close();
    std::cout << "Loaded all of the AngCor values" << std::endl;
    return result;
}

double CrossSectionCalculation(double **CrossSectionTable, double ThetaAlpha)
{
    double result = 0;
    
    if(VerboseFlag)std::cout << "CrossSection index: " << (int)(ThetaAlpha/DeltaThetaAlpha) << std::endl;
    
    result = CrossSectionTable[(int)(ThetaAlpha/DeltaThetaAlpha)][1];   
    
    if(VerboseFlag)std::cout << "Cross section value: " << result << std::endl;
    
    if(result==0)std::cout << "Not found a good cross section" << std::endl;
    return result;
}

TH1F* MakeAngCorHistogram(double ***AngCorTable, double ThetaCM, double PhiAlphaCM)
{
    TH1F *result = new TH1F("hResult","",181,0,181);
    if(VerboseFlag)std::cout << "hResult Bin Width: " << result->GetBinWidth(1) << std::endl;
    
    for(int i=0;i<=NumberThetaDecayPoints;i++)
    {
        if(VerboseFlag)std::cout << "i: " << i << "\t (int)round(ThetaCM/DeltaThetaAlpha): " << (int)round(ThetaCM/DeltaThetaAlpha) << "\t(int)PhiAlphaCM: " << (int)PhiAlphaCM << std::endl;
        result->SetBinContent(i+1,AngCorTable[(int)round(ThetaCM/DeltaThetaAlpha)][(int)PhiAlphaCM][i]);   
    }
    
    char buffer[256];
    sprintf(buffer,"figures/hAngCor%d_%d.root",(int)round(ThetaCM/DeltaThetaAlpha),(int)PhiAlphaCM);
//     result->SaveAs(buffer);
    
    return result;
}

TLorentzVector* CalculateEjectileRecoilVectors(double *Masses, double TBeam, double Ex, double ThetaAlphaCM, double PhiAlphaCM)
{
    //Function to calculate the recoil TLorentzVector in the lab frame
    TLorentzVector *result = new TLorentzVector[2];
    
    double RecoilMass = Masses[3] + Ex; //Convert the recoil mass to the right value including the excitation-energy dependence
    
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
    
    TLorentzVector CoM4Momentum0 = Lab4Momentum0;
    CoM4Momentum0.Boost(TVector3(0,0,-BetaCM));
    
    TLorentzVector CoM4Momentum1 = Lab4Momentum1;
    CoM4Momentum1.Boost(TVector3(0,0,-BetaCM));
    
    TLorentzVector CoM4Momentum2 = TLorentzVector(PCM2 * sin(ThetaAlphaCM * TMath::Pi()/180.) * cos(PhiAlphaCM * TMath::Pi()/180.),
                                                  PCM2 * sin(ThetaAlphaCM * TMath::Pi()/180.) * sin(PhiAlphaCM * TMath::Pi()/180.),
                                                  PCM2 * cos(ThetaAlphaCM * TMath::Pi()/180.),
                                                  ECM2);
    TLorentzVector CoM4Momentum3 = f4MomentumCentreOfMass - CoM4Momentum2;
    
    TLorentzVector Lab4Momentum2 = CoM4Momentum2;
    Lab4Momentum2.Boost(0,0,BetaCM);
    TLorentzVector Lab4Momentum3 = CoM4Momentum3;
    Lab4Momentum3.Boost(0,0,BetaCM);
    
    if(VerboseFlag)Lab4Momentum2.Print();
    if(VerboseFlag)Lab4Momentum3.Print();
    
    if(VerboseFlag && (Lab4Momentum0+Lab4Momentum1-Lab4Momentum2-Lab4Momentum3).Mag()>1.e-6)
    {
        std::cout << "Problem with two-body kinematics" << std::endl;
        (Lab4Momentum0+Lab4Momentum1-Lab4Momentum2-Lab4Momentum3).Print();
    }
    
    result[0] = Lab4Momentum2;
    result[1] = Lab4Momentum3;
    
    return result;
}

TLorentzVector* CalculateDecayResidualVectors(double *Masses, double Ex, TLorentzVector Recoil4Vector,double DecayTheta,double DecayPhi)
{
    TLorentzVector *result = new TLorentzVector[2];
    
    double RecoilMass = Masses[3] + Ex;
    
    DecayTheta *= TMath::Pi()/180.;//These angles relative to the beam direction
    DecayPhi *= TMath::Pi()/180.;//Relative to the beam direction
    
    //Start off in rest in the recoil frame
    
    //Energies for the breakup given by the excitation energy of the recoil in the rest frame
    TLorentzVector Initial4Momentum(0,0,0,RecoilMass);
    
    if(VerboseFlag)std::cout << "Q-value for the decay: " << Masses[4] + Masses[5] - Masses[3] << std::endl;
    if(VerboseFlag)std::cout << "Energy available for the decay: " << RecoilMass - Masses[4] - Masses[5] << std::endl;
    
    double ECM_Decay = (pow(RecoilMass,2.) + pow(Masses[4],2.) - pow(Masses[5],2.))/(2*RecoilMass);
    double ECM_Residual = (pow(RecoilMass,2.) - pow(Masses[4],2.) + pow(Masses[5],2.))/(2*RecoilMass);
    
    double PCM_Decay = sqrt(pow(ECM_Decay,2.) - pow(Masses[4],2.));
    double PCM_Residual = sqrt(pow(ECM_Residual,2.) - pow(Masses[5],2.));
    
    if(VerboseFlag)
    {
        std::cout << "ECM_Decay: " << ECM_Decay << std::endl;
        std::cout << "ECM_Residual: " << ECM_Residual << std::endl;
        std::cout << "PCM_Decay: " << PCM_Decay << std::endl;
        std::cout << "PCM_Residual: " << PCM_Residual << std::endl;
    }
    
    TVector3 RF3MomentumDecay = TVector3(PCM_Decay*sin(DecayTheta)*cos(DecayPhi),
                                         PCM_Decay*sin(DecayTheta)*sin(DecayPhi),
                                         PCM_Decay*cos(DecayTheta));
    
    TLorentzVector RF4MomentumDecay = TLorentzVector(RF3MomentumDecay,ECM_Decay);
    
    
    
    TLorentzVector RF4MomentumResidual = TLorentzVector(-RF3MomentumDecay,ECM_Residual);
    
    //Now need to boost the 4-momenta of the decay and the residual particles from the recoil's rest frame into the lab frame
    TLorentzVector Lab4MomentumDecay = RF4MomentumDecay;
    TLorentzVector Lab4MomentumResidual = RF4MomentumResidual;
    
    if(VerboseFlag)Lab4MomentumDecay.Print();
    if(VerboseFlag)Lab4MomentumResidual.Print();
    if(VerboseFlag)Recoil4Vector.Print();
    if(VerboseFlag)Recoil4Vector.Vect().Print();
    if(VerboseFlag)std::cout << "Recoil4Vector.Beta(): " << Recoil4Vector.Beta() << std::endl;
    if(VerboseFlag)Recoil4Vector.BoostVector().Print();
    
    TVector3 Recoil3Momentum = Recoil4Vector.Vect();
    
    Lab4MomentumDecay.Boost(Recoil4Vector.BoostVector());
    Lab4MomentumResidual.Boost(Recoil4Vector.BoostVector());
    
    result[0] = Lab4MomentumDecay;
    result[1] = Lab4MomentumResidual;
    
    if(VerboseFlag)Recoil4Vector.Print();
    if(VerboseFlag)Lab4MomentumDecay.Print();
    if(VerboseFlag)Lab4MomentumResidual.Print();
    
    return result;
}

bool TestKinematics(double *Masses, double TBeam, double Ex, TLorentzVector* KinematicVectors, TLorentzVector* KinematicVectorsDecay)
{
    bool result = false;
    
    double EBeam = TBeam + Masses[0];
    double PBeam = sqrt(pow(EBeam,2.) - pow(Masses[0],2.));
    
    TLorentzVector Initial4Momemtum = TLorentzVector(0,0,PBeam,EBeam);
    
    TLorentzVector Balance4Momentum = Initial4Momemtum - KinematicVectors[0] - KinematicVectorsDecay[0] - KinematicVectorsDecay[1];
    
    if(Balance4Momentum.Mag()<1.e-6)result = true;
    else
    {
        std::cout << "Final Momentum doesn't balance!" << std::endl;
        Balance4Momentum.Print();
    }
    
    return result;
}