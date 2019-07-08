//Simple code to demonastrate how to generate a simple angular correlation output
{
    TFile *fin = TFile::Open("DipoleOutput.root");
    TTree *oak = (TTree*)fin->Get("AngCorData");
 
    TCanvas *c1 = new TCanvas("c1","c1",800,800);
    
    oak->Draw("ThetaDecayCM>>hW(181,0,180)","Weight*(ThetaAlphaLab<2.)","");
}