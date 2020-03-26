//Simple code to demonastrate how to generate a simple angular correlation output - should make two PDFs with the functions drawn in them - YOU NEED TO HAVE RUN THE TESTCALCUALTION SCRIPT FIRST!
{
    TFile *fin = TFile::Open("DipoleOutput.root");
    TTree *oak = (TTree*)fin->Get("AngCorData");
 
    TCanvas *c1 = new TCanvas("c1","c1",1200,800);
    
    oak->Draw("ThetaDecayCM>>hW(181,0,180)","Weight*(ThetaAlphaLab<2.)","");
    
    TH1F *hW = (TH1F*)gROOT->FindObjectAny("hW");
    hW->GetXaxis()->SetTitle("#theta [deg]");
    hW->GetXaxis()->CenterTitle();
    hW->GetYaxis()->SetTitle("W(#theta) [a.u.]");
    hW->GetYaxis()->CenterTitle();
    hW->SetStats(0);
    hW->SetTitle("");
    hW->Draw();
    c1->Update();
    
    TFile *fin2 = TFile::Open("QuadOutput.root");
    TTree *oak2 = (TTree*)fin2->Get("AngCorData");
 
    TCanvas *c2 = new TCanvas("c2","c2",1200,800);
    
    oak2->Draw("ThetaDecayCM>>hW2(181,0,180)","Weight*(ThetaAlphaLab<2.)","");
    
    TH1F *hW2 = (TH1F*)gROOT->FindObjectAny("hW2");
    hW2->GetXaxis()->SetTitle("#theta [deg]");
    hW2->GetXaxis()->CenterTitle();
    hW2->GetYaxis()->SetTitle("W(#theta) [a.u.]");
    hW2->GetYaxis()->CenterTitle();
    hW2->SetTitle("");
    hW2->SetStats(0);
    hW2->Draw();
    c2->Update();
    
    c1->SaveAs("DipoleOutput.pdf");
    c2->SaveAs("QuadOutput.pdf");
    
    gROOT->ProcessLine(".q");
}