#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <TFile.h>
#include <TProfile.h>
#include <TMath.h>
#include <TCut.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <iostream>

// void plotEff(char *infname="/d102/dgulhan/trackingForest/CMSSW_5_3_14/src/HiForest-test.root")
// void plotEff(char *infname="/d102/dgulhan/trackAnalysis/forests/reco-100-0.4-0.4.root")
// void plotEff(char *infname="/d102/dgulhan/trackingForest/CMSSW_5_3_14/src/HiForest-test-step0_default-step1_0p4.root")

void drawText(const char *text, float xp, float yp, int size=22){
  TLatex *tex = new TLatex(xp,yp,text);
  tex->SetTextFont(43);
  tex->SetTextSize(size);
  tex->SetTextColor(kBlack);
  tex->SetLineWidth(1);
  tex->SetNDC();
  tex->Draw("same");
}

void plotEff(bool do_highpurity=true)
{
 TString var_reco[]={"trkDxy1/trkDxyError1","trkDz1/trkDzError1"};
 TString axis[]={"trkDxy1/trkDxyError1","trkDz1/trkDzError1"};
 TH1D::SetDefaultSumw2();
 
 double eta_cut=2.4;
 double pt_cut=0.5;
 double pt_cut_high=1;
 TCut cut_fake = "trkFake==1";
 TCut cut_eff = "mtrkPt>0";
 TCut cut_fake_highpurity, cut_eff_highpurity;
 
 if(do_highpurity){ 
  cut_fake_highpurity = Form("highPurity && abs(trkEta)<%.1f && trkPt>%.1f && trkPt<%.1f",eta_cut,pt_cut,pt_cut_high);
  cut_eff_highpurity = Form("mtrkQual && abs(pEta)<%.1f && pPt>%.1f && pPt<%.1f",eta_cut,pt_cut,pt_cut_high);
 }else{
  cut_fake_highpurity = Form("abs(trkEta)<%.1f && trkPt>%.1f && trkPt<%.1f",eta_cut,pt_cut,pt_cut_high);
  cut_eff_highpurity =  Form("abs(pEta)<%.1f && pPt>%.1f && pPt<%.1f",eta_cut,pt_cut,pt_cut_high);
 }
 
 const int nx=50;
 double x[nx+1];  
   
 for(int i=0;i<2;i++){
  
   double xmin =-4;
   double xmax = 4;
 
   double step = (xmax-xmin)/nx;
   x[0]=xmin; 
   for(int ix=1; ix<nx+1;ix++){
    x[ix]=x[ix-1]+step; 
   } 

 
  TFile *inf53x = new TFile("/mnt/hadoop/cms/store/user/dgulhan/PYTHIA_HYDJET_Track9_Jet30_Pyquen_DiJet_TuneZ2_Unquenched_Hydjet1p8_2760GeV_merged/HiForest_PYTHIA_HYDJET_pthat80_Track9_Jet30_matchEqR_merged_forest_0.root");
   // TFile *inf74x = new TFile("../../merge74x/HiForest_PyquenUnquenched_pthat80_74x_merged.root");
  TFile *inf74x = new TFile("/mnt/hadoop/cms/store/user/dgulhan/HiForest_Dijet_pthat80_740pre6_BS_merged/HiForest_Dijet_pthat80_740pre6_BS.root");
  TTree *ttrack53x = (TTree*)inf53x->Get("anaTrack/trackTree");
  TTree *ttrack74x = (TTree*)inf74x->Get("anaTrack/trackTree");
  TTree *tgen53x = (TTree*)inf53x->Get("HiGenParticleAna/hi");
  TTree *tgen74x = (TTree*)inf74x->Get("HiGenParticleAna/hi");
  TTree *thi53x = (TTree*)inf53x->Get("hiEvtAnalyzer/HiTree");
  TTree *thi74x = (TTree*)inf74x->Get("hiEvtAnalyzer/HiTree");
  ttrack74x->AddFriend(tgen74x);
  ttrack53x->AddFriend(tgen53x);
  ttrack74x->AddFriend(thi74x);
  ttrack53x->AddFriend(thi53x); 

  TCanvas *c = new TCanvas("c","",600,600);
  TH1D* h = new TH1D("h","",nx,x);
  TH1D* h2 = new TH1D("h2","",nx,x);
  TH1D* h_denom = new TH1D("h_denom","",nx,x);
  TH1D* h2_denom = new TH1D("h2_denom","",nx,x);

  h2->SetLineColor(2);
  h2->SetMarkerColor(2);
  h->SetXTitle(Form("%s",axis[i].Data()));
  h->SetYTitle("74x/53x");
  ttrack53x->Draw(Form("%s>>h2",var_reco[i].Data()),cut_eff&&cut_eff_highpurity,"",5000);
  ttrack74x->Draw(Form("%s>>h",var_reco[i].Data()),cut_eff&&cut_eff_highpurity,"");
  h->Scale(1/h->Integral());
  h2->Scale(1/h2->Integral());
  h->Divide(h2);
  h->Draw();
  h->SetAxisRange(0,2,"Y");


  TLegend *leg = new TLegend(0.2,0.75,0.6,0.9);
  TLegend *leg2 = new TLegend(0.17,0.9,0.54,0.95);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry("h","53x 2.76 TeV","pl");
  // leg->Draw("same");
  if(do_highpurity)drawText("high purity",0.55,0.9);
  else drawText("no quality cut",0.6,0.9);
  drawText("PYTHIA+HYDJET",0.6,0.8);
  drawText(Form("%.1f<p_{T}<%.1f",pt_cut,pt_cut_high),0.25,0.7);

  c->SaveAs(Form("plots/ratio_%d_hp%d_ptcut%.1f_%.1f.png",i,do_highpurity,pt_cut,pt_cut_high));
  c->SaveAs(Form("plots/ratio_%d_hp%d_ptcut%.1f_%.1f.pdf",i,do_highpurity,pt_cut,pt_cut_high));
  c->SaveAs(Form("plots/ratio_%d_hp%d_ptcut%.1f_%.1f.C",i,do_highpurity,pt_cut,pt_cut_high));

 
 }
}