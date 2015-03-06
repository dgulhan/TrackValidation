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
 TString var_gen[]={"pPt","pEta","pPhi","hiNtracksEtaPtCut"};
 TString var_reco[]={"trkPt","trkEta","trkPhi","hiNtracksEtaPtCut"};
 TString axis[]={"p_{T}","#eta","#phi","multiplicity (p_{T}>0.5, |#eta|<2.4)"};
 TH1D::SetDefaultSumw2();
 
 double eta_cut=2.4;
 double pt_cut=10;
 double pt_cut_high=100;
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
   
 for(int i=0;i<4;i++){
  if(i==0){
   double ptmin = pt_cut;
   double ptmax = pt_cut_high;
 
   double inix=log(ptmin)/log(10);
   double step = (log(ptmax)-log(ptmin))/(nx*log(10));

   for(int ix=0; ix<nx+1;ix++){
    x[ix]=pow(10,inix);
    inix+=step; 
   } 
  }else if(i==1){
   double etamin =-2.4;
   double etamax = 2.4;
 
   double step = (etamax-etamin)/nx;
   x[0]=etamin; 
   for(int ix=1; ix<nx+1;ix++){
    x[ix]=x[ix-1]+step; 
   } 
  }else if(i==2){
   double phimin = -TMath::Pi();
   double phimax = TMath::Pi();
 
   double step = (phimax-phimin)/nx;
   x[0]=phimin;
   for(int ix=1; ix<nx+1;ix++){
    x[ix]=x[ix-1]+step; 
   } 
  }else if(i==3){
   double ntrkmin = 0;
   double ntrkmax = 1000;
 
   double step = (ntrkmax-ntrkmin)/nx;
   x[0]=ntrkmin;
   for(int ix=1; ix<nx+1;ix++){
    x[ix]=x[ix-1]+step; 
   } 
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
  if(i==0) c->SetLogx();
  TH1D* h = new TH1D("h","",nx,x);
  TH1D* h2 = new TH1D("h2","",nx,x);
  TH1D* h_denom = new TH1D("h_denom","",nx,x);
  TH1D* h2_denom = new TH1D("h2_denom","",nx,x);

  h2->SetLineColor(2);
  h2->SetMarkerColor(2);
  h->SetXTitle(Form("%s",axis[i].Data()));
  h->SetYTitle("Reconstruction Efficiency");
  ttrack53x->Draw(Form("%s>>h",var_gen[i].Data()),cut_eff&&cut_eff_highpurity,"",5000);
  ttrack74x->Draw(Form("%s>>h2",var_gen[i].Data()),cut_eff&&cut_eff_highpurity,"");
  ttrack53x->Draw(Form("%s>>h_denom",var_gen[i].Data()),Form("abs(pEta)<%.1f && pPt>%.1f&& pPt<%.1f",eta_cut,pt_cut,pt_cut_high),"",5000);
  ttrack74x->Draw(Form("%s>>h2_denom",var_gen[i].Data()),Form("abs(pEta)<%.1f && pPt>%.1f&& pPt<%.1f",eta_cut,pt_cut,pt_cut_high),"");
  h->Divide(h_denom);
  h2->Divide(h2_denom);
  h->Draw();
  h2->Draw("same");
  h->SetAxisRange(0,1.2,"Y");


  TLegend *leg = new TLegend(0.2,0.75,0.6,0.9);
  TLegend *leg2 = new TLegend(0.17,0.9,0.54,0.95);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry("h","53x 2.76 TeV","pl");
  leg->AddEntry("h2","74x 5.02 TeV","pl");
  leg->Draw("same");
  if(do_highpurity)drawText("high purity",0.55,0.9);
  else drawText("no quality cut",0.55,0.9);
  drawText("PYTHIA+HYDJET",0.6,0.8);  
  drawText(Form("%.1f<p_{T}<%.1f",pt_cut,pt_cut_high),0.25,0.7);

  c->SaveAs(Form("plots/eff53x74x_%s_hp%d_ptcut%.1f_%.1f.png",var_gen[i].Data(),do_highpurity,pt_cut,pt_cut_high));
  c->SaveAs(Form("plots/eff53x74x_%s_hp%d_ptcut%.1f_%.1f.pdf",var_gen[i].Data(),do_highpurity,pt_cut,pt_cut_high));
  c->SaveAs(Form("plots/eff53x74x_%s_hp%d_ptcut%.1f_%.1f.C",var_gen[i].Data(),do_highpurity,pt_cut,pt_cut_high));

  TCanvas *cFake = new TCanvas("cFake","",600,600);
  if(i==0)cFake->SetLogx();
  TH1D* hF = new TH1D("hF","",nx,x);
  TH1D* hF2 = new TH1D("hF2","",nx,x);
  TH1D* hF_denom = new TH1D("hF_denom","",nx,x);
  TH1D* hF2_denom = new TH1D("hF2_denom","",nx,x);
  hF2->SetLineColor(2);
  hF2->SetMarkerColor(2);
  hF->SetXTitle(Form("%s",axis[i].Data()));
  hF->SetYTitle("Fake Rate");
  h->SetXTitle(Form("%s",axis[i].Data()));
  h->SetYTitle("Reconstruction Efficiency");
 
  ttrack53x->Draw(Form("%s>>hF",var_reco[i].Data()),cut_fake&&cut_fake_highpurity,"",5000);
  ttrack74x->Draw(Form("%s>>hF2",var_reco[i].Data()),cut_fake&&cut_fake_highpurity,"");
  ttrack53x->Draw(Form("%s>>hF_denom",var_reco[i].Data()),Form("abs(trkEta)<%.1f && trkPt>%.1f&& trkPt<%.1f",eta_cut,pt_cut,pt_cut_high),"",5000);
  ttrack74x->Draw(Form("%s>>hF2_denom",var_reco[i].Data()),Form("abs(trkEta)<%.1f && trkPt>%.1f&& trkPt<%.1f",eta_cut,pt_cut,pt_cut_high),"");
  hF->Divide(hF_denom);
  hF2->Divide(hF2_denom);
  hF->SetAxisRange(0,1.2,"Y");
  hF->Draw();
  hF2->Draw("same");
  leg->Draw("same");

  if(do_highpurity)drawText("high purity",0.55,0.9);
  else drawText("no quality cut",0.55,0.9); 
  drawText("PYTHIA+HYDJET",0.6,0.8);
  drawText(Form("%.1f<p_{T}<%.1f",pt_cut,pt_cut_high),0.25,0.7);

  cFake->SaveAs(Form("plots/fake53x74x_%s_hp%d_ptcut%.1f_%.1f.png",var_reco[i].Data(),do_highpurity,pt_cut,pt_cut_high));
  cFake->SaveAs(Form("plots/fake53x74x_%s_hp%d_ptcut%.1f_%.1f.pdf",var_reco[i].Data(),do_highpurity,pt_cut,pt_cut_high));
  cFake->SaveAs(Form("plots/fake53x74x_%s_hp%d_ptcut%.1f_%.1f.C",var_reco[i].Data(),do_highpurity,pt_cut,pt_cut_high));
  
  if(i==0){
   TFile * outf = new TFile("outf.root","recreate");
   hF_denom->Write();
   hF2_denom->Write();
   h_denom->Write();
   h2_denom->Write();
   outf->Close();
  }
 }
 
}
