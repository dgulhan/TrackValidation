#include <iostream>
#include "TCanvas.h"
#include "TError.h"
#include "TPad.h"
#include "TString.h"
#include "TRandom.h"
#include "TH1F.h"
#include "TF1.h"
 
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TProfile.h"
#include "TLatex.h"
#include "TString.h"  
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TCut.h"
#include "TNtuple.h"
#include "TProfile2D.h"
#include "TMath.h"
#include "TLine.h"

void track_ES_ER(){
TH1D::SetDefaultSumw2();
double ptmin = 0.4;
double ptmax = 100;

TFile * f = new TFile("/mnt/hadoop/cms/store/user/yenjie/HiForest_v27/Dijet80_HydjetDrum_v27_mergedV1.root");
 
TTree * t= (TTree*)f->Get("anaTrack/trackTree");

const int ny=50;
double x[ny+1];
double inix=log(ptmin)/log(10);
double step = (log(ptmax)-log(ptmin))/(ny*log(10));

TH1D* hTER[ny];
TH1D* hTER_qual[ny];

double RMS[ny];
double RMSerr[ny];
double RMS_qual[ny];
double RMSerr_qual[ny];
double xerr[ny];
for(int ix=0; ix<ny+1;ix++){ 
 x[ix]=pow(10,inix);
 
 hTER[ix]=new TH1D(Form("hTER_%d",ix),"",50,0,2);
 hTER_qual[ix]=new TH1D(Form("hTER_qual_%d",ix),"",50,0,2);
 t->Draw(Form("mtrkPt/pPt>>hTER_%d",ix),Form("mtrkPt>0 && pPt>=%.2f && pPt<%.2f",x[ix],pow(10,inix+step)));
 t->Draw(Form("mtrkPt/pPt>>hTER_qual_%d",ix),Form("mtrkPt>0 && mtrkQual && abs(pEta)<2.4 && (mtrkDxy1/mtrkDxyError1)<3.0 && (mtrkDz1/mtrkDzError1)<3 && (mtrkPtError/mtrkPt)<0.1 && pPt>=%.2f && pPt<%.2f",x[ix],pow(10,inix+step)));
 RMS[ix]=hTER[ix]->GetRMS();
 RMSerr[ix]=hTER[ix]->GetRMSError();
 RMS_qual[ix]=hTER_qual[ix]->GetRMS();
 RMSerr_qual[ix]=hTER_qual[ix]->GetRMSError();
 xerr[ix]=0;
 inix+=step;

}  

TGraphErrors *gTER= new TGraphErrors(ny,x,RMS,xerr,RMSerr);
TGraphErrors *gTER_qual= new TGraphErrors(ny,x,RMS_qual,xerr,RMSerr_qual);


TH2D * h_ES_pt = new TH2D("h_ES_pt",";p_{T}^{gen};p_{T}^{trk}/p_{T}^{gen}",ny,x,50,0.7,1.3);
t->Draw("mtrkPt/pPt:pPt>>h_ES_pt");

TH2D * h_ES_pt_qual = new TH2D("h_ES_pt_qual",";p_{T}^{gen};p_{T}^{trk}/p_{T}^{gen}",ny,x,50,0.7,1.3);
t->Draw("mtrkPt/pPt:pPt>>h_ES_pt_qual","mtrkQual && abs(pEta)<2.4 && (mtrkDxy1/mtrkDxyError1)<3.0 && (mtrkDz1/mtrkDzError1)<3 && (mtrkPtError/mtrkPt)<0.1");

TProfile *p_ES_pt = (TProfile*)h_ES_pt->ProfileX("p_ES_pt");
TProfile *p_ES_pt_qual = (TProfile*)h_ES_pt_qual->ProfileX("p_ES_pt_qual");

TCanvas * c1 = new TCanvas("c1","",600,600);
c1->SetLogx();
c1->SetRightMargin(0.1);
h_ES_pt->Draw("colz");
p_ES_pt->Draw("same");
c1->SaveAs("plots/TES.png");

TCanvas * c2 = new TCanvas("c2","",600,600);
c2->SetLogx();
c2->SetRightMargin(0.1);
h_ES_pt_qual->Draw("colz"); 
p_ES_pt_qual->Draw("same");
c2->SaveAs("plots/TES_qual.png");

TCanvas *c3 = new TCanvas("c3","",600,600);
c3->SetLogx();
gTER->SetMaximum(0.05);
gTER->SetMinimum(0);
gTER->SetTitle("");
gTER->GetYaxis()->SetTitle("sigma(p_{T}^{trk}/p_{T}^{gen})");
gTER->GetXaxis()->SetTitle("p_{T}^{gen}");

gTER->Draw("Ap");
c3->SaveAs("plots/TER.png");

TCanvas *c4 = new TCanvas("c4","",600,600);
c4->SetLogx();

gTER_qual->SetMaximum(0.05);
gTER_qual->SetMinimum(0);
gTER_qual->SetTitle("");
gTER_qual->GetYaxis()->SetTitle("sigma(p_{T}^{trk}/p_{T}^{gen})");
gTER_qual->GetXaxis()->SetTitle("p_{T}^{gen}");

gTER_qual->Draw("Ap");

c4->SaveAs("plots/TER_qual.png");

}
