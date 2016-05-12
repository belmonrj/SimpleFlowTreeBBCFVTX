#include <TStyle.h>
#include <TFile.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TLine.h>
#include <TLatex.h>
#include <TMath.h>

#include <fstream>
#include <vector>
#include <iostream>


#include "RpPar.h"

using namespace std;

float sqr(float x) {return x*x;}
double sqr(double x){return x*x;}

void make_plots()
{
  float v2_ep_resolution_mpc = 0.163;
  float v3_ep_resolution_mpc = 0.032;
  
  float v2_ep_resolution_fvtx = 0.274;
  float v3_ep_resolution_fvtx = 0.084;

  gStyle->SetOptStat(0);

  TFile *file = TFile::Open("processed/total_alltrig_allzvtx.root");
  if(!file)
    {
      cout<<"couldn't find file, breaking"<<endl;
      return;
    }
  bool draw_sine_points = true;

  TH1D *tot_zvtx =  (TH1D*)file->Get("tot_zvtx");
  TH1D *forward_zvtx =  (TH1D*)file->Get("forward_zvtx");
  TH1D *backward_zvtx =  (TH1D*)file->Get("backward_zvtx");
  
  TProfile *q_vec_2_fvtx_mb = (TProfile*)file->Get("q_vec_2_fvtx_mb");
  TProfile *q_vec_2_fvtx_ct = (TProfile*)file->Get("q_vec_2_fvtx_ct");
  TProfile *q_vec_3_fvtx_mb = (TProfile*)file->Get("q_vec_3_fvtx_mb");
  TProfile *q_vec_3_fvtx_ct = (TProfile*)file->Get("q_vec_3_fvtx_ct");

  TProfile *q_vec_2_mpc_mb = (TProfile*)file->Get("q_vec_2_mpc_mb");
  TProfile *q_vec_2_mpc_ct = (TProfile*)file->Get("q_vec_2_mpc_ct");
  TProfile *q_vec_3_mpc_mb = (TProfile*)file->Get("q_vec_3_mpc_mb");
  TProfile *q_vec_3_mpc_ct = (TProfile*)file->Get("q_vec_3_mpc_ct");

  TProfile *phi_v2_fvtx = (TProfile*)file->Get("phi_v2_fvtx");
  TProfile *phi_v3_fvtx = (TProfile*)file->Get("phi_v3_fvtx");
  TProfile *phi_v2_mpc = (TProfile*)file->Get("phi_v2_mpc");
  TProfile *phi_v3_mpc = (TProfile*)file->Get("phi_v3_mpc");

  TProfile *fvtx_v2 = (TProfile*)file->Get("hadronfvtxv2");
  TProfile *fvtx_v3 = (TProfile*)file->Get("hadronfvtxv3");

  TProfile *mpc_v2 = (TProfile*)file->Get("hadronmpcv2");
  TProfile *mpc_v3 = (TProfile*)file->Get("hadronmpcv3");

  TProfile *sfvtx_v2 = (TProfile*)file->Get("shadronfvtxv2");
  TProfile *sfvtx_v3 = (TProfile*)file->Get("shadronfvtxv3");

  TProfile *smpc_v2 = (TProfile*)file->Get("shadronmpcv2");
  TProfile *smpc_v3 = (TProfile*)file->Get("shadronmpcv3");

  fvtx_v2->Scale(1.0/v2_ep_resolution_fvtx);
  fvtx_v3->Scale(1.0/v3_ep_resolution_fvtx);
  mpc_v2->Scale(1.0/v2_ep_resolution_mpc);
  mpc_v3->Scale(1.0/v3_ep_resolution_mpc);

  sfvtx_v2->Scale(1.0/v2_ep_resolution_fvtx);
  sfvtx_v3->Scale(1.0/v3_ep_resolution_fvtx);
  smpc_v2->Scale(1.0/v2_ep_resolution_mpc);
  smpc_v3->Scale(1.0/v3_ep_resolution_mpc);

  TH1D * eta_west_dist = (TH1D*)file->Get("eta_west_dist");
  TH1D * eta_east_dist = (TH1D*)file->Get("eta_east_dist");

  vector<TProfile*> fvtx_eta_vec;
  vector<TProfile*> mpc_eta_vec;

  for(int i = 0; i < 8; i++)
    {
      fvtx_eta_vec.push_back((TProfile*)file->Get(Form("hadronfvtxv2_eta%d",i)));
      mpc_eta_vec.push_back((TProfile*)file->Get(Form("hadronmpcv2_eta%d",i)));
      fvtx_eta_vec[i]->Scale(1.0/v2_ep_resolution_fvtx);
      mpc_eta_vec[i]->Scale(1.0/v2_ep_resolution_mpc);
    }

  vector<TProfile*> sfvtx_eta_vec;
  vector<TProfile*> smpc_eta_vec;

  for(int i = 0; i < 8; i++)
    {
      sfvtx_eta_vec.push_back((TProfile*)file->Get(Form("shadronfvtxv2_eta%d",i)));
      smpc_eta_vec.push_back((TProfile*)file->Get(Form("shadronmpcv2_eta%d",i)));
      sfvtx_eta_vec[i]->Scale(1.0/v2_ep_resolution_fvtx);
      smpc_eta_vec[i]->Scale(1.0/v2_ep_resolution_mpc);
    }

  vector<TProfile*> fvtx_eta_vec_v3;
  vector<TProfile*> mpc_eta_vec_v3;

  for(int i = 0; i < 8; i++)
    {
      fvtx_eta_vec_v3.push_back((TProfile*)file->Get(Form("hadronfvtxv3_eta%d",i)));
      mpc_eta_vec_v3.push_back((TProfile*)file->Get(Form("hadronmpcv3_eta%d",i)));
      fvtx_eta_vec_v3[i]->Scale(1.0/v3_ep_resolution_fvtx);
      mpc_eta_vec_v3[i]->Scale(1.0/v3_ep_resolution_mpc);
    }

  vector<TProfile*> sfvtx_eta_vec_v3;
  vector<TProfile*> smpc_eta_vec_v3;

  for(int i = 0; i < 8; i++)
    {
      sfvtx_eta_vec_v3.push_back((TProfile*)file->Get(Form("shadronfvtxv3_eta%d",i)));
      smpc_eta_vec_v3.push_back((TProfile*)file->Get(Form("shadronmpcv3_eta%d",i)));
      sfvtx_eta_vec_v3[i]->Scale(1.0/v3_ep_resolution_fvtx);
      smpc_eta_vec_v3[i]->Scale(1.0/v3_ep_resolution_mpc);
    }

  vector<TH1D*> q_vec_fvtx_v2_bf;
  vector<TH1D*> q_vec_fvtx_v2_af;
    
  vector<TH1D*> q_vec_mpc_v2_bf;
  vector<TH1D*> q_vec_mpc_v2_af;

  for(int i = 0; i < NZPS; i++)
    {
      q_vec_fvtx_v2_bf.push_back((TH1D*)file->Get(Form("q_vec_fvtx_v2_bf_%d",i)));
      q_vec_fvtx_v2_af.push_back((TH1D*)file->Get(Form("q_vec_fvtx_v2_af_%d",i)));
      q_vec_mpc_v2_bf.push_back((TH1D*)file->Get(Form("q_vec_mpc_v2_bf_%d",i)));
      q_vec_mpc_v2_af.push_back((TH1D*)file->Get(Form("q_vec_mpc_v2_af_%d",i)));

    }
  
  TProfile *seta_fvtx_v2 = (TProfile*)file->Get("eta_shadronfvtxv2");
  TProfile *seta_fvtx_v3 = (TProfile*)file->Get("eta_shadronfvtxv3");
  seta_fvtx_v2->Scale(1.0/v2_ep_resolution_fvtx);
  seta_fvtx_v3->Scale(1.0/v3_ep_resolution_fvtx);

  TProfile *seta_mpc_v2 = (TProfile*)file->Get("eta_shadronmpcv2");
  TProfile *seta_mpc_v3 = (TProfile*)file->Get("eta_shadronmpcv3");
  seta_mpc_v2->Scale(1.0/v2_ep_resolution_mpc);
  seta_mpc_v3->Scale(1.0/v3_ep_resolution_mpc);

  TProfile *eta_fvtx_v2 = (TProfile*)file->Get("eta_hadronfvtxv2");
  TProfile *eta_fvtx_v3 = (TProfile*)file->Get("eta_hadronfvtxv3");
  eta_fvtx_v2->Scale(1.0/v2_ep_resolution_fvtx);
  eta_fvtx_v3->Scale(1.0/v3_ep_resolution_fvtx);

  TProfile *eta_mpc_v2 = (TProfile*)file->Get("eta_hadronmpcv2");
  TProfile *eta_mpc_v3 = (TProfile*)file->Get("eta_hadronmpcv3");
  eta_mpc_v2->Scale(1.0/v2_ep_resolution_mpc);
  eta_mpc_v3->Scale(1.0/v3_ep_resolution_mpc);


  TProfile *fvtx_east_v2 = (TProfile*)file->Get("east_arm_hadronfvtxv2");
  TProfile *fvtx_east_v3 = (TProfile*)file->Get("east_arm_hadronfvtxv3");
  fvtx_east_v2->Scale(1.0/v2_ep_resolution_fvtx);
  fvtx_east_v3->Scale(1.0/v3_ep_resolution_fvtx);

  TProfile *mpc_east_v2 = (TProfile*)file->Get("east_arm_hadronmpcv2");
  TProfile *mpc_east_v3 = (TProfile*)file->Get("east_arm_hadronmpcv3");
  mpc_east_v2->Scale(1.0/v2_ep_resolution_mpc);
  mpc_east_v3->Scale(1.0/v3_ep_resolution_mpc);

  TProfile *fvtx_west_v2 = (TProfile*)file->Get("west_arm_hadronfvtxv2");
  TProfile *fvtx_west_v3 = (TProfile*)file->Get("west_arm_hadronfvtxv3");
  fvtx_west_v2->Scale(1.0/v2_ep_resolution_fvtx);
  fvtx_west_v3->Scale(1.0/v3_ep_resolution_fvtx);

  TProfile *mpc_west_v2 = (TProfile*)file->Get("west_arm_hadronmpcv2");
  TProfile *mpc_west_v3 = (TProfile*)file->Get("west_arm_hadronmpcv3");
  mpc_west_v2->Scale(1.0/v2_ep_resolution_mpc);
  mpc_west_v3->Scale(1.0/v3_ep_resolution_mpc);


  TProfile *eta_west_fvtx_v2 = (TProfile*)file->Get("eta_west_arm_hadronfvtxv2");
  TProfile *eta_west_fvtx_v3 = (TProfile*)file->Get("eta_west_arm_hadronfvtxv3");

  TProfile *eta_west_mpc_v2 = (TProfile*)file->Get("eta_west_arm_hadronmpcv2");
  TProfile *eta_west_mpc_v3 = (TProfile*)file->Get("eta_west_arm_hadronmpcv3");

  cout<<"finished retrieving histograms from file"<<endl;
  
  TH1D *dNdeta = new TH1D("dNdeta","eta ",20,-5,5);
  dNdeta->SetBinContent(0,8.50273);
  dNdeta->SetBinContent(1,9.21619);
  dNdeta->SetBinContent(2,14.965);
  dNdeta->SetBinContent(3,21.7497);
  dNdeta->SetBinContent(4,28.8159);
  dNdeta->SetBinContent(5,35.1076);
  dNdeta->SetBinContent(6,39.4898);
  dNdeta->SetBinContent(7,41.2418);
  dNdeta->SetBinContent(8,40.3878);
  dNdeta->SetBinContent(9,37.7328);
  dNdeta->SetBinContent(10,34.8883);
  dNdeta->SetBinContent(11,33.2834);
  dNdeta->SetBinContent(12,32.7893);
  dNdeta->SetBinContent(13,32.0679);
  dNdeta->SetBinContent(14,30.0184);
  dNdeta->SetBinContent(15,26.3628);
  dNdeta->SetBinContent(16,21.466);
  dNdeta->SetBinContent(17,16.0255);
  dNdeta->SetBinContent(18,10.8419);
  dNdeta->SetBinContent(19,6.56325);
  dNdeta->SetBinContent(20,3.5275);
  dNdeta->SetBinContent(21,2.77855);
  dNdeta->SetEntries(6.46582e+08);

  float pt[24];
  float ept[24], spt[24];

  float v2mpcs[24], ev2mpcs[24], sv2mpcs[24];
  float v3mpcs[24], ev3mpcs[24], sv3mpcs[24];
  TLegend *zleg = new TLegend(0.1,0.7,0.4,0.9);
  TCanvas *zcan = new TCanvas("zcan","zcan",700,700);
  tot_zvtx->SetLineColor(1);
  tot_zvtx->Draw();
  forward_zvtx->SetLineColor(2);
  forward_zvtx->Draw("same");
  backward_zvtx->SetLineColor(4);
  backward_zvtx->Draw("same");
  zleg->AddEntry(tot_zvtx,"all events","L");
  zleg->AddEntry(forward_zvtx,"events with #eta > 1.4","L");
  zleg->AddEntry(backward_zvtx,"events with #eta < -1.4","L");
  zleg->Draw("same");
  cout<<"drawing q vec histograms"<<endl;

  TCanvas *fvtx_qvec = new TCanvas("fvtx_qvec","fvtx_qvec",1500,600);
  fvtx_qvec->Divide(5,2);
  
  
  for(int i = 0; i < NZPS; i++)
    {
      fvtx_qvec->cd(i+1);
      q_vec_fvtx_v2_bf[i]->SetLineColor(2);
      q_vec_fvtx_v2_bf[i]->SetTitle(Form("zbin %d",i));
      q_vec_fvtx_v2_af[i]->SetLineColor(4);
      q_vec_fvtx_v2_bf[i]->SetXTitle("#Psi_{2}");
      q_vec_fvtx_v2_bf[i]->Draw();
      q_vec_fvtx_v2_af[i]->Draw("same");

    }

  TCanvas *mpc_qvec = new TCanvas("mpc_qvec","mpc_qvec",1500,600);
  mpc_qvec->Divide(5,2);
  for(int i = 0; i < NZPS; i++)
    {
      mpc_qvec->cd(i+1);
      q_vec_mpc_v2_bf[i]->SetLineColor(2);
      q_vec_mpc_v2_bf[i]->SetTitle(Form("zbin %d",i));
      q_vec_mpc_v2_af[i]->SetLineColor(4);
      q_vec_mpc_v2_af[i]->SetXTitle("#Psi_{2}");
      q_vec_mpc_v2_bf[i]->Draw();
      q_vec_mpc_v2_af[i]->Draw("same");

    }

  cout<<"drawing eta dist"<<endl;
  
  TLegend *etaleg = new TLegend(0.7,0.8,0.9,0.9);

  TCanvas *c7 = new TCanvas("c7","c7",800,800);
  c7->SetLogy();
  eta_west_dist->SetLineColor(kRed);
  eta_west_dist->SetXTitle("#eta");
  eta_west_dist->Draw();
  eta_east_dist->SetLineColor(kBlue);
  eta_east_dist->Draw("same");
  etaleg->AddEntry(eta_west_dist,"west tracks","L");
  etaleg->AddEntry(eta_east_dist,"east tracks","L");
  etaleg->Draw("same");
  //dNdeta->SetMarkerStyle(20);
  //dNdeta->SetMarkerColor(kBlack);
  //dNdeta->SetMarkerSize(2);
  //dNdeta->Draw("same");

  //TLegend *leg5 = new TLegend(0.65,.75,.9,.85);
  //leg5->AddEntry(eta_west_dist,"West Arm Tracks","L");
  //leg5->AddEntry(eta_east_dist,"East Arm Tracks","L");
  //leg5->AddEntry(dNdeta,"Hydro Calculation from Bozek","L");
  //leg5->Draw("same");


  cout<<"drawing v2"<<endl;

  TCanvas *c0 = new TCanvas("c0","c0",700,700);

  fvtx_v2->Rebin(2);
  //fvtx_v3->Rebin(5);
  mpc_v2->Rebin(2);
  //mpc_v3->Rebin(5);

  TLegend *leg = new TLegend(0.13,0.68,0.72,0.87);
  mpc_v2->SetXTitle("p_{T}");
  mpc_v2->SetMarkerColor(kBlue);
  mpc_v2->SetMarkerStyle(22);
  mpc_v2->SetMarkerSize(2);
  mpc_v2->GetXaxis()->SetRangeUser(0.0,3.5);
  mpc_v2->GetYaxis()->SetRangeUser(-0.025,0.25);

  fvtx_v2->SetMarkerColor(kGreen);
  fvtx_v2->SetMarkerStyle(23);
  fvtx_v2->SetMarkerSize(2);


  leg->AddEntry(mpc_v2,"VTX Tracks w MPC EP v2","p");
  leg->AddEntry(fvtx_v2,"VTX Tracks w FVTX EP v2","p");
  // final PPG181 data!!!!

  ifstream finv2("He3Au_0_5_v2.dat");
  ifstream finv3("He3Au_0_5_v3.dat");

  int imax = 13;

  for(int i=0; i<imax; i++){
    ept[i]=0;
    spt[i]=0.05;
    finv2>>pt[i]>>v2mpcs[i]>>ev2mpcs[i]>>sv2mpcs[i];
  }
  finv2.close();

  TGraphErrors *grmpcs2 = new TGraphErrors(imax, pt, v2mpcs, ept, ev2mpcs);
  TGraphErrors *sgrmpcs2 = new TGraphErrors(imax, pt, v2mpcs, spt, sv2mpcs);

  for(int i=0; i<imax; i++){
    ept[i]=0;
    spt[i]=0.05;
    finv3>>pt[i]>>v3mpcs[i]>>ev3mpcs[i]>>sv3mpcs[i];
  }
  finv3.close();

  TGraphErrors *grmpcs3 = new TGraphErrors(imax, pt, v3mpcs, ept, ev3mpcs);
  TGraphErrors *sgrmpcs3 = new TGraphErrors(imax, pt, v3mpcs, spt, sv3mpcs);

  mpc_v2->Draw("p");

  grmpcs2->SetMarkerStyle(20);
  grmpcs2->SetMarkerSize(1.4);
  grmpcs2->SetMarkerColor(2);
  grmpcs2->Draw("Psame");

  sgrmpcs2->SetMarkerStyle(20);
  sgrmpcs2->SetMarkerSize(1.4);
  sgrmpcs2->SetMarkerColor(2);
  sgrmpcs2->SetFillColor(15);
  sgrmpcs2->Draw("PE2same");


  leg->AddEntry(grmpcs2,"Final PPG181 v_{2} w FVTX EP","p");
  leg->Draw("same");


  mpc_v2->Draw("psame");
  fvtx_v2->Draw("psame");

  if(draw_sine_points)
    {

      sfvtx_v2->SetMarkerColor(kGreen);
      sfvtx_v2->SetMarkerStyle(32);
      sfvtx_v2->SetMarkerSize(2);
      sfvtx_v2->Rebin(2);
      smpc_v2->SetMarkerColor(kBlue);
      smpc_v2->SetMarkerStyle(26);
      smpc_v2->SetMarkerSize(2);
      smpc_v2->Rebin(2);
      smpc_v2->Draw("psame");
      sfvtx_v2->Draw("psame");
      leg->AddEntry(sfvtx_v2,"Sin(2*dphi) w FVTX EP","p");
      leg->AddEntry(smpc_v2,"Sin(2*dphi) w MPC EP","p");

    }
  sgrmpcs3->SetMarkerStyle(21);
  sgrmpcs3->SetMarkerSize(1.4);
  sgrmpcs3->SetMarkerColor(kRed);
  sgrmpcs3->SetFillColor(15);
  //sgrmpcs3->Draw("PE2same");

  grmpcs3->SetMarkerStyle(21);
  grmpcs3->SetMarkerSize(1.4);
  grmpcs3->SetMarkerColor(kRed);
  //grmpcs3->Draw("Psame");

  cout<<"drawing v3"<<endl;

  TCanvas *c4 = new TCanvas("c4","c4",700,700);
  TLegend *leg3 = new TLegend(0.13,0.72,0.48,0.87);
  mpc_v3->SetXTitle("p_{T}");
  mpc_v3->SetMarkerColor(kBlue);
  mpc_v3->SetMarkerStyle(22);
  mpc_v3->SetMarkerSize(2);
  mpc_v3->Rebin(2);
  mpc_v3->GetXaxis()->SetRangeUser(0.0,3.5);
  mpc_v3->GetYaxis()->SetRangeUser(-0.05,0.15);
  mpc_v3->Draw("p");
  fvtx_v3->Rebin(2);
  fvtx_v3->SetMarkerColor(kGreen);
  fvtx_v3->SetMarkerStyle(23);
  fvtx_v3->SetMarkerSize(2);
  fvtx_v3->Draw("psame");
  sgrmpcs3->Draw("PE2same");
  grmpcs3->Draw("Psame");
  //sgrmpcs2->Draw("PE2same");
  //grmpcs2->Draw("Psame");
  if(draw_sine_points)
    {

      sfvtx_v3->SetMarkerColor(kGreen);
      sfvtx_v3->SetMarkerStyle(32);
      sfvtx_v3->SetMarkerSize(2);
      sfvtx_v3->Rebin(2);
      smpc_v3->SetMarkerColor(kBlue);
      smpc_v3->SetMarkerStyle(26);
      smpc_v3->SetMarkerSize(2);
      smpc_v3->Rebin(2);
      //smpc_v3->Draw("psame");
      sfvtx_v3->Draw("psame");
      leg3->AddEntry(sfvtx_v3,"Sin(3*dphi) w FVTX EP","p");
      //leg3->AddEntry(smpc_v3,"Sin(3*dphi) w MPC EP","p");

    }

  //leg3->AddEntry(mpc_v3,"VTX Tracks w MPC EP v3","p");
  leg3->AddEntry(fvtx_v3,"VTX Tracks w FVTX EP v3","p");
  leg3->AddEntry(grmpcs3,"Final PPG181 v3 w FVTX EP","p");
  leg3->Draw("same");

  cout<<"drawing eta v3"<<endl;

  
  float eta[14] = {-1.5,-1.25,-1,-.75,-.5,-.25,0,.25,.5,.75,1,1.25,1.5};

  TCanvas *c5 = new TCanvas("c5","c5",1400,800);
  c5->Divide(3,2);
  cout<<"drawing eta v2 bins"<<endl;

  vector<TLatex*> latex_vec;

  for(int i = 0; i < 6; i++)
    {
      c5->cd(i+1);
      //latex_vec.push_back(new TLatex());
      //latex_vec[i]->SetNDC();
      //latex_vec[i]->SetText(0.55, 0.175,Form("Number of Entries %d",fvtx_eta_vec[i]->GetEntries()));
      //latex_vec[i]->SetTextSize(0.03);
      fvtx_eta_vec[i+1]->Rebin(2);
      fvtx_eta_vec[i+1]->SetTitle(Form("#eta Range: %0.1f to  %0.1f",eta[2*i],eta[2*i+2]));
      fvtx_eta_vec[i+1]->SetXTitle("p_{T}");
      fvtx_eta_vec[i+1]->GetXaxis()->SetRangeUser(0.0,3.5);
      fvtx_eta_vec[i+1]->GetYaxis()->SetRangeUser(-0.05,0.25);
      fvtx_eta_vec[i+1]->SetMarkerStyle(23);
      fvtx_eta_vec[i+1]->SetMarkerSize(2);
      fvtx_eta_vec[i+1]->SetMarkerColor(kGreen);
      fvtx_eta_vec[i+1]->Draw("p");
      mpc_eta_vec[i+1]->Rebin(2);
      mpc_eta_vec[i+1]->SetMarkerStyle(22);
      mpc_eta_vec[i+1]->SetMarkerSize(2);
      mpc_eta_vec[i+1]->SetMarkerColor(kBlue);
      mpc_eta_vec[i+1]->Draw("psame");
      sgrmpcs2->Draw("PE2same");
      grmpcs2->Draw("Psame");
      
      if(draw_sine_points)
	{
	  sfvtx_eta_vec[i+1]->SetMarkerColor(kGreen);
	  sfvtx_eta_vec[i+1]->SetMarkerStyle(32);
	  sfvtx_eta_vec[i+1]->SetMarkerSize(2);
	  sfvtx_eta_vec[i+1]->Rebin(2);
	  smpc_eta_vec[i+1]->Rebin(2);
	  smpc_eta_vec[i+1]->SetMarkerColor(kBlue);
	  smpc_eta_vec[i+1]->SetMarkerStyle(26);
	  smpc_eta_vec[i+1]->SetMarkerSize(2);
	  
	  smpc_eta_vec[i+1]->Draw("psame");
	  sfvtx_eta_vec[i+1]->Draw("psame");

	}



      if(i==1) leg->Draw("same");
      //latex_vec[i]->Draw("same");
    }

  cout<<"drawing eta v3 bins"<<endl;

  TCanvas *v3_eta_bin = new TCanvas("v3_eta_bin","v3_eta_bin",1400,800);
  v3_eta_bin->Divide(3,2);
  
  for(int i = 0; i < 6; i++)
    {
      v3_eta_bin->cd(i+1);
      //latex_vec.push_back(new TLatex());
      //latex_vec[i]->SetNDC();
      //latex_vec[i]->SetText(0.55, 0.175,Form("Number of Entries %d",fvtx_eta_vec[i]->GetEntries()));
      //latex_vec[i]->SetTextSize(0.03);
      fvtx_eta_vec_v3[i+1]->Rebin(2);
      fvtx_eta_vec_v3[i+1]->SetTitle(Form("#eta Range: %0.1f to  %0.1f",eta[2*i],eta[2*i+2]));
      fvtx_eta_vec_v3[i+1]->SetXTitle("p_{T}");
      fvtx_eta_vec_v3[i+1]->GetXaxis()->SetRangeUser(0.0,3.5);
      fvtx_eta_vec_v3[i+1]->GetYaxis()->SetRangeUser(-0.05,0.1);
      fvtx_eta_vec_v3[i+1]->SetMarkerStyle(23);
      fvtx_eta_vec_v3[i+1]->SetMarkerSize(2);
      fvtx_eta_vec_v3[i+1]->SetMarkerColor(kGreen);
      fvtx_eta_vec_v3[i+1]->Draw("p");
      mpc_eta_vec_v3[i+1]->Rebin(2);
      mpc_eta_vec_v3[i+1]->SetMarkerStyle(22);
      mpc_eta_vec_v3[i+1]->SetMarkerSize(2);
      mpc_eta_vec_v3[i+1]->SetMarkerColor(kBlue);
      //mpc_eta_vec_v3[i+1]->Draw("psame");
      sgrmpcs3->Draw("PE2same");
      grmpcs3->Draw("Psame");

      if(draw_sine_points )
	{
	sfvtx_eta_vec_v3[i+1]->SetMarkerColor(kGreen);
	sfvtx_eta_vec_v3[i+1]->SetMarkerStyle(32);
	sfvtx_eta_vec_v3[i+1]->SetMarkerSize(2);
	sfvtx_eta_vec_v3[i+1]->Rebin(2);
	smpc_eta_vec_v3[i+1]->Rebin(2);
	smpc_eta_vec_v3[i+1]->SetMarkerColor(kBlue);
	smpc_eta_vec_v3[i+1]->SetMarkerStyle(26);
	smpc_eta_vec_v3[i+1]->SetMarkerSize(2);
	  
	//smpc_eta_vec_v3[i+1]->Draw("psame");
	sfvtx_eta_vec_v3[i+1]->Draw("psame");

	}
      


      if(i==1) leg3->Draw("same");
      //latex_vec[i]->Draw("same");
    }


  cout<<"drawing v2 east/west"<<endl;

  TLegend *east_west_leg = new TLegend(0.1,0.70,0.35,0.85);
  TCanvas *c2 = new TCanvas("c2","c2",700,1400);
  c2->Divide(1,2);
  c2->cd(1);
  fvtx_east_v2->Rebin(2);
  fvtx_west_v2->Rebin(2);
  fvtx_east_v2->SetMarkerStyle(26);
  fvtx_east_v2->SetMarkerColor(kBlue);
  fvtx_east_v2->SetMarkerSize(2);
  fvtx_west_v2->SetMarkerStyle(32);
  fvtx_west_v2->SetMarkerColor(kRed);
  fvtx_west_v2->SetMarkerSize(2);
  fvtx_east_v2->SetTitle("FVTX EP v_{2}");
  fvtx_east_v2->Draw("p");
  fvtx_west_v2->Draw("psame");
  fvtx_v2->Draw("psame");
  fvtx_east_v2->GetXaxis()->SetRangeUser(0.0,3.5);
  east_west_leg->AddEntry(fvtx_east_v2,"East v_{2}","p");
  east_west_leg->AddEntry(fvtx_west_v2,"West v_{2}","p");
  east_west_leg->AddEntry(fvtx_v2,"Inclusive v_{2}","p");
  east_west_leg->Draw("same");

  c2->cd(2);
  TH1D *fvtx_east_v2hist = new TH1D("fvtx_east_v2hist","fvtx_east_v2hist",35,0,7);
  for(int i = 1; i < fvtx_east_v2->GetNbinsX()+1;i++)
    {
      double cdiv = fvtx_east_v2->GetBinContent(i)/fvtx_west_v2->GetBinContent(i);
      fvtx_east_v2hist->SetBinContent(i,cdiv);
      fvtx_east_v2hist->SetBinError(i,TMath::Sqrt(cdiv*cdiv*(sqr(fvtx_east_v2->GetBinError(i)/fvtx_east_v2->GetBinContent(i))+fvtx_west_v2->GetBinError(i)/fvtx_west_v2->GetBinContent(i))));
    }
  fvtx_east_v2hist->SetMarkerStyle(20);
  fvtx_east_v2hist->SetMarkerColor(1);
  fvtx_east_v2hist->SetMarkerSize(2);
  fvtx_east_v2hist->Draw("p");
  fvtx_east_v2hist->GetXaxis()->SetRangeUser(0.0,3.5);
  TLine *line = new TLine(0.0,1.0,3.5,1.0);
  line->Draw("same");

  TCanvas *c3 = new TCanvas("c3","c3",700,1400);
  c3->Divide(1,2);
  c3->cd(1);
  mpc_east_v2->Rebin(2);
  mpc_west_v2->Rebin(2);
  mpc_east_v2->SetMarkerStyle(26);
  mpc_east_v2->SetMarkerColor(kBlue);
  mpc_east_v2->SetMarkerSize(2);
  mpc_west_v2->SetMarkerStyle(32);
  mpc_west_v2->SetMarkerColor(kRed);
  mpc_west_v2->SetMarkerSize(2);
  mpc_east_v2->Draw("p");
  mpc_west_v2->Draw("psame");
  mpc_v2->Draw("psame");
  mpc_east_v2->GetXaxis()->SetRangeUser(0.0,3.5);
  east_west_leg->Draw("same");
  c3->cd(2);
  TH1D *mpc_east_v2hist = new TH1D("mpc_east_v2hist","mpc_east_v2hist",35,0,7);
  for(int i = 1; i < mpc_east_v2->GetNbinsX()+1;i++)
    {
      double cdiv = mpc_east_v2->GetBinContent(i)/mpc_west_v2->GetBinContent(i);
      mpc_east_v2hist->SetBinContent(i,cdiv);
      mpc_east_v2hist->SetBinError(i,TMath::Sqrt(cdiv*cdiv*(sqr(mpc_east_v2->GetBinError(i)/mpc_east_v2->GetBinContent(i))+mpc_west_v2->GetBinError(i)/mpc_west_v2->GetBinContent(i))));
    }
  mpc_east_v2hist->SetMarkerStyle(20);
  mpc_east_v2hist->SetMarkerColor(1);
  mpc_east_v2hist->SetMarkerSize(2);
  mpc_east_v2hist->Draw("p");
  mpc_east_v2hist->GetXaxis()->SetRangeUser(0.0,3.5);
  line->Draw("same");

  TLegend *v3east_west_leg = new TLegend(0.1,0.70,0.35,0.85);
  TCanvas *fvtx_v3_eastwestcan = new TCanvas("fvtx_v3_eastwestcan","fvtx_v3_eastwestcan",700,1400);
  fvtx_v3_eastwestcan->Divide(1,2);
  fvtx_v3_eastwestcan->cd(1);
  fvtx_east_v3->Rebin(2);
  fvtx_west_v3->Rebin(2);
  fvtx_east_v3->SetMarkerStyle(26);
  fvtx_east_v3->SetMarkerColor(kBlue);
  fvtx_east_v3->SetMarkerSize(2);
  fvtx_west_v3->SetMarkerStyle(32);
  fvtx_west_v3->SetMarkerColor(kRed);
  fvtx_west_v3->SetMarkerSize(2);
  fvtx_east_v3->Draw("p");
  fvtx_west_v3->Draw("psame");
  fvtx_v3->Draw("psame");
  fvtx_east_v3->GetXaxis()->SetRangeUser(0.0,3.5);
  v3east_west_leg->AddEntry(fvtx_east_v3,"East v_{2}","p");
  v3east_west_leg->AddEntry(fvtx_west_v3,"West v_{2}","p");
  v3east_west_leg->AddEntry(fvtx_v3,"Inclusive v_{2}","p");
  v3east_west_leg->Draw("same");
  fvtx_v3_eastwestcan->cd(2);
  TH1D *fvtx_east_v3hist = new TH1D("fvtx_east_v3hist","fvtx_east_v3hist",35,0,7);
  for(int i = 1; i < fvtx_east_v3->GetNbinsX()+1;i++)
    {
      double cdiv = fvtx_east_v3->GetBinContent(i)/fvtx_west_v3->GetBinContent(i);
      fvtx_east_v3hist->SetBinContent(i,cdiv);
      fvtx_east_v3hist->SetBinError(i,TMath::Sqrt(cdiv*cdiv*(sqr(fvtx_east_v3->GetBinError(i)/fvtx_east_v3->GetBinContent(i))+fvtx_west_v3->GetBinError(i)/fvtx_west_v3->GetBinContent(i))));
    }
  fvtx_east_v3hist->SetMarkerStyle(20);
  fvtx_east_v3hist->SetMarkerColor(1);
  fvtx_east_v3hist->SetMarkerSize(2);
  fvtx_east_v3hist->Draw("p");
  fvtx_east_v3hist->GetXaxis()->SetRangeUser(0.0,3.5);
  line->Draw("same");

  float bins[4] = {-2.0,-0.35,0.35,2.0};
  TProfile *tmp_ppg181_v2 = new TProfile("tmp_ppg181","tmp_ppg181",3,bins);

  float ppg181_entries[13] = {2.72405e+06,1.77226e+06,1.11693e+06,694079,432750,268815,170151,107307,69351,44607,29454,19749,13405};
  float sumppg181entries = 0;
  double favgv2 = 0;
  for(int i = 0; i < imax; i++) sumppg181entries+=ppg181_entries[i];
  for(int i = 0; i <imax;i++)
    {
      tmp_ppg181_v2->Fill(0.0,(double)v2mpcs[i],(double)ppg181_entries[i]/sumppg181entries);
    }
  
  cout<<"drawing eta v2 "<<endl;

  TCanvas *c1 = new TCanvas("c1","c1",700,700);
  TLegend *leg2 = new TLegend(0.55,0.7,0.9,0.85);
  tmp_ppg181_v2->SetMarkerStyle(20);
  tmp_ppg181_v2->SetMarkerColor(kRed);
  tmp_ppg181_v2->SetMarkerSize(2);
  eta_fvtx_v2->SetMarkerStyle(22);
  eta_fvtx_v2->SetMarkerColor(kGreen);
  eta_fvtx_v2->SetMarkerSize(2);
  eta_fvtx_v2->SetXTitle("#eta");
  eta_fvtx_v2->GetYaxis()->SetRangeUser(-0.05,0.23);
  eta_fvtx_v2->Draw("p");
  eta_mpc_v2->SetMarkerStyle(23);
  eta_mpc_v2->SetMarkerColor(kBlue);
  eta_mpc_v2->SetMarkerSize(2);
  eta_mpc_v2->Draw("psame");
  tmp_ppg181_v2->Draw("psame");
  leg2->AddEntry(eta_fvtx_v2,"VTX Tracks v_{2} w FVTX EP","p");
  leg2->AddEntry(eta_mpc_v2,"VTX Tracks v_{2} w MPC EP","p");
  leg2->AddEntry(tmp_ppg181_v2,"Final PPG 181 v_{2}","p");

  if(draw_sine_points)
    {
      seta_fvtx_v2->SetMarkerColor(kGreen);
      seta_fvtx_v2->SetMarkerStyle(32);
      seta_fvtx_v2->SetMarkerSize(2);

      seta_mpc_v2->SetMarkerColor(kBlue);
      seta_mpc_v2->SetMarkerStyle(26);
      seta_mpc_v2->SetMarkerSize(2);

      seta_mpc_v2->Draw("psame");
      seta_fvtx_v2->Draw("psame");

      leg2->AddEntry(seta_fvtx_v2,"Sin(2*dphi) w FVTX EP","p");
      leg2->AddEntry(seta_mpc_v2,"Sin(2*dphi) w MPC EP","p");

    }

  
  leg2->Draw("same");

  TProfile *tmp_ppg181_v3 = new TProfile("tmp_ppg181_v3","tmp_ppg181_v3",3,bins);

  for(int i = 0; i <imax;i++)
    {
      tmp_ppg181_v3->Fill(0.0,(double)v3mpcs[i],(double)ppg181_entries[i]/sumppg181entries);
    }

  cout<<"drawing eta v3 "<<endl;

  TCanvas *eta_v3_can = new TCanvas("eta_v3_can","eta_v3_can",700,700);
  TLegend *leg3e = new TLegend(0.55,0.7,0.9,0.85);
  tmp_ppg181_v3->SetMarkerStyle(20);
  tmp_ppg181_v3->SetMarkerColor(kRed);
  tmp_ppg181_v3->SetMarkerSize(2);
  eta_fvtx_v3->SetMarkerStyle(22);
  eta_fvtx_v3->SetMarkerColor(kGreen);
  eta_fvtx_v3->SetMarkerSize(2);
  eta_fvtx_v3->SetXTitle("#eta");
  eta_fvtx_v3->GetYaxis()->SetRangeUser(-0.05,eta_fvtx_v3->GetMaximum()*1.1);
  eta_fvtx_v3->Draw("p");
  eta_mpc_v3->SetMarkerStyle(23);
  eta_mpc_v3->SetMarkerColor(kBlue);
  eta_mpc_v3->SetMarkerSize(2);
  eta_mpc_v3->Draw("psame");
  tmp_ppg181_v3->Draw("same");
  leg3e->AddEntry(eta_fvtx_v3,"VTX Tracks v_{3} w FVTX EP","p");
  //leg3e->AddEntry(eta_mpc_v3,"VTX Tracks v_{3} w MPC EP","p");
  leg3e->AddEntry(tmp_ppg181_v3,"PPG 181 v_{3}","p");
  
  if(draw_sine_points)
    {
      seta_fvtx_v3->SetMarkerColor(kGreen);
      seta_fvtx_v3->SetMarkerStyle(32);
      seta_fvtx_v3->SetMarkerSize(2);

      seta_mpc_v3->SetMarkerColor(kBlue);
      seta_mpc_v3->SetMarkerStyle(26);
      seta_mpc_v3->SetMarkerSize(2);

      //seta_mpc_v3->Draw("psame");
      seta_fvtx_v3->Draw("psame");

      leg3e->AddEntry(seta_fvtx_v3,"Sin(3*dphi) w FVTX EP","p");
      //leg2->AddEntry(seta_mpc_v3,"Sin(3*dphi) w MPC EP","p");

    }
  leg3e->Draw("same");

  cout<<"drawing q vec trig bias "<<endl;


  q_vec_2_fvtx_mb->Sumw2();
  q_vec_2_fvtx_ct->Sumw2();
  q_vec_3_fvtx_mb->Sumw2();
  q_vec_3_fvtx_ct->Sumw2();

  q_vec_2_fvtx_mb->Scale(1.0/q_vec_2_fvtx_mb->Integral());
  q_vec_2_fvtx_ct->Scale(1.0/q_vec_2_fvtx_ct->Integral());
  q_vec_3_fvtx_mb->Scale(1.0/q_vec_3_fvtx_mb->Integral());
  q_vec_3_fvtx_ct->Scale(1.0/q_vec_3_fvtx_ct->Integral());

  TLegend *bias_leg = new TLegend(0.65,0.7,0.95,0.9);
  TCanvas*q_vec_bias_fvtx = new TCanvas("q_vec_bias_fvtx","q_vec_bias_fvtx",1200,1200);
  q_vec_bias_fvtx->Divide(2,2);
  q_vec_bias_fvtx->cd(1);
  q_vec_2_fvtx_mb->SetLineColor(4);
  q_vec_2_fvtx_mb->SetXTitle("FVTX #Psi_{2}");
  q_vec_2_fvtx_mb->SetTitle("");
  q_vec_2_fvtx_mb->Draw();
  q_vec_2_fvtx_ct->SetLineColor(2);
  q_vec_2_fvtx_ct->Draw("same");
  q_vec_bias_fvtx->cd(2);
  q_vec_3_fvtx_mb->SetLineColor(4);
  q_vec_3_fvtx_mb->SetXTitle("FVTX #Psi_{3}");
  q_vec_3_fvtx_mb->SetTitle("");
  q_vec_3_fvtx_mb->Draw();
  q_vec_3_fvtx_ct->SetLineColor(2);
  q_vec_3_fvtx_ct->Draw("same");
  bias_leg->AddEntry(q_vec_2_fvtx_ct,"0-5% Central Trigger");
  bias_leg->AddEntry(q_vec_2_fvtx_mb,"MinBias No Central Trigger");
  bias_leg->Draw("same");
  
  q_vec_bias_fvtx->cd(3);
  TH1D* tmp_fvtx_2 = (TH1D*)q_vec_2_fvtx_mb->Clone();
  tmp_fvtx_2->Divide(q_vec_2_fvtx_ct);
  tmp_fvtx_2->GetYaxis()->SetRangeUser(0.97,1.03);
  tmp_fvtx_2->SetTitle(" ");
  tmp_fvtx_2->SetTitle("Minbias (no cent trig)/ Central Trig");
  tmp_fvtx_2->SetXTitle("FVTX #Psi_{2}");
  tmp_fvtx_2->Draw();
  q_vec_bias_fvtx->cd(4);
  TH1D* tmp_fvtx_3 = (TH1D*)q_vec_3_fvtx_mb->Clone();
  tmp_fvtx_3->Divide(q_vec_3_fvtx_ct);
  tmp_fvtx_3->GetYaxis()->SetRangeUser(0.97,1.03);
  tmp_fvtx_3->SetTitle("Minbias (no cent trig)/ Central Trig");
  tmp_fvtx_3->SetXTitle("FVTX #Psi_{3}");
  tmp_fvtx_3->Draw();

  q_vec_2_mpc_mb->Sumw2();
  q_vec_2_mpc_ct->Sumw2();
  q_vec_3_mpc_mb->Sumw2();
  q_vec_3_mpc_ct->Sumw2();

  q_vec_2_mpc_mb->Scale(1.0/q_vec_2_mpc_mb->Integral());
  q_vec_2_mpc_ct->Scale(1.0/q_vec_2_mpc_ct->Integral());
  q_vec_3_mpc_mb->Scale(1.0/q_vec_3_mpc_mb->Integral());
  q_vec_3_mpc_ct->Scale(1.0/q_vec_3_mpc_ct->Integral());

  TCanvas*q_vec_bias_mpc = new TCanvas("q_vec_bias_mpc","q_vec_bias_mpc",1200,1200);
  q_vec_bias_mpc->Divide(2,2);
  q_vec_bias_mpc->cd(1);
  q_vec_2_mpc_mb->SetLineColor(4);
  q_vec_2_mpc_mb->SetXTitle("MPC #Psi_{2}");
  q_vec_2_mpc_mb->SetTitle("");
  q_vec_2_mpc_mb->Draw();
  q_vec_2_mpc_ct->SetLineColor(2);
  q_vec_2_mpc_ct->Draw("same");
  q_vec_bias_mpc->cd(2);
  q_vec_3_mpc_mb->SetLineColor(4);
  q_vec_3_mpc_mb->SetXTitle("MPC #Psi_{3}");
  q_vec_3_mpc_mb->SetTitle("");
  q_vec_3_mpc_mb->Draw();
  q_vec_3_mpc_ct->SetLineColor(2);
  q_vec_3_mpc_ct->Draw("same");
  bias_leg->Draw("same");

  q_vec_bias_mpc->cd(3);
  TH1D* tmp_mpc_2 = (TH1D*)q_vec_2_mpc_mb->Clone();
  tmp_mpc_2->Divide(q_vec_2_mpc_ct);
  tmp_mpc_2->GetYaxis()->SetRangeUser(0.97,1.03);
  tmp_mpc_2->SetTitle("Minbias (no cent trig)/ Central Trig");
  tmp_mpc_2->SetXTitle("MPC #Psi_{2}");
  tmp_mpc_2->Draw();
  q_vec_bias_mpc->cd(4);
  TH1D* tmp_mpc_3 = (TH1D*)q_vec_3_mpc_mb->Clone();
  tmp_mpc_3->Divide(q_vec_3_mpc_ct);
  //tmp_mpc_3->GetYaxis()->SetRangeUser(0.97,1.03);
  tmp_fvtx_3->SetTitle("Minbias (no cent trig)/ Central Trig");
  tmp_fvtx_3->SetXTitle("MPC #Psi_{3}");
  tmp_mpc_3->Draw();

  TLegend *phi_leg = new TLegend(0.7,0.7,0.9,0.9);
  TCanvas *phi_vcan = new TCanvas("phi_vcan","phi_vcan",700,1400);
  phi_vcan->Divide(1,2);
  phi_vcan->cd(1);
  phi_v2_mpc->SetMarkerStyle(23);
  phi_v2_mpc->SetMarkerColor(kBlue);
  phi_v2_mpc->SetMarkerSize(2);
  phi_v2_fvtx->SetXTitle("#phi");
  phi_v2_fvtx->SetMarkerStyle(22);
  phi_v2_fvtx->SetMarkerColor(kGreen);
  phi_v2_fvtx->SetMarkerSize(2);
  phi_v2_fvtx->SetYTitle("v_{2}");
  phi_v2_fvtx->SetTitle("v_{2} #phi dependence");
  phi_v2_fvtx->GetYaxis()->SetRangeUser(-0.004,0.026);
  phi_v2_fvtx->Draw("p");
  phi_v2_mpc->Draw("psame");

  phi_vcan->cd(2);
  phi_v3_fvtx->SetMarkerStyle(22);
  phi_v3_fvtx->SetMarkerColor(kGreen);
  phi_v3_fvtx->SetMarkerSize(2);
  phi_v3_fvtx->GetYaxis()->SetRangeUser(-0.005,0.01);
  phi_v3_fvtx->SetXTitle("#phi");
  phi_v3_fvtx->SetTitle("");
  phi_v3_fvtx->SetYTitle("v_{3}");
  phi_v3_fvtx->Draw("p");
  phi_leg->AddEntry(phi_v2_mpc,"mpc","p");
  phi_leg->AddEntry(phi_v2_fvtx,"fvtx","p");
  phi_leg->Draw("same");
}
