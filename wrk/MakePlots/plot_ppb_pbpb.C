#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <numeric>
#include <iterator>
#include <cmath>
#include <TLegend.h>
#include <TLatex.h>
#include <TGraphErrors.h>

#include "TLorentzVector.h"
#include "TFile.h"
#include "TStyle.h"
#include "TNtuple.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TF1.h"

using namespace std;

void plot_ppb_pbpb()
{
  gStyle->SetOptStat(0);

  int x_max = 130;   // PbPb: maximum 6000     pPb: 300
  int bin_size = 50; // PbPb: maximum 100      pPb: 50

  // ---------------------------------------------------------------------
  //parton participant plane
  // TFile *fp_pplane = TFile::Open("pbpb_919_ppplane.root");
  // TFile *fp_pplane = TFile::Open("pbpb_ppplane_less10000.root");   // pbpb
  TFile *fp_pplane = TFile::Open("dAu_ppplane200.root");

  // TFile *fp_pplane = TFile::Open("dAu_ppplane_mixed.root");
  // TFile *fp_pplane = TFile::Open("ppb_ppplane_ep2.root");
  // TFile *fp_pplane = TFile::Open("dAu_mixed_ppplane.root");

  // TFile *fp_pplane = TFile::Open("pp_5t_ppplane.root");  // pp

  //TFile *fp_pplane = TFile::Open("dAu_D_ppplane.root"); //Darren

  TH1F* ppplane = (TH1F*)fp_pplane->Get("v2s");
  TProfile* tp1f_epsilon2 = (TProfile*)fp_pplane->Get("epsilon2_nch");

  TH1D* th1d_epsilon2 = tp1f_epsilon2->ProjectionX("th1d_epsilon2");
  for (int i=1; i<bin_size+1; i++)
    {
      double ep2_value = tp1f_epsilon2->GetBinContent(i);
      double ep2_entries = tp1f_epsilon2->GetBinEntries(i);
      double ep2_error = tp1f_epsilon2->GetBinError(i);

      if (ep2_value > 0)
	{
	  double sigep2_ep2 = ep2_error * TMath::Sqrt(ep2_entries) / ep2_value;
	  th1d_epsilon2->SetBinContent(i, sigep2_ep2);
	}
      else th1d_epsilon2->SetBinContent(i, -9999);
    }


  // ---------------------------------------------------------------------
  // TFile *file = TFile::Open("pbpb_919_cumulant.root");
  // TFile *file = TFile::Open("pbpb_cumulant_less10000.root");    // pbpb
  TFile *file = TFile::Open("dAu_cumulant200.root");

  // TFile *file = TFile::Open("dAu_cumulant_mixed.root");
  // TFile *file = TFile::Open("ppb_cumulant_sigma.root");
  // TFile *file = TFile::Open("dAu_mixed_cumulant.root");

  // TFile *file = TFile::Open("pp_5t_cumulant.root"); // pp

  //TFile *file = TFile::Open("dAu_eff.root"); //Darren

  TProfile* tp1f_raa4 = (TProfile*)file->Get("raa4_Ncharge");
  TProfile* tp1f_raa2 = (TProfile*)file->Get("raa2_Ncharge");
  TProfile* tp1f_daa4 = (TProfile*)file->Get("daa4_Ncharge");
  TProfile* tp1f_daa2 = (TProfile*)file->Get("daa2_Ncharge");
  TProfile* tp1f_gapp = (TProfile*)file->Get("comp_Ncharge");


  TH1D* th1d_raa4 = tp1f_raa4->ProjectionX("th1d_raa4");
  TH1D* th1d_raa2 = tp1f_raa2->ProjectionX("th1d_raa2");
  TH1D* th1d_daa4 = tp1f_daa4->ProjectionX("th1d_daa4");
  TH1D* th1d_daa2 = tp1f_daa2->ProjectionX("th1d_daa2");
  TH1D* th1d_gapp = tp1f_gapp->ProjectionX("th1d_gapp");

  // Calculate d_{2}{4}
  TH1D* th1d_d2 = (TH1D*)th1d_daa2->Clone("th1d_d2");
  th1d_d2->Multiply(th1d_raa2);
  th1d_d2->Scale(2.0);

  TH1D* th1d_d24 = (TH1D*)th1d_daa4->Clone("th1d_d24");
  th1d_d24->Add(th1d_d2, -1.0);

  // Calculate c_{2}{4}
  TH1D* th1d_c2 = (TH1D*)th1d_raa2->Clone("th1d_c2");
  th1d_c2->Multiply(th1d_raa2);
  th1d_c2->Scale(2.0);

  TH1D* th1d_c24 = (TH1D*)th1d_raa4->Clone("th1d_c24");
  th1d_c24->Add(th1d_c2, -1.0);

  TH1D* th1d_cv24 = (TH1D*)th1d_c24->Clone("th1d_cv24");

  for (int i = 1; i < bin_size+1; i++)
    {
      double temp_1 = th1d_cv24->GetBinContent(i);
      double erro_1 = th1d_cv24->GetBinError(i);

      if (temp_1 >= 0) th1d_cv24->SetBinContent(i, -9999);
      else
	{
	  double temp_2 = pow(-temp_1, 0.25);
	  double erro_2 = 0.25 * erro_1 * abs(temp_2) / abs(temp_1);
	  th1d_cv24->SetBinContent(i, temp_2);
	  th1d_cv24->SetBinError(i, erro_2);
	}
    }

  TH1D* th1d_cv22 = (TH1D*)th1d_gapp->Clone("th1d_cv22");
  for (int i = 1; i < bin_size+1; i++)
    {
      double temp_1 = th1d_cv22->GetBinContent(i);
      double erro_1 = th1d_cv22->GetBinError(i);

      if (temp_1 <= 0) th1d_cv22->SetBinContent(i, -9999);
      else
	{
	  double temp_2 = pow(temp_1, 0.5);
	  double erro_2 = 0.5 * erro_1 * abs(temp_2) / abs(temp_1);
	  th1d_cv22->SetBinContent(i, temp_2);
	  th1d_cv22->SetBinError(i, erro_2);
	}
    }

  // ------------------------------------------------------------------------------
  TH1D* th1d_v2_mid = (TH1D*)th1d_c24->Clone("th1d_v2_mid");
  TH1D* th1d_v2_sigma = (TH1D*)th1d_c24->Clone("th1d_v2_sigma");
  for (int i = 1; i < bin_size+1; i++)
    {
      double temp_1 = th1d_cv24->GetBinContent(i);
      double temp_2 = th1d_cv22->GetBinContent(i);

      if (temp_1 <= 0)
	{
	  th1d_v2_mid->SetBinContent(i, temp_1);
	  th1d_v2_sigma->SetBinContent(i, temp_1);
	}
      else
	{
	  double v2_mid = (temp_1 * temp_1 + temp_2 * temp_2) / 2;
	  double v2_sigma = (temp_2 * temp_2 - temp_1 * temp_1) / 2;
	  v2_mid = TMath::Sqrt(v2_mid);
	  if (v2_sigma > 0) v2_sigma = TMath::Sqrt(v2_sigma)/v2_mid;
	  else v2_sigma = -9999;
	  th1d_v2_mid->SetBinContent(i, v2_mid);
	  th1d_v2_sigma->SetBinContent(i, v2_sigma);
	}
    }


  // -----------------------------------------------------------------------------------
  TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
  c1->Divide(1, 2, 0, 0);

  TF1 *axis1 = new TF1("axis1", "0", 0, x_max);
  TF1 *axis2 = new TF1("axis2", "0", 0, x_max);

  //------------------------------------------------
  //Pad 1
  c1->cd(1);
  gPad->SetTicky();
  gPad->SetTickx();
  gPad->SetTopMargin(0.05);
  gPad->SetPad(.005, .3, .9, .92 );

  axis1->SetTitle("");
  axis1->Draw("AXIS");
  axis1->GetXaxis()->SetTitle("N_char");
  axis1->GetYaxis()->SetTitle("v_{2}");
  axis1->GetXaxis()->SetTitleFont(62);
  axis1->GetYaxis()->SetRangeUser(0, 0.1);
  axis1->GetYaxis()->SetLabelSize(0.05);
  axis1->GetYaxis()->SetTitleSize(0.07);
  axis1->GetXaxis()->SetLabelSize(0.05);
  axis1->GetXaxis()->SetTitleSize(0.05);
  axis1->GetYaxis()->SetTitleOffset(0.6);
  axis1->GetYaxis()->SetLabelFont(62);
  axis1->GetYaxis()->SetTitleFont(62);

  ppplane->SetMarkerStyle(33);
  ppplane->SetMarkerColor(1);
  ppplane->SetMarkerSize(1.5);

  th1d_v2_mid->SetMarkerStyle(27);
  th1d_v2_mid->SetMarkerColor(1);
  th1d_v2_mid->SetMarkerSize(1.5);

  th1d_cv22->SetMarkerStyle(24);
  th1d_cv24->SetMarkerStyle(25);
  th1d_cv22->SetMarkerColor(kRed);
  th1d_cv24->SetMarkerColor(kBlue);

  th1d_cv24->Draw("p,same");
  ppplane->Draw("p,same");
  th1d_cv22->Draw("p,same");
  th1d_v2_mid->Draw("p,same");

  TLegend *leg1 = new TLegend(0.61,0.03,0.96,0.23);
  leg1->AddEntry(th1d_cv24, "4 particle cumulant", "lep");
  leg1->AddEntry(th1d_cv22, "2 particle cumulant", "lep");
  leg1->AddEntry(ppplane, "participant plane with partons", "lep");
  leg1->AddEntry(th1d_v2_mid, "v_{2} from v_{2}{2} and v_{2}{4}", "lep");
  leg1->Draw("same");

  //------------------------------------------------
  //Pad 2
  c1->cd(2);
  gPad->SetTicky();
  gPad->SetTickx();
  gPad->SetPad(.005, .005, .9, .3);
  gPad->SetBottomMargin(0.3);

  TF1 *one = new TF1("one", "1", 0, x_max);

  axis2->SetTitle("");
  axis2->Draw("AXIS");
  axis2->GetXaxis()->SetTitle("N_char");
  axis2->GetXaxis()->SetLabelSize(0.10);
  axis2->GetXaxis()->SetTitleSize(0.12);
  axis2->GetYaxis()->SetLabelSize(0.10);
  axis2->GetYaxis()->SetTitleSize(0.10);
  axis2->GetXaxis()->SetLabelFont(62);
  axis2->GetXaxis()->SetTitleFont(62);
  axis2->GetYaxis()->SetRangeUser(0.7, 1.3);
  axis2->GetYaxis()->SetTitleOffset(0.5);
  axis2->GetYaxis()->SetLabelFont(62);
  axis2->GetYaxis()->SetTitleFont(62);

  one->SetLineWidth(1);
  one->SetLineStyle(7);
  one->Draw("same");

  TH1D* ratio_v2mid_ppp = (TH1D*)th1d_v2_mid->Clone("ratio_v2mid_ppp");
  TH1D* ratio_v22_ppp = (TH1D*)th1d_cv22->Clone("ratio_v22_ppp");
  TH1D* ratio_v24_ppp = (TH1D*)th1d_cv24->Clone("ratio_v24_ppp");
  for (int i = 1; i < bin_size+1; i++)
    {
      double temp_v2_sigma = th1d_v2_mid->GetBinContent(i);
      double temp_cv22 = th1d_cv22->GetBinContent(i);
      double temp_cv24 = th1d_cv24->GetBinContent(i);
      double temp_2 = ppplane->GetBinContent(i);

      if (temp_2 <= 0) ratio_v2mid_ppp->SetBinContent(i, -9999);
      else
	{

	  ratio_v2mid_ppp->SetBinContent(i, temp_v2_sigma / temp_2);
	  ratio_v22_ppp->SetBinContent(i, temp_cv22 / temp_2);
	  ratio_v24_ppp->SetBinContent(i, temp_cv24 / temp_2);
	}
    }

  ratio_v2mid_ppp->SetMarkerStyle(27);
  ratio_v2mid_ppp->SetMarkerColor(1);
  ratio_v2mid_ppp->SetMarkerSize(1.5);

  ratio_v22_ppp->SetMarkerStyle(24);
  ratio_v22_ppp->SetMarkerColor(kRed);

  ratio_v24_ppp->SetMarkerStyle(25);
  ratio_v24_ppp->SetMarkerColor(kBlue);

  ratio_v2mid_ppp->Draw("p,same");
  ratio_v22_ppp->Draw("p,same");
  ratio_v24_ppp->Draw("p,same");

  c1->Print("ampt_v22v24etc.png");
  c1->Print("ampt_v22v24etc.pdf");

  // ---------------------------------------------------------------------------------
  // ---------------------------------------------------------------------------------
  TCanvas *c2 = new TCanvas("c2", "c2", 800, 600);
  c2->Divide(1, 2, 0, 0);
  TF1 *axis3 = new TF1("axis3", "0", 0, x_max);
  TF1 *axis4 = new TF1("axis4", "0", 0, x_max);

  //------------------------------------------------
  //Pad 1
  c2->cd(1);
  gPad->SetTicky();
  gPad->SetTickx();
  gPad->SetTopMargin(0.05);
  gPad->SetPad(.005, .3, .9, .92 );

  axis3->SetTitle("");
  axis3->Draw("AXIS");
  axis3->GetXaxis()->SetTitle("N_char");
  axis3->GetXaxis()->SetTitleFont(62);
  axis3->GetYaxis()->SetRangeUser(0, 1);
  axis3->GetYaxis()->SetLabelSize(0.05);
  axis3->GetYaxis()->SetTitleSize(0.05);
  axis3->GetXaxis()->SetLabelSize(0.05);
  axis3->GetXaxis()->SetTitleSize(0.05);
  axis3->GetYaxis()->SetTitleOffset(0.5);
  axis3->GetYaxis()->SetLabelFont(62);
  axis3->GetYaxis()->SetTitleFont(62);

  th1d_epsilon2->SetMarkerStyle(22);
  th1d_epsilon2->SetMarkerColor(1);

  th1d_v2_sigma->SetMarkerStyle(26);
  th1d_v2_sigma->SetMarkerColor(1);

  th1d_epsilon2->Draw("p,same");
  th1d_v2_sigma->Draw("p,same");

  TLegend *leg3 = new TLegend(0.61,0.03,0.96,0.23);
  leg3->AddEntry(th1d_epsilon2, "#sigma_{#epsilon_{2}}/#epsilon_{2}", "lep");
  leg3->AddEntry(th1d_v2_sigma, "#sigma_{v_{2}}/v_{2}", "lep");
  leg3->Draw("same");

  //------------------------------------------------
  //Pad 2
  c2->cd(2);
  gPad->SetTicky();
  gPad->SetTickx();
  gPad->SetPad(.005, .005, .9, .3);
  gPad->SetBottomMargin(0.3);

  axis4->SetTitle("");
  axis4->Draw("AXIS");
  axis4->GetXaxis()->SetTitle("N_char");
  axis4->GetXaxis()->SetLabelSize(0.10);
  axis4->GetXaxis()->SetTitleSize(0.12);
  axis4->GetYaxis()->SetLabelSize(0.10);
  axis4->GetYaxis()->SetTitleSize(0.10);
  axis4->GetXaxis()->SetLabelFont(62);
  axis4->GetXaxis()->SetTitleFont(62);
  axis4->GetYaxis()->SetRangeUser(0.7, 1.3);
  axis4->GetYaxis()->SetTitleOffset(0.5);
  axis4->GetYaxis()->SetLabelFont(62);
  axis4->GetYaxis()->SetTitleFont(62);

  one->Draw("same");

  TH1D* ratio_ep2_v2 = (TH1D*)th1d_epsilon2->Clone("ratio_ep2_v2");
  for (int i = 1; i < bin_size+1; i++)
    {
      double temp_sigep2 = th1d_epsilon2->GetBinContent(i);
      double temp_sigv2  = th1d_v2_sigma->GetBinContent(i);

      ratio_ep2_v2->SetBinContent(i, temp_sigep2 / temp_sigv2);
    }

  ratio_ep2_v2->SetMarkerStyle(22);
  ratio_ep2_v2->SetMarkerColor(1);

  ratio_ep2_v2->Draw("p,same");

  // 	TFile *fout = new TFile("ep_pbpb.root", "RECREATE");

  // 	th1d_epsilon2->Write();
  // 	th1d_v2_sigma->Write();

  // 	fout->Close();


  c2->Print("ampt_sigma_epsilon_v.png");
  c2->Print("ampt_sigma_epsilon_v.pdf");

  TFile* fout = TFile::Open("epsilon_for_cumulants.root","recreate");
  th1d_epsilon2->SetName("th1d_epsilon2_200");
  th1d_epsilon2->Write();
  fout->Write();
  fout->Close();
}
