#include <iostream>
#include "TROOT.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TStyle.h"

const float pi = 3.1415926535;

TH1D* pleasefixme3h(TH1D*);
TH1D* pleasefixme2h(TH1D*);

void domain();
void dosomethingelse();
void donew();

TH1D* dosimplefix(TH1D*);

void delta()
{

  //  donew();
  domain();
  //  dosomethingelse();

}


void dosomethingelse()
{

  TCanvas* c1 = new TCanvas("c1","");

  gStyle->SetOptTitle(1);

  TFile* file = TFile::Open("input/combined_200.root");

  TH1D* th1dA = (TH1D*)file->Get("fvtxs_v3_both_phipsi");
  TH1D* th1d0 = (TH1D*)file->Get("fvtxs0_v3_both_phipsi");
  TH1D* th1d1 = (TH1D*)file->Get("fvtxs1_v3_both_phipsi");
  TH1D* th1d2 = (TH1D*)file->Get("fvtxs2_v3_both_phipsi");
  TH1D* th1d3 = (TH1D*)file->Get("fvtxs3_v3_both_phipsi");

  TH1D* th1dAE = (TH1D*)file->Get("fvtxs_v3_east_phipsi");
  TH1D* th1d0E = (TH1D*)file->Get("fvtxs0_v3_east_phipsi");
  TH1D* th1d1E = (TH1D*)file->Get("fvtxs1_v3_east_phipsi");
  TH1D* th1d2E = (TH1D*)file->Get("fvtxs2_v3_east_phipsi");
  TH1D* th1d3E = (TH1D*)file->Get("fvtxs3_v3_east_phipsi");

  TH1D* th1dAW = (TH1D*)file->Get("fvtxs_v3_west_phipsi");
  TH1D* th1d0W = (TH1D*)file->Get("fvtxs0_v3_west_phipsi");
  TH1D* th1d1W = (TH1D*)file->Get("fvtxs1_v3_west_phipsi");
  TH1D* th1d2W = (TH1D*)file->Get("fvtxs2_v3_west_phipsi");
  TH1D* th1d3W = (TH1D*)file->Get("fvtxs3_v3_west_phipsi");

  TH1D* th1dB = (TH1D*)file->Get("bbcs_v3_both_phipsi");
  TH1D* th1dBE = (TH1D*)file->Get("bbcs_v3_east_phipsi");
  TH1D* th1dBW = (TH1D*)file->Get("bbcs_v3_west_phipsi");


  th1dB->Sumw2();
  th1dB->Draw();
  th1dB->SetLineColor(kBlack);
  th1dBE->SetLineColor(kBlue);
  th1dBW->SetLineColor(kRed);
  th1dB->SetTitle("BBC South");
  th1dB->GetXaxis()->SetTitle("#phi-#psi_{3}");
  th1dBE->Draw("same");
  th1dBW->Draw("same");
  TLegend legA(0.18,0.68,0.38,0.88);
  legA.AddEntry(th1dB,"both arms","el");
  legA.AddEntry(th1dBE,"east arm","el");
  legA.AddEntry(th1dBW,"west arm","el");
  legA.SetTextSize(0.045);
  legA.Sumw2();
  legA.Draw();
  c1->Print("FigsEventPlane/superduper_bbcsAEW.png");
  c1->Print("FigsEventPlane/superduper_bbcsAEW.pdf");

  th1dA->Sumw2();
  th1dA->Draw();
  th1dA->SetLineColor(kBlack);
  th1dAE->SetLineColor(kBlue);
  th1dAW->SetLineColor(kRed);
  th1dA->SetTitle("FVTXS all layers");
  th1dA->GetXaxis()->SetTitle("#phi-#psi_{3}");
  th1dAE->Draw("same");
  th1dAW->Draw("same");
  TLegend legA(0.18,0.68,0.38,0.88);
  legA.AddEntry(th1dA,"both arms","el");
  legA.AddEntry(th1dAE,"east arm","el");
  legA.AddEntry(th1dAW,"west arm","el");
  legA.SetTextSize(0.045);
  legA.Sumw2();
  legA.Draw();
  c1->Print("FigsEventPlane/superduper_fvtxsAEW.png");
  c1->Print("FigsEventPlane/superduper_fvtxsAEW.pdf");

  th1d0->Sumw2();
  th1d0->Draw();
  th1d0->SetLineColor(kBlack);
  th1d0E->SetLineColor(kBlue);
  th1d0W->SetLineColor(kRed);
  th1d0->SetTitle("FVTXS layer 0");
  th1d0->GetXaxis()->SetTitle("#phi-#psi_{3}");
  th1d0E->Draw("same");
  th1d0W->Draw("same");
  TLegend leg1(0.18,0.68,0.38,0.88);
  leg1.AddEntry(th1d0,"both arms","el");
  leg1.AddEntry(th1d0E,"east arm","el");
  leg1.AddEntry(th1d0W,"west arm","el");
  leg1.SetTextSize(0.045);
  leg1.Sumw2();
  leg1.Draw();
  c1->Print("FigsEventPlane/superduper_fvtxs0EW.png");
  c1->Print("FigsEventPlane/superduper_fvtxs0EW.pdf");

  th1d1->Sumw2();
  th1d1->Draw();
  th1d1->SetLineColor(kBlack);
  th1d1E->SetLineColor(kBlue);
  th1d1W->SetLineColor(kRed);
  th1d1->SetTitle("FVTXS layer 1");
  th1d1->GetXaxis()->SetTitle("#phi-#psi_{3}");
  th1d1E->Draw("same");
  th1d1W->Draw("same");
  TLegend leg2(0.18,0.68,0.38,0.88);
  leg2.AddEntry(th1d1,"both arms","el");
  leg2.AddEntry(th1d1E,"east arm","el");
  leg2.AddEntry(th1d1W,"west arm","el");
  leg2.SetTextSize(0.045);
  leg2.Sumw2();
  leg2.Draw();
  c1->Print("FigsEventPlane/superduper_fvtxs1EW.png");
  c1->Print("FigsEventPlane/superduper_fvtxs1EW.pdf");

  th1d2->Sumw2();
  th1d2->Draw();
  th1d2->SetLineColor(kBlack);
  th1d2E->SetLineColor(kBlue);
  th1d2W->SetLineColor(kRed);
  th1d2->SetTitle("FVTXS layer 2");
  th1d2->GetXaxis()->SetTitle("#phi-#psi_{3}");
  th1d2E->Draw("same");
  th1d2W->Draw("same");
  TLegend leg3(0.18,0.68,0.38,0.88);
  leg3.AddEntry(th1d2,"both arms","el");
  leg3.AddEntry(th1d2E,"east arm","el");
  leg3.AddEntry(th1d2W,"west arm","el");
  leg3.SetTextSize(0.045);
  leg3.Sumw2();
  leg3.Draw();
  c1->Print("FigsEventPlane/superduper_fvtxs2EW.png");
  c1->Print("FigsEventPlane/superduper_fvtxs2EW.pdf");

  th1d3->Sumw2();
  th1d3->Draw();
  th1d3->SetLineColor(kBlack);
  th1d3E->SetLineColor(kBlue);
  th1d3W->SetLineColor(kRed);
  th1d3->SetTitle("FVTXS layer 3");
  th1d3->GetXaxis()->SetTitle("#phi-#psi_{3}");
  th1d3E->Draw("same");
  th1d3W->Draw("same");
  TLegend leg4(0.18,0.68,0.38,0.88);
  leg4.AddEntry(th1d3,"both arms","el");
  leg4.AddEntry(th1d3E,"east arm","el");
  leg4.AddEntry(th1d3W,"west arm","el");
  leg4.SetTextSize(0.045);
  leg4.Sumw2();
  leg4.Draw();
  c1->Print("FigsEventPlane/superduper_fvtxs3EW.png");
  c1->Print("FigsEventPlane/superduper_fvtxs3EW.pdf");

  cout << "triple checking here " << endl;
  TH1D* th1d_B = (TH1D*)file->Get("th1d_os_dreso3_BBC_CNT");
  TH1D* th1d_S = (TH1D*)file->Get("th1d_os_dreso3_CNT_FVTX");
  th1d_B->Rebin(2);
  cout << "integral for bbc is " << th1d_B->Integral() << endl;
  cout << "integral before is " << th1d_S->Integral() << endl;
  th1d_S->Rebin(2);
  cout << "integral after is " << th1d_S->Integral() << endl;

  th1d_S->Sumw2();
  th1d_S->Draw();
  th1dA->SetLineColor(kRed);
  th1dA->Draw("same");
  cout << "integral is " << th1dA->Integral() << endl;
  c1->Print("triple_check_binning.png");

  th1d_S->Divide(th1d3);
  th1d_S->Sumw2();
  th1d_S->Draw();
  c1->Print("triple_check_binning_ratio.png");

  cout << "check binning here " << endl;
  TH1D* th1d_B_fixed = pleasefixme3h(th1d_B);
  TH1D* th1d_S_fixed = pleasefixme3h(th1d_S);
  // th1d_B_fixed->Rebin(2);
  // th1d_S_fixed->Rebin(2);

  th1d_B->Sumw2();
  th1d_B->Draw();
  th1d_B->SetTitle("BBCS");
  th1d_B->GetXaxis()->SetTitle("#phi-#Psi_{3}");
  c1->Print("FigsEventPlane/superduper_bbcs.png");
  c1->Print("FigsEventPlane/superduper_bbcs.pdf");

  th1d_S->SetLineColor(kGreen+2);
  th1d_S->Draw("same");
  th1d_B->SetTitle("BBCS and FVTXS");
  c1->Print("FigsEventPlane/superduper_bbcsfvtxs.png");
  c1->Print("FigsEventPlane/superduper_bbcsfvtxs.pdf");

  th1d_B_fixed->Sumw2();
  th1d_B_fixed->Draw();
  double max = th1d_B_fixed->GetBinContent(th1d_B_fixed->GetMaximumBin());
  th1d_B_fixed->SetMaximum(max*1.002);
  th1d_B_fixed->SetMinimum(max*0.990);
  th1d_B_fixed->GetXaxis()->SetRangeUser(-1.2,1.2);
  th1d_B_fixed->GetXaxis()->SetTitle("#phi-#Psi_{3}");
  TF1* fun = new TF1("fun","[0]+[1]*TMath::Cos(3*(x-[2]))",-2,2);
  fun->SetParameter(0,max);
  fun->SetParameter(1,max*0.0001);
  fun->SetParameter(2,0);
  th1d_B_fixed->Fit(fun,"","",-1.1,1.1);
  TLatex tex(0.5,0.8,Form("offset %.2e #pm %.2e",fun->GetParameter(2),fun->GetParError(2)));
  tex.SetNDC();
  tex.DrawLatex(0.45,0.82,Form("offset %.2e #pm %.2e",fun->GetParameter(2),fun->GetParError(2)));
  c1->Print("FigsEventPlane/superduper_bbcs_fixed.png");
  c1->Print("FigsEventPlane/superduper_bbcs_fixed.pdf");

  th1d_S_fixed->SetLineColor(kGreen+2);
  th1d_S_fixed->Draw("same");
  c1->Print("FigsEventPlane/superduper_bbcsfvtxs_fixed.png");
  c1->Print("FigsEventPlane/superduper_bbcsfvtxs_fixed.pdf");

}

void domain()
{

  gStyle->SetOptTitle(1);

  TCanvas* c1 = new TCanvas("c1","");

  TFile* file = TFile::Open("input/combined_200.root");


  TH1D* th1dA = (TH1D*)file->Get("fvtxs_v3_both_phipsi");
  TH1D* th1d0 = (TH1D*)file->Get("fvtxs0_v3_both_phipsi");
  TH1D* th1d1 = (TH1D*)file->Get("fvtxs1_v3_both_phipsi");
  TH1D* th1d2 = (TH1D*)file->Get("fvtxs2_v3_both_phipsi");
  TH1D* th1d3 = (TH1D*)file->Get("fvtxs3_v3_both_phipsi");

  TH1D* th1dAE = (TH1D*)file->Get("fvtxs_v3_east_phipsi");
  TH1D* th1d0E = (TH1D*)file->Get("fvtxs0_v3_east_phipsi");
  TH1D* th1d1E = (TH1D*)file->Get("fvtxs1_v3_east_phipsi");
  TH1D* th1d2E = (TH1D*)file->Get("fvtxs2_v3_east_phipsi");
  TH1D* th1d3E = (TH1D*)file->Get("fvtxs3_v3_east_phipsi");

  TH1D* th1dAW = (TH1D*)file->Get("fvtxs_v3_west_phipsi");
  TH1D* th1d0W = (TH1D*)file->Get("fvtxs0_v3_west_phipsi");
  TH1D* th1d1W = (TH1D*)file->Get("fvtxs1_v3_west_phipsi");
  TH1D* th1d2W = (TH1D*)file->Get("fvtxs2_v3_west_phipsi");
  TH1D* th1d3W = (TH1D*)file->Get("fvtxs3_v3_west_phipsi");

  cout << "check binning here " << endl;
  TH1D* th1dA_fixed = pleasefixme3h(th1dA);
  th1dA_fixed->Sumw2();
  th1dA_fixed->Draw();
  double max = th1dA_fixed->GetBinContent(th1dA_fixed->GetMaximumBin());
  th1dA_fixed->SetMaximum(max*1.002);
  th1dA_fixed->SetMinimum(max*0.995);
  th1dA_fixed->GetXaxis()->SetRangeUser(-1.2,1.2);
  th1dA_fixed->GetXaxis()->SetTitle("#phi-#Psi_{3}");
  TF1* fun = new TF1("fun","[0]+[1]*TMath::Cos(3*(x-[2]))",-2,2);
  //fun->SetParLimits(2,-0.6,0.6); // root is being very stupid with these fits
  fun->SetParLimits(2,-1.1,1.1);
  TLatex tex(0.5,0.8,Form("offset %.2e #pm %.2e",fun->GetParameter(2),fun->GetParError(2)));
  tex.SetNDC();
  fun->SetParameter(0,max);
  fun->SetParameter(1,max*0.0001);
  fun->SetParameter(2,0);
  th1dA_fixed->Fit(fun,"","",-1.1,1.1);
  tex.DrawLatex(0.45,0.82,Form("offset %.2e #pm %.2e",(fun->GetParameter(2)),(fun->GetParError(2))));
  TLine* line = new TLine(max*0.995,fun->GetParameter(2),max*1.002,fun->GetParameter(2));
  line->SetLineStyle(2);
  line->DrawLine(max*0.995,fun->GetParameter(2),max*1.002,fun->GetParameter(2));
  th1dA_fixed->SetTitle("FVTX South all layers");
  c1->Print("FigsEventPlane/phipsi3_th1dA.png");
  c1->Print("FigsEventPlane/phipsi3_th1dA.pdf");
  //fun->SetParLimits(2,-1.1,1.1);

  TH1D* th1dAE_fixed = pleasefixme3h(th1dAE);
  th1dAE_fixed->Sumw2();
  th1dAE_fixed->Draw();
  max = th1dAE_fixed->GetBinContent(th1dAE_fixed->GetMaximumBin());
  th1dAE_fixed->SetMaximum(max*1.002);
  th1dAE_fixed->SetMinimum(max*0.995);
  th1dAE_fixed->GetXaxis()->SetRangeUser(-1.2,1.2);
  th1dAE_fixed->GetXaxis()->SetTitle("#phi-#Psi_{3}");
  tex.SetNDC();
  fun->SetParameter(0,max);
  fun->SetParameter(1,max*0.0001);
  fun->SetParameter(2,0);
  th1dAE_fixed->Fit(fun,"","",-1.1,1.1);
  tex.DrawLatex(0.45,0.82,Form("offset %.2e #pm %.2e",fun->GetParameter(2),fun->GetParError(2)));
  c1->Print("FigsEventPlane/phipsi3_th1dAE.png");
  c1->Print("FigsEventPlane/phipsi3_th1dAE.pdf");

  TH1D* th1dAW_fixed = pleasefixme3h(th1dAW);
  th1dAW_fixed->Sumw2();
  th1dAW_fixed->Draw();
  max = th1dAW_fixed->GetBinContent(th1dAW_fixed->GetMaximumBin());
  th1dAW_fixed->SetMaximum(max*1.002);
  th1dAW_fixed->SetMinimum(max*0.995);
  th1dAW_fixed->GetXaxis()->SetRangeUser(-1.2,1.2);
  th1dAW_fixed->GetXaxis()->SetTitle("#phi-#Psi_{3}");
  tex.SetNDC();
  fun->SetParameter(0,max);
  fun->SetParameter(1,max*0.0001);
  fun->SetParameter(2,0);
  th1dAW_fixed->Fit(fun,"","",-1.1,1.1);
  tex.DrawLatex(0.45,0.82,Form("offset %.2e #pm %.2e",fun->GetParameter(2),fun->GetParError(2)));
  c1->Print("FigsEventPlane/phipsi3_th1dAW.png");
  c1->Print("FigsEventPlane/phipsi3_th1dAW.pdf");


  TH1D* th1dB = (TH1D*)file->Get("bbcs_v3_both_phipsi");
  TH1D* th1dBE = (TH1D*)file->Get("bbcs_v3_east_phipsi");
  TH1D* th1dBW = (TH1D*)file->Get("bbcs_v3_west_phipsi");

  TH1D* th1dB_fixed = pleasefixme3h(th1dB);
  th1dB_fixed->Sumw2();
  th1dB_fixed->Draw("hist e");
  max = th1dB_fixed->GetBinContent(th1dB_fixed->GetMaximumBin());
  th1dB_fixed->SetMaximum(max*1.002);
  th1dB_fixed->SetMinimum(max*0.995);
  th1dB_fixed->GetXaxis()->SetRangeUser(-1.2,1.2);
  th1dB_fixed->GetXaxis()->SetTitle("#phi-#Psi_{3}");
  tex.SetNDC();
  fun->SetParameter(0,max);
  fun->SetParameter(1,max*0.0001);
  fun->SetParameter(2,0);
  th1dB_fixed->Fit(fun,"","",-1.1,1.1);
  tex.DrawLatex(0.45,0.82,Form("offset %.2e #pm %.2e",fun->GetParameter(2),fun->GetParError(2)));
  line->DrawLine(max*0.995,fun->GetParameter(2),max*1.002,fun->GetParameter(2));
  th1dB_fixed->SetTitle("BBC South");
  c1->Print("FigsEventPlane/phipsi3_th1dB.png");
  c1->Print("FigsEventPlane/phipsi3_th1dB.pdf");

  TH1D* th1dBE_fixed = pleasefixme3h(th1dBE);
  th1dBE_fixed->Sumw2();
  th1dBE_fixed->Draw("hist e");
  max = th1dBE_fixed->GetBinContent(th1dBE_fixed->GetMaximumBin());
  th1dBE_fixed->SetMaximum(max*1.002);
  th1dBE_fixed->SetMinimum(max*0.995);
  th1dBE_fixed->GetXaxis()->SetRangeUser(-1.2,1.2);
  th1dBE_fixed->GetXaxis()->SetTitle("#phi-#Psi_{3}");
  tex.SetNDC();
  fun->SetParameter(0,max);
  fun->SetParameter(1,max*0.0001);
  fun->SetParameter(2,0);
  th1dBE_fixed->Fit(fun,"","",-1.1,1.1);
  tex.DrawLatex(0.45,0.82,Form("offset %.2e #pm %.2e",fun->GetParameter(2),fun->GetParError(2)));
  c1->Print("FigsEventPlane/phipsi3_th1dBE.png");
  c1->Print("FigsEventPlane/phipsi3_th1dBE.pdf");

  TH1D* th1dBW_fixed = pleasefixme3h(th1dBW);
  th1dBW_fixed->Sumw2();
  th1dBW_fixed->Draw("hist e");
  max = th1dBW_fixed->GetBinContent(th1dBW_fixed->GetMaximumBin());
  th1dBW_fixed->SetMaximum(max*1.002);
  th1dBW_fixed->SetMinimum(max*0.995);
  th1dBW_fixed->GetXaxis()->SetRangeUser(-1.2,1.2);
  th1dBW_fixed->GetXaxis()->SetTitle("#phi-#Psi_{3}");
  tex.SetNDC();
  fun->SetParameter(0,max);
  fun->SetParameter(1,max*0.0001);
  fun->SetParameter(2,0);
  th1dBW_fixed->Fit(fun,"","",-1.1,1.1);
  tex.DrawLatex(0.45,0.82,Form("offset %.2e #pm %.2e",fun->GetParameter(2),fun->GetParError(2)));
  c1->Print("FigsEventPlane/phipsi3_th1dBW.png");
  c1->Print("FigsEventPlane/phipsi3_th1dBW.pdf");


  TH1D* th1d0_fixed = pleasefixme3h(th1d0);
  th1d0_fixed->Sumw2();
  th1d0_fixed->Draw();
  max = th1d0_fixed->GetBinContent(th1d0_fixed->GetMaximumBin());
  th1d0_fixed->SetMaximum(max*1.002);
  th1d0_fixed->SetMinimum(max*0.995);
  th1d0_fixed->GetXaxis()->SetRangeUser(-1.2,1.2);
  th1d0_fixed->GetXaxis()->SetTitle("#phi-#Psi_{3}");
  fun->SetParameter(0,max);
  fun->SetParameter(1,max*0.0001);
  fun->SetParameter(2,0);
  th1d0_fixed->Fit(fun,"","",-1.1,1.1);
  tex.DrawLatex(0.45,0.82,Form("offset %.2e #pm %.2e",fun->GetParameter(2),fun->GetParError(2)));
  th1d0_fixed->SetTitle("FVTX South layer 0");
  c1->Print("FigsEventPlane/phipsi3_th1d0.png");
  c1->Print("FigsEventPlane/phipsi3_th1d0.pdf");

  TH1D* th1d0E_fixed = pleasefixme3h(th1d0E);
  th1d0E_fixed->Sumw2();
  th1d0E_fixed->Draw();
  max = th1d0E_fixed->GetBinContent(th1d0E_fixed->GetMaximumBin());
  th1d0E_fixed->SetMaximum(max*1.002);
  th1d0E_fixed->SetMinimum(max*0.995);
  th1d0E_fixed->GetXaxis()->SetRangeUser(-1.2,1.2);
  th1d0E_fixed->GetXaxis()->SetTitle("#phi-#Psi_{3}");
  tex.SetNDC();
  fun->SetParameter(0,max);
  fun->SetParameter(1,max*0.0001);
  fun->SetParameter(2,0);
  th1d0E_fixed->Fit(fun,"","",-1.1,1.1);
  tex.DrawLatex(0.45,0.82,Form("offset %.2e #pm %.2e",fun->GetParameter(2),fun->GetParError(2)));
  c1->Print("FigsEventPlane/phipsi3_th1d0E.png");
  c1->Print("FigsEventPlane/phipsi3_th1d0E.pdf");

  TH1D* th1d0W_fixed = pleasefixme3h(th1d0W);
  th1d0W_fixed->Sumw2();
  th1d0W_fixed->Draw();
  max = th1d0W_fixed->GetBinContent(th1d0W_fixed->GetMaximumBin());
  th1d0W_fixed->SetMaximum(max*1.002);
  th1d0W_fixed->SetMinimum(max*0.995);
  th1d0W_fixed->GetXaxis()->SetRangeUser(-1.2,1.2);
  th1d0W_fixed->GetXaxis()->SetTitle("#phi-#Psi_{3}");
  tex.SetNDC();
  fun->SetParameter(0,max);
  fun->SetParameter(1,max*0.0001);
  fun->SetParameter(2,0);
  th1d0W_fixed->Fit(fun,"","",-1.1,1.1);
  tex.DrawLatex(0.45,0.82,Form("offset %.2e #pm %.2e",fun->GetParameter(2),fun->GetParError(2)));
  c1->Print("FigsEventPlane/phipsi3_th1d0W.png");
  c1->Print("FigsEventPlane/phipsi3_th1d0W.pdf");

  TH1D* th1d1_fixed = pleasefixme3h(th1d1);
  th1d1_fixed->Sumw2();
  th1d1_fixed->Draw();
  max = th1d1_fixed->GetBinContent(th1d1_fixed->GetMaximumBin());
  th1d1_fixed->SetMaximum(max*1.002);
  th1d1_fixed->SetMinimum(max*0.995);
  th1d1_fixed->GetXaxis()->SetRangeUser(-1.2,1.2);
  th1d1_fixed->GetXaxis()->SetTitle("#phi-#Psi_{3}");
  fun->SetParameter(0,max);
  fun->SetParameter(1,max*0.0001);
  fun->SetParameter(2,0);
  th1d1_fixed->Fit(fun,"","",-1.1,1.1);
  tex.DrawLatex(0.45,0.82,Form("offset %.2e #pm %.2e",fun->GetParameter(2),fun->GetParError(2)));
  th1d1_fixed->SetTitle("FVTX South layer 1");
  c1->Print("FigsEventPlane/phipsi3_th1d1.png");
  c1->Print("FigsEventPlane/phipsi3_th1d1.pdf");

  TH1D* th1d1E_fixed = pleasefixme3h(th1d1E);
  th1d1E_fixed->Sumw2();
  th1d1E_fixed->Draw();
  max = th1d1E_fixed->GetBinContent(th1d1E_fixed->GetMaximumBin());
  th1d1E_fixed->SetMaximum(max*1.002);
  th1d1E_fixed->SetMinimum(max*0.995);
  th1d1E_fixed->GetXaxis()->SetRangeUser(-1.2,1.2);
  th1d1E_fixed->GetXaxis()->SetTitle("#phi-#Psi_{3}");
  tex.SetNDC();
  fun->SetParameter(0,max);
  fun->SetParameter(1,max*0.0001);
  fun->SetParameter(2,0);
  th1d1E_fixed->Fit(fun,"","",-1.1,1.1);
  tex.DrawLatex(0.45,0.82,Form("offset %.2e #pm %.2e",fun->GetParameter(2),fun->GetParError(2)));
  c1->Print("FigsEventPlane/phipsi3_th1d1E.png");
  c1->Print("FigsEventPlane/phipsi3_th1d1E.pdf");

  TH1D* th1d1W_fixed = pleasefixme3h(th1d1W);
  th1d1W_fixed->Sumw2();
  th1d1W_fixed->Draw();
  max = th1d1W_fixed->GetBinContent(th1d1W_fixed->GetMaximumBin());
  th1d1W_fixed->SetMaximum(max*1.002);
  th1d1W_fixed->SetMinimum(max*0.995);
  th1d1W_fixed->GetXaxis()->SetRangeUser(-1.2,1.2);
  th1d1W_fixed->GetXaxis()->SetTitle("#phi-#Psi_{3}");
  tex.SetNDC();
  fun->SetParameter(0,max);
  fun->SetParameter(1,max*0.0001);
  fun->SetParameter(2,0);
  th1d1W_fixed->Fit(fun,"","",-1.1,1.1);
  tex.DrawLatex(0.45,0.82,Form("offset %.2e #pm %.2e",fun->GetParameter(2),fun->GetParError(2)));
  c1->Print("FigsEventPlane/phipsi3_th1d1W.png");
  c1->Print("FigsEventPlane/phipsi3_th1d1W.pdf");

  TH1D* th1d2_fixed = pleasefixme3h(th1d2);
  th1d2_fixed->Sumw2();
  th1d2_fixed->Draw();
  max = th1d2_fixed->GetBinContent(th1d2_fixed->GetMaximumBin());
  th1d2_fixed->SetMaximum(max*1.002);
  th1d2_fixed->SetMinimum(max*0.995);
  th1d2_fixed->GetXaxis()->SetRangeUser(-1.2,1.2);
  th1d2_fixed->GetXaxis()->SetTitle("#phi-#Psi_{3}");
  fun->SetParameter(0,max);
  fun->SetParameter(1,max*0.0001);
  fun->SetParameter(2,0);
  th1d2_fixed->Fit(fun,"","",-1.1,1.1);
  tex.DrawLatex(0.45,0.82,Form("offset %.2e #pm %.2e",fun->GetParameter(2),fun->GetParError(2)));
  th1d2_fixed->SetTitle("FVTX South layer 2");
  c1->Print("FigsEventPlane/phipsi3_th1d2.png");
  c1->Print("FigsEventPlane/phipsi3_th1d2.pdf");

  TH1D* th1d2E_fixed = pleasefixme3h(th1d2E);
  th1d2E_fixed->Sumw2();
  th1d2E_fixed->Draw();
  max = th1d2E_fixed->GetBinContent(th1d2E_fixed->GetMaximumBin());
  th1d2E_fixed->SetMaximum(max*1.002);
  th1d2E_fixed->SetMinimum(max*0.995);
  th1d2E_fixed->GetXaxis()->SetRangeUser(-1.2,1.2);
  th1d2E_fixed->GetXaxis()->SetTitle("#phi-#Psi_{3}");
  tex.SetNDC();
  fun->SetParameter(0,max);
  fun->SetParameter(1,max*0.0001);
  fun->SetParameter(2,0);
  th1d2E_fixed->Fit(fun,"","",-1.1,1.1);
  tex.DrawLatex(0.45,0.82,Form("offset %.2e #pm %.2e",fun->GetParameter(2),fun->GetParError(2)));
  c1->Print("FigsEventPlane/phipsi3_th1d2E.png");
  c1->Print("FigsEventPlane/phipsi3_th1d2E.pdf");

  TH1D* th1d2W_fixed = pleasefixme3h(th1d2W);
  th1d2W_fixed->Sumw2();
  th1d2W_fixed->Draw();
  max = th1d2W_fixed->GetBinContent(th1d2W_fixed->GetMaximumBin());
  th1d2W_fixed->SetMaximum(max*1.002);
  th1d2W_fixed->SetMinimum(max*0.995);
  th1d2W_fixed->GetXaxis()->SetRangeUser(-1.2,1.2);
  th1d2W_fixed->GetXaxis()->SetTitle("#phi-#Psi_{3}");
  tex.SetNDC();
  fun->SetParameter(0,max);
  fun->SetParameter(1,max*0.0001);
  fun->SetParameter(2,0);
  th1d2W_fixed->Fit(fun,"","",-1.1,1.1);
  tex.DrawLatex(0.45,0.82,Form("offset %.2e #pm %.2e",fun->GetParameter(2),fun->GetParError(2)));
  c1->Print("FigsEventPlane/phipsi3_th1d2W.png");
  c1->Print("FigsEventPlane/phipsi3_th1d2W.pdf");

  TH1D* th1d3_fixed = pleasefixme3h(th1d3);
  th1d3_fixed->Sumw2();
  th1d3_fixed->Draw();
  max = th1d3_fixed->GetBinContent(th1d3_fixed->GetMaximumBin());
  th1d3_fixed->SetMaximum(max*1.002);
  th1d3_fixed->SetMinimum(max*0.995);
  th1d3_fixed->GetXaxis()->SetRangeUser(-1.2,1.2);
  th1d3_fixed->GetXaxis()->SetTitle("#phi-#Psi_{3}");
  fun->SetParameter(0,max);
  fun->SetParameter(1,max*0.0001);
  fun->SetParameter(2,0);
  th1d3_fixed->Fit(fun,"","",-1.1,1.1);
  tex.DrawLatex(0.45,0.82,Form("offset %.2e #pm %.2e",fun->GetParameter(2),fun->GetParError(2)));
  th1d3_fixed->SetTitle("FVTX South layer 3");
  c1->Print("FigsEventPlane/phipsi3_th1d3.png");
  c1->Print("FigsEventPlane/phipsi3_th1d3.pdf");

  TH1D* th1d3E_fixed = pleasefixme3h(th1d3E);
  th1d3E_fixed->Sumw2();
  th1d3E_fixed->Draw();
  max = th1d3E_fixed->GetBinContent(th1d3E_fixed->GetMaximumBin());
  th1d3E_fixed->SetMaximum(max*1.002);
  th1d3E_fixed->SetMinimum(max*0.995);
  th1d3E_fixed->GetXaxis()->SetRangeUser(-1.2,1.2);
  th1d3E_fixed->GetXaxis()->SetTitle("#phi-#Psi_{3}");
  tex.SetNDC();
  fun->SetParameter(0,max);
  fun->SetParameter(1,max*0.0001);
  fun->SetParameter(2,0);
  th1d3E_fixed->Fit(fun,"","",-1.1,1.1);
  tex.DrawLatex(0.45,0.82,Form("offset %.2e #pm %.2e",fun->GetParameter(2),fun->GetParError(2)));
  c1->Print("FigsEventPlane/phipsi3_th1d3E.png");
  c1->Print("FigsEventPlane/phipsi3_th1d3E.pdf");

  TH1D* th1d3W_fixed = pleasefixme3h(th1d3W);
  th1d3W_fixed->Sumw2();
  th1d3W_fixed->Draw();
  max = th1d3W_fixed->GetBinContent(th1d3W_fixed->GetMaximumBin());
  th1d3W_fixed->SetMaximum(max*1.002);
  th1d3W_fixed->SetMinimum(max*0.995);
  th1d3W_fixed->GetXaxis()->SetRangeUser(-1.2,1.2);
  th1d3W_fixed->GetXaxis()->SetTitle("#phi-#Psi_{3}");
  tex.SetNDC();
  fun->SetParameter(0,max);
  fun->SetParameter(1,max*0.0001);
  fun->SetParameter(2,0);
  th1d3W_fixed->Fit(fun,"","",-1.1,1.1);
  tex.DrawLatex(0.45,0.82,Form("offset %.2e #pm %.2e",fun->GetParameter(2),fun->GetParError(2)));
  c1->Print("FigsEventPlane/phipsi3_th1d3W.png");
  c1->Print("FigsEventPlane/phipsi3_th1d3W.pdf");

  // ---

  TH1D* th1d_B = (TH1D*)file->Get("th1d_os_dreso2_BBC_CNT");
  TH1D* th1d_S = (TH1D*)file->Get("th1d_os_dreso2_CNT_FVTX");

  TF1* fun2 = new TF1("fun2","[0]+[1]*TMath::Cos(2*x-[2])",-2,2);

  TH1D* th1d_B_fixed = pleasefixme2h(th1d_B);
  TH1D* th1d_S_fixed = pleasefixme2h(th1d_S);

  th1d_B_fixed->Sumw2();
  th1d_B_fixed->Draw();
  max = th1d_B_fixed->GetBinContent(th1d_B_fixed->GetMaximumBin());
  th1d_B_fixed->SetMaximum(max*1.002);
  th1d_B_fixed->SetMinimum(max*0.95);
  th1d_B_fixed->GetXaxis()->SetRangeUser(-1.7,1.7);
  th1d_B_fixed->GetXaxis()->SetTitle("#phi-#Psi_{2}");
  th1d_B_fixed->Fit(fun2,"","",-1.6,1.6);
  tex.DrawLatex(0.45,0.82,Form("offset %.2e #pm %.2e",fun->GetParameter(2),fun->GetParError(2)));
  c1->Print("FigsEventPlane/phipsi2_th1d_B.png");
  c1->Print("FigsEventPlane/phipsi2_th1d_B.pdf");

  th1d_S_fixed->Sumw2();
  th1d_S_fixed->Draw();
  max = th1d_S_fixed->GetBinContent(th1d_S_fixed->GetMaximumBin());
  th1d_S_fixed->SetMaximum(max*1.002);
  th1d_S_fixed->SetMinimum(max*0.92);
  th1d_S_fixed->GetXaxis()->SetRangeUser(-1.7,1.7);
  th1d_S_fixed->GetXaxis()->SetTitle("#phi-#Psi_{2}");
  th1d_S_fixed->Fit(fun2,"","",-1.6,1.6);
  tex.DrawLatex(0.45,0.82,Form("offset %.2e #pm %.2e",fun->GetParameter(2),fun->GetParError(2)));
  c1->Print("FigsEventPlane/phipsi2_th1d_S.png");
  c1->Print("FigsEventPlane/phipsi2_th1d_S.pdf");

  // ----------------------------------

  TH1D* th1d_2_A = (TH1D*)file->Get("fvtxs_v2_both_phipsi");
  TH1D* th1d_2_0 = (TH1D*)file->Get("fvtxs0_v2_both_phipsi");
  TH1D* th1d_2_1 = (TH1D*)file->Get("fvtxs1_v2_both_phipsi");
  TH1D* th1d_2_2 = (TH1D*)file->Get("fvtxs2_v2_both_phipsi");
  TH1D* th1d_2_3 = (TH1D*)file->Get("fvtxs3_v2_both_phipsi");

  TH1D* th1d_2_A_fixed = pleasefixme2h(th1d_2_A);
  th1d_2_A_fixed->Sumw2();
  th1d_2_A_fixed->Draw();
  max = th1d_2_A_fixed->GetBinContent(th1d_2_A_fixed->GetMaximumBin());
  th1d_2_A_fixed->SetMaximum(max*1.002);
  th1d_2_A_fixed->SetMinimum(max*0.92);
  th1d_2_A_fixed->GetXaxis()->SetRangeUser(-1.7,1.7);
  th1d_2_A_fixed->GetXaxis()->SetTitle("#phi-#Psi_{3}");
  delete fun;
  fun = new TF1("fun","[0]+[1]*TMath::Cos(2*x-[2])",-2,2);
  th1d_2_A_fixed->Fit(fun,"","",-1.6,1.6);
  tex.DrawLatex(0.45,0.82,Form("offset %.2e #pm %.2e",fun->GetParameter(2),fun->GetParError(2)));
  c1->Print("FigsEventPlane/phipsi2_th1dA.png");
  c1->Print("FigsEventPlane/phipsi2_th1dA.pdf");

  TH1D* th1d_2_0_fixed = pleasefixme2h(th1d_2_0);
  th1d_2_0_fixed->Sumw2();
  th1d_2_0_fixed->Draw();
  max = th1d_2_0_fixed->GetBinContent(th1d_2_0_fixed->GetMaximumBin());
  th1d_2_0_fixed->SetMaximum(max*1.002);
  th1d_2_0_fixed->SetMinimum(max*0.92);
  th1d_2_0_fixed->GetXaxis()->SetRangeUser(-1.7,1.7);
  th1d_2_0_fixed->GetXaxis()->SetTitle("#phi-#Psi_{3}");
  th1d_2_0_fixed->Fit(fun,"","",-1.6,1.6);
  tex.DrawLatex(0.45,0.82,Form("offset %.2e #pm %.2e",fun->GetParameter(2),fun->GetParError(2)));
  c1->Print("FigsEventPlane/phipsi2_th1d0.png");
  c1->Print("FigsEventPlane/phipsi2_th1d0.pdf");

  TH1D* th1d_2_1_fixed = pleasefixme2h(th1d_2_1);
  th1d_2_1_fixed->Sumw2();
  th1d_2_1_fixed->Draw();
  max = th1d_2_1_fixed->GetBinContent(th1d_2_1_fixed->GetMaximumBin());
  th1d_2_1_fixed->SetMaximum(max*1.002);
  th1d_2_1_fixed->SetMinimum(max*0.92);
  th1d_2_1_fixed->GetXaxis()->SetRangeUser(-1.7,1.7);
  th1d_2_1_fixed->GetXaxis()->SetTitle("#phi-#Psi_{3}");
  th1d_2_1_fixed->Fit(fun,"","",-1.6,1.6);
  tex.DrawLatex(0.45,0.82,Form("offset %.2e #pm %.2e",fun->GetParameter(2),fun->GetParError(2)));
  c1->Print("FigsEventPlane/phipsi2_th1d1.png");
  c1->Print("FigsEventPlane/phipsi2_th1d1.pdf");

  TH1D* th1d_2_2_fixed = pleasefixme2h(th1d_2_2);
  th1d_2_2_fixed->Sumw2();
  th1d_2_2_fixed->Draw();
  max = th1d_2_2_fixed->GetBinContent(th1d_2_2_fixed->GetMaximumBin());
  th1d_2_2_fixed->SetMaximum(max*1.002);
  th1d_2_2_fixed->SetMinimum(max*0.92);
  th1d_2_2_fixed->GetXaxis()->SetRangeUser(-1.7,1.7);
  th1d_2_2_fixed->GetXaxis()->SetTitle("#phi-#Psi_{3}");
  th1d_2_2_fixed->Fit(fun,"","",-1.6,1.6);
  tex.DrawLatex(0.45,0.82,Form("offset %.2e #pm %.2e",fun->GetParameter(2),fun->GetParError(2)));
  c1->Print("FigsEventPlane/phipsi2_th1d2.png");
  c1->Print("FigsEventPlane/phipsi2_th1d2.pdf");

  TH1D* th1d_2_3_fixed = pleasefixme2h(th1d_2_3);
  th1d_2_3_fixed->Sumw2();
  th1d_2_3_fixed->Draw();
  max = th1d_2_3_fixed->GetBinContent(th1d_2_3_fixed->GetMaximumBin());
  th1d_2_3_fixed->SetMaximum(max*1.002);
  th1d_2_3_fixed->SetMinimum(max*0.92);
  th1d_2_3_fixed->GetXaxis()->SetRangeUser(-1.7,1.7);
  th1d_2_3_fixed->GetXaxis()->SetTitle("#phi-#Psi_{3}");
  th1d_2_3_fixed->Fit(fun,"","",-1.6,1.6);
  tex.DrawLatex(0.45,0.82,Form("offset %.2e #pm %.2e",fun->GetParameter(2),fun->GetParError(2)));
  c1->Print("FigsEventPlane/phipsi2_th1d3.png");
  c1->Print("FigsEventPlane/phipsi2_th1d3.pdf");

  // --- now event more fun stuff

  TH1D* th1d012 = (TH1D*)file->Get("fvtxs012_v3_both_phipsi");
  TH1D* th1d013 = (TH1D*)file->Get("fvtxs013_v3_both_phipsi");
  TH1D* th1d023 = (TH1D*)file->Get("fvtxs023_v3_both_phipsi");
  TH1D* th1d123 = (TH1D*)file->Get("fvtxs123_v3_both_phipsi");

  TH1D* th1d012_fixed = pleasefixme3h(th1d012);
  th1d012_fixed->Sumw2();
  th1d012_fixed->Draw();
  max = th1d012_fixed->GetBinContent(th1d012_fixed->GetMaximumBin());
  th1d012_fixed->SetMaximum(max*1.002);
  th1d012_fixed->SetMinimum(max*0.995);
  th1d012_fixed->GetXaxis()->SetRangeUser(-1.2,1.2);
  th1d012_fixed->GetXaxis()->SetTitle("#phi-#Psi_{3}");
  fun->SetParameter(0,max);
  fun->SetParameter(1,max*0.0001);
  fun->SetParameter(2,0);
  th1d012_fixed->Fit(fun,"","",-1.1,1.1);
  tex.DrawLatex(0.45,0.82,Form("offset %.2e #pm %.2e",fun->GetParameter(2),fun->GetParError(2)));
  c1->Print("FigsEventPlane/phipsi3_th1d012.png");
  c1->Print("FigsEventPlane/phipsi3_th1d012.pdf");

  TH1D* th1d013_fixed = pleasefixme3h(th1d013);
  th1d013_fixed->Sumw2();
  th1d013_fixed->Draw();
  max = th1d013_fixed->GetBinContent(th1d013_fixed->GetMaximumBin());
  th1d013_fixed->SetMaximum(max*1.002);
  th1d013_fixed->SetMinimum(max*0.995);
  th1d013_fixed->GetXaxis()->SetRangeUser(-1.2,1.2);
  th1d013_fixed->GetXaxis()->SetTitle("#phi-#Psi_{3}");
  fun->SetParameter(0,max);
  fun->SetParameter(1,max*0.0001);
  fun->SetParameter(2,0);
  th1d013_fixed->Fit(fun,"","",-1.1,1.1);
  tex.DrawLatex(0.45,0.82,Form("offset %.2e #pm %.2e",fun->GetParameter(2),fun->GetParError(2)));
  c1->Print("FigsEventPlane/phipsi3_th1d013.png");
  c1->Print("FigsEventPlane/phipsi3_th1d013.pdf");

  TH1D* th1d023_fixed = pleasefixme3h(th1d023);
  th1d023_fixed->Sumw2();
  th1d023_fixed->Draw();
  max = th1d023_fixed->GetBinContent(th1d023_fixed->GetMaximumBin());
  th1d023_fixed->SetMaximum(max*1.002);
  th1d023_fixed->SetMinimum(max*0.995);
  th1d023_fixed->GetXaxis()->SetRangeUser(-1.2,1.2);
  th1d023_fixed->GetXaxis()->SetTitle("#phi-#Psi_{3}");
  fun->SetParameter(0,max);
  fun->SetParameter(1,max*0.0001);
  fun->SetParameter(2,0);
  th1d023_fixed->Fit(fun,"","",-1.1,1.1);
  tex.DrawLatex(0.45,0.82,Form("offset %.2e #pm %.2e",fun->GetParameter(2),fun->GetParError(2)));
  c1->Print("FigsEventPlane/phipsi3_th1d023.png");
  c1->Print("FigsEventPlane/phipsi3_th1d023.pdf");

  TH1D* th1d123_fixed = pleasefixme3h(th1d123);
  th1d123_fixed->Sumw2();
  th1d123_fixed->Draw();
  max = th1d123_fixed->GetBinContent(th1d123_fixed->GetMaximumBin());
  th1d123_fixed->SetMaximum(max*1.002);
  th1d123_fixed->SetMinimum(max*0.995);
  th1d123_fixed->GetXaxis()->SetRangeUser(-1.2,1.2);
  th1d123_fixed->GetXaxis()->SetTitle("#phi-#Psi_{3}");
  fun->SetParameter(0,max);
  fun->SetParameter(1,max*0.0001);
  fun->SetParameter(2,0);
  th1d123_fixed->Fit(fun,"","",-1.1,1.1);
  tex.DrawLatex(0.45,0.82,Form("offset %.2e #pm %.2e",fun->GetParameter(2),fun->GetParError(2)));
  c1->Print("FigsEventPlane/phipsi3_th1d123.png");
  c1->Print("FigsEventPlane/phipsi3_th1d123.pdf");

  // --- these are all messed up for some reason, maybe come back to this some other time

  th1d_B = (TH1D*)file->Get("th1d_os_dreso3_BBC_CNT");
  th1d_S = (TH1D*)file->Get("th1d_os_dreso3_CNT_FVTX");
  th1d_B_fixed = pleasefixme3h(th1d_B);
  th1d_S_fixed = pleasefixme3h(th1d_S);

  th1d_B_fixed->Sumw2();
  th1d_B_fixed->Draw();
  max = th1d_B_fixed->GetBinContent(th1d_B_fixed->GetMaximumBin());
  th1d_B_fixed->SetMaximum(max*1.002);
  th1d_B_fixed->SetMinimum(max*0.995);
  th1d_B_fixed->GetXaxis()->SetRangeUser(-1.2,1.2);
  th1d_B_fixed->GetXaxis()->SetTitle("#phi-#Psi_{3}");
  fun->SetParameter(0,max);
  fun->SetParameter(1,max*0.0001);
  fun->SetParameter(2,0);
  th1d_B_fixed->Fit(fun,"","",-1.1,1.1);
  tex.DrawLatex(0.45,0.82,Form("offset %.2e #pm %.2e",fun->GetParameter(2),fun->GetParError(2)));
  c1->Print("FigsEventPlane/phipsi3_th1d_B.png");
  c1->Print("FigsEventPlane/phipsi3_th1d_B.pdf");

  //  th1d_S_fixed->Rebin(2);
  th1d_S_fixed->Sumw2();
  th1d_S_fixed->Draw();
  max = th1d_S_fixed->GetBinContent(th1d_S_fixed->GetMaximumBin());
  th1d_S_fixed->SetMaximum(max*1.002);
  th1d_S_fixed->SetMinimum(max*0.995);
  th1d_S_fixed->GetXaxis()->SetRangeUser(-1.2,1.2);
  th1d_S_fixed->GetXaxis()->SetTitle("#phi-#Psi_{3}");
  fun->SetParameter(0,max);
  fun->SetParameter(1,max*0.0001);
  fun->SetParameter(2,0);
  th1d_S_fixed->Fit(fun,"","",-1.1,1.1);
  tex.DrawLatex(0.45,0.82,Form("offset %.2e #pm %.2e",fun->GetParameter(2),fun->GetParError(2)));
  c1->Print("FigsEventPlane/phipsi3_th1d_S.png");
  c1->Print("FigsEventPlane/phipsi3_th1d_S.pdf");



}

TH1D* pleasefixme3h(TH1D* in)
{
  TH1D* clone = (TH1D*)in->Clone();
  int nbins = clone->GetNbinsX();
  for ( int i = 0; i < nbins; ++i )
    {
      float phi = clone->GetBinCenter(i+1);
      float amount = clone->GetBinContent(i+1);
      if ( phi < -pi )
	{
	  if ( i < 5 )
	    {
	      cout << "phi is " << phi << " " << clone->GetBinLowEdge(i+1) << " " << clone->GetBinLowEdge(i+2) << endl;
	    }
	  clone->SetBinContent(i+1,0);
	  clone->Fill(phi+2*pi,amount);
	}
      if ( phi > pi )
	{
	  clone->SetBinContent(i+1,0);
	  clone->Fill(phi-2*pi,amount);
	}
    }
  TH1D* clone2 = (TH1D*)clone->Clone();
  for ( int i = 0; i < nbins; ++i )
    {
      float phi = clone2->GetBinCenter(i+1);
      float amount = clone2->GetBinContent(i+1);
      if ( phi < -pi/3 )
	{
	  clone2->SetBinContent(i+1,0);
	  clone2->Fill(phi+(2*pi/3),amount);
	}
      if ( phi > pi/3 )
	{
	  clone2->SetBinContent(i+1,0);
	  clone2->Fill(phi-(2*pi/3),amount);
	}
    }
  return clone2;
}


TH1D* pleasefixme2h(TH1D* in)
{
  TH1D* clone = (TH1D*)in->Clone(Form("%s_clone",in->GetName()));
  int nbins = clone->GetNbinsX();
  for ( int i = 0; i < nbins; ++i )
    {
      float phi = clone->GetBinCenter(i+1);
      float amount = clone->GetBinContent(i+1);
      if ( phi < -pi )
	{
	  clone->SetBinContent(i+1,0);
	  clone->Fill(phi+2*pi,amount);
	}
      if ( phi > pi )
	{
	  clone->SetBinContent(i+1,0);
	  clone->Fill(phi-2*pi,amount);
	}
    }
  TH1D* clone2 = (TH1D*)clone->Clone(Form("%s_clone2",in->GetName()));
  for ( int i = 0; i < nbins; ++i )
    {
      float phi = clone2->GetBinCenter(i+1);
      float amount = clone2->GetBinContent(i+1);
      if ( phi < -pi/2 )
	{
	  clone2->SetBinContent(i+1,0);
	  clone2->Fill(phi+pi,amount);
	}
      if ( phi > pi/2 )
	{
	  clone2->SetBinContent(i+1,0);
	  clone2->Fill(phi-pi,amount);
	}
    }
  return clone2;
}


void donew()
{

  TCanvas* c1 = new TCanvas("c1","");

  gStyle->SetOptTitle(1);

  TFile* file = TFile::Open("input/combined_200.root");

  TH1D* th1dA = (TH1D*)file->Get("fvtxs_v3_both_phipsi3");
  TH1D* th1d0 = (TH1D*)file->Get("fvtxs0_v3_both_phipsi3");
  TH1D* th1d1 = (TH1D*)file->Get("fvtxs1_v3_both_phipsi3");
  TH1D* th1d2 = (TH1D*)file->Get("fvtxs2_v3_both_phipsi3");
  TH1D* th1d3 = (TH1D*)file->Get("fvtxs3_v3_both_phipsi3");

  TH1D* th1dAE = (TH1D*)file->Get("fvtxs_v3_east_phipsi3");
  TH1D* th1d0E = (TH1D*)file->Get("fvtxs0_v3_east_phipsi3");
  TH1D* th1d1E = (TH1D*)file->Get("fvtxs1_v3_east_phipsi3");
  TH1D* th1d2E = (TH1D*)file->Get("fvtxs2_v3_east_phipsi3");
  TH1D* th1d3E = (TH1D*)file->Get("fvtxs3_v3_east_phipsi3");

  TH1D* th1dAW = (TH1D*)file->Get("fvtxs_v3_west_phipsi3");
  TH1D* th1d0W = (TH1D*)file->Get("fvtxs0_v3_west_phipsi3");
  TH1D* th1d1W = (TH1D*)file->Get("fvtxs1_v3_west_phipsi3");
  TH1D* th1d2W = (TH1D*)file->Get("fvtxs2_v3_west_phipsi3");
  TH1D* th1d3W = (TH1D*)file->Get("fvtxs3_v3_west_phipsi3");

  double max = 1.0;
  TF1* fun = new TF1("fun","[0]+[1]*TMath::Cos(x-[2])",-3.3,3.3);
  TF1* fun3 = new TF1("fun3","[0]+[1]*TMath::Cos(3*(x-[2]))",-1.1,1.1);
  TLatex tex(0.5,0.8,Form("offset %.2e #pm %.2e",fun->GetParameter(2),fun->GetParError(2)));
  tex.SetNDC();

  //TH1D* th1dA = th1dA;
  //TH1D* th1dA = dosimplefix(th1dA);
  th1dA->SetLineColor(kBlack);
  th1dA->Draw("ep");
  max = th1dA->GetBinContent(th1dA->GetMaximumBin());
  th1dA->SetMaximum(max*1.002);
  th1dA->SetMinimum(max*0.995);
  th1dA->SetTitle("FVTXS all layers");
  th1dA->GetXaxis()->SetTitle("#phi-#psi_{3}");
  //th1dA->GetXaxis()->SetLabelSize(0);
  //th1dA->GetXaxis()->SetLimits(-2.1,2.1); // shocking how badly this works in root
  //th1dA->GetXaxis()->SetRangeUser(-1.2,1.2);
  th1dA->GetXaxis()->SetRangeUser(-3.3,3.3);
  fun->SetParameter(0,max);
  fun->SetParameter(1,max*0.0001);
  fun->SetParameter(2,0);
  // fun3->SetParameter(0,max);
  // fun3->SetParameter(1,max*0.0001);
  // fun3->SetParameter(2,0);
  //th1dA->Fit(fun3,"","",-1.1,1.1);
  th1dA->Fit(fun,"","",-3.1,3.1);
  tex.DrawLatex(0.45,0.82,Form("offset %.2e #pm %.2e",fun->GetParameter(2),fun->GetParError(2)));
  //TGaxis *xaxis0 = new TGaxis(0.15,0.15,0.9,0.15,-1.067,1.067,510,"");
  TGaxis *xaxis0 = new TGaxis(-3.2,max*0.995,3.2,max*0.995,-1.067,1.067,500,"");
  xaxis0->SetName("xaxis0");
  xaxis0->SetTitle("#phi-#psi_{3}");
  xaxis0->SetLabelSize(0.05);
  xaxis0->SetTitleSize(0.05);
  xaxis0->SetTitleOffset(1.0);
  xaxis0->SetTickSize(0);
  //xaxis0->Sumw2();
  //xaxis0->Draw();
  //th1dA->GetXaxis()->SetRange(-1.1,1.1); // shocking how badly this works in root
  c1->Print("FigsEventPlane/superduperultramega_fvtxsA.png");
  c1->Print("FigsEventPlane/superduperultramega_fvtxsA.pdf");

  //return;

  //TH1D* th1d0 = th1d0;
  //TH1D* th1d0 = dosimplefix(th1d0);
  th1d0->SetLineColor(kBlack);
  th1d0->Draw("ep");
  max = th1d0->GetBinContent(th1d0->GetMaximumBin());
  th1d0->SetMaximum(max*1.002);
  th1d0->SetMinimum(max*0.995);
  th1d0->SetTitle("FVTXS layer 0");
  th1d0->GetXaxis()->SetTitle("#phi-#psi_{3}");
  //th1d0->GetXaxis()->SetLimits(-2.1,2.1); // shocking how badly this works in root
  //th1d0->GetXaxis()->SetRangeUser(-1.2,1.2);
  th1d0->GetXaxis()->SetRangeUser(-3.2,3.2);
  fun->SetParameter(0,max);
  fun->SetParameter(1,max*0.0001);
  fun->SetParameter(2,0);
  // fun3->SetParameter(0,max);
  // fun3->SetParameter(1,max*0.0001);
  // fun3->SetParameter(2,0);
  //th1d0->Fit(fun3,"","",-1.1,1.1);
  th1d0->Fit(fun,"","",-3.1,3.1);
  tex.DrawLatex(0.45,0.82,Form("offset %.2e #pm %.2e",fun->GetParameter(2),fun->GetParError(2)));
  c1->Print("FigsEventPlane/superduperultramega_fvtxs0.png");
  c1->Print("FigsEventPlane/superduperultramega_fvtxs0.pdf");

  //TH1D* th1d1 = th1d1;
  //TH1D* th1d1 = dosimplefix(th1d1);
  th1d1->SetLineColor(kBlack);
  th1d1->Draw("ep");
  max = th1d1->GetBinContent(th1d1->GetMaximumBin());
  th1d1->SetMaximum(max*1.002);
  th1d1->SetMinimum(max*0.995);
  th1d1->SetTitle("FVTXS layer 1");
  th1d1->GetXaxis()->SetTitle("#phi-#psi_{3}");
  //th1d1->GetXaxis()->SetLimits(-2.1,2.1); // shocking how badly this works in root
  //th1d1->GetXaxis()->SetRangeUser(-1.2,1.2);
  th1d1->GetXaxis()->SetRangeUser(-3.2,3.2);
  fun->SetParameter(0,max);
  fun->SetParameter(1,max*0.0001);
  fun->SetParameter(2,0);
  // fun3->SetParameter(0,max);
  // fun3->SetParameter(1,max*0.0001);
  // fun3->SetParameter(2,0);
  //th1d1->Fit(fun3,"","",-1.1,1.1);
  th1d1->Fit(fun,"","",-3.1,3.1);
  tex.DrawLatex(0.45,0.82,Form("offset %.2e #pm %.2e",fun->GetParameter(2),fun->GetParError(2)));
  c1->Print("FigsEventPlane/superduperultramega_fvtxs1.png");
  c1->Print("FigsEventPlane/superduperultramega_fvtxs1.pdf");

  //TH1D* th1d2 = th1d2;
  //TH1D* th1d2 = dosimplefix(th1d2);
  th1d2->SetLineColor(kBlack);
  th1d2->Draw("ep");
  max = th1d2->GetBinContent(th1d2->GetMaximumBin());
  th1d2->SetMaximum(max*1.002);
  th1d2->SetMinimum(max*0.995);
  th1d2->SetTitle("FVTXS layer 2");
  th1d2->GetXaxis()->SetTitle("#phi-#psi_{3}");
  //th1d2->GetXaxis()->SetLimits(-2.1,2.1); // shocking how badly this works in root
  //th1d2->GetXaxis()->SetRangeUser(-1.2,1.2);
  th1d2->GetXaxis()->SetRangeUser(-3.2,3.2);
  fun->SetParameter(0,max);
  fun->SetParameter(1,max*0.0001);
  fun->SetParameter(2,0);
  // fun3->SetParameter(0,max);
  // fun3->SetParameter(1,max*0.0001);
  // fun3->SetParameter(2,0);
  //th1d2->Fit(fun3,"","",-1.1,1.1);
  th1d2->Fit(fun,"","",-3.1,3.1);
  tex.DrawLatex(0.45,0.82,Form("offset %.2e #pm %.2e",fun->GetParameter(2),fun->GetParError(2)));
  c1->Print("FigsEventPlane/superduperultramega_fvtxs2.png");
  c1->Print("FigsEventPlane/superduperultramega_fvtxs2.pdf");

  //TH1D* th1d3 = th1d3;
  //TH1D* th1d3 = dosimplefix(th1d3);
  th1d3->SetLineColor(kBlack);
  th1d3->Draw("ep");
  max = th1d3->GetBinContent(th1d3->GetMaximumBin());
  th1d3->SetMaximum(max*1.002);
  th1d3->SetMinimum(max*0.995);
  th1d3->SetTitle("FVTXS layer 3");
  th1d3->GetXaxis()->SetTitle("#phi-#psi_{3}");
  //th1d3->GetXaxis()->SetLimits(-2.1,2.1); // shocking how badly this works in root
  //th1d3->GetXaxis()->SetRangeUser(-1.2,1.2);
  th1d3->GetXaxis()->SetRangeUser(-3.2,3.2);
  fun->SetParameter(0,max);
  fun->SetParameter(1,max*0.0001);
  fun->SetParameter(2,0);
  // fun3->SetParameter(0,max);
  // fun3->SetParameter(1,max*0.0001);
  // fun3->SetParameter(2,0);
  //th1d3->Fit(fun3,"","",-1.1,1.1);
  th1d3->Fit(fun,"","",-3.1,3.1);
  tex.DrawLatex(0.45,0.82,Form("offset %.2e #pm %.2e",fun->GetParameter(2),fun->GetParError(2)));
  c1->Print("FigsEventPlane/superduperultramega_fvtxs3.png");
  c1->Print("FigsEventPlane/superduperultramega_fvtxs3.pdf");

  return;

}

TH1D* dosimplefix(TH1D* in)
{
  int nbins = in->GetNbinsX();
  TH1D* htemp = new TH1D("temp","",nbins,-2.1,2.1);
  for ( int i = 0; i < nbins; ++i )
    {
      htemp->SetBinContent(i+1,in->GetBinContent(i+1));
      htemp->SetBinError(i+1,in->GetBinError(i+1));
    }
  htemp->SetName(Form("%s_append",in->GetName()));
  return htemp;
}
