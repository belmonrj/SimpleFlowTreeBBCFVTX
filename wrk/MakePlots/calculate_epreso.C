#include "../RpPar.h"

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TStyle.h>
#include <TString.h>
#include <TProfile.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TF1.h>
#include <TBox.h>

#include <iostream>
#include <utility>
#include <math.h>

using namespace std;

typedef pair<double, double> ValErr;

void doenergy(int, int);
ValErr calc_epreso(ValErr AB, ValErr AC, ValErr BC);

void calculate_epreso()
{

  doenergy(200, 2);
  doenergy(62, 2);
  doenergy(39, 2);
  // doenergy(20, 2);

  // doenergy(200,3);
  // doenergy(62,3);
  // doenergy(39,3);
  // doenergy(20,3);

}

void doenergy(int energy, int harmonic)
{

  cout << endl;
  cout << " ---------- " << endl;
  cout << " -- " << energy << endl;
  cout << " ---------- " << endl;
  // gStyle->SetOptTitle(1);

  // TFile* file = TFile::Open(Form("input/combined_%d.root",energy));
  TFile* file = TFile::Open(Form("../output/hist_dAu%d.root", energy));
  if (!file)
  {
    cout << "ERROR!! Unable to find file for energy " << energy << endl;
    return;
  }

  // ---
  // ---
  // ---

  TH1D* th1d_FVTXN_nclus = (TH1D*)file->Get("th1d_FVTXN_nclus");
  TH1D* th1d_FVTXS_nclus = (TH1D*)file->Get("th1d_FVTXS_nclus");
  TH1D* th1d_BBC_charge = (TH1D*)file->Get("th1d_BBC_charge");

  float mean_fvtxn = th1d_FVTXN_nclus->GetMean();
  float mean_fvtxs = th1d_FVTXS_nclus->GetMean();
  float mean_bbcs = th1d_BBC_charge->GetMean();

  TString dm = Form("%d GeV & %.2e & %.2e & %.2e \\\\",
                    energy, mean_bbcs, mean_fvtxs, mean_fvtxn);

  // ---
  // --- Now Centrality & zvrtxdependent
  // ---
  TProfile* tp1f_BBCS_FVTXN[NMUL][NZPS];
  TProfile* tp1f_BBCS_FVTXS[NMUL][NZPS];
  TProfile* tp1f_FVTXN_FVTXS[NMUL][NZPS];
  TProfile* tp1f_BBCS_CNT[NMUL][NZPS];
  TProfile* tp1f_CNT_FVTXS[NMUL][NZPS];
  TProfile* tp1f_CNT_FVTXN[NMUL][NZPS];

  TProfile* tp1f_BBCS_FVTXSA[NMUL][NZPS];
  TProfile* tp1f_BBCS_FVTXSB[NMUL][NZPS];
  TProfile* tp1f_BBCS_FVTXSL[NMUL][NZPS][NFVTXLAY];
  TProfile* tp1f_CNT_FVTXSA[NMUL][NZPS];
  TProfile* tp1f_CNT_FVTXSB[NMUL][NZPS];
  TProfile* tp1f_CNT_FVTXSL[NMUL][NZPS][NFVTXLAY];

  TProfile* tp1f_BBCS_FVTXN_sum[NMUL];
  TProfile* tp1f_BBCS_FVTXS_sum[NMUL];
  TProfile* tp1f_FVTXN_FVTXS_sum[NMUL];
  TProfile* tp1f_BBCS_CNT_sum[NMUL];
  TProfile* tp1f_CNT_FVTXS_sum[NMUL];
  TProfile* tp1f_CNT_FVTXN_sum[NMUL];

  TProfile* tp1f_BBCS_FVTXSA_sum[NMUL];
  TProfile* tp1f_BBCS_FVTXSB_sum[NMUL];
  TProfile* tp1f_BBCS_FVTXSL_sum[NMUL][NFVTXLAY];
  TProfile* tp1f_CNT_FVTXSA_sum[NMUL];
  TProfile* tp1f_CNT_FVTXSB_sum[NMUL];
  TProfile* tp1f_CNT_FVTXSL_sum[NMUL][NFVTXLAY];


  // ---
  // --- Resolutions
  // ---
  ValErr RBBCS[NMUL];
  ValErr RFVTXS[NMUL];
  ValErr RFVTXSA[NMUL];
  ValErr RFVTXSB[NMUL];
  ValErr RFVTXSL[NMUL][NFVTXLAY];

  // ---
  // --- Graphs for plotting
  // ---
  TGraphErrors *gRBBCS[NMUL];
  TGraphErrors *gRFVTXS[NMUL];
  TGraphErrors *gRFVTXSA[NMUL];
  TGraphErrors *gRFVTXSB[NMUL];
  TGraphErrors *gRFVTXSL[NMUL][NFVTXLAY];

  int detMarker[] =
  {kFullCircle, kFullSquare, kFullDiamond, kFullStar, kFullCross, kOpenCircle, kOpenSquare};
  float detSize[] =
  {1, 1, 1.5, 1.5, 1, 1, 1};
  int detColor[] =
  {kBlack, kBlue, kRed, kGreen + 2, kMagenta + 2, kYellow + 2, kOrange + 5};

  // TString dcor[NMUL];
  // TString dcnt[NMUL];
  // TString dfb[NMUL];
  // TString davg[NMUL];
  // TString drat[NMUL];


  TString dcnt[NMUL];
  TString dfvtxsab[NMUL];
  TString dfvtxsl[NMUL];

  for (int ic = 0; ic < NMUL; ic++)
  {
    gRBBCS[ic] = new TGraphErrors();
    gRBBCS[ic]->SetMarkerStyle(kOpenCircle);
    gRBBCS[ic]->SetMarkerSize(1.0);
    gRBBCS[ic]->SetMarkerColor(kRed);
    gRBBCS[ic]->SetLineColor(kRed);

    gRFVTXS[ic] = new TGraphErrors();
    gRFVTXS[ic]->SetMarkerStyle(detMarker[0]);
    gRFVTXS[ic]->SetMarkerSize(detSize[0]);
    gRFVTXS[ic]->SetMarkerColor(detColor[0]);
    gRFVTXS[ic]->SetLineColor(detColor[0]);

    for (int il = 0; il < NFVTXLAY; il++)
    {
      gRFVTXSL[ic][il] = new TGraphErrors();
      gRFVTXSL[ic][il]->SetMarkerStyle(detMarker[1 + il]);
      gRFVTXSL[ic][il]->SetMarkerSize(detSize[1 + il]);
      gRFVTXSL[ic][il]->SetMarkerColor(detColor[1 + il]);
      gRFVTXSL[ic][il]->SetLineColor(detColor[1 + il]);
    } // il

    gRFVTXSA[ic] = new TGraphErrors();
    gRFVTXSA[ic]->SetMarkerStyle(detMarker[NFVTXLAY + 2]);
    gRFVTXSA[ic]->SetMarkerSize(detSize[NFVTXLAY + 2]);
    gRFVTXSA[ic]->SetMarkerColor(detColor[NFVTXLAY + 2]);
    gRFVTXSA[ic]->SetLineColor(detColor[NFVTXLAY + 2]);

    gRFVTXSB[ic] = new TGraphErrors();
    gRFVTXSB[ic]->SetMarkerStyle(detMarker[NFVTXLAY + 3]);
    gRFVTXSB[ic]->SetMarkerSize(detSize[NFVTXLAY + 3]);
    gRFVTXSB[ic]->SetMarkerColor(detColor[NFVTXLAY + 3]);
    gRFVTXSB[ic]->SetLineColor(detColor[NFVTXLAY + 3]);

    for (int iz = 0; iz < NZPS; iz++)
    {

      // --- get histograms from file
      tp1f_BBCS_FVTXN[ic][iz] = (TProfile*)file->Get(Form("tp1f_c%d_z%d_reso%d_BBC_FVTXN", ic, iz, harmonic));
      if (!tp1f_BBCS_FVTXN[ic][iz])
      {
        cout << "ERROR!! Unable to find " << Form("tp1f_c%d_z%d_reso%d_BBC_FVTXN", ic, iz, harmonic)
             << " in file" << endl;
        return;
      }
      tp1f_BBCS_FVTXS[ic][iz] = (TProfile*)file->Get(Form("tp1f_c%d_z%d_reso%d_BBC_FVTX", ic, iz, harmonic));
      tp1f_FVTXN_FVTXS[ic][iz] = (TProfile*)file->Get(Form("tp1f_c%d_z%d_reso%d_FVTXS_FVTXN", ic, iz, harmonic));
      tp1f_BBCS_CNT[ic][iz] = (TProfile*)file->Get(Form("tp1f_c%d_z%d_reso%d_BBC_CNT", ic, iz, harmonic));
      tp1f_CNT_FVTXS[ic][iz] = (TProfile*)file->Get(Form("tp1f_c%d_z%d_reso%d_CNT_FVTX", ic, iz, harmonic));
      tp1f_CNT_FVTXN[ic][iz] = (TProfile*)file->Get(Form("tp1f_c%d_z%d_reso%d_CNT_FVTXN", ic, iz, harmonic));

      tp1f_BBCS_FVTXSA[ic][iz] = (TProfile*)file->Get(Form("tp1f_c%d_z%d_reso%d_BBC_FVTXSA", ic, iz, harmonic));
      tp1f_BBCS_FVTXSB[ic][iz] = (TProfile*)file->Get(Form("tp1f_c%d_z%d_reso%d_BBC_FVTXSB", ic, iz, harmonic));
      for (int il = 0; il < NFVTXLAY; il++)
        tp1f_BBCS_FVTXSL[ic][iz][il] = (TProfile*)file->Get(Form("tp1f_c%d_z%d_reso%d_BBC_FVTXSL%i", ic, iz, harmonic, il));

      tp1f_CNT_FVTXSA[ic][iz] = (TProfile*)file->Get(Form("tp1f_c%d_z%d_reso%d_CNT_FVTXSA", ic, iz, harmonic));
      tp1f_CNT_FVTXSB[ic][iz] = (TProfile*)file->Get(Form("tp1f_c%d_z%d_reso%d_CNT_FVTXSB", ic, iz, harmonic));
      for (int il = 0; il < NFVTXLAY; il++)
        tp1f_CNT_FVTXSL[ic][iz][il] = (TProfile*)file->Get(Form("tp1f_c%d_z%d_reso%d_CNT_FVTXSL%i", ic, iz, harmonic, il));

      // --- fill summed histograms
      if ( iz == 0 )
      {
        tp1f_BBCS_FVTXS_sum[ic] =
          (TProfile*) tp1f_BBCS_FVTXS[ic][iz]->Clone(Form("tp1f_BBCS_FVTXS_sum_c%i", ic));
        tp1f_BBCS_CNT_sum[ic] =
          (TProfile*) tp1f_BBCS_CNT[ic][iz]->Clone(Form("tp1f_BBCS_CNT_sum_c%i", ic));
        tp1f_CNT_FVTXS_sum[ic] =
          (TProfile*) tp1f_CNT_FVTXS[ic][iz]->Clone(Form("tp1f_CNT_FVTXS_sum_c%i", ic));

        tp1f_BBCS_FVTXSA_sum[ic] =
          (TProfile*) tp1f_BBCS_FVTXSA[ic][iz]->Clone(Form("tp1f_BBCS_FVTXSA_sum_c%i", ic));
        tp1f_BBCS_FVTXSB_sum[ic] =
          (TProfile*) tp1f_BBCS_FVTXSB[ic][iz]->Clone(Form("tp1f_BBCS_FVTXSB_sum_c%i", ic));
        for (int il = 0; il < NFVTXLAY; il++)
        {
          tp1f_BBCS_FVTXSL_sum[ic][il] =
            (TProfile*) tp1f_BBCS_FVTXSL[ic][iz][il]->Clone(Form("tp1f_BBCS_FVTXSL_sum_c%i_l%i", ic, il));
        }
        tp1f_CNT_FVTXSA_sum[ic] =
          (TProfile*) tp1f_CNT_FVTXSA[ic][iz]->Clone(Form("tp1f_CNT_FVTXSA_sum_c%i", ic));
        tp1f_CNT_FVTXSB_sum[ic] =
          (TProfile*) tp1f_CNT_FVTXSB[ic][iz]->Clone(Form("tp1f_CNT_FVTXSB_sum_c%i", ic));
        for (int il = 0; il < NFVTXLAY; il++)
        {
          tp1f_CNT_FVTXSL_sum[ic][il] =
            (TProfile*) tp1f_CNT_FVTXSL[ic][iz][il]->Clone(Form("tp1f_CNT_FVTXSL_sum_c%i_l%i", ic, il));
        }

      }
      else
      {
        tp1f_BBCS_CNT_sum[ic]->Add(tp1f_BBCS_CNT[ic][iz]);
        tp1f_BBCS_FVTXS_sum[ic]->Add(tp1f_BBCS_FVTXS[ic][iz]);
        tp1f_BBCS_FVTXSA_sum[ic]->Add(tp1f_BBCS_FVTXSA[ic][iz]);
        tp1f_BBCS_FVTXSB_sum[ic]->Add(tp1f_BBCS_FVTXSB[ic][iz]);
        for (int il = 0; il < NFVTXLAY; il++)
          tp1f_BBCS_FVTXSL_sum[ic][il]->Add(tp1f_BBCS_FVTXSL[ic][iz][il]);

        tp1f_CNT_FVTXS_sum[ic]->Add(tp1f_CNT_FVTXS[ic][iz]);
        tp1f_CNT_FVTXSA_sum[ic]->Add(tp1f_CNT_FVTXSA[ic][iz]);
        tp1f_CNT_FVTXSB_sum[ic]->Add(tp1f_CNT_FVTXSB[ic][iz]);
        for (int il = 0; il < NFVTXLAY; il++)
          tp1f_CNT_FVTXSL_sum[ic][il]->Add(tp1f_CNT_FVTXSL[ic][iz][il]);

      }


      //-- get means & errors & calculate resolution for CNT-Bakcward combos
      ValErr BBCS_FVTXN = make_pair(tp1f_BBCS_FVTXN[ic][iz]->GetBinContent(1),
                                    tp1f_BBCS_FVTXN[ic][iz]->GetBinError(1));
      ValErr BBCS_FVTXS = make_pair(tp1f_BBCS_FVTXS[ic][iz]->GetBinContent(1),
                                    tp1f_BBCS_FVTXS[ic][iz]->GetBinError(1));
      ValErr FVTXN_FVTXS = make_pair(tp1f_FVTXN_FVTXS[ic][iz]->GetBinContent(1),
                                     tp1f_FVTXN_FVTXS[ic][iz]->GetBinError(1));
      ValErr CNT_BBCS = make_pair(tp1f_BBCS_CNT[ic][iz]->GetBinContent(1),
                                  tp1f_BBCS_CNT[ic][iz]->GetBinError(1));
      ValErr CNT_FVTXS = make_pair(tp1f_CNT_FVTXS[ic][iz]->GetBinContent(1),
                                   tp1f_CNT_FVTXS[ic][iz]->GetBinError(1));
      ValErr CNT_FVTXN = make_pair(tp1f_CNT_FVTXN[ic][iz]->GetBinContent(1),
                                   tp1f_CNT_FVTXN[ic][iz]->GetBinError(1));

      ValErr BBCS_FVTXSA = make_pair(tp1f_BBCS_FVTXSA[ic][iz]->GetBinContent(1),
                                     tp1f_BBCS_FVTXSA[ic][iz]->GetBinError(1));
      ValErr BBCS_FVTXSB = make_pair(tp1f_BBCS_FVTXSB[ic][iz]->GetBinContent(1),
                                     tp1f_BBCS_FVTXSB[ic][iz]->GetBinError(1));
      ValErr BBCS_FVTXSL[NFVTXLAY];
      ValErr CNT_FVTXSA = make_pair(tp1f_CNT_FVTXSA[ic][iz]->GetBinContent(1),
                                    tp1f_CNT_FVTXSA[ic][iz]->GetBinError(1));
      ValErr CNT_FVTXSB = make_pair(tp1f_CNT_FVTXSB[ic][iz]->GetBinContent(1),
                                    tp1f_CNT_FVTXSB[ic][iz]->GetBinError(1));
      ValErr CNT_FVTXSL[NFVTXLAY];
      for (int il = 0; il < NFVTXLAY; il++)
      {
        BBCS_FVTXSL[il] = make_pair(tp1f_BBCS_FVTXSL[ic][iz][il]->GetBinContent(1),
                                    tp1f_BBCS_FVTXSL[ic][iz][il]->GetBinError(1));
        CNT_FVTXSL[il] = make_pair(tp1f_CNT_FVTXSL[ic][iz][il]->GetBinContent(1),
                                   tp1f_CNT_FVTXSL[ic][iz][il]->GetBinError(1));
      }




      ValErr RBBCS = calc_epreso(CNT_BBCS, BBCS_FVTXS, CNT_FVTXS);
      ValErr RFVTXS = calc_epreso(CNT_FVTXS, BBCS_FVTXS, CNT_BBCS);
      ValErr RFVTXN = calc_epreso(CNT_FVTXN, BBCS_FVTXN, CNT_BBCS);

      ValErr RFVTXSA = calc_epreso(CNT_FVTXSA, BBCS_FVTXSA, CNT_BBCS);
      ValErr RFVTXSB = calc_epreso(CNT_FVTXSB, BBCS_FVTXSB, CNT_BBCS);
      ValErr RFVTXSL[NFVTXLAY];
      for (int il = 0; il < NFVTXLAY; il++)
        RFVTXSL[il] = calc_epreso(CNT_FVTXSL[il], BBCS_FVTXSL[il], CNT_BBCS);


      // --- Fill tgraphs
      double z = -10. + (iz + 0.5) * 20. / (float)NZPS;

      gRBBCS[ic]->SetPoint(iz, z, RBBCS.first);
      gRBBCS[ic]->SetPointError(iz, 0, RBBCS.second);

      gRFVTXS[ic]->SetPoint(iz, z, RFVTXS.first);
      gRFVTXS[ic]->SetPointError(iz, 0, RFVTXS.second);

      gRFVTXSA[ic]->SetPoint(iz, z, RFVTXSA.first);
      gRFVTXSA[ic]->SetPointError(iz, 0, RFVTXSA.second);

      gRFVTXSB[ic]->SetPoint(iz, z, RFVTXSB.first);
      gRFVTXSB[ic]->SetPointError(iz, 0, RFVTXSB.second);

      for (int il = 0; il < NFVTXLAY; il++)
      {
        gRFVTXSL[ic][il]->SetPoint(iz, z, RFVTXSL[il].first);
        gRFVTXSL[ic][il]->SetPointError(iz, 0, RFVTXSL[il].second);
      }

    } // iz

    // --- get means & errors & calculate resolution for CNT-Bakcward combos
    // --- summed over z
    ValErr BBCS_FVTXS = make_pair(tp1f_BBCS_FVTXS_sum[ic]->GetBinContent(1),
                                  tp1f_BBCS_FVTXS_sum[ic]->GetBinError(1));
    ValErr CNT_BBCS = make_pair(tp1f_BBCS_CNT_sum[ic]->GetBinContent(1),
                                tp1f_BBCS_CNT_sum[ic]->GetBinError(1));
    ValErr CNT_FVTXS = make_pair(tp1f_CNT_FVTXS_sum[ic]->GetBinContent(1),
                                 tp1f_CNT_FVTXS_sum[ic]->GetBinError(1));

    ValErr BBCS_FVTXSA = make_pair(tp1f_BBCS_FVTXSA_sum[ic]->GetBinContent(1),
                                   tp1f_BBCS_FVTXSA_sum[ic]->GetBinError(1));
    ValErr BBCS_FVTXSB = make_pair(tp1f_BBCS_FVTXSB_sum[ic]->GetBinContent(1),
                                   tp1f_BBCS_FVTXSB_sum[ic]->GetBinError(1));
    ValErr BBCS_FVTXSL[NFVTXLAY];
    ValErr CNT_FVTXSA = make_pair(tp1f_CNT_FVTXSA_sum[ic]->GetBinContent(1),
                                  tp1f_CNT_FVTXSA_sum[ic]->GetBinError(1));
    ValErr CNT_FVTXSB = make_pair(tp1f_CNT_FVTXSB_sum[ic]->GetBinContent(1),
                                  tp1f_CNT_FVTXSB_sum[ic]->GetBinError(1));
    ValErr CNT_FVTXSL[NFVTXLAY];
    for (int il = 0; il < NFVTXLAY; il++)
    {
      BBCS_FVTXSL[il] = make_pair(tp1f_BBCS_FVTXSL_sum[ic][il]->GetBinContent(1),
                                  tp1f_BBCS_FVTXSL_sum[ic][il]->GetBinError(1));
      CNT_FVTXSL[il] = make_pair(tp1f_CNT_FVTXSL_sum[ic][il]->GetBinContent(1),
                                 tp1f_CNT_FVTXSL_sum[ic][il]->GetBinError(1));
    }




    RBBCS[ic] = calc_epreso(CNT_BBCS, BBCS_FVTXS, CNT_FVTXS);
    RFVTXS[ic] = calc_epreso(CNT_FVTXS, BBCS_FVTXS, CNT_BBCS);
    // ValErr RFVTXN = calc_epreso(CNT_FVTXN, BBCS_FVTXN, CNT_BBCS);

    RFVTXSA[ic] = calc_epreso(CNT_FVTXSA, BBCS_FVTXSA, CNT_BBCS);
    RFVTXSB[ic] = calc_epreso(CNT_FVTXSB, BBCS_FVTXSB, CNT_BBCS);
    for (int il = 0; il < NFVTXLAY; il++)
      RFVTXSL[ic][il] = calc_epreso(CNT_FVTXSL[il], BBCS_FVTXSL[il], CNT_BBCS);


    dcnt[ic] = Form("%d GeV & %d & %.3e $\\pm$ %.3e & %.3e $\\pm$ %.3e \\\\",
                    energy, ic,
                    RBBCS[ic].first, RBBCS[ic].second,
                    RFVTXS[ic].first, RFVTXS[ic].second
                   );

    dfvtxsab[ic] = Form("%d GeV & %d & %.3e$\\pm$%.3e & %.3e$\\pm$%.3e & %.3e$\\pm$%.3e & %.3e$\\pm$%.3e \\\\",
                        energy, ic,
                        RBBCS[ic].first, RBBCS[ic].second,
                        RFVTXS[ic].first, RFVTXS[ic].second,
                        RFVTXSA[ic].first, RFVTXSA[ic].second,
                        RFVTXSB[ic].first, RFVTXSB[ic].second
                       );

    dfvtxsl[ic] = Form("%d GeV & %d & %.2e$\\pm$%.2e & %.2e$\\pm$%.2e & %.2e$\\pm$%.2e & %.2e$\\pm$%.2e & %.2e$\\pm$%.2e \\\\",
                       energy, ic,
                       RFVTXS[ic].first, RFVTXS[ic].second,
                       RFVTXSL[ic][0].first, RFVTXSL[ic][0].second,
                       RFVTXSL[ic][1].first, RFVTXSL[ic][1].second,
                       RFVTXSL[ic][2].first, RFVTXSL[ic][2].second,
                       RFVTXSL[ic][3].first, RFVTXSL[ic][3].second
                      );


    // // --- event planes and correlations using FVTXN-FVTXS-BBCS
    // float reso_BBCS_fn  = sqrt((float_BBCS_FVTXN * float_BBCS_FVTXS) / float_FVTXN_FVTXS); // BNBS/NS
    // float reso_FVTXS_fn = sqrt((float_FVTXN_FVTXS * float_BBCS_FVTXS) / float_BBCS_FVTXN); // NSBS/BN
    // float reso_FVTXN_fn = sqrt((float_FVTXN_FVTXS * float_BBCS_FVTXN) / float_BBCS_FVTXS); // NSBN/BS

    // float ereso_BBCS_fn  = sqrt( ( efloat_BBCS_FVTXN * efloat_BBCS_FVTXN / 4 * float_BBCS_FVTXN )
    //                              + ( efloat_BBCS_FVTXS * efloat_BBCS_FVTXS / 4 * float_BBCS_FVTXS )
    //                              + ( efloat_FVTXN_FVTXS * efloat_FVTXN_FVTXS / 4 * pow(float_FVTXN_FVTXS, 3) ) );
    // float ereso_FVTXS_fn = sqrt( ( efloat_FVTXN_FVTXS * efloat_FVTXN_FVTXS / 4 * float_FVTXN_FVTXS )
    //                              + ( efloat_BBCS_FVTXS * efloat_BBCS_FVTXS / 4 * float_BBCS_FVTXS )
    //                              + ( efloat_BBCS_FVTXN * efloat_BBCS_FVTXN / 4 * pow(float_BBCS_FVTXN, 3) ) );
    // float ereso_FVTXN_fn = sqrt( ( efloat_FVTXN_FVTXS * efloat_FVTXN_FVTXS / 4 * float_FVTXN_FVTXS )
    //                              + ( efloat_BBCS_FVTXN * efloat_BBCS_FVTXN / 4 * float_BBCS_FVTXN )
    //                              + ( efloat_BBCS_FVTXS * efloat_BBCS_FVTXS / 4 * pow(float_BBCS_FVTXS, 3) ) );


    // dfb[ic] = Form("%d GeV & %d & %.2e $\\pm$ %.2e & %.2e $\\pm$ %.2e & %.2e $\\pm$ %.2e \\\\",
    //                energy, ic,
    //                reso_BBCS_fn, ereso_BBCS_fn,
    //                reso_FVTXS_fn, ereso_FVTXS_fn,
    //                reso_FVTXN_fn, ereso_FVTXN_fn
    //               );


    // // --- resolution comparisons
    // float rat_BBCS  = reso_BBCS / reso_BBCS_fn;
    // float rat_FVTXS = reso_FVTXS / reso_FVTXS_fn;
    // float rat_FVTXN = reso_FVTXN / reso_FVTXN_fn;

    // drat[ic] = Form("%d GeV & %d & %.2e & %.2e & %.2e \\\\",
    //                 energy, ic,
    //                 rat_BBCS, rat_FVTXS, rat_FVTXN);

    // float avg_BBCS  = 0.5 * (reso_BBCS + reso_BBCS_fn);
    // float avg_FVTXS = 0.5 * (reso_FVTXS + reso_FVTXS_fn);
    // float avg_FVTXN = 0.5 * (reso_FVTXN + reso_FVTXN_fn);

    // davg[ic] = Form("%d GeV & %d & %.2e & %.2e & %.2e \\\\",
    //                 energy, ic,
    //                 avg_BBCS, avg_FVTXS, avg_FVTXN);

    // // ---
    // float reso_FVTXN_xb = sqrt ( ( float_FVTXN_FVTXS * float_CNT_FVTXN ) / float_CNT_FVTXS ) ; // NSNC/SC

    // float nreso_CNT = sqrt((float_CNT_FVTXN * float_CNT_FVTXS) / float_CNT_FVTXS); // BCBS/CS
    // float nreso_FVTXS = sqrt((float_CNT_FVTXS * float_FVTXN_FVTXS) / float_CNT_FVTXN); // CSBS/BC
    // float nreso_FVTXN = sqrt((float_CNT_FVTXN * float_FVTXN_FVTXS) / float_CNT_FVTXS); // CNBN/BC

    // // --- event planes and correlations using CNT
    // TString data4 = Form("%d GeV & %.2e & %.2e & %.2e $\\pm\$ %.2e & %.2e \$\\pm\$ %.2e & %.2e \$\\pm\$ %.2e \\\\",
    //                      energy, nreso_FVTXN, nreso_FVTXS,
    //                      float_FVTXN_FVTXS, efloat_FVTXN_FVTXS,
    //                      float_CNT_FVTXN, efloat_CNT_FVTXN,
    //                      float_CNT_FVTXS, efloat_CNT_FVTXS)

    // //  cout << data4.Data() << endl;
    // // --- very roughly 85%
    // float ratio = nreso_FVTXS / reso_FVTXS;
    // cout << nreso_FVTXS << " " << reso_FVTXS << " " << ratio << endl;

    // float ave = ( nreso_FVTXS + reso_FVTXS_fn ) / 2.0;
    // ratio = ave / reso_FVTXS;
    // cout << nreso_FVTXS << " " << reso_FVTXS << " " << reso_FVTXS_fn << " " << ave << " " << ratio << endl;


  } // ic


  // --- Print everything
  // cout << " means:" << endl;
  // cout << dm << endl;

  // cout << endl;
  // cout << " ep correlations:" << endl;
  // for (int ic = 0; ic < NMUL; ic++)
  //   cout << dcor[ic] << endl;

  cout << endl;
  cout << " resos cnt:" << endl;
  cout << " E & c & R(BBCS) & R(FVTXS) " << endl;
  for (int ic = 0; ic < NMUL; ic++)
    cout << dcnt[ic] << endl;

  cout << endl;
  cout << " resos fvtxs ab:" << endl;
  cout << " E & c & R(BBCS) & R(FVTXS) & R(FVTXSA) & R(FVTXSB)" << endl;
  for (int ic = 0; ic < NMUL; ic++)
    cout << dfvtxsab[ic] << endl;

  cout << endl;
  cout << " resos fvtxs layers:" << endl;
  cout << " E & c & R(FVTXS) & R(FVTXSL0) & R(FVTXSL1) & R(FVTXSL2) & R(FVTXSL3)" << endl;
  for (int ic = 0; ic < NMUL; ic++)
    cout << dfvtxsl[ic] << endl;



  // --- Plot
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  TF1 *fc = new TF1("fc", "[0]", -10, 10);
  fc->SetLineStyle(2);

  int NDRAW = 2;
  // int icdraw[] = {0, 1, 2, 3, 4, 5};
  int icdraw[] = {0, 2, 4};
  int drawMarker[] = {kFullCircle, kFullSquare, kFullDiamond, kFullCross, kFullStar, kOpenCircle};
  int drawSize[] = {1.0, 1.0, 1.5, 1.5, 1.5, 1.0};
  int drawColor[] = {kBlack, kBlue, kRed, kGreen + 2, kOrange + 5, kMagenta + 2};

  int ie = 0;
  if ( energy == 200)
    ie = 0;
  else if ( energy == 62 )
    ie = 1;
  else if ( energy == 39 )
    ie = 2;
  else if ( energy == 20 )
    ie = 3;
  else
    ie = 0;
  const char* centLabel[4][6] =
  {
    {"0-5%", "5-10%", "10-20%", "20-40%", "40-60%", "60-88%"},
    {"0-5%", "5-10%", "10-20%", "20-40%", "40-60%", "60-88%"},
    {"0-10%", "10-20%", "20-40%", "40-60%", "60-88%", ""},
    {"0-20%", "20-40%", "40-60%", "60-88%", "", ""},
  };

  for (int i = 0; i < NDRAW; i++)
  {
    gRBBCS[icdraw[i]]->SetMarkerStyle(drawMarker[i]);
    gRBBCS[icdraw[i]]->SetMarkerSize(drawSize[i]);
    gRBBCS[icdraw[i]]->SetMarkerColor(drawColor[i]);
    gRBBCS[icdraw[i]]->SetLineColor(drawColor[i]);

    gRFVTXS[icdraw[i]]->SetMarkerStyle(drawMarker[i]);
    gRFVTXS[icdraw[i]]->SetMarkerSize(drawSize[i]);
    gRFVTXS[icdraw[i]]->SetMarkerColor(drawColor[i]);
    gRFVTXS[icdraw[i]]->SetLineColor(drawColor[i]);

    for (int il = 0; il < NFVTXLAY; il++)
    {
      gRFVTXSL[icdraw[i]][il]->SetMarkerStyle(drawMarker[i]);
      gRFVTXSL[icdraw[i]][il]->SetMarkerSize(drawSize[i]);
      gRFVTXSL[icdraw[i]][il]->SetMarkerColor(drawColor[i]);
      gRFVTXSL[icdraw[i]][il]->SetLineColor(drawColor[i]);

    }
  }

  TLegend *legbbcs;
  TLegend *legfvtxs;
  TLegend *legfvtxsl[NFVTXLAY];

  // TLegend *legdraw = new TLegend(0.15, 0.8, 0.9, 0.9);
  // legdraw->SetFillStyle(0);
  // legdraw->SetNColumns(NDRAW);
  // for (int i = 0; i < NDRAW; i++)
  //   legdraw->AddEntry(gRBBCS[icdraw[i]], centLabel[ie][icdraw[i]], "P");

  TH1D* haxis_zvrtx = new TH1D("haxis_zvrtx",
                               ";z_{vrtx} [cm];#Psi_{2} resolution",
                               100, -10, 10);
  haxis_zvrtx->GetYaxis()->SetTitleFont(63);
  haxis_zvrtx->GetYaxis()->SetLabelFont(63);
  haxis_zvrtx->GetYaxis()->SetTitleSize(24);
  haxis_zvrtx->GetYaxis()->SetTitleOffset(1.4);
  haxis_zvrtx->GetYaxis()->SetLabelSize(20);
  haxis_zvrtx->GetXaxis()->SetTitleFont(63);
  haxis_zvrtx->GetXaxis()->SetLabelFont(63);
  haxis_zvrtx->GetXaxis()->SetTitleSize(24);
  haxis_zvrtx->GetXaxis()->SetTitleOffset(1.0);
  haxis_zvrtx->GetXaxis()->SetLabelSize(20);

  TLine lreso;
  lreso.SetLineStyle(2);

  TBox breso;
  breso.SetFillStyle(3002);

  TLatex ltitle;
  ltitle.SetNDC();
  ltitle.SetTextAlign(22);

  TCanvas *cab = new TCanvas("cab", "ab", 1000, 1000);
  cab->Divide(1, 2);

  cab->cd(1);
  haxis_zvrtx->GetYaxis()->SetRangeUser(0.0, 0.149);
  haxis_zvrtx->DrawCopy();

  legbbcs = new TLegend(0.10, 0.8, 0.9, 0.9);
  legbbcs->SetFillStyle(0);
  legbbcs->SetNColumns(NDRAW);

  for (int i = 0; i < NDRAW; i++)
  {
    // gRBBCS[icdraw[i]]->Fit(fc, "RQ0N");

    // lreso.SetLineColor(drawColor[i]);
    // lreso.DrawLine(-10, fc->GetParameter(0), 10, fc->GetParameter(0));

    // legbbcs->AddEntry(gRBBCS[icdraw[i]],
    //                   Form("%s #LTR(#Psi_{2})#GT = %.3f #pm %.3f",
    //                        centLabel[ie][icdraw[i]],
    //                        fc->GetParameter(0),
    //                        fc->GetParError(0) ),
    //                        "P");

    lreso.SetLineColor(drawColor[i]);
    lreso.DrawLine(-10, RBBCS[icdraw[i]].first, 10, RBBCS[icdraw[i]].first);

    breso.SetFillColor(drawColor[i]);
    breso.DrawBox(-10, RBBCS[icdraw[i]].first - RBBCS[icdraw[i]].second,
                  +10, RBBCS[icdraw[i]].first + RBBCS[icdraw[i]].second);


    legbbcs->AddEntry(gRBBCS[icdraw[i]],
                      Form("%s #LTR(#Psi_{2})#GT = %.3f #pm %.3f",
                           centLabel[ie][icdraw[i]],
                           RBBCS[icdraw[i]].first,
                           RBBCS[icdraw[i]].second ),
                      "P");

    gRBBCS[icdraw[i]]->Draw("P");
  }
  ltitle.DrawLatex(0.5, 0.95, Form("R(BBCS)  d+Au #sqrt{s_{_{NN}}} = %d GeV", energy));
  legbbcs->Draw("same");

  cab->cd(2);
  haxis_zvrtx->GetYaxis()->SetRangeUser(0.0, 0.34);
  haxis_zvrtx->DrawCopy();

  legfvtxs = new TLegend(0.10, 0.8, 0.9, 0.9);
  legfvtxs->SetFillStyle(0);
  legfvtxs->SetNColumns(NDRAW);

  for (int i = 0; i < NDRAW; i++)
  {
    // gRFVTXS[icdraw[i]]->Fit(fc, "RQ0N");

    // lreso.SetLineColor(drawColor[i]);
    // lreso.DrawLine(-10, fc->GetParameter(0), 10, fc->GetParameter(0));

    // legfvtxs->AddEntry(gRFVTXS[icdraw[i]],
    //                   Form("%s #LTR(#Psi_{2})#GT = %.3f #pm %.3f",
    //                        centLabel[ie][icdraw[i]],
    //                        fc->GetParameter(0),
    //                        fc->GetParError(0) ),
    //                        "P");

    lreso.SetLineColor(drawColor[i]);
    lreso.DrawLine(-10, RFVTXS[icdraw[i]].first, 10, RFVTXS[icdraw[i]].first);

    breso.SetFillColor(drawColor[i]);
    breso.DrawBox(-10, RFVTXS[icdraw[i]].first - RFVTXS[icdraw[i]].second,
                  +10, RFVTXS[icdraw[i]].first + RFVTXS[icdraw[i]].second);


    legfvtxs->AddEntry(gRFVTXS[icdraw[i]],
                       Form("%s #LTR(#Psi_{2})#GT = %.3f #pm %.3f",
                            centLabel[ie][icdraw[i]],
                            RFVTXS[icdraw[i]].first,
                            RFVTXS[icdraw[i]].second ),
                       "P");

    gRFVTXS[icdraw[i]]->Draw("P");
  }

  ltitle.DrawLatex(0.5, 0.95, Form("R(FVTXS)  d+Au #sqrt{s_{_{NN}}} = %d GeV", energy));
  legfvtxs->Draw("same");


  cab->Print(Form("pdfs/calc_epreso_zvrtx_dau%d.pdf", energy));


  TCanvas *cfvtxl = new TCanvas("cfvtxl", "fvtxl", 500, 1000);
  cfvtxl->Divide(1, 4);

  for (int il = 0; il < NFVTXLAY; il++)
  {
    cfvtxl->cd(il + 1);
    haxis_zvrtx->GetYaxis()->SetRangeUser(0.0, 0.30);
    haxis_zvrtx->GetYaxis()->SetTitleSize(16);
    haxis_zvrtx->GetYaxis()->SetTitleOffset(2.0);
    haxis_zvrtx->GetYaxis()->SetLabelSize(12);
    haxis_zvrtx->GetXaxis()->SetTitleSize(16);
    haxis_zvrtx->GetXaxis()->SetTitleOffset(1.5);
    haxis_zvrtx->GetXaxis()->SetLabelSize(12);
    haxis_zvrtx->DrawCopy();

    legfvtxsl[il] = new TLegend(0.10, 0.8, 0.9, 0.9);
    legfvtxsl[il]->SetFillStyle(0);
    legfvtxsl[il]->SetNColumns(NDRAW);

    for (int i = 0; i < NDRAW; i++)
    {
      // gRFVTXS[icdraw[i]]->Fit(fc, "RQ0N");

      // lreso.SetLineColor(drawColor[i]);
      // lreso.DrawLine(-10, fc->GetParameter(0), 10, fc->GetParameter(0));

      // legfvtxsl[il]->AddEntry(gRFVTXS[icdraw[i]],
      //                   Form("%s #LTR(#Psi_{2})#GT = %.3f #pm %.3f",
      //                        centLabel[ie][icdraw[i]],
      //                        fc->GetParameter(0),
      //                        fc->GetParError(0) ),
      //                        "P");

      lreso.SetLineColor(drawColor[i]);
      lreso.DrawLine(-10, RFVTXSL[icdraw[i]][il].first, 10, RFVTXSL[icdraw[i]][il].first);

      breso.SetFillColor(drawColor[i]);
      breso.DrawBox(-10, RFVTXSL[icdraw[i]][il].first - RFVTXSL[icdraw[i]][il].second,
                    +10, RFVTXSL[icdraw[i]][il].first + RFVTXSL[icdraw[i]][il].second);


      legfvtxsl[il]->AddEntry(gRFVTXSL[icdraw[i]][il],
                              Form("%s #LTR(#Psi_{2})#GT = %.3f #pm %.3f",
                                   centLabel[ie][icdraw[i]],
                                   RFVTXSL[icdraw[i]][il].first,
                                   RFVTXSL[icdraw[i]][il].second ),
                              "P");

      gRFVTXSL[icdraw[i]][il]->Draw("P");
    }
    legfvtxsl[il]->Draw("same");

    ltitle.DrawLatex(0.5, 0.95, Form("R(FVTXS-L%i)  d+Au #sqrt{s_{_{NN}}} = %d GeV", il, energy));

    // if ( il == 0 )
    //   legdraw->Draw("same");

  }
  cfvtxl->Print(Form("pdfs/calc_epreso_zvrtx_fvtxl_dau%d.pdf", energy));




  TCanvas *clcent = new TCanvas("clcent", "fvtxsl cent", 1200, 1000);
  clcent->SetMargin(0, 0, 0, 0);
  clcent->Divide(4, NMUL);

  for (int ic = 0; ic < NMUL; ic++)
  {
    for (int il = 0; il < NFVTXLAY; il++)
    {
      int pad = NFVTXLAY * ic + il + 1;

      // cout << " " << ic << " " << il << " " << pad << endl;

      clcent->GetPad(pad)->SetMargin(0.1, 0.02, 0.1, 0.02);
      clcent->GetPad(pad)->SetTicks(1, 1);

      clcent->cd(pad);
      haxis_zvrtx->GetYaxis()->SetRangeUser(0.0, 0.30);
      haxis_zvrtx->GetYaxis()->SetTitleSize(14);
      haxis_zvrtx->GetYaxis()->SetTitleOffset(2.0);
      haxis_zvrtx->GetYaxis()->SetLabelSize(10);
      haxis_zvrtx->GetXaxis()->SetTitleSize(14);
      haxis_zvrtx->GetXaxis()->SetTitleOffset(1.5);
      haxis_zvrtx->GetXaxis()->SetLabelSize(10);
      haxis_zvrtx->DrawCopy();

      lreso.SetLineColor(drawColor[ic]);
      lreso.DrawLine(-10, RFVTXSL[ic][il].first, 10, RFVTXSL[ic][il].first);

      breso.SetFillColor(drawColor[ic]);
      breso.DrawBox(-10, RFVTXSL[ic][il].first - RFVTXSL[ic][il].second,
                    +10, RFVTXSL[ic][il].first + RFVTXSL[ic][il].second);

      gRFVTXSL[ic][il]->SetMarkerStyle(drawMarker[ic]);
      gRFVTXSL[ic][il]->SetMarkerSize(drawSize[ic]);
      gRFVTXSL[ic][il]->SetMarkerColor(drawColor[ic]);
      gRFVTXSL[ic][il]->SetLineColor(drawColor[ic]);
      gRFVTXSL[ic][il]->Draw("P");

      ltitle.DrawLatex(0.5, 0.90,
                       Form("FVTXS-L%i  %d GeV  %s  #LTR(#Psi_{2})#GT = %.3f #pm %.3f",
                            il, energy, centLabel[ie][ic],
                            RFVTXSL[ic][il].first,
                            RFVTXSL[ic][il].second));

    } // il
  } // ic

  clcent->Print(Form("pdfs/calc_epreso_fvtxsl_all_dau%d.pdf", energy));

  // cout << endl;
  // cout << " resos fvtx:" << endl;
  // for (int ic = 0; ic < NMUL; ic++)
  //   cout << dfb[ic] << endl;

  // cout << endl;
  // cout << " resos ratio:" << endl;
  // for (int ic = 0; ic < NMUL; ic++)
  //   cout << drat[ic] << endl;

  // cout << endl;
  // cout << " resos avg:" << endl;
  // for (int ic = 0; ic < NMUL; ic++)
  //   cout << davg[ic] << endl;


}


ValErr calc_epreso(ValErr AB, ValErr AC, ValErr BC)
{

  double reso = sqrt( (AB.first * AC.first) / BC.first );
  double ereso = reso / 2. * sqrt( pow(AB.second / AB.first, 2) +
                                   pow(AC.second / AC.first, 2) +
                                   pow(BC.second / BC.first, 2) );

  return make_pair(reso, ereso);
}


