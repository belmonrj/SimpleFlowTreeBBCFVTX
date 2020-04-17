#include "../RpPar.h"

#include <TFile.h>
#include <TProfile.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLine.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>

#include <iostream>
#include <utility>

using namespace std;

typedef pair<double, double> ValErr;

void doenergy(int, int);
void diagnostic(int, int);
ValErr calc_epreso(ValErr AB, ValErr AC, ValErr BC);
ValErr int_vn(TProfile *pvn, int bl, int bh);


TFile* outfile;

void makeplot_vnEP()
{

  outfile = TFile::Open("histograms_vnEPcent.root", "recreate");

  for ( int i = 2; i < 4; ++i )
  {
    doenergy(200, i);
    doenergy(62, i);
    doenergy(39, i);
    doenergy(20, i);
  // ---
  // diagnostic(200, i);
  // diagnostic(62, i);
  // diagnostic(39, i);
  // diagnostic(20, i);
  }
  //  doenergy(200, 2);
  //  doenergy(62, 2);
  //  doenergy(39, 2);
  //  doenergy(20, 2);

  outfile->Write();
  outfile->Close();

}

void doenergy(int energy, int harmonic)
{

  gStyle->SetOptTitle(1);

  // TCanvas* c1 = new TCanvas("c1", "");

  // TFile* file = TFile::Open(Form("input/combined_%d.root",energy));
  TFile* file = TFile::Open(Form("../output/hist_dAu%d.root", energy));

  // --- Histograms for v2 vs centrality at low & high pT
  TH1D* th1d_vncent_fvtxs_lowpt = new TH1D("th1d_vncent_fvtxs_lowpt", "", NMUL, -0.5, NMUL - 0.5);
  TH1D* th1d_vncent_fvtxs_highpt = new TH1D("th1d_vncent_fvtxs_highpt", "", NMUL, -0.5, NMUL - 0.5);

  // --- graphs of integrated v2 vs centrality
  TGraphErrors *gvn_bbcs_both_cent = new TGraphErrors();
  TGraphErrors *gvn_fvtxs_both_cent = new TGraphErrors();
  TGraphErrors *gvn_fvtxsa_both_cent = new TGraphErrors();
  TGraphErrors *gvn_fvtxsb_both_cent = new TGraphErrors();
  TGraphErrors *gvn_fvtxsl_both_cent[NFVTXLAY];
  for (int il = 0; il < NFVTXLAY; il++)
    gvn_fvtxsl_both_cent[il] = new TGraphErrors();

  TGraphErrors *gvn_bbcs_east_cent = new TGraphErrors();
  TGraphErrors *gvn_fvtxs_east_cent = new TGraphErrors();
  TGraphErrors *gvn_fvtxsa_east_cent = new TGraphErrors();
  TGraphErrors *gvn_fvtxsb_east_cent = new TGraphErrors();
  TGraphErrors *gvn_fvtxsl_east_cent[NFVTXLAY];
  for (int il = 0; il < NFVTXLAY; il++)
    gvn_fvtxsl_east_cent[il] = new TGraphErrors();

  TGraphErrors *gvn_bbcs_west_cent = new TGraphErrors();
  TGraphErrors *gvn_fvtxs_west_cent = new TGraphErrors();
  TGraphErrors *gvn_fvtxsa_west_cent = new TGraphErrors();
  TGraphErrors *gvn_fvtxsb_west_cent = new TGraphErrors();
  TGraphErrors *gvn_fvtxsl_west_cent[NFVTXLAY];
  for (int il = 0; il < NFVTXLAY; il++)
    gvn_fvtxsl_west_cent[il] = new TGraphErrors();



  double lopt = 0.7;
  double hipt = 2.25;

  // --- loop over each centrality bin
  for (int ic = 0; ic < NMUL; ic++)
  {
    cout << endl;
    cout << " --- " << endl;
    cout << " --- E:" << energy << " h:" << harmonic << " c:" << ic << endl;
    cout << " --- " << endl;

    // --- Declare zvrtx summed histograms
    TProfile *tp1f_CNT_BBCS;
    TProfile *tp1f_CNT_FVTXS;
    TProfile *tp1f_CNT_FVTXSA;
    TProfile *tp1f_CNT_FVTXSB;
    TProfile *tp1f_CNT_FVTXSL[NFVTXLAY];

    TProfile *tp1f_BBCS_FVTXS;
    TProfile *tp1f_BBCS_FVTXSA;
    TProfile *tp1f_BBCS_FVTXSB;
    TProfile *tp1f_BBCS_FVTXSL[NFVTXLAY];

    TProfile* hvnpt_bbcs_east;
    TProfile* hvnpt_bbcs_west;
    TProfile* hvnpt_bbcs_both;
    TProfile* hvneta_bbcs_east;
    TProfile* hvneta_bbcs_west;
    TProfile* hvneta_bbcs_both;

    TProfile* hvnpt_fvtxs_east;
    TProfile* hvnpt_fvtxs_west;
    TProfile* hvnpt_fvtxs_both;
    TProfile* hvneta_fvtxs_east;
    TProfile* hvneta_fvtxs_west;
    TProfile* hvneta_fvtxs_both;

    TProfile* hvnpt_fvtxsa_east;
    TProfile* hvnpt_fvtxsa_west;
    TProfile* hvnpt_fvtxsa_both;
    TProfile* hvneta_fvtxsa_east;
    TProfile* hvneta_fvtxsa_west;
    TProfile* hvneta_fvtxsa_both;

    TProfile* hvnpt_fvtxsb_east;
    TProfile* hvnpt_fvtxsb_west;
    TProfile* hvnpt_fvtxsb_both;
    TProfile* hvneta_fvtxsb_east;
    TProfile* hvneta_fvtxsb_west;
    TProfile* hvneta_fvtxsb_both;

    TProfile* hvnpt_fvtxsl_east[NFVTXLAY];
    TProfile* hvnpt_fvtxsl_west[NFVTXLAY];
    TProfile* hvnpt_fvtxsl_both[NFVTXLAY];
    TProfile* hvneta_fvtxsl_east[NFVTXLAY];
    TProfile* hvneta_fvtxsl_west[NFVTXLAY];
    TProfile* hvneta_fvtxsl_both[NFVTXLAY];

    TProfile* hvnpt_fvtxs_east_zres;
    TProfile* hvnpt_fvtxs_west_zres;
    TProfile* hvnpt_fvtxs_both_zres;

    TProfile* hvnpt_bbcs_east_zres;
    TProfile* hvnpt_bbcs_west_zres;
    TProfile* hvnpt_bbcs_both_zres;

    TProfile* hvnpt_fvtxsl_east_zres[NFVTXLAY];
    TProfile* hvnpt_fvtxsl_west_zres[NFVTXLAY];
    TProfile* hvnpt_fvtxsl_both_zres[NFVTXLAY];

    TProfile* hvnpt_bbcs_east_zvtx[NZPS];
    TProfile* hvnpt_bbcs_west_zvtx[NZPS];
    TProfile* hvnpt_bbcs_both_zvtx[NZPS];

    TProfile* hvnpt_fvtxs_east_zvtx[NZPS];
    TProfile* hvnpt_fvtxs_west_zvtx[NZPS];
    TProfile* hvnpt_fvtxs_both_zvtx[NZPS];

    TProfile* hvnpt_fvtxsa_east_zvtx[NZPS];
    TProfile* hvnpt_fvtxsa_west_zvtx[NZPS];
    TProfile* hvnpt_fvtxsa_both_zvtx[NZPS];

    TProfile* hvnpt_fvtxsb_east_zvtx[NZPS];
    TProfile* hvnpt_fvtxsb_west_zvtx[NZPS];
    TProfile* hvnpt_fvtxsb_both_zvtx[NZPS];

    TProfile* hvnpt_fvtxsl_east_zvtx[NZPS][NFVTXLAY];
    TProfile* hvnpt_fvtxsl_west_zvtx[NZPS][NFVTXLAY];
    TProfile* hvnpt_fvtxsl_both_zvtx[NZPS][NFVTXLAY];

    TGraphErrors *gvn_bbcs_zdep = new TGraphErrors();
    TGraphErrors *gvn_bbcs_zdep_sys = new TGraphErrors();

    TGraphErrors *gvn_fvtxs_zdep = new TGraphErrors();
    TGraphErrors *gvn_fvtxs_zdep_sys = new TGraphErrors();


    // ---
    // --- loop over zvrtx
    // ---
    for (int iz = 0; iz < NZPS; iz++)
    {

      // --- get histograms from file
      file->cd();
      TProfile* tp1f_CNT_BBCS_z = (TProfile*)file->Get(Form("tp1f_c%d_z%d_reso%d_BBC_CNT", ic, iz, harmonic));
      TProfile* tp1f_CNT_FVTXS_z = (TProfile*)file->Get(Form("tp1f_c%d_z%d_reso%d_CNT_FVTX", ic, iz, harmonic));
      TProfile* tp1f_CNT_FVTXSA_z = (TProfile*)file->Get(Form("tp1f_c%d_z%d_reso%d_CNT_FVTXSA", ic, iz, harmonic));
      TProfile* tp1f_CNT_FVTXSB_z = (TProfile*)file->Get(Form("tp1f_c%d_z%d_reso%d_CNT_FVTXSB", ic, iz, harmonic));
      TProfile* tp1f_CNT_FVTXSL_z[NFVTXLAY];
      for (int il = 0; il < NFVTXLAY; il++)
        tp1f_CNT_FVTXSL_z[il] = (TProfile*)file->Get(Form("tp1f_c%d_z%d_reso%d_CNT_FVTXSL%d", ic, iz, harmonic, il));

      TProfile* tp1f_BBCS_FVTXS_z = (TProfile*)file->Get(Form("tp1f_c%d_z%d_reso%d_BBC_FVTX", ic, iz, harmonic));
      TProfile* tp1f_BBCS_FVTXSA_z = (TProfile*)file->Get(Form("tp1f_c%d_z%d_reso%d_BBC_FVTXSA", ic, iz, harmonic));
      TProfile* tp1f_BBCS_FVTXSB_z = (TProfile*)file->Get(Form("tp1f_c%d_z%d_reso%d_BBC_FVTXSB", ic, iz, harmonic));
      TProfile* tp1f_BBCS_FVTXSL_z[NFVTXLAY];
      for (int il = 0; il < NFVTXLAY; il++)
        tp1f_BBCS_FVTXSL_z[il] = (TProfile*)file->Get(Form("tp1f_c%d_z%d_reso%d_BBC_FVTXSL%d", ic, iz, harmonic, il));

      TProfile* hvnpt_bbcs_east_z = (TProfile*)file->Get(Form("bbcs_v%d_east_docalib_cent%d_zvtx%d", harmonic, ic, iz));
      TProfile* hvnpt_bbcs_west_z = (TProfile*)file->Get(Form("bbcs_v%d_west_docalib_cent%d_zvtx%d", harmonic, ic, iz));
      TProfile* hvnpt_bbcs_both_z = (TProfile*)file->Get(Form("bbcs_v%d_both_docalib_cent%d_zvtx%d", harmonic, ic, iz));
      TProfile* hvneta_bbcs_east_z = (TProfile*)file->Get(Form("bbcs_v%deta_east_docalib_cent%d_zvtx%d", harmonic, ic, iz));
      TProfile* hvneta_bbcs_west_z = (TProfile*)file->Get(Form("bbcs_v%deta_west_docalib_cent%d_zvtx%d", harmonic, ic, iz));
      TProfile* hvneta_bbcs_both_z = (TProfile*)file->Get(Form("bbcs_v%deta_both_docalib_cent%d_zvtx%d", harmonic, ic, iz));

      TProfile* hvnpt_fvtxs_east_z = (TProfile*)file->Get(Form("fvtxs_v%d_east_docalib_cent%d_zvtx%d", harmonic, ic, iz));
      TProfile* hvnpt_fvtxs_west_z = (TProfile*)file->Get(Form("fvtxs_v%d_west_docalib_cent%d_zvtx%d", harmonic, ic, iz));
      TProfile* hvnpt_fvtxs_both_z = (TProfile*)file->Get(Form("fvtxs_v%d_both_docalib_cent%d_zvtx%d", harmonic, ic, iz));
      TProfile* hvneta_fvtxs_east_z = (TProfile*)file->Get(Form("fvtxs_v%deta_east_docalib_cent%d_zvtx%d", harmonic, ic, iz));
      TProfile* hvneta_fvtxs_west_z = (TProfile*)file->Get(Form("fvtxs_v%deta_west_docalib_cent%d_zvtx%d", harmonic, ic, iz));
      TProfile* hvneta_fvtxs_both_z = (TProfile*)file->Get(Form("fvtxs_v%deta_both_docalib_cent%d_zvtx%d", harmonic, ic, iz));

      TProfile* hvnpt_fvtxsa_east_z = (TProfile*)file->Get(Form("fvtxsa_v%d_east_docalib_cent%d_zvtx%d", harmonic, ic, iz));
      TProfile* hvnpt_fvtxsa_west_z = (TProfile*)file->Get(Form("fvtxsa_v%d_west_docalib_cent%d_zvtx%d", harmonic, ic, iz));
      TProfile* hvnpt_fvtxsa_both_z = (TProfile*)file->Get(Form("fvtxsa_v%d_both_docalib_cent%d_zvtx%d", harmonic, ic, iz));
      // TProfile* hvneta_fvtxsa_east_z = (TProfile*)file->Get(Form("fvtxsa_v%deta_east_docalib_cent%d_zvtx%d", harmonic, ic, iz));
      // TProfile* hvneta_fvtxsa_west_z = (TProfile*)file->Get(Form("fvtxsa_v%deta_west_docalib_cent%d_zvtx%d", harmonic, ic, iz));
      // TProfile* hvneta_fvtxsa_both_z = (TProfile*)file->Get(Form("fvtxsa_v%deta_both_docalib_cent%d_zvtx%d", harmonic, ic, iz));

      TProfile* hvnpt_fvtxsb_east_z = (TProfile*)file->Get(Form("fvtxsb_v%d_east_docalib_cent%d_zvtx%d", harmonic, ic, iz));
      TProfile* hvnpt_fvtxsb_west_z = (TProfile*)file->Get(Form("fvtxsb_v%d_west_docalib_cent%d_zvtx%d", harmonic, ic, iz));
      TProfile* hvnpt_fvtxsb_both_z = (TProfile*)file->Get(Form("fvtxsb_v%d_both_docalib_cent%d_zvtx%d", harmonic, ic, iz));
      // TProfile* hvneta_fvtxsb_east_z = (TProfile*)file->Get(Form("fvtxsb_v%deta_east_docalib_cent%d_zvtx%d", harmonic, ic, iz));
      // TProfile* hvneta_fvtxsb_west_z = (TProfile*)file->Get(Form("fvtxsb_v%deta_west_docalib_cent%d_zvtx%d", harmonic, ic, iz));
      // TProfile* hvneta_fvtxsb_both_z = (TProfile*)file->Get(Form("fvtxsb_v%deta_both_docalib_cent%d_zvtx%d", harmonic, ic, iz));

      TProfile* hvnpt_fvtxsl_east_z[NFVTXLAY];
      TProfile* hvnpt_fvtxsl_west_z[NFVTXLAY];
      TProfile* hvnpt_fvtxsl_both_z[NFVTXLAY];
      // TProfile* hvneta_fvtxsl_east_z[NFVTXLAY];
      // TProfile* hvneta_fvtxsl_west_z[NFVTXLAY];
      // TProfile* hvneta_fvtxsl_both_z[NFVTXLAY];
      for (int il = 0; il < NFVTXLAY; il++)
      {
        hvnpt_fvtxsl_east_z[il] = (TProfile*)file->Get(Form("fvtxsl%d_v%d_east_docalib_cent%d_zvtx%d", il, harmonic, ic, iz));
        hvnpt_fvtxsl_west_z[il] = (TProfile*)file->Get(Form("fvtxsl%d_v%d_west_docalib_cent%d_zvtx%d", il, harmonic, ic, iz));
        hvnpt_fvtxsl_both_z[il] = (TProfile*)file->Get(Form("fvtxsl%d_v%d_both_docalib_cent%d_zvtx%d", il, harmonic, ic, iz));
        // hvneta_fvtxsl_east_z[il] = (TProfile*)file->Get(Form("fvtxsl%d_v%deta_east_docalib_cent%d_zvtx%d", il, harmonic, ic, iz));
        // hvneta_fvtxsl_west_z[il] = (TProfile*)file->Get(Form("fvtxsl%d_v%deta_west_docalib_cent%d_zvtx%d", il, harmonic, ic, iz));
        // hvneta_fvtxsl_both_z[il] = (TProfile*)file->Get(Form("fvtxsl%d_v%deta_both_docalib_cent%d_zvtx%d", il, harmonic, ic, iz));
      } // il

      // --- summ histograms over zvrtx w/o correcting
      if ( iz == 0 )
      {

        tp1f_CNT_BBCS = (TProfile*) tp1f_CNT_BBCS_z->Clone("tp1f_CNT_BBCS");
        tp1f_CNT_FVTXS = (TProfile*) tp1f_CNT_FVTXS_z->Clone("tp1f_CNT_FVTXS");
        tp1f_CNT_FVTXSA = (TProfile*) tp1f_CNT_FVTXSA_z->Clone("tp1f_CNT_FVTXSA");
        tp1f_CNT_FVTXSB = (TProfile*) tp1f_CNT_FVTXSB_z->Clone("tp1f_CNT_FVTXSB");
        for (int il = 0; il < NFVTXLAY; il++)
          tp1f_CNT_FVTXSL[il] = (TProfile*) tp1f_CNT_FVTXSL_z[il]->Clone(Form("tp1f_CNT_FVTXSL%d", il));

        tp1f_BBCS_FVTXS = (TProfile*) tp1f_BBCS_FVTXS_z->Clone("tp1f_BBCS_FVTXS");
        tp1f_BBCS_FVTXSA = (TProfile*) tp1f_BBCS_FVTXSA_z->Clone("tp1f_BBCS_FVTXSA");
        tp1f_BBCS_FVTXSB = (TProfile*) tp1f_BBCS_FVTXSB_z->Clone("tp1f_BBCS_FVTXSB");
        for (int il = 0; il < NFVTXLAY; il++)
          tp1f_BBCS_FVTXSL[il] = (TProfile*) tp1f_BBCS_FVTXSL_z[il]->Clone(Form("tp1f_BBCS_FVTXSL%d", il));


        hvnpt_bbcs_east = (TProfile*) hvnpt_bbcs_east_z->Clone("hvnpt_bbcs_east");
        hvnpt_bbcs_west = (TProfile*) hvnpt_bbcs_west_z->Clone("hvnpt_bbcs_west");
        hvnpt_bbcs_both = (TProfile*) hvnpt_bbcs_both_z->Clone("hvnpt_bbcs_both");
        hvneta_bbcs_east = (TProfile*) hvneta_bbcs_east_z->Clone("hvneta_bbcs_east");
        hvneta_bbcs_west = (TProfile*) hvneta_bbcs_west_z->Clone("hvneta_bbcs_west");
        hvneta_bbcs_both = (TProfile*) hvneta_bbcs_both_z->Clone("hvneta_bbcs_both");

        hvnpt_fvtxs_east = (TProfile*) hvnpt_fvtxs_east_z->Clone("hvnpt_fvtxs_east");
        hvnpt_fvtxs_west = (TProfile*) hvnpt_fvtxs_west_z->Clone("hvnpt_fvtxs_west");
        hvnpt_fvtxs_both = (TProfile*) hvnpt_fvtxs_both_z->Clone("hvnpt_fvtxs_both");
        hvneta_fvtxs_east = (TProfile*) hvneta_fvtxs_east_z->Clone("hvneta_fvtxs_east");
        hvneta_fvtxs_west = (TProfile*) hvneta_fvtxs_west_z->Clone("hvneta_fvtxs_west");
        hvneta_fvtxs_both = (TProfile*) hvneta_fvtxs_both_z->Clone("hvneta_fvtxs_both");

        hvnpt_fvtxsa_east = (TProfile*) hvnpt_fvtxsa_east_z->Clone("hvnpt_fvtxsa_east");
        hvnpt_fvtxsa_west = (TProfile*) hvnpt_fvtxsa_west_z->Clone("hvnpt_fvtxsa_west");
        hvnpt_fvtxsa_both = (TProfile*) hvnpt_fvtxsa_both_z->Clone("hvnpt_fvtxsa_both");
        // hvneta_fvtxsa_east = (TProfile*) hvneta_fvtxsa_east_z->Clone("hvneta_fvtxsa_east");
        // hvneta_fvtxsa_west = (TProfile*) hvneta_fvtxsa_west_z->Clone("hvneta_fvtxsa_west");
        // hvneta_fvtxsa_both = (TProfile*) hvneta_fvtxsa_both_z->Clone("hvneta_fvtxsa_both");

        hvnpt_fvtxsb_east = (TProfile*) hvnpt_fvtxsb_east_z->Clone("hvnpt_fvtxsb_east");
        hvnpt_fvtxsb_west = (TProfile*) hvnpt_fvtxsb_west_z->Clone("hvnpt_fvtxsb_west");
        hvnpt_fvtxsb_both = (TProfile*) hvnpt_fvtxsb_both_z->Clone("hvnpt_fvtxsb_both");
        // hvneta_fvtxsb_east = (TProfile*) hvneta_fvtxsb_east_z->Clone("hvneta_fvtxsb_east");
        // hvneta_fvtxsb_west = (TProfile*) hvneta_fvtxsb_west_z->Clone("hvneta_fvtxsb_west");
        // hvneta_fvtxsb_both = (TProfile*) hvneta_fvtxsb_both_z->Clone("hvneta_fvtxsb_both");

        for (int il = 0; il < NFVTXLAY; il++)
        {
          hvnpt_fvtxsl_east[il] = (TProfile*) hvnpt_fvtxsl_east_z[il]->Clone(Form("hvnpt_fvtxsl%d_east", il));
          hvnpt_fvtxsl_west[il] = (TProfile*) hvnpt_fvtxsl_west_z[il]->Clone(Form("hvnpt_fvtxsl%d_west", il));
          hvnpt_fvtxsl_both[il] = (TProfile*) hvnpt_fvtxsl_both_z[il]->Clone(Form("hvnpt_fvtxsl%d_both", il));
          // hvneta_fvtxsl_east[il] = (TProfile*) hvneta_fvtxsl_east_z[il]->Clone(Form("hvneta_fvtxsl%d_east", il));
          // hvneta_fvtxsl_west[il] = (TProfile*) hvneta_fvtxsl_west_z[il]->Clone(Form("hvneta_fvtxsl%d_west", il));
          // hvneta_fvtxsl_both[il] = (TProfile*) hvneta_fvtxsl_both_z[il]->Clone(Form("hvneta_fvtxsl%d_both", il));
        }

      }
      else
      {
        tp1f_CNT_BBCS->Add(tp1f_CNT_BBCS_z);
        tp1f_CNT_FVTXS->Add(tp1f_CNT_FVTXS_z);
        tp1f_CNT_FVTXSA->Add(tp1f_CNT_FVTXSA_z);
        tp1f_CNT_FVTXSB->Add(tp1f_CNT_FVTXSB_z);
        for (int il = 0; il < NFVTXLAY; il++)
          tp1f_CNT_FVTXSL[il]->Add(tp1f_CNT_FVTXSL_z[il]);

        tp1f_BBCS_FVTXS->Add(tp1f_BBCS_FVTXS_z);
        tp1f_BBCS_FVTXSA->Add(tp1f_BBCS_FVTXSA_z);
        tp1f_BBCS_FVTXSB->Add(tp1f_BBCS_FVTXSB_z);
        for (int il = 0; il < NFVTXLAY; il++)
          tp1f_BBCS_FVTXSL[il]->Add(tp1f_BBCS_FVTXSL_z[il]);

        hvnpt_bbcs_east->Add(hvnpt_bbcs_east_z);
        hvnpt_bbcs_west->Add(hvnpt_bbcs_west_z);
        hvnpt_bbcs_both->Add(hvnpt_bbcs_both_z);
        hvneta_bbcs_east->Add(hvneta_bbcs_east_z);
        hvneta_bbcs_west->Add(hvneta_bbcs_west_z);
        hvneta_bbcs_both->Add(hvneta_bbcs_both_z);

        hvnpt_fvtxs_east->Add(hvnpt_fvtxs_east_z);
        hvnpt_fvtxs_west->Add(hvnpt_fvtxs_west_z);
        hvnpt_fvtxs_both->Add(hvnpt_fvtxs_both_z);
        hvneta_fvtxs_east->Add(hvneta_fvtxs_east_z);
        hvneta_fvtxs_west->Add(hvneta_fvtxs_west_z);
        hvneta_fvtxs_both->Add(hvneta_fvtxs_both_z);

        hvnpt_fvtxsa_east->Add(hvnpt_fvtxsa_east_z);
        hvnpt_fvtxsa_west->Add(hvnpt_fvtxsa_west_z);
        hvnpt_fvtxsa_both->Add(hvnpt_fvtxsa_both_z);
        // hvneta_fvtxsa_east->Add(hvneta_fvtxsa_east_z);
        // hvneta_fvtxsa_west->Add(hvneta_fvtxsa_west_z);
        // hvneta_fvtxsa_both->Add(hvneta_fvtxsa_both_z);

        hvnpt_fvtxsb_east->Add(hvnpt_fvtxsb_east_z);
        hvnpt_fvtxsb_west->Add(hvnpt_fvtxsb_west_z);
        hvnpt_fvtxsb_both->Add(hvnpt_fvtxsb_both_z);
        // hvneta_fvtxsb_east->Add(hvneta_fvtxsb_east_z);
        // hvneta_fvtxsb_west->Add(hvneta_fvtxsb_west_z);
        // hvneta_fvtxsb_both->Add(hvneta_fvtxsb_both_z);

        for (int il = 0; il < NFVTXLAY; il++)
        {
          hvnpt_fvtxsl_east[il]->Add(hvnpt_fvtxsl_east_z[il]);
          hvnpt_fvtxsl_west[il]->Add(hvnpt_fvtxsl_west_z[il]);
          hvnpt_fvtxsl_both[il]->Add(hvnpt_fvtxsl_both_z[il]);
          // hvneta_fvtxsl_east[il]->Add(hvneta_fvtxsl_east_z[il]);
          // hvneta_fvtxsl_west[il]->Add(hvneta_fvtxsl_west_z[il]);
          // hvneta_fvtxsl_both[il]->Add(hvneta_fvtxsl_both_z[il]);
        } // il

      } // else

      // --- calculate ep resolution
      ValErr CNT_BBCS_z = make_pair(tp1f_CNT_BBCS_z->GetBinContent(1),
                                    tp1f_CNT_BBCS_z->GetBinError(1));
      ValErr CNT_FVTXS_z = make_pair(tp1f_CNT_FVTXS_z->GetBinContent(1),
                                     tp1f_CNT_FVTXS_z->GetBinError(1));
      ValErr CNT_FVTXSL_z[NFVTXLAY];
      for (int il = 0; il < NFVTXLAY; il++)
        CNT_FVTXSL_z[il] = make_pair(tp1f_CNT_FVTXSL_z[il]->GetBinContent(1),
                                     tp1f_CNT_FVTXSL_z[il]->GetBinError(1));

      ValErr BBCS_FVTXS_z = make_pair(tp1f_BBCS_FVTXS_z->GetBinContent(1),
                                      tp1f_BBCS_FVTXS_z->GetBinError(1));
      ValErr BBCS_FVTXSL_z[NFVTXLAY];
      for (int il = 0; il < NFVTXLAY; il++)
        BBCS_FVTXSL_z[il] = make_pair(tp1f_BBCS_FVTXSL_z[il]->GetBinContent(1),
                                      tp1f_BBCS_FVTXSL_z[il]->GetBinError(1));


      ValErr reso_BBCS_z = (calc_epreso(CNT_BBCS_z, BBCS_FVTXS_z, CNT_FVTXS_z));
      ValErr reso_FVTXS_z = (calc_epreso(CNT_FVTXS_z, BBCS_FVTXS_z, CNT_BBCS_z));
      ValErr reso_FVTXSL_z[NFVTXLAY];
      for (int il = 0; il < NFVTXLAY; il++)
        reso_FVTXSL_z[il] = (calc_epreso(CNT_FVTXSL_z[il], BBCS_FVTXSL_z[il], CNT_BBCS_z));

      // --- get the pt integrated result (before resolution correction)
      float z = -10. + (iz + 0.5) * 20. / (float)NZPS;

      int bl = hvnpt_bbcs_both_z->GetXaxis()->FindBin(0.5);
      int bh = hvnpt_bbcs_both_z->GetXaxis()->FindBin(2.0) - 1;

      // bbcs
      ValErr bbcs_int = int_vn(hvnpt_bbcs_both_z, bl, bh);
      float bbcs_sys = reso_BBCS_z.second / reso_BBCS_z.first;

      // correct for resolution
      bbcs_int.first /= reso_BBCS_z.first;
      bbcs_int.second /= reso_BBCS_z.first;

      // add resolution uncertainty in quadriture
      bbcs_int.second = bbcs_int.first * sqrt(pow(bbcs_int.second / bbcs_int.first, 2) + pow(bbcs_sys, 2));

      gvn_bbcs_zdep->SetPoint(iz, z - 0.1, bbcs_int.first);
      gvn_bbcs_zdep->SetPointError(iz, 0, bbcs_int.second);

      // fvtxs
      ValErr fvtxs_int = int_vn(hvnpt_fvtxs_both_z, bl, bh);
      float fvtxs_sys = reso_FVTXS_z.second / reso_FVTXS_z.first;

      fvtxs_int.first /= reso_FVTXS_z.first;
      fvtxs_int.second /= reso_FVTXS_z.first;

      fvtxs_int.second = fvtxs_int.first * sqrt(pow(fvtxs_int.second / fvtxs_int.first, 2) + pow(fvtxs_sys, 2));

      gvn_fvtxs_zdep->SetPoint(iz, z - 0.1, fvtxs_int.first);
      gvn_fvtxs_zdep->SetPointError(iz, 0, fvtxs_int.second);

      // --- correct ep
      hvnpt_bbcs_east_z->Scale(1. / reso_BBCS_z.first);
      hvnpt_bbcs_west_z->Scale(1. / reso_BBCS_z.first);
      hvnpt_bbcs_both_z->Scale(1. / reso_BBCS_z.first);

      hvnpt_fvtxs_east_z->Scale(1. / reso_FVTXS_z.first);
      hvnpt_fvtxs_west_z->Scale(1. / reso_FVTXS_z.first);
      hvnpt_fvtxs_both_z->Scale(1. / reso_FVTXS_z.first);

      for (int il = 0; il < NFVTXLAY; il++)
      {
        hvnpt_fvtxsl_east_z[il]->Scale(1. / reso_FVTXSL_z[il].first);
        hvnpt_fvtxsl_west_z[il]->Scale(1. / reso_FVTXSL_z[il].first);
        hvnpt_fvtxsl_both_z[il]->Scale(1. / reso_FVTXSL_z[il].first);
      }

      // --- sum corrections
      if ( iz == 0 )
      {
        hvnpt_bbcs_east_zres = (TProfile*) hvnpt_bbcs_east_z->Clone("hvnpt_bbcs_east_zres");
        hvnpt_bbcs_west_zres = (TProfile*) hvnpt_bbcs_west_z->Clone("hvnpt_bbcs_west_zres");
        hvnpt_bbcs_both_zres = (TProfile*) hvnpt_bbcs_both_z->Clone("hvnpt_bbcs_both_zres");

        hvnpt_fvtxs_east_zres = (TProfile*) hvnpt_fvtxs_east_z->Clone("hvnpt_fvtxs_east_zres");
        hvnpt_fvtxs_west_zres = (TProfile*) hvnpt_fvtxs_west_z->Clone("hvnpt_fvtxs_west_zres");
        hvnpt_fvtxs_both_zres = (TProfile*) hvnpt_fvtxs_both_z->Clone("hvnpt_fvtxs_both_zres");

        for (int il = 0; il < NFVTXLAY; il++)
        {
          hvnpt_fvtxsl_east_zres[il] = (TProfile*) hvnpt_fvtxsl_east_z[il]->Clone(Form("hvnpt_fvtxsl%d_east_zres", il));
          hvnpt_fvtxsl_west_zres[il] = (TProfile*) hvnpt_fvtxsl_west_z[il]->Clone(Form("hvnpt_fvtxsl%d_west_zres", il));
          hvnpt_fvtxsl_both_zres[il] = (TProfile*) hvnpt_fvtxsl_both_z[il]->Clone(Form("hvnpt_fvtxsl%d_both_zres", il));
        }
      }
      else
      {
        hvnpt_bbcs_east_zres->Add(hvnpt_bbcs_east_z);
        hvnpt_bbcs_west_zres->Add(hvnpt_bbcs_west_z);
        hvnpt_bbcs_both_zres->Add(hvnpt_bbcs_both_z);

        hvnpt_fvtxs_east_zres->Add(hvnpt_fvtxs_east_z);
        hvnpt_fvtxs_west_zres->Add(hvnpt_fvtxs_west_z);
        hvnpt_fvtxs_both_zres->Add(hvnpt_fvtxs_both_z);

        for (int il = 0; il < NFVTXLAY; il++)
        {
          hvnpt_fvtxsl_east_zres[il]->Add(hvnpt_fvtxsl_east_z[il]);
          hvnpt_fvtxsl_west_zres[il]->Add(hvnpt_fvtxsl_west_z[il]);
          hvnpt_fvtxsl_both_zres[il]->Add(hvnpt_fvtxsl_both_z[il]);
        }
      }


      // --- save the zvrtx dependent v2's after resolution correction
      //     note I should have just done this above, this is really a hack

      hvnpt_bbcs_east_zvtx[iz] = (TProfile*) hvnpt_bbcs_east_z->Clone(Form("hvnpt_bbcs_east_z%i", iz));
      hvnpt_bbcs_west_zvtx[iz] = (TProfile*) hvnpt_bbcs_west_z->Clone(Form("hvnpt_bbcs_west_z%i", iz));
      hvnpt_bbcs_both_zvtx[iz] = (TProfile*) hvnpt_bbcs_both_z->Clone(Form("hvnpt_bbcs_both_z%i", iz));

      hvnpt_fvtxs_east_zvtx[iz] = (TProfile*) hvnpt_fvtxs_east_z->Clone(Form("hvnpt_fvtxs_east_z%i", iz));
      hvnpt_fvtxs_west_zvtx[iz] = (TProfile*) hvnpt_fvtxs_west_z->Clone(Form("hvnpt_fvtxs_west_z%i", iz));
      hvnpt_fvtxs_both_zvtx[iz] = (TProfile*) hvnpt_fvtxs_both_z->Clone(Form("hvnpt_fvtxs_both_z%i", iz));

      for (int il = 0; il < NFVTXLAY; il++)
      {
        hvnpt_fvtxsl_east_zvtx[iz][il] = (TProfile*) hvnpt_fvtxsl_east_z[il]->Clone(Form("hvnpt_fvtxsl%i_east_z%i", il, iz));
        hvnpt_fvtxsl_west_zvtx[iz][il] = (TProfile*) hvnpt_fvtxsl_west_z[il]->Clone(Form("hvnpt_fvtxsl%i_west_z%i", il, iz));
        hvnpt_fvtxsl_both_zvtx[iz][il] = (TProfile*) hvnpt_fvtxsl_both_z[il]->Clone(Form("hvnpt_fvtxsl%i_both_z%i", il, iz));
      } // il


    } // iz

    // ---
    // --- Calculate EP resolutions
    // ---

    ValErr CNT_BBCS = make_pair(tp1f_CNT_BBCS->GetBinContent(1),
                                tp1f_CNT_BBCS->GetBinError(1));
    ValErr CNT_FVTXS = make_pair(tp1f_CNT_FVTXS->GetBinContent(1),
                                 tp1f_CNT_FVTXS->GetBinError(1));
    ValErr CNT_FVTXSA = make_pair(tp1f_CNT_FVTXSA->GetBinContent(1),
                                  tp1f_CNT_FVTXSA->GetBinError(1));
    ValErr CNT_FVTXSB = make_pair(tp1f_CNT_FVTXSB->GetBinContent(1),
                                  tp1f_CNT_FVTXSB->GetBinError(1));
    ValErr CNT_FVTXSL[NFVTXLAY];
    for (int il = 0; il < NFVTXLAY; il++)
      CNT_FVTXSL[il] = make_pair(tp1f_CNT_FVTXSL[il]->GetBinContent(1),
                                 tp1f_CNT_FVTXSL[il]->GetBinError(1));

    ValErr BBCS_FVTXS = make_pair(tp1f_BBCS_FVTXS->GetBinContent(1),
                                  tp1f_BBCS_FVTXS->GetBinError(1));
    ValErr BBCS_FVTXSA = make_pair(tp1f_BBCS_FVTXSA->GetBinContent(1),
                                   tp1f_BBCS_FVTXSA->GetBinError(1));
    ValErr BBCS_FVTXSB = make_pair(tp1f_BBCS_FVTXSB->GetBinContent(1),
                                   tp1f_BBCS_FVTXSB->GetBinError(1));
    ValErr BBCS_FVTXSL[NFVTXLAY];
    for (int il = 0; il < NFVTXLAY; il++)
      BBCS_FVTXSL[il] = make_pair(tp1f_BBCS_FVTXSL[il]->GetBinContent(1),
                                  tp1f_BBCS_FVTXSL[il]->GetBinError(1));


    ValErr reso_BBCS = calc_epreso(CNT_BBCS, BBCS_FVTXS, CNT_FVTXS);
    ValErr reso_FVTXS = calc_epreso(CNT_FVTXS, BBCS_FVTXS, CNT_BBCS);
    ValErr reso_FVTXSA = calc_epreso(CNT_FVTXSA, BBCS_FVTXSA, CNT_BBCS);
    ValErr reso_FVTXSB = calc_epreso(CNT_FVTXSB, BBCS_FVTXSB, CNT_BBCS);
    ValErr reso_FVTXSL[NFVTXLAY];
    for (int il = 0; il < NFVTXLAY; il++)
      reso_FVTXSL[il] = calc_epreso(CNT_FVTXSL[il], BBCS_FVTXSL[il], CNT_BBCS);


    if ( energy == 20 )
    {
      // --- i like chocolate and fudge and other tasty deserts
      cout << "Now making empirical adjustment for 20 GeV data" << endl;
      // reso_BBCS = 0.015;
      // reso_FVTXS = 0.04;
      // don't resolution correct the 20 GeV, so that we can rebin it later
      reso_BBCS.first = 1;
      reso_FVTXS.first = 1;
    }

    cout << "bbc resolution is     " << reso_BBCS.first << endl;
    cout << "fvtx resolution is    " << reso_FVTXS.first << endl;
    cout << "fvtxsA resolution is  " << reso_FVTXSA.first << endl;
    cout << "fvtxsB resolution is  " << reso_FVTXSB.first << endl;
    for (int il = 0; il < NFVTXLAY; il++)
      cout << "fvtxsL" << il << " resolution is " << reso_FVTXSL[il].first << endl;



    // --- get the pt integrated result (before resolution correction)
    int bl = hvnpt_bbcs_both->GetXaxis()->FindBin(0.5);
    int bh = hvnpt_bbcs_both->GetXaxis()->FindBin(2.0) - 1;

    // bbcs
    ValErr bbcs_int = int_vn(hvnpt_bbcs_both, bl, bh);
    float bbcs_sys = reso_BBCS.second / reso_BBCS.first;
    // correct for resolution
    bbcs_int.first /= reso_BBCS.first;
    bbcs_int.second /= reso_BBCS.first;
    // add resolution uncertainty in quadriture
    bbcs_int.second = bbcs_int.first * sqrt(pow(bbcs_int.second / bbcs_int.first, 2) + pow(bbcs_sys, 2));
    // fill tgraph
    gvn_bbcs_both_cent->SetPoint(ic, ic, bbcs_int.first);
    gvn_bbcs_both_cent->SetPointError(ic, 0, bbcs_int.second);

    bbcs_int = int_vn(hvnpt_bbcs_east, bl, bh);
    bbcs_sys = reso_BBCS.second / reso_BBCS.first;
    // correct for resolution
    bbcs_int.first /= reso_BBCS.first;
    bbcs_int.second /= reso_BBCS.first;
    // add resolution uncertainty in quadriture
    bbcs_int.second = bbcs_int.first * sqrt(pow(bbcs_int.second / bbcs_int.first, 2) + pow(bbcs_sys, 2));
    // fill tgraph
    gvn_bbcs_east_cent->SetPoint(ic, ic, bbcs_int.first);
    gvn_bbcs_east_cent->SetPointError(ic, 0, bbcs_int.second);

    bbcs_int = int_vn(hvnpt_bbcs_west, bl, bh);
    bbcs_sys = reso_BBCS.second / reso_BBCS.first;
    // correct for resolution
    bbcs_int.first /= reso_BBCS.first;
    bbcs_int.second /= reso_BBCS.first;
    // add resolution uncertainty in quadriture
    bbcs_int.second = bbcs_int.first * sqrt(pow(bbcs_int.second / bbcs_int.first, 2) + pow(bbcs_sys, 2));
    // fill tgraph
    gvn_bbcs_west_cent->SetPoint(ic, ic, bbcs_int.first);
    gvn_bbcs_west_cent->SetPointError(ic, 0, bbcs_int.second);


    // fvtxs
    ValErr fvtxs_int = int_vn(hvnpt_fvtxs_both, bl, bh);
    float fvtxs_sys = reso_FVTXS.second / reso_FVTXS.first;
    // correct for resolution
    fvtxs_int.first /= reso_FVTXS.first;
    fvtxs_int.second /= reso_FVTXS.first;
    // add resolution uncertainty in quadriture
    fvtxs_int.second = fvtxs_int.first * sqrt(pow(fvtxs_int.second / fvtxs_int.first, 2) + pow(fvtxs_sys, 2));
    // fill tgraph
    gvn_fvtxs_both_cent->SetPoint(ic, ic, fvtxs_int.first);
    gvn_fvtxs_both_cent->SetPointError(ic, 0, fvtxs_int.second);

    fvtxs_int = int_vn(hvnpt_fvtxs_east, bl, bh);
    fvtxs_sys = reso_FVTXS.second / reso_FVTXS.first;
    // correct for resolution
    fvtxs_int.first /= reso_FVTXS.first;
    fvtxs_int.second /= reso_FVTXS.first;
    // add resolution uncertainty in quadriture
    fvtxs_int.second = fvtxs_int.first * sqrt(pow(fvtxs_int.second / fvtxs_int.first, 2) + pow(fvtxs_sys, 2));
    // fill tgraph
    gvn_fvtxs_east_cent->SetPoint(ic, ic, fvtxs_int.first);
    gvn_fvtxs_east_cent->SetPointError(ic, 0, fvtxs_int.second);

    fvtxs_int = int_vn(hvnpt_fvtxs_west, bl, bh);
    fvtxs_sys = reso_FVTXS.second / reso_FVTXS.first;
    // correct for resolution
    fvtxs_int.first /= reso_FVTXS.first;
    fvtxs_int.second /= reso_FVTXS.first;
    // add resolution uncertainty in quadriture
    fvtxs_int.second = fvtxs_int.first * sqrt(pow(fvtxs_int.second / fvtxs_int.first, 2) + pow(fvtxs_sys, 2));
    // fill tgraph
    gvn_fvtxs_west_cent->SetPoint(ic, ic, fvtxs_int.first);
    gvn_fvtxs_west_cent->SetPointError(ic, 0, fvtxs_int.second);


    // fvtxsa
    ValErr fvtxsa_int = int_vn(hvnpt_fvtxsa_both, bl, bh);
    float fvtxsa_sys = reso_FVTXSA.second / reso_FVTXSA.first;
    // correct for resolution
    fvtxsa_int.first /= reso_FVTXSA.first;
    fvtxsa_int.second /= reso_FVTXSA.first;
    // add resolution uncertainty in quadriture
    fvtxsa_int.second = fvtxsa_int.first * sqrt(pow(fvtxsa_int.second / fvtxsa_int.first, 2) + pow(fvtxsa_sys, 2));
    // fill tgraph
    gvn_fvtxsa_both_cent->SetPoint(ic, ic, fvtxsa_int.first);
    gvn_fvtxsa_both_cent->SetPointError(ic, 0, fvtxsa_int.second);

    fvtxsa_int = int_vn(hvnpt_fvtxsa_east, bl, bh);
    fvtxsa_sys = reso_FVTXSA.second / reso_FVTXSA.first;
    // correct for resolution
    fvtxsa_int.first /= reso_FVTXSA.first;
    fvtxsa_int.second /= reso_FVTXSA.first;
    // add resolution uncertainty in quadriture
    fvtxsa_int.second = fvtxsa_int.first * sqrt(pow(fvtxsa_int.second / fvtxsa_int.first, 2) + pow(fvtxsa_sys, 2));
    // fill tgraph
    gvn_fvtxsa_east_cent->SetPoint(ic, ic, fvtxsa_int.first);
    gvn_fvtxsa_east_cent->SetPointError(ic, 0, fvtxsa_int.second);

    fvtxsa_int = int_vn(hvnpt_fvtxsa_west, bl, bh);
    fvtxsa_sys = reso_FVTXSA.second / reso_FVTXSA.first;
    // correct for resolution
    fvtxsa_int.first /= reso_FVTXSA.first;
    fvtxsa_int.second /= reso_FVTXSA.first;
    // add resolution uncertainty in quadriture
    fvtxsa_int.second = fvtxsa_int.first * sqrt(pow(fvtxsa_int.second / fvtxsa_int.first, 2) + pow(fvtxsa_sys, 2));
    // fill tgraph
    gvn_fvtxsa_west_cent->SetPoint(ic, ic, fvtxsa_int.first);
    gvn_fvtxsa_west_cent->SetPointError(ic, 0, fvtxsa_int.second);


    // fvtxsb
    ValErr fvtxsb_int = int_vn(hvnpt_fvtxsb_both, bl, bh);
    float fvtxsb_sys = reso_FVTXSB.second / reso_FVTXSB.first;
    // correct for resolution
    fvtxsb_int.first /= reso_FVTXSB.first;
    fvtxsb_int.second /= reso_FVTXSB.first;
    // add resolution uncertainty in quadriture
    fvtxsb_int.second = fvtxsb_int.first * sqrt(pow(fvtxsb_int.second / fvtxsb_int.first, 2) + pow(fvtxsb_sys, 2));
    // fill tgraph
    gvn_fvtxsb_both_cent->SetPoint(ic, ic, fvtxsb_int.first);
    gvn_fvtxsb_both_cent->SetPointError(ic, 0, fvtxsb_int.second);

    fvtxsb_int = int_vn(hvnpt_fvtxsb_east, bl, bh);
    fvtxsb_sys = reso_FVTXSB.second / reso_FVTXSB.first;
    // correct for resolution
    fvtxsb_int.first /= reso_FVTXSB.first;
    fvtxsb_int.second /= reso_FVTXSB.first;
    // add resolution uncertainty in quadriture
    fvtxsb_int.second = fvtxsb_int.first * sqrt(pow(fvtxsb_int.second / fvtxsb_int.first, 2) + pow(fvtxsb_sys, 2));
    // fill tgraph
    gvn_fvtxsb_east_cent->SetPoint(ic, ic, fvtxsb_int.first);
    gvn_fvtxsb_east_cent->SetPointError(ic, 0, fvtxsb_int.second);

    fvtxsb_int = int_vn(hvnpt_fvtxsb_west, bl, bh);
    fvtxsb_sys = reso_FVTXSB.second / reso_FVTXSB.first;
    // correct for resolution
    fvtxsb_int.first /= reso_FVTXSB.first;
    fvtxsb_int.second /= reso_FVTXSB.first;
    // add resolution uncertainty in quadriture
    fvtxsb_int.second = fvtxsb_int.first * sqrt(pow(fvtxsb_int.second / fvtxsb_int.first, 2) + pow(fvtxsb_sys, 2));
    // fill tgraph
    gvn_fvtxsb_west_cent->SetPoint(ic, ic, fvtxsb_int.first);
    gvn_fvtxsb_west_cent->SetPointError(ic, 0, fvtxsb_int.second);


    // fvtxsl
    for (int il = 0; il < NFVTXLAY; il++)
    {
      ValErr fvtxsl_int = int_vn(hvnpt_fvtxsl_both[il], bl, bh);
      float fvtxsl_sys = reso_FVTXSL[il].second / reso_FVTXSL[il].first;
      // correct for resolution
      fvtxsl_int.first /= reso_FVTXSL[il].first;
      fvtxsl_int.second /= reso_FVTXSL[il].first;
      // add resolution uncertainty in quadriture
      fvtxsl_int.second = fvtxsl_int.first * sqrt(pow(fvtxsl_int.second / fvtxsl_int.first, 2) + pow(fvtxsl_sys, 2));
      // fill tgraph
      gvn_fvtxsl_both_cent[il]->SetPoint(ic, ic, fvtxsl_int.first);
      gvn_fvtxsl_both_cent[il]->SetPointError(ic, 0, fvtxsl_int.second);

      fvtxsl_int = int_vn(hvnpt_fvtxsl_east[il], bl, bh);
      fvtxsl_sys = reso_FVTXSL[il].second / reso_FVTXSL[il].first;
      // correct for resolution
      fvtxsl_int.first /= reso_FVTXSL[il].first;
      fvtxsl_int.second /= reso_FVTXSL[il].first;
      // add resolution uncertainty in quadriture
      fvtxsl_int.second = fvtxsl_int.first * sqrt(pow(fvtxsl_int.second / fvtxsl_int.first, 2) + pow(fvtxsl_sys, 2));
      // fill tgraph
      gvn_fvtxsl_east_cent[il]->SetPoint(ic, ic, fvtxsl_int.first);
      gvn_fvtxsl_east_cent[il]->SetPointError(ic, 0, fvtxsl_int.second);

      fvtxsl_int = int_vn(hvnpt_fvtxsl_west[il], bl, bh);
      fvtxsl_sys = reso_FVTXSL[il].second / reso_FVTXSL[il].first;
      // correct for resolution
      fvtxsl_int.first /= reso_FVTXSL[il].first;
      fvtxsl_int.second /= reso_FVTXSL[il].first;
      // add resolution uncertainty in quadriture
      fvtxsl_int.second = fvtxsl_int.first * sqrt(pow(fvtxsl_int.second / fvtxsl_int.first, 2) + pow(fvtxsl_sys, 2));
      // fill tgraph
      gvn_fvtxsl_west_cent[il]->SetPoint(ic, ic, fvtxsl_int.first);
      gvn_fvtxsl_west_cent[il]->SetPointError(ic, 0, fvtxsl_int.second);
    } // il




    // --- divide by resolution
    hvnpt_bbcs_east->Scale(1.0 / reso_BBCS.first);
    hvnpt_bbcs_west->Scale(1.0 / reso_BBCS.first);
    hvnpt_bbcs_both->Scale(1.0 / reso_BBCS.first);
    hvneta_bbcs_east->Scale(1.0 / reso_BBCS.first);
    hvneta_bbcs_west->Scale(1.0 / reso_BBCS.first);
    hvneta_bbcs_both->Scale(1.0 / reso_BBCS.first);

    hvnpt_fvtxs_east->Scale(1.0 / reso_FVTXS.first);
    hvnpt_fvtxs_west->Scale(1.0 / reso_FVTXS.first);
    hvnpt_fvtxs_both->Scale(1.0 / reso_FVTXS.first);
    hvneta_fvtxs_east->Scale(1.0 / reso_FVTXS.first);
    hvneta_fvtxs_west->Scale(1.0 / reso_FVTXS.first);
    hvneta_fvtxs_both->Scale(1.0 / reso_FVTXS.first);

    hvnpt_fvtxsa_east->Scale(1.0 / reso_FVTXSA.first);
    hvnpt_fvtxsa_west->Scale(1.0 / reso_FVTXSA.first);
    hvnpt_fvtxsa_both->Scale(1.0 / reso_FVTXSA.first);
    // hvneta_fvtxsa_east->Scale(1.0 / reso_FVTXSA.first);
    // hvneta_fvtxsa_west->Scale(1.0 / reso_FVTXSA.first);
    // hvneta_fvtxsa_both->Scale(1.0 / reso_FVTXSA.first);

    hvnpt_fvtxsb_east->Scale(1.0 / reso_FVTXSB.first);
    hvnpt_fvtxsb_west->Scale(1.0 / reso_FVTXSB.first);
    hvnpt_fvtxsb_both->Scale(1.0 / reso_FVTXSB.first);
    // hvneta_fvtxsb_east->Scale(1.0 / reso_FVTXSB.first);
    // hvneta_fvtxsb_west->Scale(1.0 / reso_FVTXSB.first);
    // hvneta_fvtxsb_both->Scale(1.0 / reso_FVTXSB.first);

    for (int il = 0; il < NFVTXLAY; il++)
    {
      hvnpt_fvtxsl_east[il]->Scale(1.0 / reso_FVTXSL[il].first);
      hvnpt_fvtxsl_west[il]->Scale(1.0 / reso_FVTXSL[il].first);
      hvnpt_fvtxsl_both[il]->Scale(1.0 / reso_FVTXSL[il].first);
      // hvneta_fvtxsl_east[il]->Scale(1.0 / reso_FVTXSL[il].first);
      // hvneta_fvtxsl_west[il]->Scale(1.0 / reso_FVTXSL[il].first);
      // hvneta_fvtxsl_both[il]->Scale(1.0 / reso_FVTXSL[il].first);
    }


    float blpt = hvnpt_fvtxs_both->FindBin(lopt);
    th1d_vncent_fvtxs_lowpt->SetBinContent(ic + 1, hvnpt_fvtxs_both->GetBinContent(blpt));
    th1d_vncent_fvtxs_lowpt->SetBinError(ic + 1, hvnpt_fvtxs_both->GetBinError(blpt));

    float bhpt = hvnpt_fvtxs_both->FindBin(hipt);
    th1d_vncent_fvtxs_highpt->SetBinContent(ic + 1, hvnpt_fvtxs_both->GetBinContent(bhpt));
    th1d_vncent_fvtxs_highpt->SetBinError(ic + 1, hvnpt_fvtxs_both->GetBinError(bhpt));


    // --- write to file
    outfile->cd();

    hvnpt_bbcs_both->SetName(Form("tprofile_v%d_pT_eventplane_bbcs_c%d_%d", harmonic, ic, energy));
    hvnpt_bbcs_east->SetName(Form("tprofile_v%d_pT_east_eventplane_bbcs_c%d_%d", harmonic, ic, energy));
    hvnpt_bbcs_west->SetName(Form("tprofile_v%d_pT_west_eventplane_bbcs_c%d_%d", harmonic, ic, energy));
    hvneta_bbcs_both->SetName(Form("tprofile_v%d_eta_eventplane_bbcs_c%d_%d", harmonic, ic, energy));
    hvneta_bbcs_east->SetName(Form("tprofile_v%d_eta_east_eventplane_bbcs_c%d_%d", harmonic, ic, energy));
    hvneta_bbcs_west->SetName(Form("tprofile_v%d_eta_west_eventplane_bbcs_c%d_%d", harmonic, ic, energy));

    hvnpt_fvtxs_both->SetName(Form("tprofile_v%d_pT_eventplane_fvtxs_c%d_%d", harmonic, ic, energy));
    hvnpt_fvtxs_east->SetName(Form("tprofile_v%d_pT_east_eventplane_fvtxs_c%d_%d", harmonic, ic, energy));
    hvnpt_fvtxs_west->SetName(Form("tprofile_v%d_pT_west_eventplane_fvtxs_c%d_%d", harmonic, ic, energy));
    hvneta_fvtxs_both->SetName(Form("tprofile_v%d_eta_eventplane_fvtxs_c%d_%d", harmonic, ic, energy));
    hvneta_fvtxs_east->SetName(Form("tprofile_v%d_eta_east_eventplane_fvtxs_c%d_%d", harmonic, ic, energy));
    hvneta_fvtxs_west->SetName(Form("tprofile_v%d_eta_west_eventplane_fvtxs_c%d_%d", harmonic, ic, energy));

    hvnpt_fvtxsa_both->SetName(Form("tprofile_v%d_pT_eventplane_fvtxsa_c%d_%d", harmonic, ic, energy));
    hvnpt_fvtxsa_east->SetName(Form("tprofile_v%d_pT_east_eventplane_fvtxsa_c%d_%d", harmonic, ic, energy));
    hvnpt_fvtxsa_west->SetName(Form("tprofile_v%d_pT_west_eventplane_fvtxsa_c%d_%d", harmonic, ic, energy));
    // hvneta_fvtxsa_both->SetName(Form("tprofile_v%d_eta_eventplane_fvtxsa_c%d_%d", harmonic, ic, energy));
    // hvneta_fvtxsa_east->SetName(Form("tprofile_v%d_eta_east_eventplane_fvtxsa_c%d_%d", harmonic, ic, energy));
    // hvneta_fvtxsa_west->SetName(Form("tprofile_v%d_eta_west_eventplane_fvtxsa_c%d_%d", harmonic, ic, energy));

    hvnpt_fvtxsb_both->SetName(Form("tprofile_v%d_pT_eventplane_fvtxsb_c%d_%d", harmonic, ic, energy));
    hvnpt_fvtxsb_east->SetName(Form("tprofile_v%d_pT_east_eventplane_fvtxsb_c%d_%d", harmonic, ic, energy));
    hvnpt_fvtxsb_west->SetName(Form("tprofile_v%d_pT_west_eventplane_fvtxsb_c%d_%d", harmonic, ic, energy));
    // hvneta_fvtxsb_both->SetName(Form("tprofile_v%d_eta_eventplane_fvtxsb_c%d_%d", harmonic, ic, energy));
    // hvneta_fvtxsb_east->SetName(Form("tprofile_v%d_eta_east_eventplane_fvtxsb_c%d_%d", harmonic, ic, energy));
    // hvneta_fvtxsb_west->SetName(Form("tprofile_v%d_eta_west_eventplane_fvtxsb_c%d_%d", harmonic, ic, energy));

    for (int il = 0; il < NFVTXLAY; il++)
    {
      hvnpt_fvtxsl_both[il]->SetName(Form("tprofile_v%d_pT_eventplane_fvtxsl%d_c%d_%d", harmonic, il, ic, energy));
      hvnpt_fvtxsl_east[il]->SetName(Form("tprofile_v%d_pT_east_eventplane_fvtxsl%d_c%d_%d", harmonic, il, ic, energy));
      hvnpt_fvtxsl_west[il]->SetName(Form("tprofile_v%d_pT_west_eventplane_fvtxsl%d_c%d_%d", harmonic, il, ic, energy));
      // hvneta_fvtxsl_both[il]->SetName(Form("tprofile_v%d_eta_eventplane_fvtxsl%d_c%d_%d", harmonic, il, ic, energy));
      // hvneta_fvtxsl_east[il]->SetName(Form("tprofile_v%d_eta_east_eventplane_fvtxsl%d_c%d_%d", harmonic, il, ic, energy));
      // hvneta_fvtxsl_west[il]->SetName(Form("tprofile_v%d_eta_west_eventplane_fvtxsl%d_c%d_%d", harmonic, il, ic, energy));
    } // il

    hvnpt_bbcs_east_zres->SetName(Form("tprofile_v%d_pT_east_eventplane_bbcs_zres_c%d_%d", harmonic, ic, energy));
    hvnpt_bbcs_west_zres->SetName(Form("tprofile_v%d_pT_west_eventplane_bbcs_zres_c%d_%d", harmonic, ic, energy));
    hvnpt_bbcs_both_zres->SetName(Form("tprofile_v%d_pT_eventplane_bbcs_zres_c%d_%d", harmonic, ic, energy));

    hvnpt_fvtxs_east_zres->SetName(Form("tprofile_v%d_pT_east_eventplane_fvtxs_zres_c%d_%d", harmonic, ic, energy));
    hvnpt_fvtxs_west_zres->SetName(Form("tprofile_v%d_pT_west_eventplane_fvtxs_zres_c%d_%d", harmonic, ic, energy));
    hvnpt_fvtxs_both_zres->SetName(Form("tprofile_v%d_pT_eventplane_fvtxs_zres_c%d_%d", harmonic, ic, energy));

    for (int il = 0; il < NFVTXLAY; il++)
    {
      hvnpt_fvtxsl_east_zres[il]->SetName(Form("tprofile_v%d_pT_east_eventplane_fvtxsl%d_zres_c%d_%d", harmonic, il, ic, energy));
      hvnpt_fvtxsl_west_zres[il]->SetName(Form("tprofile_v%d_pT_west_eventplane_fvtxsl%d_zres_c%d_%d", harmonic, il, ic, energy));
      hvnpt_fvtxsl_both_zres[il]->SetName(Form("tprofile_v%d_pT_eventplane_fvtxsl%d_zres_c%d_%d", harmonic, il, ic, energy));
    } // il


    for (int iz = 0; iz < NZPS; iz++)
    {
      hvnpt_bbcs_east_zvtx[iz]->SetName(Form("tprofile_v%d_pT_east_eventplane_bbcs_z%d_c%d_%d", harmonic, iz, ic, energy));
      hvnpt_bbcs_west_zvtx[iz]->SetName(Form("tprofile_v%d_pT_west_eventplane_bbcs_z%d_c%d_%d", harmonic, iz, ic, energy));
      hvnpt_bbcs_both_zvtx[iz]->SetName(Form("tprofile_v%d_pT_eventplane_bbcs_z%d_c%d_%d", harmonic, iz, ic, energy));

      hvnpt_fvtxs_east_zvtx[iz]->SetName(Form("tprofile_v%d_pT_east_eventplane_fvtxs_z%d_c%d_%d", harmonic, iz, ic, energy));
      hvnpt_fvtxs_west_zvtx[iz]->SetName(Form("tprofile_v%d_pT_west_eventplane_fvtxs_z%d_c%d_%d", harmonic, iz, ic, energy));
      hvnpt_fvtxs_both_zvtx[iz]->SetName(Form("tprofile_v%d_pT_eventplane_fvtxs_z%d_c%d_%d", harmonic, iz, ic, energy));

      for (int il = 0; il < NFVTXLAY; il++)
      {
        hvnpt_fvtxsl_east_zvtx[iz][il]->SetName(Form("tprofile_v%d_pT_east_eventplane_fvtxsl%d_z%d_c%d_%d", harmonic, il, iz, ic, energy));
        hvnpt_fvtxsl_west_zvtx[iz][il]->SetName(Form("tprofile_v%d_pT_west_eventplane_fvtxsl%d_z%d_c%d_%d", harmonic, il, iz, ic, energy));
        hvnpt_fvtxsl_both_zvtx[iz][il]->SetName(Form("tprofile_v%d_pT_eventplane_fvtxsl%d_z%d_c%d_%d", harmonic, il, iz, ic, energy));
      } // il

    } // iz

    gvn_bbcs_zdep->SetName(Form("tgraph_v%d_zvtx_eventplane_bbcs_c%d_%d", harmonic, ic, energy));
    gvn_bbcs_zdep_sys->SetName(Form("tgraph_v%d_zvtx_sys_eventplane_bbcs_c%d_%d", harmonic, ic, energy));

    gvn_fvtxs_zdep->SetName(Form("tgraph_v%d_zvtx_eventplane_fvtxs_c%d_%d", harmonic, ic, energy));
    gvn_fvtxs_zdep_sys->SetName(Form("tgraph_v%d_zvtx_sys_eventplane_fvtxs_c%d_%d", harmonic, ic, energy));

    hvnpt_bbcs_east->Write();
    hvnpt_bbcs_west->Write();
    hvnpt_bbcs_both->Write();
    hvneta_bbcs_east->Write();
    hvneta_bbcs_west->Write();
    hvneta_bbcs_both->Write();

    hvnpt_fvtxs_east->Write();
    hvnpt_fvtxs_west->Write();
    hvnpt_fvtxs_both->Write();
    hvneta_fvtxs_east->Write();
    hvneta_fvtxs_west->Write();
    hvneta_fvtxs_both->Write();

    hvnpt_fvtxsa_east->Write();
    hvnpt_fvtxsa_west->Write();
    hvnpt_fvtxsa_both->Write();
    // hvneta_fvtxsa_east->Write();
    // hvneta_fvtxsa_west->Write();
    // hvneta_fvtxsa_both->Write();

    hvnpt_fvtxsb_east->Write();
    hvnpt_fvtxsb_west->Write();
    hvnpt_fvtxsb_both->Write();
    // hvneta_fvtxsb_east->Write();
    // hvneta_fvtxsb_west->Write();
    // hvneta_fvtxsb_both->Write();

    hvnpt_bbcs_east_zres->Write();
    hvnpt_bbcs_west_zres->Write();
    hvnpt_bbcs_both_zres->Write();

    hvnpt_fvtxs_east_zres->Write();
    hvnpt_fvtxs_west_zres->Write();
    hvnpt_fvtxs_both_zres->Write();

    for (int il = 0; il < NFVTXLAY; il++)
    {
      hvnpt_fvtxsl_east[il]->Write();
      hvnpt_fvtxsl_west[il]->Write();
      hvnpt_fvtxsl_both[il]->Write();
      // hvneta_fvtxsl_east[il]->Write();
      // hvneta_fvtxsl_west[il]->Write();
      // hvneta_fvtxsl_both[il]->Write();

      hvnpt_fvtxsl_east_zres[il]->Write();
      hvnpt_fvtxsl_west_zres[il]->Write();
      hvnpt_fvtxsl_both_zres[il]->Write();
    }

    for (int iz = 0; iz < NZPS; iz++)
    {
      hvnpt_bbcs_east_zvtx[iz]->Write();
      hvnpt_bbcs_west_zvtx[iz]->Write();
      hvnpt_bbcs_both_zvtx[iz]->Write();

      hvnpt_fvtxs_east_zvtx[iz]->Write();
      hvnpt_fvtxs_west_zvtx[iz]->Write();
      hvnpt_fvtxs_both_zvtx[iz]->Write();

      for (int il = 0; il < NFVTXLAY; il++)
      {
        hvnpt_fvtxsl_east_zvtx[iz][il]->Write();
        hvnpt_fvtxsl_west_zvtx[iz][il]->Write();
        hvnpt_fvtxsl_both_zvtx[iz][il]->Write();
      } // il
    } // iz

    gvn_bbcs_zdep->Write();
    gvn_bbcs_zdep_sys->Write();

    gvn_fvtxs_zdep->Write();
    gvn_fvtxs_zdep_sys->Write();

  } // ic


  outfile->cd();
  th1d_vncent_fvtxs_lowpt->SetName(Form("th1d_v%d_cent_lowpt_eventplane_fvtxs_%d", harmonic, energy));
  th1d_vncent_fvtxs_highpt->SetName(Form("th1d_v%d_cent_highpt_eventplane_fvtxs_%d", harmonic, energy));

  gvn_bbcs_both_cent->SetName(Form("tgraph_v%d_cent_eventplane_bbcs_%d", harmonic, energy));
  gvn_fvtxs_both_cent->SetName(Form("tgraph_v%d_cent_eventplane_fvtxs_%d", harmonic, energy));
  gvn_fvtxsa_both_cent->SetName(Form("tgraph_v%d_cent_eventplane_fvtxsa_%d", harmonic, energy));
  gvn_fvtxsb_both_cent->SetName(Form("tgraph_v%d_cent_eventplane_fvtxsb_%d", harmonic, energy));
  for (int il = 0; il < NFVTXLAY; il++)
    gvn_fvtxsl_both_cent[il]->SetName(Form("tgraph_v%d_cent_eventplane_fvtxsl%d_%d", harmonic, il, energy));

  gvn_bbcs_east_cent->SetName(Form("tgraph_v%d_cent_east_eventplane_bbcs_%d", harmonic, energy));
  gvn_fvtxs_east_cent->SetName(Form("tgraph_v%d_cent_east_eventplane_fvtxs_%d", harmonic, energy));
  gvn_fvtxsa_east_cent->SetName(Form("tgraph_v%d_cent_east_eventplane_fvtxsa_%d", harmonic, energy));
  gvn_fvtxsb_east_cent->SetName(Form("tgraph_v%d_cent_east_eventplane_fvtxsb_%d", harmonic, energy));
  for (int il = 0; il < NFVTXLAY; il++)
    gvn_fvtxsl_east_cent[il]->SetName(Form("tgraph_v%d_cent_east_eventplane_fvtxsl%d_%d", harmonic, il, energy));

  gvn_bbcs_west_cent->SetName(Form("tgraph_v%d_cent_west_eventplane_bbcs_%d", harmonic, energy));
  gvn_fvtxs_west_cent->SetName(Form("tgraph_v%d_cent_west_eventplane_fvtxs_%d", harmonic, energy));
  gvn_fvtxsa_west_cent->SetName(Form("tgraph_v%d_cent_west_eventplane_fvtxsa_%d", harmonic, energy));
  gvn_fvtxsb_west_cent->SetName(Form("tgraph_v%d_cent_west_eventplane_fvtxsb_%d", harmonic, energy));
  for (int il = 0; il < NFVTXLAY; il++)
    gvn_fvtxsl_west_cent[il]->SetName(Form("tgraph_v%d_cent_west_eventplane_fvtxsl%d_%d", harmonic, il, energy));

  th1d_vncent_fvtxs_lowpt->Write();
  th1d_vncent_fvtxs_highpt->Write();

  gvn_bbcs_both_cent->Write();
  gvn_fvtxs_both_cent->Write();
  gvn_fvtxsa_both_cent->Write();
  gvn_fvtxsb_both_cent->Write();
  for (int il = 0; il < NFVTXLAY; il++)
    gvn_fvtxsl_both_cent[il]->Write();

  gvn_bbcs_east_cent->Write();
  gvn_fvtxs_east_cent->Write();
  gvn_fvtxsa_east_cent->Write();
  gvn_fvtxsb_east_cent->Write();
  for (int il = 0; il < NFVTXLAY; il++)
    gvn_fvtxsl_east_cent[il]->Write();

  gvn_bbcs_west_cent->Write();
  gvn_fvtxs_west_cent->Write();
  gvn_fvtxsa_west_cent->Write();
  gvn_fvtxsb_west_cent->Write();
  for (int il = 0; il < NFVTXLAY; il++)
    gvn_fvtxsl_west_cent[il]->Write();



}


void diagnostic(int energy, int harmonic)
{

  gStyle->SetOptTitle(1);

  TCanvas* c1 = new TCanvas("c1", "");

  TFile* file = TFile::Open(Form("input/combined_%d.root", energy));



  float MIN = -1;
  float MAX = 1;



  // ---------
  // --- FVTXS
  // ---------

  TProfile* hvn_fvtxs_B = (TProfile*)file->Get(Form("fvtxs_v%d_both_docalib", harmonic));
  TProfile* hvn_fvtxs_E = (TProfile*)file->Get(Form("fvtxs_v%d_east_docalib", harmonic));
  TProfile* hvn_fvtxs_W = (TProfile*)file->Get(Form("fvtxs_v%d_west_docalib", harmonic));

  hvn_fvtxs_B->SetLineColor(kBlack);
  hvn_fvtxs_E->SetLineColor(kRed);
  hvn_fvtxs_W->SetLineColor(kBlue);


  if ( energy == 200 ) {MIN = -0.01; MAX = 0.05;}
  if ( energy == 62 ) {MIN = -0.01; MAX = 0.05;}
  if ( energy == 39 ) {MIN = -0.03; MAX = 0.05;}
  if ( energy == 20 ) {MIN = -0.1; MAX = 0.2;}
  if ( harmonic == 3 ) {MIN = -0.003; MAX = 0.008;}
  TLine line(0, 0, 3, 0);
  line.SetLineStyle(2);
  line.SetLineWidth(2);
  TH2D* h2dummy = new TH2D("h2dummy", "", 1, 0.0, 3.0, 1, MIN, MAX);
  h2dummy->Draw();
  line.Draw();
  h2dummy->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h2dummy->GetYaxis()->SetTitle(Form("v_{%d} not corrected for EP resolution", harmonic));
  h2dummy->GetYaxis()->SetTitleOffset(1.3);

  hvn_fvtxs_B->Draw("same");
  hvn_fvtxs_E->Draw("same");
  hvn_fvtxs_W->Draw("same");

  TLegend* leg = new TLegend(0.18, 0.68, 0.38, 0.88);
  leg->SetHeader(Form("%d GeV, FVTXS", energy));
  leg->AddEntry(hvn_fvtxs_B, "both", "el");
  leg->AddEntry(hvn_fvtxs_E, "east", "el");
  leg->AddEntry(hvn_fvtxs_W, "west", "el");
  leg->SetTextSize(0.05);
  leg->Draw();

  c1->Print(Form("FigsHarmonicCoefficient/diagnostic_fvtxs_EBW_energy%d_harm%d.png", energy, harmonic));
  c1->Print(Form("FigsHarmonicCoefficient/diagnostic_fvtxs_EBW_energy%d_harm%d.pdf", energy, harmonic));


  //if ( harmonic == 2 )
  //{
  TProfile* hvn_fvtxs0_B = (TProfile*)file->Get(Form("fvtxs0_v%d_both_docalib", harmonic));
  TProfile* hvn_fvtxs1_B = (TProfile*)file->Get(Form("fvtxs1_v%d_both_docalib", harmonic));
  TProfile* hvn_fvtxs2_B = (TProfile*)file->Get(Form("fvtxs2_v%d_both_docalib", harmonic));
  TProfile* hvn_fvtxs3_B = (TProfile*)file->Get(Form("fvtxs3_v%d_both_docalib", harmonic));

  hvn_fvtxs0_B->SetLineColor(kBlue);
  hvn_fvtxs1_B->SetLineColor(kRed);
  hvn_fvtxs2_B->SetLineColor(kGreen + 2);
  hvn_fvtxs3_B->SetLineColor(kMagenta + 2);

  if ( energy == 200 ) {MIN = -0.01; MAX = 0.05;}
  if ( energy == 62 ) {MIN = -0.01; MAX = 0.05;}
  if ( energy == 39 ) {MIN = -0.03; MAX = 0.05;}
  if ( energy == 20 ) {MIN = -0.1; MAX = 0.2;}
  if ( harmonic == 3 ) {MIN = -0.003; MAX = 0.008;}
  delete h2dummy;
  h2dummy = new TH2D("h2dummy", "", 1, 0.0, 3.0, 1, MIN, MAX);
  h2dummy->Draw();
  line.Draw();
  h2dummy->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h2dummy->GetYaxis()->SetTitle(Form("v_{%d} not corrected for EP resolution", harmonic));
  h2dummy->GetYaxis()->SetTitleOffset(1.3);

  hvn_fvtxs_B->Draw("same");
  hvn_fvtxs0_B->Draw("same");
  hvn_fvtxs1_B->Draw("same");
  hvn_fvtxs2_B->Draw("same");
  hvn_fvtxs3_B->Draw("same");

  delete leg;
  leg = new TLegend(0.18, 0.68, 0.38, 0.88);
  leg->SetHeader(Form("%d GeV, FVTXS", energy));
  leg->AddEntry(hvn_fvtxs_B, "all layers", "el");
  leg->AddEntry(hvn_fvtxs0_B, "layer 0", "el");
  leg->AddEntry(hvn_fvtxs1_B, "layer 1", "el");
  leg->AddEntry(hvn_fvtxs2_B, "layer 2", "el");
  leg->AddEntry(hvn_fvtxs3_B, "layer 3", "el");
  leg->SetTextSize(0.05);
  leg->Draw();

  c1->Print(Form("FigsHarmonicCoefficient/diagnostic_fvtxsL_B_energy%d_harm%d.png", energy, harmonic));
  c1->Print(Form("FigsHarmonicCoefficient/diagnostic_fvtxsL_B_energy%d_harm%d.pdf", energy, harmonic));
  //}


  // ---------
  // --- BBCS
  // ---------

  TProfile* hvn_bbcs_B = (TProfile*)file->Get(Form("bbcs_v%d_both_docalib", harmonic));
  TProfile* hvn_bbcs_E = (TProfile*)file->Get(Form("bbcs_v%d_east_docalib", harmonic));
  TProfile* hvn_bbcs_W = (TProfile*)file->Get(Form("bbcs_v%d_west_docalib", harmonic));

  hvn_bbcs_B->SetLineColor(kBlack);
  hvn_bbcs_E->SetLineColor(kRed);
  hvn_bbcs_W->SetLineColor(kBlue);

  if ( energy == 200 ) {MIN = -0.005; MAX = 0.03;}
  if ( energy == 62 ) {MIN = -0.005; MAX = 0.02;}
  if ( energy == 39 ) {MIN = -0.02; MAX = 0.02;}
  if ( energy == 20 ) {MIN = -0.05; MAX = 0.05;}
  if ( harmonic == 3 ) {MIN = -0.002; MAX = 0.008;}
  delete h2dummy;
  h2dummy = new TH2D("h2dummy", "", 1, 0.0, 3.0, 1, MIN, MAX);
  h2dummy->Draw();
  line.Draw();
  h2dummy->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h2dummy->GetYaxis()->SetTitle(Form("v_{%d} not corrected for EP resolution", harmonic));
  h2dummy->GetYaxis()->SetTitleOffset(1.3);

  hvn_bbcs_B->Draw("same");
  hvn_bbcs_E->Draw("same");
  hvn_bbcs_W->Draw("same");

  delete leg;
  leg = new TLegend(0.18, 0.68, 0.38, 0.88);
  leg->SetHeader(Form("%d GeV, BBCS", energy));
  leg->AddEntry(hvn_fvtxs_B, "both", "el");
  leg->AddEntry(hvn_fvtxs_E, "east", "el");
  leg->AddEntry(hvn_fvtxs_W, "west", "el");
  leg->SetTextSize(0.05);
  leg->Draw();

  c1->Print(Form("FigsHarmonicCoefficient/diagnostic_bbcs_EBW_energy%d_harm%d.png", energy, harmonic));
  c1->Print(Form("FigsHarmonicCoefficient/diagnostic_bbcs_EBW_energy%d_harm%d.pdf", energy, harmonic));



  // ---------
  // --- FVTXN
  // ---------

  TProfile* hvn_fvtxn_B = (TProfile*)file->Get(Form("fvtxn_v%d_both_docalib", harmonic));
  TProfile* hvn_fvtxn_E = (TProfile*)file->Get(Form("fvtxn_v%d_east_docalib", harmonic));
  TProfile* hvn_fvtxn_W = (TProfile*)file->Get(Form("fvtxn_v%d_west_docalib", harmonic));

  hvn_fvtxn_B->SetLineColor(kBlack);
  hvn_fvtxn_E->SetLineColor(kRed);
  hvn_fvtxn_W->SetLineColor(kBlue);

  if ( energy == 200 ) {MIN = -0.01; MAX = 0.05;}
  if ( energy == 62 ) {MIN = -0.01; MAX = 0.05;}
  if ( energy == 39 ) {MIN = -0.03; MAX = 0.05;}
  if ( energy == 20 ) {MIN = -0.1; MAX = 0.2;}
  if ( harmonic == 3 ) {MIN = -0.003; MAX = 0.008;}
  delete h2dummy;
  h2dummy = new TH2D("h2dummy", "", 1, 0.0, 3.0, 1, MIN, MAX);
  h2dummy->Draw();
  line.Draw();
  h2dummy->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h2dummy->GetYaxis()->SetTitle(Form("v_{%d} not corrected for EP resolution", harmonic));
  h2dummy->GetYaxis()->SetTitleOffset(1.3);

  hvn_fvtxn_B->Draw("same");
  hvn_fvtxn_E->Draw("same");
  hvn_fvtxn_W->Draw("same");

  TLegend* leg2 = new TLegend(0.18, 0.68, 0.38, 0.88);
  leg2->SetHeader(Form("%d GeV, FVTXN", energy));
  leg2->AddEntry(hvn_fvtxn_B, "both", "el");
  leg2->AddEntry(hvn_fvtxn_E, "east", "el");
  leg2->AddEntry(hvn_fvtxn_W, "west", "el");
  leg2->SetTextSize(0.05);
  leg2->Draw();

  c1->Print(Form("FigsHarmonicCoefficient/diagnostic_fvtxn_EBW_energy%d_harm%d.png", energy, harmonic));
  c1->Print(Form("FigsHarmonicCoefficient/diagnostic_fvtxn_EBW_energy%d_harm%d.pdf", energy, harmonic));



  delete c1;

}



ValErr calc_epreso(ValErr AB, ValErr AC, ValErr BC)
{

  double reso = sqrt( (AB.first * AC.first) / BC.first );
  double ereso = reso / 2. * sqrt( pow(AB.second / AB.first, 2) +
                                   pow(AC.second / AC.first, 2) +
                                   pow(BC.second / BC.first, 2) );

  return make_pair(reso, ereso);
}

ValErr int_vn(TProfile *pvn, int bl, int bh)
{
  float vn_int = 0;
  float vn_e = 0;
  float w = 0;
  for (int ix = bl; ix <= bh; ix++)
  {
    float bc, be;
    bc = pvn->GetBinContent(ix);
    be = pvn->GetBinError(ix);

    vn_int += bc / (be * be);
    w += 1. / (be * be);
  }
  vn_int /= w;
  vn_e = sqrt(1. / w);

  return make_pair(vn_int, vn_e);
}

