#include "../RpPar.h"

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TStyle.h>
#include <TString.h>
#include <TProfile.h>

#include <iostream>

using namespace std;

void doenergy(int, int);

void calculate_epreso()
{

  doenergy(200, 2);
  doenergy(62, 2);
  doenergy(39, 2);
  doenergy(20, 2);

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
  // --- Now Centrality dependent
  // ---
  TProfile* tp1f_BBCS_FVTXN[NMUL];
  TProfile* tp1f_BBCS_FVTXS[NMUL];
  TProfile* tp1f_FVTXN_FVTXS[NMUL];
  TProfile* tp1f_BBCS_CNT[NMUL];
  TProfile* tp1f_CNT_FVTXS[NMUL];
  TProfile* tp1f_CNT_FVTXN[NMUL];

  TString dcor[NMUL];
  TString dcnt[NMUL];
  TString dfb[NMUL];
  TString davg[NMUL];
  TString drat[NMUL];

  for (int ic = 0; ic < NMUL; ic++)
  {

    // --- get histograms from file
    tp1f_BBCS_FVTXN[ic] = (TProfile*)file->Get(Form("tp1f_c%d_reso%d_BBC_FVTXN", ic, harmonic));
    if (!tp1f_BBCS_FVTXN[ic])
    {
      cout << "ERROR!! Unable to find " << Form("tp1f_c%d_reso%d_BBC_FVTXN", ic, harmonic)
           << " in file" << endl;
      return;
    }
    tp1f_BBCS_FVTXS[ic] = (TProfile*)file->Get(Form("tp1f_c%d_reso%d_BBC_FVTX", ic, harmonic));
    tp1f_FVTXN_FVTXS[ic] = (TProfile*)file->Get(Form("tp1f_c%d_reso%d_FVTXS_FVTXN", ic, harmonic));
    tp1f_BBCS_CNT[ic] = (TProfile*)file->Get(Form("tp1f_c%d_reso%d_BBC_CNT", ic, harmonic));
    tp1f_CNT_FVTXS[ic] = (TProfile*)file->Get(Form("tp1f_c%d_reso%d_CNT_FVTX", ic, harmonic));
    tp1f_CNT_FVTXN[ic] = (TProfile*)file->Get(Form("tp1f_c%d_reso%d_CNT_FVTXN", ic, harmonic));

    //-- get means & errors
    float float_BBCS_FVTXN = tp1f_BBCS_FVTXN[ic]->GetBinContent(1);
    float float_BBCS_FVTXS = tp1f_BBCS_FVTXS[ic]->GetBinContent(1);
    float float_FVTXN_FVTXS = tp1f_FVTXN_FVTXS[ic]->GetBinContent(1);
    float float_BBCS_CNT = tp1f_BBCS_CNT[ic]->GetBinContent(1);
    float float_CNT_FVTXS = tp1f_CNT_FVTXS[ic]->GetBinContent(1);
    float float_CNT_FVTXN = tp1f_CNT_FVTXN[ic]->GetBinContent(1);

    float efloat_BBCS_FVTXN = tp1f_BBCS_FVTXN[ic]->GetBinError(1);
    float efloat_BBCS_FVTXS = tp1f_BBCS_FVTXS[ic]->GetBinError(1);
    float efloat_FVTXN_FVTXS = tp1f_FVTXN_FVTXS[ic]->GetBinError(1);
    float efloat_BBCS_CNT = tp1f_BBCS_CNT[ic]->GetBinError(1);
    float efloat_CNT_FVTXS = tp1f_CNT_FVTXS[ic]->GetBinError(1);
    float efloat_CNT_FVTXN = tp1f_CNT_FVTXN[ic]->GetBinError(1);

    dcor[ic] = Form("%d GeV & %d & %.2e $\\pm$ %.2e & %.2e $\\pm$ %.2e & %.2e $\\pm$ %.2e & %.2e $\\pm$ %.2e & %.2e $\\pm$ %.2e & %.2e $\\pm$ %.2e \\\\",
                    energy, ic,
                    float_BBCS_FVTXS, efloat_BBCS_FVTXS,
                    float_BBCS_FVTXN, efloat_BBCS_FVTXN,
                    float_BBCS_CNT, efloat_BBCS_CNT,
                    float_CNT_FVTXS, efloat_CNT_FVTXS,
                    float_CNT_FVTXN, efloat_CNT_FVTXN,
                    float_FVTXN_FVTXS, efloat_FVTXN_FVTXS
                   );


    // --- event planes and correlations using CNT-BBCS-FVTXS(N)
    float reso_BBCS = sqrt((float_BBCS_CNT * float_BBCS_FVTXS) / float_CNT_FVTXS); // BCBS/CS
    float reso_FVTXS = sqrt((float_CNT_FVTXS * float_BBCS_FVTXS) / float_BBCS_CNT); // CSBS/BC
    float reso_FVTXN = sqrt((float_CNT_FVTXN * float_BBCS_FVTXN) / float_BBCS_CNT); // CNBN/BC

    float ereso_BBCS = sqrt( ( efloat_BBCS_CNT * efloat_BBCS_CNT / 4 * float_BBCS_CNT )
                             + ( efloat_BBCS_FVTXS * efloat_BBCS_FVTXS / 4 * float_BBCS_FVTXS )
                             + ( efloat_CNT_FVTXS * efloat_CNT_FVTXS / 4 * pow(float_CNT_FVTXS, 3) ) );
    float ereso_FVTXS = sqrt( ( efloat_CNT_FVTXS * efloat_CNT_FVTXS / 4 * float_CNT_FVTXS )
                              + ( efloat_BBCS_FVTXS * efloat_BBCS_FVTXS / 4 * float_BBCS_FVTXS )
                              + ( efloat_BBCS_CNT * efloat_BBCS_CNT / 4 * pow(float_BBCS_CNT, 3) ) );
    float ereso_FVTXN = sqrt( ( efloat_CNT_FVTXN * efloat_CNT_FVTXN / 4 * float_CNT_FVTXN )
                              + ( efloat_BBCS_FVTXN * efloat_BBCS_FVTXN / 4 * float_BBCS_FVTXN )
                              + ( efloat_BBCS_CNT * efloat_BBCS_CNT / 4 * pow(float_BBCS_CNT, 3) ) );

    dcnt[ic] = Form("%d GeV & %d & %.2e $\\pm$ %.2e & %.2e $\\pm$ %.2e & %.2e $\\pm$ %.2e \\\\",
                            energy, ic,
                            reso_BBCS, ereso_BBCS,
                            reso_FVTXS, ereso_FVTXS,
                            reso_FVTXN, ereso_FVTXN
                           );



    // --- event planes and correlations using FVTXN-FVTXS-BBCS
    float reso_BBCS_fn  = sqrt((float_BBCS_FVTXN * float_BBCS_FVTXS) / float_FVTXN_FVTXS); // BNBS/NS
    float reso_FVTXS_fn = sqrt((float_FVTXN_FVTXS * float_BBCS_FVTXS) / float_BBCS_FVTXN); // NSBS/BN
    float reso_FVTXN_fn = sqrt((float_FVTXN_FVTXS * float_BBCS_FVTXN) / float_BBCS_FVTXS); // NSBN/BS

    float ereso_BBCS_fn  = sqrt( ( efloat_BBCS_FVTXN * efloat_BBCS_FVTXN / 4 * float_BBCS_FVTXN )
                                 + ( efloat_BBCS_FVTXS * efloat_BBCS_FVTXS / 4 * float_BBCS_FVTXS )
                                 + ( efloat_FVTXN_FVTXS * efloat_FVTXN_FVTXS / 4 * pow(float_FVTXN_FVTXS, 3) ) );
    float ereso_FVTXS_fn = sqrt( ( efloat_FVTXN_FVTXS * efloat_FVTXN_FVTXS / 4 * float_FVTXN_FVTXS )
                                 + ( efloat_BBCS_FVTXS * efloat_BBCS_FVTXS / 4 * float_BBCS_FVTXS )
                                 + ( efloat_BBCS_FVTXN * efloat_BBCS_FVTXN / 4 * pow(float_BBCS_FVTXN, 3) ) );
    float ereso_FVTXN_fn = sqrt( ( efloat_FVTXN_FVTXS * efloat_FVTXN_FVTXS / 4 * float_FVTXN_FVTXS )
                                 + ( efloat_BBCS_FVTXN * efloat_BBCS_FVTXN / 4 * float_BBCS_FVTXN )
                                 + ( efloat_BBCS_FVTXS * efloat_BBCS_FVTXS / 4 * pow(float_BBCS_FVTXS, 3) ) );


    dfb[ic] = Form("%d GeV & %d & %.2e $\\pm$ %.2e & %.2e $\\pm$ %.2e & %.2e $\\pm$ %.2e \\\\",
                           energy, ic,
                           reso_BBCS_fn, ereso_BBCS_fn,
                           reso_FVTXS_fn, ereso_FVTXS_fn,
                           reso_FVTXN_fn, ereso_FVTXN_fn
                          );


    // --- resolution comparisons
    float rat_BBCS  = reso_BBCS / reso_BBCS_fn;
    float rat_FVTXS = reso_FVTXS / reso_FVTXS_fn;
    float rat_FVTXN = reso_FVTXN / reso_FVTXN_fn;

    drat[ic] = Form("%d GeV & %d & %.2e & %.2e & %.2e \\\\",
                    energy, ic,
                    rat_BBCS, rat_FVTXS, rat_FVTXN);

    float avg_BBCS  = 0.5 * (reso_BBCS + reso_BBCS_fn);
    float avg_FVTXS = 0.5 * (reso_FVTXS + reso_FVTXS_fn);
    float avg_FVTXN = 0.5 * (reso_FVTXN + reso_FVTXN_fn);

    davg[ic] = Form("%d GeV & %d & %.2e & %.2e & %.2e \\\\",
                    energy, ic,
                    avg_BBCS, avg_FVTXS, avg_FVTXN);

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
  cout << " means:" << endl;
  cout << dm << endl;

  cout << endl;
  cout << " ep correlations:" << endl;
  for (int ic = 0; ic < NMUL; ic++)
    cout << dcor[ic] << endl;

  cout << endl;
  cout << " resos cnt:" << endl;
  for (int ic = 0; ic < NMUL; ic++)
    cout << dcnt[ic] << endl;

  cout << endl;
  cout << " resos fvtx:" << endl;
  for (int ic = 0; ic < NMUL; ic++)
    cout << dfb[ic] << endl;

  cout << endl;
  cout << " resos ratio:" << endl;
  for (int ic = 0; ic < NMUL; ic++)
    cout << drat[ic] << endl;

  cout << endl;
  cout << " resos avg:" << endl;
  for (int ic = 0; ic < NMUL; ic++)
    cout << davg[ic] << endl;


}

