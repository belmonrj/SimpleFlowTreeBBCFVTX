#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>

#ifndef RPPAR_H_
#include "RpPar.h"
#endif

//#ifndef __CINT__
#include "TROOT.h"
#include "TFile.h"
#include "TChain.h"
#include "TMath.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TVector2.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TRandom.h"
#include "TLorentzVector.h"
#include "TSystem.h"
#include "TComplex.h"
#include "TGraph.h"
//#endif

// global pmt variables
float d_pmt_x[64];
float d_pmt_y[64];
float d_pmt_z = -1443.5; // same for all tubes

//tree invariables
//static const int max_nh = 500; // see from ana taxi code
static const int max_nh = 50; // see from ana taxi code
//static const int max_nh = 10; // see from ana taxi code
static const int max_nf = 4000; // see from ana taxi code
//static const int max_nf = 3000; // see from ana taxi code

using namespace std;

//==================================================================//
//                                                                  //
//       The purpose of this compiled root macro is to              //
//       read in event wise ttrees, read in BBC and FVTX            //
//       event planes, perform flattening and recentering           //
//       and make v2 w/r/t vtx tracks.                              //
//                                                                  //
//       You need to run 3 iterations (rp_recal_pass = 1, 2, 3).    //
//       The first two passes are just for calibration, the third   //
//       pass is to actually create the v2. Calib text files are    //
//       created                                                    //
//                                                                  //
//       This module also has the capability to calculate the FVTX  //
//       BBC event planes from clusters and tubes respectively (and //
//       apply corrections).                                        //
//                                                                  //
//       By Theo Koblesky, May 11, 2016                             //
//       theodore.koblesky@colorado.edu                             //
//                                                                  //
//==================================================================//
//                                                                  //
//       Modify to calculate v2/v3 w/r/t CNT tracks using the EP    //
//       framework. Also calculate v2/v3 vs eta from CNT & FVTX     //
//       tracks. Remove all non-essential pieces to reduce runtime  //
//       and improve readability.                                   //
//                                                                  //
//       Darren McGlinchey                                          //
//       28 Feb 2017                                                //
//       darren.mcglinchey@colorado.edu                             //
//                                                                  //
//==================================================================//





bool DIAG = false;


int get_fvtx_layer(float);
void initialize_pmt_position();
int get_pmt_layer(int);
void flatten(int, int);
float calc2_event(float, float, float);
float calc4_event(float, float, float, float, float);
float calc4_track_flag(float, float, float, float, float, float, float, float, float, bool);


using namespace std;

int main(int argc, char *argv[])
{

  if ( argc != 2 )
  {
    cerr << "FATAL ERROR: this program takes the run number as an argument" << endl;
    return -1;
  }

  int run = atoi(argv[1]);

  cout << "Now processing with run number " << run << endl;

  // flatten(run, 1);
  // flatten(run, 2);
  flatten(run, 3);

  return 0;

}

// -----------------------------------------------------------------
void flatten(int runNumber, int rp_recal_pass)
{

  //------------------------------------------------------------//
  //                                                            //
  //       Set running conditions for this Run Number           //
  //                                                            //
  //------------------------------------------------------------//

  cout << "runNumber = " << runNumber << " rp_recal_pass = " << rp_recal_pass << endl;

  int energyflag = 0;
  // --- Run16dAu200
  if ( runNumber >= 454774 && runNumber <= 455639 ) energyflag = 200;
  // --- Run16dAu62
  if ( runNumber >= 455792 && runNumber <= 456283 ) energyflag = 62;
  // --- Run16dAu20
  if ( runNumber >= 456652 && runNumber <= 457298 ) energyflag = 20;
  // --- Run16dAu39
  if ( runNumber >= 457634 && runNumber <= 458167 ) energyflag = 39;

  int eidx = -1;
  if (energyflag == 200)
    eidx = 0;
  else if (energyflag == 62)
    eidx = 1;
  else if (energyflag == 39)
    eidx = 2;
  else if (energyflag == 20)
    eidx = 3;
  else
  {
    cout << "couldn't find valid energy index"
         << " runNumber:" << runNumber << " energyflag:" << energyflag << " eidx:" << eidx
         << endl;
    return;
  }

  int verbosity = 2;

  char calibfile[500];
  sprintf(calibfile, "output/flattening_data/flattening_%d_%d.dat", runNumber, rp_recal_pass - 1);

  cout << "calib text output file: " << calibfile << endl;

  char filename[500];

  // float fracCut = 0.95; // pile up rejection < fracCut (better place??)
  float fracCut = 0.98; // pile up rejection < fracCut (better place??)
  // float fracCut = 0.0; // pile up rejection < fracCut (better place??)

  float qxOffset = 0.0; // offset to Qx values
  // float qyOffset[] = { 0, 0, 0, 0, 0, 0}; // offset Qy values
  // centrality dependent Qy offsets
  // float qyOffset[] = { -0.001, -0.005, -0.008, -0.016, -0.029, -0.069}; // 1st iteration
  // float qyOffset[] = { -0.001, -0.007, -0.012, -0.024, -0.044, -0.104}; // 2nd iteration
  // float qyOffset[] = { -0.001, -0.008, -0.014, -0.028, -0.051, -0.121}; // 3rd iteration
  // float qyOffset[4][6] = {{0}, {0}, {0}, {0}};
  //---
  // full values (4th iteration)
  float qyOffset[4][6] =
  {
    { -0.000, -0.000, -0.000, -0.016, -0.044, -0.121}, //200 GeV

    { -0.001, -0.008, -0.014, -0.024, -0.055, -0.129}, //62 GeV

    { -0.000, -0.005, -0.014, -0.028, -0.051, -0.104}, //39 GeV

    { -0.001, -0.008, -0.014, -0.028, -0.051, -0.121}, //20 GeV
  };
  //---
  //---
  // Systematic check (x1.20)
  // float qyOffset[4][6] =
  // {
  //   { -0.000, -0.000, -0.000, -0.019, -0.053, -0.145}, //200 GeV

  //   { -0.001, -0.010, -0.017, -0.029, -0.066, -0.155}, //62 GeV

  //   { -0.000, -0.006, -0.017, -0.034, -0.061, -0.125}, //39 GeV

  //   { -0.001, -0.010, -0.017, -0.034, -0.061, -0.145}, //20 GeV
  // };
  //---
  //---
  // Systematic check (x0.80)
  // float qyOffset[4][6] =
  // {
  //   { -0.000, -0.000, -0.000, -0.013, -0.035, -0.097}, //200 GeV

  //   { -0.001, -0.006, -0.011, -0.019, -0.044, -0.103}, //62 GeV

  //   { -0.000, -0.004, -0.011, -0.022, -0.041, -0.083}, //39 GeV

  //   { -0.001, -0.006, -0.011, -0.022, -0.041, -0.097}, //20 GeV
  // };
  //---

  cout << " frac cut: " << fracCut << endl;
  cout << " Qx offset: " << qxOffset << endl;
  cout << " Qy offset: " << endl;
  for (int i = 0; i < NMUL; i++)
    cout << "    " << i << " " << qyOffset[eidx][i] << endl;

  //------------------------------------------------------------//
  //                                                            //
  //       Chain files for this Run Number                      //
  //                                                            //
  //------------------------------------------------------------//


  // --- get the number of files for this run number
  string pipe_out = (string) gSystem->GetFromPipe(Form("ls input/%d_*.root | grep -c r", runNumber));
  int nfiles = 0;
  nfiles = atoi(pipe_out.c_str());
  cout << "nfiles: " << nfiles << endl;
  if (nfiles == 0) return;

  // --- make a new TChain for the tree
  TChain *ntp_event_chain = new TChain("ntp_event");
  for ( int ifile = 0; ifile < nfiles; ++ifile )
  {
    sprintf(filename, "input/%d_%d.root", runNumber, ifile);
    cout << "adding to tchain: " << filename << endl;
    ntp_event_chain->Add(filename);
  }




  //------------------------------------------------------------//
  //                                                            //
  //       Initializing Calibration Arrays & Histograms         //
  //                                                            //
  //------------------------------------------------------------//

  // ---

  TString tubegaincorrectionfile_name = Form("SpecialProjects/WeightFiles/bbctube_run%d.root", runNumber);
  TFile* tube_gaincorrection_file = TFile::Open(tubegaincorrectionfile_name);
  float tube_gaincorrection[128];
  for ( int i = 0; i < 128; ++i ) tube_gaincorrection[i] = 1.0;
  if ( tube_gaincorrection_file )
  {
    TH1D* th1d_tube_gaincorrection = (TH1D*)tube_gaincorrection_file->Get("th1d_tubegaincorrection");
    if ( th1d_tube_gaincorrection )
    {
      for ( int i = 0; i < 64; ++i )
      {
        tube_gaincorrection[i] = th1d_tube_gaincorrection->GetBinContent(i + 1);
      }
      cout << "All tube gain gain values collected and ready for correction" << endl;
    }
    else cout << "WARNING could not get BBC tube gaincorrection histogram" << endl;
  }
  else cout << "WARNING could not get BBC tube gaincorrection file " << tubegaincorrectionfile_name.Data() << endl;
  // ---

  cout << "Initalizing PMT positions for the BBC" << endl;

  initialize_pmt_position();

  cout << "Setting flags and indecies" << endl;

  bool fvtx_clusters = true;
  bool bbc_pmts      = true;
  bool cnt_tracks    = true;
  bool fvtx_tracks   = true;


  int bbcs_index      =  0;
  int fvtxs_index     =  1;
  int fvtxn_index     =  2;


  //------------------------------------------------------------//
  //                                                            //
  //       Initializing Calibration Arrays & Histograms         //
  //                                                            //
  //------------------------------------------------------------//

  // --- problems with these array dimensions??  lots of compile errors

  cout << "Lots of arrays and stuff" << endl;



  char outFile1[300];
  sprintf(outFile1, "%s%d%s%d%s", "output/files_", energyflag, "/hist_", runNumber, ".root");
  cout << "histogram output file: " << outFile1 << endl;

  TFile *mData1 = TFile::Open(outFile1, "recreate");
  mData1->cd();


  // --- lots of comments needed here
  // --- flattening parameters output to file
  TProfile *ave[NMUL][NZPS][NHAR][NDETSHORT]; // average Psi
  TProfile *flt[NMUL][NZPS][NHAR][NDETSHORT]; // flattening parameters

  TH2D     *psi_bf[NMUL][NHAR][NDETSHORT];
  TH2D     *psi_mf[NMUL][NHAR][NDETSHORT];
  TH2D     *psi_af[NMUL][NHAR][NDETSHORT];

  // flattening parameters read in from file
  float    mean[NMUL][NZPS][NHAR][NDETSHORT][2]; // mean of qx, qy (? double check)
  float    widt[NMUL][NZPS][NHAR][NDETSHORT][2]; // width of qx, qy (?? double check)
  float    four[NMUL][NZPS][NHAR][NDETSHORT][2][NORD]; // fourier components for flattening



  // --- diagnostic histograms
  TH1D* th1d_BBC_charge = new TH1D("th1d_BBC_charge", "", 200, -0.5, 199.5);
  TH1D* th1d_FVTX_nclus = new TH1D("th1d_FVTX_nclus", "", 200, -0.5, 1999.5);
  TH1D* th1d_FVTX_ntrk = new TH1D("th1d_FVTX_ntrk", "", 100, -0.5, 99.5);
  TH2D* th2d_qBBC_nFVTX = new TH2D("th2d_qBBC_nFVTX", "", 200, -0.5, 199.5, 200, -0.5, 1999.5);
  TH1D* th1d_FVTXS_nclus = new TH1D("th1d_FVTXS_nclus", "", 200, -0.5, 1999.5);
  TH1D* th1d_FVTXN_nclus = new TH1D("th1d_FVTXN_nclus", "", 200, -0.5, 1999.5);

  TProfile* tp1f_bbc_fcharge_ring = new TProfile("tp1f_bbc_fcharge_ring", ";BBC ring; charge/total", 5, -0.5, 4.5);


  // --- event plane resolution, need to be centrality dependent
  TProfile* tp1f_reso2_BBC_CNT[NMUL];
  TProfile* tp1f_reso2_BBC_FVTX[NMUL];
  TProfile* tp1f_reso2_CNT_FVTX[NMUL];
  TProfile* tp1f_reso2_BBC_FVTXN[NMUL];
  TProfile* tp1f_reso2_CNT_FVTXN[NMUL];
  TProfile* tp1f_reso2_FVTXS_FVTXN[NMUL];

  TProfile* tp1f_reso3_BBC_CNT[NMUL];
  TProfile* tp1f_reso3_BBC_FVTX[NMUL];
  TProfile* tp1f_reso3_CNT_FVTX[NMUL];
  TProfile* tp1f_reso3_BBC_FVTXN[NMUL];
  TProfile* tp1f_reso3_CNT_FVTXN[NMUL];
  TProfile* tp1f_reso3_FVTXS_FVTXN[NMUL];

  for (int ic = 0; ic < NMUL; ic++)
  {
    tp1f_reso2_BBC_CNT[ic]  = new TProfile(Form("tp1f_c%i_reso2_BBC_CNT", ic), "", 1, -0.5, 0.5, -1e6, 1e6, "");
    tp1f_reso2_BBC_FVTX[ic] = new TProfile(Form("tp1f_c%i_reso2_BBC_FVTX", ic), "", 1, -0.5, 0.5, -1e6, 1e6, "");
    tp1f_reso2_CNT_FVTX[ic] = new TProfile(Form("tp1f_c%i_reso2_CNT_FVTX", ic), "", 1, -0.5, 0.5, -1e6, 1e6, "");
    tp1f_reso2_BBC_FVTXN[ic] = new TProfile(Form("tp1f_c%i_reso2_BBC_FVTXN", ic), "", 1, -0.5, 0.5, -1e6, 1e6, "");
    tp1f_reso2_CNT_FVTXN[ic] = new TProfile(Form("tp1f_c%i_reso2_CNT_FVTXN", ic), "", 1, -0.5, 0.5, -1e6, 1e6, "");
    tp1f_reso2_FVTXS_FVTXN[ic] = new TProfile(Form("tp1f_c%i_reso2_FVTXS_FVTXN", ic), "", 1, -0.5, 0.5, -1e6, 1e6, "");

    tp1f_reso3_BBC_CNT[ic]  = new TProfile(Form("tp1f_c%i_reso3_BBC_CNT", ic), "", 1, -0.5, 0.5, -1e6, 1e6, "");
    tp1f_reso3_BBC_FVTX[ic] = new TProfile(Form("tp1f_c%i_reso3_BBC_FVTX", ic), "", 1, -0.5, 0.5, -1e6, 1e6, "");
    tp1f_reso3_CNT_FVTX[ic] = new TProfile(Form("tp1f_c%i_reso3_CNT_FVTX", ic), "", 1, -0.5, 0.5, -1e6, 1e6, "");
    tp1f_reso3_BBC_FVTXN[ic] = new TProfile(Form("tp1f_c%i_reso3_BBC_FVTXN", ic), "", 1, -0.5, 0.5, -1e6, 1e6, "");
    tp1f_reso3_CNT_FVTXN[ic] = new TProfile(Form("tp1f_c%i_reso3_CNT_FVTXN", ic), "", 1, -0.5, 0.5, -1e6, 1e6, "");
    tp1f_reso3_FVTXS_FVTXN[ic] = new TProfile(Form("tp1f_c%i_reso3_FVTXS_FVTXN", ic), "", 1, -0.5, 0.5, -1e6, 1e6, "");
  } // ic



  cout << "Making TProfile histograms" << endl;

  // --- profile histograms for average of Psi and flattening parameters
  char name[200];
  for ( int ic = 0; ic < NMUL; ic ++ )
  {
    for ( int iz = 0; iz < NZPS; iz++ )
    {
      for ( int ih = 1; ih < NHAR; ih++ )
      {
        for ( int id = 0; id < NDETSHORT; id++ )
        {
          // --- average (of?)
          sprintf(name, "ave_%d_%d_%d_%d", ic, iz, ih, id);
          ave[ic][iz][ih][id] = new TProfile(name, name, 4, -0.5, 3.5, -10.1, 10.1, "S"); //for SMD -1.1,1.1
          // --- flattening parameter (?)
          sprintf(name, "flt_%d_%d_%d_%d", ic, iz, ih, id);
          flt[ic][iz][ih][id] = new TProfile(name, name, 4 * NORD, -0.5, NORD * 4.0 - 0.5, -1.1, 1.1);
        } // loop over NDETSHORTectors
      } // loop over harmonics
    } // loop over z_vertex bins
  } // loop over centrality bins

  // --- TH2D histograms for Q vector components
  for ( int ic = 0; ic < NMUL; ic++ )
  {
    for ( int ih = 1; ih < NHAR; ih++)
    {
      for ( int id = 0; id < NDETSHORT; id++)
      {
        // --- psi_bf (event plane before recentering and flattening)
        sprintf(name, "psi_bf_%d_%d_%d", ic, ih, id);
        psi_bf[ic][ih][id] = new TH2D(name, name, NZPS * 3, -0.5, NZPS * 3.0 - 0.5, 220, -4.1, 4.1);

        // --- psi_mf (event plane after recentering but before flattening)
        sprintf(name, "psi_mf_%d_%d_%d", ic, ih, id);
        psi_mf[ic][ih][id] = new TH2D(name, name, NZPS * 3, -0.5, NZPS * 3.0 - 0.5, 220, -4.1, 4.1);

        // --- psi_af (event plane after recentering and flattening)
        sprintf(name, "psi_af_%d_%d_%d", ic, ih, id);
        psi_af[ic][ih][id] = new TH2D(name, name, NZPS * 3, -0.5, NZPS * 3.0 - 0.5, 220, -4.1, 4.1);
      } // loop over detectors
    } // loop over harmonics
  } // loop over centrality bins

  cout << "Initalizing calibration parameters to zero" << endl;

  //Initializing the calibration parameters to be read in from the file
  for ( int ic = 0; ic < NMUL; ic ++ )
  {
    for ( int iz = 0; iz < NZPS; iz++ )
    {
      for ( int ih = 1; ih < NHAR; ih++ )
      {
        for ( int id = 0; id < NDETSHORT; id++ )
        {
          for ( int ib = 0; ib < 2; ib++ )
          {
            // --- mean (of q-vectors)
            mean[ic][iz][ih][id][ib] = 0.0;

            // --- width (of q-vectors)
            widt[ic][iz][ih][id][ib] = 1.0;

            // --- fourier components for flattening
            for ( int io = 0; io < NORD; io++ )
            {
              four[ic][iz][ih][id][ib][io] = 0.0;
            } // orders
          } // x and y
        } // detector
      } // harmonics
    } // z_vertex bins
  } // centrality bins

  //------------------------------------------------------------//
  //   Finished Initializing Calibration Arrays & Histograms    //
  //------------------------------------------------------------//

  //------------------------------------------------------------//
  //                                                            //
  //      Reading in flattening calibration parameters          //
  //                                                            //
  //------------------------------------------------------------//

  cout << "Checking for calibration pass" << endl;

  // --- dont read in for first pass
  if ( rp_recal_pass >= 2 )
  {
    cout << "reading calibration file : " << calibfile << endl;
    float f0, f1, f2, f3; //f4,f5,f6,f7;
    ifstream ifs;
    ifs.open(calibfile);
    for ( int ic = 0; ic < NMUL; ic++ )
    {
      for ( int iz = 0; iz < NZPS; iz++ )
      {
        for ( int ih = 1; ih < NHAR; ih++ )
        {
          for ( int id = 0; id < NDETSHORT; id++ )
          {
            ifs >> f0 >> f1 >> f2 >> f3;
            if ( f1 <= 0.0 ) f1 = 1.0;
            if ( f3 <= 0.0 ) f3 = 1.0;
            mean[ic][iz][ih][id][0] = f0;
            widt[ic][iz][ih][id][0] = f1;
            mean[ic][iz][ih][id][1] = f2;
            widt[ic][iz][ih][id][1] = f3;
            if ( id == 2 && ih == 1 && DIAG ) cout << f0 << " " << f1 << " " << f2 << " " << f3 << endl; //bbc psi 2 parameters
            if ( id == 2 && ih == 1 && DIAG ) cout << "---" << endl;
            if ( id == 3 && ih == 1 && DIAG ) cout << f0 << " " << f1 << " " << f2 << " " << f3 << endl; //bbc psi 2 parameters
            for ( int ib = 0; ib < 2; ib++ )
            {
              for ( int io = 0; io < NORD; io++ )
              {
                ifs >> four[ic][iz][ih][id][ib][io];
              } // orders
            } // x and y
          } // detectors
        } // harmonics
      } // z_vertex bins
    } // centrality bins
    ifs.close();
  } // check on second or third pass


  //------------------------------------------------------------//
  //  Finished Reading in flattening calibration parameters     //
  //------------------------------------------------------------//



  //------------------------------------------------------------//
  //                  Initializing histograms                   //
  //------------------------------------------------------------//

  cout << "Initializing more histograms" << endl;

  // --- docalib

  // vs pT
  const int NPTBINS = 16;
  double ptlim[NPTBINS + 1] = {
    0.0, 0.2, 0.4, 0.6, 0.8,
    1.0, 1.2, 1.4, 1.6, 1.8,
    2.0, 2.5, 3.0, 3.5, 4.0,
    4.5, 5.0
  };
  TProfile* fvtxs_v2_west_cosphi[NMUL];
  TProfile* fvtxs_v2_east_cosphi[NMUL];
  TProfile* fvtxs_v2_both_cosphi[NMUL];


  TProfile* fvtxs_v2_west_sinphi[NMUL];
  TProfile* fvtxs_v2_east_sinphi[NMUL];
  TProfile* fvtxs_v2_both_sinphi[NMUL];


  TProfile* bbcs_v2_east_docalib[NMUL];
  TProfile* bbcs_v2_west_docalib[NMUL];
  TProfile* bbcs_v2_both_docalib[NMUL];

  TProfile* fvtxs_v2_east_docalib[NMUL];
  TProfile* fvtxs_v2_west_docalib[NMUL];
  TProfile* fvtxs_v2_both_docalib[NMUL];

  TProfile* fvtxn_v2_east_docalib[NMUL];
  TProfile* fvtxn_v2_west_docalib[NMUL];
  TProfile* fvtxn_v2_both_docalib[NMUL];

  TProfile* bbcs_v3_east_docalib[NMUL];
  TProfile* bbcs_v3_west_docalib[NMUL];
  TProfile* bbcs_v3_both_docalib[NMUL];

  TProfile* fvtxs_v3_east_docalib[NMUL];
  TProfile* fvtxs_v3_west_docalib[NMUL];
  TProfile* fvtxs_v3_both_docalib[NMUL];

  TProfile* fvtxn_v3_east_docalib[NMUL];
  TProfile* fvtxn_v3_west_docalib[NMUL];
  TProfile* fvtxn_v3_both_docalib[NMUL];

  for ( int ic = 0; ic < NMUL; ++ic )
  {
    //-- cos(phi)
    fvtxs_v2_west_cosphi[ic] = new TProfile(Form("fvtxs_v2_west_cosphi_cent%d", ic), "", NPTBINS, ptlim, -1.1, 1.1);
    fvtxs_v2_east_cosphi[ic] = new TProfile(Form("fvtxs_v2_east_cosphi_cent%d", ic), "", NPTBINS, ptlim, -1.1, 1.1);
    fvtxs_v2_both_cosphi[ic] = new TProfile(Form("fvtxs_v2_both_cosphi_cent%d", ic), "", NPTBINS, ptlim, -1.1, 1.1);

    //-- sin(phi)
    fvtxs_v2_west_sinphi[ic] = new TProfile(Form("fvtxs_v2_west_sinphi_cent%d", ic), "", NPTBINS, ptlim, -1.1, 1.1);
    fvtxs_v2_east_sinphi[ic] = new TProfile(Form("fvtxs_v2_east_sinphi_cent%d", ic), "", NPTBINS, ptlim, -1.1, 1.1);
    fvtxs_v2_both_sinphi[ic] = new TProfile(Form("fvtxs_v2_both_sinphi_cent%d", ic), "", NPTBINS, ptlim, -1.1, 1.1);

    //-- v2
    bbcs_v2_both_docalib[ic] = new TProfile(Form("bbcs_v2_both_docalib_cent%d", ic), "", NPTBINS, ptlim, -1.1, 1.1);
    bbcs_v2_east_docalib[ic] = new TProfile(Form("bbcs_v2_east_docalib_cent%d", ic), "", NPTBINS, ptlim, -1.1, 1.1);
    bbcs_v2_west_docalib[ic] = new TProfile(Form("bbcs_v2_west_docalib_cent%d", ic), "", NPTBINS, ptlim, -1.1, 1.1);

    fvtxs_v2_both_docalib[ic] = new TProfile(Form("fvtxs_v2_both_docalib_cent%d", ic), "", NPTBINS, ptlim, -1.1, 1.1);
    fvtxs_v2_east_docalib[ic] = new TProfile(Form("fvtxs_v2_east_docalib_cent%d", ic), "", NPTBINS, ptlim, -1.1, 1.1);
    fvtxs_v2_west_docalib[ic] = new TProfile(Form("fvtxs_v2_west_docalib_cent%d", ic), "", NPTBINS, ptlim, -1.1, 1.1);

    fvtxn_v2_both_docalib[ic] = new TProfile(Form("fvtxn_v2_both_docalib_cent%d", ic), "", NPTBINS, ptlim, -1.1, 1.1);
    fvtxn_v2_east_docalib[ic] = new TProfile(Form("fvtxn_v2_east_docalib_cent%d", ic), "", NPTBINS, ptlim, -1.1, 1.1);
    fvtxn_v2_west_docalib[ic] = new TProfile(Form("fvtxn_v2_west_docalib_cent%d", ic), "", NPTBINS, ptlim, -1.1, 1.1);

    bbcs_v3_both_docalib[ic] = new TProfile(Form("bbcs_v3_both_docalib_cent%d", ic), "", NPTBINS, ptlim, -1.1, 1.1);
    bbcs_v3_east_docalib[ic] = new TProfile(Form("bbcs_v3_east_docalib_cent%d", ic), "", NPTBINS, ptlim, -1.1, 1.1);
    bbcs_v3_west_docalib[ic] = new TProfile(Form("bbcs_v3_west_docalib_cent%d", ic), "", NPTBINS, ptlim, -1.1, 1.1);

    fvtxs_v3_both_docalib[ic] = new TProfile(Form("fvtxs_v3_both_docalib_cent%d", ic), "", NPTBINS, ptlim, -1.1, 1.1);
    fvtxs_v3_east_docalib[ic] = new TProfile(Form("fvtxs_v3_east_docalib_cent%d", ic), "", NPTBINS, ptlim, -1.1, 1.1);
    fvtxs_v3_west_docalib[ic] = new TProfile(Form("fvtxs_v3_west_docalib_cent%d", ic), "", NPTBINS, ptlim, -1.1, 1.1);

    fvtxn_v3_both_docalib[ic] = new TProfile(Form("fvtxn_v3_both_docalib_cent%d", ic), "", NPTBINS, ptlim, -1.1, 1.1);
    fvtxn_v3_east_docalib[ic] = new TProfile(Form("fvtxn_v3_east_docalib_cent%d", ic), "", NPTBINS, ptlim, -1.1, 1.1);
    fvtxn_v3_west_docalib[ic] = new TProfile(Form("fvtxn_v3_west_docalib_cent%d", ic), "", NPTBINS, ptlim, -1.1, 1.1);
  } // ic


  // vs eta
  TProfile* bbcs_v2eta_east_docalib[NMUL];
  TProfile* bbcs_v2eta_west_docalib[NMUL];
  TProfile* bbcs_v2eta_both_docalib[NMUL];

  TProfile* fvtxs_v2eta_east_docalib[NMUL];
  TProfile* fvtxs_v2eta_west_docalib[NMUL];
  TProfile* fvtxs_v2eta_both_docalib[NMUL];

  TProfile* fvtxn_v2eta_east_docalib[NMUL];
  TProfile* fvtxn_v2eta_west_docalib[NMUL];
  TProfile* fvtxn_v2eta_both_docalib[NMUL];

  TProfile* bbcs_v3eta_east_docalib[NMUL];
  TProfile* bbcs_v3eta_west_docalib[NMUL];
  TProfile* bbcs_v3eta_both_docalib[NMUL];

  TProfile* fvtxs_v3eta_east_docalib[NMUL];
  TProfile* fvtxs_v3eta_west_docalib[NMUL];
  TProfile* fvtxs_v3eta_both_docalib[NMUL];

  TProfile* fvtxn_v3eta_east_docalib[NMUL];
  TProfile* fvtxn_v3eta_west_docalib[NMUL];
  TProfile* fvtxn_v3eta_both_docalib[NMUL];

  for ( int ic = 0; ic < NMUL; ++ic )
  {
    bbcs_v2eta_both_docalib[ic] = new TProfile(Form("bbcs_v2eta_both_docalib_cent%d", ic), "", 32, -3.2, 3.2, -1.1, 1.1);
    bbcs_v2eta_east_docalib[ic] = new TProfile(Form("bbcs_v2eta_east_docalib_cent%d", ic), "", 32, -3.2, 3.2, -1.1, 1.1);
    bbcs_v2eta_west_docalib[ic] = new TProfile(Form("bbcs_v2eta_west_docalib_cent%d", ic), "", 32, -3.2, 3.2, -1.1, 1.1);

    fvtxs_v2eta_both_docalib[ic] = new TProfile(Form("fvtxs_v2eta_both_docalib_cent%d", ic), "", 32, -3.2, 3.2, -1.1, 1.1);
    fvtxs_v2eta_east_docalib[ic] = new TProfile(Form("fvtxs_v2eta_east_docalib_cent%d", ic), "", 32, -3.2, 3.2, -1.1, 1.1);
    fvtxs_v2eta_west_docalib[ic] = new TProfile(Form("fvtxs_v2eta_west_docalib_cent%d", ic), "", 32, -3.2, 3.2, -1.1, 1.1);

    fvtxn_v2eta_both_docalib[ic] = new TProfile(Form("fvtxn_v2eta_both_docalib_cent%d", ic), "", 32, -3.2, 3.2, -1.1, 1.1);
    fvtxn_v2eta_east_docalib[ic] = new TProfile(Form("fvtxn_v2eta_east_docalib_cent%d", ic), "", 32, -3.2, 3.2, -1.1, 1.1);
    fvtxn_v2eta_west_docalib[ic] = new TProfile(Form("fvtxn_v2eta_west_docalib_cent%d", ic), "", 32, -3.2, 3.2, -1.1, 1.1);

    bbcs_v3eta_both_docalib[ic] = new TProfile(Form("bbcs_v3eta_both_docalib_cent%d", ic), "", 32, -3.2, 3.2, -1.1, 1.1);
    bbcs_v3eta_east_docalib[ic] = new TProfile(Form("bbcs_v3eta_east_docalib_cent%d", ic), "", 32, -3.2, 3.2, -1.1, 1.1);
    bbcs_v3eta_west_docalib[ic] = new TProfile(Form("bbcs_v3eta_west_docalib_cent%d", ic), "", 32, -3.2, 3.2, -1.1, 1.1);

    fvtxs_v3eta_both_docalib[ic] = new TProfile(Form("fvtxs_v3eta_both_docalib_cent%d", ic), "", 32, -3.2, 3.2, -1.1, 1.1);
    fvtxs_v3eta_east_docalib[ic] = new TProfile(Form("fvtxs_v3eta_east_docalib_cent%d", ic), "", 32, -3.2, 3.2, -1.1, 1.1);
    fvtxs_v3eta_west_docalib[ic] = new TProfile(Form("fvtxs_v3eta_west_docalib_cent%d", ic), "", 32, -3.2, 3.2, -1.1, 1.1);

    fvtxn_v3eta_both_docalib[ic] = new TProfile(Form("fvtxn_v3eta_both_docalib_cent%d", ic), "", 32, -3.2, 3.2, -1.1, 1.1);
    fvtxn_v3eta_east_docalib[ic] = new TProfile(Form("fvtxn_v3eta_east_docalib_cent%d", ic), "", 32, -3.2, 3.2, -1.1, 1.1);
    fvtxn_v3eta_west_docalib[ic] = new TProfile(Form("fvtxn_v3eta_west_docalib_cent%d", ic), "", 32, -3.2, 3.2, -1.1, 1.1);
  } // ic

  // ---------------------------------------------------------------------------------------------------------


  //------------------------------------------------------------//
  //  Finished initializing histograms                          //
  //------------------------------------------------------------//



  //------------------------------------------------------------//
  //                                                            //
  //               Initializing Tree Variables                  //
  //                                                            //
  //------------------------------------------------------------//

  cout << "Now getting ready to read in the tree branch addresses and stuff...." << endl;



  //tree variables
  float        event;
  float        d_bbcz;    // bbcz
  float        centrality; // integer but stored as float in PHGlobal etc
  float        frac;
  float        bbc_qn;
  float        bbc_qs;
  unsigned int trigger_scaled;
  unsigned int trigger_live;
  float        bc_x;
  float        bc_y;
  float        vtx_z;
  float        eventfvtx_x;
  float        eventfvtx_y;
  float        eventfvtx_z;
  float        d_Qx[9];
  float        d_Qy[9];
  float        d_Qw[9];
  float        d_BBC_charge[128];

  int          npc1;
  int          d_nFVTX_clus;
  int          d_nFVTXN_clus;
  int          d_nFVTXS_clus;
  float        d_FVTX_x[max_nf];
  float        d_FVTX_y[max_nf];
  float        d_FVTX_z[max_nf];

  int          d_ntrk;
  float        d_px[max_nh];
  float        d_py[max_nh];
  float        d_pz[max_nh];
  float        d_charge[max_nh];
  float        d_pc3sdz[max_nh];
  float        d_pc3sdphi[max_nh];

  int nfvtxt;
  int fnhits[75];
  float feta[75];
  float fphi[75];
  float fchisq[75];
  float fdcax[75];
  float fdcay[75];

  // List of branches
  TBranch* b_event;   //!
  TBranch* b_bbc_z;   //!
  TBranch* b_centrality;   //!
  TBranch* b_frac;   //!
  TBranch* b_bbc_qn;   //!
  TBranch* b_bbc_qs;   //!
  TBranch* b_npc1;   //!
  TBranch* b_trigger_scaled;   //!
  TBranch* b_trigger_live;   //!
  TBranch* b_d_Qx;   //!
  TBranch* b_d_Qy;   //!
  TBranch* b_d_Qw;   //!
  TBranch* b_bc_x;   //!
  TBranch* b_bc_y;   //!
  TBranch* b_vtx_z;   //!
  TBranch* b_fvtx_x;   //!
  TBranch* b_fvtx_y;   //!
  TBranch* b_fvtx_z;   //!
  TBranch* b_d_BBC_charge;   //!
  TBranch* b_d_nFVTX_clus;   //!
  TBranch* b_d_nFVTXN_clus;   //!
  TBranch* b_d_nFVTXS_clus;   //!
  TBranch* b_d_FVTX_x;   //!
  TBranch* b_d_FVTX_y;   //!
  TBranch* b_d_FVTX_z;   //!
  TBranch* b_ntrk;   //!
  TBranch* b_px;   //!
  TBranch* b_py;   //!
  TBranch* b_pz;   //!
  TBranch* b_charge;   //!
  TBranch* b_pc3sdz;   //!
  TBranch* b_pc3sdphi;   //!
  TBranch* b_nfvtxt;   //!
  TBranch* b_fnhits;   //!
  TBranch* b_fphi;   //!
  TBranch* b_feta;   //!
  TBranch* b_fchisq;   //!
  TBranch* b_fdcax;   //!
  TBranch* b_fdcay;   //!
  // TBranch* b_d_ntrk;   //!
  // TBranch* b_d_cntpx;   //!
  // TBranch* b_d_cntpy;   //!
  // TBranch* b_d_cntpz;   //!

  ntp_event_chain->SetBranchAddress("bbc_z", &d_bbcz, &b_bbc_z);
  ntp_event_chain->SetBranchAddress("centrality", &centrality, &b_centrality);
  ntp_event_chain->SetBranchAddress("frac", &frac, &b_frac);
  ntp_event_chain->SetBranchAddress("bbc_qn", &bbc_qn, &b_bbc_qn);
  ntp_event_chain->SetBranchAddress("bbc_qs", &bbc_qs, &b_bbc_qs);
  ntp_event_chain->SetBranchAddress("npc1", &npc1, &b_npc1);
  ntp_event_chain->SetBranchAddress("event", &event, &b_event);
  ntp_event_chain->SetBranchAddress("trigger_scaled", &trigger_scaled, &b_trigger_scaled);
  ntp_event_chain->SetBranchAddress("trigger_live", &trigger_live, &b_trigger_live);
  ntp_event_chain->SetBranchAddress("bc_x", &bc_x, &b_bc_x);
  ntp_event_chain->SetBranchAddress("bc_y", &bc_y, &b_bc_y);
  ntp_event_chain->SetBranchAddress("vtx_z", &vtx_z, &b_vtx_z);

  ntp_event_chain->SetBranchAddress("fvtx_x", &eventfvtx_x, &b_fvtx_x);
  ntp_event_chain->SetBranchAddress("fvtx_y", &eventfvtx_y, &b_fvtx_y);
  ntp_event_chain->SetBranchAddress("fvtx_z", &eventfvtx_z, &b_fvtx_z);

  ntp_event_chain->SetBranchAddress("d_BBC_charge", d_BBC_charge, &b_d_BBC_charge);
  ntp_event_chain->SetBranchAddress("d_Qx", d_Qx, &b_d_Qx);
  ntp_event_chain->SetBranchAddress("d_Qy", d_Qy, &b_d_Qy);
  ntp_event_chain->SetBranchAddress("d_Qw", d_Qw, &b_d_Qw);

  ntp_event_chain->SetBranchAddress("d_nFVTX_clus", &d_nFVTX_clus, &b_d_nFVTX_clus);
  ntp_event_chain->SetBranchAddress("d_nFVTXN_clus", &d_nFVTXN_clus, &b_d_nFVTXN_clus);
  ntp_event_chain->SetBranchAddress("d_nFVTXS_clus", &d_nFVTXS_clus, &b_d_nFVTXS_clus);
  ntp_event_chain->SetBranchAddress("d_FVTX_x", d_FVTX_x, &b_d_FVTX_x);
  ntp_event_chain->SetBranchAddress("d_FVTX_y", d_FVTX_y, &b_d_FVTX_y);
  ntp_event_chain->SetBranchAddress("d_FVTX_z", d_FVTX_z, &b_d_FVTX_z);

  ntp_event_chain->SetBranchAddress("d_ntrk", &d_ntrk, &b_ntrk);
  ntp_event_chain->SetBranchAddress("d_cntpx", d_px, &b_px);
  ntp_event_chain->SetBranchAddress("d_cntpy", d_py, &b_py);
  ntp_event_chain->SetBranchAddress("d_cntpz", d_pz, &b_pz);
  ntp_event_chain->SetBranchAddress("d_cntcharge", d_charge, &b_charge);
  ntp_event_chain->SetBranchAddress("d_cntpc3sdz", d_pc3sdz, &b_pc3sdz);
  ntp_event_chain->SetBranchAddress("d_cntpc3sdphi", d_pc3sdphi, &b_pc3sdphi);

  ntp_event_chain->SetBranchAddress("ntracklets", &nfvtxt, &b_nfvtxt);
  ntp_event_chain->SetBranchAddress("fnhits", fnhits, &b_fnhits);
  ntp_event_chain->SetBranchAddress("fphi", fphi, &b_fphi);
  ntp_event_chain->SetBranchAddress("feta", feta, &b_feta);
  ntp_event_chain->SetBranchAddress("fchisq", fchisq, &b_fchisq);
  ntp_event_chain->SetBranchAddress("fDCA_X", fdcax, &b_fdcax);
  ntp_event_chain->SetBranchAddress("fDCA_Y", fdcay, &b_fdcay);




  //------------------------------------------------------------//
  //------------------------------------------------------------//
  //                                                            //
  //                   Looping Over Event Tree                  //
  //                                                            //
  //------------------------------------------------------------//
  //------------------------------------------------------------//

  int event_counter = 0;

  int all_counter = 0;
  int bad_cent_counter = 0;

  int bad_vertex_counter = 0;

  Long64_t cluster_counter = 0;
  Long64_t gapcut_counter = 0;

  cout << "starting loop over events in the tree" << endl;

  int nentries = ntp_event_chain->GetEntries();
  cout << "total events = " << nentries << endl;

  for ( int ievt = 0; ievt < nentries; ++ievt )
  {

    // --- break and continue statements should happen much, much earlier --------------------
    if (rp_recal_pass < 1 || rp_recal_pass > 3) break; // rp_recal_pass only valid between 1 and 3

    //if ( ievt >= 100000 ) break; // just 100k events for testing, runs a little on the slow side...
    //if ( ievt >= 10000 ) break; // just 100k events for testing, runs a little on the slow side...
    ++all_counter;

    bool say_event = ( ievt % 1000 == 0 );

    if ( say_event ) cout << "event number = " << ievt << endl;

    //if ( ( say_event && verbosity > 0 ) || verbosity > 4 ) cout << "getting event level variables" << endl;
    ntp_event_chain->GetEntry(ievt);
    //if ( ( say_event && verbosity > 0 ) || verbosity > 4 ) cout << "Finished getting tree variables" << endl;

    // ---------------------
    // --- trigger selection
    // ---------------------
    unsigned int trigger_FVTXNSBBCScentral = 0x00100000;
    unsigned int trigger_FVTXNSBBCS        = 0x00400000;
    unsigned int trigger_BBCLL1narrowcent  = 0x00000008;
    unsigned int trigger_BBCLL1narrow      = 0x00000010;

    unsigned int accepted_triggers = 0;
    if ( energyflag == 200 ) accepted_triggers = trigger_BBCLL1narrowcent  | trigger_BBCLL1narrow;
    if ( energyflag == 62  ) accepted_triggers = trigger_BBCLL1narrowcent  | trigger_BBCLL1narrow;
    if ( energyflag == 20  ) accepted_triggers = trigger_FVTXNSBBCScentral | trigger_FVTXNSBBCS;
    if ( energyflag == 39  ) accepted_triggers = trigger_FVTXNSBBCScentral | trigger_FVTXNSBBCS;

    unsigned int passes_trigger = trigger_scaled & accepted_triggers;

    if ( passes_trigger == 0 )
    {
      if ( verbosity > 1 ) cout << "trigger rejected" << endl;
      continue;
    }

    // ---------------------
    // --- centrality
    // ---------------------
    int icent = -1;
    if ( centrality > -1 )
    {
      if ( energyflag == 200 || energyflag == 62 )
      {
        if ( centrality <= 5 ) icent = 0;
        else if ( centrality <= 10 ) icent = 1;
        else if ( centrality <= 20 ) icent = 2;
        else if ( centrality <= 40 ) icent = 3;
        else if ( centrality <= 60 ) icent = 4;
        else if ( centrality <= 88 ) icent = 5;
      }
      if ( energyflag == 39  )
      {
        if ( centrality <= 10 ) icent = 0;
        else if ( centrality <= 20 ) icent = 1;
        else if ( centrality <= 40 ) icent = 2;
        else if ( centrality <= 60 ) icent = 3;
        else if ( centrality <= 88 ) icent = 4;
      }
      if ( energyflag == 20  )
      {
        if ( centrality <= 20 ) icent = 0;
        else if ( centrality <= 40 ) icent = 1;
        else if ( centrality <= 60 ) icent = 2;
        else if ( centrality <= 88 ) icent = 3;
      }
    }
    else
    {
      if ( verbosity > 1 ) cout << "centrality undefined, skipping event" << endl;
      ++bad_cent_counter;
      continue;
    }
    if ( icent < 0 || icent > 5 )
    {
      if ( verbosity > 0 ) cout << "Well this is embarrassing.  icent is " << icent << " because centrality is " << centrality << endl;
      ++bad_cent_counter;
      continue;
    }

    // ---------------------
    // --- pile up rejection
    // ---------------------
    // reject pile up (double interaction) events based on the
    // frac value (< 0.95)
    // Only do this for 0-20% dAu 200 GeV
    if ( energyflag == 200 && centrality <= 20 && frac < fracCut )
      continue;


    // ---------------------
    // --- vertex binning
    // ---------------------
    double ZVTX = -9999;
    // if ( runNumber >= 454774 && runNumber <= 456283 ) ZVTX = d_bbcz;
    // if ( runNumber >= 456652 && runNumber <= 458167 ) ZVTX = eventfvtx_z;
    if ( energyflag == 200 || energyflag == 62 ) ZVTX = d_bbcz;
    if ( energyflag == 39 || energyflag == 20 ) ZVTX = eventfvtx_z;
    if ( fabs(ZVTX) > 10.0 )
    {
      if ( verbosity > 1 ) cout << "vertex rejected" << endl;
      continue;
    }

    // make sure bin number doesn't exceed number of bins
    int izvtx = NZPS * (ZVTX + 10) / 20;
    int izvtx_jamie = izvtx + 20;
    if ( izvtx_jamie > 59 )
    {
      cout << "YOU'RE GONNA DIE" << endl;
      continue;
    }
    if ( izvtx < 0 || izvtx >= NZPS )
    {
      cout << "z vertex bin count problem!!!!" << endl;
      cout << "bbcz = " << d_bbcz << endl;
      cout << "fvtx_z = " << eventfvtx_z << endl;
      cout << "bin number is " << izvtx << endl;
      continue;
    }

    // ---------------------
    // --- max clusters
    // ---------------------
    int toomanyclusters = 9999;
    if ( energyflag == 200 ) toomanyclusters = 4000;
    if ( energyflag == 62  ) toomanyclusters = 4000;
    if ( energyflag == 20  ) toomanyclusters = 300;
    if ( energyflag == 39  ) toomanyclusters = 500;

    bool is_okaync = ( d_nFVTX_clus < toomanyclusters );
    if ( !is_okaync )
    {
      if ( verbosity > 1 ) cout << "too many clusters" << endl;
      continue;
    }


    //------------------------------------------------------------//
    //                Calculating Event Planes                    //
    //------------------------------------------------------------//

    //if ( ( say_event && verbosity > 0 ) || verbosity > 4 ) cout << "Calculating event planes" << endl;

    // --- all numbers from Darren 2016-06-23
    const float x_off = 0.3;
    const float y_off = 0.02;
    const float beam_angle = 0.001;
    float vtx_z = d_bbcz;
    if ( eventfvtx_z > -999 ) vtx_z = eventfvtx_z;
    float vtx_x = x_off + atan(beam_angle) * vtx_z;
    float vtx_y = 0.02;

    // --- radius cut using FVTX coordinates
    // if ( sqrt(pow(eventfvtx_x-vtx_x,2.0) +  pow(eventfvtx_y-vtx_y,2.0)) >= 0.15 )
    //   {
    //     if ( verbosity > 1 ) cout << "rejecting event due to radius cut" << endl;
    //     continue;
    //   }

    // -------------
    // --- BBC stuff
    // -------------

    float bbc_qx2 = 0;
    float bbc_qy2 = 0;
    float bbc_qx3 = 0;
    float bbc_qy3 = 0;
    float bbc_qx4 = 0;
    float bbc_qy4 = 0;
    float bbc_qw = 0;

    float bbc_nw_qw = 0;

    float bbcq_ring[5] = {0};

    //if ( ( say_event && verbosity > 0 ) || verbosity > 4 ) cout << "Looping over BBC stuff now" << endl;

    if ( bbc_pmts )
    {
      //cout << "starting tubes loop" << endl;
      //for(int ipmt = 64; ipmt < 128; ipmt++)
      for (int ipmt = 0; ipmt < 64; ipmt++)
      {
        float bbc_charge = d_BBC_charge[ipmt];
        //cout << ipmt << " " << bbc_charge << " " << bbc_nw_qw << endl;
        if ( bbc_charge <= 0 ) continue;

        float bbc_x      = d_pmt_x[ipmt] - vtx_x * 10; //pmt location in mm
        float bbc_y      = d_pmt_y[ipmt] - vtx_y * 10;
        float bbc_z      = d_pmt_z       - vtx_z * 10;

        // --- rotation
        bbc_x = bbc_z * sin(-beam_angle) + bbc_x * cos(-beam_angle);

        float phi = TMath::ATan2(bbc_y, bbc_x);

        int ring = get_pmt_layer(ipmt);
        // float tube = ipmt;
        float bbc_charge_corrected = bbc_charge * tube_gaincorrection[ipmt]; // tube by tube gain correction from above
        //float correction = tg_jamie_bbcs_zvtx[izvtx_jamie]->Eval(phi);
        //cout << "info from Jamie histogram is " << correction << " (and by the way the jamie vertex bin is " << izvtx_jamie << ")" << endl;
        // if ( correction != correction ) correction = 1;
        // if ( correction < 0 ) correction = 1;
        // if ( correction > 10 ) correction = 1;
        // float bbc_charge_corrected = bbc_charge * correction;
        if ( bbc_charge_corrected <= 0 ) continue;

        bbc_qx2 += bbc_charge_corrected * TMath::Cos(2 * phi);
        bbc_qy2 += bbc_charge_corrected * TMath::Sin(2 * phi);
        bbc_qx3 += bbc_charge_corrected * TMath::Cos(3 * phi);
        bbc_qy3 += bbc_charge_corrected * TMath::Sin(3 * phi);
        bbc_qx4 += bbc_charge_corrected * TMath::Cos(4 * phi);
        bbc_qy4 += bbc_charge_corrected * TMath::Sin(4 * phi);
        bbc_qw += bbc_charge_corrected;

        bbc_nw_qw += bbc_charge;

        if (ring >= 0 && ring < 5)
          bbcq_ring[ring] += bbc_charge;

        //cout << ipmt << " " << bbc_charge << " " << bbc_nw_qw << endl;
      } // loop over tubes
    } // check on tubes


    if ( say_event )
    {
      //cout << "bbc charge = " << bbc_qw << endl;
      cout << "bbc charge = " << bbc_nw_qw << endl;
      //cout << "bbc charge direct sum = " << bbc_qs << "+" << bbc_qn << "=" << bbc_qs+bbc_qn << endl;
      cout << "centrality = " << centrality << endl;
    }

    //cout << "HELLO HERE I AM" << endl;

    ++event_counter;

    th1d_BBC_charge->Fill(bbc_nw_qw);
    th1d_FVTX_nclus->Fill(d_nFVTX_clus);
    th1d_FVTX_ntrk->Fill(nfvtxt);
    th2d_qBBC_nFVTX->Fill(bbc_nw_qw, d_nFVTX_clus);
    // fill charge in each ring
    for (int ir = 0; ir < 5; ir++)
    {
      if ( bbc_nw_qw > 0 )
        tp1f_bbc_fcharge_ring->Fill(ir, bbcq_ring[ir] / bbc_nw_qw);
    }

    // --------------
    // --- FVTX stuff
    // --------------

    float fvtxs_qx2[10];//all layers then 0 1 2 3
    float fvtxs_qy2[10];
    float fvtxs_qx3[10];//all layers then 0 1 2 3
    float fvtxs_qy3[10];
    float fvtxs_qx4[10];//all layers then 0 1 2 3
    float fvtxs_qy4[10];
    float fvtxs_qw[10];

    for (int ilayer = 0; ilayer < 10; ilayer++)
    {
      fvtxs_qx2[ilayer] = 0.0;
      fvtxs_qy2[ilayer] = 0.0;
      fvtxs_qx3[ilayer] = 0.0;
      fvtxs_qy3[ilayer] = 0.0;
      fvtxs_qx4[ilayer] = 0.0;
      fvtxs_qy4[ilayer] = 0.0;
      fvtxs_qw[ilayer] = 0.0;
    } // loop over layers

    // --- now FVTX North

    float fvtxn_qx2[10];//all layers then 0 1 2 3
    float fvtxn_qy2[10];
    float fvtxn_qx3[10];//all layers then 0 1 2 3
    float fvtxn_qy3[10];
    float fvtxn_qx4[10];//all layers then 0 1 2 3
    float fvtxn_qy4[10];
    float fvtxn_qw[10];

    for (int ilayer = 0; ilayer < 10; ilayer++)
    {
      fvtxn_qx2[ilayer] = 0.0;
      fvtxn_qy2[ilayer] = 0.0;
      fvtxn_qx3[ilayer] = 0.0;
      fvtxn_qy3[ilayer] = 0.0;
      fvtxn_qx4[ilayer] = 0.0;
      fvtxn_qy4[ilayer] = 0.0;
      fvtxn_qw[ilayer] = 0.0;
    } // loop over layers


    //if ( ( say_event && verbosity > 0 ) || verbosity > 4 ) cout << "Looping over FVTX cluster" << endl;
    if ( fvtx_clusters )
    {
      for (int iclus = 0; iclus < d_nFVTX_clus; iclus++)
      {
        float fvtx_x_orig = d_FVTX_x[iclus] - vtx_x;
        float fvtx_y_orig = d_FVTX_y[iclus] - vtx_y;
        float fvtx_z_orig = d_FVTX_z[iclus] - vtx_z;

        // --- rotation
        float fvtx_x, fvtx_y, fvtx_z;
        fvtx_x = fvtx_z_orig * sin(-beam_angle) + fvtx_x_orig * cos(-beam_angle);
        fvtx_y = fvtx_y_orig;
        fvtx_z = -1 * fvtx_x_orig * sin(-beam_angle) + fvtx_z_orig * cos(-beam_angle);

        double fvtx_r = sqrt(pow(fvtx_x, 2.0) + pow(fvtx_y, 2.0));
        double fvtx_the = atan2(fvtx_r, fvtx_z);
        double fvtx_eta = -log(tan(0.5 * fvtx_the));
        float phi = TMath::ATan2(fvtx_y, fvtx_x);

        // --- determine layer based on cluster z
        int fvtx_layer = get_fvtx_layer(d_FVTX_z[iclus]); // raw z to get layer

        // --- gap cut removes events where fvtx z vertex is right below fvtx south
        int igap = (fabs(fvtx_eta) - 1.0) / 0.5;
        int id_fvtx = fvtx_layer * 5 + igap;
        ++cluster_counter;
        if (!(id_fvtx >= 0 && id_fvtx < 40))
        {
          if ( verbosity > 5 )
          {
            cout << "gap cut rejecting cluster??? cluster z = "
                 << d_FVTX_z[iclus] << " fvtx z = "
                 << eventfvtx_z << " bbc z = "
                 << d_bbcz << " eta = "
                 << fvtx_eta << " "
                 << endl;
          }
          ++gapcut_counter;
          continue;
        }
        // --------------------------------------


        // cluster radius cut for lower energies (USE ENERGY FLAG!!)
        double FVTX_r = sqrt(pow(d_FVTX_x[iclus], 2.0) + pow(d_FVTX_y[iclus], 2.0));
        if ( (energyflag == 20 || energyflag == 39) && FVTX_r < 5.2 ) continue;
        // --------------------------------------


        // --- south side
        if ( d_FVTX_z[iclus] < 0 )
        {
          fvtxs_qx2[0] += TMath::Cos(2 * phi);
          fvtxs_qy2[0] += TMath::Sin(2 * phi);
          fvtxs_qx3[0] += TMath::Cos(3 * phi);
          fvtxs_qy3[0] += TMath::Sin(3 * phi);
          fvtxs_qx4[0] += TMath::Cos(4 * phi);
          fvtxs_qy4[0] += TMath::Sin(4 * phi);

          fvtxs_qw[0] += 1;
        } // check on south

        // --- north side
        if ( d_FVTX_z[iclus] > 0 )
        {
          fvtxn_qx2[0] += TMath::Cos(2 * phi);
          fvtxn_qy2[0] += TMath::Sin(2 * phi);
          fvtxn_qx3[0] += TMath::Cos(3 * phi);
          fvtxn_qy3[0] += TMath::Sin(3 * phi);
          fvtxn_qx4[0] += TMath::Cos(4 * phi);
          fvtxn_qy4[0] += TMath::Sin(4 * phi);

          fvtxn_qw[0] += 1;
        } // check on north

      } // loop over cluster
    } // check on clusters

    th1d_FVTXS_nclus->Fill(fvtxs_qw[0]);
    th1d_FVTXN_nclus->Fill(fvtxn_qw[0]);

    // --- the array that has all of the Q vectors
    float sumxy[NHAR][NDETSHORT][4];
    for (int i = 0; i < NHAR; i++)
    {
      for (int j = 0; j < NDETSHORT; j++)
      {
        for (int k = 0; k < 4; k++) //qx qy qw psi
        {
          sumxy[i][j][k] = 0; // initialize to 0
        } // x,y,w,psi
      } // detectors
    } // harmonics


    if ( bbc_pmts )
    {
      sumxy[1][bbcs_index][0] = bbc_qx2;
      sumxy[1][bbcs_index][1] = bbc_qy2;
      sumxy[1][bbcs_index][2] = bbc_qw;
      sumxy[2][bbcs_index][0] = bbc_qx3;
      sumxy[2][bbcs_index][1] = bbc_qy3;
      sumxy[2][bbcs_index][2] = bbc_qw;
      sumxy[3][bbcs_index][0] = bbc_qx4;
      sumxy[3][bbcs_index][1] = bbc_qy4;
      sumxy[3][bbcs_index][2] = bbc_qw;
    }

    if ( fvtx_clusters )
    {
      sumxy[1][fvtxs_index][0] = fvtxs_qx2[0];
      sumxy[1][fvtxs_index][1] = fvtxs_qy2[0];
      sumxy[1][fvtxs_index][2] = fvtxs_qw[0];
      sumxy[1][fvtxn_index][0] = fvtxn_qx2[0];
      sumxy[1][fvtxn_index][1] = fvtxn_qy2[0];
      sumxy[1][fvtxn_index][2] = fvtxn_qw[0];
      sumxy[2][fvtxs_index][0] = fvtxs_qx3[0];
      sumxy[2][fvtxs_index][1] = fvtxs_qy3[0];
      sumxy[2][fvtxs_index][2] = fvtxs_qw[0];
      sumxy[2][fvtxn_index][0] = fvtxn_qx3[0];
      sumxy[2][fvtxn_index][1] = fvtxn_qy3[0];
      sumxy[2][fvtxn_index][2] = fvtxn_qw[0];
      sumxy[3][fvtxs_index][0] = fvtxs_qx4[0];
      sumxy[3][fvtxs_index][1] = fvtxs_qy4[0];
      sumxy[3][fvtxs_index][2] = fvtxs_qw[0];
      sumxy[3][fvtxn_index][0] = fvtxn_qx4[0];
      sumxy[3][fvtxn_index][1] = fvtxn_qy4[0];
      sumxy[3][fvtxn_index][2] = fvtxn_qw[0];
    } // check on clusters


    if ( DIAG )
    {
      cout << "bbc from node tree: " << d_Qx[5] << " " << d_Qy[5] << " " << d_Qw[5] << endl;
      cout << "bbc from me: " << bbc_qx2 << " " << bbc_qy2 << " " << bbc_qw << endl;

      cout << "fvtx raw: " << endl;
      cout << "from node tree: " << d_Qx[4] << " " << d_Qy[4] << " " << d_Qw[4] << endl;
      cout << "from clusters: " << fvtxs_qx2[0] << " " << fvtxs_qy2[0] << " " << fvtxs_qw[0] << endl;
    }

    for ( int ih = 1; ih < NHAR; ih++ )
    {
      for (int id = 0; id < NDETSHORT; id++ )
      {
        if (sumxy[ih][id][2] > 0)
        {
          //float psi = atan2(sumxy[ih][id][1],sumxy[ih][id][0])/2.0;
          float psi = atan2(sumxy[ih][id][1], sumxy[ih][id][0]) / float(ih + 1);
          if ( DIAG ) cout << "RAW: for id: " << id << " psi: " << psi << endl;
          psi_bf[icent][ih][id]->Fill(izvtx, psi);
        } // check on weight
      } // detectors
    } // harmonics

    //------------------------------------------------------------//
    //                Flattening iteration                        //
    //------------------------------------------------------------//
    //int icent = 0;
    for ( int ih = 1; ih < NHAR; ih++ )
    {
      for ( int id = 0; id < NDETSHORT; id++ )
      {
        if ( sumxy[ih][id][2] > 0.0 )
        {
          sumxy[ih][id][3] = atan2(sumxy[ih][id][1], sumxy[ih][id][0]) / (ih + 1.0);
        }
        if ( sumxy[ih][id][2] > 0.0 ) // check on weight (x,y,w,psi)
        {
          for ( int ib = 0; ib < 2; ib++ )
          {
            sumxy[ih][id][ib] /= sumxy[ih][id][2]; // normalize to the weight

            //if(ih==1 && id==0 && ib==0 && sumxy[ih][id][ib]>1) cout<<sumxy[ih][id][ib]<<endl;
            if ( rp_recal_pass > 0 )
            {
              ave[icent][izvtx][ih][id]->Fill(ib + 0.0, sumxy[ih][id][ib]);
              if (id == 0 && DIAG) cout << "filled ave: " << ih << " " << id << " " << ib << " with: " << sumxy[ih][id][ib] << endl;
            } // pass > 0
            float sxy = sumxy[ih][id][ib];
            float mxy = mean[icent][izvtx][ih][id][ib]; // for recentering qx and qy (???)
            float wxy = widt[icent][izvtx][ih][id][ib]; // for recentering qx and qy (???)

            //if(ic==0 && izvtx==0 && ih==1 && id==0) cout<<ib<<" "<<sxy<<" "<<mxy<<" "<<wxy<<endl;
            sumxy[ih][id][ib] = (sxy - mxy) / wxy; // recentered by mean and renormalized to width
            if ( rp_recal_pass > 0 )
            {
              ave[icent][izvtx][ih][id]->Fill(ib + 2.0, sumxy[ih][id][ib]); // ib+2 to avoid overlap
              if (id == 0 && DIAG) cout << "filled ave2: " << ih << " " << id << " " << ib << " with: " << sumxy[ih][id][ib] << endl;
            } // pass > 0
          } // if weight > 0

          sumxy[ih][id][3] = atan2(sumxy[ih][id][1], sumxy[ih][id][0]) / (ih + 1.0);
          if ( rp_recal_pass > 0 )
          {
            // my own simpler version of the above histogram
            psi_mf[icent][ih][id]->Fill(izvtx, sumxy[ih][id][3]);
          }

          float psi = sumxy[ih][id][3] * (ih + 1.0);
          if ( ih == 1 && id == 0 && DIAG )  cout << "psi-1 bbc: " << psi << endl;
          float dp = 0.0;
          // --- flattening part, fourier components of psi distribution
          for (int io = 0; io < NORD; io++)
          {
            float cc = cos((io + 1.0) * psi);
            float ss = sin((io + 1.0) * psi);
            // first set of fourier components of psi
            if (rp_recal_pass > 0) flt[icent][izvtx][ih][id]->Fill(io + 0.0, cc);
            if (rp_recal_pass > 0) flt[icent][izvtx][ih][id]->Fill(io + NORD, ss);
            // --- four means fourier
            float aa = four[icent][izvtx][ih][id][0][io]; // mean cos
            float bb = four[icent][izvtx][ih][id][1][io]; // mean sin
            // dp is offset to psi, aa and bb are zero in first pass, non zero later
            dp += (aa * ss - bb * cc) * 2.0 / (io + 1.0); // ( trig identity cos(A+B) = cosAsinB - cosBsinA )
          } // orders
          psi += dp; // shift psi by...
          psi = atan2(sin(psi), cos(psi)); // trick to readjust the range
          if ( ih == 1 && id == 0 && DIAG )  cout << "psi-2 bbc: " << psi << endl;
          for (int io = 0; io < NORD; io++)
          {
            float cc = cos((io + 1.0) * psi);
            float ss = sin((io + 1.0) * psi);
            // --- fourier components of modified psi
            if (rp_recal_pass > 0) flt[icent][izvtx][ih][id]->Fill(io + NORD * 2.0, cc);
            if (rp_recal_pass > 0) flt[icent][izvtx][ih][id]->Fill(io + NORD * 3.0, ss);
          }
          sumxy[ih][id][3] = psi / (ih + 1.0);
        } // end if weight > 0
        else
        {
          sumxy[ih][id][3] = -9999.9;
        } // otherwise set psi to some crazy number
      } // detectors
    } // harmonics



    // testing Qx/Qy offset on FVTXS Psi2
    sumxy[1][fvtxs_index][0] = cos(sumxy[1][fvtxs_index][3] * 2) + qxOffset;
    sumxy[1][fvtxs_index][1] = sin(sumxy[1][fvtxs_index][3] * 2) + qyOffset[eidx][icent];
    sumxy[1][fvtxs_index][3] = atan2(sumxy[1][fvtxs_index][1], sumxy[1][fvtxs_index][0]) / (1 + 1.0);




    if ( DIAG ) cout << "bbc_rp2: " << sumxy[1][0][3] << endl;

    for ( int ih = 1; ih < NHAR; ih++ )
    {
      for ( int id = 0; id < NDETSHORT; id++ )
      {
        if ( sumxy[ih][id][2] > 0 )
        {
          psi_af[icent][ih][id]->Fill(izvtx, sumxy[ih][id][3]);
          if ( DIAG ) cout << "CORR: for id: " << id << " psi: " << sumxy[ih][id][3] << endl;
        } // check on weight
      } // detectors
    } // harmonics


    // ---
    // --- now going to calculate v2
    // ---
    if ( rp_recal_pass < 3 ) continue; // don't calculate v2 except for final pass



    // --------------------------------------
    // --- now the standard event plane stuff
    // --------------------------------------

    // ---
    // --- south
    // ---

    float bbc_south_psi2_docalib;
    float bbc_south_psi3_docalib;

    float fvtx_south_psi2_docalib;
    float fvtx_south_psi3_docalib;

    float fvtx_north_psi2_docalib;
    float fvtx_north_psi3_docalib;

    if ( bbc_pmts )
    {
      bbc_south_psi2_docalib = (sumxy[1][bbcs_index][2] > 0) ? sumxy[1][bbcs_index][3] : -9999.9;
      bbc_south_psi3_docalib = (sumxy[2][bbcs_index][2] > 0) ? sumxy[2][bbcs_index][3] : -9999.9;
    }
    if ( fvtx_clusters )
    {
      // ---
      fvtx_south_psi2_docalib  = (sumxy[1][fvtxs_index][2] > 4)  ? sumxy[1][fvtxs_index][3]  : -9999.9;
      fvtx_south_psi3_docalib  = (sumxy[2][fvtxs_index][2] > 4)  ? sumxy[2][fvtxs_index][3]  : -9999.9;

      // ---
      fvtx_north_psi2_docalib  = (sumxy[1][fvtxn_index][2] > 4)  ? sumxy[1][fvtxn_index][3]  : -9999.9;
      fvtx_north_psi3_docalib  = (sumxy[2][fvtxn_index][2] > 4)  ? sumxy[2][fvtxn_index][3]  : -9999.9;
    }


    if ( verbosity > 0 && ( fvtx_north_psi2_docalib < -999 || fvtx_south_psi2_docalib < -999 || bbc_south_psi2_docalib < -999 ) )
    {
      cout << "POSSIBLE ISSUE WITH EVENT PLANES!!!  ONE OR MORE IS -9999" << endl;
      cout << "BBC south event plane " << bbc_south_psi2_docalib << endl;
      cout << "FVTX south event plane " << fvtx_south_psi2_docalib << endl;
      cout << "FVTX north event plane " << fvtx_north_psi2_docalib << endl;
      cout << "BBC charge " << bbc_qw << endl;
      cout << "centrality " << centrality << endl;
    }



    // ---
    // --- resolution histograms
    // ---

    // --- BBC and FVTX south
    tp1f_reso2_BBC_FVTX[icent]->Fill(0.0, cos(2 * (bbc_south_psi2_docalib - fvtx_south_psi2_docalib)));
    tp1f_reso3_BBC_FVTX[icent]->Fill(0.0, cos(3 * (bbc_south_psi3_docalib - fvtx_south_psi3_docalib)));

    // --- BBC and FVTX north
    tp1f_reso2_BBC_FVTXN[icent]->Fill(0.0, cos(2 * (bbc_south_psi2_docalib - fvtx_north_psi2_docalib)));
    tp1f_reso3_BBC_FVTXN[icent]->Fill(0.0, cos(3 * (bbc_south_psi3_docalib - fvtx_north_psi3_docalib)));

    // --- FVTX south and north
    tp1f_reso2_FVTXS_FVTXN[icent]->Fill(0.0, cos(2 * (fvtx_south_psi2_docalib - fvtx_north_psi2_docalib)));
    tp1f_reso3_FVTXS_FVTXN[icent]->Fill(0.0, cos(3 * (fvtx_south_psi3_docalib - fvtx_north_psi3_docalib)));

    // ----------------------------------------------------------------------------

    //start of cnt track loop
    if ( cnt_tracks )
    {
      for (int itrk = 0; itrk < d_ntrk; itrk++)
      {
        float px    = d_px[itrk];
        float py    = d_py[itrk];
        float pz    = d_pz[itrk];
        // float charge    = d_charge[itrk];
        // float pc3sdz    = d_pc3sdz[itrk];
        // float pc3sdphi  = d_pc3sdphi[itrk];
        // if ( fabs(pc3sdz) > 2.0 || fabs(pc3sdphi) > 2.0 ) continue;

        int dcarm = 0;
        if (px > 0) dcarm = 1;

        float phi0 = TMath::ATan2(py, px);
        float pt = sqrt(px * px + py * py);

        if ( pt < 0.2 || pt > 5.0 ) continue; // pt cut added 2016-06-30

        // -------------------------------------------------------
        // --- finished with nodetree part, now doing docalib part
        // -------------------------------------------------------

        // --- rotation now done in trees
        float phi_angle = phi0;
        float pt_angle = pt;
        double theta = atan2(pt, pz);
        double eta = -log(tan(theta / 2));
        //cout << "pz is " << pz << " and eta is " << eta << endl;


        // BBC event plane
        if ( bbc_pmts )
        {
          if (-4.0 < bbc_south_psi2_docalib && bbc_south_psi2_docalib < 4.0)
          {
            // --- 2nd harmonic
            double bbc_dphi2_docalib = phi_angle - bbc_south_psi2_docalib;
            double cosbbc_dphi2_docalib = TMath::Cos(2 * bbc_dphi2_docalib);

            bbcs_v2_both_docalib[icent]->Fill(pt_angle, cosbbc_dphi2_docalib);
            if ( dcarm == 1 ) bbcs_v2_west_docalib[icent]->Fill(pt_angle, cosbbc_dphi2_docalib);
            if ( dcarm == 0 ) bbcs_v2_east_docalib[icent]->Fill(pt_angle, cosbbc_dphi2_docalib);

            bbcs_v2eta_both_docalib[icent]->Fill(eta, cosbbc_dphi2_docalib);
            if ( dcarm == 1 ) bbcs_v2eta_west_docalib[icent]->Fill(eta, cosbbc_dphi2_docalib);
            if ( dcarm == 0 ) bbcs_v2eta_east_docalib[icent]->Fill(eta, cosbbc_dphi2_docalib);


            // --- 3rd harmonic
            double bbc_dphi3_docalib = phi_angle - bbc_south_psi3_docalib;
            double cosbbc_dphi3_docalib = TMath::Cos(3 * bbc_dphi3_docalib);

            bbcs_v3_both_docalib[icent]->Fill(pt_angle, cosbbc_dphi3_docalib);
            if ( dcarm == 1 ) bbcs_v3_west_docalib[icent]->Fill(pt_angle, cosbbc_dphi3_docalib);
            if ( dcarm == 0 ) bbcs_v3_east_docalib[icent]->Fill(pt_angle, cosbbc_dphi3_docalib);

            bbcs_v3eta_both_docalib[icent]->Fill(eta, cosbbc_dphi3_docalib);
            if ( dcarm == 1 ) bbcs_v3eta_west_docalib[icent]->Fill(eta, cosbbc_dphi3_docalib);
            if ( dcarm == 0 ) bbcs_v3eta_east_docalib[icent]->Fill(eta, cosbbc_dphi3_docalib);


            // --- resolutions
            if ( pt_angle > 0.4 && pt_angle < 3.0 )
            {
              tp1f_reso2_BBC_CNT[icent]->Fill(0.0, cosbbc_dphi2_docalib);
              tp1f_reso3_BBC_CNT[icent]->Fill(0.0, cosbbc_dphi3_docalib);
            }
          } // check on bbc EP
        } // check on tubes

        // FVTX Event Plane
        if ( fvtx_clusters )
        {

          // south all layers
          if ( -4.0 < fvtx_south_psi2_docalib && fvtx_south_psi2_docalib < 4.0 )
          {

            // --- 2nd harmonic
            double fvtx_dphi2_docalib = phi_angle - fvtx_south_psi2_docalib;
            double cosfvtx_dphi2_docalib = TMath::Cos(2 * fvtx_dphi2_docalib);

            fvtxs_v2_both_docalib[icent]->Fill(pt_angle, cosfvtx_dphi2_docalib);
            if ( dcarm == 1 ) fvtxs_v2_west_docalib[icent]->Fill(pt_angle, cosfvtx_dphi2_docalib);
            if ( dcarm == 0 ) fvtxs_v2_east_docalib[icent]->Fill(pt_angle, cosfvtx_dphi2_docalib);

            fvtxs_v2eta_both_docalib[icent]->Fill(eta, cosfvtx_dphi2_docalib);
            if ( dcarm == 1 ) fvtxs_v2eta_west_docalib[icent]->Fill(eta, cosfvtx_dphi2_docalib);
            if ( dcarm == 0 ) fvtxs_v2eta_east_docalib[icent]->Fill(eta, cosfvtx_dphi2_docalib);

            fvtxs_v2_both_cosphi[icent]->Fill(pt_angle, cos(2 * phi_angle));
            if ( dcarm == 1 ) fvtxs_v2_west_cosphi[icent]->Fill(pt_angle, cos(2 * phi_angle));
            if ( dcarm == 0 ) fvtxs_v2_east_cosphi[icent]->Fill(pt_angle, cos(2 * phi_angle));

            fvtxs_v2_both_sinphi[icent]->Fill(pt_angle, sin(2 * phi_angle));
            if ( dcarm == 1 ) fvtxs_v2_west_sinphi[icent]->Fill(pt_angle, sin(2 * phi_angle));
            if ( dcarm == 0 ) fvtxs_v2_east_sinphi[icent]->Fill(pt_angle, sin(2 * phi_angle));

            // --- 3rd harmonic
            double fvtx_dphi3_docalib = phi_angle - fvtx_south_psi3_docalib;
            double cosfvtx_dphi3_docalib = TMath::Cos(3 * fvtx_dphi3_docalib);

            fvtxs_v3_both_docalib[icent]->Fill(pt_angle, cosfvtx_dphi3_docalib);
            if ( dcarm == 1 ) fvtxs_v3_west_docalib[icent]->Fill(pt_angle, cosfvtx_dphi3_docalib);
            if ( dcarm == 0 ) fvtxs_v3_east_docalib[icent]->Fill(pt_angle, cosfvtx_dphi3_docalib);

            fvtxs_v3eta_both_docalib[icent]->Fill(eta, cosfvtx_dphi3_docalib);
            if ( dcarm == 1 ) fvtxs_v3eta_west_docalib[icent]->Fill(eta, cosfvtx_dphi3_docalib);
            if ( dcarm == 0 ) fvtxs_v3eta_east_docalib[icent]->Fill(eta, cosfvtx_dphi3_docalib);


            // --- ep resolutions
            if ( pt_angle > 0.4 && pt_angle < 3.0 )
            {
              tp1f_reso2_CNT_FVTX[icent]->Fill(0.0, cosfvtx_dphi2_docalib);
              tp1f_reso3_CNT_FVTX[icent]->Fill(0.0, cosfvtx_dphi3_docalib);
            }

          } // check on ep south all layers


          // north all layers
          if ( -4.0 < fvtx_north_psi2_docalib && fvtx_north_psi2_docalib < 4.0 )
          {

            // --- 2nd harmonic
            double fvtxn_dphi2_docalib = phi_angle - fvtx_north_psi2_docalib;
            double cosfvtxn_dphi2_docalib = TMath::Cos(2 * fvtxn_dphi2_docalib);

            fvtxn_v2_both_docalib[icent]->Fill(pt_angle, cosfvtxn_dphi2_docalib);
            if ( dcarm == 1 ) fvtxn_v2_west_docalib[icent]->Fill(pt_angle, cosfvtxn_dphi2_docalib);
            if ( dcarm == 0 ) fvtxn_v2_east_docalib[icent]->Fill(pt_angle, cosfvtxn_dphi2_docalib);

            fvtxn_v2eta_both_docalib[icent]->Fill(eta, cosfvtxn_dphi2_docalib);
            if ( dcarm == 1 ) fvtxn_v2eta_west_docalib[icent]->Fill(eta, cosfvtxn_dphi2_docalib);
            if ( dcarm == 0 ) fvtxn_v2eta_east_docalib[icent]->Fill(eta, cosfvtxn_dphi2_docalib);


            // --- 3rd harmonic
            double fvtxn_dphi3_docalib = phi_angle - fvtx_north_psi3_docalib;
            double cosfvtxn_dphi3_docalib = TMath::Cos(3 * fvtxn_dphi3_docalib);

            fvtxn_v3_both_docalib[icent]->Fill(pt_angle, cosfvtxn_dphi3_docalib);
            if ( dcarm == 1 ) fvtxn_v3_west_docalib[icent]->Fill(pt_angle, cosfvtxn_dphi3_docalib);
            if ( dcarm == 0 ) fvtxn_v3_east_docalib[icent]->Fill(pt_angle, cosfvtxn_dphi3_docalib);

            fvtxn_v3eta_both_docalib[icent]->Fill(eta, cosfvtxn_dphi3_docalib);
            if ( dcarm == 1 ) fvtxn_v3eta_west_docalib[icent]->Fill(eta, cosfvtxn_dphi3_docalib);
            if ( dcarm == 0 ) fvtxn_v3eta_east_docalib[icent]->Fill(eta, cosfvtxn_dphi3_docalib);


            // --- ep resolutions
            if ( pt_angle > 0.4 && pt_angle < 3.0 )
            {
              tp1f_reso2_CNT_FVTXN[icent]->Fill(0.0, cosfvtxn_dphi2_docalib);
              tp1f_reso3_CNT_FVTXN[icent]->Fill(0.0, cosfvtxn_dphi3_docalib);
            }

          } // check on ep north all layers

        } // check on fvtx clusters

      } // loop over cnt tracks

    } // check on cnt tracks


    if ( fvtx_tracks )
    {
      // loop over fvtx tracks
      for ( int i = 0; i < nfvtxt; ++i )
      {
        float phi = fphi[i];
        float eta = feta[i];
        float dcax = fdcax[i];
        float dcay = fdcay[i];
        float chisq = fchisq[i];
        // if ( fabs(dcax) > 0.5 || fabs(dcay) > 0.5 ) continue;
        if ( fabs(dcax - x_off) > 2.0 || fabs(dcay - y_off) > 2.0 ) continue;
        if ( chisq > 5 ) continue;
        int ns = 0;
        int nn = 0;
        if ( eta > 0 ) nn = 1;
        else if ( eta < 0 ) ns = 1;
        // cout << "eta is " << eta << " and ns is " << ns << " and nn is " << nn << endl;

        bool west = (phi > -pi / 2. && phi < pi / 2.);

        // BBC event plane
        if ( bbc_pmts )
        {
          if (-4.0 < bbc_south_psi2_docalib && bbc_south_psi2_docalib < 4.0)
          {
            // --- 2nd harmonic
            double bbc_dphi2_docalib = phi - bbc_south_psi2_docalib;
            double cosbbc_dphi2_docalib = TMath::Cos(2 * bbc_dphi2_docalib);

            bbcs_v2eta_both_docalib[icent]->Fill(eta, cosbbc_dphi2_docalib);
            if ( west ) bbcs_v2eta_west_docalib[icent]->Fill(eta, cosbbc_dphi2_docalib);
            if ( !west ) bbcs_v2eta_east_docalib[icent]->Fill(eta, cosbbc_dphi2_docalib);


            // --- 3rd harmonic
            double bbc_dphi3_docalib = phi - bbc_south_psi3_docalib;
            double cosbbc_dphi3_docalib = TMath::Cos(3 * bbc_dphi3_docalib);

            bbcs_v3eta_both_docalib[icent]->Fill(eta, cosbbc_dphi3_docalib);
            if ( west ) bbcs_v3eta_west_docalib[icent]->Fill(eta, cosbbc_dphi3_docalib);
            if ( !west ) bbcs_v3eta_east_docalib[icent]->Fill(eta, cosbbc_dphi3_docalib);

          } // check on bbc EP
        } // check on tubes

        // FVTX Event Plane
        if ( fvtx_clusters )
        {

          // south all layers
          if ( -4.0 < fvtx_south_psi2_docalib && fvtx_south_psi2_docalib < 4.0 )
          {

            // --- 2nd harmonic
            double fvtx_dphi2_docalib = phi - fvtx_south_psi2_docalib;
            double cosfvtx_dphi2_docalib = TMath::Cos(2 * fvtx_dphi2_docalib);

            fvtxs_v2eta_both_docalib[icent]->Fill(eta, cosfvtx_dphi2_docalib);
            if ( west ) fvtxs_v2eta_west_docalib[icent]->Fill(eta, cosfvtx_dphi2_docalib);
            if ( !west ) fvtxs_v2eta_east_docalib[icent]->Fill(eta, cosfvtx_dphi2_docalib);


            // --- 3rd harmonic
            double fvtx_dphi3_docalib = phi - fvtx_south_psi3_docalib;
            double cosfvtx_dphi3_docalib = TMath::Cos(3 * fvtx_dphi3_docalib);

            fvtxs_v3eta_both_docalib[icent]->Fill(eta, cosfvtx_dphi3_docalib);
            if ( west ) fvtxs_v3eta_west_docalib[icent]->Fill(eta, cosfvtx_dphi3_docalib);
            if ( !west ) fvtxs_v3eta_east_docalib[icent]->Fill(eta, cosfvtx_dphi3_docalib);

          } // check on ep south all layers


          // north all layers
          if ( -4.0 < fvtx_north_psi2_docalib && fvtx_north_psi2_docalib < 4.0 )
          {

            // --- 2nd harmonic
            double fvtxn_dphi2_docalib = phi - fvtx_north_psi2_docalib;
            double cosfvtxn_dphi2_docalib = TMath::Cos(2 * fvtxn_dphi2_docalib);

            fvtxn_v2eta_both_docalib[icent]->Fill(eta, cosfvtxn_dphi2_docalib);
            if ( west ) fvtxn_v2eta_west_docalib[icent]->Fill(eta, cosfvtxn_dphi2_docalib);
            if ( !west ) fvtxn_v2eta_east_docalib[icent]->Fill(eta, cosfvtxn_dphi2_docalib);


            // --- 3rd harmonic
            double fvtxn_dphi3_docalib = phi - fvtx_north_psi3_docalib;
            double cosfvtxn_dphi3_docalib = TMath::Cos(3 * fvtxn_dphi3_docalib);

            fvtxn_v3eta_both_docalib[icent]->Fill(eta, cosfvtxn_dphi3_docalib);
            if ( west ) fvtxn_v3eta_west_docalib[icent]->Fill(eta, cosfvtxn_dphi3_docalib);
            if ( !west ) fvtxn_v3eta_east_docalib[icent]->Fill(eta, cosfvtxn_dphi3_docalib);

          } // check on ep north all layers

        } // check on fvtx clusters

      } // loop over fvtx tracks

    } // check on fvtx tracks

  }//end of event

  cout << "Processed " << event_counter << "/" << all_counter << " events (" << (float)event_counter / (float)all_counter << ")" << endl;
  cout << "Events with bad centrality = " << bad_cent_counter << " (" << (float)bad_cent_counter / (float)event_counter << ")" << endl;
  cout << "Events with vertex disagreement = " << bad_vertex_counter << " (" << (float)bad_vertex_counter / (float)event_counter << ")" << endl;
  cout << "Gap cuts applied " << gapcut_counter << "/" << cluster_counter << " (" << (float)gapcut_counter / (float)cluster_counter << ")" << endl;

  if ( rp_recal_pass < 3 && rp_recal_pass > 0 )
  {
    // --- previous pass calib file is named above, rename it here
    sprintf(calibfile, "output/flattening_data/flattening_%d_%d.dat", runNumber, rp_recal_pass);
    cout << "writing calibration file : " << calibfile << endl;
    ofstream ofs;
    ofs.open(calibfile);

    cout << "writing out flattening parameters" << endl;
    //int ic = 0;
    for ( int ic = 0; ic < NMUL; ic++ )
    {
      for ( int iz = 0; iz < NZPS; iz++ )
      {
        for ( int ih = 1; ih < NHAR; ih++ )
        {
          for ( int id = 0; id < NDETSHORT; id++ )
          {
            for ( int ib = 0; ib < 2; ib++ )
            {
              if ( DIAG ) cout << "writing ave:  " << ic << " " << iz << " " << ih << " " << id << endl;
              // write out average qx, qx error, qy, qy error
              ofs << ave[ic][iz][ih][id]->GetBinContent(ib + 1) << " ";
              ofs << ave[ic][iz][ih][id]->GetBinError  (ib + 1) << " ";
            } // x and y
            ofs << endl;
            for ( int ib = 0; ib < 2; ib++ )
            {
              for ( int io = 0; io < NORD; io++ )
              {
                // write first 12 orders for fourier fit of psi
                // we are unsure of what's being written out here
                ofs << flt[ic][iz][ih][id]->GetBinContent(ib * NORD + io + 1) << " ";
              } // orders
              ofs << endl;
            } // x and y
          } // detectors
        } // harmonics
      } // z_vertex bins
    } // centrality bins
    ofs.close();
  } // if pass 1 or 2



  cout << "Now attempting to close and write data files" << endl;

  cout << "histogram output file: " << outFile1 << endl;

  if ( rp_recal_pass == 3 )
  {
    mData1->Write();
    mData1->Close();
  }

  cout << "cleaning up" << endl;



  ntp_event_chain->Delete();

  // ---
  cout << "end of program ana" << endl;

  return;

}







void initialize_pmt_position()
{

  d_pmt_x[0] = -123;
  d_pmt_y[0] = 42.6;
  d_pmt_x[1] = -123;
  d_pmt_y[1] = 14.2;
  d_pmt_x[2] = -98.4;
  d_pmt_y[2] = 85.2;
  d_pmt_x[3] = -98.4;
  d_pmt_y[3] = 56.8;
  d_pmt_x[4] = -98.4;
  d_pmt_y[4] = 28.4;
  d_pmt_x[5] = -73.8;
  d_pmt_y[5] = 99.4;
  d_pmt_x[6] = -73.8;
  d_pmt_y[6] = 71;
  d_pmt_x[7] = -73.8;
  d_pmt_y[7] = 42.6;
  d_pmt_x[8] = -73.8;
  d_pmt_y[8] = 14.2;
  d_pmt_x[9] = -49.2;
  d_pmt_y[9] = 113.6;
  d_pmt_x[10] = -49.2;
  d_pmt_y[10] = 85.2;
  d_pmt_x[11] = -49.2;
  d_pmt_y[11] = 56.8;
  d_pmt_x[12] = -24.6;
  d_pmt_y[12] = 127.8;
  d_pmt_x[13] = -24.6;
  d_pmt_y[13] = 99.4;
  d_pmt_x[14] = -24.6;
  d_pmt_y[14] = 71;
  d_pmt_x[15] = 0;
  d_pmt_y[15] = 113.6;
  d_pmt_x[16] = 0;
  d_pmt_y[16] = 85.2;
  d_pmt_x[17] = 24.6;
  d_pmt_y[17] = 127.8;
  d_pmt_x[18] = 24.6;
  d_pmt_y[18] = 99.4;
  d_pmt_x[19] = 24.6;
  d_pmt_y[19] = 71;
  d_pmt_x[20] = 49.2;
  d_pmt_y[20] = 113.6;
  d_pmt_x[21] = 49.2;
  d_pmt_y[21] = 85.2;
  d_pmt_x[22] = 49.2;
  d_pmt_y[22] = 56.8;
  d_pmt_x[23] = 73.8;
  d_pmt_y[23] = 99.4;
  d_pmt_x[24] = 73.8;
  d_pmt_y[24] = 71;
  d_pmt_x[25] = 73.8;
  d_pmt_y[25] = 42.6;
  d_pmt_x[26] = 73.8;
  d_pmt_y[26] = 14.2;
  d_pmt_x[27] = 98.4;
  d_pmt_y[27] = 85.2;
  d_pmt_x[28] = 98.4;
  d_pmt_y[28] = 56.8;
  d_pmt_x[29] = 98.4;
  d_pmt_y[29] = 28.4;
  d_pmt_x[30] = 123;
  d_pmt_y[30] = 42.6;
  d_pmt_x[31] = 123;
  d_pmt_y[31] = 14.2;
  d_pmt_x[32] = 123;
  d_pmt_y[32] = -42.6;
  d_pmt_x[33] = 123;
  d_pmt_y[33] = -14.2;
  d_pmt_x[34] = 98.4;
  d_pmt_y[34] = -85.2;
  d_pmt_x[35] = 98.4;
  d_pmt_y[35] = -56.8;
  d_pmt_x[36] = 98.4;
  d_pmt_y[36] = -28.4;
  d_pmt_x[37] = 73.8;
  d_pmt_y[37] = -99.4;
  d_pmt_x[38] = 73.8;
  d_pmt_y[38] = -71;
  d_pmt_x[39] = 73.8;
  d_pmt_y[39] = -42.6;
  d_pmt_x[40] = 73.8;
  d_pmt_y[40] = -14.2;
  d_pmt_x[41] = 49.2;
  d_pmt_y[41] = -113.6;
  d_pmt_x[42] = 49.2;
  d_pmt_y[42] = -85.2;
  d_pmt_x[43] = 49.2;
  d_pmt_y[43] = -56.8;
  d_pmt_x[44] = 24.6;
  d_pmt_y[44] = -127.8;
  d_pmt_x[45] = 24.6;
  d_pmt_y[45] = -99.4;
  d_pmt_x[46] = 24.6;
  d_pmt_y[46] = -71;
  d_pmt_x[47] = -0;
  d_pmt_y[47] = -113.6;
  d_pmt_x[48] = -0;
  d_pmt_y[48] = -85.2;
  d_pmt_x[49] = -24.6;
  d_pmt_y[49] = -127.8;
  d_pmt_x[50] = -24.6;
  d_pmt_y[50] = -99.4;
  d_pmt_x[51] = -24.6;
  d_pmt_y[51] = -71;
  d_pmt_x[52] = -49.2;
  d_pmt_y[52] = -113.6;
  d_pmt_x[53] = -49.2;
  d_pmt_y[53] = -85.2;
  d_pmt_x[54] = -49.2;
  d_pmt_y[54] = -56.8;
  d_pmt_x[55] = -73.8;
  d_pmt_y[55] = -99.4;
  d_pmt_x[56] = -73.8;
  d_pmt_y[56] = -71;
  d_pmt_x[57] = -73.8;
  d_pmt_y[57] = -42.6;
  d_pmt_x[58] = -73.8;
  d_pmt_y[58] = -14.2;
  d_pmt_x[59] = -98.4;
  d_pmt_y[59] = -85.2;
  d_pmt_x[60] = -98.4;
  d_pmt_y[60] = -56.8;
  d_pmt_x[61] = -98.4;
  d_pmt_y[61] = -28.4;
  d_pmt_x[62] = -123;
  d_pmt_y[62] = -42.6;
  d_pmt_x[63] = -123;
  d_pmt_y[63] = -14.2;

}



int get_pmt_layer(int i)
{

  if ( i == 8 ||
       i == 11 ||
       i == 14 ||
       i == 19 ||
       i == 22 ||
       i == 26 ||
       i == 40 ||
       i == 43 ||
       i == 46 ||
       i == 51 ||
       i == 54 ||
       i == 58 ) return 0; // inner layer

  if ( i == 7 ||
       i == 16 ||
       i == 25 ||
       i == 39 ||
       i == 48 ||
       i == 57 ) return 1; // inner middle layer

  if ( i == 4 ||
       i == 6 ||
       i == 10 ||
       i == 13 ||
       i == 18 ||
       i == 21 ||
       i == 24 ||
       i == 29 ||
       i == 36 ||
       i == 38 ||
       i == 42 ||
       i == 45 ||
       i == 45 ||
       i == 50 ||
       i == 53 ||
       i == 56 ||
       i == 61 ) return 2; // middle layer

  if ( i == 3 ||
       i == 15 ||
       i == 28 ||
       i == 35 ||
       i == 47 ||
       i == 60 ) return 3; //outer middle layer

  if ( i == 0 ||
       i == 1 ||
       i == 2 ||
       i == 5 ||
       i == 9 ||
       i == 12 ||
       i == 17 ||
       i == 20 ||
       i == 23 ||
       i == 27 ||
       i == 30 ||
       i == 31 ||
       i == 32 ||
       i == 33 ||
       i == 34 ||
       i == 37 ||
       i == 41 ||
       i == 44 ||
       i == 49 ||
       i == 52 ||
       i == 55 ||
       i == 59 ||
       i == 62 ||
       i == 63 ) return 4; // outer layer

  return -1;

}


int get_fvtx_layer(float z)
{
  // --- south side
  if ( z < -18 && z > -24 ) return 0;
  if ( z < -24 && z > -30 ) return 1;
  if ( z < -30 && z > -35 ) return 2;
  if ( z < -35 )            return 3;
  // --- north side
  if ( z > 18 && z < 24 ) return 0;
  if ( z > 24 && z < 30 ) return 1;
  if ( z > 30 && z < 35 ) return 2;
  if ( z > 35 )           return 3;
  // --- invalid numbers...
  cout << "get_fvtx_layer::invalid z =  " << z << endl;
  return -1;
}


float calc2_event(float Xn, float Yn, float M)
{

  float numerator = Xn * Xn + Yn * Yn - M;
  float denominator = M * (M - 1);

  return numerator / denominator;

}


float calc4_event(float Xn, float Yn, float X2n, float Y2n, float M)
{

  if ( M < 5 ) return -9999;

  float Qn2 = Xn * Xn + Yn * Yn;
  float Qn2d = Xn * Xn - Yn * Yn;

  float one   = Qn2 * Qn2;
  float two   = X2n * X2n + Y2n * Y2n;
  float three = (2 * (X2n * Qn2d + 2 * Y2n * Xn * Yn));
  float four  = 2 * (2 * (M - 2) * Qn2);
  float five  = 2 * M * (M - 3);

  float numerator = one + two - three - four + five;
  float denominator = M * (M - 1) * (M - 2) * (M - 3);

  return numerator / denominator;

}


float calc4_track_flag(float xn, float yn, float x2n, float y2n, float Xn, float Yn, float X2n, float Y2n, float M, bool is_POI_in_RP)
{

  if ( is_POI_in_RP )
  {
    float one   = (xn * Xn + yn * Yn) * (Xn * Xn + Yn * Yn);
    float two   = x2n * Xn * Xn - x2n * Yn * Yn + 2 * y2n * Xn * Yn;
    float three = xn * Xn * X2n + xn * Yn * Y2n - yn * (X2n * Yn - Xn * Y2n);
    float four  = 2 * M * (xn * Xn + yn * Yn);
    float five  = 2 * (Xn * Xn + Yn * Yn);
    float six   = 7 * (xn * Xn + yn * Yn);
    float seven = xn * Xn + yn * Yn;
    float eight = x2n * X2n + y2n * Y2n;
    float nine = 2 * (xn * Xn + yn * Yn);
    // ---
    float numerator = one - two - three - four - five + six - seven + eight + nine + 2 * M - 6;
    float denominator = (M - 1) * (M - 2) * (M - 3);
    // ---
    return numerator / denominator;
  }
  else
  {
    float one   = (xn * Xn + yn * Yn) * (Xn * Xn + Yn * Yn);
    float three = xn * Xn * X2n + xn * Yn * Y2n - yn * (X2n * Yn - Xn * Y2n);
    float four  = 2 * M * (xn * Xn + yn * Yn);
    float nine = 2 * (xn * Xn + yn * Yn);
    // ---
    float numerator = one - three - four + nine;
    float denominator = M * (M - 1) * (M - 2);
    // ---
    return numerator / denominator;
  }

}
