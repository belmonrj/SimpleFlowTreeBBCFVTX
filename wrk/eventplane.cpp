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

  flatten(run, 1);
  flatten(run, 2);
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

  float fracCut = 0.95; // pile up rejection < fracCut (better place??)
  // float fracCut = 0.98; // pile up rejection < fracCut (better place??)

  bool tight_trkcuts = false; // flag for tight cnt & fvtx track cuts (true=tight)

  // float reso_lpt = 0.4; // low-pT cut for PHCentralTracks when calculating ep resolution
  float reso_lpt = 0.2; // low-pT cut for PHCentralTracks when calculating ep resolution
  // float reso_hpt = 3.0; // hi-pT cut for PHCentralTracks when calculating ep resolution
  float reso_hpt = 2.0; // hi-pT cut for PHCentralTracks when calculating ep resolution

  //---
  // all zeros
  float qy2_offset_fvtxs[4][6] = {{0}, {0}, {0}, {0}};
  float qy2_offset_bbcs[4][6] = {{0}, {0}, {0}, {0}};
  float qx3_offset_fvtxs[4][6] = {{0}, {0}, {0}, {0}};
  float qx3_offset_bbcs[4][6] = {{0}, {0}, {0}, {0}};
  float qy2_offset_fvtxsl[4][4][6] ={{{0}}};
  float qx3_offset_fvtxsl[4][4][6] ={{{0}}};
  //---

  //---
  // full values (4th iteration)
//   float qy2_offset_fvtxs[4][6] =
//   {
//     { -0.0026, -0.0000, -0.0000, -0.0130, -0.0407, -0.1108}, //200 GeV
//     { -0.0017, -0.0072, -0.0126, -0.0244, -0.0566, -0.1332}, //62 GeV
//     { -0.0028, -0.0045, -0.0164, -0.0333, -0.0433, 0.0}, //39 GeV
//     { -0.0154, -0.0139, -0.0110, -0.0497, 0.0, 0.0}, //20 GeV
//   };

//   // itr 4
//   float qy2_offset_fvtxsl[4][4][6] =
//   {
// //200 GeV
//     {
//       { -0.0244, -0.0338, -0.0342, -0.0595, -0.1147, -0.2237}, // layer 0
//       { +0.0037, +0.0114, +0.0145, +0.0183, +0.0186, +0.0300}, // layer 1
//       { +0.0020, +0.0143, +0.0113, -0.0036, -0.0067, -0.0319}, // layer 2
//       { +0.0063, +0.0233, +0.0201, +0.0143, +0.0156, +0.0436}, // layer 3
//     },
// //62 GeV
//     {
//       { -0.0325, -0.0450, -0.0487, -0.0865, -0.1272, -0.2329}, // layer 0
//       { +0.0137, +0.0100, +0.0101, +0.0111, +0.0207, +0.0290}, // layer 1
//       { +0.0027, -0.0023, -0.0060, -0.0091, -0.0256, -0.0623}, // layer 2
//       { +0.0084, +0.0093, +0.0041, +0.0069, +0.0158, +0.0318}, // layer 3
//     },
// //39 GeV
//     {
//       { -0.0286, -0.0343, -0.0504, -0.0704, -0.0910, +0.0000}, // layer 0
//       { +0.0168, +0.0155, +0.0151, +0.0165, +0.0240, +0.0000}, // layer 1
//       { +0.0080, +0.0029, +0.0047, -0.0032, -0.0026, +0.0000}, // layer 2
//       { -0.0002, -0.0023, -0.0066, -0.0060, -0.0163, +0.0000}, // layer 3
//     },
// //20 GeV
//     {
//       { -0.0577, -0.0432, -0.0828, -0.0411, +0.0000, +0.0000}, // layer 0
//       { +0.0062, +0.0129, +0.0297, -0.0732, +0.0000, +0.0000}, // layer 1
//       { -0.0064, -0.0006, -0.0110, -0.0990, +0.0000, +0.0000}, // layer 2
//       { -0.0058, +0.0101, +0.0170, -0.0472, +0.0000, +0.0000}, // layer 3
//     },
//   };


  //---
  //---
  // Systematic check (x1.50)
  //---
  // float qy2_offset_fvtxs[4][6] =
  // {
  //   { -0.000, -0.000, -0.000, -0.024, -0.066, -0.181}, //200 GeV

  //   { -0.002, -0.012, -0.021, -0.036, -0.082, -0.194}, //62 GeV

  //   { -0.000, -0.008, -0.021, -0.042, -0.076, -0.154}, //39 GeV

  //   { -0.002, -0.012, -0.021, -0.042, -0.076, -0.181}, //20 GeV
  // };
  //---
  // Systematic check (x0.50)
  // float qy2_offset_fvtxs[4][6] =
  // {
  //   { -0.000, -0.000, -0.000, -0.007, -0.018, -0.048}, //200 GeV

  //   { -0.000, -0.003, -0.006, -0.009, -0.022, -0.051}, //62 GeV

  //   { -0.000, -0.002, -0.006, -0.011, -0.020, -0.042}, //39 GeV

  //   { -0.000, -0.003, -0.006, -0.011, -0.020, -0.048}, //20 GeV
  // };
  //---


  //---
  // 2nd order Qy offset for BBCS
  // float qy2_offset_bbcs[4][6] =
  // {
  //   { -0.0037, -0.0096, -0.0000, -0.0000, -0.0082, -0.0000}, //200 GeV
  //   { -0.0020, -0.0020, -0.0036, -0.0041, +0.0010, -0.0009}, //62 GeV
  //   { -0.0052, -0.0038, +0.0014, -0.0019, -0.0060, -0.0000}, //39 GeV
  //   { +0.0005, -0.0181, +0.0002, -0.0451, -0.0000, -0.0000}, //20 GeV
  // };
  //---


  //--
  // 3rd order Qx offset for FVTXS
// itr 2
  // float qx3_offset_fvtxs[4][6] =
  // {
  //   { -0.0004, -0.0069, +0.0038, +0.0202, +0.0543, +0.1627 }, // 200 GeV
  //   { +0.0067, +0.0105, +0.0129, +0.0353, +0.0666, +0.1466 }, // 62 GeV
  //   { +0.0085, +0.0086, +0.0229, +0.0310, +0.0392, +0.0000 }, // 39 GeV
  //   { -0.0063, +0.0099, +0.0223, -0.0174, +0.0000, +0.0000 }, // 20 GeV
  // };
  //--

  //--
  // 3rd order Qx offset for BBCS
// itr 2
  // float qx3_offset_bbcs[4][6] =
  // {
  //   { +0.0059, -0.0040, +0.0083, +0.0103, -0.0105, -0.0012 }, // 200 GeV
  //   { +0.0021, +0.0001, -0.0019, +0.0037, -0.0024, -0.0024 }, // 62 GeV
  //   { +0.0011, +0.0007, +0.0021, +0.0017, +0.0019, +0.0000 }, // 39 GeV
  //   { +0.0189, +0.0119, -0.0079, -0.0904, +0.0000, +0.0000 }, // 20 GeV
  // };
  //--

  //-- 3rd order Qx for FVTXSL
// itr 2
//   float qx3_offset_fvtxsl[4][4][6] =
//   {
// // 200 GeV
//     {
//       { +0.0276, +0.0355, +0.0386, +0.0717, +0.1333, +0.2943}, // layer 0
//       { -0.0054, -0.0147, -0.0159, -0.0105, -0.0019, +0.0142}, // layer 1
//       { -0.0062, -0.0086, -0.0043, -0.0021, +0.0048, +0.0231}, // layer 2
//       { -0.0094, -0.0114, -0.0128, -0.0036, -0.0001, +0.0259}, // layer 3
//     },
// // 62 GeV
//     {
//       { +0.0381, +0.0425, +0.0553, +0.0887, +0.1464, +0.2891}, // layer 0
//       { -0.0078, -0.0077, -0.0108, +0.0012, +0.0020, +0.0266}, // layer 1
//       { -0.0002, -0.0006, +0.0005, +0.0042, +0.0073, +0.0228}, // layer 2
//       { -0.0059, -0.0034, -0.0036, +0.0012, -0.0019, +0.0017}, // layer 3
//     },
// // 39 GeV
//     {
//       { +0.0366, +0.0398, +0.0621, +0.0823, +0.1060, +0.0000}, // layer 0
//       { -0.0032, -0.0060, -0.0043, -0.0001, -0.0053, +0.0000}, // layer 1
//       { +0.0011, +0.0013, +0.0066, +0.0067, +0.0055, +0.0000}, // layer 2
//       { -0.0054, -0.0071, -0.0027, +0.0003, -0.0019, +0.0000}, // layer 3
//     },
// // 20 GeV
//     {
//       { +0.0246, +0.0542, +0.0758, +0.0816, +0.0000, +0.0000}, // layer 0
//       { -0.0165, -0.0097, -0.0072, -0.0960, +0.0000, +0.0000}, // layer 1
//       { +0.0027, -0.0012, +0.0029, -0.0194, +0.0000, +0.0000}, // layer 2
//       { -0.0189, -0.0062, -0.0154, -0.0036, +0.0000, +0.0000}, // layer 3
//     },
//   };
  
  //--
  cout << " frac cut: " << fracCut << endl;
  cout << " tight_trkcuts: " << tight_trkcuts << endl;
  cout << " Qy2 offset fvtxs: " << endl;
  for (int i = 0; i < NMUL; i++)
    cout << "    " << i << " " << qy2_offset_fvtxs[eidx][i] << endl;
  cout << " Qy2 offset bbcs: " << endl;
  for (int i = 0; i < NMUL; i++)
    cout << "    " << i << " " << qy2_offset_bbcs[eidx][i] << endl;
  cout << " Qx3 offset fvtxs: " << endl;
  for (int i = 0; i < NMUL; i++)
    cout << "    " << i << " " << qx3_offset_fvtxs[eidx][i] << endl;
  cout << " Qx3 offset bbcs: " << endl;
  for (int i = 0; i < NMUL; i++)
    cout << "    " << i << " " << qx3_offset_bbcs[eidx][i] << endl;

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
  int fvtxsl_index    =  3;
  int fvtxsa_index    =  7;

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

  TH1D* th1d_cent_MBonly = new TH1D("th1d_cent_MBonly", "", 100, -0.5, 99.5);
  TH1D* th1d_cent_CNonly = new TH1D("th1d_cent_CNonly", "", 100, -0.5, 99.5);
  TH1D* th1d_cent_CNetMB = new TH1D("th1d_cent_CNetMB", "", 100, -0.5, 99.5);


  TProfile* tp1f_bbc_fcharge_ring = new TProfile("tp1f_bbc_fcharge_ring", ";BBC ring; charge/total", 5, -0.5, 4.5);

  TH2D* th2d_BBCq_cent = new TH2D("th2d_bbcq_cent", "", 200, -0.5, 199.5, 101, -0.5, 100.5);
  TH2D* th2d_FVTXSntrk_cent = new TH2D("th2d_fvtxsntrk_cent", "", 100, -0.5, 99.5, 101, -0.5, 100.5);
  TH2D* th2d_FVTXnntrk_cent = new TH2D("th2d_fvtxnntrk_cent", "", 100, -0.5, 99.5, 101, -0.5, 100.5);
  TH2D* th2d_FVTXSnclus_cent = new TH2D("th2d_fvtxsnclus_cent", "", 200, -0.5, 1999.5, 101, -0.5, 100.5);
  TH2D* th2d_FVTXnnclus_cent = new TH2D("th2d_fvtxnnclus_cent", "", 200, -0.5, 1999.5, 101, -0.5, 100.5);

  TH3D* th3d_FVTXSclus_eta_zvrtx_lay =
    new TH3D("th3d_FVTXSclus_eta_zvrtx_lay", "",
             40, -4, 0,
             110, -11, 11,
             5, -0.5, 4.5);

  TH3D* th3d_FVTXSclus_eta_zvrtx_cent =
    new TH3D("th3d_FVTXSclus_eta_zvrtx_cent", "",
             40, -4, 0,
             110, -11, 11,
             101, -0.5, 100.5);

  TH1D* th1d_nevent_MB = new TH1D("th1d_nevent_mb", "", 101, -0.5, 100.5);
  TH1D* th1d_nevent_HM = new TH1D("th1d_nevent_hm", "", 101, -0.5, 100.5);
  TH1D* th1d_nevent_all = new TH1D("th1d_nevent_all", "", 101, -0.5, 100.5);

  // --- event plane resolution, need to be centrality dependent
  TProfile* tp1f_reso_BBC_CNT[NHAR][NMUL][NZPS];
  TProfile* tp1f_reso_BBC_FVTX[NHAR][NMUL][NZPS];
  TProfile* tp1f_reso_CNT_FVTX[NHAR][NMUL][NZPS];
  TProfile* tp1f_reso_BBC_FVTXN[NHAR][NMUL][NZPS];
  TProfile* tp1f_reso_CNT_FVTXN[NHAR][NMUL][NZPS];
  TProfile* tp1f_reso_FVTXS_FVTXN[NHAR][NMUL][NZPS];
  TProfile* tp1f_reso_BBC_FVTXSA[NHAR][NMUL][NZPS];
  TProfile* tp1f_reso_CNT_FVTXSA[NHAR][NMUL][NZPS];
  TProfile* tp1f_reso_BBC_FVTXSB[NHAR][NMUL][NZPS];
  TProfile* tp1f_reso_CNT_FVTXSB[NHAR][NMUL][NZPS];
  TProfile* tp1f_reso_BBC_FVTXSL[NHAR][NMUL][NZPS][NFVTXLAY];
  TProfile* tp1f_reso_CNT_FVTXSL[NHAR][NMUL][NZPS][NFVTXLAY];

  for (int in = 1; in < NHAR; in++)
  {
    for (int ic = 0; ic < NMUL; ic++)
    {
      for (int iz = 0; iz < NZPS; iz++)
      {
        tp1f_reso_BBC_CNT[in][ic][iz]  = new TProfile(Form("tp1f_c%i_z%i_reso%i_BBC_CNT", ic, iz, in + 1), "", 1, -0.5, 0.5, -1e6, 1e6, "");
        tp1f_reso_BBC_FVTX[in][ic][iz] = new TProfile(Form("tp1f_c%i_z%i_reso%i_BBC_FVTX", ic, iz, in + 1), "", 1, -0.5, 0.5, -1e6, 1e6, "");
        tp1f_reso_CNT_FVTX[in][ic][iz] = new TProfile(Form("tp1f_c%i_z%i_reso%i_CNT_FVTX", ic, iz, in + 1), "", 1, -0.5, 0.5, -1e6, 1e6, "");
        tp1f_reso_BBC_FVTXN[in][ic][iz] = new TProfile(Form("tp1f_c%i_z%i_reso%i_BBC_FVTXN", ic, iz, in + 1), "", 1, -0.5, 0.5, -1e6, 1e6, "");
        tp1f_reso_CNT_FVTXN[in][ic][iz] = new TProfile(Form("tp1f_c%i_z%i_reso%i_CNT_FVTXN", ic, iz, in + 1), "", 1, -0.5, 0.5, -1e6, 1e6, "");
        tp1f_reso_FVTXS_FVTXN[in][ic][iz] = new TProfile(Form("tp1f_c%i_z%i_reso%i_FVTXS_FVTXN", ic, iz, in + 1), "", 1, -0.5, 0.5, -1e6, 1e6, "");

        tp1f_reso_BBC_FVTXSA[in][ic][iz] = new TProfile(Form("tp1f_c%i_z%i_reso%i_BBC_FVTXSA", ic, iz, in + 1), "", 1, -0.5, 0.5, -1e6, 1e6, "");
        tp1f_reso_CNT_FVTXSA[in][ic][iz] = new TProfile(Form("tp1f_c%i_z%i_reso%i_CNT_FVTXSA", ic, iz, in + 1), "", 1, -0.5, 0.5, -1e6, 1e6, "");
        tp1f_reso_BBC_FVTXSB[in][ic][iz] = new TProfile(Form("tp1f_c%i_z%i_reso%i_BBC_FVTXSB", ic, iz, in + 1), "", 1, -0.5, 0.5, -1e6, 1e6, "");
        tp1f_reso_CNT_FVTXSB[in][ic][iz] = new TProfile(Form("tp1f_c%i_z%i_reso%i_CNT_FVTXSB", ic, iz, in + 1), "", 1, -0.5, 0.5, -1e6, 1e6, "");
        for (int il = 0; il < NFVTXLAY; il++)
        {
          tp1f_reso_BBC_FVTXSL[in][ic][iz][il] = new TProfile(Form("tp1f_c%i_z%i_reso%i_BBC_FVTXSL%i", ic, iz, in + 1, il), "", 1, -0.5, 0.5, -1e6, 1e6, "");
          tp1f_reso_CNT_FVTXSL[in][ic][iz][il] = new TProfile(Form("tp1f_c%i_z%i_reso%i_CNT_FVTXSL%i", ic, iz, in + 1, il), "", 1, -0.5, 0.5, -1e6, 1e6, "");
        }
      } // iz
    } // ic
  } // in



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
  TProfile* cnt_west_cosnphi[NHAR][NMUL][NZPS];
  TProfile* cnt_east_cosnphi[NHAR][NMUL][NZPS];
  TProfile* cnt_both_cosnphi[NHAR][NMUL][NZPS];

  TProfile* cnt_west_sinnphi[NHAR][NMUL][NZPS];
  TProfile* cnt_east_sinnphi[NHAR][NMUL][NZPS];
  TProfile* cnt_both_sinnphi[NHAR][NMUL][NZPS];

  TProfile* bbcs_vn_east_docalib[NHAR][NMUL][NZPS];
  TProfile* bbcs_vn_west_docalib[NHAR][NMUL][NZPS];
  TProfile* bbcs_vn_both_docalib[NHAR][NMUL][NZPS];

  TProfile* fvtxs_vn_east_docalib[NHAR][NMUL][NZPS];
  TProfile* fvtxs_vn_west_docalib[NHAR][NMUL][NZPS];
  TProfile* fvtxs_vn_both_docalib[NHAR][NMUL][NZPS];

  TProfile* fvtxsa_vn_east_docalib[NHAR][NMUL][NZPS];
  TProfile* fvtxsa_vn_west_docalib[NHAR][NMUL][NZPS];
  TProfile* fvtxsa_vn_both_docalib[NHAR][NMUL][NZPS];

  TProfile* fvtxsb_vn_east_docalib[NHAR][NMUL][NZPS];
  TProfile* fvtxsb_vn_west_docalib[NHAR][NMUL][NZPS];
  TProfile* fvtxsb_vn_both_docalib[NHAR][NMUL][NZPS];

  TProfile* fvtxsl_vn_east_docalib[NHAR][NMUL][NZPS][NFVTXLAY];
  TProfile* fvtxsl_vn_west_docalib[NHAR][NMUL][NZPS][NFVTXLAY];
  TProfile* fvtxsl_vn_both_docalib[NHAR][NMUL][NZPS][NFVTXLAY];

  TProfile* fvtxn_vn_east_docalib[NHAR][NMUL][NZPS];
  TProfile* fvtxn_vn_west_docalib[NHAR][NMUL][NZPS];
  TProfile* fvtxn_vn_both_docalib[NHAR][NMUL][NZPS];


  for ( int in = 1; in < NHAR; ++in )
  {
    for ( int ic = 0; ic < NMUL; ++ic )
    {
      for ( int iz = 0; iz < NZPS; ++iz )
      {
        int n = in + 1; // harmonic number

        //-- cos(nphi)
        cnt_west_cosnphi[in][ic][iz] = new TProfile(Form("cnt_west_cos%iphi_cent%d_zvtx%d", n, ic, iz), "", NPTBINS, ptlim, -1.1, 1.1);
        cnt_east_cosnphi[in][ic][iz] = new TProfile(Form("cnt_east_cos%iphi_cent%d_zvtx%d", n, ic, iz), "", NPTBINS, ptlim, -1.1, 1.1);
        cnt_both_cosnphi[in][ic][iz] = new TProfile(Form("cnt_both_cos%iphi_cent%d_zvtx%d", n, ic, iz), "", NPTBINS, ptlim, -1.1, 1.1);

        //-- sin(nphi)
        cnt_west_sinnphi[in][ic][iz] = new TProfile(Form("cnt_west_sin%iphi_cent%d_zvtx%d", n, ic, iz), "", NPTBINS, ptlim, -1.1, 1.1);
        cnt_east_sinnphi[in][ic][iz] = new TProfile(Form("cnt_east_sin%iphi_cent%d_zvtx%d", n, ic, iz), "", NPTBINS, ptlim, -1.1, 1.1);
        cnt_both_sinnphi[in][ic][iz] = new TProfile(Form("cnt_both_sin%iphi_cent%d_zvtx%d", n, ic, iz), "", NPTBINS, ptlim, -1.1, 1.1);

        //-- vn
        bbcs_vn_both_docalib[in][ic][iz] = new TProfile(Form("bbcs_v%i_both_docalib_cent%d_zvtx%d", n, ic, iz), "", NPTBINS, ptlim, -1.1, 1.1);
        bbcs_vn_east_docalib[in][ic][iz] = new TProfile(Form("bbcs_v%i_east_docalib_cent%d_zvtx%d", n, ic, iz), "", NPTBINS, ptlim, -1.1, 1.1);
        bbcs_vn_west_docalib[in][ic][iz] = new TProfile(Form("bbcs_v%i_west_docalib_cent%d_zvtx%d", n, ic, iz), "", NPTBINS, ptlim, -1.1, 1.1);

        fvtxs_vn_both_docalib[in][ic][iz] = new TProfile(Form("fvtxs_v%i_both_docalib_cent%d_zvtx%d", n, ic, iz), "", NPTBINS, ptlim, -1.1, 1.1);
        fvtxs_vn_east_docalib[in][ic][iz] = new TProfile(Form("fvtxs_v%i_east_docalib_cent%d_zvtx%d", n, ic, iz), "", NPTBINS, ptlim, -1.1, 1.1);
        fvtxs_vn_west_docalib[in][ic][iz] = new TProfile(Form("fvtxs_v%i_west_docalib_cent%d_zvtx%d", n, ic, iz), "", NPTBINS, ptlim, -1.1, 1.1);

        fvtxsa_vn_both_docalib[in][ic][iz] = new TProfile(Form("fvtxsa_v%i_both_docalib_cent%d_zvtx%d", n, ic, iz), "", NPTBINS, ptlim, -1.1, 1.1);
        fvtxsa_vn_east_docalib[in][ic][iz] = new TProfile(Form("fvtxsa_v%i_east_docalib_cent%d_zvtx%d", n, ic, iz), "", NPTBINS, ptlim, -1.1, 1.1);
        fvtxsa_vn_west_docalib[in][ic][iz] = new TProfile(Form("fvtxsa_v%i_west_docalib_cent%d_zvtx%d", n, ic, iz), "", NPTBINS, ptlim, -1.1, 1.1);

        fvtxsb_vn_both_docalib[in][ic][iz] = new TProfile(Form("fvtxsb_v%i_both_docalib_cent%d_zvtx%d", n, ic, iz), "", NPTBINS, ptlim, -1.1, 1.1);
        fvtxsb_vn_east_docalib[in][ic][iz] = new TProfile(Form("fvtxsb_v%i_east_docalib_cent%d_zvtx%d", n, ic, iz), "", NPTBINS, ptlim, -1.1, 1.1);
        fvtxsb_vn_west_docalib[in][ic][iz] = new TProfile(Form("fvtxsb_v%i_west_docalib_cent%d_zvtx%d", n, ic, iz), "", NPTBINS, ptlim, -1.1, 1.1);

        for (int il = 0; il < NFVTXLAY; il++)
        {
          fvtxsl_vn_both_docalib[in][ic][iz][il] = new TProfile(Form("fvtxsl%i_v%i_both_docalib_cent%d_zvtx%d", il, n, ic, iz), "", NPTBINS, ptlim, -1.1, 1.1);
          fvtxsl_vn_east_docalib[in][ic][iz][il] = new TProfile(Form("fvtxsl%i_v%i_east_docalib_cent%d_zvtx%d", il, n, ic, iz), "", NPTBINS, ptlim, -1.1, 1.1);
          fvtxsl_vn_west_docalib[in][ic][iz][il] = new TProfile(Form("fvtxsl%i_v%i_west_docalib_cent%d_zvtx%d", il, n, ic, iz), "", NPTBINS, ptlim, -1.1, 1.1);
        }

        fvtxn_vn_both_docalib[in][ic][iz] = new TProfile(Form("fvtxn_v%i_both_docalib_cent%d_zvtx%d", n, ic, iz), "", NPTBINS, ptlim, -1.1, 1.1);
        fvtxn_vn_east_docalib[in][ic][iz] = new TProfile(Form("fvtxn_v%i_east_docalib_cent%d_zvtx%d", n, ic, iz), "", NPTBINS, ptlim, -1.1, 1.1);
        fvtxn_vn_west_docalib[in][ic][iz] = new TProfile(Form("fvtxn_v%i_west_docalib_cent%d_zvtx%d", n, ic, iz), "", NPTBINS, ptlim, -1.1, 1.1);

      } // iz
    } // ic
  } // in


  // vs eta
  TProfile* cosnphi_eta_east[NHAR][NMUL][NZPS];
  TProfile* cosnphi_eta_west[NHAR][NMUL][NZPS];
  TProfile* cosnphi_eta_both[NHAR][NMUL][NZPS];

  TProfile* sinnphi_eta_east[NHAR][NMUL][NZPS];
  TProfile* sinnphi_eta_west[NHAR][NMUL][NZPS];
  TProfile* sinnphi_eta_both[NHAR][NMUL][NZPS];

  TProfile* bbcs_vneta_east_docalib[NHAR][NMUL][NZPS];
  TProfile* bbcs_vneta_west_docalib[NHAR][NMUL][NZPS];
  TProfile* bbcs_vneta_both_docalib[NHAR][NMUL][NZPS];

  TProfile* fvtxs_vneta_east_docalib[NHAR][NMUL][NZPS];
  TProfile* fvtxs_vneta_west_docalib[NHAR][NMUL][NZPS];
  TProfile* fvtxs_vneta_both_docalib[NHAR][NMUL][NZPS];

  TProfile* fvtxn_vneta_east_docalib[NHAR][NMUL][NZPS];
  TProfile* fvtxn_vneta_west_docalib[NHAR][NMUL][NZPS];
  TProfile* fvtxn_vneta_both_docalib[NHAR][NMUL][NZPS];

  for ( int in = 1; in < NHAR; ++in )
  {
    for ( int ic = 0; ic < NMUL; ++ic )
    {
      for ( int iz = 0; iz < NZPS; ++iz )
      {
        int n = in + 1; // harmonic number

        // <cos(2*phi)> or <sin(2*phi)>
        cosnphi_eta_east[in][ic][iz] = new TProfile(Form("cos%iphi_eta_east_cent%d_zvtx%d", n, ic, iz), "", 32, -3.2, 3.2, -1.1, 1.1);
        cosnphi_eta_west[in][ic][iz] = new TProfile(Form("cos%iphi_eta_west_cent%d_zvtx%d", n, ic, iz), "", 32, -3.2, 3.2, -1.1, 1.1);
        cosnphi_eta_both[in][ic][iz] = new TProfile(Form("cos%iphi_eta_both_cent%d_zvtx%d", n, ic, iz), "", 32, -3.2, 3.2, -1.1, 1.1);

        sinnphi_eta_east[in][ic][iz] = new TProfile(Form("sin%iphi_eta_east_cent%d_zvtx%d", n, ic, iz), "", 32, -3.2, 3.2, -1.1, 1.1);
        sinnphi_eta_west[in][ic][iz] = new TProfile(Form("sin%iphi_eta_west_cent%d_zvtx%d", n, ic, iz), "", 32, -3.2, 3.2, -1.1, 1.1);
        sinnphi_eta_both[in][ic][iz] = new TProfile(Form("sin%iphi_eta_both_cent%d_zvtx%d", n, ic, iz), "", 32, -3.2, 3.2, -1.1, 1.1);

        // v2 vs eta
        bbcs_vneta_both_docalib[in][ic][iz] = new TProfile(Form("bbcs_v%ieta_both_docalib_cent%d_zvtx%d", n, ic, iz), "", 32, -3.2, 3.2, -1.1, 1.1);
        bbcs_vneta_east_docalib[in][ic][iz] = new TProfile(Form("bbcs_v%ieta_east_docalib_cent%d_zvtx%d", n, ic, iz), "", 32, -3.2, 3.2, -1.1, 1.1);
        bbcs_vneta_west_docalib[in][ic][iz] = new TProfile(Form("bbcs_v%ieta_west_docalib_cent%d_zvtx%d", n, ic, iz), "", 32, -3.2, 3.2, -1.1, 1.1);

        fvtxs_vneta_both_docalib[in][ic][iz] = new TProfile(Form("fvtxs_v%ieta_both_docalib_cent%d_zvtx%d", n, ic, iz), "", 32, -3.2, 3.2, -1.1, 1.1);
        fvtxs_vneta_east_docalib[in][ic][iz] = new TProfile(Form("fvtxs_v%ieta_east_docalib_cent%d_zvtx%d", n, ic, iz), "", 32, -3.2, 3.2, -1.1, 1.1);
        fvtxs_vneta_west_docalib[in][ic][iz] = new TProfile(Form("fvtxs_v%ieta_west_docalib_cent%d_zvtx%d", n, ic, iz), "", 32, -3.2, 3.2, -1.1, 1.1);

        fvtxn_vneta_both_docalib[in][ic][iz] = new TProfile(Form("fvtxn_v%ieta_both_docalib_cent%d_zvtx%d", n, ic, iz), "", 32, -3.2, 3.2, -1.1, 1.1);
        fvtxn_vneta_east_docalib[in][ic][iz] = new TProfile(Form("fvtxn_v%ieta_east_docalib_cent%d_zvtx%d", n, ic, iz), "", 32, -3.2, 3.2, -1.1, 1.1);
        fvtxn_vneta_west_docalib[in][ic][iz] = new TProfile(Form("fvtxn_v%ieta_west_docalib_cent%d_zvtx%d", n, ic, iz), "", 32, -3.2, 3.2, -1.1, 1.1);

      } // iz
    } // ic
  } // in

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

    bool passes_mb = false;
    if ( energyflag == 200 ) passes_mb = trigger_scaled & trigger_BBCLL1narrow;
    if ( energyflag == 62  ) passes_mb = trigger_scaled & trigger_BBCLL1narrow;
    if ( energyflag == 20  ) passes_mb = trigger_scaled & trigger_FVTXNSBBCS;
    if ( energyflag == 39  ) passes_mb = trigger_scaled & trigger_FVTXNSBBCS;

    bool passes_hm = false;
    if ( energyflag == 200 ) passes_hm = trigger_scaled & trigger_BBCLL1narrowcent;
    if ( energyflag == 62  ) passes_hm = trigger_scaled & trigger_BBCLL1narrowcent;
    if ( energyflag == 20  ) passes_hm = trigger_scaled & trigger_FVTXNSBBCScentral;
    if ( energyflag == 39  ) passes_hm = trigger_scaled & trigger_FVTXNSBBCScentral;

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

    // --- event counting
    if ( passes_mb )
      th1d_nevent_MB->Fill(centrality);
    if ( passes_hm )
      th1d_nevent_HM->Fill(centrality);
    if ( passes_trigger )
      th1d_nevent_all->Fill(centrality);



    //if ( ( say_event && verbosity > 0 ) || verbosity > 4 ) cout << "Calculating event planes" << endl;

    // --- all numbers from Darren 2016-06-23
    const float x_off = 0.3;
    const float y_off = 0.02;
    const float beam_angle = 0.001;
    // const float beam_angle = 0.00;
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

    float bbc_qxn[NHAR] = {0};
    float bbc_qyn[NHAR] = {0};
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

        for (int ih = 1; ih < NHAR; ih++)
        {
          bbc_qxn[ih] += bbc_charge_corrected * TMath::Cos((ih + 1) * phi);
          bbc_qyn[ih] += bbc_charge_corrected * TMath::Sin((ih + 1) * phi);
        }

        // bbc_qx2 += bbc_charge_corrected * TMath::Cos(2 * phi);
        // bbc_qy2 += bbc_charge_corrected * TMath::Sin(2 * phi);
        // bbc_qx3 += bbc_charge_corrected * TMath::Cos(3 * phi);
        // bbc_qy3 += bbc_charge_corrected * TMath::Sin(3 * phi);
        // bbc_qx4 += bbc_charge_corrected * TMath::Cos(4 * phi);
        // bbc_qy4 += bbc_charge_corrected * TMath::Sin(4 * phi);
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

    // if ( isCN || isMB  ) th1d_cent_CNetMB->Fill(centrality);
    // if ( isCN && !isMB ) th1d_cent_CNonly->Fill(centrality);
    // if ( isMB && !isCN ) th1d_cent_MBonly->Fill(centrality);
    th1d_cent_CNetMB->Fill(centrality);
    if ( passes_hm ) th1d_cent_CNonly->Fill(centrality);
    if ( passes_mb ) th1d_cent_MBonly->Fill(centrality);

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
    th2d_BBCq_cent->Fill(bbc_nw_qw, centrality);


    // --------------
    // --- FVTX stuff
    // --------------

    float fvtxs_qxn[NHAR][10];//all layers then 0 1 2 3 A B
    float fvtxs_qyn[NHAR][10];
    float fvtxs_qw[10];

    float fvtxn_qxn[NHAR][10];//all layers then 0 1 2 3 A B
    float fvtxn_qyn[NHAR][10];
    float fvtxn_qw[10];

    for (int ilayer = 0; ilayer < 10; ilayer++)
    {
      for (int ih = 1; ih < NHAR; ih++)
      {
        fvtxs_qxn[ih][ilayer] = 0.0;
        fvtxs_qyn[ih][ilayer] = 0.0;

        fvtxn_qxn[ih][ilayer] = 0.0;
        fvtxn_qyn[ih][ilayer] = 0.0;
      } // loop over harmonics

      fvtxs_qw[ilayer] = 0.0;
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
          th3d_FVTXSclus_eta_zvrtx_lay->Fill(fvtx_eta, ZVTX, fvtx_layer);
          th3d_FVTXSclus_eta_zvrtx_cent->Fill(fvtx_eta, ZVTX, centrality);


          // --- all layers
          for (int ih = 1; ih < NHAR; ih++)
          {
            fvtxs_qxn[ih][0] += TMath::Cos((ih + 1) * phi);
            fvtxs_qyn[ih][0] += TMath::Sin((ih + 1) * phi);
            fvtxs_qxn[ih][fvtxsl_index + fvtx_layer] += TMath::Cos((ih + 1) * phi);
            fvtxs_qyn[ih][fvtxsl_index + fvtx_layer] += TMath::Sin((ih + 1) * phi);

          } // ih
          fvtxs_qw[fvtxsl_index + fvtx_layer] += 1;
          fvtxs_qw[0] += 1;

          // --- A or B
          if ( fvtx_eta > -2 )
          {
            for (int ih = 1; ih < NHAR; ih++)
            {
              fvtxs_qxn[ih][fvtxsa_index] += TMath::Cos((ih + 1) * phi);
              fvtxs_qyn[ih][fvtxsa_index] += TMath::Sin((ih + 1) * phi);
            }
            fvtxs_qw[fvtxsa_index] += 1;
          }
          else
          {
            for (int ih = 1; ih < NHAR; ih++)
            {
              fvtxs_qxn[ih][fvtxsa_index + 1] += TMath::Cos((ih + 1) * phi);
              fvtxs_qyn[ih][fvtxsa_index + 1] += TMath::Sin((ih + 1) * phi);
            }
            fvtxs_qw[fvtxsa_index + 1] += 1;
          }

        } // check on south

        // --- north side
        if ( d_FVTX_z[iclus] > 0 )
        {
          for ( int ih = 1; ih < NHAR; ih++ )
          {
            fvtxn_qxn[ih][0] += TMath::Cos((ih + 1) * phi);
            fvtxn_qyn[ih][0] += TMath::Sin((ih + 1) * phi);
          }
          fvtxn_qw[0] += 1;
        } // check on north

      } // loop over cluster
    } // check on clusters

    th1d_FVTXS_nclus->Fill(fvtxs_qw[0]);
    th1d_FVTXN_nclus->Fill(fvtxn_qw[0]);
    th2d_FVTXSnclus_cent->Fill(fvtxs_qw[0], centrality);
    th2d_FVTXnnclus_cent->Fill(fvtxn_qw[0], centrality);

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
      for (int ih = 1; ih < NHAR; ih++)
      {
        sumxy[ih][bbcs_index][0] = bbc_qxn[ih];
        sumxy[ih][bbcs_index][1] = bbc_qyn[ih];
        sumxy[ih][bbcs_index][2] = bbc_qw;
      }
    }

    if ( fvtx_clusters )
    {
      for (int ih = 1; ih < NHAR; ih++)
      {
        sumxy[ih][fvtxs_index][0] = fvtxs_qxn[ih][0];
        sumxy[ih][fvtxs_index][1] = fvtxs_qyn[ih][0];
        sumxy[ih][fvtxs_index][2] = fvtxs_qw[0];
        sumxy[ih][fvtxn_index][0] = fvtxn_qxn[ih][0];
        sumxy[ih][fvtxn_index][1] = fvtxn_qyn[ih][0];
        sumxy[ih][fvtxn_index][2] = fvtxn_qw[0];
        for (int il = 0; il < NFVTXLAY; il++)
        {
          sumxy[ih][fvtxsl_index + il][0] = fvtxs_qxn[ih][fvtxsl_index + il];
          sumxy[ih][fvtxsl_index + il][1] = fvtxs_qyn[ih][fvtxsl_index + il];
          sumxy[ih][fvtxsl_index + il][2] = fvtxs_qw[fvtxsl_index + il];
        }
        sumxy[ih][fvtxsa_index][0] = fvtxs_qxn[ih][fvtxsa_index];
        sumxy[ih][fvtxsa_index][1] = fvtxs_qyn[ih][fvtxsa_index];
        sumxy[ih][fvtxsa_index][2] = fvtxs_qw[fvtxsa_index];
        sumxy[ih][fvtxsa_index + 1][0] = fvtxs_qxn[ih][fvtxsa_index + 1];
        sumxy[ih][fvtxsa_index + 1][1] = fvtxs_qyn[ih][fvtxsa_index + 1];
        sumxy[ih][fvtxsa_index + 1][2] = fvtxs_qw[fvtxsa_index + 1];
      }
    } // check on clusters


    if ( DIAG )
    {
      cout << "bbc from node tree: " << d_Qx[5] << " " << d_Qy[5] << " " << d_Qw[5] << endl;
      cout << "bbc from me: " << bbc_qxn[1] << " " << bbc_qyn[1] << " " << bbc_qw << endl;

      cout << "fvtx raw: " << endl;
      cout << "from node tree: " << d_Qx[4] << " " << d_Qy[4] << " " << d_Qw[4] << endl;
      cout << "from clusters: " << fvtxs_qxn[1][0] << " " << fvtxs_qyn[1][0] << " " << fvtxs_qw[0] << endl;
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



    //
    // Add Q vector offsets.
    // Note, need to recalculate Qx and Qy to propagate flattening
    //
    for (int ih = 1; ih < NHAR; ih++)
    {
      float n = ih + 1.0;

      float qxoff, qyoff;

      // --- FVTXS all layers
      if ( ih == 1) // psi2
      {
        qxoff = 0;
        qyoff = qy2_offset_fvtxs[eidx][icent];
      }
      else if ( ih == 2) // psi3
      {
        qxoff = qx3_offset_fvtxs[eidx][icent];
        qyoff = 0;
      }
      else
      {
        qxoff = 0;
        qyoff = 0;
      }
      sumxy[ih][fvtxs_index][0] = cos(sumxy[ih][fvtxs_index][3] * n) + qxoff;
      sumxy[ih][fvtxs_index][1] = sin(sumxy[ih][fvtxs_index][3] * n) + qyoff;
      sumxy[ih][fvtxs_index][3] = atan2(sumxy[ih][fvtxs_index][1], sumxy[ih][fvtxs_index][0]) / n;

      // --- FVTXS A
      // (same offsets)
      sumxy[ih][fvtxsa_index][0] = cos(sumxy[ih][fvtxsa_index][3] * n) + qxoff;
      sumxy[ih][fvtxsa_index][1] = sin(sumxy[ih][fvtxsa_index][3] * n) + qyoff;
      sumxy[ih][fvtxsa_index][3] = atan2(sumxy[ih][fvtxsa_index][1], sumxy[ih][fvtxsa_index][0]) / n;


      // --- FVTXS B
      // (same offsets)
      sumxy[ih][fvtxsa_index + 1][0] = cos(sumxy[ih][fvtxsa_index + 1][3] * n) + qxoff;
      sumxy[ih][fvtxsa_index + 1][1] = sin(sumxy[ih][fvtxsa_index + 1][3] * n) + qyoff;
      sumxy[ih][fvtxsa_index + 1][3] = atan2(sumxy[ih][fvtxsa_index + 1][1], sumxy[ih][fvtxsa_index + 1][0]) / n;


      // --- FVTXS L
      for (int il = 0; il < NFVTXLAY; il++)
      {
        if ( ih == 1) // psi2
        {
          qxoff = 0;
          qyoff = qy2_offset_fvtxsl[eidx][il][icent];
        }
        else if ( ih == 2) // psi3
        {
          qxoff = qx3_offset_fvtxsl[eidx][il][icent];
          qyoff = 0;
        }
        else
        {
          qxoff = 0;
          qyoff = 0;
        }
        sumxy[ih][fvtxsl_index + il][0] = cos(sumxy[ih][fvtxsl_index + il][3] * n) + qxoff;
        sumxy[ih][fvtxsl_index + il][1] = sin(sumxy[ih][fvtxsl_index + il][3] * n) + qyoff;
        sumxy[ih][fvtxsl_index + il][3] = atan2(sumxy[ih][fvtxsl_index + il][1], sumxy[ih][fvtxsl_index + il][0]) / n;
      }


      // --- BBCS
      if ( ih == 1) // psi2
      {
        qxoff = 0;
        qyoff = qy2_offset_bbcs[eidx][icent];
      }
      else if ( ih == 2) // psi3
      {
        qxoff = qx3_offset_bbcs[eidx][icent];
        qyoff = 0;
      }
      else
      {
        qxoff = 0;
        qyoff = 0;
      }
      sumxy[ih][bbcs_index][0] = cos(sumxy[ih][bbcs_index][3] * n) + qxoff;
      sumxy[ih][bbcs_index][1] = sin(sumxy[ih][bbcs_index][3] * n) + qyoff;
      sumxy[ih][bbcs_index][3] = atan2(sumxy[ih][bbcs_index][1], sumxy[ih][bbcs_index][0]) / n;


    } // ih


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

    float bbcs_psin_docalib[NHAR] = {0};
    float fvtxs_psin_docalib[NHAR] = {0};
    float fvtxn_psin_docalib[NHAR] = {0};

    for (int ih = 1; ih < NHAR; ih++)
    {
      if ( bbc_pmts )
      {
        bbcs_psin_docalib[ih] = (sumxy[ih][bbcs_index][2] > 0) ? sumxy[ih][bbcs_index][3] : -9999.9;
      }
      if ( fvtx_clusters )
      {
        // ---
        fvtxs_psin_docalib[ih]  = (sumxy[ih][fvtxs_index][2] > 4)  ? sumxy[ih][fvtxs_index][3]  : -9999.9;

        // ---
        fvtxn_psin_docalib[ih]  = (sumxy[ih][fvtxn_index][2] > 4)  ? sumxy[ih][fvtxn_index][3]  : -9999.9;
      }
    }

    if ( verbosity > 0 && ( fvtxs_psin_docalib[1] < -999 || bbcs_psin_docalib[1] < -999 ) )
    {
      cout << "POSSIBLE ISSUE WITH EVENT PLANES!!!  ONE OR MORE IS -9999" << endl;
      cout << "BBC south event plane " << bbcs_psin_docalib[1] << endl;
      cout << "FVTX south event plane " << fvtxs_psin_docalib[1] << endl;
      cout << "FVTX north event plane " << fvtxn_psin_docalib[1] << endl;
      cout << "BBC charge " << bbc_qw << endl;
      cout << "centrality " << centrality << endl;
    }



    // ---
    // --- resolution histograms
    // ---
    for (int ih = 1; ih < NHAR; ih++)
    {
      // --- BBC and FVTX south
      tp1f_reso_BBC_FVTX[ih][icent][izvtx]->Fill(0.0, cos(float(ih + 1.) * (bbcs_psin_docalib[ih] - fvtxs_psin_docalib[ih])));

      // --- BBC and FVTX north
      tp1f_reso_BBC_FVTXN[ih][icent][izvtx]->Fill(0.0, cos(float(ih + 1.) * (bbcs_psin_docalib[ih] - fvtxn_psin_docalib[ih])));

      // --- FVTX south and north
      tp1f_reso_FVTXS_FVTXN[ih][icent][izvtx]->Fill(0.0, cos(float(ih + 1.) * (fvtxs_psin_docalib[ih] - fvtxn_psin_docalib[ih])));

      // --- BBC and FVTX S layers (also A, B)
      tp1f_reso_BBC_FVTXSA[ih][icent][izvtx]->Fill(0.0, cos(float(ih + 1.) * (bbcs_psin_docalib[ih] - sumxy[ih][fvtxsa_index][3])));
      tp1f_reso_BBC_FVTXSB[ih][icent][izvtx]->Fill(0.0, cos(float(ih + 1.) * (bbcs_psin_docalib[ih] - sumxy[ih][fvtxsa_index + 1][3])));
      for (int il = 0; il < NFVTXLAY; il++)
      {
        tp1f_reso_BBC_FVTXSL[ih][icent][izvtx][il]->Fill(0.0, cos(float(ih + 1.) * (bbcs_psin_docalib[ih] - sumxy[ih][fvtxsl_index + il][3])));
      }
    }
    // ----------------------------------------------------------------------------

    //start of cnt track loop
    if ( cnt_tracks )
    {
      for (int itrk = 0; itrk < d_ntrk; itrk++)
      {
        float px    = d_px[itrk];
        float py    = d_py[itrk];
        float pz    = d_pz[itrk];

        // to test effect of different beam angles
        // we first need to rotate back by the
        // 0.001 beam angle used when filling the trees
        float pxo = px;
        float pzo = pz;

        px = pxo * cos(0.001) + pzo * sin(0.001);
        pz = -1 * pxo * sin(0.001) + pzo * cos(0.001);

        pxo = px;
        pzo = pz;

        // now we can rotate based on what's in beam_angle above
        px = pxo * cos(-beam_angle) + pzo * sin(-beam_angle);
        pz = -1 * pxo * sin(-beam_angle) + pzo * cos(-beam_angle);

        // // reimplement rotation bug
        // px = pxo * cos(-beam_angle) + pzo * sin(-beam_angle);
        // pz = -1 * px * sin(-beam_angle) + pzo * cos(-beam_angle);

        // float charge    = d_charge[itrk];
        float pc3sdz    = d_pc3sdz[itrk];
        float pc3sdphi  = d_pc3sdphi[itrk];
        if ( tight_trkcuts )
          if ( fabs(pc3sdz) > 2.0 || fabs(pc3sdphi) > 2.0 ) continue;

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

        // --- Loop over harmonics
        for (int ih = 1; ih < NHAR; ih++)
        {
          float n = ih + 1.;

          cnt_both_cosnphi[ih][icent][izvtx]->Fill(pt_angle, cos(n * phi_angle));
          if ( dcarm == 1 ) cnt_west_cosnphi[ih][icent][izvtx]->Fill(pt_angle, cos(n * phi_angle));
          if ( dcarm == 0 ) cnt_east_cosnphi[ih][icent][izvtx]->Fill(pt_angle, cos(n * phi_angle));

          cnt_both_sinnphi[ih][icent][izvtx]->Fill(pt_angle, sin(n * phi_angle));
          if ( dcarm == 1 ) cnt_west_sinnphi[ih][icent][izvtx]->Fill(pt_angle, sin(n * phi_angle));
          if ( dcarm == 0 ) cnt_east_sinnphi[ih][icent][izvtx]->Fill(pt_angle, sin(n * phi_angle));


          cosnphi_eta_both[ih][icent][izvtx]->Fill(eta, cos(n * phi_angle));
          if ( dcarm == 1 ) cosnphi_eta_west[ih][icent][izvtx]->Fill(eta, cos(n * phi_angle));
          if ( dcarm == 0 ) cosnphi_eta_east[ih][icent][izvtx]->Fill(eta, cos(n * phi_angle));

          sinnphi_eta_both[ih][icent][izvtx]->Fill(eta, sin(n * phi_angle));
          if ( dcarm == 1 ) sinnphi_eta_west[ih][icent][izvtx]->Fill(eta, sin(n * phi_angle));
          if ( dcarm == 0 ) sinnphi_eta_east[ih][icent][izvtx]->Fill(eta, sin(n * phi_angle));


          // BBC event plane
          if ( bbc_pmts )
          {
            if (-4.0 < bbcs_psin_docalib[ih] && bbcs_psin_docalib[ih] < 4.0)
            {
              // --- nth harmonic
              double bbc_dphi = phi_angle - bbcs_psin_docalib[ih];
              double cosbbc_dphi = TMath::Cos(n * bbc_dphi);

              bbcs_vn_both_docalib[ih][icent][izvtx]->Fill(pt_angle, cosbbc_dphi);
              if ( dcarm == 1 ) bbcs_vn_west_docalib[ih][icent][izvtx]->Fill(pt_angle, cosbbc_dphi);
              if ( dcarm == 0 ) bbcs_vn_east_docalib[ih][icent][izvtx]->Fill(pt_angle, cosbbc_dphi);

              bbcs_vneta_both_docalib[ih][icent][izvtx]->Fill(eta, cosbbc_dphi);
              if ( dcarm == 1 ) bbcs_vneta_west_docalib[ih][icent][izvtx]->Fill(eta, cosbbc_dphi);
              if ( dcarm == 0 ) bbcs_vneta_east_docalib[ih][icent][izvtx]->Fill(eta, cosbbc_dphi);

              // --- resolutions
              if ( pt_angle > reso_lpt && pt_angle < reso_hpt )
              {
                tp1f_reso_BBC_CNT[ih][icent][izvtx]->Fill(0.0, cosbbc_dphi);
              }
            } // check on bbc EP
          } // check on tubes

          // FVTX Event Plane
          if ( fvtx_clusters )
          {

            // south all layers
            if ( -4.0 < fvtxs_psin_docalib[ih] && fvtxs_psin_docalib[ih] < 4.0 )
            {

              // --- nth harmonic
              double fvtxs_dphi = phi_angle - fvtxs_psin_docalib[ih];
              double cosfvtxs_dphi = TMath::Cos(n * fvtxs_dphi);

              fvtxs_vn_both_docalib[ih][icent][izvtx]->Fill(pt_angle, cosfvtxs_dphi);
              if ( dcarm == 1 ) fvtxs_vn_west_docalib[ih][icent][izvtx]->Fill(pt_angle, cosfvtxs_dphi);
              if ( dcarm == 0 ) fvtxs_vn_east_docalib[ih][icent][izvtx]->Fill(pt_angle, cosfvtxs_dphi);

              fvtxs_vneta_both_docalib[ih][icent][izvtx]->Fill(eta, cosfvtxs_dphi);
              if ( dcarm == 1 ) fvtxs_vneta_west_docalib[ih][icent][izvtx]->Fill(eta, cosfvtxs_dphi);
              if ( dcarm == 0 ) fvtxs_vneta_east_docalib[ih][icent][izvtx]->Fill(eta, cosfvtxs_dphi);

              // --- ep resolutions
              if ( pt_angle > reso_lpt && pt_angle < reso_hpt )
              {
                tp1f_reso_CNT_FVTX[ih][icent][izvtx]->Fill(0.0, cosfvtxs_dphi);
              }

            } // check on ep south all layers


            // north all layers
            if ( -4.0 < fvtxn_psin_docalib[ih] && fvtxn_psin_docalib[ih] < 4.0 )
            {

              // --- 2nd harmonic
              double fvtxn_dphi = phi_angle - fvtxn_psin_docalib[ih];
              double cosfvtxn_dphi = TMath::Cos(n * fvtxn_dphi);

              fvtxn_vn_both_docalib[ih][icent][izvtx]->Fill(pt_angle, cosfvtxn_dphi);
              if ( dcarm == 1 ) fvtxn_vn_west_docalib[ih][icent][izvtx]->Fill(pt_angle, cosfvtxn_dphi);
              if ( dcarm == 0 ) fvtxn_vn_east_docalib[ih][icent][izvtx]->Fill(pt_angle, cosfvtxn_dphi);

              fvtxn_vneta_both_docalib[ih][icent][izvtx]->Fill(eta, cosfvtxn_dphi);
              if ( dcarm == 1 ) fvtxn_vneta_west_docalib[ih][icent][izvtx]->Fill(eta, cosfvtxn_dphi);
              if ( dcarm == 0 ) fvtxn_vneta_east_docalib[ih][icent][izvtx]->Fill(eta, cosfvtxn_dphi);


              // --- ep resolutions
              if ( pt_angle > reso_lpt && pt_angle < reso_hpt )
              {
                tp1f_reso_CNT_FVTXN[ih][icent][izvtx]->Fill(0.0, cosfvtxn_dphi);
              }

            } // check on ep north all layers


            // south A layers
            if ( -4.0 < sumxy[ih][fvtxsa_index][3] && sumxy[ih][fvtxsa_index][3] < 4.0 )
            {

              // --- nth harmonic
              double fvtxsa_dphi = phi_angle - sumxy[ih][fvtxsa_index][3];
              double cosfvtxsa_dphi = TMath::Cos(n * fvtxsa_dphi);

              fvtxsa_vn_both_docalib[ih][icent][izvtx]->Fill(pt_angle, cosfvtxsa_dphi);
              if ( dcarm == 1 ) fvtxsa_vn_west_docalib[ih][icent][izvtx]->Fill(pt_angle, cosfvtxsa_dphi);
              if ( dcarm == 0 ) fvtxsa_vn_east_docalib[ih][icent][izvtx]->Fill(pt_angle, cosfvtxsa_dphi);

              // --- ep resolutions
              if ( pt_angle > reso_lpt && pt_angle < reso_hpt )
              {
                tp1f_reso_CNT_FVTXSA[ih][icent][izvtx]->Fill(0.0, cosfvtxsa_dphi);
              }

            } // check on ep south A layers


            // south B layers
            if ( -4.0 < sumxy[ih][fvtxsa_index + 1][3] && sumxy[ih][fvtxsa_index + 1][3] < 4.0 )
            {

              // --- nth harmonic
              double fvtxsb_dphi = phi_angle - sumxy[ih][fvtxsa_index + 1][3];
              double cosfvtxsb_dphi = TMath::Cos(n * fvtxsb_dphi);

              fvtxsb_vn_both_docalib[ih][icent][izvtx]->Fill(pt_angle, cosfvtxsb_dphi);
              if ( dcarm == 1 ) fvtxsb_vn_west_docalib[ih][icent][izvtx]->Fill(pt_angle, cosfvtxsb_dphi);
              if ( dcarm == 0 ) fvtxsb_vn_east_docalib[ih][icent][izvtx]->Fill(pt_angle, cosfvtxsb_dphi);

              // --- ep resolutions
              if ( pt_angle > reso_lpt && pt_angle < reso_hpt )
              {
                tp1f_reso_CNT_FVTXSB[ih][icent][izvtx]->Fill(0.0, cosfvtxsb_dphi);
              }

            } // check on ep south B layers

            // south layers
            for (int il = 0; il < NFVTXLAY; il++)
            {
              if ( -4.0 < sumxy[ih][fvtxsl_index + il][3] && sumxy[ih][fvtxsl_index + il][3] < 4.0 )
              {

                // --- nth harmonic
                double fvtxsl_dphi = phi_angle - sumxy[ih][fvtxsl_index + il][3];
                double cosfvtxsl_dphi = TMath::Cos(n * fvtxsl_dphi);

                fvtxsl_vn_both_docalib[ih][icent][izvtx][il]->Fill(pt_angle, cosfvtxsl_dphi);
                if ( dcarm == 1 ) fvtxsl_vn_west_docalib[ih][icent][izvtx][il]->Fill(pt_angle, cosfvtxsl_dphi);
                if ( dcarm == 0 ) fvtxsl_vn_east_docalib[ih][icent][izvtx][il]->Fill(pt_angle, cosfvtxsl_dphi);

                // --- ep resolutions
                if ( pt_angle > reso_lpt && pt_angle < reso_hpt )
                {
                  tp1f_reso_CNT_FVTXSL[ih][icent][izvtx][il]->Fill(0.0, cosfvtxsl_dphi);
                }

              } // check on ep south layers
            } // il

          } // check on fvtx clusters
        } // loop over harmonics
      } // loop over cnt tracks

    } // check on cnt tracks


    if ( fvtx_tracks )
    {
      int nfvtxts = 0;
      int nfvtxtn = 0;
      // loop over fvtx tracks
      for ( int i = 0; i < nfvtxt; ++i )
      {
        float phi = fphi[i];
        float eta = feta[i];
        float dcax = fdcax[i];
        float dcay = fdcay[i];
        float chisq = fchisq[i];

        if ( fabs(dcax - x_off) > 2.0 || fabs(dcay - y_off) > 2.0 ) continue;
        if ( tight_trkcuts )
          if ( fabs(dcax - x_off) > 0.5 || fabs(dcay - y_off) > 0.5 ) continue;
        if ( chisq > 5 ) continue;
        int ns = 0;
        int nn = 0;
        if ( eta > 0 )
        {
          nfvtxtn++;
          nn = 1;
        }
        else if ( eta < 0 )
        {
          nfvtxts++;
          ns = 1;
        }
        // cout << "eta is " << eta << " and ns is " << ns << " and nn is " << nn << endl;

        bool west = (phi > -pi / 2. && phi < pi / 2.);


        // --- loop over harmonics
        for (int ih = 1; ih < NHAR; ih++)
        {
          float n = ih + 1.;


          cosnphi_eta_both[ih][icent][izvtx]->Fill(eta, cos(n * phi));
          if ( west ) cosnphi_eta_west[ih][icent][izvtx]->Fill(eta, cos(n * phi));
          if ( !west ) cosnphi_eta_east[ih][icent][izvtx]->Fill(eta, cos(n * phi));

          sinnphi_eta_both[ih][icent][izvtx]->Fill(eta, sin(n * phi));
          if ( west ) sinnphi_eta_west[ih][icent][izvtx]->Fill(eta, sin(n * phi));
          if ( !west ) sinnphi_eta_east[ih][icent][izvtx]->Fill(eta, sin(n * phi));




          // BBC event plane
          if ( bbc_pmts )
          {
            if (-4.0 < bbcs_psin_docalib[ih] && bbcs_psin_docalib[ih] < 4.0)
            {
              // --- nth harmonic
              double bbcs_dphi = phi - bbcs_psin_docalib[ih];
              double cosbbcs_dphi = TMath::Cos(n * bbcs_dphi);

              bbcs_vneta_both_docalib[ih][icent][izvtx]->Fill(eta, cosbbcs_dphi);
              if ( west ) bbcs_vneta_west_docalib[ih][icent][izvtx]->Fill(eta, cosbbcs_dphi);
              if ( !west ) bbcs_vneta_east_docalib[ih][icent][izvtx]->Fill(eta, cosbbcs_dphi);
            } // check on bbc EP
          } // check on tubes

          // FVTX Event Plane
          if ( fvtx_clusters )
          {

            // south all layers
            if ( -4.0 < fvtxs_psin_docalib[ih] && fvtxs_psin_docalib[ih] < 4.0 )
            {

              // --- nth harmonic
              double fvtxs_dphi = phi - fvtxs_psin_docalib[ih];
              double cosfvtxs_dphi = TMath::Cos(n * fvtxs_dphi);

              fvtxs_vneta_both_docalib[ih][icent][izvtx]->Fill(eta, cosfvtxs_dphi);
              if ( west ) fvtxs_vneta_west_docalib[ih][icent][izvtx]->Fill(eta, cosfvtxs_dphi);
              if ( !west ) fvtxs_vneta_east_docalib[ih][icent][izvtx]->Fill(eta, cosfvtxs_dphi);
            } // check on ep south all layers


            // north all layers
            if ( -4.0 < fvtxn_psin_docalib[ih] && fvtxn_psin_docalib[ih] < 4.0 )
            {

              // --- nth harmonic
              double fvtxn_dphi = phi - fvtxn_psin_docalib[ih];
              double cosfvtxn_dphi = TMath::Cos(n * fvtxn_dphi);

              fvtxn_vneta_both_docalib[ih][icent][izvtx]->Fill(eta, cosfvtxn_dphi);
              if ( west ) fvtxn_vneta_west_docalib[ih][icent][izvtx]->Fill(eta, cosfvtxn_dphi);
              if ( !west ) fvtxn_vneta_east_docalib[ih][icent][izvtx]->Fill(eta, cosfvtxn_dphi);
            } // check on ep north all layers

          } // check on fvtx clusters

        } // harmonics

      } // loop over fvtx tracks

      th2d_FVTXSntrk_cent->Fill(nfvtxts, centrality);
      th2d_FVTXnntrk_cent->Fill(nfvtxtn, centrality);

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
