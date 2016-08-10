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




bool DIAG = false;


int get_fvtx_layer(float);
void initialize_pmt_position();
int get_pmt_layer(int);
void flatten(int, int);
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

  flatten(run,1);
  flatten(run,2);
  flatten(run,3);

  return 0;

}

// -----------------------------------------------------------------
void flatten(int runNumber, int rp_recal_pass)
{

  int verbosity = 0;

  char outFile1[300];
  sprintf(outFile1,"%s%d%s","output/hist_",runNumber,".root");

  char outFile2[100];
  sprintf(outFile2,"%s%d%s","output/hrp_",runNumber,".root");

  cout << "runNumber = " << runNumber << " rp_recal_pass = " << rp_recal_pass << endl;

  cout << "v2 histogram output file: " << outFile1 << endl;

  char calibfile[500];
  sprintf(calibfile,"output/flattening_%d_%d.dat",runNumber,rp_recal_pass-1);

  char filename[500];

  // --- get the number of files for this run number
  string pipe_out = (string) gSystem->GetFromPipe(Form("ls input/tree_%010d_*.root | grep -c r",runNumber));
  int nfiles = 0;
  nfiles = atoi(pipe_out.c_str());
  cout<<"nfiles: "<<nfiles<<endl;
  if(nfiles==0) return;

  // --- make a new TChain for the tree
  TChain *ntp_event_chain = new TChain("ntp_event");
  for ( int ifile = 0; ifile < nfiles; ++ifile )
    {
      sprintf(filename,"input/tree_%010d_%04d.root",runNumber,ifile);
      cout<<"adding to tchain: "<<filename<<endl;
      ntp_event_chain->Add(filename);
    }

  // ---

  TFile *phi_weight_file = TFile::Open(Form("SpecialProjects/WeightFiles/weight2d_run%d.root", runNumber)); // COME BACK HERE AND HAVE A LOOK
  if (!phi_weight_file)
    {
      cout << "ERROR could not open phi weight file: " << endl;
      return;
    }

  TH1D* th1d_fvtxs_phi_weight[10][5];
  for ( int i = 0; i < 10; ++i )
    {
      th1d_fvtxs_phi_weight[i][0] = (TH1D*)phi_weight_file->Get(Form("th1d_weight_fvtxs_zvtx%d_clus_phi_GC",i)); // COME BACK HERE AND HAVE A LOOK
      th1d_fvtxs_phi_weight[i][1] = (TH1D*)phi_weight_file->Get(Form("th1d_weight_fvtxs0_zvtx%d_clus_phi_GC",i)); // COME BACK HERE AND HAVE A LOOK
      th1d_fvtxs_phi_weight[i][2] = (TH1D*)phi_weight_file->Get(Form("th1d_weight_fvtxs1_zvtx%d_clus_phi_GC",i)); // COME BACK HERE AND HAVE A LOOK
      th1d_fvtxs_phi_weight[i][3] = (TH1D*)phi_weight_file->Get(Form("th1d_weight_fvtxs2_zvtx%d_clus_phi_GC",i)); // COME BACK HERE AND HAVE A LOOK
      th1d_fvtxs_phi_weight[i][4] = (TH1D*)phi_weight_file->Get(Form("th1d_weight_fvtxs3_zvtx%d_clus_phi_GC",i)); // COME BACK HERE AND HAVE A LOOK
      cout << "memory address of th1d_fvtxs_phi_weight[i][0] is " << th1d_fvtxs_phi_weight[i][0] << endl;
      cout << "memory address of th1d_fvtxs_phi_weight[i][1] is " << th1d_fvtxs_phi_weight[i][1] << endl;
      cout << "memory address of th1d_fvtxs_phi_weight[i][2] is " << th1d_fvtxs_phi_weight[i][2] << endl;
      cout << "memory address of th1d_fvtxs_phi_weight[i][3] is " << th1d_fvtxs_phi_weight[i][3] << endl;
      cout << "memory address of th1d_fvtxs_phi_weight[i][4] is " << th1d_fvtxs_phi_weight[i][4] << endl;
    }


  // ---

  cout << "Initalizing PMT positions for the BBC" << endl;

  initialize_pmt_position();

  bool fvtx_clusters = true;
  bool bbc_pmts      = true;
  bool cnt_tracks    = true;




  int n_angle_config = 1;
  // --- see below...
  int south_bbc_angle = 2;                     // 2
  int south_fvtx_angle = n_angle_config+2;     // 3
  int south_fvtx_0_angle = 2*n_angle_config+2; // 4
  int south_fvtx_1_angle = 3*n_angle_config+2; // 5
  int south_fvtx_2_angle = 4*n_angle_config+2; // 6
  int south_fvtx_3_angle = 5*n_angle_config+2; // 7
  int north_fvtx_angle = 6*n_angle_config+2;   // 8
  int north_fvtx_0_angle = 7*n_angle_config+2; // 9
  int north_fvtx_1_angle = 8*n_angle_config+2; // 10
  int north_fvtx_2_angle = 9*n_angle_config+2; // 11
  int north_fvtx_3_angle = 10*n_angle_config+2; // 12


  float pi = acos(-1.0);

  //------------------------------------------------------------//
  //                                                            //
  //       Initializing Calibration Arrays & Histograms         //
  //                                                            //
  //------------------------------------------------------------//

  // --- problems with these array dimensions??  lots of compile errors

  cout << "Lots of arrays and stuff" << endl;

  // --- lots of comments needed here
  // --- flattening parameters output to file
  TH2D     *qx[NMUL][NHAR][NDET]; // Q vector x component, Q vector vs z_vertex bin
  TH2D     *qy[NMUL][NHAR][NDET]; // Q vector y component, Q vector vs z_vertex bin
  TProfile *ave[NMUL][NZPS][NHAR][NDET]; // average Psi
  TProfile *flt[NMUL][NZPS][NHAR][NDET]; // flattening parameters
  TH2D     *dis[NMUL][NHAR][NDET]; // displacement?  function of z_vertex bin

  TH2D     *psi_bf[NMUL][NHAR][NDET];
  TH2D     *psi_af[NMUL][NHAR][NDET];

  // flattening parameters read in from file
  float    mean[NMUL][NZPS][NHAR][NDET][2]; // mean of Psi distribution (???)
  float    widt[NMUL][NZPS][NHAR][NDET][2]; // width of Psi distribution (???)
  float    four[NMUL][NZPS][NHAR][NDET][2][NORD]; // ?

  // ---
  TH1D* th1d_BBC_charge = new TH1D("th1d_BBC_charge","",200,-0.5,199.5);
  TH1D* th1d_FVTX_nclus = new TH1D("th1d_FVTX_nclus","",200,-0.5,1999.5);
  TH2D* th2d_qBBC_nFVTX = new TH2D("th2d_qBBC_nFVTX","",200,-0.5,199.5,200,-0.5,1999.5);
  TH1D* th1d_FVTXS_nclus = new TH1D("th1d_FVTXS_nclus","",200,-0.5,1999.5);
  TH1D* th1d_FVTXN_nclus = new TH1D("th1d_FVTXN_nclus","",200,-0.5,1999.5);

  // --- come back here to make rings
  TH1D* th1d_bbc_charge_phi = new TH1D("th1d_bbc_charge_phi","",50,-pi,pi);
  TH1D* th1d_bbc0_charge_phi = new TH1D("th1d_bbc0_charge_phi","",50,-pi,pi);
  TH1D* th1d_bbc1_charge_phi = new TH1D("th1d_bbc1_charge_phi","",50,-pi,pi);
  TH1D* th1d_bbc2_charge_phi = new TH1D("th1d_bbc2_charge_phi","",50,-pi,pi);
  TH1D* th1d_bbc3_charge_phi = new TH1D("th1d_bbc3_charge_phi","",50,-pi,pi);
  TH1D* th1d_bbc4_charge_phi = new TH1D("th1d_bbc4_charge_phi","",50,-pi,pi);

  TH1D* th1d_fvtxs_clus_phi = new TH1D("th1d_fvtxs_clus_phi","",50,-pi,pi);
  TH1D* th1d_fvtxs0_clus_phi = new TH1D("th1d_fvtxs0_clus_phi","",50,-pi,pi);
  TH1D* th1d_fvtxs1_clus_phi = new TH1D("th1d_fvtxs1_clus_phi","",50,-pi,pi);
  TH1D* th1d_fvtxs2_clus_phi = new TH1D("th1d_fvtxs2_clus_phi","",50,-pi,pi);
  TH1D* th1d_fvtxs3_clus_phi = new TH1D("th1d_fvtxs3_clus_phi","",50,-pi,pi);

  TH1D* th1d_fvtxn_clus_phi = new TH1D("th1d_fvtxn_clus_phi","",50,-pi,pi);
  TH1D* th1d_fvtxn0_clus_phi = new TH1D("th1d_fvtxn0_clus_phi","",50,-pi,pi);
  TH1D* th1d_fvtxn1_clus_phi = new TH1D("th1d_fvtxn1_clus_phi","",50,-pi,pi);
  TH1D* th1d_fvtxn2_clus_phi = new TH1D("th1d_fvtxn2_clus_phi","",50,-pi,pi);
  TH1D* th1d_fvtxn3_clus_phi = new TH1D("th1d_fvtxn3_clus_phi","",50,-pi,pi);

  // --- event plane resolution
  TProfile* tp1f_reso2_BBC_CNT = new TProfile("tp1f_reso2_BBC_CNT","",1,-0.5,0.5,-1e6,1e6,"");
  TProfile* tp1f_reso2_BBC_FVTX = new TProfile("tp1f_reso2_BBC_FVTX","",1,-0.5,0.5,-1e6,1e6,"");
  TProfile* tp1f_reso2_CNT_FVTX = new TProfile("tp1f_reso2_CNT_FVTX","",1,-0.5,0.5,-1e6,1e6,"");
  TProfile* tp1f_reso3_BBC_CNT = new TProfile("tp1f_reso3_BBC_CNT","",1,-0.5,0.5,-1e6,1e6,"");
  TProfile* tp1f_reso3_BBC_FVTX = new TProfile("tp1f_reso3_BBC_FVTX","",1,-0.5,0.5,-1e6,1e6,"");
  TProfile* tp1f_reso3_CNT_FVTX = new TProfile("tp1f_reso3_CNT_FVTX","",1,-0.5,0.5,-1e6,1e6,"");

  TH1D* th1d_reso2_BBC_CNT = new TH1D("th1d_reso2_BBC_CNT","",220,-1.1,1.1);
  TH1D* th1d_reso2_BBC_FVTX = new TH1D("th1d_reso2_BBC_FVTX","",220,-1.1,1.1);
  TH1D* th1d_reso2_CNT_FVTX = new TH1D("th1d_reso2_CNT_FVTX","",220,-1.1,1.1);
  TH1D* th1d_reso3_BBC_CNT = new TH1D("th1d_reso3_BBC_CNT","",220,-1.1,1.1);
  TH1D* th1d_reso3_BBC_FVTX = new TH1D("th1d_reso3_BBC_FVTX","",220,-1.1,1.1);
  TH1D* th1d_reso3_CNT_FVTX = new TH1D("th1d_reso3_CNT_FVTX","",220,-1.1,1.1);

  TH1D* th1d_dreso2_BBC_CNT = new TH1D("th1d_dreso2_BBC_CNT","",252,-6.3,6.3);
  TH1D* th1d_dreso2_BBC_FVTX = new TH1D("th1d_dreso2_BBC_FVTX","",252,-6.3,6.3);
  TH1D* th1d_dreso2_CNT_FVTX = new TH1D("th1d_dreso2_CNT_FVTX","",252,-6.3,6.3);
  TH1D* th1d_dreso3_BBC_CNT = new TH1D("th1d_dreso3_BBC_CNT","",252,-6.3,6.3);
  TH1D* th1d_dreso3_BBC_FVTX = new TH1D("th1d_dreso3_BBC_FVTX","",252,-6.3,6.3);
  TH1D* th1d_dreso3_CNT_FVTX = new TH1D("th1d_dreso3_CNT_FVTX","",252,-6.3,6.3);

  // --- new EP resolution histograms using FVTX north
  TProfile* tp1f_reso2_BBC_FVTXN = new TProfile("tp1f_reso2_BBC_FVTXN","",1,-0.5,0.5,-1e6,1e6,"");
  TProfile* tp1f_reso3_BBC_FVTXN = new TProfile("tp1f_reso3_BBC_FVTXN","",1,-0.5,0.5,-1e6,1e6,"");
  TProfile* tp1f_reso2_CNT_FVTXN = new TProfile("tp1f_reso2_CNT_FVTXN","",1,-0.5,0.5,-1e6,1e6,"");
  TProfile* tp1f_reso3_CNT_FVTXN = new TProfile("tp1f_reso3_CNT_FVTXN","",1,-0.5,0.5,-1e6,1e6,"");
  TProfile* tp1f_reso2_FVTXS_FVTXN = new TProfile("tp1f_reso2_FVTXS_FVTXN","",1,-0.5,0.5,-1e6,1e6,"");
  TProfile* tp1f_reso3_FVTXS_FVTXN = new TProfile("tp1f_reso3_FVTXS_FVTXN","",1,-0.5,0.5,-1e6,1e6,"");

  TH1D* th1d_reso2_BBC_FVTXN = new TH1D("th1d_reso2_BBC_FVTXN","",220,-1.1,1.1);
  TH1D* th1d_reso3_BBC_FVTXN = new TH1D("th1d_reso3_BBC_FVTXN","",220,-1.1,1.1);
  TH1D* th1d_reso2_CNT_FVTXN = new TH1D("th1d_reso2_CNT_FVTXN","",220,-1.1,1.1);
  TH1D* th1d_reso3_CNT_FVTXN = new TH1D("th1d_reso3_CNT_FVTXN","",220,-1.1,1.1);
  TH1D* th1d_reso2_FVTXS_FVTXN = new TH1D("th1d_reso2_FVTXS_FVTXN","",220,-1.1,1.1);
  TH1D* th1d_reso3_FVTXS_FVTXN = new TH1D("th1d_reso3_FVTXS_FVTXN","",220,-1.1,1.1);

  TH1D* th1d_dreso2_BBC_FVTXN = new TH1D("th1d_dreso2_BBC_FVTXN","",252,-6.3,6.3);
  TH1D* th1d_dreso3_BBC_FVTXN = new TH1D("th1d_dreso3_BBC_FVTXN","",252,-6.3,6.3);
  TH1D* th1d_dreso2_CNT_FVTXN = new TH1D("th1d_dreso2_CNT_FVTXN","",252,-6.3,6.3);
  TH1D* th1d_dreso3_CNT_FVTXN = new TH1D("th1d_dreso3_CNT_FVTXN","",252,-6.3,6.3);
  TH1D* th1d_dreso2_FVTXS_FVTXN = new TH1D("th1d_dreso2_FVTXS_FVTXN","",252,-6.3,6.3);
  TH1D* th1d_dreso3_FVTXS_FVTXN = new TH1D("th1d_dreso3_FVTXS_FVTXN","",252,-6.3,6.3);

  // --- even more stuff, for v4{psi2}

  TProfile* tp1f_reso42_BBC_CNT = new TProfile("tp1f_reso42_BBC_CNT","",1,-0.5,0.5,-1e6,1e6,"");
  TProfile* tp1f_reso42_BBC_FVTX = new TProfile("tp1f_reso42_BBC_FVTX","",1,-0.5,0.5,-1e6,1e6,"");
  TProfile* tp1f_reso42_CNT_FVTX = new TProfile("tp1f_reso42_CNT_FVTX","",1,-0.5,0.5,-1e6,1e6,"");
  TProfile* tp1f_reso42_BBC_FVTXN = new TProfile("tp1f_reso42_BBC_FVTXN","",1,-0.5,0.5,-1e6,1e6,"");
  TProfile* tp1f_reso42_CNT_FVTXN = new TProfile("tp1f_reso42_CNT_FVTXN","",1,-0.5,0.5,-1e6,1e6,"");
  TProfile* tp1f_reso42_FVTXS_FVTXN = new TProfile("tp1f_reso42_FVTXS_FVTXN","",1,-0.5,0.5,-1e6,1e6,"");

  // --- event plane resolutions with offset q-vectors

  // --- event plane resolution
  TProfile* tp1f_os_reso2_BBC_CNT = new TProfile("tp1f_os_reso2_BBC_CNT","",1,-0.5,0.5,-1e6,1e6,"");
  TProfile* tp1f_os_reso2_BBC_FVTX = new TProfile("tp1f_os_reso2_BBC_FVTX","",1,-0.5,0.5,-1e6,1e6,"");
  TProfile* tp1f_os_reso2_CNT_FVTX = new TProfile("tp1f_os_reso2_CNT_FVTX","",1,-0.5,0.5,-1e6,1e6,"");
  TProfile* tp1f_os_reso3_BBC_CNT = new TProfile("tp1f_os_reso3_BBC_CNT","",1,-0.5,0.5,-1e6,1e6,"");
  TProfile* tp1f_os_reso3_BBC_FVTX = new TProfile("tp1f_os_reso3_BBC_FVTX","",1,-0.5,0.5,-1e6,1e6,"");
  TProfile* tp1f_os_reso3_CNT_FVTX = new TProfile("tp1f_os_reso3_CNT_FVTX","",1,-0.5,0.5,-1e6,1e6,"");

  TH1D* th1d_os_reso2_BBC_CNT = new TH1D("th1d_os_reso2_BBC_CNT","",220,-1.1,1.1);
  TH1D* th1d_os_reso2_BBC_FVTX = new TH1D("th1d_os_reso2_BBC_FVTX","",220,-1.1,1.1);
  TH1D* th1d_os_reso2_CNT_FVTX = new TH1D("th1d_os_reso2_CNT_FVTX","",220,-1.1,1.1);
  TH1D* th1d_os_reso3_BBC_CNT = new TH1D("th1d_os_reso3_BBC_CNT","",220,-1.1,1.1);
  TH1D* th1d_os_reso3_BBC_FVTX = new TH1D("th1d_os_reso3_BBC_FVTX","",220,-1.1,1.1);
  TH1D* th1d_os_reso3_CNT_FVTX = new TH1D("th1d_os_reso3_CNT_FVTX","",220,-1.1,1.1);

  TH1D* th1d_os_dreso2_BBC_CNT = new TH1D("th1d_os_dreso2_BBC_CNT","",252,-6.3,6.3);
  TH1D* th1d_os_dreso2_BBC_FVTX = new TH1D("th1d_os_dreso2_BBC_FVTX","",252,-6.3,6.3);
  TH1D* th1d_os_dreso2_CNT_FVTX = new TH1D("th1d_os_dreso2_CNT_FVTX","",252,-6.3,6.3);
  TH1D* th1d_os_dreso3_BBC_CNT = new TH1D("th1d_os_dreso3_BBC_CNT","",252,-6.3,6.3);
  TH1D* th1d_os_dreso3_BBC_FVTX = new TH1D("th1d_os_dreso3_BBC_FVTX","",252,-6.3,6.3);
  TH1D* th1d_os_dreso3_CNT_FVTX = new TH1D("th1d_os_dreso3_CNT_FVTX","",252,-6.3,6.3);

  // --- new EP resolution histograms using FVTX north
  TProfile* tp1f_os_reso2_BBC_FVTXN = new TProfile("tp1f_os_reso2_BBC_FVTXN","",1,-0.5,0.5,-1e6,1e6,"");
  TProfile* tp1f_os_reso3_BBC_FVTXN = new TProfile("tp1f_os_reso3_BBC_FVTXN","",1,-0.5,0.5,-1e6,1e6,"");
  TProfile* tp1f_os_reso2_CNT_FVTXN = new TProfile("tp1f_os_reso2_CNT_FVTXN","",1,-0.5,0.5,-1e6,1e6,"");
  TProfile* tp1f_os_reso3_CNT_FVTXN = new TProfile("tp1f_os_reso3_CNT_FVTXN","",1,-0.5,0.5,-1e6,1e6,"");
  TProfile* tp1f_os_reso2_FVTXS_FVTXN = new TProfile("tp1f_os_reso2_FVTXS_FVTXN","",1,-0.5,0.5,-1e6,1e6,"");
  TProfile* tp1f_os_reso3_FVTXS_FVTXN = new TProfile("tp1f_os_reso3_FVTXS_FVTXN","",1,-0.5,0.5,-1e6,1e6,"");

  TH1D* th1d_os_reso2_BBC_FVTXN = new TH1D("th1d_os_reso2_BBC_FVTXN","",220,-1.1,1.1);
  TH1D* th1d_os_reso3_BBC_FVTXN = new TH1D("th1d_os_reso3_BBC_FVTXN","",220,-1.1,1.1);
  TH1D* th1d_os_reso2_CNT_FVTXN = new TH1D("th1d_os_reso2_CNT_FVTXN","",220,-1.1,1.1);
  TH1D* th1d_os_reso3_CNT_FVTXN = new TH1D("th1d_os_reso3_CNT_FVTXN","",220,-1.1,1.1);
  TH1D* th1d_os_reso2_FVTXS_FVTXN = new TH1D("th1d_os_reso2_FVTXS_FVTXN","",220,-1.1,1.1);
  TH1D* th1d_os_reso3_FVTXS_FVTXN = new TH1D("th1d_os_reso3_FVTXS_FVTXN","",220,-1.1,1.1);

  TH1D* th1d_os_dreso2_BBC_FVTXN = new TH1D("th1d_os_dreso2_BBC_FVTXN","",252,-6.3,6.3);
  TH1D* th1d_os_dreso3_BBC_FVTXN = new TH1D("th1d_os_dreso3_BBC_FVTXN","",252,-6.3,6.3);
  TH1D* th1d_os_dreso2_CNT_FVTXN = new TH1D("th1d_os_dreso2_CNT_FVTXN","",252,-6.3,6.3);
  TH1D* th1d_os_dreso3_CNT_FVTXN = new TH1D("th1d_os_dreso3_CNT_FVTXN","",252,-6.3,6.3);
  TH1D* th1d_os_dreso2_FVTXS_FVTXN = new TH1D("th1d_os_dreso2_FVTXS_FVTXN","",252,-6.3,6.3);
  TH1D* th1d_os_dreso3_FVTXS_FVTXN = new TH1D("th1d_os_dreso3_FVTXS_FVTXN","",252,-6.3,6.3);



  cout << "Making TProfile histograms" << endl;

  // --- profile histograms for average of Psi and flattening parameters
  char name[200];
  int icent = 0;
  int ic = icent;
  for (int iz=0; iz<NZPS; iz++)
    {
      for (int ih=1; ih<NHAR; ih++)
        {
          for (int id=0; id<NDET; id++)
            {
              sprintf(name,"ave_%d_%d_%d_%d",ic,iz,ih,id);
              ave[ic][iz][ih][id] = new TProfile(name,name,4,-0.5,3.5,-10.1,10.1,"S");//for SMD -1.1,1.1

              sprintf(name,"flt_%d_%d_%d_%d",ic,iz,ih,id);
              flt[ic][iz][ih][id] = new TProfile(name,name,4*NORD,-0.5,NORD*4.0-0.5,-1.1,1.1);
            } // loop over ndetectors
        } // loop over harmonics
    } // loop over z_vertex bins

  // --- TH2D histograms for Q vector components
  for (int ih=1; ih<NHAR; ih++)
    {
      for (int id=0; id<NDET; id++)
        {
          sprintf(name,"dis_%d_%d_%d",ic,ih,id);
          dis[ic][ih][id] = new TH2D(name,name,NZPS*3,-0.5,NZPS*3.0-0.5,50,-pi,pi);

          sprintf(name,"qx_%d_%d_%d",ic,ih,id);
          qx[ic][ih][id] = new TH2D(name,name,NZPS*3,-0.5,NZPS*3.0-0.5, 220,-4.1,4.1);

          sprintf(name,"qy_%d_%d_%d",ic,ih,id);
          qy[ic][ih][id] = new TH2D(name,name,NZPS*3,-0.5,NZPS*3.0-0.5, 220,-4.1,4.1);

          sprintf(name,"psi_bf_%d_%d_%d",ic,ih,id);
          psi_bf[ic][ih][id] = new TH2D(name,name,NZPS*3,-0.5,NZPS*3.0-0.5, 220,-4.1,4.1);

          sprintf(name,"psi_af_%d_%d_%d",ic,ih,id);
          psi_af[ic][ih][id] = new TH2D(name,name,NZPS*3,-0.5,NZPS*3.0-0.5, 220,-4.1,4.1);

        } // loop over detectors
    } // loop over harmonics

  cout << "Initalizing calibration parameters to zero" << endl;

  //Initializing the calibration parameters to be read in from the file
  for (int iz=0; iz<NZPS; iz++)
    {
      for (int ih=1; ih<NHAR; ih++)
        {
          for (int id=0; id<NDET; id++)
            {
              for (int ib=0; ib<2; ib++)
                {
                  mean[ic][iz][ih][id][ib]=0.0;
                  widt[ic][iz][ih][id][ib]=1.0;

                  for (int io=0; io<NORD; io++)
                    {
                      four[ic][iz][ih][id][ib][io]=0.0;
                    } // orders
                } // x and y
            } // detector
        } // harmonics
    } // z_vertex bins

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
  if(rp_recal_pass>=2)
    {

      cout << "reading calibration file : " << calibfile << endl;
      float f0,f1,f2,f3;//f4,f5,f6,f7;
      ifstream ifs;
      ifs.open(calibfile);
      for (int iz=0; iz<NZPS; iz++)
        {
          for (int ih=1; ih<NHAR; ih++)
            {
              for (int id=0; id<NDET; id++)
                {
                  ifs >> f0 >> f1 >> f2 >> f3;
                  if ( f1 <= 0.0 ) f1=1.0;
                  if ( f3 <= 0.0 ) f3=1.0;
                  mean[ic][iz][ih][id][0]=f0;
                  widt[ic][iz][ih][id][0]=f1;
                  mean[ic][iz][ih][id][1]=f2;
                  widt[ic][iz][ih][id][1]=f3;
                  if ( id==2 && ih == 1 && DIAG ) cout<<f0<<" "<<f1<<" "<<f2<<" "<<f3<<endl;//bbc psi 2 parameters
                  if ( id==2 && ih == 1 && DIAG ) cout << "---" << endl;
                  if ( id==3 && ih == 1 && DIAG ) cout<<f0<<" "<<f1<<" "<<f2<<" "<<f3<<endl;//bbc psi 2 parameters
                  for (int ib=0; ib<2; ib++)
                    {
                      for (int io=0; io<NORD; io++)
                        {
                          ifs>>four[ic][iz][ih][id][ib][io];
                        } // orders
                    } // x and y
                } // detectors
            } // harmonics
        } // z_vertex bins
      ifs.close();
    } // check on second or third pass


  //------------------------------------------------------------//
  //  Finished Reading in flattening calibration parameters     //
  //------------------------------------------------------------//



  //------------------------------------------------------------//
  //                  Initializing histograms                   //
  //------------------------------------------------------------//

  cout << "Initializing more histograms" << endl;

  TProfile* bbcs_v2_incl_nodetree = new TProfile("bbcs_v2_incl_nodetree","bbcs_v2_incl_nodetree",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* bbcs_v2_east_nodetree = new TProfile("bbcs_v2_east_nodetree","bbcs_v2_east_nodetree",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* bbcs_v2_west_nodetree = new TProfile("bbcs_v2_west_nodetree","bbcs_v2_west_nodetree",15, 0.0, 3.0, -1.1, 1.1);

  TProfile* fvtxs_v2_incl_nodetree = new TProfile("fvtxs_v2_incl_nodetree","fvtxs_v2_incl_nodetree",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* fvtxs_v2_east_nodetree = new TProfile("fvtxs_v2_east_nodetree","fvtxs_v2_east_nodetree",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* fvtxs_v2_west_nodetree = new TProfile("fvtxs_v2_west_nodetree","fvtxs_v2_west_nodetree",15, 0.0, 3.0, -1.1, 1.1);

  TProfile* bbcs_v2_west_docalib = new TProfile(Form("bbcs_v2_west_docalib"),Form("bbcs_v2_west_docalib"),15, 0.0, 3.0, -1.1, 1.1);
  TProfile* bbcs_v2_east_docalib = new TProfile(Form("bbcs_v2_east_docalib"),Form("bbcs_v2_east_docalib"),15, 0.0, 3.0, -1.1, 1.1);
  TProfile* bbcs_v2_both_docalib = new TProfile(Form("bbcs_v2_both_docalib"),Form("bbcs_v2_both_docalib"),15, 0.0, 3.0, -1.1, 1.1);

  TProfile* fvtxs_v2_west_docalib = new TProfile(Form("fvtxs_v2_west_docalib"),Form("fvtxs_v2_west_docalib"),15, 0.0, 3.0, -1.1, 1.1);
  TProfile* fvtxs_v2_east_docalib = new TProfile(Form("fvtxs_v2_east_docalib"),Form("fvtxs_v2_east_docalib"),15, 0.0, 3.0, -1.1, 1.1);
  TProfile* fvtxs_v2_both_docalib = new TProfile(Form("fvtxs_v2_both_docalib"),Form("fvtxs_v2_both_docalib"),15, 0.0, 3.0, -1.1, 1.1);

  TProfile* fvtxs0_v2_west_docalib = new TProfile(Form("fvtxs0_v2_west_docalib"),Form("fvtxs0_v2_west_docalib"),15, 0.0, 3.0, -1.1, 1.1);
  TProfile* fvtxs0_v2_east_docalib = new TProfile(Form("fvtxs0_v2_east_docalib"),Form("fvtxs0_v2_east_docalib"),15, 0.0, 3.0, -1.1, 1.1);
  TProfile* fvtxs0_v2_both_docalib = new TProfile(Form("fvtxs0_v2_both_docalib"),Form("fvtxs0_v2_both_docalib"),15, 0.0, 3.0, -1.1, 1.1);

  TProfile* fvtxs1_v2_west_docalib = new TProfile(Form("fvtxs1_v2_west_docalib"),Form("fvtxs1_v2_west_docalib"),15, 0.0, 3.0, -1.1, 1.1);
  TProfile* fvtxs1_v2_east_docalib = new TProfile(Form("fvtxs1_v2_east_docalib"),Form("fvtxs1_v2_east_docalib"),15, 0.0, 3.0, -1.1, 1.1);
  TProfile* fvtxs1_v2_both_docalib = new TProfile(Form("fvtxs1_v2_both_docalib"),Form("fvtxs1_v2_both_docalib"),15, 0.0, 3.0, -1.1, 1.1);

  TProfile* fvtxs2_v2_west_docalib = new TProfile(Form("fvtxs2_v2_west_docalib"),Form("fvtxs2_v2_west_docalib"),15, 0.0, 3.0, -1.1, 1.1);
  TProfile* fvtxs2_v2_east_docalib = new TProfile(Form("fvtxs2_v2_east_docalib"),Form("fvtxs2_v2_east_docalib"),15, 0.0, 3.0, -1.1, 1.1);
  TProfile* fvtxs2_v2_both_docalib = new TProfile(Form("fvtxs2_v2_both_docalib"),Form("fvtxs2_v2_both_docalib"),15, 0.0, 3.0, -1.1, 1.1);

  TProfile* fvtxs3_v2_west_docalib = new TProfile(Form("fvtxs3_v2_west_docalib"),Form("fvtxs3_v2_west_docalib"),15, 0.0, 3.0, -1.1, 1.1);
  TProfile* fvtxs3_v2_east_docalib = new TProfile(Form("fvtxs3_v2_east_docalib"),Form("fvtxs3_v2_east_docalib"),15, 0.0, 3.0, -1.1, 1.1);
  TProfile* fvtxs3_v2_both_docalib = new TProfile(Form("fvtxs3_v2_both_docalib"),Form("fvtxs3_v2_both_docalib"),15, 0.0, 3.0, -1.1, 1.1);

  // ---

  TProfile* bbcs_v3_west_docalib = new TProfile(Form("bbcs_v3_west_docalib"),Form("bbcs_v3_west_docalib"),15, 0.0, 3.0, -1.1, 1.1);
  TProfile* bbcs_v3_east_docalib = new TProfile(Form("bbcs_v3_east_docalib"),Form("bbcs_v3_east_docalib"),15, 0.0, 3.0, -1.1, 1.1);
  TProfile* bbcs_v3_both_docalib = new TProfile(Form("bbcs_v3_both_docalib"),Form("bbcs_v3_both_docalib"),15, 0.0, 3.0, -1.1, 1.1);

  TProfile* fvtxs_v3_west_docalib = new TProfile(Form("fvtxs_v3_west_docalib"),Form("fvtxs_v3_west_docalib"),15, 0.0, 3.0, -1.1, 1.1);
  TProfile* fvtxs_v3_east_docalib = new TProfile(Form("fvtxs_v3_east_docalib"),Form("fvtxs_v3_east_docalib"),15, 0.0, 3.0, -1.1, 1.1);
  TProfile* fvtxs_v3_both_docalib = new TProfile(Form("fvtxs_v3_both_docalib"),Form("fvtxs_v3_both_docalib"),15, 0.0, 3.0, -1.1, 1.1);

  // ---

  // ---

  TProfile* bbcs_v4_4Psi2_west_docalib = new TProfile(Form("bbcs_v4_4Psi2_west_docalib"),Form("bbcs_v4_4Psi2_west_docalib"),15, 0.0, 3.0, -1.1, 1.1);
  TProfile* bbcs_v4_4Psi2_east_docalib = new TProfile(Form("bbcs_v4_4Psi2_east_docalib"),Form("bbcs_v4_4Psi2_east_docalib"),15, 0.0, 3.0, -1.1, 1.1);
  TProfile* bbcs_v4_4Psi2_both_docalib = new TProfile(Form("bbcs_v4_4Psi2_both_docalib"),Form("bbcs_v4_4Psi2_both_docalib"),15, 0.0, 3.0, -1.1, 1.1);

  TProfile* fvtxs_v4_4Psi2_west_docalib = new TProfile(Form("fvtxs_v4_4Psi2_west_docalib"),Form("fvtxs_v4_4Psi2_west_docalib"),15, 0.0, 3.0, -1.1, 1.1);
  TProfile* fvtxs_v4_4Psi2_east_docalib = new TProfile(Form("fvtxs_v4_4Psi2_east_docalib"),Form("fvtxs_v4_4Psi2_east_docalib"),15, 0.0, 3.0, -1.1, 1.1);
  TProfile* fvtxs_v4_4Psi2_both_docalib = new TProfile(Form("fvtxs_v4_4Psi2_both_docalib"),Form("fvtxs_v4_4Psi2_both_docalib"),15, 0.0, 3.0, -1.1, 1.1);

  // ---

  TProfile* fvtxn_v2_west_docalib = new TProfile(Form("fvtxn_v2_west_docalib"),Form("fvtxn_v2_west_docalib"),15, 0.0, 3.0, -1.1, 1.1);
  TProfile* fvtxn_v2_east_docalib = new TProfile(Form("fvtxn_v2_east_docalib"),Form("fvtxn_v2_east_docalib"),15, 0.0, 3.0, -1.1, 1.1);
  TProfile* fvtxn_v2_both_docalib = new TProfile(Form("fvtxn_v2_both_docalib"),Form("fvtxn_v2_both_docalib"),15, 0.0, 3.0, -1.1, 1.1);

  TProfile* fvtxn_v3_west_docalib = new TProfile(Form("fvtxn_v3_west_docalib"),Form("fvtxn_v3_west_docalib"),15, 0.0, 3.0, -1.1, 1.1);
  TProfile* fvtxn_v3_east_docalib = new TProfile(Form("fvtxn_v3_east_docalib"),Form("fvtxn_v3_east_docalib"),15, 0.0, 3.0, -1.1, 1.1);
  TProfile* fvtxn_v3_both_docalib = new TProfile(Form("fvtxn_v3_both_docalib"),Form("fvtxn_v3_both_docalib"),15, 0.0, 3.0, -1.1, 1.1);

  TProfile* fvtxn_v4_4Psi2_west_docalib = new TProfile(Form("fvtxn_v4_4Psi2_west_docalib"),Form("fvtxn_v4_4Psi2_west_docalib"),15, 0.0, 3.0, -1.1, 1.1);
  TProfile* fvtxn_v4_4Psi2_east_docalib = new TProfile(Form("fvtxn_v4_4Psi2_east_docalib"),Form("fvtxn_v4_4Psi2_east_docalib"),15, 0.0, 3.0, -1.1, 1.1);
  TProfile* fvtxn_v4_4Psi2_both_docalib = new TProfile(Form("fvtxn_v4_4Psi2_both_docalib"),Form("fvtxn_v4_4Psi2_both_docalib"),15, 0.0, 3.0, -1.1, 1.1);

  // ---------------------------------------------------------------------------------------------------------
  // --- COME BACK HERE FOR 2PC HISTOGRAMS

  TProfile* bbcs_d22_west = new TProfile(Form("bbcs_d22_west"),Form("bbcs_d22_west"),15, 0.0, 3.0, -1.1, 1.1);
  TProfile* bbcs_d22_east = new TProfile(Form("bbcs_d22_east"),Form("bbcs_d22_east"),15, 0.0, 3.0, -1.1, 1.1);
  TProfile* bbcs_d22_both = new TProfile(Form("bbcs_d22_both"),Form("bbcs_d22_both"),15, 0.0, 3.0, -1.1, 1.1);

  TProfile* fvtxs_d22_west = new TProfile(Form("fvtxs_d22_west"),Form("fvtxs_d22_west"),15, 0.0, 3.0, -1.1, 1.1);
  TProfile* fvtxs_d22_east = new TProfile(Form("fvtxs_d22_east"),Form("fvtxs_d22_east"),15, 0.0, 3.0, -1.1, 1.1);
  TProfile* fvtxs_d22_both = new TProfile(Form("fvtxs_d22_both"),Form("fvtxs_d22_both"),15, 0.0, 3.0, -1.1, 1.1);

  TProfile* bbcs_d32_west = new TProfile(Form("bbcs_d32_west"),Form("bbcs_d32_west"),15, 0.0, 3.0, -1.1, 1.1);
  TProfile* bbcs_d32_east = new TProfile(Form("bbcs_d32_east"),Form("bbcs_d32_east"),15, 0.0, 3.0, -1.1, 1.1);
  TProfile* bbcs_d32_both = new TProfile(Form("bbcs_d32_both"),Form("bbcs_d32_both"),15, 0.0, 3.0, -1.1, 1.1);

  TProfile* fvtxs_d32_west = new TProfile(Form("fvtxs_d32_west"),Form("fvtxs_d32_west"),15, 0.0, 3.0, -1.1, 1.1);
  TProfile* fvtxs_d32_east = new TProfile(Form("fvtxs_d32_east"),Form("fvtxs_d32_east"),15, 0.0, 3.0, -1.1, 1.1);
  TProfile* fvtxs_d32_both = new TProfile(Form("fvtxs_d32_both"),Form("fvtxs_d32_both"),15, 0.0, 3.0, -1.1, 1.1);

  TProfile* bbcs_sin22_west = new TProfile(Form("bbcs_sin22_west"),Form("bbcs_sin22_west"),15, 0.0, 3.0, -1.1, 1.1);
  TProfile* bbcs_sin22_east = new TProfile(Form("bbcs_sin22_east"),Form("bbcs_sin22_east"),15, 0.0, 3.0, -1.1, 1.1);
  TProfile* bbcs_sin22_both = new TProfile(Form("bbcs_sin22_both"),Form("bbcs_sin22_both"),15, 0.0, 3.0, -1.1, 1.1);

  TProfile* fvtxs_sin22_west = new TProfile(Form("fvtxs_sin22_west"),Form("fvtxs_sin22_west"),15, 0.0, 3.0, -1.1, 1.1);
  TProfile* fvtxs_sin22_east = new TProfile(Form("fvtxs_sin22_east"),Form("fvtxs_sin22_east"),15, 0.0, 3.0, -1.1, 1.1);
  TProfile* fvtxs_sin22_both = new TProfile(Form("fvtxs_sin22_both"),Form("fvtxs_sin22_both"),15, 0.0, 3.0, -1.1, 1.1);

  TProfile* bbcs_sin32_west = new TProfile(Form("bbcs_sin32_west"),Form("bbcs_sin32_west"),15, 0.0, 3.0, -1.1, 1.1);
  TProfile* bbcs_sin32_east = new TProfile(Form("bbcs_sin32_east"),Form("bbcs_sin32_east"),15, 0.0, 3.0, -1.1, 1.1);
  TProfile* bbcs_sin32_both = new TProfile(Form("bbcs_sin32_both"),Form("bbcs_sin32_both"),15, 0.0, 3.0, -1.1, 1.1);

  TProfile* fvtxs_sin32_west = new TProfile(Form("fvtxs_sin32_west"),Form("fvtxs_sin32_west"),15, 0.0, 3.0, -1.1, 1.1);
  TProfile* fvtxs_sin32_east = new TProfile(Form("fvtxs_sin32_east"),Form("fvtxs_sin32_east"),15, 0.0, 3.0, -1.1, 1.1);
  TProfile* fvtxs_sin32_both = new TProfile(Form("fvtxs_sin32_both"),Form("fvtxs_sind32_both"),15, 0.0, 3.0, -1.1, 1.1);

  TProfile* bbcs_cos22_west = new TProfile(Form("bbcs_cos22_west"),Form("bbcs_cos22_west"),15, 0.0, 3.0, -1.1, 1.1);
  TProfile* bbcs_cos22_east = new TProfile(Form("bbcs_cos22_east"),Form("bbcs_cos22_east"),15, 0.0, 3.0, -1.1, 1.1);
  TProfile* bbcs_cos22_both = new TProfile(Form("bbcs_cos22_both"),Form("bbcs_cos22_both"),15, 0.0, 3.0, -1.1, 1.1);

  TProfile* fvtxs_cos22_west = new TProfile(Form("fvtxs_cos22_west"),Form("fvtxs_cos22_west"),15, 0.0, 3.0, -1.1, 1.1);
  TProfile* fvtxs_cos22_east = new TProfile(Form("fvtxs_cos22_east"),Form("fvtxs_cos22_east"),15, 0.0, 3.0, -1.1, 1.1);
  TProfile* fvtxs_cos22_both = new TProfile(Form("fvtxs_cos22_both"),Form("fvtxs_cos22_both"),15, 0.0, 3.0, -1.1, 1.1);

  TProfile* bbcs_cos32_west = new TProfile(Form("bbcs_cos32_west"),Form("bbcs_cos32_west"),15, 0.0, 3.0, -1.1, 1.1);
  TProfile* bbcs_cos32_east = new TProfile(Form("bbcs_cos32_east"),Form("bbcs_cos32_east"),15, 0.0, 3.0, -1.1, 1.1);
  TProfile* bbcs_cos32_both = new TProfile(Form("bbcs_cos32_both"),Form("bbcs_cos32_both"),15, 0.0, 3.0, -1.1, 1.1);

  TProfile* fvtxs_cos32_west = new TProfile(Form("fvtxs_cos32_west"),Form("fvtxs_cos32_west"),15, 0.0, 3.0, -1.1, 1.1);
  TProfile* fvtxs_cos32_east = new TProfile(Form("fvtxs_cos32_east"),Form("fvtxs_cos32_east"),15, 0.0, 3.0, -1.1, 1.1);
  TProfile* fvtxs_cos32_both = new TProfile(Form("fvtxs_cos32_both"),Form("fvtxs_cos32_both"),15, 0.0, 3.0, -1.1, 1.1);

  TProfile* bbcs_c22 = new TProfile(Form("bbcs_c22"),Form("bbcs_c22"),1, -0.5, 0.5, -1.1, 1.1);
  TProfile* fvtxs_c22 = new TProfile(Form("fvtxs_c22"),Form("fvtxs_c22"),1, -0.5, 0.5, -1.1, 1.1);
  TProfile* bbcs_c32 = new TProfile(Form("bbcs_c32"),Form("bbcs_c32"),1, -0.5, 0.5, -1.1, 1.1);
  TProfile* fvtxs_c32 = new TProfile(Form("fvtxs_c32"),Form("fvtxs_c32"),1, -0.5, 0.5, -1.1, 1.1);

  TProfile* bbcs_sin22 = new TProfile(Form("bbcs_sin22"),Form("bbcs_sin22"),1, -0.5, 0.5, -1.1, 1.1);
  TProfile* fvtxs_sin22 = new TProfile(Form("fvtxs_sin22"),Form("fvtxs_sin22"),1, -0.5, 0.5, -1.1, 1.1);
  TProfile* bbcs_sin32 = new TProfile(Form("bbcs_sin32"),Form("bbcs_sin32"),1, -0.5, 0.5, -1.1, 1.1);
  TProfile* fvtxs_sin32 = new TProfile(Form("fvtxs_sin32"),Form("fvtxs_sin32"),1, -0.5, 0.5, -1.1, 1.1);

  TProfile* bbcs_cos22 = new TProfile(Form("bbcs_cos22"),Form("bbcs_cos22"),1, -0.5, 0.5, -1.1, 1.1);
  TProfile* fvtxs_cos22 = new TProfile(Form("fvtxs_cos22"),Form("fvtxs_cos22"),1, -0.5, 0.5, -1.1, 1.1);
  TProfile* bbcs_cos32 = new TProfile(Form("bbcs_cos32"),Form("bbcs_cos32"),1, -0.5, 0.5, -1.1, 1.1);
  TProfile* fvtxs_cos32 = new TProfile(Form("fvtxs_cos32"),Form("fvtxs_cos32"),1, -0.5, 0.5, -1.1, 1.1);

  // ---------------------------------------------------------------------------------------------------------

  TProfile* os_bbcs_d22_west = new TProfile(Form("os_bbcs_d22_west"),Form("os_bbcs_d22_west"),15, 0.0, 3.0, -1.1, 1.1);
  TProfile* os_bbcs_d22_east = new TProfile(Form("os_bbcs_d22_east"),Form("os_bbcs_d22_east"),15, 0.0, 3.0, -1.1, 1.1);
  TProfile* os_bbcs_d22_both = new TProfile(Form("os_bbcs_d22_both"),Form("os_bbcs_d22_both"),15, 0.0, 3.0, -1.1, 1.1);

  TProfile* os_fvtxs_d22_west = new TProfile(Form("os_fvtxs_d22_west"),Form("os_fvtxs_d22_west"),15, 0.0, 3.0, -1.1, 1.1);
  TProfile* os_fvtxs_d22_east = new TProfile(Form("os_fvtxs_d22_east"),Form("os_fvtxs_d22_east"),15, 0.0, 3.0, -1.1, 1.1);
  TProfile* os_fvtxs_d22_both = new TProfile(Form("os_fvtxs_d22_both"),Form("os_fvtxs_d22_both"),15, 0.0, 3.0, -1.1, 1.1);

  TProfile* os_fvtxn_d22_west = new TProfile(Form("os_fvtxn_d22_west"),Form("os_fvtxn_d22_west"),15, 0.0, 3.0, -1.1, 1.1);
  TProfile* os_fvtxn_d22_east = new TProfile(Form("os_fvtxn_d22_east"),Form("os_fvtxn_d22_east"),15, 0.0, 3.0, -1.1, 1.1);
  TProfile* os_fvtxn_d22_both = new TProfile(Form("os_fvtxn_d22_both"),Form("os_fvtxn_d22_both"),15, 0.0, 3.0, -1.1, 1.1);

  TProfile* os_bbcs_d32_west = new TProfile(Form("os_bbcs_d32_west"),Form("os_bbcs_d32_west"),15, 0.0, 3.0, -1.1, 1.1);
  TProfile* os_bbcs_d32_east = new TProfile(Form("os_bbcs_d32_east"),Form("os_bbcs_d32_east"),15, 0.0, 3.0, -1.1, 1.1);
  TProfile* os_bbcs_d32_both = new TProfile(Form("os_bbcs_d32_both"),Form("os_bbcs_d32_both"),15, 0.0, 3.0, -1.1, 1.1);

  TProfile* os_fvtxs_d32_west = new TProfile(Form("os_fvtxs_d32_west"),Form("os_fvtxs_d32_west"),15, 0.0, 3.0, -1.1, 1.1);
  TProfile* os_fvtxs_d32_east = new TProfile(Form("os_fvtxs_d32_east"),Form("os_fvtxs_d32_east"),15, 0.0, 3.0, -1.1, 1.1);
  TProfile* os_fvtxs_d32_both = new TProfile(Form("os_fvtxs_d32_both"),Form("os_fvtxs_d32_both"),15, 0.0, 3.0, -1.1, 1.1);

  TProfile* os_fvtxn_d32_west = new TProfile(Form("os_fvtxn_d32_west"),Form("os_fvtxn_d32_west"),15, 0.0, 3.0, -1.1, 1.1);
  TProfile* os_fvtxn_d32_east = new TProfile(Form("os_fvtxn_d32_east"),Form("os_fvtxn_d32_east"),15, 0.0, 3.0, -1.1, 1.1);
  TProfile* os_fvtxn_d32_both = new TProfile(Form("os_fvtxn_d32_both"),Form("os_fvtxn_d32_both"),15, 0.0, 3.0, -1.1, 1.1);

  TProfile* os_bbcs_sin22_west = new TProfile(Form("os_bbcs_sin22_west"),Form("os_bbcs_sin22_west"),15, 0.0, 3.0, -1.1, 1.1);
  TProfile* os_bbcs_sin22_east = new TProfile(Form("os_bbcs_sin22_east"),Form("os_bbcs_sin22_east"),15, 0.0, 3.0, -1.1, 1.1);
  TProfile* os_bbcs_sin22_both = new TProfile(Form("os_bbcs_sin22_both"),Form("os_bbcs_sin22_both"),15, 0.0, 3.0, -1.1, 1.1);

  TProfile* os_fvtxs_sin22_west = new TProfile(Form("os_fvtxs_sin22_west"),Form("os_fvtxs_sin22_west"),15, 0.0, 3.0, -1.1, 1.1);
  TProfile* os_fvtxs_sin22_east = new TProfile(Form("os_fvtxs_sin22_east"),Form("os_fvtxs_sin22_east"),15, 0.0, 3.0, -1.1, 1.1);
  TProfile* os_fvtxs_sin22_both = new TProfile(Form("os_fvtxs_sin22_both"),Form("os_fvtxs_sin22_both"),15, 0.0, 3.0, -1.1, 1.1);

  TProfile* os_fvtxn_sin22_west = new TProfile(Form("os_fvtxn_sin22_west"),Form("os_fvtxn_sin22_west"),15, 0.0, 3.0, -1.1, 1.1);
  TProfile* os_fvtxn_sin22_east = new TProfile(Form("os_fvtxn_sin22_east"),Form("os_fvtxn_sin22_east"),15, 0.0, 3.0, -1.1, 1.1);
  TProfile* os_fvtxn_sin22_both = new TProfile(Form("os_fvtxn_sin22_both"),Form("os_fvtxn_sin22_both"),15, 0.0, 3.0, -1.1, 1.1);

  TProfile* os_bbcs_sin32_west = new TProfile(Form("os_bbcs_sin32_west"),Form("os_bbcs_sin32_west"),15, 0.0, 3.0, -1.1, 1.1);
  TProfile* os_bbcs_sin32_east = new TProfile(Form("os_bbcs_sin32_east"),Form("os_bbcs_sin32_east"),15, 0.0, 3.0, -1.1, 1.1);
  TProfile* os_bbcs_sin32_both = new TProfile(Form("os_bbcs_sin32_both"),Form("os_bbcs_sin32_both"),15, 0.0, 3.0, -1.1, 1.1);

  TProfile* os_fvtxs_sin32_west = new TProfile(Form("os_fvtxs_sin32_west"),Form("os_fvtxs_sin32_west"),15, 0.0, 3.0, -1.1, 1.1);
  TProfile* os_fvtxs_sin32_east = new TProfile(Form("os_fvtxs_sin32_east"),Form("os_fvtxs_sin32_east"),15, 0.0, 3.0, -1.1, 1.1);
  TProfile* os_fvtxs_sin32_both = new TProfile(Form("os_fvtxs_sin32_both"),Form("os_fvtxs_sind32_both"),15, 0.0, 3.0, -1.1, 1.1);

  TProfile* os_fvtxn_sin32_west = new TProfile(Form("os_fvtxn_sin32_west"),Form("os_fvtxn_sin32_west"),15, 0.0, 3.0, -1.1, 1.1);
  TProfile* os_fvtxn_sin32_east = new TProfile(Form("os_fvtxn_sin32_east"),Form("os_fvtxn_sin32_east"),15, 0.0, 3.0, -1.1, 1.1);
  TProfile* os_fvtxn_sin32_both = new TProfile(Form("os_fvtxn_sin32_both"),Form("os_fvtxn_sind32_both"),15, 0.0, 3.0, -1.1, 1.1);

  TProfile* os_bbcs_cos22_west = new TProfile(Form("os_bbcs_cos22_west"),Form("os_bbcs_cos22_west"),15, 0.0, 3.0, -1.1, 1.1);
  TProfile* os_bbcs_cos22_east = new TProfile(Form("os_bbcs_cos22_east"),Form("os_bbcs_cos22_east"),15, 0.0, 3.0, -1.1, 1.1);
  TProfile* os_bbcs_cos22_both = new TProfile(Form("os_bbcs_cos22_both"),Form("os_bbcs_cos22_both"),15, 0.0, 3.0, -1.1, 1.1);

  TProfile* os_fvtxs_cos22_west = new TProfile(Form("os_fvtxs_cos22_west"),Form("os_fvtxs_cos22_west"),15, 0.0, 3.0, -1.1, 1.1);
  TProfile* os_fvtxs_cos22_east = new TProfile(Form("os_fvtxs_cos22_east"),Form("os_fvtxs_cos22_east"),15, 0.0, 3.0, -1.1, 1.1);
  TProfile* os_fvtxs_cos22_both = new TProfile(Form("os_fvtxs_cos22_both"),Form("os_fvtxs_cos22_both"),15, 0.0, 3.0, -1.1, 1.1);

  TProfile* os_fvtxn_cos22_west = new TProfile(Form("os_fvtxn_cos22_west"),Form("os_fvtxn_cos22_west"),15, 0.0, 3.0, -1.1, 1.1);
  TProfile* os_fvtxn_cos22_east = new TProfile(Form("os_fvtxn_cos22_east"),Form("os_fvtxn_cos22_east"),15, 0.0, 3.0, -1.1, 1.1);
  TProfile* os_fvtxn_cos22_both = new TProfile(Form("os_fvtxn_cos22_both"),Form("os_fvtxn_cos22_both"),15, 0.0, 3.0, -1.1, 1.1);

  TProfile* os_bbcs_cos32_west = new TProfile(Form("os_bbcs_cos32_west"),Form("os_bbcs_cos32_west"),15, 0.0, 3.0, -1.1, 1.1);
  TProfile* os_bbcs_cos32_east = new TProfile(Form("os_bbcs_cos32_east"),Form("os_bbcs_cos32_east"),15, 0.0, 3.0, -1.1, 1.1);
  TProfile* os_bbcs_cos32_both = new TProfile(Form("os_bbcs_cos32_both"),Form("os_bbcs_cos32_both"),15, 0.0, 3.0, -1.1, 1.1);

  TProfile* os_fvtxs_cos32_west = new TProfile(Form("os_fvtxs_cos32_west"),Form("os_fvtxs_cos32_west"),15, 0.0, 3.0, -1.1, 1.1);
  TProfile* os_fvtxs_cos32_east = new TProfile(Form("os_fvtxs_cos32_east"),Form("os_fvtxs_cos32_east"),15, 0.0, 3.0, -1.1, 1.1);
  TProfile* os_fvtxs_cos32_both = new TProfile(Form("os_fvtxs_cos32_both"),Form("os_fvtxs_cos32_both"),15, 0.0, 3.0, -1.1, 1.1);

  TProfile* os_fvtxn_cos32_west = new TProfile(Form("os_fvtxn_cos32_west"),Form("os_fvtxn_cos32_west"),15, 0.0, 3.0, -1.1, 1.1);
  TProfile* os_fvtxn_cos32_east = new TProfile(Form("os_fvtxn_cos32_east"),Form("os_fvtxn_cos32_east"),15, 0.0, 3.0, -1.1, 1.1);
  TProfile* os_fvtxn_cos32_both = new TProfile(Form("os_fvtxn_cos32_both"),Form("os_fvtxn_cos32_both"),15, 0.0, 3.0, -1.1, 1.1);

  TProfile* os_bbcs_c22 = new TProfile(Form("os_bbcs_c22"),Form("os_bbcs_c22"),1, -0.5, 0.5, -1.1, 1.1);
  TProfile* os_fvtxs_c22 = new TProfile(Form("os_fvtxs_c22"),Form("os_fvtxs_c22"),1, -0.5, 0.5, -1.1, 1.1);
  TProfile* os_fvtxn_c22 = new TProfile(Form("os_fvtxn_c22"),Form("os_fvtxn_c22"),1, -0.5, 0.5, -1.1, 1.1);
  TProfile* os_bbcs_c32 = new TProfile(Form("os_bbcs_c32"),Form("os_bbcs_c32"),1, -0.5, 0.5, -1.1, 1.1);
  TProfile* os_fvtxs_c32 = new TProfile(Form("os_fvtxs_c32"),Form("os_fvtxs_c32"),1, -0.5, 0.5, -1.1, 1.1);
  TProfile* os_fvtxn_c32 = new TProfile(Form("os_fvtxn_c32"),Form("os_fvtxn_c32"),1, -0.5, 0.5, -1.1, 1.1);

  TProfile* os_bbcs_sin22 = new TProfile(Form("os_bbcs_sin22"),Form("os_bbcs_sin22"),1, -0.5, 0.5, -1.1, 1.1);
  TProfile* os_fvtxs_sin22 = new TProfile(Form("os_fvtxs_sin22"),Form("os_fvtxs_sin22"),1, -0.5, 0.5, -1.1, 1.1);
  TProfile* os_fvtxn_sin22 = new TProfile(Form("os_fvtxn_sin22"),Form("os_fvtxn_sin22"),1, -0.5, 0.5, -1.1, 1.1);
  TProfile* os_bbcs_sin32 = new TProfile(Form("os_bbcs_sin32"),Form("os_bbcs_sin32"),1, -0.5, 0.5, -1.1, 1.1);
  TProfile* os_fvtxs_sin32 = new TProfile(Form("os_fvtxs_sin32"),Form("os_fvtxs_sin32"),1, -0.5, 0.5, -1.1, 1.1);
  TProfile* os_fvtxn_sin32 = new TProfile(Form("os_fvtxn_sin32"),Form("os_fvtxn_sin32"),1, -0.5, 0.5, -1.1, 1.1);

  TProfile* os_bbcs_cos22 = new TProfile(Form("os_bbcs_cos22"),Form("os_bbcs_cos22"),1, -0.5, 0.5, -1.1, 1.1);
  TProfile* os_fvtxs_cos22 = new TProfile(Form("os_fvtxs_cos22"),Form("os_fvtxs_cos22"),1, -0.5, 0.5, -1.1, 1.1);
  TProfile* os_fvtxn_cos22 = new TProfile(Form("os_fvtxn_cos22"),Form("os_fvtxn_cos22"),1, -0.5, 0.5, -1.1, 1.1);
  TProfile* os_bbcs_cos32 = new TProfile(Form("os_bbcs_cos32"),Form("os_bbcs_cos32"),1, -0.5, 0.5, -1.1, 1.1);
  TProfile* os_fvtxs_cos32 = new TProfile(Form("os_fvtxs_cos32"),Form("os_fvtxs_cos32"),1, -0.5, 0.5, -1.1, 1.1);
  TProfile* os_fvtxn_cos32 = new TProfile(Form("os_fvtxn_cos32"),Form("os_fvtxn_cos32"),1, -0.5, 0.5, -1.1, 1.1);

  TProfile* os_bbcsfvtxs_c22 = new TProfile(Form("os_bbcsfvtxs_c22"),Form("os_bbcsfvtxs_c22"),1, -0.5, 0.5, -1.1, 1.1);
  TProfile* os_bbcsfvtxs_c32 = new TProfile(Form("os_bbcsfvtxs_c32"),Form("os_bbcsfvtxs_c32"),1, -0.5, 0.5, -1.1, 1.1);
  TProfile* os_bbcsfvtxn_c22 = new TProfile(Form("os_bbcsfvtxn_c22"),Form("os_bbcsfvtxn_c22"),1, -0.5, 0.5, -1.1, 1.1);
  TProfile* os_bbcsfvtxn_c32 = new TProfile(Form("os_bbcsfvtxn_c32"),Form("os_bbcsfvtxn_c32"),1, -0.5, 0.5, -1.1, 1.1);
  TProfile* os_fvtxsfvtxn_c22 = new TProfile(Form("os_fvtxsfvtxn_c22"),Form("os_fvtxsfvtxn_c22"),1, -0.5, 0.5, -1.1, 1.1);
  TProfile* os_fvtxsfvtxn_c32 = new TProfile(Form("os_fvtxsfvtxn_c32"),Form("os_fvtxsfvtxn_c32"),1, -0.5, 0.5, -1.1, 1.1);

  TH1D* os_bbcs_1dPsi2 = new TH1D(Form("os_bbcs_1dPsi2"),Form("os_bbcs_1dPsi2"),220,-4.1,4.1); // weird binning from elsewhere, leave the same to match
  TH1D* os_bbcs_1dPsi3 = new TH1D(Form("os_bbcs_1dPsi3"),Form("os_bbcs_1dPsi3"),220,-4.1,4.1); // weird binning from elsewhere, leave the same to match
  TH1D* os_fvtxs_1dPsi2 = new TH1D(Form("os_fvtxs_1dPsi2"),Form("os_fvtxs_1dPsi2"),200,-4.1,4.1); // weird binning from elsewhere, leave the same to match
  TH1D* os_fvtxs_1dPsi3 = new TH1D(Form("os_fvtxs_1dPsi3"),Form("os_fvtxs_1dPsi3"),200,-4.1,4.1); // weird binning from elsewhere, leave the same to match
  TH1D* os_fvtxn_1dPsi2 = new TH1D(Form("os_fvtxn_1dPsi2"),Form("os_fvtxn_1dPsi2"),200,-4.1,4.1); // weird binning from elsewhere, leave the same to match
  TH1D* os_fvtxn_1dPsi3 = new TH1D(Form("os_fvtxn_1dPsi3"),Form("os_fvtxn_1dPsi3"),200,-4.1,4.1); // weird binning from elsewhere, leave the same to match

  TProfile* os_bbcs_v2_west = new TProfile(Form("os_bbcs_v2_west"),Form("os_bbcs_v2_west"),15, 0.0, 3.0, -1.1, 1.1);
  TProfile* os_bbcs_v2_east = new TProfile(Form("os_bbcs_v2_east"),Form("os_bbcs_v2_east"),15, 0.0, 3.0, -1.1, 1.1);
  TProfile* os_bbcs_v2_both = new TProfile(Form("os_bbcs_v2_both"),Form("os_bbcs_v2_both"),15, 0.0, 3.0, -1.1, 1.1);

  TProfile* os_fvtxs_v2_west = new TProfile(Form("os_fvtxs_v2_west"),Form("os_fvtxs_v2_west"),15, 0.0, 3.0, -1.1, 1.1);
  TProfile* os_fvtxs_v2_east = new TProfile(Form("os_fvtxs_v2_east"),Form("os_fvtxs_v2_east"),15, 0.0, 3.0, -1.1, 1.1);
  TProfile* os_fvtxs_v2_both = new TProfile(Form("os_fvtxs_v2_both"),Form("os_fvtxs_v2_both"),15, 0.0, 3.0, -1.1, 1.1);

  TProfile* os_fvtxn_v2_west = new TProfile(Form("os_fvtxn_v2_west"),Form("os_fvtxn_v2_west"),15, 0.0, 3.0, -1.1, 1.1);
  TProfile* os_fvtxn_v2_east = new TProfile(Form("os_fvtxn_v2_east"),Form("os_fvtxn_v2_east"),15, 0.0, 3.0, -1.1, 1.1);
  TProfile* os_fvtxn_v2_both = new TProfile(Form("os_fvtxn_v2_both"),Form("os_fvtxn_v2_both"),15, 0.0, 3.0, -1.1, 1.1);

  TProfile* os_bbcs_v3_west = new TProfile(Form("os_bbcs_v3_west"),Form("os_bbcs_v3_west"),15, 0.0, 3.0, -1.1, 1.1);
  TProfile* os_bbcs_v3_east = new TProfile(Form("os_bbcs_v3_east"),Form("os_bbcs_v3_east"),15, 0.0, 3.0, -1.1, 1.1);
  TProfile* os_bbcs_v3_both = new TProfile(Form("os_bbcs_v3_both"),Form("os_bbcs_v3_both"),15, 0.0, 3.0, -1.1, 1.1);

  TProfile* os_fvtxs_v3_west = new TProfile(Form("os_fvtxs_v3_west"),Form("os_fvtxs_v3_west"),15, 0.0, 3.0, -1.1, 1.1);
  TProfile* os_fvtxs_v3_east = new TProfile(Form("os_fvtxs_v3_east"),Form("os_fvtxs_v3_east"),15, 0.0, 3.0, -1.1, 1.1);
  TProfile* os_fvtxs_v3_both = new TProfile(Form("os_fvtxs_v3_both"),Form("os_fvtxs_v3_both"),15, 0.0, 3.0, -1.1, 1.1);

  TProfile* os_fvtxn_v3_west = new TProfile(Form("os_fvtxn_v3_west"),Form("os_fvtxn_v3_west"),15, 0.0, 3.0, -1.1, 1.1);
  TProfile* os_fvtxn_v3_east = new TProfile(Form("os_fvtxn_v3_east"),Form("os_fvtxn_v3_east"),15, 0.0, 3.0, -1.1, 1.1);
  TProfile* os_fvtxn_v3_both = new TProfile(Form("os_fvtxn_v3_both"),Form("os_fvtxn_v3_both"),15, 0.0, 3.0, -1.1, 1.1);

  TProfile* os_bbcs_c24 = new TProfile(Form("os_bbcs_c24"),Form("os_bbcs_c24"),1,-0.5,0.5,-1.1,1.1);
  TProfile* os_fvtxs_c24 = new TProfile(Form("os_fvtxs_c24"),Form("os_fvtxs_c24"),1,-0.5,0.5,-1.1,1.1);
  TProfile* os_fvtxn_c24 = new TProfile(Form("os_fvtxn_c24"),Form("os_fvtxn_c24"),1,-0.5,0.5,-1.1,1.1);
  TProfile* os_bbcs_c24_vs_cent = new TProfile(Form("os_bbcs_c24_vs_cent"),Form("os_bbcs_c24_vs_cent"),101,-0.5,100.5,-1.1,1.1);
  TProfile* os_fvtxs_c24_vs_cent = new TProfile(Form("os_fvtxs_c24_vs_cent"),Form("os_fvtxs_c24_vs_cent"),101,-0.5,100.5,-1.1,1.1);
  TProfile* os_fvtxn_c24_vs_cent = new TProfile(Form("os_fvtxn_c24_vs_cent"),Form("os_fvtxn_c24_vs_cent"),101,-0.5,100.5,-1.1,1.1);
  TProfile* os_bbcs_c24_vs_bbcqs = new TProfile(Form("os_bbcs_c24_vs_bbcqs"),Form("os_bbcs_c24_vs_bbcqs"),200,0,200,-1.1,1.1);
  TProfile* os_fvtxs_c24_vs_bbcqs = new TProfile(Form("os_fvtxs_c24_vs_bbcqs"),Form("os_fvtxs_c24_vs_bbcqs"),200,0,200,-1.1,1.1);
  TProfile* os_fvtxn_c24_vs_bbcqs = new TProfile(Form("os_fvtxn_c24_vs_bbcqs"),Form("os_fvtxn_c24_vs_bbcqs"),200,0,200,-1.1,1.1);
  TProfile* os_bbcs_c24_vs_nfvtxs = new TProfile(Form("os_bbcs_c24_vs_nfvtxs"),Form("os_bbcs_c24_vs_nfvtxs"),1000,-0.5,999.5,-1.1,1.1);
  TProfile* os_fvtxs_c24_vs_nfvtxs = new TProfile(Form("os_fvtxs_c24_vs_nfvtxs"),Form("os_fvtxs_c24_vs_nfvtxs"),1000,-0.5,999.5,-1.1,1.1);
  TProfile* os_fvtxn_c24_vs_nfvtxs = new TProfile(Form("os_fvtxn_c24_vs_nfvtxs"),Form("os_fvtxn_c24_vs_nfvtxs"),1000,-0.5,999.5,-1.1,1.1);
  TProfile* os_bbcs_d24_in_both = new TProfile(Form("os_bbcs_d24_in_both"),Form("os_bbcs_d24_in_both"),15,0.0,3.0,-1.1,1.1);
  TProfile* os_fvtxs_d24_in_both = new TProfile(Form("os_fvtxs_d24_in_both"),Form("os_fvtxs_d24_in_both"),15,0.0,3.0,-1.1,1.1);
  TProfile* os_fvtxn_d24_in_both = new TProfile(Form("os_fvtxn_d24_in_both"),Form("os_fvtxn_d24_in_both"),15,0.0,3.0,-1.1,1.1);
  TProfile* os_bbcs_d24_out_both = new TProfile(Form("os_bbcs_d24_out_both"),Form("os_bbcs_d24_out_both"),15,0.0,3.0,-1.1,1.1);
  TProfile* os_fvtxs_d24_out_both = new TProfile(Form("os_fvtxs_d24_out_both"),Form("os_fvtxs_d24_out_both"),15,0.0,3.0,-1.1,1.1);
  TProfile* os_fvtxn_d24_out_both = new TProfile(Form("os_fvtxn_d24_out_both"),Form("os_fvtxn_d24_out_both"),15,0.0,3.0,-1.1,1.1);

  // ---------------------------------------------------------------------------------------------------------



  //------------------------------------------------------------//
  //               Initializing Tree Variables                  //
  //------------------------------------------------------------//

  cout << "Now getting ready to read in the tree branch addresses and stuff...." << endl;



  //tree variables
  float        event;
  float        d_bbcz;    // bbcz
  float        centrality; // integer but stored as float in PHGlobal etc
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
  float        d_BBC_charge[64];

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

  // List of branches
  TBranch* b_event;   //!
  TBranch* b_bbc_z;   //!
  TBranch* b_centrality;   //!
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
  // TBranch* b_d_ntrk;   //!
  // TBranch* b_d_cntpx;   //!
  // TBranch* b_d_cntpy;   //!
  // TBranch* b_d_cntpz;   //!

  ntp_event_chain->SetBranchAddress("bbc_z",&d_bbcz,&b_bbc_z);
  ntp_event_chain->SetBranchAddress("centrality",&centrality,&b_centrality);
  ntp_event_chain->SetBranchAddress("bbc_qn",&bbc_qn,&b_bbc_qn);
  ntp_event_chain->SetBranchAddress("bbc_qs",&bbc_qs,&b_bbc_qs);
  ntp_event_chain->SetBranchAddress("npc1",&npc1,&b_npc1);
  ntp_event_chain->SetBranchAddress("event",&event,&b_event);
  ntp_event_chain->SetBranchAddress("trigger_scaled",&trigger_scaled,&b_trigger_scaled);
  ntp_event_chain->SetBranchAddress("trigger_live",&trigger_live,&b_trigger_live);
  ntp_event_chain->SetBranchAddress("bc_x",&bc_x,&b_bc_x);
  ntp_event_chain->SetBranchAddress("bc_y",&bc_y,&b_bc_y);
  ntp_event_chain->SetBranchAddress("vtx_z",&vtx_z,&b_vtx_z);

  ntp_event_chain->SetBranchAddress("fvtx_x",&eventfvtx_x,&b_fvtx_x);
  ntp_event_chain->SetBranchAddress("fvtx_y",&eventfvtx_y,&b_fvtx_y);
  ntp_event_chain->SetBranchAddress("fvtx_z",&eventfvtx_z,&b_fvtx_z);

  ntp_event_chain->SetBranchAddress("d_BBC_charge",d_BBC_charge,&b_d_BBC_charge);
  ntp_event_chain->SetBranchAddress("d_Qx",d_Qx,&b_d_Qx);
  ntp_event_chain->SetBranchAddress("d_Qy",d_Qy,&b_d_Qy);
  ntp_event_chain->SetBranchAddress("d_Qw",d_Qw,&b_d_Qw);

  ntp_event_chain->SetBranchAddress("d_nFVTX_clus",&d_nFVTX_clus,&b_d_nFVTX_clus);
  ntp_event_chain->SetBranchAddress("d_nFVTXN_clus",&d_nFVTXN_clus,&b_d_nFVTXN_clus);
  ntp_event_chain->SetBranchAddress("d_nFVTXS_clus",&d_nFVTXS_clus,&b_d_nFVTXS_clus);
  ntp_event_chain->SetBranchAddress("d_FVTX_x",d_FVTX_x,&b_d_FVTX_x);
  ntp_event_chain->SetBranchAddress("d_FVTX_y",d_FVTX_y,&b_d_FVTX_y);
  ntp_event_chain->SetBranchAddress("d_FVTX_z",d_FVTX_z,&b_d_FVTX_z);

  ntp_event_chain->SetBranchAddress("d_ntrk",&d_ntrk,&b_ntrk);
  ntp_event_chain->SetBranchAddress("d_cntpx",d_px,&b_px);
  ntp_event_chain->SetBranchAddress("d_cntpy",d_py,&b_py);
  ntp_event_chain->SetBranchAddress("d_cntpz",d_pz,&b_pz);




  //------------------------------------------------------------//
  //          Finished Initializing Tree Variables              //
  //------------------------------------------------------------//

  //------------------------------------------------------------//
  //                   Looping Over Event Tree                  //
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
      if(rp_recal_pass<1 || rp_recal_pass > 3) break;// rp_recal_pass only valid between 1 and 3

      //if ( ievt >= 100000 ) break; // just 100k events for testing, runs a little on the slow side...
      ++all_counter;

      bool say_event = ( ievt%1000==0 );

      if ( say_event ) cout << "event number = " << ievt << endl;

      if ( ( say_event && verbosity > 0 ) || verbosity > 1 ) cout << "getting event level variables" << endl;
      ntp_event_chain->GetEntry(ievt);
      if ( ( say_event && verbosity > 0 ) || verbosity > 1 ) cout << "Finished getting tree variables" << endl;

      // ---------------------
      // --- trigger selection
      // ---------------------

      unsigned int trigger_FVTXNSBBCScentral = 0x00100000;
      unsigned int trigger_FVTXNSBBCS        = 0x00400000;
      unsigned int trigger_BBCLL1narrowcent  = 0x00000008;
      unsigned int trigger_BBCLL1narrow      = 0x00000010;

      unsigned int accepted_triggers = 0;
      // --- Run16dAu200
      if ( runNumber >= 454774 && runNumber <= 455639 ) accepted_triggers = trigger_BBCLL1narrowcent | trigger_BBCLL1narrow;
      // --- Run16dAu62
      if ( runNumber >= 455792 && runNumber <= 456283 ) accepted_triggers = trigger_BBCLL1narrowcent | trigger_BBCLL1narrow;
      // --- Run16dAu20
      if ( runNumber >= 456652 && runNumber <= 457298 ) accepted_triggers = trigger_FVTXNSBBCScentral | trigger_FVTXNSBBCS;
      // --- Run16dAu39
      if ( runNumber >= 457634 && runNumber <= 458167 ) accepted_triggers = trigger_FVTXNSBBCScentral | trigger_FVTXNSBBCS;

      unsigned int passes_trigger = trigger_scaled & accepted_triggers;
      if ( passes_trigger == 0 )
        {
          if ( verbosity > 0 ) cout << "trigger rejected" << endl;
          continue;
        }

      double ZVTX = -9999;
      if ( runNumber >= 454774 && runNumber <= 456283 ) ZVTX = d_bbcz;
      if ( runNumber >= 456652 && runNumber <= 458167 ) ZVTX = eventfvtx_z;
      if ( fabs(ZVTX) > 10.0 )
        {
          if ( verbosity > 0 ) cout << "vertex rejected" << endl;
          continue;
        }
      // --- this cut might be a good idea but it throws out too many events, further study needed
      // if ( d_bbcz > -999 && eventfvtx_z > -999 && fabs(d_bbcz-eventfvtx_z) > 5 )
      //   {
      //     if ( verbosity > 0 ) cout << "bbc and fvtx vertex exist but out of range of each other " << d_bbcz << " " << eventfvtx_z << endl;
      //     ++bad_vertex_counter;
      //     continue;
      //   }

      // make sure bin number doesn't exceed number of bins
      int izvtx = NZPS*(ZVTX+10)/20;
      if ( izvtx < 0 || izvtx >= NZPS )
        {
          cout << "z vertex bin count problem!!!!" << endl;
          cout << "bbcz = " << d_bbcz << endl;
          cout << "fvtx_z = " << eventfvtx_z << endl;
          cout << "bin number is " << izvtx << endl;
          continue;
        }

      int toomanyclusters = 9999;
      // --- Run16dAu200
      if ( runNumber >= 454774 && runNumber <= 455639 ) toomanyclusters = 1000;
      // --- Run16dAu62
      if ( runNumber >= 455792 && runNumber <= 456283 ) toomanyclusters = 1000;
      // --- Run16dAu20
      if ( runNumber >= 456652 && runNumber <= 457298 ) toomanyclusters = 500;
      // --- Run16dAu39
      if ( runNumber >= 457634 && runNumber <= 458167 ) toomanyclusters = 300;

      bool is_okaync = ( d_nFVTX_clus < toomanyclusters );
      if ( !is_okaync )
        {
          if ( verbosity > 0 ) cout << "too many clusters" << endl;
          continue;
        }

      // ---------------------------------------------------------------------------------------

      //------------------------------------------------------------//
      //                Calculating Event Planes                    //
      //------------------------------------------------------------//

      if ( ( say_event && verbosity > 0 ) || verbosity > 1 ) cout << "Calculating event planes" << endl;

      // --- all numbers from Darren 2016-06-23
      const float x_off = 0.3;
      const float beam_angle = 0.001;
      float vtx_z = d_bbcz;
      if ( eventfvtx_z > -999 ) vtx_z = eventfvtx_z;
      float vtx_x = x_off + atan(beam_angle)*vtx_z;
      float vtx_y = 0.02;

      // --- radius cut using FVTX coordinates
      if ( sqrt(pow(eventfvtx_x-vtx_x,2.0) +  pow(eventfvtx_y-vtx_y,2.0)) >= 0.15 )
        {
          if ( verbosity > 0 ) cout << "rejecting event due to radius cut" << endl;
          continue;
        }

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

      if ( ( say_event && verbosity > 0 ) || verbosity > 1 ) cout << "Looping over BBC stuff now" << endl;

      if ( bbc_pmts )
        {
          for(int ipmt = 0; ipmt < 64; ipmt++)
            {
              float bbc_charge = d_BBC_charge[ipmt];
              if(bbc_charge <= 0) continue;

              float bbc_x      = d_pmt_x[ipmt] - vtx_x*10;//pmt location in mm
              float bbc_y      = d_pmt_y[ipmt] - vtx_y*10;
              float bbc_z      = d_pmt_z       - vtx_z*10;

              // --- rotation
              bbc_x = bbc_z*sin(-beam_angle) + bbc_x*cos(-beam_angle);

              float phi = TMath::ATan2(bbc_y,bbc_x);

              int ring = get_pmt_layer(ipmt);
              th1d_bbc_charge_phi->Fill(phi,bbc_charge);
              if ( ring == 0 ) th1d_bbc0_charge_phi->Fill(phi,bbc_charge);
              if ( ring == 1 ) th1d_bbc1_charge_phi->Fill(phi,bbc_charge);
              if ( ring == 2 ) th1d_bbc2_charge_phi->Fill(phi,bbc_charge);
              if ( ring == 3 ) th1d_bbc3_charge_phi->Fill(phi,bbc_charge);
              if ( ring == 4 ) th1d_bbc4_charge_phi->Fill(phi,bbc_charge);

              bbc_qx2 += bbc_charge*TMath::Cos(2*phi);
              bbc_qy2 += bbc_charge*TMath::Sin(2*phi);
              bbc_qx3 += bbc_charge*TMath::Cos(3*phi);
              bbc_qy3 += bbc_charge*TMath::Sin(3*phi);
              bbc_qx4 += bbc_charge*TMath::Cos(4*phi);
              bbc_qy4 += bbc_charge*TMath::Sin(4*phi);
              bbc_qw += bbc_charge;
            } // loop over tubes
        } // check on tubes

      // --- do centrality cut here!!!

      if ( centrality > -999 )
        {
          if ( runNumber >= 454744 && runNumber <= 455639 && centrality > 5  ) continue; // dAu 200 GeV
          if ( runNumber >= 455792 && runNumber <= 456283 && centrality > 10 ) continue; // dAu 62 GeV
          if ( runNumber >= 456652 && runNumber <= 457298 && centrality > 20 ) continue; // dAu 20 GeV
          if ( runNumber >= 457634 && runNumber <= 458167 && centrality > 20 ) continue; // dAu 39 GeV
        }
      else
        {
          //cout << "centrality undefined, cutting on bbc charge" << endl;
          // --- revise these numbers as needed
          if ( runNumber >= 454744 && runNumber <= 455639 && bbc_qw < 60.0 ) continue; // dAu 200 GeV
          if ( runNumber >= 455792 && runNumber <= 456283 && bbc_qw < 40.0 ) continue; // dAu 62 GeV
          if ( runNumber >= 456652 && runNumber <= 457298 && bbc_qw < 25.0 ) continue; // dAu 20 GeV
          if ( runNumber >= 457634 && runNumber <= 458167 && bbc_qw < 30.0 ) continue; // dAu 39 GeV
          ++bad_cent_counter;
        }

      if ( say_event )
        {
          cout << "bbc charge = " << bbc_qw << endl;
          cout << "centrality = " << centrality << endl;
        }

      bbcs_cos22->Fill(0.0,bbc_qx2);
      bbcs_sin22->Fill(0.0,bbc_qy2);
      bbcs_cos32->Fill(0.0,bbc_qx3);
      bbcs_sin32->Fill(0.0,bbc_qy3);
      float bbc_qq2 = ( (bbc_qx2*bbc_qx2) + (bbc_qy2*bbc_qy2) - bbc_qw ) / ( (bbc_qw*bbc_qw) - bbc_qw );
      float bbc_qq3 = ( (bbc_qx3*bbc_qx3) + (bbc_qy3*bbc_qy3) - bbc_qw ) / ( (bbc_qw*bbc_qw) - bbc_qw );
      bbcs_c22->Fill(0.0,bbc_qq2);
      bbcs_c32->Fill(0.0,bbc_qq3);

      //cout << "HELLO HERE I AM" << endl;

      ++event_counter;

      th1d_BBC_charge->Fill(bbc_qw);
      th1d_FVTX_nclus->Fill(d_nFVTX_clus);
      th2d_qBBC_nFVTX->Fill(bbc_qw,d_nFVTX_clus);

      // --------------
      // --- FVTX stuff
      // --------------

      float fvtxs_qx2[5];//all layers then 0 1 2 3
      float fvtxs_qy2[5];
      float fvtxs_qx3[5];//all layers then 0 1 2 3
      float fvtxs_qy3[5];
      float fvtxs_qx4[5];//all layers then 0 1 2 3
      float fvtxs_qy4[5];
      float fvtxs_qw[5];

      for(int ilayer = 0; ilayer < 5; ilayer++)
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

      float fvtxn_qx2[5];//all layers then 0 1 2 3
      float fvtxn_qy2[5];
      float fvtxn_qx3[5];//all layers then 0 1 2 3
      float fvtxn_qy3[5];
      float fvtxn_qx4[5];//all layers then 0 1 2 3
      float fvtxn_qy4[5];
      float fvtxn_qw[5];

      for(int ilayer = 0; ilayer < 5; ilayer++)
        {
          fvtxn_qx2[ilayer] = 0.0;
          fvtxn_qy2[ilayer] = 0.0;
          fvtxn_qx3[ilayer] = 0.0;
          fvtxn_qy3[ilayer] = 0.0;
          fvtxn_qx4[ilayer] = 0.0;
          fvtxn_qy4[ilayer] = 0.0;
          fvtxn_qw[ilayer] = 0.0;
        } // loop over layers


      if ( ( say_event && verbosity > 0 ) || verbosity > 1 ) cout << "Looping over FVTX cluster" << endl;
      if ( fvtx_clusters )
        {
          for(int iclus = 0; iclus < d_nFVTX_clus; iclus++)
            {
              float fvtx_x      = d_FVTX_x[iclus] - vtx_x;
              float fvtx_y      = d_FVTX_y[iclus] - vtx_y;
              float fvtx_z      = d_FVTX_z[iclus] - vtx_z;

              // --- rotation
              fvtx_x = fvtx_z*sin(-beam_angle) + fvtx_x*cos(-beam_angle);

              double fvtx_r = sqrt(pow(fvtx_x,2.0)+pow(fvtx_y,2.0));
              double fvtx_the = atan2(fvtx_r,fvtx_z);
              double fvtx_eta = -log(tan(0.5*fvtx_the));

              // --- determine layer based on cluster z
              int fvtx_layer = get_fvtx_layer(d_FVTX_z[iclus]); // raw z to get layer

              // --- gap cut removes events where fvtx z vertex is right below fvtx south
              int igap = (fabs(fvtx_eta)-1.0)/0.5;
              int id_fvtx = fvtx_layer*5+igap;
              ++cluster_counter;
              if(!(id_fvtx>=0 && id_fvtx<40))
                {
                  if ( verbosity > 0 )
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

              double FVTX_r = sqrt(pow(d_FVTX_x[iclus],2.0)+pow(d_FVTX_y[iclus],2.0));
              if ( runNumber >= 456652 && runNumber <= 458167 && FVTX_r < 5.2 ) continue;

              float phi = TMath::ATan2(fvtx_y,fvtx_x);
              int izvtx2 = izvtx/2;
              int phi_bin = th1d_fvtxs_phi_weight[izvtx2][fvtx_layer+1]->FindBin(phi); // COME BACK HERE AND HAVE A LOOK
              float fvtx_weight = th1d_fvtxs_phi_weight[izvtx2][fvtx_layer+1]->GetBinContent(phi_bin);
              //float fvtx_weight = 1.0; // don't want to use weights right now...

              // --- south side
              if ( d_FVTX_z[iclus] < 0 )
                {
                  fvtxs_qx2[fvtx_layer+1] += fvtx_weight * TMath::Cos(2*phi);
                  fvtxs_qy2[fvtx_layer+1] += fvtx_weight * TMath::Sin(2*phi);
                  fvtxs_qx3[fvtx_layer+1] += fvtx_weight * TMath::Cos(3*phi);
                  fvtxs_qy3[fvtx_layer+1] += fvtx_weight * TMath::Sin(3*phi);

                  fvtxs_qx2[0] += fvtx_weight * TMath::Cos(2*phi);
                  fvtxs_qy2[0] += fvtx_weight * TMath::Sin(2*phi);
                  fvtxs_qx3[0] += fvtx_weight * TMath::Cos(3*phi);
                  fvtxs_qy3[0] += fvtx_weight * TMath::Sin(3*phi);
                  fvtxs_qx4[0] += fvtx_weight * TMath::Cos(4*phi);
                  fvtxs_qy4[0] += fvtx_weight * TMath::Sin(4*phi);

                  fvtxs_qw[fvtx_layer+1] += fvtx_weight;
                  fvtxs_qw[0] += fvtx_weight;

                  th1d_fvtxs_clus_phi->Fill(phi);
                  if ( fvtx_layer == 0 ) th1d_fvtxs0_clus_phi->Fill(phi);
                  if ( fvtx_layer == 1 ) th1d_fvtxs1_clus_phi->Fill(phi);
                  if ( fvtx_layer == 2 ) th1d_fvtxs2_clus_phi->Fill(phi);
                  if ( fvtx_layer == 3 ) th1d_fvtxs3_clus_phi->Fill(phi);
                } // check on south

              // --- north side
              if ( d_FVTX_z[iclus] > 0 )
                {
                  fvtxn_qx2[fvtx_layer+1] += TMath::Cos(2*phi);
                  fvtxn_qy2[fvtx_layer+1] += TMath::Sin(2*phi);
                  fvtxn_qx3[fvtx_layer+1] += TMath::Cos(3*phi);
                  fvtxn_qy3[fvtx_layer+1] += TMath::Sin(3*phi);

                  fvtxn_qx2[0] += TMath::Cos(2*phi);
                  fvtxn_qy2[0] += TMath::Sin(2*phi);
                  fvtxn_qx3[0] += TMath::Cos(3*phi);
                  fvtxn_qy3[0] += TMath::Sin(3*phi);
                  fvtxn_qx4[0] += TMath::Cos(4*phi);
                  fvtxn_qy4[0] += TMath::Sin(4*phi);

                  fvtxn_qw[fvtx_layer+1] ++;
                  fvtxn_qw[0] ++;

                  th1d_fvtxn_clus_phi->Fill(phi);
                  if ( fvtx_layer == 0 ) th1d_fvtxn0_clus_phi->Fill(phi);
                  if ( fvtx_layer == 1 ) th1d_fvtxn1_clus_phi->Fill(phi);
                  if ( fvtx_layer == 2 ) th1d_fvtxn2_clus_phi->Fill(phi);
                  if ( fvtx_layer == 3 ) th1d_fvtxn3_clus_phi->Fill(phi);
                } // check on north

            } // loop over cluster
        } // check on clusters

      th1d_FVTXS_nclus->Fill(fvtxs_qw[0]);
      th1d_FVTXN_nclus->Fill(fvtxn_qw[0]);

      fvtxs_cos22->Fill(0.0,fvtxs_qx2[0]);
      fvtxs_sin22->Fill(0.0,fvtxs_qy2[0]);
      fvtxs_cos32->Fill(0.0,fvtxs_qx3[0]);
      fvtxs_sin32->Fill(0.0,fvtxs_qy3[0]);
      float fvtxs_qq2 = ( (fvtxs_qx2[0]*fvtxs_qx2[0]) + (fvtxs_qy2[0]*fvtxs_qy2[0]) - fvtxs_qw[0] ) / ( (fvtxs_qw[0]*fvtxs_qw[0]) - fvtxs_qw[0] );
      float fvtxs_qq3 = ( (fvtxs_qx3[0]*fvtxs_qx3[0]) + (fvtxs_qy3[0]*fvtxs_qy3[0]) - fvtxs_qw[0] ) / ( (fvtxs_qw[0]*fvtxs_qw[0]) - fvtxs_qw[0] );
      fvtxs_c22->Fill(0.0,fvtxs_qq2);
      fvtxs_c32->Fill(0.0,fvtxs_qq3);


      // --- the array that has all of the Q vectors
      float sumxy[NHAR][NDET][4];
      for (int i=0; i<NHAR; i++)
        {
          for (int j=0; j<NDET; j++)
            {
              for (int k=0; k<4; k++)//qx qy qw psi
                {
                  sumxy[i][j][k]=0; // initialize to 0
                } // x,y,w,psi
            } // detectors
        } // harmonics

      //save bbc q vec 2
      //cout<<"bbc from node tree: "<<d_Qx[5]<<" "<<d_Qy[5]<<" "<<d_Qw[5]<<endl;
      // --- should add cut for only assign if Qw > 0
      sumxy[1][0][0] = d_Qx[5];
      sumxy[1][0][1] = d_Qy[5];
      sumxy[1][0][2] = d_Qw[5];

      //save fvtx q vec 2
      // --- should add cut for only assign if Qw > 4 and Qw < 1000 (or some other large number, 1000 is for pAu)
      sumxy[1][1][0] = d_Qx[4];
      sumxy[1][1][1] = d_Qy[4];
      sumxy[1][1][2] = d_Qw[4];



      if ( bbc_pmts )
        {
          sumxy[1][south_bbc_angle][0] = bbc_qx2;
          sumxy[1][south_bbc_angle][1] = bbc_qy2;
          sumxy[1][south_bbc_angle][2] = bbc_qw;
          sumxy[2][south_bbc_angle][0] = bbc_qx3;
          sumxy[2][south_bbc_angle][1] = bbc_qy3;
          sumxy[2][south_bbc_angle][2] = bbc_qw;
        }

      if ( fvtx_clusters )
        {
          // --- south side, 2nd harmonic

          sumxy[1][south_fvtx_angle][0] = fvtxs_qx2[0];
          sumxy[1][south_fvtx_angle][1] = fvtxs_qy2[0];
          sumxy[1][south_fvtx_angle][2] = fvtxs_qw[0];

          sumxy[1][south_fvtx_0_angle][0] = fvtxs_qx2[1];
          sumxy[1][south_fvtx_0_angle][1] = fvtxs_qy2[1];
          sumxy[1][south_fvtx_0_angle][2] = fvtxs_qw[1];

          sumxy[1][south_fvtx_1_angle][0] = fvtxs_qx2[2];
          sumxy[1][south_fvtx_1_angle][1] = fvtxs_qy2[2];
          sumxy[1][south_fvtx_1_angle][2] = fvtxs_qw[2];

          sumxy[1][south_fvtx_2_angle][0] = fvtxs_qx2[3];
          sumxy[1][south_fvtx_2_angle][1] = fvtxs_qy2[3];
          sumxy[1][south_fvtx_2_angle][2] = fvtxs_qw[3];

          sumxy[1][south_fvtx_3_angle][0] = fvtxs_qx2[4];
          sumxy[1][south_fvtx_3_angle][1] = fvtxs_qy2[4];
          sumxy[1][south_fvtx_3_angle][2] = fvtxs_qw[4];

          // --- north side, 2nd harmonic

          sumxy[1][north_fvtx_angle][0] = fvtxn_qx2[0];
          sumxy[1][north_fvtx_angle][1] = fvtxn_qy2[0];
          sumxy[1][north_fvtx_angle][2] = fvtxn_qw[0];

          sumxy[1][north_fvtx_0_angle][0] = fvtxn_qx2[1];
          sumxy[1][north_fvtx_0_angle][1] = fvtxn_qy2[1];
          sumxy[1][north_fvtx_0_angle][2] = fvtxn_qw[1];

          sumxy[1][north_fvtx_1_angle][0] = fvtxn_qx2[2];
          sumxy[1][north_fvtx_1_angle][1] = fvtxn_qy2[2];
          sumxy[1][north_fvtx_1_angle][2] = fvtxn_qw[2];

          sumxy[1][north_fvtx_2_angle][0] = fvtxn_qx2[3];
          sumxy[1][north_fvtx_2_angle][1] = fvtxn_qy2[3];
          sumxy[1][north_fvtx_2_angle][2] = fvtxn_qw[3];

          sumxy[1][north_fvtx_3_angle][0] = fvtxn_qx2[4];
          sumxy[1][north_fvtx_3_angle][1] = fvtxn_qy2[4];
          sumxy[1][north_fvtx_3_angle][2] = fvtxn_qw[4];

          // --- south side, 3rd harmonic

          sumxy[2][south_fvtx_angle][0] = fvtxs_qx3[0];
          sumxy[2][south_fvtx_angle][1] = fvtxs_qy3[0];
          sumxy[2][south_fvtx_angle][2] = fvtxs_qw[0];

          sumxy[2][south_fvtx_0_angle][0] = fvtxs_qx3[1];
          sumxy[2][south_fvtx_0_angle][1] = fvtxs_qy3[1];
          sumxy[2][south_fvtx_0_angle][2] = fvtxs_qw[1];

          sumxy[2][south_fvtx_1_angle][0] = fvtxs_qx3[2];
          sumxy[2][south_fvtx_1_angle][1] = fvtxs_qy3[2];
          sumxy[2][south_fvtx_1_angle][2] = fvtxs_qw[2];

          sumxy[2][south_fvtx_2_angle][0] = fvtxs_qx3[3];
          sumxy[2][south_fvtx_2_angle][1] = fvtxs_qy3[3];
          sumxy[2][south_fvtx_2_angle][2] = fvtxs_qw[3];

          sumxy[2][south_fvtx_3_angle][0] = fvtxs_qx3[4];
          sumxy[2][south_fvtx_3_angle][1] = fvtxs_qy3[4];
          sumxy[2][south_fvtx_3_angle][2] = fvtxs_qw[4];

          // --- north side, 3rd harmonic

          sumxy[2][north_fvtx_angle][0] = fvtxn_qx3[0];
          sumxy[2][north_fvtx_angle][1] = fvtxn_qy3[0];
          sumxy[2][north_fvtx_angle][2] = fvtxn_qw[0];

          sumxy[2][north_fvtx_0_angle][0] = fvtxn_qx3[1];
          sumxy[2][north_fvtx_0_angle][1] = fvtxn_qy3[1];
          sumxy[2][north_fvtx_0_angle][2] = fvtxn_qw[1];

          sumxy[2][north_fvtx_1_angle][0] = fvtxn_qx3[2];
          sumxy[2][north_fvtx_1_angle][1] = fvtxn_qy3[2];
          sumxy[2][north_fvtx_1_angle][2] = fvtxn_qw[2];

          sumxy[2][north_fvtx_2_angle][0] = fvtxn_qx3[3];
          sumxy[2][north_fvtx_2_angle][1] = fvtxn_qy3[3];
          sumxy[2][north_fvtx_2_angle][2] = fvtxn_qw[3];

          sumxy[2][north_fvtx_3_angle][0] = fvtxn_qx3[4];
          sumxy[2][north_fvtx_3_angle][1] = fvtxn_qy3[4];
          sumxy[2][north_fvtx_3_angle][2] = fvtxn_qw[4];
        }


      if ( DIAG )
        {
          cout<<"bbc from node tree: "<<d_Qx[5]<<" "<<d_Qy[5]<<" "<<d_Qw[5]<<endl;
          cout<<"bbc from me: "<<bbc_qx2<<" "<<bbc_qy2<<" "<<bbc_qw<<endl;

          cout<<"fvtx raw: "<<endl;
          cout<<"from node tree: "<<d_Qx[4]<<" "<<d_Qy[4]<<" "<<d_Qw[4]<<endl;
          cout<<"from clusters: " <<fvtxs_qx2[0]<<" "<<fvtxs_qy2[0]<<" "<<fvtxs_qw[0]<<endl;
        }

      for (int ih=1; ih<NHAR; ih++)
        {
          for (int id=0; id<NDET; id++)
            {
              if(sumxy[ih][id][2]>0)
                {
                  //float psi = atan2(sumxy[ih][id][1],sumxy[ih][id][0])/2.0;
                  float psi = atan2(sumxy[ih][id][1],sumxy[ih][id][0])/float(ih+1);
                  if ( DIAG ) cout<<"RAW: for id: "<<id<<" psi: "<<psi<<endl;
                  psi_bf[ic][ih][id]->Fill(izvtx,psi);
                }
            }
        }

      //------------------------------------------------------------//
      //                Flattening iteration                        //
      //------------------------------------------------------------//
      //int icent = 0;
      for (int ih=1; ih<NHAR; ih++)
        {
          for(int id=0;id<NDET;id++)
            {
              if (sumxy[ih][id][2]>0.0)
                {
                  sumxy[ih][id][3]=atan2(sumxy[ih][id][1],sumxy[ih][id][0])/(ih+1.0);
                  if (rp_recal_pass>0) dis[icent][ih][id]->Fill(izvtx,sumxy[ih][id][3]*(ih+1.0));
                }
              if (sumxy[ih][id][2]>0.0) // check on weight (x,y,w,psi)
                {
                  for (int ib=0; ib<2; ib++)
                    {
                      sumxy[ih][id][ib]/=sumxy[ih][id][2]; // normalize to the weight
                      //if(ih==1 && id==0 && ib==0 && sumxy[ih][id][ib]>1) cout<<sumxy[ih][id][ib]<<endl;
                      if (rp_recal_pass>0)
                        {
                          ave[icent][izvtx][ih][id]->Fill(ib+0.0,sumxy[ih][id][ib]);
                          if(id==0 && DIAG) cout<<"filled ave: "<<ih<<" "<<id<<" "<<ib<<" with: "<<sumxy[ih][id][ib]<<endl;
                          if(ib==0) qx[icent][ih][id]->Fill(izvtx,sumxy[ih][id][0]);
                          if(ib==1) qy[icent][ih][id]->Fill(izvtx,sumxy[ih][id][1]);
                        } // pass > 0
                      float sxy=sumxy[ih][id][ib];
                      float mxy=mean[icent][izvtx][ih][id][ib]; // for recentering qx and qy (???)
                      float wxy=widt[icent][izvtx][ih][id][ib]; // for recentering qx and qy (???)

                      //if(icent==0 && izvtx==0 && ih==1 && id==0) cout<<ib<<" "<<sxy<<" "<<mxy<<" "<<wxy<<endl;
                      sumxy[ih][id][ib]=(sxy-mxy)/wxy; // recentered by mean and renormalized to width
                      if (rp_recal_pass>0)
                        {
                          ave[icent][izvtx][ih][id]->Fill(ib+2.0,sumxy[ih][id][ib]);  // ib+2 to avoid overlap
                          if(id==0 && DIAG) cout<<"filled ave2: "<<ih<<" "<<id<<" "<<ib<<" with: "<<sumxy[ih][id][ib]<<endl;
                          if(ib==0) qx[icent][ih][id]->Fill(izvtx+NZPS,sumxy[ih][id][0]);
                          if(ib==1) qy[icent][ih][id]->Fill(izvtx+NZPS,sumxy[ih][id][1]);
                        } // pass > 0
                    } // if weight > 0

                  sumxy[ih][id][3]=atan2(sumxy[ih][id][1],sumxy[ih][id][0])/(ih+1.0);
                  if (rp_recal_pass>0)
                    {
                      // fill histogram with psi calculated with recenter q vectors
                      dis[icent][ih][id]->Fill(izvtx+NZPS,sumxy[ih][id][3]*(ih+1.0));
                    }

                  float psi=sumxy[ih][id][3]*(ih+1.0);
                  if(ih==1 && id==0 && DIAG)  cout<<"psi-1 bbc: "<<psi<<endl;
                  float dp=0.0;
                  // --- flattening part, fourier components of psi distribution
                  for (int io=0; io<NORD; io++)
                    {
                      float cc=cos((io+1.0)*psi);
                      float ss=sin((io+1.0)*psi);
                      // first set of fourier components of psi
                      if (rp_recal_pass>0) flt[icent][izvtx][ih][id]->Fill(io+0.0,cc);
                      if (rp_recal_pass>0) flt[icent][izvtx][ih][id]->Fill(io+NORD,ss);
                      // --- four means fourier
                      float aa=four[icent][izvtx][ih][id][0][io]; // mean cos
                      float bb=four[icent][izvtx][ih][id][1][io]; // mean sin
                      // dp is offset to psi, aa and bb are zero in first pass, non zero later
                      dp+=(aa*ss-bb*cc)*2.0/(io+1.0); // ( trig identity cos(A+B) = cosAsinB - cosBsinA )
                    } // orders
                  psi+=dp; // shift psi by...
                  psi=atan2(sin(psi),cos(psi)); // trick to readjust the range
                  if(ih==1 && id==0 && DIAG)  cout<<"psi-2 bbc: "<<psi<<endl;
                  for (int io=0; io<NORD; io++)
                    {
                      float cc=cos((io+1.0)*psi);
                      float ss=sin((io+1.0)*psi);
                      // --- fourier components of modified psi
                      if (rp_recal_pass>0) flt[icent][izvtx][ih][id]->Fill(io+NORD*2.0,cc);
                      if (rp_recal_pass>0) flt[icent][izvtx][ih][id]->Fill(io+NORD*3.0,ss);
                    }
                  sumxy[ih][id][3]=psi/(ih+1.0);
                  if (rp_recal_pass>0) dis[icent][ih][id]->Fill(izvtx+NZPS*2.0,sumxy[ih][id][3]*(ih+1.0));
                } // end if weight > 0
              else
                {
                  sumxy[ih][id][3]=-9999.9;
                } // otherwise set psi to some crazy number
            } // end of detectors
        } // end of harmonics

      if ( DIAG ) cout<<"bbc_rp2: "<<sumxy[1][0][3]<<endl;

      for (int ih=1; ih<NHAR; ih++)
        {
          for (int id=0; id<NDET; id++)
            {
              if(sumxy[ih][id][2]>0)
                {
                  psi_af[ic][ih][id]->Fill(izvtx,sumxy[ih][id][3]);
                  if ( DIAG ) cout<<"CORR: for id: "<<id<<" psi: "<<sumxy[ih][id][3]<<endl;
                }
            }
        }

      // ---
      // --- now going to calculate v2
      // ---
      if ( rp_recal_pass < 3 ) continue; // don't calculate v2 except for final pass

      // --- let's see if we can figure out the mean qx and qy from the flattening procedure above...
      // --- this one? mean[icent][izvtx][ih][id][ib] typical agreement better than 1%
      // --- or maybe this one? sumxy[ih][id][ib] actually i think these are recentered Q-vectors, but if that's the case i could use them directly...
      // --- adding fourth harmonic to the flattening procedure for the first time... 20160803-1702 EDT
      float os_bbc_qw = bbc_qw; // probably the same...
      float os_bbc_qx2 = mean[icent][izvtx][1][2][0];
      float os_bbc_qy2 = mean[icent][izvtx][1][2][1];
      float os_bbc_qx3 = mean[icent][izvtx][2][2][0];
      float os_bbc_qy3 = mean[icent][izvtx][2][2][1];
      float os_bbc_qx4 = mean[icent][izvtx][3][2][0];
      float os_bbc_qy4 = mean[icent][izvtx][3][2][1];
      os_bbcs_cos22->Fill(0.0,os_bbc_qx2);
      os_bbcs_sin22->Fill(0.0,os_bbc_qy2);
      os_bbcs_cos32->Fill(0.0,os_bbc_qx3);
      os_bbcs_sin32->Fill(0.0,os_bbc_qy3);
      os_bbc_qx2 = bbc_qx2 - os_bbc_qw*os_bbc_qx2;
      os_bbc_qy2 = bbc_qy2 - os_bbc_qw*os_bbc_qy2;
      os_bbc_qx3 = bbc_qx3 - os_bbc_qw*os_bbc_qx3;
      os_bbc_qy3 = bbc_qy3 - os_bbc_qw*os_bbc_qy3;
      os_bbc_qx4 = bbc_qx4 - os_bbc_qw*os_bbc_qx4;
      os_bbc_qy4 = bbc_qy4 - os_bbc_qw*os_bbc_qy4;
      float os_bbc_qq2 = ( (os_bbc_qx2*os_bbc_qx2) + (os_bbc_qy2*os_bbc_qy2) - os_bbc_qw ) / ( (os_bbc_qw*os_bbc_qw) - os_bbc_qw );
      float os_bbc_qq3 = ( (os_bbc_qx3*os_bbc_qx3) + (os_bbc_qy3*os_bbc_qy3) - os_bbc_qw ) / ( (os_bbc_qw*os_bbc_qw) - os_bbc_qw );
      // float os_bbc_qq4 = ( (os_bbc_qx4*os_bbc_qx4) + (os_bbc_qy4*os_bbc_qy4) - os_bbc_qw ) / ( (os_bbc_qw*os_bbc_qw) - os_bbc_qw );
      os_bbcs_c22->Fill(0.0,os_bbc_qq2);
      os_bbcs_c32->Fill(0.0,os_bbc_qq3);
      float os_bbc_psi2 = atan2(os_bbc_qy2,os_bbc_qx2)/2.0;
      float os_bbc_psi3 = atan2(os_bbc_qy3,os_bbc_qx3)/3.0;
      os_bbcs_1dPsi2->Fill(os_bbc_psi2);
      os_bbcs_1dPsi3->Fill(os_bbc_psi3);

      float os_fvtxs_qw = fvtxs_qw[0];
      float os_fvtxs_qx2 = mean[icent][izvtx][1][3][0];
      float os_fvtxs_qy2 = mean[icent][izvtx][1][3][1];
      float os_fvtxs_qx3 = mean[icent][izvtx][2][3][0];
      float os_fvtxs_qy3 = mean[icent][izvtx][2][3][1];
      float os_fvtxs_qx4 = mean[icent][izvtx][3][3][0];
      float os_fvtxs_qy4 = mean[icent][izvtx][3][3][1];
      os_fvtxs_cos22->Fill(0.0,os_fvtxs_qx2);
      os_fvtxs_sin22->Fill(0.0,os_fvtxs_qy2);
      os_fvtxs_cos32->Fill(0.0,os_fvtxs_qx3);
      os_fvtxs_sin32->Fill(0.0,os_fvtxs_qy3);
      os_fvtxs_qx2 = fvtxs_qx2[0] - os_fvtxs_qw*os_fvtxs_qx2;
      os_fvtxs_qy2 = fvtxs_qy2[0] - os_fvtxs_qw*os_fvtxs_qy2;
      os_fvtxs_qx3 = fvtxs_qx3[0] - os_fvtxs_qw*os_fvtxs_qx3;
      os_fvtxs_qy3 = fvtxs_qy3[0] - os_fvtxs_qw*os_fvtxs_qy3;
      os_fvtxs_qx4 = fvtxs_qx4[0] - os_fvtxs_qw*os_fvtxs_qx4;
      os_fvtxs_qy4 = fvtxs_qy4[0] - os_fvtxs_qw*os_fvtxs_qy4;
      float os_fvtxs_qq2 = ( (os_fvtxs_qx2*os_fvtxs_qx2) + (os_fvtxs_qy2*os_fvtxs_qy2) - os_fvtxs_qw ) / ( (os_fvtxs_qw*os_fvtxs_qw) - os_fvtxs_qw );
      float os_fvtxs_qq3 = ( (os_fvtxs_qx3*os_fvtxs_qx3) + (os_fvtxs_qy3*os_fvtxs_qy3) - os_fvtxs_qw ) / ( (os_fvtxs_qw*os_fvtxs_qw) - os_fvtxs_qw );
      // float os_fvtxs_qq4 = ( (os_fvtxs_qx4*os_fvtxs_qx4) + (os_fvtxs_qy4*os_fvtxs_qy4) - os_fvtxs_qw ) / ( (os_fvtxs_qw*os_fvtxs_qw) - os_fvtxs_qw );
      os_fvtxs_c22->Fill(0.0,os_fvtxs_qq2);
      os_fvtxs_c32->Fill(0.0,os_fvtxs_qq3);
      float os_fvtxs_psi2 = atan2(os_fvtxs_qy2,os_fvtxs_qx2)/2.0;
      float os_fvtxs_psi3 = atan2(os_fvtxs_qy3,os_fvtxs_qx3)/3.0;
      os_fvtxs_1dPsi2->Fill(os_fvtxs_psi2);
      os_fvtxs_1dPsi3->Fill(os_fvtxs_psi3);

      float os_fvtxn_qw = fvtxn_qw[0];
      float os_fvtxn_qx2 = mean[icent][izvtx][1][8][0];
      float os_fvtxn_qy2 = mean[icent][izvtx][1][8][1];
      float os_fvtxn_qx3 = mean[icent][izvtx][2][8][0];
      float os_fvtxn_qy3 = mean[icent][izvtx][2][8][1];
      float os_fvtxn_qx4 = mean[icent][izvtx][3][8][0];
      float os_fvtxn_qy4 = mean[icent][izvtx][3][8][1];
      os_fvtxn_cos22->Fill(0.0,os_fvtxn_qx2);
      os_fvtxn_sin22->Fill(0.0,os_fvtxn_qy2);
      os_fvtxn_cos32->Fill(0.0,os_fvtxn_qx3);
      os_fvtxn_sin32->Fill(0.0,os_fvtxn_qy3);
      os_fvtxn_qx2 = fvtxn_qx2[0] - os_fvtxn_qw*os_fvtxn_qx2;
      os_fvtxn_qy2 = fvtxn_qy2[0] - os_fvtxn_qw*os_fvtxn_qy2;
      os_fvtxn_qx3 = fvtxn_qx3[0] - os_fvtxn_qw*os_fvtxn_qx3;
      os_fvtxn_qy3 = fvtxn_qy3[0] - os_fvtxn_qw*os_fvtxn_qy3;
      os_fvtxn_qx4 = fvtxn_qx4[0] - os_fvtxn_qw*os_fvtxn_qx4;
      os_fvtxn_qy4 = fvtxn_qy4[0] - os_fvtxn_qw*os_fvtxn_qy4;
      float os_fvtxn_qq2 = ( (os_fvtxn_qx2*os_fvtxn_qx2) + (os_fvtxn_qy2*os_fvtxn_qy2) - os_fvtxn_qw ) / ( (os_fvtxn_qw*os_fvtxn_qw) - os_fvtxn_qw );
      float os_fvtxn_qq3 = ( (os_fvtxn_qx3*os_fvtxn_qx3) + (os_fvtxn_qy3*os_fvtxn_qy3) - os_fvtxn_qw ) / ( (os_fvtxn_qw*os_fvtxn_qw) - os_fvtxn_qw );
      // float os_fvtxn_qq4 = ( (os_fvtxn_qx4*os_fvtxn_qx4) + (os_fvtxn_qy4*os_fvtxn_qy4) - os_fvtxn_qw ) / ( (os_fvtxn_qw*os_fvtxn_qw) - os_fvtxn_qw );
      os_fvtxn_c22->Fill(0.0,os_fvtxn_qq2);
      os_fvtxn_c32->Fill(0.0,os_fvtxn_qq3);
      float os_fvtxn_psi2 = atan2(os_fvtxn_qy2,os_fvtxn_qx2)/2.0;
      float os_fvtxn_psi3 = atan2(os_fvtxn_qy3,os_fvtxn_qx3)/3.0;
      os_fvtxn_1dPsi2->Fill(os_fvtxn_psi2);
      os_fvtxn_1dPsi3->Fill(os_fvtxn_psi3);

      // --- have a look at some different correlations

      // --- bbcs-fvtxs
      float os_bbcfvtxs_qq2 = ( (os_bbc_qx2*os_fvtxs_qx2) + (os_bbc_qy2*os_fvtxs_qy2) ) / ( os_bbc_qw*os_fvtxs_qw );
      float os_bbcfvtxs_qq3 = ( (os_bbc_qx3*os_fvtxs_qx3) + (os_bbc_qy3*os_fvtxs_qy3) ) / ( os_bbc_qw*os_fvtxs_qw );
      os_bbcsfvtxs_c22->Fill(0.0,os_bbcfvtxs_qq2);
      os_bbcsfvtxs_c32->Fill(0.0,os_bbcfvtxs_qq3);
      // --- bbcs-fvtxn
      float os_bbcfvtxn_qq2 = ( (os_bbc_qx2*os_fvtxn_qx2) + (os_bbc_qy2*os_fvtxn_qy2) ) / ( os_bbc_qw*os_fvtxn_qw );
      float os_bbcfvtxn_qq3 = ( (os_bbc_qx3*os_fvtxn_qx3) + (os_bbc_qy3*os_fvtxn_qy3) ) / ( os_bbc_qw*os_fvtxn_qw );
      os_bbcsfvtxn_c22->Fill(0.0,os_bbcfvtxn_qq2);
      os_bbcsfvtxn_c32->Fill(0.0,os_bbcfvtxn_qq3);
      // --- fvtxs-fvtxn
      float os_fvtxsfvtxn_qq2 = ( (os_fvtxs_qx2*os_fvtxn_qx2) + (os_fvtxs_qy2*os_fvtxn_qy2) ) / ( os_fvtxs_qw*os_fvtxn_qw );
      float os_fvtxsfvtxn_qq3 = ( (os_fvtxs_qx3*os_fvtxn_qx3) + (os_fvtxs_qy3*os_fvtxn_qy3) ) / ( os_fvtxs_qw*os_fvtxn_qw );
      os_fvtxsfvtxn_c22->Fill(0.0,os_fvtxsfvtxn_qq2);
      os_fvtxsfvtxn_c32->Fill(0.0,os_fvtxsfvtxn_qq3);

      // --- now have a look at some 4 particle cumulants

      float os_bbc_qqqq4 = calc4_event(os_bbc_qx2,os_bbc_qy2,os_bbc_qx4,os_bbc_qy4,os_bbc_qw);
      float os_fvtxs_qqqq4 = calc4_event(os_fvtxs_qx2,os_fvtxs_qy2,os_fvtxs_qx4,os_fvtxs_qy4,os_fvtxs_qw);
      float os_fvtxn_qqqq4 = calc4_event(os_fvtxn_qx2,os_fvtxn_qy2,os_fvtxn_qx4,os_fvtxn_qy4,os_fvtxn_qw);
      os_bbcs_c24->Fill(0.0,os_bbc_qqqq4);
      os_fvtxs_c24->Fill(0.0,os_fvtxs_qqqq4);
      os_fvtxn_c24->Fill(0.0,os_fvtxn_qqqq4);
      os_bbcs_c24_vs_cent->Fill(centrality,os_bbc_qqqq4);
      os_fvtxs_c24_vs_cent->Fill(centrality,os_fvtxs_qqqq4);
      os_fvtxn_c24_vs_cent->Fill(centrality,os_fvtxn_qqqq4);
      os_bbcs_c24_vs_bbcqs->Fill(bbc_qs,os_bbc_qqqq4);
      os_fvtxs_c24_vs_bbcqs->Fill(bbc_qs,os_fvtxs_qqqq4);
      os_fvtxn_c24_vs_bbcqs->Fill(bbc_qs,os_fvtxn_qqqq4);
      os_bbcs_c24_vs_nfvtxs->Fill(d_nFVTXS_clus,os_bbc_qqqq4);
      os_fvtxs_c24_vs_nfvtxs->Fill(d_nFVTXS_clus,os_fvtxs_qqqq4);
      os_fvtxn_c24_vs_nfvtxs->Fill(d_nFVTXS_clus,os_fvtxn_qqqq4);


      // --------------------------------------
      // --- now the standard event plane stuff
      // --------------------------------------

      // --- looks like these are already done above
      float bbc_psi2 = (sumxy[1][0][2]>0)?sumxy[1][0][3]:-9999.9;
      // --- maybe 12 should be a different number (we have 4 above)
      float fvtx_psi2 = (sumxy[1][1][2]>12)?sumxy[1][1][3]:-9999.9;

      // ---
      // --- south
      // ---

      float bbc_south_psi2_docalib;
      float fvtx_south_psi2_docalib;
      float fvtx0_south_psi2_docalib;
      float fvtx1_south_psi2_docalib;
      float fvtx2_south_psi2_docalib;
      float fvtx3_south_psi2_docalib;

      float bbc_south_psi3_docalib;
      float fvtx_south_psi3_docalib;
      float fvtx0_south_psi3_docalib;
      float fvtx1_south_psi3_docalib;
      float fvtx2_south_psi3_docalib;
      float fvtx3_south_psi3_docalib;

      if ( bbc_pmts )
        {
          bbc_south_psi2_docalib = (sumxy[1][south_bbc_angle][2]>0)?sumxy[1][south_bbc_angle][3]:-9999.9;
          bbc_south_psi3_docalib = (sumxy[2][south_bbc_angle][2]>0)?sumxy[2][south_bbc_angle][3]:-9999.9;
        }
      if ( fvtx_clusters )
        {
          fvtx_south_psi2_docalib = (sumxy[1][south_fvtx_angle][2]>4)?sumxy[1][south_fvtx_angle][3]:-9999.9;
          fvtx0_south_psi2_docalib = (sumxy[1][south_fvtx_0_angle][2]>4)?sumxy[1][south_fvtx_0_angle][3]:-9999.9;
          fvtx1_south_psi2_docalib = (sumxy[1][south_fvtx_1_angle][2]>4)?sumxy[1][south_fvtx_1_angle][3]:-9999.9;
          fvtx2_south_psi2_docalib = (sumxy[1][south_fvtx_2_angle][2]>4)?sumxy[1][south_fvtx_2_angle][3]:-9999.9;
          fvtx3_south_psi2_docalib = (sumxy[1][south_fvtx_3_angle][2]>4)?sumxy[1][south_fvtx_3_angle][3]:-9999.9;
          fvtx_south_psi3_docalib = (sumxy[2][south_fvtx_angle][2]>4)?sumxy[2][south_fvtx_angle][3]:-9999.9;
          fvtx0_south_psi3_docalib = (sumxy[2][south_fvtx_0_angle][2]>4)?sumxy[2][south_fvtx_0_angle][3]:-9999.9;
          fvtx1_south_psi3_docalib = (sumxy[2][south_fvtx_1_angle][2]>4)?sumxy[2][south_fvtx_1_angle][3]:-9999.9;
          fvtx2_south_psi3_docalib = (sumxy[2][south_fvtx_2_angle][2]>4)?sumxy[2][south_fvtx_2_angle][3]:-9999.9;
          fvtx3_south_psi3_docalib = (sumxy[2][south_fvtx_3_angle][2]>4)?sumxy[2][south_fvtx_3_angle][3]:-9999.9;
        }

      // ---
      // --- north
      // ---

      float fvtx_north_psi2_docalib;
      float fvtx0_north_psi2_docalib;
      float fvtx1_north_psi2_docalib;
      float fvtx2_north_psi2_docalib;
      float fvtx3_north_psi2_docalib;

      float fvtx_north_psi3_docalib;
      float fvtx0_north_psi3_docalib;
      float fvtx1_north_psi3_docalib;
      float fvtx2_north_psi3_docalib;
      float fvtx3_north_psi3_docalib;

      if ( fvtx_clusters )
        {
          fvtx_north_psi2_docalib = (sumxy[1][north_fvtx_angle][2]>4)?sumxy[1][north_fvtx_angle][3]:-9999.9;
          fvtx0_north_psi2_docalib = (sumxy[1][north_fvtx_0_angle][2]>4)?sumxy[1][north_fvtx_0_angle][3]:-9999.9;
          fvtx1_north_psi2_docalib = (sumxy[1][north_fvtx_1_angle][2]>4)?sumxy[1][north_fvtx_1_angle][3]:-9999.9;
          fvtx2_north_psi2_docalib = (sumxy[1][north_fvtx_2_angle][2]>4)?sumxy[1][north_fvtx_2_angle][3]:-9999.9;
          fvtx3_north_psi2_docalib = (sumxy[1][north_fvtx_3_angle][2]>4)?sumxy[1][north_fvtx_3_angle][3]:-9999.9;
          fvtx_north_psi3_docalib = (sumxy[2][north_fvtx_angle][2]>4)?sumxy[2][north_fvtx_angle][3]:-9999.9;
          fvtx0_north_psi3_docalib = (sumxy[2][north_fvtx_0_angle][2]>4)?sumxy[2][north_fvtx_0_angle][3]:-9999.9;
          fvtx1_north_psi3_docalib = (sumxy[2][north_fvtx_1_angle][2]>4)?sumxy[2][north_fvtx_1_angle][3]:-9999.9;
          fvtx2_north_psi3_docalib = (sumxy[2][north_fvtx_2_angle][2]>4)?sumxy[2][north_fvtx_2_angle][3]:-9999.9;
          fvtx3_north_psi3_docalib = (sumxy[2][north_fvtx_3_angle][2]>4)?sumxy[2][north_fvtx_3_angle][3]:-9999.9;
        }


      if ( verbosity > 0 && ( fvtx_north_psi2_docalib < -999 || fvtx_south_psi2_docalib < -999 || bbc_south_psi2_docalib < -999 ) )
        {
          cout << "POSSIBLE ISSUE WITH EVENT PLANES!!!  ONE OR MORE IS -9999" << endl;
          cout << "BBC south event plane " << bbc_south_psi2_docalib << endl;
          cout << "FVTX south event plane " << fvtx_south_psi2_docalib << endl;
          cout << "FVTX north event plane " << fvtx_north_psi2_docalib << endl;
          cout << "BBC charge " << bbc_qw << endl;
          cout << "Total number of clusters is " << d_nFVTX_clus << endl;
          cout << "Number of FVTXS clusters " << fvtxs_qw[0] << endl;
          cout << "Number of FVTXN clusters " << fvtxn_qw[0] << endl;
        }



      // ---
      // --- resolution histograms
      // ---

      // --- BBC and FVTX south

      tp1f_reso2_BBC_FVTX->Fill(0.0,cos(2*(bbc_south_psi2_docalib-fvtx_south_psi2_docalib)));
      tp1f_reso3_BBC_FVTX->Fill(0.0,cos(3*(bbc_south_psi3_docalib-fvtx_south_psi3_docalib)));

      th1d_reso2_BBC_FVTX->Fill(cos(2*(bbc_south_psi2_docalib-fvtx_south_psi2_docalib)));
      th1d_reso3_BBC_FVTX->Fill(cos(3*(bbc_south_psi3_docalib-fvtx_south_psi3_docalib)));

      th1d_dreso2_BBC_FVTX->Fill(bbc_south_psi2_docalib-fvtx_south_psi2_docalib);
      th1d_dreso3_BBC_FVTX->Fill(bbc_south_psi3_docalib-fvtx_south_psi3_docalib);

      // --- BBC and FVTX north

      tp1f_reso2_BBC_FVTXN->Fill(0.0,cos(2*(bbc_south_psi2_docalib-fvtx_north_psi2_docalib)));
      tp1f_reso3_BBC_FVTXN->Fill(0.0,cos(3*(bbc_south_psi3_docalib-fvtx_north_psi3_docalib)));

      th1d_reso2_BBC_FVTXN->Fill(cos(2*(bbc_south_psi2_docalib-fvtx_north_psi2_docalib)));
      th1d_reso3_BBC_FVTXN->Fill(cos(3*(bbc_south_psi3_docalib-fvtx_north_psi3_docalib)));

      th1d_dreso2_BBC_FVTXN->Fill(bbc_south_psi2_docalib-fvtx_north_psi2_docalib);
      th1d_dreso3_BBC_FVTXN->Fill(bbc_south_psi3_docalib-fvtx_north_psi3_docalib);

      // --- FVTX south and north

      tp1f_reso2_FVTXS_FVTXN->Fill(0.0,cos(2*(fvtx_south_psi2_docalib-fvtx_north_psi2_docalib)));
      tp1f_reso3_FVTXS_FVTXN->Fill(0.0,cos(3*(fvtx_south_psi3_docalib-fvtx_north_psi3_docalib)));

      th1d_reso2_FVTXS_FVTXN->Fill(cos(2*(fvtx_south_psi2_docalib-fvtx_north_psi2_docalib)));
      th1d_reso3_FVTXS_FVTXN->Fill(cos(3*(fvtx_south_psi3_docalib-fvtx_north_psi3_docalib)));

      th1d_dreso2_FVTXS_FVTXN->Fill(fvtx_south_psi2_docalib-fvtx_north_psi2_docalib);
      th1d_dreso3_FVTXS_FVTXN->Fill(fvtx_south_psi3_docalib-fvtx_north_psi3_docalib);

      // --- for v4psi2

      tp1f_reso42_BBC_FVTX->Fill(0.0,cos(4*(bbc_south_psi2_docalib-fvtx_south_psi2_docalib)));
      tp1f_reso42_BBC_FVTXN->Fill(0.0,cos(4*(bbc_south_psi2_docalib-fvtx_north_psi2_docalib)));
      tp1f_reso42_FVTXS_FVTXN->Fill(0.0,cos(4*(fvtx_south_psi2_docalib-fvtx_north_psi2_docalib)));

      // ---
      // --- resolution histograms
      // ---

      // --- BBC and FVTX south

      tp1f_os_reso2_BBC_FVTX->Fill(0.0,cos(2*(os_bbc_psi2-os_fvtxs_psi2)));
      tp1f_os_reso3_BBC_FVTX->Fill(0.0,cos(3*(os_bbc_psi3-os_fvtxs_psi3)));

      th1d_os_reso2_BBC_FVTX->Fill(cos(2*(os_bbc_psi2-os_fvtxs_psi2)));
      th1d_os_reso3_BBC_FVTX->Fill(cos(3*(os_bbc_psi3-os_fvtxs_psi3)));

      th1d_os_dreso2_BBC_FVTX->Fill(os_bbc_psi2-os_fvtxs_psi2);
      th1d_os_dreso3_BBC_FVTX->Fill(os_bbc_psi3-os_fvtxs_psi3);

      // --- BBC and FVTX north

      tp1f_os_reso2_BBC_FVTXN->Fill(0.0,cos(2*(os_bbc_psi2-os_fvtxn_psi2)));
      tp1f_os_reso3_BBC_FVTXN->Fill(0.0,cos(3*(os_bbc_psi3-os_fvtxn_psi3)));

      th1d_os_reso2_BBC_FVTXN->Fill(cos(2*(os_bbc_psi2-os_fvtxn_psi2)));
      th1d_os_reso3_BBC_FVTXN->Fill(cos(3*(os_bbc_psi3-os_fvtxn_psi3)));

      th1d_os_dreso2_BBC_FVTXN->Fill(os_bbc_psi2-os_fvtxn_psi2);
      th1d_os_dreso3_BBC_FVTXN->Fill(os_bbc_psi3-os_fvtxn_psi3);

      // --- FVTX south and north

      tp1f_os_reso2_FVTXS_FVTXN->Fill(0.0,cos(2*(os_fvtxs_psi2-os_fvtxn_psi2)));
      tp1f_os_reso3_FVTXS_FVTXN->Fill(0.0,cos(3*(os_fvtxs_psi3-os_fvtxn_psi3)));

      th1d_os_reso2_FVTXS_FVTXN->Fill(cos(2*(os_fvtxs_psi2-os_fvtxn_psi2)));
      th1d_os_reso3_FVTXS_FVTXN->Fill(cos(3*(os_fvtxs_psi3-os_fvtxn_psi3)));

      th1d_os_dreso2_FVTXS_FVTXN->Fill(os_fvtxs_psi2-os_fvtxn_psi2);
      th1d_os_dreso3_FVTXS_FVTXN->Fill(os_fvtxs_psi3-os_fvtxn_psi3);

      // ----------------------------------------------------------------------------

      //start of vtx stand alone track loop
      if ( cnt_tracks )
        {
          for(int itrk=0; itrk< d_ntrk; itrk++)
            {
              float px    = d_px[itrk];
              float py    = d_py[itrk];
              float pz    = d_pz[itrk];

              int dcarm=0;
              if(px>0) dcarm=1;

              float phi0 = TMath::ATan2(py,px);
              float pt = sqrt(px*px+py*py);

              if ( pt < 0.2 || pt > 5.0 ) continue; // pt cut added 2016-06-30

              float bbc_dphi2 = phi0 - bbc_psi2;
              if(-4.0<bbc_psi2 && bbc_psi2<4.0 ) // why this weird cut? why not just -pi to pi? // checking against -9999 from above
                {
                  bbcs_v2_incl_nodetree->Fill(pt,cos(2*bbc_dphi2));
                  if(dcarm==0) bbcs_v2_east_nodetree->Fill(pt,cos(2*bbc_dphi2));
                  else if(dcarm==1) bbcs_v2_west_nodetree->Fill(pt,cos(2*bbc_dphi2));
                }

              float fvtx_dphi2 = phi0 - fvtx_psi2;
              if(-4.0<fvtx_psi2 && fvtx_psi2<4.0 )
                {
                  fvtxs_v2_incl_nodetree->Fill(pt,cos(2*fvtx_dphi2));
                  if(dcarm==0) fvtxs_v2_east_nodetree->Fill(pt,cos(2*fvtx_dphi2));
                  else if(dcarm==1) fvtxs_v2_west_nodetree->Fill(pt,cos(2*fvtx_dphi2));
                }

              // -------------------------------------------------------
              // --- finished with nodetree part, now doing docalib part
              // -------------------------------------------------------



              // --- rotation on single particles here
              px = pz*sin(-beam_angle) + px*cos(-beam_angle);
              float phi_angle = TMath::ATan2(py,px);
              float pt_angle = sqrt(px*px+py*py);

              if ( pt_angle < 0.2 || pt_angle > 5.0 ) continue; // pt cut added 2016-06-30

              // ------------------------------------------------------------
              // --- COME BACK HERE
              // --- doing the 2pc up front...
              float ux2 = cos(2*phi_angle);
              float uy2 = sin(2*phi_angle);
              float ux3 = cos(3*phi_angle);
              float uy3 = sin(3*phi_angle);
              float ux4 = cos(4*phi_angle);
              float uy4 = sin(4*phi_angle);

              float bbc_uq2 = ( (ux2*bbc_qx2) + (uy2*bbc_qy2) ) / ( bbc_qw );
              float bbc_uq3 = ( (ux3*bbc_qx3) + (uy3*bbc_qy3) ) / ( bbc_qw );

              bbcs_cos22_both->Fill(pt_angle,ux2);
              if ( dcarm == 1 ) bbcs_cos22_west->Fill(pt_angle,ux2);
              if ( dcarm == 0 ) bbcs_cos22_east->Fill(pt_angle,ux2);
              bbcs_sin22_both->Fill(pt_angle,uy2);
              if ( dcarm == 1 ) bbcs_sin22_west->Fill(pt_angle,uy2);
              if ( dcarm == 0 ) bbcs_sin22_east->Fill(pt_angle,uy2);
              bbcs_d22_both->Fill(pt_angle,bbc_uq2);
              if ( dcarm == 1 ) bbcs_d22_west->Fill(pt_angle,bbc_uq2);
              if ( dcarm == 0 ) bbcs_d22_east->Fill(pt_angle,bbc_uq2);

              bbcs_cos32_both->Fill(pt_angle,ux3);
              if ( dcarm == 1 ) bbcs_cos32_west->Fill(pt_angle,ux3);
              if ( dcarm == 0 ) bbcs_cos32_east->Fill(pt_angle,ux3);
              bbcs_sin32_both->Fill(pt_angle,uy3);
              if ( dcarm == 1 ) bbcs_sin32_west->Fill(pt_angle,uy3);
              if ( dcarm == 0 ) bbcs_sin32_east->Fill(pt_angle,uy3);
              bbcs_d32_both->Fill(pt_angle,bbc_uq3);
              if ( dcarm == 1 ) bbcs_d32_west->Fill(pt_angle,bbc_uq3);
              if ( dcarm == 0 ) bbcs_d32_east->Fill(pt_angle,bbc_uq3);

              // ---

              float fvtxs_uq2 = ( (ux2*fvtxs_qx2[0]) + (uy2*fvtxs_qy2[0]) ) / ( fvtxs_qw[0] );
              float fvtxs_uq3 = ( (ux3*fvtxs_qx3[0]) + (uy3*fvtxs_qy3[0]) ) / ( fvtxs_qw[0] );

              fvtxs_cos22_both->Fill(pt_angle,ux2);
              if ( dcarm == 1 ) fvtxs_cos22_west->Fill(pt_angle,ux2);
              if ( dcarm == 0 ) fvtxs_cos22_east->Fill(pt_angle,ux2);
              fvtxs_sin22_both->Fill(pt_angle,uy2);
              if ( dcarm == 1 ) fvtxs_sin22_west->Fill(pt_angle,uy2);
              if ( dcarm == 0 ) fvtxs_sin22_east->Fill(pt_angle,uy2);
              fvtxs_d22_both->Fill(pt_angle,fvtxs_uq2);
              if ( dcarm == 1 ) fvtxs_d22_west->Fill(pt_angle,fvtxs_uq2);
              if ( dcarm == 0 ) fvtxs_d22_east->Fill(pt_angle,fvtxs_uq2);

              fvtxs_cos32_both->Fill(pt_angle,ux3);
              if ( dcarm == 1 ) fvtxs_cos32_west->Fill(pt_angle,ux3);
              if ( dcarm == 0 ) fvtxs_cos32_east->Fill(pt_angle,ux3);
              fvtxs_sin32_both->Fill(pt_angle,uy3);
              if ( dcarm == 1 ) fvtxs_sin32_west->Fill(pt_angle,uy3);
              if ( dcarm == 0 ) fvtxs_sin32_east->Fill(pt_angle,uy3);
              fvtxs_d32_both->Fill(pt_angle,fvtxs_uq3);
              if ( dcarm == 1 ) fvtxs_d32_west->Fill(pt_angle,fvtxs_uq3);
              if ( dcarm == 0 ) fvtxs_d32_east->Fill(pt_angle,fvtxs_uq3);

              // ------------------------------------------------------------

              float os_bbc_uq2 = ( (ux2*os_bbc_qx2) + (uy2*os_bbc_qy2) ) / ( os_bbc_qw );
              float os_bbc_uq3 = ( (ux3*os_bbc_qx3) + (uy3*os_bbc_qy3) ) / ( os_bbc_qw );

              os_bbcs_cos22_both->Fill(pt_angle,ux2);
              if ( dcarm == 1 ) os_bbcs_cos22_west->Fill(pt_angle,ux2);
              if ( dcarm == 0 ) os_bbcs_cos22_east->Fill(pt_angle,ux2);
              os_bbcs_sin22_both->Fill(pt_angle,uy2);
              if ( dcarm == 1 ) os_bbcs_sin22_west->Fill(pt_angle,uy2);
              if ( dcarm == 0 ) os_bbcs_sin22_east->Fill(pt_angle,uy2);
              os_bbcs_d22_both->Fill(pt_angle,os_bbc_uq2);
              if ( dcarm == 1 ) os_bbcs_d22_west->Fill(pt_angle,os_bbc_uq2);
              if ( dcarm == 0 ) os_bbcs_d22_east->Fill(pt_angle,os_bbc_uq2);

              os_bbcs_cos32_both->Fill(pt_angle,ux3);
              if ( dcarm == 1 ) os_bbcs_cos32_west->Fill(pt_angle,ux3);
              if ( dcarm == 0 ) os_bbcs_cos32_east->Fill(pt_angle,ux3);
              os_bbcs_sin32_both->Fill(pt_angle,uy3);
              if ( dcarm == 1 ) os_bbcs_sin32_west->Fill(pt_angle,uy3);
              if ( dcarm == 0 ) os_bbcs_sin32_east->Fill(pt_angle,uy3);
              os_bbcs_d32_both->Fill(pt_angle,os_bbc_uq3);
              if ( dcarm == 1 ) os_bbcs_d32_west->Fill(pt_angle,os_bbc_uq3);
              if ( dcarm == 0 ) os_bbcs_d32_east->Fill(pt_angle,os_bbc_uq3);

              // ---

              float os_fvtxs_uq2 = ( (ux2*os_fvtxs_qx2) + (uy2*os_fvtxs_qy2) ) / ( os_fvtxs_qw );
              float os_fvtxs_uq3 = ( (ux3*os_fvtxs_qx3) + (uy3*os_fvtxs_qy3) ) / ( os_fvtxs_qw );

              os_fvtxs_cos22_both->Fill(pt_angle,ux2);
              if ( dcarm == 1 ) os_fvtxs_cos22_west->Fill(pt_angle,ux2);
              if ( dcarm == 0 ) os_fvtxs_cos22_east->Fill(pt_angle,ux2);
              os_fvtxs_sin22_both->Fill(pt_angle,uy2);
              if ( dcarm == 1 ) os_fvtxs_sin22_west->Fill(pt_angle,uy2);
              if ( dcarm == 0 ) os_fvtxs_sin22_east->Fill(pt_angle,uy2);
              os_fvtxs_d22_both->Fill(pt_angle,os_fvtxs_uq2);
              if ( dcarm == 1 ) os_fvtxs_d22_west->Fill(pt_angle,os_fvtxs_uq2);
              if ( dcarm == 0 ) os_fvtxs_d22_east->Fill(pt_angle,os_fvtxs_uq2);

              os_fvtxs_cos32_both->Fill(pt_angle,ux3);
              if ( dcarm == 1 ) os_fvtxs_cos32_west->Fill(pt_angle,ux3);
              if ( dcarm == 0 ) os_fvtxs_cos32_east->Fill(pt_angle,ux3);
              os_fvtxs_sin32_both->Fill(pt_angle,uy3);
              if ( dcarm == 1 ) os_fvtxs_sin32_west->Fill(pt_angle,uy3);
              if ( dcarm == 0 ) os_fvtxs_sin32_east->Fill(pt_angle,uy3);
              os_fvtxs_d32_both->Fill(pt_angle,os_fvtxs_uq3);
              if ( dcarm == 1 ) os_fvtxs_d32_west->Fill(pt_angle,os_fvtxs_uq3);
              if ( dcarm == 0 ) os_fvtxs_d32_east->Fill(pt_angle,os_fvtxs_uq3);

              // ---

              float os_fvtxn_uq2 = ( (ux2*os_fvtxn_qx2) + (uy2*os_fvtxn_qy2) ) / ( os_fvtxn_qw );
              float os_fvtxn_uq3 = ( (ux3*os_fvtxn_qx3) + (uy3*os_fvtxn_qy3) ) / ( os_fvtxn_qw );

              os_fvtxn_cos22_both->Fill(pt_angle,ux2);
              if ( dcarm == 1 ) os_fvtxn_cos22_west->Fill(pt_angle,ux2);
              if ( dcarm == 0 ) os_fvtxn_cos22_east->Fill(pt_angle,ux2);
              os_fvtxn_sin22_both->Fill(pt_angle,uy2);
              if ( dcarm == 1 ) os_fvtxn_sin22_west->Fill(pt_angle,uy2);
              if ( dcarm == 0 ) os_fvtxn_sin22_east->Fill(pt_angle,uy2);
              os_fvtxn_d22_both->Fill(pt_angle,os_fvtxn_uq2);
              if ( dcarm == 1 ) os_fvtxn_d22_west->Fill(pt_angle,os_fvtxn_uq2);
              if ( dcarm == 0 ) os_fvtxn_d22_east->Fill(pt_angle,os_fvtxn_uq2);

              os_fvtxn_cos32_both->Fill(pt_angle,ux3);
              if ( dcarm == 1 ) os_fvtxn_cos32_west->Fill(pt_angle,ux3);
              if ( dcarm == 0 ) os_fvtxn_cos32_east->Fill(pt_angle,ux3);
              os_fvtxn_sin32_both->Fill(pt_angle,uy3);
              if ( dcarm == 1 ) os_fvtxn_sin32_west->Fill(pt_angle,uy3);
              if ( dcarm == 0 ) os_fvtxn_sin32_east->Fill(pt_angle,uy3);
              os_fvtxn_d32_both->Fill(pt_angle,os_fvtxn_uq3);
              if ( dcarm == 1 ) os_fvtxn_d32_west->Fill(pt_angle,os_fvtxn_uq3);
              if ( dcarm == 0 ) os_fvtxn_d32_east->Fill(pt_angle,os_fvtxn_uq3);

              // ------------------------------------------------------------

              os_bbcs_v2_both->Fill(pt_angle,cos(2.0*(phi_angle-os_bbc_psi2)));
              if ( dcarm == 1 ) os_bbcs_v2_west->Fill(pt_angle,cos(2.0*(phi_angle-os_bbc_psi2)));
              if ( dcarm == 0 ) os_bbcs_v2_east->Fill(pt_angle,cos(2.0*(phi_angle-os_bbc_psi2)));
              os_fvtxs_v2_both->Fill(pt_angle,cos(2.0*(phi_angle-os_fvtxs_psi2)));
              if ( dcarm == 1 ) os_fvtxs_v2_west->Fill(pt_angle,cos(2.0*(phi_angle-os_fvtxs_psi2)));
              if ( dcarm == 0 ) os_fvtxs_v2_east->Fill(pt_angle,cos(2.0*(phi_angle-os_fvtxs_psi2)));
              os_fvtxn_v2_both->Fill(pt_angle,cos(2.0*(phi_angle-os_fvtxn_psi2)));
              if ( dcarm == 1 ) os_fvtxn_v2_west->Fill(pt_angle,cos(2.0*(phi_angle-os_fvtxn_psi2)));
              if ( dcarm == 0 ) os_fvtxn_v2_east->Fill(pt_angle,cos(2.0*(phi_angle-os_fvtxn_psi2)));

              os_bbcs_v3_both->Fill(pt_angle,cos(3.0*(phi_angle-os_bbc_psi3)));
              if ( dcarm == 1 ) os_bbcs_v3_west->Fill(pt_angle,cos(3.0*(phi_angle-os_bbc_psi3)));
              if ( dcarm == 0 ) os_bbcs_v3_east->Fill(pt_angle,cos(3.0*(phi_angle-os_bbc_psi3)));
              os_fvtxs_v3_both->Fill(pt_angle,cos(3.0*(phi_angle-os_fvtxs_psi3)));
              if ( dcarm == 1 ) os_fvtxs_v3_west->Fill(pt_angle,cos(3.0*(phi_angle-os_fvtxs_psi3)));
              if ( dcarm == 0 ) os_fvtxs_v3_east->Fill(pt_angle,cos(3.0*(phi_angle-os_fvtxs_psi3)));
              os_fvtxn_v3_both->Fill(pt_angle,cos(3.0*(phi_angle-os_fvtxn_psi3)));
              if ( dcarm == 1 ) os_fvtxn_v3_west->Fill(pt_angle,cos(3.0*(phi_angle-os_fvtxn_psi3)));
              if ( dcarm == 0 ) os_fvtxn_v3_east->Fill(pt_angle,cos(3.0*(phi_angle-os_fvtxn_psi3)));

              // --- 4 particle stuff

              // --- POI is in RP because we forcibly add it
              float bbc_d24_in = calc4_track_flag(ux2,uy2,ux4,uy4,ux2+os_bbc_qx2,uy2+os_bbc_qy2,ux4+os_bbc_qx4,uy4+os_bbc_qy4,1+os_bbc_qw,true);
              float fvtxs_d24_in = calc4_track_flag(ux2,uy2,ux4,uy4,ux2+os_fvtxs_qx2,uy2+os_fvtxs_qy2,ux4+os_fvtxs_qx4,uy4+os_fvtxs_qy4,1+os_fvtxs_qw,true);
              float fvtxn_d24_in = calc4_track_flag(ux2,uy2,ux4,uy4,ux2+os_fvtxn_qx2,uy2+os_fvtxn_qy2,ux4+os_fvtxn_qx4,uy4+os_fvtxn_qy4,1+os_fvtxn_qw,true);
              os_bbcs_d24_in_both->Fill(pt_angle,bbc_d24_in);
              os_fvtxs_d24_in_both->Fill(pt_angle,fvtxs_d24_in);
              os_fvtxn_d24_in_both->Fill(pt_angle,fvtxn_d24_in);
              // --- POI is not in RP, this is the natural situation because the POI and RP are in different detectors
              float bbc_d24_out = calc4_track_flag(ux2,uy2,ux4,uy4,os_bbc_qx2,os_bbc_qy2,os_bbc_qx4,os_bbc_qy4,os_bbc_qw,false);
              float fvtxs_d24_out = calc4_track_flag(ux2,uy2,ux4,uy4,os_fvtxs_qx2,os_fvtxs_qy2,os_fvtxs_qx4,os_fvtxs_qy4,os_fvtxs_qw,false);
              float fvtxn_d24_out = calc4_track_flag(ux2,uy2,ux4,uy4,os_fvtxn_qx2,os_fvtxn_qy2,os_fvtxn_qx4,os_fvtxn_qy4,os_fvtxn_qw,false);
              os_bbcs_d24_out_both->Fill(pt_angle,bbc_d24_out);
              os_fvtxs_d24_out_both->Fill(pt_angle,fvtxs_d24_out);
              os_fvtxn_d24_out_both->Fill(pt_angle,fvtxn_d24_out);


              //bbc angle
              if ( bbc_pmts )
                {
                  if(-4.0<bbc_south_psi2_docalib && bbc_south_psi2_docalib<4.0)
                    {
                      // --- 2nd harmonic
                      double bbc_dphi2_docalib = phi_angle - bbc_south_psi2_docalib;
                      double cosbbc_dphi2_docalib = TMath::Cos(2*bbc_dphi2_docalib);
                      bbcs_v2_both_docalib->Fill(pt_angle,cosbbc_dphi2_docalib);
                      if ( dcarm == 1 ) bbcs_v2_west_docalib->Fill(pt_angle,cosbbc_dphi2_docalib);
                      if ( dcarm == 0 ) bbcs_v2_east_docalib->Fill(pt_angle,cosbbc_dphi2_docalib);
                      // --- 3rd harmonic
                      double bbc_dphi3_docalib = phi_angle - bbc_south_psi3_docalib;
                      double cosbbc_dphi3_docalib = TMath::Cos(3*bbc_dphi3_docalib);
                      bbcs_v3_both_docalib->Fill(pt_angle,cosbbc_dphi3_docalib);
                      if ( dcarm == 1 ) bbcs_v3_west_docalib->Fill(pt_angle,cosbbc_dphi3_docalib);
                      if ( dcarm == 0 ) bbcs_v3_east_docalib->Fill(pt_angle,cosbbc_dphi3_docalib);
                      // --- 4th harmonic
                      double fourfour = 4*phi_angle - 4*bbc_south_psi2_docalib;
                      double cosfourfour = TMath::Cos(fourfour);
                      bbcs_v4_4Psi2_both_docalib->Fill(pt_angle,cosfourfour);
                      if ( dcarm == 1 ) bbcs_v4_4Psi2_west_docalib->Fill(pt_angle,cosfourfour);
                      if ( dcarm == 0 ) bbcs_v4_4Psi2_east_docalib->Fill(pt_angle,cosfourfour);
                      // --- ep reso
                      tp1f_reso2_BBC_CNT->Fill(0.0,cosbbc_dphi2_docalib);
                      tp1f_reso3_BBC_CNT->Fill(0.0,cosbbc_dphi3_docalib);
                      th1d_reso2_BBC_CNT->Fill(cosbbc_dphi2_docalib);
                      th1d_reso3_BBC_CNT->Fill(cosbbc_dphi3_docalib);
                      th1d_dreso2_BBC_CNT->Fill(bbc_dphi2_docalib);
                      th1d_dreso3_BBC_CNT->Fill(bbc_dphi3_docalib);
                      // ---
                      tp1f_reso42_BBC_CNT->Fill(0.0,cosfourfour);
                      // --- ep reso with offset q-vectors
                      tp1f_os_reso2_BBC_CNT->Fill(0.0, cos(2.0*(phi_angle-os_bbc_psi2)) );
                      tp1f_os_reso3_BBC_CNT->Fill(0.0, cos(3.0*(phi_angle-os_bbc_psi3)) );
                      th1d_os_reso2_BBC_CNT->Fill( cos(2.0*(phi_angle-os_bbc_psi2)) );
                      th1d_os_reso3_BBC_CNT->Fill( cos(3.0*(phi_angle-os_bbc_psi3)) );
                      th1d_os_dreso2_BBC_CNT->Fill( phi_angle-os_bbc_psi2 );
                      th1d_os_dreso3_BBC_CNT->Fill( phi_angle-os_bbc_psi3 );
                    }
                } // check on tubes

              if(!fvtx_clusters) continue;
              //fvtx all layers
              if(-4.0<fvtx_south_psi2_docalib && fvtx_south_psi2_docalib<4.0)
                {
                  // --- 2nd harmonic
                  double fvtx_dphi2_docalib = phi_angle - fvtx_south_psi2_docalib;
                  double cosfvtx_dphi2_docalib = TMath::Cos(2*fvtx_dphi2_docalib);
                  fvtxs_v2_both_docalib->Fill(pt_angle,cosfvtx_dphi2_docalib);
                  if ( dcarm == 1 ) fvtxs_v2_west_docalib->Fill(pt_angle,cosfvtx_dphi2_docalib);
                  if ( dcarm == 0 ) fvtxs_v2_east_docalib->Fill(pt_angle,cosfvtx_dphi2_docalib);
                  // --- 3rd harmonic
                  double fvtx_dphi3_docalib = phi_angle - fvtx_south_psi3_docalib;
                  double cosfvtx_dphi3_docalib = TMath::Cos(3*fvtx_dphi3_docalib);
                  fvtxs_v3_both_docalib->Fill(pt_angle,cosfvtx_dphi3_docalib);
                  if ( dcarm == 1 ) fvtxs_v3_west_docalib->Fill(pt_angle,cosfvtx_dphi3_docalib);
                  if ( dcarm == 0 ) fvtxs_v3_east_docalib->Fill(pt_angle,cosfvtx_dphi3_docalib);
                  // --- 4th harmonic
                  double fourfour = 4*phi_angle - 4*fvtx_south_psi2_docalib;
                  double cosfourfour = TMath::Cos(fourfour);
                  fvtxs_v4_4Psi2_both_docalib->Fill(pt_angle,cosfourfour);
                  if ( dcarm == 1 ) fvtxs_v4_4Psi2_west_docalib->Fill(pt_angle,cosfourfour);
                  if ( dcarm == 0 ) fvtxs_v4_4Psi2_east_docalib->Fill(pt_angle,cosfourfour);
                  // --- ep reso
                  tp1f_reso2_CNT_FVTX->Fill(0.0,cosfvtx_dphi2_docalib);
                  tp1f_reso3_CNT_FVTX->Fill(0.0,cosfvtx_dphi3_docalib);
                  th1d_reso2_CNT_FVTX->Fill(cosfvtx_dphi2_docalib);
                  th1d_reso3_CNT_FVTX->Fill(cosfvtx_dphi3_docalib);
                  th1d_dreso2_CNT_FVTX->Fill(fvtx_dphi2_docalib);
                  th1d_dreso3_CNT_FVTX->Fill(fvtx_dphi3_docalib);
                  // ---
                  tp1f_reso42_CNT_FVTX->Fill(0.0,cosfourfour);
                  // --- ep reso with offset q-vectors
                  tp1f_os_reso2_CNT_FVTX->Fill(0.0, cos(2.0*(phi_angle-os_fvtxs_psi2)) );
                  tp1f_os_reso3_CNT_FVTX->Fill(0.0, cos(3.0*(phi_angle-os_fvtxs_psi3)) );
                  th1d_os_reso2_CNT_FVTX->Fill( cos(2.0*(phi_angle-os_fvtxs_psi2)) );
                  th1d_os_reso3_CNT_FVTX->Fill( cos(3.0*(phi_angle-os_fvtxs_psi3)) );
                  th1d_os_dreso2_CNT_FVTX->Fill( phi_angle-os_fvtxs_psi2 );
                  th1d_os_dreso3_CNT_FVTX->Fill( phi_angle-os_fvtxs_psi3 );
                }

              // --- North

              //fvtx all layers
              if(-4.0<fvtx_north_psi2_docalib && fvtx_north_psi2_docalib<4.0)
                {
                  // --- 2nd harmonic
                  double fvtx_dphi2_docalib = phi_angle - fvtx_north_psi2_docalib;
                  double cosfvtx_dphi2_docalib = TMath::Cos(2*fvtx_dphi2_docalib);
                  fvtxn_v2_both_docalib->Fill(pt_angle,cosfvtx_dphi2_docalib);
                  if ( dcarm == 1 ) fvtxn_v2_west_docalib->Fill(pt_angle,cosfvtx_dphi2_docalib);
                  if ( dcarm == 0 ) fvtxn_v2_east_docalib->Fill(pt_angle,cosfvtx_dphi2_docalib);
                  // --- 3rd harmonic
                  double fvtx_dphi3_docalib = phi_angle - fvtx_north_psi3_docalib;
                  double cosfvtx_dphi3_docalib = TMath::Cos(3*fvtx_dphi3_docalib);
                  fvtxn_v3_both_docalib->Fill(pt_angle,cosfvtx_dphi3_docalib);
                  if ( dcarm == 1 ) fvtxn_v3_west_docalib->Fill(pt_angle,cosfvtx_dphi3_docalib);
                  if ( dcarm == 0 ) fvtxn_v3_east_docalib->Fill(pt_angle,cosfvtx_dphi3_docalib);
                  // --- 4th harmonic
                  double fourfour = 4*phi_angle - 4*fvtx_north_psi2_docalib;
                  double cosfourfour = TMath::Cos(fourfour);
                  fvtxn_v4_4Psi2_both_docalib->Fill(pt_angle,cosfourfour);
                  if ( dcarm == 1 ) fvtxn_v4_4Psi2_west_docalib->Fill(pt_angle,cosfourfour);
                  if ( dcarm == 0 ) fvtxn_v4_4Psi2_east_docalib->Fill(pt_angle,cosfourfour);
                  // --- ep reso
                  tp1f_reso2_CNT_FVTXN->Fill(0.0,cosfvtx_dphi2_docalib);
                  tp1f_reso3_CNT_FVTXN->Fill(0.0,cosfvtx_dphi3_docalib);
                  th1d_reso2_CNT_FVTXN->Fill(cosfvtx_dphi2_docalib);
                  th1d_reso3_CNT_FVTXN->Fill(cosfvtx_dphi3_docalib);
                  th1d_dreso2_CNT_FVTXN->Fill(fvtx_dphi2_docalib);
                  th1d_dreso3_CNT_FVTXN->Fill(fvtx_dphi3_docalib);
                  // ---
                  tp1f_reso42_CNT_FVTXN->Fill(0.0,cosfourfour);
                  // --- ep reso with offset q-vectors
                  tp1f_os_reso2_CNT_FVTXN->Fill(0.0, cos(2.0*(phi_angle-os_fvtxn_psi2)) );
                  tp1f_os_reso3_CNT_FVTXN->Fill(0.0, cos(3.0*(phi_angle-os_fvtxn_psi3)) );
                  th1d_os_reso2_CNT_FVTXN->Fill( cos(2.0*(phi_angle-os_fvtxn_psi2)) );
                  th1d_os_reso3_CNT_FVTXN->Fill( cos(3.0*(phi_angle-os_fvtxn_psi3)) );
                  th1d_os_dreso2_CNT_FVTXN->Fill( phi_angle-os_fvtxn_psi2 );
                  th1d_os_dreso3_CNT_FVTXN->Fill( phi_angle-os_fvtxn_psi3 );
                }

              // --- now fvtx layers

              //fvtx layer 0
              if(-4.0<fvtx0_south_psi2_docalib && fvtx0_south_psi2_docalib<4.0)
                {
                  double fvtx_dphi2_docalib = phi_angle - fvtx0_south_psi2_docalib;
                  double cosfvtx_dphi2_docalib = TMath::Cos(2*fvtx_dphi2_docalib);
                  if ( dcarm == 1 ) fvtxs0_v2_west_docalib->Fill(pt_angle,cosfvtx_dphi2_docalib);
                  if ( dcarm == 0 ) fvtxs0_v2_east_docalib->Fill(pt_angle,cosfvtx_dphi2_docalib);
                  fvtxs0_v2_both_docalib->Fill(pt_angle,cosfvtx_dphi2_docalib);
                }

              //fvtx layer 1
              if(-4.0<fvtx1_south_psi2_docalib && fvtx1_south_psi2_docalib<4.0)
                {
                  double fvtx_dphi2_docalib = phi_angle - fvtx1_south_psi2_docalib;
                  double cosfvtx_dphi2_docalib = TMath::Cos(2*fvtx_dphi2_docalib);
                  if ( dcarm == 1 ) fvtxs1_v2_west_docalib->Fill(pt_angle,cosfvtx_dphi2_docalib);
                  if ( dcarm == 0 ) fvtxs1_v2_east_docalib->Fill(pt_angle,cosfvtx_dphi2_docalib);
                  fvtxs1_v2_both_docalib->Fill(pt_angle,cosfvtx_dphi2_docalib);
                }

              //fvtx layer 2
              if(-4.0<fvtx2_south_psi2_docalib && fvtx2_south_psi2_docalib<4.0)
                {
                  double fvtx_dphi2_docalib = phi_angle - fvtx2_south_psi2_docalib;
                  double cosfvtx_dphi2_docalib = TMath::Cos(2*fvtx_dphi2_docalib);
                  if ( dcarm == 1 ) fvtxs2_v2_west_docalib->Fill(pt_angle,cosfvtx_dphi2_docalib);
                  if ( dcarm == 0 ) fvtxs2_v2_east_docalib->Fill(pt_angle,cosfvtx_dphi2_docalib);
                  fvtxs2_v2_both_docalib->Fill(pt_angle,cosfvtx_dphi2_docalib);
                }

              //fvtx layer 3
              if(-4.0<fvtx3_south_psi2_docalib && fvtx3_south_psi2_docalib<4.0)
                {
                  double fvtx_dphi2_docalib = phi_angle - fvtx3_south_psi2_docalib;
                  double cosfvtx_dphi2_docalib = TMath::Cos(2*fvtx_dphi2_docalib);
                  if ( dcarm == 1 ) fvtxs3_v2_west_docalib->Fill(pt_angle,cosfvtx_dphi2_docalib);
                  if ( dcarm == 0 ) fvtxs3_v2_east_docalib->Fill(pt_angle,cosfvtx_dphi2_docalib);
                  fvtxs3_v2_both_docalib->Fill(pt_angle,cosfvtx_dphi2_docalib);
                } // check on psi
            } // loop over tracks
        } // check on tracks
    }//end of event

  cout << "Processed " << event_counter << "/" << all_counter << " events (" << (float)event_counter/(float)all_counter << ")" << endl;
  cout << "Events with bad centrality = " << bad_cent_counter << " (" << (float)bad_cent_counter/(float)event_counter << ")" << endl;
  cout << "Events with vertex disagreement = " << bad_vertex_counter << " (" << (float)bad_vertex_counter/(float)event_counter << ")" << endl;
  cout << "Gap cuts applied " << gapcut_counter << "/" << cluster_counter << " (" << (float)gapcut_counter/(float)cluster_counter << ")" << endl;

  if(rp_recal_pass<3 && rp_recal_pass>0)
    {
      // --- previous pass calib file is named above, rename it here
      sprintf(calibfile,"output/flattening_%d_%d.dat",runNumber,rp_recal_pass);
      cout << "writing calibration file : " << calibfile << endl;
      ofstream ofs;
      ofs.open(calibfile);

      cout << "writing out flattening parameters" << endl;
      //int icent = 0;
      for (int iz=0; iz<NZPS; iz++)
        {
          for (int ih=1; ih<NHAR; ih++)
            {
              for (int id=0; id<NDET; id++)
                {
                  for (int ib=0; ib<2; ib++)
                    {
                      if ( DIAG ) cout<<"writing ave:  "<<ic<<" "<<iz<<" "<<ih<<" "<<id<<endl;
                      // write out average qx, qx error, qy, qy error
                      ofs << ave[icent][iz][ih][id]->GetBinContent(ib+1) << " ";
                      ofs << ave[icent][iz][ih][id]->GetBinError  (ib+1) << " ";
                    } // x and y
                  ofs << endl;
                  for (int ib=0; ib<2; ib++)
                    {
                      for (int io=0; io<NORD; io++)
                        {
                          // write first 12 orders for fourier fit of psi
                          // we are unsure of what's being written out here
                          ofs << flt[icent][iz][ih][id]->GetBinContent(ib*NORD+io+1) << " ";
                        } // orders
                      ofs << endl;
                    } // x and y
                } // detectors
            } // harmonics
        } // z_vertex bins
      ofs.close();
    } // if pass 1 or 2

  // --- QA purpose
  if(rp_recal_pass>0)
    {
      TFile *mData2=TFile::Open(outFile2,"recreate");
      mData2->cd();
      for (int iz=0; iz<NZPS; iz++)
        {
          for (int ih=1; ih<NHAR; ih++)
            {
              for (int id=0; id<NDET; id++)
                {
                  ave[ic][iz][ih][id]->Write();
                  flt[ic][iz][ih][id]->Write();
                } // detectors
            } // harmonics
        } // z_vertex bins

      for (int ih=1; ih<NHAR; ih++)
        {
          for (int id=0; id<NDET; id++)
            {
              qx[ic][ih][id]->Write();
              qy[ic][ih][id]->Write();
              dis[ic][ih][id]->Write();
              psi_bf[ic][ih][id]->Write();
              psi_af[ic][ih][id]->Write();
            } // detectors
        } // harmonics
      os_bbcs_1dPsi2->Write();
      os_bbcs_1dPsi3->Write();
      os_fvtxs_1dPsi2->Write();
      os_fvtxs_1dPsi3->Write();

      th1d_BBC_charge->Write();
      th1d_FVTX_nclus->Write();
      th1d_FVTXS_nclus->Write();
      th1d_FVTXN_nclus->Write();
      th2d_qBBC_nFVTX->Write();

      th1d_bbc_charge_phi->Write();
      th1d_bbc0_charge_phi->Write();
      th1d_bbc1_charge_phi->Write();
      th1d_bbc2_charge_phi->Write();
      th1d_bbc3_charge_phi->Write();
      th1d_bbc4_charge_phi->Write();

      th1d_fvtxs_clus_phi->Write();
      th1d_fvtxs0_clus_phi->Write();
      th1d_fvtxs1_clus_phi->Write();
      th1d_fvtxs2_clus_phi->Write();
      th1d_fvtxs3_clus_phi->Write();

      th1d_fvtxn_clus_phi->Write();
      th1d_fvtxn0_clus_phi->Write();
      th1d_fvtxn1_clus_phi->Write();
      th1d_fvtxn2_clus_phi->Write();
      th1d_fvtxn3_clus_phi->Write();

      mData2->Close();
    }

  // --- write v2 histograms if on last pass
  if(rp_recal_pass>2)
    {
      TFile *mData1=TFile::Open(outFile1,"recreate");
      mData1->cd();

      bbcs_v2_west_nodetree->Write();
      bbcs_v2_east_nodetree->Write();
      bbcs_v2_incl_nodetree->Write();

      fvtxs_v2_west_nodetree->Write();
      fvtxs_v2_east_nodetree->Write();
      fvtxs_v2_incl_nodetree->Write();

      if ( bbc_pmts )
        {
          bbcs_v2_east_docalib->Write();
          bbcs_v2_west_docalib->Write();
          bbcs_v2_both_docalib->Write();
          bbcs_v3_east_docalib->Write();
          bbcs_v3_west_docalib->Write();
          bbcs_v3_both_docalib->Write();
          bbcs_v4_4Psi2_east_docalib->Write();
          bbcs_v4_4Psi2_west_docalib->Write();
          bbcs_v4_4Psi2_both_docalib->Write();
        }
      if ( fvtx_clusters )
        {
          fvtxs_v2_east_docalib->Write();
          fvtxs_v2_west_docalib->Write();
          fvtxs_v2_both_docalib->Write();
          // ---
          fvtxs_v3_east_docalib->Write();
          fvtxs_v3_west_docalib->Write();
          fvtxs_v3_both_docalib->Write();
          // ---
          fvtxs0_v2_east_docalib->Write();
          fvtxs1_v2_east_docalib->Write();
          fvtxs2_v2_east_docalib->Write();
          fvtxs3_v2_east_docalib->Write();
          // ---
          fvtxs0_v2_west_docalib->Write();
          fvtxs1_v2_west_docalib->Write();
          fvtxs2_v2_west_docalib->Write();
          fvtxs3_v2_west_docalib->Write();
          // ---
          fvtxs0_v2_both_docalib->Write();
          fvtxs1_v2_both_docalib->Write();
          fvtxs2_v2_both_docalib->Write();
          fvtxs3_v2_both_docalib->Write();
          // ---
          fvtxs_v4_4Psi2_east_docalib->Write();
          fvtxs_v4_4Psi2_west_docalib->Write();
          fvtxs_v4_4Psi2_both_docalib->Write();
          // ---
          fvtxn_v2_east_docalib->Write();
          fvtxn_v2_west_docalib->Write();
          fvtxn_v2_both_docalib->Write();
          fvtxn_v3_east_docalib->Write();
          fvtxn_v3_west_docalib->Write();
          fvtxn_v3_both_docalib->Write();
          fvtxn_v4_4Psi2_east_docalib->Write();
          fvtxn_v4_4Psi2_west_docalib->Write();
          fvtxn_v4_4Psi2_both_docalib->Write();
        } // fvtx clsuters

      tp1f_reso2_BBC_CNT->Write();
      tp1f_reso2_BBC_FVTX->Write();
      tp1f_reso2_CNT_FVTX->Write();
      tp1f_reso3_BBC_CNT->Write();
      tp1f_reso3_BBC_FVTX->Write();
      tp1f_reso3_CNT_FVTX->Write();

      th1d_reso2_BBC_CNT->Write();
      th1d_reso2_BBC_FVTX->Write();
      th1d_reso2_CNT_FVTX->Write();
      th1d_reso3_BBC_CNT->Write();
      th1d_reso3_BBC_FVTX->Write();
      th1d_reso3_CNT_FVTX->Write();

      th1d_dreso2_BBC_CNT->Write();
      th1d_dreso2_BBC_FVTX->Write();
      th1d_dreso2_CNT_FVTX->Write();
      th1d_dreso3_BBC_CNT->Write();
      th1d_dreso3_BBC_FVTX->Write();
      th1d_dreso3_CNT_FVTX->Write();

      tp1f_reso2_BBC_FVTXN->Write();
      tp1f_reso3_BBC_FVTXN->Write();
      tp1f_reso2_CNT_FVTXN->Write();
      tp1f_reso3_CNT_FVTXN->Write();
      tp1f_reso2_FVTXS_FVTXN->Write();
      tp1f_reso3_FVTXS_FVTXN->Write();

      th1d_reso2_BBC_FVTXN->Write();
      th1d_reso3_BBC_FVTXN->Write();
      th1d_reso2_CNT_FVTXN->Write();
      th1d_reso3_CNT_FVTXN->Write();
      th1d_reso2_FVTXS_FVTXN->Write();
      th1d_reso3_FVTXS_FVTXN->Write();

      th1d_dreso2_BBC_FVTXN->Write();
      th1d_dreso3_BBC_FVTXN->Write();
      th1d_dreso2_CNT_FVTXN->Write();
      th1d_dreso3_CNT_FVTXN->Write();
      th1d_dreso2_FVTXS_FVTXN->Write();
      th1d_dreso3_FVTXS_FVTXN->Write();

      tp1f_reso42_BBC_CNT->Write();
      tp1f_reso42_BBC_FVTX->Write();
      tp1f_reso42_CNT_FVTX->Write();
      tp1f_reso42_BBC_FVTXN->Write();
      tp1f_reso42_CNT_FVTXN->Write();
      tp1f_reso42_FVTXS_FVTXN->Write();

      th1d_BBC_charge->Write();
      th1d_FVTX_nclus->Write();
      th1d_FVTXS_nclus->Write();
      th1d_FVTXN_nclus->Write();
      th2d_qBBC_nFVTX->Write();

      tp1f_reso2_CNT_FVTXN->Write();
      tp1f_os_reso2_BBC_CNT->Write();
      tp1f_os_reso2_BBC_FVTX->Write();
      tp1f_os_reso2_CNT_FVTX->Write();
      tp1f_os_reso3_BBC_CNT->Write();
      tp1f_os_reso3_BBC_FVTX->Write();
      tp1f_os_reso3_CNT_FVTX->Write();
      th1d_os_reso2_BBC_CNT->Write();
      th1d_os_reso2_BBC_FVTX->Write();
      th1d_os_reso2_CNT_FVTX->Write();
      th1d_os_reso3_BBC_CNT->Write();
      th1d_os_reso3_BBC_FVTX->Write();
      th1d_os_reso3_CNT_FVTX->Write();
      th1d_os_dreso2_BBC_CNT->Write();
      th1d_os_dreso2_BBC_FVTX->Write();
      th1d_os_dreso2_CNT_FVTX->Write();
      th1d_os_dreso3_BBC_CNT->Write();
      th1d_os_dreso3_BBC_FVTX->Write();
      th1d_os_dreso3_CNT_FVTX->Write();
      tp1f_os_reso2_BBC_FVTXN->Write();
      tp1f_os_reso3_BBC_FVTXN->Write();
      tp1f_os_reso2_CNT_FVTXN->Write();
      tp1f_os_reso3_CNT_FVTXN->Write();
      tp1f_os_reso2_FVTXS_FVTXN->Write();
      tp1f_os_reso3_FVTXS_FVTXN->Write();
      th1d_os_reso2_BBC_FVTXN->Write();
      th1d_os_reso3_BBC_FVTXN->Write();
      th1d_os_reso2_CNT_FVTXN->Write();
      th1d_os_reso3_CNT_FVTXN->Write();
      th1d_os_reso2_FVTXS_FVTXN->Write();
      th1d_os_reso3_FVTXS_FVTXN->Write();
      th1d_os_dreso2_BBC_FVTXN->Write();
      th1d_os_dreso3_BBC_FVTXN->Write();
      th1d_os_dreso2_CNT_FVTXN->Write();
      th1d_os_dreso3_CNT_FVTXN->Write();
      th1d_os_dreso2_FVTXS_FVTXN->Write();
      th1d_os_dreso3_FVTXS_FVTXN->Write();

      th1d_bbc_charge_phi->Write();
      th1d_bbc0_charge_phi->Write();
      th1d_bbc1_charge_phi->Write();
      th1d_bbc2_charge_phi->Write();
      th1d_bbc3_charge_phi->Write();
      th1d_bbc4_charge_phi->Write();

      th1d_fvtxs_clus_phi->Write();
      th1d_fvtxs0_clus_phi->Write();
      th1d_fvtxs1_clus_phi->Write();
      th1d_fvtxs2_clus_phi->Write();
      th1d_fvtxs3_clus_phi->Write();

      th1d_fvtxn_clus_phi->Write();
      th1d_fvtxn0_clus_phi->Write();
      th1d_fvtxn1_clus_phi->Write();
      th1d_fvtxn2_clus_phi->Write();
      th1d_fvtxn3_clus_phi->Write();

      // ---

      bbcs_d22_west->Write();
      bbcs_d22_east->Write();
      bbcs_d22_both->Write();

      fvtxs_d22_west->Write();
      fvtxs_d22_east->Write();
      fvtxs_d22_both->Write();

      // fvtxn_d22_west->Write();
      // fvtxn_d22_east->Write();
      // fvtxn_d22_both->Write();

      bbcs_d32_west->Write();
      bbcs_d32_east->Write();
      bbcs_d32_both->Write();

      fvtxs_d32_west->Write();
      fvtxs_d32_east->Write();
      fvtxs_d32_both->Write();

      // fvtxn_d32_west->Write();
      // fvtxn_d32_east->Write();
      // fvtxn_d32_both->Write();

      bbcs_sin22_west->Write();
      bbcs_sin22_east->Write();
      bbcs_sin22_both->Write();

      fvtxs_sin22_west->Write();
      fvtxs_sin22_east->Write();
      fvtxs_sin22_both->Write();

      // fvtxn_sin22_west->Write();
      // fvtxn_sin22_east->Write();
      // fvtxn_sin22_both->Write();

      bbcs_sin32_west->Write();
      bbcs_sin32_east->Write();
      bbcs_sin32_both->Write();

      fvtxs_sin32_west->Write();
      fvtxs_sin32_east->Write();
      fvtxs_sin32_both->Write();

      // fvtxn_sin32_west->Write();
      // fvtxn_sin32_east->Write();
      // fvtxn_sin32_both->Write();

      bbcs_cos22_west->Write();
      bbcs_cos22_east->Write();
      bbcs_cos22_both->Write();

      fvtxs_cos22_west->Write();
      fvtxs_cos22_east->Write();
      fvtxs_cos22_both->Write();

      // fvtxn_cos22_west->Write();
      // fvtxn_cos22_east->Write();
      // fvtxn_cos22_both->Write();

      bbcs_cos32_west->Write();
      bbcs_cos32_east->Write();
      bbcs_cos32_both->Write();

      fvtxs_cos32_west->Write();
      fvtxs_cos32_east->Write();
      fvtxs_cos32_both->Write();

      // fvtxn_cos32_west->Write();
      // fvtxn_cos32_east->Write();
      // fvtxn_cos32_both->Write();

      bbcs_c22->Write();
      fvtxs_c22->Write();
      // fvtxn_c22->Write();
      bbcs_c32->Write();
      fvtxs_c32->Write();
      // fvtxn_c32->Write();

      bbcs_sin22->Write();
      fvtxs_sin22->Write();
      // fvtxn_sin22->Write();
      bbcs_sin32->Write();
      fvtxs_sin32->Write();
      // fvtxn_sin32->Write();

      bbcs_cos22->Write();
      fvtxs_cos22->Write();
      // fvtxn_cos22->Write();
      bbcs_cos32->Write();
      fvtxs_cos32->Write();
      // fvtxn_cos32->Write();

      // ---

      os_bbcs_d22_west->Write();
      os_bbcs_d22_east->Write();
      os_bbcs_d22_both->Write();

      os_fvtxs_d22_west->Write();
      os_fvtxs_d22_east->Write();
      os_fvtxs_d22_both->Write();

      os_fvtxn_d22_west->Write();
      os_fvtxn_d22_east->Write();
      os_fvtxn_d22_both->Write();

      os_bbcs_d32_west->Write();
      os_bbcs_d32_east->Write();
      os_bbcs_d32_both->Write();

      os_fvtxs_d32_west->Write();
      os_fvtxs_d32_east->Write();
      os_fvtxs_d32_both->Write();

      os_fvtxn_d32_west->Write();
      os_fvtxn_d32_east->Write();
      os_fvtxn_d32_both->Write();

      os_bbcs_sin22_west->Write();
      os_bbcs_sin22_east->Write();
      os_bbcs_sin22_both->Write();

      os_fvtxs_sin22_west->Write();
      os_fvtxs_sin22_east->Write();
      os_fvtxs_sin22_both->Write();

      os_fvtxn_sin22_west->Write();
      os_fvtxn_sin22_east->Write();
      os_fvtxn_sin22_both->Write();

      os_bbcs_sin32_west->Write();
      os_bbcs_sin32_east->Write();
      os_bbcs_sin32_both->Write();

      os_fvtxs_sin32_west->Write();
      os_fvtxs_sin32_east->Write();
      os_fvtxs_sin32_both->Write();

      os_fvtxn_sin32_west->Write();
      os_fvtxn_sin32_east->Write();
      os_fvtxn_sin32_both->Write();

      os_bbcs_cos22_west->Write();
      os_bbcs_cos22_east->Write();
      os_bbcs_cos22_both->Write();

      os_fvtxs_cos22_west->Write();
      os_fvtxs_cos22_east->Write();
      os_fvtxs_cos22_both->Write();

      os_fvtxn_cos22_west->Write();
      os_fvtxn_cos22_east->Write();
      os_fvtxn_cos22_both->Write();

      os_bbcs_cos32_west->Write();
      os_bbcs_cos32_east->Write();
      os_bbcs_cos32_both->Write();

      os_fvtxs_cos32_west->Write();
      os_fvtxs_cos32_east->Write();
      os_fvtxs_cos32_both->Write();

      os_fvtxn_cos32_west->Write();
      os_fvtxn_cos32_east->Write();
      os_fvtxn_cos32_both->Write();

      os_bbcs_c22->Write();
      os_fvtxs_c22->Write();
      os_fvtxn_c22->Write();
      os_bbcs_c32->Write();
      os_fvtxs_c32->Write();
      os_fvtxn_c32->Write();

      os_bbcs_sin22->Write();
      os_fvtxs_sin22->Write();
      os_fvtxn_sin22->Write();
      os_bbcs_sin32->Write();
      os_fvtxs_sin32->Write();
      os_fvtxn_sin32->Write();

      os_bbcs_cos22->Write();
      os_fvtxs_cos22->Write();
      os_fvtxn_cos22->Write();
      os_bbcs_cos32->Write();
      os_fvtxs_cos32->Write();
      os_fvtxn_cos32->Write();

      os_bbcs_1dPsi2->Write();
      os_bbcs_1dPsi3->Write();
      os_fvtxs_1dPsi2->Write();
      os_fvtxs_1dPsi3->Write();
      os_fvtxn_1dPsi2->Write();
      os_fvtxn_1dPsi3->Write();

      os_bbcs_v2_west->Write();
      os_bbcs_v2_east->Write();
      os_bbcs_v2_both->Write();

      os_fvtxs_v2_west->Write();
      os_fvtxs_v2_east->Write();
      os_fvtxs_v2_both->Write();

      os_fvtxn_v2_west->Write();
      os_fvtxn_v2_east->Write();
      os_fvtxn_v2_both->Write();

      os_bbcs_v3_west->Write();
      os_bbcs_v3_east->Write();
      os_bbcs_v3_both->Write();

      os_fvtxs_v3_west->Write();
      os_fvtxs_v3_east->Write();
      os_fvtxs_v3_both->Write();

      os_fvtxn_v3_west->Write();
      os_fvtxn_v3_east->Write();
      os_fvtxn_v3_both->Write();

      os_bbcsfvtxs_c22->Write();
      os_bbcsfvtxs_c32->Write();
      os_bbcsfvtxn_c22->Write();
      os_bbcsfvtxn_c32->Write();
      os_fvtxsfvtxn_c22->Write();
      os_fvtxsfvtxn_c32->Write();

      os_bbcs_c24->Write();
      os_fvtxs_c24->Write();
      os_fvtxn_c24->Write();
      os_bbcs_c24_vs_cent->Write();
      os_fvtxs_c24_vs_cent->Write();
      os_fvtxn_c24_vs_cent->Write();
      os_bbcs_c24_vs_bbcqs->Write();
      os_fvtxs_c24_vs_bbcqs->Write();
      os_fvtxn_c24_vs_bbcqs->Write();
      os_bbcs_c24_vs_nfvtxs->Write();
      os_fvtxs_c24_vs_nfvtxs->Write();
      os_fvtxn_c24_vs_nfvtxs->Write();
      os_bbcs_d24_in_both->Write();
      os_fvtxs_d24_in_both->Write();
      os_fvtxn_d24_in_both->Write();
      os_bbcs_d24_out_both->Write();
      os_fvtxs_d24_out_both->Write();
      os_fvtxn_d24_out_both->Write();

      // ---

      mData1->Close();

    } // check on last pass


  //return;

  cout<<"cleaning up"<<endl;

  ntp_event_chain->Delete();
  // htree->Delete();
  // f->Close();
  // delete f;

  delete bbcs_v2_incl_nodetree;
  delete bbcs_v2_east_nodetree;
  delete bbcs_v2_west_nodetree;

  delete fvtxs_v2_incl_nodetree;
  delete fvtxs_v2_east_nodetree;
  delete fvtxs_v2_west_nodetree;

  delete bbcs_v2_west_docalib;
  delete fvtxs_v2_west_docalib;
  delete fvtxs0_v2_west_docalib;
  delete fvtxs1_v2_west_docalib;
  delete fvtxs2_v2_west_docalib;
  delete fvtxs3_v2_west_docalib;

  delete bbcs_v2_east_docalib;
  delete fvtxs_v2_east_docalib;
  delete fvtxs0_v2_east_docalib;
  delete fvtxs1_v2_east_docalib;
  delete fvtxs2_v2_east_docalib;
  delete fvtxs3_v2_east_docalib;

  delete bbcs_v2_both_docalib;
  delete fvtxs_v2_both_docalib;
  delete fvtxs0_v2_both_docalib;
  delete fvtxs1_v2_both_docalib;
  delete fvtxs2_v2_both_docalib;
  delete fvtxs3_v2_both_docalib;

  delete bbcs_v3_west_docalib;
  delete fvtxs_v3_west_docalib;
  // delete fvtxs0_v3_west_docalib;
  // delete fvtxs1_v3_west_docalib;
  // delete fvtxs2_v3_west_docalib;
  // delete fvtxs3_v3_west_docalib;

  delete bbcs_v3_east_docalib;
  delete fvtxs_v3_east_docalib;
  // delete fvtxs0_v3_east_docalib;
  // delete fvtxs1_v3_east_docalib;
  // delete fvtxs2_v3_east_docalib;
  // delete fvtxs3_v3_east_docalib;

  delete bbcs_v3_both_docalib;
  delete fvtxs_v3_both_docalib;
  // delete fvtxs0_v3_both_docalib;
  // delete fvtxs1_v3_both_docalib;
  // delete fvtxs2_v3_both_docalib;
  // delete fvtxs3_v3_both_docalib;

  for (int ic=0; ic<NMUL; ic++) {
  for (int iz=0; iz<NZPS; iz++) {
  for (int ih=1; ih<NHAR; ih++) {
  for (int id=0; id<NDET; id++) {
  delete ave[ic][iz][ih][id];
  delete flt[ic][iz][ih][id];
  }
  }
  }
  }

  for (int ic=0; ic<NMUL; ic++) {
  for (int ih=1; ih<NHAR; ih++) {
  for (int id=0; id<NDET; id++) {
  delete dis[ic][ih][id];
  delete qx[ic][ih][id];
  delete qy[ic][ih][id];
  delete psi_bf[ic][ih][id];
  delete psi_af[ic][ih][id];
  }
  }
  }

  delete th1d_BBC_charge;
  delete th1d_FVTX_nclus;
  delete th1d_FVTXN_nclus;
  delete th1d_FVTXS_nclus;
  delete th2d_qBBC_nFVTX;
  delete tp1f_reso2_BBC_CNT;
  delete tp1f_reso2_BBC_FVTX;
  delete tp1f_reso2_CNT_FVTX;
  delete tp1f_reso3_BBC_CNT;
  delete tp1f_reso3_BBC_FVTX;
  delete tp1f_reso3_CNT_FVTX;
  delete th1d_reso2_BBC_CNT;
  delete th1d_reso2_BBC_FVTX;
  delete th1d_reso2_CNT_FVTX;
  delete th1d_reso3_BBC_CNT;
  delete th1d_reso3_BBC_FVTX;
  delete th1d_reso3_CNT_FVTX;
  delete th1d_dreso2_BBC_CNT;
  delete th1d_dreso2_BBC_FVTX;
  delete th1d_dreso2_CNT_FVTX;
  delete th1d_dreso3_BBC_CNT;
  delete th1d_dreso3_BBC_FVTX;
  delete th1d_dreso3_CNT_FVTX;
  delete tp1f_reso2_BBC_FVTXN;
  delete tp1f_reso3_BBC_FVTXN;
  delete tp1f_reso2_FVTXS_FVTXN;
  delete tp1f_reso3_FVTXS_FVTXN;
  delete th1d_reso2_BBC_FVTXN;
  delete th1d_reso3_BBC_FVTXN;
  delete th1d_reso2_FVTXS_FVTXN;
  delete th1d_reso3_FVTXS_FVTXN;
  delete th1d_dreso2_BBC_FVTXN;
  delete th1d_dreso3_BBC_FVTXN;
  delete th1d_dreso2_FVTXS_FVTXN;
  delete th1d_dreso3_FVTXS_FVTXN;

  delete th1d_bbc_charge_phi;
  delete th1d_bbc0_charge_phi;
  delete th1d_bbc1_charge_phi;
  delete th1d_bbc2_charge_phi;
  delete th1d_bbc3_charge_phi;
  delete th1d_bbc4_charge_phi;

  delete th1d_fvtxs_clus_phi;
  delete th1d_fvtxs0_clus_phi;
  delete th1d_fvtxs1_clus_phi;
  delete th1d_fvtxs2_clus_phi;
  delete th1d_fvtxs3_clus_phi;

  delete th1d_fvtxn_clus_phi;
  delete th1d_fvtxn0_clus_phi;
  delete th1d_fvtxn1_clus_phi;
  delete th1d_fvtxn2_clus_phi;
  delete th1d_fvtxn3_clus_phi;

  delete tp1f_reso3_CNT_FVTXN;
  delete th1d_reso2_CNT_FVTXN;
  delete th1d_reso3_CNT_FVTXN;
  delete th1d_dreso2_CNT_FVTXN;
  delete th1d_dreso3_CNT_FVTXN;
  delete tp1f_reso42_BBC_CNT;
  delete tp1f_reso42_BBC_FVTX;
  delete tp1f_reso42_CNT_FVTX;
  delete tp1f_reso42_BBC_FVTXN;
  delete tp1f_reso42_CNT_FVTXN;
  delete tp1f_reso42_FVTXS_FVTXN;
  delete bbcs_v4_4Psi2_west_docalib;
  delete bbcs_v4_4Psi2_east_docalib;
  delete bbcs_v4_4Psi2_both_docalib;
  delete fvtxs_v4_4Psi2_west_docalib;
  delete fvtxs_v4_4Psi2_east_docalib;
  delete fvtxs_v4_4Psi2_both_docalib;
  delete fvtxn_v2_west_docalib;
  delete fvtxn_v2_east_docalib;
  delete fvtxn_v2_both_docalib;
  delete fvtxn_v3_west_docalib;
  delete fvtxn_v3_east_docalib;
  delete fvtxn_v3_both_docalib;
  delete fvtxn_v4_4Psi2_west_docalib;
  delete fvtxn_v4_4Psi2_east_docalib;
  delete fvtxn_v4_4Psi2_both_docalib;

  // ---

  delete bbcs_d22_west;
  delete bbcs_d22_east;
  delete bbcs_d22_both;

  delete fvtxs_d22_west;
  delete fvtxs_d22_east;
  delete fvtxs_d22_both;

  delete bbcs_d32_west;
  delete bbcs_d32_east;
  delete bbcs_d32_both;

  delete fvtxs_d32_west;
  delete fvtxs_d32_east;
  delete fvtxs_d32_both;

  delete bbcs_sin22_west;
  delete bbcs_sin22_east;
  delete bbcs_sin22_both;

  delete fvtxs_sin22_west;
  delete fvtxs_sin22_east;
  delete fvtxs_sin22_both;

  delete bbcs_sin32_west;
  delete bbcs_sin32_east;
  delete bbcs_sin32_both;

  delete fvtxs_sin32_west;
  delete fvtxs_sin32_east;
  delete fvtxs_sin32_both;

  delete bbcs_cos22_west;
  delete bbcs_cos22_east;
  delete bbcs_cos22_both;

  delete fvtxs_cos22_west;
  delete fvtxs_cos22_east;
  delete fvtxs_cos22_both;

  delete bbcs_cos32_west;
  delete bbcs_cos32_east;
  delete bbcs_cos32_both;

  delete fvtxs_cos32_west;
  delete fvtxs_cos32_east;
  delete fvtxs_cos32_both;

  delete bbcs_c22;
  delete fvtxs_c22;
  delete bbcs_c32;
  delete fvtxs_c32;

  delete bbcs_sin22;
  delete fvtxs_sin22;
  delete bbcs_sin32;
  delete fvtxs_sin32;

  delete bbcs_cos22;
  delete fvtxs_cos22;
  delete bbcs_cos32;
  delete fvtxs_cos32;

  // ---

  delete os_bbcs_d22_west;
  delete os_bbcs_d22_east;
  delete os_bbcs_d22_both;

  delete os_fvtxs_d22_west;
  delete os_fvtxs_d22_east;
  delete os_fvtxs_d22_both;

  delete os_bbcs_d32_west;
  delete os_bbcs_d32_east;
  delete os_bbcs_d32_both;

  delete os_fvtxs_d32_west;
  delete os_fvtxs_d32_east;
  delete os_fvtxs_d32_both;

  delete os_bbcs_sin22_west;
  delete os_bbcs_sin22_east;
  delete os_bbcs_sin22_both;

  delete os_fvtxs_sin22_west;
  delete os_fvtxs_sin22_east;
  delete os_fvtxs_sin22_both;

  delete os_bbcs_sin32_west;
  delete os_bbcs_sin32_east;
  delete os_bbcs_sin32_both;

  delete os_fvtxs_sin32_west;
  delete os_fvtxs_sin32_east;
  delete os_fvtxs_sin32_both;

  delete os_bbcs_cos22_west;
  delete os_bbcs_cos22_east;
  delete os_bbcs_cos22_both;

  delete os_fvtxs_cos22_west;
  delete os_fvtxs_cos22_east;
  delete os_fvtxs_cos22_both;

  delete os_bbcs_cos32_west;
  delete os_bbcs_cos32_east;
  delete os_bbcs_cos32_both;

  delete os_fvtxs_cos32_west;
  delete os_fvtxs_cos32_east;
  delete os_fvtxs_cos32_both;

  delete os_bbcs_c22;
  delete os_fvtxs_c22;
  delete os_bbcs_c32;
  delete os_fvtxs_c32;

  delete os_bbcs_sin22;
  delete os_fvtxs_sin22;
  delete os_bbcs_sin32;
  delete os_fvtxs_sin32;

  delete os_bbcs_cos22;
  delete os_fvtxs_cos22;
  delete os_bbcs_cos32;
  delete os_fvtxs_cos32;

  delete os_bbcs_1dPsi2;
  delete os_bbcs_1dPsi3;
  delete os_fvtxs_1dPsi2;
  delete os_fvtxs_1dPsi3;

  delete os_bbcs_v2_west;
  delete os_bbcs_v2_east;
  delete os_bbcs_v2_both;

  delete os_fvtxs_v2_west;
  delete os_fvtxs_v2_east;
  delete os_fvtxs_v2_both;

  delete os_bbcs_v3_west;
  delete os_bbcs_v3_east;
  delete os_bbcs_v3_both;

  delete os_fvtxs_v3_west;
  delete os_fvtxs_v3_east;
  delete os_fvtxs_v3_both;

  // ---

  delete tp1f_reso2_CNT_FVTXN;
  delete tp1f_os_reso2_BBC_CNT;
  delete tp1f_os_reso2_BBC_FVTX;
  delete tp1f_os_reso2_CNT_FVTX;
  delete tp1f_os_reso3_BBC_CNT;
  delete tp1f_os_reso3_BBC_FVTX;
  delete tp1f_os_reso3_CNT_FVTX;
  delete th1d_os_reso2_BBC_CNT;
  delete th1d_os_reso2_BBC_FVTX;
  delete th1d_os_reso2_CNT_FVTX;
  delete th1d_os_reso3_BBC_CNT;
  delete th1d_os_reso3_BBC_FVTX;
  delete th1d_os_reso3_CNT_FVTX;
  delete th1d_os_dreso2_BBC_CNT;
  delete th1d_os_dreso2_BBC_FVTX;
  delete th1d_os_dreso2_CNT_FVTX;
  delete th1d_os_dreso3_BBC_CNT;
  delete th1d_os_dreso3_BBC_FVTX;
  delete th1d_os_dreso3_CNT_FVTX;
  delete tp1f_os_reso2_BBC_FVTXN;
  delete tp1f_os_reso3_BBC_FVTXN;
  delete tp1f_os_reso2_CNT_FVTXN;
  delete tp1f_os_reso3_CNT_FVTXN;
  delete tp1f_os_reso2_FVTXS_FVTXN;
  delete tp1f_os_reso3_FVTXS_FVTXN;
  delete th1d_os_reso2_BBC_FVTXN;
  delete th1d_os_reso3_BBC_FVTXN;
  delete th1d_os_reso2_CNT_FVTXN;
  delete th1d_os_reso3_CNT_FVTXN;
  delete th1d_os_reso2_FVTXS_FVTXN;
  delete th1d_os_reso3_FVTXS_FVTXN;
  delete th1d_os_dreso2_BBC_FVTXN;
  delete th1d_os_dreso3_BBC_FVTXN;
  delete th1d_os_dreso2_CNT_FVTXN;
  delete th1d_os_dreso3_CNT_FVTXN;
  delete th1d_os_dreso2_FVTXS_FVTXN;
  delete th1d_os_dreso3_FVTXS_FVTXN;

  delete os_fvtxn_d22_west;
  delete os_fvtxn_d22_east;
  delete os_fvtxn_d22_both;
  delete os_fvtxn_d32_west;
  delete os_fvtxn_d32_east;
  delete os_fvtxn_d32_both;
  delete os_fvtxn_sin22_west;
  delete os_fvtxn_sin22_east;
  delete os_fvtxn_sin22_both;
  delete os_fvtxn_sin32_west;
  delete os_fvtxn_sin32_east;
  delete os_fvtxn_sin32_both;
  delete os_fvtxn_cos22_west;
  delete os_fvtxn_cos22_east;
  delete os_fvtxn_cos22_both;
  delete os_fvtxn_cos32_west;
  delete os_fvtxn_cos32_east;
  delete os_fvtxn_cos32_both;
  delete os_fvtxn_c22;
  delete os_fvtxn_c32;
  delete os_fvtxn_sin22;
  delete os_fvtxn_sin32;
  delete os_fvtxn_cos22;
  delete os_fvtxn_cos32;
  delete os_fvtxn_1dPsi2;
  delete os_fvtxn_1dPsi3;
  delete os_fvtxn_v2_west;
  delete os_fvtxn_v2_east;
  delete os_fvtxn_v2_both;
  delete os_fvtxn_v3_west;
  delete os_fvtxn_v3_east;
  delete os_fvtxn_v3_both;

  delete os_bbcsfvtxs_c22;
  delete os_bbcsfvtxs_c32;
  delete os_bbcsfvtxn_c22;
  delete os_bbcsfvtxn_c32;
  delete os_fvtxsfvtxn_c22;
  delete os_fvtxsfvtxn_c32;

  delete os_bbcs_c24;
  delete os_fvtxs_c24;
  delete os_fvtxn_c24;
  delete os_bbcs_c24_vs_cent;
  delete os_fvtxs_c24_vs_cent;
  delete os_fvtxn_c24_vs_cent;
  delete os_bbcs_c24_vs_bbcqs;
  delete os_fvtxs_c24_vs_bbcqs;
  delete os_fvtxn_c24_vs_bbcqs;
  delete os_bbcs_c24_vs_nfvtxs;
  delete os_fvtxs_c24_vs_nfvtxs;
  delete os_fvtxn_c24_vs_nfvtxs;
  delete os_bbcs_d24_in_both;
  delete os_fvtxs_d24_in_both;
  delete os_fvtxn_d24_in_both;
  delete os_bbcs_d24_out_both;
  delete os_fvtxs_d24_out_both;
  delete os_fvtxn_d24_out_both;


  // ---


  cout<<"end of program ana"<<endl;

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
       i == 11||
       i == 14||
       i == 19||
       i == 22||
       i == 26||
       i == 40||
       i == 43||
       i == 46||
       i == 51||
       i == 54||
       i == 58 ) return 0; // inner layer

  if ( i == 7 ||
       i == 16||
       i == 25||
       i == 39||
       i == 48||
       i == 57 ) return 1; // inner middle layer

  if ( i == 4 ||
       i == 6 ||
       i == 10||
       i == 13||
       i == 18||
       i == 21||
       i == 24||
       i == 29||
       i == 36||
       i == 38||
       i == 42||
       i == 45||
       i == 45||
       i == 50||
       i == 53||
       i == 56||
       i == 61 ) return 2; // middle layer

  if ( i == 3 ||
       i == 15||
       i == 28||
       i == 35||
       i == 47||
       i == 60 ) return 3; //outer middle layer

  if ( i == 0 ||
       i == 1 ||
       i == 2 ||
       i == 5 ||
       i == 9 ||
       i == 12||
       i == 17||
       i == 20||
       i == 23||
       i == 27||
       i == 30||
       i == 31||
       i == 32||
       i == 33||
       i == 34||
       i == 37||
       i == 41||
       i == 44||
       i == 49||
       i == 52||
       i == 55||
       i == 59||
       i == 62||
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
  cout<<"get_fvtx_layer::invalid z =  "<<z<<endl;
  return -1;
}


float calc4_event(float Xn, float Yn, float X2n, float Y2n, float M)
{

  float Qn2 = Xn*Xn+Yn*Yn;
  float Qn2d = Xn*Xn-Yn*Yn;

  float one   = Qn2*Qn2;
  float two   = X2n*X2n+Y2n*Y2n;
  float three = (2*(X2n*Qn2d + 2*Y2n*Xn*Yn));
  float four  = 2*(2*(M-2)*Qn2);
  float five  = 2*M*(M-3);

  float numerator = one + two - three - four + five;
  float denominator = M*(M-1)*(M-2)*(M-3);

  return numerator/denominator;

}


float calc4_track_flag(float xn, float yn, float x2n, float y2n, float Xn, float Yn, float X2n, float Y2n, float M, bool is_POI_in_RP)
{

  if ( is_POI_in_RP )
    {
      float one   = (xn*Xn + yn*Yn)*(Xn*Xn + Yn*Yn);
      float two   = x2n*Xn*Xn - x2n*Yn*Yn + 2*y2n*Xn*Yn;
      float three = xn*Xn*X2n + xn*Yn*Y2n - yn*(X2n*Yn - Xn*Y2n);
      float four  = 2*M*(xn*Xn + yn*Yn);
      float five  = 2*(Xn*Xn + Yn*Yn);
      float six   = 7*(xn*Xn + yn*Yn);
      float seven = xn*Xn + yn*Yn;
      float eight = x2n*X2n + y2n*Y2n;
      float nine = 2*(xn*Xn + yn*Yn);
      // ---
      float numerator = one - two - three - four - five + six - seven + eight + nine + 2*M - 6;
      float denominator = (M-1)*(M-2)*(M-3);
      // ---
      return numerator/denominator;
    }
  else
    {
      float one   = (xn*Xn + yn*Yn)*(Xn*Xn + Yn*Yn);
      float three = xn*Xn*X2n + xn*Yn*Y2n - yn*(X2n*Yn - Xn*Y2n);
      float four  = 2*M*(xn*Xn + yn*Yn);
      float nine = 2*(xn*Xn + yn*Yn);
      // ---
      float numerator = one - three - four + nine;
      float denominator = M*(M-1)*(M-2);
      // ---
      return numerator/denominator;
    }

}
