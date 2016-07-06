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

//===================================================================//
//                                                                   //
//       The purpose of this compiled root macro is to               //
//       read in event wise ttrees, read in BBCS and FVTXS           //
//       event planes, perform flattening and recentering            //
//       and make v2 w/r/t vtx tracks.                               //
//                                                                   //
//       You need to run 3 iterations (rp_recal_pass = 1, 2, 3).     //
//       The first two passes are just for calibration, the third    //
//       pass is to actually create the v2. Calib text files are     //
//       created                                                     //
//                                                                   //
//       This module also has the capability to calculate the FVTXS  //
//       BBCS event planes from clusters and tubes respectively (and //
//       apply corrections).                                         //
//                                                                   //
//       By Theo Koblesky, May 11, 2016                              //
//       theodore.koblesky@colorado.edu                              //
//                                                                   //
//===================================================================//




bool DIAG = false;


int get_fvtx_layer(float);
void initialize_pmt_position();
int get_pmt_layer(int);
void flatten(int, int);

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

  flatten(run,0); // just do step 0, only looking at cluster distributions for now...
  //flatten(run,1);
  //flatten(run,2);
  //flatten(run,3);

  return 0;

}

// -----------------------------------------------------------------
void flatten(int runNumber, int passNumber)
{

  cout << "Hello!  Now beginning..." << endl;
  cout << "runNumber = " << runNumber << " passNumber = " << passNumber << endl;

  int verbosity = 0;


  // char outFile1[300];
  // sprintf(outFile1,"%s%d%s","output/svrb_",runNumber,".root");
  TFile* tf_histo_inn = NULL;
  if ( passNumber > 0 ) tf_histo_inn = TFile::Open(Form("svrb_run%d_pass%d.root",runNumber,passNumber-1)); // need null for pass zero
  TFile* tf_histo_out = new TFile(Form("svrb_run%d_pass%d.root",runNumber,passNumber),"recreate");
  if ( !tf_histo_out )
    {
      cout << "FATAL ERROR: Unable to create output file!" << endl;
      return;
    }


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
      char filename[500];
      sprintf(filename,"input/tree_%010d_%04d.root",runNumber,ifile);
      cout<<"adding to tchain: "<<filename<<endl;
      ntp_event_chain->Add(filename);
    }


  cout << "Initalizing PMT positions for the BBCS" << endl;

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


  const float pi = TMath::Pi();

  //------------------------------------------------------------//
  //               Initializing Tree Variables                  //
  //------------------------------------------------------------//

  cout << "Now getting ready to read in the tree branch addresses and stuff...." << endl;



  //tree variables
  float        event;
  float        d_bbcz;    // bbcz
  float        centrality; // float because of dumb mistake in tree code :(
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

  int          d_nFVTX_clus;
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
  TBranch* b_d_FVTX_x;   //!
  TBranch* b_d_FVTX_y;   //!
  TBranch* b_d_FVTX_z;   //!
  TBranch* b_ntrk;   //!
  TBranch* b_px;   //!
  TBranch* b_py;   //!
  TBranch* b_pz;   //!
  TBranch* b_d_ntrk;   //!
  TBranch* b_d_cntpx;   //!
  TBranch* b_d_cntpy;   //!
  TBranch* b_d_cntpz;   //!

  ntp_event_chain->SetBranchAddress("bbc_z",&d_bbcz,&b_bbc_z);
  ntp_event_chain->SetBranchAddress("centrality",&centrality,&b_centrality);
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

  // --- BOAR = before offset and rotation
  TH1D* th1d_fvtxs_clus_phi_BOAR = new TH1D("th1d_fvtxs_clus_phi_BOAR","",50,-pi,pi);
  TH1D* th1d_fvtxs0_clus_phi_BOAR = new TH1D("th1d_fvtxs0_clus_phi_BOAR","",50,-pi,pi);
  TH1D* th1d_fvtxs1_clus_phi_BOAR = new TH1D("th1d_fvtxs1_clus_phi_BOAR","",50,-pi,pi);
  TH1D* th1d_fvtxs2_clus_phi_BOAR = new TH1D("th1d_fvtxs2_clus_phi_BOAR","",50,-pi,pi);
  TH1D* th1d_fvtxs3_clus_phi_BOAR = new TH1D("th1d_fvtxs3_clus_phi_BOAR","",50,-pi,pi);

  // --- IR = inside radius
  TH1D* th1d_fvtxs_clus_phi_IR = new TH1D("th1d_fvtxs_clus_phi_IR","",50,-pi,pi);
  TH1D* th1d_fvtxs0_clus_phi_IR = new TH1D("th1d_fvtxs0_clus_phi_IR","",50,-pi,pi);
  TH1D* th1d_fvtxs1_clus_phi_IR = new TH1D("th1d_fvtxs1_clus_phi_IR","",50,-pi,pi);
  TH1D* th1d_fvtxs2_clus_phi_IR = new TH1D("th1d_fvtxs2_clus_phi_IR","",50,-pi,pi);
  TH1D* th1d_fvtxs3_clus_phi_IR = new TH1D("th1d_fvtxs3_clus_phi_IR","",50,-pi,pi);

  // --- OR = outside radius
  TH1D* th1d_fvtxs_clus_phi_OR = new TH1D("th1d_fvtxs_clus_phi_OR","",50,-pi,pi);
  TH1D* th1d_fvtxs0_clus_phi_OR = new TH1D("th1d_fvtxs0_clus_phi_OR","",50,-pi,pi);
  TH1D* th1d_fvtxs1_clus_phi_OR = new TH1D("th1d_fvtxs1_clus_phi_OR","",50,-pi,pi);
  TH1D* th1d_fvtxs2_clus_phi_OR = new TH1D("th1d_fvtxs2_clus_phi_OR","",50,-pi,pi);
  TH1D* th1d_fvtxs3_clus_phi_OR = new TH1D("th1d_fvtxs3_clus_phi_OR","",50,-pi,pi);

  TH1D* th1d_FVTX_nclus = new TH1D("th1d_FVTX_nclus","",200,-0.5,1999.5);
  TH1D* th1d_FVTX_nclus_IR = new TH1D("th1d_FVTX_nclus_IR","",200,-0.5,1999.5);
  TH1D* th1d_FVTX_nclus_OR = new TH1D("th1d_FVTX_nclus_OR","",200,-0.5,1999.5);

  TH1D* th1d_BBC_charge = new TH1D("th1d_BBC_charge","",200,-0.5,199.5);
  TH1D* th1d_FVTXS_nclus = new TH1D("th1d_FVTXS_nclus","",200,-0.5,1999.5);
  TH1D* th1d_FVTXN_nclus = new TH1D("th1d_FVTXN_nclus","",200,-0.5,1999.5);

  TH1D* th1d_BBC_charge_IR = new TH1D("th1d_BBC_charge_IR","",200,-0.5,199.5);
  TH1D* th1d_FVTXS_nclus_IR = new TH1D("th1d_FVTXS_nclus_IR","",200,-0.5,1999.5);
  TH1D* th1d_FVTXN_nclus_IR = new TH1D("th1d_FVTXN_nclus_IR","",200,-0.5,1999.5);

  TH1D* th1d_BBC_charge_OR = new TH1D("th1d_BBC_charge_OR","",200,-0.5,199.5);
  TH1D* th1d_FVTXS_nclus_OR = new TH1D("th1d_FVTXS_nclus_OR","",200,-0.5,1999.5);
  TH1D* th1d_FVTXN_nclus_OR = new TH1D("th1d_FVTXN_nclus_OR","",200,-0.5,1999.5);


  // ---

  // --- q-vectors vs z for selected centrality
  // --- for the small systems, we are focusing on the central (0-5% or 0-10%)
  // --- for AuAu collisions, one would probably tune up the Q-vectors as a function
  // --- of centrality first, and then do the z-dependence
  // --- in any system, we do anticipate potentially significant distortions at large z-vertex
  // --- step0 - raw q-vectors with no corrections applied
  // --- Qx as a function of z
  TProfile* tp1f_step0_Qx2_zvtx_BBCS = new TProfile("tp1f_step0_Qx2_zvtx_BBCS","",120,-30,30,-1e10,1e10,"");
  TProfile* tp1f_step0_Qx3_zvtx_BBCS = new TProfile("tp1f_step0_Qx3_zvtx_BBCS","",120,-30,30,-1e10,1e10,"");
  TProfile* tp1f_step0_Qx2_zvtx_FVTXS = new TProfile("tp1f_step0_Qx2_zvtx_FVTXS","",120,-30,30,-1e10,1e10,"");
  TProfile* tp1f_step0_Qx3_zvtx_FVTXS = new TProfile("tp1f_step0_Qx3_zvtx_FVTXS","",120,-30,30,-1e10,1e10,"");
  TProfile* tp1f_step0_Qx2_zvtx_FVTXN = new TProfile("tp1f_step0_Qx2_zvtx_FVTXN","",120,-30,30,-1e10,1e10,"");
  TProfile* tp1f_step0_Qx3_zvtx_FVTXN = new TProfile("tp1f_step0_Qx3_zvtx_FVTXN","",120,-30,30,-1e10,1e10,"");
  // --- Qy as a function of z
  TProfile* tp1f_step0_Qy2_zvtx_BBCS = new TProfile("tp1f_step0_Qy2_zvtx_BBCS","",120,-30,30,-1e10,1e10,"");
  TProfile* tp1f_step0_Qy3_zvtx_BBCS = new TProfile("tp1f_step0_Qy3_zvtx_BBCS","",120,-30,30,-1e10,1e10,"");
  TProfile* tp1f_step0_Qy2_zvtx_FVTXS = new TProfile("tp1f_step0_Qy2_zvtx_FVTXS","",120,-30,30,-1e10,1e10,"");
  TProfile* tp1f_step0_Qy3_zvtx_FVTXS = new TProfile("tp1f_step0_Qy3_zvtx_FVTXS","",120,-30,30,-1e10,1e10,"");
  TProfile* tp1f_step0_Qy2_zvtx_FVTXN = new TProfile("tp1f_step0_Qy2_zvtx_FVTXN","",120,-30,30,-1e10,1e10,"");
  TProfile* tp1f_step0_Qy3_zvtx_FVTXN = new TProfile("tp1f_step0_Qy3_zvtx_FVTXN","",120,-30,30,-1e10,1e10,"");
  // --- Psi
  TH1D* th1d_step0_Psi2_BBCS = new TH1D("th1d_step0_Psi2_BBCS","",220,-4.1,4.1);
  TH1D* th1d_step0_Psi3_BBCS = new TH1D("th1d_step0_Psi3_BBCS","",220,-4.1,4.1);
  TH1D* th1d_step0_Psi2_FVTXS = new TH1D("th1d_step0_Psi2_FVTXS","",220,-4.1,4.1);
  TH1D* th1d_step0_Psi3_FVTXS = new TH1D("th1d_step0_Psi3_FVTXS","",220,-4.1,4.1);
  TH1D* th1d_step0_Psi2_FVTXN = new TH1D("th1d_step0_Psi2_FVTXN","",220,-4.1,4.1);
  TH1D* th1d_step0_Psi3_FVTXN = new TH1D("th1d_step0_Psi3_FVTXN","",220,-4.1,4.1);
  // --- step1 - q-vectors recentered
  // --- because the FVTX is a track/cluster detector, there's no practicable way to adjust the gain
  // --- so to treat the detectors equally, we do no gain correction to either, we just recenter the
  // --- q-vectors
  // --- Qx as a function of z
  TProfile* tp1f_step1_Qx2_zvtx_BBCS = new TProfile("tp1f_step1_Qx2_zvtx_BBCS","",120,-30,30,-1e10,1e10,"");
  TProfile* tp1f_step1_Qx3_zvtx_BBCS = new TProfile("tp1f_step1_Qx3_zvtx_BBCS","",120,-30,30,-1e10,1e10,"");
  TProfile* tp1f_step1_Qx2_zvtx_FVTXS = new TProfile("tp1f_step1_Qx2_zvtx_FVTXS","",120,-30,30,-1e10,1e10,"");
  TProfile* tp1f_step1_Qx3_zvtx_FVTXS = new TProfile("tp1f_step1_Qx3_zvtx_FVTXS","",120,-30,30,-1e10,1e10,"");
  TProfile* tp1f_step1_Qx2_zvtx_FVTXN = new TProfile("tp1f_step1_Qx2_zvtx_FVTXN","",120,-30,30,-1e10,1e10,"");
  TProfile* tp1f_step1_Qx3_zvtx_FVTXN = new TProfile("tp1f_step1_Qx3_zvtx_FVTXN","",120,-30,30,-1e10,1e10,"");
  // --- Qy as a function of z
  TProfile* tp1f_step1_Qy2_zvtx_BBCS = new TProfile("tp1f_step1_Qy2_zvtx_BBCS","",120,-30,30,-1e10,1e10,"");
  TProfile* tp1f_step1_Qy3_zvtx_BBCS = new TProfile("tp1f_step1_Qy3_zvtx_BBCS","",120,-30,30,-1e10,1e10,"");
  TProfile* tp1f_step1_Qy2_zvtx_FVTXS = new TProfile("tp1f_step1_Qy2_zvtx_FVTXS","",120,-30,30,-1e10,1e10,"");
  TProfile* tp1f_step1_Qy3_zvtx_FVTXS = new TProfile("tp1f_step1_Qy3_zvtx_FVTXS","",120,-30,30,-1e10,1e10,"");
  TProfile* tp1f_step1_Qy2_zvtx_FVTXN = new TProfile("tp1f_step1_Qy2_zvtx_FVTXN","",120,-30,30,-1e10,1e10,"");
  TProfile* tp1f_step1_Qy3_zvtx_FVTXN = new TProfile("tp1f_step1_Qy3_zvtx_FVTXN","",120,-30,30,-1e10,1e10,"");
  // --- Psi
  TH1D* th1d_step1_Psi2_BBCS = new TH1D("th1d_step1_Psi2_BBCS","",220,-4.1,4.1);
  TH1D* th1d_step1_Psi3_BBCS = new TH1D("th1d_step1_Psi3_BBCS","",220,-4.1,4.1);
  TH1D* th1d_step1_Psi2_FVTXS = new TH1D("th1d_step1_Psi2_FVTXS","",220,-4.1,4.1);
  TH1D* th1d_step1_Psi3_FVTXS = new TH1D("th1d_step1_Psi3_FVTXS","",220,-4.1,4.1);
  TH1D* th1d_step1_Psi2_FVTXN = new TH1D("th1d_step1_Psi2_FVTXN","",220,-4.1,4.1);
  TH1D* th1d_step1_Psi3_FVTXN = new TH1D("th1d_step1_Psi3_FVTXN","",220,-4.1,4.1);



  //------------------------------------------------------------//
  //                   Looping Over Event Tree                  //
  //------------------------------------------------------------//

  int all_counter = 0;
  int event_counter = 0;
  int bad_cent_counter = 0;
  int bad_vertex_counter = 0;
  int bad_radius_counter = 0;

  Long64_t cluster_counter = 0;
  Long64_t gapcut_counter = 0;

  cout << "starting loop over events in the tree" << endl;
  int nentries = ntp_event_chain->GetEntries();
  cout << "total events = " << nentries << endl;
  for ( int ievt = 0; ievt < nentries; ++ievt )
    {

      // --- break and continue statements should happen much, much earlier --------------------
      if ( passNumber < 0 || passNumber > 3 ) break;// passNumber only valid between 0 and 3... consider revising

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
      int ibbcz = NZPS*(ZVTX+10)/20;
      if ( ibbcz < 0 || ibbcz >= NZPS )
        {
          cout << "z vertex bin count problem!!!!" << endl;
          cout << "bbcz = " << d_bbcz << endl;
          cout << "fvtx_z = " << eventfvtx_z << endl;
          cout << "bin number is " << ibbcz << endl;
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

      bool inside_radius = true;
      if ( pow(eventfvtx_x-vtx_x,2.0) +  pow(eventfvtx_y-vtx_y,2.0) >= 0.15 ) inside_radius = false;


      // -------------
      // --- BBCS stuff
      // -------------

      float bbc_qx2 = 0;
      float bbc_qy2 = 0;
      float bbc_qx3 = 0;
      float bbc_qy3 = 0;
      float bbc_qw = 0;

      if ( ( say_event && verbosity > 0 ) || verbosity > 1 ) cout << "Looping over BBCS stuff now" << endl;

      if ( bbc_pmts )
        {
          for(int ipmt = 0; ipmt < 64; ipmt++)
            {
              float bbc_charge = d_BBC_charge[ipmt];
              if(bbc_charge <= 0) continue;

              // --- offset
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


      //cout << "HELLO HERE I AM" << endl;

      ++event_counter;

      // --------------
      // --- FVTX stuff
      // --------------

      float fvtxs_qx2[5];//all layers then 0 1 2 3
      float fvtxs_qy2[5];
      float fvtxs_qx3[5];//all layers then 0 1 2 3
      float fvtxs_qy3[5];
      float fvtxs_qw[5];

      for(int ilayer = 0; ilayer < 5; ilayer++)
        {
          fvtxs_qx2[ilayer] = 0.0;
          fvtxs_qy2[ilayer] = 0.0;
          fvtxs_qx3[ilayer] = 0.0;
          fvtxs_qy3[ilayer] = 0.0;
          fvtxs_qw[ilayer] = 0.0;
        } // loop over layers

      // --- now FVTX North

      float fvtxn_qx2[5];//all layers then 0 1 2 3
      float fvtxn_qy2[5];
      float fvtxn_qx3[5];//all layers then 0 1 2 3
      float fvtxn_qy3[5];
      float fvtxn_qw[5];

      for(int ilayer = 0; ilayer < 5; ilayer++)
        {
          fvtxn_qx2[ilayer] = 0.0;
          fvtxn_qy2[ilayer] = 0.0;
          fvtxn_qx3[ilayer] = 0.0;
          fvtxn_qy3[ilayer] = 0.0;
          fvtxn_qw[ilayer] = 0.0;
        } // loop over layers


      if ( ( say_event && verbosity > 0 ) || verbosity > 1 ) cout << "Looping over FVTX cluster" << endl;
      if ( fvtx_clusters )
        {
          for(int iclus = 0; iclus < d_nFVTX_clus; iclus++)
            {
              // --- offset
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

              float phib = TMath::ATan2(d_FVTX_x[iclus],d_FVTX_x[iclus]);
              float phi = TMath::ATan2(fvtx_y,fvtx_x);

              // --- south side
              if ( d_FVTX_z[iclus] < 0 )
                {
                  fvtxs_qx2[fvtx_layer+1] += TMath::Cos(2*phi);
                  fvtxs_qy2[fvtx_layer+1] += TMath::Sin(2*phi);
                  fvtxs_qx3[fvtx_layer+1] += TMath::Cos(3*phi);
                  fvtxs_qy3[fvtx_layer+1] += TMath::Sin(3*phi);

                  fvtxs_qx2[0] += TMath::Cos(2*phi);
                  fvtxs_qy2[0] += TMath::Sin(2*phi);
                  fvtxs_qx3[0] += TMath::Cos(3*phi);
                  fvtxs_qy3[0] += TMath::Sin(3*phi);

                  fvtxs_qw[fvtx_layer+1] ++;
                  fvtxs_qw[0] ++;

                  th1d_fvtxs_clus_phi_BOAR->Fill(phib);
                  if ( fvtx_layer == 0 ) th1d_fvtxs0_clus_phi_BOAR->Fill(phib);
                  if ( fvtx_layer == 1 ) th1d_fvtxs1_clus_phi_BOAR->Fill(phib);
                  if ( fvtx_layer == 2 ) th1d_fvtxs2_clus_phi_BOAR->Fill(phib);
                  if ( fvtx_layer == 3 ) th1d_fvtxs3_clus_phi_BOAR->Fill(phib);

                  if ( inside_radius )
                    {
                      th1d_fvtxs_clus_phi_IR->Fill(phi);
                      if ( fvtx_layer == 0 ) th1d_fvtxs0_clus_phi_IR->Fill(phi);
                      if ( fvtx_layer == 1 ) th1d_fvtxs1_clus_phi_IR->Fill(phi);
                      if ( fvtx_layer == 2 ) th1d_fvtxs2_clus_phi_IR->Fill(phi);
                      if ( fvtx_layer == 3 ) th1d_fvtxs3_clus_phi_IR->Fill(phi);
                    }
                  else
                    {
                      th1d_fvtxs_clus_phi_OR->Fill(phi);
                      if ( fvtx_layer == 0 ) th1d_fvtxs0_clus_phi_OR->Fill(phi);
                      if ( fvtx_layer == 1 ) th1d_fvtxs1_clus_phi_OR->Fill(phi);
                      if ( fvtx_layer == 2 ) th1d_fvtxs2_clus_phi_OR->Fill(phi);
                      if ( fvtx_layer == 3 ) th1d_fvtxs3_clus_phi_OR->Fill(phi);
                    }

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

      th1d_BBC_charge->Fill(bbc_qw);
      th1d_FVTX_nclus->Fill(d_nFVTX_clus);
      th1d_FVTXS_nclus->Fill(fvtxs_qw[0]);
      th1d_FVTXN_nclus->Fill(fvtxn_qw[0]);
      if ( !inside_radius )
        {
          ++bad_radius_counter;
          th1d_BBC_charge_OR->Fill(bbc_qw);
          th1d_FVTX_nclus_OR->Fill(d_nFVTX_clus);
          th1d_FVTXS_nclus_OR->Fill(fvtxs_qw[0]);
          th1d_FVTXN_nclus_OR->Fill(fvtxn_qw[0]);
        }
      else
        {
          th1d_BBC_charge_IR->Fill(bbc_qw);
          th1d_FVTX_nclus_IR->Fill(d_nFVTX_clus);
          th1d_FVTXS_nclus_IR->Fill(fvtxs_qw[0]);
          th1d_FVTXN_nclus_IR->Fill(fvtxn_qw[0]);
        }

      continue; // let's just look at clusters for a little while...

      // ---

      float step0_Qx2_BBCS = bbc_qx2/bbc_qw;
      float step0_Qx3_BBCS = bbc_qx3/bbc_qw;
      float step0_Qy2_BBCS = bbc_qy2/bbc_qw;
      float step0_Qy3_BBCS = bbc_qy3/bbc_qw;
      float step0_psi2_BBCS = atan2(step0_Qy2_BBCS,step0_Qx2_BBCS);
      float step0_psi3_BBCS = atan2(step0_Qy3_BBCS,step0_Qx3_BBCS);

      tp1f_step0_Qx2_zvtx_BBCS->Fill(ZVTX,step0_Qx2_BBCS);
      tp1f_step0_Qx3_zvtx_BBCS->Fill(ZVTX,step0_Qx3_BBCS);
      tp1f_step0_Qy2_zvtx_BBCS->Fill(ZVTX,step0_Qy2_BBCS);
      tp1f_step0_Qy3_zvtx_BBCS->Fill(ZVTX,step0_Qy3_BBCS);
      th1d_step0_Psi2_BBCS->Fill(step0_psi2_BBCS);
      th1d_step0_Psi3_BBCS->Fill(step0_psi3_BBCS);

      float step0_Qx2_FVTXS = fvtxs_qx2[0]/fvtxs_qw[0];
      float step0_Qx3_FVTXS = fvtxs_qx3[0]/fvtxs_qw[0];
      float step0_Qy2_FVTXS = fvtxs_qy2[0]/fvtxs_qw[0];
      float step0_Qy3_FVTXS = fvtxs_qy3[0]/fvtxs_qw[0];
      float step0_psi2_FVTXS = atan2(step0_Qy2_FVTXS,step0_Qx2_FVTXS);
      float step0_psi3_FVTXS = atan2(step0_Qy3_FVTXS,step0_Qx3_FVTXS);

      tp1f_step0_Qx2_zvtx_FVTXS->Fill(ZVTX,step0_Qx2_FVTXS);
      tp1f_step0_Qx3_zvtx_FVTXS->Fill(ZVTX,step0_Qx3_FVTXS);
      tp1f_step0_Qy2_zvtx_FVTXS->Fill(ZVTX,step0_Qy2_FVTXS);
      tp1f_step0_Qy3_zvtx_FVTXS->Fill(ZVTX,step0_Qy3_FVTXS);
      th1d_step0_Psi2_FVTXS->Fill(step0_psi2_FVTXS);
      th1d_step0_Psi3_FVTXS->Fill(step0_psi3_FVTXS);

      float step0_Qx2_FVTXN = fvtxn_qx2[0]/fvtxn_qw[0];
      float step0_Qx3_FVTXN = fvtxn_qx3[0]/fvtxn_qw[0];
      float step0_Qy2_FVTXN = fvtxn_qy2[0]/fvtxn_qw[0];
      float step0_Qy3_FVTXN = fvtxn_qy3[0]/fvtxn_qw[0];
      float step0_psi2_FVTXN = atan2(step0_Qy2_FVTXN,step0_Qx2_FVTXN);
      float step0_psi3_FVTXN = atan2(step0_Qy3_FVTXN,step0_Qx3_FVTXN);

      tp1f_step0_Qx2_zvtx_FVTXN->Fill(ZVTX,step0_Qx2_FVTXN);
      tp1f_step0_Qx3_zvtx_FVTXN->Fill(ZVTX,step0_Qx3_FVTXN);
      tp1f_step0_Qy2_zvtx_FVTXN->Fill(ZVTX,step0_Qy2_FVTXN);
      tp1f_step0_Qy3_zvtx_FVTXN->Fill(ZVTX,step0_Qy3_FVTXN);
      th1d_step0_Psi2_FVTXN->Fill(step0_psi2_FVTXN);
      th1d_step0_Psi3_FVTXN->Fill(step0_psi3_FVTXN);

      // ---
      // --- now moving on to step1 - recentering the q-vectors as a function of z
      // ---

      float step1_Qx2_BBCS = step0_Qx2_BBCS;
      float step1_Qx3_BBCS = step0_Qx3_BBCS;
      float step1_Qy2_BBCS = step0_Qy2_BBCS;
      float step1_Qy3_BBCS = step0_Qy3_BBCS;
      float step1_psi2_BBCS = step0_psi2_BBCS;
      float step1_psi3_BBCS = step0_psi3_BBCS;

      if ( tf_histo_inn )
        {
          // -----------------------
          // --- first the BBCS part
          // --- get the histograms from the file for the previous pass
          TProfile* FINN_tp1f_step0_Qx2_zvtx_BBCS = (TProfile*)tf_histo_inn->Get("tp1f_step0_Qx2_zvtx_BBCS");
          TProfile* FINN_tp1f_step0_Qx3_zvtx_BBCS = (TProfile*)tf_histo_inn->Get("tp1f_step0_Qx3_zvtx_BBCS");
          TProfile* FINN_tp1f_step0_Qy2_zvtx_BBCS = (TProfile*)tf_histo_inn->Get("tp1f_step0_Qy2_zvtx_BBCS");
          TProfile* FINN_tp1f_step0_Qy3_zvtx_BBCS = (TProfile*)tf_histo_inn->Get("tp1f_step0_Qy3_zvtx_BBCS");
          // --- get the average values for the current z vertex
          float FINN_float_step0_Qx2_zvtx_BBCS = FINN_tp1f_step0_Qx2_zvtx_BBCS->GetBinContent(FINN_tp1f_step0_Qx2_zvtx_BBCS->FindBin(ZVTX));
          float FINN_float_step0_Qx3_zvtx_BBCS = FINN_tp1f_step0_Qx3_zvtx_BBCS->GetBinContent(FINN_tp1f_step0_Qx3_zvtx_BBCS->FindBin(ZVTX));
          float FINN_float_step0_Qy2_zvtx_BBCS = FINN_tp1f_step0_Qy2_zvtx_BBCS->GetBinContent(FINN_tp1f_step0_Qy2_zvtx_BBCS->FindBin(ZVTX));
          float FINN_float_step0_Qy3_zvtx_BBCS = FINN_tp1f_step0_Qy3_zvtx_BBCS->GetBinContent(FINN_tp1f_step0_Qy3_zvtx_BBCS->FindBin(ZVTX));
          // --- calculate new q-vectors and event planes
          bbc_qx2 -= bbc_qw*FINN_float_step0_Qx2_zvtx_BBCS;
          bbc_qy2 -= bbc_qw*FINN_float_step0_Qy2_zvtx_BBCS;
          bbc_qx3 -= bbc_qw*FINN_float_step0_Qx3_zvtx_BBCS;
          bbc_qy3 -= bbc_qw*FINN_float_step0_Qy3_zvtx_BBCS;
          step1_Qx2_BBCS = bbc_qx2/bbc_qw;
          step1_Qx3_BBCS = bbc_qx3/bbc_qw;
          step1_Qy2_BBCS = bbc_qy2/bbc_qw;
          step1_Qy3_BBCS = bbc_qy3/bbc_qw;
          step1_psi2_BBCS = atan2(step0_Qy2_BBCS,step0_Qx2_BBCS);
          step1_psi3_BBCS = atan2(step0_Qy3_BBCS,step0_Qx3_BBCS);
        }

      tp1f_step1_Qx2_zvtx_BBCS->Fill(ZVTX,step1_Qx2_BBCS);
      tp1f_step1_Qx3_zvtx_BBCS->Fill(ZVTX,step1_Qx3_BBCS);
      tp1f_step1_Qy2_zvtx_BBCS->Fill(ZVTX,step1_Qy2_BBCS);
      tp1f_step1_Qy3_zvtx_BBCS->Fill(ZVTX,step1_Qy3_BBCS);
      th1d_step1_Psi2_BBCS->Fill(step1_psi2_BBCS);
      th1d_step1_Psi3_BBCS->Fill(step1_psi3_BBCS);

      // ---

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

              // --- rotation on single particles here
              px = pz*sin(-beam_angle) + px*cos(-beam_angle);
              float phi_angle = TMath::ATan2(py,px);
              float pt_angle = sqrt(px*px+py*py);

              if ( pt_angle < 0.2 || pt_angle > 5.0 ) continue; // pt cut added 2016-06-30

            } // loop over tracks

        } // check on tracks

    }//end of event

  cout << "Processed " << event_counter << "/" << all_counter << " events (" << (float)event_counter/(float)all_counter << ")" << endl;
  cout << "Events with bad centrality = " << bad_cent_counter << " (" << (float)bad_cent_counter/(float)event_counter << ")" << endl;
  //cout << "Events with vertex disagreement = " << bad_vertex_counter << " (" << (float)bad_vertex_counter/(float)event_counter << ")" << endl;
  cout << "Events with outside fvtx radius = " << bad_radius_counter << " (" << (float)bad_radius_counter/(float)event_counter << ")" << endl;
  cout << "Gap cuts applied " << gapcut_counter << "/" << cluster_counter << " (" << (float)gapcut_counter/(float)cluster_counter << ")" << endl;

  cout << "tf_histo_inn " << tf_histo_inn << endl;
  cout << "tf_histo_out " << tf_histo_out << endl;

  cout << "checking and attempting to close in file" << endl;
  if ( tf_histo_inn ) tf_histo_inn->Close();
  cout << "attempting to write out file" << endl;
  tf_histo_out->Write();
  cout << "attempting to close out file" << endl;
  tf_histo_out->Close();


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
