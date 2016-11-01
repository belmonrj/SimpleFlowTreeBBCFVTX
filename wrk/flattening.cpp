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

  flatten(run,1);
  flatten(run,2);
  flatten(run,3);

  return 0;

}

// -----------------------------------------------------------------
void flatten(int runNumber, int rp_recal_pass)
{

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

  int verbosity = 0;

  char calibfile[500];
  sprintf(calibfile,"output/flattening_data/flattening_%d_%d.dat",runNumber,rp_recal_pass-1);

  cout << "calib text output file: " << calibfile << endl;

  char filename[500];

  // // --- get the number of files for this run number
  // string pipe_out = (string) gSystem->GetFromPipe(Form("ls input/tree_%010d_*.root | grep -c r",runNumber));
  // int nfiles = 0;
  // nfiles = atoi(pipe_out.c_str());
  // cout<<"nfiles: "<<nfiles<<endl;
  // if(nfiles==0) return;

  // // --- make a new TChain for the tree
  // TChain *ntp_event_chain = new TChain("ntp_event");
  // for ( int ifile = 0; ifile < nfiles; ++ifile )
  //   {
  //     sprintf(filename,"input/tree_%010d_%04d.root",runNumber,ifile);
  //     cout<<"adding to tchain: "<<filename<<endl;
  //     ntp_event_chain->Add(filename);
  //   }

  TChain *ntp_event_chain = new TChain("ntp_event");
  sprintf(filename,"input/tree_%d.root",runNumber);
  cout << "adding to tchain: " << filename << endl;
  ntp_event_chain->Add(filename);

  // ---

  bool doweights = true; // change to false if you don't want to bother
  TString phiweightfile_name = Form("SpecialProjects/WeightFiles/newweight2d_run%d.root", runNumber);
  TFile* phi_weight_file = TFile::Open(phiweightfile_name); // COME BACK HERE AND HAVE A LOOK
  if ( !phi_weight_file )
    {
      if ( doweights ) cout << "WARNING could not open phi weight file: " << phiweightfile_name.Data() << endl;
      doweights = false;
    }

  TH2D* th2d_fvtxs_phi_weight[5];
  TH2D* th2d_fvtxn_phi_weight[5];
  if ( doweights )
    {
      TString histname0 = Form("th2d_fvtxs_clus_phi_weight" );
      TString histname1 = Form("th2d_fvtxs0_clus_phi_weight");
      TString histname2 = Form("th2d_fvtxs1_clus_phi_weight");
      TString histname3 = Form("th2d_fvtxs2_clus_phi_weight");
      TString histname4 = Form("th2d_fvtxs3_clus_phi_weight");
      th2d_fvtxs_phi_weight[0] = (TH2D*)phi_weight_file->Get(histname0);
      th2d_fvtxs_phi_weight[1] = (TH2D*)phi_weight_file->Get(histname1);
      th2d_fvtxs_phi_weight[2] = (TH2D*)phi_weight_file->Get(histname2);
      th2d_fvtxs_phi_weight[3] = (TH2D*)phi_weight_file->Get(histname3);
      th2d_fvtxs_phi_weight[4] = (TH2D*)phi_weight_file->Get(histname4);
      bool have_all_weight_histos =
        th2d_fvtxs_phi_weight[0] &&
        th2d_fvtxs_phi_weight[1] &&
        th2d_fvtxs_phi_weight[2] &&
        th2d_fvtxs_phi_weight[3] &&
        th2d_fvtxs_phi_weight[4];
      if ( verbosity > 0 || ( !have_all_weight_histos ) )
        {
          cout << "memory address of th2d_fvtxs_phi_weight[0] is " << th2d_fvtxs_phi_weight[0] << endl;
          cout << "memory address of th2d_fvtxs_phi_weight[1] is " << th2d_fvtxs_phi_weight[1] << endl;
          cout << "memory address of th2d_fvtxs_phi_weight[2] is " << th2d_fvtxs_phi_weight[2] << endl;
          cout << "memory address of th2d_fvtxs_phi_weight[3] is " << th2d_fvtxs_phi_weight[3] << endl;
          cout << "memory address of th2d_fvtxs_phi_weight[4] is " << th2d_fvtxs_phi_weight[4] << endl;
          cout << "name of th2d_fvtxs_phi_weight[0] is " << histname0.Data() << endl;
          cout << "name of th2d_fvtxs_phi_weight[1] is " << histname1.Data() << endl;
          cout << "name of th2d_fvtxs_phi_weight[2] is " << histname2.Data() << endl;
          cout << "name of th2d_fvtxs_phi_weight[3] is " << histname3.Data() << endl;
          cout << "name of th2d_fvtxs_phi_weight[4] is " << histname4.Data() << endl;
          cout << "The phi weight file read in was " << phiweightfile_name.Data() << endl;
        } // verbosity check
      if ( !have_all_weight_histos )
        {
          cout << "WARNING: not all weight histograms present" << endl;
          doweights = false;
        } // check on existence of histos
          // --- now get the weights for FVTX North
      histname0 = Form("th2d_fvtxn_clus_phi_weight" );
      histname1 = Form("th2d_fvtxn0_clus_phi_weight");
      histname2 = Form("th2d_fvtxn1_clus_phi_weight");
      histname3 = Form("th2d_fvtxn2_clus_phi_weight");
      histname4 = Form("th2d_fvtxn3_clus_phi_weight");
      th2d_fvtxn_phi_weight[0] = (TH2D*)phi_weight_file->Get(histname0);
      th2d_fvtxn_phi_weight[1] = (TH2D*)phi_weight_file->Get(histname1);
      th2d_fvtxn_phi_weight[2] = (TH2D*)phi_weight_file->Get(histname2);
      th2d_fvtxn_phi_weight[3] = (TH2D*)phi_weight_file->Get(histname3);
      th2d_fvtxn_phi_weight[4] = (TH2D*)phi_weight_file->Get(histname4);
      have_all_weight_histos =
        th2d_fvtxn_phi_weight[0] &&
        th2d_fvtxn_phi_weight[1] &&
        th2d_fvtxn_phi_weight[2] &&
        th2d_fvtxn_phi_weight[3] &&
        th2d_fvtxn_phi_weight[4];
      if ( verbosity > 0 || ( !have_all_weight_histos ) )
        {
          cout << "memory address of th2d_fvtxn_phi_weight[0] is " << th2d_fvtxn_phi_weight[0] << endl;
          cout << "memory address of th2d_fvtxn_phi_weight[1] is " << th2d_fvtxn_phi_weight[1] << endl;
          cout << "memory address of th2d_fvtxn_phi_weight[2] is " << th2d_fvtxn_phi_weight[2] << endl;
          cout << "memory address of th2d_fvtxn_phi_weight[3] is " << th2d_fvtxn_phi_weight[3] << endl;
          cout << "memory address of th2d_fvtxn_phi_weight[4] is " << th2d_fvtxn_phi_weight[4] << endl;
          cout << "name of th2d_fvtxn_phi_weight[0] is " << histname0.Data() << endl;
          cout << "name of th2d_fvtxn_phi_weight[1] is " << histname1.Data() << endl;
          cout << "name of th2d_fvtxn_phi_weight[2] is " << histname2.Data() << endl;
          cout << "name of th2d_fvtxn_phi_weight[3] is " << histname3.Data() << endl;
          cout << "name of th2d_fvtxn_phi_weight[4] is " << histname4.Data() << endl;
          cout << "The phi weight file read in was " << phiweightfile_name.Data() << endl;
        } // verbosity check
      if ( !have_all_weight_histos )
        {
          cout << "WARNING: not all weight histograms present" << endl;
          doweights = false;
        } // check on existence of histos
    } // check on doweights

  TH1D* th1d_fvtxs_phi_weight[5];
  TH1D* th1d_fvtxn_phi_weight[5];
  if ( doweights )
    {
      TString histname0 = Form("th1d_fvtxs_clus_phi_weight" );
      TString histname1 = Form("th1d_fvtxs0_clus_phi_weight");
      TString histname2 = Form("th1d_fvtxs1_clus_phi_weight");
      TString histname3 = Form("th1d_fvtxs2_clus_phi_weight");
      TString histname4 = Form("th1d_fvtxs3_clus_phi_weight");
      th1d_fvtxs_phi_weight[0] = (TH1D*)phi_weight_file->Get(histname0); // COME BACK HERE AND HAVE A LOOK
      th1d_fvtxs_phi_weight[1] = (TH1D*)phi_weight_file->Get(histname1);
      th1d_fvtxs_phi_weight[2] = (TH1D*)phi_weight_file->Get(histname2);
      th1d_fvtxs_phi_weight[3] = (TH1D*)phi_weight_file->Get(histname3);
      th1d_fvtxs_phi_weight[4] = (TH1D*)phi_weight_file->Get(histname4);
      bool have_all_weight_histos =
        th1d_fvtxs_phi_weight[0] &&
        th1d_fvtxs_phi_weight[1] &&
        th1d_fvtxs_phi_weight[2] &&
        th1d_fvtxs_phi_weight[3] &&
        th1d_fvtxs_phi_weight[4];
      if ( verbosity > 0 || ( !have_all_weight_histos ) )
        {
          cout << "memory address of th1d_fvtxs_phi_weight[0] is " << th1d_fvtxs_phi_weight[0] << endl;
          cout << "memory address of th1d_fvtxs_phi_weight[1] is " << th1d_fvtxs_phi_weight[1] << endl;
          cout << "memory address of th1d_fvtxs_phi_weight[2] is " << th1d_fvtxs_phi_weight[2] << endl;
          cout << "memory address of th1d_fvtxs_phi_weight[3] is " << th1d_fvtxs_phi_weight[3] << endl;
          cout << "memory address of th1d_fvtxs_phi_weight[4] is " << th1d_fvtxs_phi_weight[4] << endl;
          cout << "name of th1d_fvtxs_phi_weight[0] is " << histname0.Data() << endl;
          cout << "name of th1d_fvtxs_phi_weight[1] is " << histname1.Data() << endl;
          cout << "name of th1d_fvtxs_phi_weight[2] is " << histname2.Data() << endl;
          cout << "name of th1d_fvtxs_phi_weight[3] is " << histname3.Data() << endl;
          cout << "name of th1d_fvtxs_phi_weight[4] is " << histname4.Data() << endl;
          cout << "The phi weight file read in was " << phiweightfile_name.Data() << endl;
        } // verbosity check
      if ( !have_all_weight_histos )
        {
          cout << "WARNING: not all weight histograms present" << endl;
          doweights = false;
        } // check on existence of histos
          // --- now get the weights for FVTX North
      histname0 = Form("th1d_fvtxn_clus_phi_weight" );
      histname1 = Form("th1d_fvtxn0_clus_phi_weight");
      histname2 = Form("th1d_fvtxn1_clus_phi_weight");
      histname3 = Form("th1d_fvtxn2_clus_phi_weight");
      histname4 = Form("th1d_fvtxn3_clus_phi_weight");
      th1d_fvtxn_phi_weight[0] = (TH1D*)phi_weight_file->Get(histname0);
      th1d_fvtxn_phi_weight[1] = (TH1D*)phi_weight_file->Get(histname1);
      th1d_fvtxn_phi_weight[2] = (TH1D*)phi_weight_file->Get(histname2);
      th1d_fvtxn_phi_weight[3] = (TH1D*)phi_weight_file->Get(histname3);
      th1d_fvtxn_phi_weight[4] = (TH1D*)phi_weight_file->Get(histname4);
      have_all_weight_histos =
        th1d_fvtxn_phi_weight[0] &&
        th1d_fvtxn_phi_weight[1] &&
        th1d_fvtxn_phi_weight[2] &&
        th1d_fvtxn_phi_weight[3] &&
        th1d_fvtxn_phi_weight[4];
      if ( verbosity > 0 || ( !have_all_weight_histos ) )
        {
          cout << "memory address of th1d_fvtxn_phi_weight[0] is " << th1d_fvtxn_phi_weight[0] << endl;
          cout << "memory address of th1d_fvtxn_phi_weight[1] is " << th1d_fvtxn_phi_weight[1] << endl;
          cout << "memory address of th1d_fvtxn_phi_weight[2] is " << th1d_fvtxn_phi_weight[2] << endl;
          cout << "memory address of th1d_fvtxn_phi_weight[3] is " << th1d_fvtxn_phi_weight[3] << endl;
          cout << "memory address of th1d_fvtxn_phi_weight[4] is " << th1d_fvtxn_phi_weight[4] << endl;
          cout << "name of th1d_fvtxn_phi_weight[0] is " << histname0.Data() << endl;
          cout << "name of th1d_fvtxn_phi_weight[1] is " << histname1.Data() << endl;
          cout << "name of th1d_fvtxn_phi_weight[2] is " << histname2.Data() << endl;
          cout << "name of th1d_fvtxn_phi_weight[3] is " << histname3.Data() << endl;
          cout << "name of th1d_fvtxn_phi_weight[4] is " << histname4.Data() << endl;
          cout << "The phi weight file read in was " << phiweightfile_name.Data() << endl;
        } // verbosity check
      if ( !have_all_weight_histos )
        {
          cout << "WARNING: not all weight histograms present" << endl;
          doweights = false;
        } // check on existence of histos
    } // check on doweights

  // ---

  TString tubegaincorrectionfile_name = Form("SpecialProjects/WeightFiles/bbctube_run%d.root", runNumber);
  TFile* tube_gaincorrection_file = TFile::Open(tubegaincorrectionfile_name);
  float tube_gaincorrection[64] = {0};
  if ( tube_gaincorrection_file )
    {
      TH1D* th1d_tube_gaincorrection = (TH1D*)tube_gaincorrection_file->Get("th1d_tubegaincorrection");
      if ( th1d_tube_gaincorrection )
          for ( int i = 0; i < 64; ++i )
            tube_gaincorrection[i] = th1d_tube_gaincorrection->GetBinContent(i+1);
      else cout << "WARNING could not get BBC tube gaincorrection histogram" << endl;
    }
  else cout << "WARNING could not get BBC tube gaincorrection file " << tubegaincorrectionfile_name.Data() << endl;
  // ---

  cout << "Initalizing PMT positions for the BBC" << endl;

  initialize_pmt_position();

  bool fvtx_clusters = true;
  bool bbc_pmts      = true;
  bool cnt_tracks    = true;
  bool fvtx_tracks   = true;


  int bbcs_index      =  2;
  int fvtxs_index     =  3;
  int fvtxs0_index    =  4;
  int fvtxs1_index    =  5;
  int fvtxs2_index    =  6;
  int fvtxs3_index    =  7;
  int fvtxn_index     =  8;
  int fvtxn0_index    =  9;
  int fvtxn1_index    = 10;
  int fvtxn2_index    = 11;
  int fvtxn3_index    = 12;
  int fvtxs_nw_index  = 13;
  int fvtxs0_nw_index = 14;
  int fvtxs1_nw_index = 15;
  int fvtxs2_nw_index = 16;
  int fvtxs3_nw_index = 17;
  int fvtxn_nw_index  = 18;
  int fvtxn0_nw_index = 19;
  int fvtxn1_nw_index = 20;
  int fvtxn2_nw_index = 21;
  int fvtxn3_nw_index = 22;
  int bbcs_nw_index   = 23;
  int fvtxs_tracks_index  = 24; // come back here
  int fvtxs0_tracks_index = 25; // inner
  int fvtxs1_tracks_index = 26; // outer
  int fvtxn_tracks_index  = 27;
  int fvtxn0_tracks_index = 28;
  int fvtxn1_tracks_index = 29;


  float pi = acos(-1.0); // defined in RpPar.h, i wonder why no warning

  //------------------------------------------------------------//
  //                                                            //
  //       Initializing Calibration Arrays & Histograms         //
  //                                                            //
  //------------------------------------------------------------//

  // --- problems with these array dimensions??  lots of compile errors

  cout << "Lots of arrays and stuff" << endl;



  char outFile1[300];
  sprintf(outFile1,"%s%d%s%d%s","output/files_",energyflag,"/hist_",runNumber,".root");

  // char outFile2[100];
  // sprintf(outFile2,"%s%d%s","output/hrp_",runNumber,".root");

  // cout << "v2 histogram output file: " << outFile1 << endl;
  // cout << "EP histogram output file: " << outFile2 << endl;
  cout << "histogram output file: " << outFile1 << endl;

  TFile *mData1=TFile::Open(outFile1,"recreate");
  mData1->cd();

  // TFile *mData2=TFile::Open(outFile2,"recreate");
  // mData2->cd();

  // --- lots of comments needed here
  // --- flattening parameters output to file
  TH2D     *qx[NMUL][NHAR][NDET]; // Q vector x component, Q vector vs z_vertex bin
  TH2D     *qy[NMUL][NHAR][NDET]; // Q vector y component, Q vector vs z_vertex bin
  TProfile *ave[NMUL][NZPS][NHAR][NDET]; // average Psi
  TProfile *flt[NMUL][NZPS][NHAR][NDET]; // flattening parameters
  TH2D     *dis[NMUL][NHAR][NDET]; // displacement?  function of z_vertex bin

  TH2D     *psi_bf[NMUL][NHAR][NDET];
  TH2D     *psi_mf[NMUL][NHAR][NDET];
  TH2D     *psi_af[NMUL][NHAR][NDET];

  // flattening parameters read in from file
  float    mean[NMUL][NZPS][NHAR][NDET][2]; // mean of Psi distribution (???)
  float    widt[NMUL][NZPS][NHAR][NDET][2]; // width of Psi distribution (???)
  float    four[NMUL][NZPS][NHAR][NDET][2][NORD]; // ?

  // TFile *mData1=TFile::Open(outFile1,"recreate");
  // mData1->cd();


  // ---
  TH1D* th1d_BBC_charge = new TH1D("th1d_BBC_charge","",200,-0.5,199.5);
  TH1D* th1d_FVTX_nclus = new TH1D("th1d_FVTX_nclus","",200,-0.5,1999.5);
  TH2D* th2d_qBBC_nFVTX = new TH2D("th2d_qBBC_nFVTX","",200,-0.5,199.5,200,-0.5,1999.5);
  TH1D* th1d_FVTXS_nclus = new TH1D("th1d_FVTXS_nclus","",200,-0.5,1999.5);
  TH1D* th1d_FVTXN_nclus = new TH1D("th1d_FVTXN_nclus","",200,-0.5,1999.5);

  TProfile* tp1f_bbc_charge_phi = new TProfile("tp1f_bbc_charge_phi","",50,-pi,pi);
  TProfile* tp1f_bbc0_charge_phi = new TProfile("tp1f_bbc0_charge_phi","",50,-pi,pi);
  TProfile* tp1f_bbc1_charge_phi = new TProfile("tp1f_bbc1_charge_phi","",50,-pi,pi);
  TProfile* tp1f_bbc2_charge_phi = new TProfile("tp1f_bbc2_charge_phi","",50,-pi,pi);
  TProfile* tp1f_bbc3_charge_phi = new TProfile("tp1f_bbc3_charge_phi","",50,-pi,pi);
  TProfile* tp1f_bbc4_charge_phi = new TProfile("tp1f_bbc4_charge_phi","",50,-pi,pi);

  TProfile* tp1f_bbc_charge_tube = new TProfile("tp1f_bbc_charge_tube","",64,-0.5,63.5);
  TProfile* tp1f_bbc0_charge_tube = new TProfile("tp1f_bbc0_charge_tube","",64,-0.5,63.5);
  TProfile* tp1f_bbc1_charge_tube = new TProfile("tp1f_bbc1_charge_tube","",64,-0.5,63.5);
  TProfile* tp1f_bbc2_charge_tube = new TProfile("tp1f_bbc2_charge_tube","",64,-0.5,63.5);
  TProfile* tp1f_bbc3_charge_tube = new TProfile("tp1f_bbc3_charge_tube","",64,-0.5,63.5);
  TProfile* tp1f_bbc4_charge_tube = new TProfile("tp1f_bbc4_charge_tube","",64,-0.5,63.5);

  TProfile* tp1f_bbc_charge_wphi = new TProfile("tp1f_bbc_charge_wphi","",50,-pi,pi);
  TProfile* tp1f_bbc0_charge_wphi = new TProfile("tp1f_bbc0_charge_wphi","",50,-pi,pi);
  TProfile* tp1f_bbc1_charge_wphi = new TProfile("tp1f_bbc1_charge_wphi","",50,-pi,pi);
  TProfile* tp1f_bbc2_charge_wphi = new TProfile("tp1f_bbc2_charge_wphi","",50,-pi,pi);
  TProfile* tp1f_bbc3_charge_wphi = new TProfile("tp1f_bbc3_charge_wphi","",50,-pi,pi);
  TProfile* tp1f_bbc4_charge_wphi = new TProfile("tp1f_bbc4_charge_wphi","",50,-pi,pi);

  TProfile* tp1f_bbc_charge_wtube = new TProfile("tp1f_bbc_charge_wtube","",64,-0.5,63.5);
  TProfile* tp1f_bbc0_charge_wtube = new TProfile("tp1f_bbc0_charge_wtube","",64,-0.5,63.5);
  TProfile* tp1f_bbc1_charge_wtube = new TProfile("tp1f_bbc1_charge_wtube","",64,-0.5,63.5);
  TProfile* tp1f_bbc2_charge_wtube = new TProfile("tp1f_bbc2_charge_wtube","",64,-0.5,63.5);
  TProfile* tp1f_bbc3_charge_wtube = new TProfile("tp1f_bbc3_charge_wtube","",64,-0.5,63.5);
  TProfile* tp1f_bbc4_charge_wtube = new TProfile("tp1f_bbc4_charge_wtube","",64,-0.5,63.5);

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

  TH1D* th1d_fvtxs_clus_wphi = new TH1D("th1d_fvtxs_clus_wphi","",50,-pi,pi);
  TH1D* th1d_fvtxs0_clus_wphi = new TH1D("th1d_fvtxs0_clus_wphi","",50,-pi,pi);
  TH1D* th1d_fvtxs1_clus_wphi = new TH1D("th1d_fvtxs1_clus_wphi","",50,-pi,pi);
  TH1D* th1d_fvtxs2_clus_wphi = new TH1D("th1d_fvtxs2_clus_wphi","",50,-pi,pi);
  TH1D* th1d_fvtxs3_clus_wphi = new TH1D("th1d_fvtxs3_clus_wphi","",50,-pi,pi);

  TH1D* th1d_fvtxn_clus_wphi = new TH1D("th1d_fvtxn_clus_wphi","",50,-pi,pi);
  TH1D* th1d_fvtxn0_clus_wphi = new TH1D("th1d_fvtxn0_clus_wphi","",50,-pi,pi);
  TH1D* th1d_fvtxn1_clus_wphi = new TH1D("th1d_fvtxn1_clus_wphi","",50,-pi,pi);
  TH1D* th1d_fvtxn2_clus_wphi = new TH1D("th1d_fvtxn2_clus_wphi","",50,-pi,pi);
  TH1D* th1d_fvtxn3_clus_wphi = new TH1D("th1d_fvtxn3_clus_wphi","",50,-pi,pi);

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

  // --- have a quick look at fvtxs clusters and trucks, because why not
  TProfile* tp1f_reso2_FVTX_FVTX = new TProfile("tp1f_reso2_FVTX_FVTX","",1,-0.5,0.5,-1e6,1e6,"");
  TH1D* th1d_reso2_FVTX_FVTX = new TH1D("th1d_reso2_FVTX_FVTX","",220,-1.1,1.1);
  TH1D* th1d_dreso2_FVTX_FVTX = new TH1D("th1d_dreso2_FVTX_FVTX","",252,-6.3,6.3);
  TProfile* tp1f_reso3_FVTX_FVTX = new TProfile("tp1f_reso3_FVTX_FVTX","",1,-0.5,0.5,-1e6,1e6,"");
  TH1D* th1d_reso3_FVTX_FVTX = new TH1D("th1d_reso3_FVTX_FVTX","",220,-1.1,1.1);
  TH1D* th1d_dreso3_FVTX_FVTX = new TH1D("th1d_dreso3_FVTX_FVTX","",252,-6.3,6.3);


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

  TProfile* tp1f_nw_reso2_BBC_CNT = new TProfile("tp1f_nw_reso2_BBC_CNT","",1,-0.5,0.5,-1e6,1e6,"");
  TProfile* tp1f_nw_reso2_BBC_FVTX = new TProfile("tp1f_nw_reso2_BBC_FVTX","",1,-0.5,0.5,-1e6,1e6,"");
  TProfile* tp1f_nw_reso2_CNT_FVTX = new TProfile("tp1f_nw_reso2_CNT_FVTX","",1,-0.5,0.5,-1e6,1e6,"");
  TProfile* tp1f_nw_reso2_BBC_FVTXN = new TProfile("tp1f_nw_reso2_BBC_FVTXN","",1,-0.5,0.5,-1e6,1e6,"");
  TProfile* tp1f_nw_reso2_CNT_FVTXN = new TProfile("tp1f_nw_reso2_CNT_FVTXN","",1,-0.5,0.5,-1e6,1e6,"");
  TProfile* tp1f_nw_reso2_FVTXS_FVTXN = new TProfile("tp1f_nw_reso2_FVTXS_FVTXN","",1,-0.5,0.5,-1e6,1e6,"");

  TProfile* tp1f_nw_reso3_BBC_CNT = new TProfile("tp1f_nw_reso3_BBC_CNT","",1,-0.5,0.5,-1e6,1e6,"");
  TProfile* tp1f_nw_reso3_BBC_FVTX = new TProfile("tp1f_nw_reso3_BBC_FVTX","",1,-0.5,0.5,-1e6,1e6,"");
  TProfile* tp1f_nw_reso3_CNT_FVTX = new TProfile("tp1f_nw_reso3_CNT_FVTX","",1,-0.5,0.5,-1e6,1e6,"");
  TProfile* tp1f_nw_reso3_BBC_FVTXN = new TProfile("tp1f_nw_reso3_BBC_FVTXN","",1,-0.5,0.5,-1e6,1e6,"");
  TProfile* tp1f_nw_reso3_CNT_FVTXN = new TProfile("tp1f_nw_reso3_CNT_FVTXN","",1,-0.5,0.5,-1e6,1e6,"");
  TProfile* tp1f_nw_reso3_FVTXS_FVTXN = new TProfile("tp1f_nw_reso3_FVTXS_FVTXN","",1,-0.5,0.5,-1e6,1e6,"");

  TProfile* tp1f_nw_reso42_BBC_CNT = new TProfile("tp1f_nw_reso42_BBC_CNT","",1,-0.5,0.5,-1e6,1e6,"");
  TProfile* tp1f_nw_reso42_BBC_FVTX = new TProfile("tp1f_nw_reso42_BBC_FVTX","",1,-0.5,0.5,-1e6,1e6,"");
  TProfile* tp1f_nw_reso42_CNT_FVTX = new TProfile("tp1f_nw_reso42_CNT_FVTX","",1,-0.5,0.5,-1e6,1e6,"");
  TProfile* tp1f_nw_reso42_BBC_FVTXN = new TProfile("tp1f_nw_reso42_BBC_FVTXN","",1,-0.5,0.5,-1e6,1e6,"");
  TProfile* tp1f_nw_reso42_CNT_FVTXN = new TProfile("tp1f_nw_reso42_CNT_FVTXN","",1,-0.5,0.5,-1e6,1e6,"");
  TProfile* tp1f_nw_reso42_FVTXS_FVTXN = new TProfile("tp1f_nw_reso42_FVTXS_FVTXN","",1,-0.5,0.5,-1e6,1e6,"");

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

  // icent upto NMUL... 0-5, 5-10, 10-20, 20-40, 40-60, 60-88

  // --- profile histograms for average of Psi and flattening parameters
  char name[200];
  // int icent = 0;
  // int ic = icent;
  for ( int ic = 0; ic < NMUL; ic ++ )
    {
      for ( int iz = 0; iz < NZPS; iz++ )
	{
	  for ( int ih = 1; ih < NHAR; ih++ )
	    {
	      for ( int id = 0; id < NDET; id++ )
		{
		  // --- average (of?)
		  sprintf(name,"ave_%d_%d_%d_%d",ic,iz,ih,id);
		  ave[ic][iz][ih][id] = new TProfile(name,name,4,-0.5,3.5,-10.1,10.1,"S");//for SMD -1.1,1.1
		  // --- flattening parameter (?)
		  sprintf(name,"flt_%d_%d_%d_%d",ic,iz,ih,id);
		  flt[ic][iz][ih][id] = new TProfile(name,name,4*NORD,-0.5,NORD*4.0-0.5,-1.1,1.1);
		} // loop over ndetectors
	    } // loop over harmonics
	} // loop over z_vertex bins
    } // loop over centrality bins

  // --- TH2D histograms for Q vector components
  for ( int ic = 0; ic < NMUL; ic++ )
    {
      for ( int ih = 1; ih < NHAR; ih++)
	{
	  for ( int id = 0; id < NDET; id++)
	    {
	      // --- dis (?)
	      sprintf(name,"dis_%d_%d_%d",ic,ih,id);
	      dis[ic][ih][id] = new TH2D(name,name,NZPS*3,-0.5,NZPS*3.0-0.5,50,-pi,pi);
	      // --- qx (q-vector x component)
	      sprintf(name,"qx_%d_%d_%d",ic,ih,id);
	      qx[ic][ih][id] = new TH2D(name,name,NZPS*3,-0.5,NZPS*3.0-0.5, 220,-4.1,4.1);
	      // --- qy (q-vector y component)
	      sprintf(name,"qy_%d_%d_%d",ic,ih,id);
	      qy[ic][ih][id] = new TH2D(name,name,NZPS*3,-0.5,NZPS*3.0-0.5, 220,-4.1,4.1);
	      // --- psi_bf (event plane before recentering and flattening)
	      sprintf(name,"psi_bf_%d_%d_%d",ic,ih,id);
	      psi_bf[ic][ih][id] = new TH2D(name,name,NZPS*3,-0.5,NZPS*3.0-0.5, 220,-4.1,4.1);
	      // --- psi_mf (event plane after recentering but before flattening)
	      sprintf(name,"psi_mf_%d_%d_%d",ic,ih,id);
	      psi_mf[ic][ih][id] = new TH2D(name,name,NZPS*3,-0.5,NZPS*3.0-0.5, 220,-4.1,4.1);
	      // --- psi_af (event plane after recentering and flattening)
	      sprintf(name,"psi_af_%d_%d_%d",ic,ih,id);
	      psi_af[ic][ih][id] = new TH2D(name,name,NZPS*3,-0.5,NZPS*3.0-0.5, 220,-4.1,4.1);
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
	      for ( int id = 0; id < NDET; id++ )
		{
		  for ( int ib = 0; ib < 2; ib++ )
		    {
		      // --- mean (of q-vectors)
		      mean[ic][iz][ih][id][ib] = 0.0;
		      // --- width (of q-vectors)
		      widt[ic][iz][ih][id][ib] = 1.0;
		      for ( int io = 0; io < NORD; io++ )
			{
			  // --- fourier components for flattening
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
      float f0,f1,f2,f3;//f4,f5,f6,f7;
      ifstream ifs;
      ifs.open(calibfile);
      for ( int ic = 0; ic < NMUL; ic++ )
	{
	  for ( int iz = 0; iz < NZPS; iz++ )
	    {
	      for ( int ih = 1; ih < NHAR; ih++ )
		{
		  for ( int id = 0; id < NDET; id++ )
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
		      for ( int ib = 0; ib < 2; ib++ )
			{
			  for ( int io = 0; io < NORD; io++ )
			    {
			      ifs>>four[ic][iz][ih][id][ib][io];
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

  // --- nodetree

  TProfile* bbcs_v2_incl_nodetree = new TProfile("bbcs_v2_incl_nodetree","bbcs_v2_incl_nodetree",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* bbcs_v2_east_nodetree = new TProfile("bbcs_v2_east_nodetree","bbcs_v2_east_nodetree",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* bbcs_v2_west_nodetree = new TProfile("bbcs_v2_west_nodetree","bbcs_v2_west_nodetree",15, 0.0, 3.0, -1.1, 1.1);

  TProfile* fvtxs_v2_incl_nodetree = new TProfile("fvtxs_v2_incl_nodetree","fvtxs_v2_incl_nodetree",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* fvtxs_v2_east_nodetree = new TProfile("fvtxs_v2_east_nodetree","fvtxs_v2_east_nodetree",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* fvtxs_v2_west_nodetree = new TProfile("fvtxs_v2_west_nodetree","fvtxs_v2_west_nodetree",15, 0.0, 3.0, -1.1, 1.1);

  // --- docalib

  TProfile* bbcs_v2_both_docalib_cent[NMUL];
  TProfile* fvtxs_v2_both_docalib_cent[NMUL];
  for ( int ic = 0; ic < NMUL; ++ic )
    {
      bbcs_v2_both_docalib_cent[ic] = new TProfile(Form("bbcs_v2_both_docalib_cent%d",ic),"",15, 0.0, 3.0, -1.1, 1.1);
      fvtxs_v2_both_docalib_cent[ic] = new TProfile(Form("fvtxs_v2_both_docalib_cent%d",ic),"",15, 0.0, 3.0, -1.1, 1.1);
    }


  TProfile* bbcs_v2_west_docalib = new TProfile(Form("bbcs_v2_west_docalib"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* bbcs_v2_east_docalib = new TProfile(Form("bbcs_v2_east_docalib"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* bbcs_v2_both_docalib = new TProfile(Form("bbcs_v2_both_docalib"),"",15, 0.0, 3.0, -1.1, 1.1);

  TProfile* fvtxs_v2_west_docalib = new TProfile(Form("fvtxs_v2_west_docalib"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* fvtxs_v2_east_docalib = new TProfile(Form("fvtxs_v2_east_docalib"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* fvtxs_v2_both_docalib = new TProfile(Form("fvtxs_v2_both_docalib"),"",15, 0.0, 3.0, -1.1, 1.1);

  TProfile* fvtxs0_v2_west_docalib = new TProfile(Form("fvtxs0_v2_west_docalib"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* fvtxs0_v2_east_docalib = new TProfile(Form("fvtxs0_v2_east_docalib"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* fvtxs0_v2_both_docalib = new TProfile(Form("fvtxs0_v2_both_docalib"),"",15, 0.0, 3.0, -1.1, 1.1);

  TProfile* fvtxs1_v2_west_docalib = new TProfile(Form("fvtxs1_v2_west_docalib"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* fvtxs1_v2_east_docalib = new TProfile(Form("fvtxs1_v2_east_docalib"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* fvtxs1_v2_both_docalib = new TProfile(Form("fvtxs1_v2_both_docalib"),"",15, 0.0, 3.0, -1.1, 1.1);

  TProfile* fvtxs2_v2_west_docalib = new TProfile(Form("fvtxs2_v2_west_docalib"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* fvtxs2_v2_east_docalib = new TProfile(Form("fvtxs2_v2_east_docalib"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* fvtxs2_v2_both_docalib = new TProfile(Form("fvtxs2_v2_both_docalib"),"",15, 0.0, 3.0, -1.1, 1.1);

  TProfile* fvtxs3_v2_west_docalib = new TProfile(Form("fvtxs3_v2_west_docalib"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* fvtxs3_v2_east_docalib = new TProfile(Form("fvtxs3_v2_east_docalib"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* fvtxs3_v2_both_docalib = new TProfile(Form("fvtxs3_v2_both_docalib"),"",15, 0.0, 3.0, -1.1, 1.1);

  TProfile* bbcs_v3_west_docalib = new TProfile(Form("bbcs_v3_west_docalib"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* bbcs_v3_east_docalib = new TProfile(Form("bbcs_v3_east_docalib"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* bbcs_v3_both_docalib = new TProfile(Form("bbcs_v3_both_docalib"),"",15, 0.0, 3.0, -1.1, 1.1);

  TProfile* fvtxs_v3_west_docalib = new TProfile(Form("fvtxs_v3_west_docalib"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* fvtxs_v3_east_docalib = new TProfile(Form("fvtxs_v3_east_docalib"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* fvtxs_v3_both_docalib = new TProfile(Form("fvtxs_v3_both_docalib"),"",15, 0.0, 3.0, -1.1, 1.1);

  TProfile* bbcs_v4_4Psi2_west_docalib = new TProfile(Form("bbcs_v4_4Psi2_west_docalib"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* bbcs_v4_4Psi2_east_docalib = new TProfile(Form("bbcs_v4_4Psi2_east_docalib"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* bbcs_v4_4Psi2_both_docalib = new TProfile(Form("bbcs_v4_4Psi2_both_docalib"),"",15, 0.0, 3.0, -1.1, 1.1);

  TProfile* fvtxs_v4_4Psi2_west_docalib = new TProfile(Form("fvtxs_v4_4Psi2_west_docalib"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* fvtxs_v4_4Psi2_east_docalib = new TProfile(Form("fvtxs_v4_4Psi2_east_docalib"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* fvtxs_v4_4Psi2_both_docalib = new TProfile(Form("fvtxs_v4_4Psi2_both_docalib"),"",15, 0.0, 3.0, -1.1, 1.1);

  TProfile* fvtxn_v2_west_docalib = new TProfile(Form("fvtxn_v2_west_docalib"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* fvtxn_v2_east_docalib = new TProfile(Form("fvtxn_v2_east_docalib"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* fvtxn_v2_both_docalib = new TProfile(Form("fvtxn_v2_both_docalib"),"",15, 0.0, 3.0, -1.1, 1.1);

  TProfile* fvtxn_v3_west_docalib = new TProfile(Form("fvtxn_v3_west_docalib"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* fvtxn_v3_east_docalib = new TProfile(Form("fvtxn_v3_east_docalib"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* fvtxn_v3_both_docalib = new TProfile(Form("fvtxn_v3_both_docalib"),"",15, 0.0, 3.0, -1.1, 1.1);

  TProfile* fvtxn_v4_4Psi2_west_docalib = new TProfile(Form("fvtxn_v4_4Psi2_west_docalib"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* fvtxn_v4_4Psi2_east_docalib = new TProfile(Form("fvtxn_v4_4Psi2_east_docalib"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* fvtxn_v4_4Psi2_both_docalib = new TProfile(Form("fvtxn_v4_4Psi2_both_docalib"),"",15, 0.0, 3.0, -1.1, 1.1);

  // --- dcnw (docalib noweight)

  TProfile* bbcs_v2_west_dcnw = new TProfile(Form("bbcs_v2_west_dcnw"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* bbcs_v2_east_dcnw = new TProfile(Form("bbcs_v2_east_dcnw"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* bbcs_v2_both_dcnw = new TProfile(Form("bbcs_v2_both_dcnw"),"",15, 0.0, 3.0, -1.1, 1.1);

  TProfile* fvtxs_v2_west_dcnw = new TProfile(Form("fvtxs_v2_west_dcnw"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* fvtxs_v2_east_dcnw = new TProfile(Form("fvtxs_v2_east_dcnw"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* fvtxs_v2_both_dcnw = new TProfile(Form("fvtxs_v2_both_dcnw"),"",15, 0.0, 3.0, -1.1, 1.1);

  TProfile* fvtxs0_v2_west_dcnw = new TProfile(Form("fvtxs0_v2_west_dcnw"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* fvtxs0_v2_east_dcnw = new TProfile(Form("fvtxs0_v2_east_dcnw"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* fvtxs0_v2_both_dcnw = new TProfile(Form("fvtxs0_v2_both_dcnw"),"",15, 0.0, 3.0, -1.1, 1.1);

  TProfile* fvtxs1_v2_west_dcnw = new TProfile(Form("fvtxs1_v2_west_dcnw"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* fvtxs1_v2_east_dcnw = new TProfile(Form("fvtxs1_v2_east_dcnw"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* fvtxs1_v2_both_dcnw = new TProfile(Form("fvtxs1_v2_both_dcnw"),"",15, 0.0, 3.0, -1.1, 1.1);

  TProfile* fvtxs2_v2_west_dcnw = new TProfile(Form("fvtxs2_v2_west_dcnw"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* fvtxs2_v2_east_dcnw = new TProfile(Form("fvtxs2_v2_east_dcnw"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* fvtxs2_v2_both_dcnw = new TProfile(Form("fvtxs2_v2_both_dcnw"),"",15, 0.0, 3.0, -1.1, 1.1);

  TProfile* fvtxs3_v2_west_dcnw = new TProfile(Form("fvtxs3_v2_west_dcnw"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* fvtxs3_v2_east_dcnw = new TProfile(Form("fvtxs3_v2_east_dcnw"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* fvtxs3_v2_both_dcnw = new TProfile(Form("fvtxs3_v2_both_dcnw"),"",15, 0.0, 3.0, -1.1, 1.1);

  TProfile* bbcs_v3_west_dcnw = new TProfile(Form("bbcs_v3_west_dcnw"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* bbcs_v3_east_dcnw = new TProfile(Form("bbcs_v3_east_dcnw"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* bbcs_v3_both_dcnw = new TProfile(Form("bbcs_v3_both_dcnw"),"",15, 0.0, 3.0, -1.1, 1.1);

  TProfile* fvtxs_v3_west_dcnw = new TProfile(Form("fvtxs_v3_west_dcnw"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* fvtxs_v3_east_dcnw = new TProfile(Form("fvtxs_v3_east_dcnw"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* fvtxs_v3_both_dcnw = new TProfile(Form("fvtxs_v3_both_dcnw"),"",15, 0.0, 3.0, -1.1, 1.1);

  TProfile* bbcs_v4_4Psi2_west_dcnw = new TProfile(Form("bbcs_v4_4Psi2_west_dcnw"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* bbcs_v4_4Psi2_east_dcnw = new TProfile(Form("bbcs_v4_4Psi2_east_dcnw"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* bbcs_v4_4Psi2_both_dcnw = new TProfile(Form("bbcs_v4_4Psi2_both_dcnw"),"",15, 0.0, 3.0, -1.1, 1.1);

  TProfile* fvtxs_v4_4Psi2_west_dcnw = new TProfile(Form("fvtxs_v4_4Psi2_west_dcnw"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* fvtxs_v4_4Psi2_east_dcnw = new TProfile(Form("fvtxs_v4_4Psi2_east_dcnw"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* fvtxs_v4_4Psi2_both_dcnw = new TProfile(Form("fvtxs_v4_4Psi2_both_dcnw"),"",15, 0.0, 3.0, -1.1, 1.1);

  TProfile* fvtxn_v2_west_dcnw = new TProfile(Form("fvtxn_v2_west_dcnw"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* fvtxn_v2_east_dcnw = new TProfile(Form("fvtxn_v2_east_dcnw"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* fvtxn_v2_both_dcnw = new TProfile(Form("fvtxn_v2_both_dcnw"),"",15, 0.0, 3.0, -1.1, 1.1);

  TProfile* fvtxn_v3_west_dcnw = new TProfile(Form("fvtxn_v3_west_dcnw"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* fvtxn_v3_east_dcnw = new TProfile(Form("fvtxn_v3_east_dcnw"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* fvtxn_v3_both_dcnw = new TProfile(Form("fvtxn_v3_both_dcnw"),"",15, 0.0, 3.0, -1.1, 1.1);

  TProfile* fvtxn_v4_4Psi2_west_dcnw = new TProfile(Form("fvtxn_v4_4Psi2_west_dcnw"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* fvtxn_v4_4Psi2_east_dcnw = new TProfile(Form("fvtxn_v4_4Psi2_east_dcnw"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* fvtxn_v4_4Psi2_both_dcnw = new TProfile(Form("fvtxn_v4_4Psi2_both_dcnw"),"",15, 0.0, 3.0, -1.1, 1.1);

  // ---------------------------------------------------------------------------------------------------------
  // --- COME BACK HERE FOR 2PC HISTOGRAMS

  TProfile* bbcs_d22_west = new TProfile(Form("bbcs_d22_west"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* bbcs_d22_east = new TProfile(Form("bbcs_d22_east"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* bbcs_d22_both = new TProfile(Form("bbcs_d22_both"),"",15, 0.0, 3.0, -1.1, 1.1);

  TProfile* fvtxs_d22_west = new TProfile(Form("fvtxs_d22_west"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* fvtxs_d22_east = new TProfile(Form("fvtxs_d22_east"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* fvtxs_d22_both = new TProfile(Form("fvtxs_d22_both"),"",15, 0.0, 3.0, -1.1, 1.1);

  TProfile* bbcs_d32_west = new TProfile(Form("bbcs_d32_west"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* bbcs_d32_east = new TProfile(Form("bbcs_d32_east"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* bbcs_d32_both = new TProfile(Form("bbcs_d32_both"),"",15, 0.0, 3.0, -1.1, 1.1);

  TProfile* fvtxs_d32_west = new TProfile(Form("fvtxs_d32_west"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* fvtxs_d32_east = new TProfile(Form("fvtxs_d32_east"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* fvtxs_d32_both = new TProfile(Form("fvtxs_d32_both"),"",15, 0.0, 3.0, -1.1, 1.1);

  TProfile* bbcs_sin22_west = new TProfile(Form("bbcs_sin22_west"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* bbcs_sin22_east = new TProfile(Form("bbcs_sin22_east"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* bbcs_sin22_both = new TProfile(Form("bbcs_sin22_both"),"",15, 0.0, 3.0, -1.1, 1.1);

  TProfile* fvtxs_sin22_west = new TProfile(Form("fvtxs_sin22_west"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* fvtxs_sin22_east = new TProfile(Form("fvtxs_sin22_east"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* fvtxs_sin22_both = new TProfile(Form("fvtxs_sin22_both"),"",15, 0.0, 3.0, -1.1, 1.1);

  TProfile* bbcs_sin32_west = new TProfile(Form("bbcs_sin32_west"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* bbcs_sin32_east = new TProfile(Form("bbcs_sin32_east"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* bbcs_sin32_both = new TProfile(Form("bbcs_sin32_both"),"",15, 0.0, 3.0, -1.1, 1.1);

  TProfile* fvtxs_sin32_west = new TProfile(Form("fvtxs_sin32_west"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* fvtxs_sin32_east = new TProfile(Form("fvtxs_sin32_east"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* fvtxs_sin32_both = new TProfile(Form("fvtxs_sin32_both"),"",15, 0.0, 3.0, -1.1, 1.1);

  TProfile* bbcs_cos22_west = new TProfile(Form("bbcs_cos22_west"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* bbcs_cos22_east = new TProfile(Form("bbcs_cos22_east"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* bbcs_cos22_both = new TProfile(Form("bbcs_cos22_both"),"",15, 0.0, 3.0, -1.1, 1.1);

  TProfile* fvtxs_cos22_west = new TProfile(Form("fvtxs_cos22_west"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* fvtxs_cos22_east = new TProfile(Form("fvtxs_cos22_east"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* fvtxs_cos22_both = new TProfile(Form("fvtxs_cos22_both"),"",15, 0.0, 3.0, -1.1, 1.1);

  TProfile* bbcs_cos32_west = new TProfile(Form("bbcs_cos32_west"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* bbcs_cos32_east = new TProfile(Form("bbcs_cos32_east"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* bbcs_cos32_both = new TProfile(Form("bbcs_cos32_both"),"",15, 0.0, 3.0, -1.1, 1.1);

  TProfile* fvtxs_cos32_west = new TProfile(Form("fvtxs_cos32_west"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* fvtxs_cos32_east = new TProfile(Form("fvtxs_cos32_east"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* fvtxs_cos32_both = new TProfile(Form("fvtxs_cos32_both"),"",15, 0.0, 3.0, -1.1, 1.1);

  TProfile* bbcs_c22 = new TProfile(Form("bbcs_c22"),"",1, -0.5, 0.5, -1.1, 1.1);
  TProfile* fvtxs_c22 = new TProfile(Form("fvtxs_c22"),"",1, -0.5, 0.5, -1.1, 1.1);
  TProfile* bbcs_c32 = new TProfile(Form("bbcs_c32"),"",1, -0.5, 0.5, -1.1, 1.1);
  TProfile* fvtxs_c32 = new TProfile(Form("fvtxs_c32"),"",1, -0.5, 0.5, -1.1, 1.1);

  TProfile* bbcs_sin22 = new TProfile(Form("bbcs_sin22"),"",1, -0.5, 0.5, -1.1, 1.1);
  TProfile* fvtxs_sin22 = new TProfile(Form("fvtxs_sin22"),"",1, -0.5, 0.5, -1.1, 1.1);
  TProfile* bbcs_sin32 = new TProfile(Form("bbcs_sin32"),"",1, -0.5, 0.5, -1.1, 1.1);
  TProfile* fvtxs_sin32 = new TProfile(Form("fvtxs_sin32"),"",1, -0.5, 0.5, -1.1, 1.1);

  TProfile* bbcs_cos22 = new TProfile(Form("bbcs_cos22"),"",1, -0.5, 0.5, -1.1, 1.1);
  TProfile* fvtxs_cos22 = new TProfile(Form("fvtxs_cos22"),"",1, -0.5, 0.5, -1.1, 1.1);
  TProfile* bbcs_cos32 = new TProfile(Form("bbcs_cos32"),"",1, -0.5, 0.5, -1.1, 1.1);
  TProfile* fvtxs_cos32 = new TProfile(Form("fvtxs_cos32"),"",1, -0.5, 0.5, -1.1, 1.1);

  // ---------------------------------------------------------------------------------------------------------

  TProfile* os_bbcs_d22_west = new TProfile(Form("os_bbcs_d22_west"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* os_bbcs_d22_east = new TProfile(Form("os_bbcs_d22_east"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* os_bbcs_d22_both = new TProfile(Form("os_bbcs_d22_both"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* os_fvtxs_d22_west = new TProfile(Form("os_fvtxs_d22_west"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* os_fvtxs_d22_east = new TProfile(Form("os_fvtxs_d22_east"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* os_fvtxs_d22_both = new TProfile(Form("os_fvtxs_d22_both"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* os_fvtxn_d22_west = new TProfile(Form("os_fvtxn_d22_west"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* os_fvtxn_d22_east = new TProfile(Form("os_fvtxn_d22_east"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* os_fvtxn_d22_both = new TProfile(Form("os_fvtxn_d22_both"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* os_fvtxs_tracks_d22_west = new TProfile(Form("os_fvtxs_tracks_d22_west"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* os_fvtxs_tracks_d22_east = new TProfile(Form("os_fvtxs_tracks_d22_east"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* os_fvtxs_tracks_d22_both = new TProfile(Form("os_fvtxs_tracks_d22_both"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* os_fvtxn_tracks_d22_west = new TProfile(Form("os_fvtxn_tracks_d22_west"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* os_fvtxn_tracks_d22_east = new TProfile(Form("os_fvtxn_tracks_d22_east"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* os_fvtxn_tracks_d22_both = new TProfile(Form("os_fvtxn_tracks_d22_both"),"",15, 0.0, 3.0, -1.1, 1.1);

  TProfile* os_bbcs_d32_west = new TProfile(Form("os_bbcs_d32_west"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* os_bbcs_d32_east = new TProfile(Form("os_bbcs_d32_east"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* os_bbcs_d32_both = new TProfile(Form("os_bbcs_d32_both"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* os_fvtxs_d32_west = new TProfile(Form("os_fvtxs_d32_west"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* os_fvtxs_d32_east = new TProfile(Form("os_fvtxs_d32_east"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* os_fvtxs_d32_both = new TProfile(Form("os_fvtxs_d32_both"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* os_fvtxn_d32_west = new TProfile(Form("os_fvtxn_d32_west"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* os_fvtxn_d32_east = new TProfile(Form("os_fvtxn_d32_east"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* os_fvtxn_d32_both = new TProfile(Form("os_fvtxn_d32_both"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* os_fvtxs_tracks_d32_west = new TProfile(Form("os_fvtxs_tracks_d32_west"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* os_fvtxs_tracks_d32_east = new TProfile(Form("os_fvtxs_tracks_d32_east"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* os_fvtxs_tracks_d32_both = new TProfile(Form("os_fvtxs_tracks_d32_both"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* os_fvtxn_tracks_d32_west = new TProfile(Form("os_fvtxn_tracks_d32_west"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* os_fvtxn_tracks_d32_east = new TProfile(Form("os_fvtxn_tracks_d32_east"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* os_fvtxn_tracks_d32_both = new TProfile(Form("os_fvtxn_tracks_d32_both"),"",15, 0.0, 3.0, -1.1, 1.1);

  TProfile* os_bbcs_sin22_west = new TProfile(Form("os_bbcs_sin22_west"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* os_bbcs_sin22_east = new TProfile(Form("os_bbcs_sin22_east"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* os_bbcs_sin22_both = new TProfile(Form("os_bbcs_sin22_both"),"",15, 0.0, 3.0, -1.1, 1.1);

  TProfile* os_fvtxs_sin22_west = new TProfile(Form("os_fvtxs_sin22_west"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* os_fvtxs_sin22_east = new TProfile(Form("os_fvtxs_sin22_east"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* os_fvtxs_sin22_both = new TProfile(Form("os_fvtxs_sin22_both"),"",15, 0.0, 3.0, -1.1, 1.1);

  TProfile* os_fvtxn_sin22_west = new TProfile(Form("os_fvtxn_sin22_west"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* os_fvtxn_sin22_east = new TProfile(Form("os_fvtxn_sin22_east"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* os_fvtxn_sin22_both = new TProfile(Form("os_fvtxn_sin22_both"),"",15, 0.0, 3.0, -1.1, 1.1);

  TProfile* os_bbcs_sin32_west = new TProfile(Form("os_bbcs_sin32_west"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* os_bbcs_sin32_east = new TProfile(Form("os_bbcs_sin32_east"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* os_bbcs_sin32_both = new TProfile(Form("os_bbcs_sin32_both"),"",15, 0.0, 3.0, -1.1, 1.1);

  TProfile* os_fvtxs_sin32_west = new TProfile(Form("os_fvtxs_sin32_west"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* os_fvtxs_sin32_east = new TProfile(Form("os_fvtxs_sin32_east"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* os_fvtxs_sin32_both = new TProfile(Form("os_fvtxs_sin32_both"),"",15, 0.0, 3.0, -1.1, 1.1);

  TProfile* os_fvtxn_sin32_west = new TProfile(Form("os_fvtxn_sin32_west"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* os_fvtxn_sin32_east = new TProfile(Form("os_fvtxn_sin32_east"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* os_fvtxn_sin32_both = new TProfile(Form("os_fvtxn_sin32_both"),"",15, 0.0, 3.0, -1.1, 1.1);

  TProfile* os_bbcs_cos22_west = new TProfile(Form("os_bbcs_cos22_west"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* os_bbcs_cos22_east = new TProfile(Form("os_bbcs_cos22_east"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* os_bbcs_cos22_both = new TProfile(Form("os_bbcs_cos22_both"),"",15, 0.0, 3.0, -1.1, 1.1);

  TProfile* os_fvtxs_cos22_west = new TProfile(Form("os_fvtxs_cos22_west"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* os_fvtxs_cos22_east = new TProfile(Form("os_fvtxs_cos22_east"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* os_fvtxs_cos22_both = new TProfile(Form("os_fvtxs_cos22_both"),"",15, 0.0, 3.0, -1.1, 1.1);

  TProfile* os_fvtxn_cos22_west = new TProfile(Form("os_fvtxn_cos22_west"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* os_fvtxn_cos22_east = new TProfile(Form("os_fvtxn_cos22_east"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* os_fvtxn_cos22_both = new TProfile(Form("os_fvtxn_cos22_both"),"",15, 0.0, 3.0, -1.1, 1.1);

  TProfile* os_bbcs_cos32_west = new TProfile(Form("os_bbcs_cos32_west"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* os_bbcs_cos32_east = new TProfile(Form("os_bbcs_cos32_east"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* os_bbcs_cos32_both = new TProfile(Form("os_bbcs_cos32_both"),"",15, 0.0, 3.0, -1.1, 1.1);

  TProfile* os_fvtxs_cos32_west = new TProfile(Form("os_fvtxs_cos32_west"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* os_fvtxs_cos32_east = new TProfile(Form("os_fvtxs_cos32_east"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* os_fvtxs_cos32_both = new TProfile(Form("os_fvtxs_cos32_both"),"",15, 0.0, 3.0, -1.1, 1.1);

  TProfile* os_fvtxn_cos32_west = new TProfile(Form("os_fvtxn_cos32_west"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* os_fvtxn_cos32_east = new TProfile(Form("os_fvtxn_cos32_east"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* os_fvtxn_cos32_both = new TProfile(Form("os_fvtxn_cos32_both"),"",15, 0.0, 3.0, -1.1, 1.1);

  TProfile* os_bbcs_c22 = new TProfile(Form("os_bbcs_c22"),"",1, -0.5, 0.5, -1.1, 1.1);
  TProfile* os_fvtxs_c22 = new TProfile(Form("os_fvtxs_c22"),"",1, -0.5, 0.5, -1.1, 1.1);
  TProfile* os_fvtxn_c22 = new TProfile(Form("os_fvtxn_c22"),"",1, -0.5, 0.5, -1.1, 1.1);
  TProfile* os_fvtxc_c22 = new TProfile(Form("os_fvtxc_c22"),"",1, -0.5, 0.5, -1.1, 1.1);
  TProfile* os_fvtxa_c22 = new TProfile(Form("os_fvtxa_c22"),"",1, -0.5, 0.5, -1.1, 1.1);
  TProfile* os_fvtxs_ce01_c22 = new TProfile(Form("os_fvtxs_ce01_c22"),"" ,1, -0.5, 0.5, -1.1, 1.1);
  TProfile* os_fvtxn_ce01_c22 = new TProfile(Form("os_fvtxn_ce01_c22"),"" ,1, -0.5, 0.5, -1.1, 1.1);
  TProfile* os_fvtxs_tracks_c22 = new TProfile(Form("os_fvtxs_tracks_c22"),"" ,1, -0.5, 0.5, -1.1, 1.1);
  TProfile* os_fvtxn_tracks_c22 = new TProfile(Form("os_fvtxn_tracks_c22"),"" ,1, -0.5, 0.5, -1.1, 1.1);
  TProfile* os_fvtxc_tracks_c22 = new TProfile(Form("os_fvtxc_tracks_c22"),"" ,1, -0.5, 0.5, -1.1, 1.1);
  TProfile* os_fvtxa_tracks_c22 = new TProfile(Form("os_fvtxa_tracks_c22"),"" ,1, -0.5, 0.5, -1.1, 1.1);
  TProfile* os_bbcs_c32 = new TProfile(Form("os_bbcs_c32"),"",1, -0.5, 0.5, -1.1, 1.1);
  TProfile* os_fvtxs_c32 = new TProfile(Form("os_fvtxs_c32"),"",1, -0.5, 0.5, -1.1, 1.1);
  TProfile* os_fvtxn_c32 = new TProfile(Form("os_fvtxn_c32"),"",1, -0.5, 0.5, -1.1, 1.1);
  TProfile* os_fvtxc_c32 = new TProfile(Form("os_fvtxc_c32"),"",1, -0.5, 0.5, -1.1, 1.1);
  TProfile* os_fvtxa_c32 = new TProfile(Form("os_fvtxa_c32"),"",1, -0.5, 0.5, -1.1, 1.1);
  TProfile* os_fvtxs_ce01_c32 = new TProfile(Form("os_fvtxs_ce01_c32"),"" ,1, -0.5, 0.5, -1.1, 1.1);
  TProfile* os_fvtxn_ce01_c32 = new TProfile(Form("os_fvtxn_ce01_c32"),"" ,1, -0.5, 0.5, -1.1, 1.1);
  TProfile* os_fvtxs_tracks_c32 = new TProfile(Form("os_fvtxs_tracks_c32"),"" ,1, -0.5, 0.5, -1.1, 1.1);
  TProfile* os_fvtxn_tracks_c32 = new TProfile(Form("os_fvtxn_tracks_c32"),"" ,1, -0.5, 0.5, -1.1, 1.1);
  TProfile* os_fvtxc_tracks_c32 = new TProfile(Form("os_fvtxc_tracks_c32"),"" ,1, -0.5, 0.5, -1.1, 1.1);
  TProfile* os_fvtxa_tracks_c32 = new TProfile(Form("os_fvtxa_tracks_c32"),"" ,1, -0.5, 0.5, -1.1, 1.1);

  TProfile* os_bbcs_sin22 = new TProfile(Form("os_bbcs_sin22"),"",1, -0.5, 0.5, -1.1, 1.1);
  TProfile* os_fvtxs_sin22 = new TProfile(Form("os_fvtxs_sin22"),"",1, -0.5, 0.5, -1.1, 1.1);
  TProfile* os_fvtxn_sin22 = new TProfile(Form("os_fvtxn_sin22"),"",1, -0.5, 0.5, -1.1, 1.1);
  TProfile* os_fvtxs_tracks_sin22 = new TProfile(Form("os_fvtxs_tracks_sin22"),"" ,1, -0.5, 0.5, -1.1, 1.1);
  TProfile* os_fvtxn_tracks_sin22 = new TProfile(Form("os_fvtxn_tracks_sin22"),"" ,1, -0.5, 0.5, -1.1, 1.1);
  TProfile* os_bbcs_sin32 = new TProfile(Form("os_bbcs_sin32"),"",1, -0.5, 0.5, -1.1, 1.1);
  TProfile* os_fvtxs_sin32 = new TProfile(Form("os_fvtxs_sin32"),"",1, -0.5, 0.5, -1.1, 1.1);
  TProfile* os_fvtxn_sin32 = new TProfile(Form("os_fvtxn_sin32"),"",1, -0.5, 0.5, -1.1, 1.1);
  TProfile* os_fvtxs_tracks_sin32 = new TProfile(Form("os_fvtxs_tracks_sin32"),"" ,1, -0.5, 0.5, -1.1, 1.1);
  TProfile* os_fvtxn_tracks_sin32 = new TProfile(Form("os_fvtxn_tracks_sin32"),"" ,1, -0.5, 0.5, -1.1, 1.1);

  TProfile* os_bbcs_cos22 = new TProfile(Form("os_bbcs_cos22"),  "",1, -0.5, 0.5, -1.1, 1.1);
  TProfile* os_fvtxs_cos22 = new TProfile(Form("os_fvtxs_cos22"),"",1, -0.5, 0.5, -1.1, 1.1);
  TProfile* os_fvtxn_cos22 = new TProfile(Form("os_fvtxn_cos22"),"",1, -0.5, 0.5, -1.1, 1.1);
  TProfile* os_fvtxs_tracks_cos22 = new TProfile(Form("os_fvtxs_tracks_cos22"),     "" ,1, -0.5, 0.5, -1.1, 1.1);
  TProfile* os_fvtxn_tracks_cos22 = new TProfile(Form("os_fvtxn_tracks_cos22"),     "" ,1, -0.5, 0.5, -1.1, 1.1);
  TProfile* os_bbcs_cos32 = new TProfile(Form("os_bbcs_cos32"),  "",1, -0.5, 0.5, -1.1, 1.1);
  TProfile* os_fvtxs_cos32 = new TProfile(Form("os_fvtxs_cos32"),"",1, -0.5, 0.5, -1.1, 1.1);
  TProfile* os_fvtxn_cos32 = new TProfile(Form("os_fvtxn_cos32"),"",1, -0.5, 0.5, -1.1, 1.1);
  TProfile* os_fvtxs_tracks_cos32 = new TProfile(Form("os_fvtxs_tracks_cos32"),     "" ,1, -0.5, 0.5, -1.1, 1.1);
  TProfile* os_fvtxn_tracks_cos32 = new TProfile(Form("os_fvtxn_tracks_cos32"),     "" ,1, -0.5, 0.5, -1.1, 1.1);

  TProfile* os_bbcsfvtxs_c22 = new TProfile(Form("os_bbcsfvtxs_c22"),"",1, -0.5, 0.5, -1.1, 1.1);
  TProfile* os_bbcsfvtxs_c32 = new TProfile(Form("os_bbcsfvtxs_c32"),"",1, -0.5, 0.5, -1.1, 1.1);
  TProfile* os_bbcsfvtxn_c22 = new TProfile(Form("os_bbcsfvtxn_c22"),"",1, -0.5, 0.5, -1.1, 1.1);
  TProfile* os_bbcsfvtxn_c32 = new TProfile(Form("os_bbcsfvtxn_c32"),"",1, -0.5, 0.5, -1.1, 1.1);
  TProfile* os_fvtxsfvtxn_c22 = new TProfile(Form("os_fvtxsfvtxn_c22"),"",1, -0.5, 0.5, -1.1, 1.1);
  TProfile* os_fvtxsfvtxn_c32 = new TProfile(Form("os_fvtxsfvtxn_c32"),"",1, -0.5, 0.5, -1.1, 1.1);
  TProfile* os_fvtxsfvtxn_tracks_c22 = new TProfile(Form("os_fvtxsfvtxn_tracks_c22"),"",1, -0.5, 0.5, -1.1, 1.1);
  TProfile* os_fvtxsfvtxn_tracks_c32 = new TProfile(Form("os_fvtxsfvtxn_tracks_c32"),"",1, -0.5, 0.5, -1.1, 1.1);

  TH1D* os_bbcs_1dPsi2 = new TH1D(Form("os_bbcs_1dPsi2"),"",220,-4.1,4.1); // weird binning from elsewhere, leave the same to match
  TH1D* os_bbcs_1dPsi3 = new TH1D(Form("os_bbcs_1dPsi3"),"",220,-4.1,4.1); // weird binning from elsewhere, leave the same to match
  TH1D* os_fvtxs_1dPsi2 = new TH1D(Form("os_fvtxs_1dPsi2"),"",200,-4.1,4.1); // weird binning from elsewhere, leave the same to match
  TH1D* os_fvtxs_1dPsi3 = new TH1D(Form("os_fvtxs_1dPsi3"),"",200,-4.1,4.1); // weird binning from elsewhere, leave the same to match
  TH1D* os_fvtxn_1dPsi2 = new TH1D(Form("os_fvtxn_1dPsi2"),"",200,-4.1,4.1); // weird binning from elsewhere, leave the same to match
  TH1D* os_fvtxn_1dPsi3 = new TH1D(Form("os_fvtxn_1dPsi3"),"",200,-4.1,4.1); // weird binning from elsewhere, leave the same to match
  TH1D* os_fvtxs_tracks_1dPsi2 = new TH1D(Form("os_fvtxs_tracks_1dPsi2"),     "" ,200,-4.1,4.1); // weird binning from elsewhere, leave the same to match
  TH1D* os_fvtxs_tracks_1dPsi3 = new TH1D(Form("os_fvtxs_tracks_1dPsi3"),     "" ,200,-4.1,4.1); // weird binning from elsewhere, leave the same to match
  TH1D* os_fvtxn_tracks_1dPsi2 = new TH1D(Form("os_fvtxn_tracks_1dPsi2"),     "" ,200,-4.1,4.1); // weird binning from elsewhere, leave the same to match
  TH1D* os_fvtxn_tracks_1dPsi3 = new TH1D(Form("os_fvtxn_tracks_1dPsi3"),     "" ,200,-4.1,4.1); // weird binning from elsewhere, leave the same to match

  TProfile* os_bbcs_v2_west = new TProfile(Form("os_bbcs_v2_west"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* os_bbcs_v2_east = new TProfile(Form("os_bbcs_v2_east"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* os_bbcs_v2_both = new TProfile(Form("os_bbcs_v2_both"),"",15, 0.0, 3.0, -1.1, 1.1);

  TProfile* os_fvtxs_v2_west = new TProfile(Form("os_fvtxs_v2_west"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* os_fvtxs_v2_east = new TProfile(Form("os_fvtxs_v2_east"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* os_fvtxs_v2_both = new TProfile(Form("os_fvtxs_v2_both"),"",15, 0.0, 3.0, -1.1, 1.1);

  TProfile* os_fvtxn_v2_west = new TProfile(Form("os_fvtxn_v2_west"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* os_fvtxn_v2_east = new TProfile(Form("os_fvtxn_v2_east"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* os_fvtxn_v2_both = new TProfile(Form("os_fvtxn_v2_both"),"",15, 0.0, 3.0, -1.1, 1.1);

  TProfile* os_bbcs_v3_west = new TProfile(Form("os_bbcs_v3_west"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* os_bbcs_v3_east = new TProfile(Form("os_bbcs_v3_east"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* os_bbcs_v3_both = new TProfile(Form("os_bbcs_v3_both"),"",15, 0.0, 3.0, -1.1, 1.1);

  TProfile* os_fvtxs_v3_west = new TProfile(Form("os_fvtxs_v3_west"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* os_fvtxs_v3_east = new TProfile(Form("os_fvtxs_v3_east"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* os_fvtxs_v3_both = new TProfile(Form("os_fvtxs_v3_both"),"",15, 0.0, 3.0, -1.1, 1.1);

  TProfile* os_fvtxn_v3_west = new TProfile(Form("os_fvtxn_v3_west"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* os_fvtxn_v3_east = new TProfile(Form("os_fvtxn_v3_east"),"",15, 0.0, 3.0, -1.1, 1.1);
  TProfile* os_fvtxn_v3_both = new TProfile(Form("os_fvtxn_v3_both"),"",15, 0.0, 3.0, -1.1, 1.1);

  TProfile* os_bbcs_c24 = new TProfile(Form("os_bbcs_c24"),"",1,-0.5,0.5,-1.1,1.1);
  TProfile* os_fvtxs_c24 = new TProfile(Form("os_fvtxs_c24"),"",1,-0.5,0.5,-1.1,1.1);
  TProfile* os_fvtxn_c24 = new TProfile(Form("os_fvtxn_c24"),"",1,-0.5,0.5,-1.1,1.1);
  TProfile* os_bbcs_c24_vs_cent = new TProfile(Form("os_bbcs_c24_vs_cent"),"",101,-0.5,100.5,-1.1,1.1);
  TProfile* os_fvtxs_c24_vs_cent = new TProfile(Form("os_fvtxs_c24_vs_cent"),"",101,-0.5,100.5,-1.1,1.1);
  TProfile* os_fvtxn_c24_vs_cent = new TProfile(Form("os_fvtxn_c24_vs_cent"),"",101,-0.5,100.5,-1.1,1.1);
  TProfile* os_bbcs_c24_vs_bbcqs = new TProfile(Form("os_bbcs_c24_vs_bbcqs"),"",200,0,200,-1.1,1.1);
  TProfile* os_fvtxs_c24_vs_bbcqs = new TProfile(Form("os_fvtxs_c24_vs_bbcqs"),"",200,0,200,-1.1,1.1);
  TProfile* os_fvtxn_c24_vs_bbcqs = new TProfile(Form("os_fvtxn_c24_vs_bbcqs"),"",200,0,200,-1.1,1.1);
  TProfile* os_bbcs_c24_vs_nfvtxs = new TProfile(Form("os_bbcs_c24_vs_nfvtxs"),"",1000,-0.5,999.5,-1.1,1.1);
  TProfile* os_fvtxs_c24_vs_nfvtxs = new TProfile(Form("os_fvtxs_c24_vs_nfvtxs"),"",1000,-0.5,999.5,-1.1,1.1);
  TProfile* os_fvtxn_c24_vs_nfvtxs = new TProfile(Form("os_fvtxn_c24_vs_nfvtxs"),"",1000,-0.5,999.5,-1.1,1.1);
  TProfile* os_bbcs_d24_in_both = new TProfile(Form("os_bbcs_d24_in_both"),"",15,0.0,3.0,-1.1,1.1);
  TProfile* os_fvtxs_d24_in_both = new TProfile(Form("os_fvtxs_d24_in_both"),"",15,0.0,3.0,-1.1,1.1);
  TProfile* os_fvtxn_d24_in_both = new TProfile(Form("os_fvtxn_d24_in_both"),"",15,0.0,3.0,-1.1,1.1);
  TProfile* os_bbcs_d24_out_both = new TProfile(Form("os_bbcs_d24_out_both"),"",15,0.0,3.0,-1.1,1.1);
  TProfile* os_fvtxs_d24_out_both = new TProfile(Form("os_fvtxs_d24_out_both"),"",15,0.0,3.0,-1.1,1.1);
  TProfile* os_fvtxn_d24_out_both = new TProfile(Form("os_fvtxn_d24_out_both"),"",15,0.0,3.0,-1.1,1.1);

  TProfile* npc1_bbcs_v2 = new TProfile(Form("npc1_bbcs_v2"),"",61, -0.5, 60.5, -1.1, 1.1);
  TProfile* npc1_fvtxs_v2 = new TProfile(Form("npc1_fvtxs_v2"),"",61, -0.5, 60.5, -1.1, 1.1);
  TProfile* npc1_fvtxn_v2 = new TProfile(Form("npc1_fvtxn_v2"),"",61, -0.5, 60.5, -1.1, 1.1);
  //TProfile* npc1_os_cnt_c22 = new TProfile(Form("npc1_os_cnt_c22"),"",61, -0.5, 60.5, -1.1, 1.1);
  TProfile* npc1_os_bbcs_c22 = new TProfile(Form("npc1_os_bbcs_c22"),"",61, -0.5, 60.5, -1.1, 1.1);
  TProfile* npc1_os_fvtxs_c22 = new TProfile(Form("npc1_os_fvtxs_c22"),"",61, -0.5, 60.5, -1.1, 1.1);
  TProfile* npc1_os_fvtxn_c22 = new TProfile(Form("npc1_os_fvtxn_c22"),"",61, -0.5, 60.5, -1.1, 1.1);
  TProfile* npc1_os_fvtxs_tracks_c22 = new TProfile(Form("npc1_os_fvtxs_tracks_c22"),     "" ,61, -0.5, 60.5, -1.1, 1.1);
  TProfile* npc1_os_fvtxn_tracks_c22 = new TProfile(Form("npc1_os_fvtxn_tracks_c22"),     "" ,61, -0.5, 60.5, -1.1, 1.1);
  TProfile* npc1_os_cntbbcs_c22 = new TProfile(Form("npc1_os_cntbbcs_c22"),"",61, -0.5, 60.5, -1.1, 1.1);
  TProfile* npc1_os_cntfvtxs_c22 = new TProfile(Form("npc1_os_cntfvtxs_c22"),"",61, -0.5, 60.5, -1.1, 1.1);
  TProfile* npc1_os_cntfvtxn_c22 = new TProfile(Form("npc1_os_cntfvtxn_c22"),"",61, -0.5, 60.5, -1.1, 1.1);
  TProfile* npc1_os_bbcsfvtxs_c22 = new TProfile(Form("npc1_os_bbcsfvtxs_c22"),"",61, -0.5, 60.5, -1.1, 1.1);
  TProfile* npc1_os_bbcsfvtxn_c22 = new TProfile(Form("npc1_os_bbcsfvtxn_c22"),"",61, -0.5, 60.5, -1.1, 1.1);
  TProfile* npc1_os_fvtxsfvtxn_c22 = new TProfile(Form("npc1_os_fvtxsfvtxn_c22"),"",61, -0.5, 60.5, -1.1, 1.1);
  TProfile* npc1_os_fvtxsfvtxn_tracks_c22 = new TProfile(Form("npc1_os_fvtxsfvtxn_tracks_c22"),"",61, -0.5, 60.5, -1.1, 1.1);

  TProfile* npc1_bbcs_v3 = new TProfile(Form("npc1_bbcs_v3"),"",61, -0.5, 60.5, -1.1, 1.1);
  TProfile* npc1_fvtxs_v3 = new TProfile(Form("npc1_fvtxs_v3"),"",61, -0.5, 60.5, -1.1, 1.1);
  TProfile* npc1_fvtxn_v3 = new TProfile(Form("npc1_fvtxn_v3"),"",61, -0.5, 60.5, -1.1, 1.1);
  //TProfile* npc1_os_cnt_c32 = new TProfile(Form("npc1_os_cnt_c32"),"",61, -0.5, 60.5, -1.1, 1.1);
  TProfile* npc1_os_bbcs_c32 = new TProfile(Form("npc1_os_bbcs_c32"),"",61, -0.5, 60.5, -1.1, 1.1);
  TProfile* npc1_os_fvtxs_c32 = new TProfile(Form("npc1_os_fvtxs_c32"),"",61, -0.5, 60.5, -1.1, 1.1);
  TProfile* npc1_os_fvtxn_c32 = new TProfile(Form("npc1_os_fvtxn_c32"),"",61, -0.5, 60.5, -1.1, 1.1);
  TProfile* npc1_os_fvtxs_tracks_c32 = new TProfile(Form("npc1_os_fvtxs_tracks_c32"),     "" ,61, -0.5, 60.5, -1.1, 1.1);
  TProfile* npc1_os_fvtxn_tracks_c32 = new TProfile(Form("npc1_os_fvtxn_tracks_c32"),     "" ,61, -0.5, 60.5, -1.1, 1.1);
  TProfile* npc1_os_cntbbcs_c32 = new TProfile(Form("npc1_os_cntbbcs_c32"),"",61, -0.5, 60.5, -1.1, 1.1);
  TProfile* npc1_os_cntfvtxs_c32 = new TProfile(Form("npc1_os_cntfvtxs_c32"),"",61, -0.5, 60.5, -1.1, 1.1);
  TProfile* npc1_os_cntfvtxn_c32 = new TProfile(Form("npc1_os_cntfvtxn_c32"),"",61, -0.5, 60.5, -1.1, 1.1);
  TProfile* npc1_os_bbcsfvtxs_c32 = new TProfile(Form("npc1_os_bbcsfvtxs_c32"),"",61, -0.5, 60.5, -1.1, 1.1);
  TProfile* npc1_os_bbcsfvtxn_c32 = new TProfile(Form("npc1_os_bbcsfvtxn_c32"),"",61, -0.5, 60.5, -1.1, 1.1);
  TProfile* npc1_os_fvtxsfvtxn_c32 = new TProfile(Form("npc1_os_fvtxsfvtxn_c32"),"",61, -0.5, 60.5, -1.1, 1.1);
  TProfile* npc1_os_fvtxsfvtxn_tracks_c32 = new TProfile(Form("npc1_os_fvtxsfvtxn_tracks_c32"),"",61, -0.5, 60.5, -1.1, 1.1);

  TProfile* nfvtxt_bbcs_v2 = new TProfile(Form("nfvtxt_bbcs_v2"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_fvtxs_v2 = new TProfile(Form("nfvtxt_fvtxs_v2"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_fvtxn_v2 = new TProfile(Form("nfvtxt_fvtxn_v2"),"",80, -0.5, 79.5, -1.1, 1.1);
  //TProfile* nfvtxt_os_cnt_c22 = new TProfile(Form("nfvtxt_os_cnt_c22"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_os_bbcs_c22 = new TProfile(Form("nfvtxt_os_bbcs_c22"),"",80, -0.5, 79.5, -1.1, 1.1);
  //TProfile* nfvtxt_os_fvtxs_c22 = new TProfile(Form("nfvtxt_os_fvtxs_c22"),"",80, -0.5, 79.5, -1.1, 1.1);
  //TProfile* nfvtxt_os_fvtxn_c22 = new TProfile(Form("nfvtxt_os_fvtxn_c22"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_os_cntbbcs_c22 = new TProfile(Form("nfvtxt_os_cntbbcs_c22"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_os_cntfvtxs_c22 = new TProfile(Form("nfvtxt_os_cntfvtxs_c22"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_os_cntfvtxn_c22 = new TProfile(Form("nfvtxt_os_cntfvtxn_c22"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_os_bbcsfvtxs_c22 = new TProfile(Form("nfvtxt_os_bbcsfvtxs_c22"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_os_bbcsfvtxn_c22 = new TProfile(Form("nfvtxt_os_bbcsfvtxn_c22"),"",80, -0.5, 79.5, -1.1, 1.1);
  // ---
  TProfile* nfvtxt_bbcs_v3 = new TProfile(Form("nfvtxt_bbcs_v3"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_fvtxs_v3 = new TProfile(Form("nfvtxt_fvtxs_v3"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_fvtxn_v3 = new TProfile(Form("nfvtxt_fvtxn_v3"),"",80, -0.5, 79.5, -1.1, 1.1);
  //TProfile* nfvtxt_os_cnt_c32 = new TProfile(Form("nfvtxt_os_cnt_c32"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_os_bbcs_c32 = new TProfile(Form("nfvtxt_os_bbcs_c32"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_os_fvtxs_c32 = new TProfile(Form("nfvtxt_os_fvtxs_c32"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_os_fvtxn_c32 = new TProfile(Form("nfvtxt_os_fvtxn_c32"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_os_fvtxc_c32 = new TProfile(Form("nfvtxt_os_fvtxc_c32"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_os_fvtxa_c32 = new TProfile(Form("nfvtxt_os_fvtxa_c32"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_os_fvtxs_ce01_c32 = new TProfile(Form("nfvtxt_os_fvtxs_ce01_c32"),    "",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_os_fvtxn_ce01_c32 = new TProfile(Form("nfvtxt_os_fvtxn_ce01_c32"),    "",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_os_fvtxs_tracks_c32 = new TProfile(Form("nfvtxt_os_fvtxs_tracks_c32"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_os_fvtxn_tracks_c32 = new TProfile(Form("nfvtxt_os_fvtxn_tracks_c32"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_os_fvtxc_tracks_c32 = new TProfile(Form("nfvtxt_os_fvtxc_tracks_c32"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_os_fvtxa_tracks_c32 = new TProfile(Form("nfvtxt_os_fvtxa_tracks_c32"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_os_cntbbcs_c32 = new TProfile(Form("nfvtxt_os_cntbbcs_c32"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_os_cntfvtxs_c32 = new TProfile(Form("nfvtxt_os_cntfvtxs_c32"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_os_cntfvtxn_c32 = new TProfile(Form("nfvtxt_os_cntfvtxn_c32"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_os_bbcsfvtxs_c32 = new TProfile(Form("nfvtxt_os_bbcsfvtxs_c32"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_os_bbcsfvtxn_c32 = new TProfile(Form("nfvtxt_os_bbcsfvtxn_c32"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_os_fvtxsfvtxn_c32 = new TProfile(Form("nfvtxt_os_fvtxsfvtxn_c32"),"",80, -0.5, 79.5, -1.1, 1.1);
  // ---
  TProfile* nfvtxc_bbcs_v2 = new TProfile(Form("nfvtxc_bbcs_v2"),"",200, -0.5, 1999.5, -1.1, 1.1);
  TProfile* nfvtxc_fvtxs_v2 = new TProfile(Form("nfvtxc_fvtxs_v2"),"",200, -0.5, 1999.5, -1.1, 1.1);
  TProfile* nfvtxc_fvtxn_v2 = new TProfile(Form("nfvtxc_fvtxn_v2"),"",200, -0.5, 1999.5, -1.1, 1.1);
  //TProfile* nfvtxc_os_cnt_c22 = new TProfile(Form("nfvtxc_os_cnt_c22"),"",200, -0.5, 1999.5, -1.1, 1.1);
  TProfile* nfvtxc_os_bbcs_c22 = new TProfile(Form("nfvtxc_os_bbcs_c22"),"",200, -0.5, 1999.5, -1.1, 1.1);
  //TProfile* nfvtxc_os_fvtxs_c22 = new TProfile(Form("nfvtxc_os_fvtxs_c22"),"",200, -0.5, 1999.5, -1.1, 1.1);
  //TProfile* nfvtxc_os_fvtxn_c22 = new TProfile(Form("nfvtxc_os_fvtxn_c22"),"",200, -0.5, 1999.5, -1.1, 1.1);
  TProfile* nfvtxc_os_cntbbcs_c22 = new TProfile(Form("nfvtxc_os_cntbbcs_c22"),"",200, -0.5, 1999.5, -1.1, 1.1);
  TProfile* nfvtxc_os_cntfvtxs_c22 = new TProfile(Form("nfvtxc_os_cntfvtxs_c22"),"",200, -0.5, 1999.5, -1.1, 1.1);
  TProfile* nfvtxc_os_cntfvtxn_c22 = new TProfile(Form("nfvtxc_os_cntfvtxn_c22"),"",200, -0.5, 1999.5, -1.1, 1.1);
  TProfile* nfvtxc_os_bbcsfvtxs_c22 = new TProfile(Form("nfvtxc_os_bbcsfvtxs_c22"),"",200, -0.5, 1999.5, -1.1, 1.1);
  TProfile* nfvtxc_os_bbcsfvtxn_c22 = new TProfile(Form("nfvtxc_os_bbcsfvtxn_c22"),"",200, -0.5, 1999.5, -1.1, 1.1);
  // ---
  TProfile* nfvtxc_bbcs_v3 = new TProfile(Form("nfvtxc_bbcs_v3"),"",200, -0.5, 1999.5, -1.1, 1.1);
  TProfile* nfvtxc_fvtxs_v3 = new TProfile(Form("nfvtxc_fvtxs_v3"),"",200, -0.5, 1999.5, -1.1, 1.1);
  TProfile* nfvtxc_fvtxn_v3 = new TProfile(Form("nfvtxc_fvtxn_v3"),"",200, -0.5, 1999.5, -1.1, 1.1);
  //TProfile* nfvtxc_os_cnt_c32 = new TProfile(Form("nfvtxc_os_cnt_c32"),"",200, -0.5, 1999.5, -1.1, 1.1);
  TProfile* nfvtxc_os_bbcs_c32 = new TProfile(Form("nfvtxc_os_bbcs_c32"),"",200, -0.5, 1999.5, -1.1, 1.1);
  TProfile* nfvtxc_os_fvtxs_c32 = new TProfile(Form("nfvtxc_os_fvtxs_c32"),"",200, -0.5, 1999.5, -1.1, 1.1);
  TProfile* nfvtxc_os_fvtxn_c32 = new TProfile(Form("nfvtxc_os_fvtxn_c32"),"",200, -0.5, 1999.5, -1.1, 1.1);
  TProfile* nfvtxc_os_fvtxc_c32 = new TProfile(Form("nfvtxc_os_fvtxc_c32"),"",200, -0.5, 1999.5, -1.1, 1.1);
  TProfile* nfvtxc_os_fvtxa_c32 = new TProfile(Form("nfvtxc_os_fvtxa_c32"),"",200, -0.5, 1999.5, -1.1, 1.1);
  TProfile* nfvtxc_os_fvtxs_ce01_c32 = new TProfile(Form("nfvtxc_os_fvtxs_ce01_c32"),    "",200, -0.5, 1999.5, -1.1, 1.1);
  TProfile* nfvtxc_os_fvtxn_ce01_c32 = new TProfile(Form("nfvtxc_os_fvtxn_ce01_c32"),    "",200, -0.5, 1999.5, -1.1, 1.1);
  TProfile* nfvtxc_os_fvtxs_tracks_c32 = new TProfile(Form("nfvtxc_os_fvtxs_tracks_c32"),"",200, -0.5, 1999.5, -1.1, 1.1);
  TProfile* nfvtxc_os_fvtxn_tracks_c32 = new TProfile(Form("nfvtxc_os_fvtxn_tracks_c32"),"",200, -0.5, 1999.5, -1.1, 1.1);
  TProfile* nfvtxc_os_fvtxc_tracks_c32 = new TProfile(Form("nfvtxc_os_fvtxc_tracks_c32"),"",200, -0.5, 1999.5, -1.1, 1.1);
  TProfile* nfvtxc_os_fvtxa_tracks_c32 = new TProfile(Form("nfvtxc_os_fvtxa_tracks_c32"),"",200, -0.5, 1999.5, -1.1, 1.1);
  TProfile* nfvtxc_os_cntbbcs_c32 = new TProfile(Form("nfvtxc_os_cntbbcs_c32"),"",200, -0.5, 1999.5, -1.1, 1.1);
  TProfile* nfvtxc_os_cntfvtxs_c32 = new TProfile(Form("nfvtxc_os_cntfvtxs_c32"),"",200, -0.5, 1999.5, -1.1, 1.1);
  TProfile* nfvtxc_os_cntfvtxn_c32 = new TProfile(Form("nfvtxc_os_cntfvtxn_c32"),"",200, -0.5, 1999.5, -1.1, 1.1);
  TProfile* nfvtxc_os_bbcsfvtxs_c32 = new TProfile(Form("nfvtxc_os_bbcsfvtxs_c32"),"",200, -0.5, 1999.5, -1.1, 1.1);
  TProfile* nfvtxc_os_bbcsfvtxn_c32 = new TProfile(Form("nfvtxc_os_bbcsfvtxn_c32"),"",200, -0.5, 1999.5, -1.1, 1.1);
  TProfile* nfvtxc_os_fvtxsfvtxn_c32 = new TProfile(Form("nfvtxc_os_fvtxsfvtxn_c32"),"",200, -0.5, 1999.5, -1.1, 1.1);
  /// ---
  TProfile* nfvtxt_os_fvtxs_c22 = new TProfile(Form("nfvtxt_os_fvtxs_c22"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_os_fvtxn_c22 = new TProfile(Form("nfvtxt_os_fvtxn_c22"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_os_fvtxc_c22 = new TProfile(Form("nfvtxt_os_fvtxc_c22"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_os_fvtxa_c22 = new TProfile(Form("nfvtxt_os_fvtxa_c22"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_os_fvtxs_c24 = new TProfile(Form("nfvtxt_os_fvtxs_c24"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_os_fvtxn_c24 = new TProfile(Form("nfvtxt_os_fvtxn_c24"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_os_fvtxc_c24 = new TProfile(Form("nfvtxt_os_fvtxc_c24"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_os_fvtxa_c24 = new TProfile(Form("nfvtxt_os_fvtxa_c24"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_os_fvtxs_tracks_c22 = new TProfile(Form("nfvtxt_os_fvtxs_tracks_c22"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_os_fvtxn_tracks_c22 = new TProfile(Form("nfvtxt_os_fvtxn_tracks_c22"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_os_fvtxc_tracks_c22 = new TProfile(Form("nfvtxt_os_fvtxc_tracks_c22"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_os_fvtxa_tracks_c22 = new TProfile(Form("nfvtxt_os_fvtxa_tracks_c22"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_os_fvtxs_tracks_c24 = new TProfile(Form("nfvtxt_os_fvtxs_tracks_c24"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_os_fvtxn_tracks_c24 = new TProfile(Form("nfvtxt_os_fvtxn_tracks_c24"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_os_fvtxc_tracks_c24 = new TProfile(Form("nfvtxt_os_fvtxc_tracks_c24"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_os_fvtxa_tracks_c24 = new TProfile(Form("nfvtxt_os_fvtxa_tracks_c24"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_os_fvtxsfvtxn_c22  = new TProfile(Form("nfvtxt_os_fvtxsfvtxn_c22"), "",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_os_fvtxsfvtxn_c24  = new TProfile(Form("nfvtxt_os_fvtxsfvtxn_c24"), "",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_os_fvtxsfvtxn_c24a = new TProfile(Form("nfvtxt_os_fvtxsfvtxn_c24a"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_os_fvtxsfvtxn_c24b = new TProfile(Form("nfvtxt_os_fvtxsfvtxn_c24b"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_os_fvtxsfvtxn_c24c = new TProfile(Form("nfvtxt_os_fvtxsfvtxn_c24c"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_os_fvtxsfvtxn_c24d = new TProfile(Form("nfvtxt_os_fvtxsfvtxn_c24d"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_os_fvtxsfvtxn_tracks_c22  = new TProfile(Form("nfvtxt_os_fvtxsfvtxn_tracks_c22"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_os_fvtxsfvtxn_tracks_c24  = new TProfile(Form("nfvtxt_os_fvtxsfvtxn_tracks_c24"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_os_fvtxsfvtxn_tracks_c24a = new TProfile(Form("nfvtxt_os_fvtxsfvtxn_tracks_c24a"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_os_fvtxsfvtxn_tracks_c24b = new TProfile(Form("nfvtxt_os_fvtxsfvtxn_tracks_c24b"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_os_fvtxsfvtxn_tracks_c24c = new TProfile(Form("nfvtxt_os_fvtxsfvtxn_tracks_c24c"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_os_fvtxsfvtxn_tracks_c24d = new TProfile(Form("nfvtxt_os_fvtxsfvtxn_tracks_c24d"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_os_fvtxsfvtxn_tracks_c32  = new TProfile(Form("nfvtxt_os_fvtxsfvtxn_tracks_c32"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_os_fvtxs_ce01_c22 = new TProfile(Form("nfvtxt_os_fvtxs_ce01_c22"),       "",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_os_fvtxn_ce01_c22 = new TProfile(Form("nfvtxt_os_fvtxn_ce01_c22"),       "",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_os_fvtxsfvtxn_ce01_c24  = new TProfile(Form("nfvtxt_os_fvtxsfvtxn_ce01_c24"), "",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_os_fvtxsfvtxn_ce01_c24a = new TProfile(Form("nfvtxt_os_fvtxsfvtxn_ce01_c24a"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_os_fvtxsfvtxn_ce01_c24b = new TProfile(Form("nfvtxt_os_fvtxsfvtxn_ce01_c24b"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_os_fvtxsfvtxn_ce01_c24c = new TProfile(Form("nfvtxt_os_fvtxsfvtxn_ce01_c24c"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_os_fvtxsfvtxn_ce01_c24d = new TProfile(Form("nfvtxt_os_fvtxsfvtxn_ce01_c24d"),"",80, -0.5, 79.5, -1.1, 1.1);
  /// ---
  TProfile* nfvtxc_os_fvtxs_c22 = new TProfile(Form("nfvtxc_os_fvtxs_c22"),"",200, -0.5, 1999.5, -1.1, 1.1);
  TProfile* nfvtxc_os_fvtxn_c22 = new TProfile(Form("nfvtxc_os_fvtxn_c22"),"",200, -0.5, 1999.5, -1.1, 1.1);
  TProfile* nfvtxc_os_fvtxc_c22 = new TProfile(Form("nfvtxc_os_fvtxc_c22"),"",200, -0.5, 1999.5, -1.1, 1.1);
  TProfile* nfvtxc_os_fvtxa_c22 = new TProfile(Form("nfvtxc_os_fvtxa_c22"),"",200, -0.5, 1999.5, -1.1, 1.1);
  TProfile* nfvtxc_os_fvtxs_c24 = new TProfile(Form("nfvtxc_os_fvtxs_c24"),"",200, -0.5, 1999.5, -1.1, 1.1);
  TProfile* nfvtxc_os_fvtxn_c24 = new TProfile(Form("nfvtxc_os_fvtxn_c24"),"",200, -0.5, 1999.5, -1.1, 1.1);
  TProfile* nfvtxc_os_fvtxc_c24 = new TProfile(Form("nfvtxc_os_fvtxc_c24"),"",200, -0.5, 1999.5, -1.1, 1.1);
  TProfile* nfvtxc_os_fvtxa_c24 = new TProfile(Form("nfvtxc_os_fvtxa_c24"),"",200, -0.5, 1999.5, -1.1, 1.1);
  TProfile* nfvtxc_os_fvtxs_tracks_c22 = new TProfile(Form("nfvtxc_os_fvtxs_tracks_c22"),"",200, -0.5, 1999.5, -1.1, 1.1);
  TProfile* nfvtxc_os_fvtxn_tracks_c22 = new TProfile(Form("nfvtxc_os_fvtxn_tracks_c22"),"",200, -0.5, 1999.5, -1.1, 1.1);
  TProfile* nfvtxc_os_fvtxc_tracks_c22 = new TProfile(Form("nfvtxc_os_fvtxc_tracks_c22"),"",200, -0.5, 1999.5, -1.1, 1.1);
  TProfile* nfvtxc_os_fvtxa_tracks_c22 = new TProfile(Form("nfvtxc_os_fvtxa_tracks_c22"),"",200, -0.5, 1999.5, -1.1, 1.1);
  TProfile* nfvtxc_os_fvtxs_tracks_c24 = new TProfile(Form("nfvtxc_os_fvtxs_tracks_c24"),"",200, -0.5, 1999.5, -1.1, 1.1);
  TProfile* nfvtxc_os_fvtxn_tracks_c24 = new TProfile(Form("nfvtxc_os_fvtxn_tracks_c24"),"",200, -0.5, 1999.5, -1.1, 1.1);
  TProfile* nfvtxc_os_fvtxc_tracks_c24 = new TProfile(Form("nfvtxc_os_fvtxc_tracks_c24"),"",200, -0.5, 1999.5, -1.1, 1.1);
  TProfile* nfvtxc_os_fvtxa_tracks_c24 = new TProfile(Form("nfvtxc_os_fvtxa_tracks_c24"),"",200, -0.5, 1999.5, -1.1, 1.1);
  TProfile* nfvtxc_os_fvtxsfvtxn_c22  = new TProfile(Form("nfvtxc_os_fvtxsfvtxn_c22"), "",200, -0.5, 1999.5, -1.1, 1.1);
  TProfile* nfvtxc_os_fvtxsfvtxn_c24  = new TProfile(Form("nfvtxc_os_fvtxsfvtxn_c24"), "",200, -0.5, 1999.5, -1.1, 1.1);
  TProfile* nfvtxc_os_fvtxsfvtxn_c24a = new TProfile(Form("nfvtxc_os_fvtxsfvtxn_c24a"),"",200, -0.5, 1999.5, -1.1, 1.1);
  TProfile* nfvtxc_os_fvtxsfvtxn_c24b = new TProfile(Form("nfvtxc_os_fvtxsfvtxn_c24b"),"",200, -0.5, 1999.5, -1.1, 1.1);
  TProfile* nfvtxc_os_fvtxsfvtxn_c24c = new TProfile(Form("nfvtxc_os_fvtxsfvtxn_c24c"),"",200, -0.5, 1999.5, -1.1, 1.1);
  TProfile* nfvtxc_os_fvtxsfvtxn_c24d = new TProfile(Form("nfvtxc_os_fvtxsfvtxn_c24d"),"",200, -0.5, 1999.5, -1.1, 1.1);
  TProfile* nfvtxc_os_fvtxsfvtxn_tracks_c22  = new TProfile(Form("nfvtxc_os_fvtxsfvtxn_tracks_c22"), "",200, -0.5, 1999.5, -1.1, 1.1);
  TProfile* nfvtxc_os_fvtxsfvtxn_tracks_c24  = new TProfile(Form("nfvtxc_os_fvtxsfvtxn_tracks_c24"), "",200, -0.5, 1999.5, -1.1, 1.1);
  TProfile* nfvtxc_os_fvtxsfvtxn_tracks_c24a = new TProfile(Form("nfvtxc_os_fvtxsfvtxn_tracks_c24a"),"",200, -0.5, 1999.5, -1.1, 1.1);
  TProfile* nfvtxc_os_fvtxsfvtxn_tracks_c24b = new TProfile(Form("nfvtxc_os_fvtxsfvtxn_tracks_c24b"),"",200, -0.5, 1999.5, -1.1, 1.1);
  TProfile* nfvtxc_os_fvtxsfvtxn_tracks_c24c = new TProfile(Form("nfvtxc_os_fvtxsfvtxn_tracks_c24c"),"",200, -0.5, 1999.5, -1.1, 1.1);
  TProfile* nfvtxc_os_fvtxsfvtxn_tracks_c24d = new TProfile(Form("nfvtxc_os_fvtxsfvtxn_tracks_c24d"),"",200, -0.5, 1999.5, -1.1, 1.1);
  TProfile* nfvtxc_os_fvtxsfvtxn_tracks_c32  = new TProfile(Form("nfvtxc_os_fvtxsfvtxn_tracks_c32"), "",200, -0.5, 1999.5, -1.1, 1.1);
  TProfile* nfvtxc_os_fvtxs_ce01_c22 = new TProfile(Form("nfvtxc_os_fvtxs_ce01_c22"),"",200, -0.5, 1999.5, -1.1, 1.1);
  TProfile* nfvtxc_os_fvtxn_ce01_c22 = new TProfile(Form("nfvtxc_os_fvtxn_ce01_c22"),"",200, -0.5, 1999.5, -1.1, 1.1);
  TProfile* nfvtxc_os_fvtxsfvtxn_ce01_c24  = new TProfile(Form("nfvtxc_os_fvtxsfvtxn_ce01_c24"), "",200, -0.5, 1999.5, -1.1, 1.1);
  TProfile* nfvtxc_os_fvtxsfvtxn_ce01_c24a = new TProfile(Form("nfvtxc_os_fvtxsfvtxn_ce01_c24a"),"",200, -0.5, 1999.5, -1.1, 1.1);
  TProfile* nfvtxc_os_fvtxsfvtxn_ce01_c24b = new TProfile(Form("nfvtxc_os_fvtxsfvtxn_ce01_c24b"),"",200, -0.5, 1999.5, -1.1, 1.1);
  TProfile* nfvtxc_os_fvtxsfvtxn_ce01_c24c = new TProfile(Form("nfvtxc_os_fvtxsfvtxn_ce01_c24c"),"",200, -0.5, 1999.5, -1.1, 1.1);
  TProfile* nfvtxc_os_fvtxsfvtxn_ce01_c24d = new TProfile(Form("nfvtxc_os_fvtxsfvtxn_ce01_c24d"),"",200, -0.5, 1999.5, -1.1, 1.1);


  // --- now we'll do some pseudorapidity depenence
  TProfile* bbcs_v2eta_west_docalib = new TProfile(Form("bbcs_v2eta_west_docalib"),"",32, -3.2, 3.2, -1.1, 1.1);
  TProfile* bbcs_v2eta_east_docalib = new TProfile(Form("bbcs_v2eta_east_docalib"),"",32, -3.2, 3.2, -1.1, 1.1);
  TProfile* bbcs_v2eta_both_docalib = new TProfile(Form("bbcs_v2eta_both_docalib"),"",32, -3.2, 3.2, -1.1, 1.1);
  TProfile* fvtxs_v2eta_west_docalib = new TProfile(Form("fvtxs_v2eta_west_docalib"),"",32, -3.2, 3.2, -1.1, 1.1);
  TProfile* fvtxs_v2eta_east_docalib = new TProfile(Form("fvtxs_v2eta_east_docalib"),"",32, -3.2, 3.2, -1.1, 1.1);
  TProfile* fvtxs_v2eta_both_docalib = new TProfile(Form("fvtxs_v2eta_both_docalib"),"",32, -3.2, 3.2, -1.1, 1.1);
  TProfile* fvtxn_v2eta_west_docalib = new TProfile(Form("fvtxn_v2eta_west_docalib"),"",32, -3.2, 3.2, -1.1, 1.1);
  TProfile* fvtxn_v2eta_east_docalib = new TProfile(Form("fvtxn_v2eta_east_docalib"),"",32, -3.2, 3.2, -1.1, 1.1);
  TProfile* fvtxn_v2eta_both_docalib = new TProfile(Form("fvtxn_v2eta_both_docalib"),"",32, -3.2, 3.2, -1.1, 1.1);
  TProfile* bbcs_v3eta_west_docalib = new TProfile(Form("bbcs_v3eta_west_docalib"),"",32, -3.2, 3.2, -1.1, 1.1);
  TProfile* bbcs_v3eta_east_docalib = new TProfile(Form("bbcs_v3eta_east_docalib"),"",32, -3.2, 3.2, -1.1, 1.1);
  TProfile* bbcs_v3eta_both_docalib = new TProfile(Form("bbcs_v3eta_both_docalib"),"",32, -3.2, 3.2, -1.1, 1.1);
  TProfile* fvtxs_v3eta_west_docalib = new TProfile(Form("fvtxs_v3eta_west_docalib"),"",32, -3.2, 3.2, -1.1, 1.1);
  TProfile* fvtxs_v3eta_east_docalib = new TProfile(Form("fvtxs_v3eta_east_docalib"),"",32, -3.2, 3.2, -1.1, 1.1);
  TProfile* fvtxs_v3eta_both_docalib = new TProfile(Form("fvtxs_v3eta_both_docalib"),"",32, -3.2, 3.2, -1.1, 1.1);
  TProfile* fvtxn_v3eta_west_docalib = new TProfile(Form("fvtxn_v3eta_west_docalib"),"",32, -3.2, 3.2, -1.1, 1.1);
  TProfile* fvtxn_v3eta_east_docalib = new TProfile(Form("fvtxn_v3eta_east_docalib"),"",32, -3.2, 3.2, -1.1, 1.1);
  TProfile* fvtxn_v3eta_both_docalib = new TProfile(Form("fvtxn_v3eta_both_docalib"),"",32, -3.2, 3.2, -1.1, 1.1);
  // ---
  TProfile* os_bbcs_d22eta_west = new TProfile(Form("os_bbcs_d22eta_west"),"",32, -3.2, 3.2, -1.1, 1.1);
  TProfile* os_bbcs_d22eta_east = new TProfile(Form("os_bbcs_d22eta_east"),"",32, -3.2, 3.2, -1.1, 1.1);
  TProfile* os_bbcs_d22eta_both = new TProfile(Form("os_bbcs_d22eta_both"),"",32, -3.2, 3.2, -1.1, 1.1);
  TProfile* os_fvtxs_d22eta_west = new TProfile(Form("os_fvtxs_d22eta_west"),"",32, -3.2, 3.2, -1.1, 1.1);
  TProfile* os_fvtxs_d22eta_east = new TProfile(Form("os_fvtxs_d22eta_east"),"",32, -3.2, 3.2, -1.1, 1.1);
  TProfile* os_fvtxs_d22eta_both = new TProfile(Form("os_fvtxs_d22eta_both"),"",32, -3.2, 3.2, -1.1, 1.1);
  TProfile* os_fvtxn_d22eta_west = new TProfile(Form("os_fvtxn_d22eta_west"),"",32, -3.2, 3.2, -1.1, 1.1);
  TProfile* os_fvtxn_d22eta_east = new TProfile(Form("os_fvtxn_d22eta_east"),"",32, -3.2, 3.2, -1.1, 1.1);
  TProfile* os_fvtxn_d22eta_both = new TProfile(Form("os_fvtxn_d22eta_both"),"",32, -3.2, 3.2, -1.1, 1.1);
  TProfile* os_fvtxs_tracks_d22eta_west = new TProfile(Form("os_fvtxs_tracks_d22eta_west"),"",32, -3.2, 3.2, -1.1, 1.1);
  TProfile* os_fvtxs_tracks_d22eta_east = new TProfile(Form("os_fvtxs_tracks_d22eta_east"),"",32, -3.2, 3.2, -1.1, 1.1);
  TProfile* os_fvtxs_tracks_d22eta_both = new TProfile(Form("os_fvtxs_tracks_d22eta_both"),"",32, -3.2, 3.2, -1.1, 1.1);
  TProfile* os_fvtxn_tracks_d22eta_west = new TProfile(Form("os_fvtxn_tracks_d22eta_west"),"",32, -3.2, 3.2, -1.1, 1.1);
  TProfile* os_fvtxn_tracks_d22eta_east = new TProfile(Form("os_fvtxn_tracks_d22eta_east"),"",32, -3.2, 3.2, -1.1, 1.1);
  TProfile* os_fvtxn_tracks_d22eta_both = new TProfile(Form("os_fvtxn_tracks_d22eta_both"),"",32, -3.2, 3.2, -1.1, 1.1);
  TProfile* os_bbcs_d32eta_west = new TProfile(Form("os_bbcs_d32eta_west"),"",32, -3.2, 3.2, -1.1, 1.1);
  TProfile* os_bbcs_d32eta_east = new TProfile(Form("os_bbcs_d32eta_east"),"",32, -3.2, 3.2, -1.1, 1.1);
  TProfile* os_bbcs_d32eta_both = new TProfile(Form("os_bbcs_d32eta_both"),"",32, -3.2, 3.2, -1.1, 1.1);
  TProfile* os_fvtxs_d32eta_west = new TProfile(Form("os_fvtxs_d32eta_west"),"",32, -3.2, 3.2, -1.1, 1.1);
  TProfile* os_fvtxs_d32eta_east = new TProfile(Form("os_fvtxs_d32eta_east"),"",32, -3.2, 3.2, -1.1, 1.1);
  TProfile* os_fvtxs_d32eta_both = new TProfile(Form("os_fvtxs_d32eta_both"),"",32, -3.2, 3.2, -1.1, 1.1);
  TProfile* os_fvtxn_d32eta_west = new TProfile(Form("os_fvtxn_d32eta_west"),"",32, -3.2, 3.2, -1.1, 1.1);
  TProfile* os_fvtxn_d32eta_east = new TProfile(Form("os_fvtxn_d32eta_east"),"",32, -3.2, 3.2, -1.1, 1.1);
  TProfile* os_fvtxn_d32eta_both = new TProfile(Form("os_fvtxn_d32eta_both"),"",32, -3.2, 3.2, -1.1, 1.1);
  TProfile* os_fvtxs_tracks_d32eta_west = new TProfile(Form("os_fvtxs_tracks_d32eta_west"),"",32, -3.2, 3.2, -1.1, 1.1);
  TProfile* os_fvtxs_tracks_d32eta_east = new TProfile(Form("os_fvtxs_tracks_d32eta_east"),"",32, -3.2, 3.2, -1.1, 1.1);
  TProfile* os_fvtxs_tracks_d32eta_both = new TProfile(Form("os_fvtxs_tracks_d32eta_both"),"",32, -3.2, 3.2, -1.1, 1.1);
  TProfile* os_fvtxn_tracks_d32eta_west = new TProfile(Form("os_fvtxn_tracks_d32eta_west"),"",32, -3.2, 3.2, -1.1, 1.1);
  TProfile* os_fvtxn_tracks_d32eta_east = new TProfile(Form("os_fvtxn_tracks_d32eta_east"),"",32, -3.2, 3.2, -1.1, 1.1);
  TProfile* os_fvtxn_tracks_d32eta_both = new TProfile(Form("os_fvtxn_tracks_d32eta_both"),"",32, -3.2, 3.2, -1.1, 1.1);

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

  int nfvtxt;
  float feta[75];
  float fphi[75];
  float fchisq[75];
  float fdcax[75];
  float fdcay[75];

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
  TBranch* b_nfvtxt;   //!
  TBranch* b_fphi;   //!
  TBranch* b_feta;   //!
  TBranch* b_fchisq;   //!
  TBranch* b_fdcax;   //!
  TBranch* b_fdcay;   //!
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

  ntp_event_chain->SetBranchAddress("ntracklets",&nfvtxt,&b_nfvtxt);
  ntp_event_chain->SetBranchAddress("fphi",fphi,&b_fphi);
  ntp_event_chain->SetBranchAddress("feta",feta,&b_feta);
  ntp_event_chain->SetBranchAddress("fchisq",fchisq,&b_fchisq);
  ntp_event_chain->SetBranchAddress("fDCA_X",fdcax,&b_fdcax);
  ntp_event_chain->SetBranchAddress("fDCA_Y",fdcay,&b_fdcay);



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
          if ( energyflag == 39  )
	    {
	      if ( centrality <= 20 ) icent = 0;
	      else if ( centrality <= 40 ) icent = 1;
	      else if ( centrality <= 60 ) icent = 2;
	      else if ( centrality <= 88 ) icent = 3;
	    }
        }

      double ZVTX = -9999;
      if ( runNumber >= 454774 && runNumber <= 456283 ) ZVTX = d_bbcz;
      if ( runNumber >= 456652 && runNumber <= 458167 ) ZVTX = eventfvtx_z;
      if ( fabs(ZVTX) > 10.0 )
        {
          if ( verbosity > 1 ) cout << "vertex rejected" << endl;
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

      int nfvtxc = d_nFVTX_clus;
      int nfvtxc_south = d_nFVTXN_clus;
      int nfvtxc_north = d_nFVTXS_clus;

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

      float bbc_nw_qx2 = 0;
      float bbc_nw_qy2 = 0;
      float bbc_nw_qx3 = 0;
      float bbc_nw_qy3 = 0;
      float bbc_nw_qx4 = 0;
      float bbc_nw_qy4 = 0;
      float bbc_nw_qw = 0;

      if ( ( say_event && verbosity > 0 ) || verbosity > 1 ) cout << "Looping over BBC stuff now" << endl;

      if ( bbc_pmts )
        {
          for(int ipmt = 0; ipmt < 64; ipmt++)
            {
              float bbc_charge = d_BBC_charge[ipmt];
              if ( bbc_charge <= 0 ) continue;

              float bbc_x      = d_pmt_x[ipmt] - vtx_x*10;//pmt location in mm
              float bbc_y      = d_pmt_y[ipmt] - vtx_y*10;
              float bbc_z      = d_pmt_z       - vtx_z*10;

              // --- rotation
              bbc_x = bbc_z*sin(-beam_angle) + bbc_x*cos(-beam_angle);

              float phi = TMath::ATan2(bbc_y,bbc_x);

              int ring = get_pmt_layer(ipmt);
              float tube = ipmt;
              float bbc_charge_corrected = bbc_charge * tube_gaincorrection[ipmt]; // tube by tube gain correction from above
              if ( bbc_charge_corrected <= 0 ) continue;

              tp1f_bbc_charge_phi->Fill(phi,bbc_charge);
              if ( ring == 0 ) tp1f_bbc0_charge_phi->Fill(phi,bbc_charge);
              if ( ring == 1 ) tp1f_bbc1_charge_phi->Fill(phi,bbc_charge);
              if ( ring == 2 ) tp1f_bbc2_charge_phi->Fill(phi,bbc_charge);
              if ( ring == 3 ) tp1f_bbc3_charge_phi->Fill(phi,bbc_charge);
              if ( ring == 4 ) tp1f_bbc4_charge_phi->Fill(phi,bbc_charge);

              tp1f_bbc_charge_tube->Fill(tube,bbc_charge);
              if ( ring == 0 ) tp1f_bbc0_charge_tube->Fill(tube,bbc_charge);
              if ( ring == 1 ) tp1f_bbc1_charge_tube->Fill(tube,bbc_charge);
              if ( ring == 2 ) tp1f_bbc2_charge_tube->Fill(tube,bbc_charge);
              if ( ring == 3 ) tp1f_bbc3_charge_tube->Fill(tube,bbc_charge);
              if ( ring == 4 ) tp1f_bbc4_charge_tube->Fill(tube,bbc_charge);

              tp1f_bbc_charge_wphi->Fill(phi,bbc_charge_corrected);
              if ( ring == 0 ) tp1f_bbc0_charge_wphi->Fill(phi,bbc_charge_corrected);
              if ( ring == 1 ) tp1f_bbc1_charge_wphi->Fill(phi,bbc_charge_corrected);
              if ( ring == 2 ) tp1f_bbc2_charge_wphi->Fill(phi,bbc_charge_corrected);
              if ( ring == 3 ) tp1f_bbc3_charge_wphi->Fill(phi,bbc_charge_corrected);
              if ( ring == 4 ) tp1f_bbc4_charge_wphi->Fill(phi,bbc_charge_corrected);

              tp1f_bbc_charge_wtube->Fill(tube,bbc_charge_corrected);
              if ( ring == 0 ) tp1f_bbc0_charge_wtube->Fill(tube,bbc_charge_corrected);
              if ( ring == 1 ) tp1f_bbc1_charge_wtube->Fill(tube,bbc_charge_corrected);
              if ( ring == 2 ) tp1f_bbc2_charge_wtube->Fill(tube,bbc_charge_corrected);
              if ( ring == 3 ) tp1f_bbc3_charge_wtube->Fill(tube,bbc_charge_corrected);
              if ( ring == 4 ) tp1f_bbc4_charge_wtube->Fill(tube,bbc_charge_corrected);

              bbc_qx2 += bbc_charge_corrected*TMath::Cos(2*phi);
              bbc_qy2 += bbc_charge_corrected*TMath::Sin(2*phi);
              bbc_qx3 += bbc_charge_corrected*TMath::Cos(3*phi);
              bbc_qy3 += bbc_charge_corrected*TMath::Sin(3*phi);
              bbc_qx4 += bbc_charge_corrected*TMath::Cos(4*phi);
              bbc_qy4 += bbc_charge_corrected*TMath::Sin(4*phi);
              bbc_qw += bbc_charge_corrected;

              bbc_nw_qx2 += bbc_charge*TMath::Cos(2*phi);
              bbc_nw_qy2 += bbc_charge*TMath::Sin(2*phi);
              bbc_nw_qx3 += bbc_charge*TMath::Cos(3*phi);
              bbc_nw_qy3 += bbc_charge*TMath::Sin(3*phi);
              bbc_nw_qx4 += bbc_charge*TMath::Cos(4*phi);
              bbc_nw_qy4 += bbc_charge*TMath::Sin(4*phi);
              bbc_nw_qw += bbc_charge;
            } // loop over tubes
        } // check on tubes

      // --- do centrality cut here!!!

      if ( centrality < 0 )
        {
          //cout << "centrality undefined, cutting on bbc charge" << endl;
          // --- revise these numbers as needed
          if ( energyflag == 200 && bbc_nw_qw < 60.0 ) continue;
          if ( energyflag == 62  && bbc_nw_qw < 40.0 ) continue;
          if ( energyflag == 20  && bbc_nw_qw < 25.0 ) continue;
          if ( energyflag == 39  && bbc_nw_qw < 30.0 ) continue;
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

      th1d_BBC_charge->Fill(bbc_nw_qw);
      th1d_FVTX_nclus->Fill(d_nFVTX_clus);
      th2d_qBBC_nFVTX->Fill(bbc_nw_qw,d_nFVTX_clus);

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

      for(int ilayer = 0; ilayer < 10; ilayer++)
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

      for(int ilayer = 0; ilayer < 10; ilayer++)
        {
          fvtxn_qx2[ilayer] = 0.0;
          fvtxn_qy2[ilayer] = 0.0;
          fvtxn_qx3[ilayer] = 0.0;
          fvtxn_qy3[ilayer] = 0.0;
          fvtxn_qx4[ilayer] = 0.0;
          fvtxn_qy4[ilayer] = 0.0;
          fvtxn_qw[ilayer] = 0.0;
        } // loop over layers

      int nclus_south_inner = 0; // good_4_event
      int nclus_north_inner = 0;
      int nclus_south_outer = 0;
      int nclus_north_outer = 0;

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

	      // ------------------------------------------------
	      // --- determine the weights for south clusters ---
              float fvtx_weight = 1.0;
              if ( doweights )
                {
                  if ( !th2d_fvtxs_phi_weight[fvtx_layer+1] )
                    {
                      cout << "WARNING!!!  Problem with weight histograms in cluster loop..." << endl;
                      continue;
                    }
                  int phi_bin = th2d_fvtxs_phi_weight[fvtx_layer+1]->GetYaxis()->FindBin(phi); // COME BACK HERE AND HAVE A LOOK
                  int zvtx_bin = th2d_fvtxs_phi_weight[fvtx_layer+1]->GetXaxis()->FindBin(ZVTX); // COME BACK HERE AND HAVE A LOOK
                  fvtx_weight = th2d_fvtxs_phi_weight[fvtx_layer+1]->GetBinContent(zvtx_bin,phi_bin);
                  // ---
                  if ( !th1d_fvtxs_phi_weight[fvtx_layer+1] )
                    {
                      cout << "WARNING!!!  Problem with weight histograms in cluster loop..." << endl;
                      continue;
                    }
                  phi_bin = th1d_fvtxs_phi_weight[fvtx_layer+1]->FindBin(phi); // COME BACK HERE AND HAVE A LOOK
                  fvtx_weight = th1d_fvtxs_phi_weight[fvtx_layer+1]->GetBinContent(phi_bin);
                }
              if ( fvtx_weight != fvtx_weight ) fvtx_weight = 0;
              if ( fvtx_weight < 0 ) fvtx_weight = 0;
              if ( fvtx_weight > 10 ) fvtx_weight = 0;

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

		  // --- inner eta clusters
		  if ( fvtx_eta > -2.0 && fvtx_eta < -0.5 )
		    {
		      fvtxs_qx2[6] += fvtx_weight * TMath::Cos(2*phi);
		      fvtxs_qy2[6] += fvtx_weight * TMath::Sin(2*phi);
		      fvtxs_qx3[6] += fvtx_weight * TMath::Cos(3*phi);
		      fvtxs_qy3[6] += fvtx_weight * TMath::Sin(3*phi);
		      fvtxs_qx4[6] += fvtx_weight * TMath::Cos(4*phi);
		      fvtxs_qy4[6] += fvtx_weight * TMath::Sin(4*phi);
		      fvtxs_qw[6] += fvtx_weight;
		      ++nclus_south_inner; // good_4_event
		    }
		  // --- outer eta clusters
		  if ( fvtx_eta > -3.5 && fvtx_eta < -2.0 )
		    {
		      fvtxs_qx2[7] += fvtx_weight * TMath::Cos(2*phi);
		      fvtxs_qy2[7] += fvtx_weight * TMath::Sin(2*phi);
		      fvtxs_qx3[7] += fvtx_weight * TMath::Cos(3*phi);
		      fvtxs_qy3[7] += fvtx_weight * TMath::Sin(3*phi);
		      fvtxs_qx4[7] += fvtx_weight * TMath::Cos(4*phi);
		      fvtxs_qy4[7] += fvtx_weight * TMath::Sin(4*phi);
		      fvtxs_qw[7] += fvtx_weight;
		      ++nclus_south_outer; // good_4_event
		    }

                  // ------------------------------------------
                  // --- explicitly excluding weights ---------
                  fvtxs_qx2[5] += TMath::Cos(2*phi);
                  fvtxs_qy2[5] += TMath::Sin(2*phi);
                  fvtxs_qx3[5] += TMath::Cos(3*phi);
                  fvtxs_qy3[5] += TMath::Sin(3*phi);
                  fvtxs_qx4[5] += TMath::Cos(4*phi);
                  fvtxs_qy4[5] += TMath::Sin(4*phi);
                  fvtxs_qw[5] += 1;
                  // -----------------------------------------

                  th1d_fvtxs_clus_phi->Fill(phi);
                  if ( fvtx_layer == 0 ) th1d_fvtxs0_clus_phi->Fill(phi);
                  if ( fvtx_layer == 1 ) th1d_fvtxs1_clus_phi->Fill(phi);
                  if ( fvtx_layer == 2 ) th1d_fvtxs2_clus_phi->Fill(phi);
                  if ( fvtx_layer == 3 ) th1d_fvtxs3_clus_phi->Fill(phi);
                  th1d_fvtxs_clus_wphi->Fill(phi,fvtx_weight);
                  if ( fvtx_layer == 0 ) th1d_fvtxs0_clus_wphi->Fill(phi,fvtx_weight);
                  if ( fvtx_layer == 1 ) th1d_fvtxs1_clus_wphi->Fill(phi,fvtx_weight);
                  if ( fvtx_layer == 2 ) th1d_fvtxs2_clus_wphi->Fill(phi,fvtx_weight);
                  if ( fvtx_layer == 3 ) th1d_fvtxs3_clus_wphi->Fill(phi,fvtx_weight);
                } // check on south

	      // ------------------------------------------------
	      // --- determine the weights for north clusters ---
              if ( doweights )
                {
                  if ( !th2d_fvtxn_phi_weight[fvtx_layer+1] )
                    {
                      cout << "WARNING!!!  Problem with weight histograms in cluster loop..." << endl;
                      continue;
                    }
                  int phi_bin = th2d_fvtxn_phi_weight[fvtx_layer+1]->GetYaxis()->FindBin(phi); // COME BACK HERE AND HAVE A LOOK
                  int zvtx_bin = th2d_fvtxn_phi_weight[fvtx_layer+1]->GetXaxis()->FindBin(ZVTX); // COME BACK HERE AND HAVE A LOOK
                  fvtx_weight = th2d_fvtxn_phi_weight[fvtx_layer+1]->GetBinContent(zvtx_bin,phi_bin);
                  // ---
                  if ( !th1d_fvtxn_phi_weight[fvtx_layer+1] )
                    {
                      cout << "WARNING!!!  Problem with weight histograms in cluster loop..." << endl;
                      continue;
                    }
                  phi_bin = th1d_fvtxn_phi_weight[fvtx_layer+1]->FindBin(phi); // COME BACK HERE AND HAVE A LOOK
                  fvtx_weight = th1d_fvtxn_phi_weight[fvtx_layer+1]->GetBinContent(phi_bin);
                }
              if ( fvtx_weight != fvtx_weight ) fvtx_weight = 0;
              if ( fvtx_weight < 0 ) fvtx_weight = 0;
              if ( fvtx_weight > 10 ) fvtx_weight = 0;

              // --- north side
              if ( d_FVTX_z[iclus] > 0 )
                {
                  fvtxn_qx2[fvtx_layer+1] += fvtx_weight * TMath::Cos(2*phi);
                  fvtxn_qy2[fvtx_layer+1] += fvtx_weight * TMath::Sin(2*phi);
                  fvtxn_qx3[fvtx_layer+1] += fvtx_weight * TMath::Cos(3*phi);
                  fvtxn_qy3[fvtx_layer+1] += fvtx_weight * TMath::Sin(3*phi);

                  fvtxn_qx2[0] += fvtx_weight * TMath::Cos(2*phi);
                  fvtxn_qy2[0] += fvtx_weight * TMath::Sin(2*phi);
                  fvtxn_qx3[0] += fvtx_weight * TMath::Cos(3*phi);
                  fvtxn_qy3[0] += fvtx_weight * TMath::Sin(3*phi);
                  fvtxn_qx4[0] += fvtx_weight * TMath::Cos(4*phi);
                  fvtxn_qy4[0] += fvtx_weight * TMath::Sin(4*phi);

                  fvtxn_qw[fvtx_layer+1] += fvtx_weight;
                  fvtxn_qw[0] += fvtx_weight;

                  // ------------------------------------------
                  // --- explicitly excluding weights ---------
                  fvtxn_qx2[5] += TMath::Cos(2*phi);
                  fvtxn_qy2[5] += TMath::Sin(2*phi);
                  fvtxn_qx3[5] += TMath::Cos(3*phi);
                  fvtxn_qy3[5] += TMath::Sin(3*phi);
                  fvtxn_qx4[5] += TMath::Cos(4*phi);
                  fvtxn_qy4[5] += TMath::Sin(4*phi);
                  fvtxn_qw[5] += 1;
                  // -----------------------------------------

		  // --- inner eta clusters
		  if ( fvtx_eta > 0.5 && fvtx_eta < 2.0 )
		    {
		      fvtxn_qx2[6] += fvtx_weight * TMath::Cos(2*phi);
		      fvtxn_qy2[6] += fvtx_weight * TMath::Sin(2*phi);
		      fvtxn_qx3[6] += fvtx_weight * TMath::Cos(3*phi);
		      fvtxn_qy3[6] += fvtx_weight * TMath::Sin(3*phi);
		      fvtxn_qx4[6] += fvtx_weight * TMath::Cos(4*phi);
		      fvtxn_qy4[6] += fvtx_weight * TMath::Sin(4*phi);
		      fvtxn_qw[6] += fvtx_weight;
		      ++nclus_north_inner; // good_4_event
		    }
		  // --- outer eta clusters
		  if ( fvtx_eta > 2.0 && fvtx_eta < 3.5 )
		    {
		      fvtxn_qx2[7] += fvtx_weight * TMath::Cos(2*phi);
		      fvtxn_qy2[7] += fvtx_weight * TMath::Sin(2*phi);
		      fvtxn_qx3[7] += fvtx_weight * TMath::Cos(3*phi);
		      fvtxn_qy3[7] += fvtx_weight * TMath::Sin(3*phi);
		      fvtxn_qx4[7] += fvtx_weight * TMath::Cos(4*phi);
		      fvtxn_qy4[7] += fvtx_weight * TMath::Sin(4*phi);
		      fvtxn_qw[7] += fvtx_weight;
		      ++nclus_north_outer; // good_4_event
		    }

                  th1d_fvtxn_clus_phi->Fill(phi);
                  if ( fvtx_layer == 0 ) th1d_fvtxn0_clus_phi->Fill(phi);
                  if ( fvtx_layer == 1 ) th1d_fvtxn1_clus_phi->Fill(phi);
                  if ( fvtx_layer == 2 ) th1d_fvtxn2_clus_phi->Fill(phi);
                  if ( fvtx_layer == 3 ) th1d_fvtxn3_clus_phi->Fill(phi);
                  th1d_fvtxn_clus_wphi->Fill(phi,fvtx_weight);
                  if ( fvtx_layer == 0 ) th1d_fvtxn0_clus_wphi->Fill(phi,fvtx_weight);
                  if ( fvtx_layer == 1 ) th1d_fvtxn1_clus_wphi->Fill(phi,fvtx_weight);
                  if ( fvtx_layer == 2 ) th1d_fvtxn2_clus_wphi->Fill(phi,fvtx_weight);
                  if ( fvtx_layer == 3 ) th1d_fvtxn3_clus_wphi->Fill(phi,fvtx_weight);
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

      bool good_4_event = ( nclus_south_inner > 3 ) && ( nclus_south_outer > 3 ) && ( nclus_north_inner > 3 ) && ( nclus_north_outer > 3 ) ;

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


      // --- fvtx tracks
      float fvtxs_tracks_qx2[3]; // both, inner, outer
      float fvtxs_tracks_qy2[3];
      float fvtxs_tracks_qx3[3];
      float fvtxs_tracks_qy3[3];
      float fvtxs_tracks_qx4[3];
      float fvtxs_tracks_qy4[3];
      float fvtxs_tracks_qw[3];
      float fvtxn_tracks_qx2[3]; // both, inner, outer
      float fvtxn_tracks_qy2[3];
      float fvtxn_tracks_qx3[3];
      float fvtxn_tracks_qy3[3];
      float fvtxn_tracks_qx4[3];
      float fvtxn_tracks_qy4[3];
      float fvtxn_tracks_qw[3];

      for ( int i = 0; i < 3; ++i )
        {
          fvtxs_tracks_qx2[i] = 0.0;
          fvtxs_tracks_qy2[i] = 0.0;
          fvtxs_tracks_qx3[i] = 0.0;
          fvtxs_tracks_qy3[i] = 0.0;
          fvtxs_tracks_qx4[i] = 0.0;
          fvtxs_tracks_qy4[i] = 0.0;
          fvtxs_tracks_qw[i] = 0.0;
          fvtxn_tracks_qx2[i] = 0.0;
          fvtxn_tracks_qy2[i] = 0.0;
          fvtxn_tracks_qx3[i] = 0.0;
          fvtxn_tracks_qy3[i] = 0.0;
          fvtxn_tracks_qx4[i] = 0.0;
          fvtxn_tracks_qy4[i] = 0.0;
          fvtxn_tracks_qw[i] = 0.0;
        } // loop over layers

      int ntrack_south_inner = 0; // good_4_event
      int ntrack_north_inner = 0;
      int ntrack_south_outer = 0;
      int ntrack_north_outer = 0;

      // --- first fvtx track loop
      for ( int i = 0; i < nfvtxt; ++i )
	{
	  // --- rotation now done in trees
	  float phi = fphi[i];
	  float eta = feta[i];

	  bool is_south = ( eta < 0 );
	  bool is_south_inner = ( eta > -2 && eta < 0 );
	  bool is_south_outer = ( eta < -2 );
	  bool is_north = ( eta > 0 );
	  bool is_north_inner = ( eta < 2 && eta > 0 );
	  bool is_north_outer = ( eta > 2 );

	  float fvtx_weight = 1.0; // need to weight tracks at some point...

	  if ( is_south )
	    {
	      fvtxs_tracks_qx2[0] += fvtx_weight * TMath::Cos(2*phi);
	      fvtxs_tracks_qy2[0] += fvtx_weight * TMath::Sin(2*phi);
	      fvtxs_tracks_qx3[0] += fvtx_weight * TMath::Cos(3*phi);
	      fvtxs_tracks_qy3[0] += fvtx_weight * TMath::Sin(3*phi);
	      fvtxs_tracks_qx4[0] += fvtx_weight * TMath::Cos(4*phi);
	      fvtxs_tracks_qy4[0] += fvtx_weight * TMath::Sin(4*phi);
	      fvtxs_tracks_qw[0] += fvtx_weight;
	    }
	  if ( is_south_inner )
	    {
	      fvtxs_tracks_qx2[1] += fvtx_weight * TMath::Cos(2*phi);
	      fvtxs_tracks_qy2[1] += fvtx_weight * TMath::Sin(2*phi);
	      fvtxs_tracks_qx3[1] += fvtx_weight * TMath::Cos(3*phi);
	      fvtxs_tracks_qy3[1] += fvtx_weight * TMath::Sin(3*phi);
	      fvtxs_tracks_qx4[1] += fvtx_weight * TMath::Cos(4*phi);
	      fvtxs_tracks_qy4[1] += fvtx_weight * TMath::Sin(4*phi);
	      fvtxs_tracks_qw[1] += fvtx_weight;
	      ++ntrack_south_inner; // good_4_event
	    }
	  if ( is_south_outer )
	    {
	      fvtxs_tracks_qx2[2] += fvtx_weight * TMath::Cos(2*phi);
	      fvtxs_tracks_qy2[2] += fvtx_weight * TMath::Sin(2*phi);
	      fvtxs_tracks_qx3[2] += fvtx_weight * TMath::Cos(3*phi);
	      fvtxs_tracks_qy3[2] += fvtx_weight * TMath::Sin(3*phi);
	      fvtxs_tracks_qx4[2] += fvtx_weight * TMath::Cos(4*phi);
	      fvtxs_tracks_qy4[2] += fvtx_weight * TMath::Sin(4*phi);
	      fvtxs_tracks_qw[2] += fvtx_weight;
	      ++ntrack_south_outer; // good_4_event
	    }
	  if ( is_north )
	    {
	      fvtxn_tracks_qx2[0] += fvtx_weight * TMath::Cos(2*phi);
	      fvtxn_tracks_qy2[0] += fvtx_weight * TMath::Sin(2*phi);
	      fvtxn_tracks_qx3[0] += fvtx_weight * TMath::Cos(3*phi);
	      fvtxn_tracks_qy3[0] += fvtx_weight * TMath::Sin(3*phi);
	      fvtxn_tracks_qx4[0] += fvtx_weight * TMath::Cos(4*phi);
	      fvtxn_tracks_qy4[0] += fvtx_weight * TMath::Sin(4*phi);
	      fvtxn_tracks_qw[0] += fvtx_weight;
	    }
	  if ( is_north_inner )
	    {
	      fvtxn_tracks_qx2[1] += fvtx_weight * TMath::Cos(2*phi);
	      fvtxn_tracks_qy2[1] += fvtx_weight * TMath::Sin(2*phi);
	      fvtxn_tracks_qx3[1] += fvtx_weight * TMath::Cos(3*phi);
	      fvtxn_tracks_qy3[1] += fvtx_weight * TMath::Sin(3*phi);
	      fvtxn_tracks_qx4[1] += fvtx_weight * TMath::Cos(4*phi);
	      fvtxn_tracks_qy4[1] += fvtx_weight * TMath::Sin(4*phi);
	      fvtxn_tracks_qw[1] += fvtx_weight;
	      ++ntrack_north_inner; // good_4_event
	    }
	  if ( is_north_outer )
	    {
	      fvtxn_tracks_qx2[2] += fvtx_weight * TMath::Cos(2*phi);
	      fvtxn_tracks_qy2[2] += fvtx_weight * TMath::Sin(2*phi);
	      fvtxn_tracks_qx3[2] += fvtx_weight * TMath::Cos(3*phi);
	      fvtxn_tracks_qy3[2] += fvtx_weight * TMath::Sin(3*phi);
	      fvtxn_tracks_qx4[2] += fvtx_weight * TMath::Cos(4*phi);
	      fvtxn_tracks_qy4[2] += fvtx_weight * TMath::Sin(4*phi);
	      fvtxn_tracks_qw[2] += fvtx_weight;
	      ++ntrack_north_outer; // good_4_event
	    }
	}

      good_4_event = good_4_event && ( ntrack_south_inner > 0 ) && ( ntrack_south_outer > 0 ) && ( ntrack_north_inner > 0 ) && ( ntrack_north_outer > 0 ) ;
      int nfvtxt_south = ntrack_south_inner + ntrack_south_outer;
      int nfvtxt_north = ntrack_north_inner + ntrack_north_outer;

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
          sumxy[1][bbcs_nw_index][0] = bbc_nw_qx2;
          sumxy[1][bbcs_nw_index][1] = bbc_nw_qy2;
          sumxy[1][bbcs_nw_index][2] = bbc_nw_qw;
          sumxy[2][bbcs_nw_index][0] = bbc_nw_qx3;
          sumxy[2][bbcs_nw_index][1] = bbc_nw_qy3;
          sumxy[2][bbcs_nw_index][2] = bbc_nw_qw;
          sumxy[3][bbcs_nw_index][0] = bbc_nw_qx4;
          sumxy[3][bbcs_nw_index][1] = bbc_nw_qy4;
          sumxy[3][bbcs_nw_index][2] = bbc_nw_qw;
        }

      if ( fvtx_clusters )
        {
          for ( int i = 0; i < 5; ++i )
            {
              sumxy[1][fvtxs_index+i][0] = fvtxs_qx2[i];
              sumxy[1][fvtxs_index+i][1] = fvtxs_qy2[i];
              sumxy[1][fvtxs_index+i][2] = fvtxs_qw[i];
              sumxy[1][fvtxn_index+i][0] = fvtxn_qx2[i];
              sumxy[1][fvtxn_index+i][1] = fvtxn_qy2[i];
              sumxy[1][fvtxn_index+i][2] = fvtxn_qw[i];
              sumxy[2][fvtxs_index+i][0] = fvtxs_qx3[i];
              sumxy[2][fvtxs_index+i][1] = fvtxs_qy3[i];
              sumxy[2][fvtxs_index+i][2] = fvtxs_qw[i];
              sumxy[2][fvtxn_index+i][0] = fvtxn_qx3[i];
              sumxy[2][fvtxn_index+i][1] = fvtxn_qy3[i];
              sumxy[2][fvtxn_index+i][2] = fvtxn_qw[i];
              sumxy[3][fvtxs_index+i][0] = fvtxs_qx4[i];
              sumxy[3][fvtxs_index+i][1] = fvtxs_qy4[i];
              sumxy[3][fvtxs_index+i][2] = fvtxs_qw[i];
              sumxy[3][fvtxn_index+i][0] = fvtxn_qx4[i];
              sumxy[3][fvtxn_index+i][1] = fvtxn_qy4[i];
              sumxy[3][fvtxn_index+i][2] = fvtxn_qw[i];
              // ---
              sumxy[1][fvtxs_nw_index+i][0] = fvtxs_qx2[i+5];
              sumxy[1][fvtxs_nw_index+i][1] = fvtxs_qy2[i+5];
              sumxy[1][fvtxs_nw_index+i][2] = fvtxs_qw[i+5];
              sumxy[1][fvtxn_nw_index+i][0] = fvtxn_qx2[i+5];
              sumxy[1][fvtxn_nw_index+i][1] = fvtxn_qy2[i+5];
              sumxy[1][fvtxn_nw_index+i][2] = fvtxn_qw[i+5];
              sumxy[2][fvtxs_nw_index+i][0] = fvtxs_qx3[i+5];
              sumxy[2][fvtxs_nw_index+i][1] = fvtxs_qy3[i+5];
              sumxy[2][fvtxs_nw_index+i][2] = fvtxs_qw[i+5];
              sumxy[2][fvtxn_nw_index+i][0] = fvtxn_qx3[i+5];
              sumxy[2][fvtxn_nw_index+i][1] = fvtxn_qy3[i+5];
              sumxy[2][fvtxn_nw_index+i][2] = fvtxn_qw[i+5];
              sumxy[3][fvtxs_nw_index+i][0] = fvtxs_qx4[i+5];
              sumxy[3][fvtxs_nw_index+i][1] = fvtxs_qy4[i+5];
              sumxy[3][fvtxs_nw_index+i][2] = fvtxs_qw[i+5];
              sumxy[3][fvtxn_nw_index+i][0] = fvtxn_qx4[i+5];
              sumxy[3][fvtxn_nw_index+i][1] = fvtxn_qy4[i+5];
              sumxy[3][fvtxn_nw_index+i][2] = fvtxn_qw[i+5];
            } // loop over layers
        } // check on clusters
      if ( fvtx_tracks )
        {
          for ( int i = 0; i < 3; ++i )
            {
              sumxy[1][fvtxs_tracks_index+i][0] = fvtxs_tracks_qx2[i];
              sumxy[1][fvtxs_tracks_index+i][1] = fvtxs_tracks_qy2[i];
              sumxy[1][fvtxs_tracks_index+i][2] = fvtxs_tracks_qw[i];
              sumxy[1][fvtxn_tracks_index+i][0] = fvtxn_tracks_qx2[i];
              sumxy[1][fvtxn_tracks_index+i][1] = fvtxn_tracks_qy2[i];
              sumxy[1][fvtxn_tracks_index+i][2] = fvtxn_tracks_qw[i];
              sumxy[2][fvtxs_tracks_index+i][0] = fvtxs_tracks_qx3[i];
              sumxy[2][fvtxs_tracks_index+i][1] = fvtxs_tracks_qy3[i];
              sumxy[2][fvtxs_tracks_index+i][2] = fvtxs_tracks_qw[i];
              sumxy[2][fvtxn_tracks_index+i][0] = fvtxn_tracks_qx3[i];
              sumxy[2][fvtxn_tracks_index+i][1] = fvtxn_tracks_qy3[i];
              sumxy[2][fvtxn_tracks_index+i][2] = fvtxn_tracks_qw[i];
              sumxy[3][fvtxs_tracks_index+i][0] = fvtxs_tracks_qx4[i];
              sumxy[3][fvtxs_tracks_index+i][1] = fvtxs_tracks_qy4[i];
              sumxy[3][fvtxs_tracks_index+i][2] = fvtxs_tracks_qw[i];
              sumxy[3][fvtxn_tracks_index+i][0] = fvtxn_tracks_qx4[i];
              sumxy[3][fvtxn_tracks_index+i][1] = fvtxn_tracks_qy4[i];
              sumxy[3][fvtxn_tracks_index+i][2] = fvtxn_tracks_qw[i];
            } // loop over layers
        } // check on tracks


      if ( DIAG )
        {
          cout<<"bbc from node tree: "<<d_Qx[5]<<" "<<d_Qy[5]<<" "<<d_Qw[5]<<endl;
          cout<<"bbc from me: "<<bbc_nw_qx2<<" "<<bbc_nw_qy2<<" "<<bbc_nw_qw<<endl;

          cout<<"fvtx raw: "<<endl;
          cout<<"from node tree: "<<d_Qx[4]<<" "<<d_Qy[4]<<" "<<d_Qw[4]<<endl;
          cout<<"from clusters: " <<fvtxs_qx2[0]<<" "<<fvtxs_qy2[0]<<" "<<fvtxs_qw[0]<<endl;
        }

      for ( int ih = 1; ih < NHAR; ih++ )
	{
	  for (int id = 0; id < NDET; id++ )
	    {
	      if(sumxy[ih][id][2]>0)
		{
		  //float psi = atan2(sumxy[ih][id][1],sumxy[ih][id][0])/2.0;
		  float psi = atan2(sumxy[ih][id][1],sumxy[ih][id][0])/float(ih+1);
		  if ( DIAG ) cout<<"RAW: for id: "<<id<<" psi: "<<psi<<endl;
		  psi_bf[icent][ih][id]->Fill(izvtx,psi);
		} // check on weight
	    } // detectors
	} // harmonics

      //------------------------------------------------------------//
      //                Flattening iteration                        //
      //------------------------------------------------------------//
      //int icent = 0;
      for ( int ih = 1; ih < NHAR; ih++ )
	{
	  for( int id = 0; id < NDET; id++ )
	    {
	      if ( sumxy[ih][id][2] > 0.0 )
		{
		  sumxy[ih][id][3]=atan2(sumxy[ih][id][1],sumxy[ih][id][0])/(ih+1.0);
		  if (rp_recal_pass>0) dis[icent][ih][id]->Fill(izvtx,sumxy[ih][id][3]*(ih+1.0));
		}
	      if ( sumxy[ih][id][2] > 0.0 ) // check on weight (x,y,w,psi)
		{
		  for ( int ib = 0; ib < 2; ib++ )
		    {
		      sumxy[ih][id][ib]/=sumxy[ih][id][2]; // normalize to the weight
		      //if(ih==1 && id==0 && ib==0 && sumxy[ih][id][ib]>1) cout<<sumxy[ih][id][ib]<<endl;
		      if ( rp_recal_pass > 0 )
			{
			  ave[icent][izvtx][ih][id]->Fill(ib+0.0,sumxy[ih][id][ib]);
			  if(id==0 && DIAG) cout<<"filled ave: "<<ih<<" "<<id<<" "<<ib<<" with: "<<sumxy[ih][id][ib]<<endl;
			  if(ib==0) qx[icent][ih][id]->Fill(izvtx,sumxy[ih][id][0]);
			  if(ib==1) qy[icent][ih][id]->Fill(izvtx,sumxy[ih][id][1]);
			} // pass > 0
		      float sxy=sumxy[ih][id][ib];
		      float mxy=mean[icent][izvtx][ih][id][ib]; // for recentering qx and qy (???)
		      float wxy=widt[icent][izvtx][ih][id][ib]; // for recentering qx and qy (???)

		      //if(ic==0 && izvtx==0 && ih==1 && id==0) cout<<ib<<" "<<sxy<<" "<<mxy<<" "<<wxy<<endl;
		      sumxy[ih][id][ib]=(sxy-mxy)/wxy; // recentered by mean and renormalized to width
		      if ( rp_recal_pass > 0 )
			{
			  ave[icent][izvtx][ih][id]->Fill(ib+2.0,sumxy[ih][id][ib]);  // ib+2 to avoid overlap
			  if(id==0 && DIAG) cout<<"filled ave2: "<<ih<<" "<<id<<" "<<ib<<" with: "<<sumxy[ih][id][ib]<<endl;
			  if(ib==0) qx[icent][ih][id]->Fill(izvtx+NZPS,sumxy[ih][id][0]);
			  if(ib==1) qy[icent][ih][id]->Fill(izvtx+NZPS,sumxy[ih][id][1]);
			} // pass > 0
		    } // if weight > 0

		  sumxy[ih][id][3]=atan2(sumxy[ih][id][1],sumxy[ih][id][0])/(ih+1.0);
		  if ( rp_recal_pass > 0 )
		    {
		      // fill histogram with psi calculated with recenter q vectors(??)
		      dis[icent][ih][id]->Fill(izvtx+NZPS,sumxy[ih][id][3]*(ih+1.0));
		      // my own simpler version of the above histogram
		      psi_mf[icent][ih][id]->Fill(izvtx,sumxy[ih][id][3]);
		    }

		  float psi = sumxy[ih][id][3]*(ih+1.0);
		  if ( ih == 1 && id == 0 && DIAG )  cout<<"psi-1 bbc: "<<psi<<endl;
		  float dp = 0.0;
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
		  psi += dp; // shift psi by...
		  psi = atan2(sin(psi),cos(psi)); // trick to readjust the range
		  if ( ih == 1 && id == 0 && DIAG )  cout<<"psi-2 bbc: "<<psi<<endl;
		  for (int io=0; io<NORD; io++)
		    {
		      float cc=cos((io+1.0)*psi);
		      float ss=sin((io+1.0)*psi);
		      // --- fourier components of modified psi
		      if (rp_recal_pass>0) flt[icent][izvtx][ih][id]->Fill(io+NORD*2.0,cc);
		      if (rp_recal_pass>0) flt[icent][izvtx][ih][id]->Fill(io+NORD*3.0,ss);
		    }
		  sumxy[ih][id][3]=psi/(ih+1.0);
		  if ( rp_recal_pass > 0 ) dis[icent][ih][id]->Fill(izvtx+NZPS*2.0,sumxy[ih][id][3]*(ih+1.0));
		} // end if weight > 0
	      else
		{
		  sumxy[ih][id][3]=-9999.9;
		} // otherwise set psi to some crazy number
	    } // detectors
	} // harmonics



      if ( DIAG ) cout<<"bbc_rp2: "<<sumxy[1][0][3]<<endl;

      for ( int ih = 1; ih < NHAR; ih++ )
	{
	  for ( int id = 0; id < NDET; id++ )
	    {
	      if ( sumxy[ih][id][2] > 0 )
		{
		  psi_af[icent][ih][id]->Fill(izvtx,sumxy[ih][id][3]);
		  if ( DIAG ) cout<<"CORR: for id: "<<id<<" psi: "<<sumxy[ih][id][3]<<endl;
		} // check on weight
	    } // detectors
	} // harmonics
      // ---
      // --- now going to calculate v2
      // ---
      if ( rp_recal_pass < 3 ) continue; // don't calculate v2 except for final pass

      // --- let's see if we can figure out the mean qx and qy from the flattening procedure above...
      // --- this one? mean[icent][izvtx][ih][id][ib] typical agreement better than 1%
      // --- or maybe this one? sumxy[ih][id][ib] actually i think these are recentered Q-vectors, but if that's the case i could use them directly...
      // --- adding fourth harmonic to the flattening procedure for the first time... 20160803-1702 EDT
      float os_bbc_qw = bbc_qw; // the same
      float os_bbc_qx2 = bbc_qx2 - os_bbc_qw*mean[icent][izvtx][1][bbcs_index][0];
      float os_bbc_qy2 = bbc_qy2 - os_bbc_qw*mean[icent][izvtx][1][bbcs_index][1];
      float os_bbc_qx3 = bbc_qx3 - os_bbc_qw*mean[icent][izvtx][2][bbcs_index][0];
      float os_bbc_qy3 = bbc_qy3 - os_bbc_qw*mean[icent][izvtx][2][bbcs_index][1];
      float os_bbc_qx4 = bbc_qx4 - os_bbc_qw*mean[icent][izvtx][3][bbcs_index][0];
      float os_bbc_qy4 = bbc_qy4 - os_bbc_qw*mean[icent][izvtx][3][bbcs_index][1];
      float os_bbc_qq2 = calc2_event(os_bbc_qx2,os_bbc_qy2,os_bbc_qw);
      float os_bbc_qq3 = calc2_event(os_bbc_qx3,os_bbc_qy3,os_bbc_qw);
      os_bbcs_cos22->Fill(0.0,os_bbc_qx2);
      os_bbcs_sin22->Fill(0.0,os_bbc_qy2);
      os_bbcs_cos32->Fill(0.0,os_bbc_qx3);
      os_bbcs_sin32->Fill(0.0,os_bbc_qy3);
      os_bbcs_c22->Fill(0.0,os_bbc_qq2);
      os_bbcs_c32->Fill(0.0,os_bbc_qq3);
      npc1_os_bbcs_c22->Fill(npc1,os_bbc_qq2);
      npc1_os_bbcs_c32->Fill(npc1,os_bbc_qq3);
      nfvtxt_os_bbcs_c22->Fill(nfvtxt,os_bbc_qq2);
      nfvtxc_os_bbcs_c22->Fill(nfvtxc,os_bbc_qq2);
      nfvtxt_os_bbcs_c32->Fill(nfvtxt,os_bbc_qq3);
      nfvtxc_os_bbcs_c32->Fill(nfvtxc,os_bbc_qq3);
      float os_bbc_psi2 = atan2(os_bbc_qy2,os_bbc_qx2)/2.0;
      float os_bbc_psi3 = atan2(os_bbc_qy3,os_bbc_qx3)/3.0;
      os_bbcs_1dPsi2->Fill(os_bbc_psi2);
      os_bbcs_1dPsi3->Fill(os_bbc_psi3);

      float os_fvtxs_qw = fvtxs_qw[0];
      float os_fvtxs_qx2 = fvtxs_qx2[0] - os_fvtxs_qw*mean[icent][izvtx][1][fvtxs_index][0];
      float os_fvtxs_qy2 = fvtxs_qy2[0] - os_fvtxs_qw*mean[icent][izvtx][1][fvtxs_index][1];
      float os_fvtxs_qx3 = fvtxs_qx3[0] - os_fvtxs_qw*mean[icent][izvtx][2][fvtxs_index][0];
      float os_fvtxs_qy3 = fvtxs_qy3[0] - os_fvtxs_qw*mean[icent][izvtx][2][fvtxs_index][1];
      float os_fvtxs_qx4 = fvtxs_qx4[0] - os_fvtxs_qw*mean[icent][izvtx][3][fvtxs_index][0];
      float os_fvtxs_qy4 = fvtxs_qy4[0] - os_fvtxs_qw*mean[icent][izvtx][3][fvtxs_index][1];
      float os_fvtxs_qq2 = calc2_event(os_fvtxs_qx2,os_fvtxs_qy2,os_fvtxs_qw);
      float os_fvtxs_qq3 = calc2_event(os_fvtxs_qx3,os_fvtxs_qy3,os_fvtxs_qw);
      os_fvtxs_cos22->Fill(0.0,os_fvtxs_qx2);
      os_fvtxs_sin22->Fill(0.0,os_fvtxs_qy2);
      os_fvtxs_cos32->Fill(0.0,os_fvtxs_qx3);
      os_fvtxs_sin32->Fill(0.0,os_fvtxs_qy3);
      os_fvtxs_c22->Fill(0.0,os_fvtxs_qq2);
      os_fvtxs_c32->Fill(0.0,os_fvtxs_qq3);
      npc1_os_fvtxs_c22->Fill(npc1,os_fvtxs_qq2);
      npc1_os_fvtxs_c32->Fill(npc1,os_fvtxs_qq3);
      nfvtxt_os_fvtxs_c22->Fill(nfvtxt_south,os_fvtxs_qq2);
      nfvtxc_os_fvtxs_c22->Fill(nfvtxc_south,os_fvtxs_qq2);
      nfvtxt_os_fvtxs_c32->Fill(nfvtxt_south,os_fvtxs_qq3);
      nfvtxc_os_fvtxs_c32->Fill(nfvtxc_south,os_fvtxs_qq3);
      float os_fvtxs_psi2 = atan2(os_fvtxs_qy2,os_fvtxs_qx2)/2.0;
      float os_fvtxs_psi3 = atan2(os_fvtxs_qy3,os_fvtxs_qx3)/3.0;
      os_fvtxs_1dPsi2->Fill(os_fvtxs_psi2);
      os_fvtxs_1dPsi3->Fill(os_fvtxs_psi3);

      float os_fvtxn_qw = fvtxn_qw[0];
      float os_fvtxn_qx2 = fvtxn_qx2[0] - os_fvtxn_qw*mean[icent][izvtx][1][fvtxn_index][0];
      float os_fvtxn_qy2 = fvtxn_qy2[0] - os_fvtxn_qw*mean[icent][izvtx][1][fvtxn_index][1];
      float os_fvtxn_qx3 = fvtxn_qx3[0] - os_fvtxn_qw*mean[icent][izvtx][2][fvtxn_index][0];
      float os_fvtxn_qy3 = fvtxn_qy3[0] - os_fvtxn_qw*mean[icent][izvtx][2][fvtxn_index][1];
      float os_fvtxn_qx4 = fvtxn_qx4[0] - os_fvtxn_qw*mean[icent][izvtx][3][fvtxn_index][0];
      float os_fvtxn_qy4 = fvtxn_qy4[0] - os_fvtxn_qw*mean[icent][izvtx][3][fvtxn_index][1];
      float os_fvtxn_qq2 = calc2_event(os_fvtxn_qx2,os_fvtxn_qy2,os_fvtxn_qw);
      float os_fvtxn_qq3 = calc2_event(os_fvtxn_qx3,os_fvtxn_qy3,os_fvtxn_qw);
      os_fvtxn_cos22->Fill(0.0,os_fvtxn_qx2);
      os_fvtxn_sin22->Fill(0.0,os_fvtxn_qy2);
      os_fvtxn_cos32->Fill(0.0,os_fvtxn_qx3);
      os_fvtxn_sin32->Fill(0.0,os_fvtxn_qy3);
      os_fvtxn_c22->Fill(0.0,os_fvtxn_qq2);
      os_fvtxn_c32->Fill(0.0,os_fvtxn_qq3);
      npc1_os_fvtxn_c22->Fill(npc1,os_fvtxn_qq2);
      npc1_os_fvtxn_c32->Fill(npc1,os_fvtxn_qq3);
      nfvtxt_os_fvtxn_c22->Fill(nfvtxt_north,os_fvtxn_qq2);
      nfvtxc_os_fvtxn_c22->Fill(nfvtxc_north,os_fvtxn_qq2);
      nfvtxt_os_fvtxn_c32->Fill(nfvtxt_north,os_fvtxn_qq3);
      nfvtxc_os_fvtxn_c32->Fill(nfvtxc_north,os_fvtxn_qq3);
      float os_fvtxn_psi2 = atan2(os_fvtxn_qy2,os_fvtxn_qx2)/2.0;
      float os_fvtxn_psi3 = atan2(os_fvtxn_qy3,os_fvtxn_qx3)/3.0;
      os_fvtxn_1dPsi2->Fill(os_fvtxn_psi2);
      os_fvtxn_1dPsi3->Fill(os_fvtxn_psi3);

      // --- eta dependent clusters

      float os_fvtxs_ce0_qw = fvtxs_qw[6];
      float os_fvtxs_ce0_qx2 = fvtxs_qx2[6] - os_fvtxs_ce0_qw*mean[icent][izvtx][1][fvtxs0_nw_index][0];
      float os_fvtxs_ce0_qy2 = fvtxs_qy2[6] - os_fvtxs_ce0_qw*mean[icent][izvtx][1][fvtxs0_nw_index][1];
      float os_fvtxs_ce0_qx3 = fvtxs_qx3[6] - os_fvtxs_ce0_qw*mean[icent][izvtx][2][fvtxs0_nw_index][0];
      float os_fvtxs_ce0_qy3 = fvtxs_qy3[6] - os_fvtxs_ce0_qw*mean[icent][izvtx][2][fvtxs0_nw_index][1];
      // float os_fvtxs_ce0_qx4 = fvtxs_qx4[6] - os_fvtxs_ce0_qw*mean[icent][izvtx][3][fvtxs0_nw_index][0];
      // float os_fvtxs_ce0_qy4 = fvtxs_qy4[6] - os_fvtxs_ce0_qw*mean[icent][izvtx][3][fvtxs0_nw_index][1];
      float os_fvtxs_ce1_qw = fvtxs_qw[7];
      float os_fvtxs_ce1_qx2 = fvtxs_qx2[7] - os_fvtxs_ce1_qw*mean[icent][izvtx][1][fvtxs1_nw_index][0];
      float os_fvtxs_ce1_qy2 = fvtxs_qy2[7] - os_fvtxs_ce1_qw*mean[icent][izvtx][1][fvtxs1_nw_index][1];
      float os_fvtxs_ce1_qx3 = fvtxs_qx3[7] - os_fvtxs_ce1_qw*mean[icent][izvtx][2][fvtxs1_nw_index][0];
      float os_fvtxs_ce1_qy3 = fvtxs_qy3[7] - os_fvtxs_ce1_qw*mean[icent][izvtx][2][fvtxs1_nw_index][1];
      // float os_fvtxs_ce1_qx4 = fvtxs_qx4[7] - os_fvtxs_ce1_qw*mean[icent][izvtx][3][fvtxs1_nw_index][0];
      // float os_fvtxs_ce1_qy4 = fvtxs_qy4[7] - os_fvtxs_ce1_qw*mean[icent][izvtx][3][fvtxs1_nw_index][1];
      float os_fvtxn_ce0_qw = fvtxn_qw[6];
      float os_fvtxn_ce0_qx2 = fvtxn_qx2[6] - os_fvtxn_ce0_qw*mean[icent][izvtx][1][fvtxn0_nw_index][0];
      float os_fvtxn_ce0_qy2 = fvtxn_qy2[6] - os_fvtxn_ce0_qw*mean[icent][izvtx][1][fvtxn0_nw_index][1];
      float os_fvtxn_ce0_qx3 = fvtxn_qx3[6] - os_fvtxn_ce0_qw*mean[icent][izvtx][2][fvtxn0_nw_index][0];
      float os_fvtxn_ce0_qy3 = fvtxn_qy3[6] - os_fvtxn_ce0_qw*mean[icent][izvtx][2][fvtxn0_nw_index][1];
      // float os_fvtxn_ce0_qx4 = fvtxn_qx4[6] - os_fvtxn_ce0_qw*mean[icent][izvtx][3][fvtxn0_nw_index][0];
      // float os_fvtxn_ce0_qy4 = fvtxn_qy4[6] - os_fvtxn_ce0_qw*mean[icent][izvtx][3][fvtxn0_nw_index][1];
      float os_fvtxn_ce1_qw = fvtxn_qw[7];
      float os_fvtxn_ce1_qx2 = fvtxn_qx2[7] - os_fvtxn_ce1_qw*mean[icent][izvtx][1][fvtxn1_nw_index][0];
      float os_fvtxn_ce1_qy2 = fvtxn_qy2[7] - os_fvtxn_ce1_qw*mean[icent][izvtx][1][fvtxn1_nw_index][1];
      float os_fvtxn_ce1_qx3 = fvtxn_qx3[7] - os_fvtxn_ce1_qw*mean[icent][izvtx][2][fvtxn1_nw_index][0];
      float os_fvtxn_ce1_qy3 = fvtxn_qy3[7] - os_fvtxn_ce1_qw*mean[icent][izvtx][2][fvtxn1_nw_index][1];
      // float os_fvtxn_ce1_qx4 = fvtxn_qx4[7] - os_fvtxn_ce1_qw*mean[icent][izvtx][3][fvtxn1_nw_index][0];
      // float os_fvtxn_ce1_qy4 = fvtxn_qy4[7] - os_fvtxn_ce1_qw*mean[icent][izvtx][3][fvtxn1_nw_index][1];
      // ---
      float os_fvtxs_ce01_qq2 = ( os_fvtxs_ce0_qx2*os_fvtxs_ce1_qx2 + os_fvtxs_ce0_qy2*os_fvtxs_ce1_qy2 ) / ( os_fvtxs_ce0_qw*os_fvtxs_ce1_qw );
      float os_fvtxs_ce01_qq3 = ( os_fvtxs_ce0_qx3*os_fvtxs_ce1_qx3 + os_fvtxs_ce0_qy3*os_fvtxs_ce1_qy3 ) / ( os_fvtxs_ce0_qw*os_fvtxs_ce1_qw );
      float os_fvtxn_ce01_qq2 = ( os_fvtxn_ce0_qx2*os_fvtxn_ce1_qx2 + os_fvtxn_ce0_qy2*os_fvtxn_ce1_qy2 ) / ( os_fvtxn_ce0_qw*os_fvtxn_ce1_qw );
      float os_fvtxn_ce01_qq3 = ( os_fvtxn_ce0_qx3*os_fvtxn_ce1_qx3 + os_fvtxn_ce0_qy3*os_fvtxn_ce1_qy3 ) / ( os_fvtxn_ce0_qw*os_fvtxn_ce1_qw );
      // --- come back here to add some histograms and do some 4 particle stuff
      os_fvtxs_ce01_c22->Fill(0.0,os_fvtxs_ce01_qq2);
      os_fvtxs_ce01_c32->Fill(0.0,os_fvtxs_ce01_qq3);
      os_fvtxn_ce01_c22->Fill(0.0,os_fvtxn_ce01_qq2);
      os_fvtxn_ce01_c32->Fill(0.0,os_fvtxn_ce01_qq3);
      nfvtxt_os_fvtxs_ce01_c22->Fill(nfvtxt_south,os_fvtxs_ce01_qq2);
      nfvtxc_os_fvtxs_ce01_c22->Fill(nfvtxc_south,os_fvtxs_ce01_qq2);
      nfvtxt_os_fvtxs_ce01_c32->Fill(nfvtxt_south,os_fvtxs_ce01_qq3);
      nfvtxc_os_fvtxs_ce01_c32->Fill(nfvtxc_south,os_fvtxs_ce01_qq3);
      nfvtxt_os_fvtxn_ce01_c22->Fill(nfvtxt_north,os_fvtxn_ce01_qq2);
      nfvtxc_os_fvtxn_ce01_c22->Fill(nfvtxc_north,os_fvtxn_ce01_qq2);
      nfvtxt_os_fvtxn_ce01_c32->Fill(nfvtxt_north,os_fvtxn_ce01_qq3);
      nfvtxc_os_fvtxn_ce01_c32->Fill(nfvtxc_north,os_fvtxn_ce01_qq3);

      // --- now fvtx tracks

      float os_fvtxs_tracks_qw = fvtxs_tracks_qw[0];
      float os_fvtxs_tracks_qx2 = fvtxs_tracks_qx2[0] - os_fvtxs_tracks_qw*mean[icent][izvtx][1][fvtxs_tracks_index][0];
      float os_fvtxs_tracks_qy2 = fvtxs_tracks_qy2[0] - os_fvtxs_tracks_qw*mean[icent][izvtx][1][fvtxs_tracks_index][1];
      float os_fvtxs_tracks_qx3 = fvtxs_tracks_qx3[0] - os_fvtxs_tracks_qw*mean[icent][izvtx][2][fvtxs_tracks_index][0];
      float os_fvtxs_tracks_qy3 = fvtxs_tracks_qy3[0] - os_fvtxs_tracks_qw*mean[icent][izvtx][2][fvtxs_tracks_index][1];
      float os_fvtxs_tracks_qx4 = fvtxs_tracks_qx4[0] - os_fvtxs_tracks_qw*mean[icent][izvtx][3][fvtxs_tracks_index][0];
      float os_fvtxs_tracks_qy4 = fvtxs_tracks_qy4[0] - os_fvtxs_tracks_qw*mean[icent][izvtx][3][fvtxs_tracks_index][1];
      float os_fvtxs_tracks_qq2 = calc2_event(os_fvtxs_tracks_qx2,os_fvtxs_tracks_qy2,os_fvtxs_tracks_qw);
      float os_fvtxs_tracks_qq3 = calc2_event(os_fvtxs_tracks_qx3,os_fvtxs_tracks_qy3,os_fvtxs_tracks_qw);
      os_fvtxs_tracks_cos22->Fill(0.0,os_fvtxs_tracks_qx2);
      os_fvtxs_tracks_sin22->Fill(0.0,os_fvtxs_tracks_qy2);
      os_fvtxs_tracks_cos32->Fill(0.0,os_fvtxs_tracks_qx3);
      os_fvtxs_tracks_sin32->Fill(0.0,os_fvtxs_tracks_qy3);
      os_fvtxs_tracks_c22->Fill(0.0,os_fvtxs_tracks_qq2);
      os_fvtxs_tracks_c32->Fill(0.0,os_fvtxs_tracks_qq3);
      npc1_os_fvtxs_tracks_c22->Fill(npc1,os_fvtxs_tracks_qq2);
      npc1_os_fvtxs_tracks_c32->Fill(npc1,os_fvtxs_tracks_qq3);
      nfvtxt_os_fvtxs_tracks_c22->Fill(nfvtxt_south,os_fvtxs_tracks_qq2);
      nfvtxc_os_fvtxs_tracks_c22->Fill(nfvtxc_south,os_fvtxs_tracks_qq2);
      nfvtxt_os_fvtxs_tracks_c32->Fill(nfvtxt_south,os_fvtxs_tracks_qq3);
      nfvtxc_os_fvtxs_tracks_c32->Fill(nfvtxc_south,os_fvtxs_tracks_qq3);
      float os_fvtxs_tracks_psi2 = atan2(os_fvtxs_tracks_qy2,os_fvtxs_tracks_qx2)/2.0;
      float os_fvtxs_tracks_psi3 = atan2(os_fvtxs_tracks_qy3,os_fvtxs_tracks_qx3)/3.0;
      os_fvtxs_tracks_1dPsi2->Fill(os_fvtxs_tracks_psi2);
      os_fvtxs_tracks_1dPsi3->Fill(os_fvtxs_tracks_psi3);

      float os_fvtxn_tracks_qw = fvtxn_tracks_qw[0];
      float os_fvtxn_tracks_qx2 = fvtxn_tracks_qx2[0] - os_fvtxn_tracks_qw*mean[icent][izvtx][1][fvtxn_tracks_index][0];
      float os_fvtxn_tracks_qy2 = fvtxn_tracks_qy2[0] - os_fvtxn_tracks_qw*mean[icent][izvtx][1][fvtxn_tracks_index][1];
      float os_fvtxn_tracks_qx3 = fvtxn_tracks_qx3[0] - os_fvtxn_tracks_qw*mean[icent][izvtx][2][fvtxn_tracks_index][0];
      float os_fvtxn_tracks_qy3 = fvtxn_tracks_qy3[0] - os_fvtxn_tracks_qw*mean[icent][izvtx][2][fvtxn_tracks_index][1];
      float os_fvtxn_tracks_qx4 = fvtxn_tracks_qx4[0] - os_fvtxn_tracks_qw*mean[icent][izvtx][3][fvtxn_tracks_index][0];
      float os_fvtxn_tracks_qy4 = fvtxn_tracks_qy4[0] - os_fvtxn_tracks_qw*mean[icent][izvtx][3][fvtxn_tracks_index][1];
      float os_fvtxn_tracks_qq2 = calc2_event(os_fvtxn_tracks_qx2,os_fvtxn_tracks_qy2,os_fvtxn_tracks_qw);
      float os_fvtxn_tracks_qq3 = calc2_event(os_fvtxn_tracks_qx3,os_fvtxn_tracks_qy3,os_fvtxn_tracks_qw);
      os_fvtxn_tracks_cos22->Fill(0.0,os_fvtxn_tracks_qx2);
      os_fvtxn_tracks_sin22->Fill(0.0,os_fvtxn_tracks_qy2);
      os_fvtxn_tracks_cos32->Fill(0.0,os_fvtxn_tracks_qx3);
      os_fvtxn_tracks_sin32->Fill(0.0,os_fvtxn_tracks_qy3);
      os_fvtxn_tracks_c22->Fill(0.0,os_fvtxn_tracks_qq2);
      os_fvtxn_tracks_c32->Fill(0.0,os_fvtxn_tracks_qq3);
      npc1_os_fvtxn_tracks_c22->Fill(npc1,os_fvtxn_tracks_qq2);
      npc1_os_fvtxn_tracks_c32->Fill(npc1,os_fvtxn_tracks_qq3);
      nfvtxt_os_fvtxn_tracks_c22->Fill(nfvtxt_north,os_fvtxn_tracks_qq2);
      nfvtxc_os_fvtxn_tracks_c22->Fill(nfvtxc_north,os_fvtxn_tracks_qq2);
      nfvtxt_os_fvtxn_tracks_c32->Fill(nfvtxt_north,os_fvtxn_tracks_qq3);
      nfvtxc_os_fvtxn_tracks_c32->Fill(nfvtxc_north,os_fvtxn_tracks_qq3);
      float os_fvtxn_tracks_psi2 = atan2(os_fvtxn_tracks_qy2,os_fvtxn_tracks_qx2)/2.0;
      float os_fvtxn_tracks_psi3 = atan2(os_fvtxn_tracks_qy3,os_fvtxn_tracks_qx3)/3.0;
      os_fvtxn_tracks_1dPsi2->Fill(os_fvtxn_tracks_psi2);
      os_fvtxn_tracks_1dPsi3->Fill(os_fvtxn_tracks_psi3);

      // --- combined fvtx clusters

      float os_fvtxc_qx2 = os_fvtxs_qx2 + os_fvtxn_qx2;
      float os_fvtxc_qy2 = os_fvtxs_qy2 + os_fvtxn_qy2;
      float os_fvtxc_qx3 = os_fvtxs_qx3 + os_fvtxn_qx3;
      float os_fvtxc_qy3 = os_fvtxs_qy3 + os_fvtxn_qy3;
      float os_fvtxc_qx4 = os_fvtxs_qx4 + os_fvtxn_qx4;
      float os_fvtxc_qy4 = os_fvtxs_qy4 + os_fvtxn_qy4;
      float os_fvtxc_qw = os_fvtxs_qw + os_fvtxn_qw;
      float os_fvtxc_qq2 = calc2_event(os_fvtxc_qx2,os_fvtxc_qy2,os_fvtxc_qw);
      float os_fvtxc_qq3 = calc2_event(os_fvtxc_qx3,os_fvtxc_qy3,os_fvtxc_qw);
      os_fvtxc_c22->Fill(0.0,os_fvtxc_qq2);
      os_fvtxc_c32->Fill(0.0,os_fvtxc_qq3);
      nfvtxt_os_fvtxc_c22->Fill(nfvtxt,os_fvtxc_qq2);
      nfvtxc_os_fvtxc_c22->Fill(nfvtxc,os_fvtxc_qq2);
      nfvtxt_os_fvtxc_c32->Fill(nfvtxt,os_fvtxc_qq3);
      nfvtxc_os_fvtxc_c32->Fill(nfvtxc,os_fvtxc_qq3);

      // --- combined fvtx tracks

      float os_fvtxc_tracks_qx2 = os_fvtxs_tracks_qx2 + os_fvtxn_tracks_qx2;
      float os_fvtxc_tracks_qy2 = os_fvtxs_tracks_qy2 + os_fvtxn_tracks_qy2;
      float os_fvtxc_tracks_qx3 = os_fvtxs_tracks_qx3 + os_fvtxn_tracks_qx3;
      float os_fvtxc_tracks_qy3 = os_fvtxs_tracks_qy3 + os_fvtxn_tracks_qy3;
      float os_fvtxc_tracks_qx4 = os_fvtxs_tracks_qx4 + os_fvtxn_tracks_qx4;
      float os_fvtxc_tracks_qy4 = os_fvtxs_tracks_qy4 + os_fvtxn_tracks_qy4;
      float os_fvtxc_tracks_qw = os_fvtxs_tracks_qw + os_fvtxn_tracks_qw;
      float os_fvtxc_tracks_qq2 = calc2_event(os_fvtxc_tracks_qx2,os_fvtxc_tracks_qy2,os_fvtxc_tracks_qw);
      float os_fvtxc_tracks_qq3 = calc2_event(os_fvtxc_tracks_qx3,os_fvtxc_tracks_qy3,os_fvtxc_tracks_qw);
      os_fvtxc_tracks_c22->Fill(0.0,os_fvtxc_tracks_qq2);
      os_fvtxc_tracks_c32->Fill(0.0,os_fvtxc_tracks_qq3);
      nfvtxt_os_fvtxc_tracks_c22->Fill(nfvtxt,os_fvtxc_tracks_qq2);
      nfvtxc_os_fvtxc_tracks_c22->Fill(nfvtxc,os_fvtxc_tracks_qq2);
      nfvtxt_os_fvtxc_tracks_c32->Fill(nfvtxt,os_fvtxc_tracks_qq3);
      nfvtxc_os_fvtxc_tracks_c32->Fill(nfvtxc,os_fvtxc_tracks_qq3);

      // --- have a look at some different correlations

      // --- bbcs-fvtxs
      float os_bbcfvtxs_qq2 = ( (os_bbc_qx2*os_fvtxs_qx2) + (os_bbc_qy2*os_fvtxs_qy2) ) / ( os_bbc_qw*os_fvtxs_qw );
      float os_bbcfvtxs_qq3 = ( (os_bbc_qx3*os_fvtxs_qx3) + (os_bbc_qy3*os_fvtxs_qy3) ) / ( os_bbc_qw*os_fvtxs_qw );
      os_bbcsfvtxs_c22->Fill(0.0,os_bbcfvtxs_qq2);
      os_bbcsfvtxs_c32->Fill(0.0,os_bbcfvtxs_qq3);
      npc1_os_bbcsfvtxs_c22->Fill(npc1,os_bbcfvtxs_qq2);
      npc1_os_bbcsfvtxs_c32->Fill(npc1,os_bbcfvtxs_qq3);
      nfvtxt_os_bbcsfvtxs_c22->Fill(nfvtxt,os_bbcfvtxs_qq2);
      nfvtxc_os_bbcsfvtxs_c22->Fill(nfvtxc,os_bbcfvtxs_qq2);
      nfvtxt_os_bbcsfvtxs_c32->Fill(nfvtxt,os_bbcfvtxs_qq3);
      nfvtxc_os_bbcsfvtxs_c32->Fill(nfvtxc,os_bbcfvtxs_qq3);
      // --- bbcs-fvtxn
      float os_bbcfvtxn_qq2 = ( (os_bbc_qx2*os_fvtxn_qx2) + (os_bbc_qy2*os_fvtxn_qy2) ) / ( os_bbc_qw*os_fvtxn_qw );
      float os_bbcfvtxn_qq3 = ( (os_bbc_qx3*os_fvtxn_qx3) + (os_bbc_qy3*os_fvtxn_qy3) ) / ( os_bbc_qw*os_fvtxn_qw );
      os_bbcsfvtxn_c22->Fill(0.0,os_bbcfvtxn_qq2);
      os_bbcsfvtxn_c32->Fill(0.0,os_bbcfvtxn_qq3);
      npc1_os_bbcsfvtxn_c22->Fill(npc1,os_bbcfvtxn_qq2);
      npc1_os_bbcsfvtxn_c32->Fill(npc1,os_bbcfvtxn_qq3);
      nfvtxt_os_bbcsfvtxn_c22->Fill(nfvtxt,os_bbcfvtxn_qq2);
      nfvtxc_os_bbcsfvtxn_c22->Fill(nfvtxc,os_bbcfvtxn_qq2);
      nfvtxt_os_bbcsfvtxn_c32->Fill(nfvtxt,os_bbcfvtxn_qq3);
      nfvtxc_os_bbcsfvtxn_c32->Fill(nfvtxc,os_bbcfvtxn_qq3);
      // --- fvtxs-fvtxn
      float os_fvtxsfvtxn_qq2 = ( (os_fvtxs_qx2*os_fvtxn_qx2) + (os_fvtxs_qy2*os_fvtxn_qy2) ) / ( os_fvtxs_qw*os_fvtxn_qw );
      float os_fvtxsfvtxn_qq3 = ( (os_fvtxs_qx3*os_fvtxn_qx3) + (os_fvtxs_qy3*os_fvtxn_qy3) ) / ( os_fvtxs_qw*os_fvtxn_qw );
      os_fvtxsfvtxn_c22->Fill(0.0,os_fvtxsfvtxn_qq2);
      os_fvtxsfvtxn_c32->Fill(0.0,os_fvtxsfvtxn_qq3);
      npc1_os_fvtxsfvtxn_c22->Fill(npc1,os_fvtxsfvtxn_qq2);
      npc1_os_fvtxsfvtxn_c32->Fill(npc1,os_fvtxsfvtxn_qq3);
      nfvtxt_os_fvtxsfvtxn_c22->Fill(nfvtxt,os_fvtxsfvtxn_qq2); // see below
      nfvtxc_os_fvtxsfvtxn_c22->Fill(nfvtxc,os_fvtxsfvtxn_qq2); // see below
      nfvtxt_os_fvtxsfvtxn_c32->Fill(nfvtxt,os_fvtxsfvtxn_qq3);
      nfvtxc_os_fvtxsfvtxn_c32->Fill(nfvtxc,os_fvtxsfvtxn_qq3);
      // ---

      float os_fvtxsfvtxn_tracks_qq2 = ( (os_fvtxs_tracks_qx2*os_fvtxn_tracks_qx2) + (os_fvtxs_tracks_qy2*os_fvtxn_tracks_qy2) ) / ( os_fvtxs_tracks_qw*os_fvtxn_tracks_qw );
      float os_fvtxsfvtxn_tracks_qq3 = ( (os_fvtxs_tracks_qx3*os_fvtxn_tracks_qx3) + (os_fvtxs_tracks_qy3*os_fvtxn_tracks_qy3) ) / ( os_fvtxs_tracks_qw*os_fvtxn_tracks_qw );
      os_fvtxsfvtxn_tracks_c22->Fill(0.0,os_fvtxsfvtxn_tracks_qq2);
      os_fvtxsfvtxn_tracks_c32->Fill(0.0,os_fvtxsfvtxn_tracks_qq3);
      npc1_os_fvtxsfvtxn_tracks_c22->Fill(npc1,os_fvtxsfvtxn_tracks_qq2);
      npc1_os_fvtxsfvtxn_tracks_c32->Fill(npc1,os_fvtxsfvtxn_tracks_qq3);
      nfvtxt_os_fvtxsfvtxn_tracks_c22->Fill(nfvtxt,os_fvtxsfvtxn_tracks_qq2);
      nfvtxc_os_fvtxsfvtxn_tracks_c22->Fill(nfvtxc,os_fvtxsfvtxn_tracks_qq2);
      nfvtxt_os_fvtxsfvtxn_tracks_c32->Fill(nfvtxt,os_fvtxsfvtxn_tracks_qq3);
      nfvtxc_os_fvtxsfvtxn_tracks_c32->Fill(nfvtxc,os_fvtxsfvtxn_tracks_qq3);

      // --- now have a look at some 4 particle cumulants
      if ( good_4_event )
	{
	  float os_bbc_qqqq4 = calc4_event(os_bbc_qx2,os_bbc_qy2,os_bbc_qx4,os_bbc_qy4,os_bbc_qw);
	  float os_fvtxs_qqqq4 = calc4_event(os_fvtxs_qx2,os_fvtxs_qy2,os_fvtxs_qx4,os_fvtxs_qy4,os_fvtxs_qw);
	  float os_fvtxn_qqqq4 = calc4_event(os_fvtxn_qx2,os_fvtxn_qy2,os_fvtxn_qx4,os_fvtxn_qy4,os_fvtxn_qw);
	  float os_fvtxc_qqqq4 = calc4_event(os_fvtxc_qx2,os_fvtxc_qy2,os_fvtxc_qx4,os_fvtxc_qy4,os_fvtxc_qw);
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
	  // ---
	  //nfvtxt_os_fvtxs_c22->Fill(nfvtxt,os_fvtxs_qq2);
	  //nfvtxc_os_fvtxs_c22->Fill(nfvtxc,os_fvtxs_qq2);
	  //nfvtxt_os_fvtxn_c22->Fill(nfvtxt,os_fvtxn_qq2);
	  //nfvtxc_os_fvtxn_c22->Fill(nfvtxc,os_fvtxn_qq2);
	  nfvtxt_os_fvtxs_c24->Fill(nfvtxt,os_fvtxs_qqqq4);
	  nfvtxc_os_fvtxs_c24->Fill(nfvtxc,os_fvtxs_qqqq4);
	  nfvtxt_os_fvtxn_c24->Fill(nfvtxt,os_fvtxn_qqqq4);
	  nfvtxc_os_fvtxn_c24->Fill(nfvtxc,os_fvtxn_qqqq4);
	  nfvtxt_os_fvtxc_c24->Fill(nfvtxt,os_fvtxc_qqqq4);
	  nfvtxc_os_fvtxc_c24->Fill(nfvtxc,os_fvtxc_qqqq4);
	  //nfvtxt_os_fvtxsfvtxn_c22->Fill(nfvtxt,os_fvtxsfvtxn_qq2);
	  //nfvtxc_os_fvtxsfvtxn_c22->Fill(nfvtxc,os_fvtxsfvtxn_qq2);
	  nfvtxt_os_fvtxsfvtxn_c24a->Fill(nfvtxt,os_fvtxsfvtxn_qq2*os_fvtxsfvtxn_qq2); // doesn't account for cross terms
	  nfvtxc_os_fvtxsfvtxn_c24a->Fill(nfvtxc,os_fvtxsfvtxn_qq2*os_fvtxsfvtxn_qq2); // doesn't account for cross terms
	  nfvtxt_os_fvtxsfvtxn_c24b->Fill(nfvtxt,os_fvtxs_qq2*os_fvtxs_qq2*os_fvtxn_qq2*os_fvtxn_qq2); // i think this power counts to v^8
	  nfvtxc_os_fvtxsfvtxn_c24b->Fill(nfvtxc,os_fvtxs_qq2*os_fvtxs_qq2*os_fvtxn_qq2*os_fvtxn_qq2); // i think this power counts to v^8
	  nfvtxt_os_fvtxsfvtxn_c24c->Fill(nfvtxt,os_fvtxsfvtxn_qq2*os_fvtxs_qq2*os_fvtxn_qq2); // i think this power counts to v^6
	  nfvtxc_os_fvtxsfvtxn_c24c->Fill(nfvtxc,os_fvtxsfvtxn_qq2*os_fvtxs_qq2*os_fvtxn_qq2); // i think this power counts to v^6
	  nfvtxt_os_fvtxsfvtxn_c24d->Fill(nfvtxt,os_fvtxs_qq2*os_fvtxn_qq2); // i think this is right (but problem with autocorrelations)
	  nfvtxc_os_fvtxsfvtxn_c24d->Fill(nfvtxc,os_fvtxs_qq2*os_fvtxn_qq2); // i think this is right (but problem with autocorrelations)
	  nfvtxt_os_fvtxsfvtxn_c24->Fill(nfvtxt,os_fvtxs_qq2*os_fvtxn_qq2); // might as well use this as an anchor
	  nfvtxc_os_fvtxsfvtxn_c24->Fill(nfvtxc,os_fvtxs_qq2*os_fvtxn_qq2); // might as well use this as an anchor

	  float os_fvtxs_tracks_qqqq4 = calc4_event(os_fvtxs_tracks_qx2,os_fvtxs_tracks_qy2,os_fvtxs_tracks_qx4,os_fvtxs_tracks_qy4,os_fvtxs_tracks_qw);
	  float os_fvtxn_tracks_qqqq4 = calc4_event(os_fvtxn_tracks_qx2,os_fvtxn_tracks_qy2,os_fvtxn_tracks_qx4,os_fvtxn_tracks_qy4,os_fvtxn_tracks_qw);
	  float os_fvtxc_tracks_qqqq4 = calc4_event(os_fvtxc_tracks_qx2,os_fvtxc_tracks_qy2,os_fvtxc_tracks_qx4,os_fvtxc_tracks_qy4,os_fvtxc_tracks_qw);
	  nfvtxt_os_fvtxs_tracks_c24->Fill(nfvtxt_south,os_fvtxs_tracks_qqqq4);
	  nfvtxc_os_fvtxs_tracks_c24->Fill(nfvtxc_south,os_fvtxs_tracks_qqqq4);
	  nfvtxt_os_fvtxn_tracks_c24->Fill(nfvtxt_north,os_fvtxn_tracks_qqqq4);
	  nfvtxc_os_fvtxn_tracks_c24->Fill(nfvtxc_north,os_fvtxn_tracks_qqqq4);
	  nfvtxt_os_fvtxc_tracks_c24->Fill(nfvtxt,os_fvtxc_tracks_qqqq4);
	  nfvtxc_os_fvtxc_tracks_c24->Fill(nfvtxc,os_fvtxc_tracks_qqqq4);

	  // cout << "tracks cumulant " << os_fvtxs_tracks_qq2 << endl;
	  // cout << "tracks qx " << os_fvtxs_tracks_qx2 << endl;
	  // cout << "tracks qy " << os_fvtxs_tracks_qy2 << endl;
	  // cout << "tracks qw " << os_fvtxs_tracks_qw << endl;
	  // nfvtxt_os_fvtxs_tracks_c22->Fill(nfvtxt,os_fvtxs_tracks_qq2);
	  // nfvtxc_os_fvtxs_tracks_c22->Fill(nfvtxc,os_fvtxs_tracks_qq2);
	  // nfvtxt_os_fvtxn_tracks_c22->Fill(nfvtxt,os_fvtxn_tracks_qq2);
	  // nfvtxc_os_fvtxn_tracks_c22->Fill(nfvtxc,os_fvtxn_tracks_qq2);
	  // nfvtxt_os_fvtxs_tracks_c24->Fill(nfvtxt,os_fvtxs_tracks_qqqq4);
	  // nfvtxc_os_fvtxs_tracks_c24->Fill(nfvtxc,os_fvtxs_tracks_qqqq4);
	  // nfvtxt_os_fvtxn_tracks_c24->Fill(nfvtxt,os_fvtxn_tracks_qqqq4);
	  // nfvtxc_os_fvtxn_tracks_c24->Fill(nfvtxc,os_fvtxn_tracks_qqqq4);
	  //nfvtxt_os_fvtxsfvtxn_tracks_c22->Fill(nfvtxt,os_fvtxsfvtxn_tracks_qq2);
	  //nfvtxc_os_fvtxsfvtxn_tracks_c22->Fill(nfvtxc,os_fvtxsfvtxn_tracks_qq2);
	  nfvtxt_os_fvtxsfvtxn_tracks_c24a->Fill(nfvtxt,os_fvtxsfvtxn_tracks_qq2*os_fvtxsfvtxn_tracks_qq2); // doesn't account for cross terms
	  nfvtxc_os_fvtxsfvtxn_tracks_c24a->Fill(nfvtxc,os_fvtxsfvtxn_tracks_qq2*os_fvtxsfvtxn_tracks_qq2); // doesn't account for cross terms
	  nfvtxt_os_fvtxsfvtxn_tracks_c24b->Fill(nfvtxt,os_fvtxs_tracks_qq2*os_fvtxs_tracks_qq2*os_fvtxn_tracks_qq2*os_fvtxn_tracks_qq2); // v^8 ?
	  nfvtxc_os_fvtxsfvtxn_tracks_c24b->Fill(nfvtxc,os_fvtxs_tracks_qq2*os_fvtxs_tracks_qq2*os_fvtxn_tracks_qq2*os_fvtxn_tracks_qq2); // v^8 ?
	  nfvtxt_os_fvtxsfvtxn_tracks_c24c->Fill(nfvtxt,os_fvtxsfvtxn_tracks_qq2*os_fvtxs_tracks_qq2*os_fvtxn_tracks_qq2); // v^6 ?
	  nfvtxc_os_fvtxsfvtxn_tracks_c24c->Fill(nfvtxc,os_fvtxsfvtxn_tracks_qq2*os_fvtxs_tracks_qq2*os_fvtxn_tracks_qq2); // v^6 ?
	  nfvtxt_os_fvtxsfvtxn_tracks_c24d->Fill(nfvtxt,os_fvtxs_tracks_qq2*os_fvtxn_tracks_qq2); // i think this is right (and tracks shouldn't have autocorrelation issues)
	  nfvtxc_os_fvtxsfvtxn_tracks_c24d->Fill(nfvtxc,os_fvtxs_tracks_qq2*os_fvtxn_tracks_qq2); // i think this is right (and tracks shouldn't have autocorrelation issues)
	  nfvtxt_os_fvtxsfvtxn_tracks_c24->Fill(nfvtxt,os_fvtxs_tracks_qq2*os_fvtxn_tracks_qq2); // use this as an anchor
	  nfvtxc_os_fvtxsfvtxn_tracks_c24->Fill(nfvtxc,os_fvtxs_tracks_qq2*os_fvtxn_tracks_qq2); // use this as an anchor

	  // --- now let's look at some subevent 4pc stuff
	  TComplex tc_fvtxs_ce0(os_fvtxs_ce0_qx2,os_fvtxs_ce0_qy2);
	  TComplex tc_fvtxs_ce1(os_fvtxs_ce1_qx2,os_fvtxs_ce1_qy2);
	  TComplex tc_fvtxn_ce0(os_fvtxn_ce0_qx2,os_fvtxn_ce0_qy2);
	  TComplex tc_fvtxn_ce1(os_fvtxn_ce1_qx2,os_fvtxn_ce1_qy2);
	  // --- first
	  TComplex tc_fvtxns_ce01 = tc_fvtxs_ce0 * tc_fvtxs_ce1 * TComplex::Conjugate(tc_fvtxn_ce0) * TComplex::Conjugate(tc_fvtxn_ce1);
	  float os_fvtxns_ce01_qqqq4 = tc_fvtxns_ce01.Re();
	  float norm = os_fvtxs_ce0_qw*os_fvtxs_ce1_qw*os_fvtxn_ce0_qw*os_fvtxn_ce1_qw;
	  os_fvtxns_ce01_qqqq4 /= norm;
	  nfvtxt_os_fvtxsfvtxn_ce01_c24->Fill(nfvtxt,os_fvtxns_ce01_qqqq4); // seems redundant but i sort of need an anchor result
	  nfvtxc_os_fvtxsfvtxn_ce01_c24->Fill(nfvtxc,os_fvtxns_ce01_qqqq4); // seems redundant but i sort of need an anchor result
	  nfvtxt_os_fvtxsfvtxn_ce01_c24a->Fill(nfvtxt,os_fvtxns_ce01_qqqq4);
	  nfvtxc_os_fvtxsfvtxn_ce01_c24a->Fill(nfvtxc,os_fvtxns_ce01_qqqq4);
	  // --- second
	  tc_fvtxns_ce01 = tc_fvtxs_ce0 * TComplex::Conjugate(tc_fvtxs_ce1) * tc_fvtxn_ce0 * TComplex::Conjugate(tc_fvtxn_ce1);
	  os_fvtxns_ce01_qqqq4 = tc_fvtxns_ce01.Re();
	  os_fvtxns_ce01_qqqq4 /= norm;
	  nfvtxt_os_fvtxsfvtxn_ce01_c24b->Fill(nfvtxt,os_fvtxns_ce01_qqqq4);
	  nfvtxc_os_fvtxsfvtxn_ce01_c24b->Fill(nfvtxc,os_fvtxns_ce01_qqqq4);
	  // --- third
	  tc_fvtxns_ce01 = TComplex::Conjugate(tc_fvtxs_ce0) * tc_fvtxs_ce1 * TComplex::Conjugate(tc_fvtxn_ce0) * tc_fvtxn_ce1;
	  os_fvtxns_ce01_qqqq4 = tc_fvtxns_ce01.Re();
	  os_fvtxns_ce01_qqqq4 /= norm;
	  nfvtxt_os_fvtxsfvtxn_ce01_c24c->Fill(nfvtxt,os_fvtxns_ce01_qqqq4);
	  nfvtxc_os_fvtxsfvtxn_ce01_c24c->Fill(nfvtxc,os_fvtxns_ce01_qqqq4);
	  // --- fourth
	  tc_fvtxns_ce01 = TComplex::Conjugate(tc_fvtxs_ce0) * TComplex::Conjugate(tc_fvtxs_ce1) * tc_fvtxn_ce0 * tc_fvtxn_ce1;
	  os_fvtxns_ce01_qqqq4 = tc_fvtxns_ce01.Re();
	  os_fvtxns_ce01_qqqq4 /= norm;
	  nfvtxt_os_fvtxsfvtxn_ce01_c24d->Fill(nfvtxt,os_fvtxns_ce01_qqqq4);
	  nfvtxc_os_fvtxsfvtxn_ce01_c24d->Fill(nfvtxc,os_fvtxns_ce01_qqqq4);
	}

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
      float bbc_south_psi3_docalib;

      float fvtx_south_psi2_docalib;
      float fvtx0_south_psi2_docalib;
      float fvtx1_south_psi2_docalib;
      float fvtx2_south_psi2_docalib;
      float fvtx3_south_psi2_docalib;

      float fvtx_south_psi3_docalib;
      float fvtx0_south_psi3_docalib;
      float fvtx1_south_psi3_docalib;
      float fvtx2_south_psi3_docalib;
      float fvtx3_south_psi3_docalib;

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

      float bbc_south_psi2_dcnw;
      float bbc_south_psi3_dcnw;

      float fvtx_south_psi2_dcnw;
      float fvtx0_south_psi2_dcnw;
      float fvtx1_south_psi2_dcnw;
      float fvtx2_south_psi2_dcnw;
      float fvtx3_south_psi2_dcnw;

      float fvtx_south_psi3_dcnw;
      float fvtx0_south_psi3_dcnw;
      float fvtx1_south_psi3_dcnw;
      float fvtx2_south_psi3_dcnw;
      float fvtx3_south_psi3_dcnw;

      float fvtx_north_psi2_dcnw;
      float fvtx0_north_psi2_dcnw;
      float fvtx1_north_psi2_dcnw;
      float fvtx2_north_psi2_dcnw;
      float fvtx3_north_psi2_dcnw;

      float fvtx_north_psi3_dcnw;
      float fvtx0_north_psi3_dcnw;
      float fvtx1_north_psi3_dcnw;
      float fvtx2_north_psi3_dcnw;
      float fvtx3_north_psi3_dcnw;

      float fvtx_south_psi2_tracks;
      float fvtx0_south_psi2_tracks;
      float fvtx1_south_psi2_tracks;

      float fvtx_south_psi3_tracks;
      float fvtx0_south_psi3_tracks;
      float fvtx1_south_psi3_tracks;

      float fvtx_north_psi2_tracks;
      float fvtx0_north_psi2_tracks;
      float fvtx1_north_psi2_tracks;

      float fvtx_north_psi3_tracks;
      float fvtx0_north_psi3_tracks;
      float fvtx1_north_psi3_tracks;

      if ( bbc_pmts )
        {
          bbc_south_psi2_docalib = (sumxy[1][bbcs_index][2]>0)?sumxy[1][bbcs_index][3]:-9999.9;
          bbc_south_psi3_docalib = (sumxy[2][bbcs_index][2]>0)?sumxy[2][bbcs_index][3]:-9999.9;
          bbc_south_psi2_dcnw = (sumxy[1][bbcs_nw_index][2]>0)?sumxy[1][bbcs_nw_index][3]:-9999.9;
          bbc_south_psi3_dcnw = (sumxy[2][bbcs_nw_index][2]>0)?sumxy[2][bbcs_nw_index][3]:-9999.9;
        }
      if ( fvtx_clusters )
        {
          // ---
          fvtx_south_psi2_docalib  = (sumxy[1][fvtxs_index][2]>4)  ? sumxy[1][fvtxs_index][3]  : -9999.9;
          fvtx0_south_psi2_docalib = (sumxy[1][fvtxs0_index][2]>4) ? sumxy[1][fvtxs0_index][3] : -9999.9;
          fvtx1_south_psi2_docalib = (sumxy[1][fvtxs1_index][2]>4) ? sumxy[1][fvtxs1_index][3] : -9999.9;
          fvtx2_south_psi2_docalib = (sumxy[1][fvtxs2_index][2]>4) ? sumxy[1][fvtxs2_index][3] : -9999.9;
          fvtx3_south_psi2_docalib = (sumxy[1][fvtxs3_index][2]>4) ? sumxy[1][fvtxs3_index][3] : -9999.9;
          fvtx_south_psi3_docalib  = (sumxy[2][fvtxs_index][2]>4)  ? sumxy[2][fvtxs_index][3]  : -9999.9;
          fvtx0_south_psi3_docalib = (sumxy[2][fvtxs0_index][2]>4) ? sumxy[2][fvtxs0_index][3] : -9999.9;
          fvtx1_south_psi3_docalib = (sumxy[2][fvtxs1_index][2]>4) ? sumxy[2][fvtxs1_index][3] : -9999.9;
          fvtx2_south_psi3_docalib = (sumxy[2][fvtxs2_index][2]>4) ? sumxy[2][fvtxs2_index][3] : -9999.9;
          fvtx3_south_psi3_docalib = (sumxy[2][fvtxs3_index][2]>4) ? sumxy[2][fvtxs3_index][3] : -9999.9;
          // ---
          fvtx_north_psi2_docalib  = (sumxy[1][fvtxn_index][2]>4)  ? sumxy[1][fvtxn_index][3]  : -9999.9;
          fvtx0_north_psi2_docalib = (sumxy[1][fvtxn0_index][2]>4) ? sumxy[1][fvtxn0_index][3] : -9999.9;
          fvtx1_north_psi2_docalib = (sumxy[1][fvtxn1_index][2]>4) ? sumxy[1][fvtxn1_index][3] : -9999.9;
          fvtx2_north_psi2_docalib = (sumxy[1][fvtxn2_index][2]>4) ? sumxy[1][fvtxn2_index][3] : -9999.9;
          fvtx3_north_psi2_docalib = (sumxy[1][fvtxn3_index][2]>4) ? sumxy[1][fvtxn3_index][3] : -9999.9;
          fvtx_north_psi3_docalib  = (sumxy[2][fvtxn_index][2]>4)  ? sumxy[2][fvtxn_index][3]  : -9999.9;
          fvtx0_north_psi3_docalib = (sumxy[2][fvtxn0_index][2]>4) ? sumxy[2][fvtxn0_index][3] : -9999.9;
          fvtx1_north_psi3_docalib = (sumxy[2][fvtxn1_index][2]>4) ? sumxy[2][fvtxn1_index][3] : -9999.9;
          fvtx2_north_psi3_docalib = (sumxy[2][fvtxn2_index][2]>4) ? sumxy[2][fvtxn2_index][3] : -9999.9;
          fvtx3_north_psi3_docalib = (sumxy[2][fvtxn3_index][2]>4) ? sumxy[2][fvtxn3_index][3] : -9999.9;
          // ---
          fvtx_south_psi2_dcnw  = (sumxy[1][fvtxs_nw_index][2]>4)  ? sumxy[1][fvtxs_nw_index][3]  : -9999.9;
          fvtx0_south_psi2_dcnw = (sumxy[1][fvtxs0_nw_index][2]>4) ? sumxy[1][fvtxs0_nw_index][3] : -9999.9;
          fvtx1_south_psi2_dcnw = (sumxy[1][fvtxs1_nw_index][2]>4) ? sumxy[1][fvtxs1_nw_index][3] : -9999.9;
          fvtx2_south_psi2_dcnw = (sumxy[1][fvtxs2_nw_index][2]>4) ? sumxy[1][fvtxs2_nw_index][3] : -9999.9;
          fvtx3_south_psi2_dcnw = (sumxy[1][fvtxs3_nw_index][2]>4) ? sumxy[1][fvtxs3_nw_index][3] : -9999.9;
          fvtx_south_psi3_dcnw  = (sumxy[2][fvtxs_nw_index][2]>4)  ? sumxy[2][fvtxs_nw_index][3]  : -9999.9;
          fvtx0_south_psi3_dcnw = (sumxy[2][fvtxs0_nw_index][2]>4) ? sumxy[2][fvtxs0_nw_index][3] : -9999.9;
          fvtx1_south_psi3_dcnw = (sumxy[2][fvtxs1_nw_index][2]>4) ? sumxy[2][fvtxs1_nw_index][3] : -9999.9;
          fvtx2_south_psi3_dcnw = (sumxy[2][fvtxs2_nw_index][2]>4) ? sumxy[2][fvtxs2_nw_index][3] : -9999.9;
          fvtx3_south_psi3_dcnw = (sumxy[2][fvtxs3_nw_index][2]>4) ? sumxy[2][fvtxs3_nw_index][3] : -9999.9;
          // ---
          fvtx_north_psi2_dcnw  = (sumxy[1][fvtxn_nw_index][2]>4)  ? sumxy[1][fvtxn_nw_index][3]  : -9999.9;
          fvtx0_north_psi2_dcnw = (sumxy[1][fvtxn0_nw_index][2]>4) ? sumxy[1][fvtxn0_nw_index][3] : -9999.9;
          fvtx1_north_psi2_dcnw = (sumxy[1][fvtxn1_nw_index][2]>4) ? sumxy[1][fvtxn1_nw_index][3] : -9999.9;
          fvtx2_north_psi2_dcnw = (sumxy[1][fvtxn2_nw_index][2]>4) ? sumxy[1][fvtxn2_nw_index][3] : -9999.9;
          fvtx3_north_psi2_dcnw = (sumxy[1][fvtxn3_nw_index][2]>4) ? sumxy[1][fvtxn3_nw_index][3] : -9999.9;
          fvtx_north_psi3_dcnw  = (sumxy[2][fvtxn_nw_index][2]>4)  ? sumxy[2][fvtxn_nw_index][3]  : -9999.9;
          fvtx0_north_psi3_dcnw = (sumxy[2][fvtxn0_nw_index][2]>4) ? sumxy[2][fvtxn0_nw_index][3] : -9999.9;
          fvtx1_north_psi3_dcnw = (sumxy[2][fvtxn1_nw_index][2]>4) ? sumxy[2][fvtxn1_nw_index][3] : -9999.9;
          fvtx2_north_psi3_dcnw = (sumxy[2][fvtxn2_nw_index][2]>4) ? sumxy[2][fvtxn2_nw_index][3] : -9999.9;
          fvtx3_north_psi3_dcnw = (sumxy[2][fvtxn3_nw_index][2]>4) ? sumxy[2][fvtxn3_nw_index][3] : -9999.9;
          // --- tracks don't belong with clusters I guess, but...
          fvtx_south_psi2_tracks  = (sumxy[1][fvtxs_tracks_index][2]>4)  ? sumxy[1][fvtxs_tracks_index][3]  : -9999.9;
          fvtx0_south_psi2_tracks = (sumxy[1][fvtxs0_tracks_index][2]>4) ? sumxy[1][fvtxs0_tracks_index][3] : -9999.9;
          fvtx1_south_psi2_tracks = (sumxy[1][fvtxs1_tracks_index][2]>4) ? sumxy[1][fvtxs1_tracks_index][3] : -9999.9;
          fvtx_south_psi3_tracks  = (sumxy[2][fvtxs_tracks_index][2]>4)  ? sumxy[2][fvtxs_tracks_index][3]  : -9999.9;
          fvtx0_south_psi3_tracks = (sumxy[2][fvtxs0_tracks_index][2]>4) ? sumxy[2][fvtxs0_tracks_index][3] : -9999.9;
          fvtx1_south_psi3_tracks = (sumxy[2][fvtxs1_tracks_index][2]>4) ? sumxy[2][fvtxs1_tracks_index][3] : -9999.9;
          // ---
          fvtx_north_psi2_tracks  = (sumxy[1][fvtxn_tracks_index][2]>4)  ? sumxy[1][fvtxn_tracks_index][3]  : -9999.9;
          fvtx0_north_psi2_tracks = (sumxy[1][fvtxn0_tracks_index][2]>4) ? sumxy[1][fvtxn0_tracks_index][3] : -9999.9;
          fvtx1_north_psi2_tracks = (sumxy[1][fvtxn1_tracks_index][2]>4) ? sumxy[1][fvtxn1_tracks_index][3] : -9999.9;
          fvtx_north_psi3_tracks  = (sumxy[2][fvtxn_tracks_index][2]>4)  ? sumxy[2][fvtxn_tracks_index][3]  : -9999.9;
          fvtx0_north_psi3_tracks = (sumxy[2][fvtxn0_tracks_index][2]>4) ? sumxy[2][fvtxn0_tracks_index][3] : -9999.9;
          fvtx1_north_psi3_tracks = (sumxy[2][fvtxn1_tracks_index][2]>4) ? sumxy[2][fvtxn1_tracks_index][3] : -9999.9;
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

      // --- FVTX clusters and tracks
      tp1f_reso2_FVTX_FVTX->Fill(0.0,cos(2*(fvtx_south_psi2_docalib-fvtx_south_psi2_tracks)));
      tp1f_reso3_FVTX_FVTX->Fill(0.0,cos(3*(fvtx_south_psi3_docalib-fvtx_south_psi3_tracks)));
      th1d_reso2_FVTX_FVTX->Fill(cos(2*(fvtx_south_psi2_docalib-fvtx_south_psi2_tracks)));
      th1d_reso3_FVTX_FVTX->Fill(cos(3*(fvtx_south_psi3_docalib-fvtx_south_psi3_tracks)));
      th1d_dreso2_FVTX_FVTX->Fill(fvtx_south_psi2_docalib-fvtx_south_psi2_tracks);
      th1d_dreso3_FVTX_FVTX->Fill(fvtx_south_psi3_docalib-fvtx_south_psi3_tracks);
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

      // --- for no weight histograms
      tp1f_nw_reso2_BBC_FVTX->Fill(0.0,cos(2*(bbc_south_psi2_dcnw-fvtx_south_psi2_dcnw)));
      tp1f_nw_reso2_BBC_FVTXN->Fill(0.0,cos(2*(bbc_south_psi2_dcnw-fvtx_north_psi2_dcnw)));
      tp1f_nw_reso2_FVTXS_FVTXN->Fill(0.0,cos(2*(fvtx_south_psi2_dcnw-fvtx_north_psi2_dcnw)));
      tp1f_nw_reso3_BBC_FVTX->Fill(0.0,cos(3*(bbc_south_psi3_dcnw-fvtx_south_psi3_dcnw)));
      tp1f_nw_reso3_BBC_FVTXN->Fill(0.0,cos(3*(bbc_south_psi3_dcnw-fvtx_north_psi3_dcnw)));
      tp1f_nw_reso3_FVTXS_FVTXN->Fill(0.0,cos(3*(fvtx_south_psi3_dcnw-fvtx_north_psi3_dcnw)));
      tp1f_nw_reso42_BBC_FVTX->Fill(0.0,cos(4*(bbc_south_psi2_dcnw-fvtx_south_psi2_dcnw)));
      tp1f_nw_reso42_BBC_FVTXN->Fill(0.0,cos(4*(bbc_south_psi2_dcnw-fvtx_north_psi2_dcnw)));
      tp1f_nw_reso42_FVTXS_FVTXN->Fill(0.0,cos(4*(fvtx_south_psi2_dcnw-fvtx_north_psi2_dcnw)));

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

	      // --- rotation now done in trees
              float phi_angle = phi0;
              float pt_angle = pt;
              double theta = atan2(pt,pz);
              double eta = -log(tan(theta/2));
              //cout << "pz is " << pz << " and eta is " << eta << endl;

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

              float os_fvtxs_tracks_uq2 = ( (ux2*os_fvtxs_tracks_qx2) + (uy2*os_fvtxs_tracks_qy2) ) / ( os_fvtxs_tracks_qw );
              float os_fvtxs_tracks_uq3 = ( (ux3*os_fvtxs_tracks_qx3) + (uy3*os_fvtxs_tracks_qy3) ) / ( os_fvtxs_tracks_qw );
              float os_fvtxn_tracks_uq2 = ( (ux2*os_fvtxn_tracks_qx2) + (uy2*os_fvtxn_tracks_qy2) ) / ( os_fvtxn_tracks_qw );
              float os_fvtxn_tracks_uq3 = ( (ux3*os_fvtxn_tracks_qx3) + (uy3*os_fvtxn_tracks_qy3) ) / ( os_fvtxn_tracks_qw );

              os_fvtxs_tracks_d22_both->Fill(pt_angle,os_fvtxs_tracks_uq2);
              if ( dcarm == 1 ) os_fvtxs_tracks_d22_west->Fill(pt_angle,os_fvtxs_tracks_uq2);
              if ( dcarm == 0 ) os_fvtxs_tracks_d22_east->Fill(pt_angle,os_fvtxs_tracks_uq2);
              os_fvtxn_tracks_d22_both->Fill(pt_angle,os_fvtxn_tracks_uq2);
              if ( dcarm == 1 ) os_fvtxn_tracks_d22_west->Fill(pt_angle,os_fvtxn_tracks_uq2);
              if ( dcarm == 0 ) os_fvtxn_tracks_d22_east->Fill(pt_angle,os_fvtxn_tracks_uq2);

              os_fvtxs_tracks_d32_both->Fill(pt_angle,os_fvtxs_tracks_uq3);
              if ( dcarm == 1 ) os_fvtxs_tracks_d32_west->Fill(pt_angle,os_fvtxs_tracks_uq3);
              if ( dcarm == 0 ) os_fvtxs_tracks_d32_east->Fill(pt_angle,os_fvtxs_tracks_uq3);
              os_fvtxn_tracks_d32_both->Fill(pt_angle,os_fvtxn_tracks_uq3);
              if ( dcarm == 1 ) os_fvtxn_tracks_d32_west->Fill(pt_angle,os_fvtxn_tracks_uq3);
              if ( dcarm == 0 ) os_fvtxn_tracks_d32_east->Fill(pt_angle,os_fvtxn_tracks_uq3);

              // ------------------------------------------------------------

              npc1_os_cntbbcs_c22->Fill(npc1,os_bbc_uq2);
              npc1_os_cntfvtxs_c22->Fill(npc1,os_fvtxs_uq2);
              npc1_os_cntfvtxn_c22->Fill(npc1,os_fvtxn_uq2);
              npc1_os_cntbbcs_c32->Fill(npc1,os_bbc_uq3);
              npc1_os_cntfvtxs_c32->Fill(npc1,os_fvtxs_uq3);
              npc1_os_cntfvtxn_c32->Fill(npc1,os_fvtxn_uq3);
              nfvtxt_os_cntbbcs_c22->Fill(nfvtxt,os_bbc_uq2);
              nfvtxc_os_cntbbcs_c22->Fill(nfvtxc,os_bbc_uq2);
              nfvtxt_os_cntfvtxs_c22->Fill(nfvtxt,os_fvtxs_uq2);
              nfvtxc_os_cntfvtxs_c22->Fill(nfvtxc,os_fvtxs_uq2);
              nfvtxt_os_cntfvtxn_c22->Fill(nfvtxt,os_fvtxn_uq2);
              nfvtxc_os_cntfvtxn_c22->Fill(nfvtxc,os_fvtxn_uq2);
              nfvtxt_os_cntbbcs_c32->Fill(nfvtxt,os_bbc_uq3);
              nfvtxc_os_cntbbcs_c32->Fill(nfvtxc,os_bbc_uq3);
              nfvtxt_os_cntfvtxs_c32->Fill(nfvtxt,os_fvtxs_uq3);
              nfvtxc_os_cntfvtxs_c32->Fill(nfvtxc,os_fvtxs_uq3);
              nfvtxt_os_cntfvtxn_c32->Fill(nfvtxt,os_fvtxn_uq3);
              nfvtxc_os_cntfvtxn_c32->Fill(nfvtxc,os_fvtxn_uq3);

              // ------------------------------------------------------------

              os_bbcs_d22eta_both->Fill(eta,os_bbc_uq2);
              if ( dcarm == 1 ) os_bbcs_d22eta_west->Fill(eta,os_bbc_uq2);
              if ( dcarm == 0 ) os_bbcs_d22eta_east->Fill(eta,os_bbc_uq2);
              os_fvtxs_d22eta_both->Fill(eta,os_fvtxs_uq2);
              if ( dcarm == 1 ) os_fvtxs_d22eta_west->Fill(eta,os_fvtxs_uq2);
              if ( dcarm == 0 ) os_fvtxs_d22eta_east->Fill(eta,os_fvtxs_uq2);
              os_fvtxs_d22eta_both->Fill(eta,os_fvtxn_uq2);
              if ( dcarm == 1 ) os_fvtxn_d22eta_west->Fill(eta,os_fvtxn_uq2);
              if ( dcarm == 0 ) os_fvtxn_d22eta_east->Fill(eta,os_fvtxn_uq2);

              os_bbcs_d32eta_both->Fill(eta,os_bbc_uq3);
              if ( dcarm == 1 ) os_bbcs_d32eta_west->Fill(eta,os_bbc_uq3);
              if ( dcarm == 0 ) os_bbcs_d32eta_east->Fill(eta,os_bbc_uq3);
              os_fvtxs_d32eta_both->Fill(eta,os_fvtxs_uq2);
              if ( dcarm == 1 ) os_fvtxs_d32eta_west->Fill(eta,os_fvtxs_uq3);
              if ( dcarm == 0 ) os_fvtxs_d32eta_east->Fill(eta,os_fvtxs_uq3);
              os_fvtxn_d32eta_both->Fill(eta,os_fvtxn_uq2);
              if ( dcarm == 1 ) os_fvtxn_d32eta_west->Fill(eta,os_fvtxn_uq3);
              if ( dcarm == 0 ) os_fvtxn_d32eta_east->Fill(eta,os_fvtxn_uq3);

              os_fvtxs_tracks_d22eta_both->Fill(eta,os_fvtxs_tracks_uq2);
              if ( dcarm == 1 ) os_fvtxs_tracks_d22eta_west->Fill(eta,os_fvtxs_tracks_uq2);
              if ( dcarm == 0 ) os_fvtxs_tracks_d22eta_east->Fill(eta,os_fvtxs_tracks_uq2);
              os_fvtxn_tracks_d22eta_both->Fill(eta,os_fvtxn_tracks_uq2);
              if ( dcarm == 1 ) os_fvtxn_tracks_d22eta_west->Fill(eta,os_fvtxn_tracks_uq2);
              if ( dcarm == 0 ) os_fvtxn_tracks_d22eta_east->Fill(eta,os_fvtxn_tracks_uq2);

              os_fvtxs_tracks_d32eta_both->Fill(eta,os_fvtxs_tracks_uq3);
              if ( dcarm == 1 ) os_fvtxs_tracks_d32eta_west->Fill(eta,os_fvtxs_tracks_uq3);
              if ( dcarm == 0 ) os_fvtxs_tracks_d32eta_east->Fill(eta,os_fvtxs_tracks_uq3);
              os_fvtxn_tracks_d32eta_both->Fill(eta,os_fvtxn_tracks_uq3);
              if ( dcarm == 1 ) os_fvtxn_tracks_d32eta_west->Fill(eta,os_fvtxn_tracks_uq3);
              if ( dcarm == 0 ) os_fvtxn_tracks_d32eta_east->Fill(eta,os_fvtxn_tracks_uq3);

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

              // ------------------------------------------------------------

              npc1_bbcs_v2->Fill(npc1,cos(2.0*(phi_angle-os_bbc_psi2)));
              npc1_fvtxs_v2->Fill(npc1,cos(2.0*(phi_angle-os_fvtxs_psi2)));
              npc1_fvtxn_v2->Fill(npc1,cos(2.0*(phi_angle-os_fvtxn_psi2)));
              npc1_bbcs_v3->Fill(npc1,cos(3.0*(phi_angle-os_bbc_psi3)));
              npc1_fvtxs_v3->Fill(npc1,cos(3.0*(phi_angle-os_fvtxs_psi3)));
              npc1_fvtxn_v3->Fill(npc1,cos(3.0*(phi_angle-os_fvtxn_psi3)));
              nfvtxt_bbcs_v2->Fill(nfvtxt,cos(2.0*(phi_angle-os_bbc_psi2)));
              nfvtxt_fvtxs_v2->Fill(nfvtxt,cos(2.0*(phi_angle-os_fvtxs_psi2)));
              nfvtxt_fvtxn_v2->Fill(nfvtxt,cos(2.0*(phi_angle-os_fvtxn_psi2)));
              nfvtxt_bbcs_v3->Fill(nfvtxt,cos(3.0*(phi_angle-os_bbc_psi3)));
              nfvtxt_fvtxs_v3->Fill(nfvtxt,cos(3.0*(phi_angle-os_fvtxs_psi3)));
              nfvtxt_fvtxn_v3->Fill(nfvtxt,cos(3.0*(phi_angle-os_fvtxn_psi3)));
              nfvtxc_bbcs_v2->Fill(nfvtxc,cos(2.0*(phi_angle-os_bbc_psi2)));
              nfvtxc_fvtxs_v2->Fill(nfvtxc,cos(2.0*(phi_angle-os_fvtxs_psi2)));
              nfvtxc_fvtxn_v2->Fill(nfvtxc,cos(2.0*(phi_angle-os_fvtxn_psi2)));
              nfvtxc_bbcs_v3->Fill(nfvtxc,cos(3.0*(phi_angle-os_bbc_psi3)));
              nfvtxc_fvtxs_v3->Fill(nfvtxc,cos(3.0*(phi_angle-os_fvtxs_psi3)));
              nfvtxc_fvtxn_v3->Fill(nfvtxc,cos(3.0*(phi_angle-os_fvtxn_psi3)));

              // ------------------------------------------------------------

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
                      // ---
                      bbcs_v2eta_both_docalib->Fill(eta,cosbbc_dphi2_docalib);
                      if ( dcarm == 1 ) bbcs_v2eta_west_docalib->Fill(eta,cosbbc_dphi2_docalib);
                      if ( dcarm == 0 ) bbcs_v2eta_east_docalib->Fill(eta,cosbbc_dphi2_docalib);
                      bbcs_v3eta_both_docalib->Fill(eta,cosbbc_dphi3_docalib);
                      if ( dcarm == 1 ) bbcs_v3eta_west_docalib->Fill(eta,cosbbc_dphi3_docalib);
                      if ( dcarm == 0 ) bbcs_v3eta_east_docalib->Fill(eta,cosbbc_dphi3_docalib);
                    } // check on bbc EP
                  // ---
                  if(-4.0<bbc_south_psi2_dcnw && bbc_south_psi2_dcnw<4.0)
                    {
                      // --- 2nd harmonic
                      double bbc_dphi2_dcnw = phi_angle - bbc_south_psi2_dcnw;
                      double cosbbc_dphi2_dcnw = TMath::Cos(2*bbc_dphi2_dcnw);
                      bbcs_v2_both_dcnw->Fill(pt_angle,cosbbc_dphi2_dcnw);
                      if ( dcarm == 1 ) bbcs_v2_west_dcnw->Fill(pt_angle,cosbbc_dphi2_dcnw);
                      if ( dcarm == 0 ) bbcs_v2_east_dcnw->Fill(pt_angle,cosbbc_dphi2_dcnw);
                      // --- 3rd harmonic
                      double bbc_dphi3_dcnw = phi_angle - bbc_south_psi3_dcnw;
                      double cosbbc_dphi3_dcnw = TMath::Cos(3*bbc_dphi3_dcnw);
                      bbcs_v3_both_dcnw->Fill(pt_angle,cosbbc_dphi3_dcnw);
                      if ( dcarm == 1 ) bbcs_v3_west_dcnw->Fill(pt_angle,cosbbc_dphi3_dcnw);
                      if ( dcarm == 0 ) bbcs_v3_east_dcnw->Fill(pt_angle,cosbbc_dphi3_dcnw);
                      // --- 4th harmonic
                      double fourfour = 4*phi_angle - 4*bbc_south_psi2_dcnw;
                      double cosfourfour = TMath::Cos(fourfour);
                      bbcs_v4_4Psi2_both_dcnw->Fill(pt_angle,cosfourfour);
                      if ( dcarm == 1 ) bbcs_v4_4Psi2_west_dcnw->Fill(pt_angle,cosfourfour);
                      if ( dcarm == 0 ) bbcs_v4_4Psi2_east_dcnw->Fill(pt_angle,cosfourfour);
                      // --- ep reso
                      tp1f_nw_reso2_BBC_CNT->Fill(0.0,cosbbc_dphi2_dcnw);
                      tp1f_nw_reso3_BBC_CNT->Fill(0.0,cosbbc_dphi3_dcnw);
                      tp1f_nw_reso42_BBC_CNT->Fill(0.0,cosfourfour);
                    } // check on EP
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
                  // ---
                  fvtxs_v2eta_both_docalib->Fill(eta,cosfvtx_dphi2_docalib);
                  if ( dcarm == 1 ) fvtxs_v2eta_west_docalib->Fill(eta,cosfvtx_dphi2_docalib);
                  if ( dcarm == 0 ) fvtxs_v2eta_east_docalib->Fill(eta,cosfvtx_dphi2_docalib);
                  fvtxs_v3eta_both_docalib->Fill(eta,cosfvtx_dphi3_docalib);
                  if ( dcarm == 1 ) fvtxs_v3eta_west_docalib->Fill(eta,cosfvtx_dphi3_docalib);
                  if ( dcarm == 0 ) fvtxs_v3eta_east_docalib->Fill(eta,cosfvtx_dphi3_docalib);
                }
              // --- now dcnw
              if(-4.0<fvtx_south_psi2_dcnw && fvtx_south_psi2_dcnw<4.0)
                {
                  // --- 2nd harmonic
                  double fvtx_dphi2_dcnw = phi_angle - fvtx_south_psi2_dcnw;
                  double cosfvtx_dphi2_dcnw = TMath::Cos(2*fvtx_dphi2_dcnw);
                  fvtxs_v2_both_dcnw->Fill(pt_angle,cosfvtx_dphi2_dcnw);
                  if ( dcarm == 1 ) fvtxs_v2_west_dcnw->Fill(pt_angle,cosfvtx_dphi2_dcnw);
                  if ( dcarm == 0 ) fvtxs_v2_east_dcnw->Fill(pt_angle,cosfvtx_dphi2_dcnw);
                  // --- 3rd harmonic
                  double fvtx_dphi3_dcnw = phi_angle - fvtx_south_psi3_dcnw;
                  double cosfvtx_dphi3_dcnw = TMath::Cos(3*fvtx_dphi3_dcnw);
                  fvtxs_v3_both_dcnw->Fill(pt_angle,cosfvtx_dphi3_dcnw);
                  if ( dcarm == 1 ) fvtxs_v3_west_dcnw->Fill(pt_angle,cosfvtx_dphi3_dcnw);
                  if ( dcarm == 0 ) fvtxs_v3_east_dcnw->Fill(pt_angle,cosfvtx_dphi3_dcnw);
                  // --- 4th harmonic
                  double fourfour = 4*phi_angle - 4*fvtx_south_psi2_dcnw;
                  double cosfourfour = TMath::Cos(fourfour);
                  fvtxs_v4_4Psi2_both_dcnw->Fill(pt_angle,cosfourfour);
                  if ( dcarm == 1 ) fvtxs_v4_4Psi2_west_dcnw->Fill(pt_angle,cosfourfour);
                  if ( dcarm == 0 ) fvtxs_v4_4Psi2_east_dcnw->Fill(pt_angle,cosfourfour);
                  // --- ep reso
                  tp1f_nw_reso2_CNT_FVTX->Fill(0.0,cosfvtx_dphi2_dcnw);
                  tp1f_nw_reso3_CNT_FVTX->Fill(0.0,cosfvtx_dphi3_dcnw);
                  tp1f_nw_reso42_CNT_FVTX->Fill(0.0,cosfourfour);
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
                  // ---
                  fvtxn_v2eta_both_docalib->Fill(eta,cosfvtx_dphi2_docalib);
                  if ( dcarm == 1 ) fvtxn_v2eta_west_docalib->Fill(eta,cosfvtx_dphi2_docalib);
                  if ( dcarm == 0 ) fvtxn_v2eta_east_docalib->Fill(eta,cosfvtx_dphi2_docalib);
                  fvtxn_v3eta_both_docalib->Fill(eta,cosfvtx_dphi3_docalib);
                  if ( dcarm == 1 ) fvtxn_v3eta_west_docalib->Fill(eta,cosfvtx_dphi3_docalib);
                  if ( dcarm == 0 ) fvtxn_v3eta_east_docalib->Fill(eta,cosfvtx_dphi3_docalib);
                }
              // --- now dcnw
              if(-4.0<fvtx_north_psi2_dcnw && fvtx_north_psi2_dcnw<4.0)
                {
                  // --- 2nd harmonic
                  double fvtx_dphi2_dcnw = phi_angle - fvtx_north_psi2_dcnw;
                  double cosfvtx_dphi2_dcnw = TMath::Cos(2*fvtx_dphi2_dcnw);
                  fvtxn_v2_both_dcnw->Fill(pt_angle,cosfvtx_dphi2_dcnw);
                  if ( dcarm == 1 ) fvtxn_v2_west_dcnw->Fill(pt_angle,cosfvtx_dphi2_dcnw);
                  if ( dcarm == 0 ) fvtxn_v2_east_dcnw->Fill(pt_angle,cosfvtx_dphi2_dcnw);
                  // --- 3rd harmonic
                  double fvtx_dphi3_dcnw = phi_angle - fvtx_north_psi3_dcnw;
                  double cosfvtx_dphi3_dcnw = TMath::Cos(3*fvtx_dphi3_dcnw);
                  fvtxn_v3_both_dcnw->Fill(pt_angle,cosfvtx_dphi3_dcnw);
                  if ( dcarm == 1 ) fvtxn_v3_west_dcnw->Fill(pt_angle,cosfvtx_dphi3_dcnw);
                  if ( dcarm == 0 ) fvtxn_v3_east_dcnw->Fill(pt_angle,cosfvtx_dphi3_dcnw);
                  // --- 4th harmonic
                  double fourfour = 4*phi_angle - 4*fvtx_north_psi2_dcnw;
                  double cosfourfour = TMath::Cos(fourfour);
                  fvtxn_v4_4Psi2_both_dcnw->Fill(pt_angle,cosfourfour);
                  if ( dcarm == 1 ) fvtxn_v4_4Psi2_west_dcnw->Fill(pt_angle,cosfourfour);
                  if ( dcarm == 0 ) fvtxn_v4_4Psi2_east_dcnw->Fill(pt_angle,cosfourfour);
                  // --- ep reso
                  tp1f_nw_reso2_CNT_FVTXN->Fill(0.0,cosfvtx_dphi2_dcnw);
                  tp1f_nw_reso3_CNT_FVTXN->Fill(0.0,cosfvtx_dphi3_dcnw);
                  tp1f_nw_reso42_CNT_FVTXN->Fill(0.0,cosfourfour);
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

              // --- now fvtx layers

              //fvtx layer 0
              if(-4.0<fvtx0_south_psi2_dcnw && fvtx0_south_psi2_dcnw<4.0)
                {
                  double fvtx_dphi2_dcnw = phi_angle - fvtx0_south_psi2_dcnw;
                  double cosfvtx_dphi2_dcnw = TMath::Cos(2*fvtx_dphi2_dcnw);
                  if ( dcarm == 1 ) fvtxs0_v2_west_dcnw->Fill(pt_angle,cosfvtx_dphi2_dcnw);
                  if ( dcarm == 0 ) fvtxs0_v2_east_dcnw->Fill(pt_angle,cosfvtx_dphi2_dcnw);
                  fvtxs0_v2_both_dcnw->Fill(pt_angle,cosfvtx_dphi2_dcnw);
                }

              //fvtx layer 1
              if(-4.0<fvtx1_south_psi2_dcnw && fvtx1_south_psi2_dcnw<4.0)
                {
                  double fvtx_dphi2_dcnw = phi_angle - fvtx1_south_psi2_dcnw;
                  double cosfvtx_dphi2_dcnw = TMath::Cos(2*fvtx_dphi2_dcnw);
                  if ( dcarm == 1 ) fvtxs1_v2_west_dcnw->Fill(pt_angle,cosfvtx_dphi2_dcnw);
                  if ( dcarm == 0 ) fvtxs1_v2_east_dcnw->Fill(pt_angle,cosfvtx_dphi2_dcnw);
                  fvtxs1_v2_both_dcnw->Fill(pt_angle,cosfvtx_dphi2_dcnw);
                }

              //fvtx layer 2
              if(-4.0<fvtx2_south_psi2_dcnw && fvtx2_south_psi2_dcnw<4.0)
                {
                  double fvtx_dphi2_dcnw = phi_angle - fvtx2_south_psi2_dcnw;
                  double cosfvtx_dphi2_dcnw = TMath::Cos(2*fvtx_dphi2_dcnw);
                  if ( dcarm == 1 ) fvtxs2_v2_west_dcnw->Fill(pt_angle,cosfvtx_dphi2_dcnw);
                  if ( dcarm == 0 ) fvtxs2_v2_east_dcnw->Fill(pt_angle,cosfvtx_dphi2_dcnw);
                  fvtxs2_v2_both_dcnw->Fill(pt_angle,cosfvtx_dphi2_dcnw);
                }

              //fvtx layer 3
              if(-4.0<fvtx3_south_psi2_dcnw && fvtx3_south_psi2_dcnw<4.0)
                {
                  double fvtx_dphi2_dcnw = phi_angle - fvtx3_south_psi2_dcnw;
                  double cosfvtx_dphi2_dcnw = TMath::Cos(2*fvtx_dphi2_dcnw);
                  if ( dcarm == 1 ) fvtxs3_v2_west_dcnw->Fill(pt_angle,cosfvtx_dphi2_dcnw);
                  if ( dcarm == 0 ) fvtxs3_v2_east_dcnw->Fill(pt_angle,cosfvtx_dphi2_dcnw);
                  fvtxs3_v2_both_dcnw->Fill(pt_angle,cosfvtx_dphi2_dcnw);
                } // check on psi

            } // loop over tracks

          for ( int i = 0; i < nfvtxt; ++i )
            {
              float phi = fphi[i];
              float eta = feta[i];
              int ns = 0;
              int nn = 0;
              if ( eta > 0 ) nn = 1;
              else if ( eta < 0 ) ns = 1;
              // cout << "eta is " << eta << " and ns is " << ns << " and nn is " << nn << endl;

              double cosbbc_dphi2_docalib = cos(2*(phi-bbc_south_psi2_docalib));
              bbcs_v2eta_both_docalib->Fill(eta,cosbbc_dphi2_docalib);
              bbcs_v2eta_west_docalib->Fill(eta,cosbbc_dphi2_docalib);
              bbcs_v2eta_east_docalib->Fill(eta,cosbbc_dphi2_docalib);
              double cosbbc_dphi3_docalib = cos(3*(phi-bbc_south_psi3_docalib));
              bbcs_v3eta_both_docalib->Fill(eta,cosbbc_dphi3_docalib);
              bbcs_v3eta_west_docalib->Fill(eta,cosbbc_dphi3_docalib);
              bbcs_v3eta_east_docalib->Fill(eta,cosbbc_dphi3_docalib);

              // --- autocorrelations in FVTX south...
              double cosfvtx_dphi2_docalib = cos(2*(phi-fvtx_south_psi2_docalib));
              fvtxs_v2eta_both_docalib->Fill(eta,cosfvtx_dphi2_docalib);
              fvtxs_v2eta_west_docalib->Fill(eta,cosfvtx_dphi2_docalib);
              fvtxs_v2eta_east_docalib->Fill(eta,cosfvtx_dphi2_docalib);
              double cosfvtx_dphi3_docalib = cos(3*(phi-fvtx_south_psi3_docalib));
              fvtxs_v3eta_both_docalib->Fill(eta,cosfvtx_dphi3_docalib);
              fvtxs_v3eta_west_docalib->Fill(eta,cosfvtx_dphi3_docalib);
              fvtxs_v3eta_east_docalib->Fill(eta,cosfvtx_dphi3_docalib);

              // --- autocorrelations in FVTX north...
              cosfvtx_dphi2_docalib = cos(2*(phi-fvtx_north_psi2_docalib));
              fvtxn_v2eta_both_docalib->Fill(eta,cosfvtx_dphi2_docalib);
              fvtxn_v2eta_west_docalib->Fill(eta,cosfvtx_dphi2_docalib);
              fvtxn_v2eta_east_docalib->Fill(eta,cosfvtx_dphi2_docalib);
              cosfvtx_dphi3_docalib = cos(3*(phi-fvtx_north_psi3_docalib));
              fvtxn_v3eta_both_docalib->Fill(eta,cosfvtx_dphi3_docalib);
              fvtxn_v3eta_west_docalib->Fill(eta,cosfvtx_dphi3_docalib);
              fvtxn_v3eta_east_docalib->Fill(eta,cosfvtx_dphi3_docalib);

              // --- now to do scalar products
              float ux2 = cos(2*phi);
              float uy2 = sin(2*phi);
              float ux3 = cos(3*phi);
              float uy3 = sin(3*phi);
              // float ux4 = cos(4*phi);
              // float uy4 = sin(4*phi);
              float os_bbc_uq2 = ( (ux2*os_bbc_qx2) + (uy2*os_bbc_qy2) ) / ( os_bbc_qw );
              float os_bbc_uq3 = ( (ux3*os_bbc_qx3) + (uy3*os_bbc_qy3) ) / ( os_bbc_qw );
              float os_fvtxs_uq2 = ( (ux2*os_fvtxs_qx2) + (uy2*os_fvtxs_qy2) - ns ) / ( os_fvtxs_qw - ns );
              float os_fvtxs_uq3 = ( (ux3*os_fvtxs_qx3) + (uy3*os_fvtxs_qy3) - ns ) / ( os_fvtxs_qw - ns );
              float os_fvtxn_uq2 = ( (ux2*os_fvtxn_qx2) + (uy2*os_fvtxn_qy2) - nn ) / ( os_fvtxn_qw - nn );
              float os_fvtxn_uq3 = ( (ux3*os_fvtxn_qx3) + (uy3*os_fvtxn_qy3) - nn ) / ( os_fvtxn_qw - nn );
              float os_fvtxs_tracks_uq2 = ( (ux2*os_fvtxs_tracks_qx2) + (uy2*os_fvtxs_tracks_qy2) - ns ) / ( os_fvtxs_tracks_qw - ns );
              float os_fvtxs_tracks_uq3 = ( (ux3*os_fvtxs_tracks_qx3) + (uy3*os_fvtxs_tracks_qy3) - ns ) / ( os_fvtxs_tracks_qw - ns );
              float os_fvtxn_tracks_uq2 = ( (ux2*os_fvtxn_tracks_qx2) + (uy2*os_fvtxn_tracks_qy2) - nn ) / ( os_fvtxn_tracks_qw - nn );
              float os_fvtxn_tracks_uq3 = ( (ux3*os_fvtxn_tracks_qx3) + (uy3*os_fvtxn_tracks_qy3) - nn ) / ( os_fvtxn_tracks_qw - nn );

              os_bbcs_d22eta_both->Fill(eta,os_bbc_uq2);
              os_bbcs_d22eta_west->Fill(eta,os_bbc_uq2);
              os_bbcs_d22eta_east->Fill(eta,os_bbc_uq2);
              os_fvtxs_d22eta_both->Fill(eta,os_fvtxs_uq2);
              os_fvtxs_d22eta_west->Fill(eta,os_fvtxs_uq2);
              os_fvtxs_d22eta_east->Fill(eta,os_fvtxs_uq2);
              os_fvtxn_d22eta_both->Fill(eta,os_fvtxn_uq2);
              os_fvtxn_d22eta_west->Fill(eta,os_fvtxn_uq2);
              os_fvtxn_d22eta_east->Fill(eta,os_fvtxn_uq2);
              os_fvtxs_tracks_d22eta_both->Fill(eta,os_fvtxs_tracks_uq2);
              os_fvtxs_tracks_d22eta_west->Fill(eta,os_fvtxs_tracks_uq2);
              os_fvtxs_tracks_d22eta_east->Fill(eta,os_fvtxs_tracks_uq2);
              os_fvtxn_tracks_d22eta_both->Fill(eta,os_fvtxn_tracks_uq2);
              os_fvtxn_tracks_d22eta_west->Fill(eta,os_fvtxn_tracks_uq2);
              os_fvtxn_tracks_d22eta_east->Fill(eta,os_fvtxn_tracks_uq2);

              os_bbcs_d32eta_both->Fill(eta,os_bbc_uq3);
              os_bbcs_d32eta_west->Fill(eta,os_bbc_uq3);
              os_bbcs_d32eta_east->Fill(eta,os_bbc_uq3);
              os_fvtxs_d32eta_both->Fill(eta,os_fvtxs_uq3);
              os_fvtxs_d32eta_west->Fill(eta,os_fvtxs_uq3);
              os_fvtxs_d32eta_east->Fill(eta,os_fvtxs_uq3);
              os_fvtxn_d32eta_both->Fill(eta,os_fvtxn_uq3);
              os_fvtxn_d32eta_west->Fill(eta,os_fvtxn_uq3);
              os_fvtxn_d32eta_east->Fill(eta,os_fvtxn_uq3);
              os_fvtxs_tracks_d32eta_both->Fill(eta,os_fvtxs_tracks_uq3);
              os_fvtxs_tracks_d32eta_west->Fill(eta,os_fvtxs_tracks_uq3);
              os_fvtxs_tracks_d32eta_east->Fill(eta,os_fvtxs_tracks_uq3);
              os_fvtxn_tracks_d32eta_both->Fill(eta,os_fvtxn_tracks_uq3);
              os_fvtxn_tracks_d32eta_west->Fill(eta,os_fvtxn_tracks_uq3);
              os_fvtxn_tracks_d32eta_east->Fill(eta,os_fvtxn_tracks_uq3);

            }


        } // check on tracks
    }//end of event

  cout << "Processed " << event_counter << "/" << all_counter << " events (" << (float)event_counter/(float)all_counter << ")" << endl;
  cout << "Events with bad centrality = " << bad_cent_counter << " (" << (float)bad_cent_counter/(float)event_counter << ")" << endl;
  cout << "Events with vertex disagreement = " << bad_vertex_counter << " (" << (float)bad_vertex_counter/(float)event_counter << ")" << endl;
  cout << "Gap cuts applied " << gapcut_counter << "/" << cluster_counter << " (" << (float)gapcut_counter/(float)cluster_counter << ")" << endl;

  if ( rp_recal_pass < 3 && rp_recal_pass > 0 )
    {
      // --- previous pass calib file is named above, rename it here
      sprintf(calibfile,"output/flattening_data/flattening_%d_%d.dat",runNumber,rp_recal_pass);
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
		  for ( int id = 0; id < NDET; id++ )
		    {
		      for ( int ib = 0; ib < 2; ib++ )
			{
			  if ( DIAG ) cout<<"writing ave:  "<<ic<<" "<<iz<<" "<<ih<<" "<<id<<endl;
			  // write out average qx, qx error, qy, qy error
			  ofs << ave[ic][iz][ih][id]->GetBinContent(ib+1) << " ";
			  ofs << ave[ic][iz][ih][id]->GetBinError  (ib+1) << " ";
			} // x and y
		      ofs << endl;
		      for ( int ib = 0; ib < 2; ib++ )
			{
			  for ( int io = 0; io < NORD; io++ )
			    {
			      // write first 12 orders for fourier fit of psi
			      // we are unsure of what's being written out here
			      ofs << flt[ic][iz][ih][id]->GetBinContent(ib*NORD+io+1) << " ";
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

  // cout << "v2 histogram output file: " << outFile1 << endl;
  // cout << "EP histogram output file: " << outFile2 << endl;
  cout << "histogram output file: " << outFile1 << endl;

  if ( rp_recal_pass == 3 )
    {
      mData1->Write();
      //mData2->Write();
      mData1->Close();
      //mData2->Close();
    }

  cout<<"cleaning up"<<endl;



  ntp_event_chain->Delete();

  // ---

  cout << "now attempting to close weight file" << endl;

  if ( phi_weight_file ) phi_weight_file->Close();

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


float calc2_event(float Xn, float Yn, float M)
{

  float numerator = Xn*Xn + Yn*Yn - M;
  float denominator = M*(M-1);

  return numerator/denominator;

}


float calc4_event(float Xn, float Yn, float X2n, float Y2n, float M)
{

  if ( M < 5 ) return -9999;

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
