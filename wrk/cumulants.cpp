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





bool DIAG = false;


int get_fvtx_layer(float);
void dooffsets(int);
void documulants(int);
float calc2_event(float, float, float);
float calc4_event(float, float, float, float, float);
float calc6_event(TComplex&, TComplex&, TComplex&, float);

float calccossum2_event(TComplex&, TComplex&, float);
float calcsinsum2_event(TComplex&, TComplex&, float);
float calccos3_event(TComplex&, TComplex&, float);
float calcsin3_event(TComplex&, TComplex&, float);


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

  dooffsets(run);
  documulants(run);

  return 0;

}

// -----------------------------------------------------------------
void documulants(int runNumber)
{

  cout << "runNumber = " << runNumber << endl;

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


  char filename[500];

  // --- get the number of files for this run number
  string pipe_out = (string) gSystem->GetFromPipe(Form("ls input/%d_*.root | grep -c r",runNumber));
  int nfiles = 0;
  nfiles = atoi(pipe_out.c_str());
  cout<<"nfiles: "<<nfiles<<endl;
  if(nfiles==0) return;

  // --- make a new TChain for the tree
  TChain *ntp_event_chain = new TChain("ntp_event");
  for ( int ifile = 0; ifile < nfiles; ++ifile )
    {
      sprintf(filename,"input/%d_%d.root",runNumber,ifile);
      cout<<"adding to tchain: "<<filename<<endl;
      ntp_event_chain->Add(filename);
    }


  // ---
  bool doweights = true; // change to false if you don't want to bother
  TString phiweightfile_name = Form("SpecialProjects/WeightFiles/newweight2d_run%d.root", runNumber);
  TFile* phi_weight_file = TFile::Open(phiweightfile_name); // COME BACK HERE AND HAVE A LOOK
  if ( !phi_weight_file )
    {
      if ( doweights ) cout << "WARNING could not open phi weight file: " << phiweightfile_name.Data() << endl;
      doweights = false;
    }


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

  if ( doweights ) cout << "All histograms present and ready for inverse phi weighting" << endl;

  // ---

  char outFile2[300];
  sprintf(outFile2,"%s%d%s%d%s","output/files_",energyflag,"/coffsets_",runNumber,".root");
  cout << "histogram for qvector offsets input file: " << outFile2 << endl;
  TFile *mData2=TFile::Open(outFile2);
  mData2->cd();

  TProfile* nfvtxt_tracks_south_qx2 = (TProfile*)mData2->Get("nfvtxt_tracks_south_qx2");
  TProfile* nfvtxt_tracks_south_qx3 = (TProfile*)mData2->Get("nfvtxt_tracks_south_qx3");
  TProfile* nfvtxt_tracks_south_qx4 = (TProfile*)mData2->Get("nfvtxt_tracks_south_qx4");
  TProfile* nfvtxt_tracks_south_qx6 = (TProfile*)mData2->Get("nfvtxt_tracks_south_qx6");
  TProfile* nfvtxt_tracks_south_qy2 = (TProfile*)mData2->Get("nfvtxt_tracks_south_qy2");
  TProfile* nfvtxt_tracks_south_qy3 = (TProfile*)mData2->Get("nfvtxt_tracks_south_qy3");
  TProfile* nfvtxt_tracks_south_qy4 = (TProfile*)mData2->Get("nfvtxt_tracks_south_qy4");
  TProfile* nfvtxt_tracks_south_qy6 = (TProfile*)mData2->Get("nfvtxt_tracks_south_qy6");

  TProfile* nfvtxt_tracks_south_inner_qx2 = (TProfile*)mData2->Get("nfvtxt_tracks_south_inner_qx2");
  TProfile* nfvtxt_tracks_south_inner_qx3 = (TProfile*)mData2->Get("nfvtxt_tracks_south_inner_qx3");
  TProfile* nfvtxt_tracks_south_inner_qx4 = (TProfile*)mData2->Get("nfvtxt_tracks_south_inner_qx4");
  TProfile* nfvtxt_tracks_south_inner_qx6 = (TProfile*)mData2->Get("nfvtxt_tracks_south_inner_qx6");
  TProfile* nfvtxt_tracks_south_inner_qy2 = (TProfile*)mData2->Get("nfvtxt_tracks_south_inner_qy2");
  TProfile* nfvtxt_tracks_south_inner_qy3 = (TProfile*)mData2->Get("nfvtxt_tracks_south_inner_qy3");
  TProfile* nfvtxt_tracks_south_inner_qy4 = (TProfile*)mData2->Get("nfvtxt_tracks_south_inner_qy4");
  TProfile* nfvtxt_tracks_south_inner_qy6 = (TProfile*)mData2->Get("nfvtxt_tracks_south_inner_qy6");

  TProfile* nfvtxt_tracks_south_outer_qx2 = (TProfile*)mData2->Get("nfvtxt_tracks_south_outer_qx2");
  TProfile* nfvtxt_tracks_south_outer_qx3 = (TProfile*)mData2->Get("nfvtxt_tracks_south_outer_qx3");
  TProfile* nfvtxt_tracks_south_outer_qx4 = (TProfile*)mData2->Get("nfvtxt_tracks_south_outer_qx4");
  TProfile* nfvtxt_tracks_south_outer_qx6 = (TProfile*)mData2->Get("nfvtxt_tracks_south_outer_qx6");
  TProfile* nfvtxt_tracks_south_outer_qy2 = (TProfile*)mData2->Get("nfvtxt_tracks_south_outer_qy2");
  TProfile* nfvtxt_tracks_south_outer_qy3 = (TProfile*)mData2->Get("nfvtxt_tracks_south_outer_qy3");
  TProfile* nfvtxt_tracks_south_outer_qy4 = (TProfile*)mData2->Get("nfvtxt_tracks_south_outer_qy4");
  TProfile* nfvtxt_tracks_south_outer_qy6 = (TProfile*)mData2->Get("nfvtxt_tracks_south_outer_qy6");

  TProfile* nfvtxt_tracks_north_qx2 = (TProfile*)mData2->Get("nfvtxt_tracks_north_qx2");
  TProfile* nfvtxt_tracks_north_qx3 = (TProfile*)mData2->Get("nfvtxt_tracks_north_qx3");
  TProfile* nfvtxt_tracks_north_qx4 = (TProfile*)mData2->Get("nfvtxt_tracks_north_qx4");
  TProfile* nfvtxt_tracks_north_qx6 = (TProfile*)mData2->Get("nfvtxt_tracks_north_qx6");
  TProfile* nfvtxt_tracks_north_qy2 = (TProfile*)mData2->Get("nfvtxt_tracks_north_qy2");
  TProfile* nfvtxt_tracks_north_qy3 = (TProfile*)mData2->Get("nfvtxt_tracks_north_qy3");
  TProfile* nfvtxt_tracks_north_qy4 = (TProfile*)mData2->Get("nfvtxt_tracks_north_qy4");
  TProfile* nfvtxt_tracks_north_qy6 = (TProfile*)mData2->Get("nfvtxt_tracks_north_qy6");

  TProfile* nfvtxt_tracks_north_inner_qx2 = (TProfile*)mData2->Get("nfvtxt_tracks_north_inner_qx2");
  TProfile* nfvtxt_tracks_north_inner_qx3 = (TProfile*)mData2->Get("nfvtxt_tracks_north_inner_qx3");
  TProfile* nfvtxt_tracks_north_inner_qx4 = (TProfile*)mData2->Get("nfvtxt_tracks_north_inner_qx4");
  TProfile* nfvtxt_tracks_north_inner_qx6 = (TProfile*)mData2->Get("nfvtxt_tracks_north_inner_qx6");
  TProfile* nfvtxt_tracks_north_inner_qy2 = (TProfile*)mData2->Get("nfvtxt_tracks_north_inner_qy2");
  TProfile* nfvtxt_tracks_north_inner_qy3 = (TProfile*)mData2->Get("nfvtxt_tracks_north_inner_qy3");
  TProfile* nfvtxt_tracks_north_inner_qy4 = (TProfile*)mData2->Get("nfvtxt_tracks_north_inner_qy4");
  TProfile* nfvtxt_tracks_north_inner_qy6 = (TProfile*)mData2->Get("nfvtxt_tracks_north_inner_qy6");

  TProfile* nfvtxt_tracks_north_outer_qx2 = (TProfile*)mData2->Get("nfvtxt_tracks_north_outer_qx2");
  TProfile* nfvtxt_tracks_north_outer_qx3 = (TProfile*)mData2->Get("nfvtxt_tracks_north_outer_qx3");
  TProfile* nfvtxt_tracks_north_outer_qx4 = (TProfile*)mData2->Get("nfvtxt_tracks_north_outer_qx4");
  TProfile* nfvtxt_tracks_north_outer_qx6 = (TProfile*)mData2->Get("nfvtxt_tracks_north_outer_qx6");
  TProfile* nfvtxt_tracks_north_outer_qy2 = (TProfile*)mData2->Get("nfvtxt_tracks_north_outer_qy2");
  TProfile* nfvtxt_tracks_north_outer_qy3 = (TProfile*)mData2->Get("nfvtxt_tracks_north_outer_qy3");
  TProfile* nfvtxt_tracks_north_outer_qy4 = (TProfile*)mData2->Get("nfvtxt_tracks_north_outer_qy4");
  TProfile* nfvtxt_tracks_north_outer_qy6 = (TProfile*)mData2->Get("nfvtxt_tracks_north_outer_qy6");

  if ( !nfvtxt_tracks_south_qx2 ) cout << "Problem finding histogram nfvtxt_tracks_south_qx2 " << endl;
  if ( !nfvtxt_tracks_south_qx3 ) cout << "Problem finding histogram nfvtxt_tracks_south_qx3 " << endl;
  if ( !nfvtxt_tracks_south_qx4 ) cout << "Problem finding histogram nfvtxt_tracks_south_qx4 " << endl;
  if ( !nfvtxt_tracks_south_qx6 ) cout << "Problem finding histogram nfvtxt_tracks_south_qx6 " << endl;
  if ( !nfvtxt_tracks_south_qy2 ) cout << "Problem finding histogram nfvtxt_tracks_south_qy2 " << endl;
  if ( !nfvtxt_tracks_south_qy3 ) cout << "Problem finding histogram nfvtxt_tracks_south_qy3 " << endl;
  if ( !nfvtxt_tracks_south_qy4 ) cout << "Problem finding histogram nfvtxt_tracks_south_qy4 " << endl;
  if ( !nfvtxt_tracks_south_qy6 ) cout << "Problem finding histogram nfvtxt_tracks_south_qy6 " << endl;

  if ( !nfvtxt_tracks_south_inner_qx2 ) cout << "Problem finding histogram nfvtxt_tracks_south_inner_qx2 " << endl;
  if ( !nfvtxt_tracks_south_inner_qx3 ) cout << "Problem finding histogram nfvtxt_tracks_south_inner_qx3 " << endl;
  if ( !nfvtxt_tracks_south_inner_qx4 ) cout << "Problem finding histogram nfvtxt_tracks_south_inner_qx4 " << endl;
  if ( !nfvtxt_tracks_south_inner_qx6 ) cout << "Problem finding histogram nfvtxt_tracks_south_inner_qx6 " << endl;
  if ( !nfvtxt_tracks_south_inner_qy2 ) cout << "Problem finding histogram nfvtxt_tracks_south_inner_qy2 " << endl;
  if ( !nfvtxt_tracks_south_inner_qy3 ) cout << "Problem finding histogram nfvtxt_tracks_south_inner_qy3 " << endl;
  if ( !nfvtxt_tracks_south_inner_qy4 ) cout << "Problem finding histogram nfvtxt_tracks_south_inner_qy4 " << endl;
  if ( !nfvtxt_tracks_south_inner_qy6 ) cout << "Problem finding histogram nfvtxt_tracks_south_inner_qy6 " << endl;

  if ( !nfvtxt_tracks_south_outer_qx2 ) cout << "Problem finding histogram nfvtxt_tracks_south_outer_qx2 " << endl;
  if ( !nfvtxt_tracks_south_outer_qx3 ) cout << "Problem finding histogram nfvtxt_tracks_south_outer_qx3 " << endl;
  if ( !nfvtxt_tracks_south_outer_qx4 ) cout << "Problem finding histogram nfvtxt_tracks_south_outer_qx4 " << endl;
  if ( !nfvtxt_tracks_south_outer_qx6 ) cout << "Problem finding histogram nfvtxt_tracks_south_outer_qx6 " << endl;
  if ( !nfvtxt_tracks_south_outer_qy2 ) cout << "Problem finding histogram nfvtxt_tracks_south_outer_qy2 " << endl;
  if ( !nfvtxt_tracks_south_outer_qy3 ) cout << "Problem finding histogram nfvtxt_tracks_south_outer_qy3 " << endl;
  if ( !nfvtxt_tracks_south_outer_qy4 ) cout << "Problem finding histogram nfvtxt_tracks_south_outer_qy4 " << endl;
  if ( !nfvtxt_tracks_south_outer_qy6 ) cout << "Problem finding histogram nfvtxt_tracks_south_outer_qy6 " << endl;

  if ( !nfvtxt_tracks_north_qx2 ) cout << "Problem finding histogram nfvtxt_tracks_north_qx2 " << endl;
  if ( !nfvtxt_tracks_north_qx3 ) cout << "Problem finding histogram nfvtxt_tracks_north_qx3 " << endl;
  if ( !nfvtxt_tracks_north_qx4 ) cout << "Problem finding histogram nfvtxt_tracks_north_qx4 " << endl;
  if ( !nfvtxt_tracks_north_qx6 ) cout << "Problem finding histogram nfvtxt_tracks_north_qx6 " << endl;
  if ( !nfvtxt_tracks_north_qy2 ) cout << "Problem finding histogram nfvtxt_tracks_north_qy2 " << endl;
  if ( !nfvtxt_tracks_north_qy3 ) cout << "Problem finding histogram nfvtxt_tracks_north_qy3 " << endl;
  if ( !nfvtxt_tracks_north_qy4 ) cout << "Problem finding histogram nfvtxt_tracks_north_qy4 " << endl;
  if ( !nfvtxt_tracks_north_qy6 ) cout << "Problem finding histogram nfvtxt_tracks_north_qy6 " << endl;

  if ( !nfvtxt_tracks_north_inner_qx2 ) cout << "Problem finding histogram nfvtxt_tracks_north_inner_qx2 " << endl;
  if ( !nfvtxt_tracks_north_inner_qx3 ) cout << "Problem finding histogram nfvtxt_tracks_north_inner_qx3 " << endl;
  if ( !nfvtxt_tracks_north_inner_qx4 ) cout << "Problem finding histogram nfvtxt_tracks_north_inner_qx4 " << endl;
  if ( !nfvtxt_tracks_north_inner_qx6 ) cout << "Problem finding histogram nfvtxt_tracks_north_inner_qx6 " << endl;
  if ( !nfvtxt_tracks_north_inner_qy2 ) cout << "Problem finding histogram nfvtxt_tracks_north_inner_qy2 " << endl;
  if ( !nfvtxt_tracks_north_inner_qy3 ) cout << "Problem finding histogram nfvtxt_tracks_north_inner_qy3 " << endl;
  if ( !nfvtxt_tracks_north_inner_qy4 ) cout << "Problem finding histogram nfvtxt_tracks_north_inner_qy4 " << endl;
  if ( !nfvtxt_tracks_north_inner_qy6 ) cout << "Problem finding histogram nfvtxt_tracks_north_inner_qy6 " << endl;

  if ( !nfvtxt_tracks_north_outer_qx2 ) cout << "Problem finding histogram nfvtxt_tracks_north_outer_qx2 " << endl;
  if ( !nfvtxt_tracks_north_outer_qx3 ) cout << "Problem finding histogram nfvtxt_tracks_north_outer_qx3 " << endl;
  if ( !nfvtxt_tracks_north_outer_qx4 ) cout << "Problem finding histogram nfvtxt_tracks_north_outer_qx4 " << endl;
  if ( !nfvtxt_tracks_north_outer_qx6 ) cout << "Problem finding histogram nfvtxt_tracks_north_outer_qx6 " << endl;
  if ( !nfvtxt_tracks_north_outer_qy2 ) cout << "Problem finding histogram nfvtxt_tracks_north_outer_qy2 " << endl;
  if ( !nfvtxt_tracks_north_outer_qy3 ) cout << "Problem finding histogram nfvtxt_tracks_north_outer_qy3 " << endl;
  if ( !nfvtxt_tracks_north_outer_qy4 ) cout << "Problem finding histogram nfvtxt_tracks_north_outer_qy4 " << endl;
  if ( !nfvtxt_tracks_north_outer_qy6 ) cout << "Problem finding histogram nfvtxt_tracks_north_outer_qy6 " << endl;



  char outFile1[300];
  sprintf(outFile1,"%s%d%s%d%s","output/files_",energyflag,"/cumulants_",runNumber,".root");
  cout << "histogram output file: " << outFile1 << endl;
  TFile *mData1=TFile::Open(outFile1,"recreate");
  mData1->cd();

  //------------------------------------------------------------//
  //                  Initializing histograms                   //
  //------------------------------------------------------------//

  cout << "Initializing more histograms" << endl;

  // --- base histograms

  TProfile* nfvtxt_ac_fvtxs_c22 = new TProfile(Form("nfvtxt_ac_fvtxs_c22"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_ac_fvtxn_c22 = new TProfile(Form("nfvtxt_ac_fvtxn_c22"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_ac_fvtxc_c22 = new TProfile(Form("nfvtxt_ac_fvtxc_c22"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_ac_fvtxs_c24 = new TProfile(Form("nfvtxt_ac_fvtxs_c24"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_ac_fvtxn_c24 = new TProfile(Form("nfvtxt_ac_fvtxn_c24"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_ac_fvtxc_c24 = new TProfile(Form("nfvtxt_ac_fvtxc_c24"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_ac_fvtxs_tracks_c22 = new TProfile(Form("nfvtxt_ac_fvtxs_tracks_c22"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_ac_fvtxn_tracks_c22 = new TProfile(Form("nfvtxt_ac_fvtxn_tracks_c22"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_ac_fvtxc_tracks_c22 = new TProfile(Form("nfvtxt_ac_fvtxc_tracks_c22"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_ac_fvtxs_tracks_c24 = new TProfile(Form("nfvtxt_ac_fvtxs_tracks_c24"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_ac_fvtxn_tracks_c24 = new TProfile(Form("nfvtxt_ac_fvtxn_tracks_c24"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_ac_fvtxc_tracks_c24 = new TProfile(Form("nfvtxt_ac_fvtxc_tracks_c24"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_ac_fvtxs_tracks_c26 = new TProfile(Form("nfvtxt_ac_fvtxs_tracks_c26"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_ac_fvtxn_tracks_c26 = new TProfile(Form("nfvtxt_ac_fvtxn_tracks_c26"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_ac_fvtxc_tracks_c26 = new TProfile(Form("nfvtxt_ac_fvtxc_tracks_c26"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_ac_fvtxs_tracks_cov24 = new TProfile(Form("nfvtxt_ac_fvtxs_tracks_cov24"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_ac_fvtxn_tracks_cov24 = new TProfile(Form("nfvtxt_ac_fvtxn_tracks_cov24"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_ac_fvtxc_tracks_cov24 = new TProfile(Form("nfvtxt_ac_fvtxc_tracks_cov24"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_ac_fvtxsfvtxn_c22  = new TProfile(Form("nfvtxt_ac_fvtxsfvtxn_c22"), "",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_ac_fvtxsfvtxn_c24  = new TProfile(Form("nfvtxt_ac_fvtxsfvtxn_c24"), "",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_ac_fvtxsfvtxn_c24a = new TProfile(Form("nfvtxt_ac_fvtxsfvtxn_c24a"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_ac_fvtxsfvtxn_c24b = new TProfile(Form("nfvtxt_ac_fvtxsfvtxn_c24b"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_ac_fvtxsfvtxn_c24c = new TProfile(Form("nfvtxt_ac_fvtxsfvtxn_c24c"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_ac_fvtxsfvtxn_c24d = new TProfile(Form("nfvtxt_ac_fvtxsfvtxn_c24d"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_ac_fvtxsfvtxn_tracks_c22  = new TProfile(Form("nfvtxt_ac_fvtxsfvtxn_tracks_c22"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_ac_fvtxsfvtxn_tracks_c24  = new TProfile(Form("nfvtxt_ac_fvtxsfvtxn_tracks_c24"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_ac_fvtxsfvtxn_tracks_c24a = new TProfile(Form("nfvtxt_ac_fvtxsfvtxn_tracks_c24a"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_ac_fvtxsfvtxn_tracks_c24b = new TProfile(Form("nfvtxt_ac_fvtxsfvtxn_tracks_c24b"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_ac_fvtxsfvtxn_tracks_c24c = new TProfile(Form("nfvtxt_ac_fvtxsfvtxn_tracks_c24c"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_ac_fvtxsfvtxn_tracks_c24d = new TProfile(Form("nfvtxt_ac_fvtxsfvtxn_tracks_c24d"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_ac_fvtxs_ce01_c22 = new TProfile(Form("nfvtxt_ac_fvtxs_ce01_c22"),       "",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_ac_fvtxn_ce01_c22 = new TProfile(Form("nfvtxt_ac_fvtxn_ce01_c22"),       "",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_ac_fvtxsfvtxn_ce01_c24  = new TProfile(Form("nfvtxt_ac_fvtxsfvtxn_ce01_c24"), "",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_ac_fvtxsfvtxn_ce01_c24a = new TProfile(Form("nfvtxt_ac_fvtxsfvtxn_ce01_c24a"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_ac_fvtxsfvtxn_ce01_c24b = new TProfile(Form("nfvtxt_ac_fvtxsfvtxn_ce01_c24b"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_ac_fvtxsfvtxn_ce01_c24c = new TProfile(Form("nfvtxt_ac_fvtxsfvtxn_ce01_c24c"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_ac_fvtxsfvtxn_ce01_c24d = new TProfile(Form("nfvtxt_ac_fvtxsfvtxn_ce01_c24d"),"",80, -0.5, 79.5, -1.1, 1.1);

  // --- correction histograms

  // --- <<cos(n(phi1))>>
  TProfile* nfvtxt_ac_fvtxs_cos21 = new TProfile(Form("nfvtxt_ac_fvtxs_cos21"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_ac_fvtxn_cos21 = new TProfile(Form("nfvtxt_ac_fvtxn_cos21"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_ac_fvtxc_cos21 = new TProfile(Form("nfvtxt_ac_fvtxc_cos21"),"",80, -0.5, 79.5, -1.1, 1.1);
  // --- <<sin(n(phi1))>>
  TProfile* nfvtxt_ac_fvtxs_sin21 = new TProfile(Form("nfvtxt_ac_fvtxs_sin21"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_ac_fvtxn_sin21 = new TProfile(Form("nfvtxt_ac_fvtxn_sin21"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_ac_fvtxc_sin21 = new TProfile(Form("nfvtxt_ac_fvtxc_sin21"),"",80, -0.5, 79.5, -1.1, 1.1);
  // --- <<cos(n(phi1+phi2))>>
  TProfile* nfvtxt_ac_fvtxs_cossum22 = new TProfile(Form("nfvtxt_ac_fvtxs_cossum22"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_ac_fvtxn_cossum22 = new TProfile(Form("nfvtxt_ac_fvtxn_cossum22"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_ac_fvtxc_cossum22 = new TProfile(Form("nfvtxt_ac_fvtxc_cossum22"),"",80, -0.5, 79.5, -1.1, 1.1);
  // --- <<sin(n(phi1+phi2))>>
  TProfile* nfvtxt_ac_fvtxs_sinsum22 = new TProfile(Form("nfvtxt_ac_fvtxs_sinsum22"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_ac_fvtxn_sinsum22 = new TProfile(Form("nfvtxt_ac_fvtxn_sinsum22"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_ac_fvtxc_sinsum22 = new TProfile(Form("nfvtxt_ac_fvtxc_sinsum22"),"",80, -0.5, 79.5, -1.1, 1.1);
  // --- <<cos(n(phi1-phi2-phi3))>>
  TProfile* nfvtxt_ac_fvtxs_cos23 = new TProfile(Form("nfvtxt_ac_fvtxs_cos23"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_ac_fvtxn_cos23 = new TProfile(Form("nfvtxt_ac_fvtxn_cos23"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_ac_fvtxc_cos23 = new TProfile(Form("nfvtxt_ac_fvtxc_cos23"),"",80, -0.5, 79.5, -1.1, 1.1);
  // --- <<sin(n(phi1-phi2-phi3))>>
  TProfile* nfvtxt_ac_fvtxs_sin23 = new TProfile(Form("nfvtxt_ac_fvtxs_sin23"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_ac_fvtxn_sin23 = new TProfile(Form("nfvtxt_ac_fvtxn_sin23"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_ac_fvtxc_sin23 = new TProfile(Form("nfvtxt_ac_fvtxc_sin23"),"",80, -0.5, 79.5, -1.1, 1.1);

  // --- <<cos(n(phi1))>>
  TProfile* nfvtxt_ac_fvtxs_tracks_cos21 = new TProfile(Form("nfvtxt_ac_fvtxs_tracks_cos21"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_ac_fvtxn_tracks_cos21 = new TProfile(Form("nfvtxt_ac_fvtxn_tracks_cos21"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_ac_fvtxc_tracks_cos21 = new TProfile(Form("nfvtxt_ac_fvtxc_tracks_cos21"),"",80, -0.5, 79.5, -1.1, 1.1);
  // --- <<sin(n(phi1))>>
  TProfile* nfvtxt_ac_fvtxs_tracks_sin21 = new TProfile(Form("nfvtxt_ac_fvtxs_tracks_sin21"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_ac_fvtxn_tracks_sin21 = new TProfile(Form("nfvtxt_ac_fvtxn_tracks_sin21"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_ac_fvtxc_tracks_sin21 = new TProfile(Form("nfvtxt_ac_fvtxc_tracks_sin21"),"",80, -0.5, 79.5, -1.1, 1.1);
  // --- <<cos(n(phi1+phi2))>>
  TProfile* nfvtxt_ac_fvtxs_tracks_cossum22 = new TProfile(Form("nfvtxt_ac_fvtxs_tracks_cossum22"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_ac_fvtxn_tracks_cossum22 = new TProfile(Form("nfvtxt_ac_fvtxn_tracks_cossum22"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_ac_fvtxc_tracks_cossum22 = new TProfile(Form("nfvtxt_ac_fvtxc_tracks_cossum22"),"",80, -0.5, 79.5, -1.1, 1.1);
  // --- <<sin(n(phi1+phi2))>>
  TProfile* nfvtxt_ac_fvtxs_tracks_sinsum22 = new TProfile(Form("nfvtxt_ac_fvtxs_tracks_sinsum22"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_ac_fvtxn_tracks_sinsum22 = new TProfile(Form("nfvtxt_ac_fvtxn_tracks_sinsum22"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_ac_fvtxc_tracks_sinsum22 = new TProfile(Form("nfvtxt_ac_fvtxc_tracks_sinsum22"),"",80, -0.5, 79.5, -1.1, 1.1);
  // --- <<cos(n(phi1-phi2-phi3))>>
  TProfile* nfvtxt_ac_fvtxs_tracks_cos23 = new TProfile(Form("nfvtxt_ac_fvtxs_tracks_cos23"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_ac_fvtxn_tracks_cos23 = new TProfile(Form("nfvtxt_ac_fvtxn_tracks_cos23"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_ac_fvtxc_tracks_cos23 = new TProfile(Form("nfvtxt_ac_fvtxc_tracks_cos23"),"",80, -0.5, 79.5, -1.1, 1.1);
  // --- <<sin(n(phi1-phi2-phi3))>>
  TProfile* nfvtxt_ac_fvtxs_tracks_sin23 = new TProfile(Form("nfvtxt_ac_fvtxs_tracks_sin23"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_ac_fvtxn_tracks_sin23 = new TProfile(Form("nfvtxt_ac_fvtxn_tracks_sin23"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_ac_fvtxc_tracks_sin23 = new TProfile(Form("nfvtxt_ac_fvtxc_tracks_sin23"),"",80, -0.5, 79.5, -1.1, 1.1);

  // --- special stuff...
  // --- main correlations
  TProfile* nfvtxtsp_ac_fvtxs_tracks_c22 = new TProfile(Form("nfvtxtsp_ac_fvtxs_tracks_c22"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxtsp_ac_fvtxn_tracks_c22 = new TProfile(Form("nfvtxtsp_ac_fvtxn_tracks_c22"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxtsp_ac_fvtxc_tracks_c22 = new TProfile(Form("nfvtxtsp_ac_fvtxc_tracks_c22"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxtsp_ac_fvtxs_tracks_c24 = new TProfile(Form("nfvtxtsp_ac_fvtxs_tracks_c24"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxtsp_ac_fvtxn_tracks_c24 = new TProfile(Form("nfvtxtsp_ac_fvtxn_tracks_c24"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxtsp_ac_fvtxc_tracks_c24 = new TProfile(Form("nfvtxtsp_ac_fvtxc_tracks_c24"),"",80, -0.5, 79.5, -1.1, 1.1);
  // --- <<cos(n(phi1))>>
  TProfile* nfvtxtsp_ac_fvtxs_tracks_cos21 = new TProfile(Form("nfvtxtsp_ac_fvtxs_tracks_cos21"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxtsp_ac_fvtxn_tracks_cos21 = new TProfile(Form("nfvtxtsp_ac_fvtxn_tracks_cos21"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxtsp_ac_fvtxc_tracks_cos21 = new TProfile(Form("nfvtxtsp_ac_fvtxc_tracks_cos21"),"",80, -0.5, 79.5, -1.1, 1.1);
  // --- <<sin(n(phi1))>>
  TProfile* nfvtxtsp_ac_fvtxs_tracks_sin21 = new TProfile(Form("nfvtxtsp_ac_fvtxs_tracks_sin21"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxtsp_ac_fvtxn_tracks_sin21 = new TProfile(Form("nfvtxtsp_ac_fvtxn_tracks_sin21"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxtsp_ac_fvtxc_tracks_sin21 = new TProfile(Form("nfvtxtsp_ac_fvtxc_tracks_sin21"),"",80, -0.5, 79.5, -1.1, 1.1);
  // --- <<cos(n(phi1+phi2))>>
  TProfile* nfvtxtsp_ac_fvtxs_tracks_cossum22 = new TProfile(Form("nfvtxtsp_ac_fvtxs_tracks_cossum22"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxtsp_ac_fvtxn_tracks_cossum22 = new TProfile(Form("nfvtxtsp_ac_fvtxn_tracks_cossum22"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxtsp_ac_fvtxc_tracks_cossum22 = new TProfile(Form("nfvtxtsp_ac_fvtxc_tracks_cossum22"),"",80, -0.5, 79.5, -1.1, 1.1);
  // --- <<sin(n(phi1+phi2))>>
  TProfile* nfvtxtsp_ac_fvtxs_tracks_sinsum22 = new TProfile(Form("nfvtxtsp_ac_fvtxs_tracks_sinsum22"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxtsp_ac_fvtxn_tracks_sinsum22 = new TProfile(Form("nfvtxtsp_ac_fvtxn_tracks_sinsum22"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxtsp_ac_fvtxc_tracks_sinsum22 = new TProfile(Form("nfvtxtsp_ac_fvtxc_tracks_sinsum22"),"",80, -0.5, 79.5, -1.1, 1.1);
  // --- <<cos(n(phi1-phi2-phi3))>>
  TProfile* nfvtxtsp_ac_fvtxs_tracks_cos23 = new TProfile(Form("nfvtxtsp_ac_fvtxs_tracks_cos23"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxtsp_ac_fvtxn_tracks_cos23 = new TProfile(Form("nfvtxtsp_ac_fvtxn_tracks_cos23"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxtsp_ac_fvtxc_tracks_cos23 = new TProfile(Form("nfvtxtsp_ac_fvtxc_tracks_cos23"),"",80, -0.5, 79.5, -1.1, 1.1);
  // --- <<sin(n(phi1-phi2-phi3))>>
  TProfile* nfvtxtsp_ac_fvtxs_tracks_sin23 = new TProfile(Form("nfvtxtsp_ac_fvtxs_tracks_sin23"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxtsp_ac_fvtxn_tracks_sin23 = new TProfile(Form("nfvtxtsp_ac_fvtxn_tracks_sin23"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxtsp_ac_fvtxc_tracks_sin23 = new TProfile(Form("nfvtxtsp_ac_fvtxc_tracks_sin23"),"",80, -0.5, 79.5, -1.1, 1.1);

  // ---------------------------------------------------------------------------------------------------------



  // ---------------------------------------------------------------------------------------------------------

  // --- os used to mean offset but now i have to come up with something else

  // TProfile* nfvtxt_zzyzx_fvtxs_c22 = new TProfile(Form("nfvtxt_zzyzx_fvtxs_c22"),"",80, -0.5, 79.5, -1.1, 1.1);
  // TProfile* nfvtxt_zzyzx_fvtxn_c22 = new TProfile(Form("nfvtxt_zzyzx_fvtxn_c22"),"",80, -0.5, 79.5, -1.1, 1.1);
  // TProfile* nfvtxt_zzyzx_fvtxc_c22 = new TProfile(Form("nfvtxt_zzyzx_fvtxc_c22"),"",80, -0.5, 79.5, -1.1, 1.1);
  // TProfile* nfvtxt_zzyzx_fvtxs_c24 = new TProfile(Form("nfvtxt_zzyzx_fvtxs_c24"),"",80, -0.5, 79.5, -1.1, 1.1);
  // TProfile* nfvtxt_zzyzx_fvtxn_c24 = new TProfile(Form("nfvtxt_zzyzx_fvtxn_c24"),"",80, -0.5, 79.5, -1.1, 1.1);
  // TProfile* nfvtxt_zzyzx_fvtxc_c24 = new TProfile(Form("nfvtxt_zzyzx_fvtxc_c24"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_zzyzx_fvtxs_tracks_c22 = new TProfile(Form("nfvtxt_zzyzx_fvtxs_tracks_c22"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_zzyzx_fvtxn_tracks_c22 = new TProfile(Form("nfvtxt_zzyzx_fvtxn_tracks_c22"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_zzyzx_fvtxc_tracks_c22 = new TProfile(Form("nfvtxt_zzyzx_fvtxc_tracks_c22"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_zzyzx_fvtxs_tracks_c24 = new TProfile(Form("nfvtxt_zzyzx_fvtxs_tracks_c24"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_zzyzx_fvtxn_tracks_c24 = new TProfile(Form("nfvtxt_zzyzx_fvtxn_tracks_c24"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_zzyzx_fvtxc_tracks_c24 = new TProfile(Form("nfvtxt_zzyzx_fvtxc_tracks_c24"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_zzyzx_fvtxs_tracks_c26 = new TProfile(Form("nfvtxt_zzyzx_fvtxs_tracks_c26"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_zzyzx_fvtxn_tracks_c26 = new TProfile(Form("nfvtxt_zzyzx_fvtxn_tracks_c26"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_zzyzx_fvtxc_tracks_c26 = new TProfile(Form("nfvtxt_zzyzx_fvtxc_tracks_c26"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_zzyzx_fvtxs_tracks_cov24 = new TProfile(Form("nfvtxt_zzyzx_fvtxs_tracks_cov24"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_zzyzx_fvtxn_tracks_cov24 = new TProfile(Form("nfvtxt_zzyzx_fvtxn_tracks_cov24"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_zzyzx_fvtxc_tracks_cov24 = new TProfile(Form("nfvtxt_zzyzx_fvtxc_tracks_cov24"),"",80, -0.5, 79.5, -1.1, 1.1);
  // TProfile* nfvtxt_zzyzx_fvtxsfvtxn_c22  = new TProfile(Form("nfvtxt_zzyzx_fvtxsfvtxn_c22"), "",80, -0.5, 79.5, -1.1, 1.1);
  // TProfile* nfvtxt_zzyzx_fvtxsfvtxn_c24  = new TProfile(Form("nfvtxt_zzyzx_fvtxsfvtxn_c24"), "",80, -0.5, 79.5, -1.1, 1.1);
  // TProfile* nfvtxt_zzyzx_fvtxsfvtxn_c24a = new TProfile(Form("nfvtxt_zzyzx_fvtxsfvtxn_c24a"),"",80, -0.5, 79.5, -1.1, 1.1);
  // TProfile* nfvtxt_zzyzx_fvtxsfvtxn_c24b = new TProfile(Form("nfvtxt_zzyzx_fvtxsfvtxn_c24b"),"",80, -0.5, 79.5, -1.1, 1.1);
  // TProfile* nfvtxt_zzyzx_fvtxsfvtxn_c24c = new TProfile(Form("nfvtxt_zzyzx_fvtxsfvtxn_c24c"),"",80, -0.5, 79.5, -1.1, 1.1);
  // TProfile* nfvtxt_zzyzx_fvtxsfvtxn_c24d = new TProfile(Form("nfvtxt_zzyzx_fvtxsfvtxn_c24d"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_zzyzx_fvtxsfvtxn_tracks_c22  = new TProfile(Form("nfvtxt_zzyzx_fvtxsfvtxn_tracks_c22"),"",80, -0.5, 79.5, -1.1, 1.1);
  // TProfile* nfvtxt_zzyzx_fvtxsfvtxn_tracks_c24  = new TProfile(Form("nfvtxt_zzyzx_fvtxsfvtxn_tracks_c24"),"",80, -0.5, 79.5, -1.1, 1.1);
  // TProfile* nfvtxt_zzyzx_fvtxsfvtxn_tracks_c24a = new TProfile(Form("nfvtxt_zzyzx_fvtxsfvtxn_tracks_c24a"),"",80, -0.5, 79.5, -1.1, 1.1);
  // TProfile* nfvtxt_zzyzx_fvtxsfvtxn_tracks_c24b = new TProfile(Form("nfvtxt_zzyzx_fvtxsfvtxn_tracks_c24b"),"",80, -0.5, 79.5, -1.1, 1.1);
  // TProfile* nfvtxt_zzyzx_fvtxsfvtxn_tracks_c24c = new TProfile(Form("nfvtxt_zzyzx_fvtxsfvtxn_tracks_c24c"),"",80, -0.5, 79.5, -1.1, 1.1);
  // TProfile* nfvtxt_zzyzx_fvtxsfvtxn_tracks_c24d = new TProfile(Form("nfvtxt_zzyzx_fvtxsfvtxn_tracks_c24d"),"",80, -0.5, 79.5, -1.1, 1.1);
  // TProfile* nfvtxt_zzyzx_fvtxs_ce01_c22 = new TProfile(Form("nfvtxt_zzyzx_fvtxs_ce01_c22"),       "",80, -0.5, 79.5, -1.1, 1.1);
  // TProfile* nfvtxt_zzyzx_fvtxn_ce01_c22 = new TProfile(Form("nfvtxt_zzyzx_fvtxn_ce01_c22"),       "",80, -0.5, 79.5, -1.1, 1.1);
  // TProfile* nfvtxt_zzyzx_fvtxsfvtxn_ce01_c24  = new TProfile(Form("nfvtxt_zzyzx_fvtxsfvtxn_ce01_c24"), "",80, -0.5, 79.5, -1.1, 1.1);
  // TProfile* nfvtxt_zzyzx_fvtxsfvtxn_ce01_c24a = new TProfile(Form("nfvtxt_zzyzx_fvtxsfvtxn_ce01_c24a"),"",80, -0.5, 79.5, -1.1, 1.1);
  // TProfile* nfvtxt_zzyzx_fvtxsfvtxn_ce01_c24b = new TProfile(Form("nfvtxt_zzyzx_fvtxsfvtxn_ce01_c24b"),"",80, -0.5, 79.5, -1.1, 1.1);
  // TProfile* nfvtxt_zzyzx_fvtxsfvtxn_ce01_c24c = new TProfile(Form("nfvtxt_zzyzx_fvtxsfvtxn_ce01_c24c"),"",80, -0.5, 79.5, -1.1, 1.1);
  // TProfile* nfvtxt_zzyzx_fvtxsfvtxn_ce01_c24d = new TProfile(Form("nfvtxt_zzyzx_fvtxsfvtxn_ce01_c24d"),"",80, -0.5, 79.5, -1.1, 1.1);

  // --- correction histograms

  // // --- <<cos(n(phi1))>>
  // TProfile* nfvtxt_zzyzx_fvtxs_cos21 = new TProfile(Form("nfvtxt_zzyzx_fvtxs_cos21"),"",80, -0.5, 79.5, -1.1, 1.1);
  // TProfile* nfvtxt_zzyzx_fvtxn_cos21 = new TProfile(Form("nfvtxt_zzyzx_fvtxn_cos21"),"",80, -0.5, 79.5, -1.1, 1.1);
  // TProfile* nfvtxt_zzyzx_fvtxc_cos21 = new TProfile(Form("nfvtxt_zzyzx_fvtxc_cos21"),"",80, -0.5, 79.5, -1.1, 1.1);
  // // --- <<sin(n(phi1))>>
  // TProfile* nfvtxt_zzyzx_fvtxs_sin21 = new TProfile(Form("nfvtxt_zzyzx_fvtxs_sin21"),"",80, -0.5, 79.5, -1.1, 1.1);
  // TProfile* nfvtxt_zzyzx_fvtxn_sin21 = new TProfile(Form("nfvtxt_zzyzx_fvtxn_sin21"),"",80, -0.5, 79.5, -1.1, 1.1);
  // TProfile* nfvtxt_zzyzx_fvtxc_sin21 = new TProfile(Form("nfvtxt_zzyzx_fvtxc_sin21"),"",80, -0.5, 79.5, -1.1, 1.1);
  // // --- <<cos(n(phi1+phi2))>>
  // TProfile* nfvtxt_zzyzx_fvtxs_cossum22 = new TProfile(Form("nfvtxt_zzyzx_fvtxs_cossum22"),"",80, -0.5, 79.5, -1.1, 1.1);
  // TProfile* nfvtxt_zzyzx_fvtxn_cossum22 = new TProfile(Form("nfvtxt_zzyzx_fvtxn_cossum22"),"",80, -0.5, 79.5, -1.1, 1.1);
  // TProfile* nfvtxt_zzyzx_fvtxc_cossum22 = new TProfile(Form("nfvtxt_zzyzx_fvtxc_cossum22"),"",80, -0.5, 79.5, -1.1, 1.1);
  // // --- <<sin(n(phi1+phi2))>>
  // TProfile* nfvtxt_zzyzx_fvtxs_sinsum22 = new TProfile(Form("nfvtxt_zzyzx_fvtxs_sinsum22"),"",80, -0.5, 79.5, -1.1, 1.1);
  // TProfile* nfvtxt_zzyzx_fvtxn_sinsum22 = new TProfile(Form("nfvtxt_zzyzx_fvtxn_sinsum22"),"",80, -0.5, 79.5, -1.1, 1.1);
  // TProfile* nfvtxt_zzyzx_fvtxc_sinsum22 = new TProfile(Form("nfvtxt_zzyzx_fvtxc_sinsum22"),"",80, -0.5, 79.5, -1.1, 1.1);
  // // --- <<cos(n(phi1-phi2-phi3))>>
  // TProfile* nfvtxt_zzyzx_fvtxs_cos23 = new TProfile(Form("nfvtxt_zzyzx_fvtxs_cos23"),"",80, -0.5, 79.5, -1.1, 1.1);
  // TProfile* nfvtxt_zzyzx_fvtxn_cos23 = new TProfile(Form("nfvtxt_zzyzx_fvtxn_cos23"),"",80, -0.5, 79.5, -1.1, 1.1);
  // TProfile* nfvtxt_zzyzx_fvtxc_cos23 = new TProfile(Form("nfvtxt_zzyzx_fvtxc_cos23"),"",80, -0.5, 79.5, -1.1, 1.1);
  // // --- <<sin(n(phi1-phi2-phi3))>>
  // TProfile* nfvtxt_zzyzx_fvtxs_sin23 = new TProfile(Form("nfvtxt_zzyzx_fvtxs_sin23"),"",80, -0.5, 79.5, -1.1, 1.1);
  // TProfile* nfvtxt_zzyzx_fvtxn_sin23 = new TProfile(Form("nfvtxt_zzyzx_fvtxn_sin23"),"",80, -0.5, 79.5, -1.1, 1.1);
  // TProfile* nfvtxt_zzyzx_fvtxc_sin23 = new TProfile(Form("nfvtxt_zzyzx_fvtxc_sin23"),"",80, -0.5, 79.5, -1.1, 1.1);

  // --- <<cos(n(phi1))>>
  TProfile* nfvtxt_zzyzx_fvtxs_tracks_cos21 = new TProfile(Form("nfvtxt_zzyzx_fvtxs_tracks_cos21"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_zzyzx_fvtxn_tracks_cos21 = new TProfile(Form("nfvtxt_zzyzx_fvtxn_tracks_cos21"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_zzyzx_fvtxc_tracks_cos21 = new TProfile(Form("nfvtxt_zzyzx_fvtxc_tracks_cos21"),"",80, -0.5, 79.5, -1.1, 1.1);
  // --- <<sin(n(phi1))>>
  TProfile* nfvtxt_zzyzx_fvtxs_tracks_sin21 = new TProfile(Form("nfvtxt_zzyzx_fvtxs_tracks_sin21"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_zzyzx_fvtxn_tracks_sin21 = new TProfile(Form("nfvtxt_zzyzx_fvtxn_tracks_sin21"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_zzyzx_fvtxc_tracks_sin21 = new TProfile(Form("nfvtxt_zzyzx_fvtxc_tracks_sin21"),"",80, -0.5, 79.5, -1.1, 1.1);
  // --- <<cos(n(phi1+phi2))>>
  TProfile* nfvtxt_zzyzx_fvtxs_tracks_cossum22 = new TProfile(Form("nfvtxt_zzyzx_fvtxs_tracks_cossum22"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_zzyzx_fvtxn_tracks_cossum22 = new TProfile(Form("nfvtxt_zzyzx_fvtxn_tracks_cossum22"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_zzyzx_fvtxc_tracks_cossum22 = new TProfile(Form("nfvtxt_zzyzx_fvtxc_tracks_cossum22"),"",80, -0.5, 79.5, -1.1, 1.1);
  // --- <<sin(n(phi1+phi2))>>
  TProfile* nfvtxt_zzyzx_fvtxs_tracks_sinsum22 = new TProfile(Form("nfvtxt_zzyzx_fvtxs_tracks_sinsum22"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_zzyzx_fvtxn_tracks_sinsum22 = new TProfile(Form("nfvtxt_zzyzx_fvtxn_tracks_sinsum22"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_zzyzx_fvtxc_tracks_sinsum22 = new TProfile(Form("nfvtxt_zzyzx_fvtxc_tracks_sinsum22"),"",80, -0.5, 79.5, -1.1, 1.1);
  // --- <<cos(n(phi1-phi2-phi3))>>
  TProfile* nfvtxt_zzyzx_fvtxs_tracks_cos23 = new TProfile(Form("nfvtxt_zzyzx_fvtxs_tracks_cos23"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_zzyzx_fvtxn_tracks_cos23 = new TProfile(Form("nfvtxt_zzyzx_fvtxn_tracks_cos23"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_zzyzx_fvtxc_tracks_cos23 = new TProfile(Form("nfvtxt_zzyzx_fvtxc_tracks_cos23"),"",80, -0.5, 79.5, -1.1, 1.1);
  // --- <<sin(n(phi1-phi2-phi3))>>
  TProfile* nfvtxt_zzyzx_fvtxs_tracks_sin23 = new TProfile(Form("nfvtxt_zzyzx_fvtxs_tracks_sin23"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_zzyzx_fvtxn_tracks_sin23 = new TProfile(Form("nfvtxt_zzyzx_fvtxn_tracks_sin23"),"",80, -0.5, 79.5, -1.1, 1.1);
  TProfile* nfvtxt_zzyzx_fvtxc_tracks_sin23 = new TProfile(Form("nfvtxt_zzyzx_fvtxc_tracks_sin23"),"",80, -0.5, 79.5, -1.1, 1.1);

  // ---------------------------------------------------------------------------------------------------------
  TProfile* tp1f_special_fvtx_tracks_ab[8];
  for ( int i = 0; i < 8; ++i ) tp1f_special_fvtx_tracks_ab[i] = new TProfile(Form("tp1f_special_fvtx_tracks_ab%d",i),"",12,-3,3,-1.1,1.1,"");
  TProfile* tp1f_special_fvtx_tracks_aa = new TProfile("tp1f_special_fvtx_tracks_aa","",12,-3,3,-1.1,1.1,"");
  TProfile* tp1f_special_fvtx_tracks_aa_cos = new TProfile("tp1f_special_fvtx_tracks_aa_cos","",12,-3,3,-1.1,1.1,"");
  TProfile* tp1f_special_fvtx_tracks_aa_sin = new TProfile("tp1f_special_fvtx_tracks_aa_sin","",12,-3,3,-1.1,1.1,"");

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
  ntp_event_chain->SetBranchAddress("fnhits",fnhits,&b_fnhits);
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


  int all_counter = 0;

  Long64_t cluster_counter = 0;
  Long64_t gapcut_counter = 0;

  cout << "starting loop over events in the tree" << endl;
  int nentries = ntp_event_chain->GetEntries();
  cout << "total events = " << nentries << endl;
  for ( int ievt = 0; ievt < nentries; ++ievt )
    {

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

      if ( centrality > -1 )
        {
          // if ( energyflag == 200 && centrality > 5  ) continue;
          // if ( energyflag == 62  && centrality > 10 ) continue;
          // if ( energyflag == 20  && centrality > 20 ) continue;
          // if ( energyflag == 39  && centrality > 20 ) continue;
          if ( energyflag == 200 && centrality > 5  ) continue;
          if ( energyflag == 62  && centrality > 5  ) continue;
          if ( energyflag == 20  && centrality > 20 ) continue;
          if ( energyflag == 39  && centrality > 10 ) continue;
        }

      double ZVTX = -9999;
      if ( runNumber >= 454774 && runNumber <= 456283 ) ZVTX = d_bbcz;
      if ( runNumber >= 456652 && runNumber <= 458167 ) ZVTX = eventfvtx_z;
      if ( fabs(ZVTX) > 10.0 )
      //if ( fabs(ZVTX) > 5.0 )
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

      // int nfvtxc = d_nFVTX_clus;
      // int nfvtxc_south = d_nFVTXN_clus;
      // int nfvtxc_north = d_nFVTXS_clus;

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
	      if ( !th1d_fvtxs_phi_weight[fvtx_layer+1] )
		{
		  cout << "WARNING!!!  Problem with weight histograms in cluster loop..." << endl;
		  continue;
		}
	      int phi_bin = th1d_fvtxs_phi_weight[fvtx_layer+1]->FindBin(phi); // COME BACK HERE AND HAVE A LOOK
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

	    } // check on south

	      // ------------------------------------------------
	      // --- determine the weights for north clusters ---
	  if ( doweights )
	    {
	      if ( !th1d_fvtxn_phi_weight[fvtx_layer+1] )
		{
		  cout << "WARNING!!!  Problem with weight histograms in cluster loop..." << endl;
		  continue;
		}
	      int phi_bin = th1d_fvtxn_phi_weight[fvtx_layer+1]->FindBin(phi); // COME BACK HERE AND HAVE A LOOK
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

	    } // check on north

	} // loop over cluster


      //float fvtxs_qq2 = ( (fvtxs_qx2[0]*fvtxs_qx2[0]) + (fvtxs_qy2[0]*fvtxs_qy2[0]) - fvtxs_qw[0] ) / ( (fvtxs_qw[0]*fvtxs_qw[0]) - fvtxs_qw[0] );
      //float fvtxs_qq3 = ( (fvtxs_qx3[0]*fvtxs_qx3[0]) + (fvtxs_qy3[0]*fvtxs_qy3[0]) - fvtxs_qw[0] ) / ( (fvtxs_qw[0]*fvtxs_qw[0]) - fvtxs_qw[0] );

      bool good_4_event = ( nclus_south_inner > 3 ) && ( nclus_south_outer > 3 ) && ( nclus_north_inner > 3 ) && ( nclus_north_outer > 3 ) ;



      // --- fvtx tracks
      float fvtxs_tracks_qx2[3]; // both, inner, outer
      float fvtxs_tracks_qy2[3];
      float fvtxs_tracks_qx3[3];
      float fvtxs_tracks_qy3[3];
      float fvtxs_tracks_qx4[3];
      float fvtxs_tracks_qy4[3];
      float fvtxs_tracks_qx6[3];
      float fvtxs_tracks_qy6[3];
      float fvtxs_tracks_qw[3];
      float fvtxn_tracks_qx2[3]; // both, inner, outer
      float fvtxn_tracks_qy2[3];
      float fvtxn_tracks_qx3[3];
      float fvtxn_tracks_qy3[3];
      float fvtxn_tracks_qx4[3];
      float fvtxn_tracks_qy4[3];
      float fvtxn_tracks_qx6[3];
      float fvtxn_tracks_qy6[3];
      float fvtxn_tracks_qw[3];

      for ( int i = 0; i < 3; ++i )
        {
          fvtxs_tracks_qx2[i] = 0.0;
          fvtxs_tracks_qy2[i] = 0.0;
          fvtxs_tracks_qx3[i] = 0.0;
          fvtxs_tracks_qy3[i] = 0.0;
          fvtxs_tracks_qx4[i] = 0.0;
          fvtxs_tracks_qy4[i] = 0.0;
          fvtxs_tracks_qx6[i] = 0.0;
          fvtxs_tracks_qy6[i] = 0.0;
          fvtxs_tracks_qw[i] = 0.0;
          fvtxn_tracks_qx2[i] = 0.0;
          fvtxn_tracks_qy2[i] = 0.0;
          fvtxn_tracks_qx3[i] = 0.0;
          fvtxn_tracks_qy3[i] = 0.0;
          fvtxn_tracks_qx4[i] = 0.0;
          fvtxn_tracks_qy4[i] = 0.0;
          fvtxn_tracks_qx6[i] = 0.0;
          fvtxn_tracks_qy6[i] = 0.0;
          fvtxn_tracks_qw[i] = 0.0;
        } // loop over layers

      int ntrack_south_inner = 0; // good_4_event
      int ntrack_north_inner = 0;
      int ntrack_south_outer = 0;
      int ntrack_north_outer = 0;

      float special_fvtx_tracks_qx2[8];
      float special_fvtx_tracks_qy2[8];
      float special_fvtx_tracks_qw[8];
      for ( int i = 0; i < 8; ++i )
	{
	  special_fvtx_tracks_qx2[i] = 0;
	  special_fvtx_tracks_qy2[i] = 0;
	  special_fvtx_tracks_qw[i] = 0;
	}

      bool fvtx_track_passes[nfvtxt];
      //int number_of_tracks_that_fail = 0;
      int number_of_tracks_that_pass = 0;
      for ( int i = 0; i < nfvtxt; ++i )
	{
	  fvtx_track_passes[i] = true;
	}
      for ( int i = 0; i < nfvtxt; ++i )
	{
	  for ( int j = i+1; j < nfvtxt; ++j )
	    {
	      //double phi1 = fphi[i];
	      double eta1 = feta[i];
	      //double phi2 = fphi[j];
	      double eta2 = feta[j];
	      if ( fabs(eta1-eta2) < 0.01 )
		{
		  //++number_of_tracks_that_fail; // trying to figure this out
		  fvtx_track_passes[i] = false;
		  fvtx_track_passes[j] = false;
		}
	    }
	  if ( fvtx_track_passes[i] == true ) ++number_of_tracks_that_pass;
	}
      //cout << "total tracks " << nfvtxt << endl;
      //cout << "tracks that pass " << number_of_tracks_that_pass << endl;
      //cout << "tracks that fail " << number_of_tracks_that_fail << endl;

      // --- first fvtx track loop
      for ( int i = 0; i < nfvtxt; ++i )
	{
	  // --- rotation now done in trees
	  //if ( !fvtx_track_passes[i] ) continue;
	  //int nhits = fnhits[i];
	  float phi = fphi[i];
	  float eta = feta[i];
	  //float dcax = fdcax[i];
	  //float dcay = fdcay[i];

	  //if ( nhits < 4 ) continue;
	  //if ( fabs(dcax) >  1.5 || fabs(dcay) > 1.5 ) continue;

	  //cout << dcax << " " << dcay << " " << nhits << endl;

	  int special_index = -1;
	  if ( eta > -3.0 && eta < -2.5 ) special_index = 0;
	  if ( eta > -2.5 && eta < -2.0 ) special_index = 1;
	  if ( eta > -2.0 && eta < -1.5 ) special_index = 2;
	  if ( eta > -1.5 && eta < -1.0 ) special_index = 3;
	  if ( eta > 1.0 && eta < 1.5 ) special_index = 4;
	  if ( eta > 1.5 && eta < 2.0 ) special_index = 5;
	  if ( eta > 2.0 && eta < 2.5 ) special_index = 6;
	  if ( eta > 2.5 && eta < 3.0 ) special_index = 7;
	  if ( special_index > -1 && special_index < 8 )
	    {
	      special_fvtx_tracks_qx2[special_index] += cos(2*phi);
	      special_fvtx_tracks_qy2[special_index] += sin(2*phi);
	      special_fvtx_tracks_qw[special_index] += 1.0;
	    }

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
	      fvtxs_tracks_qx6[0] += fvtx_weight * TMath::Cos(6*phi);
	      fvtxs_tracks_qy6[0] += fvtx_weight * TMath::Sin(6*phi);
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
	      fvtxn_tracks_qx6[0] += fvtx_weight * TMath::Cos(6*phi);
	      fvtxn_tracks_qy6[0] += fvtx_weight * TMath::Sin(6*phi);
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
	} // fvtx track loop

      for ( int i = 0; i < 8; ++i )
	{
	  for ( int j = 0; j < 8; ++j )
	    {
	      int kd = 0; // kronecker delta
	      if ( i == j ) kd = 1;
	      float ax = special_fvtx_tracks_qx2[i];
	      float bx = special_fvtx_tracks_qx2[j];
	      float ay = special_fvtx_tracks_qy2[i];
	      float by = special_fvtx_tracks_qy2[j];
	      float aw = special_fvtx_tracks_qw[i];
	      float bw = special_fvtx_tracks_qw[j];
	      float quant  = ( ax*bx + ay*by - kd*bw) / (aw*bw - kd*bw) ;
	      float eta = -9;
	      if ( j == 0 ) eta = -2.75;
	      if ( j == 1 ) eta = -2.25;
	      if ( j == 2 ) eta = -1.75;
	      if ( j == 3 ) eta = -1.25;
	      if ( j == 4 ) eta = 1.25;
	      if ( j == 5 ) eta = 1.75;
	      if ( j == 6 ) eta = 2.25;
	      if ( j == 7 ) eta = 2.75;
	      tp1f_special_fvtx_tracks_ab[i]->Fill(eta,quant);
	      if ( kd == 1 ) tp1f_special_fvtx_tracks_aa->Fill(eta,quant);
	      if ( kd == 1 ) tp1f_special_fvtx_tracks_aa_cos->Fill(eta,ax);
	      if ( kd == 1 ) tp1f_special_fvtx_tracks_aa_sin->Fill(eta,ay);
	    }
	}



      good_4_event = good_4_event && ( ntrack_south_inner > 0 ) && ( ntrack_south_outer > 0 ) && ( ntrack_north_inner > 0 ) && ( ntrack_north_outer > 0 ) ;
      int nfvtxt_south = ntrack_south_inner + ntrack_south_outer;
      int nfvtxt_north = ntrack_north_inner + ntrack_north_outer;


      float os_fvtxs_qw = fvtxs_qw[0];
      float os_fvtxs_qx2 = fvtxs_qx2[0];
      float os_fvtxs_qy2 = fvtxs_qy2[0];
      // float os_fvtxs_qx3 = fvtxs_qx3[0];
      // float os_fvtxs_qy3 = fvtxs_qy3[0];
      float os_fvtxs_qx4 = fvtxs_qx4[0];
      float os_fvtxs_qy4 = fvtxs_qy4[0];
      float os_fvtxs_qq2 = calc2_event(os_fvtxs_qx2,os_fvtxs_qy2,os_fvtxs_qw);
      //float os_fvtxs_qq3 = calc2_event(os_fvtxs_qx3,os_fvtxs_qy3,os_fvtxs_qw);
      nfvtxt_ac_fvtxs_c22->Fill(nfvtxt,os_fvtxs_qq2);
      nfvtxt_ac_fvtxs_cos21->Fill(nfvtxt,os_fvtxs_qx2/os_fvtxs_qw);
      nfvtxt_ac_fvtxs_sin21->Fill(nfvtxt,os_fvtxs_qy2/os_fvtxs_qw);
      TComplex tc_os_fvtxs_Q2(os_fvtxs_qx2,os_fvtxs_qy2);
      TComplex tc_os_fvtxs_Q4(os_fvtxs_qx4,os_fvtxs_qy4);
      float os_fvtxs_cossum2 = calccossum2_event(tc_os_fvtxs_Q2,tc_os_fvtxs_Q4,os_fvtxs_qw);
      float os_fvtxs_sinsum2 = calcsinsum2_event(tc_os_fvtxs_Q2,tc_os_fvtxs_Q4,os_fvtxs_qw);
      nfvtxt_ac_fvtxs_cossum22->Fill(nfvtxt,os_fvtxs_cossum2);
      nfvtxt_ac_fvtxs_sinsum22->Fill(nfvtxt,os_fvtxs_sinsum2);
      float os_fvtxs_cos23 = calccos3_event(tc_os_fvtxs_Q2,tc_os_fvtxs_Q4,os_fvtxs_qw);
      float os_fvtxs_sin23 = calcsin3_event(tc_os_fvtxs_Q2,tc_os_fvtxs_Q4,os_fvtxs_qw);
      nfvtxt_ac_fvtxs_cos23->Fill(nfvtxt,os_fvtxs_cos23);
      nfvtxt_ac_fvtxs_sin23->Fill(nfvtxt,os_fvtxs_sin23);

      float os_fvtxn_qw = fvtxn_qw[0];
      float os_fvtxn_qx2 = fvtxn_qx2[0];
      float os_fvtxn_qy2 = fvtxn_qy2[0];
      // float os_fvtxn_qx3 = fvtxn_qx3[0];
      // float os_fvtxn_qy3 = fvtxn_qy3[0];
      float os_fvtxn_qx4 = fvtxn_qx4[0];
      float os_fvtxn_qy4 = fvtxn_qy4[0];
      float os_fvtxn_qq2 = calc2_event(os_fvtxn_qx2,os_fvtxn_qy2,os_fvtxn_qw);
      //float os_fvtxn_qq3 = calc2_event(os_fvtxn_qx3,os_fvtxn_qy3,os_fvtxn_qw);
      nfvtxt_ac_fvtxn_c22->Fill(nfvtxt,os_fvtxn_qq2);
      nfvtxt_ac_fvtxn_cos21->Fill(nfvtxt,os_fvtxn_qx2/os_fvtxn_qw);
      nfvtxt_ac_fvtxn_sin21->Fill(nfvtxt,os_fvtxn_qy2/os_fvtxn_qw);
      TComplex tc_os_fvtxn_Q2(os_fvtxn_qx2,os_fvtxn_qy2);
      TComplex tc_os_fvtxn_Q4(os_fvtxn_qx4,os_fvtxn_qy4);
      float os_fvtxn_cossum2 = calccossum2_event(tc_os_fvtxn_Q2,tc_os_fvtxn_Q4,os_fvtxn_qw);
      float os_fvtxn_sinsum2 = calcsinsum2_event(tc_os_fvtxn_Q2,tc_os_fvtxn_Q4,os_fvtxn_qw);
      nfvtxt_ac_fvtxn_cossum22->Fill(nfvtxt,os_fvtxn_cossum2);
      nfvtxt_ac_fvtxn_sinsum22->Fill(nfvtxt,os_fvtxn_sinsum2);
      float os_fvtxn_cos23 = calccos3_event(tc_os_fvtxn_Q2,tc_os_fvtxn_Q4,os_fvtxn_qw);
      float os_fvtxn_sin23 = calcsin3_event(tc_os_fvtxn_Q2,tc_os_fvtxn_Q4,os_fvtxn_qw);
      nfvtxt_ac_fvtxn_cos23->Fill(nfvtxt,os_fvtxn_cos23);
      nfvtxt_ac_fvtxn_sin23->Fill(nfvtxt,os_fvtxn_sin23);

      // --- eta dependent clusters

      float os_fvtxs_ce0_qw = fvtxs_qw[6];
      float os_fvtxs_ce0_qx2 = fvtxs_qx2[6];
      float os_fvtxs_ce0_qy2 = fvtxs_qy2[6];
      // float os_fvtxs_ce0_qx3 = fvtxs_qx3[6];
      // float os_fvtxs_ce0_qy3 = fvtxs_qy3[6];
      float os_fvtxs_ce1_qw = fvtxs_qw[7];
      float os_fvtxs_ce1_qx2 = fvtxs_qx2[7];
      float os_fvtxs_ce1_qy2 = fvtxs_qy2[7];
      // float os_fvtxs_ce1_qx3 = fvtxs_qx3[7];
      // float os_fvtxs_ce1_qy3 = fvtxs_qy3[7];
      float os_fvtxn_ce0_qw = fvtxn_qw[6];
      float os_fvtxn_ce0_qx2 = fvtxn_qx2[6];
      float os_fvtxn_ce0_qy2 = fvtxn_qy2[6];
      // float os_fvtxn_ce0_qx3 = fvtxn_qx3[6];
      // float os_fvtxn_ce0_qy3 = fvtxn_qy3[6];
      float os_fvtxn_ce1_qw = fvtxn_qw[7];
      float os_fvtxn_ce1_qx2 = fvtxn_qx2[7];
      float os_fvtxn_ce1_qy2 = fvtxn_qy2[7];
      // float os_fvtxn_ce1_qx3 = fvtxn_qx3[7];
      // float os_fvtxn_ce1_qy3 = fvtxn_qy3[7];
      // ---
      float os_fvtxs_ce01_qq2 = ( os_fvtxs_ce0_qx2*os_fvtxs_ce1_qx2 + os_fvtxs_ce0_qy2*os_fvtxs_ce1_qy2 ) / ( os_fvtxs_ce0_qw*os_fvtxs_ce1_qw );
      //float os_fvtxs_ce01_qq3 = ( os_fvtxs_ce0_qx3*os_fvtxs_ce1_qx3 + os_fvtxs_ce0_qy3*os_fvtxs_ce1_qy3 ) / ( os_fvtxs_ce0_qw*os_fvtxs_ce1_qw );
      float os_fvtxn_ce01_qq2 = ( os_fvtxn_ce0_qx2*os_fvtxn_ce1_qx2 + os_fvtxn_ce0_qy2*os_fvtxn_ce1_qy2 ) / ( os_fvtxn_ce0_qw*os_fvtxn_ce1_qw );
      //float os_fvtxn_ce01_qq3 = ( os_fvtxn_ce0_qx3*os_fvtxn_ce1_qx3 + os_fvtxn_ce0_qy3*os_fvtxn_ce1_qy3 ) / ( os_fvtxn_ce0_qw*os_fvtxn_ce1_qw );
      // --- come back here to add some histograms and do some 4 particle stuff
      nfvtxt_ac_fvtxs_ce01_c22->Fill(nfvtxt,os_fvtxs_ce01_qq2);
      nfvtxt_ac_fvtxn_ce01_c22->Fill(nfvtxt,os_fvtxn_ce01_qq2);

      // --- now fvtx tracks

      float os_fvtxs_tracks_qw = fvtxs_tracks_qw[0];
      float os_fvtxs_tracks_qx2 = fvtxs_tracks_qx2[0];
      float os_fvtxs_tracks_qy2 = fvtxs_tracks_qy2[0];
      // float os_fvtxs_tracks_qx3 = fvtxs_tracks_qx3[0];
      // float os_fvtxs_tracks_qy3 = fvtxs_tracks_qy3[0];
      float os_fvtxs_tracks_qx4 = fvtxs_tracks_qx4[0];
      float os_fvtxs_tracks_qy4 = fvtxs_tracks_qy4[0];
      float os_fvtxs_tracks_qx6 = fvtxs_tracks_qx6[0];
      float os_fvtxs_tracks_qy6 = fvtxs_tracks_qy6[0];
      float os_fvtxs_tracks_qq2 = calc2_event(os_fvtxs_tracks_qx2,os_fvtxs_tracks_qy2,os_fvtxs_tracks_qw);
      //float os_fvtxs_tracks_qq3 = calc2_event(os_fvtxs_tracks_qx3,os_fvtxs_tracks_qy3,os_fvtxs_tracks_qw);
      nfvtxt_ac_fvtxs_tracks_c22->Fill(nfvtxt,os_fvtxs_tracks_qq2);
      nfvtxt_ac_fvtxs_tracks_cos21->Fill(nfvtxt,os_fvtxs_tracks_qx2/os_fvtxs_tracks_qw);
      nfvtxt_ac_fvtxs_tracks_sin21->Fill(nfvtxt,os_fvtxs_tracks_qy2/os_fvtxs_tracks_qw);
      nfvtxtsp_ac_fvtxs_tracks_c22->Fill(nfvtxt_south,os_fvtxs_tracks_qq2);
      nfvtxtsp_ac_fvtxs_tracks_cos21->Fill(nfvtxt_south,os_fvtxs_tracks_qx2/os_fvtxs_tracks_qw);
      nfvtxtsp_ac_fvtxs_tracks_sin21->Fill(nfvtxt_south,os_fvtxs_tracks_qy2/os_fvtxs_tracks_qw);
      TComplex tc_os_fvtxs_tracks_Q2(os_fvtxs_tracks_qx2,os_fvtxs_tracks_qy2);
      TComplex tc_os_fvtxs_tracks_Q4(os_fvtxs_tracks_qx4,os_fvtxs_tracks_qy4);
      TComplex tc_os_fvtxs_tracks_Q6(os_fvtxs_tracks_qx6,os_fvtxs_tracks_qy6);
      float os_fvtxs_tracks_cossum2 = calccossum2_event(tc_os_fvtxs_tracks_Q2,tc_os_fvtxs_tracks_Q4,os_fvtxs_tracks_qw);
      float os_fvtxs_tracks_sinsum2 = calcsinsum2_event(tc_os_fvtxs_tracks_Q2,tc_os_fvtxs_tracks_Q4,os_fvtxs_tracks_qw);
      nfvtxt_ac_fvtxs_tracks_cossum22->Fill(nfvtxt,os_fvtxs_tracks_cossum2);
      nfvtxt_ac_fvtxs_tracks_sinsum22->Fill(nfvtxt,os_fvtxs_tracks_sinsum2);
      nfvtxtsp_ac_fvtxs_tracks_cossum22->Fill(nfvtxt_south,os_fvtxs_tracks_cossum2);
      nfvtxtsp_ac_fvtxs_tracks_sinsum22->Fill(nfvtxt_south,os_fvtxs_tracks_sinsum2);
      float os_fvtxs_tracks_cos23 = calccos3_event(tc_os_fvtxs_tracks_Q2,tc_os_fvtxs_tracks_Q4,os_fvtxs_tracks_qw);
      float os_fvtxs_tracks_sin23 = calcsin3_event(tc_os_fvtxs_tracks_Q2,tc_os_fvtxs_tracks_Q4,os_fvtxs_tracks_qw);
      nfvtxt_ac_fvtxs_tracks_cos23->Fill(nfvtxt,os_fvtxs_tracks_cos23);
      nfvtxt_ac_fvtxs_tracks_sin23->Fill(nfvtxt,os_fvtxs_tracks_sin23);
      nfvtxtsp_ac_fvtxs_tracks_cos23->Fill(nfvtxt_south,os_fvtxs_tracks_cos23);
      nfvtxtsp_ac_fvtxs_tracks_sin23->Fill(nfvtxt_south,os_fvtxs_tracks_sin23);

      float os_fvtxn_tracks_qw = fvtxn_tracks_qw[0];
      float os_fvtxn_tracks_qx2 = fvtxn_tracks_qx2[0];
      float os_fvtxn_tracks_qy2 = fvtxn_tracks_qy2[0];
      // float os_fvtxn_tracks_qx3 = fvtxn_tracks_qx3[0];
      // float os_fvtxn_tracks_qy3 = fvtxn_tracks_qy3[0];
      float os_fvtxn_tracks_qx4 = fvtxn_tracks_qx4[0];
      float os_fvtxn_tracks_qy4 = fvtxn_tracks_qy4[0];
      float os_fvtxn_tracks_qx6 = fvtxn_tracks_qx6[0];
      float os_fvtxn_tracks_qy6 = fvtxn_tracks_qy6[0];
      float os_fvtxn_tracks_qq2 = calc2_event(os_fvtxn_tracks_qx2,os_fvtxn_tracks_qy2,os_fvtxn_tracks_qw);
      //float os_fvtxn_tracks_qq3 = calc2_event(os_fvtxn_tracks_qx3,os_fvtxn_tracks_qy3,os_fvtxn_tracks_qw);
      nfvtxt_ac_fvtxn_tracks_c22->Fill(nfvtxt,os_fvtxn_tracks_qq2);
      nfvtxt_ac_fvtxn_tracks_cos21->Fill(nfvtxt,os_fvtxn_tracks_qx2/os_fvtxn_tracks_qw);
      nfvtxt_ac_fvtxn_tracks_sin21->Fill(nfvtxt,os_fvtxn_tracks_qy2/os_fvtxn_tracks_qw);
      nfvtxtsp_ac_fvtxn_tracks_c22->Fill(nfvtxt_north,os_fvtxn_tracks_qq2);
      nfvtxtsp_ac_fvtxn_tracks_cos21->Fill(nfvtxt_north,os_fvtxn_tracks_qx2/os_fvtxn_tracks_qw);
      nfvtxtsp_ac_fvtxn_tracks_sin21->Fill(nfvtxt_north,os_fvtxn_tracks_qy2/os_fvtxn_tracks_qw);
      TComplex tc_os_fvtxn_tracks_Q2(os_fvtxn_tracks_qx2,os_fvtxn_tracks_qy2);
      TComplex tc_os_fvtxn_tracks_Q4(os_fvtxn_tracks_qx4,os_fvtxn_tracks_qy4);
      TComplex tc_os_fvtxn_tracks_Q6(os_fvtxn_tracks_qx6,os_fvtxn_tracks_qy6);
      float os_fvtxn_tracks_cossum2 = calccossum2_event(tc_os_fvtxn_tracks_Q2,tc_os_fvtxn_tracks_Q4,os_fvtxn_tracks_qw);
      float os_fvtxn_tracks_sinsum2 = calcsinsum2_event(tc_os_fvtxn_tracks_Q2,tc_os_fvtxn_tracks_Q4,os_fvtxn_tracks_qw);
      nfvtxt_ac_fvtxn_tracks_cossum22->Fill(nfvtxt,os_fvtxn_tracks_cossum2);
      nfvtxt_ac_fvtxn_tracks_sinsum22->Fill(nfvtxt,os_fvtxn_tracks_sinsum2);
      nfvtxtsp_ac_fvtxn_tracks_cossum22->Fill(nfvtxt_north,os_fvtxn_tracks_cossum2);
      nfvtxtsp_ac_fvtxn_tracks_sinsum22->Fill(nfvtxt_north,os_fvtxn_tracks_sinsum2);
      float os_fvtxn_tracks_cos23 = calccos3_event(tc_os_fvtxn_tracks_Q2,tc_os_fvtxn_tracks_Q4,os_fvtxn_tracks_qw);
      float os_fvtxn_tracks_sin23 = calcsin3_event(tc_os_fvtxn_tracks_Q2,tc_os_fvtxn_tracks_Q4,os_fvtxn_tracks_qw);
      nfvtxt_ac_fvtxn_tracks_cos23->Fill(nfvtxt,os_fvtxn_tracks_cos23);
      nfvtxt_ac_fvtxn_tracks_sin23->Fill(nfvtxt,os_fvtxn_tracks_sin23);
      nfvtxtsp_ac_fvtxn_tracks_cos23->Fill(nfvtxt_north,os_fvtxn_tracks_cos23);
      nfvtxtsp_ac_fvtxn_tracks_sin23->Fill(nfvtxt_north,os_fvtxn_tracks_sin23);

      // --- combined fvtx clusters

      float os_fvtxc_qx2 = os_fvtxs_qx2 + os_fvtxn_qx2;
      float os_fvtxc_qy2 = os_fvtxs_qy2 + os_fvtxn_qy2;
      // float os_fvtxc_qx3 = os_fvtxs_qx3 + os_fvtxn_qx3;
      // float os_fvtxc_qy3 = os_fvtxs_qy3 + os_fvtxn_qy3;
      float os_fvtxc_qx4 = os_fvtxs_qx4 + os_fvtxn_qx4;
      float os_fvtxc_qy4 = os_fvtxs_qy4 + os_fvtxn_qy4;
      float os_fvtxc_qw = os_fvtxs_qw + os_fvtxn_qw;
      float os_fvtxc_qq2 = calc2_event(os_fvtxc_qx2,os_fvtxc_qy2,os_fvtxc_qw);
      //float os_fvtxc_qq3 = calc2_event(os_fvtxc_qx3,os_fvtxc_qy3,os_fvtxc_qw);
      nfvtxt_ac_fvtxc_c22->Fill(nfvtxt,os_fvtxc_qq2);
      nfvtxt_ac_fvtxc_cos21->Fill(nfvtxt,os_fvtxc_qx2/os_fvtxc_qw);
      nfvtxt_ac_fvtxc_sin21->Fill(nfvtxt,os_fvtxc_qy2/os_fvtxc_qw);
      TComplex tc_os_fvtxc_Q2(os_fvtxc_qx2,os_fvtxc_qy2);
      TComplex tc_os_fvtxc_Q4(os_fvtxc_qx4,os_fvtxc_qy4);
      float os_fvtxc_cossum2 = calccossum2_event(tc_os_fvtxc_Q2,tc_os_fvtxc_Q4,os_fvtxc_qw);
      float os_fvtxc_sinsum2 = calcsinsum2_event(tc_os_fvtxc_Q2,tc_os_fvtxc_Q4,os_fvtxc_qw);
      nfvtxt_ac_fvtxc_cossum22->Fill(nfvtxt,os_fvtxc_cossum2);
      nfvtxt_ac_fvtxc_sinsum22->Fill(nfvtxt,os_fvtxc_sinsum2);
      float os_fvtxc_cos23 = calccos3_event(tc_os_fvtxc_Q2,tc_os_fvtxc_Q4,os_fvtxc_qw);
      float os_fvtxc_sin23 = calcsin3_event(tc_os_fvtxc_Q2,tc_os_fvtxc_Q4,os_fvtxc_qw);
      nfvtxt_ac_fvtxc_cos23->Fill(nfvtxt,os_fvtxc_cos23);
      nfvtxt_ac_fvtxc_sin23->Fill(nfvtxt,os_fvtxc_sin23);

      // --- combined fvtx tracks

      float os_fvtxc_tracks_qx2 = os_fvtxs_tracks_qx2 + os_fvtxn_tracks_qx2;
      float os_fvtxc_tracks_qy2 = os_fvtxs_tracks_qy2 + os_fvtxn_tracks_qy2;
      // float os_fvtxc_tracks_qx3 = os_fvtxs_tracks_qx3 + os_fvtxn_tracks_qx3;
      // float os_fvtxc_tracks_qy3 = os_fvtxs_tracks_qy3 + os_fvtxn_tracks_qy3;
      float os_fvtxc_tracks_qx4 = os_fvtxs_tracks_qx4 + os_fvtxn_tracks_qx4;
      float os_fvtxc_tracks_qy4 = os_fvtxs_tracks_qy4 + os_fvtxn_tracks_qy4;
      float os_fvtxc_tracks_qx6 = os_fvtxs_tracks_qx6 + os_fvtxn_tracks_qx6;
      float os_fvtxc_tracks_qy6 = os_fvtxs_tracks_qy6 + os_fvtxn_tracks_qy6;
      float os_fvtxc_tracks_qw = os_fvtxs_tracks_qw + os_fvtxn_tracks_qw;
      float os_fvtxc_tracks_qq2 = calc2_event(os_fvtxc_tracks_qx2,os_fvtxc_tracks_qy2,os_fvtxc_tracks_qw);
      //float os_fvtxc_tracks_qq3 = calc2_event(os_fvtxc_tracks_qx3,os_fvtxc_tracks_qy3,os_fvtxc_tracks_qw);
      nfvtxt_ac_fvtxc_tracks_c22->Fill(nfvtxt,os_fvtxc_tracks_qq2);
      nfvtxt_ac_fvtxc_tracks_cos21->Fill(nfvtxt,os_fvtxc_tracks_qx2/os_fvtxc_tracks_qw);
      nfvtxt_ac_fvtxc_tracks_sin21->Fill(nfvtxt,os_fvtxc_tracks_qy2/os_fvtxc_tracks_qw);
      nfvtxtsp_ac_fvtxc_tracks_c22->Fill(nfvtxt,os_fvtxc_tracks_qq2);
      nfvtxtsp_ac_fvtxc_tracks_cos21->Fill(nfvtxt,os_fvtxc_tracks_qx2/os_fvtxc_tracks_qw);
      nfvtxtsp_ac_fvtxc_tracks_sin21->Fill(nfvtxt,os_fvtxc_tracks_qy2/os_fvtxc_tracks_qw);
      TComplex tc_os_fvtxc_tracks_Q2(os_fvtxc_tracks_qx2,os_fvtxc_tracks_qy2);
      TComplex tc_os_fvtxc_tracks_Q4(os_fvtxc_tracks_qx4,os_fvtxc_tracks_qy4);
      TComplex tc_os_fvtxc_tracks_Q6(os_fvtxc_tracks_qx6,os_fvtxc_tracks_qy6);
      float os_fvtxc_tracks_cossum2 = calccossum2_event(tc_os_fvtxc_tracks_Q2,tc_os_fvtxc_tracks_Q4,os_fvtxc_tracks_qw);
      float os_fvtxc_tracks_sinsum2 = calcsinsum2_event(tc_os_fvtxc_tracks_Q2,tc_os_fvtxc_tracks_Q4,os_fvtxc_tracks_qw);
      nfvtxt_ac_fvtxc_tracks_cossum22->Fill(nfvtxt,os_fvtxc_tracks_cossum2);
      nfvtxt_ac_fvtxc_tracks_sinsum22->Fill(nfvtxt,os_fvtxc_tracks_sinsum2);
      nfvtxtsp_ac_fvtxc_tracks_cossum22->Fill(nfvtxt,os_fvtxc_tracks_cossum2);
      nfvtxtsp_ac_fvtxc_tracks_sinsum22->Fill(nfvtxt,os_fvtxc_tracks_sinsum2);
      float os_fvtxc_tracks_cos23 = calccos3_event(tc_os_fvtxc_tracks_Q2,tc_os_fvtxc_tracks_Q4,os_fvtxc_tracks_qw);
      float os_fvtxc_tracks_sin23 = calcsin3_event(tc_os_fvtxc_tracks_Q2,tc_os_fvtxc_tracks_Q4,os_fvtxc_tracks_qw);
      nfvtxt_ac_fvtxc_tracks_cos23->Fill(nfvtxt,os_fvtxc_tracks_cos23);
      nfvtxt_ac_fvtxc_tracks_sin23->Fill(nfvtxt,os_fvtxc_tracks_sin23);
      nfvtxtsp_ac_fvtxc_tracks_cos23->Fill(nfvtxt,os_fvtxc_tracks_cos23);
      nfvtxtsp_ac_fvtxc_tracks_sin23->Fill(nfvtxt,os_fvtxc_tracks_sin23);

      // --- have a look at some different correlations

      // --- fvtxs-fvtxn
      float os_fvtxsfvtxn_qq2 = ( (os_fvtxs_qx2*os_fvtxn_qx2) + (os_fvtxs_qy2*os_fvtxn_qy2) ) / ( os_fvtxs_qw*os_fvtxn_qw );
      //float os_fvtxsfvtxn_qq3 = ( (os_fvtxs_qx3*os_fvtxn_qx3) + (os_fvtxs_qy3*os_fvtxn_qy3) ) / ( os_fvtxs_qw*os_fvtxn_qw );
      nfvtxt_ac_fvtxsfvtxn_c22->Fill(nfvtxt,os_fvtxsfvtxn_qq2); // see below
      // ---

      float os_fvtxsfvtxn_tracks_qq2 = ( (os_fvtxs_tracks_qx2*os_fvtxn_tracks_qx2) + (os_fvtxs_tracks_qy2*os_fvtxn_tracks_qy2) ) / ( os_fvtxs_tracks_qw*os_fvtxn_tracks_qw );
      //float os_fvtxsfvtxn_tracks_qq3 = ( (os_fvtxs_tracks_qx3*os_fvtxn_tracks_qx3) + (os_fvtxs_tracks_qy3*os_fvtxn_tracks_qy3) ) / ( os_fvtxs_tracks_qw*os_fvtxn_tracks_qw );
      nfvtxt_ac_fvtxsfvtxn_tracks_c22->Fill(nfvtxt,os_fvtxsfvtxn_tracks_qq2);



      // ------------------------------------
      // --- now for the q-vector recentering
      // ------------------------------------

      // --- need to get bin number based on nfvtxt, then get bincontent to set the offset
      int binn = nfvtxt_tracks_south_qx2->GetXaxis()->FindBin(nfvtxt);
      float offset_tracks_south_qx2 = nfvtxt_tracks_south_qx2->GetBinContent(binn);
      //float offset_tracks_south_qx3 = nfvtxt_tracks_south_qx3->GetBinContent(binn);
      float offset_tracks_south_qx4 = nfvtxt_tracks_south_qx4->GetBinContent(binn);
      float offset_tracks_south_qx6 = nfvtxt_tracks_south_qx6->GetBinContent(binn);
      float offset_tracks_south_qy2 = nfvtxt_tracks_south_qy2->GetBinContent(binn);
      //float offset_tracks_south_qy3 = nfvtxt_tracks_south_qy3->GetBinContent(binn);
      float offset_tracks_south_qy4 = nfvtxt_tracks_south_qy4->GetBinContent(binn);
      float offset_tracks_south_qy6 = nfvtxt_tracks_south_qy6->GetBinContent(binn);
      float offset_tracks_north_qx2 = nfvtxt_tracks_north_qx2->GetBinContent(binn);
      //float offset_tracks_north_qx3 = nfvtxt_tracks_north_qx3->GetBinContent(binn);
      float offset_tracks_north_qx4 = nfvtxt_tracks_north_qx4->GetBinContent(binn);
      float offset_tracks_north_qx6 = nfvtxt_tracks_north_qx6->GetBinContent(binn);
      float offset_tracks_north_qy2 = nfvtxt_tracks_north_qy2->GetBinContent(binn);
      //float offset_tracks_north_qy3 = nfvtxt_tracks_north_qy3->GetBinContent(binn);
      float offset_tracks_north_qy4 = nfvtxt_tracks_north_qy4->GetBinContent(binn);
      float offset_tracks_north_qy6 = nfvtxt_tracks_north_qy6->GetBinContent(binn);

      // if ( say_event ) cout << "bin number is " << binn << " and number of fvtx tracks is " << nfvtxt << endl;
      // if ( say_event ) cout << "qx2 south offset is " << offset_tracks_south_qx2 << endl;

      float zzyzx_fvtxs_tracks_qw = fvtxs_tracks_qw[0];
      float zzyzx_fvtxs_tracks_qx2 = fvtxs_tracks_qx2[0] - offset_tracks_south_qx2;
      float zzyzx_fvtxs_tracks_qy2 = fvtxs_tracks_qy2[0] - offset_tracks_south_qy2;
      // float zzyzx_fvtxs_tracks_qx3 = fvtxs_tracks_qx3[0] - offset_tracks_south_qx3;
      // float zzyzx_fvtxs_tracks_qy3 = fvtxs_tracks_qy3[0] - offset_tracks_south_qy3;
      float zzyzx_fvtxs_tracks_qx4 = fvtxs_tracks_qx4[0] - offset_tracks_south_qx4;
      float zzyzx_fvtxs_tracks_qy4 = fvtxs_tracks_qy4[0] - offset_tracks_south_qy4;
      float zzyzx_fvtxs_tracks_qx6 = fvtxs_tracks_qx6[0] - offset_tracks_south_qx6;
      float zzyzx_fvtxs_tracks_qy6 = fvtxs_tracks_qy6[0] - offset_tracks_south_qy6;
      float zzyzx_fvtxs_tracks_qq2 = calc2_event(zzyzx_fvtxs_tracks_qx2,zzyzx_fvtxs_tracks_qy2,zzyzx_fvtxs_tracks_qw);
      //float zzyzx_fvtxs_tracks_qq3 = calc2_event(zzyzx_fvtxs_tracks_qx3,zzyzx_fvtxs_tracks_qy3,zzyzx_fvtxs_tracks_qw);
      nfvtxt_zzyzx_fvtxs_tracks_c22->Fill(nfvtxt,zzyzx_fvtxs_tracks_qq2);
      nfvtxt_zzyzx_fvtxs_tracks_cos21->Fill(nfvtxt,zzyzx_fvtxs_tracks_qx2/zzyzx_fvtxs_tracks_qw);
      nfvtxt_zzyzx_fvtxs_tracks_sin21->Fill(nfvtxt,zzyzx_fvtxs_tracks_qy2/zzyzx_fvtxs_tracks_qw);
      TComplex tc_zzyzx_fvtxs_tracks_Q2(zzyzx_fvtxs_tracks_qx2,zzyzx_fvtxs_tracks_qy2);
      TComplex tc_zzyzx_fvtxs_tracks_Q4(zzyzx_fvtxs_tracks_qx4,zzyzx_fvtxs_tracks_qy4);
      TComplex tc_zzyzx_fvtxs_tracks_Q6(zzyzx_fvtxs_tracks_qx6,zzyzx_fvtxs_tracks_qy6);
      float zzyzx_fvtxs_tracks_cossum2 = calccossum2_event(tc_zzyzx_fvtxs_tracks_Q2,tc_zzyzx_fvtxs_tracks_Q4,zzyzx_fvtxs_tracks_qw);
      float zzyzx_fvtxs_tracks_sinsum2 = calcsinsum2_event(tc_zzyzx_fvtxs_tracks_Q2,tc_zzyzx_fvtxs_tracks_Q4,zzyzx_fvtxs_tracks_qw);
      nfvtxt_zzyzx_fvtxs_tracks_cossum22->Fill(nfvtxt,zzyzx_fvtxs_tracks_cossum2);
      nfvtxt_zzyzx_fvtxs_tracks_sinsum22->Fill(nfvtxt,zzyzx_fvtxs_tracks_sinsum2);
      float zzyzx_fvtxs_tracks_cos23 = calccos3_event(tc_zzyzx_fvtxs_tracks_Q2,tc_zzyzx_fvtxs_tracks_Q4,zzyzx_fvtxs_tracks_qw);
      float zzyzx_fvtxs_tracks_sin23 = calcsin3_event(tc_zzyzx_fvtxs_tracks_Q2,tc_zzyzx_fvtxs_tracks_Q4,zzyzx_fvtxs_tracks_qw);
      nfvtxt_zzyzx_fvtxs_tracks_cos23->Fill(nfvtxt,zzyzx_fvtxs_tracks_cos23);
      nfvtxt_zzyzx_fvtxs_tracks_sin23->Fill(nfvtxt,zzyzx_fvtxs_tracks_sin23);

      float zzyzx_fvtxn_tracks_qw = fvtxn_tracks_qw[0];
      float zzyzx_fvtxn_tracks_qx2 = fvtxn_tracks_qx2[0] - offset_tracks_north_qx2;
      float zzyzx_fvtxn_tracks_qy2 = fvtxn_tracks_qy2[0] - offset_tracks_north_qy2;
      // float zzyzx_fvtxn_tracks_qx3 = fvtxn_tracks_qx3[0] - offset_tracks_north_qx3;
      // float zzyzx_fvtxn_tracks_qy3 = fvtxn_tracks_qy3[0] - offset_tracks_north_qy3;
      float zzyzx_fvtxn_tracks_qx4 = fvtxn_tracks_qx4[0] - offset_tracks_north_qx4;
      float zzyzx_fvtxn_tracks_qy4 = fvtxn_tracks_qy4[0] - offset_tracks_north_qy4;
      float zzyzx_fvtxn_tracks_qx6 = fvtxn_tracks_qx6[0] - offset_tracks_north_qx6;
      float zzyzx_fvtxn_tracks_qy6 = fvtxn_tracks_qy6[0] - offset_tracks_north_qy6;
      float zzyzx_fvtxn_tracks_qq2 = calc2_event(zzyzx_fvtxn_tracks_qx2,zzyzx_fvtxn_tracks_qy2,zzyzx_fvtxn_tracks_qw);
      //float zzyzx_fvtxn_tracks_qq3 = calc2_event(zzyzx_fvtxn_tracks_qx3,zzyzx_fvtxn_tracks_qy3,zzyzx_fvtxn_tracks_qw);
      nfvtxt_zzyzx_fvtxn_tracks_c22->Fill(nfvtxt,zzyzx_fvtxn_tracks_qq2);
      nfvtxt_zzyzx_fvtxn_tracks_cos21->Fill(nfvtxt,zzyzx_fvtxn_tracks_qx2/zzyzx_fvtxn_tracks_qw);
      nfvtxt_zzyzx_fvtxn_tracks_sin21->Fill(nfvtxt,zzyzx_fvtxn_tracks_qy2/zzyzx_fvtxn_tracks_qw);
      TComplex tc_zzyzx_fvtxn_tracks_Q2(zzyzx_fvtxn_tracks_qx2,zzyzx_fvtxn_tracks_qy2);
      TComplex tc_zzyzx_fvtxn_tracks_Q4(zzyzx_fvtxn_tracks_qx4,zzyzx_fvtxn_tracks_qy4);
      TComplex tc_zzyzx_fvtxn_tracks_Q6(zzyzx_fvtxn_tracks_qx6,zzyzx_fvtxn_tracks_qy6);
      float zzyzx_fvtxn_tracks_cossum2 = calccossum2_event(tc_zzyzx_fvtxn_tracks_Q2,tc_zzyzx_fvtxn_tracks_Q4,zzyzx_fvtxn_tracks_qw);
      float zzyzx_fvtxn_tracks_sinsum2 = calcsinsum2_event(tc_zzyzx_fvtxn_tracks_Q2,tc_zzyzx_fvtxn_tracks_Q4,zzyzx_fvtxn_tracks_qw);
      nfvtxt_zzyzx_fvtxn_tracks_cossum22->Fill(nfvtxt,zzyzx_fvtxn_tracks_cossum2);
      nfvtxt_zzyzx_fvtxn_tracks_sinsum22->Fill(nfvtxt,zzyzx_fvtxn_tracks_sinsum2);
      float zzyzx_fvtxn_tracks_cos23 = calccos3_event(tc_zzyzx_fvtxn_tracks_Q2,tc_zzyzx_fvtxn_tracks_Q4,zzyzx_fvtxn_tracks_qw);
      float zzyzx_fvtxn_tracks_sin23 = calcsin3_event(tc_zzyzx_fvtxn_tracks_Q2,tc_zzyzx_fvtxn_tracks_Q4,zzyzx_fvtxn_tracks_qw);
      nfvtxt_zzyzx_fvtxn_tracks_cos23->Fill(nfvtxt,zzyzx_fvtxn_tracks_cos23);
      nfvtxt_zzyzx_fvtxn_tracks_sin23->Fill(nfvtxt,zzyzx_fvtxn_tracks_sin23);

      float zzyzx_fvtxc_tracks_qx2 = zzyzx_fvtxs_tracks_qx2 + zzyzx_fvtxn_tracks_qx2;
      float zzyzx_fvtxc_tracks_qy2 = zzyzx_fvtxs_tracks_qy2 + zzyzx_fvtxn_tracks_qy2;
      // float zzyzx_fvtxc_tracks_qx3 = zzyzx_fvtxs_tracks_qx3 + zzyzx_fvtxn_tracks_qx3;
      // float zzyzx_fvtxc_tracks_qy3 = zzyzx_fvtxs_tracks_qy3 + zzyzx_fvtxn_tracks_qy3;
      float zzyzx_fvtxc_tracks_qx4 = zzyzx_fvtxs_tracks_qx4 + zzyzx_fvtxn_tracks_qx4;
      float zzyzx_fvtxc_tracks_qy4 = zzyzx_fvtxs_tracks_qy4 + zzyzx_fvtxn_tracks_qy4;
      float zzyzx_fvtxc_tracks_qx6 = zzyzx_fvtxs_tracks_qx6 + zzyzx_fvtxn_tracks_qx6;
      float zzyzx_fvtxc_tracks_qy6 = zzyzx_fvtxs_tracks_qy6 + zzyzx_fvtxn_tracks_qy6;
      float zzyzx_fvtxc_tracks_qw = zzyzx_fvtxs_tracks_qw + zzyzx_fvtxn_tracks_qw;
      float zzyzx_fvtxc_tracks_qq2 = calc2_event(zzyzx_fvtxc_tracks_qx2,zzyzx_fvtxc_tracks_qy2,zzyzx_fvtxc_tracks_qw);
      //float zzyzx_fvtxc_tracks_qq3 = calc2_event(zzyzx_fvtxc_tracks_qx3,zzyzx_fvtxc_tracks_qy3,zzyzx_fvtxc_tracks_qw);
      nfvtxt_zzyzx_fvtxc_tracks_c22->Fill(nfvtxt,zzyzx_fvtxc_tracks_qq2);
      nfvtxt_zzyzx_fvtxc_tracks_cos21->Fill(nfvtxt,zzyzx_fvtxc_tracks_qx2/zzyzx_fvtxc_tracks_qw);
      nfvtxt_zzyzx_fvtxc_tracks_sin21->Fill(nfvtxt,zzyzx_fvtxc_tracks_qy2/zzyzx_fvtxc_tracks_qw);
      TComplex tc_zzyzx_fvtxc_tracks_Q2(zzyzx_fvtxc_tracks_qx2,zzyzx_fvtxc_tracks_qy2);
      TComplex tc_zzyzx_fvtxc_tracks_Q4(zzyzx_fvtxc_tracks_qx4,zzyzx_fvtxc_tracks_qy4);
      TComplex tc_zzyzx_fvtxc_tracks_Q6(zzyzx_fvtxc_tracks_qx6,zzyzx_fvtxc_tracks_qy6);
      float zzyzx_fvtxc_tracks_cossum2 = calccossum2_event(tc_zzyzx_fvtxc_tracks_Q2,tc_zzyzx_fvtxc_tracks_Q4,zzyzx_fvtxc_tracks_qw);
      float zzyzx_fvtxc_tracks_sinsum2 = calcsinsum2_event(tc_zzyzx_fvtxc_tracks_Q2,tc_zzyzx_fvtxc_tracks_Q4,zzyzx_fvtxc_tracks_qw);
      nfvtxt_zzyzx_fvtxc_tracks_cossum22->Fill(nfvtxt,zzyzx_fvtxc_tracks_cossum2);
      nfvtxt_zzyzx_fvtxc_tracks_sinsum22->Fill(nfvtxt,zzyzx_fvtxc_tracks_sinsum2);
      float zzyzx_fvtxc_tracks_cos23 = calccos3_event(tc_zzyzx_fvtxc_tracks_Q2,tc_zzyzx_fvtxc_tracks_Q4,zzyzx_fvtxc_tracks_qw);
      float zzyzx_fvtxc_tracks_sin23 = calcsin3_event(tc_zzyzx_fvtxc_tracks_Q2,tc_zzyzx_fvtxc_tracks_Q4,zzyzx_fvtxc_tracks_qw);
      nfvtxt_zzyzx_fvtxc_tracks_cos23->Fill(nfvtxt,zzyzx_fvtxc_tracks_cos23);
      nfvtxt_zzyzx_fvtxc_tracks_sin23->Fill(nfvtxt,zzyzx_fvtxc_tracks_sin23);

      // --- fvtxs-fvtxn
      float zzyzx_fvtxsfvtxn_tracks_qq2 = ( (zzyzx_fvtxs_tracks_qx2*zzyzx_fvtxn_tracks_qx2) + (zzyzx_fvtxs_tracks_qy2*zzyzx_fvtxn_tracks_qy2) ) / ( zzyzx_fvtxs_tracks_qw*zzyzx_fvtxn_tracks_qw );
      //float zzyzx_fvtxsfvtxn_tracks_qq3 = ( (zzyzx_fvtxs_tracks_qx3*zzyzx_fvtxn_tracks_qx3) + (zzyzx_fvtxs_tracks_qy3*zzyzx_fvtxn_tracks_qy3) ) / ( zzyzx_fvtxs_tracks_qw*zzyzx_fvtxn_tracks_qw );
      nfvtxt_zzyzx_fvtxsfvtxn_tracks_c22->Fill(nfvtxt,zzyzx_fvtxsfvtxn_tracks_qq2);


      // --- now have a look at some 4 particle cumulants
      if ( good_4_event )
	{
	  // ---
	  float os_fvtxs_tracks_six = calc6_event(tc_os_fvtxs_tracks_Q2,tc_os_fvtxs_tracks_Q4,tc_os_fvtxs_tracks_Q6,os_fvtxs_tracks_qw);
	  float os_fvtxn_tracks_six = calc6_event(tc_os_fvtxn_tracks_Q2,tc_os_fvtxn_tracks_Q4,tc_os_fvtxn_tracks_Q6,os_fvtxn_tracks_qw);
	  float os_fvtxc_tracks_six = calc6_event(tc_os_fvtxc_tracks_Q2,tc_os_fvtxc_tracks_Q4,tc_os_fvtxc_tracks_Q6,os_fvtxc_tracks_qw);
	  nfvtxt_ac_fvtxs_tracks_c26->Fill(nfvtxt,os_fvtxs_tracks_six);
	  nfvtxt_ac_fvtxn_tracks_c26->Fill(nfvtxt,os_fvtxn_tracks_six);
	  nfvtxt_ac_fvtxc_tracks_c26->Fill(nfvtxt,os_fvtxc_tracks_six);

	  float os_fvtxs_tracks_qqqq4 = calc4_event(os_fvtxs_tracks_qx2,os_fvtxs_tracks_qy2,os_fvtxs_tracks_qx4,os_fvtxs_tracks_qy4,os_fvtxs_tracks_qw);
	  float os_fvtxn_tracks_qqqq4 = calc4_event(os_fvtxn_tracks_qx2,os_fvtxn_tracks_qy2,os_fvtxn_tracks_qx4,os_fvtxn_tracks_qy4,os_fvtxn_tracks_qw);
	  float os_fvtxc_tracks_qqqq4 = calc4_event(os_fvtxc_tracks_qx2,os_fvtxc_tracks_qy2,os_fvtxc_tracks_qx4,os_fvtxc_tracks_qy4,os_fvtxc_tracks_qw);
	  nfvtxt_ac_fvtxs_tracks_c24->Fill(nfvtxt,os_fvtxs_tracks_qqqq4);
	  nfvtxt_ac_fvtxn_tracks_c24->Fill(nfvtxt,os_fvtxn_tracks_qqqq4);
	  nfvtxt_ac_fvtxc_tracks_c24->Fill(nfvtxt,os_fvtxc_tracks_qqqq4);

	  nfvtxtsp_ac_fvtxs_tracks_c24->Fill(nfvtxt_south,os_fvtxs_tracks_qqqq4);
	  nfvtxtsp_ac_fvtxn_tracks_c24->Fill(nfvtxt_north,os_fvtxn_tracks_qqqq4);
	  nfvtxtsp_ac_fvtxc_tracks_c24->Fill(nfvtxt,os_fvtxc_tracks_qqqq4);

	  nfvtxt_ac_fvtxs_tracks_cov24->Fill(nfvtxt,os_fvtxs_tracks_qqqq4*os_fvtxs_tracks_qq2);
	  nfvtxt_ac_fvtxn_tracks_cov24->Fill(nfvtxt,os_fvtxn_tracks_qqqq4*os_fvtxn_tracks_qq2);
	  nfvtxt_ac_fvtxc_tracks_cov24->Fill(nfvtxt,os_fvtxc_tracks_qqqq4*os_fvtxc_tracks_qq2);




	  float os_fvtxs_qqqq4 = calc4_event(os_fvtxs_qx2,os_fvtxs_qy2,os_fvtxs_qx4,os_fvtxs_qy4,os_fvtxs_qw);
	  float os_fvtxn_qqqq4 = calc4_event(os_fvtxn_qx2,os_fvtxn_qy2,os_fvtxn_qx4,os_fvtxn_qy4,os_fvtxn_qw);
	  float os_fvtxc_qqqq4 = calc4_event(os_fvtxc_qx2,os_fvtxc_qy2,os_fvtxc_qx4,os_fvtxc_qy4,os_fvtxc_qw);
	  nfvtxt_ac_fvtxs_c24->Fill(nfvtxt,os_fvtxs_qqqq4);
	  nfvtxt_ac_fvtxn_c24->Fill(nfvtxt,os_fvtxn_qqqq4);
	  nfvtxt_ac_fvtxc_c24->Fill(nfvtxt,os_fvtxc_qqqq4);
	  nfvtxt_ac_fvtxsfvtxn_c24a->Fill(nfvtxt,os_fvtxsfvtxn_qq2*os_fvtxsfvtxn_qq2); // doesn't account for cross terms
	  nfvtxt_ac_fvtxsfvtxn_c24b->Fill(nfvtxt,os_fvtxs_qq2*os_fvtxs_qq2*os_fvtxn_qq2*os_fvtxn_qq2); // i think this power counts to v^8
	  nfvtxt_ac_fvtxsfvtxn_c24c->Fill(nfvtxt,os_fvtxsfvtxn_qq2*os_fvtxs_qq2*os_fvtxn_qq2); // i think this power counts to v^6
	  nfvtxt_ac_fvtxsfvtxn_c24d->Fill(nfvtxt,os_fvtxs_qq2*os_fvtxn_qq2); // i think this is right (but problem with autocorrelations)
	  nfvtxt_ac_fvtxsfvtxn_c24->Fill(nfvtxt,os_fvtxs_qq2*os_fvtxn_qq2); // might as well use this as an anchor


	  nfvtxt_ac_fvtxsfvtxn_tracks_c24a->Fill(nfvtxt,os_fvtxsfvtxn_tracks_qq2*os_fvtxsfvtxn_tracks_qq2); // doesn't account for cross terms
	  nfvtxt_ac_fvtxsfvtxn_tracks_c24b->Fill(nfvtxt,os_fvtxs_tracks_qq2*os_fvtxs_tracks_qq2*os_fvtxn_tracks_qq2*os_fvtxn_tracks_qq2); // v^8 ?
	  nfvtxt_ac_fvtxsfvtxn_tracks_c24c->Fill(nfvtxt,os_fvtxsfvtxn_tracks_qq2*os_fvtxs_tracks_qq2*os_fvtxn_tracks_qq2); // v^6 ?
	  nfvtxt_ac_fvtxsfvtxn_tracks_c24d->Fill(nfvtxt,os_fvtxs_tracks_qq2*os_fvtxn_tracks_qq2); // i think this is right (and tracks shouldn't have autocorrelation issues)
	  nfvtxt_ac_fvtxsfvtxn_tracks_c24->Fill(nfvtxt,os_fvtxs_tracks_qq2*os_fvtxn_tracks_qq2); // use this as an anchor

	  // --- now let's look at some subevent 4pc stuff
	  TComplex tc_os_fvtxs_ce0(os_fvtxs_ce0_qx2,os_fvtxs_ce0_qy2);
	  TComplex tc_os_fvtxs_ce1(os_fvtxs_ce1_qx2,os_fvtxs_ce1_qy2);
	  TComplex tc_os_fvtxn_ce0(os_fvtxn_ce0_qx2,os_fvtxn_ce0_qy2);
	  TComplex tc_os_fvtxn_ce1(os_fvtxn_ce1_qx2,os_fvtxn_ce1_qy2);
	  // --- first
	  TComplex tc_os_fvtxns_ce01 = tc_os_fvtxs_ce0 * tc_os_fvtxs_ce1 * TComplex::Conjugate(tc_os_fvtxn_ce0) * TComplex::Conjugate(tc_os_fvtxn_ce1);
	  float os_fvtxns_ce01_qqqq4 = tc_os_fvtxns_ce01.Re();
	  float norm = os_fvtxs_ce0_qw*os_fvtxs_ce1_qw*os_fvtxn_ce0_qw*os_fvtxn_ce1_qw;
	  os_fvtxns_ce01_qqqq4 /= norm;
	  nfvtxt_ac_fvtxsfvtxn_ce01_c24->Fill(nfvtxt,os_fvtxns_ce01_qqqq4); // seems redundant but i sort of need an anchor result
	  nfvtxt_ac_fvtxsfvtxn_ce01_c24a->Fill(nfvtxt,os_fvtxns_ce01_qqqq4);
	  // --- second
	  tc_os_fvtxns_ce01 = tc_os_fvtxs_ce0 * TComplex::Conjugate(tc_os_fvtxs_ce1) * tc_os_fvtxn_ce0 * TComplex::Conjugate(tc_os_fvtxn_ce1);
	  os_fvtxns_ce01_qqqq4 = tc_os_fvtxns_ce01.Re();
	  os_fvtxns_ce01_qqqq4 /= norm;
	  nfvtxt_ac_fvtxsfvtxn_ce01_c24b->Fill(nfvtxt,os_fvtxns_ce01_qqqq4);
	  // --- third
	  tc_os_fvtxns_ce01 = TComplex::Conjugate(tc_os_fvtxs_ce0) * tc_os_fvtxs_ce1 * TComplex::Conjugate(tc_os_fvtxn_ce0) * tc_os_fvtxn_ce1;
	  os_fvtxns_ce01_qqqq4 = tc_os_fvtxns_ce01.Re();
	  os_fvtxns_ce01_qqqq4 /= norm;
	  nfvtxt_ac_fvtxsfvtxn_ce01_c24c->Fill(nfvtxt,os_fvtxns_ce01_qqqq4);
	  // --- fourth
	  tc_os_fvtxns_ce01 = TComplex::Conjugate(tc_os_fvtxs_ce0) * TComplex::Conjugate(tc_os_fvtxs_ce1) * tc_os_fvtxn_ce0 * tc_os_fvtxn_ce1;
	  os_fvtxns_ce01_qqqq4 = tc_os_fvtxns_ce01.Re();
	  os_fvtxns_ce01_qqqq4 /= norm;
	  nfvtxt_ac_fvtxsfvtxn_ce01_c24d->Fill(nfvtxt,os_fvtxns_ce01_qqqq4);

	  // ---
	  // ---

	  float zzyzx_fvtxs_tracks_qqqq4 = calc4_event(zzyzx_fvtxs_tracks_qx2,zzyzx_fvtxs_tracks_qy2,zzyzx_fvtxs_tracks_qx4,zzyzx_fvtxs_tracks_qy4,zzyzx_fvtxs_tracks_qw);
	  float zzyzx_fvtxn_tracks_qqqq4 = calc4_event(zzyzx_fvtxn_tracks_qx2,zzyzx_fvtxn_tracks_qy2,zzyzx_fvtxn_tracks_qx4,zzyzx_fvtxn_tracks_qy4,zzyzx_fvtxn_tracks_qw);
	  float zzyzx_fvtxc_tracks_qqqq4 = calc4_event(zzyzx_fvtxc_tracks_qx2,zzyzx_fvtxc_tracks_qy2,zzyzx_fvtxc_tracks_qx4,zzyzx_fvtxc_tracks_qy4,zzyzx_fvtxc_tracks_qw);
	  nfvtxt_zzyzx_fvtxs_tracks_c24->Fill(nfvtxt,zzyzx_fvtxs_tracks_qqqq4);
	  nfvtxt_zzyzx_fvtxn_tracks_c24->Fill(nfvtxt,zzyzx_fvtxn_tracks_qqqq4);
	  nfvtxt_zzyzx_fvtxc_tracks_c24->Fill(nfvtxt,zzyzx_fvtxc_tracks_qqqq4);

	  nfvtxt_zzyzx_fvtxs_tracks_cov24->Fill(nfvtxt,zzyzx_fvtxs_tracks_qqqq4*zzyzx_fvtxs_tracks_qq2);
	  nfvtxt_zzyzx_fvtxn_tracks_cov24->Fill(nfvtxt,zzyzx_fvtxn_tracks_qqqq4*zzyzx_fvtxn_tracks_qq2);
	  nfvtxt_zzyzx_fvtxc_tracks_cov24->Fill(nfvtxt,zzyzx_fvtxc_tracks_qqqq4*zzyzx_fvtxc_tracks_qq2);

	  float zzyzx_fvtxs_tracks_six = calc6_event(tc_zzyzx_fvtxs_tracks_Q2,tc_zzyzx_fvtxs_tracks_Q4,tc_zzyzx_fvtxs_tracks_Q6,zzyzx_fvtxs_tracks_qw);
	  float zzyzx_fvtxn_tracks_six = calc6_event(tc_zzyzx_fvtxn_tracks_Q2,tc_zzyzx_fvtxn_tracks_Q4,tc_zzyzx_fvtxn_tracks_Q6,zzyzx_fvtxn_tracks_qw);
	  float zzyzx_fvtxc_tracks_six = calc6_event(tc_zzyzx_fvtxc_tracks_Q2,tc_zzyzx_fvtxc_tracks_Q4,tc_zzyzx_fvtxc_tracks_Q6,zzyzx_fvtxc_tracks_qw);
	  nfvtxt_zzyzx_fvtxs_tracks_c26->Fill(nfvtxt,zzyzx_fvtxs_tracks_six);
	  nfvtxt_zzyzx_fvtxn_tracks_c26->Fill(nfvtxt,zzyzx_fvtxn_tracks_six);
	  nfvtxt_zzyzx_fvtxc_tracks_c26->Fill(nfvtxt,zzyzx_fvtxc_tracks_six);

	}





    } //end of event

  cout << "histogram output file: " << outFile1 << endl;

  mData1->Write();
  mData1->Close();

  cout<<"cleaning up"<<endl;



  ntp_event_chain->Delete();

  // ---

  cout << "now attempting to close weight file" << endl;

  if ( phi_weight_file ) phi_weight_file->Close();

  cout<<"end of program ana"<<endl;

  return;

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
  if ( M < 2 ) return -9999;
  float numerator = Xn*Xn + Yn*Yn - M;
  float denominator = M*(M-1);
  return numerator/denominator;
}


float calccossum2_event(TComplex& Qn, TComplex& Q2n, float M)
{
  if ( M < 2 ) return -9999;
  TComplex result = Qn*Qn - Q2n;
  float numerator = result.Re();
  float denominator = M*(M-1);
  return numerator/denominator;
}

float calcsinsum2_event(TComplex& Qn, TComplex& Q2n, float M)
{
  if ( M < 2 ) return -9999;
  TComplex result = Qn*Qn - Q2n;
  float numerator = result.Im();
  float denominator = M*(M-1);
  return numerator/denominator;
}

float calccos3_event(TComplex& Qn, TComplex& Q2n, float M)
{
  if ( M < 3 ) return -9999;
  TComplex result = Qn*TComplex::Conjugate(Qn)*TComplex::Conjugate(Qn) - Qn*TComplex::Conjugate(Q2n);
  float numerator = result.Re() - 2*(M-1)*TComplex::Conjugate(Qn).Re();
  float denominator = M*(M-1)*(M-2);
  return numerator/denominator;
}

float calcsin3_event(TComplex& Qn, TComplex& Q2n, float M)
{
  if ( M < 3 ) return -9999;
  TComplex result = Qn*TComplex::Conjugate(Qn)*TComplex::Conjugate(Qn) - Qn*TComplex::Conjugate(Q2n);
  float numerator = result.Im() - 2*(M-1)*TComplex::Conjugate(Qn).Im();
  float denominator = M*(M-1)*(M-2);
  return numerator/denominator;
}

float calc4_event(float Xn, float Yn, float X2n, float Y2n, float M)
{

  if ( M < 4 ) return -9999;

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

float calc6_event(TComplex& qn, TComplex& q2n, TComplex& q3n, float M)
{

  if ( M < 6 ) return -9999;

  // TComplex qn, q2n, q3n;
  // qn = TComplex(Q2x,Q2y);
  // q2n = TComplex(Q4x,Q4y);
  // q3n = TComplex(Q6x,Q6y);

  TComplex temp1;

  // first term
  // |Qn|^6 + 9*|Q2n|^2|Qn|^2 - 6 x Re[Q2n x Qn x Qn* x Qn* x Qn*] / (Mx(M-1)x(M-2)x(M-3)x(M-4)x(M-5)
  double term1a = TMath::Power((qn*TComplex::Conjugate(qn)),3);
  double term1b = 9.0 * q2n*TComplex::Conjugate(q2n) * qn*TComplex::Conjugate(qn);
  temp1 = q2n * qn * TComplex::Conjugate(qn) * TComplex::Conjugate(qn) * TComplex::Conjugate(qn);
  double term1c = -6.0 * temp1.Re();
  double term1 = (term1a+term1b+term1c)/(M*(M-1)*(M-2)*(M-3)*(M-4)*(M-5));

  // second term
  // 4 * [Re[Q3nQn*Qn*Qn*] - 3 Re[Q3nQ2n*Qn*]] / (M(M-1)(M-2)(M-3)(M-4)(M-5)
  temp1 = q3n * TComplex::Conjugate(qn) * TComplex::Conjugate(qn) * TComplex::Conjugate(qn);
  double term2a = temp1.Re();
  temp1 = q3n * TComplex::Conjugate(q2n) * TComplex::Conjugate(qn);
  double term2b = -3.0 * temp1.Re();
  double term2 = 4.0 * (term2a+term2b)/(M*(M-1)*(M-2)*(M-3)*(M-4)*(M-5));

  // third term
  // +2 * (9*(M-4)*Re[Q2nQn*qn*] + 2 |Q3n|^2) / ((M(M-1)(M-2)(M-3)(M-4)(M-5))
  temp1 = q2n*TComplex::Conjugate(qn)*TComplex::Conjugate(qn);
  double term3a = 9.0*(M-4)*temp1.Re();
  double term3b = 2.0*q3n*TComplex::Conjugate(q3n);
  double term3 = 2.0 * (term3a + term3b) / (M*(M-1)*(M-2)*(M-3)*(M-4)*(M-5));

  // fourth term
  //double term4 = -9.0 * (TMath::Power(qn*TComplex::Conjugate(qn),2)+q2n*TComplex::Conjugate(q2n)) / (M*(M-1)*(M-2)*(M-3)*(M-5));
  double term4 = -9.0 * (TMath::Power(qn*TComplex::Conjugate(qn),2)+q2n*TComplex::Conjugate(q2n)) ;
  term4 /= (M*(M-1)*(M-2)*(M-3)*(M-5));

  // fifth term
  //double term5 = 18.0 * qn*TComplex::Conjugate(qn) / (M*(M-1)*(M-3)*(M-4));
  double term5 = 18.0 * qn*TComplex::Conjugate(qn) ;
  term5 /=  (M*(M-1)*(M-3)*(M-4));

  // sixth term
  double term6 = -6.0/((M-1)*(M-2)*(M-3));

  // cos(n(phi1+phi2+phi3-phi4-phi5-phi6))
  double six = term1 + term2 + term3 + term4 + term5 + term6;

  return (float)six; // should be smarted about this at some point

}

// -----------------------------------------------------------------
void dooffsets(int runNumber)
{

  cout << "runNumber = " << runNumber << endl;

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


  char filename[500];

  // --- get the number of files for this run number
  string pipe_out = (string) gSystem->GetFromPipe(Form("ls input/%d_*.root | grep -c r",runNumber));
  int nfiles = 0;
  nfiles = atoi(pipe_out.c_str());
  cout<<"nfiles: "<<nfiles<<endl;
  if(nfiles==0) return;

  // --- make a new TChain for the tree
  TChain *ntp_event_chain = new TChain("ntp_event");
  for ( int ifile = 0; ifile < nfiles; ++ifile )
    {
      sprintf(filename,"input/%d_%d.root",runNumber,ifile);
      cout<<"adding to tchain: "<<filename<<endl;
      ntp_event_chain->Add(filename);
    }


  char outFile1[300];
  sprintf(outFile1,"%s%d%s%d%s","output/files_",energyflag,"/coffsets_",runNumber,".root");
  cout << "histogram output file: " << outFile1 << endl;
  TFile *mData1=TFile::Open(outFile1,"recreate");
  mData1->cd();

  // ---
  // --- make histograms
  // ---

  TProfile* nfvtxt_tracks_south_qx2 = new TProfile("nfvtxt_tracks_south_qx2","",80,-0.5,79.5,-1.1,1.1,"");
  TProfile* nfvtxt_tracks_south_qx3 = new TProfile("nfvtxt_tracks_south_qx3","",80,-0.5,79.5,-1.1,1.1,"");
  TProfile* nfvtxt_tracks_south_qx4 = new TProfile("nfvtxt_tracks_south_qx4","",80,-0.5,79.5,-1.1,1.1,"");
  TProfile* nfvtxt_tracks_south_qx6 = new TProfile("nfvtxt_tracks_south_qx6","",80,-0.5,79.5,-1.1,1.1,"");
  TProfile* nfvtxt_tracks_south_qy2 = new TProfile("nfvtxt_tracks_south_qy2","",80,-0.5,79.5,-1.1,1.1,"");
  TProfile* nfvtxt_tracks_south_qy3 = new TProfile("nfvtxt_tracks_south_qy3","",80,-0.5,79.5,-1.1,1.1,"");
  TProfile* nfvtxt_tracks_south_qy4 = new TProfile("nfvtxt_tracks_south_qy4","",80,-0.5,79.5,-1.1,1.1,"");
  TProfile* nfvtxt_tracks_south_qy6 = new TProfile("nfvtxt_tracks_south_qy6","",80,-0.5,79.5,-1.1,1.1,"");

  TProfile* nfvtxt_tracks_south_inner_qx2 = new TProfile("nfvtxt_tracks_south_inner_qx2","",80,-0.5,79.5,-1.1,1.1,"");
  TProfile* nfvtxt_tracks_south_inner_qx3 = new TProfile("nfvtxt_tracks_south_inner_qx3","",80,-0.5,79.5,-1.1,1.1,"");
  TProfile* nfvtxt_tracks_south_inner_qx4 = new TProfile("nfvtxt_tracks_south_inner_qx4","",80,-0.5,79.5,-1.1,1.1,"");
  TProfile* nfvtxt_tracks_south_inner_qx6 = new TProfile("nfvtxt_tracks_south_inner_qx6","",80,-0.5,79.5,-1.1,1.1,"");
  TProfile* nfvtxt_tracks_south_inner_qy2 = new TProfile("nfvtxt_tracks_south_inner_qy2","",80,-0.5,79.5,-1.1,1.1,"");
  TProfile* nfvtxt_tracks_south_inner_qy3 = new TProfile("nfvtxt_tracks_south_inner_qy3","",80,-0.5,79.5,-1.1,1.1,"");
  TProfile* nfvtxt_tracks_south_inner_qy4 = new TProfile("nfvtxt_tracks_south_inner_qy4","",80,-0.5,79.5,-1.1,1.1,"");
  TProfile* nfvtxt_tracks_south_inner_qy6 = new TProfile("nfvtxt_tracks_south_inner_qy6","",80,-0.5,79.5,-1.1,1.1,"");

  TProfile* nfvtxt_tracks_south_outer_qx2 = new TProfile("nfvtxt_tracks_south_outer_qx2","",80,-0.5,79.5,-1.1,1.1,"");
  TProfile* nfvtxt_tracks_south_outer_qx3 = new TProfile("nfvtxt_tracks_south_outer_qx3","",80,-0.5,79.5,-1.1,1.1,"");
  TProfile* nfvtxt_tracks_south_outer_qx4 = new TProfile("nfvtxt_tracks_south_outer_qx4","",80,-0.5,79.5,-1.1,1.1,"");
  TProfile* nfvtxt_tracks_south_outer_qx6 = new TProfile("nfvtxt_tracks_south_outer_qx6","",80,-0.5,79.5,-1.1,1.1,"");
  TProfile* nfvtxt_tracks_south_outer_qy2 = new TProfile("nfvtxt_tracks_south_outer_qy2","",80,-0.5,79.5,-1.1,1.1,"");
  TProfile* nfvtxt_tracks_south_outer_qy3 = new TProfile("nfvtxt_tracks_south_outer_qy3","",80,-0.5,79.5,-1.1,1.1,"");
  TProfile* nfvtxt_tracks_south_outer_qy4 = new TProfile("nfvtxt_tracks_south_outer_qy4","",80,-0.5,79.5,-1.1,1.1,"");
  TProfile* nfvtxt_tracks_south_outer_qy6 = new TProfile("nfvtxt_tracks_south_outer_qy6","",80,-0.5,79.5,-1.1,1.1,"");

  TProfile* nfvtxt_tracks_north_qx2 = new TProfile("nfvtxt_tracks_north_qx2","",80,-0.5,79.5,-1.1,1.1,"");
  TProfile* nfvtxt_tracks_north_qx3 = new TProfile("nfvtxt_tracks_north_qx3","",80,-0.5,79.5,-1.1,1.1,"");
  TProfile* nfvtxt_tracks_north_qx4 = new TProfile("nfvtxt_tracks_north_qx4","",80,-0.5,79.5,-1.1,1.1,"");
  TProfile* nfvtxt_tracks_north_qx6 = new TProfile("nfvtxt_tracks_north_qx6","",80,-0.5,79.5,-1.1,1.1,"");
  TProfile* nfvtxt_tracks_north_qy2 = new TProfile("nfvtxt_tracks_north_qy2","",80,-0.5,79.5,-1.1,1.1,"");
  TProfile* nfvtxt_tracks_north_qy3 = new TProfile("nfvtxt_tracks_north_qy3","",80,-0.5,79.5,-1.1,1.1,"");
  TProfile* nfvtxt_tracks_north_qy4 = new TProfile("nfvtxt_tracks_north_qy4","",80,-0.5,79.5,-1.1,1.1,"");
  TProfile* nfvtxt_tracks_north_qy6 = new TProfile("nfvtxt_tracks_north_qy6","",80,-0.5,79.5,-1.1,1.1,"");

  TProfile* nfvtxt_tracks_north_inner_qx2 = new TProfile("nfvtxt_tracks_north_inner_qx2","",80,-0.5,79.5,-1.1,1.1,"");
  TProfile* nfvtxt_tracks_north_inner_qx3 = new TProfile("nfvtxt_tracks_north_inner_qx3","",80,-0.5,79.5,-1.1,1.1,"");
  TProfile* nfvtxt_tracks_north_inner_qx4 = new TProfile("nfvtxt_tracks_north_inner_qx4","",80,-0.5,79.5,-1.1,1.1,"");
  TProfile* nfvtxt_tracks_north_inner_qx6 = new TProfile("nfvtxt_tracks_north_inner_qx6","",80,-0.5,79.5,-1.1,1.1,"");
  TProfile* nfvtxt_tracks_north_inner_qy2 = new TProfile("nfvtxt_tracks_north_inner_qy2","",80,-0.5,79.5,-1.1,1.1,"");
  TProfile* nfvtxt_tracks_north_inner_qy3 = new TProfile("nfvtxt_tracks_north_inner_qy3","",80,-0.5,79.5,-1.1,1.1,"");
  TProfile* nfvtxt_tracks_north_inner_qy4 = new TProfile("nfvtxt_tracks_north_inner_qy4","",80,-0.5,79.5,-1.1,1.1,"");
  TProfile* nfvtxt_tracks_north_inner_qy6 = new TProfile("nfvtxt_tracks_north_inner_qy6","",80,-0.5,79.5,-1.1,1.1,"");

  TProfile* nfvtxt_tracks_north_outer_qx2 = new TProfile("nfvtxt_tracks_north_outer_qx2","",80,-0.5,79.5,-1.1,1.1,"");
  TProfile* nfvtxt_tracks_north_outer_qx3 = new TProfile("nfvtxt_tracks_north_outer_qx3","",80,-0.5,79.5,-1.1,1.1,"");
  TProfile* nfvtxt_tracks_north_outer_qx4 = new TProfile("nfvtxt_tracks_north_outer_qx4","",80,-0.5,79.5,-1.1,1.1,"");
  TProfile* nfvtxt_tracks_north_outer_qx6 = new TProfile("nfvtxt_tracks_north_outer_qx6","",80,-0.5,79.5,-1.1,1.1,"");
  TProfile* nfvtxt_tracks_north_outer_qy2 = new TProfile("nfvtxt_tracks_north_outer_qy2","",80,-0.5,79.5,-1.1,1.1,"");
  TProfile* nfvtxt_tracks_north_outer_qy3 = new TProfile("nfvtxt_tracks_north_outer_qy3","",80,-0.5,79.5,-1.1,1.1,"");
  TProfile* nfvtxt_tracks_north_outer_qy4 = new TProfile("nfvtxt_tracks_north_outer_qy4","",80,-0.5,79.5,-1.1,1.1,"");
  TProfile* nfvtxt_tracks_north_outer_qy6 = new TProfile("nfvtxt_tracks_north_outer_qy6","",80,-0.5,79.5,-1.1,1.1,"");

  // ---
  // --- now get the tree ready
  // ---

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

  // ---
  // --- now do event loop
  // ---

  int all_counter = 0;

  cout << "starting loop over events in the tree" << endl;
  int nentries = ntp_event_chain->GetEntries();
  cout << "total events = " << nentries << endl;
  for ( int ievt = 0; ievt < nentries; ++ievt )
    {

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

      if ( centrality > -1 )
        {
          // if ( energyflag == 200 && centrality > 5  ) continue;
          // if ( energyflag == 62  && centrality > 10 ) continue;
          // if ( energyflag == 20  && centrality > 20 ) continue;
          // if ( energyflag == 39  && centrality > 20 ) continue;
          if ( energyflag == 200 && centrality > 5  ) continue;
          if ( energyflag == 62  && centrality > 5  ) continue;
          if ( energyflag == 20  && centrality > 20 ) continue;
          if ( energyflag == 39  && centrality > 10 ) continue;
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

      // int nfvtxc = d_nFVTX_clus;
      // int nfvtxc_south = d_nFVTXN_clus;
      // int nfvtxc_north = d_nFVTXS_clus;

      // ---------------------------------------------------------------------------------------


      if ( ( say_event && verbosity > 0 ) || verbosity > 1 ) cout << "Calculating event planes" << endl;

      // --- all numbers from Darren 2016-06-23
      //const float x_off = 0.3;
      //const float beam_angle = 0.001;
      float vtx_z = d_bbcz;
      if ( eventfvtx_z > -999 ) vtx_z = eventfvtx_z;
      //float vtx_x = x_off + atan(beam_angle)*vtx_z;
      //float vtx_y = 0.02;





      // --- fvtx tracks
      float fvtxs_tracks_qx2[3]; // both, inner, outer
      float fvtxs_tracks_qy2[3];
      float fvtxs_tracks_qx3[3];
      float fvtxs_tracks_qy3[3];
      float fvtxs_tracks_qx4[3];
      float fvtxs_tracks_qy4[3];
      float fvtxs_tracks_qx6[3];
      float fvtxs_tracks_qy6[3];
      float fvtxs_tracks_qw[3];
      float fvtxn_tracks_qx2[3]; // both, inner, outer
      float fvtxn_tracks_qy2[3];
      float fvtxn_tracks_qx3[3];
      float fvtxn_tracks_qy3[3];
      float fvtxn_tracks_qx4[3];
      float fvtxn_tracks_qy4[3];
      float fvtxn_tracks_qx6[3];
      float fvtxn_tracks_qy6[3];
      float fvtxn_tracks_qw[3];

      for ( int i = 0; i < 3; ++i )
        {
          fvtxs_tracks_qx2[i] = 0.0;
          fvtxs_tracks_qy2[i] = 0.0;
          fvtxs_tracks_qx3[i] = 0.0;
          fvtxs_tracks_qy3[i] = 0.0;
          fvtxs_tracks_qx4[i] = 0.0;
          fvtxs_tracks_qy4[i] = 0.0;
          fvtxs_tracks_qx6[i] = 0.0;
          fvtxs_tracks_qy6[i] = 0.0;
          fvtxs_tracks_qw[i] = 0.0;
          fvtxn_tracks_qx2[i] = 0.0;
          fvtxn_tracks_qy2[i] = 0.0;
          fvtxn_tracks_qx3[i] = 0.0;
          fvtxn_tracks_qy3[i] = 0.0;
          fvtxn_tracks_qx4[i] = 0.0;
          fvtxn_tracks_qy4[i] = 0.0;
          fvtxn_tracks_qx6[i] = 0.0;
          fvtxn_tracks_qy6[i] = 0.0;
          fvtxn_tracks_qw[i] = 0.0;
        } // loop over layers

      int ntrack_south_inner = 0;
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
	      fvtxs_tracks_qx6[0] += fvtx_weight * TMath::Cos(6*phi);
	      fvtxs_tracks_qy6[0] += fvtx_weight * TMath::Sin(6*phi);
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
	      fvtxs_tracks_qx6[1] += fvtx_weight * TMath::Cos(6*phi);
	      fvtxs_tracks_qy6[1] += fvtx_weight * TMath::Sin(6*phi);
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
	      fvtxs_tracks_qx6[2] += fvtx_weight * TMath::Cos(6*phi);
	      fvtxs_tracks_qy6[2] += fvtx_weight * TMath::Sin(6*phi);
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
	      fvtxn_tracks_qx6[0] += fvtx_weight * TMath::Cos(6*phi);
	      fvtxn_tracks_qy6[0] += fvtx_weight * TMath::Sin(6*phi);
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
	      fvtxn_tracks_qx6[1] += fvtx_weight * TMath::Cos(6*phi);
	      fvtxn_tracks_qy6[1] += fvtx_weight * TMath::Sin(6*phi);
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
	      fvtxn_tracks_qx6[2] += fvtx_weight * TMath::Cos(6*phi);
	      fvtxn_tracks_qy6[2] += fvtx_weight * TMath::Sin(6*phi);
	      fvtxn_tracks_qw[2] += fvtx_weight;
	      ++ntrack_north_outer; // good_4_event
	    }
	} // fvtx track loop

      nfvtxt_tracks_south_qx2->Fill(nfvtxt,fvtxs_tracks_qx2[0]/fvtxs_tracks_qw[0]);
      nfvtxt_tracks_south_qx3->Fill(nfvtxt,fvtxs_tracks_qx3[0]/fvtxs_tracks_qw[0]);
      nfvtxt_tracks_south_qx4->Fill(nfvtxt,fvtxs_tracks_qx4[0]/fvtxs_tracks_qw[0]);
      nfvtxt_tracks_south_qx6->Fill(nfvtxt,fvtxs_tracks_qx6[0]/fvtxs_tracks_qw[0]);
      nfvtxt_tracks_south_qy2->Fill(nfvtxt,fvtxs_tracks_qy2[0]/fvtxs_tracks_qw[0]);
      nfvtxt_tracks_south_qy3->Fill(nfvtxt,fvtxs_tracks_qy3[0]/fvtxs_tracks_qw[0]);
      nfvtxt_tracks_south_qy4->Fill(nfvtxt,fvtxs_tracks_qy4[0]/fvtxs_tracks_qw[0]);
      nfvtxt_tracks_south_qy6->Fill(nfvtxt,fvtxs_tracks_qy6[0]/fvtxs_tracks_qw[0]);

      nfvtxt_tracks_south_inner_qx2->Fill(nfvtxt,fvtxs_tracks_qx2[1]/fvtxs_tracks_qw[1]);
      nfvtxt_tracks_south_inner_qx3->Fill(nfvtxt,fvtxs_tracks_qx3[1]/fvtxs_tracks_qw[1]);
      nfvtxt_tracks_south_inner_qx4->Fill(nfvtxt,fvtxs_tracks_qx4[1]/fvtxs_tracks_qw[1]);
      nfvtxt_tracks_south_inner_qx6->Fill(nfvtxt,fvtxs_tracks_qx6[1]/fvtxs_tracks_qw[1]);
      nfvtxt_tracks_south_inner_qy2->Fill(nfvtxt,fvtxs_tracks_qy2[1]/fvtxs_tracks_qw[1]);
      nfvtxt_tracks_south_inner_qy3->Fill(nfvtxt,fvtxs_tracks_qy3[1]/fvtxs_tracks_qw[1]);
      nfvtxt_tracks_south_inner_qy4->Fill(nfvtxt,fvtxs_tracks_qy4[1]/fvtxs_tracks_qw[1]);
      nfvtxt_tracks_south_inner_qy6->Fill(nfvtxt,fvtxs_tracks_qy6[1]/fvtxs_tracks_qw[1]);

      nfvtxt_tracks_south_outer_qx2->Fill(nfvtxt,fvtxs_tracks_qx2[2]/fvtxs_tracks_qw[2]);
      nfvtxt_tracks_south_outer_qx3->Fill(nfvtxt,fvtxs_tracks_qx3[2]/fvtxs_tracks_qw[2]);
      nfvtxt_tracks_south_outer_qx4->Fill(nfvtxt,fvtxs_tracks_qx4[2]/fvtxs_tracks_qw[2]);
      nfvtxt_tracks_south_outer_qx6->Fill(nfvtxt,fvtxs_tracks_qx6[2]/fvtxs_tracks_qw[2]);
      nfvtxt_tracks_south_outer_qy2->Fill(nfvtxt,fvtxs_tracks_qy2[2]/fvtxs_tracks_qw[2]);
      nfvtxt_tracks_south_outer_qy3->Fill(nfvtxt,fvtxs_tracks_qy3[2]/fvtxs_tracks_qw[2]);
      nfvtxt_tracks_south_outer_qy4->Fill(nfvtxt,fvtxs_tracks_qy4[2]/fvtxs_tracks_qw[2]);
      nfvtxt_tracks_south_outer_qy6->Fill(nfvtxt,fvtxs_tracks_qy6[2]/fvtxs_tracks_qw[2]);

      nfvtxt_tracks_north_qx2->Fill(nfvtxt,fvtxn_tracks_qx2[0]/fvtxn_tracks_qw[0]);
      nfvtxt_tracks_north_qx3->Fill(nfvtxt,fvtxn_tracks_qx3[0]/fvtxn_tracks_qw[0]);
      nfvtxt_tracks_north_qx4->Fill(nfvtxt,fvtxn_tracks_qx4[0]/fvtxn_tracks_qw[0]);
      nfvtxt_tracks_north_qx6->Fill(nfvtxt,fvtxn_tracks_qx6[0]/fvtxn_tracks_qw[0]);
      nfvtxt_tracks_north_qy2->Fill(nfvtxt,fvtxn_tracks_qy2[0]/fvtxn_tracks_qw[0]);
      nfvtxt_tracks_north_qy3->Fill(nfvtxt,fvtxn_tracks_qy3[0]/fvtxn_tracks_qw[0]);
      nfvtxt_tracks_north_qy4->Fill(nfvtxt,fvtxn_tracks_qy4[0]/fvtxn_tracks_qw[0]);
      nfvtxt_tracks_north_qy6->Fill(nfvtxt,fvtxn_tracks_qy6[0]/fvtxn_tracks_qw[0]);

      nfvtxt_tracks_north_inner_qx2->Fill(nfvtxt,fvtxn_tracks_qx2[1]/fvtxn_tracks_qw[1]);
      nfvtxt_tracks_north_inner_qx3->Fill(nfvtxt,fvtxn_tracks_qx3[1]/fvtxn_tracks_qw[1]);
      nfvtxt_tracks_north_inner_qx4->Fill(nfvtxt,fvtxn_tracks_qx4[1]/fvtxn_tracks_qw[1]);
      nfvtxt_tracks_north_inner_qx6->Fill(nfvtxt,fvtxn_tracks_qx6[1]/fvtxn_tracks_qw[1]);
      nfvtxt_tracks_north_inner_qy2->Fill(nfvtxt,fvtxn_tracks_qy2[1]/fvtxn_tracks_qw[1]);
      nfvtxt_tracks_north_inner_qy3->Fill(nfvtxt,fvtxn_tracks_qy3[1]/fvtxn_tracks_qw[1]);
      nfvtxt_tracks_north_inner_qy4->Fill(nfvtxt,fvtxn_tracks_qy4[1]/fvtxn_tracks_qw[1]);
      nfvtxt_tracks_north_inner_qy6->Fill(nfvtxt,fvtxn_tracks_qy6[1]/fvtxn_tracks_qw[1]);

      nfvtxt_tracks_north_outer_qx2->Fill(nfvtxt,fvtxn_tracks_qx2[2]/fvtxn_tracks_qw[2]);
      nfvtxt_tracks_north_outer_qx3->Fill(nfvtxt,fvtxn_tracks_qx3[2]/fvtxn_tracks_qw[2]);
      nfvtxt_tracks_north_outer_qx4->Fill(nfvtxt,fvtxn_tracks_qx4[2]/fvtxn_tracks_qw[2]);
      nfvtxt_tracks_north_outer_qx6->Fill(nfvtxt,fvtxn_tracks_qx6[2]/fvtxn_tracks_qw[2]);
      nfvtxt_tracks_north_outer_qy2->Fill(nfvtxt,fvtxn_tracks_qy2[2]/fvtxn_tracks_qw[2]);
      nfvtxt_tracks_north_outer_qy3->Fill(nfvtxt,fvtxn_tracks_qy3[2]/fvtxn_tracks_qw[2]);
      nfvtxt_tracks_north_outer_qy4->Fill(nfvtxt,fvtxn_tracks_qy4[2]/fvtxn_tracks_qw[2]);
      nfvtxt_tracks_north_outer_qy6->Fill(nfvtxt,fvtxn_tracks_qy6[2]/fvtxn_tracks_qw[2]);

    } // end of event loop

  cout << "histogram output file: " << outFile1 << endl;
  mData1->Write();
  mData1->Close();
  cout<<"cleaning up"<<endl;
  ntp_event_chain->Delete();
  cout<<"end of program ana"<<endl;
  return;

} // end of dooffsets



