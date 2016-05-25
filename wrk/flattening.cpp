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
//#endif

// global pmt variables
float d_pmt_x[64];
float d_pmt_y[64];
float d_pmt_z = -1443.5; // same for all tubes

//tree invariables
static const int max_nh = 500; // see from ana taxi code
//static const int max_nh = 10; // see from ana taxi code
static const int max_nf = 3000; // see from ana taxi code

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

  flatten(run,1);
  flatten(run,2);
  flatten(run,3);

  return 0;

}

// -----------------------------------------------------------------
void flatten(int runNumber, int rp_recal_pass)
{

  int verbosity = 0;

  // --- need some agreed upon standard for dealing with paths

  char outFile1[300];
  sprintf(outFile1,"%s%d%s","output/hist_",runNumber,".root"); // absolute paths need to be dealt with

  char outFile2[100];
  sprintf(outFile2,"%s%d%s","output/hrp_",runNumber,".root"); // absolute paths need to be dealt with

  cout<<"runNumber = " <<runNumber<<" "
      <<"rp_recal_pass = "<<rp_recal_pass<<endl;

  char filename[500];
  sprintf(filename,"input/tree_merged_%010d.root",runNumber); // abslute paths need to be dealt with
  //sprintf(filename,"/gpfs/mnt/gpfs02/phenix/plhf/plhf1/theok/taxi/Run15pAu200FVTXClusAna503/8833/data/434905.root"); // abslute paths need to be dealt with

  cout << "tree input file: " << filename << endl;
  cout << "v2 histogram output file: " << outFile1 << endl;

  char calibfile[500];
  sprintf(calibfile,"output/flattening_%d_%d.dat",runNumber,rp_recal_pass-1);

  TFile *f=TFile::Open( filename);

  if(!f)
    {
      cout<<"ERROR, file: "<<filename<<" could not be opened"<<endl;
      return;
    }

  TTree *htree = (TTree *)f->Get("ntp_event");

  if(!htree)
    {
      cout<<"ERROR, ntp_event could not be opened"<<endl;
      return;
    }

  cout << "Initalizing PMT positions for the BBC" << endl;

  initialize_pmt_position();

  bool fvtx_clusters = true;
  bool bbc_pmts      = true;
  bool vtx_tracks    = true;



  //int n_angle_config = 64; // move below, number is 8 * 8 (nothing to do with number of bbc tubes)
  int n_angle_config = 1; // move below, number is 8 * 8 (nothing to do with number of bbc tubes)
  // --- see below...
  int first_bbc_angle = 2;                     // 2
  int first_fvtx_angle = n_angle_config+2;     // 3
  int first_fvtx_0_angle = 2*n_angle_config+2; // 4
  int first_fvtx_1_angle = 3*n_angle_config+2; // 5
  int first_fvtx_2_angle = 4*n_angle_config+2; // 6
  int first_fvtx_3_angle = 5*n_angle_config+2; // 7


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
  float    mean[NMUL][NZPS][NHAR][NDET][2]; // mean of Psi distribution
  float    widt[NMUL][NZPS][NHAR][NDET][2]; // width of Psi distribution
  float    four[NMUL][NZPS][NHAR][NDET][2][NORD]; // ?

  // ---
  TH1D* th1d_BBC_charge = new TH1D("th1d_BBC_charge","",200,-0.5,199.5);
  TH1D* th1d_FVTX_nclus = new TH1D("th1d_FVTX_nclus","",200,-0.5,1999.5);
  TH2D* th2d_qBBC_nFVTX = new TH2D("th2d_qBBC_nFVTX","",200,-0.5,199.5,200,-0.5,1999.5);


  // --- event plane resolution
  TProfile* tp1f_reso2_BBC_CNT = new TProfile("tp1f_reso2_BBC_CNT","",1,-0.5,0.5,-1e6,1e6,"");
  TProfile* tp1f_reso2_BBC_FVTX = new TProfile("tp1f_reso2_BBC_FVTX","",1,-0.5,0.5,-1e6,1e6,"");
  TProfile* tp1f_reso2_CNT_FVTX = new TProfile("tp1f_reso2_CNT_FVTX","",1,-0.5,0.5,-1e6,1e6,"");
  TProfile* tp1f_reso3_BBC_CNT = new TProfile("tp1f_reso3_BBC_CNT","",1,-0.5,0.5,-1e6,1e6,"");
  TProfile* tp1f_reso3_BBC_FVTX = new TProfile("tp1f_reso3_BBC_FVTX","",1,-0.5,0.5,-1e6,1e6,"");
  TProfile* tp1f_reso3_CNT_FVTX = new TProfile("tp1f_reso3_CNT_FVTX","",1,-0.5,0.5,-1e6,1e6,"");

  TH1D* th1d_reso2_BBC_CNT = new TH1D("th1d_reso2_BBC_CNT","",252,-6.3,6.3);
  TH1D* th1d_reso2_BBC_FVTX = new TH1D("th1d_reso2_BBC_FVTX","",252,-6.3,6.3);
  TH1D* th1d_reso2_CNT_FVTX = new TH1D("th1d_reso2_CNT_FVTX","",252,-6.3,6.3);
  TH1D* th1d_reso3_BBC_CNT = new TH1D("th1d_reso3_BBC_CNT","",252,-6.3,6.3);
  TH1D* th1d_reso3_BBC_FVTX = new TH1D("th1d_reso3_BBC_FVTX","",252,-6.3,6.3);
  TH1D* th1d_reso3_CNT_FVTX = new TH1D("th1d_reso3_CNT_FVTX","",252,-6.3,6.3);

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



  //------------------------------------------------------------//
  //               Initializing Tree Variables                  //
  //------------------------------------------------------------//

  cout << "Now getting ready to read in the tree branch addresses and stuff...." << endl;

  //tree variables
  float        event;
  float        d_bbcz;    // bbcz
  unsigned int trigger;
  float        bc_x;
  float        bc_y;
  float        vtx_z;
  float        d_Qx[9];
  float        d_Qy[9];
  float        d_Qw[9];
  float        d_BBC_charge[64];

  int          d_nFVTX_clus;
  float        d_FVTX_x[max_nf];
  float        d_FVTX_y[max_nf];
  float        d_FVTX_z[max_nf];

  int          d_nsegments;
  float        d_px[max_nh];
  float        d_py[max_nh];
  float        d_pz[max_nh];

  TBranch *b_bbcz = htree->GetBranch("bbc_z");
  TBranch *b_event = htree->GetBranch("event");
  TBranch *b_trigger = htree->GetBranch("trigger");
  TBranch *b_bc_x = htree->GetBranch("bc_x");
  TBranch *b_bc_y = htree->GetBranch("bc_y");
  TBranch *b_vtx_z = htree->GetBranch("vtx_z");
  TBranch *b_Qx = htree->GetBranch("d_Qx");
  TBranch *b_Qy = htree->GetBranch("d_Qy");
  TBranch *b_Qw = htree->GetBranch("d_Qw");
  TBranch *b_d_BBC_charge = htree->GetBranch("d_BBC_charge");

  b_bbcz->SetAddress(&d_bbcz);
  b_event->SetAddress(&event);
  b_trigger->SetAddress(&trigger);
  b_bc_x->SetAddress(&bc_x);
  b_bc_y->SetAddress(&bc_y);
  b_vtx_z->SetAddress(&vtx_z);

  b_d_BBC_charge->SetAddress(d_BBC_charge);
  b_Qx->SetAddress(d_Qx);
  b_Qy->SetAddress(d_Qy);
  b_Qw->SetAddress(d_Qw);

  TBranch *b_d_nFVTX_clus = htree->GetBranch("d_nFVTX_clus");
  TBranch *b_d_FVTX_x = htree->GetBranch("d_FVTX_x");
  TBranch *b_d_FVTX_y = htree->GetBranch("d_FVTX_y");
  TBranch *b_d_FVTX_z = htree->GetBranch("d_FVTX_z");

  b_d_nFVTX_clus->SetAddress(&d_nFVTX_clus);
  b_d_FVTX_x->SetAddress(d_FVTX_x);
  b_d_FVTX_y->SetAddress(d_FVTX_y);
  b_d_FVTX_z->SetAddress(d_FVTX_z);

  // TBranch *b_nsegments = htree->GetBranch("nsegments");
  // TBranch *b_px = htree->GetBranch("px");
  // TBranch *b_py = htree->GetBranch("py");
  // TBranch *b_pz = htree->GetBranch("pz");
  TBranch *b_nsegments = htree->GetBranch("d_ntrk");
  TBranch *b_px = htree->GetBranch("d_cntpx");
  TBranch *b_py = htree->GetBranch("d_cntpy");
  TBranch *b_pz = htree->GetBranch("d_cntpz");

  b_nsegments->SetAddress(&d_nsegments);
  b_px->SetAddress(d_px);
  b_py->SetAddress(d_py);
  b_pz->SetAddress(d_pz);

  //------------------------------------------------------------//
  //          Finished Initializing Tree Variables              //
  //------------------------------------------------------------//

  //------------------------------------------------------------//
  //                   Looping Over Event Tree                  //
  //------------------------------------------------------------//

  cout << "starting loop over events in the tree" << endl;

  int nentries = htree->GetEntries();
  cout<<"total events = " << nentries<<endl;
  for ( int ievt = 0 ; ievt < nentries ; ievt++ )
    {

      //if ( ievt >= 1000000 ) break; // just 1M events for now, runs a little on the slow side...

      bool say_event = ( ievt%1000==0 );

      if ( say_event ) cout<<"event number = "<<ievt<<endl;
      //if(ievt ==1000) break;

      if ( ( say_event && verbosity > 0 ) || verbosity > 1 ) cout << "getting event level variables" << endl;

      b_bbcz->GetEntry(ievt);
      b_event->GetEntry(ievt);
      b_trigger->GetEntry(ievt);
      b_bc_x->GetEntry(ievt);
      b_bc_y->GetEntry(ievt);
      b_bbcz->GetEntry(ievt);
      b_vtx_z->GetEntry(ievt);
      //standard q vector
      b_Qx->GetEntry(ievt);
      b_Qy->GetEntry(ievt);
      b_Qw->GetEntry(ievt);

      if ( ( say_event && verbosity > 0 ) || verbosity > 1 ) cout << "getting BBC PMT information" << endl;

      if ( bbc_pmts )  b_d_BBC_charge->GetEntry(ievt);//bbc pmts charge

      if ( ( say_event && verbosity > 0 ) || verbosity > 1 ) cout << "getting FVTX cluster information" << endl;

      //FVTX clusters
      if ( fvtx_clusters )
        {
          b_d_nFVTX_clus->GetEntry(ievt);
          b_d_FVTX_x->GetEntry(ievt);
          b_d_FVTX_y->GetEntry(ievt);
          b_d_FVTX_z->GetEntry(ievt);
        }

      if ( ( say_event && verbosity > 0 ) || verbosity > 1 ) cout << "getting track information" << endl;

      //VTX Tracks
      if ( vtx_tracks )
        {
          b_nsegments->GetEntry(ievt);
          b_px->GetEntry(ievt);
          b_py->GetEntry(ievt);
          b_pz->GetEntry(ievt);
        }

      if ( ( say_event && verbosity > 0 ) || verbosity > 1 ) cout << "Finished getting tree variables" << endl;

      //if(vtx_z!=vtx_z) continue;
      //if( fabs(vtx_z) > 100) continue;
      // --- some big questions here about the z-vertex cut
      // --- really need to double and triple check where the zvertex cuts are applied and what they are
      //int ibbcz  = NZPS*(d_bbcz+30)/60;//bbc z bin for -30 <bbc z < 30 // how do you specify the number of bins here
      if(TMath::Abs(d_bbcz) > 10.0 ) continue;
      int ibbcz  = NZPS*(d_bbcz+10)/20;//bbc z bin for -10 <bbc z < 10 // how do you specify the number of bins here

      // --- break and continue statements should happen much, much earlier --------------------
      if(rp_recal_pass<1 || rp_recal_pass > 3) break;// rp_recal_pass only valid between 1 and 3

      // make sure bin number doesn't exceed number of bins
      if ( ibbcz<0 || ibbcz >= NZPS ) continue;

      // don't do analysis if no tracks on third pass
      //if(d_nsegments==0 && rp_recal_pass > 2) continue;
      // ---------------------------------------------------------------------------------------

      //------------------------------------------------------------//
      //                Calculating Event Planes                    //
      //------------------------------------------------------------//

      if ( ( say_event && verbosity > 0 ) || verbosity > 1 ) cout << "Calculating event planes" << endl;

      //float vtx_x = bc_x + d_bbcz*0.025/10; // z dependent x position because of beam angle rotation issues // what are these numbers?  // these are really specific to p+Au
      // --- just use the beam center for now
      // --- beam center not available yet!!!
      // float vtx_x = bc_x;
      // float vtx_y = bc_y;
      float vtx_x = 0;
      float vtx_y = 0;

      float bbc_qx2 = 0;
      float bbc_qy2 = 0;
      float bbc_qx3 = 0;
      float bbc_qy3 = 0;
      float bbc_qw = 0;

      if ( ( say_event && verbosity > 0 ) || verbosity > 1 ) cout << "Looping over BBC stuff now" << endl;

      if ( bbc_pmts )
        {
          for(int ipmt = 0; ipmt < 64; ipmt++)
            {
              float bbc_x      = d_pmt_x[ipmt] - vtx_x*10;//pmt location in mm
              float bbc_y      = d_pmt_y[ipmt] - vtx_y*10;
              float bbc_z      = d_pmt_z       - d_bbcz*10;
              float bbc_charge = d_BBC_charge[ipmt];
              if(bbc_charge <= 0) continue;
              //cout<<"bbc_x: "<<bbc_x<<" bbc_y: "<<bbc_y<<" bbc_z: "<<bbc_z<<endl;

              double bbc_r = sqrt(pow(bbc_x,2.0)+pow(bbc_y,2.0));

              double bbc_the = atan2(bbc_r,bbc_z); //fvtx_z-bbcv
              double bbc_eta = -log(tan(0.5*bbc_the));

              float phi = TMath::ATan2(bbc_y,bbc_x);

              float mass = 0.1396;//assume charged pion mass
              float pT = 0.25;
              float px = pT * TMath::Cos(phi);
              float py = pT * TMath::Sin(phi);
              float pz = pT * TMath::SinH(bbc_eta);

              float energy = TMath::Sqrt(px*px+py*py+pz*pz+mass*mass);
              TLorentzVector particle_vec(px,py,pz,energy);
              phi = TMath::ATan2(particle_vec.Py(),particle_vec.Px());

              // --- need to add a harmonic loop here, this will only give psi for now
              bbc_qx2 += bbc_charge*TMath::Cos(2*phi);
              bbc_qy2 += bbc_charge*TMath::Sin(2*phi);
              bbc_qx3 += bbc_charge*TMath::Cos(3*phi);
              bbc_qy3 += bbc_charge*TMath::Sin(3*phi);
              bbc_qw += bbc_charge;
            } // loop over tubes
        }

      th1d_BBC_charge->Fill(bbc_qw);
      th1d_FVTX_nclus->Fill(d_nFVTX_clus);
      th2d_qBBC_nFVTX->Fill(bbc_qw,d_nFVTX_clus);

      //continue; // testing to get the charge distribution to make a centrality selection

      // --- do centrality cut here!!!
      //if ( bbc_qw < 61.5 ) continue; // dAu 200 GeV
      if ( bbc_qw < 30.0 ) continue; // very rough dAu 62 GeV

      float fvtx_qx2[5];//all layers then 0 1 2 3
      float fvtx_qy2[5];
      float fvtx_qx3[5];//all layers then 0 1 2 3
      float fvtx_qy3[5];
      float fvtx_qw[5];

      for(int ilayer = 0; ilayer < 5; ilayer++)
        {
          fvtx_qx2[ilayer] = 0.0;
          fvtx_qy2[ilayer] = 0.0;
          fvtx_qx3[ilayer] = 0.0;
          fvtx_qy3[ilayer] = 0.0;
          fvtx_qw[ilayer] = 0.0;
        } // loop over layers


      if ( ( say_event && verbosity > 0 ) || verbosity > 1 ) cout << "Looping over FVTX cluster" << endl;
      if ( fvtx_clusters )
        {
          for(int iclus = 0; iclus < d_nFVTX_clus; iclus++)
            {
              float fvtx_x      = d_FVTX_x[iclus] - vtx_x; // calculate for each event, function of z
              float fvtx_y      = d_FVTX_y[iclus] - vtx_y;
              float fvtx_z      = d_FVTX_z[iclus]; // need raw z to get layer

              double fvtx_r = sqrt(pow(fvtx_x,2.0)+pow(fvtx_y,2.0));

              double fvtx_the = atan2(fvtx_r,fvtx_z - d_bbcz); //fvtx_z-bbcv // add a new variable to make it clear that you're using a corrected z vertex
              double fvtx_eta = -log(tan(0.5*fvtx_the));

              // --- probaby don't need to cut on acceptance
              // --- shouldn't get any clusters outside of acceptance
              if(!(fabs(fvtx_eta)>1.0 && fabs(fvtx_eta)<3.5)) continue;
              // cout<<"fvtx_x: "<<fvtx_x<<" fvtx_y: "<<fvtx_y<<" fvtx_z"<<fvtx_z<<endl;

              int fvtx_layer    = get_fvtx_layer(fvtx_z); // raw z to get layer
              if(fvtx_layer < 0) continue;
              // --- gap cut, not sure what this does
              int igap = (fabs(fvtx_eta)-1.0)/0.5;

              int id_fvtx = fvtx_layer*5+igap;

              if(!(id_fvtx>=0 && id_fvtx<40)) continue;
              // --------------------------------------

              float phi = TMath::ATan2(fvtx_y,fvtx_x);

              float mass = 0.1396;//assume charged pion mass
              float pT = 0.25;
              float px = pT * TMath::Cos(phi);
              float py = pT * TMath::Sin(phi);
              float pz = pT * TMath::SinH(fvtx_eta);

              float energy = TMath::Sqrt(px*px+py*py+pz*pz+mass*mass);
              TLorentzVector particle_vec(px,py,pz,energy);


              phi = TMath::ATan2(particle_vec.Py(),particle_vec.Px());

              fvtx_qx2[fvtx_layer+1] += TMath::Cos(2*phi);
              fvtx_qy2[fvtx_layer+1] += TMath::Sin(2*phi);
              fvtx_qx3[fvtx_layer+1] += TMath::Cos(3*phi);
              fvtx_qy3[fvtx_layer+1] += TMath::Sin(3*phi);

              fvtx_qx2[0] += TMath::Cos(2*phi);
              fvtx_qy2[0] += TMath::Sin(2*phi);
              fvtx_qx3[0] += TMath::Cos(3*phi);
              fvtx_qy3[0] += TMath::Sin(3*phi);

              fvtx_qw[fvtx_layer+1] ++;
              fvtx_qw[0] ++;
            } // loop over cluster
        } // check on clusters



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
          sumxy[1][first_bbc_angle][0] = bbc_qx2;
          sumxy[1][first_bbc_angle][1] = bbc_qy2;
          sumxy[1][first_bbc_angle][2] = bbc_qw;
          sumxy[2][first_bbc_angle][0] = bbc_qx3;
          sumxy[2][first_bbc_angle][1] = bbc_qy3;
          sumxy[2][first_bbc_angle][2] = bbc_qw;
        }

      if ( fvtx_clusters )
        {
          sumxy[1][first_fvtx_angle][0] = fvtx_qx2[0];
          sumxy[1][first_fvtx_angle][1] = fvtx_qy2[0];
          sumxy[1][first_fvtx_angle][2] = fvtx_qw[0];

          sumxy[1][first_fvtx_0_angle][0] = fvtx_qx2[1];
          sumxy[1][first_fvtx_0_angle][1] = fvtx_qy2[1];
          sumxy[1][first_fvtx_0_angle][2] = fvtx_qw[1];

          sumxy[1][first_fvtx_1_angle][0] = fvtx_qx2[2];
          sumxy[1][first_fvtx_1_angle][1] = fvtx_qy2[2];
          sumxy[1][first_fvtx_1_angle][2] = fvtx_qw[2];

          sumxy[1][first_fvtx_2_angle][0] = fvtx_qx2[3];
          sumxy[1][first_fvtx_2_angle][1] = fvtx_qy2[3];
          sumxy[1][first_fvtx_2_angle][2] = fvtx_qw[3];

          sumxy[1][first_fvtx_3_angle][0] = fvtx_qx2[4];
          sumxy[1][first_fvtx_3_angle][1] = fvtx_qy2[4];
          sumxy[1][first_fvtx_3_angle][2] = fvtx_qw[4];

          // ---

          sumxy[2][first_fvtx_angle][0] = fvtx_qx3[0];
          sumxy[2][first_fvtx_angle][1] = fvtx_qy3[0];
          sumxy[2][first_fvtx_angle][2] = fvtx_qw[0];

          sumxy[2][first_fvtx_0_angle][0] = fvtx_qx3[1];
          sumxy[2][first_fvtx_0_angle][1] = fvtx_qy3[1];
          sumxy[2][first_fvtx_0_angle][2] = fvtx_qw[1];

          sumxy[2][first_fvtx_1_angle][0] = fvtx_qx3[2];
          sumxy[2][first_fvtx_1_angle][1] = fvtx_qy3[2];
          sumxy[2][first_fvtx_1_angle][2] = fvtx_qw[2];

          sumxy[2][first_fvtx_2_angle][0] = fvtx_qx3[3];
          sumxy[2][first_fvtx_2_angle][1] = fvtx_qy3[3];
          sumxy[2][first_fvtx_2_angle][2] = fvtx_qw[3];

          sumxy[2][first_fvtx_3_angle][0] = fvtx_qx3[4];
          sumxy[2][first_fvtx_3_angle][1] = fvtx_qy3[4];
          sumxy[2][first_fvtx_3_angle][2] = fvtx_qw[4];
        }


      if ( DIAG )
        {
          cout<<"bbc from node tree: "<<d_Qx[5]<<" "<<d_Qy[5]<<" "<<d_Qw[5]<<endl;
          cout<<"bbc from me: "<<bbc_qx2<<" "<<bbc_qy2<<" "<<bbc_qw<<endl;

          cout<<"fvtx raw: "<<endl;
          cout<<"from node tree: "<<d_Qx[4]<<" "<<d_Qy[4]<<" "<<d_Qw[4]<<endl;
          cout<<"from clusters: " <<fvtx_qx2[0]<<" "<<fvtx_qy2[0]<<" "<<fvtx_qw[0]<<endl;
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
                  psi_bf[ic][ih][id]->Fill(ibbcz,psi);
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
                  if (rp_recal_pass>0) dis[icent][ih][id]->Fill(ibbcz,sumxy[ih][id][3]*(ih+1.0));
                }
              if (sumxy[ih][id][2]>0.0) // check on weight (x,y,w,psi)
                {
                  for (int ib=0; ib<2; ib++)
                    {
                      sumxy[ih][id][ib]/=sumxy[ih][id][2]; // normalize to the weight
                      //if(ih==1 && id==0 && ib==0 && sumxy[ih][id][ib]>1) cout<<sumxy[ih][id][ib]<<endl;
                      if (rp_recal_pass>0)
                        {
                          ave[icent][ibbcz][ih][id]->Fill(ib+0.0,sumxy[ih][id][ib]);
                          if(id==0 && DIAG) cout<<"filled ave: "<<ih<<" "<<id<<" "<<ib<<" with: "<<sumxy[ih][id][ib]<<endl;
                          if(ib==0) qx[icent][ih][id]->Fill(ibbcz,sumxy[ih][id][0]);
                          if(ib==1) qy[icent][ih][id]->Fill(ibbcz,sumxy[ih][id][1]);
                        } // pass > 0
                      float sxy=sumxy[ih][id][ib];
                      float mxy=mean[icent][ibbcz][ih][id][ib]; // for recentering qx and qy
                      float wxy=widt[icent][ibbcz][ih][id][ib]; // for recentering qx and qy

                      //if(icent==0 && ibbcz==0 && ih==1 && id==0) cout<<ib<<" "<<sxy<<" "<<mxy<<" "<<wxy<<endl;
                      sumxy[ih][id][ib]=(sxy-mxy)/wxy; // recentered by mean and renormalized to width
                      if (rp_recal_pass>0)
                        {
                          ave[icent][ibbcz][ih][id]->Fill(ib+2.0,sumxy[ih][id][ib]);  // ib+2 to avoid overlap
                          if(id==0 && DIAG) cout<<"filled ave2: "<<ih<<" "<<id<<" "<<ib<<" with: "<<sumxy[ih][id][ib]<<endl;
                          if(ib==0) qx[icent][ih][id]->Fill(ibbcz+NZPS,sumxy[ih][id][0]);
                          if(ib==1) qy[icent][ih][id]->Fill(ibbcz+NZPS,sumxy[ih][id][1]);
                        } // pass > 0
                    } // if weight > 0

                  sumxy[ih][id][3]=atan2(sumxy[ih][id][1],sumxy[ih][id][0])/(ih+1.0);
                  if (rp_recal_pass>0)
                    {
                      // fill histogram with psi calculated with recenter q vectors
                      dis[icent][ih][id]->Fill(ibbcz+NZPS,sumxy[ih][id][3]*(ih+1.0));
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
                      if (rp_recal_pass>0) flt[icent][ibbcz][ih][id]->Fill(io+0.0,cc);
                      if (rp_recal_pass>0) flt[icent][ibbcz][ih][id]->Fill(io+NORD,ss);
                      // --- four means fourier
                      float aa=four[icent][ibbcz][ih][id][0][io]; // mean cos
                      float bb=four[icent][ibbcz][ih][id][1][io]; // mean sin
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
                      if (rp_recal_pass>0) flt[icent][ibbcz][ih][id]->Fill(io+NORD*2.0,cc);
                      if (rp_recal_pass>0) flt[icent][ibbcz][ih][id]->Fill(io+NORD*3.0,ss);
                    }
                  sumxy[ih][id][3]=psi/(ih+1.0);
                  if (rp_recal_pass>0) dis[icent][ih][id]->Fill(ibbcz+NZPS*2.0,sumxy[ih][id][3]*(ih+1.0));
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
                  psi_af[ic][ih][id]->Fill(ibbcz,sumxy[ih][id][3]);
                  if ( DIAG ) cout<<"CORR: for id: "<<id<<" psi: "<<sumxy[ih][id][3]<<endl;
                }
            }
        }

      // ---
      // --- now going to calculate v2
      // ---
      if(rp_recal_pass<3) continue; // don't calculate v2 except for final pass

      // --- looks like these are already done above
      float bbc_psi2 = (sumxy[1][0][2]>0)?sumxy[1][0][3]:-9999.9;
      // --- maybe 12 should be a different number (we have 4 above)
      float fvtx_psi2 = (sumxy[1][1][2]>12)?sumxy[1][1][3]:-9999.9;

      float bbc_psi2_docalib;
      float fvtx_psi2_docalib;
      float fvtx0_psi2_docalib;
      float fvtx1_psi2_docalib;
      float fvtx2_psi2_docalib;
      float fvtx3_psi2_docalib;

      float bbc_psi3_docalib;
      float fvtx_psi3_docalib;
      float fvtx0_psi3_docalib;
      float fvtx1_psi3_docalib;
      float fvtx2_psi3_docalib;
      float fvtx3_psi3_docalib;

      if ( bbc_pmts )
        {
          bbc_psi2_docalib = (sumxy[1][first_bbc_angle][2]>0)?sumxy[1][first_bbc_angle][3]:-9999.9;
          bbc_psi3_docalib = (sumxy[2][first_bbc_angle][2]>0)?sumxy[2][first_bbc_angle][3]:-9999.9;
        }
      if ( fvtx_clusters )
        {
          fvtx_psi2_docalib = (sumxy[1][first_fvtx_angle][2]>4)?sumxy[1][first_fvtx_angle][3]:-9999.9;
          fvtx0_psi2_docalib = (sumxy[1][first_fvtx_0_angle][2]>4)?sumxy[1][first_fvtx_0_angle][3]:-9999.9;
          fvtx1_psi2_docalib = (sumxy[1][first_fvtx_1_angle][2]>4)?sumxy[1][first_fvtx_1_angle][3]:-9999.9;
          fvtx2_psi2_docalib = (sumxy[1][first_fvtx_2_angle][2]>4)?sumxy[1][first_fvtx_2_angle][3]:-9999.9;
          fvtx3_psi2_docalib = (sumxy[1][first_fvtx_3_angle][2]>4)?sumxy[1][first_fvtx_3_angle][3]:-9999.9;
          fvtx_psi3_docalib = (sumxy[2][first_fvtx_angle][2]>4)?sumxy[2][first_fvtx_angle][3]:-9999.9;
          fvtx0_psi3_docalib = (sumxy[2][first_fvtx_0_angle][2]>4)?sumxy[2][first_fvtx_0_angle][3]:-9999.9;
          fvtx1_psi3_docalib = (sumxy[2][first_fvtx_1_angle][2]>4)?sumxy[2][first_fvtx_1_angle][3]:-9999.9;
          fvtx2_psi3_docalib = (sumxy[2][first_fvtx_2_angle][2]>4)?sumxy[2][first_fvtx_2_angle][3]:-9999.9;
          fvtx3_psi3_docalib = (sumxy[2][first_fvtx_3_angle][2]>4)?sumxy[2][first_fvtx_3_angle][3]:-9999.9;
        }

      tp1f_reso2_BBC_FVTX->Fill(0.0,cos(2*(bbc_psi2_docalib-fvtx_psi2_docalib)));
      tp1f_reso3_BBC_FVTX->Fill(0.0,cos(3*(bbc_psi3_docalib-fvtx_psi3_docalib)));

      th1d_reso2_BBC_FVTX->Fill(bbc_psi2_docalib-fvtx_psi2_docalib);
      th1d_reso3_BBC_FVTX->Fill(bbc_psi3_docalib-fvtx_psi3_docalib);

      //start of vtx stand alone track loop
      if ( vtx_tracks )
        {
          for(int itrk=0; itrk< d_nsegments; itrk++)
            {
              float px    = d_px[itrk];
              float py    = d_py[itrk];
              float pz    = d_pz[itrk];
              int dcarm=0;
              if(px>0) dcarm=1;

              float phi0 = TMath::ATan2(py,px);
              float pt = sqrt(px*px+py*py);

              float bbc_dphi2 = phi0 - bbc_psi2; // move this lower?

              if(-4.0<bbc_psi2 && bbc_psi2<4.0 ) // why this weird cut? why not just -pi to pi? // checking against -9999 from above
                {
                  bbcs_v2_incl_nodetree->Fill(pt,cos(2*bbc_dphi2));
                  if(dcarm==0)
                    {
                      bbcs_v2_east_nodetree->Fill(pt,cos(2*bbc_dphi2));
                    }
                  else if(dcarm==1)
                    {
                      bbcs_v2_west_nodetree->Fill(pt,cos(2*bbc_dphi2));
                    }
                }

              float fvtx_dphi2 = phi0 - fvtx_psi2;

              if(-4.0<fvtx_psi2 && fvtx_psi2<4.0 )
                {
                  fvtxs_v2_incl_nodetree->Fill(pt,cos(2*fvtx_dphi2));
                  if(dcarm==0)
                    {
                      fvtxs_v2_east_nodetree->Fill(pt,cos(2*fvtx_dphi2));
                    }
                  else if(dcarm==1)
                    {
                      fvtxs_v2_west_nodetree->Fill(pt,cos(2*fvtx_dphi2));
                    }
                }

              // --- finished with nodetree part, now doing docalib part

              float mass = 0.1396;//assume charged pion mass
              float energy = TMath::Sqrt(px*px+py*py+pz*pz+mass*mass);
              TLorentzVector particle_vec(px,py,pz,energy);

              float phi_angle = TMath::ATan2(particle_vec.Py(),particle_vec.Px()); // rotated phi // "modified"
              float pt_angle = TMath::Sqrt(particle_vec.Py()*particle_vec.Py()+particle_vec.Px()*particle_vec.Px()); // rotated pt "modified"

              //bbc angle
              if ( bbc_pmts )
                {
                  if(-4.0<bbc_psi2_docalib && bbc_psi2_docalib<4.0)
                    {
                      // --- 2nd harmonic
                      double bbc_dphi2_docalib = phi_angle - bbc_psi2_docalib;
                      double cosbbc_dphi2_docalib = TMath::Cos(2*bbc_dphi2_docalib);
                      bbcs_v2_both_docalib->Fill(pt_angle,cosbbc_dphi2_docalib);
                      if ( dcarm == 1 ) bbcs_v2_west_docalib->Fill(pt_angle,cosbbc_dphi2_docalib);
                      if ( dcarm == 0 ) bbcs_v2_east_docalib->Fill(pt_angle,cosbbc_dphi2_docalib);
                      // --- 3rd harmonic
                      double bbc_dphi3_docalib = phi_angle - bbc_psi3_docalib;
                      double cosbbc_dphi3_docalib = TMath::Cos(3*bbc_dphi3_docalib);
                      bbcs_v3_both_docalib->Fill(pt_angle,cosbbc_dphi3_docalib);
                      if ( dcarm == 1 ) bbcs_v3_west_docalib->Fill(pt_angle,cosbbc_dphi3_docalib);
                      if ( dcarm == 0 ) bbcs_v3_east_docalib->Fill(pt_angle,cosbbc_dphi3_docalib);
                      // --- ep reso
                      tp1f_reso2_BBC_CNT->Fill(0.0,cosbbc_dphi2_docalib);
                      tp1f_reso3_BBC_CNT->Fill(0.0,cosbbc_dphi3_docalib);
                      th1d_reso2_BBC_CNT->Fill(bbc_dphi2_docalib);
                      th1d_reso3_BBC_CNT->Fill(bbc_dphi3_docalib);
                    }
                } // check on tubes

              if(!fvtx_clusters) continue;
              //fvtx all layers
              if(-4.0<fvtx_psi2_docalib && fvtx_psi2_docalib<4.0)
                {
                  // --- 2nd harmonic
                  double fvtx_dphi2_docalib = phi_angle - fvtx_psi2_docalib;
                  double cosfvtx_dphi2_docalib = TMath::Cos(2*fvtx_dphi2_docalib);
                  fvtxs_v2_both_docalib->Fill(pt_angle,cosfvtx_dphi2_docalib);
                  if ( dcarm == 1 ) fvtxs_v2_west_docalib->Fill(pt_angle,cosfvtx_dphi2_docalib);
                  if ( dcarm == 0 ) fvtxs_v2_east_docalib->Fill(pt_angle,cosfvtx_dphi2_docalib);
                  // --- 3rd harmonic
                  double fvtx_dphi3_docalib = phi_angle - fvtx_psi3_docalib;
                  double cosfvtx_dphi3_docalib = TMath::Cos(3*fvtx_dphi3_docalib);
                  fvtxs_v3_both_docalib->Fill(pt_angle,cosfvtx_dphi3_docalib);
                  if ( dcarm == 1 ) fvtxs_v3_west_docalib->Fill(pt_angle,cosfvtx_dphi3_docalib);
                  if ( dcarm == 0 ) fvtxs_v3_east_docalib->Fill(pt_angle,cosfvtx_dphi3_docalib);
                  // --- ep reso
                  tp1f_reso2_CNT_FVTX->Fill(0.0,cosfvtx_dphi2_docalib);
                  tp1f_reso3_CNT_FVTX->Fill(0.0,cosfvtx_dphi3_docalib);
                  th1d_reso2_CNT_FVTX->Fill(fvtx_dphi2_docalib);
                  th1d_reso3_CNT_FVTX->Fill(fvtx_dphi3_docalib);
                }

              // --- now fvtx layers

              //fvtx layer 0
              if(-4.0<fvtx0_psi2_docalib && fvtx0_psi2_docalib<4.0)
                {
                  double fvtx_dphi2_docalib = phi_angle - fvtx0_psi2_docalib;
                  double cosfvtx_dphi2_docalib = TMath::Cos(2*fvtx_dphi2_docalib);

                  if(dcarm==1)
                    {
                      fvtxs0_v2_west_docalib->Fill(pt_angle,cosfvtx_dphi2_docalib);
                    }
                  else if( dcarm==0)
                    {
                      fvtxs0_v2_east_docalib->Fill(pt_angle,cosfvtx_dphi2_docalib);
                    }
                  fvtxs0_v2_both_docalib->Fill(pt_angle,cosfvtx_dphi2_docalib);
                }

              //fvtx layer 1
              if(-4.0<fvtx1_psi2_docalib && fvtx1_psi2_docalib<4.0)
                {
                  double fvtx_dphi2_docalib = phi_angle - fvtx1_psi2_docalib;
                  double cosfvtx_dphi2_docalib = TMath::Cos(2*fvtx_dphi2_docalib);

                  if(dcarm==1)
                    {
                      fvtxs1_v2_west_docalib->Fill(pt_angle,cosfvtx_dphi2_docalib);
                    }
                  else if( dcarm==0)
                    {
                      fvtxs1_v2_east_docalib->Fill(pt_angle,cosfvtx_dphi2_docalib);
                    }
                  fvtxs1_v2_both_docalib->Fill(pt_angle,cosfvtx_dphi2_docalib);
                }

              //fvtx layer 2
              if(-4.0<fvtx2_psi2_docalib && fvtx2_psi2_docalib<4.0)
                {
                  double fvtx_dphi2_docalib = phi_angle - fvtx2_psi2_docalib;
                  double cosfvtx_dphi2_docalib = TMath::Cos(2*fvtx_dphi2_docalib);

                  if(dcarm==1)
                    {
                      fvtxs2_v2_west_docalib->Fill(pt_angle,cosfvtx_dphi2_docalib);
                    }
                  else if( dcarm==0)
                    {
                      fvtxs2_v2_east_docalib->Fill(pt_angle,cosfvtx_dphi2_docalib);
                    }
                  fvtxs2_v2_both_docalib->Fill(pt_angle,cosfvtx_dphi2_docalib);
                }

              //fvtx layer 3
              if(-4.0<fvtx3_psi2_docalib && fvtx3_psi2_docalib<4.0)
                {
                  double fvtx_dphi2_docalib = phi_angle - fvtx3_psi2_docalib;
                  double cosfvtx_dphi2_docalib = TMath::Cos(2*fvtx_dphi2_docalib);

                  if(dcarm==1)
                    {
                      fvtxs3_v2_west_docalib->Fill(pt_angle,cosfvtx_dphi2_docalib);
                    }
                  else if( dcarm==0)
                    {
                      fvtxs3_v2_east_docalib->Fill(pt_angle,cosfvtx_dphi2_docalib);
                    }
                  fvtxs3_v2_both_docalib->Fill(pt_angle,cosfvtx_dphi2_docalib);
                } // check on psi
            } // loop over tracks
        } // check on tracks
    }//end of event



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
      th1d_BBC_charge->Write();
      th1d_FVTX_nclus->Write();
      th2d_qBBC_nFVTX->Write();
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

      th1d_BBC_charge->Write();
      th1d_FVTX_nclus->Write();
      th2d_qBBC_nFVTX->Write();

      mData1->Close();

    } // check on last pass


  //return;

  cout<<"cleaning up"<<endl;

  htree->Delete();
  f->Close();
  delete f;
  /*
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
    }
    }
    }
  */
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



int get_fvtx_layer(float z)
{
  if(z < -18 && z > -24) return 0;
  if(z < -24 && z > -30) return 1;
  if(z < -30 && z > -35) return 2;
  if(z < -35)            return 3;

  cout<<"get_fvtx_layer::invalid z =  "<<z<<endl;
  return -1;
}

