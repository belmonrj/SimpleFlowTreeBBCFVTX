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
static const int max_nh = 10; // see from ana taxi code
static const int max_nf = 750; // see from ana taxi code

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



// Original code sent by Theo K. on May 05, 2016
// J. Nagle - modifications (comments and diagnostics added - see tag "// JLN"

bool DIAG = false;

// JLN - 1) this function takes as 2nd, 3rd arguments the beam vectors in the Lab Frame
//       2) it then determines the needed boost and rotatation (in that order) into the center-of-mass frame 
//          with the beams along the z-axis
//       3) it then uses these boost and rotation values to either change the vector "vec" 
//          from the Lab to the CoM frame (option 1) or from the CoM fram to the Lab frame (option 2)
//       4) JLN added an option flag for this function
//       5) Note that in this original version the function also returns the relative phi angle between frames... 
float boost_and_rotate(TLorentzVector & vec, TLorentzVector input1, TLorentzVector input2, int option)
{

  // JLN - just define a special vector for the z-axis
  TVector3 z(0,0,1);

  if (DIAG) {
    cout << "===================================================" << endl;
    cout << "boost_and_rotate:   input beam vectors (1,2) and then particle of interest" << endl;
    cout << "px1 = " << input1.Px() << ", py1 = " << input1.Py() << ", pz1 = " << input1.Pz() << ", E = " << input1.E() << endl;
    cout << "px2 = " << input2.Px() << ", py2 = " << input2.Py() << ", pz2 = " << input2.Pz() << ", E = " << input2.E() << endl;
    cout << "px  = " << vec.Px()    << ", py  = " << vec.Py()    << ", pz  = " << vec.Pz()    << ", E = " << vec.E()    << endl;
    cout << " " << endl;
  }

  // calculate the total four vector and then boost both beam vectors into CoM frame (not yet along z-axis)
  TLorentzVector cms = input1 + input2;
  input1.Boost(-cms.BoostVector()); // blue
  input2.Boost(-cms.BoostVector()); // yellow


  if (DIAG) {
    cout << "boost_and_rotate:   input beam vectors after BOOST TO CoM ONLY" << endl;
    cout << "px1 = " << input1.Px() << ", py1 = " << input1.Py() << ", pz1 = " << input1.Pz() << ", E = " << input1.E() << endl;
    cout << "px2 = " << input2.Px() << ", py2 = " << input2.Py() << ", pz2 = " << input2.Pz() << ", E = " << input2.E() << endl;
    cout << "Boost pxcms = " << cms.Px() << ", pycms = " << cms.Py() << ", pzcms = " << cms.Pz() << ", E = " << cms.E() << endl;
    // JLN - as a check, boost cms vector into cvs frame (should have identically zero momentum)
    // problem later on !!!!!
    //    cms.Boost(-cms.BoostVector());
    //    cout << "Check pxcms = " << cms.Px() << ", pycms = " << cms.Py() << ", pzcms = " << cms.Pz() << ", E = " << cms.E() << endl;
    cout << " " << endl;
  }

  // Now rotate the beams about x to align them with z axis
  double rotAngleY1 = -input1.Angle(z);
  double rotAngleY2 = -input2.Angle(z);   // this assumes rotation is just interchanging x and z !!!!
 
  if (DIAG) {
    cout << "angle around Y = " << rotAngleY1 << " (and alternate check = " << rotAngleY2 << " )" << endl;
  }

  // JLN - assuming check works...
  double rotAngleY = rotAngleY2;
  // JLN - current check by Theo about whether rotation is in wrong direction....
  input1.RotateY(rotAngleY);
  if(TMath::Abs(input1.Px()) >  0.00001 ) {
    // wrong rotation direction
    // rotate back and set the new angle
    rotAngleY = -1.0 * rotAngleY;
    input1.RotateY(rotAngleY);
  }
  // JLN - do the correct rotation of the beam vectors
  input1.RotateY(rotAngleY);
  input2.RotateY(rotAngleY);

  if (DIAG) {
    cout << "boost_and_rotate:   input beam vectors after BOOST TO CoM and ROTATION TOO" << endl;
    cout << "px1 = " << input1.Px() << ", py1 = " << input1.Py() << ", pz1 = " << input1.Pz() << ", E = " << input1.E() << endl;
    cout << "px2 = " << input2.Px() << ", py2 = " << input2.Py() << ", pz2 = " << input2.Pz() << ", E = " << input2.E() << endl;
    cout << " " << endl;
  }


  // --- Only now apply the boost + rotation to the particle vector (see option part)
  // now change vec either (option==1 from the Lab to Com) or (option==2 from the Com to Lab)
  if (option == 1) {
    // Lab to CoM - order is boost and then rotate (identical to beam manipulations)
    vec.Boost(-cms.BoostVector()); 
    vec.RotateY(rotAngleY);
  } else if (option == 2) {
    // CoM to Lab - order is rotate and boost (opposite sign of beam manipulations)
    vec.RotateY(-rotAngleY);
    vec.Boost(cms.BoostVector()); 
  } else {
    cout << "boost_and_rotate:   Error with option not equal to 1 or 2:  option = " << option << endl;
  }

  // not really used here.... as return value
  return rotAngleY;

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

void initialize_pmt_position();

void post_ana_bbc_fvtx(int runNumber = 435823, int rp_recal_pass = 1){

  // --- need some agreed upon standard for dealing with paths

  char outFile1[300];
  sprintf(outFile1,"%s%d%s","/phenix/plhf/theok/taxi/Run15pAu200FVTXClusAna503/8819/processed/hist_",runNumber,".root"); // absolute paths need to be dealt with

  char outFile2[100];
  sprintf(outFile2,"%s%d%s","vtx_ep_calib/shengli/rp/hrp_",runNumber,".root"); // absolute paths need to be dealt with

  cout<<"runNumber = " <<runNumber<<" "
      <<"rp_recal_pass = "<<rp_recal_pass<<endl;

  char filename[500];
  sprintf(filename,"/gpfs/mnt/gpfs02/phenix/plhf/plhf1/theok/taxi/Run15pAu200FVTXClusAna503/8819/data/%d.root",runNumber); // abslute paths need to be dealt with

  cout << "v2 input file: " << filename << endl;
  cout << "v2 output file: " << outFile1 << endl;

  char calibfile[500];
  sprintf(calibfile,"vtx_ep_calib/new_wide_z_flattening/flattening_%d_%d.dat",runNumber,rp_recal_pass-1);

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

  initialize_pmt_position();

  bool fvtx_clusters = true;
  bool bbc_pmts      = true;
  bool vtx_tracks    = true;

  int n_angle_config = 64; // move below, number is 8 * 8 (nothing to do with number of bbc tubes)
  int n_side_angle = sqrt(n_angle_config);
  // --- these things for bookkeeping...???
  int first_bbc_angle = 2;
  int first_fvtx_angle = 64+2;
  int first_fvtx_0_angle = 2*64+2;
  int first_fvtx_1_angle = 3*64+2;
  int first_fvtx_2_angle = 4*64+2;
  int first_fvtx_3_angle = 5*64+2;

  float max_blue_angle   = 0.0032;
  float min_blue_angle   = 0.000;
  float max_yellow_angle = 0.0072;
  float min_yellow_angle = 0.000;

  float pi = acos(-1.0);

  //------------------------------------------------------------//
  //                                                            //
  //       Initializing Calibration Arrays & Histograms         //
  //                                                            //
  //------------------------------------------------------------//

  // --- lots of comments needed here
  TH2D     *qx[NMUL][NHAR][NDET];
  TH2D     *qy[NMUL][NHAR][NDET];
  TProfile *ave[NMUL][NZPS][NHAR][NDET];
  TProfile *flt[NMUL][NZPS][NHAR][NDET];
  TH2D     *dis[NMUL][NHAR][NDET];

  float    mean[NMUL][NZPS][NHAR][NDET][2];
  float    widt[NMUL][NZPS][NHAR][NDET][2];
  float    four[NMUL][NZPS][NHAR][NDET][2][NORD];


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

            }
        }
    }

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
        }
    }


  //Initializing the calibration parameters
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
                    }
                }
            }
        }
    }

  //------------------------------------------------------------//
  //   Finished Initializing Calibration Arrays & Histograms    //
  //------------------------------------------------------------//

  //------------------------------------------------------------//
  //                                                            //
  //      Reading in flattening calibration parameters          //
  //                                                            //
  //------------------------------------------------------------//

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
                  if (f1==0.0 || f1<0.0) f1=1.0;
                  if (f3==0.0 || f3<0.0) f3=1.0;
                  mean[ic][iz][ih][id][0]=f0;
                  widt[ic][iz][ih][id][0]=f1;
                  mean[ic][iz][ih][id][1]=f2;
                  widt[ic][iz][ih][id][1]=f3;
                  if(id==0 && ih == 1 && DIAG) cout<<f0<<" "<<f1<<" "<<f2<<" "<<f3<<endl;//bbc psi 2 parameters
                  for (int ib=0; ib<2; ib++) 
                    {
                      for (int io=0; io<NORD; io++) 
                        {
                          ifs>>four[ic][iz][ih][id][ib][io];
                        }
                    }
                }
            }
        }
      ifs.close();
    }


  //------------------------------------------------------------//
  //  Finished Reading in flattening calibration parameters     //
  //------------------------------------------------------------//

  //------------------------------------------------------------//
  //                  Initializing boost vectors                //
  //------------------------------------------------------------//

  vector<TLorentzVector > blue_vecs;
  vector<TLorentzVector > yellow_vecs;
  float proton_mass = 0.938;

  for(int iangle1 = 0; iangle1 < n_side_angle; iangle1++)
    {
      float blue_angle = min_blue_angle+iangle1*(max_blue_angle-min_blue_angle)/n_side_angle; // always a good idea to cast ints as floats
      float blue_px = 100*TMath::Sin(blue_angle);
      float blue_py = 0.0;
      float blue_pz = 100*TMath::Cos(blue_angle);
      float blue_energy = TMath::Sqrt(100*100+proton_mass*proton_mass);
      blue_vecs.push_back(TLorentzVector(blue_px,blue_py,blue_pz,blue_energy));
    }
  for(int iangle2 = 0; iangle2 < n_side_angle; iangle2++)
    {
      float yellow_angle = TMath::Pi()+min_yellow_angle+iangle2*(max_yellow_angle-min_yellow_angle)/n_side_angle; // always a good idea to cast ints as floats
      float yellow_px = 100*TMath::Sin(yellow_angle);
      float yellow_py = 0.0;
      float yellow_pz = 100*TMath::Cos(yellow_angle);
      float yellow_energy = TMath::Sqrt(100*100+proton_mass*proton_mass);
      yellow_vecs.push_back(TLorentzVector(yellow_px,yellow_py,yellow_pz,yellow_energy));

    }

  //------------------------------------------------------------//
  //                  Initializing histograms                   //
  //------------------------------------------------------------//

  TProfile *bbcs_v2_incl_default = new TProfile("bbcs_v2_incl_default","bbcs_v2_incl_default",15, 0.0, 3.0,-1.1,1.1);
  TProfile *bbcs_v2_east_default = new TProfile("bbcs_v2_east_default","bbcs_v2_east_default",15, 0.0, 3.0,-1.1,1.1);
  TProfile *bbcs_v2_west_default = new TProfile("bbcs_v2_west_default","bbcs_v2_west_default",15, 0.0, 3.0,-1.1,1.1);

  vector<TProfile * > bbcs_v2_west_angle;
  vector<TProfile * > bbcs_v2_east_angle;
  vector<TProfile * > bbcs_v2_incl_angle;

  TProfile *fvtxs_v2_incl_default = new TProfile("fvtxs_v2_incl_default","fvtxs_v2_incl_default",15, 0.0, 3.0,-1.1,1.1);
  TProfile *fvtxs_v2_east_default = new TProfile("fvtxs_v2_east_default","fvtxs_v2_east_default",15, 0.0, 3.0,-1.1,1.1);
  TProfile *fvtxs_v2_west_default = new TProfile("fvtxs_v2_west_default","fvtxs_v2_west_default",15, 0.0, 3.0,-1.1,1.1);

  vector<TProfile * > fvtxs_v2_west_angle;
  vector<TProfile * > fvtxs_v2_east_angle;
  vector<TProfile * > fvtxs_v2_incl_angle;

  vector<TProfile * > fvtxs0_v2_west_angle;
  vector<TProfile * > fvtxs0_v2_east_angle;
  vector<TProfile * > fvtxs0_v2_incl_angle;

  vector<TProfile * > fvtxs1_v2_west_angle;
  vector<TProfile * > fvtxs1_v2_east_angle;
  vector<TProfile * > fvtxs1_v2_incl_angle;

  vector<TProfile * > fvtxs2_v2_west_angle;
  vector<TProfile * > fvtxs2_v2_east_angle;
  vector<TProfile * > fvtxs2_v2_incl_angle;

  vector<TProfile * > fvtxs3_v2_west_angle;
  vector<TProfile * > fvtxs3_v2_east_angle;
  vector<TProfile * > fvtxs3_v2_incl_angle;

  for(int iangle = 0; iangle < n_angle_config; iangle++)
    {
      bbcs_v2_west_angle.push_back(new TProfile(Form("bbcs_v2_west_angle_%d",iangle),Form("bbcs_v2_west_angle_%d",iangle),15, 0.0, 3.0,-1.1,1.1));
      bbcs_v2_east_angle.push_back(new TProfile(Form("bbcs_v2_east_angle_%d",iangle),Form("bbcs_v2_east_angle_%d",iangle),15, 0.0, 3.0,-1.1,1.1));
      bbcs_v2_incl_angle.push_back(new TProfile(Form("bbcs_v2_incl_angle_%d",iangle),Form("bbcs_v2_incl_angle_%d",iangle),15, 0.0, 3.0,-1.1,1.1));

      fvtxs_v2_west_angle.push_back(new TProfile(Form("fvtxs_v2_west_angle_%d",iangle),Form("fvtxs_v2_west_angle_%d",iangle),15, 0.0, 3.0,-1.1,1.1));
      fvtxs_v2_east_angle.push_back(new TProfile(Form("fvtxs_v2_east_angle_%d",iangle),Form("fvtxs_v2_east_angle_%d",iangle),15, 0.0, 3.0,-1.1,1.1));
      fvtxs_v2_incl_angle.push_back(new TProfile(Form("fvtxs_v2_incl_angle_%d",iangle),Form("fvtxs_v2_incl_angle_%d",iangle),15, 0.0, 3.0,-1.1,1.1));

      fvtxs0_v2_west_angle.push_back(new TProfile(Form("fvtxs0_v2_west_angle_%d",iangle),Form("fvtxs0_v2_west_angle_%d",iangle),15, 0.0, 3.0,-1.1,1.1));
      fvtxs0_v2_east_angle.push_back(new TProfile(Form("fvtxs0_v2_east_angle_%d",iangle),Form("fvtxs0_v2_east_angle_%d",iangle),15, 0.0, 3.0,-1.1,1.1));
      fvtxs0_v2_incl_angle.push_back(new TProfile(Form("fvtxs0_v2_incl_angle_%d",iangle),Form("fvtxs0_v2_incl_angle_%d",iangle),15, 0.0, 3.0,-1.1,1.1));

      fvtxs1_v2_west_angle.push_back(new TProfile(Form("fvtxs1_v2_west_angle_%d",iangle),Form("fvtxs1_v2_west_angle_%d",iangle),15, 0.0, 3.0,-1.1,1.1));
      fvtxs1_v2_east_angle.push_back(new TProfile(Form("fvtxs1_v2_east_angle_%d",iangle),Form("fvtxs1_v2_east_angle_%d",iangle),15, 0.0, 3.0,-1.1,1.1));
      fvtxs1_v2_incl_angle.push_back(new TProfile(Form("fvtxs1_v2_incl_angle_%d",iangle),Form("fvtxs1_v2_incl_angle_%d",iangle),15, 0.0, 3.0,-1.1,1.1));

      fvtxs2_v2_west_angle.push_back(new TProfile(Form("fvtxs2_v2_west_angle_%d",iangle),Form("fvtxs2_v2_west_angle_%d",iangle),15, 0.0, 3.0,-1.1,1.1));
      fvtxs2_v2_east_angle.push_back(new TProfile(Form("fvtxs2_v2_east_angle_%d",iangle),Form("fvtxs2_v2_east_angle_%d",iangle),15, 0.0, 3.0,-1.1,1.1));
      fvtxs2_v2_incl_angle.push_back(new TProfile(Form("fvtxs2_v2_incl_angle_%d",iangle),Form("fvtxs2_v2_incl_angle_%d",iangle),15, 0.0, 3.0,-1.1,1.1));

      fvtxs3_v2_west_angle.push_back(new TProfile(Form("fvtxs3_v2_west_angle_%d",iangle),Form("fvtxs3_v2_west_angle_%d",iangle),15, 0.0, 3.0,-1.1,1.1));
      fvtxs3_v2_east_angle.push_back(new TProfile(Form("fvtxs3_v2_east_angle_%d",iangle),Form("fvtxs3_v2_east_angle_%d",iangle),15, 0.0, 3.0,-1.1,1.1));
      fvtxs3_v2_incl_angle.push_back(new TProfile(Form("fvtxs3_v2_incl_angle_%d",iangle),Form("fvtxs3_v2_incl_angle_%d",iangle),15, 0.0, 3.0,-1.1,1.1));

    }

  //------------------------------------------------------------//
  //               Initializing Tree Variables                  //
  //------------------------------------------------------------//

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

  TBranch *b_nsegments = htree->GetBranch("nsegments");
  TBranch *b_px = htree->GetBranch("px");
  TBranch *b_py = htree->GetBranch("py");
  TBranch *b_pz = htree->GetBranch("pz");

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

  int nentries = htree->GetEntries();
  cout<<"total events = " << nentries<<endl;
  for ( int ievt = 0 ; ievt < nentries ; ievt++ ) {
    
    if(ievt%1000==0) cout<<"event number = "<<ievt<<endl;
    //if(ievt ==1000) break;

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

    if(bbc_pmts)
      b_d_BBC_charge->GetEntry(ievt);//bbc pmts charge

    //FVTX clusters    
    if(fvtx_clusters)
      {
        b_d_nFVTX_clus->GetEntry(ievt);
        b_d_FVTX_x->GetEntry(ievt);
        b_d_FVTX_y->GetEntry(ievt);
        b_d_FVTX_z->GetEntry(ievt);
      }

    //VTX Tracks
    if(vtx_tracks)
      {
        b_nsegments->GetEntry(ievt);
        b_px->GetEntry(ievt);
        b_py->GetEntry(ievt);
        b_pz->GetEntry(ievt);
      }

    //if(vtx_z!=vtx_z) continue;
    //if( fabs(vtx_z) > 100) continue;
    // --- some big questions here about the z-vertex cut
    // --- really need to double and triple check where the zvertex cuts are applied and what they are
    int ibbcz  = NZPS*(d_bbcz+30)/60;//bbc z bin for -30 <bbc z < 30 // how do you specify the number of bins here

    // --- break and continue statements should happen much, much earlier --------------------
    if(rp_recal_pass<1 || rp_recal_pass > 3) break;// rp_recal_pass only valid between 1 and 3

    if(ibbcz<0||ibbcz>=NZPS) continue;

    if(d_nsegments==0 && rp_recal_pass > 2) continue;
    // ---------------------------------------------------------------------------------------

    //------------------------------------------------------------//
    //                Calculating Event Planes                    //
    //------------------------------------------------------------//

    float vtx_x = bc_x + d_bbcz*0.025/10; // z dependent x position because of beam angle rotation issues // what are these numbers?  // these are really specific to p+Au
    float vtx_y = bc_y;

    float bbc_qx[n_angle_config];
    float bbc_qy[n_angle_config];
    float bbc_qw             = 0;

    for(int iangle = 0; iangle < n_angle_config; iangle++)
      {
        bbc_qx[iangle] = 0;
        bbc_qy[iangle] = 0;
      }

    if(bbc_pmts)
      for(int iangle = 0; iangle < n_angle_config; iangle++)
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

              // --- need explanation of these two things
              int iangle_blue = iangle/n_side_angle;
              int iangle_yellow = iangle%n_side_angle;

              boost_and_rotate(particle_vec,blue_vecs[iangle_blue],yellow_vecs[iangle_yellow],1);// option 1 Lab to CoM

              phi = TMath::ATan2(particle_vec.Py(),particle_vec.Px());

              bbc_qx[iangle] += bbc_charge*TMath::Cos(2*phi);
              bbc_qy[iangle] += bbc_charge*TMath::Sin(2*phi);
              if(iangle == 0)
                bbc_qw += bbc_charge;
            }
        }

    float fvtx_qx[5][n_angle_config];//all layers then 0 1 2 3
    float fvtx_qy[5][n_angle_config];
    float fvtx_qw[5]                ;

    for(int iangle = 0; iangle < n_angle_config; iangle++)
      {
        for(int ilayer = 0; ilayer < 5; ilayer++)
          {
            fvtx_qx[ilayer][iangle] = 0.0;
            fvtx_qy[ilayer][iangle] = 0.0;
            if(iangle==0)
              fvtx_qw[ilayer] = 0.0;
          }
      }

    if(fvtx_clusters)
      for(int iangle = 0; iangle< n_angle_config; iangle++)
        {
          for(int iclus = 0; iclus < d_nFVTX_clus; iclus++)
            {
              float fvtx_x      = d_FVTX_x[iclus] - vtx_x; // calculate for each event, function of z
              float fvtx_y      = d_FVTX_y[iclus] - vtx_y;
              float fvtx_z      = d_FVTX_z[iclus];

              double fvtx_r = sqrt(pow(fvtx_x,2.0)+pow(fvtx_y,2.0));

              double fvtx_the = atan2(fvtx_r,fvtx_z - d_bbcz); //fvtx_z-bbcv // add a new variable to make it clear that you're using a corrected z vertex
              double fvtx_eta = -log(tan(0.5*fvtx_the));

              // --- z dependent eta cut for the track, but not for the clusters???
              if(!(fabs(fvtx_eta)>1.0 && fabs(fvtx_eta)<3.5)) continue;
              // cout<<"fvtx_x: "<<fvtx_x<<" fvtx_y: "<<fvtx_y<<" fvtx_z"<<fvtx_z<<endl;

              int fvtx_layer    = get_fvtx_layer(fvtx_z);

              int igap = (fabs(fvtx_eta)-1.0)/0.5;

              int id_fvtx = fvtx_layer*5+igap;

              if(!(id_fvtx>=0 && id_fvtx<40)) continue;

              float phi = TMath::ATan2(fvtx_y,fvtx_x);

              float mass = 0.1396;//assume charged pion mass
              float pT = 0.25;
              float px = pT * TMath::Cos(phi);
              float py = pT * TMath::Sin(phi);
              float pz = pT * TMath::SinH(fvtx_eta);

              float energy = TMath::Sqrt(px*px+py*py+pz*pz+mass*mass);
              TLorentzVector particle_vec(px,py,pz,energy);

              int iangle_blue = iangle/n_side_angle;
              int iangle_yellow = iangle%n_side_angle;

              boost_and_rotate(particle_vec,blue_vecs[iangle_blue],yellow_vecs[iangle_yellow],1);// option 1 Lab to CoM

              phi = TMath::ATan2(particle_vec.Py(),particle_vec.Px());

              fvtx_qx[fvtx_layer+1][iangle] += TMath::Cos(2*phi);
              fvtx_qy[fvtx_layer+1][iangle] += TMath::Sin(2*phi);

              fvtx_qx[0][iangle] += TMath::Cos(2*phi);
              fvtx_qy[0][iangle] += TMath::Sin(2*phi);

              if(iangle==0)
                { 
                  fvtx_qw[fvtx_layer+1] ++;
                  fvtx_qw[0] ++;
                }

            }
        }

    float sumxy[NHAR][NDET][4];
    for (int i=0; i<NHAR; i++) 
      {
        for (int j=0; j<NDET; j++) 
          {
            for (int k=0; k<4; k++)//qx qy qw psi
              {
                sumxy[i][j][k]=0;
              }
          }
      }

    //save bbc q vec 2
    //cout<<"bbc from node tree: "<<d_Qx[5]<<" "<<d_Qy[5]<<" "<<d_Qw[5]<<endl;
    sumxy[1][0][0] = d_Qx[5];
    sumxy[1][0][1] = d_Qy[5];
    sumxy[1][0][2] = d_Qw[5];

    //save fvtx q vec 2
    sumxy[1][1][0] = d_Qx[4];
    sumxy[1][1][1] = d_Qy[4];
    sumxy[1][1][2] = d_Qw[4];

    if(fvtx_clusters || bbc_pmts)
      for(int iangle = 0; iangle < n_angle_config; iangle++)
        {
          if(bbc_pmts)
            {
              sumxy[1][iangle+first_bbc_angle][0] = bbc_qx[iangle];
              sumxy[1][iangle+first_bbc_angle][1] = bbc_qy[iangle];
              sumxy[1][iangle+first_bbc_angle][2] = bbc_qw;
            }

          if(fvtx_clusters)
            {
              sumxy[1][iangle+first_fvtx_angle][0] = fvtx_qx[0][iangle];
              sumxy[1][iangle+first_fvtx_angle][1] = fvtx_qy[0][iangle];
              sumxy[1][iangle+first_fvtx_angle][2] = fvtx_qw[0];

              sumxy[1][iangle+first_fvtx_0_angle][0] = fvtx_qx[1][iangle];
              sumxy[1][iangle+first_fvtx_0_angle][1] = fvtx_qy[1][iangle];
              sumxy[1][iangle+first_fvtx_0_angle][2] = fvtx_qw[1];

              sumxy[1][iangle+first_fvtx_1_angle][0] = fvtx_qx[2][iangle];
              sumxy[1][iangle+first_fvtx_1_angle][1] = fvtx_qy[2][iangle];
              sumxy[1][iangle+first_fvtx_1_angle][2] = fvtx_qw[2];

              sumxy[1][iangle+first_fvtx_2_angle][0] = fvtx_qx[3][iangle];
              sumxy[1][iangle+first_fvtx_2_angle][1] = fvtx_qy[3][iangle];
              sumxy[1][iangle+first_fvtx_2_angle][2] = fvtx_qw[3];

              sumxy[1][iangle+first_fvtx_3_angle][0] = fvtx_qx[4][iangle];
              sumxy[1][iangle+first_fvtx_3_angle][1] = fvtx_qy[4][iangle];
              sumxy[1][iangle+first_fvtx_3_angle][2] = fvtx_qw[4];
            }
        }

    if(DIAG)
      {
        cout<<"bbc from node tree: "<<d_Qx[5]<<" "<<d_Qy[5]<<" "<<d_Qw[5]<<endl;
        cout<<"bbc from me: "<<bbc_qx[0]<<" "<<bbc_qy[0]<<" "<<bbc_qw<<endl;

        cout<<"fvtx raw: "<<endl;
        cout<<"from node tree: "<<d_Qx[4]<<" "<<d_Qy[4]<<" "<<d_Qw[4]<<endl;
        cout<<"from clusters: " <<fvtx_qx[0][0]<<" "<<fvtx_qy[0][0]<<" "<<fvtx_qw[0]<<endl;
      }
    //------------------------------------------------------------//
    //                Flattening iteration                        //
    //------------------------------------------------------------//
    //int icent = 0;
    for (int ih=1; ih<NHAR; ih++) {
      for(int id=0;id<NDET;id++){
        if (sumxy[ih][id][2]>0.0) {
          sumxy[ih][id][3]=atan2(sumxy[ih][id][1],sumxy[ih][id][0])/(ih+1.0);
          if (rp_recal_pass>0) dis[icent][ih][id]->Fill(ibbcz,sumxy[ih][id][3]*(ih+1.0));
        }
        if (sumxy[ih][id][2]>0.0) {
          for (int ib=0; ib<2; ib++) {
            sumxy[ih][id][ib]/=sumxy[ih][id][2];
            //if(ih==1 && id==0 && ib==0 && sumxy[ih][id][ib]>1) cout<<sumxy[ih][id][ib]<<endl;
            if (rp_recal_pass>0) {
              ave[icent][ibbcz][ih][id]->Fill(ib+0.0,sumxy[ih][id][ib]);
              if(id==0 && DIAG) cout<<"filled ave: "<<ih<<" "<<id<<" "<<ib<<" with: "<<sumxy[ih][id][ib]<<endl;
              if(ib==0) qx[icent][ih][id]->Fill(ibbcz,sumxy[ih][id][0]);
              if(ib==1) qy[icent][ih][id]->Fill(ibbcz,sumxy[ih][id][1]);
            }
            float sxy=sumxy[ih][id][ib];
            float mxy=mean[icent][ibbcz][ih][id][ib];
            float wxy=widt[icent][ibbcz][ih][id][ib];

            //if(icent==0 && ibbcz==0 && ih==1 && id==0) cout<<ib<<" "<<sxy<<" "<<mxy<<" "<<wxy<<endl;
            sumxy[ih][id][ib]=(sxy-mxy)/wxy;
            if (rp_recal_pass>0) {
              ave[icent][ibbcz][ih][id]->Fill(ib+2.0,sumxy[ih][id][ib]);
              if(id==0 && DIAG) cout<<"filled ave2: "<<ih<<" "<<id<<" "<<ib<<" with: "<<sumxy[ih][id][ib]<<endl;
              if(ib==0) qx[icent][ih][id]->Fill(ibbcz+NZPS,sumxy[ih][id][0]);
              if(ib==1) qy[icent][ih][id]->Fill(ibbcz+NZPS,sumxy[ih][id][1]);
            }
          }

          sumxy[ih][id][3]=atan2(sumxy[ih][id][1],sumxy[ih][id][0])/(ih+1.0);
          if (rp_recal_pass>0) {
            dis[icent][ih][id]->Fill(ibbcz+NZPS,sumxy[ih][id][3]*(ih+1.0));
          }

          float psi=sumxy[ih][id][3]*(ih+1.0);
          if(ih==1 && id==0 && DIAG)  cout<<"psi2-1 bbc: "<<psi<<endl;
          float dp=0.0;
          for (int io=0; io<NORD; io++) {
            float cc=cos((io+1.0)*psi);
            float ss=sin((io+1.0)*psi);
            if (rp_recal_pass>0) flt[icent][ibbcz][ih][id]->Fill(io+0.0,cc);
            if (rp_recal_pass>0) flt[icent][ibbcz][ih][id]->Fill(io+NORD,ss);
            float aa=four[icent][ibbcz][ih][id][0][io]; // mean cos
            float bb=four[icent][ibbcz][ih][id][1][io]; // mean sin
            dp+=(aa*ss-bb*cc)*2.0/(io+1.0);
          }
          psi+=dp;
          psi=atan2(sin(psi),cos(psi));
          if(ih==1 && id==0 && DIAG)  cout<<"psi2-2 bbc: "<<psi<<endl;
          for (int io=0; io<NORD; io++) {
            float cc=cos((io+1.0)*psi);
            float ss=sin((io+1.0)*psi);
            if (rp_recal_pass>0) flt[icent][ibbcz][ih][id]->Fill(io+NORD*2.0,cc);
            if (rp_recal_pass>0) flt[icent][ibbcz][ih][id]->Fill(io+NORD*3.0,ss);
          }
          sumxy[ih][id][3]=psi/(ih+1.0);
          if (rp_recal_pass>0) dis[icent][ih][id]->Fill(ibbcz+NZPS*2.0,sumxy[ih][id][3]*(ih+1.0));
        } else {
          sumxy[ih][id][3]=-9999.9;
        }
      }//end of id
    }//end of ih
      
    if(rp_recal_pass<3) continue;

    float bbc_psi2 = (sumxy[1][0][2]>0)?sumxy[1][0][3]:-9999.9;
    float fvtx_psi2 = (sumxy[1][1][2]>12)?sumxy[1][1][3]:-9999.9;

    float bbc_psi2_angle[n_angle_config];

    float fvtx_psi2_angle[n_angle_config];
    float fvtx0_psi2_angle[n_angle_config];
    float fvtx1_psi2_angle[n_angle_config];
    float fvtx2_psi2_angle[n_angle_config];
    float fvtx3_psi2_angle[n_angle_config];

    if(fvtx_clusters || bbc_pmts)
      for(int iangle = 0; iangle < n_angle_config; iangle++)
        {
          if(bbc_pmts)
            bbc_psi2_angle[iangle] = (sumxy[1][iangle+first_bbc_angle][2]>0)?sumxy[1][iangle+first_bbc_angle][3]:-9999.9;

          if(fvtx_clusters)
            {
              fvtx_psi2_angle[iangle] = (sumxy[1][iangle+first_fvtx_angle][2]>4)?sumxy[1][iangle+first_fvtx_angle][3]:-9999.9;
              fvtx0_psi2_angle[iangle] = (sumxy[1][iangle+first_fvtx_0_angle][2]>4)?sumxy[1][iangle+first_fvtx_0_angle][3]:-9999.9;
              fvtx1_psi2_angle[iangle] = (sumxy[1][iangle+first_fvtx_1_angle][2]>4)?sumxy[1][iangle+first_fvtx_1_angle][3]:-9999.9;
              fvtx2_psi2_angle[iangle] = (sumxy[1][iangle+first_fvtx_2_angle][2]>4)?sumxy[1][iangle+first_fvtx_2_angle][3]:-9999.9;
              fvtx3_psi2_angle[iangle] = (sumxy[1][iangle+first_fvtx_3_angle][2]>4)?sumxy[1][iangle+first_fvtx_3_angle][3]:-9999.9;
            }
        }


    //start of vtx stand alone track loop
    if(vtx_tracks)
      for(int itrk=0; itrk< d_nsegments; itrk++)
        {
          float px    = d_px[itrk];
          float py    = d_py[itrk];
          float pz    = d_pz[itrk];
          int dcarm=0;
          if(px>0) dcarm=1;

          float phi0 = TMath::ATan2(py,px);
          float pt = sqrt(px*px+py*py);

          float bbc_dphi_2 = phi0 - bbc_psi2; // move this lower?

          if(-4.0<bbc_psi2 && bbc_psi2<4.0 ) // why this weird cut? why not just -pi to pi? // checking against -9999 from above
            {
              bbcs_v2_incl_default->Fill(pt,cos(2*bbc_dphi_2));
              if(dcarm==0) 
                {
                  bbcs_v2_east_default->Fill(pt,cos(2*bbc_dphi_2));
                }
              else if(dcarm==1)
                {
                  bbcs_v2_west_default->Fill(pt,cos(2*bbc_dphi_2));
                }
            }

          float fvtx_dphi_2 = phi0 - fvtx_psi2;

          if(-4.0<fvtx_psi2 && fvtx_psi2<4.0 )
            {
              fvtxs_v2_incl_default->Fill(pt,cos(2*fvtx_dphi_2));
              if(dcarm==0) 
                {
                  fvtxs_v2_east_default->Fill(pt,cos(2*fvtx_dphi_2));
                }
              else if(dcarm==1)
                {
                  fvtxs_v2_west_default->Fill(pt,cos(2*fvtx_dphi_2));
                }
            }

          if(fvtx_clusters || bbc_pmts)
            for(int iangle = 0; iangle < n_angle_config; iangle++)
              {

                float mass = 0.1396;//assume charged pion mass
                float energy = TMath::Sqrt(px*px+py*py+pz*pz+mass*mass);
                TLorentzVector particle_vec(px,py,pz,energy);

                int iangle_blue = iangle/n_side_angle;
                int iangle_yellow = iangle%n_side_angle;

                boost_and_rotate(particle_vec,blue_vecs[iangle_blue],yellow_vecs[iangle_yellow],1);// option 1 Lab to CoM

                float phi_angle = TMath::ATan2(particle_vec.Py(),particle_vec.Px()); // rotated phi // "modified"
                float pt_angle = TMath::Sqrt(particle_vec.Py()*particle_vec.Py()+particle_vec.Px()*particle_vec.Px()); // rotated pt "modified"

                //bbc angle
                if(bbc_pmts)
                  if(-4.0<bbc_psi2_angle[iangle] && bbc_psi2_angle[iangle]<4.0)
                    {
                      double bbc_dphi2_angle = phi_angle - bbc_psi2_angle[iangle];
                      double cosbbc_dphi2_angle = TMath::Cos(2*bbc_dphi2_angle);

                      if(dcarm==1)
                        {
                          bbcs_v2_west_angle[iangle]->Fill(pt_angle,cosbbc_dphi2_angle);
                        }
                      else if( dcarm==0)
                        {
                          bbcs_v2_east_angle[iangle]->Fill(pt_angle,cosbbc_dphi2_angle);
                        }
                      bbcs_v2_incl_angle[iangle]->Fill(pt_angle,cosbbc_dphi2_angle);
                    }

                if(!fvtx_clusters) continue;
                //fvtx all layers
                if(-4.0<fvtx_psi2_angle[iangle] && fvtx_psi2_angle[iangle]<4.0)
                  {
                    double fvtx_dphi2_angle = phi_angle - fvtx_psi2_angle[iangle];
                    double cosfvtx_dphi2_angle = TMath::Cos(2*fvtx_dphi2_angle);

                    if(dcarm==1)
                      {
                        fvtxs_v2_west_angle[iangle]->Fill(pt_angle,cosfvtx_dphi2_angle);
                      }
                    else if( dcarm==0)
                      {
                        fvtxs_v2_east_angle[iangle]->Fill(pt_angle,cosfvtx_dphi2_angle);
                      }
                    fvtxs_v2_incl_angle[iangle]->Fill(pt_angle,cosfvtx_dphi2_angle);
                  }

                //fvtx layer 0
                if(-4.0<fvtx0_psi2_angle[iangle] && fvtx0_psi2_angle[iangle]<4.0)
                  {
                    double fvtx_dphi2_angle = phi_angle - fvtx0_psi2_angle[iangle];
                    double cosfvtx_dphi2_angle = TMath::Cos(2*fvtx_dphi2_angle);

                    if(dcarm==1)
                      {
                        fvtxs0_v2_west_angle[iangle]->Fill(pt_angle,cosfvtx_dphi2_angle);
                      }
                    else if( dcarm==0)
                      {
                        fvtxs0_v2_east_angle[iangle]->Fill(pt_angle,cosfvtx_dphi2_angle);
                      }
                    fvtxs0_v2_incl_angle[iangle]->Fill(pt_angle,cosfvtx_dphi2_angle);
                  }

                //fvtx layer 1
                if(-4.0<fvtx1_psi2_angle[iangle] && fvtx1_psi2_angle[iangle]<4.0)
                  {
                    double fvtx_dphi2_angle = phi_angle - fvtx1_psi2_angle[iangle];
                    double cosfvtx_dphi2_angle = TMath::Cos(2*fvtx_dphi2_angle);

                    if(dcarm==1)
                      {
                        fvtxs1_v2_west_angle[iangle]->Fill(pt_angle,cosfvtx_dphi2_angle);
                      }
                    else if( dcarm==0)
                      {
                        fvtxs1_v2_east_angle[iangle]->Fill(pt_angle,cosfvtx_dphi2_angle);
                      }
                    fvtxs1_v2_incl_angle[iangle]->Fill(pt_angle,cosfvtx_dphi2_angle);
                  }

                //fvtx layer 2
                if(-4.0<fvtx2_psi2_angle[iangle] && fvtx2_psi2_angle[iangle]<4.0)
                  {
                    double fvtx_dphi2_angle = phi_angle - fvtx2_psi2_angle[iangle];
                    double cosfvtx_dphi2_angle = TMath::Cos(2*fvtx_dphi2_angle);

                    if(dcarm==1)
                      {
                        fvtxs2_v2_west_angle[iangle]->Fill(pt_angle,cosfvtx_dphi2_angle);
                      }
                    else if( dcarm==0)
                      {
                        fvtxs2_v2_east_angle[iangle]->Fill(pt_angle,cosfvtx_dphi2_angle);
                      }
                    fvtxs2_v2_incl_angle[iangle]->Fill(pt_angle,cosfvtx_dphi2_angle);
                  }

                //fvtx layer 3
                if(-4.0<fvtx3_psi2_angle[iangle] && fvtx3_psi2_angle[iangle]<4.0)
                  {
                    double fvtx_dphi2_angle = phi_angle - fvtx3_psi2_angle[iangle];
                    double cosfvtx_dphi2_angle = TMath::Cos(2*fvtx_dphi2_angle);

                    if(dcarm==1)
                      {
                        fvtxs3_v2_west_angle[iangle]->Fill(pt_angle,cosfvtx_dphi2_angle);
                      }
                    else if( dcarm==0)
                      {
                        fvtxs3_v2_east_angle[iangle]->Fill(pt_angle,cosfvtx_dphi2_angle);
                      }
                    fvtxs3_v2_incl_angle[iangle]->Fill(pt_angle,cosfvtx_dphi2_angle);
                  }
              }
        }
  }//end of event 
    
  //}

  if(rp_recal_pass<3 && rp_recal_pass>0)
    {
      //sprintf(calibfile,"vtx_ep_calib/new_wide_z_flattening/flattening_%d_%d.dat",runNumber,rp_recal_pass);
      sprintf(calibfile,"vtx_ep_calib/new_wide_z_flattening/flattening_%d_%d.dat",runNumber,rp_recal_pass);
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
                      if(DIAG)
                        cout<<"writing ave:  "<<ic<<" "<<iz<<" "<<ih<<" "<<id<<endl;

                      ofs << ave[icent][iz][ih][id]->GetBinContent(ib+1) << " ";
                      ofs << ave[icent][iz][ih][id]->GetBinError  (ib+1) << " ";
                    }
                  ofs << endl;
                  for (int ib=0; ib<2; ib++) 
                    {
                      for (int io=0; io<NORD; io++) 
                        {
                          ofs << flt[icent][iz][ih][id]->GetBinContent(ib*NORD+io+1) << " ";
                        }
                      ofs << endl;
                    }
                }
            }
        }
      ofs.close();
    }


  if(rp_recal_pass>0){
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
              }
          }
      }

    for (int ih=1; ih<NHAR; ih++) 
      {
        for (int id=0; id<NDET; id++) 
          {
            //if(ih==1 && id==4)
            //{
            qx[ic][ih][id]->Write();
            qy[ic][ih][id]->Write();

            dis[ic][ih][id]->Write();
            // }
          }
      }

    mData2->Close();
  }

  if(rp_recal_pass>2)
    {
      TFile *mData1=TFile::Open(outFile1,"recreate");
      mData1->cd();

      bbcs_v2_west_default->Write();
      bbcs_v2_east_default->Write();
      bbcs_v2_incl_default->Write();

      fvtxs_v2_west_default->Write();
      fvtxs_v2_east_default->Write();
      fvtxs_v2_incl_default->Write();

      if(fvtx_clusters || bbc_pmts)
        for(int iangle = 0; iangle < n_angle_config; iangle++)
          {
            if(bbc_pmts)
              bbcs_v2_west_angle[iangle]->Write();
            if(fvtx_clusters)
              {
                fvtxs_v2_west_angle[iangle]->Write();
                fvtxs0_v2_west_angle[iangle]->Write();
                fvtxs1_v2_west_angle[iangle]->Write();
                fvtxs2_v2_west_angle[iangle]->Write();
                fvtxs3_v2_west_angle[iangle]->Write();
              }
          }

      mData1->Close();
    }


  //return;

  cout<<"cleaning up"<<endl;

  htree->Delete(); 
  f->Close(); 
  delete f; 
  /*
    delete bbcs_v2_incl_default;
    delete bbcs_v2_east_default;
    delete bbcs_v2_west_default;

    delete fvtxs_v2_incl_default;
    delete fvtxs_v2_east_default;
    delete fvtxs_v2_west_default;

    for(int iangle = 0; iangle < n_angle_config; iangle++)
    {
    delete bbcs_v2_west_angle[iangle];

    delete fvtxs_v2_west_angle[iangle];
    delete fvtxs0_v2_west_angle[iangle];
    delete fvtxs1_v2_west_angle[iangle];
    delete fvtxs2_v2_west_angle[iangle];
    delete fvtxs3_v2_west_angle[iangle];
    }

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
