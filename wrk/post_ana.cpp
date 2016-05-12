#include <TTree.h>
#include <TProfile.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TBranch.h>
#include <TMath.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLorentzVector.h>
#include <TVector3.h>

#include <fstream>
#include <iostream>

#include "RpPar.h"

using namespace std;

void boost_and_rotate(TLorentzVector & vec, TLorentzVector input1, TLorentzVector input2)
{

  TLorentzVector cms = input1 + input2;

  TVector3 z(0,0,1);

  input1.Boost(-cms.BoostVector());
  input2.Boost(-cms.BoostVector());
  vec.Boost(-cms.BoostVector());

  float rotAngleY = -input1.Angle(z);

  //Now rotate the beams about x to align them with z axis
  vec.RotateY(rotAngleY);

}

void post_ana(int runnumber, int calFlag = 3)
{

  cout<<"run post_ana.cpp runnumber: "<<runnumber<<" and calFlag: "<<calFlag<<endl;

  bool trimmed_tree = true;

  //int calFlag = 1;
  float    mean[NMUL][NZPS][NHAR][NDET][2];
  float    widt[NMUL][NZPS][NHAR][NDET][2];
  float    four[NMUL][NZPS][NHAR][NDET][2][NORD];

  static const int n_bbc_angles = 63;

  //bool include_fvtx_tracks = true;
  bool include_fvtx_tracks = false;
  //TFile *file = TFile::Open(Form("/phenix/u/theok/hhj2/taxi/Run14HeAu200MinBias/6321/data/%d.root",runnumber));
  //TFile *file = TFile::Open(Form("~/plhf/taxi/Run14HeAu200MinBias/6506/data/%d.root",runnumber));
  TFile *file = TFile::Open(Form("/gpfs/mnt/gpfs02/phenix/plhf/plhf1/theok/taxi/Run15pAu200CAMBPro104/8650/data/%d.root",runnumber));
  //TFile *file = TFile::Open(Form("data/%d.root",runnumber));
  //TFile *file = TFile::Open("/phenix/u/theok/test/VTX_event_plane/out/fvtx_test.root");
  if(!file) return;
  TTree *ntp_event = (TTree*)file->Get("ntp_event");
  if(!ntp_event) return;


  //shift calibration parameters from shengli
  // float Qy2mpcs_shift=0.005;
  // float Qy2fvtxs_shift=-0.007;

  // float Qx3mpcs_shift=-0.010;
  // float Qx3fvtxs_shift=0.006;


  float event;
  float centrality;
  float trigger;
  float bbc_qn;
  float bbc_qs;
  float bbc_z;
  float vtx_x;
  float vtx_y;
  float vtx_z;
  float fvtx_z;
  float d_Qx[9];
  float d_Qy[9];
  float d_Qw[9];
  float d_BBCs_Q[189];
  int nsegments;
  int eventok;

  TBranch *b_bbc_z = ntp_event->GetBranch("bbc_z");
  b_bbc_z->SetAddress(&bbc_z);

  TBranch *b_d_Qx = ntp_event->GetBranch("d_Qx");
  b_d_Qx->SetAddress(&d_Qx);


  TBranch *b_d_Qy = ntp_event->GetBranch("d_Qy");
  b_d_Qy->SetAddress(&d_Qy);


  TBranch *b_d_Qw = ntp_event->GetBranch("d_Qw");
  b_d_Qw->SetAddress(&d_Qw);


  TBranch *b_d_BBCs_Q = ntp_event->GetBranch("d_BBCs_Q");
  b_d_BBCs_Q->SetAddress(&d_BBCs_Q);

  TBranch *b_event = ntp_event->GetBranch("event");
  TBranch *b_centrality = ntp_event->GetBranch("centrality");
  TBranch *b_trigger = ntp_event->GetBranch("trigger");
  TBranch *b_bbc_qn = ntp_event->GetBranch("bbc_qn");
  TBranch *b_bbc_qs = ntp_event->GetBranch("bbc_qs");
  TBranch *b_vtx_x = ntp_event->GetBranch("vtx_x");
  TBranch *b_vtx_y = ntp_event->GetBranch("vtx_y");
  TBranch *b_vtx_z = ntp_event->GetBranch("vtx_z");
  TBranch *b_fvtx_z = ntp_event->GetBranch("fvtx_z");

  if(!trimmed_tree)
  {
    b_event->SetAddress(&event);

    b_centrality->SetAddress(&centrality);

    b_trigger->SetAddress(&trigger);

    b_bbc_qn->SetAddress(&bbc_qn);

    b_bbc_qs->SetAddress(&bbc_qs);

    b_vtx_x->SetAddress(&vtx_x);


    b_vtx_y->SetAddress(&vtx_y);


    b_vtx_z->SetAddress(&vtx_z);


    b_fvtx_z->SetAddress(&fvtx_z);

  }


  TBranch *b_nsegments = ntp_event->GetBranch("nsegments");
  b_nsegments->SetAddress(&nsegments);


  int trackID[10];
  int charge[10];
  float chisq[10];
  int ndf[10];
  int nhit0[10];
  int nhit1[10];
  int nhit2[10];
  int nhit3[10];
  float px[10];
  float py[10];
  float pz[10];
  float dca[10];
  float dca2d[10];
  int segmentok[10];

  TBranch *b_px = ntp_event->GetBranch("px");
  b_px->SetAddress(&px);


  TBranch *b_py = ntp_event->GetBranch("py");
  b_py->SetAddress(&py);


  TBranch *b_pz = ntp_event->GetBranch("pz");
  b_pz->SetAddress(&pz);

  TBranch *b_eventok = ntp_event->GetBranch("eventok");
  TBranch *b_trackID = ntp_event->GetBranch("trackID");
  TBranch *b_charge = ntp_event->GetBranch("charge");
  TBranch *b_chisq = ntp_event->GetBranch("chisq");
  TBranch *b_ndf = ntp_event->GetBranch("ndf");
  TBranch *b_nhit0 = ntp_event->GetBranch("nhit0");
  TBranch *b_nhit1 = ntp_event->GetBranch("nhit1");
  TBranch *b_nhit2 = ntp_event->GetBranch("nhit2");
  TBranch *b_nhit3 = ntp_event->GetBranch("nhit3");
  TBranch *b_dca = ntp_event->GetBranch("dca");
  TBranch *b_dca2d = ntp_event->GetBranch("dca2d");
  TBranch *b_segmentok = ntp_event->GetBranch("segmentok");

  if(!trimmed_tree)
  {
    b_eventok->SetAddress(&eventok);

  //if(!include_fvtx_tracks)
  //{
    b_trackID->SetAddress(&trackID);

    b_charge->SetAddress(&charge);

    b_chisq->SetAddress(&chisq);

    b_ndf->SetAddress(&ndf);

    b_nhit0->SetAddress(&nhit0);

    b_nhit1->SetAddress(&nhit1);

    b_nhit2->SetAddress(&nhit2);

    b_nhit3->SetAddress(&nhit3);

    b_dca->SetAddress(&dca);

    b_dca2d->SetAddress(&dca2d);

    b_segmentok->SetAddress(&segmentok);
  }
  //}
  /*
  int ntracklets;
  float feta[75];
  float fphi[75];
  float fchisq[75];
  int farm[75];
  int fnhits[75];
  float fDCA_X[75];
  float fDCA_Y[75];
*/

  //if(include_fvtx_tracks)
  //{
/*
  TBranch *b_ntracklets = ntp_event->GetBranch("ntracklets");
  b_ntracklets->SetAddress(&ntracklets);


  TBranch *b_feta = ntp_event->GetBranch("feta");
  b_feta->SetAddress(&feta);


  TBranch *b_fphi = ntp_event->GetBranch("fphi");
  b_fphi->SetAddress(&fphi);


  TBranch *b_fchisq = ntp_event->GetBranch("fchisq");
  b_fchisq->SetAddress(&fchisq);

  TBranch *b_fnhits = ntp_event->GetBranch("fnhits");
  b_fnhits->SetAddress(&fnhits);
  TBranch *b_farm = ntp_event->GetBranch("farm");
  b_farm->SetAddress(&farm);

  TBranch *b_fDCA_X = ntp_event->GetBranch("fDCA_X");
  b_fDCA_X->SetAddress(&fDCA_X);

  TBranch *b_fDCA_Y = ntp_event->GetBranch("fDCA_Y");
  b_fDCA_Y->SetAddress(&fDCA_Y);
*/

  //}

  for (int ic=0; ic<NMUL; ic++) {
    for (int iz=0; iz<NZPS; iz++) {
      for (int ih=0; ih<NHAR; ih++) {
       for (int id=0; id<NDET; id++) {
         for (int ib=0; ib<2; ib++) {
           mean[ic][iz][ih][id][ib]=0.0;
           widt[ic][iz][ih][id][ib]=1.0;

           for (int io=0; io<NORD; io++) {
             four[ic][iz][ih][id][ib][io]=0.0;
           }
         }
       }
     }
   }
 }

 float    bbc_angle_mean[NZPS][n_bbc_angles][2];
 float    bbc_angle_widt[NZPS][n_bbc_angles][2];
 float    bbc_angle_four[NZPS][n_bbc_angles][2][NORD];

 for(int iz = 0; iz < NZPS; iz++)
 {
  for(int iangle = 0; iangle < n_bbc_angles; iangle++)
  {
    for(int ib = 0; ib < 2; ib++)
    {
      bbc_angle_mean[iz][iangle][ib]=0.0;
      bbc_angle_widt[iz][iangle][ib]=1.0;
      for(int io = 0; io < NORD; io++)
      {
        bbc_angle_four[iz][iangle][ib][io]=0.0;
      }
    }
  }
}

if(calFlag>1)
{
      float f0,f1,f2,f3;//f4,f5,f6,f7;
      ifstream ifs;
      //ifs.open(Form("tsukuba/flattening_%d_%d.dat",runnumber,calFlag-1));
      ifs.open(Form("vtx_ep_calib/new_wide_z_flattening/flattening_%d_%d.dat",runnumber,calFlag-1));

      if(!ifs)
      {
    //cout<<"ERROR failed to open file: "<<Form("vtx_ep_calib/flattening_%d_2.dat",runnumber)<<endl;
       cout<<"ERROR failed to open file: "<<Form("vtx_ep_calib/new_wide_z_flattening/flattening_%d_2.dat",runnumber)<<endl;
       cout<<"EXITING"<<endl;
       return;
     }
     for (int ic=0; ic<NMUL; ic++) {
       for (int iz=0; iz<NZPS; iz++) {
         for (int ih=0; ih<NHAR; ih++) {
           for (int id=0; id<NDET; id++) {
        //for (int ib=0; ib<2; ib++) {
        //mean[ic][iz][ih][id][ib]=0.0;
        //widt[ic][iz][ih][id][ib]=1.0;
        //}
             ifs >> f0 >> f1 >> f2 >> f3;
             if (f1==0.0 || f1<0.0) f1=1.0;
             if (f3==0.0 || f3<0.0) f3=1.0;
             mean[ic][iz][ih][id][0]=f0;
             widt[ic][iz][ih][id][0]=f1;
             mean[ic][iz][ih][id][1]=f2;
             widt[ic][iz][ih][id][1]=f3;

             if(ih==1 && id==4 && iz == 0)
             {
              cout<<"for iz: "<<iz<<" f0: "<<f0<<" f1: "<<f1<<" f2: "<<f2<<" f3: "<<f3<<endl;
            }

            for (int ib=0; ib<2; ib++) {
    //for (int io=0; io<8; io++) {
    //four[ic][iz][ih][id][ib][io]=0.0;
    //}
              for (int io=0; io<NORD; io++) {
                ifs>>four[ic][iz][ih][id][ib][io];
              }
    //ifs >> f0 >> f1 >> f2 >> f3 >> f4 >> f5 >> f6 >> f7;
    //four[ic][iz][ih][id][ib][0]=f0;
    //four[ic][iz][ih][id][ib][1]=f1;
    //four[ic][iz][ih][id][ib][2]=f2;
    //four[ic][iz][ih][id][ib][3]=f3;
    //four[ic][iz][ih][id][ib][4]=f4;
    //four[ic][iz][ih][id][ib][5]=f5;
    //four[ic][iz][ih][id][ib][6]=f6;
    //four[ic][iz][ih][id][ib][7]=f7;
            }
          }
        }
      }
    }

    ifs.close();

    ifs.open(Form("vtx_ep_calib/new_wide_z_flattening/flattening_bbc_angle_%d_%d.dat",runnumber,calFlag-1));
    for(int iz = 0; iz < NZPS; iz++)
    {
      for(int iangle = 0; iangle < n_bbc_angles; iangle++)
      {
       ifs >> f0 >> f1 >> f2 >> f3;
       if (f1==0.0 || f1<0.0) f1=1.0;
       if (f3==0.0 || f3<0.0) f3=1.0;
       bbc_angle_mean[iz][iangle][0]=f0;
       bbc_angle_widt[iz][iangle][0]=f1;
       bbc_angle_mean[iz][iangle][1]=f2;
       bbc_angle_widt[iz][iangle][1]=f3;

       if(iangle == 0 && iz == 0)
       {
        cout<<"for iz: "<<iz<<" f0: "<<f0<<" f1: "<<f1<<" f2: "<<f2<<" f3: "<<f3<<endl;
      }

      for(int ib=0; ib<2; ib++)
      {
        for(int io = 0; io < NORD; io++)
        {
          ifs>>bbc_angle_four[iz][iangle][ib][io];
        }
      }
    }
  }
  ifs.close();

}

  //TFile *outfile = TFile::Open(Form("processed/%d.root",runnumber),"RECREATE");
  //TFile *outfile = TFile::Open(Form("/phenix/u/theok/hhj2/taxi/Run14HeAu200MinBias/6321/processed/%d.root",runnumber),"RECREATE");
  //TFile *outfile = TFile::Open(Form("~/plhf/taxi/Run14HeAu200MinBias/6506/processed/%d.root",runnumber),"RECREATE");
TFile *outfile = TFile::Open(Form("/gpfs/mnt/gpfs02/phenix/plhf/plhf1/theok/taxi/Run15pAu200CAMBPro104/8650/processed/%d.root",runnumber),"RECREATE");
  //TFile *outfile = TFile::Open("/phenix/u/theok/test/VTX_event_plane/out/fvtx_out.root","RECREATE");

TH2D     *qx[NHAR][NDET];
TH2D     *qy[NHAR][NDET];
TProfile *ave[NZPS][NHAR][NDET];
TProfile *flt[NZPS][NHAR][NDET];
TH2D     *dis[NHAR][NDET];

TH2D     *bbc_angle_qx[n_bbc_angles];
TH2D     *bbc_angle_qy[n_bbc_angles];
TProfile *bbc_angle_ave[NZPS][n_bbc_angles];
TProfile *bbc_angle_flt[NZPS][n_bbc_angles];
TH2D     *bbc_angle_dis[n_bbc_angles];

char name[200];

for(int iz = 0; iz < NZPS; iz++)
{
  for(int iangle = 0; iangle < n_bbc_angles; iangle++)
  {
    sprintf(name,"bbc_angle_ave_%d_%d",iz,iangle);
    bbc_angle_ave[iz][iangle] = new TProfile(name,name,4,-0.5,3.5,-10.1,10.1,"S");
    sprintf(name,"bbc_angle_flt_%d_%d",iz,iangle);
    bbc_angle_flt[iz][iangle] = new TProfile(name,name,4*NORD,-0.5,NORD*4.0-0.5,-1.1,1.1);
  }
}

for (int iz=0; iz<NZPS; iz++) {
  for (int ih=0; ih<NHAR; ih++) {
    for (int id=0; id<NDET; id++) {
     sprintf(name,"ave_%d_%d_%d",iz,ih,id);
       ave[iz][ih][id] = new TProfile(name,name,4,-0.5,3.5,-10.1,10.1,"S");//for SMD -1.1,1.1

       sprintf(name,"flt_%d_%d_%d",iz,ih,id);
       flt[iz][ih][id] = new TProfile(name,name,4*NORD,-0.5,NORD*4.0-0.5,-1.1,1.1);

     }
   }
 }

 for(int iangle = 0; iangle < n_bbc_angles; iangle++)
 {
  sprintf(name,"bbc_angle_dis_%d",iangle);
  bbc_angle_dis[iangle] = new TH2D(name,name,NZPS*3,-0.5,NZPS*3.0-0.5,50,-pi,pi);

  sprintf(name,"bbc_angle_qx_%d",iangle);
  bbc_angle_qx[iangle] = new TH2D(name,name,NZPS*3,-0.5,NZPS*3.0-0.5, 220,-4.1,4.1);

  sprintf(name,"bbc_angle_qy_%d",iangle);
  bbc_angle_qy[iangle] = new TH2D(name,name,NZPS*3,-0.5,NZPS*3.0-0.5, 220,-4.1,4.1);

}
for (int ih=0; ih<NHAR; ih++) {
  for (int id=0; id<NDET; id++) {
    sprintf(name,"dis_%d_%d",ih,id);
    dis[ih][id] = new TH2D(name,name,NZPS*3,-0.5,NZPS*3.0-0.5,50,-pi,pi);

    sprintf(name,"qx_%d_%d",ih,id);
    qx[ih][id] = new TH2D(name,name,NZPS*3,-0.5,NZPS*3.0-0.5, 220,-4.1,4.1);

    sprintf(name,"qy_%d_%d",ih,id);
    qy[ih][id] = new TH2D(name,name,NZPS*3,-0.5,NZPS*3.0-0.5, 220,-4.1,4.1);
  }
}


  //TH1D *resolution_hist = new TH1D("resolution_hist","resolution_hist",3,0,3);

TH1D *eta_east_top_dist = new TH1D("eta_east_top_dist","eta_east_top_dist",70,-3.5,3.5);
TH1D *eta_east_bottom_dist = new TH1D("eta_east_bottom_dist","eta_east_bottom_dist",70,-3.5,3.5);

TH1D *eta_west_top_dist = new TH1D("eta_west_top_dist","eta_west_top_dist",70,-3.5,3.5);
TH1D *eta_west_bottom_dist = new TH1D("eta_west_bottom_dist","eta_west_bottom_dist",70,-3.5,3.5);

TH1D *vtx_phi_dist = new TH1D("vtx_phi_dist","vtx_phi_dist",50,-3.14159,3.14159);

TH1D  *eta_east_dist = new TH1D("eta_east_dist","eta_east_dist",80,-2.0,2.0);
TH1D  *eta_west_dist = new TH1D("eta_west_dist","eta_west_dist",80,-2.0,2.0);

TH1D *eta_fvtx_dist  = new TH1D("eta_fvtx_dist","eta_fvtx_dist",140,-3.5,3.5);
TH1D *eta_tot_dist  = new TH1D("eta_tot_dist","eta_tot_dist",140,-3.5,3.5);

TH1D *forward_zvtx = new TH1D("forward_zvtx","forward_zvtx",150,-30,30);
TH1D *backward_zvtx = new TH1D("backward_zvtx","backward_zvtx",150,-30,30);
TH1D *tot_zvtx = new TH1D("tot_zvtx","tot_zvtx",150,-30,30);

  //TH1D *ngoodtracks_per_eve = new TH1D("ngoodtracks_per_eve","ngoodtracks_per_eve",0,);


TH2D *hphi_zvtx = new TH2D("hphi_zvtx","hphi_zvtx",50,-3.14159/2.0,1.5*3.14159,40,-2,2);
TH2D *hphi_theta = new TH2D("hphi_theta","hphi_theta",50,-3.14159/2.0,1.5*3.14159,50,0,3.14159);


TProfile *phi_v2_fvtx = new TProfile("phi_v2_fvtx","phi_v2_fvtx",24,-pi/2,3*pi/2);
TProfile *phi_v3_fvtx = new TProfile("phi_v3_fvtx","phi_v3_fvtx",24,-pi/2,3*pi/2);
TProfile *phi_v2_mpc = new TProfile("phi_v2_mpc","phi_v2_mpc",24,-pi/2,3*pi/2);
TProfile *phi_v3_mpc = new TProfile("phi_v3_mpc","phi_v3_mpc",24,-pi/2,3*pi/2);
TProfile *phi_v2_bbc = new TProfile("phi_v2_bbc","phi_v2_bbc",24,-pi/2,3*pi/2);
TProfile *phi_v3_bbc = new TProfile("phi_v3_bbc","phi_v3_bbc",24,-pi/2,3*pi/2);

vector< TH1D* > dummy1;
vector< TH1D* > dummy2;
vector <vector<TH1D*> > q_vec_bbc_angle_bf;
vector <vector<TH1D*> > q_vec_bbc_angle_af;

vector<TH1D*> q_vec_fvtx_v2_bf;
vector<TH1D*> q_vec_fvtx_v2_af;

vector<TH1D*> q_vec_vtx_v2_bf;
vector<TH1D*> q_vec_vtx_v2_af;

vector<TH1D*> q_vec_mpc_v2_bf;
vector<TH1D*> q_vec_mpc_v2_af;

vector<TH1D*> q_vec_fvtx_v3_bf;
vector<TH1D*> q_vec_fvtx_v3_af;

vector<TH1D*> q_vec_vtx_v3_bf;
vector<TH1D*> q_vec_vtx_v3_af;

vector<TH1D*> q_vec_mpc_v3_bf;
vector<TH1D*> q_vec_mpc_v3_af;

vector<TH1D*> q_vec_bbc_v3_bf;
vector<TH1D*> q_vec_bbc_v3_af;

vector<TH1D*> q_vec_bbc_v2_bf;
vector<TH1D*> q_vec_bbc_v2_af;

for(int iangle = 0; iangle < n_bbc_angles; iangle++)
{
  for(int izbin = 0; izbin < NZPS; izbin++)
  {
    dummy1.push_back(new TH1D(Form("q_vec_bbc_angle_bf_angle%d_zbin_%d",iangle,izbin),Form("q_vec_bbc_angle_bf_angle%d_zbin_%d",iangle,izbin),50,-3.14159,3.14159));
    dummy2.push_back(new TH1D(Form("q_vec_bbc_angle_af_angle%d_zbin_%d",iangle,izbin),Form("q_vec_bbc_angle_af_angle%d_zbin_%d",iangle,izbin),50,-3.14159,3.14159));
  }
  q_vec_bbc_angle_bf.push_back(dummy1);
  q_vec_bbc_angle_af.push_back(dummy2);
  dummy1.clear();
  dummy2.clear();
}

for(int izbin = 0; izbin < NZPS; izbin++)
{
  q_vec_fvtx_v2_bf.push_back(new TH1D(Form("q_vec_fvtx_v2_bf_%d",izbin),Form("q_vec_fvtx_v2_bf_%d",izbin),50,-3.14159,3.14159));
  q_vec_fvtx_v2_af.push_back(new TH1D(Form("q_vec_fvtx_v2_af_%d",izbin),Form("q_vec_fvtx_v2_af_%d",izbin),50,-3.14159,3.14159));
  q_vec_vtx_v2_bf.push_back(new TH1D(Form("q_vec_vtx_v2_bf_%d",izbin),Form("q_vec_vtx_v2_bf_%d",izbin),50,-3.14159,3.14159));
  q_vec_vtx_v2_af.push_back(new TH1D(Form("q_vec_vtx_v2_af_%d",izbin),Form("q_vec_vtx_v2_af_%d",izbin),50,-3.14159,3.14159));
  q_vec_mpc_v2_bf.push_back(new TH1D(Form("q_vec_mpc_v2_bf_%d",izbin),Form("q_vec_mpc_v2_bf_%d",izbin),50,-3.14159,3.14159));
  q_vec_mpc_v2_af.push_back(new TH1D(Form("q_vec_mpc_v2_af_%d",izbin),Form("q_vec_mpc_v2_af_%d",izbin),50,-3.14159,3.14159));
  q_vec_bbc_v2_bf.push_back(new TH1D(Form("q_vec_bbc_v2_bf_%d",izbin),Form("q_vec_bbc_v2_bf_%d",izbin),50,-3.14159,3.14159));
  q_vec_bbc_v2_af.push_back(new TH1D(Form("q_vec_bbc_v2_af_%d",izbin),Form("q_vec_bbc_v2_af_%d",izbin),50,-3.14159,3.14159));

  q_vec_fvtx_v3_bf.push_back(new TH1D(Form("q_vec_fvtx_v3_bf_%d",izbin),Form("q_vec_fvtx_v3_bf_%d",izbin),50,-3.14159,3.14159));
  q_vec_fvtx_v3_af.push_back(new TH1D(Form("q_vec_fvtx_v3_af_%d",izbin),Form("q_vec_fvtx_v3_af_%d",izbin),50,-3.14159,3.14159));
  q_vec_vtx_v3_bf.push_back(new TH1D(Form("q_vec_vtx_v3_bf_%d",izbin),Form("q_vec_vtx_v3_bf_%d",izbin),50,-3.14159,3.14159));
  q_vec_vtx_v3_af.push_back(new TH1D(Form("q_vec_vtx_v3_af_%d",izbin),Form("q_vec_vtx_v3_af_%d",izbin),50,-3.14159,3.14159));
  q_vec_mpc_v3_bf.push_back(new TH1D(Form("q_vec_mpc_v3_bf_%d",izbin),Form("q_vec_mpc_v3_bf_%d",izbin),50,-3.14159,3.14159));
  q_vec_mpc_v3_af.push_back(new TH1D(Form("q_vec_mpc_v3_af_%d",izbin),Form("q_vec_mpc_v3_af_%d",izbin),50,-3.14159,3.14159));
  q_vec_bbc_v3_bf.push_back(new TH1D(Form("q_vec_bbc_v3_bf_%d",izbin),Form("q_vec_bbc_v3_bf_%d",izbin),50,-3.14159,3.14159));
  q_vec_bbc_v3_af.push_back(new TH1D(Form("q_vec_bbc_v3_af_%d",izbin),Form("q_vec_bbc_v3_af_%d",izbin),50,-3.14159,3.14159));



}

vector<TProfile*> ptbin_phi_v2_fvtx;
vector<TProfile*> ptbin_phi_v3_fvtx;
vector<TProfile*> ptbin_phi_v2_mpc;
vector<TProfile*> ptbin_phi_v3_mpc;
vector<TProfile*> ptbin_phi_v2_bbc;
vector<TProfile*> ptbin_phi_v3_bbc;

for(int i = 0; i < 6; i++)
{
  ptbin_phi_v2_fvtx.push_back(new TProfile(Form("ptbin_phi_v2_fvtx_%d",i),Form("ptbin_phi_v2_fvtx_%d",i),24,-pi/2,3*pi/2));
  ptbin_phi_v3_fvtx.push_back(new TProfile(Form("ptbin_phi_v3_fvtx_%d",i),Form("ptbin_phi_v3_fvtx_%d",i),24,-pi/2,3*pi/2));
  ptbin_phi_v2_mpc.push_back(new TProfile(Form("ptbin_phi_v2_mpc_%d",i),Form("ptbin_phi_v2_mpc_%d",i),24,-pi/2,3*pi/2));
  ptbin_phi_v3_mpc.push_back(new TProfile(Form("ptbin_phi_v3_mpc_%d",i),Form("ptbin_phi_v3_mpc_%d",i),24,-pi/2,3*pi/2));
  ptbin_phi_v2_bbc.push_back(new TProfile(Form("ptbin_phi_v2_bbc_%d",i),Form("ptbin_phi_v2_bbc_%d",i),24,-pi/2,3*pi/2));
  ptbin_phi_v3_bbc.push_back(new TProfile(Form("ptbin_phi_v3_bbc_%d",i),Form("ptbin_phi_v3_bbc_%d",i),24,-pi/2,3*pi/2));

}

  //TProfile *bfcali_ep_resolution = new TProfile("bfcali_ep_resolution","bfcali_ep_resolution",NZPS,0,NZPS);
TProfile *cos_bbc_vtx_psi2 = new TProfile("cos_bbc_vtx_psi2","cos_bbc_vtx_psi2",NZPS,0,NZPS,-10,10);
TProfile *cos_bbc_fvtx_psi2 = new TProfile("cos_bbc_fvtx_psi2","cos_bbc_fvtx_psi2",NZPS,0,NZPS,-10,10);
TProfile *cos_vtx_fvtx_psi2 = new TProfile("cos_vtx_fvtx_psi2","cos_vtx_fvtx_psi2",NZPS,0,NZPS,-10,10);
TProfile *cos_mpc_fvtx_psi2 = new TProfile("cos_mpc_fvtx_psi2","cos_mpc_fvtx_psi2",NZPS,0,NZPS,-10,10);
TProfile *cos_mpc_vtx_psi2 = new TProfile("cos_mpc_vtx_psi2","cos_mpc_vtx_psi2",NZPS,0,NZPS,-10,10);
TProfile *cos_mpc_bbc_psi2 = new TProfile("cos_mpc_bbc_psi2","cos_mpc_bbc_psi2",NZPS,0,NZPS,-10,10);

TProfile *cos_bbc_vtx_psi3 = new TProfile("cos_bbc_vtx_psi3","cos_bbc_vtx_psi3",NZPS,0,NZPS,-10,10);
TProfile *cos_bbc_fvtx_psi3 = new TProfile("cos_bbc_fvtx_psi3","cos_bbc_fvtx_psi3",NZPS,0,NZPS,-10,10);
TProfile *cos_vtx_fvtx_psi3 = new TProfile("cos_vtx_fvtx_psi3","cos_vtx_fvtx_psi3",NZPS,0,NZPS,-10,10);
TProfile *cos_mpc_fvtx_psi3 = new TProfile("cos_mpc_fvtx_psi3","cos_mpc_fvtx_psi3",NZPS,0,NZPS,-10,10);
TProfile *cos_mpc_vtx_psi3 = new TProfile("cos_mpc_vtx_psi3","cos_mpc_vtx_psi3",NZPS,0,NZPS,-10,10);
TProfile *cos_mpc_bbc_psi3 = new TProfile("cos_mpc_bbc_psi3","cos_mpc_bbc_psi3",NZPS,0,NZPS,-10,10);


TH1D* q_vec_2_mpc_mb = new TH1D("q_vec_2_mpc_mb","q_vec_2_mpc_mb",50,-3.14159,3.14159);
TH1D* q_vec_2_mpc_ct = new TH1D("q_vec_2_mpc_ct","q_vec_2_mpc_ct",50,-3.14159,3.14159);
TH1D* q_vec_3_mpc_mb = new TH1D("q_vec_3_mpc_mb","q_vec_3_mpc_mb",50,-3.14159/2,3.14159/2);
TH1D* q_vec_3_mpc_ct = new TH1D("q_vec_3_mpc_ct","q_vec_3_mpc_ct",50,-3.14159/2,3.14159/2);

TH1D* q_vec_2_fvtx_mb = new TH1D("q_vec_2_fvtx_mb","q_vec_2_fvtx_mb",50,-3.14159,3.14159);
TH1D* q_vec_2_fvtx_ct = new TH1D("q_vec_2_fvtx_ct","q_vec_2_fvtx_ct",50,-3.14159,3.14159);
TH1D* q_vec_3_fvtx_mb = new TH1D("q_vec_3_fvtx_mb","q_vec_3_fvtx_mb",50,-3.14159/2,3.14159/2);
TH1D* q_vec_3_fvtx_ct = new TH1D("q_vec_3_fvtx_ct","q_vec_3_fvtx_ct",50,-3.14159/2,3.14159/2);

TProfile *hadronmpcv2 = new TProfile("hadronmpcv2","hadronmpcv2", 70, 0.0, 7.0, -1.1, 1.1);
TProfile *hadronmpcv3 =  new TProfile("hadronmpcv3","hadronmpcv3", 70, 0.0, 7.0, -1.1, 1.1);

TProfile *hadronbbcv2 = new TProfile("hadronbbcv2","hadronbbcv2", 70, 0.0, 7.0, -1.1, 1.1);
TProfile *hadronbbcv3 =  new TProfile("hadronbbcv3","hadronbbcv3", 70, 0.0, 7.0, -1.1, 1.1);

TProfile *hadronfvtxv2 = new TProfile("hadronfvtxv2","hadronfvtxv2", 70, 0.0, 7.0, -1.1, 1.1);
TProfile *hadronfvtxv3 =  new TProfile("hadronfvtxv3","hadronfvtxv3", 70, 0.0, 7.0, -1.1, 1.1);

TProfile *hadronmpcv2_sqr = new TProfile("hadronmpcv2_sqr","hadronmpcv2_sqr", 70, 0.0, 7.0, -1.1, 1.1);
TProfile *hadronmpcv3_sqr =  new TProfile("hadronmpcv3_sqr","hadronmpcv3_sqr", 70, 0.0, 7.0, -1.1, 1.1);

TProfile *hadronbbcv2_sqr = new TProfile("hadronbbcv2_sqr","hadronbbcv2_sqr", 70, 0.0, 7.0, -1.1, 1.1);
TProfile *hadronbbcv3_sqr =  new TProfile("hadronbbcv3_sqr","hadronbbcv3_sqr", 70, 0.0, 7.0, -1.1, 1.1);

TProfile *hadronfvtxv2_sqr = new TProfile("hadronfvtxv2_sqr","hadronfvtxv2_sqr", 70, 0.0, 7.0, -1.1, 1.1);
TProfile *hadronfvtxv3_sqr =  new TProfile("hadronfvtxv3_sqr","hadronfvtxv3_sqr", 70, 0.0, 7.0, -1.1, 1.1);

TProfile *eta_hadronmpcv2 = new TProfile("eta_hadronmpcv2","eta_hadronmpcv2", 28, -3.5, 3.5, -1.1, 1.1);
TProfile *eta_hadronmpcv3 =  new TProfile("eta_hadronmpcv3","eta_hadronmpcv3", 28, -3.5, 3.5, -1.1, 1.1);

TProfile *eta_hadronbbcv2 = new TProfile("eta_hadronbbcv2","eta_hadronbbcv2", 28, -3.5, 3.5, -1.1, 1.1);
TProfile *eta_hadronbbcv3 =  new TProfile("eta_hadronbbcv3","eta_hadronbbcv3", 28, -3.5, 3.5, -1.1, 1.1);

TProfile *eta_hadronfvtxv2 = new TProfile("eta_hadronfvtxv2","eta_hadronfvtxv2", 28, -3.5, 3.5, -1.1, 1.1);
TProfile *eta_hadronfvtxv3 =  new TProfile("eta_hadronfvtxv3","eta_hadronfvtxv3", 28, -3.5, 3.5, -1.1, 1.1);

TProfile *feta_hadronmpcv2 = new TProfile("feta_hadronmpcv2","feta_hadronmpcv2", 28, -3.5, 3.5, -1.1, 1.1);
TProfile *feta_hadronmpcv3 =  new TProfile("feta_hadronmpcv3","feta_hadronmpcv3", 28, -3.5, 3.5, -1.1, 1.1);

TProfile *feta_hadronbbcv2 = new TProfile("feta_hadronbbcv2","feta_hadronbbcv2", 28, -3.5, 3.5, -1.1, 1.1);
TProfile *feta_hadronbbcv3 =  new TProfile("feta_hadronbbcv3","feta_hadronbbcv3", 28, -3.5, 3.5, -1.1, 1.1);

TProfile *feta_hadronfvtxv2 = new TProfile("feta_hadronfvtxv2","feta_hadronfvtxv2", 28, -3.5, 3.5, -1.1, 1.1);
TProfile *feta_hadronfvtxv3 =  new TProfile("feta_hadronfvtxv3","feta_hadronfvtxv3", 28, -3.5, 3.5, -1.1, 1.1);

vector<TProfile*> hadronfvtxv2_eta;
vector<TProfile*> hadronmpcv2_eta;
vector<TProfile*> hadronbbcv2_eta;

vector<TProfile*> hadronfvtxv2_eta_sqr;
vector<TProfile*> hadronmpcv2_eta_sqr;
vector<TProfile*> hadronbbcv2_eta_sqr;


vector<TProfile*> shadronfvtxv2_eta;
vector<TProfile*> shadronmpcv2_eta;
vector<TProfile*> shadronbbcv2_eta;

for(int i = 0; i < NETA; i++)
{
  hadronfvtxv2_eta.push_back(new TProfile(Form("hadronfvtxv2_eta%d",i),Form("hadronfvtxv2_eta%d",i), 70, 0.0, 7.0, -1.1, 1.1));
  hadronmpcv2_eta.push_back(new TProfile(Form("hadronmpcv2_eta%d",i),Form("hadronmpcv2_eta%d",i), 70, 0.0, 7.0, -1.1, 1.1));
  hadronbbcv2_eta.push_back(new TProfile(Form("hadronbbcv2_eta%d",i),Form("hadronbbcv2_eta%d",i), 70, 0.0, 7.0, -1.1, 1.1));
  shadronfvtxv2_eta.push_back(new TProfile(Form("shadronfvtxv2_eta%d",i),Form("shadronfvtxv2_eta%d",i), 70, 0.0, 7.0, -1.1, 1.1));
  shadronmpcv2_eta.push_back(new TProfile(Form("shadronmpcv2_eta%d",i),Form("shadronmpcv2_eta%d",i), 70, 0.0, 7.0, -1.1, 1.1));
  shadronbbcv2_eta.push_back(new TProfile(Form("shadronbbcv2_eta%d",i),Form("shadronbbcv2_eta%d",i), 70, 0.0, 7.0, -1.1, 1.1));

  hadronfvtxv2_eta_sqr.push_back(new TProfile(Form("hadronfvtxv2_eta%d_sqr",i),Form("hadronfvtxv2_eta%d_sqr",i), 70, 0.0, 7.0, -1.1, 1.1));
  hadronmpcv2_eta_sqr.push_back(new TProfile(Form("hadronmpcv2_eta%d_sqr",i),Form("hadronmpcv2_eta%d_sqr",i), 70, 0.0, 7.0, -1.1, 1.1));
  hadronbbcv2_eta_sqr.push_back(new TProfile(Form("hadronbbcv2_eta%d_sqr",i),Form("hadronbbcv2_eta%d_sqr",i), 70, 0.0, 7.0, -1.1, 1.1));
}

vector<TProfile *> hadronfvtxv3_eta;
vector<TProfile *> hadronmpcv3_eta;
vector<TProfile *> hadronbbcv3_eta;
vector<TProfile *> shadronfvtxv3_eta;
vector<TProfile *> shadronbbcv3_eta;
vector<TProfile *> shadronmpcv3_eta;
for(int i = 0; i < NETA; i++)
{
  hadronfvtxv3_eta.push_back(new TProfile(Form("hadronfvtxv3_eta%d",i),Form("hadronfvtxv3_eta%d",i), 70, 0.0, 7.0, -1.1, 1.1));
  hadronmpcv3_eta.push_back(new TProfile(Form("hadronmpcv3_eta%d",i),Form("hadronmpcv3_eta%d",i), 70, 0.0, 7.0, -1.1, 1.1));
  hadronbbcv3_eta.push_back(new TProfile(Form("hadronbbcv3_eta%d",i),Form("hadronbbcv3_eta%d",i), 70, 0.0, 7.0, -1.1, 1.1));

  shadronfvtxv3_eta.push_back(new TProfile(Form("shadronfvtxv3_eta%d",i),Form("shadronfvtxv3_eta%d",i), 70, 0.0, 7.0, -1.1, 1.1));
  shadronmpcv3_eta.push_back(new TProfile(Form("shadronmpcv3_eta%d",i),Form("shadronmpcv3_eta%d",i), 70, 0.0, 7.0, -1.1, 1.1));
  shadronbbcv3_eta.push_back(new TProfile(Form("shadronbbcv3_eta%d",i),Form("shadronbbcv3_eta%d",i), 70, 0.0, 7.0, -1.1, 1.1));


}

vector<TH1D*> vtx_pt_eta;

for(int i = 0; i < NETA; i++)
{
  vtx_pt_eta.push_back(new TH1D(Form("vtx_pt_eta_%d",i),Form("vtx_pt_eta_%d",i),100,0.0,5.0));
}

vector<TProfile * > east_arm_bbcv2_angle;
vector<TProfile * > west_arm_bbcv2_angle;
vector<TProfile * > both_arm_bbcv2_angle;

for(int iangle = 0; iangle < n_bbc_angles; iangle++)
{
  east_arm_bbcv2_angle.push_back(new TProfile(Form("east_arm_bbcv2_angle%d",iangle),Form("east_arm_bbcv2_angle%d",iangle),70,0.0,7.0,-1.1,1.1));
  west_arm_bbcv2_angle.push_back(new TProfile(Form("west_arm_bbcv2_angle%d",iangle),Form("west_arm_bbcv2_angle%d",iangle),70,0.0,7.0,-1.1,1.1));
  both_arm_bbcv2_angle.push_back(new TProfile(Form("both_arm_bbcv2_angle%d",iangle),Form("both_arm_bbcv2_angle%d",iangle),70,0.0,7.0,-1.1,1.1));

}

TProfile *east_arm_hadronmpcv2 = new TProfile("east_arm_hadronmpcv2","east_arm_hadronmpcv2", 70, 0.0, 7.0, -1.1, 1.1);
TProfile *east_arm_hadronmpcv3 =  new TProfile("east_arm_hadronmpcv3","east_arm_hadronmpcv3", 70, 0.0, 7.0, -1.1, 1.1);

TProfile *east_arm_hadronbbcv2 = new TProfile("east_arm_hadronbbcv2","east_arm_hadronbbcv2", 70, 0.0, 7.0, -1.1, 1.1);
TProfile *east_arm_hadronbbcv3 =  new TProfile("east_arm_hadronbbcv3","east_arm_hadronbbcv3", 70, 0.0, 7.0, -1.1, 1.1);

TProfile *east_arm_hadronfvtxv2 = new TProfile("east_arm_hadronfvtxv2","east_arm_hadronfvtxv2", 70, 0.0, 7.0, -1.1, 1.1);
TProfile *east_arm_hadronfvtxv3 =  new TProfile("east_arm_hadronfvtxv3","east_arm_hadronfvtxv3", 70, 0.0, 7.0, -1.1, 1.1);

TProfile *west_arm_hadronmpcv2 = new TProfile("west_arm_hadronmpcv2","west_arm_hadronmpcv2", 70, 0.0, 7.0, -1.1, 1.1);
TProfile *west_arm_hadronmpcv3 =  new TProfile("west_arm_hadronmpcv3","west_arm_hadronmpcv3", 70, 0.0, 7.0, -1.1, 1.1);

TProfile *west_arm_hadronbbcv2 = new TProfile("west_arm_hadronbbcv2","west_arm_hadronbbcv2", 70, 0.0, 7.0, -1.1, 1.1);
TProfile *west_arm_hadronbbcv3 =  new TProfile("west_arm_hadronbbcv3","west_arm_hadronbbcv3", 70, 0.0, 7.0, -1.1, 1.1);

TProfile *west_arm_hadronfvtxv2 = new TProfile("west_arm_hadronfvtxv2","west_arm_hadronfvtxv2", 70, 0.0, 7.0, -1.1, 1.1);
TProfile *west_arm_hadronfvtxv3 =  new TProfile("west_arm_hadronfvtxv3","west_arm_hadronfvtxv3", 70, 0.0, 7.0, -1.1, 1.1);


TProfile *q_weight_mpc_pt = new TProfile("q_weight_mpc_pt","q_weight_mpc_pt",70,0.0,7,0.0,1200);
TProfile *q_weight_fvtx_pt = new TProfile("q_weight_fvtx_pt","q_weight_fvtx_pt",70,0.0,7,0.0,1200);
TProfile *q_weight_bbc_pt = new TProfile("q_weight_bbc_pt","q_weight_bbc_pt",70,0.0,7,0.0,1200);


  // sin hadron distributions
TProfile *shadronmpcv2 = new TProfile("shadronmpcv2","shadronmpcv2", 70, 0.0, 7.0, -1.1, 1.1);
TProfile *shadronmpcv3 =  new TProfile("shadronmpcv3","shadronmpcv3", 70, 0.0, 7.0, -1.1, 1.1);

TProfile *shadronbbcv2 = new TProfile("shadronbbcv2","shadronbbcv2", 70, 0.0, 7.0, -1.1, 1.1);
TProfile *shadronbbcv3 =  new TProfile("shadronbbcv3","shadronbbcv3", 70, 0.0, 7.0, -1.1, 1.1);

TProfile *shadronfvtxv2 = new TProfile("shadronfvtxv2","shadronfvtxv2", 70, 0.0, 7.0, -1.1, 1.1);
TProfile *shadronfvtxv3 =  new TProfile("shadronfvtxv3","shadronfvtxv3", 70, 0.0, 7.0, -1.1, 1.1);

TProfile *eta_shadronmpcv2 = new TProfile("eta_shadronmpcv2","eta_shadronmpcv2", 28, -3.5, 3.5, -1.1, 1.1);
TProfile *eta_shadronmpcv3 =  new TProfile("eta_shadronmpcv3","eta_shadronmpcv3", 28, -3.5, 3.5, -1.1, 1.1);

TProfile *eta_shadronbbcv2 = new TProfile("eta_shadronbbcv2","eta_shadronbbcv2", 28, -3.5, 3.5, -1.1, 1.1);
TProfile *eta_shadronbbcv3 =  new TProfile("eta_shadronbbcv3","eta_shadronbbcv3", 28, -3.5, 3.5, -1.1, 1.1);

TProfile *eta_shadronfvtxv2 = new TProfile("eta_shadronfvtxv2","eta_shadronfvtxv2", 28, -3.5, 3.5, -1.1, 1.1);
TProfile *eta_shadronfvtxv3 =  new TProfile("eta_shadronfvtxv3","eta_shadronfvtxv3", 28, -3.5, 3.5, -1.1, 1.1);

TProfile *feta_shadronmpcv2 = new TProfile("feta_shadronmpcv2","feta_shadronmpcv2", 28, -3.5, 3.5, -1.1, 1.1);
TProfile *feta_shadronmpcv3 =  new TProfile("feta_shadronmpcv3","feta_shadronmpcv3", 28, -3.5, 3.5, -1.1, 1.1);

TProfile *feta_shadronbbcv2 = new TProfile("feta_shadronbbcv2","feta_shadronbbcv2", 28, -3.5, 3.5, -1.1, 1.1);
TProfile *feta_shadronbbcv3 =  new TProfile("feta_shadronbbcv3","feta_shadronbbcv3", 28, -3.5, 3.5, -1.1, 1.1);

TProfile *feta_shadronfvtxv2 = new TProfile("feta_shadronfvtxv2","feta_shadronfvtxv2", 28, -3.5, 3.5, -1.1, 1.1);
TProfile *feta_shadronfvtxv3 =  new TProfile("feta_shadronfvtxv3","feta_shadronfvtxv3", 28, -3.5, 3.5, -1.1, 1.1);

TProfile *eta_hadronmpcv2_sqr = new TProfile("eta_hadronmpcv2_sqr","eta_hadronmpcv2_sqr", 28, -3.5, 3.5, -1.1, 1.1);
TProfile *eta_hadronbbcv2_sqr = new TProfile("eta_hadronbbcv2_sqr","eta_hadronbbcv2_sqr", 28, -3.5, 3.5, -1.1, 1.1);
TProfile *eta_hadronfvtxv2_sqr = new TProfile("eta_hadronfvtxv2_sqr","eta_hadronfvtxv2_sqr", 28, -3.5, 3.5, -1.1, 1.1);
TProfile *feta_hadronmpcv2_sqr = new TProfile("feta_hadronmpcv2_sqr","feta_hadronmpcv2_sqr", 28, -3.5, 3.5, -1.1, 1.1);
TProfile *feta_hadronbbcv2_sqr = new TProfile("feta_hadronbbcv2_sqr","feta_hadronbbcv2_sqr", 28, -3.5, 3.5, -1.1, 1.1);
TProfile *feta_hadronfvtxv2_sqr = new TProfile("feta_hadronfvtxv2_sqr","feta_hadronfvtxv2_sqr", 28, -3.5, 3.5, -1.1, 1.1);


TH1D *all_vtx_pt = new TH1D("all_vtx_pt","all_vtx_pt",100,0.0,7.0);
TH1D *fvtx_vtx_pt = new TH1D("fvtx_vtx_pt","fvtx_vtx_pt",100,0.0,7.0);
TH1D *bbc_vtx_pt = new TH1D("bbc_vtx_pt","bbc_vtx_pt",100,0.0,7.0);
TH1D *mpc_vtx_pt = new TH1D("mpc_vtx_pt","mpc_vtx_pt",100,0.0,7.0);


int nentries = ntp_event->GetEntries();
cout<<"total events = " << nentries<<endl;
for ( Int_t ievt = 0 ; ievt < nentries ; ievt++ ) {

  if(ievt%100000==0)
    cout<<"procesed "<<ievt*100.0/nentries<<" percentage of events"<<endl;

    //if(ievt>5) break;

    //event variables
  if(!trimmed_tree)
  {
    b_event->GetEntry(ievt);
    b_centrality->GetEntry(ievt);
    b_trigger->GetEntry(ievt);
    b_bbc_qn->GetEntry(ievt);
    b_bbc_qs->GetEntry(ievt);
    b_vtx_x->GetEntry(ievt);
    b_vtx_y->GetEntry(ievt);
    b_vtx_z->GetEntry(ievt);
    b_fvtx_z->GetEntry(ievt);
    b_eventok->GetEntry(ievt);
  }
  b_bbc_z->GetEntry(ievt);

  b_d_Qx->GetEntry(ievt);
  b_d_Qy->GetEntry(ievt);
  b_d_Qw->GetEntry(ievt);

  b_d_BBCs_Q->GetEntry(ievt);

  b_nsegments->GetEntry(ievt);



    //track variables

  b_px->GetEntry(ievt);
  b_py->GetEntry(ievt);
  b_pz->GetEntry(ievt);

  if(!trimmed_tree)
  {
    b_trackID->GetEntry(ievt);
    b_charge->GetEntry(ievt);
    b_chisq->GetEntry(ievt);
    b_ndf->GetEntry(ievt);
    b_nhit0->GetEntry(ievt);
    b_nhit1->GetEntry(ievt);
    b_nhit2->GetEntry(ievt);
    b_nhit3->GetEntry(ievt);
    b_dca->GetEntry(ievt);
    b_dca2d->GetEntry(ievt);
    b_segmentok->GetEntry(ievt);
  }
    //if(include_fvtx_tracks)
    //{
  /*
  b_ntracklets->GetEntry(ievt);
  b_feta->GetEntry(ievt);
  b_fphi->GetEntry(ievt);
  b_fchisq->GetEntry(ievt);
  b_farm->GetEntry(ievt);
  */
    //}
    //cout<<"trigger "<<trigger<<endl;

    //if(trigger!=4) continue;
  if(centrality > 5 && !trimmed_tree) continue;

    //int ibbcz  = NZPS*(vtx_z+10)/20;
  int ibbcz  = NZPS*(bbc_z+30)/60;


  if(ibbcz<0||ibbcz>=NZPS)
  {
   if(ibbcz < 0)
     ibbcz  = 0;
   if(ibbcz >= NZPS)
    ibbcz = NZPS-1;
  continue;
}

    //if(fabs(bbc_z ) > 1.0) continue;

//if(vtx_z!=vtx_z) continue;
//if(fabs(vtx_z) > 10.0) continue;

float sumxy[NHAR][NDET][4];
for (int i=0; i<NHAR; i++) {
  for (int j=0; j<NDET; j++) {
   for (int k=0; k<4; k++) {
     sumxy[i][j][k]=0;
   }
 }
}

float sumxy_bbc_angle[n_bbc_angles][4];
for(int i = 0; i < n_bbc_angles;i++)
  for(int j = 0; j<4; j++)
    sumxy_bbc_angle[i][j]=0;
  bool good_mpc_ep = false;
  bool good_fvtx_ep = false;
  bool good_bbc_ep = false;

  for(int ihar=0; ihar<NHAR; ihar++){
      //MPC
    if(d_Qw[ihar*3+0]>0)
    {
    sumxy[ihar][0][0] = d_Qx[ihar*3+0]; //south
    sumxy[ihar][0][1] = d_Qy[ihar*3+0];
    sumxy[ihar][0][2] = d_Qw[ihar*3+0];

    if(ihar==1) good_mpc_ep = true;

  }
      if(d_Qw[ihar*3+1]<1000 && d_Qw[ihar*3+1]>5)//FVTXs
      {
       sumxy[ihar][6][0] = d_Qx[ihar*3+1];
       sumxy[ihar][6][1] = d_Qy[ihar*3+1];
       sumxy[ihar][6][2] = d_Qw[ihar*3+1];

       if(ihar==1) good_fvtx_ep = true;
     }
      if(d_Qw[ihar*3+2]>0)//bbc
      {
       sumxy[ihar][4][0] = d_Qx[ihar*3+2];
       sumxy[ihar][4][1] = d_Qy[ihar*3+2];
       sumxy[ihar][4][2] = d_Qw[ihar*3+2];

       if(ihar==1)
       {
         for(int iangle = 0; iangle < n_bbc_angles; iangle++)
         {
           for(int icomp = 0; icomp < 3; icomp++)
             sumxy_bbc_angle[iangle][icomp] = d_BBCs_Q[iangle*3 + icomp];
         }
       }


       if(ihar==1) good_bbc_ep = true;
     }
    }//end of ihar


    for(int itrk=0; itrk< nsegments; itrk++){
      if(itrk>=50) break;

      float phi = TMath::ATan2(py[itrk],px[itrk]);

      for(int ihar = 0; ihar<NHAR;ihar++)
      {
       float qx = TMath::Cos(phi*(ihar+1));
       float qy = TMath::Sin(phi*(ihar+1));
       sumxy[ihar][5][0] += qx;
       sumxy[ihar][5][1] += qy;
       sumxy[ihar][5][2] += 1;
     }

     if(itrk>=50) break;
     if(!trimmed_tree)
     {
      if(segmentok[itrk]!=1) continue;
      if(dca[itrk]*dca[itrk] > 0.125 || dca2d[itrk]*dca2d[itrk] > 0.0075) continue;
      if(nhit0[itrk]!=4) continue;
    }
    float pt = TMath::Sqrt(px[itrk]*px[itrk]+py[itrk]*py[itrk]);

    all_vtx_pt->Fill(pt);


    if(good_bbc_ep) bbc_vtx_pt->Fill(pt);
    if(good_fvtx_ep) fvtx_vtx_pt->Fill(pt);
    if(good_mpc_ep) mpc_vtx_pt->Fill(pt);
  }


  for(int iangle = 0; iangle < n_bbc_angles; iangle++)
  {
    float bbc_angle_psi2 = TMath::ATan2(sumxy_bbc_angle[iangle][1],sumxy_bbc_angle[iangle][0])/2.0;
    if(sumxy_bbc_angle[iangle][2] > 0.0)
      q_vec_bbc_angle_bf[iangle][ibbcz]->Fill(bbc_angle_psi2,1);
  }

  double mpc_psi2 = TMath::ATan2(sumxy[1][0][1],sumxy[1][0][0])/2.0;
  double mpc_psi3 = TMath::ATan2(sumxy[2][0][1],sumxy[2][0][0])/3.0;

  double fvtx_psi2 = TMath::ATan2(sumxy[1][6][1],sumxy[1][6][0])/2.0;
  double fvtx_psi3 = TMath::ATan2(sumxy[2][6][1],sumxy[2][6][0])/3.0;

  double bbc_psi2 = TMath::ATan2(sumxy[1][4][1],sumxy[1][4][0])/2.0;
  double bbc_psi3 = TMath::ATan2(sumxy[2][4][1],sumxy[2][4][0])/3.0;

  double vtx_psi2 = TMath::ATan2(sumxy[1][5][1],sumxy[1][5][0])/2.0;
  double vtx_psi3 = TMath::ATan2(sumxy[2][5][1],sumxy[2][5][0])/3.0;

  q_vec_bbc_v2_bf[ibbcz]->Fill(bbc_psi2,1);
  q_vec_bbc_v3_bf[ibbcz]->Fill(bbc_psi3,1);

  if(sumxy[1][0][2]>0.0)
  {
   q_vec_mpc_v2_bf[ibbcz]->Fill(mpc_psi2,1);
   if(trigger!=4 || trimmed_tree){
     q_vec_2_mpc_mb->Fill(mpc_psi2,1);
     if(sumxy[2][0][2]>0.0)
       q_vec_3_mpc_mb->Fill(mpc_psi3,1);
   }
   else{
     q_vec_2_mpc_ct->Fill(mpc_psi2,1);
     if(sumxy[2][0][2]>0.0)
       q_vec_3_mpc_ct->Fill(mpc_psi3,1);
   }
 }
 if(sumxy[1][6][2]>0.0)
 {
   q_vec_fvtx_v2_bf[ibbcz]->Fill(fvtx_psi2,1);
   if(trigger!=4 || trimmed_tree){
     q_vec_2_fvtx_mb->Fill(fvtx_psi2,1);
     if(sumxy[2][6][2]>0.0)
       q_vec_3_fvtx_mb->Fill(fvtx_psi3,1);
   }
   else{
     q_vec_2_fvtx_ct->Fill(fvtx_psi2,1);
     if(sumxy[2][6][2]>0.0)
       q_vec_3_fvtx_ct->Fill(fvtx_psi3,1);
   }
 }

 if(sumxy[2][0][2]>0.0)
 {
   q_vec_mpc_v3_bf[ibbcz]->Fill(mpc_psi3,1);
 }
 if(sumxy[2][6][2]>0.0)
 {
   q_vec_fvtx_v3_bf[ibbcz]->Fill(fvtx_psi3,1);
 }

 q_vec_vtx_v2_bf[ibbcz]->Fill(vtx_psi2,1);
 q_vec_vtx_v3_bf[ibbcz]->Fill(vtx_psi3,1);

   //bbc angle flattening
 for(int iangle = 0; iangle<n_bbc_angles; iangle++)
 {
  if(sumxy_bbc_angle[iangle][2]>0.0)
  {
    sumxy_bbc_angle[iangle][3]=TMath::ATan2(sumxy_bbc_angle[iangle][1],sumxy_bbc_angle[iangle][0])/(2);
       if (calFlag < 3 && calFlag>0) bbc_angle_dis[iangle]->Fill(ibbcz,sumxy_bbc_angle[iangle][3]*(2));//fill the PSI2
     }
     if (sumxy_bbc_angle[iangle][2]>0.0) {

       for (int ib=0; ib<2; ib++) {
         sumxy_bbc_angle[iangle][ib]/=sumxy_bbc_angle[iangle][2];//dividing by weight?

         if (calFlag < 3 && calFlag>0) {
           bbc_angle_ave[ibbcz][iangle]->Fill(ib+0.0,sumxy_bbc_angle[iangle][ib]);//fill the components after normalizing
           if(ib==0) bbc_angle_qx[iangle]->Fill(ibbcz,sumxy_bbc_angle[iangle][0]);
           if(ib==1) bbc_angle_qy[iangle]->Fill(ibbcz,sumxy_bbc_angle[iangle][1]);
         }

         float sxy=sumxy_bbc_angle[iangle][ib];
         float mxy=bbc_angle_mean[ibbcz][iangle][ib];//for calflag 1: = 0
         float wxy=bbc_angle_widt[ibbcz][iangle][ib];// for calflag 1: = 1

         sumxy_bbc_angle[iangle][ib]=(sxy-mxy)/wxy;//event (angle mean - mean )/ width
         if (calFlag < 3 && calFlag>0) {
           bbc_angle_ave[ibbcz][iangle]->Fill(ib+2.0,sumxy_bbc_angle[iangle][ib]);
           if(ib==0) bbc_angle_qx[iangle]->Fill(ibbcz+NZPS,sumxy_bbc_angle[iangle][0]);
           if(ib==1) bbc_angle_qy[iangle]->Fill(ibbcz+NZPS,sumxy_bbc_angle[iangle][1]);
         }

       }

       sumxy_bbc_angle[iangle][3]=TMath::ATan2(sumxy_bbc_angle[iangle][1],sumxy_bbc_angle[iangle][0])/(2);
       if (calFlag < 3 && calFlag>0) {
         bbc_angle_dis[iangle]->Fill(ibbcz+NZPS,sumxy_bbc_angle[iangle][3]*(2));
       }

       float psi=sumxy_bbc_angle[iangle][3]*(2);
       float dp=0.0;
       for (int io=0; io<NORD; io++) {
         float cc=cos((io+1.0)*psi);
         float ss=sin((io+1.0)*psi);
         if (calFlag < 3 && calFlag>0) bbc_angle_flt[ibbcz][iangle]->Fill(io+0.0,cc);
         if (calFlag < 3 && calFlag>0) bbc_angle_flt[ibbcz][iangle]->Fill(io+NORD,ss);
        float aa= bbc_angle_four[ibbcz][iangle][0][io]; // mean cos
        float bb= bbc_angle_four[ibbcz][iangle][1][io]; // mean sin
        dp+=(aa*ss-bb*cc)*2.0/(io+1.0);
      }

      psi+=dp;
      psi=TMath::ATan2(sin(psi),cos(psi));
      for (int io=0; io<NORD; io++) {
       float cc=cos((io+1.0)*psi);
       float ss=sin((io+1.0)*psi);
       if (calFlag < 3 && calFlag>0) bbc_angle_flt[ibbcz][iangle]->Fill(io+NORD*2.0,cc);
       if (calFlag < 3 && calFlag>0) bbc_angle_flt[ibbcz][iangle]->Fill(io+NORD*3.0,ss);
     }
     sumxy_bbc_angle[iangle][3]=psi/(2);
     if(iangle==0 && ievt%500000==0)
     {
      cout<<"for event: "<<ievt<<" match psi/2: "<<psi/(2)<<endl;
    }

    if (calFlag < 3 && calFlag>0) bbc_angle_dis[iangle]->Fill(ibbcz+NZPS*2.0,sumxy_bbc_angle[iangle][3]*(2));

  }
  else sumxy_bbc_angle[iangle][3]=-9999.9;
}

for (int ih=0; ih<NHAR; ih++) {
  for(int id=0;id<NDET;id++){
   if (sumxy[ih][id][2]>0.0) {
     sumxy[ih][id][3]=TMath::ATan2(sumxy[ih][id][1],sumxy[ih][id][0])/(ih+1.0);
     if (calFlag < 3 && calFlag>0) dis[ih][id]->Fill(ibbcz,sumxy[ih][id][3]*(ih+1.0));
   }
   if (sumxy[ih][id][2]>0.0) {
     for (int ib=0; ib<2; ib++) {
       sumxy[ih][id][ib]/=sumxy[ih][id][2];
      //if(ih==1 && id==0 && ib==0 && sumxy[ih][id][ib]>1) cout<<sumxy[ih][id][ib]<<endl;
       if (calFlag < 3 && calFlag>0) {
         ave[ibbcz][ih][id]->Fill(ib+0.0,sumxy[ih][id][ib]);
         if(ib==0) qx[ih][id]->Fill(ibbcz,sumxy[ih][id][0]);
         if(ib==1) qy[ih][id]->Fill(ibbcz,sumxy[ih][id][1]);
       }

       float sxy=sumxy[ih][id][ib];
       float mxy=mean[0][ibbcz][ih][id][ib];
       float wxy=widt[0][ibbcz][ih][id][ib];

      //if(icent==0 && ibbcz==0 && ih==1 && id==0) cout<<ib<<" "<<sxy<<" "<<mxy<<" "<<wxy<<endl;
       sumxy[ih][id][ib]=(sxy-mxy)/wxy;
       if (calFlag < 3 && calFlag>0) {
         ave[ibbcz][ih][id]->Fill(ib+2.0,sumxy[ih][id][ib]);
         if(ib==0) qx[ih][id]->Fill(ibbcz+NZPS,sumxy[ih][id][0]);
         if(ib==1) qy[ih][id]->Fill(ibbcz+NZPS,sumxy[ih][id][1]);
       }

     }

     sumxy[ih][id][3]=TMath::ATan2(sumxy[ih][id][1],sumxy[ih][id][0])/(ih+1.0);
     if (calFlag < 3 && calFlag>0) {
       dis[ih][id]->Fill(ibbcz+NZPS,sumxy[ih][id][3]*(ih+1.0));
     }

     float psi=sumxy[ih][id][3]*(ih+1.0);
     float dp=0.0;
     for (int io=0; io<NORD; io++) {
       float cc=cos((io+1.0)*psi);
       float ss=sin((io+1.0)*psi);
       if (calFlag < 3 && calFlag>0) flt[ibbcz][ih][id]->Fill(io+0.0,cc);
       if (calFlag < 3 && calFlag>0) flt[ibbcz][ih][id]->Fill(io+NORD,ss);
      float aa=four[0][ibbcz][ih][id][0][io]; // mean cos
      float bb=four[0][ibbcz][ih][id][1][io]; // mean sin
      dp+=(aa*ss-bb*cc)*2.0/(io+1.0);
    }
    psi+=dp;
    psi=TMath::ATan2(sin(psi),cos(psi));
    for (int io=0; io<NORD; io++) {
     float cc=cos((io+1.0)*psi);
     float ss=sin((io+1.0)*psi);
     if (calFlag < 3 && calFlag>0) flt[ibbcz][ih][id]->Fill(io+NORD*2.0,cc);
     if (calFlag < 3 && calFlag>0) flt[ibbcz][ih][id]->Fill(io+NORD*3.0,ss);
   }
   sumxy[ih][id][3]=psi/(ih+1.0);
   if(ih==1 && id==4 && ievt%500000==0)
   {
     cout<<"for event: "<<ievt<<" correct psi/2: "<<psi/(2)<<endl;
   }
   if (calFlag < 3 && calFlag>0) dis[ih][id]->Fill(ibbcz+NZPS*2.0,sumxy[ih][id][3]*(ih+1.0));
 }
 else sumxy[ih][id][3]=-9999.9;

      }//end of id

    }  //end of ih

    if(calFlag <3) continue;


    float         mpcsRP2  = (sumxy[1][0][2]>0)?sumxy[1][0][3]:-9999.9;//south
    float         fvtx0sRP2  = (sumxy[1][6][2]>0)?sumxy[1][6][3]:-9999.9;//south
    float         vtx0sRP2  = (sumxy[1][5][2]>0)?sumxy[1][5][3]:-9999.9;//south
    float         bbcsRP2  = (sumxy[1][4][2]>0)?sumxy[1][4][3]:-9999.9;//south

    //PSi3
    float         mpcsRP3  = (sumxy[2][0][2]>0)?sumxy[2][0][3]:-9999.9;//south
    float         fvtx0sRP3  = (sumxy[2][6][2]>0)?sumxy[2][6][3]:-9999.9;//south
    float         vtx0sRP3  = (sumxy[2][5][2]>0)?sumxy[2][5][3]:-9999.9;//south
    float         bbcsRP3  = (sumxy[2][4][2]>0)?sumxy[2][4][3]:-9999.9;//south
    //cout<<"nsegments: "<<nsegments<<endl;

    float bbcsRP2_angle[n_bbc_angles];
    for(int iangle = 0; iangle < n_bbc_angles; iangle++)
    {
      bbcsRP2_angle[iangle] = (sumxy_bbc_angle[iangle][2]>0)?sumxy_bbc_angle[iangle][3]:-9999.9;
    }

    q_vec_mpc_v2_af[ibbcz]->Fill(mpcsRP2);
    q_vec_fvtx_v2_af[ibbcz]->Fill(fvtx0sRP2);
    q_vec_vtx_v2_af[ibbcz]->Fill(vtx0sRP2);
    q_vec_bbc_v2_af[ibbcz]->Fill(bbcsRP2);

    q_vec_mpc_v3_af[ibbcz]->Fill(mpcsRP3);
    q_vec_fvtx_v3_af[ibbcz]->Fill(fvtx0sRP3);
    q_vec_vtx_v3_af[ibbcz]->Fill(vtx0sRP3);
    q_vec_bbc_v3_af[ibbcz]->Fill(bbcsRP3);

    for(int iangle = 0; iangle < n_bbc_angles; iangle++)
    {
      float bbc_angle_psi2_af = (sumxy_bbc_angle[iangle][2]>0)?sumxy_bbc_angle[iangle][3]:-9999.9;
      q_vec_bbc_angle_af[iangle][ibbcz]->Fill(bbc_angle_psi2_af,1);
    }


    //if(mpcsRP2 < -1000 || fvtx0sRP2 < -1000 || mpcsRP3 < -1000 || fvtx0sRP3< -1000 || bbcsRP2 < -1000 || bbcsRP3 < -1000 || vtx0sRP2 < -1000 || vtx0sRP3 < -1000) continue;


    cos_bbc_vtx_psi2->Fill(ibbcz,TMath::Cos(2*(bbcsRP2-vtx0sRP2)));
    cos_bbc_fvtx_psi2->Fill(ibbcz,TMath::Cos(2*(bbcsRP2-fvtx0sRP2)));
    cos_vtx_fvtx_psi2->Fill(ibbcz,TMath::Cos(2*(vtx0sRP2-fvtx0sRP2)));
    cos_mpc_fvtx_psi2->Fill(ibbcz,TMath::Cos(2*(mpcsRP2-fvtx0sRP2)));
    cos_mpc_vtx_psi2->Fill(ibbcz,TMath::Cos(2*(mpcsRP2-vtx0sRP2)));
    cos_mpc_bbc_psi2->Fill(ibbcz,TMath::Cos(2*(mpcsRP2-bbcsRP2)));

    cos_bbc_vtx_psi3->Fill(ibbcz,TMath::Cos(3*(bbcsRP3-vtx0sRP3)));
    cos_bbc_fvtx_psi3->Fill(ibbcz,TMath::Cos(3*(bbcsRP3-fvtx0sRP3)));
    cos_vtx_fvtx_psi3->Fill(ibbcz,TMath::Cos(3*(vtx0sRP3-fvtx0sRP3)));
    cos_mpc_fvtx_psi3->Fill(ibbcz,TMath::Cos(3*(mpcsRP3-fvtx0sRP3)));
    cos_mpc_vtx_psi3->Fill(ibbcz,TMath::Cos(3*(mpcsRP3-vtx0sRP3)));
    cos_mpc_bbc_psi3->Fill(ibbcz,TMath::Cos(3*(mpcsRP3-bbcsRP3)));


    // float Qx2mpcs = TMath::Cos(2*mpcsRP2);
    // float Qy2mpcs = TMath::Sin(2*mpcsRP2)/*-Qy2mpcs_shift*/;
    // float Qw2mpcs = sumxy[1][0][2];

    // float Qx3mpcs = TMath::Cos(3*mpcsRP3)/*-Qx3mpcs_shift*/;
    // float Qy3mpcs = TMath::Sin(3*mpcsRP3);
    // float Qw3mpcs = sumxy[2][0][2];

    // float Qx2fvtx0s = TMath::Cos(2*fvtx0sRP2);
    // float Qy2fvtx0s = TMath::Sin(2*fvtx0sRP2)/*-Qy2fvtxs_shift*/;
    // float Qw2fvtx0s = sumxy[1][6][2];

    // float Qx3fvtx0s = TMath::Cos(3*fvtx0sRP3)/*-Qx3fvtxs_shift*/;
    // float Qy3fvtx0s = TMath::Sin(3*fvtx0sRP3);
    // float Qw3fvtx0s = sumxy[2][6][2];



    //if(Qw2mpcs<=0 || Qw3mpcs<=0) continue;
    //if(!(Qw2fvtx0s>5 && Qw2fvtx0s<1000)) continue;
    //if(!(Qw3fvtx0s>5 && Qw3fvtx0s<1000)) continue;


    //cout<<"filling segments"<<endl;
    bool first_track = true;
    if(!trimmed_tree)
      tot_zvtx->Fill(vtx_z);
    for(int itrk=0; itrk< nsegments; itrk++){
      if(itrk>=50) break;
      if(!trimmed_tree)
      {
        if(segmentok[itrk]!=1) continue;
        if(dca[itrk]*dca[itrk] > 0.125 || dca2d[itrk]*dca2d[itrk] > 0.0075) continue;

      //int nhits = nhit0[itrk]+nhit1[itrk]+nhit2[itrk]+nhit3[itrk];
      //cout<<"nhits: "<<nhits<<endl;
        if(nhit0[itrk]!=4) continue;
      }
      float phi = TMath::ATan2(py[itrk],px[itrk]);

      //float pi = 3.1419265;

      float pt = TMath::Sqrt(px[itrk]*px[itrk]+py[itrk]*py[itrk]);
      float eta = TMath::ASinH(pz[itrk]/pt);

      if(phi < -1.5)
        eta_east_top_dist->Fill(eta);
      else if(phi < 0.0)
        eta_west_bottom_dist->Fill(eta);
      else if(phi < 1.5)
        eta_west_top_dist->Fill(eta);
      else 
        eta_east_bottom_dist->Fill(eta);

      vtx_phi_dist->Fill(phi);

      //perform phi folding
      if(phi < -pi/2.)
      {
       phi+=2.*pi;
     }
     else if(phi > 3.*pi/2.)
     {
       phi-=2.*pi;
     }


     float theta = TMath::ACos(pz[itrk]/TMath::Sqrt(pt*pt+pz[itrk]*pz[itrk]));

     if(!trimmed_tree)
       hphi_zvtx->Fill(phi,vtx_z);
     hphi_theta->Fill(phi,theta);


     q_weight_mpc_pt->Fill(pt,sumxy[1][0][2]);
     q_weight_fvtx_pt->Fill(pt,sumxy[1][6][2]);
     q_weight_bbc_pt->Fill(pt,sumxy[1][4][2]);

     if(eta > 1.8 && first_track && !trimmed_tree)
     {
       forward_zvtx->Fill(vtx_z);
       first_track = false;
     }

     if(eta < -1.8 && first_track && !trimmed_tree)
     {
       backward_zvtx->Fill(vtx_z);
       first_track = false;
     }

     eta_tot_dist->Fill(eta);

     double mpc_dphi2 = phi - mpcsRP2;
     double bbc_dphi2 = phi - bbcsRP2;
     double fvtx_dphi2 = phi - fvtx0sRP2;

     double mpc_dphi3 = phi - mpcsRP3;
     double bbc_dphi3 = phi - bbcsRP3;
     double fvtx_dphi3 = phi - fvtx0sRP3;


      /*    float cosmpc_dphi2 = TMath::Cos(2*phi)*Qx2mpcs+TMath::Sin(2*phi)*Qy2mpcs;
      float cosmpc_dphi3 = TMath::Cos(3*phi)*Qx3mpcs+TMath::Sin(3*phi)*Qy3mpcs;

      float cosfvtx_dphi2 = TMath::Cos(2*phi)*Qx2fvtx0s+TMath::Sin(2*phi)*Qy2fvtx0s;
      float cosfvtx_dphi3 = TMath::Cos(3*phi)*Qx3fvtx0s+TMath::Sin(3*phi)*Qy3fvtx0s;
      */

      double cosmpc_dphi2 = TMath::Cos(2*mpc_dphi2);
      double cosmpc_dphi3 = TMath::Cos(3*mpc_dphi3);
      double cosbbc_dphi2 = TMath::Cos(2*bbc_dphi2);
      double cosbbc_dphi3 = TMath::Cos(3*bbc_dphi3);
      double cosfvtx_dphi2 = TMath::Cos(2*fvtx_dphi2);
      double cosfvtx_dphi3 = TMath::Cos(3*fvtx_dphi3);

      double sinmpc_dphi2 = TMath::Sin(2*mpc_dphi2);
      double sinmpc_dphi3 = TMath::Sin(3*mpc_dphi3);
      double sinbbc_dphi2 = TMath::Sin(2*bbc_dphi2);
      double sinbbc_dphi3 = TMath::Sin(3*bbc_dphi3);
      double sinfvtx_dphi2 = TMath::Sin(2*fvtx_dphi2);
      double sinfvtx_dphi3 = TMath::Sin(3*fvtx_dphi3);

      for(int iangle1 = 0; iangle1 < 7; iangle1++)
      {
        float blue_angle = 0.0016+0.0008-iangle1*2*0.0024/18;
        float blue_px = 100*TMath::Sin(blue_angle);
        float blue_py = 0.0;
        float blue_pz = 100*TMath::Cos(blue_angle);
        float proton_mass = 0.938;
        float blue_energy = TMath::Sqrt(100*100+proton_mass*proton_mass);
        TLorentzVector blue_beam(blue_px,blue_py,blue_pz,blue_energy);

        for(int iangle2 = 0; iangle2 < 9; iangle2++)
        {
          float yellow_angle = TMath::Pi()+0.0036+0.0016-iangle2*2*0.0048/24;
          float yellow_px = 100*TMath::Sin(yellow_angle);
          float yellow_py = 0.0;
          float yellow_pz = 100*TMath::Cos(yellow_angle);
          float au_mass = 197*0.939;
          float yellow_energy = TMath::Sqrt(100*100+au_mass*au_mass);
          TLorentzVector yellow_beam(yellow_px,yellow_py,yellow_pz,yellow_energy);
          int iangle = iangle1*9+iangle2;
          //if(iangle >72) continue;

          float mass = 0.1396;//assume charged pion mass
          float energy = TMath::Sqrt(px[itrk]*px[itrk]+py[itrk]*py[itrk]+pz[itrk]*pz[itrk]+mass*mass);
          TLorentzVector particle_vec(px[itrk],py[itrk],pz[itrk],energy);

          boost_and_rotate(particle_vec,blue_beam,yellow_beam);

          float phi_angle = TMath::ATan2(particle_vec.Py(),particle_vec.Px());
          float pt_angle = TMath::Sqrt(particle_vec.Py()*particle_vec.Py()+particle_vec.Px()*particle_vec.Px());
          float eta_angle = TMath::ASinH(particle_vec.Pz()/pt_angle);
          if(TMath::Abs(eta_angle) <= 0.35)
          {
            double bbc_dphi2_angle = bbcsRP2_angle[iangle] - phi_angle;
            double cosbbc_dphi2_angle = TMath::Cos(2*bbc_dphi2_angle);
            if(phi_angle < pi/3.0)
            {
             west_arm_bbcv2_angle[iangle]->Fill(pt_angle,cosbbc_dphi2_angle);
           }
           else if( phi_angle > pi/3.0)
           {
            east_arm_bbcv2_angle[iangle]->Fill(pt_angle,cosbbc_dphi2_angle);
          }
          both_arm_bbcv2_angle[iangle]->Fill(pt_angle,cosbbc_dphi2_angle);
        }
      }
    }

    if(TMath::Abs(eta) <= 0.35)
    {
     hadronmpcv2->Fill(pt,cosmpc_dphi2);
     hadronmpcv3->Fill(pt,cosmpc_dphi3);

     hadronbbcv2->Fill(pt,cosbbc_dphi2);
     hadronbbcv3->Fill(pt,cosbbc_dphi3);

     hadronmpcv2_sqr->Fill(pt,cosmpc_dphi2*cosmpc_dphi2);
     hadronmpcv3_sqr->Fill(pt,cosmpc_dphi3*cosmpc_dphi3);

     hadronbbcv2_sqr->Fill(pt,cosbbc_dphi2*cosbbc_dphi2);
     hadronbbcv3_sqr->Fill(pt,cosbbc_dphi3*cosbbc_dphi3);


     hadronfvtxv2->Fill(pt,cosfvtx_dphi2);
     hadronfvtxv3->Fill(pt,cosfvtx_dphi3);

     hadronfvtxv2_sqr->Fill(pt,cosfvtx_dphi2*cosfvtx_dphi2);
     hadronfvtxv3_sqr->Fill(pt,cosfvtx_dphi3*cosfvtx_dphi3);

     shadronmpcv2->Fill(pt,sinmpc_dphi2);
     shadronmpcv3->Fill(pt,sinmpc_dphi3);

     shadronbbcv2->Fill(pt,sinbbc_dphi2);
     shadronbbcv3->Fill(pt,sinbbc_dphi3);


     shadronfvtxv2->Fill(pt,sinfvtx_dphi2);
     shadronfvtxv3->Fill(pt,sinfvtx_dphi3);

     if(phi < pi/3.0)
     {
       eta_west_dist->Fill(eta);
       west_arm_hadronmpcv2->Fill(pt,cosmpc_dphi2);
       west_arm_hadronmpcv3->Fill(pt,cosmpc_dphi3);

       west_arm_hadronbbcv2->Fill(pt,cosbbc_dphi2);
       west_arm_hadronbbcv3->Fill(pt,cosbbc_dphi3);

       west_arm_hadronfvtxv2->Fill(pt,cosfvtx_dphi2);
       west_arm_hadronfvtxv3->Fill(pt,cosfvtx_dphi3);

     }
     else if( phi > pi/3.0)
     {

       eta_east_dist->Fill(eta);
       east_arm_hadronmpcv2->Fill(pt,cosmpc_dphi2);
       east_arm_hadronmpcv3->Fill(pt,cosmpc_dphi3);

       east_arm_hadronbbcv2->Fill(pt,cosbbc_dphi2);
       east_arm_hadronbbcv3->Fill(pt,cosbbc_dphi3);

       east_arm_hadronfvtxv2->Fill(pt,cosfvtx_dphi2);
       east_arm_hadronfvtxv3->Fill(pt,cosfvtx_dphi3);
     }

   }

   phi_v2_fvtx->Fill(phi,cosfvtx_dphi2);
   phi_v3_fvtx->Fill(phi,cosfvtx_dphi3);
   phi_v2_mpc->Fill(phi,cosmpc_dphi2);
   phi_v3_mpc->Fill(phi,cosmpc_dphi3);

   phi_v2_bbc->Fill(phi,cosbbc_dphi2);
   phi_v3_bbc->Fill(phi,cosbbc_dphi3);

   int iptbin = 6*pt/3.0;
   if(iptbin >=0 && iptbin < 6)
   {
     ptbin_phi_v2_fvtx[iptbin]->Fill(phi,cosfvtx_dphi2);
     ptbin_phi_v3_fvtx[iptbin]->Fill(phi,cosfvtx_dphi3);
     ptbin_phi_v2_mpc[iptbin]->Fill(phi,cosmpc_dphi2);
     ptbin_phi_v3_mpc[iptbin]->Fill(phi,cosmpc_dphi3);

     ptbin_phi_v2_bbc[iptbin]->Fill(phi,cosbbc_dphi2);
     ptbin_phi_v3_bbc[iptbin]->Fill(phi,cosbbc_dphi3);

   }
     /*
     if(phi < pi/3.0)
     {
       eta_west_dist->Fill(eta);
       west_arm_hadronmpcv2->Fill(pt,cosmpc_dphi2);
       west_arm_hadronmpcv3->Fill(pt,cosmpc_dphi3);

       west_arm_hadronbbcv2->Fill(pt,cosbbc_dphi2);
       west_arm_hadronbbcv3->Fill(pt,cosbbc_dphi3);

       west_arm_hadronfvtxv2->Fill(pt,cosfvtx_dphi2);
       west_arm_hadronfvtxv3->Fill(pt,cosfvtx_dphi3);

     }
     else if( phi > pi/3.0)
     {
       eta_east_dist->Fill(eta);
       east_arm_hadronmpcv2->Fill(pt,cosmpc_dphi2);
       east_arm_hadronmpcv3->Fill(pt,cosmpc_dphi3);

       east_arm_hadronbbcv2->Fill(pt,cosbbc_dphi2);
       east_arm_hadronbbcv3->Fill(pt,cosbbc_dphi3);

       east_arm_hadronfvtxv2->Fill(pt,cosfvtx_dphi2);
       east_arm_hadronfvtxv3->Fill(pt,cosfvtx_dphi3);
     }
*/
      //if(eta > -2.0 && eta < 2.0)
      //{
      //if(pt > 0.49 && pt < 3.01)
      //{
     eta_hadronfvtxv2->Fill(eta,cosfvtx_dphi2);
     eta_hadronfvtxv3->Fill(eta,cosfvtx_dphi3);

     eta_hadronmpcv2->Fill(eta,cosmpc_dphi2);
     eta_hadronmpcv3->Fill(eta,cosmpc_dphi3);

     eta_hadronbbcv2->Fill(eta,cosbbc_dphi2);
     eta_hadronbbcv3->Fill(eta,cosbbc_dphi3);

     eta_shadronfvtxv2->Fill(eta,sinfvtx_dphi2);
     eta_shadronfvtxv3->Fill(eta,sinfvtx_dphi3);

     eta_shadronmpcv2->Fill(eta,sinmpc_dphi2);
     eta_shadronmpcv3->Fill(eta,sinmpc_dphi3);

     eta_shadronbbcv2->Fill(eta,sinbbc_dphi2);
     eta_shadronbbcv3->Fill(eta,sinbbc_dphi3);

     eta_hadronfvtxv2_sqr->Fill(eta,cosfvtx_dphi2*cosfvtx_dphi2);
     eta_hadronmpcv2_sqr->Fill(eta,cosmpc_dphi2*cosmpc_dphi2);
     eta_hadronbbcv2_sqr->Fill(eta,cosbbc_dphi2*cosbbc_dphi2);


      //}

     int eta_bin = NETA*(eta+2.0)/4.0;
     if(eta_bin >= NETA || eta_bin < 0)
       continue;


     vtx_pt_eta[eta_bin]->Fill(pt);

     hadronfvtxv2_eta[eta_bin]->Fill(pt,cosfvtx_dphi2);
     hadronmpcv2_eta[eta_bin]->Fill(pt,cosmpc_dphi2);
     hadronbbcv2_eta[eta_bin]->Fill(pt,cosbbc_dphi2);

     shadronfvtxv2_eta[eta_bin]->Fill(pt,sinfvtx_dphi2);
     shadronmpcv2_eta[eta_bin]->Fill(pt,sinmpc_dphi2);
     shadronbbcv2_eta[eta_bin]->Fill(pt,sinbbc_dphi2);

     hadronfvtxv3_eta[eta_bin]->Fill(pt,cosfvtx_dphi3);
     hadronmpcv3_eta[eta_bin]->Fill(pt,cosmpc_dphi3);
     hadronbbcv3_eta[eta_bin]->Fill(pt,cosbbc_dphi3);

     shadronfvtxv3_eta[eta_bin]->Fill(pt,sinfvtx_dphi3);
     shadronmpcv3_eta[eta_bin]->Fill(pt,sinmpc_dphi3);
     shadronbbcv3_eta[eta_bin]->Fill(pt,sinbbc_dphi3);


     hadronfvtxv2_eta_sqr[eta_bin]->Fill(pt,cosfvtx_dphi2*cosfvtx_dphi2);
     hadronmpcv2_eta_sqr[eta_bin]->Fill(pt,cosmpc_dphi2*cosmpc_dphi2);
     hadronbbcv2_eta_sqr[eta_bin]->Fill(pt,cosbbc_dphi2*cosbbc_dphi2);

      //}

   }

    //cout<<"ntracklets: "<<ntracklets<<endl;

/*
    //cout<<"starting fvtx track loop"<<endl;
   for(int trk = 0; trk < ntracklets; trk++)
   {
     if(trk >=75) continue;
  //if(trk%10!=0) continue;
  //if(fphi[trk] < -1000 || feta[trk] < -1000 || fchisq[trk] < -1000) continue;
  //if(fchisq[trk] > 5 ) continue;
  //if(fnhits[trk]<3) continue;
     eta_fvtx_dist->Fill(feta[trk]);
     eta_tot_dist->Fill(feta[trk]);
     double mpc_dphi2 = fphi[trk] - mpcsRP2;
     double bbc_dphi2 = fphi[trk] - bbcsRP2;
     double fvtx_dphi2 = fphi[trk] - fvtx0sRP2;

     double mpc_dphi3 = fphi[trk] - mpcsRP3;
     double bbc_dphi3 = fphi[trk] - bbcsRP3;
     double fvtx_dphi3 = fphi[trk] - fvtx0sRP3;

     double cosmpc_dphi2 = TMath::Cos(2*mpc_dphi2);
     double cosmpc_dphi3 = TMath::Cos(3*mpc_dphi3);
     double cosbbc_dphi2 = TMath::Cos(2*bbc_dphi2);
     double cosbbc_dphi3 = TMath::Cos(3*bbc_dphi3);
     double cosfvtx_dphi2 = TMath::Cos(2*fvtx_dphi2);
     double cosfvtx_dphi3 = TMath::Cos(3*fvtx_dphi3);

     double sinmpc_dphi2 = TMath::Sin(2*mpc_dphi2);
     double sinmpc_dphi3 = TMath::Sin(3*mpc_dphi3);
     double sinbbc_dphi2 = TMath::Sin(2*bbc_dphi2);
     double sinbbc_dphi3 = TMath::Sin(3*bbc_dphi3);
     double sinfvtx_dphi2 = TMath::Sin(2*fvtx_dphi2);
     double sinfvtx_dphi3 = TMath::Sin(3*fvtx_dphi3);

     feta_hadronmpcv2->Fill(feta[trk],cosmpc_dphi2);
     feta_hadronmpcv3->Fill(feta[trk],cosmpc_dphi3);

     feta_hadronbbcv2->Fill(feta[trk],cosbbc_dphi2);
     feta_hadronbbcv3->Fill(feta[trk],cosbbc_dphi3);

     feta_hadronfvtxv2->Fill(feta[trk],cosfvtx_dphi2);
     feta_hadronfvtxv3->Fill(feta[trk],cosfvtx_dphi3);

     feta_shadronmpcv2->Fill(feta[trk],sinmpc_dphi2);
     feta_shadronmpcv3->Fill(feta[trk],sinmpc_dphi3);

     feta_shadronbbcv2->Fill(feta[trk],sinbbc_dphi2);
     feta_shadronbbcv3->Fill(feta[trk],sinbbc_dphi3);

     feta_shadronfvtxv2->Fill(feta[trk],sinfvtx_dphi2);
     feta_shadronfvtxv3->Fill(feta[trk],sinfvtx_dphi3);

     feta_hadronmpcv2_sqr->Fill(feta[trk],cosmpc_dphi2*cosmpc_dphi2);
     feta_hadronbbcv2_sqr->Fill(feta[trk],cosbbc_dphi2*cosbbc_dphi2);
     feta_hadronfvtxv2_sqr->Fill(feta[trk],cosfvtx_dphi2*cosfvtx_dphi2);


   }
*/
 }


 hadronmpcv2->SetTitle("MPC VTX v2");
 hadronmpcv3->SetTitle("MPC VTX v3");
 hadronfvtxv2->SetTitle("FVTX VTX v2");
 hadronfvtxv3->SetTitle("FVTX VTX v3");

 hadronmpcv2->SetXTitle("pT");
 hadronmpcv3->SetXTitle("pT");
 hadronfvtxv2->SetXTitle("pT");
 hadronfvtxv3->SetXTitle("pT");

 cout<<"finished loop, writing out histograms to file"<<endl;
  //eta_east_dist->Scale(1.0/nentries);
  //eta_west_dist->Scale(1.0/nentries);

 feta_hadronmpcv2->Write();
 feta_hadronmpcv3->Write();

 feta_hadronbbcv2->Write();
 feta_hadronbbcv3->Write();

 feta_hadronfvtxv2->Write();
 feta_hadronfvtxv3->Write();

 feta_shadronmpcv2->Write();
 feta_shadronmpcv3->Write();

 feta_shadronbbcv2->Write();
 feta_shadronbbcv3->Write();

 feta_shadronfvtxv2->Write();
 feta_shadronfvtxv3->Write();

 feta_hadronmpcv2_sqr->Write();
 feta_hadronbbcv2_sqr->Write();
 feta_hadronfvtxv2_sqr->Write();


 phi_v2_fvtx->Write();
 phi_v3_fvtx->Write();
 phi_v2_mpc->Write();
 phi_v3_mpc->Write();
 phi_v2_bbc->Write();
 phi_v3_bbc->Write();

 for(int izbin = 0; izbin < NZPS; izbin++)
 {
  q_vec_fvtx_v2_bf[izbin]->Write();
  q_vec_vtx_v2_bf[izbin]->Write();
  q_vec_mpc_v2_bf[izbin]->Write();
  q_vec_fvtx_v2_af[izbin]->Write();
  q_vec_vtx_v2_af[izbin]->Write();
  q_vec_mpc_v2_af[izbin]->Write();
  q_vec_bbc_v2_bf[izbin]->Write();
  q_vec_bbc_v2_af[izbin]->Write();

  q_vec_fvtx_v3_bf[izbin]->Write();
  q_vec_vtx_v3_bf[izbin]->Write();
  q_vec_mpc_v3_bf[izbin]->Write();
  q_vec_fvtx_v3_af[izbin]->Write();
  q_vec_vtx_v3_af[izbin]->Write();
  q_vec_mpc_v3_af[izbin]->Write();
  q_vec_bbc_v3_bf[izbin]->Write();
  q_vec_bbc_v3_af[izbin]->Write();
}

cos_bbc_vtx_psi2->Write();
cos_bbc_fvtx_psi2->Write();
cos_vtx_fvtx_psi2->Write();
cos_mpc_fvtx_psi2->Write();
cos_mpc_vtx_psi2->Write();
cos_mpc_bbc_psi2->Write();

all_vtx_pt->Write();


bbc_vtx_pt->Write();
fvtx_vtx_pt->Write();
mpc_vtx_pt->Write();


cos_bbc_vtx_psi3->Write();
cos_bbc_fvtx_psi3->Write();
cos_vtx_fvtx_psi3->Write();
cos_mpc_fvtx_psi3->Write();
cos_mpc_vtx_psi3->Write();
cos_mpc_bbc_psi3->Write();

q_vec_2_mpc_mb->Write();
q_vec_2_mpc_ct->Write();
q_vec_3_mpc_mb->Write();
q_vec_3_mpc_ct->Write();

q_vec_2_fvtx_mb->Write();
q_vec_2_fvtx_ct->Write();
q_vec_3_fvtx_mb->Write();
q_vec_3_fvtx_ct->Write();

backward_zvtx->Write();
forward_zvtx->Write();
tot_zvtx->Write();

hadronmpcv2->Write();
hadronmpcv3->Write();
hadronbbcv2->Write();
hadronbbcv3->Write();
hadronfvtxv2->Write();
hadronfvtxv3->Write();

hadronmpcv2_sqr->Write();
hadronmpcv3_sqr->Write();
hadronbbcv2_sqr->Write();
hadronbbcv3_sqr->Write();
hadronfvtxv2_sqr->Write();
hadronfvtxv3_sqr->Write();


hphi_zvtx->Write();
hphi_theta->Write();

eta_east_dist->Write();
eta_west_dist->Write();

for(int i = 0; i < 16; i++)
{
  hadronfvtxv2_eta[i]->Write();
  hadronmpcv2_eta[i]->Write();
  hadronbbcv2_eta[i]->Write();

  shadronfvtxv2_eta[i]->Write();
  shadronmpcv2_eta[i]->Write();
  shadronbbcv2_eta[i]->Write();

  hadronfvtxv3_eta[i]->Write();
  hadronmpcv3_eta[i]->Write();
  hadronbbcv3_eta[i]->Write();

  shadronfvtxv3_eta[i]->Write();
  shadronmpcv3_eta[i]->Write();
  shadronbbcv3_eta[i]->Write();

  vtx_pt_eta[i]->Write();

  hadronfvtxv2_eta_sqr[i]->Write();
  hadronmpcv2_eta_sqr[i]->Write();
  hadronbbcv2_eta_sqr[i]->Write();
}

for(int i = 0; i < 6; i++)
{
  ptbin_phi_v2_fvtx[i]->Write();
  ptbin_phi_v3_fvtx[i]->Write();
  ptbin_phi_v2_mpc[i]->Write();
  ptbin_phi_v3_mpc[i]->Write();
  ptbin_phi_v2_bbc[i]->Write();
  ptbin_phi_v3_bbc[i]->Write();

}
eta_hadronmpcv2->Write();
eta_hadronmpcv3->Write();
eta_hadronbbcv2->Write();
eta_hadronbbcv3->Write();
eta_hadronfvtxv2->Write();
eta_hadronfvtxv3->Write();

eta_hadronmpcv2_sqr->Write();
  //eta_hadronmpcv3_sqr->Write();
eta_hadronbbcv2_sqr->Write();
  //eta_hadronbbcv3_sqr->Write();
eta_hadronfvtxv2_sqr->Write();
  //eta_hadronfvtxv3_sqr->Write();


eta_fvtx_dist->Write();
eta_tot_dist->Write();


  //sin histograms
shadronmpcv2->Write();
shadronmpcv3->Write();
shadronbbcv2->Write();
shadronbbcv3->Write();
shadronfvtxv2->Write();
shadronfvtxv3->Write();

eta_shadronmpcv2->Write();
eta_shadronmpcv3->Write();
eta_shadronbbcv2->Write();
eta_shadronbbcv3->Write();
eta_shadronfvtxv2->Write();
eta_shadronfvtxv3->Write();

  //east west arm
west_arm_hadronmpcv2->Write();
west_arm_hadronmpcv3->Write();
west_arm_hadronbbcv2->Write();
west_arm_hadronbbcv3->Write();
west_arm_hadronfvtxv2->Write();
west_arm_hadronfvtxv3->Write();

east_arm_hadronmpcv2->Write();
east_arm_hadronmpcv3->Write();
east_arm_hadronbbcv2->Write();
east_arm_hadronbbcv3->Write();
east_arm_hadronfvtxv2->Write();
east_arm_hadronfvtxv3->Write();

eta_east_bottom_dist->Write();
eta_east_top_dist->Write();
eta_west_bottom_dist->Write();
eta_west_top_dist->Write();

for(int iangle = 0; iangle < n_bbc_angles; iangle++)
{
  east_arm_bbcv2_angle[iangle]->Write();
  west_arm_bbcv2_angle[iangle]->Write();
  both_arm_bbcv2_angle[iangle]->Write();
  for(int izbin = 0; izbin < NZPS; izbin++)
  {
    q_vec_bbc_angle_bf[iangle][izbin]->Write();
    q_vec_bbc_angle_af[iangle][izbin]->Write();
  }
}

vtx_phi_dist->Write();

q_weight_mpc_pt->Write();
q_weight_fvtx_pt->Write();
q_weight_bbc_pt->Write();

if(calFlag < 3 && calFlag>0)
{
  for (int iz=0; iz<NZPS; iz++) {
    for (int ih=0; ih<NHAR; ih++) {
     for (int id=0; id<NDET; id++) {
       if(iz==0)
       {
         dis[ih][id]->Write();
         qx[ih][id]->Write();
         qy[ih][id]->Write();
       }

       flt[iz][ih][id]->Write();
       ave[iz][ih][id]->Write();
     }
   }
 }
 for(int iangle = 0; iangle < n_bbc_angles; iangle++)
 {
  for (int iz=0; iz<NZPS; iz++) {
   if(iz==0)
   {
     bbc_angle_dis[iangle]->Write();
     bbc_angle_qx[iangle]->Write();
     bbc_angle_qy[iangle]->Write();
   }

   bbc_angle_flt[iz][iangle]->Write();
   bbc_angle_ave[iz][iangle]->Write();
 }
}
}

if(calFlag < 3 && calFlag>0)
{
  sprintf(name,"vtx_ep_calib/new_wide_z_flattening/flattening_%d_%d.dat",runnumber,calFlag);
      //sprintf(name,"tsukuba/flattening_%d_%d.dat",runnumber,calFlag);
  cout << "writing calibration file : " << name << endl;
  ofstream ofs;
  ofs.open(name);

  //cout << "flt para" << endl;
  for (int iz=0; iz<NZPS; iz++) {
   for (int ih=0; ih<NHAR; ih++) {
     for (int id=0; id<NDET; id++) {
       for (int ib=0; ib<2; ib++) {
        //cout<<"ave["<<iz<<"]["<<ih<<"]["<<id<<"]->Get("<<ib+1<<")"<<endl;
         ofs << ave[iz][ih][id]->GetBinContent(ib+1) << " ";
         ofs << ave[iz][ih][id]->GetBinError  (ib+1) << " ";
       }
       ofs << endl;
       for (int ib=0; ib<2; ib++) {
         for (int io=0; io<NORD; io++) {
    //cout<<"flt["<<iz<<"]["<<ih<<"]["<<id<<"]->Get("<<ib*NORD+io+1<<")"<<endl;
          ofs << flt[iz][ih][id]->GetBinContent(ib*NORD+io+1) << " ";
        }
        ofs << endl;
      }
    }
  }
}

ofs.close();

sprintf(name,"vtx_ep_calib/new_wide_z_flattening/flattening_bbc_angle_%d_%d.dat",runnumber,calFlag);
      //sprintf(name,"tsukuba/flattening_%d_%d.dat",runnumber,calFlag);
cout << "writing calibration 2nd file : " << name << endl;
ofs.open(name);

for(int iz = 0; iz < NZPS; iz++)
{
  for(int iangle = 0; iangle < n_bbc_angles; iangle++)
  {
    for(int ib = 0; ib < 2; ib++)
    {
      ofs << bbc_angle_ave[iz][iangle]->GetBinContent(ib+1) << " ";
      ofs << bbc_angle_ave[iz][iangle]->GetBinError(ib+1) << " ";
    }
    ofs << endl;
    for(int ib=0; ib<2; ib++)
    {
      for(int io = 0; io < NORD; io++)
      {
       ofs << bbc_angle_flt[iz][iangle]->GetBinContent(ib*NORD+io+1) << " ";
     }
     ofs << endl;
   }
 }
}

}

outfile->Close();


}
