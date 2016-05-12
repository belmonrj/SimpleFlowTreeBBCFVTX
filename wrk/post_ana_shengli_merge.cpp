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

//#include "Run15pAupc3dphidzcalibsmoothpass1.h"


using namespace std;

void boost_and_rotate(TLorentzVector & vec, TLorentzVector input1, TLorentzVector input2);

void post_ana_shengli_merge(int runNumber = 435823, int rp_recal_pass = 1){

  char outFile1[300];
  sprintf(outFile1,"%s%d%s","/phenix/plhf/theok/taxi/Run15pAu200CAMBPro104/8691/processed/shengli/hist_",runNumber,".root");

  char outFile2[100];
  sprintf(outFile2,"%s%d%s","vtx_ep_calib/shengli/rp/hrp_",runNumber,".root");

  cout<<"runNumber = " <<runNumber<<" "
  <<"rp_recal_pass = "<<rp_recal_pass<<endl;

  static const int n_bbc_angles = 63;

   //histogram
  //char cc[100];

  float pi = acos(-1.0);

  TH2D     *qx[NMUL][NHAR][NDET];
  TH2D     *qy[NMUL][NHAR][NDET];
  TProfile *ave[NMUL][NZPS][NHAR][NDET];
  TProfile *flt[NMUL][NZPS][NHAR][NDET];
  TH2D     *dis[NMUL][NHAR][NDET];
  float    mean[NMUL][NZPS][NHAR][NDET][2];
  float    widt[NMUL][NZPS][NHAR][NDET][2];
  float    four[NMUL][NZPS][NHAR][NDET][2][NORD];

  float    bbc_angle_mean[NZPS][n_bbc_angles][2];
  float    bbc_angle_widt[NZPS][n_bbc_angles][2];
  float    bbc_angle_four[NZPS][n_bbc_angles][2][NORD];

  TH2D     *bbc_angle_qx[n_bbc_angles];
  TH2D     *bbc_angle_qy[n_bbc_angles];
  TProfile *bbc_angle_ave[NZPS][n_bbc_angles];
  TProfile *bbc_angle_flt[NZPS][n_bbc_angles];
  TH2D     *bbc_angle_dis[n_bbc_angles];


  char name[200];
  for (int ic=0; ic<NMUL; ic++) {
    for (int iz=0; iz<NZPS; iz++) {
      for (int ih=0; ih<NHAR; ih++) {
       for (int id=0; id<NDET; id++) {
         sprintf(name,"ave_%d_%d_%d_%d",ic,iz,ih,id);
	       ave[ic][iz][ih][id] = new TProfile(name,name,4,-0.5,3.5,-10.1,10.1,"S");//for SMD -1.1,1.1

	       sprintf(name,"flt_%d_%d_%d_%d",ic,iz,ih,id);
	       flt[ic][iz][ih][id] = new TProfile(name,name,4*NORD,-0.5,NORD*4.0-0.5,-1.1,1.1);

      }
    }
  }
}

for (int ic=0; ic<NMUL; ic++) {
  for (int ih=0; ih<NHAR; ih++) {
    for (int id=0; id<NDET; id++) {
     sprintf(name,"dis_%d_%d_%d",ic,ih,id);
     dis[ic][ih][id] = new TH2D(name,name,NZPS*3,-0.5,NZPS*3.0-0.5,50,-pi,pi);

     sprintf(name,"qx_%d_%d_%d",ic,ih,id);
     qx[ic][ih][id] = new TH2D(name,name,NZPS*3,-0.5,NZPS*3.0-0.5, 220,-4.1,4.1);

     sprintf(name,"qy_%d_%d_%d",ic,ih,id);
     qy[ic][ih][id] = new TH2D(name,name,NZPS*3,-0.5,NZPS*3.0-0.5, 220,-4.1,4.1);
   }
 }
}

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

for(int iangle = 0; iangle < n_bbc_angles; iangle++)
{
  sprintf(name,"bbc_angle_dis_%d",iangle);
  bbc_angle_dis[iangle] = new TH2D(name,name,NZPS*3,-0.5,NZPS*3.0-0.5,50,-pi,pi);

  sprintf(name,"bbc_angle_qx_%d",iangle);
  bbc_angle_qx[iangle] = new TH2D(name,name,NZPS*3,-0.5,NZPS*3.0-0.5, 220,-4.1,4.1);

  sprintf(name,"bbc_angle_qy_%d",iangle);
  bbc_angle_qy[iangle] = new TH2D(name,name,NZPS*3,-0.5,NZPS*3.0-0.5, 220,-4.1,4.1);

}

int calFlag = rp_recal_pass;
cout << "calibration flag : " << calFlag << endl;

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

if(calFlag>=2){
  sprintf(name,"vtx_ep_calib/new_wide_z_flattening/flattening_%d_%d.dat",runNumber,calFlag-1);

  cout << "reading calibration file : " << name << endl;
    float f0,f1,f2,f3;//f4,f5,f6,f7;
    ifstream ifs;
    ifs.open(name);
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

  ifs.open(Form("vtx_ep_calib/new_wide_z_flattening/flattening_bbc_angle_%d_%d.dat",runNumber,calFlag-1));
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


  //centrality

  //psi2 dis
TH1F *disnew[10][2];
TH1F *disnew2[10][2];
for (int ic=0; ic<10; ic++) {
  for (int id=0; id<2; id++) {
    sprintf(name,"disnew_%d_%d",ic,id);
    disnew[ic][id] = new TH1F(name, name, 100, -pi/2.0, pi/2.0);

    sprintf(name,"disnew2_%d_%d",ic,id);
    disnew2[ic][id] = new TH1F(name, name, 100, -pi/2.0, pi/2.0);
  }
}


  //add by shengli for east and west resolution
  /*
  TH2F *uhadronv1bbcs = new TH2F("uhadronv1bbcs","uhadronv1bbcs", 15, 0.0, 3.0, 2200, -1.1, 1.1);
  TH2F *uhadwestv1bbcs = new TH2F("uhadwestv1bbcs","uhadwestv1bbcs", 15, 0.0, 3.0, 2200, -1.1, 1.1);
  TH2F *uhadeastv1bbcs = new TH2F("uhadeastv1bbcs","uhadeastv1bbcs", 15, 0.0, 3.0, 2200, -1.1, 1.1);

  TH2F *uhadronv2bbcs = new TH2F("uhadronv2bbcs","uhadronv2bbcs", 15, 0.0, 3.0, 2200, -1.1, 1.1);
  TH2F *uhadwestv2bbcs = new TH2F("uhadwestv2bbcs","uhadwestv2bbcs", 15, 0.0, 3.0, 2200, -1.1, 1.1);
  TH2F *uhadeastv2bbcs = new TH2F("uhadeastv2bbcs","uhadeastv2bbcs", 15, 0.0, 3.0, 2200, -1.1, 1.1);
  */


  TProfile *bbcs_v2_incl = new TProfile("bbcs_v2_incl","bbcs_v2_incl",15, 0.0, 3.0,-1.1,1.1);
  TProfile *bbcs_v2_east = new TProfile("bbcs_v2_east","bbcs_v2_east",15, 0.0, 3.0,-1.1,1.1);
  TProfile *bbcs_v2_west = new TProfile("bbcs_v2_west","bbcs_v2_west",15, 0.0, 3.0,-1.1,1.1);

  vector<TProfile * > bbcs_v2_west_angle;
  vector<TProfile * > bbcs_v2_east_angle;
  vector<TProfile * > bbcs_v2_incl_angle;

  for(int i = 0; i < n_bbc_angles; i++)
  {
    bbcs_v2_west_angle.push_back(new TProfile(Form("bbcs_v2_west_angle_%d",i),Form("bbcs_v2_west_angle_%d",i),15, 0.0, 3.0,-1.1,1.1));
    bbcs_v2_east_angle.push_back(new TProfile(Form("bbcs_v2_east_angle_%d",i),Form("bbcs_v2_east_angle_%d",i),15, 0.0, 3.0,-1.1,1.1));
    bbcs_v2_incl_angle.push_back(new TProfile(Form("bbcs_v2_incl_angle_%d",i),Form("bbcs_v2_incl_angle_%d",i),15, 0.0, 3.0,-1.1,1.1));
  }

  //tree invariables
  static const int max_nh = 10;
  
  // theo tree structure

  float    d_bbcz;    // bbcz

  float d_Qx[9];
  float d_Qy[9];
  float d_Qw[9];
  float d_BBCs_Q[189];

  int   d_nsegments;
  float d_px[max_nh];
  float d_py[max_nh];
  float d_pz[max_nh];

  int nevent=0;

  /*
  Char_t inFile[100];
  
  Char_t inFile1[100];

  sprintf(inFile,"%s%d%s","lst/run_",runNumber,".lst");
  cout<<"*********************************** "<<inFile<<endl;
  char filename[500];

  int ifile=0; 
  ifstream fin(inFile);
  ifstream fin1(inFile1);
  
  while(fin.getline(filename, 500)){
    ifile++;
  */
    char filename[500];
    sprintf(filename,"/gpfs/mnt/gpfs02/phenix/plhf/plhf1/theok/taxi/Run15pAu200CAMBPro104/8691/data/%d.root",runNumber);
  //TFile *f = TFile::Open(Form("/gpfs/mnt/gpfs02/phenix/plhf/plhf1/theok/taxi/Run15pAu200CAMBPro104/8691/data/%d.root",runNumber));
    cout << " ana file " << filename << endl;

    TFile *f=new TFile( filename);

    if(!f) return;

    
    TTree *htree = (TTree *)f->Get("ntp_event");

    if(!htree) return;


    TBranch *b_bbcz = htree->GetBranch("bbc_z");
    TBranch *b_Qx = htree->GetBranch("d_Qx");
    TBranch *b_Qy = htree->GetBranch("d_Qy");
    TBranch *b_Qw = htree->GetBranch("d_Qw");

    TBranch *b_d_BBCs_Q = htree->GetBranch("d_BBCs_Q");
    b_d_BBCs_Q->SetAddress(&d_BBCs_Q);


    b_bbcz->SetAddress(&d_bbcz);
    b_Qx->SetAddress(d_Qx);
    b_Qy->SetAddress(d_Qy);
    b_Qw->SetAddress(d_Qw);

    cout<<"finish Qw"<<endl;
    

    TBranch *b_nsegments = htree->GetBranch("nsegments");
    TBranch *b_px = htree->GetBranch("px");
    TBranch *b_py = htree->GetBranch("py");
    TBranch *b_pz = htree->GetBranch("pz");

    b_nsegments->SetAddress(&d_nsegments);
    b_px->SetAddress(d_px);
    b_py->SetAddress(d_py);
    b_pz->SetAddress(d_pz);

    int nentries = htree->GetEntries();
    //int nentries1 = htree1->GetEntries();
    cout<<"total events = " << nentries<<endl;
    for ( Int_t ievt = 0 ; ievt < nentries ; ievt++ ) {
    //for ( Int_t ievt = 0 ; ievt < 1e5 ; ievt++ ) {
      nevent++;
      
      //for ( Int_t ievt = 0 ; ievt < 200 ; ievt++ ) {
      if(ievt%10000==0) cout<<"events number = "<<ievt<<" "<<nevent<<endl;
      //if(nevent ==5) break;
      b_bbcz->GetEntry(ievt);

      b_Qx->GetEntry(ievt);
      b_Qy->GetEntry(ievt);
      b_Qw->GetEntry(ievt);
      b_d_BBCs_Q->GetEntry(ievt);
      //track
      b_nsegments->GetEntry(ievt);
      b_px->GetEntry(ievt);
      b_py->GetEntry(ievt);
      b_pz->GetEntry(ievt);
      
      
      int ibbcz  = NZPS*(d_bbcz+15)/30;

      if(ibbcz<0||ibbcz>=NZPS) continue;

      int icent=0;      
      if(icent<0||icent>=NMUL) continue;

      if(rp_recal_pass<1) continue;

      //cout<<"faint 1 "<<d_nh<<" "<<d_nvtxtrk<<endl;

      float sumxy[NHAR][NDET][4];
      for (int i=0; i<NHAR; i++) {
       for (int j=0; j<NDET; j++) {
         for (int k=0; k<4; k++) {
           sumxy[i][j][k]=0;
         }
       }
     }

     float         bbcsRP1[NDET]; 
     float         bbcsRP2[NDET]; 
     float         bbcsRP2_angle[n_bbc_angles];

     sumxy[0][0][0] = d_Qx[2];
     sumxy[0][0][1] = d_Qy[2];
     sumxy[0][0][2] = d_Qw[2];

     sumxy[1][0][0] = d_Qx[5];
     sumxy[1][0][1] = d_Qy[5];
     sumxy[1][0][2] = d_Qw[5];

     float sumxy_bbc_angle[n_bbc_angles][4];
     for(int i = 0; i < n_bbc_angles;i++)
      for(int j = 0; j<4; j++)
        sumxy_bbc_angle[i][j]=0;

      for(int iangle = 0; iangle < n_bbc_angles; iangle++)
      {
       for(int icomp = 0; icomp < 3; icomp++)
         sumxy_bbc_angle[iangle][icomp] = d_BBCs_Q[iangle*3 + icomp];
     }

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

    if (calFlag < 3 && calFlag>0) bbc_angle_dis[iangle]->Fill(ibbcz+NZPS*2.0,sumxy_bbc_angle[iangle][3]*(2));

  }
  else sumxy_bbc_angle[iangle][3]=-9999.9;
}


for (int ih=0; ih<NHAR; ih++) {
 for(int id=0;id<NDET;id++){
   if (sumxy[ih][id][2]>0.0) {
     sumxy[ih][id][3]=atan2(sumxy[ih][id][1],sumxy[ih][id][0])/(ih+1.0);
     if (calFlag>0) dis[icent][ih][id]->Fill(ibbcz,sumxy[ih][id][3]*(ih+1.0));
   }
   if (sumxy[ih][id][2]>0.0) {
     for (int ib=0; ib<2; ib++) {
       sumxy[ih][id][ib]/=sumxy[ih][id][2];
	      //if(ih==1 && id==0 && ib==0 && sumxy[ih][id][ib]>1) cout<<sumxy[ih][id][ib]<<endl;
       if (calFlag>0) {
        ave[icent][ibbcz][ih][id]->Fill(ib+0.0,sumxy[ih][id][ib]);
        if(ib==0) qx[icent][ih][id]->Fill(ibbcz,sumxy[ih][id][0]);
        if(ib==1) qy[icent][ih][id]->Fill(ibbcz,sumxy[ih][id][1]);
      }
      float sxy=sumxy[ih][id][ib];
      float mxy=mean[icent][ibbcz][ih][id][ib];
      float wxy=widt[icent][ibbcz][ih][id][ib];

	      //if(icent==0 && ibbcz==0 && ih==1 && id==0) cout<<ib<<" "<<sxy<<" "<<mxy<<" "<<wxy<<endl;
      sumxy[ih][id][ib]=(sxy-mxy)/wxy;
      if (calFlag>0) {
        ave[icent][ibbcz][ih][id]->Fill(ib+2.0,sumxy[ih][id][ib]);
        if(ib==0) qx[icent][ih][id]->Fill(ibbcz+NZPS,sumxy[ih][id][0]);
        if(ib==1) qy[icent][ih][id]->Fill(ibbcz+NZPS,sumxy[ih][id][1]);
      }
    }

    sumxy[ih][id][3]=atan2(sumxy[ih][id][1],sumxy[ih][id][0])/(ih+1.0);
    if (calFlag>0) {
     dis[icent][ih][id]->Fill(ibbcz+NZPS,sumxy[ih][id][3]*(ih+1.0));
   }

   float psi=sumxy[ih][id][3]*(ih+1.0);
   float dp=0.0;
   for (int io=0; io<NORD; io++) {
     float cc=cos((io+1.0)*psi);
     float ss=sin((io+1.0)*psi);
     if (calFlag>0) flt[icent][ibbcz][ih][id]->Fill(io+0.0,cc);
     if (calFlag>0) flt[icent][ibbcz][ih][id]->Fill(io+NORD,ss);
	      float aa=four[icent][ibbcz][ih][id][0][io]; // mean cos
	      float bb=four[icent][ibbcz][ih][id][1][io]; // mean sin
	      dp+=(aa*ss-bb*cc)*2.0/(io+1.0);
	    }
	    psi+=dp;
	    psi=atan2(sin(psi),cos(psi));
	    for (int io=0; io<NORD; io++) {
       float cc=cos((io+1.0)*psi);
       float ss=sin((io+1.0)*psi);
       if (calFlag>0) flt[icent][ibbcz][ih][id]->Fill(io+NORD*2.0,cc);
       if (calFlag>0) flt[icent][ibbcz][ih][id]->Fill(io+NORD*3.0,ss);
     }
     sumxy[ih][id][3]=psi/(ih+1.0);
     if (calFlag>0) dis[icent][ih][id]->Fill(ibbcz+NZPS*2.0,sumxy[ih][id][3]*(ih+1.0));
   } else {
     sumxy[ih][id][3]=-9999.9;
   }

	}//end of id
      }//end of ih

      
      //if(icent!=0) continue;
      
      
      for(int id=0; id<NDET; id++)
      {
	     bbcsRP1[id] = (sumxy[0][id][2]>0)?sumxy[0][id][3]:-9999.9;//south
	     bbcsRP2[id] = (sumxy[1][id][2]>0)?sumxy[1][id][3]:-9999.9;//south
      }

      for(int iangle = 0; iangle < n_bbc_angles; iangle++)
      {
        bbcsRP2_angle[iangle] = (sumxy_bbc_angle[iangle][2]>0)?sumxy_bbc_angle[iangle][3]:-9999.9;
        //cout<<"bbcsrp2_angle["<<iangle<<"]: "<<bbcsRP2_angle[iangle]<<endl;
      }

      //-----------------------------------------v2------------------------------------------------------------

      //continue;

      if(rp_recal_pass<3) continue;
      //cout<<d_nh<<" "<<d_nvtxtrk<<endl;


      int id = 0;//sqrt(NDET)*iblue+iyellow;
        //cout<<"looping over segments"<<endl;

      //start of vtx stand alone trk
      for(int itrk=0; itrk< d_nsegments; itrk++){
                //cout<<"segment: "<<itrk<<endl;
       float px    = d_px[itrk];
       float py    = d_py[itrk];
       float pz    = d_pz[itrk];
       int dcarm=0;
       if(px>0) dcarm=1;


       float phi0 = TMath::ATan2(py,px);
       float pt = sqrt(px*px+py*py);

	     //v2
       float dphi_2 = -9999;
       dphi_2 = phi0 - bbcsRP2[id];
	     //cout<<"bbcsRP2[id]: "<<bbcsRP2[id]<<endl;
       if(-4.0<bbcsRP2[id] && bbcsRP2[id]<4.0 ){
         bbcs_v2_incl->Fill(pt,cos(2*dphi_2));
           //cout<<"filling histos"<<endl;
         if(dcarm==0) 
         {
          bbcs_v2_east->Fill(pt,cos(2*dphi_2));
        }
        else if(dcarm==1)
        {
          bbcs_v2_west->Fill(pt,cos(2*dphi_2));
        }
      }

      for(int iangle1 = 0; iangle1 < 7; iangle1++)
      {
        float blue_angle = 0.0016+0.0024-iangle1*2*0.0024/6;
        float blue_px = 100*TMath::Sin(blue_angle);
        float blue_py = 0.0;
        float blue_pz = 100*TMath::Cos(blue_angle);
        float proton_mass = 0.938;
        float blue_energy = TMath::Sqrt(100*100+proton_mass*proton_mass);
        TLorentzVector blue_beam(blue_px,blue_py,blue_pz,blue_energy);

        for(int iangle2 = 0; iangle2 < 9; iangle2++)
        {
          float yellow_angle = TMath::Pi()+0.0036+0.0048-iangle2*2*0.0048/8;
          float yellow_px = 100*TMath::Sin(yellow_angle);
          float yellow_py = 0.0;
          float yellow_pz = 100*TMath::Cos(yellow_angle);
          float yellow_energy = TMath::Sqrt(100*100+proton_mass*proton_mass);
          TLorentzVector yellow_beam(yellow_px,yellow_py,yellow_pz,yellow_energy);
          int iangle = iangle1*9+iangle2;
          //if(iangle >72) continue;

          if(!(-4.0<bbcsRP2_angle[iangle] && bbcsRP2_angle[iangle]<4.0 )) continue;

          float mass = 0.1396;//assume charged pion mass
          float energy = TMath::Sqrt(px*px+py*py+pz*pz+mass*mass);
          TLorentzVector particle_vec(px,py,pz,energy);

          boost_and_rotate(particle_vec,blue_beam,yellow_beam);

          float phi_angle = TMath::ATan2(particle_vec.Py(),particle_vec.Px());
          float pt_angle = TMath::Sqrt(particle_vec.Py()*particle_vec.Py()+particle_vec.Px()*particle_vec.Px());
          //float eta_angle = TMath::ASinH(particle_vec.Pz()/pt_angle);
          //if(TMath::Abs(eta_angle) <= 1.0)
          //{
          double bbc_dphi2_angle = bbcsRP2_angle[iangle] - phi_angle;
          double cosbbc_dphi2_angle = TMath::Cos(2*bbc_dphi2_angle);
          //cout<<"filling bbc angle v2 for iangle: "<<iangle<<" is:  "<<cosbbc_dphi2_angle<<endl;
          if(dcarm==1)
          {
           bbcs_v2_west_angle[iangle]->Fill(pt_angle,cosbbc_dphi2_angle);
         }
         else if( dcarm==0)
         {
           bbcs_v2_east_angle[iangle]->Fill(pt_angle,cosbbc_dphi2_angle);
         }
         bbcs_v2_incl_angle[iangle]->Fill(pt_angle,cosbbc_dphi2_angle);
        //}
       }
     }

   }
    }//end of event 
    
    htree->Delete(); 
    f->Close(); 
    delete f; 
    //}

    if(calFlag<3 && calFlag>0){
    //if(calFlag==1) sprintf(name,"tsukuba/flattening_%d.dat",runNumber);
    //if(calFlag==2) 
      sprintf(name,"vtx_ep_calib/new_wide_z_flattening/flattening_%d_%d.dat",runNumber,calFlag);
      cout << "writing calibration file : " << name << endl;
      ofstream ofs;
      ofs.open(name);

      cout << "flt para" << endl;
      for (int ic=0; ic<NMUL; ic++) {
        for (int iz=0; iz<NZPS; iz++) {
         for (int ih=0; ih<NHAR; ih++) {
           for (int id=0; id<NDET; id++) {
             for (int ib=0; ib<2; ib++) {
               ofs << ave[ic][iz][ih][id]->GetBinContent(ib+1) << " ";
               ofs << ave[ic][iz][ih][id]->GetBinError  (ib+1) << " ";
             }
             ofs << endl;
             for (int ib=0; ib<2; ib++) {
               for (int io=0; io<NORD; io++) {
                ofs << flt[ic][iz][ih][id]->GetBinContent(ib*NORD+io+1) << " ";
              }
              ofs << endl;
            }
          }
        }
      }
    }
    ofs.close();

    sprintf(name,"vtx_ep_calib/new_wide_z_flattening/flattening_bbc_angle_%d_%d.dat",runNumber,calFlag);
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

 if(calFlag>0){
  TFile *mData2=new TFile(outFile2,"recreate");
  cout<<outFile2<<endl;
  mData2->cd();
  for (int ic=0; ic<NMUL; ic++) {
    for (int iz=0; iz<NZPS; iz++) {
     for (int ih=0; ih<NHAR; ih++) {
       for (int id=0; id<NDET; id++) {
         ave[ic][iz][ih][id]->Write();
         flt[ic][iz][ih][id]->Write();

       }
     }
   }
 }

 for (int ic=0; ic<NMUL; ic++) {
  for (int ih=0; ih<NHAR; ih++) {
   for (int id=0; id<NDET; id++) {
     qx[ic][ih][id]->Write();
     qy[ic][ih][id]->Write();

     dis[ic][ih][id]->Write();
   }
 }
}

mData2->Close();
}

if(calFlag>0){
  cout<<outFile1<<endl;
  TFile *mData1=new TFile(outFile1,"recreate");
  mData1->cd();

  bbcs_v2_west->Write();
  bbcs_v2_east->Write();
  bbcs_v2_incl->Write();

  for(int iangle = 0; iangle < n_bbc_angles; iangle++)
  {
    bbcs_v2_incl_angle[iangle]->Write();
    bbcs_v2_west_angle[iangle]->Write();
    bbcs_v2_east_angle[iangle]->Write();
  }
    //end of adding


  if(calFlag< 3)
  {
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


  mData1->Close();


}


cout<<"cleaning up"<<endl;

delete bbcs_v2_incl;
delete bbcs_v2_east;
delete bbcs_v2_west;

for (int ic=0; ic<NMUL; ic++) {
  for (int iz=0; iz<NZPS; iz++) {
    for (int ih=0; ih<NHAR; ih++) {
     for (int id=0; id<NDET; id++) {
       delete ave[ic][iz][ih][id];

       delete flt[ic][iz][ih][id];


     }
   }
 }
}

for (int ic=0; ic<NMUL; ic++) {
  for (int ih=0; ih<NHAR; ih++) {
    for (int id=0; id<NDET; id++) {
     delete dis[ic][ih][id];

     delete qx[ic][ih][id];

     delete qy[ic][ih][id];
   }
 }
}

for (int ic=0; ic<10; ic++) {
  for (int id=0; id<2; id++) {
    delete disnew[ic][id];

    delete disnew2[ic][id];
  }
}

for(int iangle = 0; iangle < n_bbc_angles; iangle++)
{
  delete bbcs_v2_incl_angle[iangle];
  delete bbcs_v2_west_angle[iangle];
  delete bbcs_v2_east_angle[iangle];
}

for(int iangle = 0; iangle < n_bbc_angles; iangle++)
{
  for (int iz=0; iz<NZPS; iz++) {
   if(iz==0)
   {
     delete bbc_angle_dis[iangle];
     delete bbc_angle_qx[iangle];
     delete bbc_angle_qy[iangle];
   }

   delete bbc_angle_flt[iz][iangle];
   delete bbc_angle_ave[iz][iangle];
 }
}

cout<<"end of program ana"<<endl;
return;

}



/*
float Pc3Sdphi(int dcarm, int charge, float pt, float dphi){
  if(pt<=0.4||pt>=4.0) return -9999;
  float mean = pc3dphimean[dcarm][charge][0]+
    pc3dphimean[dcarm][charge][1]*pt+
    pc3dphimean[dcarm][charge][2]/sqrt(pt)+
    pc3dphimean[dcarm][charge][3]/pt/pt;

  float sigma = pc3dphisigma[dcarm][charge][0]+
    pc3dphisigma[dcarm][charge][1]*pt+
    pc3dphisigma[dcarm][charge][2]/sqrt(pt)+
    pc3dphisigma[dcarm][charge][3]/pt/pt;

  return fabs(dphi-mean)/sigma;
  
  else if(pt<0.5) return (dphi-0.000195674)/0.00958994;
  else if(pt<1.0) return (dphi-8.71664e-05)/0.00346979;
  else if(pt<1.5) return (dphi-0.000120889)/0.00248171;
  else if(pt<2.0) return (dphi-0.000111251)/0.00219466;
  else if(pt<2.5) return (dphi-0.000113817)/0.00210342;
  else if(pt<3.0) return (dphi-9.21377e-05)/0.00211282;
  else if(pt<3.5) return (dphi-0.000101712)/0.0021572;
  else if(pt<4.0) return (dphi-4.31527e-05)/0.00227459;
  else if(pt<4.5) return (dphi+1.95281e-06)/0.00269572;
  else if(pt<5.0) return (dphi+0.000112822)/0.00325696;
  else return -9999;
  
}

float Pc3Sdz(int dcarm, int charge, float pt, float dz){
  if(pt<=0.4||pt>=4.0) return -9999;
  float mean = pc3dzmean[dcarm][charge][0]+
    pc3dzmean[dcarm][charge][1]*pt+
    pc3dzmean[dcarm][charge][2]/sqrt(pt)+
    pc3dzmean[dcarm][charge][3]/pt/pt;

  float sigma= pc3dzsigma[dcarm][charge][0]+
    pc3dzsigma[dcarm][charge][1]*pt+
    pc3dzsigma[dcarm][charge][2]/sqrt(pt)+
    pc3dzsigma[dcarm][charge][3]/pt/pt;

  return fabs(dz-mean)/sigma;
  
  else if(pt<0.5) return (dz-0.0120409)/3.60411;
  else if(pt<1.0) return (dz+0.0163131)/2.45785;
  else if(pt<1.5) return (dz+0.0116410)/2.16330;
  else if(pt<2.0) return (dz+0.0027346)/1.97138;
  else if(pt<2.5) return (dz+0.0197661)/1.96495;
  else if(pt<3.0) return (dz-0.0777760)/1.86931;
  else if(pt<5.0) return (dz+0.0147089)/2.03005;
  else return -9999;
  
}
*/

void boost_and_rotate(TLorentzVector & vec, TLorentzVector input1, TLorentzVector input2){
	
  TLorentzVector cms = input1 + input2;

  TVector3 z(0,0,1);
  //TLorentzVector z(0,0,1,1);

  input1.Boost(-cms.BoostVector());
  input2.Boost(-cms.BoostVector());
  vec.Boost(-cms.BoostVector());
  
  float rotAngleY = -input1.Angle(z);
  
  //float rotAngleY = -input1.Angle(z.Vect());

  //Now rotate the beams about x to align them with z axis
  vec.RotateY(rotAngleY);
  
}

