#ifndef __SIMPLEFLOWTREEBBCFVTX_H__
#define __SIMPLEFLOWTREEBBCFVTX_H__

// standard includes
#include <string>
#include <vector>

// ROOT includes
#include <TFile.h>
#include <TNtuple.h>
//#include <TH1.h>
//#include <TH1D.h>
//#include <TH2.h>
//#include <TH2D.h>
#include <TProfile2D.h>

// PHENIX includes
#include <SubsysReco.h>
#include <PHParticle.h>
//#include <UltraLight/UltraLight.h>
//#include <UltraLight/UltraLightTrack.h>
#include <SvxSegment.h>


// global variables for array sizes
//static const int N_FVTX_CLUSTER_MAX = 750; // maximum number of fvtx clusters, may need to be changed for larger systems
static const int N_FVTX_CLUSTER_MAX = 4000; // what's the largest possible number?
static const int N_CTRK_MAX = 75; // max number of PHCentralTracks
static const int N_FTRK_MAX = 75; // max number of FVTX tracks

class Fun4AllHistoManager;
class BbcCalib;
class BbcGeo;
class TVector3;
class TLorentzVector;
class dAuBES_utils;


class SimpleFlowTreeBBCFVTX: public SubsysReco
{
 public:
  SimpleFlowTreeBBCFVTX();
  virtual ~SimpleFlowTreeBBCFVTX();

  /// Fun4All calls...
  int  Init         (PHCompositeNode *topNode);
  int  InitRun      (PHCompositeNode *topNode);
  int  process_event(PHCompositeNode *topNode);
  int  ResetEvent   (PHCompositeNode *topNode);
  int  End          (PHCompositeNode *topNode);
  int  EndRun       (PHCompositeNode *topNode);
  void Verbosity    (int verbosity) {_verbosity = verbosity;}

  /// Single particle ntuple output...
  void set_output_filename(std::string filename) {_output_filename = filename;} // select output file name externally
  void set_use_runlist(bool b){ _use_runlist = b;} // text file for runs to analyze
  void set_runlist_file(std::string filename) {_runlist_filename = filename;} // name of file for above
  //void set_calib_file(std::string calibfile){_calib_file = calibfile;}
  void set_create_ttree(bool b){_create_ttree = b;} // ??


  // --- set methods for variousthings to write out to trees/ntuples
  void set_write_bbc(bool b){_write_bbc = b;}
  void set_write_cnt(bool b){_write_cnt = b;}
  void set_write_fvtx(bool b){_write_fvtx = b;}
  void set_write_fvtx_clusters(bool b){_write_fvtx_clusters = b;}
  bool is_run_in_list(int runnumber);

 protected:

  // ---

  bool pass_eta_cut(float eta, int bbcz_bin);

  /// current event
  unsigned long _ievent;

  /// verbosity level
  int _verbosity;

  //changes what is being written out
  bool _create_ttree;
  //histogram vector
  Fun4AllHistoManager* HistoManager;

  bool _write_bbc;
  bool _write_cnt;
  bool _write_fvtx;
  bool _write_fvtx_clusters;
  /// module output filename
  std::string _output_filename;
  TFile *_output_file;
  //std::string _calib_file;

  //run list stuff
  bool _use_runlist;
  std::string _runlist_filename;

  /// module output ntuple
  TTree *_ntp_event;
  TNtuple *_ntp_cluster;

  //BBC objects to recalibrate BBCs event plane
  BbcCalib *m_bbccalib;
  BbcGeo   *m_bbcgeo;

  dAuBES_utils* _utils;            ///< Utilities class


  int tmp_evt;


  //-- ntp_event variables
  int event;
  float bbc_z;
  float centrality;
  int npc1;
  unsigned int trigger_scaled;
  unsigned int trigger_live;
  float d_Qx[9]; // x comp of Q vec, 3 detectors times 3 harmonics
  float d_Qy[9]; // y comp of Q vec, 3 detectors times 3 harmonics
  float d_Qw[9]; // weight for q-vector (not harmonic dependent? but the option is available)
  float bc_x;
  float bc_y;
  float vtx_z;
  float FVTX_X;
  float FVTX_Y;
  float FVTX_Z;

  //float d_BBCs_Qy[221];
  //float d_BBCs_Qw[221];

  // if (_write_bbc)
  float bbc_qn;
  float bbc_qs;
  float d_BBC_charge[128]; // charge from each bbc tube for south (first 64) and north (second 64)

  // if (_write_fvtx_clusters)
  int   d_nFVTX_clus;
  int   d_nFVTXN_clus;
  int   d_nFVTXS_clus;
  float d_FVTX_x[N_FVTX_CLUSTER_MAX];
  float d_FVTX_y[N_FVTX_CLUSTER_MAX];
  float d_FVTX_z[N_FVTX_CLUSTER_MAX];

  //  
  int d_ntrk;
  float d_cntpx[N_CTRK_MAX];
  float d_cntpy[N_CTRK_MAX];
  float d_cntpz[N_CTRK_MAX];
  float d_cntpc3sdz[N_CTRK_MAX];
  float d_cntpc3sdphi[N_CTRK_MAX];

  // if (_write_fvtx)
  int ntracklets;
  float fphi[N_FTRK_MAX];
  float feta[N_FTRK_MAX];
  float fthe[N_FTRK_MAX];
  float fchisq[N_FTRK_MAX];
  int farm[N_FTRK_MAX];
  int fnhits[N_FTRK_MAX];
  float fDCA_X[N_FTRK_MAX];
  float fDCA_Y[N_FTRK_MAX];


  //-- Other variables



};

#endif /* __VTXEVENTPLANERECO_H__ */



