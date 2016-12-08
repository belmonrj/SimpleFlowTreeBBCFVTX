#include "SimpleFlowTreeBBCFVTX.h"

#include <cmath>
#include <algorithm>
#include <stdlib.h>
#include <iostream>
#include <fcntl.h>
#include <fstream>


#include "TOAD.h"

#include <TMath.h>
#include <TString.h>

#include <TLorentzVector.h>
#include <TVector3.h>

// Fun4All includes...
#include <Fun4AllReturnCodes.h>
#include <Fun4AllServer.h>
#include <Fun4AllHistoManager.h>
#include <PHCompositeNode.h>
#include <PHIODataNode.h>
#include <phool.h>
#include <getClass.h>
#include <PHGlobalv9.h>
#include <EventHeader.h>
#include <TrigLvl1.h>
#include <SvxSegmentList.h>
#include <SvxSegment.h>
#include <VtxOut.h>
#include <PHPoint.h>
#include <PHAngle.h>
#include <RunHeader.h>
//#include <BbcLL1Out.h>
#include <TFvtxCompactTrkMap.h>
#include <BbcRaw.h>
#include <BbcOut.h>
#include <BbcCalib.hh>
#include <BbcGeo.hh>

#include "SvxClusterList.h"
#include "SvxCluster.h"

#include "RpSnglSumXY.h"
#include "RpSumXYObject.h"
#include "ReactionPlaneObject.h"
#include "RpConst.h"
#include "recoConsts.h"
#include <TFvtxCompactCoordMap.h>

#include <PreviousEvent.h>

#include "dAuBES_utils.h"


// ------------------------
#include "PHCentralTrack.h"
#include "PHSnglCentralTrack.h"



using namespace std;


// --- class constructor
SimpleFlowTreeBBCFVTX::SimpleFlowTreeBBCFVTX():
  SubsysReco("SIMPLEFLOWTREEBBCFVTX"),
  _ievent(0),
  _verbosity(0),
  _create_ttree(true),
  _write_bbc(false),
  _write_cnt(false),
  _write_fvtx(false),
  _write_fvtx_clusters(false),
  _output_filename("NULL"),
  _output_file(NULL),
  _use_runlist(false),
  _runlist_filename(""),
  _ntp_event(NULL),
  _utils(NULL),
  tmp_evt(0)
{
  ResetEvent(NULL);

  m_bbccalib = new BbcCalib();
  m_bbcgeo   = new BbcGeo();
  return;
}


// --- class destructor
SimpleFlowTreeBBCFVTX::~SimpleFlowTreeBBCFVTX()
{
  delete m_bbccalib;
  delete m_bbcgeo;
  if ( _utils ) delete _utils;
}


// --- Init method, part of Fun4All inheriance
int SimpleFlowTreeBBCFVTX::Init(PHCompositeNode *topNode)
{

  ResetEvent(topNode); // is this needed?

  if (_verbosity > 1) cout << PHWHERE << "::Init() - entered." << endl;

  _output_file = new TFile(_output_filename.c_str(), "RECREATE");

  if (_create_ttree)
  {
    _ntp_event = new TTree("ntp_event", "event-wise ntuple");
    _ntp_event->SetAutoFlush(1000);
    _ntp_event->SetMaxTreeSize(100000000000LL);
    _ntp_event -> Branch("event", &event, "event/F");
    _ntp_event -> Branch("bbc_z", &bbc_z, "bbc_z/F");
    _ntp_event -> Branch("centrality", &centrality, "centrality/F");
    _ntp_event -> Branch("npc1", &npc1, "npc1/I");
    _ntp_event -> Branch("trigger_scaled", &trigger_scaled, "trigger_scaled/i");
    _ntp_event -> Branch("trigger_live", &trigger_live, "trigger_live/i");
    _ntp_event -> Branch("d_Qx", &d_Qx, "d_Qx[9]/F");
    _ntp_event -> Branch("d_Qy", &d_Qy, "d_Qy[9]/F");
    _ntp_event -> Branch("d_Qw", &d_Qw, "d_Qw[9]/F");
    _ntp_event -> Branch("bc_x", &bc_x, "bc_x/F");
    _ntp_event -> Branch("bc_y", &bc_y, "bc_y/F");
    _ntp_event -> Branch("vtx_z", &vtx_z, "vtx_z/F");
    _ntp_event -> Branch("fvtx_x", &FVTX_X, "fvtx_x/F");
    _ntp_event -> Branch("fvtx_y", &FVTX_Y, "fvtx_y/F");
    _ntp_event -> Branch("fvtx_z", &FVTX_Z, "fvtx_z/F");
    if (_write_bbc)
    {
      _ntp_event -> Branch("bbc_qn", &bbc_qn, "bbc_qn/F");
      _ntp_event -> Branch("bbc_qs", &bbc_qs, "bbc_qs/F");
      _ntp_event -> Branch("d_BBC_charge", &d_BBC_charge, "d_BBC_charge[128]/F");
      //_ntp_event -> Branch("d_BBCs_Qy",&d_BBCs_Qy,"d_BBCs_Qy[221]/F");
      //_ntp_event -> Branch("d_BBCs_Qw",&d_BBCs_Qw,"d_BBCs_Qw[221]/F");
    }
    if (_write_fvtx_clusters)
    {
      _ntp_event -> Branch("d_nFVTX_clus", &d_nFVTX_clus, "d_nFVTX_clus/I");
      _ntp_event -> Branch("d_FVTX_x", &d_FVTX_x, "d_FVTX_x[d_nFVTX_clus]/F");
      _ntp_event -> Branch("d_FVTX_y", &d_FVTX_y, "d_FVTX_y[d_nFVTX_clus]/F");
      _ntp_event -> Branch("d_FVTX_z", &d_FVTX_z, "d_FVTX_z[d_nFVTX_clus]/F");
      _ntp_event -> Branch("d_nFVTXN_clus", &d_nFVTXN_clus, "d_nFVTXN_clus/I");
      _ntp_event -> Branch("d_nFVTXS_clus", &d_nFVTXS_clus, "d_nFVTXS_clus/I");
    }
    //now track based arrays
    if (_write_cnt)
    {
      _ntp_event -> Branch("d_ntrk", &d_ntrk, "d_ntrk/I");
      _ntp_event -> Branch("d_cntpx", &d_cntpx, "d_cntpx[d_ntrk]/F");
      _ntp_event -> Branch("d_cntpy", &d_cntpy, "d_cntpy[d_ntrk]/F");
      _ntp_event -> Branch("d_cntpz", &d_cntpz, "d_cntpz[d_ntrk]/F");
      _ntp_event -> Branch("d_cntcharge", &d_cntcharge, "d_cntcharge[d_ntrk]/F");
      _ntp_event -> Branch("d_cntpc3sdz", &d_cntpc3sdz, "d_cntpc3sdz[d_ntrk]/F");
      _ntp_event -> Branch("d_cntpc3sdphi", &d_cntpc3sdphi, "d_cntpc3sdphi[d_ntrk]/F");
    }
    if (_write_fvtx)
    {
      //fvtx tracking parameters
      //_ntp_event -> Branch("fvtx_z",&fvtx_z,"fvtx_z/F");
      _ntp_event -> Branch("ntracklets", &ntracklets, "ntracklets/I");
      _ntp_event -> Branch("feta", &feta, "feta[ntracklets]/F");
      _ntp_event -> Branch("fthe", &fthe, "fthe[ntracklets]/F");
      _ntp_event -> Branch("fphi", &fphi, "fphi[ntracklets]/F");
      _ntp_event -> Branch("fchisq", &fchisq, "fchisq[ntracklets]/F");
      _ntp_event -> Branch("farm", &farm, "farm[ntracklets]/I");
      _ntp_event -> Branch("fnhits", &fnhits, "fnhits[ntracklets]/I");
      _ntp_event -> Branch("fDCA_X", &fDCA_X, "fDCA_X[ntracklets]/F");
      _ntp_event -> Branch("fDCA_Y", &fDCA_Y, "fDCA_Y[ntracklets]/F");
    }

  }




  return EVENT_OK;
}


// --- InitRun, part of Fun4All inheritance
int SimpleFlowTreeBBCFVTX::InitRun(PHCompositeNode *topNode)
{

  int runnumber = 0;

  RunHeader *rh = findNode::getClass<RunHeader>(topNode, "RunHeader");
  if (!rh)
  {
    cout << PHWHERE << " ERROR::RunHeader not found" << endl;
    return ABORTEVENT;
  }
  else
  {
    runnumber = rh->get_RunNumber();
  }

  recoConsts *rc = recoConsts::instance();
  int icalibversion = 4002; // after Run5 and Field_ON

  rc->set_IntFlag("BBCCALIBVERSION", icalibversion);

  if (m_bbccalib) {
    PHTimeStamp TimeStp = rc->get_TimeStamp();
    int BBCCALIBVERSION = rc->get_IntFlag("BBCCALIBVERSION");
    cout << "SvxRpSumXYReco::InitRun - restored constants are for " << TimeStp << endl;
    m_bbccalib->restore(TimeStp, BBCCALIBVERSION);
  }


  // Setup the utility class
  // This is done in init run so that the collision system can be
  // determined from the run number
  TString _collsys = "Run16dAu200"; // default to 200 GeV
  // --- Run15pAu200
  if ( runnumber >= 432637 && runnumber <= 436647 )
    _collsys = "Run15pAu200";
  // --- Run15pAl200
  if ( runnumber >= 436759 && runnumber <= 438422 )
    _collsys = "Run15pAl200";
  // --- Run16dAu200
  if ( runnumber >= 454774 && runnumber <= 455639 )
    _collsys = "Run16dAu200";
  // --- Run16dAu62
  if ( runnumber >= 455792 && runnumber <= 456283 )
    _collsys = "Run16dAu62";
  // --- Run16dAu20
  if ( runnumber >= 456652 && runnumber <= 457298 )
    _collsys = "Run16dAu20";
  // --- Run16dAu39
  if ( runnumber >= 457634 && runnumber <= 458167 )
    _collsys = "Run16dAu39";

  // --- delete this pointer in EndRun
  _utils = new dAuBES_utils(_collsys, true);
  // _utils->is_sim(_is_sim);


  return EVENT_OK;
}


// --- ResetEvent, part of Fun4All inheritance, called after every event by Fun4All
int SimpleFlowTreeBBCFVTX::ResetEvent(PHCompositeNode *topNode)
{
  if (_verbosity > 1) cout << PHWHERE << "::ResetEvent() - entered." << endl;

  event         = -9999;
  centrality    = -9999;
  npc1          = -9999;
  trigger_scaled = -9999;
  trigger_live   = -9999;
  bbc_qn        = -9999;
  bbc_qs        = -9999;
  bbc_z         = -9999;
  vtx_z         = -9999;
  bc_x          = -9999;
  bc_y          = -9999;
  ntracklets    = -9999;

  for (int i = 0; i < 9; i++)
  {
    d_Qx[i] = -9999;
    d_Qy[i] = -9999;
    d_Qw[i] = 0;
  }

  for (int i = 0; i < 128; i++)
  {
    d_BBC_charge[i] = 0.0;

  }

  for (int i = 0; i < N_FVTX_CLUSTER_MAX; i++)
  {
    d_FVTX_x[i] = -9999;
    d_FVTX_y[i] = -9999;
    d_FVTX_z[i] = -9999;
  }

  if (_write_fvtx)
    for (int i = 0; i < 75; i++)
    {
      farm[i]      = -9999;
      fphi[i]      = -9999;
      feta[i]      = -9999;
      fthe[i]      = -9999;
      fchisq[i]    = -9999;
      fnhits[i]    = -9999;
      fDCA_X[i]    = -9999;
      fDCA_Y[i]    = -9999;
    }

  return EVENT_OK;

}



// --- process_event, inherited from Fun4All, does the main part of the analysis on an event by event basis
int SimpleFlowTreeBBCFVTX::process_event(PHCompositeNode *topNode)
{

  if (_verbosity > 1) cout << PHWHERE << "::process_event() - entered on event #" << _ievent << endl;
  else if ((_ievent) % 10000 == 0) cout << PHWHERE << "::process_event() - event #" << _ievent << endl;

  // advance event counter
  _ievent++;

  event = _ievent;



  //---------------------------------
  // Run Number Selection Capability
  //---------------------------------

  int runnumber = 0;

  RunHeader *rh = findNode::getClass<RunHeader>(topNode, "RunHeader");
  if (!rh)
  {
    cout << PHWHERE << " ERROR::RunHeader not found" << endl;
    return ABORTEVENT;
  }
  else
  {
    runnumber = rh->get_RunNumber();
  }

  if (_use_runlist)
  {
    if (is_run_in_list(runnumber) == false)
    {
      if (_verbosity > 1) {
        cout << endl << "ABORTING RUN, Run number: " << runnumber << " is not in the run list" << endl << endl;
      }
      return ABORTEVENT;
    }
  }


  //-------------------------------
  // Grab info off of the node tree
  //-------------------------------

  TrigLvl1 *triggers = findNode::getClass<TrigLvl1>(topNode, "TrigLvl1");
  if (!triggers)
  {
    cout << PHWHERE << " ERROR::TrigLvl1 not found" << endl;
    return ABORTEVENT;
  }
  PHGlobal *global = findNode::getClass<PHGlobal>(topNode, "PHGlobal");
  if (!global)
  {
    cout << PHWHERE << " ERROR::PHGlobal not found" << endl;
    return ABORTEVENT;
  }
  EventHeader *evthead = findNode::getClass<EventHeader>(topNode, "EventHeader");
  if (!evthead)
  {
    cout << PHWHERE << " ERROR::EventHeader not found" << endl;
    return ABORTEVENT;
  }
  VtxOut *vertexes = findNode::getClass<VtxOut>(topNode, "VtxOut");
  if (!vertexes)
  {
    cout << PHWHERE << " ERROR::VtxOut not found" << endl;
    return ABORTEVENT;
  }

  TFvtxCompactTrkMap* trkfvtx_map = NULL;
  if (_write_fvtx)
  {
    trkfvtx_map = findNode::getClass<TFvtxCompactTrkMap>(topNode, "TFvtxCompactTrkMap");
    if (!trkfvtx_map) {
      cerr << PHWHERE << " No TFvtxCompactTrkMap object !" << endl;
      return ABORTEVENT;
    }
  }

  TFvtxCompactCoordMap* fvtx_coord_map = findNode::getClass<TFvtxCompactCoordMap>(topNode, "TFvtxCompactCoordMap");
  if (!fvtx_coord_map && _write_fvtx_clusters)
  {
    cout << PHWHERE << " ERROR::TFvtxCompactCoordMap not found" << endl;
    return ABORTEVENT;
  }

  RpSumXYObject* d_rp = findNode::getClass<RpSumXYObject>(topNode, "RpSumXYObject");
  if (!d_rp)
  {
    if ( _verbosity > 4 ) cout << PHWHERE << "Could not find the RPSumXYObject" << endl;
    //return ABORTEVENT;
  }

  BbcRaw *bbcraw = findNode::getClass<BbcRaw>(topNode, "BbcRaw");
  if (!bbcraw)
  {
    cout << PHWHERE << "Could not find Bbcraw!" << endl;
    if (_write_bbc)
      return ABORTEVENT;
  }


  //---------------------------------------------------------//
  //
  //         Make Event Selection
  //
  //---------------------------------------------------------//

  if (!_utils->is_event_ok(topNode))
    return EVENT_OK;



  //---------------------------------------------------------//
  //
  //         Reading in Global Event Information into Tree
  //
  //---------------------------------------------------------//

  PHPoint precise_vertex1 = vertexes->get_Vertex("SVX_PRECISE");
  vtx_z = precise_vertex1.getZ();
  if (vtx_z != vtx_z) vtx_z = -9999; //NAN check

  PHPoint svx_fast = vertexes->get_Vertex("SVX");//seed vertex
  bc_x = svx_fast.getX();//these are actually the beam center
  bc_y = svx_fast.getY();

  // phglobal fields
  centrality  = global->getCentrality();
  bbc_qn      = global->getBbcChargeN();
  bbc_qs      = global->getBbcChargeS();
  npc1        = global->getNumberPC1Hits();
  event = evthead->get_EvtSequence();
  trigger_scaled = triggers->get_lvl1_trigscaled();
  trigger_live = triggers->get_lvl1_triglive();

  // // --- these numbers taken from run 16 run control log
  // unsigned int trigger_FVTXNSBBCScentral = 0x00100000;
  // unsigned int trigger_FVTXNSBBCS        = 0x00400000;
  // unsigned int trigger_BBCLL1narrowcent  = 0x00000008;
  // unsigned int trigger_BBCLL1narrow      = 0x00000010;

  // unsigned int accepted_triggers = 0;
  // // accepted_triggers = trigger_FVTXNSBBCScentral | trigger_FVTXNSBBCS | trigger_BBCLL1narrowcent | trigger_BBCLL1narrow ;
  // // --- Run16dAu200
  // if ( runnumber >= 454774 && runnumber <= 455639 ) accepted_triggers = trigger_BBCLL1narrowcent | trigger_BBCLL1narrow;
  // // --- Run16dAu62
  // if ( runnumber >= 455792 && runnumber <= 456283 ) accepted_triggers = trigger_BBCLL1narrowcent | trigger_BBCLL1narrow;
  // // --- Run16dAu20
  // if ( runnumber >= 456652 && runnumber <= 457298 ) accepted_triggers = trigger_FVTXNSBBCScentral | trigger_FVTXNSBBCS;
  // // --- Run16dAu39
  // if ( runnumber >= 457634 && runnumber <= 458167 ) accepted_triggers = trigger_FVTXNSBBCScentral | trigger_FVTXNSBBCS;



  // --- bbc_z...
  PHPoint vertex1 = vertexes->get_Vertex("BBC");
  bbc_z = vertex1.getZ();
  if ( bbc_z != bbc_z ) bbc_z = -9999; // reassign nan

  PHPoint fvtx_vertex = vertexes->get_Vertex("FVTX");
  FVTX_X = fvtx_vertex.getX();
  FVTX_Y = fvtx_vertex.getY();
  FVTX_Z = fvtx_vertex.getZ();
  if ( FVTX_Z != FVTX_Z ) FVTX_Z = -9999; // reassign nan

  // cout << endl;
  // cout << "--- starting vertex checking ---" << endl;
  float zvtx = _utils->get_vrtx(topNode);


  if ( _verbosity > 1 ) cout << "FVTX vertex points: " << FVTX_X << " " << FVTX_Y << " " << FVTX_Z << endl;




  //int ibbcz_bin = (bbc_z+30.0)/10;//for fvtx eta cuts


  //---------------------------------------------------------//
  //  Finished Reading in Global Event Information into Tree
  //---------------------------------------------------------//



  //if(_write_fvtx && !trkfvtx_map) return ABORTEVENT;

  for (int i = 0; i < 9; i++) { //fvtx bbc cnt smd vtx 3+3+1+3+%d
    d_Qx[i] = -9999;
    d_Qy[i] = -9999;
    d_Qw[i] = 0.0;
  }

  float sumxy_fvtx[3][8][3];
  for (int ihar = 0; ihar < 3; ihar++) { //fvtx bbc cnt smd vtx 3+3+1+3+%d
    for (int idet = 0; idet < 8; idet++) {
      for (int ixy = 0; ixy < 3; ixy++) {
        sumxy_fvtx[ihar][idet][ixy] = 0.0;
      }
    }
  }

  //---------------------------------------------------------//
  //
  //                  Retrieving Event Planes
  //
  //---------------------------------------------------------//

  for (int idet = 62; idet < 68; idet++) {
    for (int ihar = 0; ihar < 3; ihar++) {
      int idcode = -9999;

      if (idet < 65) idcode = RP::calcIdCode(RP::ID_MPC, idet - 62, ihar);
      else idcode = RP::calcIdCode(RP::ID_BBC, idet - 65, ihar);

      if (idcode >= 0 && d_rp) {
        RpSnglSumXY *s_rp = d_rp->getRpSumXY(idcode);

        if (idet == 62) { //MPC
          d_Qx[ihar * 3 + idet - 62] = s_rp->QVector(0);
          d_Qy[ihar * 3 + idet - 62] = s_rp->QVector(1);
          d_Qw[ihar * 3 + idet - 62] = s_rp->Weight();
        }
        if (idet == 65) { //BBC
          d_Qx[ihar * 3 + idet - 65 + 2] = s_rp->QVector(0);
          d_Qy[ihar * 3 + idet - 65 + 2] = s_rp->QVector(1);
          d_Qw[ihar * 3 + idet - 65 + 2] = s_rp->Weight();
        }
      }
    }
  }
  //FVTX
  //cout<<RP::ID_FVT<<endl;
  //south

  for (int idet = 0; idet < 20; idet++) {
    for (int ihar = 0; ihar < 3; ihar++) {
      int idcode = RP::calcIdCode(RP::ID_FVT, idet, ihar);

      if (idcode >= 0 && d_rp) {
        RpSnglSumXY *s_rp = d_rp->getRpSumXY(idcode);
        int ifvtx = idet % 5;
        if (ifvtx < 4) { //-1.0-3.0
          sumxy_fvtx[ihar][0][0] += s_rp->QVector(0); //compiling the 20 FVTX event planes
          sumxy_fvtx[ihar][0][1] += s_rp->QVector(1);
          sumxy_fvtx[ihar][0][2] += s_rp->Weight();
        }
      }
    }
  }

  for (int ihar = 0; ihar < 3; ihar++) {
    //for(int idet=0; idet<8; idet++){
    d_Qx[ihar * 3 + 1] = sumxy_fvtx[ihar][0][0];
    d_Qy[ihar * 3 + 1] = sumxy_fvtx[ihar][0][1];
    d_Qw[ihar * 3 + 1] = sumxy_fvtx[ihar][0][2];
    //}
  }


  //---------------------------------------------------------//
  //             Finished retrieving Event Planes
  //---------------------------------------------------------//

  //---------------------------------------------------------//
  //
  //            Writing out BBC Raw Objects
  //
  //---------------------------------------------------------//

  if (_write_bbc)
  {
    if (bbcraw)
    {
      for ( int ipmt = 0; ipmt < 128; ++ipmt )
      {

        short iadc    = bbcraw->get_Adc(ipmt);
        short itdc    = bbcraw->get_Tdc0(ipmt);
        //float tdc1=bbcraw->get_Tdc1(ipmt);
        float itime0  = m_bbccalib->getHitTime0(ipmt, itdc, iadc);
        float icharge = m_bbccalib->getCharge(ipmt, iadc);
        if ( ( itime0 <= 0 || icharge <= 0 ) && verbosity > 8 )
        {
          // --- this happens A LOT
          cout << "EITHER SIDE!!!" << endl;
          cout << "For event number " << _ievent - 1
               << " and BBC tube number " << ipmt
               << " the time is " << itime0
               << " and the charge is " << icharge
               << endl;
        }
        if ( icharge <= 0 && itime0 >= 0 )
        {
          // --- this never happens
          cout << "THIS SHOULDN'T BE HAPPENING" << endl;
          cout << "For event number " << _ievent - 1
               << " and BBC tube number " << ipmt
               << " the time is " << itime0
               << " and the charge is " << icharge
               << endl;
        }
        if ( itime0 <= 0 ) icharge -= 10000.0;
        //float ibbc_x  = m_bbcgeo->getX(ipmt);
        //float ibbc_y  = m_bbcgeo->getY(ipmt);
        //float ibbc_z  = m_bbcgeo->getZ(ipmt);
        //cout<<"d_pmt_x["<<ipmt<<"] = "<<ibbc_x<<";"<<endl;
        //cout<<"d_pmt_y["<<ipmt<<"] = "<<ibbc_y<<";"<<endl;
        //cout<<"d_pmt_z["<<ipmt<<"] = "<<ibbc_z<<";"<<endl;

        d_BBC_charge[ipmt] = icharge;

      }//end of ipmt loop
    }
  }

  //---------------------------------------------------------//
  //                 finished Get BBC Raw
  //---------------------------------------------------------//

  //---------------------------------------------------------//
  //
  //            Get FVTX clusters
  //
  //---------------------------------------------------------//

  int nfvtx_raw_clus = 0;
  int nfvtxn_raw_clus = 0;
  int nfvtxs_raw_clus = 0;
  if ( fvtx_coord_map )
  {
    TFvtxCompactCoordMap::iterator _iter( fvtx_coord_map->range() );
    while ( TFvtxCompactCoordMap::const_pointer fvtx_ptr = _iter.next() )
    {
      TFvtxCompactCoord* fvtxcoord = fvtx_ptr->get();
      PHPoint fvtx_coord_point = fvtxcoord->get_coord_midpoint();
      //int iarm = fvtxcoord->get_arm();
      float fvtx_x = fvtx_coord_point.getX();
      float fvtx_y = fvtx_coord_point.getY();
      float fvtx_z = fvtx_coord_point.getZ();
      //float fvtx_r = sqrt(pow(fvtx_x,2.0)+pow(fvtx_y,2.0));
      if ( (fabs(fvtx_x) > 999) || (fabs(fvtx_y) > 999) || (fabs(fvtx_z) > 999)) continue;
      //float fvtx_the = atan2(fvtx_r,fvtx_z-vtx_z);
      //float fvtx_phi = atan2(fvtx_y,fvtx_x);
      //float fvtx_eta = -log(tan(0.5*fvtx_the));
      //if(fvtx_z < 0)
      //{
      if (nfvtx_raw_clus >= N_FVTX_CLUSTER_MAX)
      {
        cout << "butting against the max fvtx cluster size " << nfvtx_raw_clus << "/" << N_FVTX_CLUSTER_MAX << ", breaking" << endl;
        break;
      }
      d_FVTX_x[nfvtx_raw_clus] = fvtx_x;
      d_FVTX_y[nfvtx_raw_clus] = fvtx_y;
      d_FVTX_z[nfvtx_raw_clus] = fvtx_z;
      ++nfvtx_raw_clus;
      if ( fvtx_z > 0 ) ++nfvtxn_raw_clus;
      if ( fvtx_z < 0 ) ++nfvtxs_raw_clus;
      //cout<<"fvtx_eta: "<<fvtx_eta<<endl;
      //} // south

    } // while loop over iterator
  } // check on fvtx_coord_map

  d_nFVTX_clus = nfvtx_raw_clus;
  d_nFVTXN_clus = nfvtxn_raw_clus;
  d_nFVTXS_clus = nfvtxs_raw_clus;
  //cout<<"nfvtxs_raw_clus: "<<nfvtxs_raw_clus <<endl;

  //---------------------------------------------------------//
  //                 finished Get FVTX Clusters
  //---------------------------------------------------------//

  //---------------------------------------------------------//
  //
  //                 Get FVTX Tracks
  //
  //---------------------------------------------------------//

  //int ntr = -1;
  int ntr = 0;

  if ( trkfvtx_map && _write_fvtx )
  {
    TFvtxCompactTrkMap::const_iterator trk_iter = trkfvtx_map->range();
    while ( TFvtxCompactTrkMap::const_pointer trk_ptr = trk_iter.next() )
    {

      TFvtxCompactTrk* fvtx_trk = trk_ptr->get();

      //-- Only write out good fvtx tracks
      if ( !_utils->is_fvtx_track_ok(fvtx_trk, zvtx) )
        continue;

      float the = fvtx_trk->get_fvtx_theta();
      float eta = fvtx_trk->get_fvtx_eta();
      float phi = fvtx_trk->get_fvtx_phi();
      int   arm = (int)fvtx_trk->get_arm();
      float fvtx_x      = fvtx_trk->get_fvtx_vtx().getX();
      float fvtx_y      = fvtx_trk->get_fvtx_vtx().getY();
      float fvtx_z      = fvtx_trk->get_fvtx_vtx().getZ();
      int   nfhits      = (int)fvtx_trk->get_nhits();

      // fix total momentum to 1.0 (for rotating due to beamtilt)
      double px = 1.0 * TMath::Sin(the) * TMath::Cos(phi);
      double py = 1.0 * TMath::Sin(the) * TMath::Sin(phi);
      double pz = 1.0 * TMath::Cos(the);

      // rotate based on beamtilt
      px = _utils->rotate_x(px, pz);
      pz = _utils->rotate_z(px, pz);
      phi = TMath::ATan2(py, px);
      the = TMath::ACos(pz / TMath::Sqrt(px * px + py * py + pz * pz));


      float vertex_z = bbc_z;
      if ( FVTX_Z > -999 ) vertex_z = FVTX_Z;
      float DCA_x      = fvtx_x + tan(the) * cos(phi) * (vertex_z - fvtx_z);
      float DCA_y      = fvtx_y + tan(the) * sin(phi) * (vertex_z - fvtx_z);
      //cout<<"nfhits: "<<nfhits<<endl;

      //use this code of you want to figure out which layers have the hits
      /*      short nhits = 0;
              for (int i=0; i<8; i++)
              nhits += has_hit(i);
              for (int i=0; i<4; i++)
              nhits += has_svxhit(i);
              return nhits;
      */
      //if(the==0 || phi==0 || fvtx_x==0 || fvtx_y==0 || fvtx_z==0) continue;
      //if(the==0) continue;

      //if ( nfhits < 3 ) continue;
      //if ( !pass_eta_cut(eta,ibbcz_bin) ) continue;
      if ( fvtx_trk->get_chi2_ndf() > 5 ) continue;
      if ( fabs(DCA_x - 0.3) > 2.0 || fabs(DCA_y - 0.02) > 2.0 ) continue;


      //float DCA_R      = sqrt((DCA_x*DCA_x) + (DCA_y*DCA_y));

      if (ntr < N_FTRK_MAX)
      {
        feta[ntr]   = eta;
        fthe[ntr]   = the;
        fphi[ntr]   = phi;
        fchisq[ntr] = fvtx_trk->get_chi2_ndf();
        farm[ntr]   = arm;
        fnhits[ntr] = nfhits;
        fDCA_X[ntr] = DCA_x;
        fDCA_Y[ntr] = DCA_y;
      }
      else
      {
        cout << "butting up against the boundary of fvtx tracks" << endl;
        break;
      }
      //short_chi2 = fvtx_trk->get_short_chi2_ndf();
      //chi2       = fvtx_trk->get_chi2_ndf();

      ++ntr;

    } // end while loop over tracks
  } // check on fvtx track map

  ntracklets = ntr;

  //---------------------------------------------------------//
  //                 finished Get FVTX Tracks
  //---------------------------------------------------------//


  //---------------------------------------------------------//
  //
  //                 Get CNT Tracks
  //
  //---------------------------------------------------------//

  d_ntrk = 0;
  PHCentralTrack *ctrk = findNode::getClass<PHCentralTrack>(topNode, "PHCentralTrack");
  if ( ctrk )
  {
    int ntrk = ctrk->get_npart();
    if ( ntrk > N_CTRK_MAX )
    {
      cout << PHWHERE << " WARNING: too many tracks, skipping event" << endl;
      return ABORTEVENT;
    }
    int counter = 0;
    for ( int itrk = 0; itrk < ntrk; ++itrk)
    {

      PHSnglCentralTrack *strk = ctrk->get_track(itrk);

      //-- Only write out good PHCentralTracks
      if ( !_utils->is_cnt_track_ok(strk) )
        continue;


      float mom         = strk->get_mom();
      float zed         = strk->get_zed();
      int quality       = strk->get_quality();
      if ( mom < 0.0 || mom > 50.0 ) continue;
      if ( fabs(zed) < 3.0 || fabs(zed) > 70.0 ) continue;
      if ( quality != 63 && quality != 31 ) continue;

      float px = strk->get_px();
      float py = strk->get_py();
      float pz = strk->get_pz();

      // int arm = 0;
      // if ( px > 0 ) arm = 1;

      // rotate based on beamtilt
      px = _utils->rotate_x(px, pz);
      pz = _utils->rotate_z(px, pz);


      int charge = strk->get_charge();

      float pc3dphi = strk->get_pc3dphi();
      float pc3dz = strk->get_pc3dz();

      if ( pc3dphi < -9990 || pc3dz < -9990 ) continue;

      // float pc3sdphi = calcsdphi(pc3dphi, arm, charge, mom);
      // float pc3sdz = calcsdz(pc3dz, arm, charge, mom);
      float pc3sdphi = strk->get_pc3sdphi();
      float pc3sdz = strk->get_pc3sdz();

      if ( fabs(pc3sdphi) > 3.0 || fabs(pc3sdz) > 3.0 ) continue;

      d_cntpx[counter] = px;
      d_cntpy[counter] = py;
      d_cntpz[counter] = pz;
      d_cntcharge[counter] = charge;
      d_cntpc3sdz[counter] = pc3sdz;
      d_cntpc3sdphi[counter] = pc3sdphi;

      ++counter;

    } // loop over tracks

    d_ntrk = counter;

  } // check on track pointer

  //---------------------------------------------------------//
  //                 finished Get CNT Tracks
  //---------------------------------------------------------//


  if ( _create_ttree ) _ntp_event->Fill();

  if ( _verbosity > 0 ) cout << "sucessfully processed this event" << endl;

  tmp_evt++;//to keep track of how many events pass event cuts
  return EVENT_OK;

} // end of process_event



int SimpleFlowTreeBBCFVTX::EndRun(PHCompositeNode *topNode)
{
  if ( _utils ) delete _utils;
  return EVENT_OK;
}

int SimpleFlowTreeBBCFVTX::End(PHCompositeNode *topNode)
{
  if (_verbosity > 1) cout << PHWHERE << "::End() - entered." << endl;

  cout << "total events: " << _ievent << " fraction passing vtx cut: " << tmp_evt * 1.0 / _ievent << endl;

  _output_file->cd();

  if (_create_ttree)
    _ntp_event->Write();

  _output_file->Close();
  delete _output_file;

  //delete  m_bbccalib;
  //delete  m_bbcgeo;

  //cout << PHWHERE << "::End() - ntuple output written to file: " << _output_filename << endl;


  return EVENT_OK;
}




bool SimpleFlowTreeBBCFVTX::is_run_in_list(int runnumber)
{
  ifstream runlist;
  runlist.open(_runlist_filename.c_str());
  int run_num;
  while (1) {
    if (!runlist.good()) break;
    runlist >> run_num;
    if (run_num == runnumber)
      return true;
  }
  return false;
}

//these vaules were obtained from ana note 1162
bool SimpleFlowTreeBBCFVTX::pass_eta_cut(float eta, int bbcz_bin)
{
  if (bbcz_bin == 0)
  {
    if (eta > 1.9 && eta < 3.2)
      return true;
  }
  if (bbcz_bin == 1)
  {
    if (eta > -1.7 && eta < -0.8)
      return true;
    if (eta > 1.7 && eta < 3.2)
      return true;
  }
  if (bbcz_bin == 2)
  {
    if (eta > -2.1 && eta < -0.7)
      return true;
    if (eta > 1.5 && eta < 3.1)
      return true;
  }
  if (bbcz_bin == 3)
  {
    if (eta > -2.4 && eta < -0.9)
      return true;
    if (eta > 1.3 && eta < 3.0)
      return true;
  }
  if (bbcz_bin == 4)
  {
    if (eta > -2.6 && eta < -1.0)
      return true;
    if (eta > 1.2 && eta < 2.8)
      return true;
  }
  if (bbcz_bin == 5)
  {
    if (eta > -2.8 && eta < -1.2)
      return true;
    if (eta > 1.0 && eta < 2.6)
      return true;
  }
  if (bbcz_bin == 6)
  {
    if (eta > -3.0 && eta < -1.3)
      return true;
    if (eta > 0.9 && eta < 2.4)
      return true;
  }
  if (bbcz_bin == 7)
  {
    if (eta > -3.1 && eta < -1.4)
      return true;
    if (eta > 0.7 && eta < 2.1)
      return true;
  }
  if (bbcz_bin == 8)
  {
    if (eta > -3.2 && eta < -1.5)
      return true;
    if (eta > 0.7 && eta < 1.6)
      return true;
  }
  if (bbcz_bin == 9)
  {
    if (eta > -3.3 && eta < -1.8)
      return true;
  }

  return false;
}

