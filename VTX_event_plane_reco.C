#include "VTX_event_plane_reco.h"
#include "Run16RoughMatchingPC3.h"

#include <cmath>
#include <algorithm>
#include <stdlib.h>
#include <iostream>
#include <fcntl.h>
#include <fstream>


#include "TOAD.h"

#include <TMath.h>

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

// ------------------------
#include "PHCentralTrack.h"
#include "PHSnglCentralTrack.h"



using namespace std;


// --- class constructor
VTX_event_plane_reco::VTX_event_plane_reco(): SubsysReco("VTXULTRALIGHTRECO")
{
  _ievent = 0;
  tmp_evt = 0;
  _verbosity = 0;
  _z_vertex_range = 30;//cm
  _ntp_event = NULL;
  _output_filename = "NULL";
  _output_file = NULL;
  _create_ttree = true;
  _write_bbc = false;
  _trimmed_tree = false;
  _write_vtx = true;
  _write_fvtx_clusters = false;
  int file = open("/dev/urandom", O_RDONLY);
  int seed;
  read(file, &seed, sizeof seed);
  close(file);
  _write_clusters = false;
  srand(seed);
  cout << PHWHERE << " random seed = " << seed << endl;
  _write_fvtx = false;

  m_bbccalib = new BbcCalib();
  m_bbcgeo   = new BbcGeo();
  return;
}


// --- class destructor
VTX_event_plane_reco::~VTX_event_plane_reco()
{
  delete m_bbccalib;
  delete m_bbcgeo;
}


// --- Init method, part of Fun4All inheriance
int VTX_event_plane_reco::Init(PHCompositeNode *topNode)
{

  ResetEvent(topNode); // is this needed?

  if (_verbosity > 1) cout << PHWHERE << "::Init() - entered." << endl;

  _output_file = new TFile(_output_filename.c_str(), "RECREATE");

  if(_create_ttree)
    {
      _ntp_event = new TTree("ntp_event","event-wise ntuple");
      _ntp_event->SetAutoFlush(1000);
      _ntp_event -> Branch("event",&event,"event/F");
      _ntp_event -> Branch("bbc_z",&bbc_z,"bbc_z/F");
      _ntp_event -> Branch("centrality",&centrality,"centrality/F");
      _ntp_event -> Branch("trigger_scaled",&trigger_scaled,"trigger_scaled/i");
      _ntp_event -> Branch("trigger_live",&trigger_live,"trigger_live/i");
      _ntp_event -> Branch("d_Qx",&d_Qx,"d_Qx[9]/F");
      _ntp_event -> Branch("d_Qy",&d_Qy,"d_Qy[9]/F");
      _ntp_event -> Branch("d_Qw",&d_Qw,"d_Qw[9]/F");
      _ntp_event -> Branch("bc_x",&bc_x,"bc_x/F");
      _ntp_event -> Branch("bc_y",&bc_y,"bc_y/F");
      _ntp_event -> Branch("vtx_z",&vtx_z,"vtx_z/F");
      _ntp_event -> Branch("fvtx_x",&FVTX_X,"fvtx_x/F");
      _ntp_event -> Branch("fvtx_y",&FVTX_Y,"fvtx_y/F");
      _ntp_event -> Branch("fvtx_z",&FVTX_Z,"fvtx_z/F");
      if(_write_bbc)
        _ntp_event -> Branch("d_BBC_charge",&d_BBC_charge,"d_BBC_charge[64]/F");
      if(_write_fvtx_clusters)
        {
          _ntp_event -> Branch("d_nFVTX_clus",&d_nFVTX_clus,"d_nFVTX_clus/I");
          _ntp_event -> Branch("d_FVTX_x",&d_FVTX_x,"d_FVTX_x[d_nFVTX_clus]/F");
          _ntp_event -> Branch("d_FVTX_y",&d_FVTX_y,"d_FVTX_y[d_nFVTX_clus]/F");
          _ntp_event -> Branch("d_FVTX_z",&d_FVTX_z,"d_FVTX_z[d_nFVTX_clus]/F");
        }
      _ntp_event -> Branch("d_ntrk",&d_ntrk,"d_ntrk/I");
      _ntp_event -> Branch("d_cntpx",&d_cntpx,"d_cntpx[d_ntrk]/F");
      _ntp_event -> Branch("d_cntpy",&d_cntpy,"d_cntpy[d_ntrk]/F");
      _ntp_event -> Branch("d_cntpz",&d_cntpz,"d_cntpz[d_ntrk]/F");
      _ntp_event -> Branch("d_cntpc3sdz",&d_cntpc3sdz,"d_cntpc3sdz[d_ntrk]/F");
      _ntp_event -> Branch("d_cntpc3sdphi",&d_cntpc3sdphi,"d_cntpc3sdphi[d_ntrk]/F");
      //_ntp_event -> Branch("d_BBCs_Qy",&d_BBCs_Qy,"d_BBCs_Qy[221]/F");
      //_ntp_event -> Branch("d_BBCs_Qw",&d_BBCs_Qw,"d_BBCs_Qw[221]/F");
      //now track based arrays
      if(_write_vtx)
        {
          _ntp_event -> Branch("nsegments",&nsegments,"nsegments/I");
          _ntp_event -> Branch("px",&px,"px[nsegments]/F");
          _ntp_event -> Branch("py",&py,"py[nsegments]/F");
          _ntp_event -> Branch("pz",&pz,"pz[nsegments]/F");
        }
      if(!_trimmed_tree)
        {
          _ntp_event -> Branch("vtx_x",&vtx_x,"vtx_x/F");
          _ntp_event -> Branch("vtx_y",&vtx_y,"vtx_y/F");
          //          _ntp_event -> Branch("centrality",&centrality,"centrality/F");
          _ntp_event -> Branch("bbc_qn",&bbc_qn,"bbc_qn/F");
          _ntp_event -> Branch("bbc_qs",&bbc_qs,"bbc_qs/F");
          _ntp_event -> Branch("eventok",&eventok,"eventok/I");
          _ntp_event -> Branch("trackID",&trackID,"trackID[nsegments]/I");
          _ntp_event -> Branch("charge",&charge,"charge[nsegments]/I");
          _ntp_event -> Branch("chisq",&chisq,"chisq[nsegments]/F");
          _ntp_event -> Branch("ndf",&ndf,"ndf[nsegments]/I");
          _ntp_event -> Branch("nhit0",&nhit0,"nhit0[nsegments]/I");
          _ntp_event -> Branch("nhit1",&nhit0,"nhit0[nsegments]/I");
          _ntp_event -> Branch("nhit2",&nhit0,"nhit0[nsegments]/I");
          _ntp_event -> Branch("nhit3",&nhit0,"nhit0[nsegments]/I");
          _ntp_event -> Branch("dca",&dca,"dca[nsegments]/F");
          _ntp_event -> Branch("dca2d",&dca2d,"dca2d[nsegments]/F");
          _ntp_event -> Branch("segmentok",&segmentok,"segmentok[nsegments]/I");
          _ntp_event -> Branch("svtx_z",&svtx_z,"svtx_z/F");
          _ntp_event -> Branch("vtxposE_x",&vtxposE_x,"vtxposE_x/F");
          _ntp_event -> Branch("vtxposE_y",&vtxposE_y,"vtxposE_y/F");
          _ntp_event -> Branch("vtxposE_z",&vtxposE_z,"vtxposE_z/F");
          _ntp_event -> Branch("vtxposW_x",&vtxposW_x,"vtxposW_x/F");
          _ntp_event -> Branch("vtxposW_y",&vtxposW_y,"vtxposW_y/F");
          _ntp_event -> Branch("vtxposW_z",&vtxposW_z,"vtxposW_z/F");

        }
      if(_write_fvtx)
        {
          //fvtx tracking parameters
          //_ntp_event -> Branch("fvtx_z",&fvtx_z,"fvtx_z/F");
          _ntp_event -> Branch("ntracklets",&ntracklets,"ntracklets/I");
          _ntp_event -> Branch("feta",&feta,"feta[75]/F");
          _ntp_event -> Branch("fphi",&fphi,"fphi[75]/F");
          _ntp_event -> Branch("fchisq",&fchisq,"fchisq[75]/F");
          _ntp_event -> Branch("farm",&farm,"farm[75]/I");
          _ntp_event -> Branch("fnhits",&fnhits,"fnhits[75]/I");
          _ntp_event -> Branch("fDCA_X",&fDCA_X,"fDCA_X[75]/F");
          _ntp_event -> Branch("fDCA_Y",&fDCA_Y,"fDCA_Y[75]/F");
        }

      if(_write_clusters)
        _ntp_cluster = new TNtuple("ntp_cluster","cluster-wise ntuple",
                                   "event:layer:x:y:z:adc1:adc2:section:ladder:sensor:lx:ly:lz:size:sizex:sizez:vz");

    }

  return EVENT_OK;
}


// --- InitRun, part of Fun4All inheritance
int VTX_event_plane_reco::InitRun(PHCompositeNode *topNode)
{

  int runnumber = 0;

  RunHeader *rh = findNode::getClass<RunHeader>(topNode, "RunHeader");
  if(!rh)
    {
      cout<<PHWHERE<<" ERROR::RunHeader not found"<<endl;
      return ABORTEVENT;
    }
  else
    {
      runnumber = rh->get_RunNumber();
    }

  recoConsts *rc = recoConsts::instance();
  int icalibversion = 4002; // after Run5 and Field_ON

  rc->set_IntFlag("BBCCALIBVERSION", icalibversion);

  if(m_bbccalib){
    PHTimeStamp TimeStp = rc->get_TimeStamp();
    int BBCCALIBVERSION = rc->get_IntFlag("BBCCALIBVERSION");
    cout << "SvxRpSumXYReco::InitRun - restored constants are for " << TimeStp << endl;
    m_bbccalib->restore(TimeStp, BBCCALIBVERSION);
  }

  return EVENT_OK;
}


// --- ResetEvent, part of Fun4All inheritance, called after every event by Fun4All
int VTX_event_plane_reco::ResetEvent(PHCompositeNode *topNode)
{
  if (_verbosity > 1) cout << PHWHERE << "::ResetEvent() - entered." << endl;

  event         = -9999;
  centrality    = -9999;
  trigger_scaled = -9999;
  trigger_live   = -9999;
  bbc_qn        = -9999;
  bbc_qs        = -9999;
  bbc_z         = -9999;
  vtx_x         = -9999;
  vtx_y         = -9999;
  vtx_z         = -9999;
  bc_x          = -9999;
  bc_y          = -9999;
  svtx_z        = -9999;
  vtxposE_x     = -9999;
  vtxposE_x     = -9999;
  vtxposE_y     = -9999;
  vtxposE_z     = -9999;
  vtxposW_x     = -9999;
  vtxposW_y     = -9999;
  vtxposW_z     = -9999;
  fvtx_z        = -9999;
  nsegments     = -9999;
  ntracklets    = -9999;
  eventok       = -9999;

  for(int i = 0; i < 9; i++)
    {
      d_Qx[i] = -9999;
      d_Qy[i] = -9999;
      d_Qw[i] = 0;
    }

  for(int i = 0; i < 64; i++)
    {
      d_BBC_charge[i] = 0.0;

    }

  for(int i = 0; i <N_SVX_TRACK_MAX; i++)
    {
      trackID[i]   = -9999;
      charge[i]    = -9999;
      chisq[i]     = -9999;
      ndf[i]       = -9999;
      nhit0[i]     = -9999;
      nhit1[i]     = -9999;
      nhit2[i]     = -9999;
      nhit3[i]     = -9999;
      px[i]        = -9999;
      py[i]        = -9999;
      pz[i]        = -9999;
      dca[i]       = -9999;
      dca2d[i]     = -9999;
      segmentok[i] = -9999;
    }

  for(int i = 0; i < N_FVTX_CLUSTER_MAX; i++)
    {
      d_FVTX_x[i] = -9999;
      d_FVTX_y[i] = -9999;
      d_FVTX_z[i] = -9999;
    }

  if(_write_fvtx)
    for(int i = 0; i < 75; i++)
      {
        farm[i]      = -9999;
        fphi[i]      = -9999;
        feta[i]      = -9999;
        fchisq[i]    = -9999;
        fnhits[i]    = -9999;
        fDCA_X[i]    = -9999;
        fDCA_Y[i]    = -9999;
      }

  return EVENT_OK;

}



// --- process_event, inherited from Fun4All, does the main part of the analysis on an event by event basis
int VTX_event_plane_reco::process_event(PHCompositeNode *topNode)
{

  if (_verbosity > 1) cout << PHWHERE << "::process_event() - entered on event #" << _ievent << endl;
  else if ((_ievent)%10000 == 0) cout << PHWHERE << "::process_event() - event #" << _ievent << endl;

  // advance event counter
  _ievent++;

  event = _ievent;



  //---------------------------------
  // Run Number Selection Capability
  //---------------------------------

  int runnumber = 0;

  RunHeader *rh = findNode::getClass<RunHeader>(topNode, "RunHeader");
  if(!rh)
    {
      cout<<PHWHERE<<" ERROR::RunHeader not found"<<endl;
      return ABORTEVENT;
    }
  else
    {
      runnumber = rh->get_RunNumber();
    }

  if(_use_runlist)
    {
      if(is_run_in_list(runnumber)==false)
        {
          if(_verbosity > 1){
            cout<<endl<<"ABORTING RUN, Run number: "<<runnumber<<" is not in the run list"<<endl<<endl;
          }
          return ABORTEVENT;
        }
    }


  //-------------------------------
  // Grab info off of the node tree
  //-------------------------------

  //------------------PreviousEvent--------------------------------
  PreviousEvent *d_peve = NULL;
  d_peve    = findNode::getClass<PreviousEvent>(topNode, "PreviousEvent"); // for Tick cut
  if (d_peve == NULL)
    {

      std::cout << "SvxQAEventSelection::GetNodes -"
                << " No PreviousEvent object !" << std::endl;

      return ABORTEVENT;
    }

  TrigLvl1 *triggers = findNode::getClass<TrigLvl1>(topNode, "TrigLvl1");
  if(!triggers)
    {
      cout<<PHWHERE<<" ERROR::TrigLvl1 not found"<<endl;
      return ABORTEVENT;
    }
  PHGlobal *global = findNode::getClass<PHGlobal>(topNode, "PHGlobal");
  if(!global)
    {
      cout<<PHWHERE<<" ERROR::PHGlobal not found"<<endl;
      return ABORTEVENT;
    }
  EventHeader *evthead = findNode::getClass<EventHeader>(topNode, "EventHeader");
  if(!evthead)
    {
      cout<<PHWHERE<<" ERROR::EventHeader not found"<<endl;
      return ABORTEVENT;
    }
  VtxOut *vertexes = findNode::getClass<VtxOut>(topNode, "VtxOut");
  if(!vertexes)
    {
      cout<<PHWHERE<<" ERROR::VtxOut not found"<<endl;
      return ABORTEVENT;
    }
  SvxSegmentList *segments = findNode::getClass<SvxSegmentList>(topNode, "SvxSegmentList");
  if(!segments && _write_vtx)
    {
      cout<<PHWHERE<<" ERROR::SvxSegmentList not found"<<endl;
      return ABORTEVENT;
    }

  TFvtxCompactTrkMap* trkfvtx_map = NULL;
  if(_write_fvtx)
    {
      trkfvtx_map = findNode::getClass<TFvtxCompactTrkMap>(topNode,"TFvtxCompactTrkMap");
      if(!trkfvtx_map){
        cerr << PHWHERE << " No TFvtxCompactTrkMap object !" << endl;
        return ABORTEVENT;
      }
    }

  TFvtxCompactCoordMap* fvtx_coord_map = findNode::getClass<TFvtxCompactCoordMap>(topNode,"TFvtxCompactCoordMap");
  if(!fvtx_coord_map && _write_fvtx_clusters)
    {
      cout<<PHWHERE<<" ERROR::TFvtxCompactCoordMap not found"<<endl;
      return ABORTEVENT;
    }




  SvxClusterList *d_svxcls = NULL;
  d_svxcls = findNode::getClass<SvxClusterList>(topNode, "SvxClusterList");
  if (d_svxcls == NULL && _write_clusters)
    {
      std::cerr << "SvxClusterList node not found.Register SvxReco module" << std::endl;
      return ABORTEVENT;
    }

  RpSumXYObject* d_rp = findNode::getClass<RpSumXYObject>(topNode, "RpSumXYObject");

  if(!d_rp){
    if ( _verbosity > 4 ) cout<< PHWHERE << "Could not find the RPSumXYObject"<< endl;
    //return ABORTEVENT;
  }

  BbcRaw *bbcraw=findNode::getClass<BbcRaw>(topNode,"BbcRaw");
  if(!bbcraw){
    cout << "could not find Bbcraw!" << endl;
    if(_write_bbc)
      return ABORTEVENT;
  }


  int pticks[3] = {0};
  for ( int i = 0; i < 3; i++ ) pticks[i] = (d_peve != NULL) ? d_peve->get_clockticks(i) : -999;

  bool tick =  ( ( 20 < pticks[0] && pticks[0] < 120) ||
                 (700 < pticks[1] && pticks[1] < 780) );


  bool failed_tick_cut = false;

  if (tick)
    {
      if(_verbosity > 1)
        {
          std::cout << "SvxQAEventSelection::EventSelection - "
                    << "Failed tick cut. Not a good event!" << std::endl;
          std::cout << "                                      pticks: "
                    << pticks[0] << " " << pticks[1] << " " << pticks[2] << std::endl;
        }
      //return ABORTEVENT;
      failed_tick_cut = true;// skip SVX Segments if true
    }


  //---------------------------------------------------------//
  //
  //         Reading in Global Event Information into Tree
  //
  //---------------------------------------------------------//

  PHPoint precise_vertex1 = vertexes->get_Vertex("SVX_PRECISE");
  vtx_x = precise_vertex1.getX();
  vtx_y = precise_vertex1.getY();
  vtx_z = precise_vertex1.getZ();
  if(vtx_z!=vtx_z) vtx_z = -9999;//NAN check

  PHPoint svx_fast = vertexes->get_Vertex("SVX");//seed vertex
  bc_x = svx_fast.getX();//these are actually the beam center
  bc_y = svx_fast.getY();
  svtx_z = svx_fast.getZ();

  // phglobal fields
  centrality  = global->getCentrality();
  bbc_qn      = global->getBbcChargeN();
  bbc_qs      = global->getBbcChargeS();
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



  // unsigned int passes_trigger = trigger_scaled & accepted_triggers;
  // if ( passes_trigger == 0 )
  //   {
  //     if ( _verbosity > 0 ) cout << "trigger rejected" << endl;
  //     return ABORTEVENT;
  //   }
  // else if ( _verbosity > 0 ) cout << "trigger accepted" << endl;




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
  // float zvtx = -9999;
  // if ( runnumber >= 454744 && runnumber <= 456283 ) zvtx = bbc_z;
  // if ( runnumber >= 456652 && runnumber <= 458167 ) zvtx = FVTX_Z;
  // if ( fabs(zvtx) > _z_vertex_range )
  //   {
  //     if ( _verbosity > 0 ) cout << "rejecting event because of bad vertex " << zvtx << " cm" << endl;
  //     return ABORTEVENT;
  //   }
  // else if ( _verbosity > 0 ) cout << "event accepted vertex is " << zvtx << " cm" << endl;

  if ( _verbosity > 1 ) cout << "FVTX vertex points: " << FVTX_X << " " << FVTX_Y << " " << FVTX_Z << endl;
  // cout << "FVTX_Z is " << FVTX_Z << endl;




  int ibbcz_bin = (bbc_z+30.0)/10;//for fvtx eta cuts

  //get the East vertex (if available)
  PHPoint vtxposE;
  vtxposE = vertexes->get_Vertex("SVX_PRECISEE");
  vtxposE_x = vtxposE.getX();
  vtxposE_y = vtxposE.getY();
  vtxposE_z = vtxposE.getZ();

  //get the West vertex (if available)
  PHPoint vtxposW;
  vtxposW = vertexes->get_Vertex("SVX_PRECISEW");
  vtxposW_x = vtxposW.getX();
  vtxposW_y = vtxposW.getY();
  vtxposW_z = vtxposW.getZ();

  //---------------------------------------------------------//
  //  Finished Reading in Global Event Information into Tree
  //---------------------------------------------------------//


  //---------------------------------------------------------//
  //
  //           Writing out VTX Clusters in nTuple
  //
  //---------------------------------------------------------//


  if(_write_clusters)
    {
      PHPoint default_vertex1 = vertexes->get_Vertex();
      float default_z = default_vertex1.getZ();//

      for (int iclus = 0; iclus < d_svxcls->get_nClusters(); iclus++)
        {
          SvxCluster *clus = d_svxcls->get_Cluster(iclus);
          if (!clus) continue;
          float layer = clus->get_layer();
          float adc1 = clus->get_adc(0);
          float adc2 = clus->get_adc(1);
          float x = clus->get_xyz_global(0);
          float y = clus->get_xyz_global(1);
          float z = clus->get_xyz_global(2);
          float lx = clus->get_xyz_local(0);
          float ly = clus->get_xyz_local(1);
          float lz = clus->get_xyz_local(2);
          float section = clus->get_svxSection();
          float ladder = clus->get_ladder();
          float sensor = clus->get_sensor();
          float size = clus->get_size();
          float sizex = clus->get_xz_size(0);
          float sizez = clus->get_xz_size(1);

          float clus_data[17] =
            {
              _ievent,
              layer,
              x,
              y,
              z,
              adc1,
              adc2,
              section,
              ladder,
              sensor,
              lx,
              ly,
              lz,
              size,
              sizex,
              sizez,
              default_z
            };

          _ntp_cluster->Fill(clus_data);

        }
    }
  //---------------------------------------------------------//
  //         Finished writing out VTX Clusters in nTuple
  //---------------------------------------------------------//

  //if(_write_fvtx && !trkfvtx_map) return ABORTEVENT;

  //if( failed_tick_cut && _write_vtx) return ABORTEVENT;

  for(int i=0; i<9; i++){//fvtx bbc cnt smd vtx 3+3+1+3+%d
    d_Qx[i] = -9999;
    d_Qy[i] = -9999;
    d_Qw[i] = 0.0;
  }

  float sumxy_fvtx[3][8][3];
  for(int ihar=0; ihar<3; ihar++){//fvtx bbc cnt smd vtx 3+3+1+3+%d
    for(int idet=0; idet<8; idet++){
      for(int ixy=0; ixy<3; ixy++){
        sumxy_fvtx[ihar][idet][ixy] = 0.0;
      }
    }
  }

  eventok = is_event_ok(topNode);//only z vertex cut and centrality cut

  if (!eventok)
    {
      return ABORTEVENT;
    }

  //  tmp_evt++;//to keep track of how many events pass event cuts

  //---------------------------------------------------------//
  //
  //                  Retrieving Event Planes
  //
  //---------------------------------------------------------//

  for(int idet=62; idet<68; idet++){
    for(int ihar=0; ihar<3; ihar++){
      int idcode = -9999;

      if(idet<65) idcode = RP::calcIdCode(RP::ID_MPC,idet-62, ihar);
      else idcode = RP::calcIdCode(RP::ID_BBC,idet-65, ihar);

      if(idcode>=0 && d_rp){
        RpSnglSumXY *s_rp = d_rp->getRpSumXY(idcode);

        if(idet==62){//MPC
          d_Qx[ihar*3+idet-62]=s_rp->QVector(0);
          d_Qy[ihar*3+idet-62]=s_rp->QVector(1);
          d_Qw[ihar*3+idet-62]=s_rp->Weight();
        }
        if(idet==65){//BBC
          d_Qx[ihar*3+idet-65+2]=s_rp->QVector(0);
          d_Qy[ihar*3+idet-65+2]=s_rp->QVector(1);
          d_Qw[ihar*3+idet-65+2]=s_rp->Weight();
        }
      }
    }
  }
  //FVTX
  //cout<<RP::ID_FVT<<endl;
  //south

  for(int idet=0; idet<20; idet++){
    for(int ihar=0; ihar<3; ihar++){
      int idcode = RP::calcIdCode(RP::ID_FVT,idet, ihar);

      if(idcode>=0 && d_rp){
        RpSnglSumXY *s_rp = d_rp->getRpSumXY(idcode);
        int ifvtx = idet%5;
        if(ifvtx<4){//-1.0-3.0
          sumxy_fvtx[ihar][0][0]+=s_rp->QVector(0);//compiling the 20 FVTX event planes
          sumxy_fvtx[ihar][0][1]+=s_rp->QVector(1);
          sumxy_fvtx[ihar][0][2]+=s_rp->Weight();
        }
      }
    }
  }

  for(int ihar=0; ihar<3; ihar++){
    //for(int idet=0; idet<8; idet++){
    d_Qx[ihar*3+1]=sumxy_fvtx[ihar][0][0];
    d_Qy[ihar*3+1]=sumxy_fvtx[ihar][0][1];
    d_Qw[ihar*3+1]=sumxy_fvtx[ihar][0][2];
    //}
  }


  //---------------------------------------------------------//
  //             Finished retrieving Event Planes
  //---------------------------------------------------------//

  //---------------------------------------------------------//
  //
  //            Writing out BBC and FVTX Raw Objects
  //
  //---------------------------------------------------------//

  if(_write_bbc)
    {
      if(bbcraw)
        {
          for(int ipmt=0; ipmt<64; ipmt++)
            {

              short iadc    = bbcraw->get_Adc(ipmt);
              short itdc    = bbcraw->get_Tdc0(ipmt);
              //float tdc1=bbcraw->get_Tdc1(ipmt);
              float itime0  = m_bbccalib->getHitTime0(ipmt,itdc,iadc);
              float icharge = m_bbccalib->getCharge(ipmt,iadc);
              if ( ( itime0 <=0 || icharge <=0 ) && verbosity > 8 )
                {
                  // --- this happens A LOT
                  cout << "EITHER SIDE!!!" << endl;
                  cout << "For event number " << _ievent-1
                       << " and BBC tube number " << ipmt
                       << " the time is " << itime0
                       << " and the charge is " << icharge
                       << endl;
                }
              if ( icharge <= 0 && itime0 >=0 )
                {
                  // --- this never happens
                  cout << "Holy shit!" << endl;
                  cout << "For event number " << _ievent-1
                       << " and BBC tube number " << ipmt
                       << " the time is " << itime0
                       << " and the charge is " << icharge
                       << endl;
                }
              if ( itime0 <= 0 ) icharge -= 10000.0;
              //float ibbc_x  = m_bbcgeo->getX(ipmt);
              //float ibbc_y  = m_bbcgeo->getY(ipmt);
              float ibbc_z  = m_bbcgeo->getZ(ipmt);
              //cout<<"d_pmt_x["<<ipmt<<"] = "<<ibbc_x<<";"<<endl;
              //cout<<"d_pmt_y["<<ipmt<<"] = "<<ibbc_y<<";"<<endl;
              //cout<<"d_pmt_z["<<ipmt<<"] = "<<ibbc_z<<";"<<endl;
              // --- i think this cut is unnecessary...
              if(ibbc_z > 0) continue;//only select on south bbc

              if ( ( itime0 <=0 || icharge <=0 ) && verbosity > 8 )
                {
                  // --- this happens A LOT
                  cout << "SOUTH SIDE!!!" << endl;
                  cout << "For event number " << _ievent-1
                       << " and BBC tube number " << ipmt
                       << " the time is " << itime0
                       << " and the charge is " << icharge
                       << endl;
                }

              d_BBC_charge[ipmt] = icharge;

            }//end of ipmt loop
        }
    }


  int nfvtxs_raw_clus = 0;
  if ( fvtx_coord_map )
    {
      TFvtxCompactCoordMap::iterator _iter( fvtx_coord_map->range() );
      while( TFvtxCompactCoordMap::const_pointer fvtx_ptr = _iter.next() )
        {
          TFvtxCompactCoord* fvtxcoord = fvtx_ptr->get();
          PHPoint fvtx_coord_point = fvtxcoord->get_coord_midpoint();
          //int iarm = fvtxcoord->get_arm();
          float fvtx_x = fvtx_coord_point.getX();
          float fvtx_y = fvtx_coord_point.getY();
          float fvtx_z = fvtx_coord_point.getZ();
          //float fvtx_r = sqrt(pow(fvtx_x,2.0)+pow(fvtx_y,2.0));
          if( (fabs(fvtx_x)>999) ||(fabs(fvtx_y)>999) || (fabs(fvtx_z)>999)) continue;
          //float fvtx_the = atan2(fvtx_r,fvtx_z-vtx_z);
          //float fvtx_phi = atan2(fvtx_y,fvtx_x);
          //float fvtx_eta = -log(tan(0.5*fvtx_the));
          //if(fvtx_z < 0)
          //{
          if(nfvtxs_raw_clus >= N_FVTX_CLUSTER_MAX)
            {
              cout<<"butting against the max fvtx cluster size " << nfvtxs_raw_clus << "/" << N_FVTX_CLUSTER_MAX << ", breaking"<<endl;
              break;
            }
          d_FVTX_x[nfvtxs_raw_clus] = fvtx_x;
          d_FVTX_y[nfvtxs_raw_clus] = fvtx_y;
          d_FVTX_z[nfvtxs_raw_clus] = fvtx_z;
          nfvtxs_raw_clus++;
          //cout<<"fvtx_eta: "<<fvtx_eta<<endl;
          //} // south

        } // while loop over iterator
    } // check on fvtx_coord_map

  d_nFVTX_clus = nfvtxs_raw_clus;
  //cout<<"nfvtxs_raw_clus: "<<nfvtxs_raw_clus <<endl;

  //---------------------------------------------------------//
  //
  //                 Get FVTX Tracks
  //
  //---------------------------------------------------------//

  int ntr = -1;

  if (trkfvtx_map && _write_fvtx){
    TFvtxCompactTrkMap::const_iterator trk_iter = trkfvtx_map->range();
    while( TFvtxCompactTrkMap::const_pointer trk_ptr = trk_iter.next() ){

      TFvtxCompactTrk* fvtx_trk = trk_ptr->get();

      float the = fvtx_trk->get_fvtx_theta();
      float eta = fvtx_trk->get_fvtx_eta();
      float phi = fvtx_trk->get_fvtx_phi();
      int   arm = (int)fvtx_trk->get_arm();
      float fvtx_x      = fvtx_trk->get_fvtx_vtx().getX();
      float fvtx_y      = fvtx_trk->get_fvtx_vtx().getY();
      float fvtx_z      = fvtx_trk->get_fvtx_vtx().getZ();
      int   nfhits      = (int)fvtx_trk->get_nhits();

      float DCA_x      = fvtx_x + tan(the)*cos(phi)*(bbc_z - fvtx_z);
      float DCA_y      = fvtx_y + tan(the)*sin(phi)*(bbc_z - fvtx_z);
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
      if(the==0) continue;

      if(nfhits<3) continue;
      if(!pass_eta_cut(eta,ibbcz_bin)) continue;
      if(fvtx_trk->get_chi2_ndf() > 5) continue;
      if(fabs(DCA_x) > 2.0 || fabs(DCA_y) > 2.0) continue;
      ntr++;


      //float DCA_R      = sqrt((DCA_x*DCA_x) + (DCA_y*DCA_y));

      if(ntr < 75)
        {
          feta[ntr]   = eta;
          fphi[ntr]   = phi;
          fchisq[ntr] = fvtx_trk->get_chi2_ndf();
          farm[ntr]   = arm;
          fnhits[ntr] = nfhits;
          fDCA_X[ntr] = DCA_x;
          fDCA_Y[ntr] = DCA_y;
        }
      else
        {
          cout<<"butting up against the boundary of fvtx tracks"<<endl;
          break;
        }
      //short_chi2 = fvtx_trk->get_short_chi2_ndf();
      //chi2       = fvtx_trk->get_chi2_ndf();

    }
  }

  ntracklets = ntr;

  //---------------------------------------------------------//
  //                 finished Get FVTX Tracks
  //---------------------------------------------------------//

  //---------------------------------------------------------//
  //
  //                 Get VTX Tracks
  //
  //---------------------------------------------------------//

  if(_write_vtx)
    nsegments = segments->get_nSegments();
  //  cout<<"nsegments in this event: "<<nsegments<<endl;

  int igoodseg = 0;

  if(segments && !failed_tick_cut && _write_vtx )
    {
      //vector<int> indexes;
      //for(int isegment = 0; isegment < nsegments; isegment++)
      //{
      //indexes.push_back(isegment);
      //}
      //random_shuffle( indexes.begin(), indexes.end() );

      for(int isegment = 0; isegment < nsegments; isegment++)
        {

          SvxSegment *segment = segments->get_segment(isegment);

          if(!is_segment_ok(segment)) continue;
          if(igoodseg >= N_SVX_TRACK_MAX)
            {
              cout<<"bumping up againt track limit"<<endl;
              break;
            }
          charge[igoodseg] = 0;
          if(segment->IsPositive())
            {
              charge[igoodseg] = +1.0;
            }
          else
            {
              charge[igoodseg] = -1.0;
            }

          float ipx = segment->get3MomentumAtPrimaryVertex(0);
          float ipy = segment->get3MomentumAtPrimaryVertex(1);
          float ipz = segment->get3MomentumAtPrimaryVertex(2);

          trackID[igoodseg]   = segment->getSegmentID();
          chisq[igoodseg]     = calc_chisq_fromquality(segment->getSegmentQuality(),segment->getSegmentScore());
          ndf[igoodseg]       = segment->getNDF();
          nhit0[igoodseg]     = segment->getNhits(0);
          nhit1[igoodseg]     = segment->getNhits(1);
          nhit2[igoodseg]     = segment->getNhits(2);
          nhit3[igoodseg]     = segment->getNhits(3);
          dca[igoodseg]       = segment->getDCA();
          dca2d[igoodseg]     = segment->getDCA2D();
          px[igoodseg]       = ipx;
          py[igoodseg]       = ipy;
          pz[igoodseg]       = ipz;
          segmentok[igoodseg] = is_segment_ok(segment);

          igoodseg++;
	  //if(!is_segment_ok(segment)) continue;

        }
    }

  nsegments = igoodseg;

  //---------------------------------------------------------//
  //                 finished Get VTX Tracks
  //---------------------------------------------------------//

  //if(nsegments==0) return EVENT_OK; // BAD CUT

  d_ntrk = 0;
  PHCentralTrack *ctrk = findNode::getClass<PHCentralTrack>(topNode,"PHCentralTrack");
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

          float mom         = strk->get_mom();
          float zed         = strk->get_zed();
          int quality       = strk->get_quality();
          if ( mom < 0.0 || mom > 50.0 ) continue;
          if ( fabs(zed) < 3.0 || fabs(zed) > 70.0 ) continue;
          if ( quality != 63 && quality != 31 ) continue;

          float px = strk->get_px();
          float py = strk->get_py();
          float pz = strk->get_pz();

          int arm = 0;
          if ( px > 0 ) arm = 1;

          int charge = strk->get_charge();

          float pc3dphi = strk->get_pc3dphi();
          float pc3dz = strk->get_pc3dz();

          if ( pc3dphi < -9990 || pc3dz < -9990 ) continue;

          float pc3sdphi = calcsdphi(pc3dphi, arm, charge, mom);
          float pc3sdz = calcsdz(pc3dz, arm, charge, mom);

          if ( fabs(pc3sdphi) > 3.0 || fabs(pc3sdz) > 3.0 ) continue;

          d_cntpx[counter] = px;
          d_cntpy[counter] = py;
          d_cntpz[counter] = pz;
          d_cntpc3sdz[counter] = pc3sdz;
          d_cntpc3sdphi[counter] = pc3sdphi;

          ++counter;

        } // loop over tracks

      d_ntrk = counter;

    } // check on track pointer

  if ( _create_ttree ) _ntp_event->Fill();

  ResetEvent(topNode);

  if ( _verbosity > 0 ) cout << "sucessfully processed this event" << endl;

  tmp_evt++;//to keep track of how many events pass event cuts
  return EVENT_OK;

} // end of process_event




int VTX_event_plane_reco::End(PHCompositeNode *topNode)
{
  if (_verbosity > 1) cout << PHWHERE << "::End() - entered." << endl;

  cout<<"total events: "<<_ievent<<" fraction passing vtx cut: "<<tmp_evt*1.0/_ievent<<endl;

  _output_file->cd();

  if(_create_ttree)
    {
      if(!_write_clusters)
        _ntp_event->Write();
      else if(_write_clusters)
        _ntp_cluster->Write();
    }
  _output_file->Close();
  delete _output_file;

  //delete  m_bbccalib;
  //delete  m_bbcgeo;

  //cout << PHWHERE << "::End() - ntuple output written to file: " << _output_filename << endl;


  return EVENT_OK;
}


bool VTX_event_plane_reco::is_event_ok(PHCompositeNode *topNode)
{

  if (_verbosity > 1) cout << PHWHERE << "::is_event_ok() - entered." << endl;

  // grab from the node tree
  //TrigLvl1 *triggers = findNode::getClass<TrigLvl1>(topNode, "TrigLvl1");
  //if(!triggers)
  //{
  //cout<<PHWHERE<<" ERROR::TrigLvl1 not found"<<endl;
  //return ABORTEVENT;
  //}
  PHGlobal *global = findNode::getClass<PHGlobal>(topNode, "PHGlobal");
  if(!global)
    {
      cout<<PHWHERE<<" ERROR::PHGlobal not found"<<endl;
      return ABORTEVENT;
    }
  VtxOut *vertexes = findNode::getClass<VtxOut>(topNode, "VtxOut");
  if(!vertexes)
    {
      cout<<PHWHERE<<" ERROR::VtxOut not found"<<endl;
      return ABORTEVENT;
    }

  // --- why not use the class variable?  cut now applied below with vtxout object
  // --- bbc_z
  // double bbcz = global->getBbcZVertex();
  // if(fabs(bbcz) > 30)
  //   {
  //     if(_verbosity > 0)
  //       cout<<"event rejected because bbc z vertex outside of 30 cm"<<endl;
  //     return false;
  //   }

  return true;

  // bail on bad centrality
  float tmpcentrality = global->getCentrality();
  if((tmpcentrality < 0.0)||(tmpcentrality > 100.0))
    {
      if(_verbosity > 0)
        cout<<"event rejected because of outside of sensible centrality range"<<endl;
      return false;
    }
  //if(centrality < 20.0) return false; // peripheral cut
  if(tmpcentrality > 5.0)
    {
      if(_verbosity > 0)
        cout<<"event rejected because outside of 5%% centrality"<<endl;
      return false;
    }
  return true;

} // is_event_ok




bool VTX_event_plane_reco::is_segment_ok(SvxSegment *segment)
{
  // single particle cuts go here...
  if (_verbosity > 1) cout << PHWHERE << "::is_segment_ok() - entered." << endl;

  // reject tracks with too few hits...
  unsigned int nhits = segment->getNhits(0)
    + segment->getNhits(1)
    + segment->getNhits(2)
    + segment->getNhits(3);
  //cout<<"nhit0: "<<segment->getNhits(0)<<" nhit1: "<<segment->getNhits(1)<<" nhit2: "<<segment->getNhits(2)<<" nhit3: "<<segment->getNhits(3)<<endl;
  if(nhits < 4)
    {
      if(_verbosity > 1)
        cout<<"track rejected because not enough hits"<<endl;
      return false;
    }
  //if(segment->getNhits(0)+segment->getNhits(1)!=2)
  //{
  //if(_verbosity > 1)
  //cout<<"track rejected because no hit in b0 and b1"<<endl;
  //return false;
  //}
  // reject tracks with off-vertex dca...
  if(fabs(segment->getDCA2D()) > 0.03 || fabs(segment->getDCA())> 0.5)
    {
      if(_verbosity > 1)
        cout<<"track rejected because of DCA cut"<<endl;
      return false;
    }
  // reject tracks with chisq/ndf<%d...
  //int ndf = 2*nhits - 5;
  double px = segment->get3MomentumAtPrimaryVertex(0);
  double py = segment->get3MomentumAtPrimaryVertex(1);
  double pt = sqrt(px*px+py*py);
  float chisqndf = calc_chisq_fromquality(segment->getSegmentQuality(),segment->getSegmentScore());
  if(!pass_chisq_cut(chisqndf,pt,nhits))
    {
      if(_verbosity > 1)
        cout<<"track rejected because chisq cut"<<endl;
      return false;
    }

  return true;
}


bool VTX_event_plane_reco::pass_chisq_cut(double chisq, double pt, int nhits)
{
  if(nhits == 3)
    {
      if((pt < 1.0)&&(chisq<3)) return true;
      if((pt > 1.0)&&(chisq<2)) return true;
    }
  else if(nhits == 4)
    {
      if(chisq<2) return true;
    }

  return false;
}

bool VTX_event_plane_reco::is_run_in_list(int runnumber)
{
  ifstream runlist;
  runlist.open(_runlist_filename.c_str());
  int run_num;
  while(1){
    if(!runlist.good()) break;
    runlist >> run_num;
    if(run_num==runnumber)
      return true;
  }
  return false;
}

//these vaules were obtained from ana note 1162
bool VTX_event_plane_reco::pass_eta_cut(float eta, int bbcz_bin)
{
  if(bbcz_bin==0)
    {
      if(eta > 1.9 && eta < 3.2)
        return true;
    }
  if(bbcz_bin==1)
    {
      if(eta > -1.7 && eta < -0.8)
        return true;
      if(eta > 1.7 && eta < 3.2)
        return true;
    }
  if(bbcz_bin==2)
    {
      if(eta > -2.1 && eta < -0.7)
        return true;
      if(eta > 1.5 && eta < 3.1)
        return true;
    }
  if(bbcz_bin==3)
    {
      if(eta > -2.4 && eta < -0.9)
        return true;
      if(eta > 1.3 && eta < 3.0)
        return true;
    }
  if(bbcz_bin==4)
    {
      if(eta > -2.6 && eta < -1.0)
        return true;
      if(eta > 1.2 && eta < 2.8)
        return true;
    }
  if(bbcz_bin==5)
    {
      if(eta > -2.8 && eta < -1.2)
        return true;
      if(eta > 1.0 && eta < 2.6)
        return true;
    }
  if(bbcz_bin==6)
    {
      if(eta > -3.0 && eta < -1.3)
        return true;
      if(eta > 0.9 && eta < 2.4)
        return true;
    }
  if(bbcz_bin==7)
    {
      if(eta > -3.1 && eta < -1.4)
        return true;
      if(eta > 0.7 && eta < 2.1)
        return true;
    }
  if(bbcz_bin==8)
    {
      if(eta > -3.2 && eta < -1.5)
        return true;
      if(eta > 0.7 && eta < 1.6)
        return true;
    }
  if(bbcz_bin==9)
    {
      if(eta > -3.3 && eta < -1.8)
        return true;
    }

  return false;
}


float VTX_event_plane_reco::calc_chisq_fromquality(float quality, float score)
{
  float chisq = quality - score/100.;
  chisq = 1/chisq-2.0;
  if (chisq < 0)
    {
      //if(_verbosity > 0)
      //cout << "WARNING!! calc_chisq_fromquality(" << quality << "," << score <<") gives chisq/ndf=" << chisq << endl;
      return 99999;
    }
  return chisq;
}
