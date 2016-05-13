void Run_VTX_event_plane(const char *outFile = "test_train_output.root")
{

  int z_vertex_range = 30;

  //----------
  // Libraries
  //----------
  //gSystem->Load("libfvtx_subsysreco.so");
  //gSystem->Load("libfun4allfuncs_muons");
  //gSystem->Load("libMWGOO");
  gSystem->Load("libVTX_event_plane");
  //gSystem->Load("libcompactCNT.so");
  //gSystem->Load("libTOAD");

  //--------------------------
  // TOAD absolute path finder
  //--------------------------

  //TOAD *toad_loader = new TOAD("theok");
  //string file_location = toad_loader->location("runlist.txt");
  //delete toad_loader;


  //---------------
  // Fun4All Server
  //---------------

  Fun4AllServer *se = Fun4AllServer::instance();


  //DeathToMemoryHogs *dt = new DeathToMemoryHogs();
  //se->registerSubsystem(dt);

  //to get the FVTX tracks
  FvtxReadbackDST *fvtxoo = new FvtxReadbackDST();
  fvtxoo->Verbosity( 0 );
  se->registerSubsystem( fvtxoo );

  //SvxCompactToDST* svxcomptodst = new SvxCompactToDST();
  //se->registerSubsystem(svxcomptodst);

  //SubsysReco *svxprimvtxfinder_west = new SvxPrimVertexFinder("SVXPRIMVTXFINDERW", 1);
  //se->registerSubsystem(svxprimvtxfinder_west);

  //SubsysReco *svxprimvtxfinder_east = new SvxPrimVertexFinder("SVXPRIMVTXFINDERE", 2);
  //se->registerSubsystem(svxprimvtxfinder_east);

  //FvtxPrimVertex* fvtxprimvtx = new FvtxPrimVertex();
  //fvtxprimvtx->set_fvtx_Rres(0.4); //crossing window (0.1 for CuAu, 0.5 for pp)
  //fvtxprimvtx->set_source( FvtxPrimVertex::Coordinate ); // or FvtxPrimVertex::Tracks
  //fvtxprimvtx->set_bbcz_window(2.0); // FVTX-BBC vertex matching window
  //se->registerSubsystem(fvtxprimvtx);

  //----------------------
  // VTX_event_plane_reco
  //----------------------
  VTX_event_plane_reco *vtx_ep = new VTX_event_plane_reco();
  vtx_ep->set_z_vertex_range(z_vertex_range);         //set z_vertex range for accepted tracks
  vtx_ep->set_use_runlist(false);                     //use the runlist to select for good runs
  vtx_ep->set_create_ttree(true);
  vtx_ep->set_write_clusters(false);
  vtx_ep->set_write_bbc(true);
  vtx_ep->set_write_fvtx_clusters(true);
  vtx_ep->set_write_fvtx(false);
  vtx_ep->set_trimmed_tree(true);
  vtx_ep->set_write_vtx(true);
  //vtx_ep->set_runlist_file(string);
  vtx_ep->set_output_filename(outFile);//set the ntuple output file
  vtx_ep->Verbosity(0);
  se->registerSubsystem(vtx_ep);

  //return;
}

void InputData(vector<string> &indata)
{
  indata.push_back("CNT");
  indata.push_back("DST_SVX");
  //indata.push_back("MWG");
  indata.push_back("dst_fvtx");
  indata.push_back("dst_eve");
  return;
}
