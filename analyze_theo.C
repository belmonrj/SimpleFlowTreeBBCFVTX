void analyze_theo(int runnumber, int segmentnumber = 0)
{


  //----------
  // Libraries
  //----------

  //gSystem->Load("libfun4all.so");
  //gSystem->Load("libfun4allfuncs.so");
  gSystem->Load("libCNT.so");
  gSystem->Load("libcompactCNT.so");
  //gSystem->Load("librecal.so");

  gSystem->Load("libfvtx_subsysreco.so");
  gSystem->Load("libmutoo_subsysreco");
  gSystem->Load("libfun4all.so");
  gSystem->Load("libfun4allfuncs.so");
  gSystem->Load("libfun4allfuncs_muons");
  gSystem->Load("libMWGOO");
  gSystem->Load("libmutrg");
  gSystem->Load("librpc_subsysreco");
  gSystem->Load("librpc_muotrackreco");
  gSystem->Load("libcompactCNT.so");
  gSystem->Load("librecal.so");
  gSystem->Load("libBbcMultipleVtx.so" );


  gSystem->Load("libVTX_event_plane.so");

  //gSystem->Load("libReadnDST.so");
  gROOT->SetBatch(true);

  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity(0);
  //recoConsts *rc = recoConsts::instance();

  MasterRecalibratorManager *mr = new MasterRecalibratorManager();
  mr->FillHistos(0);
  se->registerSubsystem(mr);

  FvtxReadbackDST *fvtxoo = new FvtxReadbackDST();
  fvtxoo->Verbosity( 0 );
  se->registerSubsystem( fvtxoo );

  //----------------------
  // VTX_event_plane_reco
  //----------------------
  int z_vertex_range = 10;
  VTX_event_plane_reco *vtx_ep = new VTX_event_plane_reco();
  vtx_ep->set_z_vertex_range(z_vertex_range);         //set z_vertex range for accepted tracks
  vtx_ep->set_use_runlist(false);                     //use the runlist to select for good runs
  vtx_ep->set_create_ttree(true);
  vtx_ep->set_write_clusters(false); // svx clusters not available
  vtx_ep->set_write_bbc(true);
  vtx_ep->set_write_fvtx_clusters(true);
  vtx_ep->set_write_fvtx(true);
  vtx_ep->set_trimmed_tree(true);
  vtx_ep->set_write_vtx(false); // svx tracks not available
  //vtx_ep->set_runlist_file(string);
  vtx_ep->set_output_filename(Form("out/tree_%010d_%04d.root",runnumber,segmentnumber));//set the ntuple output file
  vtx_ep->Verbosity(0);
  se->registerSubsystem(vtx_ep);

  Fun4AllInputManager *in1 = new Fun4AllDstInputManager("DSTin1", "DST");
  in1->Verbosity(0);
  se->registerInputManager(in1);

  cout << "Analysis started " << endl;
  gSystem->Exec("date");
  if ( runnumber >= 454000 & runnumber < 455000 ) se->fileopen("DSTin1", Form("/phenix/prod/online_production/run16_online_ca/run_0000454000_0000455000/CNT/CNT_run16_online_ca-%010d-%04d.root",runnumber,segmentnumber));
  if ( runnumber >= 455000 & runnumber < 456000 ) se->fileopen("DSTin1", Form("/phenix/prod/online_production/run16_online_ca/run_0000455000_0000456000/CNT/CNT_run16_online_ca-%010d-%04d.root",runnumber,segmentnumber));
  if ( runnumber >= 456000 & runnumber < 457000 ) se->fileopen("DSTin1", Form("/phenix/prod/online_production/run16_online_ca/run_0000456000_0000457000/CNT/CNT_run16_online_ca-%010d-%04d.root",runnumber,segmentnumber));
  if ( runnumber >= 457000 & runnumber < 458000 ) se->fileopen("DSTin1", Form("/phenix/prod/online_production/run16_online_ca/run_0000457000_0000458000/CNT/CNT_run16_online_ca-%010d-%04d.root",runnumber,segmentnumber));
  if ( runnumber >= 458000 & runnumber < 459000 ) se->fileopen("DSTin1", Form("/phenix/prod/online_production/run16_online_ca/run_0000458000_0000459000/CNT/CNT_run16_online_ca-%010d-%04d.root",runnumber,segmentnumber));


  Fun4AllInputManager *in2 = new Fun4AllDstInputManager("DSTin2", "DST");
  in2->Verbosity(0);
  se->registerInputManager(in2);

  gSystem->Exec("date");
  if ( runnumber >= 454000 & runnumber < 455000 ) se->fileopen("DSTin2", Form("/phenix/prod/online_production/run16_online_ca/run_0000454000_0000455000/DST_EVE/DST_EVE_run16_online_ca-%010d-%04d.root",runnumber,segmentnumber));
  if ( runnumber >= 455000 & runnumber < 456000 ) se->fileopen("DSTin2", Form("/phenix/prod/online_production/run16_online_ca/run_0000455000_0000456000/DST_EVE/DST_EVE_run16_online_ca-%010d-%04d.root",runnumber,segmentnumber));
  if ( runnumber >= 456000 & runnumber < 457000 ) se->fileopen("DSTin2", Form("/phenix/prod/online_production/run16_online_ca/run_0000456000_0000457000/DST_EVE/DST_EVE_run16_online_ca-%010d-%04d.root",runnumber,segmentnumber));
  if ( runnumber >= 457000 & runnumber < 458000 ) se->fileopen("DSTin2", Form("/phenix/prod/online_production/run16_online_ca/run_0000457000_0000458000/DST_EVE/DST_EVE_run16_online_ca-%010d-%04d.root",runnumber,segmentnumber));
  if ( runnumber >= 458000 & runnumber < 459000 ) se->fileopen("DSTin2", Form("/phenix/prod/online_production/run16_online_ca/run_0000458000_0000459000/DST_EVE/DST_EVE_run16_online_ca-%010d-%04d.root",runnumber,segmentnumber));

  se->run(0);

  se->End();
  //se->dumpHistos(outFile);
  //Fun4AllHistoManager *hm;
  //hm = se->getHistoManager("Run11AuAuSvxEventHists");
  //hm->dumpHistos(histfile);

  cout << "Analysis finished " << endl;

}
