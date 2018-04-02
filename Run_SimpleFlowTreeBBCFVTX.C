void Run_SimpleFlowTreeBBCFVTX(
  const char *outFile = "test_train_output.root")
{

  //----------
  // Libraries
  //----------

  // gSystem->Load("libUltraLight");
  // gSystem->Load("libCabanaBoy");

  //-- required for Run15pAu200FvtxMBPro107
  // not for Run16dAu
  gSystem->Load("libfvtx_subsysreco.so");

  gSystem->Load("libdAuBES_utils.so");
  gSystem->Load("libDoubleInteractionUtil.so");
  gSystem->Load("libSimpleFlowTreeBBCFVTX.so");

  gSystem->ListLibraries();

  //------------------//
  // Fun4All Server --//
  //------------------//

  Fun4AllServer *se = Fun4AllServer::instance();

  // To get the FVTX tracks
  se->registerSubsystem(new FvtxReadbackDST());

  //---------------//
  // User Module --//
  //---------------//

  SimpleFlowTreeBBCFVTX *sflow = new SimpleFlowTreeBBCFVTX();
  sflow->set_create_ttree(true);
  sflow->set_write_bbc(true);
  sflow->set_write_cnt(true);
  sflow->set_write_fvtx_clusters(false);
  sflow->set_write_fvtx(true);
  sflow->set_use_runlist(false);
  //sflow->set_runlist_file(string);
  sflow->set_output_filename(outFile);//set the ntuple output file
  sflow->Verbosity(0);
  se->registerSubsystem(sflow);

}

void InputData(vector<string> &indata)
{
  indata.push_back("CNT");
  indata.push_back("DST_EVE");
  // indata.push_back("DST_FVTX"); // MWG file is here for run15pau
  indata.push_back("MWG"); // MWG file is here for run15pal
  return;
}
