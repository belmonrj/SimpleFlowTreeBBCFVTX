void cospsidiff_run(int, int);
void psidiff_run(int, int);

void beforeafterflat_run(int, int);


void event_plane_testing()
{

  // psidiff_run(455050,2);
  // cospsidiff_run(455050,2);

  beforeafterflat_run(454774,2); // 200
  return;
  beforeafterflat_run(455050,2); // 200
  beforeafterflat_run(455355,2); // 200
  beforeafterflat_run(456015,2); // 62
  beforeafterflat_run(456201,2); // 62
  beforeafterflat_run(457015,2); // 20
  beforeafterflat_run(457213,2); // 20
  beforeafterflat_run(458014,2); // 39
  beforeafterflat_run(458167,2); // 39

  beforeafterflat_run(455050,3); // 200
  beforeafterflat_run(455355,3); // 200
  beforeafterflat_run(456015,3); // 62
  beforeafterflat_run(456201,3); // 62
  beforeafterflat_run(457015,3); // 20
  beforeafterflat_run(457213,3); // 20
  beforeafterflat_run(458014,3); // 39
  beforeafterflat_run(458167,3); // 39

}


void beforeafterflat_run(int run, int harmonic)
{

  gStyle->SetOptTitle(0);

  TCanvas* c1 = new TCanvas("c1","");

  //TFile* file = TFile::Open(Form("input/hist_%d.root",run));
  TFile* file = TFile::Open(Form("input/files_200/hist_%d.root",run)); // uh oh
  if ( !file ) { cout << "no file" << endl; return; }
  TList* list = file->GetListOfKeys();
  if ( list->GetSize() < 1 ) { cout << "no keys" << endl; return; }

  TH1D* hbbcW = (TH1D*) ((TH2D*)file->Get(Form("psi_bf_0_%d_23",harmonic-1)))->ProjectionY();
  hbbcW->Draw();
  hbbcW->SetLineColor(kRed);
  hbbcW->GetXaxis()->SetTitle(Form("#Psi_{%d}, BBCS",harmonic));
  TH1D* hbbc = (TH1D*) ((TH2D*)file->Get(Form("psi_bf_0_%d_2",harmonic-1)))->ProjectionY();
  hbbc->Draw("same");
  hbbc->SetLineColor(kGreen+2);
  hbbc->GetXaxis()->SetTitle(Form("#Psi_{%d}, BBCS",harmonic));
  TH1D* hbbcM = (TH1D*) ((TH2D*)file->Get(Form("psi_mf_0_%d_2",harmonic-1)))->ProjectionY();
  hbbcM->Draw("same");
  hbbcM->SetLineColor(kBlue);
  TH1D* hbbcA = (TH1D*) ((TH2D*)file->Get(Form("psi_af_0_%d_2",harmonic-1)))->ProjectionY();
  hbbcA->Draw("same");
  hbbcA->SetLineColor(kBlack);
  TLegend* leg = new TLegend(0.68,0.68,0.88,0.88);
  leg->AddEntry(hbbcW,"Raw","l");
  leg->AddEntry(hbbc,"Gain correction","l");
  leg->AddEntry(hbbcM,"Recentering","l");
  leg->AddEntry(hbbcA,"Flattening","l");
  leg->SetTextSize(0.045);
  leg->Draw();
  c1->Print(Form("FigsEventPlane/after_psi%d_bbcs_%d.png",harmonic,run));
  c1->Print(Form("FigsEventPlane/after_psi%d_bbcs_%d.pdf",harmonic,run));
  delete leg;

  TH1D* hfvtxsW = (TH1D*) ((TH2D*)file->Get(Form("psi_bf_0_%d_13",harmonic-1)))->ProjectionY();
  hfvtxsW->Draw();
  hfvtxsW->SetLineColor(kRed);
  hfvtxsW->GetXaxis()->SetTitle(Form("#Psi_{%d}, FVTXS",harmonic));
  TH1D* hfvtxs = (TH1D*) ((TH2D*)file->Get(Form("psi_bf_0_%d_3",harmonic-1)))->ProjectionY();
  hfvtxs->Draw("same");
  hfvtxs->SetLineColor(kGreen+2);
  hfvtxs->GetXaxis()->SetTitle(Form("#Psi_{%d}, FVTXS",harmonic));
  TH1D* hfvtxsM = (TH1D*) ((TH2D*)file->Get(Form("psi_mf_0_%d_3",harmonic-1)))->ProjectionY();
  hfvtxsM->Draw("same");
  hfvtxsM->SetLineColor(kBlue);
  TH1D* hfvtxsA = (TH1D*) ((TH2D*)file->Get(Form("psi_af_0_%d_3",harmonic-1)))->ProjectionY();
  hfvtxsA->Draw("same");
  hfvtxsA->SetLineColor(kBlack);
  leg = new TLegend(0.68,0.68,0.88,0.88);
  leg->AddEntry(hfvtxsW,"Raw","l");
  leg->AddEntry(hfvtxs,"Weighting","l");
  leg->AddEntry(hfvtxsM,"Recentering","l");
  leg->AddEntry(hfvtxsA,"Flattening","l");
  leg->SetTextSize(0.045);
  leg->Draw();
  c1->Print(Form("FigsEventPlane/after_psi%d_fvtxs_%d.png",harmonic,run));
  c1->Print(Form("FigsEventPlane/after_psi%d_fvtxs_%d.pdf",harmonic,run));

  TH1D* hfvtxs0W = (TH1D*) ((TH2D*)file->Get(Form("psi_bf_0_%d_14",harmonic-1)))->ProjectionY();
  hfvtxs0W->Draw();
  hfvtxs0W->SetLineColor(kRed);
  hfvtxs0W->GetXaxis()->SetTitle(Form("#Psi_{%d}, FVTXS0",harmonic));
  TH1D* hfvtxs0 = (TH1D*) ((TH2D*)file->Get(Form("psi_bf_0_%d_4",harmonic-1)))->ProjectionY();
  hfvtxs0->Draw("same");
  hfvtxs0->SetLineColor(kGreen+2);
  hfvtxs0->GetXaxis()->SetTitle(Form("#Psi_{%d}, FVTXS0",harmonic));
  TH1D* hfvtxs0M = (TH1D*) ((TH2D*)file->Get(Form("psi_mf_0_%d_4",harmonic-1)))->ProjectionY();
  hfvtxs0M->Draw("same");
  hfvtxs0M->SetLineColor(kBlue);
  TH1D* hfvtxs0A = (TH1D*) ((TH2D*)file->Get(Form("psi_af_0_%d_4",harmonic-1)))->ProjectionY();
  hfvtxs0A->Draw("same");
  hfvtxs0A->SetLineColor(kBlack);
  leg = new TLegend(0.68,0.68,0.88,0.88);
  leg->AddEntry(hfvtxs0W,"Raw","l");
  leg->AddEntry(hfvtxs0,"Weighting","l");
  leg->AddEntry(hfvtxs0M,"Recentering","l");
  leg->AddEntry(hfvtxs0A,"Flattening","l");
  leg->SetTextSize(0.045);
  leg->Draw();
  c1->Print(Form("FigsEventPlane/after_psi%d_fvtxs0_%d.png",harmonic,run));
  c1->Print(Form("FigsEventPlane/after_psi%d_fvtxs0_%d.pdf",harmonic,run));

  TH1D* hfvtxs1W = (TH1D*) ((TH2D*)file->Get(Form("psi_bf_0_%d_15",harmonic-1)))->ProjectionY();
  hfvtxs1W->Draw();
  hfvtxs1W->SetLineColor(kRed);
  hfvtxs1W->GetXaxis()->SetTitle(Form("#Psi_{%d}, FVTXS1",harmonic));
  TH1D* hfvtxs1 = (TH1D*) ((TH2D*)file->Get(Form("psi_bf_0_%d_5",harmonic-1)))->ProjectionY();
  hfvtxs1->Draw("same");
  hfvtxs1->SetLineColor(kGreen+2);
  hfvtxs1->GetXaxis()->SetTitle(Form("#Psi_{%d}, FVTXS1",harmonic));
  TH1D* hfvtxs1M = (TH1D*) ((TH2D*)file->Get(Form("psi_mf_0_%d_5",harmonic-1)))->ProjectionY();
  hfvtxs1M->Draw("same");
  hfvtxs1M->SetLineColor(kBlue);
  TH1D* hfvtxs1A = (TH1D*) ((TH2D*)file->Get(Form("psi_af_0_%d_5",harmonic-1)))->ProjectionY();
  hfvtxs1A->Draw("same");
  hfvtxs1A->SetLineColor(kBlack);
  leg = new TLegend(0.68,0.68,0.88,0.88);
  leg->AddEntry(hfvtxs1W,"Raw","l");
  leg->AddEntry(hfvtxs1,"Weighting","l");
  leg->AddEntry(hfvtxs1M,"Recentering","l");
  leg->AddEntry(hfvtxs1A,"Flattening","l");
  leg->SetTextSize(0.045);
  leg->Draw();
  c1->Print(Form("FigsEventPlane/after_psi%d_fvtxs1_%d.png",harmonic,run));
  c1->Print(Form("FigsEventPlane/after_psi%d_fvtxs1_%d.pdf",harmonic,run));

  TH1D* hfvtxs2W = (TH1D*) ((TH2D*)file->Get(Form("psi_bf_0_%d_15",harmonic-1)))->ProjectionY();
  hfvtxs2W->Draw();
  hfvtxs2W->SetLineColor(kRed);
  hfvtxs2W->GetXaxis()->SetTitle(Form("#Psi_{%d}, FVTXS2",harmonic));
  TH1D* hfvtxs2 = (TH1D*) ((TH2D*)file->Get(Form("psi_bf_0_%d_5",harmonic-1)))->ProjectionY();
  hfvtxs2->Draw("same");
  hfvtxs2->SetLineColor(kGreen+2);
  hfvtxs2->GetXaxis()->SetTitle(Form("#Psi_{%d}, FVTXS2",harmonic));
  TH1D* hfvtxs2M = (TH1D*) ((TH2D*)file->Get(Form("psi_mf_0_%d_5",harmonic-1)))->ProjectionY();
  hfvtxs2M->Draw("same");
  hfvtxs2M->SetLineColor(kBlue);
  TH1D* hfvtxs2A = (TH1D*) ((TH2D*)file->Get(Form("psi_af_0_%d_5",harmonic-1)))->ProjectionY();
  hfvtxs2A->Draw("same");
  hfvtxs2A->SetLineColor(kBlack);
  leg = new TLegend(0.68,0.68,0.88,0.88);
  leg->AddEntry(hfvtxs2W,"Raw","l");
  leg->AddEntry(hfvtxs2,"Weighting","l");
  leg->AddEntry(hfvtxs2M,"Recentering","l");
  leg->AddEntry(hfvtxs2A,"Flattening","l");
  leg->SetTextSize(0.045);
  leg->Draw();
  c1->Print(Form("FigsEventPlane/after_psi%d_fvtxs2_%d.png",harmonic,run));
  c1->Print(Form("FigsEventPlane/after_psi%d_fvtxs2_%d.pdf",harmonic,run));

  TH1D* hfvtxs3W = (TH1D*) ((TH2D*)file->Get(Form("psi_bf_0_%d_16",harmonic-1)))->ProjectionY();
  hfvtxs3W->Draw();
  hfvtxs3W->SetLineColor(kRed);
  hfvtxs3W->GetXaxis()->SetTitle(Form("#Psi_{%d}, FVTXS3",harmonic));
  TH1D* hfvtxs3 = (TH1D*) ((TH2D*)file->Get(Form("psi_bf_0_%d_6",harmonic-1)))->ProjectionY();
  hfvtxs3->Draw("same");
  hfvtxs3->SetLineColor(kGreen+2);
  hfvtxs3->GetXaxis()->SetTitle(Form("#Psi_{%d}, FVTXS3",harmonic));
  TH1D* hfvtxs3M = (TH1D*) ((TH2D*)file->Get(Form("psi_mf_0_%d_6",harmonic-1)))->ProjectionY();
  hfvtxs3M->Draw("same");
  hfvtxs3M->SetLineColor(kBlue);
  TH1D* hfvtxs3A = (TH1D*) ((TH2D*)file->Get(Form("psi_af_0_%d_6",harmonic-1)))->ProjectionY();
  hfvtxs3A->Draw("same");
  hfvtxs3A->SetLineColor(kBlack);
  leg = new TLegend(0.68,0.68,0.88,0.88);
  leg->AddEntry(hfvtxs3W,"Raw","l");
  leg->AddEntry(hfvtxs3,"Weighting","l");
  leg->AddEntry(hfvtxs3M,"Recentering","l");
  leg->AddEntry(hfvtxs3A,"Flattening","l");
  leg->SetTextSize(0.045);
  leg->Draw();
  c1->Print(Form("FigsEventPlane/after_psi%d_fvtxs3_%d.png",harmonic,run));
  c1->Print(Form("FigsEventPlane/after_psi%d_fvtxs3_%d.pdf",harmonic,run));



  // ---

  TH1D* hfvtxnW = (TH1D*) ((TH2D*)file->Get(Form("psi_bf_0_%d_18",harmonic-1)))->ProjectionY();
  hfvtxnW->Draw();
  hfvtxnW->SetLineColor(kRed);
  hfvtxnW->GetXaxis()->SetTitle(Form("#Psi_{%d}, FVTXN",harmonic));
  TH1D* hfvtxn = (TH1D*) ((TH2D*)file->Get(Form("psi_bf_0_%d_8",harmonic-1)))->ProjectionY();
  hfvtxn->Draw("same");
  hfvtxn->SetLineColor(kGreen+2);
  hfvtxn->GetXaxis()->SetTitle(Form("#Psi_{%d}, FVTXN",harmonic));
  TH1D* hfvtxnM = (TH1D*) ((TH2D*)file->Get(Form("psi_mf_0_%d_8",harmonic-1)))->ProjectionY();
  hfvtxnM->Draw("same");
  hfvtxnM->SetLineColor(kBlue);
  TH1D* hfvtxnA = (TH1D*) ((TH2D*)file->Get(Form("psi_af_0_%d_8",harmonic-1)))->ProjectionY();
  hfvtxnA->Draw("same");
  hfvtxnA->SetLineColor(kBlack);
  leg->Draw();
  c1->Print(Form("FigsEventPlane/after_psi%d_fvtxn_%d.png",harmonic,run));
  c1->Print(Form("FigsEventPlane/after_psi%d_fvtxn_%d.pdf",harmonic,run));


  // --- and now for something related but different

  TH1D* shbbcW = (TH1D*) ((TH2D*)file->Get(Form("psi_bf_0_%d_23",harmonic-1)))->ProjectionY();
  shbbcW->Draw();
  shbbcW->SetLineColor(kRed);
  shbbcW->GetXaxis()->SetTitle(Form("#Psi_{%d}, BBCS",harmonic));
  TH1D* shbbcM = (TH1D*) ((TH2D*)file->Get(Form("psi_mf_0_%d_23",harmonic-1)))->ProjectionY();
  shbbcM->Draw("same");
  shbbcM->SetLineColor(kBlue);
  TH1D* shbbcA = (TH1D*) ((TH2D*)file->Get(Form("psi_af_0_%d_23",harmonic-1)))->ProjectionY();
  shbbcA->Draw("same");
  shbbcA->SetLineColor(kBlack);
  leg = new TLegend(0.68,0.68,0.88,0.88);
  leg->AddEntry(shbbcW,"Raw","l");
  leg->AddEntry(shbbcM,"Recentering","l");
  leg->AddEntry(shbbcA,"Flattening","l");
  leg->SetTextSize(0.045);
  leg->Draw();
  c1->Print(Form("FigsEventPlane/specafter_psi%d_bbcs_%d.png",harmonic,run));
  c1->Print(Form("FigsEventPlane/specafter_psi%d_bbcs_%d.pdf",harmonic,run));
  delete leg;

  TH1D* shfvtxsW = (TH1D*) ((TH2D*)file->Get(Form("psi_bf_0_%d_13",harmonic-1)))->ProjectionY();
  shfvtxsW->Draw();
  shfvtxsW->SetLineColor(kRed);
  shfvtxsW->GetXaxis()->SetTitle(Form("#Psi_{%d}, FVTXS",harmonic));
  TH1D* shfvtxsM = (TH1D*) ((TH2D*)file->Get(Form("psi_mf_0_%d_13",harmonic-1)))->ProjectionY();
  shfvtxsM->Draw("same");
  shfvtxsM->SetLineColor(kBlue);
  TH1D* shfvtxsA = (TH1D*) ((TH2D*)file->Get(Form("psi_af_0_%d_13",harmonic-1)))->ProjectionY();
  shfvtxsA->Draw("same");
  shfvtxsA->SetLineColor(kBlack);
  leg = new TLegend(0.68,0.68,0.88,0.88);
  leg->AddEntry(shfvtxsW,"Raw","l");
  leg->AddEntry(shfvtxsM,"Recentering","l");
  leg->AddEntry(shfvtxsA,"Flattening","l");
  leg->SetTextSize(0.045);
  leg->Draw();
  c1->Print(Form("FigsEventPlane/specafter_psi%d_fvtxs_%d.png",harmonic,run));
  c1->Print(Form("FigsEventPlane/specafter_psi%d_fvtxs_%d.pdf",harmonic,run));

  TH1D* shfvtxnW = (TH1D*) ((TH2D*)file->Get(Form("psi_bf_0_%d_18",harmonic-1)))->ProjectionY();
  shfvtxnW->Draw();
  shfvtxnW->SetLineColor(kRed);
  shfvtxnW->GetXaxis()->SetTitle(Form("#Psi_{%d}, FVTXN",harmonic));
  TH1D* shfvtxnM = (TH1D*) ((TH2D*)file->Get(Form("psi_mf_0_%d_18",harmonic-1)))->ProjectionY();
  shfvtxnM->Draw("same");
  shfvtxnM->SetLineColor(kBlue);
  TH1D* shfvtxnA = (TH1D*) ((TH2D*)file->Get(Form("psi_af_0_%d_18",harmonic-1)))->ProjectionY();
  shfvtxnA->Draw("same");
  shfvtxnA->SetLineColor(kBlack);
  leg->Draw();
  c1->Print(Form("FigsEventPlane/specafter_psi%d_fvtxn_%d.png",harmonic,run));
  c1->Print(Form("FigsEventPlane/specafter_psi%d_fvtxn_%d.pdf",harmonic,run));

  // ---

  delete c1;

  delete leg;

}

void psidiff_run(int run, int harmonic)
{

  gStyle->SetOptTitle(1);

  TCanvas* c1 = new TCanvas("c1","");

  TFile* file = TFile::Open(Form("input/hist_%d.root",run));

  // ---

  TH1D* th1d_dreso_bbc_cnt = (TH1D*)file->Get(Form("th1d_dreso%d_BBC_CNT",harmonic));
  TH1D* th1d_dreso_bbc_fvtx = (TH1D*)file->Get(Form("th1d_dreso%d_BBC_FVTX",harmonic));
  TH1D* th1d_dreso_cnt_fvtx = (TH1D*)file->Get(Form("th1d_dreso%d_CNT_FVTX",harmonic));
  TH1D* th1d_dreso_bbc_fvtxn = (TH1D*)file->Get(Form("th1d_dreso%d_BBC_FVTXN",harmonic));
  TH1D* th1d_dreso_cnt_fvtxn = (TH1D*)file->Get(Form("th1d_dreso%d_CNT_FVTXN",harmonic));
  TH1D* th1d_dreso_fvtxs_fvtxn = (TH1D*)file->Get(Form("th1d_dreso%d_FVTXS_FVTXN",harmonic));

  th1d_dreso_bbc_cnt->SetLineColor(kBlue);
  th1d_dreso_bbc_fvtx->SetLineColor(kBlack);
  th1d_dreso_cnt_fvtx->SetLineColor(kRed);
  th1d_dreso_bbc_fvtxn->SetLineColor(kGreen+2);
  th1d_dreso_cnt_fvtxn->SetLineColor(kOrange+2);
  th1d_dreso_fvtxs_fvtxn->SetLineColor(kGray);

  double max = 0;
  double max_BC = 0;
  double max_BN = 0;
  double max_BS = 0;
  double max_CN = 0;
  double max_CS = 0;
  double max_NS = 0;

  max_BC = th1d_dreso_bbc_cnt->GetMaximum();
  max_BN = th1d_dreso_bbc_fvtxn->GetMaximum();
  max_BS = th1d_dreso_bbc_fvtx->GetMaximum();
  max_CN = th1d_dreso_cnt_fvtxn->GetMaximum();
  max_CS = th1d_dreso_cnt_fvtx->GetMaximum();
  max_NS = th1d_dreso_fvtxs_fvtxn->GetMaximum();

  max = 1;
  if ( max_BC > max ) max = max_BC;
  if ( max_BN > max ) max = max_BN;
  if ( max_BS > max ) max = max_BS;
  if ( max_CN > max ) max = max_CN;
  if ( max_CS > max ) max = max_CS;
  if ( max_NS > max ) max = max_NS;

  TH2* th2 = new TH2D("th2","",1,-6.0,6.0,1,0.0,1.1*max);
  th2->Draw();
  th1d_dreso_bbc_cnt->Draw("same");
  th1d_dreso_bbc_fvtx->Draw("same");
  th1d_dreso_cnt_fvtx->Draw("same");
  th1d_dreso_bbc_fvtxn->Draw("same");
  th1d_dreso_cnt_fvtxn->Draw("same");
  th1d_dreso_fvtxs_fvtxn->Draw("same");
  TLegend* leg = new TLegend();
  c1->Print(Form("FigsEventPlane/psidiff_dreso%d_%d.png",harmonic,run));
  c1->Print(Form("FigsEventPlane/psidiff_dreso%d_%d.pdf",harmonic,run));
  delete th2;
  delete leg;

  // ---

  th1d_dreso_bbc_fvtx->SetLineColor(kRed);
  th1d_dreso_bbc_fvtxn->SetLineColor(kGreen+2);
  th1d_dreso_fvtxs_fvtxn->SetLineColor(kBlue);

  th1d_dreso_bbc_cnt->SetLineColor(kBlue);
  th1d_dreso_cnt_fvtx->SetLineColor(kRed);
  th1d_dreso_cnt_fvtxn->SetLineColor(kBlack);

  // ---

  max = 1;
  if ( max_BN > max ) max = max_BN;
  if ( max_BS > max ) max = max_BS;
  if ( max_NS > max ) max = max_NS;

  th2 = new TH2D("th2","",1,-6.0,6.0,1,0.0,1.1*max);
  th2->Draw();
  th1d_dreso_bbc_fvtx->Draw("same");
  th1d_dreso_bbc_fvtxn->Draw("same");
  th1d_dreso_fvtxs_fvtxn->Draw("same");
  leg = new TLegend(0.18,0.68,0.38,0.88);
  leg->SetHeader(Form("Run %d #Psi_{%d}",run,harmonic));
  leg->AddEntry(th1d_dreso_bbc_fvtx,"BBCS-FVTXS","l");
  leg->AddEntry(th1d_dreso_bbc_fvtxn,"BBCS-FVTXN","l");
  leg->AddEntry(th1d_dreso_fvtxs_fvtxn,"FVTXS-FVTXN","l");
  leg->SetTextSize(0.045);
  leg->Draw();
  c1->Print(Form("FigsEventPlane/psidiff_dreso%d_%d_SANSCNT.png",harmonic,run));
  c1->Print(Form("FigsEventPlane/psidiff_dreso%d_%d_SANSCNT.pdf",harmonic,run));
  delete th2;

  max = 1;
  if ( max_BC > max ) max = max_BC;
  if ( max_CN > max ) max = max_CN;
  if ( max_CS > max ) max = max_CS;

  th2 = new TH2D("th2","",1,-6.0,6.0,1,0.0,1.1*max);
  th2->Draw();
  th1d_dreso_bbc_cnt->Draw("same");
  th1d_dreso_cnt_fvtx->Draw("same");
  th1d_dreso_cnt_fvtxn->Draw("same");
  leg = new TLegend(0.18,0.68,0.38,0.88);
  leg->SetHeader(Form("Run %d #Psi_{%d}",run,harmonic));
  leg->AddEntry(th1d_dreso_bbc_cnt,"BBCS-CNT","l");
  leg->AddEntry(th1d_dreso_cnt_fvtx,"CNT-FVTXS","l");
  leg->AddEntry(th1d_dreso_cnt_fvtxn,"CNT-FVTXN","l");
  leg->SetTextSize(0.045);
  leg->Draw();
  c1->Print(Form("FigsEventPlane/psidiff_dreso%d_%d_AVECCNT.png",harmonic,run));
  c1->Print(Form("FigsEventPlane/psidiff_dreso%d_%d_AVECCNT.pdf",harmonic,run));
  delete th2;



  delete c1;

}



void cospsidiff_run(int run, int harmonic)
{

  gStyle->SetOptTitle(1);

  TCanvas* c1 = new TCanvas("c1","");

  TFile* file = TFile::Open(Form("input/hist_%d.root",run));

  // ---

  TH1D* th1d_reso_bbc_cnt = (TH1D*)file->Get(Form("th1d_reso%d_BBC_CNT",harmonic));
  TH1D* th1d_reso_bbc_fvtx = (TH1D*)file->Get(Form("th1d_reso%d_BBC_FVTX",harmonic));
  TH1D* th1d_reso_cnt_fvtx = (TH1D*)file->Get(Form("th1d_reso%d_CNT_FVTX",harmonic));
  TH1D* th1d_reso_bbc_fvtxn = (TH1D*)file->Get(Form("th1d_reso%d_BBC_FVTXN",harmonic));
  TH1D* th1d_reso_cnt_fvtxn = (TH1D*)file->Get(Form("th1d_reso%d_CNT_FVTXN",harmonic));
  TH1D* th1d_reso_fvtxs_fvtxn = (TH1D*)file->Get(Form("th1d_reso%d_FVTXS_FVTXN",harmonic));

  th1d_reso_bbc_cnt->SetLineColor(kBlue);
  th1d_reso_bbc_fvtx->SetLineColor(kBlack);
  th1d_reso_cnt_fvtx->SetLineColor(kRed);
  th1d_reso_bbc_fvtxn->SetLineColor(kGreen+2);
  th1d_reso_cnt_fvtxn->SetLineColor(kOrange+2);
  th1d_reso_fvtxs_fvtxn->SetLineColor(kGray);

  double max = 0;
  double max_BC = 0;
  double max_BN = 0;
  double max_BS = 0;
  double max_CN = 0;
  double max_CS = 0;
  double max_NS = 0;

  max_BC = th1d_reso_bbc_cnt->GetMaximum();
  max_BN = th1d_reso_bbc_fvtxn->GetMaximum();
  max_BS = th1d_reso_bbc_fvtx->GetMaximum();
  max_CN = th1d_reso_cnt_fvtxn->GetMaximum();
  max_CS = th1d_reso_cnt_fvtx->GetMaximum();
  max_NS = th1d_reso_fvtxs_fvtxn->GetMaximum();

  max = 1;
  if ( max_BC > max ) max = max_BC;
  if ( max_BN > max ) max = max_BN;
  if ( max_BS > max ) max = max_BS;
  if ( max_CN > max ) max = max_CN;
  if ( max_CS > max ) max = max_CS;
  if ( max_NS > max ) max = max_NS;

  TH2* th2 = new TH2D("th2","",1,-1.1,1.1,1,0.0,1.1*max);
  th2->Draw();
  th1d_reso_bbc_cnt->Draw("same");
  th1d_reso_bbc_fvtx->Draw("same");
  th1d_reso_cnt_fvtx->Draw("same");
  th1d_reso_bbc_fvtxn->Draw("same");
  th1d_reso_cnt_fvtxn->Draw("same");
  th1d_reso_fvtxs_fvtxn->Draw("same");
  TLegend* leg = new TLegend();
  c1->Print(Form("FigsEventPlane/psidiff_reso%d_%d.png",harmonic,run));
  c1->Print(Form("FigsEventPlane/psidiff_reso%d_%d.pdf",harmonic,run));
  delete th2;
  delete leg;

  // ---

  th1d_reso_bbc_fvtx->SetLineColor(kRed);
  th1d_reso_bbc_fvtxn->SetLineColor(kGreen+2);
  th1d_reso_fvtxs_fvtxn->SetLineColor(kBlue);

  th1d_reso_bbc_cnt->SetLineColor(kBlue);
  th1d_reso_cnt_fvtx->SetLineColor(kRed);
  th1d_reso_cnt_fvtxn->SetLineColor(kBlack);

  // ---

  max = 1;
  if ( max_BN > max ) max = max_BN;
  if ( max_BS > max ) max = max_BS;
  if ( max_NS > max ) max = max_NS;

  th2 = new TH2D("th2","",1,-1.1,1.1,1,0.0,1.1*max);
  th2->Draw();
  th1d_reso_bbc_fvtx->Draw("same");
  th1d_reso_bbc_fvtxn->Draw("same");
  th1d_reso_fvtxs_fvtxn->Draw("same");
  leg = new TLegend(0.25,0.68,0.38,0.88);
  leg->SetHeader(Form("Run %d #Psi_{%d}",run,harmonic));
  leg->AddEntry(th1d_reso_bbc_fvtx,"BBCS-FVTXS","l");
  leg->AddEntry(th1d_reso_bbc_fvtxn,"BBCS-FVTXN","l");
  leg->AddEntry(th1d_reso_fvtxs_fvtxn,"FVTXS-FVTXN","l");
  leg->SetTextSize(0.045);
  leg->Draw();
  c1->Print(Form("FigsEventPlane/psidiff_reso%d_%d_SANSCNT.png",harmonic,run));
  c1->Print(Form("FigsEventPlane/psidiff_reso%d_%d_SANSCNT.pdf",harmonic,run));
  delete th2;

  max = 1;
  if ( max_BC > max ) max = max_BC;
  if ( max_CN > max ) max = max_CN;
  if ( max_CS > max ) max = max_CS;

  th2 = new TH2D("th2","",1,-1.1,1.1,1,0.0,1.1*max);
  th2->Draw();
  th1d_reso_bbc_cnt->Draw("same");
  th1d_reso_cnt_fvtx->Draw("same");
  th1d_reso_cnt_fvtxn->Draw("same");
  leg = new TLegend(0.25,0.68,0.38,0.88);
  leg->SetHeader(Form("Run %d #Psi_{%d}",run,harmonic));
  leg->AddEntry(th1d_reso_bbc_cnt,"BBCS-CNT","l");
  leg->AddEntry(th1d_reso_cnt_fvtx,"CNT-FVTXS","l");
  leg->AddEntry(th1d_reso_cnt_fvtxn,"CNT-FVTXN","l");
  leg->SetTextSize(0.045);
  leg->Draw();
  c1->Print(Form("FigsEventPlane/psidiff_reso%d_%d_AVECCNT.png",harmonic,run));
  c1->Print(Form("FigsEventPlane/psidiff_reso%d_%d_AVECCNT.pdf",harmonic,run));
  delete th2;



  delete c1;

}

