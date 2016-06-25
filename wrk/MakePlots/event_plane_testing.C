void psidiff_run(int);


void event_plane_testing()
{

  psidiff_run(455050);
  psidiff_run(456015);
  psidiff_run(457015);
  psidiff_run(458014);

}


void psidiff_run(int run)
{

  gStyle->SetOptTitle(1);

  TCanvas* c1 = new TCanvas("c1","");

  TFile* file = TFile::Open(Form("input/hist_%d.root",run));

  // ---
  // ---
  // ---

  TH1D* th1d_reso2bbc_cnt = (TH1D*)file->Get("th1d_reso2_BBC_CNT");
  TH1D* th1d_reso2bbc_fvtx = (TH1D*)file->Get("th1d_reso2_BBC_FVTX");
  TH1D* th1d_reso2cnt_fvtx = (TH1D*)file->Get("th1d_reso2_CNT_FVTX");

  th1d_reso2bbc_cnt->SetLineColor(kBlue);
  th1d_reso2bbc_fvtx->SetLineColor(kBlack);
  th1d_reso2cnt_fvtx->SetLineColor(kRed);

  double max_BC = th1d_reso2bbc_cnt->GetMaximum();
  double max_BS = th1d_reso2bbc_fvtx->GetMaximum();
  double max_CS = th1d_reso2cnt_fvtx->GetMaximum();

  double max = 1.0;
  if ( max_BC > max ) max = max_BC;
  if ( max_BS > max ) max = max_BS;
  if ( max_CS > max ) max = max_CS;

  double min_BC = th1d_reso2bbc_cnt->GetMinimum(1);
  double min_BS = th1d_reso2bbc_fvtx->GetMinimum(1);
  double min_CS = th1d_reso2cnt_fvtx->GetMinimum(1);

  double min = min_BC;
  if ( min_BS < min ) min = min_BS;
  if ( min_CS < min ) min = min_CS;

  TH2D* th2 = new TH2D("th2","",1,-1.1,1.1,1,0.9*min,1.1*max);
  th2->Draw();
  th1d_reso2bbc_cnt->Draw("same");
  th1d_reso2bbc_fvtx->Draw("same");
  th1d_reso2cnt_fvtx->Draw("same");


  c1->SetLogy(1);
  c1->Print(Form("psidiff_reso2_%d.png",run));
  c1->Print(Form("psidiff_reso2_%d.pdf",run));
  c1->SetLogy(0);

  // ---

  TH1D* th1d_dreso2bbc_cnt = (TH1D*)file->Get("th1d_dreso2_BBC_CNT");
  TH1D* th1d_dreso2bbc_fvtx = (TH1D*)file->Get("th1d_dreso2_BBC_FVTX");
  TH1D* th1d_dreso2cnt_fvtx = (TH1D*)file->Get("th1d_dreso2_CNT_FVTX");

  th1d_dreso2bbc_cnt->SetLineColor(kBlue);
  th1d_dreso2bbc_fvtx->SetLineColor(kBlack);
  th1d_dreso2cnt_fvtx->SetLineColor(kRed);

  max_BC = th1d_dreso2bbc_cnt->GetMaximum();
  max_BS = th1d_dreso2bbc_fvtx->GetMaximum();
  max_CS = th1d_dreso2cnt_fvtx->GetMaximum();

  max = 1;
  if ( max_BC > max ) max = max_BC;
  if ( max_BS > max ) max = max_BS;
  if ( max_CS > max ) max = max_CS;

  delete th2;
  th2 = new TH2D("th2","",1,-6.0,6.0,1,0.0,1.1*max);
  th2->Draw();
  th1d_dreso2bbc_cnt->Draw("same");
  th1d_dreso2bbc_fvtx->Draw("same");
  th1d_dreso2cnt_fvtx->Draw("same");

  c1->SetLogy(0);
  c1->Print(Form("psidiff_dreso2_%d.png",run));
  c1->Print(Form("psidiff_dreso2_%d.pdf",run));
  c1->SetLogy(1);

  // ---
  // ---
  // ---

  TH1D* th1d_reso3bbc_cnt = (TH1D*)file->Get("th1d_reso3_BBC_CNT");
  TH1D* th1d_reso3bbc_fvtx = (TH1D*)file->Get("th1d_reso3_BBC_FVTX");
  TH1D* th1d_reso3cnt_fvtx = (TH1D*)file->Get("th1d_reso3_CNT_FVTX");

  th1d_reso3bbc_cnt->SetLineColor(kBlue);
  th1d_reso3bbc_fvtx->SetLineColor(kBlack);
  th1d_reso3cnt_fvtx->SetLineColor(kRed);

  max_BC = th1d_reso3bbc_cnt->GetMaximum();
  max_BS = th1d_reso3bbc_fvtx->GetMaximum();
  max_CS = th1d_reso3cnt_fvtx->GetMaximum();

  max = 1;
  if ( max_BC > max ) max = max_BC;
  if ( max_BS > max ) max = max_BS;
  if ( max_CS > max ) max = max_CS;

  min_BC = th1d_reso3bbc_cnt->GetMinimum(1);
  min_BS = th1d_reso3bbc_fvtx->GetMinimum(1);
  min_CS = th1d_reso3cnt_fvtx->GetMinimum(1);

  min = min_BC;
  if ( min_BS < min ) min = min_BS;
  if ( min_CS < min ) min = min_CS;

  delete th2;
  th2 = new TH2D("th2","",1,-1.1,1.1,1,0.9*min,1.1*max);
  th2->Draw();
  th1d_reso3bbc_cnt->Draw("same");
  th1d_reso3bbc_fvtx->Draw("same");
  th1d_reso3cnt_fvtx->Draw("same");
  th1d_reso3bbc_cnt->SetMinimum(th1d_reso3bbc_fvtx->GetMinimum(1.0));

  c1->SetLogy(1);
  c1->Print(Form("psidiff_reso3_%d.png",run));
  c1->Print(Form("psidiff_reso3_%d.pdf",run));
  c1->SetLogy(0);

  // ---

  TH1D* th1d_dreso3bbc_cnt = (TH1D*)file->Get("th1d_dreso3_BBC_CNT");
  TH1D* th1d_dreso3bbc_fvtx = (TH1D*)file->Get("th1d_dreso3_BBC_FVTX");
  TH1D* th1d_dreso3cnt_fvtx = (TH1D*)file->Get("th1d_dreso3_CNT_FVTX");

  th1d_dreso3bbc_cnt->SetLineColor(kBlue);
  th1d_dreso3bbc_fvtx->SetLineColor(kBlack);
  th1d_dreso3cnt_fvtx->SetLineColor(kRed);

  max_BC = th1d_dreso3bbc_cnt->GetMaximum();
  max_BS = th1d_dreso3bbc_fvtx->GetMaximum();
  max_CS = th1d_dreso3cnt_fvtx->GetMaximum();

  max = 1;
  if ( max_BC > max ) max = max_BC;
  if ( max_BS > max ) max = max_BS;
  if ( max_CS > max ) max = max_CS;

  delete th2;
  th2 = new TH2D("th2","",1,-5.0,5.0,1,0.0,1.1*max);
  th2->Draw();
  th1d_dreso3bbc_cnt->Draw("same");
  th1d_dreso3bbc_fvtx->Draw("same");
  th1d_dreso3cnt_fvtx->Draw("same");

  c1->SetLogy(0);
  c1->Print(Form("psidiff_dreso3_%d.png",run));
  c1->Print(Form("psidiff_dreso3_%d.pdf",run));
  //c1->SetLogy(1); //???

  delete c1;

}

