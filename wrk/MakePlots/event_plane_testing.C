void psidiff_run(int);


void event_plane_testing()
{

  psidiff_run(455050);

}


void psidiff_run(int run)
{

  gStyle->SetOptTitle(1);

  TCanvas* c1 = new TCanvas("c1","");

  TFile* file = TFile::Open(Form("input/hist_%d.root",run));

  // ---
  // ---
  // ---

  TProfile* tp1f_reso2bbc_cnt = (TProfile*)file->Get("tp1f_reso2_BBC_CNT");
  TProfile* tp1f_reso2bbc_fvtx = (TProfile*)file->Get("tp1f_reso2_BBC_FVTX");
  TProfile* tp1f_reso2cnt_fvtx = (TProfile*)file->Get("tp1f_reso2_CNT_FVTX");

  double double_reso2_bbc_cnt = tp1f_reso2bbc_cnt->GetBinContent(1);
  double double_reso2_bbc_fvtx = tp1f_reso2bbc_fvtx->GetBinContent(1);
  double double_reso2_cnt_fvtx = tp1f_reso2cnt_fvtx->GetBinContent(1);

  cout << "BBC CNT correlation " << double_reso2_bbc_cnt << endl;
  cout << "BBC FVTX correlation " << double_reso2_bbc_fvtx << endl;
  cout << "CNT FVTX correlation " << double_reso2_cnt_fvtx << endl;

  double reso_bbc = sqrt((double_reso2_bbc_cnt*double_reso2_bbc_fvtx)/double_reso2_cnt_fvtx);
  double reso_fvtx = sqrt((double_reso2_cnt_fvtx*double_reso2_bbc_fvtx)/double_reso2_bbc_cnt);

  cout << "bbc resolution is " << reso_bbc << endl;
  cout << "fvtx resolution is " << reso_fvtx << endl;

  // ---

  TH1D* th1d_reso2bbc_cnt = (TH1D*)file->Get("th1d_reso2_BBC_CNT");
  TH1D* th1d_reso2bbc_fvtx = (TH1D*)file->Get("th1d_reso2_BBC_FVTX");
  TH1D* th1d_reso2cnt_fvtx = (TH1D*)file->Get("th1d_reso2_CNT_FVTX");

  th1d_reso2bbc_cnt->SetLineColor(kBlue);
  th1d_reso2bbc_fvtx->SetLineColor(kBlack);
  th1d_reso2cnt_fvtx->SetLineColor(kRed);

  th1d_reso2bbc_cnt->Draw();
  th1d_reso2bbc_fvtx->Draw("same");
  th1d_reso2cnt_fvtx->Draw("same");
  th1d_reso2bbc_cnt->SetMinimum(th1d_reso2bbc_fvtx->GetMinimum(1.0));

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

  th1d_dreso2bbc_cnt->Draw();
  th1d_dreso2bbc_fvtx->Draw("same");
  th1d_dreso2cnt_fvtx->Draw("same");

  c1->SetLogy(0);
  c1->Print(Form("psidiff_dreso2_%d.png",run));
  c1->Print(Form("psidiff_dreso2_%d.pdf",run));
  c1->SetLogy(1);

  // ---
  // ---
  // ---

  TProfile* tp1f_reso3bbc_cnt = (TProfile*)file->Get("tp1f_reso3_BBC_CNT");
  TProfile* tp1f_reso3bbc_fvtx = (TProfile*)file->Get("tp1f_reso3_BBC_FVTX");
  TProfile* tp1f_reso3cnt_fvtx = (TProfile*)file->Get("tp1f_reso3_CNT_FVTX");

  double double_reso3_bbc_cnt = tp1f_reso3bbc_cnt->GetBinContent(1);
  double double_reso3_bbc_fvtx = tp1f_reso3bbc_fvtx->GetBinContent(1);
  double double_reso3_cnt_fvtx = tp1f_reso3cnt_fvtx->GetBinContent(1);

  cout << "BBC CNT correlation " << double_reso3_bbc_cnt << endl;
  cout << "BBC FVTX correlation " << double_reso3_bbc_fvtx << endl;
  cout << "CNT FVTX correlation " << double_reso3_cnt_fvtx << endl;

  double reso_bbc = sqrt((double_reso3_bbc_cnt*double_reso3_bbc_fvtx)/double_reso3_cnt_fvtx);
  double reso_fvtx = sqrt((double_reso3_cnt_fvtx*double_reso3_bbc_fvtx)/double_reso3_bbc_cnt);

  cout << "bbc resolution is " << reso_bbc << endl;
  cout << "fvtx resolution is " << reso_fvtx << endl;

  // ---

  TH1D* th1d_reso3bbc_cnt = (TH1D*)file->Get("th1d_reso3_BBC_CNT");
  TH1D* th1d_reso3bbc_fvtx = (TH1D*)file->Get("th1d_reso3_BBC_FVTX");
  TH1D* th1d_reso3cnt_fvtx = (TH1D*)file->Get("th1d_reso3_CNT_FVTX");

  th1d_reso3bbc_cnt->SetLineColor(kBlue);
  th1d_reso3bbc_fvtx->SetLineColor(kBlack);
  th1d_reso3cnt_fvtx->SetLineColor(kRed);

  th1d_reso3bbc_cnt->Draw();
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

  th1d_dreso3bbc_cnt->Draw();
  th1d_dreso3bbc_fvtx->Draw("same");
  th1d_dreso3cnt_fvtx->Draw("same");

  c1->SetLogy(0);
  c1->Print(Form("psidiff_dreso3_%d.png",run));
  c1->Print(Form("psidiff_dreso3_%d.pdf",run));
  c1->SetLogy(1);

}

