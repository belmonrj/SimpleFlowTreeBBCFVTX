void doit(int);

void temp_2pc()
{

  //doit(456652);
  doit(200);
  doit(62);
  doit(39);
  doit(20);

}

void doit(int handle)
{

  gStyle->SetOptTitle(1);

  TCanvas* c1 = new TCanvas("c1","");

  TFile* file = NULL;
  if ( handle <= 200 ) file = TFile::Open(Form("input/combined_%d.root",handle));
  else if ( handle > 454000 ) file = TFile::Open(Form("input/hist_%d.root",handle));
  else
    {
      cout << "YOU'RE GONNA DIE" << endl;
      return;
    }

  // ---
  // ---
  // ---

  // ------------
  // --- BBCS ---
  // ------------

  TProfile* tp1f_os_bbcs_d22_both = (TProfile*)file->Get("os_bbcs_d22_both");
  TProfile* tp1f_os_bbcs_cos22_both = (TProfile*)file->Get("os_bbcs_cos22_both");
  TProfile* tp1f_os_bbcs_sin22_both = (TProfile*)file->Get("os_bbcs_sin22_both");

  TProfile* tp1f_os_bbcs_c22 = (TProfile*)file->Get("os_bbcs_c22");
  TProfile* tp1f_os_bbcs_cos22 = (TProfile*)file->Get("os_bbcs_cos22");
  TProfile* tp1f_os_bbcs_sin22 = (TProfile*)file->Get("os_bbcs_sin22");

  double os_bbcs_c22_raw = tp1f_os_bbcs_c22->GetBinContent(1);
  double os_bbcs_cos22_raw = tp1f_os_bbcs_cos22->GetBinContent(1);
  double os_bbcs_sin22_raw = tp1f_os_bbcs_sin22->GetBinContent(1);

  double os_bbcs_c22_corr = os_bbcs_c22_raw;
  if ( os_bbcs_c22_corr < 0 ) cout << "YOU'RE GONNA DIE" << endl;
  double os_bbcs_v22_corr = sqrt(fabs(os_bbcs_c22_corr));
  cout << os_bbcs_v22_corr << endl;

  TH1D* th1d_os_bbcs_d22_both = (TH1D*)tp1f_os_bbcs_d22_both->ProjectionX();
  TH1D* th1d_os_bbcs_d22_both_corr = (TH1D*)th1d_os_bbcs_d22_both->Clone();
  th1d_os_bbcs_d22_both_corr->Draw();
  c1->Print(Form("FigsTwo/corr_os_bbcs_d22_%d.png",handle));
  c1->Print(Form("FigsTwo/corr_os_bbcs_d22_%d.pdf",handle));

  th1d_os_bbcs_d22_both_corr->Scale(1.0/os_bbcs_v22_corr);
  if ( handle <= 200 ) th1d_os_bbcs_d22_both_corr->SetTitle(Form("d+Au #sqrt{s_{NN}} = %d GeV",handle));
  th1d_os_bbcs_d22_both_corr->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  th1d_os_bbcs_d22_both_corr->GetYaxis()->SetTitle("v_{2}{2}");
  th1d_os_bbcs_d22_both_corr->Draw();
  c1->Print(Form("FigsTwo/corr_os_bbcs_v22_%d.png",handle));
  c1->Print(Form("FigsTwo/corr_os_bbcs_v22_%d.pdf",handle));



  // -------------
  // --- FVTXS ---
  // -------------

  TProfile* tp1f_os_fvtxs_d22_both = (TProfile*)file->Get("os_fvtxs_d22_both");
  TProfile* tp1f_os_fvtxs_cos22_both = (TProfile*)file->Get("os_fvtxs_cos22_both");
  TProfile* tp1f_os_fvtxs_sin22_both = (TProfile*)file->Get("os_fvtxs_sin22_both");

  TProfile* tp1f_os_fvtxs_c22 = (TProfile*)file->Get("os_fvtxs_c22");
  TProfile* tp1f_os_fvtxs_cos22 = (TProfile*)file->Get("os_fvtxs_cos22");
  TProfile* tp1f_os_fvtxs_sin22 = (TProfile*)file->Get("os_fvtxs_sin22");

  double os_fvtxs_c22_raw = tp1f_os_fvtxs_c22->GetBinContent(1);
  double os_fvtxs_cos22_raw = tp1f_os_fvtxs_cos22->GetBinContent(1);
  double os_fvtxs_sin22_raw = tp1f_os_fvtxs_sin22->GetBinContent(1);

  double os_fvtxs_c22_corr = os_fvtxs_c22_raw;
  if ( os_fvtxs_c22_corr < 0 ) cout << "YOU'RE GONNA DIE" << endl;
  double os_fvtxs_v22_corr = sqrt(fabs(os_fvtxs_c22_corr));
  cout << os_fvtxs_v22_corr << endl;

  TH1D* th1d_os_fvtxs_d22_both = (TH1D*)tp1f_os_fvtxs_d22_both->ProjectionX();
  TH1D* th1d_os_fvtxs_d22_both_corr = (TH1D*)th1d_os_fvtxs_d22_both->Clone();
  th1d_os_fvtxs_d22_both_corr->Draw();
  c1->Print(Form("FigsTwo/corr_os_fvtxs_d22_%d.png",handle));
  c1->Print(Form("FigsTwo/corr_os_fvtxs_d22_%d.pdf",handle));

  th1d_os_fvtxs_d22_both_corr->Scale(1.0/os_fvtxs_v22_corr);
  if ( handle <= 200 ) th1d_os_fvtxs_d22_both_corr->SetTitle(Form("d+Au #sqrt{s_{NN}} = %d GeV",handle));
  th1d_os_fvtxs_d22_both_corr->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  th1d_os_fvtxs_d22_both_corr->GetYaxis()->SetTitle("v_{2}{2}");
  th1d_os_fvtxs_d22_both_corr->Draw();
  c1->Print(Form("FigsTwo/corr_os_fvtxs_v22_%d.png",handle));
  c1->Print(Form("FigsTwo/corr_os_fvtxs_v22_%d.pdf",handle));

  delete c1;

}
