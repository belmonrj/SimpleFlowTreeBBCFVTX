void doit(int);

void temp_4pc()
{

  doit(200);

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

  TProfile* tp1f_os_bbcs_d22_both = (TProfile*)file->Get("os_bbcs_d22_both");
  TProfile* tp1f_os_bbcs_c22 = (TProfile*)file->Get("os_bbcs_c22");
  double os_bbcs_c22 = tp1f_os_bbcs_c22->GetBinContent(1);
  if ( os_bbcs_c22 < 0 ) cout << "YOU'RE GONNA DIE" << endl;
  double os_bbcs_v22 = sqrt(fabs(os_bbcs_c22));
  cout << os_bbcs_v22 << endl;

  TProfile* tp1f_os_bbcs_d24_out_both = (TProfile*)file->Get("os_bbcs_d24_out_both");
  TProfile* tp1f_os_bbcs_c24 = (TProfile*)file->Get("os_bbcs_c24");
  double os_bbcs_c24 = tp1f_os_bbcs_c24->GetBinContent(1);

  double os_bbcs_cumulant = 2*os_bbcs_c22*os_bbcs_c22 - os_bbcs_c24;
  cout << os_bbcs_c24 << endl;
  cout << 2*os_bbcs_c22*os_bbcs_c22 << endl;
  if ( os_bbcs_cumulant < 0 ) cout << "YOU'RE GONNA DIE" << endl;
  double os_bbcs_v24 = pow(fabs(os_bbcs_c24),0.25);
  cout << os_bbcs_v24 << endl;

  tp1f_os_bbcs_d24_out_both->Draw();
  c1->Print("fourpc_bbcs_pt.png");

  TH1D* th1d_os_bbcs_d22_both = tp1f_os_bbcs_d22_both->ProjectionX();
  TH1D* th1d_os_bbcs_d24_out_both = tp1f_os_bbcs_d24_out_both->ProjectionX();
  th1d_os_bbcs_d22_both->Scale(2*os_bbcs_c22);
  th1d_os_bbcs_d22_both->Add(th1d_os_bbcs_d24_out_both,-1.0);
  th1d_os_bbcs_d22_both->Draw();
  c1->Print("cumulant4pc_bbcs_pt.png");

  // ---

  TProfile* tp1f_os_fvtxs_d22_both = (TProfile*)file->Get("os_fvtxs_d22_both");
  TProfile* tp1f_os_fvtxs_c22 = (TProfile*)file->Get("os_fvtxs_c22");
  double os_fvtxs_c22 = tp1f_os_fvtxs_c22->GetBinContent(1);
  if ( os_fvtxs_c22 < 0 ) cout << "YOU'RE GONNA DIE" << endl;
  double os_fvtxs_v22 = sqrt(fabs(os_fvtxs_c22));
  cout << os_fvtxs_v22 << endl;

  TProfile* tp1f_os_fvtxs_d24_out_both = (TProfile*)file->Get("os_fvtxs_d24_out_both");
  TProfile* tp1f_os_fvtxs_c24 = (TProfile*)file->Get("os_fvtxs_c24");
  double os_fvtxs_c24 = tp1f_os_fvtxs_c24->GetBinContent(1);

  double os_fvtxs_cumulant = 2*os_fvtxs_c22*os_fvtxs_c22 - os_fvtxs_c24;
  cout << os_fvtxs_c24 << endl;
  cout << 2*os_fvtxs_c22*os_fvtxs_c22 << endl;
  if ( os_fvtxs_cumulant < 0 ) cout << "YOU'RE GONNA DIE" << endl;
  double os_fvtxs_v24 = pow(fabs(os_fvtxs_c24),0.25);
  cout << os_fvtxs_v24 << endl;

  tp1f_os_fvtxs_d24_out_both->Draw();
  c1->Print("fourpc_fvtxs_pt.png");

  TH1D* th1d_os_fvtxs_d22_both = tp1f_os_fvtxs_d22_both->ProjectionX();
  TH1D* th1d_os_fvtxs_d24_out_both = tp1f_os_fvtxs_d24_out_both->ProjectionX();
  th1d_os_fvtxs_d22_both->Scale(2*os_fvtxs_c22);
  th1d_os_fvtxs_d22_both->Add(th1d_os_fvtxs_d24_out_both,-1.0);
  th1d_os_fvtxs_d22_both->Draw();
  c1->Print("cumulant4pc_fvtxs_pt.png");

  // ---

  delete c1;

}
