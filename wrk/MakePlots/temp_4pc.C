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
  // --- now let's have a look at the c24 vs nfvtxt
  // ---

  TProfile* tp1f_nfvtxt_fvtxns_c22 = (TProfile*)file->Get("nfvtxt_os_fvtxsfvtxn_c22");
  TProfile* tp1f_nfvtxt_fvtxns_c24a = (TProfile*)file->Get("nfvtxt_os_fvtxsfvtxn_c24a");
  TProfile* tp1f_nfvtxt_fvtxns_c24b = (TProfile*)file->Get("nfvtxt_os_fvtxsfvtxn_c24b");
  TProfile* tp1f_nfvtxt_fvtxns_c24c = (TProfile*)file->Get("nfvtxt_os_fvtxsfvtxn_c24c");
  TH1D* th1d_nfvtxt_fvtxns_c24a = tp1f_nfvtxt_fvtxns_c24a->ProjectionX();
  TH1D* th1d_nfvtxt_fvtxns_c24b = tp1f_nfvtxt_fvtxns_c24b->ProjectionX();
  TH1D* th1d_nfvtxt_fvtxns_c24c = tp1f_nfvtxt_fvtxns_c24c->ProjectionX();

  TH1D* th1d_nfvtxt_fvtxns_c22 = tp1f_nfvtxt_fvtxns_c22->ProjectionX();
  TH1D* th1d_nfvtxt_fvtxns_c22a = (TH1D*)th1d_nfvtxt_fvtxns_c22->Clone("th1d_nfvtxt_fvtxns_c22a");
  TH1D* th1d_nfvtxt_fvtxns_c22b = (TH1D*)th1d_nfvtxt_fvtxns_c22->Clone("th1d_nfvtxt_fvtxns_c22b");
  TH1D* th1d_nfvtxt_fvtxns_c22c = (TH1D*)th1d_nfvtxt_fvtxns_c22->Clone("th1d_nfvtxt_fvtxns_c22c");

  th1d_nfvtxt_fvtxns_c22a->Multiply(th1d_nfvtxt_fvtxns_c22);
  th1d_nfvtxt_fvtxns_c22a->Scale(-2.0);
  th1d_nfvtxt_fvtxns_c22a->Add(th1d_nfvtxt_fvtxns_c24a,1.0);
  th1d_nfvtxt_fvtxns_c22a->Draw();
  TLine line(0.0,0.0,75.0,0.0);
  line.SetLineWidth(2);
  line.SetLineStyle(2);
  line.Draw();
  c1->Print("c24a_nfvtxt.png");

  th1d_nfvtxt_fvtxns_c22b->Multiply(th1d_nfvtxt_fvtxns_c22);
  th1d_nfvtxt_fvtxns_c22b->Scale(-2.0);
  th1d_nfvtxt_fvtxns_c22b->Add(th1d_nfvtxt_fvtxns_c24b,1.0);
  th1d_nfvtxt_fvtxns_c22b->SetMaximum(1e-5);
  th1d_nfvtxt_fvtxns_c22b->SetMinimum(-1e-5);
  th1d_nfvtxt_fvtxns_c22b->Draw();
  line.Draw();
  c1->Print("c24b_nfvtxt.png");

  th1d_nfvtxt_fvtxns_c22c->Multiply(th1d_nfvtxt_fvtxns_c22);
  th1d_nfvtxt_fvtxns_c22c->Scale(-2.0);
  th1d_nfvtxt_fvtxns_c22c->Add(th1d_nfvtxt_fvtxns_c24c,1.0);
  th1d_nfvtxt_fvtxns_c22c->SetMaximum(1e-5);
  th1d_nfvtxt_fvtxns_c22c->SetMinimum(-1e-5);
  th1d_nfvtxt_fvtxns_c22c->Draw();
  line.Draw();
  c1->Print("c24c_nfvtxt.png");

  const int nbins = tp1f_nfvtxt_fvtxns_c22->GetNbinsX();
  TH1D* th1d_nfvtxt_fvtxns_v24a = tp1f_nfvtxt_fvtxns_c24a->ProjectionX("th1d_nfvtxt_fvtxns_v24a");
  TH1D* th1d_nfvtxt_fvtxns_v24b = tp1f_nfvtxt_fvtxns_c24b->ProjectionX("th1d_nfvtxt_fvtxns_v24b");
  TH1D* th1d_nfvtxt_fvtxns_v24c = tp1f_nfvtxt_fvtxns_c24c->ProjectionX("th1d_nfvtxt_fvtxns_v24c");
  for ( int i = 0; i < nbins; ++i )
    {
      // --- v24a
      float c24 = th1d_nfvtxt_fvtxns_c22a->GetBinContent(i+1);
      float ec24 = th1d_nfvtxt_fvtxns_c22a->GetBinError(i+1);
      float v24 = 0;
      if ( c24 < 0 ) v24 = sqrt(sqrt(-c24));
      float ev24 = fabs(ec24*v24/c24); // not right, get back to it soon
      th1d_nfvtxt_fvtxns_v24a->SetBinContent(i+1,v24);
      th1d_nfvtxt_fvtxns_v24a->SetBinError(i+1,ev24);
      // --- v24b
      c24 = th1d_nfvtxt_fvtxns_c22b->GetBinContent(i+1);
      ec24 = th1d_nfvtxt_fvtxns_c22b->GetBinError(i+1);
      v24 = 0;
      if ( c24 < 0 ) v24 = sqrt(sqrt(-c24));
      ev24 = fabs(ec24*v24/c24); // not right, get back to it soon
      th1d_nfvtxt_fvtxns_v24b->SetBinContent(i+1,v24);
      th1d_nfvtxt_fvtxns_v24b->SetBinError(i+1,ev24);
      cout << c24 << " " << v24 << " " << ev24 << endl;
      // --- v24c
      c24 = th1d_nfvtxt_fvtxns_c22c->GetBinContent(i+1);
      ec24 = th1d_nfvtxt_fvtxns_c22c->GetBinError(i+1);
      v24 = 0;
      if ( c24 < 0 ) v24 = sqrt(sqrt(-c24));
      ev24 = fabs(ec24*v24/c24); // not right, get back to it soon
      th1d_nfvtxt_fvtxns_v24c->SetBinContent(i+1,v24);
      th1d_nfvtxt_fvtxns_v24c->SetBinError(i+1,ev24);
      // cout << c24 << " " << v24 << " " << ev24 << endl;
    }


  th1d_nfvtxt_fvtxns_v24a->Draw();
  c1->Print("v24a_nfvtxt.png");
  th1d_nfvtxt_fvtxns_v24b->Draw();
  c1->Print("v24b_nfvtxt.png");
  th1d_nfvtxt_fvtxns_v24c->Draw();
  c1->Print("v24c_nfvtxt.png");

  delete c1;

}
