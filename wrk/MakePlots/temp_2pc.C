void doit(int);

void temp_2pc()
{

  doit(456652);
  // doit(200);
  // doit(62);
  // doit(39);
  // doit(20);

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

  // -------------
  // --- FVTXN ---
  // -------------

  TProfile* tp1f_os_fvtxn_d22_both = (TProfile*)file->Get("os_fvtxn_d22_both");
  TProfile* tp1f_os_fvtxn_cos22_both = (TProfile*)file->Get("os_fvtxn_cos22_both");
  TProfile* tp1f_os_fvtxn_sin22_both = (TProfile*)file->Get("os_fvtxn_sin22_both");

  TProfile* tp1f_os_fvtxn_c22 = (TProfile*)file->Get("os_fvtxn_c22");
  TProfile* tp1f_os_fvtxn_cos22 = (TProfile*)file->Get("os_fvtxn_cos22");
  TProfile* tp1f_os_fvtxn_sin22 = (TProfile*)file->Get("os_fvtxn_sin22");

  double os_fvtxn_c22_raw = tp1f_os_fvtxn_c22->GetBinContent(1);
  double os_fvtxn_cos22_raw = tp1f_os_fvtxn_cos22->GetBinContent(1);
  double os_fvtxn_sin22_raw = tp1f_os_fvtxn_sin22->GetBinContent(1);

  double os_fvtxn_c22_corr = os_fvtxn_c22_raw;
  if ( os_fvtxn_c22_corr < 0 ) cout << "YOU'RE GONNA DIE" << endl;
  double os_fvtxn_v22_corr = sqrt(fabs(os_fvtxn_c22_corr));
  cout << os_fvtxn_v22_corr << endl;

  TH1D* th1d_os_fvtxn_d22_both = (TH1D*)tp1f_os_fvtxn_d22_both->ProjectionX();
  TH1D* th1d_os_fvtxn_d22_both_corr = (TH1D*)th1d_os_fvtxn_d22_both->Clone();
  th1d_os_fvtxn_d22_both_corr->Draw();
  c1->Print(Form("FigsTwo/corr_os_fvtxn_d22_%d.png",handle));
  c1->Print(Form("FigsTwo/corr_os_fvtxn_d22_%d.pdf",handle));

  th1d_os_fvtxn_d22_both_corr->Scale(1.0/os_fvtxn_v22_corr);
  if ( handle <= 200 ) th1d_os_fvtxn_d22_both_corr->SetTitle(Form("d+Au #sqrt{s_{NN}} = %d GeV",handle));
  th1d_os_fvtxn_d22_both_corr->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  th1d_os_fvtxn_d22_both_corr->GetYaxis()->SetTitle("v_{2}{2}");
  th1d_os_fvtxn_d22_both_corr->Draw();
  c1->Print(Form("FigsTwo/corr_os_fvtxn_v22_%d.png",handle));
  c1->Print(Form("FigsTwo/corr_os_fvtxn_v22_%d.pdf",handle));

  // ---

  // -------------------------------------------------
  // --- different method of estimating reference flow
  // -------------------------------------------------

  // --- note that the bbcs-fvtxn correlation is so weak in the lower energies that this method might not work

  TProfile* tp1f_os_bbcsfvtxs_c22 = (TProfile*)file->Get("os_bbcsfvtxs_c22");
  TProfile* tp1f_os_bbcsfvtxn_c22 = (TProfile*)file->Get("os_bbcsfvtxn_c22");
  TProfile* tp1f_os_fvtxsfvtxn_c22 = (TProfile*)file->Get("os_fvtxsfvtxn_c22");

  double os_bbcsfvtxs_c22_raw = tp1f_os_bbcsfvtxs_c22->GetBinContent(1);
  double os_bbcsfvtxn_c22_raw = tp1f_os_bbcsfvtxn_c22->GetBinContent(1);
  double os_fvtxsfvtxn_c22_raw = tp1f_os_fvtxsfvtxn_c22->GetBinContent(1);

  double os_newref_c22_bbcs = (os_bbcsfvtxs_c22_raw*os_bbcsfvtxn_c22_raw)/os_fvtxsfvtxn_c22_raw;
  double os_newref_c22_fvtxs = (os_bbcsfvtxs_c22_raw*os_fvtxsfvtxn_c22_raw)/os_bbcsfvtxn_c22_raw;
  double os_newref_c22_fvtxn = (os_bbcsfvtxn_c22_raw*os_fvtxsfvtxn_c22_raw)/os_bbcsfvtxs_c22_raw;

  if ( os_newref_c22_bbcs < 0 ) cout << "YOU'RE GONNA DIE" << endl;
  if ( os_newref_c22_fvtxs < 0 ) cout << "YOU'RE GONNA DIE" << endl;
  if ( os_newref_c22_fvtxn < 0 ) cout << "YOU'RE GONNA DIE" << endl;

  double os_newref_v22_bbcs =  sqrt(fabs(os_newref_c22_bbcs));
  double os_newref_v22_fvtxs = sqrt(fabs(os_newref_c22_fvtxs));
  double os_newref_v22_fvtxn = sqrt(fabs(os_newref_c22_fvtxn));

  cout << os_newref_v22_bbcs << endl;
  cout << os_newref_v22_fvtxs << endl;
  cout << os_newref_v22_fvtxn << endl;

  // --- rescale to remove the old v2 and apply the new v2
  th1d_os_bbcs_d22_both_corr->Scale(os_bbcs_v22_corr/os_newref_v22_bbcs);
  th1d_os_fvtxs_d22_both_corr->Scale(os_fvtxs_v22_corr/os_newref_v22_fvtxs);
  th1d_os_fvtxn_d22_both_corr->Scale(os_fvtxn_v22_corr/os_newref_v22_fvtxn);

  th1d_os_bbcs_d22_both_corr->Draw();
  c1->Print(Form("FigsTwo/newref_os_bbcs_v22_%d.png",handle));
  c1->Print(Form("FigsTwo/newref_os_bbcs_v22_%d.pdf",handle));

  th1d_os_fvtxs_d22_both_corr->Draw();
  c1->Print(Form("FigsTwo/newref_os_fvtxs_v22_%d.png",handle));
  c1->Print(Form("FigsTwo/newref_os_fvtxs_v22_%d.pdf",handle));

  th1d_os_fvtxn_d22_both_corr->Draw();
  c1->Print(Form("FigsTwo/newref_os_fvtxn_v22_%d.png",handle));
  c1->Print(Form("FigsTwo/newref_os_fvtxn_v22_%d.pdf",handle));

  // ---

  // --------------------------------------------------------------------------------------------
  // --- 3 component differential flow scalar product (makes me contemplate 3 particle cumulants)
  // --------------------------------------------------------------------------------------------

  TH1D* th1d_os_bbcsfvtxs_v22_3csp = (TH1D*)th1d_os_bbcs_d22_both_corr->Clone();
  TH1D* th1d_os_bbcsfvtxn_v22_3csp = (TH1D*)th1d_os_bbcs_d22_both_corr->Clone();
  TH1D* th1d_os_fvtxsfvtxn_v22_3csp = (TH1D*)th1d_os_bbcs_d22_both_corr->Clone();

  int nbins = th1d_os_bbcsfvtxs_v22_3csp->GetNbinsX();
  for ( int i = 0; i < nbins; ++i )
    {
      double d22_bbcs = th1d_os_bbcs_d22_both->GetBinContent(i+1);
      double d22_fvtxs = th1d_os_fvtxs_d22_both->GetBinContent(i+1);
      double d22_fvtxn = th1d_os_fvtxn_d22_both->GetBinContent(i+1);

      double answer_bbcsfvtxs = (d22_bbcs*d22_fvtxs)/os_bbcsfvtxs_c22_raw;
      double answer_bbcsfvtxn = (d22_bbcs*d22_fvtxn)/os_bbcsfvtxn_c22_raw;
      double answer_fvtxsfvtxn = (d22_fvtxs*d22_fvtxs)/os_fvtxsfvtxn_c22_raw;

      if ( answer_bbcsfvtxs != answer_bbcsfvtxs ) answer_bbcsfvtxs = 0;
      if ( answer_bbcsfvtxn != answer_bbcsfvtxn ) answer_bbcsfvtxn = 0;
      if ( answer_fvtxsfvtxn != answer_fvtxsfvtxn ) answer_fvtxsfvtxn = 0;

      th1d_os_bbcsfvtxs_v22_3csp->SetBinContent(i+1,answer_bbcsfvtxs);
      th1d_os_bbcsfvtxn_v22_3csp->SetBinContent(i+1,answer_bbcsfvtxn);
      th1d_os_fvtxsfvtxn_v22_3csp->SetBinContent(i+1,answer_fvtxsfvtxn);
    }

  th1d_os_bbcsfvtxs_v22_3csp->Draw();
  c1->Print(Form("FigsTwo/threesp_os_bbcsfvtxs_v22_%d.png",handle));
  c1->Print(Form("FigsTwo/threesp_os_bbcsfvtxs_v22_%d.pdf",handle));

  th1d_os_bbcsfvtxn_v22_3csp->Draw();
  c1->Print(Form("FigsTwo/threesp_os_bbcsfvtxn_v22_%d.png",handle));
  c1->Print(Form("FigsTwo/threesp_os_bbcsfvtxn_v22_%d.pdf",handle));

  th1d_os_fvtxsfvtxn_v22_3csp->Draw();
  c1->Print(Form("FigsTwo/threesp_os_fvtxsfvtxn_v22_%d.png",handle));
  c1->Print(Form("FigsTwo/threesp_os_fvtxsfvtxn_v22_%d.pdf",handle));

  // ---

  delete c1;

}
