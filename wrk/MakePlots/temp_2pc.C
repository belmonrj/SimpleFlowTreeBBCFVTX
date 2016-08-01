void doit(int);

void temp_2pc()
{

  //doit(456652);
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

  TProfile* tp1f_bbcs_d22_both = (TProfile*)file->Get("bbcs_d22_both");
  TProfile* tp1f_bbcs_cos22_both = (TProfile*)file->Get("bbcs_cos22_both");
  TProfile* tp1f_bbcs_sin22_both = (TProfile*)file->Get("bbcs_sin22_both");

  TProfile* tp1f_bbcs_c22 = (TProfile*)file->Get("bbcs_c22");
  TProfile* tp1f_bbcs_cos22 = (TProfile*)file->Get("bbcs_cos22");
  TProfile* tp1f_bbcs_sin22 = (TProfile*)file->Get("bbcs_sin22");

  double bbcs_c22_raw = tp1f_bbcs_c22->GetBinContent(1);
  double bbcs_cos22_raw = tp1f_bbcs_cos22->GetBinContent(1);
  double bbcs_sin22_raw = tp1f_bbcs_sin22->GetBinContent(1);

  // --- see equation C1 in 1010.0233
  double bbcs_c22_corr = bbcs_c22_raw - pow(bbcs_cos22_raw,2.0) - pow(bbcs_sin22_raw,2.0);
  if ( bbcs_c22_corr < 0 ) cout << "YOU'RE GONNA DIE" << endl;
  cout << bbcs_c22_corr << " " << bbcs_c22_raw << " " << bbcs_cos22_raw << " " << bbcs_sin22_raw << " " << pow(bbcs_cos22_raw,2.0) << " " << pow(bbcs_sin22_raw,2.0) << endl;
  double bbcs_v22_corr = sqrt(fabs(bbcs_c22_corr));
  cout << bbcs_v22_corr << endl;

  TH1D* th1d_bbcs_d22_both = (TH1D*)tp1f_bbcs_d22_both->ProjectionX();
  TH1D* th1d_bbcs_cos22_both = (TH1D*)tp1f_bbcs_cos22_both->ProjectionX();
  TH1D* th1d_bbcs_sin22_both = (TH1D*)tp1f_bbcs_sin22_both->ProjectionX();

  th1d_bbcs_d22_both->Draw();
  c1->Print(Form("FigsTwo/raw_bbcs_d22_%d.png",handle));
  c1->Print(Form("FigsTwo/raw_bbcs_d22_%d.pdf",handle));

  th1d_bbcs_cos22_both->Draw();
  c1->Print(Form("FigsTwo/raw_bbcs_cos22_%d.png",handle));
  c1->Print(Form("FigsTwo/raw_bbcs_cos22_%d.pdf",handle));

  th1d_bbcs_sin22_both->Draw();
  c1->Print(Form("FigsTwo/raw_bbcs_sin22_%d.png",handle));
  c1->Print(Form("FigsTwo/raw_bbcs_sin22_%d.pdf",handle));

  TH1D* th1d_bbcs_d22_both_corr = (TH1D*)th1d_bbcs_d22_both->Clone();
  int nbinsx = th1d_bbcs_d22_both_corr->GetNbinsX();
  for ( int i = 0; i < nbinsx; ++i )
    {
      double pt = th1d_bbcs_d22_both->GetBinCenter(i+1);
      double d22 = th1d_bbcs_d22_both->GetBinContent(i+1);
      double cos22 = th1d_bbcs_cos22_both->GetBinContent(i+1);
      double sin22 = th1d_bbcs_sin22_both->GetBinContent(i+1);
      // --- see equation C13 in 1010.0233
      double answer = d22 - cos22*bbcs_cos22_raw - sin22*bbcs_sin22_raw;
      th1d_bbcs_d22_both_corr->SetBinContent(i+1,answer);
      if ( pt >= 1.0 && pt <= 1.5 )
        {
          cout << "pt " << pt << " "
               << "corrected d22 " << answer << " "
               << "d22 " << d22 << " "
               << "cntcos " << cos22 << " "
               << "cntsin " << sin22 << " "
               << "bbcscos " << bbcs_cos22_raw << " "
               << "bbcssin " << bbcs_sin22_raw << " "
               << "coscos " << cos22*bbcs_cos22_raw << " "
               << "sinsin " << sin22*bbcs_sin22_raw << " "
               << endl;
        }
    }

  th1d_bbcs_d22_both_corr->Draw();
  c1->Print(Form("FigsTwo/corr_bbcs_d22_%d.png",handle));
  c1->Print(Form("FigsTwo/corr_bbcs_d22_%d.pdf",handle));

  th1d_bbcs_d22_both_corr->Scale(1.0/bbcs_v22_corr);
  th1d_bbcs_d22_both_corr->Draw();
  c1->Print(Form("FigsTwo/corr_bbcs_v22_%d.png",handle));
  c1->Print(Form("FigsTwo/corr_bbcs_v22_%d.pdf",handle));

  // ---
  // ---
  // ---

  TProfile* tp1f_fvtxs_d22_both = (TProfile*)file->Get("fvtxs_d22_both");
  TProfile* tp1f_fvtxs_cos22_both = (TProfile*)file->Get("fvtxs_cos22_both");
  TProfile* tp1f_fvtxs_sin22_both = (TProfile*)file->Get("fvtxs_sin22_both");

  TProfile* tp1f_fvtxs_c22 = (TProfile*)file->Get("fvtxs_c22");
  TProfile* tp1f_fvtxs_cos22 = (TProfile*)file->Get("fvtxs_cos22");
  TProfile* tp1f_fvtxs_sin22 = (TProfile*)file->Get("fvtxs_sin22");

  double fvtxs_c22_raw = tp1f_fvtxs_c22->GetBinContent(1);
  double fvtxs_cos22_raw = tp1f_fvtxs_cos22->GetBinContent(1);
  double fvtxs_sin22_raw = tp1f_fvtxs_sin22->GetBinContent(1);

  // --- see equation C1 in 1010.0233
  double fvtxs_c22_corr = fvtxs_c22_raw - pow(fvtxs_cos22_raw,2.0) - pow(fvtxs_sin22_raw,2.0);
  if ( fvtxs_c22_corr < 0 ) cout << "YOU'RE GONNA DIE" << endl;
  cout << fvtxs_c22_corr << " " << fvtxs_c22_raw << " " << fvtxs_cos22_raw << " " << fvtxs_sin22_raw << " " << pow(fvtxs_cos22_raw,2.0) << " " << pow(fvtxs_sin22_raw,2.0) << endl;
  double fvtxs_v22_corr = sqrt(fabs(fvtxs_c22_corr));
  cout << fvtxs_v22_corr << endl;

  TH1D* th1d_fvtxs_d22_both = (TH1D*)tp1f_fvtxs_d22_both->ProjectionX();
  TH1D* th1d_fvtxs_cos22_both = (TH1D*)tp1f_fvtxs_cos22_both->ProjectionX();
  TH1D* th1d_fvtxs_sin22_both = (TH1D*)tp1f_fvtxs_sin22_both->ProjectionX();

  th1d_fvtxs_d22_both->Draw();
  c1->Print(Form("FigsTwo/raw_fvtxs_d22_%d.png",handle));
  c1->Print(Form("FigsTwo/raw_fvtxs_d22_%d.pdf",handle));

  th1d_fvtxs_cos22_both->Draw();
  c1->Print(Form("FigsTwo/raw_fvtxs_cos22_%d.png",handle));
  c1->Print(Form("FigsTwo/raw_fvtxs_cos22_%d.pdf",handle));

  th1d_fvtxs_sin22_both->Draw();
  c1->Print(Form("FigsTwo/raw_fvtxs_sin22_%d.png",handle));
  c1->Print(Form("FigsTwo/raw_fvtxs_sin22_%d.pdf",handle));

  TH1D* th1d_fvtxs_d22_both_corr = (TH1D*)th1d_fvtxs_d22_both->Clone();
  int nbinsx = th1d_fvtxs_d22_both_corr->GetNbinsX();
  for ( int i = 0; i < nbinsx; ++i )
    {
      double pt = th1d_fvtxs_d22_both->GetBinCenter(i+1);
      double d22 = th1d_fvtxs_d22_both->GetBinContent(i+1);
      double cos22 = th1d_fvtxs_cos22_both->GetBinContent(i+1);
      double sin22 = th1d_fvtxs_sin22_both->GetBinContent(i+1);
      // --- see equation C13 in 1010.0233
      double answer = d22 - cos22*fvtxs_cos22_raw - sin22*fvtxs_sin22_raw;
      th1d_fvtxs_d22_both_corr->SetBinContent(i+1,answer);
      if ( pt >= 1.0 && pt <= 1.5 )
        {
          cout << "pt " << pt << " "
               << "corrected d22 " << answer << " "
               << "d22 " << d22 << " "
               << "cntcos " << cos22 << " "
               << "cntsin " << sin22 << " "
               << "fvtxscos " << fvtxs_cos22_raw << " "
               << "fvtxssin " << fvtxs_sin22_raw << " "
               << "coscos " << cos22*fvtxs_cos22_raw << " "
               << "sinsin " << sin22*fvtxs_sin22_raw << " "
               << endl;
        }
    }

  th1d_fvtxs_d22_both_corr->Draw();
  c1->Print(Form("FigsTwo/corr_fvtxs_d22_%d.png",handle));
  c1->Print(Form("FigsTwo/corr_fvtxs_d22_%d.pdf",handle));

  th1d_fvtxs_d22_both_corr->Scale(1.0/fvtxs_v22_corr);
  th1d_fvtxs_d22_both_corr->Draw();
  c1->Print(Form("FigsTwo/corr_fvtxs_v22_%d.png",handle));
  c1->Print(Form("FigsTwo/corr_fvtxs_v22_%d.pdf",handle));





}
