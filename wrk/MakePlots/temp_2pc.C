void doit(int);

void temp_2pc()
{

  doit(456652);

}

void doit(int handle)
{

  gStyle->SetOptTitle(1);

  TCanvas* c1 = new TCanvas("c1","");

  TFile* file = TFile::Open(Form("input/hist_%d.root",handle));

  TProfile* tp1f_bbcs_d22_both = (TProfile*)file->Get("bbcs_d22_both");
  TProfile* tp1f_bbcs_cos22_both = (TProfile*)file->Get("bbcs_cos22_both");
  TProfile* tp1f_bbcs_sin22_both = (TProfile*)file->Get("bbcs_sin22_both");

  TProfile* tp1f_bbcs_c22 = (TProfile*)file->Get("bbcs_c22");
  TProfile* tp1f_bbcs_cos22 = (TProfile*)file->Get("bbcs_cos22");
  TProfile* tp1f_bbcs_sin22 = (TProfile*)file->Get("bbcs_sin22");

  double bbcs_c22_raw = tp1f_bbcs_c22->GetBinContent(1);
  double bbcs_cos22_raw = tp1f_bbcs_cos22->GetBinContent(1);
  double bbcs_sin22_raw = tp1f_bbcs_sin22->GetBinContent(1);

  double bbcs_c22_corr = bbcs_c22_raw - pow(bbcs_cos22_raw,2.0) - pow(bbcs_sin22_raw,2.0);
  if ( bbcs_c22_corr < 0 ) cout << "YOU'RE GONNA DIE" << endl;
  cout << bbcs_c22_corr << " " << bbcs_c22_raw << " " << bbcs_cos22_raw << " " << bbcs_sin22_raw << endl;

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

      double d22 = th1d_bbcs_d22_both->GetBinContent(i+1);
      double cos22 = th1d_bbcs_cos22_both->GetBinContent(i+1);
      double sin22 = th1d_bbcs_sin22_both->GetBinContent(i+1);
      double answer = d22 - pow(cos22,2.0) - pow(sin22,2.0);
      th1d_bbcs_d22_both_corr->SetBinContent(i+1,answer);
    }

  th1d_bbcs_d22_both_corr->Draw();
  c1->Print(Form("FigsTwo/corr_bbcs_d22_%d.png",handle));
  c1->Print(Form("FigsTwo/corr_bbcs_d22_%d.pdf",handle));

  th1d_bbcs_d22_both_corr->Scale(1.0/sqrt(fabs(bbcs_c22_corr)));
  th1d_bbcs_d22_both_corr->Draw();
  c1->Print(Form("FigsTwo/corr_bbcs_v22_%d.png",handle));
  c1->Print(Form("FigsTwo/corr_bbcs_v22_%d.pdf",handle));


}
