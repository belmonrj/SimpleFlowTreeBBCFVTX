TH1D* sqrt(const TH1D*);

void go_frac_yourself()
{

  TCanvas* c1 = new TCanvas();

  TFile* fbase = TFile::Open("cumulants_pAu_base.root");
  TProfile* tp1f_c22gap_base = (TProfile*)fbase->Get("nfvtxt_ac_fvtxsfvtxn_tracks_c22");
  TProfile* tp1f_c22_base = (TProfile*)fbase->Get("nfvtxt_ac_fvtxc_tracks_c22");

  TFile* ffrac = TFile::Open("cumulants_pAu_frac.root");
  TProfile* tp1f_c22gap_frac = (TProfile*)ffrac->Get("nfvtxt_ac_fvtxsfvtxn_tracks_c22");
  TProfile* tp1f_c22_frac = (TProfile*)ffrac->Get("nfvtxt_ac_fvtxc_tracks_c22");

  TH1D* th1d_c22gap_base = tp1f_c22gap_base->ProjectionX("th1d_c22gap_base");
  TH1D* th1d_c22_base = tp1f_c22_base->ProjectionX("th1d_c22_base");
  TH1D* th1d_c22gap_frac = tp1f_c22gap_frac->ProjectionX("th1d_c22gap_frac");
  TH1D* th1d_c22_frac = tp1f_c22_frac->ProjectionX("th1d_c22_frac");

  TH1D* th1d_v22gap_base = sqrt(th1d_c22gap_base);
  TH1D* th1d_v22_base = sqrt(th1d_c22_base);
  TH1D* th1d_v22gap_frac = sqrt(th1d_c22gap_frac);
  TH1D* th1d_v22_frac = sqrt(th1d_c22_frac);

  th1d_v22gap_base->SetMarkerStyle(kFullDiamond);
  th1d_v22gap_frac->SetMarkerStyle(kOpenDiamond);
  th1d_v22_base->SetMarkerStyle(kFullCircle);
  th1d_v22_frac->SetMarkerStyle(kOpenCircle);

  th1d_v22gap_base->SetMarkerColor(kMagenta+2);
  th1d_v22gap_frac->SetMarkerColor(kMagenta+2);
  th1d_v22_base->SetMarkerColor(kRed);
  th1d_v22_frac->SetMarkerColor(kRed);

  th1d_v22gap_base->Draw("ex0p");
  th1d_v22gap_base->SetMaximum(0.18);
  th1d_v22gap_base->SetMinimum(0.0);
  th1d_v22gap_frac->Draw("ex0p same");
  th1d_v22_base->Draw("ex0p same");
  th1d_v22_frac->Draw("ex0p same");

  TLegend* leg = new TLegend(0.18,0.68,0.38,0.88);
  leg->AddEntry(th1d_v22_base,"v22, no gap, no frac cut","p");
  leg->AddEntry(th1d_v22_frac,"v22, no gap, frac cut","p");
  leg->AddEntry(th1d_v22gap_base,"v22, gap, no frac cut","p");
  leg->AddEntry(th1d_v22gap_frac,"v22, gap, frac cut","p");
  leg->SetTextSize(0.045);
  leg->Draw();

  c1->Print("what_the_frac.png");

}

TH1D* sqrt(const TH1D* hold)
{
  TH1D* hnew = (TH1D*)hold->Clone();
  int nbins = hold->GetNbinsX();
  for ( int i = 0; i < nbins; ++i )
    {
      double oldcontent = hold->GetBinContent(i+1);
      double olderror = hold->GetBinError(i+1);
      double newcontent = -9;
      double newerror = 999;
      if ( oldcontent >= 0 )
	{
	  newcontent = sqrt(oldcontent);
	  newerror = olderror/sqrt(oldcontent);
	}
      hnew->SetBinContent(i+1,newcontent);
      hnew->SetBinError(i+1,newerror);
    }
  hnew->SetName(Form("%s_sqrt",hold->GetName()));
  return hnew;
}
