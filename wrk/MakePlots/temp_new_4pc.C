void takehistograms(TProfile*,TProfile*,TProfile*,TProfile*,TProfile*,TProfile*,TProfile*,TProfile*);

void calc_corr_four(double,double,double,double,double,double,double,double);

void temp_new_4pc()
{

  TFile* fin = TFile::Open("input/cumulants_200.root");

  TProfile* tp1f_four = (TProfile*)fin->Get("nfvtxt_os_fvtxc_tracks_c24");
  TProfile* tp1f_two = (TProfile*)fin->Get("nfvtxt_os_fvtxc_tracks_c22");
  TProfile* tp1f_cos1 = (TProfile*)fin->Get("nfvtxt_os_fvtxc_tracks_cos21");
  TProfile* tp1f_sin1 = (TProfile*)fin->Get("nfvtxt_os_fvtxc_tracks_sin21");
  TProfile* tp1f_cossum2 = (TProfile*)fin->Get("nfvtxt_os_fvtxc_tracks_cossum22");
  TProfile* tp1f_sinsum2 = (TProfile*)fin->Get("nfvtxt_os_fvtxc_tracks_sinsum22");
  TProfile* tp1f_cos3 = (TProfile*)fin->Get("nfvtxt_os_fvtxc_tracks_cos23");
  TProfile* tp1f_sin3 = (TProfile*)fin->Get("nfvtxt_os_fvtxc_tracks_sin23");

  takehistograms(tp1f_four,tp1f_two,tp1f_cos1,tp1f_sin1,tp1f_cossum2,tp1f_sinsum2,tp1f_cos3,tp1f_sin3);

}

void takehistograms
(
 TProfile* tp1f_four,
 TProfile* tp1f_two,
 TProfile* tp1f_cos1,
 TProfile* tp1f_sin1,
 TProfile* tp1f_cossum2,
 TProfile* tp1f_sinsum2,
 TProfile* tp1f_cos3,
 TProfile* tp1f_sin3
)
{

  int nbinsx = tp1f_four->GetNbinsX();
  const int nbins = nbinsx;

  double four[nbins];
  double two[nbins];
  double cos1[nbins];
  double sin1[nbins];
  double cossum2[nbins];
  double sinsum2[nbins];
  double cos3[nbins];
  double sin3[nbins];

  double x[nbins]; // generic x-axis, multiplicity in our case for now

  double corr_c2[nbins];
  double corr_c4[nbins];
  double uncorr_c4[nbins];

  for ( int i = 0; i < nbins; ++i )
    {
      // --- xaxis
      x[i] = tp1f_four->GetBinCenter(i+1);
      // --- components
      four[i]    = tp1f_four->GetBinContent(i+1);
      two[i]     = tp1f_two->GetBinContent(i+1);
      cos1[i]    = tp1f_cos1->GetBinContent(i+1);
      sin1[i]    = tp1f_sin1->GetBinContent(i+1);
      cossum2[i] = tp1f_cossum2->GetBinContent(i+1);
      sinsum2[i] = tp1f_sinsum2->GetBinContent(i+1);
      cos3[i]    = tp1f_cos3->GetBinContent(i+1);
      sin3[i]    = tp1f_sin3->GetBinContent(i+1);
      // --- some corrections...
      uncorr_c4[i] = four[i] - 2*two[i]*two[i]; // not very useful but have for completeness
      corr_c2[i] = two[i] - cos1[i]*cos1[i] - sin1[i]*sin1[i]; // simple enough it doesn't need it's own function
      corr_c4[i] = calc_corr_four(); // come back to this asap
    }

}
