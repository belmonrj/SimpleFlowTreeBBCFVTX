void takehistograms(TProfile*,TProfile*,TProfile*,TProfile*,TProfile*,TProfile*,TProfile*,TProfile*);

double calc_corr_four(double,double,double,double,double,double,double,double);

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
  double uncorr_c2[nbins];
  double uncorr_c4[nbins];
  double corr_four[nbins];
  double corr_222[nbins];
  double uncorr_four[nbins];
  double uncorr_222[nbins];

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
      corr_c4[i] = calc_corr_four(four[i],two[i],cos1[i],sin1[i],cossum2[i],sinsum2[i],cos3[i],sin3[i]);
      corr_222[i] = 2*corr_c2[i]*corr_c2[i];
      corr_four[i] = corr_c4[i] + corr_222[i];
      uncorr_four[i] = four[i];
      uncorr_222[i] = 2*two[i]*two[i];
    }

  TGraphErrors* tge_corr_c24 = new TGraphErrors(nbins,x,corr_c4,0,0);
  tge_corr_c24->SetMarkerStyle(kOpenCircle);
  tge_corr_c24->SetMarkerColor(kBlack);
  tge_corr_c24->Draw("ap");
  tge_corr_c24->GetXaxis()->SetLimits(0,80);
  tge_corr_c24->GetXaxis()->SetTitle("N^{FVTX}_{tracks}");
  c1->Print("testcorr4pc.png");
  tge_corr_c24->SetMaximum(1e-4);
  tge_corr_c24->SetMinimum(-1e-4);
  c1->Print("testcorr4pcfz.png");

  TGraphErrors* tge_corr_c22 = new TGraphErrors(nbins,x,corr_c2,0,0);
  tge_corr_c22->SetMarkerStyle(kOpenCircle);
  tge_corr_c22->SetMarkerColor(kBlack);
  tge_corr_c22->Draw("ap");
  tge_corr_c22->GetXaxis()->SetLimits(0,80);
  tge_corr_c22->GetXaxis()->SetTitle("N^{FVTX}_{tracks}");
  c1->Print("testcorr2pc.png");
  tge_corr_c22->SetMaximum(2e-2);
  tge_corr_c22->SetMinimum(0);
  c1->Print("testcorr2pcfz.png");

  TGraphErrors* tge_corr_222 = new TGraphErrors(nbins,x,corr_222,0,0);
  tge_corr_222->SetMarkerStyle(kOpenCircle);
  tge_corr_222->SetMarkerColor(kBlack);
  tge_corr_222->Draw("ap");
  tge_corr_222->GetXaxis()->SetLimits(0,80);
  tge_corr_222->GetXaxis()->SetTitle("N^{FVTX}_{tracks}");
  TGraphErrors* tge_corr_four = new TGraphErrors(nbins,x,corr_four,0,0);
  tge_corr_four->SetMarkerStyle(kOpenSquare);
  tge_corr_four->SetMarkerColor(kRed);
  tge_corr_four->Draw("p");
  TLegend* leg = new TLegend(0.62,0.68,0.88,0.88);
  leg->AddEntry(tge_corr_222,"2#LT#LT2#GT#GT^{2} - a.c","p");
  leg->AddEntry(tge_corr_four,"#LT#LT4#GT#GT - a.c.","p");
  leg->SetTextSize(0.05);
  leg->Draw();
  c1->Print("testcorrcomponents.png");
  tge_corr_222->SetMaximum(1e-3);
  tge_corr_222->SetMinimum(-1e-4);
  c1->Print("testcorrcomponentsfz.png");
  tge_corr_222->SetMaximum(2e-4);
  tge_corr_222->SetMinimum(-2e-5);
  c1->Print("testcorrcomponentsfzz.png");

  TGraphErrors* tge_uncorr_222 = new TGraphErrors(nbins,x,uncorr_222,0,0);
  tge_uncorr_222->SetMarkerStyle(kOpenCircle);
  tge_uncorr_222->SetMarkerColor(kBlack);
  tge_uncorr_222->Draw("ap");
  tge_uncorr_222->GetXaxis()->SetLimits(0,80);
  tge_uncorr_222->GetXaxis()->SetTitle("N^{FVTX}_{tracks}");
  TGraphErrors* tge_uncorr_four = new TGraphErrors(nbins,x,uncorr_four,0,0);
  tge_uncorr_four->SetMarkerStyle(kOpenSquare);
  tge_uncorr_four->SetMarkerColor(kRed);
  tge_uncorr_four->Draw("p");
  TLegend* leg = new TLegend(0.62,0.68,0.88,0.88);
  leg->AddEntry(tge_uncorr_222,"2#LT#LT2#GT#GT^{2}","p");
  leg->AddEntry(tge_uncorr_four,"#LT#LT4#GT#GT","p");
  leg->SetTextSize(0.05);
  leg->Draw();
  c1->Print("testuncorrcomponents.png");
  tge_uncorr_222->SetMaximum(1e-3);
  tge_uncorr_222->SetMinimum(-1e-4);
  c1->Print("testuncorrcomponentsfz.png");
  tge_uncorr_222->SetMaximum(2e-4);
  tge_uncorr_222->SetMinimum(-2e-5);
  c1->Print("testuncorrcomponentsfzz.png");

}

double calc_corr_four(double four, double two, double cos1, double sin1, double cossum2, double sinsum2, double cos3, double sin3)
{
  double uncorr = four - 2*two*two;
  double corr_term1 = 4*cos1*cos3;
  double corr_term2 = 4*sin1*sin3;
  double corr_term3 = cossum2*cossum2;
  double corr_term4 = sinsum2*sinsum2;
  double corr_term5 = 4*cossum2*(cos1*cos1 - sin1*sin1);
  double corr_term6 = 8*sinsum2*sin1*cos1;
  double corr_term7 = 8*two*(cos1*cos1 + sin1*sin1);
  double corr_term8 = 6*(cos1*cos1 + sin1*sin1)*(cos1*cos1 + sin1*sin1);
  double result = uncorr - corr_term1 + corr_term2 - corr_term3 - corr_term4 + corr_term5 + corr_term6 + corr_term7 - corr_term8;
  return result;
}
