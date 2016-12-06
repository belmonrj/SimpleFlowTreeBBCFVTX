void takehistograms(TProfile*,TProfile*,TProfile*,TProfile*,TProfile*,TProfile*,TProfile*,TProfile*,TProfile*,TProfile*,TProfile*);

double calc_corr_four(double,double,double,double,double,double,double,double);

void temp_new_4pc()
{

  TFile* fin = TFile::Open("input/cumulants_200.root");

  TProfile* tp1f_four = (TProfile*)fin->Get("nfvtxt_os_fvtxc_tracks_c24");
  TProfile* tp1f_two = (TProfile*)fin->Get("nfvtxt_os_fvtxc_tracks_c22");
  TProfile* tp1f_S_four = (TProfile*)fin->Get("nfvtxt_zzyzx_fvtxc_tracks_c24");
  TProfile* tp1f_S_two = (TProfile*)fin->Get("nfvtxt_zzyzx_fvtxc_tracks_c22");
  TProfile* tp1f_cos1 = (TProfile*)fin->Get("nfvtxt_os_fvtxc_tracks_cos21");
  TProfile* tp1f_sin1 = (TProfile*)fin->Get("nfvtxt_os_fvtxc_tracks_sin21");
  TProfile* tp1f_cossum2 = (TProfile*)fin->Get("nfvtxt_os_fvtxc_tracks_cossum22");
  TProfile* tp1f_sinsum2 = (TProfile*)fin->Get("nfvtxt_os_fvtxc_tracks_sinsum22");
  TProfile* tp1f_cos3 = (TProfile*)fin->Get("nfvtxt_os_fvtxc_tracks_cos23");
  TProfile* tp1f_sin3 = (TProfile*)fin->Get("nfvtxt_os_fvtxc_tracks_sin23");
  TProfile* tp1f_SS_two = (TProfile*)fin->Get("nfvtxt_zzyzx_fvtxsfvtxn_tracks_c22");

  takehistograms(tp1f_four,tp1f_two,tp1f_cos1,tp1f_sin1,tp1f_cossum2,tp1f_sinsum2,tp1f_cos3,tp1f_sin3,tp1f_S_four,tp1f_S_two,tp1f_SS_two);

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
 TProfile* tp1f_sin3,
 TProfile* tp1f_S_four,
 TProfile* tp1f_S_two,
 TProfile* tp1f_SS_two
)
{

  int nbinsx = tp1f_four->GetNbinsX();
  const int nbins = nbinsx;

  double four[nbins];
  double two[nbins];
  double sfour[nbins];
  double stwo[nbins];
  double cos1[nbins];
  double sin1[nbins];
  double cossum2[nbins];
  double sinsum2[nbins];
  double cos3[nbins];
  double sin3[nbins];

  double x[nbins]; // generic x-axis, multiplicity in our case for now
  double ex[nbins];

  double corr_c2[nbins];
  double corr_c4[nbins];
  double uncorr_c2[nbins];
  double uncorr_c4[nbins];
  double corr_four[nbins];
  double corr_222[nbins];
  double uncorr_four[nbins];
  double uncorr_222[nbins];

  double corr_ssc2[nbins];
  double corr_sc2[nbins];
  double corr_sc4[nbins];
  double corr_sfour[nbins];
  double corr_s222[nbins];

  double corr_v24[nbins];
  double corr_sv24[nbins];
  double corr_sv22[nbins];
  double corr_ssv22[nbins];

  TH1D* th1d_222 = new TH1D("th1d_222","",80,0,80);
  TH1D* th1d_four = new TH1D("th1d_four","",80,0,80);
  TH1D* th1d_s222 = new TH1D("th1d_s222","",80,0,80);
  TH1D* th1d_sfour = new TH1D("th1d_sfour","",80,0,80);

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
      // --- make some histograms
      th1d_222->SetBinContent(i+1,corr_222[i]);
      th1d_four->SetBinContent(i+1,corr_four[i]);
      // --- california
      sfour[i]    = tp1f_S_four->GetBinContent(i+1);
      stwo[i]     = tp1f_S_two->GetBinContent(i+1);
      corr_sfour[i] = sfour[i];
      corr_s222[i] = 2*stwo[i]*stwo[i];
      corr_sc4[i] = sfour[i] - 2*stwo[i]*stwo[i];
      corr_sc2[i] = stwo[i];
      corr_ssc2[i] = tp1f_SS_two->GetBinContent(i+1);
      // ---
      corr_v24[i] = -9;
      if ( corr_c4[i] < 0 && x[i] > 19 && x[i] < 60 ) corr_v24[i] = sqrt(sqrt(-corr_c4[i]));
      corr_sv24[i] = -9;
      if ( corr_sc4[i] < 0 && x[i] > 19 && x[i] < 60 ) corr_sv24[i] = sqrt(sqrt(-corr_sc4[i]));
      corr_sv22[i] = -9;
      if ( corr_sc2[i] > 0 && x[i] > 10 && x[i] < 70 ) corr_sv22[i] = sqrt(corr_sc2[i]);
      corr_ssv22[i] = -9;
      if ( corr_ssc2[i] > 0 && x[i] > 10 && x[i] < 70 ) corr_ssv22[i] = sqrt(corr_ssc2[i]);
    }

  TGraphErrors* tge_corr_c24 = new TGraphErrors(nbins,x,corr_c4,0,0);
  tge_corr_c24->SetMarkerStyle(kOpenCircle);
  tge_corr_c24->SetMarkerColor(kBlack);
  tge_corr_c24->Draw("ap");
  tge_corr_c24->GetXaxis()->SetLimits(0,80);
  tge_corr_c24->GetXaxis()->SetTitle("N^{FVTX}_{tracks}");
  c1->Print("FigsFour/testcorr4pc.png");
  tge_corr_c24->SetMaximum(1e-4);
  tge_corr_c24->SetMinimum(-1e-4);
  c1->Print("FigsFour/testcorr4pcfz.png");

  TGraphErrors* tge_corr_c22 = new TGraphErrors(nbins,x,corr_c2,0,0);
  tge_corr_c22->SetMarkerStyle(kOpenCircle);
  tge_corr_c22->SetMarkerColor(kBlack);
  tge_corr_c22->Draw("ap");
  tge_corr_c22->GetXaxis()->SetLimits(0,80);
  tge_corr_c22->GetXaxis()->SetTitle("N^{FVTX}_{tracks}");
  c1->Print("FigsFour/testcorr2pc.png");
  tge_corr_c22->SetMaximum(2e-2);
  tge_corr_c22->SetMinimum(0);
  c1->Print("FigsFour/testcorr2pcfz.png");

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
  leg->AddEntry(tge_corr_222,"2#LT#LT2#GT#GT^{2} - a.c.","p");
  leg->AddEntry(tge_corr_four,"#LT#LT4#GT#GT - a.c.","p");
  leg->SetTextSize(0.05);
  leg->Draw();
  c1->Print("FigsFour/testcorrcomponents.png");
  tge_corr_222->SetMaximum(1e-3);
  tge_corr_222->SetMinimum(-1e-4);
  c1->Print("FigsFour/testcorrcomponentsfz.png");
  tge_corr_222->SetMaximum(2e-4);
  tge_corr_222->SetMinimum(-2e-5);
  c1->Print("FigsFour/testcorrcomponentsfzz.png");

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
  c1->Print("FigsFour/testuncorrcomponents.png");
  tge_uncorr_222->SetMaximum(1e-3);
  tge_uncorr_222->SetMinimum(-1e-4);
  c1->Print("FigsFour/testuncorrcomponentsfz.png");
  tge_uncorr_222->SetMaximum(2e-4);
  tge_uncorr_222->SetMinimum(-2e-5);
  c1->Print("FigsFour/testuncorrcomponentsfzz.png");

  // ---

  TH1D* th1d_222_clone = (TH1D*)th1d_222->Clone();
  TH1D* th1d_four_clone = (TH1D*)th1d_four->Clone();

  TFile* file = TFile::Open(Form("input/combined_200.root"));
  TProfile* tp1f_c22 = (TProfile*)file->Get("nfvtxt_os_fvtxc_tracks_c22");
  TProfile* tp1f_c24a = (TProfile*)file->Get("nfvtxt_os_fvtxc_tracks_c24");
  TH1D* th1d_c24a = tp1f_c24a->ProjectionX(Form("th1d_c24a"));
  TH1D* th1d_c22 = tp1f_c22->ProjectionX(Form("th1d_c22"));
  TH1D* th1d_c22a = (TH1D*)th1d_c22->Clone(Form("th1d_c22a"));
  th1d_c22a->Multiply(th1d_c22);
  th1d_c22a->Scale(2.0);

  // th1d_222_clone->Divide(th1d_c22a);
  // th1d_222_clone->Draw();
  // c1->Print("FigsFour/testratio222.png");
  // th1d_four_clone->Divide(th1d_c24a);
  // th1d_four_clone->Draw();
  // c1->Print("FigsFour/testratiofour.png");
  // th1d_222_clone->Draw("same");
  // c1->Print("FigsFour/testratiocomponents.png");

  float ratio_222[nbins];
  float ratio_four[nbins];
  float eratio_222[nbins];
  float eratio_four[nbins];
  float n[nbins];
  float cumulant_ac[nbins];
  float cumulant_os[nbins];
  for ( int i = 0; i < nbins; ++i )
    {
      float value_222_a = corr_222[i];
      float value_222_b = th1d_c22a->GetBinContent(i+1);
      float value_four_a = corr_four[i];
      float value_four_b = th1d_c24a->GetBinContent(i+1);
      ratio_222[i] = -9;
      ratio_four[i] = -9;
      if ( value_222_b > 0 ) ratio_222[i] = value_222_a/value_222_b;
      if ( value_four_b > 0 ) ratio_four[i] = value_four_a/value_four_b;
      eratio_222[i] = 0;
      eratio_four[i] = 0;
      ex[i] = 0;
      n[i] = th1d_c22a->GetBinCenter(i+1);
      if ( i%5 == 0 ) cout << x[i] << " " << n[i] << " " << ratio_222[i] << endl;
      cumulant_ac[i] = corr_c4[i];
      cumulant_os[i] = value_four_b - value_222_b;
    }

  TGraphErrors* tge_cumulant_ac = new TGraphErrors(nbins,n,cumulant_ac,0,0);
  TGraphErrors* tge_cumulant_os = new TGraphErrors(nbins,n,cumulant_os,0,0);
  tge_cumulant_ac->SetMarkerStyle(kOpenCircle);
  tge_cumulant_ac->SetMarkerColor(kBlack);
  tge_cumulant_os->SetMarkerStyle(kOpenSquare);
  tge_cumulant_os->SetMarkerColor(kBlue); // all the deepest blues are black
  tge_cumulant_ac->Draw("ap");
  tge_cumulant_os->Draw("p");
  tge_cumulant_ac->SetMaximum(5e-5);
  tge_cumulant_ac->SetMinimum(-5e-5);
  tge_cumulant_ac->GetXaxis()->SetLimits(0,80);
  tge_cumulant_ac->GetXaxis()->SetTitle("N^{FVTX}_{tracks}");
  TLine line0(0,0,80,0);
  line0.SetLineStyle(2);
  line0.SetLineWidth(2);
  line0.Draw();
  c1->Print("FigsFour/testcomparecumulants.png");

  TGraphErrors* tge_ratio_222 = new TGraphErrors(nbins,n,ratio_222,0,0);
  tge_ratio_222->SetMarkerStyle(kOpenCircle);
  tge_ratio_222->SetMarkerColor(kBlack);
  tge_ratio_222->Draw("ap");
  tge_ratio_222->SetMaximum(2.0);
  tge_ratio_222->SetMinimum(0.0);
  tge_ratio_222->GetXaxis()->SetLimits(0,80);
  tge_ratio_222->GetXaxis()->SetTitle("N^{FVTX}_{tracks}");
  TGraphErrors* tge_ratio_four = new TGraphErrors(nbins,n,ratio_four,0,0);
  tge_ratio_four->SetMarkerStyle(kOpenSquare);
  tge_ratio_four->SetMarkerColor(kRed);
  tge_ratio_four->Draw("p");
  TLegend* leg = new TLegend(0.62,0.68,0.88,0.88);
  leg->AddEntry(tge_corr_222,"2#LT#LT2#GT#GT^{2} - a.c.","p");
  leg->AddEntry(tge_corr_four,"#LT#LT4#GT#GT - a.c.","p");
  leg->SetTextSize(0.05);
  leg->Draw();
  TLine line1(0,1,80,1);
  line1.SetLineStyle(2);
  line1.SetLineWidth(2);
  line1.Draw();
  TLine line2(0,1.1,80,1.1);
  line2.SetLineStyle(2);
  line2.SetLineWidth(2);
  line2.Draw();
  c1->Print("FigsFour/testcorrcomponentsratio.png");

  // --- now go back to california

  TGraphErrors* tge_corr_sc24 = new TGraphErrors(nbins,x,corr_sc4,0,0);
  tge_corr_sc24->SetMarkerStyle(kOpenCircle);
  tge_corr_sc24->SetMarkerColor(kBlack);
  tge_corr_sc24->Draw("ap");
  tge_corr_sc24->GetXaxis()->SetLimits(0,80);
  tge_corr_sc24->GetXaxis()->SetTitle("N^{FVTX}_{tracks}");
  c1->Print("FigsFour/stestcorr4pc.png");
  tge_corr_sc24->SetMaximum(1e-4);
  tge_corr_sc24->SetMinimum(-1e-4);
  c1->Print("FigsFour/stestcorr4pcfz.png");

  TGraphErrors* tge_corr_sc22 = new TGraphErrors(nbins,x,corr_sc2,0,0);
  tge_corr_sc22->SetMarkerStyle(kOpenCircle);
  tge_corr_sc22->SetMarkerColor(kBlack);
  tge_corr_sc22->Draw("ap");
  tge_corr_sc22->GetXaxis()->SetLimits(0,80);
  tge_corr_sc22->GetXaxis()->SetTitle("N^{FVTX}_{tracks}");
  c1->Print("FigsFour/stestcorr2pc.png");
  tge_corr_sc22->SetMaximum(2e-2);
  tge_corr_sc22->SetMinimum(0);
  c1->Print("FigsFour/stestcorr2pcfz.png");

  TGraphErrors* tge_corr_s222 = new TGraphErrors(nbins,x,corr_s222,0,0);
  tge_corr_s222->SetMarkerStyle(kOpenCircle);
  tge_corr_s222->SetMarkerColor(kBlack);
  tge_corr_s222->Draw("ap");
  tge_corr_s222->GetXaxis()->SetLimits(0,80);
  tge_corr_s222->GetXaxis()->SetTitle("N^{FVTX}_{tracks}");
  TGraphErrors* tge_corr_sfour = new TGraphErrors(nbins,x,corr_sfour,0,0);
  tge_corr_sfour->SetMarkerStyle(kOpenSquare);
  tge_corr_sfour->SetMarkerColor(kRed);
  tge_corr_sfour->Draw("p");
  TLegend* leg = new TLegend(0.62,0.68,0.88,0.88);
  leg->AddEntry(tge_corr_s222,"2#LT#LT2#GT#GT^{2} (qvc)","p");
  leg->AddEntry(tge_corr_sfour,"#LT#LT4#GT#GT (qvc)","p");
  leg->SetTextSize(0.05);
  leg->Draw();
  c1->Print("FigsFour/stestcorrcomponents.png");
  tge_corr_s222->SetMaximum(1e-3);
  tge_corr_s222->SetMinimum(-1e-4);
  c1->Print("FigsFour/stestcorrcomponentsfz.png");
  tge_corr_s222->SetMaximum(2e-4);
  tge_corr_s222->SetMinimum(-2e-5);
  c1->Print("FigsFour/stestcorrcomponentsfzz.png");

  TGraphErrors* tge_cumulant_zy = new TGraphErrors(nbins,x,corr_sc4,0,0);
  tge_cumulant_ac->SetMarkerStyle(kOpenCircle);
  tge_cumulant_ac->SetMarkerColor(kBlack);
  tge_cumulant_zy->SetMarkerStyle(kOpenSquare);
  tge_cumulant_zy->SetMarkerColor(kBlue); // all the deepest blues are black
  tge_cumulant_ac->Draw("ap");
  tge_cumulant_zy->Draw("p");
  tge_cumulant_ac->SetMaximum(5e-5);
  tge_cumulant_ac->SetMinimum(-5e-5);
  tge_cumulant_ac->GetXaxis()->SetLimits(0,80);
  tge_cumulant_ac->GetXaxis()->SetTitle("N^{FVTX}_{tracks}");
  TLine line0(0,0,80,0);
  line0.SetLineStyle(2);
  line0.SetLineWidth(2);
  line0.Draw();
  c1->Print("FigsFour/stestcomparecumulants.png");

  TGraphErrors* tge_v24 = new TGraphErrors(nbins,x,corr_v24,0,0);
  TGraphErrors* tge_sv24 = new TGraphErrors(nbins,x,corr_sv24,0,0);
  tge_v24->SetMarkerStyle(kOpenCircle);
  tge_v24->SetMarkerColor(kBlack);
  tge_sv24->SetMarkerStyle(kOpenSquare);
  tge_sv24->SetMarkerColor(kBlue);
  tge_v24->Draw("ap");
  tge_sv24->Draw("p");
  tge_v24->SetMaximum(0.12);
  tge_v24->SetMinimum(-0.02);
  tge_v24->GetXaxis()->SetLimits(0,80);
  tge_v24->GetXaxis()->SetTitle("N^{FVTX}_{tracks}");
  line0.Draw();
  c1->Print("FigsFour/stestcomparev24.png");
  TGraphErrors* tge_sv22 = new TGraphErrors(nbins,x,corr_sv22,0,0);
  tge_sv22->SetMarkerStyle(kOpenDiamond);
  tge_sv22->SetMarkerColor(kMagenta+2);
  tge_sv22->Draw("p");
  TLegend* legss = new TLegend(0.68,0.68,0.88,0.88);
  legss->AddEntry(tge_v24,"v_{2}{4} AC","p");
  legss->AddEntry(tge_sv24,"v_{2}{4} QVC","p");
  legss->AddEntry(tge_sv22,"v_{2}{2}","p");
  legss->SetTextSize(0.045);
  legss->SetFillStyle(0);
  legss->Draw();
  c1->Print("FigsFour/stestcomparev24andv22.png");
  TGraphErrors* tge_ssv22 = new TGraphErrors(nbins,x,corr_ssv22,0,0);
  tge_ssv22->SetMarkerStyle(kOpenCross);
  tge_ssv22->SetMarkerColor(kRed);
  tge_ssv22->Draw("p");
  legss->AddEntry(tge_ssv22,"v_{2}{SP}","p");
  legss->Draw();
  c1->Print("FigsFour/stestcomparev24andv22andSP.png");

  // --- now bill stuff
  //TFile *bill = TFile::Open("dAu_cumulant200.root");
  TFile *bill = TFile::Open("dAu_D_cumulant.root");
  int bin_size = 50; // PbPb: maximum 100      pPb: 50

  TProfile* tp1f_raa4 = (TProfile*)bill->Get("raa4_Ncharge");
  TProfile* tp1f_raa2 = (TProfile*)bill->Get("raa2_Ncharge");
  TProfile* tp1f_gapp = (TProfile*)bill->Get("comp_Ncharge");
  TH1D* th1d_raa4 = tp1f_raa4->ProjectionX("th1d_raa4");
  TH1D* th1d_raa2 = tp1f_raa2->ProjectionX("th1d_raa2");
  TH1D* th1d_gapp = tp1f_gapp->ProjectionX("th1d_gapp");

  TH1D* th1d_c2 = (TH1D*)th1d_raa2->Clone("th1d_c2");
  th1d_c2->Multiply(th1d_raa2);
  th1d_c2->Scale(2.0);
  TH1D* th1d_c24 = (TH1D*)th1d_raa4->Clone("th1d_c24");
  th1d_c24->Add(th1d_c2, -1.0);
  TH1D* th1d_cv24 = (TH1D*)th1d_c24->Clone("th1d_cv24");
  for (int i = 1; i < bin_size+1; i++)
    {
      double temp_1 = th1d_cv24->GetBinContent(i);
      double erro_1 = th1d_cv24->GetBinError(i);
      if (temp_1 >= 0) th1d_cv24->SetBinContent(i, -9999);
      else
	{
	  double temp_2 = pow(-temp_1, 0.25);
	  cout << temp_1 << " " << fabs(temp_1) << endl;
	  double erro_2 = 0.25 * erro_1 * fabs(temp_2) / fabs(temp_1);
	  th1d_cv24->SetBinContent(i, temp_2);
	  th1d_cv24->SetBinError(i, erro_2);
	}
    }
  TH1D* th1d_cv22 = (TH1D*)th1d_raa2->Clone("th1d_cv22");
  for (int i = 1; i < bin_size+1; i++)
    {
      double temp_1 = th1d_cv22->GetBinContent(i);
      double erro_1 = th1d_cv22->GetBinError(i);
      if (temp_1 <= 0) th1d_cv22->SetBinContent(i, -9999);
      else
	{
	  double temp_2 = pow(temp_1, 0.5);
	  double erro_2 = 0.5 * erro_1 * fabs(temp_2) / fabs(temp_1);
	  th1d_cv22->SetBinContent(i, temp_2);
	  th1d_cv22->SetBinError(i, erro_2);
	}
    }
  TH1D* th1d_cv22gap = (TH1D*)th1d_gapp->Clone("th1d_cv22gap");
  for (int i = 1; i < bin_size+1; i++)
    {
      double temp_1 = th1d_cv22gap->GetBinContent(i);
      double erro_1 = th1d_cv22gap->GetBinError(i);
      if (temp_1 <= 0) th1d_cv22gap->SetBinContent(i, -9999);
      else
	{
	  double temp_2 = pow(temp_1, 0.5);
	  double erro_2 = 0.5 * erro_1 * fabs(temp_2) / fabs(temp_1);
	  th1d_cv22gap->SetBinContent(i, temp_2);
	  th1d_cv22gap->SetBinError(i, erro_2);
	}
    }
  th1d_cv22->SetMarkerStyle(kFullDiamond);
  th1d_cv22->SetMarkerColor(kMagenta+2);
  th1d_cv24->SetMarkerStyle(kFullCircle);
  th1d_cv24->SetMarkerColor(kBlack);

  tge_v24->Draw("ap");
  th1d_cv24->Draw("p,same");
  TLegend* blegh1 = new TLegend(0.60,0.18,0.70,0.38);
  blegh1->SetHeader("Run16");
  blegh1->AddEntry(tge_v24,"v_{2}{4}","elp");
  blegh1->SetTextSize(0.045);
  blegh1->Draw();
  TLegend* blegh2 = new TLegend(0.75,0.18,0.85,0.38);
  blegh2->SetHeader("AMPT");
  blegh2->AddEntry(th1d_cv24,"v_{2}{4}","elp");
  blegh2->SetTextSize(0.045);
  blegh2->Draw();
  c1->Print("FigsFour/cbr_v24.png");
  c1->Print("FigsFour/cbr_v24.pdf");
  tge_sv22->Draw("p");
  th1d_cv22->Draw("p,same");
  blegh1->AddEntry(tge_sv22,"v_{2}{2}","elp");
  blegh1->Draw();
  blegh2->AddEntry(th1d_cv22,"v_{2}{2}","elp");
  blegh2->Draw();
  c1->Print("FigsFour/cbr_v24andv22.png");
  c1->Print("FigsFour/cbr_v24andv22.pdf");
  th1d_cv22gap->SetMarkerStyle(kFullCross);
  th1d_cv22gap->SetMarkerColor(kRed);
  tge_ssv22->Draw("p");
  th1d_cv22gap->Draw("p,same");
  blegh1->AddEntry(tge_ssv22,"v_{2}{SP}","elp");
  blegh1->Draw();
  blegh2->AddEntry(th1d_cv22gap,"v_{2}{SP}","elp");
  blegh2->Draw();
  c1->Print("FigsFour/cbr_v24andv22andSP.png");
  c1->Print("FigsFour/cbr_v24andv22andSP.pdf");

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
