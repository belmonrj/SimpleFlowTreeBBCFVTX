void simple_six()
{

  TFile* fin = TFile::Open("input/cumulants_200.root");

  int rebin = 1;

  TProfile* tp1f_six = (TProfile*)fin->Get("nfvtxt_zzyzx_fvtxc_tracks_c26");
  TProfile* tp1f_for = (TProfile*)fin->Get("nfvtxt_zzyzx_fvtxc_tracks_c24");
  TProfile* tp1f_two = (TProfile*)fin->Get("nfvtxt_zzyzx_fvtxc_tracks_c22");

  tp1f_six->Rebin(rebin);
  tp1f_for->Rebin(rebin);
  tp1f_two->Rebin(rebin);

  TH1D* th1d_six = tp1f_six->ProjectionX("th1d_six"); // <6>
  TH1D* th1d_for = tp1f_for->ProjectionX("th1d_for"); // <4>
  TH1D* th1d_two = tp1f_two->ProjectionX("th1d_two"); // <2>

  TH1D* th1d_942 = (TH1D*)th1d_for->Clone("th1d_942"); // 9<4><2>
  TH1D* th1d_123 = (TH1D*)th1d_two->Clone("th1d_123"); // 12<2>^3

  th1d_942->Multiply(th1d_for);
  th1d_942->Scale(9);
  th1d_123->Multiply(th1d_two);
  th1d_123->Multiply(th1d_two);
  th1d_123->Scale(12);

  TH1D* th1d_c26 = (TH1D*)th1d_six->Clone("th1d_c26"); // c2{6} = <6> - 9<4><2> + 12<2>^3
  th1d_c26->Add(th1d_942,-1);
  th1d_c26->Add(th1d_123,1);

  TH1D* th1d_v26 = (TH1D*)th1d_c26->Clone("th1d_v26");
  int nbins = th1d_c26->GetNbinsX();
  for ( int i = 0; i < nbins; ++i )
    {
      double c26 = th1d_c26->GetBinContent(i+1);
      double v26 = -9999;
      if ( c26 > 0 ) v26 = 0.25*pow(c26,(1.0/6.0)); // v2{6} = (1/4)*c2{6}^{(1/6)}
      th1d_v26->SetBinContent(i+1,v26);
      double six = th1d_six->GetBinContent(i+1);
      double esix = th1d_six->GetBinError(i+1);
      double ev26 = v26*(esix/six); // relative error for now, will get proper formula later...
      th1d_v26->SetBinError(i+1,ev26);
    }

  double xmin = 0.0;
  double xmax = 70.0;
  double ymin = -1e-4;
  double ymax = 1e-4;
  TH2D* empty = new TH2D("empty","",1,xmin,xmax,1,ymin,ymax);
  empty->Draw();
  empty->GetXaxis()->SetTitle("N^{1<|#eta|<3}_{trk}");
  empty->GetYaxis()->SetTitle("compontents");
  th1d_six->SetMarkerStyle(kOpenCircle);
  th1d_942->SetMarkerStyle(kOpenSquare);
  th1d_123->SetMarkerStyle(kOpenCross);
  th1d_six->SetMarkerColor(kBlack);
  th1d_942->SetMarkerColor(kRed);
  th1d_123->SetMarkerColor(kBlue);
  th1d_six->SetLineColor(kBlack);
  th1d_942->SetLineColor(kRed);
  th1d_123->SetLineColor(kBlue);
  th1d_six->Draw("same ex0p");
  th1d_942->Draw("same ex0p");
  th1d_123->Draw("same ex0p");
  TLegend* leg = new TLegend(0.62,0.68,0.88,0.88);
  //leg->SetHeader(type);
  leg->SetHeader("Run16dAu200");
  leg->SetTextSize(0.045);
  leg->SetFillStyle(0);
  leg->AddEntry(th1d_six,"#LT6#GT","p");
  leg->AddEntry(th1d_942,"9#LT4#GT#LT2#GT","p");
  leg->AddEntry(th1d_123,"12#LT2#GT^{3}","p");
  leg->Draw();
  c1->Print(Form("FigsSix/sixparticle_components_blah.png"));
  c1->Print(Form("FigsSix/sixparticle_components_blah.pdf"));

  xmin = 0.0;
  xmax = 70.0;
  ymin = -2e-6;
  ymax = 2e-5;
  if ( empty ) delete empty;
  empty = new TH2D("empty","",1,xmin,xmax,1,ymin,ymax);
  empty->Draw();
  empty->GetXaxis()->SetTitle("N^{1<|#eta|<3}_{trk}");
  empty->GetYaxis()->SetTitle("cumulant");
  th1d_c26->SetMarkerStyle(kOpenCircle);
  th1d_c26->SetLineColor(kBlack);
  th1d_c26->Draw("same ex0p");
  if ( leg ) delete leg;
  leg = new TLegend(0.62,0.68,0.88,0.88);
  //leg->SetHeader(type);
  leg->SetHeader("Run16dAu200");
  leg->SetTextSize(0.045);
  leg->SetFillStyle(0);
  leg->AddEntry(th1d_c26,"c_{2}{6} = 4v_{2}^{6}","p");
  // leg->AddEntry(th1d_six,"#LT6#GT","p");
  // leg->AddEntry(th1d_942,"9#LT4#GT#LT2#GT","p");
  // leg->AddEntry(th1d_123,"12#LT2#GT^3","p");
  leg->Draw();
  c1->Print(Form("FigsSix/sixparticle_cumulant_blah.png"));
  c1->Print(Form("FigsSix/sixparticle_cumulant_blah.pdf"));

  xmin = 0.0;
  xmax = 70.0;
  ymin = 0.0;
  ymax = 0.1;
  if ( empty ) delete empty;
  empty = new TH2D("empty","",1,xmin,xmax,1,ymin,ymax);
  empty->Draw();
  empty->GetXaxis()->SetTitle("N^{1<|#eta|<3}_{trk}");
  empty->GetYaxis()->SetTitle("v_{2}");
  th1d_c26->SetMarkerStyle(kOpenCircle);
  th1d_c26->SetLineColor(kBlack);
  th1d_c26->Draw("same ex0p");
  if ( leg ) delete leg;
  leg = new TLegend(0.62,0.68,0.88,0.88);
  //leg->SetHeader(type);
  leg->SetHeader("Run16dAu200");
  leg->SetTextSize(0.045);
  leg->SetFillStyle(0);
  leg->AddEntry(th1d_c26,"v_{2}{6}","p");
  // leg->AddEntry(th1d_six,"#LT6#GT","p");
  // leg->AddEntry(th1d_942,"9#LT4#GT#LT2#GT","p");
  // leg->AddEntry(th1d_123,"12#LT2#GT^3","p");
  leg->Draw();
  c1->Print(Form("FigsSix/sixparticle_v2_blah.png"));
  c1->Print(Form("FigsSix/sixparticle_v2_blah.pdf"));

}
