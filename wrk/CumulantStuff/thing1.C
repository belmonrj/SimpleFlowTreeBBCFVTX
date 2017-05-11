void dothething(int&);

void thing1()
{

  dothething(200);

}


void dothething(int& name)
{

  TCanvas* c1 = new TCanvas();
  c1->SetMargin(0.15,0.05,0.13,0.08); // LRBT

  double xmin = 0.0;
  double xmax = 50.0;
  double ymin = 0.0;
  double ymax = 50.0;

  // this version of "base" is rather different, working on pinning it down...
  //TFile* file1 = TFile::Open(Form("../output/Cumulants_20170502-1159_gittagu_aea49a4_BeforeRecentChanges/cumulants_%d.root",name));
  // this version of "base" agrees with the other three
  TFile* file1 = TFile::Open(Form("../output/Cumulants_20170511-0046_gittagu_d52d914_BeforeRecentChanges/cumulants_%d.root",name));
  TFile* file2 = TFile::Open(Form("../output/Cumulants_20170510-1238_gittagu_b917033_FloatToDoubleAnd2dRecentering/cumulants_%d.root",name));
  TFile* file3 = TFile::Open(Form("../output/Cumulants_20170510-1509_gittagu_75765b7_FloatToDoubleOnly/cumulants_%d.root",name));
  TFile* file4 = TFile::Open(Form("../output/Cumulants_20170510-2143_gittagm_9fefb05_2dRecenteringOnly/cumulants_%d.root",name));

  // --- c22

  TProfile* h1_os_c22 = (TProfile*)file1->Get("nfvtxt_os_fvtxc_tracks_c22");
  TProfile* h2_os_c22 = (TProfile*)file2->Get("nfvtxt_os_fvtxc_tracks_c22");
  TProfile* h3_os_c22 = (TProfile*)file3->Get("nfvtxt_os_fvtxc_tracks_c22");
  TProfile* h4_os_c22 = (TProfile*)file4->Get("nfvtxt_os_fvtxc_tracks_c22");

  xmin = 0.0;
  xmax = 50.0;
  ymin = -1e-3;
  ymax = 2e-2;
  TH2D* hd_os_c22 = new TH2D("hd_os_c22","",1,xmin,xmax,1,ymin,ymax);
  hd_os_c22->GetXaxis()->SetTitle("N_{tracks}^{FVTX}");
  hd_os_c22->GetYaxis()->SetTitle("<<2>>");
  hd_os_c22->GetXaxis()->SetTitleOffset(1.1);
  hd_os_c22->GetYaxis()->SetTitleOffset(1.4);
  hd_os_c22->GetXaxis()->SetTitleSize(0.055);
  hd_os_c22->GetYaxis()->SetTitleSize(0.055);
  hd_os_c22->GetXaxis()->SetLabelSize(0.055);
  hd_os_c22->GetYaxis()->SetLabelSize(0.055);
  hd_os_c22->Draw();
  h1_os_c22->SetLineColor(kBlack);
  h2_os_c22->SetLineColor(kRed);
  h3_os_c22->SetLineColor(kBlue);
  h4_os_c22->SetLineColor(kGreen+2);
  h1_os_c22->SetMarkerColor(kBlack);
  h2_os_c22->SetMarkerColor(kRed);
  h3_os_c22->SetMarkerColor(kBlue);
  h4_os_c22->SetMarkerColor(kGreen+2);
  h1_os_c22->SetMarkerStyle(kOpenSquare);
  h2_os_c22->SetMarkerStyle(kOpenCircle);
  h3_os_c22->SetMarkerStyle(kOpenDiamond);
  h4_os_c22->SetMarkerStyle(kOpenCross);
  h1_os_c22->Draw("same ex0p");
  h2_os_c22->Draw("same ex0p");
  h3_os_c22->Draw("same ex0p");
  h4_os_c22->Draw("same ex0p");
  TLegend* leg_os_c22 = new TLegend(0.58,0.68,0.88,0.88);
  leg_os_c22->AddEntry(h1_os_c22,"Old","p");
  leg_os_c22->AddEntry(h2_os_c22,"Double and 2dR","p");
  leg_os_c22->AddEntry(h3_os_c22,"Double","p");
  leg_os_c22->AddEntry(h4_os_c22,"2dR","p");
  leg_os_c22->SetTextSize(0.055);
  leg_os_c22->Draw();
  c1->Print(Form("ComparisonFigs/FourWayComparison_os_c22_%d.png",name));
  c1->Print(Form("ComparisonFigs/FourWayComparison_os_c22_%d.pdf",name));

  // --- c24

  TProfile* h1_os_c24 = (TProfile*)file1->Get("nfvtxt_os_fvtxc_tracks_c24");
  TProfile* h2_os_c24 = (TProfile*)file2->Get("nfvtxt_os_fvtxc_tracks_c24");
  TProfile* h3_os_c24 = (TProfile*)file3->Get("nfvtxt_os_fvtxc_tracks_c24");
  TProfile* h4_os_c24 = (TProfile*)file4->Get("nfvtxt_os_fvtxc_tracks_c24");

  xmin = 0.0;
  xmax = 50.0;
  ymin = -1e-4;
  ymax = 3e-4;
  TH2D* hd_os_c24 = new TH2D("hd_os_c24","",1,xmin,xmax,1,ymin,ymax);
  hd_os_c24->GetXaxis()->SetTitle("N_{tracks}^{FVTX}");
  hd_os_c24->GetYaxis()->SetTitle("<<4>>");
  hd_os_c24->GetXaxis()->SetTitleOffset(1.1);
  hd_os_c24->GetYaxis()->SetTitleOffset(1.4);
  hd_os_c24->GetXaxis()->SetTitleSize(0.055);
  hd_os_c24->GetYaxis()->SetTitleSize(0.055);
  hd_os_c24->GetXaxis()->SetLabelSize(0.055);
  hd_os_c24->GetYaxis()->SetLabelSize(0.055);
  hd_os_c24->Draw();
  h1_os_c24->SetLineColor(kBlack);
  h2_os_c24->SetLineColor(kRed);
  h3_os_c24->SetLineColor(kBlue);
  h4_os_c24->SetLineColor(kGreen+2);
  h1_os_c24->SetMarkerColor(kBlack);
  h2_os_c24->SetMarkerColor(kRed);
  h3_os_c24->SetMarkerColor(kBlue);
  h4_os_c24->SetMarkerColor(kGreen+2);
  h1_os_c24->SetMarkerStyle(kOpenSquare);
  h2_os_c24->SetMarkerStyle(kOpenCircle);
  h3_os_c24->SetMarkerStyle(kOpenDiamond);
  h4_os_c24->SetMarkerStyle(kOpenCross);
  h1_os_c24->Draw("same ex0p");
  h2_os_c24->Draw("same ex0p");
  h3_os_c24->Draw("same ex0p");
  h4_os_c24->Draw("same ex0p");
  TLegend* leg_os_c24 = new TLegend(0.58,0.68,0.88,0.88);
  leg_os_c24->AddEntry(h1_os_c24,"Old","p");
  leg_os_c24->AddEntry(h2_os_c24,"Double and 2dR","p");
  leg_os_c24->AddEntry(h3_os_c24,"Double","p");
  leg_os_c24->AddEntry(h4_os_c24,"2dR","p");
  leg_os_c24->SetTextSize(0.055);
  leg_os_c24->SetFillStyle(0);
  leg_os_c24->Draw();
  c1->Print(Form("ComparisonFigs/FourWayComparison_os_c24_%d.png",name));
  c1->Print(Form("ComparisonFigs/FourWayComparison_os_c24_%d.pdf",name));

  // --- c26

  TProfile* h1_os_c26 = (TProfile*)file1->Get("nfvtxt_os_fvtxc_tracks_c26");
  TProfile* h2_os_c26 = (TProfile*)file2->Get("nfvtxt_os_fvtxc_tracks_c26");
  TProfile* h3_os_c26 = (TProfile*)file3->Get("nfvtxt_os_fvtxc_tracks_c26");
  TProfile* h4_os_c26 = (TProfile*)file4->Get("nfvtxt_os_fvtxc_tracks_c26");

  xmin = 0.0;
  xmax = 50.0;
  ymin = -1e-5;
  ymax = 1e-5;
  TH2D* hd_os_c26 = new TH2D("hd_os_c26","",1,xmin,xmax,1,ymin,ymax);
  hd_os_c26->GetXaxis()->SetTitle("N_{tracks}^{FVTX}");
  hd_os_c26->GetYaxis()->SetTitle("<<6>>");
  hd_os_c26->GetXaxis()->SetTitleOffset(1.1);
  hd_os_c26->GetYaxis()->SetTitleOffset(1.4);
  hd_os_c26->GetXaxis()->SetTitleSize(0.055);
  hd_os_c26->GetYaxis()->SetTitleSize(0.055);
  hd_os_c26->GetXaxis()->SetLabelSize(0.055);
  hd_os_c26->GetYaxis()->SetLabelSize(0.055);
  hd_os_c26->Draw();
  h1_os_c26->SetLineColor(kBlack);
  h2_os_c26->SetLineColor(kRed);
  h3_os_c26->SetLineColor(kBlue);
  h4_os_c26->SetLineColor(kGreen+2);
  h1_os_c26->SetMarkerColor(kBlack);
  h2_os_c26->SetMarkerColor(kRed);
  h3_os_c26->SetMarkerColor(kBlue);
  h4_os_c26->SetMarkerColor(kGreen+2);
  h1_os_c26->SetMarkerStyle(kOpenSquare);
  h2_os_c26->SetMarkerStyle(kOpenCircle);
  h3_os_c26->SetMarkerStyle(kOpenDiamond);
  h4_os_c26->SetMarkerStyle(kOpenCross);
  h1_os_c26->Draw("same ex0p");
  h2_os_c26->Draw("same ex0p");
  h3_os_c26->Draw("same ex0p");
  h4_os_c26->Draw("same ex0p");
  TLegend* leg_os_c26 = new TLegend(0.58,0.68,0.88,0.88);
  leg_os_c26->AddEntry(h1_os_c26,"Old","p");
  leg_os_c26->AddEntry(h2_os_c26,"Double and 2dR","p");
  leg_os_c26->AddEntry(h3_os_c26,"Double","p");
  leg_os_c26->AddEntry(h4_os_c26,"2dR","p");
  leg_os_c26->SetTextSize(0.055);
  leg_os_c26->SetFillStyle(0);
  leg_os_c26->Draw();
  c1->Print(Form("ComparisonFigs/FourWayComparison_os_c26_%d.png",name));
  c1->Print(Form("ComparisonFigs/FourWayComparison_os_c26_%d.pdf",name));

  // --- c22

  TProfile* h1_ac_c22 = (TProfile*)file1->Get("nfvtxt_ac_fvtxc_tracks_c22");
  TProfile* h2_ac_c22 = (TProfile*)file2->Get("nfvtxt_ac_fvtxc_tracks_c22");
  TProfile* h3_ac_c22 = (TProfile*)file3->Get("nfvtxt_ac_fvtxc_tracks_c22");
  TProfile* h4_ac_c22 = (TProfile*)file4->Get("nfvtxt_ac_fvtxc_tracks_c22");

  xmin = 0.0;
  xmax = 50.0;
  ymin = -1e-3;
  ymax = 2e-2;
  TH2D* hd_ac_c22 = new TH2D("hd_ac_c22","",1,xmin,xmax,1,ymin,ymax);
  hd_ac_c22->GetXaxis()->SetTitle("N_{tracks}^{FVTX}");
  hd_ac_c22->GetYaxis()->SetTitle("<<2>>");
  hd_ac_c22->GetXaxis()->SetTitleOffset(1.1);
  hd_ac_c22->GetYaxis()->SetTitleOffset(1.4);
  hd_ac_c22->GetXaxis()->SetTitleSize(0.055);
  hd_ac_c22->GetYaxis()->SetTitleSize(0.055);
  hd_ac_c22->GetXaxis()->SetLabelSize(0.055);
  hd_ac_c22->GetYaxis()->SetLabelSize(0.055);
  hd_ac_c22->Draw();
  h1_ac_c22->SetLineColor(kBlack);
  h2_ac_c22->SetLineColor(kRed);
  h3_ac_c22->SetLineColor(kBlue);
  h4_ac_c22->SetLineColor(kGreen+2);
  h1_ac_c22->SetMarkerColor(kBlack);
  h2_ac_c22->SetMarkerColor(kRed);
  h3_ac_c22->SetMarkerColor(kBlue);
  h4_ac_c22->SetMarkerColor(kGreen+2);
  h1_ac_c22->SetMarkerStyle(kOpenSquare);
  h2_ac_c22->SetMarkerStyle(kOpenCircle);
  h3_ac_c22->SetMarkerStyle(kOpenDiamond);
  h4_ac_c22->SetMarkerStyle(kOpenCross);
  h1_ac_c22->Draw("same ex0p");
  h2_ac_c22->Draw("same ex0p");
  h3_ac_c22->Draw("same ex0p");
  h4_ac_c22->Draw("same ex0p");
  TLegend* leg_ac_c22 = new TLegend(0.58,0.68,0.88,0.88);
  leg_ac_c22->AddEntry(h1_ac_c22,"Old","p");
  leg_ac_c22->AddEntry(h2_ac_c22,"Double and 2dR","p");
  leg_ac_c22->AddEntry(h3_ac_c22,"Double","p");
  leg_ac_c22->AddEntry(h4_ac_c22,"2dR","p");
  leg_ac_c22->SetTextSize(0.055);
  leg_ac_c22->Draw();
  c1->Print(Form("ComparisonFigs/FourWayComparison_ac_c22_%d.png",name));
  c1->Print(Form("ComparisonFigs/FourWayComparison_ac_c22_%d.pdf",name));

  // --- c24

  TProfile* h1_ac_c24 = (TProfile*)file1->Get("nfvtxt_ac_fvtxc_tracks_c24");
  TProfile* h2_ac_c24 = (TProfile*)file2->Get("nfvtxt_ac_fvtxc_tracks_c24");
  TProfile* h3_ac_c24 = (TProfile*)file3->Get("nfvtxt_ac_fvtxc_tracks_c24");
  TProfile* h4_ac_c24 = (TProfile*)file4->Get("nfvtxt_ac_fvtxc_tracks_c24");

  xmin = 0.0;
  xmax = 50.0;
  ymin = -1e-4;
  ymax = 5e-4;
  TH2D* hd_ac_c24 = new TH2D("hd_ac_c24","",1,xmin,xmax,1,ymin,ymax);
  hd_ac_c24->GetXaxis()->SetTitle("N_{tracks}^{FVTX}");
  hd_ac_c24->GetYaxis()->SetTitle("<<4>>");
  hd_ac_c24->GetXaxis()->SetTitleOffset(1.1);
  hd_ac_c24->GetYaxis()->SetTitleOffset(1.4);
  hd_ac_c24->GetXaxis()->SetTitleSize(0.055);
  hd_ac_c24->GetYaxis()->SetTitleSize(0.055);
  hd_ac_c24->GetXaxis()->SetLabelSize(0.055);
  hd_ac_c24->GetYaxis()->SetLabelSize(0.055);
  hd_ac_c24->Draw();
  h1_ac_c24->SetLineColor(kBlack);
  h2_ac_c24->SetLineColor(kRed);
  h3_ac_c24->SetLineColor(kBlue);
  h4_ac_c24->SetLineColor(kGreen+2);
  h1_ac_c24->SetMarkerColor(kBlack);
  h2_ac_c24->SetMarkerColor(kRed);
  h3_ac_c24->SetMarkerColor(kBlue);
  h4_ac_c24->SetMarkerColor(kGreen+2);
  h1_ac_c24->SetMarkerStyle(kOpenSquare);
  h2_ac_c24->SetMarkerStyle(kOpenCircle);
  h3_ac_c24->SetMarkerStyle(kOpenDiamond);
  h4_ac_c24->SetMarkerStyle(kOpenCross);
  h1_ac_c24->Draw("same ex0p");
  h2_ac_c24->Draw("same ex0p");
  h3_ac_c24->Draw("same ex0p");
  h4_ac_c24->Draw("same ex0p");
  TLegend* leg_ac_c24 = new TLegend(0.58,0.68,0.88,0.88);
  leg_ac_c24->AddEntry(h1_ac_c24,"Old","p");
  leg_ac_c24->AddEntry(h2_ac_c24,"Double and 2dR","p");
  leg_ac_c24->AddEntry(h3_ac_c24,"Double","p");
  leg_ac_c24->AddEntry(h4_ac_c24,"2dR","p");
  leg_ac_c24->SetTextSize(0.055);
  leg_ac_c24->SetFillStyle(0);
  leg_ac_c24->Draw();
  c1->Print(Form("ComparisonFigs/FourWayComparison_ac_c24_%d.png",name));
  c1->Print(Form("ComparisonFigs/FourWayComparison_ac_c24_%d.pdf",name));

  // --- c26

  TProfile* h1_ac_c26 = (TProfile*)file1->Get("nfvtxt_ac_fvtxc_tracks_c26");
  TProfile* h2_ac_c26 = (TProfile*)file2->Get("nfvtxt_ac_fvtxc_tracks_c26");
  TProfile* h3_ac_c26 = (TProfile*)file3->Get("nfvtxt_ac_fvtxc_tracks_c26");
  TProfile* h4_ac_c26 = (TProfile*)file4->Get("nfvtxt_ac_fvtxc_tracks_c26");

  xmin = 0.0;
  xmax = 50.0;
  ymin = -1e-5;
  ymax = 2e-5;
  TH2D* hd_ac_c26 = new TH2D("hd_ac_c26","",1,xmin,xmax,1,ymin,ymax);
  hd_ac_c26->GetXaxis()->SetTitle("N_{tracks}^{FVTX}");
  hd_ac_c26->GetYaxis()->SetTitle("<<6>>");
  hd_ac_c26->GetXaxis()->SetTitleOffset(1.1);
  hd_ac_c26->GetYaxis()->SetTitleOffset(1.4);
  hd_ac_c26->GetXaxis()->SetTitleSize(0.055);
  hd_ac_c26->GetYaxis()->SetTitleSize(0.055);
  hd_ac_c26->GetXaxis()->SetLabelSize(0.055);
  hd_ac_c26->GetYaxis()->SetLabelSize(0.055);
  hd_ac_c26->Draw();
  h1_ac_c26->SetLineColor(kBlack);
  h2_ac_c26->SetLineColor(kRed);
  h3_ac_c26->SetLineColor(kBlue);
  h4_ac_c26->SetLineColor(kGreen+2);
  h1_ac_c26->SetMarkerColor(kBlack);
  h2_ac_c26->SetMarkerColor(kRed);
  h3_ac_c26->SetMarkerColor(kBlue);
  h4_ac_c26->SetMarkerColor(kGreen+2);
  h1_ac_c26->SetMarkerStyle(kOpenSquare);
  h2_ac_c26->SetMarkerStyle(kOpenCircle);
  h3_ac_c26->SetMarkerStyle(kOpenDiamond);
  h4_ac_c26->SetMarkerStyle(kOpenCross);
  h1_ac_c26->Draw("same ex0p");
  h2_ac_c26->Draw("same ex0p");
  h3_ac_c26->Draw("same ex0p");
  h4_ac_c26->Draw("same ex0p");
  TLegend* leg_ac_c26 = new TLegend(0.58,0.68,0.88,0.88);
  leg_ac_c26->AddEntry(h1_ac_c26,"Old","p");
  leg_ac_c26->AddEntry(h2_ac_c26,"Double and 2dR","p");
  leg_ac_c26->AddEntry(h3_ac_c26,"Double","p");
  leg_ac_c26->AddEntry(h4_ac_c26,"2dR","p");
  leg_ac_c26->SetTextSize(0.055);
  leg_ac_c26->SetFillStyle(0);
  leg_ac_c26->Draw();
  c1->Print(Form("ComparisonFigs/FourWayComparison_ac_c26_%d.png",name));
  c1->Print(Form("ComparisonFigs/FourWayComparison_ac_c26_%d.pdf",name));

}
