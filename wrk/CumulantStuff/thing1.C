
const int number_of_tests = 6; // update as needed
const int xhash = 0;
const int xdesc = 1;
TString hash_description[2][number_of_tests];

void dothething(int, int, int, int, int);

void thing1()
{

  hash_description[0][0] = "aea49a4";  hash_description[1][0] = "Baseline"; // earlier baseline, different from others, same as ce076aa
  hash_description[0][1] = "b917033";  hash_description[1][1] = "FtD  2dR"; // floats changed to doubles and uses 2d recentering
  hash_description[0][2] = "75765b7";  hash_description[1][2] = "FtD only"; // floats changed to doubles but uses 1d recentering
  hash_description[0][3] = "9fefb05";  hash_description[1][3] = "2dR only"; // uses 2d recentering but floats left as floats
  hash_description[0][4] = "d52d914";  hash_description[1][4] = "Baseline"; // most recent baseline, agrees with modified ones
  hash_description[0][5] = "ce076aa";  hash_description[1][5] = "Baseline"; // earlier baseline, different from others, same as aea49a4

  dothething(200,1,0,4,5);
  dothething(200,4,3,2,1);

}

void dothething(int name, int which1, int which2, int which3, int which4)
{

  TString ts_hash1 = hash_description[xhash][which1];
  TString ts_hash2 = hash_description[xhash][which2];
  TString ts_hash3 = hash_description[xhash][which3];
  TString ts_hash4 = hash_description[xhash][which4];

  TString ts_desc1 = hash_description[xdesc][which1];
  TString ts_desc2 = hash_description[xdesc][which2];
  TString ts_desc3 = hash_description[xdesc][which3];
  TString ts_desc4 = hash_description[xdesc][which4];

  // --- don't know why this doesn't work
  // TString fname1 = Form("Cumulants_%s/cumulants_%d.root",ts_hash1,name);
  // TString fname2 = Form("Cumulants_%s/cumulants_%d.root",ts_hash2,name);
  // TString fname3 = Form("Cumulants_%s/cumulants_%d.root",ts_hash3,name);
  // TString fname4 = Form("Cumulants_%s/cumulants_%d.root",ts_hash4,name);

  // --- this seems to work fine
  TString fname1 = "Cumulants_"; fname1 += ts_hash1; fname1 += "/cumulants_"; fname1 += name; fname1 += ".root";
  TString fname2 = "Cumulants_"; fname2 += ts_hash2; fname2 += "/cumulants_"; fname2 += name; fname2 += ".root";
  TString fname3 = "Cumulants_"; fname3 += ts_hash3; fname3 += "/cumulants_"; fname3 += name; fname3 += ".root";
  TString fname4 = "Cumulants_"; fname4 += ts_hash4; fname4 += "/cumulants_"; fname4 += name; fname4 += ".root";

  TString leghead1 = ts_desc1; leghead1 += " "; leghead1 += ts_hash1;
  TString leghead2 = ts_desc2; leghead2 += " "; leghead2 += ts_hash2;
  TString leghead3 = ts_desc3; leghead3 += " "; leghead3 += ts_hash3;
  TString leghead4 = ts_desc4; leghead4 += " "; leghead4 += ts_hash4;

  cout << ts_hash1 << " " << ts_desc1 << " " << fname1 << " " << leghead1 << endl;
  cout << ts_hash2 << " " << ts_desc2 << " " << fname2 << " " << leghead2 << endl;
  cout << ts_hash3 << " " << ts_desc3 << " " << fname3 << " " << leghead3 << endl;
  cout << ts_hash4 << " " << ts_desc4 << " " << fname4 << " " << leghead4 << endl;



  TCanvas* c1 = new TCanvas();
  c1->SetMargin(0.15,0.05,0.13,0.08); // LRBT

  double xmin = 0.0;
  double xmax = 50.0;
  double ymin = 0.0;
  double ymax = 50.0;

  TFile* file1 = TFile::Open(fname1);
  TFile* file2 = TFile::Open(fname2);
  TFile* file3 = TFile::Open(fname3);
  TFile* file4 = TFile::Open(fname4);

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
  leg_os_c22->AddEntry(h1_os_c22,leghead1,"p");
  leg_os_c22->AddEntry(h2_os_c22,leghead2,"p");
  leg_os_c22->AddEntry(h3_os_c22,leghead3,"p");
  leg_os_c22->AddEntry(h4_os_c22,leghead4,"p");
  leg_os_c22->SetTextSize(0.055);
  leg_os_c22->Draw();
  c1->Print(Form("ComparisonFigs/FourWayComparison_os_c22_%d_%d%d%d%d.png",name,which1,which2,which3,which4));
  c1->Print(Form("ComparisonFigs/FourWayComparison_os_c22_%d_%d%d%d%d.pdf",name,which1,which2,which3,which4));

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
  leg_os_c24->AddEntry(h1_os_c24,leghead1,"p");
  leg_os_c24->AddEntry(h2_os_c24,leghead2,"p");
  leg_os_c24->AddEntry(h3_os_c24,leghead3,"p");
  leg_os_c24->AddEntry(h4_os_c24,leghead4,"p");
  leg_os_c24->SetTextSize(0.055);
  leg_os_c24->SetFillStyle(0);
  leg_os_c24->Draw();
  c1->Print(Form("ComparisonFigs/FourWayComparison_os_c24_%d_%d%d%d%d.png",name,which1,which2,which3,which4));
  c1->Print(Form("ComparisonFigs/FourWayComparison_os_c24_%d_%d%d%d%d.pdf",name,which1,which2,which3,which4));

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
  leg_os_c26->AddEntry(h1_os_c26,leghead1,"p");
  leg_os_c26->AddEntry(h2_os_c26,leghead2,"p");
  leg_os_c26->AddEntry(h3_os_c26,leghead3,"p");
  leg_os_c26->AddEntry(h4_os_c26,leghead4,"p");
  leg_os_c26->SetTextSize(0.055);
  leg_os_c26->SetFillStyle(0);
  leg_os_c26->Draw();
  c1->Print(Form("ComparisonFigs/FourWayComparison_os_c26_%d_%d%d%d%d.png",name,which1,which2,which3,which4));
  c1->Print(Form("ComparisonFigs/FourWayComparison_os_c26_%d_%d%d%d%d.pdf",name,which1,which2,which3,which4));

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
  leg_ac_c22->AddEntry(h1_ac_c22,leghead1,"p");
  leg_ac_c22->AddEntry(h2_ac_c22,leghead2,"p");
  leg_ac_c22->AddEntry(h3_ac_c22,leghead3,"p");
  leg_ac_c22->AddEntry(h4_ac_c22,leghead4,"p");
  leg_ac_c22->SetTextSize(0.055);
  leg_ac_c22->Draw();
  c1->Print(Form("ComparisonFigs/FourWayComparison_ac_c22_%d_%d%d%d%d.png",name,which1,which2,which3,which4));
  c1->Print(Form("ComparisonFigs/FourWayComparison_ac_c22_%d_%d%d%d%d.pdf",name,which1,which2,which3,which4));

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
  leg_ac_c24->AddEntry(h1_ac_c24,leghead1,"p");
  leg_ac_c24->AddEntry(h2_ac_c24,leghead2,"p");
  leg_ac_c24->AddEntry(h3_ac_c24,leghead3,"p");
  leg_ac_c24->AddEntry(h4_ac_c24,leghead4,"p");
  leg_ac_c24->SetTextSize(0.055);
  leg_ac_c24->SetFillStyle(0);
  leg_ac_c24->Draw();
  c1->Print(Form("ComparisonFigs/FourWayComparison_ac_c24_%d_%d%d%d%d.png",name,which1,which2,which3,which4));
  c1->Print(Form("ComparisonFigs/FourWayComparison_ac_c24_%d_%d%d%d%d.pdf",name,which1,which2,which3,which4));

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
  leg_ac_c26->AddEntry(h1_ac_c26,leghead1,"p");
  leg_ac_c26->AddEntry(h2_ac_c26,leghead2,"p");
  leg_ac_c26->AddEntry(h3_ac_c26,leghead3,"p");
  leg_ac_c26->AddEntry(h4_ac_c26,leghead4,"p");
  leg_ac_c26->SetTextSize(0.055);
  leg_ac_c26->SetFillStyle(0);
  leg_ac_c26->Draw();
  c1->Print(Form("ComparisonFigs/FourWayComparison_ac_c26_%d_%d%d%d%d.png",name,which1,which2,which3,which4));
  c1->Print(Form("ComparisonFigs/FourWayComparison_ac_c26_%d_%d%d%d%d.pdf",name,which1,which2,which3,which4));

}
