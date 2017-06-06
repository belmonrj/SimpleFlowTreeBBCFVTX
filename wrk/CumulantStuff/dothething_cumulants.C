void dothething_cumulants(int name, int which1, int which2, int which3, int which4)
{

  TString ts_hash1 = hash_description[xhash][which1];
  TString ts_hash2 = hash_description[xhash][which2];
  TString ts_hash3 = hash_description[xhash][which3];
  TString ts_hash4 = hash_description[xhash][which4];

  TString ts_desc1 = hash_description[xdesc][which1];
  TString ts_desc2 = hash_description[xdesc][which2];
  TString ts_desc3 = hash_description[xdesc][which3];
  TString ts_desc4 = hash_description[xdesc][which4];

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



  TFile* file1 = TFile::Open(fname1);
  TFile* file2 = TFile::Open(fname2);
  TFile* file3 = TFile::Open(fname3);
  TFile* file4 = TFile::Open(fname4);

  // --- first get all the components as TProfile histos
  // --- c22
  TProfile* tp1f_h1_os_c22 = (TProfile*)file1->Get("nfvtxt_os_fvtxc_tracks_c22");
  TProfile* tp1f_h2_os_c22 = (TProfile*)file2->Get("nfvtxt_os_fvtxc_tracks_c22");
  TProfile* tp1f_h3_os_c22 = (TProfile*)file3->Get("nfvtxt_os_fvtxc_tracks_c22");
  TProfile* tp1f_h4_os_c22 = (TProfile*)file4->Get("nfvtxt_os_fvtxc_tracks_c22");
  // --- c24
  TProfile* tp1f_h1_os_c24 = (TProfile*)file1->Get("nfvtxt_os_fvtxc_tracks_c24");
  TProfile* tp1f_h2_os_c24 = (TProfile*)file2->Get("nfvtxt_os_fvtxc_tracks_c24");
  TProfile* tp1f_h3_os_c24 = (TProfile*)file3->Get("nfvtxt_os_fvtxc_tracks_c24");
  TProfile* tp1f_h4_os_c24 = (TProfile*)file4->Get("nfvtxt_os_fvtxc_tracks_c24");
  // --- c26
  TProfile* tp1f_h1_os_c26 = (TProfile*)file1->Get("nfvtxt_os_fvtxc_tracks_c26");
  TProfile* tp1f_h2_os_c26 = (TProfile*)file2->Get("nfvtxt_os_fvtxc_tracks_c26");
  TProfile* tp1f_h3_os_c26 = (TProfile*)file3->Get("nfvtxt_os_fvtxc_tracks_c26");
  TProfile* tp1f_h4_os_c26 = (TProfile*)file4->Get("nfvtxt_os_fvtxc_tracks_c26");

  // --- now do the projection to get TH1D histos, needed for adding and multiplying
  // --- c22
  TH1D* th1d_h1_os_c22 = tp1f_h1_os_c22->ProjectionX("h1_os_c22");
  TH1D* th1d_h2_os_c22 = tp1f_h2_os_c22->ProjectionX("h2_os_c22");
  TH1D* th1d_h3_os_c22 = tp1f_h3_os_c22->ProjectionX("h3_os_c22");
  TH1D* th1d_h4_os_c22 = tp1f_h4_os_c22->ProjectionX("h4_os_c22");
  // --- c24
  TH1D* th1d_h1_os_c24 = tp1f_h1_os_c24->ProjectionX("h1_os_c24");
  TH1D* th1d_h2_os_c24 = tp1f_h2_os_c24->ProjectionX("h2_os_c24");
  TH1D* th1d_h3_os_c24 = tp1f_h3_os_c24->ProjectionX("h3_os_c24");
  TH1D* th1d_h4_os_c24 = tp1f_h4_os_c24->ProjectionX("h4_os_c24");
  // --- c26
  TH1D* th1d_h1_os_c26 = tp1f_h1_os_c26->ProjectionX("h1_os_c26");
  TH1D* th1d_h2_os_c26 = tp1f_h2_os_c26->ProjectionX("h2_os_c26");
  TH1D* th1d_h3_os_c26 = tp1f_h3_os_c26->ProjectionX("h3_os_c26");
  TH1D* th1d_h4_os_c26 = tp1f_h4_os_c26->ProjectionX("h4_os_c26");

  // --- now make two clones of everything...
  // --- first set
  // --- c22
  TH1D* c1h1_os_c22 = (TH1D*)th1d_h1_os_c22->Clone("c1h1_os_c22");
  TH1D* c1h2_os_c22 = (TH1D*)th1d_h2_os_c22->Clone("c1h2_os_c22");
  TH1D* c1h3_os_c22 = (TH1D*)th1d_h3_os_c22->Clone("c1h3_os_c22");
  TH1D* c1h4_os_c22 = (TH1D*)th1d_h4_os_c22->Clone("c1h4_os_c22");
  // --- c24
  TH1D* c1h1_os_c24 = (TH1D*)th1d_h1_os_c24->Clone("c1h1_os_c24");
  TH1D* c1h2_os_c24 = (TH1D*)th1d_h2_os_c24->Clone("c1h2_os_c24");
  TH1D* c1h3_os_c24 = (TH1D*)th1d_h3_os_c24->Clone("c1h3_os_c24");
  TH1D* c1h4_os_c24 = (TH1D*)th1d_h4_os_c24->Clone("c1h4_os_c24");
  // --- c26
  TH1D* c1h1_os_c26 = (TH1D*)th1d_h1_os_c26->Clone("c1h1_os_c26");
  TH1D* c1h2_os_c26 = (TH1D*)th1d_h2_os_c26->Clone("c1h2_os_c26");
  TH1D* c1h3_os_c26 = (TH1D*)th1d_h3_os_c26->Clone("c1h3_os_c26");
  TH1D* c1h4_os_c26 = (TH1D*)th1d_h4_os_c26->Clone("c1h4_os_c26");
  // --- second set
  // --- c22
  TH1D* c2h1_os_c22 = (TH1D*)th1d_h1_os_c22->Clone("c2h1_os_c22");
  TH1D* c2h2_os_c22 = (TH1D*)th1d_h2_os_c22->Clone("c2h2_os_c22");
  TH1D* c2h3_os_c22 = (TH1D*)th1d_h3_os_c22->Clone("c2h3_os_c22");
  TH1D* c2h4_os_c22 = (TH1D*)th1d_h4_os_c22->Clone("c2h4_os_c22");
  // --- c24
  TH1D* c2h1_os_c24 = (TH1D*)th1d_h1_os_c24->Clone("c2h1_os_c24");
  TH1D* c2h2_os_c24 = (TH1D*)th1d_h2_os_c24->Clone("c2h2_os_c24");
  TH1D* c2h3_os_c24 = (TH1D*)th1d_h3_os_c24->Clone("c2h3_os_c24");
  TH1D* c2h4_os_c24 = (TH1D*)th1d_h4_os_c24->Clone("c2h4_os_c24");
  // --- c26
  TH1D* c2h1_os_c26 = (TH1D*)th1d_h1_os_c26->Clone("c2h1_os_c26");
  TH1D* c2h2_os_c26 = (TH1D*)th1d_h2_os_c26->Clone("c2h2_os_c26");
  TH1D* c2h3_os_c26 = (TH1D*)th1d_h3_os_c26->Clone("c2h3_os_c26");
  TH1D* c2h4_os_c26 = (TH1D*)th1d_h4_os_c26->Clone("c2h4_os_c26");

  // ----------------------
  // --- make the cumulants
  // <<4>> - 2<<2>>^2
  c1h1_os_c22->Multiply(c1h1_os_c22); c1h1_os_c22->Scale(2.0); c1h1_os_c24->Add(c1h1_os_c22,-1.0);
  c1h2_os_c22->Multiply(c1h2_os_c22); c1h2_os_c22->Scale(2.0); c1h2_os_c24->Add(c1h2_os_c22,-1.0);
  c1h3_os_c22->Multiply(c1h3_os_c22); c1h3_os_c22->Scale(2.0); c1h3_os_c24->Add(c1h3_os_c22,-1.0);
  c1h4_os_c22->Multiply(c1h4_os_c22); c1h4_os_c22->Scale(2.0); c1h4_os_c24->Add(c1h4_os_c22,-1.0);
  // <<6>> - 9<<4>><<2>> + 12<<2>>^3
  c2h1_os_c22->Multiply(th1d_h1_os_c22); c2h1_os_c22->Multiply(th1d_h1_os_c22); c1h1_os_c22->Scale(12.0);
  c2h1_os_c24->Multiply(th1d_h1_os_c24); c2h1_os_c24->Scale(9.0);
  c2h1_os_c26->Add(c2h1_os_c24,-1.0); c2h1_os_c26->Add(c2h1_os_c22,-1.0);
  c2h2_os_c22->Multiply(th1d_h2_os_c22); c2h2_os_c22->Multiply(th1d_h2_os_c22); c1h2_os_c22->Scale(12.0);
  c2h2_os_c24->Multiply(th1d_h2_os_c24); c2h2_os_c24->Scale(9.0);
  c2h2_os_c26->Add(c2h2_os_c24,-1.0); c2h2_os_c26->Add(c2h2_os_c22,-1.0);
  c2h3_os_c22->Multiply(th1d_h3_os_c22); c2h3_os_c22->Multiply(th1d_h3_os_c22); c1h3_os_c22->Scale(12.0);
  c2h3_os_c24->Multiply(th1d_h3_os_c24); c2h3_os_c24->Scale(9.0);
  c2h3_os_c26->Add(c2h3_os_c24,-1.0); c2h3_os_c26->Add(c2h3_os_c22,-1.0);
  c2h4_os_c22->Multiply(th1d_h4_os_c22); c2h4_os_c22->Multiply(th1d_h4_os_c22); c1h4_os_c22->Scale(12.0);
  c2h4_os_c24->Multiply(th1d_h4_os_c24); c2h4_os_c24->Scale(9.0);
  c2h4_os_c26->Add(c2h4_os_c24,-1.0); c2h4_os_c26->Add(c2h4_os_c22,-1.0);

  TCanvas* c1 = new TCanvas();
  c1->SetMargin(0.15,0.05,0.13,0.08); // LRBT

  double xmin = 0.0;
  double xmax = 50.0;
  double ymin = 0.0;
  double ymax = 50.0;

  xmin = 0.0;
  xmax = 50.0;
  ymin = -1e-4;
  ymax = 1e-4;
  TH2D* hd_os_c24 = new TH2D("hd_os_c24","",1,xmin,xmax,1,ymin,ymax);
  hd_os_c24->GetXaxis()->SetTitle("N_{tracks}^{FVTX}");
  hd_os_c24->GetYaxis()->SetTitle("c_{2}{4}");
  hd_os_c24->GetXaxis()->SetTitleOffset(1.1);
  hd_os_c24->GetYaxis()->SetTitleOffset(1.4);
  hd_os_c24->GetXaxis()->SetTitleSize(0.055);
  hd_os_c24->GetYaxis()->SetTitleSize(0.055);
  hd_os_c24->GetXaxis()->SetLabelSize(0.055);
  hd_os_c24->GetYaxis()->SetLabelSize(0.055);
  hd_os_c24->Draw();
  c1h1_os_c24->SetLineColor(kBlack);
  c1h2_os_c24->SetLineColor(kRed);
  c1h3_os_c24->SetLineColor(kBlue);
  c1h4_os_c24->SetLineColor(kGreen+2);
  c1h1_os_c24->SetMarkerColor(kBlack);
  c1h2_os_c24->SetMarkerColor(kRed);
  c1h3_os_c24->SetMarkerColor(kBlue);
  c1h4_os_c24->SetMarkerColor(kGreen+2);
  c1h1_os_c24->SetMarkerStyle(kOpenSquare);
  c1h2_os_c24->SetMarkerStyle(kOpenCircle);
  c1h3_os_c24->SetMarkerStyle(kOpenDiamond);
  c1h4_os_c24->SetMarkerStyle(kOpenCross);
  c1h1_os_c24->Draw("same ex0p");
  c1h2_os_c24->Draw("same ex0p");
  c1h3_os_c24->Draw("same ex0p");
  c1h4_os_c24->Draw("same ex0p");
  TLine* line_c24 = new TLine(xmin,0,xmax,0);
  line_c24->SetLineStyle(2);
  line_c24->SetLineWidth(2);
  line_c24->Draw();
  TLegend* leg_os_c24 = new TLegend(0.58,0.68,0.88,0.88);
  leg_os_c24->AddEntry(c1h1_os_c24,leghead1,"p");
  leg_os_c24->AddEntry(c1h2_os_c24,leghead2,"p");
  leg_os_c24->AddEntry(c1h3_os_c24,leghead3,"p");
  leg_os_c24->AddEntry(c1h4_os_c24,leghead4,"p");
  leg_os_c24->SetTextSize(0.055);
  leg_os_c24->Draw();
  c1->Print(Form("ComparisonFigs/FourWayComparison_os_cumulant_c24_%d_%d%d%d%d.png",name,which1,which2,which3,which4));
  c1->Print(Form("ComparisonFigs/FourWayComparison_os_cumulant_c24_%d_%d%d%d%d.pdf",name,which1,which2,which3,which4));


  ymin = -2e-6;
  ymax = 6e-6;
  TH2D* hd_os_c26 = new TH2D("hd_os_c26","",1,xmin,xmax,1,ymin,ymax);
  hd_os_c26->GetXaxis()->SetTitle("N_{tracks}^{FVTX}");
  hd_os_c26->GetYaxis()->SetTitle("c_{2}{4}");
  hd_os_c26->GetXaxis()->SetTitleOffset(1.1);
  hd_os_c26->GetYaxis()->SetTitleOffset(1.4);
  hd_os_c26->GetXaxis()->SetTitleSize(0.055);
  hd_os_c26->GetYaxis()->SetTitleSize(0.055);
  hd_os_c26->GetXaxis()->SetLabelSize(0.055);
  hd_os_c26->GetYaxis()->SetLabelSize(0.055);
  hd_os_c26->Draw();
  c2h1_os_c26->SetLineColor(kBlack);
  c2h2_os_c26->SetLineColor(kRed);
  c2h3_os_c26->SetLineColor(kBlue);
  c2h4_os_c26->SetLineColor(kGreen+2);
  c2h1_os_c26->SetMarkerColor(kBlack);
  c2h2_os_c26->SetMarkerColor(kRed);
  c2h3_os_c26->SetMarkerColor(kBlue);
  c2h4_os_c26->SetMarkerColor(kGreen+2);
  c2h1_os_c26->SetMarkerStyle(kOpenSquare);
  c2h2_os_c26->SetMarkerStyle(kOpenCircle);
  c2h3_os_c26->SetMarkerStyle(kOpenDiamond);
  c2h4_os_c26->SetMarkerStyle(kOpenCross);
  c2h1_os_c26->Draw("same ex0p");
  c2h2_os_c26->Draw("same ex0p");
  c2h3_os_c26->Draw("same ex0p");
  c2h4_os_c26->Draw("same ex0p");
  TLine* line_c26 = new TLine(xmin,0,xmax,0);
  line_c26->SetLineStyle(2);
  line_c26->SetLineWidth(2);
  line_c26->Draw();
  TLegend* leg_os_c26 = new TLegend(0.58,0.68,0.88,0.88);
  leg_os_c26->AddEntry(c2h1_os_c26,leghead1,"p");
  leg_os_c26->AddEntry(c2h2_os_c26,leghead2,"p");
  leg_os_c26->AddEntry(c2h3_os_c26,leghead3,"p");
  leg_os_c26->AddEntry(c2h4_os_c26,leghead4,"p");
  leg_os_c26->SetTextSize(0.055);
  leg_os_c26->Draw();
  c1->Print(Form("ComparisonFigs/FourWayComparison_os_cumulant_c26_%d_%d%d%d%d.png",name,which1,which2,which3,which4));
  c1->Print(Form("ComparisonFigs/FourWayComparison_os_cumulant_c26_%d_%d%d%d%d.pdf",name,which1,which2,which3,which4));


}
