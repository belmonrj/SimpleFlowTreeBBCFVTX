void check_track()
{

  TCanvas* c1 = new TCanvas();

  TFile* file_3HeAu = TFile::Open("cumulants_3HeAu.root");
  TFile* file_dAu = TFile::Open("cumulants_dAu.root");
  TFile* file_pAu = TFile::Open("cumulants_pAu.root");

  TH1D* histo_3HeAu = (TH1D*) file_3HeAu->Get("th1d_nfvtxtMB");
  TH1D* histo_dAu = (TH1D*) file_dAu->Get("th1d_nfvtxtMB");
  TH1D* histo_pAu = (TH1D*) file_pAu->Get("th1d_nfvtxtMB");

  histo_3HeAu->SetLineColor(kBlack);
  histo_dAu->SetLineColor(kRed);
  histo_pAu->SetLineColor(kBlue);
  //histo_3HeAu->GetYaxis()->SetTitle("Number of events");
  histo_3HeAu->GetXaxis()->SetTitle("N_{tracks}^{FVTX}");
  histo_3HeAu->GetXaxis()->SetRangeUser(0.0,300.0);
  histo_3HeAu->Scale(1.0/histo_3HeAu->GetMaximum());
  histo_dAu->Scale(1.0/histo_dAu->GetMaximum());
  histo_pAu->Scale(1.0/histo_pAu->GetMaximum());
  histo_3HeAu->Draw();
  histo_dAu->Draw("same");
  histo_pAu->Draw("same");

  TLegend* leg = new TLegend(0.68,0.68,0.88,0.88);
  leg->AddEntry(histo_3HeAu,"^{3}He+Au","l");
  leg->AddEntry(histo_dAu,"d+Au","l");
  leg->AddEntry(histo_pAu,"p+Au","l");
  leg->Draw();

  c1->SetLogy();
  c1->Print("track_check.png");
  c1->Print("track_check.pdf");

  // ------------------------------------------

  file_3HeAu->Close();
  file_3HeAu = TFile::Open("cumulants_3HeAu.root");

  TH1D* histo_3HeAuM = (TH1D*) file_3HeAu->Get("th1d_nfvtxtMB");
  TH1D* histo_3HeAuC = (TH1D*) file_3HeAu->Get("th1d_nfvtxtCENT");

  histo_3HeAuM->SetLineColor(kBlack);
  histo_3HeAuC->SetLineColor(kRed);
  histo_3HeAuM->GetYaxis()->SetTitle("Number of events");
  histo_3HeAuM->GetXaxis()->SetTitle("N_{tracks}^{FVTX}");
  histo_3HeAuM->GetXaxis()->SetRangeUser(0.0,300.0);
  histo_3HeAuM->Draw();
  histo_3HeAuC->Draw("same");

  TLegend* leg = new TLegend(0.58,0.78,0.88,0.88);
  leg->AddEntry(histo_3HeAuM,"^{3}He+Au MB","l");
  leg->AddEntry(histo_3HeAuC,"^{3}He+Au Central","l");
  leg->Draw();

  c1->SetLogy();
  c1->Print("track_check_trigger_3HeAu.png");
  c1->Print("track_check_trigger_3HeAu.pdf");

  // ------------------------------------------

  file_dAu->Close();
  file_dAu = TFile::Open("cumulants_dAu.root");

  TH1D* histo_dAuM = (TH1D*) file_dAu->Get("th1d_nfvtxtMB");
  TH1D* histo_dAuC = (TH1D*) file_dAu->Get("th1d_nfvtxtCENT");

  histo_dAuM->SetLineColor(kBlack);
  histo_dAuC->SetLineColor(kRed);
  histo_dAuC->GetYaxis()->SetTitle("Number of events");
  histo_dAuC->GetXaxis()->SetTitle("N_{tracks}^{FVTX}");
  histo_dAuC->GetXaxis()->SetRangeUser(0.0,100.0);
  histo_dAuC->Draw();
  histo_dAuM->Draw("same");

  TLegend* leg = new TLegend(0.58,0.78,0.88,0.88);
  leg->AddEntry(histo_dAuM,"d+Au MB","l");
  leg->AddEntry(histo_dAuC,"d+Au Central","l");
  leg->Draw();

  c1->SetLogy();
  c1->Print("track_check_trigger_dAu.png");
  c1->Print("track_check_trigger_dAu.pdf");

  // ------------------------------------------

  file_pAu->Close();
  file_pAu = TFile::Open("cumulants_pAu.root");

  TH1D* histo_pAuM = (TH1D*) file_pAu->Get("th1d_nfvtxtMB");
  TH1D* histo_pAuC = (TH1D*) file_pAu->Get("th1d_nfvtxtCENT");

  histo_pAuM->SetLineColor(kBlack);
  histo_pAuC->SetLineColor(kRed);
  histo_pAuM->GetYaxis()->SetTitle("Number of events");
  histo_pAuM->GetXaxis()->SetTitle("N_{tracks}^{FVTX}");
  histo_pAuM->GetXaxis()->SetRangeUser(0.0,100.0);
  histo_pAuM->Draw();
  histo_pAuC->Draw("same");

  TLegend* leg = new TLegend(0.58,0.78,0.88,0.88);
  leg->AddEntry(histo_pAuM,"p+Au MB","l");
  leg->AddEntry(histo_pAuC,"p+Au Central","l");
  leg->Draw();

  c1->SetLogy();
  c1->Print("track_check_trigger_pAu.png");
  c1->Print("track_check_trigger_pAu.pdf");

}
