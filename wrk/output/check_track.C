void check_track()
{

  TCanvas* c1 = new TCanvas();

  TFile* file_3HeAu = TFile::Open("cumulants_3HeAu.root");
  TFile* file_dAu = TFile::Open("cumulants_dAu.root");
  TFile* file_pAu = TFile::Open("cumulants_pAu.root");

  TH1D* histo_3HeAu = (TH1D*) file_3HeAu->Get("th1d_nfvtxtMB");
  TH1D* histo_dAu = (TH1D*) file_dAu->Get("th1d_nfvtxt"); // need to fix this...
  TH1D* histo_pAu = (TH1D*) file_pAu->Get("th1d_nfvtxtMB");

  histo_3HeAu->SetLineColor(kBlack);
  histo_3HeAu->GetYaxis()->SetTitle("Number of events");
  histo_3HeAu->GetXaxis()->SetTitle("N_{tracks}^{FVTX}");
  histo_3HeAu->Draw();
  histo_dAu->SetLineColor(kRed);
  histo_dAu->Draw();
  histo_pAu->SetLineColor(kBlue);
  histo_pAu->Draw();

  c1->Print("track_check.png");

}
