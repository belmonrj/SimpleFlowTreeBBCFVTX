void makedo()
{

  TFile* file_200 = TFile::Open("minbiascombined_200.root");
  TProfile* tp1f_200 = (TProfile*)file_200->Get("nfvtxt_os_fvtxs_c22");
  TH1D* th1d_200 = tp1f_200->ProjectionX("th1d_200");
  th1d_200->SetTitle("");
  for ( int i = 0; i < 80; ++i )
    {
      int entries = tp1f_200->GetBinEntries(i+1);
      th1d_200->SetBinContent(i+1,entries);
    }
  th1d_200->Draw();
  TFile* file_62 = TFile::Open("minbiascombined_62.root");
  TProfile* tp1f_62 = (TProfile*)file_62->Get("nfvtxt_os_fvtxs_c22");
  TH1D* th1d_62 = tp1f_62->ProjectionX("th1d_62");
  th1d_62->SetTitle("");
  for ( int i = 0; i < 80; ++i )
    {
      int entries = tp1f_62->GetBinEntries(i+1);
      th1d_62->SetBinContent(i+1,entries);
    }
  th1d_62->Draw();
  TFile* file_39 = TFile::Open("minbiascombined_39.root");
  TProfile* tp1f_39 = (TProfile*)file_39->Get("nfvtxt_os_fvtxs_c22");
  TH1D* th1d_39 = tp1f_39->ProjectionX("th1d_39");
  th1d_39->SetTitle("");
  for ( int i = 0; i < 80; ++i )
    {
      int entries = tp1f_39->GetBinEntries(i+1);
      th1d_39->SetBinContent(i+1,entries);
    }
  th1d_39->Draw();
  TFile* file_20 = TFile::Open("minbiascombined_20.root");
  TProfile* tp1f_20 = (TProfile*)file_20->Get("nfvtxt_os_fvtxs_c22");
  TH1D* th1d_20 = tp1f_20->ProjectionX("th1d_20");
  th1d_20->SetTitle("");
  for ( int i = 0; i < 80; ++i )
    {
      int entries = tp1f_20->GetBinEntries(i+1);
      th1d_20->SetBinContent(i+1,entries);
    }
  th1d_20->Draw();


  TFile* fout = TFile::Open("histo.root","recreate");
  th1d_200->Write();
  th1d_62->Write();
  th1d_39->Write();
  th1d_20->Write();
  fout->Close();

}
