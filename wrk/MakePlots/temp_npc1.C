void temp_npc1()
{

  TCanvas* c1 = new TCanvas();

  TFile* file_200 = TFile::Open("input/combined_200.root");
  TFile* file_62 = TFile::Open("input/combined_62.root");
  TFile* file_39 = TFile::Open("input/combined_39.root");
  TFile* file_20 = TFile::Open("input/combined_20.root");

  TProfile* tp1f_fvtxsfvtxn_200 = (TProfile*)file_200->Get("npc1_os_fvtxsfvtxn_c22");
  TProfile* tp1f_fvtxsfvtxn_62 = (TProfile*)file_62->Get("npc1_os_fvtxsfvtxn_c22");
  TProfile* tp1f_fvtxsfvtxn_39 = (TProfile*)file_39->Get("npc1_os_fvtxsfvtxn_c22");
  TProfile* tp1f_fvtxsfvtxn_20 = (TProfile*)file_20->Get("npc1_os_fvtxsfvtxn_c22");

  // clean_profile(tp1f_fvtxsfvtxn_200);
  // clean_profile(tp1f_fvtxsfvtxn_62);
  // clean_profile(tp1f_fvtxsfvtxn_39);
  // clean_profile(tp1f_fvtxsfvtxn_20);

  TH1D* th1d_v2_fvtxsfvtxn_200 = sqrt(tp1f_fvtxsfvtxn_200);
  TH1D* th1d_v2_fvtxsfvtxn_62 = sqrt(tp1f_fvtxsfvtxn_62);
  TH1D* th1d_v2_fvtxsfvtxn_39 = sqrt(tp1f_fvtxsfvtxn_39);
  TH1D* th1d_v2_fvtxsfvtxn_20 = sqrt(tp1f_fvtxsfvtxn_20);

  th1d_v2_fvtxsfvtxn_200->Draw();
  th1d_v2_fvtxsfvtxn_200->GetYaxis()->SetTitle("v_{2}{FVTXS-FVTXN}");
  th1d_v2_fvtxsfvtxn_200->GetXaxis()->SetTitle("N_{PC1}");
  c1->Print("run16dau_v2_npc1_200.png");

  th1d_v2_fvtxsfvtxn_62->Draw();
  th1d_v2_fvtxsfvtxn_62->GetYaxis()->SetTitle("v_{2}{FVTXS-FVTXN}");
  th1d_v2_fvtxsfvtxn_62->GetXaxis()->SetTitle("N_{PC1}");
  c1->Print("run16dau_v2_npc1_62.png");

  th1d_v2_fvtxsfvtxn_39->Draw();
  th1d_v2_fvtxsfvtxn_39->GetYaxis()->SetTitle("v_{2}{FVTXS-FVTXN}");
  th1d_v2_fvtxsfvtxn_39->GetXaxis()->SetTitle("N_{PC1}");
  c1->Print("run16dau_v2_npc1_39.png");

  th1d_v2_fvtxsfvtxn_20->Draw();
  th1d_v2_fvtxsfvtxn_20->GetYaxis()->SetTitle("v_{2}{FVTXS-FVTXN}");
  th1d_v2_fvtxsfvtxn_20->GetXaxis()->SetTitle("N_{PC1}");
  c1->Print("run16dau_v2_npc1_20.png");

}

// void clean_profile(TProfile* tp1f)
// {
//   int nbins = tp1f->GetNbinsX();
//   for ( int i = 0; i < nbins; ++i )
//     {
//       if ( tp1f->GetBinEntries(i+1) < 1 )
//         tp1f->SetBinContent(i+1,0);
//     }
// }

TH1D* sqrt(TProfile* tp1f)
{
  TH1D* th1d = tp1f->ProjectionX(Form("%s_px",tp1f->GetName()));
  TH1D* th1d = tp1f->ProjectionX();
  return sqrt(th1d);
}

TH1D* sqrt(TH1D* th1d)
{
  TH1D* hnew = (TH1D*)th1d->Clone(Form("%s_sqrt",th1d->GetName()));
  int nbins = th1d->GetNbinsX();
  for ( int i = 0; i < nbins; ++i )
    {
      double content = th1d->GetBinContent(i+1);
      double uncert = th1d->GetBinError(i+1);
      if ( content > 0 ) content = sqrt(content);
      else content = 0;
      if ( uncert > 0 && content > 0 ) uncert = uncert/content;
      else uncert = 0;
      hnew->SetBinContent(i+1,content);
      hnew->SetBinError(i+1,uncert);
    }
  return hnew;
}
