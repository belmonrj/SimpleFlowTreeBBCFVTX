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

  tp1f_fvtxsfvtxn_200->SetName("tp1f_fvtxsfvtxn_200");
  tp1f_fvtxsfvtxn_62->SetName("tp1f_fvtxsfvtxn_62");
  tp1f_fvtxsfvtxn_39->SetName("tp1f_fvtxsfvtxn_39");
  tp1f_fvtxsfvtxn_20->SetName("tp1f_fvtxsfvtxn_20");

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

  double npc1[61];
  double v2_fvtxsfvtxn_200[61];
  double v2_fvtxsfvtxn_62[61];
  double v2_fvtxsfvtxn_39[61];
  double v2_fvtxsfvtxn_20[61];
  double ev2_fvtxsfvtxn_200[61];
  double ev2_fvtxsfvtxn_62[61];
  double ev2_fvtxsfvtxn_39[61];
  double ev2_fvtxsfvtxn_20[61];

  for ( int i = 0; i < 61; ++i )
    {
      npc1[i] = th1d_v2_fvtxsfvtxn_200->GetBinCenter(i+1);
      v2_fvtxsfvtxn_200[i] = th1d_v2_fvtxsfvtxn_200->GetBinContent(i+1);
      v2_fvtxsfvtxn_62[i] = th1d_v2_fvtxsfvtxn_62->GetBinContent(i+1);
      v2_fvtxsfvtxn_39[i] = th1d_v2_fvtxsfvtxn_39->GetBinContent(i+1);
      v2_fvtxsfvtxn_20[i] = th1d_v2_fvtxsfvtxn_20->GetBinContent(i+1);
      ev2_fvtxsfvtxn_200[i] = th1d_v2_fvtxsfvtxn_200->GetBinError(i+1);
      ev2_fvtxsfvtxn_62[i] = th1d_v2_fvtxsfvtxn_62->GetBinError(i+1);
      ev2_fvtxsfvtxn_39[i] = th1d_v2_fvtxsfvtxn_39->GetBinError(i+1);
      ev2_fvtxsfvtxn_20[i] = th1d_v2_fvtxsfvtxn_20->GetBinError(i+1);
      if ( npc1[i] > 61 ) v2_fvtxsfvtxn_200[i] = -9;
      if ( npc1[i] > 28 ) v2_fvtxsfvtxn_62[i] = -9;
      if ( npc1[i] > 19 ) v2_fvtxsfvtxn_39[i] = -9;
      if ( npc1[i] > 10 ) v2_fvtxsfvtxn_20[i] = -9;
    }

  TH2D* hdummy = new TH2D("hdummy","",1,-0.5,60.5,1,0,0.08);
  hdummy->GetYaxis()->SetTitle("v_{2}{FVTXS-FVTXN}");
  hdummy->GetXaxis()->SetTitle("N_{PC1}");
  TGraphErrors* tge_fvtxsfvtxn_200 = new TGraphErrors(61,npc1,v2_fvtxsfvtxn_200,0,ev2_fvtxsfvtxn_200);
  TGraphErrors* tge_fvtxsfvtxn_62 = new TGraphErrors(61,npc1,v2_fvtxsfvtxn_62,0,ev2_fvtxsfvtxn_62);
  TGraphErrors* tge_fvtxsfvtxn_39 = new TGraphErrors(61,npc1,v2_fvtxsfvtxn_39,0,ev2_fvtxsfvtxn_39);
  TGraphErrors* tge_fvtxsfvtxn_20 = new TGraphErrors(61,npc1,v2_fvtxsfvtxn_20,0,ev2_fvtxsfvtxn_20);
  tge_fvtxsfvtxn_200->SetMarkerStyle(kFullCircle);
  tge_fvtxsfvtxn_62->SetMarkerStyle(kFullSquare);
  tge_fvtxsfvtxn_39->SetMarkerStyle(kFullCross);
  tge_fvtxsfvtxn_20->SetMarkerStyle(kFullDiamond);
  tge_fvtxsfvtxn_200->SetMarkerColor(kBlack);
  tge_fvtxsfvtxn_62->SetMarkerColor(kBlue);
  tge_fvtxsfvtxn_39->SetMarkerColor(kRed);
  tge_fvtxsfvtxn_20->SetMarkerColor(kGreen+2);
  hdummy->Draw();
  tge_fvtxsfvtxn_200->Draw("p");
  tge_fvtxsfvtxn_62->Draw("p");
  tge_fvtxsfvtxn_39->Draw("p");
  tge_fvtxsfvtxn_20->Draw("p");
  c1->Print("run16dau_v2_npc1_all.png");

}

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
