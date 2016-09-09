void temp_npc1()
{

  takename("bbcsfvtxs",2);
  takename("bbcsfvtxn",2);
  takename("fvtxsfvtxn",2);
  takename("bbcsfvtxs",3);
  takename("bbcsfvtxn",3);
  takename("fvtxsfvtxn",3);

}


void takename(const char* name, int harmonic)
{

  TCanvas* c1 = new TCanvas();

  TFile* file_200 = TFile::Open("input/minbiascombined_200.root");
  TFile* file_62 = TFile::Open("input/minbiascombined_62.root");
  TFile* file_39 = TFile::Open("input/minbiascombined_39.root");
  TFile* file_20 = TFile::Open("input/minbiascombined_20.root");

  TProfile* tp1f_200 = (TProfile*)file_200->Get(Form("npc1_os_%s_c%d2",name,harmonic));
  TProfile* tp1f_62 = (TProfile*)file_62->Get(Form("npc1_os_%s_c%d2",name,harmonic));
  TProfile* tp1f_39 = (TProfile*)file_39->Get(Form("npc1_os_%s_c%d2",name,harmonic));
  TProfile* tp1f_20 = (TProfile*)file_20->Get(Form("npc1_os_%s_c%d2",name,harmonic));

  tp1f_200->SetName("tp1f_200");
  tp1f_62->SetName("tp1f_62");
  tp1f_39->SetName("tp1f_39");
  tp1f_20->SetName("tp1f_20");

  TH1D* th1d_vn_200 = sqrt(tp1f_200);
  TH1D* th1d_vn_62 = sqrt(tp1f_62);
  TH1D* th1d_vn_39 = sqrt(tp1f_39);
  TH1D* th1d_vn_20 = sqrt(tp1f_20);

  th1d_vn_200->Draw();
  th1d_vn_200->GetYaxis()->SetTitle(Form("v_{%d}{2}",harmonic));
  th1d_vn_200->GetXaxis()->SetTitle("N_{PC1}");
  c1->Print(Form("FigsMultiplicity/run16dau_v%d_%s_npc1_200.png",harmonic,name));

  th1d_vn_62->Draw();
  th1d_vn_62->GetYaxis()->SetTitle(Form("v_{%d}{2}",harmonic));
  th1d_vn_62->GetXaxis()->SetTitle("N_{PC1}");
  c1->Print(Form("FigsMultiplicity/run16dau_v%d_%s_npc1_62.png",harmonic,name));

  th1d_vn_39->Draw();
  th1d_vn_39->GetYaxis()->SetTitle(Form("v_{%d}{2}",harmonic));
  th1d_vn_39->GetXaxis()->SetTitle("N_{PC1}");
  c1->Print(Form("FigsMultiplicity/run16dau_v%d_%s_npc1_39.png",harmonic,name));

  th1d_vn_20->Draw();
  th1d_vn_20->GetYaxis()->SetTitle(Form("v_{%d}{2}",harmonic));
  th1d_vn_20->GetXaxis()->SetTitle("N_{PC1}");
  c1->Print(Form("FigsMultiplicity/run16dau_v%d_%s_npc1_20.png",harmonic,name));

  double npc1[61];
  double v2_200[61];
  double v2_62[61];
  double v2_39[61];
  double v2_20[61];
  double ev2_200[61];
  double ev2_62[61];
  double ev2_39[61];
  double ev2_20[61];

  for ( int i = 0; i < 61; ++i )
    {
      npc1[i] = th1d_vn_200->GetBinCenter(i+1);
      v2_200[i] = th1d_vn_200->GetBinContent(i+1);
      v2_62[i] = th1d_vn_62->GetBinContent(i+1);
      v2_39[i] = th1d_vn_39->GetBinContent(i+1);
      v2_20[i] = th1d_vn_20->GetBinContent(i+1);
      ev2_200[i] = th1d_vn_200->GetBinError(i+1);
      ev2_62[i] = th1d_vn_62->GetBinError(i+1);
      ev2_39[i] = th1d_vn_39->GetBinError(i+1);
      ev2_20[i] = th1d_vn_20->GetBinError(i+1);
      if ( npc1[i] > 61 ) v2_200[i] = -9;
      if ( npc1[i] > 28 ) v2_62[i] = -9;
      if ( npc1[i] > 19 ) v2_39[i] = -9;
      if ( npc1[i] > 10 ) v2_20[i] = -9;
    }

  TH2D* hdummy = new TH2D("hdummy","",1,-0.5,60.5,1,0,0.08);
  hdummy->GetYaxis()->SetTitle(Form("v_{%d}{2}",harmonic));
  hdummy->GetXaxis()->SetTitle("N_{PC1}");
  TGraphErrors* tge_200 = new TGraphErrors(61,npc1,v2_200,0,ev2_200);
  TGraphErrors* tge_62 = new TGraphErrors(61,npc1,v2_62,0,ev2_62);
  TGraphErrors* tge_39 = new TGraphErrors(61,npc1,v2_39,0,ev2_39);
  TGraphErrors* tge_20 = new TGraphErrors(61,npc1,v2_20,0,ev2_20);
  tge_200->SetMarkerStyle(kFullCircle);
  tge_62->SetMarkerStyle(kFullSquare);
  tge_39->SetMarkerStyle(kFullCross);
  tge_20->SetMarkerStyle(kFullDiamond);
  tge_200->SetMarkerColor(kBlack);
  tge_62->SetMarkerColor(kBlue);
  tge_39->SetMarkerColor(kRed);
  tge_20->SetMarkerColor(kGreen+2);
  hdummy->Draw();
  tge_200->Draw("p");
  tge_62->Draw("p");
  tge_39->Draw("p");
  tge_20->Draw("p");
  TLegend* leg = new TLegend(0.68,0.68,0.88,0.88);
  leg->AddEntry(tge_200,"200 GeV","p");
  leg->AddEntry(tge_62,"62.4 GeV","p");
  leg->AddEntry(tge_39,"39 GeV","p");
  leg->AddEntry(tge_20,"19.6 GeV","p");
  leg->SetTextSize(0.045);
  leg->SetFillStyle(0);
  leg->Draw();
  c1->Print(Form("FigsMultiplicity/run16dau_v%d_%s_npc1_all.png",harmonic,name));

  // ---

  TFile* file_200C = TFile::Open("input/combined_200.root");
  TFile* file_62C = TFile::Open("input/combined_62.root");
  TFile* file_39C = TFile::Open("input/combined_39.root");
  TFile* file_20C = TFile::Open("input/combined_20.root");

  TProfile* tp1f_200C = (TProfile*)file_200C->Get(Form("npc1_os_%s_c%d2",name,harmonic));
  TProfile* tp1f_62C = (TProfile*)file_62C->Get(Form("npc1_os_%s_c%d2",name,harmonic));
  TProfile* tp1f_39C = (TProfile*)file_39C->Get(Form("npc1_os_%s_c%d2",name,harmonic));
  TProfile* tp1f_20C = (TProfile*)file_20C->Get(Form("npc1_os_%s_c%d2",name,harmonic));

  tp1f_200C->SetName("tp1f_200C");
  tp1f_62C->SetName("tp1f_62C");
  tp1f_39C->SetName("tp1f_39C");
  tp1f_20C->SetName("tp1f_20C");

  TH1D* th1d_vn_200C = sqrt(tp1f_200C);
  TH1D* th1d_vn_62C = sqrt(tp1f_62C);
  TH1D* th1d_vn_39C = sqrt(tp1f_39C);
  TH1D* th1d_vn_20C = sqrt(tp1f_20C);

  // ---

  th1d_vn_200C->SetLineColor(kRed);
  th1d_vn_200->SetLineColor(kBlack);
  th1d_vn_200->Draw();
  th1d_vn_200->SetMinimum(0);
  th1d_vn_200->SetMaximum(0.08);
  th1d_vn_200C->Draw("same");
  c1->Print(Form("FigsMultiplicity/run16dau_cmbc_v%d_%s_npc1_200.png",harmonic,name));

  th1d_vn_62C->SetLineColor(kRed);
  th1d_vn_62->SetLineColor(kBlack);
  th1d_vn_62->Draw();
  th1d_vn_62->SetMinimum(0);
  th1d_vn_62->SetMaximum(0.08);
  th1d_vn_62->GetXaxis()->SetRangeUser(-0.5,30.5);
  th1d_vn_62C->Draw("same");
  c1->Print(Form("FigsMultiplicity/run16dau_cmbc_v%d_%s_npc1_62.png",harmonic,name));

  th1d_vn_39C->SetLineColor(kRed);
  th1d_vn_39->SetLineColor(kBlack);
  th1d_vn_39->Draw();
  th1d_vn_39->SetMinimum(0);
  th1d_vn_39->SetMaximum(0.08);
  th1d_vn_39->GetXaxis()->SetRangeUser(-0.5,20.5);
  th1d_vn_39C->Draw("same");
  c1->Print(Form("FigsMultiplicity/run16dau_cmbc_v%d_%s_npc1_39.png",harmonic,name));

  th1d_vn_20C->SetLineColor(kRed);
  th1d_vn_20->SetLineColor(kBlack);
  th1d_vn_20->Draw();
  th1d_vn_20->SetMinimum(0);
  th1d_vn_20->SetMaximum(0.08);
  th1d_vn_20->GetXaxis()->SetRangeUser(-0.5,10.5);
  th1d_vn_20C->Draw("same");
  c1->Print(Form("FigsMultiplicity/run16dau_cmbc_v%d_%s_npc1_20.png",harmonic,name));

}

TH1D* sqrt(TProfile* tp1f)
{
  TH1D* th1d = tp1f->ProjectionX(Form("%s_px",tp1f->GetName()));
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
