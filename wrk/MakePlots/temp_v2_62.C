void temp_v2_62()
{

  TCanvas* c1 = new TCanvas("c1","");

  //TFile* file = TFile::Open("input/combined_10.root");
  //TFile* file = TFile::Open("input/combined_F5.root");
  //TFile* file = TFile::Open("input/combined_S5.root");
  TFile* file = TFile::Open("input/hist_455795.root");
  //TFile* file = TFile::Open("../output/hist_454811.root");

  // ---

  TProfile* tp1f_bbc_cnt = (TProfile*)file->Get("tp1f_reso2_BBC_CNT");
  TProfile* tp1f_bbc_fvtx = (TProfile*)file->Get("tp1f_reso2_BBC_FVTX");
  TProfile* tp1f_cnt_fvtx = (TProfile*)file->Get("tp1f_reso2_CNT_FVTX");

  float float_bbc_cnt = tp1f_bbc_cnt->GetBinContent(1);
  float float_bbc_fvtx = tp1f_bbc_fvtx->GetBinContent(1);
  float float_cnt_fvtx = tp1f_cnt_fvtx->GetBinContent(1);

  float reso_bbc = sqrt((float_bbc_cnt*float_bbc_fvtx)/float_cnt_fvtx);
  float reso_fvtx = sqrt((float_cnt_fvtx*float_bbc_fvtx)/float_bbc_cnt);

  cout << "bbc resolution is " << reso_bbc << endl;
  cout << "fvtx resolution is " << reso_fvtx << endl;

  // ---

  TProfile* hv2_fvtxs = (TProfile*)file->Get("fvtxs_v2_both_docalib");

  //hv2_fvtxs->Scale(1.0/0.237392);
  hv2_fvtxs->Scale(1.0/reso_fvtx);
  hv2_fvtxs->Draw();
  hv2_fvtxs->SetMaximum(0.15);
  hv2_fvtxs->SetMinimum(0.0);
  hv2_fvtxs->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hv2_fvtxs->GetYaxis()->SetTitle("v_{2}{EP}");

  ifstream finpub("ppg161.dat");
  float pt[13], pubv2[13], epubv2[13], esyspubv2[13];
  for ( int i = 0; i < 13; ++i )
    {
      finpub>>pt[i]>>pubv2[i]>>epubv2[i]>>esyspubv2[i];
    }

  TGraphErrors* tge_pub = new TGraphErrors(13,pt,pubv2,0,epubv2);
  tge_pub->SetMarkerStyle(kFullCircle);
  tge_pub->Draw("p");

  TLegend *leg = new TLegend(0.18,0.68,0.38,0.88);
  leg->AddEntry(hv2_fvtxs,"Run16 FVTXS","el");
  leg->AddEntry(tge_pub,"PPG161","p");
  leg->SetTextSize(0.05);
  leg->Draw();

  c1->Print("run16dau62_v2_fvtxs.pdf");
  c1->Print("run16dau62_v2_fvtxs.png");

  TProfile* hv2_bbcs = (TProfile*)file->Get("bbcs_v2_both_docalib");
  hv2_bbcs->SetLineColor(kRed);
  //hv2_bbcs->Scale(1.0/0.104519);
  hv2_bbcs->Scale(1.0/reso_bbc);
  hv2_bbcs->Draw("same");

  delete leg;

  leg = new TLegend(0.18,0.68,0.38,0.88);
  leg->AddEntry(hv2_fvtxs,"Run16 FVTXS","el");
  leg->AddEntry(hv2_bbcs,"Run16 BBCS","el");
  leg->AddEntry(tge_pub,"PPG161","p");
  leg->SetTextSize(0.05);
  leg->Draw();

  c1->Print("run16dau62_v2_fvtxsbbcs.pdf");
  c1->Print("run16dau62_v2_fvtxsbbcs.png");

  hv2_fvtxs->SetMaximum(1.0);

  c1->Print("run16dau62_v2_fvtxsbbcs_fullscale.pdf");
  c1->Print("run16dau62_v2_fvtxsbbcs_fullscale.png");


}

