void temp_v3()
{

  TCanvas* c1 = new TCanvas("c1","");

  TFile* file = TFile::Open("input/combined_10.root");
  //TFile* file = TFile::Open("combined.root");
  //TFile* file = TFile::Open("../output/hist_454811.root");

  // ---

  TProfile* tp1f_bbc_cnt = (TProfile*)file->Get("tp1f_reso3_BBC_CNT");
  TProfile* tp1f_bbc_fvtx = (TProfile*)file->Get("tp1f_reso3_BBC_FVTX");
  TProfile* tp1f_cnt_fvtx = (TProfile*)file->Get("tp1f_reso3_CNT_FVTX");

  float float_bbc_cnt = tp1f_bbc_cnt->GetBinContent(1);
  float float_bbc_fvtx = tp1f_bbc_fvtx->GetBinContent(1);
  float float_cnt_fvtx = tp1f_cnt_fvtx->GetBinContent(1);

  float reso_bbc = sqrt((float_bbc_cnt*float_bbc_fvtx)/float_cnt_fvtx);
  float reso_fvtx = sqrt((float_cnt_fvtx*float_bbc_fvtx)/float_bbc_cnt);

  cout << "bbc resolution is " << reso_bbc << endl;
  cout << "fvtx resolution is " << reso_fvtx << endl;

  // ---

  TProfile* hv3_fvtxs = (TProfile*)file->Get("fvtxs_v3_both_docalib");

  hv3_fvtxs->Scale(1.0/reso_fvtx);
  hv3_fvtxs->Draw();
  hv3_fvtxs->SetMaximum(0.15);
  hv3_fvtxs->SetMinimum(0.0);
  hv3_fvtxs->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hv3_fvtxs->GetYaxis()->SetTitle("v_{3}{EP}");

  TLegend *leg = new TLegend(0.18,0.68,0.38,0.88);
  leg->AddEntry(hv3_fvtxs,"Run16 FVTXS","el");
  leg->SetTextSize(0.05);
  leg->Draw();

  c1->Print("run16dau200_v3_fvtxs.pdf");
  c1->Print("run16dau200_v3_fvtxs.png");

  TProfile* hv3_bbcs = (TProfile*)file->Get("bbcs_v3_both_docalib");
  hv3_bbcs->SetLineColor(kRed);
  hv3_bbcs->Scale(1.0/reso_bbc);
  hv3_bbcs->Draw("same");

  delete leg;

  leg = new TLegend(0.18,0.68,0.38,0.88);
  leg->AddEntry(hv3_fvtxs,"Run16 FVTXS","el");
  leg->AddEntry(hv3_bbcs,"Run16 BBCS","el");
  leg->SetTextSize(0.05);
  leg->Draw();

  c1->Print("run16dau200_v3_fvtxsbbcs.pdf");
  c1->Print("run16dau200_v3_fvtxsbbcs.png");

  hv3_fvtxs->SetMinimum(-0.1);
  hv3_fvtxs->SetMaximum(0.1);

  c1->Print("run16dau200_v3_fvtxsbbcs_fullscale.pdf");
  c1->Print("run16dau200_v3_fvtxsbbcs_fullscale.png");


}

