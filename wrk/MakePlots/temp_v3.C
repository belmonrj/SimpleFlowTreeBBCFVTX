void doenergy(int);
void doenergy_eastwest(int);


void temp_v3()
{

  doenergy(200);
  doenergy(62);
  doenergy(20);

}


void doenergy(int energy)
{

  gStyle->SetOptTitle(1);

  TCanvas* c1 = new TCanvas("c1","");

  TFile* file = TFile::Open(Form("input/combined_%d.root",energy));

  // ---

  TProfile* tp1f_bbc_fvtxN = (TProfile*)file->Get("tp1f_reso3_BBC_FVTXN");
  TProfile* tp1f_bbc_fvtx = (TProfile*)file->Get("tp1f_reso3_BBC_FVTX");
  TProfile* tp1f_fvtxN_fvtx = (TProfile*)file->Get("tp1f_reso3_FVTXS_FVTXN");

  float float_bbc_fvtxN = tp1f_bbc_fvtxN->GetBinContent(1);
  float float_bbc_fvtx = tp1f_bbc_fvtx->GetBinContent(1);
  float float_fvtxN_fvtx = tp1f_fvtxN_fvtx->GetBinContent(1);

  cout <<  float_bbc_fvtxN << endl;
  cout <<  float_bbc_fvtx << endl;
  cout <<  float_fvtxN_fvtx << endl;

  float reso_bbc = sqrt(fabs((float_bbc_fvtxN*float_bbc_fvtx)/float_fvtxN_fvtx));
  float reso_fvtx = sqrt(fabs((float_fvtxN_fvtx*float_bbc_fvtx)/float_bbc_fvtxN));

  cout << "bbc resolution is " << reso_bbc << endl;
  cout << "fvtx resolution is " << reso_fvtx << endl;

  // ---

  TProfile* tp1f_bbc_cnt = (TProfile*)file->Get("tp1f_reso3_BBC_CNT");
  tp1f_bbc_fvtx = (TProfile*)file->Get("tp1f_reso3_BBC_FVTX");
  TProfile* tp1f_cnt_fvtx = (TProfile*)file->Get("tp1f_reso3_CNT_FVTX");

  float float_bbc_cnt = tp1f_bbc_cnt->GetBinContent(1);
  float_bbc_fvtx = tp1f_bbc_fvtx->GetBinContent(1);
  float float_cnt_fvtx = tp1f_cnt_fvtx->GetBinContent(1);

  cout <<  float_bbc_cnt << endl;
  cout <<  float_bbc_fvtx << endl;
  cout <<  float_cnt_fvtx << endl;

  reso_bbc = sqrt(fabs((float_bbc_cnt*float_bbc_fvtx)/float_cnt_fvtx));
  reso_fvtx = sqrt(fabs((float_cnt_fvtx*float_bbc_fvtx)/float_bbc_cnt));

  cout << "bbc resolution is " << reso_bbc << endl;
  cout << "fvtx resolution is " << reso_fvtx << endl;

  // ---

  TProfile* hv3_fvtxs = (TProfile*)file->Get("fvtxs_v3_both_docalib");

  hv3_fvtxs->Scale(1.0/reso_fvtx);
  hv3_fvtxs->Draw();
  hv3_fvtxs->SetTitle(Form("d+Au collisions at #sqrt{s_{NN}} = %d GeV",energy));
  if ( energy == 62 ) hv3_fvtxs->SetTitle(Form("d+Au collisions at #sqrt{s_{NN}} = 62.4 GeV"));
  if ( energy == 20 ) hv3_fvtxs->SetTitle(Form("d+Au collisions at #sqrt{s_{NN}} = 19.6 GeV"));
  hv3_fvtxs->SetMaximum(0.15);
  hv3_fvtxs->SetMinimum(0.0);
  hv3_fvtxs->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hv3_fvtxs->GetYaxis()->SetTitle("v_{3}{EP}");
  hv3_fvtxs->GetYaxis()->SetTitleOffset(1.25);

  TLegend *leg = new TLegend(0.18,0.68,0.38,0.88);
  leg->AddEntry(hv3_fvtxs,"Run16 FVTXS","el");
  leg->SetTextSize(0.05);
  leg->Draw();

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

  // c1->Print("run16dau200_v3_fvtxsbbcs.pdf");
  // c1->Print("run16dau200_v3_fvtxsbbcs.png");

  hv3_fvtxs->SetMinimum(-0.05);
  hv3_fvtxs->SetMaximum(0.15);

  TLine line(0.0,0.0,3.0,0.0);
  line.SetLineStyle(2);
  line.SetLineWidth(2);
  line.Draw();

  // c1->Print("run16dau200_v3_fvtxsbbcs_fullscale.pdf");
  // c1->Print("run16dau200_v3_fvtxsbbcs_fullscale.png");
  c1->Print(Form("run16dau%d_v3_fvtxsbbcs.pdf",energy));
  c1->Print(Form("run16dau%d_v3_fvtxsbbcs.png",energy));

  delete c1;

}


void doenergy_eastwest(int energy)
{

  cout << "East-West code non-functional for now, will revisit later" << endl;
  return;

  TCanvas* c1 = new TCanvas();
  // ---

  TProfile* hv3_fvtxs_east = (TProfile*)file->Get("fvtxs_v3_east_docalib");
  TProfile* hv3_fvtxs_west = (TProfile*)file->Get("fvtxs_v3_west_docalib");
  hv3_fvtxs_east->Scale(1.0/reso_fvtx);
  hv3_fvtxs_west->Scale(1.0/reso_fvtx);
  hv3_fvtxs->SetLineColor(kBlack);
  hv3_fvtxs_east->SetLineColor(kBlue);
  hv3_fvtxs_west->SetLineColor(kRed);
  hv3_fvtxs->SetMinimum(-1);
  hv3_fvtxs->SetMaximum(1);
  hv3_fvtxs->Draw();
  hv3_fvtxs_east->Draw("same");
  hv3_fvtxs_west->Draw("same");

  TProfile* hv3_bbcs_east = (TProfile*)file->Get("bbcs_v3_east_docalib");
  TProfile* hv3_bbcs_west = (TProfile*)file->Get("bbcs_v3_west_docalib");
  hv3_bbcs_east->Scale(1.0/reso_bbc);
  hv3_bbcs_west->Scale(1.0/reso_bbc);
  hv3_bbcs->SetLineColor(kBlack);
  hv3_bbcs_east->SetLineColor(kBlue);
  hv3_bbcs_west->SetLineColor(kRed);
  hv3_bbcs->SetMinimum(-1);
  hv3_bbcs->SetMaximum(1);
  hv3_bbcs->Draw("same");
  hv3_bbcs_east->Draw("same");
  hv3_bbcs_west->Draw("same");

  c1->Print("run16dau200_v3_fvtxsbbcs_fullscale_eastwest.pdf");
  c1->Print("run16dau200_v3_fvtxsbbcs_fullscale_eastwest.png");


}

