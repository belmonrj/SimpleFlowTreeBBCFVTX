void doenergy(int);

void temp_v2v4()
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

  TProfile* tp1f_bbc_fvtxN = (TProfile*)file->Get("tp1f_reso2_BBC_FVTXN");
  TProfile* tp1f_bbc_fvtx = (TProfile*)file->Get("tp1f_reso2_BBC_FVTX");
  TProfile* tp1f_fvtxN_fvtx = (TProfile*)file->Get("tp1f_reso2_FVTXS_FVTXN");

  float float_bbc_fvtxN = tp1f_bbc_fvtxN->GetBinContent(1);
  float float_bbc_fvtx = tp1f_bbc_fvtx->GetBinContent(1);
  float float_fvtxN_fvtx = tp1f_fvtxN_fvtx->GetBinContent(1);

  cout <<  float_bbc_fvtxN << endl;
  cout <<  float_bbc_fvtx << endl;
  cout <<  float_fvtxN_fvtx << endl;

  float reso_bbc = sqrt((float_bbc_fvtxN*float_bbc_fvtx)/float_fvtxN_fvtx);
  float reso_fvtx = sqrt((float_fvtxN_fvtx*float_bbc_fvtx)/float_bbc_fvtxN);

  cout << "bbc resolution is " << reso_bbc << endl;
  cout << "fvtx resolution is " << reso_fvtx << endl;

  // ---

  TProfile* tp1f_bbc_cnt = (TProfile*)file->Get("tp1f_reso2_BBC_CNT");
  tp1f_bbc_fvtx = (TProfile*)file->Get("tp1f_reso2_BBC_FVTX");
  TProfile* tp1f_cnt_fvtx = (TProfile*)file->Get("tp1f_reso2_CNT_FVTX");

  float float_bbc_cnt = tp1f_bbc_cnt->GetBinContent(1);
  float_bbc_fvtx = tp1f_bbc_fvtx->GetBinContent(1);
  float float_cnt_fvtx = tp1f_cnt_fvtx->GetBinContent(1);

  cout <<  float_bbc_cnt << endl;
  cout <<  float_bbc_fvtx << endl;
  cout <<  float_cnt_fvtx << endl;

  reso_bbc = sqrt((float_bbc_cnt*float_bbc_fvtx)/float_cnt_fvtx);
  reso_fvtx = sqrt((float_cnt_fvtx*float_bbc_fvtx)/float_bbc_cnt);

  cout << "bbc resolution is " << reso_bbc << endl;
  cout << "fvtx resolution is " << reso_fvtx << endl;

  // ---

  TProfile* hv2_fvtxs = (TProfile*)file->Get("fvtxs_v2_both_docalib");

  hv2_fvtxs->Scale(1.0/reso_fvtx);
  hv2_fvtxs->Draw();
  // the 62 GeV is actually 62.4 and the 20 GeV is actually 19.6, so need to modify
  hv2_fvtxs->SetTitle(Form("d+Au collisions at #sqrt{s_{NN}} = %d GeV",energy));
  if ( energy == 62 ) hv2_fvtxs->SetTitle(Form("d+Au collisions at #sqrt{s_{NN}} = 62.4 GeV"));
  if ( energy == 20 ) hv2_fvtxs->SetTitle(Form("d+Au collisions at #sqrt{s_{NN}} = 19.6 GeV"));
  hv2_fvtxs->SetMaximum(0.17);
  hv2_fvtxs->SetMinimum(0.0);
  hv2_fvtxs->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hv2_fvtxs->GetYaxis()->SetTitle("v_{2}{EP}");
  hv2_fvtxs->GetYaxis()->SetTitleOffset(1.25);

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
  leg->AddEntry(tge_pub,"Run8 (ppg161)","p");
  leg->SetTextSize(0.05);
  leg->Draw();

  TProfile* hv2_bbcs = (TProfile*)file->Get("bbcs_v2_both_docalib");
  hv2_bbcs->SetLineColor(kRed);
  hv2_bbcs->Scale(1.0/reso_bbc);
  hv2_bbcs->Draw("same");

  delete leg;

  leg = new TLegend(0.18,0.68,0.38,0.88);
  leg->AddEntry(hv2_fvtxs,"Run16 FVTXS","el");
  leg->AddEntry(hv2_bbcs,"Run16 BBCS","el");
  leg->AddEntry(tge_pub,"Run8 (ppg161)","p");
  leg->SetTextSize(0.05);
  leg->Draw();

  c1->Print(Form("run16dau%d_v2_fvtxsbbcs.pdf",energy));
  c1->Print(Form("run16dau%d_v2_fvtxsbbcs.png",energy));


  // ---
  // --- 4Psi2
  // ---

  TProfile* hv4_4Psi2_fvtxs = (TProfile*)file->Get("fvtxs_v4_4Psi2_both_docalib");

  hv4_4Psi2_fvtxs->Scale(1.0/reso_fvtx);
  hv4_4Psi2_fvtxs->Draw();
  // the 62 GeV is actually 62.4 and the 20 GeV is actually 19.6, so need to modify
  hv4_4Psi2_fvtxs->SetTitle(Form("d+Au collisions at #sqrt{s_{NN}} = %d GeV",energy));
  if ( energy == 62 ) hv4_4Psi2_fvtxs->SetTitle(Form("d+Au collisions at #sqrt{s_{NN}} = 62.4 GeV"));
  if ( energy == 20 ) hv4_4Psi2_fvtxs->SetTitle(Form("d+Au collisions at #sqrt{s_{NN}} = 19.6 GeV"));
  hv4_4Psi2_fvtxs->SetMaximum(0.05);
  hv4_4Psi2_fvtxs->SetMinimum(-0.05);
  hv4_4Psi2_fvtxs->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hv4_4Psi2_fvtxs->GetYaxis()->SetTitle("v_{4}{#Psi_{2}}");
  hv4_4Psi2_fvtxs->GetYaxis()->SetTitleOffset(1.25);

  if ( leg ) delete leg;

  TLegend *leg = new TLegend(0.18,0.68,0.38,0.88);
  leg->AddEntry(hv4_4Psi2_fvtxs,"Run16 FVTXS","el");
  leg->SetTextSize(0.05);
  leg->Draw();

  TProfile* hv4_4Psi2_bbcs = (TProfile*)file->Get("bbcs_v4_4Psi2_both_docalib");
  hv4_4Psi2_bbcs->SetLineColor(kRed);
  //hv4_4Psi2_bbcs->Scale(1.0/0.104519);
  hv4_4Psi2_bbcs->Scale(1.0/reso_bbc);
  hv4_4Psi2_bbcs->Draw("same");

  delete leg;

  leg = new TLegend(0.18,0.68,0.38,0.88);
  leg->AddEntry(hv4_4Psi2_fvtxs,"Run16 FVTXS","el");
  leg->AddEntry(hv4_4Psi2_bbcs,"Run16 BBCS","el");
  leg->SetTextSize(0.05);
  leg->Draw();

  TLine line(0.0,0.0,3.0,0.0);
  line.SetLineStyle(2);
  line.SetLineWidth(2);
  line.Draw();

  c1->Print(Form("run16dau%d_v4_4Psi2_fvtxsbbcs.pdf",energy));
  c1->Print(Form("run16dau%d_v4_4Psi2_fvtxsbbcs.png",energy));



  delete c1;

}

