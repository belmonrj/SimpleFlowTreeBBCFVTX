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

  cout << "bbc-fvtxn " << float_bbc_fvtxN << endl;
  cout << "bbc-fvtxs " << float_bbc_fvtx << endl;
  cout << "fvtxn-fvtxs " << float_fvtxN_fvtx << endl;

  float reso_bbc_fn = sqrt((float_bbc_fvtxN*float_bbc_fvtx)/float_fvtxN_fvtx); // BNBS/NS
  float reso_fvtx_fn = sqrt((float_fvtxN_fvtx*float_bbc_fvtx)/float_bbc_fvtxN); // NSBS/BN
  float reso_fvtxN_fn = sqrt((float_fvtxN_fvtx*float_bbc_fvtxN)/float_bbc_fvtx); // NSBN/BS

  cout << "bbc resolution is " << reso_bbc_fn << endl;
  cout << "fvtx resolution is " << reso_fvtx_fn << endl;
  cout << "fvtxN resolution is " << reso_fvtxN_fn << endl;

  // ---

  TProfile* tp1f_bbc_cnt = (TProfile*)file->Get("tp1f_reso2_BBC_CNT");
  TProfile* tp1f_cnt_fvtx = (TProfile*)file->Get("tp1f_reso2_CNT_FVTX");
  TProfile* tp1f_cnt_fvtxN = (TProfile*)file->Get("tp1f_reso2_CNT_FVTXN");

  float float_bbc_cnt = tp1f_bbc_cnt->GetBinContent(1);
  float float_cnt_fvtx = tp1f_cnt_fvtx->GetBinContent(1);
  float float_cnt_fvtxN = tp1f_cnt_fvtxN->GetBinContent(1);

  cout << "bbc-cnt " << float_bbc_cnt << endl;
  cout << "cnt-fvtxs " << float_cnt_fvtx << endl;
  cout << "cbt-fvtxn " << float_cnt_fvtxN << endl;

  float reso_bbc = sqrt((float_bbc_cnt*float_bbc_fvtx)/float_cnt_fvtx); // BCBS/CS
  float reso_fvtx = sqrt((float_cnt_fvtx*float_bbc_fvtx)/float_bbc_cnt); // CSBS/BC
  float reso_fvtxN = sqrt((float_cnt_fvtxN*float_bbc_fvtxN)/float_bbc_cnt); // CNBN/BC

  cout << "bbc resolution is " << reso_bbc << endl;
  cout << "fvtx resolution is " << reso_fvtx << endl;
  cout << "fvtxN resolution is " << reso_fvtxN << endl;

  // ---

  float reso_fvtxN_xb = sqrt ( ( float_fvtxN_fvtx * float_cnt_fvtxN ) / float_cnt_fvtx ) ; // NSNC/SC

  cout << "fvtxN resolution is " << reso_fvtxN_xb << endl;

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

  TProfile* hv2_bbcs = (TProfile*)file->Get("bbcs_v2_both_docalib");
  hv2_bbcs->SetLineColor(kRed);
  hv2_bbcs->Scale(1.0/reso_bbc);
  hv2_bbcs->Draw("same");




  TLegend* leg = new TLegend(0.18,0.68,0.38,0.88);
  leg->AddEntry(hv2_fvtxs,"Run16 FVTXS","el");
  leg->AddEntry(hv2_bbcs,"Run16 BBCS","el");
  if ( energy == 200 ) leg->AddEntry(tge_pub,"Run8 (ppg161)","p");
  else leg->AddEntry(tge_pub,"Run8 (200 GeV)","p");
  leg->SetTextSize(0.05);
  leg->Draw();

  TLatex* latex = new TLatex(0.2,0.2,"EP resolutions calculated with CNT");
  latex->SetNDC();
  latex->SetTextSize(0.05);
  //latex->Draw();

  c1->Print(Form("run16dau%d_v2_fvtxsbbcs.pdf",energy));
  c1->Print(Form("run16dau%d_v2_fvtxsbbcs.png",energy));

  TProfile* hv2_fvtxn = (TProfile*)file->Get("fvtxn_v2_both_docalib");
  hv2_fvtxn->SetLineColor(kGreen+2);
  hv2_fvtxn->Scale(1.0/reso_fvtxN);
  hv2_fvtxn->Draw("same");

  delete leg;

  TLegend* leg3 = new TLegend(0.18,0.68,0.38,0.88);
  leg3->AddEntry(hv2_fvtxs,"Run16 FVTXS","el");
  leg3->AddEntry(hv2_fvtxn,"Run16 FVTXN","el");
  leg3->AddEntry(hv2_bbcs,"Run16 BBCS","el");
  if ( energy == 200 ) leg3->AddEntry(tge_pub,"Run8 (ppg161)","p");
  else leg3->AddEntry(tge_pub,"Run8 (200 GeV)","p");
  leg3->SetTextSize(0.05);
  leg3->Draw();

  //latex->Draw();

  c1->Print(Form("run16dau%d_v2_fvtxsnbbcs.pdf",energy));
  c1->Print(Form("run16dau%d_v2_fvtxsnbbcs.png",energy));

  hv2_fvtxn->Scale(reso_fvtxN/reso_fvtxN_fn);
  hv2_fvtxs->Scale(reso_fvtx/reso_fvtx_fn);
  hv2_bbcs->Scale(reso_bbc/reso_bbc_fn);

  delete latex;
  latex = new TLatex(0.2,0.2,"EP resolutions calculated without CNT");
  latex->SetNDC();
  latex->SetTextSize(0.05);
  latex->Draw();

  //  latex->DrawLatex(0.2,0.2,"EP resolutions calculated without CNT");

  c1->Print(Form("run16dau%d_v2_fvtxsnbbcs_fn.pdf",energy));
  c1->Print(Form("run16dau%d_v2_fvtxsnbbcs_fn.png",energy));

  c1->Clear();
  hv2_fvtxs->Draw();
  hv2_bbcs->Draw("same");
  tge_pub->Draw("p");
  leg = new TLegend(0.18,0.68,0.38,0.88);
  leg->AddEntry(hv2_fvtxs,"Run16 FVTXS","el");
  leg->AddEntry(hv2_bbcs,"Run16 BBCS","el");
  if ( energy == 200 ) leg->AddEntry(tge_pub,"Run8 (ppg161)","p");
  else leg->AddEntry(tge_pub,"Run8 (200 GeV)","p");
  leg->SetTextSize(0.05);
  leg->Draw();

  latex->Draw();

  c1->Print(Form("run16dau%d_v2_fvtxsbbcs_fn.pdf",energy));
  c1->Print(Form("run16dau%d_v2_fvtxsbbcs_fn.png",energy));

  hv2_fvtxn->Scale(2*reso_fvtxN_fn/(reso_fvtxN_fn+reso_fvtxN));
  hv2_fvtxs->Scale(2*reso_fvtx_fn/(reso_fvtx_fn+reso_fvtx));
  hv2_bbcs->Scale(2*reso_bbc_fn/(reso_bbc_fn+reso_bbc));

  delete latex;
  latex = new TLatex(0.2,0.2,"Average EP resolutions with and without CNT");
  latex->SetNDC();
  latex->SetTextSize(0.05);
  latex->Draw();

  //  latex->DrawLatex(0.2,0.2,"Average of EP resolutions with and without CNT");

  c1->Print(Form("run16dau%d_v2_fvtxsbbcs_ave.pdf",energy));
  c1->Print(Form("run16dau%d_v2_fvtxsbbcs_ave.png",energy));

  delete leg;
  hv2_fvtxn->Draw("same");
  leg3->Draw();

  //latex->Draw();

  c1->Print(Form("run16dau%d_v2_fvtxsnbbcs_ave.pdf",energy));
  c1->Print(Form("run16dau%d_v2_fvtxsnbbcs_ave.png",energy));

  c1->Clear();

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

