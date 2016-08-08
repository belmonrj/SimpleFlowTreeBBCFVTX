void doenergy(int);

void temp_split_v2()
{

  doenergy(200);
  doenergy(62);
  doenergy(39);
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
  ofstream fout((const char*)Form("v2fvtxs_%d.dat",energy));
  for ( int i = 0; i < hv2_fvtxs->GetNbinsX(); ++i ) fout << hv2_fvtxs->GetBinCenter(i+1) << " " << hv2_fvtxs->GetBinContent(i+1) << " " << hv2_fvtxs->GetBinError(i+1) << " " << endl;
  fout.close();
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

  TLegend *leg = new TLegend(0.18,0.72,0.38,0.88);
  leg->AddEntry(tge_pub,"Run8 200 GeV (ppg161)","p");
  leg->AddEntry(hv2_fvtxs,"Run16 FVTXS","el");
  leg->SetTextSize(0.045);
  leg->Draw();
  float tlx = 0.55;
  float tly = 0.23;
  TLatex* tlref = new TLatex(tlx,tly,Form("FVTXS Res = %.3f",reso_fvtx));
  tlref->SetTextSize(0.05);
  tlref->SetNDC();
  tlref->Draw();
  hv2_fvtxs->SetMaximum(0.17);
  if ( energy <= 39 ) hv2_fvtxs->SetMaximum(0.35);
  c1->Print(Form("FigsHarmonicCoefficient/run16dau%d_v2_fvtxs.pdf",energy));
  c1->Print(Form("FigsHarmonicCoefficient/run16dau%d_v2_fvtxs.png",energy));

  TProfile* hv2_bbcs = (TProfile*)file->Get("bbcs_v2_both_docalib");
  hv2_bbcs->SetLineColor(kRed);
  hv2_bbcs->Scale(1.0/reso_bbc);
  hv2_bbcs->Draw("same");
  delete leg;
  leg = new TLegend(0.18,0.72,0.38,0.88);
  leg->AddEntry(tge_pub,"Run8 200 GeV (ppg161)","p");
  leg->AddEntry(hv2_fvtxs,"Run16 FVTXS","el");
  leg->AddEntry(hv2_bbcs,"Run16 BBCS","el");
  leg->SetTextSize(0.045);
  leg->Draw();
  float tlx2 = 0.55;
  float tly2 = 0.18;
  TLatex* tlref2 = new TLatex(tlx2,tly2,Form("BBCS Res = %.3f",reso_bbc));
  tlref2->SetTextSize(0.05);
  tlref2->SetNDC();
  tlref2->Draw();
  c1->Print(Form("FigsHarmonicCoefficient/run16dau%d_v2_fvtxsbbcs.pdf",energy));
  c1->Print(Form("FigsHarmonicCoefficient/run16dau%d_v2_fvtxsbbcs.png",energy));

  c1->Clear();

  delete c1;

  c1 = new TCanvas("c1","",700,900);
  c1->Divide(1,2);

  c1->cd(1);

  //gPad->SetTopMargin(0.05);
  gPad->SetBottomMargin(0.0);
  gPad->SetRightMargin(0.0);
  gPad->SetPad(.005, .3, .98, .98 );

  hv2_fvtxs->SetMinimum(0.00001);
  hv2_fvtxs->GetXaxis()->SetLabelSize(0.0);
  hv2_fvtxs->GetXaxis()->SetTitleSize(0.0);
  hv2_fvtxs->Draw();
  hv2_bbcs->Draw("same");
  if ( energy == 39 ) hv2_fvtxs->SetMaximum(0.5);
  tge_pub->Draw("p");
  leg->Draw();

  c1->cd(2);

  gPad->SetPad(.005, .005, .98, .3);
  gPad->SetTopMargin(0.0);
  gPad->SetRightMargin(0.0);
  gPad->SetBottomMargin(0.3);

  TH1D* hpv2_fvtxs = hv2_fvtxs->ProjectionX("hpv2_fvtxs");
  TH1D* hpv2_bbcs = hv2_bbcs->ProjectionX("hpv2_bbcs");
  hpv2_fvtxs->Divide(hpv2_bbcs);
  hpv2_fvtxs->SetLineColor(kBlack);
  hpv2_fvtxs->Draw();
  hpv2_fvtxs->SetTitle("");
  hpv2_fvtxs->SetMinimum(0.75);
  hpv2_fvtxs->SetMaximum(1.25);
  if ( energy == 62 ) hpv2_fvtxs->SetMinimum(0.5);
  if ( energy == 62 ) hpv2_fvtxs->SetMaximum(1.5);
  if ( energy == 39 ) hpv2_fvtxs->SetMinimum(0.0);
  if ( energy == 39 ) hpv2_fvtxs->SetMaximum(6.0);
  hpv2_fvtxs->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hpv2_fvtxs->GetXaxis()->SetTitleSize(0.10);
  hpv2_fvtxs->GetXaxis()->SetLabelSize(0.10);
  hpv2_fvtxs->GetYaxis()->SetTitle("FVTXS/BBC");
  hpv2_fvtxs->GetYaxis()->SetTitleOffset(0.6); // very surpised this is the right number???
  hpv2_fvtxs->GetYaxis()->SetTitleSize(0.10);
  hpv2_fvtxs->GetYaxis()->SetLabelSize(0.10);
  hpv2_fvtxs->GetYaxis()->SetNdivisions(505);

  TLine* line1 = new TLine(0.0,1.0,3.0,1.0);
  TLine* line2 = new TLine(0.0,1.15,3.0,1.15);
  line1->SetLineStyle(2);
  line2->SetLineStyle(2);
  line1->Draw();
  line2->Draw();

  c1->Print(Form("FigsHarmonicCoefficient/run16dau%d_v2_fvtxsbbcs_split.pdf",energy));
  c1->Print(Form("FigsHarmonicCoefficient/run16dau%d_v2_fvtxsbbcs_split.png",energy));

  delete c1;

}

