void doenergy(int, int);
void doenergy3(int, int);

TFile* outfile;

void show_layers_vn()
{

  doenergy(200,2);
  doenergy(200,3);
  doenergy3(200,2);
  doenergy3(200,3);

}

void doenergy(int energy, int harmonic)
{

  gStyle->SetOptTitle(1);

  TCanvas* c1 = new TCanvas("c1","");

  TFile* file = TFile::Open(Form("input/combined_%d.root",energy));

  // ---

  TProfile* tp1f_bbc_cnt = (TProfile*)file->Get(Form("tp1f_reso%d_BBC_CNT",harmonic));
  TProfile* tp1f_bbc_fvtx = (TProfile*)file->Get(Form("tp1f_reso%d_BBC_FVTX",harmonic));
  TProfile* tp1f_cnt_fvtx = (TProfile*)file->Get(Form("tp1f_reso%d_CNT_FVTX",harmonic));
  TProfile* tp1f_bbc_fvtx0 = (TProfile*)file->Get(Form("tp1f_reso%d_BBC_FVTX0",harmonic));
  TProfile* tp1f_cnt_fvtx0 = (TProfile*)file->Get(Form("tp1f_reso%d_CNT_FVTX0",harmonic));
  TProfile* tp1f_bbc_fvtx1 = (TProfile*)file->Get(Form("tp1f_reso%d_BBC_FVTX1",harmonic));
  TProfile* tp1f_cnt_fvtx1 = (TProfile*)file->Get(Form("tp1f_reso%d_CNT_FVTX1",harmonic));
  TProfile* tp1f_bbc_fvtx2 = (TProfile*)file->Get(Form("tp1f_reso%d_BBC_FVTX2",harmonic));
  TProfile* tp1f_cnt_fvtx2 = (TProfile*)file->Get(Form("tp1f_reso%d_CNT_FVTX2",harmonic));
  TProfile* tp1f_bbc_fvtx3 = (TProfile*)file->Get(Form("tp1f_reso%d_BBC_FVTX3",harmonic));
  TProfile* tp1f_cnt_fvtx3 = (TProfile*)file->Get(Form("tp1f_reso%d_CNT_FVTX3",harmonic));

  float float_bbc_cnt = tp1f_bbc_cnt->GetBinContent(1);
  float float_bbc_fvtx = tp1f_bbc_fvtx->GetBinContent(1);
  float float_cnt_fvtx = tp1f_cnt_fvtx->GetBinContent(1);
  float float_bbc_fvtx0 = tp1f_bbc_fvtx0->GetBinContent(1);
  float float_cnt_fvtx0 = tp1f_cnt_fvtx0->GetBinContent(1);
  float float_bbc_fvtx1 = tp1f_bbc_fvtx1->GetBinContent(1);
  float float_cnt_fvtx1 = tp1f_cnt_fvtx1->GetBinContent(1);
  float float_bbc_fvtx2 = tp1f_bbc_fvtx2->GetBinContent(1);
  float float_cnt_fvtx2 = tp1f_cnt_fvtx2->GetBinContent(1);
  float float_bbc_fvtx3 = tp1f_bbc_fvtx3->GetBinContent(1);
  float float_cnt_fvtx3 = tp1f_cnt_fvtx3->GetBinContent(1);

  cout << "bbc-cnt " << float_bbc_cnt << endl;
  cout << "bbc-fvtxs " << float_bbc_fvtx << endl;
  cout << "cnt-fvtxs " << float_cnt_fvtx << endl;


  float reso_bbc = sqrt((float_bbc_cnt*float_bbc_fvtx)/float_cnt_fvtx); // BCBS/CS
  float reso_fvtx = sqrt((float_cnt_fvtx*float_bbc_fvtx)/float_bbc_cnt); // CSBS/BC
  float reso_fvtx0 = sqrt((float_cnt_fvtx0*float_bbc_fvtx0)/float_bbc_cnt); // CSBS/BC
  float reso_fvtx1 = sqrt((float_cnt_fvtx1*float_bbc_fvtx1)/float_bbc_cnt); // CSBS/BC
  float reso_fvtx2 = sqrt((float_cnt_fvtx2*float_bbc_fvtx2)/float_bbc_cnt); // CSBS/BC
  float reso_fvtx3 = sqrt((float_cnt_fvtx3*float_bbc_fvtx3)/float_bbc_cnt); // CSBS/BC

  cout << "bbc resolution is " << reso_bbc << endl;
  cout << "fvtx resolution is " << reso_fvtx << endl;
  cout << "fvtx0 resolution is " << reso_fvtx0 << endl;
  cout << "fvtx1 resolution is " << reso_fvtx1 << endl;
  cout << "fvtx2 resolution is " << reso_fvtx2 << endl;
  cout << "fvtx3 resolution is " << reso_fvtx3 << endl;

  // ---

  //return;

  TProfile* hvn_bbcs_B = (TProfile*)file->Get(Form("bbcs_v%d_both_docalib",harmonic));
  hvn_bbcs_B->SetMarkerStyle(kOpenSquare);
  hvn_bbcs_B->SetMarkerColor(kBlack);
  hvn_bbcs_B->SetLineColor(kBlack);

  TProfile* hvn_fvtxs_B = (TProfile*)file->Get(Form("fvtxs_v%d_both_docalib",harmonic));
  TProfile* hvn_fvtxs0_B = (TProfile*)file->Get(Form("fvtxs0_v%d_both_docalib",harmonic));
  TProfile* hvn_fvtxs1_B = (TProfile*)file->Get(Form("fvtxs1_v%d_both_docalib",harmonic));
  TProfile* hvn_fvtxs2_B = (TProfile*)file->Get(Form("fvtxs2_v%d_both_docalib",harmonic));
  TProfile* hvn_fvtxs3_B = (TProfile*)file->Get(Form("fvtxs3_v%d_both_docalib",harmonic));

  hvn_fvtxs_B->SetMarkerStyle(kOpenCircle);
  hvn_fvtxs_B->SetMarkerColor(kBlack);
  hvn_fvtxs_B->SetLineColor(kBlack);
  hvn_fvtxs0_B->SetLineColor(kBlue);
  hvn_fvtxs1_B->SetLineColor(kRed);
  hvn_fvtxs2_B->SetLineColor(kGreen+2);
  hvn_fvtxs3_B->SetLineColor(kMagenta+2);

  double MIN = 0;
  double MAX = 1;
  if ( harmonic == 2 ) {MIN = -0.01; MAX = 0.05;}
  if ( harmonic == 3 ) {MIN = -0.002; MAX = 0.005;}
  TH2D* h2dummy = new TH2D("h2dummy","",1,0.0,3.0,1,MIN,MAX);
  h2dummy->Draw();
  TLine line(0,0,3,0);
  line.SetLineStyle(2);
  line.SetLineWidth(2);
  line.Draw();
  h2dummy->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h2dummy->GetYaxis()->SetTitle(Form("v_{%d} not corrected for EP resolution",harmonic));
  h2dummy->GetYaxis()->SetTitleOffset(1.3);

  hvn_fvtxs_B->Draw("same");
  hvn_fvtxs0_B->Draw("same");
  hvn_fvtxs1_B->Draw("same");
  hvn_fvtxs2_B->Draw("same");
  hvn_fvtxs3_B->Draw("same");

  TLegend* leg = new TLegend(0.18,0.68,0.38,0.88);
  leg->SetHeader(Form("%d GeV, FVTXS",energy));
  leg->AddEntry(hvn_fvtxs_B,"all layers","elp");
  leg->AddEntry(hvn_fvtxs0_B,"layer 0","el");
  leg->AddEntry(hvn_fvtxs1_B,"layer 1","el");
  leg->AddEntry(hvn_fvtxs2_B,"layer 2","el");
  leg->AddEntry(hvn_fvtxs3_B,"layer 3","el");
  leg->SetTextSize(0.05);
  leg->Draw();

  c1->Print(Form("FigsCheckThree/special_fvtxsL_B_energy%d_harm%d.png",energy,harmonic));
  c1->Print(Form("FigsCheckThree/special_fvtxsL_B_energy%d_harm%d.pdf",energy,harmonic));

  hvn_fvtxs_B->Scale(1.0/reso_fvtx);
  hvn_fvtxs0_B->Scale(1.0/reso_fvtx0);
  hvn_fvtxs1_B->Scale(1.0/reso_fvtx1);
  hvn_fvtxs2_B->Scale(1.0/reso_fvtx2);
  hvn_fvtxs3_B->Scale(1.0/reso_fvtx3);

  if ( harmonic == 2 ) {MIN = -0.02; MAX = 0.2;}
  if ( harmonic == 3 ) {MIN = -0.02; MAX = 0.08;}
  delete h2dummy;
  h2dummy = new TH2D("h2dummy","",1,0.0,3.0,1,MIN,MAX);
  h2dummy->Draw();
  h2dummy->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h2dummy->GetYaxis()->SetTitle(Form("v_{%d}",harmonic));
  h2dummy->GetYaxis()->SetTitleOffset(1.3);
  hvn_fvtxs_B->Draw("same");
  hvn_fvtxs0_B->Draw("same");
  hvn_fvtxs1_B->Draw("same");
  hvn_fvtxs2_B->Draw("same");
  hvn_fvtxs3_B->Draw("same");
  line.Draw();
  leg->Draw();

  c1->Print(Form("FigsCheckThree/special_eprc_fvtxsL_B_energy%d_harm%d.png",energy,harmonic));
  c1->Print(Form("FigsCheckThree/special_eprc_fvtxsL_B_energy%d_harm%d.pdf",energy,harmonic));

  hvn_fvtxs_B->Rebin(5);
  hvn_fvtxs0_B->Rebin(5);
  hvn_fvtxs1_B->Rebin(5);
  hvn_fvtxs2_B->Rebin(5);
  hvn_fvtxs3_B->Rebin(5);

  h2dummy->Draw();
  hvn_fvtxs_B->Draw("same");
  hvn_fvtxs0_B->Draw("same");
  hvn_fvtxs1_B->Draw("same");
  hvn_fvtxs2_B->Draw("same");
  hvn_fvtxs3_B->Draw("same");
  line.Draw();
  leg->Draw();

  c1->Print(Form("FigsCheckThree/special_eprc_fvtxsL_B_energy%d_harm%d_REBIN.png",energy,harmonic));
  c1->Print(Form("FigsCheckThree/special_eprc_fvtxsL_B_energy%d_harm%d_REBIN.pdf",energy,harmonic));

  hvn_bbcs_B->Scale(1.0/reso_bbc);
  hvn_bbcs_B->Draw("same");
  leg->AddEntry(hvn_bbcs_B,"BBCS","elp");
  leg->Draw();

  c1->Print(Form("FigsCheckThree/special_eprcab_fvtxsL_B_energy%d_harm%d.png",energy,harmonic));
  c1->Print(Form("FigsCheckThree/special_eprcab_fvtxsL_B_energy%d_harm%d.pdf",energy,harmonic));



  // TH1D* th1d_hvn_fvtxs_B = hvn_fvtxs_B->ProjectionX();
  // for ( int i = 0; i < 15; ++i )
  //   {
  //     double value0 = hvn_fvtxs0_B->GetBinContent(i+1);
  //     double error0 = hvn_fvtxs0_B->GetBinError(i+1);
  //     double sigma0 = error0*error0;
  //     double value1 = hvn_fvtxs1_B->GetBinContent(i+1);
  //     double error1 = hvn_fvtxs1_B->GetBinError(i+1);
  //     double sigma1 = error1*error1;
  //     double value2 = hvn_fvtxs2_B->GetBinContent(i+1);
  //     double error2 = hvn_fvtxs2_B->GetBinError(i+1);
  //     double sigma2 = error2*error2;
  //     double value3 = hvn_fvtxs3_B->GetBinContent(i+1);
  //     double error3 = hvn_fvtxs3_B->GetBinError(i+1);
  //     double sigma3 = error3*error3;
  //     // ---
  //     sigma0 = hvn_fvtxs0_B->GetBinEntries(i+1);
  //     sigma1 = hvn_fvtxs1_B->GetBinEntries(i+1);
  //     sigma2 = hvn_fvtxs2_B->GetBinEntries(i+1);
  //     sigma3 = hvn_fvtxs3_B->GetBinEntries(i+1);
  //     // ---
  //     double numer0 = value0*sigma0;
  //     double numer1 = value1*sigma0;
  //     double numer2 = value2*sigma0;
  //     double numer3 = value3*sigma0;
  //     // ---
  //     double numerator = numer0 + numer1 + numer2 + numer3;
  //     double denominator = sigma0 + sigma1+ sigma2 + sigma3;
  //     numerator = numer0 + numer1 + numer2;
  //     denominator = sigma0 + sigma1+ sigma2;
  //     double finalvalue = -9; if ( denominator > 0 ) finalvalue = numerator/denominator;
  //     cout << "final value is " << finalvalue << endl;
  //     // ---
  //     th1d_hvn_fvtxs_B->SetBinContent(i+1,finalvalue);
  //   }

  // h2dummy->Draw();
  // hvn_bbcs_B->Draw("same ex0p");
  // hvn_fvtxs_B->Draw("same ex0p");
  // th1d_hvn_fvtxs_B->SetMarkerStyle(kOpenCross);
  // th1d_hvn_fvtxs_B->SetMarkerColor(kBlack);
  // th1d_hvn_fvtxs_B->SetLineColor(kBlack);
  // th1d_hvn_fvtxs_B->Draw("same ex0p");
  // line.Draw();

  // TLegend* leg2 = new TLegend(0.18,0.68,0.38,0.88);
  // leg2->AddEntry(hvn_bbcs_B,"BBCS","ep");
  // leg2->AddEntry(hvn_fvtxs_B,"FVTXS all layers","ep");
  // leg2->AddEntry(th1d_hvn_fvtxs_B,"FVTXS layers 0 1 2 ESTIMATE ONLY","ep");
  // leg2->SetTextSize(0.05);
  // leg2->Draw();

  // c1->Print(Form("FigsCheckThree/special_calcave_fvtxsL_B_energy%d_harm%d.png",energy,harmonic));
  // c1->Print(Form("FigsCheckThree/special_calcave_fvtxsL_B_energy%d_harm%d.pdf",energy,harmonic));

  delete c1;

}



void doenergy3(int energy, int harmonic)
{

  gStyle->SetOptTitle(1);

  TCanvas* c1 = new TCanvas("c1","");

  TFile* file = TFile::Open(Form("input/combined_%d.root",energy));

  // ---

  TProfile* tp1f_bbc_cnt = (TProfile*)file->Get(Form("tp1f_reso%d_BBC_CNT",harmonic));
  TProfile* tp1f_bbc_fvtx = (TProfile*)file->Get(Form("tp1f_reso%d_BBC_FVTX",harmonic));
  TProfile* tp1f_cnt_fvtx = (TProfile*)file->Get(Form("tp1f_reso%d_CNT_FVTX",harmonic));
  TProfile* tp1f_bbc_fvtx012 = (TProfile*)file->Get(Form("tp1f_reso%d_BBC_FVTX012",harmonic));
  TProfile* tp1f_cnt_fvtx012 = (TProfile*)file->Get(Form("tp1f_reso%d_CNT_FVTX012",harmonic));
  TProfile* tp1f_bbc_fvtx013 = (TProfile*)file->Get(Form("tp1f_reso%d_BBC_FVTX013",harmonic));
  TProfile* tp1f_cnt_fvtx013 = (TProfile*)file->Get(Form("tp1f_reso%d_CNT_FVTX013",harmonic));
  TProfile* tp1f_bbc_fvtx023 = (TProfile*)file->Get(Form("tp1f_reso%d_BBC_FVTX023",harmonic));
  TProfile* tp1f_cnt_fvtx023 = (TProfile*)file->Get(Form("tp1f_reso%d_CNT_FVTX023",harmonic));
  TProfile* tp1f_bbc_fvtx123 = (TProfile*)file->Get(Form("tp1f_reso%d_BBC_FVTX123",harmonic));
  TProfile* tp1f_cnt_fvtx123 = (TProfile*)file->Get(Form("tp1f_reso%d_CNT_FVTX123",harmonic));

  float float_bbc_cnt = tp1f_bbc_cnt->GetBinContent(1);
  float float_bbc_fvtx = tp1f_bbc_fvtx->GetBinContent(1);
  float float_cnt_fvtx = tp1f_cnt_fvtx->GetBinContent(1);
  float float_bbc_fvtx012 = tp1f_bbc_fvtx012->GetBinContent(1);
  float float_cnt_fvtx012 = tp1f_cnt_fvtx012->GetBinContent(1);
  float float_bbc_fvtx013 = tp1f_bbc_fvtx013->GetBinContent(1);
  float float_cnt_fvtx013 = tp1f_cnt_fvtx013->GetBinContent(1);
  float float_bbc_fvtx023 = tp1f_bbc_fvtx023->GetBinContent(1);
  float float_cnt_fvtx023 = tp1f_cnt_fvtx023->GetBinContent(1);
  float float_bbc_fvtx123 = tp1f_bbc_fvtx123->GetBinContent(1);
  float float_cnt_fvtx123 = tp1f_cnt_fvtx123->GetBinContent(1);

  cout << "bbc-cnt " << float_bbc_cnt << endl;
  cout << "bbc-fvtxs " << float_bbc_fvtx << endl;
  cout << "cnt-fvtxs " << float_cnt_fvtx << endl;


  float reso_bbc = sqrt((float_bbc_cnt*float_bbc_fvtx)/float_cnt_fvtx); // BCBS/CS
  float reso_fvtx = sqrt((float_cnt_fvtx*float_bbc_fvtx)/float_bbc_cnt); // CSBS/BC
  float reso_fvtx012 = sqrt((float_cnt_fvtx012*float_bbc_fvtx012)/float_bbc_cnt); // CSBS/BC
  float reso_fvtx013 = sqrt((float_cnt_fvtx013*float_bbc_fvtx013)/float_bbc_cnt); // CSBS/BC
  float reso_fvtx023 = sqrt((float_cnt_fvtx023*float_bbc_fvtx023)/float_bbc_cnt); // CSBS/BC
  float reso_fvtx123 = sqrt((float_cnt_fvtx123*float_bbc_fvtx123)/float_bbc_cnt); // CSBS/BC

  cout << "bbc resolution is " << reso_bbc << endl;
  cout << "fvtx resolution is " << reso_fvtx << endl;
  cout << "fvtx012 resolution is " << reso_fvtx012 << endl;
  cout << "fvtx013 resolution is " << reso_fvtx013 << endl;
  cout << "fvtx023 resolution is " << reso_fvtx023 << endl;
  cout << "fvtx123 resolution is " << reso_fvtx123 << endl;

  // ---

  //return;

  TProfile* hvn_bbcs_B = (TProfile*)file->Get(Form("bbcs_v%d_both_docalib",harmonic));
  hvn_bbcs_B->SetMarkerStyle(kOpenSquare);
  hvn_bbcs_B->SetMarkerColor(kBlack);
  hvn_bbcs_B->SetLineColor(kBlack);

  TProfile* hvn_fvtxs_B = (TProfile*)file->Get(Form("fvtxs_v%d_both_docalib",harmonic));
  TProfile* hvn_fvtxs012_B = (TProfile*)file->Get(Form("fvtxs012_v%d_both_docalib",harmonic));
  TProfile* hvn_fvtxs013_B = (TProfile*)file->Get(Form("fvtxs013_v%d_both_docalib",harmonic));
  TProfile* hvn_fvtxs023_B = (TProfile*)file->Get(Form("fvtxs023_v%d_both_docalib",harmonic));
  TProfile* hvn_fvtxs123_B = (TProfile*)file->Get(Form("fvtxs123_v%d_both_docalib",harmonic));

  hvn_fvtxs_B->SetMarkerStyle(kOpenCircle);
  hvn_fvtxs_B->SetMarkerColor(kBlack);
  hvn_fvtxs_B->SetLineColor(kBlack);
  hvn_fvtxs012_B->SetLineColor(kBlue);
  hvn_fvtxs013_B->SetLineColor(kRed);
  hvn_fvtxs023_B->SetLineColor(kGreen+2);
  hvn_fvtxs123_B->SetLineColor(kMagenta+2);

  double MIN = 0;
  double MAX = 1;
  if ( harmonic == 2 ) {MIN = -0.01; MAX = 0.05;}
  if ( harmonic == 3 ) {MIN = -0.002; MAX = 0.005;}
  TH2D* h2dummy = new TH2D("h2dummy","",1,0.0,3.0,1,MIN,MAX);
  h2dummy->Draw();
  TLine line(0,0,3,0);
  line.SetLineStyle(2);
  line.SetLineWidth(2);
  line.Draw();
  h2dummy->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h2dummy->GetYaxis()->SetTitle(Form("v_{%d} not corrected for EP resolution",harmonic));
  h2dummy->GetYaxis()->SetTitleOffset(1.3);

  hvn_fvtxs_B->Draw("same");
  hvn_fvtxs012_B->Draw("same");
  hvn_fvtxs013_B->Draw("same");
  hvn_fvtxs023_B->Draw("same");
  hvn_fvtxs123_B->Draw("same");

  TLegend* leg = new TLegend(0.18,0.68,0.38,0.88);
  leg->SetHeader(Form("%d GeV, FVTXS",energy));
  leg->AddEntry(hvn_fvtxs_B,"all layers","elp");
  leg->AddEntry(hvn_fvtxs012_B,"layer 012","el");
  leg->AddEntry(hvn_fvtxs013_B,"layer 013","el");
  leg->AddEntry(hvn_fvtxs023_B,"layer 023","el");
  leg->AddEntry(hvn_fvtxs123_B,"layer 123","el");
  leg->SetTextSize(0.05);
  leg->Draw();

  c1->Print(Form("FigsCheckThree/superspecial_fvtxsL_B_energy%d_harm%d.png",energy,harmonic));
  c1->Print(Form("FigsCheckThree/superspecial_fvtxsL_B_energy%d_harm%d.pdf",energy,harmonic));

  hvn_fvtxs_B->Scale(1.0/reso_fvtx);
  hvn_fvtxs012_B->Scale(1.0/reso_fvtx012);
  hvn_fvtxs013_B->Scale(1.0/reso_fvtx013);
  hvn_fvtxs023_B->Scale(1.0/reso_fvtx023);
  hvn_fvtxs123_B->Scale(1.0/reso_fvtx123);

  if ( harmonic == 2 ) {MIN = -0.02; MAX = 0.2;}
  if ( harmonic == 3 ) {MIN = -0.02; MAX = 0.08;}
  delete h2dummy;
  h2dummy = new TH2D("h2dummy","",1,0.0,3.0,1,MIN,MAX);
  h2dummy->Draw();
  h2dummy->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h2dummy->GetYaxis()->SetTitle(Form("v_{%d}",harmonic));
  h2dummy->GetYaxis()->SetTitleOffset(1.3);
  hvn_fvtxs_B->Draw("same");
  hvn_fvtxs012_B->Draw("same");
  hvn_fvtxs013_B->Draw("same");
  hvn_fvtxs023_B->Draw("same");
  hvn_fvtxs123_B->Draw("same");
  line.Draw();
  leg->Draw();

  c1->Print(Form("FigsCheckThree/superspecial_eprc_fvtxsL_B_energy%d_harm%d.png",energy,harmonic));
  c1->Print(Form("FigsCheckThree/superspecial_eprc_fvtxsL_B_energy%d_harm%d.pdf",energy,harmonic));

  hvn_bbcs_B->Scale(1.0/reso_bbc);
  hvn_bbcs_B->Draw("same");
  leg->AddEntry(hvn_bbcs_B,"BBCS","elp");
  leg->Draw();

  c1->Print(Form("FigsCheckThree/superspecial_eprcab_fvtxsL_B_energy%d_harm%d.png",energy,harmonic));
  c1->Print(Form("FigsCheckThree/superspecial_eprcab_fvtxsL_B_energy%d_harm%d.pdf",energy,harmonic));



  delete c1;

}
