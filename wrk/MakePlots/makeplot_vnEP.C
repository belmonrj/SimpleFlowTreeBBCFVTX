void doenergy(int, int);
void diagnostic(int, int);

TFile* outfile;

void makeplot_vnEP()
{

  outfile = TFile::Open("histograms_vnEP.root","recreate");

  for ( int i = 2; i < 4; ++i )
    {
      doenergy(200,i);
      doenergy(62,i);
      doenergy(39,i);
      doenergy(20,i);
      // ---
      diagnostic(200,i);
      diagnostic(62,i);
      diagnostic(39,i);
      diagnostic(20,i);
    }

  outfile->Write();
  outfile->Close();

}

void doenergy(int energy, int harmonic)
{

  gStyle->SetOptTitle(1);

  TCanvas* c1 = new TCanvas("c1","");

  TFile* file = TFile::Open(Form("input/combined_%d.root",energy));

  // ---

  TProfile* tp1f_bbc_fvtxN = (TProfile*)file->Get(Form("tp1f_reso%d_BBC_FVTXN",harmonic));
  TProfile* tp1f_bbc_fvtx = (TProfile*)file->Get(Form("tp1f_reso%d_BBC_FVTX",harmonic));
  TProfile* tp1f_fvtxN_fvtx = (TProfile*)file->Get(Form("tp1f_reso%d_FVTXS_FVTXN",harmonic));

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

  TProfile* tp1f_bbc_cnt = (TProfile*)file->Get(Form("tp1f_reso%d_BBC_CNT",harmonic));
  TProfile* tp1f_cnt_fvtx = (TProfile*)file->Get(Form("tp1f_reso%d_CNT_FVTX",harmonic));
  TProfile* tp1f_cnt_fvtxN = (TProfile*)file->Get(Form("tp1f_reso%d_CNT_FVTXN",harmonic));

  float float_bbc_cnt = tp1f_bbc_cnt->GetBinContent(1);
  float float_cnt_fvtx = tp1f_cnt_fvtx->GetBinContent(1);
  float float_cnt_fvtxN = tp1f_cnt_fvtxN->GetBinContent(1);

  cout << "bbc-cnt " << float_bbc_cnt << endl;
  cout << "cnt-fvtxs " << float_cnt_fvtx << endl;
  cout << "cbt-fvtxn " << float_cnt_fvtxN << endl;

  if ( float_bbc_cnt < 0 ) { cout << "YOU'RE GONNA DIE (bbcs-cnt)" << endl; }
  if ( float_cnt_fvtx < 0 ) { cout << "YOU'RE GONNA DIE (cnt-fvtxs)" << endl; }
  if ( float_cnt_fvtxN < 0 ) { cout << "YOU'RE GONNA DIE (cnt-fvtxn)" << endl; }

  float reso_bbc = sqrt((float_bbc_cnt*float_bbc_fvtx)/float_cnt_fvtx); // BCBS/CS
  float reso_fvtx = sqrt((float_cnt_fvtx*float_bbc_fvtx)/float_bbc_cnt); // CSBS/BC
  float reso_fvtxN = sqrt((float_cnt_fvtxN*float_bbc_fvtxN)/float_bbc_cnt); // CNBN/BC

  cout << "bbc resolution is " << reso_bbc << endl;
  cout << "fvtx resolution is " << reso_fvtx << endl;
  cout << "fvtxN resolution is " << reso_fvtxN << endl;

  // ---

  float reso_fvtxN_xb = sqrt ( ( float_fvtxN_fvtx * float_cnt_fvtxN ) / float_cnt_fvtx ) ; // NSNC/SC

  cout << "fvtxN resolution is " << reso_fvtxN_xb << endl;

  // ---

  if ( energy == 20 )
    {
      // --- i like chocolate and fudge and other tasty deserts
      cout << "Now making empirical adjustment for 20 GeV data" << endl;
      reso_bbc = 0.015;
      reso_fvtx = 0.04;
    }

  // ---

  TProfile* hvn_fvtxs = (TProfile*)file->Get(Form("fvtxs_v%d_both_docalib",harmonic));
  hvn_fvtxs->Scale(1.0/reso_fvtx);
  ofstream fout((const char*)Form("DataTextFiles/run16dau%d_v%dfvtxs.dat",energy,harmonic));
  for ( int i = 0; i < hvn_fvtxs->GetNbinsX(); ++i ) fout << hvn_fvtxs->GetBinCenter(i+1) << " " << hvn_fvtxs->GetBinContent(i+1) << " " << hvn_fvtxs->GetBinError(i+1) << " " << endl;
  fout.close();
  hvn_fvtxs->Draw();
  // the 62 GeV is actually 62.4 and the 20 GeV is actually 19.6, so need to modify
  hvn_fvtxs->SetTitle(Form("d+Au collisions at #sqrt{s_{NN}} = %d GeV",energy));
  if ( energy == 62 ) hvn_fvtxs->SetTitle(Form("d+Au collisions at #sqrt{s_{NN}} = 62.4 GeV"));
  if ( energy == 20 ) hvn_fvtxs->SetTitle(Form("d+Au collisions at #sqrt{s_{NN}} = 19.6 GeV"));
  hvn_fvtxs->SetMaximum(0.17);
  if ( energy == 39 ) hvn_fvtxs->SetMaximum(0.27);
  hvn_fvtxs->SetMinimum(0.0);
  TLine line(0,0,3,0);
  line.SetLineStyle(2);
  line.SetLineWidth(2);
  if ( energy == 20 )
    {
      hvn_fvtxs->SetMaximum(0.5);
      hvn_fvtxs->SetMinimum(-0.2);
      line.Draw();
    }
  if ( harmonic == 3 )
    {
      hvn_fvtxs->SetMaximum(0.08);
      hvn_fvtxs->SetMinimum(-0.02);
      line.Draw();
    }
  hvn_fvtxs->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hvn_fvtxs->GetYaxis()->SetTitle(Form("v_{%d}{EP}",harmonic));
  hvn_fvtxs->GetYaxis()->SetTitleOffset(1.25);

  ifstream finpub("ppg161.dat");
  float pt[13], pubvn[13], epubvn[13], esyspubvn[13];
  for ( int i = 0; i < 13; ++i )
    {
      finpub>>pt[i]>>pubvn[i]>>epubvn[i]>>esyspubvn[i];
    }

  TGraphErrors* tge_pub = new TGraphErrors(13,pt,pubvn,0,epubvn);
  tge_pub->SetMarkerStyle(kFullCircle);
  if ( harmonic == 2 ) tge_pub->Draw("p");

  c1->Print(Form("FigsHarmonicCoefficient/run16dau%d_v%d_fvtxs.pdf",energy,harmonic));
  c1->Print(Form("FigsHarmonicCoefficient/run16dau%d_v%d_fvtxs.png",energy,harmonic));

  TProfile* hvn_bbcs = (TProfile*)file->Get(Form("bbcs_v%d_both_docalib",harmonic));
  hvn_bbcs->SetLineColor(kRed);
  hvn_bbcs->Scale(1.0/reso_bbc);
  hvn_bbcs->Draw("same");

  if ( energy == 20 )
    {
      hvn_fvtxs->SetMaximum(1.0);
      hvn_fvtxs->SetMinimum(-0.5);
      line.Draw();
    }

  TLegend* leg = new TLegend(0.18,0.68,0.38,0.88);
  leg->AddEntry(hvn_fvtxs,"Run16 FVTXS","el");
  leg->AddEntry(hvn_bbcs,"Run16 BBCS","el");
  if ( energy == 200 && harmonic == 2 ) leg->AddEntry(tge_pub,"Run8 (ppg161)","p");
  else if ( harmonic == 2 ) leg->AddEntry(tge_pub,"Run8 (200 GeV)","p");
  leg->SetTextSize(0.05);
  leg->Draw();

  TLatex* latex = new TLatex(0.2,0.2,"EP resolutions calculated with CNT");
  latex->SetNDC();
  latex->SetTextSize(0.05);
  //latex->Draw();

  c1->Print(Form("FigsHarmonicCoefficient/run16dau%d_v%d_fvtxsbbcs.pdf",energy,harmonic));
  c1->Print(Form("FigsHarmonicCoefficient/run16dau%d_v%d_fvtxsbbcs.png",energy,harmonic));

  TProfile* xhvn_bbcs = (TProfile*)hvn_bbcs->Clone();
  TProfile* xhvn_fvtxs = (TProfile*)hvn_fvtxs->Clone();
  xhvn_fvtxs->Add(xhvn_bbcs);
  xhvn_fvtxs->SetLineColor(kBlack);
  xhvn_fvtxs->SetMarkerColor(kBlack);
  xhvn_fvtxs->SetMarkerStyle(kOpenCircle);
  xhvn_fvtxs->Draw("same");

  double x[15], y[15], exh[15], eyh[15], exl[15], eyl[15];
  for ( int i = 0; i < 15; ++i )
    {
      x[i] = xhvn_fvtxs->GetBinCenter(i+1);
      y[i] = xhvn_fvtxs->GetBinContent(i+1);
      exl[i] = 0;
      exh[i] = 0;
      double a = hvn_fvtxs->GetBinContent(i+1) - y[i];
      double b = hvn_bbcs->GetBinContent(i+1) - y[i];
      if ( a > b )
        {
          eyl[i] = fabs(b);
          eyh[i] = fabs(a);
        }
      else
        {
          eyl[i] = fabs(a);
          eyh[i] = fabs(b);
        }
    }
  TGraphAsymmErrors* tgae = new TGraphAsymmErrors(15,x,y,exl,exh,eyl,eyh);
  tgae->SetMarkerStyle(1);
  tgae->SetMarkerColor(kGray);
  tgae->SetLineColor(kGray);
  tgae->SetLineWidth(15);
  tgae->Draw("pz");
  xhvn_fvtxs->Draw("same");
  if ( harmonic == 3 ) line.Draw();

  c1->Print(Form("FigsHarmonicCoefficient/run16dau%d_v%d_fvtxsbbcsA.pdf",energy,harmonic));
  c1->Print(Form("FigsHarmonicCoefficient/run16dau%d_v%d_fvtxsbbcsA.png",energy,harmonic));

  xhvn_fvtxs->Draw();
  if ( harmonic == 3 ) line.Draw();

  c1->Print(Form("FigsHarmonicCoefficient/run16dau%d_v%d_fvtxsbbcsAO.pdf",energy,harmonic));
  c1->Print(Form("FigsHarmonicCoefficient/run16dau%d_v%d_fvtxsbbcsAO.png",energy,harmonic));



  TProfile* hvn_fvtxn = (TProfile*)file->Get(Form("fvtxn_v%d_both_docalib",harmonic));
  hvn_fvtxn->SetLineColor(kGreen+2);
  hvn_fvtxn->Scale(1.0/reso_fvtxN);
  hvn_fvtxn->Draw("same");

  delete leg;

  TLegend* leg3 = new TLegend(0.18,0.68,0.38,0.88);
  leg3->AddEntry(hvn_fvtxs,"Run16 FVTXS","el");
  leg3->AddEntry(hvn_fvtxn,"Run16 FVTXN","el");
  leg3->AddEntry(hvn_bbcs,"Run16 BBCS","el");
  if ( energy == 200 && harmonic == 2 ) leg3->AddEntry(tge_pub,"Run8 (ppg161)","p");
  else if ( harmonic == 2 ) leg3->AddEntry(tge_pub,"Run8 (200 GeV)","p");
  leg3->SetTextSize(0.05);
  leg3->Draw();

  //latex->Draw();

  c1->Print(Form("FigsHarmonicCoefficient/run16dau%d_v%d_fvtxsnbbcs.pdf",energy,harmonic));
  c1->Print(Form("FigsHarmonicCoefficient/run16dau%d_v%d_fvtxsnbbcs.png",energy,harmonic));

  // hvn_fvtxn->Scale(reso_fvtxN/reso_fvtxN_fn);
  // hvn_fvtxs->Scale(reso_fvtx/reso_fvtx_fn);
  // hvn_bbcs->Scale(reso_bbc/reso_bbc_fn);

  // delete latex;
  // latex = new TLatex(0.2,0.2,"EP resolutions calculated without CNT");
  // latex->SetNDC();
  // latex->SetTextSize(0.05);
  // latex->Draw();

  // //  latex->DrawLatex(0.2,0.2,"EP resolutions calculated without CNT");

  // c1->Print(Form("FigsHarmonicCoefficient/run16dau%d_v%d_fvtxsnbbcs_fn.pdf",energy,harmonic));
  // c1->Print(Form("FigsHarmonicCoefficient/run16dau%d_v%d_fvtxsnbbcs_fn.png",energy,harmonic));

  // c1->Clear();
  // hvn_fvtxs->Draw();
  // hvn_bbcs->Draw("same");
  // tge_pub->Draw("p");
  // leg = new TLegend(0.18,0.68,0.38,0.88);
  // leg->AddEntry(hvn_fvtxs,"Run16 FVTXS","el");
  // leg->AddEntry(hvn_bbcs,"Run16 BBCS","el");
  // if ( energy == 200 && harmonic == 2 ) leg->AddEntry(tge_pub,"Run8 (ppg161)","p");
  // else if ( harmonic == 2 ) leg->AddEntry(tge_pub,"Run8 (200 GeV)","p");
  // leg->SetTextSize(0.05);
  // leg->Draw();

  // latex->Draw();

  // c1->Print(Form("FigsHarmonicCoefficient/run16dau%d_v%d_fvtxsbbcs_fn.pdf",energy,harmonic));
  // c1->Print(Form("FigsHarmonicCoefficient/run16dau%d_v%d_fvtxsbbcs_fn.png",energy,harmonic));

  // hvn_fvtxn->Scale(2*reso_fvtxN_fn/(reso_fvtxN_fn+reso_fvtxN));
  // hvn_fvtxs->Scale(2*reso_fvtx_fn/(reso_fvtx_fn+reso_fvtx));
  // hvn_bbcs->Scale(2*reso_bbc_fn/(reso_bbc_fn+reso_bbc));

  // delete latex;
  // latex = new TLatex(0.2,0.2,"Average EP resolutions with and without CNT");
  // latex->SetNDC();
  // latex->SetTextSize(0.05);
  // latex->Draw();

  // //  latex->DrawLatex(0.2,0.2,"Average of EP resolutions with and without CNT");

  // c1->Print(Form("FigsHarmonicCoefficient/run16dau%d_v%d_fvtxsbbcs_ave.pdf",energy,harmonic));
  // c1->Print(Form("FigsHarmonicCoefficient/run16dau%d_v%d_fvtxsbbcs_ave.png",energy,harmonic));

  // delete leg;
  // hvn_fvtxn->Draw("same");
  // leg3->Draw();

  //latex->Draw();

  c1->Print(Form("FigsHarmonicCoefficient/run16dau%d_v%d_fvtxsnbbcs_ave.pdf",energy,harmonic));
  c1->Print(Form("FigsHarmonicCoefficient/run16dau%d_v%d_fvtxsnbbcs_ave.png",energy,harmonic));

  c1->Clear();



  // --- now eta plots

  TProfile* hvneta_bbcs = (TProfile*)file->Get(Form("bbcs_v%deta_both_docalib",harmonic));
  hvneta_bbcs->SetLineColor(kBlack);
  hvneta_bbcs->Scale(1.0/reso_bbc);
  hvneta_bbcs->Draw();
  hvneta_bbcs->SetTitle(Form("d+Au collisions at #sqrt{s_{NN}} = %d GeV",energy));
  if ( energy == 62 ) hvneta_bbcs->SetTitle(Form("d+Au collisions at #sqrt{s_{NN}} = 62.4 GeV"));
  if ( energy == 20 ) hvneta_bbcs->SetTitle(Form("d+Au collisions at #sqrt{s_{NN}} = 19.6 GeV"));
  hvneta_bbcs->SetMaximum(0.1);
  hvneta_bbcs->SetMinimum(0.0);
  if ( energy <= 39 || harmonic == 3 ) hvneta_bbcs->SetMinimum(-0.02);
  if ( energy == 20 )
    {
      hvneta_bbcs->SetMaximum(0.3);
      hvneta_bbcs->SetMinimum(-0.2);
      TLine line(-3.1,0,3.1,0);
      line.SetLineStyle(2);
      line.SetLineWidth(2);
      line.Draw();
    }
  hvneta_bbcs->GetXaxis()->SetTitle("#eta");
  hvneta_bbcs->GetYaxis()->SetTitle(Form("v_{%d}{EP}",harmonic));
  hvneta_bbcs->GetYaxis()->SetTitleOffset(1.25);
  TLegend* leta = new TLegend(0.68,0.68,0.88,0.88);
  leta->AddEntry(hvneta_bbcs,"BBCS","el");
  leta->SetTextSize(0.05);
  leta->Draw();

  c1->Print(Form("FigsHarmonicCoefficient/run16dau%d_v%deta_bbcs.pdf",energy,harmonic));
  c1->Print(Form("FigsHarmonicCoefficient/run16dau%d_v%deta_bbcs.png",energy,harmonic));

  TProfile* hvneta_fvtxs = (TProfile*)file->Get(Form("fvtxs_v%deta_both_docalib",harmonic));
  hvneta_fvtxs->SetLineColor(kBlue);
  hvneta_fvtxs->Scale(1.0/reso_fvtx);
  hvneta_fvtxs->Draw("same");
  leta->AddEntry(hvneta_fvtxs,"FVTXS","el");
  leta->Draw();

  c1->Print(Form("FigsHarmonicCoefficient/run16dau%d_v%deta_fvtxsbbcs.pdf",energy,harmonic));
  c1->Print(Form("FigsHarmonicCoefficient/run16dau%d_v%deta_fvtxsbbcs.png",energy,harmonic));

  // --- clean
  hvneta_bbcs->SetBinContent(1,-9);
  hvneta_bbcs->SetBinContent(12,-9);
  hvneta_bbcs->SetBinContent(21,-9);
  hvneta_bbcs->SetBinContent(32,-9);
  hvneta_fvtxs->SetBinContent(21,-9);
  hvneta_fvtxs->SetBinContent(32,-9);

  c1->Print(Form("FigsHarmonicCoefficient/run16dau%d_v%deta_cleaned_fvtxsbbcs.pdf",energy,harmonic));
  c1->Print(Form("FigsHarmonicCoefficient/run16dau%d_v%deta_cleaned_fvtxsbbcs.png",energy,harmonic));

  TH1D* hvneta_bbcs_m = (TH1D*)hvneta_bbcs->Clone();
  TH1D* hvneta_fvtxs_m = (TH1D*)hvneta_fvtxs->Clone();
  hvneta_fvtxs_m->Divide(hvneta_bbcs_m);
  hvneta_fvtxs_m->Draw();
  hvneta_fvtxs_m->SetMinimum(0.0);
  hvneta_fvtxs_m->SetMaximum(2.0);
  c1->Print(Form("FigsHarmonicCoefficient/run16dau%d_v%deta_ratio_fvtxsbbcs.pdf",energy,harmonic));
  c1->Print(Form("FigsHarmonicCoefficient/run16dau%d_v%deta_ratio_fvtxsbbcs.png",energy,harmonic));

  // TFile* fcorr = TFile::Open("ampt_vneta_correction.root");
  // //TProfile* tp1f_corr = (TProfile*)fcorr->Get("hv2eta_corr_sys");
  // TProfile* tp1f_corr = (TProfile*)fcorr->Get("v2eta_correction_dAu200");
  // TH1D* th1d_corr = tp1f_corr->ProjectionX();
  // hvneta_bbcs->Divide(th1d_corr);
  // hvneta_fvtxs->Divide(th1d_corr);
  // hvneta_bbcs->Draw();
  // hvneta_fvtxs->Draw("same");
  // leta->Draw();

  // c1->Print(Form("FigsHarmonicCoefficient/run16dau%d_v%deta_corrected_fvtxsbbcs.pdf",energy,harmonic));
  // c1->Print(Form("FigsHarmonicCoefficient/run16dau%d_v%deta_corrected_fvtxsbbcs.png",energy,harmonic));

  delete c1;


  if ( harmonic == 3 ) hvn_fvtxs = xhvn_fvtxs;

  hvneta_fvtxs->SetTitle(Form("d+Au collisions at #sqrt{s_{NN}} = %d GeV",energy));
  if ( energy == 62 ) hvneta_fvtxs->SetTitle(Form("d+Au collisions at #sqrt{s_{NN}} = 62.4 GeV"));
  if ( energy == 20 ) hvneta_fvtxs->SetTitle(Form("d+Au collisions at #sqrt{s_{NN}} = 19.6 GeV"));

  outfile->cd();
  hvn_fvtxs->SetName(Form("tprofile_v%d_pT_eventplane_fvtxs_%d",harmonic,energy));
  hvneta_bbcs->SetName(Form("tprofile_v%d_eta_eventplane_bbcs_%d",harmonic,energy));
  hvneta_fvtxs->SetName(Form("tprofile_v%d_eta_eventplane_fvtxs_%d",harmonic,energy));
  hvn_fvtxs->Write();
  hvneta_bbcs->Write();
  hvneta_fvtxs->Write();

}


void diagnostic(int energy, int harmonic)
{

  gStyle->SetOptTitle(1);

  TCanvas* c1 = new TCanvas("c1","");

  TFile* file = TFile::Open(Form("input/combined_%d.root",energy));



  float MIN = -1;
  float MAX = 1;



  // ---------
  // --- FVTXS
  // ---------

  TProfile* hvn_fvtxs_B = (TProfile*)file->Get(Form("fvtxs_v%d_both_docalib",harmonic));
  TProfile* hvn_fvtxs_E = (TProfile*)file->Get(Form("fvtxs_v%d_east_docalib",harmonic));
  TProfile* hvn_fvtxs_W = (TProfile*)file->Get(Form("fvtxs_v%d_west_docalib",harmonic));

  hvn_fvtxs_B->SetLineColor(kBlack);
  hvn_fvtxs_E->SetLineColor(kRed);
  hvn_fvtxs_W->SetLineColor(kBlue);


  if ( energy == 200 ) {MIN = -0.01; MAX = 0.05;}
  if ( energy == 62 ) {MIN = -0.01; MAX = 0.05;}
  if ( energy == 39 ) {MIN = -0.03; MAX = 0.05;}
  if ( energy == 20 ) {MIN = -0.1; MAX = 0.2;}
  if ( harmonic == 3 ) {MIN = -0.003; MAX = 0.008;}
  TLine line(0,0,3,0);
  line.SetLineStyle(2);
  line.SetLineWidth(2);
  TH2D* h2dummy = new TH2D("h2dummy","",1,0.0,3.0,1,MIN,MAX);
  h2dummy->Draw();
  line.Draw();
  h2dummy->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h2dummy->GetYaxis()->SetTitle(Form("v_{%d} not corrected for EP resolution",harmonic));
  h2dummy->GetYaxis()->SetTitleOffset(1.3);

  hvn_fvtxs_B->Draw("same");
  hvn_fvtxs_E->Draw("same");
  hvn_fvtxs_W->Draw("same");

  TLegend* leg = new TLegend(0.18,0.68,0.38,0.88);
  leg->SetHeader(Form("%d GeV, FVTXS",energy));
  leg->AddEntry(hvn_fvtxs_B,"both","el");
  leg->AddEntry(hvn_fvtxs_E,"east","el");
  leg->AddEntry(hvn_fvtxs_W,"west","el");
  leg->SetTextSize(0.05);
  leg->Draw();

  c1->Print(Form("FigsHarmonicCoefficient/diagnostic_fvtxs_EBW_energy%d_harm%d.png",energy,harmonic));
  c1->Print(Form("FigsHarmonicCoefficient/diagnostic_fvtxs_EBW_energy%d_harm%d.pdf",energy,harmonic));


  if ( harmonic == 2 )
  {
  TProfile* hvn_fvtxs0_B = (TProfile*)file->Get(Form("fvtxs0_v%d_both_docalib",harmonic));
  TProfile* hvn_fvtxs1_B = (TProfile*)file->Get(Form("fvtxs1_v%d_both_docalib",harmonic));
  TProfile* hvn_fvtxs2_B = (TProfile*)file->Get(Form("fvtxs2_v%d_both_docalib",harmonic));
  TProfile* hvn_fvtxs3_B = (TProfile*)file->Get(Form("fvtxs3_v%d_both_docalib",harmonic));

  hvn_fvtxs0_B->SetLineColor(kBlue);
  hvn_fvtxs1_B->SetLineColor(kRed);
  hvn_fvtxs2_B->SetLineColor(kGreen+2);
  hvn_fvtxs3_B->SetLineColor(kMagenta+2);

  if ( energy == 200 ) {MIN = -0.01; MAX = 0.05;}
  if ( energy == 62 ) {MIN = -0.01; MAX = 0.05;}
  if ( energy == 39 ) {MIN = -0.03; MAX = 0.05;}
  if ( energy == 20 ) {MIN = -0.1; MAX = 0.2;}
  if ( harmonic == 3 ) {MIN = -0.003; MAX = 0.008;}
  delete h2dummy;
  h2dummy = new TH2D("h2dummy","",1,0.0,3.0,1,MIN,MAX);
  h2dummy->Draw();
  line.Draw();
  h2dummy->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h2dummy->GetYaxis()->SetTitle(Form("v_{%d} not corrected for EP resolution",harmonic));
  h2dummy->GetYaxis()->SetTitleOffset(1.3);

  hvn_fvtxs_B->Draw("same");
  hvn_fvtxs0_B->Draw("same");
  hvn_fvtxs1_B->Draw("same");
  hvn_fvtxs2_B->Draw("same");
  hvn_fvtxs3_B->Draw("same");

  delete leg;
  leg = new TLegend(0.18,0.68,0.38,0.88);
  leg->SetHeader(Form("%d GeV, FVTXS",energy));
  leg->AddEntry(hvn_fvtxs_B,"all layers","el");
  leg->AddEntry(hvn_fvtxs0_B,"layer 0","el");
  leg->AddEntry(hvn_fvtxs1_B,"layer 1","el");
  leg->AddEntry(hvn_fvtxs2_B,"layer 2","el");
  leg->AddEntry(hvn_fvtxs3_B,"layer 3","el");
  leg->SetTextSize(0.05);
  leg->Draw();

  c1->Print(Form("FigsHarmonicCoefficient/diagnostic_fvtxsL_B_energy%d_harm%d.png",energy,harmonic));
  c1->Print(Form("FigsHarmonicCoefficient/diagnostic_fvtxsL_B_energy%d_harm%d.pdf",energy,harmonic));
  }


  // ---------
  // --- BBCS
  // ---------

  TProfile* hvn_bbcs_B = (TProfile*)file->Get(Form("bbcs_v%d_both_docalib",harmonic));
  TProfile* hvn_bbcs_E = (TProfile*)file->Get(Form("bbcs_v%d_east_docalib",harmonic));
  TProfile* hvn_bbcs_W = (TProfile*)file->Get(Form("bbcs_v%d_west_docalib",harmonic));

  hvn_bbcs_B->SetLineColor(kBlack);
  hvn_bbcs_E->SetLineColor(kRed);
  hvn_bbcs_W->SetLineColor(kBlue);

  if ( energy == 200 ) {MIN = -0.005; MAX = 0.03;}
  if ( energy == 62 ) {MIN = -0.005; MAX = 0.02;}
  if ( energy == 39 ) {MIN = -0.02; MAX = 0.02;}
  if ( energy == 20 ) {MIN = -0.05; MAX = 0.05;}
  if ( harmonic == 3 ) {MIN = -0.002; MAX = 0.008;}
  delete h2dummy;
  h2dummy = new TH2D("h2dummy","",1,0.0,3.0,1,MIN,MAX);
  h2dummy->Draw();
  line.Draw();
  h2dummy->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h2dummy->GetYaxis()->SetTitle(Form("v_{%d} not corrected for EP resolution",harmonic));
  h2dummy->GetYaxis()->SetTitleOffset(1.3);

  hvn_bbcs_B->Draw("same");
  hvn_bbcs_E->Draw("same");
  hvn_bbcs_W->Draw("same");

  delete leg;
  leg = new TLegend(0.18,0.68,0.38,0.88);
  leg->SetHeader(Form("%d GeV, BBCS",energy));
  leg->AddEntry(hvn_fvtxs_B,"both","el");
  leg->AddEntry(hvn_fvtxs_E,"east","el");
  leg->AddEntry(hvn_fvtxs_W,"west","el");
  leg->SetTextSize(0.05);
  leg->Draw();

  c1->Print(Form("FigsHarmonicCoefficient/diagnostic_bbcs_EBW_energy%d_harm%d.png",energy,harmonic));
  c1->Print(Form("FigsHarmonicCoefficient/diagnostic_bbcs_EBW_energy%d_harm%d.pdf",energy,harmonic));



  // ---------
  // --- FVTXN
  // ---------

  TProfile* hvn_fvtxn_B = (TProfile*)file->Get(Form("fvtxn_v%d_both_docalib",harmonic));
  TProfile* hvn_fvtxn_E = (TProfile*)file->Get(Form("fvtxn_v%d_east_docalib",harmonic));
  TProfile* hvn_fvtxn_W = (TProfile*)file->Get(Form("fvtxn_v%d_west_docalib",harmonic));

  hvn_fvtxn_B->SetLineColor(kBlack);
  hvn_fvtxn_E->SetLineColor(kRed);
  hvn_fvtxn_W->SetLineColor(kBlue);

  if ( energy == 200 ) {MIN = -0.01; MAX = 0.05;}
  if ( energy == 62 ) {MIN = -0.01; MAX = 0.05;}
  if ( energy == 39 ) {MIN = -0.03; MAX = 0.05;}
  if ( energy == 20 ) {MIN = -0.1; MAX = 0.2;}
  if ( harmonic == 3 ) {MIN = -0.003; MAX = 0.008;}
  delete h2dummy;
  h2dummy = new TH2D("h2dummy","",1,0.0,3.0,1,MIN,MAX);
  h2dummy->Draw();
  line.Draw();
  h2dummy->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h2dummy->GetYaxis()->SetTitle(Form("v_{%d} not corrected for EP resolution",harmonic));
  h2dummy->GetYaxis()->SetTitleOffset(1.3);

  hvn_fvtxn_B->Draw("same");
  hvn_fvtxn_E->Draw("same");
  hvn_fvtxn_W->Draw("same");

  TLegend* leg = new TLegend(0.18,0.68,0.38,0.88);
  leg->SetHeader(Form("%d GeV, FVTXN",energy));
  leg->AddEntry(hvn_fvtxn_B,"both","el");
  leg->AddEntry(hvn_fvtxn_E,"east","el");
  leg->AddEntry(hvn_fvtxn_W,"west","el");
  leg->SetTextSize(0.05);
  leg->Draw();

  c1->Print(Form("FigsHarmonicCoefficient/diagnostic_fvtxn_EBW_energy%d_harm%d.png",energy,harmonic));
  c1->Print(Form("FigsHarmonicCoefficient/diagnostic_fvtxn_EBW_energy%d_harm%d.pdf",energy,harmonic));



  delete c1;

}
