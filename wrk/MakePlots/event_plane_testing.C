void cospsidiff_run(int, int);
void psidiff_run(int, int);

void beforeafterflat_run(int, int);


void event_plane_testing()
{

  psidiff_run(455050,2);
  psidiff_run(456015,2);
  psidiff_run(457015,2);
  psidiff_run(457043,2);
  psidiff_run(457213,2);
  psidiff_run(458014,2);
  psidiff_run(458153,2);
  psidiff_run(458167,2);

  cospsidiff_run(455050,2);
  cospsidiff_run(456015,2);
  cospsidiff_run(457015,2);
  cospsidiff_run(457043,2);
  cospsidiff_run(457213,2);
  cospsidiff_run(458014,2);
  cospsidiff_run(458153,2);
  cospsidiff_run(458167,2);

  beforeafterflat_run(455050,2);
  beforeafterflat_run(456015,2);
  beforeafterflat_run(457015,2);
  beforeafterflat_run(457043,2);
  beforeafterflat_run(457213,2);
  beforeafterflat_run(458014,2);
  beforeafterflat_run(458153,2);
  beforeafterflat_run(458167,2);

}


void beforeafterflat_run(int run, int harmonic)
{

  gStyle->SetOptTitle(1);

  TCanvas* c1 = new TCanvas("c1","");

  TFile* file = TFile::Open(Form("input/hrp_%d.root",run));

  TH1D* hbbc = (TH1D*) ((TH2D*)file->Get("psi_bf_0_1_2"))->ProjectionY();
  hbbc->Draw();
  hbbc->GetXaxis()->SetTitle("#Psi_{2}, BBCS");
  c1->Print(Form("FigsEventPlane/before_psi2_bbcs_%d.png",run));
  c1->Print(Form("FigsEventPlane/before_psi2_bbcs_%d.pdf",run));
  TH1D* hbbcA = (TH1D*) ((TH2D*)file->Get("psi_af_0_1_2"))->ProjectionY();
  hbbcA->Draw("same");
  hbbcA->GetXaxis()->SetTitle("#Psi_{2}, BBCS");
  c1->Print(Form("FigsEventPlane/after_psi2_bbcs_%d.png",run));
  c1->Print(Form("FigsEventPlane/after_psi2_bbcs_%d.pdf",run));

  TH1D* hfvtxs = (TH1D*) ((TH2D*)file->Get("psi_bf_0_1_3"))->ProjectionY();
  hfvtxs->Draw();
  hfvtxs->GetXaxis()->SetTitle("#Psi_{2}, FVTXS");
  c1->Print(Form("FigsEventPlane/before_psi2_fvtxs_%d.png",run));
  c1->Print(Form("FigsEventPlane/before_psi2_fvtxs_%d.pdf",run));
  TH1D* hfvtxsA = (TH1D*) ((TH2D*)file->Get("psi_af_0_1_3"))->ProjectionY();
  hfvtxsA->Draw("same");
  hfvtxsA->GetXaxis()->SetTitle("#Psi_{2}, FVTXS");
  c1->Print(Form("FigsEventPlane/after_psi2_fvtxs_%d.png",run));
  c1->Print(Form("FigsEventPlane/after_psi2_fvtxs_%d.pdf",run));

  TH1D* hfvtxn = (TH1D*) ((TH2D*)file->Get("psi_bf_0_1_8"))->ProjectionY();
  hfvtxn->Draw();
  hfvtxn->GetXaxis()->SetTitle("#Psi_{2}, FVTXN");
  c1->Print(Form("FigsEventPlane/before_psi2_fvtxn_%d.png",run));
  c1->Print(Form("FigsEventPlane/before_psi2_fvtxn_%d.pdf",run));
  TH1D* hfvtxnA = (TH1D*) ((TH2D*)file->Get("psi_af_0_1_8"))->ProjectionY();
  hfvtxnA->Draw("same");
  hfvtxnA->GetXaxis()->SetTitle("#Psi_{2}, FVTXN");
  c1->Print(Form("FigsEventPlane/after_psi2_fvtxn_%d.png",run));
  c1->Print(Form("FigsEventPlane/after_psi2_fvtxn_%d.pdf",run));

  delete c1;

}

void psidiff_run(int run, int harmonic)
{

  gStyle->SetOptTitle(1);

  TCanvas* c1 = new TCanvas("c1","");

  TFile* file = TFile::Open(Form("input/hist_%d.root",run));

  // ---

  TH1D* th1d_dreso_bbc_cnt = (TH1D*)file->Get(Form("th1d_dreso%d_BBC_CNT",harmonic));
  TH1D* th1d_dreso_bbc_fvtx = (TH1D*)file->Get(Form("th1d_dreso%d_BBC_FVTX",harmonic));
  TH1D* th1d_dreso_cnt_fvtx = (TH1D*)file->Get(Form("th1d_dreso%d_CNT_FVTX",harmonic));
  TH1D* th1d_dreso_bbc_fvtxn = (TH1D*)file->Get(Form("th1d_dreso%d_BBC_FVTXN",harmonic));
  TH1D* th1d_dreso_cnt_fvtxn = (TH1D*)file->Get(Form("th1d_dreso%d_CNT_FVTXN",harmonic));
  TH1D* th1d_dreso_fvtxs_fvtxn = (TH1D*)file->Get(Form("th1d_dreso%d_FVTXS_FVTXN",harmonic));

  th1d_dreso_bbc_cnt->SetLineColor(kBlue);
  th1d_dreso_bbc_fvtx->SetLineColor(kBlack);
  th1d_dreso_cnt_fvtx->SetLineColor(kRed);
  th1d_dreso_bbc_fvtxn->SetLineColor(kGreen+2);
  th1d_dreso_cnt_fvtxn->SetLineColor(kOrange+2);
  th1d_dreso_fvtxs_fvtxn->SetLineColor(kGray);

  double max = 0;
  double max_BC = 0;
  double max_BN = 0;
  double max_BS = 0;
  double max_CN = 0;
  double max_CS = 0;
  double max_NS = 0;

  max_BC = th1d_dreso_bbc_cnt->GetMaximum();
  max_BN = th1d_dreso_bbc_fvtxn->GetMaximum();
  max_BS = th1d_dreso_bbc_fvtx->GetMaximum();
  max_CN = th1d_dreso_cnt_fvtxn->GetMaximum();
  max_CS = th1d_dreso_cnt_fvtx->GetMaximum();
  max_NS = th1d_dreso_fvtxs_fvtxn->GetMaximum();

  max = 1;
  if ( max_BC > max ) max = max_BC;
  if ( max_BN > max ) max = max_BN;
  if ( max_BS > max ) max = max_BS;
  if ( max_CN > max ) max = max_CN;
  if ( max_CS > max ) max = max_CS;
  if ( max_NS > max ) max = max_NS;

  TH2* th2 = new TH2D("th2","",1,-6.0,6.0,1,0.0,1.1*max);
  th2->Draw();
  th1d_dreso_bbc_cnt->Draw("same");
  th1d_dreso_bbc_fvtx->Draw("same");
  th1d_dreso_cnt_fvtx->Draw("same");
  th1d_dreso_bbc_fvtxn->Draw("same");
  th1d_dreso_cnt_fvtxn->Draw("same");
  th1d_dreso_fvtxs_fvtxn->Draw("same");
  TLegend* leg = new TLegend();
  c1->Print(Form("FigsEventPlane/psidiff_dreso%d_%d.png",harmonic,run));
  c1->Print(Form("FigsEventPlane/psidiff_dreso%d_%d.pdf",harmonic,run));
  delete th2;
  delete leg;

  // ---

  th1d_dreso_bbc_fvtx->SetLineColor(kRed);
  th1d_dreso_bbc_fvtxn->SetLineColor(kGreen+2);
  th1d_dreso_fvtxs_fvtxn->SetLineColor(kBlue);

  th1d_dreso_bbc_cnt->SetLineColor(kBlue);
  th1d_dreso_cnt_fvtx->SetLineColor(kRed);
  th1d_dreso_cnt_fvtxn->SetLineColor(kBlack);

  // ---

  max = 1;
  if ( max_BN > max ) max = max_BN;
  if ( max_BS > max ) max = max_BS;
  if ( max_NS > max ) max = max_NS;

  th2 = new TH2D("th2","",1,-6.0,6.0,1,0.0,1.1*max);
  th2->Draw();
  th1d_dreso_bbc_fvtx->Draw("same");
  th1d_dreso_bbc_fvtxn->Draw("same");
  th1d_dreso_fvtxs_fvtxn->Draw("same");
  leg = new TLegend(0.18,0.68,0.38,0.88);
  leg->SetHeader(Form("Run %d #Psi_{%d}",run,harmonic));
  leg->AddEntry(th1d_dreso_bbc_fvtx,"BBCS-FVTXS","l");
  leg->AddEntry(th1d_dreso_bbc_fvtxn,"BBCS-FVTXN","l");
  leg->AddEntry(th1d_dreso_fvtxs_fvtxn,"FVTXS-FVTXN","l");
  leg->SetTextSize(0.045);
  leg->Draw();
  c1->Print(Form("FigsEventPlane/psidiff_dreso%d_%d_SANSCNT.png",harmonic,run));
  c1->Print(Form("FigsEventPlane/psidiff_dreso%d_%d_SANSCNT.pdf",harmonic,run));
  delete th2;

  max = 1;
  if ( max_BC > max ) max = max_BC;
  if ( max_CN > max ) max = max_CN;
  if ( max_CS > max ) max = max_CS;

  th2 = new TH2D("th2","",1,-6.0,6.0,1,0.0,1.1*max);
  th2->Draw();
  th1d_dreso_bbc_cnt->Draw("same");
  th1d_dreso_cnt_fvtx->Draw("same");
  th1d_dreso_cnt_fvtxn->Draw("same");
  leg = new TLegend(0.18,0.68,0.38,0.88);
  leg->SetHeader(Form("Run %d #Psi_{%d}",run,harmonic));
  leg->AddEntry(th1d_dreso_bbc_cnt,"BBCS-CNT","l");
  leg->AddEntry(th1d_dreso_cnt_fvtx,"CNT-FVTXS","l");
  leg->AddEntry(th1d_dreso_cnt_fvtxn,"CNT-FVTXN","l");
  leg->SetTextSize(0.045);
  leg->Draw();
  c1->Print(Form("FigsEventPlane/psidiff_dreso%d_%d_AVECCNT.png",harmonic,run));
  c1->Print(Form("FigsEventPlane/psidiff_dreso%d_%d_AVECCNT.pdf",harmonic,run));
  delete th2;



  delete c1;

}



void cospsidiff_run(int run, int harmonic)
{

  gStyle->SetOptTitle(1);

  TCanvas* c1 = new TCanvas("c1","");

  TFile* file = TFile::Open(Form("input/hist_%d.root",run));

  // ---

  TH1D* th1d_reso_bbc_cnt = (TH1D*)file->Get(Form("th1d_reso%d_BBC_CNT",harmonic));
  TH1D* th1d_reso_bbc_fvtx = (TH1D*)file->Get(Form("th1d_reso%d_BBC_FVTX",harmonic));
  TH1D* th1d_reso_cnt_fvtx = (TH1D*)file->Get(Form("th1d_reso%d_CNT_FVTX",harmonic));
  TH1D* th1d_reso_bbc_fvtxn = (TH1D*)file->Get(Form("th1d_reso%d_BBC_FVTXN",harmonic));
  TH1D* th1d_reso_cnt_fvtxn = (TH1D*)file->Get(Form("th1d_reso%d_CNT_FVTXN",harmonic));
  TH1D* th1d_reso_fvtxs_fvtxn = (TH1D*)file->Get(Form("th1d_reso%d_FVTXS_FVTXN",harmonic));

  th1d_reso_bbc_cnt->SetLineColor(kBlue);
  th1d_reso_bbc_fvtx->SetLineColor(kBlack);
  th1d_reso_cnt_fvtx->SetLineColor(kRed);
  th1d_reso_bbc_fvtxn->SetLineColor(kGreen+2);
  th1d_reso_cnt_fvtxn->SetLineColor(kOrange+2);
  th1d_reso_fvtxs_fvtxn->SetLineColor(kGray);

  double max = 0;
  double max_BC = 0;
  double max_BN = 0;
  double max_BS = 0;
  double max_CN = 0;
  double max_CS = 0;
  double max_NS = 0;

  max_BC = th1d_reso_bbc_cnt->GetMaximum();
  max_BN = th1d_reso_bbc_fvtxn->GetMaximum();
  max_BS = th1d_reso_bbc_fvtx->GetMaximum();
  max_CN = th1d_reso_cnt_fvtxn->GetMaximum();
  max_CS = th1d_reso_cnt_fvtx->GetMaximum();
  max_NS = th1d_reso_fvtxs_fvtxn->GetMaximum();

  max = 1;
  if ( max_BC > max ) max = max_BC;
  if ( max_BN > max ) max = max_BN;
  if ( max_BS > max ) max = max_BS;
  if ( max_CN > max ) max = max_CN;
  if ( max_CS > max ) max = max_CS;
  if ( max_NS > max ) max = max_NS;

  TH2* th2 = new TH2D("th2","",1,-1.1,1.1,1,0.0,1.1*max);
  th2->Draw();
  th1d_reso_bbc_cnt->Draw("same");
  th1d_reso_bbc_fvtx->Draw("same");
  th1d_reso_cnt_fvtx->Draw("same");
  th1d_reso_bbc_fvtxn->Draw("same");
  th1d_reso_cnt_fvtxn->Draw("same");
  th1d_reso_fvtxs_fvtxn->Draw("same");
  TLegend* leg = new TLegend();
  c1->Print(Form("FigsEventPlane/psidiff_reso%d_%d.png",harmonic,run));
  c1->Print(Form("FigsEventPlane/psidiff_reso%d_%d.pdf",harmonic,run));
  delete th2;
  delete leg;

  // ---

  th1d_reso_bbc_fvtx->SetLineColor(kRed);
  th1d_reso_bbc_fvtxn->SetLineColor(kGreen+2);
  th1d_reso_fvtxs_fvtxn->SetLineColor(kBlue);

  th1d_reso_bbc_cnt->SetLineColor(kBlue);
  th1d_reso_cnt_fvtx->SetLineColor(kRed);
  th1d_reso_cnt_fvtxn->SetLineColor(kBlack);

  // ---

  max = 1;
  if ( max_BN > max ) max = max_BN;
  if ( max_BS > max ) max = max_BS;
  if ( max_NS > max ) max = max_NS;

  th2 = new TH2D("th2","",1,-1.1,1.1,1,0.0,1.1*max);
  th2->Draw();
  th1d_reso_bbc_fvtx->Draw("same");
  th1d_reso_bbc_fvtxn->Draw("same");
  th1d_reso_fvtxs_fvtxn->Draw("same");
  leg = new TLegend(0.25,0.68,0.38,0.88);
  leg->SetHeader(Form("Run %d #Psi_{%d}",run,harmonic));
  leg->AddEntry(th1d_reso_bbc_fvtx,"BBCS-FVTXS","l");
  leg->AddEntry(th1d_reso_bbc_fvtxn,"BBCS-FVTXN","l");
  leg->AddEntry(th1d_reso_fvtxs_fvtxn,"FVTXS-FVTXN","l");
  leg->SetTextSize(0.045);
  leg->Draw();
  c1->Print(Form("FigsEventPlane/psidiff_reso%d_%d_SANSCNT.png",harmonic,run));
  c1->Print(Form("FigsEventPlane/psidiff_reso%d_%d_SANSCNT.pdf",harmonic,run));
  delete th2;

  max = 1;
  if ( max_BC > max ) max = max_BC;
  if ( max_CN > max ) max = max_CN;
  if ( max_CS > max ) max = max_CS;

  th2 = new TH2D("th2","",1,-1.1,1.1,1,0.0,1.1*max);
  th2->Draw();
  th1d_reso_bbc_cnt->Draw("same");
  th1d_reso_cnt_fvtx->Draw("same");
  th1d_reso_cnt_fvtxn->Draw("same");
  leg = new TLegend(0.25,0.68,0.38,0.88);
  leg->SetHeader(Form("Run %d #Psi_{%d}",run,harmonic));
  leg->AddEntry(th1d_reso_bbc_cnt,"BBCS-CNT","l");
  leg->AddEntry(th1d_reso_cnt_fvtx,"CNT-FVTXS","l");
  leg->AddEntry(th1d_reso_cnt_fvtxn,"CNT-FVTXN","l");
  leg->SetTextSize(0.045);
  leg->Draw();
  c1->Print(Form("FigsEventPlane/psidiff_reso%d_%d_AVECCNT.png",harmonic,run));
  c1->Print(Form("FigsEventPlane/psidiff_reso%d_%d_AVECCNT.pdf",harmonic,run));
  delete th2;



  delete c1;

}
