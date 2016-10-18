void doit(int, int);



void temp_clustrack()
{



  // doit(456652,2);
  doit(200,2);
  return;
  for ( int i = 2; i < 4; ++i )
    {
      doit(200,i);
      doit(62,i);
      doit(39,i);
      doit(20,i);
    }


}

void doit(int handle, int harmonic)
{

  gStyle->SetOptTitle(1);

  TCanvas* c1 = new TCanvas("c1","");

  TFile* file = NULL;
  if ( handle <= 200 ) file = TFile::Open(Form("input/combined_%d.root",handle));
  else if ( handle > 454000 ) file = TFile::Open(Form("input/hist_%d.root",handle));
  else
    {
      cout << "YOU'RE GONNA DIE" << endl;
      return;
    }

  ifstream finep;
  if ( handle <= 200 ) finep.open((const char*)Form("DataTextFiles/run16dau%d_v%dfvtxs.dat",handle,harmonic));
  float pt[15], epv2[15], eepv2[15], esysepv2[15];
  for ( int i = 0; i < 15; ++i )
    {
      if ( finep.is_open() ) finep>>pt[i]>>epv2[i]>>eepv2[i];
      if ( epv2[i] <= 0 ) epv2[i] = -9;
    }
  TGraphErrors* tge_ep = new TGraphErrors(15,pt,epv2,0,eepv2);
  tge_ep->SetMarkerStyle(kFullCircle);

  // ---
  // ---
  // ---

  // ------------
  // --- BBCS ---
  // ------------

  TProfile* tp1f_os_bbcs_dn2_both = (TProfile*)file->Get(Form("os_bbcs_d%d2_both",harmonic));
  TProfile* tp1f_os_bbcs_cn2 = (TProfile*)file->Get(Form("os_bbcs_c%d2",harmonic));
  double os_bbcs_cn2_raw = tp1f_os_bbcs_cn2->GetBinContent(1);
  double os_bbcs_cn2_corr = os_bbcs_cn2_raw;
  if ( os_bbcs_cn2_corr < 0 ) cout << "YOU'RE GONNA DIE" << endl;
  double os_bbcs_vn2_corr = sqrt(os_bbcs_cn2_corr);
  cout << os_bbcs_vn2_corr << endl;

  TH1D* th1d_os_bbcs_dn2_both = (TH1D*)tp1f_os_bbcs_dn2_both->ProjectionX();
  TH1D* th1d_os_bbcs_dn2_both_corr = (TH1D*)th1d_os_bbcs_dn2_both->Clone();
  th1d_os_bbcs_dn2_both_corr->Draw();
  c1->Print(Form("FigsClusTrack/corr_os_bbcs_d%d2_%d.png",harmonic,handle));
  c1->Print(Form("FigsClusTrack/corr_os_bbcs_d%d2_%d.pdf",harmonic,handle));

  th1d_os_bbcs_dn2_both_corr->Scale(1.0/os_bbcs_vn2_corr);
  if ( handle <= 200 ) th1d_os_bbcs_dn2_both_corr->SetTitle(Form("d+Au #sqrt{s_{NN}} = %d GeV",handle));
  if ( handle == 62 ) th1d_os_bbcs_dn2_both_corr->SetTitle(Form("d+Au #sqrt{s_{NN}} = 62.4 GeV"));
  if ( handle == 20 ) th1d_os_bbcs_dn2_both_corr->SetTitle(Form("d+Au #sqrt{s_{NN}} = 19.6 GeV"));
  th1d_os_bbcs_dn2_both_corr->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  th1d_os_bbcs_dn2_both_corr->GetYaxis()->SetTitle(Form("v_{%d}{2}",harmonic));
  //  th1d_os_bbcs_dn2_both_corr->SetMarkerStyle(kOpenCircle);
  th1d_os_bbcs_dn2_both_corr->SetMaximum(th1d_os_bbcs_dn2_both_corr->GetMaximum()*1.2);
  th1d_os_bbcs_dn2_both_corr->Draw();
  float tlx = 0.65;
  float tly = 0.2;
  TLatex* tlref = new TLatex(tlx,tly,Form("#sqrt{c_{%d}{2}} = %.3f",harmonic,os_bbcs_vn2_corr));
  tlref->SetTextSize(0.05);
  tlref->SetNDC();
  tlref->Draw();
  float tlx2 = 0.17;
  float tly2 = 0.85;
  TLatex* tlmeth = new TLatex(tlx2,tly2,Form("scalar product with BBCS"));
  tlmeth->SetTextSize(0.05);
  tlmeth->SetNDC();
  tlmeth->Draw();
  float tlx3 = 0.17;
  float tly3 = 0.8;
  TLatex* tlmeth3 = new TLatex(tlx3,tly3,Form("reference flow using BBCS"));
  tlmeth3->SetTextSize(0.05);
  tlmeth3->SetNDC();
  tlmeth3->Draw();
  c1->Print(Form("FigsClusTrack/corr_os_bbcs_v%d2_%d.png",harmonic,handle));
  c1->Print(Form("FigsClusTrack/corr_os_bbcs_v%d2_%d.pdf",harmonic,handle));



  // -------------
  // --- FVTXS ---
  // -------------

  TProfile* tp1f_os_fvtxs_dn2_both = (TProfile*)file->Get(Form("os_fvtxs_d%d2_both",harmonic));
  TProfile* tp1f_os_fvtxs_cn2 = (TProfile*)file->Get(Form("os_fvtxs_c%d2",harmonic));
  double os_fvtxs_cn2_raw = tp1f_os_fvtxs_cn2->GetBinContent(1);
  double os_fvtxs_cn2_corr = os_fvtxs_cn2_raw;
  if ( os_fvtxs_cn2_corr < 0 ) cout << "YOU'RE GONNA DIE" << endl;
  double os_fvtxs_vn2_corr = sqrt(os_fvtxs_cn2_corr);
  cout << os_fvtxs_vn2_corr << endl;

  TH1D* th1d_os_fvtxs_dn2_both = (TH1D*)tp1f_os_fvtxs_dn2_both->ProjectionX();
  TH1D* th1d_os_fvtxs_dn2_both_corr = (TH1D*)th1d_os_fvtxs_dn2_both->Clone();
  //  th1d_os_fvtxs_dn2_both_corr->SetMarkerStyle(kOpenCircle);
  th1d_os_fvtxs_dn2_both_corr->Draw();
  c1->Print(Form("FigsClusTrack/corr_os_fvtxs_d%d2_%d.png",harmonic,handle));
  c1->Print(Form("FigsClusTrack/corr_os_fvtxs_d%d2_%d.pdf",harmonic,handle));

  th1d_os_fvtxs_dn2_both_corr->Scale(1.0/os_fvtxs_vn2_corr);
  if ( handle <= 200 ) th1d_os_fvtxs_dn2_both_corr->SetTitle(Form("d+Au #sqrt{s_{NN}} = %d GeV",handle));
  th1d_os_fvtxs_dn2_both_corr->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  th1d_os_fvtxs_dn2_both_corr->GetYaxis()->SetTitle("v_{2}{2}");
  //  th1d_os_fvtxs_dn2_both_corr->SetMarkerStyle(kOpenCircle);
  th1d_os_fvtxs_dn2_both_corr->SetMaximum(th1d_os_fvtxs_dn2_both_corr->GetMaximum()*1.2);
  th1d_os_fvtxs_dn2_both_corr->Draw();
  tlref->DrawLatex(tlx,tly,Form("#sqrt{c_{%d}{2}} = %.3f",harmonic,os_fvtxs_vn2_corr));
  tlmeth->DrawLatex(tlx2,tly2,Form("scalar product with FVTXS"));
  tlmeth3->DrawLatex(tlx3,tly3,Form("reference flow using FVTXS"));
  c1->Print(Form("FigsClusTrack/corr_os_fvtxs_v%d2_%d.png",harmonic,handle));
  c1->Print(Form("FigsClusTrack/corr_os_fvtxs_v%d2_%d.pdf",harmonic,handle));
  th1d_os_fvtxs_dn2_both_corr->Scale(os_fvtxs_vn2_corr); // remove old scaling
  TProfile* tp1f_os_fvtxs_ce01_cn2 = (TProfile*)file->Get(Form("os_fvtxs_ce01_c%d2",harmonic));
  double os_fvtxs_ce01_cn2_raw = tp1f_os_fvtxs_ce01_cn2->GetBinContent(1);
  double os_fvtxs_ce01_cn2_corr = os_fvtxs_ce01_cn2_raw;
  if ( os_fvtxs_ce01_cn2_corr < 0 ) cout << "YOU'RE GONNA DIE" << endl;
  double os_fvtxs_ce01_vn2_corr = sqrt(os_fvtxs_ce01_cn2_corr);
  cout << os_fvtxs_ce01_vn2_corr << endl;
  th1d_os_fvtxs_dn2_both_corr->Scale(1.0/os_fvtxs_ce01_vn2_corr); // apply new scaling
  th1d_os_fvtxs_dn2_both_corr->Draw();
  tlref->DrawLatex(tlx,tly,Form("#sqrt{c_{%d}{2}} = %.3f",harmonic,os_fvtxs_ce01_vn2_corr));
  tlmeth->DrawLatex(tlx2,tly2,Form("scalar product with FVTXS"));
  tlmeth3->DrawLatex(tlx3,tly3,Form("reference flow using FVTXS"));
  c1->Print(Form("FigsClusTrack/corr_os_fvtxs_ce01_v%d2_%d.png",harmonic,handle));
  c1->Print(Form("FigsClusTrack/corr_os_fvtxs_ce01_v%d2_%d.pdf",harmonic,handle));



  // -------------
  // --- FVTXS_TRACKS ---
  // -------------

  TProfile* tp1f_os_fvtxs_tracks_dn2_both = (TProfile*)file->Get(Form("os_fvtxs_tracks_d%d2_both",harmonic));
  TProfile* tp1f_os_fvtxs_tracks_cn2 = (TProfile*)file->Get(Form("os_fvtxs_tracks_c%d2",harmonic));
  double os_fvtxs_tracks_cn2_raw = tp1f_os_fvtxs_tracks_cn2->GetBinContent(1);
  double os_fvtxs_tracks_cn2_corr = os_fvtxs_tracks_cn2_raw;
  if ( os_fvtxs_tracks_cn2_corr < 0 ) cout << "YOU'RE GONNA DIE" << endl;
  double os_fvtxs_tracks_vn2_corr = sqrt(os_fvtxs_tracks_cn2_corr);
  cout << os_fvtxs_tracks_vn2_corr << endl;

  TH1D* th1d_os_fvtxs_tracks_dn2_both = (TH1D*)tp1f_os_fvtxs_tracks_dn2_both->ProjectionX();
  TH1D* th1d_os_fvtxs_tracks_dn2_both_corr = (TH1D*)th1d_os_fvtxs_tracks_dn2_both->Clone();
  //  th1d_os_fvtxs_tracks_dn2_both_corr->SetMarkerStyle(kOpenCircle);
  th1d_os_fvtxs_tracks_dn2_both_corr->Draw();
  c1->Print(Form("FigsClusTrack/corr_os_fvtxs_tracks_d%d2_%d.png",harmonic,handle));
  c1->Print(Form("FigsClusTrack/corr_os_fvtxs_tracks_d%d2_%d.pdf",harmonic,handle));

  th1d_os_fvtxs_tracks_dn2_both_corr->Scale(1.0/os_fvtxs_tracks_vn2_corr);
  if ( handle <= 200 ) th1d_os_fvtxs_tracks_dn2_both_corr->SetTitle(Form("d+Au #sqrt{s_{NN}} = %d GeV",handle));
  th1d_os_fvtxs_tracks_dn2_both_corr->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  th1d_os_fvtxs_tracks_dn2_both_corr->GetYaxis()->SetTitle("v_{2}{2}");
  //  th1d_os_fvtxs_tracks_dn2_both_corr->SetMarkerStyle(kOpenCircle);
  th1d_os_fvtxs_tracks_dn2_both_corr->SetMaximum(th1d_os_fvtxs_tracks_dn2_both_corr->GetMaximum()*1.2);
  th1d_os_fvtxs_tracks_dn2_both_corr->Draw();
  tlref->DrawLatex(tlx,tly,Form("#sqrt{c_{%d}{2}} = %.3f",harmonic,os_fvtxs_tracks_vn2_corr));
  tlmeth->DrawLatex(tlx2,tly2,Form("scalar product with FVTXS_TRACKS"));
  tlmeth3->DrawLatex(tlx3,tly3,Form("reference flow using FVTXS_TRACKS"));
  c1->Print(Form("FigsClusTrack/corr_os_fvtxs_tracks_v%d2_%d.png",harmonic,handle));
  c1->Print(Form("FigsClusTrack/corr_os_fvtxs_tracks_v%d2_%d.pdf",harmonic,handle));

  // -------------
  // --- FVTXN ---
  // -------------

  TProfile* tp1f_os_fvtxn_dn2_both = (TProfile*)file->Get(Form("os_fvtxn_d%d2_both",harmonic));
  TProfile* tp1f_os_fvtxn_cn2 = (TProfile*)file->Get(Form("os_fvtxn_c%d2",harmonic));
  double os_fvtxn_cn2_raw = tp1f_os_fvtxn_cn2->GetBinContent(1);
  double os_fvtxn_cn2_corr = os_fvtxn_cn2_raw;
  if ( os_fvtxn_cn2_corr < 0 ) cout << "YOU'RE GONNA DIE" << endl;
  double os_fvtxn_vn2_corr = sqrt(os_fvtxn_cn2_corr);
  cout << os_fvtxn_vn2_corr << endl;

  TH1D* th1d_os_fvtxn_dn2_both = (TH1D*)tp1f_os_fvtxn_dn2_both->ProjectionX();
  TH1D* th1d_os_fvtxn_dn2_both_corr = (TH1D*)th1d_os_fvtxn_dn2_both->Clone();
  //  th1d_os_fvtxn_dn2_both_corr->SetMarkerStyle(kOpenCircle);
  th1d_os_fvtxn_dn2_both_corr->SetMaximum(th1d_os_fvtxn_dn2_both_corr->GetMaximum()*1.2);
  th1d_os_fvtxn_dn2_both_corr->Draw();
  c1->Print(Form("FigsClusTrack/corr_os_fvtxn_d%d2_%d.png",harmonic,handle));
  c1->Print(Form("FigsClusTrack/corr_os_fvtxn_d%d2_%d.pdf",harmonic,handle));

  th1d_os_fvtxn_dn2_both_corr->Scale(1.0/os_fvtxn_vn2_corr);
  if ( handle <= 200 ) th1d_os_fvtxn_dn2_both_corr->SetTitle(Form("d+Au #sqrt{s_{NN}} = %d GeV",handle));
  th1d_os_fvtxn_dn2_both_corr->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  th1d_os_fvtxn_dn2_both_corr->GetYaxis()->SetTitle("v_{2}{2}");
  //  th1d_os_fvtxn_dn2_both_corr->SetMarkerStyle(kOpenCircle);
  th1d_os_fvtxn_dn2_both_corr->SetMaximum(th1d_os_fvtxn_dn2_both_corr->GetMaximum()*1.2);
  th1d_os_fvtxn_dn2_both_corr->Draw();
  tlref->DrawLatex(tlx,tly,Form("#sqrt{c_{%d}{2}} = %.3f",harmonic,os_fvtxn_vn2_corr));
  tlmeth->DrawLatex(tlx2,tly2,Form("scalar product with FVTXN"));
  tlmeth3->DrawLatex(tlx3,tly3,Form("reference flow using FVTXN"));
  c1->Print(Form("FigsClusTrack/corr_os_fvtxn_v%d2_%d.png",harmonic,handle));
  c1->Print(Form("FigsClusTrack/corr_os_fvtxn_v%d2_%d.pdf",harmonic,handle));

  // -------------
  // --- FVTXN_TRACKS ---
  // -------------

  TProfile* tp1f_os_fvtxn_tracks_dn2_both = (TProfile*)file->Get(Form("os_fvtxn_tracks_d%d2_both",harmonic));
  TProfile* tp1f_os_fvtxn_tracks_cn2 = (TProfile*)file->Get(Form("os_fvtxn_tracks_c%d2",harmonic));
  double os_fvtxn_tracks_cn2_raw = tp1f_os_fvtxn_tracks_cn2->GetBinContent(1);
  double os_fvtxn_tracks_cn2_corr = os_fvtxn_tracks_cn2_raw;
  if ( os_fvtxn_tracks_cn2_corr < 0 ) cout << "YOU'RE GONNA DIE" << endl;
  double os_fvtxn_tracks_vn2_corr = sqrt(os_fvtxn_tracks_cn2_corr);
  cout << os_fvtxn_tracks_vn2_corr << endl;

  TH1D* th1d_os_fvtxn_tracks_dn2_both = (TH1D*)tp1f_os_fvtxn_tracks_dn2_both->ProjectionX();
  TH1D* th1d_os_fvtxn_tracks_dn2_both_corr = (TH1D*)th1d_os_fvtxn_tracks_dn2_both->Clone();
  //  th1d_os_fvtxn_tracks_dn2_both_corr->SetMarkerStyle(kOpenCircle);
  th1d_os_fvtxn_tracks_dn2_both_corr->SetMaximum(th1d_os_fvtxn_tracks_dn2_both_corr->GetMaximum()*1.2);
  th1d_os_fvtxn_tracks_dn2_both_corr->Draw();
  c1->Print(Form("FigsClusTrack/corr_os_fvtxn_tracks_d%d2_%d.png",harmonic,handle));
  c1->Print(Form("FigsClusTrack/corr_os_fvtxn_tracks_d%d2_%d.pdf",harmonic,handle));

  th1d_os_fvtxn_tracks_dn2_both_corr->Scale(1.0/os_fvtxn_tracks_vn2_corr);
  if ( handle <= 200 ) th1d_os_fvtxn_tracks_dn2_both_corr->SetTitle(Form("d+Au #sqrt{s_{NN}} = %d GeV",handle));
  th1d_os_fvtxn_tracks_dn2_both_corr->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  th1d_os_fvtxn_tracks_dn2_both_corr->GetYaxis()->SetTitle("v_{2}{2}");
  //  th1d_os_fvtxn_tracks_dn2_both_corr->SetMarkerStyle(kOpenCircle);
  th1d_os_fvtxn_tracks_dn2_both_corr->SetMaximum(th1d_os_fvtxn_tracks_dn2_both_corr->GetMaximum()*1.2);
  th1d_os_fvtxn_tracks_dn2_both_corr->Draw();
  tlref->DrawLatex(tlx,tly,Form("#sqrt{c_{%d}{2}} = %.3f",harmonic,os_fvtxn_tracks_vn2_corr));
  tlmeth->DrawLatex(tlx2,tly2,Form("scalar product with FVTXN_TRACKS"));
  tlmeth3->DrawLatex(tlx3,tly3,Form("reference flow using FVTXN_TRACKS"));
  c1->Print(Form("FigsClusTrack/corr_os_fvtxn_tracks_v%d2_%d.png",harmonic,handle));
  c1->Print(Form("FigsClusTrack/corr_os_fvtxn_tracks_v%d2_%d.pdf",harmonic,handle));

  // ---

  // -------------------------------------------------
  // --- different method of estimating reference flow
  // -------------------------------------------------

  // --- note that the bbcs-fvtxn correlation is so weak in the lower energies that this method might not work

  TProfile* tp1f_os_bbcsfvtxs_cn2 = (TProfile*)file->Get(Form("os_bbcsfvtxs_c%d2",harmonic));
  TProfile* tp1f_os_bbcsfvtxn_cn2 = (TProfile*)file->Get(Form("os_bbcsfvtxn_c%d2",harmonic));
  TProfile* tp1f_os_fvtxsfvtxn_cn2 = (TProfile*)file->Get(Form("os_fvtxsfvtxn_c%d2",harmonic));

  double os_bbcsfvtxs_cn2_raw = tp1f_os_bbcsfvtxs_cn2->GetBinContent(1);
  double os_bbcsfvtxn_cn2_raw = tp1f_os_bbcsfvtxn_cn2->GetBinContent(1);
  double os_fvtxsfvtxn_cn2_raw = tp1f_os_fvtxsfvtxn_cn2->GetBinContent(1);

  double os_newref_cn2_bbcs = (os_bbcsfvtxs_cn2_raw*os_bbcsfvtxn_cn2_raw)/os_fvtxsfvtxn_cn2_raw;
  double os_newref_cn2_fvtxs = (os_bbcsfvtxs_cn2_raw*os_fvtxsfvtxn_cn2_raw)/os_bbcsfvtxn_cn2_raw;
  double os_newref_cn2_fvtxn = (os_bbcsfvtxn_cn2_raw*os_fvtxsfvtxn_cn2_raw)/os_bbcsfvtxs_cn2_raw;

  if ( os_newref_cn2_bbcs < 0 ) cout << "YOU'RE GONNA DIE" << endl;
  if ( os_newref_cn2_fvtxs < 0 ) cout << "YOU'RE GONNA DIE" << endl;
  if ( os_newref_cn2_fvtxn < 0 ) cout << "YOU'RE GONNA DIE" << endl;

  double os_newref_vn2_bbcs =  sqrt(os_newref_cn2_bbcs);
  double os_newref_vn2_fvtxs = sqrt(os_newref_cn2_fvtxs);
  double os_newref_vn2_fvtxn = sqrt(os_newref_cn2_fvtxn);

  cout << os_newref_vn2_bbcs << endl;
  cout << os_newref_vn2_fvtxs << endl;
  cout << os_newref_vn2_fvtxn << endl;

  // --- rescale to remove the old v2 and apply the new v2
  th1d_os_bbcs_dn2_both_corr->Scale(os_bbcs_vn2_corr/os_newref_vn2_bbcs);
  th1d_os_fvtxs_dn2_both_corr->Scale(os_fvtxs_vn2_corr/os_newref_vn2_fvtxs);
  th1d_os_fvtxn_dn2_both_corr->Scale(os_fvtxn_vn2_corr/os_newref_vn2_fvtxn);

  //  th1d_os_bbcs_dn2_both_corr->SetMarkerStyle(kOpenCircle);
  th1d_os_bbcs_dn2_both_corr->SetMaximum(th1d_os_bbcs_dn2_both_corr->GetMaximum()*1.2);
  th1d_os_bbcs_dn2_both_corr->Draw();
  tlref->DrawLatex(tlx,tly,Form("#sqrt{c_{%d}{2}} = %.3f",harmonic,os_newref_vn2_bbcs));
  tlmeth->DrawLatex(tlx2,tly2,Form("scalar product with BBCS"));
  tlmeth3->DrawLatex(tlx3,tly3,Form("reference flow using BBCS, FVTXS, FVTXN"));
  c1->Print(Form("FigsClusTrack/newref_os_bbcs_v%d2_%d.png",harmonic,handle));
  c1->Print(Form("FigsClusTrack/newref_os_bbcs_v%d2_%d.pdf",harmonic,handle));

  //  th1d_os_fvtxs_dn2_both_corr->SetMarkerStyle(kOpenCircle);
  th1d_os_fvtxs_dn2_both_corr->SetMaximum(th1d_os_fvtxs_dn2_both_corr->GetMaximum()*1.2);
  th1d_os_fvtxs_dn2_both_corr->Draw();
  tlref->DrawLatex(tlx,tly,Form("#sqrt{c_{%d}{2}} = %.3f",harmonic,os_newref_vn2_fvtxs));
  tlmeth->DrawLatex(tlx2,tly2,Form("scalar product with FVTXS"));
  tlmeth3->DrawLatex(tlx3,tly3,Form("reference flow using BBCS, FVTXS, FVTXN"));
  c1->Print(Form("FigsClusTrack/newref_os_fvtxs_v%d2_%d.png",harmonic,handle));
  c1->Print(Form("FigsClusTrack/newref_os_fvtxs_v%d2_%d.pdf",harmonic,handle));

  //  th1d_os_fvtxn_dn2_both_corr->SetMarkerStyle(kOpenCircle);
  th1d_os_fvtxn_dn2_both_corr->SetMaximum(th1d_os_fvtxn_dn2_both_corr->GetMaximum()*1.2);
  th1d_os_fvtxn_dn2_both_corr->Draw();
  tlref->DrawLatex(tlx,tly,Form("#sqrt{c_{%d}{2}} = %.3f",harmonic,os_newref_vn2_fvtxn));
  tlmeth->DrawLatex(tlx2,tly2,Form("scalar product with FVTXN"));
  tlmeth3->DrawLatex(tlx3,tly3,Form("reference flow using BBCS, FVTXS, FVTXN"));
  c1->Print(Form("FigsClusTrack/newref_os_fvtxn_v%d2_%d.png",harmonic,handle));
  c1->Print(Form("FigsClusTrack/newref_os_fvtxn_v%d2_%d.pdf",harmonic,handle));

  // ---

  // --------------------------------------------------------------------------------------------
  // --- 3 component differential flow scalar product (makes me contemplate 3 particle cumulants)
  // --------------------------------------------------------------------------------------------

  TH1D* th1d_os_bbcsfvtxs_vn2_3csp = (TH1D*)th1d_os_bbcs_dn2_both_corr->Clone();
  TH1D* th1d_os_bbcsfvtxn_vn2_3csp = (TH1D*)th1d_os_bbcs_dn2_both_corr->Clone();
  TH1D* th1d_os_fvtxsfvtxn_vn2_3csp = (TH1D*)th1d_os_bbcs_dn2_both_corr->Clone();

  int nbins = th1d_os_bbcsfvtxs_vn2_3csp->GetNbinsX();
  for ( int i = 0; i < nbins; ++i )
    {
      double dn2_bbcs = th1d_os_bbcs_dn2_both->GetBinContent(i+1);
      double dn2_fvtxs = th1d_os_fvtxs_dn2_both->GetBinContent(i+1);
      double dn2_fvtxn = th1d_os_fvtxn_dn2_both->GetBinContent(i+1);

      double answer_bbcsfvtxs = sqrt((dn2_bbcs*dn2_fvtxs)/os_bbcsfvtxs_cn2_raw);
      double answer_bbcsfvtxn = sqrt((dn2_bbcs*dn2_fvtxn)/os_bbcsfvtxn_cn2_raw);
      double answer_fvtxsfvtxn = sqrt((dn2_fvtxs*dn2_fvtxs)/os_fvtxsfvtxn_cn2_raw);

      if ( answer_bbcsfvtxs != answer_bbcsfvtxs ) answer_bbcsfvtxs = 0;
      if ( answer_bbcsfvtxn != answer_bbcsfvtxn ) answer_bbcsfvtxn = 0;
      if ( answer_fvtxsfvtxn != answer_fvtxsfvtxn ) answer_fvtxsfvtxn = 0;

      th1d_os_bbcsfvtxs_vn2_3csp->SetBinContent(i+1,answer_bbcsfvtxs);
      th1d_os_bbcsfvtxn_vn2_3csp->SetBinContent(i+1,answer_bbcsfvtxn);
      th1d_os_fvtxsfvtxn_vn2_3csp->SetBinContent(i+1,answer_fvtxsfvtxn);
    }

  //  th1d_os_bbcsfvtxs_vn2_3csp->SetMarkerStyle(kOpenCircle);
  th1d_os_bbcsfvtxs_vn2_3csp->SetMaximum(th1d_os_bbcsfvtxs_vn2_3csp->GetBinContent(th1d_os_bbcsfvtxs_vn2_3csp->GetMaximumBin())*1.2);
  th1d_os_bbcsfvtxs_vn2_3csp->Draw();
  tlref->DrawLatex(tlx,tly,Form("#sqrt{c_{%d}{2}} = %.3f",harmonic,sqrt(os_bbcsfvtxs_cn2_raw)));
  tlmeth->DrawLatex(tlx2,tly2,Form("double scalar product with BBCS and FVTXS"));
  tlmeth3->DrawLatex(tlx3,tly3,Form("reference flow using BBCS and FVTXS"));
  c1->Print(Form("FigsClusTrack/threesp_os_bbcsfvtxs_v%d2_%d.png",harmonic,handle));
  c1->Print(Form("FigsClusTrack/threesp_os_bbcsfvtxs_v%d2_%d.pdf",harmonic,handle));

  //  th1d_os_bbcsfvtxn_vn2_3csp->SetMarkerStyle(kOpenCircle);
  th1d_os_bbcsfvtxn_vn2_3csp->SetMaximum(th1d_os_bbcsfvtxn_vn2_3csp->GetBinContent(th1d_os_bbcsfvtxn_vn2_3csp->GetMaximumBin())*1.2);
  th1d_os_bbcsfvtxn_vn2_3csp->Draw();
  tlref->DrawLatex(tlx,tly,Form("#sqrt{c_{%d}{2}} = %.3f",harmonic,sqrt(os_bbcsfvtxn_cn2_raw)));
  tlmeth->DrawLatex(tlx2,tly2,Form("double scalar product with BBCS and FVTXN"));
  tlmeth3->DrawLatex(tlx3,tly3,Form("reference flow using BBCS and FVTXN"));
  c1->Print(Form("FigsClusTrack/threesp_os_bbcsfvtxn_v%d2_%d.png",harmonic,handle));
  c1->Print(Form("FigsClusTrack/threesp_os_bbcsfvtxn_v%d2_%d.pdf",harmonic,handle));

  //  th1d_os_fvtxsfvtxn_vn2_3csp->SetMarkerStyle(kOpenCircle);
  th1d_os_fvtxsfvtxn_vn2_3csp->SetMaximum(th1d_os_fvtxsfvtxn_vn2_3csp->GetBinContent(th1d_os_fvtxsfvtxn_vn2_3csp->GetMaximumBin())*1.2);
  th1d_os_fvtxsfvtxn_vn2_3csp->Draw();
  tlref->DrawLatex(tlx,tly,Form("#sqrt{c_{%d}{2}} = %.3f",harmonic,sqrt(os_fvtxsfvtxn_cn2_raw)));
  tlmeth->DrawLatex(tlx2,tly2,Form("double scalar product with FVTXS and FVTXN"));
  tlmeth3->DrawLatex(tlx3,tly3,Form("reference flow using FVTXS and FVTXN"));
  c1->Print(Form("FigsClusTrack/threesp_os_fvtxsfvtxn_v%d2_%d.png",harmonic,handle));
  c1->Print(Form("FigsClusTrack/threesp_os_fvtxsfvtxn_v%d2_%d.pdf",harmonic,handle));

  // ---

  // ----------------------------
  // --- a few simple comparisons
  // ----------------------------

  th1d_os_bbcsfvtxs_vn2_3csp->SetMinimum(0);
  th1d_os_bbcsfvtxs_vn2_3csp->SetMaximum(0.2);
  if ( handle <= 39 ) th1d_os_bbcsfvtxs_vn2_3csp->SetMaximum(0.35);
  th1d_os_bbcsfvtxs_vn2_3csp->Draw();
  th1d_os_bbcsfvtxn_vn2_3csp->Draw("same");
  th1d_os_fvtxsfvtxn_vn2_3csp->Draw("same");
  tge_ep->Draw("p");
  th1d_os_bbcsfvtxs_vn2_3csp->SetLineColor(kBlack);
  th1d_os_bbcsfvtxn_vn2_3csp->SetLineColor(kRed);
  th1d_os_fvtxsfvtxn_vn2_3csp->SetLineColor(kBlue);
  TLegend* leg = new TLegend(0.18,0.68,0.38,0.88);
  leg->SetHeader("double scalar product");
  leg->AddEntry(th1d_os_bbcsfvtxs_vn2_3csp,"BBCS and FVTXS","el");
  leg->AddEntry(th1d_os_bbcsfvtxn_vn2_3csp,"BBCS and FVTXN","el");
  leg->AddEntry(th1d_os_fvtxsfvtxn_vn2_3csp,"FVTXS and FVTXN","el");
  leg->AddEntry(tge_ep,"FVTXS event plane","ep");
  leg->SetTextSize(0.04);
  leg->SetFillStyle(0);
  leg->Draw();
  c1->Print(Form("FigsClusTrack/threesp_os_compare_v%d2_%d.png",harmonic,handle));
  c1->Print(Form("FigsClusTrack/threesp_os_compare_v%d2_%d.pdf",harmonic,handle));

  th1d_os_bbcs_dn2_both_corr->SetMinimum(0);
  th1d_os_bbcs_dn2_both_corr->SetMaximum(0.2);
  if ( handle <= 39 ) th1d_os_bbcs_dn2_both_corr->SetMaximum(0.35);
  th1d_os_bbcs_dn2_both_corr->Draw();
  th1d_os_fvtxs_dn2_both_corr->Draw("same");
  th1d_os_fvtxn_dn2_both_corr->Draw("same");
  tge_ep->Draw("p");
  th1d_os_bbcs_dn2_both_corr->SetLineColor(kBlack);
  th1d_os_fvtxs_dn2_both_corr->SetLineColor(kRed);
  th1d_os_fvtxn_dn2_both_corr->SetLineColor(kBlue);
  if ( leg ) delete leg;
  leg = new TLegend(0.18,0.68,0.38,0.88);
  leg->SetHeader("scalar product");
  leg->AddEntry(th1d_os_bbcs_dn2_both_corr,"BBCS","el");
  leg->AddEntry(th1d_os_fvtxs_dn2_both_corr,"FVTXS","el");
  leg->AddEntry(th1d_os_fvtxn_dn2_both_corr,"FVTXN","el");
  leg->AddEntry(tge_ep,"FVTXS event plane","ep");
  leg->SetTextSize(0.04);
  leg->SetFillStyle(0);
  leg->Draw();
  c1->Print(Form("FigsClusTrack/newref_os_compare_v%d2_%d.png",harmonic,handle));
  c1->Print(Form("FigsClusTrack/newref_os_compare_v%d2_%d.pdf",harmonic,handle));

  th1d_os_bbcs_dn2_both_corr->SetLineColor(kBlue);
  th1d_os_fvtxs_dn2_both_corr->SetLineColor(kRed);
  th1d_os_bbcsfvtxs_vn2_3csp->SetLineColor(kBlack);
  th1d_os_bbcsfvtxs_vn2_3csp->Draw();
  th1d_os_bbcs_dn2_both_corr->Draw("same");
  th1d_os_fvtxs_dn2_both_corr->Draw("same");
  tge_ep->Draw("p");
  if ( leg ) delete leg;
  leg = new TLegend(0.18,0.68,0.38,0.88);
  //leg->SetHeader("scalar product");
  leg->AddEntry(tge_ep,"FVTXS event plane","ep");
  leg->AddEntry(th1d_os_bbcsfvtxs_vn2_3csp,"BBCS and FVTXS","el");
  leg->AddEntry(th1d_os_bbcs_dn2_both_corr,"BBCS","el");
  leg->AddEntry(th1d_os_fvtxs_dn2_both_corr,"FVTXS","el");
  leg->SetTextSize(0.04);
  leg->SetFillStyle(0);
  leg->Draw();
  c1->Print(Form("FigsClusTrack/nr3c_os_compare_v%d2_%d.png",harmonic,handle));
  c1->Print(Form("FigsClusTrack/nr3c_os_compare_v%d2_%d.pdf",harmonic,handle));

  // ---

  th1d_os_bbcsfvtxs_vn2_3csp->Draw();
  th1d_os_bbcsfvtxs_vn2_3csp->GetYaxis()->SetTitle(Form("v_{%d}",harmonic));
  th1d_os_bbcsfvtxs_vn2_3csp->GetYaxis()->SetTitleOffset(1.2);
  if ( harmonic == 3 )
    {
      th1d_os_bbcsfvtxs_vn2_3csp->SetMaximum(0.08);
      th1d_os_bbcsfvtxs_vn2_3csp->SetMinimum(-0.02);
    }
  tge_ep->Draw("p");
  if ( leg ) delete leg;
  leg = new TLegend(0.18,0.68,0.38,0.88);
  leg->AddEntry(tge_ep,"event plane (FVTXS)","ep");
  leg->AddEntry(th1d_os_bbcsfvtxs_vn2_3csp,"scalar product (3sub)","el");
  leg->SetTextSize(0.04);
  leg->SetFillStyle(0);
  leg->Draw();
  c1->Print(Form("FigsClusTrack/run16dau_2pcEPcompare_v%d2_%d.png",harmonic,handle));
  c1->Print(Form("FigsClusTrack/run16dau_2pcEPcompare_v%d2_%d.pdf",harmonic,handle));

  delete c1;


}
