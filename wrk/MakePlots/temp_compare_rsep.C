void doit(int);

void temp_compare_rsep()
{

  doit(456652);
  // doit(200);
  // doit(62);
  // doit(39);
  // doit(20);

}

void doit(int handle)
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

  // ---
  // ---
  // ---

  // --- BBCS

  TProfile* tp1f_os_bbcs_v2_both = (TProfile*)file->Get("os_bbcs_v2_both");
  TProfile* hv2_bbcs = (TProfile*)file->Get("bbcs_v2_both_docalib");

  if ( handle <= 200 ) tp1f_os_bbcs_v2_both->SetTitle(Form("d+Au #sqrt{s_{NN}} = %d GeV",handle));
  tp1f_os_bbcs_v2_both->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  tp1f_os_bbcs_v2_both->GetYaxis()->SetTitle("uncorrected v_{2}{EP}");
  tp1f_os_bbcs_v2_both->GetYaxis()->SetTitleOffset(1.4);

  if ( handle <= 200 ) hv2_bbcs->SetTitle(Form("d+Au #sqrt{s_{NN}} = %d GeV",handle));
  hv2_bbcs->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hv2_bbcs->GetYaxis()->SetTitle("uncorrected v_{2}{EP}");
  hv2_bbcs->GetYaxis()->SetTitleOffset(1.4);

  hv2_bbcs->SetLineColor(kBlack);
  hv2_bbcs->Draw();
  c1->Print(Form("FigsEventPlane/comp_s_bbcs_v2raw_%d.png",handle));
  c1->Print(Form("FigsEventPlane/comp_s_bbcs_v2raw_%d.pdf",handle));

  tp1f_os_bbcs_v2_both->SetLineColor(kBlue);
  tp1f_os_bbcs_v2_both->Draw();
  c1->Print(Form("FigsEventPlane/comp_r_bbcs_v2raw_%d.png",handle));
  c1->Print(Form("FigsEventPlane/comp_r_bbcs_v2raw_%d.pdf",handle));

  hv2_bbcs->Draw("same");
  TLegend* leg = new TLegend(0.16,0.78,0.38,0.88);
  leg->AddEntry(hv2_bbcs,"Fourier decomposition","el");
  leg->AddEntry(tp1f_os_bbcs_v2_both,"Q-vector recentering","el");
  leg->SetTextSize(0.045);
  leg->SetFillStyle(0);
  leg->Draw();
  c1->Print(Form("FigsEventPlane/comp_rs_bbcs_v2raw_%d.png",handle));
  c1->Print(Form("FigsEventPlane/comp_rs_bbcs_v2raw_%d.pdf",handle));

  // --- FVTXS

  TProfile* tp1f_os_fvtxs_v2_both = (TProfile*)file->Get("os_fvtxs_v2_both");
  TProfile* hv2_fvtxs = (TProfile*)file->Get("fvtxs_v2_both_docalib");


  if ( handle <= 200 ) tp1f_os_fvtxs_v2_both->SetTitle(Form("d+Au #sqrt{s_{NN}} = %d GeV",handle));
  tp1f_os_fvtxs_v2_both->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  tp1f_os_fvtxs_v2_both->GetYaxis()->SetTitle("uncorrected v_{2}{EP}");
  tp1f_os_fvtxs_v2_both->GetYaxis()->SetTitleOffset(1.4);

  if ( handle <= 200 ) hv2_fvtxs->SetTitle(Form("d+Au #sqrt{s_{NN}} = %d GeV",handle));
  hv2_fvtxs->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hv2_fvtxs->GetYaxis()->SetTitle("uncorrected v_{2}{EP}");
  hv2_fvtxs->GetYaxis()->SetTitleOffset(1.4);

  hv2_fvtxs->SetLineColor(kBlack);
  hv2_fvtxs->Draw();
  c1->Print(Form("FigsEventPlane/comp_s_fvtxs_v2raw_%d.png",handle));
  c1->Print(Form("FigsEventPlane/comp_s_fvtxs_v2raw_%d.pdf",handle));

  tp1f_os_fvtxs_v2_both->SetLineColor(kBlue);
  tp1f_os_fvtxs_v2_both->Draw();
  c1->Print(Form("FigsEventPlane/comp_r_fvtxs_v2raw_%d.png",handle));
  c1->Print(Form("FigsEventPlane/comp_r_fvtxs_v2raw_%d.pdf",handle));

  hv2_fvtxs->Draw("same");
  leg->Draw();
  c1->Print(Form("FigsEventPlane/comp_rs_fvtxs_v2raw_%d.png",handle));
  c1->Print(Form("FigsEventPlane/comp_rs_fvtxs_v2raw_%d.pdf",handle));

  // --- FVTXN

  TProfile* tp1f_os_fvtxn_v2_both = (TProfile*)file->Get("os_fvtxn_v2_both");
  TProfile* hv2_fvtxn = (TProfile*)file->Get("fvtxn_v2_both_docalib");


  if ( handle <= 200 ) tp1f_os_fvtxn_v2_both->SetTitle(Form("d+Au #sqrt{s_{NN}} = %d GeV",handle));
  tp1f_os_fvtxn_v2_both->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  tp1f_os_fvtxn_v2_both->GetYaxis()->SetTitle("uncorrected v_{2}{EP}");
  tp1f_os_fvtxn_v2_both->GetYaxis()->SetTitleOffset(1.4);

  if ( handle <= 200 ) hv2_fvtxn->SetTitle(Form("d+Au #sqrt{s_{NN}} = %d GeV",handle));
  hv2_fvtxn->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hv2_fvtxn->GetYaxis()->SetTitle("uncorrected v_{2}{EP}");
  hv2_fvtxn->GetYaxis()->SetTitleOffset(1.4);

  hv2_fvtxn->SetLineColor(kBlack);
  hv2_fvtxn->Draw();
  c1->Print(Form("FigsEventPlane/comp_s_fvtxn_v2raw_%d.png",handle));
  c1->Print(Form("FigsEventPlane/comp_s_fvtxn_v2raw_%d.pdf",handle));

  tp1f_os_fvtxn_v2_both->SetLineColor(kBlue);
  tp1f_os_fvtxn_v2_both->Draw();
  c1->Print(Form("FigsEventPlane/comp_r_fvtxn_v2raw_%d.png",handle));
  c1->Print(Form("FigsEventPlane/comp_r_fvtxn_v2raw_%d.pdf",handle));

  hv2_fvtxn->Draw("same");
  leg->Draw();
  c1->Print(Form("FigsEventPlane/comp_rs_fvtxn_v2raw_%d.png",handle));
  c1->Print(Form("FigsEventPlane/comp_rs_fvtxn_v2raw_%d.pdf",handle));

  // ---

  delete leg;
  delete c1;

}
