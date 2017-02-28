TH1D* dueprocess(TProfile*,TProfile*,TProfile*,TProfile*,TProfile*);

TH1D* dueprocess(TH1D*,TH1D*,TH1D*,TH1D*,TH1D*);

void diehanddieverletzt()
{

  TCanvas* c1 = new TCanvas();
  c1->SetGrid();

  TFile* file = TFile::Open("input/combined_200.root");

  TProfile* tp1f_base;
  TProfile* tp1f_cosphi;
  TProfile* tp1f_sinphi;
  TProfile* tp1f_cospsi;
  TProfile* tp1f_sinpsi;

  // --- all layers
  tp1f_base = (TProfile*)file->Get("bbcs_v3_west_docalib");
  tp1f_cosphi = (TProfile*)file->Get("bbcs_v3_west_cosphi");
  tp1f_sinphi = (TProfile*)file->Get("bbcs_v3_west_sinphi");
  tp1f_cospsi = (TProfile*)file->Get("bbcs_v3_west_cospsi3");
  tp1f_sinpsi = (TProfile*)file->Get("bbcs_v3_west_sinpsi3");
  TH1D* th1d_sB_west = dueprocess(tp1f_base,tp1f_cosphi,tp1f_sinphi,tp1f_cospsi,tp1f_sinpsi);
  tp1f_base = (TProfile*)file->Get("bbcs_v3_east_docalib");
  tp1f_cosphi = (TProfile*)file->Get("bbcs_v3_east_cosphi");
  tp1f_sinphi = (TProfile*)file->Get("bbcs_v3_east_sinphi");
  tp1f_cospsi = (TProfile*)file->Get("bbcs_v3_east_cospsi3");
  tp1f_sinpsi = (TProfile*)file->Get("bbcs_v3_east_sinpsi3");
  TH1D* th1d_sB_east = dueprocess(tp1f_base,tp1f_cosphi,tp1f_sinphi,tp1f_cospsi,tp1f_sinpsi);
  tp1f_base = (TProfile*)file->Get("bbcs_v3_both_docalib");
  tp1f_cosphi = (TProfile*)file->Get("bbcs_v3_both_cosphi");
  tp1f_sinphi = (TProfile*)file->Get("bbcs_v3_both_sinphi");
  tp1f_cospsi = (TProfile*)file->Get("bbcs_v3_both_cospsi3");
  tp1f_sinpsi = (TProfile*)file->Get("bbcs_v3_both_sinpsi3");
  TH1D* th1d_sB_both = dueprocess(tp1f_base,tp1f_cosphi,tp1f_sinphi,tp1f_cospsi,tp1f_sinpsi);
  // --- all layers
  tp1f_base = (TProfile*)file->Get("fvtxs_v3_west_docalib");
  tp1f_cosphi = (TProfile*)file->Get("fvtxs_v3_west_cosphi");
  tp1f_sinphi = (TProfile*)file->Get("fvtxs_v3_west_sinphi");
  tp1f_cospsi = (TProfile*)file->Get("fvtxs_v3_west_cospsi3");
  tp1f_sinpsi = (TProfile*)file->Get("fvtxs_v3_west_sinpsi3");
  TH1D* th1d_sA_west = dueprocess(tp1f_base,tp1f_cosphi,tp1f_sinphi,tp1f_cospsi,tp1f_sinpsi);
  tp1f_base = (TProfile*)file->Get("fvtxs_v3_east_docalib");
  tp1f_cosphi = (TProfile*)file->Get("fvtxs_v3_east_cosphi");
  tp1f_sinphi = (TProfile*)file->Get("fvtxs_v3_east_sinphi");
  tp1f_cospsi = (TProfile*)file->Get("fvtxs_v3_east_cospsi3");
  tp1f_sinpsi = (TProfile*)file->Get("fvtxs_v3_east_sinpsi3");
  TH1D* th1d_sA_east = dueprocess(tp1f_base,tp1f_cosphi,tp1f_sinphi,tp1f_cospsi,tp1f_sinpsi);
  tp1f_base = (TProfile*)file->Get("fvtxs_v3_both_docalib");
  tp1f_cosphi = (TProfile*)file->Get("fvtxs_v3_both_cosphi");
  tp1f_sinphi = (TProfile*)file->Get("fvtxs_v3_both_sinphi");
  tp1f_cospsi = (TProfile*)file->Get("fvtxs_v3_both_cospsi3");
  tp1f_sinpsi = (TProfile*)file->Get("fvtxs_v3_both_sinpsi3");
  TH1D* th1d_sA_both = dueprocess(tp1f_base,tp1f_cosphi,tp1f_sinphi,tp1f_cospsi,tp1f_sinpsi);
  // --- layer s0
  tp1f_base = (TProfile*)file->Get("fvtxs0_v3_west_docalib");
  tp1f_cosphi = (TProfile*)file->Get("fvtxs0_v3_west_cosphi");
  tp1f_sinphi = (TProfile*)file->Get("fvtxs0_v3_west_sinphi");
  tp1f_cospsi = (TProfile*)file->Get("fvtxs0_v3_west_cospsi3");
  tp1f_sinpsi = (TProfile*)file->Get("fvtxs0_v3_west_sinpsi3");
  TH1D* th1d_s0_west = dueprocess(tp1f_base,tp1f_cosphi,tp1f_sinphi,tp1f_cospsi,tp1f_sinpsi);
  tp1f_base = (TProfile*)file->Get("fvtxs0_v3_east_docalib");
  tp1f_cosphi = (TProfile*)file->Get("fvtxs0_v3_east_cosphi");
  tp1f_sinphi = (TProfile*)file->Get("fvtxs0_v3_east_sinphi");
  tp1f_cospsi = (TProfile*)file->Get("fvtxs0_v3_east_cospsi3");
  tp1f_sinpsi = (TProfile*)file->Get("fvtxs0_v3_east_sinpsi3");
  TH1D* th1d_s0_east = dueprocess(tp1f_base,tp1f_cosphi,tp1f_sinphi,tp1f_cospsi,tp1f_sinpsi);
  tp1f_base = (TProfile*)file->Get("fvtxs0_v3_both_docalib");
  tp1f_cosphi = (TProfile*)file->Get("fvtxs0_v3_both_cosphi");
  tp1f_sinphi = (TProfile*)file->Get("fvtxs0_v3_both_sinphi");
  tp1f_cospsi = (TProfile*)file->Get("fvtxs0_v3_both_cospsi3");
  tp1f_sinpsi = (TProfile*)file->Get("fvtxs0_v3_both_sinpsi3");
  TH1D* th1d_s0_both = dueprocess(tp1f_base,tp1f_cosphi,tp1f_sinphi,tp1f_cospsi,tp1f_sinpsi);
  // --- layer s1
  tp1f_base = (TProfile*)file->Get("fvtxs1_v3_west_docalib");
  tp1f_cosphi = (TProfile*)file->Get("fvtxs1_v3_west_cosphi");
  tp1f_sinphi = (TProfile*)file->Get("fvtxs1_v3_west_sinphi");
  tp1f_cospsi = (TProfile*)file->Get("fvtxs1_v3_west_cospsi3");
  tp1f_sinpsi = (TProfile*)file->Get("fvtxs1_v3_west_sinpsi3");
  TH1D* th1d_s1_west = dueprocess(tp1f_base,tp1f_cosphi,tp1f_sinphi,tp1f_cospsi,tp1f_sinpsi);
  tp1f_base = (TProfile*)file->Get("fvtxs1_v3_east_docalib");
  tp1f_cosphi = (TProfile*)file->Get("fvtxs1_v3_east_cosphi");
  tp1f_sinphi = (TProfile*)file->Get("fvtxs1_v3_east_sinphi");
  tp1f_cospsi = (TProfile*)file->Get("fvtxs1_v3_east_cospsi3");
  tp1f_sinpsi = (TProfile*)file->Get("fvtxs1_v3_east_sinpsi3");
  TH1D* th1d_s1_east = dueprocess(tp1f_base,tp1f_cosphi,tp1f_sinphi,tp1f_cospsi,tp1f_sinpsi);
  tp1f_base = (TProfile*)file->Get("fvtxs1_v3_both_docalib");
  tp1f_cosphi = (TProfile*)file->Get("fvtxs1_v3_both_cosphi");
  tp1f_sinphi = (TProfile*)file->Get("fvtxs1_v3_both_sinphi");
  tp1f_cospsi = (TProfile*)file->Get("fvtxs1_v3_both_cospsi3");
  tp1f_sinpsi = (TProfile*)file->Get("fvtxs1_v3_both_sinpsi3");
  TH1D* th1d_s1_both = dueprocess(tp1f_base,tp1f_cosphi,tp1f_sinphi,tp1f_cospsi,tp1f_sinpsi);
  // --- layer s2
  tp1f_base = (TProfile*)file->Get("fvtxs2_v3_west_docalib");
  tp1f_cosphi = (TProfile*)file->Get("fvtxs2_v3_west_cosphi");
  tp1f_sinphi = (TProfile*)file->Get("fvtxs2_v3_west_sinphi");
  tp1f_cospsi = (TProfile*)file->Get("fvtxs2_v3_west_cospsi3");
  tp1f_sinpsi = (TProfile*)file->Get("fvtxs2_v3_west_sinpsi3");
  TH1D* th1d_s2_west = dueprocess(tp1f_base,tp1f_cosphi,tp1f_sinphi,tp1f_cospsi,tp1f_sinpsi);
  tp1f_base = (TProfile*)file->Get("fvtxs2_v3_east_docalib");
  tp1f_cosphi = (TProfile*)file->Get("fvtxs2_v3_east_cosphi");
  tp1f_sinphi = (TProfile*)file->Get("fvtxs2_v3_east_sinphi");
  tp1f_cospsi = (TProfile*)file->Get("fvtxs2_v3_east_cospsi3");
  tp1f_sinpsi = (TProfile*)file->Get("fvtxs2_v3_east_sinpsi3");
  TH1D* th1d_s2_east = dueprocess(tp1f_base,tp1f_cosphi,tp1f_sinphi,tp1f_cospsi,tp1f_sinpsi);
  tp1f_base = (TProfile*)file->Get("fvtxs2_v3_both_docalib");
  tp1f_cosphi = (TProfile*)file->Get("fvtxs2_v3_both_cosphi");
  tp1f_sinphi = (TProfile*)file->Get("fvtxs2_v3_both_sinphi");
  tp1f_cospsi = (TProfile*)file->Get("fvtxs2_v3_both_cospsi3");
  tp1f_sinpsi = (TProfile*)file->Get("fvtxs2_v3_both_sinpsi3");
  TH1D* th1d_s2_both = dueprocess(tp1f_base,tp1f_cosphi,tp1f_sinphi,tp1f_cospsi,tp1f_sinpsi);
  // --- layer s3
  tp1f_base = (TProfile*)file->Get("fvtxs3_v3_west_docalib");
  tp1f_cosphi = (TProfile*)file->Get("fvtxs3_v3_west_cosphi");
  tp1f_sinphi = (TProfile*)file->Get("fvtxs3_v3_west_sinphi");
  tp1f_cospsi = (TProfile*)file->Get("fvtxs3_v3_west_cospsi3");
  tp1f_sinpsi = (TProfile*)file->Get("fvtxs3_v3_west_sinpsi3");
  TH1D* th1d_s3_west = dueprocess(tp1f_base,tp1f_cosphi,tp1f_sinphi,tp1f_cospsi,tp1f_sinpsi);
  tp1f_base = (TProfile*)file->Get("fvtxs3_v3_east_docalib");
  tp1f_cosphi = (TProfile*)file->Get("fvtxs3_v3_east_cosphi");
  tp1f_sinphi = (TProfile*)file->Get("fvtxs3_v3_east_sinphi");
  tp1f_cospsi = (TProfile*)file->Get("fvtxs3_v3_east_cospsi3");
  tp1f_sinpsi = (TProfile*)file->Get("fvtxs3_v3_east_sinpsi3");
  TH1D* th1d_s3_east = dueprocess(tp1f_base,tp1f_cosphi,tp1f_sinphi,tp1f_cospsi,tp1f_sinpsi);
  tp1f_base = (TProfile*)file->Get("fvtxs3_v3_both_docalib");
  tp1f_cosphi = (TProfile*)file->Get("fvtxs3_v3_both_cosphi");
  tp1f_sinphi = (TProfile*)file->Get("fvtxs3_v3_both_sinphi");
  tp1f_cospsi = (TProfile*)file->Get("fvtxs3_v3_both_cospsi3");
  tp1f_sinpsi = (TProfile*)file->Get("fvtxs3_v3_both_sinpsi3");
  TH1D* th1d_s3_both = dueprocess(tp1f_base,tp1f_cosphi,tp1f_sinphi,tp1f_cospsi,tp1f_sinpsi);



  TLegend* leg = NULL;

  th1d_sB_both->SetLineColor(kBlack);
  th1d_sB_west->SetLineColor(kRed);
  th1d_sB_east->SetLineColor(kBlue);
  th1d_sB_both->Draw();
  th1d_sB_both->SetMaximum(0.005);
  th1d_sB_both->SetMinimum(-0.005);
  th1d_sB_both->GetYaxis()->SetTitle("v_{3} (no EP corr)");
  th1d_sB_both->GetYaxis()->SetTitleOffset(1.5);
  th1d_sB_both->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  th1d_sB_west->Draw("same");
  th1d_sB_east->Draw("same");
  leg = new TLegend(0.18,0.68,0.38,0.88);
  leg->AddEntry(th1d_sB_both,"both arms","el");
  leg->AddEntry(th1d_sB_east,"east arm","el");
  leg->AddEntry(th1d_sB_west,"west arm","el");
  leg->SetTextSize(0.045);
  leg->Draw();
  c1->Print("FigsCheckThree/fig_bbcs_allarms.png");
  c1->Print("FigsCheckThree/fig_bbcs_allarms.pdf");

  th1d_sA_both->SetLineColor(kBlack);
  th1d_sA_west->SetLineColor(kRed);
  th1d_sA_east->SetLineColor(kBlue);
  th1d_sA_both->Draw();
  th1d_sA_both->SetMaximum(0.005);
  th1d_sA_both->SetMinimum(-0.005);
  th1d_sA_both->GetYaxis()->SetTitle("v_{3} (no EP corr)");
  th1d_sA_both->GetYaxis()->SetTitleOffset(1.5);
  th1d_sA_both->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  th1d_sA_west->Draw("same");
  th1d_sA_east->Draw("same");
  leg = new TLegend(0.18,0.68,0.38,0.88);
  leg->AddEntry(th1d_sA_both,"both arms","el");
  leg->AddEntry(th1d_sA_east,"east arm","el");
  leg->AddEntry(th1d_sA_west,"west arm","el");
  leg->SetTextSize(0.045);
  leg->Draw();
  c1->Print("FigsCheckThree/fig_layer_sA_allarms.png");
  c1->Print("FigsCheckThree/fig_layer_sA_allarms.pdf");

  th1d_s0_both->SetLineColor(kBlack);
  th1d_s0_west->SetLineColor(kRed);
  th1d_s0_east->SetLineColor(kBlue);
  th1d_s0_both->Draw();
  th1d_s0_both->SetMaximum(0.005);
  th1d_s0_both->SetMinimum(-0.005);
  th1d_s0_both->GetYaxis()->SetTitle("v_{3} (no EP corr)");
  th1d_s0_both->GetYaxis()->SetTitleOffset(1.5);
  th1d_s0_both->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  th1d_s0_west->Draw("same");
  th1d_s0_east->Draw("same");
  leg->Draw();
  c1->Print("FigsCheckThree/fig_layer_s0_allarms.png");
  c1->Print("FigsCheckThree/fig_layer_s0_allarms.pdf");

  th1d_s1_both->SetLineColor(kBlack);
  th1d_s1_west->SetLineColor(kRed);
  th1d_s1_east->SetLineColor(kBlue);
  th1d_s1_both->Draw();
  th1d_s1_both->SetMaximum(0.005);
  th1d_s1_both->SetMinimum(-0.005);
  th1d_s1_both->GetYaxis()->SetTitle("v_{3} (no EP corr)");
  th1d_s1_both->GetYaxis()->SetTitleOffset(1.5);
  th1d_s1_both->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  th1d_s1_west->Draw("same");
  th1d_s1_east->Draw("same");
  leg->Draw();
  c1->Print("FigsCheckThree/fig_layer_s1_allarms.png");
  c1->Print("FigsCheckThree/fig_layer_s1_allarms.pdf");

  th1d_s2_both->SetLineColor(kBlack);
  th1d_s2_west->SetLineColor(kRed);
  th1d_s2_east->SetLineColor(kBlue);
  th1d_s2_both->Draw();
  th1d_s2_both->SetMaximum(0.005);
  th1d_s2_both->SetMinimum(-0.005);
  th1d_s2_both->GetYaxis()->SetTitle("v_{3} (no EP corr)");
  th1d_s2_both->GetYaxis()->SetTitleOffset(1.5);
  th1d_s2_both->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  th1d_s2_west->Draw("same");
  th1d_s2_east->Draw("same");
  leg->Draw();
  c1->Print("FigsCheckThree/fig_layer_s2_allarms.png");
  c1->Print("FigsCheckThree/fig_layer_s2_allarms.pdf");

  th1d_s3_both->SetLineColor(kBlack);
  th1d_s3_west->SetLineColor(kRed);
  th1d_s3_east->SetLineColor(kBlue);
  th1d_s3_both->Draw();
  th1d_s3_both->SetMaximum(0.005);
  th1d_s3_both->SetMinimum(-0.005);
  th1d_s3_both->GetYaxis()->SetTitle("v_{3} (no EP corr)");
  th1d_s3_both->GetYaxis()->SetTitleOffset(1.5);
  th1d_s3_both->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  th1d_s3_west->Draw("same");
  th1d_s3_east->Draw("same");
  leg->Draw();
  c1->Print("FigsCheckThree/fig_layer_s3_allarms.png");
  c1->Print("FigsCheckThree/fig_layer_s3_allarms.pdf");

  // ---

  TH1D* th1d_nc_sB_both = ( (TProfile*) file->Get("bbcs_v3_both_docalib") )->ProjectionX("bbcs_v3_both_docalib_px2nc");
  TH1D* th1d_nc_sB_west = ( (TProfile*) file->Get("bbcs_v3_west_docalib") )->ProjectionX("bbcs_v3_west_docalib_px2nc");
  TH1D* th1d_nc_sB_east = ( (TProfile*) file->Get("bbcs_v3_east_docalib") )->ProjectionX("bbcs_v3_east_docalib_px2nc");
  th1d_nc_sB_both->SetLineColor(kBlack);
  th1d_nc_sB_west->SetLineColor(kRed);
  th1d_nc_sB_east->SetLineColor(kBlue);
  th1d_nc_sB_both->Draw();
  th1d_nc_sB_both->SetMaximum(0.005);
  th1d_nc_sB_both->SetMinimum(-0.005);
  th1d_nc_sB_both->GetYaxis()->SetTitle("v_{3} (no EP corr)");
  th1d_nc_sB_both->GetYaxis()->SetTitleOffset(1.5);
  th1d_nc_sB_both->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  th1d_nc_sB_west->Draw("same");
  th1d_nc_sB_east->Draw("same");
  leg->Draw();
  c1->Print("FigsCheckThree/fig_nocorrterms_bbcs_allarms.png");
  c1->Print("FigsCheckThree/fig_nocorrterms_bbcs_allarms.pdf");

  TH1D* th1d_nc_sA_both = ( (TProfile*) file->Get("fvtxs_v3_both_docalib") )->ProjectionX("fvtxs_v3_both_docalib_px2nc");
  TH1D* th1d_nc_sA_west = ( (TProfile*) file->Get("fvtxs_v3_west_docalib") )->ProjectionX("fvtxs_v3_west_docalib_px2nc");
  TH1D* th1d_nc_sA_east = ( (TProfile*) file->Get("fvtxs_v3_east_docalib") )->ProjectionX("fvtxs_v3_east_docalib_px2nc");
  th1d_nc_sA_both->SetLineColor(kBlack);
  th1d_nc_sA_west->SetLineColor(kRed);
  th1d_nc_sA_east->SetLineColor(kBlue);
  th1d_nc_sA_both->Draw();
  th1d_nc_sA_both->SetMaximum(0.005);
  th1d_nc_sA_both->SetMinimum(-0.005);
  th1d_nc_sA_both->GetYaxis()->SetTitle("v_{3} (no EP corr)");
  th1d_nc_sA_both->GetYaxis()->SetTitleOffset(1.5);
  th1d_nc_sA_both->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  th1d_nc_sA_west->Draw("same");
  th1d_nc_sA_east->Draw("same");
  leg->Draw();
  c1->Print("FigsCheckThree/fig_nocorrterms_layer_sA_allarms.png");
  c1->Print("FigsCheckThree/fig_nocorrterms_layer_sA_allarms.pdf");

  TH1D* th1d_nc_s0_both = ( (TProfile*) file->Get("fvtxs0_v3_both_docalib") )->ProjectionX("fvtxs0_v3_both_docalib_px2nc");
  TH1D* th1d_nc_s0_west = ( (TProfile*) file->Get("fvtxs0_v3_west_docalib") )->ProjectionX("fvtxs0_v3_west_docalib_px2nc");
  TH1D* th1d_nc_s0_east = ( (TProfile*) file->Get("fvtxs0_v3_east_docalib") )->ProjectionX("fvtxs0_v3_east_docalib_px2nc");
  th1d_nc_s0_both->SetLineColor(kBlack);
  th1d_nc_s0_west->SetLineColor(kRed);
  th1d_nc_s0_east->SetLineColor(kBlue);
  th1d_nc_s0_both->Draw();
  th1d_nc_s0_both->SetMaximum(0.005);
  th1d_nc_s0_both->SetMinimum(-0.005);
  th1d_nc_s0_both->GetYaxis()->SetTitle("v_{3} (no EP corr)");
  th1d_nc_s0_both->GetYaxis()->SetTitleOffset(1.5);
  th1d_nc_s0_both->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  th1d_nc_s0_west->Draw("same");
  th1d_nc_s0_east->Draw("same");
  leg->Draw();
  c1->Print("FigsCheckThree/fig_nocorrterms_layer_s0_allarms.png");
  c1->Print("FigsCheckThree/fig_nocorrterms_layer_s0_allarms.pdf");

  TH1D* th1d_nc_s1_both = ( (TProfile*) file->Get("fvtxs1_v3_both_docalib") )->ProjectionX("fvtxs1_v3_both_docalib_px2nc");
  TH1D* th1d_nc_s1_west = ( (TProfile*) file->Get("fvtxs1_v3_west_docalib") )->ProjectionX("fvtxs1_v3_west_docalib_px2nc");
  TH1D* th1d_nc_s1_east = ( (TProfile*) file->Get("fvtxs1_v3_east_docalib") )->ProjectionX("fvtxs1_v3_east_docalib_px2nc");
  th1d_nc_s1_both->SetLineColor(kBlack);
  th1d_nc_s1_west->SetLineColor(kRed);
  th1d_nc_s1_east->SetLineColor(kBlue);
  th1d_nc_s1_both->Draw();
  th1d_nc_s1_both->SetMaximum(0.005);
  th1d_nc_s1_both->SetMinimum(-0.005);
  th1d_nc_s1_both->GetYaxis()->SetTitle("v_{3} (no EP corr)");
  th1d_nc_s1_both->GetYaxis()->SetTitleOffset(1.5);
  th1d_nc_s1_both->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  th1d_nc_s1_west->Draw("same");
  th1d_nc_s1_east->Draw("same");
  leg->Draw();
  c1->Print("FigsCheckThree/fig_nocorrterms_layer_s1_allarms.png");
  c1->Print("FigsCheckThree/fig_nocorrterms_layer_s1_allarms.pdf");

  TH1D* th1d_nc_s2_both = ( (TProfile*) file->Get("fvtxs2_v3_both_docalib") )->ProjectionX("fvtxs2_v3_both_docalib_px2nc");
  TH1D* th1d_nc_s2_west = ( (TProfile*) file->Get("fvtxs2_v3_west_docalib") )->ProjectionX("fvtxs2_v3_west_docalib_px2nc");
  TH1D* th1d_nc_s2_east = ( (TProfile*) file->Get("fvtxs2_v3_east_docalib") )->ProjectionX("fvtxs2_v3_east_docalib_px2nc");
  th1d_nc_s2_both->SetLineColor(kBlack);
  th1d_nc_s2_west->SetLineColor(kRed);
  th1d_nc_s2_east->SetLineColor(kBlue);
  th1d_nc_s2_both->Draw();
  th1d_nc_s2_both->SetMaximum(0.005);
  th1d_nc_s2_both->SetMinimum(-0.005);
  th1d_nc_s2_both->GetYaxis()->SetTitle("v_{3} (no EP corr)");
  th1d_nc_s2_both->GetYaxis()->SetTitleOffset(1.5);
  th1d_nc_s2_both->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  th1d_nc_s2_west->Draw("same");
  th1d_nc_s2_east->Draw("same");
  leg->Draw();
  c1->Print("FigsCheckThree/fig_nocorrterms_layer_s2_allarms.png");
  c1->Print("FigsCheckThree/fig_nocorrterms_layer_s2_allarms.pdf");

  TH1D* th1d_nc_s3_both = ( (TProfile*) file->Get("fvtxs3_v3_both_docalib") )->ProjectionX("fvtxs3_v3_both_docalib_px2nc");
  TH1D* th1d_nc_s3_west = ( (TProfile*) file->Get("fvtxs3_v3_west_docalib") )->ProjectionX("fvtxs3_v3_west_docalib_px2nc");
  TH1D* th1d_nc_s3_east = ( (TProfile*) file->Get("fvtxs3_v3_east_docalib") )->ProjectionX("fvtxs3_v3_east_docalib_px2nc");
  th1d_nc_s3_both->SetLineColor(kBlack);
  th1d_nc_s3_west->SetLineColor(kRed);
  th1d_nc_s3_east->SetLineColor(kBlue);
  th1d_nc_s3_both->Draw();
  th1d_nc_s3_both->SetMaximum(0.005);
  th1d_nc_s3_both->SetMinimum(-0.005);
  th1d_nc_s3_both->GetYaxis()->SetTitle("v_{3} (no EP corr)");
  th1d_nc_s3_both->GetYaxis()->SetTitleOffset(1.5);
  th1d_nc_s3_both->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  th1d_nc_s3_west->Draw("same");
  th1d_nc_s3_east->Draw("same");
  leg->Draw();
  c1->Print("FigsCheckThree/fig_nocorrterms_layer_s3_allarms.png");
  c1->Print("FigsCheckThree/fig_nocorrterms_layer_s3_allarms.pdf");

  // --- Sein ist die hand die verletzt...

  int harmonic = 3;
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

  th1d_sB_both->Scale(1.0/reso_bbc);
  th1d_sB_west->Scale(1.0/reso_bbc);
  th1d_sB_east->Scale(1.0/reso_bbc);

  th1d_sA_both->Scale(1.0/reso_fvtx);
  th1d_sA_west->Scale(1.0/reso_fvtx);
  th1d_sA_east->Scale(1.0/reso_fvtx);

  th1d_sB_both->SetLineColor(kBlue+2);
  th1d_sA_both->SetLineColor(kRed);
  th1d_sB_both->Draw();
  th1d_sB_both->SetMaximum(0.05);
  th1d_sB_both->SetMinimum(-0.01);
  th1d_sA_both->Draw("same");
  c1->Print("FigsCheckThree/reso_corrected_botharms_bbcsfvtxs.png");
  c1->Print("FigsCheckThree/reso_corrected_botharms_bbcsfvtxs.pdf");

  th1d_nc_sB_both->Scale(1.0/reso_bbc);
  th1d_nc_sB_west->Scale(1.0/reso_bbc);
  th1d_nc_sB_east->Scale(1.0/reso_bbc);

  th1d_nc_sA_both->Scale(1.0/reso_fvtx);
  th1d_nc_sA_west->Scale(1.0/reso_fvtx);
  th1d_nc_sA_east->Scale(1.0/reso_fvtx);

  th1d_nc_sB_both->SetLineColor(kBlue+2);
  th1d_nc_sA_both->SetLineColor(kRed);
  th1d_nc_sB_both->Draw();
  th1d_nc_sB_both->SetMaximum(0.05);
  th1d_nc_sB_both->SetMinimum(-0.01);
  th1d_nc_sA_both->Draw("same");
  c1->Print("FigsCheckThree/reso_corrected_noterms_botharms_bbcsfvtxs.png");
  c1->Print("FigsCheckThree/reso_corrected_noterms_botharms_bbcsfvtxs.pdf");


  th1d_sB_both->SetLineColor(kBlue+2);
  th1d_sA_both->SetLineColor(kRed+2);
  th1d_nc_sB_both->SetLineColor(kBlue);
  th1d_nc_sA_both->SetLineColor(kRed);
  th1d_sB_both->Draw();
  th1d_sA_both->Draw("same");
  th1d_nc_sB_both->Draw("same");
  th1d_nc_sA_both->Draw("same");
  c1->Print("FigsCheckThree/reso_corrected_compterms_botharms_bbcsfvtxs.png");
  c1->Print("FigsCheckThree/reso_corrected_compterms_botharms_bbcsfvtxs.pdf");

  th1d_sB_east->SetLineColor(kBlue+2);
  th1d_sA_east->SetLineColor(kRed+2);
  th1d_nc_sB_east->SetLineColor(kBlue);
  th1d_nc_sA_east->SetLineColor(kRed);
  th1d_sB_east->Draw();
  th1d_sA_east->Draw("same");
  th1d_nc_sB_east->Draw("same");
  th1d_nc_sA_east->Draw("same");
  c1->Print("FigsCheckThree/reso_corrected_compterms_eastarm_bbcsfvtxs.png");
  c1->Print("FigsCheckThree/reso_corrected_compterms_eastarm_bbcsfvtxs.pdf");

  th1d_sB_west->SetLineColor(kBlue+2);
  th1d_sA_west->SetLineColor(kRed+2);
  th1d_nc_sB_west->SetLineColor(kBlue);
  th1d_nc_sA_west->SetLineColor(kRed);
  th1d_sB_west->Draw();
  th1d_sA_west->Draw("same");
  th1d_nc_sB_west->Draw("same");
  th1d_nc_sA_west->Draw("same");
  c1->Print("FigsCheckThree/reso_corrected_compterms_westarm_bbcsfvtxs.png");
  c1->Print("FigsCheckThree/reso_corrected_compterms_westarm_bbcsfvtxs.pdf");

  TFile* fout = TFile::Open("shengli_stuff.root","recreate");
  fout->cd();
  th1d_sB_both->SetName("hvn_withcorrterms_bbcs_B");
  th1d_sB_east->SetName("hvn_withcorrterms_bbcs_E");
  th1d_sB_west->SetName("hvn_withcorrterms_bbcs_W");
  th1d_sA_both->SetName("hvn_withcorrterms_fvtxs_B");
  th1d_sA_east->SetName("hvn_withcorrterms_fvtxs_E");
  th1d_sA_west->SetName("hvn_withcorrterms_fvtxs_W");
  th1d_s0_both->SetName("hvn_withcorrterms_fvtxs0_B");
  th1d_s0_east->SetName("hvn_withcorrterms_fvtxs0_E");
  th1d_s0_west->SetName("hvn_withcorrterms_fvtxs0_W");
  th1d_s1_both->SetName("hvn_withcorrterms_fvtxs1_B");
  th1d_s1_east->SetName("hvn_withcorrterms_fvtxs1_E");
  th1d_s1_west->SetName("hvn_withcorrterms_fvtxs1_W");
  th1d_s2_both->SetName("hvn_withcorrterms_fvtxs2_B");
  th1d_s2_east->SetName("hvn_withcorrterms_fvtxs2_E");
  th1d_s2_west->SetName("hvn_withcorrterms_fvtxs2_W");
  th1d_s3_both->SetName("hvn_withcorrterms_fvtxs3_B");
  th1d_s3_east->SetName("hvn_withcorrterms_fvtxs3_E");
  th1d_s3_west->SetName("hvn_withcorrterms_fvtxs3_W");
  th1d_sB_both->Write();
  th1d_sB_east->Write();
  th1d_sB_west->Write();
  th1d_sA_both->Write();
  th1d_sA_east->Write();
  th1d_sA_west->Write();
  th1d_s0_both->Write();
  th1d_s0_east->Write();
  th1d_s0_west->Write();
  th1d_s1_both->Write();
  th1d_s1_east->Write();
  th1d_s1_west->Write();
  th1d_s2_both->Write();
  th1d_s2_east->Write();
  th1d_s2_west->Write();
  th1d_s3_both->Write();
  th1d_s3_east->Write();
  th1d_s3_west->Write();
  th1d_nc_sB_both->SetName("hvn_nocorrterms_bbcs_B");
  th1d_nc_sB_east->SetName("hvn_nocorrterms_bbcs_E");
  th1d_nc_sB_west->SetName("hvn_nocorrterms_bbcs_W");
  th1d_nc_sA_both->SetName("hvn_nocorrterms_fvtxs_B");
  th1d_nc_sA_east->SetName("hvn_nocorrterms_fvtxs_E");
  th1d_nc_sA_west->SetName("hvn_nocorrterms_fvtxs_W");
  th1d_nc_s0_both->SetName("hvn_nocorrterms_fvtxs0_B");
  th1d_nc_s0_east->SetName("hvn_nocorrterms_fvtxs0_E");
  th1d_nc_s0_west->SetName("hvn_nocorrterms_fvtxs0_W");
  th1d_nc_s1_both->SetName("hvn_nocorrterms_fvtxs1_B");
  th1d_nc_s1_east->SetName("hvn_nocorrterms_fvtxs1_E");
  th1d_nc_s1_west->SetName("hvn_nocorrterms_fvtxs1_W");
  th1d_nc_s2_both->SetName("hvn_nocorrterms_fvtxs2_B");
  th1d_nc_s2_east->SetName("hvn_nocorrterms_fvtxs2_E");
  th1d_nc_s2_west->SetName("hvn_nocorrterms_fvtxs2_W");
  th1d_nc_s3_both->SetName("hvn_nocorrterms_fvtxs3_B");
  th1d_nc_s3_east->SetName("hvn_nocorrterms_fvtxs3_E");
  th1d_nc_s3_west->SetName("hvn_nocorrterms_fvtxs3_W");
  th1d_nc_sB_both->Write();
  th1d_nc_sB_east->Write();
  th1d_nc_sB_west->Write();
  th1d_nc_sA_both->Write();
  th1d_nc_sA_east->Write();
  th1d_nc_sA_west->Write();
  th1d_nc_s0_both->Write();
  th1d_nc_s0_east->Write();
  th1d_nc_s0_west->Write();
  th1d_nc_s1_both->Write();
  th1d_nc_s1_east->Write();
  th1d_nc_s1_west->Write();
  th1d_nc_s2_both->Write();
  th1d_nc_s2_east->Write();
  th1d_nc_s2_west->Write();
  th1d_nc_s3_both->Write();
  th1d_nc_s3_east->Write();
  th1d_nc_s3_west->Write();
  fout->Write();
  fout->Close();

}


TH1D* dueprocess(TProfile* tp1f_base, TProfile* tp1f_cosphi, TProfile* tp1f_sinphi, TProfile* tp1f_cospsi, TProfile* tp1f_sinpsi)
{
  TH1D* th1d_base = tp1f_base->ProjectionX();
  TH1D* th1d_cosphi = tp1f_cosphi->ProjectionX();
  TH1D* th1d_sinphi = tp1f_sinphi->ProjectionX();
  TH1D* th1d_cospsi = tp1f_cospsi->ProjectionX();
  TH1D* th1d_sinpsi = tp1f_sinpsi->ProjectionX();
  TH1D* result = dueprocess(th1d_base,th1d_cosphi,th1d_sinphi,th1d_cospsi,th1d_sinpsi);
  cout << result << endl;
  return result;
}


TH1D* dueprocess(TH1D* th1d_base, TH1D* th1d_cosphi, TH1D* th1d_sinphi, TH1D* th1d_cospsi, TH1D* th1d_sinpsi)
{
  TH1D* th1d_combined_cos = (TH1D*)th1d_cosphi->Clone();
  TH1D* th1d_combined_sin = (TH1D*)th1d_sinphi->Clone();
  th1d_combined_cos->Multiply(th1d_cospsi);
  th1d_combined_sin->Multiply(th1d_sinpsi);
  TH1D* th1d_return = (TH1D*)th1d_base->Clone();
  th1d_return->Add(th1d_combined_cos,-1.0);
  th1d_return->Add(th1d_combined_sin,-1.0);
  cout << th1d_return << endl;
  return th1d_return;
}
