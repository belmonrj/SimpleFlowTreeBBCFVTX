TH1D* dueprocess(TProfile*,TProfile*,TProfile*,TProfile*,TProfile*);

TH1D* dueprocess(TH1D*,TH1D*,TH1D*,TH1D*,TH1D*);

void diehanddieverletzt()
{

  TFile* file = TFile::Open("input/combined_200.root");

  TProfile* tp1f_base;
  TProfile* tp1f_cosphi;
  TProfile* tp1f_sinphi;
  TProfile* tp1f_cospsi;
  TProfile* tp1f_sinpsi;

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

  th1d_sA_both->SetLineColor(kBlack);
  th1d_sA_west->SetLineColor(kRed);
  th1d_sA_east->SetLineColor(kBlue);
  th1d_sA_both->Draw();
  th1d_sA_both->SetMaximum(0.002);
  th1d_sA_both->SetMinimum(-0.002);
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
  th1d_s0_both->SetMaximum(0.002);
  th1d_s0_both->SetMinimum(-0.002);
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
  th1d_s1_both->SetMaximum(0.002);
  th1d_s1_both->SetMinimum(-0.002);
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
  th1d_s2_both->SetMaximum(0.002);
  th1d_s2_both->SetMinimum(-0.002);
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
  th1d_s3_both->SetMaximum(0.002);
  th1d_s3_both->SetMinimum(-0.002);
  th1d_s3_both->GetYaxis()->SetTitle("v_{3} (no EP corr)");
  th1d_s3_both->GetYaxis()->SetTitleOffset(1.5);
  th1d_s3_both->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  th1d_s3_west->Draw("same");
  th1d_s3_east->Draw("same");
  leg->Draw();
  c1->Print("FigsCheckThree/fig_layer_s3_allarms.png");
  c1->Print("FigsCheckThree/fig_layer_s3_allarms.pdf");

  // ---

  th1d_sA_both = ( (TProfile*) file->Get("fvtxs_v3_both_docalib") )->ProjectionX("fvtxs_v3_both_docalib_px2nc");
  th1d_sA_west = ( (TProfile*) file->Get("fvtxs_v3_west_docalib") )->ProjectionX("fvtxs_v3_west_docalib_px2nc");
  th1d_sA_east = ( (TProfile*) file->Get("fvtxs_v3_east_docalib") )->ProjectionX("fvtxs_v3_east_docalib_px2nc");
  th1d_sA_both->SetLineColor(kBlack);
  th1d_sA_west->SetLineColor(kRed);
  th1d_sA_east->SetLineColor(kBlue);
  th1d_sA_both->Draw();
  th1d_sA_both->SetMaximum(0.002);
  th1d_sA_both->SetMinimum(-0.002);
  th1d_sA_both->GetYaxis()->SetTitle("v_{3} (no EP corr)");
  th1d_sA_both->GetYaxis()->SetTitleOffset(1.5);
  th1d_sA_both->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  th1d_sA_west->Draw("same");
  th1d_sA_east->Draw("same");
  leg->Draw();
  c1->Print("FigsCheckThree/fig_nocorrterms_layer_sA_allarms.png");
  c1->Print("FigsCheckThree/fig_nocorrterms_layer_sA_allarms.pdf");

  th1d_s0_both = ( (TProfile*) file->Get("fvtxs0_v3_both_docalib") )->ProjectionX("fvtxs0_v3_both_docalib_px2nc");
  th1d_s0_west = ( (TProfile*) file->Get("fvtxs0_v3_west_docalib") )->ProjectionX("fvtxs0_v3_west_docalib_px2nc");
  th1d_s0_east = ( (TProfile*) file->Get("fvtxs0_v3_east_docalib") )->ProjectionX("fvtxs0_v3_east_docalib_px2nc");
  th1d_s0_both->SetLineColor(kBlack);
  th1d_s0_west->SetLineColor(kRed);
  th1d_s0_east->SetLineColor(kBlue);
  th1d_s0_both->Draw();
  th1d_s0_both->SetMaximum(0.002);
  th1d_s0_both->SetMinimum(-0.002);
  th1d_s0_both->GetYaxis()->SetTitle("v_{3} (no EP corr)");
  th1d_s0_both->GetYaxis()->SetTitleOffset(1.5);
  th1d_s0_both->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  th1d_s0_west->Draw("same");
  th1d_s0_east->Draw("same");
  leg->Draw();
  c1->Print("FigsCheckThree/fig_nocorrterms_layer_s0_allarms.png");
  c1->Print("FigsCheckThree/fig_nocorrterms_layer_s0_allarms.pdf");

  th1d_s1_both = ( (TProfile*) file->Get("fvtxs1_v3_both_docalib") )->ProjectionX("fvtxs1_v3_both_docalib_px2nc");
  th1d_s1_west = ( (TProfile*) file->Get("fvtxs1_v3_west_docalib") )->ProjectionX("fvtxs1_v3_west_docalib_px2nc");
  th1d_s1_east = ( (TProfile*) file->Get("fvtxs1_v3_east_docalib") )->ProjectionX("fvtxs1_v3_east_docalib_px2nc");
  th1d_s1_both->SetLineColor(kBlack);
  th1d_s1_west->SetLineColor(kRed);
  th1d_s1_east->SetLineColor(kBlue);
  th1d_s1_both->Draw();
  th1d_s1_both->SetMaximum(0.002);
  th1d_s1_both->SetMinimum(-0.002);
  th1d_s1_both->GetYaxis()->SetTitle("v_{3} (no EP corr)");
  th1d_s1_both->GetYaxis()->SetTitleOffset(1.5);
  th1d_s1_both->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  th1d_s1_west->Draw("same");
  th1d_s1_east->Draw("same");
  leg->Draw();
  c1->Print("FigsCheckThree/fig_nocorrterms_layer_s1_allarms.png");
  c1->Print("FigsCheckThree/fig_nocorrterms_layer_s1_allarms.pdf");

  th1d_s2_both = ( (TProfile*) file->Get("fvtxs2_v3_both_docalib") )->ProjectionX("fvtxs2_v3_both_docalib_px2nc");
  th1d_s2_west = ( (TProfile*) file->Get("fvtxs2_v3_west_docalib") )->ProjectionX("fvtxs2_v3_west_docalib_px2nc");
  th1d_s2_east = ( (TProfile*) file->Get("fvtxs2_v3_east_docalib") )->ProjectionX("fvtxs2_v3_east_docalib_px2nc");
  th1d_s2_both->SetLineColor(kBlack);
  th1d_s2_west->SetLineColor(kRed);
  th1d_s2_east->SetLineColor(kBlue);
  th1d_s2_both->Draw();
  th1d_s2_both->SetMaximum(0.002);
  th1d_s2_both->SetMinimum(-0.002);
  th1d_s2_both->GetYaxis()->SetTitle("v_{3} (no EP corr)");
  th1d_s2_both->GetYaxis()->SetTitleOffset(1.5);
  th1d_s2_both->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  th1d_s2_west->Draw("same");
  th1d_s2_east->Draw("same");
  leg->Draw();
  c1->Print("FigsCheckThree/fig_nocorrterms_layer_s2_allarms.png");
  c1->Print("FigsCheckThree/fig_nocorrterms_layer_s2_allarms.pdf");

  th1d_s3_both = ( (TProfile*) file->Get("fvtxs3_v3_both_docalib") )->ProjectionX("fvtxs3_v3_both_docalib_px2nc");
  th1d_s3_west = ( (TProfile*) file->Get("fvtxs3_v3_west_docalib") )->ProjectionX("fvtxs3_v3_west_docalib_px2nc");
  th1d_s3_east = ( (TProfile*) file->Get("fvtxs3_v3_east_docalib") )->ProjectionX("fvtxs3_v3_east_docalib_px2nc");
  th1d_s3_both->SetLineColor(kBlack);
  th1d_s3_west->SetLineColor(kRed);
  th1d_s3_east->SetLineColor(kBlue);
  th1d_s3_both->Draw();
  th1d_s3_both->SetMaximum(0.002);
  th1d_s3_both->SetMinimum(-0.002);
  th1d_s3_both->GetYaxis()->SetTitle("v_{3} (no EP corr)");
  th1d_s3_both->GetYaxis()->SetTitleOffset(1.5);
  th1d_s3_both->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  th1d_s3_west->Draw("same");
  th1d_s3_east->Draw("same");
  leg->Draw();
  c1->Print("FigsCheckThree/fig_nocorrterms_layer_s3_allarms.png");
  c1->Print("FigsCheckThree/fig_nocorrterms_layer_s3_allarms.pdf");


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
