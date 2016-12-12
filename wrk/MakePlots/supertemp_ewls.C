void supertemp_ewls()
{

  gStyle->SetOptTitle(1);

  TFile* file = new TFile("input/combined_200.root");

  TProfile* hv3_bbcs_both = (TProfile*)file->Get("bbcs_v3_both_docalib");
  TProfile* hv3_fvtxs_both = (TProfile*)file->Get("fvtxs_v3_both_docalib");
  TProfile* hv3_layer0_both = (TProfile*)file->Get("fvtxs0_v3_both_docalib");
  TProfile* hv3_layer1_both = (TProfile*)file->Get("fvtxs1_v3_both_docalib");
  TProfile* hv3_layer2_both = (TProfile*)file->Get("fvtxs2_v3_both_docalib");
  TProfile* hv3_layer3_both = (TProfile*)file->Get("fvtxs3_v3_both_docalib");

  TProfile* hv3_bbcs_east = (TProfile*)file->Get("bbcs_v3_east_docalib");
  TProfile* hv3_fvtxs_east = (TProfile*)file->Get("fvtxs_v3_east_docalib");
  TProfile* hv3_layer0_east = (TProfile*)file->Get("fvtxs0_v3_east_docalib");
  TProfile* hv3_layer1_east = (TProfile*)file->Get("fvtxs1_v3_east_docalib");
  TProfile* hv3_layer2_east = (TProfile*)file->Get("fvtxs2_v3_east_docalib");
  TProfile* hv3_layer3_east = (TProfile*)file->Get("fvtxs3_v3_east_docalib");

  TProfile* hv3_bbcs_west = (TProfile*)file->Get("bbcs_v3_west_docalib");
  TProfile* hv3_fvtxs_west = (TProfile*)file->Get("fvtxs_v3_west_docalib");
  TProfile* hv3_layer0_west = (TProfile*)file->Get("fvtxs0_v3_west_docalib");
  TProfile* hv3_layer1_west = (TProfile*)file->Get("fvtxs1_v3_west_docalib");
  TProfile* hv3_layer2_west = (TProfile*)file->Get("fvtxs2_v3_west_docalib");
  TProfile* hv3_layer3_west = (TProfile*)file->Get("fvtxs3_v3_west_docalib");

  hv3_bbcs_both->SetLineColor(kBlack);
  hv3_fvtxs_both->SetLineColor(kBlack);
  hv3_layer0_both->SetLineColor(kBlack);
  hv3_layer1_both->SetLineColor(kBlack);
  hv3_layer2_both->SetLineColor(kBlack);
  hv3_layer3_both->SetLineColor(kBlack);

  hv3_bbcs_east->SetLineColor(kRed);
  hv3_fvtxs_east->SetLineColor(kRed);
  hv3_layer0_east->SetLineColor(kRed);
  hv3_layer1_east->SetLineColor(kRed);
  hv3_layer2_east->SetLineColor(kRed);
  hv3_layer3_east->SetLineColor(kRed);

  hv3_bbcs_west->SetLineColor(kBlue);
  hv3_fvtxs_west->SetLineColor(kBlue);
  hv3_layer0_west->SetLineColor(kBlue);
  hv3_layer1_west->SetLineColor(kBlue);
  hv3_layer2_west->SetLineColor(kBlue);
  hv3_layer3_west->SetLineColor(kBlue);

  hv3_bbcs_both->Draw();
  hv3_bbcs_both->SetTitle("BBC South");
  hv3_bbcs_both->SetMaximum(0.005);
  hv3_bbcs_both->SetMinimum(-0.005);
  hv3_bbcs_both->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hv3_bbcs_both->GetYaxis()->SetTitle("v_{3} (no EP res corr)");
  hv3_bbcs_both->GetYaxis()->SetTitleOffset(1.25);
  hv3_bbcs_east->Draw("same");
  hv3_bbcs_west->Draw("same");
  TLegend* leg = new TLegend(0.18,0.68,0.38,0.88);
  leg->AddEntry(hv3_bbcs_both,"both arms","el");
  leg->AddEntry(hv3_bbcs_east,"east arms","el");
  leg->AddEntry(hv3_bbcs_west,"west arms","el");
  leg->Draw();
  TLine* line = new TLine(0.0,0.0,3.0,0.0);
  line->SetLineStyle(2);
  line->SetLineWidth(2);
  line->Draw();
  c1->Print("FigsCheckThree/supertemp_bbcs_bew.png");
  c1->Print("FigsCheckThree/supertemp_bbcs_bew.pdf");

  hv3_fvtxs_both->Draw();
  hv3_fvtxs_both->SetTitle("FVTX South, all layers");
  hv3_fvtxs_both->SetMaximum(0.005);
  hv3_fvtxs_both->SetMinimum(-0.005);
  hv3_fvtxs_both->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hv3_fvtxs_both->GetYaxis()->SetTitle("v_{3} (no EP res corr)");
  hv3_fvtxs_both->GetYaxis()->SetTitleOffset(1.25);
  hv3_fvtxs_east->Draw("same");
  hv3_fvtxs_west->Draw("same");
  leg->Draw();
  line->Draw();
  c1->Print("FigsCheckThree/supertemp_fvtxs_bew.png");
  c1->Print("FigsCheckThree/supertemp_fvtxs_bew.pdf");

  hv3_layer0_both->Draw();
  hv3_layer0_both->SetTitle("FVTX South, layer 0");
  hv3_layer0_both->SetMaximum(0.005);
  hv3_layer0_both->SetMinimum(-0.005);
  hv3_layer0_both->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hv3_layer0_both->GetYaxis()->SetTitle("v_{3} (no EP res corr)");
  hv3_layer0_both->GetYaxis()->SetTitleOffset(1.25);
  hv3_layer0_east->Draw("same");
  hv3_layer0_west->Draw("same");
  leg->Draw();
  line->Draw();
  c1->Print("FigsCheckThree/supertemp_layer0_bew.png");
  c1->Print("FigsCheckThree/supertemp_layer0_bew.pdf");

  hv3_layer1_both->Draw();
  hv3_layer1_both->SetTitle("FVTX South, layer 1");
  hv3_layer1_both->SetMaximum(0.005);
  hv3_layer1_both->SetMinimum(-0.005);
  hv3_layer1_both->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hv3_layer1_both->GetYaxis()->SetTitle("v_{3} (no EP res corr)");
  hv3_layer1_both->GetYaxis()->SetTitleOffset(1.25);
  hv3_layer1_east->Draw("same");
  hv3_layer1_west->Draw("same");
  leg->Draw();
  line->Draw();
  c1->Print("FigsCheckThree/supertemp_layer1_bew.png");
  c1->Print("FigsCheckThree/supertemp_layer1_bew.pdf");

  hv3_layer2_both->Draw();
  hv3_layer2_both->SetTitle("FVTX South, layer 2");
  hv3_layer2_both->SetMaximum(0.005);
  hv3_layer2_both->SetMinimum(-0.005);
  hv3_layer2_both->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hv3_layer2_both->GetYaxis()->SetTitle("v_{3} (no EP res corr)");
  hv3_layer2_both->GetYaxis()->SetTitleOffset(1.25);
  hv3_layer2_east->Draw("same");
  hv3_layer2_west->Draw("same");
  leg->Draw();
  line->Draw();
  c1->Print("FigsCheckThree/supertemp_layer2_bew.png");
  c1->Print("FigsCheckThree/supertemp_layer2_bew.pdf");

  hv3_layer3_both->Draw();
  hv3_layer3_both->SetTitle("FVTX South, layer 3");
  hv3_layer3_both->SetMaximum(0.005);
  hv3_layer3_both->SetMinimum(-0.005);
  hv3_layer3_both->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hv3_layer3_both->GetYaxis()->SetTitle("v_{3} (no EP res corr)");
  hv3_layer3_both->GetYaxis()->SetTitleOffset(1.25);
  hv3_layer3_east->Draw("same");
  hv3_layer3_west->Draw("same");
  leg->Draw();
  line->Draw();
  c1->Print("FigsCheckThree/supertemp_layer3_bew.png");
  c1->Print("FigsCheckThree/supertemp_layer3_bew.pdf");

}
