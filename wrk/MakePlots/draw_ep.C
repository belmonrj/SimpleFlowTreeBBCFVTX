void draw_ep()
{

  TCanvas* c1 = new TCanvas("c1","");

  TFile* file = TFile::Open("input/EverythingIsCopasetic/hrp_455050.root");

  TH1D* psi2_BBC_before = (TH1D*)((TH2D*)file->Get("psi_bf_0_1_2"))->ProjectionY();
  TH1D* psi2_BBC_after = (TH1D*)((TH2D*)file->Get("psi_af_0_1_2"))->ProjectionY();
  psi2_BBC_before->Draw();
  psi2_BBC_before->GetXaxis()->SetRangeUser(-2.0,2.0);
  psi2_BBC_after->SetLineColor(kRed);
  psi2_BBC_after->Draw("same");
  c1->Print("flattening_psi2_BBC_455050.png");
  c1->Print("flattening_psi2_BBC_455050.pdf");

  TH1D* psi2_FVTX_before = (TH1D*)((TH2D*)file->Get("psi_bf_0_1_3"))->ProjectionY();
  TH1D* psi2_FVTX_after = (TH1D*)((TH2D*)file->Get("psi_af_0_1_3"))->ProjectionY();
  psi2_FVTX_before->Draw();
  psi2_FVTX_before->GetXaxis()->SetRangeUser(-2.0,2.0);
  psi2_FVTX_after->SetLineColor(kRed);
  psi2_FVTX_after->Draw("same");
  c1->Print("flattening_psi2_FVTX_455050.png");
  c1->Print("flattening_psi2_FVTX_455050.pdf");

  // ---

  TH1D* psi3_BBC_before = (TH1D*)((TH2D*)file->Get("psi_bf_0_2_2"))->ProjectionY();
  TH1D* psi3_BBC_after = (TH1D*)((TH2D*)file->Get("psi_af_0_2_2"))->ProjectionY();
  psi3_BBC_before->Draw();
  psi3_BBC_before->GetXaxis()->SetRangeUser(-2.0,2.0);
  psi3_BBC_after->SetLineColor(kRed);
  psi3_BBC_after->Draw("same");
  c1->Print("flattening_psi3_BBC_455050.png");
  c1->Print("flattening_psi3_BBC_455050.pdf");

  TH1D* psi3_FVTX_before = (TH1D*)((TH2D*)file->Get("psi_bf_0_2_3"))->ProjectionY();
  TH1D* psi3_FVTX_after = (TH1D*)((TH2D*)file->Get("psi_af_0_2_3"))->ProjectionY();
  psi3_FVTX_before->Draw();
  psi3_FVTX_before->GetXaxis()->SetRangeUser(-2.0,2.0);
  psi3_FVTX_after->SetLineColor(kRed);
  psi3_FVTX_after->Draw("same");
  c1->Print("flattening_psi3_FVTX_455050.png");
  c1->Print("flattening_psi3_FVTX_455050.pdf");

}
