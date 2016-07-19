void takerun(int);
void takeenergy(int);
void takefile(TFile*, int);

void compgb();
void compE();

void special_x1()
{

  gStyle->SetOptTitle(1);

  // compbg();
  // compE();
  //return;

  // takerun(456652);
  // takerun(456943);
  // takerun(457013);
  // takerun(457298);

  takeenergy(20);
  takeenergy(39);
  takeenergy(62);

}


void takerun(int run)
{
  TFile* file = TFile::Open(Form("RootFiles/svrb_run%d_pass0.root",run));
  takefile(file,run);
}


void takeenergy(int energy)
{
  TFile* file = TFile::Open(Form("RootFiles/sum%d.root",energy));
  takefile(file,energy);
}


void takefile(TFile* file, int handle)
{

  TCanvas* c1 = new TCanvas("c1","");

  TH1D* th1d_counter = (TH1D*)file->Get("th1d_FVTX_nclus");
  TH1D* th1d_counter_IR = (TH1D*)file->Get("th1d_FVTX_nclus_IR");
  TH1D* th1d_counter_NCIR = (TH1D*)file->Get("th1d_FVTX_nclus_NCIR");
  TH1D* th1d_counter_OR = (TH1D*)file->Get("th1d_FVTX_nclus_OR");
  int nevents = th1d_counter->GetEntries();
  int nevents_IR = th1d_counter_IR->GetEntries();
  int nevents_NCIR = th1d_counter_NCIR->GetEntries();
  int nevents_OR = th1d_counter_OR->GetEntries();

  cout << nevents << endl;
  cout << nevents_IR << endl;
  cout << nevents_NCIR << endl;

  TH1D* th1d_fvtxs_clus_phi = (TH1D*)file->Get("th1d_fvtxs_clus_phi");
  TH1D* th1d_fvtxs_clus_phi_IR = (TH1D*)file->Get("th1d_fvtxs_clus_phi_IR");
  TH1D* th1d_fvtxs_clus_phi_NCIR = (TH1D*)file->Get("th1d_fvtxs_clus_phi_NCIR");
  th1d_fvtxs_clus_phi->Scale(1.0/nevents);
  th1d_fvtxs_clus_phi_IR->Scale(1.0/nevents_IR);
  th1d_fvtxs_clus_phi_NCIR->Scale(1.0/nevents_NCIR);
  th1d_fvtxs_clus_phi->SetMinimum(0);
  th1d_fvtxs_clus_phi->SetTitle("all layers, FVTX South");
  th1d_fvtxs_clus_phi->GetXaxis()->SetTitle("cluster #phi (rad)");
  th1d_fvtxs_clus_phi->GetYaxis()->SetTitle("counts/event");
  th1d_fvtxs_clus_phi->SetLineColor(kBlack);
  th1d_fvtxs_clus_phi->Draw();
  th1d_fvtxs_clus_phi_NCIR->SetLineColor(kRed);
  th1d_fvtxs_clus_phi_NCIR->Draw("same");
  th1d_fvtxs_clus_phi_IR->SetLineColor(kBlue);
  th1d_fvtxs_clus_phi_IR->Draw("same");
  TLegend* leg = new TLegend(0.18,0.68,0.38,0.88);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.045);
  leg->AddEntry(th1d_fvtxs_clus_phi,"all","l");
  leg->AddEntry(th1d_fvtxs_clus_phi_IR,"inside radius","l");
  leg->AddEntry(th1d_fvtxs_clus_phi_NCIR,"inside radius and cluster cut","l");
  leg->Draw();
  c1->Print(Form("FigsPhi/radiuscut_fvtxs_clus_phi_%d.png",handle));
  c1->Print(Form("FigsPhi/radiuscut_fvtxs_clus_phi_%d.pdf",handle));

  // ---

  TH1D* th1d_fvtxs0_clus_phi = (TH1D*)file->Get("th1d_fvtxs0_clus_phi");
  TH1D* th1d_fvtxs0_clus_phi_IR = (TH1D*)file->Get("th1d_fvtxs0_clus_phi_IR");
  TH1D* th1d_fvtxs0_clus_phi_NCIR = (TH1D*)file->Get("th1d_fvtxs0_clus_phi_NCIR");
  th1d_fvtxs0_clus_phi->Scale(1.0/nevents);
  th1d_fvtxs0_clus_phi_IR->Scale(1.0/nevents_IR);
  th1d_fvtxs0_clus_phi_NCIR->Scale(1.0/nevents_NCIR);
  th1d_fvtxs0_clus_phi->SetMinimum(0);
  th1d_fvtxs0_clus_phi->SetTitle("layer 0, FVTX South");
  th1d_fvtxs0_clus_phi->GetXaxis()->SetTitle("cluster #phi (rad)");
  th1d_fvtxs0_clus_phi->GetYaxis()->SetTitle("counts/event");
  th1d_fvtxs0_clus_phi->SetLineColor(kBlack);
  th1d_fvtxs0_clus_phi->Draw();
  th1d_fvtxs0_clus_phi_NCIR->SetLineColor(kRed);
  th1d_fvtxs0_clus_phi_NCIR->Draw("same");
  th1d_fvtxs0_clus_phi_IR->SetLineColor(kBlue);
  th1d_fvtxs0_clus_phi_IR->Draw("same");
  leg->Draw();
  c1->Print(Form("FigsPhi/radiuscut_fvtxs0_clus_phi_%d.png",handle));
  c1->Print(Form("FigsPhi/radiuscut_fvtxs0_clus_phi_%d.pdf",handle));

  // ---

  TH1D* th1d_fvtxs1_clus_phi = (TH1D*)file->Get("th1d_fvtxs1_clus_phi");
  TH1D* th1d_fvtxs1_clus_phi_IR = (TH1D*)file->Get("th1d_fvtxs1_clus_phi_IR");
  TH1D* th1d_fvtxs1_clus_phi_NCIR = (TH1D*)file->Get("th1d_fvtxs1_clus_phi_NCIR");
  th1d_fvtxs1_clus_phi->Scale(1.0/nevents);
  th1d_fvtxs1_clus_phi_IR->Scale(1.0/nevents_IR);
  th1d_fvtxs1_clus_phi_NCIR->Scale(1.0/nevents_NCIR);
  th1d_fvtxs1_clus_phi->SetMinimum(0);
  th1d_fvtxs1_clus_phi->SetTitle("layer 1, FVTX South");
  th1d_fvtxs1_clus_phi->GetXaxis()->SetTitle("cluster #phi (rad)");
  th1d_fvtxs1_clus_phi->GetYaxis()->SetTitle("counts/event");
  th1d_fvtxs1_clus_phi->SetLineColor(kBlack);
  th1d_fvtxs1_clus_phi->Draw();
  th1d_fvtxs1_clus_phi_NCIR->SetLineColor(kRed);
  th1d_fvtxs1_clus_phi_NCIR->Draw("same");
  th1d_fvtxs1_clus_phi_IR->SetLineColor(kBlue);
  th1d_fvtxs1_clus_phi_IR->Draw("same");
  leg->Draw();
  c1->Print(Form("FigsPhi/radiuscut_fvtxs1_clus_phi_%d.png",handle));
  c1->Print(Form("FigsPhi/radiuscut_fvtxs1_clus_phi_%d.pdf",handle));

  // ---

  TH1D* th1d_fvtxs2_clus_phi = (TH1D*)file->Get("th1d_fvtxs2_clus_phi");
  TH1D* th1d_fvtxs2_clus_phi_IR = (TH1D*)file->Get("th1d_fvtxs2_clus_phi_IR");
  TH1D* th1d_fvtxs2_clus_phi_NCIR = (TH1D*)file->Get("th1d_fvtxs2_clus_phi_NCIR");
  th1d_fvtxs2_clus_phi->Scale(1.0/nevents);
  th1d_fvtxs2_clus_phi_IR->Scale(1.0/nevents_IR);
  th1d_fvtxs2_clus_phi_NCIR->Scale(1.0/nevents_NCIR);
  th1d_fvtxs2_clus_phi->SetMinimum(0);
  th1d_fvtxs2_clus_phi->SetTitle("layer 2, FVTX South");
  th1d_fvtxs2_clus_phi->GetXaxis()->SetTitle("cluster #phi (rad)");
  th1d_fvtxs2_clus_phi->GetYaxis()->SetTitle("counts/event");
  th1d_fvtxs2_clus_phi->SetLineColor(kBlack);
  th1d_fvtxs2_clus_phi->Draw();
  th1d_fvtxs2_clus_phi_NCIR->SetLineColor(kRed);
  th1d_fvtxs2_clus_phi_NCIR->Draw("same");
  th1d_fvtxs2_clus_phi_IR->SetLineColor(kBlue);
  th1d_fvtxs2_clus_phi_IR->Draw("same");
  leg->Draw();
  c1->Print(Form("FigsPhi/radiuscut_fvtxs2_clus_phi_%d.png",handle));
  c1->Print(Form("FigsPhi/radiuscut_fvtxs2_clus_phi_%d.pdf",handle));

  // ---

  TH1D* th1d_fvtxs3_clus_phi = (TH1D*)file->Get("th1d_fvtxs3_clus_phi");
  TH1D* th1d_fvtxs3_clus_phi_IR = (TH1D*)file->Get("th1d_fvtxs3_clus_phi_IR");
  TH1D* th1d_fvtxs3_clus_phi_NCIR = (TH1D*)file->Get("th1d_fvtxs3_clus_phi_NCIR");
  th1d_fvtxs3_clus_phi->Scale(1.0/nevents);
  th1d_fvtxs3_clus_phi_IR->Scale(1.0/nevents_IR);
  th1d_fvtxs3_clus_phi_NCIR->Scale(1.0/nevents_NCIR);
  th1d_fvtxs3_clus_phi->SetMinimum(0);
  th1d_fvtxs3_clus_phi->SetTitle("layer 3, FVTX South");
  th1d_fvtxs3_clus_phi->GetXaxis()->SetTitle("cluster #phi (rad)");
  th1d_fvtxs3_clus_phi->GetYaxis()->SetTitle("counts/event");
  th1d_fvtxs3_clus_phi->SetLineColor(kBlack);
  th1d_fvtxs3_clus_phi->Draw();
  th1d_fvtxs3_clus_phi_NCIR->SetLineColor(kRed);
  th1d_fvtxs3_clus_phi_NCIR->Draw("same");
  th1d_fvtxs3_clus_phi_IR->SetLineColor(kBlue);
  th1d_fvtxs3_clus_phi_IR->Draw("same");
  leg->Draw();
  c1->Print(Form("FigsPhi/radiuscut_fvtxs3_clus_phi_%d.png",handle));
  c1->Print(Form("FigsPhi/radiuscut_fvtxs3_clus_phi_%d.pdf",handle));



  // ---
  // --- come back here for north

  TH1D* th1d_fvtxn_clus_phi = (TH1D*)file->Get("th1d_fvtxn_clus_phi");
  TH1D* th1d_fvtxn_clus_phi_IR = (TH1D*)file->Get("th1d_fvtxn_clus_phi_IR");
  TH1D* th1d_fvtxn_clus_phi_NCIR = (TH1D*)file->Get("th1d_fvtxn_clus_phi_NCIR");
  th1d_fvtxn_clus_phi->Scale(1.0/nevents);
  th1d_fvtxn_clus_phi_IR->Scale(1.0/nevents_IR);
  th1d_fvtxn_clus_phi_NCIR->Scale(1.0/nevents_NCIR);
  th1d_fvtxn_clus_phi->SetMinimum(0);
  th1d_fvtxn_clus_phi->SetTitle("all layers, FVTX North");
  th1d_fvtxn_clus_phi->GetXaxis()->SetTitle("cluster #phi (rad)");
  th1d_fvtxn_clus_phi->GetYaxis()->SetTitle("counts/event");
  th1d_fvtxn_clus_phi->SetLineColor(kBlack);
  th1d_fvtxn_clus_phi->Draw();
  th1d_fvtxn_clus_phi_NCIR->SetLineColor(kRed);
  th1d_fvtxn_clus_phi_NCIR->Draw("same");
  th1d_fvtxn_clus_phi_IR->SetLineColor(kBlue);
  th1d_fvtxn_clus_phi_IR->Draw("same");
  TLegend* leg = new TLegend(0.18,0.68,0.38,0.88);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.045);
  leg->AddEntry(th1d_fvtxn_clus_phi,"all","l");
  leg->AddEntry(th1d_fvtxn_clus_phi_IR,"inside radius","l");
  leg->AddEntry(th1d_fvtxn_clus_phi_NCIR,"outside radius","l");
  leg->Draw();
  c1->Print(Form("FigsPhi/radiuscut_fvtxn_clus_phi_%d.png",handle));
  c1->Print(Form("FigsPhi/radiuscut_fvtxn_clus_phi_%d.pdf",handle));

  // ---

  TH1D* th1d_fvtxn0_clus_phi = (TH1D*)file->Get("th1d_fvtxn0_clus_phi");
  TH1D* th1d_fvtxn0_clus_phi_IR = (TH1D*)file->Get("th1d_fvtxn0_clus_phi_IR");
  TH1D* th1d_fvtxn0_clus_phi_NCIR = (TH1D*)file->Get("th1d_fvtxn0_clus_phi_NCIR");
  th1d_fvtxn0_clus_phi->Scale(1.0/nevents);
  th1d_fvtxn0_clus_phi_IR->Scale(1.0/nevents_IR);
  th1d_fvtxn0_clus_phi_NCIR->Scale(1.0/nevents_NCIR);
  th1d_fvtxn0_clus_phi->SetMinimum(0);
  th1d_fvtxn0_clus_phi->SetTitle("layer 0, FVTX North");
  th1d_fvtxn0_clus_phi->GetXaxis()->SetTitle("cluster #phi (rad)");
  th1d_fvtxn0_clus_phi->GetYaxis()->SetTitle("counts/event");
  th1d_fvtxn0_clus_phi->SetLineColor(kBlack);
  th1d_fvtxn0_clus_phi->Draw();
  th1d_fvtxn0_clus_phi_NCIR->SetLineColor(kRed);
  th1d_fvtxn0_clus_phi_NCIR->Draw("same");
  th1d_fvtxn0_clus_phi_IR->SetLineColor(kBlue);
  th1d_fvtxn0_clus_phi_IR->Draw("same");
  leg->Draw();
  c1->Print(Form("FigsPhi/radiuscut_fvtxn0_clus_phi_%d.png",handle));
  c1->Print(Form("FigsPhi/radiuscut_fvtxn0_clus_phi_%d.pdf",handle));

  // ---

  TH1D* th1d_fvtxn1_clus_phi = (TH1D*)file->Get("th1d_fvtxn1_clus_phi");
  TH1D* th1d_fvtxn1_clus_phi_IR = (TH1D*)file->Get("th1d_fvtxn1_clus_phi_IR");
  TH1D* th1d_fvtxn1_clus_phi_NCIR = (TH1D*)file->Get("th1d_fvtxn1_clus_phi_NCIR");
  th1d_fvtxn1_clus_phi->Scale(1.0/nevents);
  th1d_fvtxn1_clus_phi_IR->Scale(1.0/nevents_IR);
  th1d_fvtxn1_clus_phi_NCIR->Scale(1.0/nevents_NCIR);
  th1d_fvtxn1_clus_phi->SetMinimum(0);
  th1d_fvtxn1_clus_phi->SetTitle("layer 1, FVTX North");
  th1d_fvtxn1_clus_phi->GetXaxis()->SetTitle("cluster #phi (rad)");
  th1d_fvtxn1_clus_phi->GetYaxis()->SetTitle("counts/event");
  th1d_fvtxn1_clus_phi->SetLineColor(kBlack);
  th1d_fvtxn1_clus_phi->Draw();
  th1d_fvtxn1_clus_phi_NCIR->SetLineColor(kRed);
  th1d_fvtxn1_clus_phi_NCIR->Draw("same");
  th1d_fvtxn1_clus_phi_IR->SetLineColor(kBlue);
  th1d_fvtxn1_clus_phi_IR->Draw("same");
  leg->Draw();
  c1->Print(Form("FigsPhi/radiuscut_fvtxn1_clus_phi_%d.png",handle));
  c1->Print(Form("FigsPhi/radiuscut_fvtxn1_clus_phi_%d.pdf",handle));

  // ---

  TH1D* th1d_fvtxn2_clus_phi = (TH1D*)file->Get("th1d_fvtxn2_clus_phi");
  TH1D* th1d_fvtxn2_clus_phi_IR = (TH1D*)file->Get("th1d_fvtxn2_clus_phi_IR");
  TH1D* th1d_fvtxn2_clus_phi_NCIR = (TH1D*)file->Get("th1d_fvtxn2_clus_phi_NCIR");
  th1d_fvtxn2_clus_phi->Scale(1.0/nevents);
  th1d_fvtxn2_clus_phi_IR->Scale(1.0/nevents_IR);
  th1d_fvtxn2_clus_phi_NCIR->Scale(1.0/nevents_NCIR);
  th1d_fvtxn2_clus_phi->SetMinimum(0);
  th1d_fvtxn2_clus_phi->SetTitle("layer 2, FVTX North");
  th1d_fvtxn2_clus_phi->GetXaxis()->SetTitle("cluster #phi (rad)");
  th1d_fvtxn2_clus_phi->GetYaxis()->SetTitle("counts/event");
  th1d_fvtxn2_clus_phi->SetLineColor(kBlack);
  th1d_fvtxn2_clus_phi->Draw();
  th1d_fvtxn2_clus_phi_NCIR->SetLineColor(kRed);
  th1d_fvtxn2_clus_phi_NCIR->Draw("same");
  th1d_fvtxn2_clus_phi_IR->SetLineColor(kBlue);
  th1d_fvtxn2_clus_phi_IR->Draw("same");
  leg->Draw();
  c1->Print(Form("FigsPhi/radiuscut_fvtxn2_clus_phi_%d.png",handle));
  c1->Print(Form("FigsPhi/radiuscut_fvtxn2_clus_phi_%d.pdf",handle));

  // ---

  TH1D* th1d_fvtxn3_clus_phi = (TH1D*)file->Get("th1d_fvtxn3_clus_phi");
  TH1D* th1d_fvtxn3_clus_phi_IR = (TH1D*)file->Get("th1d_fvtxn3_clus_phi_IR");
  TH1D* th1d_fvtxn3_clus_phi_NCIR = (TH1D*)file->Get("th1d_fvtxn3_clus_phi_NCIR");
  th1d_fvtxn3_clus_phi->Scale(1.0/nevents);
  th1d_fvtxn3_clus_phi_IR->Scale(1.0/nevents_IR);
  th1d_fvtxn3_clus_phi_NCIR->Scale(1.0/nevents_NCIR);
  th1d_fvtxn3_clus_phi->SetMinimum(0);
  th1d_fvtxn3_clus_phi->SetTitle("layer 3, FVTX North");
  th1d_fvtxn3_clus_phi->GetXaxis()->SetTitle("cluster #phi (rad)");
  th1d_fvtxn3_clus_phi->GetYaxis()->SetTitle("counts/event");
  th1d_fvtxn3_clus_phi->SetLineColor(kBlack);
  th1d_fvtxn3_clus_phi->Draw();
  th1d_fvtxn3_clus_phi_NCIR->SetLineColor(kRed);
  th1d_fvtxn3_clus_phi_NCIR->Draw("same");
  th1d_fvtxn3_clus_phi_IR->SetLineColor(kBlue);
  th1d_fvtxn3_clus_phi_IR->Draw("same");
  leg->Draw();
  c1->Print(Form("FigsPhi/radiuscut_fvtxn3_clus_phi_%d.png",handle));
  c1->Print(Form("FigsPhi/radiuscut_fvtxn3_clus_phi_%d.pdf",handle));




  // ---
  // --- come back here for xy, and maybe eventually etaphi
  // ---

  TCanvas* c2 = new TCanvas("c2","",800,800);
  c2->cd();

  TH2D* th2d_fvtxs_clus_xy = (TH2D*)file->Get("th2d_fvtxs_clus_xy");
  TH2D* th2d_fvtxs_clus_xy_IR = (TH2D*)file->Get("th2d_fvtxs_clus_xy_IR");
  TH2D* th2d_fvtxs_clus_xy_OR = (TH2D*)file->Get("th2d_fvtxs_clus_xy_OR");
  th2d_fvtxs_clus_xy->Scale(1.0/nevents);
  th2d_fvtxs_clus_xy_IR->Scale(1.0/nevents_IR);
  th2d_fvtxs_clus_xy_OR->Scale(1.0/nevents_OR);
  th2d_fvtxs_clus_xy_IR->Draw("colz");
  th2d_fvtxs_clus_xy_IR->SetTitle("FVTX South, All Layers");
  th2d_fvtxs_clus_xy_IR->GetXaxis()->SetTitle("cluster x (cm)");
  th2d_fvtxs_clus_xy_IR->GetYaxis()->SetTitle("cluster y (cm)");
  TEllipse circlecutA1(0.0,0.0,5.2,5.2,180,60);
  circlecutA1.SetFillStyle(0);
  circlecutA1.SetLineWidth(2);
  TEllipse circlecutA2(0.0,0.0,9.2,9.2,119,112);
  circlecutA2.SetFillStyle(0);
  circlecutA2.SetLineWidth(2);
  TEllipse circlecutA3(0.0,0.0,5.2,5.2,200,60);
  circlecutA3.SetFillStyle(0);
  circlecutA3.SetLineWidth(2);
  TEllipse circlecutA4(0.0,0.0,5.2,5.2,20,-20);
  circlecutA4.SetFillStyle(0);
  circlecutA4.SetLineWidth(2);
  circlecutA1.Draw();
  circlecutA2.Draw();
  circlecutA3.Draw();
  circlecutA4.Draw();
  c2->SetLogz(0);
  c2->Print(Form("FigsPhi/plot2d_fvtxs_clus_xy_%d.png",handle));
  c2->Print(Form("FigsPhi/plot2d_fvtxs_clus_xy_%d.pdf",handle));
  c2->SetLogz(1);
  c2->Print(Form("FigsPhi/plot2d_fvtxs_clus_xylog_%d.png",handle));
  c2->Print(Form("FigsPhi/plot2d_fvtxs_clus_xylog_%d.pdf",handle));

  TH2D* th2d_fvtxs0_clus_xy = (TH2D*)file->Get("th2d_fvtxs0_clus_xy");
  TH2D* th2d_fvtxs0_clus_xy_IR = (TH2D*)file->Get("th2d_fvtxs0_clus_xy_IR");
  TH2D* th2d_fvtxs0_clus_xy_OR = (TH2D*)file->Get("th2d_fvtxs0_clus_xy_OR");
  th2d_fvtxs0_clus_xy->Scale(1.0/nevents);
  th2d_fvtxs0_clus_xy_IR->Scale(1.0/nevents_IR);
  th2d_fvtxs0_clus_xy_OR->Scale(1.0/nevents_OR);
  th2d_fvtxs0_clus_xy_IR->Draw("colz");

  th2d_fvtxs0_clus_xy_IR->SetTitle("FVTX South, Layer 0");
  th2d_fvtxs0_clus_xy_IR->GetXaxis()->SetTitle("cluster x (cm)");
  th2d_fvtxs0_clus_xy_IR->GetYaxis()->SetTitle("cluster y (cm)");
  circlecutA1.Draw();
  circlecutA2.Draw();
  c2->SetLogz(0);
  c2->Print(Form("FigsPhi/plot2d_fvtxs0_clus_xy_%d.png",handle));
  c2->Print(Form("FigsPhi/plot2d_fvtxs0_clus_xy_%d.pdf",handle));
  c2->SetLogz(1);
  c2->Print(Form("FigsPhi/plot2d_fvtxs0_clus_xylog_%d.png",handle));
  c2->Print(Form("FigsPhi/plot2d_fvtxs0_clus_xylog_%d.pdf",handle));

  TH2D* th2d_fvtxs1_clus_xy = (TH2D*)file->Get("th2d_fvtxs1_clus_xy");
  TH2D* th2d_fvtxs1_clus_xy_IR = (TH2D*)file->Get("th2d_fvtxs1_clus_xy_IR");
  TH2D* th2d_fvtxs1_clus_xy_OR = (TH2D*)file->Get("th2d_fvtxs1_clus_xy_OR");
  th2d_fvtxs1_clus_xy->Scale(1.0/nevents);
  th2d_fvtxs1_clus_xy_IR->Scale(1.0/nevents_IR);
  th2d_fvtxs1_clus_xy_OR->Scale(1.0/nevents_OR);
  th2d_fvtxs1_clus_xy_IR->Draw("colz");
  th2d_fvtxs1_clus_xy_IR->SetTitle("FVTX South, Layer 1");
  th2d_fvtxs1_clus_xy_IR->GetXaxis()->SetTitle("cluster x (cm)");
  th2d_fvtxs1_clus_xy_IR->GetYaxis()->SetTitle("cluster y (cm)");
  circlecutA1.Draw();
  c2->SetLogz(0);
  c2->Print(Form("FigsPhi/plot2d_fvtxs1_clus_xy_%d.png",handle));
  c2->Print(Form("FigsPhi/plot2d_fvtxs1_clus_xy_%d.pdf",handle));
  c2->SetLogz(1);
  c2->Print(Form("FigsPhi/plot2d_fvtxs1_clus_xylog_%d.png",handle));
  c2->Print(Form("FigsPhi/plot2d_fvtxs1_clus_xylog_%d.pdf",handle));

  TH2D* th2d_fvtxs2_clus_xy = (TH2D*)file->Get("th2d_fvtxs2_clus_xy");
  TH2D* th2d_fvtxs2_clus_xy_IR = (TH2D*)file->Get("th2d_fvtxs2_clus_xy_IR");
  TH2D* th2d_fvtxs2_clus_xy_OR = (TH2D*)file->Get("th2d_fvtxs2_clus_xy_OR");
  th2d_fvtxs2_clus_xy->Scale(1.0/nevents);
  th2d_fvtxs2_clus_xy_IR->Scale(1.0/nevents_IR);
  th2d_fvtxs2_clus_xy_OR->Scale(1.0/nevents_OR);
  th2d_fvtxs2_clus_xy_IR->Draw("colz");
  th2d_fvtxs2_clus_xy_IR->SetTitle("FVTX South, Layer 2");
  th2d_fvtxs2_clus_xy_IR->GetXaxis()->SetTitle("cluster x (cm)");
  th2d_fvtxs2_clus_xy_IR->GetYaxis()->SetTitle("cluster y (cm)");
  circlecutA1.Draw();
  c2->SetLogz(0);
  c2->Print(Form("FigsPhi/plot2d_fvtxs2_clus_xy_%d.png",handle));
  c2->Print(Form("FigsPhi/plot2d_fvtxs2_clus_xy_%d.pdf",handle));
  c2->SetLogz(1);
  c2->Print(Form("FigsPhi/plot2d_fvtxs2_clus_xylog_%d.png",handle));
  c2->Print(Form("FigsPhi/plot2d_fvtxs2_clus_xylog_%d.pdf",handle));

  TH2D* th2d_fvtxs3_clus_xy = (TH2D*)file->Get("th2d_fvtxs3_clus_xy");
  TH2D* th2d_fvtxs3_clus_xy_IR = (TH2D*)file->Get("th2d_fvtxs3_clus_xy_IR");
  TH2D* th2d_fvtxs3_clus_xy_OR = (TH2D*)file->Get("th2d_fvtxs3_clus_xy_OR");
  th2d_fvtxs3_clus_xy->Scale(1.0/nevents);
  th2d_fvtxs3_clus_xy_IR->Scale(1.0/nevents_IR);
  th2d_fvtxs3_clus_xy_OR->Scale(1.0/nevents_OR);
  th2d_fvtxs3_clus_xy_IR->Draw("colz");
  th2d_fvtxs3_clus_xy_IR->SetTitle("FVTX South, Layer 3");
  th2d_fvtxs3_clus_xy_IR->GetXaxis()->SetTitle("cluster x (cm)");
  th2d_fvtxs3_clus_xy_IR->GetYaxis()->SetTitle("cluster y (cm)");
  circlecutA3.Draw();
  circlecutA4.Draw();
  c2->SetLogz(0);
  c2->Print(Form("FigsPhi/plot2d_fvtxs3_clus_xy_%d.png",handle));
  c2->Print(Form("FigsPhi/plot2d_fvtxs3_clus_xy_%d.pdf",handle));
  c2->SetLogz(1);
  c2->Print(Form("FigsPhi/plot2d_fvtxs3_clus_xylog_%d.png",handle));
  c2->Print(Form("FigsPhi/plot2d_fvtxs3_clus_xylog_%d.pdf",handle));

  // ---

  TH2D* th2d_fvtxn_clus_xy = (TH2D*)file->Get("th2d_fvtxn_clus_xy");
  TH2D* th2d_fvtxn_clus_xy_IR = (TH2D*)file->Get("th2d_fvtxn_clus_xy_IR");
  TH2D* th2d_fvtxn_clus_xy_OR = (TH2D*)file->Get("th2d_fvtxn_clus_xy_OR");
  th2d_fvtxn_clus_xy->Scale(1.0/nevents);
  th2d_fvtxn_clus_xy_IR->Scale(1.0/nevents_IR);
  th2d_fvtxn_clus_xy_OR->Scale(1.0/nevents_OR);
  th2d_fvtxn_clus_xy_IR->Draw("colz");
  th2d_fvtxn_clus_xy_IR->SetTitle("FVTX North, All Layers");
  th2d_fvtxn_clus_xy_IR->GetXaxis()->SetTitle("cluster x (cm)");
  th2d_fvtxn_clus_xy_IR->GetYaxis()->SetTitle("cluster y (cm)");
  TEllipse circlecutA5(0.0,0.0,5.2,5.2,170,90);
  circlecutA5.SetFillStyle(0);
  circlecutA5.SetLineWidth(2);
  TEllipse circlecutA6(0.0,0.0,5.2,5.2,170,60);
  circlecutA6.SetFillStyle(0);
  circlecutA6.SetLineWidth(2);
  TEllipse circlecutA7(0.0,0.0,17.0,17.0,163,156);
  circlecutA7.SetFillStyle(0);
  circlecutA7.SetLineWidth(2);
  TEllipse circlecutA8(0.0,0.0,17.0,17.0,196,194);
  circlecutA8.SetFillStyle(0);
  circlecutA8.SetLineWidth(2);
  //circlecutA3.Draw();
  circlecutA5.Draw();
  circlecutA6.Draw();
  circlecutA7.Draw();
  circlecutA8.Draw();
  c2->SetLogz(0);
  c2->Print(Form("FigsPhi/plot2d_fvtxn_clus_xy_%d.png",handle));
  c2->Print(Form("FigsPhi/plot2d_fvtxn_clus_xy_%d.pdf",handle));
  c2->SetLogz(1);
  c2->Print(Form("FigsPhi/plot2d_fvtxn_clus_xylog_%d.png",handle));
  c2->Print(Form("FigsPhi/plot2d_fvtxn_clus_xylog_%d.pdf",handle));

  TH2D* th2d_fvtxn0_clus_xy = (TH2D*)file->Get("th2d_fvtxn0_clus_xy");
  TH2D* th2d_fvtxn0_clus_xy_IR = (TH2D*)file->Get("th2d_fvtxn0_clus_xy_IR");
  TH2D* th2d_fvtxn0_clus_xy_OR = (TH2D*)file->Get("th2d_fvtxn0_clus_xy_OR");
  th2d_fvtxn0_clus_xy->Scale(1.0/nevents);
  th2d_fvtxn0_clus_xy_IR->Scale(1.0/nevents_IR);
  th2d_fvtxn0_clus_xy_OR->Scale(1.0/nevents_OR);
  th2d_fvtxn0_clus_xy_IR->Draw("colz");
  th2d_fvtxn0_clus_xy_IR->SetTitle("FVTX North, Layer 0");
  th2d_fvtxn0_clus_xy_IR->GetXaxis()->SetTitle("cluster x (cm)");
  th2d_fvtxn0_clus_xy_IR->GetYaxis()->SetTitle("cluster y (cm)");
  circlecutA5.Draw();
  c2->SetLogz(0);
  c2->Print(Form("FigsPhi/plot2d_fvtxn0_clus_xy_%d.png",handle));
  c2->Print(Form("FigsPhi/plot2d_fvtxn0_clus_xy_%d.pdf",handle));
  c2->SetLogz(1);
  c2->Print(Form("FigsPhi/plot2d_fvtxn0_clus_xylog_%d.png",handle));
  c2->Print(Form("FigsPhi/plot2d_fvtxn0_clus_xylog_%d.pdf",handle));

  TH2D* th2d_fvtxn1_clus_xy = (TH2D*)file->Get("th2d_fvtxn1_clus_xy");
  TH2D* th2d_fvtxn1_clus_xy_IR = (TH2D*)file->Get("th2d_fvtxn1_clus_xy_IR");
  TH2D* th2d_fvtxn1_clus_xy_OR = (TH2D*)file->Get("th2d_fvtxn1_clus_xy_OR");
  th2d_fvtxn1_clus_xy->Scale(1.0/nevents);
  th2d_fvtxn1_clus_xy_IR->Scale(1.0/nevents_IR);
  th2d_fvtxn1_clus_xy_OR->Scale(1.0/nevents_OR);
  th2d_fvtxn1_clus_xy_IR->Draw("colz");
  th2d_fvtxn1_clus_xy_IR->SetTitle("FVTX North, Layer 1");
  th2d_fvtxn1_clus_xy_IR->GetXaxis()->SetTitle("cluster x (cm)");
  th2d_fvtxn1_clus_xy_IR->GetYaxis()->SetTitle("cluster y (cm)");
  circlecutA5.Draw();
  c2->SetLogz(0);
  c2->Print(Form("FigsPhi/plot2d_fvtxn1_clus_xy_%d.png",handle));
  c2->Print(Form("FigsPhi/plot2d_fvtxn1_clus_xy_%d.pdf",handle));
  c2->SetLogz(1);
  c2->Print(Form("FigsPhi/plot2d_fvtxn1_clus_xylog_%d.png",handle));
  c2->Print(Form("FigsPhi/plot2d_fvtxn1_clus_xylog_%d.pdf",handle));

  TH2D* th2d_fvtxn2_clus_xy = (TH2D*)file->Get("th2d_fvtxn2_clus_xy");
  TH2D* th2d_fvtxn2_clus_xy_IR = (TH2D*)file->Get("th2d_fvtxn2_clus_xy_IR");
  TH2D* th2d_fvtxn2_clus_xy_OR = (TH2D*)file->Get("th2d_fvtxn2_clus_xy_OR");
  th2d_fvtxn2_clus_xy->Scale(1.0/nevents);
  th2d_fvtxn2_clus_xy_IR->Scale(1.0/nevents_IR);
  th2d_fvtxn2_clus_xy_OR->Scale(1.0/nevents_OR);
  th2d_fvtxn2_clus_xy_IR->Draw("colz");
  th2d_fvtxn2_clus_xy_IR->SetTitle("FVTX North, Layer 2");
  th2d_fvtxn2_clus_xy_IR->GetXaxis()->SetTitle("cluster x (cm)");
  th2d_fvtxn2_clus_xy_IR->GetYaxis()->SetTitle("cluster y (cm)");
  circlecutA5.Draw();
  c2->SetLogz(0);
  c2->Print(Form("FigsPhi/plot2d_fvtxn2_clus_xy_%d.png",handle));
  c2->Print(Form("FigsPhi/plot2d_fvtxn2_clus_xy_%d.pdf",handle));
  c2->SetLogz(1);
  c2->Print(Form("FigsPhi/plot2d_fvtxn2_clus_xylog_%d.png",handle));
  c2->Print(Form("FigsPhi/plot2d_fvtxn2_clus_xylog_%d.pdf",handle));

  TH2D* th2d_fvtxn3_clus_xy = (TH2D*)file->Get("th2d_fvtxn3_clus_xy");
  TH2D* th2d_fvtxn3_clus_xy_IR = (TH2D*)file->Get("th2d_fvtxn3_clus_xy_IR");
  TH2D* th2d_fvtxn3_clus_xy_OR = (TH2D*)file->Get("th2d_fvtxn3_clus_xy_OR");
  th2d_fvtxn3_clus_xy->Scale(1.0/nevents);
  th2d_fvtxn3_clus_xy_IR->Scale(1.0/nevents_IR);
  th2d_fvtxn3_clus_xy_OR->Scale(1.0/nevents_OR);
  th2d_fvtxn3_clus_xy_IR->Draw("colz");
  th2d_fvtxn3_clus_xy_IR->SetTitle("FVTX North, Layer 3");
  th2d_fvtxn3_clus_xy_IR->GetXaxis()->SetTitle("cluster x (cm)");
  th2d_fvtxn3_clus_xy_IR->GetYaxis()->SetTitle("cluster y (cm)");
  circlecutA6.Draw();
  circlecutA7.Draw();
  circlecutA8.Draw();
  c2->SetLogz(0);
  c2->Print(Form("FigsPhi/plot2d_fvtxn3_clus_xy_%d.png",handle));
  c2->Print(Form("FigsPhi/plot2d_fvtxn3_clus_xy_%d.pdf",handle));
  c2->SetLogz(1);
  c2->Print(Form("FigsPhi/plot2d_fvtxn3_clus_xylog_%d.png",handle));
  c2->Print(Form("FigsPhi/plot2d_fvtxn3_clus_xylog_%d.pdf",handle));

  // --- come back here for periph

  TH2D* th2d_fvtxs_clus_periphxy = (TH2D*)file->Get("th2d_fvtxs_clus_periphxy");
  TH2D* th2d_fvtxs_clus_periphxy_IR = (TH2D*)file->Get("th2d_fvtxs_clus_periphxy_IR");
  TH2D* th2d_fvtxs_clus_periphxy_OR = (TH2D*)file->Get("th2d_fvtxs_clus_periphxy_OR");
  th2d_fvtxs_clus_periphxy->Scale(1.0/nevents);
  th2d_fvtxs_clus_periphxy_IR->Scale(1.0/nevents_IR);
  th2d_fvtxs_clus_periphxy_OR->Scale(1.0/nevents_OR);
  th2d_fvtxs_clus_periphxy_IR->Draw("colz");
  th2d_fvtxs_clus_periphxy_IR->SetTitle("FVTX South, All Layers");
  th2d_fvtxs_clus_periphxy_IR->GetXaxis()->SetTitle("cluster x (cm)");
  th2d_fvtxs_clus_periphxy_IR->GetYaxis()->SetTitle("cluster y (cm)");
  circlecutA1.Draw();
  circlecutA2.Draw();
  circlecutA3.Draw();
  circlecutA4.Draw();
  c2->SetLogz(0);
  c2->Print(Form("FigsPhi/plot2d_fvtxs_clus_periphxy_%d.png",handle));
  c2->Print(Form("FigsPhi/plot2d_fvtxs_clus_periphxy_%d.pdf",handle));
  c2->SetLogz(1);
  c2->Print(Form("FigsPhi/plot2d_fvtxs_clus_periphxylog_%d.png",handle));
  c2->Print(Form("FigsPhi/plot2d_fvtxs_clus_periphxylog_%d.pdf",handle));

  TH2D* th2d_fvtxs0_clus_periphxy = (TH2D*)file->Get("th2d_fvtxs0_clus_periphxy");
  TH2D* th2d_fvtxs0_clus_periphxy_IR = (TH2D*)file->Get("th2d_fvtxs0_clus_periphxy_IR");
  TH2D* th2d_fvtxs0_clus_periphxy_OR = (TH2D*)file->Get("th2d_fvtxs0_clus_periphxy_OR");
  th2d_fvtxs0_clus_periphxy->Scale(1.0/nevents);
  th2d_fvtxs0_clus_periphxy_IR->Scale(1.0/nevents_IR);
  th2d_fvtxs0_clus_periphxy_OR->Scale(1.0/nevents_OR);
  th2d_fvtxs0_clus_periphxy_IR->Draw("colz");
  th2d_fvtxs0_clus_periphxy_IR->SetTitle("FVTX South, Layer 0");
  th2d_fvtxs0_clus_periphxy_IR->GetXaxis()->SetTitle("cluster x (cm)");
  th2d_fvtxs0_clus_periphxy_IR->GetYaxis()->SetTitle("cluster y (cm)");
  circlecutA1.Draw();
  circlecutA2.Draw();
  c2->SetLogz(0);
  c2->Print(Form("FigsPhi/plot2d_fvtxs0_clus_periphxy_%d.png",handle));
  c2->Print(Form("FigsPhi/plot2d_fvtxs0_clus_periphxy_%d.pdf",handle));
  c2->SetLogz(1);
  c2->Print(Form("FigsPhi/plot2d_fvtxs0_clus_periphxylog_%d.png",handle));
  c2->Print(Form("FigsPhi/plot2d_fvtxs0_clus_periphxylog_%d.pdf",handle));

  TH2D* th2d_fvtxs1_clus_periphxy = (TH2D*)file->Get("th2d_fvtxs1_clus_periphxy");
  TH2D* th2d_fvtxs1_clus_periphxy_IR = (TH2D*)file->Get("th2d_fvtxs1_clus_periphxy_IR");
  TH2D* th2d_fvtxs1_clus_periphxy_OR = (TH2D*)file->Get("th2d_fvtxs1_clus_periphxy_OR");
  th2d_fvtxs1_clus_periphxy->Scale(1.0/nevents);
  th2d_fvtxs1_clus_periphxy_IR->Scale(1.0/nevents_IR);
  th2d_fvtxs1_clus_periphxy_OR->Scale(1.0/nevents_OR);
  th2d_fvtxs1_clus_periphxy_IR->Draw("colz");
  th2d_fvtxs1_clus_periphxy_IR->SetTitle("FVTX South, Layer 1");
  th2d_fvtxs1_clus_periphxy_IR->GetXaxis()->SetTitle("cluster x (cm)");
  th2d_fvtxs1_clus_periphxy_IR->GetYaxis()->SetTitle("cluster y (cm)");
  circlecutA1.Draw();
  c2->SetLogz(0);
  c2->Print(Form("FigsPhi/plot2d_fvtxs1_clus_periphxy_%d.png",handle));
  c2->Print(Form("FigsPhi/plot2d_fvtxs1_clus_periphxy_%d.pdf",handle));
  c2->SetLogz(1);
  c2->Print(Form("FigsPhi/plot2d_fvtxs1_clus_periphxylog_%d.png",handle));
  c2->Print(Form("FigsPhi/plot2d_fvtxs1_clus_periphxylog_%d.pdf",handle));

  TH2D* th2d_fvtxs2_clus_periphxy = (TH2D*)file->Get("th2d_fvtxs2_clus_periphxy");
  TH2D* th2d_fvtxs2_clus_periphxy_IR = (TH2D*)file->Get("th2d_fvtxs2_clus_periphxy_IR");
  TH2D* th2d_fvtxs2_clus_periphxy_OR = (TH2D*)file->Get("th2d_fvtxs2_clus_periphxy_OR");
  th2d_fvtxs2_clus_periphxy->Scale(1.0/nevents);
  th2d_fvtxs2_clus_periphxy_IR->Scale(1.0/nevents_IR);
  th2d_fvtxs2_clus_periphxy_OR->Scale(1.0/nevents_OR);
  th2d_fvtxs2_clus_periphxy_IR->Draw("colz");
  th2d_fvtxs2_clus_periphxy_IR->SetTitle("FVTX South, Layer 2");
  th2d_fvtxs2_clus_periphxy_IR->GetXaxis()->SetTitle("cluster x (cm)");
  th2d_fvtxs2_clus_periphxy_IR->GetYaxis()->SetTitle("cluster y (cm)");
  circlecutA1.Draw();
  c2->SetLogz(0);
  c2->Print(Form("FigsPhi/plot2d_fvtxs2_clus_periphxy_%d.png",handle));
  c2->Print(Form("FigsPhi/plot2d_fvtxs2_clus_periphxy_%d.pdf",handle));
  c2->SetLogz(1);
  c2->Print(Form("FigsPhi/plot2d_fvtxs2_clus_periphxylog_%d.png",handle));
  c2->Print(Form("FigsPhi/plot2d_fvtxs2_clus_periphxylog_%d.pdf",handle));

  TH2D* th2d_fvtxs3_clus_periphxy = (TH2D*)file->Get("th2d_fvtxs3_clus_periphxy");
  TH2D* th2d_fvtxs3_clus_periphxy_IR = (TH2D*)file->Get("th2d_fvtxs3_clus_periphxy_IR");
  TH2D* th2d_fvtxs3_clus_periphxy_OR = (TH2D*)file->Get("th2d_fvtxs3_clus_periphxy_OR");
  th2d_fvtxs3_clus_periphxy->Scale(1.0/nevents);
  th2d_fvtxs3_clus_periphxy_IR->Scale(1.0/nevents_IR);
  th2d_fvtxs3_clus_periphxy_OR->Scale(1.0/nevents_OR);
  th2d_fvtxs3_clus_periphxy_IR->Draw("colz");
  th2d_fvtxs3_clus_periphxy_IR->SetTitle("FVTX South, Layer 3");
  th2d_fvtxs3_clus_periphxy_IR->GetXaxis()->SetTitle("cluster x (cm)");
  th2d_fvtxs3_clus_periphxy_IR->GetYaxis()->SetTitle("cluster y (cm)");
  circlecutA3.Draw();
  circlecutA4.Draw();
  c2->SetLogz(0);
  c2->Print(Form("FigsPhi/plot2d_fvtxs3_clus_periphxy_%d.png",handle));
  c2->Print(Form("FigsPhi/plot2d_fvtxs3_clus_periphxy_%d.pdf",handle));
  c2->SetLogz(1);
  c2->Print(Form("FigsPhi/plot2d_fvtxs3_clus_periphxylog_%d.png",handle));
  c2->Print(Form("FigsPhi/plot2d_fvtxs3_clus_periphxylog_%d.pdf",handle));

  // ---

  TH2D* th2d_fvtxn_clus_periphxy = (TH2D*)file->Get("th2d_fvtxn_clus_periphxy");
  TH2D* th2d_fvtxn_clus_periphxy_IR = (TH2D*)file->Get("th2d_fvtxn_clus_periphxy_IR");
  TH2D* th2d_fvtxn_clus_periphxy_OR = (TH2D*)file->Get("th2d_fvtxn_clus_periphxy_OR");
  th2d_fvtxn_clus_periphxy->Scale(1.0/nevents);
  th2d_fvtxn_clus_periphxy_IR->Scale(1.0/nevents_IR);
  th2d_fvtxn_clus_periphxy_OR->Scale(1.0/nevents_OR);
  th2d_fvtxn_clus_periphxy_IR->Draw("colz");
  th2d_fvtxn_clus_periphxy_IR->SetTitle("FVTX North, All Layers");
  th2d_fvtxn_clus_periphxy_IR->GetXaxis()->SetTitle("cluster x (cm)");
  th2d_fvtxn_clus_periphxy_IR->GetYaxis()->SetTitle("cluster y (cm)");
  circlecutA5.Draw();
  circlecutA6.Draw();
  circlecutA7.Draw();
  circlecutA8.Draw();
  c2->SetLogz(0);
  c2->Print(Form("FigsPhi/plot2d_fvtxn_clus_periphxy_%d.png",handle));
  c2->Print(Form("FigsPhi/plot2d_fvtxn_clus_periphxy_%d.pdf",handle));
  c2->SetLogz(1);
  c2->Print(Form("FigsPhi/plot2d_fvtxn_clus_periphxylog_%d.png",handle));
  c2->Print(Form("FigsPhi/plot2d_fvtxn_clus_periphxylog_%d.pdf",handle));

  TH2D* th2d_fvtxn0_clus_periphxy = (TH2D*)file->Get("th2d_fvtxn0_clus_periphxy");
  TH2D* th2d_fvtxn0_clus_periphxy_IR = (TH2D*)file->Get("th2d_fvtxn0_clus_periphxy_IR");
  TH2D* th2d_fvtxn0_clus_periphxy_OR = (TH2D*)file->Get("th2d_fvtxn0_clus_periphxy_OR");
  th2d_fvtxn0_clus_periphxy->Scale(1.0/nevents);
  th2d_fvtxn0_clus_periphxy_IR->Scale(1.0/nevents_IR);
  th2d_fvtxn0_clus_periphxy_OR->Scale(1.0/nevents_OR);
  th2d_fvtxn0_clus_periphxy_IR->Draw("colz");
  th2d_fvtxn0_clus_periphxy_IR->SetTitle("FVTX North, Layer 0");
  th2d_fvtxn0_clus_periphxy_IR->GetXaxis()->SetTitle("cluster x (cm)");
  th2d_fvtxn0_clus_periphxy_IR->GetYaxis()->SetTitle("cluster y (cm)");
  circlecutA5.Draw();
  c2->SetLogz(0);
  c2->Print(Form("FigsPhi/plot2d_fvtxn0_clus_periphxy_%d.png",handle));
  c2->Print(Form("FigsPhi/plot2d_fvtxn0_clus_periphxy_%d.pdf",handle));
  c2->SetLogz(1);
  c2->Print(Form("FigsPhi/plot2d_fvtxn0_clus_periphxylog_%d.png",handle));
  c2->Print(Form("FigsPhi/plot2d_fvtxn0_clus_periphxylog_%d.pdf",handle));

  TH2D* th2d_fvtxn1_clus_periphxy = (TH2D*)file->Get("th2d_fvtxn1_clus_periphxy");
  TH2D* th2d_fvtxn1_clus_periphxy_IR = (TH2D*)file->Get("th2d_fvtxn1_clus_periphxy_IR");
  TH2D* th2d_fvtxn1_clus_periphxy_OR = (TH2D*)file->Get("th2d_fvtxn1_clus_periphxy_OR");
  th2d_fvtxn1_clus_periphxy->Scale(1.0/nevents);
  th2d_fvtxn1_clus_periphxy_IR->Scale(1.0/nevents_IR);
  th2d_fvtxn1_clus_periphxy_OR->Scale(1.0/nevents_OR);
  th2d_fvtxn1_clus_periphxy_IR->Draw("colz");
  th2d_fvtxn1_clus_periphxy_IR->SetTitle("FVTX North, Layer 1");
  th2d_fvtxn1_clus_periphxy_IR->GetXaxis()->SetTitle("cluster x (cm)");
  th2d_fvtxn1_clus_periphxy_IR->GetYaxis()->SetTitle("cluster y (cm)");
  circlecutA5.Draw();
  c2->SetLogz(0);
  c2->Print(Form("FigsPhi/plot2d_fvtxn1_clus_periphxy_%d.png",handle));
  c2->Print(Form("FigsPhi/plot2d_fvtxn1_clus_periphxy_%d.pdf",handle));
  c2->SetLogz(1);
  c2->Print(Form("FigsPhi/plot2d_fvtxn1_clus_periphxylog_%d.png",handle));
  c2->Print(Form("FigsPhi/plot2d_fvtxn1_clus_periphxylog_%d.pdf",handle));

  TH2D* th2d_fvtxn2_clus_periphxy = (TH2D*)file->Get("th2d_fvtxn2_clus_periphxy");
  TH2D* th2d_fvtxn2_clus_periphxy_IR = (TH2D*)file->Get("th2d_fvtxn2_clus_periphxy_IR");
  TH2D* th2d_fvtxn2_clus_periphxy_OR = (TH2D*)file->Get("th2d_fvtxn2_clus_periphxy_OR");
  th2d_fvtxn2_clus_periphxy->Scale(1.0/nevents);
  th2d_fvtxn2_clus_periphxy_IR->Scale(1.0/nevents_IR);
  th2d_fvtxn2_clus_periphxy_OR->Scale(1.0/nevents_OR);
  th2d_fvtxn2_clus_periphxy_IR->Draw("colz");
  th2d_fvtxn2_clus_periphxy_IR->SetTitle("FVTX North, Layer 2");
  th2d_fvtxn2_clus_periphxy_IR->GetXaxis()->SetTitle("cluster x (cm)");
  th2d_fvtxn2_clus_periphxy_IR->GetYaxis()->SetTitle("cluster y (cm)");
  circlecutA5.Draw();
  c2->SetLogz(0);
  c2->Print(Form("FigsPhi/plot2d_fvtxn2_clus_periphxy_%d.png",handle));
  c2->Print(Form("FigsPhi/plot2d_fvtxn2_clus_periphxy_%d.pdf",handle));
  c2->SetLogz(1);
  c2->Print(Form("FigsPhi/plot2d_fvtxn2_clus_periphxylog_%d.png",handle));
  c2->Print(Form("FigsPhi/plot2d_fvtxn2_clus_periphxylog_%d.pdf",handle));

  TH2D* th2d_fvtxn3_clus_periphxy = (TH2D*)file->Get("th2d_fvtxn3_clus_periphxy");
  TH2D* th2d_fvtxn3_clus_periphxy_IR = (TH2D*)file->Get("th2d_fvtxn3_clus_periphxy_IR");
  TH2D* th2d_fvtxn3_clus_periphxy_OR = (TH2D*)file->Get("th2d_fvtxn3_clus_periphxy_OR");
  th2d_fvtxn3_clus_periphxy->Scale(1.0/nevents);
  th2d_fvtxn3_clus_periphxy_IR->Scale(1.0/nevents_IR);
  th2d_fvtxn3_clus_periphxy_OR->Scale(1.0/nevents_OR);
  th2d_fvtxn3_clus_periphxy_IR->Draw("colz");
  th2d_fvtxn3_clus_periphxy_IR->SetTitle("FVTX North, Layer 3");
  th2d_fvtxn3_clus_periphxy_IR->GetXaxis()->SetTitle("cluster x (cm)");
  th2d_fvtxn3_clus_periphxy_IR->GetYaxis()->SetTitle("cluster y (cm)");
  circlecutA6.Draw();
  circlecutA7.Draw();
  circlecutA8.Draw();
  c2->SetLogz(0);
  c2->Print(Form("FigsPhi/plot2d_fvtxn3_clus_periphxy_%d.png",handle));
  c2->Print(Form("FigsPhi/plot2d_fvtxn3_clus_periphxy_%d.pdf",handle));
  c2->SetLogz(1);
  c2->Print(Form("FigsPhi/plot2d_fvtxn3_clus_periphxylog_%d.png",handle));
  c2->Print(Form("FigsPhi/plot2d_fvtxn3_clus_periphxylog_%d.pdf",handle));



  // --- come back here for NCIR

  TH2D* th2d_fvtxs_clus_xy_NCIR = (TH2D*)file->Get("th2d_fvtxs_clus_xy_NCIR");
  th2d_fvtxs_clus_xy_NCIR->Scale(1.0/nevents_NCIR);
  th2d_fvtxs_clus_xy_NCIR->Draw("colz");
  th2d_fvtxs_clus_xy_NCIR->SetTitle("FVTX South, All Layers");
  th2d_fvtxs_clus_xy_NCIR->GetXaxis()->SetTitle("cluster x (cm)");
  th2d_fvtxs_clus_xy_NCIR->GetYaxis()->SetTitle("cluster y (cm)");
  circlecutA1.Draw();
  circlecutA2.Draw();
  circlecutA3.Draw();
  circlecutA4.Draw();
  c2->SetLogz(0);
  c2->Print(Form("FigsPhi/plot2d_fvtxs_clus_xy_NCIR_%d.png",handle));
  c2->Print(Form("FigsPhi/plot2d_fvtxs_clus_xy_NCIR_%d.pdf",handle));
  c2->SetLogz(1);
  c2->Print(Form("FigsPhi/plot2d_fvtxs_clus_xylog_NCIR_%d.png",handle));
  c2->Print(Form("FigsPhi/plot2d_fvtxs_clus_xylog_NCIR_%d.pdf",handle));

  TH2D* th2d_fvtxs0_clus_xy_NCIR = (TH2D*)file->Get("th2d_fvtxs0_clus_xy_NCIR");
  th2d_fvtxs0_clus_xy_NCIR->Scale(1.0/nevents_NCIR);
  th2d_fvtxs0_clus_xy_NCIR->Draw("colz");
  th2d_fvtxs0_clus_xy_NCIR->SetTitle("FVTX South, Layer 0");
  th2d_fvtxs0_clus_xy_NCIR->GetXaxis()->SetTitle("cluster x (cm)");
  th2d_fvtxs0_clus_xy_NCIR->GetYaxis()->SetTitle("cluster y (cm)");
  circlecutA1.Draw();
  circlecutA2.Draw();
  c2->SetLogz(0);
  c2->Print(Form("FigsPhi/plot2d_fvtxs0_clus_xy_NCIR_%d.png",handle));
  c2->Print(Form("FigsPhi/plot2d_fvtxs0_clus_xy_NCIR_%d.pdf",handle));
  c2->SetLogz(1);
  c2->Print(Form("FigsPhi/plot2d_fvtxs0_clus_xylog_NCIR_%d.png",handle));
  c2->Print(Form("FigsPhi/plot2d_fvtxs0_clus_xylog_NCIR_%d.pdf",handle));

  TH2D* th2d_fvtxs1_clus_xy_NCIR = (TH2D*)file->Get("th2d_fvtxs1_clus_xy_NCIR");
  th2d_fvtxs1_clus_xy_NCIR->Scale(1.0/nevents_NCIR);
  th2d_fvtxs1_clus_xy_NCIR->Draw("colz");
  th2d_fvtxs1_clus_xy_NCIR->SetTitle("FVTX South, Layer 1");
  th2d_fvtxs1_clus_xy_NCIR->GetXaxis()->SetTitle("cluster x (cm)");
  th2d_fvtxs1_clus_xy_NCIR->GetYaxis()->SetTitle("cluster y (cm)");
  circlecutA1.Draw();
  c2->SetLogz(0);
  c2->Print(Form("FigsPhi/plot2d_fvtxs1_clus_xy_NCIR_%d.png",handle));
  c2->Print(Form("FigsPhi/plot2d_fvtxs1_clus_xy_NCIR_%d.pdf",handle));
  c2->SetLogz(1);
  c2->Print(Form("FigsPhi/plot2d_fvtxs1_clus_xylog_NCIR_%d.png",handle));
  c2->Print(Form("FigsPhi/plot2d_fvtxs1_clus_xylog_NCIR_%d.pdf",handle));

  TH2D* th2d_fvtxs2_clus_xy_NCIR = (TH2D*)file->Get("th2d_fvtxs2_clus_xy_NCIR");
  th2d_fvtxs2_clus_xy_NCIR->Scale(1.0/nevents_NCIR);
  th2d_fvtxs2_clus_xy_NCIR->Draw("colz");
  th2d_fvtxs2_clus_xy_NCIR->SetTitle("FVTX South, Layer 2");
  th2d_fvtxs2_clus_xy_NCIR->GetXaxis()->SetTitle("cluster x (cm)");
  th2d_fvtxs2_clus_xy_NCIR->GetYaxis()->SetTitle("cluster y (cm)");
  circlecutA1.Draw();
  c2->SetLogz(0);
  c2->Print(Form("FigsPhi/plot2d_fvtxs2_clus_xy_NCIR_%d.png",handle));
  c2->Print(Form("FigsPhi/plot2d_fvtxs2_clus_xy_NCIR_%d.pdf",handle));
  c2->SetLogz(1);
  c2->Print(Form("FigsPhi/plot2d_fvtxs2_clus_xylog_NCIR_%d.png",handle));
  c2->Print(Form("FigsPhi/plot2d_fvtxs2_clus_xylog_NCIR_%d.pdf",handle));

  TH2D* th2d_fvtxs3_clus_xy_NCIR = (TH2D*)file->Get("th2d_fvtxs3_clus_xy_NCIR");
  th2d_fvtxs3_clus_xy_NCIR->Scale(1.0/nevents_NCIR);
  th2d_fvtxs3_clus_xy_NCIR->Draw("colz");
  th2d_fvtxs3_clus_xy_NCIR->SetTitle("FVTX South, Layer 3");
  th2d_fvtxs3_clus_xy_NCIR->GetXaxis()->SetTitle("cluster x (cm)");
  th2d_fvtxs3_clus_xy_NCIR->GetYaxis()->SetTitle("cluster y (cm)");
  circlecutA3.Draw();
  circlecutA4.Draw();
  c2->SetLogz(0);
  c2->Print(Form("FigsPhi/plot2d_fvtxs3_clus_xy_NCIR_%d.png",handle));
  c2->Print(Form("FigsPhi/plot2d_fvtxs3_clus_xy_NCIR_%d.pdf",handle));
  c2->SetLogz(1);
  c2->Print(Form("FigsPhi/plot2d_fvtxs3_clus_xylog_NCIR_%d.png",handle));
  c2->Print(Form("FigsPhi/plot2d_fvtxs3_clus_xylog_NCIR_%d.pdf",handle));

  // ---

  TH2D* th2d_fvtxn_clus_xy_NCIR = (TH2D*)file->Get("th2d_fvtxn_clus_xy_NCIR");
  th2d_fvtxn_clus_xy_NCIR->Scale(1.0/nevents_NCIR);
  th2d_fvtxn_clus_xy_NCIR->Draw("colz");
  th2d_fvtxn_clus_xy_NCIR->SetTitle("FVTX North, All Layers");
  th2d_fvtxn_clus_xy_NCIR->GetXaxis()->SetTitle("cluster x (cm)");
  th2d_fvtxn_clus_xy_NCIR->GetYaxis()->SetTitle("cluster y (cm)");
  circlecutA5.Draw();
  circlecutA6.Draw();
  circlecutA7.Draw();
  circlecutA8.Draw();
  c2->SetLogz(0);
  c2->Print(Form("FigsPhi/plot2d_fvtxn_clus_xy_NCIR_%d.png",handle));
  c2->Print(Form("FigsPhi/plot2d_fvtxn_clus_xy_NCIR_%d.pdf",handle));
  c2->SetLogz(1);
  c2->Print(Form("FigsPhi/plot2d_fvtxn_clus_xylog_NCIR_%d.png",handle));
  c2->Print(Form("FigsPhi/plot2d_fvtxn_clus_xylog_NCIR_%d.pdf",handle));

  TH2D* th2d_fvtxn0_clus_xy_NCIR = (TH2D*)file->Get("th2d_fvtxn0_clus_xy_NCIR");
  th2d_fvtxn0_clus_xy_NCIR->Scale(1.0/nevents_NCIR);
  th2d_fvtxn0_clus_xy_NCIR->Draw("colz");
  th2d_fvtxn0_clus_xy_NCIR->SetTitle("FVTX North, Layer 0");
  th2d_fvtxn0_clus_xy_NCIR->GetXaxis()->SetTitle("cluster x (cm)");
  th2d_fvtxn0_clus_xy_NCIR->GetYaxis()->SetTitle("cluster y (cm)");
  circlecutA5.Draw();
  c2->SetLogz(0);
  c2->Print(Form("FigsPhi/plot2d_fvtxn0_clus_xy_NCIR_%d.png",handle));
  c2->Print(Form("FigsPhi/plot2d_fvtxn0_clus_xy_NCIR_%d.pdf",handle));
  c2->SetLogz(1);
  c2->Print(Form("FigsPhi/plot2d_fvtxn0_clus_xylog_NCIR_%d.png",handle));
  c2->Print(Form("FigsPhi/plot2d_fvtxn0_clus_xylog_NCIR_%d.pdf",handle));

  TH2D* th2d_fvtxn1_clus_xy_NCIR = (TH2D*)file->Get("th2d_fvtxn1_clus_xy_NCIR");
  th2d_fvtxn1_clus_xy_NCIR->Scale(1.0/nevents_NCIR);
  th2d_fvtxn1_clus_xy_NCIR->Draw("colz");
  th2d_fvtxn1_clus_xy_NCIR->SetTitle("FVTX North, Layer 1");
  th2d_fvtxn1_clus_xy_NCIR->GetXaxis()->SetTitle("cluster x (cm)");
  th2d_fvtxn1_clus_xy_NCIR->GetYaxis()->SetTitle("cluster y (cm)");
  circlecutA5.Draw();
  c2->SetLogz(0);
  c2->Print(Form("FigsPhi/plot2d_fvtxn1_clus_xy_NCIR_%d.png",handle));
  c2->Print(Form("FigsPhi/plot2d_fvtxn1_clus_xy_NCIR_%d.pdf",handle));
  c2->SetLogz(1);
  c2->Print(Form("FigsPhi/plot2d_fvtxn1_clus_xylog_NCIR_%d.png",handle));
  c2->Print(Form("FigsPhi/plot2d_fvtxn1_clus_xylog_NCIR_%d.pdf",handle));

  TH2D* th2d_fvtxn2_clus_xy_NCIR = (TH2D*)file->Get("th2d_fvtxn2_clus_xy_NCIR");
  th2d_fvtxn2_clus_xy_NCIR->Scale(1.0/nevents_NCIR);
  th2d_fvtxn2_clus_xy_NCIR->Draw("colz");
  th2d_fvtxn2_clus_xy_NCIR->SetTitle("FVTX North, Layer 2");
  th2d_fvtxn2_clus_xy_NCIR->GetXaxis()->SetTitle("cluster x (cm)");
  th2d_fvtxn2_clus_xy_NCIR->GetYaxis()->SetTitle("cluster y (cm)");
  circlecutA5.Draw();
  c2->SetLogz(0);
  c2->Print(Form("FigsPhi/plot2d_fvtxn2_clus_xy_NCIR_%d.png",handle));
  c2->Print(Form("FigsPhi/plot2d_fvtxn2_clus_xy_NCIR_%d.pdf",handle));
  c2->SetLogz(1);
  c2->Print(Form("FigsPhi/plot2d_fvtxn2_clus_xylog_NCIR_%d.png",handle));
  c2->Print(Form("FigsPhi/plot2d_fvtxn2_clus_xylog_NCIR_%d.pdf",handle));

  TH2D* th2d_fvtxn3_clus_xy_NCIR = (TH2D*)file->Get("th2d_fvtxn3_clus_xy_NCIR");
  th2d_fvtxn3_clus_xy_NCIR->Scale(1.0/nevents_NCIR);
  th2d_fvtxn3_clus_xy_NCIR->Draw("colz");
  th2d_fvtxn3_clus_xy_NCIR->SetTitle("FVTX North, Layer 3");
  th2d_fvtxn3_clus_xy_NCIR->GetXaxis()->SetTitle("cluster x (cm)");
  th2d_fvtxn3_clus_xy_NCIR->GetYaxis()->SetTitle("cluster y (cm)");
  circlecutA6.Draw();
  circlecutA7.Draw();
  circlecutA8.Draw();
  c2->SetLogz(0);
  c2->Print(Form("FigsPhi/plot2d_fvtxn3_clus_xy_NCIR_%d.png",handle));
  c2->Print(Form("FigsPhi/plot2d_fvtxn3_clus_xy_NCIR_%d.pdf",handle));
  c2->SetLogz(1);
  c2->Print(Form("FigsPhi/plot2d_fvtxn3_clus_xylog_NCIR_%d.png",handle));
  c2->Print(Form("FigsPhi/plot2d_fvtxn3_clus_xylog_NCIR_%d.pdf",handle));



  c1->cd();

  // ---
  // --- come back here for phieta
  // ---


  TH2D* th2d_fvtxs_clus_phieta = (TH2D*)file->Get("th2d_fvtxs_clus_phieta");
  TH2D* th2d_fvtxs_clus_phieta_IR = (TH2D*)file->Get("th2d_fvtxs_clus_phieta_IR");
  TH2D* th2d_fvtxs_clus_phieta_OR = (TH2D*)file->Get("th2d_fvtxs_clus_phieta_OR");
  th2d_fvtxs_clus_phieta->Scale(1.0/nevents);
  th2d_fvtxs_clus_phieta_IR->Scale(1.0/nevents_IR);
  th2d_fvtxs_clus_phieta_OR->Scale(1.0/nevents_OR);
  th2d_fvtxs_clus_phieta_IR->Draw("colz");
  c1->SetLogz(0);
  c1->Print(Form("FigsPhi/plot2d_fvtxs_clus_phieta_%d.png",handle));
  c1->Print(Form("FigsPhi/plot2d_fvtxs_clus_phieta_%d.pdf",handle));
  c1->SetLogz(1);
  c1->Print(Form("FigsPhi/plot2d_fvtxs_clus_phietalog_%d.png",handle));
  c1->Print(Form("FigsPhi/plot2d_fvtxs_clus_phietalog_%d.pdf",handle));

  TH2D* th2d_fvtxs0_clus_phieta = (TH2D*)file->Get("th2d_fvtxs0_clus_phieta");
  TH2D* th2d_fvtxs0_clus_phieta_IR = (TH2D*)file->Get("th2d_fvtxs0_clus_phieta_IR");
  TH2D* th2d_fvtxs0_clus_phieta_OR = (TH2D*)file->Get("th2d_fvtxs0_clus_phieta_OR");
  th2d_fvtxs0_clus_phieta->Scale(1.0/nevents);
  th2d_fvtxs0_clus_phieta_IR->Scale(1.0/nevents_IR);
  th2d_fvtxs0_clus_phieta_OR->Scale(1.0/nevents_OR);
  th2d_fvtxs0_clus_phieta_IR->Draw("colz");
  c1->SetLogz(0);
  c1->Print(Form("FigsPhi/plot2d_fvtxs0_clus_phieta_%d.png",handle));
  c1->Print(Form("FigsPhi/plot2d_fvtxs0_clus_phieta_%d.pdf",handle));
  c1->SetLogz(1);
  c1->Print(Form("FigsPhi/plot2d_fvtxs0_clus_phietalog_%d.png",handle));
  c1->Print(Form("FigsPhi/plot2d_fvtxs0_clus_phietalog_%d.pdf",handle));

  TH2D* th2d_fvtxs1_clus_phieta = (TH2D*)file->Get("th2d_fvtxs1_clus_phieta");
  TH2D* th2d_fvtxs1_clus_phieta_IR = (TH2D*)file->Get("th2d_fvtxs1_clus_phieta_IR");
  TH2D* th2d_fvtxs1_clus_phieta_OR = (TH2D*)file->Get("th2d_fvtxs1_clus_phieta_OR");
  th2d_fvtxs1_clus_phieta->Scale(1.0/nevents);
  th2d_fvtxs1_clus_phieta_IR->Scale(1.0/nevents_IR);
  th2d_fvtxs1_clus_phieta_OR->Scale(1.0/nevents_OR);
  th2d_fvtxs1_clus_phieta_IR->Draw("colz");
  c1->SetLogz(0);
  c1->Print(Form("FigsPhi/plot2d_fvtxs1_clus_phieta_%d.png",handle));
  c1->Print(Form("FigsPhi/plot2d_fvtxs1_clus_phieta_%d.pdf",handle));
  c1->SetLogz(1);
  c1->Print(Form("FigsPhi/plot2d_fvtxs1_clus_phietalog_%d.png",handle));
  c1->Print(Form("FigsPhi/plot2d_fvtxs1_clus_phietalog_%d.pdf",handle));

  TH2D* th2d_fvtxs2_clus_phieta = (TH2D*)file->Get("th2d_fvtxs2_clus_phieta");
  TH2D* th2d_fvtxs2_clus_phieta_IR = (TH2D*)file->Get("th2d_fvtxs2_clus_phieta_IR");
  TH2D* th2d_fvtxs2_clus_phieta_OR = (TH2D*)file->Get("th2d_fvtxs2_clus_phieta_OR");
  th2d_fvtxs2_clus_phieta->Scale(1.0/nevents);
  th2d_fvtxs2_clus_phieta_IR->Scale(1.0/nevents_IR);
  th2d_fvtxs2_clus_phieta_OR->Scale(1.0/nevents_OR);
  th2d_fvtxs2_clus_phieta_IR->Draw("colz");
  c1->SetLogz(0);
  c1->Print(Form("FigsPhi/plot2d_fvtxs2_clus_phieta_%d.png",handle));
  c1->Print(Form("FigsPhi/plot2d_fvtxs2_clus_phieta_%d.pdf",handle));
  c1->SetLogz(1);
  c1->Print(Form("FigsPhi/plot2d_fvtxs2_clus_phietalog_%d.png",handle));
  c1->Print(Form("FigsPhi/plot2d_fvtxs2_clus_phietalog_%d.pdf",handle));

  TH2D* th2d_fvtxs3_clus_phieta = (TH2D*)file->Get("th2d_fvtxs3_clus_phieta");
  TH2D* th2d_fvtxs3_clus_phieta_IR = (TH2D*)file->Get("th2d_fvtxs3_clus_phieta_IR");
  TH2D* th2d_fvtxs3_clus_phieta_OR = (TH2D*)file->Get("th2d_fvtxs3_clus_phieta_OR");
  th2d_fvtxs3_clus_phieta->Scale(1.0/nevents);
  th2d_fvtxs3_clus_phieta_IR->Scale(1.0/nevents_IR);
  th2d_fvtxs3_clus_phieta_OR->Scale(1.0/nevents_OR);
  th2d_fvtxs3_clus_phieta_IR->Draw("colz");
  c1->SetLogz(0);
  c1->Print(Form("FigsPhi/plot2d_fvtxs3_clus_phieta_%d.png",handle));
  c1->Print(Form("FigsPhi/plot2d_fvtxs3_clus_phieta_%d.pdf",handle));
  c1->SetLogz(1);
  c1->Print(Form("FigsPhi/plot2d_fvtxs3_clus_phietalog_%d.png",handle));
  c1->Print(Form("FigsPhi/plot2d_fvtxs3_clus_phietalog_%d.pdf",handle));

  // ---

  TH2D* th2d_fvtxn_clus_phieta = (TH2D*)file->Get("th2d_fvtxn_clus_phieta");
  TH2D* th2d_fvtxn_clus_phieta_IR = (TH2D*)file->Get("th2d_fvtxn_clus_phieta_IR");
  TH2D* th2d_fvtxn_clus_phieta_OR = (TH2D*)file->Get("th2d_fvtxn_clus_phieta_OR");
  th2d_fvtxn_clus_phieta->Scale(1.0/nevents);
  th2d_fvtxn_clus_phieta_IR->Scale(1.0/nevents_IR);
  th2d_fvtxn_clus_phieta_OR->Scale(1.0/nevents_OR);
  th2d_fvtxn_clus_phieta_IR->Draw("colz");
  c1->SetLogz(0);
  c1->Print(Form("FigsPhi/plot2d_fvtxn_clus_phieta_%d.png",handle));
  c1->Print(Form("FigsPhi/plot2d_fvtxn_clus_phieta_%d.pdf",handle));
  c1->SetLogz(1);
  c1->Print(Form("FigsPhi/plot2d_fvtxn_clus_phietalog_%d.png",handle));
  c1->Print(Form("FigsPhi/plot2d_fvtxn_clus_phietalog_%d.pdf",handle));

  TH2D* th2d_fvtxn0_clus_phieta = (TH2D*)file->Get("th2d_fvtxn0_clus_phieta");
  TH2D* th2d_fvtxn0_clus_phieta_IR = (TH2D*)file->Get("th2d_fvtxn0_clus_phieta_IR");
  TH2D* th2d_fvtxn0_clus_phieta_OR = (TH2D*)file->Get("th2d_fvtxn0_clus_phieta_OR");
  th2d_fvtxn0_clus_phieta->Scale(1.0/nevents);
  th2d_fvtxn0_clus_phieta_IR->Scale(1.0/nevents_IR);
  th2d_fvtxn0_clus_phieta_OR->Scale(1.0/nevents_OR);
  th2d_fvtxn0_clus_phieta_IR->Draw("colz");
  c1->SetLogz(0);
  c1->Print(Form("FigsPhi/plot2d_fvtxn0_clus_phieta_%d.png",handle));
  c1->Print(Form("FigsPhi/plot2d_fvtxn0_clus_phieta_%d.pdf",handle));
  c1->SetLogz(1);
  c1->Print(Form("FigsPhi/plot2d_fvtxn0_clus_phietalog_%d.png",handle));
  c1->Print(Form("FigsPhi/plot2d_fvtxn0_clus_phietalog_%d.pdf",handle));

  TH2D* th2d_fvtxn1_clus_phieta = (TH2D*)file->Get("th2d_fvtxn1_clus_phieta");
  TH2D* th2d_fvtxn1_clus_phieta_IR = (TH2D*)file->Get("th2d_fvtxn1_clus_phieta_IR");
  TH2D* th2d_fvtxn1_clus_phieta_OR = (TH2D*)file->Get("th2d_fvtxn1_clus_phieta_OR");
  th2d_fvtxn1_clus_phieta->Scale(1.0/nevents);
  th2d_fvtxn1_clus_phieta_IR->Scale(1.0/nevents_IR);
  th2d_fvtxn1_clus_phieta_OR->Scale(1.0/nevents_OR);
  th2d_fvtxn1_clus_phieta_IR->Draw("colz");
  c1->SetLogz(0);
  c1->Print(Form("FigsPhi/plot2d_fvtxn1_clus_phieta_%d.png",handle));
  c1->Print(Form("FigsPhi/plot2d_fvtxn1_clus_phieta_%d.pdf",handle));
  c1->SetLogz(1);
  c1->Print(Form("FigsPhi/plot2d_fvtxn1_clus_phietalog_%d.png",handle));
  c1->Print(Form("FigsPhi/plot2d_fvtxn1_clus_phietalog_%d.pdf",handle));

  TH2D* th2d_fvtxn2_clus_phieta = (TH2D*)file->Get("th2d_fvtxn2_clus_phieta");
  TH2D* th2d_fvtxn2_clus_phieta_IR = (TH2D*)file->Get("th2d_fvtxn2_clus_phieta_IR");
  TH2D* th2d_fvtxn2_clus_phieta_OR = (TH2D*)file->Get("th2d_fvtxn2_clus_phieta_OR");
  th2d_fvtxn2_clus_phieta->Scale(1.0/nevents);
  th2d_fvtxn2_clus_phieta_IR->Scale(1.0/nevents_IR);
  th2d_fvtxn2_clus_phieta_OR->Scale(1.0/nevents_OR);
  th2d_fvtxn2_clus_phieta_IR->Draw("colz");
  c1->SetLogz(0);
  c1->Print(Form("FigsPhi/plot2d_fvtxn2_clus_phieta_%d.png",handle));
  c1->Print(Form("FigsPhi/plot2d_fvtxn2_clus_phieta_%d.pdf",handle));
  c1->SetLogz(1);
  c1->Print(Form("FigsPhi/plot2d_fvtxn2_clus_phietalog_%d.png",handle));
  c1->Print(Form("FigsPhi/plot2d_fvtxn2_clus_phietalog_%d.pdf",handle));

  TH2D* th2d_fvtxn3_clus_phieta = (TH2D*)file->Get("th2d_fvtxn3_clus_phieta");
  TH2D* th2d_fvtxn3_clus_phieta_IR = (TH2D*)file->Get("th2d_fvtxn3_clus_phieta_IR");
  TH2D* th2d_fvtxn3_clus_phieta_OR = (TH2D*)file->Get("th2d_fvtxn3_clus_phieta_OR");
  th2d_fvtxn3_clus_phieta->Scale(1.0/nevents);
  th2d_fvtxn3_clus_phieta_IR->Scale(1.0/nevents_IR);
  th2d_fvtxn3_clus_phieta_OR->Scale(1.0/nevents_OR);
  th2d_fvtxn3_clus_phieta_IR->Draw("colz");
  c1->SetLogz(0);
  c1->Print(Form("FigsPhi/plot2d_fvtxn3_clus_phieta_%d.png",handle));
  c1->Print(Form("FigsPhi/plot2d_fvtxn3_clus_phieta_%d.pdf",handle));
  c1->SetLogz(1);
  c1->Print(Form("FigsPhi/plot2d_fvtxn3_clus_phietalog_%d.png",handle));
  c1->Print(Form("FigsPhi/plot2d_fvtxn3_clus_phietalog_%d.pdf",handle));



  delete c1;
  delete c2;

}




// COME BACK HERE


void compbg()
{

  TCanvas* c1 = new TCanvas("c1","");

  TFile* fileB = TFile::Open(Form("RootFiles/sum20_bad.root"));
  TFile* fileG = TFile::Open(Form("RootFiles/sum20_better.root"));

  TH1D* th1dG_counter = (TH1D*)fileG->Get("th1d_FVTX_nclus");
  TH1D* th1dG_counter_IR = (TH1D*)fileG->Get("th1d_FVTX_nclus_IR");
  TH1D* th1dG_counter_OR = (TH1D*)fileG->Get("th1d_FVTX_nclus_OR");
  int neventsG = th1dG_counter->GetEntries();
  int neventsG_IR = th1dG_counter_IR->GetEntries();
  int neventsG_OR = th1dG_counter_OR->GetEntries();

  TH1D* th1dB_counter = (TH1D*)fileB->Get("th1d_FVTX_nclus");
  TH1D* th1dB_counter_IR = (TH1D*)fileB->Get("th1d_FVTX_nclus_IR");
  TH1D* th1dB_counter_OR = (TH1D*)fileB->Get("th1d_FVTX_nclus_OR");
  int neventsB = th1dB_counter->GetEntries();
  int neventsB_IR = th1dB_counter_IR->GetEntries();
  int neventsB_OR = th1dB_counter_OR->GetEntries();

  // ---
  // ---
  // ---

  // ---
  TH1D* th1dG_fvtxs_clus_phi = (TH1D*)fileG->Get("th1d_fvtxs_clus_phi");
  TH1D* th1dG_fvtxs_clus_phi_IR = (TH1D*)fileG->Get("th1d_fvtxs_clus_phi_IR");
  TH1D* th1dG_fvtxs_clus_phi_OR = (TH1D*)fileG->Get("th1d_fvtxs_clus_phi_OR");
  th1dG_fvtxs_clus_phi->Scale(1.0/neventsG);
  th1dG_fvtxs_clus_phi_IR->Scale(1.0/neventsG_IR);
  th1dG_fvtxs_clus_phi_OR->Scale(1.0/neventsG_OR);

  TH1D* th1dB_fvtxs_clus_phi = (TH1D*)fileB->Get("th1d_fvtxs_clus_phi");
  TH1D* th1dB_fvtxs_clus_phi_IR = (TH1D*)fileB->Get("th1d_fvtxs_clus_phi_IR");
  TH1D* th1dB_fvtxs_clus_phi_OR = (TH1D*)fileB->Get("th1d_fvtxs_clus_phi_OR");
  th1dB_fvtxs_clus_phi->Scale(1.0/neventsB);
  th1dB_fvtxs_clus_phi_IR->Scale(1.0/neventsB_IR);
  th1dB_fvtxs_clus_phi_OR->Scale(1.0/neventsB_OR);

  // ---
  th1dB_fvtxs_clus_phi->Draw();
  th1dB_fvtxs_clus_phi->SetMinimum(0);
  th1dB_fvtxs_clus_phi->SetTitle(Form("d+Au 20 GeV, all x,y, FVTX South"));
  th1dB_fvtxs_clus_phi->GetXaxis()->SetTitle("cluster #phi (rad)");
  th1dB_fvtxs_clus_phi->GetYaxis()->SetTitle("counts/event");
  th1dG_fvtxs_clus_phi->Draw("same");
  th1dG_fvtxs_clus_phi->SetLineColor(kRed);
  TLegend* leg = new TLegend(0.18,0.68,0.38,0.88);
  leg->AddEntry(th1dB_fvtxs_clus_phi,"all runs","l");
  leg->AddEntry(th1dG_fvtxs_clus_phi,"good runs","l");
  leg->SetTextSize(0.045);
  leg->Draw();
  c1->Print(Form("FigsPhi/compbg_fvtxs_clus_phi_20.png"));
  c1->Print(Form("FigsPhi/compbg_fvtxs_clus_phi_20.pdf"));

  th1dB_fvtxs_clus_phi_IR->Draw();
  th1dB_fvtxs_clus_phi_IR->SetMinimum(0);
  th1dB_fvtxs_clus_phi_IR->SetTitle(Form("d+Au 20 GeV, x,y within 0.15 cm of beam center, FVTX South"));
  th1dB_fvtxs_clus_phi_IR->GetXaxis()->SetTitle("cluster #phi (rad)");
  th1dB_fvtxs_clus_phi_IR->GetYaxis()->SetTitle("counts/event");
  th1dG_fvtxs_clus_phi_IR->Draw("same");
  th1dG_fvtxs_clus_phi_IR->SetLineColor(kRed);
  leg->Draw();
  c1->Print(Form("FigsPhi/compbg_fvtxs_clus_phi_20_IR.png"));
  c1->Print(Form("FigsPhi/compbg_fvtxs_clus_phi_20_IR.pdf"));

  th1dB_fvtxs_clus_phi->Draw();
  th1dB_fvtxs_clus_phi->SetMinimum(0);
  th1dB_fvtxs_clus_phi->SetTitle(Form("d+Au 20 GeV, FVTX South"));
  th1dB_fvtxs_clus_phi->GetXaxis()->SetTitle("cluster #phi (rad)");
  th1dB_fvtxs_clus_phi->GetYaxis()->SetTitle("counts/event");
  th1dG_fvtxs_clus_phi_IR->Draw("same");
  th1dG_fvtxs_clus_phi_IR->SetLineColor(kRed);
  delete leg;
  leg = new TLegend(0.18,0.68,0.38,0.88);
  leg->AddEntry(th1dB_fvtxs_clus_phi,"all runs, no radius cut","l");
  leg->AddEntry(th1dG_fvtxs_clus_phi_IR,"good runs, with radius cut","l");
  leg->SetTextSize(0.045);
  leg->Draw();
  c1->Print(Form("FigsPhi/compbg_fvtxs_clus_phi_20_BGIR.png"));
  c1->Print(Form("FigsPhi/compbg_fvtxs_clus_phi_20_BGIR.pdf"));



  // ---
  // --- fvtxs0
  // ---

  // ---
  TH1D* th1dG_fvtxs0_clus_phi = (TH1D*)fileG->Get("th1d_fvtxs0_clus_phi");
  TH1D* th1dG_fvtxs0_clus_phi_IR = (TH1D*)fileG->Get("th1d_fvtxs0_clus_phi_IR");
  TH1D* th1dG_fvtxs0_clus_phi_OR = (TH1D*)fileG->Get("th1d_fvtxs0_clus_phi_OR");
  th1dG_fvtxs0_clus_phi->Scale(1.0/neventsG);
  th1dG_fvtxs0_clus_phi_IR->Scale(1.0/neventsG_IR);
  th1dG_fvtxs0_clus_phi_OR->Scale(1.0/neventsG_OR);

  TH1D* th1dB_fvtxs0_clus_phi = (TH1D*)fileB->Get("th1d_fvtxs0_clus_phi");
  TH1D* th1dB_fvtxs0_clus_phi_IR = (TH1D*)fileB->Get("th1d_fvtxs0_clus_phi_IR");
  TH1D* th1dB_fvtxs0_clus_phi_OR = (TH1D*)fileB->Get("th1d_fvtxs0_clus_phi_OR");
  th1dB_fvtxs0_clus_phi->Scale(1.0/neventsB);
  th1dB_fvtxs0_clus_phi_IR->Scale(1.0/neventsB_IR);
  th1dB_fvtxs0_clus_phi_OR->Scale(1.0/neventsB_OR);

  // ---
  th1dB_fvtxs0_clus_phi->Draw();
  th1dB_fvtxs0_clus_phi->SetMinimum(0);
  th1dB_fvtxs0_clus_phi->SetTitle(Form("d+Au 20 GeV, all x,y, FVTX South"));
  th1dB_fvtxs0_clus_phi->GetXaxis()->SetTitle("cluster #phi (rad)");
  th1dB_fvtxs0_clus_phi->GetYaxis()->SetTitle("counts/event");
  th1dG_fvtxs0_clus_phi->Draw("same");
  th1dG_fvtxs0_clus_phi->SetLineColor(kRed);
  delete leg;
  leg = new TLegend(0.18,0.68,0.38,0.88);
  leg->AddEntry(th1dB_fvtxs0_clus_phi,"all runs","l");
  leg->AddEntry(th1dG_fvtxs0_clus_phi,"good runs","l");
  leg->SetTextSize(0.045);
  leg->Draw();
  c1->Print(Form("FigsPhi/compbg_fvtxs0_clus_phi_20.png"));
  c1->Print(Form("FigsPhi/compbg_fvtxs0_clus_phi_20.pdf"));

  th1dB_fvtxs0_clus_phi_IR->Draw();
  th1dB_fvtxs0_clus_phi_IR->SetMinimum(0);
  th1dB_fvtxs0_clus_phi_IR->SetTitle(Form("d+Au 20 GeV, x,y within 0.15 cm of beam center, FVTX South"));
  th1dB_fvtxs0_clus_phi_IR->GetXaxis()->SetTitle("cluster #phi (rad)");
  th1dB_fvtxs0_clus_phi_IR->GetYaxis()->SetTitle("counts/event");
  th1dG_fvtxs0_clus_phi_IR->Draw("same");
  th1dG_fvtxs0_clus_phi_IR->SetLineColor(kRed);
  leg->Draw();
  c1->Print(Form("FigsPhi/compbg_fvtxs0_clus_phi_20_IR.png"));
  c1->Print(Form("FigsPhi/compbg_fvtxs0_clus_phi_20_IR.pdf"));

  th1dB_fvtxs0_clus_phi->Draw();
  th1dB_fvtxs0_clus_phi->SetMinimum(0);
  th1dB_fvtxs0_clus_phi->SetTitle(Form("d+Au 20 GeV, FVTX South"));
  th1dB_fvtxs0_clus_phi->GetXaxis()->SetTitle("cluster #phi (rad)");
  th1dB_fvtxs0_clus_phi->GetYaxis()->SetTitle("counts/event");
  th1dG_fvtxs0_clus_phi_IR->Draw("same");
  th1dG_fvtxs0_clus_phi_IR->SetLineColor(kRed);
  delete leg;
  leg = new TLegend(0.18,0.68,0.38,0.88);
  leg->AddEntry(th1dB_fvtxs0_clus_phi,"all runs, no radius cut","l");
  leg->AddEntry(th1dG_fvtxs0_clus_phi_IR,"good runs, with radius cut","l");
  leg->SetTextSize(0.045);
  leg->Draw();
  c1->Print(Form("FigsPhi/compbg_fvtxs0_clus_phi_20_BGIR.png"));
  c1->Print(Form("FigsPhi/compbg_fvtxs0_clus_phi_20_BGIR.pdf"));



  // ---
  // --- fvtxs1
  // ---

  // ---
  TH1D* th1dG_fvtxs1_clus_phi = (TH1D*)fileG->Get("th1d_fvtxs1_clus_phi");
  TH1D* th1dG_fvtxs1_clus_phi_IR = (TH1D*)fileG->Get("th1d_fvtxs1_clus_phi_IR");
  TH1D* th1dG_fvtxs1_clus_phi_OR = (TH1D*)fileG->Get("th1d_fvtxs1_clus_phi_OR");
  th1dG_fvtxs1_clus_phi->Scale(1.0/neventsG);
  th1dG_fvtxs1_clus_phi_IR->Scale(1.0/neventsG_IR);
  th1dG_fvtxs1_clus_phi_OR->Scale(1.0/neventsG_OR);

  TH1D* th1dB_fvtxs1_clus_phi = (TH1D*)fileB->Get("th1d_fvtxs1_clus_phi");
  TH1D* th1dB_fvtxs1_clus_phi_IR = (TH1D*)fileB->Get("th1d_fvtxs1_clus_phi_IR");
  TH1D* th1dB_fvtxs1_clus_phi_OR = (TH1D*)fileB->Get("th1d_fvtxs1_clus_phi_OR");
  th1dB_fvtxs1_clus_phi->Scale(1.0/neventsB);
  th1dB_fvtxs1_clus_phi_IR->Scale(1.0/neventsB_IR);
  th1dB_fvtxs1_clus_phi_OR->Scale(1.0/neventsB_OR);

  // ---
  th1dB_fvtxs1_clus_phi->Draw();
  th1dB_fvtxs1_clus_phi->SetMinimum(0);
  th1dB_fvtxs1_clus_phi->SetTitle(Form("d+Au 20 GeV, all x,y, FVTX South"));
  th1dB_fvtxs1_clus_phi->GetXaxis()->SetTitle("cluster #phi (rad)");
  th1dB_fvtxs1_clus_phi->GetYaxis()->SetTitle("counts/event");
  th1dG_fvtxs1_clus_phi->Draw("same");
  th1dG_fvtxs1_clus_phi->SetLineColor(kRed);
  delete leg;
  leg = new TLegend(0.18,0.68,0.38,0.88);
  leg->AddEntry(th1dB_fvtxs1_clus_phi,"all runs","l");
  leg->AddEntry(th1dG_fvtxs1_clus_phi,"good runs","l");
  leg->SetTextSize(0.045);
  leg->Draw();
  c1->Print(Form("FigsPhi/compbg_fvtxs1_clus_phi_20.png"));
  c1->Print(Form("FigsPhi/compbg_fvtxs1_clus_phi_20.pdf"));

  th1dB_fvtxs1_clus_phi_IR->Draw();
  th1dB_fvtxs1_clus_phi_IR->SetMinimum(0);
  th1dB_fvtxs1_clus_phi_IR->SetTitle(Form("d+Au 20 GeV, x,y within 0.15 cm of beam center, FVTX South"));
  th1dB_fvtxs1_clus_phi_IR->GetXaxis()->SetTitle("cluster #phi (rad)");
  th1dB_fvtxs1_clus_phi_IR->GetYaxis()->SetTitle("counts/event");
  th1dG_fvtxs1_clus_phi_IR->Draw("same");
  th1dG_fvtxs1_clus_phi_IR->SetLineColor(kRed);
  leg->Draw();
  c1->Print(Form("FigsPhi/compbg_fvtxs1_clus_phi_20_IR.png"));
  c1->Print(Form("FigsPhi/compbg_fvtxs1_clus_phi_20_IR.pdf"));

  th1dB_fvtxs1_clus_phi->Draw();
  th1dB_fvtxs1_clus_phi->SetMinimum(0);
  th1dB_fvtxs1_clus_phi->SetTitle(Form("d+Au 20 GeV, FVTX South"));
  th1dB_fvtxs1_clus_phi->GetXaxis()->SetTitle("cluster #phi (rad)");
  th1dB_fvtxs1_clus_phi->GetYaxis()->SetTitle("counts/event");
  th1dG_fvtxs1_clus_phi_IR->Draw("same");
  th1dG_fvtxs1_clus_phi_IR->SetLineColor(kRed);
  delete leg;
  leg = new TLegend(0.18,0.68,0.38,0.88);
  leg->AddEntry(th1dB_fvtxs1_clus_phi,"all runs, no radius cut","l");
  leg->AddEntry(th1dG_fvtxs1_clus_phi_IR,"good runs, with radius cut","l");
  leg->SetTextSize(0.045);
  leg->Draw();
  c1->Print(Form("FigsPhi/compbg_fvtxs1_clus_phi_20_BGIR.png"));
  c1->Print(Form("FigsPhi/compbg_fvtxs1_clus_phi_20_BGIR.pdf"));



  // ---
  // --- fvtxs2
  // ---

  // ---
  TH1D* th1dG_fvtxs2_clus_phi = (TH1D*)fileG->Get("th1d_fvtxs2_clus_phi");
  TH1D* th1dG_fvtxs2_clus_phi_IR = (TH1D*)fileG->Get("th1d_fvtxs2_clus_phi_IR");
  TH1D* th1dG_fvtxs2_clus_phi_OR = (TH1D*)fileG->Get("th1d_fvtxs2_clus_phi_OR");
  th1dG_fvtxs2_clus_phi->Scale(1.0/neventsG);
  th1dG_fvtxs2_clus_phi_IR->Scale(1.0/neventsG_IR);
  th1dG_fvtxs2_clus_phi_OR->Scale(1.0/neventsG_OR);

  TH1D* th1dB_fvtxs2_clus_phi = (TH1D*)fileB->Get("th1d_fvtxs2_clus_phi");
  TH1D* th1dB_fvtxs2_clus_phi_IR = (TH1D*)fileB->Get("th1d_fvtxs2_clus_phi_IR");
  TH1D* th1dB_fvtxs2_clus_phi_OR = (TH1D*)fileB->Get("th1d_fvtxs2_clus_phi_OR");
  th1dB_fvtxs2_clus_phi->Scale(1.0/neventsB);
  th1dB_fvtxs2_clus_phi_IR->Scale(1.0/neventsB_IR);
  th1dB_fvtxs2_clus_phi_OR->Scale(1.0/neventsB_OR);

  // ---
  th1dB_fvtxs2_clus_phi->Draw();
  th1dB_fvtxs2_clus_phi->SetMinimum(0);
  th1dB_fvtxs2_clus_phi->SetTitle(Form("d+Au 20 GeV, all x,y, FVTX South"));
  th1dB_fvtxs2_clus_phi->GetXaxis()->SetTitle("cluster #phi (rad)");
  th1dB_fvtxs2_clus_phi->GetYaxis()->SetTitle("counts/event");
  th1dG_fvtxs2_clus_phi->Draw("same");
  th1dG_fvtxs2_clus_phi->SetLineColor(kRed);
  delete leg;
  leg = new TLegend(0.18,0.68,0.38,0.88);
  leg->AddEntry(th1dB_fvtxs2_clus_phi,"all runs","l");
  leg->AddEntry(th1dG_fvtxs2_clus_phi,"good runs","l");
  leg->SetTextSize(0.045);
  leg->Draw();
  c1->Print(Form("FigsPhi/compbg_fvtxs2_clus_phi_20.png"));
  c1->Print(Form("FigsPhi/compbg_fvtxs2_clus_phi_20.pdf"));

  th1dB_fvtxs2_clus_phi_IR->Draw();
  th1dB_fvtxs2_clus_phi_IR->SetMinimum(0);
  th1dB_fvtxs2_clus_phi_IR->SetTitle(Form("d+Au 20 GeV, x,y within 0.15 cm of beam center, FVTX South"));
  th1dB_fvtxs2_clus_phi_IR->GetXaxis()->SetTitle("cluster #phi (rad)");
  th1dB_fvtxs2_clus_phi_IR->GetYaxis()->SetTitle("counts/event");
  th1dG_fvtxs2_clus_phi_IR->Draw("same");
  th1dG_fvtxs2_clus_phi_IR->SetLineColor(kRed);
  leg->Draw();
  c1->Print(Form("FigsPhi/compbg_fvtxs2_clus_phi_20_IR.png"));
  c1->Print(Form("FigsPhi/compbg_fvtxs2_clus_phi_20_IR.pdf"));

  th1dB_fvtxs2_clus_phi->Draw();
  th1dB_fvtxs2_clus_phi->SetMinimum(0);
  th1dB_fvtxs2_clus_phi->SetTitle(Form("d+Au 20 GeV, FVTX South"));
  th1dB_fvtxs2_clus_phi->GetXaxis()->SetTitle("cluster #phi (rad)");
  th1dB_fvtxs2_clus_phi->GetYaxis()->SetTitle("counts/event");
  th1dG_fvtxs2_clus_phi_IR->Draw("same");
  th1dG_fvtxs2_clus_phi_IR->SetLineColor(kRed);
  delete leg;
  leg = new TLegend(0.18,0.68,0.38,0.88);
  leg->AddEntry(th1dB_fvtxs2_clus_phi,"all runs, no radius cut","l");
  leg->AddEntry(th1dG_fvtxs2_clus_phi_IR,"good runs, with radius cut","l");
  leg->SetTextSize(0.045);
  leg->Draw();
  c1->Print(Form("FigsPhi/compbg_fvtxs2_clus_phi_20_BGIR.png"));
  c1->Print(Form("FigsPhi/compbg_fvtxs2_clus_phi_20_BGIR.pdf"));



  // ---
  // --- fvtxs3
  // ---

  // ---
  TH1D* th1dG_fvtxs3_clus_phi = (TH1D*)fileG->Get("th1d_fvtxs3_clus_phi");
  TH1D* th1dG_fvtxs3_clus_phi_IR = (TH1D*)fileG->Get("th1d_fvtxs3_clus_phi_IR");
  TH1D* th1dG_fvtxs3_clus_phi_OR = (TH1D*)fileG->Get("th1d_fvtxs3_clus_phi_OR");
  th1dG_fvtxs3_clus_phi->Scale(1.0/neventsG);
  th1dG_fvtxs3_clus_phi_IR->Scale(1.0/neventsG_IR);
  th1dG_fvtxs3_clus_phi_OR->Scale(1.0/neventsG_OR);

  TH1D* th1dB_fvtxs3_clus_phi = (TH1D*)fileB->Get("th1d_fvtxs3_clus_phi");
  TH1D* th1dB_fvtxs3_clus_phi_IR = (TH1D*)fileB->Get("th1d_fvtxs3_clus_phi_IR");
  TH1D* th1dB_fvtxs3_clus_phi_OR = (TH1D*)fileB->Get("th1d_fvtxs3_clus_phi_OR");
  th1dB_fvtxs3_clus_phi->Scale(1.0/neventsB);
  th1dB_fvtxs3_clus_phi_IR->Scale(1.0/neventsB_IR);
  th1dB_fvtxs3_clus_phi_OR->Scale(1.0/neventsB_OR);

  // ---
  th1dB_fvtxs3_clus_phi->Draw();
  th1dB_fvtxs3_clus_phi->SetMinimum(0);
  th1dB_fvtxs3_clus_phi->SetTitle(Form("d+Au 20 GeV, all x,y, FVTX South"));
  th1dB_fvtxs3_clus_phi->GetXaxis()->SetTitle("cluster #phi (rad)");
  th1dB_fvtxs3_clus_phi->GetYaxis()->SetTitle("counts/event");
  th1dG_fvtxs3_clus_phi->Draw("same");
  th1dG_fvtxs3_clus_phi->SetLineColor(kRed);
  delete leg;
  leg = new TLegend(0.18,0.68,0.38,0.88);
  leg->AddEntry(th1dB_fvtxs3_clus_phi,"all runs","l");
  leg->AddEntry(th1dG_fvtxs3_clus_phi,"good runs","l");
  leg->SetTextSize(0.045);
  leg->Draw();
  c1->Print(Form("FigsPhi/compbg_fvtxs3_clus_phi_20.png"));
  c1->Print(Form("FigsPhi/compbg_fvtxs3_clus_phi_20.pdf"));

  th1dB_fvtxs3_clus_phi_IR->Draw();
  th1dB_fvtxs3_clus_phi_IR->SetMinimum(0);
  th1dB_fvtxs3_clus_phi_IR->SetTitle(Form("d+Au 20 GeV, x,y within 0.15 cm of beam center, FVTX South"));
  th1dB_fvtxs3_clus_phi_IR->GetXaxis()->SetTitle("cluster #phi (rad)");
  th1dB_fvtxs3_clus_phi_IR->GetYaxis()->SetTitle("counts/event");
  th1dG_fvtxs3_clus_phi_IR->Draw("same");
  th1dG_fvtxs3_clus_phi_IR->SetLineColor(kRed);
  leg->Draw();
  c1->Print(Form("FigsPhi/compbg_fvtxs3_clus_phi_20_IR.png"));
  c1->Print(Form("FigsPhi/compbg_fvtxs3_clus_phi_20_IR.pdf"));

  th1dB_fvtxs3_clus_phi->Draw();
  th1dB_fvtxs3_clus_phi->SetMinimum(0);
  th1dB_fvtxs3_clus_phi->SetTitle(Form("d+Au 20 GeV, FVTX South"));
  th1dB_fvtxs3_clus_phi->GetXaxis()->SetTitle("cluster #phi (rad)");
  th1dB_fvtxs3_clus_phi->GetYaxis()->SetTitle("counts/event");
  th1dG_fvtxs3_clus_phi_IR->Draw("same");
  th1dG_fvtxs3_clus_phi_IR->SetLineColor(kRed);
  delete leg;
  leg = new TLegend(0.18,0.68,0.38,0.88);
  leg->AddEntry(th1dB_fvtxs3_clus_phi,"all runs, no radius cut","l");
  leg->AddEntry(th1dG_fvtxs3_clus_phi_IR,"good runs, with radius cut","l");
  leg->SetTextSize(0.045);
  leg->Draw();
  c1->Print(Form("FigsPhi/compbg_fvtxs3_clus_phi_20_BGIR.png"));
  c1->Print(Form("FigsPhi/compbg_fvtxs3_clus_phi_20_BGIR.pdf"));



  // ---
  // --- come back here for north

  // ---
  // ---
  // ---

  // ---
  TH1D* th1dG_fvtxn_clus_phi = (TH1D*)fileG->Get("th1d_fvtxn_clus_phi");
  TH1D* th1dG_fvtxn_clus_phi_IR = (TH1D*)fileG->Get("th1d_fvtxn_clus_phi_IR");
  TH1D* th1dG_fvtxn_clus_phi_OR = (TH1D*)fileG->Get("th1d_fvtxn_clus_phi_OR");
  th1dG_fvtxn_clus_phi->Scale(1.0/neventsG);
  th1dG_fvtxn_clus_phi_IR->Scale(1.0/neventsG_IR);
  th1dG_fvtxn_clus_phi_OR->Scale(1.0/neventsG_OR);

  TH1D* th1dB_fvtxn_clus_phi = (TH1D*)fileB->Get("th1d_fvtxn_clus_phi");
  TH1D* th1dB_fvtxn_clus_phi_IR = (TH1D*)fileB->Get("th1d_fvtxn_clus_phi_IR");
  TH1D* th1dB_fvtxn_clus_phi_OR = (TH1D*)fileB->Get("th1d_fvtxn_clus_phi_OR");
  th1dB_fvtxn_clus_phi->Scale(1.0/neventsB);
  th1dB_fvtxn_clus_phi_IR->Scale(1.0/neventsB_IR);
  th1dB_fvtxn_clus_phi_OR->Scale(1.0/neventsB_OR);

  // ---
  th1dB_fvtxn_clus_phi->Draw();
  th1dB_fvtxn_clus_phi->SetMinimum(0);
  th1dB_fvtxn_clus_phi->SetTitle(Form("d+Au 20 GeV, all x,y, FVTX North"));
  th1dB_fvtxn_clus_phi->GetXaxis()->SetTitle("cluster #phi (rad)");
  th1dB_fvtxn_clus_phi->GetYaxis()->SetTitle("counts/event");
  th1dG_fvtxn_clus_phi->Draw("same");
  th1dG_fvtxn_clus_phi->SetLineColor(kRed);
  TLegend* leg = new TLegend(0.18,0.68,0.38,0.88);
  leg->AddEntry(th1dB_fvtxn_clus_phi,"all runs","l");
  leg->AddEntry(th1dG_fvtxn_clus_phi,"good runs","l");
  leg->SetTextSize(0.045);
  leg->Draw();
  c1->Print(Form("FigsPhi/compbg_fvtxn_clus_phi_20.png"));
  c1->Print(Form("FigsPhi/compbg_fvtxn_clus_phi_20.pdf"));

  th1dB_fvtxn_clus_phi_IR->Draw();
  th1dB_fvtxn_clus_phi_IR->SetMinimum(0);
  th1dB_fvtxn_clus_phi_IR->SetTitle(Form("d+Au 20 GeV, x,y within 0.15 cm of beam center, FVTX North"));
  th1dB_fvtxn_clus_phi_IR->GetXaxis()->SetTitle("cluster #phi (rad)");
  th1dB_fvtxn_clus_phi_IR->GetYaxis()->SetTitle("counts/event");
  th1dG_fvtxn_clus_phi_IR->Draw("same");
  th1dG_fvtxn_clus_phi_IR->SetLineColor(kRed);
  leg->Draw();
  c1->Print(Form("FigsPhi/compbg_fvtxn_clus_phi_20_IR.png"));
  c1->Print(Form("FigsPhi/compbg_fvtxn_clus_phi_20_IR.pdf"));

  th1dB_fvtxn_clus_phi->Draw();
  th1dB_fvtxn_clus_phi->SetMinimum(0);
  th1dB_fvtxn_clus_phi->SetTitle(Form("d+Au 20 GeV, FVTX North"));
  th1dB_fvtxn_clus_phi->GetXaxis()->SetTitle("cluster #phi (rad)");
  th1dB_fvtxn_clus_phi->GetYaxis()->SetTitle("counts/event");
  th1dG_fvtxn_clus_phi_IR->Draw("same");
  th1dG_fvtxn_clus_phi_IR->SetLineColor(kRed);
  delete leg;
  leg = new TLegend(0.18,0.68,0.38,0.88);
  leg->AddEntry(th1dB_fvtxn_clus_phi,"all runs, no radius cut","l");
  leg->AddEntry(th1dG_fvtxn_clus_phi_IR,"good runs, with radius cut","l");
  leg->SetTextSize(0.045);
  leg->Draw();
  c1->Print(Form("FigsPhi/compbg_fvtxn_clus_phi_20_BGIR.png"));
  c1->Print(Form("FigsPhi/compbg_fvtxn_clus_phi_20_BGIR.pdf"));



  // ---
  // --- fvtxn0
  // ---

  // ---
  TH1D* th1dG_fvtxn0_clus_phi = (TH1D*)fileG->Get("th1d_fvtxn0_clus_phi");
  TH1D* th1dG_fvtxn0_clus_phi_IR = (TH1D*)fileG->Get("th1d_fvtxn0_clus_phi_IR");
  TH1D* th1dG_fvtxn0_clus_phi_OR = (TH1D*)fileG->Get("th1d_fvtxn0_clus_phi_OR");
  th1dG_fvtxn0_clus_phi->Scale(1.0/neventsG);
  th1dG_fvtxn0_clus_phi_IR->Scale(1.0/neventsG_IR);
  th1dG_fvtxn0_clus_phi_OR->Scale(1.0/neventsG_OR);

  TH1D* th1dB_fvtxn0_clus_phi = (TH1D*)fileB->Get("th1d_fvtxn0_clus_phi");
  TH1D* th1dB_fvtxn0_clus_phi_IR = (TH1D*)fileB->Get("th1d_fvtxn0_clus_phi_IR");
  TH1D* th1dB_fvtxn0_clus_phi_OR = (TH1D*)fileB->Get("th1d_fvtxn0_clus_phi_OR");
  th1dB_fvtxn0_clus_phi->Scale(1.0/neventsB);
  th1dB_fvtxn0_clus_phi_IR->Scale(1.0/neventsB_IR);
  th1dB_fvtxn0_clus_phi_OR->Scale(1.0/neventsB_OR);

  // ---
  th1dB_fvtxn0_clus_phi->Draw();
  th1dB_fvtxn0_clus_phi->SetMinimum(0);
  th1dB_fvtxn0_clus_phi->SetTitle(Form("d+Au 20 GeV, all x,y, FVTX North"));
  th1dB_fvtxn0_clus_phi->GetXaxis()->SetTitle("cluster #phi (rad)");
  th1dB_fvtxn0_clus_phi->GetYaxis()->SetTitle("counts/event");
  th1dG_fvtxn0_clus_phi->Draw("same");
  th1dG_fvtxn0_clus_phi->SetLineColor(kRed);
  delete leg;
  leg = new TLegend(0.18,0.68,0.38,0.88);
  leg->AddEntry(th1dB_fvtxn0_clus_phi,"all runs","l");
  leg->AddEntry(th1dG_fvtxn0_clus_phi,"good runs","l");
  leg->SetTextSize(0.045);
  leg->Draw();
  c1->Print(Form("FigsPhi/compbg_fvtxn0_clus_phi_20.png"));
  c1->Print(Form("FigsPhi/compbg_fvtxn0_clus_phi_20.pdf"));

  th1dB_fvtxn0_clus_phi_IR->Draw();
  th1dB_fvtxn0_clus_phi_IR->SetMinimum(0);
  th1dB_fvtxn0_clus_phi_IR->SetTitle(Form("d+Au 20 GeV, x,y within 0.15 cm of beam center, FVTX North"));
  th1dB_fvtxn0_clus_phi_IR->GetXaxis()->SetTitle("cluster #phi (rad)");
  th1dB_fvtxn0_clus_phi_IR->GetYaxis()->SetTitle("counts/event");
  th1dG_fvtxn0_clus_phi_IR->Draw("same");
  th1dG_fvtxn0_clus_phi_IR->SetLineColor(kRed);
  leg->Draw();
  c1->Print(Form("FigsPhi/compbg_fvtxn0_clus_phi_20_IR.png"));
  c1->Print(Form("FigsPhi/compbg_fvtxn0_clus_phi_20_IR.pdf"));

  th1dB_fvtxn0_clus_phi->Draw();
  th1dB_fvtxn0_clus_phi->SetMinimum(0);
  th1dB_fvtxn0_clus_phi->SetTitle(Form("d+Au 20 GeV, FVTX North"));
  th1dB_fvtxn0_clus_phi->GetXaxis()->SetTitle("cluster #phi (rad)");
  th1dB_fvtxn0_clus_phi->GetYaxis()->SetTitle("counts/event");
  th1dG_fvtxn0_clus_phi_IR->Draw("same");
  th1dG_fvtxn0_clus_phi_IR->SetLineColor(kRed);
  delete leg;
  leg = new TLegend(0.18,0.68,0.38,0.88);
  leg->AddEntry(th1dB_fvtxn0_clus_phi,"all runs, no radius cut","l");
  leg->AddEntry(th1dG_fvtxn0_clus_phi_IR,"good runs, with radius cut","l");
  leg->SetTextSize(0.045);
  leg->Draw();
  c1->Print(Form("FigsPhi/compbg_fvtxn0_clus_phi_20_BGIR.png"));
  c1->Print(Form("FigsPhi/compbg_fvtxn0_clus_phi_20_BGIR.pdf"));



  // ---
  // --- fvtxn1
  // ---

  // ---
  TH1D* th1dG_fvtxn1_clus_phi = (TH1D*)fileG->Get("th1d_fvtxn1_clus_phi");
  TH1D* th1dG_fvtxn1_clus_phi_IR = (TH1D*)fileG->Get("th1d_fvtxn1_clus_phi_IR");
  TH1D* th1dG_fvtxn1_clus_phi_OR = (TH1D*)fileG->Get("th1d_fvtxn1_clus_phi_OR");
  th1dG_fvtxn1_clus_phi->Scale(1.0/neventsG);
  th1dG_fvtxn1_clus_phi_IR->Scale(1.0/neventsG_IR);
  th1dG_fvtxn1_clus_phi_OR->Scale(1.0/neventsG_OR);

  TH1D* th1dB_fvtxn1_clus_phi = (TH1D*)fileB->Get("th1d_fvtxn1_clus_phi");
  TH1D* th1dB_fvtxn1_clus_phi_IR = (TH1D*)fileB->Get("th1d_fvtxn1_clus_phi_IR");
  TH1D* th1dB_fvtxn1_clus_phi_OR = (TH1D*)fileB->Get("th1d_fvtxn1_clus_phi_OR");
  th1dB_fvtxn1_clus_phi->Scale(1.0/neventsB);
  th1dB_fvtxn1_clus_phi_IR->Scale(1.0/neventsB_IR);
  th1dB_fvtxn1_clus_phi_OR->Scale(1.0/neventsB_OR);

  // ---
  th1dB_fvtxn1_clus_phi->Draw();
  th1dB_fvtxn1_clus_phi->SetMinimum(0);
  th1dB_fvtxn1_clus_phi->SetTitle(Form("d+Au 20 GeV, all x,y, FVTX North"));
  th1dB_fvtxn1_clus_phi->GetXaxis()->SetTitle("cluster #phi (rad)");
  th1dB_fvtxn1_clus_phi->GetYaxis()->SetTitle("counts/event");
  th1dG_fvtxn1_clus_phi->Draw("same");
  th1dG_fvtxn1_clus_phi->SetLineColor(kRed);
  delete leg;
  leg = new TLegend(0.18,0.68,0.38,0.88);
  leg->AddEntry(th1dB_fvtxn1_clus_phi,"all runs","l");
  leg->AddEntry(th1dG_fvtxn1_clus_phi,"good runs","l");
  leg->SetTextSize(0.045);
  leg->Draw();
  c1->Print(Form("FigsPhi/compbg_fvtxn1_clus_phi_20.png"));
  c1->Print(Form("FigsPhi/compbg_fvtxn1_clus_phi_20.pdf"));

  th1dB_fvtxn1_clus_phi_IR->Draw();
  th1dB_fvtxn1_clus_phi_IR->SetMinimum(0);
  th1dB_fvtxn1_clus_phi_IR->SetTitle(Form("d+Au 20 GeV, x,y within 0.15 cm of beam center, FVTX North"));
  th1dB_fvtxn1_clus_phi_IR->GetXaxis()->SetTitle("cluster #phi (rad)");
  th1dB_fvtxn1_clus_phi_IR->GetYaxis()->SetTitle("counts/event");
  th1dG_fvtxn1_clus_phi_IR->Draw("same");
  th1dG_fvtxn1_clus_phi_IR->SetLineColor(kRed);
  leg->Draw();
  c1->Print(Form("FigsPhi/compbg_fvtxn1_clus_phi_20_IR.png"));
  c1->Print(Form("FigsPhi/compbg_fvtxn1_clus_phi_20_IR.pdf"));

  th1dB_fvtxn1_clus_phi->Draw();
  th1dB_fvtxn1_clus_phi->SetMinimum(0);
  th1dB_fvtxn1_clus_phi->SetTitle(Form("d+Au 20 GeV, FVTX North"));
  th1dB_fvtxn1_clus_phi->GetXaxis()->SetTitle("cluster #phi (rad)");
  th1dB_fvtxn1_clus_phi->GetYaxis()->SetTitle("counts/event");
  th1dG_fvtxn1_clus_phi_IR->Draw("same");
  th1dG_fvtxn1_clus_phi_IR->SetLineColor(kRed);
  delete leg;
  leg = new TLegend(0.18,0.68,0.38,0.88);
  leg->AddEntry(th1dB_fvtxn1_clus_phi,"all runs, no radius cut","l");
  leg->AddEntry(th1dG_fvtxn1_clus_phi_IR,"good runs, with radius cut","l");
  leg->SetTextSize(0.045);
  leg->Draw();
  c1->Print(Form("FigsPhi/compbg_fvtxn1_clus_phi_20_BGIR.png"));
  c1->Print(Form("FigsPhi/compbg_fvtxn1_clus_phi_20_BGIR.pdf"));



  // ---
  // --- fvtxn2
  // ---

  // ---
  TH1D* th1dG_fvtxn2_clus_phi = (TH1D*)fileG->Get("th1d_fvtxn2_clus_phi");
  TH1D* th1dG_fvtxn2_clus_phi_IR = (TH1D*)fileG->Get("th1d_fvtxn2_clus_phi_IR");
  TH1D* th1dG_fvtxn2_clus_phi_OR = (TH1D*)fileG->Get("th1d_fvtxn2_clus_phi_OR");
  th1dG_fvtxn2_clus_phi->Scale(1.0/neventsG);
  th1dG_fvtxn2_clus_phi_IR->Scale(1.0/neventsG_IR);
  th1dG_fvtxn2_clus_phi_OR->Scale(1.0/neventsG_OR);

  TH1D* th1dB_fvtxn2_clus_phi = (TH1D*)fileB->Get("th1d_fvtxn2_clus_phi");
  TH1D* th1dB_fvtxn2_clus_phi_IR = (TH1D*)fileB->Get("th1d_fvtxn2_clus_phi_IR");
  TH1D* th1dB_fvtxn2_clus_phi_OR = (TH1D*)fileB->Get("th1d_fvtxn2_clus_phi_OR");
  th1dB_fvtxn2_clus_phi->Scale(1.0/neventsB);
  th1dB_fvtxn2_clus_phi_IR->Scale(1.0/neventsB_IR);
  th1dB_fvtxn2_clus_phi_OR->Scale(1.0/neventsB_OR);

  // ---
  th1dB_fvtxn2_clus_phi->Draw();
  th1dB_fvtxn2_clus_phi->SetMinimum(0);
  th1dB_fvtxn2_clus_phi->SetTitle(Form("d+Au 20 GeV, all x,y, FVTX North"));
  th1dB_fvtxn2_clus_phi->GetXaxis()->SetTitle("cluster #phi (rad)");
  th1dB_fvtxn2_clus_phi->GetYaxis()->SetTitle("counts/event");
  th1dG_fvtxn2_clus_phi->Draw("same");
  th1dG_fvtxn2_clus_phi->SetLineColor(kRed);
  delete leg;
  leg = new TLegend(0.18,0.68,0.38,0.88);
  leg->AddEntry(th1dB_fvtxn2_clus_phi,"all runs","l");
  leg->AddEntry(th1dG_fvtxn2_clus_phi,"good runs","l");
  leg->SetTextSize(0.045);
  leg->Draw();
  c1->Print(Form("FigsPhi/compbg_fvtxn2_clus_phi_20.png"));
  c1->Print(Form("FigsPhi/compbg_fvtxn2_clus_phi_20.pdf"));

  th1dB_fvtxn2_clus_phi_IR->Draw();
  th1dB_fvtxn2_clus_phi_IR->SetMinimum(0);
  th1dB_fvtxn2_clus_phi_IR->SetTitle(Form("d+Au 20 GeV, x,y within 0.15 cm of beam center, FVTX North"));
  th1dB_fvtxn2_clus_phi_IR->GetXaxis()->SetTitle("cluster #phi (rad)");
  th1dB_fvtxn2_clus_phi_IR->GetYaxis()->SetTitle("counts/event");
  th1dG_fvtxn2_clus_phi_IR->Draw("same");
  th1dG_fvtxn2_clus_phi_IR->SetLineColor(kRed);
  leg->Draw();
  c1->Print(Form("FigsPhi/compbg_fvtxn2_clus_phi_20_IR.png"));
  c1->Print(Form("FigsPhi/compbg_fvtxn2_clus_phi_20_IR.pdf"));

  th1dB_fvtxn2_clus_phi->Draw();
  th1dB_fvtxn2_clus_phi->SetMinimum(0);
  th1dB_fvtxn2_clus_phi->SetTitle(Form("d+Au 20 GeV, FVTX North"));
  th1dB_fvtxn2_clus_phi->GetXaxis()->SetTitle("cluster #phi (rad)");
  th1dB_fvtxn2_clus_phi->GetYaxis()->SetTitle("counts/event");
  th1dG_fvtxn2_clus_phi_IR->Draw("same");
  th1dG_fvtxn2_clus_phi_IR->SetLineColor(kRed);
  delete leg;
  leg = new TLegend(0.18,0.68,0.38,0.88);
  leg->AddEntry(th1dB_fvtxn2_clus_phi,"all runs, no radius cut","l");
  leg->AddEntry(th1dG_fvtxn2_clus_phi_IR,"good runs, with radius cut","l");
  leg->SetTextSize(0.045);
  leg->Draw();
  c1->Print(Form("FigsPhi/compbg_fvtxn2_clus_phi_20_BGIR.png"));
  c1->Print(Form("FigsPhi/compbg_fvtxn2_clus_phi_20_BGIR.pdf"));



  // ---
  // --- fvtxn3
  // ---

  // ---
  TH1D* th1dG_fvtxn3_clus_phi = (TH1D*)fileG->Get("th1d_fvtxn3_clus_phi");
  TH1D* th1dG_fvtxn3_clus_phi_IR = (TH1D*)fileG->Get("th1d_fvtxn3_clus_phi_IR");
  TH1D* th1dG_fvtxn3_clus_phi_OR = (TH1D*)fileG->Get("th1d_fvtxn3_clus_phi_OR");
  th1dG_fvtxn3_clus_phi->Scale(1.0/neventsG);
  th1dG_fvtxn3_clus_phi_IR->Scale(1.0/neventsG_IR);
  th1dG_fvtxn3_clus_phi_OR->Scale(1.0/neventsG_OR);

  TH1D* th1dB_fvtxn3_clus_phi = (TH1D*)fileB->Get("th1d_fvtxn3_clus_phi");
  TH1D* th1dB_fvtxn3_clus_phi_IR = (TH1D*)fileB->Get("th1d_fvtxn3_clus_phi_IR");
  TH1D* th1dB_fvtxn3_clus_phi_OR = (TH1D*)fileB->Get("th1d_fvtxn3_clus_phi_OR");
  th1dB_fvtxn3_clus_phi->Scale(1.0/neventsB);
  th1dB_fvtxn3_clus_phi_IR->Scale(1.0/neventsB_IR);
  th1dB_fvtxn3_clus_phi_OR->Scale(1.0/neventsB_OR);

  // ---
  th1dB_fvtxn3_clus_phi->Draw();
  th1dB_fvtxn3_clus_phi->SetMinimum(0);
  th1dB_fvtxn3_clus_phi->SetTitle(Form("d+Au 20 GeV, all x,y, FVTX North"));
  th1dB_fvtxn3_clus_phi->GetXaxis()->SetTitle("cluster #phi (rad)");
  th1dB_fvtxn3_clus_phi->GetYaxis()->SetTitle("counts/event");
  th1dG_fvtxn3_clus_phi->Draw("same");
  th1dG_fvtxn3_clus_phi->SetLineColor(kRed);
  delete leg;
  leg = new TLegend(0.18,0.68,0.38,0.88);
  leg->AddEntry(th1dB_fvtxn3_clus_phi,"all runs","l");
  leg->AddEntry(th1dG_fvtxn3_clus_phi,"good runs","l");
  leg->SetTextSize(0.045);
  leg->Draw();
  c1->Print(Form("FigsPhi/compbg_fvtxn3_clus_phi_20.png"));
  c1->Print(Form("FigsPhi/compbg_fvtxn3_clus_phi_20.pdf"));

  th1dB_fvtxn3_clus_phi_IR->Draw();
  th1dB_fvtxn3_clus_phi_IR->SetMinimum(0);
  th1dB_fvtxn3_clus_phi_IR->SetTitle(Form("d+Au 20 GeV, x,y within 0.15 cm of beam center, FVTX North"));
  th1dB_fvtxn3_clus_phi_IR->GetXaxis()->SetTitle("cluster #phi (rad)");
  th1dB_fvtxn3_clus_phi_IR->GetYaxis()->SetTitle("counts/event");
  th1dG_fvtxn3_clus_phi_IR->Draw("same");
  th1dG_fvtxn3_clus_phi_IR->SetLineColor(kRed);
  leg->Draw();
  c1->Print(Form("FigsPhi/compbg_fvtxn3_clus_phi_20_IR.png"));
  c1->Print(Form("FigsPhi/compbg_fvtxn3_clus_phi_20_IR.pdf"));

  th1dB_fvtxn3_clus_phi->Draw();
  th1dB_fvtxn3_clus_phi->SetMinimum(0);
  th1dB_fvtxn3_clus_phi->SetTitle(Form("d+Au 20 GeV, FVTX North"));
  th1dB_fvtxn3_clus_phi->GetXaxis()->SetTitle("cluster #phi (rad)");
  th1dB_fvtxn3_clus_phi->GetYaxis()->SetTitle("counts/event");
  th1dG_fvtxn3_clus_phi_IR->Draw("same");
  th1dG_fvtxn3_clus_phi_IR->SetLineColor(kRed);
  delete leg;
  leg = new TLegend(0.18,0.68,0.38,0.88);
  leg->AddEntry(th1dB_fvtxn3_clus_phi,"all runs, no radius cut","l");
  leg->AddEntry(th1dG_fvtxn3_clus_phi_IR,"good runs, with radius cut","l");
  leg->SetTextSize(0.045);
  leg->Draw();
  c1->Print(Form("FigsPhi/compbg_fvtxn3_clus_phi_20_BGIR.png"));
  c1->Print(Form("FigsPhi/compbg_fvtxn3_clus_phi_20_BGIR.pdf"));





  delete c1;

}


void compE()
{

  TCanvas* c1 = new TCanvas("c1","");

  TFile* fileB = TFile::Open(Form("RootFiles/sum20_better.root"));
  TFile* fileG = TFile::Open(Form("RootFiles/sum39.root"));

  TH1D* th1dG_counter = (TH1D*)fileG->Get("th1d_FVTX_nclus");
  TH1D* th1dG_counter_IR = (TH1D*)fileG->Get("th1d_FVTX_nclus_IR");
  TH1D* th1dG_counter_OR = (TH1D*)fileG->Get("th1d_FVTX_nclus_OR");
  int neventsG = th1dG_counter->GetEntries();
  int neventsG_IR = th1dG_counter_IR->GetEntries();
  int neventsG_OR = th1dG_counter_OR->GetEntries();
  TH1D* th1dB_counter = (TH1D*)fileB->Get("th1d_FVTX_nclus");
  TH1D* th1dB_counter_IR = (TH1D*)fileB->Get("th1d_FVTX_nclus_IR");
  TH1D* th1dB_counter_OR = (TH1D*)fileB->Get("th1d_FVTX_nclus_OR");
  int neventsB = th1dB_counter->GetEntries();
  int neventsB_IR = th1dB_counter_IR->GetEntries();
  int neventsB_OR = th1dB_counter_OR->GetEntries();

  // ---
  TH1D* th1dB_fvtxs_clus_phi = (TH1D*)fileB->Get("th1d_fvtxs_clus_phi");
  TH1D* th1dB_fvtxs_clus_phi_IR = (TH1D*)fileB->Get("th1d_fvtxs_clus_phi_IR");
  TH1D* th1dB_fvtxs_clus_phi_OR = (TH1D*)fileB->Get("th1d_fvtxs_clus_phi_OR");
  th1dB_fvtxs_clus_phi->Scale(1.0/neventsB);
  th1dB_fvtxs_clus_phi_IR->Scale(1.0/neventsB_IR);
  th1dB_fvtxs_clus_phi_OR->Scale(1.0/neventsB_OR);

  TH1D* th1dG_fvtxs_clus_phi = (TH1D*)fileG->Get("th1d_fvtxs_clus_phi");
  TH1D* th1dG_fvtxs_clus_phi_IR = (TH1D*)fileG->Get("th1d_fvtxs_clus_phi_IR");
  TH1D* th1dG_fvtxs_clus_phi_OR = (TH1D*)fileG->Get("th1d_fvtxs_clus_phi_OR");
  th1dG_fvtxs_clus_phi->Scale(1.0/neventsG);
  th1dG_fvtxs_clus_phi_IR->Scale(1.0/neventsG_IR);
  th1dG_fvtxs_clus_phi_OR->Scale(1.0/neventsG_OR);

  // ---
  th1dB_fvtxs_clus_phi->Draw();
  th1dB_fvtxs_clus_phi->SetMinimum(0);
  th1dB_fvtxs_clus_phi->SetTitle(Form("d+Au 20 and 39 GeV, all x,y, FVTX South"));
  th1dB_fvtxs_clus_phi->GetXaxis()->SetTitle("cluster #phi (rad)");
  th1dB_fvtxs_clus_phi->GetYaxis()->SetTitle("counts/event");
  th1dG_fvtxs_clus_phi->Draw("same");
  th1dG_fvtxs_clus_phi->SetLineColor(kRed);
  TLegend* leg = new TLegend(0.18,0.68,0.38,0.88);
  leg->AddEntry(th1dB_fvtxs_clus_phi,"20 GeV","l");
  leg->AddEntry(th1dG_fvtxs_clus_phi,"39 GeV","l");
  leg->SetTextSize(0.045);
  leg->Draw();
  c1->Print(Form("FigsPhi/compE_fvtxs_clus_phi_2039.png"));
  c1->Print(Form("FigsPhi/compE_fvtxs_clus_phi_2039.pdf"));

  // ---
  th1dB_fvtxs_clus_phi_IR->Draw();
  th1dB_fvtxs_clus_phi_IR->SetMinimum(0);
  th1dB_fvtxs_clus_phi_IR->SetTitle(Form("d+Au 20 and 39 GeV, x,y within 0.15 cm of beam center, FVTX South"));
  th1dB_fvtxs_clus_phi_IR->GetXaxis()->SetTitle("cluster #phi (rad)");
  th1dB_fvtxs_clus_phi_IR->GetYaxis()->SetTitle("counts/event");
  th1dG_fvtxs_clus_phi_IR->Draw("same");
  th1dG_fvtxs_clus_phi_IR->SetLineColor(kRed);
  leg->Draw();
  c1->Print(Form("FigsPhi/compE_fvtxs_clus_phi_2039_IR.png"));
  c1->Print(Form("FigsPhi/compE_fvtxs_clus_phi_2039_IR.pdf"));

  // ---
  // --- come back here for layer
  // ---

  // ---
  TH1D* th1dB_fvtxs0_clus_phi = (TH1D*)fileB->Get("th1d_fvtxs0_clus_phi");
  TH1D* th1dB_fvtxs0_clus_phi_IR = (TH1D*)fileB->Get("th1d_fvtxs0_clus_phi_IR");
  TH1D* th1dB_fvtxs0_clus_phi_OR = (TH1D*)fileB->Get("th1d_fvtxs0_clus_phi_OR");
  th1dB_fvtxs0_clus_phi->Scale(1.0/neventsB);
  th1dB_fvtxs0_clus_phi_IR->Scale(1.0/neventsB_IR);
  th1dB_fvtxs0_clus_phi_OR->Scale(1.0/neventsB_OR);

  TH1D* th1dG_fvtxs0_clus_phi = (TH1D*)fileG->Get("th1d_fvtxs0_clus_phi");
  TH1D* th1dG_fvtxs0_clus_phi_IR = (TH1D*)fileG->Get("th1d_fvtxs0_clus_phi_IR");
  TH1D* th1dG_fvtxs0_clus_phi_OR = (TH1D*)fileG->Get("th1d_fvtxs0_clus_phi_OR");
  th1dG_fvtxs0_clus_phi->Scale(1.0/neventsG);
  th1dG_fvtxs0_clus_phi_IR->Scale(1.0/neventsG_IR);
  th1dG_fvtxs0_clus_phi_OR->Scale(1.0/neventsG_OR);

  // ---
  th1dB_fvtxs0_clus_phi->Draw();
  th1dB_fvtxs0_clus_phi->SetMinimum(0);
  th1dB_fvtxs0_clus_phi->SetTitle(Form("d+Au 20 and 39 GeV, all x,y, FVTX South"));
  th1dB_fvtxs0_clus_phi->GetXaxis()->SetTitle("cluster #phi (rad)");
  th1dB_fvtxs0_clus_phi->GetYaxis()->SetTitle("counts/event");
  th1dG_fvtxs0_clus_phi->Draw("same");
  th1dG_fvtxs0_clus_phi->SetLineColor(kRed);
  TLegend* leg = new TLegend(0.18,0.68,0.38,0.88);
  leg->AddEntry(th1dB_fvtxs0_clus_phi,"20 GeV","l");
  leg->AddEntry(th1dG_fvtxs0_clus_phi,"39 GeV","l");
  leg->SetTextSize(0.045);
  leg->Draw();
  c1->Print(Form("FigsPhi/compE_fvtxs0_clus_phi_2039.png"));
  c1->Print(Form("FigsPhi/compE_fvtxs0_clus_phi_2039.pdf"));

  // ---
  th1dB_fvtxs0_clus_phi_IR->Draw();
  th1dB_fvtxs0_clus_phi_IR->SetMinimum(0);
  th1dB_fvtxs0_clus_phi_IR->SetTitle(Form("d+Au 20 and 39 GeV, x,y within 0.15 cm of beam center, FVTX South"));
  th1dB_fvtxs0_clus_phi_IR->GetXaxis()->SetTitle("cluster #phi (rad)");
  th1dB_fvtxs0_clus_phi_IR->GetYaxis()->SetTitle("counts/event");
  th1dG_fvtxs0_clus_phi_IR->Draw("same");
  th1dG_fvtxs0_clus_phi_IR->SetLineColor(kRed);
  leg->Draw();
  c1->Print(Form("FigsPhi/compE_fvtxs0_clus_phi_2039_IR.png"));
  c1->Print(Form("FigsPhi/compE_fvtxs0_clus_phi_2039_IR.pdf"));



  // ---
  // --- come back here for layer
  // ---

  // ---
  TH1D* th1dB_fvtxs1_clus_phi = (TH1D*)fileB->Get("th1d_fvtxs1_clus_phi");
  TH1D* th1dB_fvtxs1_clus_phi_IR = (TH1D*)fileB->Get("th1d_fvtxs1_clus_phi_IR");
  TH1D* th1dB_fvtxs1_clus_phi_OR = (TH1D*)fileB->Get("th1d_fvtxs1_clus_phi_OR");
  th1dB_fvtxs1_clus_phi->Scale(1.0/neventsB);
  th1dB_fvtxs1_clus_phi_IR->Scale(1.0/neventsB_IR);
  th1dB_fvtxs1_clus_phi_OR->Scale(1.0/neventsB_OR);

  TH1D* th1dG_fvtxs1_clus_phi = (TH1D*)fileG->Get("th1d_fvtxs1_clus_phi");
  TH1D* th1dG_fvtxs1_clus_phi_IR = (TH1D*)fileG->Get("th1d_fvtxs1_clus_phi_IR");
  TH1D* th1dG_fvtxs1_clus_phi_OR = (TH1D*)fileG->Get("th1d_fvtxs1_clus_phi_OR");
  th1dG_fvtxs1_clus_phi->Scale(1.0/neventsG);
  th1dG_fvtxs1_clus_phi_IR->Scale(1.0/neventsG_IR);
  th1dG_fvtxs1_clus_phi_OR->Scale(1.0/neventsG_OR);

  // ---
  th1dB_fvtxs1_clus_phi->Draw();
  th1dB_fvtxs1_clus_phi->SetMinimum(0);
  th1dB_fvtxs1_clus_phi->SetTitle(Form("d+Au 20 and 39 GeV, all x,y, FVTX South"));
  th1dB_fvtxs1_clus_phi->GetXaxis()->SetTitle("cluster #phi (rad)");
  th1dB_fvtxs1_clus_phi->GetYaxis()->SetTitle("counts/event");
  th1dG_fvtxs1_clus_phi->Draw("same");
  th1dG_fvtxs1_clus_phi->SetLineColor(kRed);
  TLegend* leg = new TLegend(0.18,0.68,0.38,0.88);
  leg->AddEntry(th1dB_fvtxs1_clus_phi,"20 GeV","l");
  leg->AddEntry(th1dG_fvtxs1_clus_phi,"39 GeV","l");
  leg->SetTextSize(0.045);
  leg->Draw();
  c1->Print(Form("FigsPhi/compE_fvtxs1_clus_phi_2039.png"));
  c1->Print(Form("FigsPhi/compE_fvtxs1_clus_phi_2039.pdf"));

  // ---
  th1dB_fvtxs1_clus_phi_IR->Draw();
  th1dB_fvtxs1_clus_phi_IR->SetMinimum(0);
  th1dB_fvtxs1_clus_phi_IR->SetTitle(Form("d+Au 20 and 39 GeV, x,y within 0.15 cm of beam center, FVTX South"));
  th1dB_fvtxs1_clus_phi_IR->GetXaxis()->SetTitle("cluster #phi (rad)");
  th1dB_fvtxs1_clus_phi_IR->GetYaxis()->SetTitle("counts/event");
  th1dG_fvtxs1_clus_phi_IR->Draw("same");
  th1dG_fvtxs1_clus_phi_IR->SetLineColor(kRed);
  leg->Draw();
  c1->Print(Form("FigsPhi/compE_fvtxs1_clus_phi_2039_IR.png"));
  c1->Print(Form("FigsPhi/compE_fvtxs1_clus_phi_2039_IR.pdf"));


  // ---
  // --- come back here for layer
  // ---

  // ---
  TH1D* th1dB_fvtxs2_clus_phi = (TH1D*)fileB->Get("th1d_fvtxs2_clus_phi");
  TH1D* th1dB_fvtxs2_clus_phi_IR = (TH1D*)fileB->Get("th1d_fvtxs2_clus_phi_IR");
  TH1D* th1dB_fvtxs2_clus_phi_OR = (TH1D*)fileB->Get("th1d_fvtxs2_clus_phi_OR");
  th1dB_fvtxs2_clus_phi->Scale(1.0/neventsB);
  th1dB_fvtxs2_clus_phi_IR->Scale(1.0/neventsB_IR);
  th1dB_fvtxs2_clus_phi_OR->Scale(1.0/neventsB_OR);

  TH1D* th1dG_fvtxs2_clus_phi = (TH1D*)fileG->Get("th1d_fvtxs2_clus_phi");
  TH1D* th1dG_fvtxs2_clus_phi_IR = (TH1D*)fileG->Get("th1d_fvtxs2_clus_phi_IR");
  TH1D* th1dG_fvtxs2_clus_phi_OR = (TH1D*)fileG->Get("th1d_fvtxs2_clus_phi_OR");
  th1dG_fvtxs2_clus_phi->Scale(1.0/neventsG);
  th1dG_fvtxs2_clus_phi_IR->Scale(1.0/neventsG_IR);
  th1dG_fvtxs2_clus_phi_OR->Scale(1.0/neventsG_OR);

  // ---
  th1dB_fvtxs2_clus_phi->Draw();
  th1dB_fvtxs2_clus_phi->SetMinimum(0);
  th1dB_fvtxs2_clus_phi->SetTitle(Form("d+Au 20 and 39 GeV, all x,y, FVTX South"));
  th1dB_fvtxs2_clus_phi->GetXaxis()->SetTitle("cluster #phi (rad)");
  th1dB_fvtxs2_clus_phi->GetYaxis()->SetTitle("counts/event");
  th1dG_fvtxs2_clus_phi->Draw("same");
  th1dG_fvtxs2_clus_phi->SetLineColor(kRed);
  TLegend* leg = new TLegend(0.18,0.68,0.38,0.88);
  leg->AddEntry(th1dB_fvtxs2_clus_phi,"20 GeV","l");
  leg->AddEntry(th1dG_fvtxs2_clus_phi,"39 GeV","l");
  leg->SetTextSize(0.045);
  leg->Draw();
  c1->Print(Form("FigsPhi/compE_fvtxs2_clus_phi_2039.png"));
  c1->Print(Form("FigsPhi/compE_fvtxs2_clus_phi_2039.pdf"));

  // ---
  th1dB_fvtxs2_clus_phi_IR->Draw();
  th1dB_fvtxs2_clus_phi_IR->SetMinimum(0);
  th1dB_fvtxs2_clus_phi_IR->SetTitle(Form("d+Au 20 and 39 GeV, x,y within 0.15 cm of beam center, FVTX South"));
  th1dB_fvtxs2_clus_phi_IR->GetXaxis()->SetTitle("cluster #phi (rad)");
  th1dB_fvtxs2_clus_phi_IR->GetYaxis()->SetTitle("counts/event");
  th1dG_fvtxs2_clus_phi_IR->Draw("same");
  th1dG_fvtxs2_clus_phi_IR->SetLineColor(kRed);
  leg->Draw();
  c1->Print(Form("FigsPhi/compE_fvtxs2_clus_phi_2039_IR.png"));
  c1->Print(Form("FigsPhi/compE_fvtxs2_clus_phi_2039_IR.pdf"));

  // ---
  // --- come back here for layer
  // ---

  // ---
  TH1D* th1dB_fvtxs3_clus_phi = (TH1D*)fileB->Get("th1d_fvtxs3_clus_phi");
  TH1D* th1dB_fvtxs3_clus_phi_IR = (TH1D*)fileB->Get("th1d_fvtxs3_clus_phi_IR");
  TH1D* th1dB_fvtxs3_clus_phi_OR = (TH1D*)fileB->Get("th1d_fvtxs3_clus_phi_OR");
  th1dB_fvtxs3_clus_phi->Scale(1.0/neventsB);
  th1dB_fvtxs3_clus_phi_IR->Scale(1.0/neventsB_IR);
  th1dB_fvtxs3_clus_phi_OR->Scale(1.0/neventsB_OR);

  TH1D* th1dG_fvtxs3_clus_phi = (TH1D*)fileG->Get("th1d_fvtxs3_clus_phi");
  TH1D* th1dG_fvtxs3_clus_phi_IR = (TH1D*)fileG->Get("th1d_fvtxs3_clus_phi_IR");
  TH1D* th1dG_fvtxs3_clus_phi_OR = (TH1D*)fileG->Get("th1d_fvtxs3_clus_phi_OR");
  th1dG_fvtxs3_clus_phi->Scale(1.0/neventsG);
  th1dG_fvtxs3_clus_phi_IR->Scale(1.0/neventsG_IR);
  th1dG_fvtxs3_clus_phi_OR->Scale(1.0/neventsG_OR);

  // ---
  th1dB_fvtxs3_clus_phi->Draw();
  th1dB_fvtxs3_clus_phi->SetMinimum(0);
  th1dB_fvtxs3_clus_phi->SetTitle(Form("d+Au 20 and 39 GeV, all x,y, FVTX South"));
  th1dB_fvtxs3_clus_phi->GetXaxis()->SetTitle("cluster #phi (rad)");
  th1dB_fvtxs3_clus_phi->GetYaxis()->SetTitle("counts/event");
  th1dG_fvtxs3_clus_phi->Draw("same");
  th1dG_fvtxs3_clus_phi->SetLineColor(kRed);
  TLegend* leg = new TLegend(0.18,0.68,0.38,0.88);
  leg->AddEntry(th1dB_fvtxs3_clus_phi,"20 GeV","l");
  leg->AddEntry(th1dG_fvtxs3_clus_phi,"39 GeV","l");
  leg->SetTextSize(0.045);
  leg->Draw();
  c1->Print(Form("FigsPhi/compE_fvtxs3_clus_phi_2039.png"));
  c1->Print(Form("FigsPhi/compE_fvtxs3_clus_phi_2039.pdf"));

  // ---
  th1dB_fvtxs3_clus_phi_IR->Draw();
  th1dB_fvtxs3_clus_phi_IR->SetMinimum(0);
  th1dB_fvtxs3_clus_phi_IR->SetTitle(Form("d+Au 20 and 39 GeV, x,y within 0.15 cm of beam center, FVTX South"));
  th1dB_fvtxs3_clus_phi_IR->GetXaxis()->SetTitle("cluster #phi (rad)");
  th1dB_fvtxs3_clus_phi_IR->GetYaxis()->SetTitle("counts/event");
  th1dG_fvtxs3_clus_phi_IR->Draw("same");
  th1dG_fvtxs3_clus_phi_IR->SetLineColor(kRed);
  leg->Draw();
  c1->Print(Form("FigsPhi/compE_fvtxs3_clus_phi_2039_IR.png"));
  c1->Print(Form("FigsPhi/compE_fvtxs3_clus_phi_2039_IR.pdf"));




  // ---
  // --- come back here for north

  // ---
  TH1D* th1dB_fvtxn_clus_phi = (TH1D*)fileB->Get("th1d_fvtxn_clus_phi");
  TH1D* th1dB_fvtxn_clus_phi_IR = (TH1D*)fileB->Get("th1d_fvtxn_clus_phi_IR");
  TH1D* th1dB_fvtxn_clus_phi_OR = (TH1D*)fileB->Get("th1d_fvtxn_clus_phi_OR");
  th1dB_fvtxn_clus_phi->Scale(1.0/neventsB);
  th1dB_fvtxn_clus_phi_IR->Scale(1.0/neventsB_IR);
  th1dB_fvtxn_clus_phi_OR->Scale(1.0/neventsB_OR);

  TH1D* th1dG_fvtxn_clus_phi = (TH1D*)fileG->Get("th1d_fvtxn_clus_phi");
  TH1D* th1dG_fvtxn_clus_phi_IR = (TH1D*)fileG->Get("th1d_fvtxn_clus_phi_IR");
  TH1D* th1dG_fvtxn_clus_phi_OR = (TH1D*)fileG->Get("th1d_fvtxn_clus_phi_OR");
  th1dG_fvtxn_clus_phi->Scale(1.0/neventsG);
  th1dG_fvtxn_clus_phi_IR->Scale(1.0/neventsG_IR);
  th1dG_fvtxn_clus_phi_OR->Scale(1.0/neventsG_OR);

  // ---
  th1dB_fvtxn_clus_phi->Draw();
  th1dB_fvtxn_clus_phi->SetMinimum(0);
  th1dB_fvtxn_clus_phi->SetTitle(Form("d+Au 20 and 39 GeV, all x,y, FVTX North"));
  th1dB_fvtxn_clus_phi->GetXaxis()->SetTitle("cluster #phi (rad)");
  th1dB_fvtxn_clus_phi->GetYaxis()->SetTitle("counts/event");
  th1dG_fvtxn_clus_phi->Draw("same");
  th1dG_fvtxn_clus_phi->SetLineColor(kRed);
  TLegend* leg = new TLegend(0.18,0.68,0.38,0.88);
  leg->AddEntry(th1dB_fvtxn_clus_phi,"20 GeV","l");
  leg->AddEntry(th1dG_fvtxn_clus_phi,"39 GeV","l");
  leg->SetTextSize(0.045);
  leg->Draw();
  c1->Print(Form("FigsPhi/compE_fvtxn_clus_phi_2039.png"));
  c1->Print(Form("FigsPhi/compE_fvtxn_clus_phi_2039.pdf"));

  // ---
  th1dB_fvtxn_clus_phi_IR->Draw();
  th1dB_fvtxn_clus_phi_IR->SetMinimum(0);
  th1dB_fvtxn_clus_phi_IR->SetTitle(Form("d+Au 20 and 39 GeV, x,y within 0.15 cm of beam center, FVTX North"));
  th1dB_fvtxn_clus_phi_IR->GetXaxis()->SetTitle("cluster #phi (rad)");
  th1dB_fvtxn_clus_phi_IR->GetYaxis()->SetTitle("counts/event");
  th1dG_fvtxn_clus_phi_IR->Draw("same");
  th1dG_fvtxn_clus_phi_IR->SetLineColor(kRed);
  leg->Draw();
  c1->Print(Form("FigsPhi/compE_fvtxn_clus_phi_2039_IR.png"));
  c1->Print(Form("FigsPhi/compE_fvtxn_clus_phi_2039_IR.pdf"));

  // ---
  // --- come back here for layer
  // ---

  // ---
  TH1D* th1dB_fvtxn0_clus_phi = (TH1D*)fileB->Get("th1d_fvtxn0_clus_phi");
  TH1D* th1dB_fvtxn0_clus_phi_IR = (TH1D*)fileB->Get("th1d_fvtxn0_clus_phi_IR");
  TH1D* th1dB_fvtxn0_clus_phi_OR = (TH1D*)fileB->Get("th1d_fvtxn0_clus_phi_OR");
  th1dB_fvtxn0_clus_phi->Scale(1.0/neventsB);
  th1dB_fvtxn0_clus_phi_IR->Scale(1.0/neventsB_IR);
  th1dB_fvtxn0_clus_phi_OR->Scale(1.0/neventsB_OR);

  TH1D* th1dG_fvtxn0_clus_phi = (TH1D*)fileG->Get("th1d_fvtxn0_clus_phi");
  TH1D* th1dG_fvtxn0_clus_phi_IR = (TH1D*)fileG->Get("th1d_fvtxn0_clus_phi_IR");
  TH1D* th1dG_fvtxn0_clus_phi_OR = (TH1D*)fileG->Get("th1d_fvtxn0_clus_phi_OR");
  th1dG_fvtxn0_clus_phi->Scale(1.0/neventsG);
  th1dG_fvtxn0_clus_phi_IR->Scale(1.0/neventsG_IR);
  th1dG_fvtxn0_clus_phi_OR->Scale(1.0/neventsG_OR);

  // ---
  th1dB_fvtxn0_clus_phi->Draw();
  th1dB_fvtxn0_clus_phi->SetMinimum(0);
  th1dB_fvtxn0_clus_phi->SetTitle(Form("d+Au 20 and 39 GeV, all x,y, FVTX North"));
  th1dB_fvtxn0_clus_phi->GetXaxis()->SetTitle("cluster #phi (rad)");
  th1dB_fvtxn0_clus_phi->GetYaxis()->SetTitle("counts/event");
  th1dG_fvtxn0_clus_phi->Draw("same");
  th1dG_fvtxn0_clus_phi->SetLineColor(kRed);
  TLegend* leg = new TLegend(0.18,0.68,0.38,0.88);
  leg->AddEntry(th1dB_fvtxn0_clus_phi,"20 GeV","l");
  leg->AddEntry(th1dG_fvtxn0_clus_phi,"39 GeV","l");
  leg->SetTextSize(0.045);
  leg->Draw();
  c1->Print(Form("FigsPhi/compE_fvtxn0_clus_phi_2039.png"));
  c1->Print(Form("FigsPhi/compE_fvtxn0_clus_phi_2039.pdf"));

  // ---
  th1dB_fvtxn0_clus_phi_IR->Draw();
  th1dB_fvtxn0_clus_phi_IR->SetMinimum(0);
  th1dB_fvtxn0_clus_phi_IR->SetTitle(Form("d+Au 20 and 39 GeV, x,y within 0.15 cm of beam center, FVTX North"));
  th1dB_fvtxn0_clus_phi_IR->GetXaxis()->SetTitle("cluster #phi (rad)");
  th1dB_fvtxn0_clus_phi_IR->GetYaxis()->SetTitle("counts/event");
  th1dG_fvtxn0_clus_phi_IR->Draw("same");
  th1dG_fvtxn0_clus_phi_IR->SetLineColor(kRed);
  leg->Draw();
  c1->Print(Form("FigsPhi/compE_fvtxn0_clus_phi_2039_IR.png"));
  c1->Print(Form("FigsPhi/compE_fvtxn0_clus_phi_2039_IR.pdf"));



  // ---
  // --- come back here for layer
  // ---

  // ---
  TH1D* th1dB_fvtxn1_clus_phi = (TH1D*)fileB->Get("th1d_fvtxn1_clus_phi");
  TH1D* th1dB_fvtxn1_clus_phi_IR = (TH1D*)fileB->Get("th1d_fvtxn1_clus_phi_IR");
  TH1D* th1dB_fvtxn1_clus_phi_OR = (TH1D*)fileB->Get("th1d_fvtxn1_clus_phi_OR");
  th1dB_fvtxn1_clus_phi->Scale(1.0/neventsB);
  th1dB_fvtxn1_clus_phi_IR->Scale(1.0/neventsB_IR);
  th1dB_fvtxn1_clus_phi_OR->Scale(1.0/neventsB_OR);

  TH1D* th1dG_fvtxn1_clus_phi = (TH1D*)fileG->Get("th1d_fvtxn1_clus_phi");
  TH1D* th1dG_fvtxn1_clus_phi_IR = (TH1D*)fileG->Get("th1d_fvtxn1_clus_phi_IR");
  TH1D* th1dG_fvtxn1_clus_phi_OR = (TH1D*)fileG->Get("th1d_fvtxn1_clus_phi_OR");
  th1dG_fvtxn1_clus_phi->Scale(1.0/neventsG);
  th1dG_fvtxn1_clus_phi_IR->Scale(1.0/neventsG_IR);
  th1dG_fvtxn1_clus_phi_OR->Scale(1.0/neventsG_OR);

  // ---
  th1dB_fvtxn1_clus_phi->Draw();
  th1dB_fvtxn1_clus_phi->SetMinimum(0);
  th1dB_fvtxn1_clus_phi->SetTitle(Form("d+Au 20 and 39 GeV, all x,y, FVTX North"));
  th1dB_fvtxn1_clus_phi->GetXaxis()->SetTitle("cluster #phi (rad)");
  th1dB_fvtxn1_clus_phi->GetYaxis()->SetTitle("counts/event");
  th1dG_fvtxn1_clus_phi->Draw("same");
  th1dG_fvtxn1_clus_phi->SetLineColor(kRed);
  TLegend* leg = new TLegend(0.18,0.68,0.38,0.88);
  leg->AddEntry(th1dB_fvtxn1_clus_phi,"20 GeV","l");
  leg->AddEntry(th1dG_fvtxn1_clus_phi,"39 GeV","l");
  leg->SetTextSize(0.045);
  leg->Draw();
  c1->Print(Form("FigsPhi/compE_fvtxn1_clus_phi_2039.png"));
  c1->Print(Form("FigsPhi/compE_fvtxn1_clus_phi_2039.pdf"));

  // ---
  th1dB_fvtxn1_clus_phi_IR->Draw();
  th1dB_fvtxn1_clus_phi_IR->SetMinimum(0);
  th1dB_fvtxn1_clus_phi_IR->SetTitle(Form("d+Au 20 and 39 GeV, x,y within 0.15 cm of beam center, FVTX North"));
  th1dB_fvtxn1_clus_phi_IR->GetXaxis()->SetTitle("cluster #phi (rad)");
  th1dB_fvtxn1_clus_phi_IR->GetYaxis()->SetTitle("counts/event");
  th1dG_fvtxn1_clus_phi_IR->Draw("same");
  th1dG_fvtxn1_clus_phi_IR->SetLineColor(kRed);
  leg->Draw();
  c1->Print(Form("FigsPhi/compE_fvtxn1_clus_phi_2039_IR.png"));
  c1->Print(Form("FigsPhi/compE_fvtxn1_clus_phi_2039_IR.pdf"));


  // ---
  // --- come back here for layer
  // ---

  // ---
  TH1D* th1dB_fvtxn2_clus_phi = (TH1D*)fileB->Get("th1d_fvtxn2_clus_phi");
  TH1D* th1dB_fvtxn2_clus_phi_IR = (TH1D*)fileB->Get("th1d_fvtxn2_clus_phi_IR");
  TH1D* th1dB_fvtxn2_clus_phi_OR = (TH1D*)fileB->Get("th1d_fvtxn2_clus_phi_OR");
  th1dB_fvtxn2_clus_phi->Scale(1.0/neventsB);
  th1dB_fvtxn2_clus_phi_IR->Scale(1.0/neventsB_IR);
  th1dB_fvtxn2_clus_phi_OR->Scale(1.0/neventsB_OR);

  TH1D* th1dG_fvtxn2_clus_phi = (TH1D*)fileG->Get("th1d_fvtxn2_clus_phi");
  TH1D* th1dG_fvtxn2_clus_phi_IR = (TH1D*)fileG->Get("th1d_fvtxn2_clus_phi_IR");
  TH1D* th1dG_fvtxn2_clus_phi_OR = (TH1D*)fileG->Get("th1d_fvtxn2_clus_phi_OR");
  th1dG_fvtxn2_clus_phi->Scale(1.0/neventsG);
  th1dG_fvtxn2_clus_phi_IR->Scale(1.0/neventsG_IR);
  th1dG_fvtxn2_clus_phi_OR->Scale(1.0/neventsG_OR);

  // ---
  th1dB_fvtxn2_clus_phi->Draw();
  th1dB_fvtxn2_clus_phi->SetMinimum(0);
  th1dB_fvtxn2_clus_phi->SetTitle(Form("d+Au 20 and 39 GeV, all x,y, FVTX North"));
  th1dB_fvtxn2_clus_phi->GetXaxis()->SetTitle("cluster #phi (rad)");
  th1dB_fvtxn2_clus_phi->GetYaxis()->SetTitle("counts/event");
  th1dG_fvtxn2_clus_phi->Draw("same");
  th1dG_fvtxn2_clus_phi->SetLineColor(kRed);
  TLegend* leg = new TLegend(0.18,0.68,0.38,0.88);
  leg->AddEntry(th1dB_fvtxn2_clus_phi,"20 GeV","l");
  leg->AddEntry(th1dG_fvtxn2_clus_phi,"39 GeV","l");
  leg->SetTextSize(0.045);
  leg->Draw();
  c1->Print(Form("FigsPhi/compE_fvtxn2_clus_phi_2039.png"));
  c1->Print(Form("FigsPhi/compE_fvtxn2_clus_phi_2039.pdf"));

  // ---
  th1dB_fvtxn2_clus_phi_IR->Draw();
  th1dB_fvtxn2_clus_phi_IR->SetMinimum(0);
  th1dB_fvtxn2_clus_phi_IR->SetTitle(Form("d+Au 20 and 39 GeV, x,y within 0.15 cm of beam center, FVTX North"));
  th1dB_fvtxn2_clus_phi_IR->GetXaxis()->SetTitle("cluster #phi (rad)");
  th1dB_fvtxn2_clus_phi_IR->GetYaxis()->SetTitle("counts/event");
  th1dG_fvtxn2_clus_phi_IR->Draw("same");
  th1dG_fvtxn2_clus_phi_IR->SetLineColor(kRed);
  leg->Draw();
  c1->Print(Form("FigsPhi/compE_fvtxn2_clus_phi_2039_IR.png"));
  c1->Print(Form("FigsPhi/compE_fvtxn2_clus_phi_2039_IR.pdf"));

  // ---
  // --- come back here for layer
  // ---

  // ---
  TH1D* th1dB_fvtxn3_clus_phi = (TH1D*)fileB->Get("th1d_fvtxn3_clus_phi");
  TH1D* th1dB_fvtxn3_clus_phi_IR = (TH1D*)fileB->Get("th1d_fvtxn3_clus_phi_IR");
  TH1D* th1dB_fvtxn3_clus_phi_OR = (TH1D*)fileB->Get("th1d_fvtxn3_clus_phi_OR");
  th1dB_fvtxn3_clus_phi->Scale(1.0/neventsB);
  th1dB_fvtxn3_clus_phi_IR->Scale(1.0/neventsB_IR);
  th1dB_fvtxn3_clus_phi_OR->Scale(1.0/neventsB_OR);

  TH1D* th1dG_fvtxn3_clus_phi = (TH1D*)fileG->Get("th1d_fvtxn3_clus_phi");
  TH1D* th1dG_fvtxn3_clus_phi_IR = (TH1D*)fileG->Get("th1d_fvtxn3_clus_phi_IR");
  TH1D* th1dG_fvtxn3_clus_phi_OR = (TH1D*)fileG->Get("th1d_fvtxn3_clus_phi_OR");
  th1dG_fvtxn3_clus_phi->Scale(1.0/neventsG);
  th1dG_fvtxn3_clus_phi_IR->Scale(1.0/neventsG_IR);
  th1dG_fvtxn3_clus_phi_OR->Scale(1.0/neventsG_OR);

  // ---
  th1dB_fvtxn3_clus_phi->Draw();
  th1dB_fvtxn3_clus_phi->SetMinimum(0);
  th1dB_fvtxn3_clus_phi->SetTitle(Form("d+Au 20 and 39 GeV, all x,y, FVTX North"));
  th1dB_fvtxn3_clus_phi->GetXaxis()->SetTitle("cluster #phi (rad)");
  th1dB_fvtxn3_clus_phi->GetYaxis()->SetTitle("counts/event");
  th1dG_fvtxn3_clus_phi->Draw("same");
  th1dG_fvtxn3_clus_phi->SetLineColor(kRed);
  TLegend* leg = new TLegend(0.18,0.68,0.38,0.88);
  leg->AddEntry(th1dB_fvtxn3_clus_phi,"20 GeV","l");
  leg->AddEntry(th1dG_fvtxn3_clus_phi,"39 GeV","l");
  leg->SetTextSize(0.045);
  leg->Draw();
  c1->Print(Form("FigsPhi/compE_fvtxn3_clus_phi_2039.png"));
  c1->Print(Form("FigsPhi/compE_fvtxn3_clus_phi_2039.pdf"));

  // ---
  th1dB_fvtxn3_clus_phi_IR->Draw();
  th1dB_fvtxn3_clus_phi_IR->SetMinimum(0);
  th1dB_fvtxn3_clus_phi_IR->SetTitle(Form("d+Au 20 and 39 GeV, x,y within 0.15 cm of beam center, FVTX North"));
  th1dB_fvtxn3_clus_phi_IR->GetXaxis()->SetTitle("cluster #phi (rad)");
  th1dB_fvtxn3_clus_phi_IR->GetYaxis()->SetTitle("counts/event");
  th1dG_fvtxn3_clus_phi_IR->Draw("same");
  th1dG_fvtxn3_clus_phi_IR->SetLineColor(kRed);
  leg->Draw();
  c1->Print(Form("FigsPhi/compE_fvtxn3_clus_phi_2039_IR.png"));
  c1->Print(Form("FigsPhi/compE_fvtxn3_clus_phi_2039_IR.pdf"));



  delete c1;

}
