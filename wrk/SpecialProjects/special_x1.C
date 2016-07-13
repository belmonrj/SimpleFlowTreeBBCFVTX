void takerun(int);
void takeenergy(int);

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
  //  takeenergy(39);

}

void takerun(int run)
{

  TCanvas* c1 = new TCanvas("c1","");

  TFile* file = TFile::Open(Form("RootFiles/svrb_run%d_pass0.root",run));

  TH1D* th1d_counter = (TH1D*)file->Get("th1d_FVTX_nclus");
  TH1D* th1d_counter_IR = (TH1D*)file->Get("th1d_FVTX_nclus_IR");
  TH1D* th1d_counter_OR = (TH1D*)file->Get("th1d_FVTX_nclus_OR");
  int nevents = th1d_counter->GetEntries();
  int nevents_IR = th1d_counter_IR->GetEntries();
  int nevents_OR = th1d_counter_OR->GetEntries();

  TH1D* th1d_fvtxs_clus_phi = (TH1D*)file->Get("th1d_fvtxs_clus_phi");
  TH1D* th1d_fvtxs_clus_phi_IR = (TH1D*)file->Get("th1d_fvtxs_clus_phi_IR");
  TH1D* th1d_fvtxs_clus_phi_OR = (TH1D*)file->Get("th1d_fvtxs_clus_phi_OR");
  th1d_fvtxs_clus_phi->Scale(1.0/nevents);
  th1d_fvtxs_clus_phi_IR->Scale(1.0/nevents_IR);
  th1d_fvtxs_clus_phi_OR->Scale(1.0/nevents_OR);
  th1d_fvtxs_clus_phi->SetMinimum(0);
  th1d_fvtxs_clus_phi->SetLineColor(kBlack);
  th1d_fvtxs_clus_phi->Draw();
  th1d_fvtxs_clus_phi_OR->SetLineColor(kRed);
  th1d_fvtxs_clus_phi_OR->Draw("same");
  th1d_fvtxs_clus_phi_IR->SetLineColor(kBlue);
  th1d_fvtxs_clus_phi_IR->Draw("same");
  TLegend* leg = new TLegend(0.18,0.68,0.38,0.88);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.045);
  leg->AddEntry(th1d_fvtxs_clus_phi,"all","l");
  leg->AddEntry(th1d_fvtxs_clus_phi_IR,"inside radius","l");
  leg->AddEntry(th1d_fvtxs_clus_phi_OR,"outside radius","l");
  leg->Draw();
  c1->Print(Form("FigsPhi/radiuscut_fvtxs_clus_phi_run%d.png",run));
  c1->Print(Form("FigsPhi/radiuscut_fvtxs_clus_phi_run%d.pdf",run));

  // ---

  TH1D* th1d_fvtxs0_clus_phi = (TH1D*)file->Get("th1d_fvtxs0_clus_phi");
  TH1D* th1d_fvtxs0_clus_phi_IR = (TH1D*)file->Get("th1d_fvtxs0_clus_phi_IR");
  TH1D* th1d_fvtxs0_clus_phi_OR = (TH1D*)file->Get("th1d_fvtxs0_clus_phi_OR");
  th1d_fvtxs0_clus_phi->Scale(1.0/nevents);
  th1d_fvtxs0_clus_phi_IR->Scale(1.0/nevents_IR);
  th1d_fvtxs0_clus_phi_OR->Scale(1.0/nevents_OR);
  th1d_fvtxs0_clus_phi->SetMinimum(0);
  th1d_fvtxs0_clus_phi->SetLineColor(kBlack);
  th1d_fvtxs0_clus_phi->Draw();
  th1d_fvtxs0_clus_phi_OR->SetLineColor(kRed);
  th1d_fvtxs0_clus_phi_OR->Draw("same");
  th1d_fvtxs0_clus_phi_IR->SetLineColor(kBlue);
  th1d_fvtxs0_clus_phi_IR->Draw("same");
  leg->Draw();
  c1->Print(Form("FigsPhi/radiuscut_fvtxs0_clus_phi_run%d.png",run));
  c1->Print(Form("FigsPhi/radiuscut_fvtxs0_clus_phi_run%d.pdf",run));

  // ---

  TH1D* th1d_fvtxs1_clus_phi = (TH1D*)file->Get("th1d_fvtxs1_clus_phi");
  TH1D* th1d_fvtxs1_clus_phi_IR = (TH1D*)file->Get("th1d_fvtxs1_clus_phi_IR");
  TH1D* th1d_fvtxs1_clus_phi_OR = (TH1D*)file->Get("th1d_fvtxs1_clus_phi_OR");
  th1d_fvtxs1_clus_phi->Scale(1.0/nevents);
  th1d_fvtxs1_clus_phi_IR->Scale(1.0/nevents_IR);
  th1d_fvtxs1_clus_phi_OR->Scale(1.0/nevents_OR);
  th1d_fvtxs1_clus_phi->SetMinimum(0);
  th1d_fvtxs1_clus_phi->SetLineColor(kBlack);
  th1d_fvtxs1_clus_phi->Draw();
  th1d_fvtxs1_clus_phi_OR->SetLineColor(kRed);
  th1d_fvtxs1_clus_phi_OR->Draw("same");
  th1d_fvtxs1_clus_phi_IR->SetLineColor(kBlue);
  th1d_fvtxs1_clus_phi_IR->Draw("same");
  leg->Draw();
  c1->Print(Form("FigsPhi/radiuscut_fvtxs1_clus_phi_run%d.png",run));
  c1->Print(Form("FigsPhi/radiuscut_fvtxs1_clus_phi_run%d.pdf",run));

  // ---

  TH1D* th1d_fvtxs2_clus_phi = (TH1D*)file->Get("th1d_fvtxs2_clus_phi");
  TH1D* th1d_fvtxs2_clus_phi_IR = (TH1D*)file->Get("th1d_fvtxs2_clus_phi_IR");
  TH1D* th1d_fvtxs2_clus_phi_OR = (TH1D*)file->Get("th1d_fvtxs2_clus_phi_OR");
  th1d_fvtxs2_clus_phi->Scale(1.0/nevents);
  th1d_fvtxs2_clus_phi_IR->Scale(1.0/nevents_IR);
  th1d_fvtxs2_clus_phi_OR->Scale(1.0/nevents_OR);
  th1d_fvtxs2_clus_phi->SetMinimum(0);
  th1d_fvtxs2_clus_phi->SetLineColor(kBlack);
  th1d_fvtxs2_clus_phi->Draw();
  th1d_fvtxs2_clus_phi_OR->SetLineColor(kRed);
  th1d_fvtxs2_clus_phi_OR->Draw("same");
  th1d_fvtxs2_clus_phi_IR->SetLineColor(kBlue);
  th1d_fvtxs2_clus_phi_IR->Draw("same");
  leg->Draw();
  c1->Print(Form("FigsPhi/radiuscut_fvtxs2_clus_phi_run%d.png",run));
  c1->Print(Form("FigsPhi/radiuscut_fvtxs2_clus_phi_run%d.pdf",run));

  // ---

  TH1D* th1d_fvtxs3_clus_phi = (TH1D*)file->Get("th1d_fvtxs3_clus_phi");
  TH1D* th1d_fvtxs3_clus_phi_IR = (TH1D*)file->Get("th1d_fvtxs3_clus_phi_IR");
  TH1D* th1d_fvtxs3_clus_phi_OR = (TH1D*)file->Get("th1d_fvtxs3_clus_phi_OR");
  th1d_fvtxs3_clus_phi->Scale(1.0/nevents);
  th1d_fvtxs3_clus_phi_IR->Scale(1.0/nevents_IR);
  th1d_fvtxs3_clus_phi_OR->Scale(1.0/nevents_OR);
  th1d_fvtxs3_clus_phi->SetMinimum(0);
  th1d_fvtxs3_clus_phi->SetLineColor(kBlack);
  th1d_fvtxs3_clus_phi->Draw();
  th1d_fvtxs3_clus_phi_OR->SetLineColor(kRed);
  th1d_fvtxs3_clus_phi_OR->Draw("same");
  th1d_fvtxs3_clus_phi_IR->SetLineColor(kBlue);
  th1d_fvtxs3_clus_phi_IR->Draw("same");
  leg->Draw();
  c1->Print(Form("FigsPhi/radiuscut_fvtxs3_clus_phi_run%d.png",run));
  c1->Print(Form("FigsPhi/radiuscut_fvtxs3_clus_phi_run%d.pdf",run));



  // ---
  // --- come back here for north

  TH1D* th1d_fvtxn_clus_phi = (TH1D*)file->Get("th1d_fvtxn_clus_phi");
  TH1D* th1d_fvtxn_clus_phi_IR = (TH1D*)file->Get("th1d_fvtxn_clus_phi_IR");
  TH1D* th1d_fvtxn_clus_phi_OR = (TH1D*)file->Get("th1d_fvtxn_clus_phi_OR");
  th1d_fvtxn_clus_phi->Scale(1.0/nevents);
  th1d_fvtxn_clus_phi_IR->Scale(1.0/nevents_IR);
  th1d_fvtxn_clus_phi_OR->Scale(1.0/nevents_OR);
  th1d_fvtxn_clus_phi->SetMinimum(0);
  th1d_fvtxn_clus_phi->SetLineColor(kBlack);
  th1d_fvtxn_clus_phi->Draw();
  th1d_fvtxn_clus_phi_OR->SetLineColor(kRed);
  th1d_fvtxn_clus_phi_OR->Draw("same");
  th1d_fvtxn_clus_phi_IR->SetLineColor(kBlue);
  th1d_fvtxn_clus_phi_IR->Draw("same");
  TLegend* leg = new TLegend(0.18,0.68,0.38,0.88);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.045);
  leg->AddEntry(th1d_fvtxn_clus_phi,"all","l");
  leg->AddEntry(th1d_fvtxn_clus_phi_IR,"inside radius","l");
  leg->AddEntry(th1d_fvtxn_clus_phi_OR,"outside radius","l");
  leg->Draw();
  c1->Print(Form("FigsPhi/radiuscut_fvtxn_clus_phi_run%d.png",run));
  c1->Print(Form("FigsPhi/radiuscut_fvtxn_clus_phi_run%d.pdf",run));

  // ---

  TH1D* th1d_fvtxn0_clus_phi = (TH1D*)file->Get("th1d_fvtxn0_clus_phi");
  TH1D* th1d_fvtxn0_clus_phi_IR = (TH1D*)file->Get("th1d_fvtxn0_clus_phi_IR");
  TH1D* th1d_fvtxn0_clus_phi_OR = (TH1D*)file->Get("th1d_fvtxn0_clus_phi_OR");
  th1d_fvtxn0_clus_phi->Scale(1.0/nevents);
  th1d_fvtxn0_clus_phi_IR->Scale(1.0/nevents_IR);
  th1d_fvtxn0_clus_phi_OR->Scale(1.0/nevents_OR);
  th1d_fvtxn0_clus_phi->SetMinimum(0);
  th1d_fvtxn0_clus_phi->SetLineColor(kBlack);
  th1d_fvtxn0_clus_phi->Draw();
  th1d_fvtxn0_clus_phi_OR->SetLineColor(kRed);
  th1d_fvtxn0_clus_phi_OR->Draw("same");
  th1d_fvtxn0_clus_phi_IR->SetLineColor(kBlue);
  th1d_fvtxn0_clus_phi_IR->Draw("same");
  leg->Draw();
  c1->Print(Form("FigsPhi/radiuscut_fvtxn0_clus_phi_run%d.png",run));
  c1->Print(Form("FigsPhi/radiuscut_fvtxn0_clus_phi_run%d.pdf",run));

  // ---

  TH1D* th1d_fvtxn1_clus_phi = (TH1D*)file->Get("th1d_fvtxn1_clus_phi");
  TH1D* th1d_fvtxn1_clus_phi_IR = (TH1D*)file->Get("th1d_fvtxn1_clus_phi_IR");
  TH1D* th1d_fvtxn1_clus_phi_OR = (TH1D*)file->Get("th1d_fvtxn1_clus_phi_OR");
  th1d_fvtxn1_clus_phi->Scale(1.0/nevents);
  th1d_fvtxn1_clus_phi_IR->Scale(1.0/nevents_IR);
  th1d_fvtxn1_clus_phi_OR->Scale(1.0/nevents_OR);
  th1d_fvtxn1_clus_phi->SetMinimum(0);
  th1d_fvtxn1_clus_phi->SetLineColor(kBlack);
  th1d_fvtxn1_clus_phi->Draw();
  th1d_fvtxn1_clus_phi_OR->SetLineColor(kRed);
  th1d_fvtxn1_clus_phi_OR->Draw("same");
  th1d_fvtxn1_clus_phi_IR->SetLineColor(kBlue);
  th1d_fvtxn1_clus_phi_IR->Draw("same");
  leg->Draw();
  c1->Print(Form("FigsPhi/radiuscut_fvtxn1_clus_phi_run%d.png",run));
  c1->Print(Form("FigsPhi/radiuscut_fvtxn1_clus_phi_run%d.pdf",run));

  // ---

  TH1D* th1d_fvtxn2_clus_phi = (TH1D*)file->Get("th1d_fvtxn2_clus_phi");
  TH1D* th1d_fvtxn2_clus_phi_IR = (TH1D*)file->Get("th1d_fvtxn2_clus_phi_IR");
  TH1D* th1d_fvtxn2_clus_phi_OR = (TH1D*)file->Get("th1d_fvtxn2_clus_phi_OR");
  th1d_fvtxn2_clus_phi->Scale(1.0/nevents);
  th1d_fvtxn2_clus_phi_IR->Scale(1.0/nevents_IR);
  th1d_fvtxn2_clus_phi_OR->Scale(1.0/nevents_OR);
  th1d_fvtxn2_clus_phi->SetMinimum(0);
  th1d_fvtxn2_clus_phi->SetLineColor(kBlack);
  th1d_fvtxn2_clus_phi->Draw();
  th1d_fvtxn2_clus_phi_OR->SetLineColor(kRed);
  th1d_fvtxn2_clus_phi_OR->Draw("same");
  th1d_fvtxn2_clus_phi_IR->SetLineColor(kBlue);
  th1d_fvtxn2_clus_phi_IR->Draw("same");
  leg->Draw();
  c1->Print(Form("FigsPhi/radiuscut_fvtxn2_clus_phi_run%d.png",run));
  c1->Print(Form("FigsPhi/radiuscut_fvtxn2_clus_phi_run%d.pdf",run));

  // ---

  TH1D* th1d_fvtxn3_clus_phi = (TH1D*)file->Get("th1d_fvtxn3_clus_phi");
  TH1D* th1d_fvtxn3_clus_phi_IR = (TH1D*)file->Get("th1d_fvtxn3_clus_phi_IR");
  TH1D* th1d_fvtxn3_clus_phi_OR = (TH1D*)file->Get("th1d_fvtxn3_clus_phi_OR");
  th1d_fvtxn3_clus_phi->Scale(1.0/nevents);
  th1d_fvtxn3_clus_phi_IR->Scale(1.0/nevents_IR);
  th1d_fvtxn3_clus_phi_OR->Scale(1.0/nevents_OR);
  th1d_fvtxn3_clus_phi->SetMinimum(0);
  th1d_fvtxn3_clus_phi->SetLineColor(kBlack);
  th1d_fvtxn3_clus_phi->Draw();
  th1d_fvtxn3_clus_phi_OR->SetLineColor(kRed);
  th1d_fvtxn3_clus_phi_OR->Draw("same");
  th1d_fvtxn3_clus_phi_IR->SetLineColor(kBlue);
  th1d_fvtxn3_clus_phi_IR->Draw("same");
  leg->Draw();
  c1->Print(Form("FigsPhi/radiuscut_fvtxn3_clus_phi_run%d.png",run));
  c1->Print(Form("FigsPhi/radiuscut_fvtxn3_clus_phi_run%d.pdf",run));





  delete c1;

}



void takeenergy(int energy)
{

  TCanvas* c1 = new TCanvas("c1","");

  TFile* file = TFile::Open(Form("RootFiles/sum%d.root",energy));

  TH1D* th1d_counter = (TH1D*)file->Get("th1d_FVTX_nclus");
  TH1D* th1d_counter_IR = (TH1D*)file->Get("th1d_FVTX_nclus_IR");
  TH1D* th1d_counter_OR = (TH1D*)file->Get("th1d_FVTX_nclus_OR");
  int nevents = th1d_counter->GetEntries();
  int nevents_IR = th1d_counter_IR->GetEntries();
  int nevents_OR = th1d_counter_OR->GetEntries();

  TH1D* th1d_fvtxs_clus_phi = (TH1D*)file->Get("th1d_fvtxs_clus_phi");
  TH1D* th1d_fvtxs_clus_phi_IR = (TH1D*)file->Get("th1d_fvtxs_clus_phi_IR");
  TH1D* th1d_fvtxs_clus_phi_OR = (TH1D*)file->Get("th1d_fvtxs_clus_phi_OR");
  th1d_fvtxs_clus_phi->Scale(1.0/nevents);
  th1d_fvtxs_clus_phi_IR->Scale(1.0/nevents_IR);
  th1d_fvtxs_clus_phi_OR->Scale(1.0/nevents_OR);
  th1d_fvtxs_clus_phi->SetMinimum(0);
  th1d_fvtxs_clus_phi->SetTitle("all layers, FVTX South");
  th1d_fvtxs_clus_phi->GetXaxis()->SetTitle("cluster #phi (rad)");
  th1d_fvtxs_clus_phi->GetYaxis()->SetTitle("counts/event");
  th1d_fvtxs_clus_phi->SetLineColor(kBlack);
  th1d_fvtxs_clus_phi->Draw();
  th1d_fvtxs_clus_phi_OR->SetLineColor(kRed);
  th1d_fvtxs_clus_phi_OR->Draw("same");
  th1d_fvtxs_clus_phi_IR->SetLineColor(kBlue);
  th1d_fvtxs_clus_phi_IR->Draw("same");
  TLegend* leg = new TLegend(0.18,0.68,0.38,0.88);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.045);
  leg->AddEntry(th1d_fvtxs_clus_phi,"all","l");
  leg->AddEntry(th1d_fvtxs_clus_phi_IR,"inside radius","l");
  leg->AddEntry(th1d_fvtxs_clus_phi_OR,"outside radius","l");
  leg->Draw();
  c1->Print(Form("FigsPhi/radiuscut_fvtxs_clus_phi_energy%d.png",energy));
  c1->Print(Form("FigsPhi/radiuscut_fvtxs_clus_phi_energy%d.pdf",energy));

  // ---

  TH1D* th1d_fvtxs0_clus_phi = (TH1D*)file->Get("th1d_fvtxs0_clus_phi");
  TH1D* th1d_fvtxs0_clus_phi_IR = (TH1D*)file->Get("th1d_fvtxs0_clus_phi_IR");
  TH1D* th1d_fvtxs0_clus_phi_OR = (TH1D*)file->Get("th1d_fvtxs0_clus_phi_OR");
  th1d_fvtxs0_clus_phi->Scale(1.0/nevents);
  th1d_fvtxs0_clus_phi_IR->Scale(1.0/nevents_IR);
  th1d_fvtxs0_clus_phi_OR->Scale(1.0/nevents_OR);
  th1d_fvtxs0_clus_phi->SetMinimum(0);
  th1d_fvtxs0_clus_phi->SetTitle("layer 0, FVTX South");
  th1d_fvtxs0_clus_phi->GetXaxis()->SetTitle("cluster #phi (rad)");
  th1d_fvtxs0_clus_phi->GetYaxis()->SetTitle("counts/event");
  th1d_fvtxs0_clus_phi->SetLineColor(kBlack);
  th1d_fvtxs0_clus_phi->Draw();
  th1d_fvtxs0_clus_phi_OR->SetLineColor(kRed);
  th1d_fvtxs0_clus_phi_OR->Draw("same");
  th1d_fvtxs0_clus_phi_IR->SetLineColor(kBlue);
  th1d_fvtxs0_clus_phi_IR->Draw("same");
  leg->Draw();
  c1->Print(Form("FigsPhi/radiuscut_fvtxs0_clus_phi_energy%d.png",energy));
  c1->Print(Form("FigsPhi/radiuscut_fvtxs0_clus_phi_energy%d.pdf",energy));

  // ---

  TH1D* th1d_fvtxs1_clus_phi = (TH1D*)file->Get("th1d_fvtxs1_clus_phi");
  TH1D* th1d_fvtxs1_clus_phi_IR = (TH1D*)file->Get("th1d_fvtxs1_clus_phi_IR");
  TH1D* th1d_fvtxs1_clus_phi_OR = (TH1D*)file->Get("th1d_fvtxs1_clus_phi_OR");
  th1d_fvtxs1_clus_phi->Scale(1.0/nevents);
  th1d_fvtxs1_clus_phi_IR->Scale(1.0/nevents_IR);
  th1d_fvtxs1_clus_phi_OR->Scale(1.0/nevents_OR);
  th1d_fvtxs1_clus_phi->SetMinimum(0);
  th1d_fvtxs1_clus_phi->SetTitle("layer 1, FVTX South");
  th1d_fvtxs1_clus_phi->GetXaxis()->SetTitle("cluster #phi (rad)");
  th1d_fvtxs1_clus_phi->GetYaxis()->SetTitle("counts/event");
  th1d_fvtxs1_clus_phi->SetLineColor(kBlack);
  th1d_fvtxs1_clus_phi->Draw();
  th1d_fvtxs1_clus_phi_OR->SetLineColor(kRed);
  th1d_fvtxs1_clus_phi_OR->Draw("same");
  th1d_fvtxs1_clus_phi_IR->SetLineColor(kBlue);
  th1d_fvtxs1_clus_phi_IR->Draw("same");
  leg->Draw();
  c1->Print(Form("FigsPhi/radiuscut_fvtxs1_clus_phi_energy%d.png",energy));
  c1->Print(Form("FigsPhi/radiuscut_fvtxs1_clus_phi_energy%d.pdf",energy));

  // ---

  TH1D* th1d_fvtxs2_clus_phi = (TH1D*)file->Get("th1d_fvtxs2_clus_phi");
  TH1D* th1d_fvtxs2_clus_phi_IR = (TH1D*)file->Get("th1d_fvtxs2_clus_phi_IR");
  TH1D* th1d_fvtxs2_clus_phi_OR = (TH1D*)file->Get("th1d_fvtxs2_clus_phi_OR");
  th1d_fvtxs2_clus_phi->Scale(1.0/nevents);
  th1d_fvtxs2_clus_phi_IR->Scale(1.0/nevents_IR);
  th1d_fvtxs2_clus_phi_OR->Scale(1.0/nevents_OR);
  th1d_fvtxs2_clus_phi->SetMinimum(0);
  th1d_fvtxs2_clus_phi->SetTitle("layer 2, FVTX South");
  th1d_fvtxs2_clus_phi->GetXaxis()->SetTitle("cluster #phi (rad)");
  th1d_fvtxs2_clus_phi->GetYaxis()->SetTitle("counts/event");
  th1d_fvtxs2_clus_phi->SetLineColor(kBlack);
  th1d_fvtxs2_clus_phi->Draw();
  th1d_fvtxs2_clus_phi_OR->SetLineColor(kRed);
  th1d_fvtxs2_clus_phi_OR->Draw("same");
  th1d_fvtxs2_clus_phi_IR->SetLineColor(kBlue);
  th1d_fvtxs2_clus_phi_IR->Draw("same");
  leg->Draw();
  c1->Print(Form("FigsPhi/radiuscut_fvtxs2_clus_phi_energy%d.png",energy));
  c1->Print(Form("FigsPhi/radiuscut_fvtxs2_clus_phi_energy%d.pdf",energy));

  // ---

  TH1D* th1d_fvtxs3_clus_phi = (TH1D*)file->Get("th1d_fvtxs3_clus_phi");
  TH1D* th1d_fvtxs3_clus_phi_IR = (TH1D*)file->Get("th1d_fvtxs3_clus_phi_IR");
  TH1D* th1d_fvtxs3_clus_phi_OR = (TH1D*)file->Get("th1d_fvtxs3_clus_phi_OR");
  th1d_fvtxs3_clus_phi->Scale(1.0/nevents);
  th1d_fvtxs3_clus_phi_IR->Scale(1.0/nevents_IR);
  th1d_fvtxs3_clus_phi_OR->Scale(1.0/nevents_OR);
  th1d_fvtxs3_clus_phi->SetMinimum(0);
  th1d_fvtxs3_clus_phi->SetTitle("layer 3, FVTX South");
  th1d_fvtxs3_clus_phi->GetXaxis()->SetTitle("cluster #phi (rad)");
  th1d_fvtxs3_clus_phi->GetYaxis()->SetTitle("counts/event");
  th1d_fvtxs3_clus_phi->SetLineColor(kBlack);
  th1d_fvtxs3_clus_phi->Draw();
  th1d_fvtxs3_clus_phi_OR->SetLineColor(kRed);
  th1d_fvtxs3_clus_phi_OR->Draw("same");
  th1d_fvtxs3_clus_phi_IR->SetLineColor(kBlue);
  th1d_fvtxs3_clus_phi_IR->Draw("same");
  leg->Draw();
  c1->Print(Form("FigsPhi/radiuscut_fvtxs3_clus_phi_energy%d.png",energy));
  c1->Print(Form("FigsPhi/radiuscut_fvtxs3_clus_phi_energy%d.pdf",energy));



  // ---
  // --- come back here for north

  TH1D* th1d_fvtxn_clus_phi = (TH1D*)file->Get("th1d_fvtxn_clus_phi");
  TH1D* th1d_fvtxn_clus_phi_IR = (TH1D*)file->Get("th1d_fvtxn_clus_phi_IR");
  TH1D* th1d_fvtxn_clus_phi_OR = (TH1D*)file->Get("th1d_fvtxn_clus_phi_OR");
  th1d_fvtxn_clus_phi->Scale(1.0/nevents);
  th1d_fvtxn_clus_phi_IR->Scale(1.0/nevents_IR);
  th1d_fvtxn_clus_phi_OR->Scale(1.0/nevents_OR);
  th1d_fvtxn_clus_phi->SetMinimum(0);
  th1d_fvtxn_clus_phi->SetTitle("all layers, FVTX North");
  th1d_fvtxn_clus_phi->GetXaxis()->SetTitle("cluster #phi (rad)");
  th1d_fvtxn_clus_phi->GetYaxis()->SetTitle("counts/event");
  th1d_fvtxn_clus_phi->SetLineColor(kBlack);
  th1d_fvtxn_clus_phi->Draw();
  th1d_fvtxn_clus_phi_OR->SetLineColor(kRed);
  th1d_fvtxn_clus_phi_OR->Draw("same");
  th1d_fvtxn_clus_phi_IR->SetLineColor(kBlue);
  th1d_fvtxn_clus_phi_IR->Draw("same");
  TLegend* leg = new TLegend(0.18,0.68,0.38,0.88);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.045);
  leg->AddEntry(th1d_fvtxn_clus_phi,"all","l");
  leg->AddEntry(th1d_fvtxn_clus_phi_IR,"inside radius","l");
  leg->AddEntry(th1d_fvtxn_clus_phi_OR,"outside radius","l");
  leg->Draw();
  c1->Print(Form("FigsPhi/radiuscut_fvtxn_clus_phi_energy%d.png",energy));
  c1->Print(Form("FigsPhi/radiuscut_fvtxn_clus_phi_energy%d.pdf",energy));

  // ---

  TH1D* th1d_fvtxn0_clus_phi = (TH1D*)file->Get("th1d_fvtxn0_clus_phi");
  TH1D* th1d_fvtxn0_clus_phi_IR = (TH1D*)file->Get("th1d_fvtxn0_clus_phi_IR");
  TH1D* th1d_fvtxn0_clus_phi_OR = (TH1D*)file->Get("th1d_fvtxn0_clus_phi_OR");
  th1d_fvtxn0_clus_phi->Scale(1.0/nevents);
  th1d_fvtxn0_clus_phi_IR->Scale(1.0/nevents_IR);
  th1d_fvtxn0_clus_phi_OR->Scale(1.0/nevents_OR);
  th1d_fvtxn0_clus_phi->SetMinimum(0);
  th1d_fvtxn0_clus_phi->SetTitle("layer 0, FVTX North");
  th1d_fvtxn0_clus_phi->GetXaxis()->SetTitle("cluster #phi (rad)");
  th1d_fvtxn0_clus_phi->GetYaxis()->SetTitle("counts/event");
  th1d_fvtxn0_clus_phi->SetLineColor(kBlack);
  th1d_fvtxn0_clus_phi->Draw();
  th1d_fvtxn0_clus_phi_OR->SetLineColor(kRed);
  th1d_fvtxn0_clus_phi_OR->Draw("same");
  th1d_fvtxn0_clus_phi_IR->SetLineColor(kBlue);
  th1d_fvtxn0_clus_phi_IR->Draw("same");
  leg->Draw();
  c1->Print(Form("FigsPhi/radiuscut_fvtxn0_clus_phi_energy%d.png",energy));
  c1->Print(Form("FigsPhi/radiuscut_fvtxn0_clus_phi_energy%d.pdf",energy));

  // ---

  TH1D* th1d_fvtxn1_clus_phi = (TH1D*)file->Get("th1d_fvtxn1_clus_phi");
  TH1D* th1d_fvtxn1_clus_phi_IR = (TH1D*)file->Get("th1d_fvtxn1_clus_phi_IR");
  TH1D* th1d_fvtxn1_clus_phi_OR = (TH1D*)file->Get("th1d_fvtxn1_clus_phi_OR");
  th1d_fvtxn1_clus_phi->Scale(1.0/nevents);
  th1d_fvtxn1_clus_phi_IR->Scale(1.0/nevents_IR);
  th1d_fvtxn1_clus_phi_OR->Scale(1.0/nevents_OR);
  th1d_fvtxn1_clus_phi->SetMinimum(0);
  th1d_fvtxn1_clus_phi->SetTitle("layer 1, FVTX North");
  th1d_fvtxn1_clus_phi->GetXaxis()->SetTitle("cluster #phi (rad)");
  th1d_fvtxn1_clus_phi->GetYaxis()->SetTitle("counts/event");
  th1d_fvtxn1_clus_phi->SetLineColor(kBlack);
  th1d_fvtxn1_clus_phi->Draw();
  th1d_fvtxn1_clus_phi_OR->SetLineColor(kRed);
  th1d_fvtxn1_clus_phi_OR->Draw("same");
  th1d_fvtxn1_clus_phi_IR->SetLineColor(kBlue);
  th1d_fvtxn1_clus_phi_IR->Draw("same");
  leg->Draw();
  c1->Print(Form("FigsPhi/radiuscut_fvtxn1_clus_phi_energy%d.png",energy));
  c1->Print(Form("FigsPhi/radiuscut_fvtxn1_clus_phi_energy%d.pdf",energy));

  // ---

  TH1D* th1d_fvtxn2_clus_phi = (TH1D*)file->Get("th1d_fvtxn2_clus_phi");
  TH1D* th1d_fvtxn2_clus_phi_IR = (TH1D*)file->Get("th1d_fvtxn2_clus_phi_IR");
  TH1D* th1d_fvtxn2_clus_phi_OR = (TH1D*)file->Get("th1d_fvtxn2_clus_phi_OR");
  th1d_fvtxn2_clus_phi->Scale(1.0/nevents);
  th1d_fvtxn2_clus_phi_IR->Scale(1.0/nevents_IR);
  th1d_fvtxn2_clus_phi_OR->Scale(1.0/nevents_OR);
  th1d_fvtxn2_clus_phi->SetMinimum(0);
  th1d_fvtxn2_clus_phi->SetTitle("layer 2, FVTX North");
  th1d_fvtxn2_clus_phi->GetXaxis()->SetTitle("cluster #phi (rad)");
  th1d_fvtxn2_clus_phi->GetYaxis()->SetTitle("counts/event");
  th1d_fvtxn2_clus_phi->SetLineColor(kBlack);
  th1d_fvtxn2_clus_phi->Draw();
  th1d_fvtxn2_clus_phi_OR->SetLineColor(kRed);
  th1d_fvtxn2_clus_phi_OR->Draw("same");
  th1d_fvtxn2_clus_phi_IR->SetLineColor(kBlue);
  th1d_fvtxn2_clus_phi_IR->Draw("same");
  leg->Draw();
  c1->Print(Form("FigsPhi/radiuscut_fvtxn2_clus_phi_energy%d.png",energy));
  c1->Print(Form("FigsPhi/radiuscut_fvtxn2_clus_phi_energy%d.pdf",energy));

  // ---

  TH1D* th1d_fvtxn3_clus_phi = (TH1D*)file->Get("th1d_fvtxn3_clus_phi");
  TH1D* th1d_fvtxn3_clus_phi_IR = (TH1D*)file->Get("th1d_fvtxn3_clus_phi_IR");
  TH1D* th1d_fvtxn3_clus_phi_OR = (TH1D*)file->Get("th1d_fvtxn3_clus_phi_OR");
  th1d_fvtxn3_clus_phi->Scale(1.0/nevents);
  th1d_fvtxn3_clus_phi_IR->Scale(1.0/nevents_IR);
  th1d_fvtxn3_clus_phi_OR->Scale(1.0/nevents_OR);
  th1d_fvtxn3_clus_phi->SetMinimum(0);
  th1d_fvtxn3_clus_phi->SetTitle("layer 3, FVTX North");
  th1d_fvtxn3_clus_phi->GetXaxis()->SetTitle("cluster #phi (rad)");
  th1d_fvtxn3_clus_phi->GetYaxis()->SetTitle("counts/event");
  th1d_fvtxn3_clus_phi->SetLineColor(kBlack);
  th1d_fvtxn3_clus_phi->Draw();
  th1d_fvtxn3_clus_phi_OR->SetLineColor(kRed);
  th1d_fvtxn3_clus_phi_OR->Draw("same");
  th1d_fvtxn3_clus_phi_IR->SetLineColor(kBlue);
  th1d_fvtxn3_clus_phi_IR->Draw("same");
  leg->Draw();
  c1->Print(Form("FigsPhi/radiuscut_fvtxn3_clus_phi_energy%d.png",energy));
  c1->Print(Form("FigsPhi/radiuscut_fvtxn3_clus_phi_energy%d.pdf",energy));




  // ---
  // --- come back here for xy, and maybe eventually etaphi
  // ---


  TH2D* th2d_fvtxs_clus_xy = (TH2D*)file->Get("th2d_fvtxs_clus_xy");
  TH2D* th2d_fvtxs_clus_xy_IR = (TH2D*)file->Get("th2d_fvtxs_clus_xy_IR");
  TH2D* th2d_fvtxs_clus_xy_OR = (TH2D*)file->Get("th2d_fvtxs_clus_xy_OR");
  th2d_fvtxs_clus_xy->Scale(1.0/nevents);
  th2d_fvtxs_clus_xy_IR->Scale(1.0/nevents_IR);
  th2d_fvtxs_clus_xy_OR->Scale(1.0/nevents_OR);
  th2d_fvtxs_clus_xy_IR->Draw("colz");
  c1->SetLogz(0);
  c1->Print(Form("FigsPhi/plot2d_fvtxs_clus_xy_energy%d.png",energy));
  c1->Print(Form("FigsPhi/plot2d_fvtxs_clus_xy_energy%d.pdf",energy));
  c1->SetLogz(1);
  c1->Print(Form("FigsPhi/plot2d_fvtxs_clus_xylog_energy%d.png",energy));
  c1->Print(Form("FigsPhi/plot2d_fvtxs_clus_xylog_energy%d.pdf",energy));

  TH2D* th2d_fvtxs0_clus_xy = (TH2D*)file->Get("th2d_fvtxs0_clus_xy");
  TH2D* th2d_fvtxs0_clus_xy_IR = (TH2D*)file->Get("th2d_fvtxs0_clus_xy_IR");
  TH2D* th2d_fvtxs0_clus_xy_OR = (TH2D*)file->Get("th2d_fvtxs0_clus_xy_OR");
  th2d_fvtxs0_clus_xy->Scale(1.0/nevents);
  th2d_fvtxs0_clus_xy_IR->Scale(1.0/nevents_IR);
  th2d_fvtxs0_clus_xy_OR->Scale(1.0/nevents_OR);
  th2d_fvtxs0_clus_xy_IR->Draw("colz");
  c1->SetLogz(0);
  c1->Print(Form("FigsPhi/plot2d_fvtxs0_clus_xy_energy%d.png",energy));
  c1->Print(Form("FigsPhi/plot2d_fvtxs0_clus_xy_energy%d.pdf",energy));
  c1->SetLogz(1);
  c1->Print(Form("FigsPhi/plot2d_fvtxs0_clus_xylog_energy%d.png",energy));
  c1->Print(Form("FigsPhi/plot2d_fvtxs0_clus_xylog_energy%d.pdf",energy));

  TH2D* th2d_fvtxs1_clus_xy = (TH2D*)file->Get("th2d_fvtxs1_clus_xy");
  TH2D* th2d_fvtxs1_clus_xy_IR = (TH2D*)file->Get("th2d_fvtxs1_clus_xy_IR");
  TH2D* th2d_fvtxs1_clus_xy_OR = (TH2D*)file->Get("th2d_fvtxs1_clus_xy_OR");
  th2d_fvtxs1_clus_xy->Scale(1.0/nevents);
  th2d_fvtxs1_clus_xy_IR->Scale(1.0/nevents_IR);
  th2d_fvtxs1_clus_xy_OR->Scale(1.0/nevents_OR);
  th2d_fvtxs1_clus_xy_IR->Draw("colz");
  c1->SetLogz(0);
  c1->Print(Form("FigsPhi/plot2d_fvtxs1_clus_xy_energy%d.png",energy));
  c1->Print(Form("FigsPhi/plot2d_fvtxs1_clus_xy_energy%d.pdf",energy));
  c1->SetLogz(1);
  c1->Print(Form("FigsPhi/plot2d_fvtxs1_clus_xylog_energy%d.png",energy));
  c1->Print(Form("FigsPhi/plot2d_fvtxs1_clus_xylog_energy%d.pdf",energy));

  TH2D* th2d_fvtxs2_clus_xy = (TH2D*)file->Get("th2d_fvtxs2_clus_xy");
  TH2D* th2d_fvtxs2_clus_xy_IR = (TH2D*)file->Get("th2d_fvtxs2_clus_xy_IR");
  TH2D* th2d_fvtxs2_clus_xy_OR = (TH2D*)file->Get("th2d_fvtxs2_clus_xy_OR");
  th2d_fvtxs2_clus_xy->Scale(1.0/nevents);
  th2d_fvtxs2_clus_xy_IR->Scale(1.0/nevents_IR);
  th2d_fvtxs2_clus_xy_OR->Scale(1.0/nevents_OR);
  th2d_fvtxs2_clus_xy_IR->Draw("colz");
  c1->SetLogz(0);
  c1->Print(Form("FigsPhi/plot2d_fvtxs2_clus_xy_energy%d.png",energy));
  c1->Print(Form("FigsPhi/plot2d_fvtxs2_clus_xy_energy%d.pdf",energy));
  c1->SetLogz(1);
  c1->Print(Form("FigsPhi/plot2d_fvtxs2_clus_xylog_energy%d.png",energy));
  c1->Print(Form("FigsPhi/plot2d_fvtxs2_clus_xylog_energy%d.pdf",energy));

  TH2D* th2d_fvtxs3_clus_xy = (TH2D*)file->Get("th2d_fvtxs3_clus_xy");
  TH2D* th2d_fvtxs3_clus_xy_IR = (TH2D*)file->Get("th2d_fvtxs3_clus_xy_IR");
  TH2D* th2d_fvtxs3_clus_xy_OR = (TH2D*)file->Get("th2d_fvtxs3_clus_xy_OR");
  th2d_fvtxs3_clus_xy->Scale(1.0/nevents);
  th2d_fvtxs3_clus_xy_IR->Scale(1.0/nevents_IR);
  th2d_fvtxs3_clus_xy_OR->Scale(1.0/nevents_OR);
  th2d_fvtxs3_clus_xy_IR->Draw("colz");
  c1->SetLogz(0);
  c1->Print(Form("FigsPhi/plot2d_fvtxs3_clus_xy_energy%d.png",energy));
  c1->Print(Form("FigsPhi/plot2d_fvtxs3_clus_xy_energy%d.pdf",energy));
  c1->SetLogz(1);
  c1->Print(Form("FigsPhi/plot2d_fvtxs3_clus_xylog_energy%d.png",energy));
  c1->Print(Form("FigsPhi/plot2d_fvtxs3_clus_xylog_energy%d.pdf",energy));

  // ---

  TH2D* th2d_fvtxn_clus_xy = (TH2D*)file->Get("th2d_fvtxn_clus_xy");
  TH2D* th2d_fvtxn_clus_xy_IR = (TH2D*)file->Get("th2d_fvtxn_clus_xy_IR");
  TH2D* th2d_fvtxn_clus_xy_OR = (TH2D*)file->Get("th2d_fvtxn_clus_xy_OR");
  th2d_fvtxn_clus_xy->Scale(1.0/nevents);
  th2d_fvtxn_clus_xy_IR->Scale(1.0/nevents_IR);
  th2d_fvtxn_clus_xy_OR->Scale(1.0/nevents_OR);
  th2d_fvtxn_clus_xy_IR->Draw("colz");
  c1->SetLogz(0);
  c1->Print(Form("FigsPhi/plot2d_fvtxn_clus_xy_energy%d.png",energy));
  c1->Print(Form("FigsPhi/plot2d_fvtxn_clus_xy_energy%d.pdf",energy));
  c1->SetLogz(1);
  c1->Print(Form("FigsPhi/plot2d_fvtxn_clus_xylog_energy%d.png",energy));
  c1->Print(Form("FigsPhi/plot2d_fvtxn_clus_xylog_energy%d.pdf",energy));

  TH2D* th2d_fvtxn0_clus_xy = (TH2D*)file->Get("th2d_fvtxn0_clus_xy");
  TH2D* th2d_fvtxn0_clus_xy_IR = (TH2D*)file->Get("th2d_fvtxn0_clus_xy_IR");
  TH2D* th2d_fvtxn0_clus_xy_OR = (TH2D*)file->Get("th2d_fvtxn0_clus_xy_OR");
  th2d_fvtxn0_clus_xy->Scale(1.0/nevents);
  th2d_fvtxn0_clus_xy_IR->Scale(1.0/nevents_IR);
  th2d_fvtxn0_clus_xy_OR->Scale(1.0/nevents_OR);
  th2d_fvtxn0_clus_xy_IR->Draw("colz");
  c1->SetLogz(0);
  c1->Print(Form("FigsPhi/plot2d_fvtxn0_clus_xy_energy%d.png",energy));
  c1->Print(Form("FigsPhi/plot2d_fvtxn0_clus_xy_energy%d.pdf",energy));
  c1->SetLogz(1);
  c1->Print(Form("FigsPhi/plot2d_fvtxn0_clus_xylog_energy%d.png",energy));
  c1->Print(Form("FigsPhi/plot2d_fvtxn0_clus_xylog_energy%d.pdf",energy));

  TH2D* th2d_fvtxn1_clus_xy = (TH2D*)file->Get("th2d_fvtxn1_clus_xy");
  TH2D* th2d_fvtxn1_clus_xy_IR = (TH2D*)file->Get("th2d_fvtxn1_clus_xy_IR");
  TH2D* th2d_fvtxn1_clus_xy_OR = (TH2D*)file->Get("th2d_fvtxn1_clus_xy_OR");
  th2d_fvtxn1_clus_xy->Scale(1.0/nevents);
  th2d_fvtxn1_clus_xy_IR->Scale(1.0/nevents_IR);
  th2d_fvtxn1_clus_xy_OR->Scale(1.0/nevents_OR);
  th2d_fvtxn1_clus_xy_IR->Draw("colz");
  c1->SetLogz(0);
  c1->Print(Form("FigsPhi/plot2d_fvtxn1_clus_xy_energy%d.png",energy));
  c1->Print(Form("FigsPhi/plot2d_fvtxn1_clus_xy_energy%d.pdf",energy));
  c1->SetLogz(1);
  c1->Print(Form("FigsPhi/plot2d_fvtxn1_clus_xylog_energy%d.png",energy));
  c1->Print(Form("FigsPhi/plot2d_fvtxn1_clus_xylog_energy%d.pdf",energy));

  TH2D* th2d_fvtxn2_clus_xy = (TH2D*)file->Get("th2d_fvtxn2_clus_xy");
  TH2D* th2d_fvtxn2_clus_xy_IR = (TH2D*)file->Get("th2d_fvtxn2_clus_xy_IR");
  TH2D* th2d_fvtxn2_clus_xy_OR = (TH2D*)file->Get("th2d_fvtxn2_clus_xy_OR");
  th2d_fvtxn2_clus_xy->Scale(1.0/nevents);
  th2d_fvtxn2_clus_xy_IR->Scale(1.0/nevents_IR);
  th2d_fvtxn2_clus_xy_OR->Scale(1.0/nevents_OR);
  th2d_fvtxn2_clus_xy_IR->Draw("colz");
  c1->SetLogz(0);
  c1->Print(Form("FigsPhi/plot2d_fvtxn2_clus_xy_energy%d.png",energy));
  c1->Print(Form("FigsPhi/plot2d_fvtxn2_clus_xy_energy%d.pdf",energy));
  c1->SetLogz(1);
  c1->Print(Form("FigsPhi/plot2d_fvtxn2_clus_xylog_energy%d.png",energy));
  c1->Print(Form("FigsPhi/plot2d_fvtxn2_clus_xylog_energy%d.pdf",energy));

  TH2D* th2d_fvtxn3_clus_xy = (TH2D*)file->Get("th2d_fvtxn3_clus_xy");
  TH2D* th2d_fvtxn3_clus_xy_IR = (TH2D*)file->Get("th2d_fvtxn3_clus_xy_IR");
  TH2D* th2d_fvtxn3_clus_xy_OR = (TH2D*)file->Get("th2d_fvtxn3_clus_xy_OR");
  th2d_fvtxn3_clus_xy->Scale(1.0/nevents);
  th2d_fvtxn3_clus_xy_IR->Scale(1.0/nevents_IR);
  th2d_fvtxn3_clus_xy_OR->Scale(1.0/nevents_OR);
  th2d_fvtxn3_clus_xy_IR->Draw("colz");
  c1->SetLogz(0);
  c1->Print(Form("FigsPhi/plot2d_fvtxn3_clus_xy_energy%d.png",energy));
  c1->Print(Form("FigsPhi/plot2d_fvtxn3_clus_xy_energy%d.pdf",energy));
  c1->SetLogz(1);
  c1->Print(Form("FigsPhi/plot2d_fvtxn3_clus_xylog_energy%d.png",energy));
  c1->Print(Form("FigsPhi/plot2d_fvtxn3_clus_xylog_energy%d.pdf",energy));



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
  c1->Print(Form("FigsPhi/plot2d_fvtxs_clus_phieta_energy%d.png",energy));
  c1->Print(Form("FigsPhi/plot2d_fvtxs_clus_phieta_energy%d.pdf",energy));
  c1->SetLogz(1);
  c1->Print(Form("FigsPhi/plot2d_fvtxs_clus_phietalog_energy%d.png",energy));
  c1->Print(Form("FigsPhi/plot2d_fvtxs_clus_phietalog_energy%d.pdf",energy));

  TH2D* th2d_fvtxs0_clus_phieta = (TH2D*)file->Get("th2d_fvtxs0_clus_phieta");
  TH2D* th2d_fvtxs0_clus_phieta_IR = (TH2D*)file->Get("th2d_fvtxs0_clus_phieta_IR");
  TH2D* th2d_fvtxs0_clus_phieta_OR = (TH2D*)file->Get("th2d_fvtxs0_clus_phieta_OR");
  th2d_fvtxs0_clus_phieta->Scale(1.0/nevents);
  th2d_fvtxs0_clus_phieta_IR->Scale(1.0/nevents_IR);
  th2d_fvtxs0_clus_phieta_OR->Scale(1.0/nevents_OR);
  th2d_fvtxs0_clus_phieta_IR->Draw("colz");
  c1->SetLogz(0);
  c1->Print(Form("FigsPhi/plot2d_fvtxs0_clus_phieta_energy%d.png",energy));
  c1->Print(Form("FigsPhi/plot2d_fvtxs0_clus_phieta_energy%d.pdf",energy));
  c1->SetLogz(1);
  c1->Print(Form("FigsPhi/plot2d_fvtxs0_clus_phietalog_energy%d.png",energy));
  c1->Print(Form("FigsPhi/plot2d_fvtxs0_clus_phietalog_energy%d.pdf",energy));

  TH2D* th2d_fvtxs1_clus_phieta = (TH2D*)file->Get("th2d_fvtxs1_clus_phieta");
  TH2D* th2d_fvtxs1_clus_phieta_IR = (TH2D*)file->Get("th2d_fvtxs1_clus_phieta_IR");
  TH2D* th2d_fvtxs1_clus_phieta_OR = (TH2D*)file->Get("th2d_fvtxs1_clus_phieta_OR");
  th2d_fvtxs1_clus_phieta->Scale(1.0/nevents);
  th2d_fvtxs1_clus_phieta_IR->Scale(1.0/nevents_IR);
  th2d_fvtxs1_clus_phieta_OR->Scale(1.0/nevents_OR);
  th2d_fvtxs1_clus_phieta_IR->Draw("colz");
  c1->SetLogz(0);
  c1->Print(Form("FigsPhi/plot2d_fvtxs1_clus_phieta_energy%d.png",energy));
  c1->Print(Form("FigsPhi/plot2d_fvtxs1_clus_phieta_energy%d.pdf",energy));
  c1->SetLogz(1);
  c1->Print(Form("FigsPhi/plot2d_fvtxs1_clus_phietalog_energy%d.png",energy));
  c1->Print(Form("FigsPhi/plot2d_fvtxs1_clus_phietalog_energy%d.pdf",energy));

  TH2D* th2d_fvtxs2_clus_phieta = (TH2D*)file->Get("th2d_fvtxs2_clus_phieta");
  TH2D* th2d_fvtxs2_clus_phieta_IR = (TH2D*)file->Get("th2d_fvtxs2_clus_phieta_IR");
  TH2D* th2d_fvtxs2_clus_phieta_OR = (TH2D*)file->Get("th2d_fvtxs2_clus_phieta_OR");
  th2d_fvtxs2_clus_phieta->Scale(1.0/nevents);
  th2d_fvtxs2_clus_phieta_IR->Scale(1.0/nevents_IR);
  th2d_fvtxs2_clus_phieta_OR->Scale(1.0/nevents_OR);
  th2d_fvtxs2_clus_phieta_IR->Draw("colz");
  c1->SetLogz(0);
  c1->Print(Form("FigsPhi/plot2d_fvtxs2_clus_phieta_energy%d.png",energy));
  c1->Print(Form("FigsPhi/plot2d_fvtxs2_clus_phieta_energy%d.pdf",energy));
  c1->SetLogz(1);
  c1->Print(Form("FigsPhi/plot2d_fvtxs2_clus_phietalog_energy%d.png",energy));
  c1->Print(Form("FigsPhi/plot2d_fvtxs2_clus_phietalog_energy%d.pdf",energy));

  TH2D* th2d_fvtxs3_clus_phieta = (TH2D*)file->Get("th2d_fvtxs3_clus_phieta");
  TH2D* th2d_fvtxs3_clus_phieta_IR = (TH2D*)file->Get("th2d_fvtxs3_clus_phieta_IR");
  TH2D* th2d_fvtxs3_clus_phieta_OR = (TH2D*)file->Get("th2d_fvtxs3_clus_phieta_OR");
  th2d_fvtxs3_clus_phieta->Scale(1.0/nevents);
  th2d_fvtxs3_clus_phieta_IR->Scale(1.0/nevents_IR);
  th2d_fvtxs3_clus_phieta_OR->Scale(1.0/nevents_OR);
  th2d_fvtxs3_clus_phieta_IR->Draw("colz");
  c1->SetLogz(0);
  c1->Print(Form("FigsPhi/plot2d_fvtxs3_clus_phieta_energy%d.png",energy));
  c1->Print(Form("FigsPhi/plot2d_fvtxs3_clus_phieta_energy%d.pdf",energy));
  c1->SetLogz(1);
  c1->Print(Form("FigsPhi/plot2d_fvtxs3_clus_phietalog_energy%d.png",energy));
  c1->Print(Form("FigsPhi/plot2d_fvtxs3_clus_phietalog_energy%d.pdf",energy));

  // ---

  TH2D* th2d_fvtxn_clus_phieta = (TH2D*)file->Get("th2d_fvtxn_clus_phieta");
  TH2D* th2d_fvtxn_clus_phieta_IR = (TH2D*)file->Get("th2d_fvtxn_clus_phieta_IR");
  TH2D* th2d_fvtxn_clus_phieta_OR = (TH2D*)file->Get("th2d_fvtxn_clus_phieta_OR");
  th2d_fvtxn_clus_phieta->Scale(1.0/nevents);
  th2d_fvtxn_clus_phieta_IR->Scale(1.0/nevents_IR);
  th2d_fvtxn_clus_phieta_OR->Scale(1.0/nevents_OR);
  th2d_fvtxn_clus_phieta_IR->Draw("colz");
  c1->SetLogz(0);
  c1->Print(Form("FigsPhi/plot2d_fvtxn_clus_phieta_energy%d.png",energy));
  c1->Print(Form("FigsPhi/plot2d_fvtxn_clus_phieta_energy%d.pdf",energy));
  c1->SetLogz(1);
  c1->Print(Form("FigsPhi/plot2d_fvtxn_clus_phietalog_energy%d.png",energy));
  c1->Print(Form("FigsPhi/plot2d_fvtxn_clus_phietalog_energy%d.pdf",energy));

  TH2D* th2d_fvtxn0_clus_phieta = (TH2D*)file->Get("th2d_fvtxn0_clus_phieta");
  TH2D* th2d_fvtxn0_clus_phieta_IR = (TH2D*)file->Get("th2d_fvtxn0_clus_phieta_IR");
  TH2D* th2d_fvtxn0_clus_phieta_OR = (TH2D*)file->Get("th2d_fvtxn0_clus_phieta_OR");
  th2d_fvtxn0_clus_phieta->Scale(1.0/nevents);
  th2d_fvtxn0_clus_phieta_IR->Scale(1.0/nevents_IR);
  th2d_fvtxn0_clus_phieta_OR->Scale(1.0/nevents_OR);
  th2d_fvtxn0_clus_phieta_IR->Draw("colz");
  c1->SetLogz(0);
  c1->Print(Form("FigsPhi/plot2d_fvtxn0_clus_phieta_energy%d.png",energy));
  c1->Print(Form("FigsPhi/plot2d_fvtxn0_clus_phieta_energy%d.pdf",energy));
  c1->SetLogz(1);
  c1->Print(Form("FigsPhi/plot2d_fvtxn0_clus_phietalog_energy%d.png",energy));
  c1->Print(Form("FigsPhi/plot2d_fvtxn0_clus_phietalog_energy%d.pdf",energy));

  TH2D* th2d_fvtxn1_clus_phieta = (TH2D*)file->Get("th2d_fvtxn1_clus_phieta");
  TH2D* th2d_fvtxn1_clus_phieta_IR = (TH2D*)file->Get("th2d_fvtxn1_clus_phieta_IR");
  TH2D* th2d_fvtxn1_clus_phieta_OR = (TH2D*)file->Get("th2d_fvtxn1_clus_phieta_OR");
  th2d_fvtxn1_clus_phieta->Scale(1.0/nevents);
  th2d_fvtxn1_clus_phieta_IR->Scale(1.0/nevents_IR);
  th2d_fvtxn1_clus_phieta_OR->Scale(1.0/nevents_OR);
  th2d_fvtxn1_clus_phieta_IR->Draw("colz");
  c1->SetLogz(0);
  c1->Print(Form("FigsPhi/plot2d_fvtxn1_clus_phieta_energy%d.png",energy));
  c1->Print(Form("FigsPhi/plot2d_fvtxn1_clus_phieta_energy%d.pdf",energy));
  c1->SetLogz(1);
  c1->Print(Form("FigsPhi/plot2d_fvtxn1_clus_phietalog_energy%d.png",energy));
  c1->Print(Form("FigsPhi/plot2d_fvtxn1_clus_phietalog_energy%d.pdf",energy));

  TH2D* th2d_fvtxn2_clus_phieta = (TH2D*)file->Get("th2d_fvtxn2_clus_phieta");
  TH2D* th2d_fvtxn2_clus_phieta_IR = (TH2D*)file->Get("th2d_fvtxn2_clus_phieta_IR");
  TH2D* th2d_fvtxn2_clus_phieta_OR = (TH2D*)file->Get("th2d_fvtxn2_clus_phieta_OR");
  th2d_fvtxn2_clus_phieta->Scale(1.0/nevents);
  th2d_fvtxn2_clus_phieta_IR->Scale(1.0/nevents_IR);
  th2d_fvtxn2_clus_phieta_OR->Scale(1.0/nevents_OR);
  th2d_fvtxn2_clus_phieta_IR->Draw("colz");
  c1->SetLogz(0);
  c1->Print(Form("FigsPhi/plot2d_fvtxn2_clus_phieta_energy%d.png",energy));
  c1->Print(Form("FigsPhi/plot2d_fvtxn2_clus_phieta_energy%d.pdf",energy));
  c1->SetLogz(1);
  c1->Print(Form("FigsPhi/plot2d_fvtxn2_clus_phietalog_energy%d.png",energy));
  c1->Print(Form("FigsPhi/plot2d_fvtxn2_clus_phietalog_energy%d.pdf",energy));

  TH2D* th2d_fvtxn3_clus_phieta = (TH2D*)file->Get("th2d_fvtxn3_clus_phieta");
  TH2D* th2d_fvtxn3_clus_phieta_IR = (TH2D*)file->Get("th2d_fvtxn3_clus_phieta_IR");
  TH2D* th2d_fvtxn3_clus_phieta_OR = (TH2D*)file->Get("th2d_fvtxn3_clus_phieta_OR");
  th2d_fvtxn3_clus_phieta->Scale(1.0/nevents);
  th2d_fvtxn3_clus_phieta_IR->Scale(1.0/nevents_IR);
  th2d_fvtxn3_clus_phieta_OR->Scale(1.0/nevents_OR);
  th2d_fvtxn3_clus_phieta_IR->Draw("colz");
  c1->SetLogz(0);
  c1->Print(Form("FigsPhi/plot2d_fvtxn3_clus_phieta_energy%d.png",energy));
  c1->Print(Form("FigsPhi/plot2d_fvtxn3_clus_phieta_energy%d.pdf",energy));
  c1->SetLogz(1);
  c1->Print(Form("FigsPhi/plot2d_fvtxn3_clus_phietalog_energy%d.png",energy));
  c1->Print(Form("FigsPhi/plot2d_fvtxn3_clus_phietalog_energy%d.pdf",energy));



  delete c1;

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
