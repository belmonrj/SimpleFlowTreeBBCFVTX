void doe(int);


void phidist_fvtx_layers()
{

  doe(62);
  doe(39);
  doe(20);

}


void doe(int energy)
{

  TCanvas* c1 = new TCanvas("c1","");

  TFile* file = TFile::Open(Form("input/combined_%d.root",energy));

  TH1D* th1d_south_layerA = (TH1D*)file->Get("th1d_fvtxs_clus_phi");
  TH1D* th1d_south_layer0 = (TH1D*)file->Get("th1d_fvtxs0_clus_phi");
  TH1D* th1d_south_layer1 = (TH1D*)file->Get("th1d_fvtxs1_clus_phi");
  TH1D* th1d_south_layer2 = (TH1D*)file->Get("th1d_fvtxs2_clus_phi");
  TH1D* th1d_south_layer3 = (TH1D*)file->Get("th1d_fvtxs3_clus_phi");

  th1d_south_layerA->Draw();
  th1d_south_layerA->SetMinimum(0);
  th1d_south_layerA->GetXaxis()->SetTitle("cluster #phi (rad)");
  th1d_south_layerA->GetYaxis()->SetTitle("counts");
  c1->Print(Form("cluster_phidist_fvtxs_energy%d.png",energy));
  c1->Print(Form("cluster_phidist_fvtxs_energy%d.pdf",energy));

  th1d_south_layer0->Draw();
  th1d_south_layer0->SetMinimum(0);
  th1d_south_layer0->GetXaxis()->SetTitle("cluster #phi (rad)");
  th1d_south_layer0->GetYaxis()->SetTitle("counts");
  c1->Print(Form("cluster_phidist_fvtxs0_energy%d.png",energy));
  c1->Print(Form("cluster_phidist_fvtxs0_energy%d.pdf",energy));

  th1d_south_layer1->Draw();
  th1d_south_layer1->SetMinimum(0);
  th1d_south_layer1->GetXaxis()->SetTitle("cluster #phi (rad)");
  th1d_south_layer1->GetYaxis()->SetTitle("counts");
  c1->Print(Form("cluster_phidist_fvtxs1_energy%d.png",energy));
  c1->Print(Form("cluster_phidist_fvtxs1_energy%d.pdf",energy));

  th1d_south_layer2->Draw();
  th1d_south_layer2->SetMinimum(0);
  th1d_south_layer2->GetXaxis()->SetTitle("cluster #phi (rad)");
  th1d_south_layer2->GetYaxis()->SetTitle("counts");
  c1->Print(Form("cluster_phidist_fvtxs2_energy%d.png",energy));
  c1->Print(Form("cluster_phidist_fvtxs2_energy%d.pdf",energy));

  th1d_south_layer3->Draw();
  th1d_south_layer3->SetMinimum(0);
  th1d_south_layer3->GetXaxis()->SetTitle("cluster #phi (rad)");
  th1d_south_layer3->GetYaxis()->SetTitle("counts");
  c1->Print(Form("cluster_phidist_fvtxs3_energy%d.png",energy));
  c1->Print(Form("cluster_phidist_fvtxs3_energy%d.pdf",energy));

  // ---

  TH1D* th1d_north_layerA = (TH1D*)file->Get("th1d_fvtxn_clus_phi");
  TH1D* th1d_north_layer0 = (TH1D*)file->Get("th1d_fvtxn0_clus_phi");
  TH1D* th1d_north_layer1 = (TH1D*)file->Get("th1d_fvtxn1_clus_phi");
  TH1D* th1d_north_layer2 = (TH1D*)file->Get("th1d_fvtxn2_clus_phi");
  TH1D* th1d_north_layer3 = (TH1D*)file->Get("th1d_fvtxn3_clus_phi");

  th1d_north_layerA->Draw();
  th1d_north_layerA->SetMinimum(0);
  th1d_north_layerA->GetXaxis()->SetTitle("cluster #phi (rad)");
  th1d_north_layerA->GetYaxis()->SetTitle("counts");
  c1->Print(Form("cluster_phidist_fvtxn_energy%d.png",energy));
  c1->Print(Form("cluster_phidist_fvtxn_energy%d.pdf",energy));

  th1d_north_layer0->Draw();
  th1d_north_layer0->SetMinimum(0);
  th1d_north_layer0->GetXaxis()->SetTitle("cluster #phi (rad)");
  th1d_north_layer0->GetYaxis()->SetTitle("counts");
  c1->Print(Form("cluster_phidist_fvtxn0_energy%d.png",energy));
  c1->Print(Form("cluster_phidist_fvtxn0_energy%d.pdf",energy));

  th1d_north_layer1->Draw();
  th1d_north_layer1->SetMinimum(0);
  th1d_north_layer1->GetXaxis()->SetTitle("cluster #phi (rad)");
  th1d_north_layer1->GetYaxis()->SetTitle("counts");
  c1->Print(Form("cluster_phidist_fvtxn1_energy%d.png",energy));
  c1->Print(Form("cluster_phidist_fvtxn1_energy%d.pdf",energy));

  th1d_north_layer2->Draw();
  th1d_north_layer2->SetMinimum(0);
  th1d_north_layer2->GetXaxis()->SetTitle("cluster #phi (rad)");
  th1d_north_layer2->GetYaxis()->SetTitle("counts");
  c1->Print(Form("cluster_phidist_fvtxn2_energy%d.png",energy));
  c1->Print(Form("cluster_phidist_fvtxn2_energy%d.pdf",energy));

  th1d_north_layer3->Draw();
  th1d_north_layer3->SetMinimum(0);
  th1d_north_layer3->GetXaxis()->SetTitle("cluster #phi (rad)");
  th1d_north_layer3->GetYaxis()->SetTitle("counts");
  c1->Print(Form("cluster_phidist_fvtxn3_energy%d.png",energy));
  c1->Print(Form("cluster_phidist_fvtxn3_energy%d.pdf",energy));



  delete c1;

}
