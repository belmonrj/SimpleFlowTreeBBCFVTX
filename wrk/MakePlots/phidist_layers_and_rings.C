void doe(int);


void phidist_layers_and_rings()
{

  doe(62);
  doe(39);
  doe(20);

}


void doe(int energy)
{

  gStyle->SetOptTitle(1);

  TCanvas* c1 = new TCanvas("c1","");

  TFile* file = TFile::Open(Form("input/combined_%d.root",energy));

  TH1D* th1d_bbc_charge = (TH1D*)file->Get("th1d_BBC_charge");
  int nevents = th1d_bbc_charge->GetEntries();

  TH1D* th1d_south_layerA = (TH1D*)file->Get("th1d_fvtxs_clus_phi");
  TH1D* th1d_south_layer0 = (TH1D*)file->Get("th1d_fvtxs0_clus_phi");
  TH1D* th1d_south_layer1 = (TH1D*)file->Get("th1d_fvtxs1_clus_phi");
  TH1D* th1d_south_layer2 = (TH1D*)file->Get("th1d_fvtxs2_clus_phi");
  TH1D* th1d_south_layer3 = (TH1D*)file->Get("th1d_fvtxs3_clus_phi");

  th1d_south_layerA->Draw();
  th1d_south_layerA->Scale(1.0/nevents);
  th1d_south_layerA->SetMinimum(0);
  th1d_south_layerA->SetTitle(Form("FVTX South, all layers, %d GeV",energy));
  th1d_south_layerA->GetXaxis()->SetTitle("cluster #phi (rad)");
  th1d_south_layerA->GetYaxis()->SetTitle("counts/event");
  c1->Print(Form("FigsOther/cluster_phidist_fvtxs_energy%d.png",energy));
  c1->Print(Form("FigsOther/cluster_phidist_fvtxs_energy%d.pdf",energy));
  if ( energy == 62 ) th1d_south_layerA->SetMinimum(2.0);
  if ( energy == 62 ) th1d_south_layerA->SetMaximum(4.0);
  if ( energy == 39 ) th1d_south_layerA->SetMinimum(2.0);
  if ( energy == 39 ) th1d_south_layerA->SetMaximum(4.0);
  if ( energy == 20 ) th1d_south_layerA->SetMinimum(3.0);
  if ( energy == 20 ) th1d_south_layerA->SetMaximum(6.0);
  c1->Print(Form("FigsOther/cluster_phidist_fvtxs_energy%d_zoom.png",energy));
  c1->Print(Form("FigsOther/cluster_phidist_fvtxs_energy%d_zoom.pdf",energy));

  th1d_south_layer0->Draw();
  th1d_south_layer0->Scale(1.0/nevents);
  th1d_south_layer0->SetMinimum(0);
  th1d_south_layer0->SetTitle(Form("FVTX South, layer 0, %d GeV",energy));
  th1d_south_layer0->GetXaxis()->SetTitle("cluster #phi (rad)");
  th1d_south_layer0->GetYaxis()->SetTitle("counts/event");
  c1->Print(Form("FigsOther/cluster_phidist_fvtxs0_energy%d.png",energy));
  c1->Print(Form("FigsOther/cluster_phidist_fvtxs0_energy%d.pdf",energy));

  th1d_south_layer1->Draw();
  th1d_south_layer1->Scale(1.0/nevents);
  th1d_south_layer1->SetMinimum(0);
  th1d_south_layer1->SetTitle(Form("FVTX South, layer 1, %d GeV",energy));
  th1d_south_layer1->GetXaxis()->SetTitle("cluster #phi (rad)");
  th1d_south_layer1->GetYaxis()->SetTitle("counts/event");
  c1->Print(Form("FigsOther/cluster_phidist_fvtxs1_energy%d.png",energy));
  c1->Print(Form("FigsOther/cluster_phidist_fvtxs1_energy%d.pdf",energy));

  th1d_south_layer2->Draw();
  th1d_south_layer2->Scale(1.0/nevents);
  th1d_south_layer2->SetMinimum(0);
  th1d_south_layer2->SetTitle(Form("FVTX South, layer 2, %d GeV",energy));
  th1d_south_layer2->GetXaxis()->SetTitle("cluster #phi (rad)");
  th1d_south_layer2->GetYaxis()->SetTitle("counts/event");
  c1->Print(Form("FigsOther/cluster_phidist_fvtxs2_energy%d.png",energy));
  c1->Print(Form("FigsOther/cluster_phidist_fvtxs2_energy%d.pdf",energy));

  th1d_south_layer3->Draw();
  th1d_south_layer3->Scale(1.0/nevents);
  th1d_south_layer3->SetMinimum(0);
  th1d_south_layer3->SetTitle(Form("FVTX South, layer 3, %d GeV",energy));
  th1d_south_layer3->GetXaxis()->SetTitle("cluster #phi (rad)");
  th1d_south_layer3->GetYaxis()->SetTitle("counts/event");
  c1->Print(Form("FigsOther/cluster_phidist_fvtxs3_energy%d.png",energy));
  c1->Print(Form("FigsOther/cluster_phidist_fvtxs3_energy%d.pdf",energy));

  // ---

  TH1D* th1d_north_layerA = (TH1D*)file->Get("th1d_fvtxn_clus_phi");
  TH1D* th1d_north_layer0 = (TH1D*)file->Get("th1d_fvtxn0_clus_phi");
  TH1D* th1d_north_layer1 = (TH1D*)file->Get("th1d_fvtxn1_clus_phi");
  TH1D* th1d_north_layer2 = (TH1D*)file->Get("th1d_fvtxn2_clus_phi");
  TH1D* th1d_north_layer3 = (TH1D*)file->Get("th1d_fvtxn3_clus_phi");

  th1d_north_layerA->Draw();
  th1d_north_layerA->Scale(1.0/nevents);
  th1d_north_layerA->SetMinimum(0);
  th1d_north_layerA->SetTitle(Form("FVTX North, all layers, %d GeV",energy));
  th1d_north_layerA->GetXaxis()->SetTitle("cluster #phi (rad)");
  th1d_north_layerA->GetYaxis()->SetTitle("counts/event");
  c1->Print(Form("FigsOther/cluster_phidist_fvtxn_energy%d.png",energy));
  c1->Print(Form("FigsOther/cluster_phidist_fvtxn_energy%d.pdf",energy));

  th1d_north_layer0->Draw();
  th1d_north_layer0->Scale(1.0/nevents);
  th1d_north_layer0->SetMinimum(0);
  th1d_north_layer0->SetTitle(Form("FVTX North, layer 0, %d GeV",energy));
  th1d_north_layer0->GetXaxis()->SetTitle("cluster #phi (rad)");
  th1d_north_layer0->GetYaxis()->SetTitle("counts/event");
  c1->Print(Form("FigsOther/cluster_phidist_fvtxn0_energy%d.png",energy));
  c1->Print(Form("FigsOther/cluster_phidist_fvtxn0_energy%d.pdf",energy));

  th1d_north_layer1->Draw();
  th1d_north_layer1->Scale(1.0/nevents);
  th1d_north_layer1->SetMinimum(0);
  th1d_north_layer1->SetTitle(Form("FVTX North, layer 1, %d GeV",energy));
  th1d_north_layer1->GetXaxis()->SetTitle("cluster #phi (rad)");
  th1d_north_layer1->GetYaxis()->SetTitle("counts/event");
  c1->Print(Form("FigsOther/cluster_phidist_fvtxn1_energy%d.png",energy));
  c1->Print(Form("FigsOther/cluster_phidist_fvtxn1_energy%d.pdf",energy));

  th1d_north_layer2->Draw();
  th1d_north_layer2->Scale(1.0/nevents);
  th1d_north_layer2->SetMinimum(0);
  th1d_north_layer2->SetTitle(Form("FVTX North, layer 2, %d GeV",energy));
  th1d_north_layer2->GetXaxis()->SetTitle("cluster #phi (rad)");
  th1d_north_layer2->GetYaxis()->SetTitle("counts/event");
  c1->Print(Form("FigsOther/cluster_phidist_fvtxn2_energy%d.png",energy));
  c1->Print(Form("FigsOther/cluster_phidist_fvtxn2_energy%d.pdf",energy));

  th1d_north_layer3->Draw();
  th1d_north_layer3->Scale(1.0/nevents);
  th1d_north_layer3->SetMinimum(0);
  th1d_north_layer3->SetTitle(Form("FVTX North, layer 3, %d GeV",energy));
  th1d_north_layer3->GetXaxis()->SetTitle("cluster #phi (rad)");
  th1d_north_layer3->GetYaxis()->SetTitle("counts/event");
  c1->Print(Form("FigsOther/cluster_phidist_fvtxn3_energy%d.png",energy));
  c1->Print(Form("FigsOther/cluster_phidist_fvtxn3_energy%d.pdf",energy));

  // --- now bbc

  TH1D* th1d_south_ringA = (TH1D*)file->Get("th1d_bbc_charge_phi");
  TH1D* th1d_south_ring0 = (TH1D*)file->Get("th1d_bbc0_charge_phi");
  TH1D* th1d_south_ring1 = (TH1D*)file->Get("th1d_bbc1_charge_phi");
  TH1D* th1d_south_ring2 = (TH1D*)file->Get("th1d_bbc2_charge_phi");
  TH1D* th1d_south_ring3 = (TH1D*)file->Get("th1d_bbc3_charge_phi");
  TH1D* th1d_south_ring4 = (TH1D*)file->Get("th1d_bbc4_charge_phi");

  th1d_south_ringA->Draw();
  th1d_south_ringA->Scale(1.0/nevents);
  th1d_south_ringA->SetMinimum(0);
  th1d_south_ringA->SetTitle(Form("BBC South, all rings, %d GeV",energy));
  th1d_south_ringA->GetXaxis()->SetTitle("tube #phi (rad)");
  th1d_south_ringA->GetYaxis()->SetTitle("(counts*charge)/event");
  c1->Print(Form("FigsOther/charge_phidist_bbc_energy%d.png",energy));
  c1->Print(Form("FigsOther/charge_phidist_bbc_energy%d.pdf",energy));
  // if ( energy == 62 ) th1d_south_ringA->SetMinimum(5.0e8);
  // if ( energy == 62 ) th1d_south_ringA->SetMaximum(9.0e8);
  // if ( energy == 39 ) th1d_south_ringA->SetMinimum(3.0e8);
  // if ( energy == 39 ) th1d_south_ringA->SetMaximum(6.6e8);
  // if ( energy == 20 ) th1d_south_ringA->SetMinimum(6.0e7);
  // if ( energy == 20 ) th1d_south_ringA->SetMaximum(9.5e7);
  // c1->Print(Form("FigsOther/charge_phidist_bbc_energy%d_zoom.png",energy));
  // c1->Print(Form("FigsOther/charge_phidist_bbc_energy%d_zoom.pdf",energy));

  th1d_south_ring0->Draw();
  th1d_south_ring0->Scale(1.0/nevents);
  th1d_south_ring0->SetMinimum(0);
  th1d_south_ring0->SetTitle(Form("BBC South, ring 0, %d GeV",energy));
  th1d_south_ring0->GetXaxis()->SetTitle("tube #phi (rad)");
  th1d_south_ring0->GetYaxis()->SetTitle("(counts*charge)/event");
  c1->Print(Form("FigsOther/charge_phidist_bbc0_energy%d.png",energy));
  c1->Print(Form("FigsOther/charge_phidist_bbc0_energy%d.pdf",energy));

  th1d_south_ring1->Draw();
  th1d_south_ring1->Scale(1.0/nevents);
  th1d_south_ring1->SetMinimum(0);
  th1d_south_ring1->SetTitle(Form("BBC South, ring 1, %d GeV",energy));
  th1d_south_ring1->GetXaxis()->SetTitle("tube #phi (rad)");
  th1d_south_ring1->GetYaxis()->SetTitle("(counts*charge)/event");
  c1->Print(Form("FigsOther/charge_phidist_bbc1_energy%d.png",energy));
  c1->Print(Form("FigsOther/charge_phidist_bbc1_energy%d.pdf",energy));

  th1d_south_ring2->Draw();
  th1d_south_ring2->Scale(1.0/nevents);
  th1d_south_ring2->SetMinimum(0);
  th1d_south_ring2->SetTitle(Form("BBC South, ring 2, %d GeV",energy));
  th1d_south_ring2->GetXaxis()->SetTitle("tube #phi (rad)");
  th1d_south_ring2->GetYaxis()->SetTitle("(counts*charge)/event");
  c1->Print(Form("FigsOther/charge_phidist_bbc2_energy%d.png",energy));
  c1->Print(Form("FigsOther/charge_phidist_bbc2_energy%d.pdf",energy));

  th1d_south_ring3->Draw();
  th1d_south_ring3->Scale(1.0/nevents);
  th1d_south_ring3->SetMinimum(0);
  th1d_south_ring3->SetTitle(Form("BBC South, ring 3, %d GeV",energy));
  th1d_south_ring3->GetXaxis()->SetTitle("tube #phi (rad)");
  th1d_south_ring3->GetYaxis()->SetTitle("(counts*charge)/event");
  c1->Print(Form("FigsOther/charge_phidist_bbc3_energy%d.png",energy));
  c1->Print(Form("FigsOther/charge_phidist_bbc3_energy%d.pdf",energy));

  th1d_south_ring4->Draw();
  th1d_south_ring4->Scale(1.0/nevents);
  th1d_south_ring4->SetMinimum(0);
  th1d_south_ring4->SetTitle(Form("BBC South, ring 4, %d GeV",energy));
  th1d_south_ring4->GetXaxis()->SetTitle("tube #phi (rad)");
  th1d_south_ring4->GetYaxis()->SetTitle("(counts*charge)/event");
  c1->Print(Form("FigsOther/charge_phidist_bbc4_energy%d.png",energy));
  c1->Print(Form("FigsOther/charge_phidist_bbc4_energy%d.pdf",energy));

  // ---



  delete c1;

}
