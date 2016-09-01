void beforeafterflat_run(int);


void checkweights()
{

  beforeafterflat_run(455050); // 200
  beforeafterflat_run(455355); // 200
  beforeafterflat_run(456015); // 62
  beforeafterflat_run(456201); // 62
  beforeafterflat_run(457015); // 20
  beforeafterflat_run(457213); // 20
  beforeafterflat_run(458014); // 39
  beforeafterflat_run(458167); // 39

}


void beforeafterflat_run(int run)
{

  gStyle->SetOptTitle(1);

  TCanvas* c1 = new TCanvas("c1","");

  TFile* file = TFile::Open(Form("../output/hist_%d.root",run));
  if ( !file ) { cout << "no file" << endl; return; }
  TList* list = file->GetListOfKeys();
  if ( list->GetSize() < 1 ) { cout << "no keys" << endl; return; }

  TProfile* tp1f_bbc_before = (TProfile*)file->Get("tp1f_bbc_charge_tube");
  TProfile* tp1f_bbc_after = (TProfile*)file->Get("tp1f_bbc_charge_wtube");

  tp1f_bbc_before->SetLineColor(kRed);
  tp1f_bbc_after->SetLineColor(kBlack);
  tp1f_bbc_before->SetMarkerColor(kRed);
  tp1f_bbc_after->SetMarkerColor(kBlack);
  tp1f_bbc_before->SetMarkerStyle(kFullCircle);
  tp1f_bbc_after->SetMarkerStyle(kFullCircle);
  tp1f_bbc_before->Draw();
  tp1f_bbc_before->SetTitle("BBC South");
  tp1f_bbc_before->GetXaxis()->SetTitle("tube number");
  tp1f_bbc_before->GetYaxis()->SetTitle("tube charge");
  tp1f_bbc_after->Draw("same");

  c1->Print(Form("FigsWeight/check_bbc_run%d.png",run));
  c1->Print(Form("FigsWeight/check_bbc_run%d.pdf",run));

  TH1D* th1d_fvtxs_before = (TH1D*)file->Get("th1d_fvtxs_clus_phi");
  TH1D* th1d_fvtxs_after = (TH1D*)file->Get("th1d_fvtxs_clus_wphi");

  th1d_fvtxs_before->SetLineColor(kRed);
  th1d_fvtxs_after->SetLineColor(kBlack);
  th1d_fvtxs_before->Draw();
  th1d_fvtxs_before->SetTitle("FVTX South");
  th1d_fvtxs_before->GetXaxis()->SetTitle("cluster #phi");
  th1d_fvtxs_before->GetYaxis()->SetTitle("number of clusters");
  th1d_fvtxs_after->Draw("same");

  c1->Print(Form("FigsWeight/check_fvtxs_run%d.png",run));
  c1->Print(Form("FigsWeight/check_fvtxs_run%d.pdf",run));

  TH1D* th1d_fvtxs0_before = (TH1D*)file->Get("th1d_fvtxs0_clus_phi");
  TH1D* th1d_fvtxs0_after = (TH1D*)file->Get("th1d_fvtxs0_clus_wphi");

  th1d_fvtxs0_before->SetLineColor(kRed);
  th1d_fvtxs0_after->SetLineColor(kBlack);
  th1d_fvtxs0_before->Draw();
  th1d_fvtxs0_before->SetTitle("FVTX South Layer 0");
  th1d_fvtxs0_before->GetXaxis()->SetTitle("cluster #phi");
  th1d_fvtxs0_before->GetYaxis()->SetTitle("number of clusters");
  th1d_fvtxs0_after->Draw("same");

  c1->Print(Form("FigsWeight/check_fvtxs0_run%d.png",run));
  c1->Print(Form("FigsWeight/check_fvtxs0_run%d.pdf",run));

  TH1D* th1d_fvtxs1_before = (TH1D*)file->Get("th1d_fvtxs1_clus_phi");
  TH1D* th1d_fvtxs1_after = (TH1D*)file->Get("th1d_fvtxs1_clus_wphi");

  th1d_fvtxs1_before->SetLineColor(kRed);
  th1d_fvtxs1_after->SetLineColor(kBlack);
  th1d_fvtxs1_before->Draw();
  th1d_fvtxs1_before->SetTitle("FVTX South Layer 1");
  th1d_fvtxs1_before->GetXaxis()->SetTitle("cluster #phi");
  th1d_fvtxs1_before->GetYaxis()->SetTitle("number of clusters");
  th1d_fvtxs1_after->Draw("same");

  c1->Print(Form("FigsWeight/check_fvtxs1_run%d.png",run));
  c1->Print(Form("FigsWeight/check_fvtxs1_run%d.pdf",run));

  TH1D* th1d_fvtxs2_before = (TH1D*)file->Get("th1d_fvtxs2_clus_phi");
  TH1D* th1d_fvtxs2_after = (TH1D*)file->Get("th1d_fvtxs2_clus_wphi");

  th1d_fvtxs2_before->SetLineColor(kRed);
  th1d_fvtxs2_after->SetLineColor(kBlack);
  th1d_fvtxs2_before->Draw();
  th1d_fvtxs2_before->SetTitle("FVTX South Layer 2");
  th1d_fvtxs2_before->GetXaxis()->SetTitle("cluster #phi");
  th1d_fvtxs2_before->GetYaxis()->SetTitle("number of clusters");
  th1d_fvtxs2_after->Draw("same");

  c1->Print(Form("FigsWeight/check_fvtxs2_run%d.png",run));
  c1->Print(Form("FigsWeight/check_fvtxs2_run%d.pdf",run));

  TH1D* th1d_fvtxs3_before = (TH1D*)file->Get("th1d_fvtxs3_clus_phi");
  TH1D* th1d_fvtxs3_after = (TH1D*)file->Get("th1d_fvtxs3_clus_wphi");

  th1d_fvtxs3_before->SetLineColor(kRed);
  th1d_fvtxs3_after->SetLineColor(kBlack);
  th1d_fvtxs3_before->Draw();
  th1d_fvtxs3_before->SetTitle("FVTX South Layer 3");
  th1d_fvtxs3_before->GetXaxis()->SetTitle("cluster #phi");
  th1d_fvtxs3_before->GetYaxis()->SetTitle("number of clusters");
  th1d_fvtxs3_after->Draw("same");

  c1->Print(Form("FigsWeight/check_fvtxs3_run%d.png",run));
  c1->Print(Form("FigsWeight/check_fvtxs3_run%d.pdf",run));

  delete c1;

}

