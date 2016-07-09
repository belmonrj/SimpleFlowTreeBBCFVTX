void weight()
{

  TCanvas* c1 = new TCanvas("c1","");

  TFile* file = TFile::Open(Form("RootFiles/sum20_better.root"));
  TH1D* th1d_fvtxs_clus_phi_IR = (TH1D*)file->Get("th1d_fvtxs_clus_phi_IR");
  TH1D* th1d_fvtxs0_clus_phi_IR = (TH1D*)file->Get("th1d_fvtxs0_clus_phi_IR");
  TH1D* th1d_fvtxs1_clus_phi_IR = (TH1D*)file->Get("th1d_fvtxs1_clus_phi_IR");
  TH1D* th1d_fvtxs2_clus_phi_IR = (TH1D*)file->Get("th1d_fvtxs2_clus_phi_IR");
  TH1D* th1d_fvtxs3_clus_phi_IR = (TH1D*)file->Get("th1d_fvtxs3_clus_phi_IR");


  const int nbins = 50;
  if ( th1d_fvtxs_clus_phi_IR->GetNbinsX() != nbins )
    {
      cout << "YOU'RE GONNA DIE" << endl;
      return;
    }

  double ave = th1d_fvtxs_clus_phi_IR->Integral(1,nbins); // use 1 and nbins to exclude underflow (0) and overflow (nbins+1)
  double ave0 = th1d_fvtxs0_clus_phi_IR->Integral(1,nbins);
  double ave1 = th1d_fvtxs1_clus_phi_IR->Integral(1,nbins);
  double ave2 = th1d_fvtxs2_clus_phi_IR->Integral(1,nbins);
  double ave3 = th1d_fvtxs3_clus_phi_IR->Integral(1,nbins);
  ave /= nbins;
  ave0 /= nbins;
  ave1 /= nbins;
  ave2 /= nbins;
  ave3 /= nbins;
  double pphi[nbins];
  double pweight[nbins];
  double pweight0[nbins];
  double pweight1[nbins];
  double pweight2[nbins];
  double pweight3[nbins];
  ofstream fout((const char*)Form("weight.dat"));
  for ( int i = 0; i < nbins; ++i )
    {
      double phi = th1d_fvtxs_clus_phi_IR->GetBinCenter(i+1);
      double weight = ave/th1d_fvtxs_clus_phi_IR->GetBinContent(i+1);
      double weight0 = ave0/th1d_fvtxs0_clus_phi_IR->GetBinContent(i+1);
      double weight1 = ave1/th1d_fvtxs1_clus_phi_IR->GetBinContent(i+1);
      double weight2 = ave2/th1d_fvtxs2_clus_phi_IR->GetBinContent(i+1);
      double weight3 = ave3/th1d_fvtxs3_clus_phi_IR->GetBinContent(i+1);
      if ( !TMath::Finite(weight) ) weight = 0;
      if ( !TMath::Finite(weight0) ) weight0 = 0;
      if ( !TMath::Finite(weight1) ) weight1 = 0;
      if ( !TMath::Finite(weight2) ) weight2 = 0;
      if ( !TMath::Finite(weight3) ) weight3 = 0;
      pphi[i] = phi;
      pweight[i] = weight;
      pweight0[i] = weight0;
      pweight1[i] = weight1;
      pweight2[i] = weight2;
      pweight3[i] = weight3;
      fout << phi << " "
           << weight << " "
           << weight0 << " "
           << weight1 << " "
           << weight2 << " "
           << weight3 << endl;
    }
  fout.close();

  TGraph* tg_weight = new TGraph(nbins,pphi,pweight);
  tg_weight->SetMarkerStyle(kFullCircle);
  tg_weight->Draw("ap");
  tg_weight->SetMinimum(0.0);
  tg_weight->SetMaximum(4.0);
  tg_weight->GetXaxis()->SetLimits(-TMath::Pi(),TMath::Pi());
  c1->Print(Form("weight.png"));

  TGraph* tg_weight0 = new TGraph(nbins,pphi,pweight0);
  tg_weight0->SetMarkerStyle(kFullCircle);
  tg_weight0->Draw("ap");
  tg_weight0->SetMinimum(0.0);
  tg_weight0->SetMaximum(4.0);
  tg_weight0->GetXaxis()->SetLimits(-TMath::Pi(),TMath::Pi());
  c1->Print(Form("weight0.png"));

  TGraph* tg_weight1 = new TGraph(nbins,pphi,pweight1);
  tg_weight1->SetMarkerStyle(kFullCircle);
  tg_weight1->Draw("ap");
  tg_weight1->SetMinimum(0.0);
  tg_weight1->SetMaximum(4.0);
  tg_weight1->GetXaxis()->SetLimits(-TMath::Pi(),TMath::Pi());
  c1->Print(Form("weight1.png"));

  TGraph* tg_weight2 = new TGraph(nbins,pphi,pweight2);
  tg_weight2->SetMarkerStyle(kFullCircle);
  tg_weight2->Draw("ap");
  tg_weight2->SetMinimum(0.0);
  tg_weight2->SetMaximum(4.0);
  tg_weight2->GetXaxis()->SetLimits(-TMath::Pi(),TMath::Pi());
  c1->Print(Form("weight2.png"));

  TGraph* tg_weight3 = new TGraph(nbins,pphi,pweight3);
  tg_weight3->SetMarkerStyle(kFullCircle);
  tg_weight3->Draw("ap");
  tg_weight3->SetMinimum(0.0);
  tg_weight3->SetMaximum(4.0);
  tg_weight3->GetXaxis()->SetLimits(-TMath::Pi(),TMath::Pi());
  c1->Print(Form("weight3.png"));

  fout.open((const char*)Form("fweight.dat"));
  for ( int i = 0; i < nbins; ++i )
    {
      double phi = th1d_fvtxs_clus_phi_IR->GetBinCenter(i+1);
      double weight = ave/th1d_fvtxs_clus_phi_IR->GetBinContent(i+1);
      double weight0 = ave0/th1d_fvtxs0_clus_phi_IR->GetBinContent(i+1);
      double weight1 = ave1/th1d_fvtxs1_clus_phi_IR->GetBinContent(i+1);
      double weight2 = ave2/th1d_fvtxs2_clus_phi_IR->GetBinContent(i+1);
      double weight3 = ave3/th1d_fvtxs3_clus_phi_IR->GetBinContent(i+1);
      if ( weight > 1.5 || weight < 0.5 ) weight = 0;
      if ( weight0 > 1.5 || weight0 < 0.5 ) weight0 = 0;
      if ( weight1 > 1.5 || weight1 < 0.5 ) weight1 = 0;
      if ( weight2 > 1.5 || weight2 < 0.5 ) weight2 = 0;
      if ( weight3 > 1.5 || weight3 < 0.5 ) weight3 = 0;
      pphi[i] = phi;
      pweight[i] = weight;
      pweight0[i] = weight0;
      pweight1[i] = weight1;
      pweight2[i] = weight2;
      pweight3[i] = weight3;
      fout << phi << " "
           << weight << " "
           << weight0 << " "
           << weight1 << " "
           << weight2 << " "
           << weight3 << endl;
    }
  fout.close();

  TGraph* tg_fweight = new TGraph(nbins,pphi,pweight);
  tg_fweight->SetMarkerStyle(kFullCircle);
  tg_fweight->Draw("ap");
  tg_fweight->SetMinimum(0.0);
  tg_fweight->SetMaximum(2.0);
  tg_fweight->GetXaxis()->SetLimits(-TMath::Pi(),TMath::Pi());
  c1->Print(Form("fweight.png"));

  TGraph* tg_fweight0 = new TGraph(nbins,pphi,pweight0);
  tg_fweight0->SetMarkerStyle(kFullCircle);
  tg_fweight0->Draw("ap");
  tg_fweight0->SetMinimum(0.0);
  tg_fweight0->SetMaximum(2.0);
  tg_fweight0->GetXaxis()->SetLimits(-TMath::Pi(),TMath::Pi());
  c1->Print(Form("fweight0.png"));

  TGraph* tg_fweight1 = new TGraph(nbins,pphi,pweight1);
  tg_fweight1->SetMarkerStyle(kFullCircle);
  tg_fweight1->Draw("ap");
  tg_fweight1->SetMinimum(0.0);
  tg_fweight1->SetMaximum(2.0);
  tg_fweight1->GetXaxis()->SetLimits(-TMath::Pi(),TMath::Pi());
  c1->Print(Form("fweight1.png"));

  TGraph* tg_fweight2 = new TGraph(nbins,pphi,pweight2);
  tg_fweight2->SetMarkerStyle(kFullCircle);
  tg_fweight2->Draw("ap");
  tg_fweight2->SetMinimum(0.0);
  tg_fweight2->SetMaximum(2.0);
  tg_fweight2->GetXaxis()->SetLimits(-TMath::Pi(),TMath::Pi());
  c1->Print(Form("fweight2.png"));

  TGraph* tg_fweight3 = new TGraph(nbins,pphi,pweight3);
  tg_fweight3->SetMarkerStyle(kFullCircle);
  tg_fweight3->Draw("ap");
  tg_fweight3->SetMinimum(0.0);
  tg_fweight3->SetMaximum(2.0);
  tg_fweight3->GetXaxis()->SetLimits(-TMath::Pi(),TMath::Pi());
  c1->Print(Form("fweight3.png"));

}
