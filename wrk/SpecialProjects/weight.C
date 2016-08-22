void runweight2d(int);
void only2d(int);


void weight()
{

  //runweight2d(456652); // remember when this was needed?
  only2d(455355); // now this is needed, and that could change again in a little while (must be some weird environment thing)
  //runweight2d(455355);

  //  return;

  int run;
  ifstream fin;

  fin.open("list_20.short");
  while ( fin >> run ) only2d(run);
  fin.close();

  // fin.open("list_39.short");
  // while ( fin >> run ) only2d(run);
  // fin.close();

  // fin.open("list_62.short");
  // while ( fin >> run ) only2d(run);
  // fin.close();

  // fin.open("list_200.short");
  // while ( fin >> run ) only2d(run);
  // fin.close();

  // fin.open("list_20.short");
  // while ( fin >> run ) runweight2d(run);
  // fin.close();

  // fin.open("list_39.short");
  // while ( fin >> run ) runweight2d(run);
  // fin.close();

  // fin.open("list_62.short");
  // while ( fin >> run ) runweight2d(run);
  // fin.close();

  // fin.open("list_200.short");
  // while ( fin >> run ) runweight2d(run);
  // fin.close();

}


void only2d(int run)
{

  cout << "Now processing run " << run << endl;

  TCanvas* c1 = new TCanvas("c1","");

  TFile* file = TFile::Open(Form("RootFiles/svrb_run%d_pass0.root",run));
  if ( !file )
    {
      cout << "WARNING: file does not exist for run " << run << endl;
      return;
    }

  // ------------- //
  // --- south --- //
  // ------------- //

  TH2D* th2d_fvtxs_clus_phi = (TH2D*)file->Get("th2d_fvtxs_clus_phi_GC");
  TH2D* th2d_fvtxs0_clus_phi = (TH2D*)file->Get("th2d_fvtxs0_clus_phi_GC");
  TH2D* th2d_fvtxs1_clus_phi = (TH2D*)file->Get("th2d_fvtxs1_clus_phi_GC");
  TH2D* th2d_fvtxs2_clus_phi = (TH2D*)file->Get("th2d_fvtxs2_clus_phi_GC");
  TH2D* th2d_fvtxs3_clus_phi = (TH2D*)file->Get("th2d_fvtxs3_clus_phi_GC");

  bool ihavethehistos = th2d_fvtxs_clus_phi && th2d_fvtxs0_clus_phi && th2d_fvtxs1_clus_phi && th2d_fvtxs2_clus_phi && th2d_fvtxs3_clus_phi;

  if ( !ihavethehistos )
    {
      cout << "YOU'RE GONNA DIE (missing histograms)" << endl;
      return;
    }

  const int nbinsx = 20;
  if ( th2d_fvtxs_clus_phi->GetNbinsX() != nbinsx )
    {
      cout << "YOU'RE GONNA DIE " << nbinsx << " " << th2d_fvtxs_clus_phi->GetNbinsX() << endl;
      return;
    }
  const int nbinsy = 50;
  if ( th2d_fvtxs_clus_phi->GetNbinsY() != nbinsy )
    {
      cout << "YOU'RE GONNA DIE " << nbinsy << " " << th2d_fvtxs_clus_phi->GetNbinsX() << endl;
      return;
    }

  double integral = th2d_fvtxs_clus_phi->Integral(1,nbinsx,1,nbinsy);
  double totalbins = (double)nbinsx*(double)nbinsy;
  double average = integral/totalbins;

  cout << "average is " << average << endl;

  TFile* fout = TFile::Open(Form("WeightFiles/only2d_run%d.root",run),"recreate");

  TH2D* th2d_fvtxs_clus_phi_weight = (TH2D*)th2d_fvtxs_clus_phi->Clone("th2d_fvtxs_clus_phi_weight");
  TH2D* th2d_fvtxs0_clus_phi_weight = (TH2D*)th2d_fvtxs0_clus_phi->Clone("th2d_fvtxs0_clus_phi_weight");
  TH2D* th2d_fvtxs1_clus_phi_weight = (TH2D*)th2d_fvtxs1_clus_phi->Clone("th2d_fvtxs1_clus_phi_weight");
  TH2D* th2d_fvtxs2_clus_phi_weight = (TH2D*)th2d_fvtxs2_clus_phi->Clone("th2d_fvtxs2_clus_phi_weight");
  TH2D* th2d_fvtxs3_clus_phi_weight = (TH2D*)th2d_fvtxs3_clus_phi->Clone("th2d_fvtxs3_clus_phi_weight");
  for ( int i = 0; i < nbinsx; ++i )
    {
      for ( int j = 0; j < nbinsy; ++j )
        {
          double bincontent = th2d_fvtxs_clus_phi->GetBinContent(i+1,j+1);
          double newvalue = 0;
          double temp = average/bincontent;
          if ( temp == temp ) newvalue = temp;
          th2d_fvtxs_clus_phi_weight->SetBinContent(i+1,j+1,newvalue);
          // ---
          bincontent = th2d_fvtxs0_clus_phi->GetBinContent(i+1,j+1);
          newvalue = 0;
          temp = average/bincontent;
          if ( temp == temp ) newvalue = temp;
          th2d_fvtxs0_clus_phi_weight->SetBinContent(i+1,j+1,newvalue);
          // ---
          bincontent = th2d_fvtxs1_clus_phi->GetBinContent(i+1,j+1);
          newvalue = 0;
          temp = average/bincontent;
          if ( temp == temp ) newvalue = temp;
          th2d_fvtxs1_clus_phi_weight->SetBinContent(i+1,j+1,newvalue);
          // ---
          bincontent = th2d_fvtxs2_clus_phi->GetBinContent(i+1,j+1);
          newvalue = 0;
          temp = average/bincontent;
          if ( temp == temp ) newvalue = temp;
          th2d_fvtxs2_clus_phi_weight->SetBinContent(i+1,j+1,newvalue);
          // ---
          bincontent = th2d_fvtxs3_clus_phi->GetBinContent(i+1,j+1);
          newvalue = 0;
          temp = average/bincontent;
          if ( temp == temp ) newvalue = temp;
          th2d_fvtxs3_clus_phi_weight->SetBinContent(i+1,j+1,newvalue);
        }
    }

  th2d_fvtxs_clus_phi->Draw("colz");
  c1->Print(Form("FigsWeight/testing_only2ddist_%d.png",run));
  th2d_fvtxs_clus_phi_weight->Draw("colz");
  c1->Print(Form("FigsWeight/testing_only2dweight_%d.png",run));

  th2d_fvtxs0_clus_phi->Draw("colz");
  c1->Print(Form("FigsWeight/only2ddist_fvtxs0_%d.png",run));
  th2d_fvtxs0_clus_phi_weight->Draw("colz");
  c1->Print(Form("FigsWeight/only2dweight_fvtxs0_%d.png",run));

  th2d_fvtxs1_clus_phi->Draw("colz");
  c1->Print(Form("FigsWeight/only2ddist_fvtxs1_%d.png",run));
  th2d_fvtxs1_clus_phi_weight->Draw("colz");
  c1->Print(Form("FigsWeight/only2dweight_fvtxs1_%d.png",run));

  th2d_fvtxs2_clus_phi->Draw("colz");
  c1->Print(Form("FigsWeight/only2ddist_fvtxs2_%d.png",run));
  th2d_fvtxs2_clus_phi_weight->Draw("colz");
  c1->Print(Form("FigsWeight/only2dweight_fvtxs2_%d.png",run));

  th2d_fvtxs3_clus_phi->Draw("colz");
  c1->Print(Form("FigsWeight/only2ddist_fvtxs3_%d.png",run));
  th2d_fvtxs3_clus_phi_weight->Draw("colz");
  c1->Print(Form("FigsWeight/only2dweight_fvtxs3_%d.png",run));

  fout->Write();
  fout->Close();

  delete c1;

} // only2d



void runweight2d(int run)
{

  cout << "Now processing run " << run << endl;

  TCanvas* c1 = new TCanvas("c1","");

  cout << "Attempting to open the file" << endl;

  TFile* file = TFile::Open(Form("RootFiles/svrb_run%d_pass0.root",run));
  if ( !file )
    {
      cout << "WARNING: file does not exist for run " << run << endl;
      return;
    }

  // ------------- //
  // --- south --- //
  // ------------- //

  TH2D* th2d_fvtxs_clus_phi_GC = (TH2D*)file->Get("th2d_fvtxs_clus_phi_GC");
  TH2D* th2d_fvtxs0_clus_phi_GC = (TH2D*)file->Get("th2d_fvtxs0_clus_phi_GC");
  TH2D* th2d_fvtxs1_clus_phi_GC = (TH2D*)file->Get("th2d_fvtxs1_clus_phi_GC");
  TH2D* th2d_fvtxs2_clus_phi_GC = (TH2D*)file->Get("th2d_fvtxs2_clus_phi_GC");
  TH2D* th2d_fvtxs3_clus_phi_GC = (TH2D*)file->Get("th2d_fvtxs3_clus_phi_GC");

  const int nbinsx = 20;
  if ( th2d_fvtxs_clus_phi_GC->GetNbinsX() != nbinsx )
    {
      cout << "YOU'RE GONNA DIE " << nbinsx << " " << th2d_fvtxs_clus_phi_GC->GetNbinsX() << endl;
      return;
    }
  const int nbinsy = 50;
  if ( th2d_fvtxs_clus_phi_GC->GetNbinsY() != nbinsy )
    {
      cout << "YOU'RE GONNA DIE " << nbinsy << " " << th2d_fvtxs_clus_phi_GC->GetNbinsX() << endl;
      return;
    }

  TFile* fout = TFile::Open(Form("WeightFiles/weight2d_run%d.root",run),"recreate");

  // -- loop over z bins
  for ( int i = 0; i < nbinsx; ++i )
    {
      // --- get the 1d projections
      TH1D* th1d_fvtxs_zvtxbin_clus_phi_GC = (TH1D*)th2d_fvtxs_clus_phi_GC->ProjectionY(Form("th1d_fvtxs_zvtx%d_clus_phi_GC",i),i+1,i+1);
      th1d_fvtxs_zvtxbin_clus_phi_GC->Draw();
      c1->Print(Form("FigsWeight/dist2d_fvtxs_zvtx%d_clus_phi_GC_run%d.png",i,run));

      TH1D* th1d_fvtxs0_zvtxbin_clus_phi_GC = (TH1D*)th2d_fvtxs0_clus_phi_GC->ProjectionY(Form("th1d_fvtxs0_zvtx%d_clus_phi_GC",i),i+1,i+1);
      th1d_fvtxs0_zvtxbin_clus_phi_GC->Draw();
      c1->Print(Form("FigsWeight/dist2d_fvtxs0_zvtx%d_clus_phi_GC_run%d.png",i,run));

      TH1D* th1d_fvtxs1_zvtxbin_clus_phi_GC = (TH1D*)th2d_fvtxs1_clus_phi_GC->ProjectionY(Form("th1d_fvtxs1_zvtx%d_clus_phi_GC",i),i+1,i+1);
      th1d_fvtxs1_zvtxbin_clus_phi_GC->Draw();
      c1->Print(Form("FigsWeight/dist2d_fvtxs1_zvtx%d_clus_phi_GC_run%d.png",i,run));

      TH1D* th1d_fvtxs2_zvtxbin_clus_phi_GC = (TH1D*)th2d_fvtxs2_clus_phi_GC->ProjectionY(Form("th1d_fvtxs2_zvtx%d_clus_phi_GC",i),i+1,i+1);
      th1d_fvtxs2_zvtxbin_clus_phi_GC->Draw();
      c1->Print(Form("FigsWeight/dist2d_fvtxs2_zvtx%d_clus_phi_GC_run%d.png",i,run));

      TH1D* th1d_fvtxs3_zvtxbin_clus_phi_GC = (TH1D*)th2d_fvtxs3_clus_phi_GC->ProjectionY(Form("th1d_fvtxs3_zvtx%d_clus_phi_GC",i),i+1,i+1);
      th1d_fvtxs3_zvtxbin_clus_phi_GC->Draw();
      c1->Print(Form("FigsWeight/dist2d_fvtxs3_zvtx%d_clus_phi_GC_run%d.png",i,run));

      // --- make 1d clones
      TH1D* th1d_weight_fvtxs_zvtxbin_clus_phi_GC = (TH1D*)th1d_fvtxs_zvtxbin_clus_phi_GC->Clone(Form("th1d_weight_fvtxs_zvtx%d_clus_phi_GC",i));
      TH1D* th1d_weight_fvtxs0_zvtxbin_clus_phi_GC = (TH1D*)th1d_fvtxs0_zvtxbin_clus_phi_GC->Clone(Form("th1d_weight_fvtxs0_zvtx%d_clus_phi_GC",i));
      TH1D* th1d_weight_fvtxs1_zvtxbin_clus_phi_GC = (TH1D*)th1d_fvtxs1_zvtxbin_clus_phi_GC->Clone(Form("th1d_weight_fvtxs1_zvtx%d_clus_phi_GC",i));
      TH1D* th1d_weight_fvtxs2_zvtxbin_clus_phi_GC = (TH1D*)th1d_fvtxs2_zvtxbin_clus_phi_GC->Clone(Form("th1d_weight_fvtxs2_zvtx%d_clus_phi_GC",i));
      TH1D* th1d_weight_fvtxs3_zvtxbin_clus_phi_GC = (TH1D*)th1d_fvtxs3_zvtxbin_clus_phi_GC->Clone(Form("th1d_weight_fvtxs3_zvtx%d_clus_phi_GC",i));
      double ave = th1d_fvtxs_zvtxbin_clus_phi_GC->Integral(1,nbinsy); // use 1 and nbins to exclude underflow (0) and overflow (nbins+1)
      double ave0 = th1d_fvtxs0_zvtxbin_clus_phi_GC->Integral(1,nbinsy);
      double ave1 = th1d_fvtxs1_zvtxbin_clus_phi_GC->Integral(1,nbinsy);
      double ave2 = th1d_fvtxs2_zvtxbin_clus_phi_GC->Integral(1,nbinsy);
      double ave3 = th1d_fvtxs3_zvtxbin_clus_phi_GC->Integral(1,nbinsy);
      ave /= nbinsy;
      ave0 /= nbinsy;
      ave1 /= nbinsy;
      ave2 /= nbinsy;
      ave3 /= nbinsy;
      for ( int j = 0; j < nbinsy; ++j )
	{
	  double phi = th1d_fvtxs_zvtxbin_clus_phi_GC->GetBinCenter(j+1);
	  double weight = ave/th1d_fvtxs_zvtxbin_clus_phi_GC->GetBinContent(j+1);
	  double weight0 = ave0/th1d_fvtxs0_zvtxbin_clus_phi_GC->GetBinContent(j+1);
	  double weight1 = ave1/th1d_fvtxs1_zvtxbin_clus_phi_GC->GetBinContent(j+1);
	  double weight2 = ave2/th1d_fvtxs2_zvtxbin_clus_phi_GC->GetBinContent(j+1);
	  double weight3 = ave3/th1d_fvtxs3_zvtxbin_clus_phi_GC->GetBinContent(j+1);
	  if ( !TMath::Finite(weight) ) weight = 0;
	  if ( !TMath::Finite(weight0) ) weight0 = 0;
	  if ( !TMath::Finite(weight1) ) weight1 = 0;
	  if ( !TMath::Finite(weight2) ) weight2 = 0;
	  if ( !TMath::Finite(weight3) ) weight3 = 0;
	  if ( weight  > 1.5 || weight  < 0.5 ) weight = 0;
	  if ( weight0 > 1.5 || weight0 < 0.5 ) weight0 = 0;
	  if ( weight1 > 1.5 || weight1 < 0.5 ) weight1 = 0;
	  if ( weight2 > 1.5 || weight2 < 0.5 ) weight2 = 0;
	  if ( weight3 > 1.5 || weight3 < 0.5 ) weight3 = 0;
	  // ---
	  th1d_weight_fvtxs_zvtxbin_clus_phi_GC->SetBinContent(j+1,weight);
	  th1d_weight_fvtxs0_zvtxbin_clus_phi_GC->SetBinContent(j+1,weight0);
	  th1d_weight_fvtxs1_zvtxbin_clus_phi_GC->SetBinContent(j+1,weight1);
	  th1d_weight_fvtxs2_zvtxbin_clus_phi_GC->SetBinContent(j+1,weight2);
	  th1d_weight_fvtxs3_zvtxbin_clus_phi_GC->SetBinContent(j+1,weight3);
	}

      th1d_weight_fvtxs_zvtxbin_clus_phi_GC->Draw();
      c1->Print(Form("FigsWeight/weight2d_fvtxs_zvtx%d_clus_phi_GC_run%d.png",i,run));
      th1d_fvtxs_zvtxbin_clus_phi_GC->Multiply(th1d_weight_fvtxs_zvtxbin_clus_phi_GC);
      th1d_fvtxs_zvtxbin_clus_phi_GC->Draw();
      c1->Print(Form("FigsWeight/check2d_fvtxs_zvtx%d_clus_phi_GC_run%d.png",i,run));

      th1d_weight_fvtxs0_zvtxbin_clus_phi_GC->Draw();
      c1->Print(Form("FigsWeight/weight2d_fvtxs0_zvtx%d_clus_phi_GC_run%d.png",i,run));
      th1d_fvtxs0_zvtxbin_clus_phi_GC->Multiply(th1d_weight_fvtxs0_zvtxbin_clus_phi_GC);
      th1d_fvtxs0_zvtxbin_clus_phi_GC->Draw();
      c1->Print(Form("FigsWeight/check2d_fvtxs0_zvtx%d_clus_phi_GC_run%d.png",i,run));

      th1d_weight_fvtxs1_zvtxbin_clus_phi_GC->Draw();
      c1->Print(Form("FigsWeight/weight2d_fvtxs1_zvtx%d_clus_phi_GC_run%d.png",i,run));
      th1d_fvtxs1_zvtxbin_clus_phi_GC->Multiply(th1d_weight_fvtxs1_zvtxbin_clus_phi_GC);
      th1d_fvtxs1_zvtxbin_clus_phi_GC->Draw();
      c1->Print(Form("FigsWeight/check2d_fvtxs1_zvtx%d_clus_phi_GC_run%d.png",i,run));

      th1d_weight_fvtxs2_zvtxbin_clus_phi_GC->Draw();
      c1->Print(Form("FigsWeight/weight2d_fvtxs2_zvtx%d_clus_phi_GC_run%d.png",i,run));
      th1d_fvtxs2_zvtxbin_clus_phi_GC->Multiply(th1d_weight_fvtxs2_zvtxbin_clus_phi_GC);
      th1d_fvtxs2_zvtxbin_clus_phi_GC->Draw();
      c1->Print(Form("FigsWeight/check2d_fvtxs2_zvtx%d_clus_phi_GC_run%d.png",i,run));

      th1d_weight_fvtxs3_zvtxbin_clus_phi_GC->Draw();
      c1->Print(Form("FigsWeight/weight2d_fvtxs3_zvtx%d_clus_phi_GC_run%d.png",i,run));
      th1d_fvtxs3_zvtxbin_clus_phi_GC->Multiply(th1d_weight_fvtxs3_zvtxbin_clus_phi_GC);
      th1d_fvtxs3_zvtxbin_clus_phi_GC->Draw();
      c1->Print(Form("FigsWeight/check2d_fvtxs3_zvtx%d_clus_phi_GC_run%d.png",i,run));

    }

  // ------------- //
  // --- north --- //
  // ------------- //

  TH2D* th2d_fvtxn_clus_phi_GC = (TH2D*)file->Get("th2d_fvtxn_clus_phi_GC");
  TH2D* th2d_fvtxn0_clus_phi_GC = (TH2D*)file->Get("th2d_fvtxn0_clus_phi_GC");
  TH2D* th2d_fvtxn1_clus_phi_GC = (TH2D*)file->Get("th2d_fvtxn1_clus_phi_GC");
  TH2D* th2d_fvtxn2_clus_phi_GC = (TH2D*)file->Get("th2d_fvtxn2_clus_phi_GC");
  TH2D* th2d_fvtxn3_clus_phi_GC = (TH2D*)file->Get("th2d_fvtxn3_clus_phi_GC");

  const int nbinsx = 20;
  if ( th2d_fvtxn_clus_phi_GC->GetNbinsX() != nbinsx )
    {
      cout << "YOU'RE GONNA DIE " << nbinsx << " " << th2d_fvtxn_clus_phi_GC->GetNbinsX() << endl;
      return;
    }
  const int nbinsy = 50;
  if ( th2d_fvtxn_clus_phi_GC->GetNbinsY() != nbinsy )
    {
      cout << "YOU'RE GONNA DIE " << nbinsy << " " << th2d_fvtxn_clus_phi_GC->GetNbinsX() << endl;
      return;
    }

  // TFile* fout = TFile::Open(Form("WeightFiles/weight2d_run%d.root",run),"recreate");

  // -- loop over z bins
  for ( int i = 0; i < nbinsx; ++i )
    {
      // --- get the 1d projections
      TH1D* th1d_fvtxn_zvtxbin_clus_phi_GC = (TH1D*)th2d_fvtxn_clus_phi_GC->ProjectionY(Form("th1d_fvtxn_zvtx%d_clus_phi_GC",i),i+1,i+1);
      th1d_fvtxn_zvtxbin_clus_phi_GC->Draw();
      c1->Print(Form("FigsWeight/weight2d_fvtxn_zvtx%d_clus_phi_GC_run%d.png",i,run));

      TH1D* th1d_fvtxn0_zvtxbin_clus_phi_GC = (TH1D*)th2d_fvtxn0_clus_phi_GC->ProjectionY(Form("th1d_fvtxn0_zvtx%d_clus_phi_GC",i),i+1,i+1);
      th1d_fvtxn0_zvtxbin_clus_phi_GC->Draw();
      c1->Print(Form("FigsWeight/weight2d_fvtxn0_zvtx%d_clus_phi_GC_run%d.png",i,run));

      TH1D* th1d_fvtxn1_zvtxbin_clus_phi_GC = (TH1D*)th2d_fvtxn1_clus_phi_GC->ProjectionY(Form("th1d_fvtxn1_zvtx%d_clus_phi_GC",i),i+1,i+1);
      th1d_fvtxn1_zvtxbin_clus_phi_GC->Draw();
      c1->Print(Form("FigsWeight/weight2d_fvtxn1_zvtx%d_clus_phi_GC_run%d.png",i,run));

      TH1D* th1d_fvtxn2_zvtxbin_clus_phi_GC = (TH1D*)th2d_fvtxn2_clus_phi_GC->ProjectionY(Form("th1d_fvtxn2_zvtx%d_clus_phi_GC",i),i+1,i+1);
      th1d_fvtxn2_zvtxbin_clus_phi_GC->Draw();
      c1->Print(Form("FigsWeight/weight2d_fvtxn2_zvtx%d_clus_phi_GC_run%d.png",i,run));

      TH1D* th1d_fvtxn3_zvtxbin_clus_phi_GC = (TH1D*)th2d_fvtxn3_clus_phi_GC->ProjectionY(Form("th1d_fvtxn3_zvtx%d_clus_phi_GC",i),i+1,i+1);
      th1d_fvtxn3_zvtxbin_clus_phi_GC->Draw();
      c1->Print(Form("FigsWeight/weight2d_fvtxn3_zvtx%d_clus_phi_GC_run%d.png",i,run));

      // --- make 1d clones
      TH1D* th1d_weight_fvtxn_zvtxbin_clus_phi_GC = (TH1D*)th1d_fvtxn_zvtxbin_clus_phi_GC->Clone(Form("th1d_weight_fvtxn_zvtx%d_clus_phi_GC",i));
      TH1D* th1d_weight_fvtxn0_zvtxbin_clus_phi_GC = (TH1D*)th1d_fvtxn0_zvtxbin_clus_phi_GC->Clone(Form("th1d_weight_fvtxn0_zvtx%d_clus_phi_GC",i));
      TH1D* th1d_weight_fvtxn1_zvtxbin_clus_phi_GC = (TH1D*)th1d_fvtxn1_zvtxbin_clus_phi_GC->Clone(Form("th1d_weight_fvtxn1_zvtx%d_clus_phi_GC",i));
      TH1D* th1d_weight_fvtxn2_zvtxbin_clus_phi_GC = (TH1D*)th1d_fvtxn2_zvtxbin_clus_phi_GC->Clone(Form("th1d_weight_fvtxn2_zvtx%d_clus_phi_GC",i));
      TH1D* th1d_weight_fvtxn3_zvtxbin_clus_phi_GC = (TH1D*)th1d_fvtxn3_zvtxbin_clus_phi_GC->Clone(Form("th1d_weight_fvtxn3_zvtx%d_clus_phi_GC",i));
      double ave = th1d_fvtxn_zvtxbin_clus_phi_GC->Integral(1,nbinsy); // use 1 and nbins to exclude underflow (0) and overflow (nbins+1)
      double ave0 = th1d_fvtxn0_zvtxbin_clus_phi_GC->Integral(1,nbinsy);
      double ave1 = th1d_fvtxn1_zvtxbin_clus_phi_GC->Integral(1,nbinsy);
      double ave2 = th1d_fvtxn2_zvtxbin_clus_phi_GC->Integral(1,nbinsy);
      double ave3 = th1d_fvtxn3_zvtxbin_clus_phi_GC->Integral(1,nbinsy);
      ave /= nbinsy;
      ave0 /= nbinsy;
      ave1 /= nbinsy;
      ave2 /= nbinsy;
      ave3 /= nbinsy;
      for ( int j = 0; j < nbinsy; ++j )
	{
	  double phi = th1d_fvtxn_zvtxbin_clus_phi_GC->GetBinCenter(j+1);
	  double weight = ave/th1d_fvtxn_zvtxbin_clus_phi_GC->GetBinContent(j+1);
	  double weight0 = ave0/th1d_fvtxn0_zvtxbin_clus_phi_GC->GetBinContent(j+1);
	  double weight1 = ave1/th1d_fvtxn1_zvtxbin_clus_phi_GC->GetBinContent(j+1);
	  double weight2 = ave2/th1d_fvtxn2_zvtxbin_clus_phi_GC->GetBinContent(j+1);
	  double weight3 = ave3/th1d_fvtxn3_zvtxbin_clus_phi_GC->GetBinContent(j+1);
	  if ( !TMath::Finite(weight) ) weight = 0;
	  if ( !TMath::Finite(weight0) ) weight0 = 0;
	  if ( !TMath::Finite(weight1) ) weight1 = 0;
	  if ( !TMath::Finite(weight2) ) weight2 = 0;
	  if ( !TMath::Finite(weight3) ) weight3 = 0;
	  if ( weight  > 1.5 || weight  < 0.5 ) weight = 0;
	  if ( weight0 > 1.5 || weight0 < 0.5 ) weight0 = 0;
	  if ( weight1 > 1.5 || weight1 < 0.5 ) weight1 = 0;
	  if ( weight2 > 1.5 || weight2 < 0.5 ) weight2 = 0;
	  if ( weight3 > 1.5 || weight3 < 0.5 ) weight3 = 0;
	  // ---
	  th1d_weight_fvtxn_zvtxbin_clus_phi_GC->SetBinContent(j+1,weight);
	  th1d_weight_fvtxn0_zvtxbin_clus_phi_GC->SetBinContent(j+1,weight0);
	  th1d_weight_fvtxn1_zvtxbin_clus_phi_GC->SetBinContent(j+1,weight1);
	  th1d_weight_fvtxn2_zvtxbin_clus_phi_GC->SetBinContent(j+1,weight2);
	  th1d_weight_fvtxn3_zvtxbin_clus_phi_GC->SetBinContent(j+1,weight3);
	}

      th1d_weight_fvtxn_zvtxbin_clus_phi_GC->Draw();
      c1->Print(Form("FigsWeight/weight2d_fvtxn_zvtx%d_clus_phi_GC_run%d.png",i,run));
      th1d_fvtxn_zvtxbin_clus_phi_GC->Multiply(th1d_weight_fvtxn_zvtxbin_clus_phi_GC);
      th1d_fvtxn_zvtxbin_clus_phi_GC->Draw();
      c1->Print(Form("FigsWeight/check2d_fvtxn_zvtx%d_clus_phi_GC_run%d.png",i,run));

      th1d_weight_fvtxn0_zvtxbin_clus_phi_GC->Draw();
      c1->Print(Form("FigsWeight/weight2d_fvtxn0_zvtx%d_clus_phi_GC_run%d.png",i,run));
      th1d_fvtxn0_zvtxbin_clus_phi_GC->Multiply(th1d_weight_fvtxn0_zvtxbin_clus_phi_GC);
      th1d_fvtxn0_zvtxbin_clus_phi_GC->Draw();
      c1->Print(Form("FigsWeight/check2d_fvtxn0_zvtx%d_clus_phi_GC_run%d.png",i,run));

      th1d_weight_fvtxn1_zvtxbin_clus_phi_GC->Draw();
      c1->Print(Form("FigsWeight/weight2d_fvtxn1_zvtx%d_clus_phi_GC_run%d.png",i,run));
      th1d_fvtxn1_zvtxbin_clus_phi_GC->Multiply(th1d_weight_fvtxn1_zvtxbin_clus_phi_GC);
      th1d_fvtxn1_zvtxbin_clus_phi_GC->Draw();
      c1->Print(Form("FigsWeight/check2d_fvtxn1_zvtx%d_clus_phi_GC_run%d.png",i,run));

      th1d_weight_fvtxn2_zvtxbin_clus_phi_GC->Draw();
      c1->Print(Form("FigsWeight/weight2d_fvtxn2_zvtx%d_clus_phi_GC_run%d.png",i,run));
      th1d_fvtxn2_zvtxbin_clus_phi_GC->Multiply(th1d_weight_fvtxn2_zvtxbin_clus_phi_GC);
      th1d_fvtxn2_zvtxbin_clus_phi_GC->Draw();
      c1->Print(Form("FigsWeight/check2d_fvtxn2_zvtx%d_clus_phi_GC_run%d.png",i,run));

      th1d_weight_fvtxn3_zvtxbin_clus_phi_GC->Draw();
      c1->Print(Form("FigsWeight/weight2d_fvtxn3_zvtx%d_clus_phi_GC_run%d.png",i,run));
      th1d_fvtxn3_zvtxbin_clus_phi_GC->Multiply(th1d_weight_fvtxn3_zvtxbin_clus_phi_GC);
      th1d_fvtxn3_zvtxbin_clus_phi_GC->Draw();
      c1->Print(Form("FigsWeight/check2d_fvtxn3_zvtx%d_clus_phi_GC_run%d.png",i,run));

    }

  fout->Write();
  fout->Close();

  file->Close();

  delete c1;

}


