void runweight2d(int);
void newweight2d(int);


void newweight()
{

  //runweight2d(456652); // remember when this was needed?
  //newweight2d(455355); // now this is needed, and that could change again in a little while (must be some weird environment thing)
  newweight2d(456652); // now this is needed, and that could change again in a little while (must be some weird environment thing)
  //runweight2d(455355);

  //  return;

  int run;
  ifstream fin;

  // fin.open("list_20.short");
  // while ( fin >> run ) newweight2d(run);
  // fin.close();

  // fin.open("list_39.short");
  // while ( fin >> run ) newweight2d(run);
  // fin.close();

  // fin.open("list_62.short");
  // while ( fin >> run ) newweight2d(run);
  // fin.close();

  // fin.open("list_200.short");
  // while ( fin >> run ) newweight2d(run);
  // fin.close();

}


void newweight2d(int run)
{

  cout << "Now processing run " << run << endl;

  TCanvas* c1 = new TCanvas("c1","");

  TFile* file = TFile::Open(Form("RootFiles/weight_run%d_pass0.root",run));
  if ( !file )
    {
      cout << "WARNING: file does not exist for run " << run << endl;
      return;
    }

  // ------------- //
  // --- south --- //
  // ------------- //

  TFile* fout = TFile::Open(Form("WeightFiles/newweight2d_run%d.root",run),"recreate");

  TH2D* th2d_fvtxs_clus_phi = (TH2D*)((TH2D*)file->Get("th2d_fvtxs_clus_phi"))->Clone("th2d_fvtxs_clus_phi");
  TH2D* th2d_fvtxs0_clus_phi = (TH2D*)((TH2D*)file->Get("th2d_fvtxs0_clus_phi"))->Clone("th2d_fvtxs0_clus_phi");
  TH2D* th2d_fvtxs1_clus_phi = (TH2D*)((TH2D*)file->Get("th2d_fvtxs1_clus_phi"))->Clone("th2d_fvtxs1_clus_phi");
  TH2D* th2d_fvtxs2_clus_phi = (TH2D*)((TH2D*)file->Get("th2d_fvtxs2_clus_phi"))->Clone("th2d_fvtxs2_clus_phi");
  TH2D* th2d_fvtxs3_clus_phi = (TH2D*)((TH2D*)file->Get("th2d_fvtxs3_clus_phi"))->Clone("th2d_fvtxs3_clus_phi");

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
          double temp = average/bincontent;
          if ( temp != temp ) temp = 0;
          if ( bincontent < 2 ) temp = 0;
          if ( i == 10 )
            {
              cout << "bincontent " << bincontent << endl;
              cout << "average " << average << endl;
              cout << "weight " << temp <<  " " << i+1 << " " << j+1 << endl;
            }
          th2d_fvtxs_clus_phi_weight->SetBinContent(i+1,j+1,temp);
          if ( i == 10 )cout << th2d_fvtxs_clus_phi_weight->GetBinContent(i+1,j+1,temp) << endl;
          // ---
          bincontent = th2d_fvtxs0_clus_phi->GetBinContent(i+1,j+1);
          temp = average/bincontent;
          if ( temp != temp ) temp = 0;
          if ( bincontent < 2 ) temp = 0;
          if ( i == 10 )
            {
              cout << "bincontent " << bincontent << endl;
              cout << "average " << average << endl;
              cout << "weight " << temp <<  " " << i+1 << " " << j+1 << endl;
            }
          th2d_fvtxs0_clus_phi_weight->SetBinContent(i+1,j+1,temp);
          if ( i == 10 )cout << th2d_fvtxs0_clus_phi_weight->GetBinContent(i+1,j+1,temp) << endl;
          // ---
          bincontent = th2d_fvtxs1_clus_phi->GetBinContent(i+1,j+1);
          temp = average/bincontent;
          if ( temp != temp ) temp = 0;
          if ( bincontent < 2 ) temp = 0;
          if ( i == 10 )
            {
              cout << "bincontent " << bincontent << endl;
              cout << "average " << average << endl;
              cout << "weight " << temp <<  " " << i+1 << " " << j+1 << endl;
            }
          th2d_fvtxs1_clus_phi_weight->SetBinContent(i+1,j+1,temp);
          if ( i == 10 )cout << th2d_fvtxs1_clus_phi_weight->GetBinContent(i+1,j+1,temp) << endl;
          // ---
          bincontent = th2d_fvtxs2_clus_phi->GetBinContent(i+1,j+1);
          temp = average/bincontent;
          if ( temp != temp ) temp = 0;
          if ( bincontent < 2 ) temp = 0;
          if ( i == 10 )
            {
              cout << "bincontent " << bincontent << endl;
              cout << "average " << average << endl;
              cout << "weight " << temp <<  " " << i+1 << " " << j+1 << endl;
            }
          th2d_fvtxs2_clus_phi_weight->SetBinContent(i+1,j+1,temp);
          if ( i == 10 )cout << th2d_fvtxs2_clus_phi_weight->GetBinContent(i+1,j+1,temp) << endl;
          // ---
          bincontent = th2d_fvtxs3_clus_phi->GetBinContent(i+1,j+1);
          temp = average/bincontent;
          if ( temp != temp ) temp = 0;
          if ( bincontent < 2 ) temp = 0;
          if ( i == 10 )
            {
              cout << "bincontent " << bincontent << endl;
              cout << "average " << average << endl;
              cout << "weight " << temp <<  " " << i+1 << " " << j+1 << endl;
            }
          th2d_fvtxs3_clus_phi_weight->SetBinContent(i+1,j+1,temp);
          if ( i == 10 )cout << th2d_fvtxs3_clus_phi_weight->GetBinContent(i+1,j+1,temp) << endl;
        }
    }

  // ------------- //
  // --- north --- //
  // ------------- //

  TH2D* th2d_fvtxn_clus_phi = (TH2D*)file->Get("th2d_fvtxn_clus_phi");
  TH2D* th2d_fvtxn0_clus_phi = (TH2D*)file->Get("th2d_fvtxn0_clus_phi");
  TH2D* th2d_fvtxn1_clus_phi = (TH2D*)file->Get("th2d_fvtxn1_clus_phi");
  TH2D* th2d_fvtxn2_clus_phi = (TH2D*)file->Get("th2d_fvtxn2_clus_phi");
  TH2D* th2d_fvtxn3_clus_phi = (TH2D*)file->Get("th2d_fvtxn3_clus_phi");

  ihavethehistos = th2d_fvtxn_clus_phi && th2d_fvtxn0_clus_phi && th2d_fvtxn1_clus_phi && th2d_fvtxn2_clus_phi && th2d_fvtxn3_clus_phi;

  if ( !ihavethehistos )
    {
      cout << "YOU'RE GONNA DIE (missing histograms)" << endl;
      return;
    }

  if ( th2d_fvtxn_clus_phi->GetNbinsX() != nbinsx )
    {
      cout << "YOU'RE GONNA DIE " << nbinsx << " " << th2d_fvtxn_clus_phi->GetNbinsX() << endl;
      return;
    }
  if ( th2d_fvtxn_clus_phi->GetNbinsY() != nbinsy )
    {
      cout << "YOU'RE GONNA DIE " << nbinsy << " " << th2d_fvtxn_clus_phi->GetNbinsX() << endl;
      return;
    }

  integral = th2d_fvtxn_clus_phi->Integral(1,nbinsx,1,nbinsy);
  totalbins = (double)nbinsx*(double)nbinsy;
  average = integral/totalbins;

  cout << "average is " << average << endl;

  TH2D* th2d_fvtxn_clus_phi_weight = (TH2D*)th2d_fvtxn_clus_phi->Clone("th2d_fvtxn_clus_phi_weight");
  TH2D* th2d_fvtxn0_clus_phi_weight = (TH2D*)th2d_fvtxn0_clus_phi->Clone("th2d_fvtxn0_clus_phi_weight");
  TH2D* th2d_fvtxn1_clus_phi_weight = (TH2D*)th2d_fvtxn1_clus_phi->Clone("th2d_fvtxn1_clus_phi_weight");
  TH2D* th2d_fvtxn2_clus_phi_weight = (TH2D*)th2d_fvtxn2_clus_phi->Clone("th2d_fvtxn2_clus_phi_weight");
  TH2D* th2d_fvtxn3_clus_phi_weight = (TH2D*)th2d_fvtxn3_clus_phi->Clone("th2d_fvtxn3_clus_phi_weight");
  for ( int i = 0; i < nbinsx; ++i )
    {
      for ( int j = 0; j < nbinsy; ++j )
        {
          double bincontent = th2d_fvtxn_clus_phi->GetBinContent(i+1,j+1);
          double temp = average/bincontent;
          if ( temp != temp ) temp = 0;
          if ( bincontent < 2 ) temp = 0;
          th2d_fvtxn_clus_phi_weight->SetBinContent(i+1,j+1,temp);
          // ---
          bincontent = th2d_fvtxn0_clus_phi->GetBinContent(i+1,j+1);
          temp = average/bincontent;
          if ( temp != temp ) temp = 0;
          if ( bincontent < 2 ) temp = 0;
          th2d_fvtxn0_clus_phi_weight->SetBinContent(i+1,j+1,temp);
          // ---
          bincontent = th2d_fvtxn1_clus_phi->GetBinContent(i+1,j+1);
          temp = average/bincontent;
          if ( temp != temp ) temp = 0;
          if ( bincontent < 2 ) temp = 0;
          th2d_fvtxn1_clus_phi_weight->SetBinContent(i+1,j+1,temp);
          // ---
          bincontent = th2d_fvtxn2_clus_phi->GetBinContent(i+1,j+1);
          temp = average/bincontent;
          if ( temp != temp ) temp = 0;
          if ( bincontent < 2 ) temp = 0;
          th2d_fvtxn2_clus_phi_weight->SetBinContent(i+1,j+1,temp);
          // ---
          bincontent = th2d_fvtxn3_clus_phi->GetBinContent(i+1,j+1);
          temp = average/bincontent;
          if ( temp != temp ) temp = 0;
          if ( bincontent < 2 ) temp = 0;
          th2d_fvtxn3_clus_phi_weight->SetBinContent(i+1,j+1,temp);
        }
    }

  fout->Write();
  //  fout->Close();

  // ---

  th2d_fvtxs_clus_phi->Draw("colz");
  c1->Print(Form("FigsWeight/testing_newweight2ddist_fvtxs_%d.png",run));
  th2d_fvtxs_clus_phi_weight->Draw("colz");
  c1->Print(Form("FigsWeight/testing_newweight2dweight_fvtxs_%d.png",run));

  th2d_fvtxs0_clus_phi->Draw("colz");
  c1->Print(Form("FigsWeight/newweight2ddist_fvtxs0_%d.png",run));
  th2d_fvtxs0_clus_phi_weight->Draw("colz");
  c1->Print(Form("FigsWeight/newweight2dweight_fvtxs0_%d.png",run));

  th2d_fvtxs1_clus_phi->Draw("colz");
  c1->Print(Form("FigsWeight/newweight2ddist_fvtxs1_%d.png",run));
  th2d_fvtxs1_clus_phi_weight->Draw("colz");
  c1->Print(Form("FigsWeight/newweight2dweight_fvtxs1_%d.png",run));

  th2d_fvtxs2_clus_phi->Draw("colz");
  c1->Print(Form("FigsWeight/newweight2ddist_fvtxs2_%d.png",run));
  th2d_fvtxs2_clus_phi_weight->Draw("colz");
  c1->Print(Form("FigsWeight/newweight2dweight_fvtxs2_%d.png",run));

  th2d_fvtxs3_clus_phi->Draw("colz");
  c1->Print(Form("FigsWeight/newweight2ddist_fvtxs3_%d.png",run));
  th2d_fvtxs3_clus_phi_weight->Draw("colz");
  c1->Print(Form("FigsWeight/newweight2dweight_fvtxs3_%d.png",run));

  // ---

  th2d_fvtxn_clus_phi->Draw("colz");
  c1->Print(Form("FigsWeight/testing_newweight2ddist_fvtxn_%d.png",run));
  th2d_fvtxn_clus_phi_weight->Draw("colz");
  c1->Print(Form("FigsWeight/testing_newweight2dweight_fvtxn_%d.png",run));

  th2d_fvtxn0_clus_phi->Draw("colz");
  c1->Print(Form("FigsWeight/newweight2ddist_fvtxn0_%d.png",run));
  th2d_fvtxn0_clus_phi_weight->Draw("colz");
  c1->Print(Form("FigsWeight/newweight2dweight_fvtxn0_%d.png",run));

  th2d_fvtxn1_clus_phi->Draw("colz");
  c1->Print(Form("FigsWeight/newweight2ddist_fvtxn1_%d.png",run));
  th2d_fvtxn1_clus_phi_weight->Draw("colz");
  c1->Print(Form("FigsWeight/newweight2dweight_fvtxn1_%d.png",run));

  th2d_fvtxn2_clus_phi->Draw("colz");
  c1->Print(Form("FigsWeight/newweight2ddist_fvtxn2_%d.png",run));
  th2d_fvtxn2_clus_phi_weight->Draw("colz");
  c1->Print(Form("FigsWeight/newweight2dweight_fvtxn2_%d.png",run));

  th2d_fvtxn3_clus_phi->Draw("colz");
  c1->Print(Form("FigsWeight/newweight2ddist_fvtxn3_%d.png",run));
  th2d_fvtxn3_clus_phi_weight->Draw("colz");
  c1->Print(Form("FigsWeight/newweight2dweight_fvtxn3_%d.png",run));

  fout->Close();

  delete c1;

} // newweight2d



