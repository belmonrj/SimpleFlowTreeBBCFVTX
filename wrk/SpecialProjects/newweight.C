void newweight2d(int);

void bbctube(int);

void fvtx_track_weight1d(int);

void newweight()
{

  //fvtx_track_weight1d(455355);

  //return;

  //bbctube(456652);

  //  return;

  //runweight2d(456652); // remember when this was needed?
  //newweight2d(455355); // now this is needed, and that could change again in a little while (must be some weird environment thing)
  //newweight2d(456652); // now this is needed, and that could change again in a little while (must be some weird environment thing)
  //runweight2d(455355);

  // newweight2d(454774);
  // bbctube(454774);
  // fvtx_track_weight1d(454774);

  newweight2d(456652);
  bbctube(456652);
  fvtx_track_weight1d(456652);

  //return;

  int run;
  ifstream fin;

  // fin.open("list200.txt");
  // while ( fin >> run )
  //   {
  //     newweight2d(run);
  //     bbctube(run);
  //     fvtx_track_weight1d(run);
  //   }
  // fin.close();

  // fin.open("list62.txt");
  // while ( fin >> run )
  //   {
  //     newweight2d(run);
  //     bbctube(run);
  //     fvtx_track_weight1d(run);
  //   }
  // fin.close();

  // fin.open("list39.txt");
  // while ( fin >> run )
  //   {
  //     newweight2d(run);
  //     bbctube(run);
  //     fvtx_track_weight1d(run);
  //   }
  // fin.close();

  // fin.open("list20.txt");
  // while ( fin >> run )
  //   {
  //     newweight2d(run);
  //     bbctube(run);
  //     fvtx_track_weight1d(run);
  //   }
  // fin.close();

}


void newweight2d(int run)
{

  cout << "Now processing run " << run << " to do cluster weights" << endl;

  TCanvas* c1 = new TCanvas("c1","");

  TFile* file = TFile::Open(Form("RootFiles/hist_%d.root",run));
  if ( !file )
    {
      cout << "WARNING: file does not exist for run " << run << endl;
      return;
    }
  TList* list = (TList*)file->GetListOfKeys();
  int size = list->GetSize();
  if ( size < 1 )
    {
      cout << "WARNING: file does not have enough keys, run " << run << endl;
      return;
    }

  // ------------- //
  // --- south --- //
  // ------------- //

  TFile* fout = TFile::Open(Form("WeightFiles/newweight2d_run%d.root",run),"recreate");

  // TH2D* th2d_fvtxs_clus_phi = (TH2D*)((TH2D*)file->Get("th2d_fvtxs_clus_phi"))->Clone("th2d_fvtxs_clus_phi");
  // TH2D* th2d_fvtxs0_clus_phi = (TH2D*)((TH2D*)file->Get("th2d_fvtxs0_clus_phi"))->Clone("th2d_fvtxs0_clus_phi");
  // TH2D* th2d_fvtxs1_clus_phi = (TH2D*)((TH2D*)file->Get("th2d_fvtxs1_clus_phi"))->Clone("th2d_fvtxs1_clus_phi");
  // TH2D* th2d_fvtxs2_clus_phi = (TH2D*)((TH2D*)file->Get("th2d_fvtxs2_clus_phi"))->Clone("th2d_fvtxs2_clus_phi");
  // TH2D* th2d_fvtxs3_clus_phi = (TH2D*)((TH2D*)file->Get("th2d_fvtxs3_clus_phi"))->Clone("th2d_fvtxs3_clus_phi");

  // bool ihavethehistos = th2d_fvtxs_clus_phi && th2d_fvtxs0_clus_phi && th2d_fvtxs1_clus_phi && th2d_fvtxs2_clus_phi && th2d_fvtxs3_clus_phi;

  // if ( !ihavethehistos )
  //   {
  //     cout << "YOU'RE GONNA DIE (missing histograms)" << endl;
  //     return;
  //   }

  // const int nbinsx = 20;
  // if ( th2d_fvtxs_clus_phi->GetNbinsX() != nbinsx )
  //   {
  //     cout << "YOU'RE GONNA DIE " << nbinsx << " " << th2d_fvtxs_clus_phi->GetNbinsX() << endl;
  //     return;
  //   }
  // const int nbinsy = 50;
  // if ( th2d_fvtxs_clus_phi->GetNbinsY() != nbinsy )
  //   {
  //     cout << "YOU'RE GONNA DIE " << nbinsy << " " << th2d_fvtxs_clus_phi->GetNbinsX() << endl;
  //     return;
  //   }

  // double integral = th2d_fvtxs_clus_phi->Integral(1,nbinsx,1,nbinsy);
  // double totalbins = (double)nbinsx*(double)nbinsy;
  // double average = integral/totalbins;

  // cout << "average is " << average << endl;

  // TH2D* th2d_fvtxs_clus_phi_weight = (TH2D*)th2d_fvtxs_clus_phi->Clone("th2d_fvtxs_clus_phi_weight");
  // TH2D* th2d_fvtxs0_clus_phi_weight = (TH2D*)th2d_fvtxs0_clus_phi->Clone("th2d_fvtxs0_clus_phi_weight");
  // TH2D* th2d_fvtxs1_clus_phi_weight = (TH2D*)th2d_fvtxs1_clus_phi->Clone("th2d_fvtxs1_clus_phi_weight");
  // TH2D* th2d_fvtxs2_clus_phi_weight = (TH2D*)th2d_fvtxs2_clus_phi->Clone("th2d_fvtxs2_clus_phi_weight");
  // TH2D* th2d_fvtxs3_clus_phi_weight = (TH2D*)th2d_fvtxs3_clus_phi->Clone("th2d_fvtxs3_clus_phi_weight");
  // for ( int i = 0; i < nbinsx; ++i )
  //   {
  //     for ( int j = 0; j < nbinsy; ++j )
  //       {
  //         double bincontent = th2d_fvtxs_clus_phi->GetBinContent(i+1,j+1);
  //         double temp = average/bincontent;
  //         if ( temp != temp ) temp = 0;
  //         if ( bincontent < 2 ) temp = 0;
  //         if ( i == 99 )
  //           {
  //             cout << "bincontent " << bincontent << endl;
  //             cout << "average " << average << endl;
  //             cout << "weight " << temp <<  " " << i+1 << " " << j+1 << endl;
  //           }
  //         th2d_fvtxs_clus_phi_weight->SetBinContent(i+1,j+1,temp);
  //         if ( i == 99 )cout << th2d_fvtxs_clus_phi_weight->GetBinContent(i+1,j+1,temp) << endl;
  //         // ---
  //         bincontent = th2d_fvtxs0_clus_phi->GetBinContent(i+1,j+1);
  //         temp = average/bincontent;
  //         if ( temp != temp ) temp = 0;
  //         if ( bincontent < 2 ) temp = 0;
  //         if ( i == 99 )
  //           {
  //             cout << "bincontent " << bincontent << endl;
  //             cout << "average " << average << endl;
  //             cout << "weight " << temp <<  " " << i+1 << " " << j+1 << endl;
  //           }
  //         th2d_fvtxs0_clus_phi_weight->SetBinContent(i+1,j+1,temp);
  //         if ( i == 99 )cout << th2d_fvtxs0_clus_phi_weight->GetBinContent(i+1,j+1,temp) << endl;
  //         // ---
  //         bincontent = th2d_fvtxs1_clus_phi->GetBinContent(i+1,j+1);
  //         temp = average/bincontent;
  //         if ( temp != temp ) temp = 0;
  //         if ( bincontent < 2 ) temp = 0;
  //         if ( i == 99 )
  //           {
  //             cout << "bincontent " << bincontent << endl;
  //             cout << "average " << average << endl;
  //             cout << "weight " << temp <<  " " << i+1 << " " << j+1 << endl;
  //           }
  //         th2d_fvtxs1_clus_phi_weight->SetBinContent(i+1,j+1,temp);
  //         if ( i == 99 )cout << th2d_fvtxs1_clus_phi_weight->GetBinContent(i+1,j+1,temp) << endl;
  //         // ---
  //         bincontent = th2d_fvtxs2_clus_phi->GetBinContent(i+1,j+1);
  //         temp = average/bincontent;
  //         if ( temp != temp ) temp = 0;
  //         if ( bincontent < 2 ) temp = 0;
  //         if ( i == 99 )
  //           {
  //             cout << "bincontent " << bincontent << endl;
  //             cout << "average " << average << endl;
  //             cout << "weight " << temp <<  " " << i+1 << " " << j+1 << endl;
  //           }
  //         th2d_fvtxs2_clus_phi_weight->SetBinContent(i+1,j+1,temp);
  //         if ( i == 99 )cout << th2d_fvtxs2_clus_phi_weight->GetBinContent(i+1,j+1,temp) << endl;
  //         // ---
  //         bincontent = th2d_fvtxs3_clus_phi->GetBinContent(i+1,j+1);
  //         temp = average/bincontent;
  //         if ( temp != temp ) temp = 0;
  //         if ( bincontent < 2 ) temp = 0;
  //         if ( i == 99 )
  //           {
  //             cout << "bincontent " << bincontent << endl;
  //             cout << "average " << average << endl;
  //             cout << "weight " << temp <<  " " << i+1 << " " << j+1 << endl;
  //           }
  //         th2d_fvtxs3_clus_phi_weight->SetBinContent(i+1,j+1,temp);
  //         if ( i == 99 )cout << th2d_fvtxs3_clus_phi_weight->GetBinContent(i+1,j+1,temp) << endl;
  //       }
  //   }

  // // ------------- //
  // // --- north --- //
  // // ------------- //

  // TH2D* th2d_fvtxn_clus_phi = (TH2D*)((TH2D*)file->Get("th2d_fvtxn_clus_phi"))->Clone("th2d_fvtxn_clus_phi");
  // TH2D* th2d_fvtxn0_clus_phi = (TH2D*)((TH2D*)file->Get("th2d_fvtxn0_clus_phi"))->Clone("th2d_fvtxn0_clus_phi");
  // TH2D* th2d_fvtxn1_clus_phi = (TH2D*)((TH2D*)file->Get("th2d_fvtxn1_clus_phi"))->Clone("th2d_fvtxn1_clus_phi");
  // TH2D* th2d_fvtxn2_clus_phi = (TH2D*)((TH2D*)file->Get("th2d_fvtxn2_clus_phi"))->Clone("th2d_fvtxn2_clus_phi");
  // TH2D* th2d_fvtxn3_clus_phi = (TH2D*)((TH2D*)file->Get("th2d_fvtxn3_clus_phi"))->Clone("th2d_fvtxn3_clus_phi");

  // ihavethehistos = th2d_fvtxn_clus_phi && th2d_fvtxn0_clus_phi && th2d_fvtxn1_clus_phi && th2d_fvtxn2_clus_phi && th2d_fvtxn3_clus_phi;

  // if ( !ihavethehistos )
  //   {
  //     cout << "YOU'RE GONNA DIE (missing histograms)" << endl;
  //     return;
  //   }

  // if ( th2d_fvtxn_clus_phi->GetNbinsX() != nbinsx )
  //   {
  //     cout << "YOU'RE GONNA DIE " << nbinsx << " " << th2d_fvtxn_clus_phi->GetNbinsX() << endl;
  //     return;
  //   }
  // if ( th2d_fvtxn_clus_phi->GetNbinsY() != nbinsy )
  //   {
  //     cout << "YOU'RE GONNA DIE " << nbinsy << " " << th2d_fvtxn_clus_phi->GetNbinsX() << endl;
  //     return;
  //   }

  // integral = th2d_fvtxn_clus_phi->Integral(1,nbinsx,1,nbinsy);
  // totalbins = (double)nbinsx*(double)nbinsy;
  // average = integral/totalbins;

  // cout << "average is " << average << endl;

  // TH2D* th2d_fvtxn_clus_phi_weight = (TH2D*)th2d_fvtxn_clus_phi  ->Clone("th2d_fvtxn_clus_phi_weight");
  // TH2D* th2d_fvtxn0_clus_phi_weight = (TH2D*)th2d_fvtxn0_clus_phi->Clone("th2d_fvtxn0_clus_phi_weight");
  // TH2D* th2d_fvtxn1_clus_phi_weight = (TH2D*)th2d_fvtxn1_clus_phi->Clone("th2d_fvtxn1_clus_phi_weight");
  // TH2D* th2d_fvtxn2_clus_phi_weight = (TH2D*)th2d_fvtxn2_clus_phi->Clone("th2d_fvtxn2_clus_phi_weight");
  // TH2D* th2d_fvtxn3_clus_phi_weight = (TH2D*)th2d_fvtxn3_clus_phi->Clone("th2d_fvtxn3_clus_phi_weight");
  // for ( int i = 0; i < nbinsx; ++i )
  //   {
  //     for ( int j = 0; j < nbinsy; ++j )
  //       {
  //         double bincontent = th2d_fvtxn_clus_phi->GetBinContent(i+1,j+1);
  //         double temp = average/bincontent;
  //         if ( temp != temp ) temp = 0;
  //         if ( bincontent < 2 ) temp = 0;
  //         th2d_fvtxn_clus_phi_weight->SetBinContent(i+1,j+1,temp);
  //         // ---
  //         bincontent = th2d_fvtxn0_clus_phi->GetBinContent(i+1,j+1);
  //         temp = average/bincontent;
  //         if ( temp != temp ) temp = 0;
  //         if ( bincontent < 2 ) temp = 0;
  //         th2d_fvtxn0_clus_phi_weight->SetBinContent(i+1,j+1,temp);
  //         // ---
  //         bincontent = th2d_fvtxn1_clus_phi->GetBinContent(i+1,j+1);
  //         temp = average/bincontent;
  //         if ( temp != temp ) temp = 0;
  //         if ( bincontent < 2 ) temp = 0;
  //         th2d_fvtxn1_clus_phi_weight->SetBinContent(i+1,j+1,temp);
  //         // ---
  //         bincontent = th2d_fvtxn2_clus_phi->GetBinContent(i+1,j+1);
  //         temp = average/bincontent;
  //         if ( temp != temp ) temp = 0;
  //         if ( bincontent < 2 ) temp = 0;
  //         th2d_fvtxn2_clus_phi_weight->SetBinContent(i+1,j+1,temp);
  //         // ---
  //         bincontent = th2d_fvtxn3_clus_phi->GetBinContent(i+1,j+1);
  //         temp = average/bincontent;
  //         if ( temp != temp ) temp = 0;
  //         if ( bincontent < 2 ) temp = 0;
  //         th2d_fvtxn3_clus_phi_weight->SetBinContent(i+1,j+1,temp);
  //       }
  //   }

  // ---
  // ---
  // ---

  TH1D* th1d_fvtxs_clus_phi  = (TH1D*)((TH1D*)file->Get("th1d_fvtxs_clus_phi" ))->Clone();
  TH1D* th1d_fvtxs0_clus_phi = (TH1D*)((TH1D*)file->Get("th1d_fvtxs0_clus_phi"))->Clone();
  TH1D* th1d_fvtxs1_clus_phi = (TH1D*)((TH1D*)file->Get("th1d_fvtxs1_clus_phi"))->Clone();
  TH1D* th1d_fvtxs2_clus_phi = (TH1D*)((TH1D*)file->Get("th1d_fvtxs2_clus_phi"))->Clone();
  TH1D* th1d_fvtxs3_clus_phi = (TH1D*)((TH1D*)file->Get("th1d_fvtxs3_clus_phi"))->Clone();

  TH1D* th1d_weight_fvtxs  = (TH1D*)th1d_fvtxs_clus_phi ->Clone("th1d_fvtxs_clus_phi_weight");
  TH1D* th1d_weight_fvtxs0 = (TH1D*)th1d_fvtxs0_clus_phi->Clone("th1d_fvtxs0_clus_phi_weight");
  TH1D* th1d_weight_fvtxs1 = (TH1D*)th1d_fvtxs1_clus_phi->Clone("th1d_fvtxs1_clus_phi_weight");
  TH1D* th1d_weight_fvtxs2 = (TH1D*)th1d_fvtxs2_clus_phi->Clone("th1d_fvtxs2_clus_phi_weight");
  TH1D* th1d_weight_fvtxs3 = (TH1D*)th1d_fvtxs3_clus_phi->Clone("th1d_fvtxs3_clus_phi_weight");

  const int nbins = 50;
  if ( th1d_fvtxs_clus_phi->GetNbinsX() != nbins )
    {
      cout << "YOU'RE GONNA DIE" << endl;
      return;
    }

  double ave = th1d_fvtxs_clus_phi->Integral(1,nbins); // use 1 and nbins to exclude underflow (0) and overflow (nbins+1)
  double ave0 = th1d_fvtxs0_clus_phi->Integral(1,nbins);
  double ave1 = th1d_fvtxs1_clus_phi->Integral(1,nbins);
  double ave2 = th1d_fvtxs2_clus_phi->Integral(1,nbins);
  double ave3 = th1d_fvtxs3_clus_phi->Integral(1,nbins);
  ave /= nbins;
  ave0 /= nbins;
  ave1 /= nbins;
  ave2 /= nbins;
  ave3 /= nbins;
  for ( int i = 0; i < nbins; ++i )
    {
      double phi = th1d_fvtxs_clus_phi->GetBinCenter(i+1);
      double weight = ave/th1d_fvtxs_clus_phi->GetBinContent(i+1);
      double weight0 = ave0/th1d_fvtxs0_clus_phi->GetBinContent(i+1);
      double weight1 = ave1/th1d_fvtxs1_clus_phi->GetBinContent(i+1);
      double weight2 = ave2/th1d_fvtxs2_clus_phi->GetBinContent(i+1);
      double weight3 = ave3/th1d_fvtxs3_clus_phi->GetBinContent(i+1);
      if ( !TMath::Finite(weight) ) weight = 0;
      if ( !TMath::Finite(weight0) ) weight0 = 0;
      if ( !TMath::Finite(weight1) ) weight1 = 0;
      if ( !TMath::Finite(weight2) ) weight2 = 0;
      if ( !TMath::Finite(weight3) ) weight3 = 0;
      th1d_weight_fvtxs->SetBinContent(i+1,weight);
      th1d_weight_fvtxs0->SetBinContent(i+1,weight0);
      th1d_weight_fvtxs1->SetBinContent(i+1,weight1);
      th1d_weight_fvtxs2->SetBinContent(i+1,weight2);
      th1d_weight_fvtxs3->SetBinContent(i+1,weight3);
    }


  // ---
  // ---
  // ---

  TH1D* th1d_fvtxn_clus_phi  = (TH1D*)((TH1D*)file->Get("th1d_fvtxn_clus_phi" ))->Clone();
  TH1D* th1d_fvtxn0_clus_phi = (TH1D*)((TH1D*)file->Get("th1d_fvtxn0_clus_phi"))->Clone();
  TH1D* th1d_fvtxn1_clus_phi = (TH1D*)((TH1D*)file->Get("th1d_fvtxn1_clus_phi"))->Clone();
  TH1D* th1d_fvtxn2_clus_phi = (TH1D*)((TH1D*)file->Get("th1d_fvtxn2_clus_phi"))->Clone();
  TH1D* th1d_fvtxn3_clus_phi = (TH1D*)((TH1D*)file->Get("th1d_fvtxn3_clus_phi"))->Clone();

  TH1D* th1d_weight_fvtxn  = (TH1D*)th1d_fvtxn_clus_phi ->Clone("th1d_fvtxn_clus_phi_weight");
  TH1D* th1d_weight_fvtxn0 = (TH1D*)th1d_fvtxn0_clus_phi->Clone("th1d_fvtxn0_clus_phi_weight");
  TH1D* th1d_weight_fvtxn1 = (TH1D*)th1d_fvtxn1_clus_phi->Clone("th1d_fvtxn1_clus_phi_weight");
  TH1D* th1d_weight_fvtxn2 = (TH1D*)th1d_fvtxn2_clus_phi->Clone("th1d_fvtxn2_clus_phi_weight");
  TH1D* th1d_weight_fvtxn3 = (TH1D*)th1d_fvtxn3_clus_phi->Clone("th1d_fvtxn3_clus_phi_weight");

  if ( th1d_fvtxn_clus_phi->GetNbinsX() != nbins )
    {
      cout << "YOU'RE GONNA DIE" << endl;
      return;
    }

  double ave = th1d_fvtxn_clus_phi->Integral(1,nbins); // use 1 and nbins to exclude underflow (0) and overflow (nbins+1)
  double ave0 = th1d_fvtxn0_clus_phi->Integral(1,nbins);
  double ave1 = th1d_fvtxn1_clus_phi->Integral(1,nbins);
  double ave2 = th1d_fvtxn2_clus_phi->Integral(1,nbins);
  double ave3 = th1d_fvtxn3_clus_phi->Integral(1,nbins);
  ave /= nbins;
  ave0 /= nbins;
  ave1 /= nbins;
  ave2 /= nbins;
  ave3 /= nbins;
  for ( int i = 0; i < nbins; ++i )
    {
      double phi = th1d_fvtxn_clus_phi->GetBinCenter(i+1);
      double weight = ave/th1d_fvtxn_clus_phi->GetBinContent(i+1);
      double weight0 = ave0/th1d_fvtxn0_clus_phi->GetBinContent(i+1);
      double weight1 = ave1/th1d_fvtxn1_clus_phi->GetBinContent(i+1);
      double weight2 = ave2/th1d_fvtxn2_clus_phi->GetBinContent(i+1);
      double weight3 = ave3/th1d_fvtxn3_clus_phi->GetBinContent(i+1);
      if ( !TMath::Finite(weight) ) weight = 0;
      if ( !TMath::Finite(weight0) ) weight0 = 0;
      if ( !TMath::Finite(weight1) ) weight1 = 0;
      if ( !TMath::Finite(weight2) ) weight2 = 0;
      if ( !TMath::Finite(weight3) ) weight3 = 0;
      th1d_weight_fvtxn->SetBinContent(i+1,weight);
      th1d_weight_fvtxn0->SetBinContent(i+1,weight0);
      th1d_weight_fvtxn1->SetBinContent(i+1,weight1);
      th1d_weight_fvtxn2->SetBinContent(i+1,weight2);
      th1d_weight_fvtxn3->SetBinContent(i+1,weight3);
    }





  fout->Write();
  //  fout->Close();

  // ---

  // th2d_fvtxs_clus_phi->Draw("colz");
  // c1->Print(Form("FigsWeight/testing_newweight2ddist_fvtxs_%d.png",run));
  // th2d_fvtxs_clus_phi_weight->Draw("colz");
  // c1->Print(Form("FigsWeight/testing_newweight2dweight_fvtxs_%d.png",run));

  // th2d_fvtxs0_clus_phi->Draw("colz");
  // c1->Print(Form("FigsWeight/newweight2ddist_fvtxs0_%d.png",run));
  // th2d_fvtxs0_clus_phi_weight->Draw("colz");
  // c1->Print(Form("FigsWeight/newweight2dweight_fvtxs0_%d.png",run));

  // th2d_fvtxs1_clus_phi->Draw("colz");
  // c1->Print(Form("FigsWeight/newweight2ddist_fvtxs1_%d.png",run));
  // th2d_fvtxs1_clus_phi_weight->Draw("colz");
  // c1->Print(Form("FigsWeight/newweight2dweight_fvtxs1_%d.png",run));

  // th2d_fvtxs2_clus_phi->Draw("colz");
  // c1->Print(Form("FigsWeight/newweight2ddist_fvtxs2_%d.png",run));
  // th2d_fvtxs2_clus_phi_weight->Draw("colz");
  // c1->Print(Form("FigsWeight/newweight2dweight_fvtxs2_%d.png",run));

  // th2d_fvtxs3_clus_phi->Draw("colz");
  // c1->Print(Form("FigsWeight/newweight2ddist_fvtxs3_%d.png",run));
  // th2d_fvtxs3_clus_phi_weight->Draw("colz");
  // c1->Print(Form("FigsWeight/newweight2dweight_fvtxs3_%d.png",run));

  // // ---

  // th2d_fvtxn_clus_phi->Draw("colz");
  // c1->Print(Form("FigsWeight/testing_newweight2ddist_fvtxn_%d.png",run));
  // th2d_fvtxn_clus_phi_weight->Draw("colz");
  // c1->Print(Form("FigsWeight/testing_newweight2dweight_fvtxn_%d.png",run));

  // th2d_fvtxn0_clus_phi->Draw("colz");
  // c1->Print(Form("FigsWeight/newweight2ddist_fvtxn0_%d.png",run));
  // th2d_fvtxn0_clus_phi_weight->Draw("colz");
  // c1->Print(Form("FigsWeight/newweight2dweight_fvtxn0_%d.png",run));

  // th2d_fvtxn1_clus_phi->Draw("colz");
  // c1->Print(Form("FigsWeight/newweight2ddist_fvtxn1_%d.png",run));
  // th2d_fvtxn1_clus_phi_weight->Draw("colz");
  // c1->Print(Form("FigsWeight/newweight2dweight_fvtxn1_%d.png",run));

  // th2d_fvtxn2_clus_phi->Draw("colz");
  // c1->Print(Form("FigsWeight/newweight2ddist_fvtxn2_%d.png",run));
  // th2d_fvtxn2_clus_phi_weight->Draw("colz");
  // c1->Print(Form("FigsWeight/newweight2dweight_fvtxn2_%d.png",run));

  // th2d_fvtxn3_clus_phi->Draw("colz");
  // c1->Print(Form("FigsWeight/newweight2ddist_fvtxn3_%d.png",run));
  // th2d_fvtxn3_clus_phi_weight->Draw("colz");
  // c1->Print(Form("FigsWeight/newweight2dweight_fvtxn3_%d.png",run));

  fout->Close();

  delete c1;

} // newweight2d



void bbctube(int run)
{

  cout << "Now processing run " << run << " to do tube gain corrections " << endl;

  TCanvas* c1 = new TCanvas("c1","");

  TFile* file = TFile::Open(Form("RootFiles/hist_%d.root",run));
  if ( !file )
    {
      cout << "WARNING: file does not exist for run " << run << endl;
      return;
    }
  TList* list = (TList*)file->GetListOfKeys();
  int size = list->GetSize();
  if ( size < 1 )
    {
      cout << "WARNING: file does not have enough keys, run " << run << endl;
      return;
    }

  gROOT->ProcessLine(".L bbcrings.C");

  TProfile* tp1f_bbc_charge_tube = (TProfile*)file->Get("tp1f_bbc_charge_tube");
  TProfile* tp1f_bbc0_charge_tube = (TProfile*)file->Get("tp1f_bbc0_charge_tube");
  TProfile* tp1f_bbc1_charge_tube = (TProfile*)file->Get("tp1f_bbc1_charge_tube");
  TProfile* tp1f_bbc2_charge_tube = (TProfile*)file->Get("tp1f_bbc2_charge_tube");
  TProfile* tp1f_bbc3_charge_tube = (TProfile*)file->Get("tp1f_bbc3_charge_tube");
  TProfile* tp1f_bbc4_charge_tube = (TProfile*)file->Get("tp1f_bbc4_charge_tube");

  double mean0 = tp1f_bbc0_charge_tube->GetMean(2);
  double mean1 = tp1f_bbc1_charge_tube->GetMean(2);
  double mean2 = tp1f_bbc2_charge_tube->GetMean(2);
  double mean3 = tp1f_bbc3_charge_tube->GetMean(2);
  double mean4 = tp1f_bbc4_charge_tube->GetMean(2);

  cout << mean0 << endl;
  cout << mean1 << endl;
  cout << mean2 << endl;
  cout << mean3 << endl;
  cout << mean4 << endl;

  tp1f_bbc0_charge_tube->Draw();
  TLine line0(0,mean0,63,mean0);
  line0.Draw();
  c1->Print(Form("FigsWeight/bbctube0_charge_run%d.png",run));
  tp1f_bbc1_charge_tube->Draw();
  TLine line1(0,mean1,63,mean1);
  line1.Draw();
  c1->Print(Form("FigsWeight/bbctube1_charge_run%d.png",run));
  tp1f_bbc2_charge_tube->Draw();
  TLine line2(0,mean2,63,mean2);
  line2.Draw();
  c1->Print(Form("FigsWeight/bbctube2_charge_run%d.png",run));
  tp1f_bbc3_charge_tube->Draw();
  TLine line3(0,mean3,63,mean3);
  line3.Draw();
  c1->Print(Form("FigsWeight/bbctube3_charge_run%d.png",run));
  tp1f_bbc4_charge_tube->Draw();
  TLine line4(0,mean4,63,mean4);
  line4.Draw();
  c1->Print(Form("FigsWeight/bbctube4_charge_run%d.png",run));

  TFile* fout = TFile::Open(Form("WeightFiles/bbctube_run%d.root",run),"recreate");
  TH1D* th1d_tubegaincorrection = new TH1D("th1d_tubegaincorrection","",64,-0.5,63.5);
  for ( int i = 0; i < 64; ++i )
    {
      double rawcharge = 0;
      double gaincorrection = 0;
      int layer = bbcrings(i);
      // cout << i << " " << layer << endl;
      if ( layer == 0 ) { rawcharge = tp1f_bbc0_charge_tube->GetBinContent(i+1); gaincorrection = mean0 / rawcharge; }
      if ( layer == 1 ) { rawcharge = tp1f_bbc1_charge_tube->GetBinContent(i+1); gaincorrection = mean1 / rawcharge; }
      if ( layer == 2 ) { rawcharge = tp1f_bbc2_charge_tube->GetBinContent(i+1); gaincorrection = mean2 / rawcharge; }
      if ( layer == 3 ) { rawcharge = tp1f_bbc3_charge_tube->GetBinContent(i+1); gaincorrection = mean3 / rawcharge; }
      if ( layer == 4 ) { rawcharge = tp1f_bbc4_charge_tube->GetBinContent(i+1); gaincorrection = mean4 / rawcharge; }
      th1d_tubegaincorrection->SetBinContent(i+1,gaincorrection);
    }

  th1d_tubegaincorrection->Draw();
  c1->Print(Form("FigsWeight/bbctube_gaincorrection_run%d.png",run));

  fout->Write();
  fout->Close();

  delete c1;

}


// --- come back here

void fvtx_track_weight1d(int run)
{

  cout << "Now processing run " << run << " to do track weights " << endl;
  cout << "Haha no just kidding I'm not doing anything until we put those histograms into the flattening code" << endl;
  return;

  TCanvas* c1 = new TCanvas("c1","");

  TFile* file = TFile::Open(Form("RootFiles/hist_%d.root",run));
  if ( !file )
    {
      cout << "WARNING: file does not exist for run " << run << endl;
      return;
    }
  TList* list = (TList*)file->GetListOfKeys();
  int size = list->GetSize();
  if ( size < 1 )
    {
      cout << "WARNING: file does not have enough keys, run " << run << endl;
      return;
    }

  cout << "File name is " << file->GetName() << " and address is " << file << endl;

  // ---
  // ---
  // ---

  TFile* fout = TFile::Open(Form("WeightFiles/weight1d_run%d.root",run),"recreate");

  TH1D* th1d_fvtxs_track_phi  = ((TH1D*)file->Get("th1d_fvtxs_track_phi" ));
  TH1D* th1d_fvtxs0_track_phi = ((TH1D*)file->Get("th1d_fvtxs0_track_phi"));
  TH1D* th1d_fvtxs1_track_phi = ((TH1D*)file->Get("th1d_fvtxs1_track_phi"));

  if ( !th1d_fvtxs_track_phi ||
       !th1d_fvtxs0_track_phi ||
       !th1d_fvtxs1_track_phi )
    {
      cout << "no track weight histograms " << endl;
      return;
    }

  TH1D* th1d_weight_fvtxs  = (TH1D*)th1d_fvtxs_track_phi ->Clone("th1d_fvtxs_track_phi_weight");
  TH1D* th1d_weight_fvtxs0 = (TH1D*)th1d_fvtxs0_track_phi->Clone("th1d_fvtxs0_track_phi_weight");
  TH1D* th1d_weight_fvtxs1 = (TH1D*)th1d_fvtxs1_track_phi->Clone("th1d_fvtxs1_track_phi_weight");


  const int nbins = 50;
  if ( th1d_fvtxs_track_phi->GetNbinsX() != nbins )
    {
      cout << "YOU'RE GONNA DIE" << endl;
      return;
    }

  double ave = th1d_fvtxs_track_phi->Integral(1,nbins); // use 1 and nbins to exclude underflow (0) and overflow (nbins+1)
  double ave0 = th1d_fvtxs0_track_phi->Integral(1,nbins);
  double ave1 = th1d_fvtxs1_track_phi->Integral(1,nbins);
  // double ave2 = th1d_fvtxs2_track_phi->Integral(1,nbins);
  // double ave3 = th1d_fvtxs3_track_phi->Integral(1,nbins);
  ave /= nbins;
  ave0 /= nbins;
  ave1 /= nbins;
  // ave2 /= nbins;
  // ave3 /= nbins;
  for ( int i = 0; i < nbins; ++i )
    {
      double phi = th1d_fvtxs_track_phi->GetBinCenter(i+1);
      double weight = ave/th1d_fvtxs_track_phi->GetBinContent(i+1);
      double weight0 = ave0/th1d_fvtxs0_track_phi->GetBinContent(i+1);
      double weight1 = ave1/th1d_fvtxs1_track_phi->GetBinContent(i+1);
      // double weight2 = ave2/th1d_fvtxs2_track_phi->GetBinContent(i+1);
      // double weight3 = ave3/th1d_fvtxs3_track_phi->GetBinContent(i+1);
      if ( !TMath::Finite(weight) ) weight = 0;
      if ( !TMath::Finite(weight0) ) weight0 = 0;
      if ( !TMath::Finite(weight1) ) weight1 = 0;
      // if ( !TMath::Finite(weight2) ) weight2 = 0;
      // if ( !TMath::Finite(weight3) ) weight3 = 0;
      th1d_weight_fvtxs->SetBinContent(i+1,weight);
      th1d_weight_fvtxs0->SetBinContent(i+1,weight0);
      th1d_weight_fvtxs1->SetBinContent(i+1,weight1);
      // th1d_weight_fvtxs2->SetBinContent(i+1,weight2);
      // th1d_weight_fvtxs3->SetBinContent(i+1,weight3);
    }

  // ---
  // ---
  // ---

  TH1D* th1d_fvtxn_track_phi  = (TH1D*)((TH1D*)file->Get("th1d_fvtxn_track_phi" ))->Clone();
  TH1D* th1d_fvtxn0_track_phi = (TH1D*)((TH1D*)file->Get("th1d_fvtxn0_track_phi"))->Clone();
  TH1D* th1d_fvtxn1_track_phi = (TH1D*)((TH1D*)file->Get("th1d_fvtxn1_track_phi"))->Clone();
  // TH1D* th1d_fvtxn2_track_phi = (TH1D*)((TH1D*)file->Get("th1d_fvtxn2_track_phi"))->Clone();
  // TH1D* th1d_fvtxn3_track_phi = (TH1D*)((TH1D*)file->Get("th1d_fvtxn3_track_phi"))->Clone();

  TH1D* th1d_weight_fvtxn  = (TH1D*)th1d_fvtxn_track_phi ->Clone("th1d_fvtxn_track_phi_weight");
  TH1D* th1d_weight_fvtxn0 = (TH1D*)th1d_fvtxn0_track_phi->Clone("th1d_fvtxn0_track_phi_weight");
  TH1D* th1d_weight_fvtxn1 = (TH1D*)th1d_fvtxn1_track_phi->Clone("th1d_fvtxn1_track_phi_weight");
  // TH1D* th1d_weight_fvtxn2 = (TH1D*)th1d_fvtxn2_track_phi->Clone("th1d_fvtxn2_track_phi_weight");
  // TH1D* th1d_weight_fvtxn3 = (TH1D*)th1d_fvtxn3_track_phi->Clone("th1d_fvtxn3_track_phi_weight");

  if ( th1d_fvtxn_track_phi->GetNbinsX() != nbins )
    {
      cout << "YOU'RE GONNA DIE" << endl;
      return;
    }

  double ave = th1d_fvtxn_track_phi->Integral(1,nbins); // use 1 and nbins to exclude underflow (0) and overflow (nbins+1)
  double ave0 = th1d_fvtxn0_track_phi->Integral(1,nbins);
  double ave1 = th1d_fvtxn1_track_phi->Integral(1,nbins);
  // double ave2 = th1d_fvtxn2_track_phi->Integral(1,nbins);
  // double ave3 = th1d_fvtxn3_track_phi->Integral(1,nbins);
  ave /= nbins;
  ave0 /= nbins;
  ave1 /= nbins;
  // ave2 /= nbins;
  // ave3 /= nbins;
  for ( int i = 0; i < nbins; ++i )
    {
      double phi = th1d_fvtxn_track_phi->GetBinCenter(i+1);
      double weight = ave/th1d_fvtxn_track_phi->GetBinContent(i+1);
      double weight0 = ave0/th1d_fvtxn0_track_phi->GetBinContent(i+1);
      double weight1 = ave1/th1d_fvtxn1_track_phi->GetBinContent(i+1);
      // double weight2 = ave2/th1d_fvtxn2_track_phi->GetBinContent(i+1);
      // double weight3 = ave3/th1d_fvtxn3_track_phi->GetBinContent(i+1);
      if ( !TMath::Finite(weight) ) weight = 0;
      if ( !TMath::Finite(weight0) ) weight0 = 0;
      if ( !TMath::Finite(weight1) ) weight1 = 0;
      // if ( !TMath::Finite(weight2) ) weight2 = 0;
      // if ( !TMath::Finite(weight3) ) weight3 = 0;
      th1d_weight_fvtxn->SetBinContent(i+1,weight);
      th1d_weight_fvtxn0->SetBinContent(i+1,weight0);
      th1d_weight_fvtxn1->SetBinContent(i+1,weight1);
      // th1d_weight_fvtxn2->SetBinContent(i+1,weight2);
      // th1d_weight_fvtxn3->SetBinContent(i+1,weight3);
    }

  fout->Write();
  fout->Close();

  cout << "Hopefully wrote out file " << fout->GetName() << endl;

  delete c1;

} // newweight2d
