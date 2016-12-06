void doenergy(int);

void temp_decorrelation()
{

  doenergy(200);
  // doenergy(62);
  // doenergy(20);
  // doenergy(39);

}

void doenergy(int energy)
{

  gStyle->SetOptTitle(1);

  TCanvas* c1 = new TCanvas("c1","");

  TFile* file = TFile::Open(Form("input/cumulants_%d.root",energy));

  TProfile* tp1f_special_fvtx_tracks_aa_cos = (TProfile*)file->Get("tp1f_special_fvtx_tracks_aa_cos"); // acceptance correction terms
  TProfile* tp1f_special_fvtx_tracks_aa_sin = (TProfile*)file->Get("tp1f_special_fvtx_tracks_aa_sin"); // acceptance correction terms
  TProfile* tp1f_special_fvtx_tracks_aa = (TProfile*)file->Get("tp1f_special_fvtx_tracks_aa");
  TProfile* tp1f_special_fvtx_tracks_ab0 = (TProfile*)file->Get("tp1f_special_fvtx_tracks_ab0");
  TProfile* tp1f_special_fvtx_tracks_ab1 = (TProfile*)file->Get("tp1f_special_fvtx_tracks_ab1");
  TProfile* tp1f_special_fvtx_tracks_ab2 = (TProfile*)file->Get("tp1f_special_fvtx_tracks_ab2");
  TProfile* tp1f_special_fvtx_tracks_ab3 = (TProfile*)file->Get("tp1f_special_fvtx_tracks_ab3");
  TProfile* tp1f_special_fvtx_tracks_ab4 = (TProfile*)file->Get("tp1f_special_fvtx_tracks_ab4");
  TProfile* tp1f_special_fvtx_tracks_ab5 = (TProfile*)file->Get("tp1f_special_fvtx_tracks_ab5");
  TProfile* tp1f_special_fvtx_tracks_ab6 = (TProfile*)file->Get("tp1f_special_fvtx_tracks_ab6");
  TProfile* tp1f_special_fvtx_tracks_ab7 = (TProfile*)file->Get("tp1f_special_fvtx_tracks_ab7");

  // -------------------------------------------------------------------------------------------
  // --- bin 0
  // ---
  double content_ab0[12];
  double content_aa0[12];
  double content_bb0[12];
  double content_final0[12];
  double econtent_ab0[12];
  double econtent_aa0[12];
  double econtent_bb0[12];
  double econtent_final0[12];
  double x[12];
  // -------------------------------------------------------------------------------------------
  for ( int i = 0; i < 12; ++i )
    {
      // --- set the indices for particles a and b
      int index_a = 1;
      int index_b = i+1;
      // --- get the acceptance correction terms
      double cosb = tp1f_special_fvtx_tracks_aa_cos->GetBinContent(index_b);
      double cosa = tp1f_special_fvtx_tracks_aa_cos->GetBinContent(index_a);
      double sinb = tp1f_special_fvtx_tracks_aa_sin->GetBinContent(index_b);
      double sina = tp1f_special_fvtx_tracks_aa_sin->GetBinContent(index_a);
      // --- get eta_b
      x[i] = tp1f_special_fvtx_tracks_ab0->GetBinCenter(i+1);
      // --- get the components
      content_ab0[i] = tp1f_special_fvtx_tracks_ab0->GetBinContent(index_b);
      content_bb0[i] = tp1f_special_fvtx_tracks_aa->GetBinContent(index_b);
      content_aa0[i] = tp1f_special_fvtx_tracks_aa->GetBinContent(index_a);
      econtent_ab0[i] = tp1f_special_fvtx_tracks_ab0->GetBinError(index_b);
      econtent_bb0[i] = tp1f_special_fvtx_tracks_aa->GetBinError(index_b);
      econtent_aa0[i] = tp1f_special_fvtx_tracks_aa->GetBinError(index_a);
      // --- apply the acceptance corrections
      content_ab0[i] -= (cosa*cosb + sina*sinb);
      content_bb0[i] -= (cosb*cosb + sinb*sinb);
      content_aa0[i] -= (cosa*cosa + sina*sina);
      // --- calculate the final ratio
      content_final0[i] = -9;
      if ( content_aa0[i] > 0 && content_bb0[i] > 0 )
	{
	  content_final0[i] = content_ab0[i]/sqrt(content_aa0[i]*content_bb0[i]);
	  econtent_final0[i] = content_final0[i] * sqrt(
							(econtent_ab0[i]/content_ab0[i])*(econtent_ab0[i]/content_ab0[i]) +
							(econtent_aa0[i]/content_aa0[i])*(econtent_aa0[i]/content_aa0[i]) +
							(econtent_bb0[i]/content_bb0[i])*(econtent_bb0[i]/content_bb0[i])
							);
	  cout << content_final0[i] << " " << econtent_final0[i] << endl;
	}
    }
  // -------------------------------------------------------------------------------------------
  TGraphErrors* tge_ab0 = new TGraphErrors(12,x,content_ab0,0,econtent_ab0);
  tge_ab0->SetTitle(Form("d+Au collisions at #sqrt{s_{NN}} = %d GeV",energy));
  tge_ab0->SetMarkerStyle(kOpenSquare);
  tge_ab0->SetMarkerColor(kBlack);
  tge_ab0->Draw("ap");
  tge_ab0->SetMaximum(0.05);
  tge_ab0->SetMinimum(0.0);
  tge_ab0->GetXaxis()->SetLimits(-3,3);
  tge_ab0->GetXaxis()->SetTitle("#eta_{b}");
  tge_ab0->GetYaxis()->SetTitle("v_{a,b}(#eta_{a} = -2.75)");
  tge_ab0->GetYaxis()->SetTitleOffset(1.25);
  c1->Print("FigsDecorrelation/early_ab0.png");
  c1->Print("FigsDecorrelation/early_ab0.pdf");
  TGraphErrors* tge_aa0 = new TGraphErrors(12,x,content_aa0,0,econtent_aa0);
  tge_aa0->SetTitle(Form("d+Au collisions at #sqrt{s_{NN}} = %d GeV",energy));
  tge_aa0->SetMarkerStyle(kOpenSquare);
  tge_aa0->SetMarkerColor(kBlack);
  tge_aa0->Draw("ap");
  tge_aa0->SetMaximum(0.05);
  tge_aa0->SetMinimum(0.0);
  tge_aa0->GetXaxis()->SetLimits(-3,3);
  tge_aa0->GetXaxis()->SetTitle("#eta_{b}");
  tge_aa0->GetYaxis()->SetTitle("v_{a,a}(#eta_{a} = -2.75)");
  tge_aa0->GetYaxis()->SetTitleOffset(1.25);
  c1->Print("FigsDecorrelation/early_aa0.png");
  c1->Print("FigsDecorrelation/early_aa0.pdf");
  TGraphErrors* tge_bb0 = new TGraphErrors(12,x,content_bb0,0,econtent_bb0);
  tge_bb0->SetTitle(Form("d+Au collisions at #sqrt{s_{NN}} = %d GeV",energy));
  tge_bb0->SetMarkerStyle(kOpenSquare);
  tge_bb0->SetMarkerColor(kBlack);
  tge_bb0->Draw("ap");
  tge_bb0->SetMaximum(0.05);
  tge_bb0->SetMinimum(0.0);
  tge_bb0->GetXaxis()->SetLimits(-3,3);
  tge_bb0->GetXaxis()->SetTitle("#eta_{b}");
  tge_bb0->GetYaxis()->SetTitle("v_{b,b}(#eta_{a} = -2.75)");
  tge_bb0->GetYaxis()->SetTitleOffset(1.25);
  c1->Print("FigsDecorrelation/early_bb0.png");
  c1->Print("FigsDecorrelation/early_bb0.pdf");
  TGraphErrors* tge_final0 = new TGraphErrors(12,x,content_final0,0,econtent_final0);
  tge_final0->SetTitle(Form("d+Au collisions at #sqrt{s_{NN}} = %d GeV",energy));
  tge_final0->SetMarkerStyle(kOpenSquare);
  tge_final0->SetMarkerColor(kBlack);
  tge_final0->Draw("ap");
  tge_final0->SetMaximum(1.5);
  tge_final0->SetMinimum(0.0);
  tge_final0->GetXaxis()->SetLimits(-3,3);
  tge_final0->GetXaxis()->SetTitle("#eta_{b}");
  tge_final0->GetYaxis()->SetTitle("r(#eta_{a} = -2.75)");
  c1->Print("FigsDecorrelation/early_final0.png");
  c1->Print("FigsDecorrelation/early_final0.pdf");
  // -------------------------------------------------------------------------------------------



  // -------------------------------------------------------------------------------------------
  // --- come back here for bin 1
  // ---
  double content_ab1[12];
  double content_aa1[12];
  double content_bb1[12];
  double content_final1[12];
  double econtent_ab1[12];
  double econtent_aa1[12];
  double econtent_bb1[12];
  double econtent_final1[12];
  double x[12];
  // -------------------------------------------------------------------------------------------
  for ( int i = 0; i < 12; ++i )
    {
      // --- set the indices for particles a and b
      int index_a = 2;
      int index_b = i+1;
      // --- get the acceptance correction terms
      double cosb = tp1f_special_fvtx_tracks_aa_cos->GetBinContent(index_b);
      double cosa = tp1f_special_fvtx_tracks_aa_cos->GetBinContent(index_a);
      double sinb = tp1f_special_fvtx_tracks_aa_sin->GetBinContent(index_b);
      double sina = tp1f_special_fvtx_tracks_aa_sin->GetBinContent(index_a);
      // --- get eta_b
      x[i] = tp1f_special_fvtx_tracks_ab1->GetBinCenter(i+1);
      // --- get the components
      content_ab1[i] = tp1f_special_fvtx_tracks_ab1->GetBinContent(index_b);
      content_bb1[i] = tp1f_special_fvtx_tracks_aa->GetBinContent(index_b);
      content_aa1[i] = tp1f_special_fvtx_tracks_aa->GetBinContent(index_a);
      econtent_ab1[i] = tp1f_special_fvtx_tracks_ab1->GetBinError(index_b);
      econtent_bb1[i] = tp1f_special_fvtx_tracks_aa->GetBinError(index_b);
      econtent_aa1[i] = tp1f_special_fvtx_tracks_aa->GetBinError(index_a);
      // --- apply the acceptance corrections
      content_ab1[i] -= (cosa*cosb + sina*sinb);
      content_bb1[i] -= (cosb*cosb + sinb*sinb);
      content_aa1[i] -= (cosa*cosa + sina*sina);
      // --- calculate the final ratio
      content_final1[i] = -9;
      if ( content_aa1[i] > 0 && content_bb1[i] > 0 )
	{
	  content_final1[i] = content_ab1[i]/sqrt(content_aa1[i]*content_bb1[i]);
	  econtent_final1[i] = content_final1[i] * sqrt(
							(econtent_ab1[i]/content_ab1[i])*(econtent_ab1[i]/content_ab1[i]) +
							(econtent_aa1[i]/content_aa1[i])*(econtent_aa1[i]/content_aa1[i]) +
							(econtent_bb1[i]/content_bb1[i])*(econtent_bb1[i]/content_bb1[i])
							);
	  cout << content_final1[i] << " " << econtent_final1[i] << endl;
	}
    }
  // -------------------------------------------------------------------------------------------
  TGraphErrors* tge_ab1 = new TGraphErrors(12,x,content_ab1,0,econtent_ab1);
  tge_ab1->SetTitle(Form("d+Au collisions at #sqrt{s_{NN}} = %d GeV",energy));
  tge_ab1->SetMarkerStyle(kOpenSquare);
  tge_ab1->SetMarkerColor(kBlack);
  tge_ab1->Draw("ap");
  tge_ab1->SetMaximum(0.05);
  tge_ab1->SetMinimum(0.0);
  tge_ab1->GetXaxis()->SetLimits(-3,3);
  tge_ab1->GetXaxis()->SetTitle("#eta_{b}");
  tge_ab1->GetYaxis()->SetTitle("v_{a,b}(#eta_{a} = -2.25)");
  tge_ab1->GetYaxis()->SetTitleOffset(1.25);
  c1->Print("FigsDecorrelation/early_ab1.png");
  c1->Print("FigsDecorrelation/early_ab1.pdf");
  TGraphErrors* tge_aa1 = new TGraphErrors(12,x,content_aa1,0,econtent_aa1);
  tge_aa1->SetTitle(Form("d+Au collisions at #sqrt{s_{NN}} = %d GeV",energy));
  tge_aa1->SetMarkerStyle(kOpenSquare);
  tge_aa1->SetMarkerColor(kBlack);
  tge_aa1->Draw("ap");
  tge_aa1->SetMaximum(0.05);
  tge_aa1->SetMinimum(0.0);
  tge_aa1->GetXaxis()->SetLimits(-3,3);
  tge_aa1->GetXaxis()->SetTitle("#eta_{b}");
  tge_aa1->GetYaxis()->SetTitle("v_{a,a}(#eta_{a} = -2.25)");
  tge_aa1->GetYaxis()->SetTitleOffset(1.25);
  c1->Print("FigsDecorrelation/early_aa1.png");
  c1->Print("FigsDecorrelation/early_aa1.pdf");
  TGraphErrors* tge_bb1 = new TGraphErrors(12,x,content_bb1,0,econtent_bb1);
  tge_bb1->SetTitle(Form("d+Au collisions at #sqrt{s_{NN}} = %d GeV",energy));
  tge_bb1->SetMarkerStyle(kOpenSquare);
  tge_bb1->SetMarkerColor(kBlack);
  tge_bb1->Draw("ap");
  tge_bb1->SetMaximum(0.05);
  tge_bb1->SetMinimum(0.0);
  tge_bb1->GetXaxis()->SetLimits(-3,3);
  tge_bb1->GetXaxis()->SetTitle("#eta_{b}");
  tge_bb1->GetYaxis()->SetTitle("v_{b,b}(#eta_{a} = -2.25)");
  tge_bb1->GetYaxis()->SetTitleOffset(1.25);
  c1->Print("FigsDecorrelation/early_bb1.png");
  c1->Print("FigsDecorrelation/early_bb1.pdf");
  TGraphErrors* tge_final1 = new TGraphErrors(12,x,content_final1,0,econtent_final1);
  tge_final1->SetTitle(Form("d+Au collisions at #sqrt{s_{NN}} = %d GeV",energy));
  tge_final1->SetMarkerStyle(kOpenSquare);
  tge_final1->SetMarkerColor(kBlack);
  tge_final1->Draw("ap");
  tge_final1->SetMaximum(1.5);
  tge_final1->SetMinimum(0.0);
  tge_final1->GetXaxis()->SetLimits(-3,3);
  tge_final1->GetXaxis()->SetTitle("#eta_{b}");
  tge_final1->GetYaxis()->SetTitle("r(#eta_{a} = -2.25)");
  c1->Print("FigsDecorrelation/early_final1.png");
  c1->Print("FigsDecorrelation/early_final1.pdf");
  // -------------------------------------------------------------------------------------------



  // -------------------------------------------------------------------------------------------
  // --- come back here for bin 1
  // ---
  double content_ab2[12];
  double content_aa2[12];
  double content_bb2[12];
  double content_final2[12];
  double econtent_ab2[12];
  double econtent_aa2[12];
  double econtent_bb2[12];
  double econtent_final2[12];
  double x[12];
  // -------------------------------------------------------------------------------------------
  for ( int i = 0; i < 12; ++i )
    {
      // --- set the indices for particles a and b
      int index_a = 3;
      int index_b = i+1;
      // --- get the acceptance correction terms
      double cosb = tp1f_special_fvtx_tracks_aa_cos->GetBinContent(index_b);
      double cosa = tp1f_special_fvtx_tracks_aa_cos->GetBinContent(index_a);
      double sinb = tp1f_special_fvtx_tracks_aa_sin->GetBinContent(index_b);
      double sina = tp1f_special_fvtx_tracks_aa_sin->GetBinContent(index_a);
      // --- get eta_b
      x[i] = tp1f_special_fvtx_tracks_ab2->GetBinCenter(i+1);
      // --- get the components
      content_ab2[i] = tp1f_special_fvtx_tracks_ab2->GetBinContent(index_b);
      content_bb2[i] = tp1f_special_fvtx_tracks_aa->GetBinContent(index_b);
      content_aa2[i] = tp1f_special_fvtx_tracks_aa->GetBinContent(index_a);
      econtent_ab2[i] = tp1f_special_fvtx_tracks_ab2->GetBinError(index_b);
      econtent_bb2[i] = tp1f_special_fvtx_tracks_aa->GetBinError(index_b);
      econtent_aa2[i] = tp1f_special_fvtx_tracks_aa->GetBinError(index_a);
      // --- apply the acceptance corrections
      content_ab2[i] -= (cosa*cosb + sina*sinb);
      content_bb2[i] -= (cosb*cosb + sinb*sinb);
      content_aa2[i] -= (cosa*cosa + sina*sina);
      // --- calculate the final ratio
      content_final2[i] = -9;
      if ( content_aa2[i] > 0 && content_bb2[i] > 0 )
	{
	  content_final2[i] = content_ab2[i]/sqrt(content_aa2[i]*content_bb2[i]);
	  econtent_final2[i] = content_final2[i] * sqrt(
							(econtent_ab2[i]/content_ab2[i])*(econtent_ab2[i]/content_ab2[i]) +
							(econtent_aa2[i]/content_aa2[i])*(econtent_aa2[i]/content_aa2[i]) +
							(econtent_bb2[i]/content_bb2[i])*(econtent_bb2[i]/content_bb2[i])
							);
	  cout << content_final2[i] << " " << econtent_final2[i] << endl;
	}
    }
  // -------------------------------------------------------------------------------------------
  TGraphErrors* tge_ab2 = new TGraphErrors(12,x,content_ab2,0,econtent_ab2);
  tge_ab2->SetTitle(Form("d+Au collisions at #sqrt{s_{NN}} = %d GeV",energy));
  tge_ab2->SetMarkerStyle(kOpenSquare);
  tge_ab2->SetMarkerColor(kBlack);
  tge_ab2->Draw("ap");
  tge_ab2->SetMaximum(0.05);
  tge_ab2->SetMinimum(0.0);
  tge_ab2->GetXaxis()->SetLimits(-3,3);
  tge_ab2->GetXaxis()->SetTitle("#eta_{b}");
  tge_ab2->GetYaxis()->SetTitle("v_{a,b}(#eta_{a} = -1.75)");
  tge_ab2->GetYaxis()->SetTitleOffset(1.25);
  c1->Print("FigsDecorrelation/early_ab2.png");
  c1->Print("FigsDecorrelation/early_ab2.pdf");
  TGraphErrors* tge_aa2 = new TGraphErrors(12,x,content_aa2,0,econtent_aa2);
  tge_aa2->SetTitle(Form("d+Au collisions at #sqrt{s_{NN}} = %d GeV",energy));
  tge_aa2->SetMarkerStyle(kOpenSquare);
  tge_aa2->SetMarkerColor(kBlack);
  tge_aa2->Draw("ap");
  tge_aa2->SetMaximum(0.05);
  tge_aa2->SetMinimum(0.0);
  tge_aa2->GetXaxis()->SetLimits(-3,3);
  tge_aa2->GetXaxis()->SetTitle("#eta_{b}");
  tge_aa2->GetYaxis()->SetTitle("v_{a,a}(#eta_{a} = -1.75)");
  tge_aa2->GetYaxis()->SetTitleOffset(1.25);
  c1->Print("FigsDecorrelation/early_aa2.png");
  c1->Print("FigsDecorrelation/early_aa2.pdf");
  TGraphErrors* tge_bb2 = new TGraphErrors(12,x,content_bb2,0,econtent_bb2);
  tge_bb2->SetTitle(Form("d+Au collisions at #sqrt{s_{NN}} = %d GeV",energy));
  tge_bb2->SetMarkerStyle(kOpenSquare);
  tge_bb2->SetMarkerColor(kBlack);
  tge_bb2->Draw("ap");
  tge_bb2->SetMaximum(0.05);
  tge_bb2->SetMinimum(0.0);
  tge_bb2->GetXaxis()->SetLimits(-3,3);
  tge_bb2->GetXaxis()->SetTitle("#eta_{b}");
  tge_bb2->GetYaxis()->SetTitle("v_{b,b}(#eta_{a} = -1.75)");
  tge_bb2->GetYaxis()->SetTitleOffset(1.25);
  c1->Print("FigsDecorrelation/early_bb2.png");
  c1->Print("FigsDecorrelation/early_bb2.pdf");
  TGraphErrors* tge_final2 = new TGraphErrors(12,x,content_final2,0,econtent_final2);
  tge_final2->SetTitle(Form("d+Au collisions at #sqrt{s_{NN}} = %d GeV",energy));
  tge_final2->SetMarkerStyle(kOpenSquare);
  tge_final2->SetMarkerColor(kBlack);
  tge_final2->Draw("ap");
  tge_final2->SetMaximum(1.5);
  tge_final2->SetMinimum(0.0);
  tge_final2->GetXaxis()->SetLimits(-3,3);
  tge_final2->GetXaxis()->SetTitle("#eta_{b}");
  tge_final2->GetYaxis()->SetTitle("r(#eta_{a} = -1.75)");
  c1->Print("FigsDecorrelation/early_final2.png");
  c1->Print("FigsDecorrelation/early_final2.pdf");
  // -------------------------------------------------------------------------------------------



  // -------------------------------------------------------------------------------------------
  // --- come back here for bin 3
  // ---
  double content_ab3[12];
  double content_aa3[12];
  double content_bb3[12];
  double content_final3[12];
  double econtent_ab3[12];
  double econtent_aa3[12];
  double econtent_bb3[12];
  double econtent_final3[12];
  double x[12];
  // -------------------------------------------------------------------------------------------
  for ( int i = 0; i < 12; ++i )
    {
      // --- set the indices for particles a and b
      int index_a = 4;
      int index_b = i+1;
      // --- get the acceptance correction terms
      double cosb = tp1f_special_fvtx_tracks_aa_cos->GetBinContent(index_b);
      double cosa = tp1f_special_fvtx_tracks_aa_cos->GetBinContent(index_a);
      double sinb = tp1f_special_fvtx_tracks_aa_sin->GetBinContent(index_b);
      double sina = tp1f_special_fvtx_tracks_aa_sin->GetBinContent(index_a);
      // --- get eta_b
      x[i] = tp1f_special_fvtx_tracks_ab3->GetBinCenter(i+1);
      // --- get the components
      content_ab3[i] = tp1f_special_fvtx_tracks_ab3->GetBinContent(index_b);
      content_bb3[i] = tp1f_special_fvtx_tracks_aa->GetBinContent(index_b);
      content_aa3[i] = tp1f_special_fvtx_tracks_aa->GetBinContent(index_a);
      econtent_ab3[i] = tp1f_special_fvtx_tracks_ab3->GetBinError(index_b);
      econtent_bb3[i] = tp1f_special_fvtx_tracks_aa->GetBinError(index_b);
      econtent_aa3[i] = tp1f_special_fvtx_tracks_aa->GetBinError(index_a);
      // --- apply the acceptance corrections
      content_ab3[i] -= (cosa*cosb + sina*sinb);
      content_bb3[i] -= (cosb*cosb + sinb*sinb);
      content_aa3[i] -= (cosa*cosa + sina*sina);
      // --- calculate the final ratio
      content_final3[i] = -9;
      if ( content_aa3[i] > 0 && content_bb3[i] > 0 )
	{
	  content_final3[i] = content_ab3[i]/sqrt(content_aa3[i]*content_bb3[i]);
	  econtent_final3[i] = content_final3[i] * sqrt(
							(econtent_ab3[i]/content_ab3[i])*(econtent_ab3[i]/content_ab3[i]) +
							(econtent_aa3[i]/content_aa3[i])*(econtent_aa3[i]/content_aa3[i]) +
							(econtent_bb3[i]/content_bb3[i])*(econtent_bb3[i]/content_bb3[i])
							);
	  cout << content_final3[i] << " " << econtent_final3[i] << endl;
	}
    }
  // -------------------------------------------------------------------------------------------
  TGraphErrors* tge_ab3 = new TGraphErrors(12,x,content_ab3,0,econtent_ab3);
  tge_ab3->SetTitle(Form("d+Au collisions at #sqrt{s_{NN}} = %d GeV",energy));
  tge_ab3->SetMarkerStyle(kOpenSquare);
  tge_ab3->SetMarkerColor(kBlack);
  tge_ab3->Draw("ap");
  tge_ab3->SetMaximum(0.05);
  tge_ab3->SetMinimum(0.0);
  tge_ab3->GetXaxis()->SetLimits(-3,3);
  tge_ab3->GetXaxis()->SetTitle("#eta_{b}");
  tge_ab3->GetYaxis()->SetTitle("v_{a,b}(#eta_{a} = -1.25)");
  tge_ab3->GetYaxis()->SetTitleOffset(1.25);
  c1->Print("FigsDecorrelation/early_ab3.png");
  c1->Print("FigsDecorrelation/early_ab3.pdf");
  TGraphErrors* tge_aa3 = new TGraphErrors(12,x,content_aa3,0,econtent_aa3);
  tge_aa3->SetTitle(Form("d+Au collisions at #sqrt{s_{NN}} = %d GeV",energy));
  tge_aa3->SetMarkerStyle(kOpenSquare);
  tge_aa3->SetMarkerColor(kBlack);
  tge_aa3->Draw("ap");
  tge_aa3->SetMaximum(0.05);
  tge_aa3->SetMinimum(0.0);
  tge_aa3->GetXaxis()->SetLimits(-3,3);
  tge_aa3->GetXaxis()->SetTitle("#eta_{b}");
  tge_aa3->GetYaxis()->SetTitle("v_{a,a}(#eta_{a} = -1.25)");
  tge_aa3->GetYaxis()->SetTitleOffset(1.25);
  c1->Print("FigsDecorrelation/early_aa3.png");
  c1->Print("FigsDecorrelation/early_aa3.pdf");
  TGraphErrors* tge_bb3 = new TGraphErrors(12,x,content_bb3,0,econtent_bb3);
  tge_bb3->SetTitle(Form("d+Au collisions at #sqrt{s_{NN}} = %d GeV",energy));
  tge_bb3->SetMarkerStyle(kOpenSquare);
  tge_bb3->SetMarkerColor(kBlack);
  tge_bb3->Draw("ap");
  tge_bb3->SetMaximum(0.05);
  tge_bb3->SetMinimum(0.0);
  tge_bb3->GetXaxis()->SetLimits(-3,3);
  tge_bb3->GetXaxis()->SetTitle("#eta_{b}");
  tge_bb3->GetYaxis()->SetTitle("v_{b,b}(#eta_{a} = -1.25)");
  tge_bb3->GetYaxis()->SetTitleOffset(1.25);
  c1->Print("FigsDecorrelation/early_bb3.png");
  c1->Print("FigsDecorrelation/early_bb3.pdf");
  TGraphErrors* tge_final3 = new TGraphErrors(12,x,content_final3,0,econtent_final3);
  tge_final3->SetTitle(Form("d+Au collisions at #sqrt{s_{NN}} = %d GeV",energy));
  tge_final3->SetMarkerStyle(kOpenSquare);
  tge_final3->SetMarkerColor(kBlack);
  tge_final3->Draw("ap");
  tge_final3->SetMaximum(1.5);
  tge_final3->SetMinimum(0.0);
  tge_final3->GetXaxis()->SetLimits(-3,3);
  tge_final3->GetXaxis()->SetTitle("#eta_{b}");
  tge_final3->GetYaxis()->SetTitle("r(#eta_{a} = -1.25)");
  c1->Print("FigsDecorrelation/early_final3.png");
  c1->Print("FigsDecorrelation/early_final3.pdf");
  // -------------------------------------------------------------------------------------------



  // -------------------------------------------------------------------------------------------
  // --- come back here for bin 4 (north side)
  // ---
  double content_ab4[12];
  double content_aa4[12];
  double content_bb4[12];
  double content_final4[12];
  double econtent_ab4[12];
  double econtent_aa4[12];
  double econtent_bb4[12];
  double econtent_final4[12];
  double x[12];
  // -------------------------------------------------------------------------------------------
  for ( int i = 0; i < 12; ++i )
    {
      // --- set the indices for particles a and b
      int index_a = 9; // offset yby 4 because of the 4 empty bins outside the acceptance
      int index_b = i+1;
      // --- get the acceptance correction terms
      double cosb = tp1f_special_fvtx_tracks_aa_cos->GetBinContent(index_b);
      double cosa = tp1f_special_fvtx_tracks_aa_cos->GetBinContent(index_a);
      double sinb = tp1f_special_fvtx_tracks_aa_sin->GetBinContent(index_b);
      double sina = tp1f_special_fvtx_tracks_aa_sin->GetBinContent(index_a);
      // --- get eta_b
      x[i] = tp1f_special_fvtx_tracks_ab4->GetBinCenter(i+1);
      // --- get the components
      content_ab4[i] = tp1f_special_fvtx_tracks_ab4->GetBinContent(index_b);
      content_bb4[i] = tp1f_special_fvtx_tracks_aa->GetBinContent(index_b);
      content_aa4[i] = tp1f_special_fvtx_tracks_aa->GetBinContent(index_a);
      econtent_ab4[i] = tp1f_special_fvtx_tracks_ab4->GetBinError(index_b);
      econtent_bb4[i] = tp1f_special_fvtx_tracks_aa->GetBinError(index_b);
      econtent_aa4[i] = tp1f_special_fvtx_tracks_aa->GetBinError(index_a);
      // --- apply the acceptance corrections
      content_ab4[i] -= (cosa*cosb + sina*sinb);
      content_bb4[i] -= (cosb*cosb + sinb*sinb);
      content_aa4[i] -= (cosa*cosa + sina*sina);
      // --- calculate the final ratio
      content_final4[i] = -9;
      if ( content_aa4[i] > 0 && content_bb4[i] > 0 )
	{
	  content_final4[i] = content_ab4[i]/sqrt(content_aa4[i]*content_bb4[i]);
	  econtent_final4[i] = content_final4[i] * sqrt(
							(econtent_ab4[i]/content_ab4[i])*(econtent_ab4[i]/content_ab4[i]) +
							(econtent_aa4[i]/content_aa4[i])*(econtent_aa4[i]/content_aa4[i]) +
							(econtent_bb4[i]/content_bb4[i])*(econtent_bb4[i]/content_bb4[i])
							);
	}
      cout << content_final4[i] << " " << econtent_final4[i] << endl;
      // cout << content_ab4[i] << " " << econtent_ab4[i] << endl;
      // cout << content_aa4[i] << " " << econtent_aa4[i] << endl;
      // cout << content_bb4[i] << " " << econtent_bb4[i] << endl;
    }
  // -------------------------------------------------------------------------------------------
  TGraphErrors* tge_ab4 = new TGraphErrors(12,x,content_ab4,0,econtent_ab4);
  tge_ab4->SetTitle(Form("d+Au collisions at #sqrt{s_{NN}} = %d GeV",energy));
  tge_ab4->SetMarkerStyle(kOpenSquare);
  tge_ab4->SetMarkerColor(kBlack);
  tge_ab4->Draw("ap");
  tge_ab4->SetMaximum(0.05);
  tge_ab4->SetMinimum(0.0);
  tge_ab4->GetXaxis()->SetLimits(-3,3);
  tge_ab4->GetXaxis()->SetTitle("#eta_{b}");
  tge_ab4->GetYaxis()->SetTitle("v_{a,b}(#eta_{a} = 1.25)");
  tge_ab4->GetYaxis()->SetTitleOffset(1.25);
  c1->Print("FigsDecorrelation/early_ab4.png");
  c1->Print("FigsDecorrelation/early_ab4.pdf");
  TGraphErrors* tge_aa4 = new TGraphErrors(12,x,content_aa4,0,econtent_aa4);
  tge_aa4->SetTitle(Form("d+Au collisions at #sqrt{s_{NN}} = %d GeV",energy));
  tge_aa4->SetMarkerStyle(kOpenSquare);
  tge_aa4->SetMarkerColor(kBlack);
  tge_aa4->Draw("ap");
  tge_aa4->SetMaximum(0.05);
  tge_aa4->SetMinimum(0.0);
  tge_aa4->GetXaxis()->SetLimits(-3,3);
  tge_aa4->GetXaxis()->SetTitle("#eta_{b}");
  tge_aa4->GetYaxis()->SetTitle("v_{a,a}(#eta_{a} = 1.25)");
  tge_aa4->GetYaxis()->SetTitleOffset(1.25);
  c1->Print("FigsDecorrelation/early_aa4.png");
  c1->Print("FigsDecorrelation/early_aa4.pdf");
  TGraphErrors* tge_bb4 = new TGraphErrors(12,x,content_bb4,0,econtent_bb4);
  tge_bb4->SetTitle(Form("d+Au collisions at #sqrt{s_{NN}} = %d GeV",energy));
  tge_bb4->SetMarkerStyle(kOpenSquare);
  tge_bb4->SetMarkerColor(kBlack);
  tge_bb4->Draw("ap");
  tge_bb4->SetMaximum(0.05);
  tge_bb4->SetMinimum(0.0);
  tge_bb4->GetXaxis()->SetLimits(-3,3);
  tge_bb4->GetXaxis()->SetTitle("#eta_{b}");
  tge_bb4->GetYaxis()->SetTitle("v_{b,b}(#eta_{a} = 1.25)");
  tge_bb4->GetYaxis()->SetTitleOffset(1.25);
  c1->Print("FigsDecorrelation/early_bb4.png");
  c1->Print("FigsDecorrelation/early_bb4.pdf");
  TGraphErrors* tge_final4 = new TGraphErrors(12,x,content_final4,0,econtent_final4);
  tge_final4->SetTitle(Form("d+Au collisions at #sqrt{s_{NN}} = %d GeV",energy));
  tge_final4->SetMarkerStyle(kOpenSquare);
  tge_final4->SetMarkerColor(kBlack);
  tge_final4->Draw("ap");
  tge_final4->SetMaximum(1.5);
  tge_final4->SetMinimum(0.0);
  tge_final4->GetXaxis()->SetLimits(-3,3);
  tge_final4->GetXaxis()->SetTitle("#eta_{b}");
  tge_final4->GetYaxis()->SetTitle("r(#eta_{a} = 1.25)");
  c1->Print("FigsDecorrelation/early_final4.png");
  c1->Print("FigsDecorrelation/early_final4.pdf");
  // -------------------------------------------------------------------------------------------



  // -------------------------------------------------------------------------------------------
  // --- come back here for bin 5 (north side)
  // ---
  double content_ab5[12];
  double content_aa5[12];
  double content_bb5[12];
  double content_final5[12];
  double econtent_ab5[12];
  double econtent_aa5[12];
  double econtent_bb5[12];
  double econtent_final5[12];
  double x[12];
  // -------------------------------------------------------------------------------------------
  for ( int i = 0; i < 12; ++i )
    {
      // --- set the indices for particles a and b
      int index_a = 10; // offset yby 4 because of the 4 empty bins outside the acceptance
      int index_b = i+1;
      // --- get the acceptance correction terms
      double cosb = tp1f_special_fvtx_tracks_aa_cos->GetBinContent(index_b);
      double cosa = tp1f_special_fvtx_tracks_aa_cos->GetBinContent(index_a);
      double sinb = tp1f_special_fvtx_tracks_aa_sin->GetBinContent(index_b);
      double sina = tp1f_special_fvtx_tracks_aa_sin->GetBinContent(index_a);
      // --- get eta_b
      x[i] = tp1f_special_fvtx_tracks_ab5->GetBinCenter(i+1);
      // --- get the components
      content_ab5[i] = tp1f_special_fvtx_tracks_ab5->GetBinContent(index_b);
      content_bb5[i] = tp1f_special_fvtx_tracks_aa->GetBinContent(index_b);
      content_aa5[i] = tp1f_special_fvtx_tracks_aa->GetBinContent(index_a);
      econtent_ab5[i] = tp1f_special_fvtx_tracks_ab5->GetBinError(index_b);
      econtent_bb5[i] = tp1f_special_fvtx_tracks_aa->GetBinError(index_b);
      econtent_aa5[i] = tp1f_special_fvtx_tracks_aa->GetBinError(index_a);
      // --- apply the acceptance corrections
      content_ab5[i] -= (cosa*cosb + sina*sinb);
      content_bb5[i] -= (cosb*cosb + sinb*sinb);
      content_aa5[i] -= (cosa*cosa + sina*sina);
      // --- calculate the final ratio
      content_final5[i] = -9;
      if ( content_aa5[i] > 0 && content_bb5[i] > 0 )
	{
	  content_final5[i] = content_ab5[i]/sqrt(content_aa5[i]*content_bb5[i]);
	  econtent_final5[i] = content_final5[i] * sqrt(
							(econtent_ab5[i]/content_ab5[i])*(econtent_ab5[i]/content_ab5[i]) +
							(econtent_aa5[i]/content_aa5[i])*(econtent_aa5[i]/content_aa5[i]) +
							(econtent_bb5[i]/content_bb5[i])*(econtent_bb5[i]/content_bb5[i])
							);
	}
      cout << content_final5[i] << " " << econtent_final5[i] << endl;
      // cout << content_ab5[i] << " " << econtent_ab5[i] << endl;
      // cout << content_aa5[i] << " " << econtent_aa5[i] << endl;
      // cout << content_bb5[i] << " " << econtent_bb5[i] << endl;
    }
  // -------------------------------------------------------------------------------------------
  TGraphErrors* tge_ab5 = new TGraphErrors(12,x,content_ab5,0,econtent_ab5);
  tge_ab5->SetTitle(Form("d+Au collisions at #sqrt{s_{NN}} = %d GeV",energy));
  tge_ab5->SetMarkerStyle(kOpenSquare);
  tge_ab5->SetMarkerColor(kBlack);
  tge_ab5->Draw("ap");
  tge_ab5->SetMaximum(0.05);
  tge_ab5->SetMinimum(0.0);
  tge_ab5->GetXaxis()->SetLimits(-3,3);
  tge_ab5->GetXaxis()->SetTitle("#eta_{b}");
  tge_ab5->GetYaxis()->SetTitle("v_{a,b}(#eta_{a} = 1.75)");
  tge_ab5->GetYaxis()->SetTitleOffset(1.25);
  c1->Print("FigsDecorrelation/early_ab5.png");
  c1->Print("FigsDecorrelation/early_ab5.pdf");
  TGraphErrors* tge_aa5 = new TGraphErrors(12,x,content_aa5,0,econtent_aa5);
  tge_aa5->SetTitle(Form("d+Au collisions at #sqrt{s_{NN}} = %d GeV",energy));
  tge_aa5->SetMarkerStyle(kOpenSquare);
  tge_aa5->SetMarkerColor(kBlack);
  tge_aa5->Draw("ap");
  tge_aa5->SetMaximum(0.05);
  tge_aa5->SetMinimum(0.0);
  tge_aa5->GetXaxis()->SetLimits(-3,3);
  tge_aa5->GetXaxis()->SetTitle("#eta_{b}");
  tge_aa5->GetYaxis()->SetTitle("v_{a,a}(#eta_{a} = 1.75)");
  tge_aa5->GetYaxis()->SetTitleOffset(1.25);
  c1->Print("FigsDecorrelation/early_aa5.png");
  c1->Print("FigsDecorrelation/early_aa5.pdf");
  TGraphErrors* tge_bb5 = new TGraphErrors(12,x,content_bb5,0,econtent_bb5);
  tge_bb5->SetTitle(Form("d+Au collisions at #sqrt{s_{NN}} = %d GeV",energy));
  tge_bb5->SetMarkerStyle(kOpenSquare);
  tge_bb5->SetMarkerColor(kBlack);
  tge_bb5->Draw("ap");
  tge_bb5->SetMaximum(0.05);
  tge_bb5->SetMinimum(0.0);
  tge_bb5->GetXaxis()->SetLimits(-3,3);
  tge_bb5->GetXaxis()->SetTitle("#eta_{b}");
  tge_bb5->GetYaxis()->SetTitle("v_{b,b}(#eta_{a} = 1.75)");
  tge_bb5->GetYaxis()->SetTitleOffset(1.25);
  c1->Print("FigsDecorrelation/early_bb5.png");
  c1->Print("FigsDecorrelation/early_bb5.pdf");
  TGraphErrors* tge_final5 = new TGraphErrors(12,x,content_final5,0,econtent_final5);
  tge_final5->SetTitle(Form("d+Au collisions at #sqrt{s_{NN}} = %d GeV",energy));
  tge_final5->SetMarkerStyle(kOpenSquare);
  tge_final5->SetMarkerColor(kBlack);
  tge_final5->Draw("ap");
  tge_final5->SetMaximum(1.5);
  tge_final5->SetMinimum(0.0);
  tge_final5->GetXaxis()->SetLimits(-3,3);
  tge_final5->GetXaxis()->SetTitle("#eta_{b}");
  tge_final5->GetYaxis()->SetTitle("r(#eta_{a} = 1.75)");
  c1->Print("FigsDecorrelation/early_final5.png");
  c1->Print("FigsDecorrelation/early_final5.pdf");
  // -------------------------------------------------------------------------------------------



  // -------------------------------------------------------------------------------------------
  // --- come back here for bin 6 (north side)
  // ---
  double content_ab6[12];
  double content_aa6[12];
  double content_bb6[12];
  double content_final6[12];
  double econtent_ab6[12];
  double econtent_aa6[12];
  double econtent_bb6[12];
  double econtent_final6[12];
  double x[12];
  // -------------------------------------------------------------------------------------------
  for ( int i = 0; i < 12; ++i )
    {
      // --- set the indices for particles a and b
      int index_a = 11; // offset yby 4 because of the 4 empty bins outside the acceptance
      int index_b = i+1;
      // --- get the acceptance correction terms
      double cosb = tp1f_special_fvtx_tracks_aa_cos->GetBinContent(index_b);
      double cosa = tp1f_special_fvtx_tracks_aa_cos->GetBinContent(index_a);
      double sinb = tp1f_special_fvtx_tracks_aa_sin->GetBinContent(index_b);
      double sina = tp1f_special_fvtx_tracks_aa_sin->GetBinContent(index_a);
      // --- get eta_b
      x[i] = tp1f_special_fvtx_tracks_ab6->GetBinCenter(i+1);
      // --- get the components
      content_ab6[i] = tp1f_special_fvtx_tracks_ab6->GetBinContent(index_b);
      content_bb6[i] = tp1f_special_fvtx_tracks_aa->GetBinContent(index_b);
      content_aa6[i] = tp1f_special_fvtx_tracks_aa->GetBinContent(index_a);
      econtent_ab6[i] = tp1f_special_fvtx_tracks_ab6->GetBinError(index_b);
      econtent_bb6[i] = tp1f_special_fvtx_tracks_aa->GetBinError(index_b);
      econtent_aa6[i] = tp1f_special_fvtx_tracks_aa->GetBinError(index_a);
      // --- apply the acceptance corrections
      content_ab6[i] -= (cosa*cosb + sina*sinb);
      content_bb6[i] -= (cosb*cosb + sinb*sinb);
      content_aa6[i] -= (cosa*cosa + sina*sina);
      // --- calculate the final ratio
      content_final6[i] = -9;
      if ( content_aa6[i] > 0 && content_bb6[i] > 0 )
	{
	  content_final6[i] = content_ab6[i]/sqrt(content_aa6[i]*content_bb6[i]);
	  econtent_final6[i] = content_final6[i] * sqrt(
							(econtent_ab6[i]/content_ab6[i])*(econtent_ab6[i]/content_ab6[i]) +
							(econtent_aa6[i]/content_aa6[i])*(econtent_aa6[i]/content_aa6[i]) +
							(econtent_bb6[i]/content_bb6[i])*(econtent_bb6[i]/content_bb6[i])
							);
	}
      cout << content_final6[i] << " " << econtent_final6[i] << endl;
      // cout << content_ab6[i] << " " << econtent_ab6[i] << endl;
      // cout << content_aa6[i] << " " << econtent_aa6[i] << endl;
      // cout << content_bb6[i] << " " << econtent_bb6[i] << endl;
    }
  // -------------------------------------------------------------------------------------------
  TGraphErrors* tge_ab6 = new TGraphErrors(12,x,content_ab6,0,econtent_ab6);
  tge_ab6->SetTitle(Form("d+Au collisions at #sqrt{s_{NN}} = %d GeV",energy));
  tge_ab6->SetMarkerStyle(kOpenSquare);
  tge_ab6->SetMarkerColor(kBlack);
  tge_ab6->Draw("ap");
  tge_ab6->SetMaximum(0.05);
  tge_ab6->SetMinimum(0.0);
  tge_ab6->GetXaxis()->SetLimits(-3,3);
  tge_ab6->GetXaxis()->SetTitle("#eta_{b}");
  tge_ab6->GetYaxis()->SetTitle("v_{a,b}(#eta_{a} = 2.25)");
  tge_ab6->GetYaxis()->SetTitleOffset(1.25);
  c1->Print("FigsDecorrelation/early_ab6.png");
  c1->Print("FigsDecorrelation/early_ab6.pdf");
  TGraphErrors* tge_aa6 = new TGraphErrors(12,x,content_aa6,0,econtent_aa6);
  tge_aa6->SetTitle(Form("d+Au collisions at #sqrt{s_{NN}} = %d GeV",energy));
  tge_aa6->SetMarkerStyle(kOpenSquare);
  tge_aa6->SetMarkerColor(kBlack);
  tge_aa6->Draw("ap");
  tge_aa6->SetMaximum(0.05);
  tge_aa6->SetMinimum(0.0);
  tge_aa6->GetXaxis()->SetLimits(-3,3);
  tge_aa6->GetXaxis()->SetTitle("#eta_{b}");
  tge_aa6->GetYaxis()->SetTitle("v_{a,a}(#eta_{a} = 2.25)");
  tge_aa6->GetYaxis()->SetTitleOffset(1.25);
  c1->Print("FigsDecorrelation/early_aa6.png");
  c1->Print("FigsDecorrelation/early_aa6.pdf");
  TGraphErrors* tge_bb6 = new TGraphErrors(12,x,content_bb6,0,econtent_bb6);
  tge_bb6->SetTitle(Form("d+Au collisions at #sqrt{s_{NN}} = %d GeV",energy));
  tge_bb6->SetMarkerStyle(kOpenSquare);
  tge_bb6->SetMarkerColor(kBlack);
  tge_bb6->Draw("ap");
  tge_bb6->SetMaximum(0.05);
  tge_bb6->SetMinimum(0.0);
  tge_bb6->GetXaxis()->SetLimits(-3,3);
  tge_bb6->GetXaxis()->SetTitle("#eta_{b}");
  tge_bb6->GetYaxis()->SetTitle("v_{b,b}(#eta_{a} = 2.25)");
  tge_bb6->GetYaxis()->SetTitleOffset(1.25);
  c1->Print("FigsDecorrelation/early_bb6.png");
  c1->Print("FigsDecorrelation/early_bb6.pdf");
  TGraphErrors* tge_final6 = new TGraphErrors(12,x,content_final6,0,econtent_final6);
  tge_final6->SetTitle(Form("d+Au collisions at #sqrt{s_{NN}} = %d GeV",energy));
  tge_final6->SetMarkerStyle(kOpenSquare);
  tge_final6->SetMarkerColor(kBlack);
  tge_final6->Draw("ap");
  tge_final6->SetMaximum(1.5);
  tge_final6->SetMinimum(0.0);
  tge_final6->GetXaxis()->SetLimits(-3,3);
  tge_final6->GetXaxis()->SetTitle("#eta_{b}");
  tge_final6->GetYaxis()->SetTitle("r(#eta_{a} = 2.25)");
  c1->Print("FigsDecorrelation/early_final6.png");
  c1->Print("FigsDecorrelation/early_final6.pdf");
  // -------------------------------------------------------------------------------------------



  // -------------------------------------------------------------------------------------------
  // --- come back here for bin 7 (last bin, north side)
  // ---
  double content_ab7[12];
  double content_aa7[12];
  double content_bb7[12];
  double content_final7[12];
  double econtent_ab7[12];
  double econtent_aa7[12];
  double econtent_bb7[12];
  double econtent_final7[12];
  double x[12];
  // -------------------------------------------------------------------------------------------
  for ( int i = 0; i < 12; ++i )
    {
      // --- set the indices for particles a and b
      int index_a = 12; // offset yby 4 because of the 4 empty bins outside the acceptance
      int index_b = i+1;
      // --- get the acceptance correction terms
      double cosb = tp1f_special_fvtx_tracks_aa_cos->GetBinContent(index_b);
      double cosa = tp1f_special_fvtx_tracks_aa_cos->GetBinContent(index_a);
      double sinb = tp1f_special_fvtx_tracks_aa_sin->GetBinContent(index_b);
      double sina = tp1f_special_fvtx_tracks_aa_sin->GetBinContent(index_a);
      // --- get eta_b
      x[i] = tp1f_special_fvtx_tracks_ab7->GetBinCenter(i+1);
      // --- get the components
      content_ab7[i] = tp1f_special_fvtx_tracks_ab7->GetBinContent(index_b);
      content_bb7[i] = tp1f_special_fvtx_tracks_aa->GetBinContent(index_b);
      content_aa7[i] = tp1f_special_fvtx_tracks_aa->GetBinContent(index_a);
      econtent_ab7[i] = tp1f_special_fvtx_tracks_ab7->GetBinError(index_b);
      econtent_bb7[i] = tp1f_special_fvtx_tracks_aa->GetBinError(index_b);
      econtent_aa7[i] = tp1f_special_fvtx_tracks_aa->GetBinError(index_a);
      // --- apply the acceptance corrections
      content_ab7[i] -= (cosa*cosb + sina*sinb);
      content_bb7[i] -= (cosb*cosb + sinb*sinb);
      content_aa7[i] -= (cosa*cosa + sina*sina);
      // --- calculate the final ratio
      content_final7[i] = -9;
      if ( content_aa7[i] > 0 && content_bb7[i] > 0 )
	{
	  content_final7[i] = content_ab7[i]/sqrt(content_aa7[i]*content_bb7[i]);
	  econtent_final7[i] = content_final7[i] * sqrt(
							(econtent_ab7[i]/content_ab7[i])*(econtent_ab7[i]/content_ab7[i]) +
							(econtent_aa7[i]/content_aa7[i])*(econtent_aa7[i]/content_aa7[i]) +
							(econtent_bb7[i]/content_bb7[i])*(econtent_bb7[i]/content_bb7[i])
							);
	}
      cout << content_final7[i] << " " << econtent_final7[i] << endl;
      // cout << content_ab7[i] << " " << econtent_ab7[i] << endl;
      // cout << content_aa7[i] << " " << econtent_aa7[i] << endl;
      // cout << content_bb7[i] << " " << econtent_bb7[i] << endl;
    }
  // -------------------------------------------------------------------------------------------
  TGraphErrors* tge_ab7 = new TGraphErrors(12,x,content_ab7,0,econtent_ab7);
  tge_ab7->SetTitle(Form("d+Au collisions at #sqrt{s_{NN}} = %d GeV",energy));
  tge_ab7->SetMarkerStyle(kOpenSquare);
  tge_ab7->SetMarkerColor(kBlack);
  tge_ab7->Draw("ap");
  tge_ab7->SetMaximum(0.05);
  tge_ab7->SetMinimum(0.0);
  tge_ab7->GetXaxis()->SetLimits(-3,3);
  tge_ab7->GetXaxis()->SetTitle("#eta_{b}");
  tge_ab7->GetYaxis()->SetTitle("v_{a,b}(#eta_{a} = 2.75)");
  tge_ab7->GetYaxis()->SetTitleOffset(1.25);
  c1->Print("FigsDecorrelation/early_ab7.png");
  c1->Print("FigsDecorrelation/early_ab7.pdf");
  TGraphErrors* tge_aa7 = new TGraphErrors(12,x,content_aa7,0,econtent_aa7);
  tge_aa7->SetTitle(Form("d+Au collisions at #sqrt{s_{NN}} = %d GeV",energy));
  tge_aa7->SetMarkerStyle(kOpenSquare);
  tge_aa7->SetMarkerColor(kBlack);
  tge_aa7->Draw("ap");
  tge_aa7->SetMaximum(0.05);
  tge_aa7->SetMinimum(0.0);
  tge_aa7->GetXaxis()->SetLimits(-3,3);
  tge_aa7->GetXaxis()->SetTitle("#eta_{b}");
  tge_aa7->GetYaxis()->SetTitle("v_{a,a}(#eta_{a} = 2.75)");
  tge_aa7->GetYaxis()->SetTitleOffset(1.25);
  c1->Print("FigsDecorrelation/early_aa7.png");
  c1->Print("FigsDecorrelation/early_aa7.pdf");
  TGraphErrors* tge_bb7 = new TGraphErrors(12,x,content_bb7,0,econtent_bb7);
  tge_bb7->SetTitle(Form("d+Au collisions at #sqrt{s_{NN}} = %d GeV",energy));
  tge_bb7->SetMarkerStyle(kOpenSquare);
  tge_bb7->SetMarkerColor(kBlack);
  tge_bb7->Draw("ap");
  tge_bb7->SetMaximum(0.05);
  tge_bb7->SetMinimum(0.0);
  tge_bb7->GetXaxis()->SetLimits(-3,3);
  tge_bb7->GetXaxis()->SetTitle("#eta_{b}");
  tge_bb7->GetYaxis()->SetTitle("v_{b,b}(#eta_{a} = 2.75)");
  tge_bb7->GetYaxis()->SetTitleOffset(1.25);
  c1->Print("FigsDecorrelation/early_bb7.png");
  c1->Print("FigsDecorrelation/early_bb7.pdf");
  TGraphErrors* tge_final7 = new TGraphErrors(12,x,content_final7,0,econtent_final7);
  tge_final7->SetTitle(Form("d+Au collisions at #sqrt{s_{NN}} = %d GeV",energy));
  tge_final7->SetMarkerStyle(kOpenSquare);
  tge_final7->SetMarkerColor(kBlack);
  tge_final7->Draw("ap");
  tge_final7->SetMaximum(1.5);
  tge_final7->SetMinimum(0.0);
  tge_final7->GetXaxis()->SetLimits(-3,3);
  tge_final7->GetXaxis()->SetTitle("#eta_{b}");
  tge_final7->GetYaxis()->SetTitle("r(#eta_{a} = 2.75)");
  c1->Print("FigsDecorrelation/early_final7.png");
  c1->Print("FigsDecorrelation/early_final7.pdf");
  // -------------------------------------------------------------------------------------------

}

