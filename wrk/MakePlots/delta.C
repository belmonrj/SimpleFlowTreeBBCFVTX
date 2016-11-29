const float pi = 3.1415926535;

TH1D* pleasefixme3h(TH1D*);


void delta()
{

  TFile* file = TFile::Open("input/combined_200.root");


  TH1D* th1dA = (TH1D*)file->Get("fvtxs_v3_both_phipsi");
  TH1D* th1d0 = (TH1D*)file->Get("fvtxs0_v3_both_phipsi");
  TH1D* th1d1 = (TH1D*)file->Get("fvtxs1_v3_both_phipsi");
  TH1D* th1d2 = (TH1D*)file->Get("fvtxs2_v3_both_phipsi");
  TH1D* th1d3 = (TH1D*)file->Get("fvtxs3_v3_both_phipsi");

  TH1D* th1dA_fixed = pleasefixme3h(th1dA);
  th1dA_fixed->Draw();
  double max = th1dA_fixed->GetBinContent(th1dA_fixed->GetMaximumBin());
  th1dA_fixed->SetMaximum(max*1.002);
  th1dA_fixed->SetMinimum(max*0.995);
  th1dA_fixed->GetXaxis()->SetRangeUser(-1.2,1.2);
  th1dA_fixed->GetXaxis()->SetTitle("#phi-#Psi_{3}");
  TF1* fun = new TF1("fun","[0]+[1]*TMath::Cos(3*x-[2])",-2,2);
  TLatex tex(0.5,0.8,Form("offset %.2e #pm %.2e",fun->GetParameter(2),fun->GetParError(2)));
  tex.SetNDC();
  fun->SetParLimits(2,-0.9,0.9);
  th1dA_fixed->Fit(fun,"","",-1.1,1.1);
  tex.DrawLatex(0.45,0.82,Form("offset %.2e #pm %.2e",fun->GetParameter(2),fun->GetParError(2)));
  c1->Print("FigsEventPlane/phipsi3_th1dA.png");
  c1->Print("FigsEventPlane/phipsi3_th1dA.pdf");

  TH1D* th1d0_fixed = pleasefixme3h(th1d0);
  th1d0_fixed->Draw();
  max = th1d0_fixed->GetBinContent(th1d0_fixed->GetMaximumBin());
  th1d0_fixed->SetMaximum(max*1.002);
  th1d0_fixed->SetMinimum(max*0.995);
  th1d0_fixed->GetXaxis()->SetRangeUser(-1.2,1.2);
  th1d0_fixed->GetXaxis()->SetTitle("#phi-#Psi_{3}");
  th1d0_fixed->Fit(fun,"","",-1.1,1.1);
  tex.DrawLatex(0.45,0.82,Form("offset %.2e #pm %.2e",fun->GetParameter(2),fun->GetParError(2)));
  c1->Print("FigsEventPlane/phipsi3_th1d0.png");
  c1->Print("FigsEventPlane/phipsi3_th1d0.pdf");

  TH1D* th1d1_fixed = pleasefixme3h(th1d1);
  th1d1_fixed->Draw();
  max = th1d1_fixed->GetBinContent(th1d1_fixed->GetMaximumBin());
  th1d1_fixed->SetMaximum(max*1.002);
  th1d1_fixed->SetMinimum(max*0.995);
  th1d1_fixed->GetXaxis()->SetRangeUser(-1.2,1.2);
  th1d1_fixed->GetXaxis()->SetTitle("#phi-#Psi_{3}");
  th1d1_fixed->Fit(fun,"","",-1.1,1.1);
  tex.DrawLatex(0.45,0.82,Form("offset %.2e #pm %.2e",fun->GetParameter(2),fun->GetParError(2)));
  c1->Print("FigsEventPlane/phipsi3_th1d1.png");
  c1->Print("FigsEventPlane/phipsi3_th1d1.pdf");

  TH1D* th1d2_fixed = pleasefixme3h(th1d2);
  th1d2_fixed->Draw();
  max = th1d2_fixed->GetBinContent(th1d2_fixed->GetMaximumBin());
  th1d2_fixed->SetMaximum(max*1.002);
  th1d2_fixed->SetMinimum(max*0.995);
  th1d2_fixed->GetXaxis()->SetRangeUser(-1.2,1.2);
  th1d2_fixed->GetXaxis()->SetTitle("#phi-#Psi_{3}");
  th1d2_fixed->Fit(fun,"","",-1.1,1.1);
  tex.DrawLatex(0.45,0.82,Form("offset %.2e #pm %.2e",fun->GetParameter(2),fun->GetParError(2)));
  c1->Print("FigsEventPlane/phipsi3_th1d2.png");
  c1->Print("FigsEventPlane/phipsi3_th1d2.pdf");

  TH1D* th1d3_fixed = pleasefixme3h(th1d3);
  th1d3_fixed->Draw();
  max = th1d3_fixed->GetBinContent(th1d3_fixed->GetMaximumBin());
  th1d3_fixed->SetMaximum(max*1.002);
  th1d3_fixed->SetMinimum(max*0.995);
  th1d3_fixed->GetXaxis()->SetRangeUser(-1.2,1.2);
  th1d3_fixed->GetXaxis()->SetTitle("#phi-#Psi_{3}");
  th1d3_fixed->Fit(fun,"","",-1.1,1.1);
  tex.DrawLatex(0.45,0.82,Form("offset %.2e #pm %.2e",fun->GetParameter(2),fun->GetParError(2)));
  c1->Print("FigsEventPlane/phipsi3_th1d3.png");
  c1->Print("FigsEventPlane/phipsi3_th1d3.pdf");

  // ---

  TH1D* th1d_B = (TH1D*)file->Get("th1d_os_dreso2_BBC_CNT");
  TH1D* th1d_S = (TH1D*)file->Get("th1d_os_dreso2_CNT_FVTX");

  TF1* fun2 = new TF1("fun2","[0]+[1]*TMath::Cos(2*x-[2])",-2,2);

  TH1D* th1d_B_fixed = pleasefixme2h(th1d_B);
  TH1D* th1d_S_fixed = pleasefixme2h(th1d_S);

  th1d_B_fixed->Draw();
  max = th1d_B_fixed->GetBinContent(th1d_B_fixed->GetMaximumBin());
  th1d_B_fixed->SetMaximum(max*1.002);
  th1d_B_fixed->SetMinimum(max*0.95);
  th1d_B_fixed->GetXaxis()->SetRangeUser(-1.7,1.7);
  th1d_B_fixed->GetXaxis()->SetTitle("#phi-#Psi_{2}");
  th1d_B_fixed->Fit(fun2,"","",-1.6,1.6);
  tex.DrawLatex(0.45,0.82,Form("offset %.2e #pm %.2e",fun->GetParameter(2),fun->GetParError(2)));
  c1->Print("FigsEventPlane/phipsi2_th1d_B.png");
  c1->Print("FigsEventPlane/phipsi2_th1d_B.pdf");

  th1d_S_fixed->Draw();
  max = th1d_S_fixed->GetBinContent(th1d_S_fixed->GetMaximumBin());
  th1d_S_fixed->SetMaximum(max*1.002);
  th1d_S_fixed->SetMinimum(max*0.92);
  th1d_S_fixed->GetXaxis()->SetRangeUser(-1.7,1.7);
  th1d_S_fixed->GetXaxis()->SetTitle("#phi-#Psi_{2}");
  th1d_S_fixed->Fit(fun2,"","",-1.6,1.6);
  tex.DrawLatex(0.45,0.82,Form("offset %.2e #pm %.2e",fun->GetParameter(2),fun->GetParError(2)));
  c1->Print("FigsEventPlane/phipsi2_th1d_S.png");
  c1->Print("FigsEventPlane/phipsi2_th1d_S.pdf");

  // ----------------------------------

  TH1D* th1d_2_A = (TH1D*)file->Get("fvtxs_v2_both_phipsi");
  TH1D* th1d_2_0 = (TH1D*)file->Get("fvtxs0_v2_both_phipsi");
  TH1D* th1d_2_1 = (TH1D*)file->Get("fvtxs1_v2_both_phipsi");
  TH1D* th1d_2_2 = (TH1D*)file->Get("fvtxs2_v2_both_phipsi");
  TH1D* th1d_2_3 = (TH1D*)file->Get("fvtxs3_v2_both_phipsi");

  TH1D* th1d_2_A_fixed = pleasefixme2h(th1d_2_A);
  th1d_2_A_fixed->Draw();
  max = th1d_2_A_fixed->GetBinContent(th1d_2_A_fixed->GetMaximumBin());
  th1d_2_A_fixed->SetMaximum(max*1.002);
  th1d_2_A_fixed->SetMinimum(max*0.92);
  th1d_2_A_fixed->GetXaxis()->SetRangeUser(-1.7,1.7);
  th1d_2_A_fixed->GetXaxis()->SetTitle("#phi-#Psi_{3}");
  delete fun;
  fun = new TF1("fun","[0]+[1]*TMath::Cos(2*x-[2])",-2,2);
  fun->SetParLimits(2,-0.9,0.9);
  th1d_2_A_fixed->Fit(fun,"","",-1.6,1.6);
  tex.DrawLatex(0.45,0.82,Form("offset %.2e #pm %.2e",fun->GetParameter(2),fun->GetParError(2)));
  c1->Print("FigsEventPlane/phipsi2_th1dA.png");
  c1->Print("FigsEventPlane/phipsi2_th1dA.pdf");

  TH1D* th1d_2_0_fixed = pleasefixme2h(th1d_2_0);
  th1d_2_0_fixed->Draw();
  max = th1d_2_0_fixed->GetBinContent(th1d_2_0_fixed->GetMaximumBin());
  th1d_2_0_fixed->SetMaximum(max*1.002);
  th1d_2_0_fixed->SetMinimum(max*0.92);
  th1d_2_0_fixed->GetXaxis()->SetRangeUser(-1.7,1.7);
  th1d_2_0_fixed->GetXaxis()->SetTitle("#phi-#Psi_{3}");
  th1d_2_0_fixed->Fit(fun,"","",-1.6,1.6);
  tex.DrawLatex(0.45,0.82,Form("offset %.2e #pm %.2e",fun->GetParameter(2),fun->GetParError(2)));
  c1->Print("FigsEventPlane/phipsi2_th1d0.png");
  c1->Print("FigsEventPlane/phipsi2_th1d0.pdf");

  TH1D* th1d_2_1_fixed = pleasefixme2h(th1d_2_1);
  th1d_2_1_fixed->Draw();
  max = th1d_2_1_fixed->GetBinContent(th1d_2_1_fixed->GetMaximumBin());
  th1d_2_1_fixed->SetMaximum(max*1.002);
  th1d_2_1_fixed->SetMinimum(max*0.92);
  th1d_2_1_fixed->GetXaxis()->SetRangeUser(-1.7,1.7);
  th1d_2_1_fixed->GetXaxis()->SetTitle("#phi-#Psi_{3}");
  th1d_2_1_fixed->Fit(fun,"","",-1.6,1.6);
  tex.DrawLatex(0.45,0.82,Form("offset %.2e #pm %.2e",fun->GetParameter(2),fun->GetParError(2)));
  c1->Print("FigsEventPlane/phipsi2_th1d1.png");
  c1->Print("FigsEventPlane/phipsi2_th1d1.pdf");

  TH1D* th1d_2_2_fixed = pleasefixme2h(th1d_2_2);
  th1d_2_2_fixed->Draw();
  max = th1d_2_2_fixed->GetBinContent(th1d_2_2_fixed->GetMaximumBin());
  th1d_2_2_fixed->SetMaximum(max*1.002);
  th1d_2_2_fixed->SetMinimum(max*0.92);
  th1d_2_2_fixed->GetXaxis()->SetRangeUser(-1.7,1.7);
  th1d_2_2_fixed->GetXaxis()->SetTitle("#phi-#Psi_{3}");
  th1d_2_2_fixed->Fit(fun,"","",-1.6,1.6);
  tex.DrawLatex(0.45,0.82,Form("offset %.2e #pm %.2e",fun->GetParameter(2),fun->GetParError(2)));
  c1->Print("FigsEventPlane/phipsi2_th1d2.png");
  c1->Print("FigsEventPlane/phipsi2_th1d2.pdf");

  TH1D* th1d_2_3_fixed = pleasefixme2h(th1d_2_3);
  th1d_2_3_fixed->Draw();
  max = th1d_2_3_fixed->GetBinContent(th1d_2_3_fixed->GetMaximumBin());
  th1d_2_3_fixed->SetMaximum(max*1.002);
  th1d_2_3_fixed->SetMinimum(max*0.92);
  th1d_2_3_fixed->GetXaxis()->SetRangeUser(-1.7,1.7);
  th1d_2_3_fixed->GetXaxis()->SetTitle("#phi-#Psi_{3}");
  th1d_2_3_fixed->Fit(fun,"","",-1.6,1.6);
  tex.DrawLatex(0.45,0.82,Form("offset %.2e #pm %.2e",fun->GetParameter(2),fun->GetParError(2)));
  c1->Print("FigsEventPlane/phipsi2_th1d3.png");
  c1->Print("FigsEventPlane/phipsi2_th1d3.pdf");

  // --- now event more fun stuff

  TH1D* th1d012 = (TH1D*)file->Get("fvtxs012_v3_both_phipsi");
  TH1D* th1d013 = (TH1D*)file->Get("fvtxs013_v3_both_phipsi");
  TH1D* th1d023 = (TH1D*)file->Get("fvtxs023_v3_both_phipsi");
  TH1D* th1d123 = (TH1D*)file->Get("fvtxs123_v3_both_phipsi");

  TH1D* th1d012_fixed = pleasefixme3h(th1d012);
  th1d012_fixed->Draw();
  max = th1d012_fixed->GetBinContent(th1d012_fixed->GetMaximumBin());
  th1d012_fixed->SetMaximum(max*1.002);
  th1d012_fixed->SetMinimum(max*0.995);
  th1d012_fixed->GetXaxis()->SetRangeUser(-1.2,1.2);
  th1d012_fixed->GetXaxis()->SetTitle("#phi-#Psi_{3}");
  th1d012_fixed->Fit(fun,"","",-1.1,1.1);
  tex.DrawLatex(0.45,0.82,Form("offset %.2e #pm %.2e",fun->GetParameter(2),fun->GetParError(2)));
  c1->Print("FigsEventPlane/phipsi3_th1d012.png");
  c1->Print("FigsEventPlane/phipsi3_th1d012.pdf");

  TH1D* th1d013_fixed = pleasefixme3h(th1d013);
  th1d013_fixed->Draw();
  max = th1d013_fixed->GetBinContent(th1d013_fixed->GetMaximumBin());
  th1d013_fixed->SetMaximum(max*1.002);
  th1d013_fixed->SetMinimum(max*0.995);
  th1d013_fixed->GetXaxis()->SetRangeUser(-1.2,1.2);
  th1d013_fixed->GetXaxis()->SetTitle("#phi-#Psi_{3}");
  th1d013_fixed->Fit(fun,"","",-1.1,1.1);
  tex.DrawLatex(0.45,0.82,Form("offset %.2e #pm %.2e",fun->GetParameter(2),fun->GetParError(2)));
  c1->Print("FigsEventPlane/phipsi3_th1d013.png");
  c1->Print("FigsEventPlane/phipsi3_th1d013.pdf");

  TH1D* th1d023_fixed = pleasefixme3h(th1d023);
  th1d023_fixed->Draw();
  max = th1d023_fixed->GetBinContent(th1d023_fixed->GetMaximumBin());
  th1d023_fixed->SetMaximum(max*1.002);
  th1d023_fixed->SetMinimum(max*0.995);
  th1d023_fixed->GetXaxis()->SetRangeUser(-1.2,1.2);
  th1d023_fixed->GetXaxis()->SetTitle("#phi-#Psi_{3}");
  th1d023_fixed->Fit(fun,"","",-1.1,1.1);
  tex.DrawLatex(0.45,0.82,Form("offset %.2e #pm %.2e",fun->GetParameter(2),fun->GetParError(2)));
  c1->Print("FigsEventPlane/phipsi3_th1d023.png");
  c1->Print("FigsEventPlane/phipsi3_th1d023.pdf");

  TH1D* th1d123_fixed = pleasefixme3h(th1d123);
  th1d123_fixed->Draw();
  max = th1d123_fixed->GetBinContent(th1d123_fixed->GetMaximumBin());
  th1d123_fixed->SetMaximum(max*1.002);
  th1d123_fixed->SetMinimum(max*0.995);
  th1d123_fixed->GetXaxis()->SetRangeUser(-1.2,1.2);
  th1d123_fixed->GetXaxis()->SetTitle("#phi-#Psi_{3}");
  th1d123_fixed->Fit(fun,"","",-1.1,1.1);
  tex.DrawLatex(0.45,0.82,Form("offset %.2e #pm %.2e",fun->GetParameter(2),fun->GetParError(2)));
  c1->Print("FigsEventPlane/phipsi3_th1d123.png");
  c1->Print("FigsEventPlane/phipsi3_th1d123.pdf");

  // --- these are all messed up for some reason, maybe come back to this some other time

  th1d_B = (TH1D*)file->Get("th1d_os_dreso3_BBC_CNT");
  th1d_S = (TH1D*)file->Get("th1d_os_dreso3_CNT_FVTX");
  th1d_B_fixed = pleasefixme3h(th1d_B);
  th1d_S_fixed = pleasefixme3h(th1d_S);

  th1d_B_fixed->Draw();
  max = th1d_B_fixed->GetBinContent(th1d_B_fixed->GetMaximumBin());
  th1d_B_fixed->SetMaximum(max*1.002);
  th1d_B_fixed->SetMinimum(max*0.995);
  th1d_B_fixed->GetXaxis()->SetRangeUser(-1.2,1.2);
  th1d_B_fixed->GetXaxis()->SetTitle("#phi-#Psi_{3}");
  th1d_B_fixed->Fit(fun,"","",-1.1,1.1);
  tex.DrawLatex(0.45,0.82,Form("offset %.2e #pm %.2e",fun->GetParameter(2),fun->GetParError(2)));
  c1->Print("FigsEventPlane/phipsi3_th1d_B.png");
  c1->Print("FigsEventPlane/phipsi3_th1d_B.pdf");

  th1d_S_fixed->Draw();
  max = th1d_S_fixed->GetBinContent(th1d_S_fixed->GetMaximumBin());
  th1d_S_fixed->SetMaximum(max*1.002);
  th1d_S_fixed->SetMinimum(max*0.995);
  th1d_S_fixed->GetXaxis()->SetRangeUser(-1.2,1.2);
  th1d_S_fixed->GetXaxis()->SetTitle("#phi-#Psi_{3}");
  th1d_S_fixed->Fit(fun,"","",-1.1,1.1);
  tex.DrawLatex(0.45,0.82,Form("offset %.2e #pm %.2e",fun->GetParameter(2),fun->GetParError(2)));
  c1->Print("FigsEventPlane/phipsi3_th1d_S.png");
  c1->Print("FigsEventPlane/phipsi3_th1d_S.pdf");



}

TH1D* pleasefixme3h(TH1D* in)
{
  TH1D* clone = (TH1D*)in->Clone();
  int nbins = clone->GetNbinsX();
  for ( int i = 0; i < nbins; ++i )
    {
      float phi = clone->GetBinCenter(i+1);
      float amount = clone->GetBinContent(i+1);
      if ( phi < -pi )
	{
	  clone->SetBinContent(i+1,0);
	  clone->Fill(phi+2*pi,amount);
	}
      if ( phi > pi )
	{
	  clone->SetBinContent(i+1,0);
	  clone->Fill(phi-2*pi,amount);
	}
    }
  TH1D* clone2 = (TH1D*)clone->Clone();
  for ( int i = 0; i < nbins; ++i )
    {
      float phi = clone2->GetBinCenter(i+1);
      float amount = clone2->GetBinContent(i+1);
      if ( phi < -pi/3 )
	{
	  clone2->SetBinContent(i+1,0);
	  clone2->Fill(phi+(2*pi/3),amount);
	}
      if ( phi > pi/3 )
	{
	  clone2->SetBinContent(i+1,0);
	  clone2->Fill(phi-(2*pi/3),amount);
	}
    }
  return clone2;
}


TH1D* pleasefixme2h(TH1D* in)
{
  TH1D* clone = (TH1D*)in->Clone();
  int nbins = clone->GetNbinsX();
  for ( int i = 0; i < nbins; ++i )
    {
      float phi = clone->GetBinCenter(i+1);
      float amount = clone->GetBinContent(i+1);
      if ( phi < -pi )
	{
	  clone->SetBinContent(i+1,0);
	  clone->Fill(phi+2*pi,amount);
	}
      if ( phi > pi )
	{
	  clone->SetBinContent(i+1,0);
	  clone->Fill(phi-2*pi,amount);
	}
    }
  TH1D* clone2 = (TH1D*)clone->Clone();
  for ( int i = 0; i < nbins; ++i )
    {
      float phi = clone2->GetBinCenter(i+1);
      float amount = clone2->GetBinContent(i+1);
      if ( phi < -pi/2 )
	{
	  clone2->SetBinContent(i+1,0);
	  clone2->Fill(phi+pi,amount);
	}
      if ( phi > pi/2 )
	{
	  clone2->SetBinContent(i+1,0);
	  clone2->Fill(phi-pi,amount);
	}
    }
  return clone2;
}
