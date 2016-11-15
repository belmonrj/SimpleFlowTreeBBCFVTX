const float pi = 3.1415926535;

TH1D* pleasefixme(TH1D*);


void delta()
{

  TFile* file = TFile::Open("input/combined_200.root");
  cout << file << endl;

  TH1D* th1dA = (TH1D*)file->Get("fvtxs_v3_both_phipsi");
  cout << th1dA << endl;
  TH1D* th1d0 = (TH1D*)file->Get("fvtxs0_v3_both_phipsi");
  TH1D* th1d1 = (TH1D*)file->Get("fvtxs1_v3_both_phipsi");
  TH1D* th1d2 = (TH1D*)file->Get("fvtxs2_v3_both_phipsi");
  TH1D* th1d3 = (TH1D*)file->Get("fvtxs3_v3_both_phipsi");

  TH1D* th1dA_fixed = pleasefixme(th1dA);
  th1dA_fixed->Draw();
  double max = th1dA_fixed->GetBinContent(th1dA_fixed->GetMaximumBin());
  th1dA_fixed->SetMaximum(max*1.002);
  th1dA_fixed->SetMinimum(max*0.995);
  th1dA_fixed->GetXaxis()->SetRangeUser(-1.2,1.2);
  th1dA_fixed->GetXaxis()->SetTitle("#phi-#Psi_{3}");
  TF1* fun = new TF1("fun","[0]+[1]*TMath::Cos(3*x-[2])",-1.2,1.2);
  TLatex tex(0.5,0.8,Form("offset %.2e",fun->GetParameter(2)));
  tex.SetNDC();
  fun->SetParLimits(2,-0.9,0.9);
  th1dA_fixed->Fit(fun,"R");
  tex.DrawLatex(0.65,0.82,Form("offset %.2e",fun->GetParameter(2)));
  c1->Print("phipsi3_th1dA.png");
  c1->Print("phipsi3_th1dA.pdf");

  TH1D* th1d0_fixed = pleasefixme(th1d0);
  th1d0_fixed->Draw();
  max = th1d0_fixed->GetBinContent(th1d0_fixed->GetMaximumBin());
  th1d0_fixed->SetMaximum(max*1.002);
  th1d0_fixed->SetMinimum(max*0.995);
  th1d0_fixed->GetXaxis()->SetRangeUser(-1.2,1.2);
  th1d0_fixed->GetXaxis()->SetTitle("#phi-#Psi_{3}");
  th1d0_fixed->Fit(fun,"R");
  tex.DrawLatex(0.65,0.82,Form("offset %.2e",fun->GetParameter(2)));
  c1->Print("phipsi3_th1d0.png");
  c1->Print("phipsi3_th1d0.pdf");

  TH1D* th1d1_fixed = pleasefixme(th1d1);
  th1d1_fixed->Draw();
  max = th1d1_fixed->GetBinContent(th1d1_fixed->GetMaximumBin());
  th1d1_fixed->SetMaximum(max*1.002);
  th1d1_fixed->SetMinimum(max*0.995);
  th1d1_fixed->GetXaxis()->SetRangeUser(-1.2,1.2);
  th1d1_fixed->GetXaxis()->SetTitle("#phi-#Psi_{3}");
  th1d1_fixed->Fit(fun,"R");
  tex.DrawLatex(0.65,0.82,Form("offset %.2e",fun->GetParameter(2)));
  c1->Print("phipsi3_th1d1.png");
  c1->Print("phipsi3_th1d1.pdf");

  TH1D* th1d2_fixed = pleasefixme(th1d2);
  th1d2_fixed->Draw();
  max = th1d2_fixed->GetBinContent(th1d2_fixed->GetMaximumBin());
  th1d2_fixed->SetMaximum(max*1.002);
  th1d2_fixed->SetMinimum(max*0.995);
  th1d2_fixed->GetXaxis()->SetRangeUser(-1.2,1.2);
  th1d2_fixed->GetXaxis()->SetTitle("#phi-#Psi_{3}");
  th1d2_fixed->Fit(fun,"R");
  tex.DrawLatex(0.65,0.82,Form("offset %.2e",fun->GetParameter(2)));
  c1->Print("phipsi3_th1d2.png");
  c1->Print("phipsi3_th1d2.pdf");

  TH1D* th1d3_fixed = pleasefixme(th1d3);
  th1d3_fixed->Draw();
  max = th1d3_fixed->GetBinContent(th1d3_fixed->GetMaximumBin());
  th1d3_fixed->SetMaximum(max*1.002);
  th1d3_fixed->SetMinimum(max*0.995);
  th1d3_fixed->GetXaxis()->SetRangeUser(-1.2,1.2);
  th1d3_fixed->GetXaxis()->SetTitle("#phi-#Psi_{3}");
  th1d3_fixed->Fit(fun,"R");
  tex.DrawLatex(0.65,0.82,Form("offset %.2e",fun->GetParameter(2)));
  c1->Print("phipsi3_th1d3.png");
  c1->Print("phipsi3_th1d3.pdf");

}

TH1D* pleasefixme(TH1D* in)
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
