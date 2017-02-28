void epguess()
{

  TCanvas* c1 = new TCanvas();

  float bbcs_200 = 83.8;
  float bbcs_62 = 56.1;
  float bbcs_39 = 43.0;
  float bbcs_20 = 32.5;

  float fvtxs_num_200 = 242.0;
  float fvtxs_num_62 = 153.0;
  float fvtxs_num_39 = 110.0;
  float fvtxs_num_20 = 94.5;

  float fvtxs_res_200 = 0.227;
  float fvtxs_res_62 = 0.126;
  float fvtxs_res_39 = 0.064;
  float fvtxs_res_20 = 0.040;

  float fvtxs_sqt_200 = sqrt(fvtxs_num_200);
  float fvtxs_sqt_62  = sqrt(fvtxs_num_62 );
  float fvtxs_sqt_39  = sqrt(fvtxs_num_39 );
  float fvtxs_sqt_20  = sqrt(fvtxs_num_20 );

  float fvtxs_num[4] = {fvtxs_num_200,fvtxs_num_62,fvtxs_num_39,fvtxs_num_20};
  float fvtxs_sqt[4] = {fvtxs_sqt_200,fvtxs_sqt_62,fvtxs_sqt_39,fvtxs_sqt_20};
  float fvtxs_res[4] = {fvtxs_res_200,fvtxs_res_62,fvtxs_res_39,fvtxs_res_20};

  //TF1* funk = new TF1("funk","pol1",0,250);
  TF1* funk = new TF1("funk","pol6",0,250);
  funk->SetParameter(0,-0.074);
  funk->SetParameter(1,0.001257);
  funk->SetParameter(2,0.0);
  funk->SetParameter(3,0.0);
  funk->SetParameter(4,0.0);
  funk->SetParameter(5,0.0);
  funk->SetParameter(6,0.0);
  TF1* fun = new TF1("fun","pol1",0,250);

  TGraph* tg_resvssqt = new TGraph(4,fvtxs_sqt,fvtxs_res);
  tg_resvssqt->SetMarkerStyle(kFullCircle);
  tg_resvssqt->SetMarkerColor(kBlack);
  tg_resvssqt->GetXaxis()->SetTitle("#sqrt{N^{clusters}_{FVTXS}}");
  tg_resvssqt->GetYaxis()->SetTitle("Res(#Psi_{2}^{FVTXS})");
  tg_resvssqt->GetYaxis()->SetTitleOffset(1.2);
  tg_resvssqt->Draw("ap");
  tg_resvssqt->Fit(fun,"","",10,49);
  fun->Draw("same");
  TLatex tex1;
  //tex1.SetTextSize(0.05)
  tex1.DrawLatex(fvtxs_sqt_20+0.2,fvtxs_res_20-0.005,"19.6 GeV (extrapolated)");
  tex1.DrawLatex(fvtxs_sqt_39+0.2,fvtxs_res_39-0.005,"39 GeV");
  tex1.DrawLatex(fvtxs_sqt_62+0.2,fvtxs_res_62-0.005,"62.4 GeV");
  tex1.DrawLatex(fvtxs_sqt_200-1.2,fvtxs_res_200-0.005,"200 GeV");
  c1->Print("figsqt.png");
  c1->Print("figsqt.pdf");


  // cout << "guessed point is " << fun->Eval(fvtxs_sqt_20) << endl;

}
