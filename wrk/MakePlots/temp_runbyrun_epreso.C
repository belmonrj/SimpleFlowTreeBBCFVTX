void dostuff(int, double&, double&);

void temp_runbyrun_epreso()
{

  gStyle->SetOptTitle(0);

  TCanvas* c1 = new TCanvas("c1","");

  ifstream fin("list.short");

  double index[67];
  double bbc_reso[67];
  double fvtx_reso[67];
  int run;
  for ( int i = 0; i < 67; ++i )
    {
      fin>>run;
      double bbc_temp = 0;
      double fvtx_temp = 0;
      dostuff(run,bbc_temp,fvtx_temp);
      // if ( bbc_temp != bbc_temp ) bbc_temp = 0;
      // if ( fvtx_temp != fvtx_temp ) fvtx_temp = 0;
      index[i] = i+0.5;
      bbc_reso[i] = bbc_temp;
      fvtx_reso[i] = fvtx_temp;
    }

  TGraph* tg_bbc = new TGraph(67,index,bbc_reso);
  tg_bbc->SetMarkerStyle(kFullCircle);
  tg_bbc->SetMarkerColor(kRed);
  tg_bbc->Draw("ap");
  tg_bbc->SetTitle("");
  tg_bbc->GetXaxis()->SetLimits(-2,69);
  tg_bbc->GetXaxis()->SetTitle("Run Index");
  tg_bbc->GetYaxis()->SetTitle("Resolution");
  tg_bbc->SetMinimum(0);
  tg_bbc->SetMaximum(0.5);

  TGraph* tg_fvtx = new TGraph(67,index,fvtx_reso);
  tg_fvtx->SetMarkerStyle(kFullCircle);
  tg_fvtx->SetMarkerColor(kBlue);
  tg_fvtx->Draw("p");

  TLegend* leg = new TLegend(0.18,0.68,0.28,0.88);
  leg->SetHeader("Event plane resolution run by run");
  leg->AddEntry(tg_bbc,"BBC","p");
  leg->AddEntry(tg_fvtx,"FVTX","p");
  leg->SetTextSize(0.05);
  leg->Draw();

  // tg_bbc->Fit("pol0","","",0,67);
  // tg_fvtx->Fit("pol0","","",0,67);
  TF1* fun_bbc = new TF1("fun_bbc","pol0",0,67);
  TF1* fun_fvtx = new TF1("fun_fvtx","pol0",0,67);
  fun_bbc->SetLineColor(kBlack);
  fun_fvtx->SetLineColor(kBlack);
  // tg_bbc->Fit(fun_bbc,"R");
  // tg_fvtx->Fit(fun_fvtx,"R");
  fun_bbc->SetParameter(0,0.0969358);
  fun_fvtx->SetParameter(0,0.224709);
  fun_bbc->Draw("same");
  fun_fvtx->Draw("same");

  c1->Print("runbyrun_epreso.png");
  c1->Print("runbyrun_epreso.pdf");

}

void dostuff(int run, double& bbc, double& fvtx)
{

  int verbosity  = 0;

  gStyle->SetOptTitle(1);

  //TCanvas* c1 = new TCanvas("c1","");

  TFile* file = TFile::Open(Form("input/hist_%d.root",run));

  // ---

  TProfile* tp1f_bbc_fvtxN = (TProfile*)file->Get("tp1f_reso2_BBC_FVTXN");
  TProfile* tp1f_bbc_fvtx = (TProfile*)file->Get("tp1f_reso2_BBC_FVTX");
  TProfile* tp1f_fvtxN_fvtx = (TProfile*)file->Get("tp1f_reso2_FVTXS_FVTXN");

  float float_bbc_fvtxN = tp1f_bbc_fvtxN->GetBinContent(1);
  float float_bbc_fvtx = tp1f_bbc_fvtx->GetBinContent(1);
  float float_fvtxN_fvtx = tp1f_fvtxN_fvtx->GetBinContent(1);

  if ( verbosity > 1 ) cout <<  float_bbc_fvtxN << endl;
  if ( verbosity > 1 ) cout <<  float_bbc_fvtx << endl;
  if ( verbosity > 1 ) cout <<  float_fvtxN_fvtx << endl;

  float reso_bbc = sqrt((float_bbc_fvtxN*float_bbc_fvtx)/float_fvtxN_fvtx);
  float reso_fvtx = sqrt((float_fvtxN_fvtx*float_bbc_fvtx)/float_bbc_fvtxN);

  if ( verbosity > 1 ) cout << "bbc resolution is " << reso_bbc << endl;
  if ( verbosity > 1 ) cout << "fvtx resolution is " << reso_fvtx << endl;

  // ---

  TProfile* tp1f_bbc_cnt = (TProfile*)file->Get("tp1f_reso2_BBC_CNT");
  tp1f_bbc_fvtx = (TProfile*)file->Get("tp1f_reso2_BBC_FVTX");
  TProfile* tp1f_cnt_fvtx = (TProfile*)file->Get("tp1f_reso2_CNT_FVTX");

  float float_bbc_cnt = tp1f_bbc_cnt->GetBinContent(1);
  float_bbc_fvtx = tp1f_bbc_fvtx->GetBinContent(1);
  float float_cnt_fvtx = tp1f_cnt_fvtx->GetBinContent(1);

  if ( verbosity > 1 ) cout <<  float_bbc_cnt << endl;
  if ( verbosity > 1 ) cout <<  float_bbc_fvtx << endl;
  if ( verbosity > 1 ) cout <<  float_cnt_fvtx << endl;

  reso_bbc = sqrt((float_bbc_cnt*float_bbc_fvtx)/float_cnt_fvtx);
  reso_fvtx = sqrt((float_cnt_fvtx*float_bbc_fvtx)/float_bbc_cnt);

  if ( verbosity > 0 ) cout << "bbc resolution is " << reso_bbc << endl;
  if ( verbosity > 0 ) cout << "fvtx resolution is " << reso_fvtx << endl;

  bbc = reso_bbc;
  fvtx = reso_fvtx;

  file->Close();

}

