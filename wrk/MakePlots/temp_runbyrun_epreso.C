void dostuff(int, int, double&, double&, double&, double&);

void makeplots(int, int);



void temp_runbyrun_epreso()
{

  makeplots(200,2);
  makeplots(200,3);
  makeplots(200,42);

  makeplots(62,2);
  makeplots(62,3);
  makeplots(62,42);

  makeplots(20,2);
  makeplots(20,3);
  makeplots(20,42);

  makeplots(39,2);
  makeplots(39,3);
  makeplots(39,42);

}

void makeplots(int energy, int harmonic)
{

  gStyle->SetOptTitle(0);

  TCanvas* c1 = new TCanvas("c1","");

  ifstream fin((const char*)Form("list_%d.short",energy));

  double index[100];
  double bbc_reso[100];
  double fvtx_reso[100];
  double bbc_reso_fn[100];
  double fvtx_reso_fn[100];
  int run;
  int counter = 0;
  for ( int i = 0; i < 100; ++i )
    {
      if ( fin.eof() ) break;
      fin >> run;
      //cout << run << endl;
      ++counter;
      double bbc_temp = 0;
      double fvtx_temp = 0;
      double bbc_temp_fn = 0;
      double fvtx_temp_fn = 0;
      dostuff(run,harmonic,bbc_temp,fvtx_temp,bbc_temp_fn,fvtx_temp_fn);
      // if ( bbc_temp != bbc_temp ) bbc_temp = 0;
      // if ( fvtx_temp != fvtx_temp ) fvtx_temp = 0;
      index[i] = i+0.5;
      bbc_reso[i] = bbc_temp;
      fvtx_reso[i] = fvtx_temp;
      bbc_reso_fn[i] = bbc_temp_fn;
      fvtx_reso_fn[i] = fvtx_temp_fn;
    }

  TGraph* tg_bbc = new TGraph(counter,index,bbc_reso);
  tg_bbc->SetMarkerStyle(kFullCircle);
  tg_bbc->SetMarkerColor(kRed);
  tg_bbc->Draw("ap");
  tg_bbc->SetTitle("");
  tg_bbc->GetXaxis()->SetLimits(-2,counter+2);
  tg_bbc->GetXaxis()->SetTitle("Run Index");
  tg_bbc->GetYaxis()->SetTitle("Resolution");
  tg_bbc->SetMinimum(0);
  tg_bbc->SetMaximum(0.5);

  TGraph* tg_fvtx = new TGraph(counter,index,fvtx_reso);
  tg_fvtx->SetMarkerStyle(kFullCircle);
  tg_fvtx->SetMarkerColor(kBlue);
  tg_fvtx->Draw("p");

  TLegend* leg = new TLegend(0.18,0.68,0.28,0.88);
  leg->SetHeader("Event plane resolution run by run");
  leg->AddEntry(tg_bbc,"BBC","p");
  leg->AddEntry(tg_fvtx,"FVTX","p");
  leg->SetTextSize(0.05);
  leg->Draw();

  TF1* fun_bbc = new TF1("fun_bbc","pol0",0,counter);
  TF1* fun_fvtx = new TF1("fun_fvtx","pol0",0,counter);
  fun_bbc->SetLineColor(kBlack);
  fun_fvtx->SetLineColor(kBlack);
  if ( energy == 200 && harmonic == 2 )
    {
      fun_bbc->SetParameter(0,0.0969358);
      fun_fvtx->SetParameter(0,0.224709);
      fun_bbc->Draw("same");
      fun_fvtx->Draw("same");
    }
  // tg_bbc->Fit(fun_bbc,"R");
  // tg_fvtx->Fit(fun_fvtx,"R");

  c1->Print(Form("FigsEventPlane/runbyrun_%d_epreso_%d.png",energy,harmonic));
  c1->Print(Form("FigsEventPlane/runbyrun_%d_epreso_%d.pdf",energy,harmonic));

  TGraph* tg_bbc_fn = new TGraph(67,index,bbc_reso_fn);
  tg_bbc_fn->SetMarkerStyle(kOpenCircle);
  tg_bbc_fn->SetMarkerColor(kRed);
  tg_bbc_fn->Draw("p");

  TGraph* tg_fvtx_fn = new TGraph(67,index,fvtx_reso_fn);
  tg_fvtx_fn->SetMarkerStyle(kOpenCircle);
  tg_fvtx_fn->SetMarkerColor(kBlue);
  tg_fvtx_fn->Draw("p");

  c1->Print(Form("FigsEventPlane/runbyrun_%d_epreso_%d_fn.png",energy,harmonic));
  c1->Print(Form("FigsEventPlane/runbyrun_%d_epreso_%d_fn.pdf",energy,harmonic));

  delete c1;

}



void dostuff(int run, int hh, double& bbc, double& fvtx, double& bbc_fn, double& fvtx_fn)
{

  if ( run <= 0 )
    {
      cout << "FATAL: bad run number" << endl;
      exit(-1);
    }

  int verbosity  = 0;

  gStyle->SetOptTitle(1);

  //TCanvas* c1 = new TCanvas("c1","");

  TFile* file = TFile::Open(Form("input/hist_%d.root",run));

  // ---

  TProfile* tp1f_bbc_fvtxN = (TProfile*)file->Get(Form("tp1f_reso%d_BBC_FVTXN",hh));
  TProfile* tp1f_bbc_fvtx = (TProfile*)file->Get(Form("tp1f_reso%d_BBC_FVTX",hh));
  TProfile* tp1f_fvtxN_fvtx = (TProfile*)file->Get(Form("tp1f_reso%d_FVTXS_FVTXN",hh));

  float float_bbc_fvtxN = tp1f_bbc_fvtxN->GetBinContent(1);
  float float_bbc_fvtx = tp1f_bbc_fvtx->GetBinContent(1);
  float float_fvtxN_fvtx = tp1f_fvtxN_fvtx->GetBinContent(1);

  if ( verbosity > 1 ) cout <<  float_bbc_fvtxN << endl;
  if ( verbosity > 1 ) cout <<  float_bbc_fvtx << endl;
  if ( verbosity > 1 ) cout <<  float_fvtxN_fvtx << endl;

  float reso_bbc = sqrt((float_bbc_fvtxN*float_bbc_fvtx)/float_fvtxN_fvtx);
  float reso_fvtx = sqrt((float_fvtxN_fvtx*float_bbc_fvtx)/float_bbc_fvtxN);

  bbc_fn = reso_bbc;
  fvtx_fn = reso_fvtx;

  if ( verbosity > 1 ) cout << "bbc resolution is " << reso_bbc << endl;
  if ( verbosity > 1 ) cout << "fvtx resolution is " << reso_fvtx << endl;

  // ---

  TProfile*  tp1f_bbc_cnt = (TProfile*)file->Get(Form("tp1f_reso2_BBC_CNT",hh));
  tp1f_bbc_fvtx = (TProfile*)file->Get(Form("tp1f_reso2_BBC_FVTX",hh));
  TProfile* tp1f_cnt_fvtx = (TProfile*)file->Get(Form("tp1f_reso2_CNT_FVTX",hh));

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

