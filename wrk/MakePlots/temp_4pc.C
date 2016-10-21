void dohandle(int);

void docalc(TProfile*, TProfile*, const char*, const char*, const char*, int);

void temp_4pc()
{

  dohandle(200);
  // dohandle(62);
  // dohandle(39);
  // dohandle(20);

}

void dohandle(int handle)
{

  gStyle->SetOptTitle(0);

  TCanvas* c1 = new TCanvas("c1","");

  TFile* file = NULL;
  if ( handle <= 200 ) file = TFile::Open(Form("input/combined_%d.root",handle));
  else if ( handle > 454000 ) file = TFile::Open(Form("input/hist_%d.root",handle));
  else
    {
      cout << "YOU'RE GONNA DIE" << endl;
      return;
    }

  TProfile* two;
  TProfile* four;

  // --- type A clusters (v^4 but no autocorrelation corrections)
  two = (TProfile*)file->Get("nfvtxc_os_fvtxsfvtxn_c22");
  four = (TProfile*)file->Get("nfvtxc_os_fvtxsfvtxn_c24a");
  docalc(two,four,"clusters","mixed","2subA",handle);

  // --- type D clusters (v^4 with autocorrelation corrections, but clusters cause problems for that)
  two = (TProfile*)file->Get("nfvtxc_os_fvtxsfvtxn_c22");
  four = (TProfile*)file->Get("nfvtxc_os_fvtxsfvtxn_c24d");
  docalc(two,four,"clusters","mixed","2subD",handle);

  // ---

  // --- type A clusters with 4 subevents (++--)
  two = (TProfile*)file->Get("nfvtxt_os_fvtxsfvtxn_c22");
  four = (TProfile*)file->Get("nfvtxt_os_fvtxsfvtxn_ce01_c24a");
  docalc(two,four,"clusters","mixed","4subA",handle);

  // --- type D clusters with 4 subevents (--++)
  two = (TProfile*)file->Get("nfvtxt_os_fvtxsfvtxn_c22");
  four = (TProfile*)file->Get("nfvtxt_os_fvtxsfvtxn_ce01_c24d");
  docalc(two,four,"clusters","mixed","4subD",handle);

  // ---

  // --- type A tracks (v^4 but no autocorrelation corrections)
  two = (TProfile*)file->Get("nfvtxt_os_fvtxsfvtxn_tracks_c22");
  four = (TProfile*)file->Get("nfvtxt_os_fvtxsfvtxn_tracks_c24a");
  docalc(two,four,"tracks","mixed","2subA",handle);

  // --- type D tracks (v^4 with autocorrelation corrections, should be okay for tracks???)
  two = (TProfile*)file->Get("nfvtxt_os_fvtxsfvtxn_tracks_c22");
  four = (TProfile*)file->Get("nfvtxt_os_fvtxsfvtxn_tracks_c24d");
  docalc(two,four,"tracks","mixed","2subD",handle);

  // ---
  // ---
  // ---

  // --- south tracks only
  two = (TProfile*)file->Get("nfvtxt_os_fvtxs_tracks_c22");
  four = (TProfile*)file->Get("nfvtxt_os_fvtxs_tracks_c24");
  docalc(two,four,"tracks","south","1sub",handle);

  // --- north tracks only
  two = (TProfile*)file->Get("nfvtxt_os_fvtxn_tracks_c22");
  four = (TProfile*)file->Get("nfvtxt_os_fvtxn_tracks_c24");
  docalc(two,four,"tracks","north","1sub",handle);

  // --- combined tracks only
  two = (TProfile*)file->Get("nfvtxt_os_fvtxc_tracks_c22");
  four = (TProfile*)file->Get("nfvtxt_os_fvtxc_tracks_c24");
  docalc(two,four,"tracks","combined","1sub",handle);

  // --- south clusters only
  two = (TProfile*)file->Get("nfvtxc_os_fvtxs_c22");
  four = (TProfile*)file->Get("nfvtxc_os_fvtxs_c24");
  docalc(two,four,"clusters","south","1sub",handle);

  // --- north clusters only
  two = (TProfile*)file->Get("nfvtxc_os_fvtxn_c22");
  four = (TProfile*)file->Get("nfvtxc_os_fvtxn_c24");
  docalc(two,four,"clusters","north","1sub",handle);

  // --- combined clusters only
  two = (TProfile*)file->Get("nfvtxc_os_fvtxc_c22");
  four = (TProfile*)file->Get("nfvtxc_os_fvtxc_c24");
  docalc(two,four,"clusters","combined","1sub",handle);


}

void docalc(TProfile* tp1f_c22, TProfile* tp1f_c24a, const char* type, const char* side, const char* detail, int handle)
{

  char description[100];
  sprintf(description,"%s_%s_%s_%d",type,side,detail,handle);

  cout << "description is " << description << endl;

  TH1D* th1d_c24a = tp1f_c24a->ProjectionX(Form("th1d_c24a_%s",description));
  TH1D* th1d_c22 = tp1f_c22->ProjectionX(Form("th1d_c22_%s",description));
  TH1D* th1d_c22a = (TH1D*)th1d_c22->Clone(Form("th1d_c22a_%s",description));

  th1d_c22a->Multiply(th1d_c22);
  th1d_c22a->Scale(2.0);
  th1d_c22a->SetMaximum(2e-05);
  th1d_c22a->SetMinimum(-4e-06);
  th1d_c22a->Draw();
  th1d_c22a->GetXaxis()->SetTitle(Form("N_{%s}^{%s}",type,side));
  th1d_c24a->SetLineColor(kRed);
  th1d_c24a->Draw("same");
  TLegend leg1(0.68,0.68,0.88,0.88);
  leg1.AddEntry(th1d_c22a,"2#LT#LT2#GT#GT^{2}","el");
  leg1.AddEntry(th1d_c24a,"#LT#LT4#GT#GT","el");
  leg1.SetTextSize(0.045);
  leg1.Draw();
  TLine line(0.0,0.0,75.0,0.0);
  line.SetLineWidth(2);
  line.SetLineStyle(2);
  line.Draw();
  c1->Print(Form("FigsFour/components_222and4_nfvtxt_%s.png",description));
  c1->Print(Form("FigsFour/components_222and4_nfvtxt_%s.pdf",description));
  th1d_c22a->SetMaximum(1e-03);
  th1d_c22a->SetMinimum(-1e-03);
  c1->Print(Form("FigsFour/components_222and4_nfvtxt_%s_zout.png",description));
  c1->Print(Form("FigsFour/components_222and4_nfvtxt_%s_zout.pdf",description));

  th1d_c22a->Scale(-1.0);
  th1d_c22a->Add(th1d_c24a,1.0);
  th1d_c22a->SetMaximum(2e-05);
  th1d_c22a->SetMinimum(-2e-05);
  th1d_c22a->Draw();
  line.Draw();
  c1->Print(Form("FigsFour/c24_nfvtxt_%s.png",description));
  c1->Print(Form("FigsFour/c24_nfvtxt_%s.pdf",description));
  th1d_c22a->SetMaximum(1e-03);
  th1d_c22a->SetMinimum(-1e-03);
  c1->Print(Form("FigsFour/c24_nfvtxt_%s_zout.png",description));
  c1->Print(Form("FigsFour/c24_nfvtxt_%s_zout.pdf",description));

  const int nbins = tp1f_c22->GetNbinsX();
  TH1D* th1d_v24a = tp1f_c24a->ProjectionX(Form("th1d_v24a_%s",description));
  for ( int i = 0; i < nbins; ++i )
    {
      float q22 = th1d_c22->GetBinContent(i+1);
      float eq22 = th1d_c22->GetBinError(i+1);
      float eq24 = 0;
      // --- v24a
      float c24 = th1d_c22a->GetBinContent(i+1);
      float ec24 = th1d_c22a->GetBinError(i+1);
      float v24 = 0;
      float ev24 = 0;
      if ( c24 < 0 )
	{
	  v24 = sqrt(sqrt(-c24));
	  eq24 = th1d_c24a->GetBinError(i+1);
	  ev24 = (1.0/pow(fabs(-c24),0.75))*sqrt((q22*q22*eq22*eq22)+(0.0625*eq24*eq24));
	}
      th1d_v24a->SetBinContent(i+1,v24);
      th1d_v24a->SetBinError(i+1,ev24);
    }

  TH1D* th1d_vn_22 = sqrt(th1d_c22);
  th1d_vn_22->SetLineColor(kRed);

  th1d_v24a->Draw();
  th1d_v24a->GetXaxis()->SetTitle(Form("N_{%s}^{%s}",type,side));
  th1d_v24a->GetYaxis()->SetTitle("v_{2}{}");
  th1d_v24a->SetMaximum(0.2);
  th1d_v24a->SetMinimum(-0.02);
  TLegend leg(0.68,0.68,0.88,0.88);
  leg.AddEntry(th1d_vn_22,"v_{2}{2}","el");
  leg.AddEntry(th1d_v24a,"v_{2}{4}","el");
  leg.SetTextSize(0.045);
  leg.Draw();
  line.Draw();
  c1->Print(Form("FigsFour/v24_nfvtxt_%s.png",description));
  c1->Print(Form("FigsFour/v24_nfvtxt_%s.pdf",description));
  th1d_vn_22->Draw("same");
  c1->Print(Form("FigsFour/v24_and_v22_nfvtxt_%s.png",description));
  c1->Print(Form("FigsFour/v24_and_v22_nfvtxt_%s.pdf",description));

  delete c1;

}

TH1D* sqrt(TProfile* tp1f)
{
  TH1D* th1d = tp1f->ProjectionX(Form("%s_px",tp1f->GetName()));
  return sqrt(th1d);
}

TH1D* sqrt(TH1D* th1d)
{
  TH1D* hnew = (TH1D*)th1d->Clone(Form("%s_sqrt",th1d->GetName()));
  int nbins = th1d->GetNbinsX();
  for ( int i = 0; i < nbins; ++i )
    {
      double content = th1d->GetBinContent(i+1);
      double uncert = th1d->GetBinError(i+1);
      if ( content > 0 ) content = sqrt(content);
      else content = 0;
      if ( uncert > 0 && content > 0 ) uncert = uncert/content;
      else uncert = 0;
      hnew->SetBinContent(i+1,content);
      hnew->SetBinError(i+1,uncert);
    }
  return hnew;
}

