void doit(int);

void temp_4pc()
{

  doit(200);

}

void doit(int handle)
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

  // ---
  // --- now let's have a look at the c24 vs nfvtxt
  // ---

  TProfile* tp1f_nfvtxt_fvtxns_c22 = (TProfile*)file->Get("nfvtxt_os_fvtxsfvtxn_c22");
  TProfile* tp1f_nfvtxt_fvtxns_c24a = (TProfile*)file->Get("nfvtxt_os_fvtxsfvtxn_c24a");
  TProfile* tp1f_nfvtxt_fvtxns_c24b = (TProfile*)file->Get("nfvtxt_os_fvtxsfvtxn_c24b");
  TProfile* tp1f_nfvtxt_fvtxns_c24c = (TProfile*)file->Get("nfvtxt_os_fvtxsfvtxn_c24c");
  TH1D* th1d_nfvtxt_fvtxns_c24a = tp1f_nfvtxt_fvtxns_c24a->ProjectionX();
  TH1D* th1d_nfvtxt_fvtxns_c24b = tp1f_nfvtxt_fvtxns_c24b->ProjectionX();
  TH1D* th1d_nfvtxt_fvtxns_c24c = tp1f_nfvtxt_fvtxns_c24c->ProjectionX();

  TH1D* th1d_nfvtxt_fvtxns_c22 = tp1f_nfvtxt_fvtxns_c22->ProjectionX();
  TH1D* th1d_nfvtxt_fvtxns_c22a = (TH1D*)th1d_nfvtxt_fvtxns_c22->Clone("th1d_nfvtxt_fvtxns_c22a");
  TH1D* th1d_nfvtxt_fvtxns_c22b = (TH1D*)th1d_nfvtxt_fvtxns_c22->Clone("th1d_nfvtxt_fvtxns_c22b");
  TH1D* th1d_nfvtxt_fvtxns_c22c = (TH1D*)th1d_nfvtxt_fvtxns_c22->Clone("th1d_nfvtxt_fvtxns_c22c");

  th1d_nfvtxt_fvtxns_c22a->Multiply(th1d_nfvtxt_fvtxns_c22);
  th1d_nfvtxt_fvtxns_c22a->Scale(-2.0);
  th1d_nfvtxt_fvtxns_c22a->Add(th1d_nfvtxt_fvtxns_c24a,1.0);
  th1d_nfvtxt_fvtxns_c22a->Draw();
  TLine line(0.0,0.0,75.0,0.0);
  line.SetLineWidth(2);
  line.SetLineStyle(2);
  line.Draw();
  c1->Print("FigsFour/c24a_nfvtxt.png");
  c1->Print("FigsFour/c24a_nfvtxt.pdf");

  th1d_nfvtxt_fvtxns_c22b->Multiply(th1d_nfvtxt_fvtxns_c22);
  th1d_nfvtxt_fvtxns_c22b->Scale(-2.0);
  th1d_nfvtxt_fvtxns_c22b->Add(th1d_nfvtxt_fvtxns_c24b,1.0);
  th1d_nfvtxt_fvtxns_c22b->SetMaximum(1e-5);
  th1d_nfvtxt_fvtxns_c22b->SetMinimum(-1e-5);
  th1d_nfvtxt_fvtxns_c22b->Draw();
  line.Draw();
  c1->Print("FigsFour/c24b_nfvtxt.png");
  c1->Print("FigsFour/c24b_nfvtxt.pdf");

  th1d_nfvtxt_fvtxns_c22c->Multiply(th1d_nfvtxt_fvtxns_c22);
  th1d_nfvtxt_fvtxns_c22c->Scale(-2.0);
  th1d_nfvtxt_fvtxns_c22c->Add(th1d_nfvtxt_fvtxns_c24c,1.0);
  th1d_nfvtxt_fvtxns_c22c->SetMaximum(1e-5);
  th1d_nfvtxt_fvtxns_c22c->SetMinimum(-1e-5);
  th1d_nfvtxt_fvtxns_c22c->Draw();
  line.Draw();
  c1->Print("FigsFour/c24c_nfvtxt.png");
  c1->Print("FigsFour/c24c_nfvtxt.pdf");

  const int nbins = tp1f_nfvtxt_fvtxns_c22->GetNbinsX();
  TH1D* th1d_nfvtxt_fvtxns_v24a = tp1f_nfvtxt_fvtxns_c24a->ProjectionX("th1d_nfvtxt_fvtxns_v24a");
  TH1D* th1d_nfvtxt_fvtxns_v24b = tp1f_nfvtxt_fvtxns_c24b->ProjectionX("th1d_nfvtxt_fvtxns_v24b");
  TH1D* th1d_nfvtxt_fvtxns_v24c = tp1f_nfvtxt_fvtxns_c24c->ProjectionX("th1d_nfvtxt_fvtxns_v24c");
  for ( int i = 0; i < nbins; ++i )
    {
      float q22 = th1d_nfvtxt_fvtxns_c22->GetBinContent(i+1);
      float eq22 = th1d_nfvtxt_fvtxns_c22->GetBinError(i+1);
      float eq24 = 0;
      // --- v24a
      float c24 = th1d_nfvtxt_fvtxns_c22a->GetBinContent(i+1);
      float ec24 = th1d_nfvtxt_fvtxns_c22a->GetBinError(i+1);
      float v24 = 0;
      float ev24 = 0;
      if ( c24 < 0 )
	{
	  v24 = sqrt(sqrt(-c24));
	  eq24 = th1d_nfvtxt_fvtxns_c24a->GetBinError(i+1);
	  //ev24 = fabs(ec24*v24/c24); // not right, get back to it soon
	  ev24 = (1.0/pow(fabs(-c24),0.75))*sqrt((q22*q22*eq22*eq22)+(0.0625*eq24*eq24));
	}
      th1d_nfvtxt_fvtxns_v24a->SetBinContent(i+1,v24);
      th1d_nfvtxt_fvtxns_v24a->SetBinError(i+1,ev24);
      // --- v24b
      c24 = th1d_nfvtxt_fvtxns_c22b->GetBinContent(i+1);
      ec24 = th1d_nfvtxt_fvtxns_c22b->GetBinError(i+1);
      v24 = 0;
      ev24 = 0;
      if ( c24 < 0 )
	{
	  v24 = sqrt(sqrt(-c24));
	  eq24 = th1d_nfvtxt_fvtxns_c24b->GetBinError(i+1);
	  //ev24 = fabs(ec24*v24/c24); // not right, get back to it soon
	  ev24 = (1.0/pow(fabs(-c24),0.75))*sqrt((q22*q22*eq22*eq22)+(0.0625*eq24*eq24));
	}
      th1d_nfvtxt_fvtxns_v24b->SetBinContent(i+1,v24);
      th1d_nfvtxt_fvtxns_v24b->SetBinError(i+1,ev24);
      // cout << c24 << " " << v24 << " " << ev24 << endl;
      // cout << c24 << " " << th1d_nfvtxt_fvtxns_v24b->GetBinContent(i+1) << " " << th1d_nfvtxt_fvtxns_v24b->GetBinError(i+1) << endl;
      // --- v24c
      c24 = th1d_nfvtxt_fvtxns_c22c->GetBinContent(i+1);
      ec24 = th1d_nfvtxt_fvtxns_c22c->GetBinError(i+1);
      v24 = 0;
      ev24 = 0;
      if ( c24 < 0 )
	{
	  v24 = sqrt(sqrt(-c24));
	  eq24 = th1d_nfvtxt_fvtxns_c24c->GetBinError(i+1);
	  //ev24 = fabs(ec24*v24/c24); // not right, get back to it soon
	  ev24 = (1.0/pow(fabs(-c24),0.75))*sqrt((q22*q22*eq22*eq22)+(0.0625*eq24*eq24));
	}
      th1d_nfvtxt_fvtxns_v24c->SetBinContent(i+1,v24);
      th1d_nfvtxt_fvtxns_v24c->SetBinError(i+1,ev24);
      // cout << c24 << " " << v24 << " " << ev24 << endl;
    }

  //TProfile* tp1f_22 = (TProfile*)file_200->Get(Form("nfvtxt_os_%s_c%d2",name,harmonic));
  TProfile* tp1f_22 = (TProfile*)file->Get(Form("nfvtxt_os_fvtxsfvtxn_c22"));
  tp1f_22->SetName("tp1f_22");
  TH1D* th1d_vn_22 = sqrt(tp1f_22);
  th1d_vn_22->SetLineColor(kRed);

  th1d_nfvtxt_fvtxns_v24a->Draw();
  th1d_nfvtxt_fvtxns_v24a->GetXaxis()->SetTitle("N_{tracks}^{FVTX}");
  th1d_nfvtxt_fvtxns_v24a->GetYaxis()->SetTitle("v_{2}{}");
  th1d_nfvtxt_fvtxns_v24a->SetMaximum(0.1);
  th1d_nfvtxt_fvtxns_v24a->SetMinimum(-0.01);
  TLegend leg(0.68,0.68,0.88,0.88);
  leg.AddEntry(th1d_vn_22,"v_{2}{2}","el");
  leg.AddEntry(th1d_nfvtxt_fvtxns_v24a,"v_{2}{4}","el");
  leg.SetTextSize(0.045);
  leg.Draw();
  line.Draw();
  c1->Print("FigsFour/v24a_nfvtxt.png");
  c1->Print("FigsFour/v24a_nfvtxt.pdf");
  th1d_vn_22->Draw("same");
  c1->Print("FigsFour/v24a_and_v22_nfvtxt.png");
  c1->Print("FigsFour/v24a_and_v22_nfvtxt.pdf");
  th1d_nfvtxt_fvtxns_v24b->Draw();
  th1d_nfvtxt_fvtxns_v24b->GetXaxis()->SetTitle("N_{tracks}^{FVTX}");
  th1d_nfvtxt_fvtxns_v24b->GetYaxis()->SetTitle("v_{2}{}");
  th1d_nfvtxt_fvtxns_v24b->SetMaximum(0.1);
  th1d_nfvtxt_fvtxns_v24b->SetMinimum(-0.01);
  leg.Draw();
  line.Draw();
  c1->Print("FigsFour/v24b_nfvtxt.png");
  c1->Print("FigsFour/v24b_nfvtxt.pdf");
  th1d_vn_22->Draw("same");
  c1->Print("FigsFour/v24b_and_v22_nfvtxt.png");
  c1->Print("FigsFour/v24b_and_v22_nfvtxt.pdf");
  th1d_nfvtxt_fvtxns_v24c->Draw();
  th1d_nfvtxt_fvtxns_v24c->GetXaxis()->SetTitle("N_{tracks}^{FVTX}");
  th1d_nfvtxt_fvtxns_v24c->GetYaxis()->SetTitle("v_{2}{}");
  th1d_nfvtxt_fvtxns_v24c->SetMaximum(0.1);
  th1d_nfvtxt_fvtxns_v24c->SetMinimum(-0.01);
  leg.Draw();
  line.Draw();
  c1->Print("FigsFour/v24c_nfvtxt.png");
  c1->Print("FigsFour/v24c_nfvtxt.pdf");
  th1d_vn_22->Draw("same");
  c1->Print("FigsFour/v24c_and_v22_nfvtxt.png");
  c1->Print("FigsFour/v24c_and_v22_nfvtxt.pdf");

  // ---

  TProfile* tp1f_nfvtxt_fvtxs_c22 = (TProfile*)file->Get("nfvtxt_os_fvtxs_c22");
  TProfile* tp1f_nfvtxt_fvtxs_c24a = (TProfile*)file->Get("nfvtxt_os_fvtxs_c24");
  TH1D* th1d_nfvtxt_fvtxs_c24a = tp1f_nfvtxt_fvtxs_c24a->ProjectionX();
  TH1D* th1d_nfvtxt_fvtxs_c22 = tp1f_nfvtxt_fvtxs_c22->ProjectionX();
  TH1D* th1d_nfvtxt_fvtxs_c22a = (TH1D*)th1d_nfvtxt_fvtxs_c22->Clone("th1d_nfvtxt_fvtxs_c22a");
  th1d_nfvtxt_fvtxs_c22a->Multiply(th1d_nfvtxt_fvtxs_c22);
  th1d_nfvtxt_fvtxs_c22a->Scale(-2.0);
  th1d_nfvtxt_fvtxs_c22a->Add(th1d_nfvtxt_fvtxs_c24a,1.0);
  th1d_nfvtxt_fvtxs_c22a->Draw();
  th1d_nfvtxt_fvtxs_c22a->SetMaximum(1e-4);
  th1d_nfvtxt_fvtxs_c22a->SetMinimum(-1e-4);
  TLine line(0.0,0.0,75.0,0.0);
  line.SetLineWidth(2);
  line.SetLineStyle(2);
  line.Draw();
  c1->Print("FigsFour/c24_fvtxs_nfvtxt.png");
  c1->Print("FigsFour/c24_fvtxs_nfvtxt.pdf");

  TProfile* tp1f_nfvtxt_fvtxn_c22 = (TProfile*)file->Get("nfvtxt_os_fvtxn_c22");
  TProfile* tp1f_nfvtxt_fvtxn_c24a = (TProfile*)file->Get("nfvtxt_os_fvtxn_c24");
  TH1D* th1d_nfvtxt_fvtxn_c24a = tp1f_nfvtxt_fvtxn_c24a->ProjectionX();
  TH1D* th1d_nfvtxt_fvtxn_c22 = tp1f_nfvtxt_fvtxn_c22->ProjectionX();
  TH1D* th1d_nfvtxt_fvtxn_c22a = (TH1D*)th1d_nfvtxt_fvtxn_c22->Clone("th1d_nfvtxt_fvtxn_c22a");
  th1d_nfvtxt_fvtxn_c22a->Multiply(th1d_nfvtxt_fvtxn_c22);
  th1d_nfvtxt_fvtxn_c22a->Scale(-2.0);
  th1d_nfvtxt_fvtxn_c22a->Add(th1d_nfvtxt_fvtxn_c24a,1.0);
  th1d_nfvtxt_fvtxn_c22a->Draw();
  th1d_nfvtxt_fvtxn_c22a->SetMaximum(1e-4);
  th1d_nfvtxt_fvtxn_c22a->SetMinimum(-1e-4);
  TLine line(0.0,0.0,75.0,0.0);
  line.SetLineWidth(2);
  line.SetLineStyle(2);
  line.Draw();
  c1->Print("FigsFour/c24_fvtxn_nfvtxt.png");
  c1->Print("FigsFour/c24_fvtxn_nfvtxt.pdf");

  TProfile* tp1f_nfvtxt_fvtxs_tracks_c22 = (TProfile*)file->Get("nfvtxt_os_fvtxs_tracks_c22");
  TProfile* tp1f_nfvtxt_fvtxs_tracks_c24a = (TProfile*)file->Get("nfvtxt_os_fvtxs_tracks_c24");
  TH1D* th1d_nfvtxt_fvtxs_tracks_c24a = tp1f_nfvtxt_fvtxs_tracks_c24a->ProjectionX();
  TH1D* th1d_nfvtxt_fvtxs_tracks_c22 = tp1f_nfvtxt_fvtxs_tracks_c22->ProjectionX();
  TH1D* th1d_nfvtxt_fvtxs_tracks_c22a = (TH1D*)th1d_nfvtxt_fvtxs_tracks_c22->Clone("th1d_nfvtxt_fvtxs_tracks_c22a");
  th1d_nfvtxt_fvtxs_tracks_c22a->Multiply(th1d_nfvtxt_fvtxs_tracks_c22);
  th1d_nfvtxt_fvtxs_tracks_c22a->Scale(-2.0);
  th1d_nfvtxt_fvtxs_tracks_c22a->Add(th1d_nfvtxt_fvtxs_tracks_c24a,1.0);
  th1d_nfvtxt_fvtxs_tracks_c22a->Draw();
  th1d_nfvtxt_fvtxs_tracks_c22a->SetMaximum(1e-4);
  th1d_nfvtxt_fvtxs_tracks_c22a->SetMinimum(-1e-4);
  TLine line(0.0,0.0,75.0,0.0);
  line.SetLineWidth(2);
  line.SetLineStyle(2);
  line.Draw();
  c1->Print("FigsFour/c24_fvtxs_tracks_nfvtxt.png");
  c1->Print("FigsFour/c24_fvtxs_tracks_nfvtxt.pdf");

  TProfile* tp1f_nfvtxt_fvtxn_tracks_c22 = (TProfile*)file->Get("nfvtxt_os_fvtxn_tracks_c22");
  TProfile* tp1f_nfvtxt_fvtxn_tracks_c24a = (TProfile*)file->Get("nfvtxt_os_fvtxn_tracks_c24");
  TH1D* th1d_nfvtxt_fvtxn_tracks_c24a = tp1f_nfvtxt_fvtxn_tracks_c24a->ProjectionX();
  TH1D* th1d_nfvtxt_fvtxn_tracks_c22 = tp1f_nfvtxt_fvtxn_tracks_c22->ProjectionX();
  TH1D* th1d_nfvtxt_fvtxn_tracks_c22a = (TH1D*)th1d_nfvtxt_fvtxn_tracks_c22->Clone("th1d_nfvtxt_fvtxn_tracks_c22a");
  th1d_nfvtxt_fvtxn_tracks_c22a->Multiply(th1d_nfvtxt_fvtxn_tracks_c22);
  th1d_nfvtxt_fvtxn_tracks_c22a->Scale(-2.0);
  th1d_nfvtxt_fvtxn_tracks_c22a->Add(th1d_nfvtxt_fvtxn_tracks_c24a,1.0);
  th1d_nfvtxt_fvtxn_tracks_c22a->Draw();
  th1d_nfvtxt_fvtxn_tracks_c22a->SetMaximum(1e-4);
  th1d_nfvtxt_fvtxn_tracks_c22a->SetMinimum(-1e-4);
  TLine line(0.0,0.0,75.0,0.0);
  line.SetLineWidth(2);
  line.SetLineStyle(2);
  line.Draw();
  c1->Print("FigsFour/c24_fvtxn_tracks_nfvtxt.png");
  c1->Print("FigsFour/c24_fvtxn_tracks_nfvtxt.pdf");

  TProfile* tp1f_nfvtxt_fvtxc_tracks_c22 = (TProfile*)file->Get("nfvtxt_os_fvtxc_tracks_c22");
  TProfile* tp1f_nfvtxt_fvtxc_tracks_c24a = (TProfile*)file->Get("nfvtxt_os_fvtxc_tracks_c24");
  TH1D* th1d_nfvtxt_fvtxc_tracks_c24a = tp1f_nfvtxt_fvtxc_tracks_c24a->ProjectionX();
  TH1D* th1d_nfvtxt_fvtxc_tracks_c22 = tp1f_nfvtxt_fvtxc_tracks_c22->ProjectionX();
  TH1D* th1d_nfvtxt_fvtxc_tracks_c22a = (TH1D*)th1d_nfvtxt_fvtxc_tracks_c22->Clone("th1d_nfvtxt_fvtxc_tracks_c22a");
  th1d_nfvtxt_fvtxc_tracks_c22a->Multiply(th1d_nfvtxt_fvtxc_tracks_c22);
  th1d_nfvtxt_fvtxc_tracks_c22a->Scale(-2.0);
  th1d_nfvtxt_fvtxc_tracks_c22a->Add(th1d_nfvtxt_fvtxc_tracks_c24a,1.0);
  th1d_nfvtxt_fvtxc_tracks_c22a->Draw();
  th1d_nfvtxt_fvtxc_tracks_c22a->SetMaximum(1e-4);
  th1d_nfvtxt_fvtxc_tracks_c22a->SetMinimum(-1e-4);
  TLine line(0.0,0.0,75.0,0.0);
  line.SetLineWidth(2);
  line.SetLineStyle(2);
  line.Draw();
  c1->Print("FigsFour/c24_fvtxc_tracks_nfvtxt.png");
  c1->Print("FigsFour/c24_fvtxc_tracks_nfvtxt.pdf");



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

