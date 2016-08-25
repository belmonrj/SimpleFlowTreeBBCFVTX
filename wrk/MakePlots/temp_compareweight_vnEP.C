void temp_doit();


void temp_compareweight_vnEP()
{

  temp_doit(200,2);
  temp_doit(200,3);

}


void temp_doit(int energy, int harmonic)
{

  gStyle->SetOptTitle(1);

  TCanvas* c1 = new TCanvas("c1","");

  TFile* file = TFile::Open(Form("input/combined_%d.root",energy));

  TProfile* hvn_fvtxs = (TProfile*)file->Get(Form("fvtxs_v%d_both_docalib",harmonic));
  ofstream fout((const char*)Form("DataTextFiles/run16dau%d_v%dfvtxs.dat",energy,harmonic));
  for ( int i = 0; i < hvn_fvtxs->GetNbinsX(); ++i ) fout << hvn_fvtxs->GetBinCenter(i+1) << " " << hvn_fvtxs->GetBinContent(i+1) << " " << hvn_fvtxs->GetBinError(i+1) << " " << endl;
  fout.close();
  hvn_fvtxs->Draw();
  // the 62 GeV is actually 62.4 and the 20 GeV is actually 19.6, so need to modify
  hvn_fvtxs->SetTitle(Form("d+Au collisions at #sqrt{s_{NN}} = %d GeV",energy));
  if ( energy == 62 ) hvn_fvtxs->SetTitle(Form("d+Au collisions at #sqrt{s_{NN}} = 62.4 GeV"));
  if ( energy == 20 ) hvn_fvtxs->SetTitle(Form("d+Au collisions at #sqrt{s_{NN}} = 19.6 GeV"));
  hvn_fvtxs->SetMaximum(0.04);
  hvn_fvtxs->SetMinimum(0.0);
  if ( harmonic == 3 )
  {
  hvn_fvtxs->SetMaximum(0.004);
  hvn_fvtxs->SetMinimum(-0.002);
  }
  hvn_fvtxs->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hvn_fvtxs->GetYaxis()->SetTitle(Form("v_{%d}{EP}",harmonic));
  hvn_fvtxs->GetYaxis()->SetTitleOffset(1.25);

  TProfile* hvn_fvtxsNW = (TProfile*)file->Get(Form("fvtxs_v%d_both_dcnw",harmonic));
  hvn_fvtxsNW->SetLineColor(kRed);
  hvn_fvtxsNW->Draw("same");

  TLegend* leg = new TLegend(0.18,0.68,0.38,0.88);
  leg->AddEntry(hvn_fvtxsNW,"Run16 FVTXS no weight","el");
  leg->AddEntry(hvn_fvtxs,"Run16 FVTXS","el");
  leg->SetTextSize(0.05);
  leg->Draw();

  c1->Print(Form("FigsHarmonicCoefficient/weightcompare_run16dau%d_uncorrv%d_fvtxs.pdf",energy,harmonic));
  c1->Print(Form("FigsHarmonicCoefficient/weightcompare_run16dau%d_uncorrv%d_fvtxs.png",energy,harmonic));





  TProfile* hvn_fvtxn = (TProfile*)file->Get(Form("fvtxn_v%d_both_docalib",harmonic));
  ofstream fout((const char*)Form("DataTextFiles/run16dau%d_v%dfvtxn.dat",energy,harmonic));
  for ( int i = 0; i < hvn_fvtxn->GetNbinsX(); ++i ) fout << hvn_fvtxn->GetBinCenter(i+1) << " " << hvn_fvtxn->GetBinContent(i+1) << " " << hvn_fvtxn->GetBinError(i+1) << " " << endl;
  fout.close();
  hvn_fvtxn->Draw();
  // the 62 GeV is actually 62.4 and the 20 GeV is actually 19.6, so need to modify
  hvn_fvtxn->SetTitle(Form("d+Au collisions at #sqrt{s_{NN}} = %d GeV",energy));
  if ( energy == 62 ) hvn_fvtxn->SetTitle(Form("d+Au collisions at #sqrt{s_{NN}} = 62.4 GeV"));
  if ( energy == 20 ) hvn_fvtxn->SetTitle(Form("d+Au collisions at #sqrt{s_{NN}} = 19.6 GeV"));
  hvn_fvtxn->SetMaximum(0.04);
  hvn_fvtxn->SetMinimum(0.0);
  if ( harmonic == 3 )
  {
  hvn_fvtxn->SetMaximum(0.004);
  hvn_fvtxn->SetMinimum(-0.002);
  }
  hvn_fvtxn->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hvn_fvtxn->GetYaxis()->SetTitle(Form("v_{%d}{EP}",harmonic));
  hvn_fvtxn->GetYaxis()->SetTitleOffset(1.25);

  TProfile* hvn_fvtxnNW = (TProfile*)file->Get(Form("fvtxn_v%d_both_dcnw",harmonic));
  hvn_fvtxnNW->SetLineColor(kRed);
  hvn_fvtxnNW->Draw("same");

  TLegend* leg = new TLegend(0.18,0.68,0.38,0.88);
  leg->AddEntry(hvn_fvtxnNW,"Run16 FVTXN no weight","el");
  leg->AddEntry(hvn_fvtxn,"Run16 FVTXN","el");
  leg->SetTextSize(0.05);
  leg->Draw();

  c1->Print(Form("FigsHarmonicCoefficient/weightcompare_run16dau%d_uncorrv%d_fvtxn.pdf",energy,harmonic));
  c1->Print(Form("FigsHarmonicCoefficient/weightcompare_run16dau%d_uncorrv%d_fvtxn.png",energy,harmonic));





  TProfile* hvn_bbcs = (TProfile*)file->Get(Form("bbcs_v%d_both_docalib",harmonic));
  ofstream fout((const char*)Form("DataTextFiles/run16dau%d_v%dbbcs.dat",energy,harmonic));
  for ( int i = 0; i < hvn_bbcs->GetNbinsX(); ++i ) fout << hvn_bbcs->GetBinCenter(i+1) << " " << hvn_bbcs->GetBinContent(i+1) << " " << hvn_bbcs->GetBinError(i+1) << " " << endl;
  fout.close();
  hvn_bbcs->Draw();
  // the 62 GeV is actually 62.4 and the 20 GeV is actually 19.6, so need to modify
  hvn_bbcs->SetTitle(Form("d+Au collisions at #sqrt{s_{NN}} = %d GeV",energy));
  if ( energy == 62 ) hvn_bbcs->SetTitle(Form("d+Au collisions at #sqrt{s_{NN}} = 62.4 GeV"));
  if ( energy == 20 ) hvn_bbcs->SetTitle(Form("d+Au collisions at #sqrt{s_{NN}} = 19.6 GeV"));
  hvn_bbcs->SetMaximum(0.02);
  hvn_bbcs->SetMinimum(0.0);
  if ( harmonic == 3 )
  {
  hvn_bbcs->SetMaximum(0.004);
  hvn_bbcs->SetMinimum(-0.002);
  }
  hvn_bbcs->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hvn_bbcs->GetYaxis()->SetTitle(Form("v_{%d}{EP}",harmonic));
  hvn_bbcs->GetYaxis()->SetTitleOffset(1.25);

  TProfile* hvn_bbcsNW = (TProfile*)file->Get(Form("bbcs_v%d_both_dcnw",harmonic));
  hvn_bbcsNW->SetLineColor(kRed);
  hvn_bbcsNW->Draw("same");

  TLegend* leg = new TLegend(0.18,0.68,0.38,0.88);
  leg->AddEntry(hvn_bbcsNW,"Run16 BBCS no weight","el");
  leg->AddEntry(hvn_bbcs,"Run16 BBCS","el");
  leg->SetTextSize(0.05);
  leg->Draw();

  c1->Print(Form("FigsHarmonicCoefficient/weightcompare_run16dau%d_uncorrv%d_bbcs.pdf",energy,harmonic));
  c1->Print(Form("FigsHarmonicCoefficient/weightcompare_run16dau%d_uncorrv%d_bbcs.png",energy,harmonic));




  // --- need to look at event plane resolutions...

  delete c1;


}
