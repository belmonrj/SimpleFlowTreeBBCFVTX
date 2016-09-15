void temp_stats()
{
  ts_energy(200);
  ts_energy(62);
  ts_energy(39);
  ts_energy(20);
}



void ts_energy(int energy)
{

  TCanvas* c1 = new TCanvas();
  c1->SetLogy();

  gStyle->SetOptTitle(1);

  TFile* file = TFile::Open(Form("input/combined_%d.root",energy));
  TProfile* tp1f_fvtxsv2 = (TProfile*)file->Get("fvtxs_v2_both_docalib");
  tp1f_fvtxsv2->SetName("tp1f_fvtxsv2");
  TH1D* th1d_fvtxsv2 = tp1f_fvtxsv2->ProjectionX("th1d_fvtxsv2");

  TH1D* th1d_numtracks = (TH1D*)th1d_fvtxsv2->Clone("th1d_numtracks");
  TH1D* th1d_counterr = (TH1D*)th1d_fvtxsv2->Clone("th1d_counterr");
  TH1D* th1d_projv2err = (TH1D*)th1d_fvtxsv2->Clone("th1d_projv2err");
  TH1D* th1d_v2err = (TH1D*)th1d_fvtxsv2->Clone("th1d_v2err");

  const int nbins = tp1f_fvtxsv2->GetNbinsX();

  for ( int i = 0; i < nbins; ++i )
    {
      double numtracks = tp1f_fvtxsv2->GetBinEntries(i+1);
      double counterr = 0;
      if ( numtracks > 0 ) counterr = 1.0/sqrt(numtracks);
      double v2 = tp1f_fvtxsv2->GetBinContent(i+1);
      double v2err = tp1f_fvtxsv2->GetBinError(i+1);
      double projv2err = v2*counterr;
      th1d_numtracks->SetBinContent(i+1,numtracks);
      th1d_numtracks->SetBinError(i+1,numtracks*counterr);
      th1d_counterr->SetBinContent(i+1,counterr);
      th1d_counterr->SetBinError(i+1,0);
      th1d_v2err->SetBinContent(i+1,v2err);
      th1d_v2err->SetBinError(i+1,0);
      th1d_projv2err->SetBinContent(i+1,projv2err);
      th1d_projv2err->SetBinError(i+1,0);
    }

  th1d_numtracks->Draw();
  c1->Print(Form("fig_numtracks_%d.png",energy));
  th1d_counterr->Draw();
  c1->Print(Form("fig_counterr_%d.png",energy));
  th1d_v2err->Draw();
  c1->Print(Form("fig_v2err_%d.png",energy));
  th1d_projv2err->Draw();
  c1->Print(Form("fig_projv2err_%d.png",energy));

  th1d_counterr->SetLineColor(kBlue);
  th1d_counterr->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  th1d_counterr->GetYaxis()->SetTitle("Statistical Uncertainty");
  th1d_counterr->SetTitle(Form("d+Au collisions at #sqrt{s_{NN}} = %d GeV",energy));
  th1d_counterr->Draw();
  th1d_v2err->SetLineColor(kBlack);
  th1d_v2err->Draw("same");
  TLegend* leg = new TLegend(0.18,0.68,0.38,0.88);
  leg->AddEntry(th1d_counterr,"1/#sqrt{N_{tracks}}","l");
  leg->AddEntry(th1d_v2err,"From v_{2} histogram","l");
  leg->SetTextSize(0.045);
  leg->Draw();
  c1->Print(Form("fig_countandv2err_%d.png",energy));

  c1->SetLogy(0);
  th1d_v2err->Divide(th1d_counterr);
  th1d_v2err->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  th1d_v2err->GetYaxis()->SetTitle("Ratio of Statistical Uncertainties");
  th1d_v2err->SetTitle(Form("d+Au collisions at #sqrt{s_{NN}} = %d GeV",energy));
  th1d_v2err->Draw();
  c1->Print(Form("fig_countandv2err_ratio_%d.png",energy));


}
