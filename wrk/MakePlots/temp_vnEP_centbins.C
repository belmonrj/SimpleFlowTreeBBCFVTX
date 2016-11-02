void doenergy(int, int);
void diagnostic(int, int);



void temp_vnEP_centbins()
{

  doenergy(62,2,0);

}

void doenergy(int energy, int harmonic, int centbin)
{

  gStyle->SetOptTitle(1);

  TCanvas* c1 = new TCanvas("c1","");

  //TFile* file = TFile::Open(Form("input/combined_%d.root",energy));
  TFile* file = TFile::Open(Form("cent_input/hist_455792_centbins.root"));

  // ---

  TProfile* tp1f_bbc_fvtx = (TProfile*)file->Get(Form("tp1f_reso%d_BBC_FVTX_cent%d",harmonic,centbin));
  TProfile* tp1f_bbc_cnt = (TProfile*)file->Get(Form("tp1f_reso%d_BBC_CNT_cent%d",harmonic,centbin));
  TProfile* tp1f_cnt_fvtx = (TProfile*)file->Get(Form("tp1f_reso%d_CNT_FVTX_cent%d",harmonic,centbin));

  float float_bbc_fvtx = tp1f_bbc_fvtx->GetBinContent(1);
  float float_bbc_cnt = tp1f_bbc_cnt->GetBinContent(1);
  float float_cnt_fvtx = tp1f_cnt_fvtx->GetBinContent(1);

  cout << "bbc-fvtxs " << float_bbc_fvtx << endl;
  cout << "bbc-cnt " << float_bbc_cnt << endl;
  cout << "cnt-fvtxs " << float_cnt_fvtx << endl;

  if ( float_bbc_cnt < 0 ) { cout << "YOU'RE GONNA DIE (bbcs-cnt)" << endl; }
  if ( float_cnt_fvtx < 0 ) { cout << "YOU'RE GONNA DIE (cnt-fvtxs)" << endl; }

  float reso_bbc = sqrt((float_bbc_cnt*float_bbc_fvtx)/float_cnt_fvtx); // BCBS/CS
  float reso_fvtx = sqrt((float_cnt_fvtx*float_bbc_fvtx)/float_bbc_cnt); // CSBS/BC

  cout << "bbc resolution is " << reso_bbc << endl;
  cout << "fvtx resolution is " << reso_fvtx << endl;

  // ---



  // ---

  if ( energy == 20 )
    {
      // --- i like chocolate and fudge and other tasty deserts
      cout << "Now making empirical adjustment for 20 GeV data" << endl;
      reso_bbc = 0.015;
      reso_fvtx = 0.04;
    }

  // ---

  TProfile* hvn_fvtxs = (TProfile*)file->Get(Form("fvtxs_v%d_both_docalib_cent%d",harmonic,centbin));
  hvn_fvtxs->Scale(1.0/reso_fvtx);
  hvn_fvtxs->Draw();
  // the 62 GeV is actually 62.4 and the 20 GeV is actually 19.6, so need to modify
  hvn_fvtxs->SetTitle(Form("d+Au collisions at #sqrt{s_{NN}} = %d GeV",energy));
  if ( energy == 62 ) hvn_fvtxs->SetTitle(Form("d+Au collisions at #sqrt{s_{NN}} = 62.4 GeV"));
  if ( energy == 20 ) hvn_fvtxs->SetTitle(Form("d+Au collisions at #sqrt{s_{NN}} = 19.6 GeV"));
  hvn_fvtxs->SetMaximum(0.17);
  if ( energy == 39 ) hvn_fvtxs->SetMaximum(0.27);
  hvn_fvtxs->SetMinimum(0.0);
  TLine line(0,0,3,0);
  line.SetLineStyle(2);
  line.SetLineWidth(2);
  if ( energy == 20 )
    {
      hvn_fvtxs->SetMaximum(0.5);
      hvn_fvtxs->SetMinimum(-0.2);
      line.Draw();
    }
  if ( harmonic == 3 )
    {
      hvn_fvtxs->SetMaximum(0.08);
      hvn_fvtxs->SetMinimum(-0.02);
      line.Draw();
    }
  hvn_fvtxs->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hvn_fvtxs->GetYaxis()->SetTitle(Form("v_{%d}{EP}",harmonic));
  hvn_fvtxs->GetYaxis()->SetTitleOffset(1.25);


  c1->Print(Form("FigsHarmonicCoefficient/centbins_run16dau%d_v%d_fvtxs_cent%d.pdf",energy,harmonic,centbin));
  c1->Print(Form("FigsHarmonicCoefficient/centbins_run16dau%d_v%d_fvtxs_cent%d.png",energy,harmonic,centbin));

  TProfile* hvn_bbcs = (TProfile*)file->Get(Form("bbcs_v%d_both_docalib_cent%d",harmonic,centbin));
  hvn_bbcs->SetLineColor(kRed);
  hvn_bbcs->Scale(1.0/reso_bbc);
  hvn_bbcs->Draw("same");

  if ( energy == 20 )
    {
      hvn_fvtxs->SetMaximum(1.0);
      hvn_fvtxs->SetMinimum(-0.5);
      line.Draw();
    }

  TLegend* leg = new TLegend(0.18,0.68,0.38,0.88);
  leg->AddEntry(hvn_fvtxs,"Run16 FVTXS","el");
  leg->AddEntry(hvn_bbcs,"Run16 BBCS","el");
  leg->SetTextSize(0.05);
  leg->Draw();

  c1->Print(Form("FigsHarmonicCoefficient/centbins_run16dau%d_v%d_fvtxsbbcs_cent%d.pdf",energy,harmonic,centbin));
  c1->Print(Form("FigsHarmonicCoefficient/centbins_run16dau%d_v%d_fvtxsbbcs_cent%d.png",energy,harmonic,centbin));





  c1->Clear();



  // --- now eta plots

  TProfile* hvneta_bbcs = (TProfile*)file->Get(Form("bbcs_v%deta_both_docalib_cent%d",harmonic,centbin));
  hvneta_bbcs->SetLineColor(kBlack);
  hvneta_bbcs->Scale(1.0/reso_bbc);
  hvneta_bbcs->Draw();
  hvneta_bbcs->SetTitle(Form("d+Au collisions at #sqrt{s_{NN}} = %d GeV",energy));
  if ( energy == 62 ) hvneta_bbcs->SetTitle(Form("d+Au collisions at #sqrt{s_{NN}} = 62.4 GeV"));
  if ( energy == 20 ) hvneta_bbcs->SetTitle(Form("d+Au collisions at #sqrt{s_{NN}} = 19.6 GeV"));
  hvneta_bbcs->SetMaximum(0.1);
  hvneta_bbcs->SetMinimum(0.0);
  if ( energy <= 39 || harmonic == 3 ) hvneta_bbcs->SetMinimum(-0.02);
  if ( energy == 20 )
    {
      hvneta_bbcs->SetMaximum(0.3);
      hvneta_bbcs->SetMinimum(-0.2);
      TLine line(-3.1,0,3.1,0);
      line.SetLineStyle(2);
      line.SetLineWidth(2);
      line.Draw();
    }
  hvneta_bbcs->GetXaxis()->SetTitle("#eta");
  hvneta_bbcs->GetYaxis()->SetTitle(Form("v_{%d}{EP}",harmonic));
  hvneta_bbcs->GetYaxis()->SetTitleOffset(1.25);
  TLegend* leta = new TLegend(0.68,0.68,0.88,0.88);
  leta->AddEntry(hvneta_bbcs,"BBCS","el");
  leta->SetTextSize(0.05);
  leta->Draw();

  c1->Print(Form("FigsHarmonicCoefficient/centbins_run16dau%d_v%deta_bbcs_cent%d.pdf",energy,harmonic,centbin));
  c1->Print(Form("FigsHarmonicCoefficient/centbins_run16dau%d_v%deta_bbcs_cent%d.png",energy,harmonic,centbin));

  TProfile* hvneta_fvtxs = (TProfile*)file->Get(Form("fvtxs_v%deta_both_docalib_cent%d",harmonic,centbin));
  hvneta_fvtxs->SetLineColor(kBlue);
  hvneta_fvtxs->Scale(1.0/reso_fvtx);
  hvneta_fvtxs->Draw("same");
  leta->AddEntry(hvneta_fvtxs,"FVTXS","el");
  leta->Draw();

  c1->Print(Form("FigsHarmonicCoefficient/centbins_run16dau%d_v%deta_fvtxsbbcs_cent%d.pdf",energy,harmonic,centbin));
  c1->Print(Form("FigsHarmonicCoefficient/centbins_run16dau%d_v%deta_fvtxsbbcs_cent%d.png",energy,harmonic,centbin));

  // --- clean
  hvneta_bbcs->SetBinContent(1,-9);
  hvneta_bbcs->SetBinContent(12,-9);
  hvneta_bbcs->SetBinContent(21,-9);
  hvneta_bbcs->SetBinContent(32,-9);
  hvneta_fvtxs->SetBinContent(21,-9);
  hvneta_fvtxs->SetBinContent(32,-9);

  c1->Print(Form("FigsHarmonicCoefficient/centbins_run16dau%d_v%deta_cleaned_fvtxsbbcs_cent%d.pdf",energy,harmonic,centbin));
  c1->Print(Form("FigsHarmonicCoefficient/centbins_run16dau%d_v%deta_cleaned_fvtxsbbcs_cent%d.png",energy,harmonic,centbin));

  delete c1;

  if ( harmonic == 3 ) hvn_fvtxs = xhvn_fvtxs;

  hvneta_fvtxs->SetTitle(Form("d+Au collisions at #sqrt{s_{NN}} = %d GeV",energy));
  if ( energy == 62 ) hvneta_fvtxs->SetTitle(Form("d+Au collisions at #sqrt{s_{NN}} = 62.4 GeV"));
  if ( energy == 20 ) hvneta_fvtxs->SetTitle(Form("d+Au collisions at #sqrt{s_{NN}} = 19.6 GeV"));

}


