void getstuff(int, int, double&, double&, double&, double&, double&, double&, double&, double&, double&, double&, double&, double&);

void getmult(int, double&, double&, double&, double&, double&, double&);

void getmult_new(int, double&, double&, double&, double&); // bbcs charge, cnt tracks, fvtx tracks, fvtx clusters

void makeplots(int, int);

void makemult(int);

void makemult_new(int);

void someratios(int);

void someplots(int);



void temp_runbyrun_fullsuite()
{

  // makemult_new(200);
  // makemult_new(62);
  makemult_new(39);
  makemult_new(20);

  return;

  // makeplots(200,2);
  // makeplots(200,3);
  // makeplots(200,42);

  // makeplots(62,2);
  // makeplots(62,3);
  // makeplots(62,42);

  // makeplots(20,2);
  // makeplots(20,3);
  // makeplots(20,42);

  // makeplots(39,2);
  // makeplots(39,3);
  // makeplots(39,42);

  // someratios(200);
  // someratios(62);
  someratios(20);
  someratios(39);

  // someplots(200);
  // someplots(62);
  // someplots(20);
  // someplots(39);

  // makemult(200);
  // makemult(62);
  // makemult(20);
  // makemult(39);

}



void someplots(int energy)
{

  TCanvas* c1 = new TCanvas("c1","");
  //c1->SetLogz();

  ifstream fin((const char*)Form("list_%d.short",energy));
  int run;
  while ( fin >> run )
    {
      TFile* file = TFile::Open(Form("RootFiles/svrb_run%d_pass0.root",run));
      if ( !file )
        {
          cout << "YOU'RE GONNA DIE" << endl;
          continue;
        }

      TH2D* th2d_corr_bbcqn_bbcqs = (TH2D*)file->Get("th2d_corr_bbcqn_bbcqs");
      TH2D* th2d_corr_bbcqs_fvtxs = (TH2D*)file->Get("th2d_corr_bbcqs_fvtxs");
      TH2D* th2d_corr_bbcqn_fvtxn = (TH2D*)file->Get("th2d_corr_bbcqn_fvtxn");
      TH2D* th2d_corr_bbcqn_fvtxs = (TH2D*)file->Get("th2d_corr_bbcqn_fvtxs");
      TH2D* th2d_corr_bbcqs_fvtxn = (TH2D*)file->Get("th2d_corr_bbcqs_fvtxn");
      TH2D* th2d_corr_fvtxn_fvtxs = (TH2D*)file->Get("th2d_corr_fvtxn_fvtxs");

      TH2D* th2d_corr_npc1_bbcqn = (TH2D*)file->Get("th2d_corr_npc1_bbcqn");
      TH2D* th2d_corr_npc1_bbcqs = (TH2D*)file->Get("th2d_corr_npc1_bbcqs");
      TH2D* th2d_corr_npc1_fvtxn = (TH2D*)file->Get("th2d_corr_npc1_fvtxn");
      TH2D* th2d_corr_npc1_fvtxs = (TH2D*)file->Get("th2d_corr_npc1_fvtxs");

      int nevents = th2d_corr_bbcqn_bbcqs->GetEntries();
      if ( nevents < 1 ) nevents = 1;

      th2d_corr_bbcqn_bbcqs->Draw("colz");
      th2d_corr_bbcqn_bbcqs->Scale(1.0/nevents);
      th2d_corr_bbcqn_bbcqs->GetXaxis()->SetTitle("BBC North charge");
      th2d_corr_bbcqn_bbcqs->GetYaxis()->SetTitle("BBC South charge");
      c1->Print(Form("FigsOther/plot2dcorr_bbcqn_bbcqs_%d.png",run));
      c1->Print(Form("FigsOther/plot2dcorr_bbcqn_bbcqs_%d.pdf",run));
      th2d_corr_bbcqs_fvtxs->Draw("colz");
      th2d_corr_bbcqs_fvtxs->Scale(1.0/nevents);
      th2d_corr_bbcqs_fvtxs->GetXaxis()->SetTitle("BBC South charge");
      th2d_corr_bbcqs_fvtxs->GetYaxis()->SetTitle("FVTX South cluster multiplicity");
      c1->Print(Form("FigsOther/plot2dcorr_bbcqs_fvtxs_%d.png",run));
      c1->Print(Form("FigsOther/plot2dcorr_bbcqs_fvtxs_%d.pdf",run));
      th2d_corr_bbcqn_fvtxn->Draw("colz");
      th2d_corr_bbcqn_fvtxn->Scale(1.0/nevents);
      th2d_corr_bbcqn_fvtxn->GetXaxis()->SetTitle("BBC North charge");
      th2d_corr_bbcqn_fvtxn->GetYaxis()->SetTitle("FVTX North cluster multiplicity");
      c1->Print(Form("FigsOther/plot2dcorr_bbcqn_fvtxn_%d.png",run));
      c1->Print(Form("FigsOther/plot2dcorr_bbcqn_fvtxn_%d.pdf",run));
      th2d_corr_bbcqn_fvtxs->Draw("colz");
      th2d_corr_bbcqn_fvtxs->Scale(1.0/nevents);
      th2d_corr_bbcqn_fvtxs->GetXaxis()->SetTitle("BBC North charge");
      th2d_corr_bbcqn_fvtxs->GetYaxis()->SetTitle("FVTX South cluster multiplicity");
      c1->Print(Form("FigsOther/plot2dcorr_bbcqn_fvtxs_%d.png",run));
      c1->Print(Form("FigsOther/plot2dcorr_bbcqn_fvtxs_%d.pdf",run));
      th2d_corr_bbcqs_fvtxn->Draw("colz");
      th2d_corr_bbcqs_fvtxn->Scale(1.0/nevents);
      th2d_corr_bbcqs_fvtxn->GetXaxis()->SetTitle("BBC South charge");
      th2d_corr_bbcqs_fvtxn->GetYaxis()->SetTitle("FVTX North cluster multiplicity");
      c1->Print(Form("FigsOther/plot2dcorr_bbcqs_fvtxn_%d.png",run));
      c1->Print(Form("FigsOther/plot2dcorr_bbcqs_fvtxn_%d.pdf",run));
      th2d_corr_fvtxn_fvtxs->Draw("colz");
      th2d_corr_fvtxn_fvtxs->Scale(1.0/nevents);
      th2d_corr_fvtxn_fvtxs->GetXaxis()->SetTitle("FVTX North cluster multiplicity");
      th2d_corr_fvtxn_fvtxs->GetYaxis()->SetTitle("FVTX South cluster multiplicity");
      c1->Print(Form("FigsOther/plot2dcorr_fvtxn_fvtxs_%d.png",run));
      c1->Print(Form("FigsOther/plot2dcorr_fvtxn_fvtxs_%d.pdf",run));

      th2d_corr_npc1_bbcqn->Draw("colz");
      th2d_corr_npc1_bbcqn->Scale(1.0/nevents);
      th2d_corr_npc1_bbcqn->GetXaxis()->SetTitle("PC1 cluster multiplicity");
      th2d_corr_npc1_bbcqn->GetYaxis()->SetTitle("BBC North charge");
      c1->Print(Form("FigsOther/plot2dcorr_npc1_bbcqn_%d.png",run));
      c1->Print(Form("FigsOther/plot2dcorr_npc1_bbcqn_%d.pdf",run));
      th2d_corr_npc1_bbcqs->Draw("colz");
      th2d_corr_npc1_bbcqs->Scale(1.0/nevents);
      th2d_corr_npc1_bbcqs->GetXaxis()->SetTitle("PC1 cluster multiplicity");
      th2d_corr_npc1_bbcqs->GetYaxis()->SetTitle("BBC South charge");
      c1->Print(Form("FigsOther/plot2dcorr_npc1_bbcqs_%d.png",run));
      c1->Print(Form("FigsOther/plot2dcorr_npc1_bbcqs_%d.pdf",run));
      th2d_corr_npc1_fvtxn->Draw("colz");
      th2d_corr_npc1_fvtxn->Scale(1.0/nevents);
      th2d_corr_npc1_fvtxn->GetXaxis()->SetTitle("PC1 cluster multiplicity");
      th2d_corr_npc1_fvtxn->GetYaxis()->SetTitle("FVTX North cluster multiplicity");
      c1->Print(Form("FigsOther/plot2dcorr_npc1_fvtxn_%d.png",run));
      c1->Print(Form("FigsOther/plot2dcorr_npc1_fvtxn_%d.pdf",run));
      th2d_corr_npc1_fvtxs->Draw("colz");
      th2d_corr_npc1_fvtxs->Scale(1.0/nevents);
      th2d_corr_npc1_fvtxs->GetXaxis()->SetTitle("PC1 cluster multiplicity");
      th2d_corr_npc1_fvtxs->GetYaxis()->SetTitle("FVTX South cluster multiplicity");
      c1->Print(Form("FigsOther/plot2dcorr_npc1_fvtxs_%d.png",run));
      c1->Print(Form("FigsOther/plot2dcorr_npc1_fvtxs_%d.pdf",run));

      file->Close();
    }

}



void someratios(int energy)
{

  TCanvas* c1 = new TCanvas("c1","");

  ifstream fin((const char*)Form("runbyrun_fromlogfiles_%d.dat",energy));
  int run[115];
  int counter = 0;
  float index[115];
  float cratio[115];
  float vratio[115];
  float eratio[115];
  while ( fin >> run[counter] >> cratio[counter] >> vratio[counter] >> eratio[counter] )
    {
      index[counter] = counter + 0.5;
      eratio[counter] = 1 - eratio[counter];
      ++counter;
    }

  TGraph* tg_cratio = new TGraph(counter,index,cratio);
  TGraph* tg_vratio = new TGraph(counter,index,vratio);
  TGraph* tg_eratio = new TGraph(counter,index,eratio);

  tg_cratio->SetMarkerStyle(kOpenCircle);
  tg_vratio->SetMarkerStyle(kOpenSquare);
  tg_eratio->SetMarkerStyle(kOpenCross);

  tg_cratio->SetMarkerColor(kBlack);
  tg_vratio->SetMarkerColor(kBlue);
  tg_eratio->SetMarkerColor(kRed);

  tg_cratio->GetXaxis()->SetLimits(-1,counter);
  tg_vratio->GetXaxis()->SetLimits(-1,counter);
  tg_eratio->GetXaxis()->SetLimits(-1,counter);

  tg_cratio->SetMaximum(1.0);
  tg_vratio->SetMaximum(1.0);
  tg_eratio->SetMaximum(1.0);

  tg_cratio->SetMinimum(0.0);
  tg_vratio->SetMinimum(0.0);
  tg_eratio->SetMinimum(0.0);

  tg_cratio->Draw("ap");
  tg_vratio->Draw("p");
  tg_eratio->Draw("p");
  TLegend* leg3 = new TLegend(0.48,0.18,0.88,0.38);
  leg3->AddEntry(tg_cratio,"Ratio of clusters within 5.2 cm","p");
  leg3->AddEntry(tg_vratio,"Ratio of events with vertex outside/inside 5 cm","p");
  leg3->AddEntry(tg_eratio,"Ratio of events passing cluster cut","p");
  leg3->SetFillStyle(0);
  leg3->Draw();
  c1->Print(Form("FigsOther/ratios_runbyrun_fromlogfiles_%d.png",energy));
  c1->Print(Form("FigsOther/ratios_runbyrun_fromlogfiles_%d.pdf",energy));

  tg_cratio->Draw("ap");
  c1->Print(Form("FigsOther/ratio_cluster_runbyrun_fromlogfiles_%d.png",energy));
  c1->Print(Form("FigsOther/ratio_cluster_runbyrun_fromlogfiles_%d.pdf",energy));

  tg_vratio->Draw("ap");
  c1->Print(Form("FigsOther/ratio_vertex_runbyrun_fromlogfiles_%d.png",energy));
  c1->Print(Form("FigsOther/ratio_vertex_runbyrun_fromlogfiles_%d.pdf",energy));

  tg_eratio->Draw("ap");
  c1->Print(Form("FigsOther/ratio_pass_runbyrun_fromlogfiles_%d.png",energy));
  c1->Print(Form("FigsOther/ratio_pass_runbyrun_fromlogfiles_%d.pdf",energy));

  tg_eratio->Draw("ap");
  tg_vratio->Draw("p");
  TLegend* leg2 = new TLegend(0.48,0.18,0.88,0.38);
  leg2->AddEntry(tg_vratio,"Ratio of events with vertex outside/inside 5 cm","p");
  leg2->AddEntry(tg_eratio,"Ratio of events passing cluster cut","p");
  leg2->SetFillStyle(0);
  leg2->Draw();
  c1->Print(Form("FigsOther/ratio_vertexandpass_runbyrun_fromlogfiles_%d.png",energy));
  c1->Print(Form("FigsOther/ratio_vertexandpass_runbyrun_fromlogfiles_%d.pdf",energy));

  delete c1;

  delete leg3;
  delete leg2;

}



void makeplots(int energy, int harmonic)
{

  if ( harmonic != 2 ) return;

  gStyle->SetOptTitle(0);

  TCanvas* c1 = new TCanvas("c1","");

  // -----------------------------------------
  // --- first get the individual correlations
  // -----------------------------------------
  double index[110];
  double
    array_BC[110],
    array_BN[110],
    array_BS[110],
    array_CN[110],
    array_CS[110],
    array_NS[110],
    array_eBC[110],
    array_eBN[110],
    array_eBS[110],
    array_eCN[110],
    array_eCS[110],
    array_eNS[110];
  int run;
  int counter = 0;
  ifstream fin((const char*)Form("list_%d.short",energy));
  for ( int i = 0; i < 110; ++i )
    {
      index[i] = i + 0.5;
      double
        temp_BC,
        temp_BN,
        temp_BS,
        temp_CN,
        temp_CS,
        temp_NS,
        temp_eBC,
        temp_eBN,
        temp_eBS,
        temp_eCN,
        temp_eCS,
        temp_eNS;
      if ( fin.eof() ) break;
      fin >> run;
      //cout << run << endl;
      ++counter;
      getstuff(run,harmonic,
               temp_BC,
               temp_BN,
               temp_BS,
               temp_CN,
               temp_CS,
               temp_NS,
               temp_eBC,
               temp_eBN,
               temp_eBS,
               temp_eCN,
               temp_eCS,
               temp_eNS
               ); // end getstuff
        array_BC[i] = temp_BC;
        array_BN[i] = temp_BN;
        array_BS[i] = temp_BS;
        array_CN[i] = temp_CN;
        array_CS[i] = temp_CS;
        array_NS[i] = temp_NS;
        array_eBC[i] = temp_eBC;
        array_eBN[i] = temp_eBN;
        array_eBS[i] = temp_eBS;
        array_eCN[i] = temp_eCN;
        array_eCS[i] = temp_eCS;
        array_eNS[i] = temp_eNS;

        // cout << "Testing the array values..." << endl;
        // cout << array_BC[i] << endl;
        // cout << array_BN[i] << endl;
        // cout << array_BS[i] << endl;
        // cout << array_CN[i] << endl;
        // cout << array_CS[i] << endl;
        // cout << array_NS[i] << endl;
        // cout << "Testing the error values..." << endl;
        // cout << array_eBC[i] << endl;
        // cout << array_eBN[i] << endl;
        // cout << array_eBS[i] << endl;
        // cout << array_eCN[i] << endl;
        // cout << array_eCS[i] << endl;
        // cout << array_eNS[i] << endl;

    } // end for loop



  // -------------------------------------------------------
  // --- now get the resolutions built from the correlations
  // -------------------------------------------------------

  double resoB_CN[110]; // BBCS(B) resolution using CNT(C) and FVTXN(N)
  double resoB_CS[110]; // BBCS(B) resolution using CNT(C) and FVTXS(S)
  double resoB_NS[110]; // BBCS(B) resolution using FVTXN(N) and FVTXS(S)

  double resoN_BC[110]; // FVTXN(N) resolution using BBCS(B) and CNT(C)
  double resoN_BS[110]; // FVTXN(N) resolution using BBCS(B) and FVTXS(S)
  double resoN_CS[110]; // FVTXN(N) resolution using CNT(C) and FVTXS(S)

  double resoS_BC[110]; // FVTXS(S) resolution using BBCS(B) and CNT(C)
  double resoS_BN[110]; // FVTXS(S) resolution using BBCS(B) and CNT(C)
  double resoS_CN[110]; // FVTXS(S) resolution using CNT(C) and FVTXN(N)

  double eresoB_CN[110]; // BBCS(B) resolution using CNT(C) and FVTXN(N)
  double eresoB_CS[110]; // BBCS(B) resolution using CNT(C) and FVTXS(S)
  double eresoB_NS[110]; // BBCS(B) resolution using FVTXN(N) and FVTXS(S)

  double eresoN_BC[110]; // FVTXN(N) resolution using BBCS(B) and CNT(C)
  double eresoN_BS[110]; // FVTXN(N) resolution using BBCS(B) and FVTXS(S)
  double eresoN_CS[110]; // FVTXN(N) resolution using CNT(C) and FVTXS(S)

  double eresoS_BC[110]; // FVTXS(S) resolution using BBCS(B) and CNT(C)
  double eresoS_BN[110]; // FVTXS(S) resolution using BBCS(B) and CNT(C)
  double eresoS_CN[110]; // FVTXS(S) resolution using CNT(C) and FVTXN(N)

  for ( int i = 0; i < counter-1; ++i )
    {
      resoB_CN[i] = sqrt ( ( array_BC[i] * array_BN[i] ) / array_CN[i] );
      resoB_CS[i] = sqrt ( ( array_BC[i] * array_BS[i] ) / array_CS[i] );
      resoB_NS[i] = sqrt ( ( array_BN[i] * array_BS[i] ) / array_NS[i] );

      resoN_BC[i] = sqrt ( ( array_BN[i] * array_CN[i] ) / array_BC[i] );
      resoN_BS[i] = sqrt ( ( array_BN[i] * array_NS[i] ) / array_BS[i] );
      resoN_CS[i] = sqrt ( ( array_BC[i] * array_BS[i] ) / array_CS[i] );

      resoS_BC[i] = sqrt ( ( array_BS[i] * array_CS[i] ) / array_BC[i] );
      resoS_BN[i] = sqrt ( ( array_BS[i] * array_NS[i] ) / array_BN[i] );
      resoS_CN[i] = sqrt ( ( array_CS[i] * array_NS[i] ) / array_CN[i] );

      // --- propagation of uncertainties, wee fun...

      eresoB_CN[i] = resoB_CN[i] * sqrt ( pow(array_eBC[i]/array_BC[i],2) + pow(array_eBN[i]/array_BN[i],2) + pow(array_eCN[i]/array_CN[i],2) );
      eresoB_CS[i] = resoB_CS[i] * sqrt ( pow(array_eBC[i]/array_BC[i],2) + pow(array_eBS[i]/array_BS[i],2) + pow(array_eCS[i]/array_CS[i],2) );
      eresoB_NS[i] = resoB_NS[i] * sqrt ( pow(array_eBN[i]/array_BN[i],2) + pow(array_eBS[i]/array_BS[i],2) + pow(array_eNS[i]/array_NS[i],2) );

      eresoN_BC[i] = resoN_BC[i] * sqrt ( pow(array_eBN[i]/array_BN[i],2) + pow(array_eCN[i]/array_CN[i],2) + pow(array_eBC[i]/array_BC[i],2) );
      eresoN_BS[i] = resoN_BS[i] * sqrt ( pow(array_eBN[i]/array_BN[i],2) + pow(array_eNS[i]/array_NS[i],2) + pow(array_eBS[i]/array_BS[i],2) );
      eresoN_CS[i] = resoN_CS[i] * sqrt ( pow(array_eBC[i]/array_BC[i],2) + pow(array_eBS[i]/array_BS[i],2) + pow(array_eCS[i]/array_CS[i],2) );

      eresoS_BC[i] = resoS_BC[i] * sqrt ( pow(array_eBS[i]/array_BS[i],2) + pow(array_eCS[i]/array_CS[i],2) + pow(array_eBC[i]/array_BC[i],2) );
      eresoS_BN[i] = resoS_BN[i] * sqrt ( pow(array_eBS[i]/array_BS[i],2) + pow(array_eNS[i]/array_NS[i],2) + pow(array_eBN[i]/array_BN[i],2) );
      eresoS_CN[i] = resoS_CN[i] * sqrt ( pow(array_eCS[i]/array_CS[i],2) + pow(array_eNS[i]/array_NS[i],2) + pow(array_eCN[i]/array_CN[i],2) );

      // ---

      // cout << "Resolutions before resetting nans" << endl;

      // cout << resoB_CN[i] << endl;
      // cout << resoB_CS[i] << endl;
      // cout << resoB_NS[i] << endl;

      if ( TMath::IsNaN(resoB_CN[i]) ) { resoB_CN[i] = -9; eresoB_CN[i] = 999; }
      if ( TMath::IsNaN(resoB_CS[i]) ) { resoB_CS[i] = -9; eresoB_CS[i] = 999; }
      if ( TMath::IsNaN(resoB_NS[i]) ) { resoB_NS[i] = -9; eresoB_NS[i] = 999; }

      if ( TMath::IsNaN(resoN_BC[i]) ) { resoN_BC[i] = -9; eresoN_BC[i] = 999; }
      if ( TMath::IsNaN(resoN_BS[i]) ) { resoN_BS[i] = -9; eresoN_BS[i] = 999; }
      if ( TMath::IsNaN(resoN_CS[i]) ) { resoN_CS[i] = -9; eresoN_CS[i] = 999; }

      if ( TMath::IsNaN(resoS_BC[i]) ) { resoS_BC[i] = -9; eresoS_BC[i] = 999; }
      if ( TMath::IsNaN(resoS_BN[i]) ) { resoS_BN[i] = -9; eresoS_BN[i] = 999; }
      if ( TMath::IsNaN(resoS_CN[i]) ) { resoS_CN[i] = -9; eresoS_CN[i] = 999; }

      // cout << "Resolutions after resetting nans" << endl;

      // cout << resoB_CN[i] << endl;
      // cout << resoB_CS[i] << endl;
      // cout << resoB_NS[i] << endl;

    }

  TGraphErrors* tgeresoB_CN = new TGraphErrors(counter-1,index,resoB_CN,0,eresoB_CN);
  TGraphErrors* tgeresoB_CS = new TGraphErrors(counter-1,index,resoB_CS,0,eresoB_CS);
  TGraphErrors* tgeresoB_NS = new TGraphErrors(counter-1,index,resoB_NS,0,eresoB_NS);

  TGraphErrors* tgeresoN_BC = new TGraphErrors(counter-1,index,resoN_BC,0,eresoN_BC);
  TGraphErrors* tgeresoN_BS = new TGraphErrors(counter-1,index,resoN_BS,0,eresoN_BS);
  TGraphErrors* tgeresoN_CS = new TGraphErrors(counter-1,index,resoN_CS,0,eresoN_CS);

  TGraphErrors* tgeresoS_BC = new TGraphErrors(counter-1,index,resoS_BC,0,eresoS_BC);
  TGraphErrors* tgeresoS_BN = new TGraphErrors(counter-1,index,resoS_BN,0,eresoS_BN);
  TGraphErrors* tgeresoS_CN = new TGraphErrors(counter-1,index,resoS_CN,0,eresoS_CN);

  // ---

  tgeresoB_CN->SetMarkerColor(kRed);
  tgeresoB_CS->SetMarkerColor(kGreen+2);
  tgeresoB_NS->SetMarkerColor(kBlue);

  tgeresoB_CN->SetMarkerStyle(kOpenCircle);
  tgeresoB_CS->SetMarkerStyle(kOpenSquare);
  tgeresoB_NS->SetMarkerStyle(kOpenDiamond);

  tgeresoB_CN->Draw("ap");
  tgeresoB_CS->Draw("p");
  tgeresoB_NS->Draw("p");

  tgeresoB_CN->GetXaxis()->SetLimits(-2,counter+1);
  tgeresoB_CN->GetXaxis()->SetTitle("run index");
  tgeresoB_CN->GetYaxis()->SetTitle("BBCS EP resolution");
  tgeresoB_CN->SetTitle("BBCS EP resolution");
  tgeresoB_CN->SetMaximum(0.2);
  tgeresoB_CN->SetMinimum(0.0);
  if ( energy != 200 || harmonic != 2 ) tgeresoB_CN->SetMinimum(-0.05);

  TF1* funB_CN = new TF1("funB_CN","pol0",0,110);
  TF1* funB_CS = new TF1("funB_CS","pol0",0,110);
  TF1* funB_NS = new TF1("funB_NS","pol0",0,110);

  // funB_CN->SetParLimits(0,0,1);
  // funB_CS->SetParLimits(0,0,1);
  // funB_NS->SetParLimits(0,0,1);

  funB_CN->SetLineColor(kRed);
  funB_CS->SetLineColor(kGreen+2);
  funB_NS->SetLineColor(kBlue);

  tgeresoB_CN->Fit(funB_CN,"Q","",0,counter-1);
  tgeresoB_CS->Fit(funB_CS,"Q","",0,counter-1);
  tgeresoB_NS->Fit(funB_NS,"Q","",0,counter-1);

  TLegend* legB = new TLegend(0.18,0.68,0.28,0.88);
  TLegendEntry* tleB_CN = legB->AddEntry(tgeresoB_CN,Form("BBCS, CNT, FVTXN,   average = %f",funB_CN->GetParameter(0)),"p");
  TLegendEntry* tleB_CS = legB->AddEntry(tgeresoB_CS,Form("BBCS, CNT, FVTXS,   average = %f",funB_CS->GetParameter(0)),"p");
  TLegendEntry* tleB_NS = legB->AddEntry(tgeresoB_NS,Form("BBCS, FVTXN, FVTXS, average = %f",funB_NS->GetParameter(0)),"p");
  tleB_CN->SetTextColor(kRed);
  tleB_CS->SetTextColor(kGreen+2);
  tleB_NS->SetTextColor(kBlue);
  legB->SetTextSize(0.05);
  legB->Draw();

  c1->Print(Form("FigsEventPlane/figreso_bbcs_energy%d_harmonic%d.png",energy,harmonic));
  c1->Print(Form("FigsEventPlane/figreso_bbcs_energy%d_harmonic%d.pdf",energy,harmonic));

  // ---

  tgeresoN_BC->SetMarkerColor(kRed);
  tgeresoN_BS->SetMarkerColor(kGreen+2);
  tgeresoN_CS->SetMarkerColor(kBlue);

  tgeresoN_BC->SetMarkerStyle(kOpenCircle);
  tgeresoN_BS->SetMarkerStyle(kOpenSquare);
  tgeresoN_CS->SetMarkerStyle(kOpenDiamond);

  tgeresoN_BC->Draw("ap");
  tgeresoN_BS->Draw("p");
  tgeresoN_CS->Draw("p");

  tgeresoN_BC->GetXaxis()->SetLimits(-2,counter+1);
  tgeresoN_BC->GetXaxis()->SetTitle("run index");
  tgeresoN_BC->GetYaxis()->SetTitle("FVTXN EP resolution");
  tgeresoN_BC->SetTitle("FVTXN EP resolution");
  tgeresoN_BC->SetMaximum(0.2);
  tgeresoN_BC->SetMinimum(0.0);
  if ( energy != 200 || harmonic != 2 ) tgeresoN_BC->SetMinimum(-0.05);

  // TLegend* legB = new TLegend(0.18,0.68,0.28,0.88);
  // legB->AddEntry(tgeresoN_BC,"FVTXN, BBCS, CNT","p");
  // legB->AddEntry(tgeresoN_BS,"FVTXN, BBCS, FVTXS","p");
  // legB->AddEntry(tgeresoN_CS,"FVTXN, CNT, FVTXS","p");
  // legB->SetTextSize(0.05);
  // legB->Draw();

  TF1* funN_BC = new TF1("funN_BC","pol0",0,110);
  TF1* funN_BS = new TF1("funN_BS","pol0",0,110);
  TF1* funN_CS = new TF1("funN_CS","pol0",0,110);

  // funN_BC->SetParLimits(0,0,1);
  // funN_BS->SetParLimits(0,0,1);
  // funN_CS->SetParLimits(0,0,1);

  funN_BC->SetLineColor(kRed);
  funN_BS->SetLineColor(kGreen+2);
  funN_CS->SetLineColor(kBlue);

  tgeresoN_BC->Fit(funN_BC,"Q","",0,counter-1);
  tgeresoN_BS->Fit(funN_BS,"Q","",0,counter-1);
  tgeresoN_CS->Fit(funN_CS,"Q","",0,counter-1);

  TLegend* legN = new TLegend(0.18,0.68,0.28,0.88);
  TLegendEntry* tleN_BC = legN->AddEntry(tgeresoN_BC,Form("FVTXN, BBCS, CNT,   average = %f",funN_BC->GetParameter(0)),"p");
  TLegendEntry* tleN_BS = legN->AddEntry(tgeresoN_BS,Form("FVTXN, BBCS, FVTXS, average = %f",funN_BS->GetParameter(0)),"p");
  TLegendEntry* tleN_CS = legN->AddEntry(tgeresoN_CS,Form("FVTXN, CNT, FVTXS, average = %f",funN_CS->GetParameter(0)),"p");
  tleN_BC->SetTextColor(kRed);
  tleN_BS->SetTextColor(kGreen+2);
  tleN_CS->SetTextColor(kBlue);
  legN->SetTextSize(0.05);
  legN->Draw();

  c1->Print(Form("FigsEventPlane/figreso_fvtxn_energy%d_harmonic%d.png",energy,harmonic));
  c1->Print(Form("FigsEventPlane/figreso_fvtxn_energy%d_harmonic%d.pdf",energy,harmonic));

  // ---

  tgeresoS_BC->SetMarkerColor(kRed);
  tgeresoS_BN->SetMarkerColor(kGreen+2);
  tgeresoS_CN->SetMarkerColor(kBlue);

  tgeresoS_BC->SetMarkerStyle(kOpenCircle);
  tgeresoS_BN->SetMarkerStyle(kOpenSquare);
  tgeresoS_CN->SetMarkerStyle(kOpenDiamond);

  tgeresoS_BC->Draw("ap");
  tgeresoS_BN->Draw("p");
  tgeresoS_CN->Draw("p");

  tgeresoS_BC->GetXaxis()->SetLimits(-2,counter+1);
  tgeresoS_BC->GetXaxis()->SetTitle("run index");
  tgeresoS_BC->GetYaxis()->SetTitle("FVTXS EP resolution");
  tgeresoS_BC->SetTitle("FVTXS EP resolution");
  tgeresoS_BC->SetMaximum(0.4);
  tgeresoS_BC->SetMinimum(0.0);
  if ( energy != 200 || harmonic != 2 ) tgeresoS_BC->SetMinimum(-0.05);

  TF1* funS_BC = new TF1("funS_BC","pol0",0,110);
  TF1* funS_BN = new TF1("funS_BN","pol0",0,110);
  TF1* funS_CN = new TF1("funS_CN","pol0",0,110);

  // funS_BC->SetParLimits(0,0,1);
  // funS_BN->SetParLimits(0,0,1);
  // funS_CN->SetParLimits(0,0,1);

  funS_BC->SetLineColor(kRed);
  funS_BN->SetLineColor(kGreen+2);
  funS_CN->SetLineColor(kBlue);

  tgeresoS_BC->Fit(funS_BC,"Q","",0,counter-1);
  tgeresoS_BN->Fit(funS_BN,"Q","",0,counter-1);
  tgeresoS_CN->Fit(funS_CN,"Q","",0,counter-1);

  TLegend* legS = new TLegend(0.18,0.68,0.28,0.88);
  TLegendEntry* tleS_BC = legS->AddEntry(tgeresoS_BC,Form("FVTXS, BBCS, CNT,   average = %f",funS_BC->GetParameter(0)),"p");
  TLegendEntry* tleS_BN = legS->AddEntry(tgeresoS_BN,Form("FVTXS, BBCS, FVTXN, average = %f",funS_BN->GetParameter(0)),"p");
  TLegendEntry* tleS_CN = legS->AddEntry(tgeresoS_CN,Form("FVTXS, CNT, FVTXN, average = %f",funS_CN->GetParameter(0)),"p");
  tleS_BC->SetTextColor(kRed);
  tleS_BN->SetTextColor(kGreen+2);
  tleS_CN->SetTextColor(kBlue);
  legS->SetTextSize(0.05);
  legS->Draw();

  c1->Print(Form("FigsEventPlane/figreso_fvtxs_energy%d_harmonic%d.png",energy,harmonic));
  c1->Print(Form("FigsEventPlane/figreso_fvtxs_energy%d_harmonic%d.pdf",energy,harmonic));

  // ---

  delete c1;

}



void getstuff(int run, int hh,
              double& double_BC,
              double& double_BN,
              double& double_BS,
              double& double_CN,
              double& double_CS,
              double& double_NS,
              double& double_eBC,
              double& double_eBN,
              double& double_eBS,
              double& double_eCN,
              double& double_eCS,
              double& double_eNS
              )
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

  TProfile* tp1f_BC = (TProfile*)file->Get(Form("tp1f_reso%d_BBC_CNT",hh));
  TProfile* tp1f_BN = (TProfile*)file->Get(Form("tp1f_reso%d_BBC_FVTXN",hh));
  TProfile* tp1f_BS = (TProfile*)file->Get(Form("tp1f_reso%d_BBC_FVTX",hh));
  TProfile* tp1f_CN = (TProfile*)file->Get(Form("tp1f_reso%d_CNT_FVTXN",hh));
  TProfile* tp1f_CS = (TProfile*)file->Get(Form("tp1f_reso%d_CNT_FVTX",hh));
  TProfile* tp1f_NS = (TProfile*)file->Get(Form("tp1f_reso%d_FVTXS_FVTXN",hh));

  double_BC = tp1f_BC->GetBinContent(1);
  double_BN = tp1f_BN->GetBinContent(1);
  double_BS = tp1f_BS->GetBinContent(1);
  double_CN = tp1f_CN->GetBinContent(1);
  double_CS = tp1f_CS->GetBinContent(1);
  double_NS = tp1f_NS->GetBinContent(1);

  double_eBC = tp1f_BC->GetBinError(1);
  double_eBN = tp1f_BN->GetBinError(1);
  double_eBS = tp1f_BS->GetBinError(1);
  double_eCN = tp1f_CN->GetBinError(1);
  double_eCS = tp1f_CS->GetBinError(1);
  double_eNS = tp1f_NS->GetBinError(1);

  file->Close();

}


void makemult(int energy)
{

  gStyle->SetOptTitle(0);

  TCanvas* c1 = new TCanvas("c1","");

  // -----------------------------------------
  // --- first get the individual correlations
  // -----------------------------------------
  int runn[110];
  double index[110];
  double
    array_BBCS[110],
    array_FVTXS[110],
    array_FVTXN[110],
    array_eBBCS[110],
    array_eFVTXS[110],
    array_eFVTXN[110];
  int run;
  int counter = 0;
  ifstream fin((const char*)Form("list_%d.short",energy));
  for ( int i = 0; i < 110; ++i )
    {
      if ( fin.eof() ) break;
      fin >> run;
      // cout << "energy is " << energy << endl;
      // cout << "list name is " << (const char*)Form("list_%d.short",energy) << endl;
      // cout << "run is " << run << endl;
      double temp_bbcs, temp_fvtxs, temp_fvtxn;
      double temp_ebbcs, temp_efvtxs, temp_efvtxn;
      index[i] = i+0.5;
      runn[i] = run;
      if ( !getmult(run,temp_bbcs,temp_fvtxs,temp_fvtxn,temp_ebbcs,temp_efvtxs,temp_efvtxn) ) continue;
      ++counter;
      array_BBCS[i] = temp_bbcs;
      array_FVTXS[i] = temp_fvtxs;
      array_FVTXN[i] = temp_fvtxn;
      array_eBBCS[i] = temp_ebbcs;
      array_eFVTXS[i] = temp_efvtxs;
      array_eFVTXN[i] = temp_efvtxn;
    }
  fin.close();

  TGraphErrors* tge_fvtxs = new TGraphErrors(counter,index,array_FVTXS,0,array_eFVTXS);
  TGraphErrors* tge_fvtxn = new TGraphErrors(counter,index,array_FVTXN,0,array_eFVTXN);
  TGraphErrors* tge_bbcs = new TGraphErrors(counter,index,array_BBCS,0,array_eBBCS);

  tge_bbcs->SetMarkerColor(kBlack);
  tge_fvtxs->SetMarkerColor(kRed);
  tge_fvtxn->SetMarkerColor(kBlue);

  tge_bbcs->SetMarkerStyle(kOpenCircle);
  tge_fvtxs->SetMarkerStyle(kOpenCircle);
  tge_fvtxn->SetMarkerStyle(kOpenCircle);

  tge_fvtxs->Draw("ap");
  tge_fvtxs->SetMaximum(650);
  if ( energy == 39 ) tge_fvtxs->SetMaximum(300);
  tge_fvtxs->SetMinimum(0);
  tge_fvtxs->GetXaxis()->SetLimits(-2,counter+1);
  tge_fvtxs->GetXaxis()->SetTitle("run index");
  tge_fvtxs->GetYaxis()->SetTitle("Mean Multiplicity");
  tge_fvtxn->Draw("p");
  tge_bbcs->Draw("p");

  TLegend* leg = new TLegend(0.78,0.78,0.98,0.98);
  leg->AddEntry(tge_bbcs,"BBCS","p");
  leg->AddEntry(tge_fvtxs,"FVTXS","p");
  leg->AddEntry(tge_fvtxn,"FVTXN","p");
  leg->SetFillStyle(0);
  leg->SetTextSize(0.045);
  leg->Draw();

  TF1* fun_bbcs = new TF1("fun_bbcs","pol0",0,110);
  TF1* fun_fvtxs = new TF1("fun_fvtxs","pol0",0,110);
  TF1* fun_fvtxn = new TF1("fun_fvtxn","pol0",0,110);
  fun_bbcs->SetLineColor(kBlack);
  fun_fvtxs->SetLineColor(kRed);
  fun_fvtxn->SetLineColor(kBlue);

  tge_bbcs->Fit(fun_bbcs,"","",0,counter);
  tge_fvtxs->Fit(fun_fvtxs,"","",0,counter);
  tge_fvtxn->Fit(fun_fvtxn,"","",0,counter);

  // c1->Print(Form("FigsOther/mult_runbyrun_energy%d.png",energy));
  // c1->Print(Form("FigsOther/mult_runbyrun_energy%d.pdf",energy));
  c1->Print(Form("FigsOther/mult_runbyrun_energy%d_IR.png",energy));
  c1->Print(Form("FigsOther/mult_runbyrun_energy%d_IR.pdf",energy));

  // -------------------------------------------------
  // --- now attempting to do simplistic run by run QA
  // -------------------------------------------------

  double ave_BBCS = 0;
  double rms_BBCS = 0;
  double ave_FVTXS = 0;
  double rms_FVTXS = 0;
  double ave_FVTXN = 0;
  double rms_FVTXN = 0;
  for ( int i = 0; i < counter; ++i )
    {
      ave_BBCS += array_BBCS[i];
      rms_BBCS += pow(array_BBCS[i],2);
      ave_FVTXS += array_FVTXS[i];
      rms_FVTXS += pow(array_FVTXS[i],2);
      ave_FVTXN += array_FVTXN[i];
      rms_FVTXN += pow(array_FVTXN[i],2);
    }
  ave_BBCS /= counter;
  rms_BBCS /= counter;
  ave_FVTXS /= counter;
  rms_FVTXS /= counter;
  ave_FVTXN /= counter;
  rms_FVTXN /= counter;
  double fit_BBCS = fun_bbcs->GetParameter(0);
  double fit_FVTXS = fun_fvtxs->GetParameter(0);
  double fit_FVTXN = fun_fvtxn->GetParameter(0);
  for ( int i = 0; i < counter; ++i )
    {
      //if ( fabs(array_BBCS[i]-ave_BBCS) > (rms_BBCS-ave_BBCS) )
      //if ( fabs(array_BBCS[i]-ave_BBCS) > (1.2*ave_BBCS) )
      //if ( fabs(array_BBCS[i]-fit_BBCS) > (1.05*fit_BBCS) )
      if ( array_BBCS[i] > 999 )
        {
          cout << "BBCS out of bounds for runnumber " << runn[i] << endl;
          array_BBCS[i] = 999; // just playing around for visualizations...
          array_eBBCS[i] = 999; // just playing around for visualizations...
        }
      //if ( fabs(array_FVTXS[i]-ave_FVTXS) > (rms_FVTXS-ave_FVTXS) )
      //if ( fabs(array_FVTXS[i]-ave_FVTXS) > (1.2*ave_FVTXS) )
      if ( array_FVTXS[i] > 999 )
        {
          cout << "FVTXS out of bounds for runnumber " << runn[i] << endl;
          array_FVTXS[i] = 999; // just playing around for visualizations...
          array_eFVTXS[i] = 999; // just playing around for visualizations...
        }
      //if ( fabs(array_FVTXN[i]-ave_FVTXN) > (rms_FVTXN-ave_FVTXN) )
      //if ( fabs(array_FVTXN[i]-ave_FVTXN) > (1.2*ave_FVTXN) )
      if ( array_FVTXN[i] > 999 )
        {
          cout << "FVTXN out of bounds for runnumber " << runn[i] << endl;
          array_FVTXN[i] = 999; // just playing around for visualizations...
          array_eFVTXN[i] = 999; // just playing around for visualizations...
        }
    }

  TGraphErrors* tge_fvtxs = new TGraphErrors(counter,index,array_FVTXS,0,array_eFVTXS);
  TGraphErrors* tge_fvtxn = new TGraphErrors(counter,index,array_FVTXN,0,array_eFVTXN);
  TGraphErrors* tge_bbcs = new TGraphErrors(counter,index,array_BBCS,0,array_eBBCS);

  tge_bbcs->SetMarkerColor(kBlack);
  tge_fvtxs->SetMarkerColor(kRed);
  tge_fvtxn->SetMarkerColor(kBlue);

  tge_bbcs->SetMarkerStyle(kOpenCircle);
  tge_fvtxs->SetMarkerStyle(kOpenCircle);
  tge_fvtxn->SetMarkerStyle(kOpenCircle);

  tge_fvtxs->Draw("ap");
  tge_fvtxs->SetMaximum(650);
  if ( energy == 39 ) tge_fvtxs->SetMaximum(300);
  tge_fvtxs->SetMinimum(0);
  tge_fvtxs->GetXaxis()->SetLimits(-2,counter+1);
  tge_fvtxs->GetXaxis()->SetTitle("run index");
  tge_fvtxs->GetYaxis()->SetTitle("Mean Multiplicity");
  tge_fvtxn->Draw("p");
  tge_bbcs->Draw("p");
  leg->Draw();
  TLine line1(0,220,counter,220);
  line1.SetLineStyle(2);
  line1.Draw();
  TLine line2(0,150,counter,150);
  line2.SetLineStyle(2);
  line2.Draw();
  tge_bbcs->Fit(fun_bbcs,"","",0,counter);
  tge_fvtxs->Fit(fun_fvtxs,"","",0,counter);
  tge_fvtxn->Fit(fun_fvtxn,"","",0,counter);
  // c1->Print(Form("FigsOther/mult_runbyrun_energy%d_rejects.png",energy));
  // c1->Print(Form("FigsOther/mult_runbyrun_energy%d_rejects.pdf",energy));
  c1->Print(Form("FigsOther/mult_runbyrun_energy%d_rejects_IR.png",energy));
  c1->Print(Form("FigsOther/mult_runbyrun_energy%d_rejects_IR.pdf",energy));

  int runj;
  double* ratio = new double[counter];
  ifstream finj((const char*)Form("runratio_%d.dat",energy));
  for ( int i = 0; i < counter-1; ++i )
    {
      finj >> runj >> ratio[i];
      //if ( ratio[i] > 30.0 || array_FVTXS[i] > 200.0 ) cout << runj << " " << runn[i] << " " << array_FVTXS[i] << " " << ratio[i] << endl;
      //if ( ratio[i] > 30.0 ) cout << runj << " " << runn[i] << " " << array_FVTXS[i] << " " << ratio[i] << endl;
      cout << runj << " " << runn[i] << " " << array_FVTXS[i] << " " << ratio[i] << endl;
    }
  cout << "counter is " << counter << endl;

  TGraph* tgj = new TGraph(84,array_FVTXS,ratio);
  tgj->SetMarkerColor(kBlack);
  tgj->SetMarkerStyle(kOpenCircle);
  tgj->Draw("ap");
  tgj->GetXaxis()->SetTitle("FVTXS mean mult");
  tgj->GetYaxis()->SetTitle("ratio of triggers");
  if ( energy == 20 )
    {
      TLine line1(0.0,30.0,200.0,30.0);
      line1.SetLineStyle(2);
      line1.Draw();
      TLine line2(200.0,0.0,200.0,30.0);
      line2.SetLineStyle(2);
      line2.Draw();
    }
  if ( energy == 39 )
    {
      TLine line1(0.0,6.0,151.0,6.0);
      line1.SetLineStyle(2);
      line1.Draw();
      TLine line2(151.0,0.0,151.0,6.0);
      line2.SetLineStyle(2);
      line2.Draw();
    }
  c1->Print(Form("FigsOther/mult_runbyrun_energy%d_ratioj.png",energy));
  c1->Print(Form("FigsOther/mult_runbyrun_energy%d_ratioj.pdf",energy));
  if ( energy == 20 )
    {
      tgj->GetXaxis()->SetLimits(0,200);
      tgj->SetMaximum(30);
      c1->Print(Form("FigsOther/mult_runbyrun_energy%d_ratioj_zoom.png",energy));
      c1->Print(Form("FigsOther/mult_runbyrun_energy%d_ratioj_zoom.pdf",energy));
    }

  delete c1;

}



bool getmult(int run, double& bbcs, double& fvtxs, double& fvtxn, double& ebbcs, double& efvtxs, double& efvtxn)
{

  if ( run <= 0 )
    {
      cout << "FATAL: bad run number" << endl;
      exit(-1);
    }

  //TFile* file = TFile::Open(Form("input/hist_%d.root",run));
  //TFile* file = TFile::Open(Form("svrb_run%d_pass0.root",run));
  TFile* file = TFile::Open(Form("RootFiles/svrb_run%d_pass0.root",run));
  if ( !file )
    {
      cout << "missing run " << run << endl;
      return false;
    }

  TH1D* th1d_bbcs = (TH1D*)file->Get("th1d_BBC_charge");
  TH1D* th1d_fvtxs = (TH1D*)file->Get("th1d_FVTXS_nclus");
  TH1D* th1d_fvtxn = (TH1D*)file->Get("th1d_FVTXN_nclus");
  // TH1D* th1d_bbcs = (TH1D*)file->Get("th1d_BBC_charge_IR");
  // TH1D* th1d_fvtxs = (TH1D*)file->Get("th1d_FVTXS_nclus_IR");
  // TH1D* th1d_fvtxn = (TH1D*)file->Get("th1d_FVTXN_nclus_IR");

  if ( !th1d_bbcs || !th1d_fvtxs || !th1d_fvtxn )
    {
      cout << "mising histograms " << run << " " << th1d_bbcs << " " << th1d_fvtxs << " " << th1d_fvtxn << endl;
      return false;
    }

  bbcs = th1d_bbcs->GetMean();
  fvtxs = th1d_fvtxs->GetMean();
  fvtxn = th1d_fvtxn->GetMean();

  ebbcs = th1d_bbcs->GetMeanError();
  efvtxs = th1d_fvtxs->GetMeanError();
  efvtxn = th1d_fvtxn->GetMeanError();

  return true;

}



bool getmult_new(int run, double& bbcs_charge, double& cnt_tracks, double& fvtx_tracks, double& fvtx_clusters)
{

  if ( run <= 0 )
    {
      cout << "FATAL: bad run number" << endl;
      exit(-1);
    }

  //TFile* file = TFile::Open(Form("input/hist_%d.root",run));
  //TFile* file = TFile::Open(Form("svrb_run%d_pass0.root",run));
  TFile* file = TFile::Open(Form("RootFiles/svrb_run%d_pass0.root",run));
  if ( !file )
    {
      cout << "missing run " << run << endl;
      return false;
    }

  TH1D* th1d_cnt = (TH1D*)file->Get("th1d_CNT_ntrk_GC");
  TH1D* th1d_bbcs = (TH1D*)file->Get("th1d_BBC_charge_GC");
  TH1D* th1d_fvtxc = (TH1D*)file->Get("th1d_FVTX_nclus_GC");
  TH1D* th1d_fvtxt = (TH1D*)file->Get("th1d_FVTX_ntrk_GC");

  if ( !th1d_bbcs || !th1d_cnt || !th1d_fvtxc || !th1d_fvtxt )
    {
      cout << "mising histograms " << run << " " << th1d_bbcs << " " << th1d_cnt << " " << th1d_fvtxc << " " << th1d_fvtxt << endl;
      return false;
    }

  cnt_tracks = th1d_cnt->GetMean();
  bbcs_charge = th1d_bbcs->GetMean();
  fvtx_tracks = th1d_fvtxt->GetMean();
  fvtx_clusters = th1d_fvtxc->GetMean();



  return true;

}




void makemult_new(int energy)
{

  gStyle->SetOptTitle(0);

  TCanvas* c1 = new TCanvas("c1","");

  // ------------------------------------------------------------------------------
  // --- these are just first guess, use additional info below to get better values
  double cut_bbc = 0;
  if ( energy == 200 ) cut_bbc = 68.0;
  if ( energy == 62 )  cut_bbc = 46.0;
  if ( energy == 39 )  cut_bbc = 25.0;
  if ( energy == 20 )  cut_bbc = 15.0;
  double cut_cnt = 0;
  if ( energy == 200 ) cut_cnt = 0.70;
  if ( energy == 62 )  cut_cnt = 0.47;
  if ( energy == 39 )  cut_cnt = 0.32;
  if ( energy == 20 )  cut_cnt = 0.18; // also 0.3 from above?
  double cut_fxt = 0;
  if ( energy == 200 ) cut_fxt = 15.0;
  if ( energy == 62 )  cut_fxt = 11.4;
  if ( energy == 39 )  cut_fxt =  6.0;
  if ( energy == 20 )  cut_fxt =  5.0; // also 6.0 from above?
  double cut_fxc = 0;
  if ( energy == 200 ) cut_fxc = 320.0;
  if ( energy == 62 )  cut_fxc = 210.0;
  if ( energy == 39 )  cut_fxc = 210.0; // same as 62 GeV???
  if ( energy == 20 )  cut_fxc = 150.0; // also 200.0 from above?
  // ---
  double cut_bbc_above = 30.0;
  double cut_cnt_above = 0.3;
  double cut_fxt_above = 6.0;
  double cut_fxc_above = 200;

  double array_run[110]; // dumb to make it a double but need it for TGraph...
  double index[110];
  double
    array_cnt[110],
    array_bbc[110],
    array_fxc[110],
    array_fxt[110];
  int run = 0;
  int counter = 0;
  int mcounter = 0;
  double mean_bbc = 0;
  double mean_cnt = 0;
  double mean_fxt = 0;
  double mean_fxc = 0;
  double mean2_bbc = 0;
  double mean2_cnt = 0;
  double mean2_fxt = 0;
  double mean2_fxc = 0;
  ifstream fin((const char*)Form("list_%d.short",energy));
  for ( int i = 0; i < 110; ++i )
    {
      // cout << "energy is " << energy << endl;
      // cout << "list name is " << (const char*)Form("list_%d.short",energy) << endl;
      if ( fin.eof() ) break;
      fin >> run;
      ++counter;
      array_cnt[i] = 0;
      array_bbc[i] = 0;
      array_fxc[i] = 0;
      array_fxt[i] = 0;
      // cout << "run is " << run << endl;
      double temp_bbc, temp_cnt, temp_fxt, temp_fxc;
      index[i] = i+0.5;
      array_run[i] = run;
      if ( !getmult_new(run,temp_bbc,temp_cnt,temp_fxt,temp_fxc) ) continue;
      ++mcounter;
      array_cnt[i] = temp_cnt;
      array_bbc[i] = temp_bbc;
      array_fxc[i] = temp_fxc;
      array_fxt[i] = temp_fxt;
      mean_cnt += temp_cnt;
      mean_bbc += temp_bbc;
      mean_fxt += temp_fxt;
      mean_fxc += temp_fxc;
      mean2_cnt += temp_cnt*temp_cnt;
      mean2_bbc += temp_bbc*temp_bbc;
      mean2_fxt += temp_fxt*temp_fxt;
      mean2_fxc += temp_fxc*temp_fxc;
    }
  fin.close();

  mean_cnt /= mcounter;
  mean_bbc /= mcounter;
  mean_fxt /= mcounter;
  mean_fxc /= mcounter;

  mean2_cnt /= mcounter;
  mean2_bbc /= mcounter;
  mean2_fxt /= mcounter;
  mean2_fxc /= mcounter;

  double vari_bbc = mean2_bbc - mean_bbc*mean_bbc;
  double vari_cnt = mean2_cnt - mean_cnt*mean_cnt;
  double vari_fxt = mean2_fxt - mean_fxt*mean_fxt;
  double vari_fxc = mean2_fxc - mean_fxc*mean_fxc;

  double sigm_bbc = sqrt(vari_bbc);
  double sigm_cnt = sqrt(vari_cnt);
  double sigm_fxt = sqrt(vari_fxt);
  double sigm_fxc = sqrt(vari_fxc);

  cut_bbc = mean_bbc - 3*sigm_bbc;
  cut_cnt = mean_cnt - 3*sigm_cnt;
  cut_fxt = mean_fxt - 3*sigm_fxt;
  cut_fxc = mean_fxc - 3*sigm_fxc;

  cut_bbc_above = mean_bbc + 3*sigm_bbc;
  cut_cnt_above = mean_cnt + 3*sigm_cnt;
  cut_fxt_above = mean_fxt + 3*sigm_fxt;
  cut_fxc_above = mean_fxc + 3*sigm_fxc;

  // ---

  TGraph* tg_run = new TGraph(counter,index,array_run);
  TGraph* tg_cnt = new TGraph(counter,index,array_cnt);
  TGraph* tg_bbc = new TGraph(counter,index,array_bbc);
  TGraph* tg_fxc = new TGraph(counter,index,array_fxc);
  TGraph* tg_fxt = new TGraph(counter,index,array_fxt);

  // ---

  TLine* line = NULL;
  TLine* aline = NULL;
  TLine* mline = NULL;

  tg_cnt->Draw("ap");
  tg_cnt->SetMinimum(0);
  tg_cnt->SetMaximum(1.17*cut_cnt_above);
  if ( line ) delete line;
  line = new TLine(0,cut_cnt,counter,cut_cnt);
  line->SetLineStyle(2);
  line->SetLineWidth(2);
  line->Draw();
  if ( aline ) delete aline;
  aline = new TLine(0,cut_cnt_above,counter,cut_cnt_above);
  aline->SetLineStyle(2);
  aline->SetLineWidth(2);
  aline->Draw();
  if ( mline ) delete mline;
  mline = new TLine(0,mean_cnt,counter,mean_cnt);
  mline->SetLineStyle(2);
  mline->SetLineWidth(2);
  mline->Draw();
  tg_cnt->GetXaxis()->SetLimits(-1,counter+1);
  tg_cnt->SetMarkerStyle(kFullCircle);
  tg_cnt->GetXaxis()->SetTitle("run index");
  tg_cnt->GetYaxis()->SetTitle("cnt tracks per event");
  c1->Print(Form("FigsRun/runindex_cnt_energy%d.png",energy));
  c1->Print(Form("FigsRun/runindex_cnt_energy%d.pdf",energy));

  tg_bbc->Draw("ap");
  tg_bbc->SetMinimum(0);
  tg_bbc->SetMaximum(1.17*cut_bbc_above);
  if ( line ) delete line;
  line = new TLine(0,cut_bbc,counter,cut_bbc);
  line->SetLineStyle(2);
  line->SetLineWidth(2);
  line->Draw();
  line->Draw();
  if ( aline ) delete aline;
  aline = new TLine(0,cut_bbc_above,counter,cut_bbc_above);
  aline->SetLineStyle(2);
  aline->SetLineWidth(2);
  aline->Draw();
  if ( mline ) delete mline;
  mline = new TLine(0,mean_bbc,counter,mean_bbc);
  mline->SetLineStyle(2);
  mline->SetLineWidth(2);
  mline->Draw();
  tg_bbc->GetXaxis()->SetLimits(-1,counter+1);
  tg_bbc->SetMarkerStyle(kFullCircle);
  tg_bbc->GetXaxis()->SetTitle("run index");
  tg_bbc->GetYaxis()->SetTitle("bbc charge");
  c1->Print(Form("FigsRun/runindex_bbc_energy%d.png",energy));
  c1->Print(Form("FigsRun/runindex_bbc_energy%d.pdf",energy));

  tg_fxt->Draw("ap");
  tg_fxt->SetMinimum(0);
  tg_fxt->SetMaximum(1.17*cut_fxt_above);
  if ( line ) delete line;
  line = new TLine(0,cut_fxt,counter,cut_fxt);
  line->SetLineStyle(2);
  line->SetLineWidth(2);
  line->Draw();
  if ( aline ) delete aline;
  aline = new TLine(0,cut_fxt_above,counter,cut_fxt_above);
  aline->SetLineStyle(2);
  aline->SetLineWidth(2);
  aline->Draw();
  if ( mline ) delete mline;
  mline = new TLine(0,mean_fxt,counter,mean_fxt);
  mline->SetLineStyle(2);
  mline->SetLineWidth(2);
  mline->Draw();
  tg_fxt->GetXaxis()->SetLimits(-1,counter+1);
  tg_fxt->SetMarkerStyle(kFullCircle);
  tg_fxt->GetXaxis()->SetTitle("run index");
  tg_fxt->GetYaxis()->SetTitle("FVTX tracks per event");
  c1->Print(Form("FigsRun/runindex_fxt_energy%d.png",energy));
  c1->Print(Form("FigsRun/runindex_fxt_energy%d.pdf",energy));

  tg_fxc->Draw("ap");
  tg_fxc->SetMinimum(0);
  tg_fxc->SetMaximum(1.17*cut_fxc_above);
  if ( line ) delete line;
  line = new TLine(0,cut_fxc,counter,cut_fxc);
  line->SetLineStyle(2);
  line->SetLineWidth(2);
  line->Draw();
  if ( aline ) delete aline;
  aline = new TLine(0,cut_fxc_above,counter,cut_fxc_above);
  aline->SetLineStyle(2);
  aline->SetLineWidth(2);
  aline->Draw();
  if ( mline ) delete mline;
  mline = new TLine(0,mean_fxc,counter,mean_fxc);
  mline->SetLineStyle(2);
  mline->SetLineWidth(2);
  mline->Draw();
  tg_fxc->GetXaxis()->SetLimits(-1,counter+1);
  tg_fxc->SetMarkerStyle(kFullCircle);
  tg_fxc->GetXaxis()->SetTitle("run index");
  tg_fxc->GetYaxis()->SetTitle("fvtx clusters per event");
  c1->Print(Form("FigsRun/runindex_fxc_energy%d.png",energy));
  c1->Print(Form("FigsRun/runindex_fxc_energy%d.pdf",energy));

  tg_run->Draw("ap");
  tg_run->GetXaxis()->SetLimits(-1,counter+1);
  tg_run->SetMarkerStyle(kFullCircle);
  tg_run->GetXaxis()->SetTitle("run index");
  tg_run->GetYaxis()->SetTitle("run number");
  c1->Print(Form("FigsRun/runindex_run_energy%d.png",energy));
  c1->Print(Form("FigsRun/runindex_run_energy%d.pdf",energy));

  // ---

  // ---

  TGraph* tgr_run = new TGraph(counter,array_run,array_run);
  TGraph* tgr_cnt = new TGraph(counter,array_run,array_cnt);
  TGraph* tgr_bbc = new TGraph(counter,array_run,array_bbc);
  TGraph* tgr_fxc = new TGraph(counter,array_run,array_fxc);
  TGraph* tgr_fxt = new TGraph(counter,array_run,array_fxt);

  tgr_cnt->Draw("ap");
  tgr_cnt->SetMinimum(0);
  tgr_cnt->SetMaximum(1.17*cut_cnt_above);
  if ( line ) delete line;
  line = new TLine(0,cut_cnt,counter,cut_cnt);
  line->SetLineStyle(2);
  line->SetLineWidth(2);
  line->Draw();
  if ( aline ) delete aline;
  aline = new TLine(0,cut_cnt_above,counter,cut_cnt_above);
  aline->SetLineStyle(2);
  aline->SetLineWidth(2);
  aline->Draw();
  if ( mline ) delete mline;
  mline = new TLine(0,mean_cnt,counter,mean_cnt);
  mline->SetLineStyle(2);
  mline->SetLineWidth(2);
  mline->Draw();
  tgr_cnt->GetXaxis()->SetLimits(array_run[0]-100,array_run[counter-1]+100);
  tgr_cnt->SetMarkerStyle(kFullCircle);
  tgr_cnt->GetXaxis()->SetTitle("run number");
  tgr_cnt->GetYaxis()->SetTitle("cnt tracks per event");
  c1->Print(Form("FigsRun/run_cnt_energy%d.png",energy));
  c1->Print(Form("FigsRun/run_cnt_energy%d.pdf",energy));

  tgr_bbc->Draw("ap");
  tgr_bbc->SetMinimum(0);
  tgr_bbc->SetMaximum(1.17*cut_bbc_above);
  if ( line ) delete line;
  line = new TLine(0,cut_bbc,counter,cut_bbc);
  line->SetLineStyle(2);
  line->SetLineWidth(2);
  line->Draw();
  line->Draw();
  if ( aline ) delete aline;
  aline = new TLine(0,cut_bbc_above,counter,cut_bbc_above);
  aline->SetLineStyle(2);
  aline->SetLineWidth(2);
  aline->Draw();
  if ( mline ) delete mline;
  mline = new TLine(0,mean_bbc,counter,mean_bbc);
  mline->SetLineStyle(2);
  mline->SetLineWidth(2);
  mline->Draw();
  tgr_bbc->GetXaxis()->SetLimits(array_run[0]-100,array_run[counter-1]+100);
  tgr_bbc->SetMarkerStyle(kFullCircle);
  tgr_bbc->GetXaxis()->SetTitle("run number");
  tgr_bbc->GetYaxis()->SetTitle("bbc charge");
  c1->Print(Form("FigsRun/run_bbc_energy%d.png",energy));
  c1->Print(Form("FigsRun/run_bbc_energy%d.pdf",energy));

  tgr_fxt->Draw("ap");
  tgr_fxt->SetMinimum(0);
  tgr_fxt->SetMaximum(1.17*cut_fxt_above);
  if ( line ) delete line;
  line = new TLine(0,cut_fxt,counter,cut_fxt);
  line->SetLineStyle(2);
  line->SetLineWidth(2);
  line->Draw();
  if ( aline ) delete aline;
  aline = new TLine(0,cut_fxt_above,counter,cut_fxt_above);
  aline->SetLineStyle(2);
  aline->SetLineWidth(2);
  aline->Draw();
  if ( mline ) delete mline;
  mline = new TLine(0,mean_fxt,counter,mean_fxt);
  mline->SetLineStyle(2);
  mline->SetLineWidth(2);
  mline->Draw();
  tgr_fxt->GetXaxis()->SetLimits(array_run[0]-100,array_run[counter-1]+100);
  tgr_fxt->SetMarkerStyle(kFullCircle);
  tgr_fxt->GetXaxis()->SetTitle("run number");
  tgr_fxt->GetYaxis()->SetTitle("FVTX tracks per event");
  c1->Print(Form("FigsRun/run_fxt_energy%d.png",energy));
  c1->Print(Form("FigsRun/run_fxt_energy%d.pdf",energy));

  tgr_fxc->Draw("ap");
  tgr_fxc->SetMinimum(0);
  tgr_fxc->SetMaximum(1.17*cut_fxc_above);
  if ( line ) delete line;
  line = new TLine(0,cut_fxc,counter,cut_fxc);
  line->SetLineStyle(2);
  line->SetLineWidth(2);
  line->Draw();
  if ( aline ) delete aline;
  aline = new TLine(0,cut_fxc_above,counter,cut_fxc_above);
  aline->SetLineStyle(2);
  aline->SetLineWidth(2);
  aline->Draw();
  if ( mline ) delete mline;
  mline = new TLine(0,mean_fxc,counter,mean_fxc);
  mline->SetLineStyle(2);
  mline->SetLineWidth(2);
  mline->Draw();
  tgr_fxc->GetXaxis()->SetLimits(array_run[0]-100,array_run[counter-1]+100);
  tgr_fxc->SetMarkerStyle(kFullCircle);
  tgr_fxc->GetXaxis()->SetTitle("run number");
  tgr_fxc->GetYaxis()->SetTitle("fvtx clusters per event");
  c1->Print(Form("FigsRun/run_fxc_energy%d.png",energy));
  c1->Print(Form("FigsRun/run_fxc_energy%d.pdf",energy));



  if ( line ) delete line;
  if ( aline ) delete aline;
  if ( mline ) delete mline;
  delete c1;



  TFile* fout = TFile::Open(Form("runindex_histograms_energy%d.root",energy),"recreate");
  fout->cd();
  TH1D* th1d_bbc = new TH1D("th1d_bbc","",counter,-0.5,counter-0.5);
  TH1D* th1d_cnt = new TH1D("th1d_cnt","",counter,-0.5,counter-0.5);
  TH1D* th1d_fxt = new TH1D("th1d_fxt","",counter,-0.5,counter-0.5);
  TH1D* th1d_fxc = new TH1D("th1d_fxc","",counter,-0.5,counter-0.5);
  TH1D* th1d_run = new TH1D("th1d_run","",counter,-0.5,counter-0.5);
  for ( int i = 0; i < counter; ++i )
    {
      th1d_bbc->SetBinContent(i+1,array_bbc[i]);
      th1d_cnt->SetBinContent(i+1,array_cnt[i]);
      th1d_fxt->SetBinContent(i+1,array_fxt[i]);
      th1d_fxc->SetBinContent(i+1,array_fxc[i]);
      th1d_run->SetBinContent(i+1,array_run[i]);
    }
  fout->Write();
  fout->Close();


}

