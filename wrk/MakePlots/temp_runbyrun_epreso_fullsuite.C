void getstuff(int, int, double&, double&, double&, double&, double&, double&, double&, double&, double&, double&, double&, double&);


void makeplots(int, int);


void temp_runbyrun_epreso_fullsuite()
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

}


void makeplots(int energy, int harmonic)
{

  gStyle->SetOptTitle(0);

  TCanvas* c1 = new TCanvas("c1","");

  // -----------------------------------------
  // --- first get the individual correlations
  // -----------------------------------------
  double index[100];
  double
    array_BC[100],
    array_BN[100],
    array_BS[100],
    array_CN[100],
    array_CS[100],
    array_NS[100],
    array_eBC[100],
    array_eBN[100],
    array_eBS[100],
    array_eCN[100],
    array_eCS[100],
    array_eNS[100];
  int run;
  int counter = 0;
  ifstream fin((const char*)Form("list_%d.short",energy));
  for ( int i = 0; i < 100; ++i )
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

  double resoB_CN[100]; // BBC(B) resolution using CNT(C) and FVTXN(N)
  double resoB_CS[100]; // BBC(B) resolution using CNT(C) and FVTXS(S)
  double resoB_NS[100]; // BBC(B) resolution using FVTXN(N) and FVTXS(S)

  double resoN_BC[100]; // FVTXN(N) resolution using BBC(B) and CNT(C)
  double resoN_BS[100]; // FVTXN(N) resolution using BBC(B) and FVTXS(S)
  double resoN_CS[100]; // FVTXN(N) resolution using CNT(C) and FVTXS(S)

  double resoS_BC[100]; // FVTXS(S) resolution using BBC(B) and CNT(C)
  double resoS_BN[100]; // FVTXS(S) resolution using BBC(B) and CNT(C)
  double resoS_CN[100]; // FVTXS(S) resolution using CNT(C) and FVTXN(N)

  double eresoB_CN[100]; // BBC(B) resolution using CNT(C) and FVTXN(N)
  double eresoB_CS[100]; // BBC(B) resolution using CNT(C) and FVTXS(S)
  double eresoB_NS[100]; // BBC(B) resolution using FVTXN(N) and FVTXS(S)

  double eresoN_BC[100]; // FVTXN(N) resolution using BBC(B) and CNT(C)
  double eresoN_BS[100]; // FVTXN(N) resolution using BBC(B) and FVTXS(S)
  double eresoN_CS[100]; // FVTXN(N) resolution using CNT(C) and FVTXS(S)

  double eresoS_BC[100]; // FVTXS(S) resolution using BBC(B) and CNT(C)
  double eresoS_BN[100]; // FVTXS(S) resolution using BBC(B) and CNT(C)
  double eresoS_CN[100]; // FVTXS(S) resolution using CNT(C) and FVTXN(N)

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
  tgeresoB_CN->GetXaxis()->SetTitle("Run Index");
  tgeresoB_CN->GetYaxis()->SetTitle("BBC EP resolution");
  tgeresoB_CN->SetTitle("BBC EP resolution");
  tgeresoB_CN->SetMaximum(0.2);
  tgeresoB_CN->SetMinimum(0.0);

  TLegend* legB = new TLegend(0.18,0.68,0.28,0.88);
  legB->AddEntry(tgeresoB_CN,"BBC, CNT, FVTXN","p");
  legB->AddEntry(tgeresoB_CS,"BBC, CNT, FVTXS","p");
  legB->AddEntry(tgeresoB_NS,"BBC, FVTXN, FVTXS","p");
  legB->SetTextSize(0.05);
  legB->Draw();

  c1->Print(Form("figreso_bbc_energy%d_harmonic%d.png",energy,harmonic));
  c1->Print(Form("figreso_bbc_energy%d_harmonic%d.pdf",energy,harmonic));

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
  tgeresoN_BC->GetXaxis()->SetTitle("Run Index");
  tgeresoN_BC->GetYaxis()->SetTitle("FVTXN EP resolution");
  tgeresoN_BC->SetTitle("FVTXN EP resolution");
  tgeresoN_BC->SetMaximum(0.2);
  tgeresoN_BC->SetMinimum(0.0);

  TLegend* legB = new TLegend(0.18,0.68,0.28,0.88);
  legB->AddEntry(tgeresoN_BC,"FVTXN, BBC, CNT","p");
  legB->AddEntry(tgeresoN_BS,"FVTXN, BBC, FVTXS","p");
  legB->AddEntry(tgeresoN_CS,"FVTXN, CNT, FVTXS","p");
  legB->SetTextSize(0.05);
  legB->Draw();

  c1->Print(Form("figreso_fvtxn_energy%d_harmonic%d.png",energy,harmonic));
  c1->Print(Form("figreso_fvtxn_energy%d_harmonic%d.pdf",energy,harmonic));

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
  tgeresoS_BC->GetXaxis()->SetTitle("Run Index");
  tgeresoS_BC->GetYaxis()->SetTitle("FVTXS EP resolution");
  tgeresoS_BC->SetTitle("FVTXS EP resolution");
  tgeresoS_BC->SetMaximum(0.4);
  tgeresoS_BC->SetMinimum(0.0);

  TLegend* legB = new TLegend(0.18,0.68,0.28,0.88);
  legB->AddEntry(tgeresoS_BC,"FVTXS, BBC, CNT","p");
  legB->AddEntry(tgeresoS_BN,"FVTXS, BBC, FVTXN","p");
  legB->AddEntry(tgeresoS_CN,"FVTXS, CNT, FVTXN","p");
  legB->SetTextSize(0.05);
  legB->Draw();

  c1->Print(Form("figreso_fvtxs_energy%d_harmonic%d.png",energy,harmonic));
  c1->Print(Form("figreso_fvtxs_energy%d_harmonic%d.pdf",energy,harmonic));

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

  TProfile* tp1f_BC = (TProfile*)file->Get(Form("tp1f_reso2_BBC_CNT",hh));
  TProfile* tp1f_BN = (TProfile*)file->Get(Form("tp1f_reso%d_BBC_FVTXN",hh));
  TProfile* tp1f_BS = (TProfile*)file->Get(Form("tp1f_reso%d_BBC_FVTX",hh));
  TProfile* tp1f_CN = (TProfile*)file->Get(Form("tp1f_reso2_CNT_FVTXN",hh));
  TProfile* tp1f_CS = (TProfile*)file->Get(Form("tp1f_reso2_CNT_FVTX",hh));
  TProfile* tp1f_NS = (TProfile*)file->Get(Form("tp1f_reso%d_FVTXS_FVTXN",hh));

  // cout << tp1f_BC << endl;
  // cout << tp1f_BC << endl;
  // cout << tp1f_BS << endl;
  // cout << tp1f_CN << endl;
  // cout << tp1f_CS << endl;
  // cout << tp1f_NS << endl;

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

  // cout << "Values assigned..." << endl;
  // cout << double_BC << endl;
  // cout << double_BN << endl;
  // cout << double_BS << endl;
  // cout << double_CN << endl;
  // cout << double_CS << endl;
  // cout << double_NS << endl;

  file->Close();

}

