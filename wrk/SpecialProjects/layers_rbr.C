int counter = 0;

float psi3_offset_fvtxs0[111];
float psi3_offset_fvtxs1[111];
float psi3_offset_fvtxs2[111];
float psi3_offset_fvtxs3[111];

float psi3_magnitude_fvtxs0[111];
float psi3_magnitude_fvtxs1[111];
float psi3_magnitude_fvtxs2[111];
float psi3_magnitude_fvtxs3[111];

float psi3_difference_fvtxs0[111];
float psi3_difference_fvtxs1[111];
float psi3_difference_fvtxs2[111];
float psi3_difference_fvtxs3[111];

float psi3_center_fvtxs0[111];
float psi3_center_fvtxs1[111];
float psi3_center_fvtxs2[111];
float psi3_center_fvtxs3[111];

float e_psi3_offset_fvtxs0[111];
float e_psi3_offset_fvtxs1[111];
float e_psi3_offset_fvtxs2[111];
float e_psi3_offset_fvtxs3[111];

float e_psi3_magnitude_fvtxs0[111];
float e_psi3_magnitude_fvtxs1[111];
float e_psi3_magnitude_fvtxs2[111];
float e_psi3_magnitude_fvtxs3[111];

float psi2_offset_fvtxs0[111];
float psi2_offset_fvtxs1[111];
float psi2_offset_fvtxs2[111];
float psi2_offset_fvtxs3[111];

float psi2_magnitude_fvtxs0[111];
float psi2_magnitude_fvtxs1[111];
float psi2_magnitude_fvtxs2[111];
float psi2_magnitude_fvtxs3[111];

float psi2_difference_fvtxs0[111];
float psi2_difference_fvtxs1[111];
float psi2_difference_fvtxs2[111];
float psi2_difference_fvtxs3[111];

float psi2_center_fvtxs0[111];
float psi2_center_fvtxs1[111];
float psi2_center_fvtxs2[111];
float psi2_center_fvtxs3[111];

float e_psi2_offset_fvtxs0[111];
float e_psi2_offset_fvtxs1[111];
float e_psi2_offset_fvtxs2[111];
float e_psi2_offset_fvtxs3[111];

float e_psi2_magnitude_fvtxs0[111];
float e_psi2_magnitude_fvtxs1[111];
float e_psi2_magnitude_fvtxs2[111];
float e_psi2_magnitude_fvtxs3[111];

float runnumber[111]; // actually an int but convenient to have a float
float runindex[111];

void doit(int);

void layers_rbr()
{

  gStyle->SetOptFit(0);
  gStyle->SetOptTitle(1);

  int run;
  ifstream fin;

  counter = 0;
  fin.open("list200.txt");
  while ( fin >> run )
    {
      doit(run);
      runnumber[counter] = run;
      runindex[counter] = counter;
      ++counter;
    }
  fin.close();

  TCanvas* c1 = new TCanvas();

  for ( int i = 0; i < counter; ++i )
    {
      cout << "LOOK HERE!!! " << i << " " << runindex[i] << " " << psi3_difference_fvtxs0[i] << endl;
    }



  TGraph* tg_psi3_index_difference_fvtxs0 = new TGraph(counter,runindex,psi3_difference_fvtxs0);
  tg_psi3_index_difference_fvtxs0->SetMarkerStyle(kFullCircle);
  tg_psi3_index_difference_fvtxs0->SetMarkerColor(kBlack);
  tg_psi3_index_difference_fvtxs0->SetTitle("FVTXS layer 0");
  tg_psi3_index_difference_fvtxs0->Draw("ap");
  tg_psi3_index_difference_fvtxs0->SetMaximum(1.0);
  tg_psi3_index_difference_fvtxs0->SetMinimum(0.0);
  tg_psi3_index_difference_fvtxs0->GetXaxis()->SetLimits(-1,counter+1);
  tg_psi3_index_difference_fvtxs0->GetXaxis()->SetTitle("run index");
  tg_psi3_index_difference_fvtxs0->GetYaxis()->SetTitle("relative difference");
  c1->Print("Jamie/Figs/final_psi3_difference_fvtxs0.png");
  c1->Print("Jamie/Figs/final_psi3_difference_fvtxs0.pdf");

  TGraph* tg_psi3_index_center_fvtxs0 = new TGraph(counter,runindex,psi3_center_fvtxs0);
  tg_psi3_index_center_fvtxs0->SetMarkerStyle(kFullCircle);
  tg_psi3_index_center_fvtxs0->SetMarkerColor(kBlack);
  tg_psi3_index_center_fvtxs0->SetTitle("FVTXS layer 0");
  tg_psi3_index_center_fvtxs0->Draw("ap");
  tg_psi3_index_center_fvtxs0->SetMaximum(1.0);
  tg_psi3_index_center_fvtxs0->SetMinimum(-1.0);
  tg_psi3_index_center_fvtxs0->GetXaxis()->SetLimits(-1,counter+1);
  tg_psi3_index_center_fvtxs0->GetXaxis()->SetTitle("run index");
  tg_psi3_index_center_fvtxs0->GetYaxis()->SetTitle("center");
  c1->Print("Jamie/Figs/final_psi3_center_fvtxs0.png");
  c1->Print("Jamie/Figs/final_psi3_center_fvtxs0.pdf");


  TGraph* tg_psi3_index_difference_fvtxs1 = new TGraph(counter,runindex,psi3_difference_fvtxs1);
  tg_psi3_index_difference_fvtxs1->SetMarkerStyle(kFullCircle);
  tg_psi3_index_difference_fvtxs1->SetMarkerColor(kBlack);
  tg_psi3_index_difference_fvtxs1->SetTitle("FVTXS layer 1");
  tg_psi3_index_difference_fvtxs1->Draw("ap");
  tg_psi3_index_difference_fvtxs1->SetMaximum(1.0);
  tg_psi3_index_difference_fvtxs1->SetMinimum(0.0);
  tg_psi3_index_difference_fvtxs1->GetXaxis()->SetLimits(-1,counter+1);
  tg_psi3_index_difference_fvtxs1->GetXaxis()->SetTitle("run index");
  tg_psi3_index_difference_fvtxs1->GetYaxis()->SetTitle("relative difference");
  c1->Print("Jamie/Figs/final_psi3_difference_fvtxs1.png");
  c1->Print("Jamie/Figs/final_psi3_difference_fvtxs1.pdf");

  TGraph* tg_psi3_index_center_fvtxs1 = new TGraph(counter,runindex,psi3_center_fvtxs1);
  tg_psi3_index_center_fvtxs1->SetMarkerStyle(kFullCircle);
  tg_psi3_index_center_fvtxs1->SetMarkerColor(kBlack);
  tg_psi3_index_center_fvtxs1->SetTitle("FVTXS layer 1");
  tg_psi3_index_center_fvtxs1->Draw("ap");
  tg_psi3_index_center_fvtxs1->SetMaximum(1.0);
  tg_psi3_index_center_fvtxs1->SetMinimum(-1.0);
  tg_psi3_index_center_fvtxs1->GetXaxis()->SetLimits(-1,counter+1);
  tg_psi3_index_center_fvtxs1->GetXaxis()->SetTitle("run index");
  tg_psi3_index_center_fvtxs1->GetYaxis()->SetTitle("center");
  c1->Print("Jamie/Figs/final_psi3_center_fvtxs1.png");
  c1->Print("Jamie/Figs/final_psi3_center_fvtxs1.pdf");


  TGraph* tg_psi3_index_difference_fvtxs2 = new TGraph(counter,runindex,psi3_difference_fvtxs2);
  tg_psi3_index_difference_fvtxs2->SetMarkerStyle(kFullCircle);
  tg_psi3_index_difference_fvtxs2->SetMarkerColor(kBlack);
  tg_psi3_index_difference_fvtxs2->SetTitle("FVTXS layer 2");
  tg_psi3_index_difference_fvtxs2->Draw("ap");
  tg_psi3_index_difference_fvtxs2->SetMaximum(1.0);
  tg_psi3_index_difference_fvtxs2->SetMinimum(0.0);
  tg_psi3_index_difference_fvtxs2->GetXaxis()->SetLimits(-1,counter+1);
  tg_psi3_index_difference_fvtxs2->GetXaxis()->SetTitle("run index");
  tg_psi3_index_difference_fvtxs2->GetYaxis()->SetTitle("relative difference");
  c1->Print("Jamie/Figs/final_psi3_difference_fvtxs2.png");
  c1->Print("Jamie/Figs/final_psi3_difference_fvtxs2.pdf");

  TGraph* tg_psi3_index_center_fvtxs2 = new TGraph(counter,runindex,psi3_center_fvtxs2);
  tg_psi3_index_center_fvtxs2->SetMarkerStyle(kFullCircle);
  tg_psi3_index_center_fvtxs2->SetMarkerColor(kBlack);
  tg_psi3_index_center_fvtxs2->SetTitle("FVTXS layer 2");
  tg_psi3_index_center_fvtxs2->Draw("ap");
  tg_psi3_index_center_fvtxs2->SetMaximum(1.0);
  tg_psi3_index_center_fvtxs2->SetMinimum(-1.0);
  tg_psi3_index_center_fvtxs2->GetXaxis()->SetLimits(-1,counter+1);
  tg_psi3_index_center_fvtxs2->GetXaxis()->SetTitle("run index");
  tg_psi3_index_center_fvtxs2->GetYaxis()->SetTitle("center");
  c1->Print("Jamie/Figs/final_psi3_center_fvtxs2.png");
  c1->Print("Jamie/Figs/final_psi3_center_fvtxs2.pdf");


  TGraph* tg_psi3_index_difference_fvtxs3 = new TGraph(counter,runindex,psi3_difference_fvtxs3);
  tg_psi3_index_difference_fvtxs3->SetMarkerStyle(kFullCircle);
  tg_psi3_index_difference_fvtxs3->SetMarkerColor(kBlack);
  tg_psi3_index_difference_fvtxs3->SetTitle("FVTXS layer 3");
  tg_psi3_index_difference_fvtxs3->Draw("ap");
  tg_psi3_index_difference_fvtxs3->SetMaximum(1.0);
  tg_psi3_index_difference_fvtxs3->SetMinimum(0.0);
  tg_psi3_index_difference_fvtxs3->GetXaxis()->SetLimits(-1,counter+1);
  tg_psi3_index_difference_fvtxs3->GetXaxis()->SetTitle("run index");
  tg_psi3_index_difference_fvtxs3->GetYaxis()->SetTitle("relative difference");
  c1->Print("Jamie/Figs/final_psi3_difference_fvtxs3.png");
  c1->Print("Jamie/Figs/final_psi3_difference_fvtxs3.pdf");

  TGraph* tg_psi3_index_center_fvtxs3 = new TGraph(counter,runindex,psi3_center_fvtxs3);
  tg_psi3_index_center_fvtxs3->SetMarkerStyle(kFullCircle);
  tg_psi3_index_center_fvtxs3->SetMarkerColor(kBlack);
  tg_psi3_index_center_fvtxs3->SetTitle("FVTXS layer 3");
  tg_psi3_index_center_fvtxs3->Draw("ap");
  tg_psi3_index_center_fvtxs3->SetMaximum(1.0);
  tg_psi3_index_center_fvtxs3->SetMinimum(-1.0);
  tg_psi3_index_center_fvtxs3->GetXaxis()->SetLimits(-1,counter+1);
  tg_psi3_index_center_fvtxs3->GetXaxis()->SetTitle("run index");
  tg_psi3_index_center_fvtxs3->GetYaxis()->SetTitle("center");
  c1->Print("Jamie/Figs/final_psi3_center_fvtxs3.png");
  c1->Print("Jamie/Figs/final_psi3_center_fvtxs3.pdf");



  // ---
  // --- come back here for psi2
  // ---

  TGraph* tg_psi2_index_difference_fvtxs0 = new TGraph(counter,runindex,psi2_difference_fvtxs0);
  tg_psi2_index_difference_fvtxs0->SetMarkerStyle(kFullCircle);
  tg_psi2_index_difference_fvtxs0->SetMarkerColor(kBlack);
  tg_psi2_index_difference_fvtxs0->SetTitle("FVTXS layer 0");
  tg_psi2_index_difference_fvtxs0->Draw("ap");
  tg_psi2_index_difference_fvtxs0->SetMaximum(1.0);
  tg_psi2_index_difference_fvtxs0->SetMinimum(0.0);
  tg_psi2_index_difference_fvtxs0->GetXaxis()->SetLimits(-1,counter+1);
  tg_psi2_index_difference_fvtxs0->GetXaxis()->SetTitle("run index");
  tg_psi2_index_difference_fvtxs0->GetYaxis()->SetTitle("relative difference");
  c1->Print("Jamie/Figs/final_psi2_difference_fvtxs0.png");
  c1->Print("Jamie/Figs/final_psi2_difference_fvtxs0.pdf");

  TGraph* tg_psi2_index_center_fvtxs0 = new TGraph(counter,runindex,psi2_center_fvtxs0);
  tg_psi2_index_center_fvtxs0->SetMarkerStyle(kFullCircle);
  tg_psi2_index_center_fvtxs0->SetMarkerColor(kBlack);
  tg_psi2_index_center_fvtxs0->SetTitle("FVTXS layer 0");
  tg_psi2_index_center_fvtxs0->Draw("ap");
  tg_psi2_index_center_fvtxs0->SetMaximum(1.0);
  tg_psi2_index_center_fvtxs0->SetMinimum(-1.0);
  tg_psi2_index_center_fvtxs0->GetXaxis()->SetLimits(-1,counter+1);
  tg_psi2_index_center_fvtxs0->GetXaxis()->SetTitle("run index");
  tg_psi2_index_center_fvtxs0->GetYaxis()->SetTitle("center");
  c1->Print("Jamie/Figs/final_psi2_center_fvtxs0.png");
  c1->Print("Jamie/Figs/final_psi2_center_fvtxs0.pdf");


  TGraph* tg_psi2_index_difference_fvtxs1 = new TGraph(counter,runindex,psi2_difference_fvtxs1);
  tg_psi2_index_difference_fvtxs1->SetMarkerStyle(kFullCircle);
  tg_psi2_index_difference_fvtxs1->SetMarkerColor(kBlack);
  tg_psi2_index_difference_fvtxs1->SetTitle("FVTXS layer 1");
  tg_psi2_index_difference_fvtxs1->Draw("ap");
  tg_psi2_index_difference_fvtxs1->SetMaximum(1.0);
  tg_psi2_index_difference_fvtxs1->SetMinimum(0.0);
  tg_psi2_index_difference_fvtxs1->GetXaxis()->SetLimits(-1,counter+1);
  tg_psi2_index_difference_fvtxs1->GetXaxis()->SetTitle("run index");
  tg_psi2_index_difference_fvtxs1->GetYaxis()->SetTitle("relative difference");
  c1->Print("Jamie/Figs/final_psi2_difference_fvtxs1.png");
  c1->Print("Jamie/Figs/final_psi2_difference_fvtxs1.pdf");

  TGraph* tg_psi2_index_center_fvtxs1 = new TGraph(counter,runindex,psi2_center_fvtxs1);
  tg_psi2_index_center_fvtxs1->SetMarkerStyle(kFullCircle);
  tg_psi2_index_center_fvtxs1->SetMarkerColor(kBlack);
  tg_psi2_index_center_fvtxs1->SetTitle("FVTXS layer 1");
  tg_psi2_index_center_fvtxs1->Draw("ap");
  tg_psi2_index_center_fvtxs1->SetMaximum(1.0);
  tg_psi2_index_center_fvtxs1->SetMinimum(-1.0);
  tg_psi2_index_center_fvtxs1->GetXaxis()->SetLimits(-1,counter+1);
  tg_psi2_index_center_fvtxs1->GetXaxis()->SetTitle("run index");
  tg_psi2_index_center_fvtxs1->GetYaxis()->SetTitle("center");
  c1->Print("Jamie/Figs/final_psi2_center_fvtxs1.png");
  c1->Print("Jamie/Figs/final_psi2_center_fvtxs1.pdf");


  TGraph* tg_psi2_index_difference_fvtxs2 = new TGraph(counter,runindex,psi2_difference_fvtxs2);
  tg_psi2_index_difference_fvtxs2->SetMarkerStyle(kFullCircle);
  tg_psi2_index_difference_fvtxs2->SetMarkerColor(kBlack);
  tg_psi2_index_difference_fvtxs2->SetTitle("FVTXS layer 2");
  tg_psi2_index_difference_fvtxs2->Draw("ap");
  tg_psi2_index_difference_fvtxs2->SetMaximum(1.0);
  tg_psi2_index_difference_fvtxs2->SetMinimum(0.0);
  tg_psi2_index_difference_fvtxs2->GetXaxis()->SetLimits(-1,counter+1);
  tg_psi2_index_difference_fvtxs2->GetXaxis()->SetTitle("run index");
  tg_psi2_index_difference_fvtxs2->GetYaxis()->SetTitle("relative difference");
  c1->Print("Jamie/Figs/final_psi2_difference_fvtxs2.png");
  c1->Print("Jamie/Figs/final_psi2_difference_fvtxs2.pdf");

  TGraph* tg_psi2_index_center_fvtxs2 = new TGraph(counter,runindex,psi2_center_fvtxs2);
  tg_psi2_index_center_fvtxs2->SetMarkerStyle(kFullCircle);
  tg_psi2_index_center_fvtxs2->SetMarkerColor(kBlack);
  tg_psi2_index_center_fvtxs2->SetTitle("FVTXS layer 2");
  tg_psi2_index_center_fvtxs2->Draw("ap");
  tg_psi2_index_center_fvtxs2->SetMaximum(1.0);
  tg_psi2_index_center_fvtxs2->SetMinimum(-1.0);
  tg_psi2_index_center_fvtxs2->GetXaxis()->SetLimits(-1,counter+1);
  tg_psi2_index_center_fvtxs2->GetXaxis()->SetTitle("run index");
  tg_psi2_index_center_fvtxs2->GetYaxis()->SetTitle("center");
  c1->Print("Jamie/Figs/final_psi2_center_fvtxs2.png");
  c1->Print("Jamie/Figs/final_psi2_center_fvtxs2.pdf");


  TGraph* tg_psi2_index_difference_fvtxs3 = new TGraph(counter,runindex,psi2_difference_fvtxs3);
  tg_psi2_index_difference_fvtxs3->SetMarkerStyle(kFullCircle);
  tg_psi2_index_difference_fvtxs3->SetMarkerColor(kBlack);
  tg_psi2_index_difference_fvtxs3->SetTitle("FVTXS layer 3");
  tg_psi2_index_difference_fvtxs3->Draw("ap");
  tg_psi2_index_difference_fvtxs3->SetMaximum(1.0);
  tg_psi2_index_difference_fvtxs3->SetMinimum(0.0);
  tg_psi2_index_difference_fvtxs3->GetXaxis()->SetLimits(-1,counter+1);
  tg_psi2_index_difference_fvtxs3->GetXaxis()->SetTitle("run index");
  tg_psi2_index_difference_fvtxs3->GetYaxis()->SetTitle("relative difference");
  c1->Print("Jamie/Figs/final_psi2_difference_fvtxs3.png");
  c1->Print("Jamie/Figs/final_psi2_difference_fvtxs3.pdf");

  TGraph* tg_psi2_index_center_fvtxs3 = new TGraph(counter,runindex,psi2_center_fvtxs3);
  tg_psi2_index_center_fvtxs3->SetMarkerStyle(kFullCircle);
  tg_psi2_index_center_fvtxs3->SetMarkerColor(kBlack);
  tg_psi2_index_center_fvtxs3->SetTitle("FVTXS layer 3");
  tg_psi2_index_center_fvtxs3->Draw("ap");
  tg_psi2_index_center_fvtxs3->SetMaximum(1.0);
  tg_psi2_index_center_fvtxs3->SetMinimum(-1.0);
  tg_psi2_index_center_fvtxs3->GetXaxis()->SetLimits(-1,counter+1);
  tg_psi2_index_center_fvtxs3->GetXaxis()->SetTitle("run index");
  tg_psi2_index_center_fvtxs3->GetYaxis()->SetTitle("center");
  c1->Print("Jamie/Figs/final_psi2_center_fvtxs3.png");
  c1->Print("Jamie/Figs/final_psi2_center_fvtxs3.pdf");




}


void doit(int run)
{

  //TFile* file = TFile::Open(Form("RootFiles/hist_%d.root",run));
  TFile* file = TFile::Open(Form("../output/files_200_mostlyprevious/hist_%d.root",run));
  if ( !file )
    {
      cout << "WARNING: file does not exist for run " << run << endl;
      return;
    }
  TList* list = (TList*)file->GetListOfKeys();
  int size = list->GetSize();
  if ( size < 1 )
    {
      cout << "WARNING: file does not have enough keys, run " << run << endl;
      return;
    }

  TCanvas* c1 = new TCanvas();

  TH2D* th2d_psi3_fvtxs0 = (TH2D*)file->Get("psi_bf_0_2_4");
  TH2D* th2d_psi3_fvtxs1 = (TH2D*)file->Get("psi_bf_0_2_5");
  TH2D* th2d_psi3_fvtxs2 = (TH2D*)file->Get("psi_bf_0_2_6");
  TH2D* th2d_psi3_fvtxs3 = (TH2D*)file->Get("psi_bf_0_2_7");
  TH1D* th1d_psi3_fvtxs0 = (TH1D*)th2d_psi3_fvtxs0->ProjectionY("th1d_psi3_fvtxs0");
  TH1D* th1d_psi3_fvtxs1 = (TH1D*)th2d_psi3_fvtxs1->ProjectionY("th1d_psi3_fvtxs1");
  TH1D* th1d_psi3_fvtxs2 = (TH1D*)th2d_psi3_fvtxs2->ProjectionY("th1d_psi3_fvtxs2");
  TH1D* th1d_psi3_fvtxs3 = (TH1D*)th2d_psi3_fvtxs3->ProjectionY("th1d_psi3_fvtxs3");

  TF1* fun = new TF1("fun","[0]+[1]*cos(3*(x-[2]))",-1.0,1.0);
  fun->SetParLimits(2,-1.0,1.0);
  float par0, par1, par2;
  float epar0, epar1, epar2;
  float max;
  float min;
  float cen;
  float off;
  float mag;
  TLine* line_min = NULL;
  TLine* line_max = NULL;
  TLine* line_cen = NULL;

  max = th1d_psi3_fvtxs0->GetMaximum();
  cen = th1d_psi3_fvtxs0->GetBinCenter(th1d_psi3_fvtxs0->GetMaximumBin());
  min = th1d_psi3_fvtxs0->GetMinimum(max*0.1); // minimum bin with at least 1 entry
  fun->SetParameter(0,max);
  fun->SetParameter(1,max);
  fun->SetParameter(2,0.0);
  th1d_psi3_fvtxs0->SetTitle("FVTXS layer 0");
  th1d_psi3_fvtxs0->Draw();
  th1d_psi3_fvtxs0->Fit(fun,"","",-1.0,1.0);
  if ( line_max ) delete line_max;
  line_max = new TLine(-4,max,4,max);
  line_max->SetLineStyle(2);
  line_max->SetLineWidth(2);
  line_max->Draw();
  if ( line_min ) delete line_min;
  line_min = new TLine(-4,min,4,min);
  line_min->SetLineStyle(2);
  line_min->SetLineWidth(2);
  line_min->Draw();
  if ( line_cen ) delete line_cen;
  line_cen = new TLine(cen,0,cen,max);
  line_cen->SetLineStyle(2);
  line_cen->SetLineWidth(2);
  line_cen->Draw();
  c1->Print(Form("Jamie/Figs/psi3fit_fvtxs0_run%d.png",run));
  c1->Print(Form("Jamie/Figs/psi3fit_fvtxs0_run%d.pdf",run));
  par0 = fun->GetParameter(0);
  par1 = fun->GetParameter(1);
  par2 = fun->GetParameter(2);
  epar0 = fun->GetParError(0);
  epar1 = fun->GetParError(1);
  epar2 = fun->GetParError(2);
  if ( par0 > 0 ) mag = par1/par0;
  off = par2;
  psi3_offset_fvtxs0[counter] = off;
  psi3_magnitude_fvtxs0[counter] = mag;
  e_psi3_offset_fvtxs0[counter] = epar2;
  e_psi3_magnitude_fvtxs0[counter] = (epar1/par1)*mag;
  psi3_difference_fvtxs0[counter] = (max-min)/(max);
  psi3_center_fvtxs0[counter] = cen;

  max = th1d_psi3_fvtxs1->GetMaximum();
  cen = th1d_psi3_fvtxs1->GetBinCenter(th1d_psi3_fvtxs1->GetMaximumBin());
  min = th1d_psi3_fvtxs1->GetMinimum(max*0.1); // minimum bin with at least 1 entry
  fun->SetParameter(0,max);
  fun->SetParameter(1,max);
  fun->SetParameter(2,0.0);
  th1d_psi3_fvtxs1->SetTitle("FVTXS layer 1");
  th1d_psi3_fvtxs1->Draw();
  th1d_psi3_fvtxs1->Fit(fun,"","",-1.0,1.0);
  if ( line_max ) delete line_max;
  line_max = new TLine(-4,max,4,max);
  line_max->SetLineStyle(2);
  line_max->SetLineWidth(2);
  line_max->Draw();
  if ( line_min ) delete line_min;
  line_min = new TLine(-4,min,4,min);
  line_min->SetLineStyle(2);
  line_min->SetLineWidth(2);
  line_min->Draw();
  if ( line_cen ) delete line_cen;
  line_cen = new TLine(cen,0,cen,max);
  line_cen->SetLineStyle(2);
  line_cen->SetLineWidth(2);
  line_cen->Draw();
  c1->Print(Form("Jamie/Figs/psi3fit_fvtxs1_run%d.png",run));
  c1->Print(Form("Jamie/Figs/psi3fit_fvtxs1_run%d.pdf",run));
  par0 = fun->GetParameter(0);
  par1 = fun->GetParameter(1);
  par2 = fun->GetParameter(2);
  epar0 = fun->GetParError(0);
  epar1 = fun->GetParError(1);
  epar2 = fun->GetParError(2);
  if ( par0 > 0 ) mag = par1/par0;
  off = par2;
  psi3_offset_fvtxs1[counter] = off;
  psi3_magnitude_fvtxs1[counter] = mag;
  e_psi3_offset_fvtxs1[counter] = epar2;
  e_psi3_magnitude_fvtxs1[counter] = (epar1/par1)*mag;
  psi3_difference_fvtxs1[counter] = (max-min)/(max);
  psi3_center_fvtxs1[counter] = cen;

  max = th1d_psi3_fvtxs2->GetMaximum();
  cen = th1d_psi3_fvtxs2->GetBinCenter(th1d_psi3_fvtxs2->GetMaximumBin());
  min = th1d_psi3_fvtxs2->GetMinimum(max*0.1); // minimum bin with at least 1 entry
  fun->SetParameter(0,max);
  fun->SetParameter(1,max);
  fun->SetParameter(2,0.0);
  th1d_psi3_fvtxs2->SetTitle("FVTXS layer 2");
  th1d_psi3_fvtxs2->Draw();
  th1d_psi3_fvtxs2->Fit(fun,"","",-1.0,1.0);
  if ( line_max ) delete line_max;
  line_max = new TLine(-4,max,4,max);
  line_max->SetLineStyle(2);
  line_max->SetLineWidth(2);
  line_max->Draw();
  if ( line_min ) delete line_min;
  line_min = new TLine(-4,min,4,min);
  line_min->SetLineStyle(2);
  line_min->SetLineWidth(2);
  line_min->Draw();
  if ( line_cen ) delete line_cen;
  line_cen = new TLine(cen,0,cen,max);
  line_cen->SetLineStyle(2);
  line_cen->SetLineWidth(2);
  line_cen->Draw();
  c1->Print(Form("Jamie/Figs/psi3fit_fvtxs2_run%d.png",run));
  c1->Print(Form("Jamie/Figs/psi3fit_fvtxs2_run%d.pdf",run));
  par0 = fun->GetParameter(0);
  par1 = fun->GetParameter(1);
  par2 = fun->GetParameter(2);
  epar0 = fun->GetParError(0);
  epar1 = fun->GetParError(1);
  epar2 = fun->GetParError(2);
  if ( par0 > 0 ) mag = par1/par0;
  off = par2;
  psi3_offset_fvtxs2[counter] = off;
  psi3_magnitude_fvtxs2[counter] = mag;
  e_psi3_offset_fvtxs2[counter] = epar2;
  e_psi3_magnitude_fvtxs2[counter] = (epar1/par1)*mag;
  psi3_difference_fvtxs2[counter] = (max-min)/(max);
  psi3_center_fvtxs2[counter] = cen;

  max = th1d_psi3_fvtxs3->GetMaximum();
  cen = th1d_psi3_fvtxs3->GetBinCenter(th1d_psi3_fvtxs3->GetMaximumBin());
  min = th1d_psi3_fvtxs3->GetMinimum(max*0.1); // minimum bin with at least 1 entry
  fun->SetParameter(0,max);
  fun->SetParameter(1,max);
  fun->SetParameter(2,0.0);
  th1d_psi3_fvtxs3->SetTitle("FVTXS layer 3");
  th1d_psi3_fvtxs3->Draw();
  th1d_psi3_fvtxs3->Fit(fun,"","",-1.0,1.0);
  if ( line_max ) delete line_max;
  line_max = new TLine(-4,max,4,max);
  line_max->SetLineStyle(2);
  line_max->SetLineWidth(2);
  line_max->Draw();
  if ( line_min ) delete line_min;
  line_min = new TLine(-4,min,4,min);
  line_min->SetLineStyle(2);
  line_min->SetLineWidth(2);
  line_min->Draw();
  if ( line_cen ) delete line_cen;
  line_cen = new TLine(cen,0,cen,max);
  line_cen->SetLineStyle(2);
  line_cen->SetLineWidth(2);
  line_cen->Draw();
  c1->Print(Form("Jamie/Figs/psi3fit_fvtxs3_run%d.png",run));
  c1->Print(Form("Jamie/Figs/psi3fit_fvtxs3_run%d.pdf",run));
  par0 = fun->GetParameter(0);
  par1 = fun->GetParameter(1);
  par2 = fun->GetParameter(2);
  epar0 = fun->GetParError(0);
  epar1 = fun->GetParError(1);
  epar2 = fun->GetParError(2);
  if ( par0 > 0 ) mag = par1/par0;
  off = par2;
  psi3_offset_fvtxs3[counter] = off;
  psi3_magnitude_fvtxs3[counter] = mag;
  e_psi3_offset_fvtxs3[counter] = epar2;
  e_psi3_magnitude_fvtxs3[counter] = (epar1/par1)*mag;
  psi3_difference_fvtxs3[counter] = (max-min)/(max);
  psi3_center_fvtxs3[counter] = cen;



  // ---
  // --- come back here for psi2
  // ---

  TH2D* th2d_psi2_fvtxs0 = (TH2D*)file->Get("psi_bf_0_1_4");
  TH2D* th2d_psi2_fvtxs1 = (TH2D*)file->Get("psi_bf_0_1_5");
  TH2D* th2d_psi2_fvtxs2 = (TH2D*)file->Get("psi_bf_0_1_6");
  TH2D* th2d_psi2_fvtxs3 = (TH2D*)file->Get("psi_bf_0_1_7");
  TH1D* th1d_psi2_fvtxs0 = (TH1D*)th2d_psi2_fvtxs0->ProjectionY("th1d_psi2_fvtxs0");
  TH1D* th1d_psi2_fvtxs1 = (TH1D*)th2d_psi2_fvtxs1->ProjectionY("th1d_psi2_fvtxs1");
  TH1D* th1d_psi2_fvtxs2 = (TH1D*)th2d_psi2_fvtxs2->ProjectionY("th1d_psi2_fvtxs2");
  TH1D* th1d_psi2_fvtxs3 = (TH1D*)th2d_psi2_fvtxs3->ProjectionY("th1d_psi2_fvtxs3");

  if ( fun ) delete fun;
  fun = new TF1("fun","[0]+[1]*cos(2*(x-[2]))",-1.6,1.6);
  fun->SetParLimits(2,-1.57,1.57);

  max = th1d_psi2_fvtxs0->GetMaximum();
  cen = th1d_psi2_fvtxs0->GetBinCenter(th1d_psi2_fvtxs0->GetMaximumBin());
  min = th1d_psi2_fvtxs0->GetMinimum(max*0.15); // minimum bin with at least 1 entry
  fun->SetParameter(0,max);
  fun->SetParameter(1,max);
  fun->SetParameter(2,0.0);
  th1d_psi2_fvtxs0->SetTitle("FVTXS layer 0");
  th1d_psi2_fvtxs0->Draw();
  th1d_psi2_fvtxs0->Fit(fun,"","",-1.55,1.55);
  if ( line_max ) delete line_max;
  line_max = new TLine(-4,max,4,max);
  line_max->SetLineStyle(2);
  line_max->SetLineWidth(2);
  line_max->Draw();
  if ( line_min ) delete line_min;
  line_min = new TLine(-4,min,4,min);
  line_min->SetLineStyle(2);
  line_min->SetLineWidth(2);
  line_min->Draw();
  if ( line_cen ) delete line_cen;
  line_cen = new TLine(cen,0,cen,max);
  line_cen->SetLineStyle(2);
  line_cen->SetLineWidth(2);
  line_cen->Draw();
  c1->Print(Form("Jamie/Figs/psi2fit_fvtxs0_run%d.png",run));
  c1->Print(Form("Jamie/Figs/psi2fit_fvtxs0_run%d.pdf",run));
  par0 = fun->GetParameter(0);
  par1 = fun->GetParameter(1);
  par2 = fun->GetParameter(2);
  epar0 = fun->GetParError(0);
  epar1 = fun->GetParError(1);
  epar2 = fun->GetParError(2);
  if ( par0 > 0 ) mag = par1/par0;
  off = par2;
  psi2_offset_fvtxs0[counter] = off;
  psi2_magnitude_fvtxs0[counter] = mag;
  e_psi2_offset_fvtxs0[counter] = epar2;
  e_psi2_magnitude_fvtxs0[counter] = (epar1/par1)*mag;
  psi2_difference_fvtxs0[counter] = (max-min)/(max);
  psi2_center_fvtxs0[counter] = cen;

  max = th1d_psi2_fvtxs1->GetMaximum();
  cen = th1d_psi2_fvtxs1->GetBinCenter(th1d_psi2_fvtxs1->GetMaximumBin());
  min = th1d_psi2_fvtxs1->GetMinimum(max*0.15); // minimum bin with at least 1 entry
  fun->SetParameter(0,max);
  fun->SetParameter(1,max);
  fun->SetParameter(2,0.0);
  th1d_psi2_fvtxs1->SetTitle("FVTXS layer 1");
  th1d_psi2_fvtxs1->Draw();
  th1d_psi2_fvtxs1->Fit(fun,"","",-1.55,1.55);
  if ( line_max ) delete line_max;
  line_max = new TLine(-4,max,4,max);
  line_max->SetLineStyle(2);
  line_max->SetLineWidth(2);
  line_max->Draw();
  if ( line_min ) delete line_min;
  line_min = new TLine(-4,min,4,min);
  line_min->SetLineStyle(2);
  line_min->SetLineWidth(2);
  line_min->Draw();
  if ( line_cen ) delete line_cen;
  line_cen = new TLine(cen,0,cen,max);
  line_cen->SetLineStyle(2);
  line_cen->SetLineWidth(2);
  line_cen->Draw();
  c1->Print(Form("Jamie/Figs/psi2fit_fvtxs1_run%d.png",run));
  c1->Print(Form("Jamie/Figs/psi2fit_fvtxs1_run%d.pdf",run));
  par0 = fun->GetParameter(0);
  par1 = fun->GetParameter(1);
  par2 = fun->GetParameter(2);
  epar0 = fun->GetParError(0);
  epar1 = fun->GetParError(1);
  epar2 = fun->GetParError(2);
  if ( par0 > 0 ) mag = par1/par0;
  off = par2;
  psi2_offset_fvtxs1[counter] = off;
  psi2_magnitude_fvtxs1[counter] = mag;
  e_psi2_offset_fvtxs1[counter] = epar2;
  e_psi2_magnitude_fvtxs1[counter] = (epar1/par1)*mag;
  psi2_difference_fvtxs1[counter] = (max-min)/(max);
  psi2_center_fvtxs1[counter] = cen;

  max = th1d_psi2_fvtxs2->GetMaximum();
  cen = th1d_psi2_fvtxs2->GetBinCenter(th1d_psi2_fvtxs2->GetMaximumBin());
  min = th1d_psi2_fvtxs2->GetMinimum(max*0.15); // minimum bin with at least 1 entry
  fun->SetParameter(0,max);
  fun->SetParameter(1,max);
  fun->SetParameter(2,0.0);
  th1d_psi2_fvtxs2->SetTitle("FVTXS layer 2");
  th1d_psi2_fvtxs2->Draw();
  th1d_psi2_fvtxs2->Fit(fun,"","",-1.55,1.55);
  if ( line_max ) delete line_max;
  line_max = new TLine(-4,max,4,max);
  line_max->SetLineStyle(2);
  line_max->SetLineWidth(2);
  line_max->Draw();
  if ( line_min ) delete line_min;
  line_min = new TLine(-4,min,4,min);
  line_min->SetLineStyle(2);
  line_min->SetLineWidth(2);
  line_min->Draw();
  if ( line_cen ) delete line_cen;
  line_cen = new TLine(cen,0,cen,max);
  line_cen->SetLineStyle(2);
  line_cen->SetLineWidth(2);
  line_cen->Draw();
  c1->Print(Form("Jamie/Figs/psi2fit_fvtxs2_run%d.png",run));
  c1->Print(Form("Jamie/Figs/psi2fit_fvtxs2_run%d.pdf",run));
  par0 = fun->GetParameter(0);
  par1 = fun->GetParameter(1);
  par2 = fun->GetParameter(2);
  epar0 = fun->GetParError(0);
  epar1 = fun->GetParError(1);
  epar2 = fun->GetParError(2);
  if ( par0 > 0 ) mag = par1/par0;
  off = par2;
  psi2_offset_fvtxs2[counter] = off;
  psi2_magnitude_fvtxs2[counter] = mag;
  e_psi2_offset_fvtxs2[counter] = epar2;
  e_psi2_magnitude_fvtxs2[counter] = (epar1/par1)*mag;
  psi2_difference_fvtxs2[counter] = (max-min)/(max);
  psi2_center_fvtxs2[counter] = cen;

  max = th1d_psi2_fvtxs3->GetMaximum();
  cen = th1d_psi2_fvtxs3->GetBinCenter(th1d_psi2_fvtxs3->GetMaximumBin());
  min = th1d_psi2_fvtxs3->GetMinimum(max*0.15); // minimum bin with at least 1 entry
  fun->SetParameter(0,max);
  fun->SetParameter(1,max);
  fun->SetParameter(2,0.0);
  th1d_psi2_fvtxs3->SetTitle("FVTXS layer 3");
  th1d_psi2_fvtxs3->Draw();
  th1d_psi2_fvtxs3->Fit(fun,"","",-1.55,1.55);
  if ( line_max ) delete line_max;
  line_max = new TLine(-4,max,4,max);
  line_max->SetLineStyle(2);
  line_max->SetLineWidth(2);
  line_max->Draw();
  if ( line_min ) delete line_min;
  line_min = new TLine(-4,min,4,min);
  line_min->SetLineStyle(2);
  line_min->SetLineWidth(2);
  line_min->Draw();
  if ( line_cen ) delete line_cen;
  line_cen = new TLine(cen,0,cen,max);
  line_cen->SetLineStyle(2);
  line_cen->SetLineWidth(2);
  line_cen->Draw();
  c1->Print(Form("Jamie/Figs/psi2fit_fvtxs3_run%d.png",run));
  c1->Print(Form("Jamie/Figs/psi2fit_fvtxs3_run%d.pdf",run));
  par0 = fun->GetParameter(0);
  par1 = fun->GetParameter(1);
  par2 = fun->GetParameter(2);
  epar0 = fun->GetParError(0);
  epar1 = fun->GetParError(1);
  epar2 = fun->GetParError(2);
  if ( par0 > 0 ) mag = par1/par0;
  off = par2;
  psi2_offset_fvtxs3[counter] = off;
  psi2_magnitude_fvtxs3[counter] = mag;
  e_psi2_offset_fvtxs3[counter] = epar2;
  e_psi2_magnitude_fvtxs3[counter] = (epar1/par1)*mag;
  psi2_difference_fvtxs3[counter] = (max-min)/(max);
  psi2_center_fvtxs3[counter] = cen;



  delete fun;
  delete c1;

}
