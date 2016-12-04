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

float runnumber[111]; // actually an int but convenient to have a float

void doit(int);

void layers_rbr()
{

  gStyle->SetOptFit(1);

  int run;
  ifstream fin;

  counter = 0;
  fin.open("list200.txt");
  while ( fin >> run )
    {
      doit(run);
      ++counter;
      runnumber[counter] = run;
      if ( counter > 2 ) break;
    }
  fin.close();




}


void doit(int run)
{

  TFile* file = TFile::Open(Form("RootFiles/hist_%d.root",run));
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

  TF1* fun = new TF1("fun","[0]+[1]*cos(3*(x-[2]))",-1.1,1.1);
  fun->SetParLimits(2,-1.1,1.1);
  float par0, par1, par2;
  float epar0, epar1, epar2;
  float max;
  float off;
  float mag;

  max = th1d_psi3_fvtxs0->GetMaximum();
  fun->SetParameter(0,max);
  fun->SetParameter(1,max);
  fun->SetParameter(2,0.0);
  th1d_psi3_fvtxs0->Draw();
  th1d_psi3_fvtxs0->Fit(fun,"","",-1.1,1.1);
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

  max = th1d_psi3_fvtxs1->GetMaximum();
  fun->SetParameter(0,max);
  fun->SetParameter(1,max);
  fun->SetParameter(2,0.0);
  th1d_psi3_fvtxs1->Draw();
  th1d_psi3_fvtxs1->Fit(fun,"","",-1.1,1.1);
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

  max = th1d_psi3_fvtxs2->GetMaximum();
  fun->SetParameter(0,max);
  fun->SetParameter(1,max);
  fun->SetParameter(2,0.0);
  th1d_psi3_fvtxs2->Draw();
  th1d_psi3_fvtxs2->Fit(fun,"","",-1.1,1.1);
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

  max = th1d_psi3_fvtxs3->GetMaximum();
  fun->SetParameter(0,max);
  fun->SetParameter(1,max);
  fun->SetParameter(2,0.0);
  th1d_psi3_fvtxs3->Draw();
  th1d_psi3_fvtxs3->Fit(fun,"","",-1.1,1.1);
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


  delete fun;
  delete c1;

}
