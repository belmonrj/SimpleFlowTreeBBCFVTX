void doenergy(int, int);

void calculate_epreso()
{

  doenergy(200,2);
  doenergy(62,2);
  doenergy(39,2);
  doenergy(20,2);

  // doenergy(200,3);
  // doenergy(62,3);
  // doenergy(39,3);
  // doenergy(20,3);

}

void doenergy(int energy, int harmonic)
{

  gStyle->SetOptTitle(1);

  TFile* file = TFile::Open(Form("input/combined_%d.root",energy));

  // ---
  // ---
  // ---

  TH1D* th1d_FVTXN_nclus = (TH1D*)file->Get("th1d_FVTXN_nclus");
  TH1D* th1d_FVTXS_nclus = (TH1D*)file->Get("th1d_FVTXS_nclus");
  TH1D* th1d_BBC_charge = (TH1D*)file->Get("th1d_BBC_charge");

  float mean_fvtxn = th1d_FVTXN_nclus->GetMean();
  float mean_fvtxs = th1d_FVTXS_nclus->GetMean();
  float mean_bbcs = th1d_BBC_charge->GetMean();

  // ---
  // ---
  // ---


  TProfile* tp1f_BBCS_FVTXN = (TProfile*)file->Get(Form("tp1f_reso%d_BBC_FVTXN",harmonic));
  TProfile* tp1f_BBCS_FVTXS = (TProfile*)file->Get(Form("tp1f_reso%d_BBC_FVTX",harmonic));
  TProfile* tp1f_FVTXN_FVTXS = (TProfile*)file->Get(Form("tp1f_reso%d_FVTXS_FVTXN",harmonic));
  // ---
  float float_BBCS_FVTXN = tp1f_BBCS_FVTXN->GetBinContent(1);
  float float_BBCS_FVTXS = tp1f_BBCS_FVTXS->GetBinContent(1);
  float float_FVTXN_FVTXS = tp1f_FVTXN_FVTXS->GetBinContent(1);
  float efloat_BBCS_FVTXN = tp1f_BBCS_FVTXN->GetBinError(1);
  float efloat_BBCS_FVTXS = tp1f_BBCS_FVTXS->GetBinError(1);
  float efloat_FVTXN_FVTXS = tp1f_FVTXN_FVTXS->GetBinError(1);
  // ---
  float reso_BBCS_fn = sqrt((float_BBCS_FVTXN*float_BBCS_FVTXS)/float_FVTXN_FVTXS); // BNBS/NS
  float reso_FVTXS_fn = sqrt((float_FVTXN_FVTXS*float_BBCS_FVTXS)/float_BBCS_FVTXN); // NSBS/BN
  float reso_FVTXN_fn = sqrt((float_FVTXN_FVTXS*float_BBCS_FVTXN)/float_BBCS_FVTXS); // NSBN/BS

  TProfile* tp1f_BBCS_CNT = (TProfile*)file->Get(Form("tp1f_reso%d_BBC_CNT",harmonic));
  TProfile* tp1f_CNT_FVTXS = (TProfile*)file->Get(Form("tp1f_reso%d_CNT_FVTX",harmonic));
  TProfile* tp1f_CNT_FVTXN = (TProfile*)file->Get(Form("tp1f_reso%d_CNT_FVTXN",harmonic));
  // ---
  float float_BBCS_CNT = tp1f_BBCS_CNT->GetBinContent(1);
  float float_CNT_FVTXS = tp1f_CNT_FVTXS->GetBinContent(1);
  float float_CNT_FVTXN = tp1f_CNT_FVTXN->GetBinContent(1);
  float efloat_BBCS_CNT = tp1f_BBCS_CNT->GetBinError(1);
  float efloat_CNT_FVTXS = tp1f_CNT_FVTXS->GetBinError(1);
  float efloat_CNT_FVTXN = tp1f_CNT_FVTXN->GetBinError(1);
  // ---
  float reso_BBCS = sqrt((float_BBCS_CNT*float_BBCS_FVTXS)/float_CNT_FVTXS); // BCBS/CS
  float reso_FVTXS = sqrt((float_CNT_FVTXS*float_BBCS_FVTXS)/float_BBCS_CNT); // CSBS/BC
  float reso_FVTXN = sqrt((float_CNT_FVTXN*float_BBCS_FVTXN)/float_BBCS_CNT); // CNBN/BC
  float reso_FVTXN_xb = sqrt ( ( float_FVTXN_FVTXS * float_CNT_FVTXN ) / float_CNT_FVTXS ) ; // NSNC/SC

  // ---
  float ereso_BBCS = sqrt( ( efloat_BBCS_CNT*efloat_BBCS_CNT/4*float_BBCS_CNT )
                           + ( efloat_BBCS_FVTXS*efloat_BBCS_FVTXS/4*float_BBCS_FVTXS )
                           + ( efloat_CNT_FVTXS*efloat_CNT_FVTXS/4*pow(float_CNT_FVTXS,3) ) );

  // --- event planes and correlations using CNT
  TString data1 = Form("%d GeV & %.2e & %.2e & %.2e \$\\pm\$ %.2e & %.2e \$\\pm\$ %.2e & %.2e \$\\pm\$ %.2e \\\\",
                       energy,reso_BBCS,reso_FVTXS,
                       float_BBCS_FVTXS,efloat_BBCS_FVTXS,
                       float_BBCS_CNT,efloat_BBCS_CNT,
                       float_CNT_FVTXS,efloat_CNT_FVTXS);

  // --- event planes and correlations using FVTXN
  TString data2 = Form("%d GeV & %.2e & %.2e & %.2e & %.2e & %.2e & %.2e \\\\",
                       energy, reso_BBCS_fn, reso_FVTXS_fn, reso_FVTXN_fn,
                       float_BBCS_FVTXS, float_BBCS_FVTXN, float_FVTXN_FVTXS);

  TString data3 = Form("%d GeV & %.2e & %.2e & %.2e & %.2e \\\\",
                       energy, mean_bbcs, mean_fvtxs, reso_BBCS, reso_FVTXS);

  //cout << data1.Data() << endl;
  //cout << data2.Data() << endl;
  //cout << data3.Data() << endl;

  float nreso_CNT = sqrt((float_CNT_FVTXN*float_CNT_FVTXS)/float_CNT_FVTXS); // BCBS/CS
  float nreso_FVTXS = sqrt((float_CNT_FVTXS*float_FVTXN_FVTXS)/float_CNT_FVTXN); // CSBS/BC
  float nreso_FVTXN = sqrt((float_CNT_FVTXN*float_FVTXN_FVTXS)/float_CNT_FVTXS); // CNBN/BC

  // --- event planes and correlations using CNT
  TString data4 = Form("%d GeV & %.2e & %.2e & %.2e \$\\pm\$ %.2e & %.2e \$\\pm\$ %.2e & %.2e \$\\pm\$ %.2e \\\\",
                       energy,nreso_FVTXN,nreso_FVTXS,
                       float_FVTXN_FVTXS,efloat_FVTXN_FVTXS,
                       float_CNT_FVTXN,efloat_CNT_FVTXN,
                       float_CNT_FVTXS,efloat_CNT_FVTXS);

  //  cout << data4.Data() << endl;
  // --- very roughly 85%
  float ratio = nreso_FVTXS/reso_FVTXS;
  cout << nreso_FVTXS << " " << reso_FVTXS << " " << ratio << endl;

  float ave = ( nreso_FVTXS + reso_FVTXS_fn ) / 2.0;
  ratio = ave/reso_FVTXS;
  cout << nreso_FVTXS << " " << reso_FVTXS << " " << reso_FVTXS_fn << " " << ave << " " << ratio << endl;



}

