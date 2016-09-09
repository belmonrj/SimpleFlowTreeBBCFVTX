void doenergy(int, int);

void calculate_epreso()
{

  doenergy(200,2);
  doenergy(62,2);
  doenergy(39,2);
  doenergy(20,2);

  doenergy(200,3);
  doenergy(62,3);
  doenergy(39,3);
  doenergy(20,3);

}

void doenergy(int energy, int harmonic)
{

  gStyle->SetOptTitle(1);

  TFile* file = TFile::Open(Form("input/combined_%d.root",energy));

  // ---

  TProfile* tp1f_bbc_fvtxN = (TProfile*)file->Get(Form("tp1f_reso%d_BBC_FVTXN",harmonic));
  TProfile* tp1f_bbc_fvtx = (TProfile*)file->Get(Form("tp1f_reso%d_BBC_FVTX",harmonic));
  TProfile* tp1f_fvtxN_fvtx = (TProfile*)file->Get(Form("tp1f_reso%d_FVTXS_FVTXN",harmonic));

  float float_bbc_fvtxN = tp1f_bbc_fvtxN->GetBinContent(1);
  float float_bbc_fvtx = tp1f_bbc_fvtx->GetBinContent(1);
  float float_fvtxN_fvtx = tp1f_fvtxN_fvtx->GetBinContent(1);

  // cout << "bbc-fvtxn " << float_bbc_fvtxN << endl;
  // cout << "bbc-fvtxs " << float_bbc_fvtx << endl;
  // cout << "fvtxn-fvtxs " << float_fvtxN_fvtx << endl;

  float reso_bbc_fn = sqrt((float_bbc_fvtxN*float_bbc_fvtx)/float_fvtxN_fvtx); // BNBS/NS
  float reso_fvtx_fn = sqrt((float_fvtxN_fvtx*float_bbc_fvtx)/float_bbc_fvtxN); // NSBS/BN
  float reso_fvtxN_fn = sqrt((float_fvtxN_fvtx*float_bbc_fvtxN)/float_bbc_fvtx); // NSBN/BS

  // cout << "bbc resolution is " << reso_bbc_fn << endl;
  // cout << "fvtx resolution is " << reso_fvtx_fn << endl;
  // cout << "fvtxN resolution is " << reso_fvtxN_fn << endl;

  // ---

  TProfile* tp1f_bbc_cnt = (TProfile*)file->Get(Form("tp1f_reso%d_BBC_CNT",harmonic));
  TProfile* tp1f_cnt_fvtx = (TProfile*)file->Get(Form("tp1f_reso%d_CNT_FVTX",harmonic));
  TProfile* tp1f_cnt_fvtxN = (TProfile*)file->Get(Form("tp1f_reso%d_CNT_FVTXN",harmonic));

  float float_bbc_cnt = tp1f_bbc_cnt->GetBinContent(1);
  float float_cnt_fvtx = tp1f_cnt_fvtx->GetBinContent(1);
  float float_cnt_fvtxN = tp1f_cnt_fvtxN->GetBinContent(1);

  // cout << "bbc-cnt " << float_bbc_cnt << endl;
  // cout << "cnt-fvtxs " << float_cnt_fvtx << endl;
  // cout << "cbt-fvtxn " << float_cnt_fvtxN << endl;

  if ( float_bbc_cnt < 0 ) { /*float_bbc_cnt *= -1;*/ cout << "YOU'RE GONNA DIE (bbcs-cnt)" << endl; }
  if ( float_bbc_fvtx < 0 ) { /*float_bbc_fvtx *= -1;*/ cout << "YOU'RE GONNA DIE (cnt-fvtxs)" << endl; }
  if ( float_cnt_fvtx < 0 ) { /*float_cnt_fvtx *= -1;*/ cout << "YOU'RE GONNA DIE (cnt-fvtxs)" << endl; }
  if ( float_cnt_fvtxN < 0 ) { /*float_cnt_fvtxN *= -1;*/ cout << "YOU'RE GONNA DIE (cnt-fvtxn)" << endl; }

  float reso_bbc = sqrt((float_bbc_cnt*float_bbc_fvtx)/float_cnt_fvtx); // BCBS/CS
  float reso_fvtx = sqrt((float_cnt_fvtx*float_bbc_fvtx)/float_bbc_cnt); // CSBS/BC
  float reso_fvtxN = sqrt((float_cnt_fvtxN*float_bbc_fvtxN)/float_bbc_cnt); // CNBN/BC

  //  cout << "bbc resolution is " << reso_bbc << endl;
  //  cout << "fvtx resolution is " << reso_fvtx << endl;
  //  cout << "fvtxN resolution is " << reso_fvtxN << endl;

  // ---

  float reso_fvtxN_xb = sqrt ( ( float_fvtxN_fvtx * float_cnt_fvtxN ) / float_cnt_fvtx ) ; // NSNC/SC

  // cout << "fvtxN resolution is " << reso_fvtxN_xb << endl;

  cout << energy << " GeV & " << reso_bbc << " & " << reso_fvtx << " & " << float_bbc_fvtx << " & " << float_bbc_cnt << " & " << float_cnt_fvtx << " \\\\ " << endl;
  cout << energy << " GeV & " << reso_bbc_fn << " & " << reso_fvtx_fn << " & " << reso_fvtxN_fn << " & " << float_bbc_fvtx << " & " << float_bbc_fvtxN << " & " << float_fvtxN_fvtx << " \\\\ " << endl;

}

