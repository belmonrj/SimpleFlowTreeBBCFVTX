void run_flattening(int run = 454936)
{

  gROOT->ProcessLine(".L post_ana_bbc_fvtx.cpp++");

  post_ana_bbc_fvtx(run,1);
  gROOT->Reset();
  post_ana_bbc_fvtx(run,2);
  gROOT->Reset();
  post_ana_bbc_fvtx(run,3);

}
