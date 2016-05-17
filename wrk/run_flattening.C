void run_flattening()
{

  gROOT->ProcessLine(".L post_ana_bbc_fvtx.cpp++");

  post_ana_bbc_fvtx(454811,1);
  gROOT->Reset();
  post_ana_bbc_fvtx(454811,2);
  gROOT->Reset();
  post_ana_bbc_fvtx(454811,3);

}
