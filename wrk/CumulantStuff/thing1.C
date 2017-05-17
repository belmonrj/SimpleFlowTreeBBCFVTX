
const int number_of_tests = 6; // update as needed
const int xhash = 0;
const int xdesc = 1;
TString hash_description[2][number_of_tests];

//void dothething(int, int, int, int, int);

void thing1()
{

  hash_description[0][0] = "aea49a4";  hash_description[1][0] = "Baseline"; // earlier baseline, different from others, same as ce076aa
  hash_description[0][1] = "b917033";  hash_description[1][1] = "FtD  2dR"; // floats changed to doubles and uses 2d recentering
  hash_description[0][2] = "75765b7";  hash_description[1][2] = "FtD only"; // floats changed to doubles but uses 1d recentering
  hash_description[0][3] = "9fefb05";  hash_description[1][3] = "2dR only"; // uses 2d recentering but floats left as floats
  hash_description[0][4] = "d52d914";  hash_description[1][4] = "Baseline"; // most recent baseline, agrees with modified ones
  hash_description[0][5] = "ce076aa";  hash_description[1][5] = "Baseline"; // earlier baseline, different from others, same as aea49a4

  gROOT->ProcessLine(".L dothething_components.C");
  dothething_components(200,1,0,4,5);
  dothething_components(200,4,3,2,1);

  gROOT->ProcessLine(".L dothething_cumulants.C");
  dothething_cumulants(200,1,0,4,5);
  dothething_cumulants(200,4,3,2,1);

}

