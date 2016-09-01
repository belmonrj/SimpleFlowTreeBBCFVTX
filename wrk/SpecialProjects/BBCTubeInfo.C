float d_pmt_x[64];
float d_pmt_y[64];
float d_pmt_z = -1443.5; // same for all tubes

float radius_layer_inner = 85.2-10.0;
float radius_layer_innermiddle = 85.2;
float radius_layer_middle = 113.6-13.0;
float radius_layer_outermiddle = 113.6;
float radius_layer_outer = 142.0-16.0;

void initialize_pmt_position();

int get_pmt_layer(int);

void BBCTubeInfo()
{

  TCanvas* c1 = new TCanvas("c1","",800,800);

  initialize_pmt_position();

  float d_pmt_r[64];
  float d_pmt_phi[64];

  float tube_number[64];

  for ( int i = 0; i < 64; ++i )
    {
      d_pmt_r[i] = sqrt ( ( d_pmt_x[i]*d_pmt_x[i] ) + ( d_pmt_y[i]*d_pmt_y[i] ) );
      d_pmt_phi[i] = atan2 ( d_pmt_y[i], d_pmt_x[i] );
      tube_number[i] = i;
    }

  TGraph* tg_r = new TGraph(64,tube_number,d_pmt_r);
  TGraph* tg_phi = new TGraph(64,tube_number,d_pmt_phi);

  tg_r->SetMarkerStyle(kFullCircle);
  tg_r->SetMinimum(0.0);
  tg_r->SetMaximum(140.0);
  tg_r->GetXaxis()->SetLimits(-1,64);
  tg_r->Draw("ap");
  c1->Print("bbc_tube_radius.png");
  c1->Print("bbc_tube_radius.pdf");

  tg_phi->SetMarkerStyle(kFullCircle);
  tg_phi->SetMinimum(-3.2);
  tg_phi->SetMaximum(3.2);
  tg_phi->GetXaxis()->SetLimits(-1,64);
  tg_phi->Draw("ap");
  c1->Print("bbc_tube_phi.png");
  c1->Print("bbc_tube_phi.pdf");

  TGraph* tg_xy = new TGraph(64,d_pmt_x,d_pmt_y);
  tg_xy->SetMarkerStyle(kFullCircle);
  tg_xy->SetMinimum(-140);
  tg_xy->SetMaximum(140);
  tg_xy->GetXaxis()->SetLimits(-140,140);
  tg_xy->Draw("ap");
  TEllipse* el0 = new TEllipse(0,0,85.2-10,85.2-10);
  el0->SetFillStyle(0);
  el0->SetLineColor(kRed);
  el0->Draw();
  TEllipse* el1 = new TEllipse(0,0,113.6-11,113.6-13);
  el1->SetFillStyle(0);
  el1->SetLineColor(kGreen+2);
  el1->Draw();
  TEllipse* el2 = new TEllipse(0,0,142.0-12,142.0-16);
  el2->SetFillStyle(0);
  el2->SetLineColor(kBlue);
  el2->Draw();
  c1->Print("bbc_tube_xy.png");
  c1->Print("bbc_tube_xy.pdf");

  TMarker* tmtube[64];
  TLatex* tex = new TLatex();
  for ( int i = 0; i < 64; ++i )
    {
      float x = d_pmt_x[i];
      float y = d_pmt_y[i];
      tex->DrawLatex(x,y,Form("%d",i));
      int layer = get_pmt_layer(i);
      if ( layer == -1 ) cout << "no layer info for tube number " << i << endl;
      tmtube[i] = new TMarker(d_pmt_x[i],d_pmt_y[i],1);
      tmtube[i]->SetMarkerStyle(kFullCircle);
      if ( layer == 0 ) tmtube[i]->SetMarkerColor(kRed);
      if ( layer == 2 ) tmtube[i]->SetMarkerColor(kGreen+2);
      if ( layer == 4 ) tmtube[i]->SetMarkerColor(kBlue);
      if ( layer == 1 ) tmtube[i]->SetMarkerColor(kOrange+2);
      if ( layer == 3 ) tmtube[i]->SetMarkerColor(kMagenta+2);
      tmtube[i]->Draw();
    }
  c1->Print("bbc_tubenumber_xy.png");
  c1->Print("bbc_tubenumber_xy.pdf");


}


int get_pmt_layer(int i)
{

  if ( i == 8 ||
       i == 11||
       i == 14||
       i == 19||
       i == 22||
       i == 26||
       i == 40||
       i == 43||
       i == 46||
       i == 51||
       i == 54||
       i == 58 ) return 0; // inner layer

  if ( i == 7 ||
       i == 16||
       i == 25||
       i == 39||
       i == 48||
       i == 57 ) return 1; // inner middle layer

  if ( i == 4 ||
       i == 6 ||
       i == 10||
       i == 13||
       i == 18||
       i == 21||
       i == 24||
       i == 29||
       i == 36||
       i == 38||
       i == 42||
       i == 45||
       i == 45||
       i == 50||
       i == 53||
       i == 56||
       i == 61 ) return 2; // middle layer

  if ( i == 3 ||
       i == 15||
       i == 28||
       i == 35||
       i == 47||
       i == 60 ) return 3; //outer middle layer

  if ( i == 0 ||
       i == 1 ||
       i == 2 ||
       i == 5 ||
       i == 9 ||
       i == 12||
       i == 17||
       i == 20||
       i == 23||
       i == 27||
       i == 30||
       i == 31||
       i == 32||
       i == 33||
       i == 34||
       i == 37||
       i == 41||
       i == 44||
       i == 49||
       i == 52||
       i == 55||
       i == 59||
       i == 62||
       i == 63 ) return 4; // outer layer

  return -1;

}


void initialize_pmt_position()
{

  d_pmt_x[0] = -123;
  d_pmt_y[0] = 42.6;
  d_pmt_x[1] = -123;
  d_pmt_y[1] = 14.2;
  d_pmt_x[2] = -98.4;
  d_pmt_y[2] = 85.2;
  d_pmt_x[3] = -98.4;
  d_pmt_y[3] = 56.8;
  d_pmt_x[4] = -98.4;
  d_pmt_y[4] = 28.4;
  d_pmt_x[5] = -73.8;
  d_pmt_y[5] = 99.4;
  d_pmt_x[6] = -73.8;
  d_pmt_y[6] = 71;
  d_pmt_x[7] = -73.8;
  d_pmt_y[7] = 42.6;
  d_pmt_x[8] = -73.8;
  d_pmt_y[8] = 14.2;
  d_pmt_x[9] = -49.2;
  d_pmt_y[9] = 113.6;
  d_pmt_x[10] = -49.2;
  d_pmt_y[10] = 85.2;
  d_pmt_x[11] = -49.2;
  d_pmt_y[11] = 56.8;
  d_pmt_x[12] = -24.6;
  d_pmt_y[12] = 127.8;
  d_pmt_x[13] = -24.6;
  d_pmt_y[13] = 99.4;
  d_pmt_x[14] = -24.6;
  d_pmt_y[14] = 71;
  d_pmt_x[15] = 0;
  d_pmt_y[15] = 113.6;
  d_pmt_x[16] = 0;
  d_pmt_y[16] = 85.2;
  d_pmt_x[17] = 24.6;
  d_pmt_y[17] = 127.8;
  d_pmt_x[18] = 24.6;
  d_pmt_y[18] = 99.4;
  d_pmt_x[19] = 24.6;
  d_pmt_y[19] = 71;
  d_pmt_x[20] = 49.2;
  d_pmt_y[20] = 113.6;
  d_pmt_x[21] = 49.2;
  d_pmt_y[21] = 85.2;
  d_pmt_x[22] = 49.2;
  d_pmt_y[22] = 56.8;
  d_pmt_x[23] = 73.8;
  d_pmt_y[23] = 99.4;
  d_pmt_x[24] = 73.8;
  d_pmt_y[24] = 71;
  d_pmt_x[25] = 73.8;
  d_pmt_y[25] = 42.6;
  d_pmt_x[26] = 73.8;
  d_pmt_y[26] = 14.2;
  d_pmt_x[27] = 98.4;
  d_pmt_y[27] = 85.2;
  d_pmt_x[28] = 98.4;
  d_pmt_y[28] = 56.8;
  d_pmt_x[29] = 98.4;
  d_pmt_y[29] = 28.4;
  d_pmt_x[30] = 123;
  d_pmt_y[30] = 42.6;
  d_pmt_x[31] = 123;
  d_pmt_y[31] = 14.2;
  d_pmt_x[32] = 123;
  d_pmt_y[32] = -42.6;
  d_pmt_x[33] = 123;
  d_pmt_y[33] = -14.2;
  d_pmt_x[34] = 98.4;
  d_pmt_y[34] = -85.2;
  d_pmt_x[35] = 98.4;
  d_pmt_y[35] = -56.8;
  d_pmt_x[36] = 98.4;
  d_pmt_y[36] = -28.4;
  d_pmt_x[37] = 73.8;
  d_pmt_y[37] = -99.4;
  d_pmt_x[38] = 73.8;
  d_pmt_y[38] = -71;
  d_pmt_x[39] = 73.8;
  d_pmt_y[39] = -42.6;
  d_pmt_x[40] = 73.8;
  d_pmt_y[40] = -14.2;
  d_pmt_x[41] = 49.2;
  d_pmt_y[41] = -113.6;
  d_pmt_x[42] = 49.2;
  d_pmt_y[42] = -85.2;
  d_pmt_x[43] = 49.2;
  d_pmt_y[43] = -56.8;
  d_pmt_x[44] = 24.6;
  d_pmt_y[44] = -127.8;
  d_pmt_x[45] = 24.6;
  d_pmt_y[45] = -99.4;
  d_pmt_x[46] = 24.6;
  d_pmt_y[46] = -71;
  d_pmt_x[47] = -0;
  d_pmt_y[47] = -113.6;
  d_pmt_x[48] = -0;
  d_pmt_y[48] = -85.2;
  d_pmt_x[49] = -24.6;
  d_pmt_y[49] = -127.8;
  d_pmt_x[50] = -24.6;
  d_pmt_y[50] = -99.4;
  d_pmt_x[51] = -24.6;
  d_pmt_y[51] = -71;
  d_pmt_x[52] = -49.2;
  d_pmt_y[52] = -113.6;
  d_pmt_x[53] = -49.2;
  d_pmt_y[53] = -85.2;
  d_pmt_x[54] = -49.2;
  d_pmt_y[54] = -56.8;
  d_pmt_x[55] = -73.8;
  d_pmt_y[55] = -99.4;
  d_pmt_x[56] = -73.8;
  d_pmt_y[56] = -71;
  d_pmt_x[57] = -73.8;
  d_pmt_y[57] = -42.6;
  d_pmt_x[58] = -73.8;
  d_pmt_y[58] = -14.2;
  d_pmt_x[59] = -98.4;
  d_pmt_y[59] = -85.2;
  d_pmt_x[60] = -98.4;
  d_pmt_y[60] = -56.8;
  d_pmt_x[61] = -98.4;
  d_pmt_y[61] = -28.4;
  d_pmt_x[62] = -123;
  d_pmt_y[62] = -42.6;
  d_pmt_x[63] = -123;
  d_pmt_y[63] = -14.2;

}



