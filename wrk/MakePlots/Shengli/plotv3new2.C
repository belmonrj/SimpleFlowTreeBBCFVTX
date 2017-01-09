void plotv3new2(){
  TFile *f = new TFile("combined_200.root");
  fvtxs_v3_both_docalib->Draw();
  bbcs_v3_both_docalib->Draw();

  float fvtx4_cnt=0; float con_fvtx4=0;//for FVTX resolution
  float fvtx5_cnt=0; float con_fvtx5=0;//for BBC resolution
  for(int i=2; i<15; i++){//pt>0.4 and pt<3 GeV/c
    fvtx4_cnt+=fvtxs_v3_both_docalib->GetBinContent(i+1)*fvtxs_v3_both_docalib->GetBinEntries(i+1);
    con_fvtx4+=fvtxs_v3_both_docalib->GetBinEntries(i+1);

    fvtx5_cnt+=bbcs_v3_both_docalib->GetBinContent(i+1)*bbcs_v3_both_docalib->GetBinEntries(i+1);
    con_fvtx5+=bbcs_v3_both_docalib->GetBinEntries(i+1);
  }
  
  fvtx4_cnt/=con_fvtx4;
  fvtx5_cnt/=con_fvtx5;

  float bbc_fvtx =  tp1f_reso3_BBC_FVTX->GetBinContent(1);

  float reso4 = sqrt(bbc_fvtx*fvtx4_cnt/fvtx5_cnt);//fvtx resolution
  float reso5 = sqrt(bbc_fvtx/fvtx4_cnt*fvtx5_cnt);//bbc resolution

  cout<<reso4<<" "<<reso5<<endl;

  TProfile *v3fvtx4 = fvtxs_v3_both_docalib->Clone();//fvtx v3
  TProfile *v3fvtx5 = bbcs_v3_both_docalib->Clone();//bbc v3

  
  v3fvtx4->Scale(1.0/reso4);
  v3fvtx5->Scale(1.0/reso5);

  cc = new TCanvas("cc","cc");
  cc->SetFillColor(10);
  //cc->Divide(2,3);
  cc->Draw();


  TH1F *h0=new TH1F("h0","Layer0", 30, 0, 3.0);
  h0->SetMaximum(0.1);
  h0->SetMinimum(0.0);
  h0->GetXaxis()->SetTitle("pt");
  h0->GetYaxis()->SetTitle("v3");
  h0->GetYaxis()->CenterTitle(kTRUE);
  h0->Draw();


  v3fvtx4->SetMarkerStyle(20);
  v3fvtx4->SetMarkerSize(1.2);
  v3fvtx4->SetMarkerColor(4);
  v3fvtx4->Draw("Psame");

  v3fvtx5->SetMarkerStyle(24);
  v3fvtx5->SetMarkerSize(1.2);
  v3fvtx5->SetMarkerColor(6);
  v3fvtx5->Draw("Psame");

  TLegend *leg1 = new TLegend(0.65,0.68,0.90,0.88);
  leg1->SetFillColor(10);
  leg1->SetLineStyle(4000);
  leg1->SetLineColor(10);
  leg1->SetLineWidth(0.);
  leg1->SetTextSize(0.04);
  leg1->SetBorderSize(0);

  leg1->AddEntry(v3fvtx4,"FVTX","P");
  leg1->AddEntry(v3fvtx5,"BBC","P");

leg1->Draw();
}
