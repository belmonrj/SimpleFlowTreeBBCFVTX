void plotv3ew(){
  TFile *f = new TFile("combined_200.root");
  
  TH1F *hcosw=bbcs_v3_west_cosphi->Clone();
  TH1F *hcose=bbcs_v3_east_cosphi->Clone();

  TH1F *v3allbbc=bbcs_v3_both_docalib->Clone();
  TH1F *v3westbbc=bbcs_v3_west_docalib->Clone();
  TH1F *v3eastbbc=bbcs_v3_east_docalib->Clone();

  TH1F *v3allfvtx=fvtxs_v3_both_docalib->Clone();
  TH1F *v3westfvtx=fvtxs_v3_west_docalib->Clone();
  TH1F *v3eastfvtx=fvtxs_v3_east_docalib->Clone();
  /*
  uhadronvn_2_4->ProfileX("v3allfvtx");
  uhadeastvn_2_4->ProfileX("v3eastfvtx");
  uhadwestvn_2_4->ProfileX("v3westfvtx");

  uhadronvn_2_5->ProfileX("v3allbbc");
  uhadeastvn_2_5->ProfileX("v3eastbbc");
  uhadwestvn_2_5->ProfileX("v3westbbc");
  */
  TH1F *v3eastfvtxnew = new TH1F("v3eastfvtxnew","v3eastfvtxnew",15,0,3.0);
  TH1F *v3westfvtxnew = new TH1F("v3westfvtxnew","v3westfvtxnew",15,0,3.0);

  TH1F *v3eastbbcnew = new TH1F("v3eastbbcnew","v3eastbbcnew",15,0,3.0);
  TH1F *v3westbbcnew = new TH1F("v3westbbcnew","v3westbbcnew",15,0,3.0);

  /*
  TH1F *v3eastfvtxnew=static_cast<TH1F *>v3eastfvtx->Clone();
  TH1F *v3westfvtxnew=static_cast<TH1F *>v3westfvtx->Clone();

  TH1F *v3eastbbcnew=static_cast<TH1F *>v3eastbbc->Clone();
  TH1F *v3westbbcnew=static_cast<TH1F *>v3westbbc->Clone();
  
  TH1F *hcoswfvtx=static_cast<TH1F *>hcosw->Clone();
  TH1F *hcosefvtx=static_cast<TH1F *>hcose->Clone();
  
  TH1F *hcoswbbc=static_cast<TH1F *>hcosw->Clone();
  TH1F *hcosebbc=static_cast<TH1F *>hcose->Clone();
  */

  float qxfvtx=0;
  //hcosefvtx->Scale(qxfvtx);
  //hcoswfvtx->Scale(qxfvtx);

  float qxbbc=0.002;
  //hcosebbc->Scale(qxbbc);
  //hcoswbbc->Scale(qxbbc);
  //FVTX
  for(int i=0; i<15; i++){

    float cone = v3eastfvtx->GetBinContent(i+1)+qxfvtx*hcose->GetBinContent(i+1);
    float econe = sqrt(v3eastfvtx->GetBinError(i+1)**2+(qxfvtx*hcose->GetBinError(i+1))**2);
    v3eastfvtxnew->SetBinContent(i+1,cone);
    v3eastfvtxnew->SetBinError(i+1,econe);

    float conw = v3westfvtx->GetBinContent(i+1)+qxfvtx*hcosw->GetBinContent(i+1);
    float econw = sqrt(v3westfvtx->GetBinError(i+1)**2+(qxfvtx*hcosw->GetBinError(i+1))**2);
    v3westfvtxnew->SetBinContent(i+1,conw);
    v3westfvtxnew->SetBinError(i+1,econw);

    //cout<<cone<<" "<<conw<<endl;
  }

  //BBC
  for(int i=0; i<15; i++){

    float cone = v3eastbbc->GetBinContent(i+1)+qxbbc*hcose->GetBinContent(i+1);
    float econe = sqrt(v3eastbbc->GetBinError(i+1)**2+(qxbbc*hcose->GetBinError(i+1))**2);
    v3eastbbcnew->SetBinContent(i+1,cone);
    v3eastbbcnew->SetBinError(i+1,econe);

    float conw = v3westbbc->GetBinContent(i+1)+qxbbc*hcosw->GetBinContent(i+1);
    float econw = sqrt(v3westbbc->GetBinError(i+1)**2+(qxbbc*hcosw->GetBinError(i+1))**2);
    v3westbbcnew->SetBinContent(i+1,conw);
    v3westbbcnew->SetBinError(i+1,econw);

    //cout<<cone<<" "<<conw<<endl;
  }

  //for(int i=0; i<15; i++){
  //  
  //}
  //v3eastfvtxnew->Add(hcosefvtx,1.0);
  //v3westfvtxnew->Add(hcoswfvtx,1.0);

  

 

  //v3eastbbcnew->Add(hcosebbc,1.0);
  //v3westbbcnew->Add(hcoswbbc,1.0);
  
  float resofvtx=0.050;//0.0388085;
  float resobbc=0.040;//0.0329477;

  v3allfvtx->Scale(1.0/resofvtx);
  v3eastfvtxnew->Scale(1.0/resofvtx);
  v3westfvtxnew->Scale(1.0/resofvtx);

  
  v3allbbc->Scale(1.0/resobbc);
  v3eastbbcnew->Scale(1.0/resobbc);
  v3westbbcnew->Scale(1.0/resobbc);

  cc = new TCanvas("cc","cc",800,400);
  cc->SetFillColor(10);
  cc->Divide(2,1);
  cc->Draw();

  cc->cd(1);
  TH1F *h0=new TH1F("h0","FVTX", 30, 0, 3.0);
  h0->SetMaximum(0.1);
  h0->SetMinimum(0.0);
  h0->GetXaxis()->SetTitle("pt");
  h0->GetYaxis()->SetTitle("v3");
  h0->GetYaxis()->CenterTitle(kTRUE);
  h0->Draw();

  v3allfvtx->SetMarkerStyle(22);
  v3allfvtx->SetMarkerSize(1.2);
  v3allfvtx->SetMarkerColor(6);
  v3allfvtx->Draw("Psame");

  v3eastfvtxnew->SetMarkerStyle(20);
  v3eastfvtxnew->SetMarkerSize(1.2);
  v3eastfvtxnew->SetMarkerColor(2);
  v3eastfvtxnew->Draw("Psame");

  v3westfvtxnew->SetMarkerStyle(24);
  v3westfvtxnew->SetMarkerSize(1.2);
  v3westfvtxnew->SetMarkerColor(4);
  v3westfvtxnew->Draw("Psame");
  
  cc->cd(2);
  TH1F *h1=new TH1F("h1","BBC", 30, 0, 3.0);
  h1->SetMaximum(0.1);
  h1->SetMinimum(0.0);
  h1->GetXaxis()->SetTitle("pt");
  h1->GetYaxis()->SetTitle("v3");
  h1->GetYaxis()->CenterTitle(kTRUE);
  h1->Draw();

  v3allbbc->SetMarkerStyle(22);
  v3allbbc->SetMarkerSize(1.2);
  v3allbbc->SetMarkerColor(6);
  v3allbbc->Draw("Psame");

  v3eastbbcnew->SetMarkerStyle(20);
  v3eastbbcnew->SetMarkerSize(1.2);
  v3eastbbcnew->SetMarkerColor(2);
  v3eastbbcnew->Draw("Psame");

  v3westbbcnew->SetMarkerStyle(24);
  v3westbbcnew->SetMarkerSize(1.2);
  v3westbbcnew->SetMarkerColor(4);
  v3westbbcnew->Draw("Psame");

  TLegend *leg1 = new TLegend(0.65,0.68,0.90,0.88);
  leg1->SetFillColor(10);
  leg1->SetLineStyle(4000);
  leg1->SetLineColor(10);
  leg1->SetLineWidth(0.);
  leg1->SetTextSize(0.04);
  leg1->SetBorderSize(0);

  leg1->AddEntry(v3allbbc,"ALL","P");
  leg1->AddEntry(v3eastbbcnew,"East","P");
  leg1->AddEntry(v3westbbcnew,"West","P");
  leg1->Draw();
}
