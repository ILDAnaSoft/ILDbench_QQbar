// cross section : https://ild.ngt.ndu.ac.jp/CDS/mc-dbd.log/generated/metainfo-id/

TH1F* getHist(float epol = 0., float ppol = 0.)
{
  float luminosity = 4000.*0.4; // fb-1
  float w_L = (1.-epol)*(1.+ppol)/4.;
  float w_R = (1.+epol)*(1.-ppol)/4.;
  double xsec_Pyyuyyc_L = 1.6444286E+02; 
  double xsec_Pyyuyyc_R = 6.4372089E+01; 
  double xsec_Pyycyyu_L = 1.6547670E+02;
  double xsec_Pyycyyu_R = 6.4060427E+01; 
  double xsec_Pyyuyyu_L = 1.6656864E+02;
  double xsec_Pyyuyyu_R = 6.4540275E+01; 
  double xsec_Pyycyyc_L = 1.6332002E+02; 
  double xsec_Pyycyyc_R = 6.0651686E+01; 

  TString dirPath = gSystem->GetFromPipe("echo $INITDIR") + "/run_l5/test/root/";

  TString fpath_Pyyuyyu_L = dirPath + "l5_500GeV.Pyyuyyu.eL.pR_sum.root";
  TString fpath_Pyyuyyu_R = dirPath + "l5_500GeV.Pyyuyyu.eR.pL_sum.root";
  TString fpath_Pyycyyc_L = dirPath + "l5_500GeV.Pyycyyc.eL.pR_sum.root";
  TString fpath_Pyycyyc_R = dirPath + "l5_500GeV.Pyycyyc.eR.pL_sum.root";
  TString fpath_Pyyuyyc_L = dirPath + "l5_500GeV.Pyyuyyc.eL.pR_sum.root";
  TString fpath_Pyyuyyc_R = dirPath + "l5_500GeV.Pyyuyyc.eR.pL_sum.root";
  TString fpath_Pyycyyu_L = dirPath + "l5_500GeV.Pyycyyu.eL.pR_sum.root";
  TString fpath_Pyycyyu_R = dirPath + "l5_500GeV.Pyycyyu.eR.pL_sum.root";

  TFile* fin_Pyyuyyu_L = new TFile(fpath_Pyyuyyu_L.Data());
  TFile* fin_Pyyuyyu_R = new TFile(fpath_Pyyuyyu_R.Data());
  TFile* fin_Pyycyyc_L = new TFile(fpath_Pyycyyc_L.Data());
  TFile* fin_Pyycyyc_R = new TFile(fpath_Pyycyyc_R.Data());
  TFile* fin_Pyyuyyc_L = new TFile(fpath_Pyyuyyc_L.Data());
  TFile* fin_Pyyuyyc_R = new TFile(fpath_Pyyuyyc_R.Data());
  TFile* fin_Pyycyyu_L = new TFile(fpath_Pyycyyu_L.Data());
  TFile* fin_Pyycyyu_R = new TFile(fpath_Pyycyyu_R.Data());

  // histogram operation seems to require non-pointer object.
  TH1F h_Pyyuyyu_L = (TH1F)(*static_cast<TH1F*>(fin_Pyyuyyu_L->Get("MCCosThetaT")));
  TH1F h_Pyyuyyu_R = (TH1F)(*static_cast<TH1F*>(fin_Pyyuyyu_R->Get("MCCosThetaT")));
  TH1F h_Pyycyyc_L = (TH1F)(*static_cast<TH1F*>(fin_Pyycyyc_L->Get("MCCosThetaT")));
  TH1F h_Pyycyyc_R = (TH1F)(*static_cast<TH1F*>(fin_Pyycyyc_R->Get("MCCosThetaT")));
  TH1F h_Pyyuyyc_L = (TH1F)(*static_cast<TH1F*>(fin_Pyyuyyc_L->Get("MCCosThetaT")));
  TH1F h_Pyyuyyc_R = (TH1F)(*static_cast<TH1F*>(fin_Pyyuyyc_R->Get("MCCosThetaT")));
  TH1F h_Pyycyyu_L = (TH1F)(*static_cast<TH1F*>(fin_Pyycyyu_L->Get("MCCosThetaT")));
  TH1F h_Pyycyyu_R = (TH1F)(*static_cast<TH1F*>(fin_Pyycyyu_R->Get("MCCosThetaT")));

  int n_Pyyuyyu_L = h_Pyyuyyu_L.Integral();
  int n_Pyyuyyu_R = h_Pyyuyyu_R.Integral();
  int n_Pyycyyc_L = h_Pyycyyc_L.Integral();
  int n_Pyycyyc_R = h_Pyycyyc_R.Integral();
  int n_Pyyuyyc_L = h_Pyyuyyc_L.Integral();
  int n_Pyyuyyc_R = h_Pyyuyyc_R.Integral();
  int n_Pyycyyu_L = h_Pyycyyu_L.Integral();
  int n_Pyycyyu_R = h_Pyycyyu_R.Integral();
  float k_Pyyuyyu_L = w_L*xsec_Pyyuyyu_L*luminosity/n_Pyyuyyu_L;
  float k_Pyyuyyu_R = w_R*xsec_Pyyuyyu_R*luminosity/n_Pyyuyyu_R;
  float k_Pyycyyc_L = w_L*xsec_Pyycyyc_L*luminosity/n_Pyycyyc_L;
  float k_Pyycyyc_R = w_R*xsec_Pyycyyc_R*luminosity/n_Pyycyyc_R;
  float k_Pyyuyyc_L = w_L*xsec_Pyyuyyc_L*luminosity/n_Pyyuyyc_L;
  float k_Pyyuyyc_R = w_R*xsec_Pyyuyyc_R*luminosity/n_Pyyuyyc_R;
  float k_Pyycyyu_L = w_L*xsec_Pyycyyu_L*luminosity/n_Pyycyyu_L;
  float k_Pyycyyu_R = w_R*xsec_Pyycyyu_R*luminosity/n_Pyycyyu_R;
  TH1F* h_sum = new TH1F( k_Pyyuyyu_L*h_Pyyuyyu_L + k_Pyyuyyu_R*h_Pyyuyyu_R +
                          k_Pyycyyc_L*h_Pyycyyc_L + k_Pyycyyc_R*h_Pyycyyc_R +
                          k_Pyyuyyc_L*h_Pyyuyyc_L + k_Pyyuyyc_R*h_Pyyuyyc_R +
                          k_Pyycyyu_L*h_Pyycyyu_L + k_Pyycyyu_R*h_Pyycyyu_R 
                        );
  return h_sum;
  //h_em8pp3->Add(h_Pyyuyyu_R);
  //h_R->Draw("same");
}

void test()
{
  float epol1 = -0.8;
  float ppol1 = +0.3;
  float epol2 = +0.8;
  float ppol2 = -0.3;
  TH1F* h_em8pp3 = getHist(epol1,ppol1);
  h_em8pp3->SetMinimum(0.);
  h_em8pp3->SetLineColor(2);
  h_em8pp3->Draw();
  TH1F* h_ep8pm3 = getHist(epol2,ppol2);
  h_ep8pm3->SetLineColor(4);
  h_ep8pm3->Draw("same");

  // AFB computation
  // for h_em8pp3
  double toterr,nerr,perr;
  float total    = h_em8pp3->IntegralAndError(0,h_em8pp3->GetNbinsX(),toterr);
  float negative = h_em8pp3->IntegralAndError(0,h_em8pp3->GetNbinsX()/2,nerr);
  float positive = h_em8pp3->IntegralAndError(h_em8pp3->GetNbinsX()/2,h_em8pp3->GetNbinsX(),perr);
  float afb_em8pp3    = (positive-negative)/total;
  float afberr_em8pp3 = TMath::Sqrt(TMath::Power((positive-negative)/(total*total)*toterr,2)+
                                    TMath::Power((1/total)*perr,2)+
                                    TMath::Power((1/total)*nerr,2));
  //cerr << total << " " << positive << " " << negative << endl;
  // for h_ep8pm3
  total    = h_ep8pm3->IntegralAndError(0,h_ep8pm3->GetNbinsX(),toterr);
  negative = h_ep8pm3->IntegralAndError(0,h_ep8pm3->GetNbinsX()/2,nerr);
  positive = h_ep8pm3->IntegralAndError(h_ep8pm3->GetNbinsX()/2,h_ep8pm3->GetNbinsX(),perr);
  //cerr << total << " " << positive << " " << negative << endl;
  float afb_ep8pm3    = (positive-negative)/total;
  float afberr_ep8pm3 = TMath::Sqrt(TMath::Power((positive-negative)/(total*total)*toterr,2)+
                                    TMath::Power((1/total)*perr,2)+
                                    TMath::Power((1/total)*nerr,2));

  TPaveText* pt = new TPaveText(0.2,0.65,0.58,0.85,"brNDC");
  stringstream label_em8pp3, label_ep8pm3;
  label_em8pp3 << fixed << setprecision(1) << "AFB (e" << epol1 << ", p" << ppol1 << ") = " << setprecision(4) << afb_em8pp3 << " #pm " << afberr_em8pp3 << ends;
  label_ep8pm3 << fixed << setprecision(1) << "AFB (e" << epol2 << ", p" << ppol2 << ") = " << setprecision(4) << afb_ep8pm3 << " #pm " << afberr_ep8pm3 << ends;
  pt->AddText(label_em8pp3.str().data());
  ((TText*)(pt->GetListOfLines()->Last()))->SetTextColor(2);
  pt->AddText(label_ep8pm3.str().data());
  ((TText*)(pt->GetListOfLines()->Last()))->SetTextColor(4);
  pt->SetFillStyle(0);
  pt->SetBorderSize(1);
  pt->SetLineColor(0);
  pt->Draw();
}
