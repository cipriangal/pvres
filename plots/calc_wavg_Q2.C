//
#include <TMath.h>
#include <TString.h>
#include <iomanip.h>

void calc_wavg_Q2()
{

  gROOT->Reset();

  Double_t W[20],Q2[20],Apv[20],ApvErr[20],rate[20],Abeam[20],AbeamErr[20],AbeamRErr[20];
  Double_t sW[20],sQ2[20],sApv[20],sApvErr[20],srate[20],sAbeam[20],sAbeamErr[20],sAbeamRErr[20];
  Double_t Wavg[20],Apvavg[20],ApvErravg[20],Abeamavg[20],AbeamErravg[20];
  Double_t AbeamAerr[20],werr[20];
  Int_t jj,kk;
  Int_t n(0);
  Int_t debug(0);

  TString file1,file2,temp;

  //   example of macro to read data from an ascii file and
  //   create a root file with an histogram and an ntuple.
  
  
  //c1 = new TCanvas("c1","A plot",550,100,650,500);
  c1 = new TCanvas("c1","A plot",5,5,600,600);

  c1->SetFillColor(10);
  c1->SetGrid(0,0);
  c1->Clear();
  c1->DrawFrame(0.9,0.0,2.5,150.0);
  //c1->Divide(1,2);
  c1->cd(1);
  
  file1 = "datafiles/hms_errors_k4.txt";
  file2 = "datafiles/shms_errors_k3.txt";
  ifstream  infile(file1.Data());
  ofstream outfile(Form("datafiles/prot_Q2_avg_errors.txt"),ios::out); 

  jj=0;
  infile >> temp;
  infile >> temp;
  infile >> temp;
  infile >> temp;
  infile >> temp;
  infile >> temp;
  infile >> temp;
  while (!infile.eof()) {
   
    infile >> W[jj];
    infile >> Q2[jj];
    infile >> Apv[jj];
    infile >> ApvErr[jj];
    infile >> rate[jj];
    infile >> Abeam[jj];
    infile >> AbeamRErr[jj];

    AbeamErr[jj]   = AbeamRErr[jj]/100.*Abeam[jj];

    //cout << W[jj] << "\t" << Q2[jj] << "\t" << AbeamErr[jj] << endl;

    /*Apv[jj]      /= Q2[jj];
      ApvErr[jj]   /= Q2[jj];*/
    Abeam[jj]    /= Q2[jj];
    AbeamErr[jj] /= Q2[jj];

    jj++;
  }
  infile.close();
  jj--;
  cout << jj << "\t" << W[0] << "\t" << W[jj-1] << endl;

  infile.open(file2.Data());

  kk=0;
  infile >> temp;
  infile >> temp;
  infile >> temp;
  infile >> temp;
  infile >> temp;
  infile >> temp;
  infile >> temp;
  while (!infile.eof()) {
   
    infile >> sW[kk];
    infile >> sQ2[kk];
    infile >> sApv[kk];
    infile >> sApvErr[kk];
    infile >> srate[kk];
    infile >> sAbeam[kk];
    infile >> sAbeamRErr[kk];
    werr[kk] = 0.0;

    sAbeamErr[kk]   = sAbeamRErr[kk]/100.*sAbeam[kk];

    //cout << sW[kk] << "\t" << sQ2[kk] << "\t" << sAbeamErr[kk] << endl;

    /*sApv[kk]      /= sQ2[kk];
      sApvErr[kk]   /= sQ2[kk];*/
    sAbeam[kk]    /= sQ2[kk];
    sAbeamErr[kk] /= sQ2[kk];

    kk++;
  }
  infile.close();
  kk--;
  cout << kk << "\t" << sW[0] << "\t" << sW[kk-1] << endl << endl;

  for (int ii=0; ii<kk; ii++) {
    if (sW[ii] == W[ii]) {
      Wavg[n]      = sW[ii];
      Apvavg[n]    = (sApv[ii]/TMath::Power(sApvErr[ii],2.) + Apv[ii]/TMath::Power(ApvErr[ii],2.))
	/(1./TMath::Power(sApvErr[ii],2.) + 1./TMath::Power(ApvErr[ii],2.))*1e6;
      ApvErravg[n] = TMath::Sqrt(1./(1./TMath::Power(sApvErr[ii],2.) + 1./TMath::Power(ApvErr[ii],2.)));
      Abeamavg[n]  = -(sAbeam[ii]/TMath::Power(sAbeamErr[ii],2.) + Abeam[ii]/TMath::Power(AbeamErr[ii],2.))
	/(1./TMath::Power(sAbeamErr[ii],2.) + 1./TMath::Power(AbeamErr[ii],2.))*1e6;
      AbeamErravg[n] = TMath::Sqrt(1./(1./TMath::Power(sAbeamErr[ii],2.) 
				       + 1./TMath::Power(AbeamErr[ii],2.)))*1e6;
      AbeamAerr[n]   = 100.*AbeamErravg[n]/Abeamavg[n];
      
      cout << Wavg[n] << "\t" << Abeamavg[n] << "\t" <<  AbeamErravg[n] << "\t" << AbeamAerr[n] << endl;
      outfile << Wavg[n] << "\t" << Apvavg[n] << "\t" <<  ApvErravg[n] << "\t" << Abeamavg[n] 
	      << "\t" << AbeamErravg[n] << "\t" << AbeamAerr[n] << endl;
      n++;
    } else {
      Wavg[n]        = sW[ii];
      Apvavg[n]      = sApv[ii]*1e6;
      ApvErravg[n]   = sApvErr[ii];
      Abeamavg[n]    = sAbeam[ii]*1e6;
      AbeamErravg[n] = sAbeamErr[ii];
      AbeamAerr[n]   = 100.*AbeamErravg[n]/Abeamavg[n];
      cout << Wavg[n] << "\t" << Abeamavg[n] << "\t" << AbeamErravg[n] << "\t" << AbeamAerr[n] << endl;
      outfile << Wavg[n] << "\t" << Apvavg[n] << "\t" <<  ApvErravg[n] << "\t" << Abeamavg[n] 
	      << "\t" << AbeamErravg[n] << "\t" << AbeamAerr[n] << endl;
      n++;
    }
  }
  cout << endl << n << endl;
  outfile.close();
  //return;
  
  //TF1 *fitf = new TF1("fitf",fit.Data(),1490,2755);
  //TF1 *fitf = new TF1("fitf",fit.Data(),2896,2903);
  //TF1 *f1 = new TF1("f1","[0]*x+[1]*x*x",-0.1,0.11);
  //f1->SetParameters(23.0,0.0229,0.00626);

  gr1 = new TGraphErrors(n,Wavg,Abeamavg,werr,AbeamErravg);  
  gr1->SetTitle("Deuteron (HMS and SHMS)");
  gr1->SetFillColor(10);
  gr1->SetMarkerColor(1);
  gr1->SetLineColor(1);
  gr1->SetLineWidth(3);
  gr1->SetMarkerStyle(20);
  gr1->SetMarkerSize(2.0);

  gr1->GetXaxis()->SetTitle("W (GeV)");
  gr1->GetXaxis()->SetLabelSize(0.029);
  gr1->GetXaxis()->SetTickLength(-0.015);
  gr1->GetXaxis()->SetLabelOffset(0.01);
  gr1->GetYaxis()->SetTitle("A_{PV}(ppm)/Q^{2}");
  gr1->GetYaxis()->SetTitleOffset(1.05);
  gr1->GetXaxis()->SetTitleSize(0.04);
  gr1->GetXaxis()->CenterTitle();
  gr1->GetYaxis()->SetLabelSize(0.029);
  gr1->GetYaxis()->SetTitleSize(0.04);
  gr1->GetYaxis()->SetLabelOffset(0.01);
  gr1->GetYaxis()->SetTickLength(-0.015);
  gr1->GetYaxis()->CenterTitle();

  //gr1->Fit("fitf","rem");
  gr1->GetXaxis()->SetRangeUser(1.1,2.4);
  gr1->GetYaxis()->SetRangeUser(-120,-40);
  gr1->Draw("AP");

  return;

  Double_t par[4],epar[4],ycalc[6],diff(0.),sum(0.),sdev(0.);
  fitf->GetParameters(&par[0]);

  outfile << "Yield Fit Parameters:\n";
  for(int i=0; i<3; i++) {
    epar[i] = fitf->GetParError(i);
    if (debug) {
      cout << setprecision(8) << par[i]<< "\t";
      cout << setprecision(5) << epar[i] << endl;
    }
    outfile << Form("P%d\t",i) << setprecision(8) << par[i]<< "\t";
    outfile << setprecision(5) << epar[i] << endl;
  }
  outfile << endl;
  outfile << "Nubin \t Y-Ratio  \t Yield_calc: \n";
      
  for (int i=0; i<6; i++) {
    for(int j=0; j<3; j++) {
      ycalc[i] += par[j]*TMath::Power(nu[i],float(j)); 
    }
    if(debug) {
      cout << nu[i] << "\t";
      cout << setprecision(8) << yield[i]/ycalc[i] << "\t" << ycalc[i] << endl;
    }
    outfile << nu[i] << "\t";
    outfile << setprecision(8) << yield[i]/ycalc[i] << "\t" << ycalc[i] << endl;
    diff =  yield[i] - ycalc[i];
    sum += TMath::Power(diff,2.);
  }
  sdev = TMath::Sqrt(sum/5.);
  cout << endl << "Fit Standard Deviation:\t" << sum << "\t" 
       <<  sdev << "\t" << sdev/TMath::Sqrt(6.) << endl;
  outfile << endl << "Fit Standard Deviation:\t" << sum << "\t" 
       <<  sdev << "\t" << sdev/TMath::Sqrt(6.) << endl;

  outfile.close();
  
  
  return;
  
  TMultiGraph *mg = new TMultiGraph();
  mg->SetTitle("4.209 GeV, 3.16 GeV/c");

  mg->Add(gr1);
  mg->Draw("AP");
  gr1->Fit("fitf","R");
  mg->GetXaxis()->SetTitle("Run number");
  mg->GetXaxis()->SetLabelSize(0.029);
  mg->GetXaxis()->SetTickLength(-0.015);
  mg->GetXaxis()->SetLabelOffset(0.01);
  mg->GetYaxis()->SetTitle("Yield");
  mg->GetYaxis()->SetTitleOffset(1.2);
  mg->GetXaxis()->CenterTitle();
  mg->GetYaxis()->SetLabelSize(0.029);
  mg->GetYaxis()->SetLabelOffset(0.01);
  mg->GetYaxis()->SetTickLength(-0.015);
  mg->GetYaxis()->CenterTitle();

  /*TLine *a1 = new TLine(2841.5,228.96,2841.5,232.28);
  TLine *a2 = new TLine(2838.5,228.96,2838.5,232.28);
  TLine *b1 = new TLine(2848.5,228.96,2848.5,231.65);
  TLine *b2 = new TLine(2844.5,228.96,2844.5,232.28);
  TLine *c2 = new TLine(2850.5,228.96,2850.5,231.65);

  a1->SetLineColor(2);
  a1->SetLineWidth(2);
  a1->Draw();
  a2->SetLineColor(8);
  a2->SetLineWidth(2);
  a2->SetLineStyle(2);
  a2->Draw();
  b1->SetLineColor(2);
  b1->SetLineWidth(2);
  b1->Draw();
  b2->SetLineColor(8);
  b2->SetLineWidth(2);
  b2->SetLineStyle(2);
  b2->Draw();
  c2->SetLineColor(8);
  c2->SetLineWidth(2);
  c2->SetLineStyle(2);
  c2->Draw();*/

  c1->Update();

}













