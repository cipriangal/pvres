Double_t W[20],Q2[20],Apv[20],ApvErr[20],rate[20],Abeam[20],AbeamErr[20],AbeamRErr[20];
Double_t sW[20],sQ2[20],sApv[20],sApvErr[20],srate[20],sAbeam[20],sAbeamErr[20],sAbeamRErr[20];
Double_t AbeamAerr[20],werr[20];

void calc_wavg_Q2(){

  //structre [nPoint][data types]
  // data types: 0=W ; 1=Q2 ; 2=A_PV ; 3= A_PV error ; 4=rate
  //             5=A_beam ; 6=A_beam Relative Error ; 7=A_beam Absolute Error
  double h_data_hms[20][8],h_data_shms[20][8];
  double d_data_hms[20][8],d_data_shms[20][8];

  int h_n_hms  = readData("./datafiles/hyd_hms_errors_k4.txt",h_data_hms);
  int h_n_shms = readData("./datafiles/hyd_shms_errors_k3.txt",h_data_shms);
  if(h_n_hms!=h_n_shms){
    cout<<"h # data points doesn't match "<<h_n_shms<<" "<<h_n_hms<<endl;
    return;
  }
  
  int d_n_hms  = readData("./datafiles/deut_hms_errors_k4.txt",d_data_hms);
  int d_n_shms = readData("./datafiles/deut_shms_errors_k3.txt",d_data_shms);
  if(d_n_hms!=d_n_shms){
    cout<<"d # data points doesn't match "<<d_n_shms<<" "<<d_n_hms<<endl;
    return;
  }

  TCanvas *c1=new TCanvas("c1","c1");
  drawAverage(h_n_hms,h_data_hms,h_data_shms,"Hydrogen");
  TCanvas *c2=new TCanvas("c2","c2");
  drawAverage(d_n_hms,d_data_hms,d_data_shms,"Deuterium");
}

void drawAverage(int n, double dhms[20][8],double dshms[20][8],string grNm){

  Double_t Wavg[20],Apvavg[20],ApvErravg[20],Abeamavg[20],AbeamErravg[20];

  double avg(0);
  double davg(0);
  cout << "#\tW\tA_beam\td(A_beam)\td(A_beam)/A_beam" << endl;

  for(int i=0;i<n;i++){
    if(dhms[i][0]!=dshms[i][0]){
      cout<<" W does not match: i hms shms "<<i<<"\t"<<dhms[i][0]<<"\t"<<dshms[i][0]<<endl;
    }

    Wavg[i] = dhms[i][0];

    /* Not currently used
    Apvavg[i] = ( dhms[i][2]/pow(dhms[i][3],2) + dshms[i][2]/pow(dshms[i][3],2)  ) /
      ( 1/pow(dhms[i][3],2) + 1/pow(dshms[i][3],2) ) *1e6;

    ApvErravg[i] = sqrt( 1 / ( 1/pow(dhms[i][3],2) + 1/pow(dshms[i][3],2) ) *1e6;
    */

    Abeamavg[i] =  - ( dhms[i][5]/pow(dhms[i][7],2) + dshms[i][5]/pow(dshms[i][7],2)  ) /
      ( 1/pow(dhms[i][7],2) + 1/pow(dshms[i][7],2) ) *1e6;
        
    AbeamErravg[i] =  sqrt( 1 / ( 1/pow(dhms[i][7],2) + 1/pow(dshms[i][7],2) ) ) *1e6;

    cout<< i << "\t" << Wavg[i] << "\t" << Abeamavg[i] << "\t" << AbeamErravg[i]
	<< "\t" << fabs(AbeamErravg[i]/Abeamavg[i]) << endl;
    avg += Abeamavg[i]/pow(AbeamErravg[i],2);
    davg += 1/pow(AbeamErravg[i],2);
  }

  double avg = avg/davg;
  double davg = sqrt( 1 / davg );
  cout<<"Total: avg = "<<avg<<" d(avg) = "<<davg<<" d(avg)/avg = "<<fabs(davg/avg)<<endl;
  
  double werr[20]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  double fixA[20]={-110,-110,-110,-110,-110,-110,-110,-110,-110,-110,
		   -110,-110,-110,-110,-110,-110,-110,-110,-110,-110};
  TGraphErrors *gr=new TGraphErrors(n,Wavg,Abeamavg,werr,AbeamErravg);
  //TGraphErrors *gr=new TGraphErrors(n,Wavg,fixA,werr,AbeamErravg);
  gr->SetTitle(";W [GeV];A_{PV}(ppm)/Q^{2}");
  gr->SetName(grNm.c_str());

  gr->SetMarkerColor(1);
  gr->SetMarkerSize(2);
  gr->SetMarkerStyle(20);

  gr->GetXaxis()->CenterTitle();
  gr->GetXaxis()->SetTitleOffset(1.05);
  gr->GetXaxis()->SetTitleSize(0.04);

  gr->GetYaxis()->CenterTitle();
  gr->GetYaxis()->SetTitleOffset(1.05);
  gr->GetYaxis()->SetTitleSize(0.04);
  
  gr->GetXaxis()->SetRangeUser(1.1,2.4);
  gr->GetYaxis()->SetRangeUser(-120,-40);
  gr->Draw("AP");
}

int readData(string fnm,double data[20][8]){
  ifstream fin(fnm.c_str());
  string dummy;
  std::getline(fin,dummy);
  double dm[8];
  int n=0;
  while(fin>>dm[0]>>dm[1]>>dm[2]>>dm[3]>>dm[4]>>dm[5]>>dm[6]){
    for(int i=0;i<7;i++)
      data[n][i]=dm[i];
    //absolute error = relErr / 100 * A_beam
    data[n][7] = data[n][6]/100 * data[n][5];
    //absolute error /= Q2
    data[n][7] /= data[n][1];
    //A_beam /= Q2
    data[n][5] /= data[n][1];
    n++;
  }
  fin.close();
  return n;
}








