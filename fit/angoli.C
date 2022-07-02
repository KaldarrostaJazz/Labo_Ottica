
  double const pi = 3.141592653589793;

Double_t gausBeam(Double_t* x, Double_t* par) {
  Float_t xx = x[0];
  Double_t z = (xx - par[1]) / par[2];
  return par[1] * exp(-2 * z * z);
}
Double_t cosSquare(Double_t* x, Double_t* par) {
  Float_t xx = x[0];
  return par[0] * std::cos(xx * pi / 180 + par[1]) *
         std::cos(xx * pi / 180 + par[1]);
}
Double_t d_gaus(Double_t* x, Double_t* par) {
  Float_t xx = x[0];
  Double_t z1 = (xx - par[1]) / par[2];
  Double_t z2 = (xx - par[4]) / par[5];
  return par[0] * exp(-0.5 * z1 * z1) + par[3] * exp(-0.5 * z2 * z2) + par[6];
}
Double_t bueno(Double_t* x, Double_t* par) {
  Float_t xx = x[0];
  Double_t f = std::cos(xx + par[0]) * std::cos(xx + par[0]) / (std::cos(xx) * std::cos(xx));
  Double_t ac = std::acos(std::sqrt(f * std::cos(xx) * std::cos(xx))) - xx;
  return ac;
}

void angoli() {
  gROOT->SetStyle("Modern");
  gStyle->SetOptFit(0);


  const char* titles[24]{
      "Filtro polarizzante a 0#circ;Distanza (m);Intensit#grave{a} (un. arb.)",
      "Filtro polarizzante a 10#circ;Distanza (m);Intensit#grave{a} (un. arb.)",
      "Filtro polarizzante a 20#circ;Distanza (m);Intensit#grave{a} (un. arb.)",
      "Filtro polarizzante a 30#circ;Distanza (m);Intensit#grave{a} (un. arb.)",
      "Filtro polarizzante a 36#circ;Distanza (m);Intensit#grave{a} (un. arb.)",
      "Filtro polarizzante a 38#circ;Distanza (m);Intensit#grave{a} (un. arb.)",
      "Filtro polarizzante a 39#circ;Distanza (m);Intensit#grave{a} (un. arb.)",
      "Filtro polarizzante a 40#circ;Distanza (m);Intensit#grave{a} (un. arb.)",
      "Filtro polarizzante a 42#circ;Distanza (m);Intensit#grave{a} (un. arb.)",
      "Filtro polarizzante a 50#circ;Distanza (m);Intensit#grave{a} (un. arb.)",
      "Filtro polarizzante a 60#circ;Distanza (m);Intensit#grave{a} (un. arb.)",
      "Filtro polarizzante a 70#circ;Distanza (m);Intensit#grave{a} (un. arb.)",
      "Filtro polarizzante a 80#circ;Distanza (m);Intensit#grave{a} (un. arb.)",
      "Filtro polarizzante a 90#circ;Distanza (m);Intensit#grave{a} (un. arb.)",
      "Filtro polarizzante a 100#circ;Distanza (m);Intensit#grave{a} (un. "
      "arb.)",
      "Filtro polarizzante a 110#circ;Distanza (m);Intensit#grave{a} (un. "
      "arb.)",
      "Filtro polarizzante a 120#circ;Distanza (m);Intensit#grave{a} (un. "
      "arb.)",
      "Filtro polarizzante a 128#circ;Distanza (m);Intensit#grave{a} (un. "
      "arb.)",
      "Filtro polarizzante a 130#circ;Distanza (m);Intensit#grave{a} (un. "
      "arb.)",
      "Filtro polarizzante a 140#circ;Distanza (m);Intensit#grave{a} (un. "
      "arb.)",
      "Filtro polarizzante a 150#circ;Distanza (m);Intensit#grave{a} (un. "
      "arb.)",
      "Filtro polarizzante a 160#circ;Distanza (m);Intensit#grave{a} (un. "
      "arb.)",
      "Filtro polarizzante a 170#circ;Distanza (m);Intensit#grave{a} (un. "
      "arb.)",
      "Filtro polarizzante a 180#circ;Distanza (m);Intensit#grave{a} (un. "
      "arb.)"};
  const char* filenames[24]{"../osservazioni/angoli/g1_68-72_0ang",
                            "../osservazioni/angoli/g1_68-72_10ang",
                            "../osservazioni/angoli/g1_68-72_20ang",
                            "../osservazioni/angoli/g10_68-72_30ang",
                            "../osservazioni/angoli/g100_68-72_36ang",
                            "../osservazioni/angoli/g100_68-72_38ang",
                            "../osservazioni/angoli/g100_68-72_39ang",
                            "../osservazioni/angoli/g100_68-72_40ang",
                            "../osservazioni/angoli/g100_68-72_42ang",
                            "../osservazioni/angoli/g10_68-72_50ang",
                            "../osservazioni/angoli/g10_68-72_60ang",
                            "../osservazioni/angoli/g10_68-72_70ang",
                            "../osservazioni/angoli/g1_68-72_80ang",
                            "../osservazioni/angoli/g1_68-72_90ang",
                            "../osservazioni/angoli/g1_68-72_100ang",
                            "../osservazioni/angoli/g1_68-72_110ang",
                            "../osservazioni/angoli/g1_68-72_120ang",
                            "../osservazioni/angoli/g1_68-72_128ang",
                            "../osservazioni/angoli/g1_68-72_130ang",
                            "../osservazioni/angoli/g1_68-72_140ang",
                            "../osservazioni/angoli/g1_68-72_150ang",
                            "../osservazioni/angoli/g1_68-72_160ang",
                            "../osservazioni/angoli/g1_68-72_170ang",
                            "../osservazioni/angoli/g1_68-72_180ang"};
  TGraphErrors* grafici[24];
  TF1* funcs[24];
  TF1* fitFuncs[24];
  fitFuncs[0] = new TF1("fitFunc0", "d_gaus", 0.067, 0.073, 7);
  // fitFuncs[0]->SetParameters(75, 0.005, 0.0708, 75, 0.005, 0.0708, 0);

  fitFuncs[1] = new TF1("fitFunc1", "d_gaus", 0.067, 0.073, 7);
  // fitFuncs[1]->SetParameters(45, 0.005, 0.0708, 45, 0.005, 0.0708, 0);

  fitFuncs[2] = new TF1("fitFunc2", "d_gaus", 0.067, 0.073, 7);
  // fitFuncs[2]->SetParameters(20, 0.005, 0.0708, 20, 0.005, 0.0708, 0);

  fitFuncs[3] = new TF1("fitFunc3", "d_gaus", 0.067, 0.073, 7);
  // fitFuncs[3]->SetParameters(35, 0.005, 0.0708, 35, 0.005, 0.0708, 0.5);

  fitFuncs[4] = new TF1("fitFunc4", "d_gaus", 0.067, 0.073, 7);
  // fitFuncs[4]->SetParameters(35, 0.005, 0.0704, 35, 0.005, 0.0709, 0.5);

  fitFuncs[5] = new TF1("fitFunc5", "d_gaus", 0.067, 0.073, 7);
  // fitFuncs[5]->SetParameters(60, 0.005, 0.0706, 60, 0.005, 0.0708, 5);

  fitFuncs[6] = new TF1("fitFunc6", "d_gaus", 0.067, 0.073, 7);
  // fitFuncs[6]->SetParameters(30, 0.002, 0.0704, 36, 0.002, 0.0709, 4.5);

  fitFuncs[7] = new TF1("fitFunc7", "d_gaus", 0.067, 0.073, 7);
  // fitFuncs[7]->SetParameters(37, 0.0025, 0.0704, 51, 0.0025, 0.0709, 4.5);

  fitFuncs[8] = new TF1("fitFunc8", "d_gaus", 0.067, 0.073, 7);
  // fitFuncs[8]->SetParameters(32, 0.0025, 0.0704, 50, 0.0025, 0.0709, 4.5);

  fitFuncs[9] = new TF1("fitFunc9", "d_gaus", 0.067, 0.073, 7);
  // fitFuncs[9]->SetParameters(35, 0.005, 0.07, 35, 0.005, 0.071, 5);

  fitFuncs[10] = new TF1("fitFunc10", "d_gaus", 0.067, 0.073, 7);
  // fitFuncs[10]->SetParameters(35, 0.005, 0.07, 35, 0.005, 0.071, 0.5);

  fitFuncs[11] = new TF1("fitFunc11", "d_gaus", 0.067, 0.073, 7);
  // fitFuncs[11]->SetParameters(35, 0.005, 0.07, 35, 0.005, 0.071, 0.5);

  fitFuncs[12] = new TF1("fitFunc12", "d_gaus", 0.067, 0.073, 7);
  // fitFuncs[12]->SetParameters(35, 0.005, 0.07, 35, 0.005, 0.071, 0.5);

  fitFuncs[13] = new TF1("fitFunc13", "d_gaus", 0.067, 0.073, 7);
  // fitFuncs[13]->SetParameters(35, 0.005, 0.07, 35, 0.005, 0.071, 0.5);

  fitFuncs[14] = new TF1("fitFunc14", "d_gaus", 0.067, 0.073, 7);
  // fitFuncs[14]->SetParameters(35, 0.005, 0.07, 35, 0.005, 0.071, 0.5);

  fitFuncs[15] = new TF1("fitFunc15", "d_gaus", 0.067, 0.073, 7);
  // fitFuncs[15]->SetParameters(35, 0.005, 0.07, 35, 0.005, 0.071, 0.5);

  fitFuncs[16] = new TF1("fitFunc16", "d_gaus", 0.067, 0.073, 7);
  // fitFuncs[16]->SetParameters(35, 0.005, 0.07, 35, 0.005, 0.071, 0.5);

  fitFuncs[17] = new TF1("fitFunc17", "d_gaus", 0.067, 0.073, 7);
  // fitFuncs[17]->SetParameters(35, 0.005, 0.07, 35, 0.005, 0.071, 0.5);

  fitFuncs[18] = new TF1("fitFunc18", "d_gaus", 0.067, 0.073, 7);
  // fitFuncs[18]->SetParameters(35, 0.005, 0.07, 35, 0.005, 0.071, 0.5);

  fitFuncs[19] = new TF1("fitFunc19", "d_gaus", 0.067, 0.073, 7);
  // fitFuncs[19]->SetParameters(35, 0.005, 0.07, 35, 0.005, 0.071, 0.5);

  fitFuncs[20] = new TF1("fitFunc20", "d_gaus", 0.067, 0.073, 7);
  // fitFuncs[20]->SetParameters(35, 0.005, 0.07, 35, 0.005, 0.071, 0.5);

  fitFuncs[21] = new TF1("fitFunc21", "d_gaus", 0.067, 0.073, 7);
  // fitFuncs[21]->SetParameters(35, 0.005, 0.07, 35, 0.005, 0.071, 0.5);

  fitFuncs[22] = new TF1("fitFunc22", "d_gaus", 0.067, 0.073, 7);
  // fitFuncs[22]->SetParameters(35, 0.005, 0.07, 35, 0.005, 0.071, 0.5);

  fitFuncs[23] = new TF1("fitFunc23", "d_gaus", 0.067, 0.073, 7);
  // fitFuncs[23]->SetParameters(35, 0.005, 0.07, 35, 0.005, 0.071, 0.5);

  const char* name[24] = {"fitFunc0",  "fitFunc1",  "fitFunc2",  "fitFunc3",
                          "fitFunc4",  "fitFunc5",  "fitFunc6",  "fitFunc7",
                          "fitFunc8",  "fitFunc9",  "fitFunc10", "fitFunc11",
                          "fitFunc12", "fitFunc13", "fitFunc14", "fitFunc15",
                          "fitFunc16", "fitFunc17", "fitFunc18", "fitFunc19",
                          "fitFunc20", "fitFunc21", "fitFunc22", "fitFunc23"};

  for (int i = 0; i < 24; i++) {
    grafici[i] = new TGraphErrors(filenames[i], "%lg %lg %lg");
    grafici[i]->SetTitle(titles[i]);
    grafici[i]->SetMarkerStyle(2);
    grafici[i]->SetMarkerColor(12);
    grafici[i]->Fit("gaus", "Q0+");
    TF1* gauss = grafici[i]->GetFunction("gaus");
    fitFuncs[i]->SetParameters(
        gauss->GetParameter(0) / 2, gauss->GetParameter(1),
        gauss->GetParameter(2) - 0.0005, gauss->GetParameter(0) / 2,
        gauss->GetParameter(1), gauss->GetParameter(2) + 0.0005, 3);
    grafici[i]->Fit(name[i], "Q0+");
    funcs[i] = grafici[i]->GetFunction(name[i]);
    funcs[i]->SetLineWidth(1);
    funcs[i]->SetLineColor(42);
    funcs[i]->SetParNames("Constant_1", "Mean_1", "Sigma_1", "Constant_2",
                          "Mean_2", "Sigma_2", "Offset");
  }

  double gains[23]{1, 1, 1, 10, 100, 100, /*100, */100, 100, 10, 10, 10,
                   1, 1, 1, 1,  1,   1,   1,   1,   1,   1,  1,  1};
  double intensities[23];
  double intensities_err[23];
  double centers[23];
  double centers_err[23];
  for (int i = 0; i < 24; i++) {
    if (i < 6) {
    intensities[i] = (funcs[i]->GetParameter("Constant_1") + funcs[i]->GetParameter("Constant_2")) / gains[i];
    intensities_err[i] = 0.3;
    centers[i] = (funcs[i]->GetParameter("Mean_1") + funcs[i]->GetParameter("Mean_2")) / 2;
    centers_err[i] = 2.5e-5;
    }
    else if (i == 6) {std::cout << "6" << '\n';}
    else {
      intensities[i-1] = (funcs[i]->GetParameter("Constant_1") + funcs[i]->GetParameter("Constant_2")) / gains[i-1];
    intensities_err[i-1] = 0.3;
    centers[i-1] = (funcs[i]->GetParameter("Mean_1") + funcs[i]->GetParameter("Mean_2")) / 2;
    centers_err[i-1] = 2.5e-5;
    }
  }


  double x_vec[23]{0,  10, 20,  30,  36,  38,  /*39,  */40,  42,  50,  60,  70,
                   80, 90, 100, 110, 120, 128, 130, 140, 150, 160, 170, 180};

  TGraphErrors* graph_1 =
      new TGraphErrors(23, x_vec, centers, nullptr, centers_err);
  graph_1->SetTitle(
      "Posizione centro in funzione dell'angolo di polarizzazione;Angolo "
      "(#circ);Distanza (m)");
  graph_1->RemovePoint(4);
  graph_1->RemovePoint(5);
  graph_1->RemovePoint(6);
  graph_1->RemovePoint(7);
  graph_1->SetLineColor(12);
  graph_1->SetMarkerStyle(2);
  graph_1->SetMarkerColor(kBlack);
  graph_1->Fit("pol1");
  TF1* retta = graph_1->GetFunction("pol1");
  retta->SetLineWidth(1);
  retta->SetLineColor(42);

  TGraphErrors* graph_2 =
      new TGraphErrors(23, x_vec, intensities, nullptr, intensities_err);
  graph_2->SetTitle(
      "Intensit#grave{a} massima in funzione dell'angolo di "
      "polarizzazione;Angolo "
      "(#circ);Distanza (m)");
  graph_2->SetLineColor(12);
  graph_2->SetMarkerStyle(2);
  graph_2->SetMarkerColor(kBlack);

  TF1* coseno = new TF1("coseno", "cosSquare", 0, 180, 2);
  coseno->SetLineColor(42);
  coseno->SetLineWidth(1);
  coseno->SetParameters(3.68359e+02, 8.81928e-01);

  graph_2->Fit("coseno");


std::cout << "Reduced Chi Square:\n";
std::cout << retta->GetChisquare() / retta->GetNDF() << '\n' << coseno->GetChisquare() / coseno->GetNDF() << "\n\n";
std::cout << funcs[6]->GetParameter(0) << " " << funcs[6]->GetParameter(3) << '\n';


  double rapporti[16];
  double coseni[16];
  double alpha[16];
  double ax_vec[16];
  for (int i = 0; i < 16; i++) {
    if (i < 6) {
    rapporti[i] = funcs[i]->GetParameter("Constant_2") / funcs[i]->GetParameter("Constant_1");
    coseni[i] = TMath::Cos((129.49 - x_vec[i]) * pi / 180);
    ax_vec[i] = (129.49 - x_vec[i]) * pi / 180;
    } else if ( i == 6) {std::cout << "6\n";}
    else {
      rapporti[i-1] = funcs[i]->GetParameter("Constant_2")/ funcs[i]->GetParameter("Constant_1");
      coseni[i-1] = TMath::Cos((129.49 - x_vec[i-1]) * pi / 180);
      ax_vec[i-1] = (129.49 - x_vec[i-1]) * pi / 180;
    }
  }
  for (int i = 0; i < 16; i++) {
    alpha[i] = acos(sqrt(rapporti[i] * coseni[i] * coseni[i])) - ax_vec[i];
  }

  TGraph* io = new TGraph(16, ax_vec, alpha);
  io->SetLineColor(kRed);
  io->SetMarkerStyle(kOpenCircle);

io->Fit("pol1");
io->Fit("pol2", "+");

TF1* f_1 = io->GetFunction("pol1");
TF1* f_2 = io->GetFunction("pol2");

  TCanvas* c_3 = new TCanvas();
  c_3->SetGrid();
  io->Draw("APL");
  f_1->Draw("SAME");
  f_2->Draw("SAME");
}