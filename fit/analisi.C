Double_t DoppiaFen(Double_t* x, Double_t* par) {
  Float_t xx = x[0];
  Double_t sin_th = std::sin(std::atan((xx - par[0]) / par[1]));
  Double_t interferenza = par[2] * std::cos(0.5 * par[3] * par[4] * sin_th) *
                          std::cos(0.5 * par[3] * par[4] * sin_th);
  Double_t rifrazione = std::sin(0.5 * par[3] * par[5] * sin_th) /
                        (0.5 * par[3] * par[5] * sin_th);
  return interferenza * rifrazione * rifrazione + par[6];
}
Double_t Interferenza(Double_t* x, Double_t* par) {
	Float_t xx = x[0];
  Double_t sin_th = std::sin(std::atan((xx - par[0]) / par[1]));
  Double_t interferenza = par[2] * std::cos(0.5 * par[3] * par[4] * sin_th) *
                          std::cos(0.5 * par[3] * par[4] * sin_th);
  return interferenza + par[5];
}
Double_t Diffrazione(Double_t* x, Double_t* par) {
	Float_t xx = x[0];
  Double_t sin_th = std::sin(std::atan((xx - par[0]) / par[1]));
  Double_t rifrazione = std::sin(0.5 * par[3] * par[4] * sin_th) /
                        (0.5 * par[3] * par[4] * sin_th);
  return par[2] * rifrazione * rifrazione + par[5];
}

void analisi(const char* file) {
  gROOT->SetStyle("Modern");
  gStyle->SetOptFit(0);

  double const k = 2 * TMath::Pi() / 632.8e-9;

  TGraphErrors* graph = new TGraphErrors(file, "%lg %lg %lg");
  graph->SetTitle(
      "Diffrazione e intereferenza da doppia fenditura;Posizione "
      "(m);Intensit#grave{a} (un. arb.)");
  graph->GetXaxis()->SetLabelSize(0.028);
  graph->GetXaxis()->SetTitleSize(0.028);
  graph->GetXaxis()->SetTitleOffset(1.2);
  graph->GetXaxis()->SetNdivisions(513, kTRUE);
  graph->GetXaxis()->SetDecimals();
  graph->GetYaxis()->SetLabelSize(0.028);
  graph->GetYaxis()->SetTitleSize(0.028);
  graph->SetMarkerStyle(2);
  graph->SetMarkerColor(12);
  graph->SetLineColor(12);
  graph->GetHistogram()->SetMinimum(0);

  int N = graph->GetN();
  for (int i = 0; i < N; i++) {
    // double new_err = graph->GetErrorY(i) + 3.68124;
    double new_err = graph->GetErrorY(i) + 0.27;
    // double new_err = graph->GetErrorY(i) * 1.004;
    graph->SetPointError(i, 0.000025, new_err);
  }
  double x_min = graph->GetPointX(0);
  double x_max = graph->GetPointX(N - 1);

  TF1* func = new TF1("func", "DoppiaFen", x_min, x_max, 7);
  func->SetLineWidth(1);
  func->SetLineColor(46);
  func->SetParNames("Center (c)", "Slits/PhCell Distance (L)",
                    "Max. Intensity (I)", "Wavenumber (k)",
                    "Slits Distance (D)", "Slits Width (d)", "Offset");
  func->SetParameters(0.0601618, 1.1025, 94, k, 0.00100, 0.00015, 1);
  // func->SetParLimits(3, k, k);

  //graph->Fit("func", "Q0R" /*, "", 0.053, 0.068*/);

  TF1* interf = new TF1("interf", "Interferenza", x_min, x_max, 6);
  interf->SetLineColor(30);
  //interf->SetLineStyle(kDashed);
  interf->SetLineWidth(1);
  interf->SetParameters(func->GetParameter(0), func->GetParameter(1), func->GetParameter(2), func->GetParameter(3), func->GetParameter(4), func->GetParameter(6));

  TF1* diffr = new TF1("diffr", "Diffrazione", x_min, x_max, 6);
  diffr->SetLineColor(38);
  //diffr->SetLineStyle(kDashed);
  diffr->SetLineWidth(1);
  diffr->SetParameters(func->GetParameter(0), func->GetParameter(1), func->GetParameter(2), func->GetParameter(3), func->GetParameter(5), func->GetParameter(6));

  double x_vec[N];
  double y_vec1[N];
  for (int i = 0; i < N; i++) {
    x_vec[i] = graph->GetPointX(i);
    y_vec1[i] = graph->GetPointY(i) - func->Eval(graph->GetPointX(i));
  }

  TGraph* res = new TGraph(N, x_vec, y_vec1);
  res->SetLineWidth(1);
  res->SetLineColor(46);
  res->SetTitle(";Posizione (m);");
  res->GetXaxis()->SetLabelSize(0.08);
  res->GetXaxis()->SetTitleSize(0.1);
  res->GetXaxis()->SetTitleOffset(-0.2);
  res->GetXaxis()->SetLimits(graph->GetXaxis()->GetXmin(),
                             graph->GetXaxis()->GetXmax());
  res->GetYaxis()->SetLabelSize(0.08);

  TCanvas* canva = new TCanvas("", "", 800, 600);
  TLegend* leg = new TLegend(0.71, 0.79, 0.89, 0.89);
  leg->AddEntry(graph, "Dati Sperimentali", "pe");
  leg->AddEntry(func, "Valori attesi", "l");
  leg->AddEntry(interf, "Interferenza", "l");
  //leg->AddEntry(diffr, "Diffrazione", "l");
  // canva->Divide(1,2);
  // canva->cd(1)->SetPad(0, 0.25, 1, 1);
  // canva->cd(2)->SetPad(0, 0, 1, 0.25);
  // canva->cd(1)->SetBottomMargin(0.01);
  // canva->cd(2)->SetTopMargin(-0);
  // canva->cd(2)->SetBottomMargin(0.18);
  //canva->/*cd(1)->*/ SetGrid();
  graph->Draw("APLE");
  interf->Draw("SAME");
  //diffr->Draw("SAME");
  func->Draw("SAME");
  leg->Draw("SAME");
  // canva->cd(2)->SetGrid();
  // res->Draw("AL");

  std::cout << std::setfill('=') << std::setw(76) << '\n';
  std::cout << std::setfill('.') << std::setw(37) << "SINGLE FIT"
            << std::setw(39) << '\n';
  std::cout << std::setfill('.') << "DOF: " << std::setw(70) << func->GetNDF()
            << '\n';
  std::cout << std::setfill('.') << "Chi Square: " << std::setw(63)
            << func->GetChisquare() << '\n';
  std::cout << std::setfill('.') << "Reduced Chi Square: " << std::setw(55)
            << func->GetChisquare() / func->GetNDF() << '\n';
  std::cout << std::setfill('.') << "Probability: " << std::setw(62)
            << func->GetProb() << '\n';
  std::cout << std::setfill('-') << std::setw(76) << '\n';
  std::string p_tot[7];
  for (int i = 0; i < 7; i++) {
    std::string p = std::to_string(func->GetParameter(i)) + " +/- ";
    p_tot[i] = p + std::to_string(func->GetParError(i));
  }
  std::cout << std::setfill('.') << func->GetParName(0) << std::setw(65)
            << p_tot[0] << '\n';
  std::cout << std::setfill('.') << func->GetParName(1) << std::setw(50)
            << p_tot[1] << '\n';
  std::cout << std::setfill('.') << func->GetParName(2) << std::setw(57)
            << p_tot[2] << '\n';
  std::cout << std::setfill('.') << func->GetParName(3) << std::setw(61)
            << p_tot[3] << '\n';
  std::cout << std::setfill('.') << func->GetParName(4) << std::setw(57)
            << p_tot[4] << '\n';
  std::cout << std::setfill('.') << func->GetParName(5) << std::setw(60)
            << p_tot[5] << '\n';
  std::cout << std::setfill('.') << func->GetParName(6) << std::setw(69)
            << p_tot[6] << '\n';

  std::cout << std::setfill('=') << std::setw(76) << '\n';

  double center = func->GetParameter(0);
  double interval_d = 0.0606;
  double interval_u = 0.0611;
  double max_position1 = func->GetMaximumX(interval_d, interval_u);
  double min_position1 = func->GetMinimumX(center, max_position1);
  double max_position2 =
      func->GetMaximumX(interval_u, graph->GetXaxis()->GetXmax());
  double min_position2 = func->GetMinimumX(max_position1, max_position2);
  double max1 = func->GetMaximum(interval_d, interval_u);
  double max2 = func->GetMaximum(interval_u, graph->GetXaxis()->GetXmax());

  std::cout << "Massimo n.1\t" << max_position1 - center << "\t" << max1
            << '\n';
  std::cout << "Minimo n.1 \t" << min_position1 - center << '\n';
  std::cout << "Massimo n.2\t" << max_position2 - center << "\t" << max2
            << '\n';
  std::cout << "Minimo n.2 \t" << min_position2 - center << '\n';
}
