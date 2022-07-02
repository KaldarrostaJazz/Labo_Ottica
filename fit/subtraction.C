Double_t DoppiaFen(Double_t* x, Double_t* par) {
  Float_t xx = x[0];
  Double_t sin_th = std::sin(std::atan((xx - par[0]) / par[1]));
  Double_t interferenza = par[2] * std::cos(0.5 * par[3] * par[4] * sin_th) *
                          std::cos(0.5 * par[3] * par[4] * sin_th);
  Double_t rifrazione = std::sin(0.5 * par[3] * par[5] * sin_th) /
                        (0.5 * par[3] * par[5] * sin_th);
  return interferenza * rifrazione * rifrazione + par[6];
}

void subtraction(const char* file1, const char* file2) {
  gBenchmark->Start("Analisi dati con sottrazione");
  gROOT->SetStyle("Modern");
  gStyle->SetOptFit(0);

  double const k = 2 * TMath::Pi() / 632.8e-9;

  TGraphErrors graph1{file1, "%lg %lg %lg"};
  int N1 = graph1.GetN();
  double x1_min = graph1.GetPointX(0);
  double x1_max = graph1.GetPointX(N1 - 1);
  TGraphErrors graph2{file2, "%lg %lg %lg"};
  int N2 = graph2.GetN();
  double x2_min = graph2.GetPointX(0);
  double x2_max = graph2.GetPointX(N2 - 1);
  try {
    if (N1 != N2) {
      throw std::runtime_error(
          "Graph1 and Graph2 have different number of points.\n");
    } else if (x1_min != x2_min || x1_max != x2_max) {
      throw std::runtime_error("Graph1 and Graph2 have different ranges.\n");
    } else {
      int N = N1;
      double x_vec[N];
      double y_vec[N];
      double y_err[N];
      for (int i = 0; i < N; i++) {
        x_vec[i] = graph1.GetPointX(i);
        y_vec[i] = (graph1.GetPointY(i) + (graph2.GetPointY(i) / 10.));
        y_err[i] = (graph1.GetErrorY(i) + graph2.GetErrorY(i));
      }
      TGraphErrors* graph = new TGraphErrors(N, x_vec, y_vec, nullptr, y_err);
      graph->SetTitle(
          "Diffrazione e intereferenza da doppia fenditura;;Intensit#grave{a} "
          "(un. arb.)");
      graph->GetXaxis()->SetLabelSize(0);
      graph->GetYaxis()->SetTitleSize(0.04);
      graph->SetMarkerStyle(2);
      graph->SetMarkerColor(12);

      /*for (int i = 0; i < N; i++) {
          //double new_err = graph->GetErrorY(i) + 3.68124;
          double new_err = graph->GetErrorY(i) + 0.27;
          //double new_err = graph->GetErrorY(i) * 1.004;
          graph->SetPointError(i, graph->GetErrorX(i), new_err);
      }*/
      double x_min = graph->GetPointX(0);
      double x_max = graph->GetPointX(N - 1);

      TF1* func = new TF1("func", "DoppiaFen", x_min, x_max, 7);
      func->SetLineWidth(1);
      func->SetLineColor(46);
      func->SetParNames("Center (c)", "Slits/PhCell Distance (L)",
                        "Max. Intensity (I)", "Wavenumber (k)",
                        "Slits Distance (D)", "Slits Width (d)", "Offset");
      func->SetParameters(0.0602572, 1.4325, 161, k, 0.00025, 0.00015, 1);
      func->SetParLimits(3, k, k);

      graph->Fit("func", "Q0R", "", 0.054, 0.066);

      double y_res[N];
      for (int i = 0; i < N; i++) {
        x_vec[i] = graph->GetPointX(i);
        y_res[i] = graph->GetPointY(i) - func->Eval(graph->GetPointX(i));
      }

      TGraph* res = new TGraph(N, x_vec, y_res);
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
      leg->AddEntry(graph, "dati sperimentali", "pe");
      leg->AddEntry(func, "fit singolo", "l");
      canva->Divide(1, 2);
      canva->cd(1)->SetPad(0, 0.25, 1, 1);
      canva->cd(2)->SetPad(0, 0, 1, 0.25);
      canva->cd(1)->SetBottomMargin(0.01);
      canva->cd(2)->SetTopMargin(-0);
      canva->cd(2)->SetBottomMargin(0.18);
      canva->cd(1)->SetGrid();
      graph->Draw("APE");
      func->Draw("SAME");
      leg->Draw("SAME");
      canva->cd(2)->SetGrid();
      res->Draw("AL");

      std::cout << std::setfill('=') << std::setw(76) << '\n';
      std::cout << std::setfill('.') << std::setw(37) << "SINGLE FIT"
                << std::setw(39) << '\n';
      std::cout << std::setfill('.') << "DOF: " << std::setw(70)
                << func->GetNDF() << '\n';
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
    }
  } catch (std::exception& e) {
    std::cout << e.what();
  }
  gBenchmark->Show("Analisi dati con sottrazione");
}
