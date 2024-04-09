#include "DrawFitResult.h"
#include "RooMath.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TLine.h"
#include <complex>
#include <iomanip>
#include <iostream>

using namespace std;

DrawFitResult::DrawFitResult(std::function<double(const double *, const double *)> background,
                             std::function<double(const double *, const double *)> signal,
                             std::function<double(const double *, const double *)> total,
                             std::vector<TGraphErrors *> graph, Fit &fit, int NSignal, int SigTyp)
    : Background(background), Signal(signal), Total(total), Graph(graph), FitRes(fit),
      _NSignal(NSignal), _SigTyp(SigTyp) {}

DrawFitResult::DrawFitResult(std::function<double(const double *, const double *)> background,
                             std::function<double(const double *, const double *)> signal,
                             std::function<double(const double *, const double *)> total,
                             TGraphErrors *graph, int NSignal, int SigTyp)
    : Background(background), Signal(signal), Total(total), GraphSingle(graph), _NSignal(NSignal),
      _SigTyp(SigTyp) {}

void DrawFitResult::DrawResult(std::vector<double> resPara, double x_first, double x_last,
                               double fit_x_low, double fit_x_high) {
  TCanvas *c1 = new TCanvas("Intensitaeten", "", 4000, 4000);
  c1->Divide(2, 2);
  std::vector<std::vector<TF1 *>> SignalFunctions;

  int NSigPara = 2 + _SigTyp; // for BW: 3, Voigt: 4

  for (int bin = 0; bin < 4; bin++) {
    std::vector<TF1 *> SignalForThisBin;
    stringstream ss0;
    ss0 << "Total_" << bin;

    TF1 *FunctionModel = new TF1(ss0.str().c_str(), Total, x_first, x_last, 5, 1);
    FunctionModel->SetNpx(1000);
    FunctionModel->SetLineColor(kRed);

    for (int i = 0; i < _NSignal; i++) {
      stringstream ss;
      ss << "f_" << i << "_" << bin;
      TF1 *FunctionSignal = new TF1(ss.str().c_str(), Signal, x_first, x_last, NSigPara, 1);
      // for (int k = 0; k < NSigPara - 1 - 1; k++) {
      //   // FunctionModel->SetParameter(k, FunctionTotal->GetParameter(k + NSigPara * i));
      //   // FunctionSignal->SetParameter(k, FunctionTotal->GetParameter(k + NSigPara * i));
      // }
      if (bin == 0)
        FunctionSignal->SetParameter(0, resPara.at(0));
      else
        FunctionSignal->SetParameter(0, resPara.at(8 * bin - 2 * (bin - 1)));
      FunctionSignal->SetParameter(1, resPara.at(1));
      FunctionSignal->SetParameter(2, resPara.at(2));
      FunctionSignal->SetParameter(3, resPara.at(3 + (bin * 6))); // magnitude

      FunctionSignal->SetNpx(1000);
      FunctionSignal->SetLineColor(kAzure + 2 + 500 * i);
      FunctionSignal->SetLineStyle(2);
      SignalForThisBin.push_back(FunctionSignal);
    }
    SignalFunctions.push_back(SignalForThisBin);
    // Model
    if (bin == 0)
      FunctionModel->SetParameter(0, resPara.at(0));
    else
      FunctionModel->SetParameter(0, resPara.at(8 * bin - 2 * (bin - 1)));
    FunctionModel->SetParameter(1, resPara.at(1));
    FunctionModel->SetParameter(2, resPara.at(2));
    FunctionModel->SetParameter(3, resPara.at(3 + (bin * 6))); // magnitude
    FunctionModel->SetParameter(4, resPara.at(4 + (bin * 6)));
    // FunctionModel->SetParameter(5, resPara.at(5 + (bin * 6)));
    // FunctionModel->SetParameter(6, resPara.at(6 + (bin * 6)));
    // FunctionModel->SetParameter(7, resPara.at(7 + (bin * 6)));

    FunctionBackground = new TF1("Background", Background, x_first, x_last, 1, 1);
    // Background
    FunctionBackground->SetParameter(0, resPara.at(4 + (bin * 6)));
    // FunctionBackground->SetParameter(1, resPara.at(5 + (bin * 6)));
    // FunctionBackground->SetParameter(2, resPara.at(6 + (bin * 6)));
    // FunctionBackground->SetParameter(3, resPara.at(7 + (bin * 6)));

    FunctionBackground->SetNpx(1000);
    FunctionBackground->SetLineColor(kOrange + 7);
    FunctionBackground->SetLineStyle(2);

    // Legend
    TLegend *leg = new TLegend(0.15, 0.48, 0.45, 0.88);
    leg->SetTextFont(132);
    leg->AddEntry(Graph[bin], "Datenpunkte", "p");
    leg->AddEntry(FunctionModel, "Total", "l");

    for (int i = 0; i < _NSignal; i++) {
      std::stringstream ss;
      ss << "FunctionSignal " << i + 1;
      leg->AddEntry(SignalFunctions.at(bin).at(i), ss.str().c_str(), "l");

      if (_SigTyp == 1) { // BW
        std::stringstream ss1;
        ss1 << "m_{0} " << i + 1 << "= " << std::setprecision(4) << std::fixed
            << FitRes.GetResultValue().at(0 + 4 * i) << " #pm "
            << FitRes.GetResultValue().at(1 + 4 * i) << " GeV/#it{c}^{2}";
        stringstream ss2;
        ss2 << "#Gamma " << i + 1 << "= " << std::setprecision(4) << std::fixed
            << FitRes.GetResultValue().at(2 + 4 * i) << " #pm "
            << FitRes.GetResultValue().at(3 + 4 * i) << " GeV/#it{c}^{2}";

        leg->AddEntry((TObject *)0, ss1.str().c_str(), "");
        leg->AddEntry((TObject *)0, ss2.str().c_str(), "");
      } else if (_SigTyp == 2) { // Voigt
        std::stringstream ss1;
        ss1 << "m_{0} " << i + 1 << "= " << std::setprecision(4) << std::fixed
            << FitRes.GetResultValue().at(0 + 6 * i) << " #pm "
            << FitRes.GetResultValue().at(1 + 6 * i) << " GeV/#it{c}^{2}";

        stringstream ss2;
        ss2 << "#Gamma " << i + 1 << "= " << std::setprecision(4) << std::fixed
            << FitRes.GetResultValue().at(2 + 6 * i) << " #pm "
            << FitRes.GetResultValue().at(3 + 6 * i) << " GeV/#it{c}^{2}";

        stringstream ss3;
        ss3 << "#sigma " << i + 1 << "= " << std::setprecision(4) << std::fixed
            << FitRes.GetResultValue().at(4 + bin * 2) << " #pm "
            << FitRes.GetResultValue().at(5 + bin * 2) << " GeV/#it{c}^{2}";
        // ss3 << "#sigma " << i + 1 << "= " << std::setprecision(4) << std::fixed
        //     << FitRes.GetResultValue().at(4 + bin * 2) << " #pm "
        //     << FitRes.GetResultValue().at(5 + bin * 2) << " GeV/#it{c}^{2}";

        leg->AddEntry((TObject *)0, ss1.str().c_str(), "");
        leg->AddEntry((TObject *)0, ss2.str().c_str(), "");
        leg->AddEntry((TObject *)0, ss3.str().c_str(), "");
      }
      leg->SetFillStyle(0);
    }

    leg->AddEntry(FunctionBackground, "FunctionBackground", "l");
    std::stringstream ss;
    ss << "#chi^{2}_{red.}=" << std::setprecision(2) << std::fixed << FitRes.GetRedChi();
    leg->AddEntry((TObject *)0, ss.str().c_str(), "");
    leg->SetFillStyle(0);

    // Draw
    FunctionModel->GetYaxis()->SetRangeUser(0, 1.1 * FunctionModel->GetMaximum(x_first, x_last));
    FunctionModel->SetTitle(";#it{m_{#pi^{#minus}#pi^{+}#eta}} [GeV/#it{c}^{2}];Yields");
    // std::cout << "Draw" << __LINE__ << endl;
    c1->cd(bin + 1);
    Graph[bin]->GetXaxis()->SetRangeUser(x_first, x_last);
    Graph[bin]->GetYaxis()->SetRangeUser(0, 1.1 * (Graph[bin]->GetYaxis()->GetXmax()));
    Graph[bin]->Draw("apz");
    FunctionModel->Draw("same");
    FunctionBackground->Draw("same");
    for (int i = 0; i < _NSignal; i++) {
      SignalFunctions.at(bin).at(i)->Draw("same");
    }

    std::stringstream ss5;
    ss5 << "FitBereich= " << std::setprecision(1) << std::fixed << fit_x_low << " bis "
        << fit_x_high;
    leg->AddEntry((TObject *)0, ss5.str().c_str(), "");
    leg->SetFillStyle(0);
    leg->Draw("same");
  }
  c1->Print("combinedIntensitaet.pdf", "EmbedFonts");
}

// without mass bin
void DrawFitResult::DrawResultSingle(/*std::vector<double> resPara, */ double x_first,
                                     double x_last, double fit_x_low, double fit_x_high) {
  TCanvas *c1 = new TCanvas("Intensitaeten", "", 4000, 4000);
  std::vector<TF1 *> SignalFunctions;

  int NSigPara = 2 + _SigTyp; // for BW: 3, Voigt: 4

  TF1 *FunctionModel = new TF1("Total", Total, x_first, x_last, 5, 1);
  FunctionModel->SetNpx(1000);
  FunctionModel->SetLineColor(kRed);

  for (int i = 0; i < _NSignal; i++) {
    TF1 *FunctionSignal = new TF1("Signal", Signal, x_first, x_last, NSigPara, 1);
    // for (int k = 0; k < NSigPara - 1 - 1; k++) {
    //   // FunctionModel->SetParameter(k, FunctionTotal->GetParameter(k + NSigPara * i));
    //   // FunctionSignal->SetParameter(k, FunctionTotal->GetParameter(k + NSigPara * i));
    // }
    FunctionSignal->SetParameter(0, 0.0096);
    FunctionSignal->SetParameter(1, 1.288);
    FunctionSignal->SetParameter(2, 0.0183);
    FunctionSignal->SetParameter(3, 400); // magnitude

    FunctionSignal->SetNpx(1000);
    FunctionSignal->SetLineColor(kAzure + 2 + 500 * i);
    FunctionSignal->SetLineStyle(2);
    SignalFunctions.push_back(FunctionSignal);
  }
  // Model
  FunctionModel->SetParameter(0, 0.0096);
  FunctionModel->SetParameter(1, 1.288);
  FunctionModel->SetParameter(2, 0.0183);
  FunctionModel->SetParameter(3, 400); // magnitude
  FunctionModel->SetParameter(4, 0.232);
  // FunctionModel->SetParameter(5, resPara.at(5 + (bin * 6)));
  // FunctionModel->SetParameter(6, resPara.at(6 + (bin * 6)));
  // FunctionModel->SetParameter(7, resPara.at(7 + (bin * 6)));

  FunctionBackground = new TF1("Background", Background, x_first, x_last, 1, 1);
  // Background
  FunctionBackground->SetParameter(0, FunctionModel->GetParameter(4));
  // FunctionBackground->SetParameter(1, resPara.at(5 + (bin * 6)));
  // FunctionBackground->SetParameter(2, resPara.at(6 + (bin * 6)));
  // FunctionBackground->SetParameter(3, resPara.at(7 + (bin * 6)));

  FunctionBackground->SetNpx(1000);
  FunctionBackground->SetLineColor(kOrange + 7);
  FunctionBackground->SetLineStyle(2);

  // Legend
  TLegend *leg = new TLegend(0.15, 0.48, 0.45, 0.88);
  leg->SetTextFont(132);
  leg->AddEntry(GraphSingle, "Datenpunkte", "p");

  leg->AddEntry(FunctionModel, "Total", "l");

  // for (int i = 0; i < _NSignal; i++) {
  //   std::stringstream ss;
  //   ss << "FunctionSignal " << i + 1;
  //   leg->AddEntry(SignalFunctions.at(i), ss.str().c_str(), "l");

  //   if (_SigTyp == 1) { // BW
  //     std::stringstream ss1;
  //     ss1 << "m_{0} " << i + 1 << "= " << std::setprecision(4) << std::fixed
  //         << FitRes.GetResultValue().at(0 + 4 * i) << " #pm "
  //         << FitRes.GetResultValue().at(1 + 4 * i) << " GeV/#it{c}^{2}";
  //     stringstream ss2;
  //     ss2 << "#Gamma " << i + 1 << "= " << std::setprecision(4) << std::fixed
  //         << FitRes.GetResultValue().at(2 + 4 * i) << " #pm "
  //         << FitRes.GetResultValue().at(3 + 4 * i) << " GeV/#it{c}^{2}";

  //     leg->AddEntry((TObject *)0, ss1.str().c_str(), "");
  //     leg->AddEntry((TObject *)0, ss2.str().c_str(), "");
  //   } else if (_SigTyp == 2) { // Voigt
  //     std::stringstream ss1;
  //     ss1 << "m_{0} " << i + 1 << "= " << std::setprecision(4) << std::fixed
  //         << FitRes.GetResultValue().at(0 + 6 * i) << " #pm "
  //         << FitRes.GetResultValue().at(1 + 6 * i) << " GeV/#it{c}^{2}";
  // std::cout << "Draw Fit Result  " << __LINE__ << std::endl;

  //     stringstream ss2;
  //     ss2 << "#Gamma " << i + 1 << "= " << std::setprecision(4) << std::fixed
  //         << FitRes.GetResultValue().at(2 + 6 * i) << " #pm "
  //         << FitRes.GetResultValue().at(3 + 6 * i) << " GeV/#it{c}^{2}";

  //     stringstream ss3;
  //     ss3 << "#sigma " << i + 1 << "= " << std::setprecision(4) << std::fixed
  //         << FitRes.GetResultValue().at(4) << " #pm " << FitRes.GetResultValue().at(5)
  //         << " GeV/#it{c}^{2}";
  //     // ss3 << "#sigma " << i + 1 << "= " << std::setprecision(4) << std::fixed
  //     //     << FitRes.GetResultValue().at(4 + bin * 2) << " #pm "
  //     //     << FitRes.GetResultValue().at(5 + bin * 2) << " GeV/#it{c}^{2}";

  //     leg->AddEntry((TObject *)0, ss1.str().c_str(), "");
  //     leg->AddEntry((TObject *)0, ss2.str().c_str(), "");
  //     leg->AddEntry((TObject *)0, ss3.str().c_str(), "");
  //   }
  //   leg->SetFillStyle(0);
  // }

  leg->AddEntry(FunctionBackground, "FunctionBackground", "l");

  // std::stringstream ss;
  // ss << "#chi^{2}_{red.}=" << std::setprecision(2) << std::fixed << FitRes.GetRedChi();
  // leg->AddEntry((TObject *)0, ss.str().c_str(), "");
  // leg->SetFillStyle(0);

  // Draw
  // FunctionModel->GetYaxis()->SetRangeUser(0, 1.1 * FunctionModel->GetMaximum(x_first, x_last));

  FunctionModel->SetTitle(";#it{m_{#pi^{#minus}#pi^{+}#eta}} [GeV/#it{c}^{2}];Yields");
  GraphSingle->GetXaxis()->SetRangeUser(x_first, x_last);
  GraphSingle->GetYaxis()->SetRangeUser(0, 1.1 * (GraphSingle->GetYaxis()->GetXmax()));

  GraphSingle->Draw("apz");

  // FunctionModel->Draw("same");
  FunctionBackground->Draw("same");

  // for (int i = 0; i < _NSignal; i++) {
  //   SignalFunctions.at(i)->Draw("same");
  // }

  std::stringstream ss5;
  ss5 << "FitBereich= " << std::setprecision(1) << std::fixed << fit_x_low << " bis " << fit_x_high;

  leg->AddEntry((TObject *)0, ss5.str().c_str(), "");

  leg->SetFillStyle(0);
  // leg->Draw("same");
  c1->cd();

  c1->Print("TotalIntensitaet.pdf", "EmbedFonts");
}

// draw bg result for all Bin
void DrawFitResult::DrawBGResult(std::vector<double> *resPara, double x_first, double x_last,
                                 double fit_x_low, double fit_x_high) {

  TCanvas *c2 = new TCanvas("IntensitaetenBG", "", 4000, 4000);
  // c2->Divide(2, 2);

  TLegend *leg = new TLegend(0.15, 0.48, 0.45, 0.88);
  leg->SetTextFont(132);
  // leg->AddEntry(FunctionModel, "Total", "l");
  for (int bin = 0; bin < 4; bin++) {

    stringstream ss;
    // ss << "Background_" << bin;
    FunctionBackground = new TF1(ss.str().c_str(), Background, x_first, x_last, 4, 1);
    // Background
    FunctionBackground->SetParameter(0, resPara[bin].at((0)));
    FunctionBackground->SetParameter(1, resPara[bin].at(1));
    FunctionBackground->SetParameter(2, resPara[bin].at(2));
    FunctionBackground->SetParameter(3, resPara[bin].at(3));

    FunctionBackground->SetNpx(1000);
    FunctionBackground->SetLineColor(kOrange + 7);
    FunctionBackground->SetLineStyle(2);
    leg->AddEntry(FunctionBackground, "Fit", "l");
    Graph[bin]->GetXaxis()->SetRangeUser(x_first, x_last);
    Graph[bin]->GetYaxis()->SetRangeUser(0, 1.1 * (Graph[1]->GetYaxis()->GetXmax()));
    Graph[bin]->Draw("apz");
    leg->AddEntry(Graph[bin], "Datenpunkte", "p");
    // std::stringstream ss;
    ss << "#chi^{2}_{red.}=" << std::setprecision(2) << std::fixed << FitRes.GetRedChi();
    leg->AddEntry((TObject *)0, ss.str().c_str(), "");
    leg->SetFillStyle(0);
    // std::cout << "Draw" << __LINE__ << endl;
    c2->cd(bin + 1);

    FunctionBackground->Draw("same");
  }
  c2->Print("CombinedBackground.pdf", "EmbedFonts");
}
void DrawFitResult::DrawBGSGResult(std::vector<double> resPara, double x_first, double x_last,
                                   double fit_x_low, double fit_x_high) {

  TCanvas *c2 = new TCanvas("IntensitaetenBG", "", 4000, 4000);
  // c2->Divide(2, 2);

  TLegend *leg = new TLegend(0.15, 0.65, 0.35, 0.88);
  leg->SetTextFont(132);
  // leg->AddEntry(FunctionModel, "Total", "l");
  int bin = 0;
  Graph[bin]->GetXaxis()->SetRangeUser(x_first, x_last);
  Graph[bin]->GetYaxis()->SetRangeUser(0, 1.1 * (Graph[bin]->GetYaxis()->GetXmax()));
  Graph[bin]->Draw("apz");
  leg->AddEntry(Graph[bin], "Datenpunkte", "p");

  stringstream ss;
  // ss << "Background_" << bin;
  FunctionBackground = new TF1(ss.str().c_str(), Background, x_first, x_last, 4, 1);
  // Background
  FunctionBackground->SetParameter(0, resPara[0]);
  FunctionBackground->SetParameter(1, resPara[1]);
  FunctionBackground->SetParameter(2, resPara[2]);
  FunctionBackground->SetParameter(3, resPara[3]);

  FunctionBackground->SetNpx(1000);
  FunctionBackground->SetLineColor(kOrange + 7);
  FunctionBackground->SetLineStyle(2);
  leg->AddEntry(FunctionBackground, "Fit", "l");

  // std::stringstream ss;
  ss << "#chi^{2}_{red.}=" << std::setprecision(2) << std::fixed << FitRes.GetRedChi();
  leg->AddEntry((TObject *)0, ss.str().c_str(), "");
  leg->SetFillStyle(0);
  leg->Draw("same");

  // std::cout << "Draw" << __LINE__ << endl;
  // c2->cd(bin + 1);
  c2->cd();

  FunctionBackground->Draw("same");
  c2->Print("Background.pdf", "EmbedFonts");
}

void DrawFitResult::DrawSGResult(std::vector<double> *resPara, double x_first, double x_last,
                                 double fit_x_low, double fit_x_high) {

  TCanvas *c2 = new TCanvas("IntensitaetenSG", "", 4000, 4000);
  c2->Divide(2, 2);

  for (int bin = 0; bin < 4; bin++) {

    // Create TF1 to draw the total function
    stringstream ss0;
    ss0 << "Total_" << bin;
    TF1 *FunctionTotalSG =
        new TF1(ss0.str().c_str(), Total, x_first, x_last, resPara[bin].size(), 1);

    for (int i = 0; i < 8; i++) {
      FunctionTotalSG->SetParameter(i, resPara[bin].at(i));
    }
    FunctionTotalSG->SetNpx(1000);
    FunctionTotalSG->SetLineColor(kRed);

    stringstream ss10;
    ss10 << "Signal_" << bin;
    TF1 *FunctionSignalSG = new TF1(ss10.str().c_str(), Signal, x_first, x_last, 4, 1);
    FunctionSignalSG->SetParameter(0, resPara[bin].at((0)));
    FunctionSignalSG->SetParameter(1, resPara[bin].at((1)));
    FunctionSignalSG->SetParameter(2, resPara[bin].at((2)));
    FunctionSignalSG->SetParameter(3, resPara[bin].at((3)));
    FunctionSignalSG->SetNpx(1000);
    FunctionSignalSG->SetLineColor(kAzure + 2);
    FunctionSignalSG->SetLineStyle(2);
    stringstream ss11;
    ss11 << "Background_" << bin;
    FunctionBackground = new TF1(ss11.str().c_str(), Background, x_first, x_last, 4, 1);
    // Background
    FunctionBackground->SetParameter(0, resPara[bin].at(4));
    FunctionBackground->SetParameter(1, resPara[bin].at(5));
    FunctionBackground->SetParameter(2, resPara[bin].at(6));
    FunctionBackground->SetParameter(3, resPara[bin].at(7));

    FunctionBackground->SetNpx(1000);
    FunctionBackground->SetLineColor(kOrange + 7);
    FunctionBackground->SetLineStyle(2);

    TLegend *leg = new TLegend(0.15, 0.48, 0.45, 0.88);
    leg->SetTextFont(132);
    leg->AddEntry(Graph[bin], "Datenpunkte", "p");
    leg->AddEntry(FunctionTotalSG, "Total", "l");
    std::stringstream ss1;
    ss1 << "m_{0} "
        << "= " << std::setprecision(4) << std::fixed << FitRes.GetResultValueSG()[bin].at(0)
        << " #pm " << FitRes.GetResultValueSG()[bin].at(1) << " GeV/#it{c}^{2}";
    cout << "draw" << __LINE__ << endl;

    stringstream ss2;
    ss2 << "#Gamma "
        << "= " << std::setprecision(4) << std::fixed << FitRes.GetResultValueSG()[bin].at(2)
        << " #pm " << FitRes.GetResultValueSG()[bin].at(3) << " GeV/#it{c}^{2}";
    cout << "draw" << __LINE__ << endl;

    stringstream ss3;
    ss3 << "#sigma "
        << "= " << std::setprecision(4) << std::fixed << FitRes.GetResultValueSG()[bin].at(4)
        << " #pm " << FitRes.GetResultValueSG()[bin].at(5) << " GeV/#it{c}^{2}";

    leg->AddEntry((TObject *)0, ss1.str().c_str(), "");
    leg->AddEntry((TObject *)0, ss2.str().c_str(), "");
    leg->AddEntry((TObject *)0, ss3.str().c_str(), "");

    leg->SetFillStyle(0);
    std::stringstream ss;
    ss << "#chi^{2}_{red.}=" << std::setprecision(2) << std::fixed << FitRes.GetRedChiSG()[bin];
    leg->AddEntry((TObject *)0, ss.str().c_str(), "");
    cout << "draw" << __LINE__ << endl;

    leg->SetFillStyle(0);
    std::cout << "Draw" << __LINE__ << endl;
    c2->cd(bin + 1);

    Graph[bin]->GetXaxis()->SetRangeUser(x_first, x_last);
    Graph[bin]->GetYaxis()->SetRangeUser(0, 1.1 * (Graph[bin]->GetYaxis()->GetXmax()));
    Graph[bin]->Draw("apz");
    FunctionTotalSG->Draw("same");
    FunctionSignalSG->Draw("same");
    FunctionBackground->Draw("same");
    std::stringstream ss5;
    ss5 << "FitBereich= " << std::setprecision(1) << std::fixed << fit_x_low << " bis "
        << fit_x_high;
    leg->AddEntry((TObject *)0, ss5.str().c_str(), "");
    leg->AddEntry(FunctionBackground, "FunctionBackground", "l");

    leg->SetFillStyle(0);
    leg->Draw("same");
    c2->Print("SingleIntensitaet.pdf", "EmbedFonts");
  }
}

void DrawFitResult::Plot(std::function<double(const double *, const double *)> voigt) {
  // Create TF1 to draw the total function
  TF1 *Function = new TF1("Total", voigt, 1.0, 1.4, 4, 1);
  // Parameter 0 -2: Signal, 3-6:Background
  Function->SetNpx(1000);
  Function->SetLineColor(kBlack);
  Function->SetParameter(0, 0.003);
  Function->SetParameter(1, 1.282);
  Function->SetParameter(2, 0.75);
  Function->SetParameter(3, 0.68);

  TCanvas *c2 = new TCanvas();
  c2->cd();
  Function->Draw();
  c2->Print("Voigt.pdf");
}