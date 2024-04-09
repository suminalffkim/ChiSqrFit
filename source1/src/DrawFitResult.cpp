#include "DrawFitResult.h"
#include "RooMath.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TLine.h"
#include <complex>

using namespace std;

DrawFitResult::DrawFitResult(std::function<double(const double *, const double *)> background,
                             std::function<double(const double *, const double *)> signal,
                             std::function<double(const double *, const double *)> total,
                             std::vector<TGraphErrors *> graph, Fit &fit, int NSignal, int SigTyp)
    : Background(background), Signal(signal), Total(total), Graph(graph), FitRes(fit),
      _NSignal(NSignal), _SigTyp(SigTyp) {}

//
void DrawFitResult::DrawResult(std::vector<double> resPara, double x_first, double x_last,
                               double fit_x_low, double fit_x_high) {
  TCanvas *c1 = new TCanvas("Intensitaeten", "", 4000, 4000);
  c1->Divide(2, 2);

  for (int bin = 0; bin < 4; bin++) {
    // Create TF1 to draw the total function
    FunctionTotal = new TF1("Total", Total, x_first, x_last, resPara.size(), 1);
    // Parameter 0 -2: Signal, 3-6:Background

    for (int i = 0; i < int(resPara.size()); i++) {
      FunctionTotal->SetParameter(i, resPara.at(i));
    }

    FunctionTotal->SetNpx(1000);
    FunctionTotal->SetLineColor(kRed);
    std::vector<TF1 *> SignalFunctions;

    int NSigPara = 2 + _SigTyp; // for BW: 3, Voigt: 4

    for (int i = 0; i < _NSignal; i++) {
      stringstream ss;
      ss << "f_" << i;

      TF1 *FunctionSignal = new TF1(ss.str().c_str(), Signal, x_first, x_last, NSigPara, 1);
      for (int k = 0; k < NSigPara; k++) {
        FunctionSignal->SetParameter(k, FunctionTotal->GetParameter(k + NSigPara * i));
      }
      FunctionSignal->SetNpx(1000);
      FunctionSignal->SetLineColor(kAzure + 2 + 500 * i);
      FunctionSignal->SetLineStyle(2);
      SignalFunctions.push_back(FunctionSignal);
    }

    FunctionBackground = new TF1("Background", Background, x_first, x_last, 4, 1);
    // Background
    FunctionBackground->SetParameter(0, FunctionTotal->GetParameter(4 + (bin * 4) + bin));
    FunctionBackground->SetParameter(1, FunctionTotal->GetParameter(1 + 4 + (bin * 4) + bin));
    FunctionBackground->SetParameter(2, FunctionTotal->GetParameter(2 + 4 + (bin * 4) + bin));
    FunctionBackground->SetParameter(3, FunctionTotal->GetParameter(3 + 4 + (bin * 4) + bin));

    FunctionBackground->SetNpx(1000);
    FunctionBackground->SetLineColor(kOrange + 7);
    FunctionBackground->SetLineStyle(2);

    // Legend
    TLegend *leg = new TLegend(0.15, 0.48, 0.45, 0.88);
    leg->SetTextFont(132);
    leg->AddEntry(Graph[bin], "Datenpunkte", "p");
    leg->AddEntry(FunctionTotal, "Total", "l");
    for (int i = 0; i < _NSignal; i++) {
      std::stringstream ss;
      ss << "FunctionSignal " << i + 1;
      leg->AddEntry(SignalFunctions.at(i), ss.str().c_str(), "l");

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
            << FitRes.GetResultValue().at(4 + 6 * i) << " #pm "
            << FitRes.GetResultValue().at(5 + 6 * i) << " GeV/#it{c}^{2}";

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
    FunctionTotal->GetYaxis()->SetRangeUser(0, 1.1 * FunctionTotal->GetMaximum(x_first, x_last));
    FunctionTotal->SetTitle(";#it{m_{#pi^{#minus}#pi^{+}#eta}} [GeV/#it{c}^{2}];Yields");
    // std::cout << "Draw" << __LINE__ << endl;
    c1->cd(bin + 1);

    Graph[bin]->GetXaxis()->SetRangeUser(x_first, x_last);
    Graph[bin]->GetYaxis()->SetRangeUser(0, 1.1 * (Graph[bin]->GetYaxis()->GetXmax()));
    Graph[bin]->Draw("apz");
    FunctionTotal->Draw("same");
    FunctionBackground->Draw("same");
    for (int i = 0; i < _NSignal; i++) {
      SignalFunctions.at(i)->Draw("same");
    }

    // double y_high = Graph->GetYaxis()->GetXmax();
    // // Draw vertikal Line for Fit Range
    // TLine *firstLine = new TLine(fit_x_low, 0, fit_x_low, 3000);
    // TLine *secondLine = new TLine(fit_x_high, 0, fit_x_high, 1.1 * y_high);
    // firstLine->SetLineColor(kRed);
    // secondLine->SetLineColor(kRed);
    // firstLine->SetLineWidth(2);
    // secondLine->SetLineWidth(2);
    // firstLine->SetLineStyle(3);
    // secondLine->SetLineStyle(3);
    // firstLine->Draw("same");
    // secondLine->Draw("same");
    // leg->AddEntry(firstLine, "Fit Bereich", "l");
    std::stringstream ss5;
    ss5 << "FitBereich= " << std::setprecision(1) << std::fixed << fit_x_low << " bis "
        << fit_x_high;
    leg->AddEntry((TObject *)0, ss5.str().c_str(), "");
    leg->SetFillStyle(0);
    leg->Draw("same");
    c1->Print("combinedIntensitaet.pdf", "EmbedFonts");
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
  Function->SetParameter(2, 0.024);

  TCanvas *c2 = new TCanvas();
  c2->cd();
  Function->Draw();
  c2->Print("Voigt.pdf");
}