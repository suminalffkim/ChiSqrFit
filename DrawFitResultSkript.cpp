#include "TAxis.h"
#include "TF1.h"
#include "TLegend.h"
#include <complex>
#include <iostream>

void DrawFitResult() {

  // Pure signal function
  std::function<double(double *, double *)> Signal = [](double *x, double *par) -> double {
    //[0] magnitude
    //[1] position
    //[2] width
    std::complex<double> result =
        par[0] * par[1] * par[2] /
        (par[1] * par[1] - x[0] * x[0] - std::complex<double>(0, 1) * par[1] * par[2]);
    return std::pow(std::fabs(result), 2);
  };

  // Pure background function
  std::function<double(double *, double *)> Background = [](double *x, double *par) -> double {
    // return pol3 with 4 parameter
    return par[0] * std::pow(x[0], 3) + par[1] * std::pow(x[0], 2) + par[2] * x[0] + par[3];
  };

  // Combine  signal&background to total
  std::function<double(double *, double *)> Total = [Signal, Background](double *var,
                                                                         double *par) -> double {
    return Signal(&var[0], &par[0]) + Background(&var[0], &par[3]);
  };

  // Create TF1 to draw the total function
  TF1 *FunctionTotal = new TF1("Intensity & Phase", Total, 1.18, 1.4, 7, 1);
  // Signal
  FunctionTotal->SetParameter(0, 80);
  FunctionTotal->SetParameter(1, 1.281);
  FunctionTotal->SetParameter(2, 0.022);
  // Background
  FunctionTotal->SetParameter(3, -5.12936e+04);
  FunctionTotal->SetParameter(4, 1.98654e+05);
  FunctionTotal->SetParameter(5, -2.45677e+05);
  FunctionTotal->SetParameter(6, 9.89954e+04);
  FunctionTotal->SetNpx(1000);
  FunctionTotal->SetLineColor(kBlack);

  // Create TF1 to draw signal part
  TF1 *FunctionSignal = new TF1("Intensity & Phase", Signal, 1.18, 1.4, 3, 1);
  // Signal
  FunctionSignal->SetParameter(0, FunctionTotal->GetParameter(0));
  FunctionSignal->SetParameter(1, FunctionTotal->GetParameter(1));
  FunctionSignal->SetParameter(2, FunctionTotal->GetParameter(2));
  FunctionSignal->SetNpx(1000);
  FunctionSignal->SetLineColor(kAzure + 2);

  // Create TF1 to draw signal part
  TF1 *FunctionBackground = new TF1("Intensity & Phase", Background, 1.18, 1.4, 4, 1);
  // Background
  FunctionBackground->SetParameter(0, FunctionTotal->GetParameter(3));
  FunctionBackground->SetParameter(1, FunctionTotal->GetParameter(4));
  FunctionBackground->SetParameter(2, FunctionTotal->GetParameter(5));
  FunctionBackground->SetParameter(3, FunctionTotal->GetParameter(6));
  FunctionBackground->SetNpx(1000);
  FunctionBackground->SetLineColor(kOrange + 7);

  // Draw
  FunctionTotal->GetYaxis()->SetRangeUser(0, 1.1 * FunctionTotal->GetMaximum(1.18, 1.4));
  FunctionTotal->SetTitle(";#it{m_{#pi^{#minus}#pi^{+}#eta}} [GeV/#it{c}^{2}];Yields");
  FunctionTotal->Draw();
  FunctionBackground->Draw("same");
  FunctionSignal->Draw("same");

  // Legend
  TLegend *leg = new TLegend(0.6, 0.7, 0.85, 0.85);
  leg->SetTextFont(132);
  leg->AddEntry(FunctionTotal, "Total", "l");
  leg->AddEntry(FunctionSignal, "FunctionSignal", "l");
  leg->AddEntry(FunctionBackground, "FunctionBackground", "l");
  leg->AddEntry((TObject *)0, "#chi^{2}_{red.}=");
  leg->SetFillStyle(0);
  leg->Draw("same");
}