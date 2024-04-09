#include "DrawFitResult.h"
#include "Fit.h"
#include "RealData.h"
#include "RooMath.h"
#include "TInterpreter.h"
#include <complex>
#include <iostream>
using namespace std;
int main(/*int argc, char **argv*/) {
  gInterpreter->ExecuteMacro("~/.rootlogon.C");

  RealData rd;
  double fit_x_low = 1.0; // 1.0;
  double fit_x_high = 3;  // 1.6;

  int NSignal = 1;
  // Signal Typ- 1: BW, 2: Voigt
  int SignalTyp = 2;

  // Pure signal function BW
  std::function<double(const double *, const double *)> Signal = [](const double *x,
                                                                    const double *par) -> double {
    //[0] magnitude
    //[1] position
    //[2] width
    std::complex<double> result =
        par[0] * par[1] * par[2] /
        (par[1] * par[1] - x[0] * x[0] - std::complex<double>(0, 1) * par[1] * par[2]);
    return pow(fabs(result), 2);
  };

  // return complex fucntion BW
  std::function<std::complex<double>(const double *, const double *)> Signal2 =
      [&](const double *x, const double *par) -> std::complex<double> {
    //[0] magnitude
    //[1] position
    //   //[2] width
    std::complex<double> result =
        par[0] * par[1] * par[2] /
        (par[1] * par[1] - x[0] * x[0] - std::complex<double>(0, 1) * par[1] * par[2]);
    return result;
  };

  std::function<double(const double *, const double *)> Voigt = [](const double *x,
                                                                   const double *par) -> double {
    // Pos=par[1], Width=par[2] par[3]: magnitude
    double sigma = par[0]; // std::abs(par[0]);
    double mean = par[1];
    double gamma = par[2]; // std::abs(par[2]);
    return par[3] * TMath::Voigt(x[0] - mean, sigma, gamma, 4);
  };

  // double Voigt(double x, double mean, double gamma, double sigma) {
  //   std::complex<double> z;
  //   std::complex<double> i(0, 1);
  //   z = ((x - mean) + i * gamma) / (sigma * std::sqrt(2.0));

  //   std::complex<double> faddeeva = RooMath::faddeeva(z);
  //   // std::complex<double> w = std::exp(-z * z) * RooMath::erfc(-i * z);
  //   return faddeeva.real() / (sigma * std::sqrt(2.0 * TMath::Pi()));
  // }
  // Pure background function
  std::function<double(const double *, const double *)> Background =
      [](const double *x, const double *par) -> double {
    // double binWidth = rd.GetYPiMinus().at(1).first - rd.GetYPiMinus().at(0).first;
    // int bin = int(x[0] / binWidth);
    // return par[0] * rd.GetYPiMinus().at(bin).second;
    return par[0] * std::pow(x[0], 3) + par[1] * std::pow(x[0], 2) + par[2] * x[0] + par[3];
  };

  std::function<double(const double *, const double *)> BGHist = [&](const double *x,
                                                                     const double *par) -> double {
    double binWidth = rd.GetXPiMinusSingle()->at(1) - rd.GetXPiMinusSingle()->at(0);
    int bin = int(x[0] / binWidth);
    return par[0] * rd.GetYPiMinusSingle()->at(bin);
  };

  // Combine  BW signal&background to total
  std::function<double(const double *, const double *)> model = [&](const double *var,
                                                                    const double *par) -> double {
    if (NSignal == 1) {
      return Signal(&var[0], &par[0]) + Background(&var[0], &par[3]);
    } else if (NSignal == 2) {
      return Signal(&var[0], &par[0]) + Signal(&var[0], &par[3]) + Background(&var[0], &par[6]);
    } else {
      abort();
    }
  };

  // combine Voigt + Background
  std::function<double(const double *, const double *)> model2 =
      [&, Voigt, Background](const double *var, const double *par) -> double {
    if (NSignal == 1) {
      // std::cout << "Voigt=" << Voigt(&var[0], &par[0]) << ", BG=" << Background(&var[0], &par[4])
      //           << std::endl;
      return Voigt(&var[0], &par[0]) + Background(&var[0], &par[4]);
    } else if (NSignal == 2) {
      return Voigt(&var[0], &par[0]) + Voigt(&var[0], &par[4]) + Background(&var[0], &par[8]);
    } else {
      abort();
    }
  };

  std::function<double(const double *, const double *)> modelWithHistogram =
      [&](const double *var, const double *par) -> double {
    return Voigt(&var[0], &par[0]) + BGHist(&var[0], &par[4]);
  };

  // BG
  // Fit fitBGSG(rd.GetXBG(), rd.GetYBG(), rd.GetDeltaYBG(), fit_x_low, fit_x_high, NSignal,
  //             SignalTyp);
  // fitBGSG.SetChiSquareFunction(Background, 4);
  // fitBGSG.DoFitBGSG();
  // DrawFitResult draw(Background, Voigt, model2, rd.GetGraphBGSG(), fitBGSG, NSignal, SignalTyp);
  // draw.DrawBGSGResult(fitBGSG.GetResultPara(), 1, 1.5, fit_x_low, fit_x_high);

  // Fit fitBG(rd.GetXPiMinus(), rd.GetYPiMinus(), rd.GetDeltaYPiMinus(), fit_x_low, fit_x_high,
  //           NSignal, SignalTyp);
  // fitBG.SetChiSquareFunction(Background, 4);
  // fitBG.DoFitBG();
  // std::vector<double> *BGPara = fitBG.GetBackgroundParameters();

  // // DrawFitResult drawBG(Background, Voigt, model2, rd.GetGraphBG(), fitBG, NSignal, SignalTyp);
  // // drawBG.DrawBGResult(BGPara, 1, 1.5, fit_x_low, fit_x_high);

  // Fit fit(rd.GetX(), rd.GetY(), rd.GetDeltaY(), fit_x_low, fit_x_high, NSignal, SignalTyp);
  // fit.SetBackgroundParameters(BGPara);
  // int binsize = 4;
  // fit.SetChiSquareFunction(model2, 26); // (1 + SignalTyp) * NSignal + (4 + NSignal) * binsize
  // // fit.DoFitSingle();

  // Fit fit(rd.GetX(), rd.GetY(), rd.GetDeltaY(), fit_x_low, fit_x_high, NSignal, SignalTyp, true);
  // fit.SetChiSquareFunction(modelWithHistogram, 5);
  // fit.DoFitSingle();
  DrawFitResult draw(BGHist, Voigt, modelWithHistogram, rd.GetGraphSingle(), NSignal, SignalTyp);

  double x_first = 1.0; // rd.GetX().front();
  double x_last = 3.0;  // rd.GetX().back();

  draw.DrawResultSingle(x_first, x_last, fit_x_low, fit_x_high);

  // draw.DrawSGResult(fit.GetResultSGPara(), x_first, x_last, fit_x_low, fit_x_high);

  // draw.Plot(Voigt);
  // std::cout << "main" << __LINE__ << endl;
}
