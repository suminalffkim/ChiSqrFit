#pragma once
#include "Fit.h"
#include "RealData.h"
#include "RooMath.h"
#include "TAxis.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include <complex>
#include <functional>

#include <iostream>

class DrawFitResult {
public:
  DrawFitResult();
  DrawFitResult(std::function<double(const double *, const double *)> background,
                std::function<double(const double *, const double *)> signal,
                std::function<double(const double *, const double *)> total,
                std::vector<TGraphErrors *> graph, Fit &fit, int NSignal, int SigTyp);
  void DrawResult(std::vector<double> resPara, double x_first, double x_last, double fit_x_low,
                  double fit_x_high);
  void Plot(std::function<double(const double *, const double *)> voigt);

private:
  std::function<double(const double *, const double *)> Background;
  std::function<double(const double *, const double *)> Signal;
  std::function<double(const double *, const double *)> Total;

  std::vector<TGraphErrors *> Graph;
  TF1 *FunctionTotal;
  TF1 *FunctionSignal;
  TF1 *FunctionSignal2;
  std::vector<TF1 *> SignalFunctions;

  TF1 *FunctionBackground;
  Fit FitRes;
  // RealData rd;
  int _NSignal, _SigTyp;
};
