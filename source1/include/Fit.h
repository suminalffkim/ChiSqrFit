#pragma once
#include "TGraphErrors.h"
#include "TH1D.h"
#include <functional>

class Fit {
public:
  Fit(){};
  Fit(std::vector<double> *x, std::vector<double> *y, std::vector<double> *delta_y, double x_low,
      double x_high, int NSignal, int SigTyp);

  void SetChiSquareFunction(std::function<double(const double *, const double *)> model,
                            int Nparam);
  void DoFit();
  void DoFitBG();
  std::vector<double> DrawParameters(size_t NParameters);
  std::vector<double> GetResultPara() { return _resultPara; }
  double GetRedChi() { return redChiSqu; }
  std::vector<double> GetResultValue() { return _resultValue; }
  void SetBackgroundParameters(const std::vector<double> *par);
  std::vector<double> *GetBackgroundParameters() { return _BGPara; };

private:
  std::vector<std::function<double(const double *)>> chiSquare;

  //    = [I_ii_acc_pos, PositionOfFlatWave,
  //                                         doCauchyRegularization, gamma_inv_sq, b] {

  //     // create your function
  //     double result = 5.1;
  //     return result;
  //   }
  std::vector<double> _x[4], _y[4], _delta_y[4];
  double _xlow, _xhigh, _minValue, _degFree, peak, width, peakErr, widthErr, sigma, sigErr;
  int _Nparam, _NSignal, _SigTyp;
  std::vector<double> _startPara, _BGPara[4], _deltaBGPara[4], _resultPara, _deltaResPara,
      _resultValue;
  double redChiSqu;
};