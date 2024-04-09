#pragma once
#include "TGraphErrors.h"
#include "TH1D.h"
#include <functional>

class Fit {
public:
  Fit(){};

  Fit(std::vector<double> *x, std::vector<double> *y, std::vector<double> *delta_y, double x_low,
      double x_high, int NSignal, int SigTyp, bool Single);
  void SetChiSquareFunction(std::function<double(const double *, const double *)> model,
                            int Nparam);
  void DoFit();
  void DoFitBG();
  void DoFitBGSG();
  void DoFitSingle();
  std::vector<double> DrawParameters(size_t NParameters);
  std::vector<double> GetResultPara() { return _resultPara; }
  std::vector<double> *GetResultSGPara() { return _resultParaSG; }

  double GetRedChi() { return redChiSqu; }
  double *GetRedChiSG() { return redChiSquSG; }

  std::vector<double> GetResultValue() { return _resultValue; }
  std::vector<double> *GetResultValueSG() { return _resultValueSG; }

  void SetBackgroundParameters(const std::vector<double> *par);
  std::vector<double> *GetBackgroundParameters() { return _BGPara; };

private:
  std::vector<std::function<double(const double *)>> chiSquare;
  std::function<double(const double *)> chiSquareSingle;

  // Vector hat nicht funktuniert

  //    = [I_ii_acc_pos, PositionOfFlatWave,
  //                                         doCauchyRegularization, gamma_inv_sq, b] {

  //     // create your function
  //     double result = 5.1;
  //     return result;
  //   }
  bool _IsSingle;
  std::vector<double> _x[4], _y[4], _delta_y[4];
  std::vector<double> _xSingle, _ySingle, _delta_ySingle;
  double _xlow, _xhigh, _minValue, _degFree, peak, width, peakErr, widthErr, sigma, sigErr;
  int _Nparam, _NSignal, _SigTyp;
  std::vector<double> _startPara, _BGPara[4], _deltaBGPara[4], _resultPara, _deltaResPara,
      _resultParaSG[4], _deltaResParaSG[4], _resultValueSG[4], _resultValue;
  double redChiSqu, redChiSquSG[4];
};