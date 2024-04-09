#include "Fit.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/Minimizer.h"
#include "Math/MinimizerOptions.h"
#include "TMath.h"
#include <algorithm>
#include <iostream>
#include <iterator>
#include <random>
using namespace std;

Fit::Fit(std::vector<double> *x, std::vector<double> *y, std::vector<double> *delta_y, double x_low,
         double x_high, int NSignal, int SigTyp,
         bool IsSingle) //_x(x), _y(y), _delta_y(delta_y), _xlow(x_low), _xhigh(x_high)
{
  _IsSingle = IsSingle;

  _xlow = x_low;
  _xhigh = x_high;
  _NSignal = NSignal;
  _SigTyp = SigTyp;
  if (!_IsSingle) {

    std::vector<double> x_tem[4], y_tem[4], deltay_temp[4];
    for (int i = 0; i < 4; i++) {
      _x[i] = x[i];
      // cout << _x[i].size() << endl;
      _y[i] = y[i];
      _delta_y[i] = delta_y[i];
      double x_max = 0;
      for (size_t j = 0; j < _x[i].size(); j++) {
        if (x_low <= _x[i].at(j) && x_high >= x[i].at(j)) {
          x_tem[i].push_back(_x[i].at(j));
          y_tem[i].push_back(_y[i].at(j));
          deltay_temp[i].push_back(_delta_y[i].at(j));
          // cout << "bin= " << i << endl;
          // cout << "x= " << _x[i].at(j) << endl;

          if (_x[i].at(j) > x_max) {
            x_max = _x[i].at(j);
          }
        }
      }
      _x[i].swap(x_tem[i]);
      _y[i].swap(y_tem[i]);
      _delta_y[i].swap(deltay_temp[i]);
      // cout << "bin= " << i << ", x_max=" << x_max << endl;
    }
  } else {
    std::vector<double> x_tem, y_tem, deltay_temp;
    _xSingle = *x;
    _ySingle = *y;
    _delta_ySingle = *delta_y;

    double x_max = 0;
    for (size_t j = 0; j < _xSingle.size(); j++) {
      if (x_low <= _xSingle.at(j) && x_high >= _xSingle.at(j)) {
        x_tem.push_back(_xSingle.at(j));
        y_tem.push_back(_ySingle.at(j));
        deltay_temp.push_back(_delta_ySingle.at(j));

        if (_xSingle.at(j) > x_max) {
          x_max = _xSingle.at(j);
        }
      }
    }
    _xSingle.swap(x_tem);
    _ySingle.swap(y_tem);
    _delta_ySingle.swap(deltay_temp);
  }
}

void Fit::SetChiSquareFunction(std::function<double(const double *, const double *)> model,
                               int Nparam) {
  _Nparam = Nparam;
  chiSquare = {};
  // chiSquare.resize(4);
  std::vector<double> xCopy, yCopy, deltayCopy;
  if (!_IsSingle) {

    for (int j = 0; j < 1; j++) {
      xCopy = _x[j];
      yCopy = _y[j];
      deltayCopy = _delta_y[j];

      std::function<double(const double *)> OneChiSquare =
          [xCopy = xCopy, yCopy = yCopy, deltayCopy = deltayCopy,
           model = model](const double *par) -> double {
        double result = 0;
        //   create your function
        for (size_t i = 0; i < xCopy.size(); i++) {
          if (yCopy[i] != 0 && deltayCopy[i] != 0)
            result += pow(((yCopy[i] - model(&xCopy[i], &par[0])) / deltayCopy[i]), 2);
        }
        return result;
      };
      chiSquare.push_back(OneChiSquare);
    }
  } else {
    xCopy = _xSingle;
    yCopy = _ySingle;
    deltayCopy = _delta_ySingle;
    chiSquareSingle = [xCopy = xCopy, yCopy = yCopy, deltayCopy = deltayCopy,
                       model = model](const double *par) -> double {
      double result = 0;
      //   create your function
      for (size_t i = 0; i < xCopy.size(); i++) {
        if (yCopy[i] != 0 && deltayCopy[i] != 0)
          result += pow(((yCopy[i] - model(&xCopy[i], &par[0])) / deltayCopy[i]), 2);
      }
      return result;
    };
  }
}

void Fit::DoFit() {

  std::function<double(const double *)> TotalChiSquare =
      [&, chiSquare = chiSquare](const double *par) -> double {
    double result = 0;
    for (size_t bin = 0; bin < 4; bin++) { // chiSquare.size()
      std::vector<double> par_per_bin = {};
      if (bin == 0)
        par_per_bin.push_back(par[0]);
      else
        par_per_bin.push_back(par[8 * bin - 2 * (bin - 1)]);
      for (size_t same_par = 1; same_par < 3; same_par++)
        par_per_bin.push_back(par[same_par]);
      for (size_t individual_par = 0; individual_par < 5; individual_par++)
        par_per_bin.push_back(par[individual_par + 3 + 6 * bin]);
      result += chiSquare[bin](&par_per_bin[0]);
    }
    return result;
  };

  ROOT::Math::Functor f = ROOT::Math::Functor(TotalChiSquare, _Nparam);
  ROOT::Math::Minimizer *min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
  min->SetMaxFunctionCalls(1000000); // default from root example 1000000
  min->SetMaxIterations(10000);      // default from root example 10000
  min->SetTolerance(0.01);
  // Tolerance of how small the slope should be at minimum. Default is
  // 0.01 see https://root-forum.cern.ch/t/minuit2-edm-question/18470
  min->SetPrintLevel(0);
  min->SetFunction(f);
  min->SetStrategy(2);
  //   min->SetErrorDef(0.5);

  // Generate the random Start parameters
  std::vector<double> startParas = DrawParameters(_Nparam);

  // Set backgound
  for (size_t bin = 0; bin < 4; bin++)
    for (size_t par = 0; par < 4; par++) {
      startParas[par + 4 + (bin * 6)] = _BGPara[bin].at(par);
    }
  // Set signal
  startParas.at(0) = 0.005;
  startParas.at(8) = 0.005;
  startParas.at(14) = 0.005;
  startParas.at(20) = 0.005;

  startParas.at(1) = 1.282;
  startParas.at(2) = 0.022;

  startParas.at(3) = 1000;
  startParas.at(9) = 1000;
  startParas.at(15) = 1000;
  startParas.at(21) = 1000;

  for (size_t i = 0; i < startParas.size(); i++) {
    std::stringstream Name;
    Name << "Var " << i;
    min->SetVariable(i, Name.str(), startParas.at(i), 1e-5);
    cout << "Start Para " << i << " = " << startParas.at(i) << endl;
    min->SetVariableStepSize(i, 0.001);
  }

  // 1. Signal

  if (_SigTyp == 1) { // BW
    min->SetVariableValue(0, 1);
  } else if (_SigTyp == 2) { // Voigt
    // min->SetVariableValue(0, 0.005);
    // min->SetVariableLimits(0, 0.001, 0.01);
    // min->SetVariableLowerLimit(3, 0.01);

  } else {
    abort();
  }
  min->SetVariableLimits(0, 0.0, 0.03);
  min->SetVariableLimits(8, 0.0, 0.03);
  min->SetVariableLimits(14, 0.0, 0.03);
  min->SetVariableLimits(20, 0.0, 0.03);

  // min->SetVariableValue(1, 1.282); // mu
  // min->SetVariableLimits(1, 1.27, 1.29);
  // min->SetVariableValue(2, 0.024); // gamma
  min->SetVariableLimits(2, 0.005, 0.1);

  // 2.Signal
  if (_NSignal == 2) {
    if (_SigTyp == 1) { // BW
      min->SetVariableValue(4, 1.294);
      // min->FixVariable(4);
      min->SetVariableValue(5, 0.055);

    } else if (_SigTyp == 2) {         // Voigt
      min->SetVariableValue(4, 0.003); // sigma
      min->SetVariableValue(5, 1.294); // Peak Pos
      // min->FixVariable(5);
      // min->SetVariableLimits(5, 1.270, 1.310);
      min->SetVariableValue(6, 0.055); // Peak Breite
      // min->FixVariable(6);
      min->SetVariableLowerLimit(7, 0.01);

    } else {
      abort();
    }
  }
  bool success = false;
  cout << "Start Minimization" << endl;

  while (!success) {
    success = min->Minimize();
    for (size_t p = 0; p < startParas.size(); p++)
      cout << "p=" << p << ", X(p)=" << min->X()[p] << endl;
    cout << "chi=" << min->MinValue() << endl;
    // cout << "fit" << __LINE__ << endl;
  }
  cout << "End Minimization" << endl;

  const double *x = min->X();
  const double *errors = min->Errors();

  for (int i = 0; i < _Nparam; i++) {
    _resultPara.push_back(x[i]);
    _deltaResPara.push_back(errors[i]);
    cout << "Para " << i << " = " << x[i] << " +- " << errors[i] << endl;
  }
  _minValue = min->MinValue();
  _degFree = _x[0].size() + _x[1].size() + _x[2].size() + _x[3].size() - _Nparam;
  redChiSqu = _minValue / _degFree;
  cout << "red. Chi^2 = " << redChiSqu << endl;

  // std::cout << "Tolerance: " << min->Tolerance() << std::endl;
  // std::cout << "ErrorDef: " << min->ErrorDef() << std::endl;
  // std::cout << "MaxFunctionCalls: " << min->MaxFunctionCalls() << std::endl;

  if (_SigTyp == 1) {
    for (int i = 0; i < _NSignal; i++) {
      sigma = _resultPara.at(0 + i * 4);
      sigErr = _deltaResPara.at(0 + i * 4);
      peak = _resultPara.at(1 + i * 3);
      peakErr = _deltaResPara.at(1 + i * 3);
      width = _resultPara.at(2 + i * 3);
      widthErr = _deltaResPara.at(2 + i * 3);
      _resultValue.push_back(peak);
      _resultValue.push_back(peakErr);
      _resultValue.push_back(width);
      _resultValue.push_back(widthErr);
      _resultValue.push_back(sigma);
      _resultValue.push_back(sigErr);
    }
  } else if (_SigTyp == 2) {
    for (int i = 0; i < _NSignal; i++) {
      sigma = _resultPara.at(0 + i * 4);
      sigErr = _deltaResPara.at(0 + i * 4);
      peak = _resultPara.at(1 + i * 4);
      peakErr = _deltaResPara.at(1 + i * 4);
      width = _resultPara.at(2 + i * 4);
      widthErr = _deltaResPara.at(2 + i * 4);
      _resultValue.push_back(peak);
      _resultValue.push_back(peakErr);
      _resultValue.push_back(width);
      _resultValue.push_back(widthErr);
      _resultValue.push_back(sigma);
      _resultValue.push_back(sigErr);
      for (int bin = 1; bin < 4; bin++) {
        sigma = _resultPara.at(8 * bin - 2 * (bin - 1));
        sigErr = _deltaResPara.at(8 * bin - 2 * (bin - 1));
        _resultValue.push_back(sigma);
        _resultValue.push_back(sigErr);
      }
    }
  }
}

void Fit::DoFitBG() {
  for (int j = 0; j < 4; j++) {
    ROOT::Math::Functor f = ROOT::Math::Functor(chiSquare[j], _Nparam);
    ROOT::Math::Minimizer *min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
    min->SetMaxFunctionCalls(1000000); // default from root example 1000000
    min->SetMaxIterations(10000);      // default from root example 10000
    min->SetTolerance(0.01);
    // Tolerance of how small the slope should be at minimum. Default is
    // 0.01 see https://root-forum.cern.ch/t/minuit2-edm-question/18470
    min->SetPrintLevel(0);
    min->SetFunction(f);
    //   min->SetStrategy(0);
    //   min->SetErrorDef(0.5);

    // Generate the random Start parameters
    std::vector<double> startParas = DrawParameters(_Nparam);
    for (size_t i = 0; i < startParas.size(); ++i) {
      std::stringstream Name;
      Name << "Var " << i;
      min->SetVariable(i, Name.str(), startParas.at(i), 1e-5);
      cout << "Start Para BG" << i << " = " << startParas.at(i) << endl;
    }
    bool success = false;
    cout << "Start Minimization BG" << endl;

    cout << j << endl;

    while (!success) {
      success = min->Minimize();
    }

    cout << "End Minimization BG" << endl;

    const double *x = min->X();
    const double *errors = min->Errors();

    _BGPara[j] = {};
    for (int i = 0; i < _Nparam; i++) {
      _BGPara[j].push_back(x[i]);
      _deltaBGPara[j].push_back(errors[i]);
      cout << "Para Background" << i << " = " << x[i] << " +- " << errors[i] << endl;
    }
    _minValue = min->MinValue();
    _degFree = _x[j].size() - _Nparam;
    double redChiSquBG = _minValue / _degFree;
    cout << "red. Chi^2 = " << redChiSquBG << endl;
  }
}

void Fit::DoFitSingle() {
  // for (int bin = 0; bin < 4; bin++) {
  int bin = 0;
  ROOT::Math::Functor f = ROOT::Math::Functor(chiSquare[bin], 5);
  ROOT::Math::Minimizer *min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
  min->SetMaxFunctionCalls(1000000); // default from root example 1000000
  min->SetMaxIterations(10000);      // default from root example 10000
  min->SetTolerance(0.01);
  // Tolerance of how small the slope should be at minimum. Default is
  // 0.01 see https://root-forum.cern.ch/t/minuit2-edm-question/18470
  min->SetPrintLevel(0);
  min->SetFunction(f);
  min->SetStrategy(2);
  //   min->SetErrorDef(0.5);

  // Generate the random Start parameters
  std::vector<double> startParas = DrawParameters(5);

  // // Set backgound
  // for (size_t par = 0; par < 4; par++) {
  //   startParas.at(par + 4) = _BGPara[bin].at(par);
  // }
  // Set signal

  startParas.at(0) = 0.005;
  startParas.at(1) = 1.282;
  startParas.at(2) = 0.022;
  startParas.at(3) = 1000;

  for (size_t i = 0; i < startParas.size(); i++) {
    std::stringstream Name;
    Name << "Var " << i;
    min->SetVariable(i, Name.str(), startParas.at(i), 1e-5);
    cout << "Start Para " << i << " = " << startParas.at(i) << endl;
    min->SetVariableStepSize(i, 0.001);
  }

  min->SetVariableLimits(0, 0.0, 0.03);

  // min->SetVariableValue(1, 1.282); // mu
  min->SetVariableLimits(1, 1.27, 1.29);
  // min->SetVariableValue(2, 0.024); // gamma
  min->SetVariableLimits(2, 0.005, 0.1);
  bool success = false;
  cout << "Start Minimization" << endl;

  while (!success) {
    success = min->Minimize();
    // for (size_t p = 0; p < startParas.size(); p++)
    //   cout << "p=" << p << ", X(p)=" << min->X()[p] << endl;
    cout << "chi=" << min->MinValue() << endl;
    // cout << "fit" << __LINE__ << endl;
  }
  cout << "End Minimization" << endl;

  const double *x = min->X();
  const double *errors = min->Errors();

  for (int i = 0; i < 5; i++) {
    _resultParaSG[bin].push_back(x[i]);
    _deltaResParaSG[bin].push_back(errors[i]);
    cout << "Para " << i << " = " << x[i] << " +- " << errors[i] << endl;
  }
  _minValue = min->MinValue();
  _degFree = _x[bin].size() - 5;
  redChiSquSG[bin] = _minValue / _degFree;
  cout << "red. Chi^2 = " << redChiSquSG[bin] << endl;

  // std::cout << "Tolerance: " << min->Tolerance() << std::endl;
  // std::cout << "ErrorDef: " << min->ErrorDef() << std::endl;
  // std::cout << "MaxFunctionCalls: " << min->MaxFunctionCalls() << std::endl;

  if (_SigTyp == 2) {
    for (int i = 0; i < _NSignal; i++) {
      sigma = _resultParaSG[bin].at(0 + i * 4);
      sigErr = _deltaResParaSG[bin].at(0 + i * 4);
      peak = _resultParaSG[bin].at(1 + i * 4);
      peakErr = _deltaResParaSG[bin].at(1 + i * 4);
      width = _resultParaSG[bin].at(2 + i * 4);
      widthErr = _deltaResParaSG[bin].at(2 + i * 4);
      _resultValueSG[bin].push_back(peak);
      _resultValueSG[bin].push_back(peakErr);
      _resultValueSG[bin].push_back(width);
      _resultValueSG[bin].push_back(widthErr);
      _resultValueSG[bin].push_back(sigma);
      _resultValueSG[bin].push_back(sigErr);
    }
  }
}
// }
void Fit::DoFitBGSG() {

  ROOT::Math::Functor f = ROOT::Math::Functor(chiSquare[0], 4);
  ROOT::Math::Minimizer *min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
  min->SetMaxFunctionCalls(1000000); // default from root example 1000000
  min->SetMaxIterations(10000);      // default from root example 10000
  min->SetTolerance(0.01);
  // Tolerance of how small the slope should be at minimum. Default is
  // 0.01 see https://root-forum.cern.ch/t/minuit2-edm-question/18470
  min->SetPrintLevel(0);
  min->SetFunction(f);
  //   min->SetStrategy(0);
  //   min->SetErrorDef(0.5);

  // Generate the random Start parameters
  std::vector<double> startParas = DrawParameters(4);
  for (size_t i = 0; i < startParas.size(); ++i) {
    std::stringstream Name;
    Name << "Var " << i;
    min->SetVariable(i, Name.str(), startParas.at(i), 1e-5);
    cout << "Start Para BG" << i << " = " << startParas.at(i) << endl;
  }
  bool success = false;
  cout << "Start Minimization BG" << endl;

  while (!success) {
    success = min->Minimize();
  }

  cout << "End Minimization BG" << endl;

  const double *x = min->X();
  const double *errors = min->Errors();
  _resultPara = {};
  for (int i = 0; i < 4; i++) {
    _resultPara.push_back(x[i]);
    _deltaResPara.push_back(errors[i]);
    cout << "Para " << i << " = " << x[i] << " +- " << errors[i] << endl;
  }
  _minValue = min->MinValue();
  _degFree = _x[0].size() + _x[1].size() + _x[2].size() + _x[3].size() - 4;
  redChiSqu = _minValue / _degFree;
  cout << "red. Chi^2 = " << redChiSqu << endl;
}

void Fit::SetBackgroundParameters(const std::vector<double> *par) {
  for (int i = 0; i < 4; i++) {
    _BGPara[i] = {};
    for (int j = 0; j < 4; j++) {
      _BGPara[i].push_back(par[i].at(j));
    }
  }
}

std::vector<double> Fit::DrawParameters(size_t NParameters) {
  // create random engine
  std::random_device rd;
  std::mt19937 engine{rd()};
  // vector containing the start parameters
  std::vector<double> out;
  std::uniform_real_distribution<double> dist{0, 1};
  for (size_t i = 0; i < NParameters; ++i) {
    out.push_back(dist(engine));
  }
  _startPara = out;
  return out;
}