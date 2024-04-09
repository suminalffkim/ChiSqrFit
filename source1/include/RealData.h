#pragma once
#include "TGraphErrors.h"
#include "TH1D.h"

class RealData {
public:
  RealData();
  std::vector<double> *GetX() { return _x; }
  std::vector<double> *GetY() { return _y; }
  std::vector<double> *GetDeltaY() { return _delta_y; }
  std::vector<double> *GetXPiMinus() { return _xPiMiPiMiEta; }
  double GetXBegin() { return x_begin; }
  std::vector<double> *GetYPiMinus() { return _yPiMiPiMiEta; }
  std::vector<double> *GetDeltaYPiMinus() { return _delta_y_PiMiPiMiEta; }
  std::vector<TGraphErrors *> GetGraph() { return GData; }
  std::vector<TGraphErrors *> GetBGHGraph() { return GDataBGSG; };

private:
  void PlotData();
  void GetData();
  TH1D *Hpipiminuseta_m;
  // TH1D *Hpipluspiminuseta_m;
  std::vector<TH1D *> Hpipluspiminuseta_m, Hpiminuspiminuseta_m;
  std::vector<TGraphErrors *> GData, GDataBGSG;
  std::vector<double> _x[4], _xPiMiPiMiEta[4], _y[4], _yPiMiPiMiEta[4];
  std::vector<double> _delta_x[4], _delta_y[4], _delta_y_PiMiPiMiEta[4];
  double x_begin;
};
