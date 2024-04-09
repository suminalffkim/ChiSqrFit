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
  std::vector<double> *GetXPiMinusSingle() { return &_xPiMiPiMiEtaSingle; }
  std::vector<double> *GetYPiMinusSingle() { return &_yPiMiPiMiEtaSingle; }

  double GetXBegin() { return x_begin; }
  std::vector<double> *GetYPiMinus() { return _yPiMiPiMiEta; }
  std::vector<double> *GetDeltaYPiMinus() { return _delta_y_PiMiPiMiEta; }
  std::vector<TGraphErrors *> GetGraph() { return GData; }
  TGraphErrors *GetGraphSingle() { return GDataSingle; }

  std::vector<TGraphErrors *> GetGraphBG() { return GDataBG; }
  std::vector<TGraphErrors *> GetGraphBGSG() { return GDataBGSG; }
  std::vector<double> *GetXBG() { return _xBG; }
  std::vector<double> *GetYBG() { return _yBG; }
  std::vector<double> *GetDeltaYBG() { return _delta_yBG; }

private:
  // with mass bin
  void PlotData();
  void GetData();
  // without mass binning
  void PlotDataSingle();
  void GetDataSG();
  TH1D *Hpipiminuseta_m;
  // TH1D *Hpipluspiminuseta_m;
  std::vector<TH1D *> Hpipluspiminuseta_m, Hpiminuspiminuseta_m, HistBGSG;
  TH1D *Hpipluspiminuseta_mSingle, *Hpiminuspiminuseta_mSingle;
  std::vector<TGraphErrors *> GData, GDataBG, GDataBGSG;
  TGraphErrors *GDataSingle;
  std::vector<double> _x[4], _xPiMiPiMiEta[4], _y[4], _yPiMiPiMiEta[4];
  std::vector<double> _delta_ySingle, _delta_xSingle, _delta_x2Single, _xSingle,
      _xPiMiPiMiEtaSingle, _ySingle, _yPiMiPiMiEtaSingle, _delta_y_PiMiPiMiEtaSingle;
  std::vector<double> _delta_x[4], _delta_y[4], _delta_y_PiMiPiMiEta[4], _delta_x2[4];
  double x_begin;
  int Ninterval;
  std::vector<double> _xBG[4], _yBG[4], _delta_yBG[4], _delta_xBG[4];
};
