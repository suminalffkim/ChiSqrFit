#include "RealData.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH2D.h"
#include "TLegend.h"
#include "TLorentzVector.h"
#include "TTreeReader.h"
#include <iostream>
#include <sstream>
using namespace std;

void RealData::PlotData() {
  // Read File
  std::string File = "/mnt/data/compass/2009/3Pi2G_Eventselection_Slot5/EventSelection_2023/"
                     "PiPiPiEta/Total_08_09_June29.root";
  TFile f(File.c_str(), "read");

  if (!f.IsOpen()) {
    std::cout << "Could not open " << File << "!" << endl;
    abort();
  }

  TTreeReader TreeReader("AfterSelection", &f); // !!!AfterSelection PseudoData_MC
  if (TreeReader.GetEntries(true) == -1) {
    std::cout << "There is no tree called \"AfterSelection\" in file \"" << File << "\"!" << endl;
    abort();
  }

  TTreeReaderValue<std::vector<double>> BeamPion(TreeReader, "BeamPion");
  TTreeReaderValue<std::vector<double>> PiMinus_1(TreeReader, "PiMinus1"); // PiMinus1
  TTreeReaderValue<std::vector<double>> PiMinus_2(TreeReader, "PiMinus2"); // PiMinus2
  TTreeReaderValue<std::vector<double>> PiPlus_(TreeReader, "PiPlus");
  TTreeReaderValue<std::vector<double>> Eta_(TreeReader, "Neutral"); // Eta
  TTreeReaderValue<double> mx(TreeReader, "m_x");
  TTreeReaderValue<double> t_prime(TreeReader, "tprime");

  int counter = 1;
  double FitInvtervall[] = {0.8, 1.5, 2.0, 2.5, 3.0};

  for (int intv = 0; intv < 4; intv++) {
    // temp Histogram for piplpimieta + pimipipleta
    stringstream ss2;
    ss2 << "Hpipluspiminuseta" << intv;
    TH1D *Hist = new TH1D(
        ss2.str().c_str(),
        "Title; #it{m_{#pi^{#minus}#pi^{#plus}#eta}} [GeV/#it{c}^{2}]; Anzahl / (2MeV/#it{c}^{2})",
        1000, 0, 3);

    stringstream ss3;
    ss3 << "Hpiminuspiminuseta" << intv;
    TH1D *HistBG = new TH1D(
        ss3.str().c_str(),
        "Title; #it{m_{#pi^{#minus}#pi^{#minus}#eta}} [GeV/#it{c}^{2}]; Anzahl / (2MeV/#it{c}^{2})",
        1000, 0, 3);

    for (auto i = TreeReader.begin(); i != TreeReader.end(); ++i) {

      // if (counter % (TreeReader.GetEntries() / 100) == 0) {
      //   cout << "... " << int(100.0 * double(counter) / double(TreeReader.GetEntries()))
      //        << " % of the real data have been read ..." << endl;
      // }

      TreeReader.SetEntry(*i);

      std::pair<double, double> tp_mx = std::make_pair(*t_prime, *mx);
      // cout << "m_x= " << tp_mx.second << endl;
      // Hmx->Fill(tp_mx.second); // m_x
      // Htp->Fill(tp_mx.first);  // t_prime
      // H_mx_tp->Fill(tp_mx.second, tp_mx.first);

      TLorentzVector Beam(BeamPion->at(0), BeamPion->at(1), BeamPion->at(2), BeamPion->at(3));
      TLorentzVector PiMinus1(PiMinus_1->at(0), PiMinus_1->at(1), PiMinus_1->at(2),
                              PiMinus_1->at(3));
      TLorentzVector PiMinus2(PiMinus_2->at(0), PiMinus_2->at(1), PiMinus_2->at(2),
                              PiMinus_2->at(3));
      TLorentzVector PiPlus(PiPlus_->at(0), PiPlus_->at(1), PiPlus_->at(2), PiPlus_->at(3));
      TLorentzVector Eta(Eta_->at(0), Eta_->at(1), Eta_->at(2), Eta_->at(3));

      // Pi+Pi-Eta
      TLorentzVector PiplusPiminusEta = PiMinus1 + PiPlus + Eta;

      if ((tp_mx.second > FitInvtervall[intv]) && (tp_mx.second <= FitInvtervall[intv + 1])) {
        Hist->Fill(PiplusPiminusEta.M());
        TLorentzVector PiMinusPiMinusEta = PiMinus1 + PiMinus2 + Eta;

        // pi-pi-Eta
        HistBG->Fill(PiMinusPiMinusEta.M());
      }
      counter++;
    }

    Hpipluspiminuseta_m.push_back(Hist);
    Hpiminuspiminuseta_m.push_back(HistBG);
    Hpiminuspiminuseta_m.at(intv)->Rebin(4);
  }

  for (int intv = 0; intv < 4; intv++) {

    std::cout << "Draw Histogram" << endl;
    TCanvas *c1 = new TCanvas("c1", "pi+pi-eta", 1000, 1000);
    c1->cd();

    TH1D *HClone = dynamic_cast<TH1D *>(Hpipluspiminuseta_m.at(intv)->Clone());
    HClone->SetDirectory(0);
    Hpipluspiminuseta_m.at(intv)->SetDirectory(0);
    Hpiminuspiminuseta_m.at(intv)->SetDirectory(0);

    HClone->SetName("HPiPlPiMiEtaClone");
    HClone->Rebin(10);

    HClone->GetXaxis()->SetRangeUser(0.8, 3);
    HClone->Draw();

    //   c1->Modified();
    //   c1->Update();
    std::stringstream ss;
    ss << "PiplusPiminusEtaPlot_" << intv << ".pdf";
    c1->Print(ss.str().c_str(), "EmbedFonts");
  }
  f.Close();

  // TCanvas *c2 = new TCanvas("c2", "pi-pi-eta", 1000, 1000);
  // c2->cd();
  // Hpiminuspiminuseta_m->GetXaxis()->SetRangeUser(1, 2);
  // Hpiminuspiminuseta_m->Draw();
  // c2->Print("PiminusPiminusEtaPlot.pdf", "EmbedFonts");
  std::cout << "Finish" << endl;
  return;
}

void RealData::PlotDataSingle() {
  // Read File
  std::string File = "/mnt/data/compass/2009/3Pi2G_Eventselection_Slot5/EventSelection_2023/"
                     "PiPiPiEta/Total_08_09_June29.root";
  TFile f(File.c_str(), "read");

  if (!f.IsOpen()) {
    std::cout << "Could not open " << File << "!" << endl;
    abort();
  }

  TTreeReader TreeReader("AfterSelection", &f); // !!!AfterSelection PseudoData_MC
  if (TreeReader.GetEntries(true) == -1) {
    std::cout << "There is no tree called \"AfterSelection\" in file \"" << File << "\"!" << endl;
    abort();
  }

  TTreeReaderValue<std::vector<double>> BeamPion(TreeReader, "BeamPion");
  TTreeReaderValue<std::vector<double>> PiMinus_1(TreeReader, "PiMinus1"); // PiMinus1
  TTreeReaderValue<std::vector<double>> PiMinus_2(TreeReader, "PiMinus2"); // PiMinus2
  TTreeReaderValue<std::vector<double>> PiPlus_(TreeReader, "PiPlus");
  TTreeReaderValue<std::vector<double>> Eta_(TreeReader, "Neutral"); // Eta
  TTreeReaderValue<double> mx(TreeReader, "m_x");
  TTreeReaderValue<double> t_prime(TreeReader, "tprime");

  // temp Histogram for piplpimieta + pimipipleta

  TH1D *HistSingle = new TH1D(
      "HpipluspiminusetaGesamt",
      "Title; #it{m_{#pi^{#minus}#pi^{#plus}#eta}} [GeV/#it{c}^{2}]; Anzahl / (2MeV/#it{c}^{2})",
      1000, 0, 3);

  TH1D *HistBGSingle = new TH1D(
      "HpiminuspiminusetaGesamt",
      "Title; #it{m_{#pi^{#minus}#pi^{#minus}#eta}} [GeV/#it{c}^{2}]; Anzahl / (2MeV/#it{c}^{2})",
      1000, 0, 3);

  for (auto i = TreeReader.begin(); i != TreeReader.end(); ++i) {

    // if (counter % (TreeReader.GetEntries() / 100) == 0) {
    //   cout << "... " << int(100.0 * double(counter) / double(TreeReader.GetEntries()))
    //        << " % of the real data have been read ..." << endl;
    // }

    TreeReader.SetEntry(*i);

    std::pair<double, double> tp_mx = std::make_pair(*t_prime, *mx);
    // cout << "m_x= " << tp_mx.second << endl;
    // Hmx->Fill(tp_mx.second); // m_x
    // Htp->Fill(tp_mx.first);  // t_prime
    // H_mx_tp->Fill(tp_mx.second, tp_mx.first);

    TLorentzVector Beam(BeamPion->at(0), BeamPion->at(1), BeamPion->at(2), BeamPion->at(3));
    TLorentzVector PiMinus1(PiMinus_1->at(0), PiMinus_1->at(1), PiMinus_1->at(2), PiMinus_1->at(3));
    TLorentzVector PiMinus2(PiMinus_2->at(0), PiMinus_2->at(1), PiMinus_2->at(2), PiMinus_2->at(3));
    TLorentzVector PiPlus(PiPlus_->at(0), PiPlus_->at(1), PiPlus_->at(2), PiPlus_->at(3));
    TLorentzVector Eta(Eta_->at(0), Eta_->at(1), Eta_->at(2), Eta_->at(3));

    // Pi+Pi-Eta
    TLorentzVector PiplusPiminusEta = PiMinus1 + PiPlus + Eta;
    TLorentzVector PiMinusPiMinusEta = PiMinus1 + PiMinus2 + Eta;
    // for single total histogram
    HistSingle->Fill(PiplusPiminusEta.M());
    HistBGSingle->Fill(PiMinusPiMinusEta.M());
  }

  Hpipluspiminuseta_mSingle = HistSingle;
  Hpiminuspiminuseta_mSingle = HistBGSingle;

  std::cout << "Draw Histogram" << endl;
  TCanvas *c2 = new TCanvas("c2", "pi+pi-eta", 1000, 1000);
  c2->cd();
  Hpiminuspiminuseta_mSingle->Rebin(4);
  TH1D *HClone = dynamic_cast<TH1D *>(Hpipluspiminuseta_mSingle->Clone());
  HClone->SetDirectory(0);
  Hpipluspiminuseta_mSingle->SetDirectory(0);
  Hpiminuspiminuseta_mSingle->SetDirectory(0);

  HClone->SetName("HPiPlPiMiEtaClone");
  HClone->Rebin(10);

  HClone->GetXaxis()->SetRangeUser(0.8, 3);
  HClone->Draw();

  //   c1->Modified();
  //   c1->Update();
  std::stringstream ss;
  ss << "PiplusPiminusEtaPlot_total"
     << ".pdf";
  c2->Print(ss.str().c_str(), "EmbedFonts");

  f.Close();

  // TCanvas *c2 = new TCanvas("c2", "pi-pi-eta", 1000, 1000);
  // c2->cd();
  // Hpiminuspiminuseta_m->GetXaxis()->SetRangeUser(1, 2);
  // Hpiminuspiminuseta_m->Draw();
  // c2->Print("PiminusPiminusEtaPlot.pdf", "EmbedFonts");
  cout << "Finish" << endl;
  return;
}

void RealData::GetData() {
  // Hpipluspiminuseta_m->Rebin(2);

  for (int i = 0; i < 4; i++) {
    _x[i] = {};
    _y[i] = {};
    _delta_y[i] = {};

    for (int bin = 0; bin < Hpipluspiminuseta_m.at(i)->GetNbinsX(); bin++) {
      _x[i].push_back(Hpipluspiminuseta_m.at(i)->GetBinCenter(bin));
      _y[i].push_back(Hpipluspiminuseta_m.at(i)->GetBinContent(bin));
      _delta_y[i].push_back(std::sqrt(Hpipluspiminuseta_m.at(i)->GetBinContent(bin)));
      _delta_x[i].push_back(0);
    }

    // Hpiminuspiminuseta_m->Rebin(2);
    std::vector<double> delta_x2;
    for (int bin = 0; bin < Hpiminuspiminuseta_m.at(i)->GetNbinsX(); bin++) {
      _xPiMiPiMiEta[i].push_back(Hpiminuspiminuseta_m.at(i)->GetBinCenter(bin));
      _yPiMiPiMiEta[i].push_back(Hpiminuspiminuseta_m.at(i)->GetBinContent(bin));
      _delta_y_PiMiPiMiEta[i].push_back(std::sqrt(Hpiminuspiminuseta_m.at(i)->GetBinContent(bin)));
      delta_x2.push_back(0);
    }

    std::cout << "We just stored the data we need for Chi Square.." << std::endl;
    TGraphErrors *Graph = new TGraphErrors(_x[i].size(), &_x[i].front(), &_y[i].front(),
                                           &_delta_x[i].front(), &_delta_y[i].front());
    GData.push_back(Graph);
    GData.at(i)->SetTitle(
        ";#it{m_{#pi^{#plus}#pi^{#minus}#eta}} [GeV/#it{c}^{2}]; Anzahl / (2MeV/#it{c}^{2})");
    GData.at(i)->SetMarkerStyle(20);
    GData.at(i)->SetMarkerSize(1.5);
    GData.at(i)->SetLineWidth(2);

    // TCanvas *c1 = new TCanvas("GData", "Graph Data", 1000, 1000);
    // c1->cd();
    // GData->Draw("APZ");
    // c1->Print("GraphOfData.pdf", "EmbedFonts");

    // for (int i = 0; i < _y.size(); i++) {
    //   cout << "y= " << _y[i] << endl;
    // }
  }
}

void RealData::GetDataSG() {
  // Hpipluspiminuseta_m->Rebin(2);
  _xSingle = {};
  _ySingle = {};
  _delta_ySingle = {};

  for (int bin = 0; bin < Hpipluspiminuseta_mSingle->GetNbinsX(); bin++) {
    _xSingle.push_back(Hpipluspiminuseta_mSingle->GetBinCenter(bin));
    _ySingle.push_back(Hpipluspiminuseta_mSingle->GetBinContent(bin));
    _delta_ySingle.push_back(std::sqrt(Hpipluspiminuseta_mSingle->GetBinContent(bin)));
    _delta_xSingle.push_back(0);
  }

  // Hpiminuspiminuseta_m->Rebin(2);
  std::vector<double> delta_x2;

  for (int bin = 0; bin < Hpiminuspiminuseta_mSingle->GetNbinsX(); bin++) {
    _xPiMiPiMiEtaSingle.push_back(Hpiminuspiminuseta_mSingle->GetBinCenter(bin));
    _yPiMiPiMiEtaSingle.push_back(Hpiminuspiminuseta_mSingle->GetBinContent(bin));
    _delta_y_PiMiPiMiEtaSingle.push_back(std::sqrt(Hpiminuspiminuseta_mSingle->GetBinContent(bin)));
    _delta_x2Single.push_back(0);
  }

  std::cout << "We just stored the data we need for Chi Square without binning.." << std::endl;
  GDataSingle = new TGraphErrors(_xSingle.size(), &_xSingle.front(), &_ySingle.front(),
                                 &_delta_xSingle.front(), &_delta_ySingle.front());
  GDataSingle->SetTitle(
      ";#it{m_{#pi^{#plus}#pi^{#minus}#eta}} [GeV/#it{c}^{2}]; Anzahl / (2MeV/#it{c}^{2})");
  GDataSingle->SetMarkerStyle(20);
  GDataSingle->SetMarkerSize(1.5);
  GDataSingle->SetLineWidth(2);

  // TCanvas *c1 = new TCanvas("GData", "Graph Data", 1000, 1000);
  // c1->cd();
  // GData->Draw("APZ");
  // c1->Print("GraphOfData.pdf", "EmbedFonts");

  // for (int i = 0; i < _y.size(); i++) {
  //   cout << "y= " << _y[i] << endl;
  // }
}
RealData::RealData() {
  cout << "Get Real Data" << endl;
  // PlotData();
  PlotDataSingle();
  GetDataSG();
  // GetData();
}
