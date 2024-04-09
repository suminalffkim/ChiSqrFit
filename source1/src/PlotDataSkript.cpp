#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLegend.h"
#include "TLorentzVector.h"
#include "TTreeReader.h"
#include <iostream>
using namespace std;

void PlotData() {
  cout << "Test" << endl;
  TH1D *Htp = new TH1D("Htp", "Title; t' [(GeV/c)^{2}]; Anzahl", 100, 0, 5);
  TH1D *Hmx = new TH1D("Hmx", "Title; m_{x} [GeV/c^{2}]; Anzahl", 100, 0, 10);
  TH2D *H_mx_tp =
      new TH2D("H_mx_tp", "2D; m_{x} [GeV/c^{2}]; t' [(GeV/c)^{2}]", 50, 1, 3, 50, 0, 1);
  TH1D *Hpipipi_mx = new TH1D("Hpipipi_mx", "Title; m_{x} [GeV/c^{2}]; Anzahl", 300, 0, 3);
  TH1D *Hpipipi_p = new TH1D("Hpipipi", "Title; p [GeV/c]; Anzahl", 100, 0, 200);
  TH1D *Heta_p = new TH1D("HEtap", "Titel; p [GeV/c]; Anzahl", 100, 0, 200);
  TH1D *Hpipi_m = new TH1D("Hpipi_m", "Title; m_{x} [GeV/c^{2}]; Anzahl", 400, 0, 2);
  TH1D *Hpiminuseta_m = new TH1D("Hpiminuseta_m", "Title; m_{x} [GeV/c^{2}]; Anzahl", 400, 0, 2);
  TH1D *Hpipluseta_m = new TH1D("Hpipluseta_m", "Title; m_{x} [GeV/c^{2}]; Anzahl", 400, 0, 2);
  TH1D *Hpiminuseta_p = new TH1D("HPiminuseta_p", "Titel; p [GeV/c]; Anzahl", 100, 0, 200);
  TH1D *Hpipluseta_p = new TH1D("HPipluseta_p", "Titel; p [GeV/c]; Anzahl", 100, 0, 200);
  TH1D *Hpipiminuseta_m =
      new TH1D("Hpipiminuseta_m", "Title; m_{x} [GeV/c^{2}]; Anzahl", 500, 0, 3);
  TH1D *Hpipluspiminuseta_m =
      new TH1D("Hpipluspiminuseta_m", "Title; m_{x} [GeV/c^{2}]; Anzahl", 500, 0, 3);

  // Read in real data
  std::string File = "/mnt/data/compass/2009/3Pi2G_Eventselection_Slot5/EventSelection_2023/"
                     "PiPiPiEta/Total_08_09_July05.root";
  TFile f(File.c_str(), "read");

  if (!f.IsOpen()) {
    cout << "Could not open " << File << "!" << endl;
    abort();
  }

  TTreeReader TreeReader("AfterSelection", &f); // !!!AfterSelection PseudoData_MC
  if (TreeReader.GetEntries(true) == -1) {
    cout << "There is no tree called \"AfterSelection\" in file \"" << File << "\"!" << endl;
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

  // TH1D *Htp = new TH1D("Htp", "Title; t'; Anzahl", 10, 0, 5);

  for (auto i = TreeReader.begin(); i != TreeReader.end(); ++i) {

    // if (counter % (TreeReader.GetEntries() / 100) == 0) {
    //   cout << "... " << int(100.0 * double(counter) / double(TreeReader.GetEntries()))
    //        << " % of the real data have been read ..." << endl;
    // }

    TreeReader.SetEntry(*i);

    std::pair<double, double> tp_mx = std::make_pair(*t_prime, *mx);
    // cout << "m_x= " << tp_mx.second << endl;
    Hmx->Fill(tp_mx.second); // m_x
    Htp->Fill(tp_mx.first);  // t_prime
    H_mx_tp->Fill(tp_mx.second, tp_mx.first);

    TLorentzVector Beam(BeamPion->at(0), BeamPion->at(1), BeamPion->at(2), BeamPion->at(3));
    TLorentzVector PiMinus1(PiMinus_1->at(0), PiMinus_1->at(1), PiMinus_1->at(2), PiMinus_1->at(3));
    TLorentzVector PiMinus2(PiMinus_2->at(0), PiMinus_2->at(1), PiMinus_2->at(2), PiMinus_2->at(3));
    TLorentzVector PiPlus(PiPlus_->at(0), PiPlus_->at(1), PiPlus_->at(2), PiPlus_->at(3));
    TLorentzVector Eta(Eta_->at(0), Eta_->at(1), Eta_->at(2), Eta_->at(3));

    TLorentzVector PiPiPi = PiMinus1 + PiMinus2 + PiPlus;

    Hpipipi_mx->Fill(PiPiPi.M());
    Hpipipi_p->Fill(PiPiPi.P());
    Heta_p->Fill(Eta.P());

    TLorentzVector PiPi_1 = PiMinus1 + PiPlus;
    TLorentzVector PiPi_2 = PiMinus2 + PiPlus;

    Hpipi_m->Fill(PiPi_1.M());
    Hpipi_m->Fill(PiPi_2.M());

    TLorentzVector PiminusEta = PiMinus1 + Eta;
    Hpiminuseta_m->Fill(PiminusEta.M());
    Hpiminuseta_p->Fill(PiminusEta.P());

    TLorentzVector PiplusEta = PiPlus + Eta;
    Hpipluseta_m->Fill(PiplusEta.M());
    Hpipluseta_p->Fill(PiplusEta.P());

    TLorentzVector PiPiminusEta = PiMinus1 + PiMinus2 + Eta;
    Hpipiminuseta_m->Fill(PiPiminusEta.M());
    TLorentzVector PiplusPiminusEta = PiMinus1 + PiPlus + Eta;
    Hpipluspiminuseta_m->Fill(PiplusPiminusEta.M());

    counter++;
  }
  f.Close();
  // TCanvas *C1 = new TCanvas("C1", "", 1000, 1000);
  // C1->cd();
  // Htp->Draw();

  // TCanvas *c2 = new TCanvas("c2", "", 1000, 1000);
  // c2->cd();
  // Hmx->Draw();

  // TCanvas *c3 = new TCanvas("c3", "", 1000, 1000);
  // c3->cd();
  // c3->SetLogz();
  // H_mx_tp->Draw("colz");

  // TCanvas *c4 = new TCanvas("c4", "PiPiPi_m", 1000, 1000);
  // c4->cd();
  // Hpipipi_mx->Draw();

  // TCanvas *c5 = new TCanvas("c5", "PiPiPi_Impuls", 1000, 1000);
  // c5->cd();
  // Hpipipi_p->Draw();

  // TCanvas *c6 = new TCanvas("c6", "Eta_p", 1000, 1000);
  // c6->cd();
  // Heta_p->Draw();

  TCanvas *c7 = new TCanvas("c7", "pi_plus_pi_minus", 1000, 1000);
  c7->cd();
  Hpipi_m->Draw();

  // TCanvas *c8 = new TCanvas("c8", "pi_minus/plus_eta", 1000, 1000);
  // c8->cd();
  // Hpiminuseta_m->Draw();
  // Hpipluseta_m->SetLineColor(kRed);
  // Hpipluseta_m->Draw("same");
  // TLegend *legend = new TLegend(0.6, 0.7, 0.85, 0.85);
  // legend->SetTextFont(132);
  // legend->AddEntry(Hpiminuseta_m, "#pi^{#minus}#eta", "l");
  // legend->AddEntry(Hpipluseta_m, "#pi^{#plus}#eta", "l");
  // legend->Draw("same");

  TCanvas *c9 = new TCanvas("c8", "pi_minus/plus_pi_minus_eta", 1000, 1000);
  c9->cd();
  Hpipluspiminuseta_m->SetLineColor(kRed);

  Hpipluspiminuseta_m->Draw();

  Hpipiminuseta_m->Draw("same");
  TLegend *legend2 = new TLegend(0.6, 0.7, 0.8, 0.8);
  legend2->SetTextFont(132);
  legend2->AddEntry(Hpipiminuseta_m, "#pi^{#minus}#pi^{#minus}#eta", "l");
  legend2->AddEntry(Hpipluspiminuseta_m, "#pi^{#plus}#pi^{#minus}#eta", "l");
  legend2->Draw();

  return;
}