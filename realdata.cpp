
#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TLorentzVector.h"
#include "TTreeReader.h"
#include <iostream>
#include <sstream>
using namespace std;

void realdata() {
  // Read File
  std::string File = "/mnt/data/compass/2009/3Pi2G_Eventselection_Slot5/EventSelection_2023/"
                     "PiPiPiEta/Total_08_09_June29.root";
  TFile f(File.c_str(), "read");
  std::cout << "open file" << std::endl;

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

  TH1D *Hist1 = new TH1D(
      "pimpip",
      "Title; #it{m_{#pi^{#minus}#pi^{#plus}}} [GeV/#it{c}^{2}]; Anzahl / (4MeV/#it{c}^{2})", 1000,
      0, 3);

  TH1D *Hist2 = new TH1D(
      "pimeta", "Title; #it{m_{#pi^{#minus}#eta}} [GeV/#it{c}^{2}]; Anzahl / (4MeV/#it{c}^{2})",
      1000, 0, 3);
  TH1D *Hist3 = new TH1D(
      "pimpipeta",
      "Title; #it{m_{#pi^{#minus}#pi^{#plus}#eta}} [GeV/#it{c}^{2}]; Anzahl / (4MeV/#it{c}^{2})",
      1000, 0, 3);
  TH1D *Hist4 = new TH1D("pimpippim",
                         "Title; #it{m_{#pi^{#minus}#pi^{#plus}#pi^{#minus}}} [GeV/#it{c}^{2}]; "
                         "Anzahl / (4MeV/#it{c}^{2})",
                         1000, 0, 3);

  std::vector<TH1D *> mHist1, mHist2, mHist3, mHist4;

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
    if (tp_mx.second <= 3.0) {
      if (tp_mx.first < 1.0 && tp_mx.first > 0.1) {

        TLorentzVector Beam(BeamPion->at(0), BeamPion->at(1), BeamPion->at(2), BeamPion->at(3));
        TLorentzVector PiMinus1(PiMinus_1->at(0), PiMinus_1->at(1), PiMinus_1->at(2),
                                PiMinus_1->at(3));
        TLorentzVector PiMinus2(PiMinus_2->at(0), PiMinus_2->at(1), PiMinus_2->at(2),
                                PiMinus_2->at(3));
        TLorentzVector PiPlus(PiPlus_->at(0), PiPlus_->at(1), PiPlus_->at(2), PiPlus_->at(3));
        TLorentzVector Eta(Eta_->at(0), Eta_->at(1), Eta_->at(2), Eta_->at(3));

        // Pi+Pi-Eta
        TLorentzVector PiminusPiplus = PiMinus1 + PiPlus;
        TLorentzVector PiminusEta = PiMinus1 + Eta;

        TLorentzVector PiminusPiplusEta = PiMinus1 + PiPlus + Eta;
        TLorentzVector PiminusPiplusPiminus = PiMinus1 + PiPlus + PiMinus2;

        Hist1->Fill(PiminusPiplus.M());
        Hist2->Fill(PiminusEta.M());
        Hist3->Fill(PiminusPiplusEta.M());
        Hist4->Fill(PiminusPiplusPiminus.M());
      }
    }
    counter++;
  }
  std::cout << "Draw Histogram" << endl;
  TCanvas *c1 = new TCanvas("c1", "pi-pi+", 1000, 1000);
  c1->cd();

  Hist1->Rebin(4);
  Hist2->Rebin(4);
  Hist3->Rebin(4);
  Hist4->Rebin(4);
  c1->cd();

  Hist1->Draw();
  c1->Print("piminuspiplus.pdf", "EmbedFonts");

  TCanvas *c2 = new TCanvas("c2", "pi-eta", 1000, 1000);
  c2->cd();
  Hist2->Draw();
  c2->Print("piminuseta.pdf", "EmbedFonts");

  TCanvas *c3 = new TCanvas("c3", "pi-pi+eta", 1000, 1000);
  c3->cd();
  Hist3->Draw();
  c3->Print("piminuspipluseta.pdf", "EmbedFonts");

  TCanvas *c4 = new TCanvas("c4", "pi-pi+pi-", 1000, 1000);
  c4->cd();
  Hist4->Draw();
  c4->Print("piminuspipluspiminus.pdf", "EmbedFonts");

  f.Close();
}
