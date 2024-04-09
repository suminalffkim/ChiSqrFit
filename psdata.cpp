
#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TLorentzVector.h"
#include "TTreeReader.h"
#include <iostream>
#include <sstream>
using namespace std;

void psdata() {
  // Read File
  std::string File = "/home/kim/Desktop/PiPiPiGG_PWA/rootfiles/PseudoDataFile_Resonances.root";
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

  TTreeReaderValue<double> mx(TreeReader, "m_x");
  TTreeReaderValue<double> t_prime(TreeReader, "tprime");

  int counter = 1;

  TH1D *Hist1 = new TH1D("mX", "Title; #it{m_{X}} [GeV/#it{c}^{2}]; Intensity / (40 MeV/c^{2})",
                         200, 1.2, 2.52);

  TH1D *Hist2 =
      new TH1D("tpirme", "Title; #it{t'} [(GeV/#it{c})^{2}]; Anzahl / (40 MeV/c)^{2}", 500, 0, 1.1);

  std::vector<TH1D *> mHist1, mHist2;

  for (auto i = TreeReader.begin(); i != TreeReader.end(); ++i) {

    // if (counter % (TreeReader.GetEntries() / 100) == 0) {
    //   cout << "... " << int(100.0 * double(counter) / double(TreeReader.GetEntries()))
    //        << " % of the real data have been read ..." << endl;
    // }

    TreeReader.SetEntry(*i);

    std::pair<double, double> tp_mx = std::make_pair(*t_prime, *mx);
    // cout << "m_x= " << tp_mx.second << endl;
    Hist1->Fill(tp_mx.second); // m_x
    Hist2->Fill(tp_mx.first);  // t_prime
    // H_mx_tp->Fill(tp_mx.second, tp_mx.first);

    counter++;
  }

  Hist1->Rebin(8);
  Hist2->Rebin(4);

  std::cout << "Draw Histogram" << endl;
  TCanvas *c1 = new TCanvas("c1", "mX", 1000, 1000);
  c1->cd();

  Hist1->Draw();
  c1->Print("mXPhaseSpace.pdf", "EmbedFonts");

  TCanvas *c2 = new TCanvas("c2", "tp", 1000, 1000);
  c2->cd();
  Hist2->Draw();
  c2->Print("tpPhaseSpace.pdf", "EmbedFonts");

  f.Close();
}
