// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include <fstream>
#include "TCanvas.h"
#include "TH1F.h"
#include "TFile.h"
#include "TApplication.h"
#include "TGraphErrors.h"

#include "MCHBase/Digit.h"
#include "MCHBase/PreCluster.h"
#include "TBDigitsFileReader.h"
#include "MCHPreClustering/PreClusterFinder.h"
#include "MCHClustering/ClusteringForTest.h"

using namespace o2::mch;
using namespace std;

int main(int argc, char** argv)
{

  TBDigitsFileReader digitsReader;
  digitsReader.init(argv[1]);
  PreClusterFinder preClusterFinder;
  Clustering clustering;
  float residualsMoi[7000];
  float residualsAlberto[7000];
  float residualsDifference[7000];
  float ytrks[7000];
  int count = 0;

  //   TFile *f = new TFile("ComparisonSimuTB/ReferenceTB.root", "NEW");

  TApplication app("app", &argc, argv);

  TCanvas* cResDistrib = new TCanvas("cResDistrib", "Residuals Distribution", 0, 0, 600, 600);
  TH1F* hResDist_Seb = new TH1F("hResDist_Seb", "Residuals distribution from TB data - My Clustering", 50, -0.1, 0.1);
  TH1F* hResDist_Alberto = new TH1F("hResDist_Alberto", "Residuals distribution from TB data - Alberto Clustering", 50, -0.1, 0.1);
  TH1F* hTB_ClusterChg = new TH1F("hTB_ClusterChg", "TB data cluster charge distribution", 300, 0, 3000);
  TH1F* hTB_MeanYNHits = new TH1F("hTB_MeanYNHits", "Mean YNHits as a function of charge", 300, 0, 3000);
  TH1F* hTB_MeanNHits = new TH1F("hTB_MeanNHits", "Mean NHits as a function of charge", 300, 0, 3000);

  cResDistrib->Divide(2, 1);

  TCanvas* cYNhits = new TCanvas("cYNhits", "cYNhits", 0, 0, 600, 600);
  TH1F* hYNhitsTB = new TH1F("hYNhitsTB", "TB data NHitsY distribution", 10, 0.5, 10.5);

  TCanvas* cXNhits = new TCanvas("cXNhits", "cXNhits", 0, 0, 600, 600);
  TH1F* hXNhitsTB = new TH1F("hXNhitsTB", "TB data NHitsX distribution", 10, 0.5, 10.5);

  TCanvas* cNhits = new TCanvas("cNhits", "cNhits", 0, 0, 600, 600);
  TH1F* hNhitsTB = new TH1F("hNhitsTB", "TB data NHits distribution", 12, 0.5, 12.5);

  TCanvas* cchmax = new TCanvas("cchmax", "cchmax", 0, 0, 600, 600);
  TH1F* hchmax = new TH1F("hchmax", "TB data charge max on a pad distribution", 300, 0, 3000);

  TCanvas* cytrkspread = new TCanvas("cytrkspread", "cytrkspread", 0, 0, 600, 600);
  TH1F* hytrkspread = new TH1F("hytrkspread", "Spread of y trk", 500, 10, 15);

  TCanvas* cChargeIntervals = new TCanvas("cChargeIntervals", "NHits wrt K3 for different charge intervals", 0, 0, 600, 600);
  TH1F* hNhitsTB0_300 = new TH1F("hNhitsTB0_300", "TB data NHits distribution - Charge 0 to 300 ADC", 8, 0.5, 8.5);
  TH1F* hNhitsTB300_600 = new TH1F("hNhitsTB300_600", "TB data NHits distribution - Charge 300 to 600 ADC", 8, 0.5, 8.5);
  TH1F* hNhitsTB600_1000 = new TH1F("hNhitsTB600_1000", "TB data NHits distribution - Charge 600 to 1000 ADC", 8, 0.5, 8.5);
  TH1F* hNhitsTB1000_3000 = new TH1F("hNhitsTB1000_3000", "TB data NHits distribution - Charge 1000 to 3000 ADC", 8, 0.5, 8.5);

  TCanvas* cChargeIntervals_Y = new TCanvas("cChargeIntervals_Y", "NYHits wrt K3 for different charge intervals", 0, 0, 600, 600);
  TH1F* hYNhitsTB0_300 = new TH1F("hYNhitsTB0_300", "TB data NHitsY distribution - Charge 0 to 300 ADC", 8, 0.5, 8.5);
  TH1F* hYNhitsTB300_600 = new TH1F("hYNhitsTB300_600", "TB data NHitsY distribution - Charge 300 to 600 ADC", 8, 0.5, 8.5);
  TH1F* hYNhitsTB600_1000 = new TH1F("hYNhitsTB600_1000", "TB data NHitsY distribution - Charge 600 to 1000 ADC", 8, 0.5, 8.5);
  TH1F* hYNhitsTB1000_3000 = new TH1F("hYNhitsTB1000_3000", "TB data NHitsY distribution - Charge 1000 to 3000 ADC", 8, 0.5, 8.5);

  preClusterFinder.init();

  Digit* digitsBuffer = NULL;
  std::vector<Digit> digits(0);
  std::vector<PreCluster> preClusters(0);
  std::vector<Clustering::Cluster> clusters(0);

  // load digits from binary input file, block-by-block
  while (digitsReader.readDigitsFromFile()) {
    clusters.clear();

    // get number of loaded digits and store them into a memory buffer
    auto nDigits = digitsReader.getNumberOfDigits();
    float xtrk;
    float ytrk;
    digitsReader.get_trk_pos(819, xtrk, ytrk);
    cout << "xtrk: " << xtrk << endl;
    cout << "ytrk: " << ytrk << endl;
    printf("nDigits: %d\n", (int)nDigits);
    hytrkspread->Fill(ytrk);

    std::vector<TBCluster>& tbclusters = digitsReader.getClusters();
    std::cout << "tbclusters.size(): " << tbclusters.size() << std::endl;
    for (size_t ic = 0; ic < tbclusters.size(); ic++) {
      std::cout << "TBCluster " << ic << ":\n"
                << tbclusters[ic] << std::endl;
    }
    //continue;

    digitsBuffer = (Digit*)realloc(digitsBuffer, sizeof(Digit) * nDigits);
    digitsReader.storeDigits(digitsBuffer);

    // load the digits from the memory buffer and run the pre-clustering phase
    preClusterFinder.reset();
    preClusterFinder.loadDigits({digitsBuffer, nDigits});
    preClusterFinder.run();

    // get the preclusters and associated digits
    preClusters.clear();
    digits.clear();
    preClusterFinder.getPreClusters(preClusters, digits);

    printf("\n\n==========\nRunning Clustering\n\n");

    // Fit Mathieson
    //clustering.runFinderSimpleFit(preClusters, digits, clusters);

    //Uncomment the method you wish to use to clusterize

    //Runs the clustering of preClusters following a CenterOfGravity algorithm. Fills clusters.
    //clustering.runFinderCOG(preClusters, digits, clusters);
    //          printf("Number of clusters obtained and saved: %lu\n", clusters.size());

    // Fit Mathieson
    //       clustering.runFinderSimpleFit(preClusters, clusters);

    // Fit Simple Gaussienne
    //     clustering.runFinderGaussianFit(preClusters, clusters);

    // Fit Double Gaussienne
    //     clustering.runFinderDoubleGaussianFit(preClusters, clusters);

    // Choose to only look at nice events, where we have one precluster reconstructed as a single cluster, which are a majority (972 out of roughly 1000). This helps us get rid of pathological events.

    if (preClusters.size() == 1 && clusters.size() == 1 && tbclusters.size() == 1) {
      hNhitsTB->Fill(tbclusters[0].fNhits);
      hYNhitsTB->Fill(tbclusters[0].fYNhits);
      hXNhitsTB->Fill(tbclusters[0].fXNhits);
      hchmax->Fill(tbclusters[0].fChargemax);
      float yobtenuMoi = clusters[0].gety();
      float yobtenuAlberto = tbclusters[0].fYclus;
      float chargeAlberto = tbclusters[0].fCharge;
      float differenceMoi = ytrk - yobtenuMoi;
      float differenceAlberto = ytrk - yobtenuAlberto;
      ytrks[count] = ytrk;
      residualsMoi[count] = differenceMoi;
      residualsAlberto[count] = differenceAlberto;
      residualsDifference[count] = differenceMoi - differenceAlberto;

      if (chargeAlberto < 300) {
        hNhitsTB0_300->Fill(tbclusters[0].fNhits);
        hYNhitsTB0_300->Fill(tbclusters[0].fYNhits);
      } else if (300 <= chargeAlberto && chargeAlberto < 600) {
        hNhitsTB300_600->Fill(tbclusters[0].fNhits);
        hYNhitsTB300_600->Fill(tbclusters[0].fYNhits);
      } else if (600 <= chargeAlberto && chargeAlberto < 1000) {
        hNhitsTB600_1000->Fill(tbclusters[0].fNhits);
        hYNhitsTB600_1000->Fill(tbclusters[0].fYNhits);
      } else if (1000 <= chargeAlberto && chargeAlberto < 3000) {
        hNhitsTB1000_3000->Fill(tbclusters[0].fNhits);
        hYNhitsTB1000_3000->Fill(tbclusters[0].fYNhits);
      }

      // Look at the distribution of residuals from this TestBeam data

      hResDist_Seb->Fill(differenceMoi);
      hResDist_Seb->GetXaxis()->SetTitle("Residual y (cm)");
      hResDist_Seb->GetYaxis()->SetTitle("Count");
      hResDist_Alberto->Fill(differenceAlberto);
      hResDist_Alberto->GetXaxis()->SetTitle("Residual y (cm)");
      hResDist_Alberto->GetYaxis()->SetTitle("Count");
      hTB_ClusterChg->Fill(chargeAlberto);
      hTB_ClusterChg->GetXaxis()->SetTitle("Cluster charge (ADC)");
      hTB_ClusterChg->GetYaxis()->SetTitle("Count");
      hTB_MeanYNHits->Fill(chargeAlberto, tbclusters[0].fYNhits);
      hTB_MeanYNHits->GetXaxis()->SetTitle("Cluster charge (ADC)");
      hTB_MeanYNHits->GetYaxis()->SetTitle("Mean number of pads fired on Y");
      hTB_MeanNHits->Fill(chargeAlberto, tbclusters[0].fNhits);
      hTB_MeanNHits->GetXaxis()->SetTitle("Cluster charge (ADC)");
      hTB_MeanNHits->GetYaxis()->SetTitle("Mean number of pads fired");

      cout << "RESIDUAL y found by me: " << differenceMoi << endl;
      cout << "RESIDUAL y found by Alberto: " << differenceAlberto << endl;
      if (residualsDifference[count] != 0) {
        cout << "DISAGREE ON RESIDUALS" << endl;
        cout << "Number of preclusters found by Alberto: " << tbclusters.size() << endl;
      }
      cout << "It was cluster number " << count << endl;
      count++;
    }

    break;
  }
  cResDistrib->cd(1);
  hResDist_Seb->Draw();
  cResDistrib->cd(2);
  hResDist_Alberto->Draw();
  cResDistrib->Update();
  cResDistrib->Draw();

  cYNhits->cd();
  hYNhitsTB->Draw();
  cYNhits->Update();
  cYNhits->Draw();

  cXNhits->cd();
  hXNhitsTB->Draw();
  cXNhits->Update();
  cXNhits->Draw();

  cNhits->cd();
  hNhitsTB->Draw();
  cNhits->Update();
  cNhits->Draw();

  cchmax->cd();
  hchmax->Draw();
  cchmax->Update();
  cchmax->Draw();

  cytrkspread->cd();
  hytrkspread->Draw();
  cytrkspread->Update();
  cytrkspread->Draw();

  TCanvas* cMeanHits_Y = new TCanvas("cMeanHits_Y", "As a function of charge", 0, 0, 600, 600);
  cMeanHits_Y->Divide(2, 1);
  cMeanHits_Y->cd(1);
  hTB_ClusterChg->Draw();

  double custombins[5]{0, 300, 600, 1000, 3000};
  TH1F* hTB_ClusterChg_60bins = (TH1F*)hTB_ClusterChg->Rebin(5, "hTB_ClusterChg_60bins");
  TH1F* hTB_ClusterChg_custombins = (TH1F*)hTB_ClusterChg->Rebin(4, "hTB_ClusterChg_custombins", custombins);

  cMeanHits_Y->cd(2);
  TH1F* hTB_MeanYNHits_60bins = (TH1F*)hTB_MeanYNHits->Rebin(5, "hTB_MeanYNHits_60bins");
  TH1F* hTB_MeanYNHits_custombins = (TH1F*)hTB_MeanYNHits->Rebin(4, "hhTB_MeanYNHits_custombins", custombins);
  hTB_MeanYNHits->Divide(hTB_ClusterChg);
  hTB_MeanYNHits->Draw();
  cMeanHits_Y->Update();
  cMeanHits_Y->Draw();

  TCanvas* cMeanHits = new TCanvas("cMeanHits", "As a function of charge", 0, 0, 600, 600);
  cMeanHits->Divide(2, 1);
  cMeanHits->cd(1);
  hTB_ClusterChg->Draw();
  cMeanHits->cd(2);
  TH1F* hTB_MeanNHits_60bins = (TH1F*)hTB_MeanNHits->Rebin(5, "hTB_MeanNHits_60bins");
  TH1F* hTB_MeanNHits_custombins = (TH1F*)hTB_MeanNHits->Rebin(4, "hhTB_MeanNHits_custombins", custombins);
  hTB_MeanNHits->Divide(hTB_ClusterChg);
  hTB_MeanNHits->Draw();
  cMeanHits->Update();
  cMeanHits->Draw();

  TCanvas* cMeanHits_Y_60bins = new TCanvas("cMeanHits_Y_60bins", "As a function of charge", 0, 0, 600, 600);
  TCanvas* cMeanHits_60bins = new TCanvas("cMeanHits_60bins", "As a function of charge", 0, 0, 600, 600);

  cMeanHits_Y_60bins->Divide(2, 1);
  cMeanHits_Y_60bins->cd(1);
  hTB_ClusterChg_60bins->Draw();
  cMeanHits_Y_60bins->cd(2);
  hTB_MeanYNHits_60bins->Divide(hTB_ClusterChg_60bins);
  hTB_MeanYNHits_60bins->Draw();
  cMeanHits_Y_60bins->Update();
  cMeanHits_Y_60bins->Draw();

  cMeanHits_60bins->Divide(2, 1);
  cMeanHits_60bins->cd(1);
  hTB_ClusterChg_60bins->Draw();
  cMeanHits_60bins->cd(2);
  hTB_MeanNHits_60bins->Divide(hTB_ClusterChg_60bins);
  hTB_MeanNHits_60bins->Draw();
  cMeanHits_60bins->Update();
  cMeanHits_60bins->Draw();

  TCanvas* cMeanHits_Y_custombins = new TCanvas("cMeanHits_Y_custombins", "As a function of charge", 0, 0, 600, 600);
  TCanvas* cMeanHits_custombins = new TCanvas("cMeanHits_custombins", "As a function of charge", 0, 0, 600, 600);

  cMeanHits_Y_custombins->Divide(2, 1);
  cMeanHits_Y_custombins->cd(1);
  hTB_ClusterChg_custombins->Draw();
  cMeanHits_Y_custombins->cd(2);
  hTB_MeanYNHits_custombins->Divide(hTB_ClusterChg_custombins);
  hTB_MeanYNHits_custombins->Draw();
  cMeanHits_Y_custombins->Update();
  cMeanHits_Y_custombins->Draw();

  cMeanHits_custombins->Divide(2, 1);
  cMeanHits_custombins->cd(1);
  hTB_ClusterChg_custombins->Draw();
  cMeanHits_custombins->cd(2);
  hTB_MeanNHits_custombins->Divide(hTB_ClusterChg_custombins);
  hTB_MeanNHits_custombins->Draw();
  cMeanHits_custombins->Update();
  cMeanHits_custombins->Draw();

  cChargeIntervals->Divide(2, 2);
  cChargeIntervals->cd(1);
  hNhitsTB0_300->GetXaxis()->SetTitle("Number of pads fired");
  hNhitsTB0_300->GetYaxis()->SetTitle("Count");
  hNhitsTB0_300->Draw();
  cChargeIntervals->cd(2);
  hNhitsTB300_600->GetXaxis()->SetTitle("Number of pads fired");
  hNhitsTB300_600->GetYaxis()->SetTitle("Count");
  hNhitsTB300_600->Draw();
  cChargeIntervals->cd(3);
  hNhitsTB600_1000->GetXaxis()->SetTitle("Number of pads fired");
  hNhitsTB600_1000->GetYaxis()->SetTitle("Count");
  hNhitsTB600_1000->Draw();
  cChargeIntervals->cd(4);
  hNhitsTB1000_3000->GetXaxis()->SetTitle("Number of pads fired");
  hNhitsTB1000_3000->GetYaxis()->SetTitle("Count");
  hNhitsTB1000_3000->Draw();
  cChargeIntervals->Update();
  cChargeIntervals->Draw();

  cChargeIntervals_Y->Divide(2, 2);
  cChargeIntervals_Y->cd(1);
  hYNhitsTB0_300->GetXaxis()->SetTitle("Number of pads fired on Y");
  hYNhitsTB0_300->GetYaxis()->SetTitle("Count");
  hYNhitsTB0_300->Draw();
  cChargeIntervals_Y->cd(2);
  hYNhitsTB300_600->GetXaxis()->SetTitle("Number of pads fired on Y");
  hYNhitsTB300_600->GetYaxis()->SetTitle("Count");
  hYNhitsTB300_600->Draw();
  cChargeIntervals_Y->cd(3);
  hYNhitsTB600_1000->GetXaxis()->SetTitle("Number of pads fired on Y");
  hYNhitsTB600_1000->GetYaxis()->SetTitle("Count");
  hYNhitsTB600_1000->Draw();
  cChargeIntervals_Y->cd(4);
  hYNhitsTB1000_3000->GetXaxis()->SetTitle("Number of pads fired on Y");
  hYNhitsTB1000_3000->GetYaxis()->SetTitle("Count");
  hYNhitsTB1000_3000->Draw();
  cChargeIntervals_Y->Update();
  cChargeIntervals_Y->Draw();

  // Look at the residuals wrt coordinate y of the hit, to see if there is a bias

  TCanvas* cerr = new TCanvas("cerr", "Graph example", 0, 0, 600, 600);
  cerr->Divide(2, 1);

  cerr->cd(1);
  TGraph* gr1 = new TGraph(count, ytrks, residualsMoi);
  gr1->SetTitle("Residuals wrt input y - My Clusetring");
  gr1->SetMarkerColor(4);
  gr1->SetLineColor(4);
  gr1->SetMarkerStyle(8);
  gr1->GetXaxis()->SetTitle("y input (cm)");
  gr1->GetYaxis()->SetTitle("Residual (cm)");
  gr1->Draw("AP");

  cerr->cd(2);
  TGraph* gr2 = new TGraph(count, ytrks, residualsAlberto);
  gr2->SetTitle("Residuals wrt input y - Alberto's Clustering");
  gr2->SetMarkerColor(4);
  gr2->SetLineColor(4);
  gr2->SetMarkerStyle(8);
  gr2->GetXaxis()->SetTitle("y input (cm)");
  gr2->GetYaxis()->SetTitle("Residual (cm)");
  gr2->Draw("AP");

  cerr->Update();
  cerr->Draw();

  TCanvas* cdiff = new TCanvas("cdiff", "Graph example", 0, 0, 600, 600);

  TGraph* gr3 = new TGraph(count, ytrks, residualsDifference);
  gr3->SetTitle("Difference in residuals");
  gr3->SetMarkerColor(4);
  gr3->SetLineColor(4);
  gr3->SetMarkerStyle(8);
  gr3->GetXaxis()->SetTitle("y input (cm)");
  gr3->GetYaxis()->SetTitle("Residual difference (Mine-Alberto's) (cm)");
  gr3->Draw("AP");

  cdiff->Update();
  cdiff->Draw();

  //    f->Write();
  //    f->Close();

  app.Run(kTRUE);

  return 0;
}
