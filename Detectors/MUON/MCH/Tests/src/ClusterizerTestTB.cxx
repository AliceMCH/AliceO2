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
#include "TApplication.h"
#include "TGraphErrors.h"

#include "MCHBase/Digit.h"
#include "MCHBase/PreClusterBlock.h"
#include "TBDigitsFileReader.h"
#include "../../PreClustering/src/PreClusterFinder.h"
#include "MCHClustering/ClusteringForTest.h"

using namespace o2::mch;
using namespace std;

int main(int argc, char** argv)
{
    
  TBDigitsFileReader digitsReader(argv[1]);
  PreClusterFinder preClusterFinder;
  Clustering clustering;
  float residualsMoi[7000];
    float residualsAlberto[7000];
    float residualsDifference[7000];
    float ytrks[7000];
    int count = 0;
    
    TApplication app ("app",&argc,argv);
    
    TCanvas *cbell = new TCanvas("cbell","Bell",0,0,600,600);
    TH1F *h1 = new TH1F("h1", "Residuals distribution from TB data - My Clustering", 50, -0.1, 0.1);
    TH1F *h2 = new TH1F("h2", "Residuals distribution from TB data - Alberto Clustering", 50, -0.1, 0.1);
    TH1F *h3 = new TH1F("h3", "TB data cluster charge distribution", 60, 0, 3000);
    TH1F *h4 = new TH1F("h4", "Mean YNHits as a function of charge", 60, 0, 3000);
    TH1F *h5 = new TH1F("h5", "Mean NHits as a function of charge", 60, 0, 3000);
    
    cbell->Divide(2,1);
    
    TCanvas *cYNhits = new TCanvas("cYNhits","cYNhits",0,0,600,600);
    TH1F *hYNhitsTB = new TH1F("hYNhitsTB", "TB data NHitsY distribution", 10, 0.5, 10.5);
    
    TCanvas *cXNhits = new TCanvas("cXNhits","cXNhits",0,0,600,600);
    TH1F *hXNhitsTB = new TH1F("hXNhitsTB", "TB data NHitsX distribution", 10, 0.5, 10.5);
    
    TCanvas *cNhits = new TCanvas("cNhits","cNhits",0,0,600,600);
    TH1F *hNhitsTB = new TH1F("hNhitsTB", "TB data NHits distribution", 12, 0.5, 12.5);
    
    TCanvas *cchmax = new TCanvas("cchmax","cchmax",0,0,600,600);
    TH1F *hchmax = new TH1F("hchmax", "TB data charge max on a pad distribution", 300, 0, 3000);
    
    TCanvas *cytrkspread = new TCanvas("cytrkspread","cytrkspread",0,0,600,600);
    TH1F *hytrkspread = new TH1F("hytrkspread", "Spread of y trk", 500, 10, 15);
  
  preClusterFinder.init();

  Digit* digitsBuffer = NULL;
  std::vector<Digit> digits(0);
  std::vector<PreClusterStruct> preClusters(0);
  std::vector<Clustering::Cluster> clusters(0);
    
  // load digits from binary input file, block-by-block
  while(digitsReader.readDigitsFromFile()) {
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
    std::cout<<"tbclusters.size(): "<<tbclusters.size()<<std::endl;
    for(size_t ic = 0; ic < tbclusters.size(); ic++) {
      std::cout << "TBCluster " << ic << ":\n" << tbclusters[ic] << std::endl;
    }
    //continue;

    digitsBuffer = (Digit*)realloc(digitsBuffer, sizeof(Digit) * nDigits);
    digitsReader.storeDigits(digitsBuffer);

    // load the digits from the memory buffer and run the pre-clustering phase
    preClusterFinder.reset();
    preClusterFinder.loadDigits(digitsBuffer, nDigits);
    preClusterFinder.run();

    // get the preclusters and associated digits
    preClusterFinder.getPreClusters(preClusters, digits);
      
      printf("\n\n==========\nRunning Clustering\n\n");

      
      
       //Uncomment the method you wish to use to clusterize
            
          //Runs the clustering of preClusters following a CenterOfGravity algorithm. Fills clusters.
          clustering.runFinderCOG(preClusters, clusters);
//          printf("Number of clusters obtained and saved: %lu\n", clusters.size());
            
            // Fit Mathieson
  //       clustering.runFinderSimpleFit(preClusters, clusters);
            
            // Fit Simple Gaussienne
       //     clustering.runFinderGaussianFit(preClusters, clusters);
            
            // Fit Double Gaussienne
      //     clustering.runFinderDoubleGaussianFit(preClusters, clusters);

      
      
      
      
      // Choose to only look at nice events, where we have one precluster reconstructed as a single cluster, which are a majority (972 out of roughly 1000). This helps us get rid of pathological events.
      
      if(preClusters.size()==1 && clusters.size() == 1 && tbclusters.size() ==1){
          hNhitsTB->Fill(tbclusters[0].fNhits);
          hYNhitsTB->Fill(tbclusters[0].fYNhits);
          hXNhitsTB->Fill(tbclusters[0].fXNhits);
          hchmax->Fill(tbclusters[0].fChargemax);
          float yobtenuMoi = clusters[0].gety();
          float yobtenuAlberto = tbclusters[0].fYclus;
          float chargeAlberto = tbclusters[0].fCharge;
          float differenceMoi = ytrk-yobtenuMoi;
          float differenceAlberto = ytrk - yobtenuAlberto;
          ytrks[count] = ytrk;
          residualsMoi[count] = differenceMoi;
          residualsAlberto[count] = differenceAlberto;
          residualsDifference[count] = differenceMoi-differenceAlberto;
          
          
          // Look at the distribution of residuals from this TestBeam data
          
            h1->Fill(differenceMoi);
            h1->GetXaxis()->SetTitle("Residual y (cm)");
            h1->GetYaxis()->SetTitle("Count");
              h2->Fill(differenceAlberto);
              h2->GetXaxis()->SetTitle("Residual y (cm)");
              h2->GetYaxis()->SetTitle("Count");
              h3->Fill(chargeAlberto);
              h3->GetXaxis()->SetTitle("Cluster charge (ADC)");
              h3->GetYaxis()->SetTitle("Count");
              h4->Fill(chargeAlberto, tbclusters[0].fYNhits);
              h4->GetXaxis()->SetTitle("Cluster charge (ADC)");
              h4->GetYaxis()->SetTitle("Mean number of pads fired on Y");
          h5->Fill(chargeAlberto, tbclusters[0].fNhits);
          h5->GetXaxis()->SetTitle("Cluster charge (ADC)");
          h5->GetYaxis()->SetTitle("Mean number of pads fired");
        
          cout << "RESIDUAL y found by me: " << differenceMoi <<endl;
          cout << "RESIDUAL y found by Alberto: " << differenceAlberto <<endl;
          if(residualsDifference[count] != 0){
              cout << "DISAGREE ON RESIDUALS" <<endl;
              cout << "Number of preclusters found by Alberto: " << tbclusters.size() <<endl;
          }
          cout << "It was cluster number " << count <<endl;
          count++;
      }
      
    //break;
  }
    cbell->cd(1);
    h1->Draw();
    cbell->cd(2);
    h2->Draw();
    cbell->Update();
    cbell->Draw();
    
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
    
    TCanvas *cchg = new TCanvas("cchg","As a function of charge",0,0,600,600);
    cchg->Divide(2,1);
    cchg->cd(1);
    h3->Draw();
    cchg->cd(2);
    h4->Divide(h3);
    h4->Draw();
    cchg->Update();
    cchg->Draw();
    
    TCanvas *cchg2 = new TCanvas("cchg2","As a function of charge",0,0,600,600);
    cchg2->Divide(2,1);
    cchg2->cd(1);
    h3->Draw();
    cchg2->cd(2);
    h5->Divide(h3);
    h5->Draw();
    cchg2->Update();
    cchg2->Draw();
    
    // Look at the residuals wrt coordinate y of the hit, to see if there is a bias
    
    TCanvas *cerr = new TCanvas("cerr","Graph example",0,0,600,600);
    cerr->Divide(2,1);
       
    cerr->cd(1);
       TGraph *gr1 = new TGraph(count, ytrks, residualsMoi);
       gr1->SetTitle("Residuals wrt input y - My Clusetring");
       gr1->SetMarkerColor(4);
       gr1->SetLineColor(4);
       gr1->SetMarkerStyle(8);
       gr1->GetXaxis()->SetTitle("y input (cm)");
       gr1->GetYaxis()->SetTitle("Residual (cm)");
       gr1->Draw("AP");
    
    cerr->cd(2);
      TGraph *gr2 = new TGraph(count, ytrks, residualsAlberto);
      gr2->SetTitle("Residuals wrt input y - Alberto's Clustering");
      gr2->SetMarkerColor(4);
      gr2->SetLineColor(4);
      gr2->SetMarkerStyle(8);
      gr2->GetXaxis()->SetTitle("y input (cm)");
      gr2->GetYaxis()->SetTitle("Residual (cm)");
      gr2->Draw("AP");
    
    cerr->Update();
    cerr->Draw();
    
    TCanvas *cdiff = new TCanvas("cdiff","Graph example",0,0,600,600);
          
          TGraph *gr3 = new TGraph(count, ytrks, residualsDifference);
          gr3->SetTitle("Difference in residuals");
          gr3->SetMarkerColor(4);
          gr3->SetLineColor(4);
          gr3->SetMarkerStyle(8);
          gr3->GetXaxis()->SetTitle("y input (cm)");
          gr3->GetYaxis()->SetTitle("Residual difference (Mine-Alberto's) (cm)");
          gr3->Draw("AP");
       
       cdiff->Update();
       cdiff->Draw();
    
              app.Run(kTRUE);

  return 0;
}
