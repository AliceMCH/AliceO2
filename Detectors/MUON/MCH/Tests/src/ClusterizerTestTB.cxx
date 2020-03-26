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
//#include "MCHBase/DigitBlock.h"
#include "MCHPreClustering/PreClusterBlock.h"
#include "MCHPreClustering/PreClusterFinder.h"
#include "TBDigitsFileReader.h"
#include "MCHClustering/ClusteringForTest.h"

using namespace o2::mch;
using namespace std;

int main(int argc, char** argv)
{
    
  TBDigitsFileReader digitsReader(argv[1]);
  PreClusterFinder preClusterFinder;
  PreClusterBlock preClusterBlock;
  Clustering clustering;
  float residuals[1000];
    float ytrks[1000];
    int count = 0;
    
    TApplication app ("app",&argc,argv);
    
    TCanvas *cbell = new TCanvas("cbell","Bell",0,0,600,600);
    TH1F *h1 = new TH1F("h1", "Residuals distribution from TB data", 50, -0.1, 0.1);
  
  std::string fname;
  preClusterFinder.init(fname);

  Digit* digitsBuffer = NULL;
  char* preClustersBuffer = NULL;
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
    //continue;
    digitsBuffer = (Digit*)realloc(digitsBuffer, sizeof(Digit) * nDigits);
    digitsReader.storeDigits(digitsBuffer);

    // load the digits from the memory buffer and run the pre-clustering phase
    preClusterFinder.reset();
    preClusterFinder.loadDigits(digitsBuffer, nDigits);
    preClusterFinder.run();

    // get number of pre-clusters and store them into a memory buffer
    auto preClustersSize = preClusterBlock.getPreClustersBufferSize(preClusterFinder);
    printf("preClustersSize: %d\n", (int)preClustersSize);
    preClustersBuffer = (char*)realloc(preClustersBuffer, preClustersSize);
    preClusterBlock.storePreClusters(preClusterFinder, preClustersBuffer);

    //continue;
    printf("\n\n==========\nReading clusters\n\n");

    std::vector<PreClusterStruct> preClusters;
    preClusterBlock.readPreClusters(preClusters, preClustersBuffer, preClustersSize);
      
      printf("\n\n==========\nRunning Clustering\n\n");

      
      
       //Uncomment the method you wish to use to clusterize
            
          //Runs the clustering of preClusters following a CenterOfGravity algorithm. Fills clusters.
//          clustering.runFinderCOG(preClusters, clusters);
//          printf("Number of clusters obtained and saved: %lu\n", clusters.size());
            
            // Fit Mathieson
         clustering.runFinderSimpleFit(preClusters, clusters);
            
            // Fit Simple Gaussienne
       //     clustering.runFinderGaussianFit(preClusters, clusters);
            
            // Fit Double Gaussienne
      //     clustering.runFinderDoubleGaussianFit(preClusters, clusters);

      
      
      
      
      // Choose to only look at nice events, where we have one precluster reconstructed as a single cluster, which are a majority (972 out of roughly 1000). This helps us get rid of pathological events.
      
      if(preClusters.size()==1){
          float yobtenu = clusters[0].gety();
          float difference = ytrk-yobtenu;
          ytrks[count] = ytrk;
          residuals[count] = difference;
          
          // Look at the distribution of residuals from this TestBeam data
          
          h1->Fill(difference);
            h1->GetXaxis()->SetTitle("Residual y (cm)");
            h1->GetYaxis()->SetTitle("Count");
            h1->Draw();
          cout << "RESIDUAL y: " << difference <<endl;
          count++;
      }
      
    //break;
  }
    
    cbell->Update();
    cbell->Draw();
    
    // Look at the residuals wrt coordinate y of the hit, to see if there is a bias
    
    TCanvas *cerr = new TCanvas("cerr","Graph example",0,0,600,600);
       
       TGraph *gr1 = new TGraph(count, ytrks, residuals);
       gr1->SetTitle("Residuals wrt input y");
       gr1->SetMarkerColor(4);
       gr1->SetLineColor(4);
       gr1->SetMarkerStyle(8);
       gr1->GetXaxis()->SetTitle("y input (cm)");
       gr1->GetYaxis()->SetTitle("Residual (cm)");
       gr1->Draw("AP");
    
    cerr->Update();
    cerr->Draw();
              app.Run(kTRUE);

  return 0;
}
