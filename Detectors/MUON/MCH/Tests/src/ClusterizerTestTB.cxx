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

#include "MCHBase/Digit.h"
#include "MCHBase/PreCluster.h"
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
  std::vector<float> residuals;
    
    TApplication app ("app",&argc,argv);
    
    TCanvas *cbell = new TCanvas("cbell","Bell",0,0,600,600);
    TH1F *h1 = new TH1F("h1", "Residuals distribution from TB data", 50, -0.1, 0.1);
  
  preClusterFinder.init();

  Digit* digitsBuffer = NULL;
  std::vector<Digit> digits(0);
  std::vector<PreCluster> preClusters(0);
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

    std::vector<TBCluster>& tbclusters = digitsReader.getClusters();
    std::cout<<"tbclusters.size(): "<<tbclusters.size()<<std::endl;
    for(size_t ic = 0; ic < tbclusters.size(); ic++) {
      std::cout<<"TBCluster "<<ic<<":\n"<<tbclusters[ic];
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
      clustering.runFinderSimpleFit(preClusters, digits, clusters);

      
      if(preClusters.size()==1){
          float yobtenu = clusters[0].gety();
          float difference = ytrk-yobtenu;
          h1->Fill(difference);
            h1->GetXaxis()->SetTitle("Residual y (cm)");
            h1->GetYaxis()->SetTitle("Count");
            h1->Draw();
          cout << "RESIDUAL y: " << difference <<endl;

      }

    //break;
  }
    
    cbell->Update();
    cbell->Draw();
              app.Run(kTRUE);

  return 0;
}
