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

#include "MCHBase/Digit.h"
#include "MCHBase/PreCluster.h"
#include "DigitsFileReader.h"
#include "../../PreClustering/src/PreClusterFinder.h"
#include "MCHClustering/ClusteringForTest.h"

using namespace o2::mch;
using namespace std;

int main(int argc, char** argv)
{
  DigitsFileReader digitsReader(argv[1]);
  PreClusterFinder preClusterFinder;
  Clustering clustering;
  
  preClusterFinder.init();

  Digit* digitsBuffer = NULL;
  std::vector<Digit> digits(0);
  std::vector<PreCluster> preClusters(0);
  std::vector<Clustering::Cluster> clusters(0);

  // load digits from binary input file, block-by-block
  while(digitsReader.readDigitsFromFile()) {

    // get number of loaded digits and store them into a memory buffer
    auto nDigits = digitsReader.getNumberOfDigits();
    printf("nDigits: %d\n", (int)nDigits);
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
      
      //Uncomment the method you wish to use to clusterize
      
    //Runs the clustering of preClusters following a CenterOfGravity algorithm. Fills clusters.
//    clustering.runFinderCOG(preClusters, digits, clusters);
//    printf("Number of clusters obtained and saved: %lu\n", clusters.size());
      
      // Fit Mathieson
   clustering.runFinderSimpleFit(preClusters, digits, clusters);
      
      // Fit Simple Gaussienne
 //     clustering.runFinderGaussianFit(preClusters, digits, clusters);
      
      // Fit Double Gaussienne
//     clustering.runFinderDoubleGaussianFit(preClusters, digits, clusters);

    break;
  }

  return 0;
}
