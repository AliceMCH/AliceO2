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
//#include "MCHBase/DigitBlock.h"
#include "MCHPreClustering/PreClusterBlock.h"
#include "MCHPreClustering/PreClusterFinder.h"
#include "TBDigitsFileReader.h"

using namespace o2::mch;
using namespace std;

int main(int argc, char** argv)
{
  TBDigitsFileReader digitsReader(argv[1]);
  PreClusterFinder preClusterFinder;
  PreClusterBlock preClusterBlock;
  
  std::string fname;
  preClusterFinder.init(fname);

  Digit* digitsBuffer = NULL;
  char* preClustersBuffer = NULL;

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

    //break;
  }

  return 0;
}
