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

//#include "MCHBase/Digit.h"
//#include "MCHBase/PreClusterBlock.h"
#include "TBDigitsFileReader.h"
#include "../../PreClustering/src/PreClusterFinder.h"

using namespace o2::mch;
using namespace std;

int main(int argc, char** argv)
{
    
  TBDigitsFileReader digitsReader(argv[1]);
  ofstream outFile(argv[2],ios::out);

  PreClusterFinder preClusterFinder;
  preClusterFinder.init();

  Digit* digitsBuffer = NULL;
    
  // load digits from binary input file, block-by-block
  while(digitsReader.readDigitsFromFile()) {

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

    outFile<<preClusterFinder<<std::endl;
  }

  return 0;
}
