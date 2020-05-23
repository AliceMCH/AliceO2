// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include <memory>
#include <fstream>

#include "MCHBase/Digit.h"

using namespace o2::mch;
using namespace std;

namespace o2
{

namespace mch
{

class DigitsFileReader
{
 public:
  DigitsFileReader(std::string inputFileName);

  bool readDigitsFromFile();

  ssize_t getNumberOfDigits();
  void storeDigits(void* bufferPtr);

 private:
  std::ifstream mInputFile;
  std::vector<std::unique_ptr<Digit>> digits;
};

} // namespace mch
} // namespace o2
