// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include <iostream>

#include "DigitsFileReader.h"
#include "MCHMappingInterface/Segmentation.h"

using namespace o2::mch;
using namespace std;

struct DataBlockHeader {
  uint16_t fType;        // The type of the data block. Must contain a value defined by DataBlockType.
  uint16_t fRecordWidth; // The number of bytes each record uses.
  uint32_t fNrecords;    // Number of records in this data block.
};

struct DigitBlock {
  DataBlockHeader header; // Common data block header
};

struct DigitStruct {
  uint32_t uid;   // Digit ID in the current mapping (from OCDB)
  uint16_t index; // Digit index in the new mapping (produced internally)
  uint16_t adc;   // ADC value of signal
};

DigitsFileReader::DigitsFileReader(std::string inputFileName)
{
  mInputFile.open(inputFileName, ios::binary);
  if (!mInputFile.is_open()) {
    throw invalid_argument("Cannot open input file " + inputFileName);
  }
}

bool DigitsFileReader::readDigitsFromFile()
{
  const uint32_t SSizeOfDigitBlock = sizeof(DigitBlock);
  const uint32_t SSizeOfDigitStruct = sizeof(DigitStruct);

  digits.clear();

  /// send the digits of the current event

  DigitBlock digitBlock{};

  mInputFile.read(reinterpret_cast<char*>(&digitBlock), SSizeOfDigitBlock);
  if (mInputFile.fail()) {
    return false; // probably reached eof
  }

  if (digitBlock.header.fRecordWidth != SSizeOfDigitStruct) {
    throw length_error("incorrect size of digits. Content changed?");
  }

  // create the output message
  auto size = digitBlock.header.fNrecords * SSizeOfDigitStruct;
  //printf("fNrecords: %d  size: %d\n", (int)(digitBlock.header.fNrecords), (int)size);

  DigitStruct* bufferPtr = (DigitStruct*)malloc(size);
  if (!bufferPtr) {
    throw length_error("cannot allocate digits buffer");
  }

  // fill digits info
  if (size > 0) {
    mInputFile.read((char*)bufferPtr, size);
  } else {
    std::cout << "event is empty\n";
  }

  digits.reserve(digitBlock.header.fNrecords);
  for (unsigned int di = 0; di < digitBlock.header.fNrecords; di++) {
    const DigitStruct& digit = bufferPtr[di];
    uint32_t uid = digit.uid;
    uint16_t adc = digit.adc;
    uint32_t time = 0;

    uint16_t detId = uid & 0xFFF;
    uint16_t dualSampaId = (uid >> 12) & 0xFFF;
    uint16_t dualSampaChannel = (uid >> 24) & 0x3F;

    try {
      const mapping::Segmentation& segment = mapping::segmentation(detId);
      //mapping::Segmentation segment(detId);
      int padId = segment.findPadByFEE(dualSampaId, dualSampaChannel);
      if (padId < 0)
        continue;
      digits.push_back(std::make_unique<Digit>(time, detId, padId, adc));
      /*digits.push_back( std::make_unique<Digit>() );
      Digit* mchdigit = digits.back().get();
      mchdigit->setDetID(detId);
      mchdigit->setPadID(padId);
      mchdigit->setADC(adc);
      mchdigit->setTimeStamp(time);*/
    } catch (std::exception& e) {
      continue;
    }
  }

  return true;
}

ssize_t DigitsFileReader::getNumberOfDigits()
{
  return digits.size();
}

void DigitsFileReader::storeDigits(void* bufferPtr)
{
  Digit* ptr = (Digit*)bufferPtr;
  for (unsigned int di = 0; di < digits.size(); di++) {

    memcpy(ptr, digits[di].get(), sizeof(Digit));
    ptr += 1;
  }
}
