// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include "MCHRawDecoder/ROFFinder.h"

#include <chrono>
#include <memory>
#include <stdexcept>
#include <vector>
#include <map>
#include <algorithm>
#include <fmt/format.h>

#include <fairmq/Tools.h>
#include <FairMQLogger.h>

#include "MCHMappingInterface/Segmentation.h"

//#define ROFDEBUG 1

namespace o2
{
namespace mch
{
namespace raw
{

using namespace std;

bool operator<(const ROFFinder::RawDigitRef& d1, const ROFFinder::RawDigitRef& d2)
{
  if ((*d1) < (*d2))
    return true;
  return false;
}

bool operator<(ROFFinder::RawDigitIter& d1, ROFFinder::RawDigitIter& d2)
{
  if ((*(*d1)) < (*(*d2)))
    return true;
  return false;
}

//_________________________________________________________________________________________________
ROFFinder::ROFFinder()
{
}

//_________________________________________________________________________________________________
ROFFinder::~ROFFinder() = default;

//_________________________________________________________________________________________________
void ROFFinder::process(const DataDecoder::RawDigitVector& digits, bool dummyROFs)
{
  if (dummyROFs) {
    for (size_t i = 0; i < digits.size(); i++) {
      RawDigitRef digit{&(digits[i])};
      mOutputDigits.push_back(digit);
    }

    if (!digits.empty()) {
      mOutputROFs.emplace_back(DigitTime2IR(digits[0]), 0, mOutputDigits.size());
    }
    return;
  }

  fillDigitsArray(digits);

  resetROF();

  // initialize the references to the first digit in each orbit, which wil be used as seeds for the ROF search
  for (int iOrbit = 0; iOrbit < maxObitsInTF; iOrbit++) {
    mSeeds[iOrbit] = mOrderedDigits[iOrbit].begin();
  }

  for (int i = 0; i < maxObitsInTF; i++) {
    // loop over digits until the current orbit is completely processed
    while (std::optional<int> iOrbit = getOrbitWithErlierDigit(mOrderedDigits, i)) {
      addOneDigit(iOrbit.value());
    }
  }

  // close the last ROF
  startNewROF();
}

//_________________________________________________________________________________________________
void ROFFinder::fillDigitsArray(const DataDecoder::RawDigitVector& digits)
{
  for (size_t i = 0; i < digits.size(); i++) {
    RawDigitRef digit{&(digits[i])};
    if (!digit->timeValid()) {
      LOG(ERROR) << "Digit with invalid time, DS " << digit->info.solar << "," << digit->info.ds << "," << digit->info.chip
                 << "  pad " << digit->getDetID() << "," << digit->getPadID() << "  "
                 << digit->getOrbit() << " " << digit->getTime() << " " << digit->getBXTime();
      continue;
    }

    auto orbit = digit->getOrbit();
    if (orbit < mFirstTForbit) {
      LOG(ERROR) << "[ROFFinder::fillDigitsArray] orbit smaller than first TF orbit: " << orbit << ", " << mFirstTForbit;
      continue;
    }

    int iOrbit = orbit - mFirstTForbit;
    if (iOrbit >= maxObitsInTF) {
      LOG(ERROR) << "[ROFFinder::fillDigitsArray] orbit larger than max orbits in TF: " << orbit << ", " << mFirstTForbit;
      continue;
    }

#ifdef ROFDEBUG
    std::cout << "Inserting digit into orbit " << iOrbit << std::endl;
    std::cout << " pad " << digit->getDetID() << "," << digit->getPadID() << " "
              << digit->getOrbit() << " " << digit->getTime() << " " << digit->getBXTime() << std::endl;
#endif

    mOrderedDigits[iOrbit].insert(digit);
  }

#ifdef ROFDEBUG
  std::cout << "ORDERED DIGITS:\n";
  dumpOrderedDigits();
#endif

  size_t orderedSize = 0;
  for (int i = 0; i < maxObitsInTF; i++) {
    orderedSize += mOrderedDigits[i].size();
  }
  if (orderedSize != digits.size()) {
    LOG(ERROR) << "[ROFFinder::fillDigitsArray] digits array size mismatch: " << orderedSize << " instead of " << digits.size();
  }
}

//_________________________________________________________________________________________________
std::optional<int> ROFFinder::getOrbitWithErlierDigit(const RawDigitArray& digits, int iOrbit)
{
  auto iterIsValid = [&](RawDigitIter iter, const RawDigitSet& digits) {
    return (iter != digits.end());
  };

  // initialize references to iterators and digit sets of current orbit
  auto& iter1 = mSeeds[iOrbit];
  auto& set1 = digits[iOrbit];
  // if we reached the end of the current orbit, return a null result. This will stop the loop over the digits
  if (!iterIsValid(iter1, set1)) {
    return std::nullopt;
  }

  int iOrbitEarlier = iOrbit;

  // check the first unprocessed digit of the next orbit, unless we are at the last one in the TF
  if (iOrbit != (digits.size() - 1)) {
    auto& iter2 = mSeeds[iOrbit + 1];
    auto& set2 = digits[iOrbit + 1];
    if (iterIsValid(iter2, set2)) {
#ifdef ROFDEBUG
      std::cout << " checking if " << iOrbit + 1 << " is earlier than " << iOrbit << std::endl;
      std::cout << "  digit 1: "
                << " DE# " << (*iter1)->getDetID() << " PadId " << (*iter1)->getPadID() << "  ADC " << (*iter1)->getADC()
                << " time " << (*iter1)->getTime() << "," << (*iter1)->getSampaTime() << std::endl;
      std::cout << "  digit 2: "
                << " DE# " << (*iter2)->getDetID() << " PadId " << (*iter2)->getPadID() << "  ADC " << (*iter2)->getADC()
                << " time " << (*iter2)->getTime() << "," << (*iter2)->getSampaTime() << std::endl;
#endif
      if (iter2 < iter1) {
        iOrbitEarlier += 1;
#ifdef ROFDEBUG
        std::cout << "  YES!!!\n";
      } else {
        std::cout << "  no\n";
#endif
      }
    }
  }

  return std::optional<int>(iOrbitEarlier);
}

//_________________________________________________________________________________________________
void ROFFinder::addOneDigit(int iOrbit)
{
  auto& digitIter = mSeeds[iOrbit];
  auto digitPtr = *digitIter;
  mOutputDigits.push_back(digitPtr);
  ++digitIter;
#ifdef ROFDEBUG
  std::cout << "DIGIT ADDED TO ROF:\n ORBIT " << iOrbit
            << " DE# " << digitPtr->getDetID() << " PadId " << digitPtr->getPadID() << "  ADC " << digitPtr->getADC()
            << " time " << digitPtr->getTime() << "," << digitPtr->getSampaTime() << std::endl;
#endif

  // if no ROF has been created yet, initialize the index to the first digit and set the number of entries to 1
  // otherwise compute the dime difference with respect to the first digit in the ROF. If it smaller than 4 bunch crossings
  // add the digit to the current ROF, otherwose close the ROF and open a new one
  if (mEntries == 0) {
    startNewROF();
    return;
  }

  constexpr int oneADCclock = 4;
  const auto d1 = mOutputDigits[mFirstIdx];
  const auto d2 = digitPtr;
  auto tdiff = d2->getTime() - d1->getTime();
  if (tdiff < oneADCclock) {
    extendROF();
  } else {
    startNewROF();
  }
}

//_________________________________________________________________________________________________
void ROFFinder::reset()
{
  resetROF();
  mOutputDigits.clear();
  mOutputROFs.clear();
  for (int i = 0; i < maxObitsInTF; i++) {
    mOrderedDigits[i].clear();
  }
}

//_________________________________________________________________________________________________
void ROFFinder::resetROF()
{
  mFirstIdx = -1;
  mEntries = 0;
}

//_________________________________________________________________________________________________
void ROFFinder::extendROF()
{
  mEntries += 1;
}

//_________________________________________________________________________________________________
o2::InteractionRecord ROFFinder::DigitTime2IR(const RawDigit& digit)
{
  constexpr int bcInOneOrbit = 3564;
  uint32_t orbit = digit.getTime() / bcInOneOrbit + mFirstTForbit;
  int32_t bc = digit.getTime() % bcInOneOrbit;
  return o2::InteractionRecord(bc, orbit);
}

//_________________________________________________________________________________________________
void ROFFinder::startNewROF()
{
  constexpr int bcInOneOrbit = 3564;
  if (mEntries > 0) {
    const auto digit = mOutputDigits[mFirstIdx];
    mOutputROFs.emplace_back(DigitTime2IR(*digit), mFirstIdx, mEntries);
  }
  mFirstIdx = mOutputDigits.size() - 1;
  mEntries = 1;
}

//_________________________________________________________________________________________________
char* ROFFinder::storeDigits(size_t& bufSize)
{
  static constexpr size_t sizeOfDigit = sizeof(o2::mch::Digit);

#ifdef ROFDEBUG
  dumpOutputDigits();
#endif

  bufSize = mOutputDigits.size() * sizeOfDigit;
  o2::mch::Digit* buf = reinterpret_cast<o2::mch::Digit*>(malloc(bufSize));
  if (!buf) {
    bufSize = 0;
    return nullptr;
  }

  o2::mch::Digit* p = buf;
  for (size_t i = 0; i < mOutputDigits.size(); i++) {
    auto d = mOutputDigits[i];
    memcpy(p, &(d->digit), sizeOfDigit);
    p += 1;
  }

  return reinterpret_cast<char*>(buf);
}

//_________________________________________________________________________________________________
char* ROFFinder::storeROFRecords(size_t& bufSize)
{
  static constexpr size_t sizeOfROFRecord = sizeof(o2::mch::ROFRecord);

#ifdef ROFDEBUG
  dumpOutputROFs();
#endif

  bufSize = mOutputROFs.size() * sizeOfROFRecord;
  o2::mch::ROFRecord* buf = reinterpret_cast<o2::mch::ROFRecord*>(malloc(bufSize));
  if (!buf) {
    bufSize = 0;
    return nullptr;
  }

  o2::mch::ROFRecord* p = buf;
  for (size_t i = 0; i < mOutputROFs.size(); i++) {
    auto& rof = mOutputROFs[i];
    memcpy(p, &(rof), sizeOfROFRecord);
    p += 1;
  }

  return reinterpret_cast<char*>(buf);
}

//_________________________________________________________________________________________________
void ROFFinder::dumpOrderedDigits(int iOrbit)
{
  for (auto& digit : mOrderedDigits[iOrbit]) {
    auto& d = digit->digit;
    auto& t = digit->info;

    if (d.getPadID() < 0) {
      continue;
    }
    const o2::mch::mapping::Segmentation& segment = o2::mch::mapping::segmentation(d.getDetID());
    bool bending = segment.isBendingPad(d.getPadID());
    float X = segment.padPositionX(d.getPadID());
    float Y = segment.padPositionY(d.getPadID());
    uint32_t orbit = t.orbit;
    uint32_t bunchCrossing = t.bunchCrossing;
    uint32_t sampaTime = t.sampaTime;
    std::cout << fmt::format("    DE {:4d}  PAD {:5d}  ADC {:6d}  TIME {} ({} {} {:4d} {})",
                             d.getDetID(), d.getPadID(), d.getADC(), d.getTime(), orbit, bunchCrossing, sampaTime, t.getBXTime());
    std::cout << fmt::format("\tC {}  PAD_XY {:+2.2f} , {:+2.2f}", (bending ? (int)0 : (int)1), X, Y);
    std::cout << std::endl;
  }
}

void ROFFinder::dumpOrderedDigits()
{
  for (int i = 0; i < maxObitsInTF; i++) {
    if (mOrderedDigits[i].empty()) {
      continue;
    }
    std::cout << "  ORBIT " << i << std::endl;
    dumpOrderedDigits(i);
  }
}

//_________________________________________________________________________________________________
void ROFFinder::dumpOutputDigits()
{
  std::cout << "OUTPUT DIGITS:\n";
  for (size_t i = 0; i < mOutputDigits.size(); i++) {
    auto d = mOutputDigits[i];

    int iROF = -1;
    for (size_t j = 0; j < mOutputROFs.size(); j++) {
      const auto& rof = mOutputROFs[j];
      if (rof.getFirstIdx() <= i && rof.getLastIdx() >= i) {
        iROF = j;
      }
    }
    std::cout << fmt::format("    DIGIT {}  ROF {}  DE {:4d}  PAD {:5d}  ADC {:6d}  TIME {} ({} {} {:4d} {})\n",
                             i, iROF, d->getDetID(), d->getPadID(), d->getADC(),
                             d->getTime(), d->getOrbit(), d->getBunchCrossing(), d->getSampaTime(), d->getBXTime());
  }
}

//_________________________________________________________________________________________________
void ROFFinder::dumpOutputROFs()
{
  std::cout << "OUTPUT ROFs:\n";
  for (size_t i = 0; i < mOutputROFs.size(); i++) {
    auto& rof = mOutputROFs[i];
    std::cout << fmt::format("    ROF {} {}-{} {},{}\n",
                             i, rof.getFirstIdx(), rof.getLastIdx(), rof.getBCData().orbit, rof.getBCData().bc);
  }
}

//_________________________________________________________________________________________________
bool ROFFinder::isRofTimeMonotonic()
{
  for (size_t i = 1; i < mOutputROFs.size(); i++) {
    const auto& rof = mOutputROFs[i];
    const auto& rofPrev = mOutputROFs[i - 1];
    if (rof.getBCData() < rofPrev.getBCData()) {
      return false;
    }
  }
  return true;
}

//_________________________________________________________________________________________________
bool ROFFinder::isDigitsTimeAligned()
{
  for (size_t i = 0; i < mOutputROFs.size(); i++) {
    auto& rof = mOutputROFs[i];
    for (int di = rof.getFirstIdx() + 1; di <= rof.getLastIdx(); di++) {
      const auto& digit = *(mOutputDigits[di]);
      const auto& digitPrev = *(mOutputDigits[di - 1]);
      if (digit.getTime() != digitPrev.getTime()) {
        return false;
      }
    }
  }
  return true;
}

} // namespace raw
} // namespace mch
} // namespace o2
