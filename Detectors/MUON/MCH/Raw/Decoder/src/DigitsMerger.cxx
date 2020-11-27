// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

///
/// \file    DataDecoder.cxx
/// \author  Andrea Ferrero
///
/// \brief Implementation of a data processor to run the raw decoding
///

#include "DigitsMerger.h"

#include <iostream>
#include <fmt/format.h>

//static bool mPrint = false;

namespace o2
{
namespace mch
{
namespace raw
{


static std::ostream& operator<<(std::ostream& os, const o2::mch::Digit& d)
{
  auto tend = d.getTime().sampaTime + d.nofSamples() - 1;
  os << fmt::format("PAD ({:04d} {:04d})\tADC {:06d}  TIME ({} {} {:02d})  SIZE {}  END {}",
      d.getDetID(), d.getPadID(), d.getADC(), d.getTime().orbit, d.getTime().bunchCrossing,
      d.getTime().sampaTime, d.nofSamples(), tend)
            << ((tend >= 98) ? " *" : "");
  return os;
}

static int64_t digitsTimeDiff(const Digit& d1, const Digit& d2)
{
  const uint32_t bxCounterRollover = 0x100000;
  const uint32_t bxCounterHalfRollover = bxCounterRollover / 2;

  // compute time difference
  Digit::Time t1 = d1.getTime();
  Digit::Time t2 = d2.getTime();
  uint32_t bx1 = t1.bunchCrossing;
  uint32_t bx2 = t2.bunchCrossing;
  // correct for value rollover
  if ((bx2 + bxCounterHalfRollover) < bx1) {
    bx2 += bxCounterRollover;
  }

  int64_t t1full = bx1 + (t1.sampaTime << 2);
  int64_t t2full = bx2 + (t2.sampaTime << 2);
  int64_t timeDiff = t2full - t1full;

  return timeDiff;
}

static uint32_t digitsTimeGap(const Digit& d1, const Digit& d2)
{
  const uint32_t bxCounterRollover = 0x100000;

  // compute time difference
  Digit::Time startTime = d2.getTime();
  uint32_t bxStart = startTime.bunchCrossing;
  Digit::Time stopTime = d1.getTime();
  stopTime.sampaTime += d1.nofSamples() - 1;
  uint32_t bxStop = stopTime.bunchCrossing;
  // correct for value rollover
  if (bxStart < bxStop) {
    bxStart += bxCounterRollover;
  }

  uint32_t stopTimeFull = bxStop + (stopTime.sampaTime << 2);
  uint32_t startTimeFull = bxStart + (startTime.sampaTime << 2);
  uint32_t timeDiff = startTimeFull - stopTimeFull;

  return timeDiff;
}

static int64_t digitTime(const Digit& d)
{
  Digit::Time t1 = d.getTime();
  uint32_t bx1 = t1.bunchCrossing;
  int64_t t1full = bx1 + (t1.sampaTime << 2);

  return t1full;
}

void FeeIdMerger::setOrbit(uint32_t orbit, bool stop)
{
  // perform the merging and send digits of previous orbit if either the stop RDH is received
  // or a new orbit is started
  // the code is written in a way that is robust against the presence of a stop RDH. if missing, the
  // merging gets triggered by a change in the orbit number, otherwise the stopRDH will trigger the
  // merging earlier
  if ((orbit == currentBuffer.orbit) && (!stop)) {
    return;
  }
  if (mDebug) {
    std::cout << "[FeeIdMerger::setOrbit] changing orbit from " << currentBuffer.orbit << " to " << orbit
        << " (previous is " << previousBuffer.orbit << ")" << std::endl;
  }

  // do the merging
  mergeDigits();

  // send the digits that are not merged into others
  int nSent = 0;
  int64_t digitTimeLast = 0;
  Digit* lastDigit = nullptr;
  for (auto& d : previousBuffer.digits) {
    int64_t t = digitTime(d.first);
    if (t > digitTimeLast) {
      digitTimeLast = t;
      lastDigit = &(d.first);
    }
    if ((d.second == DIGIT_STATE_UNCHECKED) && (d.first.getPadID() >= 0)) {
      sendDigit(d.first);
      nSent += 1;
    }
  }
  if (mDebug) {
    std::cout << "[FeeIdMerger] sent " << nSent << " digits for orbit " << previousBuffer.orbit << "  current orbit is " << orbit << std::endl;
  }

  if (lastDigit) {
    previousBuffer.lastDigit = *lastDigit;
  }

  if (mDebug) {
    std::cout << "[FeeIdMerger::setOrbit] last digit of orbit " << previousBuffer.orbit << ": ";
    if (previousBuffer.lastDigit.has_value()) {
      std::cout << previousBuffer.lastDigit.value();
    } else {
      std::cout << "null";
    }
    std::cout << std::endl;

    std::cout<<"[FeeIdMerger::setOrbit] number of digits for orbit " << currentBuffer.orbit
        << ": " << currentBuffer.digits.size() << std::endl;
  }

  if (previousBuffer.lastDigit.has_value()) {
    const uint32_t cBxDiffMax = 1000;
    for (auto& d : currentBuffer.digits) {
      if ((d.second == DIGIT_STATE_UNCHECKED) && (d.first.getPadID() >= 0)) {
        int64_t tdiff = digitsTimeDiff(previousBuffer.lastDigit.value(), d.first);
        if (mDebug) {
          std::cout << "  checking digit " << d.first << " - tdiff: " << tdiff <<std::endl;
        }
        if (tdiff <= cBxDiffMax) {
          d.second = DIGIT_STATE_COMPLETED;
          sendDigit(d.first);
          nSent += 1;
        }
      }
    }
  }

  if (stop) {
    // cleat the previous buffer if a stop RDH is received
    // this way, the following change in orbit number will not trigger a new merging
    previousBuffer.digits.clear();
  }

  if (orbit != currentBuffer.orbit) {
    // clear the contents of the buffer from the previous orbit, and swap the vectors
    previousBuffer.digits.clear();
    previousBuffer.lastDigit.reset();
    previousBuffer.orbit = currentBuffer.orbit;
    std::swap(previousBuffer.digits, currentBuffer.digits);
    currentBuffer.orbit = orbit;
  }
}

// helper function to check if two digits correspond to the same pad;
static bool areSamePad(const Digit& d1, const Digit& d2)
{
  if (d1.getDetID() != d2.getDetID())
    return false;
  if (d1.getPadID() != d2.getPadID())
    return false;

  return true;
}


static void mergeTwoDigits(Digit& d1, const Digit& d2, bool debug)
{
  if (debug) {
    std::cout << "Merging digits\n  " << d1 << std::endl << "and\n  " << d2 <<std::endl;
  }
  d1.setADC(d1.getADC() + d2.getADC());
  d1.setNofSamples(d1.nofSamples() + d2.nofSamples());
}


static bool updateDigits(Digit& d1, Digit& d2, bool debug)
{
  const uint32_t oneADCclockCycle = 4;

  // skip all digits not matching the detId/padId
  if (!areSamePad(d1, d2)) {
    return false;
  }

  // compute time difference
  uint32_t timeGap = digitsTimeGap(d1, d2);

  // skip if the time difference is not equal to 1 ADC clock cycle
  if (timeGap != oneADCclockCycle) {
    return false;
  }

  // merge digits
  mergeTwoDigits(d1, d2, debug);

  return true;
}


static void mergeBuffers(MergerBuffer& buf1, MergerBuffer& buf2, bool debug)
{
  for (size_t i = 0; i < buf2.digits.size(); i++) {
    auto& d2 = buf2.digits[i];

    // skip already merged digits
    if (d2.second != DIGIT_STATE_UNCHECKED) {
      continue;
    }

    // skip digits that do not start at the beginning of the time window
    Digit::Time startTime = d2.first.getTime();
    if (startTime.sampaTime != 0) {
      continue;
    }

    for (size_t j = 0; j < buf1.digits.size(); j++) {
      auto& d1 = buf1.digits[j];

      // skip already merged digits
      if (d1.second != d2.second != DIGIT_STATE_UNCHECKED) {
        continue;
      }

      if (updateDigits(d1.first, d2.first, debug)) {
        // mark d2 as merged
        d2.second = DIGIT_STATE_MERGED;
        break;
      }
    }
  }
}


void FeeIdMerger::mergeDigits()
{
  auto& currentBuffer = getCurrentBuffer();
  auto& previousBuffer = getPreviousBuffer();

  uint32_t orbit1 = previousBuffer.orbit;
  uint32_t orbit2 = currentBuffer.orbit;
  if (mDebug) {
    std::cout << "[FeeIdMerger::mergeDigits] merging orbit " << orbit2 << std::endl;
  }
  mergeBuffers(currentBuffer, currentBuffer, mDebug);

  // only merge digits from consecutive orbits
  if ((orbit2 >= orbit1) && ((orbit2 - orbit1) > 1)) {
    return;
  }

  if (mDebug) {
    std::cout << "[FeeIdMerger::mergeDigits] merging orbits " << orbit1 << " and " << orbit2 <<std::endl;
  }
  mergeBuffers(previousBuffer, currentBuffer, mDebug);
}

Merger::Merger(bool debug): mDebug(debug)
{
  for (int i = 0; i <= MCH_MERGER_FEEID_MAX; i++) {
    mergers[i].setDebug(mDebug);
  }
}

void Merger::setOrbit(int feeId, uint32_t orbit, bool stop)
{
  if (feeId < 0 || feeId > MCH_MERGER_FEEID_MAX) {
    return;
  }

  mergers[feeId].setOrbit(orbit, stop);
}

void Merger::setDigitHandler(std::function<void(const Digit&)> h)
{
  for (int feeId = 0; feeId <= MCH_MERGER_FEEID_MAX; feeId++) {
    mergers[feeId].setDigitHandler(h);
  }
}

void Merger::addDigit(int feeId, int solarId, int dsAddr, int chAddr,
                      int deId, int padId, unsigned long adc, Digit::Time time, uint16_t nSamples)
{
  if (feeId < 0 || feeId > MCH_MERGER_FEEID_MAX) {
    return;
  }

  mergers[feeId].getCurrentBuffer().digits.emplace_back(std::make_pair(o2::mch::Digit{deId, padId, adc, time, nSamples}, DIGIT_STATE_UNCHECKED));
}

void Merger::mergeDigits(int feeId)
{
  if (feeId < 0 || feeId > MCH_MERGER_FEEID_MAX) {
    return;
  }

  mergers[feeId].mergeDigits();
}

} // namespace raw
} // namespace mch
} // end namespace o2
