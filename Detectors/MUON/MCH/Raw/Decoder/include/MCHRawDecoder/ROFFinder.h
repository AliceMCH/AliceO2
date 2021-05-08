// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file ROFFinder.h
/// \brief Class to group the fired pads according to their time stamp
///
/// \author Andrea Ferrero, CEA

#ifndef ALICEO2_MCH_ROFFINDER_H_
#define ALICEO2_MCH_ROFFINDER_H_

#include <cassert>
#include <cstdint>
#include <unordered_map>
#include <vector>

#include <gsl/span>

#include "DataFormatsMCH/Digit.h"
#include "MCHRawDecoder/DataDecoder.h"
#include "DataFormatsMCH/ROFRecord.h"

namespace o2
{
namespace mch
{
namespace raw
{

class ROFFinder
{
 public:
  static constexpr int maxObitsInTF = 256;

  struct RawDigitRef {
    const DataDecoder::RawDigit* ptr;

    RawDigitRef() : ptr(nullptr) {}
    RawDigitRef(const DataDecoder::RawDigit* p) : ptr(p) {}

    const DataDecoder::RawDigit& operator*() const
    {
      return (*ptr);
    }
    const DataDecoder::RawDigit* operator->() const
    {
      return (ptr);
    }
  };
  using RawDigit = DataDecoder::RawDigit;
  using RawDigitRef_ = const RawDigit*;
  using RawDigitSet = std::multiset<RawDigitRef>;
  using RawDigitArray = std::array<RawDigitSet, maxObitsInTF>;
  using RawDigitIter = RawDigitSet::const_iterator;
  using RawDigitIterArray = std::array<RawDigitIter, maxObitsInTF>;
  ROFFinder();
  ~ROFFinder();

  void setFirstTForbit(uint32_t orbit) { mFirstTForbit = orbit; }
  void process(const DataDecoder::RawDigitVector& digits, bool dummyROFs = false);
  char* storeDigits(size_t& bufSize);
  char* storeROFRecords(size_t& bufSize);

  o2::InteractionRecord DigitTime2IR(const RawDigit& digit);

  std::vector<RawDigitRef> getDigits() { return mOutputDigits; }
  std::vector<o2::mch::ROFRecord> getROFRecords() { return mOutputROFs; }

  bool isRofTimeMonotonic();
  bool isDigitsTimeAligned();

  void reset();

 private:
  void fillDigitsArray(const DataDecoder::RawDigitVector& digits);
  std::optional<int> getOrbitWithErlierDigit(const RawDigitArray& digits, int iOrbit);
  void addOneDigit(int iOrbit);
  void resetROF();
  void extendROF();
  void startNewROF();

  void dumpOrderedDigits(int iOrbit);
  void dumpOrderedDigits();
  void dumpOutputDigits();
  void dumpOutputROFs();

  int mFirstIdx = -1;
  int mEntries = 0;
  uint32_t mFirstTForbit;

  RawDigitIterArray mSeeds{}; ///< indexes of the digits used as seeds for starting new ROFs

  // references to decoded digits, sorted in time ascending order separately for each orbit;
  RawDigitArray mOrderedDigits;

  // vector of references to decoded digits, globally sorted in time ascending order
  std::vector<RawDigitRef> mOutputDigits{};
  // vector of ROFRecords, indexed with respect to the output digits vector
  std::vector<o2::mch::ROFRecord> mOutputROFs{};
};

} // namespace raw
} // namespace mch
} // namespace o2

#endif // ALICEO2_MCH_ROFFINDERSIMPLE_H_
