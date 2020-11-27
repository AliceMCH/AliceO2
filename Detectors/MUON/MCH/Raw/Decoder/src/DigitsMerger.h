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
/// \file    DatDecoder.cxx
/// \author  Andrea Ferrero
///
/// \brief Implementation of a data processor to run the raw decoding
///

#include <functional>
#include "MCHBase/Digit.h"

#define MCH_MERGER_FEEID_MAX 63

namespace o2
{
namespace mch
{
namespace raw
{

enum MergerDigitState
{
  DIGIT_STATE_UNCHECKED,
  DIGIT_STATE_MERGED,
  DIGIT_STATE_COMPLETED
};

//_________________________________________________________________
struct MergerBuffer {
  std::vector<std::pair<Digit,MergerDigitState>> digits;
  uint32_t orbit = {0};
  std::optional<Digit> lastDigit;
};

//_________________________________________________________________
class FeeIdMerger
{
 public:
  FeeIdMerger() = default;
  ~FeeIdMerger() = default;

  void setDebug(bool debug)
  {
    mDebug = debug;
  }

  void setId(int id)
  {
    feeId = id;
  }

  void setDigitHandler(std::function<void(const Digit&)> h)
  {
    sendDigit = h;
  }

  void setOrbit(uint32_t orbit, bool stop);

  MergerBuffer& getCurrentBuffer()
  {
    return currentBuffer;
  }

  MergerBuffer& getPreviousBuffer()
  {
    return previousBuffer;
  }

  void mergeDigits();

 private:
  MergerBuffer currentBuffer, previousBuffer;
  int feeId = {0};
  std::function<void(const Digit&)> sendDigit;
  bool mDebug = {false};
};

//_________________________________________________________________
class BaseMerger
{
 public:
  BaseMerger() = default;
  ~BaseMerger() = default;

  virtual void setDigitHandler(std::function<void(const Digit&)> h) = 0;
  virtual void setOrbit(int feeId, uint32_t orbit, bool stop) = 0;
  virtual void addDigit(int feeId, int solarId, int dsAddr, int chAddr,
                        int deId, int padId, unsigned long adc, Digit::Time time, uint16_t nSamples) = 0;
  virtual void mergeDigits(int feeId) = 0;
};

class NoOpMerger : public BaseMerger
{
 public:
  NoOpMerger(bool debug): mDebug(debug) { }
  ~NoOpMerger() = default;

  void setOrbit(int feeId, uint32_t orbit, bool stop)
  {
  }

  void setDigitHandler(std::function<void(const Digit&)> h)
  {
    sendDigit = h;
  }

  void addDigit(int feeId, int solarId, int dsAddr, int chAddr,
                int deId, int padId, unsigned long adc, Digit::Time time, uint16_t nSamples)
  {
    sendDigit(o2::mch::Digit(deId, padId, adc, time, nSamples));
  }

  void mergeDigits(int feeId) {}

 private:
  std::function<void(const Digit&)> sendDigit;
  bool mDebug = {false};
};

//_________________________________________________________________
class Merger : public BaseMerger
{
 public:
  Merger(bool debug);
  ~Merger() = default;

  void setOrbit(int feeId, uint32_t orbit, bool stop);

  void setDigitHandler(std::function<void(const Digit&)> h);

  void addDigit(int feeId, int solarId, int dsAddr, int chAddr,
                int deId, int padId, unsigned long adc, Digit::Time time, uint16_t nSamples);

  void mergeDigits(int feeId);

 private:
  int32_t feeId = {-1};
  FeeIdMerger mergers[MCH_MERGER_FEEID_MAX + 1];
  bool mDebug = {false};
};

} // namespace raw
} // namespace mch
} // end namespace o2
