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
/// \file    DataDecoder.h
/// \author  Andrea Ferrero
///
/// \brief Definition of the decoder for the MCH data
///

#ifndef O2_MCH_DATADECODER_H_
#define O2_MCH_DATADECODER_H_

#include <gsl/span>
#include <unordered_set>
#include <unordered_map>
#include <queue>
#include <mutex>
#include <condition_variable>
#include <thread>

#include "Headers/RDHAny.h"
#include "DataFormatsMCH/Digit.h"
#include "MCHRawDecoder/OrbitInfo.h"
#include "MCHRawDecoder/PageDecoder.h"

namespace o2
{
namespace mch
{
namespace raw
{

enum DigitsMappingMode
{
  eDigitsMappingStandard = 0,
  eDigitsMappingFast1 = 1,
  eDigitsMappingFast2 = 2
};

using RdhHandler = std::function<void(o2::header::RDHAny*)>;

// custom hash for OrbitInfo objects
struct OrbitInfoHash {
  std::size_t operator()(const OrbitInfo& info) const noexcept
  {
    return std::hash<uint64_t>{}(info.get());
  }
};

void dumpOrbits(const std::unordered_set<OrbitInfo, OrbitInfoHash>& mOrbits);


template<class T>
class SafeQueue
{
public:
  SafeQueue() = default;

  void stop()
  {
    std::lock_guard<std::mutex> lock(mQueueMutex);
    mStop = true;
    mQueueCondition.notify_one();
  }

  void push(const T val)
  {
    std::lock_guard<std::mutex> lock(mQueueMutex);
    mQueue.push(val);
    mQueueCondition.notify_one();
  }

  bool get(T val)
  {
    std::unique_lock<std::mutex> lock(mQueueMutex);
    mQueueCondition.wait(lock, [&]{ return (!mQueue.empty() || mStop); });

    if (mStop) { return false; }

    val = mQueue.front();
    return true;
  }

  void pop(ssize_t& remaining)
  {
    std::lock_guard<std::mutex> lock(mQueueMutex);
    mQueue.pop();
    remaining = mQueue.size();
  }

  size_t size() {
    std::lock_guard<std::mutex> lock(mQueueMutex);
    return mQueue.size();
  }

  bool empty() {
    std::lock_guard<std::mutex> lock(mQueueMutex);
    return mQueue.empty();
  }

private:
  std::mutex mQueueMutex;
  std::condition_variable mQueueCondition;
  std::queue<T> mQueue;
  bool mStop{false};
};

//_________________________________________________________________
//
// Data decoder
//_________________________________________________________________
class DataDecoder
{
 public:
  static constexpr int32_t tfTimeMax{0x7FFFFFFF};
  static constexpr int32_t tfTimeInvalid{-tfTimeMax};

  struct SampaInfo {
    union {
      uint32_t id = 0;
      struct {
        uint32_t chip : 1;
        uint32_t ds : 6;
        uint32_t solar : 16;
        uint32_t unused : 9;
      };
    };

    union {
      // default value
      uint64_t time = 0x0000000000000000;
      struct {                       ///
        uint32_t sampaTime : 10;     /// bit 0 to 9: sampa time
        uint32_t bunchCrossing : 20; /// bit 10 to 29: bunch crossing counter
        uint32_t reserved : 2;       /// bit 30 to 31: reserved
        uint32_t orbit;              /// bit 32 to 63: orbit
      };                             ///
    };
    uint32_t getBXTime() const
    {
      return (bunchCrossing + (sampaTime * 4));
    }
    int32_t tfTime;
    bool timeValid() const
    {
      return (tfTime != tfTimeInvalid);
    }

    bool operator==(const SampaInfo&) const;
    bool operator<(const SampaInfo& rhs) const
    {
      if (id < rhs.id) {
        return true;
      } else if (time < rhs.time) {
        return true;
      }
      return false;
    }
  };

  struct SampaTimeFrameStart {
    SampaTimeFrameStart() = default;
    SampaTimeFrameStart(uint32_t orbit, uint32_t bunchCrossing) : mOrbit(orbit), mBunchCrossing(bunchCrossing) {}

    uint32_t mOrbit{0};
    uint32_t mBunchCrossing{0};
  };

  struct RawDigit {
    o2::mch::Digit digit;
    SampaInfo info;
    auto getDetID() const { return digit.getDetID(); }
    auto getPadID() const { return digit.getPadID(); }
    uint32_t getADC() const { return digit.getADC(); }
    auto getTime() const { return info.tfTime; }
    bool timeValid() const { return info.timeValid(); }
    auto getOrbit() const { return info.orbit; }
    auto getBunchCrossing() const { return info.bunchCrossing; }
    auto getSampaTime() const { return info.sampaTime; }
    auto getBXTime() const { return info.getBXTime(); }

    bool operator==(const RawDigit&) const;
  };

  using RawDigitVector = std::vector<RawDigit>;

  DataDecoder(SampaChannelHandler channelHandler, RdhHandler rdhHandler,
              uint32_t sampaBcOffset,
              std::string mapCRUfile, std::string mapFECfile,
              bool ds2manu, bool verbose, bool useDummyElecMap,
              DigitsMappingMode mappingMode);

  void reset();
  void decodeBuffer(gsl::span<const std::byte> buf);

  void setFirstOrbitInRun(uint32_t orbit) { mFirstOrbitInRun = orbit; }
  std::optional<uint32_t> getFirstOrbitInRun() { return mFirstOrbitInRun; }
  void setFirstOrbitInTF(uint32_t orbit);

  void setSampaBcOffset(uint32_t offset) { mSampaTimeOffset = offset; }
  uint32_t getSampaBcOffset() const { return mSampaTimeOffset; }

  static int32_t digitsTimeDiff(uint32_t orbit1, uint32_t bc1, uint32_t orbit2, uint32_t bc2);
  static void computeDigitsTime(RawDigitVector& digits, SampaTimeFrameStart& sampaTimeFrameStart, bool debug);
  void computeDigitsTime()
  {
    computeDigitsTime(mDigits, mSampaTimeFrameStart, mDebug);
  }

  const RawDigitVector& getDigits() const { return mDigits; }
  const std::unordered_set<OrbitInfo, OrbitInfoHash>& getOrbits() const { return mOrbits; }

 private:
  void initElec2DetMapper(std::string filename);
  void initFee2SolarMapper(std::string filename);
  void init();
  void decodePage(gsl::span<const std::byte> page);
  void dumpDigits();
  bool getPadMapping(const DsElecId& dsElecId, DualSampaChannelId channel, int& deId, int& dsIddet, int& padId);
  bool addDigit(const DsElecId& dsElecId, DualSampaChannelId channel, const o2::mch::raw::SampaCluster& sc);

  Elec2DetMapper mElec2Det{nullptr};       ///< front-end electronics mapping
  FeeLink2SolarMapper mFee2Solar{nullptr}; ///< CRU electronics mapping
  std::string mMapFECfile;                 ///< optional text file with custom front-end electronics mapping
  std::string mMapCRUfile;                 ///< optional text file with custom CRU mapping

  struct MapDS
  {
    MapDS()
    {
      std::fill(padIds, padIds + 64, -1);
    }

    int deId{-1};
    int dsIddet{-1};
    int padIds[64];
  };
  using MapSolar = std::array<MapDS, 40>;
  using MapFee = std::array<MapSolar, 12>;
  std::vector<MapFee> mMapMCH{64};
  std::vector<MapSolar> mMapSolar{1024};
  DigitsMappingMode mMappingMode;

  struct Condition
  {
    std::mutex mMutex;
    std::condition_variable mCondition;
  };

  struct PageProcessor
  {
    PageProcessor(Elec2DetMapper elec2Det, std::vector<MapSolar>& mapSolar, RawDigitVector& digits, std::mutex& digitsMutex, Condition& processEnd, bool debug):
      mElec2Det(elec2Det), mMapSolar(mapSolar), mDigits(digits), mDigitsMutex(digitsMutex), mProcessEndCondition(processEnd), mDebug(debug)
    {

    }

    void spawn()
    {
      mThread = std::thread(&PageProcessor::run, this);
    }

    void run();
    void decodePage(gsl::span<const std::byte> page);
    void appendPage(gsl::span<const std::byte> page)
    {
      mPageQueue.push(page);
    }

    using PageBuffer = gsl::span<const std::byte>;

    bool mDebug;
    o2::mch::raw::PageDecoder mDecoder; ///< CRU page decoder
    SafeQueue<PageBuffer> mPageQueue;

    Elec2DetMapper mElec2Det{nullptr};       ///< front-end electronics mapping
    std::vector<MapSolar>& mMapSolar;

    std::thread mThread;
    std::mutex& mDigitsMutex;
    RawDigitVector& mDigits;
    Condition& mProcessEndCondition;
  };
  std::vector<PageProcessor> mDecoderThreadPool;
  Condition mProcessEndCondition;
  int mThreadsNum{2};

  o2::mch::raw::PageDecoder mDecoder; ///< CRU page decoder

  std::mutex mDigitsMutex;
  RawDigitVector mDigits;                               ///< vector of decoded digits
  std::unordered_set<OrbitInfo, OrbitInfoHash> mOrbits; ///< list of orbits in the processed buffer

  std::optional<uint32_t> mFirstOrbitInRun; ///< first orbit in the processed run
  SampaTimeFrameStart mSampaTimeFrameStart; ///< SAMPA bunch-crossing counter at the beiginning of the TF

  uint32_t mSampaTimeOffset{339986}; ///< SAMPA BC counter value at the beginning of the first orbit in the run

  SampaChannelHandler mChannelHandler;                  ///< optional user function to be called for each decoded SAMPA hit
  std::function<void(o2::header::RDHAny*)> mRdhHandler; ///< optional user function to be called for each RDH

  bool mDebug{false};
  bool mDs2manu{false};
  uint32_t mOrbit{0};
  bool mUseDummyElecMap{false};
};

bool operator<(const DataDecoder::RawDigit& d1, const DataDecoder::RawDigit& d2);

std::ostream& operator<<(std::ostream& os, const DataDecoder::RawDigit& d);

std::string asString(const DataDecoder::RawDigit& d);

} // namespace raw
} // namespace mch
} // end namespace o2
#endif //O2_MCH_DATADECODER_H_
