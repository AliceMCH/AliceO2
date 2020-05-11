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
/// \file    DatDecoderSpec.cxx
/// \author  Andrea Ferrero
///
/// \brief Implementation of a data processor to run the raw decoding
///

#include <random>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include "Framework/CallbackService.h"
#include "Framework/ConfigParamRegistry.h"
#include "Framework/ControlService.h"
#include "Framework/DataProcessorSpec.h"
#include "Framework/Lifetime.h"
#include "Framework/Output.h"
#include "Framework/Task.h"
#include "Framework/WorkflowSpec.h"

#include "DPLUtils/DPLRawParser.h"
#include "MCHBase/Digit.h"
#include "Headers/RAWDataHeader.h"
#include "MCHRawCommon/DataFormats.h"
#include "MCHRawDecoder/PageDecoder.h"
#include "MCHRawElecMap/Mapper.h"
#include "MCHRawCommon/RDHManip.h"
#include "MCHMappingInterface/Segmentation.h"
#include "MCHWorkflow/DataDecoderSpec.h"
#include <array>

namespace o2::header
{
extern std::ostream& operator<<(std::ostream&, const o2::header::RAWDataHeaderV4&);
}

namespace o2
{
namespace mch
{
namespace raw
{

using namespace o2;
using namespace o2::framework;
using namespace o2::mch::mapping;
using RDHv4 = o2::header::RAWDataHeaderV4;

std::array<int, 64> refManu2ds_st345_v5 = {
  63, 62, 61, 60, 59, 57, 56, 53, 51, 50, 47, 45, 44, 41, 38, 35,
  36, 33, 34, 37, 32, 39, 40, 42, 43, 46, 48, 49, 52, 54, 55, 58,
  7, 8, 5, 2, 6, 1, 3, 0, 4, 9, 10, 15, 17, 18, 22, 25,
  31, 30, 29, 28, 27, 26, 24, 23, 20, 21, 16, 19, 12, 14, 11, 13};

std::array<int, 64> refManu2ds_st345_v2 = {
  62, 61, 63, 60, 59, 55, 58, 57, 56, 54, 50, 46, 42, 39, 37, 41,
  35, 36, 33, 34, 32, 38, 43, 40, 45, 44, 47, 48, 49, 52, 51, 53,
  7, 6, 5, 4, 2, 3, 1, 0, 9, 11, 13, 15, 17, 19, 21, 23,
  31, 30, 29, 28, 27, 26, 25, 24, 22, 20, 18, 16, 14, 12, 10, 8};

#define refManu2ds_st345 refManu2ds_st345_v5

std::array<int, 64> refDs2manu_st345;

int manu2ds(int i)
{
  return refManu2ds_st345[i];
}

int ds2manu(int i)
{
  return refDs2manu_st345[i];
}

bool mDs2manu = false; ///< convert channel numbering from Run3 to Run1-2 order
bool mPrint = false;   ///< print digits

class MergerDigit
{
 public:
  MergerDigit() = default;
  ~MergerDigit() = default;

  o2::mch::Digit digit;
  o2::mch::HitTime stopTime;
  bool merged = {false};
  int solarId = {-1};
  int dsAddr = {-1};
  int chAddr = {-1};
};

struct MergerBuffer {
  std::vector<MergerDigit> digits;
  uint32_t orbit = {0};
};

class FeeIdMerger
{
  MergerBuffer buffers[2];
  int currentBufId = {1};
  int previousBufId = {0};
  int feeId = {0};

  std::function<void(const Digit&)> sendDigit;

 public:
  FeeIdMerger() = default;
  ~FeeIdMerger() = default;

  void setId(int id)
  {
    feeId = id;
  }

  void setDigitHandler(std::function<void(const Digit&)> h)
  {
    sendDigit = h;
  }

  void setOrbit(uint32_t orbit);

  MergerBuffer& getCurrentBuffer()
  {
    return buffers[currentBufId];
  }

  MergerBuffer& getPreviousBuffer()
  {
    return buffers[previousBufId];
  }

  int getCurrentBufId()
  {
    return currentBufId;
  }

  int getPreviousBufId()
  {
    return previousBufId;
  }

  void mergeDigits();

  void storeDigits();
};

void FeeIdMerger::setOrbit(uint32_t orbit)
{
  if (orbit == buffers[currentBufId].orbit) {
    return;
  }

  int nSent = 0;
  for (auto& d : buffers[previousBufId].digits) {
    if (!d.merged) {
    sendDigit(d.digit);
    nSent += 1;
    }
  }
  if (mPrint) {
    std::cout<<"[FeeIdMerger] sent "<<nSent<<" digits for orbit "<<buffers[previousBufId].orbit<<"  current orbit is "<<orbit<<std::endl;
  }

  currentBufId = 1 - currentBufId;
  previousBufId = 1 - previousBufId;
  buffers[currentBufId].digits.clear();
  buffers[currentBufId].orbit = orbit;
}


void FeeIdMerger::mergeDigits()
{
  auto checkDigits = [](MergerDigit& d1, MergerDigit& d2) -> bool {
    // skip digits that are already merged
    //std::cout<<"merged: "<<d1.merged<<","<<d2.merged<<std::endl;
    if (d1.merged || d2.merged)
      return false;

    // skip all digits not matching the detId/padId
    //std::cout<<"detId: "<<d1.digit.getDetID()<<","<<d2.digit.getDetID()<<std::endl;
    if (d1.solarId != d2.solarId)
      return false;
    if (d1.dsAddr != d2.dsAddr)
      return false;
    if (d1.chAddr != d2.chAddr)
      return false;
    if (d1.digit.getDetID() != d2.digit.getDetID())
      return false;
    //std::cout<<"padId: "<<d1.digit.getPadID()<<","<<d2.digit.getPadID()<<std::endl;
    if (d1.digit.getPadID() != d2.digit.getPadID())
      return false;

    // compute time difference
    HitTime startTime = d1.digit.getTime();
    uint32_t bxStart = startTime.bunchCrossing;
    HitTime stopTime = d2.stopTime;
    uint32_t bxStop = stopTime.bunchCrossing;
    // correct for value rollover
    if (bxStart < bxStop)
      bxStart += 0x100000;

    uint32_t stopTimeFull = bxStop + (stopTime.sampaTime << 2);
    uint32_t startTimeFull = bxStart + (startTime.sampaTime << 2);
    uint32_t timeDiff = startTimeFull - stopTimeFull;

    if (mPrint) {
      std::cout << "checkDigits: " << d1.digit.getDetID() << "  " << d1.digit.getPadID() << "  " << bxStop << "  " << stopTime.sampaTime
                << "  " << bxStart << "  " << startTime.sampaTime << "  " << timeDiff << std::endl;
    }

    // skip if the time difference is not equal to 1 ADC clock cycle
    if (timeDiff > 8)
      return false;

    // merge digits
    d2.digit.setADC(d1.digit.getADC() + d2.digit.getADC());
    d2.stopTime = d1.stopTime;
    d1.merged = true;
    return true;
  };

  auto& currentBuffer = getCurrentBuffer();
  auto currentBufId = getCurrentBufId();
  auto& previousBuffer = getPreviousBuffer();
  auto previousBufId = getPreviousBufId();

  if (mPrint) {
    std::cout << "Merging digits in " << previousBufId << " (orbit=" << previousBuffer.orbit << ")\n";
  }
  for (size_t i = 0; i < previousBuffer.digits.size(); i++) {
    MergerDigit& d1 = previousBuffer.digits[i];

    // skip digits that do not start at the beginning of the time window
    HitTime startTime = d1.digit.getTime();
    if (startTime.sampaTime != 0) {
      continue;
    }

    for (size_t j = 0; j < previousBuffer.digits.size(); j++) {
      if (i == j) {
        continue;
      }
      MergerDigit& d2 = previousBuffer.digits[j];
      if (checkDigits(d1, d2)) {
        break;
      }
    }
  }

  //if(mPrint) {
  //  std::cout<<"digitsBuffer[currentBufId].orbit: "<<digitsBuffer[currentBufId].orbit<<"\ndigitsBuffer[previousBufId].orbit: "<<digitsBuffer[previousBufId].orbit<<std::endl;
  //}

  // only merge digits from consecutive orbits
  uint32_t orbit_p = previousBuffer.orbit;
  uint32_t orbit_c = currentBuffer.orbit;
  if (mPrint) {
    std::cout << "orbit_c: " << orbit_c << "  orbit_p: " << orbit_p << std::endl;
  }
  if ((orbit_c >= orbit_p) && ((orbit_c - orbit_p) > 1)) {
    return;
  }

  if (mPrint) {
    std::cout << "Merging digits from " << currentBufId << " (orbit=" << currentBuffer.orbit << ") into "
              << previousBufId << " (orbit=" << previousBuffer.orbit << ")\n";
  }
  for (size_t i = 0; i < currentBuffer.digits.size(); i++) {
    MergerDigit& d1 = currentBuffer.digits[i];

    // skip digits that do not start at the beginning of the time window
    HitTime startTime = d1.digit.getTime();
    if (startTime.sampaTime != 0) {
      continue;
    }

    for (size_t j = 0; j < previousBuffer.digits.size(); j++) {
      MergerDigit& d2 = previousBuffer.digits[j];
      if (checkDigits(d1, d2)) {
        break;
      }
    }
  }
}



void FeeIdMerger::storeDigits()
{
  auto& digits = getPreviousBuffer().digits;
  if (mPrint) {
    std::cout << "[FeeIdMerger] storing " << digits.size() << " digits from " << feeId << "-" << getPreviousBufId()
        << " (orbit=" << getPreviousBuffer().orbit << ")" << std::endl;
  }

  for (auto& d : digits) {
    if (!d.merged && (d.digit.getPadID() >= 0)) {
      sendDigit(d.digit);
    }
  }
}


#define MCH_MERGER_FEEID_MAX 63

class MergerBase
{
 public:
  MergerBase() = default;
  ~MergerBase() = default;

  virtual void setOrbit(int feeId, uint32_t orbit) = 0;

  virtual void addDigit(int feeId, int solarId, int dsAddr, int chAddr,
                        int deId, int padId, int adc, HitTime time, HitTime stopTime) = 0;

  virtual void mergeDigits(int feeId) = 0;
};

class MergerSimple : public MergerBase
{
  //std::vector<Digit> digits;
  std::function<void(const Digit&)> sendDigit;

 public:
  MergerSimple() = default;
  ~MergerSimple() = default;

  void setOrbit(int feeId, uint32_t orbit)
  {
    //digits.clear();
  }

  void setDigitHandler(std::function<void(const Digit&)> h)
  {
    sendDigit = h;
  }

  void addDigit(int feeId, int solarId, int dsAddr, int chAddr,
                int deId, int padId, int adc, HitTime time, HitTime stopTime)
  {
    //digits.emplace_back(o2::mch::Digit(time, deId, padId, adc));
    sendDigit(o2::mch::Digit(time, deId, padId, adc));
  }

  void mergeDigits(int feeId) {}
};


class Merger : public MergerBase
{
  int32_t feeId = {-1};
  FeeIdMerger mergers[MCH_MERGER_FEEID_MAX + 1];

 public:
  Merger() = default;
  ~Merger() = default;

  void setOrbit(int feeId, uint32_t orbit);

  void setDigitHandler(std::function<void(const Digit&)> h);

  void addDigit(int feeId, int solarId, int dsAddr, int chAddr,
                int deId, int padId, int adc, HitTime time, HitTime stopTime);

  void mergeDigits(int feeId);
};


void Merger::setOrbit(int feeId, uint32_t orbit)
{
  if (feeId < 0 || feeId > MCH_MERGER_FEEID_MAX) {
    return;
  }

  mergers[feeId].setOrbit(orbit);
}

void Merger::setDigitHandler(std::function<void(const Digit&)> h)
{
  for (int feeId = 0; feeId <= MCH_MERGER_FEEID_MAX; feeId++) {
    mergers[feeId].setDigitHandler(h);
  }
}

void Merger::addDigit(int feeId, int solarId, int dsAddr, int chAddr,
                      int deId, int padId, int adc, HitTime time, HitTime stopTime)
{
  if (feeId < 0 || feeId > MCH_MERGER_FEEID_MAX) {
    return;
  }

  mergers[feeId].getCurrentBuffer().digits.emplace_back(MergerDigit{o2::mch::Digit(time, deId, padId, adc),
                                                                    stopTime, false, solarId, dsAddr, chAddr});
}

void Merger::mergeDigits(int feeId)
{
  if (feeId < 0 || feeId > MCH_MERGER_FEEID_MAX) {
    return;
  }

  mergers[feeId].mergeDigits();
}

/*void Merger::storeDigits(int feeId, std::vector<Digit>& outputDigits)
{
  if (feeId < 0 || feeId > MCH_MERGER_FEEID_MAX) {
    return;
  }

  auto& digits = mergers[feeId].getPreviousBuffer().digits;
  if (mPrint) {
    std::cout << "[Merger] storing " << digits.size() << " digits from " << feeId << "-" << mergers[feeId].getPreviousBufId()
        << " (orbit=" << mergers[feeId].getPreviousBuffer().orbit << ")" << std::endl;
  }
  int nDigits = {0};
  for (auto& d : digits) {
    if (!d.merged && (d.digit.getPadID() >= 0)) {
      outputDigits.emplace_back(d.digit);
    }
  }
}*/

//=======================
// Data decoder
class DataDecoderTask
{
private:
  void decodeBuffer(gsl::span<const std::byte> page)
  {
    size_t ndigits{0};

    uint8_t isStopRDH = 0;
    uint32_t orbit;
    uint32_t feeId;
    uint32_t linkId;

    auto channelHandler = [&](DsElecId dsElecId, uint8_t channel, o2::mch::raw::SampaCluster sc) {
      if (mDs2manu) {
        channel = ds2manu(int(channel));
      }
      if (mPrint) {
        auto s = asString(dsElecId);
        auto ch = fmt::format("{}-CH{:02d}", s, channel);
        std::cout << ch << std::endl;
      }
      double digitadc(0);
      //for (auto d = 0; d < sc.nofSamples(); d++) {
      for (auto d = 0; d < sc.samples.size(); d++) {
        digitadc += sc.samples[d];
      }

      int deId{-1};
      int dsIddet{-1};
      if (auto opt = mElec2Det(dsElecId); opt.has_value()) {
        DsDetId dsDetId = opt.value();
        dsIddet = dsDetId.dsId();
        deId = dsDetId.deId();
      }

      int padId = -1;
      const Segmentation& segment = segmentation(deId);
      if ((&segment) == nullptr) {
        return;
      }
      padId = segment.findPadByFEE(dsIddet, int(channel));
      if (mPrint) {
        auto s = asString(dsElecId);
        auto ch = fmt::format("{}-CH{:02d}", s, channel);
        std::cout << ch << "  "
            << fmt::format("PAD ({:04d} {:04d} {:04d})\tADC {:5.0f}\tTIME ({} {} {})\tSIZE {}\tEND {}", deId, dsIddet, padId, digitadc, orbit, sc.bunchCrossing, sc.timestamp, sc.nofSamples(), (sc.timestamp + sc.nofSamples() - 1))
        << (((sc.timestamp + sc.nofSamples() - 1) >= 98) ? " *" : "") << std::endl;
        //std::cout << "DS " << (int)dsElecId.elinkId() << "  CHIP " << ((int)channel) / 32 << "  CH " << ((int)channel) % 32 << "  ADC " << digitadc << "  DE# " << deId << "  DSid " << dsIddet << "  PadId " << padId << std::endl;
        if (false && padId == 994) {
          for (auto d = 0; d < sc.samples.size(); d++) {
            std::cout<<"  sample "<<d<<"  "<<sc.samples[d]<<std::endl;
          }
        }
      }

      HitTime time;
      time.sampaTime = sc.timestamp;
      time.bunchCrossing = sc.bunchCrossing;
      time.orbit = orbit;
      HitTime stopTime;
      stopTime.sampaTime = sc.timestamp + sc.nofSamples() - 1;
      stopTime.bunchCrossing = sc.bunchCrossing;
      stopTime.orbit = orbit;

      mMerger->addDigit(feeId, static_cast<int>(dsElecId.solarId()), static_cast<int>(dsElecId.elinkId()), static_cast<int>(channel),
                        deId, padId, digitadc, time, stopTime);

      if (false && mPrint)
        std::cout << "DIGIT STORED:\nADC " << digitadc << " DE# " << (int)deId << " PadId " << (int)padId
                  << " time " << sc.timestamp << " size " << sc.nofSamples() << std::endl;
      ++ndigits;
    };

    const auto patchPage = [&](gsl::span<const std::byte> rdhBuffer) {
      auto rdhPtr = reinterpret_cast<o2::header::RAWDataHeaderV4*>(const_cast<std::byte*>(&rdhBuffer[0]));
      auto& rdh = *rdhPtr;
      mNrdhs++;
      auto cruId = rdhCruId(rdh);
      rdhFeeId(rdh, cruId * 2 + rdhEndpoint(rdh));
      isStopRDH = rdhStop(rdh);
      if (true && mPrint) {
        std::cout << std::endl
                  << mNrdhs << "--" << rdh << "\n";
      }

      feeId = rdhFeeId(rdh);
      orbit = rdhOrbit(rdh);
      linkId = rdhLinkId(rdh);

      if (linkId == 15) {
        mMerger = &mMergerFull;
      } else {
        mMerger = &mMergerSimple;
      }

      if (isStopRDH == 0 && (feeId < 64)) {
        mMerger->setOrbit(feeId, orbit);
      }
    };

    if (!mDecoder) {
      mDecoder = mFee2Solar ? o2::mch::raw::createPageDecoder(page, channelHandler, mFee2Solar)
      : o2::mch::raw::createPageDecoder(page, channelHandler);
    }

    patchPage(page);

    // skip stop RDHs
    if (isStopRDH != 0)
      return;

    mDecoder(page);

    mMerger->mergeDigits(feeId);

    if (mPrint) {
      for (auto d : mOutputDigits) {
        if(d.getPadID() < 0) continue;
        const Segmentation& segment = segmentation(d.getDetID());
        float X = segment.padPositionX(d.getPadID());
        float Y = segment.padPositionY(d.getPadID());
        bool bend = !segment.isBendingPad(d.getPadID());
        if(bend) continue;
        std::cout << fmt::format("  DE {:4d}  PAD {:5d}  ADC {:6d}  TIME ({} {} {:4d})",
            d.getDetID(), d.getPadID(), d.getADC(), d.getTime().orbit, d.getTime().bunchCrossing, d.getTimeStamp());
        std::cout << fmt::format("\tC {}  PAD_XY {:+2.2f} , {:+2.2f}", (int)bend, X, Y);
        std::cout << std::endl;
      }
      for (auto d : mOutputDigits) {
        if(d.getPadID() < 0) continue;
        const Segmentation& segment = segmentation(d.getDetID());
        float X = segment.padPositionX(d.getPadID());
        float Y = segment.padPositionY(d.getPadID());
        bool bend = !segment.isBendingPad(d.getPadID());
        if(!bend) continue;
        std::cout << fmt::format("  DE {:4d}  PAD {:5d}  ADC {:6d}  TIME ({} {} {:4d})",
            d.getDetID(), d.getPadID(), d.getADC(), d.getTime().orbit, d.getTime().bunchCrossing, d.getTimeStamp());
        std::cout << fmt::format("\tC {}  PAD_XY {:+2.2f} , {:+2.2f}", (int)bend, X, Y);
        std::cout << std::endl;
      }
      for (auto d : mOutputDigits) {
        if(d.getPadID() >= 0) continue;
        std::cout << fmt::format("  DE {:4d}  PAD {:5d}  ADC {:6d}  TIME ({} {} {:4d})",
            d.getDetID(), d.getPadID(), d.getADC(), d.getTime().orbit, d.getTime().bunchCrossing, d.getTimeStamp());
        std::cout << std::endl;
      }
    }
  }

 private:
  std::string readFileContent(std::string& filename)
  {
    std::string content;
    std::string s;
    std::ifstream in(filename);
    while (std::getline(in, s)) {
      content += s;
      content += " ";
    }
    std::cout << "readFileContent(" << filename << "):" << std::endl
              << content << std::endl;
    return content;
  }

  void initElec2DetMapper(std::string filename)
  {
    std::cout << "[initElec2DetMapper] filename=" << filename << std::endl;
    if (filename.empty()) {
      mElec2Det = createElec2DetMapper<ElectronicMapperGenerated>();
    } else {
      ElectronicMapperString::sFecMap = readFileContent(filename);
      mElec2Det = createElec2DetMapper<ElectronicMapperString>();
    }
  }

  void initFee2SolarMapper(std::string filename)
  {
    std::cout << "[initFee2SolarMapper] filename=" << filename << std::endl;
    if (filename.empty()) {
      mFee2Solar = createFeeLink2SolarMapper<ElectronicMapperGenerated>();
    } else {
      ElectronicMapperString::sCruMap = readFileContent(filename);
      mFee2Solar = createFeeLink2SolarMapper<ElectronicMapperString>();
    }
  }

 public:
  //_________________________________________________________________________________________________
  void init(framework::InitContext& ic)
  {
    mNrdhs = 0;

    for (int i = 0; i < 64; i++) {
      for (int j = 0; j < 64; j++) {
        if (refManu2ds_st345[j] != i) {
          continue;
        }
        refDs2manu_st345[i] = j;
        break;
      }
    }

    mDs2manu = ic.options().get<bool>("ds2manu");
    mPrint = ic.options().get<bool>("print");

    auto mapCRUfile = ic.options().get<std::string>("cru-map");
    auto mapFECfile = ic.options().get<std::string>("fec-map");

    initElec2DetMapper(mapFECfile);
    initFee2SolarMapper(mapCRUfile);

    const auto storeDigit = [&](const Digit& d) {
      mOutputDigits.emplace_back(d);
    };

    mMergerSimple.setDigitHandler(storeDigit);
    mMergerFull.setDigitHandler(storeDigit);
  }

  //_________________________________________________________________________________________________
  void decodeTF(framework::ProcessingContext& pc)
  {
    // get the input buffer
    auto& inputs = pc.inputs();
    DPLRawParser parser(inputs, o2::framework::select("TF:MCH/RAWDATA"));

    for (auto it = parser.begin(), end = parser.end(); it != end; ++it) {
      // retrieving RDH v4
      auto const* rdh = it.get_if<o2::header::RAWDataHeaderV4>();
      // retrieving the raw pointer of the page
      auto const* raw = it.raw();
      // size of payload
      size_t payloadSize = it.size();

      if (payloadSize == 0) {
        continue;
      }

      gsl::span<const std::byte> buffer(reinterpret_cast<const std::byte*>(raw), sizeof(o2::header::RAWDataHeaderV4) + payloadSize);
      decodeBuffer(buffer);
    }
  }

  //_________________________________________________________________________________________________
  void decodeReadout(const o2::framework::DataRef& input)
  {
    static int nFrame = 1;
    // get the input buffer
    if (input.spec->binding != "readout") {
      return;
    }

    const auto* header = o2::header::get<header::DataHeader*>(input.header);
    if (!header) {
      return;
    }

    auto const* raw = input.payload;
    // size of payload
    size_t payloadSize = header->payloadSize;

    if (mPrint) {
      std::cout << nFrame << "  payloadSize=" << payloadSize << std::endl;
    }
    nFrame += 1;
    if (payloadSize == 0) {
      return;
    }

    gsl::span<const std::byte> buffer(reinterpret_cast<const std::byte*>(raw), payloadSize);
    decodeBuffer(buffer);
  }

  //_________________________________________________________________________________________________
  void run(framework::ProcessingContext& pc)
  {
    mOutputDigits.clear();
    //std::cout << "DatDecoderTask::run()" << std::endl;
    decodeTF(pc);
    for (auto&& input : pc.inputs()) {
      if (input.spec->binding == "readout")
        decodeReadout(input);
    }

    const size_t OUT_SIZE = sizeof(o2::mch::Digit) * mOutputDigits.size();

    // send the output buffer via DPL
    char* outbuffer = nullptr;
    outbuffer = (char*)realloc(outbuffer, OUT_SIZE);
    memcpy(outbuffer, mOutputDigits.data(), OUT_SIZE);

    // create the output message
    auto freefct = [](void* data, void*) { free(data); };
    pc.outputs().adoptChunk(Output{"MCH", "DIGITS", 0}, outbuffer, OUT_SIZE, freefct, nullptr);
  }

 private:
  Elec2DetMapper mElec2Det{nullptr};
  FeeLink2SolarMapper mFee2Solar{nullptr};
  o2::mch::raw::PageDecoder mDecoder;
  size_t mNrdhs{0};
  std::vector<o2::mch::Digit> mOutputDigits;

  std::ifstream mInputFile{}; ///< input file

  MergerBase* mMerger = {nullptr};
  MergerSimple mMergerSimple;
  Merger mMergerFull;
};

//_________________________________________________________________________________________________
o2::framework::DataProcessorSpec getDecodingSpec()
{
  return DataProcessorSpec{
    "DataDecoder",
    //o2::framework::select("TF:MCH/RAWDATA, re:ROUT/RAWDATA"),
    o2::framework::select("readout:ROUT/RAWDATA"),
    //o2::framework::select("TF:MCH/RAWDATA"),
    Outputs{OutputSpec{"MCH", "DIGITS", 0, Lifetime::Timeframe}},
    AlgorithmSpec{adaptFromTask<DataDecoderTask>()},
    Options{{"print", VariantType::Bool, false, {"print digits"}},
            {"cru-map", VariantType::String, "", {"custom CRU mapping"}},
            {"fec-map", VariantType::String, "", {"custom FEC mapping"}},
            {"ds2manu", VariantType::Bool, false, {"convert channel numbering from Run3 to Run1-2 order"}}}};
}

} // namespace raw
} // namespace mch
} // end namespace o2
