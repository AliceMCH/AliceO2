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

std::array<int, 64> refManu2ds_st345 = {
  63, 62, 61, 60, 59, 57, 56, 53, 51, 50, 47, 45, 44, 41, 38, 35,
  36, 33, 34, 37, 32, 39, 40, 42, 43, 46, 48, 49, 52, 54, 55, 58,
  7, 8, 5, 2, 6, 1, 3, 0, 4, 9, 10, 15, 17, 18, 22, 25,
  31, 30, 29, 28, 27, 26, 24, 23, 20, 21, 16, 19, 12, 14, 11, 13};

std::array<int, 64> refManu2ds_st345_v2 = {
  62, 61, 63, 60, 59, 55, 58, 57, 56, 54, 50, 46, 42, 39, 37, 41,
  35, 36, 33, 34, 32, 38, 43, 40, 45, 44, 47, 48, 49, 52, 51, 53,
  7, 6, 5, 4, 2, 3, 1, 0, 9, 11, 13, 15, 17, 19, 21, 23,
  31, 30, 29, 28, 27, 26, 25, 24, 22, 20, 18, 16, 14, 12, 10, 8};

std::array<int, 64> refDs2manu_st345;

int manu2ds(int i)
{
  return refManu2ds_st345[i];
}

int ds2manu(int i)
{
  return refDs2manu_st345[i];
}

class DigitInfo
{
 public:
  DigitInfo() = default;
  ~DigitInfo() = default;

  o2::mch::Digit digit;
  o2::mch::HitTime stopTime;
  bool merged = {false};
  int solarId = {-1};
  int dsAddr = {-1};
  int chAddr = {-1};
};

struct DigitsBuffer {
  std::vector<DigitInfo> digits;
  uint32_t orbit;
};

class DigitsMerger
{
  DigitsBuffer buffers[2];
  int currentBufId = {1};
  int previousBufId = {0};

 public:
  DigitsMerger() = default;
  ~DigitsMerger() = default;

  void newPage(uint32_t orbit)
  {
    currentBufId = 1 - currentBufId;
    previousBufId = 1 - previousBufId;
    buffers[currentBufId].digits.clear();
    buffers[currentBufId].orbit = orbit;
  }

  DigitsBuffer& getCurrentBuffer()
  {
    return buffers[currentBufId];
  }

  DigitsBuffer& getPreviousBuffer()
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
};

//=======================
// Data decoder
class DataDecoderTask
{
  void mergeDigits()
  {
    auto checkDigits = [](DigitInfo& d1, DigitInfo& d2) -> bool {
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

      if (false) {
        std::cout << "checkDigits: " << d1.digit.getDetID() << "  " << d1.digit.getPadID() << "  " << bxStop << "  " << stopTime.sampaTime
                  << "  " << bxStart << "  " << startTime.sampaTime << "  " << timeDiff << std::endl;
      }

      // skip if the time difference is not equal to 1 ADC clock cycle
      if (timeDiff > 8)
        return false;

      // merge digits
      d2.digit.setADC(d1.digit.getADC() + d2.digit.getADC());
      d1.merged = true;
      return true;
    };

    auto& currentBuffer = mergers[feeId].getCurrentBuffer();
    auto currentBufId = mergers[feeId].getCurrentBufId();
    auto& previousBuffer = mergers[feeId].getPreviousBuffer();
    auto previousBufId = mergers[feeId].getPreviousBufId();

    if (mPrint) {
      std::cout << "Merging digits in " << previousBufId << " (orbit=" << previousBuffer.orbit << ")\n";
    }
    for (size_t i = 0; i < previousBuffer.digits.size(); i++) {
      DigitInfo& d1 = previousBuffer.digits[i];

      // skip digits that do not start at the beginning of the time window
      HitTime startTime = d1.digit.getTime();
      if (startTime.sampaTime != 0) {
        continue;
      }

      for (size_t j = 0; j < previousBuffer.digits.size(); j++) {
        if (i == j) {
          continue;
        }
        DigitInfo& d2 = previousBuffer.digits[j];
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
      DigitInfo& d1 = currentBuffer.digits[i];

      // skip digits that do not start at the beginning of the time window
      HitTime startTime = d1.digit.getTime();
      if (startTime.sampaTime != 0) {
        continue;
      }

      for (size_t j = 0; j < previousBuffer.digits.size(); j++) {
        DigitInfo& d2 = previousBuffer.digits[j];
        if (checkDigits(d1, d2)) {
          break;
        }
      }
    }
  }

  void decodeBuffer(gsl::span<const std::byte> page)
  {
    size_t ndigits{0};

    uint8_t isStopRDH = 0;

    auto channelHandler = [&](DsElecId dsElecId, uint8_t channel, o2::mch::raw::SampaCluster sc) {
      auto& digits = mergers[feeId].getCurrentBuffer().digits;
      if (mDs2manu) {
        channel = ds2manu(int(channel));
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

      if (mPrint) {
        auto s = asString(dsElecId);
        auto ch = fmt::format("{}-CH{:02d}", s, channel);
        std::cout << ch << "  "
                  << fmt::format("PAD {:04d}-{:04d}\tADC {:5.0f}\tTIME {}-{}-{}\tSIZE {}\tEND {}", deId, dsIddet, digitadc, mergers[feeId].getCurrentBuffer().orbit, sc.bunchCrossing, sc.timestamp, sc.nofSamples(), (sc.timestamp + sc.nofSamples() - 1))
                  << (((sc.timestamp + sc.nofSamples() - 1) >= 98) ? " *" : "") << std::endl;
        //std::cout << "DS " << (int)dsElecId.elinkId() << "  CHIP " << ((int)channel) / 32 << "  CH " << ((int)channel) % 32 << "  ADC " << digitadc << "  DE# " << deId << "  DSid " << dsIddet << "  PadId " << padId << std::endl;
      }

      int padId = -1;
      try {
        const Segmentation& segment = segmentation(deId);
        padId = segment.findPadByFEE(dsIddet, int(channel));
        if (mPrint) {
          auto s = asString(dsElecId);
          auto ch = fmt::format("{}-CH{:02d}", s, channel);
          std::cout << ch << "  "
                    << fmt::format("PAD {:04d}-{:04d}-{:04d}\tADC {:5.0f}\tTIME {}-{}-{}\tSIZE {}\tEND {}", deId, dsIddet, padId, digitadc, mergers[feeId].getCurrentBuffer().orbit, sc.bunchCrossing, sc.timestamp, sc.nofSamples(), (sc.timestamp + sc.nofSamples() - 1))
                    << (((sc.timestamp + sc.nofSamples() - 1) >= 98) ? " *" : "") << std::endl;
          //std::cout << "DS " << (int)dsElecId.elinkId() << "  CHIP " << ((int)channel) / 32 << "  CH " << ((int)channel) % 32 << "  ADC " << digitadc << "  DE# " << deId << "  DSid " << dsIddet << "  PadId " << padId << std::endl;
        }
      } catch (const std::exception& e) {
        std::cout << "Failed to get padId: " << e.what() << std::endl;
        return;
      }

      HitTime time;
      time.sampaTime = sc.timestamp;
      time.bunchCrossing = sc.bunchCrossing;
      HitTime stopTime;
      stopTime.sampaTime = sc.timestamp + sc.nofSamples() - 1;
      stopTime.bunchCrossing = sc.bunchCrossing;

      digits.emplace_back(DigitInfo{o2::mch::Digit(time, deId, padId, digitadc),
                                    stopTime, false, static_cast<int>(dsElecId.solarId()), static_cast<int>(dsElecId.elinkId()), static_cast<int>(channel)});

      if (false && mPrint)
        std::cout << "DIGIT STORED:\nADC " << digits.back().digit.getADC() << " DE# " << digits.back().digit.getDetID() << " PadId " << digits.back().digit.getPadID()
                  << " time " << digits.back().digit.getTimeStamp() << " size " << sc.nofSamples() << std::endl;
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
      if (isStopRDH == 0 && (feeId < 64)) {
        mergers[feeId].newPage(rdhOrbit(rdh));
      }
    };

    o2::mch::raw::PageDecoder decode =
      mFee2Solar ? o2::mch::raw::createPageDecoder(page, channelHandler, mFee2Solar)
                 : o2::mch::raw::createPageDecoder(page, channelHandler);

    patchPage(page);
    //std::cout<<"isStopRDH: "<<(int)isStopRDH<<std::endl;
    // skip stop RDHs
    if (isStopRDH != 0)
      return;

    auto& currentBuffer = mergers[feeId].getCurrentBuffer();
    auto currentBufId = mergers[feeId].getCurrentBufId();
    auto& previousBuffer = mergers[feeId].getPreviousBuffer();
    auto previousBufId = mergers[feeId].getPreviousBufId();

    if (mPrint) {
      std::cout << "decoding page: previousBufId=" << previousBufId << "  previous orbit: "
                << previousBuffer.orbit << "  digits.size(): " << previousBuffer.digits.size() << std::endl;
      std::cout << "               currentBufId=" << currentBufId << "  current orbit: "
                << currentBuffer.orbit << "  digits.size(): " << currentBuffer.digits.size() << std::endl;
    }
    decode(page);
    if (mPrint) {
      std::cout << "page decoded: previousBufId=" << previousBufId << "  previous orbit: "
                << previousBuffer.orbit << "  digits.size(): " << previousBuffer.digits.size() << std::endl;
      std::cout << "               currentBufId=" << currentBufId << "  current orbit: "
                << currentBuffer.orbit << "  digits.size(): " << currentBuffer.digits.size() << std::endl;
    }

    mergeDigits();

    for (size_t i = 0; i < previousBuffer.digits.size(); i++) {
      Digit& d = previousBuffer.digits[i].digit;
      if (previousBuffer.digits[i].merged)
        continue;
      outputDigits.emplace_back(d);
    }

    if (mPrint) {
      std::cout << "previousBufId: " << previousBufId << "  digits.size(): " << previousBuffer.digits.size()
                << "  orbit: " << previousBuffer.orbit << std::endl;
      for (auto d : outputDigits) {
        std::cout << "  DE# " << d.getDetID() << " PadId " << d.getPadID() << " ADC " << d.getADC() << " time " << d.getTimeStamp() << std::endl;
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

    if (mPrint) {
      try {
        int deId = 819;
        int dsIddet = 1133;
        const Segmentation& segment = segmentation(deId);
        for (int channel = 0; channel < 64; channel++) {
          int padId = segment.findPadByFEE(dsIddet, int(channel));
          std::cout << "DE# " << deId << "  DSid " << dsIddet << "  channel " << channel << "  PadId " << padId << std::endl;
        }
      } catch (const std::exception& e) {
        std::cout << "Failed to get padId: " << e.what() << std::endl;
      }
    }
  }

  //_________________________________________________________________________________________________
  void
    decodeTF(framework::ProcessingContext& pc)
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
    outputDigits.clear();
    //std::cout << "DatDecoderTask::run()" << std::endl;
    decodeTF(pc);
    for (auto&& input : pc.inputs()) {
      if (input.spec->binding == "readout")
        decodeReadout(input);
    }

    const size_t OUT_SIZE = sizeof(o2::mch::Digit) * outputDigits.size();

    // send the output buffer via DPL
    char* outbuffer = nullptr;
    outbuffer = (char*)realloc(outbuffer, OUT_SIZE);
    memcpy(outbuffer, outputDigits.data(), OUT_SIZE);

    // create the output message
    auto freefct = [](void* data, void*) { free(data); };
    pc.outputs().adoptChunk(Output{"MCH", "DIGITS", 0}, outbuffer, OUT_SIZE, freefct, nullptr);
  }

 private:
  Elec2DetMapper mElec2Det{nullptr};
  FeeLink2SolarMapper mFee2Solar{nullptr};
  size_t mNrdhs{0};
  std::vector<o2::mch::Digit> outputDigits;

  std::ifstream mInputFile{}; ///< input file
  bool mDs2manu = false;      ///< print convert channel numbering from Run3 to Run1-2 order
  bool mPrint = false;        ///< print digits

  uint32_t feeId = 0;
  DigitsMerger mergers[64];
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
