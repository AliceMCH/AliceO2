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
std::array<int, 64> refDs2manu_st345;

int manu2ds(int i)
{
  return refManu2ds_st345[i];
}

int ds2manu(int i)
{
  return refDs2manu_st345[i];
}

//=============
// Classes for custom mapping implementation
#define MCH_MAX_FEEID 64
#define MCH_MAX_CRU_LINK 12
#define LINKID_MAX 0x7FF

class MapSolar
{
 public:
  int mLink = {-1}; // link ID

  MapSolar() = default;
  ~MapSolar() = default;
};


// CRU mapping
class MapCRU
{
  bool mInitialized = {false};
  MapSolar mSolarMap[MCH_MAX_FEEID][MCH_MAX_CRU_LINK];

 public:
  MapCRU() = default;
  ~MapCRU() = default;
  bool load(std::string mapFile);
  bool initialized() { return mInitialized; }
  int32_t getLink(int32_t c, int32_t l);
};


bool MapCRU::load(std::string mapFile)
{
  std::ifstream file;
  file.open(mapFile);
  if (!file) {
    std::cout << "[MapCRU::readMapping] can't open file " << mapFile << std::endl;
    return false;
  }

  int c, l, link_id;
  char tstr[500];
  while (file.getline(tstr, 499)) {
    std::string s(tstr);
    std::istringstream line(s);
    line >> link_id >> c >> l;
    if (c < 0 || c >= MCH_MAX_FEEID)
      continue;
    if (l < 0 || l >= MCH_MAX_CRU_LINK)
      continue;
    mSolarMap[c][l].mLink = link_id;
  }
  mInitialized = true;
  return true;
}

int32_t MapCRU::getLink(int32_t c, int32_t l)
{
  if (!initialized()) return -1;

  int32_t result = -1;
  if (c < 0 || c >= MCH_MAX_FEEID)
    return result;
  if (l < 0 || l >= MCH_MAX_CRU_LINK)
    return result;
  return mSolarMap[c][l].mLink;
}


class MapDualSampa
{
 public:
  int mDE = {-1};    // detector element
  int mIndex = {-1}; // DS index
  int mBad = {-1};   // if = 1 bad pad (not used for analysis)

  MapDualSampa() = default;
  ~MapDualSampa() = default;
};


// Electronics mapping
class MapFEC
{
  bool mInitialized = {false};
  MapDualSampa mDsMap[LINKID_MAX + 1][40];

 public:
  MapFEC() = default;
  ~MapFEC() = default;
  bool load(std::string mapFile);
  bool initialized() { return mInitialized; }
  bool getDsId(uint32_t link_id, uint32_t ds_addr, int& de, int& dsid);
};

bool MapFEC::load(std::string mapFile)
{
  std::ifstream file;
  file.open(mapFile);
  if (!file) {
    std::cout << "[MapFEC::readDSMapping] can't open file " << mapFile << std::endl;
    return false;
  }

  int link_id, group_id, de, ds_id[5];
  while (!file.eof()) {
    file >> link_id >> group_id >> de >> ds_id[0] >> ds_id[1] >> ds_id[2] >> ds_id[3] >> ds_id[4];
    if (link_id < 0 || link_id > LINKID_MAX)
      continue;
    for (int i = 0; i < 5; i++) {
      if (ds_id[i] <= 0)
        continue;
      int ds_addr = group_id * 5 + i;
      if (ds_addr < 0 || ds_addr >= 40)
        continue;
      mDsMap[link_id][ds_addr].mDE = de;
      mDsMap[link_id][ds_addr].mIndex = ds_id[i];
      mDsMap[link_id][ds_addr].mBad = 0;
    }
  }
  mInitialized = true;
  return true;
}

bool MapFEC::getDsId(uint32_t link_id, uint32_t ds_addr, int& de, int& dsid)
{
  if (!initialized()) return false;

  if (mDsMap[link_id][ds_addr].mBad == 1)
    return false;
  de = mDsMap[link_id][ds_addr].mDE;
  dsid = mDsMap[link_id][ds_addr].mIndex;
  return true;
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


struct DigitsBuffer
{
  std::vector<DigitInfo> digits;
  uint32_t orbit;
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
      if(d1.merged || d2.merged)
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
       if(bxStart < bxStop)
         bxStart += 0x100000;

       uint32_t stopTimeFull  = bxStop + (stopTime.sampaTime << 2);
       uint32_t startTimeFull = bxStart + (startTime.sampaTime << 2);
       uint32_t timeDiff = startTimeFull - stopTimeFull;

       std::cout<<d1.digit.getDetID()<<"  "<<d1.digit.getPadID()<<"  "<<bxStop<<"  "<<stopTime.sampaTime
           <<"  "<<bxStart<<"  "<<startTime.sampaTime<<"  "<<timeDiff<<std::endl;

       // skip if the time difference is not equal to 1 ADC clock cycle
       if(timeDiff > 8)
         return false;

       // merge digits
       d2.digit.setADC(d1.digit.getADC() + d2.digit.getADC());
       d1.merged = true;
       return true;
    };

    for(size_t i = 0; i < digitsBuffer[previousBufId].digits.size(); i++) {
      DigitInfo& d1 = digitsBuffer[previousBufId].digits[i];

      // skip digits that do not start at the beginning of the time window
      HitTime startTime = d1.digit.getTime();
      if(startTime.sampaTime != 0)
        continue;

      for(size_t j = 0; j < digitsBuffer[previousBufId].digits.size(); j++) {
        if(i == j)
          continue;

        DigitInfo& d2 = digitsBuffer[previousBufId].digits[j];
        if(checkDigits(d1, d2))
          break;
      }
    }

    // only merge digits from the same orbit
    if(digitsBuffer[currentBufId].orbit != digitsBuffer[previousBufId].orbit)
      return;

    for(size_t i = 0; i < digitsBuffer[currentBufId].digits.size(); i++) {
      DigitInfo& d1 = digitsBuffer[currentBufId].digits[i];

      // skip digits that do not start at the beginning of the time window
      HitTime startTime = d1.digit.getTime();
      if(startTime.sampaTime != 0)
        continue;

      for(size_t j = 0; j < digitsBuffer[previousBufId].digits.size(); j++) {
        DigitInfo& d2 = digitsBuffer[previousBufId].digits[j];
        if(checkDigits(d1, d2))
          break;
      }
    }
  }

  void decodeBuffer(gsl::span<const std::byte> page)
  {
    size_t ndigits{0};
    auto& digits = digitsBuffer[currentBufId].digits;

    auto linkHandler = [&](FeeLinkId feeLinkId) -> std::optional<uint16_t> {
      std::optional<uint16_t> result;
      uint16_t link = mMapCRU.getLink(feeLinkId.feeId(), feeLinkId.linkId());
      result = link;
      if (mPrint)
        std::cout << "[linkHandler] (" << (int)feeLinkId.feeId() << "," << (int)feeLinkId.linkId() << ") -> " << result.value() << std::endl;
      return result;
    };

    auto channelHandler = [&](DsElecId dsElecId, uint8_t channel, o2::mch::raw::SampaCluster sc) {
      if (mDs2manu)
        channel = ds2manu(int(channel));
      if (mPrint) {
        auto s = asString(dsElecId);
        auto ch = fmt::format("{}-CH{}", s, channel);
        std::cout << ch << std::endl;
      }
      double digitadc(0);
      //for (auto d = 0; d < sc.nofSamples(); d++) {
      for (auto d = 0; d < sc.samples.size(); d++) {
        digitadc += sc.samples[d];
      }

      int deId;
      int dsIddet;
      if (mMapFEC.initialized()) {
        if (!mMapFEC.getDsId(dsElecId.solarId(), dsElecId.elinkId(), deId, dsIddet)) {
          deId = dsIddet = -1;
        }
      } else if (auto opt = mElec2Det(dsElecId); opt.has_value()) {
        DsDetId dsDetId = opt.value();
        dsIddet = dsDetId.dsId();
        deId = dsDetId.deId();
      }

      int padId = -1;
      try {
        const Segmentation& segment = segmentation(deId);
        padId = segment.findPadByFEE(dsIddet, int(channel));
        if (mPrint)
          std::cout << "DS " << (int)dsElecId.elinkId() << "  CHIP " << ((int)channel) / 32 << "  CH " << ((int)channel) % 32 << "  ADC " << digitadc << "  DE# " << deId << "  DSid " << dsIddet << "  PadId " << padId << std::endl;
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

      if (mPrint)
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
      if (mPrint) {
        std::cout << std::endl << mNrdhs << "--" << rdh << "\n";
      }
      digitsBuffer[currentBufId].orbit = rdhOrbit(rdh);
    };

    o2::mch::raw::PageDecoder decode =
      mMapCRU.initialized() ? o2::mch::raw::createPageDecoder(page, channelHandler, linkHandler) : o2::mch::raw::createPageDecoder(page, channelHandler);
    patchPage(page);
    decode(page);

    mergeDigits();
  }

 public:
  //_________________________________________________________________________________________________
  void init(framework::InitContext& ic)
  {
    mElec2Det = createElec2DetMapper<ElectronicMapperGenerated>();
    mNrdhs = 0;

    for (int i = 0; i < 64; i++) {
      for (int j = 0; j < 64; j++) {
        if (refManu2ds_st345[j] != i)
          continue;
        refDs2manu_st345[i] = j;
        break;
      }
    }

    mDs2manu = ic.options().get<bool>("ds2manu");
    mPrint = ic.options().get<bool>("print");

    auto mapCRUfile = ic.options().get<std::string>("cru-map");
    auto mapFECfile = ic.options().get<std::string>("fec-map");

    if (!mapCRUfile.empty())
      mMapCRU.load(mapCRUfile);
    if (!mapFECfile.empty())
      mMapFEC.load(mapFECfile);


    if (mPrint) {
      try {
        int deId = 819;
        int dsIddet = 1134;
        const Segmentation& segment = segmentation(deId);
        for(int channel = 0; channel < 64; channel++) {
          int padId = segment.findPadByFEE(dsIddet, int(channel));
          std::cout << "DE# " << deId << "  DSid " << dsIddet << "  channel " << channel << "  PadId " << padId << std::endl;
        }
      } catch (const std::exception& e) {
        std::cout << "Failed to get padId: " << e.what() << std::endl;
      }
    }
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

      if (payloadSize == 0)
        continue;

      gsl::span<const std::byte> buffer(reinterpret_cast<const std::byte*>(raw), sizeof(o2::header::RAWDataHeaderV4) + payloadSize);
      decodeBuffer(buffer);
    }
  }

  //_________________________________________________________________________________________________
  void decodeReadout(const o2::framework::DataRef& input)
  {
    static int nFrame = 1;
    // get the input buffer
    if (input.spec->binding != "readout")
      return;

    const auto* header = o2::header::get<header::DataHeader*>(input.header);
    if (false)
      printf("Header: %p\n", (void*)header);
    if (!header)
      return;

    if (false)
      printf("payloadSize: %d\n", (int)header->payloadSize);
    if (false)
      printf("payload: %p\n", input.payload);

    auto const* raw = input.payload;
    // size of payload
    size_t payloadSize = header->payloadSize;

    if (mPrint)
      std::cout << nFrame << "  payloadSize=" << payloadSize << std::endl;
    nFrame += 1;
    if (payloadSize == 0)
      return;

    gsl::span<const std::byte> buffer(reinterpret_cast<const std::byte*>(raw), payloadSize);
    decodeBuffer(buffer);
  }

  //_________________________________________________________________________________________________
  void run(framework::ProcessingContext& pc)
  {
    currentBufId  = 1 - currentBufId;
    previousBufId = 1 - previousBufId;
    digitsBuffer[currentBufId].digits.clear();

    decodeTF(pc);
    for (auto&& input : pc.inputs()) {
      if (input.spec->binding == "readout")
        decodeReadout(input);
    }

    std::vector<o2::mch::Digit> digits;
    for(size_t i = 0; i < digitsBuffer[previousBufId].digits.size(); i++) {
      Digit& d = digitsBuffer[previousBufId].digits[i].digit;
      if(digitsBuffer[previousBufId].digits[i].merged)
        continue;
      digits.emplace_back(d);
    }

    if (mPrint) {
      std::cout << "previousBufId: " << previousBufId << "  digitsBuffer[previousBufId].digits.size(): " << digitsBuffer[previousBufId].digits.size()
          << "  digits.size(): " << digits.size() << std::endl;
      for (auto d : digits) {
        std::cout << "  DE# " << d.getDetID() << " PadId " << d.getPadID() << " ADC " << d.getADC() << " time " << d.getTimeStamp() << std::endl;
      }
    }

    const size_t OUT_SIZE = sizeof(o2::mch::Digit) * digits.size();

    // send the output buffer via DPL
    char* outbuffer = nullptr;
    outbuffer = (char*)realloc(outbuffer, OUT_SIZE);
    memcpy(outbuffer, digits.data(), OUT_SIZE);

    // create the output message
    auto freefct = [](void* data, void*) { free(data); };
    pc.outputs().adoptChunk(Output{"MCH", "DIGITS", 0}, outbuffer, OUT_SIZE, freefct, nullptr);
  }

 private:
  std::function<std::optional<DsDetId>(DsElecId)> mElec2Det;
  std::function<std::optional<uint16_t>(FeeLinkId)> mFee2LinkId{nullptr};
  size_t mNrdhs{0};

  std::ifstream mInputFile{}; ///< input file
  bool mDs2manu = false;      ///< print convert channel numbering from Run3 to Run1-2 order
  bool mPrint = false;        ///< print digits
  MapCRU mMapCRU;
  MapFEC mMapFEC;

  DigitsBuffer digitsBuffer[2];
  int currentBufId = {1};
  int previousBufId = {0};
};

//_________________________________________________________________________________________________
o2::framework::DataProcessorSpec getDecodingSpec()
{
  return DataProcessorSpec{
    "DataDecoder",
    //o2::framework::select("TF:MCH/RAWDATA, re:ROUT/RAWDATA"),
    o2::framework::select("readout:ROUT/RAWDATA"),
    Outputs{OutputSpec{"MCH", "DIGITS", 0, Lifetime::Timeframe}},
    AlgorithmSpec{adaptFromTask<DataDecoderTask>()},
    Options{{"print", VariantType::Bool, false, {"print digits"}},
            {"cru-map", VariantType::String, "", {"custom CRU mapping"}},
            {"fec-map", VariantType::String, "", {"custom FEC mapping"}},
            {"ds2manu", VariantType::Bool, false, {"convert channel numbering from Run3 to Run1-2 order"}}}};
}

} // end namespace raw
} // end namespace mch
} // end namespace o2
