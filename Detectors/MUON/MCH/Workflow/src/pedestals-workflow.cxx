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
#include <array>
#include <functional>

#include "Framework/CallbackService.h"
#include "Framework/ConfigParamRegistry.h"
#include "Framework/ControlService.h"
#include "Framework/DataProcessorSpec.h"
#include "Framework/Lifetime.h"
#include "Framework/Output.h"
#include "Framework/Task.h"
#include "Framework/WorkflowSpec.h"
#include "Framework/ConfigParamSpec.h"
#include "Framework/CompletionPolicyHelpers.h"

#include "Headers/RDHAny.h"
#include "MCHRawDecoder/PageDecoder.h"

#include "Headers/RAWDataHeader.h"
#include "DetectorsRaw/RDHUtils.h"
#include "DPLUtils/DPLRawParser.h"

#include "MCHRawCommon/DataFormats.h"
#include "MCHRawElecMap/Mapper.h"
#include "MCHMappingInterface/Segmentation.h"

static const size_t SOLAR_ID_MAX = 100 * 8;

namespace o2
{
namespace mch
{
namespace raw
{

using namespace o2;
using namespace o2::framework;
using namespace o2::mch::mapping;
using RDH = o2::header::RDHAny;

static std::string readFileContent(std::string& filename)
{
  std::string content;
  std::string s;
  std::ifstream in(filename);
  while (std::getline(in, s)) {
    content += s;
    content += "\n";
  }
  std::cout << "readFileContent(" << filename << "):" << std::endl
            << content << std::endl;
  return content;
};

static void patchPage(gsl::span<const std::byte> rdhBuffer, bool verbose)
{
  static int mNrdhs = 0;
  auto& rdhAny = *reinterpret_cast<RDH*>(const_cast<std::byte*>(&(rdhBuffer[0])));
  mNrdhs++;

  auto cruId = o2::raw::RDHUtils::getCRUID(rdhAny) & 0xFF;
  auto flags = o2::raw::RDHUtils::getCRUID(rdhAny) & 0xFF00;
  auto endpoint = o2::raw::RDHUtils::getEndPointID(rdhAny);
  uint32_t feeId = cruId * 2 + endpoint + flags;
  o2::raw::RDHUtils::setFEEID(rdhAny, feeId);

  if (verbose) {
    std::cout << mNrdhs << "--\n";
    o2::raw::RDHUtils::printRDH(rdhAny);
  }
};

static bool isValidDeID(int deId)
{
  for (auto id : deIdsForAllMCH) {
    if (id == deId) {
      return true;
    }
  }

  return false;
}

//=======================
// Data decoder
class PedestalsTask
{
 public:
  //PedestalsTask() : mInputSpec() {}
  //PedestalsTask(PedestalsTask& t) : mInputSpec() {}
  PedestalsTask(std::string spec) : mInputSpec(spec) {}

  void initElec2DetMapper(std::string filename)
  {
    std::cout << "[initElec2DetMapper] filename=" << filename << std::endl;
    if (filename.empty()) {
      mElec2Det = createElec2DetMapper<ElectronicMapperGenerated>();
    } else {
      ElectronicMapperString::sFecMap = readFileContent(filename);
      mElec2Det = createElec2DetMapper<ElectronicMapperString>();
    }
  };

  void initFee2SolarMapper(std::string filename)
  {
    std::cout << "[initFee2SolarMapper] filename=" << filename << std::endl;
    if (filename.empty()) {
      mFee2Solar = createFeeLink2SolarMapper<ElectronicMapperGenerated>();
    } else {
      ElectronicMapperString::sCruMap = readFileContent(filename);
      mFee2Solar = createFeeLink2SolarMapper<ElectronicMapperString>();
    }
  };

  //_________________________________________________________________________________________________
  void init(framework::InitContext& ic)
  {
    mDebug = ic.options().get<bool>("debug");

    auto mapCRUfile = ic.options().get<std::string>("cru-map");
    auto mapFECfile = ic.options().get<std::string>("fec-map");
    initFee2SolarMapper(mMapCRUfile);
    initElec2DetMapper(mMapFECfile);

    mNoiseThreshold = ic.options().get<float>("noise-threshold");
    mPedestalThreshold = ic.options().get<float>("pedestal-threshold");

    ic.services().get<CallbackService>().set(CallbackService::Id::Stop, [this]() { stop(); });
    ic.services().get<CallbackService>().set(CallbackService::Id::Reset, [this]() { reset(); });

    for (int s = 0; s <= SOLAR_ID_MAX; s++) {
      for (int i = 0; i < 40; i++) {
        for (int j = 0; j < 64; j++) {
          nhits[s][i][j] = 0;
          pedestal[s][i][j] = noise[s][i][j] = 0;
        }
      }
    }
  }

  //_________________________________________________________________________________________________
  void reset()
  {

  }

  //_________________________________________________________________________________________________
  void stop()
  {
    //if (mDebug) {
      std::cout << "\n\n============================\nStop called\n";
    //}
    for (int s = 0; s <= SOLAR_ID_MAX; s++) {
      for (int i = 0; i < 40; i++) {
        for (int j = 0; j < 64; j++) {
          //std::cout << "SOLAR " << s << "  DS " << i << "  CH " << j << "  nhits " << nhits[s][i][j] << std::endl;
          if (nhits[s][i][j] == 0) continue;

          bool ok = true;
          if (pedestal[s][i][j] > mPedestalThreshold) {
            //std::cout << "SOLAR " << s << "  DS " << i << "  CH " << j << "  excluded, pedestal=" << pedestal[s][i][j] << std::endl;
            ok = false;
          }

          double rms = std::sqrt(noise[s][i][j] / nhits[s][i][j]);
          if (rms > mNoiseThreshold) {
            //std::cout << "SOLAR " << s << "  DS " << i << "  CH " << j << "  excluded, rms=" << rms << std::endl;
            ok = false;
          }
          if (!ok) {
            std::cout << "SOLAR " << s << "  DS " << i << "  CH " << j << "  nhits " << nhits[s][i][j]
                << "  pedestal " << pedestal[s][i][j] << "  RMS " << rms << "  OK " << ok << std::endl;
          }
        }
      }
    }
  }

  //_________________________________________________________________________________________________
  void decodePage(gsl::span<const std::byte> page)
  {
    size_t ndigits{0};

    uint32_t orbit;

    auto channelHandler = [&](DsElecId dsElecId, uint8_t channel, o2::mch::raw::SampaCluster sc) {

      auto solarId = dsElecId.solarId();
      auto dsId = dsElecId.elinkId();

      for (auto s: sc.samples) {
        nhits[solarId][dsId][channel] += 1;
        uint64_t N = nhits[solarId][dsId][channel];

        double p0 = pedestal[solarId][dsId][channel];
        double p = p0 + (s - p0) / N;
        pedestal[solarId][dsId][channel] = p;

        double M0 = noise[solarId][dsId][channel];
        double M = M0 + (s - p0) * (s - p);
        noise[solarId][dsId][channel] = M;
      }
      if (mDebug) {
        std::cout << "solarId " << (int)solarId << "  dsId " << (int)dsId << "  ch " << (int)channel << "  nsamples " << sc.samples.size()
              << "  nhits "<< nhits[solarId][dsId][channel] << "  ped "<< pedestal[solarId][dsId][channel] << "  noise " << noise[solarId][dsId][channel] << std::endl;
      }
      ++ndigits;
    };

    patchPage(page, mDebug);

    if (!mDecoder) {
      DecodedDataHandlers handlers;
      handlers.sampaChannelHandler = channelHandler;
      mDecoder = mFee2Solar ? o2::mch::raw::createPageDecoder(page, handlers, mFee2Solar)
                            : o2::mch::raw::createPageDecoder(page, handlers);
    }
    mDecoder(page);
  }

  //_________________________________________________________________________________________________
  void decodeBuffer(gsl::span<const std::byte> buf)
  {
    if (mDebug) {
      std::cout << "\n\n============================\nStart of new buffer\n";
    }
    size_t bufSize = buf.size();
    size_t pageStart = 0;
    while (bufSize > pageStart) {
      RDH* rdh = reinterpret_cast<RDH*>(const_cast<std::byte*>(&(buf[pageStart])));
      auto rdhVersion = o2::raw::RDHUtils::getVersion(rdh);
      auto rdhHeaderSize = o2::raw::RDHUtils::getHeaderSize(rdh);
      if (rdhHeaderSize != 64) {
        break;
      }
      auto pageSize = o2::raw::RDHUtils::getOffsetToNext(rdh);

      gsl::span<const std::byte> page(reinterpret_cast<const std::byte*>(rdh), pageSize);
      decodePage(page);

      pageStart += pageSize;
    }
  }

  //_________________________________________________________________________________________________
  // the decodeTF() function processes the messages generated by the (sub)TimeFrame builder
  void decodeTF(framework::ProcessingContext& pc)
  {
    // get the input buffer
    auto& inputs = pc.inputs();
    DPLRawParser parser(inputs, o2::framework::select(mInputSpec.c_str()));

    for (auto it = parser.begin(), end = parser.end(); it != end; ++it) {
      auto const* raw = it.raw();
      if (!raw) {
        continue;
      }
      size_t payloadSize = it.size();

      gsl::span<const std::byte> buffer(reinterpret_cast<const std::byte*>(raw), sizeof(RDH) + payloadSize);
      decodeBuffer(buffer);
    }
  }

  //_________________________________________________________________________________________________
  // the decodeReadout() function processes the messages generated by o2-mch-cru-page-reader-workflow
  void decodeReadout(const o2::framework::DataRef& input)
  {
    const auto* header = o2::header::get<header::DataHeader*>(input.header);
    if (!header) {
      return;
    }

    auto const* raw = input.payload;
    size_t payloadSize = header->payloadSize;

    gsl::span<const std::byte> buffer(reinterpret_cast<const std::byte*>(raw), payloadSize);
    decodeBuffer(buffer);
  }

  //_________________________________________________________________________________________________
  void run(framework::ProcessingContext& pc)
  {
    reset();
    for (auto&& input : pc.inputs()) {
      if (input.spec->binding == "TF") {
        decodeTF(pc);
      }
      if (input.spec->binding == "readout") {
        decodeReadout(input);
      }
    }
  }

 private:
  o2::mch::raw::PageDecoder mDecoder;
  SampaChannelHandler mChannelHandler;

  Elec2DetMapper mElec2Det{nullptr};
  FeeLink2SolarMapper mFee2Solar{nullptr};
  std::string mMapCRUfile;
  std::string mMapFECfile;

  uint64_t nhits[SOLAR_ID_MAX+1][40][64];
  double pedestal[SOLAR_ID_MAX+1][40][64];
  double noise[SOLAR_ID_MAX+1][40][64];

  float mNoiseThreshold;
  float mPedestalThreshold;

  std::string mInputSpec;     /// selection string for the input data
  bool mDebug = {false};      /// flag to enable verbose output
};

} // namespace raw
} // namespace mch
} // end namespace o2


using namespace o2::framework;

void customize(std::vector<ConfigParamSpec>& workflowOptions)
{
  //workflowOptions.push_back(ConfigParamSpec{"dataspec", VariantType::String, "TF:MCH/RAWDATA", {"selection string for the input data"}});
  workflowOptions.push_back(ConfigParamSpec{"dataspec", VariantType::String, "readout:ROUT/RAWDATA", {"selection string for the input data"}});
}

#include "Framework/runDataProcessing.h"

using namespace o2;
using namespace o2::framework;

//_________________________________________________________________________________________________
o2::framework::DataProcessorSpec getPedestalsSpec(std::string inputSpec)
{
  //o2::mch::raw::PedestalsTask task();
  return DataProcessorSpec{
    "Pedestals",
    o2::framework::select(inputSpec.c_str()),
    Outputs{},
    AlgorithmSpec{adaptFromTask<o2::mch::raw::PedestalsTask>(inputSpec)},
    Options{{"debug", VariantType::Bool, false, {"enable verbose output"}},
            {"noise-threshold", VariantType::Float, (float)2.0, {"maximum acceptable noise value"}},
            {"pedestal-threshold", VariantType::Float, (float)150, {"maximum acceptable pedestal value"}},
            {"cru-map", VariantType::String, "", {"custom CRU mapping"}},
            {"fec-map", VariantType::String, "", {"custom FEC mapping"}}}};
}

WorkflowSpec defineDataProcessing(const ConfigContext& config)
{
  auto inputSpec = config.options().get<std::string>("dataspec");

  WorkflowSpec specs;

  DataProcessorSpec producer = getPedestalsSpec("readout:ROUT/RAWDATA");
  specs.push_back(producer);

  return specs;
}
