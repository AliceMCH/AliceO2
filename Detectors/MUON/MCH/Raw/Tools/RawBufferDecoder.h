// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#ifndef O2_MCH_RAW_HANDLERS_H
#define O2_MCH_RAW_HANDLERS_H

#include "DumpBuffer.h"
#include "Headers/RAWDataHeader.h"
#include "MCHRawCommon/DataFormats.h"
#include "MCHRawDecoder/Decoder.h"
#include "MCHRawElecMap/Mapper.h"
#include "MCHMappingFactory/CreateSegmentation.h"
#include "MCHBase/DigitBlock.h"
#include "MCHBase/Digit.h"
#include "boost/program_options.hpp"
#include <chrono>
#include <vector>
#include <fmt/format.h>
#include <fstream>
#include <gsl/span>
#include <iostream>
#include <rapidjson/document.h>
#include <rapidjson/stringbuffer.h>
#include <rapidjson/writer.h>
#include <optional>
#include <cstdint>
#include "DumpOptionsAndStat.h"


std::array<int, 64> refManu2ds_st345 = {
    63, 62, 61, 60, 59, 57, 56, 53, 51, 50, 47, 45, 44, 41, 38, 35,
    36, 33, 34, 37, 32, 39, 40, 42, 43, 46, 48, 49, 52, 54, 55, 58,
    7, 8, 5, 2, 6, 1, 3, 0, 4, 9, 10, 15, 17, 18, 22, 25,
    31, 30, 29, 28, 27, 26, 24, 23, 20, 21, 16, 19, 12, 14, 11, 13};

int manu2ds(int i){
  return refManu2ds_st345[i];
}

extern std::ostream& operator<<(std::ostream&, const o2::header::RAWDataHeaderV4&);

using namespace o2::framework;
using namespace o2::mch::mapping;
using namespace o2::mch::raw;
using RDHv4 = o2::header::RAWDataHeaderV4;

std::ostream& operator<<(std::ostream& os, const Stat& s)
{
  os << fmt::format("MEAN {:7.3f} NSAMPLES {:5d} ", s.adc, s.n);
  return os;
}


template <typename CHARGESUM, typename RDH>
class RawBufferDecoder
{
  std::function<std::optional<DsDetId>(DsElecId)> Elec2Det;
  std::function<std::optional<uint16_t>(CruLinkId)> cruLink2solar;
  bool verbose;
public:

  RawBufferDecoder(bool v=false): verbose(v)
  {
    std::vector<int> deidspan;
    for(int i=0; i<deIdsForAllMCH.size(); i++){
      deidspan.push_back(deIdsForAllMCH[i]);
    }

    Elec2Det = createElec2DetMapper<ElectronicMapperGenerated>(deidspan);
    cruLink2solar = o2::mch::raw::createCruLink2SolarMapper<ElectronicMapperGenerated>();
  }

  void decodeBuffer(gsl::span<uint8_t> sbuffer, std::vector<o2::mch::Digit> &digits)
  {
    //bool verbose = false;

    size_t ndigits{0};

    //std::vector< std::unique_ptr<o2::mch::Digit> > digits;
    //if(verbose) std::cout << "On nettoie le vector digits" << std::endl;


    auto channelHandler = [&](DsCruId dsId, uint8_t channel, o2::mch::raw::SampaCluster sc) {
      auto s = asString(dsId);
      channel = manu2ds(int(channel));
      if(verbose) {
        auto ch = fmt::format("{}-CH{}", s, channel);
        std::cout << ch << std::endl;
      }
      double digitadc(0);
      for (auto d = 0; d < sc.nofSamples(); d++) {
        digitadc += sc.samples[d];
      }


      int deId;
      int dsIddet;
      auto cruId = o2::mch::raw::decodeCruLinkId(dsId.linkId());
      if(auto solar = cruLink2solar(cruId); solar.has_value()) {
        DsElecId dsElecId(solar.value(), dsId.elinkGroupId(), dsId.elinkIndexInGroup());
        //Attention on est dans une boucle
        if(auto opt = Elec2Det(dsElecId); opt.has_value()) {
          DsDetId dsDetId = opt.value();
          dsIddet = dsDetId.dsId();
          deId = dsDetId.deId();
        }
      } else {
        dsIddet = 9999;
        deId = 819;
      }

      int padId = -1;
      try {
        const Segmentation& segment = segmentation(deId);
        //Segmentation segment(deId);

        padId = segment.findPadByFEE(dsIddet, int(channel));
        if(verbose)
          std::cout << "DS "<<(int)dsId.elinkId()<<"  CHIP "<<((int)channel)/32<<"  CH "<<((int)channel)%32<<"  ADC " << digitadc << "  DE# " << deId << "  DSid " << dsIddet << "  PadId " << padId << std::endl;
      } catch (const std::exception& e) { return; }


      int time = 0;

      digits.emplace_back( o2::mch::Digit() );
      o2::mch::Digit& mchdigit = digits.back();
      mchdigit.setDetID(deId);
      mchdigit.setPadID(padId);
      mchdigit.setADC(digitadc);
      mchdigit.setTimeStamp(time);

      if(verbose) std::cout << "DIGIT STORED:\nADC " << digits.back().getADC() << " DE# " << digits.back().getDetID() << " PadId " << digits.back().getPadID() << " time "<< digits.back().getTimeStamp() << std::endl;

      // std::cout << "For this digit we obtained a padId of " << padId << std::endl;
      ++ndigits;
    };


    size_t nrdhs{0};

    auto rdhHandler = [&](const RDH& rdh) -> std::optional<RDH> {
      nrdhs++;
      auto r = rdh;
      auto cruId = r.cruID;
      auto linkId = rdhLinkId(r);
      //auto solar = cruLink2solar(o2::mch::raw::CruLinkId(cruId, linkId));
      //if (!solar.has_value()) {
      //  std::cout << fmt::format("ERROR - Could not get solarUID from CRU,LINK=({},{})\n",
      //      cruId, linkId);
      //  return std::nullopt;
      //}
      r.feeId = cruId; //solar.value();
      if(verbose)
        std::cout << "RDH INFO: CRUID " << cruId << " LINKID " << int(linkId) << " ORBIT "<< rdhOrbit(rdh) << std::endl;
      return r;
    };

    o2::mch::raw::Decoder decode = o2::mch::raw::createDecoder<CHARGESUM, RDH>(rdhHandler, channelHandler);

    std::vector<std::chrono::microseconds> timers;

    DecoderStat decStat;

    decStat = decode(sbuffer);

    /*
    if(verbose) std::cout << "Filling buffer of digits ["<<ndigits<<"]..." << std::endl;

    o2::mch::Digit* digitsBuffer = NULL;
    digitsBuffer = (o2::mch::Digit*)realloc(digitsBuffer, sizeof(o2::mch::Digit) * ndigits);

    o2::mch::Digit* ptr = (o2::mch::Digit*)digitsBuffer;
    for(unsigned int di = 0; di < ndigits; di++) {

      memcpy(ptr, digits[di].get(), sizeof(o2::mch::Digit));
      if(verbose)
        if(di % 1 == 0){
          std::cout << "Added digit # " << di << " to buffer " << std::endl;
          std::cout<< *(digits[di].get())<<std::endl;;
        }
      ptr += 1;
      outsize += sizeof(o2::mch::Digit);
    }

    return (char*)digitsBuffer;
    */
  }
};




#endif
