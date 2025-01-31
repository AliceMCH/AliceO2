// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include "MCHDigitFiltering/DigitModifier.h"

#include "DataFormatsMCH/Digit.h"
#include "MCHRawElecMap/Mapper.h"
#include "MCHMappingInterface/Segmentation.h"
#include <functional>
#include <array>
#include <unordered_map>

namespace
{
/** initialization of the pad remapping table for Station 2 DEs
 */
void initST2PadsRemappingTable(std::unordered_map<int, std::unordered_map<int, int>>& padsRemapping)
{
  // Remapping of ST2 DS boards near the rounded part
  {
    std::array<int, 8> deToRemap{300, 301, 302, 303, 400, 401, 402, 403};
    std::array<int, 5> dsToRemap{99, 100, 101, 102, 103};

    for (auto deId : deToRemap) {

      const o2::mch::mapping::Segmentation& segment = o2::mch::mapping::segmentation(deId);
      for (auto dsId : dsToRemap) {
        // double loop on DS channels
        // 1. find the minimum pad index of the DS board
        int padIdMin = -1;
        int channelForPadIdMin = -1;
        for (int channel = 0; channel < 64; channel++) {
          auto padId = segment.findPadByFEE(dsId, int(channel));
          if (padIdMin < 0 || padId < padIdMin) {
            padIdMin = padId;
            channelForPadIdMin = channel;
          }
        }

        // 2. build the re-mapping table
        for (int channel = 0; channel < 64; channel++) {
          auto padId = segment.findPadByFEE(dsId, int(channel));
          if (padId < padIdMin) {
            // something is wrong here...
            continue;
          }
          int padIdInDS = padId - padIdMin;
          int padColumn = padIdInDS / 16;
          int padRow = padIdInDS % 16;

          int padIdRemapped = -1;

          switch (padColumn) {
            case 0:
              // shift right by 3 columns
              padIdRemapped = padId + 16 * 3;
              break;
            case 1:
              // shift right by 1 column
              padIdRemapped = padId + 16;
              break;
            case 2:
              // shift left by 1 column
              padIdRemapped = padId - 16;
              break;
            case 3:
              // shift left by 3 columns
              padIdRemapped = padId - 16 * 3;
              break;
          }

          padsRemapping[deId][padId] = padIdRemapped;
        }
      }
    }
  }
}

o2::mch::DigitModifier createST1MappingCorrector(int /*runNumber*/)
{
  return [](o2::mch::Digit& digit) {
    return;
  };
}
} // namespace

o2::mch::DigitModifier createST2MappingCorrector(int /*runNumber*/)
{
  static std::unordered_map<int, std::unordered_map<int, int>> padsRemapping;

  if (padsRemapping.empty()) {
    initST2PadsRemappingTable(padsRemapping);
  }

  return [](o2::mch::Digit& digit) {
    // Only consider DEs from ST2
    if (digit.getDetID() >= 300 && digit.getDetID() < 500) {
      // check if the current DE needs some remapping
      if (padsRemapping.count(digit.getDetID()) > 0) {
        // check if the current padID needs to be remapped
        auto& padsRemappingForDe = padsRemapping[digit.getDetID()];
        if (padsRemappingForDe.count(digit.getPadID()) > 0) {
          // get the corrected padID
          int padIDRemapped = padsRemappingForDe[digit.getPadID()];
          // update the digit
          digit.setPadID(padIDRemapped);
        }
      }
    }
  };
}

namespace o2::mch
{
DigitModifier createDigitModifier(int runNumber,
                                  bool correctST1Mapping,
                                  bool correctST2Mapping)
{
  std::vector<DigitModifier> parts;

  if (correctST1Mapping) {
    parts.emplace_back(createST1MappingCorrector(runNumber));
  }
  if (correctST2Mapping) {
    parts.emplace_back(createST2MappingCorrector(runNumber));
  }
  return [parts](Digit& digit) {
    for (const auto& p : parts) {
      p(digit);
    }
  };
}

} // namespace o2::mch
