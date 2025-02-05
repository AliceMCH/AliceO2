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

o2::mch::DigitModifier createST1MappingCorrector(int /*runNumber*/)
{
  return {};
}

/** initialization of the pad remapping table for Station 2 DEs
 */
void initST2PadsRemappingTable(std::unordered_map<int, std::unordered_map<int, int>>& padsRemapping)
{
  // Remapping of ST2 DS boards near the rounded part
  {
    std::array<int, 8> deToRemap{ 300, 301, 302, 303, 400, 401, 402, 403 };
    std::array<int, 5> dsToRemap{ 99, 100, 101, 102, 103 };

    for (auto deId : deToRemap) {

      const o2::mch::mapping::Segmentation& segment = o2::mch::mapping::segmentation(deId);
      for (auto dsId : dsToRemap) {
        // double loop on DS channels
        // 1. find the minimum pad index of the DS board
        int padIdMin = -1;
        int channelForPadIdMin = -1;
        for (int channel = 0; channel < 64; channel++) {
          auto padId = segment.findPadByFEE(dsId, int(channel));
          if (padId < 0) {
            // skip non-connected channels
            // this should never occur in this specific case, should we rise an exception instead?
            continue;
          }
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

o2::mch::DigitModifier createST2MappingCorrector(int runNumber)
{
  constexpr int lastRunToBeFixed = 560402;
  // ST2 mapping needs to be corrected only for data collected up to the end of 2024 Pb-Pb
  if (runNumber > lastRunToBeFixed) {
    // do not modify digits collected after 2024 Pb-Pb
    return {};
  }

  static std::unordered_map<int, std::unordered_map<int, int>> padsRemapping;

  if (padsRemapping.empty()) {
    initST2PadsRemappingTable(padsRemapping);
  }

  return [](o2::mch::Digit& digit) {
    // check if the current DE needs some remapping
    auto padsRemappingForDe = padsRemapping.find(digit.getDetID());
    if (padsRemappingForDe != padsRemapping.end()) {
      // check if the current padID needs to be remapped
      auto padIDRemapped = padsRemappingForDe->second.find(digit.getPadID());
      if (padIDRemapped != padsRemappingForDe->second.end()) {
        // update the digit
        digit.setPadID(padIDRemapped->second);
      }
    }
  };
}
} // namespace

namespace o2::mch
{
DigitModifier createDigitModifier(int runNumber,
                                  bool updateST1,
                                  bool updateST2)
{
  DigitModifier modifierST1 = updateST1 ? createST1MappingCorrector(runNumber) : DigitModifier{};
  DigitModifier modifierST2 = updateST2 ? createST2MappingCorrector(runNumber) : DigitModifier{};

  if (modifierST1 || modifierST2) {
    return [modifierST1, modifierST2](Digit& digit) {
      // the ST1/ST2 modifiers are mutually exclusive, depending on the DeID associated to the digit
      auto detID = digit.getDetID();
      if (modifierST1 && detID >= 100 && detID < 300) {
        modifierST1(digit);
      }
      if (modifierST2 && detID >= 300 && detID < 500) {
        modifierST2(digit);
      }
    };
  } else {
    // return an empty function if none of the modifiers is set
    return {};
  }
}

} // namespace o2::mch
