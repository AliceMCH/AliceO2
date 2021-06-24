// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file ROFTimeClusterFinder.h
/// \brief Class to group the fired pads according to their time stamp
///
/// \author Andrea Ferrero, CEA

#ifndef ALICEO2_MCH_ROFTIMECLUSTERFINDER_H_
#define ALICEO2_MCH_ROFTIMECLUSTERFINDER_H_

#include <cassert>
#include <cstdint>
#include <vector>
#include <gsl/span>

#include "DataFormatsMCH/Digit.h"
#include "DataFormatsMCH/ROFRecord.h"

namespace o2
{
namespace mch
{

class ROFTimeClusterFinder
{
 public:
  using ROFVector = std::vector<o2::mch::ROFRecord>;

  ROFTimeClusterFinder(gsl::span<const o2::mch::ROFRecord> rof, uint32_t firstTForbit, int debug);
  ~ROFTimeClusterFinder() = default;

  void process();

  ROFVector getROFRecords() { return mOutputROFs; }

  char* saveROFRsToBuffer(size_t& bufSize);
  void dumpInputROFs();
  void dumpOutputROFs();

 private:
  static constexpr uint32_t sBCinOneMicrosec = 1000 / 25;
  static constexpr uint32_t sROFSearchWin = sBCinOneMicrosec * 2;
  static constexpr uint32_t sROFSearchSteps = 3;
  static constexpr uint32_t sROFWinSteps = sROFSearchSteps - 1;
  static constexpr uint32_t sROFWinSteps2 = sROFWinSteps * 2;
  static constexpr uint32_t sROFWinStepSize = sBCinOneMicrosec / sROFWinSteps;

  void findFirstValidROF(int64_t winStart, size_t& iRof);
  void fillStepArrays(int64_t winStart, size_t& iRof, std::vector<uint64_t>& stepEntries, std::vector<uint64_t>& rofIds);
  void findMax(int& stepMax, uint32_t& maxEntries, const std::vector<uint64_t>& stepEntries);
  void storeROFs(int stepMax, uint32_t maxEntries,
      std::vector<uint64_t>& stepEntries, std::vector<uint64_t>& rofIds);
  void advanceSearchWindow(int stepMax, int64_t& winStart, size_t& iRof, std::vector<uint64_t>& rofIds);

  gsl::span<const o2::mch::ROFRecord> mInputROFs;
  o2::InteractionRecord mTFstart;
  int mDebug;

  ROFVector mOutputROFs{};
};

} // namespace mch
} // namespace o2

#endif // ALICEO2_MCH_ROFTIMECLUSTERFINDER_H_
