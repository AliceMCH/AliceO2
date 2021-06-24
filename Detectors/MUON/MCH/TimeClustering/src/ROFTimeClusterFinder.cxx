// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include "MCHTimeClustering/ROFTimeClusterFinder.h"

#include <iostream>
#include <fmt/format.h>

//#define ROFDEBUG 1

namespace o2
{
namespace mch
{

using namespace std;

//_________________________________________________________________________________________________

ROFTimeClusterFinder::ROFTimeClusterFinder(gsl::span<const o2::mch::ROFRecord> rofs, uint32_t firstTForbit, int debug) :
        mInputROFs(rofs), mTFstart(0, firstTForbit), mDebug(debug)
{
}

void ROFTimeClusterFinder::findFirstValidROF(int64_t winStart, size_t& iRof)
{
  // skip ROFs that belong to the previous TF
  for (iRof = 0; iRof < mInputROFs.size(); iRof++) {
    const auto& rof = mInputROFs[iRof];
    const auto ir = rof.getBCData();
    auto rofBc = ir.differenceInBC(mTFstart);

    if (mDebug) {
      std::cout << "  checking ROF #" << iRof << std::endl;
      std::cout << "    IR " << ir << std::endl;
      std::cout << "    TF " << mTFstart << std::endl;
      std::cout << "    diff " << rofBc << std::endl;
    }

    if (rofBc >= winStart) {
      break;
    }
    if (mDebug) {
      std::cout << "  skipping ROF #" << iRof << std::endl;
    }
  }
}

//_________________________________________________________________________________________________

void ROFTimeClusterFinder::fillStepArrays(int64_t winStart, size_t& iRof, std::vector<uint64_t>& stepEntries, std::vector<uint64_t>& rofIds)
{
  int64_t stepStart = winStart;
  for (uint32_t step = 0; step < sROFWinSteps2; step++) {
    size_t entries = 0;
    int64_t stepEnd = stepStart + sROFWinStepSize;
    rofIds[step] = iRof;
    if (mDebug) {
      std::cout << fmt::format("    step {}  start {}  end {}  iRof {}\n", step, stepStart, stepEnd, iRof);
    }

    // compute the cumulative number of digits
    for (; iRof < mInputROFs.size(); iRof++) {
      const auto& rof = mInputROFs[iRof];
      const auto& ir = rof.getBCData();
      auto rofBc = ir.differenceInBC(mTFstart);
      if (rofBc >= stepEnd) break;

      entries += rof.getNEntries();
      if (mDebug) {
        std::cout << fmt::format("      adding ROF #{} with {} digits\n", iRof, rof.getNEntries());
      }
    }
    stepEntries[step] = entries;
    if (mDebug) {
      std::cout << fmt::format("  step entries: {}\n", entries);
    }

    stepStart = stepEnd;
  }
}

//_________________________________________________________________________________________________

void ROFTimeClusterFinder::findMax(int& stepMax, uint32_t& maxEntries, const std::vector<uint64_t>& stepEntries)
{
  // search for the window with the maximum number of entries
  for (size_t i = 0; i < sROFSearchSteps; i++) {
    int64_t totEntries = 0;
    // get the total number of entries in a time window of 1 us starting at time step i
    for (int step = 0; step < sROFWinSteps; step++) {
      totEntries += stepEntries[step + i];
    }

    if (totEntries > maxEntries) {
      maxEntries = totEntries;
      stepMax = i;
    }
  }
  if (mDebug) {
    std::cout << fmt::format("  stepMax {}  entries {}\n", stepMax, maxEntries);
  }
}

//_________________________________________________________________________________________________

void ROFTimeClusterFinder::storeROFs(int stepMax, uint32_t maxEntries,
    std::vector<uint64_t>& stepEntries, std::vector<uint64_t>& rofIds)
{
  if (stepMax < 0) { return; }

  // 1. store a ROF corresponding to the time interval before the maximum
  uint32_t rofNEntries = 0;
  for (int64_t i = 0; i < stepMax; i++) {
    const auto& rof = mInputROFs[i];
    rofNEntries += stepEntries[i];
  }
  if (rofNEntries > 0) {
    const auto& firstRof = mInputROFs[rofIds[0]];
    mOutputROFs.emplace_back(firstRof.getBCData(), firstRof.getFirstIdx(), rofNEntries);
  }

  // 2. store a ROF corresponding to the time interval with the maximum, unless it is at the
  //    upper edge of the 2 us of the search range. In this case we continue the search of the maximum
  //    in the next iteration
  if (stepMax < (sROFSearchSteps - 1)) {
    const auto& firstRof = mInputROFs[rofIds[stepMax]];
    mOutputROFs.emplace_back(firstRof.getBCData(), firstRof.getFirstIdx(), maxEntries);
  }
}

//_________________________________________________________________________________________________

void ROFTimeClusterFinder::advanceSearchWindow(int stepMax, int64_t& winStart, size_t& iRof, std::vector<uint64_t>& rofIds)
{
  if (stepMax < 0) {
    winStart += sROFSearchWin;
  } else {
    // 2. store a ROF corresponding to the time interval with the maximum, unless it is at the
    //    upper edge of the 2 us of the search range. In this case we continue the search of the maximum
    //    in the next iteration
    if (stepMax < (sROFSearchSteps - 1)) {
      winStart += (stepMax + sROFWinSteps) * sROFWinStepSize;
      iRof = rofIds[stepMax + sROFWinSteps];
    } else {
      // the 1 us region with the maximum number of entries was at the upper edge of the search window
      // in this case we do not have a local maximum, and we resume the search from the beginning of this window
      winStart += stepMax * sROFWinStepSize;
      iRof = rofIds[stepMax];
    }
  }
}

//_________________________________________________________________________________________________

void ROFTimeClusterFinder::process()
{
  // sliding window algorithm
  // the algorithm counts the number of digits in a window of 1 us and searches for a local maximum
  // at each iteration the sliding 1 us window spans a time interval of 2 us and keeps the track of the
  // position that maximizes the number of digits
  int64_t winStart{0};
  size_t iRof;

  if (mDebug) {
    std::cout << fmt::format("sBCinOneMicrosec {}  sROFWinSteps {}  sROFWinStepSize {}\n", sBCinOneMicrosec, sROFWinSteps, sROFWinStepSize);
  }

  std::vector<uint64_t> stepEntries(sROFWinSteps2);
  std::vector<uint64_t> rofIds(sROFWinSteps2);

  findFirstValidROF(winStart, iRof);

  while (true) {

    if (iRof == mInputROFs.size()) { break; }

    if (mDebug) {
      std::cout << fmt::format("\n  winStart {}\n", winStart);
    }

    // fill the arrays with the first ROF index and the total number of entries for each search step
    fillStepArrays(winStart, iRof, stepEntries, rofIds);

    // search for the step that gives the maximum number of entries
    int stepMax{-1};
    uint32_t maxEntries{0};
    findMax(stepMax, maxEntries, stepEntries);

    // store the output ROFs
    storeROFs(stepMax, maxEntries, stepEntries, rofIds);

    // move the search window after the last stored ROF
    advanceSearchWindow(stepMax, winStart, iRof, rofIds);
  }

  if (mDebug) {
    std::cout << "ROF processing completed" << std::endl;
  }
}

//_________________________________________________________________________________________________

char* ROFTimeClusterFinder::saveROFRsToBuffer(size_t& bufSize)
{
  static constexpr size_t sizeOfROFRecord = sizeof(o2::mch::ROFRecord);

#ifdef ROFDEBUG
  dumpOutputROFs();
#endif

  bufSize = mOutputROFs.size() * sizeOfROFRecord;
  o2::mch::ROFRecord* buf = reinterpret_cast<o2::mch::ROFRecord*>(malloc(bufSize));
  if (!buf) {
    bufSize = 0;
    return nullptr;
  }

  o2::mch::ROFRecord* p = buf;
  for (size_t i = 0; i < mOutputROFs.size(); i++) {
    auto& rof = mOutputROFs[i];
    memcpy(p, &(rof), sizeOfROFRecord);
    p += 1;
  }

  return reinterpret_cast<char*>(buf);
}

//_________________________________________________________________________________________________

void ROFTimeClusterFinder::dumpInputROFs()
{
  for (size_t i = 0; i < mInputROFs.size(); i++) {
    auto& rof = mInputROFs[i];
    const auto ir = rof.getBCData();
    auto rofBc = ir.differenceInBC(mTFstart);
    std::cout << fmt::format("    ROF {}  RANGE {}-{}  SIZE {}  IR {}-{},{}  DIFF {}\n", i, rof.getFirstIdx(), rof.getLastIdx(),
        rof.getNEntries(), rof.getBCData().orbit, rof.getBCData().bc, rof.getBCData().toLong(), rofBc);
  }
}

//_________________________________________________________________________________________________

void ROFTimeClusterFinder::dumpOutputROFs()
{
  for (size_t i = 0; i < mOutputROFs.size(); i++) {
    auto& rof = mOutputROFs[i];
    std::cout << fmt::format("    ROF {}  RANGE {}-{}  SIZE {}  IR {}-{},{}\n", i, rof.getFirstIdx(), rof.getLastIdx(),
        rof.getNEntries(), rof.getBCData().orbit, rof.getBCData().bc, rof.getBCData().toLong());
  }
}

} // namespace mch
} // namespace o2
