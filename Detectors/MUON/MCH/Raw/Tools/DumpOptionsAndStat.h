#ifndef O2_MCH_RAW_DUMPOPTIONSANDSTAT_H
#define O2_MCH_RAW_DUMPOPTIONSANDSTAT_H

#include "DumpBuffer.h"
#include "Headers/RAWDataHeader.h"
#include "MCHRawCommon/DataFormats.h"
#include "MCHRawDecoder/Decoder.h"
#include "MCHRawElecMap/Mapper.h"
#include "MCHMappingInterface/Segmentation.h"
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

class DumpOptions
{
 public:
  DumpOptions(unsigned int deId, unsigned int maxNofRDHs, bool showRDHs, bool jsonOutput)
    : mDeId{deId}, mMaxNofRDHs{maxNofRDHs == 0 ? std::numeric_limits<unsigned int>::max() : maxNofRDHs}, mShowRDHs{showRDHs}, mJSON{jsonOutput} {}

  unsigned int deId() const
  {
    return mDeId;
  }
  unsigned int maxNofRDHs() const
  {
    return mMaxNofRDHs;
  }

  bool showRDHs() const
  {
    return mShowRDHs;
  }

  bool json() const
  {
    return mJSON;
  }

  std::optional<uint16_t> cruId() const
  {
    return mCruId;
  }

  void cruId(uint16_t c) { mCruId = c; }

 private:
  unsigned int mDeId;
  unsigned int mMaxNofRDHs;
  bool mShowRDHs;
  bool mJSON;
  std::optional<uint16_t> mCruId{std::nullopt};
};

struct Stat {
  double adc{0};
  double rms{0};
  double q{0};
  int n{0};
  void incr(int v)
   {
    n++;
       auto newAdc = adc + v;
    adc = newAdc;
  }
};
    
#endif
