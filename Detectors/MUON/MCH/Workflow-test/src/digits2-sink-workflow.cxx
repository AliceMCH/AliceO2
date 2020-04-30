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
/// \file    digits-sink-workflow.cxx
/// \author  Andrea Ferrero
///
/// \brief This is an executable that dumps to a file on disk the digits received via DPL.
///
/// This is an executable that dumps to a file on disk the digits received via the Data Processing Layer.
/// It can be used to debug the raw decoding step. For example, one can do:
/// \code{.sh}
/// o2-mch-file-to-digits-workflow --infile=some_data_file | o2-mch-digits-sink-workflow --outfile digits.txt
/// \endcode
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
#include "Framework/DataProcessorSpec.h"
#include "Framework/runDataProcessing.h"

#include "DPLUtils/DPLRawParser.h"
#include "MCHBase/Digit.h"
#include "MCHMappingFactory/CreateSegmentation.h"

using namespace o2;
using namespace o2::framework;

namespace o2
{
namespace mch
{

// \class Digit2
/// \brief MCH digit implementation, with pad coordinates
class Digit2
{
 public:
  Digit2() = default;

  ~Digit2() = default;

 public:
  uint64_t mTime;
  uint32_t mDetID;
  uint32_t mPadID;         /// PadIndex to which the digit corresponds to
  uint64_t mADC; /// Amplitude of signal
  uint8_t mCathode;
  float mX, mY, mSizeX, mSizeY;

  ClassDefNV(Digit2, 1);
}; //class Digit2


namespace raw
{

using namespace o2;
using namespace o2::framework;

class DigitsSinkTask
{
 public:
  //_________________________________________________________________________________________________
  void init(framework::InitContext& ic)
  {
    /// Get the input file and other options from the context
    LOG(INFO) << "initializing digits2 sink";

    auto outputFileName = ic.options().get<std::string>("outfile");
    mOutputFile.open(outputFileName, std::ios::out | std::ios::binary);
    if (!mOutputFile.is_open()) {
      throw std::invalid_argument("Cannot open output file" + outputFileName);
    }

    auto stop = [this]() {
      /// close the input file
      LOG(INFO) << "stop file reader";
      this->mOutputFile.close();
    };
    ic.services().get<CallbackService>().set(CallbackService::Id::Stop, stop);
  }

  //_________________________________________________________________________________________________
  void run(framework::ProcessingContext& pc)
  {
    // get the input digits
    auto digits = pc.inputs().get<gsl::span<Digit>>("digits");

    std::vector<Digit2> digits2;
    for (auto d : digits) {
      Digit2 d2;
      d2.mTime = static_cast<uint64_t>(d.getTimeStamp());
      d2.mDetID = static_cast<uint32_t>(d.getDetID());
      d2.mPadID = static_cast<uint32_t>(d.getPadID());
      d2.mADC = static_cast<uint64_t>(d.getADC());
      const mapping::Segmentation& segment = mapping::segmentation(d.getDetID());

      d2.mCathode = (segment.isBendingPad(d.getPadID()) ? 0 : 1);
      d2.mX = segment.padPositionX(d.getPadID());
      d2.mY = segment.padPositionY(d.getPadID());
      d2.mSizeX = segment.padSizeX(d.getPadID());
      d2.mSizeY = segment.padSizeY(d.getPadID());

      digits2.push_back(d2);
    }
    gsl::span<Digit2> sdigits2 = digits2;
    float xtrk = 0, ytrk = 0;
    mOutputFile.write(reinterpret_cast<char*>(&xtrk), sizeof(float));
    mOutputFile.write(reinterpret_cast<char*>(&ytrk), sizeof(float));

    int nDigits = sdigits2.length();
    mOutputFile.write(reinterpret_cast<char*>(&nDigits), sizeof(int));
    mOutputFile.write(reinterpret_cast<const char*>(sdigits2.data()), sdigits2.size_bytes());

    int nClus = 0;
    mOutputFile.write(reinterpret_cast<char*>(&nClus), sizeof(int));
  }

 private:
  std::ofstream mOutputFile{}; ///< output file
};

} // end namespace raw
} // end namespace mch
} // end namespace o2

// clang-format off
WorkflowSpec defineDataProcessing(const ConfigContext&)
{
  WorkflowSpec specs;

  // The producer to generate some data in the workflow
  DataProcessorSpec producer{
    "Digits2Sink",
    Inputs{InputSpec{"digits", "MCH", "DIGITS", 0, Lifetime::Timeframe}},
    Outputs{},
    AlgorithmSpec{adaptFromTask<o2::mch::raw::DigitsSinkTask>()},
    Options{ { "outfile", VariantType::String, "digits.out", { "output file name" } },
      {"txt", VariantType::Bool, false, {"output digits in text format"}}}
  };
  specs.push_back(producer);

  return specs;
}
// clang-format on
