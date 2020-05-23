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
/// \file    DatDecoder.cxx
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
#include "Headers/RDHAny.h"
#include "MCHRawCommon/DataFormats.h"
#include "MCHRawDecoder/PageDecoder.h"
#include "MCHRawElecMap/Mapper.h"
//#include "MCHRawCommon/RDHManip.h"
#include "MCHMappingInterface/Segmentation.h"
#include "MCHWorkflow/DataDecoderSpec.h"
//#include <array>

#define MCH_MERGER_FEEID_MAX 63

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

class MergerBase;

//_________________________________________________________________
//
// Data decoder
//_________________________________________________________________
class DataDecoder
{
 private:
  void initElec2DetMapper(std::string filename);
  void initFee2SolarMapper(std::string filename);

 public:
  void setPrint(bool print)
  {
    mPrint = print;
  }

  void setDs2manu(bool ds2manu)
  {
    mDs2manu = ds2manu;
  }

  void setMapCRUfile(const std::string& mapCRUfile)
  {
    mMapCRUfile = mapCRUfile;
  }

  void setMapFECfile(const std::string& mapFECfile)
  {
    mMapFECfile = mapFECfile;
  }

  void setChannelHandler(SampaChannelHandler channelHandler)
  {
    mChannelHandler = channelHandler;
  }

  void setRdhHandler(std::function<void(o2::header::RDHAny*)> rdhHandler)
  {
    mRdhHandler = rdhHandler;
  }

  void init();
  void reset();
  void decodeBuffer(gsl::span<const std::byte> page);

  std::vector<o2::mch::Digit>& getOutputDigits() { return mOutputDigits; }

 private:
  Elec2DetMapper mElec2Det{nullptr};
  FeeLink2SolarMapper mFee2Solar{nullptr};
  o2::mch::raw::PageDecoder mDecoder;
  size_t mNrdhs{0};
  std::vector<o2::mch::Digit> mOutputDigits;

  SampaChannelHandler mChannelHandler;
  std::function<void(o2::header::RDHAny*)> mRdhHandler;

  std::string mMapCRUfile;
  std::string mMapFECfile;

  bool mPrint{false};
  bool mDs2manu{false};

  MergerBase* mMerger = {nullptr};
};

} // namespace raw
} // namespace mch
} // end namespace o2
