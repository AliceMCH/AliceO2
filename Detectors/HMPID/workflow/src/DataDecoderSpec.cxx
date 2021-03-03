// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file   DataDecoderSpec.cxx
/// \author Antonio Franco - INFN Bari
/// \version 1.0
/// \date 01 feb 2021
/// \brief Implementation of a data processor to run the HMPID raw decoding
///

#include <random>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <array>
#include <functional>
#include <vector>

#include "TTree.h"
#include "TFile.h"

#include <gsl/span>

#include "Framework/CallbackService.h"
#include "Framework/ConfigParamRegistry.h"
#include "Framework/ControlService.h"
#include "Framework/DataProcessorSpec.h"
#include "Framework/Lifetime.h"
#include "Framework/Output.h"
#include "Framework/Task.h"
#include "Framework/WorkflowSpec.h"
#include "Framework/Logger.h"

#include "Headers/RAWDataHeader.h"
#include "DetectorsRaw/RDHUtils.h"
#include "DPLUtils/DPLRawParser.h"

#include "HMPIDBase/Digit.h"
#include "HMPIDBase/Geo.h"
#include "HMPIDReconstruction/HmpidDecodeRawMem.h"
#include "HMPIDWorkflow/DataDecoderSpec.h"

namespace o2
{
namespace hmpid
{

using namespace o2;
using namespace o2::framework;
using RDH = o2::header::RDHAny;

//=======================
// Data decoder
void DataDecoderTask::init(framework::InitContext& ic)
{

  LOG(INFO) << "[HMPID Data Decoder - Init] ( create Decoder for " << Geo::MAXEQUIPMENTS << " equipments !";

  mRootStatFile = ic.options().get<std::string>("root-file");
  mDeco = new o2::hmpid::HmpidDecodeRawDigit(Geo::MAXEQUIPMENTS);
  mDeco->init();
  mTotalDigits = 0;
  mTotalFrames = 0;

  mExTimer.start();
  return;
}

void DataDecoderTask::run(framework::ProcessingContext& pc)
{
  mDeco->mDigits.clear();
  decodeTF(pc);
  //  TODO: accept other types of Raw Streams ...
  //  decodeReadout(pc);
  // decodeRawFile(pc);

  mExTimer.elapseMes("... Decoding... Digits decoded = " + std::to_string(mTotalDigits) + " Frames received = " + std::to_string(mTotalFrames));
  return;
}

void DataDecoderTask::endOfStream(framework::EndOfStreamContext& ec)
{
  // Records the statistics
  float avgEventSize;    //[o2::hmpid::Geo::MAXEQUIPMENTS];
  float avgBusyTime;     //[o2::hmpid::Geo::MAXEQUIPMENTS];
  float numOfSamples;    //[o2::hmpid::Geo::N_MODULES][o2::hmpid::Geo::N_YCOLS][o2::hmpid::Geo::N_XROWS];
  float sumOfCharges;    //[o2::hmpid::Geo::N_MODULES][o2::hmpid::Geo::N_YCOLS][o2::hmpid::Geo::N_XROWS];
  float squareOfCharges; //[o2::hmpid::Geo::N_MODULES][o2::hmpid::Geo::N_YCOLS][o2::hmpid::Geo::N_XROWS];
  float xb;
  float yb;

  TString filename = TString::Format("%s_stat.root", mRootStatFile.c_str());
  LOG(INFO) << "Create the stat file " << filename.Data();
  TFile mfileOut(TString::Format("%s", filename.Data()), "RECREATE");
  TTree* theObj[Geo::N_MODULES + 1];
  for (int i = 0; i < Geo::N_MODULES; i++) { // Create the TTree array
    TString tit = TString::Format("HMPID Data Decoding Statistic results Mod=%d", i);
    theObj[i] = new TTree("o2hmp", tit);
    theObj[i]->Branch("x", &xb, "s");
    theObj[i]->Branch("y", &yb, "s");
    theObj[i]->Branch("Samples", &numOfSamples, "i");
    theObj[i]->Branch("Sum_of_charges", &sumOfCharges, "l");
    theObj[i]->Branch("Sum_of_square", &squareOfCharges, "l");
  }
  theObj[Geo::N_MODULES] = new TTree("o2hmp", "HMPID Data Decoding Statistic results");
  theObj[Geo::N_MODULES]->Branch("Average_Event_Size", &avgEventSize, "F");
  theObj[Geo::N_MODULES]->Branch("Average_Busy_Time", &avgBusyTime, "F");

  char summaryFileName[254];
  sprintf(summaryFileName, "%s_stat.txt", mRootStatFile.c_str());
  mDeco->writeSummaryFile(summaryFileName);
  int numEqui = mDeco->getNumberOfEquipments();
  for (int e = 0; e < numEqui; e++) {
    avgEventSize = mDeco->getAverageEventSize(e);
    avgBusyTime = mDeco->getAverageBusyTime(e);
    theObj[Geo::N_MODULES]->Fill();
  }
  for (int m = 0; m < o2::hmpid::Geo::N_MODULES; m++) {
    for (int y = 0; y < o2::hmpid::Geo::N_YCOLS; y++) {
      for (int x = 0; x < o2::hmpid::Geo::N_XROWS; x++) {
        xb = x;
        yb = y;
        numOfSamples = mDeco->getPadSamples(m, x, y);
        sumOfCharges = mDeco->getPadSum(m, x, y);
        squareOfCharges = mDeco->getPadSquares(m, x, y);
        theObj[m]->Fill();
      }
    }
  }
  for (int i = 0; i <= Geo::N_MODULES; i++) {
    theObj[i]->Write();
  }

  mExTimer.logMes("End the Decoding ! Digits decoded = " + std::to_string(mTotalDigits) + " Frames received = " + std::to_string(mTotalFrames));
  mExTimer.stop();
  //ec.services().get<ControlService>().endOfStream();
  // ec.services().get<o2::framework::ControlService>().readyToQuit(framework::QuitRequest::Me);
  return;
}
//_________________________________________________________________________________________________
// the decodeTF() function processes the the messages generated by the (sub)TimeFrame builder
void DataDecoderTask::decodeTF(framework::ProcessingContext& pc)
{
  LOG(DEBUG) << "*********** In decodeTF **************";

  // get the input buffer
  auto& inputs = pc.inputs();
  DPLRawParser parser(inputs, o2::framework::select("TF:HMP/RAWDATA"));
  for (auto it = parser.begin(), end = parser.end(); it != end; ++it) {
    mDeco->mDigits.clear();
    uint32_t* theBuffer = (uint32_t*)it.raw();
    mDeco->setUpStream(theBuffer, it.size() + it.offset());
    mDeco->decodePageFast(&theBuffer);
    mTotalFrames++;
    pc.outputs().snapshot(o2::framework::Output{"HMP", "DIGITS", 0, o2::framework::Lifetime::Timeframe}, mDeco->mDigits); //
    mTotalDigits += mDeco->mDigits.size();
    LOG(DEBUG) << "Writing " << mDeco->mDigits.size() << "/" << mTotalDigits << " Digits ...";
  }
  return;
}

//_________________________________________________________________________________________________
// the decodeReadout() function processes the messages generated by o2-mch-cru-page-reader-workflow
// TODO: rearrange, test
void DataDecoderTask::decodeReadout(framework::ProcessingContext& pc)
{
  LOG(INFO) << "*********** In decode readout **************";

  // get the input buffer
  auto& inputs = pc.inputs();
  DPLRawParser parser(inputs, o2::framework::select("readout:HMP/RAWDATA"));
  //  DPLRawParser parser(inputs, o2::framework::select("HMP/readout"));

  for (auto it = parser.begin(), end = parser.end(); it != end; ++it) {
    uint32_t* theBuffer = (uint32_t*)it.raw();
    mDeco->setUpStream(theBuffer, it.size() + it.offset());
    mDeco->decodePageFast(&theBuffer);
  }
  return;
}

// the decodeReadout() function processes the messages generated by o2-mch-cru-page-reader-workflow
// TODO: rearrange, test
void DataDecoderTask::decodeRawFile(framework::ProcessingContext& pc)
{
  LOG(INFO) << "*********** In decode rawfile **************";

  for (auto&& input : pc.inputs()) {
    if (input.spec->binding == "file") {
      const header::DataHeader* header = o2::header::get<header::DataHeader*>(input.header);
      if (!header) {
        return;
      }

      auto const* raw = input.payload;
      size_t payloadSize = header->payloadSize;

      LOG(INFO) << "  payloadSize=" << payloadSize;
      if (payloadSize == 0) {
        return;
      }

      uint32_t* theBuffer = (uint32_t*)input.payload;
      int pagesize = header->payloadSize;
      mDeco->setUpStream(theBuffer, pagesize);
      mDeco->decodePageFast(&theBuffer);
    }
  }
  return;
}

//_________________________________________________________________________________________________
o2::framework::DataProcessorSpec getDecodingSpec(std::string inputSpec)
{
  std::vector<o2::framework::InputSpec> inputs;
  inputs.emplace_back("TF", o2::framework::ConcreteDataTypeMatcher{"HMP", "RAWDATA"}, o2::framework::Lifetime::Timeframe);
  inputs.emplace_back("file", o2::framework::ConcreteDataTypeMatcher{"ROUT", "RAWDATA"}, o2::framework::Lifetime::Timeframe);
  //  inputs.emplace_back("readout", o2::header::gDataOriginHMP, "RAWDATA", 0, Lifetime::Timeframe);
  //  inputs.emplace_back("readout", o2::header::gDataOriginHMP, "RAWDATA", 0, Lifetime::Timeframe);
  //  inputs.emplace_back("rawfile", o2::header::gDataOriginHMP, "RAWDATA", 0, Lifetime::Timeframe);

  std::vector<o2::framework::OutputSpec> outputs;
  outputs.emplace_back("HMP", "DIGITS", 0, o2::framework::Lifetime::Timeframe);

  return DataProcessorSpec{
    "HMP-DataDecoder",
    o2::framework::select(inputSpec.c_str()),
    outputs,
    AlgorithmSpec{adaptFromTask<DataDecoderTask>()},
    Options{{"root-file", VariantType::String, "/tmp/hmpRawDecodeResults", {"Name of the Root file with the decoding results."}}}};
}

} // namespace hmpid
} // end namespace o2
