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
/// \file    runFileReader.cxx
/// \author  Andrea Ferrero
///
/// \brief This is an executable that reads a data file from disk and sends the data to QC via DPL.
///
/// This is an executable that reads a data file from disk and sends the data to QC via the Data Processing Layer.
/// It can be used as a data source for QC development. For example, one can do:
/// \code{.sh}
/// o2-qc-run-file-reader --infile=some_data_file | o2-qc --config json://${QUALITYCONTROL_ROOT}/etc/your_config.json
/// \endcode
///

#include <random>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include "Framework/CallbackService.h"
#include "Framework/ControlService.h"
#include "Framework/Task.h"
#include "Framework/runDataProcessing.h"
#include "DPLUtils/DPLRawParser.h"
#include "MCHBase/Digit.h"
#include "DataDecoderSpec.h"
#include "RawBufferDecoder.h"

// Dans ce code, on récupère un infut aui est un message avec le buffer, on fait tourner le code de base decodeBuffer qui est dans Handlers, et on renvoir un message de sortie (inspiré de FileReader de Andrea)

using namespace o2;
using namespace o2::framework;

class DataDecoderTask
{
  RawBufferDecoder<SampleMode, RDHv4> decoder;

public:
  //_________________________________________________________________________________________________
  void init(framework::InitContext& ic)
  {
    // Rien à initialiser
  }

  //_________________________________________________________________________________________________
  void run(framework::ProcessingContext& pc)
  {
    bool verbose = false;
    // get the input buffer

    std::vector<o2::mch::Digit> digits;

    auto& inputs = pc.inputs();

    DPLRawParser parser(inputs, o2::framework::select("TF:MCH/RAWDATA"));

    for (auto it = parser.begin(), end = parser.end(); it != end; ++it) {
      // retrieving RDH v4
      auto const* rdh = it.get_if<o2::header::RAWDataHeaderV4>();
      // retrieving the raw pointer of the page
      auto const* raw = it.raw();
      // retrieving payload pointer of the page
      auto const* payload = it.data();
      // size of payload
      size_t payloadSize = it.size();
      // offset of payload in the raw page
      size_t offset = it.offset();

      if( payloadSize == 0 ) continue;


      //std::cout<<"\n\npayloadSize: "<<payloadSize<<std::endl;
      //std::cout<<"raw:     "<<(void*)raw<<std::endl;
      //std::cout<<"payload: "<<(void*)payload<<std::endl;

     //  std::vector<uint8_t> buffer((uint8_t*)raw, ((uint8_t*)raw)+payloadSize+sizeof(o2::header::RAWDataHeaderV4));
      gsl::span<uint8_t> buffer(const_cast<uint8_t*>(reinterpret_cast<const uint8_t*>(raw)), payloadSize+sizeof(o2::header::RAWDataHeaderV4));
      decoder.decodeBuffer(buffer, digits);
    }

    for (auto&& input : pc.inputs()) {
      //QcInfoLogger::GetInstance() << "run RawDataProcessor: input " << input.spec->binding << AliceO2::InfoLogger::InfoLogger::endm;

      if (input.spec->binding != "readout")
        continue;

      const auto* header = o2::header::get<header::DataHeader*>(input.header);
      if (!header)
        continue;

      //std::cout<<"payloadSize: "<<header->payloadSize<<std::endl;

      std::vector<uint8_t> buffer((uint8_t*)input.payload, ((uint8_t*)input.payload)+header->payloadSize);
      decoder.decodeBuffer(buffer, digits);
    }

    if(false) {
    for (auto d : digits) {
      std::cout <<
          " DE# " << d.getDetID() <<
          " PadId " << d.getPadID() <<
          " ADC " << d.getADC() <<
          " time "<< d.getTimeStamp() <<
          std::endl;
    }
    }

    const size_t OUT_SIZE = sizeof(o2::mch::Digit) * digits.size();

    /// send the output buffer via DPL
    char* outbuffer = NULL;
    outbuffer = (char*)realloc(outbuffer, OUT_SIZE);
    memcpy(outbuffer, digits.data(), OUT_SIZE);

    // create the output message
    auto freefct = [](void* data, void*) { free(data); };
    pc.outputs().adoptChunk(Output{ "MCH", "DIGITS", 0 }, outbuffer, OUT_SIZE, freefct, nullptr);

    /*
    auto msgIn = pc.inputs().get<gsl::span<char>>("readout");
    auto bufferPtrIn = msgIn.data();
    auto sizeIn = msgIn.size();

    //if(verbose)
      std::cout << "We got the input buffer, sizeIn: " << sizeIn << std::endl;

    //std::vector<uint8_t> buffer(sizeIn);
    std::vector<uint8_t> buffer((uint8_t*)bufferPtrIn, ((uint8_t*)bufferPtrIn)+sizeIn);

    // Decode the buffer
    size_t outsize;
    char* outbuffer = decoder.decodeBuffer(buffer, outsize);



    /// send the output buffer via DPL

    const int OUT_SIZE = outsize;

    // create the output message
    auto freefct = [](void* data, void*) { free(data); };
    pc.outputs().adoptChunk(Output{ "MCH", "DIGITS", 0 }, (char*)outbuffer, OUT_SIZE, freefct, nullptr);
    //exit(0);
    */
  }

private:
  std::ifstream mInputFile{}; ///< input file
  bool mPrint = false;        ///< print digits

};

// clang-format off
WorkflowSpec defineDataProcessing(const ConfigContext&)
{
  WorkflowSpec specs;

  DataProcessorSpec producer{
    "DataDecoder",
    //Inputs{InputSpec{"readout", "ROUT", "RAWDATA", Lifetime::Timeframe}},
    o2::framework::select("TF:MCH/RAWDATA"),
    //o2::framework::select("readout:ROUT/RAWDATA"),
    Outputs{OutputSpec{"MCH", "DIGITS", 0, Lifetime::Timeframe}},
    AlgorithmSpec{adaptFromTask<DataDecoderTask>()},
    Options{}

  };
  specs.push_back(producer);

  return specs;
}
// clang-format on
