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
/// @author  Laurent Aphecetche, Sebastien Perrin

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
#include "Framework/DataProcessorSpec.h"
#include "RawBufferDecoder.h"
#include "DumpOptionsAndStat.h"

namespace po = boost::program_options;


using namespace o2::framework;
using namespace o2::mch::mapping;
using namespace o2::mch::raw;
using RDHv4 = o2::header::RAWDataHeaderV4;

    //  Récupère un fichier sur l'ordinateur et utilise le bout de code de base decodeBuffer pour décoder
    
    
template <typename CHARGESUM, typename RDH>
void digitdump(std::string input, DumpOptions opt)
{
  std::ifstream in(input.c_str(), std::ios::binary);
  if (!in.good()) {
    std::cout << "could not open file " << input << "\n";
    return;
  }
  constexpr size_t pageSize = 8192;

    std::vector<uint8_t> buffer(pageSize);

  size_t npages{0};
  size_t outsize;

  RawBufferDecoder<CHARGESUM, RDH> decoder(false);
  while (npages < opt.maxNofRDHs() && in.read(reinterpret_cast<char*>(&buffer[0]), pageSize)) {
    npages++;
    //char* outbuffer = decoder.decodeBuffer(buffer, outsize);
  }
    
  return;
}
    
    
int main(int argc, char* argv[])
{
  std::string prefix;
  std::vector<int> detElemIds;
  std::string inputFile;
  po::variables_map vm;
  po::options_description generic("Generic options");
  unsigned int nrdhs{0};
  unsigned int deId{0};
  bool showRDHs{false};
  bool userLogic{false}; // default format is bareformat...
  bool chargeSum{false}; //... in sample mode
  bool jsonOutput{false};

  // clang-format off
  generic.add_options()
      ("help,h", "produce help message")
      ("input-file,i", po::value<std::string>(&inputFile)->required(), "input file name")
      ("nrdhs,n", po::value<unsigned int>(&nrdhs), "number of RDHs to go through")
      ("showRDHs,s",po::bool_switch(&showRDHs),"show RDHs")
      ("userLogic,u",po::bool_switch(&userLogic),"user logic format")
      ("chargeSum,c",po::bool_switch(&chargeSum),"charge sum format")
      ("json,j",po::bool_switch(&jsonOutput),"output means and rms in json format")
      ("de,d",po::value<unsigned int>(&deId),"detection element id of the data to be decoded")
      ("cru",po::value<uint16_t>(),"force cruId")
      ;
  // clang-format on

  po::options_description cmdline;
  cmdline.add(generic);

  po::store(po::command_line_parser(argc, argv).options(cmdline).run(), vm);

  if (vm.count("help")) {
    std::cout << generic << "\n";
    return 2;
  }

  try {
    po::notify(vm);
  } catch (boost::program_options::error& e) {
    std::cout << "Error: " << e.what() << "\n";
    exit(1);
  }

  DumpOptions opt(deId, nrdhs, showRDHs, jsonOutput);

  if (vm.count("cru")) {
    opt.cruId(vm["cru"].as<uint16_t>());
  }
  if (chargeSum) {
    digitdump<ChargeSumMode, RDHv4>(inputFile, opt);
  } else {
    digitdump<SampleMode, RDHv4>(inputFile, opt);
  }

  return 0;
}

