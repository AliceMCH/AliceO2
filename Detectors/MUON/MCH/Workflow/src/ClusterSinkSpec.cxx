// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file ClusterSinkSpec.cxx
/// \brief Implementation of a data processor to write clusters
///
/// \author Philippe Pillot, Subatech
/// \author Andrea Ferrero, CEA

#include "ClusterSinkSpec.h"

#include <iostream>
#include <fstream>

#include <array>
#include <stdexcept>
#include <vector>

#include "Framework/CallbackService.h"
#include "Framework/ConfigParamRegistry.h"
#include "Framework/ControlService.h"
#include "Framework/DataProcessorSpec.h"
#include "Framework/Lifetime.h"
#include "Framework/Task.h"

#include "MCHBase/Digit.h"
#include "MCHBase/PreCluster.h"
#include "MCHBase/Cluster.h"
#include "MCHMappingInterface/Segmentation.h"

namespace o2
{
namespace mch
{

using namespace std;
using namespace o2::framework;

class ClusterSinkTask
{
 public:
  //_________________________________________________________________________________________________
  void init(framework::InitContext& ic)
  {
    /// Get the output file from the context
    LOG(INFO) << "initializing cluster sink";

    mText = ic.options().get<bool>("txt");

    auto outputFileName = ic.options().get<std::string>("outfile");
    mOutputFile.open(outputFileName, (mText ? ios::out : (ios::out | ios::binary)));
    if (!mOutputFile.is_open()) {
      throw invalid_argument("Cannot open output file" + outputFileName);
    }

    auto stop = [this]() {
      /// close the output file
      LOG(INFO) << "stop cluster sink";
      this->mOutputFile.close();
    };
    ic.services().get<CallbackService>().set(CallbackService::Id::Stop, stop);
  }

  //_________________________________________________________________________________________________
  void run(framework::ProcessingContext& pc)
  {
    /// dump the clusters with associated digits of the current event

    // get the input clusters and associated digits
    auto clusters = pc.inputs().get<gsl::span<Cluster>>("clusters");

    if (mText) {
      mOutputFile << clusters.size() << " clusters:" << endl;
      for (const auto& cluster : clusters) {
        mOutputFile << cluster.getDetID() << "  " << cluster.getCharge() << "  " << cluster.getTimeStamp()
                    << "  " << cluster.getX() << "," << cluster.getY() << " +/- " << cluster.getSigmaX() << "," << cluster.getSigmaY() << endl;
      }
    } else {
      // write the number of clusters
      int nClusters = clusters.size();
      mOutputFile.write(reinterpret_cast<char*>(&nClusters), sizeof(int));

      // write the clusters
      mOutputFile.write(reinterpret_cast<const char*>(clusters.data()), clusters.size_bytes());
    }
  }

 private:
  std::ofstream mOutputFile{}; ///< output file
  bool mText = false;          ///< output clusters in text format
};

//_________________________________________________________________________________________________
o2::framework::DataProcessorSpec getClusterSinkSpec()
{
  return DataProcessorSpec{
    "ClusterSink",
    Inputs{InputSpec{"clusters", "MCH", "CLUSTERS", 0, Lifetime::Timeframe}},
    Outputs{},
    AlgorithmSpec{adaptFromTask<ClusterSinkTask>()},
    Options{{"outfile", VariantType::String, "clusters.out", {"output filename"}},
            {"txt", VariantType::Bool, false, {"output clusters in text format"}}}};
}

} // end namespace mch
} // end namespace o2
