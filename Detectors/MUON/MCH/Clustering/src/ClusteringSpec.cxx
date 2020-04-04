// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file TBDigitsFileReaderSpec.cxx
/// \brief Implementation of a data processor to run the preclusterizer
///
/// \author Philippe Pillot, Subatech

#include "ClusteringSpec.h"

#include <string>
#include <chrono>

#include <stdexcept>

#include "Framework/CallbackService.h"
#include "Framework/ConfigParamRegistry.h"
#include "Framework/ControlService.h"
#include "Framework/DataProcessorSpec.h"
#include "Framework/Lifetime.h"
#include "Framework/Output.h"
#include "Framework/Task.h"

#include "MCHClustering/ClusteringForTest.h"

namespace o2
{
namespace mch
{

using namespace std;
using namespace o2::framework;

class ClusteringTask
{
public:
  //_________________________________________________________________________________________________
  void init(framework::InitContext& ic)
  {
    /// Prepare the preclusterizer
    LOG(INFO) << "[MCH] initializing clustering";

    mMethod = ic.options().get<int>("method");
    mPrint = ic.options().get<bool>("print");
  }
  //_________________________________________________________________________________________________
  void run(framework::ProcessingContext& pc)
  {
    // get the input buffers
    auto preClusters = pc.inputs().get<gsl::span<PreCluster>>("preclusters");
    auto digits = pc.inputs().get<gsl::span<Digit>>("preclusterdigits");

    if( mPrint ) {
      std::cout<<"Number of pre-clusters: "<<preClusters.length()<<std::endl;
    }

    Clustering clustering;
    std::vector<Clustering::Cluster> clusters(0);
    switch(mMethod) {
    case 0:
      clustering.runFinderCOG(preClusters, digits, clusters);
      break;
    case 1:
      clustering.runFinderGaussianFit(preClusters, digits, clusters);
      break;
    case 2:
      clustering.runFinderSimpleFit(preClusters, digits, clusters);
      break;
    }

    if( mPrint ) {
      std::cout<<"Number of clusters: "<<clusters.size()<<std::endl;
      int ci{0};
      for (const auto& cluster : clusters) {
        std::cout<<"Cluster "<<ci<<std::endl;
        std::cout<<"  x: "<<cluster.x<<"\n  y: "<<cluster.y<<"\n";
      }
    }

    size_t bufSize = sizeof(Clustering::Cluster) * clusters.size();
    char* clustersBuffer = (char*)malloc(bufSize);

    // create the output message
    auto freefct = [](void* data, void* /*hint*/) { free(data); };
    pc.outputs().adoptChunk(Output{ "MCH", "CLUSTERS" }, reinterpret_cast<char*>(clustersBuffer), bufSize, freefct, nullptr);
  }

private:
  int mMethod = 2;
  bool mPrint = false;                    ///< print preclusters
};

//_________________________________________________________________________________________________
DataProcessorSpec getClusteringSpec()
{
  return DataProcessorSpec{
    "Clustering",
    Inputs{InputSpec{"preclusters", "MCH", "PRECLUSTERS", 0, Lifetime::Timeframe},
      InputSpec{"preclusterdigits", "MCH", "PRECLUSTERDIGITS", 0, Lifetime::Timeframe}
    },
    Outputs{OutputSpec{"MCH", "CLUSTERS", 0, Lifetime::Timeframe}},
    AlgorithmSpec{adaptFromTask<ClusteringTask>()},
    Options{{"print", VariantType::Bool, false, {"print clusters"}},
      {"method", VariantType::Int, 2, {"clustering method [0=CoG, 1=Gauss, 2=Mathieson simple"}}}};
}

} // end namespace mch
} // end namespace o2
