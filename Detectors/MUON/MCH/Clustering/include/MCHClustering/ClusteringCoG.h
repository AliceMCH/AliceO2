// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// @since 2016-10-20
/// @author P. Pillot
/// @brief Class to fit a preclusters with the center-of-gravity method

#ifndef ALICEO2_MCH_CLUSTERINGCOG_H_
#define ALICEO2_MCH_CLUSTERINGCOG_H_

#include <vector>

#include <gsl/span>

#include "MCHBase/Digit.h"
#include "MCHBase/PreCluster.h"
#include "MCHBase/Cluster.h"

namespace o2
{
namespace mch
{

// Classes Cluster

class ClusteringCoG
{
 public:
  ClusteringCoG() = default;
  ~ClusteringCoG() = default;

  ClusteringCoG(const ClusteringCoG&) = delete;
  ClusteringCoG& operator=(const ClusteringCoG&) = delete;
  ClusteringCoG(ClusteringCoG&&) = delete;
  ClusteringCoG& operator=(ClusteringCoG&&) = delete;

  void run(gsl::span<const PreCluster> preClusters, gsl::span<const Digit> digits, std::vector<Cluster>& clusters);
  Cluster run(gsl::span<const Digit> precluster);
};

} // namespace mch
} // namespace o2

#endif /* ALICEO2_MCH_CLUSTERINGCOG_H_ */
