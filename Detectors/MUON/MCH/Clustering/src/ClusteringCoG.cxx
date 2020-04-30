// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.


#include <chrono>
#include <memory>
#include <stdexcept>
#include <iostream>
#include <vector>
#include <TMath.h>
#include <stdio.h>
#include <math.h>

#include "MCHClustering/ClusteringCoG.h"
#include "MCHMappingInterface/Segmentation.h"

namespace o2
{
namespace mch
{

using namespace std;


//_________________________________________________________________________________________________
void ClusteringCoG::run(gsl::span<const PreCluster> preClusters, gsl::span<const Digit> digits, std::vector<Cluster>& clusters)
{
  // loop on pre-clusters
  for (const auto& preCluster : preClusters) {
    // we remove the preclusters with only one hit
    if(preCluster.nDigits < 2) {
      continue;
    }

    // get the digits of this precluster
    auto preClusterDigits = digits.subspan(preCluster.firstDigit, preCluster.nDigits);

    // CoG algorithm
    auto cluster = run(preClusterDigits);
    clusters.push_back(cluster);
  }
}

//_________________________________________________________________________________________________
Cluster ClusteringCoG::run(gsl::span<const Digit> precluster)
{
  double xmin = 1E9;
  double ymin = 1E9;
  double xmax = -1E9;
  double ymax = -1E9;
  double charge[] = { 0.0, 0.0 };
  int multiplicity[] = { 0, 0 };

  double x[] = { 0.0, 0.0 };
  double y[] = { 0.0, 0.0 };

  double xsize[] = { 0.0, 0.0 };
  double ysize[] = { 0.0, 0.0 };

  int detid = precluster[0].getDetID();
  const mapping::Segmentation& segment = mapping::segmentation(detid);

  for ( size_t i = 0; i < precluster.length(); ++i ) {
    const Digit& digit = precluster[i];
    int padid = digit.getPadID();

    // position and size of current pad
    double padPosition[2] = {segment.padPositionX(padid), segment.padPositionY(padid)};
    double padSize[2] = {segment.padSizeX(padid), segment.padSizeY(padid)};

    // update of xmin/max et ymin/max
    xmin = std::min(padPosition[0]-0.5*padSize[0],xmin);
    xmax = std::max(padPosition[0]+0.5*padSize[0],xmax);
    ymin = std::min(padPosition[1]-0.5*padSize[1],ymin);
    ymax = std::max(padPosition[1]+0.5*padSize[1],ymax);

    // cathode index
    int cathode = segment.isBendingPad(padid) ? 0 : 1;

    // update of the cluster position, size, charge and multiplicity
    x[cathode] += padPosition[0] * digit.getADC();
    y[cathode] += padPosition[1] * digit.getADC();
    xsize[cathode] += padSize[0];
    ysize[cathode] += padSize[1];
    charge[cathode] += digit.getADC();
    multiplicity[cathode] += 1;
  }

  // Computation of the CoG coordinates for the two cathodes
  for ( int cathode = 0; cathode < 2; ++cathode ) {
    if ( charge[cathode] != 0 ) {
      x[cathode] /= charge[cathode];
      y[cathode] /= charge[cathode];
    }
    if ( multiplicity[cathode] != 0 ) {
      double sqrtCharge = sqrt(charge[cathode]);
      xsize[cathode] /= (multiplicity[cathode] * sqrtCharge);
      ysize[cathode] /= (multiplicity[cathode] * sqrtCharge);
    } else {
      xsize[cathode] = 1E9;
      ysize[cathode] = 1E9;
    }
  }

  // each CoG coordinate is taken from the cathode with the best precision
  double xCOG = ( xsize[0] < xsize[1] ) ? x[0] : x[1];
  double yCOG = ( ysize[0] < ysize[1] ) ? y[0] : y[1];
  double ex = ( xsize[0] < xsize[1] ) ? xsize[0] : xsize[1];
  double ey = ( ysize[0] < ysize[1] ) ? ysize[0] : ysize[1];
  // For the moment we use the time of the pad with the highest charge
  double timestamp = precluster[0].getTimeStamp();

  return Cluster{timestamp, precluster[0].getDetID(), xCOG, yCOG, ex, ey, charge[0] + charge[1]};
}

} // namespace mch
} // namespace o2
