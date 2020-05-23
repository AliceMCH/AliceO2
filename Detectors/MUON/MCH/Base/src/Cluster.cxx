// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include "MCHBase/Cluster.h"
#include <cmath>

namespace o2::mch
{

static bool closeEnough(double x, double y, double eps = 1E-6)
{
  return std::fabs(x - y) <= eps * std::max(1.0, std::max(std::fabs(x), std::fabs(y)));
}

Cluster::Cluster(float time, int detid, float x, float y, float ex, float ey, float charge)
  : mTime(time), mDetID(detid), mX(x), mY(y), mSigmaX(ex), mSigmaY(ey), mCharge(charge)
{
}

bool Cluster::operator==(const Cluster& other) const
{
  return mDetID == other.mDetID &&
         closeEnough(mX, other.mX) &&
         closeEnough(mY, other.mY) &&
         closeEnough(mSigmaX, other.mSigmaX) &&
         closeEnough(mSigmaY, other.mSigmaY) &&
         closeEnough(mCharge, other.mCharge) &&
         closeEnough(mTime, other.mTime);
}

} // namespace o2::mch
