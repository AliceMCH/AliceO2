// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/** @file Digit.h
 * C++ simple Muon MCH cluster.
 * @author  Andrea Ferrero, CEA
 * @author  Sebastien Perrin, CEA
 */

#ifndef ALICEO2_MCH_BASE_CLUSTER_H_
#define ALICEO2_MCH_BASE_CLUSTER_H_

#include "Rtypes.h"

namespace o2
{
namespace mch
{

// \class Cluster
/// \brief MCH cluster implementation
class Cluster
{
 public:
  Cluster() = default;

  Cluster(float time, int detid, float x, float y, float ex, float ey, float charge);
  ~Cluster() = default;

  bool operator==(const Cluster&) const;

  float getTimeStamp() const { return mTime; }

  int getDetID() const { return mDetID; }

  double getX() const { return mX; }
  double getY() const { return mY; }
  double getSigmaX() const { return mSigmaX; }
  double getSigmaY() const { return mSigmaY; }

  float getCharge() const { return mCharge; }

 private:
  int mDetID;     /// Detection element the cluster belongs to
  float mTime;    /// Time stamp of the charge cluster
  double mX;      /// cluster X coordinate
  double mY;      /// cluster Y coordinate
  double mSigmaX; /// error on X coordinate
  double mSigmaY; /// error on Y coordinate
  float mCharge;  /// Total charge in the cluster

  ClassDefNV(Cluster, 1);
}; //class Cluster

} //namespace mch
} //namespace o2
#endif // ALICEO2_MCH_BASE_CLUSTER_H_
