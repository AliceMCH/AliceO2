// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include "MapCRU.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <limits>

namespace o2::mch
{

MapCRU::MapCRU()
{
  mFeeLink2Solar.fill(std::numeric_limits<uint16_t>::max());
}

int MapCRU::indexFeeLink(int feeid, int linkid) const
{
  if (feeid < 0 || feeid >= sMaxFeeId) {
    return -1;
  }
  if (linkid < 0 || linkid >= sMaxLinkId) {
    return -1;
  }
  return feeid * sMaxLinkId + linkid;
}

int MapCRU::load(std::istream& in)
{
  int f, l, link_id, dummy;
  std::string a1, a2;
  std::string s;
  while (std::getline(in, s)) {
    if (s.empty()) {
      continue;
    }
    std::istringstream line(s);
    line >> link_id >> f >> l >> dummy >> a1 >> a2;
    auto ix = indexFeeLink(f, l);
    if (ix < 0) {
      continue;
    }
    mFeeLink2Solar.at(ix) = link_id;
  }
  return size();
}

size_t MapCRU::size() const
{
  return std::count_if(mFeeLink2Solar.begin(), mFeeLink2Solar.end(), [](uint16_t a) { return a != std::numeric_limits<uint16_t>::max(); });
}

std::optional<uint16_t> MapCRU::operator()(const o2::mch::raw::FeeLinkId& feeLinkId) const
{
  if (!size()) {
    return std::nullopt;
  }
  auto ix = indexFeeLink(feeLinkId.feeId(), feeLinkId.linkId());
  if (ix < 0) {
    return std::nullopt;
  }
  return mFeeLink2Solar.at(ix);
}
} // namespace o2::mch
