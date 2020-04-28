// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include "MapFEC.h"
#include <iostream>
#include <fstream>
#include <algorithm>

namespace o2::mch
{

MapFEC::MapFEC()
{
  mDsMap.fill({-1, -1, -1});
}

int MapFEC::index(uint32_t linkId, uint32_t dsAddr) const
{
  if (linkId < 0 || linkId > sMaxLinkId) {
    return -1;
  }
  if (dsAddr < 0 || dsAddr >= sMaxDs) {
    return -1;
  }
  return linkId * sMaxDs + dsAddr;
}

int MapFEC::load(std::istream& in)
{
  int link_id, group_id, de, ds_id[5];
  while (in >> link_id >> group_id >> de >> ds_id[0] >> ds_id[1] >> ds_id[2] >> ds_id[3] >> ds_id[4]) {
    for (int i = 0; i < 5; i++) {
      if (ds_id[i] <= 0) {
        continue;
      }
      int ds_addr = group_id * 5 + i;
      int ix = index(link_id, ds_addr);
      if (ix < 0) {
        continue;
      }
      mDsMap.at(ix) = {de, ds_id[i], 0};
    }
  }
  return size();
}

size_t MapFEC::size() const
{
  return std::count_if(mDsMap.begin(), mDsMap.end(), [](const MapDualSampa& m) {
    return m.deId >= 0 && m.dsId >= 0 && m.bad == 0;
  });
}

bool MapFEC::getDsId(uint32_t link_id, uint32_t ds_addr, int& de, int& dsid) const
{
  if (!size()) {
    return false;
  }

  int ix = index(link_id, ds_addr);
  std::cout << "link_id=" << link_id << " ds_addr=" << ds_addr
            << " ix=" << ix << "\n";
  if (ix < 0) {
    return false;
  }
  if (mDsMap.at(ix).bad == 1) {
    return false;
  }
  de = mDsMap.at(ix).deId;
  dsid = mDsMap.at(ix).dsId;
  return de >= 0 && dsid >= 0;
}
} // namespace o2::mch
