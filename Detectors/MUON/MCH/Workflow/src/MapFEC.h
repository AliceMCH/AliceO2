// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#ifndef O2_MCH_WORKFLOW_MAPFEC_H
#define O2_MCH_WORKFLOW_MAPFEC_H

#include <string>
#include <array>
#include <istream>

namespace o2::mch
{
class MapFEC
{
 public:
  MapFEC();
  int load(std::istream& in);
  bool getDsId(uint32_t link_id, uint32_t ds_addr, int& de, int& dsid) const;
  size_t size() const;

 private:
  struct MapDualSampa {
    int deId = -1; // detector element
    int dsId = -1; // DS index
    int bad = -1;  // if = 1 bad pad (not used for analysis)
  };
  int index(uint32_t linkId, uint32_t dsAddr) const;
  static constexpr int sMaxLinkId = 0x7ff;
  static constexpr int sMaxDs = 40;
  std::array<MapDualSampa, (sMaxLinkId + 1) * sMaxDs> mDsMap;
};

} // namespace o2::mch

#endif
