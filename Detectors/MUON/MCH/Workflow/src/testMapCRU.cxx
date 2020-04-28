// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#define BOOST_TEST_MODULE Test MCHWorkflow MapCRU
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>
#include "MapCRU.h"
#include <sstream>
#include <string>

struct STREAM {
  // cruId cruLink s/n endpoint endpoint0_addr endpoint1_addr
  //   #    0..11  #    0..1     xx:#.#           yy:#.#
  std::stringstream refStream;
  STREAM()
  {
    refStream << R"(
1       0       0       0       af:0.0  bf:0.0
2       0      11       1       af:0.0  bf:0.0
)";
  }
};

BOOST_FIXTURE_TEST_SUITE(mapcru, STREAM)

BOOST_AUTO_TEST_CASE(DefautCtorCreatesAnEmptyMap)
{
  o2::mch::MapCRU cru;
  BOOST_CHECK_EQUAL(cru.size(), 0);
}

BOOST_AUTO_TEST_CASE(RequestingInvalidLinkMustYieldInvalidDeDs)
{
  o2::mch::MapCRU cru;
  BOOST_CHECK_EQUAL(cru.load(refStream), 2);
}

BOOST_AUTO_TEST_SUITE_END()
