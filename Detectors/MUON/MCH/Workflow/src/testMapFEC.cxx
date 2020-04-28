// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#define BOOST_TEST_MODULE Test MCHWorkflow MapFEC
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>
#include "MapFEC.h"
#include <sstream>
#include <string>

struct STREAM {
  // link group   de      dsid0   dsid1   dsid2   dsid3   dsid4
  std::stringstream refStream;
  STREAM()
  {
    refStream << R"(
1       0       819     108     0       107     0       106
1       1       919     1133    0       1134    0       0
)";
  }
};

BOOST_FIXTURE_TEST_SUITE(mapfec, STREAM)

BOOST_AUTO_TEST_CASE(RequestingInvalidLinkMustYieldInvalidDeDs)
{
  o2::mch::MapFEC fec;
  BOOST_CHECK_EQUAL(fec.load(refStream), 5);
  int deId, dsId;
  bool ok = fec.getDsId(2, 0, deId, dsId);
  BOOST_CHECK_EQUAL(ok, false);
  BOOST_CHECK_EQUAL(deId, -1);
  BOOST_CHECK_EQUAL(dsId, -1);
}

BOOST_AUTO_TEST_CASE(RequestingInvalidDsAddressMustYieldInvalidDeDs)
{
  o2::mch::MapFEC fec;
  BOOST_CHECK_EQUAL(fec.load(refStream), 5);
  int deId, dsId;
  bool ok = fec.getDsId(1, 6, deId, dsId);
  BOOST_CHECK_EQUAL(ok, false);
  BOOST_CHECK_EQUAL(deId, -1);
  BOOST_CHECK_EQUAL(dsId, -1);
}

BOOST_AUTO_TEST_CASE(GetValidDeDs)
{
  o2::mch::MapFEC fec;
  BOOST_CHECK_EQUAL(fec.load(refStream), 5);
  int deId, dsId;
  bool ok = fec.getDsId(1, 7, deId, dsId);
  BOOST_CHECK_EQUAL(ok, true);
  BOOST_CHECK_EQUAL(deId, 919);
  BOOST_CHECK_EQUAL(dsId, 1134);
}

BOOST_AUTO_TEST_SUITE_END()
