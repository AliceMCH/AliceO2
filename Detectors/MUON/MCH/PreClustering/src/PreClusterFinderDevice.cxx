// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file PreClusterFinderDevice.cxx
/// \brief Implementation of a DPL device to run the preclusterizer
///
/// \author Philippe Pillot, Subatech

#include "PreClusterFinderSpec.h"
#include "PreClusterFinderDevice.h"

// clang-format off
WorkflowSpec defineDataProcessing(const ConfigContext&)
{
  WorkflowSpec specs;

  DataProcessorSpec producer = o2::mch::getPreClusterFinderSpec();
  specs.push_back(producer);

  return specs;
}
// clang-format on
