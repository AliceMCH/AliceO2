// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file PreClusterFinderDevice.h
/// \brief Definition of a DPL device to run the preclusterizer
///
/// \author Philippe Pillot, Subatech

#ifndef O2_MCH_PRECLUSTERFINDERDEVICE_H_
#define O2_MCH_PRECLUSTERFINDERDEVICE_H_

#include "Framework/CallbackService.h"
#include "Framework/ControlService.h"
#include "Framework/Task.h"
#include "Framework/runDataProcessing.h"

using namespace o2;
using namespace o2::framework;

WorkflowSpec defineDataProcessing(const ConfigContext&);

#endif // O2_MCH_PRECLUSTERFINDERSPEC_H_
