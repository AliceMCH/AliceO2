// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file   clusters-sink-workflow.cxx
/// \author Philippe Pillot, Subatech
/// \author  Andrea Ferrero, CEA
///
/// \brief This is an executable that dumps to a file on disk the clusters received via DPL.
///
/// This is an executable that dumps to a file on disk the clusters received via the Data Processing Layer.
/// It can be used to debug the clustering step. For example, one can do:
/// \code{.sh}
/// o2-mch-file-to-digits-workflow --infile=some_data_file | o2-mch-digits-to-preclusters-workflow | o2-mch-preclusters-to-clusters-workflow --method 0 | o2-mch-clusters-sink-workflow --outfile clusters.out
/// \endcode
///

#include "Framework/CallbackService.h"
#include "Framework/ControlService.h"
#include "Framework/Task.h"
#include "Framework/runDataProcessing.h"
#include "ClusterSinkSpec.h"

using namespace o2;
using namespace o2::framework;

WorkflowSpec defineDataProcessing(const ConfigContext&)
{
  WorkflowSpec specs;

  DataProcessorSpec producer = o2::mch::getClusterSinkSpec();
  specs.push_back(producer);

  return specs;
}
