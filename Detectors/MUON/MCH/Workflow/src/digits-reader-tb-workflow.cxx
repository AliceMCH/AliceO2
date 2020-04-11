// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

///
/// \file    runFileReader.cxx
/// \author  Andrea Ferrero
///
/// \brief This is an executable that reads digits from disk and sends the data to QC via DPL.
///
/// This is an executable that reads digits from disk and sends the data to QC via the Data Processing Layer.
/// The file with digits is a specially crafted binary file derived from the ROOT trees generated by the MCH test-beam analysis code
/// It can be used as a data source for pre-clustering development. For example, one can do:
/// \code{.sh}
/// o2-mch-digits-reader-tb-workflow --infile=TB_digits_data_file | o2-mch-preclustering-workflow | o2-mch-digits-to-preclusters-workflow | o2-mch-preclusters-sink-workflow
/// \endcode
///

#include "Framework/CallbackService.h"
#include "Framework/ControlService.h"
#include "Framework/Task.h"
#include "Framework/runDataProcessing.h"
#include "TBDigitsFileReaderSpec.h"

using namespace o2;
using namespace o2::framework;

WorkflowSpec defineDataProcessing(const ConfigContext&)
{
  WorkflowSpec specs;

  DataProcessorSpec producer = o2::mch::getTBDigitsFileReaderSpec();
  specs.push_back(producer);

  return specs;
}
