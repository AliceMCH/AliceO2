# Copyright CERN and copyright holders of ALICE O2. This software is distributed
# under the terms of the GNU General Public License v3 (GPL Version 3), copied
# verbatim in the file "COPYING".
#
# See http://alice-o2.web.cern.ch/license for full licensing information.
#
# In applying this license CERN does not waive the privileges and immunities
# granted to it by virtue of its status as an Intergovernmental Organization or
# submit itself to any jurisdiction.

o2_add_library(MCHCalibration
         SOURCES
	       src/PedestalDigit.cxx
	       PUBLIC_LINK_LIBRARIES ROOT::Core FairRoot::Base FairMQ::FairMQ ms_gsl::ms_gsl O2::MCHBase O2::MCHRawDecoder)

o2_target_root_dictionary(MCHCalibration
  HEADERS include/MCHCalibration/PedestalDigit.h)

o2_add_executable(
        pedestals-workflow
        SOURCES src/pedestals-workflow.cxx
        COMPONENT_NAME mch
        PUBLIC_LINK_LIBRARIES Boost::program_options O2::Framework O2::DPLUtils O2::MCHCalibration)
  