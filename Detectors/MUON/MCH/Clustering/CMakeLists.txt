# Copyright CERN and copyright holders of ALICE O2. This software is distributed
# under the terms of the GNU General Public License v3 (GPL Version 3), copied
# verbatim in the file "COPYING".
#
# See http://alice-o2.web.cern.ch/license for full licensing information.
#
# In applying this license CERN does not waive the privileges and immunities
# granted to it by virtue of its status as an Intergovernmental Organization or
# submit itself to any jurisdiction.

o2_add_library(MCHClustering
               SOURCES src/ClusteringForTest.cxx src/MyObjectDigits.cxx src/MyObjectParametersMathieson.cxx src/MyObjectParametersGauss.cxx src/MyObjectParametersDoubleGauss.cxx
               PUBLIC_LINK_LIBRARIES O2::MCHMappingInterface O2::SimulationDataFormat O2::MCHBase O2::Framework O2::MCHPreClustering O2::MCHMappingInterface O2::MCHMappingFactory O2::MCHMappingImpl3)
               
o2_target_root_dictionary(MCHClustering
               HEADERS include/MCHClustering/MyObjectDigits.h
                       include/MCHClustering/MyObjectParametersMathieson.h
                       include/MCHClustering/MyObjectParametersGauss.h
                       include/MCHClustering/MyObjectParametersDoubleGauss.h
               LINKDEF include/MCHClustering/LinkDef.h)

o2_add_executable(
        clustering
        SOURCES src/ClusteringSpec.cxx src/ClusteringDevice.cxx
        COMPONENT_NAME mch
        PUBLIC_LINK_LIBRARIES O2::Framework O2::DPLUtils O2::MCHClustering)
