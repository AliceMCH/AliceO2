# Copyright CERN and copyright holders of ALICE O2. This software is distributed
# under the terms of the GNU General Public License v3 (GPL Version 3), copied
# verbatim in the file "COPYING".
#
# See http://alice-o2.web.cern.ch/license for full licensing information.
#
# In applying this license CERN does not waive the privileges and immunities
# granted to it by virtue of its status as an Intergovernmental Organization or
# submit itself to any jurisdiction.

o2_add_executable(print-preclusters-run2
                  COMPONENT_NAME mch
                  SOURCES src/PrintPreclustersRun2.cxx src/DigitsFileReader.cxx
                  PUBLIC_LINK_LIBRARIES O2::MCHBase O2::MCHMappingInterface O2::MCHMappingImpl3 
                  O2::MCHPreClustering O2::MCHClustering)


o2_add_executable(
        digits-reader-tb-workflow
        SOURCES src/TBDigitsFileReader.cxx src/TBDigitsFileReaderSpec.cxx src/digits-reader-tb-workflow.cxx
        COMPONENT_NAME mch
        PUBLIC_LINK_LIBRARIES O2::Framework O2::DPLUtils Boost::program_options O2::MCHBase O2::MCHMappingInterface O2::MCHMappingImpl3)
                  
o2_add_executable(clusterizer-test-tb
                  COMPONENT_NAME mch
                  SOURCES src/ClusterizerTestTB.cxx src/TBDigitsFileReader.cxx
                  PUBLIC_LINK_LIBRARIES O2::MCHBase O2::MCHPreClustering O2::MCHClustering)

o2_add_executable(
        digits2-sink-workflow
        SOURCES src/digits2-sink-workflow.cxx
        COMPONENT_NAME mch
        PUBLIC_LINK_LIBRARIES O2::Framework O2::DPLUtils Boost::program_options O2::MCHBase O2::MCHMappingInterface O2::MCHMappingImpl3)

o2_add_executable(
        preclusters-sink-alt-workflow
        SOURCES src/PreClusterSinkSpec.cxx src/preclusters-sink-workflow.cxx
        COMPONENT_NAME mch
        PUBLIC_LINK_LIBRARIES  O2::MCHBase O2::MCHMappingInterface O2::MCHMappingImpl3 O2::MCHPreClustering O2::MCHClustering)

#o2_add_executable(print-preclusters-tb
#                  COMPONENT_NAME mch
#                  SOURCES src/PrintPreclustersTB.cxx src/TBDigitsFileReader.cxx
#                  PUBLIC_LINK_LIBRARIES O2::MCHBase O2::MCHPreClustering)
