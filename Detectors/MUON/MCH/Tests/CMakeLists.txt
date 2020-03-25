# Copyright CERN and copyright holders of ALICE O2. This software is distributed
# under the terms of the GNU General Public License v3 (GPL Version 3), copied
# verbatim in the file "COPYING".
#
# See http://alice-o2.web.cern.ch/license for full licensing information.
#
# In applying this license CERN does not waive the privileges and immunities
# granted to it by virtue of its status as an Intergovernmental Organization or
# submit itself to any jurisdiction.

o2_add_executable(clusterizer-test
                  COMPONENT_NAME mch
                  SOURCES src/ClusterizerTest.cxx src/PreClusterSerializer.cxx src/DigitsFileReader.cxx
                  PUBLIC_LINK_LIBRARIES O2::MCHBase O2::MCHPreClustering O2::MCHClustering)
                  
o2_add_executable(clusterizer-test-tb
COMPONENT_NAME mch
SOURCES src/ClusterizerTestTB.cxx src/TBDigitsFileReader.cxx
PUBLIC_LINK_LIBRARIES O2::MCHBase O2::MCHPreClustering O2::MCHClustering)

o2_add_executable(clusterizer-validation
                  COMPONENT_NAME mch
                  SOURCES src/ClusterizerValidation.cxx src/Validation.cxx
                  PUBLIC_LINK_LIBRARIES O2::MCHBase O2::MCHPreClustering O2::MCHClustering)
