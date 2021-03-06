# Copyright CERN and copyright holders of ALICE O2. This software is distributed
# under the terms of the GNU General Public License v3 (GPL Version 3), copied
# verbatim in the file "COPYING".
#
# See http://alice-o2.web.cern.ch/license for full licensing information.
#
# In applying this license CERN does not waive the privileges and immunities
# granted to it by virtue of its status as an Intergovernmental Organization or
# submit itself to any jurisdiction.

o2_add_library(MCHRawEncoderPayload
        SOURCES
        BitSet.cxx
        DataBlock.cxx
        PayloadEncoder.cxx
        CruBufferCreator.cxx
        PayloadPaginator.cxx
        EncoderImplHelper.cxx
        PUBLIC_LINK_LIBRARIES
        O2::CommonDataFormat
        O2::DetectorsRaw
        O2::MCHBase
        O2::MCHMappingImpl4
        O2::MCHRawCommon
        O2::MCHRawElecMap
        PRIVATE_LINK_LIBRARIES O2::MCHRawImplHelpers)

o2_add_executable(raw-create-test-buffers
        SOURCES createTestBuffers.cxx
        COMPONENT_NAME mch
        PUBLIC_LINK_LIBRARIES O2::MCHRawImplHelpers O2::MCHRawEncoderPayload
        TARGETVARNAME createTestBuffers)

if(BUILD_TESTING)

        o2_add_test(elink-encoder
                SOURCES testElinkEncoder.cxx
                COMPONENT_NAME mchraw
                LABELS "muon;mch;raw"
                PUBLIC_LINK_LIBRARIES O2::MCHRawEncoderPayload O2::MCHRawImplHelpers)

        o2_add_test(gbt-encoder
                SOURCES testGBTEncoder.cxx
                COMPONENT_NAME mchraw
                LABELS "muon;mch;raw"
                PUBLIC_LINK_LIBRARIES O2::MCHRawEncoderPayload O2::MCHTestRawRefBuffers
                O2::MCHRawImplHelpers)

        o2_add_test(payload-encoder
                SOURCES testPayloadEncoder.cxx 
                COMPONENT_NAME mchraw
                LABELS "muon;mch;raw"
                PUBLIC_LINK_LIBRARIES O2::MCHRawEncoderPayload O2::MCHRawImplHelpers)

        o2_add_test(bitset
                SOURCES testBitSet.cxx
                COMPONENT_NAME mchraw
                LABELS "muon;mch;raw"
                PUBLIC_LINK_LIBRARIES O2::MCHRawEncoderPayload)

        o2_add_test(bare-elink-encoder
                SOURCES testBareElinkEncoder.cxx
                COMPONENT_NAME mchraw
                LABELS "muon;mch;raw"
                PUBLIC_LINK_LIBRARIES O2::MCHRawEncoderPayload O2::MCHRawImplHelpers)

        add_test(NAME o2-test-mch-create-test-buffers-userlogic COMMAND ${createTestBuffers} -u)
        set_property(TEST o2-test-mch-create-test-buffers-userlogic APPEND PROPERTY LABELS "muon;mch;raw")
        add_test(NAME o2-test-mch-create-test-buffers-bare COMMAND ${createTestBuffers})
        set_property(TEST o2-test-mch-create-test-buffers-bare APPEND PROPERTY LABELS "muon;mch;raw")
                
endif()

if(benchmark_FOUND)
        o2_add_executable(bitset
                COMPONENT_NAME mchraw
                SOURCES benchBitSet.cxx
                PUBLIC_LINK_LIBRARIES O2::MCHRawEncoderPayload benchmark::benchmark
                IS_BENCHMARK)
endif()

