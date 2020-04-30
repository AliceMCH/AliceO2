// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file PreClusterSinkSpec.cxx
/// \brief Implementation of a data processor to write preclusters
///
/// \author Philippe Pillot, Subatech

#include "PreClusterSinkSpec.h"

#include <iostream>
#include <fstream>

#include <array>
#include <stdexcept>
#include <vector>

#include "Framework/CallbackService.h"
#include "Framework/ConfigParamRegistry.h"
#include "Framework/ControlService.h"
#include "Framework/DataProcessorSpec.h"
#include "Framework/Lifetime.h"
#include "Framework/Task.h"

#include "MCHBase/Digit.h"
#include "MCHBase/PreCluster.h"
#include "MCHMappingFactory/CreateSegmentation.h"
#include "MCHClustering/ClusteringForTest.h"

namespace o2
{
namespace mch
{

using namespace std;
using namespace o2::framework;


// \class Digit2
/// \brief MCH digit implementation, with pad coordinates
struct Digit2
{
 public:
  uint64_t mTime;
  uint32_t mDetID;
  uint32_t mPadID;         /// PadIndex to which the digit corresponds to
  uint64_t mADC; /// Amplitude of signal
  uint8_t mCathode;
  float mX, mY, mSizeX, mSizeY;
}; //class Digit2



class PreClusterSinkTask
{
public:
  //_________________________________________________________________________________________________
  void init(framework::InitContext& ic)
  {
    /// Get the output file from the context
    LOG(INFO) << "initializing precluster sink";

    mText = ic.options().get<bool>("txt");

    auto outputFileName = ic.options().get<std::string>("outfile");
    mOutputFile.open(outputFileName, (mText ? ios::out : (ios::out | ios::binary)));
    if (!mOutputFile.is_open()) {
      throw invalid_argument("Cannot open output file" + outputFileName);
    }

    mUseRun2DigitUID = ic.options().get<bool>("useRun2DigitUID");

    auto stop = [this]() {
      /// close the output file
      LOG(INFO) << "stop precluster sink";
      this->mOutputFile.close();
    };
    ic.services().get<CallbackService>().set(CallbackService::Id::Stop, stop);
  }

  //_________________________________________________________________________________________________
  void run(framework::ProcessingContext& pc)
  {
    /// dump the preclusters with associated digits of the current event

    // get the input preclusters and associated digits
    auto preClusters = pc.inputs().get<gsl::span<PreCluster>>("preclusters");
    auto digits = pc.inputs().get<gsl::span<Digit>>("digits");
    Clustering clustering;

    if (mText) {
      mOutputFile << preClusters.size() << " preclusters:" << endl;
      if (mUseRun2DigitUID) {
        std::vector<Digit> digitsCopy(digits.begin(), digits.end());
        convertPadID2DigitUID(digitsCopy);
        for (const auto& precluster : preClusters) {
          precluster.print(mOutputFile, digitsCopy);
        }
      } else {
        for (const auto& precluster : preClusters) {
          precluster.print(mOutputFile, digits);
        }
      }
    } else {
      float xtrk = 0, ytrk = 0;
      mOutputFile.write(reinterpret_cast<char*>(&xtrk), sizeof(float));
      mOutputFile.write(reinterpret_cast<char*>(&ytrk), sizeof(float));

      // write the number of preclusters
      int nPreClusters = preClusters.size();
      mOutputFile.write(reinterpret_cast<char*>(&nPreClusters), sizeof(int));
      cout << "nPreclusters: " << preClusters.size() << endl;

      for (auto& preCluster : preClusters) {
        gsl::span<const PreCluster> sPreCluster = {&preCluster, 1};
        //std::vector<Clustering::Cluster> clusters(0);
        // Fit Mathieson
        //clustering.runFinderSimpleFit(sPreCluster, digits, clusters);
        //float clusPosition[2] = {
        //    clusters.empty() ? 0 : clusters[0].getx(),
        //        clusters.empty() ? 0 : clusters[0].gety()
        //};
        //mOutputFile.write(reinterpret_cast<char*>(clusPosition), sizeof(clusPosition));

        auto preClusterDigits = digits.subspan(preCluster.firstDigit, preCluster.nDigits);
        std::vector<Digit> dv(preClusterDigits.begin(), preClusterDigits.end());

        // write the total number of digits in this precluster
        int nDigits = dv.size();
        mOutputFile.write(reinterpret_cast<char*>(&nDigits), sizeof(int));
        //cout<<"sizeof(Digit2): "<<sizeof(Digit2)<<endl;

        for (auto& digit : dv) {
          if(false){
            cout << "\nDetID:" << digit.getDetID() << " PadID:" << digit.getPadID() << endl;
          }

          Digit2 d2;
          d2.mTime = static_cast<uint64_t>(digit.getTimeStamp());
          d2.mDetID = static_cast<uint32_t>(digit.getDetID());
          d2.mPadID = static_cast<uint32_t>(digit.getPadID());
          d2.mADC = static_cast<uint64_t>(digit.getADC());
          const mapping::Segmentation& segment = mapping::segmentation(digit.getDetID());

          d2.mCathode = (segment.isBendingPad(digit.getPadID()) ? 0 : 1);
          d2.mX = segment.padPositionX(digit.getPadID());
          d2.mY = segment.padPositionY(digit.getPadID());
          d2.mSizeX = segment.padSizeX(digit.getPadID());
          d2.mSizeY = segment.padSizeY(digit.getPadID());

          mOutputFile.write(reinterpret_cast<char*>(&d2), sizeof(Digit2));
        }
      }

      // write the number of clusters
      //int nClusters = 0;
      //mOutputFile.write(reinterpret_cast<char*>(&nClusters), sizeof(int));
    }
  }

private:
  //_________________________________________________________________________________________________
  void convertPadID2DigitUID(std::vector<Digit>& digits)
  {
    /// convert the pad ID (i.e. index) in O2 mapping into a digit UID in run2 format

    // cathode number of the bending plane for each DE
    static const std::array<std::vector<int>, 10> bendingCathodes{
      {{0, 1, 0, 1},
        {0, 1, 0, 1},
        {0, 1, 0, 1},
        {0, 1, 0, 1},
        {0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1},
        {0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1},
        {0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1},
        {0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1},
        {0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1},
        {0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1}}};

    for (auto& digit : digits) {

      int deID = digit.getDetID();
      auto& segmentation = mapping::segmentation(deID);
      int padID = digit.getPadID();
      int cathode = bendingCathodes[deID / 100 - 1][deID % 100];
      if (!segmentation.isBendingPad(padID)) {
        cathode = 1 - cathode;
      }
      int manuID = segmentation.padDualSampaId(padID);
      int manuCh = segmentation.padDualSampaChannel(padID);

      int digitID = (deID) | (manuID << 12) | (manuCh << 24) | (cathode << 30);
      digit.setPadID(digitID);
    }
  }

  std::ofstream mOutputFile{};   ///< output file
  bool mText = false;            ///< output preclusters in text format
  bool mUseRun2DigitUID = false; ///< true if Digit.mPadID = digit UID in run2 format
};

//_________________________________________________________________________________________________
o2::framework::DataProcessorSpec getPreClusterSinkSpec()
{
  return DataProcessorSpec{
    "PreClusterSink",
    Inputs{InputSpec{"preclusters", "MCH", "PRECLUSTERS", 0, Lifetime::Timeframe},
      InputSpec{"digits", "MCH", "PRECLUSTERDIGITS", 0, Lifetime::Timeframe}},
      Outputs{},
      AlgorithmSpec{adaptFromTask<PreClusterSinkTask>()},
      Options{{"outfile", VariantType::String, "preclusters.out", {"output filename"}},
        {"txt", VariantType::Bool, false, {"output preclusters in text format"}},
        {"useRun2DigitUID", VariantType::Bool, false, {"mPadID = digit UID in run2 format"}}}};
}

} // end namespace mch
} // end namespace o2
