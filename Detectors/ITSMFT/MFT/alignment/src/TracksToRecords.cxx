// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// @file TracksToRecords.cxx

#include <TTreeReader.h>
#include <TTreeReaderValue.h>

#include "Framework/InputSpec.h"
#include "Framework/Logger.h"
#include <Framework/InputRecord.h>
#include "MFTBase/Geometry.h"
#include <MFTTracking/Constants.h>
#include "ForwardAlign/MillePedeRecord.h"

#include "ReconstructionDataFormats/PrimaryVertex.h"

#include "MFTAlignment/TracksToRecords.h"

using namespace o2::mft;

ClassImp(o2::mft::TracksToRecords);


static void FwdtoMCH(const o2::track::TrackParCovFwd& fwdtrack, std::array<double, 5>& mchPar, std::array<double, 15>& mchCov)
{
  using SMatrix55Std = ROOT::Math::SMatrix<double, 5>;
  using SMatrix55Sym = ROOT::Math::SMatrix<double, 5, 5, ROOT::Math::MatRepSym<double, 5>>;

  // Convert Forward Track parameters and covariances matrix to the MCH track format.

  // Parameter conversion
  double alpha1, alpha3, alpha4, x2, x3, x4;

  x2 = fwdtrack.getPhi();
  x3 = fwdtrack.getTanl();
  x4 = fwdtrack.getInvQPt();

  auto sinx2 = TMath::Sin(x2);
  auto cosx2 = TMath::Cos(x2);

  alpha1 = cosx2 / x3;
  alpha3 = sinx2 / x3;
  alpha4 = x4 / TMath::Sqrt(x3 * x3 + sinx2 * sinx2);

  auto K = TMath::Sqrt(x3 * x3 + sinx2 * sinx2);
  auto K3 = K * K * K;

  // Covariances matrix conversion
  SMatrix55Std jacobian;
  SMatrix55Sym covariances;

  covariances(0, 0) = fwdtrack.getCovariances()(0, 0);
  covariances(0, 1) = fwdtrack.getCovariances()(0, 1);
  covariances(0, 2) = fwdtrack.getCovariances()(0, 2);
  covariances(0, 3) = fwdtrack.getCovariances()(0, 3);
  covariances(0, 4) = fwdtrack.getCovariances()(0, 4);

  covariances(1, 1) = fwdtrack.getCovariances()(1, 1);
  covariances(1, 2) = fwdtrack.getCovariances()(1, 2);
  covariances(1, 3) = fwdtrack.getCovariances()(1, 3);
  covariances(1, 4) = fwdtrack.getCovariances()(1, 4);

  covariances(2, 2) = fwdtrack.getCovariances()(2, 2);
  covariances(2, 3) = fwdtrack.getCovariances()(2, 3);
  covariances(2, 4) = fwdtrack.getCovariances()(2, 4);

  covariances(3, 3) = fwdtrack.getCovariances()(3, 3);
  covariances(3, 4) = fwdtrack.getCovariances()(3, 4);

  covariances(4, 4) = fwdtrack.getCovariances()(4, 4);

  jacobian(0, 0) = 1;

  jacobian(1, 2) = -sinx2 / x3;
  jacobian(1, 3) = -cosx2 / (x3 * x3);

  jacobian(2, 1) = 1;

  jacobian(3, 2) = cosx2 / x3;
  jacobian(3, 3) = -sinx2 / (x3 * x3);

  jacobian(4, 2) = -x4 * sinx2 * cosx2 / K3;
  jacobian(4, 3) = -x3 * x4 / K3;
  jacobian(4, 4) = 1 / K;
  // jacobian*covariances*jacobian^T
  covariances = ROOT::Math::Similarity(jacobian, covariances);

  //mchCov = {covariances(0, 0), covariances(1, 0), covariances(1, 1), covariances(2, 0), covariances(2, 1), covariances(2, 2), covariances(3, 0), covariances(3, 1), covariances(3, 2), covariances(3, 3), covariances(4, 0), covariances(4, 1), covariances(4, 2), covariances(4, 3), covariances(4, 4)};
  mchCov = {
      covariances(0, 0),
      covariances(1, 0),
      covariances(2, 0),
      covariances(3, 0),
      covariances(4, 0),
      covariances(1, 1),
      covariances(2, 1),
      covariances(3, 1),
      covariances(4, 1),
      covariances(2, 2),
      covariances(2, 3),
      covariances(2, 4),
      covariances(3, 3),
      covariances(3, 4),
      covariances(4, 4)};
  mchPar = {fwdtrack.getX(), alpha1, fwdtrack.getY(), alpha3, alpha4};
}

static void MCHtoFwd(o2::track::TrackParCovFwd& fwdtrack, std::array<double, 5>& mchPar, std::array<double, 15>& mchCov)
{
  using SMatrix55Std = ROOT::Math::SMatrix<double, 5>;
  using SMatrix55Sym = ROOT::Math::SMatrix<double, 5, 5, ROOT::Math::MatRepSym<double, 5>>;

  // Convert a MCH Track parameters and covariances matrix to the
  // Forward track format. Must be called after propagation though the absorber

  // Parameter conversion
  double alpha1, alpha3, alpha4, x2, x3, x4;

  alpha1 = mchPar[1];
  alpha3 = mchPar[3];
  alpha4 = mchPar[4];

  x2 = TMath::ATan2(-alpha3, -alpha1);
  x3 = -1. / TMath::Sqrt(alpha3 * alpha3 + alpha1 * alpha1);
  x4 = alpha4 * -x3 * TMath::Sqrt(1 + alpha3 * alpha3);

  auto K = alpha1 * alpha1 + alpha3 * alpha3;
  auto K32 = K * TMath::Sqrt(K);
  auto L = TMath::Sqrt(alpha3 * alpha3 + 1);

  // Covariances matrix conversion
  SMatrix55Std jacobian;
  SMatrix55Sym covariances;

  covariances(0, 0) = mchCov[0]; //mchParam.getCovariances()(0, 0);
  covariances(0, 1) = mchCov[1]; //mchParam.getCovariances()(0, 1);
  covariances(0, 2) = mchCov[2]; //mchParam.getCovariances()(0, 2);
  covariances(0, 3) = mchCov[3]; //mchParam.getCovariances()(0, 3);
  covariances(0, 4) = mchCov[4]; //mchParam.getCovariances()(0, 4);

  covariances(1, 1) = mchCov[5]; //mchParam.getCovariances()(1, 1);
  covariances(1, 2) = mchCov[6]; //mchParam.getCovariances()(1, 2);
  covariances(1, 3) = mchCov[7]; //mchParam.getCovariances()(1, 3);
  covariances(1, 4) = mchCov[8]; //mchParam.getCovariances()(1, 4);

  covariances(2, 2) = mchCov[9]; //mchParam.getCovariances()(2, 2);
  covariances(2, 3) = mchCov[10]; //mchParam.getCovariances()(2, 3);
  covariances(2, 4) = mchCov[11]; //mchParam.getCovariances()(2, 4);

  covariances(3, 3) = mchCov[12]; //mchParam.getCovariances()(3, 3);
  covariances(3, 4) = mchCov[13]; //mchParam.getCovariances()(3, 4);

  covariances(4, 4) = mchCov[14]; //mchParam.getCovariances()(4, 4);

  jacobian(0, 0) = 1;

  jacobian(1, 2) = 1;

  jacobian(2, 1) = -alpha3 / K;
  jacobian(2, 3) = alpha1 / K;

  jacobian(3, 1) = alpha1 / K32;
  jacobian(3, 3) = alpha3 / K32;

  jacobian(4, 1) = -alpha1 * alpha4 * L / K32;
  jacobian(4, 3) = alpha3 * alpha4 * (1 / (TMath::Sqrt(K) * L) - L / K32);
  jacobian(4, 4) = L / TMath::Sqrt(K);

  // jacobian*covariances*jacobian^T
  covariances = ROOT::Math::Similarity(jacobian, covariances);

  // Set output
  //fwdtrack.setX(mchParam.getNonBendingCoor());
  //fwdtrack.setY(mchParam.getBendingCoor());
  //fwdtrack.setZ(mchParam.getZ());
  fwdtrack.setPhi(x2);
  fwdtrack.setTanl(x3);
  fwdtrack.setInvQPt(x4);
  //fwdtrack.setCharge(mchParam.getCharge());
  fwdtrack.setCovariances(covariances);
}

//__________________________________________________________________________
TracksToRecords::TracksToRecords()
  : Aligner(),
    mRunNumber(0),
    mBz(0),
    mNumberTFs(0),
    mNumberOfClusterChainROFs(0),
    mNumberOfTrackChainROFs(0),
    mCounterLocalEquationFailed(0),
    mCounterSkippedTracks(0),
    mCounterUsedTracks(0),
    mGlobalDerivatives(std::vector<double>(mNumberOfGlobalParam)),
    mLocalDerivatives(std::vector<double>(mNumberOfTrackParam)),
    mMinNumberClusterCut(6),
    mWeightRecord(1.),
    mDictionary(nullptr),
    mAlignPoint(new AlignPointHelper()),
    mWithControl(false),
    mNEntriesAutoSave(10000),
    mRecordWriter(new o2::fwdalign::MilleRecordWriter()),
    mWithConstraintsRecWriter(false),
    mConstraintsRecWriter(nullptr),
    mMillepede(new o2::fwdalign::MillePede2())
{
  if (mWithConstraintsRecWriter) {
    mConstraintsRecWriter = new o2::fwdalign::MilleRecordWriter();
  }
  // initialise the content of each array
  resetGlocalDerivative();
  resetLocalDerivative();
  LOGF(debug, "TracksToRecords instantiated");
}

//__________________________________________________________________________
TracksToRecords::~TracksToRecords()
{
  if (mConstraintsRecWriter) {
    delete mConstraintsRecWriter;
  }
  if (mMillepede) {
    delete mMillepede;
  }
  if (mRecordWriter) {
    delete mRecordWriter;
  }
  if (mAlignPoint) {
    delete mAlignPoint;
  }
  if (mDictionary) {
    mDictionary = nullptr;
  }
  LOGF(debug, "TracksToRecords destroyed");
}

//__________________________________________________________________________
void TracksToRecords::init()
{
  if (mIsInitDone) {
    return;
  }
  if (mDictionary == nullptr) {
    LOGF(fatal, "TracksToRecords::init() failed because no cluster dictionary is defined");
    mIsInitDone = false;
    return;
  }

  mRecordWriter->setCyclicAutoSave(mNEntriesAutoSave);
  mRecordWriter->setDataFileName(mMilleRecordsFileName);
  mMillepede->SetRecordWriter(mRecordWriter);

  if (mWithConstraintsRecWriter) {
    mConstraintsRecWriter->setCyclicAutoSave(mNEntriesAutoSave);
    mConstraintsRecWriter->setDataFileName(mMilleConstraintsRecFileName);
    mMillepede->SetConstraintsRecWriter(mConstraintsRecWriter);
  }

  mAlignPoint->setClusterDictionary(mDictionary);

  mMillepede->InitMille(mNumberOfGlobalParam,
                        mNumberOfTrackParam,
                        mChi2CutNStdDev,
                        mResCut,
                        mResCutInitial);

  LOG(info) << "-------------- TracksToRecords configured with -----------------";
  LOGF(info, "Chi2CutNStdDev = %d", mChi2CutNStdDev);
  LOGF(info, "ResidualCutInitial = %.3f", mResCutInitial);
  LOGF(info, "ResidualCut = %.3f", mResCut);
  LOGF(info, "MinNumberClusterCut = %d", mMinNumberClusterCut);
  LOGF(info, "mStartFac = %.3f", mStartFac);
  LOGF(info,
       "Allowed variation: dx = %.3f, dy = %.3f, dz = %.3f, dRz = %.4f",
       mAllowVar[0], mAllowVar[1], mAllowVar[3], mAllowVar[2]);
  LOG(info) << "-----------------------------------------------------------";

  // set allowed variations for all parameters
  for (int chipId = 0; chipId < mNumberOfSensors; ++chipId) {
    for (Int_t iPar = 0; iPar < mNDofPerSensor; ++iPar) {
      mMillepede->SetParSigma(chipId * mNDofPerSensor + iPar, mAllowVar[iPar]);
    }
  }

  // set iterations
  if (mStartFac > 1) {
    mMillepede->SetIterations(mStartFac);
  }

  mIsInitDone = true;
  LOGF(info, "TracksToRecords init done");
}

//__________________________________________________________________________
void TracksToRecords::processTimeFrame(o2::framework::ProcessingContext& ctx)
{
  mNumberTFs++; // TF Counter

  // get tracks
  mMFTTracks = ctx.inputs().get<gsl::span<o2::mft::TrackMFT>>("tracks");
  mMFTTracksROF = ctx.inputs().get<gsl::span<o2::itsmft::ROFRecord>>("tracksrofs");
  mMFTTrackClusIdx = ctx.inputs().get<gsl::span<int>>("trackClIdx");

  // get clusters
  mMFTClusters = ctx.inputs().get<gsl::span<o2::itsmft::CompClusterExt>>("compClusters");
  mMFTClustersROF = ctx.inputs().get<gsl::span<o2::itsmft::ROFRecord>>("clustersrofs");
  mMFTClusterPatterns = ctx.inputs().get<gsl::span<unsigned char>>("patterns");
  mPattIt = mMFTClusterPatterns.begin();
  mAlignPoint->convertCompactClusters(
    mMFTClusters, mPattIt, mMFTClustersLocal, mMFTClustersGlobal);
}

//__________________________________________________________________________
void TracksToRecords::processRecoTracks()
{
  if (!mIsInitDone) {
    LOGF(fatal, "TracksToRecords::processRecoTracks() aborted because init was not done !");
    return;
  }
  if (!mRecordWriter || !mRecordWriter->isInitOk()) {
    LOGF(fatal, "TracksToRecords::processRecoTracks() aborted because uninitialised mRecordWriter !");
    return;
  }

  LOG(info) << "TracksToRecords::processRecoTracks() - start";

  int nCounterAllTracks = 0;

  for (auto& oneTrack : mMFTTracks) { // track loop

    LOGF(debug, "Processing track # %5d", nCounterAllTracks);

    // Skip the track if not enough clusters
    auto ncls = oneTrack.getNumberOfPoints();
    if (ncls < mMinNumberClusterCut) {
      nCounterAllTracks++;
      mCounterSkippedTracks++;
      continue;
    }

    // Skip presumably quite low momentum track
    if (!oneTrack.isLTF()) {
      nCounterAllTracks++;
      mCounterSkippedTracks++;
      continue;
    }

    auto offset = oneTrack.getExternalClusterIndexOffset();

    mRecordWriter->getRecord()->Reset();

    // Store the initial track parameters
    auto track = oneTrack;
    mAlignPoint->resetTrackInitialParam();
    mAlignPoint->recordTrackInitialParam(track);

    bool isTrackUsed = true;

    for (int icls = 0; icls < ncls; ++icls) { // cluster loop

      mAlignPoint->resetAlignPoint();

      // Store measured positions
      auto clsEntry = mMFTTrackClusIdx[offset + icls];
      auto localCluster = mMFTClustersLocal[clsEntry];
      auto globalCluster = mMFTClustersGlobal[clsEntry];
      mAlignPoint->setMeasuredPosition(localCluster, globalCluster);
      if (!mAlignPoint->isClusterOk()) {
        LOGF(warning, "TracksToRecords::processRecoTracks() - will not use track # %5d with at least a bad cluster", nCounterAllTracks);
        mCounterSkippedTracks++;
        isTrackUsed = false;
        break;
      }

      // Propagate track to the current z plane of this cluster
      track.propagateParamToZlinear(mAlignPoint->getGlobalMeasuredPosition().Z());

      // Store reco positions
      mAlignPoint->setGlobalRecoPosition(track);

      // compute residuals
      mAlignPoint->setLocalResidual();
      mAlignPoint->setGlobalResidual();

      // Compute derivatives
      mAlignPoint->computeLocalDerivatives();
      mAlignPoint->computeGlobalDerivatives();

      // Set local equations
      bool success = true;
      success &= setLocalEquationX();
      success &= setLocalEquationY();
      success &= setLocalEquationZ();
      isTrackUsed &= success;
      if (mWithControl && success) {
        mPointControl.fill(mAlignPoint, mCounterUsedTracks);
      }
      if (!success) {
        LOGF(error, "TracksToRecords::processRecoTracks() - track %i h %d d %d l %d s %4d lMpos x %.2e y %.2e z %.2e gMpos x %.2e y %.2e z %.2e gRpos x %.2e y %.2e z %.2e",
             mCounterUsedTracks, mAlignPoint->half(), mAlignPoint->disk(), mAlignPoint->layer(), mAlignPoint->getSensorId(),
             mAlignPoint->getLocalMeasuredPosition().X(), mAlignPoint->getLocalMeasuredPosition().Y(), mAlignPoint->getLocalMeasuredPosition().Z(),
             mAlignPoint->getGlobalMeasuredPosition().X(), mAlignPoint->getGlobalMeasuredPosition().Y(), mAlignPoint->getGlobalMeasuredPosition().Z(),
             mAlignPoint->getGlobalRecoPosition().X(), mAlignPoint->getGlobalRecoPosition().Y(), mAlignPoint->getGlobalRecoPosition().Z());
      }

    } // end of loop on clusters

    if (isTrackUsed) {
      mRecordWriter->setRecordRun(mRunNumber);
      mRecordWriter->setRecordWeight(mWeightRecord);
      const bool doPrint = false;
      mRecordWriter->fillRecordTree(doPrint); // save record data
      mCounterUsedTracks++;
    }
    nCounterAllTracks++;
  } // end of loop on tracks
  LOG(info) << "TracksToRecords::processRecoTracks() - end";
}

//__________________________________________________________________________
void TracksToRecords::processROFs(TChain* mfttrackChain, TChain* mftclusterChain)
{
  if (!mIsInitDone) {
    LOGF(fatal, "TracksToRecords::processROFs() aborted because init was not done !");
    return;
  }

  if (!mRecordWriter || !mRecordWriter->isInitOk()) {
    LOGF(fatal, "TracksToRecords::processROFs() aborted because uninitialised mRecordWriter !");
    return;
  }

  LOG(info) << "TracksToRecords::processROFs() - start";

  TTreeReader mftTrackChainReader(mfttrackChain);
  TTreeReader mftClusterChainReader(mftclusterChain);
  std::vector<unsigned char>::iterator pattIterator;

  TTreeReaderValue<std::vector<o2::mft::TrackMFT>> mftTracks =
    {mftTrackChainReader, "MFTTrack"};
  TTreeReaderValue<std::vector<o2::itsmft::ROFRecord>> mftTracksROF =
    {mftTrackChainReader, "MFTTracksROF"};
  TTreeReaderValue<std::vector<int>> mftTrackClusIdx =
    {mftTrackChainReader, "MFTTrackClusIdx"};

  TTreeReaderValue<std::vector<o2::itsmft::CompClusterExt>> mftClusters =
    {mftClusterChainReader, "MFTClusterComp"};
  TTreeReaderValue<std::vector<o2::itsmft::ROFRecord>> mftClustersROF =
    {mftClusterChainReader, "MFTClustersROF"};
  TTreeReaderValue<std::vector<unsigned char>> mftClusterPatterns =
    {mftClusterChainReader, "MFTClusterPatt"};

  int nCounterAllTracks = 0;

  while (mftTrackChainReader.Next() && mftClusterChainReader.Next()) {

    mNumberOfTrackChainROFs += (*mftTracksROF).size();
    mNumberOfClusterChainROFs += (*mftClustersROF).size();
    assert(mNumberOfTrackChainROFs == mNumberOfClusterChainROFs);

    pattIterator = (*mftClusterPatterns).begin();
    mAlignPoint->convertCompactClusters(
      *mftClusters, pattIterator, mMFTClustersLocal, mMFTClustersGlobal);

    //______________________________________________________
    for (auto& oneTrack : *mftTracks) { // track loop

      LOGF(debug, "Processing track # %5d", nCounterAllTracks);

      // Skip the track if not enough clusters
      auto ncls = oneTrack.getNumberOfPoints();
      if (ncls < mMinNumberClusterCut) {
        nCounterAllTracks++;
        mCounterSkippedTracks++;
        continue;
      }

      // Skip presumably quite low momentum track
      if (!oneTrack.isLTF()) {
        nCounterAllTracks++;
        mCounterSkippedTracks++;
        continue;
      }

      auto offset = oneTrack.getExternalClusterIndexOffset();

      mRecordWriter->getRecord()->Reset();

      // Store the initial track parameters
      mAlignPoint->resetTrackInitialParam();
      mAlignPoint->recordTrackInitialParam(oneTrack);

      bool isTrackUsed = true;

      for (int icls = 0; icls < ncls; ++icls) { // cluster loop

        mAlignPoint->resetAlignPoint();

        // Store measured positions
        auto clsEntry = (*mftTrackClusIdx)[offset + icls];
        auto localCluster = mMFTClustersLocal[clsEntry];
        auto globalCluster = mMFTClustersGlobal[clsEntry];
        mAlignPoint->setMeasuredPosition(localCluster, globalCluster);
        if (!mAlignPoint->isClusterOk()) {
          LOGF(warning, "TracksToRecords::processROFs() - will not use track # %5d with at least a bad cluster", nCounterAllTracks);
          mCounterSkippedTracks++;
          isTrackUsed = false;
          break;
        }

        // Propagate track to the current z plane of this cluster
        oneTrack.propagateParamToZlinear(mAlignPoint->getGlobalMeasuredPosition().Z());

        // Store reco positions
        mAlignPoint->setGlobalRecoPosition(oneTrack);

        // compute residuals
        mAlignPoint->setLocalResidual();
        mAlignPoint->setGlobalResidual();

        // Compute derivatives
        mAlignPoint->computeLocalDerivatives();
        mAlignPoint->computeGlobalDerivatives();

        // Set local equations
        bool success = true;
        success &= setLocalEquationX();
        success &= setLocalEquationY();
        success &= setLocalEquationZ();
        isTrackUsed &= success;
        if (mWithControl && success) {
          mPointControl.fill(mAlignPoint, mCounterUsedTracks);
        }
        if (!success) {
          LOGF(error, "TracksToRecords::processROFs() - track %i h %d d %d l %d s %4d lMpos x %.2e y %.2e z %.2e gMpos x %.2e y %.2e z %.2e gRpos x %.2e y %.2e z %.2e",
               mCounterUsedTracks, mAlignPoint->half(), mAlignPoint->disk(), mAlignPoint->layer(), mAlignPoint->getSensorId(),
               mAlignPoint->getLocalMeasuredPosition().X(), mAlignPoint->getLocalMeasuredPosition().Y(), mAlignPoint->getLocalMeasuredPosition().Z(),
               mAlignPoint->getGlobalMeasuredPosition().X(), mAlignPoint->getGlobalMeasuredPosition().Y(), mAlignPoint->getGlobalMeasuredPosition().Z(),
               mAlignPoint->getGlobalRecoPosition().X(), mAlignPoint->getGlobalRecoPosition().Y(), mAlignPoint->getGlobalRecoPosition().Z());
        }

      } // end of loop on clusters

      if (isTrackUsed) {
        // copy track record
        mRecordWriter->setRecordRun(mRunNumber);
        mRecordWriter->setRecordWeight(mWeightRecord);
        const bool doPrint = false;
        mRecordWriter->fillRecordTree(doPrint); // save record data
        mCounterUsedTracks++;
      }
      nCounterAllTracks++;
    } // end of loop on tracks

  } // end of loop on TChain reader

  LOG(info) << "TracksToRecords::processROFs() - end";
}

//__________________________________________________________________________
void TracksToRecords::processROFs(TChain* itsvertexChain, TChain* mfttrackChain, TChain* mftclusterChain)
{
  //const float vertexZShift = 0.225;
  const float vertexZShift = 0.25;

  if (!mIsInitDone) {
    LOGF(fatal, "TracksToRecords::processROFs() aborted because init was not done !");
    return;
  }

  if (!mRecordWriter || !mRecordWriter->isInitOk()) {
    LOGF(fatal, "TracksToRecords::processROFs() aborted because uninitialised mRecordWriter !");
    return;
  }

  LOG(info) << "TracksToRecords::processROFs() with ITS vertex - start";

  TTreeReader vertexChainReader(itsvertexChain);
  TTreeReader mftTrackChainReader(mfttrackChain);
  TTreeReader mftClusterChainReader(mftclusterChain);
  std::vector<unsigned char>::iterator pattIterator;

  TTreeReaderValue<std::vector<o2::dataformats::PrimaryVertex>> vertices =
    {vertexChainReader, "PrimaryVertex"};

  TTreeReaderValue<std::vector<o2::mft::TrackMFT>> mftTracks =
    {mftTrackChainReader, "MFTTrack"};
  TTreeReaderValue<std::vector<o2::itsmft::ROFRecord>> mftTracksROF =
    {mftTrackChainReader, "MFTTracksROF"};
  TTreeReaderValue<std::vector<int>> mftTrackClusIdx =
    {mftTrackChainReader, "MFTTrackClusIdx"};

  TTreeReaderValue<std::vector<o2::itsmft::CompClusterExt>> mftClusters =
    {mftClusterChainReader, "MFTClusterComp"};
  TTreeReaderValue<std::vector<o2::itsmft::ROFRecord>> mftClustersROF =
    {mftClusterChainReader, "MFTClustersROF"};
  TTreeReaderValue<std::vector<unsigned char>> mftClusterPatterns =
    {mftClusterChainReader, "MFTClusterPatt"};

  int nCounterAllTracks = 0;

  size_t entry = 0;
  auto nEntries = mftTrackChainReader.GetEntries(kTRUE);
  size_t nMFTtracks = 0;
  while (vertexChainReader.Next() && mftTrackChainReader.Next() && mftClusterChainReader.Next()) {

    entry +=1;
    if ((entry % 100) == 0) {
      std::cout << std::format("Reading entry {} / {}", entry, nEntries) << std::endl;
    }

    mNumberOfTrackChainROFs += (*mftTracksROF).size();
    mNumberOfClusterChainROFs += (*mftClustersROF).size();
    assert(mNumberOfTrackChainROFs == mNumberOfClusterChainROFs);

    pattIterator = (*mftClusterPatterns).begin();
    mAlignPoint->convertCompactClusters(
      *mftClusters, pattIterator, mMFTClustersLocal, mMFTClustersGlobal);

    nMFTtracks += mftTracks->size();
    if (mftTracks->size() > 0) {
    //  std::cout << std::format("Entry #{} with {} MFT tracks (total = {})", entry-1, mftTracks->size(), nMFTtracks) << std::endl;
    }

    //std::cout << std::format("Number of ITS vertices: {}", (*vertices).size()) << std::endl;

    //______________________________________________________
    for (const auto& oneRof : *mftTracksROF) { // track ROF loop
      const auto& rofStart = oneRof.getBCData();
      auto rofEnd = rofStart + 198;

      int vertexId = -1;
      int vertexCount = -1;
      int nVerticesInRof = 0;
      //______________________________________________________
      for (const auto& oneVertex : *vertices) { // vertex loop
        vertexCount += 1;
        if (oneVertex.getIRMin() < rofStart) {
          continue;
        }
        if (oneVertex.getIRMax() > rofEnd) {
          continue;
        }
        nVerticesInRof += 1;
        vertexId = vertexCount;

        if (false) {
          std::cout << std::format("[TOTO] correlated vertex found with IR=({},{} -> {},{})  MFT IR=({},{} -> {},{})",
              oneVertex.getIRMin().orbit, oneVertex.getIRMin().bc,
              oneVertex.getIRMax().orbit, oneVertex.getIRMax().bc,
              rofStart.orbit, rofStart.bc,
              rofEnd.orbit, rofEnd.bc) << std::endl;
        }
      }
      if (false && nVerticesInRof > 0) {
        std::cout << std::format("[TOTO] number of vertices in MFT ROF: {}", nVerticesInRof) << std::endl;
        std::cout << std::format("[TOTO] selected vertex ID: {}  ", vertexId) << std::endl;
      }

      // only select ROFs with a single ITS vertex
      if (nVerticesInRof != 1) {
        continue;
      }

      const auto& theVertex = (*vertices)[vertexId];

      // discard vertices in the +/- 5cm region
      auto pvx = theVertex.getX();
      auto pvy = theVertex.getY();
      auto pvz = theVertex.getZ();
      if (std::fabs(pvz) < 5) {
        continue;
      }
      //std::cout << std::format("[TOTO] vertex z={:0.3f} IR=({},{} -> {},{})  MFT IR=({},{} -> {},{})",
      //    pvz,
      //    theVertex.getIRMin().orbit, theVertex.getIRMin().bc,
      //    theVertex.getIRMax().orbit, theVertex.getIRMax().bc,
      //    rofStart.orbit, rofStart.bc,
      //    rofEnd.orbit, rofEnd.bc) << std::endl;

      pvz -= vertexZShift;;

      std::array<double, 3> pv{ pvx, pvy, pvz };

      int firstTrackIndex = oneRof.getFirstEntry();
      int lastTrackIndex = oneRof.getFirstEntry() + oneRof.getNEntries() - 1;
      for (int iTrack = firstTrackIndex; iTrack <= lastTrackIndex; iTrack++) {
        auto oneTrack = (*mftTracks)[iTrack];

        // Skip the track if not enough clusters
        auto ncls = oneTrack.getNumberOfPoints();
        if (ncls < mMinNumberClusterCut) {
          //std::cout << "Track has " << ncls << " clusters -> skipped" << std::endl;
          nCounterAllTracks++;
          mCounterSkippedTracks++;
          continue;
        }

        // Skip presumably quite low momentum track
        if (!oneTrack.isLTF()) {
          //std::cout << "Track is not LTF -> skipped" << std::endl;
          nCounterAllTracks++;
          mCounterSkippedTracks++;
          continue;
        }

        // rotate & shift the track such that it goeas through the primary vertex
        // the rotation occurs at the middle plane of the MFT
        double zRefPlane = (o2::mft::constants::mft::LayerZCoordinate()[9] + o2::mft::constants::mft::LayerZCoordinate()[0]) / 2.f;

        // extrapolate track to reference plane
        oneTrack.propagateParamToZlinear(zRefPlane);

        //std::cout << std::format("MFT track at vertex: z={:0.3f}  dx={:0.3f}  dy={:0.3f}", pvz, (oneTrack.getX() - pvx), (oneTrack.getY() - pvy)) << std::endl;

        // compute new track slopes
        double dz = zRefPlane - pvz;
        double newSlopeX = (oneTrack.getX() - pvx) / dz;
        double newSlopeY = (oneTrack.getY() - pvy) / dz;

        // modify current track
        std::array<double, 5> mchPar{ 0.0 };
        std::array<double, 15> mchCov{ 0.0 };
        FwdtoMCH(oneTrack, mchPar, mchCov);
        mchPar[1] = newSlopeX;
        mchPar[3] = newSlopeY;
        MCHtoFwd(oneTrack, mchPar, mchCov);

        // get the track parameters at the primary vertex
        oneTrack.propagateParamToZlinear(pvz);
        //std::cout << std::format("Vertex:              x={:0.3f} y={:0.3f} z={:0.3f}", pvx, pvy, pvz) << std::endl;
        //std::cout << std::format("MFT track at vertex: x={:0.3f} y={:0.3f} z={:0.3f}",
        //    oneTrack.getX(), oneTrack.getY(), oneTrack.getZ()) << std::endl;

        auto offset = oneTrack.getExternalClusterIndexOffset();

        mRecordWriter->getRecord()->Reset();

        // Store the initial track parameters
        mAlignPoint->resetTrackInitialParam();
        mAlignPoint->recordTrackInitialParam(oneTrack);

        bool isTrackUsed = true;

        for (int icls = 0; icls < ncls; ++icls) { // cluster loop

          mAlignPoint->resetAlignPoint();

          // Store measured positions
          auto clsEntry = (*mftTrackClusIdx)[offset + icls];
          auto localCluster = mMFTClustersLocal[clsEntry];
          auto globalCluster = mMFTClustersGlobal[clsEntry];
          mAlignPoint->setMeasuredPosition(localCluster, globalCluster);
          if (!mAlignPoint->isClusterOk()) {
            LOGF(warning, "TracksToRecords::processROFs() - will not use track # %5d with at least a bad cluster", nCounterAllTracks);
            //std::cout << "Track has at least one bad cluster -> skipped" << std::endl;
            mCounterSkippedTracks++;
            isTrackUsed = false;
            break;
          }

          // Propagate track to the current z plane of this cluster
          oneTrack.propagateParamToZlinear(mAlignPoint->getGlobalMeasuredPosition().Z());

          // Store reco positions
          mAlignPoint->setGlobalRecoPosition(oneTrack);

          // compute residuals
          mAlignPoint->setLocalResidual();
          mAlignPoint->setGlobalResidual();

          // Compute derivatives
          mAlignPoint->computeLocalDerivatives();
          mAlignPoint->computeGlobalDerivatives();

          // Set local equations
          bool success = true;
          success &= setLocalEquationX();
          success &= setLocalEquationY();
          success &= setLocalEquationZ();
          isTrackUsed &= success;
          if (mWithControl && success) {
            mPointControl.fill(mAlignPoint, mCounterUsedTracks);
          }
          if (!success) {
            LOGF(error, "TracksToRecords::processROFs() - track %i h %d d %d l %d s %4d lMpos x %.2e y %.2e z %.2e gMpos x %.2e y %.2e z %.2e gRpos x %.2e y %.2e z %.2e",
                mCounterUsedTracks, mAlignPoint->half(), mAlignPoint->disk(), mAlignPoint->layer(), mAlignPoint->getSensorId(),
                mAlignPoint->getLocalMeasuredPosition().X(), mAlignPoint->getLocalMeasuredPosition().Y(), mAlignPoint->getLocalMeasuredPosition().Z(),
                mAlignPoint->getGlobalMeasuredPosition().X(), mAlignPoint->getGlobalMeasuredPosition().Y(), mAlignPoint->getGlobalMeasuredPosition().Z(),
                mAlignPoint->getGlobalRecoPosition().X(), mAlignPoint->getGlobalRecoPosition().Y(), mAlignPoint->getGlobalRecoPosition().Z());
          }

        } // end of loop on clusters

        if (isTrackUsed) {
          // copy track record
          mRecordWriter->setRecordRun(mRunNumber);
          mRecordWriter->setRecordWeight(mWeightRecord);
          const bool doPrint = false;
          mRecordWriter->fillRecordTree(doPrint); // save record data
          mCounterUsedTracks++;
        }
        nCounterAllTracks++;
      } // end of loop on tracks
    } // end of loop on ROFs

  } // end of loop on TChain reader

  LOG(info) << "TracksToRecords::processROFs() - end";
}

//__________________________________________________________________________
void TracksToRecords::printProcessTrackSummary()
{
  LOGF(info, "TracksToRecords processRecoTracks() summary: ");
  if (mNumberOfTrackChainROFs) {
    LOGF(info,
         "n ROFs = %d, used tracks = %d, skipped tracks = %d, local equations failed = %d",
         mNumberOfTrackChainROFs, mCounterUsedTracks,
         mCounterSkippedTracks, mCounterLocalEquationFailed);
  } else {
    LOGF(info,
         "n TFs = %d, used tracks = %d, skipped tracks = %d, local equations failed = %d",
         mNumberTFs, mCounterUsedTracks,
         mCounterSkippedTracks, mCounterLocalEquationFailed);
  }
}

//__________________________________________________________________________
void TracksToRecords::startRecordWriter()
{
  if (mRecordWriter) {
    mRecordWriter->init();
  }
  if (mWithControl) {
    mPointControl.setCyclicAutoSave(mNEntriesAutoSave);
    mPointControl.init();
  }
}

//__________________________________________________________________________
void TracksToRecords::endRecordWriter()
{
  if (mRecordWriter) {
    mRecordWriter->terminate(); // write record tree and close output file
  }
  if (mWithControl) {
    mPointControl.terminate();
  }
}

//__________________________________________________________________________
void TracksToRecords::startConstraintsRecWriter()
{
  if (!mWithConstraintsRecWriter) {
    return;
  }
  if (mConstraintsRecWriter) {
    mConstraintsRecWriter->changeDataBranchName();
    mConstraintsRecWriter->init();
  }
}

//__________________________________________________________________________
void TracksToRecords::endConstraintsRecWriter()
{
  if (!mWithConstraintsRecWriter) {
    return;
  }
  if (mConstraintsRecWriter) {
    mConstraintsRecWriter->terminate();
  }
}

//__________________________________________________________________________
bool TracksToRecords::setLocalDerivative(Int_t index, Double_t value)
{
  // index [0 .. 3] for {dX0, dTx, dY0, dTz}

  bool success = false;
  if (index < mNumberOfTrackParam) {
    mLocalDerivatives[index] = value;
    success = true;
  } else {
    LOGF(error,
         "AlignHelper::setLocalDerivative() - index %d >= %d",
         index, mNumberOfTrackParam);
  }
  return success;
}

//__________________________________________________________________________
bool TracksToRecords::setGlobalDerivative(Int_t index, Double_t value)
{
  // index [0 .. 3] for {dDeltaX, dDeltaY, dDeltaRz, dDeltaZ}

  bool success = false;
  if (index < mNumberOfGlobalParam) {
    mGlobalDerivatives[index] = value;
    success = true;
  } else {
    LOGF(error,
         "AlignHelper::setGlobalDerivative() - index %d >= %d",
         index, mNumberOfGlobalParam);
  }
  return success;
}

//__________________________________________________________________________
bool TracksToRecords::resetLocalDerivative()
{
  bool success = false;
  std::fill(mLocalDerivatives.begin(), mLocalDerivatives.end(), 0.);
  success = true;
  return success;
}

//__________________________________________________________________________
bool TracksToRecords::resetGlocalDerivative()
{
  bool success = false;
  std::fill(mGlobalDerivatives.begin(), mGlobalDerivatives.end(), 0.);
  success = true;
  return success;
}

//__________________________________________________________________________
bool TracksToRecords::setLocalEquationX()
{

  if (!mAlignPoint->isAlignPointSet()) {
    LOGF(error,
         "TracksToRecords::setLocalEquationX() - no align point coordinates set !");
    return false;
  }
  if (!mAlignPoint->isGlobalDerivativeDone()) {
    return false;
  }
  if (!mAlignPoint->isLocalDerivativeDone()) {
    return false;
  }

  bool success = true;

  // clean slate for the local equation for this measurement

  success &= resetGlocalDerivative();
  success &= resetLocalDerivative();

  // local derivatives
  // index [0 .. 3] for {dX0, dTx, dY0, dTz}

  success &= setLocalDerivative(0, 0.0 /*mAlignPoint->localDerivativeX().dX0()*/);
  success &= setLocalDerivative(1, mAlignPoint->localDerivativeX().dTx());
  success &= setLocalDerivative(2, 0.0 /*mAlignPoint->localDerivativeX().dY0()*/);
  success &= setLocalDerivative(3, mAlignPoint->localDerivativeX().dTy());

  // global derivatives
  // index [0 .. 3] for {dDeltaX, dDeltaY, dDeltaRz, dDeltaZ}

  Int_t chipId = mAlignPoint->getSensorId();
  success &= setGlobalDerivative(chipId * mNDofPerSensor + 0, mAlignPoint->globalDerivativeX().dDeltaX());
  success &= setGlobalDerivative(chipId * mNDofPerSensor + 1, mAlignPoint->globalDerivativeX().dDeltaY());
  success &= setGlobalDerivative(chipId * mNDofPerSensor + 2, mAlignPoint->globalDerivativeX().dDeltaRz());
  success &= setGlobalDerivative(chipId * mNDofPerSensor + 3, mAlignPoint->globalDerivativeX().dDeltaZ());

  if (success) {
    if (mCounterUsedTracks < 5) {
      LOGF(debug,
           "TracksToRecords::setLocalEquationX(): track %i sr %4d local %.3e %.3e %.3e %.3e, global %.3e %.3e %.3e %.3e X %.3e sigma %.3e",
           mCounterUsedTracks, chipId,
           mLocalDerivatives[0], mLocalDerivatives[1], mLocalDerivatives[2], mLocalDerivatives[3],
           mGlobalDerivatives[chipId * mNDofPerSensor + 0],
           mGlobalDerivatives[chipId * mNDofPerSensor + 1],
           mGlobalDerivatives[chipId * mNDofPerSensor + 2],
           mGlobalDerivatives[chipId * mNDofPerSensor + 3],
           mAlignPoint->getLocalResidual().X(),
           mAlignPoint->getLocalMeasuredPositionSigma().X());
    }
    mMillepede->SetLocalEquation(
      mGlobalDerivatives,
      mLocalDerivatives,
      mAlignPoint->getLocalResidual().X(),
      mAlignPoint->getLocalMeasuredPositionSigma().X());
  } else {
    mCounterLocalEquationFailed++;
  }

  return success;
}

//__________________________________________________________________________
bool TracksToRecords::setLocalEquationY()
{
  if (!mAlignPoint->isAlignPointSet()) {
    LOGF(error,
         "TracksToRecords::setLocalEquationY() - no align point coordinates set !");
    return false;
  }
  if (!mAlignPoint->isGlobalDerivativeDone()) {
    return false;
  }
  if (!mAlignPoint->isLocalDerivativeDone()) {
    return false;
  }

  bool success = true;

  // clean slate for the local equation for this measurement

  success &= resetGlocalDerivative();
  success &= resetLocalDerivative();

  // local derivatives
  // index [0 .. 3] for {dX0, dTx, dY0, dTz}

  success &= setLocalDerivative(0, 0.0 /*mAlignPoint->localDerivativeY().dX0()*/);
  success &= setLocalDerivative(1, mAlignPoint->localDerivativeY().dTx());
  success &= setLocalDerivative(2, 0.0 /*mAlignPoint->localDerivativeY().dY0()*/);
  success &= setLocalDerivative(3, mAlignPoint->localDerivativeY().dTy());

  // global derivatives
  // index [0 .. 3] for {dDeltaX, dDeltaY, dDeltaRz, dDeltaZ}

  Int_t chipId = mAlignPoint->getSensorId();
  success &= setGlobalDerivative(chipId * mNDofPerSensor + 0, mAlignPoint->globalDerivativeY().dDeltaX());
  success &= setGlobalDerivative(chipId * mNDofPerSensor + 1, mAlignPoint->globalDerivativeY().dDeltaY());
  success &= setGlobalDerivative(chipId * mNDofPerSensor + 2, mAlignPoint->globalDerivativeY().dDeltaRz());
  success &= setGlobalDerivative(chipId * mNDofPerSensor + 3, mAlignPoint->globalDerivativeY().dDeltaZ());

  if (success) {
    if (mCounterUsedTracks < 5) {
      LOGF(debug,
           "TracksToRecords::setLocalEquationY(): track %i sr %4d local %.3e %.3e %.3e %.3e, global %.3e %.3e %.3e %.3e Y %.3e sigma %.3e",
           mCounterUsedTracks, chipId,
           mLocalDerivatives[0], mLocalDerivatives[1], mLocalDerivatives[2], mLocalDerivatives[3],
           mGlobalDerivatives[chipId * mNDofPerSensor + 0],
           mGlobalDerivatives[chipId * mNDofPerSensor + 1],
           mGlobalDerivatives[chipId * mNDofPerSensor + 2],
           mGlobalDerivatives[chipId * mNDofPerSensor + 3],
           mAlignPoint->getLocalResidual().Y(),
           mAlignPoint->getLocalMeasuredPositionSigma().Y());
    }
    mMillepede->SetLocalEquation(
      mGlobalDerivatives,
      mLocalDerivatives,
      mAlignPoint->getLocalResidual().Y(),
      mAlignPoint->getLocalMeasuredPositionSigma().Y());
  } else {
    mCounterLocalEquationFailed++;
  }

  return success;
}

//__________________________________________________________________________
bool TracksToRecords::setLocalEquationZ()
{
  if (!mAlignPoint->isAlignPointSet()) {
    LOGF(error,
         "TracksToRecords::setLocalEquationZ() - no align point coordinates set !");
    return false;
  }
  if (!mAlignPoint->isGlobalDerivativeDone()) {
    return false;
  }
  if (!mAlignPoint->isLocalDerivativeDone()) {
    return false;
  }

  bool success = true;

  // clean slate for the local equation for this measurement

  success &= resetGlocalDerivative();
  success &= resetLocalDerivative();

  // local derivatives
  // index [0 .. 3] for {dX0, dTx, dY0, dTz}

  success &= setLocalDerivative(0, 0.0 /*mAlignPoint->localDerivativeZ().dX0()*/);
  success &= setLocalDerivative(1, mAlignPoint->localDerivativeZ().dTx());
  success &= setLocalDerivative(2, 0.0 /*mAlignPoint->localDerivativeZ().dY0()*/);
  success &= setLocalDerivative(3, mAlignPoint->localDerivativeZ().dTy());

  // global derivatives
  // index [0 .. 3] for {dDeltaX, dDeltaY, dDeltaRz, dDeltaZ}

  Int_t chipId = mAlignPoint->getSensorId();
  success &= setGlobalDerivative(chipId * mNDofPerSensor + 0, mAlignPoint->globalDerivativeZ().dDeltaX());
  success &= setGlobalDerivative(chipId * mNDofPerSensor + 1, mAlignPoint->globalDerivativeZ().dDeltaY());
  success &= setGlobalDerivative(chipId * mNDofPerSensor + 2, mAlignPoint->globalDerivativeZ().dDeltaRz());
  success &= setGlobalDerivative(chipId * mNDofPerSensor + 3, mAlignPoint->globalDerivativeZ().dDeltaZ());

  if (success) {
    if (mCounterUsedTracks < 5) {
      LOGF(debug,
           "TracksToRecords::setLocalEquationZ(): track %i sr %4d local %.3e %.3e %.3e %.3e, global %.3e %.3e %.3e %.3e Z %.3e sigma %.3e",
           mCounterUsedTracks, chipId,
           mLocalDerivatives[0], mLocalDerivatives[1], mLocalDerivatives[2], mLocalDerivatives[3],
           mGlobalDerivatives[chipId * mNDofPerSensor + 0],
           mGlobalDerivatives[chipId * mNDofPerSensor + 1],
           mGlobalDerivatives[chipId * mNDofPerSensor + 2],
           mGlobalDerivatives[chipId * mNDofPerSensor + 3],
           mAlignPoint->getLocalResidual().Z(),
           mAlignPoint->getLocalMeasuredPositionSigma().Z());
    }
    mMillepede->SetLocalEquation(
      mGlobalDerivatives,
      mLocalDerivatives,
      mAlignPoint->getLocalResidual().Z(),
      mAlignPoint->getLocalMeasuredPositionSigma().Z());
  } else {
    mCounterLocalEquationFailed++;
  }

  return success;
}
