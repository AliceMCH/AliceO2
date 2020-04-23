// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file GPUTPCGMMerger.h
/// \author Sergey Gorbunov, David Rohr

#ifndef GPUTPCGMMERGER_H
#define GPUTPCGMMERGER_H

#include "GPUParam.h"
#include "GPUTPCDef.h"
#include "GPUTPCGMBorderTrack.h"
#include "GPUTPCGMMergedTrack.h"
#include "GPUTPCGMSliceTrack.h"
#include "GPUCommonDef.h"
#include "GPUProcessor.h"
#include "GPUTPCGMMergerTypes.h"

#if !defined(GPUCA_GPUCODE)
#include <cmath>
#include <iostream>
#endif // GPUCA_GPUCODE

namespace o2
{
namespace base
{
class MatLayerCylSet;
}
} // namespace o2

namespace GPUCA_NAMESPACE
{
namespace gpu
{
class GPUTPCSliceTrack;
class GPUTPCSliceOutput;
class GPUTPCGMTrackParam;
class GPUTPCTracker;
class GPUChainTracking;
class GPUTPCGMPolynomialField;
struct GPUTPCGMLoopData;

/**
 * @class GPUTPCGMMerger
 *
 */
class GPUTPCGMMerger : public GPUProcessor
{
 public:
  GPUTPCGMMerger();
  ~GPUTPCGMMerger() CON_DEFAULT;
  GPUTPCGMMerger(const GPUTPCGMMerger&) CON_DELETE;
  const GPUTPCGMMerger& operator=(const GPUTPCGMMerger&) const CON_DELETE;
  static CONSTEXPR int NSLICES = GPUCA_NSLICES; //* N slices

  struct memory {
    GPUAtomic(unsigned int) nRetryRefit;
    GPUAtomic(unsigned int) nLoopData;
    const GPUTPCTrack* firstGlobalTracks[NSLICES];
  };

  struct trackCluster {
    unsigned int id;
    unsigned char row;
    unsigned char slice;
    unsigned char leg;
  };

  void InitializeProcessor();
  void RegisterMemoryAllocation();
  void SetMaxData(const GPUTrackingInOutPointers& io);
  void* SetPointersHostOnly(void* mem);
  void* SetPointersGPURefit(void* mem);

  void OverrideSliceTracker(GPUTPCTracker* trk) { mSliceTrackers = trk; }
  void SetTrackingChain(GPUChainTracking* c) { mChainTracking = c; }
  const GPUChainTracking* GetTrackingChain() const { return mChainTracking; }

  void SetSliceData(int index, const GPUTPCSliceOutput* SliceData);

  GPUhd() int NOutputTracks() const { return mNOutputTracks; }
  GPUhd() const GPUTPCGMMergedTrack* OutputTracks() const { return mOutputTracks; }
  GPUhd() GPUTPCGMMergedTrack* OutputTracks()
  {
    return mOutputTracks;
  }

  GPUhd() void SetMatLUT(const o2::base::MatLayerCylSet* lut) { mMatLUT = lut; }
  GPUhd() const o2::base::MatLayerCylSet* MatLUT() const { return mMatLUT; }

  GPUhd() unsigned int NClusters() const { return (mNClusters); }
  GPUhd() unsigned int NMaxClusters() const { return (mNMaxClusters); }
  GPUhd() unsigned int NMaxTracks() const { return (mNMaxTracks); }
  GPUhd() unsigned int NMaxOutputTrackClusters() const { return (mNMaxOutputTrackClusters); }
  GPUhd() unsigned int NOutputTrackClusters() const { return (mNOutputTrackClusters); }
  GPUhd() const GPUTPCGMMergedTrackHit* Clusters() const { return (mClusters); }
  GPUhd() GPUTPCGMMergedTrackHit* Clusters()
  {
    return (mClusters);
  }
  GPUhd() const GPUTPCTracker* SliceTrackers() const { return (mSliceTrackers); }
  GPUhd() GPUAtomic(unsigned int) * ClusterAttachment() const { return (mClusterAttachment); }
  GPUhd() unsigned int* TrackOrderAttach() const { return mTrackOrderAttach; }
  GPUhd() unsigned int* TrackOrderProcess() const { return mTrackOrderProcess; }
  GPUd() unsigned int NSlowTracks() const { return mNSlowTracks; }
  GPUd() unsigned int* RetryRefitIds() const { return mRetryRefitIds; }
  GPUd() GPUTPCGMLoopData* LoopData() const { return mLoopData; }
  GPUd() memory* Memory() const { return mMemory; }

  short MemoryResRefit() { return mMemoryResRefit; }

  int RefitSliceTrack(GPUTPCGMSliceTrack& sliceTrack, const GPUTPCTrack* inTrack, float alpha, int slice);
  void SetTrackClusterZT(GPUTPCGMSliceTrack& track, int iSlice, const GPUTPCTrack* sliceTr);

  void UnpackSlices();
  void MergeCEInit();
  void MergeCE();
  void MergeWithingSlices();
  void MergeSlices();
  void PrepareClustersForFit();
  void CollectMergedTracks();
  void Finalize();

#ifndef GPUCA_GPUCODE
  void DumpSliceTracks(std::ostream& out);
  void DumpMergedWithinSlices(std::ostream& out);
  void DumpMergedBetweenSlices(std::ostream& out);
  void DumpCollected(std::ostream& out);
  void DumpMergeCE(std::ostream& out);
  void DumpFitPrepare(std::ostream& out);
  void DumpRefit(std::ostream& out);
  void DumpFinal(std::ostream& out);
#endif

 private:
  void MakeBorderTracks(int iSlice, int iBorder, GPUTPCGMBorderTrack B[], int& nB, bool useOrigTrackParam = false);
  void MergeBorderTracks(int iSlice1, GPUTPCGMBorderTrack B1[], int N1, int iSlice2, GPUTPCGMBorderTrack B2[], int N2, int mergeMode = 0);

  void MergeCEFill(const GPUTPCGMSliceTrack* track, const GPUTPCGMMergedTrackHit& cls, int itr);
  void ResolveMergeSlices(bool useOrigTrackParam, bool mergeAll);
  void MergeSlicesStep(int border0, int border1, bool useOrigTrackParam);
  void ClearTrackLinks(int n);

  void CheckMergedTracks();
#ifndef GPUCA_GPUCODE
  void PrintMergeGraph(const GPUTPCGMSliceTrack* trk, std::ostream& out);
  int GetTrackLabel(const GPUTPCGMBorderTrack& trk);
#endif

  int SliceTrackInfoFirst(int iSlice) { return mSliceTrackInfoIndex[iSlice]; }
  int SliceTrackInfoLast(int iSlice) { return mSliceTrackInfoIndex[iSlice + 1]; }
  int SliceTrackInfoGlobalFirst(int iSlice) { return mSliceTrackInfoIndex[NSLICES + iSlice]; }
  int SliceTrackInfoGlobalLast(int iSlice) { return mSliceTrackInfoIndex[NSLICES + iSlice + 1]; }
  int SliceTrackInfoLocalTotal() { return mSliceTrackInfoIndex[NSLICES]; }
  int SliceTrackInfoTotal() { return mSliceTrackInfoIndex[2 * NSLICES]; }

  int mNextSliceInd[NSLICES];
  int mPrevSliceInd[NSLICES];

  const GPUTPCSliceOutput* mkSlices[NSLICES]; //* array of input slice tracks

  int* mTrackLinks;

  unsigned int mNMaxSliceTracks;         // maximum number of incoming slice tracks
  unsigned int mNMaxTracks;              // maximum number of output tracks
  unsigned int mNMaxSingleSliceTracks;   // max N tracks in one slice
  unsigned int mNMaxOutputTrackClusters; // max number of clusters in output tracks (double-counting shared clusters)
  unsigned int mNMaxClusters;            // max total unique clusters (in event)

  short mMemoryResRefit;

  int mMaxID;
  int mNClusters; // Total number of incoming clusters (from slice tracks)
  int mNOutputTracks;
  int mNOutputTrackClusters;
  GPUTPCGMMergedTrack* mOutputTracks; //* array of output merged tracks

  GPUTPCGMSliceTrack* mSliceTrackInfos; //* additional information for slice tracks
  int mSliceTrackInfoIndex[NSLICES * 2 + 1];
  GPUTPCGMMergedTrackHit* mClusters;
  int* mGlobalClusterIDs;
  GPUAtomic(unsigned int) * mClusterAttachment;
  unsigned int* mTrackOrderAttach;
  unsigned int* mTrackOrderProcess;
  unsigned int mNSlowTracks;
  char* mTmpMem;
  GPUTPCGMBorderTrack* mBorderMemory; // memory for border tracks
  GPUTPCGMBorderTrack* mBorder[NSLICES];
  GPUTPCGMBorderTrack::Range* mBorderRangeMemory;    // memory for border tracks
  GPUTPCGMBorderTrack::Range* mBorderRange[NSLICES]; // memory for border tracks
  int mBorderCETracks[2][NSLICES];
  memory* mMemory;
  unsigned int* mRetryRefitIds;
  GPUTPCGMLoopData* mLoopData;

  const GPUTPCTracker* mSliceTrackers;
  const o2::base::MatLayerCylSet* mMatLUT;
  GPUChainTracking* mChainTracking; // Tracking chain with access to input data / parameters
};
} // namespace gpu
} // namespace GPUCA_NAMESPACE

#endif // GPUTPCGMMERGER_H
