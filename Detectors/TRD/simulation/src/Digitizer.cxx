// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include <TGeoManager.h>
#include <TRandom.h>

#include "FairLogger.h"
#include "DetectorsBase/GeometryManager.h"

#include "TRDBase/TRDGeometry.h"
#include "TRDBase/TRDSimParam.h"
#include "TRDBase/TRDPadPlane.h"
#include "TRDBase/PadResponse.h"

#include "TRDSimulation/Digitizer.h"
#include "TRDSimulation/TRDSimParams.h"
#include <cmath>

#ifdef WITH_OPENMP
#include <omp.h>
#endif

using namespace o2::trd;
using namespace o2::math_utils;

// init method for late initialization
void Digitizer::init()
{
  mGeo = TRDGeometry::instance();
  mGeo->createClusterMatrixArray();          // Requiered for chamberInGeometry()
  mPRF = new PadResponse();                  // Pad response function initialization
  mSimParam = TRDSimParam::Instance();       // Instance for simulation parameters
  mCommonParam = TRDCommonParam::Instance(); // Instance for common parameters
  if (!mSimParam) {
  }
  if (!mCommonParam) {
  } else {
    if (!mCommonParam->cacheMagField()) {
    }
  }

  // obtain the number of threads from configuration
#ifdef WITH_OPENMP
  int askedthreads = TRDSimParams::Instance().digithreads;
  int maxthreads = omp_get_max_threads();
  if (askedthreads < 0) {
    mNumThreads = maxthreads;
  } else {
    mNumThreads = std::min(maxthreads, askedthreads);
  }
  LOG(INFO) << "TRD: Digitizing with " << mNumThreads << " threads ";
#endif

  // initialize structures that we need per thread
  mGausRandomRings.resize(mNumThreads);
  mFlatRandomRings.resize(mNumThreads);
  mLogRandomRings.resize(mNumThreads);
  for (int i = 0; i < mNumThreads; ++i) {
    mGausRandomRings[i].initialize(RandomRing<>::RandomType::Gaus);
    mFlatRandomRings[i].initialize(RandomRing<>::RandomType::Flat);
    mLogRandomRings[i].initialize([]() -> float { return std::log(gRandom->Rndm()); });
    mDriftEstimators.emplace_back();
  }

  setSimulationParameters();
}

void Digitizer::setSimulationParameters()
{
  mNpad = mSimParam->getNumberOfPadsInPadResponse(); // Number of pads included in the pad response
  if (mSimParam->TRFOn()) {
    mTimeBinTRFend = ((int)(mSimParam->GetTRFhi() * mCommonParam->GetSamplingFrequency())) - 1;
  }
  mMaxTimeBins = kTimeBins;     // for signals, usually set at 30 tb = 3 microseconds
  mMaxTimeBinsTRAP = kTimeBins; // for adcs; should be read from the CCDB or the TRAP config
  mSamplingRate = mCommonParam->GetSamplingFrequency();
  mElAttachProp = mSimParam->GetElAttachProp() / 100;
}

void Digitizer::clearCollections()
{
  for (int det = 0; det < kNdet; ++det) {
    mSignalsMapCollection[det].clear();
    mDigitsCollection[det].clear();
  }
  mMergedLabels.clear();
}

void Digitizer::flush(DigitContainer& digits, o2::dataformats::MCTruthContainer<MCLabel>& labels)
{
  // Convert Signals to Digits over all TRD detectors (in a parallel fashion)
#ifdef WITH_OPENMP
  omp_set_num_threads(mNumThreads);
#pragma omp parallel for schedule(dynamic)
#endif
  for (int det = 0; det < kNdet; ++det) {
#ifdef WITH_OPENMP
    const int threadid = omp_get_thread_num();
#else
    const int threadid = 0;
#endif
    if (mSignalsMapCollection[det].size() == 0) {
      continue;
    }
    bool status = convertSignalsToADC(det, mSignalsMapCollection[det], mDigitsCollection[det], threadid);
    if (!status) {
      LOG(WARN) << "TRD conversion of signals to digits failed for detector " << det;
      continue; // go to the next chamber
    }
  }

  // Finalize
  for (int det = 0; det < kNdet; ++det) {
    std::move(mDigitsCollection[det].begin(), mDigitsCollection[det].end(), std::back_inserter(digits));
    for (const auto& it : mSignalsMapCollection[det]) {
      std::move(it.second.labels.begin(), it.second.labels.end(), std::back_inserter(mMergedLabels));
    }
  }
  labels.addElements(labels.getIndexedSize(), mMergedLabels);
  clearCollections();
}

void Digitizer::process(std::vector<HitType> const& hits, DigitContainer& digits, o2::dataformats::MCTruthContainer<MCLabel>& labels)
{
  if (!mCalib) {
    LOG(FATAL) << "TRD Calibration database not available";
  }

  // if there are digits for longer than 3 microseconds, flush them
  // todo: check for pileup
  // constexpr double totalSamplingTime = kTimeBins / mCommonParam->GetSamplingFrequency();
  if ((mTime - mLastTime) > 0) {
    flush(digits, labels);
    if (digits.size() == 0) {
      LOG(WARN) << "Digit container is empty after flushing signals";
    }
  }
  mLastTime = mTime; // always update the last time to the curret time for the next iteration

  // Get the a hit container for all the hits in a given detector then call convertHits for a given detector (0 - 539)
  std::array<std::vector<HitType>, kNdet> hitsPerDetector;
  getHitContainerPerDetector(hits, hitsPerDetector);

#ifdef WITH_OPENMP
  omp_set_num_threads(mNumThreads);
// Loop over all TRD detectors (in a parallel fashion)
#pragma omp parallel for schedule(dynamic)
#endif
  for (int det = 0; det < kNdet; ++det) {
#ifdef WITH_OPENMP
    const int threadid = omp_get_thread_num();
#else
    const int threadid = 0;
#endif
    auto& signalsMap = mSignalsMapCollection[det];
    // Jump to the next detector if the detector is
    // switched off, not installed, etc
    if (mCalib->isChamberNoData(det)) {
      continue;
    }
    if (!mGeo->chamberInGeometry(det)) {
      continue;
    }

    // Go to the next detector if there are no hits
    if (hitsPerDetector[det].size() == 0) {
      continue;
    }

    if (!convertHits(det, hitsPerDetector[det], signalsMap, threadid)) {
      LOG(WARN) << "TRD conversion of hits failed for detector " << det;
      continue; // go to the next chamber
    }
  }
}

void Digitizer::getHitContainerPerDetector(const std::vector<HitType>& hits, std::array<std::vector<HitType>, kNdet>& hitsPerDetector)
{
  //
  // Fill an array of size kNdet (540)
  // The i-element of the array contains the hit collection for the i-detector
  // To be called once, before doing the loop over all detectors and process the hits
  //
  for (const auto& hit : hits) {
    hitsPerDetector[hit.GetDetectorID()].push_back(hit);
  }
}

bool Digitizer::convertHits(const int det, const std::vector<HitType>& hits, SignalContainer& signalMapCont, int thread)
{
  //
  // Convert the detector-wise sorted hits to detector signals
  //

  double padSignal[mNpad];

  const double calExBDetValue = mCalib->getExB(det); // T * V/cm (check units)
  const TRDPadPlane* padPlane = mGeo->getPadPlane(det);
  const int layer = mGeo->getLayer(det);
  const float rowEndROC = padPlane->getRowEndROC();
  const float row0 = padPlane->getRow0ROC();
  const int nRowMax = padPlane->getNrows();
  const int nColMax = padPlane->getNcols();

  // Loop over hits
  for (const auto& hit : hits) {
    const int qTotal = hit.GetCharge();
    /*
      Now the real local coordinate system of the ROC
      column direction: locC
      row direction:    locR
      time direction:   locT
      locR and locC are identical to the coordinates of the corresponding
      volumina of the drift or amplification region.
      locT is defined relative to the wire plane (i.e. middle of amplification
      region), meaning locT = 0, and is negative for hits coming from the
      drift region.
    */
    double locC = hit.getLocalC();        // col direction in amplification or drift volume
    double locR = hit.getLocalR();        // row direction in amplification or drift volume
    double locT = hit.getLocalT();        // time direction in amplification or drift volume
    const double driftLength = -1 * locT; // The drift length in cm without diffusion
    // Patch to take care of TR photons that are absorbed
    // outside the chamber volume. A real fix would actually need
    // a more clever implementation of the TR hit generation
    if (qTotal < 0) {
      if ((locR < rowEndROC) || (locR > row0)) {
        continue;
      }
      if ((driftLength < DrMin) || (driftLength > DrMax)) {
        continue;
      }
    }

    int rowE = padPlane->getPadRowNumberROC(locR);
    if (rowE < 0) {
      continue;
    }

    double rowOffset = padPlane->getPadRowOffsetROC(rowE, locR);
    double offsetTilt = padPlane->getTiltOffset(rowOffset);
    int colE = padPlane->getPadColNumber(locC + offsetTilt);
    if (colE < 0) {
      continue;
    }

    double absDriftLength = std::fabs(driftLength); // Normalized drift length
    if (mCommonParam->ExBOn()) {
      absDriftLength /= std::sqrt(1 / (1 + calExBDetValue * calExBDetValue));
    }

    float driftVelocity = mCalib->getVDrift(det, colE, rowE); // The drift velocity
    float t0 = mCalib->getT0(det, colE, rowE);                // The T0 velocity

    // Loop over all created electrons
    const int nElectrons = std::fabs(qTotal);
    for (int el = 0; el < nElectrons; ++el) {
      // Electron attachment
      if (mSimParam->ElAttachOn()) {
        if (mFlatRandomRings[thread].getNextValue() < absDriftLength * mElAttachProp) {
          continue;
        }
      }
      // scoped diffused coordinates for each electron
      double locRd{locR}, locCd{locC}, locTd{locT};

      // Apply diffusion smearing
      if (mSimParam->DiffusionOn()) {
        if (!diffusion(driftVelocity, absDriftLength, calExBDetValue, locR, locC, locT, locRd, locCd, locTd, thread)) {
          continue;
        }
      }

      // Apply E x B effects
      if (mCommonParam->ExBOn()) {
        locCd = locCd + calExBDetValue * driftLength;
      }
      // The electron position after diffusion and ExB in pad coordinates.
      rowE = padPlane->getPadRowNumberROC(locRd);
      if (rowE < 0) {
        continue;
      }
      rowOffset = padPlane->getPadRowOffsetROC(rowE, locRd);
      // The pad column (rphi-direction)
      offsetTilt = padPlane->getTiltOffset(rowOffset);
      colE = padPlane->getPadColNumber(locCd + offsetTilt);
      if (colE < 0) {
        continue;
      }
      const double colOffset = padPlane->getPadColOffset(colE, locCd + offsetTilt);
      driftVelocity = mCalib->getVDrift(det, colE, rowE); // The drift velocity for the updated col and row
      t0 = mCalib->getT0(det, colE, rowE);                // The T0 velocity for the updated col and row
      // Convert the position to drift time [mus], using either constant drift velocity or
      // time structure of drift cells (non-isochronity, GARFIELD calculation).
      // Also add absolute time of hits to take pile-up events into account properly
      double driftTime;
      if (mSimParam->TimeStructOn()) {
        // Get z-position with respect to anode wire
        double zz = row0 - locRd + padPlane->getAnodeWireOffset();
        zz -= ((int)(2 * zz)) * 0.5;
        if (zz > 0.25) {
          zz = 0.5 - zz;
        }
        // Use drift time map (GARFIELD)
        driftTime = mDriftEstimators[thread].TimeStruct(driftVelocity, 0.5 * AmWidth - 1.0 * locTd, zz) + hit.GetTime();
      } else {
        // Use constant drift velocity
        driftTime = std::fabs(locTd) / driftVelocity + hit.GetTime(); // drift time in microseconds
      }

      // Apply the gas gain including fluctuations
      const double signal = -(mSimParam->GetGasGain()) * mLogRandomRings[thread].getNextValue();

      // Apply the pad response
      if (mSimParam->PRFOn()) {
        // The distance of the electron to the center of the pad in units of pad width
        double dist = (colOffset - 0.5 * padPlane->getColSize(colE)) / padPlane->getColSize(colE);
        // ********************************************************************************
        // This is a fixed parametrization, i.e. not dependent on calibration values !
        // ********************************************************************************
        if (!(mPRF->getPRF(signal, dist, layer, padSignal))) {
          continue;
        }
      } else {
        padSignal[0] = 0;
        padSignal[1] = signal;
        padSignal[2] = 0;
      }

      // The time bin (always positive), with t0 distortion
      double timeBinIdeal = driftTime * mSamplingRate + t0;
      // Protection
      if (std::fabs(timeBinIdeal) > (2 * mMaxTimeBins)) {
        timeBinIdeal = 2 * mMaxTimeBins;
      }
      int timeBinTruncated = ((int)timeBinIdeal);
      // The distance of the position to the middle of the timebin
      double timeOffset = ((float)timeBinTruncated + 0.5 - timeBinIdeal) / mSamplingRate;
      // Sample the time response inside the drift region + additional time bins before and after.
      // The sampling is done always in the middle of the time bin
      const int firstTimeBin = std::max(timeBinTruncated, 0);
      const int lastTimeBin = std::min(timeBinTruncated + mTimeBinTRFend, mMaxTimeBins);

      // loop over pads first then over timebins for better cache friendliness
      // and less access to signalMapCont
      for (int pad = 0; pad < mNpad; ++pad) {
        int colPos = colE + pad - 1;
        if (colPos < 0) {
          continue;
        }
        if (colPos >= nColMax) {
          break;
        }

        const int key = calculateKey(det, rowE, colPos);
        auto& currentSignalData = signalMapCont[key]; // Get the old signal or make a new one if it doesn't exist
        auto& currentSignal = currentSignalData.signals;
        auto& trackIds = currentSignalData.trackIds;
        auto& labels = currentSignalData.labels;

        // add a label record
        addLabel(hit, labels, trackIds);

        // add signal with crosstalk for the non-central pads only
        if (colPos != colE) {
          for (int tb = firstTimeBin; tb < lastTimeBin; ++tb) {
            const double t = (tb - timeBinTruncated) / mSamplingRate + timeOffset;
            const double timeResponse = mSimParam->TRFOn() ? mSimParam->TimeResponse(t) : 1;
            const double crossTalk = mSimParam->CTOn() ? mSimParam->CrossTalk(t) : 0;
            currentSignal[tb] += padSignal[pad] * (timeResponse + crossTalk);
          } // end of loop time bins
        } else {
          for (int tb = firstTimeBin; tb < lastTimeBin; ++tb) {
            const double t = (tb - timeBinTruncated) / mSamplingRate + timeOffset;
            const double timeResponse = mSimParam->TRFOn() ? mSimParam->TimeResponse(t) : 1;
            currentSignal[tb] += padSignal[pad] * timeResponse;
          } // end of loop time bins
        }
      } // end of loop over pads
    }   // end of loop over electrons
  }     // end of loop over hits
  return true;
}

void Digitizer::addLabel(const o2::trd::HitType& hit, std::vector<o2::trd::MCLabel>& labels, std::unordered_map<int, int>& trackIds)
{
  if (trackIds[hit.GetTrackID()] == 0) {
    trackIds[hit.GetTrackID()] = 1;
    MCLabel label(hit.GetTrackID(), getEventID(), getSrcID());
    labels.push_back(label);
  }
}

float drawGaus(o2::math_utils::RandomRing<>& normaldistRing, float mu, float sigma)
{
  // this is using standard normally distributed random numbers and rescaling to make
  // them gaussian distributed with general mu and sigma
  return mu + sigma * normaldistRing.getNextValue();
}

bool Digitizer::convertSignalsToADC(const int det, SignalContainer& signalMapCont, DigitContainer& digits, int thread)
{
  //
  // Converts the sampled electron signals to ADC values for a given chamber
  //
  if (signalMapCont.size() == 0) {
    return false;
  }

  constexpr double kEl2fC = 1.602e-19 * 1.0e15;                                 // Converts number of electrons to fC
  double coupling = mSimParam->GetPadCoupling() * mSimParam->GetTimeCoupling(); // Coupling factor
  double convert = kEl2fC * mSimParam->GetChipGain();                           // Electronics conversion factor
  double adcConvert = mSimParam->GetADCoutRange() / mSimParam->GetADCinRange(); // ADC conversion factor
  double baseline = mSimParam->GetADCbaseline() / adcConvert;                   // The electronics baseline in mV
  double baselineEl = baseline / convert;                                       // The electronics baseline in electrons

  for (auto& signalMapIter : signalMapCont) {
    const auto key = signalMapIter.first;
    const int det = getDetectorFromKey(key);
    const int row = getRowFromKey(key);
    const int col = getColFromKey(key);
    // halfchamber masking
    int mcm = (int)(col / 18);               // current group of 18 col pads
    int halfchamberside = (mcm > 3 ? 1 : 0); // 0=Aside, 1=Bside

    // Halfchambers that are switched off, masked by mCalib
    /* Something is wrong with isHalfChamberNoData - deactivated for now
    if (mCalib->isHalfChamberNoData(det, halfchamberside)) {
      continue;
    }
    */

    // Check whether pad is masked
    // Bridged pads are not considered yet!!!
    if (mCalib->isPadMasked(det, col, row) || mCalib->isPadNotConnected(det, col, row)) {
      continue;
    }

    float padgain = mCalib->getPadGainFactor(det, row, col); // The gain factor
    if (padgain <= 0) {
      LOG(FATAL) << "Not a valid gain " << padgain << ", " << det << ", " << col << ", " << row;
    }

    // Loop over the all timebins in the ADC array
    SignalArray& signalData = signalMapIter.second;
    auto& signalArray = signalData.signals;
    ArrayADC adcs{};
    for (int tb = 0; tb < mMaxTimeBinsTRAP; ++tb) {
      float signalAmp = (float)signalArray[tb]; // The signal amplitude
      signalAmp *= coupling;                    // Pad and time coupling
      signalAmp *= padgain;                     // Gain factors
      // Add the noise, starting from minus ADC baseline in electrons
      signalAmp = std::max((double)drawGaus(mGausRandomRings[thread], signalAmp, mSimParam->GetNoise()), -baselineEl);
      signalAmp *= convert;  // Convert to mV
      signalAmp += baseline; // Add ADC baseline in mV
      // Convert to ADC counts
      // Set the overflow-bit fADCoutRange if the signal is larger than fADCinRange
      ADC_t adc = 0;
      if (signalAmp >= mSimParam->GetADCinRange()) {
        adc = ((ADC_t)mSimParam->GetADCoutRange());
      } else {
        adc = std::lround(signalAmp * adcConvert);
      }
      // update the adc array value
      adcs[tb] = adc;
    } // loop over timebins
    // Convert the map to digits here, and push them to the container
    digits.emplace_back(det, row, col, adcs, getEventTime());
  } // loop over digits
  return true;
}

bool Digitizer::diffusion(float vdrift, float absdriftlength, float exbvalue,
                          float lRow0, float lCol0, float lTime0,
                          double& lRow, double& lCol, double& lTime, int thread)
{
  //
  // Applies the diffusion smearing to the position of a single electron.
  // Depends on absolute drift length.
  //
  float diffL = 0.0;
  float diffT = 0.0;
  if (mDriftEstimators[thread].GetDiffCoeff(diffL, diffT, vdrift)) {
    float driftSqrt = std::sqrt(absdriftlength);
    float sigmaT = driftSqrt * diffT;
    float sigmaL = driftSqrt * diffL;
    lRow = drawGaus(mGausRandomRings[thread], lRow0, sigmaT);
    if (mCommonParam->ExBOn()) {
      const float exbfactor = 1.f / (1.f + exbvalue * exbvalue);
      lCol = drawGaus(mGausRandomRings[thread], lCol0, sigmaT * exbfactor);
      lTime = drawGaus(mGausRandomRings[thread], lTime0, sigmaL * exbfactor);
    } else {
      lCol = drawGaus(mGausRandomRings[thread], lCol0, sigmaT);
      lTime = drawGaus(mGausRandomRings[thread], lTime0, sigmaL);
    }
    return true;
  } else {
    return false;
  }
}
