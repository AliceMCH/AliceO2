// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include "MCHCalibration/PedestalProcessor.h"
#include <cmath>
#include <iostream>

namespace o2::mch::calibration
{

  PedestalProcessor::PedestalProcessor()
  {
    reset();
  }

  void PedestalProcessor::reset()
  {
    for (int s = 0; s < MCH_NUMBER_OF_SOLAR; s++) {
      for (int i = 0; i < 40; i++) {
        for (int j = 0; j < 64; j++) {
          mNhits[s][i][j] = 0;
          mPedestal[s][i][j] = 0;
          mNoise[s][i][j] = 0;
        }
      }
    }
  }

  void PedestalProcessor::processDigits(gsl::span<const PedestalDigit> digits)
  {
    bool mDebug = false;
    for (auto& d : digits) {
      auto solarId = d.getSolarId();
      auto dsId = d.getDsId();
      auto channel = d.getChannel();

      for (uint16_t i = 0; i < d.nofSamples(); i++) {
        auto s = d.getSample(i);

        mNhits[solarId][dsId][channel] += 1;
        uint64_t N = mNhits[solarId][dsId][channel];

        double p0 = mPedestal[solarId][dsId][channel];
        double p = p0 + (s - p0) / N;
        mPedestal[solarId][dsId][channel] = p;

        double M0 = mNoise[solarId][dsId][channel];
        double M = M0 + (s - p0) * (s - p);
        mNoise[solarId][dsId][channel] = M;
      }

      if (mDebug) {
        std::cout << "solarId " << (int)solarId << "  dsId " << (int)dsId << "  ch " << (int)channel << "  nsamples " << d.nofSamples()
                    << "  nhits "<< mNhits[solarId][dsId][channel] << "  ped "<< mPedestal[solarId][dsId][channel] << "  noise " << mNoise[solarId][dsId][channel] << std::endl;
      }
    }
  }

  double PedestalProcessor::getPedestal(uint32_t solarId, uint32_t dsId, uint32_t channel) const
  {
    if (solarId >= MCH_NUMBER_OF_SOLAR || dsId >= 40 || channel >= 64) {
      return 0;
    }

    return mPedestal[solarId][dsId][channel];
  }

  double PedestalProcessor::getRms(uint32_t solarId, uint32_t dsId, uint32_t channel) const
  {
    if (solarId >= MCH_NUMBER_OF_SOLAR || dsId >= 40 || channel >= 64) {
      return 0;
    }

    double rms = std::sqrt(mNoise[solarId][dsId][channel] / mNhits[solarId][dsId][channel]);
    return rms;
  }

} // namespace o2::mch::calibration
