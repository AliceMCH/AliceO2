// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#ifndef Validation_h
#define Validation_h

#include <memory>
#include <fstream>
#include <stdio.h>

#include "MCHBase/Digit.h"
#include "MCHPreClustering/PreClusterFinder.h"
#include "MCHClustering/ClusteringForTest.h"
#include "TH1.h"

using namespace o2::mch;
using namespace std;

namespace o2 {

namespace mch {

 Double_t myMathieson2D(Double_t *x, Double_t *par);
Double_t myMathieson2D2hits(Double_t *x, Double_t *par);
 void myMath1hit(Double_t x, Double_t y);
void myMath2hits(Double_t x1, Double_t y1, Double_t x2, Double_t y2, Double_t chg1, Double_t chg2);
void ResidualsCOG();
void ResidualsCompare();
void ResidualsPlot(double yarray[], double resyfound[], double eyfound[], const int size);
void PlotWidthWrtCharge();
bool GradualAcceptance(int charge, double dice);
void PowFitHitsWrtChg();

class Validation
{
public:
  Validation();
  void PlotMathieson2D(TH1F* hchgpads, TH1F* hchgafter, TH1F* hchmaxafter, TH1F* hNbinsafter, TH1F* hNbinsX, TH1F* hNbinsXafter, TH1F* hNbinsY, TH1F* hNbinsYafter, TH1F* hMeanYbins, TH1F* hMeanbins, TH1F* hNhits0_300, TH1F* hNhits300_600, TH1F* hNhits600_1000, TH1F* hNhits1000_3000, TH1F* hYNhits0_300, TH1F* hYNhits300_600, TH1F* hYNhits600_1000, TH1F* hYNhits1000_3000, Double_t x, Double_t y, int nsamples, int SeedMath = 0);
  void InfoDE819b();
  void InfoDE819nb();
  std::vector<Clustering::Cluster> TestClustering();
    ssize_t getNumberOfDigits();
    void storeDigits(void* bufferPtr);
    
private:
    
    vector<double> lowxsb;
    vector<double> lowysb;
    vector<double> lowxsnb;
    vector<double> lowysnb;

    PreClusterFinder preClusterFinder;
    Clustering clustering;

    
    std::vector< std::unique_ptr<mch::Digit> > digits;
    mch::Digit* digitsBuffer;
    int nDigits;

};

}
}

#endif /* Validation_h */
