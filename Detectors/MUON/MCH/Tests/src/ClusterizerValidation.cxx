// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

// This program serves as a closure test of clustering methods
// It calls function in the validation class to generate ideal hits and reconstruct them

#include <fstream>
#include <cstdio>

#include "Validation.h"

#include "TMath.h"
#include "TRandom.h"
#include "TH1.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TH2.h"
#include "TF2.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include <cmath>
#include <TApplication.h>


using namespace o2::mch;
using namespace std;

int main(int argc, char** argv){
    
    double xarray[200]{0};
    double yarray[200]{0};
    double chg[200]{0};
    double resyfound[200]{0};
    double eyfound[200]{0};
    
    Validation validation;
    std::vector<Clustering::Cluster> clusters;
    
    TApplication app ("app",&argc,argv);
    cout << "\n\n==========\nRunning the Validation procedure of (pre)clustering" << endl;
    
    TRandom *ygen = new TRandom(12345);
    TRandom *xgen = new TRandom(123456);
    TRandom *chggen = new TRandom(123);

    for(int i=0; i<200; i++){
//        yarray[i] = ygen->Uniform(0,0.5);
//        xarray[i] = 0;
        yarray[i] = ygen->Uniform(-20,20);
        xarray[i] = xgen->Uniform(-40,40);
     //   chg[i] = chggen->Uniform(20,2000);
        chg[i] = 1000;
    }

    cout << "\n\n==========\nGetting info for Bending plane\n\n" << endl;
   validation.InfoDE819b();
    cout << "\n\n==========\nGetting info for Non-Bending plane\n\n" << endl;
   validation.InfoDE819nb();

    for(int i=0; i<200 ; i++){
    cout << "\n\n==========\nHit generation, histograms plotting and digitization\n\n" << endl;
   validation.PlotMathieson2D(xarray[i], yarray[i], chg[i]);
    cout << "\n\n==========\nTesting the (pre)clustering\n\n" << endl;
        cout << "EVENT # " << i << endl;
   clusters = validation.TestClustering();
        resyfound[i] = yarray[i]-clusters[0].gety();
        eyfound[i] = clusters[0].getey();
    }

    cout << "\n\n==========\nValidation procedure terminated\n\n" << endl;


   ResidualsPlot(yarray, resyfound, eyfound, 200);

   // ResidualsCompare();

 //   PlotWidthWrtCharge();
    
   app.Run(kTRUE);
    
    return 0;
}

