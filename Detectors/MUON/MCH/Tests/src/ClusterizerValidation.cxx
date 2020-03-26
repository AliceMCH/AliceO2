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

#include "MCHBase/Digit.h"
//#include "MCHBase/DigitBlock.h"
#include "MCHPreClustering/PreClusterBlock.h"
#include "MCHPreClustering/PreClusterFinder.h"
#include "DigitsFileReader.h"
#include "MCHClustering/ClusteringForTest.h"
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
    
    // Declaration of arrays. The size of the arrays is hardcoded (it wouldn't accept an initialisation using a variable)
    // The size of the arrays should be the number of events you wish to simulate
    
    double xarray[50]{0};
    double yarray[50]{0};
    double chg[50]{0};
    double resyfound[50]{0};
    double eyfound[50]{0};
    
    Validation validation;
    std::vector<Clustering::Cluster> clusters;
    
    TApplication app ("app",&argc,argv);
    cout << "\n\n==========\nRunning the Validation procedure of (pre)clustering" << endl;
    
    //We generate the position of the hits and the charge of the clusters.
    //Depending on what you want to do you may need to fix some of these variables.
    
    /*
     USE CASES
     
     To check the dependency of residuals with respect to y (to see if there is a bias in the residuals):
     I would fix x to 0, get y uniformly between 0 and 0.5, get the charge uniformly between 20 and 2000.
     
     To see the residuals distribution on the detector 819:
     I would make x and y uniform on the detector 819 (so y from -20 to 20 and x from -40 to 40) and the charge uniformly between 20 and 2000.
     
     To see the evolution of residuals distribution width with respect to cluster charge:
     I would make x and y uniform on the detector 819 (so y from -20 to 20 and x from -40 to 40) and fix the charge to any given value I wish to look at.
     Then I would run the code, make a gaussian fit of the residuals distrubtion obtained, and hardcode the result in the dedicated function at the end of Validation.cxx
     */
    
    TRandom *ygen = new TRandom(12345);
    TRandom *xgen = new TRandom(123456);
    TRandom *chggen = new TRandom(123);

    // Dependency of residuals with respect to y
        for(int i=0; i<50; i++){
            xarray[i] = 0;
            yarray[i] = ygen->Uniform(0,0.5);
            chg[i] = chggen->Uniform(20,2000);
        }
    
    // Residuals distribution on the detector 819
//        for(int i=0; i<200; i++){
//            yarray[i] = ygen->Uniform(-20,20);
//            xarray[i] = xgen->Uniform(-40,40);
//            chg[i] = chggen->Uniform(20,2000);
//        }
    
   //  Residuals distribution width with respect to fixed cluster charge, not really useful, since it is the same process as Residuals distribution but at fixed cluster charge
//        for(int i=0; i<200; i++){
//            yarray[i] = ygen->Uniform(-20,20);
//            xarray[i] = xgen->Uniform(-40,40);
//            chg[i] = 5000;
//        }

    cout << "\n\n==========\nGetting info for Bending plane\n\n" << endl;
   validation.InfoDE819b();
    cout << "\n\n==========\nGetting info for Non-Bending plane\n\n" << endl;
   validation.InfoDE819nb();

    for(int i=0; i<50 ; i++){
    cout << "\n\n==========\nHit generation, histograms plotting and digitization\n\n" << endl;
   validation.PlotMathieson2D(xarray[i], yarray[i], chg[i]);
    cout << "\n\n==========\nTesting the (pre)clustering\n\n" << endl;
        cout << "EVENT # " << i << endl;
   clusters = validation.TestClustering();
        resyfound[i] = yarray[i]-clusters[0].gety();
        eyfound[i] = clusters[0].getey();
    }

    cout << "\n\n==========\nValidation procedure terminated\n\n" << endl;

    // If want to plot the residual dependency wrt y and/or the residuals distribution for this run of events
   ResidualsPlot(yarray, resyfound, eyfound, 50);

    
    //Useless
//  ResidualsCompare();


    // If want to plot the evolution of width with respect to cluster charge (values hardcoded in Validation.cxx)
 //   PlotWidthWrtCharge();
    
   app.Run(kTRUE);
    
    return 0;
}

