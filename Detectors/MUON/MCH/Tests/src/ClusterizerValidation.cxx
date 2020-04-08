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
#include "MCHBase/PreClusterBlock.h"
//#include "MCHPreClustering/PreClusterFinder.h"
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
    
    double xarray[6500]{0};
    double yarray[6500]{0};
    double chg[6500]{0};
    double resyfound[6500]{0};
    double eyfound[6500]{0};
    int count = 0;
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
    
    TCanvas *cchg = new TCanvas("cchg","Charge",0,0,600,600);
      TH1F *hchg = new TH1F("hchg", "Validation data cluster charge distribution without cuts or noise", 150, 0, 3000);
      TH1F *hchgafter = new TH1F("hchgafter", "Validation data cluster charge distribution with cuts and noise on bending plane, events with > 1 digit", 60, 0, 3000);
    hchgafter->SetLineColor(3);
    
    TCanvas *cNbinsY = new TCanvas("cNbinsY","NBinsY",0,0,600,600);
      TH1F *hNbinsY = new TH1F("hNbinsY", "Validation data cluster NBinsY distribution without cuts or noise", 20, 0.5, 20.5);
      TH1F *hNbinsYafter = new TH1F("hNbinsYafter", "Validation data NBinsY distribution with cuts and noise on bending plane, all events", 20, 0.5, 20.5);
    hNbinsYafter->SetLineColor(3);
    
    TCanvas *cMeanYbins = new TCanvas("cMeanYbins","cMeanYbins",0,0,600,600);
      TH1F *hMeanYbins = new TH1F("hMeanYbins", "Validation data mean Y pads wrt charge with cuts and noise on bending plane, events with > 1 digit", 60, 0, 3000);
    hMeanYbins->SetLineColor(3);
    
    TCanvas *cMeanbins = new TCanvas("cMeanbins","cMeanbins",0,0,600,600);
      TH1F *hMeanbins = new TH1F("hMeanbins", "Validation data mean pads wrt charge with cuts and noise on bending plane, events with > 1 digit", 60, 0, 3000);
    hMeanbins->SetLineColor(3);
    
    TCanvas *cNbinsX = new TCanvas("cNbinsX","NBinsX",0,0,600,600);
      TH1F *hNbinsX = new TH1F("hNbinsX", "Validation data cluster NBinsX distribution without cuts or noise", 20, 0.5, 20.5);
      TH1F *hNbinsXafter = new TH1F("hNbinsXafter", "Validation data NBinsX distribution with cuts and noise on bending plane, all events", 20, 0.5, 20.5);
    hNbinsXafter->SetLineColor(3);
    
    TCanvas *cNbins = new TCanvas("cNbins","NBins",0,0,600,600);
      TH1F *hNbinsafter = new TH1F("hNbinsafter", "Validation data NBins distribution with cuts and noise on bending plane, all events", 20, 0.5, 20.5);
    hNbinsafter->SetLineColor(3);
    
    TCanvas *cchmax = new TCanvas("cchmax","cchmax",0,0,600,600);
      TH1F *hchmaxafter = new TH1F("hchmaxafter", "Validation data cluster max pad charge distribution with cuts or noise, events with > 1 digit", 300, 0, 3000);
    hchmaxafter->SetLineColor(3);
    
    TCanvas *cchgpads = new TCanvas("cchgpads","cchgpads",0,0,600,600);
      TH1F *hchgpads = new TH1F("hchgpads", "Validation data pad charge distribution with cuts or noise, all events", 3000, 0, 3000);
    hchgpads->SetLineColor(3);
    
    cchg->Divide(2,1);
    cNbinsY->Divide(2,1);
    cNbinsX->Divide(2,1);
    cchg->cd(1);

    // Dependency of residuals with respect to y
//        for(int i=0; i<50; i++){
//            xarray[i] = 0;
//            yarray[i] = ygen->Uniform(0,0.5);
//            chg[i] = chggen->Uniform(20,2000);
//
//        }
    
    
    // Residuals distribution on the detector 819
//        for(int i=0; i<6500; i++){
//            yarray[i] = ygen->Uniform(-20,20);
//            xarray[i] = xgen->Uniform(-40,40);
//         //   chg[i] = chggen->Uniform(0,2000);
//            chg[i] = chggen->Landau(550,180);  //Distribution deduced from the fit on TB data 550 180 after
//            hchg->Fill(chg[i]);
//        }
    
      // Residuals distribution on the detector 819 like TB
            for(int i=0; i<6500; i++){
                yarray[i] = ygen->Uniform(11.7,12.9);
                xarray[i] = xgen->Uniform(8.75,11.25);
             //   chg[i] = chggen->Uniform(0,2000);
                chg[i] = chggen->Landau(565,195);  //Distribution deduced from the fit on TB data 550 180 after
                hchg->Fill(chg[i]);
            }

    
   //  Residuals distribution width with respect to fixed cluster charge, not really useful, since it is the same process as Residuals distribution but at fixed cluster charge
//        for(int i=0; i<2000; i++){
//            yarray[i] = ygen->Uniform(-20,20);
//            xarray[i] = xgen->Uniform(-40,40);
//            chg[i] = 20;
//        }
    
    
    hchg->GetXaxis()->SetTitle("Charge (ADC)");
    hchg->GetYaxis()->SetTitle("Count");
    hchg->Draw();

    cout << "\n\n==========\nGetting info for Bending plane\n\n" << endl;
   validation.InfoDE819b();
    cout << "\n\n==========\nGetting info for Non-Bending plane\n\n" << endl;
   validation.InfoDE819nb();
    for(int i=0; i<6500 ; i++){
    cout << "\n\n==========\nHit generation, histograms plotting and digitization\n\n" << endl;
   validation.PlotMathieson2D(hchgpads, hchgafter, hchmaxafter, hNbinsafter, hNbinsX, hNbinsXafter, hNbinsY, hNbinsYafter, hMeanYbins, hMeanbins, xarray[i], yarray[i], chg[i]);
    cout << "\n\n==========\nTesting the (pre)clustering\n\n" << endl;
        cout << "EVENT # " << i << endl;
   clusters = validation.TestClustering();
        if(clusters.size() == 1){
        resyfound[i] = yarray[i]-clusters[0].gety();
        eyfound[i] = clusters[0].getey();
            count++;
        }
    }
    
    hchgafter->Draw("SAME");
    cchg->cd(2);
    hchgafter->GetXaxis()->SetTitle("Charge (ADC)");
    hchgafter->GetYaxis()->SetTitle("Count");
    hchgafter->Draw();
    cchg->Update();
    cchg->Draw();
    
    cNbinsY->cd(1);
    hNbinsY->GetXaxis()->SetTitle("Number of pads on y");
    hNbinsY->GetYaxis()->SetTitle("Count");
    hNbinsY->Draw();
    hNbinsYafter->Draw("SAME");
    cNbinsY->cd(2);
    hNbinsYafter->GetXaxis()->SetTitle("Number of pads on y");
    hNbinsYafter->GetYaxis()->SetTitle("Count");
    hNbinsYafter->Draw();
    cNbinsY->Update();
    cNbinsY->Draw();
    
    cNbinsX->cd(1);
    hNbinsX->GetXaxis()->SetTitle("Number of pads on y");
    hNbinsX->GetYaxis()->SetTitle("Count");
    hNbinsX->Draw();
    hNbinsXafter->Draw("SAME");
    cNbinsX->cd(2);
    hNbinsXafter->GetXaxis()->SetTitle("Number of pads on x");
    hNbinsXafter->GetYaxis()->SetTitle("Count");
    hNbinsXafter->Draw();
    cNbinsX->Update();
    cNbinsX->Draw();
    
    cNbins->cd();
    hNbinsafter->GetXaxis()->SetTitle("Number of pads");
    hNbinsafter->GetYaxis()->SetTitle("Count");
    hNbinsafter->Draw();
    cNbins->Update();
    cNbins->Draw();
    
    cchmax->cd();
    hchmaxafter->GetXaxis()->SetTitle("Charge (ADC)");
    hchmaxafter->GetYaxis()->SetTitle("Count");
    hchmaxafter->Draw();
    cchmax->Update();
    cchmax->Draw();
    
    cchgpads->cd();
    hchgpads->GetXaxis()->SetTitle("Charge (ADC)");
    hchgpads->GetYaxis()->SetTitle("Count");
    hchgpads->Draw();
    cchgpads->Update();
    cchgpads->Draw();
    
    cMeanYbins->cd();
    hMeanYbins->GetXaxis()->SetTitle("Charge (ADC)");
    hMeanYbins->GetYaxis()->SetTitle("Mean number of pads on y");
    hMeanYbins->Divide(hchgafter);
    hMeanYbins->Draw();
    cMeanYbins->Update();
    cMeanYbins->Draw();
    
    cMeanbins->cd();
    hMeanbins->GetXaxis()->SetTitle("Charge (ADC)");
    hMeanbins->GetYaxis()->SetTitle("Mean number of pads");
    hMeanbins->Divide(hchgafter);
    hMeanbins->Draw();
    cMeanbins->Update();
    cMeanbins->Draw();
    

    cout << "\n\n==========\nValidation procedure terminated\n\n" << endl;

    // If want to plot the residual dependency wrt y and/or the residuals distribution for this run of events
   ResidualsPlot(yarray, resyfound, eyfound, count);

    
    //Useless
//  ResidualsCompare();


    // If want to plot the evolution of width with respect to cluster charge (values hardcoded in Validation.cxx)
    PlotWidthWrtCharge();
    
    PowFitHitsWrtChg();
    
   app.Run(kTRUE);
    
    return 0;
}

