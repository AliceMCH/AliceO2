#include <fstream>
#include <vector>
#include <algorithm>
#include <functional>


#include "Validation.h"
#include "MCHMappingInterface/Segmentation.h"
#include "MCHMappingInterface/CathodeSegmentation.h"
#include "MCHMappingFactory/CreateSegmentation.h"
#include "DigitsFileReader.h"

#include "TMath.h"
#include "TRandom.h"
#include "TCanvas.h"
#include "TH2.h"
#include "TH1.h"
#include "TF2.h"
#include "TF1.h"
#include "TGraphPainter.h"
#include "TGraphErrors.h"
#include "TGraph2DErrors.h"
#include "TLegend.h"
#include <cmath>

using namespace o2::mch;
using namespace std;

namespace o2 {

namespace mch {

Validation::Validation()
{
  preClusterFinder.init();
};


//PART WITH DEFINITION OF THE MATHIESON AND THE IMPACT POINTS

// A single Mathieson distribution

Double_t myMathieson2D(Double_t *x, Double_t *par){
    Double_t pitch = 0.25;
    Float_t xx = x[0]-par[2];
    Float_t yy = x[1]-par[3];
    Double_t K2x = M_PI*0.5*(1-0.5*sqrt(par[0]));
    Double_t K1x = K2x * sqrt(par[0]) / 4 / atan(sqrt(par[0]));
    Double_t fx = K1x * ((1-pow(tanh(K2x*xx/pitch),2))/(1+ par[0]*pow(tanh(K2x*xx/pitch),2)));
    Double_t K2y = M_PI*0.5*(1-0.5*sqrt(par[1]));
    Double_t K1y = K2y * sqrt(par[1]) / 4 / atan(sqrt(par[1]));
    Double_t fy = K1y * ((1-pow(tanh(K2y*yy/pitch),2))/(1+ par[1]*pow(tanh(K2y*yy/pitch),2)));
    Double_t f = fx*fy;
    return f;
}

//________________________________________________________________________________________

//2 superimposed Mathieson distributions, defined by their respective centers.

Double_t myMathieson2D2hits(Double_t *x, Double_t *par){
    Double_t pitch = 0.25;
    Float_t xx1 = x[0]-par[2];
    Float_t yy1 = x[1]-par[3];
    Float_t xx2 = x[0]-par[4];
    Float_t yy2 = x[1]-par[5];
    Double_t K2x = M_PI*0.5*(1-0.5*sqrt(par[0]));
    Double_t K1x = K2x * sqrt(par[0]) / 4 / atan(sqrt(par[0]));
    Double_t fx1 = K1x * ((1-pow(tanh(K2x*xx1/pitch),2))/(1+ par[0]*pow(tanh(K2x*xx1/pitch),2)));
    Double_t K2y = M_PI*0.5*(1-0.5*sqrt(par[1]));
    Double_t K1y = K2y * sqrt(par[1]) / 4 / atan(sqrt(par[1]));
    Double_t fy1 = K1y * ((1-pow(tanh(K2y*yy1/pitch),2))/(1+ par[1]*pow(tanh(K2y*yy1/pitch),2)));
    Double_t fx2 = K1x * ((1-pow(tanh(K2x*xx2/pitch),2))/(1+ par[0]*pow(tanh(K2x*xx2/pitch),2)));
    Double_t fy2 = K1y * ((1-pow(tanh(K2y*yy2/pitch),2))/(1+ par[1]*pow(tanh(K2y*yy2/pitch),2)));
    Double_t f = par[6]*fx1*fy1 + par[7]*fx2*fy2;
    return f;
}

//________________________________________________________________________________________
void myMath1hit(Double_t x, Double_t y){

    // Function generating a ONE HIT Mathieson histogram
    
    TF2 *f1 = new TF2("myMath",myMathieson2D,-10,10,-10,10,4);
    f1->SetParameters(0.5840, 0.5085, x, y);
    f1->SetParNames("K3x", "K3y", "Mean x", "Mean y");
}
//0.5840 0.5085  x y
    
//________________________________________________________________________________________
void myMath2hits(Double_t x1, Double_t y1, Double_t x2, Double_t y2, Double_t chg1, Double_t chg2){
    
    // Function generating a TWO HITS Mathieson histogram
    
    TF2 *f1 = new TF2("myMath",myMathieson2D2hits,-10,10,-10,10,8);
    f1->SetParameters(0.5840, 0.5085, x1, y1, x2, y2, chg1, chg2);
    f1->SetParNames("K3x", "K3y", "Mean1 x", "Mean1 y", "Mean2 x", "Mean2 y", "Charge1", "Charge2");
    //f1->Draw("colz");
}





//PART WITH: - PLOTTING OF THE MATHIESON FOLLOWING THE DEFINITION WE USE (ONE OR TWO HITS) AND MAPPING OF A GIVEN DETECTOR (DE819)
//           - DIGIT CREATION

//________________________________________________________________________________________
void Validation::PlotMathieson2D(TH1F* hchgpads, TH1F* hchgafter, TH1F* hchmaxafter, TH1F* hNbinsafter, TH1F* hNbinsX, TH1F* hNbinsXafter, TH1F* hNbinsY, TH1F* hNbinsYafter, TH1F* hMeanYbins, TH1F* hMeanbins, Double_t x, Double_t y, int nsamples, int SeedMath){
    
    digits.clear();
    
    TH2F* hb(NULL);
    TH2F* hdice(NULL);
    TH2F* hnb(NULL);
    myMath1hit(x, y);
    int gPrintLevel = 0;
    
    o2::mch::mapping::CathodeSegmentation catsegb(819, kTRUE);
    o2::mch::mapping::CathodeSegmentation catsegnb(819, kFALSE);

    int nopadsb = catsegb.nofPads();
    int nopadsnb = catsegnb.nofPads();
    int padid;
    bool keepdigit = kTRUE;
    
    
    //Conversion des vecteurs en arrays bending
    
    double xlowsb[lowxsb.size()];
    double ylowsb[lowysb.size()];
    std::copy(lowxsb.begin(), lowxsb.end(), xlowsb);
    std::copy(lowysb.begin(), lowysb.end(), ylowsb);
    
    //Conversion des vecteurs en arrays non-bending
    
    double xlowsnb[lowxsnb.size()];
    double ylowsnb[lowysnb.size()];
    std::copy(lowxsnb.begin(), lowxsnb.end(), xlowsnb);
    std::copy(lowysnb.begin(), lowysnb.end(), ylowsnb);
    
    double noise;
    double dice;
    
    TRandom *keep = new TRandom(0);
    hdice = new TH2F("hdice","hdice",lowxsb.size()-1,xlowsb,lowysb.size()-1,ylowsb);
    for(int bin=0; bin<(lowxsb.size()+1)*(lowysb.size()+1); bin++){
        dice = keep->Uniform(0,1);
     //   cout << "[dicegen] dice = " << dice <<endl;
        hdice->AddBinContent(bin,dice);
   //     cout << "Saved in hdice: " << hdice->GetBinContent(bin)<<endl;
        
    }
    
    //Création et remplissage histogrammes bending
    
    cout << "Generating histograms bending and non-bending..."<< endl;

    hb = new TH2F("hb","hb",lowxsb.size()-1,xlowsb,lowysb.size()-1,ylowsb);
    gRandom->SetSeed(SeedMath);
    hb->FillRandom("myMath", /*nsamples=*/nsamples);
    
    //Plot to find the NBinsY
    std::vector<int> fullbins;
    int appearances = 0;
    int appearancesmax = 0;
    int column = 0;
    int columnmax = 0;
    int chgbin = 0;
    int chgsum = 0;
    for(int i=1; i<lowxsb.size()+1; i++){
        for(int j=1; j<lowysb.size()+1; j++){
            chgbin = hb->GetBinContent(i,j);
            if(chgbin != 0){
                cout << "Chgbin = " << chgbin << " in bin (" << i << "," << j <<")"<<endl;
                fullbins.push_back(i);
                fullbins.push_back(j);
                
            }
        }
    }
    for(int k; k < fullbins.size(); k += 2){
        if(fullbins[k] != column){
            column = fullbins[k];
            appearances = 0;
        }
        appearances++;
        if(appearances>appearancesmax){
            appearancesmax = appearances;
            columnmax = column;
        }
    }
    cout << "Column " << columnmax << " appears the most: " << appearancesmax << " times."<<endl;
    hNbinsY->Fill(appearancesmax);
    
    cout << "On a rempli hNbinsY" <<endl;
    
    //Plot to find the NBinsX
       std::vector<int> fullbinsX;
       int appearancesX = 0;
       int appearancesmaxX = 0;
       int line = 0;
       int linemax = 0;
       int chgbinX = 0;
       for(int i=1; i<lowysb.size()+1; i++){
           for(int j=1; j<lowxsb.size()+1; j++){
               chgbinX = hb->GetBinContent(j,i);
               if(chgbinX != 0){
                   cout << "Chgbin = " << chgbinX << " in bin (" << j << "," << i <<")"<<endl;
                   fullbinsX.push_back(i);
                   fullbinsX.push_back(j);
                   
               }
           }
       }
       for(int k; k < fullbinsX.size(); k += 2){
           if(fullbinsX[k] != line){
               line = fullbinsX[k];
               appearancesX = 0;
           }
           appearancesX++;
           if(appearancesX>appearancesmaxX){
               appearancesmaxX = appearancesX;
               linemax = line;
           }
       }
       cout << "Line " << linemax << " appears the most: " << appearancesmaxX << " times."<<endl;
       hNbinsX->Fill(appearancesmaxX);
       
       cout << "On a rempli hNbinsX" <<endl;
    
    
    TRandom *noisegen = new TRandom(321);
    for(int bin=0; bin<(lowxsb.size()+1)*(lowysb.size()+1); bin++){
        noise = noisegen->Gaus(0,1);
        hb->AddBinContent(bin, abs(noise));
    }
    
    //Plot to find the NBinsY after noise and cut
       std::vector<int> fullbinsafter;
        appearances = 0;
        appearancesmax = 0;
        column = 0;
        columnmax = 0;
        chgbin = 0;
    double dice1;
       for(int i=1; i<lowxsb.size()+1; i++){
           for(int j=1; j<lowysb.size()+1; j++){
               chgbin = hb->GetBinContent(i,j);
               dice1 = hdice->GetBinContent(i,j);
               keepdigit = GradualAcceptance(chgbin, dice1);
               if(chgbin > 10 && keepdigit){
                   cout << "Chgbin = " << chgbin << " in bin (" << i << "," << j <<")"<<endl;
                   fullbinsafter.push_back(i);
                   fullbinsafter.push_back(j);
                   
               }
           }
       }
       for(int k; k < fullbinsafter.size(); k += 2){
           if(fullbinsafter[k] != column){
               column = fullbinsafter[k];
               appearances = 0;
           }
           appearances++;
           if(appearances>appearancesmax){
               appearancesmax = appearances;
               columnmax = column;
           }
       }
    cout << "Column " << columnmax << " appears the most after noise and cut: " << appearancesmax << " times." <<endl;
       hNbinsYafter->Fill(appearancesmax);
    
    cout << "On a rempli hNbinsYafter" <<endl;
    
    
    //Plot to find the NBinsX after noise and cut
          std::vector<int> fullbinsafterX;
           appearancesX = 0;
           appearancesmaxX = 0;
           line = 0;
           linemax = 0;
           chgbinX = 0;
          for(int i=1; i<lowysb.size()+1; i++){
              for(int j=1; j<lowxsb.size()+1; j++){
                  chgbinX = hb->GetBinContent(j,i);
                  dice1 = hdice->GetBinContent(j,i);
                  keepdigit = GradualAcceptance(chgbinX, dice1);
                  if(chgbinX > 10 && keepdigit){
                      cout << "ChgbinX = " << chgbinX << " in bin (" << j << "," << i <<")"<<endl;
                      fullbinsafterX.push_back(i);
                      fullbinsafterX.push_back(j);
                      
                  }
              }
          }
          for(int k; k < fullbinsafterX.size(); k += 2){
              if(fullbinsafterX[k] != line){
                  line = fullbinsafterX[k];
                  appearancesX = 0;
              }
              appearancesX++;
              if(appearancesX>appearancesmaxX){
                  appearancesmaxX = appearancesX;
                  linemax = line;
              }
          }
       cout << "Line " << linemax << " appears the most after noise and cut: " << appearancesmaxX << " times." <<endl;
          hNbinsXafter->Fill(appearancesmaxX);
       
       cout << "On a rempli hNbinsXafter" <<endl;
       
    
//    TCanvas *c2 = new TCanvas("c2", "The road to PlotDorado",
//                              200,10,600,400);
//    c2->Divide(1,2);
//    c2->cd(1);
//
//    hb->SetTitle("Mathieson 2D Distribution - Bending plane - DE819");
//    hb->Draw("LEGO2");
//
//    c2->Update();
//    c2->Draw();
//
//    TCanvas *c3 = new TCanvas("c3", "The road to PlotDorado",
//    200,10,600,400);
//    c3->Divide(1,2);
//    c3->cd(1);
//
//    hb->Draw("colz");
//
//    c3->Update();
//    c3->Draw();
    
    //Création et remplissage histogrammes non-bending
    
    hnb = new TH2F("hnb","hnb",lowxsnb.size()-1,xlowsnb,lowysnb.size()-1,ylowsnb);
    hnb->FillRandom("myMath", /*nsamples=*/nsamples);
    
    for(int bin=0; bin<(lowxsnb.size()+1)*(lowysnb.size()+1); bin++){
        noise = noisegen->Gaus(0,1);
        hnb->AddBinContent(bin, abs(noise));
    }
    
//    c2->cd(2);
//    hnb->SetTitle("Mathieson 2D Distribution - Non-Bending plane - DE819");
//    hnb->Draw("LEGO2");
//
//    c2->Update();
//    c2->Draw();
//
//    c3->cd(2);
//    hnb->Draw("colz");
//
//    c3->Update();
//    c3->Draw();
    
    
    int nbbinsb = (lowxsb.size()-1)*(lowysb.size()-1);
    int nbbinsnb = (lowxsnb.size()-1)*(lowysnb.size()-1);
    double bincontent;
    double binxcent;
    double binycent;
    double charge[nopadsb + nopadsnb];
    double dicecontent;
    double dice2[nopadsb + nopadsnb];
    
    for(int i=0; i<nopadsb + nopadsnb; i++ ){
        charge[i] = 0.;
    }
    for(int i=0; i<nopadsb + nopadsnb; i++ ){
        dice2[i] = 0.;
    }
    
    if(gPrintLevel > 1){
        cout << "Number of bending bins = " << lowxsb.size()-1 << " x " << lowysb.size()-1 << " = " << nbbinsb << endl;
        cout << "Number of non-bending bins = " << lowxsnb.size()-1 << " x " << lowysnb.size()-1 << " = " << nbbinsnb << endl;
    }
    
    
    //Remplissage digits bending plane
    
    cout << "Filling digits from bending plane..." << endl;
    
    for(int binx = 1; binx<lowxsb.size(); binx++){
        for(int biny = 1; biny<lowysb.size(); biny++){
            bincontent = hb->GetBinContent(binx,biny);
            dicecontent = hdice->GetBinContent(binx,biny);
            binxcent = (hb->GetXaxis())->GetBinCenter(binx);
            binycent = (hb->GetYaxis())->GetBinCenter(biny);
            
//            if((binx*biny)%4 == 0){
//                cout << "Content of bending bin (binx= " << binx << ", biny= " << biny << "): " << bincontent << endl;
//                cout << "Center of bending bin (binx= " << binx << ", biny= " << biny << "): (" << binxcent << "," << binycent << ")" << endl;
//            }
            
            padid = catsegb.findPadByPosition(binxcent,binycent);
            charge[padid] += bincontent;
            dice2[padid] += dicecontent;
            
//            if(charge[padid] > 2){
//               cout << "Content of bending bin (binx= " << binx << ", biny= " << biny << "): " << bincontent << endl;
//               cout << "Center of bending bin (binx= " << binx << ", biny= " << biny << "): (" << binxcent << "," << binycent << ")" << endl;
//               cout << "Padid: " << padid << "  Charge: " << charge[padid] <<endl;
//           }
        }
    }
    
    //Remplissage digits non-bending plane
    
    cout << "Filling digits from non-bending plane..." << endl;
        
        for(int binx = 1; binx<lowxsnb.size(); binx++){
            for(int biny = 1; biny<lowysnb.size(); biny++){
                bincontent = hnb->GetBinContent(binx,biny);
                binxcent = (hnb->GetXaxis())->GetBinCenter(binx);
                binycent = (hnb->GetYaxis())->GetBinCenter(biny);
                
                padid = catsegnb.findPadByPosition(binxcent,binycent) + nopadsb;
                charge[padid] += bincontent;
                if(gPrintLevel > 2){
                    if(charge[padid] > 2){
                        cout << "Content of non-bending bin (binx= " << binx << ", biny= " << biny << "): " << bincontent << endl;
                        cout << "Center of non-bending bin (binx= " << binx << ", biny= " << biny << "): (" << binxcent << "," << binycent << ")" << endl;
                        cout << "Padid: " << padid << "  Charge: " << charge[padid] <<endl;
                    }
                }
            }
        }
    
    if(gPrintLevel > 1){
        cout << "Nombre de bins bending : " << nbbinsb << " Nombre de pads bending : " << nopadsb << endl;
        cout << "Nombre de bins non-bending : " << nbbinsnb << " Nombre de pads non-bending : " << nopadsnb << endl;
    }
    
    int chargeevent = 0;
    int chpadmax = 0;
    int nbleftdigits = 0;
    
    for(int i=0; i<nopadsb + nopadsnb; i++){
        keepdigit = GradualAcceptance(charge[i], dice2[i]);
        if(charge[i] > 10 && keepdigit){  //Couper le bruit type (6ADC = 2*seuil de 3ADC, seuil de 3adc vient de 3*bruit/0.8)
            if(i<nopadsb){
                chargeevent += charge[i];
                hchgpads->Fill(charge[i]);
                if(charge[i]>chpadmax){
                    chpadmax = charge[i];
                }
                nbleftdigits++;
        
              //  cout << "Un digit est en cours de création..." << endl;

                digits.push_back( std::make_unique<mch::Digit>() );
                mch::Digit* digit = digits.back().get();

                   digit->setDetID(819);
                   digit->setPadID(i);
                   digit->setADC(charge[i]);
                if(gPrintLevel > 1){
                    cout << "Digit padid " << digit->getPadID() << " et ADC " << digit->getADC() << endl;
                }
            }
        }
    }
    
    cout << "This event had " << nbleftdigits << " digits left after noise and cuts" <<endl;
    hNbinsafter->Fill(nbleftdigits);
    
    if(nbleftdigits > 1){
        hchgafter->Fill(chargeevent);
        hMeanYbins->Fill(chargeevent, appearancesmax);
        hMeanbins->Fill(chargeevent, nbleftdigits);
        hchmaxafter->Fill(chpadmax);
    }
    
    delete hb;
    delete hnb;
    delete hdice;
}


//________________________________________________________________________________________
bool GradualAcceptance(int charge, double dice){
    //Following the charge, we will keep random events to get the pad charge distribution we want as in TB. To replicate the cutting process.
    
    double probabilitiestokeep[9] = {0.002, 0.0047, 0.0125, 0.0718, 0.1875, 0.25, 0.4333, 0.72, 0.84};
    
    
    if(charge <= 10){
           return kFALSE;
       }
    
    else if(charge >= 20){
        return kTRUE;
    }
    
    else if(dice < probabilitiestokeep[charge-11]){
        cout << "[GradualAcceptance] The charge is " << charge << " and the dice tells " << dice << " compared to " << probabilitiestokeep[charge-11] << "... KEEP" << endl;
        return kTRUE;
    }
    
    else{
        cout << "[GradualAcceptance] The charge is " << charge << " and the dice tells " << dice << " compared to " << probabilitiestokeep[charge-11] << "... THROW" <<endl;
        return kFALSE;
    }

}

//PART GETTING THE MAPPING INFO FOR THE TWO CATHODES OF DETECTOR 819 TO MAKE NICE HISTOGRAMS BASED ON SIZE AND POSITION OF PADS

//________________________________________________________________________________________
void Validation::InfoDE819b(){
    
    int detElemId =819;
    bool isbending = kTRUE;
    try {
      const o2::mch::mapping::Segmentation& segment = o2::mch::mapping::segmentation(detElemId);
      const o2::mch::mapping::CathodeSegmentation& catseg = segment.bending();
    //o2::mch::mapping::CathodeSegmentation catseg(detElemId, isbending);

    int nopads = catseg.nofPads();
    int gPrintLevel = 0;
    
    double padposx;
    double padposy;
    double padsizex;
    double padsizey;
    double lowerx;
    double lowery;
    lowxsb.clear();
    lowysb.clear();
    int valuesx = 0;
    int valuesy = 0;
    
    if(gPrintLevel > 1){
        cout << "There are " << nopads << " pads on bending plane" << endl;
    }
    
    for(int catPadindex = 0; catPadindex<nopads; catPadindex++){
        padposx = catseg.padPositionX(catPadindex);
        padposy = catseg.padPositionY(catPadindex);
        padsizex = catseg.padSizeX(catPadindex);
        padsizey = catseg.padSizeY(catPadindex);
        lowerx = padposx - (padsizex/2);
        lowery = padposy - (padsizey/2);
        
        if(!(find(lowxsb.begin(), lowxsb.end(), lowerx) != lowxsb.end())){
            lowxsb.push_back(lowerx);
            valuesx++;
        }
        if(!(find(lowysb.begin(), lowysb.end(), lowery) != lowysb.end())){
            lowysb.push_back(lowery);
            valuesy++;
        }
        
    //    if(catPadindex%10==0){
    //        cout << "catPadindex " << catPadindex << "\n"
    //        << "Position: (" << padposx << "," << padposy
    //        << ")\n"
    //        << "Size: (" << padsizex << "," << padsizey
    //        << ")\n";
    //        << "LowerLimits: (" << lowerx << "," << lowery
    //        << ")" << endl;
    //    }
        
    }
    
    if(gPrintLevel > 1){
        cout << "Il doit y avoir " << valuesx << " valeurs dans le vecteur lowxsb"
        << "et " << valuesy << " valeurs dans le vecteur lowysb" << endl;
        cout << "Il y a " << lowxsb.size() << " valeurs dans le vecteur lowxsb"
        << "et " << lowysb.size() << " valeurs dans le vecteur lowysb" << endl;
    }
    
    cout << "Rangement des bords de pads par ordre croissant ..." << endl;
    
    sort(lowxsb.begin(), lowxsb.end());
    sort(lowysb.begin(), lowysb.end());
    
    cout << "Suppression des doublons par erreur d'arrondi ..." << endl;
    
    auto lastxb = std::unique(lowxsb.begin(), lowxsb.end(),
    [](double l, double r) { return std::abs(l - r) < 0.001; });
    lowxsb.erase(lastxb, lowxsb.end());
    
    auto lastyb = std::unique(lowysb.begin(), lowysb.end(),
    [](double l, double r) { return std::abs(l - r) < 0.001; });
    lowysb.erase(lastyb, lowysb.end());
    
    if(gPrintLevel > 1){
        cout << "Il y a maintenant " << lowxsb.size() << " valeurs dans le vecteur lowxsb"
        << " et " << lowysb.size() << " valeurs dans le vecteur lowysb" << endl;
    }
    
    cout << "Ajout de la borne sup" << endl;
    lowxsb.push_back(2*lowxsb[lowxsb.size()-1]-lowxsb[lowxsb.size()-2]);
    lowysb.push_back(2*lowysb[lowysb.size()-1]-lowysb[lowysb.size()-2]);
    
    if(gPrintLevel > 1){
        cout << "xlowsb:" << lowxsb[0] << "," << lowxsb[1] << "..." << lowxsb[valuesx-1] << "," << lowxsb[valuesx] << endl;
        cout << "ylowsb:" << lowysb[0] << "," << lowysb[1] << "..." << lowysb[valuesy-1] << "," << lowysb[valuesy] << endl;
    }
    
//    for(int i=0; i<lowxsb.size(); i++){
//        cout << lowxsb[i] << endl;
//    }
//    for(int i=0; i<lowysb.size(); i++){
//        cout << lowysb[i] << endl;
//    }
    }
    catch (const std::exception& e) {
      std::cout<<"[InfoDE819b] cannot create segmentation for DE "<<detElemId<<std::endl;
        return;
    }
}

//________________________________________________________________________________________
void Validation::InfoDE819nb(){
    
    int detElemId =819;
    bool isbending = kFALSE;
    try {
      const o2::mch::mapping::Segmentation& segment = o2::mch::mapping::segmentation(detElemId);
      const o2::mch::mapping::CathodeSegmentation& catseg = segment.nonBending();
    //o2::mch::mapping::CathodeSegmentation catseg(detElemId, isbending);

    int nopads = catseg.nofPads();
    int gPrintLevel = 0;
    
    double padposx;
    double padposy;
    double padsizex;
    double padsizey;
    double lowerx;
    double lowery;
    lowxsnb.clear();
    lowysnb.clear();
    int valuesx = 0;
    int valuesy = 0;
    
    if(gPrintLevel > 1){
        cout << "There are " << nopads << " pads on non-bending plane" << endl;
    }
    
    for(int catPadindex = 0; catPadindex<nopads; catPadindex++){
        padposx = catseg.padPositionX(catPadindex);
        padposy = catseg.padPositionY(catPadindex);
        padsizex = catseg.padSizeX(catPadindex);
        padsizey = catseg.padSizeY(catPadindex);
        lowerx = padposx - (padsizex/2);
        lowery = padposy - (padsizey/2);
        
        if(!(find(lowxsnb.begin(), lowxsnb.end(), lowerx) != lowxsnb.end())){
            lowxsnb.push_back(lowerx);
            valuesx++;
        }
        if(!(find(lowysnb.begin(), lowysnb.end(), lowery) != lowysnb.end())){
            lowysnb.push_back(lowery);
            valuesy++;
        }
        
//        if(catPadindex%10==0){
//            cout << "catPadindex " << catPadindex << "\n"
//            << "Position: (" << padposx << "," << padposy
//            << ")\n"
//            << "Size: (" << padsizex << "," << padsizey
//            << ")\n"
//            << "LowerLimits: (" << lowerx << "," << lowery
//            << ")" << endl;
//        }
        
    }
    
    if(gPrintLevel > 1){
        cout << "Il doit y avoir " << valuesx << " valeurs dans le vecteur lowxsnb"
        << "et " << valuesy << " valeurs dans le vecteur lowysnb" << endl;
        cout << "Il y a " << lowxsnb.size() << " valeurs dans le vecteur lowxsnb"
        << "et " << lowysnb.size() << " valeurs dans le vecteur lowysnb" << endl;
    }
    
    cout << "Rangement des bords de pads par ordre croissant ..." << endl;
    
    sort(lowxsnb.begin(), lowxsnb.end());
    sort(lowysnb.begin(), lowysnb.end());
    
    cout << "Suppression des doublons par erreur d'arrondi ..." << endl;
    
    auto lastxnb = std::unique(lowxsnb.begin(), lowxsnb.end(),
    [](double l, double r) { return std::abs(l - r) < 0.001; });
    lowxsnb.erase(lastxnb, lowxsnb.end());
    
    auto lastynb = std::unique(lowysnb.begin(), lowysnb.end(),
    [](double l, double r) { return std::abs(l - r) < 0.001; });
    lowysnb.erase(lastynb, lowysnb.end());
    
    if(gPrintLevel > 1){
        cout << "Il y a maintenant " << lowxsnb.size() << " valeurs dans le vecteur lowxsnb"
        << " et " << lowysnb.size() << " valeurs dans le vecteur lowysnb" << endl;
    }
    
    cout << "Ajout de la borne sup" << endl;
    lowxsnb.push_back(2*lowxsnb[lowxsnb.size()-1]-lowxsnb[lowxsnb.size()-2]);
    lowysnb.push_back(2*lowysnb[lowysnb.size()-1]-lowysnb[lowysnb.size()-2]);
    
//    cout << "xlowsnb:" << lowxsnb[0] << "," << lowxsnb[1] << "..." << lowxsnb[lowxsnb.size()-2] << "," << lowxsnb[lowxsnb.size()-1] << endl;
//    cout << "ylowsnb:" << lowysnb[0] << "," << lowysnb[1] << "..." << lowysnb[lowysnb.size()-2] << "," << lowysnb[lowysnb.size()-1] << endl;
    
//    for(int i=0; i<lowxsnb.size(); i++){
//        cout << lowxsnb[i] << endl;
//    }
//    for(int i=0; i<lowysnb.size(); i++){
//        cout << lowysnb[i] << endl;
//    }
    }
    catch (const std::exception& e) {
      std::cout<<"[InfoDE819b] cannot create segmentation for DE "<<detElemId<<std::endl;
        return;
    }
}








//PART RUNNING THE PRECLUSTERING AND CLUSETRING FOLLOWING CHOSEN METHOD, BASED ON DIGITS OBTAINED FROM ABOVE

//________________________________________________________________________________________
std::vector<Clustering::Cluster> Validation::TestClustering(){
    
    cout << "Filling buffer of digits..." << endl;
    
    Digit* digitsBuffer = NULL;
    nDigits = getNumberOfDigits();
    digitsBuffer = (mch::Digit*)realloc(digitsBuffer, sizeof(mch::Digit) * nDigits);
    storeDigits(digitsBuffer);

      std::vector<mch::Clustering::Cluster> clusters(0);

        // load the digits from the memory buffer and run the pre-clustering phase
        preClusterFinder.reset();
        preClusterFinder.loadDigits(digitsBuffer, nDigits);
        preClusterFinder.run();

        // get the preclusters and associated digits
        std::vector<Digit> digits(0);
        std::vector<PreClusterStruct> preClusters(0);
        preClusterFinder.getPreClusters(preClusters, digits);

          printf("\n\n==========\nRunning Clustering\n\n");
    
    
    
    // Uncomment the desired clustering method
          
    
    //To run COG Clustering
       clustering.runFinderCOG(preClusters, clusters);
//        printf("Number of clusters obtained and saved: %lu\n", clusters.size());
    
    
    //To run Mathieson fit Clustering
 //      clustering.runFinderSimpleFit(preClusters, clusters);
          
    
    //To run Gaussian fit Clustering
 //         clustering.runFinderGaussianFit(preClusters, clusters);
    
    
    //To run Double Gaussian fit Clustering
  //  clustering.runFinderDoubleGaussianFit(preClusters, clusters);
    
    delete digitsBuffer;
    
    return clusters;

}

//________________________________________________________________________________________
ssize_t Validation::getNumberOfDigits()
{
  return digits.size();
}

//________________________________________________________________________________________
void Validation::storeDigits(void* bufferPtr)
{
  mch::Digit* ptr = (mch::Digit*)bufferPtr;
  for(unsigned int di = 0; di < digits.size(); di++) {

    memcpy(ptr, digits[di].get(), sizeof(mch::Digit));
    ptr += 1;
  }
}










//PART PLOTTING RESULTS OF RESIDUALS MEASUREMENTS OBTAINED MANUALLY (NOT IMPORTANT HERE)

//________________________________________________________________________________________
void ResidualsCOG(){
    
    const Int_t n = 15;
    
    Double_t xinput[n]  = {0, 2.5, 5, 7.5, 10, 0, 2.5, 5, 7.5, 10, 0, 2.5, 5, 7.5, 10};
    Double_t exinput[n] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    Double_t yinput[n]  = {0, 0, 0, 0, 0, 0.25, 0.25, 0.25, 0.25, 0.25, 0.5, 0.5, 0.5, 0.5, 0.5};
    Double_t eyinput[n] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    Double_t restot[n]  = {0.007840347,  0.002234084,  0.006410746,  0.003586554,  0.004374705,  0.006152688,  0.00139723,  0.005710998,  0.002378442,  0.000750866,  0.002745256,  0.003005745,  0.007499117,  0.002986651,  0.006803757};
    Double_t erestot[n] = {0.005222712,  0.004080179,  0.002880783,  0.0037952,  0.00523312,  0.005242154,  0.00214019,  0.006411414,  0.002580706,  0.005788876,  0.005806547,  0.00273805,  0.004339903,  0.003174934,  0.005813516};
    
    TCanvas *cerr = new TCanvas("cerr","Graph2DErrors example",0,0,600,600);
    
    TGraph2DErrors *dte = new TGraph2DErrors(n, xinput, yinput, restot, exinput, eyinput, erestot);
    dte->SetTitle("Residuals COG - DE819 on a pad of Bending plane");
    dte->SetFillColor(29);
    dte->SetMarkerSize(0.8);
    dte->SetMarkerStyle(20);
    dte->SetMarkerColor(kRed);
    dte->SetLineColor(kBlue-3);
    dte->SetLineWidth(2);
    dte->Draw("err p0");
    
    cerr->Update();
    cerr->Draw();
    
}

void ResidualsCompare(){
    const Int_t n = 11;
    
    //Residuals on Y for x=0 y between 0 and 0.5 cm on DE 819
    
    TCanvas *c10 = new TCanvas("c10", "The road to PlotDorado",
    200,10,600,400);
        
                Float_t meaninput[n]  = {0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5};
                Float_t emeaninput[n] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

                Float_t resgaus[n]  = {-0.00413, -0.00140, -0.01742, -0.01522, 0.00244, -0.00323, 0.00756, 0.00685, 0.01445, 0.00716, 0.00232};
                Float_t eresgaus[n] = {0.00570, 0.00577, 0.00561, 0.00559, 0.00551, 0.00560, 0.00549, 0.00605, 0.00567, 0.00569, 0.00568};

//                Float_t rescog[n]  = {-0.0114, -0.0548, -0.0840, -0.0834, -0.0432, 0.0012, 0.0350, 0.0810, 0.0775, 0.0491, 0.0018};
//                Float_t erescog[n] = {0.0061, 0.0101, 0.0053, 0.0057, 0.0068, 0.0094, 0.0080, 0.0074, 0.0042, 0.0116, 0.0061};
    
                Float_t resdoublegaus[n]  = {-0.00046, -0.00139, 0.00185, -0.00710, -0.00106, 0.00444, -0.00398, -0.00471, -0.00572, -0.00414, -0.00381};
                Float_t eresdoublegaus[n] = {0.00494, 0.00503, 0.00535, 0.00580, 0.00672, 0.00698, 0.00636, 0.00570, 0.00535, 0.00495, 0.00493};
    
                Float_t resmathieson[n]  = {0.00636, 0.00010, 0.00113, 0.00003, -0.00110, 0.00513, 0.00400, -0.00413, 0.00327, 0.00645, -0.00574};
                Float_t eresmathieson[n] = {0.00495, 0.00504, 0.00551, 0.00589, 0.00637, 0.00633, 0.00649, 0.00605, 0.00535, 0.00499, 0.00479};

                TGraphErrors *gr1 = new TGraphErrors(n,meaninput,resgaus,emeaninput,eresgaus);
                TGraphErrors *gr2 = new TGraphErrors(n,meaninput,resdoublegaus,emeaninput,eresdoublegaus);
             //   TGraphErrors *gr3 = new TGraphErrors(n,meaninput,rescog,emeaninput,erescog);
                TGraphErrors *gr4 = new TGraphErrors(n,meaninput,resmathieson,emeaninput,eresmathieson);
                TF1 *gr0 = new TF1("gr0","0",-1,1);

                gr1->SetTitle("Residuals wrt input y");
                gr1->SetMarkerColor(4);
                gr1->SetLineColor(4);
                gr1->SetMarkerStyle(8);
                gr1->GetXaxis()->SetTitle("y input (cm)");
                gr1->GetYaxis()->SetTitle("Residual (cm)");
                gr1->Draw("AP");
                gr2->SetMarkerColor(3);
                gr2->SetLineColor(3);
                gr2->SetMarkerStyle(8);
                gr2->Draw("P SAME");
//                gr3->SetMarkerColor(2);
//                gr3->SetLineColor(2);
//                gr3->SetMarkerStyle(8);
//                gr3->Draw("P SAME");
                gr4->SetMarkerColor(1);
                gr4->SetLineColor(1);
                gr4->SetMarkerStyle(8);
                gr4->Draw("P SAME");
                gr0->SetLineStyle(2);
                gr0->SetLineColor(1);
                gr0->Draw("SAME");

                auto legend = new TLegend(0.1,0.7,0.3,0.9);
                legend->SetHeader("Residuals"); // option "C" allows to center the header
                legend->AddEntry(gr1,"Residuals Gaussian fit","lep");
                legend->AddEntry(gr2,"Residuals Double Gaussian fit","lep");
//                legend->AddEntry(gr3,"Residuals COG","lep");
                legend->AddEntry(gr4,"Residuals Mathieson","lep");
                legend->Draw();
    
    c10->Update();
    c10->Draw();
    
    
}

// Function plotting the residuals as a function of the position y of the hit and the residuals distribution

void ResidualsPlot(double yarray[], double resyfound[], double eyfound[], int size){
    
    const int n = size;
    
    Double_t eyinput[n];
    for(int i=0; i<n; i++){
        eyinput[i]=0;
    }
    
    
    TCanvas *cerr = new TCanvas("cerr","GraphErrors example",0,0,600,600);
    
    TGraphErrors *gr1 = new TGraphErrors(n, yarray, resyfound, eyinput, eyfound);
    gr1->SetTitle("Residuals wrt input y");
    gr1->SetMarkerColor(4);
    gr1->SetLineColor(4);
    gr1->SetMarkerStyle(8);
    gr1->GetXaxis()->SetTitle("y input (cm)");
    gr1->GetYaxis()->SetTitle("Residual (cm)");
    gr1->Draw("AP");
    cerr->Update();
    cerr->Draw();
    
    
    TCanvas *cbell = new TCanvas("cbell","Bell",0,0,600,600);
    TH1F *h1 = new TH1F("h1", "Residuals distribution", 50, -0.1, 0.1);
    for(int i=0; i<n; i++){
        h1->Fill(resyfound[i]);
    }
    h1->GetXaxis()->SetTitle("Residual y (cm)");
    h1->GetYaxis()->SetTitle("Count");
    h1->Draw();
    cbell->Update();
    cbell->Draw();
    
}

// Function plotting the width of the residuals distribution with respect to charge of the clusters for different clustering methods. Hardcoded vectors found by simulations.

void PlotWidthWrtCharge(){
    
    const Int_t n = 8;
    
    Double_t chginput[n]  = {20, 50, 100, 200, 500, 1000, 2000, 5000};
    Double_t echginput[n] = {0, 0, 0, 0, 0, 0, 0, 0};
    
    // RESULTATS AVEC 50 EVTS
    
//    Double_t widthMathieson[n]  = {0.0724216, 0.0473039, 0.0237125, 0.0195968, 0.0136211, 0.00922341, 0.00550195, 0.00424860};
//    Double_t ewidthMathieson[n] = {0.00823728, 0.00253988, 0.000958042, 0.000588419, 0.000349494, 0.000168508, 0.0000935699, 0.0000546348};
//    Double_t widthDoubleGauss[n]  = {0.0809867, 0.0405142, 0.0245105, 0.0198315, 0.0130453, 0.00867564, 0.00585070, 0.00422690};
//    Double_t ewidthDoubleGauss[n] = {0.0110478, 0.00196554, 0.000973368, 0.000605342, 0.000327329, 0.000196655, 0.0000988013, 0.0000564364};
//    Double_t widthGauss[n]  = {0.115848, 0.0386573, 0.0268345, 0.0240388, 0.0184959, 0.0164656, 0.0157473, 0.0171971};
//    Double_t ewidthGauss[n] = {0.0261366, 0.00189652, 0.00106311, 0.000851333, 0.000518058, 0.000461726, 0.000396683, 0.000548232};
//    Double_t widthCOG[n]  = {0.111975, 0.0470066, 0.0272162, 0.0258470, 0.0250620, 0.0212521, 0.0182168, 0.0266299};
//    Double_t ewidthCOG[n] = {0.0234304, 0.00250890, 0.00103235, 0.000912423, 0.000890655, 0.000643082, 0.000500731, 0.00143591};

    // RESULTATS AVEC 2000 EVTS
    
        Double_t widthCOG[n]  = {7.60227E-02, 3.91475E-02, 2.70075E-02, 2.47138E-02, 2.36353E-02, 1.95029E-02, 1.81684E-02, 1.80579E-02};
        Double_t ewidthCOG[n] = {8.74033E-04, 1.77986E-04, 9.88354E-05, 8.31474E-05, 6.91609E-05, 5.04850E-05, 4.32215E-05, 4.11040E-05};
    
        Double_t widthGauss[n]  = {8.04075E-02, 4.18631E-02, 2.88616E-02, 2.23938E-02, 1.83937E-02, 1.64864E-02, 1.49507E-02, 1.41666E-02};
        Double_t ewidthGauss[n] = {1.00424E-03, 1.93914E-04, 1.03890E-04, 7.12706E-05, 5.00895E-05, 3.87129E-05, 3.38022E-05, 2.98527E-05};
    
        Double_t widthDoubleGauss[n]  = {7.82095E-02, 4.09391E-02, 2.75198E-02, 2.12291E-02, 1.26785E-02, 8.64849E-03, 5.82560E-03, 3.99766E-03};
        Double_t ewidthDoubleGauss[n] = {9.27416E-04, 1.86571E-04, 9.97959E-05, 6.53720E-05, 3.18873E-05, 1.71273E-05, 9.97586E-06, 5.41483E-06};
    
        Double_t widthMathieson[n]  = {8.36689E-02, 3.94728E-02, 2.70987E-02, 1.85719E-02, 1.22286E-02, 8.33011E-03, 5.73253E-03, 3.71460E-03};
        Double_t ewidthMathieson[n] = {1.12688E-03, 1.76920E-04, 9.55805E-05, 5.20617E-05, 2.77115E-05, 1.67013E-05, 9.26347E-06, 5.06832E-06};
    
        Double_t widthMathieson01[n]  = {8.05518E-02, 4.62265E-02, 3.63358E-02, 2.84678E-02, 2.41877E-02, 2.34258E-02, 2.32844E-02, 2.41844E-02};
        Double_t ewidthMathieson01[n] = {1.01157E-03, 2.47385E-04, 1.51532E-04, 1.01928E-04, 7.55642E-05, 6.69759E-05, 6.37582E-05, 6.61309E-05};
    
        Double_t widthMathieson02[n]  = {8.04581E-02, 4.53607E-02, 3.16927E-02, 2.26150E-02, 1.61596E-02, 1.51937E-02, 1.44840E-02, 1.45222E-02};
        Double_t ewidthMathieson02[n] = {1.01779E-03, 2.27123E-04, 1.22195E-04, 7.31616E-05, 4.62418E-05, 3.82145E-05, 3.24152E-05, 3.08703E-05};
    
        Double_t widthMathieson03[n]  = {7.28293E-02, 4.22823E-02, 2.89688E-02, 1.96324E-02, 1.38605E-02, 1.12013E-02, 1.00848E-02, 8.96331E-03};
        Double_t ewidthMathieson03[n] = {7.66027E-04, 2.02240E-04, 1.07416E-04, 6.20914E-05, 3.54202E-05, 2.51841E-05, 2.03993E-05, 1.59101E-05};
    
        Double_t widthMathieson04[n]  = {8.37835E-02, 4.21793E-02, 2.76214E-02, 1.92385E-02, 1.25567E-02, 9.18683E-03, 6.93955E-03, 5.47619E-03};
        Double_t ewidthMathieson04[n] = {1.10881E-03, 1.99743E-04, 1.01470E-04, 6.04768E-05, 3.06376E-05, 1.90336E-05, 1.27466E-05, 8.57512E-06};
    
    
    
    
    TCanvas *cwid = new TCanvas("cwid","width wrt charge",0,0,600,600);
    
    TGraphErrors *gr01 = new TGraphErrors(n, chginput, widthCOG, echginput, ewidthCOG);
    gr01->SetTitle("Residuals width wrt charge of clusters");
    gr01->SetMarkerColor(1);
    gr01->SetLineColor(1);
    gr01->SetMarkerStyle(8);
    gr01->GetXaxis()->SetTitle("charge of a cluster (ADC)");
    gr01->GetYaxis()->SetTitle("Residual width (cm)");
    gr01->GetHistogram()->SetMinimum(0.001);
    gr01->GetHistogram()->SetMaximum(0.1);
    gr01->Draw("AP");
    
        TGraphErrors *gr02 = new TGraphErrors(n, chginput, widthGauss, echginput, ewidthGauss);
         gr02->SetMarkerColor(2);
         gr02->SetLineColor(2);
         gr02->SetMarkerStyle(8);
        gr02->Draw("P SAME");
    
        TGraphErrors *gr03 = new TGraphErrors(n, chginput, widthDoubleGauss, echginput, ewidthDoubleGauss);
         gr03->SetMarkerColor(3);
         gr03->SetLineColor(3);
         gr03->SetMarkerStyle(8);
        gr03->Draw("P SAME");
    
    TGraphErrors *gr04 = new TGraphErrors(n, chginput, widthMathieson, echginput, ewidthMathieson);
    gr04->SetMarkerColor(4);
    gr04->SetLineColor(4);
    gr04->SetMarkerStyle(8);
    gr04->Draw("P SAME");
    
    
    
    auto legend = new TLegend(0.7,0.7,0.9,0.9);
                    legend->SetHeader("Clustering procedures"); // option "C" allows to center the header
                    legend->AddEntry(gr01, "Center of Gravity","lep");
                    legend->AddEntry(gr02, "Simple Gaussian fit, K3 AliRoot","lep");
                    legend->AddEntry(gr03, "Double Gaussian fit, K3 AliRoot","lep");
                    legend->AddEntry(gr04,"Mathieson fit, K3 AliRoot","lep");
                    legend->Draw();
    
    cwid->SetLogx(true);
    cwid->SetLogy(true);
    cwid->Update();
    cwid->Draw();
    
    // K3 related study
    
    TCanvas *cerrk = new TCanvas("cerrk","width wrt charge",0,0,600,600);
    
    TGraphErrors *gr010 = new TGraphErrors(n, chginput, widthMathieson01, echginput, ewidthMathieson01);
    gr010->SetTitle("Residuals width wrt charge of clusters for different K3");
    gr010->SetMarkerColor(1);
    gr010->SetLineColor(1);
    gr010->SetMarkerStyle(8);
    gr010->GetXaxis()->SetTitle("charge of a cluster (ADC)");
    gr010->GetYaxis()->SetTitle("Residual width (cm)");
    gr010->GetHistogram()->SetMinimum(0.001);
    gr010->GetHistogram()->SetMaximum(0.1);
    gr010->Draw("AP");
    
    TGraphErrors *gr020 = new TGraphErrors(n, chginput, widthMathieson02, echginput, ewidthMathieson02);
    gr020->SetMarkerColor(2);
    gr020->SetLineColor(2);
    gr020->SetMarkerStyle(8);
    gr020->Draw("P SAME");
    
    TGraphErrors *gr030 = new TGraphErrors(n, chginput, widthMathieson03, echginput, ewidthMathieson03);
    gr030->SetMarkerColor(3);
    gr030->SetLineColor(3);
    gr030->SetMarkerStyle(8);
    gr030->Draw("P SAME");
    
    TGraphErrors *gr040 = new TGraphErrors(n, chginput, widthMathieson04, echginput, ewidthMathieson04);
    gr040->SetMarkerColor(5);
    gr040->SetLineColor(5);
    gr040->SetMarkerStyle(8);
    gr040->Draw("P SAME");
    
    TGraphErrors *gr058 = new TGraphErrors(n, chginput, widthMathieson, echginput, ewidthMathieson);
    gr058->SetMarkerColor(4);
    gr058->SetLineColor(4);
    gr058->SetMarkerStyle(8);
    gr058->Draw("P SAME");
    
    
    
    auto legend2 = new TLegend(0.7,0.7,0.9,0.9);
                    legend2->SetHeader("Clustering procedures"); // option "C" allows to center the header
                    legend2->AddEntry(gr010, "Mathieson fit, K3 = 0,1","lep");
                    legend2->AddEntry(gr020, "Mathieson fit, K3 = 0,2","lep");
                    legend2->AddEntry(gr030, "Mathieson fit, K3 = 0,3","lep");
                    legend2->AddEntry(gr040, "Mathieson fit, K3 = 0,4","lep");
                    legend2->AddEntry(gr058,"Mathieson fit, K3 AliRoot = 0,584","lep");
                    legend2->Draw();
    
    cerrk->SetLogx(true);
    cerrk->SetLogy(true);
    cerrk->Update();
    cerrk->Draw();
    
    
    
}


void PowFitHitsWrtChg(){

    const Int_t n = 11;

    Double_t Kyinput[n]  = {0.1, 0.2, 0.3, 0.4, 0.5085, 0.5840, 0.6, 0.7, 0.8, 0.9, 1.0};
    Double_t eKyinput[n] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

    // RESULTATS Nhit

        Double_t A[n]  = {-4.21090E-01, -4.38204E-01, -4.37657E-01, -4.38090E-01, -4.50430E-01, -4.44181E-01, -4.37458E-01, -4.40858E-01, -4.61754E-01, -4.83430e-01, -5.14002E-01};
        Double_t eA[n] = {8.06906E-02, 8.07751E-02, 8.12839E-02, 8.16296E-02, 8.22543E-02, 8.28223E-02, 8.30776E-02, 8.38515E-02, 8.41772E-02, 8.46338e-02, 8.54575E-02};

        Double_t B[n]  = {1.82640E-01, 1.87841E-01, 1.91275E-01, 1.94537E-01, 1.98234E-01, 2.00095E-01, 2.00262E-01, 2.03291E-01, 2.06923E-01, 2.10539e-01, 2.14624E-01};
        Double_t eB[n] = {3.67935E-03, 3.56855E-03, 3.51521E-03, 3.46050E-03, 3.40752E-03, 3.39232E-03, 3.39789E-03, 3.36618E-03, 3.30577E-03, 3.25156e-03, 3.20159E-03};
    
    // RESULTATS NYhit
    
        Double_t Ay[n]  = {-2.04643E-01, -2.12065E-01, -1.95381E-01, -1.83955E-01, -1.88283E-01, -1.80734E-01, -1.76481E-01, -1.76784E-01, -1.89669E-01, -2.10244e-01, -2.42355E-01};
           Double_t eAy[n] = {7.43735E-02, 7.44994E-02, 7.52579E-02, 7.57952E-02, 7.64351E-02, 7.70176E-02, 7.73119E-02, 7.80585E-02, 7.85110E-02, 7.89559e-02, 7.96241E-02};

           Double_t By[n]  = {1.45580E-01, 1.51095E-01, 1.54191E-01, 1.57273E-01, 1.61149E-01, 1.63162E-01, 1.63482E-01, 1.66773E-01, 1.70706E-01, 1.74885e-01, 1.79758E-01};
           Double_t eBy[n] = {4.29345E-03, 4.15430E-03, 4.11402E-03, 4.06387E-03, 3.99849E-03, 3.97858E-03, 3.98475E-03, 3.94130E-03, 3.86845E-03, 3.79078e-03, 3.70815E-03};


    TCanvas *cpowA = new TCanvas("cpowA","Power Fit Coefficient A wrt K3",0,0,600,600);
    TCanvas *cpowB = new TCanvas("cpowB","Power Fit Coefficient B wrt K3",0,0,600,600);

    cpowA->cd();
    TGraphErrors *grA = new TGraphErrors(n, Kyinput, A, eKyinput, eA);
    grA->SetTitle("Power Fit Coefficient A wrt K3");
    grA->SetMarkerColor(1);
    grA->SetLineColor(1);
    grA->SetMarkerStyle(8);
    grA->GetXaxis()->SetTitle("K3y");
    grA->GetYaxis()->SetTitle("Values of Power fit coefficients");
    grA->Draw("AP");
    grA->GetYaxis()->SetRangeUser(-1,0.3);
    grA->Draw("AP");
    
    TF1 *fA = new TF1("f1","-6.20309E-01",0,1.2);
    fA->SetLineColor(1);
    fA->SetLineStyle(1);
    fA->Draw("SAME");
    
    TF1 *fAup = new TF1("f1","-6.20309E-01 + 2.64759E-01",0,1.2);
    fAup->SetLineColor(1);
    fAup->SetLineStyle(2);
    fAup->Draw("SAME");
    
    TF1 *fAdown = new TF1("f1","-6.20309E-01 - 2.64759E-01",0,1.2);
    fAdown->SetLineColor(1);
    fAdown->SetLineStyle(2);
    fAdown->Draw("SAME");

    
    TGraphErrors *grAy = new TGraphErrors(n, Kyinput, Ay, eKyinput, eAy);
    grAy->SetMarkerColor(2);
    grAy->SetLineColor(2);
    grAy->SetMarkerStyle(8);
    grAy->Draw("P SAME");
    
    TF1 *fAy = new TF1("f1","-3.45061E-01",0,1.2);
    fAy->SetLineColor(2);
    fAy->SetLineStyle(1);
    fAy->Draw("SAME");
    
    TF1 *fAyup = new TF1("f1","-3.45061E-01 + 2.38040E-01",0,1.2);
    fAyup->SetLineColor(2);
    fAyup->SetLineStyle(2);
    fAyup->Draw("SAME");
    
    TF1 *fAydown = new TF1("f1","-3.45061E-01 - 2.38040E-01",0,1.2);
    fAydown->SetLineColor(2);
    fAydown->SetLineStyle(2);
    fAydown->Draw("SAME");
    
    
    
    auto legendA = new TLegend(0.7,0.7,0.9,0.9);
                    legendA->SetHeader("Fit coefficients"); // option "C" allows to center the header
                    legendA->AddEntry(grA, "Constant A - Nhits","lep");
                    legendA->AddEntry(grAy, "Constant A - NYhits","lep");
                    legendA->Draw();
    
    cpowB->cd();
    TGraphErrors *grB = new TGraphErrors(n, Kyinput, B, eKyinput, eB);
     grB->SetTitle("Power Fit Coefficient B wrt K3");
     grB->SetMarkerColor(1);
     grB->SetLineColor(1);
     grB->SetMarkerStyle(5);
    grB->GetXaxis()->SetTitle("K3y");
    grB->GetYaxis()->SetTitle("Values of Power fit coefficients");
    grB->Draw("AP");
    grB->GetYaxis()->SetRangeUser(0.1,0.3);
    grB->Draw("AP");
    
    TF1 *fB = new TF1("f1","2.01511E-01",0,1.2);
       fB->SetLineColor(1);
       fB->SetLineStyle(1);
       fB->Draw("SAME");
       
       TF1 *fBup = new TF1("f1","2.01511E-01 + 1.07484E-02",0,1.2);
       fBup->SetLineColor(1);
       fBup->SetLineStyle(2);
       fBup->Draw("SAME");
       
       TF1 *fBdown = new TF1("f1","2.01511E-01 - 1.07484E-02",0,1.2);
       fBdown->SetLineColor(1);
       fBdown->SetLineStyle(2);
       fBdown->Draw("SAME");

        TGraphErrors *grBy = new TGraphErrors(n, Kyinput, By, eKyinput, eBy);
         grBy->SetMarkerColor(2);
         grBy->SetLineColor(2);
         grBy->SetMarkerStyle(5);
        grBy->Draw("P SAME");
    
    TF1 *fBy = new TF1("f1","1.63366E-01",0,1.2);
    fBy->SetLineColor(2);
    fBy->SetLineStyle(1);
    fBy->Draw("SAME");
    
    TF1 *fByup = new TF1("f1","1.63366E-01 + 1.22941E-02",0,1.2);
    fByup->SetLineColor(2);
    fByup->SetLineStyle(2);
    fByup->Draw("SAME");
    
    TF1 *fBydown = new TF1("f1","1.63366E-01 - 1.22941E-02",0,1.2);
    fBydown->SetLineColor(2);
    fBydown->SetLineStyle(2);
    fBydown->Draw("SAME");

    auto legendB = new TLegend(0.7,0.7,0.9,0.9);
                    legendB->SetHeader("Fit coefficients"); // option "C" allows to center the header
                    legendB->AddEntry(grB, "Constant B - Nhits","lep");
                    legendB->AddEntry(grBy, "Constant B - NYhits","lep");
                    legendB->Draw();

    cpowA->Update();
    cpowA->Draw();

    cpowB->Update();
    cpowB->Draw();


}


}
}
