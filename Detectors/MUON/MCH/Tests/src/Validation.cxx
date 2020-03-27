#include <fstream>
#include <vector>
#include <algorithm>
#include <functional>


#include "Validation.h"
#include "MCHMappingInterface/Segmentation.h"
#include "MCHMappingInterface/CathodeSegmentation.h"
#include "MCHMappingFactory/CreateSegmentation.h"
#include "MCHBase/Digit.h"
#include "MCHPreClustering/PreClusterBlock.h"
#include "MCHPreClustering/PreClusterFinder.h"
#include "DigitsFileReader.h"
#include "MCHClustering/ClusteringForTest.h"

#include "TMath.h"
#include "TRandom.h"
#include "TCanvas.h"
#include "TH2.h"
#include "TF2.h"
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
  std::string fname;
  preClusterFinder.init(fname);
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
    f1->SetParameters(0.5085, 0.5840, x, y);
    f1->SetParNames("K3x", "K3y", "Mean x", "Mean y");
}
    
//________________________________________________________________________________________
void myMath2hits(Double_t x1, Double_t y1, Double_t x2, Double_t y2, Double_t chg1, Double_t chg2){
    
    // Function generating a TWO HITS Mathieson histogram
    
    TF2 *f1 = new TF2("myMath",myMathieson2D2hits,-10,10,-10,10,8);
    f1->SetParameters(0.5085, 0.5840, x1, y1, x2, y2, chg1, chg2);
    f1->SetParNames("K3x", "K3y", "Mean1 x", "Mean1 y", "Mean2 x", "Mean2 y", "Charge1", "Charge2");
    //f1->Draw("colz");
}





//PART WITH: - PLOTTING OF THE MATHIESON FOLLOWING THE DEFINITION WE USE (ONE OR TWO HITS) AND MAPPING OF A GIVEN DETECTOR (DE819)
//           - DIGIT CREATION

//________________________________________________________________________________________
void Validation::PlotMathieson2D(Double_t x, Double_t y, int nsamples){
    
    digits.clear();
    
    TH2F* hb(NULL);
    TH2F* hnb(NULL);
    myMath1hit(x, y);
    int gPrintLevel = 0;
    
    o2::mch::mapping::CathodeSegmentation catsegb(819, kTRUE);
    o2::mch::mapping::CathodeSegmentation catsegnb(819, kFALSE);

    int nopadsb = catsegb.nofPads();
    int nopadsnb = catsegnb.nofPads();
    int padid;
    
    
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
    
    //Création et remplissage histogrammes bending
    
    cout << "Generating histograms bending and non-bending..."<< endl;

    hb = new TH2F("hb","hb",lowxsb.size()-1,xlowsb,lowysb.size()-1,ylowsb);
    gRandom->SetSeed(0);
    hb->FillRandom("myMath", /*nsamples=*/nsamples);
    
    TRandom *noisegen = new TRandom(321);
    for(int bin=0; bin<(lowxsb.size()+1)*(lowysb.size()+1); bin++){
        noise = noisegen->Gaus(1,0);
        hb->AddBinContent(bin, noise);
    }
    
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
        noise = noisegen->Gaus(1,0);
        hnb->AddBinContent(bin, noise);
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
    
    for(int i=0; i<nopadsb + nopadsnb; i++ ){
        charge[i] = 0.;
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
            binxcent = (hb->GetXaxis())->GetBinCenter(binx);
            binycent = (hb->GetYaxis())->GetBinCenter(biny);
            
//            if((binx*biny)%4 == 0){
//                cout << "Content of bending bin (binx= " << binx << ", biny= " << biny << "): " << bincontent << endl;
//                cout << "Center of bending bin (binx= " << binx << ", biny= " << biny << "): (" << binxcent << "," << binycent << ")" << endl;
//            }
            
            padid = catsegb.findPadByPosition(binxcent,binycent);
            charge[padid] += bincontent;
            
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
    
    for(int i=0; i<nopadsb + nopadsnb; i++){
        if(charge[i] > 6){  //Couper le bruit type (6ADC = 2*seuil de 3ADC, seuil de 3adc vient de 3*bruit/0.8)

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
    
    delete hb;
    delete hnb;
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

      char* preClustersBuffer = NULL;
      std::vector<mch::Clustering::Cluster> clusters(0);

        // load the digits from the memory buffer and run the pre-clustering phase
        preClusterFinder.reset();
        preClusterFinder.loadDigits(digitsBuffer, nDigits);
        preClusterFinder.run();

        // get number of pre-clusters and store them into a memory buffer
        auto preClustersSize = preClusterBlock.getPreClustersBufferSize(preClusterFinder);
        printf("preClustersSize: %d\n", (int)preClustersSize);
        preClustersBuffer = (char*)realloc(preClustersBuffer, preClustersSize);
        preClusterBlock.storePreClusters(preClusterFinder, preClustersBuffer);

        //continue;
        printf("\n\n==========\nReading clusters\n\n");

        std::vector<mch::PreClusterStruct> preClusters;
        preClusterBlock.readPreClusters(preClusters, preClustersBuffer, preClustersSize);
          
          printf("\n\n==========\nRunning Clustering\n\n");
    
    
    
    // Uncomment the desired clustering method
          
    
    //To run COG Clustering
//       clustering.runFinderCOG(preClusters, clusters);
//        printf("Number of clusters obtained and saved: %lu\n", clusters.size());
    
    
    //To run Mathieson fit Clustering
       clustering.runFinderSimpleFit(preClusters, clusters);
          
    
    //To run Gaussian fit Clustering
  //        clustering.runFinderGaussianFit(preClusters, clusters);
    
    
    //To run Double Gaussian fit Clustering
 //   clustering.runFinderDoubleGaussianFit(preClusters, clusters);
    
    delete digitsBuffer;
    delete preClustersBuffer;
    
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
    Double_t widthMathieson[n]  = {0.0724216, 0.0473039, 0.0237125, 0.0195968, 0.0136211, 0.00922341, 0.00550195, 0.00424860};
    Double_t ewidthMathieson[n] = {0.00823728, 0.00253988, 0.000958042, 0.000588419, 0.000349494, 0.000168508, 0.0000935699, 0.0000546348};
    Double_t widthDoubleGauss[n]  = {0.0809867, 0.0405142, 0.0245105, 0.0198315, 0.0130453, 0.00867564, 0.00585070, 0.00422690};
    Double_t ewidthDoubleGauss[n] = {0.0110478, 0.00196554, 0.000973368, 0.000605342, 0.000327329, 0.000196655, 0.0000988013, 0.0000564364};
    Double_t widthGauss[n]  = {0.115848, 0.0386573, 0.0268345, 0.0240388, 0.0184959, 0.0164656, 0.0157473, 0.0171971};
    Double_t ewidthGauss[n] = {0.0261366, 0.00189652, 0.00106311, 0.000851333, 0.000518058, 0.000461726, 0.000396683, 0.000548232};
    Double_t widthCOG[n]  = {0.111975, 0.0470066, 0.0272162, 0.0258470, 0.0250620, 0.0212521, 0.0182168, 0.0266299};
    Double_t ewidthCOG[n] = {0.0234304, 0.00250890, 0.00103235, 0.000912423, 0.000890655, 0.000643082, 0.000500731, 0.00143591};
    
    TCanvas *cerr = new TCanvas("cerr","width wrt charge",0,0,600,600);
    
    TGraphErrors *gr2 = new TGraphErrors(n, chginput, widthMathieson, echginput, ewidthMathieson);
    gr2->SetTitle("Residuals width wrt charge of clusters");
    gr2->SetMarkerColor(4);
    gr2->SetLineColor(4);
    gr2->SetMarkerStyle(8);
    gr2->GetXaxis()->SetTitle("charge of a cluster (ADC)");
    gr2->GetYaxis()->SetTitle("Residual width (cm)");
    gr2->Draw("AP");
    
    TGraphErrors *gr3 = new TGraphErrors(n, chginput, widthDoubleGauss, echginput, ewidthDoubleGauss);
     gr3->SetMarkerColor(2);
     gr3->SetLineColor(2);
     gr3->SetMarkerStyle(8);
    gr3->Draw("P SAME");
    
    TGraphErrors *gr4 = new TGraphErrors(n, chginput, widthCOG, echginput, ewidthCOG);
     gr4->SetMarkerColor(3);
     gr4->SetLineColor(3);
     gr4->SetMarkerStyle(8);
    gr4->Draw("P SAME");
    
    TGraphErrors *gr5 = new TGraphErrors(n, chginput, widthGauss, echginput, ewidthGauss);
     gr5->SetMarkerColor(1);
     gr5->SetLineColor(1);
     gr5->SetMarkerStyle(8);
    gr5->Draw("P SAME");
    
    auto legend = new TLegend(0.7,0.7,0.9,0.9);
                    legend->SetHeader("Clustering procedures"); // option "C" allows to center the header
                    legend->AddEntry(gr2,"Mathieson fit, K3 AliRoot","lep");
                    legend->AddEntry(gr3, "Double Gaussian fit, K3 AliRoot","lep");
                    legend->AddEntry(gr5, "Simple Gaussian fit, K3 AliRoot","lep");
                    legend->AddEntry(gr4, "Center of Gravity","lep");
                    legend->Draw();
    
    cerr->SetLogx(true);
    cerr->SetLogy(true);
    cerr->Update();
    cerr->Draw();
    
    
}


}
}
