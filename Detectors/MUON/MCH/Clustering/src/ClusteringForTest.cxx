//
//  ClusteringForTest.cxx
//
//  Algo de clustring repris de AliRoot, ClusterCOG
//  This program implements the clustering part
//  There are currently 4 methods implemented:
//  CenterOfGravity, (Mathieson)SimpleFit, Gaussian Fit and Double Gaussian Fit
//
//  Each fit exixts in 2 variations. An old hand-wavey one where the chi2 minimisation is done manually
//  A better one using MIGRAD
//
// FixMe: Find a better way (using TObject) to send digits data to the chi2 minimisation function instead of defining parameters by hand and keeping them fixed
//
//  Created by PERRIN Sébastien on 31/10/2019.
//

#include <chrono>
#include <memory>
#include <stdexcept>
#include <iostream>
#include <vector>
#include <TMath.h>
#include <stdio.h>
#include <math.h>

#include "MCHClustering/ClusteringForTest.h"
#include "MCHBase/Digit.h"
#include "MCHMappingInterface/Segmentation.h"
#include "MCHPreClustering/PreClusterFinder.h"
#include "MCHPreClustering/PreClusterBlock.h"
#include "TVirtualFitter.h"
#include "TObject.h"
#include "TObjArray.h"
#include "TVectorT.h"
#include "MCHClustering/MyObjectDigits.h"
#include "MCHClustering/MyObjectParametersMathieson.h"
#include "MCHClustering/MyObjectParametersGauss.h"
#include "MCHClustering/MyObjectParametersDoubleGauss.h"

#include "TRandom2.h"
#include "TError.h"


//#include "TVector2.h"

#include <fairmq/Tools.h>
#include <FairMQLogger.h>

namespace o2
{
namespace mch
{

// templateClassImp(TVectorT);

using namespace std;


// FONCTIONS POUR CLUSTERING COG


//_________________________________________________________________________________________________
void Clustering::runFinderCOG(std::vector<PreClusterStruct>& preClusters, std::vector<Cluster>& clusters)
{
    cout << "preClusters.size():" << preClusters.size() << endl;
    int sizedigit;
    int jete = 0;
    Digit digittmp;
    Cluster clustertmp;
    int gPrintLevel = 2;
    
    for (Int_t i = 0; i < preClusters.size(); i++)   //On boucle sur les preclusters identifiés
    {
        if(preClusters[i].nDigits == 1){         // Je vire les preclusters qui n'ont qu'un seul digit
            jete++;
            continue;
        }
        cout << "Precluster avec " << preClusters[i].nDigits << " digits." << endl;  //Nombre de digits dans le precluster aue l'on regarde
        
        std::vector<Digit> preclustertmp;
        const Digit* ptrdigit = preClusters[i].digits;
        
        for (int j=0; j<preClusters[i].nDigits; j++){          //Je récupère les digits et les mets dans un   vecteur de digits, preclustertmp
            digittmp = *ptrdigit;
            if(gPrintLevel > 1){
                cout << "Digit numero:" << j << endl;
                cout << "digittmp.getADC():" << digittmp.getADC() << endl;
                cout << "digittmp.getDetID():" << digittmp.getDetID() << endl;
                cout << "digittmp.getPadID():" << digittmp.getPadID() << endl;
    //            cout << "type digittmp:" << typeid(digittmp).name() << endl;
    //            cout << "type ptrdigit:" << typeid(ptrdigit).name() << endl;
    //            cout << "ptrdigit:" << ptrdigit << endl;
            }
            preclustertmp.push_back(digittmp);
            sizedigit = sizeof(digittmp);
            ptrdigit++;
        }
        
        //Je clusterise mon vecteur de digits, preclustertmp et l'ajoute au vecteur de clusters, nommé clusters.
        
       // cout << "preclustertmp[0].getADC():" << preclustertmp[0].getADC() << endl;
        clustertmp = FinderCOG(preclustertmp);
      //  cout << "COORD X DU CLUSTER AJOUTÉ A CLUSTERS: " << clustertmp.getx() << endl;
        clusters.push_back(clustertmp);
    }
    
    cout << "Nombre de preclusters au départ:" << preClusters.size() << endl;
    cout << "Nombre de preclusters de taille 1, abandonnés:" << jete << endl;
    
}

//_________________________________________________________________________________________________
Clustering::Cluster Clustering::FinderCOG(std::vector<Digit> &precluster)
{
    cout << "\n\n==========\nRunning COG Algorithm\n\n" << endl;

    Double_t xmin = 1E9;
    Double_t ymin = 1E9;
    Double_t xmax = -1E9;
    Double_t ymax = -1E9;
    Double_t charge[2] = { 0.0, 0.0 };
    Double_t multiplicity[2] = { 0.0, 0.0 };
    Cluster cluster;
    int gPrintLevel = 2;
    
    Double_t x[] = { 0.0, 0.0 };
    Double_t y[] = { 0.0, 0.0 };
    
    Double_t xsize[] = { 0.0, 0.0 };
    Double_t ysize[] = { 0.0, 0.0 };
    
         
    for ( Int_t i = 0; i < precluster.size(); ++i ) //On boucle sur les digits de notre precluster
    {
       int detid = precluster[i].getDetID();
       int padid = precluster[i].getPadID();
        
        if(gPrintLevel > 1){
            cout << "\nDetID:" << detid << " PadID:" << padid << endl;
        }
       mapping::Segmentation pad(detid); // = mapping::Segmentation(detid);
       
       
       
       //On définit le vecteur position et taille du pad en question
     std::vector<Double_t> padPosition = {pad.padPositionX(padid), pad.padPositionY(padid)};
     std::vector<Double_t> padSize = {pad.padSizeX(padid), pad.padSizeY(padid)};
       
        if(gPrintLevel > 1){
           cout << "PadPosition: "<< padPosition[0] << " " << padPosition[1] << endl;
           cout << "PadSize: "<< padSize[0] << " " << padSize[1] << endl;
        }
        
       //On met à jour les xminmax et yminmax
     xmin = TMath::Min(padPosition[0]-0.5*padSize[0],xmin);
     xmax = TMath::Max(padPosition[0]+0.5*padSize[0],xmax);
     ymin = TMath::Min(padPosition[1]-0.5*padSize[1],ymin);
     ymax = TMath::Max(padPosition[1]+0.5*padSize[1],ymax);
        
        for ( Int_t cathode = 0; cathode < 2; ++cathode ) //On boucle sur les deux plans de cathodes
         {
           
             if ( cathode == pad.isBendingPad(padid) ) //On regarde à quelle cathode appartient le pad
             {
                 if(gPrintLevel > 1){
                     cout << "Cathode:" << cathode << endl;
                 }
                 
                 //On ajoute un terme au numerateur du COG et de la précision
               x[cathode] += padPosition[0]*precluster[i].getADC();
               y[cathode] += padPosition[1]*precluster[i].getADC();
               xsize[cathode] += padSize[0];
               ysize[cathode] += padSize[1];
                 
                 //On ajoute un terme au dénominateur du COG et de la précision
               charge[cathode] += precluster[i].getADC();
               multiplicity[cathode] += 1;
             }
             
         }
                
    }
           
    for ( Int_t cathode = 0; cathode < 2; ++cathode ) //On boucle sur les deux plans de cathodes
    {
//        cout << "\nCathode:" << cathode << endl;
//        cout << "x:" << x[cathode] << endl;
//        cout << "y:" << y[cathode] << endl;
//        cout << "Charge:" << charge[cathode] << endl;
//        cout << "xsize:" << xsize[cathode] << endl;
//        cout << "ysize:" << ysize[cathode] << endl;
//        cout << "Multiplicity:" << multiplicity[cathode] << endl;
        
                 //On fait les calculs de COG et précision
               if ( charge[cathode] != 0 )
               {
                 x[cathode] /= charge[cathode];
                 y[cathode] /= charge[cathode];
               }
               if ( multiplicity[cathode] != 0 )
               {
                 xsize[cathode] /= (multiplicity[cathode]*sqrt(charge[cathode]));
                 ysize[cathode] /= (multiplicity[cathode]*sqrt(charge[cathode]));
               }
                else if( multiplicity[cathode] == 0)
                {
                    xsize[cathode] = 1E9;
                    ysize[cathode] = 1E9;
                }
        
        if(gPrintLevel > 1){
            cout << "\nPost-calcul cathode:" << cathode << endl;
            cout << "x:" << x[cathode] << endl;
            cout << "y:" << y[cathode] << endl;
            cout << "Charge:" << charge[cathode] << endl;
            cout << "xsize:" << xsize[cathode] << endl;
            cout << "ysize:" << ysize[cathode] << endl;
            cout << "Multiplicity:" << multiplicity[cathode] << endl;
        }
        
    }
    

  Double_t xCOG = 0.0;
  Double_t yCOG = 0.0;
  Double_t ex = 0.0;
  Double_t ey = 0.0;

  // On regarde sur quelle cathode la précision est la meilleure pour définir x et y du COG.
  xCOG = ( xsize[0] < xsize[1] ) ? x[0] : x[1];
  yCOG = ( ysize[0] < ysize[1] ) ? y[0] : y[1];
    
    ex = ( xsize[0] < xsize[1] ) ? xsize[0] : xsize[1];
    ey = ( ysize[0] < ysize[1] ) ? ysize[0] : ysize[1];
    
//    cout << "ex:" << ex << endl;
//    cout << "ey:" << ey << endl;
    
    double timestamp = precluster[0].getTimeStamp();  // On récupère le timestamp du premier digit du precluster
//
//    // On print en console les infos importantes pour pouvoir checker les problèmes
//    printf("\n\nCluster multiplicity %ld (x,y)=(%e,%e) boundaries=(xmin,ymin,xmax,ymax)=(%e,%e,%e,%e)"
//                  " (x0,y0,x1,y1)=(%e,%e,%e,%e) timestamp = %lf \n\n",
//                  precluster.size(),xCOG,yCOG,xmin,ymin,xmax,ymax,
//                  x[0],y[0],x[1],y[1],timestamp);
//
    //On remplit un cluster avec les infos nécessaires
    cluster.setx(xCOG);
    cluster.sety(yCOG);
    cluster.setex(ex);
    cluster.setey(ey);
    cluster.settimestamp(timestamp);    // En première approx. on peut affecter le timestamp du pre;ier digit au cluster
    
    printf("\n\nCluster multiplicity %ld \n(x,y)=(%e,%e) \nboundaries=(xmin,ymin,xmax,ymax)=(%e,%e,%e,%e) \nprecision=(ex,ey)=(%e,%e) \ntimestamp = %lf \n\n",precluster.size(),cluster.getx(),cluster.gety(),xmin,ymin,xmax,ymax,cluster.getex(),cluster.getey(),cluster.gettimestamp());
    
    return cluster;
    
    
//     PRINCIPE GENERAL COG
//     Pour chaque cathode
//        Pour chaque pad dans le precluster
//            Mettre à jour les valeurs des limites du cluster
//            Ajouter les termes au numerateur de x,y COG et au calcul de la precision
//        Diviser par les denominateurs
//     Obtenir les parametres du cluster
     
    
}

//_________________________________________________________________________________________________





// FONCTIONS POUR CLUSTERING MATHIESON (SIMPLE) FIT




//_________________________________________________________________________________________________
void Clustering::runFinderSimpleFit(std::vector<PreClusterStruct>& preClusters, std::vector<Cluster>& clusters)
{
    cout << "\n\n==========\nRunning runFinderSimpleFit Algorithm\n\n" << endl;
    cout << "preClusters.size():" << preClusters.size() << endl;
    int sizedigit;
    int jete = 0;
    Digit digittmp;
    Cluster clustertmpCOG;
    Cluster clustertmp;
    int gPrintLevel = 2;
    
    for (Int_t i = 0; i < preClusters.size(); i++)   //On boucle sur les preclusters identifiés
    {
        if(preClusters[i].nDigits == 1){         // Je vire les preclusters qui n'ont qu'un seul digit
            jete++;
            continue;
        }
        cout << "Precluster avec " << preClusters[i].nDigits << " digits." << endl;  //Nombre de digits dans le precluster aue l'on regarde
        
        std::vector<Digit> preclustertmp;
        const Digit* ptrdigit = preClusters[i].digits;
        
        for (int j=0; j<preClusters[i].nDigits; j++){          //Je récupère les digits et les mets dans un   vecteur de digits, preclustertmp
            digittmp = *ptrdigit;
            if(gPrintLevel > 1){
                cout << "Digit numero:" << j << endl;
                cout << "digittmp.getADC():" << digittmp.getADC() << endl;
                cout << "digittmp.getDetID():" << digittmp.getDetID() << endl;
                cout << "digittmp.getPadID():" << digittmp.getPadID() << endl;
    //            cout << "type digittmp:" << typeid(digittmp).name() << endl;
    //            cout << "type ptrdigit:" << typeid(ptrdigit).name() << endl;
    //            cout << "ptrdigit:" << ptrdigit << endl;
            }
            preclustertmp.push_back(digittmp);
            sizedigit = sizeof(digittmp);
            ptrdigit++;
        }
        
        //Je clusterise mon vecteur de digits, preclustertmp et l'ajoute au vecteur de clusters, nommé clusters.
        
       // cout << "preclustertmp[0].getADC():" << preclustertmp[0].getADC() << endl;
        clustertmpCOG = FinderCOG(preclustertmp);
      //  cout << "COORD X DU CLUSTER AJOUTÉ A CLUSTERS: " << clustertmp.getx() << endl;
        if(preclustertmp.size() < 3){
            clustertmp = clustertmpCOG;
        }
        else{
    //    clustertmp = FinderSimpleFit(preclustertmp, clustertmpCOG);    //Fit manuel
        clustertmp = ComputePositionClean(preclustertmp, clustertmpCOG);      //Fit Minuit
        }
        
        clusters.push_back(clustertmp);
    }
    
    cout << "Nombre de preclusters au départ:" << preClusters.size() << endl;
    cout << "Nombre de preclusters de taille 1, abandonnés:" << jete << endl;
    
}

//_________________________________________________________________________________________________
Clustering::Cluster Clustering::ComputePositionClean(std::vector<Digit> &precluster, Clustering::Cluster clustertmpCOG)
{
  /// Compute the position of the given cluster, by fitting a Mathieson
  /// charge distribution to it

    Double_t vecxmin[2] = {1E3, 1E3};
    Double_t vecxmax[2] = {-1E3, -1E3};
    Double_t vecymin[2] = {1E3, 1E3};
    Double_t vecymax[2] = {-1E3, -1E3};
    
  TVirtualFitter* fitter = TVirtualFitter::Fitter(0,2);
  fitter->SetFCN(FitFunctionClean);

  cout << "\n\n==========\nRunning MinuitMathieson Algorithm.\n\n" << endl;

  Cluster cluster;

  //Création des vecteurs position du hit sur chaque cqthode
  Double_t xhit[2] = { clustertmpCOG.getx(), clustertmpCOG.getx() };
  Double_t yhit[2] = { clustertmpCOG.gety(), clustertmpCOG.gety() };

    //Création des vecteurs de charge totale, multiplicité et tailles sur chaque cathode
       Double_t chargetot[2] = { 0.0, 0.0 };
       Double_t multiplicity[2] = { 0.0, 0.0 };
       Double_t xsize[2] = { 0.0, 0.0 };
       Double_t ysize[2] = { 0.0, 0.0 };

   //Détermination des vecteurs de charge totale, multiplicité et taille qui devient précision
   for ( Int_t i = 0; i < precluster.size(); ++i ){
       int detid = precluster[i].getDetID();
       int padid = precluster[i].getPadID();
       
       mapping::Segmentation pad(detid); // = mapping::Segmentation(detid);

       for ( int cathode = 0; cathode < 2; ++cathode ) //On boucle sur les deux plans de cathodes
       {

           if ( cathode == pad.isBendingPad(padid) ) //On regarde à quelle cathode appartient le pad
                   {
                       std::vector<Double_t> padPosition = {pad.padPositionX(padid), pad.padPositionY(padid)};
                       std::vector<Double_t> padSize = {pad.padSizeX(padid), pad.padSizeY(padid)};
                       
                       vecxmin[cathode] = TMath::Min(padPosition[0]-0.5*padSize[0],vecxmin[cathode]);
                       vecxmax[cathode] = TMath::Max(padPosition[0]+0.5*padSize[0],vecxmax[cathode]);
                       vecymin[cathode] = TMath::Min(padPosition[1]-0.5*padSize[1],vecymin[cathode]);
                       vecymax[cathode] = TMath::Max(padPosition[1]+0.5*padSize[1],vecymax[cathode]);

                       chargetot[cathode] += precluster[i].getADC();
                       multiplicity[cathode] += 1;
                   }
       }
   }
    

  Float_t stepX = 0.00001; // cm
  Float_t stepY = 0.00001; // cm
    Double_t Kx3 = 0.5085;
    Double_t Ky3 = 0.5840;
  //  Double_t Ky3 = 0.2; //Compare with Alberto
    
    MyObjectDigits digits;

//    Digit* digitptr = &precluster[0];
//    TVectorT<Digit> digits = TVectorT<Digit>(precluster.size(), digitptr);
    
    cout << "On va set digits" << endl;
    digits.setfDigits(precluster);
    cout << (digits.getfDigits())[0].getADC() << " est ADC du premier digit" << endl;
  //  Digit loule = digits.operator()(0);
   // cout << "ADC du digit 0: " << loule.getADC() << endl;
    
    Double_t arg(-2); // disable printout

    for ( int cathode = 0; cathode < 2; ++cathode ) //On boucle sur les deux plans de cathodes
    {
        fitter->Clear();
        fitter->ExecuteCommand("SET PRINT",&arg,1);
        if(vecxmin[cathode]>vecxmax[cathode]){
            vecxmin[cathode] = -40.;
            vecxmax[cathode] = +40.;
        }
        if(vecymin[cathode]>vecymax[cathode]){
            vecymin[cathode] = -20.;
            vecymax[cathode] = +20.;
        }
        cout << "On va set les params" << endl;
        
        MyObjectParametersMathieson parameters;
        parameters.setfCathode(cathode);
        parameters.setfChargetot(chargetot[cathode]);
        parameters.setfKx3(Kx3);
        parameters.setfKy3(Ky3);
        parameters.setfNbdigits(precluster.size());
        
//        double params[5] = {static_cast<double>(cathode), chargetot[cathode], Kx3, Ky3, static_cast<double>(precluster.size())};
//        double* paramptr = &params[0];
//        TVectorT<double> parameters = TVectorT<double>(5, paramptr);
        
        
        TObjArray userObjects;
        userObjects.Add(&digits);
        userObjects.Add(&parameters);
        fitter->SetObjectFit(&userObjects);
        
        fitter->SetParameter(0,"cluster X position",xhit[cathode],stepX,vecxmin[cathode],vecxmax[cathode]);
        fitter->SetParameter(1,"cluster Y position",yhit[cathode],stepY,vecymin[cathode],vecymax[cathode]);
        

        Double_t stratArg(2);
        fitter->ExecuteCommand("SET STR",&stratArg,1);
        Int_t val = fitter->ExecuteCommand("MIGRAD",0,0);
        if ( val && chargetot[cathode] != 0 )
        {
            //Fit failed with robust strategy 2 try balanced strategy 1
            Double_t stratArg(1);
            fitter->ExecuteCommand("SET STR",&stratArg,1);
            Int_t val2 = fitter->ExecuteCommand("MIGRAD",0,0);
            if ( val2 ){
                // fit failed. Using COG results, with big errors
                cout << "Fit failed on viable cathode. Using COG results for cluster." << endl;
                    return clustertmpCOG;
            }
        }

        xhit[cathode] = fitter->GetParameter(0);
        yhit[cathode] = fitter->GetParameter(1);
        xsize[cathode] = fitter->GetParError(0);
        ysize[cathode] = fitter->GetParError(1);

        Double_t amin, edm, errdef;
        Int_t nvpar, nparx;

        fitter->GetStats(amin, edm, errdef, nvpar, nparx);

        Double_t chi2 = amin;

        printf("Cluster fitted to (x,y)=(%e,%e) (xerr,yerr)=(%e,%e) \n chi2=%e ndf=%d",
                        xhit[cathode],yhit[cathode],
                        xsize[cathode],ysize[cathode],chi2,fitter->GetNumberFreeParameters());
    }

    //Détermination de la position finale du hit à l'aide de la résolution sur chaque cathode
    Double_t xhitfinal = xhit[0];
    Double_t yhitfinal = yhit[1];
    Double_t ex = xsize[0];
    Double_t ey = ysize[1];

    double timestamp = precluster[0].getTimeStamp();  // On récupère le timestamp du premier digit du precluster


    //On remplit un cluster avec les infos nécessaires
    cluster.setx(xhitfinal);
    cluster.sety(yhitfinal);
    cluster.setex(ex);
    cluster.setey(ey);
    cluster.settimestamp(timestamp);    // En première approx. on peut affecter le timestamp du pre;ier digit au cluster

            printf("\n\nCluster par simple fit:\nCluster multiplicity %ld \n(x,y)=(%e,%e) \nprecision=(ex,ey)=(%e,%e) \ntimestamp = %lf \n\n",precluster.size(),cluster.getx(),cluster.gety(),cluster.getex(),cluster.getey(),cluster.gettimestamp());

    return cluster;
}

//_________________________________________________________________________________________________
void FitFunctionClean(Int_t& /*notused*/, Double_t* /*notused*/,
                            Double_t& sum, Double_t* par,
                             Int_t /*notused*/)
  {
    /// Chi2 Function to minimize: Mathieson charge distribution in 2 dimensions
      
          //Méthode propre. On amène un pointeur vers les digits aue l'on peut lire un par un. Et un TArray avec les autres parametres (cathode, chgtot, Kx3, Ky3, nbdigits)

    TObjArray* userObjects = static_cast<TObjArray*>(TVirtualFitter::GetFitter()->GetObjectFit());
      
      MyObjectDigits* digits = static_cast<MyObjectDigits*>(userObjects->At(0));
      MyObjectParametersMathieson* parameters = static_cast<MyObjectParametersMathieson*>(userObjects->At(1));
      
//    TVectorT<Digit>* digits = static_cast<TVectorT<Digit>*>(userObjects->At(0));
//    TVectorT<double>* parameters = static_cast<TVectorT<double>*>(userObjects->At(1));
      
      Clustering cluster;

    sum = 0.0;
      
      
      int cathode = parameters->getfCathode();
      double chargetot = parameters->getfChargetot();
      double Kx3 = parameters->getfKx3();
      double Ky3 = parameters->getfKy3();
      int nbdigits = parameters->getfNbdigits();
      
//      int cathode = parameters->operator()(0);
//      double chargetot = parameters->operator()(1);
//      double Kx3 = parameters->operator()(2);
//      double Ky3 = parameters->operator()(3);
//      int nbdigits = parameters->operator()(4);
      
    for(Int_t i = 0; i < nbdigits; ++i){

//        Digit digitbeingread = *digitsptr;
//        digitsptr++;
      
        Digit digitbeingread = (digits->getfDigits())[i];
        
    //   Digit digitbeingread = digits->operator()(i);
        
        int detid = digitbeingread.getDetID();
        int padid = digitbeingread.getPadID();
        int ADC = digitbeingread.getADC();
        
        mapping::Segmentation pad(detid); // = mapping::Segmentation(detid);

        if ( cathode == pad.isBendingPad(padid) ) //On regarde à quelle cathode appartient le pad
        {
              //On définit le vecteur position et taille du pad en question
            std::vector<Double_t> padPosition = {pad.padPositionX(padid), pad.padPositionY(padid)};
            std::vector<Double_t> padSize = {pad.padSizeX(padid), pad.padSizeY(padid)};

            std::vector<Double_t> lowerLeft = {par[0] - padPosition[0] - 0.5*padSize[0], par[1] - padPosition[1] - 0.5*padSize[1]};
            std::vector<Double_t> upperRight= {lowerLeft[0] + padSize[0], lowerLeft[1] + padSize[1]};

            Double_t qfit = chargetot*cluster.IntMathiesonXY(lowerLeft[0], lowerLeft[1], upperRight[0], upperRight[1], Kx3, Ky3);

            Double_t qtrue = ADC;

            Double_t delta = qtrue-qfit;

            sum += pow(delta,2)/qtrue;
        }
    }
}

//_________________________________________________________________________________________________
Double_t Clustering::IntMathiesonXY(Double_t x1, Double_t y1, Double_t x2, Double_t y2, Double_t Kx3, Double_t Ky3)
{
    Double_t pitch = 0.25;                                          //Adapted fron AliMuonMathieson
    Double_t Kx2 = M_PI / 2. * (1. - 0.5 * sqrt(Kx3));
    Double_t cx1 = Kx2 * sqrt(Kx3) / 4. / atan(sqrt(Kx3));
    Double_t Kx4 = cx1 / Kx2 / sqrt(Kx3);
    
    Double_t Ky2 = M_PI / 2. * (1. - 0.5 * sqrt(Ky3));
    Double_t cy1 = Ky2 * sqrt(Ky3) / 4. / atan(sqrt(Ky3));
    Double_t Ky4 = cy1 / Ky2 / sqrt(Ky3);
    
    x1 /= pitch;
    x2 /= pitch;
    y1 /= pitch;
    y2 /= pitch;

    Double_t ux1 = sqrt(Kx3)*tanh(Kx2*x1);
    Double_t ux2 = sqrt(Kx3)*tanh(Kx2*x2);
    
    Double_t uy1 = sqrt(Ky3)*tanh(Ky2*y1);
    Double_t uy2 = sqrt(Ky3)*tanh(Ky2*y2);
    
    return Double_t(4.*Kx4*(atan(ux2)-atan(ux1))*
                   Ky4*(atan(uy2)-atan(uy1)));
}

//_________________________________________________________________________________________________







// FONCTIONS POUR CLUSTERING GAUSS FIT





//_________________________________________________________________________________________________
void Clustering::runFinderGaussianFit(std::vector<PreClusterStruct>& preClusters, std::vector<Cluster>& clusters)
{
    cout << "preClusters.size():" << preClusters.size() << endl;
    int sizedigit;
    int jete = 0;
    Digit digittmp;
    Cluster clustertmpCOG;
    Cluster clustertmp;
    
    for (Int_t i = 0; i < preClusters.size(); i++)   //On boucle sur les preclusters identifiés
    {
        if(preClusters[i].nDigits == 1){         // Je vire les preclusters qui n'ont qu'un seul digit
            jete++;
            continue;
        }
        cout << "Precluster avec " << preClusters[i].nDigits << " digits." << endl;  //Nombre de digits dans le precluster aue l'on regarde
        
        std::vector<Digit> preclustertmp;
        const Digit* ptrdigit = preClusters[i].digits;
        
        for (int j=0; j<preClusters[i].nDigits; j++){          //Je récupère les digits et les mets dans un   vecteur de digits, preclustertmp
            digittmp = *ptrdigit;
            cout << "Digit numero:" << j << endl;
            cout << "digittmp.getADC():" << digittmp.getADC() << endl;
            cout << "digittmp.getDetID():" << digittmp.getDetID() << endl;
            cout << "digittmp.getPadID():" << digittmp.getPadID() << endl;
//            cout << "type digittmp:" << typeid(digittmp).name() << endl;
//            cout << "type ptrdigit:" << typeid(ptrdigit).name() << endl;
//            cout << "ptrdigit:" << ptrdigit << endl;
            preclustertmp.push_back(digittmp);
            sizedigit = sizeof(digittmp);
            ptrdigit++;
        }
        
        //Je clusterise mon vecteur de digits, preclustertmp et l'ajoute au vecteur de clusters, nommé clusters.
        
       // cout << "preclustertmp[0].getADC():" << preclustertmp[0].getADC() << endl;
        clustertmpCOG = FinderCOG(preclustertmp);
      //  cout << "COORD X DU CLUSTER AJOUTÉ A CLUSTERS: " << clustertmp.getx() << endl;
        if(preclustertmp.size() < 3){
            clustertmp = clustertmpCOG;
        }
        else{
   //     clustertmp = FinderGaussianFit(preclustertmp, clustertmpCOG);         //Fit manuel
        clustertmp = ComputePositionGaussianFitClean(preclustertmp, clustertmpCOG);  //Fit Minuit
        }
        
        clusters.push_back(clustertmp);
    }
    
    cout << "Nombre de preclusters au départ:" << preClusters.size() << endl;
    cout << "Nombre de preclusters de taille 1, abandonnés:" << jete << endl;
    
}

//_________________________________________________________________________________________________
Clustering::Cluster Clustering::ComputePositionGaussianFitClean(std::vector<Digit> &precluster, Clustering::Cluster clustertmpCOG)
{
  /// Compute the position of the given cluster, by fitting a Gaussian
  /// charge distribution to it

        Double_t vecxmin[2] = {1E3, 1E3};
        Double_t vecxmax[2] = {-1E3, -1E3};
        Double_t vecymin[2] = {1E3, 1E3};
        Double_t vecymax[2] = {-1E3, -1E3};
        
      TVirtualFitter* fitter = TVirtualFitter::Fitter(0,2);
      fitter->SetFCN(FitFunctionGaussianFitClean);

      cout << "\n\n==========\nRunning MinuitGaussian Algorithm.\n\n" << endl;

      Cluster cluster;

      //Création des vecteurs position du hit sur chaque cqthode
      Double_t xhit[2] = { clustertmpCOG.getx(), clustertmpCOG.getx() };
      Double_t yhit[2] = { clustertmpCOG.gety(), clustertmpCOG.gety() };

        //Création des vecteurs de charge totale, multiplicité et tailles sur chaque cathode
           Double_t chargetot[2] = { 0.0, 0.0 };
           Double_t multiplicity[2] = { 0.0, 0.0 };
           Double_t xsize[2] = { 0.0, 0.0 };
           Double_t ysize[2] = { 0.0, 0.0 };

       //Détermination des vecteurs de charge totale, multiplicité et taille qui devient précision
       for ( Int_t i = 0; i < precluster.size(); ++i ){
           int detid = precluster[i].getDetID();
           int padid = precluster[i].getPadID();
           
           mapping::Segmentation pad(detid); // = mapping::Segmentation(detid);

           for ( int cathode = 0; cathode < 2; ++cathode ) //On boucle sur les deux plans de cathodes
           {

               if ( cathode == pad.isBendingPad(padid) ) //On regarde à quelle cathode appartient le pad
                       {
                           std::vector<Double_t> padPosition = {pad.padPositionX(padid), pad.padPositionY(padid)};
                           std::vector<Double_t> padSize = {pad.padSizeX(padid), pad.padSizeY(padid)};
                           
                           vecxmin[cathode] = TMath::Min(padPosition[0]-0.5*padSize[0],vecxmin[cathode]);
                           vecxmax[cathode] = TMath::Max(padPosition[0]+0.5*padSize[0],vecxmax[cathode]);
                           vecymin[cathode] = TMath::Min(padPosition[1]-0.5*padSize[1],vecymin[cathode]);
                           vecymax[cathode] = TMath::Max(padPosition[1]+0.5*padSize[1],vecymax[cathode]);

                           chargetot[cathode] += precluster[i].getADC();
                           multiplicity[cathode] += 1;
                       }
           }
       }
        

      Float_t stepX = 0.00001; // cm
      Float_t stepY = 0.00001; // cm
        Double_t sigx = 0.7604;
        Double_t sigy = 0.7808;
        
        
        MyObjectDigits digits;
        
        cout << "On va set digits" << endl;
        digits.setfDigits(precluster);
        cout << (digits.getfDigits())[0].getADC() << " est ADC du premier digit" << endl;
        
        Double_t arg(-2); // disable printout

        for ( int cathode = 0; cathode < 2; ++cathode ) //On boucle sur les deux plans de cathodes
        {
            fitter->Clear();
            fitter->ExecuteCommand("SET PRINT",&arg,1);
            if(vecxmin[cathode]>vecxmax[cathode]){
                vecxmin[cathode] = -40.;
                vecxmax[cathode] = +40.;
            }
            if(vecymin[cathode]>vecymax[cathode]){
                vecymin[cathode] = -20.;
                vecymax[cathode] = +20.;
            }
            cout << "On va set les params" << endl;
            
            MyObjectParametersGauss parameters;
            parameters.setfCathode(cathode);
            parameters.setfChargetot(chargetot[cathode]);
            parameters.setfSigx(sigx);
            parameters.setfSigy(sigy);
            parameters.setfNbdigits(precluster.size());
            
            
            TObjArray userObjects;
            userObjects.Add(&digits);
            userObjects.Add(&parameters);
            fitter->SetObjectFit(&userObjects);
            
            fitter->SetParameter(0,"cluster X position",xhit[cathode],stepX,vecxmin[cathode],vecxmax[cathode]);
            fitter->SetParameter(1,"cluster Y position",yhit[cathode],stepY,vecymin[cathode],vecymax[cathode]);
            

            Double_t stratArg(2);
            fitter->ExecuteCommand("SET STR",&stratArg,1);
            Int_t val = fitter->ExecuteCommand("MIGRAD",0,0);
            if ( val && chargetot[cathode] != 0 )
            {
                //Fit failed with robust strategy 2 try balanced strategy 1
                Double_t stratArg(1);
                fitter->ExecuteCommand("SET STR",&stratArg,1);
                Int_t val2 = fitter->ExecuteCommand("MIGRAD",0,0);
                if ( val2 ){
                    // fit failed. Using COG results, with big errors
                    cout << "Fit failed on viable cathode. Using COG results for cluster." << endl;
                        return clustertmpCOG;
                }
            }

            xhit[cathode] = fitter->GetParameter(0);
            yhit[cathode] = fitter->GetParameter(1);
            xsize[cathode] = fitter->GetParError(0);
            ysize[cathode] = fitter->GetParError(1);

            Double_t amin, edm, errdef;
            Int_t nvpar, nparx;

            fitter->GetStats(amin, edm, errdef, nvpar, nparx);

            Double_t chi2 = amin;

            printf("Cluster fitted to (x,y)=(%e,%e) (xerr,yerr)=(%e,%e) \n chi2=%e ndf=%d",
                            xhit[cathode],yhit[cathode],
                            xsize[cathode],ysize[cathode],chi2,fitter->GetNumberFreeParameters());
        }

        //Détermination de la position finale du hit à l'aide de la résolution sur chaque cathode
        Double_t xhitfinal = xhit[0];
        Double_t yhitfinal = yhit[1];
        Double_t ex = xsize[0];
        Double_t ey = ysize[1];

        double timestamp = precluster[0].getTimeStamp();  // On récupère le timestamp du premier digit du precluster


        //On remplit un cluster avec les infos nécessaires
        cluster.setx(xhitfinal);
        cluster.sety(yhitfinal);
        cluster.setex(ex);
        cluster.setey(ey);
        cluster.settimestamp(timestamp);    // En première approx. on peut affecter le timestamp du pre;ier digit au cluster

                printf("\n\nCluster par gaussian fit:\nCluster multiplicity %ld \n(x,y)=(%e,%e) \nprecision=(ex,ey)=(%e,%e) \ntimestamp = %lf \n\n",precluster.size(),cluster.getx(),cluster.gety(),cluster.getex(),cluster.getey(),cluster.gettimestamp());

        return cluster;
}


//_________________________________________________________________________________________________
void FitFunctionGaussianFitClean(Int_t& /*notused*/, Double_t* /*notused*/,
                            Double_t& sum, Double_t* par,
                             Int_t /*notused*/)
  {
    /// Chi2 Function to minimize: Gaussian charge distribution in 2 dimensions
          
              //Méthode propre. On amène un pointeur vers les digits aue l'on peut lire un par un. Et un TArray avec les autres parametres (cathode, chgtot, sigx, sigy, nbdigits)

        TObjArray* userObjects = static_cast<TObjArray*>(TVirtualFitter::GetFitter()->GetObjectFit());
          
          MyObjectDigits* digits = static_cast<MyObjectDigits*>(userObjects->At(0));
          MyObjectParametersGauss* parameters = static_cast<MyObjectParametersGauss*>(userObjects->At(1));
          
          
          Clustering cluster;

        sum = 0.0;
          
          
          int cathode = parameters->getfCathode();
          double chargetot = parameters->getfChargetot();
          double sigx = parameters->getfSigx();
          double sigy = parameters->getfSigy();
          int nbdigits = parameters->getfNbdigits();
          
        for(Int_t i = 0; i < nbdigits; ++i){
          
            Digit digitbeingread = (digits->getfDigits())[i];
            
            int detid = digitbeingread.getDetID();
            int padid = digitbeingread.getPadID();
            int ADC = digitbeingread.getADC();
            
            mapping::Segmentation pad(detid);

            if ( cathode == pad.isBendingPad(padid) ) //On regarde à quelle cathode appartient le pad
            {
                  //On définit le vecteur position et taille du pad en question
                std::vector<Double_t> padPosition = {pad.padPositionX(padid), pad.padPositionY(padid)};
                std::vector<Double_t> padSize = {pad.padSizeX(padid), pad.padSizeY(padid)};

                std::vector<Double_t> lowerLeft = {par[0] - padPosition[0] - 0.5*padSize[0], par[1] - padPosition[1] - 0.5*padSize[1]};
                std::vector<Double_t> upperRight= {lowerLeft[0] + padSize[0], lowerLeft[1] + padSize[1]};

                Double_t qfit = chargetot*cluster.IntGaussXY(lowerLeft[0], lowerLeft[1], upperRight[0], upperRight[1], sigx, sigy);

                Double_t qtrue = ADC;

                Double_t delta = qtrue-qfit;

                sum += pow(delta,2)/qtrue;
            }
        }
}

//_________________________________________________________________________________________________
Double_t Clustering::IntGaussXY(Double_t x1, Double_t y1, Double_t x2, Double_t y2, Double_t sigx, Double_t sigy)
{
    Double_t pitch = 0.25;
    x1 /= pitch;
    x2 /= pitch;
    y1 /= pitch;
    y2 /= pitch;
    
    Double_t erfx = 0.5*( erf( (x2)/(sigx*sqrt(2)) ) - erf( (x1)/(sigx*sqrt(2)) ) );
    Double_t erfy = 0.5*( erf( (y2)/(sigy*sqrt(2)) ) - erf( (y1)/(sigy*sqrt(2)) ) );
    Double_t integral = erfx*erfy;
    
    return integral;
}

//_________________________________________________________________________________________________






// FONCTIONS POUR CLUSTERING DOUBLE GAUSS FIT






//_________________________________________________________________________________________________
void Clustering::runFinderDoubleGaussianFit(std::vector<PreClusterStruct>& preClusters, std::vector<Cluster>& clusters)
{
    cout << "preClusters.size():" << preClusters.size() << endl;
    int sizedigit;
    int jete = 0;
    Digit digittmp;
    Cluster clustertmpCOG;
    Cluster clustertmp;
    
    for (Int_t i = 0; i < preClusters.size(); i++)   //On boucle sur les preclusters identifiés
    {
        if(preClusters[i].nDigits == 1){         // Je vire les preclusters qui n'ont qu'un seul digit
            jete++;
            continue;
        }
        cout << "Precluster avec " << preClusters[i].nDigits << " digits." << endl;  //Nombre de digits dans le precluster aue l'on regarde
        
        std::vector<Digit> preclustertmp;
        const Digit* ptrdigit = preClusters[i].digits;
        
        for (int j=0; j<preClusters[i].nDigits; j++){          //Je récupère les digits et les mets dans un   vecteur de digits, preclustertmp
            digittmp = *ptrdigit;
            cout << "Digit numero:" << j << endl;
            cout << "digittmp.getADC():" << digittmp.getADC() << endl;
            cout << "digittmp.getDetID():" << digittmp.getDetID() << endl;
            cout << "digittmp.getPadID():" << digittmp.getPadID() << endl;
//            cout << "type digittmp:" << typeid(digittmp).name() << endl;
//            cout << "type ptrdigit:" << typeid(ptrdigit).name() << endl;
//            cout << "ptrdigit:" << ptrdigit << endl;
            preclustertmp.push_back(digittmp);
            sizedigit = sizeof(digittmp);
            ptrdigit++;
        }
        
        //Je clusterise mon vecteur de digits, preclustertmp et l'ajoute au vecteur de clusters, nommé clusters.
        
       // cout << "preclustertmp[0].getADC():" << preclustertmp[0].getADC() << endl;
        clustertmpCOG = FinderCOG(preclustertmp);
      //  cout << "COORD X DU CLUSTER AJOUTÉ A CLUSTERS: " << clustertmp.getx() << endl;
        if(preclustertmp.size() < 3){
            clustertmp = clustertmpCOG;
        }
        else{
        clustertmp = ComputePositionDoubleGaussianFitClean(preclustertmp, clustertmpCOG);
        }
        
        clusters.push_back(clustertmp);
    }
    
    cout << "Nombre de preclusters au départ:" << preClusters.size() << endl;
    cout << "Nombre de preclusters de taille 1, abandonnés:" << jete << endl;
    
}

//_________________________________________________________________________________________________
Clustering::Cluster Clustering::ComputePositionDoubleGaussianFitClean(std::vector<Digit> &precluster, Clustering::Cluster clustertmpCOG)
{
  /// Compute the position of the given cluster, by fitting a Double Gaussian
  /// charge distribution to it

        Double_t vecxmin[2] = {1E3, 1E3};
        Double_t vecxmax[2] = {-1E3, -1E3};
        Double_t vecymin[2] = {1E3, 1E3};
        Double_t vecymax[2] = {-1E3, -1E3};
        
      TVirtualFitter* fitter = TVirtualFitter::Fitter(0,2);
      fitter->SetFCN(FitFunctionDoubleGaussianFitClean);

      cout << "\n\n==========\nRunning MinuitDoubleGaussian Algorithm.\n\n" << endl;

      Cluster cluster;

      //Création des vecteurs position du hit sur chaque cqthode
      Double_t xhit[2] = { clustertmpCOG.getx(), clustertmpCOG.getx() };
      Double_t yhit[2] = { clustertmpCOG.gety(), clustertmpCOG.gety() };

        //Création des vecteurs de charge totale, multiplicité et tailles sur chaque cathode
           Double_t chargetot[2] = { 0.0, 0.0 };
           Double_t multiplicity[2] = { 0.0, 0.0 };
           Double_t xsize[2] = { 0.0, 0.0 };
           Double_t ysize[2] = { 0.0, 0.0 };

       //Détermination des vecteurs de charge totale, multiplicité et taille qui devient précision
       for ( Int_t i = 0; i < precluster.size(); ++i ){
           int detid = precluster[i].getDetID();
           int padid = precluster[i].getPadID();
           
           mapping::Segmentation pad(detid); // = mapping::Segmentation(detid);

           for ( int cathode = 0; cathode < 2; ++cathode ) //On boucle sur les deux plans de cathodes
           {

               if ( cathode == pad.isBendingPad(padid) ) //On regarde à quelle cathode appartient le pad
                       {
                           std::vector<Double_t> padPosition = {pad.padPositionX(padid), pad.padPositionY(padid)};
                           std::vector<Double_t> padSize = {pad.padSizeX(padid), pad.padSizeY(padid)};
                           
                           vecxmin[cathode] = TMath::Min(padPosition[0]-0.5*padSize[0],vecxmin[cathode]);
                           vecxmax[cathode] = TMath::Max(padPosition[0]+0.5*padSize[0],vecxmax[cathode]);
                           vecymin[cathode] = TMath::Min(padPosition[1]-0.5*padSize[1],vecymin[cathode]);
                           vecymax[cathode] = TMath::Max(padPosition[1]+0.5*padSize[1],vecymax[cathode]);

                           chargetot[cathode] += precluster[i].getADC();
                           multiplicity[cathode] += 1;
                       }
           }
       }
        

      Float_t stepX = 0.00001; // cm
      Float_t stepY = 0.00001; // cm
        Double_t sig1x = 1.1059;
        Double_t sig1y = 1.1350;
        Double_t sig2x = 0.5702;
        Double_t sig2y = 0.5785;
        Double_t chgfracx = 0.6233;
        Double_t chgfracy = 0.6519;
        
        
        MyObjectDigits digits;
        
        cout << "On va set digits" << endl;
        digits.setfDigits(precluster);
        cout << (digits.getfDigits())[0].getADC() << " est ADC du premier digit" << endl;
        
        Double_t arg(-2); // disable printout

        for ( int cathode = 0; cathode < 2; ++cathode ) //On boucle sur les deux plans de cathodes
        {
            fitter->Clear();
            fitter->ExecuteCommand("SET PRINT",&arg,1);
            if(vecxmin[cathode]>vecxmax[cathode]){
                vecxmin[cathode] = -40.;
                vecxmax[cathode] = +40.;
            }
            if(vecymin[cathode]>vecymax[cathode]){
                vecymin[cathode] = -20.;
                vecymax[cathode] = +20.;
            }
            cout << "On va set les params" << endl;
            
            MyObjectParametersDoubleGauss parameters;
            parameters.setfCathode(cathode);
            parameters.setfChargetot(chargetot[cathode]);
            parameters.setfSig1x(sig1x);
            parameters.setfSig1y(sig1y);
            parameters.setfSig2x(sig2x);
            parameters.setfSig2y(sig2y);
            parameters.setfChgfracx(chgfracx);
            parameters.setfChgfracy(chgfracy);
            parameters.setfNbdigits(precluster.size());
            
            
            TObjArray userObjects;
            userObjects.Add(&digits);
            userObjects.Add(&parameters);
            fitter->SetObjectFit(&userObjects);
            
            fitter->SetParameter(0,"cluster X position",xhit[cathode],stepX,vecxmin[cathode],vecxmax[cathode]);
            fitter->SetParameter(1,"cluster Y position",yhit[cathode],stepY,vecymin[cathode],vecymax[cathode]);
            

            Double_t stratArg(2);
            fitter->ExecuteCommand("SET STR",&stratArg,1);
            Int_t val = fitter->ExecuteCommand("MIGRAD",0,0);
            if ( val && chargetot[cathode] != 0 )
            {
                //Fit failed with robust strategy 2 try balanced strategy 1
                Double_t stratArg(1);
                fitter->ExecuteCommand("SET STR",&stratArg,1);
                Int_t val2 = fitter->ExecuteCommand("MIGRAD",0,0);
                if ( val2 ){
                    // fit failed. Using COG results, with big errors
                    cout << "Fit failed on viable cathode. Using COG results for cluster." << endl;
                        return clustertmpCOG;
                }
            }

            xhit[cathode] = fitter->GetParameter(0);
            yhit[cathode] = fitter->GetParameter(1);
            xsize[cathode] = fitter->GetParError(0);
            ysize[cathode] = fitter->GetParError(1);

            Double_t amin, edm, errdef;
            Int_t nvpar, nparx;

            fitter->GetStats(amin, edm, errdef, nvpar, nparx);

            Double_t chi2 = amin;

            printf("Cluster fitted to (x,y)=(%e,%e) (xerr,yerr)=(%e,%e) \n chi2=%e ndf=%d",
                            xhit[cathode],yhit[cathode],
                            xsize[cathode],ysize[cathode],chi2,fitter->GetNumberFreeParameters());
        }

        //Détermination de la position finale du hit à l'aide de la résolution sur chaque cathode
        Double_t xhitfinal = xhit[0];
        Double_t yhitfinal = yhit[1];
        Double_t ex = xsize[0];
        Double_t ey = ysize[1];

        double timestamp = precluster[0].getTimeStamp();  // On récupère le timestamp du premier digit du precluster


        //On remplit un cluster avec les infos nécessaires
        cluster.setx(xhitfinal);
        cluster.sety(yhitfinal);
        cluster.setex(ex);
        cluster.setey(ey);
        cluster.settimestamp(timestamp);    // En première approx. on peut affecter le timestamp du pre;ier digit au cluster

                printf("\n\nCluster par double gaussian fit:\nCluster multiplicity %ld \n(x,y)=(%e,%e) \nprecision=(ex,ey)=(%e,%e) \ntimestamp = %lf \n\n",precluster.size(),cluster.getx(),cluster.gety(),cluster.getex(),cluster.getey(),cluster.gettimestamp());

        return cluster;
}


//_________________________________________________________________________________________________
void FitFunctionDoubleGaussianFitClean(Int_t& /*notused*/, Double_t* /*notused*/,
                            Double_t& sum, Double_t* par,
                             Int_t /*notused*/)
  {
    /// Chi2 Function to minimize: Mathieson charge distribution in 2 dimensions
      
          //Méthode bourrine. Définir un vecteur de paramètres énorme:
          // x   y   cathode   chgtot   sig1s   sig2s   chgfracs   digit par digit (Det Pad ADC)


      
     TObjArray* userObjects = static_cast<TObjArray*>(TVirtualFitter::GetFitter()->GetObjectFit());
        
        MyObjectDigits* digits = static_cast<MyObjectDigits*>(userObjects->At(0));
        MyObjectParametersDoubleGauss* parameters = static_cast<MyObjectParametersDoubleGauss*>(userObjects->At(1));
        
        
        Clustering cluster;

      sum = 0.0;
        
        
        int cathode = parameters->getfCathode();
        double chargetot = parameters->getfChargetot();
        double sig1x = parameters->getfSig1x();
        double sig1y = parameters->getfSig1y();
        double sig2x = parameters->getfSig2x();
        double sig2y = parameters->getfSig2y();
        double chgfracx = parameters->getfChgfracx();
        double chgfracy = parameters->getfChgfracy();
        int nbdigits = parameters->getfNbdigits();
        
      for(Int_t i = 0; i < nbdigits; ++i){
        
          Digit digitbeingread = (digits->getfDigits())[i];
          
          int detid = digitbeingread.getDetID();
          int padid = digitbeingread.getPadID();
          int ADC = digitbeingread.getADC();
        
    mapping::Segmentation pad(detid); // = mapping::Segmentation(detid);

    if ( cathode == pad.isBendingPad(padid) ) //On regarde à quelle cathode appartient le pad
    {
          //On définit le vecteur position et taille du pad en question
        std::vector<Double_t> padPosition = {pad.padPositionX(padid), pad.padPositionY(padid)};
        std::vector<Double_t> padSize = {pad.padSizeX(padid), pad.padSizeY(padid)};

        std::vector<Double_t> lowerLeft = {par[0] - padPosition[0] - 0.5*padSize[0], par[1] - padPosition[1] - 0.5*padSize[1]};
        std::vector<Double_t> upperRight= {lowerLeft[0] + padSize[0], lowerLeft[1] + padSize[1]};

        Double_t qfit = chargetot*cluster.IntDoubleGaussXY(lowerLeft[0], lowerLeft[1], upperRight[0], upperRight[1], sig1x, sig1y, sig2x, sig2y, chgfracx, chgfracy);

        Double_t qtrue = ADC;

        Double_t delta = qtrue-qfit;

        sum += pow(delta,2)/qtrue;
    }
  }
}



//_________________________________________________________________________________________________
Double_t Clustering::IntDoubleGaussXY(Double_t x1, Double_t y1, Double_t x2, Double_t y2, Double_t sig1x, Double_t sig1y, Double_t sig2x, Double_t sig2y, Double_t chgfracx, Double_t chgfracy)
{
    Double_t pitch = 0.25;
    x1 /= pitch;
    x2 /= pitch;
    y1 /= pitch;
    y2 /= pitch;
    
    Double_t erf1x = 0.5*( erf( (x2)/(sig1x*sqrt(2)) ) - erf( (x1)/(sig1x*sqrt(2)) ) );
    Double_t erf1y = 0.5*( erf( (y2)/(sig1y*sqrt(2)) ) - erf( (y1)/(sig1y*sqrt(2)) ) );
    Double_t erf2x = 0.5*( erf( (x2)/(sig2x*sqrt(2)) ) - erf( (x1)/(sig2x*sqrt(2)) ) );
    Double_t erf2y = 0.5*( erf( (y2)/(sig2y*sqrt(2)) ) - erf( (y1)/(sig2y*sqrt(2)) ) );
    Double_t erfx = chgfracx*erf1x/(chgfracx+1) + erf2x/(chgfracx+1);
    Double_t erfy = chgfracy*erf1y/(chgfracy+1) + erf2y/(chgfracy+1);
    Double_t integral = erfx*erfy;
    
    return integral;
}













// ---------------------------------------
// FONCTIONS OBSOLETES - GARDEES AU CAS OÙ
// ---------------------------------------















//_________________________________________________________________________________________________
Clustering::Cluster Clustering::FinderSimpleFit(std::vector<Digit> &precluster, Clustering::Cluster clustertmpCOG)
{
    
// FIT MANUEL - POUR AVOIR LE FIT SOUS MINUIT VOIR ComputePosition
    cout << "\n\n==========\nRunning SimpleFit Algorithm\n\n" << endl;
    
    Cluster cluster;

    //Création des vecteurs position du hit sur chaque cqthode
    Double_t xhit[2] = { clustertmpCOG.getx(), clustertmpCOG.getx() };
    Double_t yhit[2] = { clustertmpCOG.gety(), clustertmpCOG.gety() };
    
  //  cout << "clustertmpCOG.getx()" << clustertmpCOG.getx() << endl;
    
    //Création du float porteur de la coordonnée la plus précise
    Double_t xhitfinal = 0;
    Double_t yhitfinal = 0;
    
    //Création des vecteurs de charge totale, multiplicité et tailles sur chaque cathode
    Double_t chargetot[2] = { 0.0, 0.0 };
    Double_t multiplicity[2] = { 0.0, 0.0 };
    Double_t xsize[2] = { 0.0, 0.0 };
    Double_t ysize[2] = { 0.0, 0.0 };
    
    //Détermination des vecteurs de charge totale, multiplicité et taille qui devient précision
    for ( Int_t i = 0; i < precluster.size(); ++i ){
        int detid = precluster[i].getDetID();
        int padid = precluster[i].getPadID();
        mapping::Segmentation pad(detid); // = mapping::Segmentation(detid);
            
        for ( int cathode = 0; cathode < 2; ++cathode ) //On boucle sur les deux plans de cathodes
        {
                
            if ( cathode == pad.isBendingPad(padid) ) //On regarde à quelle cathode appartient le pad
                    {
                        std::vector<Double_t> padPosition = {pad.padPositionX(padid), pad.padPositionY(padid)};
                        std::vector<Double_t> padSize = {pad.padSizeX(padid), pad.padSizeY(padid)};

                        chargetot[cathode] += precluster[i].getADC();
                        multiplicity[cathode] += 1;
                        xsize[cathode] += padSize[0];
                        ysize[cathode] += padSize[1];
                    }
        }
    }
    for ( int cathode = 0; cathode < 2; ++cathode ) //On boucle sur les deux plans de cathodes
    {
        if ( multiplicity[cathode] != 0 )
        {
          xsize[cathode] /= multiplicity[cathode];
          ysize[cathode] /= multiplicity[cathode];
        }
    }
    
    //Minimisation du Chi2 en 2D sur chaque cathode
    for ( int cathode = 0; cathode < 2; ++cathode ) //On boucle sur les deux plans de cathodes
    {
    
        Double_t Kx3 = 0.5085;
        Double_t Ky3 = 0.5840;
        
        Double_t pasxhit = 0.00001;
        Double_t pasyhit = 0.00001;
        
        int compteurboucle = 0;
        
        Double_t ChiDeux = Chi2Mathieson(cathode, chargetot[cathode], precluster, clustertmpCOG, xhit[cathode], yhit[cathode], Kx3, Ky3);
        Double_t OldChiDeux = 1;
        Double_t NewChiDeux = 100;
        Double_t deltaxhit = 100;
        Double_t deltayhit = 100;
        Double_t limitstrict = 0.000001;
        Double_t limitlean = 0.0001;
        Double_t limitx = 0.;
        Double_t limity = 0.;
        
        if(cathode == 0){
            limitx = limitstrict;
            limity = limitlean;
        }
        else{
            limity = limitstrict;
            limitx = limitlean;
        }
        
        while ((sqrt(pow((NewChiDeux-OldChiDeux)/OldChiDeux,2)) > 0.001)||(abs(deltaxhit)>limitx)||(abs(deltayhit)>limity)){
            compteurboucle++;
            
            //Iteration and improvement of estimate
            Double_t dChidxhit = (Chi2Mathieson(cathode, chargetot[cathode], precluster, clustertmpCOG, xhit[cathode]+pasxhit, yhit[cathode], Kx3, Ky3)-Chi2Mathieson(cathode, chargetot[cathode], precluster, clustertmpCOG, xhit[cathode], yhit[cathode], Kx3, Ky3))/pasxhit;
            Double_t dChidyhit = (Chi2Mathieson(cathode, chargetot[cathode], precluster, clustertmpCOG, xhit[cathode], yhit[cathode]+pasyhit, Kx3, Ky3)-Chi2Mathieson(cathode, chargetot[cathode], precluster, clustertmpCOG, xhit[cathode], yhit[cathode], Kx3, Ky3))/pasyhit;
            
//            float dChidKx3 = (Chi2Mathieson(cathode, precluster, clustertmp, xhit[cathode], yhit[cathode], Kx3+pasKx3, Ky3)-Chi2Mathieson(cathode, precluster, clustertmp, xhit[cathode], yhit[cathode], Kx3, Ky3))/pasKx3;


//            float dChidKy3 = (Chi2Mathieson(cathode, precluster, clustertmp, xhit[cathode], yhit[cathode], Kx3, Ky3+pasKy3)-Chi2Mathieson(cathode, precluster, clustertmp, xhit[cathode], yhit[cathode], Kx3, Ky3))/pasKy3;

            deltaxhit = pasxhit*dChidxhit*1;
            deltayhit = pasyhit*dChidyhit*1;
            xhit[cathode] = xhit[cathode] - pasxhit*dChidxhit*1;
            yhit[cathode] = yhit[cathode] - pasyhit*dChidyhit*1;
//            Kx3 = Kx3 - pasKx3*dChidKx3;
//            Ky3 = Ky3 - pasKy3*dChidKy3;
            
            
            // Update ChiDeux
            OldChiDeux = ChiDeux;
            ChiDeux = Chi2Mathieson(cathode, chargetot[cathode], precluster, clustertmpCOG, xhit[cathode], yhit[cathode], Kx3, Ky3);
            NewChiDeux = ChiDeux;
            cout << "\nBoucle While numero " << compteurboucle << endl;
            cout << "xhit: " << xhit[cathode] << endl;
            cout << "yhit: " << yhit[cathode] << endl;
            cout << "deltaxhit: " << deltaxhit << endl;
            cout << "deltayhit: " << deltayhit << endl;
            cout << "OldChiDeux: " << OldChiDeux << endl;
            cout << "NewChiDeux: " << NewChiDeux << endl;
            
            if(compteurboucle > 10000){
                cout << "Déjà 10.000 boucles pour trouver Chi2, utilisation résultats COG..." << endl;
                return clustertmpCOG;
            }
        }
        
        cout << "\n\ncathode " << cathode << endl;
        cout << "Boucle Finale numero " << compteurboucle << endl;
        cout << "xhit final: " << xhit[cathode] << endl;
        cout << "yhit final: " << yhit[cathode] << endl;
        cout << "Kx3: " << Kx3 << endl;
        cout << "Ky3: " << Ky3 << endl;
        cout << "Chi2 Final: " << NewChiDeux << endl;
    
    }
    
    
    //Détermination de la position finale du hit à l'aide de la résolution sur chaque cathode
    xhitfinal = xhit[0];
    yhitfinal = yhit[1];
    Double_t ex = xsize[0];
    Double_t ey = ysize[1];
    
    double timestamp = precluster[0].getTimeStamp();  // On récupère le timestamp du premier digit du precluster
    
    
    //On remplit un cluster avec les infos nécessaires
    cluster.setx(xhitfinal);
    cluster.sety(yhitfinal);
    cluster.setex(ex);
    cluster.setey(ey);
    cluster.settimestamp(timestamp);    // En première approx. on peut affecter le timestamp du pre;ier digit au cluster
    
    printf("\n\nCluster par COG:\nCluster multiplicity %ld \n(x,y)=(%e,%e) \nprecision=(ex,ey)=(%e,%e) \ntimestamp = %lf \n\n",precluster.size(),clustertmpCOG.getx(),clustertmpCOG.gety(),clustertmpCOG.getex(),clustertmpCOG.getey(),clustertmpCOG.gettimestamp());
    
            printf("\n\nCluster par simple fit:\nCluster multiplicity %ld \n(x,y)=(%e,%e) \nprecision=(ex,ey)=(%e,%e) \ntimestamp = %lf \n\n",precluster.size(),cluster.getx(),cluster.gety(),cluster.getex(),cluster.getey(),cluster.gettimestamp());
    
    return cluster;

}


//_________________________________________________________________________________________________
Double_t Clustering::Chi2Mathieson(int cathode, Double_t chargetot, std::vector<Digit> &precluster, Clustering::Cluster clustertmp, Double_t xhit, Double_t yhit, Double_t Kx3, Double_t Ky3)
{
    
    Double_t sum = 0;
    for(Int_t i = 0; i < precluster.size(); ++i){
        
        int detid = precluster[i].getDetID();
        int padid = precluster[i].getPadID();
        mapping::Segmentation pad(detid); // = mapping::Segmentation(detid);
        
        if ( cathode == pad.isBendingPad(padid) ) //On regarde à quelle cathode appartient le pad
        {
              //On définit le vecteur position et taille du pad en question
            std::vector<Double_t> padPosition = {pad.padPositionX(padid), pad.padPositionY(padid)};
            std::vector<Double_t> padSize = {pad.padSizeX(padid), pad.padSizeY(padid)};
            
            std::vector<Double_t> lowerLeft = {xhit - padPosition[0] - 0.5*padSize[0], yhit - padPosition[1] - 0.5*padSize[1]};
            std::vector<Double_t> upperRight= {lowerLeft[0] + padSize[0], lowerLeft[1] + padSize[1]};
            
            Double_t qfit = chargetot*IntMathiesonXY(lowerLeft[0], lowerLeft[1], upperRight[0], upperRight[1], Kx3, Ky3);

            Double_t qtrue = precluster[i].getADC();
            
            Double_t delta = qtrue-qfit;
            
            sum += pow(delta,2)/qtrue;
        }
    }
    
    return sum;
}
//_________________________________________________________________________________________________
void FitFunction(Int_t& /*notused*/, Double_t* /*notused*/,
                            Double_t& sum, Double_t* par,
                             Int_t /*notused*/)
  {
    /// Chi2 Function to minimize: Mathieson charge distribution in 2 dimensions
      
          //Méthode bourrine. Définir un vecteur de paramètres énorme:
          // x   y   cathode   chgtot   Kx3   Ky3   digit par digit (Det Pad ADC)

//    TObjArray* userObjects = static_cast<TObjArray*>(TVirtualFitter::GetFitter()->GetObjectFit());
//
//    std::vector<Digit> precluster = static_cast<std::vector<Digit>>(userObjects->At(0));
//    int cathode = static_cast<int>(userObjects->At(1));
//    Double_t chargetot = static_cast<Double_t>(userObjects->At(2));
//    Double_t Kx3 = static_cast<Double_t>(userObjects->At(3));
//    Double_t Ky3 = static_cast<Double_t>(userObjects->At(4));
      
      Clustering cluster;

    sum = 0.0;
      
      int cathode = par[2];
      int chargetot = par[3];
      double Kx3 = par[4];
      double Ky3 = par[5];
      int nbdigits = par[6];

    for(Int_t i = 0; i < nbdigits; ++i){

    int detid = par[7+(3*i)];
    int padid = par[8+(3*i)];
    int ADC = par[9+(3*i)];
        
    mapping::Segmentation pad(detid); // = mapping::Segmentation(detid);

    if ( cathode == pad.isBendingPad(padid) ) //On regarde à quelle cathode appartient le pad
    {
          //On définit le vecteur position et taille du pad en question
        std::vector<Double_t> padPosition = {pad.padPositionX(padid), pad.padPositionY(padid)};
        std::vector<Double_t> padSize = {pad.padSizeX(padid), pad.padSizeY(padid)};

        std::vector<Double_t> lowerLeft = {par[0] - padPosition[0] - 0.5*padSize[0], par[1] - padPosition[1] - 0.5*padSize[1]};
        std::vector<Double_t> upperRight= {lowerLeft[0] + padSize[0], lowerLeft[1] + padSize[1]};

        Double_t qfit = chargetot*cluster.IntMathiesonXY(lowerLeft[0], lowerLeft[1], upperRight[0], upperRight[1], Kx3, Ky3);

        Double_t qtrue = ADC;

        Double_t delta = qtrue-qfit;

        sum += pow(delta,2)/qtrue;
    }
  }
}

//_________________________________________________________________________________________________
Clustering::Cluster Clustering::ComputePosition(std::vector<Digit> &precluster, Clustering::Cluster clustertmpCOG)
{
  /// Compute the position of the given cluster, by fitting a Mathieson
  /// charge distribution to it

    int nparams = 6+(3*precluster.size());
    Double_t vecxmin[2] = {1E3, 1E3};
    Double_t vecxmax[2] = {-1E3, -1E3};
    Double_t vecymin[2] = {1E3, 1E3};
    Double_t vecymax[2] = {-1E3, -1E3};
    
  TVirtualFitter* fitter = TVirtualFitter::Fitter(0,nparams);
  fitter->SetFCN(FitFunction);

  cout << "\n\n==========\nRunning MinuitMathieson Algorithm.\n\n" << endl;

  Cluster cluster;

  //Création des vecteurs position du hit sur chaque cqthode
  Double_t xhit[2] = { clustertmpCOG.getx(), clustertmpCOG.getx() };
  Double_t yhit[2] = { clustertmpCOG.gety(), clustertmpCOG.gety() };

    //Création des vecteurs de charge totale, multiplicité et tailles sur chaque cathode
       Double_t chargetot[2] = { 0.0, 0.0 };
       Double_t multiplicity[2] = { 0.0, 0.0 };
       Double_t xsize[2] = { 0.0, 0.0 };
       Double_t ysize[2] = { 0.0, 0.0 };

   //Détermination des vecteurs de charge totale, multiplicité et taille qui devient précision
   for ( Int_t i = 0; i < precluster.size(); ++i ){
       int detid = precluster[i].getDetID();
       int padid = precluster[i].getPadID();
       
       mapping::Segmentation pad(detid); // = mapping::Segmentation(detid);

       for ( int cathode = 0; cathode < 2; ++cathode ) //On boucle sur les deux plans de cathodes
       {

           if ( cathode == pad.isBendingPad(padid) ) //On regarde à quelle cathode appartient le pad
                   {
                       std::vector<Double_t> padPosition = {pad.padPositionX(padid), pad.padPositionY(padid)};
                       std::vector<Double_t> padSize = {pad.padSizeX(padid), pad.padSizeY(padid)};
                       
                       vecxmin[cathode] = TMath::Min(padPosition[0]-0.5*padSize[0],vecxmin[cathode]);
                       vecxmax[cathode] = TMath::Max(padPosition[0]+0.5*padSize[0],vecxmax[cathode]);
                       vecymin[cathode] = TMath::Min(padPosition[1]-0.5*padSize[1],vecymin[cathode]);
                       vecymax[cathode] = TMath::Max(padPosition[1]+0.5*padSize[1],vecymax[cathode]);

                       chargetot[cathode] += precluster[i].getADC();
                       multiplicity[cathode] += 1;
                   }
       }
   }
    

  Float_t stepX = 0.00001; // cm
  Float_t stepY = 0.00001; // cm
    Double_t Kx3 = 0.5085;
    Double_t Ky3 = 0.5840;

    Double_t arg(-1); // disable printout

    for ( int cathode = 0; cathode < 2; ++cathode ) //On boucle sur les deux plans de cathodes
    {
        fitter->Clear();
        fitter->ExecuteCommand("SET PRINT",&arg,1);
        if(vecxmin[cathode]>vecxmax[cathode]){
            vecxmin[cathode] = -40.;
            vecxmax[cathode] = +40.;
        }
        if(vecymin[cathode]>vecymax[cathode]){
            vecymin[cathode] = -20.;
            vecymax[cathode] = +20.;
        }
        
        fitter->SetParameter(0,"cluster X position",xhit[cathode],stepX,vecxmin[cathode],vecxmax[cathode]);
        fitter->SetParameter(1,"cluster Y position",yhit[cathode],stepY,vecymin[cathode],vecymax[cathode]);
        fitter->SetParameter(2,"cathode",cathode,0,0,0);
        fitter->SetParameter(3,"chargetot",chargetot[cathode],0,0,0);
        fitter->SetParameter(4,"Kx3",Kx3,0,0,0);
        fitter->SetParameter(5,"Ky3",Ky3,0,0,0);
        fitter->SetParameter(6,"NbDigits",precluster.size(),0,0,0);
        
        for(int j=0; j<precluster.size(); j++){
            fitter->SetParameter(7+(3*j),"DetID",precluster[j].getDetID(),0,0,0);
            fitter->SetParameter(8+(3*j),"PadID",precluster[j].getPadID(),0,0,0);
            fitter->SetParameter(9+(3*j),"ADC",precluster[j].getADC(),0,0,0);
        }

        Double_t stratArg(2);
        fitter->ExecuteCommand("SET STR",&stratArg,1);
        Int_t val = fitter->ExecuteCommand("MIGRAD",0,0);
        if ( val && chargetot[cathode] != 0 )
        {
            //Fit failed with robust strategy 2 try balanced strategy 1
            Double_t stratArg(1);
            fitter->ExecuteCommand("SET STR",&stratArg,1);
            Int_t val2 = fitter->ExecuteCommand("MIGRAD",0,0);
            if ( val2 ){
                // fit failed. Using COG results, with big errors
                cout << "Fit failed on viable cathode. Using COG results for cluster." << endl;
                    return clustertmpCOG;
            }
        }

        xhit[cathode] = fitter->GetParameter(0);
        yhit[cathode] = fitter->GetParameter(1);
        xsize[cathode] = fitter->GetParError(0);
        ysize[cathode] = fitter->GetParError(1);

        Double_t amin, edm, errdef;
        Int_t nvpar, nparx;

        fitter->GetStats(amin, edm, errdef, nvpar, nparx);

        Double_t chi2 = amin;

        printf("Cluster fitted to (x,y)=(%e,%e) (xerr,yerr)=(%e,%e) \n chi2=%e ndf=%d",
                        xhit[cathode],yhit[cathode],
                        xsize[cathode],ysize[cathode],chi2,fitter->GetNumberFreeParameters());
    }

    //Détermination de la position finale du hit à l'aide de la résolution sur chaque cathode
    Double_t xhitfinal = xhit[0];
    Double_t yhitfinal = yhit[1];
    Double_t ex = xsize[0];
    Double_t ey = ysize[1];

    double timestamp = precluster[0].getTimeStamp();  // On récupère le timestamp du premier digit du precluster


    //On remplit un cluster avec les infos nécessaires
    cluster.setx(xhitfinal);
    cluster.sety(yhitfinal);
    cluster.setex(ex);
    cluster.setey(ey);
    cluster.settimestamp(timestamp);    // En première approx. on peut affecter le timestamp du pre;ier digit au cluster

            printf("\n\nCluster par simple fit:\nCluster multiplicity %ld \n(x,y)=(%e,%e) \nprecision=(ex,ey)=(%e,%e) \ntimestamp = %lf \n\n",precluster.size(),cluster.getx(),cluster.gety(),cluster.getex(),cluster.getey(),cluster.gettimestamp());

    return cluster;
}

//_________________________________________________________________________________________________

//_________________________________________________________________________________________________
Clustering::Cluster Clustering::FinderGaussianFit(std::vector<Digit> &precluster, Clustering::Cluster clustertmpCOG)
{  //Fit Gaussien Manuel - voir ComputePositionGaussianFit pour le fit Minuit
    cout << "\n\n==========\nRunning GaussianFit Algorithm\n\n" << endl;
    
    Cluster cluster;

    //Création des vecteurs position du hit sur chaque cqthode
    Double_t xhit[2] = { clustertmpCOG.getx(), clustertmpCOG.getx() };
    Double_t yhit[2] = { clustertmpCOG.gety(), clustertmpCOG.gety() };
    
  //  cout << "clustertmpCOG.getx()" << clustertmpCOG.getx() << endl;
    
    //Création du float porteur de la coordonnée la plus précise
    Double_t xhitfinal = 0;
    Double_t yhitfinal = 0;
    
    //Création des vecteurs de charge totale, multiplicité et tailles sur chaque cathode
    Double_t chargetot[2] = { 0.0, 0.0 };
    Double_t multiplicity[2] = { 0.0, 0.0 };
    Double_t xsize[2] = { 0.0, 0.0 };
    Double_t ysize[2] = { 0.0, 0.0 };
    
    //Détermination des vecteurs de charge totale, multiplicité et taille qui devient précision
    for ( Int_t i = 0; i < precluster.size(); ++i ){
        int detid = precluster[i].getDetID();
        int padid = precluster[i].getPadID();
        mapping::Segmentation pad(detid); // = mapping::Segmentation(detid);
            
        for ( int cathode = 0; cathode < 2; ++cathode ) //On boucle sur les deux plans de cathodes
        {
                
            if ( cathode == pad.isBendingPad(padid) ) //On regarde à quelle cathode appartient le pad
                    {
                        std::vector<Double_t> padPosition = {pad.padPositionX(padid), pad.padPositionY(padid)};
                        std::vector<Double_t> padSize = {pad.padSizeX(padid), pad.padSizeY(padid)};

                        chargetot[cathode] += precluster[i].getADC();
                        multiplicity[cathode] += 1;
                        xsize[cathode] += padSize[0];
                        ysize[cathode] += padSize[1];
                    }
        }
    }
    for ( int cathode = 0; cathode < 2; ++cathode ) //On boucle sur les deux plans de cathodes
    {
        if ( multiplicity[cathode] != 0 )
        {
          xsize[cathode] /= multiplicity[cathode];
          ysize[cathode] /= multiplicity[cathode];
        }
    }
    
    //Minimisation du Chi2 en 2D sur chaque cathode
    for ( int cathode = 0; cathode < 2; ++cathode ) //On boucle sur les deux plans de cathodes
    {
    
        Double_t sigx = 0.7604;
        Double_t sigy = 0.7808;
        
        Double_t pasxhit = 0.00001;
        Double_t pasyhit = 0.00001;
        
        int compteurboucle = 0;
        
        Double_t ChiDeux = Chi2Gauss(cathode, chargetot[cathode], precluster, clustertmpCOG, xhit[cathode], yhit[cathode], sigx, sigy);
        Double_t OldChiDeux = 1;
        Double_t NewChiDeux = 100;
        Double_t deltaxhit = 100;
        Double_t deltayhit = 100;
        Double_t limitstrict = 0.000001;
        Double_t limitlean = 0.0001;
        Double_t limitx = 0.;
        Double_t limity = 0.;
        
        if(cathode == 0){
            limitx = limitstrict;
            limity = limitlean;
        }
        else{
            limity = limitstrict;
            limitx = limitlean;
        }
        
        while ((sqrt(pow((NewChiDeux-OldChiDeux)/OldChiDeux,2)) > 0.001)||(abs(deltaxhit)>limitx)||(abs(deltayhit)>limity)){
            compteurboucle++;
            
            //Iteration and improvement of estimate
            Double_t dChidxhit = (Chi2Gauss(cathode, chargetot[cathode], precluster, clustertmpCOG, xhit[cathode]+pasxhit, yhit[cathode], sigx, sigy)-Chi2Gauss(cathode, chargetot[cathode], precluster, clustertmpCOG, xhit[cathode], yhit[cathode], sigx, sigy))/pasxhit;
            Double_t dChidyhit = (Chi2Gauss(cathode, chargetot[cathode], precluster, clustertmpCOG, xhit[cathode], yhit[cathode]+pasyhit, sigx, sigy)-Chi2Gauss(cathode, chargetot[cathode], precluster, clustertmpCOG, xhit[cathode], yhit[cathode], sigx, sigy))/pasyhit;
            

            deltaxhit = pasxhit*dChidxhit*1;
            deltayhit = pasyhit*dChidyhit*1;
            xhit[cathode] = xhit[cathode] - pasxhit*dChidxhit*1;
            yhit[cathode] = yhit[cathode] - pasyhit*dChidyhit*1;
            
            
            // Update ChiDeux
            OldChiDeux = ChiDeux;
            ChiDeux = Chi2Gauss(cathode, chargetot[cathode], precluster, clustertmpCOG, xhit[cathode], yhit[cathode], sigx, sigy);
            NewChiDeux = ChiDeux;
            cout << "\nBoucle While numero " << compteurboucle << endl;
            cout << "xhit: " << xhit[cathode] << endl;
            cout << "yhit: " << yhit[cathode] << endl;
            cout << "deltaxhit: " << deltaxhit << endl;
            cout << "deltayhit: " << deltayhit << endl;
            cout << "OldChiDeux: " << OldChiDeux << endl;
            cout << "NewChiDeux: " << NewChiDeux << endl;
            
            if(compteurboucle > 10000){
                cout << "Déjà 10.000 boucles pour trouver Chi2, utilisation résultats COG..." << endl;
                return clustertmpCOG;
            }
        }
        
        cout << "\n\ncathode " << cathode << endl;
        cout << "Boucle Finale numero " << compteurboucle << endl;
        cout << "xhit final: " << xhit[cathode] << endl;
        cout << "yhit final: " << yhit[cathode] << endl;
        cout << "sigx: " << sigx << endl;
        cout << "sigy: " << sigy << endl;
        cout << "Chi2 Final: " << NewChiDeux << endl;
    
    }
    
    
    //Détermination de la position finale du hit à l'aide de la résolution sur chaque cathode
    xhitfinal = xhit[0];
    yhitfinal = yhit[1];
    Double_t ex = xsize[0];
    Double_t ey = ysize[1];
    
    double timestamp = precluster[0].getTimeStamp();  // On récupère le timestamp du premier digit du precluster
    
    
    //On remplit un cluster avec les infos nécessaires
    cluster.setx(xhitfinal);
    cluster.sety(yhitfinal);
    cluster.setex(ex);
    cluster.setey(ey);
    cluster.settimestamp(timestamp);    // En première approx. on peut affecter le timestamp du pre;ier digit au cluster
    
     printf("\n\nCluster par COG:\nCluster multiplicity %ld \n(x,y)=(%e,%e) \nprecision=(ex,ey)=(%e,%e) \ntimestamp = %lf \n\n",precluster.size(),clustertmpCOG.getx(),clustertmpCOG.gety(),clustertmpCOG.getex(),clustertmpCOG.getey(),clustertmpCOG.gettimestamp());
    
            printf("\n\nCluster par gaussian fit:\nCluster multiplicity %ld \n(x,y)=(%e,%e) \nprecision=(ex,ey)=(%e,%e) \ntimestamp = %lf \n\n",precluster.size(),cluster.getx(),cluster.gety(),cluster.getex(),cluster.getey(),cluster.gettimestamp());
    
    return cluster;

}

//_________________________________________________________________________________________________
Double_t Clustering::Chi2Gauss(int cathode, Double_t chargetot, std::vector<Digit> &precluster, Clustering::Cluster clustertmp, Double_t xhit, Double_t yhit, Double_t sigx, Double_t sigy)
{
    
    Double_t sum = 0;
    for(Int_t i = 0; i < precluster.size(); ++i){
        
        int detid = precluster[i].getDetID();
        int padid = precluster[i].getPadID();
        mapping::Segmentation pad(detid); // = mapping::Segmentation(detid);
        
        if ( cathode == pad.isBendingPad(padid) ) //On regarde à quelle cathode appartient le pad
        {
              //On définit le vecteur position et taille du pad en question
            std::vector<double> padPosition = {pad.padPositionX(padid), pad.padPositionY(padid)};
            std::vector<double> padSize = {pad.padSizeX(padid), pad.padSizeY(padid)};
            
            std::vector<double> lowerLeft = {xhit - padPosition[0] - 0.5*padSize[0], yhit - padPosition[1] - 0.5*padSize[1]};
            std::vector<double> upperRight= {lowerLeft[0] + padSize[0], lowerLeft[1] + padSize[1]};
            
            Double_t qfit = chargetot*IntGaussXY(lowerLeft[0], lowerLeft[1], upperRight[0], upperRight[1], sigx, sigy);

            Double_t qtrue = precluster[i].getADC();
            
            Double_t delta = qtrue-qfit;
            
            sum += pow(delta,2)/qtrue;
        }
    }
    
    return sum;
}

//_________________________________________________________________________________________________
void FitFunctionGaussianFit(Int_t& /*notused*/, Double_t* /*notused*/,
                            Double_t& sum, Double_t* par,
                             Int_t /*notused*/)
  {
    /// Chi2 Function to minimize: Gaussian charge distribution in 2 dimensions
      
          //Méthode bourrine. Définir un vecteur de paramètres énorme:
          // x   y   cathode   chgtot   sigx   sigy   digit par digit (Det Pad ADC)

//    TObjArray* userObjects = static_cast<TObjArray*>(TVirtualFitter::GetFitter()->GetObjectFit());
//
//    std::vector<Digit> precluster = static_cast<std::vector<Digit>>(userObjects->At(0));
//    int cathode = static_cast<int>(userObjects->At(1));
//    Double_t chargetot = static_cast<Double_t>(userObjects->At(2));
//    Double_t Kx3 = static_cast<Double_t>(userObjects->At(3));
//    Double_t Ky3 = static_cast<Double_t>(userObjects->At(4));
      
      Clustering cluster;

    sum = 0.0;
      
      int cathode = par[2];
      int chargetot = par[3];
      double sigx = par[4];
      double sigy = par[5];
      int nbdigits = par[6];

    for(Int_t i = 0; i < nbdigits; ++i){

    int detid = par[7+(3*i)];
    int padid = par[8+(3*i)];
    int ADC = par[9+(3*i)];
        
    mapping::Segmentation pad(detid); // = mapping::Segmentation(detid);

    if ( cathode == pad.isBendingPad(padid) ) //On regarde à quelle cathode appartient le pad
    {
          //On définit le vecteur position et taille du pad en question
        std::vector<Double_t> padPosition = {pad.padPositionX(padid), pad.padPositionY(padid)};
        std::vector<Double_t> padSize = {pad.padSizeX(padid), pad.padSizeY(padid)};

        std::vector<Double_t> lowerLeft = {par[0] - padPosition[0] - 0.5*padSize[0], par[1] - padPosition[1] - 0.5*padSize[1]};
        std::vector<Double_t> upperRight= {lowerLeft[0] + padSize[0], lowerLeft[1] + padSize[1]};

        Double_t qfit = chargetot*cluster.IntGaussXY(lowerLeft[0], lowerLeft[1], upperRight[0], upperRight[1], sigx, sigy);

        Double_t qtrue = ADC;

        Double_t delta = qtrue-qfit;

        sum += pow(delta,2)/qtrue;
    }
  }
}

//_________________________________________________________________________________________________
Clustering::Cluster Clustering::ComputePositionGaussianFit(std::vector<Digit> &precluster, Clustering::Cluster clustertmpCOG)
{
  /// Compute the position of the given cluster, by fitting a Gaussian
  /// charge distribution to it

    int nparams = 6+(3*precluster.size());
    Double_t vecxmin[2] = {1E3, 1E3};
    Double_t vecxmax[2] = {-1E3, -1E3};
    Double_t vecymin[2] = {1E3, 1E3};
    Double_t vecymax[2] = {-1E3, -1E3};
    
  TVirtualFitter* fitter = TVirtualFitter::Fitter(0,nparams);
  fitter->SetFCN(FitFunctionGaussianFit);

  cout << "\n\n==========\nRunning MinuitGaussian Algorithm.\n\n" << endl;

  Cluster cluster;

  //Création des vecteurs position du hit sur chaque cqthode
  Double_t xhit[2] = { clustertmpCOG.getx(), clustertmpCOG.getx() };
  Double_t yhit[2] = { clustertmpCOG.gety(), clustertmpCOG.gety() };

    //Création des vecteurs de charge totale, multiplicité et tailles sur chaque cathode
       Double_t chargetot[2] = { 0.0, 0.0 };
       Double_t multiplicity[2] = { 0.0, 0.0 };
       Double_t xsize[2] = { 0.0, 0.0 };
       Double_t ysize[2] = { 0.0, 0.0 };

   //Détermination des vecteurs de charge totale, multiplicité et taille qui devient précision
   for ( Int_t i = 0; i < precluster.size(); ++i ){
       int detid = precluster[i].getDetID();
       int padid = precluster[i].getPadID();
       
       mapping::Segmentation pad(detid); // = mapping::Segmentation(detid);

       for ( int cathode = 0; cathode < 2; ++cathode ) //On boucle sur les deux plans de cathodes
       {

           if ( cathode == pad.isBendingPad(padid) ) //On regarde à quelle cathode appartient le pad
                   {
                       std::vector<Double_t> padPosition = {pad.padPositionX(padid), pad.padPositionY(padid)};
                       std::vector<Double_t> padSize = {pad.padSizeX(padid), pad.padSizeY(padid)};
                       
                       vecxmin[cathode] = TMath::Min(padPosition[0]-0.5*padSize[0],vecxmin[cathode]);
                       vecxmax[cathode] = TMath::Max(padPosition[0]+0.5*padSize[0],vecxmax[cathode]);
                       vecymin[cathode] = TMath::Min(padPosition[1]-0.5*padSize[1],vecymin[cathode]);
                       vecymax[cathode] = TMath::Max(padPosition[1]+0.5*padSize[1],vecymax[cathode]);

                       chargetot[cathode] += precluster[i].getADC();
                       multiplicity[cathode] += 1;
                   }
       }
   }
    

  Float_t stepX = 0.00001; // cm
  Float_t stepY = 0.00001; // cm
    Double_t sigx = 0.7604;
    Double_t sigy = 0.7808;

    Double_t arg(-1); // disable printout

    for ( int cathode = 0; cathode < 2; ++cathode ) //On boucle sur les deux plans de cathodes
    {
        fitter->Clear();
        fitter->ExecuteCommand("SET PRINT",&arg,1);
        if(vecxmin[cathode]>vecxmax[cathode]){
            vecxmin[cathode] = -40.;
            vecxmax[cathode] = +40.;
        }
        if(vecymin[cathode]>vecymax[cathode]){
            vecymin[cathode] = -20.;
            vecymax[cathode] = +20.;
        }

        fitter->SetParameter(0,"cluster X position",xhit[cathode],stepX,vecxmin[cathode],vecxmax[cathode]);
        fitter->SetParameter(1,"cluster Y position",yhit[cathode],stepY,vecymin[cathode],vecymax[cathode]);
        fitter->SetParameter(2,"cathode",cathode,0,0,0);
        fitter->SetParameter(3,"chargetot",chargetot[cathode],0,0,0);
        fitter->SetParameter(4,"sigx",sigx,0,0,0);
        fitter->SetParameter(5,"sigy",sigy,0,0,0);
        fitter->SetParameter(6,"NbDigits",precluster.size(),0,0,0);
        
        for(int j=0; j<precluster.size(); j++){
            fitter->SetParameter(7+(3*j),"DetID",precluster[j].getDetID(),0,0,0);
            fitter->SetParameter(8+(3*j),"PadID",precluster[j].getPadID(),0,0,0);
            fitter->SetParameter(9+(3*j),"ADC",precluster[j].getADC(),0,0,0);
        }

        Double_t stratArg(2);
        fitter->ExecuteCommand("SET STR",&stratArg,1);
        Int_t val = fitter->ExecuteCommand("MIGRAD",0,0);
        if ( val && chargetot[cathode] != 0)
        {
           //Fit failed with robust strategy 2 try balanced strategy 1
            Double_t stratArg(1);
            fitter->ExecuteCommand("SET STR",&stratArg,1);
            Int_t val2 = fitter->ExecuteCommand("MIGRAD",0,0);
            if ( val2 ){
                // fit failed. Using COG results, with big errors
                cout << "Fit failed on viable cathode. Using COG results for cluster." << endl;
                    return clustertmpCOG;
            }
        }

        xhit[cathode] = fitter->GetParameter(0);
        yhit[cathode] = fitter->GetParameter(1);
        xsize[cathode] = fitter->GetParError(0);
        ysize[cathode] = fitter->GetParError(1);

        Double_t amin, edm, errdef;
        Int_t nvpar, nparx;

        fitter->GetStats(amin, edm, errdef, nvpar, nparx);

        Double_t chi2 = amin;

        printf("Cluster fitted to (x,y)=(%e,%e) (xerr,yerr)=(%e,%e) \n chi2=%e ndf=%d",
                        xhit[cathode],yhit[cathode],
                        xsize[cathode],ysize[cathode],chi2,fitter->GetNumberFreeParameters());
    }

    //Détermination de la position finale du hit à l'aide de la résolution sur chaque cathode
    Double_t xhitfinal = xhit[0];
    Double_t yhitfinal = yhit[1];
    Double_t ex = xsize[0];
    Double_t ey = ysize[1];

    double timestamp = precluster[0].getTimeStamp();  // On récupère le timestamp du premier digit du precluster
    
     printf("\n\nCluster par gaussian fit:\nCluster multiplicity %ld \n(x,y)=(%e,%e) \nprecision=(ex,ey)=(%e,%e) \ntimestamp = %lf \n\n",precluster.size(),xhitfinal,yhitfinal,ex,ey,timestamp);

    //On remplit un cluster avec les infos nécessaires
    cluster.setx(xhitfinal);
    cluster.sety(yhitfinal);
    cluster.setex(ex);
    cluster.setey(ey);
    cluster.settimestamp(timestamp);    // En première approx. on peut affecter le timestamp du pre;ier digit au cluster

    return cluster;
}


//_________________________________________________________________________________________________
void FitFunctionDoubleGaussianFit(Int_t& /*notused*/, Double_t* /*notused*/,
                            Double_t& sum, Double_t* par,
                             Int_t /*notused*/)
  {
    /// Chi2 Function to minimize: Mathieson charge distribution in 2 dimensions
      
          //Méthode bourrine. Définir un vecteur de paramètres énorme:
          // x   y   cathode   chgtot   sig1s   sig2s   chgfracs   digit par digit (Det Pad ADC)

//    TObjArray* userObjects = static_cast<TObjArray*>(TVirtualFitter::GetFitter()->GetObjectFit());
//
//    std::vector<Digit> precluster = static_cast<std::vector<Digit>>(userObjects->At(0));
//    int cathode = static_cast<int>(userObjects->At(1));
//    Double_t chargetot = static_cast<Double_t>(userObjects->At(2));
//    Double_t Kx3 = static_cast<Double_t>(userObjects->At(3));
//    Double_t Ky3 = static_cast<Double_t>(userObjects->At(4));
      
      Clustering cluster;

    sum = 0.0;
      
      int cathode = par[2];
      int chargetot = par[3];
      double sig1x = par[4];
      double sig1y = par[5];
      double sig2x = par[6];
      double sig2y = par[7];
      double chgfracx = par[8];
      double chgfracy = par[9];
      int nbdigits = par[10];

    for(Int_t i = 0; i < nbdigits; ++i){

    int detid = par[11+(3*i)];
    int padid = par[12+(3*i)];
    int ADC = par[13+(3*i)];
        
    mapping::Segmentation pad(detid); // = mapping::Segmentation(detid);

    if ( cathode == pad.isBendingPad(padid) ) //On regarde à quelle cathode appartient le pad
    {
          //On définit le vecteur position et taille du pad en question
        std::vector<Double_t> padPosition = {pad.padPositionX(padid), pad.padPositionY(padid)};
        std::vector<Double_t> padSize = {pad.padSizeX(padid), pad.padSizeY(padid)};

        std::vector<Double_t> lowerLeft = {par[0] - padPosition[0] - 0.5*padSize[0], par[1] - padPosition[1] - 0.5*padSize[1]};
        std::vector<Double_t> upperRight= {lowerLeft[0] + padSize[0], lowerLeft[1] + padSize[1]};

        Double_t qfit = chargetot*cluster.IntDoubleGaussXY(lowerLeft[0], lowerLeft[1], upperRight[0], upperRight[1], sig1x, sig1y, sig2x, sig2y, chgfracx, chgfracy);

        Double_t qtrue = ADC;

        Double_t delta = qtrue-qfit;

        sum += pow(delta,2)/qtrue;
    }
  }
}

//_________________________________________________________________________________________________
Clustering::Cluster Clustering::ComputePositionDoubleGaussianFit(std::vector<Digit> &precluster, Clustering::Cluster clustertmpCOG)
{
  /// Compute the position of the given cluster, by fitting a Gaussian
  /// charge distribution to it

    int nparams = 10+(3*precluster.size());
    Double_t vecxmin[2] = {1E3, 1E3};
    Double_t vecxmax[2] = {-1E3, -1E3};
    Double_t vecymin[2] = {1E3, 1E3};
    Double_t vecymax[2] = {-1E3, -1E3};
    
  TVirtualFitter* fitter = TVirtualFitter::Fitter(0,nparams);
  fitter->SetFCN(FitFunctionDoubleGaussianFit);

  cout << "\n\n==========\nRunning MinuitDoubleGaussian Algorithm.\n\n" << endl;

  Cluster cluster;

  //Création des vecteurs position du hit sur chaque cqthode
  Double_t xhit[2] = { clustertmpCOG.getx(), clustertmpCOG.getx() };
  Double_t yhit[2] = { clustertmpCOG.gety(), clustertmpCOG.gety() };

    //Création des vecteurs de charge totale, multiplicité et tailles sur chaque cathode
       Double_t chargetot[2] = { 0.0, 0.0 };
       Double_t multiplicity[2] = { 0.0, 0.0 };
       Double_t xsize[2] = { 0.0, 0.0 };
       Double_t ysize[2] = { 0.0, 0.0 };

   //Détermination des vecteurs de charge totale, multiplicité et taille qui devient précision
   for ( Int_t i = 0; i < precluster.size(); ++i ){
       int detid = precluster[i].getDetID();
       int padid = precluster[i].getPadID();
       
       mapping::Segmentation pad(detid); // = mapping::Segmentation(detid);

       for ( int cathode = 0; cathode < 2; ++cathode ) //On boucle sur les deux plans de cathodes
       {

           if ( cathode == pad.isBendingPad(padid) ) //On regarde à quelle cathode appartient le pad
                   {
                       std::vector<Double_t> padPosition = {pad.padPositionX(padid), pad.padPositionY(padid)};
                       std::vector<Double_t> padSize = {pad.padSizeX(padid), pad.padSizeY(padid)};
                       
                       vecxmin[cathode] = TMath::Min(padPosition[0]-0.5*padSize[0],vecxmin[cathode]);
                       vecxmax[cathode] = TMath::Max(padPosition[0]+0.5*padSize[0],vecxmax[cathode]);
                       vecymin[cathode] = TMath::Min(padPosition[1]-0.5*padSize[1],vecymin[cathode]);
                       vecymax[cathode] = TMath::Max(padPosition[1]+0.5*padSize[1],vecymax[cathode]);

                       chargetot[cathode] += precluster[i].getADC();
                       multiplicity[cathode] += 1;
                   }
       }
   }
    

  Float_t stepX = 0.00001; // cm
  Float_t stepY = 0.00001; // cm
    Double_t sig1x = 1.1059;
    Double_t sig1y = 1.1350;
    Double_t sig2x = 0.5702;
    Double_t sig2y = 0.5785;
    Double_t chgfracx = 0.6233;
    Double_t chgfracy = 0.6519;

    Double_t arg(-1); // disable printout

    for ( int cathode = 0; cathode < 2; ++cathode ) //On boucle sur les deux plans de cathodes
    {
        fitter->Clear();
        fitter->ExecuteCommand("SET PRINT",&arg,1);
        if(vecxmin[cathode]>vecxmax[cathode]){
            vecxmin[cathode] = -40.;
            vecxmax[cathode] = +40.;
        }
        if(vecymin[cathode]>vecymax[cathode]){
            vecymin[cathode] = -20.;
            vecymax[cathode] = +20.;
        }

        fitter->SetParameter(0,"cluster X position",xhit[cathode],stepX,vecxmin[cathode],vecxmax[cathode]);
        fitter->SetParameter(1,"cluster Y position",yhit[cathode],stepY,vecymin[cathode],vecymax[cathode]);
        fitter->SetParameter(2,"cathode",cathode,0,0,0);
        fitter->SetParameter(3,"chargetot",chargetot[cathode],0,0,0);
        fitter->SetParameter(4,"sig1x",sig1x,0,0,0);
        fitter->SetParameter(5,"sig1y",sig1y,0,0,0);
        fitter->SetParameter(6,"sig2x",sig2x,0,0,0);
        fitter->SetParameter(7,"sig2y",sig2y,0,0,0);
        fitter->SetParameter(8,"chgfracx",chgfracx,0,0,0);
        fitter->SetParameter(9,"chgfracy",chgfracy,0,0,0);
        fitter->SetParameter(10,"NbDigits",precluster.size(),0,0,0);
        
        for(int j=0; j<precluster.size(); j++){
            fitter->SetParameter(11+(3*j),"DetID",precluster[j].getDetID(),0,0,0);
            fitter->SetParameter(12+(3*j),"PadID",precluster[j].getPadID(),0,0,0);
            fitter->SetParameter(13+(3*j),"ADC",precluster[j].getADC(),0,0,0);
        }

        Double_t stratArg(2);
        fitter->ExecuteCommand("SET STR",&stratArg,1);
        Int_t val = fitter->ExecuteCommand("MIGRAD",0,0);
        if ( val && chargetot[cathode] != 0)
        {
            //Fit failed with robust strategy 2 try balanced strategy 1
            Double_t stratArg(1);
            fitter->ExecuteCommand("SET STR",&stratArg,1);
            Int_t val2 = fitter->ExecuteCommand("MIGRAD",0,0);
            if ( val2 ){
                // fit failed. Using COG results, with big errors
                cout << "Fit failed on viable cathode. Using COG results for cluster." << endl;
                    return clustertmpCOG;
            }
        }

        xhit[cathode] = fitter->GetParameter(0);
        yhit[cathode] = fitter->GetParameter(1);
        xsize[cathode] = fitter->GetParError(0);
        ysize[cathode] = fitter->GetParError(1);

        Double_t amin, edm, errdef;
        Int_t nvpar, nparx;

        fitter->GetStats(amin, edm, errdef, nvpar, nparx);

        Double_t chi2 = amin;

        printf("Cluster fitted to (x,y)=(%e,%e) (xerr,yerr)=(%e,%e) \n chi2=%e ndf=%d",
                        xhit[cathode],yhit[cathode],
                        xsize[cathode],ysize[cathode],chi2,fitter->GetNumberFreeParameters());
    }

    //Détermination de la position finale du hit à l'aide de la résolution sur chaque cathode
    Double_t xhitfinal = xhit[0];
    Double_t yhitfinal = yhit[1];
    Double_t ex = xsize[0];
    Double_t ey = ysize[1];

    double timestamp = precluster[0].getTimeStamp();  // On récupère le timestamp du premier digit du precluster
    
     printf("\n\nCluster par double gaussian fit:\nCluster multiplicity %ld \n(x,y)=(%e,%e) \nprecision=(ex,ey)=(%e,%e) \ntimestamp = %lf \n\n",precluster.size(),xhitfinal,yhitfinal,ex,ey,timestamp);

    //On remplit un cluster avec les infos nécessaires
    cluster.setx(xhitfinal);
    cluster.sety(yhitfinal);
    cluster.setex(ex);
    cluster.setey(ey);
    cluster.settimestamp(timestamp);    // En première approx. on peut affecter le timestamp du pre;ier digit au cluster

    return cluster;
}

/*

//_________________________________________________________________________________________________
void Clustering::ComparisonMathiesonGauss(std::string& buffer)
{
    Cluster* clusterCOG;
    Cluster* clusterMathieson;
    Cluster* clusterGauss;
    
    clusterCOG = runFinderCOG(buffer);
    clusterMathieson = runFinderSimpleFit(buffer);
    clusterGauss = runFinderGaussianFit(buffer);

}

*/

} // namespace mch
} // namespace o2
