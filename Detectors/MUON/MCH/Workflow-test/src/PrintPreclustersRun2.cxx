// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include <fstream>
#include <chrono>
#include <TH2.h>
#include <TFile.h>

#include "DigitsFileReader.h"
#include "MCHPreClustering/PreClusterFinder.h"
#include "MCHMappingInterface/Segmentation.h"
#include "MCHClustering/ClusteringForTest.h"

using namespace o2::mch;
using namespace std;

int iClus;

#define DEID 919

TH2F* plotPrecuster(PreCluster& preCluster, gsl::span<Digit> digits, std::vector<Clustering::Cluster>& clusters)
{
  int maxPadID = -1;
  int maxADC = 0;
  for (const auto& digit : digits) {
    if(digit.getADC() > maxADC && digit.getPadID() < 1024) {
      maxADC = digit.getADC();
      maxPadID = digit.getPadID();
    }
  }

  if(maxPadID < 0) return nullptr;
  const mapping::Segmentation& segment = mapping::segmentation(digits[0].getDetID());

  double padMaxPosition[2] = {segment.padPositionX(maxPadID), segment.padPositionY(maxPadID)};
  double padMaxSize[2] = {segment.padSizeX(maxPadID), segment.padSizeY(maxPadID)};

  TH2F* h2 = new TH2F(TString::Format("h2_%d",iClus),
      "Pad amplitude", 11, padMaxPosition[0]-padMaxSize[0]*5-padMaxSize[0]/2, padMaxPosition[0]+padMaxSize[0]*5+padMaxSize[0]/2,
      11, padMaxPosition[1]-padMaxSize[1]*5-padMaxSize[1]/2, padMaxPosition[1]+padMaxSize[1]*5+padMaxSize[1]/2);
  h2->SetDrawOption("colz");

  for (const auto& digit : digits) {
    if(digit.getPadID() < 1024) {
      float pX = segment.padPositionX(digit.getPadID());
      float pY = segment.padPositionY(digit.getPadID());
      int bx = h2->GetXaxis()->FindBin(pX);
      int by = h2->GetYaxis()->FindBin(pY);
      h2->SetBinContent(bx, by, digit.getADC());
      std::cout<<"pX="<<pX<<"  pY="<<pY<<"  bx="<<bx<<"  by="<<by<<"  ADC="<<digit.getADC()<<endl;
    }
  }
  return h2;
}

int main(int argc, char** argv)
{
  std::chrono::duration<double, std::milli> mTimeDigitsRead;
  std::chrono::duration<double, std::milli> mTimeDigitsStore;
  std::chrono::duration<double, std::milli> mTimeDigitsLoad;
  std::chrono::duration<double, std::milli> mTimePreClustering;

  DigitsFileReader digitsReader(argv[1]);
  ofstream outFile(argv[2],ios::out|ios::binary);

  PreClusterFinder preClusterFinder;
  preClusterFinder.init();

  TFile fHits("hits.root","RECREATE");

  TH1F* hPadCharge = new TH1F("hPadCharge", "Pad charge distribution - DE 919", 420, 0, 4200);

  iClus = 0;

  Clustering clustering;

  Digit* digitsBuffer = NULL;
  int nevent = 1;

  // load digits from binary input file, block-by-block
  while(true) {

    auto tStart = std::chrono::high_resolution_clock::now();
    auto ok = digitsReader.readDigitsFromFile();
    if(!ok) break;
    auto tEnd = std::chrono::high_resolution_clock::now();
    mTimeDigitsRead += tEnd - tStart;

    // get number of loaded digits and store them into a memory buffer
    auto nDigits = digitsReader.getNumberOfDigits();
    //printf("nEvent: %d  nDigits: %d\n", nevent, (int)nDigits);
    //continue;
    digitsBuffer = (Digit*)realloc(digitsBuffer, sizeof(Digit) * nDigits);
    tStart = std::chrono::high_resolution_clock::now();
    digitsReader.storeDigits(digitsBuffer);
    tEnd = std::chrono::high_resolution_clock::now();
    mTimeDigitsStore += tEnd - tStart;

    gsl::span<Digit> digits{digitsBuffer, nDigits};

    printf("\n\n\n====================\n\n");

    for (const auto& digit : digits) {
      if(/*digit.getPadID() < 1024 &&*/ digit.getDetID() == DEID) {
        const mapping::Segmentation& segment = mapping::segmentation(digit.getDetID());
        float pX = segment.padPositionX(digit.getPadID());
        float pY = segment.padPositionY(digit.getPadID());
        std::cout<<digit.getDetID()<<"  pX="<<pX<<"  pY="<<pY<<"  ADC="<<digit.getADC()<<endl;
        if(digit.getPadID() < 1024) hPadCharge->Fill(digit.getADC());
      }
    }

    // load the digits from the memory buffer and run the pre-clustering phase
    preClusterFinder.reset();
    tStart = std::chrono::high_resolution_clock::now();
    preClusterFinder.loadDigits({digitsBuffer, nDigits});
    tEnd = std::chrono::high_resolution_clock::now();
    mTimeDigitsLoad += tEnd - tStart;

    tStart = std::chrono::high_resolution_clock::now();
    int nPreClusters = preClusterFinder.run();
    tEnd = std::chrono::high_resolution_clock::now();
    mTimePreClustering += tEnd - tStart;

    std::vector<PreCluster> preClusters{}; ///< vector of preclusters
    std::vector<Digit> usedDigits{};       ///< vector of digits in the preclusters
    preClusters.reserve(nPreClusters); // to avoid reallocation if
    usedDigits.reserve(nDigits); // the capacity is exceeded
    preClusterFinder.getPreClusters(preClusters, usedDigits);
    if (usedDigits.size() != nDigits) {
      throw runtime_error("some digits have been lost during the preclustering");
    }



    //cout<<preClusterFinder<<std::endl;

    // write the number of preclusters
    outFile.write(reinterpret_cast<char*>(&nPreClusters), sizeof(int));

    for (auto& preCluster : preClusters) {

      auto preClusterDigits = gsl::span<Digit>(usedDigits).subspan(preCluster.firstDigit, preCluster.nDigits);
      //std::vector<Digit> dv(preClusterDigits.begin(), preClusterDigits.end());

      gsl::span<const PreCluster> sPreCluster = {&preCluster, 1};
      std::vector<Clustering::Cluster> clusters(0);

      //cout<<"DetID: "<<preClusterDigits[0].getDetID()<<endl;
      if(preClusterDigits[0].getDetID() == DEID && preClusterDigits.size()>1) {
        preCluster.print(cout, usedDigits);
        // Fit Mathieson
        //clustering.runFinderSimpleFit(sPreCluster, usedDigits, clusters);
        TH2F* h2 = plotPrecuster(preCluster, preClusterDigits, clusters);
        if(h2) {
          fHits.cd();
          h2->Write();
          TH1D* h1 = h2->ProjectionY();
          h1->Write();
        }
        iClus += 1;
        //getchar();
      }
      float clusPosition[2] = {
          clusters.empty() ? 0 : static_cast<float>(clusters[0].getx()),
              clusters.empty() ? 0 : static_cast<float>(clusters[0].gety())
      };
      outFile.write(reinterpret_cast<char*>(clusPosition), sizeof(clusPosition));

      // write the total number of digits in this precluster
      int nDigits = preClusterDigits.size();
      outFile.write(reinterpret_cast<char*>(&nDigits), sizeof(int));


      for (auto& digit : preClusterDigits) {
        int detid = digit.getDetID();
        int padid = digit.getPadID();

        if(false){
          cout << "\nDetID:" << detid << " PadID:" << padid << endl;
        }
        //mapping::Segmentation pad(detid); // = mapping::Segmentation(detid);
        const mapping::Segmentation& pad = mapping::segmentation(detid);

        int digitInfo[4] = {detid, padid, static_cast<int>(digit.getTimeStamp()), static_cast<int>(digit.getADC())};
        //On définit le vecteur position et taille du pad en question
        float padPosition[2] = {static_cast<float>(pad.padPositionX(padid)), static_cast<float>(pad.padPositionY(padid))};
        float padSize[2] = {static_cast<float>(pad.padSizeX(padid)), static_cast<float>(pad.padSizeY(padid))};

        outFile.write(reinterpret_cast<char*>(digitInfo), sizeof(digitInfo));
        outFile.write(reinterpret_cast<char*>(padPosition), sizeof(padPosition));
        outFile.write(reinterpret_cast<char*>(padSize), sizeof(padSize));
      }
    }


    nevent += 1;
    if( (nevent%10) == 0) cout<<nevent<<endl;
    if( nevent > 10000 ) break;
  }

  fHits.cd();
  hPadCharge->Write();

  std::cout<<nevent<<"  "<<mTimeDigitsRead.count()<<" ms  "<<mTimeDigitsStore.count()<<" ms  "<<mTimeDigitsLoad.count()<<" ms  "<<mTimePreClustering.count()<<" ms"<<std::endl;

  return 0;
}
