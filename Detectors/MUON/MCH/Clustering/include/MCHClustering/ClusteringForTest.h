//
//  ClusteringForTest.h
//  
//
//  Created by PERRIN Sébastien on 31/10/2019.
//

#ifndef ClusteringForTest_h
#define ClusteringForTest_h

#include <stdio.h>
#include <cassert>
#include <cstdint>
#include <unordered_map>
#include <vector>

#include <gsl/span>

//include michael
//include segmentation.h

#include "MCHBase/Digit.h"
#include "MCHMappingInterface/Segmentation.h"
#include "MCHBase/PreCluster.h"

namespace o2
{
namespace mch
{

// Classes Cluster

class Clustering
{
    public:
    struct Cluster {
    Double_t x;            // x position du cluster
    Double_t y;            // y position du cluster
    // float z;            // z position du cluster
    Double_t ex;           // resolution en x
    Double_t ey;           // resolution en y
    //  float ez;           // resolution en z
    Double_t timestamp;    // Timestamp associé au cluster
    //  int idPreCluster;   // Id du PreCluster associé à ce cluster
        
        
    Double_t getx() const { return mx; }
    void setx(Double_t x) { mx = x; }
        
    Double_t gety() const { return my; }
    void sety(Double_t y) { my = y; }
        
    Double_t getex() const { return mex; }
    void setex(Double_t ex) { mex = ex; }
        
    Double_t getey() const { return mey; }
    void setey(Double_t ey) { mey = ey; }

    Double_t gettimestamp() const { return mtimestamp; }
    void settimestamp(Double_t timestamp) { mtimestamp = timestamp; }

    private:
    Double_t mx;
    Double_t my;
    Double_t mex;
    Double_t mey;
    Double_t mtimestamp;
    };
    
    //   struct PreClusterStruct {
    //       std::vector<Digit> digits;
    //   };
    //
    //   struct PreClusterSizes {
    //       std::vector<int> sizes;
    //   };

  Clustering() = default;
  ~Clustering() = default;

  Clustering(const Clustering&) = delete;
  Clustering& operator=(const Clustering&) = delete;
  Clustering(Clustering&&) = delete;
  Clustering& operator=(Clustering&&) = delete;

    
    void runFinderCOG(gsl::span<const PreCluster> preClusters, gsl::span<const Digit> digits, std::vector<Cluster>& clusters);
    Clustering::Cluster FinderCOG(std::vector<Digit> &precluster);

    
    void runFinderSimpleFit(gsl::span<const PreCluster> preClusters, gsl::span<const Digit> digits, std::vector<Cluster>& clusters);
    Clustering::Cluster ComputePositionClean(std::vector<Digit> &precluster, Clustering::Cluster clustertmpCOG);
    Double_t IntMathiesonXY(Double_t x1, Double_t y1, Double_t x2, Double_t y2, Double_t Kx3, Double_t Ky3);
    
    
    void runFinderGaussianFit(gsl::span<const PreCluster> preClusters, gsl::span<const Digit> digits, std::vector<Cluster>& clusters);
    Clustering::Cluster ComputePositionGaussianFitClean(std::vector<Digit> &precluster, Clustering::Cluster clustertmpCOG);
    Double_t IntGaussXY(Double_t x1, Double_t y1, Double_t x2, Double_t y2, Double_t sigx, Double_t sigy);
    
    
    void runFinderDoubleGaussianFit(gsl::span<const PreCluster> preClusters, gsl::span<const Digit> digits, std::vector<Cluster>& clusters);
    Clustering::Cluster ComputePositionDoubleGaussianFitClean(std::vector<Digit> &precluster, Clustering::Cluster clustertmpCOG);
    Double_t IntDoubleGaussXY(Double_t x1, Double_t y1, Double_t x2, Double_t y2, Double_t sig1x, Double_t sig1y, Double_t sig2x, Double_t sig2y, Double_t chgfracx, Double_t chgfracy);

    
//    void ComparisonMathiesonGauss(std::string& buffer);
    
    
    
    // FONCTIONS OBSOLETES
    Clustering::Cluster FinderSimpleFit(std::vector<Digit> &precluster, Clustering::Cluster clustertmpCOG);
    Double_t Chi2Mathieson(int cathode, Double_t chargetot, std::vector<Digit> &precluster, Clustering::Cluster clustertmp, Double_t xhit, Double_t yhit, Double_t Kx3, Double_t Ky3);
    Clustering::Cluster ComputePosition(std::vector<Digit> &precluster, Clustering::Cluster clustertmpCOG);
    Clustering::Cluster FinderGaussianFit(std::vector<Digit> &precluster, Clustering::Cluster clustertmpCOG);
    Double_t Chi2Gauss(int cathode, Double_t chargetot, std::vector<Digit> &precluster, Clustering::Cluster clustertmp, Double_t xhit, Double_t yhit, Double_t sigx, Double_t sigy);
    Clustering::Cluster ComputePositionGaussianFit(std::vector<Digit> &precluster, Clustering::Cluster clustertmpCOG);
    
    Clustering::Cluster ComputePositionDoubleGaussianFit(std::vector<Digit> &precluster, Clustering::Cluster clustertmpCOG);
    
};

void FitFunctionClean(Int_t& /*notused*/, Double_t* /*notused*/,
Double_t& sum, Double_t* par,
 Int_t /*notused*/);
void FitFunctionGaussianFitClean(Int_t& /*notused*/, Double_t* /*notused*/,
Double_t& sum, Double_t* par,
 Int_t /*notused*/);
void FitFunctionDoubleGaussianFitClean(Int_t& /*notused*/, Double_t* /*notused*/,
Double_t& sum, Double_t* par,
Int_t /*notused*/);

//FONCTIONS OBSOLETES
void FitFunction(Int_t& /*notused*/, Double_t* /*notused*/,
Double_t& sum, Double_t* par,
 Int_t /*notused*/);
void FitFunctionGaussianFit(Int_t& /*notused*/, Double_t* /*notused*/,
Double_t& sum, Double_t* par,
 Int_t /*notused*/);
void FitFunctionDoubleGaussianFit(Int_t& /*notused*/, Double_t* /*notused*/,
Double_t& sum, Double_t* par,
Int_t /*notused*/);

} // namespace mch
} // namespace o2

#endif /* ClusteringForTest_h */
