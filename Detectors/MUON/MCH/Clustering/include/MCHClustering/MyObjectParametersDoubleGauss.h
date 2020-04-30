#ifndef ALICEO2_MCH_MYOBJECTPARAMETERSDOUBLEGAUSS_H_
#define ALICEO2_MCH_MYOBJECTPARAMETERSDOUBLEGAUSS_H_

#include <TROOT.h>
#include <iostream>
#include <vector>
#include <TObject.h>

namespace o2
{
namespace mch
{

using namespace std;

class MyObjectParametersDoubleGauss : public TObject
{
 public:
    MyObjectParametersDoubleGauss();
  //MyObjectParameters(int cathode, double chargetot, double Ky3, double Kx3, int Nbdigits);
    virtual ~MyObjectParametersDoubleGauss();

  int getfCathode() const { return fCathode; }
  void setfCathode(int cathode) { fCathode = cathode; }
    double getfChargetot() const { return fChargetot; }
    void setfChargetot(double chargetot) { fChargetot = chargetot; }
    double getfSig1y() const { return fSig1y; }
    void setfSig1y(double sig1y) { fSig1y = sig1y; }
    double getfSig1x() const { return fSig1x; }
    void setfSig1x(double sig1x) { fSig1x = sig1x; }
    double getfSig2y() const { return fSig2y; }
    void setfSig2y(double sig2y) { fSig2y = sig2y; }
    double getfSig2x() const { return fSig2x; }
    void setfSig2x(double sig2x) { fSig2x = sig2x; }
    double getfChgfracx() const { return fChgfracx; }
    void setfChgfracx(double chgfracx) { fChgfracx = chgfracx; }
    double getfChgfracy() const { return fChgfracy; }
    void setfChgfracy(double chgfracy) { fChgfracy = chgfracy; }
    int getfNbdigits() const { return fNbdigits; }
    void setfNbdigits(int Nbdigits) { fNbdigits = Nbdigits; }
    
    ClassDef(MyObjectParametersDoubleGauss, 1);

 private:
    int fCathode;
    double fChargetot;
    double fSig1y;
    double fSig1x;
    double fSig2y;
    double fSig2x;
    double fChgfracx;
    double fChgfracy;
    int fNbdigits;

};

}
}

#endif // ALICEO2_MCH_MYOBJECTPARAMETERSDOUBLEGAUSS_H_

