#ifndef ALICEO2_MCH_MYOBJECTPARAMETERSGAUSS_H_
#define ALICEO2_MCH_MYOBJECTPARAMETERSGAUSS_H_

#include <TROOT.h>
#include <iostream>
#include <vector>
#include <TObject.h>

namespace o2
{
namespace mch
{

using namespace std;

class MyObjectParametersGauss : public TObject
{
 public:
    MyObjectParametersGauss();
  //MyObjectParameters(int cathode, double chargetot, double Ky3, double Kx3, int Nbdigits);
    virtual ~MyObjectParametersGauss();

  int getfCathode() const { return fCathode; }
  void setfCathode(int cathode) { fCathode = cathode; }
    double getfChargetot() const { return fChargetot; }
    void setfChargetot(double chargetot) { fChargetot = chargetot; }
    double getfSigy() const { return fSigy; }
    void setfSigy(double sigy) { fSigy = sigy; }
    double getfSigx() const { return fSigx; }
    void setfSigx(double sigx) { fSigx = sigx; }
    int getfNbdigits() const { return fNbdigits; }
    void setfNbdigits(int Nbdigits) { fNbdigits = Nbdigits; }
    
    ClassDef(MyObjectParametersGauss, 1);

 private:
    int fCathode;
    double fChargetot;
    double fSigy;
    double fSigx;
    int fNbdigits;

};

}
}

#endif // ALICEO2_MCH_MYOBJECTPARAMETERSGAUSS_H_
