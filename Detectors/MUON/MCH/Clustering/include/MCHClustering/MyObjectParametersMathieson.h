#ifndef ALICEO2_MCH_MYOBJECTPARAMETERSMATHIESON_H_
#define ALICEO2_MCH_MYOBJECTPARAMETERSMATHIESON_H_

#include <TROOT.h>
#include <iostream>
#include <vector>
#include <TObject.h>

namespace o2
{
namespace mch
{

using namespace std;

class MyObjectParametersMathieson : public TObject
{
 public:
    MyObjectParametersMathieson();
  //MyObjectParameters(int cathode, double chargetot, double Ky3, double Kx3, int Nbdigits);
    virtual ~MyObjectParametersMathieson();

  int getfCathode() const { return fCathode; }
  void setfCathode(int cathode) { fCathode = cathode; }
    double getfChargetot() const { return fChargetot; }
    void setfChargetot(double chargetot) { fChargetot = chargetot; }
    double getfKy3() const { return fKy3; }
    void setfKy3(double Ky3) { fKy3 = Ky3; }
    double getfKx3() const { return fKx3; }
    void setfKx3(double Kx3) { fKx3 = Kx3; }
    int getfNbdigits() const { return fNbdigits; }
    void setfNbdigits(int Nbdigits) { fNbdigits = Nbdigits; }
    
    ClassDef(MyObjectParametersMathieson, 1);

 private:
    int fCathode;
    double fChargetot;
    double fKy3;
    double fKx3;
    int fNbdigits;

};

}
}

#endif // ALICEO2_MCH_MYOBJECTPARAMETERSMATHIESON_H_
