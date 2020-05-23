#ifndef ALICEO2_MCH_MYOBJECTDIGITS_H_
#define ALICEO2_MCH_MYOBJECTDIGITS_H_

#include <TROOT.h>
#include "TObject.h"
#include <iostream>
#include <vector>
#include "MCHBase/Digit.h"

namespace o2
{
namespace mch
{

using namespace std;

//class MyObjectDigits : public TObject
//{
//
//public:
//
//    MyObjectDigits() = default;
//    virtual ~MyObjectDigits() = default;
//
//    std::vector<Digit> fDigits;
//    void setfDigits(std::vector<Digit> &digits) { fDigits = digits; }
//
//
//    ClassDef(MyObjectDigits,1); // A value of 1 indicates that this class is streamable.
//
//};

class MyObjectDigits : public TObject
{
 public:
  MyObjectDigits();
  virtual ~MyObjectDigits();

  std::vector<Digit> getfDigits() const { return fDigits; }
  void setfDigits(std::vector<Digit> vectdigits) { fDigits = vectdigits; }

  ClassDef(MyObjectDigits, 1);

 private:
  std::vector<Digit> fDigits;
};

} // namespace mch
} // namespace o2

#endif // ALICEO2_MCH_MYOBJECTDIGITS_H_
