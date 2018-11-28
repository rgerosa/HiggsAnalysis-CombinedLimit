#ifndef ROO_PDF_INTEGRAL
#define ROO_PDF_INTEGRAL

#include "RooAbsReal.h"
#include "RooRealVar.h"
#include <string>

class RooPdfIntegral : public RooAbsReal {

 public:
  RooPdfIntegral() {
  }
  
 RooPdfIntegral(const char *name, // name for the variable
		const char *title, 
		RooAbsReal& p,  // pdf to be integrated
		RooRealVar& v,  // observable for the PDF
		const char* binstr, // bin-range to integrate
		const char* totstr  // total-range
		):
  RooAbsReal(name, title),
    pdf("pdf", "pdf", this, p),
    var(v),
    integrals("integrals","integrals",this),
    binname(binstr),
    totname(totstr)
    {
      var.setConstant(kFALSE); // take-out constness when object is created
      RooAbsReal* integral = pdf.arg().createIntegral(var, var, binname.c_str());
      integrals.add(*integral);
      integral =  pdf.arg().createIntegral(var, var, totname.c_str());
      integrals.add(*integral);            
    }
  
 RooPdfIntegral(const RooPdfIntegral& other, const char* name=0):
  RooAbsReal(other, name),
    pdf("pdf", this, other.pdf),
    var(other.var),
    integrals("integrals",this,RooListProxy()),
    binname(other.binname),
    totname(other.totname)
      {
	var.setConstant(kFALSE);
	TIterator *intIter = other.integrals.createIterator(); 
	RooAbsReal *fInt;
	while ( (fInt = (RooAbsReal*) intIter->Next()) ){
	  integrals.add(*fInt);
	}
      }
  
  virtual TObject* clone(const char* newname) const {
    return new RooPdfIntegral(*this, newname);
  }
  
  virtual Double_t evaluate() const {
    double retval =  ((RooAbsReal*) integrals.at(0))->getVal()/((RooAbsReal*)integrals.at(1))->getVal();
    return retval;
  }
  
 private:

  RooRealProxy pdf;
  RooRealVar var;
  RooListProxy integrals;
  std::string binname;
  std::string totname;
  
  ClassDef(RooPdfIntegral, 1)
};

#endif
