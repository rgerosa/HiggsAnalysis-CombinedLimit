#ifndef ROO_PDF_WITH_UNCERTAINTIES
#define ROO_PDF_WITH_UNCERTAINTIES

#include "RooAbsReal.h"
#include "RooRealVar.h"
#include <string>

class RooAbsPdfWithUncertainties : public RooAbsPdf {

 public:
  RooAbsPdfWithUncertainties() {
  }
  
 RooAbsPdfWithUncertainties(const char *name, // name for the variable
			    const char *title, 
			    RooAbsPdf &  pdf,  // pdf to be integrated
			    RooRealVar & var  // observable for the PDF
			    ):
  RooAbsPdf(name, title),
    pdf("pdf", "pdf", this, pdf),
    var("var","var",this,var),
    nuisances("nuisances","nuisances",this)
      {
	nuisanceNBins.clear();
	nuisanceRange.clear();
	nuisanceShiftUp.clear();
	nuisanceShiftDw.clear();
      }
  
 RooAbsPdfWithUncertainties(const RooAbsPdfWithUncertainties& other, const char* name=0):
  RooAbsPdf(other, name),
    pdf("pdf",this,other.pdf),
    var("var",this,other.var),
    nuisances("nuisances",this,other.nuisances)
      {
	nuisanceNBins = other.nuisanceNBins;
	nuisanceRange = other.nuisanceRange;
	nuisanceShiftUp = other.nuisanceShiftUp;
	nuisanceShiftDw = other.nuisanceShiftDw;
      }

  virtual TObject* clone(const char* newname) const {
    return new RooAbsPdfWithUncertainties(*this, newname);
  }

  // add the uncertainty
  void addUncertainty(RooRealVar& nuisance, const int & nbins, const double & xmin, const double & xmax, const std::vector<double> & nuisUp, const std::vector<double> & nuisDw){

    if(int(nuisUp.size()) != nbins or int(nuisDw.size()) != nbins){
      std::cout<<"Error --> incoherent size of vectors for nuisance value and Nbins: up "<<nuisUp.size()<<" down "<<nuisDw.size()<<" bins "<<nbins<<endl;
      return;
    }

    nuisances.add(nuisance);
    nuisanceNBins.push_back(nbins);
    nuisanceRange.push_back(std::pair<double,double> (xmin,xmax));
    nuisanceShiftUp.push_back(nuisUp);
    nuisanceShiftDw.push_back(nuisDw);            

  }

  void addUncertainty(RooRealVar& nuisance, const int & nbins, const double & xmin, const double & xmax, const std::vector<double> & nuisVal){

    if(int(nuisVal.size()) != nbins){
      std::cout<<"Error --> incoherent size of vectors for nuisance value and Nbins: "<<nuisVal.size()<<" bins "<<nbins<<endl;
      return;
    }
    addUncertainty(nuisance,nbins,xmin,xmax,nuisVal,nuisVal);
  }


  virtual Double_t evaluate() const { // evaluate the pdf

    // take the bare pdf value    
    double pdfVal = ((RooAbsPdf*) pdf.absArg())->getVal();

    // loop on nuisances
    if(nuisances.getSize() != 0){
      double pdfCorrection = 1;
      for(int inuis = 0; inuis < nuisances.getSize(); inuis++){
	double nuisVal = ((RooRealVar*) nuisances.at(inuis))->getVal();
	if(var < nuisanceRange[inuis].first or
	   var > nuisanceRange[inuis].second) // if x is outside the validity range of the nuisance, pdf value unchanged
	  pdfCorrection *= 1;
	else{	   // correct by the nuisance
	  double shift = 0;
	  if(var == nuisanceRange[inuis].second) shift  = -0.000001;
	  if(var == nuisanceRange[inuis].first)  shift  = 0.000001;
	  int ibin = int((var-nuisanceRange[inuis].first+shift)/(nuisanceRange[inuis].second-nuisanceRange[inuis].first)*nuisanceNBins[inuis]);
	  if(nuisVal >= 0) pdfCorrection *= (1+nuisVal*nuisanceShiftUp[inuis].at(ibin));
	  else pdfCorrection *= (1-nuisVal*nuisanceShiftDw[inuis].at(ibin));	  
	}
      }
      return pdfCorrection*pdfVal;
    }
    else 
      return pdfVal;    
  }

 private:

  RooRealProxy pdf;
  RooRealProxy var;
  RooListProxy nuisances;// centered around zero by construction

  std::vector<int> nuisanceNBins;
  std::vector<std::pair<double,double> > nuisanceRange;
  std::vector<std::vector<double> > nuisanceShiftUp;
  std::vector<std::vector<double> > nuisanceShiftDw;

  ClassDef(RooAbsPdfWithUncertainties, 1)
};

#endif
