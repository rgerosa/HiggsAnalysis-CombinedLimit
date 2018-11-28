#ifndef HMUMUROOPDFS
#define HMUMUROOPDFS

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooAbsReal.h"
#include "RooListProxy.h"
#include "RooDataSet.h"

class RooModZPdf : public RooAbsPdf {
    public:
        RooModZPdf() {
        };

        RooModZPdf(const char *name, const char *title, RooAbsReal& _x, RooAbsReal& _a);
        RooModZPdf(const char *name, const char *title, RooAbsReal& _x, RooAbsReal& _a, RooAbsReal& _b);
        RooModZPdf(const char *name, const char *title, RooAbsReal& _x, RooAbsReal& _a, RooAbsReal& _b, RooAbsReal& _c);
        RooModZPdf(const char *name, const char *title, RooAbsReal& _x, RooAbsReal& _a, RooAbsReal& _b, RooAbsReal& _c, RooAbsReal& _m);
        RooModZPdf(const char *name, const char *title, RooAbsReal& _x, RooAbsReal& _a, RooAbsReal& _b, RooAbsReal& _c, RooAbsReal& _m, RooAbsReal& _w);
        RooModZPdf(const char *name, const char *title, RooAbsReal& _x, RooAbsReal& _a, RooAbsReal& _b, RooAbsReal& _c, const RooArgList& _coef);

        RooModZPdf(const RooModZPdf& other, const char* name=0) ;

        virtual TObject* clone(const char* newname) const {
            return new RooModZPdf(*this,newname);
        }

        inline virtual ~RooModZPdf() {
        }

    protected:
        RooRealProxy x;
        RooRealProxy a;
        RooRealProxy b;
        RooRealProxy c;
        RooRealProxy m;
        RooRealProxy w;

        RooListProxy bernCoef;

        bool fixb;	
        bool fixc;	
        bool fixZMass;	
        bool fixZWidth;	

        Double_t evaluate() const ;

    private:
        ClassDef(RooModZPdf,1)
};

class RooExpPowPdf : public RooAbsPdf {

    public:
        RooExpPowPdf() {
        };

        RooExpPowPdf(const char *name, const char *title, RooAbsReal& _x, RooAbsReal& _a, RooAbsReal& _r);

        RooExpPowPdf(const RooExpPowPdf& other, const char* name=0) ;

        virtual TObject* clone(const char* newname) const {
            return new RooExpPowPdf(*this,newname);
        }

        inline virtual ~RooExpPowPdf() {
        }

    protected:
        RooRealProxy x;
        RooRealProxy a;
        RooRealProxy r;

        Double_t evaluate() const ;

    private:

        ClassDef(RooExpPowPdf,1)

};

class RooExpBWPdf : public RooAbsPdf {

    public:
        RooExpBWPdf() {
        };

        RooExpBWPdf(const char *name, const char *title, RooAbsReal& _x, RooAbsReal& _a, RooAbsReal& _m0, RooAbsReal& _w);

        RooExpBWPdf(const RooExpBWPdf& other, const char* name=0) ;

        virtual TObject* clone(const char* newname) const {
            return new RooExpBWPdf(*this,newname);
        }

        inline virtual ~RooExpBWPdf() {
        }

    protected:
        RooRealProxy x;
        RooRealProxy a;
        RooRealProxy m0;
        RooRealProxy w;

        Double_t evaluate() const ;

    private:

        ClassDef(RooExpBWPdf,1)

};

class RooPowerLawPdf : public RooAbsPdf {

    public:
        RooPowerLawPdf() {
        };

        RooPowerLawPdf(const char *name, const char *title, RooAbsReal& _x, RooAbsReal& _m0, const RooArgList &_coefList);

        RooPowerLawPdf(const RooPowerLawPdf& other, const char* name=0) ;

        virtual TObject* clone(const char* newname) const {
            return new RooPowerLawPdf(*this, newname);
        }

        inline virtual ~RooPowerLawPdf() {
        }

    protected:
        RooRealProxy x;
        RooRealProxy m0;
        RooListProxy coefList;

        Double_t evaluate() const ;

    private:

        ClassDef(RooPowerLawPdf,1)

};

class RooSumTwoExpPdf : public RooAbsPdf {

    public:
        RooSumTwoExpPdf() {
        };

        RooSumTwoExpPdf(const char *name, const char *title, RooAbsReal& _x, RooAbsReal& _a1, RooAbsReal& _a2, RooAbsReal& _f);

        RooSumTwoExpPdf(const RooSumTwoExpPdf& other, const char* name=0) ;

        virtual TObject* clone(const char* newname) const {
            return new RooSumTwoExpPdf(*this, newname);
        }

        inline virtual ~RooSumTwoExpPdf() {
        }

        Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const;

        Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const;

    protected:
        RooRealProxy x;
        RooRealProxy a1;
        RooRealProxy a2;
        RooRealProxy f;

        Double_t evaluate() const ;

    private:

        ClassDef(RooSumTwoExpPdf,1)

};

class RooTwoPowerLawPdf : public RooAbsPdf {

    public:
        RooTwoPowerLawPdf() {
        };

        RooTwoPowerLawPdf(const char *name, const char *title, RooAbsReal& _x, RooAbsReal& _a1, RooAbsReal& _m0);

        RooTwoPowerLawPdf(const RooTwoPowerLawPdf& other, const char* name=0) ;

        virtual TObject* clone(const char* newname) const {
            return new RooTwoPowerLawPdf(*this, newname);
        }

        inline virtual ~RooTwoPowerLawPdf() {
        }

        Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const;

        Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const;

    protected:
        RooRealProxy x;
        RooRealProxy a1;
        RooRealProxy m0;

        Double_t evaluate() const ;

    private:

        ClassDef(RooTwoPowerLawPdf,1)

};

class RooAKeysPdf : public RooAbsPdf {

 public:
  
  enum Mirror { NoMirror, MirrorLeft, MirrorRight, MirrorBoth,
		MirrorAsymLeft, MirrorAsymLeftRight,
		MirrorAsymRight, MirrorLeftAsymRight,
		MirrorAsymBoth };
  RooAKeysPdf();
  RooAKeysPdf(const char *name, const char *title,
	      RooAbsReal& x, RooDataSet& data, Mirror mirror= NoMirror,
	      Double_t rho=1);
  RooAKeysPdf(const RooAKeysPdf& other, const char* name=0);
  virtual TObject* clone(const char* newname) const {return new RooAKeysPdf(*this,newname); }
  virtual ~RooAKeysPdf();
  
  void LoadDataSet( RooDataSet& data);

  Int_t _verbosity;

 protected:
  
  RooRealProxy _x ;
  Double_t evaluate() const;

 private:
  
  Double_t evaluateFull(Double_t x) const;
  Double_t evaluate_weight(Double_t x0) const;
  
  Int_t _nEvents;
  Double_t *_dataPts;  //[_nEvents]
  Double_t *_dataWgts; //[_nEvents]
  Double_t *_weights;  //[_nEvents]
  Double_t _sumWgt ;
  
  enum { _nPoints = 1000 };
  Double_t _lookupTable[_nPoints+1];
  
  Double_t g(Double_t x,Double_t sigma) const;

  Bool_t _mirrorLeft, _mirrorRight;
  Bool_t _asymLeft, _asymRight;

  // cached info on variable
  Char_t _varName[128];
  Double_t _lo, _hi, _binWidth;
  Double_t _rho;
  
  ClassDef(RooAKeysPdf,1) // One-dimensional non-parametric kernel estimation p.d.f.
};

#endif
