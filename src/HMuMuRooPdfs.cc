#include "HiggsAnalysis/CombinedLimit/interface/HMuMuRooPdfs.h"
#include "RooRealVar.h"
#include "TMath.h"
#include "RooMath.h"
#include "TError.h"
#include <cmath>
#include <complex>
#include <iostream>

ClassImp(RooModZPdf)

RooModZPdf::RooModZPdf(const char *name, const char *title, RooAbsReal& _x, RooAbsReal& _a):
    RooAbsPdf(name, title),
    x("x", "x", this, _x),
    a("a", "a", this, _a),
    b("b", "b", this, 0.0),
    c("c", "c", this, 2.0),
    m("m", "m", this, 91.2),
    w("w", "w", this, 2.5),
    bernCoef("coefficients", "List of Bernstein coefficients", this),
    fixb(true),
    fixc(true),
    fixZMass(true),
    fixZWidth(true)
{
}

RooModZPdf::RooModZPdf(const char *name, const char *title, RooAbsReal& _x, RooAbsReal& _a, RooAbsReal& _b):
    RooAbsPdf(name, title),
    x("x", "x", this, _x),
    a("a", "a", this, _a),
    b("b", "b", this, _b),
    c("c", "c", this, 2.0),
    m("m", "m", this, 91.2),
    w("w", "w", this, 2.5),
    bernCoef("coefficients", "List of Bernstein coefficients", this),
    fixb(false),
    fixc(true),
    fixZMass(true),
    fixZWidth(true)
{
}

RooModZPdf::RooModZPdf(const char *name, const char *title, RooAbsReal& _x, RooAbsReal& _a, RooAbsReal& _b, RooAbsReal& _c):
    RooAbsPdf(name, title),
    x("x", "x", this, _x),
    a("a", "a", this, _a),
    b("b", "b", this, _b),
    c("c", "c", this, _c),
    m("m", "m", this, 91.2),
    w("w", "w", this, 2.5),
    bernCoef("coefficients", "List of Bernstein coefficients", this),
    fixb(false),
    fixc(false),
    fixZMass(true),
    fixZWidth(true)
{
}

RooModZPdf::RooModZPdf(const char *name, const char *title, RooAbsReal& _x, RooAbsReal& _a, RooAbsReal& _b, RooAbsReal& _c, const RooArgList& _coef):
    RooAbsPdf(name, title),
    x("x", "x", this, _x),
    a("a", "a", this, _a),
    b("b", "b", this, _b),
    c("c", "c", this, _c),
    m("m", "m", this, 91.2),
    w("w", "w", this, 2.5),
    bernCoef("coefficients", "List of Bernstein coefficients", this),
    fixb(false),
    fixc(false),
    fixZMass(true),
    fixZWidth(true)
{
    TIterator* coefIter = _coef.createIterator() ;
    RooAbsArg* coef ;
    while((coef = (RooAbsArg*)coefIter->Next())) {
        if (!dynamic_cast<RooAbsReal*>(coef)) {
          std::cout << "RooBernstein::ctor(" << GetName() << ") ERROR: coefficient " << coef->GetName() << " is not of type RooAbsReal" << std::endl ;
          R__ASSERT(0) ;
        }
        bernCoef.add(*coef);
    }
    delete coefIter ;
}

RooModZPdf::RooModZPdf(const char *name, const char *title, RooAbsReal& _x, RooAbsReal& _a, RooAbsReal& _b, RooAbsReal& _c, RooAbsReal& _m):
    RooAbsPdf(name, title),
    x("x", "x", this, _x),
    a("a", "a", this, _a),
    b("b", "b", this, _b),
    c("c", "c", this, _c),
    m("m", "m", this, _m),
    w("w", "w", this, 2.5),
    bernCoef("coefficients", "List of Bernstein coefficients", this),
    fixb(false),
    fixc(false),
    fixZMass(false),
    fixZWidth(true)
{
}

RooModZPdf::RooModZPdf(const char *name, const char *title, RooAbsReal& _x, RooAbsReal& _a, RooAbsReal& _b, RooAbsReal& _c, RooAbsReal& _m, RooAbsReal& _w):
    RooAbsPdf(name, title),
    x("x", "x", this, _x),
    a("a", "a", this, _a),
    b("b", "b", this, _b),
    c("c", "c", this, _c),
    m("m", "m", this, _m),
    w("w", "w", this, _w),
    bernCoef("coefficients", "List of Bernstein coefficients", this),
    fixb(false),
    fixc(false),
    fixZMass(false),
    fixZWidth(false)
{
}

RooModZPdf::RooModZPdf(const RooModZPdf& other, const char* name):
    RooAbsPdf(other, name),
    x("x", this, other.x),
    a("a", this, other.a),
    b("b", this, other.b),
    c("c", this, other.c),
    m("m", this, other.m),
    w("w", this, other.w),
    bernCoef("coefficients", this, other.bernCoef),
    fixb(other.fixb),
    fixc(other.fixc),
    fixZMass(other.fixZMass),
    fixZWidth(other.fixZWidth)
{
}

double RooModZPdf::evaluate() const {
    double zm = 91.2;
    double zw =  2.5;
    double bv =  0.0;
    double cv =  2.0;

    if (!fixZMass ) zm = m;
    if (!fixZWidth) zw = w;
    if (!fixb)      bv = b;
    if (!fixc)      cv = c;

    double val = 0.0;
    val += exp(a*x + bv*x*x);

    // non-relativistic BW
    // val /= (pow(x-zm, cv) + pow(zw/2.0, cv));

    // relativistic BW
    val *= x*x;
    val /= (pow(x*x-zm*zm,cv) + pow(x*x*zw*zw/(zm*zm),cv));

    Int_t degree = bernCoef.getSize();
    if (degree <= 0) return val;

    Double_t xmin = x.min();
    Double_t xv = (x - xmin) / (x.max() - xmin);
    RooFIter iter = bernCoef.fwdIterator();

    Double_t bernval = 1.;
    Double_t coefsum = 0.;
    Double_t coef    = 0.;
    for (Int_t i = 0; i < degree; i++) {
        coef = ((RooAbsReal *)iter.next())->getVal();
        coefsum -= coef;
        bernval += (degree+1.) * TMath::Binomial(degree, i) * pow(xv, degree-i) * pow(1.-xv, i) * coef;
    }
    bernval += coefsum * pow(1.-xv, degree);

    return val * bernval;

}

ClassImp(RooExpPowPdf)

RooExpPowPdf::RooExpPowPdf(const char *name, const char *title, RooAbsReal& _x, RooAbsReal& _a, RooAbsReal& _r):
    RooAbsPdf(name, title),
    x("x" , "x" , this, _x),
    a("a" , "a" , this, _a),
    r("r" , "r" , this, _r)
{
}

RooExpPowPdf::RooExpPowPdf(const RooExpPowPdf& other, const char* name):
    RooAbsPdf(other, name),
    x("x", this, other.x),
    a("a", this, other.a),
    r("r", this, other.r)
{
}

Double_t RooExpPowPdf::evaluate() const {

    return exp(a*x)/pow(x, r);

}

ClassImp(RooExpBWPdf)

RooExpBWPdf::RooExpBWPdf(const char *name, const char *title, RooAbsReal& _x, RooAbsReal& _a, RooAbsReal& _m0, RooAbsReal& _w):
    RooAbsPdf(name, title),
    x ("x" , "x" , this, _x),
    a ("a" , "a" , this, _a),
    m0("m0", "m0", this, _m0),
    w ("w" , "w" , this, _w)
{
}

RooExpBWPdf::RooExpBWPdf(const RooExpBWPdf& other, const char* name):
    RooAbsPdf(other, name),
    x ("x" , this, other.x ),
    a ("a" , this, other.a ),
    m0("m0", this, other.m0),
    w ("w" , this, other.w )
{
}

Double_t RooExpBWPdf::evaluate() const {

    return exp(a*x)/((x-m0)*(x-m0) + w*w/4.);

}

ClassImp(RooPowerLawPdf)

RooPowerLawPdf::RooPowerLawPdf(const char *name, const char *title, RooAbsReal& _x, RooAbsReal& _m0, const RooArgList& _coefList):
    RooAbsPdf(name, title),
    x ("x" , "x" , this, _x),
    m0("m0", "m0", this, _m0),
    coefList("coefs" , "coefs" , this)
{
    TIterator* coefIter = _coefList.createIterator() ;
    RooAbsArg* coef ;
    while((coef = (RooAbsArg*)coefIter->Next())) {
        if (!dynamic_cast<RooAbsReal*>(coef)) {
            std::cout << "RooPowerLawPdf::ctor(" << GetName() << ") ERROR: coefficient " << coef->GetName()
            << " is not of type RooAbsReal" << std::endl ;
            R__ASSERT(0);
        }
        coefList.add(*coef) ;
    }
    delete coefIter ;
}

RooPowerLawPdf::RooPowerLawPdf(const RooPowerLawPdf& other, const char* name):
    RooAbsPdf(other, name),
    x ("x" , this, other.x ),
    m0("m0", this, other.m0),
    coefList("coefList",this,other.coefList)
{
}

Double_t RooPowerLawPdf::evaluate() const {

    RooFIter iter = coefList.fwdIterator();

    if (x-m0 == 0.) return TMath::SignalingNaN(); 

    double val = 0.;
    for (int i = 0; i < coefList.getSize(); i++) {
        val += ((RooAbsReal *)iter.next())->getVal() / pow(x-m0, i+1);
    }

    return val;
}

ClassImp(RooSumTwoExpPdf)

RooSumTwoExpPdf::RooSumTwoExpPdf(const char *name, const char *title, RooAbsReal& _x, RooAbsReal& _a1, RooAbsReal& _a2, RooAbsReal& _f):
    RooAbsPdf(name, title),
    x ("x" , "x" , this, _x),
    a1("a1", "a1", this, _a1),
    a2("a2", "a2", this, _a2),
    f ("f" , "f" , this, _f )
{
}

RooSumTwoExpPdf::RooSumTwoExpPdf(const RooSumTwoExpPdf& other, const char* name):
    RooAbsPdf(other, name),
    x ("x" , this, other.x ),
    a1("a1", this, other.a1),
    a2("a2", this, other.a2),
    f ("f" , this, other.f )
{
}

Double_t RooSumTwoExpPdf::evaluate() const {

    const Double_t xmin = 110.;
    const Double_t xmax = 150.;

    Double_t retval1 = a1 / (exp(a1*xmax) - exp(a1*xmin));
    Double_t retval2 = a2 / (exp(a2*xmax) - exp(a2*xmin));

    if (a1 == 0.) retval1 = 1.0 / (xmax - xmin);
    if (a2 == 0.) retval2 = 1.0 / (xmax - xmin);

    return f * exp(a1*x) * retval1 + (1.-f) * exp(a2*x) * retval2;     

}

Int_t RooSumTwoExpPdf::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* /*rangeName*/) const {

    if (matchArgs(allVars, analVars, x)) return 1;
    else return 0;
}

Double_t RooSumTwoExpPdf::analyticalIntegral(Int_t code, const char* rangeName) const {

    R__ASSERT(code == 1);

    const Double_t xmin = x.min(rangeName);
    const Double_t xmax = x.max(rangeName);

    if (xmin == 110. && xmax == 150.) return 1.;

    Double_t integral_part1 = f;
    if (a1 != 0.) {
        integral_part1 /= exp(a1*150.) - exp(a1*110.);
        integral_part1 *= exp(a1*xmax) - exp(a1*xmin);
    }
    else {
        integral_part1 *= (xmax - xmin)/(150. - 110.);
    }

    Double_t integral_part2 = (1. - f);
    if (a2 != 0.) {
        integral_part2 /= exp(a2*150.) - exp(a2*110.);
        integral_part2 *= exp(a2*xmax) - exp(a2*xmin);
    }
    else {
        integral_part2 *= (xmax - xmin)/(150. - 110.);
    }

    return integral_part1 + integral_part2;

}

ClassImp(RooTwoPowerLawPdf)

RooTwoPowerLawPdf::RooTwoPowerLawPdf(const char *name, const char *title, RooAbsReal& _x, RooAbsReal& _a1, RooAbsReal& _m0):
    RooAbsPdf(name, title),
    x ("x" , "x" , this, _x),
    a1("a1", "a1", this, _a1),
    m0("m0", "m0", this, _m0)
{
}

RooTwoPowerLawPdf::RooTwoPowerLawPdf(const RooTwoPowerLawPdf& other, const char* name):
    RooAbsPdf(other, name),
    x ("x" , this, other.x ),
    a1("a1", this, other.a1),
    m0("m0", this, other.m0)
{
}

Double_t RooTwoPowerLawPdf::evaluate() const {

    const Double_t xmin = 110.;
    const Double_t xmax = 150.;

    Double_t retval1 = 1.0 / log((xmax-m0)/(xmin-m0));
    retval1 *= 1.0/(x-m0);
    retval1 *= a1;

    Double_t retval2 = 1.0 / (1./(xmin-m0) - 1./(xmax-m0));
    retval2 *= 1.0/((x-m0)*(x-m0));
    retval2 *= 1.0 - a1;

    return retval1 + retval2;

}

Int_t RooTwoPowerLawPdf::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* /*rangeName*/) const {

    if (matchArgs(allVars, analVars, x)) return 1;
    else return 0;
}

Double_t RooTwoPowerLawPdf::analyticalIntegral(Int_t code, const char* rangeName) const {

    R__ASSERT(code == 1);

    const Double_t xmin = x.min(rangeName);
    const Double_t xmax = x.max(rangeName);

    if (xmin == 110. && xmax == 150.) return 1.;

    Double_t integral_part1 = a1;
    integral_part1 *= log((xmax-m0)/(xmin-m0)) / log((150.-m0)/(110.-m0));

    Double_t integral_part2 = 1. - a1;
    integral_part2 *= 1./(xmin-m0) - 1./(xmax-m0);
    integral_part2 *= 1. / (1./(110.-m0) - 1./(150.-m0));

    return integral_part1 + integral_part2;

}

ClassImp(RooAKeysPdf)

RooAKeysPdf::RooAKeysPdf() : _verbosity(0), _nEvents(0), _dataPts(0), _dataWgts(0), _weights(0), _sumWgt(0),
	       _mirrorLeft(kFALSE), _mirrorRight(kFALSE), 
	       _asymLeft(kFALSE), _asymRight(kFALSE)
{ 
  // coverity[UNINIT_CTOR]
}


//_____________________________________________________________________________
RooAKeysPdf::RooAKeysPdf(const char *name, const char *title,
			 RooAbsReal& x, RooDataSet& data,
			 Mirror mirror, Double_t rho) :
  RooAbsPdf(name,title),
  _verbosity(0),
  _x("x","Dependent",this,x),
  _nEvents(0),
  _dataPts(0),
  _dataWgts(0),
  _weights(0),
  _mirrorLeft(mirror==MirrorLeft || mirror==MirrorBoth || mirror==MirrorLeftAsymRight),
  _mirrorRight(mirror==MirrorRight || mirror==MirrorBoth || mirror==MirrorAsymLeftRight),
  _asymLeft(mirror==MirrorAsymLeft || mirror==MirrorAsymLeftRight || mirror==MirrorAsymBoth),
  _asymRight(mirror==MirrorAsymRight || mirror==MirrorLeftAsymRight || mirror==MirrorAsymBoth),
  _rho(rho)
{
  // cache stuff about x
  snprintf(_varName, 128,"%s", x.GetName());
  RooRealVar real= (RooRealVar&)(_x.arg());
  _lo = real.getMin();
  _hi = real.getMax();
  _binWidth = (_hi-_lo)/(_nPoints-1);

  // form the lookup table
  LoadDataSet(data);
}



//_____________________________________________________________________________
RooAKeysPdf::RooAKeysPdf(const RooAKeysPdf& other, const char* name):
  RooAbsPdf(other,name),  _verbosity(other._verbosity),
  _x("x",this,other._x), _nEvents(other._nEvents),
  _dataPts(0), _dataWgts(0), _weights(0), _sumWgt(0),
  _mirrorLeft( other._mirrorLeft ), _mirrorRight( other._mirrorRight ),
  _asymLeft(other._asymLeft), _asymRight(other._asymRight),
  _rho( other._rho ) {

  // cache stuff about x
  snprintf(_varName, 128, "%s", other._varName );
  _lo = other._lo;
  _hi = other._hi;
  _binWidth = other._binWidth;

  // copy over data and weights... not necessary, commented out for speed
  //    _dataPts = new Double_t[_nEvents];
  //    _weights = new Double_t[_nEvents];  
  //    for (Int_t i= 0; i<_nEvents; i++) {
  //      _dataPts[i]= other._dataPts[i];
  //      _weights[i]= other._weights[i];
  //    }

  // copy over the lookup table
  for (Int_t i= 0; i<_nPoints+1; i++)
    _lookupTable[i]= other._lookupTable[i];
  
}


//_____________________________________________________________________________
RooAKeysPdf::~RooAKeysPdf() {
  delete[] _dataPts;
  delete[] _dataWgts;
  delete[] _weights;
}


void

//_____________________________________________________________________________
RooAKeysPdf::LoadDataSet( RooDataSet& data) {
  delete[] _dataPts;
  delete[] _dataWgts;
  delete[] _weights;

  // make new arrays for data and weights to fill
  _nEvents= (Int_t) data.numEntries();
  if (_mirrorLeft) _nEvents += data.numEntries();
  if (_mirrorRight) _nEvents += data.numEntries();

  _dataPts  = new Double_t[_nEvents];
  _dataWgts = new Double_t[_nEvents];
  _weights  = 0 ;
  _sumWgt = 0 ;

  Double_t x0(0);
  Double_t x1(0);
  Double_t x2(0);

  Int_t i, idata=0;
  for (i=0; i<data.numEntries(); i++) {
    const RooArgSet *values= data.get(i);
    RooRealVar real= (RooRealVar&)(values->operator[](_varName));

    _dataPts[idata]= real.getVal();
    _dataWgts[idata] = data.weight() ;
    x0 += _dataWgts[idata] ; x1+=_dataWgts[idata]*_dataPts[idata]; x2+=_dataWgts[idata]*_dataPts[idata]*_dataPts[idata];
    idata++;
    _sumWgt+= data.weight() ;

    if (_mirrorLeft) {
      _dataPts[idata]= 2*_lo - real.getVal();
      _dataWgts[idata]= data.weight() ;
      _sumWgt+= data.weight() ;
      idata++;
    }

    if (_mirrorRight) {
      _dataPts[idata]  = 2*_hi - real.getVal();
      _dataWgts[idata] = data.weight() ;
      _sumWgt+= data.weight() ;
      idata++;
    }
  }

  Double_t meanv=x1/x0;
  Double_t sigmav=sqrt(x2/x0-meanv*meanv);
  Double_t h=TMath::Power(Double_t(4)/Double_t(3),0.2)*TMath::Power(_nEvents,-0.2)*_rho;
  Double_t hmin=h*sigmav*sqrt(2.)/10;
  Double_t norm=h*sqrt(sigmav)/(2.0*sqrt(3.0));

  if (_nEvents < _nPoints ) {
    _weights=new Double_t[_nEvents];
    for(Int_t j=0;j<_nEvents;++j) {
      _weights[j]=norm/sqrt(g(_dataPts[j],h*sigmav));
      if (_weights[j]<hmin) _weights[j]=hmin;
    }
  } else {
    _weights=new Double_t[_nPoints];
    for(Int_t k=0;k<_nPoints;++k) {
      _weights[k]=norm/sqrt(g(_lo+k*_binWidth,h*sigmav));
      if (_weights[k]<hmin) _weights[k]=hmin;
    }
  } 
 
  for (i=0;i<_nPoints+1;++i) 
    _lookupTable[i]=evaluateFull( _lo+Double_t(i)*_binWidth );  
}

//_____________________________________________________________________________
Double_t RooAKeysPdf::evaluate() const {
  Int_t i = (Int_t)floor((Double_t(_x)-_lo)/_binWidth);
  if (i<0) {
    std::cerr << "got point below lower bound:"
	 << Double_t(_x) << " < " << _lo
	      << " -- performing linear extrapolation..." << std::endl;
    i=0;
  }
  if (i>_nPoints-1) {
    std::cerr << "got point above upper bound:"
	 << Double_t(_x) << " > " << _hi
	 << " -- performing linear extrapolation..." << std::endl;
    i=_nPoints-1;
  }
  Double_t dx = (Double_t(_x)-(_lo+i*_binWidth))/_binWidth;
  
  // for now do simple linear interpolation.
  // one day replace by splines...
  return (_lookupTable[i]+dx*(_lookupTable[i+1]-_lookupTable[i]));
}

//_____________________________________________________________________________
Double_t RooAKeysPdf::evaluate_weight (Double_t x0) const {
  Int_t i = (Int_t)floor((Double_t(x0)-_lo)/_binWidth);
  if (i<0) {
    if (_verbosity) std::cerr << "got point below lower bound:"
			 << Double_t(x0) << " < " << _lo
			 << " -- performing linear extrapolation..." << std::endl;
    i=0;
  }
  if (i>_nPoints-1) {
    if (_verbosity) std::cerr << "got point above upper bound:"
			      << Double_t(x0) << " > " << _hi
			      << " -- performing linear extrapolation..." << std::endl;
    i=_nPoints-1;
  }
  Double_t dx = (Double_t(x0)-(_lo+i*_binWidth))/_binWidth;
  
  // for now do simple linear interpolation.
  // one day replace by splines...
  return (_weights[i]+dx*(_weights[i+1]-_weights[i]));
}

//_____________________________________________________________________________
Double_t RooAKeysPdf::evaluateFull( Double_t x ) const {
  Double_t y=0;

  if (_nEvents < _nPoints ) {
    for (Int_t i=0;i<_nEvents;++i) {
      Double_t chi=(x-_dataPts[i])/_weights[i];
      y+=_dataWgts[i]*exp(-0.5*chi*chi)/_weights[i];

      if (_asymLeft) {
	chi=(x-(2*_lo-_dataPts[i]))/_weights[i];
	y-=_dataWgts[i]*exp(-0.5*chi*chi)/_weights[i];
      }
      if (_asymRight) {
	chi=(x-(2*_hi-_dataPts[i]))/_weights[i];
	y-=_dataWgts[i]*exp(-0.5*chi*chi)/_weights[i];
      }
    }
  } else {
    for (Int_t i=0;i<_nEvents;++i) {
      Double_t weight = evaluate_weight( _dataPts[i] );
      Double_t chi=(x-_dataPts[i])/weight;
      y+=_dataWgts[i]*exp(-0.5*chi*chi)/weight;
      if (_asymLeft) {
	weight = evaluate_weight( _dataPts[i] );
	chi=(x-(2*_lo-_dataPts[i]))/weight;
	y-=_dataWgts[i]*exp(-0.5*chi*chi)/weight;
      }
      if (_asymRight) {
	weight = evaluate_weight( _dataPts[i] );
	chi=(x-(2*_hi-_dataPts[i]))/weight;
	y-=_dataWgts[i]*exp(-0.5*chi*chi)/weight;
      }
    }
  }
  static const Double_t sqrt2pi(sqrt(2*TMath::Pi()));  
  return y/(sqrt2pi*_sumWgt);
}


//_____________________________________________________________________________
Double_t RooAKeysPdf::g(Double_t x,Double_t sigmav) const {
  
  Double_t c=Double_t(1)/(2*sigmav*sigmav);

  Double_t y=0;
  for (Int_t i=0;i<_nEvents;++i) {
    Double_t r=x-_dataPts[i];
    y+=exp(-c*r*r);
  }
  
  static const Double_t sqrt2pi(sqrt(2*TMath::Pi()));  
  return y/(sigmav*sqrt2pi*_nEvents);
}
