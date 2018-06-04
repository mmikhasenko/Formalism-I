// Copyright [06/2018] Misha Mikhasenko, mikhail.mikhasenko@gmail.com

#include <iostream>
#include <vector>
#include <complex>
#include <functional>

#define LAMBDA(x,y,z) ((x)*(x)+(y)*(y)+(z)*(z)-2*(x)*(y)-2*(y)*(z)-2*(z)*(x))

typedef std::complex<double> cd;

#define mJpsi 2.98
#define mB   4.2 
#define mPi  0.14
#define mK   0.549

//------------------------------------------------------------------------------//
//------------------------------------------------------------------------------//
//----------------------THE EQUATIONS FROM THE APPENDIX B-----------------------//
//-----------------------------------------------------------EPC ?? ??? ??------//
//------------------------------------------------------------------------------//

namespace FormalismI {

  const double sqrt2 = sqrt(2);
  const double m1sq = mJpsi*mJpsi;
  const double m2sq = mB*mB; 
  const double m3sq = mPi*mPi;
  const double m4sq = mK*mK;

  typedef struct {
    char ch; // either 's': s-channel or 't': t-channel
    double j;     // spin for the resonance
    double l;     // partial l-wave
    std::function<cd(double)> amp;  // lineshape (BW) function
  } isobar;

  //------------------------------------------------------------------------------//

  void check(char channel) {
    if ((channel != 's') && (channel != 't')) {
      std::cerr << "Error: 'channel' must be either 's' or 't', while it is " << channel << "!\n";
    }
  }
  
  double mxsq(char ch) {
    check(ch);
    return (ch=='s') ? m2sq : m3sq;
  }

  
  auto lambda12 = [](double x)->double{return LAMBDA(x,m1sq,m2sq);};
  auto lambda13 = [](double x)->double{return LAMBDA(x,m1sq,m3sq);};
  std::function<double(double)> lambda1x(char ch) {
    check(ch);
    return (ch=='s') ? lambda12 : lambda13;
  }
  auto lambda34 = [](double x)->double{return LAMBDA(x,m3sq,m4sq);};
  auto lambda24 = [](double x)->double{return LAMBDA(x,m2sq,m4sq);};
  
  auto p2 = [](double x)->double{return sqrt(LAMBDA(x,m1sq,m2sq))/(2*sqrt(x));};
  auto p3 = [](double x)->double{return sqrt(LAMBDA(x,m1sq,m3sq))/(2*sqrt(x));};
  std::function<double(double)> px(char ch) {
    check(ch);
    return (ch=='s') ? p2 : p3;
  }


  auto q2 = [](double x)->double{return sqrt(LAMBDA(x,m3sq,m4sq))/(2*sqrt(x));};
  auto q3 = [](double x)->double{return sqrt(LAMBDA(x,m2sq,m4sq))/(2*sqrt(x));};
  std::function<double(double)> qx(char ch) {
    check(ch);
    return (ch=='s') ? q2 : q3;
  }
  
  double x(char ch, double s, double t) {
    check(ch);
    return ch == 's' ? s : t;
  }

  double zx(char ch, double s, double t) {
    check(ch);
    double u = m1sq+m2sq+m3sq+m4sq - s - t;
    return ch == 's' ?
      (s*(t-u)+(m1sq-m2sq)*(m3sq-m4sq))/sqrt(lambda12(s)*lambda34(s))  :
      (t*(s-u)+(m1sq-m3sq)*(m2sq-m4sq))/sqrt(lambda13(t)*lambda24(t));
  }

  //------------------------------------------------------------------------------//

  double clebsch_gordon(uint j1, uint m1, uint j2, uint m2, uint j, uint m) {
    return 1.0;
  }
  double wignerd_hat(uint j, uint m1, uint m2, double z) {
    return 1.0;
  }
  
  double xi_minu(char ch, uint j, double x, double zx) {
    double pq = px(ch)(x)*qx(ch)(x);
    double val = pow(pq,j-1)/(4*M_PI*sqrt2)*sqrt((2*j+1)*(2*j-1)) *
      clebsch_gordon(j-1,0,1,1,j,1) *
      wignerd_hat(j,1,0,zx);
    return val;
  }

  double xi_zero(char ch, uint j, double x, double zx) {
    double val = 0.0;
    return val;
  }

  double xi_plus(char ch, uint j, double x, double zx) {
    double pq = px(ch)(x)*qx(ch)(x);
    double val = pow(pq,j-1)/(4*M_PI*sqrt2)*sqrt((2*j+1)*(2*j-1)) *
      clebsch_gordon(j-1,0,1,1,j,1) *
      wignerd_hat(j,1,0,zx);
    return val;
  }

  std::vector<std::function<double(char,uint,double,double)> > xi =
    {xi_minu, xi_zero, xi_plus};
  
  //------------------------------------------------------------------------------//

  double beta_minu(char ch, uint j, double x, double zx) {
    double pq = px(ch)(x)*qx(ch)(x);
    double val = 4*m1sq*pow(pq,j)/(4*M_PI*lambda1x(ch)(x)) * (x+m1sq-mxsq(ch))/(sqrt2*m1sq)*
      sqrt((2*j+1)*(2*j-1))*(
                             clebsch_gordon(j-1,0,1,0,j,0) *
                             wignerd_hat(j,0,0,zx)/sqrt2 + 
                             clebsch_gordon(j-1,0,1,1,j,1) *
                             wignerd_hat(j,1,0,zx)*zx
                             );
    return val;
  }

  double beta_zero(char ch, uint j, double x, double zx) {
    double val = 0.0;
    return val;
  }

  double beta_plus(char ch, uint j, double x, double zx) {
    double pq = px(ch)(x)*qx(ch)(x);
    double val = 4*m1sq*pow(pq,j)/(4*M_PI*lambda1x(ch)(x)) * (x+m1sq-mxsq(ch))/(sqrt2*m1sq)*
      sqrt((2*j+1)*(2*j+3))*(
                             clebsch_gordon(j+1,0,1,0,j,0) *
                             wignerd_hat(j,0,0,zx)/sqrt2 + 
                             clebsch_gordon(j+1,0,1,1,j,1) *
                             wignerd_hat(j,1,0,zx)*zx
                             );
    return val;
  }

  std::vector<std::function<double(char,uint,double,double)> > beta =
    {xi_minu, xi_zero, xi_plus};

  //------------------------------------------------------------------------------//

  double delta_minu(char ch, uint j, double x, double zx) {
    double val = 0.0;
    return val;
  }

  double delta_zero(char ch, uint j, double x, double zx) {
    double pq = px(ch)(x)*qx(ch)(x);
    double val = sqrt2/(4*M_PI)*(2*j+1)*pow(pq, j-1)*wignerd_hat(j,1,0,zx);
    return val;
  }

  double delta_plus(char ch, uint j, double x, double zx) {
    double val = 0.0;
    return val;
  }

  std::vector<std::function<double(char,uint,double,double)> > delta =
    {delta_minu, delta_zero, delta_plus};

  //------------------------------------------------------------------------------//

  // to be completed later
  cd four_vectors_contraction(char Zi, char Zj, double s, double t) {
    return 1.0;
  }

  //------------------------------------------------------------------------------//

  // combined function for the xi, beta, delta
  double xi_beta_delta(char Zi, char ch, uint j, uint l, double x, double zx) {
    if (abs(j-l)>1) return 0.0;
    uint f_index = l-j+1;
  
    switch(Zi) {
    case 'C': return xi.at(f_index)(ch, j, x, zx);
    case 'B': return beta.at(f_index)(ch, j, x, zx);
    case 'D': return delta.at(f_index)(ch, j, x, zx);
    }
  }

  cd V_form(char ch_i, uint j_i, uint l_i,
            char ch_j, uint j_j, uint l_j,
            double s, double t) {
    cd val = 0;
    for (char Zi : {'B', 'C', 'D'}) {
      for (char Zj : {'B', 'C', 'D'}) {
        // this block can be speeded up
        cd fvc_ij = four_vectors_contraction(Zi, Zj, s, t); //! CAN BE PRECALCULATED
        double xbd_i = xi_beta_delta(Zi, ch_i, j_i, l_i, x(ch_i,s,t), zx(ch_i,s,t)); //! CAN BE PRECALCULATED
        double xbd_j = xi_beta_delta(Zj, ch_j, j_j, l_j, x(ch_i,s,t), zx(ch_i,s,t)); //! CAN BE PRECALCULATED
        fvc_ij *= xbd_i*xbd_j;
        // --------------------------
        val += fvc_ij;
      }
    }
    return val;
  }

  //------------------------------------------------------------------------------//

  double density(std::vector<isobar> &isobars, std::vector<cd> &couplings,
                 double s, double t) {
    // calculate the amplitude
    cd density = 0.0;
    for (uint i = 0; i < isobars.size(); i++) {
      for (uint j = 0; j < isobars.size(); j++) {
        // get couplings
        cd ci = couplings[i];
        cd cj = couplings[j];

        // get isobars description
        auto iso_i = isobars[i];
        auto iso_j = isobars[j];
      
        // get particular matrix elements
        cd heavy_part_ij = FormalismI::V_form(iso_i.ch, iso_i.j, iso_i.l,
                                              iso_j.ch, iso_j.j, iso_j.l,
                                              s,t);
        cd density_ij =
          iso_i.amp(FormalismI::x(iso_i.ch,s,t)) *
          heavy_part_ij *
          std::conj(iso_j.amp(FormalismI::x(iso_j.ch,s,t)));
      
        //-----the-bilinear-form-for-couplings-----//
        density += ci * density_ij * conj(cj);
        //-----------------------------------------//
      }
    }
    std::cout << "Imag part suppose to be zero! Im@Rho = " << density << "\n";
    return real(density);
  }

}  // namespace FormalismI

// not important sugar
template<class T> std::vector<std::pair<uint, T> > enumerate(std::vector<T> input) {
  std::vector<std::pair<uint, T> > output(input.size());
  for (uint i = 0; i < input.size(); i++) output[i] = std::make_pair(i,input[i]);
  return output;
}
  
//------------------------------------------------------------------------------//
//------------------------------------------------------------------------------//
//------------------------THE PARAMETRIZATION OF ISOBARS------------------------//
//------------------------------------------------------------------------------//
//------------------------------------------------------------------------------//

namespace IsobarParams {
  cd Ks_bw(double s) {
    double mKs = 0.875;
    double GammaKs = 0.05;
    cd amp = 1.0/(mKs*mKs - s - cd(0,1.0)*mKs*GammaKs);
    return amp;
  };
}

//------------------------------------------------------------------------------//
//------------------------------------------------------------------------------//
//------------------THE MAIN: CREATE ISOBARS, CALL DENSITY----------------------//
//------------------------------------------------------------------------------//
//------------------------------------------------------------------------------//

int main() {

  // amplitudes for
  std::vector<FormalismI::isobar> isobars;
  // fill up `isobars`-vector here

  // 1) K* S-wave
  FormalismI::isobar Ks_Swave;
  Ks_Swave.ch ='s';
  Ks_Swave.j = 1;
  Ks_Swave.l = 0;
  Ks_Swave.amp = IsobarParams::Ks_bw;
  isobars.push_back(Ks_Swave);

  // 2) K* D-wave
  FormalismI::isobar Ks_Dwave;
  Ks_Dwave.ch ='s';
  Ks_Dwave.j = 1;
  Ks_Dwave.l = 2;
  Ks_Dwave.amp = IsobarParams::Ks_bw;
  isobars.push_back(Ks_Dwave);

  // vector of couplings
  std::vector<cd> couplings(isobars.size());  // the values are to be set

  // call density!
  double s = 1.1, t = 1.2;  
  double density = FormalismI::density(isobars, couplings, s, t);
  std::cout << "Final result is " << density << "\n";

  return 0.0;
}
