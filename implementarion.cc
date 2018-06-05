// Copyright [06/2018] Misha Mikhasenko, mikhail.mikhasenko@gmail.com

#include <iostream>
#include <vector>
#include <complex>
#include <functional>
#include <map>

#define LAMBDA(x,y,z) ((x)*(x)+(y)*(y)+(z)*(z)-2*(x)*(y)-2*(y)*(z)-2*(z)*(x))

typedef std::complex<double> cd;

#define mJpsi 3.0969
#define mB   5.279
#define mPi  0.13957
#define mK   0.49368
#define mMu  0.1134

//------------------------------------------------------------------------------//
//------------------------------------------------------------------------------//
//------------------------------------------------------------------------------//

double fourv_prod(const std::vector<double> &p1, const std::vector<double> &p2) {
        if (p1.size() != 4 || p2.size() != 4) std::cerr << "Error: p.size() != 4, something is wrong!\n";
        return p1[3]*p2[3]-p1[0]*p2[0]-p1[1]*p2[1]-p1[2]*p2[2];
}
double inv_masssq(const std::vector<double> &p) {
        return fourv_prod(p,p);
}
double inv_masssq_of_sum(const std::vector<double> &p1, const std::vector<double> &p2) {
        return fourv_prod(p1,p1) + fourv_prod(p2,p2) + 2*fourv_prod(p1,p2);
}

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
        double j; // spin for the resonance
        double l; // partial l-wave
        std::function<cd(double)> amp; // lineshape (BW) function
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


auto lambda12 = [](double x) -> double {
                        return LAMBDA(x,m1sq,m2sq);
                };
auto lambda13 = [](double x) -> double {
                        return LAMBDA(x,m1sq,m3sq);
                };
std::function<double(double)> lambda1x(char ch) {
        check(ch);
        return (ch=='s') ? lambda12 : lambda13;
}
auto lambda34 = [](double x) -> double {
                        return LAMBDA(x,m3sq,m4sq);
                };
auto lambda24 = [](double x) -> double {
                        return LAMBDA(x,m2sq,m4sq);
                };

auto p2 = [](double x) -> double {
                  return sqrt(LAMBDA(x,m1sq,m2sq))/(2*sqrt(x));
          };
auto p3 = [](double x) -> double {
                  return sqrt(LAMBDA(x,m1sq,m3sq))/(2*sqrt(x));
          };
std::function<double(double)> px(char ch) {
        check(ch);
        return (ch=='s') ? p2 : p3;
}


auto q2 = [](double x) -> double {
                  return sqrt(LAMBDA(x,m3sq,m4sq))/(2*sqrt(x));
          };
auto q3 = [](double x) -> double {
                  return sqrt(LAMBDA(x,m2sq,m4sq))/(2*sqrt(x));
          };
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

void construct_C(std::vector<double> &C, const std::vector<double> &pi, const std::vector<double> &pj) {
        double misq = inv_masssq(pi);
        double mjsq = inv_masssq(pj);
        for (uint i = 0; i < 4; i++) C[i] = pi[i]+pj[i];
        double s_t = inv_masssq(C);
        for (uint i = 0; i < 4; i++) C[i] = pi[i]-pj[i] + (misq-mjsq)/s_t * C[i];
}

void construct_B(std::vector<double> &B, const std::vector<double> &pi, const std::vector<double> &pj) {
        for (uint i = 0; i < 4; i++) B[i] = pi[i]+pj[i];
}

// the imaginary i might be omitted, the sign to be checked
void construct_D(std::vector<double> &D, const std::vector<double> &pi, const std::vector<double> &pj, const std::vector<double> &pk) {
        // empty the array
        for (uint mu = 0; mu < 4; mu++) D[mu] = 0.0;
        // fill up with levi-civita
        for (uint si = 0; si < 4; si++) {
                for (uint nu = 0; nu < 4; nu++) {
                        if (nu == si) continue;
                        for (uint rh = 0; rh < 4; rh++) {
                                if (rh == si || rh == nu) continue;
                                for (uint mu = 0; mu < 4; mu++) {
                                        if (mu == si || mu == nu || mu == rh) continue;
                                        int pm = (2*(nu>si)-1)*(2*(rh>si)-1)*(2*(mu>si)-1)*(2*(rh>nu)-1)*(2*(mu>nu)-1)*(2*(mu>rh)-1);
                                        D[mu] += pm*pi[si]*pj[nu]*pk[rh];
                                }
                        }
                }
        }
        // that is all
}

double attac_leptons(const std::vector<double> &L, const std::vector<double> &R, const std::vector<double> &q1, const std::vector<double> &q2) {
        // sum two vectors to get the scalar product
        std::vector<double> q1q2(4);
        for (uint i = 0; i < 4; i++) q1q2[i] = q1[i]+q2[i];
        double msq = inv_masssq(q1q2);
        // products
        double Lq1 = fourv_prod(L,q1), Lq2 = fourv_prod(L,q2);
        double Rq1 = fourv_prod(R,q1), Rq2 = fourv_prod(R,q2);
        double LR = fourv_prod(L,R);
        return Lq1*Rq2 + Lq2*Rq1 - msq/2.0*LR;
}

typedef std::map<std::pair<std::pair<char,char>,std::pair<char,char> >, double > struct_map;

void construct_structres(struct_map &fvmap, const std::vector<double> vecs[6]) {
        // get notations from the paper
        auto p1 = vecs[0], p2 = vecs[1], p3 = vecs[2], p4 = vecs[3];
        auto q1 = vecs[4], q2 = vecs[5];
        std::vector<double> v(4);
        std::map<std::pair<char,char>, std::vector<double> > fvstruct;
        construct_C(v,p3,p4);    fvstruct[std::make_pair('s','C')] = v;
        construct_B(v,p3,p4);    fvstruct[std::make_pair('s','B')] = v;
        construct_D(v,p1,p2,p3); fvstruct[std::make_pair('s','D')] = v;
        construct_C(v,p2,p4);    fvstruct[std::make_pair('t','C')] = v;
        construct_B(v,p2,p4);    fvstruct[std::make_pair('t','B')] = v;
        construct_D(v,p1,p3,p2); fvstruct[std::make_pair('t','D')] = v;

        for (auto & i : fvstruct)
                for (auto & j : fvstruct)
                        fvmap[std::make_pair(i.first,j.first)] = attac_leptons(i.second,j.second, q1, q2);
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
          double s, double t,
          const struct_map &st_map) {

        // calculate the sum
        cd val = 0;
        for (char Zi : {'B', 'C', 'D'}) {
                for (char Zj : {'B', 'C', 'D'}) {
                        double fvc_ij = st_map.at(std::make_pair(std::make_pair(ch_i, Zi),std::make_pair(ch_j, Zj)));
                        // this block can be speeded up
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
               double s, double t,
               const struct_map &st_map) {
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
                                                              s,t, st_map);
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
//---------------------GENERATOR OF ARTIGICIAL EVENTS---------------------------//
//------------------------------------------------------------------------------//
//------------------------------------------------------------------------------//

namespace Generator {

std::vector<std::vector<double> > decay_p(const std::vector<double> &p,
                                          double mu1sq, double mu2sq,
                                          double cosTh, double phi = 0.0) {
        double msq = inv_masssq(p);
        double pstar = sqrt(LAMBDA(msq, mu1sq, mu2sq)/(4*msq));
        double sinTh = sqrt(1.0-cosTh*cosTh);
        std::vector<double> rest_frame_q1 = {pstar*sinTh*cos(phi),
                                             pstar*sinTh*sin(phi),
                                             pstar*cosTh,
                                             sqrt(pstar*pstar+mu1sq)};
        std::vector<double> rest_frame_q2 = {-rest_frame_q1[0],
                                             -rest_frame_q1[1],
                                             -rest_frame_q1[2],
                                             sqrt(pstar*pstar+mu2sq)};
        // boost
        auto boost = [](std::vector<double> &p, double gamma) -> void {
                             double gb = sqrt(gamma*gamma-1.0);
                             double pp = p[2]*gamma+p[3]*gb;
                             double E  = p[3]*gamma+p[2]*gb;
                             p[2] = pp; p[3] = E;
                     };
        boost(rest_frame_q1, p[3]/sqrt(msq));
        boost(rest_frame_q2, p[3]/sqrt(msq));
        return {rest_frame_q1, rest_frame_q2};
}

std::vector<std::vector<double> > make_up_some_vectors(double s, double cosTh, double phi,
                                                       double m1sq, double m2sq, double m3sq, double m4sq) {
        std::vector<double> p1 = {0.0, 0.0, -sqrt(LAMBDA(s,m1sq,m2sq))/(2*sqrt(s)), (m2sq-s-m1sq)/(2*sqrt(s))};
        std::vector<double> p2 = {0.0, 0.0, -sqrt(LAMBDA(s,m1sq,m2sq))/(2*sqrt(s)), (s-m1sq+m2sq)/(2*sqrt(s))};
        double pstar = sqrt(LAMBDA(s,m3sq,m4sq))/(2*sqrt(s));
        std::vector<double> tot = {0.0, 0.0, 0.0, sqrt(s)};
        std::vector<std::vector<double> > p3p4 = decay_p(tot, m3sq, m4sq, cosTh, phi);
        return {p1,p2,p3p4[0],p3p4[1]};
}

}


//------------------------------------------------------------------------------//
//------------------------------------------------------------------------------//
//------------------THE MAIN: CREATE ISOBARS, CALL DENSITY----------------------//
//------------------------------------------------------------------------------//
//------------------------------------------------------------------------------//

int main() {

        // ----------------------------------------
        // here is my amplitude event
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
        std::vector<cd> couplings(isobars.size()); // the values are to be set
        std::cout << "----> The isobars are created!\n";

        // ----------------------------------------
        // here is my artificial event!
        // I make up some vectors
        double s_input = 2.25;
        std::vector<std::vector<double> > ps = Generator::make_up_some_vectors(s_input, 0.3, 0.0,
                FormalismI::m1sq, FormalismI::m2sq, FormalismI::m3sq, FormalismI::m4sq);
        std::cout << "----> p1, p2, p3, p4 are generated!\n";
        std::vector<std::vector<double> > qs = Generator::decay_p(ps[0], mMu*mMu, mMu*mMu, 0.3);
        std::vector<double> input_vectors[] = {ps[0],ps[1],ps[2],ps[3],qs[0],qs[1]};
        std::cout << "----> An artificial event is generated!\n";

        // ----------------------------------------
        // call density!
        FormalismI::struct_map stmap; FormalismI::construct_structres(stmap, input_vectors);
        double s = inv_masssq_of_sum(ps[2],ps[3]), t = inv_masssq_of_sum(ps[0],ps[3]);
        double density = FormalismI::density(isobars, couplings, s, t, stmap);
        std::cout << "----> Final result is " << density << "\n";

        return 0.0;
}
