// [06/2018] Misha Mikhasenko, mikhail.mikhasenko@gmail.com

#include <functional>
#include <iostream>
#include <complex>
#include <vector>
#include <random>
#include <map>
#include <iomanip>

#define _USE_MATH_DEFINES
#include <cmath>

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
//----------------------THE EQUATIONS FROM THE APPENDIX B-----------------------//
//-----------------------------------------------------------EPC ?? ??? ??------//
#include "FormalismI.cc"

//------------------------------------------------------------------------------//
//---------------------GENERATOR OF ARTIGICIAL EVENTS---------------------------//
//------------------------------------------------------------------------------//
#include "Generator.cc"

//------------------------------------------------------------------------------//
//------------------------THE PARAMETRIZATION OF ISOBARS------------------------//
//------------------------------------------------------------------------------//

namespace IsobarParams {
cd Ks_bw(double s) {
        double mKs = 0.875;
        double GammaKs = 0.05;
        cd amp = 1.0/(mKs*mKs - s - cd(0,1.0)*mKs*GammaKs);
        return amp;
};

cd Zc_bw(double t) {
        double mZc = 4.43;
        double GammaZc = 0.107;
        cd amp = 1.0/(mZc*mZc - t - cd(0,1.0)*mZc*GammaZc);
        return amp;
};

}

//------------------------------------------------------------------------------//
//------------------THE MAIN: CREATE ISOBARS, CALL DENSITY----------------------//
//------------------------------------------------------------------------------//

#define mMu  0.1134

int main() {
        double m1sq = FormalismI::m1sq, m2sq = FormalismI::m2sq,
               m3sq = FormalismI::m3sq, m4sq = FormalismI::m4sq;
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

        // 2) Zc S-wave
        FormalismI::isobar Zc_Swave;
        Zc_Swave.ch ='t';
        Zc_Swave.j = 1;
        Zc_Swave.l = 0;
        Zc_Swave.amp = IsobarParams::Zc_bw;
        isobars.push_back(Zc_Swave);

        // vector of couplings
        std::vector<cd> couplings(isobars.size());
        // the values are to be set
        couplings[0] = cd(1.0, 0.0);
        couplings[1] = cd(-3.0, -0.1);
        // std::cout << "----> The isobars are created!\n";

        // ----------------------------------------
        std::random_device rd;
        std::mt19937 en(rd());
        std::uniform_real_distribution<> dist_cosT(-1.0, 1.0);
        std::uniform_real_distribution<> dist_phi(-M_PI, M_PI);

        // create s-distibutrion according to the projected phase space
        auto density = [](double s) -> double {
                               return sqrt(FormalismI::lambda12(s)*FormalismI::lambda34(s))/s;
                       };
        std::vector<double> sv(100);
        sv[0] = pow(sqrt(m3sq)+sqrt(m4sq),2);  sv[sv.size()-1] = pow(sqrt(m2sq)-sqrt(m1sq),2);
        std::vector<double> w(sv.size());
        w[0] = 0; w[w.size()-1] = 0;
        for (uint i = 1; i < sv.size()-1; i++) {
                sv[i] = sv[0] + (sv[sv.size()-1]-sv[0])/(sv.size()-1)*i;
                w[i] = density(sv[i]);
        }
        std::piecewise_linear_distribution<> dist_s(sv.begin(), sv.end(), w.begin());


        // loop over the events
        const uint Nev = 10000;
        for (uint i = 0; i < Nev; i++) {
                double s_input = dist_s(en);
                double cosTh = dist_cosT(en);
                double phi = dist_phi(en);

                std::vector<std::vector<double> > ps = Generator::make_up_some_vectors(s_input, cosTh, phi,
                                                                                       m1sq, m2sq, m3sq, m4sq);
                // this function makes the decay to the muons
                double cosTh_mu = dist_cosT(en);
                double phi_mu = dist_phi(en);
                std::vector<std::vector<double> > qs = Generator::decay_p(ps[0], mMu*mMu, mMu*mMu, cosTh_mu, phi_mu);
                std::vector<std::vector<double> > input_vectors{ps[0],ps[1],ps[2],ps[3],qs[0],qs[1]};
                // for (auto & v : input_vectors) {
                //         for (auto j : v) std::cout << j << ", ";
                //         std::cout << "\n";
                // }
                double density = FormalismI::density(isobars, couplings, input_vectors);
                double t_input = inv_masssq_of_sum(ps[0], ps[2]);
                std::cout << "s, t, density: " << s_input << " " << t_input << " " << density << "\n";
        }

        return 0.0;
}
