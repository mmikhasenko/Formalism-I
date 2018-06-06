// [06/2018] Misha Mikhasenko, mikhail.mikhasenko@gmail.com

#include <iostream>
#include <vector>
#include <complex>
#include <functional>
#include <map>

#define LAMBDA(x,y,z) ((x)*(x)+(y)*(y)+(z)*(z)-2*(x)*(y)-2*(y)*(z)-2*(z)*(x))

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
}

//------------------------------------------------------------------------------//
//------------------THE MAIN: CREATE ISOBARS, CALL DENSITY----------------------//
//------------------------------------------------------------------------------//

#define mMu  0.1134

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
        std::vector<cd> couplings(isobars.size());
        // the values are to be set
        couplings[0] = cd(3.0, 0.1);
        couplings[1] = cd(5.0, 0.1);
        std::cout << "----> The isobars are created!\n";

        // ----------------------------------------
        // here is my artificial event!
        // I make up some vectors
        double s_input = 2.25;
        // this function returns p1, p2, p3, p4
        std::vector<std::vector<double> > ps = Generator::make_up_some_vectors(s_input, 0.3, 0.0,
                                                                               FormalismI::m1sq, FormalismI::m2sq, FormalismI::m3sq, FormalismI::m4sq);
        std::cout << "----> p1, p2, p3, p4 are generated!";
        // this function makes the decay to the muons
        std::vector<std::vector<double> > qs = Generator::decay_p(ps[0], mMu*mMu, mMu*mMu, 0.3);
        std::vector<double> input_vectors[] = {ps[0],ps[1],ps[2],ps[3],qs[0],qs[1]};
        std::cout << " q1, q2 are generated!\n";

        double density = FormalismI::density(isobars, couplings, input_vectors);
        std::cout << "----> Final result is " << density << "\n";

        /* in order to speed the calculations, one can do:

           1) precalculate the basis functions for the event
           FormalismI::struct_map stmap; FormalismI::construct_structres(stmap, input_vectors);

           2) get s and t
           double s = inv_masssq_of_sum(ps[2],ps[3]), t = inv_masssq_of_sum(ps[0],ps[3]);

           3) call density
           double density = FormalismI::density(isobars, couplings, s,t, stmap);

         */

        return 0.0;
}
