// [06/2018] Misha Mikhasenko, mikhail.mikhasenko@gmail.com

#include <iostream>
#include <cmath>
#include <vector>
#include <gsl/gsl_sf_coupling.h>
#include <gsl/gsl_sf_gamma.h>

namespace SpecialFunc {
//------------------------------------------------------------------------------//

double clebsch_gordon(uint j1, uint m1, uint j2, uint m2, uint j, uint m) {
        return (j1-j2+m % 2 == 1 ? -1.0 : 1.0) * sqrt(2*j+1) *
               gsl_sf_coupling_3j(2*j1, 2*j2, 2*j, 2*m1, 2*m2, -2*m);
}

// The implemented expression can be find in Wikipedia
//     https://en.wikipedia.org/wiki/Jacobi_polynomials
double jacobi_pols(uint n, uint a, uint b, double x) {
        // here comes the calculation of log(n!)
        const std::vector<double> logfact = {0.0, 0.0,
                                             6.93147180559945309e-1, 1.79175946922805500e00,
                                             3.17805383034794562e00, 4.78749174278204599e00,
                                             6.57925121201010100e00, 8.52516136106541430e00,
                                             1.06046029027452502e01, 1.28018274800814696e01,
                                             1.51044125730755153e01, 1.75023078458738858e01,
                                             1.99872144956618861e01, 2.25521638531234229e01,
                                             2.51912211827386815e01, 2.78992713838408916e01,
                                             3.06718601060806728e01, 3.35050734501368889e01,
                                             3.63954452080330536e01, 3.93398841871994940e01,
                                             4.23356164607534850e01, 4.53801388984769080e01,
                                             4.84711813518352239e01, 5.16066755677643736e01};
        if (n+a >= logfact.size() || n+b >= logfact.size()) {
                std::cerr << "Error: j is too high, please check the implementation of jacobi polynomials!\n";
                return 0.0;
        }

        double ls = log((1.0-x)/2.0);
        double lc = log((1.0+x)/2.0);
        double res = 0.0;
        for (uint s = 0; s <= n; s++) {
                double logs = logfact[n+a] + logfact[n+b]-logfact[n-s]-logfact[a+s]-logfact[s]-logfact[n+b-s];
                double args = s*ls + (n-s)*lc;
                res += (s % 2 == 0 ? 1.0 : -1.0) * exp(logs+args);
        }
        return res;
}

// The reference  for the relation between WignerD and Jacobi polynomials is
// Eq. (3.74) of L. Biedenharn, J. Louck, and P. Carruthers, Angular Momentum in Quantum Physics: Theory and Application
// see also (B1) of the FormalismII paper.
double wignerd_hat(uint j, uint m1, uint m2, double z) {
        double factor = (abs(m1-m2)+m1-m2)/2 % 2 == 0 ? 1.0 : -1.0;
        int am1 = abs(m1), am2 = abs(m2);
        double M = (am1 > am2) ? am1 : am2;
        double N = (am1 < am2) ? am2 : am1;
        return factor/pow(2,M)*
               sqrt(gsl_sf_gamma(j-M+1)*gsl_sf_gamma(j+M+1)/gsl_sf_gamma(j-N+1)/gsl_sf_gamma(j+N+1))*
               jacobi_pols(j-M, abs(m1-m2),abs(m1+m2), z);
}

}
