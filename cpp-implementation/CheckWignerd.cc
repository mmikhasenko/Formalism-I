
#include <random>
#include <iostream>
#include <functional>
#include "SpecialFunctions.cc"

int main() {

        std::vector<std::function<double(double)> > fs = {
                [](double z) -> double {
                        return SpecialFunc::wignerd(2,2,2,z);
                },
                [](double z) -> double {
                        return SpecialFunc::wignerd(2,2,-1,z);
                },
                [](double z) -> double {
                        return SpecialFunc::wignerd(1,-1,0,z);
                }
        };
        std::vector<std::function<double(double)> > gs = {
                [](double z) -> double {
                        return (1+z)*(1+z)/4.0;
                },
                [](double z) -> double {
                        double theta = acos(z);
                        return -sin(theta)*(1-z)/2;
                },
                [](double z) -> double {
                        double theta = acos(z);
                        return sin(theta)/sqrt(2);
                }
        };


        std::random_device rd;
        std::default_random_engine en(rd());
        std::uniform_real_distribution<> dist(-1, 1);
        for (uint n = 0; n < fs.size(); n++) {
                std::cout << "The next function\n";
                for (uint i = 0; i < 10; i++) {
                        double z = dist(en);
                        std::cout << "fs["<<n<<"] vs gs["<<n<<"]"
                                  << ", cosTh = " << z << " :"
                                  << fs[n](z) << ", " << gs[n](z)
                                  << ", Delta = " << fs[n](z) - gs[n](z)
                                  << ", Ratio = " << fs[n](z)/gs[n](z)
                                  << "\n";
                }
        }

        return 0.0;
}
