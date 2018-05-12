#ifndef _CONSTANT_H_
#define _CONSTANT_H_
#include <complex>

namespace mom
{

    const double PI = 3.1415926535897932384626433832795028841971693993751;
    const double pi = 3.1415926535897932384626433832795028841971693993751;
    const double pi2 = (pi * 2.0);
    const double pi4 = (pi * 4.0);
    const double pi_2 = pi / 2.;
    const std::complex<double> complex_j = std::complex<double>(0., 1.);
    const std::complex<double> p2j = pi2 * complex_j;

    const double _eta = (120 * pi);
    const double miu0 = 1.2566370614359172953850573533118e-6;
    const double epsilon0 = 8.841941282883074209382431298473e-12;
    const double C_light = 1. / sqrt(miu0 * epsilon0);

    const double c_wave_speed = 1. / sqrt(miu0 * epsilon0);
}
#endif
