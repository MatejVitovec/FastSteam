#include <cmath>

#include "HrubyEquation.hpp"


double HrubyEquation::eta0(double delta, double psi) const
{
    double out = 0.0;
    for (size_t i = 0; i < 10; i++)
    {
        out += ceoffC[i]*std::pow(calcKsi(psi), i);
    }
    return out;
}

double HrubyEquation::eta0k(double delta, double psi) const
{
    double out = 0.0;
    for (size_t i = 1; i < 10; i++)
    {
        out += i*ceoffC[i]*std::pow(calcKsi(psi), i-1);
    }
    return out;
}

double HrubyEquation::eta0kk(double delta, double psi) const
{
    double out = 0.0;
    for (size_t i = 2; i < 10; i++)
    {
        out += i*(i-1)*ceoffC[i]*std::pow(calcKsi(psi), i-2);
    }
    return out;
}


double HrubyEquation::etar(double delta, double psi) const
{
    return -calcA1(calcKsi(psi))*delta;
}

double HrubyEquation::etard(double delta, double psi) const
{
    return -calcA1(calcKsi(psi));
}

double HrubyEquation::etardd(double delta, double psi) const
{
    return 0.0;
}

double HrubyEquation::etark(double delta, double psi) const
{
    return -calcA1k(calcKsi(psi))*delta;
}

double HrubyEquation::etarkk(double delta, double psi) const
{
    return -calcA1kk(calcKsi(psi))*delta;
}

double HrubyEquation::etardk(double delta, double psi) const
{
    return -calcA1k(calcKsi(psi));
}


double HrubyEquation::eta(double delta, double psi) const
{
    return eta0(delta, psi) + etar(delta, psi);
}

double HrubyEquation::etad(double delta, double psi) const
{
    return -1.0/delta + etard(delta, psi);
}

double HrubyEquation::etadd(double delta, double psi) const
{
    return 1.0/(delta*delta) + etardd(delta, psi);
}

double HrubyEquation::etap(double delta, double psi) const
{
    return std::pow(calcKsi(psi), 2)*(eta0k(delta, psi) + etark(delta, psi));
}

double HrubyEquation::etapp(double delta, double psi) const
{
    return -2.0*calcKsi(psi)*etap(delta, psi) + std::pow(calcKsi(psi), 4)*(eta0kk(delta, psi) + etarkk(delta, psi));
}

double HrubyEquation::etadp(double delta, double psi) const
{
    return std::pow(calcKsi(psi), 2)*etardk(delta, psi);
}
