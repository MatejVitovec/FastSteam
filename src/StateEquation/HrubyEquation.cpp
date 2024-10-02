#include <cmath>

#include "HrubyEquation.hpp"

double HrubyEquation::p(double rho, double e) const
{
    return pFunc(calcDelta(rho), calcPsi(e));
}

double HrubyEquation::T(double rho, double e) const
{
    return TFunc(calcDelta(rho), calcPsi(e));
}

double HrubyEquation::s(double rho, double e) const
{

}

double HrubyEquation::h(double rho, double e) const
{

}

double HrubyEquation::a2(double rho, double e) const
{
    double delta = calcDelta(rho);
    double psi = calcPsi(e);
    double Z = compressFactorFunc(delta, psi);
    double etapVal = etap(delta, psi);

    return specGasConst*critT*(-(2.0*delta*Z*etadp(delta, psi))/(std::pow(etapp(delta, psi), 2)) - (Z*Z*etapp(delta, psi))/(std::pow(etapVal, 3)) + WFunc(delta, psi)/etapVal);
}

double HrubyEquation::a(double rho, double e) const
{
    std::sqrt(a2(rho, e));
}


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

double HrubyEquation::compressFactorFunc(double delta, double psi) const
{
    return 1.0 - delta*etard(delta, psi);
}

double HrubyEquation::WFunc(double delta, double psi) const
{
    return 1.0 - 2.0*delta*etard(delta, psi) - delta*delta*etardd(delta, psi);
}

double HrubyEquation::pFunc(double delta, double psi) const
{
    return compressFactorFunc(delta, psi)*TFunc(delta, psi)*specGasConst*delta*critRho;
}

double HrubyEquation::TFunc(double delta, double psi) const
{
    return critT/etap(delta, psi);
}