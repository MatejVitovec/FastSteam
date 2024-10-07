#include <cmath>
#include <vector>

#include "HrubyEquation.hpp"

HrubyEquation::HrubyEquation()
{
    critT = 647.096;
    critRho = 322.0;
    specGasConst = 461.51805;

    beta = {3.0136, 0.21825};
    ceoffC = {7.7419423096, 69.178357300, -474.69167775, 6104.0685828, -71467.778757, 794338.56008, -6808358.4642, 40191412.675, -140094025.97, 218677291.14};
    coeffB = {0.5125, -0.00999, 0.000000322};
}

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
    return sFunc(calcDelta(rho), calcPsi(e));
}

double HrubyEquation::h(double rho, double e) const
{
    return e + p(rho, e)/rho;
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

double HrubyEquation::eFromRhoP(double rho, double p, double guessE) const
{
    double delta = calcDelta(rho);

    double psi = nonLinearSolver.solve([=](double val) { return pFunc(delta, val) - p; },
                                       [=](double val) { return pDDeltaFunc(delta, val); },
                                       calcPsi(guessE));
    
    return psi*(critT*specGasConst);
}

double HrubyEquation::eFromRhoT(double rho, double T, double guessE) const
{
    double delta = calcDelta(rho);

    double psi = nonLinearSolver.solve([=](double val) { return TFunc(delta, val) - T; },
                                       [=](double val) { return TDDeltaFunc(delta, val); },
                                       calcPsi(guessE));
    
    return psi*(critT*specGasConst);
}

double HrubyEquation::eFromRhoS(double rho, double s, double guessE) const
{
    double delta = calcDelta(rho);

    double psi = nonLinearSolver.solve([=](double val) { return sFunc(delta, val) - s; },
                                       [=](double val) { return sDDeltaFunc(delta, val); },
                                       calcPsi(guessE));
    
    return psi*(critT*specGasConst);
}

double HrubyEquation::eFromRhoH(double rho, double h, double guessE) const
{
    double delta = calcDelta(rho);

    double psi = nonLinearSolver.solve([=](double val) { return hFunc(delta, val) - h; },
                                       [=](double val) { return sDDeltaFunc(delta, val); },
                                       calcPsi(guessE));
    
    return psi*(critT*specGasConst);
}

double HrubyEquation::rhoFromEP(double e, double p, double guessRho) const
{
    double psi = calcPsi(e);

    double delta = nonLinearSolver.solve([=](double val) { return pFunc(val, psi) - p; },
                                         [=](double val) { return pDPsiFunc(val, psi); },
                                         calcPsi(guessRho));
    
    return delta*critRho;
}

double HrubyEquation::rhoFromET(double e, double T, double guessRho) const
{
    double psi = calcPsi(e);

    double delta = nonLinearSolver.solve([=](double val) { return TFunc(val, psi) - T; },
                                         [=](double val) { return TDPsiFunc(val, psi); },
                                         calcPsi(guessRho));
    
    return delta*critRho;
}

double HrubyEquation::rhoFromES(double e, double s, double guessRho) const
{
    double psi = calcPsi(e);

    double delta = nonLinearSolver.solve([=](double val) { return sFunc(val, psi) - s; },
                                         [=](double val) { return sDPsiFunc(val, psi); },
                                         calcPsi(guessRho));
    
    return delta*critRho;
}

double HrubyEquation::rhoFromEH(double e, double h, double guessRho) const
{
    double psi = calcPsi(e);

    double delta = nonLinearSolver.solve([=](double val) { return hFunc(val, psi) - h; },
                                         [=](double val) { return hDPsiFunc(val, psi); },
                                         calcPsi(guessRho));
    
    return delta*critRho;
}

std::pair<double, double> HrubyEquation::RhoTFromSP(double s, double p, double guessRho, double guessE) const
{
    std::pair<double, double> result = nonLinearSolver.solve([=](double val1, double val2) { return sFunc(val1, val2) - s; },
                                                             [=](double val1, double val2) { return sDDeltaFunc(val1, val2); },
                                                             [=](double val1, double val2) { return sDPsiFunc(val1, val2); },
                                                             [=](double val1, double val2) { return pFunc(val1, val2) - p; },
                                                             [=](double val1, double val2) { return pDDeltaFunc(val1, val2); },
                                                             [=](double val1, double val2) { return pDPsiFunc(val1, val2); },
                                                             calcDelta(guessRho),
                                                             calcPsi(guessE));

    return std::make_pair(result.first*critRho, result.second*(critT*specGasConst));
}


//private functions

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

double HrubyEquation::sFunc(double delta, double psi) const
{
    return eta(delta, psi)*specGasConst;
}

double HrubyEquation::hFunc(double delta, double psi) const
{
    return psi*(critT*specGasConst) + pFunc(delta, psi)/(delta*critRho);
}

double HrubyEquation::pDDeltaFunc(double delta, double psi) const
{
    return 0.0; //TODO
}

double HrubyEquation::TDDeltaFunc(double delta, double psi) const
{
    return 0.0; //TODO
}

double HrubyEquation::sDDeltaFunc(double delta, double psi) const
{
    return 0.0; //TODO
}

double HrubyEquation::hDDeltaFunc(double delta, double psi) const
{
    return 0.0; //TODO
}


double HrubyEquation::pDPsiFunc(double delta, double psi) const
{
    return 0.0; //TODO
}

double HrubyEquation::TDPsiFunc(double delta, double psi) const
{
    return 0.0; //TODO
}

double HrubyEquation::sDPsiFunc(double delta, double psi) const
{
    return 0.0; //TODO
}

double HrubyEquation::hDPsiFunc(double delta, double psi) const
{
    return 0.0; //TODO
}