#include <cmath>
#include <iostream>
#include <fstream>

#include "Helmholtz.hpp"

Helmholtz::Helmholtz()
{
    nonLinearSolver = NewtonMethod();
}

double Helmholtz::calcDelta(double rho) const
{
    return rho/critRho;
}

double Helmholtz::calcTau(double T) const
{
    return critT/T;
}

double Helmholtz::p(double rho, double T) const
{
    return pFunc(calcDelta(rho), calcTau(T));
}

double Helmholtz::e(double rho, double T) const
{
    return eFunc(calcDelta(rho), calcTau(T));
}

double Helmholtz::s(double rho, double T) const
{
    return sFunc(calcDelta(rho), calcTau(T));
}

double Helmholtz::h(double rho, double T) const
{
    double delta = calcDelta(rho);
    double tau = calcTau(T);

    return specGasConst*(critT/tau)*(1.0 + tau*(phi0t(delta, tau) + phirt(delta, tau)) + delta*phird(delta, tau));
}

double Helmholtz::cv(double rho, double T) const
{
    double delta = calcDelta(rho);
    double tau = calcTau(T);

    return specGasConst*(-tau*tau*(phi0tt(delta, tau) + phirtt(delta, tau)));
}

double Helmholtz::cp(double rho, double T) const
{
    double delta = calcDelta(rho);
    double tau = calcTau(T);

    return specGasConst*(-tau*tau*(phi0tt(delta, tau) + phirtt(delta, tau)) + std::pow(1.0 + delta*phird(delta, tau) - delta*tau*phirdt(delta, tau), 2)/(1 + 2*delta*phird(delta, tau) + delta*delta*phirdd(delta, tau)));
}

double Helmholtz::a2(double rho, double T) const
{
    double delta = calcDelta(rho);
    double tau = calcTau(T);

    double phirdAux = phird(delta, tau);
    return specGasConst*T*(1.0 + 2.0*delta*phirdAux + delta*delta*phirdd(delta, tau) - pow(1.0 + delta*(phirdAux - tau*phirdt(delta, tau)), 2)/(tau*tau*(phi0tt(delta, tau) + phirtt(delta, tau))));
}

double Helmholtz::a(double rho, double T) const
{
    return std::sqrt(a2(rho, T));
}

/*double Helmholtz::vaporPressure(double T) const
{
    double theta = calcTheta(T);
    return critP*std::exp(calcTau(T)*(satPCoeffs[0]*theta
                                    + satPCoeffs[1]*std::pow(theta, 1.5)
                                    + satPCoeffs[2]*std::pow(theta, 3)
                                    + satPCoeffs[3]*std::pow(theta, 3.5)
                                    + satPCoeffs[4]*std::pow(theta, 4)
                                    + satPCoeffs[5]*std::pow(theta, 7.5)));
}

double Helmholtz::saturatedVaporDensity(double T) const
{
    double theta = calcTheta(T);
    return critRho*(1.0 + satLiqRhoCoeffs[0]*std::pow(theta, 1.0/3.0)
                        + satLiqRhoCoeffs[1]*std::pow(theta, 2.0/3.0)
                        + satLiqRhoCoeffs[2]*std::pow(theta, 5.0/3.0)
                        + satLiqRhoCoeffs[3]*std::pow(theta, 16.0/3.0)
                        + satLiqRhoCoeffs[4]*std::pow(theta, 43.0/3.0)
                        + satLiqRhoCoeffs[5]*std::pow(theta, 110.0/3.0));
}

double Helmholtz::saturatedLiquidDensity(double T) const
{
    double theta = calcTheta(T);
    return critRho*std::exp(satVapRhoCoeffs[0]*std::pow(theta, 2.0/6.0)
                          + satVapRhoCoeffs[1]*std::pow(theta, 4.0/6.0)
                          + satVapRhoCoeffs[2]*std::pow(theta, 8.0/6.0)
                          + satVapRhoCoeffs[3]*std::pow(theta, 18.0/6.0)
                          + satVapRhoCoeffs[4]*std::pow(theta, 37.0/6.0)
                          + satVapRhoCoeffs[5]*std::pow(theta, 71.0/6.0));
}*/

double Helmholtz::tFromRhoP(double rho, double p, double guessT) const
{
    double delta = calcDelta(rho);

    double tau = nonLinearSolver.solve([=](double val) { return pFunc(delta, val) - p; },
                                       [=](double val) { return pDTauFunc(delta, val); },
                                       critT/guessT);
    
    return critT/tau;
}

double Helmholtz::tFromRhoE(double rho, double e, double guessT) const
{
    double delta = calcDelta(rho);

    double tau = nonLinearSolver.solve([=](double val) { return eFunc(delta, val) - e; },
                                       [=](double val) { return eDTauFunc(delta, val);},
                                       critT/guessT);
    
    return critT/tau;
}

double Helmholtz::tFromRhoS(double rho, double s, double guessT) const
{
    double delta = calcDelta(rho);

    double tau = nonLinearSolver.solve([=](double val) { return sFunc(delta, val) - s; },
                                       [=](double val) { return sDTauFunc(delta, val); },
                                       critT/guessT);
    
    return critT/tau;
}

double Helmholtz::tFromRhoH(double T, double h, double guessRho) const
{
    double tau = calcTau(T);

    double delta = nonLinearSolver.solve([=](double val) { return hFunc(val, tau) - h; },
                                         [=](double val) { return hDDeltaFunc(val, tau); },
                                         guessRho/critRho);

    return delta*critRho;
}

double Helmholtz::rhoFromTP(double T, double p, double guessRho) const
{
    double tau = calcTau(T);

    double delta = nonLinearSolver.solve([=](double val) { return pFunc(val, tau) - p; },
                                         [=](double val) { return pDDeltaFunc(val, tau); },
                                         guessRho/critRho);

    return delta*critRho;
}

double Helmholtz::rhoFromTE(double T, double e, double guessRho) const
{
    double tau = calcTau(T);

    double delta = nonLinearSolver.solve([=](double val) { return eFunc(val, tau) - e; },
                                         [=](double val) { return eDDeltaFunc(val, tau); },
                                         guessRho/critRho);

    return delta*critRho;
}

double Helmholtz::rhoFromTS(double T, double s, double guessRho) const
{
    double tau = calcTau(T);

    double delta = nonLinearSolver.solve([=](double val) { return sFunc(val, tau) - s; },
                                         [=](double val) { return sDDeltaFunc(val, tau); },
                                         guessRho/critRho);

    return delta*critRho;
}

double Helmholtz::rhoFromTH(double T, double h, double guessRho) const
{
    double tau = calcTau(T);

    double delta = nonLinearSolver.solve([=](double val) { return hFunc(val, tau) - h; },
                                         [=](double val) { return hDDeltaFunc(val, tau); },
                                         guessRho/critRho);

    return delta*critRho;
}


std::pair<double, double> Helmholtz::RhoTFromSP(double s, double p, double guessRho, double guessT) const
{
    std::pair<double, double> result = nonLinearSolver.solve([=](double val1, double val2) { return sFunc(val1, val2) - s; },
                                                             [=](double val1, double val2) { return sDDeltaFunc(val1, val2); },
                                                             [=](double val1, double val2) { return sDTauFunc(val1, val2); },
                                                             [=](double val1, double val2) { return pFunc(val1, val2) - p; },
                                                             [=](double val1, double val2) { return pDDeltaFunc(val1, val2); },
                                                             [=](double val1, double val2) { return pDTauFunc(val1, val2); },
                                                             guessRho/critRho,
                                                             critT/guessT);

    return std::make_pair(result.first*critRho, critT/result.second);
}


double Helmholtz::pFunc(double delta, double tau) const
{
    return critT*critRho*specGasConst*(delta/tau)*(1.0 + delta*phird(delta, tau));
}

double Helmholtz::eFunc(double delta, double tau) const
{
    return specGasConst*critT*(phi0t(delta, tau) + phirt(delta, tau));
}

double Helmholtz::sFunc(double delta, double tau) const
{
    return specGasConst*(tau*(phi0t(delta, tau) + phirt(delta, tau)) - phi0(delta, tau) - phir(delta, tau));
}

double Helmholtz::hFunc(double delta, double tau) const
{
    return specGasConst*critT*tau*(1.0 + tau*(phi0t(delta, tau) + phirt(delta, tau)) + delta*phird(delta, tau));
}

double Helmholtz::pDDeltaFunc(double delta, double tau) const
{
    return critRho*critT*(specGasConst/tau)*(1.0 + 2.0*delta*phird(delta, tau) + delta*delta*phirdd(delta, tau));
}

double Helmholtz::eDDeltaFunc(double delta, double tau) const
{
    return 0.0;
}

double Helmholtz::sDDeltaFunc(double delta, double tau) const
{
    return specGasConst*(tau*(phi0dt(delta, tau) + phirdt(delta, tau)) - phi0d(delta, tau) - phird(delta, tau));
}

double Helmholtz::hDDeltaFunc(double delta, double tau) const
{
    return 0.0;
}

double Helmholtz::pDTauFunc(double delta, double tau) const
{
    return specGasConst*critRho*critT*(delta/tau)*((-1.0/tau) - (delta/tau)*phird(delta, tau) + delta*phirdt(delta, tau));
}

double Helmholtz::eDTauFunc(double delta, double tau) const
{
    return specGasConst*critT*(phi0tt(delta, tau) + phirtt(delta, tau));
}

double Helmholtz::sDTauFunc(double delta, double tau) const
{
    return specGasConst*tau*(phi0tt(delta, tau) + phirtt(delta, tau));
}

double Helmholtz::hDTauFunc(double delta, double tau) const
{
    return 0.0;
}