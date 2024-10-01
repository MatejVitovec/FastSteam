#ifndef HELMHOLTZ
#define HELMHOLTZ

#include <array>
#include <vector>

#include "NonLinearSolver/NewtonMethod.hpp"

class Helmholtz
{
    public:

        Helmholtz();

        virtual ~Helmholtz() {}

        double p(double rho, double T) const;
        double e(double rho, double T) const;
        double s(double rho, double T) const;
        double h(double rho, double T) const;
        double cv(double rho, double T) const;
        double cp(double rho, double T) const;
        double a2(double rho, double T) const;
        double a(double rho, double T) const;

        double vaporPressure(double T) const;
        double saturatedVaporDensity(double T) const;
        double saturatedLiquidDensity(double T) const;

        double tFromRhoP(double rho, double p, double guessT) const;
        double tFromRhoE(double rho, double e, double guessT) const;
        double tFromRhoS(double rho, double s, double guessT) const;
        double tFromRhoH(double rho, double h, double guessT) const;

        double rhoFromTP(double T, double p, double guessRho) const;
        double rhoFromTE(double T, double e, double guessRho) const;
        double rhoFromTS(double T, double s, double guessRho) const;
        double rhoFromTH(double T, double h, double guessRho) const;

        std::pair<double, double> RhoTFromSP(double s, double p, double guessRho, double guessT) const;

    protected:
        double critT;
        double critRho;
        double specGasConst;
        double critP;

        std::array<double, 6> satPCoeffs = {-7.85951783, 1.84408259, -11.7866497, 22.6807411, -15.9618719, -1.80122502};
        std::array<double, 6> satLiqRhoCoeffs = {1.99274064, 1.09965342, -0.510839303, -1.75493479, -45.5170352, -674694.450};
        std::array<double, 6> satVapRhoCoeffs = {-2.03150240, -2.68302940, -5.38626492, -17.2991605, -44.7586581, -63.9201063};

        double calcDelta(double rho) const;
        double calcTau(double T) const;
        double calcTheta(double T) const;

        NewtonMethod nonLinearSolver;

        virtual double phi0(double delta, double tau) const = 0;
        virtual double phi0d(double delta, double tau) const = 0;
        virtual double phi0dd(double delta, double tau) const = 0;
        virtual double phi0t(double delta, double tau) const = 0;
        virtual double phi0tt(double delta, double tau) const = 0;
        virtual double phi0dt(double delta, double tau) const = 0;

        virtual double phir(double delta, double tau) const = 0;
        virtual double phird(double delta, double tau) const = 0;
        virtual double phirdd(double delta, double tau) const = 0;
        virtual double phirt(double delta, double tau) const = 0;
        virtual double phirtt(double delta, double tau) const = 0;
        virtual double phirdt(double delta, double tau) const = 0;

        double pFunc(double delta, double tau) const;
        double eFunc(double delta, double tau) const;
        double sFunc(double delta, double tau) const;
        double hFunc(double delta, double tau) const;

        double pDDeltaFunc(double delta, double tau) const;
        double eDDeltaFunc(double delta, double tau) const; //TODO
        double sDDeltaFunc(double delta, double tau) const;
        double hDDeltaFunc(double delta, double tau) const; //TODO

        double pDTauFunc(double delta, double tau) const;
        double eDTauFunc(double delta, double tau) const;
        double sDTauFunc(double delta, double tau) const;
        double hDTauFunc(double delta, double tau) const; //TODO

};

#endif // HELMHOLTZ