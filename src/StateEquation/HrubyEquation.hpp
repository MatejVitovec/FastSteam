#ifndef HRUBYEQUATION
#define HRUBYEQUATION

#include <array>
#include <vector>

#include "NonLinearSolver/NewtonMethod.hpp"

class HrubyEquation
{
    public:

        HrubyEquation();

        virtual ~HrubyEquation() {}

        double p(double rho, double e) const;
        double T(double rho, double e) const;
        double s(double rho, double e) const;
        double h(double rho, double e) const;
        double a2(double rho, double e) const;
        double a(double rho, double e) const;

        double eFromRhoP(double rho, double p, double guessE);

        std::pair<double, double> RhoEFromSP(double s, double p, double guessRho, double guessE) const;

    private:
        double critT;
        double critRho;
        double specGasConst;

        std::array<double, 2> beta;
        std::array<double, 10> ceoffC;
        std::array<double, 3> coeffB;

        NewtonMethod nonLinearSolver;

        inline double calcEta(double s) const {return s/specGasConst; }
        inline double calcDelta(double rho) const {return rho/critRho; }
        inline double calcPsi(double e) const {return e/(critT*specGasConst); }
        inline double calcChi(double psi) const {return 1.0/(psi - beta[0]); }
        inline double calcKsi(double psi) const {return beta[1] - calcChi(psi); }
        inline double calcA1(double ksi) const { return coeffB[0] + coeffB[1]/std::pow(ksi, 2) + coeffB[2]/std::pow(ksi, 4); }
        inline double calcA1k(double ksi) const { return -2.0*coeffB[1]/std::pow(ksi, 3) - 4.0*coeffB[2]/std::pow(ksi, 5); }
        inline double calcA1kk(double ksi) const { return 6.0*coeffB[1]/std::pow(ksi, 4) + 20.0*coeffB[2]/std::pow(ksi, 6); }


        double eta0(double delta, double psi) const;
        double eta0d(double delta, double psi) const;
        double eta0dd(double delta, double psi) const;
        double eta0k(double delta, double psi) const;
        double eta0kk(double delta, double psi) const;
        double eta0dk(double delta, double psi) const;

        double etar(double delta, double psi) const;
        double etard(double delta, double psi) const;
        double etardd(double delta, double psi) const;
        double etark(double delta, double psi) const;
        double etarkk(double delta, double psi) const;
        double etardk(double delta, double psi) const;

        double eta(double delta, double psi) const;
        double etad(double delta, double psi) const;
        double etadd(double delta, double psi) const;
        double etap(double delta, double psi) const;
        double etapp(double delta, double psi) const;
        double etadp(double delta, double psi) const;

        double pFunc(double delta, double psi) const;
        double eFunc(double delta, double psi) const;
        double sFunc(double delta, double psi) const;

        double pDDeltaFunc(double delta, double psi) const;
        double sDDeltaFunc(double delta, double psi) const;

        double pDPsiFunc(double delta, double psi) const;
        double eDPsiFunc(double delta, double psi) const;
        double sDPsiFunc(double delta, double psi) const;

};

#endif // HRUBYEQUATION