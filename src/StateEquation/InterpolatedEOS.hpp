#ifndef INTERPOLATEDEOS
#define INTERPOLATEDEOS

#include <array>
#include <vector>

#include "NonLinearSolver/NewtonMethod.hpp"
#include "Interpolation/Interpolation.hpp"

template <typename INTERPOLATION>
class InterpolatedEOS
{
    public:

        InterpolatedEOS()
        {
            entropyInterpolation = INTERPOLATION(std::vector<int>({200}),
                                                 std::vector<int>({50, 100, 100}),
                                                 std::vector<double>({1/1188.87, 1/0.008}),
                                                 std::vector<double>({1800000.0, 2005000.0, 2650000.0, 4085000.27}),
                                                 Interpolation::Transformation::LOGINV,
                                                 Interpolation::Transformation::NONE);
        }

        InterpolatedEOS(INTERPOLATION _entropyInterpolation)
        {
            entropyInterpolation = _entropyInterpolation;
        }

        void init();

        double s(double rho, double e) const
        {
            return entropyInterpolation.calc(rho, e);
        }

        double p(double rho, double e) const
        {
            return - T(rho, e)*rho*rho*sDr(rho, e);
        }

        double T(double rho, double e) const
        {
            return 1.0/sDe(rho, e);
        }

        double h(double rho, double e) const
        {
            return e + p(rho, e)/rho;
        }

        double cv(double rho, double T) const
        {
            return -1.0/(T*sDee(rho, e));
        }

        double cp(double rho, double T) const
        {
            return 0.0 //TODO
        }

        double a2(double rho, double e) const
        {
            double pVal = p(rho, e);
            double TVal = T(rho, e);

            return -2.0*pVal*TVal*sDre(rho, e) - ((pVal*pVal*TVal)/(rho*rho))*sDee(rho, e) - TVal*rho*rho*sDrr(rho, e) + (2.0*pVal)/rho;
        }

        double a(double rho, double e) const
        {
            return std::sqrt(a2(rho, e));
        }

        double eFromRhoP(double rho, double p, double guessE) const;
        //double eFromRhoT(double rho, double T, double guessE) const;
        //double eFromRhoS(double rho, double s, double guessE) const;
        //double eFromRhoH(double rho, double h, double guessE) const;

        //double rhoFromEP(double e, double p, double guessRho) const;
        //double rhoFromET(double e, double T, double guessRho) const;
        //double rhoFromES(double e, double s, double guessRho) const;
        //double rhoFromEH(double e, double h, double guessRho) const;

        std::pair<double, double> RhoTFromSP(double s, double p, double guessRho, double guessE) const;


    private:
        double critT;
        double critRho;
        double specGasConst;

        INTERPOLATION entropyInterpolation;

        double sDr(double rho, double e) const
        {
            return entropyInterpolation.calcDiffX(rho, e);
        }

        double sDrr(double rho, double e) const
        {
            return entropyInterpolation.calcDiff2X(rho, e);
        }

        double sDe(double rho, double e) const
        {
            return entropyInterpolation.calcDiffY(rho, e);
        }

        double sDee(double rho, double e) const
        {
            return entropyInterpolation.calcDiff2Y(rho, e);
        }

        double sDre(double rho, double e) const
        {
            return entropyInterpolation.calcDiffXY(rho, e);
        }
};

template <typename INTERPOLATION>
double InterpolatedEOS<INTERPOLATION>::eFromRhoP(double rho, double p, double guessE) const
{
    double guessE = 0.0;

    return entropyInterpolation.calcYForDiffXbyDiffYEqualConst(rho, guessE, -p/(rho*rho));
}

template <typename INTERPOLATION>
std::pair<double, double> InterpolatedEOS<INTERPOLATION>::RhoTFromSP(double s, double p, double guessRho, double guessE) const
{
    return 0.0; //TODO
}

#endif // INTERPOLATEDEOS