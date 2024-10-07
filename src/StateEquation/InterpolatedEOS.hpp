#ifndef HRUBYEQUATION
#define HRUBYEQUATION

#include <array>
#include <vector>

#include "NonLinearSolver/NewtonMethod.hpp"

template <typename INTERPOLATION>
class HrubyEquation
{
    public:

        HrubyEquation();

        double s(double rho, double e) const;

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
        double eFromRhoT(double rho, double T, double guessE) const;
        double eFromRhoS(double rho, double s, double guessE) const;
        double eFromRhoH(double rho, double h, double guessE) const;

        double rhoFromEP(double e, double p, double guessRho) const;
        double rhoFromET(double e, double T, double guessRho) const;
        double rhoFromES(double e, double s, double guessRho) const;
        double rhoFromEH(double e, double h, double guessRho) const;

        std::pair<double, double> RhoTFromSP(double s, double p, double guessRho, double guessE) const;


    private:
        double critT;
        double critRho;
        double specGasConst;

        INTERPOLATION entropyInterpolation;

        double sDr(double rho, double e) const;
        double sDrr(double rho, double e) const;
        double sDe(double rho, double e) const;
        double sDee(double rho, double e) const;
        double sDre(double rho, double e) const;
};


template <typename INTERPOLATION>
double HrubyEquation<INTERPOLATION>::s(double rho, double e) const
{
    
}

template <typename INTERPOLATION>
double HrubyEquation<INTERPOLATION>::sDr(double rho, double e) const
{
    
}

template <typename INTERPOLATION>
double HrubyEquation<INTERPOLATION>::sDrr(double rho, double e) const
{
    
}

template <typename INTERPOLATION>
double HrubyEquation<INTERPOLATION>::sDe(double rho, double e) const
{
    
}

template <typename INTERPOLATION>
double HrubyEquation<INTERPOLATION>::sDee(double rho, double e) const
{
    
}

template <typename INTERPOLATION>
double HrubyEquation<INTERPOLATION>::sDre(double rho, double e) const
{
    
}


#endif // HRUBYEQUATION