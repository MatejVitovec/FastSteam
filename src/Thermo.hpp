#ifndef THERMO_HPP
#define THERMO_HPP

#include <array>

class Thermo
{
    public:
    
        Thermo() {}

        virtual ~Thermo() {}

        //p T a
        std::array<double, 3> internalState(double rho, double e, double TOld, double pOld);

        //rho e T a
        std::array<double, 4> inlet(double pIn, double s0, double TOld, double pOld, double rhoOld);

        //e T a
        std::array<double, 3> outlet(double rhoIn, double p2, double TOld);


        
};

#endif // THERMO_HPP