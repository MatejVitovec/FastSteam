#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <memory>
#include <chrono>
#include <fenv.h>

#include "StateEquation/Iapws95.hpp"
#include "StateEquation/HrubyEquation.hpp"

double auxPressure(double rho, double e)
{
    return (1.33263 - 1.0)*rho*(e - 1995917.326);
}

int main(int argc, char** argv)
{
    feenableexcept(FE_INVALID | FE_OVERFLOW);

    int rhoSize = 401;
    double rhoFirst = 0.00084113 ;
    double rhoLast = 125;

    int eSize = 201;
    double eFirst = 2005000.0;
    double eLast = 2650000.0;

    ////

    std::vector<double> rhoRange(rhoSize); //transformation:  log(1/x)   back - 1.0/exp(x)
    std::vector<double> eRange(eSize);

    double rhoDx = (log(1/rhoLast) - log(1/rhoFirst))/(rhoSize-1);
    for (size_t i = 0; i < rhoRange.size(); i++)
    {
        rhoRange[i] = 1.0/std::exp(log(1/rhoFirst) + i*rhoDx);
    }

    double eDx = (eLast - eFirst)/(eSize-1);
    for (size_t i = 0; i < eRange.size(); i++)
    {
        eRange[i] = eFirst + i*eDx;
    }
    

    Iapws95 i95 = Iapws95();
    HrubyEquation hruby = HrubyEquation();

    std::vector<std::vector<double>> i95Results = std::vector<std::vector<double>>(eRange.size());
    std::vector<std::vector<double>> hrubyresults = std::vector<std::vector<double>>(eRange.size());


    for (size_t j = 0; j < eRange.size(); j++)
    {
        i95Results[j] = std::vector<double>(rhoRange.size());
        hrubyresults[j] = std::vector<double>(rhoRange.size());

        for (size_t i = 0; i < rhoRange.size(); i++)
        {
            double rho = rhoRange[i];
            double e = eRange[j];

            i95Results[j][i] = i95.tFromRhoE(rho, e, auxPressure(rho, e)/(rho*461.51805));
            hrubyresults[j][i] = hruby.T(rho, e);
        }
    }

    
    

    
    return 0;
}