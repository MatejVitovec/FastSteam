#ifndef CONSISTENTEOS_HPP
#define CONSISTENTEOS_HPP

#include <array>

template <typename INTERPOLATION>
class ConsistentEoS
{
    public:    
        ConsistentEoS() {}


    private:
        INTERPOLATION entropyInterpolation;

        
};

#endif // CONSISTENTEOS_HPP