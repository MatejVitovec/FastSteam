#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <memory>
#include <chrono>
#include <fenv.h>



int main(int argc, char** argv)
{
    feenableexcept(FE_INVALID | FE_OVERFLOW);

    
    return 0;
}