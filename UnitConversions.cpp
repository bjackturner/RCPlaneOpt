#include <iostream>
#include <cmath>

double mag(double x[]) {
    
    double magnitude = 0;

    for (int i = 0; i < std::size(x); i++) {
        magnitude += x[i] * x[i];

    }

    return sqrt(magnitude); // Return mach number
}

