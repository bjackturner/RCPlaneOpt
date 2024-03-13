#include <iostream>
#include<cmath>

/*
Flow Props include functions that are related to the fluid flow such as mach number and Reynolds.
*/

// Mach number function
double Mach(double u, double gamma, double R, double T) {
    
    return u / sqrt(gamma * R * T); // Return mach number
}

// Reynolds number function
double Reynolds(double rho, double u, double l, double mu) {

    return (rho * u * l) / mu; // Return reynolds number
}