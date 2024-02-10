#include <iostream>
#include<cmath>

/*
Mach is a function that takes fluid flow velocity and the 
properties of air to find the Mach number of the flow. Mach 
number is the speed relative to the speed of sound.
*/

double Mach(double u, double gamma, double R, double T) {

    double a = sqrt(gamma * R * T); // Speed of sound
    
    return u / a; // Return mach number
}