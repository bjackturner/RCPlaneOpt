#include <iostream>
#include<cmath>

double Reynolds(double rho, double u, double l, double mu);

int main() {
    double Re = Reynolds(1.225, 20, 0.5, 1.7654e-5);

    std::cout << "Re = " << Re << std::endl;

    return 0;
}

double Reynolds(double rho, double u, double l, double mu) {

    return (rho * u * l) / mu; // Return reynolds number
}
