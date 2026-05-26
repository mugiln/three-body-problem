#include <iostream>
#include <fstream>
#include <cmath>

const double G      = 1.0;
const double M_sun  = 1.0;
const double M_jup  = 0.001;   // Jupiter ≈ 1/1000th of the Sun
const double a      = 1.0;     // orbital radius (Sun-Jupiter distance)
const double omega  = 1.0;     // angular velocity of rotating frame

double effective_potential(double x, double y) {
    // Sun sits at (-mu, 0), Jupiter at (1-mu, 0)
    // where mu = M_jup / (M_sun + M_jup)
    double mu  = M_jup / (M_sun + M_jup);

    double x_sun = -mu;          // Sun position on x-axis
    double x_jup = 1.0 - mu;    // Jupiter position on x-axis

    // Distance from (x,y) to each body — Pythagoras in 2D
    double r1 = std::sqrt((x - x_sun)*(x - x_sun) + y*y);
    double r2 = std::sqrt((x - x_jup)*(x - x_jup) + y*y);

    // Effective potential = gravity wells + centrifugal term
    return -(1.0 - mu)/r1 - mu/r2 - 0.5*(x*x + y*y);
}

int main() {
    std::ofstream file("potential.csv");
    file << "x,y,potential\n";

    int N    = 400;          // 400x400 grid = 160,000 points
    double range = 2.0;      // covers x ∈ [-2.0, 2.0], y ∈ [-2.0, 2.0]

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            double x = -range + 2.0 * range * i / (N - 1);
            double y = -range + 2.0 * range * j / (N - 1);

            double phi = effective_potential(x, y);

            // Clamp extreme values near the singularities (very close to Sun/Jupiter)
            if (phi < -10.0) phi = -10.0;

            file << x << "," << y << "," << phi << "\n";
        }
    }

    file.close();
    std::cout << "potential.csv done!\n";
    return 0;
}