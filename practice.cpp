#include <iostream>
#include <vector>
#include <array>
#include <fstream>
#include <cmath>
#include <iomanip>

struct Body
{
    std::string name;
    double mass;
    std::array<double, 3> position;
    std::array<double, 3> velocity;
    std::array<double, 3> acceleration;

};


std::array<double,3> subtract( const std::array<double,3> &a , const std::array<double,3> &b){
    return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
}

std::array<double,3> add( const std::array<double,3> &a , const std::array<double,3> &b){
    return {a[0] + b[0], a[1] + b[1], a[2] + b[2]};
}

std::array<double ,3> scale(std::array<double,3> v, double s){
    return {v[0] * s, v[1] * s, v[2] * s};
}

double norm(const std::array<double,3> &v){
    return std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

const double G = 1.0;
const double dt = 0.001;
const double total_time = 20.0;
const int num_steps = total_time / dt;


std::array<std::array<double, 3>, 3> calculate_acceleration(const std::array<Body, 3>& bodies) {
    std::array<std::array<double, 3>, 3> acc = {{{0,0,0},{0,0,0},{0,0,0}}};

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            if (i == j) continue;

            auto r_ij = subtract(bodies[j].position, bodies[i].position);
            double r = norm(r_ij);
            double factor = G * bodies[j].mass / (r * r * r);
            auto contrib = scale(r_ij, factor);
            acc[i] = add(acc[i], contrib);
        }
    }
    return acc;
}

void velocity_verlet_step(std::array<Body, 3>& bodies, double dt) {
    
    // Step 1: DRIFT
    auto acc1 = calculate_acceleration(bodies);

    // Step 2: DRIFT 
    for (int i = 0; i < 3; i++)
        bodies[i].position = add(bodies[i].position,
                              add(scale(bodies[i].velocity, dt),
                                  scale(acc1[i], 0.5 * dt * dt)));

    // Step 3: KICK 
    auto acc2 = calculate_acceleration(bodies);

    // Step 4: DRIFT
    for (int i = 0; i < 3; i++)
        bodies[i].velocity = add(bodies[i].velocity,
                              scale(add(acc1[i], acc2[i]), 0.5 * dt));
}

std::array<Body, 3> initialize_bodies(std::string scenario) {
    std::array<Body, 3> bodies;

    if (scenario == "figure-eight") {
        bodies[0] = {"Body 1", 1.0, {-0.97000436,  0.24308753, 0.0}, { 0.466203685,  0.43236573, 0.0}, {0,0,0}};
        bodies[1] = {"Body 2", 1.0, { 0.0,          0.0,        0.0}, {-0.93240737,  -0.86473146, 0.0}, {0,0,0}};
        bodies[2] = {"Body 3", 1.0, { 0.97000436,  -0.24308753, 0.0}, { 0.466203685,  0.43236573, 0.0}, {0,0,0}};
    }   
    else if (scenario == "euler") {
   
    bodies[0] = {"Body 1", 1.0, {-1.0, 0.0, 0.0}, {0.0,  0.5, 0.0}, {0,0,0}};
    bodies[1] = {"Body 2", 1.0, { 0.0, 0.0, 0.0}, {0.0,  0.0, 0.0}, {0,0,0}};
    bodies[2] = {"Body 3", 1.0, { 1.0, 0.0, 0.0}, {0.0, -0.5, 0.0}, {0,0,0}};
    }


    else if (scenario == "lagrange") {
    // Montgomery 2001 - exact stable equilateral triangle
    double v = 0.5 * std::sqrt(3.0 * G);
    bodies[0] = {"Body 1", 1.0, { 1.0,  0.0, 0.0}, { 0.0,  v, 0.0}, {0,0,0}};
    bodies[1] = {"Body 2", 1.0, {-0.5,  0.8660254, 0.0}, {-v*0.8660254,  -v*0.5, 0.0}, {0,0,0}};
    bodies[2] = {"Body 3", 1.0, {-0.5, -0.8660254, 0.0}, { v*0.8660254,  -v*0.5, 0.0}, {0,0,0}};
    }   
    else if (scenario == "chaotic") {
    bodies[0] = {"Body 1", 1.0, {-0.97000436,  0.24308753,  0.0 }, { 0.466203685,  0.43236573,  0.08}, {0,0,0}};
    bodies[1] = {"Body 2", 1.0, { 0.0,          0.0,         0.06}, {-0.93240737,  -0.86473146,  0.04}, {0,0,0}};
    bodies[2] = {"Body 3", 1.0, { 0.97000436,  -0.24308753,  0.0 }, { 0.506203685,  0.47236573,  0.0 }, {0,0,0}};
    }
    return bodies;
}

int main() {
    std::vector<std::string> scenarios = {"figure-eight", "lagrange", "chaotic", "euler"};

    for (auto& scenario : scenarios) {
        auto bodies = initialize_bodies(scenario);

        std::ofstream file(scenario + ".csv");
        file << "x0,y0,z0,x1,y1,z1,x2,y2,z2\n";

        
        int steps = (scenario == "lagrange") ? 500000 :
            (scenario == "chaotic")  ?  40000 :
            (scenario == "euler")    ? 200000 : 50000;

        for (int step = 0; step < steps; step++) {
            velocity_verlet_step(bodies, dt);
            if (step % 10 == 0) {
                file << bodies[0].position[0] << "," << bodies[0].position[1] << "," << bodies[0].position[2] << ","
                     << bodies[1].position[0] << "," << bodies[1].position[1] << "," << bodies[1].position[2] << ","
                     << bodies[2].position[0] << "," << bodies[2].position[1] << "," << bodies[2].position[2] << "\n";
            }
        }
        file.close();
        std::cout << scenario << " done!\n";
    }
    return 0;
}