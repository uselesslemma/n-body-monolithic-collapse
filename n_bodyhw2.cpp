/*
*  N-body code
*
*  Copyright © 2017 Trevor Picard. All rights reserved.
*/

///////////////////////header//////////////////////////////////
#include <iostream>
#include <cmath>
#include <math.h>
#include <random>
#include <stdlib.h>
#include <fstream>
#include <algorithm>
using namespace std;



///////////////////////global variables//////////////////////////////////
const int n_dim = 3;      // number of spatial dimensions (x, y, z)
double t = 0;             // initial time
double dt = 0.0001;       // initial time step
double t_max = 100;       // initial end time
int n_orbits = 10;        // initial number of orbits
int steps_orb = 1000;     // initial number of time steps per orbit
double epsilon = 0;       // initial softening



/////////////////////////////////////////////////////////////////////
/*
 * Requires: a mean velocity and a standard deviation.
 * Modifies: nothing.
 * Effects : returns a Gaussian random velocity given the mean velocity and
 *           standard deviation of a velocity distribution.
 * Used In : generate_sphere()
 */
double randn(double mu, double sigma) {
    
    double U1, U2, W, mult;
    static double X1, X2;
    static int call = 0;
    
    if (call == 1) {
        
        call = !call;
        return (mu + sigma * (double) X2);
    }
    
    do {
        U1 = -1 + ((double) rand () / RAND_MAX) * 2;
        U2 = -1 + ((double) rand () / RAND_MAX) * 2;
        W = pow (U1, 2) + pow (U2, 2);
    }
    
    while (W >= 1 || W == 0);
    
    mult = sqrt ((-2 * log (W)) / W);
    X1 = U1 * mult;
    X2 = U2 * mult;
    
    call = !call;
    
    return (mu + sigma * (double) X1);
}



//////////////////////////////////////////////////////////////////////////
/*
 * Requires: nothing.
 * Modifies: nothing.
 * Effects : defines a particle, which is a data structure with a mass
 *           and with 3D position, velocity, and acceleration components.
 */
struct particle {
    double m;
    double r[n_dim];
    double v[n_dim];
    double g[n_dim];
};



//////////////////////////////////////////////////////////////////////////
/*
 * Requires: an array of particles.
 * Modifies: the members of each particle in array p.
 * Effects : assigns masses and 3D initial coniditions to the particles in the
 *           simulation.
 * Used In : main()
 */
void initial_conditions(double m,
                        double rx,
                        double ry,
                        double rz,
                        double vx,
                        double vy,
                        double vz,
                        int index,
                        particle *p) {
    
    p[index].m    = m;
    p[index].r[0] = rx;
    p[index].r[1] = ry;
    p[index].r[2] = rz;
    p[index].v[0] = vx;
    p[index].v[1] = vy;
    p[index].v[2] = vz;
}



//////////////////////////////////////////////////////////////////////////
/*
 * Requires: a particle array and the number of particles.
 * Modifies: particle accelerations.
 * Effects : resets all acceleration components to 0.
 * Used In : acceleration()
 */
void reset_g(particle *p, int num_p) {
    
    // loops through particles
    for (int i = 0; i < num_p; ++i) {
        
        // sets all components to 0
        for (int j = 0; j < n_dim; ++j) {
            p[i].g[j] = 0;
        }
    }
}



//////////////////////////////////////////////////////////////////////////
/*
 * Requires: specified initial positions of the N particles in array p.
 * Modifies: particle accelerations.
 * Effects : updates the total accelerations of each particle due to every
 *           other particle in the simulation.
 * Used In : leap_frog(), euler(), sphere_leap_frog()
 */
void acceleration(particle *p, int num_p) {

    // sets the acceleration components to 0
    reset_g(p, num_p);
    
    // for each particle i...
    for (int i = 0; i < num_p; ++i) {
        
        // sum accelerations from every other particle j...
        for (int j = 0; j < num_p; ++j) {
            
            // but not due to the particle itself.
            if (i != j) {
                
                // distance between particles i & j
                double d  = pow(sqrt(pow(p[j].r[0] - p[i].r[0], 2) +
                                     pow(p[j].r[1] - p[i].r[1], 2) +
                                     pow(p[j].r[2] - p[i].r[2], 2) +
                                     pow(epsilon, 2)), 3);
                
                // acceleration components of particle i due to particle j
                double gx = p[j].m * (p[j].r[0] - p[i].r[0]) / d;
                double gy = p[j].m * (p[j].r[1] - p[i].r[1]) / d;
                double gz = p[j].m * (p[j].r[2] - p[i].r[2]) / d;
                
                // sets NaN accelerations to 0
                if (isnan(gx)) {
                    gx = 0;
                }
                if (isnan(gy)) {
                    gy = 0;
                }
                if (isnan(gz)) {
                    gz = 0;
                }
                
                // adds accelerations to the total for each component
                p[i].g[0] += gx;
                p[i].g[1] += gy;
                p[i].g[2] += gz;
            }
        }
    }
}



//////////////////////////////////////////////////////////////////////////
/*
 * Requires: N particles in array p that have specified masses, positions, & 
 *           velocities, and which type of energy to return.
 * Modifies: nothing.
 * Effects : calculates and returns the kinetic, potiential, ratio of kinetic 
 *           to potential, and/or total energy of the system.
 * Used In : leap_frog(), euler(), sphere_leap_frog()
 */
double energy(particle *p, int num_p, string e_return) {

    
/////////////total kinetic energy////////////////////
    
    double K = 0;
    
    // if the function call asks for K or K/U...
    if (e_return != "U") {
        
        // adds the KE of each particle to the total of the system
        for (int i = 0; i < num_p; ++i) {
            K += p[i].m * sqrt(pow(p[i].v[0], 2) +
                               pow(p[i].v[1], 2) +
                               pow(p[i].v[2], 2));
        }
    }
    
/////////////total potential energy//////////////////
    
    double U = 0;
    
    // if function call asks for U or K/U...
    if (e_return != "K") {
        
        // for each particle i...
        for (int i = 0; i < num_p; ++i) {
            
            // sum potential energy due to every other particle j
            for (int j = i + 1; j < num_p; ++j) {
                
                // PE between particles i & j=
                double U_temp = -(p[i].m * p[j].m) /
                sqrt(pow(p[j].r[0] - p[i].r[0], 2) +
                     pow(p[j].r[1] - p[i].r[1], 2) +
                     pow(p[j].r[2] - p[i].r[2], 2) +
                     pow(epsilon, 2));
                
                // set NaN PEs to 0
                if (isnan(U_temp)) {
                    U_temp = 0;
                }
                
                // adds PE to the total of the system
                U += U_temp;
            }
        }
    }
    
    // returns the energy component specified in the function call
    if (e_return == "U") {
        return U;
    }
    else if (e_return == "K") {
        return K;
    }
    else if (e_return == "K/U") {
        return K / U;
    }
    else if (e_return == "total") {
        return K + U;
    }
    
    // returns 0 for erroneous function calls
    else {
        return 0;
    }
}



//////////////////////////////////////////////////////////////////////////
/*
 * Requires: an open output file stream, specified initial conditions of the
 *           N particles in array p, and how many orbits should be stored.
 * Modifies: the output file stream.
 * Effects : evolves the system of particles in time using 2D leap-frog
 *           integration, keeping track of particle properties and the
 *           evolution of the system's energy.
 * Used In : check_and_run()
 */
void leap_frog(particle *p, int num_p, ofstream &file, int o_step) {
    
    // defines total energies of adjacent time steps and their difference
    double E1 = 0;
    double E2 = 0;
    double dE = 0;
    
    // defines the # of steps and orbits
    int step_count = 0;
    int orb_count  = 0;
    
    // loops through each time step
    while (t <= t_max) {
    
        // finds accelerations of all particles
        acceleration(p, num_p);
        
        
        // calculates the current total energy of the system
        if (t > 0) {
            E2 = energy(p, num_p, "total");
            dE = abs(E2 - E1);
        }
        
        // calculates the current total energy at the previous time step (> t_0)
        E1 = energy(p, num_p, "total");

    
        // loops through the particles
        for (int i = 0; i < num_p; ++i) {
        
            // writes current positions to file every 'o_step' orbits
            if (orb_count % o_step == 0) {
                file << p[i].r[0] << " " << p[i].r[1] << " ";
                
                // writes energy data once all particle data is written
                if (i == num_p - 1) {
                    file << (t / t_max) * n_orbits << " " << dE / E1 <<
                     " " << E1 << " " << endl;
                }
            }
        
            // updates positions/velocities for the next time step
            for (int j = 0; j < n_dim; ++j) {
                
                // updates velocity components only after t_0
                if (t > 0) {
                    p[i].v[j] += p[i].g[j] * dt;
                }
                p[i].r[j] += p[i].v[j] * dt;
            }
        }
    
        // moves to the next time step
        t += dt;
        ++step_count;
        
        // counts and prints # of orbits
        if (step_count % steps_orb == 0) {
            ++orb_count;
            cout << orb_count << endl;
        }
    }
}



//////////////////////////////////////////////////////////////////////////
/*
 * Requires: an open output file stream, specified initial conditions of the
 *           N particles in array p, and how many orbits should be stored.
 * Modifies: the output file stream.
 * Effects : evolves the system of particles in time using 2D Euler
 *           integration, keeping track of particle properties and the
 *           evolution of the system's energy.
 * Used In : check_and_run()
 */
void euler(particle *p, int num_p, ofstream &file, int o_step) {
    
    // defines total energies of adjacent time steps and their difference
    double E1 = 0;
    double E2 = 0;
    double dE = 0;
    
    // defines the # of steps and orbits
    int step_count = 0;
    int orb_count  = 0;
    
    // loops through each time step
    while (t <= t_max) {
        
        // finds accelerations of all particles
        acceleration(p, num_p);
        
        // calculates the current total energy of the system
        if (t > 0) {
            E2 = energy(p, num_p, "total");
            dE = abs(E2 - E1);
        }
        
        // calculates the current total energy at the previous time step (> t_0)
        E1 = energy(p, num_p, "total");

        
        // prints initial system energy
        if (t == 0) {
            cout << E1 << " ";
        }

        // loops through the particles
        for (int i = 0; i < num_p; ++i) {
            
            
            // writes current positions to file every 'o_step' orbits
            if (orb_count % o_step == 0) {
                file << p[i].r[0] << " " << p[i].r[1] << " ";
                
                // writes energy data once all particle data is written
                if (i == num_p - 1) {
                    file << (t / t_max) * n_orbits << " " << dE / E1
                         << " " << E1 << endl;
                }
            }

            
            // updates positions/velocities for the next time step
            for (int j = 0; j < n_dim; ++j) {
                p[i].r[j] += p[i].v[j] * dt + 0.5 * p[i].g[j] * pow(dt, 2);
                p[i].v[j] += p[i].g[j] * dt;
            }
        }
        
        // moves to the next time step
        t += dt;
        ++step_count;
        
        // counts and prints # of orbits
        if (step_count % n_orbits == 0) {
            ++orb_count;
            cout << orb_count << endl;
        }
    }
    
    // prints final system energy
    cout << E1 << endl;
}



//////////////////////////////////////////////////////////////////////////
/*
 * Requires: an output file stream & N particles in array p.
 * Modifies: the output file stream and the initial conditions of the system.
 * Effects : assigns random positions to N particles within a sphere of radius
 *           1, calculates the mean velocity required for an energy ratio of
 *           K/U = -0.1, then generates random velocities using randn().
 * Used In : main()
 */
void generate_sphere(particle *p, int num_p) {
    
    // initial ratio of KE to PE
    double KU = -0.1;
    
    default_random_engine generator;
    uniform_real_distribution<double> distribution(-1.0,1.0);
    
    
    for (int i = 0; i < num_p; ++i) {
        
        bool done = false;
        
        while (!done) {
            
            double x0 = distribution(generator);
            double y0 = distribution(generator);
            double z0 = distribution(generator);
            
            bool in_sphere = (sqrt(pow(x0, 2) + pow(y0, 2) + pow(z0, 2)) <= 1);
            
            if (in_sphere) {
                initial_conditions(1,                               // m
                                   x0,  // x_0
                                   y0,  // y_0
                                   z0,  // z_0
                                   0,                               // v_x,0
                                   0,                               // v_y,0
                                   0,                               // v_z,0
                                   i,
                                   p);
                done = true;
            }
        }
    }

    
    // finds the total KE & corresponding <v> for all particles
    double K = KU * energy(p, num_p, "U");
    double v_ave = sqrt((2 * K) / (n_dim * num_p));
    
    // assigns Gaussian random velocities using mu = <v> and sigma = 0.2
    double sig = 0.2;
    for (int i = 0; i < num_p; ++i) {
        for (int j = 0; j < n_dim; ++j) {
            p[i].v[j] = randn(v_ave, sig);
        }
    }
    
}



//////////////////////////////////////////////////////////////////////////
/*
 * Requires: an open output file stream & specified initial conditions
 * Modifies: the output file stream
 * Effects : Evolves the system of particles in time using leap-frog
 *           integration, keeping track of particle properties and the total
 *           energy of the system.
 * Used In : check_and_run()
 */
void sphere_leap_frog(particle *p, int num_p,
                      ofstream &efile, ofstream &pfile) {
    
    // defines the volume of a sphere with radius 1
    int R = 1;
    double vol = (4 / 3) * M_PI * pow(R, 3);
    double t_cross = 1 / sqrt(num_p / vol);
    int steps_cross = 1000;
    int t_cross_max = 100;
    t_max = t_cross * t_cross_max;
    dt = t_cross / steps_cross;
    
    // counts the steps, # of crossing times, and half-mass radii
    int step_count   = 0;
    int cross_count  = 0;
    bool one_step = true;
    double radius[num_p];
    
    // loops through each time step
    while (t <= t_max) {

        
        // finds accelerations of all particles
        acceleration(p, num_p);
        
        // writes KE and K/U ratio to the energy file
        efile << energy(p, num_p, "K") << " " << energy(p, num_p, "U") << " "
              << energy(p, num_p, "K/U") << " ";
        
        // keeps track of when to record position data
        int maxKE = 686;
        bool snapshot = (t == 0) || ((cross_count == 5) && one_step) ||
                        (step_count == maxKE) || (cross_count == 100);
        
        // loops through each particle
        for (int i = 0; i < num_p; ++i) {
            
            radius[i] = sqrt(pow(p[i].r[0], 2) + pow(p[i].r[1], 2) +
                             pow(p[i].r[2], 2));
            
            // writes requested snapshot data to the file
            if (snapshot) {
                
                // loops through and writes x, y, z dimensions to the file
                for (int j = 0; j < n_dim; ++j) {
                    pfile << p[i].r[j] << " ";
                    
                    // starts a new line once all components are written
                    if (j == n_dim - 1) {
                        pfile << endl;
                    }
                }
            }
            
            // updates positions/velocities for the next time step
            for (int j = 0; j < n_dim; ++j) {
                
                // updates velocity components only after t_0
                if (t > 0) {
                    p[i].v[j] += p[i].g[j] * dt;
                }
                p[i].r[j] += p[i].v[j] * dt;
            }
        }
        
        // skips a line after each recorded snapshot
        if (snapshot) {
            pfile << endl;
        }
        
        // only allows one snapshot to be recorded at the 5th crossing time
        if (cross_count == 5) {
            one_step = false;
        }
        
        // moves to the next time step
        t += dt;
        ++step_count;
        
        // counts # of crossing times
        if (step_count % steps_cross == 0) {
            ++cross_count;
            cout << cross_count << endl;
        }
        
        // writes the half mass radius to the file
        sort(radius, radius + num_p);
        efile << radius[num_p / 2] << endl;
    }
}



//////////////////////////////////////////////////////////////////////////
/*
 * Requires: an output file stream & an integration method
 * Modifies: the output file stream
 * Effects : Checks if a file was created and opened properly, then runs the
 *           simulation using the specified integration method.
 * Used In : main()
 */
void check_and_run(particle *p, int num_p,
                   ofstream &file, string method, int o_step) {
    
    
    // checks if the file was opened correctly
    if (file.is_open()) {
        
        // chooses a method of integration
        if (method == "leap frog") {
            leap_frog(p, num_p, file, o_step);
        }
        else if (method == "euler") {
            euler(p, num_p, file, o_step);
        }
        else {
            cout << "something's wrong boi..." << endl;
        }
        
        // closes the file and resets the time
        file.close();
        t = 0;
    }
    else {
        cout << "something's wrong boi..." << endl;
    }
}



//////////////////////////////////////////////////////////////////////////
/*
 *
 *
 * BEGIN MAIN PROGRAM
 *
 *
 */
int main() {

    // makes an array of 2 particles
    int N = 2;
    particle p2[N];
    
    /*
     * (1)
     * leap-frog method
     *
     */
    

    double T = 4 * M_PI;
    dt = T / steps_orb;
    t_max = steps_orb * n_orbits * dt;
    double vi = 0.5;
    
    // assigns masses and initial conditions to each particle
    for (int i = 0; i < N; ++i) {
        initial_conditions(1,                          // m
                           pow(-1, i + 2),             // x_0
                           0,                          // y_0
                           0,                          // z_0
                           0,                          // v_x,0
                           vi * pow(-1, i + 2),        // v_y,0
                           0,                          // v_z,0
                           i,
                           p2);
    }
    
    
    
    // runs the leap-frog method and stores position data
    ofstream circles;
    circles.open("circles.txt");
    check_and_run(p2, N, circles, "leap frog", 1);
    
    
//////////////////////////////////////////////////////////////////////////
    /*
     * (2)
     * leap-frog vs. euler method
     *
     */
    
    double r_max = 10;
    double r_min = 1;
    double a = (r_max + r_min);
    T = 2 * M_PI * sqrt(pow(a, 3) / 2);
    double e = (r_max - r_min) / (r_max + r_min);
    int M = 2;
    vi = 0.5 * sqrt(M / a) * sqrt((1 - e) / (1 + e));
    n_orbits = 10000;
    
    
    //////////200 steps per orbit///////////
    
    steps_orb = 200;
    dt = T / steps_orb;
    t_max = steps_orb * n_orbits * dt;
    
    //////////////////leap-frog/////////////////////////////////
    for (int i = 0; i < N; ++i) {
        initial_conditions(1,                          // m
                           10 * pow(-1, i + 2),        // x_0
                           0,                          // y_0
                           0,                          // z_0
                           0,                          // v_x,0
                           vi * pow(-1, i + 2),        // v_y,0
                           0,                          // v_z,0
                           i,
                           p2);
    }
    ofstream lf200;
    lf200.open("lf200.txt");
    check_and_run(p2, N, lf200, "leap frog", 1000);
    
    
    //////////////////euler//////////////////////////////////
    for (int i = 0; i < N; ++i) {
        initial_conditions(1,                          // m
                           10 * pow(-1, i + 2),        // x_0
                           0,                          // y_0
                           0,                          // z_0
                           0,                          // v_x,0
                           vi * pow(-1, i + 2),        // v_y,0
                           0,                          // v_z,0
                           i,
                           p2);
    }
    ofstream euler200;
    euler200.open("euler200.txt");
    check_and_run(p2, N, euler200, "euler", 1);
    
    
    
    
    //////////3000 steps per orbit///////////
    
    steps_orb = 3000;
    dt = T / steps_orb;
    t_max = steps_orb * n_orbits * dt;
    
    // assigns masses and initial conditions to each particle
    for (int i = 0; i < N; ++i) {
        initial_conditions(1,                          // m
                           10 * pow(-1, i + 2),        // x_0
                           0,                          // y_0
                           0,                          // z_0
                           0,                          // v_x,0
                           vi * pow(-1, i + 2),        // v_y,0
                           0,                          // v_z,0
                           i,
                           p2);
    }
    
    // runs the leap-frog method and stores position data
    ofstream lf3000;
    lf3000.open("lf3000.txt");
    check_and_run(p2, N, lf3000, "leap frog", 1000);
    
    
    
    //////////////////euler//////////////////////////////////
    
    // assigns masses and initial conditions to each particle
    for (int i = 0; i < N; ++i) {
        initial_conditions(1,                          // m
                           10 * pow(-1, i + 2),        // x_0
                           0,                          // y_0
                           0,                          // z_0
                           0,                          // v_x,0
                           vi * pow(-1, i + 2),        // v_y,0
                           0,                          // v_z,0
                           i,
                           p2);
    }
    
    // runs the euler method and stores position data
    ofstream euler3000;
    euler3000.open("euler3000.txt");
    check_and_run(p2, N, euler3000, "euler", 1);
    
    
    
//////////////////////////////////////////////////////////////////////////
/*
 * (3)
 * random spherical distribution
 *
 */
    
    N = 500;
    particle p500[N];
    epsilon = 0.001;
    generate_sphere(p500, N);
    
     
    // runs the euler method and stores position data
    ofstream p3energy;
    ofstream p3particles;
    p3energy.open("p3energy.txt");
    p3particles.open("p3particles.txt");
    
    // checks if the file was opened correctly
    if (p3energy.is_open() && p3particles.is_open()) {
        
        sphere_leap_frog(p500, N, p3energy, p3particles);
        p3energy.close();
        p3particles.close();
        t = 0;
    }
    else {
        cout << "something's wrong boi..." << endl;
    }
    
}
