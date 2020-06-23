/*
 *  Monolithic Collapse: N-body code
 *
 *  Copyright © 2017 Trevor Picard. All rights reserved.
 *
 */

///////////////////////header//////////////////////////////////
#include <iostream>
#include <cmath>
#include <math.h>
#include <random>
#include <stdlib.h>
#include <fstream>
#include <algorithm>
//#include <omp.h>
using namespace std;



///////////////////////global variables//////////////////////////////////
const int n_dim = 3;                 // spatial dimensions (x, y, z)
const double G = 6.67408e-11;        // gravitational constant
const double pc_m = 3.085677581e+16; // pc to m conversion
const double sol_kg = 1.989e+30;     // solar masses to kg conversion

// physical parameters of the simulation
const int N = 300;                         // number of particles
const double M = sol_kg * 1e+12;           // total mass of the system [kg]
const double Ri0 = 4000 * pc_m;            // initial radius of the system [m]
double vol = (4 / 3) * M_PI * pow(Ri0, 3); // initial volume of the system [m^3]
const double vi0 = 10000 / sqrt(n_dim);     // initial velocity [m/s]
const double epsilon = 0.001;              // force softening parameter

// times and time-dependent parameters
double t = 0;                                      // initial time
const double t0 = pow(G * M / pow(Ri0, 3), -0.5);  // time scale
const int n_dyn = 25;                              // number of crossing times
const int steps_dyn = 30000;                       // steps per crossing time
const double t_dyn = Ri0 / (vi0 * t0);             // dynamical timescale
const double t_relax = (N * t_dyn) / (8 * log(N)); // relaxation timescale
const double dt = t_dyn / steps_dyn;               // time step
const double t_max = n_dyn * t_dyn;                // end time


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
} p[N];



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
                        int index) {
    
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
 * Requires: a particle array.
 * Modifies: particle accelerations.
 * Effects : resets all acceleration components to 0.
 * Used In : acceleration()
 */
void reset_g() {

    // loops through particles
    #pragma omp parallel for
    for (int i = 0; i < N; ++i) {
        
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
 * Used In : leap_frog()
 */
void acceleration() {

    // sets the acceleration components to 0
    reset_g();
    
    //cout << "finding the center of mass..." << endl << endl;
    
    if (t > 0) {
        // loops through each particle and finds the center of mass of the system
        double cm[n_dim] = {0, 0, 0};
        #pragma omp parallel for
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < n_dim; ++j) {
                cm[j] += p[i].r[j] / N;
            }
        }
        
        
        //cout << "updating positions..." << endl << endl;
        
        // update the positions relative to the center of mass
        #pragma omp parallel for
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < n_dim; ++j) {
                p[i].r[j] -= cm[j];
            }
        }
    }

    

    //cout << "calculating accelerations..." << endl << endl;
    
    // for each particle i...
    #pragma omp parallel for
    for (int i = 0; i < N; ++i) {
        
        
        // sum accelerations from every other particle j...
        //#pragma omp parallel for num_threads(16)
        for (int j = 0; j < N; ++j) {
            
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
                if (std::isnan(gx)) {
                    gx = 0;
                }
                if (std::isnan(gy)) {
                    gy = 0;
                }
                if (std::isnan(gz)) {
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
 * Modifies: the density & velocity output file streams
 * Effects : calculates and returns the kinetic, potiential, ratio of kinetic 
 *           to potential, and/or total energy of the system.
 * Used In : leap_frog()
 */
void density_v_rms(ofstream &dfile, ofstream &vfile) {

    // defines a volume density array
    struct pair {
        double r;
        double v;
    };

    pair r_vrms[N] = {};
    double rad[N]  = {};
    
    
    cout << "determining density & rms velocity profiles...";
        
    // calculates the corresponding spherical volume at each particle's radius
    #pragma omp parallel for
    for (int i = 0; i < N; ++i) {
        r_vrms[i].r = sqrt(pow(p[i].r[0], 2) +
                           pow(p[i].r[1], 2) +
                           pow(p[i].r[2], 2));
        
        rad[i] = r_vrms[i].r;
        
        r_vrms[i].v = sqrt(pow(p[i].v[0], 2) +
                           pow(p[i].v[1], 2) +
                           pow(p[i].v[2], 2));
    }

    
    // sorts the array by radius
    sort(rad, rad + N);
    
    for (int i = 0; i < N; ++i) {
        vfile << r_vrms[i].r * Ri0 / t0 / 1000 << " "
              << r_vrms[i].v * Ri0 / pc_m / 1000 << endl;
    }
    vfile << endl;

 
    // calculates the density of particles starting at the center
    const int n_bins = 10;
    int bin_count = 0;
    double  rho[n_bins] = {};
    double r_step = rad[N - 1] / n_bins;
    double r_in = 0;
    double r_out = r_step;
    double V_in = 0;
    double V_out = 0;
    int i = 0;
    while (bin_count < n_bins) {
        
        // calculates the total volume enclosed within the shell's outer radius
        V_out = (4 / 3) * M_PI * pow(r_out, 3);

        // sums the masses within the shell volume
        double m_tot = 0;
        while ((rad[i] <= r_out) && (i < N)) {
            m_tot += p[i].m;
            ++i;
        }
        
        rho[bin_count] = m_tot / (V_out - V_in);

        // moves to the next volume shell and writes the density
        r_in = r_out;
        r_out += r_step;
        V_in = V_out;
        if (m_tot != 0) {
            
            dfile << (rho[bin_count] * M * pow(pc_m / Ri0, 3)) / sol_kg << " "
                  << rad[i] * Ri0 / (pc_m * 1000) << endl;
        }
        
        ++bin_count;
    }
    dfile << endl;
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
void generate_sphere() {

    
    default_random_engine generator;
    uniform_real_distribution<double> distribution(-1, 1);
    
    
    for (int i = 0; i < N; ++i) {
        
        bool done = false;
        
        while (!done) {
            
            double x0 = distribution(generator);
            double y0 = distribution(generator);
            double z0 = distribution(generator);
            
            bool in_sphere = (sqrt(pow(x0, 2) + pow(y0, 2) + pow(z0, 2)) <= 1);
            
            if (in_sphere) {
                initial_conditions(1. / N,                          // m
                                   x0,                              // x_0
                                   y0,                              // y_0
                                   z0,                              // z_0
                                   0,                               // v_x,0
                                   0,                               // v_y,0
                                   0,                               // v_z,0
                                   i);
                done = true;
            }
        }
        
        // assigns Gaussian random velocities
        double sig = 0.1 * vi0 * t0 / Ri0;
        for (int j = 0; j < n_dim; ++j) {
            p[i].v[j] = randn(vi0 * t0 / Ri0, sig);
        }
    }

    cout << "initial configuration set." << endl << endl;
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
void leap_frog(ofstream &dfile, ofstream &vfile, ofstream &pfile) {

    // counts the steps, # of crossing times, and half-mass radii
    int step_count = 0;
    int dyn_count  = 0;
    int rho_v_step = 5;
    
    bool pending_relaxed = true;
    bool print_one = true;
    bool relax_print = false;
    
    // loops through each time step
    while (t <= t_max) {
        
      //cout << step_count << endl;        

        if (pending_relaxed && (t >= t_relax)) {
            cout << endl << "aaaand the system is relaxed." << endl << endl;
            pending_relaxed = false;
            print_one = true;
            relax_print = true;
        }
        
        // finds accelerations of all particles
        acceleration();
        
        
        // keeps track of when to record position data
        bool snapshot = (step_count % 5 == 0) || (print_one && ((t == 0) || ((!pending_relaxed &&
									      (dyn_count % rho_v_step == 0)))));
        
        // loops through each particle
        for (int i = 0; i < N; ++i) {
            
            // writes requested snapshot data to the file
            if (snapshot || relax_print) {
                
                // loops through and writes x, y, z dimensions to the file
                for (int j = 0; j < n_dim; ++j) {
                    pfile << (p[i].r[j] * Ri0) / (pc_m * 1000) << " ";
                    
                    // starts a new line once all components are written
                    if (j == n_dim - 1) {
                        pfile << endl;
                    }
                }
            }
            
            // updates positions/velocities for the next time step
	        #pragma omp parallel for
            for (int j = 0; j < n_dim; ++j) {
                
                // updates velocity components only after t_0
                if (t > 0) {
                    p[i].v[j] += p[i].g[j] * dt;
                }
                p[i].r[j] += p[i].v[j] * dt;
            }
        }
        
        // skips a line after each recorded snapshot
        if (snapshot || relax_print) {
            
            // writes the density & v_rms profiles to files
            density_v_rms(dfile, vfile);
            cout << "done." << endl << endl;
            pfile << endl;
            
            print_one = false;
            relax_print = false;
        }
        
        // moves to the next time step
        t += dt;
        ++step_count;
        
        // counts # of crossing times
        if (step_count % steps_dyn == 0) {
            ++dyn_count;
            print_one = true;
            cout << endl << dyn_count << " dynamical times complete."
                 << endl << endl;
        }
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
void check_and_run(ofstream &dfile, ofstream &vfile, ofstream &pfile) {
    
    
    // checks if the file was opened correctly
    if (dfile.is_open() && vfile.is_open() && pfile.is_open()) {
        
        leap_frog(dfile, vfile, pfile);
        
        // closes the file and resets the time
        dfile.close();
        vfile.close();
        pfile.close();
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


/*
 * random spherical distribution
 *
 */
  //int N_threads = 64;
  //omp_set_num_threads(N_threads);
    
    generate_sphere();
     
    // runs the euler method and stores position data
    ofstream density;
    ofstream v_rms;
    ofstream particles;
    density.open("density.txt");
    v_rms.open("v_rms.txt");
    particles.open("particles.txt");
    
    check_and_run(density, v_rms, particles);
    
}
