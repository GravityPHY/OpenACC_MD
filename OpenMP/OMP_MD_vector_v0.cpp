/**
 * Author: Michael Graziano mjgrazia@bu.edu
 * Author: Sahan Bandara sahanb@bu.edu
 * Adapted by: Tianyi Xu tyx@bu.edu
 */

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <chrono>
#include <vector>
#include <cmath>

#define R_cut 2.0
#define R_skin 4.0


#define PRINTRESULT 1

using namespace std;

// Function Prototypes:
int main(int argc, char **argv );
void setPositions (int np, int nd, int L, double *pos);
void createNeighborhood (int np, int nd, int L, double *pos, vector<int> *neighbor);
void calcAccel (int np, int nd, int L, double mass, double *pos, double *acc, double &Epot, vector<int> *neighbor);
void boundCond (int np, int nd, int L, double *pos);
void calcVerlet (int np, int nd, double mass, double dt, int L, double &Ekin, double &Epot,
                 double *pos, double *vel, double *acc, vector<int> *neighbor);


int main(int argc, char **argv) {
    /**
     * Main function of the Molecular Dynamics Simulator
     */

    // Variable Declartions
    double *pos;              // Holds the position of the particles
    vector<int> *neighbor;    // Holds the neighborhood table for the simulation
    double *vel;              // Holds the velocity of the particles
    double *acc;              // Holds the accleration of the particles
    int L = 4096;             // length of the simulation box
    double Ekin, Epot;        // Holds the kinetic & potential energy of the system.
    double mass = 1.0;        // Holds the mass of each atom
    double dt = 0.005;        // Holds the time step
    int steps = 5000;         // Holds the maximum number of time steps to do
    int np = 1024;            // Holds the number of particles
    int nd = 2;               // Holds the number of dimensions
    int i, j, k, dim;         // General iterators

    int print_gap = steps / 100;

    if(argc != 3)
    {
        printf("Please provide L and np, \n");
        return -1;
    }
    if(argc == 3){
        L = atoi(argv[1]);
        np  = atoi(argv[2]);
        if(L*L < 4*np){
            printf("Too many particles");
            return -1;
        }
        else{
            printf("L: %d, np: %d\n", L, np);
        }
    }

    // Chrono Time Variables
    chrono::time_point<chrono::steady_clock> begin_time, end_time;
    chrono::duration<double> difference_in_time;
    double difference_in_seconds;

    // Initialize memory
    pos = new double[nd*np];
    vel = new double[nd*np];
    acc = new double[nd*np];
    neighbor = new vector<int>[np]; // array of vectors

    // Initialize Ekin
    Ekin = 0;

    // detect_threads_setting();

    // Initialize velocities
    for (i = 0; i < np; i++) {
        for (dim = 0; dim < nd; dim++) {
            vel[dim+i*nd] = 0.0;
        }
    }

    // Initialize positions
    setPositions(np, nd, L, pos);

    // Initialize neighborhood table 
    createNeighborhood(np, nd, L, pos, neighbor);

    // Initialize potential energy
    calcAccel(np, nd, L, mass, pos, acc, Epot, neighbor);

    // Initialize acceleration to zero
    for (i = 0; i < np; i++) {
        for (dim = 0; dim < nd; dim++) {
            acc[dim+i*nd] = 0.0;
        }
    }

    begin_time = chrono::steady_clock::now();

    // Perform Verlet Velocity Calculations
#ifdef PRINTRESULT
    printf("Particle | Position (X) | Position (Y) | Velocity (X) | Velocity (Y) | Accel. (X) | Accel (Y) |\n");
#endif

    for (i = 0; i < steps; i++) {
#ifdef PRINTRESULT
        if (i == 0 || i == steps-1) {
        //if ((i+1)%print_gap == 0 || i == 0 || i == steps-1) {
            for (j = 0; j < np; j++) {
                printf(" %4d   ", j);
                printf("%14.8f       %14.8f       %14.8f       %14.8f     %14.8f       %14.8f\n",
                       pos[0+j*nd],pos[1+j*nd],vel[0+j*nd],vel[1+j*nd],acc[0+j*nd],acc[1+j*nd]);
            }
        }
#endif
        calcVerlet (np, nd, mass, dt, L, Ekin, Epot, pos, vel, acc, neighbor);
    }

    end_time = chrono::steady_clock::now();
    difference_in_time = end_time - begin_time;
    difference_in_seconds = difference_in_time.count();

    cout << np << "\t" << L << "\t"
         << scientific << setprecision(3) << dt << "\t"
         << resetiosflags(ios::scientific) << setprecision(15) << difference_in_seconds << endl;

    // Clean memory
    delete [] pos;
    delete [] vel;
    delete [] acc;
    delete [] neighbor;
    return 0;
}

void setPositions (int np, int nd, int L, double *pos) {
    /**
     * Used to randomly intialize the positions of the particles in the system.
     * Code adapted from csm_ch08.pdf.
     */

    // Variable Declarations
    double min_r = pow(2.0, 1.0/3.0); // Used to prevent overlap
    double dist_diff[nd];             // Used to hold the difference in distance
    double overlap_check;             // Used to hold the computation for overlap checking
    bool overlap;                     // Used to indicate overlap occured
    int i, j, dim;                    // General iterators

    for (i = 0; i < np; i++) { // Scans through all particles
        do {
            overlap = false;
            for (dim = 0; dim < nd; dim++) { // Scans through all dimensions of particle
                pos[dim+i*nd] = L/4 + (L*static_cast<double>(rand()))/(static_cast<double>(RAND_MAX)*2);
            }
            j = 0;
            while (j < i && !overlap) {
                overlap_check = 0;
                for (dim = 0; dim < nd; dim++) {
                    overlap_check += (pos[dim+i*nd]-pos[dim+j*nd])
                                     * (pos[dim+i*nd]-pos[dim+j*nd]);
                }
                if (overlap_check < min_r) { overlap = true; }
                j++;
            }
        } while (overlap);
    }
}

void createNeighborhood (int np, int nd, int L, double *pos, vector<int> *neighbor) {
    /**
     * Used to create the neighborhood table for the simulation.
     */

    // Variable Declarations
    int i, j, dim;               // General iterators
    double dist_diff[nd];            // Difference in distance
    bool in_range;               // indicates if a particle is a neighbor

    // iterated over np*(np-1)/2 ordered pairs, 
    // add to neighbor list if distance <= skin 

    for (i = 0; i < np; i++) {
        for(j = 0; j < np; j++) {
            if (i == j) {
                continue;
            }
            for (dim = 0; dim < nd; dim++) {
                dist_diff[dim] = pos[dim + i*nd] - pos[dim + j*nd];
                if (abs(dist_diff[dim]) > (L - R_skin)) { // periodic boundary 
                    dist_diff[dim] = dist_diff[dim] >= 0 ? dist_diff[dim] - L : dist_diff[dim] + L;
                }
                in_range = (abs(dist_diff[dim]) <= R_skin);// only particles in R_skin
                if (!in_range) break;
            }
            if (in_range) {
                neighbor[i].push_back(j);
                // neighbor[j].push_back(i);
            }
        }
    }
}



void calcAccel (int np, int nd, int L, double mass, double *pos, double *acc, double &Epot, vector<int> *neighbor) {
    /**
     * Used to calculate the forces (aka acceleration) based on the interations
     * between particles.
     */

    // Variable Declarations
    double dist_diff[nd];        // Holds the difference in distance
    double r2, r2_over, r6_over; // Holds the r^2, r2_over, r6_over for force calculations
    double fp;                   // Holds the coefficient for the force equation
    double finst[nd];            // Holds the instantanous force (in each direction)
    int i, j, dim;               // General iterators
    double R_cut_2 = R_cut * R_cut;
    // Zero out force & Epot
    double local_Epot = 0.0;

    for (i = 0; i < np; i++) {
        for (dim = 0; dim < nd; dim++) {
            acc[dim+i*nd] = 0.0;
        }
    }

    // Perform force & Epot calculation based on proximity
    for (i = 0; i < np; i++){
        while(!neighbor[i].empty()){
            j = neighbor[i].back();
            r2 = 0.0;
            for (dim = 0; dim < nd; dim++){
                dist_diff[dim] = pos[dim + i*nd] - pos[dim + j*nd];
                r2 += dist_diff[dim] * dist_diff[dim];
            }
            if (R_cut_2 > r2){ // particle within R_cut
                r2_over = 1.0/r2;
                r6_over = r2_over*r2_over*r2_over; // 1/((r2)^3)
                fp = 48.0*r6_over*(r6_over-0.5)*r2_over;
                for (dim = 0; dim < nd; dim++) {
                    finst[dim] = fp * dist_diff[dim];
                    acc[dim+i*nd] += finst[dim]/mass;
                }
                local_Epot += 4.0*(r6_over*r6_over-r6_over);
            }
            neighbor[i].pop_back();
        }
    }

    Epot = local_Epot;
    
}

void boundCond (int np, int nd, int L, double *pos) {
    /**
     * Used to reposition particles that exit the simulation area.
     */

    // Variable Declaration
    int i, dim; // General Iterators

    for (i = 0; i < np; i++) {
        for (dim = 0; dim < nd; dim++) {
            if (pos[dim+i*nd] < 0) { // Particle has exited the lower boundary.
                pos[dim+i*nd] = fmod(pos[dim+i*nd], L)+L;
            } else if (pos[dim+i*nd] > L) { // Particle has exited the upper boundary.
                pos[dim+i*nd] = fmod(pos[dim+i*nd], L);
            } else { // Particle is right on a boundary or in the simulation region.
                continue;
            }
        }
    }
}

void calcVerlet (int np, int nd, double mass, double dt, int L, double &Ekin, double &Epot,
                 double *pos, double *vel, double *acc, vector<int> *neighbor) {
    /**
     * Used to perform the Verlet velocity calculation (based off the leap frog method)
     */

    // Variable Declarations
    double dt2 = 0.5 * dt; // Holds the half time step for the velocity calculation
    int i, dim;            // General iterators

    for (i = 0; i < np; i++) {
        for (dim = 0; dim < nd; dim++) {
            vel[dim+i*nd] += acc[dim+i*nd] * dt2;
            pos[dim+i*nd] += vel[dim+i*nd] * dt;
        }
    }

    boundCond (np, nd, L, pos);
    createNeighborhood(np, nd, L, pos, neighbor);
    calcAccel (np, nd, L, mass, pos, acc, Epot, neighbor);
    double local_Ekin = 0.0;

    for (i = 0; i < np; i++) {
        for (dim = 0; dim < nd; dim++) {
            vel[dim+i*nd] += acc[dim+i*nd] * dt2;
            local_Ekin += vel[dim+i*nd] * vel[dim+i*nd];
        }
    }
    Ekin = 0.5 * local_Ekin/np;
}
