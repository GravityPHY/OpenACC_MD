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
#include <cstring>

#define R_cut 2.0
#define R_skin 4.0

using namespace std;

// Function Prototypes:
int main(int argc, char **argv );
void setPositions (int np, int nd, int L, double *pos);
void boundCond (int np, int nd, int L, double *pos);

// vector version
void createNeighborhood_vector (int np, int nd, int L, double *pos, vector<int> *neighbor);
void calcAccel_vector (int np, int nd, int L, double mass, double *pos, double *acc, double &Epot, vector<int> *neighbor);
void calcVerlet_vector (int np, int nd, double mass, double dt, int L, double &Ekin, double &Epot,
                 double *pos, double *vel, double *acc, vector<int> *neighbor);
// array version
void createNeighborhood_array (int np, int nd, int L, double *pos, int *First, int *Edge);
void calcAccel_array (int np, int nd, int L, double mass, double *pos, double *acc, double &Epot, int *First, int *Edge);
void calcVerlet_array (int np, int nd, double mass, double dt, int L, double &Ekin, double &Epot,
                 double *pos, double *vel, double *acc, int *First, int *Edge);




int main(int argc, char **argv) {
    /**
     * Main function of the Molecular Dynamics Simulator
     */

    // Variable Declartions
    double *pos, *pos_array;              // Holds the position of the particles
    double *vel, *vel_array;              // Holds the velocity of the particles
    double *acc, *acc_array;              // Holds the accleration of the particles
    int L = 4096;             // length of the simulation box
    double Ekin, Epot, Ekin_array, Epot_array;        // Holds the kinetic & potential energy of the system.
    double mass = 1.0;        // Holds the mass of each atom
    double dt = 1E-12;        // Holds the time step
    int steps = 5000;         // Holds the maximum number of time steps to do
    int np = 1024;            // Holds the number of particles
    int nd = 2;               // Holds the number of dimensions
    int i, j, k, dim;         // General iterators

    // vector
    vector<int> *neighbor;    // Holds the neighborhood table for the simulation
    // array
    int *First;               // Holds the index of Edge
    int *Edge;                // Holds the neighborhood table

    // Chrono Time Variables
    chrono::time_point<chrono::steady_clock> begin_time, end_time;
    chrono::duration<double> difference_in_time;
    double difference_in_seconds;

    // Initialize memory
    // vector
    pos = new double[nd*np];
    vel = new double[nd*np];
    acc = new double[nd*np];
    neighbor = new vector<int>[np]; // array of vectors

    // array
    pos_array = new double[nd*np];
    vel_array = new double[nd*np];
    acc_array = new double[nd*np];
    First = new int[np+1];
    Edge = new int[np*np];

    // Initialize Ekin
    Ekin = 0;
    Ekin_array = 0;

    // Initialize velocities
    for (i = 0; i < np; i++) {
        for (dim = 0; dim < nd; dim++) {
            vel[dim+i*nd] = 0.0;
            vel_array[dim+i*nd] = 0.0;
        }
    }

    // Initialize positions
    setPositions(np, nd, L, pos);
    memcpy(pos_array, pos, sizeof(pos));

    // Initialize neighborhood table 
    createNeighborhood_vector(np, nd, L, pos, neighbor);
    createNeighborhood_array(np, nd, L, pos_array, First, Edge);

    // Initialize potential energy
    calcAccel_vector(np, nd, L, mass, pos, acc, Epot, neighbor);
    calcAccel_array(np, nd, L, mass, pos_array, acc_array, Epot_array, First, Edge);

    // Initialize acceleration to zero
    for (i = 0; i < np; i++) {
        for (dim = 0; dim < nd; dim++) {
            acc[dim+i*nd] = 0.0;
            acc_array[dim+i*nd] = 0.0;
        }
    }

    begin_time = chrono::steady_clock::now();

    // Perform Verlet Velocity Calculations

    for (i = 0; i < steps; i++) {
        calcVerlet_vector(np, nd, mass, dt, L, Ekin, Epot, pos, vel, acc, neighbor);
    }

    end_time = chrono::steady_clock::now();
    difference_in_time = end_time - begin_time;
    difference_in_seconds = difference_in_time.count();

    cout << "Vector: " << np << "\t" << L << "\t"
         << scientific << setprecision(3) << dt << "\t"
         << resetiosflags(ios::scientific) << setprecision(15) << difference_in_seconds << endl;


    begin_time = chrono::steady_clock::now();

    // Perform Verlet Velocity Calculations

    for (i = 0; i < steps; i++) {
        calcVerlet_array(np, nd, mass, dt, L, Ekin_array, Epot_array, pos_array, vel_array, acc_array, First, Edge);
    }

    end_time = chrono::steady_clock::now();
    difference_in_time = end_time - begin_time;
    difference_in_seconds = difference_in_time.count();

    cout << "Array:  " << np << "\t" << L << "\t"
         << scientific << setprecision(3) << dt << "\t"
         << resetiosflags(ios::scientific) << setprecision(15) << difference_in_seconds << endl;


    // Clean memory
    delete [] pos;
    delete [] vel;
    delete [] acc;
    delete [] neighbor;
    delete [] First;
    delete [] Edge;
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
                pos[dim+i*nd] = L*static_cast<double>(rand())/static_cast<double>(RAND_MAX);
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

void createNeighborhood_vector (int np, int nd, int L, double *pos, vector<int> *neighbor) {
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
        for(j = i+1; j < np; j++) {
            for (dim = 0; dim < nd; dim++) {
                dist_diff[dim] = pos[dim + i*nd] - pos[dim + j*nd];
                if (abs(dist_diff[dim]) > (L - R_skin - R_cut)) { // periodic boundary 
                    dist_diff[dim] = dist_diff[dim] >= 0 ? dist_diff[dim] - L : dist_diff[dim] + L;
                }
                in_range = (abs(dist_diff[dim]) <= R_skin);// ignore out range
                if (!in_range) break;
            }
            if (in_range) {
                neighbor[i].push_back(j);
                neighbor[j].push_back(i);
            }
        }
    }
}



void calcAccel_vector (int np, int nd, int L, double mass, double *pos, double *acc, double &Epot, vector<int> *neighbor) {
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
    // Zero out force & Epot
    Epot = 0.0;
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
            r2_over = 1.0/r2;
            r6_over = r2_over*r2_over*r2_over; // 1/((r2)^3)
            fp = 48.0*r6_over*(r6_over-0.5)*r2_over;
            for (dim = 0; dim < nd; dim++) {
                finst[dim] = fp * dist_diff[dim];
                acc[dim+i*nd] += finst[dim]/mass;
                acc[dim+j*nd] -= finst[dim]/mass;
            }
            Epot += 4.0*(r6_over*r6_over-r6_over);
            neighbor[i].pop_back();
        }
    }
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

void calcVerlet_vector (int np, int nd, double mass, double dt, int L, double &Ekin, double &Epot,
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
    createNeighborhood_vector(np, nd, L, pos, neighbor);
    calcAccel_vector(np, nd, L, mass, pos, acc, Epot, neighbor);
    Ekin = 0.0;
    for (i = 0; i < np; i++) {
        for (dim = 0; dim < nd; dim++) {
            vel[dim+i*nd] += acc[dim+i*nd] * dt2;
            Ekin += vel[dim+i*nd] * vel[dim+i*nd];
        }
    }
    Ekin = 0.5 * Ekin/np;
}


void createNeighborhood_array (int np, int nd, int L, double *pos, int *First, int *Edge) {
    /**
     * Used to create the neighborhood table for the simulation.
     */

    // Variable Declarations
    int i, j, dim;               // General iterators
    double dist_diff[nd];            // Difference in distance
    bool in_range;               // indicates if a particle is a neighbor
    int EdgeIndex = 0;          
    
    First[0] = 0; 

    // iterated over np*np ordered pairs, 
    // add index of particle to Edge array if distance <= skin 
    for(i=0; i<np; i++){
        First[i+1] = First[i];
        for (j=0; j< np; j++){
            for (dim = 0; dim < nd; dim++) {
                dist_diff[dim] = pos[dim + i*nd] - pos[dim + j*nd];
                if (abs(dist_diff[dim]) > (L - R_skin - R_cut)) { // periodic boundary 
                    dist_diff[dim] = dist_diff[dim] >= 0 ? dist_diff[dim] - L : dist_diff[dim] + L;
                }
                in_range = (abs(dist_diff[dim]) <= R_skin);// ignore out range
                if (!in_range) break;
            }
            if (in_range) {
                Edge[EdgeIndex] = j;
                EdgeIndex++;
                First[i+1] +=1;
            }
        }
    }
    Edge[EdgeIndex] = -1; // end indicator
}


void calcAccel_array (int np, int nd, int L, double mass, double *pos, double *acc, double &Epot, int *First, int *Edge) {
    /**
     * Used to calculate the forces (aka acceleration) based on the interations
     * between particles.
     */

    // Variable Declarations
    double dist_diff[nd];        // Holds the difference in distance
    double r2, r2_over, r6_over; // Holds the r^2, r2_over, r6_over for force calculations
    double fp;                   // Holds the coefficient for the force equation
    double finst[nd];            // Holds the instantanous force (in each direction)
    int i, j, e, dim;               // General iterators
    // Zero out force & Epot
    Epot = 0.0;
    for (i = 0; i < np; i++) {
        for (dim = 0; dim < nd; dim++) {
            acc[dim+i*nd] = 0.0;
        }
    }

    // Perform force & Epot calculation based on proximity
    for (i=0; i<np; i++){
        for (e = First[i]; e<First[i+1]; e++){ // neighbors of i
            j = Edge[e]; // content of Edge array == index of neighbor 
            r2 = 0.0;
            for (dim=0; dim<nd; dim++){
                dist_diff[dim] = pos[dim + i*nd] - pos[dim + j*nd];
                r2 += dist_diff[dim] * dist_diff[dim];
            }
            r2_over = 1.0/r2;
            r6_over = r2_over*r2_over*r2_over; // 1/((r2)^3)
            fp = 48.0*r6_over*(r6_over-0.5)*r2_over;
            for (dim = 0; dim < nd; dim++) {
                finst[dim] = fp * dist_diff[dim];
                acc[dim+i*nd] += finst[dim]/mass;
                acc[dim+j*nd] -= finst[dim]/mass;
            }
            Epot += 4.0*(r6_over*r6_over-r6_over);
        }
    }
}

void calcVerlet_array (int np, int nd, double mass, double dt, int L, double &Ekin, double &Epot,
                 double *pos, double *vel, double *acc, int *First, int *Edge) {
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
    createNeighborhood_array(np, nd, L, pos, First, Edge);
    calcAccel_array (np, nd, L, mass, pos, acc, Epot, First, Edge);
    Ekin = 0.0;
    for (i = 0; i < np; i++) {
        for (dim = 0; dim < nd; dim++) {
            vel[dim+i*nd] += acc[dim+i*nd] * dt2;
            Ekin += vel[dim+i*nd] * vel[dim+i*nd];
        }
    }
    Ekin = 0.5 * Ekin/np;
}
