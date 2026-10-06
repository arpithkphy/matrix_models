// main.cpp -- run the Langevin simulation: set up, thermalise, then generate configurations and measure them.
//
// How the model grows out of simpler ones (each step adds one ingredient):
//   1. matrix harmonic oscillator       one Hermitian matrix X(t),     S = (N/lambda) int dt tr[ (1/2) (dX/dt)^2 + (1/2) m^2 X^2 ]
//   2. + commutator potential           several matrices X^i(t),       add  -(1/4) sum_{i,j} tr [X^i,X^j]^2
//   3. + gauge field (holonomy D)       the time derivative becomes covariant: in the static gauge only the holonomy
//                                       D = diag(exp(i theta)) on the last link and the Faddeev-Popov term remain.
//      With 3 matrices this is THIS program (the bosonic D=4 model). Adding a mass and a Myers cubic term gives bosonic mini-BMN.
#include "model.h"      // everything

// ---- definitions of the variables declared extern in model.h ----
Matrix matrices[D][L];          // X^i_t
double angles[NCOLOR];          // holonomy angles
Matrix holonomy;                // D = diag(exp(i theta))
double base_step;               // Langevin step
int thermalization_steps;       // steps before measuring
int production_steps;           // steps with measurements
int measure_every;              // measure every this many production steps
double adaptive_factor;         // adaptive step factor (0 = off)
int seed;                       // random seed
double langevin_time;           // Langevin time of the production run

int main() {   // program entry point
    read_parameters();                                  // read params.txt and print the set-up
    filesystem::create_directories("meas");             // output folder (measure.cpp writes meas/*.txt)
    srand(seed);                                        // seed the random number generator

    // ---- starting configuration: random ("hot") matrices and random angles ----
    for (int site = 0; site < L; site++) {              // every time site
        for (int dim = 0; dim < D; dim++) {             // every matrix
            matrices[dim][site] = hot_start_matrix();   // random real symmetric traceless matrix
        }
    }
    double angle_sum = 0.0;                             // the angles must sum to zero (SU(N))
    for (int color = 0; color < NCOLOR; color++) {   // for every colour
        angles[color] = PI * (uniform_random() - 0.5);  // random angle in (-pi/2, pi/2)
        angle_sum += angles[color];   // add up the angles
    }
    for (int color = 0; color < NCOLOR; color++) {   // for every colour
        angles[color] = angles[color] - (1.0 / NCOLOR) * angle_sum;      // shift so that the sum is zero
        holonomy.set(color, color, exp(IOTA_COMPLEX * angles[color]));   // D = diag(exp(i theta))
    }

    // ---- thermalisation: evolve without measuring ----
    for (int step = 1; step <= thermalization_steps; step++) {   // thermalisation loop
        langevin_step();                                // one Langevin update
    }

    // ---- production: evolve and measure ----
    langevin_time = 0.0;                                // Langevin time in the output counts from here
    for (int step = 1; step <= production_steps; step++) {   // production loop
        langevin_step();                                // one Langevin update
        if (step % measure_every == 0) {                // every measure_every steps ...
            measure();                                  // ... write the observables of this configuration
        }
    }

    print_adaptive_report();                            // summary of the adaptive step
    return 0;   // exit code 0 = success
}
