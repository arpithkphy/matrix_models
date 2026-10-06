// measure.cpp -- the observables of the current configuration, appended to three text files in meas/.
// Each call writes one line per file:  "<Langevin time> <value>"  (the extent file also has the three per-matrix extents).
//   Polyakov loop  P  = |tr D| / N                                     0 <= P <= 1;  D = holonomy = diag(exp(i theta))
//   extent         R2 = (lambda/(N L)) sum_{t,i} tr (X^i_t)^2          summed over the three matrices
//   energy         E  = -3 lambda/(4 N L) sum_{t,i,j} tr [X^i_t, X^j_t]^2    (energy per N^2; virial estimator, exact at lambda = 1)
#include "model.h"      // matrices, holonomy, langevin_time, commutator, trace_of, ...

void measure() {   // write the observables of the current configuration
    static ofstream extent_file("meas/extent.txt");       // opened on the first call, kept open afterwards
    static ofstream polyakov_file("meas/polyakov.txt");   // polyakov file, opened once
    static ofstream energy_file("meas/energy.txt");   // energy file, opened once
    extent_file.precision(10);                             // write 10 significant digits
    polyakov_file.precision(10);   // 10 significant digits
    energy_file.precision(10);   // 10 significant digits

    // ---- extent of space: one value per matrix, and their sum ----
    complex<double> extent[D];                             // extent[i] = (lambda/(N L)) sum_t tr (X^i_t)^2
    complex<double> extent_total = 0.0;                    // sum over the three matrices
    for (int dim = 0; dim < D; dim++) {   // for every matrix
        extent[dim] = 0.0;   // start the sum at zero
        for (int site = 0; site < L; site++) {   // sum over the time sites
            extent[dim] += (LAMBDA / (1.0 * NCOLOR * L)) * trace_of(matrices[dim][site] * matrices[dim][site]);   // (lambda/(N L)) tr (X^i_t)^2
        }
        extent_total += extent[dim];   // add to the total over the three matrices
    }
    extent_file << langevin_time << "\t" << extent_total.real() << "\t";   // time, then the total
    for (int dim = 0; dim < D; dim++) {   // for every matrix
        extent_file << extent[dim].real() << "\t";                          // then the three single extents
    }
    extent_file << endl;   // end of the line

    // ---- Polyakov loop ----
    double polyakov = (1.0 / (1.0 * NCOLOR)) * abs(trace_of(holonomy));    // |tr D| / N
    polyakov_file << langevin_time << "\t" << polyakov << endl;   // time and Polyakov loop on one line

    // ---- energy per N^2 ----
    complex<double> energy = 0.0;   // running sum of tr [X^c,X^d]^2
    for (int site = 0; site < L; site++) {   // every time site
        for (int c = 0; c < D; c++) {   // every matrix c
            for (int d = 0; d < D; d++) {   // every matrix d
                Matrix comm = commutator(matrices[c][site], matrices[d][site]);    // [X^c, X^d]
                energy += trace_of(comm * comm);                                    // tr [X^c, X^d]^2
            }
        }
    }
    energy = energy * (-3.0 * LAMBDA / (4.0 * NCOLOR * L));   // multiply by -3 lambda / (4 N L): energy per N^2
    energy_file << langevin_time << "\t" << energy.real() << endl;   // time and energy on one line
}
