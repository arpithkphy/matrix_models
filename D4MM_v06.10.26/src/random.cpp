// random.cpp -- random numbers: uniform, Gaussian, the starting matrices and the Langevin noise matrices.
// All of them draw from the C library generator rand(); main.cpp seeds it once with srand(seed).
#include "model.h"      // Matrix, NCOLOR, gaussian_random, ...

double uniform_random() {   // uniform random number in [0, 1]
    return (rand() / (double) RAND_MAX);              // uniform number in [0, 1]
}

// Marsaglia polar method: two uniform numbers in the unit disc give two independent N(0,1) numbers.
// One is returned at once, the other is kept (static) for the next call.
double gaussian_random(const double mean, const double deviation) {   // Gaussian random number with the given mean and standard deviation
    static double kept = 0.0;                         // the second number of the last pair
    static int have_kept = 0;                         // 1 if `kept` has not been used yet
    if (!have_kept) {                                 // no stored number: make a new pair
        double x, y, r;                               // point in the plane and its squared distance from the origin
        do {   // repeat until the point is accepted
            x = 2.0 * rand() / RAND_MAX - 1;          // uniform in [-1, 1]
            y = 2.0 * rand() / RAND_MAX - 1;          // uniform in [-1, 1]
            r = x * x + y * y;                        // squared radius
        } while (r == 0.0 || r > 1.0);                // repeat until the point is inside the unit disc (and not the origin)
        double factor = sqrt(-2.0 * log(r) / r);      // polar-method scale factor
        kept = y * factor;                            // second Gaussian number, kept for the next call
        have_kept = 1;                                // remember that one is stored
        return x * factor * deviation + mean;         // first Gaussian number with the requested mean and deviation
    } else {
        have_kept = 0;                                // use up the stored number
        return kept * deviation + mean;   // return the stored number, scaled
    }
}

// Starting configuration: real symmetric matrix with Gaussian entries (variance 2), trace removed.
Matrix hot_start_matrix() {   // random starting matrix: real symmetric, traceless
    Matrix m;                                         // result, starts as zero
    for (int i = 0; i < NCOLOR; i++) {                // upper triangle including the diagonal
        for (int j = i; j < NCOLOR; j++) {   // upper triangle: j from i to N-1
            if (i != j) {   // off-diagonal entries
                m.set(i, j, gaussian_random(0.0, sqrt(2.0)));   // off-diagonal entry
                m.set(j, i, m.get(i, j));                       // symmetric partner
            } else {
                m.set(i, i, gaussian_random(0.0, sqrt(2.0)));   // diagonal entry
            }
        }
    }
    complex<double> trace = trace_of(m);              // the trace to remove
    for (int c = 0; c < NCOLOR; c++) {   // every diagonal entry
        m.set(c, c, m.get(c, c) - (1.0 / NCOLOR) * trace);      // subtract trace/N from every diagonal entry
    }
    return m;   // return the finished matrix
}

// Langevin noise for a Hermitian matrix field: <eta_ij eta_kl> = 2 delta_il delta_jk.
//   diagonal entries   : real, variance 2
//   off-diagonal entries: complex, real and imaginary part each with variance 1, eta_ji = conj(eta_ij)
Matrix hermitian_noise() {   // random Hermitian matrix with <eta_ij eta_kl> = 2 delta_il delta_jk
    Matrix noise;                                     // result, starts as zero
    for (int i = 0; i < NCOLOR; i++) {   // upper triangle: i is the row ...
        for (int j = i; j < NCOLOR; j++) {   // ... j >= i is the column
            if (i != j) {   // off-diagonal entries
                double real_part = gaussian_random(0.0, 1.0);        // variance 1
                double imaginary_part = gaussian_random(0.0, 1.0);   // variance 1
                noise.set(i, j, complex<double>(real_part, imaginary_part));            // upper entry
                noise.set(j, i, conj(noise.get(i, j)));                                 // lower entry: Hermitian partner
            } else {
                noise.set(i, j, gaussian_random(0.0, sqrt(2.0)));                       // diagonal: real, variance 2
            }
        }
    }
    return noise;   // the trace is removed by langevin_step after the update
}
