// evolve.cpp -- one Langevin update of all variables (the heart of the program), with an optional adaptive step.
//
// The update is the Euler-Maruyama discretisation of the Langevin equation
//     X     <- X     - eps * dS/dX     + sqrt(eps) * eta ,   <eta_ij eta_kl> = 2 delta_il delta_jk   (matrices, then trace removed)
//     theta <- theta - eps * dS/dtheta + sqrt(eps) * xi  ,   <xi xi> = 2                              (angles, then sum removed)
// All new values are computed from the OLD values and copied in at the end (synchronous update).
//
// ADAPTIVE STEP (switch off with adaptive_factor = 0 in params.txt): before an update the largest force
//     Kmax = max over all matrices (Frobenius norm of dS/dX) and angles (|dS/dtheta|)
// is compared with  adaptive_factor x (running average of Kmax).  If Kmax is larger, this one step is shortened to
//     eps_n = eps * adaptive_factor * average / Kmax ,
// so that a rare very large force (for example two angles almost coinciding, where the Faddeev-Popov force ~ 1/distance is huge)
// cannot throw the configuration far away. The noise uses the same shortened step. In all other steps eps_n = eps exactly, so the
// systematic Euler error of the base step is NOT removed by this device. Langevin time advances by eps_n.
#include "model.h"      // Matrix, matrices, angles, holonomy, base_step, adaptive_factor, ...

const double AVERAGE_WEIGHT = 1.0e-3;            // weight of the newest Kmax in the running average (average over about 1000 steps)

static double running_mean = 0.0;                // running average of Kmax; 0 before the first step
static long steps_taken = 0;                     // number of steps so far
static long steps_shortened = 0;                 // number of steps with eps_n < eps
static double smallest_ratio = 1.0;              // smallest eps_n / eps that occurred

// Frobenius norm of a matrix: sqrt( sum_ij |m_ij|^2 ).
static double frobenius_norm(const Matrix &m) {   // helper: Frobenius norm of a matrix
    double sum = 0.0;                            // running sum of |m_ij|^2
    for (int i = 0; i < NCOLOR; i++) {   // rows
        for (int j = 0; j < NCOLOR; j++) {   // columns
            sum += norm(m.get(i, j));            // norm(z) = |z|^2
        }
    }
    return sqrt(sum);   // square root of the sum of squares
}

void print_adaptive_report() {   // prints the summary written at the end of run.log
    cout << "-- ADAPTIVE STEP : factor " << adaptive_factor << ", steps shortened " << steps_shortened << " of " << steps_taken   // how many steps were shortened of all steps
         << ", smallest eps_n/eps " << smallest_ratio << ", final Langevin time " << langevin_time << "\n";   // smallest step ratio and the final Langevin time
}

void langevin_step() {   // one Langevin update of everything (matrices and angles)
    // ---- 1. all forces from the old configuration, and the largest one ----
    Matrix force_matrix[D][L];                   // dS/dX^dim_site for all dim, site
    double force_angle[NCOLOR];                  // dS/dtheta_color for all colours
    double largest_force = 0.0;                  // Kmax
    for (int site = 0; site < L; site++) {   // every time site
        for (int dim = 0; dim < D; dim++) {   // every matrix
            force_matrix[dim][site] = force_on_matrix(dim, site);                 // force on one matrix
            largest_force = max(largest_force, frobenius_norm(force_matrix[dim][site]));   // keep the largest force seen so far
        }
    }
    for (int color = 0; color < NCOLOR; color++) {   // every colour
        force_angle[color] = force_on_angle(color);                              // force on one angle
        largest_force = max(largest_force, fabs(force_angle[color]));   // keep the largest force seen so far
    }

    // ---- 2. the step of this update ----
    double step = base_step;                     // normally the base step
    if (adaptive_factor > 0.0) {                 // adaptive step switched on
        if (running_mean == 0.0) {   // first step ever: start the running average with this force
            running_mean = largest_force;        // first step: start the running average
        }
        double allowed = adaptive_factor * running_mean;      // largest force that is still accepted with the full step
        if (largest_force > allowed) {           // a force spike: shorten this step
            step = base_step * allowed / largest_force;   // shortened step: base step times (allowed force / actual largest force)
            steps_shortened++;   // count this shortened step
            smallest_ratio = min(smallest_ratio, step / base_step);   // remember the smallest ratio eps_n / eps
        }
        running_mean += AVERAGE_WEIGHT * (largest_force - running_mean);   // update the running average
    }
    steps_taken++;   // count the step
    langevin_time += step;                       // Langevin time advances by the step actually taken

    // ---- 3. new matrices: X - step * force + sqrt(step) * noise, then remove the trace ----
    Matrix new_matrix[D][L];                     // the updated matrices
    for (int site = 0; site < L; site++) {   // every time site
        for (int dim = 0; dim < D; dim++) {   // every matrix
            Matrix noise = hermitian_noise();    // random Hermitian matrix
            new_matrix[dim][site] = matrices[dim][site] - step * force_matrix[dim][site] + noise * sqrt(step);   // X - step * dS/dX + sqrt(step) * noise (Euler-Maruyama)
            complex<double> trace = trace_of(new_matrix[dim][site]);     // keep the matrices traceless (SU(N))
            for (int c = 0; c < NCOLOR; c++) {   // subtract trace/N from the diagonal ...
                new_matrix[dim][site].set(c, c, new_matrix[dim][site].get(c, c) - (1.0 / NCOLOR) * trace);   // ... so that the new matrix is traceless
            }
        }
    }

    // ---- 4. new angles: theta - step * force + sqrt(step) * noise, then make them sum to zero ----
    double new_angle[NCOLOR];   // the updated angles
    for (int color = 0; color < NCOLOR; color++) {   // every colour
        double noise = gaussian_random(0.0, sqrt(2.0));            // real Gaussian noise with variance 2
        new_angle[color] = angles[color] - step * force_angle[color] + noise * sqrt(step);   // theta - step * dS/dtheta + sqrt(step) * noise
    }
    double angle_sum = 0.0;                      // projection to sum(theta) = 0, done once after ALL angles are updated
    for (int color = 0; color < NCOLOR; color++) {   // every colour
        angle_sum += new_angle[color];   // add up the new angles
    }
    for (int color = 0; color < NCOLOR; color++) {   // every colour
        new_angle[color] = new_angle[color] - (1.0 / NCOLOR) * angle_sum;   // subtract the mean: now the angles sum to zero
    }

    // ---- 5. copy everything in and rebuild the holonomy D = diag(exp(i theta)) ----
    for (int site = 0; site < L; site++) {   // every time site
        for (int dim = 0; dim < D; dim++) {   // every matrix
            matrices[dim][site] = new_matrix[dim][site];   // copy the new matrix in
        }
    }
    for (int color = 0; color < NCOLOR; color++) {   // every colour
        angles[color] = new_angle[color];   // copy the new angle in
        holonomy.set(color, color, exp(IOTA_COMPLEX * angles[color]));   // rebuild the holonomy entry exp(i theta)
    }
}
