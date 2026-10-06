// forces.cpp -- the derivatives of the action S (see model.h) with respect to the variables.
// The Langevin update uses minus these "forces": variable <- variable - step * dS/dvariable + noise.
#include "model.h"      // Matrix, matrices, angles, holonomy, NCOLOR, L, D, SPACING, LAMBDA

// dS/dX^dim_site : the derivative of S with respect to the matrix X^dim at time site `site`.
//   S contains  (N/lambda) a tr[ (1/(2 a^2)) (X_{t+1} - X_t)^2 ]  -  (N/lambda) a (1/4) sum_{i,j} tr [X^i,X^j]^2.
//   The first term gives (N/(lambda a)) (2 X_t - X_{t+1} - X_{t-1}); the neighbour that crosses the last link is rotated by the holonomy.
//   The second term gives -(N a/lambda) sum_j [X^j, [X^dim, X^j]]   (the j = dim term vanishes).
Matrix force_on_matrix(const int &dim, const int &site) {   // force on matrix X^dim at time site `site`
    Matrix neighbours;                                // will hold X_{t+1} - 2 X_t + X_{t-1}  (with the holonomy on the last link)
    Matrix force;                                     // the result

    if (site == 0) {                                  // first site: the previous neighbour is the last site, across the last link
        neighbours = matrices[dim][(site + 1) % L] - 2.0 * matrices[dim][site]   // site 0: X_1 - 2 X_0 + (previous neighbour, rotated by D) ...
                     + dagger(holonomy) * matrices[dim][(L - 1) % L] * holonomy;      // D^dagger X_{L-1} D
    } else if (site == (L - 1) % L) {                 // last site: the next neighbour is site 0, across the last link
        neighbours = matrices[dim][(site - 1) % L] - 2.0 * matrices[dim][site]   // last site: X_{L-2} - 2 X_{L-1} + (next neighbour X_0 rotated by D) ...
                     + (holonomy) * matrices[dim][0] * dagger(holonomy);              // D X_0 D^dagger
    } else {                                          // interior site: plain neighbours
        neighbours = matrices[dim][(site + 1) % L] - 2.0 * matrices[dim][site]   // interior site: X_{t+1} - 2 X_t + X_{t-1}
                     + matrices[dim][(site - 1) % L];   // (second line of the same expression)
    }

    force = -(1.0 * NCOLOR / (SPACING * LAMBDA)) * neighbours;     // kinetic part: (N/(lambda a)) (2 X - X_+ - X_-)

    for (int i = 0; i < D; i++) {                     // commutator potential: sum over the other matrices j (called i here)
        force = force - 1.0 * (1.0 * NCOLOR * SPACING / LAMBDA)   // subtract (N a/lambda) [X^j,[X^dim,X^j]] ...
                        * commutator(matrices[i][site], commutator(matrices[dim][site], matrices[i][site]));   // -(N a/lambda) [X^j,[X^dim,X^j]]
    }

    return force;   // dS/dX^dim_site
}

// dS/dtheta_color : the derivative of S with respect to the holonomy angle theta_color.
//   Kinetic part (only the last link depends on the angles):
//        -(2N/(lambda a)) sum_{dim,m} Re[ i (X^dim_{L-1})_{m,color} (X^dim_0)_{color,m} exp(i (theta_color - theta_m)) ]
//   Faddeev-Popov part:  -sum_{m != color} cot((theta_color - theta_m)/2)   (it diverges when two angles coincide).
//   Below, `sum` is accumulated with the opposite sign (+Re[...] and +cot) and the sign is flipped in the last line.
double force_on_angle(const int &color) {   // force on the angle theta_color
    double sum = 0.0;                                 // accumulates minus the derivative, as in the update formula

    for (int d = 0; d < D; d++) {                     // every matrix X^d
        for (int m = 0; m < NCOLOR; m++) {            // every colour index m
            sum = sum + (IOTA_COMPLEX * matrices[d][(L - 1) % L].get(m, color)   // i * X_{L-1}(m,c) * X_0(c,m) * exp(i (theta_c - theta_m)), real part taken ...
                                      * matrices[d][0].get(color, m)   // (second line of the same expression)
                                      * exp(IOTA_COMPLEX * (angles[color] - angles[m]))).real();   // Re[ i X_{L-1}(m,c) X_0(c,m) e^{i(theta_c - theta_m)} ]
        }
    }
    sum = sum * (2.0 * NCOLOR / (SPACING * LAMBDA));  // factor 2N/(lambda a)

    for (int m = 0; m < NCOLOR; m++) {                // Faddeev-Popov term
        if (m != color) {   // skip m = color (the Faddeev-Popov term needs two different angles)
            sum = sum + (cos(0.5 * (angles[color] - angles[m])))   // + cos / sin = cot of half the angle difference (continued on the next line)
                      / (sin(0.5 * (angles[color] - angles[m])));   // + cot((theta_c - theta_m)/2)
        }
    }

    return -1.0 * sum;                                // sum was minus dS/dtheta, so return dS/dtheta
}
