// model.h -- everything the program shares: the physics constants, the field variables, the matrix class,
// and the declarations of all functions. Every .cpp file includes only this header.
//
// THE MODEL (bosonic D=4 matrix quantum mechanics, N x N Hermitian traceless matrices X^i, i = 1..3, on a periodic time lattice):
//   S = (N/lambda) sum_t a tr[ (1/(2 a^2)) (U_t X^i_{t+1} U_t^dag - X^i_t)^2  -  (1/4) sum_{i,j} [X^i_t, X^j_t]^2 ]
//       - sum_{j!=k} ln| sin((theta_j - theta_k)/2) |
//   with a = beta/L the lattice spacing, U_t = 1 except on the last link U_{L-1} = D = diag(exp(i theta_1), ..., exp(i theta_N)),
//   the holonomy of the gauge field, and sum_j theta_j = 0. The last term is the Faddeev-Popov (Vandermonde) term of the angles.
//   Observables: Polyakov loop |tr D|/N, extent of space R2, energy per N^2 (see measure.cpp).
#ifndef MODEL_H
#define MODEL_H

#include <cmath>      // sqrt, sin, cos, exp, log
#include <complex>    // complex<double>
#include <cstdlib>    // rand, srand, exit
#include <filesystem> // create_directories (the output folder)
#include <fstream>    // ofstream, ifstream
#include <iostream>   // cout

using namespace std;   // so that complex, cout, ... need no "std::" in front

// ---------- physics constants: edit here and recompile (run.py does this for every run) ----------
const int D = 3;                        // number of matrices X^i (3 for the D=4 model: 3 scalars + 1 time direction)
const int L = 20;                       // number of lattice sites in the time direction
const int NCOLOR = 4;                   // matrix size N (gauge group SU(N))
const double TEMP = 2.5;                // temperature T = 1/beta
const double SPACING = 1.0 / (TEMP * L);   // lattice spacing a = beta / L
const double LAMBDA = 1.0;              // 't Hooft coupling; the action carries the weight N / lambda

const complex<double> IOTA_COMPLEX(0.0, 1.0);   // the imaginary unit i
const double PI = 3.1415926535;                 // pi (only used for the random start of the angles)

// ---------- N x N complex matrix ----------
class Matrix {   // the N x N complex matrix
private:   // the entries are hidden from the rest of the program ...
    complex<double> entry[NCOLOR][NCOLOR];      // the NCOLOR x NCOLOR entries
public:   // ... and reached only through get / set
    Matrix();                                   // the zero matrix
    complex<double> get(int row, int column) const;                 // entry (row, column)
    void set(int row, int column, const complex<double> value);     // set entry (row, column)
};

Matrix operator+(const Matrix &a, const Matrix &b);      // a + b
Matrix operator-(const Matrix &a, const Matrix &b);      // a - b
Matrix operator*(const Matrix &a, const Matrix &b);      // matrix product
Matrix operator*(const Matrix &a, const double s);       // matrix times a real number
Matrix operator*(const double s, const Matrix &a);       // real number times matrix

Matrix dagger(const Matrix &a);                          // Hermitian conjugate (conjugate transpose)
Matrix commutator(const Matrix &a, const Matrix &b);     // [a, b] = a b - b a
complex<double> trace_of(const Matrix &a);               // tr a

// ---------- random numbers (random.cpp) ----------
double uniform_random();                                       // uniform in [0, 1]
double gaussian_random(const double mean, const double deviation);   // Gaussian with the given mean and standard deviation
Matrix hot_start_matrix();                                     // random real symmetric traceless matrix (starting configuration)
Matrix hermitian_noise();                                      // Langevin noise: Hermitian, <eta_ij eta_kl> = 2 delta_il delta_jk, traceless is imposed later

// ---------- the field variables (defined in main.cpp, used everywhere) ----------
extern Matrix matrices[D][L];            // X^i_t : matrices[i][t]
extern double angles[NCOLOR];            // holonomy angles theta_1 .. theta_N (real, sum to zero)
extern Matrix holonomy;                  // D = diag(exp(i theta)), rebuilt after every step

// ---------- run parameters (read from params.txt by read_parameters) ----------
extern double base_step;                 // Langevin step epsilon
extern int thermalization_steps;         // steps before the first measurement
extern int production_steps;             // steps after thermalisation
extern int measure_every;                // measure every this many production steps
extern double adaptive_factor;           // adaptive step: 0 = off, otherwise shorten a step whose largest force exceeds factor x its running mean
extern int seed;                         // seed of the random number generator
extern double langevin_time;             // Langevin time since the start of production (sum of the steps taken)

// ---------- the steps of the program ----------
void read_parameters();                  // params.cpp : read params.txt and print the run set-up
Matrix force_on_matrix(const int &dim, const int &site);   // forces.cpp : dS/dX^dim_site
double force_on_angle(const int &color);                   // forces.cpp : dS/dtheta_color
void langevin_step();                    // evolve.cpp : one Langevin update of all variables
void print_adaptive_report();            // evolve.cpp : how often the adaptive step shortened a step
void measure();                          // measure.cpp : write the observables of the current configuration

#endif
