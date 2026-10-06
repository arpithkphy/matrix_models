// matrix.cpp -- the N x N complex matrix class and the few operations the program needs.
#include "model.h"      // Matrix, NCOLOR, complex

Matrix::Matrix() {                                    // default matrix: all entries zero
    for (int i = 0; i < NCOLOR; i++) {                // loop over rows
        for (int j = 0; j < NCOLOR; j++) {            // loop over columns
            entry[i][j] = 0.0;                        // set the entry to zero
        }
    }
}

complex<double> Matrix::get(int row, int column) const {   // read one entry
    return entry[row][column];                        // read one entry
}

void Matrix::set(int row, int column, const complex<double> value) {   // write one entry
    entry[row][column] = value;                       // write one entry
}

Matrix operator+(const Matrix &a, const Matrix &b) {   // sum of two matrices
    Matrix sum;                                       // result, starts as zero
    for (int i = 0; i < NCOLOR; i++)                  // rows
        for (int j = 0; j < NCOLOR; j++)              // columns
            sum.set(i, j, a.get(i, j) + b.get(i, j)); // add entry by entry
    return sum;   // return the sum
}

Matrix operator-(const Matrix &a, const Matrix &b) {   // difference of two matrices
    Matrix difference;                                // result, starts as zero
    for (int i = 0; i < NCOLOR; i++)                  // rows
        for (int j = 0; j < NCOLOR; j++)              // columns
            difference.set(i, j, a.get(i, j) - b.get(i, j));   // subtract entry by entry
    return difference;   // return the difference
}

Matrix operator*(const Matrix &a, const Matrix &b) {   // matrix product
    Matrix product;                                   // result, starts as zero
    for (int i = 0; i < NCOLOR; i++) {                // row of the result
        for (int j = 0; j < NCOLOR; j++) {            // column of the result
            complex<double> total = 0.0;              // running sum over the inner index
            for (int k = 0; k < NCOLOR; k++) {        // inner index
                total = total + a.get(i, k) * b.get(k, j);   // (a b)_ij = sum_k a_ik b_kj
            }
            product.set(i, j, total);                 // store the finished entry
        }
    }
    return product;   // return the product
}

Matrix operator*(const Matrix &a, const double s) {   // matrix times a real number
    Matrix scaled;                                    // result, starts as zero
    for (int i = 0; i < NCOLOR; i++)                  // rows
        for (int j = 0; j < NCOLOR; j++)              // columns
            scaled.set(i, j, a.get(i, j) * s);        // multiply every entry by s
    return scaled;   // return the scaled matrix
}

Matrix operator*(const double s, const Matrix &a) {   // real number times matrix
    Matrix scaled;                                    // result, starts as zero
    for (int i = 0; i < NCOLOR; i++)                  // rows
        for (int j = 0; j < NCOLOR; j++)              // columns
            scaled.set(i, j, a.get(i, j) * s);        // multiply every entry by s
    return scaled;   // return the scaled matrix
}

Matrix dagger(const Matrix &a) {   // Hermitian conjugate
    Matrix result;                                    // result, starts as zero
    for (int i = 0; i < NCOLOR; i++)                  // rows
        for (int j = 0; j < NCOLOR; j++)              // columns
            result.set(i, j, conj(a.get(j, i)));      // (a^dagger)_ij = conjugate of a_ji
    return result;   // return the conjugate transpose
}

Matrix commutator(const Matrix &a, const Matrix &b) {   // commutator [a, b]
    return a * b - b * a;                             // [a, b] = a b - b a
}

complex<double> trace_of(const Matrix &a) {   // trace of a matrix
    complex<double> total = 0.0;                      // running sum
    for (int i = 0; i < NCOLOR; i++)                  // diagonal entries
        total = total + a.get(i, i);                  // tr a = sum_i a_ii
    return total;   // return the trace
}
