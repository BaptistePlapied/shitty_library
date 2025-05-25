#include "general.h"
#include <math.h>

#ifndef SEQUENTIAL_H_
#define SEQUENTIAL_H_

// OPERATION
// complex
static inline complex c_add(complex a, complex b) {
    return (complex){a.Re + b.Re, a.Im + b.Im};
}
static inline complex c_sub(complex a, complex b) {
    return (complex){a.Re - b.Re, a.Im - b.Im};
}
static inline complex c_scale(complex a, double b) {
    return (complex){a.Re * b, a.Im * b};
}
static inline complex c_conj(complex a) { return (complex){a.Re, -a.Im}; }
static inline double c_norm2(complex a) { return a.Re * a.Re + a.Im * a.Im; }
static inline double c_norm(complex a) { return sqrt(c_norm2(a)); }
static inline double c_arg(complex a) { return atan2(a.Re, a.Im); }
static inline complex c_mult(complex a, complex b) {
    return (complex){a.Re * b.Re - a.Im * b.Im, a.Re * b.Im + a.Im * b.Re};
}

complex c_div(complex a, complex b);
complex c_exp(complex a);
complex c_log(complex a);
uint8_t c_equal(complex a, complex b);
uint8_t c_compare(complex a, complex b); // the norm2

// vector
vector *v_add(vector *a, vector *b, vector *result);
vector *v_sub(vector *a, vector *b, vector *result);
vector *v_mult_e(vector *a, vector *b, vector *result);
vector *v_div_e(vector *a, vector *b, vector *result);
vector *v_scale(vector *a, complex alpha, vector *result);
vector *v_scale_r(vector *a, double alpha, vector *result);
vector *v_conj(vector *a, vector *result);
vector *v_Re(vector *a, vector *result);
vector *v_Im(vector *a, vector *result);
vector *v_map(vector *a, complex (*f)(complex), vector *result);
uint8_t v_equal(vector *a, vector *b, vector *result);
complex v_dot_prod(vector *a, vector *b);
double v_norm(vector *a);
double v_norm2(vector *a);

// matrix

matrix *m_add(matrix *A, matrix *B, matrix *Result);
matrix *m_sub(matrix *A, matrix *B, matrix *Result);
matrix *m_mult_e(matrix *A, matrix *B, matrix *Result);
matrix *m_div_e(matrix *A, matrix *B, matrix *Result);
matrix *m_mult(matrix *A, matrix *B, matrix *Result);
matrix *m_scale(matrix *A, complex alpha, matrix *Result);
matrix *m_scale_r(matrix *A, double alpha, matrix *Result);
matrix *m_conj(matrix *A, matrix *Result); // Element-wise conjugate
matrix *m_transpose(matrix *A, matrix *Result);
matrix *m_hermitian(matrix *A, matrix *Result); // Conjugate transpose
double m_norm(matrix *A);                       // sqrt(sum |mᵢⱼ|²)
double m_norm2(matrix *A);                      // sum |mᵢⱼ|²
matrix *m_mult(matrix *A, matrix *B, matrix *Result);

// tensor

// algo

typedef struct {
    uint64_t row_target; // Row being modified
    uint64_t row_source; // Pivot row used
    complex multiplier;  // Multiplier used to eliminate the column
} row_op;

void row_op_printf(row_op ope);

typedef struct {
    uint64_t num_ops;
    row_op *operations;
    vector *pivots;
    vector *permutation; // optional, to store swaps
} echelon_expr;

void echelon_expr_printf(echelon_expr echelon);

vector *backsub(matrix *U, vector *b, vector *result);
matrix *m_echelon(matrix *A, echelon_expr *expression, matrix *Result);

#endif // !SEQUENTIAL_H_
