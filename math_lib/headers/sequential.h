#include "general.h"

#ifndef SEQUENTIAL_H_
#define SEQUENTIAL_H_

// OPERATION
// complex
complex c_add(complex a, complex b);
complex c_sub(complex a, complex b);
complex c_mult(complex a, complex b);
complex c_scale(complex a, double b);
complex c_div(complex a, complex b);
complex c_exp(complex a);
complex c_log(complex a);
complex c_conj(complex a);
uint8_t c_equal(complex a, complex b);
double c_norm(complex a);
double c_norm2(complex a);

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
uint8_t *v_equal(vector *a, vector *b, vector *result);
double v_dot_prod(vector *a, vector *b);
double v_norm(vector *a, vector *b);
double v_norm2(vector *a, vector *b);

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

// tensor

#endif // !SEQUENTIAL_H_
