#include "../headers/general.h"
#include "../headers/sequential.h"
#include <math.h>
#include <stdint.h>
#include <stdio.h>

/*
 *
 * SWAP STRUCT
 *
 */

vector *complex_to_vector(complex *data, uint64_t m);

matrix *complex_to_matrix(complex *data, uint64_t m, uint64_t n) {
    return m_init_l(m, n, data);
}
tensor *complex_to_tensor(complex *data, uint64_t *shape, uint64_t ndim);

/*
 *
 * OPERATION FUNCTION
 *
 */

inline complex c_add(complex a, complex b) { return (complex){a.Re + b.Re, a.Im + b.Im}; }

inline complex c_sub(complex a, complex b) { return (complex){a.Re - b.Re, a.Im - b.Im}; }

inline complex c_mult(complex a, complex b) {
    return (complex){a.Re * b.Re - a.Im * b.Im, a.Re * b.Im + a.Im * b.Re};
}

inline complex c_scale(complex a, double b) { return (complex){a.Re * b, a.Im * b}; }

complex c_div(complex a, complex b) {
    double denominator = b.Re * b.Re + b.Im * b.Im;
    if (denominator == 0.0) {
        fprintf(stderr, "ERROR: Division by zero in c_div\n");
        return (complex){0.0, 0.0}; // Or use NANs or handle differently
    }
    return (complex){(a.Re * b.Re + a.Im * b.Im) / denominator,
                     (a.Im * b.Re - a.Re * b.Im) / denominator};
}
inline complex c_conj(complex a) { return (complex){a.Re, -a.Im}; }

inline double c_norm2(complex a) { return a.Re * a.Re + a.Im * a.Im; }
inline double c_norm(complex a) { return sqrt(c_norm2(a)); }
inline double c_arg(complex a) { return atan2(a.Re, a.Im); }

complex c_exp(complex a) {
    double exp_re = exp(a.Re);
    return (complex){exp_re * cos(a.Im), exp_re * sin(a.Im)};
}
complex c_log(complex a) {
    double modulus = c_norm(a);
    double angle = c_arg(a);
    return (complex){log(modulus), angle};
}
uint8_t c_equal(complex a, complex b) { return (a.Re == b.Re && a.Im == b.Im) ? 1 : 0; }
// Return: 0 = a < b, 1 = a == b, 2 = a > b
uint8_t c_compare(complex a, complex b) {
    double mag_a = c_norm2(a);
    double mag_b = c_norm2(b);
    if (fabs(mag_a - mag_b) < 1e-12)
        return 1; // Equal (with tolerance)
    return (mag_a > mag_b) ? 2 : 0;
}
