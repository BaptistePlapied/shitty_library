#include "general.h"
#include <math.h>
#include <stdint.h>

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
static inline complex c_zero() { return (complex){0.0, 0.0}; }
static inline complex c_real(double a) { return (complex){a, 0.0}; }
static inline complex c_Re(complex a) { return (complex){a.Re, 0.0}; }
static inline complex c_Im(complex a) { return (complex){0.0, a.Im}; }
static inline complex c_conj(complex a) { return (complex){a.Re, -a.Im}; }
static inline complex c_abs(complex a) { return (complex){fabs(a.Re), fabs(a.Im)}; }
static inline double c_norm2(complex a) { return a.Re * a.Re + a.Im * a.Im; }
static inline double c_norm(complex a) { return sqrt(c_norm2(a)); }
static inline double c_arg(complex a) { return atan2(a.Re, a.Im); }
static inline complex c_mult(complex a, complex b) {
    return (complex){a.Re * b.Re - a.Im * b.Im, a.Re * b.Im + a.Im * b.Re};
}

complex c_div(complex a, complex b);
complex c_sqrt(complex z);
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
matrix *v_outer(vector *a, vector *b, matrix *Result);

// tensor

// algo

typedef struct {
    uint64_t row_target; // Row being modified
    uint64_t row_source; // Pivot row used
    complex multiplier;  // Multiplier used to eliminate the column
} row_op;

void row_op_printf(row_op ope);

typedef struct {
    uint64_t v_size;
    uint64_t num_ops;
    row_op *operations;
    vector *pivots;
    vector *permutation; // optional, to store swaps
} echelon_expr;

void echelon_expr_printf(echelon_expr echelon);

vector *backsub(matrix *U, vector *b, vector *result);
matrix *m_echelon(matrix *A, echelon_expr *expression, matrix *Result);
vector *v_apply_echelon(vector *a, echelon_expr *expr, vector *result);
void lu_pp(matrix *A, matrix *L, matrix *U, matrix *P, uint64_t *swap);
complex m_tr_det(matrix *A);
complex m_det(matrix *A);
void m_mgs_qr(matrix *A, matrix *Q, matrix *R);

// TRACKER

#define MAX_TRACKED_GROUPS 8
#define MAX_TRACKED_IN_GROUPS 128
#define MAX_TRACKED 1028

typedef enum { TRACK_MATRIX, TRACK_VECTOR, TRACK_TENSOR } TrackType;

typedef struct {
    void *ptr;
    TrackType type;
} TrackedObject;

typedef struct {
    TrackedObject item[MAX_TRACKED_IN_GROUPS];
    uint64_t count;
} TrackedGroup;

static TrackedObject tracked[MAX_TRACKED];
static uint64_t tracked_count = 0;
static TrackedGroup groups[MAX_TRACKED_GROUPS];

void track_group_init();
void *track_generic_group_in(void *ptr, uint64_t group_id, TrackType type);
void *track_generic_group_get(uint64_t group_id, uint64_t index);
void track_generic_group_replace(uint64_t group_id, uint64_t index, void *ptr,
                                 TrackType type);
void track_group_clear(uint64_t group_id);
uint64_t track_group_count(uint64_t group_id);
uint64_t find_empty_group_id();

void *track_generic(void *ptr, TrackType type);
void track_clear_all();

// put in a group
#define m_track_in(A, group_id)                                                          \
    (matrix *)track_generic_group_in((void *)(A), (group_id), TRACK_MATRIX)
#define v_track_in(a, group_id)                                                          \
    (vector *)track_generic_group_in((void *)(a), (group_id), TRACK_VECTOR)
#define t_track_in(T, group_id)                                                          \
    (tensor *)track_generic_group_in((void *)(T), (group_id), TRACK_TENSOR)

// get
#define m_track_get(group_id, index) (matrix *)track_generic_group_get((group_id), index)
#define v_track_get(group_id, index) (vector *)track_generic_group_get((group_id), index)
#define t_track_get(group_id, index) (tensor *)track_generic_group_get((group_id), index)

// replace
#define m_track_replace(group_id, index, A)                                              \
    track_generic_group_replace((group_id), (index), (void *)(A), TRACK_MATRIX)
#define v_track_replace(group_id, index, a)                                              \
    track_generic_group_replace((group_id), (index), (void *)(a), TRACK_VECTOR)
#define t_track_replace(group_id, index, T)                                              \
    track_generic_group_replace((group_id), (index), (void *)(T), TRACK_TENSOR)

// global

#define m_track(A) (matrix *)track_generic((void *)(A), TRACK_MATRIX)
#define v_track(a) (vector *)track_generic((void *)(a), TRACK_VECTOR)
#define t_track(T) (tensor *)track_generic((void *)(T), TRACK_TENSOR)

#endif // !SEQUENTIAL_H_
