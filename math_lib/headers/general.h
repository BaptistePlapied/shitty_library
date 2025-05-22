#include <stdbool.h>
#include <stdint.h>

#ifndef GENERAL_H_
#define GENERAL_H_

typedef struct {
    double Re, Im;
} complex;

typedef struct {
    uint64_t m;
    complex *data;
} vector;
typedef struct {
    uint64_t m, n;
    complex *data;
} matrix;

typedef struct {
    uint64_t ndim;   // number of dimensions
    uint64_t *shape; // size of each dimension
    complex *data;   // flattened data array
} tensor;

// SPINAL FUNCTION
// complex

// vector
vector *v_init(uint64_t m);
vector *v_init_0(uint64_t m);
vector *v_init_l(uint64_t m, complex *c);
vector *v_resize(vector *a, uint64_t m);
vector *v_copy(vector *a);
void v_free(vector *a);
void v_printf(vector *a);

// matrix
matrix *m_init(uint64_t m, uint64_t n);
matrix *m_init_0(uint64_t m, uint64_t n);
matrix *m_init_l(uint64_t m, uint64_t n, complex *c);
matrix *m_reshape(matrix *A, uint64_t m, uint64_t n);
matrix *m_copy(matrix *A);
void m_free(matrix *A);
void m_printf(matrix *A);

// tensor
tensor *t_init(uint64_t *shape, uint64_t ndim);
tensor *t_init_0(uint64_t *shape, uint64_t ndim);
tensor *t_init_l(uint64_t *shape, uint64_t ndim, complex *c);
tensor *t_reshape(tensor *T, uint64_t *shape, uint64_t ndim);
tensor *t_copy(tensor *T);
void t_free(tensor *T);
void t_printf(tensor *T);

// SWAP STRUCT
// complex
vector *complex_to_vector(complex *data, uint64_t m);
matrix *complex_to_matrix(complex *data, uint64_t m, uint64_t n);
tensor *complex_to_tensor(complex *data, uint64_t *shape, uint64_t ndim);

// vector
void vector_to_complex(vector *a, complex *c, uint64_t *m);
matrix *vector_to_matrix(vector *a, bool as_col);
tensor *vector_to_tensor(vector *a);

// matrix
void matrix_to_complex(matrix *A, complex *c, uint64_t *m);
vector *matrix_to_vector(matrix *A);
tensor *matrix_to_tensor(matrix *A);

// tensor
void tensor_to_complex(tensor *T, complex *c, uint64_t *m);
matrix *tensor_to_matrix(tensor *T);
vector *tensor_to_vector(tensor *T);

#endif // GENERAL_H_
