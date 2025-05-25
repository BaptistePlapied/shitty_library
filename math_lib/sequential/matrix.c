#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../headers/general.h"
#include "../headers/sequential.h"

/*
 *
 * SPINAL FUNCTION
 *
 */

matrix *m_init(uint64_t m, uint64_t n) {
    matrix *A = (matrix *)malloc(sizeof(matrix));
    if (!A) {
        fprintf(stderr, "ERROR: Failed to allocate matrix struct\n");
        return NULL;
    }
    A->m = m;
    A->n = n;
    A->data = (complex *)malloc(sizeof(complex) * (m * n));
    if (!A->data) {
        fprintf(stderr, "ERROR: Failed to allocate matrix data\n");
        free(A);
        return NULL;
    }
    return A;
}

matrix *m_init_0(uint64_t m, uint64_t n) {
    matrix *A = (matrix *)malloc(sizeof(matrix));
    if (!A) {
        fprintf(stderr, "ERROR: Failed to allocate matrix struct\n");
        return NULL;
    }
    A->m = m;
    A->n = n;
    A->data = (complex *)calloc(m * n, sizeof(complex));
    if (!A->data) {
        fprintf(stderr, "ERROR: Failed to allocate matrix data\n");
        free(A);
        return NULL;
    }
    return A;
}
matrix *m_init_l(uint64_t m, uint64_t n, complex *c) {
    if (!c) {
        fprintf(stderr, "ERROR: Input complex array is NULL\n");
        return NULL;
    }
    matrix *A = (matrix *)malloc(sizeof(matrix));
    if (!A) {
        fprintf(stderr, "ERROR: Failed to allocate matrix struct\n");
        return NULL;
    }
    A->m = m;
    A->n = n;
    A->data = (complex *)malloc(sizeof(complex) * (m * n));
    if (!A->data) {
        fprintf(stderr, "ERROR: Failed to allocate matrix data\n"); // Not "complex array
                                                                    // does not exist"
        free(A);
        return NULL;
    }
    memcpy(A->data, c, sizeof(complex) * m * n);
    return A;
}

void m_free(matrix *A) {
    if (A) {
        free(A->data);
        free(A);
    }
}

matrix *m_reshape(matrix *A, uint64_t m, uint64_t n) {
    if (!A) {
        fprintf(stderr, "ERROR: Input matrix array is NULL\n");
        return NULL;
    }
    if (A->m * A->n == m * n) {
        matrix *R = m_init_l(m, n, A->data);
        return R;
    } else {
        fprintf(stderr, "ERROR: new dim for matrix uncompatible");
        return NULL;
    }
}

matrix *m_resize(matrix *A, uint64_t m, uint64_t n) {
    if (!A) {
        fprintf(stderr, "ERROR: Input matrix array is NULL\n");
        return NULL;
    }
    matrix *copy_A = m_init(m, n);
    if (!copy_A) {
        fprintf(stderr, "ERROR: Failed to allocate resized matrix\n");
        return NULL;
    }
    for (uint64_t i = 0; i < m; i++) {
        for (uint64_t j = 0; j < n; j++) {
            if (i < A->m && j < A->n) {
                copy_A->data[i * n + j] = A->data[i * A->n + j];
            } else {
                copy_A->data[i * n + j] = (complex){0.0, 0.0};
            }
        }
    }
    return copy_A;
}

matrix *m_copy(matrix *A) {
    if (!A) {
        fprintf(stderr, "ERROR: Input matrix array is NULL\n");
        return NULL;
    }
    matrix *R = m_init_l(A->m, A->n, A->data);
    return R;
}

void m_printf(matrix *A) {
    if (!A) {
        fprintf(stderr, "ERROR: Input matrix array is NULL\n");
        return;
    }
    for (uint64_t i = 0; i < A->m; i++) {
        printf("[ ");
        for (uint64_t j = 0; j < A->n; j++) {
            complex z = A->data[i * A->n + j];
            c_printf(z);
            if (j < A->n - 1)
                printf(", ");
        }
        printf(" ]\n");
    }
}

/*
 *
 *   SWAP FUNCTION
 *
 */

void matrix_to_complex(matrix *A, complex *c, uint64_t *m);
vector *matrix_to_vector(matrix *A);
tensor *matrix_to_tensor(matrix *A);

/*
 *
 *   OPERATION FUNCTION
 *
 */

matrix *m_add(matrix *A, matrix *B, matrix *Result) {
    if (!A || !B) {
        fprintf(stderr, "ERROR: Input matrix is NULL\n");
        return NULL;
    }
    if (A->m != B->m || A->n != B->n) {
        fprintf(stderr, "ERROR: Matrix dimensions are incompatible for addition\n");
        return NULL;
    }

    if (Result) {
        if (Result->m != A->m || Result->n != A->n) {
            fprintf(stderr, "ERROR: Result matrix has incompatible dimensions\n");
            return NULL;
        }
    } else {
        Result = m_init(A->m, A->n);
        if (!Result)
            return NULL;
    }

    for (uint64_t i = 0; i < A->m * A->n; i++) {
        Result->data[i] = c_add(A->data[i], B->data[i]);
    }
    return Result;
}

matrix *m_sub(matrix *A, matrix *B, matrix *Result) {
    if (!A || !B) {
        fprintf(stderr, "ERROR: Input matrix is NULL\n");
        return NULL;
    }
    if (A->m != B->m || A->n != B->n) {
        fprintf(stderr, "ERROR: Matrix dimensions are incompatible for substraction\n");
        return NULL;
    }

    if (Result) {
        if (Result->m != A->m || Result->n != A->n) {
            fprintf(stderr, "ERROR: Result matrix has incompatible dimensions\n");
            return NULL;
        }
    } else {
        Result = m_init(A->m, A->n);
        if (!Result)
            return NULL;
    }

    for (uint64_t i = 0; i < A->m * A->n; i++) {
        Result->data[i] = c_sub(A->data[i], B->data[i]);
    }
    return Result;
}

matrix *m_mult_e(matrix *A, matrix *B, matrix *Result) {
    if (!A || !B) {
        fprintf(stderr, "ERROR: Input matrix is NULL\n");
        return NULL;
    }
    if (A->m != B->m || A->n != B->n) {
        fprintf(stderr,
                "ERROR: Matrix dimensions are incompatible for multiplication_E\n");
        return NULL;
    }

    if (Result) {
        if (Result->m != A->m || Result->n != A->n) {
            fprintf(stderr, "ERROR: Result matrix has incompatible dimensions\n");
            return NULL;
        }
    } else {
        Result = m_init(A->m, A->n);
        if (!Result)
            return NULL;
    }

    for (uint64_t i = 0; i < A->m * A->n; i++) {
        Result->data[i] = c_mult(A->data[i], B->data[i]);
    }
    return Result;
}

matrix *m_div_e(matrix *A, matrix *B, matrix *Result) {
    if (!A || !B) {
        fprintf(stderr, "ERROR: Input matrix is NULL\n");
        return NULL;
    }
    if (A->m != B->m || A->n != B->n) {
        fprintf(stderr, "ERROR: Matrix dimensions are incompatible for division_E\n");
        return NULL;
    }

    if (Result) {
        if (Result->m != A->m || Result->n != A->n) {
            fprintf(stderr, "ERROR: Result matrix has incompatible dimensions\n");
            return NULL;
        }
    } else {
        Result = m_init(A->m, A->n);
        if (!Result)
            return NULL;
    }

    for (uint64_t i = 0; i < A->m * A->n; i++) {
        Result->data[i] = c_div(A->data[i], B->data[i]);
    }
    return Result;
}

matrix *m_scale_r(matrix *A, double alpha, matrix *Result) {
    if (!A) {
        fprintf(stderr, "ERROR: Input matrix is NULL\n");
        return NULL;
    }
    if (Result) {
        if (Result->m != A->m || Result->n != A->n) {
            fprintf(stderr, "ERROR: Result matrix has incompatible dimensions\n");
            return NULL;
        }
    } else {
        Result = m_init(A->m, A->n);
        if (!Result)
            return NULL;
    }
    for (uint64_t i = 0; i < A->m * A->n; i++) {
        Result->data[i] = c_scale(A->data[i], alpha);
    }
    return Result;
}

matrix *m_scale(matrix *A, complex alpha, matrix *Result) {
    if (!A) {
        fprintf(stderr, "ERROR: Input matrix is NULL\n");
        return NULL;
    }
    if (Result) {
        if (Result->m != A->m || Result->n != A->n) {
            fprintf(stderr, "ERROR: Result matrix has incompatible dimensions\n");
            return NULL;
        }
    } else {
        Result = m_init(A->m, A->n);
        if (!Result)
            return NULL;
    }
    for (uint64_t i = 0; i < A->m * A->n; i++) {
        Result->data[i] = c_mult(A->data[i], alpha);
    }
    return Result;
}

matrix *m_conj(matrix *A, matrix *Result) {
    if (!A) {
        fprintf(stderr, "ERROR: Input matrix is NULL\n");
        return NULL;
    }
    if (Result) {
        if (Result->m != A->m || Result->n != A->n) {
            fprintf(stderr, "ERROR: Result matrix has incompatible dimensions\n");
            return NULL;
        }
    } else {
        Result = m_init(A->m, A->n);
        if (!Result)
            return NULL;
    }
    for (uint64_t i = 0; i < A->m * A->n; i++) {
        Result->data[i] = c_conj(A->data[i]);
    }
    return Result;
}

matrix *m_transpose(matrix *A, matrix *Result) {
    if (!A) {
        fprintf(stderr, "ERROR: Input matrix is NULL\n");
        return NULL;
    }
    if (Result) {
        if (Result->m != A->n || Result->n != A->m) {
            fprintf(stderr, "ERROR: Result matrix has incompatible dimensions\n");
            return NULL;
        }
    } else {
        Result = m_init(A->n, A->m);
        if (!Result)
            return NULL;
    }
    for (uint64_t i = 0; i < A->m; i++) {
        for (uint64_t j = 0; j < A->n; j++) {
            Result->data[j * A->m + i] = A->data[i * A->n + j];
        }
    }
    return Result;
}

matrix *m_hermitian(matrix *A, matrix *Result) {
    if (!A) {
        fprintf(stderr, "ERROR: Input matrix is NULL\n");
        return NULL;
    }
    if (Result) {
        if (Result->m != A->n || Result->n != A->m) {
            fprintf(stderr, "ERROR: Result matrix has incompatible dimensions\n");
            return NULL;
        }
    } else {
        Result = m_init(A->n, A->m);
        if (!Result)
            return NULL;
    }
    matrix *temp = m_conj(A, NULL);
    if (!temp)
        return NULL;
    m_transpose(temp, Result);
    m_free(temp);
    return Result;
}

double m_norm2(matrix *A) {
    if (!A) {
        fprintf(stderr, "ERROR: Input matrix is NULL\n");
        return -1.0;
    }
    double norm2 = 0;
    for (uint64_t i = 0; i < A->m * A->n; i++) {
        norm2 += c_norm2(A->data[i]);
    }
    return norm2;
}

double m_norm(matrix *A) {
    if (!A) {
        fprintf(stderr, "ERROR: Input matrix is NULL\n");
        return -1.0;
    }
    return sqrt(m_norm2(A));
}

matrix *m_mult(matrix *A, matrix *B, matrix *Result) {
    if (!A || !B) {
        fprintf(stderr, "ERROR: Input matrix is NULL\n");
        return NULL;
    }
    if (A->n != B->m) {
        fprintf(stderr, "ERROR: Matrix dimensions are incompatible for multiplication\n");
        return NULL;
    }
    if (Result) {
        if (Result->m != A->m || Result->n != B->n) {
            fprintf(stderr, "ERROR: Result matrix has incompatible dimensions\n");
            return NULL;
        }
    } else {
        Result = m_init(A->m, B->n);
        if (!Result)
            return NULL;
    }
    for (uint64_t i = 0; i < A->m; i++) {
        for (uint64_t j = 0; j < B->n; j++) {
            Result->data[i * B->n + j] = (complex){0.0, 0.0};
            for (uint64_t k = 0; k < A->n; k++) {
                Result->data[i * B->n + j] =
                    c_add(Result->data[i * B->n + j],
                          c_mult(A->data[i * A->n + k], B->data[k * B->n + j]));
            }
        }
    }
    return Result;
}
