#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../headers/general.h"
#include "../headers/sequential.h"

/*
 *
 * OPERATION FUNCTION
 *
 */

vector *v_init(uint64_t m) {
    vector *a = (vector *)malloc(sizeof(vector));
    if (!a) {
        fprintf(stderr, "ERROR: Failed to allocate vector struct\n");
        return NULL;
    }
    a->m = m;
    a->data = (complex *)malloc(sizeof(complex) * m);
    if (!a->data) {
        fprintf(stderr, "ERROR: Failed to allocate vector data\n");
        free(a);
        return NULL;
    }
    return a;
}

vector *v_init_0(uint64_t m) {
    vector *a = (vector *)malloc(sizeof(vector));
    if (!a) {
        fprintf(stderr, "ERROR: Failed to allocate vector struct\n");
        return NULL;
    }
    a->m = m;
    a->data = (complex *)calloc(m, sizeof(complex));
    if (!a->data) {
        fprintf(stderr, "ERROR: Failed to allocate vector data\n");
        free(a);
        return NULL;
    }
    return a;
}

vector *v_init_l(uint64_t m, complex *c) {
    if (!c) {
        fprintf(stderr, "ERROR: Input complex array is NULL\n");
        return NULL;
    }
    vector *a = (vector *)malloc(sizeof(vector));
    if (!a) {
        fprintf(stderr, "ERROR: Failed to allocate vector struct\n");
        return NULL;
    }
    a->m = m;
    a->data = (complex *)malloc(sizeof(complex) * m);
    if (!a->data) {
        fprintf(stderr, "ERROR: Failed to allocate vector data\n"); // Not "complex array
                                                                    // // does not exist"
        free(a);
        return NULL;
    }
    memcpy(a->data, c, sizeof(complex) * m);
    return a;
}

void v_free(vector *a) {
    if (a) {
        free(a->data);
        free(a);
    }
}
vector *v_copy(vector *a) {
    if (!a) {
        fprintf(stderr, "ERROR: Input vector is NULL\n");
        return NULL;
    }
    return v_init_l(a->m, a->data);
}

vector *v_resize(vector *a, uint64_t new_m) {
    if (!a) {
        fprintf(stderr, "ERROR: Input vector is NULL\n");
        return NULL;
    }
    complex *new_data = (complex *)realloc(a->data, sizeof(complex) * new_m);
    if (!new_data) {
        fprintf(stderr, "ERROR: Failed to allocate memory during vector resize\n");
        return NULL;
    }
    if (new_m > a->m) {
        for (uint64_t i = a->m; i < new_m; ++i) {
            new_data[i] = (complex){0.0, 0.0};
        }
    }
    a->data = new_data;
    a->m = new_m;
    return a;
}

void v_printf(vector *a) {
    if (!a) {
        fprintf(stderr, "ERROR: Input vector is NULL\n");
        return;
    }
    printf("[ ");
    for (uint64_t i = 0; i < a->m; i++) {
        complex z = a->data[i];
        c_printf(z);
        if (i < a->m - 1)
            printf(", ");
    }
    printf(" ]\n");
}

/*
 *
 * SWAP STRUCT
 *
 */
void vector_to_complex(vector *a, complex *c, uint64_t *m);
matrix *vector_to_matrix(vector *a, bool as_col);
tensor *vector_to_tensor(vector *a);

/*
 *
 * OPERATION FUNCTION
 *
 */

vector *v_add(vector *a, vector *b, vector *result) {
    if (!a || !b) {
        fprintf(stderr, "ERROR: Input vector is NULL\n");
        return NULL;
    }
    if (a->m != b->m) {
        fprintf(stderr, "ERROR: Vectors dimensions are incompatible for addition\n");
        return NULL;
    }
    if (result) {
        if (result->m != a->m) {
            fprintf(stderr, "ERROR: Result vector has incompatible dimension\n");
            return NULL;
        }
    } else {
        result = v_init(a->m);
        if (!result)
            return NULL;
    }
    for (uint64_t i = 0; i < a->m; i++) {
        result->data[i] = c_add(a->data[i], b->data[i]);
    }
    return result;
}

vector *v_sub(vector *a, vector *b, vector *result) {
    if (!a || !b) {
        fprintf(stderr, "ERROR: Input vector is NULL\n");
        return NULL;
    }
    if (a->m != b->m) {
        fprintf(stderr, "ERROR: Vectors dimensions are incompatible for substraction\n");
        return NULL;
    }
    if (result) {
        if (result->m != a->m) {
            fprintf(stderr, "ERROR: Result vector has incompatible dimension\n");
            return NULL;
        }
    } else {
        result = v_init(a->m);
        if (!result)
            return NULL;
    }
    for (uint64_t i = 0; i < a->m; i++) {
        result->data[i] = c_sub(a->data[i], b->data[i]);
    }
    return result;
}

vector *v_mult_e(vector *a, vector *b, vector *result) {
    if (!a || !b) {
        fprintf(stderr, "ERROR: Input vector is NULL\n");
        return NULL;
    }
    if (a->m != b->m) {
        fprintf(stderr,
                "ERROR: Vectors dimensions are incompatible for multiplication_E\n");
        return NULL;
    }
    if (result) {
        if (result->m != a->m) {
            fprintf(stderr, "ERROR: Result vector has incompatible dimension\n");
            return NULL;
        }
    } else {
        result = v_init(a->m);
        if (!result)
            return NULL;
    }
    for (uint64_t i = 0; i < a->m; i++) {
        result->data[i] = c_mult(a->data[i], b->data[i]);
    }
    return result;
}

vector *v_div_e(vector *a, vector *b, vector *result) {
    if (!a || !b) {
        fprintf(stderr, "ERROR: Input vector is NULL\n");
        return NULL;
    }
    if (a->m != b->m) {
        fprintf(stderr, "ERROR: Vectors dimensions are incompatible for division_E\n");
        return NULL;
    }
    if (result) {
        if (result->m != a->m) {
            fprintf(stderr, "ERROR: Result vector has incompatible dimension\n");
            return NULL;
        }
    } else {
        result = v_init(a->m);
        if (!result)
            return NULL;
    }
    for (uint64_t i = 0; i < a->m; i++) {
        result->data[i] = c_div(a->data[i], b->data[i]);
    }
    return result;
}

vector *v_scale(vector *a, complex alpha, vector *result) {
    if (!a) {
        fprintf(stderr, "ERROR: Input vector is NULL\n");
        return NULL;
    }
    if (result) {
        if (result->m != a->m) {
            fprintf(stderr, "ERROR: Result vector has incompatible dimension\n");
            return NULL;
        }
    } else {
        result = v_init(a->m);
        if (!result)
            return NULL;
    }
    for (uint64_t i = 0; i < a->m; i++) {
        result->data[i] = c_mult(a->data[i], alpha);
    }
    return result;
}

vector *v_scale_r(vector *a, double alpha, vector *result) {
    if (!a) {
        fprintf(stderr, "ERROR: Input vector is NULL\n");
        return NULL;
    }
    if (result) {
        if (result->m != a->m) {
            fprintf(stderr, "ERROR: Result vector has incompatible dimension\n");
            return NULL;
        }
    } else {
        result = v_init(a->m);
        if (!result)
            return NULL;
    }
    for (uint64_t i = 0; i < a->m; i++) {
        result->data[i] = c_scale(a->data[i], alpha);
    }
    return result;
}

vector *v_conj(vector *a, vector *result) {
    if (!a) {
        fprintf(stderr, "ERROR: Input vector is NULL\n");
        return NULL;
    }
    if (result) {
        if (result->m != a->m) {
            fprintf(stderr, "ERROR: Result vector has incompatible dimension\n");
            return NULL;
        }
    } else {
        result = v_init(a->m);
        if (!result)
            return NULL;
    }
    for (uint64_t i = 0; i < a->m; i++) {
        result->data[i] = c_conj(a->data[i]);
    }
    return result;
}

vector *v_Re(vector *a, vector *result) {
    if (!a) {
        fprintf(stderr, "ERROR: Input vector is NULL\n");
        return NULL;
    }
    if (result) {
        if (result->m != a->m) {
            fprintf(stderr, "ERROR: Result vector has incompatible dimension\n");
            return NULL;
        }
    } else {
        result = v_init(a->m);
        if (!result)
            return NULL;
    }
    for (uint64_t i = 0; i < a->m; i++) {
        result->data[i] = (complex){a->data[i].Re, 0.0};
    }
    return result;
}

vector *v_Im(vector *a, vector *result) {
    if (!a) {
        fprintf(stderr, "ERROR: Input vector is NULL\n");
        return NULL;
    }
    if (result) {
        if (result->m != a->m) {
            fprintf(stderr, "ERROR: Result vector has incompatible dimension\n");
            return NULL;
        }
    } else {
        result = v_init(a->m);
        if (!result)
            return NULL;
    }
    for (uint64_t i = 0; i < a->m; i++) {
        result->data[i] = (complex){0.0, a->data[i].Im};
    }
    return result;
}

vector *v_map(vector *a, complex (*f)(complex), vector *result) {
    if (!a) {
        fprintf(stderr, "ERROR: Input vector is NULL\n");
        return NULL;
    }
    if (result) {
        if (result->m != a->m) {
            fprintf(stderr, "ERROR: Result vector has incompatible dimension\n");
            return NULL;
        }
    } else {
        result = v_init(a->m);
        if (!result)
            return NULL;
    }
    for (uint64_t i = 0; i < a->m; i++) {
        result->data[i] = f(a->data[i]);
    }
    return result;
}

uint8_t v_equal(vector *a, vector *b, vector *result) {
    if (!a || !b) {
        fprintf(stderr, "ERROR: Input vector is NULL\n");
        return 255;
    }
    if (a->m != b->m) {
        fprintf(stderr, "ERROR: Vectors dimensions are incompatible for equality\n");
        return 255;
    }
    for (uint64_t i = 0; i < a->m; i++) {
        if (c_equal(a->data[i], b->data[i]) != 1) {
            return 0;
        }
    }
    return 1;
}
complex v_dot_prod(vector *a, vector *b) {
    if (!a || !b) {
        fprintf(stderr, "ERROR: Input vector is NULL\n");
        return (complex){-1, -1};
    }
    if (a->m != b->m) {
        fprintf(stderr, "ERROR: Vectors dimensions are incompatible for equality\n");
        return (complex){-1, -1};
    }
    complex result = {0, 0};
    for (uint64_t i = 0; i < a->m; i++) {
        result = c_add(c_mult(a->data[i], b->data[i]), result);
    }
    return result;
}

double v_norm(vector *a) {
    if (!a) {
        fprintf(stderr, "ERROR: Input vector is NULL\n");
        return -1;
    }
    return v_norm2(a);
}
double v_norm2(vector *a) {
    if (!a) {
        fprintf(stderr, "ERROR: Input vector is NULL\n");
        return -1;
    }
    double result = 0;
    for (uint64_t i = 0; i < a->m; i++) {
        result += c_norm2(a->data[i]);
    }
    return result;
}
