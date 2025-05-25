
#include "../headers/general.h"
#include "../headers/sequential.h"
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

void row_op_printf(row_op ope) {
    printf(" [%llu -> %llu] : ", ope.row_source, ope.row_target);
    c_printf(ope.multiplier);
}

void echelon_expr_printf(echelon_expr echelon) {
    printf(" echelon : \n\tnum_ops : %llu\n\tpivots :\n\t\t", echelon.num_ops);
    v_printf(echelon.pivots);
    printf("\tpermutation :\n\t\t");
    v_printf(echelon.permutation);
    printf("\toperations :\n\t\t");
    for (uint64_t i = 0; i < echelon.num_ops; i++) {
        row_op_printf(echelon.operations[i]);
        printf("\n\t\t");
    }
    printf("\n");
}

vector *backsub(matrix *U, vector *b, vector *result) {
    if (!U || !b) {
        fprintf(stderr, "ERROR: Input matrix or vector is NULL\n");
        return NULL;
    }
    if (U->m != b->m || U->n != U->m) {
        fprintf(stderr,
                "ERROR: Matrix and vector dimensions are incompatible for backsub\n");
        return NULL;
    }
    for (uint64_t i = 1; i < U->m; i++) {
        for (uint64_t j = 0; j < i; j++) {
            if (fabs(U->data[i * U->n + j].Re) > DBL_EPSILON ||
                fabs(U->data[i * U->n + j].Im) > DBL_EPSILON) {
                fprintf(stderr, "Erreur : la matrice n'est pas triangulaire supérieure.\n");
                return NULL;
            }
        }
    }
    if (result) {
        if (result->m != b->m) {
            fprintf(stderr, "ERROR: Result vector has incompatible dimensions\n");
            return NULL;
        }
    } else {
        result = v_init(b->m);
        if (!result)
            return NULL;
    }
    vector *x = v_copy(b);
    for (int64_t i = (int64_t)U->m - 1; i >= 0; i--) {
        for (uint64_t j = i + 1; j < U->n; j++) {
            x->data[i] = c_sub(x->data[i], c_mult(U->data[i * U->n + j], x->data[j]));
        }
        complex pivot = U->data[i * U->n + i];
        if (fabs(pivot.Re) < DBL_EPSILON && fabs(pivot.Im) < DBL_EPSILON) {
            if (fabs(x->data[i].Re) > DBL_EPSILON || fabs(x->data[i].Im) > DBL_EPSILON) {
                fprintf(stderr, "Aucune solution\n");
                v_free(x);
                return NULL;
            } else {
                fprintf(stderr, "Infinité de solutions\n");
                v_free(x);
                return NULL;
            }
        }
        x->data[i] = c_div(x->data[i], pivot);
    }
    return x;
}

matrix *m_echelon(matrix *A, echelon_expr *expr, matrix *Result) {
    if (!A || !expr) {
        fprintf(stderr, "ERROR: NULL input to m_echelon\n");
        return NULL;
    }

    if (Result) {
        if (Result->m != A->m || Result->n != A->n) {
            fprintf(stderr, "ERROR: Result matrix has incompatible dimensions\n");
            return NULL;
        }
        Result = m_copy(A);
    } else {
        Result = m_init_l(A->m, A->n, A->data);
        if (!Result) return NULL;
    }

    if (expr) {
        expr->pivots = v_init(A->m);
        expr->permutation = v_init(A->m);
        expr->operations = malloc((A->m * A->m) / 2 * sizeof(row_op));
    } else {
        fprintf(stderr, "ERROR: Input expression is NULL\n");
        return NULL;
    }

    uint64_t op_index = 0;
    for (uint64_t i = 0; i < A->m; ++i) {
        complex pivot = Result->data[i * A->n + i];
        expr->pivots->data[i] = pivot;
        for (uint64_t j = i; j < A->n; ++j) {
            Result->data[i * A->n + j] = c_div(Result->data[i * A->n + j], pivot);
        }
        for (uint64_t k = i + 1; k < A->m; ++k) {
            complex factor = Result->data[k * A->n + i];
            expr->operations[op_index++] = (row_op){
                .row_target = k,
                .row_source = i,
                .multiplier = factor
            };
            for (uint64_t j = i; j < A->n; ++j) {
                complex sub = c_mult(factor, Result->data[i * A->n + j]);
                Result->data[k * A->n + j] = c_sub(Result->data[k * A->n + j], sub);
            }
        }
    }

    expr->num_ops = op_index;
    return Result;
}
