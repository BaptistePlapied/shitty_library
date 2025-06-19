
#include "../headers/general.h"
#include "../headers/sequential.h"
#include <float.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void row_op_printf(row_op ope) {
    printf(" [%llu -> %llu] : ", ope.row_source, ope.row_target);
    c_printf(ope.multiplier);
}

void echelon_expr_printf(echelon_expr echelon) {
    printf(" echelon : \n\tv_size : %llu\n\tnum_ops : %llu\n\tpivots :\n\t\t",
           echelon.v_size, echelon.num_ops);
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
                fprintf(stderr,
                        "Erreur : la matrice n'est pas triangulaire supérieure.\n");
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
        if (!Result)
            return NULL;
    }

    if (expr) {
        expr->pivots = v_init(A->m);
        expr->permutation = v_init(A->m);
        expr->operations = malloc((A->m * A->m) / 2 * sizeof(row_op));
    } else {
        fprintf(stderr, "ERROR: Input expression is NULL\n");
        return NULL;
    }
    expr->v_size = A->m;
    uint64_t op_index = 0;
    for (uint64_t i = 0; i < A->m; ++i) {
        complex pivot = Result->data[i * A->n + i];
        expr->pivots->data[i] = pivot;
        for (uint64_t j = i; j < A->n; ++j) {
            Result->data[i * A->n + j] = c_div(Result->data[i * A->n + j], pivot);
        }
        for (uint64_t k = i + 1; k < A->m; ++k) {
            complex factor = Result->data[k * A->n + i];
            expr->operations[op_index++] =
                (row_op){.row_target = k, .row_source = i, .multiplier = factor};
            for (uint64_t j = i; j < A->n; ++j) {
                complex sub = c_mult(factor, Result->data[i * A->n + j]);
                Result->data[k * A->n + j] = c_sub(Result->data[k * A->n + j], sub);
            }
        }
    }
    expr->num_ops = op_index;
    return Result;
}

vector *v_apply_echelon(vector *a, echelon_expr *expr, vector *result) {
    if (!a || !expr) {
        fprintf(stderr, "ERROR: Input vector or expression is NULL\n");
        return NULL;
    }
    if (result) {
        if (result->m != a->m) {
            fprintf(stderr, "ERROR: Result vector has incompatible dimensions\n");
            return NULL;
        }
        result = v_copy(a);
    } else {
        result = v_init_l(a->m, a->data);
        if (!result)
            return NULL;
    }
    for (uint64_t i = 0; i < expr->v_size; i++) {
        result->data[i] = c_div(result->data[i], expr->pivots->data[i]);
        for (uint64_t j = 0; j < expr->num_ops; j++) {
            if (expr->operations[j].row_source == i) {
                uint64_t k = expr->operations[j].row_target;
                complex factor = expr->operations[j].multiplier;
                result->data[k] = c_sub(result->data[k], c_mult(factor, result->data[i]));
            } else if (expr->operations[j].row_source > i) {
                break;
            }
        }
    }
    return result;
}

matrix *m_id(matrix *Result) {
    if (!Result) {
        fprintf(stderr, "ERROR: Input matrix is NULL\n");
        return NULL;
    }
    for (uint64_t i = 0; i < Result->m; i++) {
        for (uint64_t j = 0; j < Result->m; j++) {
            if (i == j) {
                Result->data[i * Result->m + j] = (complex){1.0, 0.0};
            } else {
                Result->data[i * Result->m + j] = (complex){0.0, 0.0};
            }
        }
    }
    return Result;
}

void lu_pp(matrix *A, matrix *L, matrix *U, matrix *P, uint64_t *swap) {
    if (!A || !U || !L || !P) {
        fprintf(stderr, "ERROR: Input matrix is NULL\n");
        return;
    }
    if (U->m != A->m || U->n != A->n) {
        fprintf(stderr, "ERROR: U matrix has incompatible dimensions\n");
        return;
    } else if (L->m != A->m || L->n != A->n) {
        fprintf(stderr, "ERROR: L matrix has incompatible dimensions\n");
        return;
    } else if (P->m != A->m || P->n != A->n) {
        fprintf(stderr, "ERROR: P matrix has incompatible dimensions\n");
        return;
    }
    uint64_t n = A->m;
    memcpy(U->data, A->data, sizeof(complex) * A->m * A->n);
    m_id(P);
    uint64_t swap_count = 0;
    for (uint64_t i = 0; i < n; i++) {
        uint64_t max_row = i;
        double max_val = c_norm2(U->data[i * n + i]);
        for (uint64_t k = i + 1; k < n; k++) {
            double val = c_norm2(U->data[k * n + i]);
            if (val > max_val) {
                max_val = val;
                max_row = k;
            }
        }
        if (max_row != i) {
            m_swap_rows(U, i, max_row);
            m_swap_rows(P, i, max_row);
            for (uint64_t j = 0; j < i; j++) {
                complex temp = L->data[i * n + j];
                L->data[i * n + j] = L->data[max_row * n + j];
                L->data[max_row * n + j] = temp;
            }
            swap_count++;
        }
        for (uint64_t j = i + 1; j < n; j++) {
            complex mult = c_div(U->data[j * n + i], U->data[i * n + i]);
            L->data[j * n + i] = mult;
            for (uint64_t k = i; k < n; k++) {
                complex prod = c_mult(mult, U->data[i * n + k]);
                U->data[j * n + k] = c_sub(U->data[j * n + k], prod);
            }
        }
    }
    for (uint64_t i = 0; i < n; i++) {
        L->data[i * n + i] = (complex){1.0, 0.0};
    }
    if (swap) {
        *swap = swap_count;
    }
}

complex m_tr_det(matrix *A) {
    if (!A) {
        fprintf(stderr, "ERROR: Input matrix is NULL\n");
        return (complex){-1.0, -1.0};
    }
    if (A->m != A->n) {
        fprintf(stderr, "ERROR: matrix must be squaric for det\n");
        return (complex){-1.0, -1.0};
    }
    complex result = (complex){1.0, 0.0};
    for (uint64_t i = 0; i < A->m; i++) {
        result = c_mult(result, A->data[i * A->n + i]);
    }
    return result;
}

complex m_det(matrix *A) {
    if (!A) {
        fprintf(stderr, "ERROR: Input matrix is NULL\n");
        return (complex){-1.0, -1.0};
    }
    if (A->m != A->n) {
        fprintf(stderr, "ERROR: matrix must be squaric for det\n");
        return (complex){-1.0, -1.0};
    }
    matrix *L = m_init(A->m, A->n);
    matrix *U = m_init(A->m, A->n);
    matrix *P = m_init(A->m, A->n);
    uint64_t swap;
    lu_pp(A, L, U, P, &swap);
    complex det = m_tr_det(U);
    if (swap % 2 == 1) {
        det = c_mult((complex){-1.0, 0.0}, det);
    }
    m_free(L);
    m_free(U);
    m_free(P);
    return det;
}

void m_mgs_qr(matrix *A, matrix *Q, matrix *R) {
    if (!A || !Q || !R) {
        fprintf(stderr, "ERROR: Input matrix is NULL\n");
        return;
    }
    if (A->m != A->n) {
        fprintf(stderr, "ERROR: A matrix has bad dimensions for QR\n");
        return;
    }
    if (Q->m != A->m || Q->n != A->n) {
        fprintf(stderr, "ERROR: Q matrix has incompatible dimensions\n");
        return;
    }
    if (R->m != A->n || R->n != A->n) {
        fprintf(stderr, "ERROR: R matrix must be square with size A->n\n");
        return;
    }

    uint64_t m = A->m;
    uint64_t n = A->n;

    vector *vj = v_init(m);
    vector *qi = v_init(m);
    vector *qj = v_init(m);

    for (uint64_t j = 0; j < n; ++j) {
        // Copy column j of A into vj
        for (uint64_t i = 0; i < m; ++i)
            vj->data[i] = A->data[i * n + j];
        for (uint64_t i = 0; i < j; ++i) {
            // Copy column i of Q into qi
            for (uint64_t k = 0; k < m; ++k)
                qi->data[k] = Q->data[k * n + i];
            // R[i][j] = <q_i, v_j>
            R->data[i * n + j] = v_dot_prod(qi, vj);
            // vj -= R[i][j] * q_i
            vj = v_sub(vj, v_scale(qi, R->data[i * n + j], NULL), vj);
        }
        // R[j][j] = ||vj||
        complex rjj = c_sqrt(v_dot_prod(vj, vj));
        if (c_norm2(rjj) == 0.0) {
            fprintf(stderr, "ERROR: Linearly dependent columns at column %llu\n", j);
            v_free(vj);
            v_free(qi);
            v_free(qj);
            return;
        }
        R->data[j * n + j] = rjj;
        v_scale(vj, c_div((complex){1.0, 0.0}, rjj), qj);
        // Store qj into Q
        for (uint64_t i = 0; i < m; ++i)
            Q->data[i * n + j] = qj->data[i];
    }
    v_free(vj);
    v_free(qi);
    v_free(qj);
}

vector *v_householder(vector *a, vector *result) {
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
    vector *e1 = v_init_0(a->m);
    e1->data[0] = (complex){1.0, 0.0};
    double norm_a = v_norm(a);
    for (uint64_t i = 0; i < a->m; i++) {
        result->data[i] = c_add(a->data[i], c_mult(e1->data[i], (complex){norm_a, 0.0}));
    }
    double norm_res = v_norm(result);
    for (uint64_t i = 0; i < a->m; i++) {
        result->data[i] = c_div(result->data[i], (complex){norm_res, 0.0});
    }
    v_free(e1);
    return result;
}

void m_hessen_reduc(matrix *A, matrix *A_res, matrix *Q_tot) {
    if (!A || !A_res || !Q_tot) {
        fprintf(stderr, "ERROR: One or more input matrices are NULL\n");
        return;
    }

    if (A_res->m != A->m || A_res->n != A->n) {
        fprintf(stderr, "ERROR: A_res matrix has incompatible dimensions\n");
        return;
    }

    if (Q_tot->m != A->n || Q_tot->n != A->n) {
        fprintf(stderr, "ERROR: Q_tot matrix has incompatible dimensions\n");
        return;
    }

    uint64_t group_id = find_empty_group_id();
    m_copy_on(A, A_res);

    // Initialize Q_tot as identity
    for (uint64_t i = 0; i < A->m; i++) {
        for (uint64_t j = 0; j < A->m; j++) {
            Q_tot->data[i * A->n + j] = (i == j) ? c_real(1.0) : c_zero();
        }
    }

    for (uint64_t k = 0; k < A->n - 2; k++) {
        uint64_t sub_n = A->n - k - 1;

        // x = A_res[k+1:, k]
        vector *x = v_init(sub_n);
        for (uint64_t i = 0; i < sub_n; i++) {
            x->data[i] = A_res->data[(k + 1 + i) * A_res->n + k];
        }

        // v = Householder(x)
        vector *v = v_householder(x, NULL);
        matrix *vv2 = m_scale_r(m_track_in(v_outer(v, v, NULL), group_id), 2.0, NULL);

        // Construct H_sub = I - 2vvᵗ
        matrix *H_sub = m_init(sub_n, sub_n);
        for (uint64_t i = 0; i < sub_n; i++) {
            for (uint64_t j = 0; j < sub_n; j++) {
                H_sub->data[i * sub_n + j] = (i == j) ? c_real(1.0) : c_zero();
                H_sub->data[i * sub_n + j] =
                    c_sub(H_sub->data[i * sub_n + j], vv2->data[i * vv2->n + j]);
            }
        }

        // Embed H_sub into H_k (full-size)
        matrix *H_k = m_init(A->n, A->n);
        for (uint64_t i = 0; i < A->n; i++) {
            for (uint64_t j = 0; j < A->n; j++) {
                H_k->data[i * A->n + j] = (i == j) ? c_real(1.0) : c_zero();
            }
        }

        for (uint64_t i = 0; i < sub_n; i++) {
            for (uint64_t j = 0; j < sub_n; j++) {
                H_k->data[(k + 1 + i) * A->n + (k + 1 + j)] =
                    H_sub->data[i * H_sub->n + j];
            }
        }

        // A_res ← H_k A_res H_kᵗ
        matrix *H_k_T = m_track_in(m_transpose(H_k, NULL), group_id);
        matrix *A_temp =
            m_mult(H_k, m_track_in(m_mult(A_res, H_k_T, NULL), group_id), NULL);
        m_copy_on(A_temp, A_res);

        // Q_tot ← Q_tot × H_k
        matrix *Q_temp = m_mult(Q_tot, H_k, NULL);
        m_copy_on(Q_temp, Q_tot);

        // Cleanup
        m_free(A_temp);
        m_free(Q_temp);
        v_free(x);
        v_free(v);
        m_free(H_sub);
        m_free(H_k);
    }

    track_group_clear(group_id);
}
