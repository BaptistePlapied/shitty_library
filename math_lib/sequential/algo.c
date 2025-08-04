
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

double m_norm2_offdiag(matrix *M) {
    if (!M || M->m != M->n) {
        fprintf(stderr, "ERROR: Matrix is NULL or not square in m_norm2_offdiag\n");
        return -1.0;
    }

    double norm2 = 0.0;
    uint64_t n = M->n;

    for (uint64_t i = 0; i < n; i++) {
        for (uint64_t j = 0; j < n; j++) {
            if (i != j) {
                complex z = M->data[i * n + j];
                norm2 += z.Re * z.Re + z.Im * z.Im;
            }
        }
    }

    return norm2;
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

vector *forward_sub(matrix *L, vector *b, vector *result) {
    if (!L || !b) {
        fprintf(stderr, "ERROR: Input matrix or vector is NULL\n");
        return NULL;
    }
    if (L->m != b->m || L->n != L->m) {
        fprintf(stderr,
                "ERROR: Matrix and vector dimensions are incompatible for forward_sub\n");
        return NULL;
    }

    for (uint64_t i = 0; i < L->m; i++) {
        for (uint64_t j = i + 1; j < L->n; j++) {
            if (fabs(L->data[i * L->n + j].Re) > DBL_EPSILON ||
                fabs(L->data[i * L->n + j].Im) > DBL_EPSILON) {
                fprintf(stderr,
                        "Erreur : la matrice n'est pas triangulaire inférieure.\n");
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

    vector *y = v_copy(b);

    for (uint64_t i = 0; i < L->m; i++) {
        for (uint64_t j = 0; j < i; j++) {
            y->data[i] = c_sub(y->data[i], c_mult(L->data[i * L->n + j], y->data[j]));
        }

        complex pivot = L->data[i * L->n + i];
        if (fabs(pivot.Re) < DBL_EPSILON && fabs(pivot.Im) < DBL_EPSILON) {
            if (fabs(y->data[i].Re) > DBL_EPSILON || fabs(y->data[i].Im) > DBL_EPSILON) {
                fprintf(stderr, "Aucune solution\n");
                v_free(y);
                return NULL;
            } else {
                fprintf(stderr, "Infinité de solutions\n");
                v_free(y);
                return NULL;
            }
        }
        y->data[i] = c_div(y->data[i], pivot);
    }

    return y;
}

vector *back_sub(matrix *U, vector *b, vector *result) {
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

matrix *m_inv_lu(matrix *A, matrix *Result) {
    if (!A || A->m != A->n) {
        fprintf(stderr, "ERROR: Matrix must be square for inversion.\n");
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
    uint64_t n = A->n;
    matrix *L = m_init(n, n);
    matrix *U = m_init(n, n);
    matrix *P = m_init(n, n);
    if (!L || !U || !P) {
        fprintf(stderr, "ERROR: Memory allocation failed.\n");
        return NULL;
    }
    lu_pp(A, L, U, P, NULL);
    vector *b_perm = v_init(n);
    vector *y = v_init(n);
    vector *x = v_init(n);
    for (uint64_t col = 0; col < n; col++) {
        for (uint64_t i = 0; i < n; i++) {
            b_perm->data[i] =
                P->data[i * n + col]; // P * e_j donne la col-ème colonne de P
        }
        forward_sub(L, b_perm, y);
        back_sub(U, y, x);
        for (uint64_t i = 0; i < n; i++)
            Result->data[i * n + col] = x->data[i];
    }
    v_free(b_perm);
    v_free(y);
    v_free(x);
    m_free(L);
    m_free(U);
    m_free(P);
    return Result;
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

void m_mgs(matrix *A, matrix *Q, matrix *R) {
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
    vector *v_j = v_init(m);
    for (uint64_t j = 0; j < n; j++) {
        for (uint64_t i = 0; i < m; i++) {
            v_j->data[i] = A->data[i * n + j];
        }
        for (uint64_t i = 0; i < j; i++) {
            complex r_ij = c_zero();
            for (uint64_t k = 0; k < m; k++) {
                r_ij = c_add(c_mult(Q->data[k * n + i], v_j->data[k]), r_ij);
            }
            R->data[i * n + j] = r_ij;
            for (uint64_t k = 0; k < m; k++) {
                v_j->data[k] = c_sub(v_j->data[k], c_mult(r_ij, Q->data[k * n + i]));
            }
        }
        complex r_jj = c_zero();
        for (uint64_t i = 0; i < m; i++) {
            r_jj = c_add(c_mult(v_j->data[i], v_j->data[i]), r_jj);
        }
        r_jj = c_sqrt(r_jj);
        R->data[j * n + j] = r_jj;
        for (uint64_t i = 0; i < m; i++) {
            Q->data[i * n + j] = c_div(v_j->data[i], r_jj);
        }
    }
    v_free(v_j);
}

vector *v_householder(vector *a, vector *result) {
    if (!a)
        return NULL;
    if (result) {
        if (result->m != a->m)
            return NULL;
    } else {
        result = v_init(a->m);
        if (!result)
            return NULL;
    }
    double norm_a = v_norm(a);
    if (norm_a < 1e-14) {
        v_free(result);
        result = v_init_0(a->m);
        result->data[0] = c_real(1.0);
        return result;
    }
    complex x0 = a->data[0];
    double sign_real = (x0.Re >= 0) ? 1.0 : -1.0;
    vector *e1 = v_init_0(a->m);
    e1->data[0] = (complex){sign_real * norm_a, 0.0};
    for (uint64_t i = 0; i < a->m; i++) {
        result->data[i] = c_add(a->data[i], e1->data[i]);
    }
    double norm_v = v_norm(result);
    if (norm_v < 1e-14)
        norm_v = 1.0;
    for (uint64_t i = 0; i < a->m; i++) {
        result->data[i] = c_div(result->data[i], c_real(norm_v));
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
        matrix *VVᴴ = v_outer(v, v, NULL);
        m_printf(VVᴴ);
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
        matrix *H_k_T = m_track_in(m_hermitian(H_k, NULL), group_id);
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

vector *qr_eigenvalues(matrix *A, double tol, int max_iter, vector *result) {
    if (!A || A->m != A->n) {
        fprintf(stderr, "ERROR : Wrong input for qr_eigenvalues\n");
        return NULL;
    }
    if (result) {
        if (result->m != A->m) {
            fprintf(stderr, "ERROR : incompatible dim for result eigen\n");
            return NULL;
        }
    } else {
        result = v_init(A->m);
        if (!result)
            return NULL;
    }
    matrix *R = m_init(A->m, A->n);
    matrix *Q = m_init(A->m, A->n);
    matrix *H = m_init(A->m, A->n);
    m_copy_on(A, H);
    matrix *I = m_init(A->m, A->n);
    m_id(I);
    matrix *H_shifted = m_init(A->m, A->n);

    for (uint64_t i = 0; i < max_iter; i++) {
        complex shift = H->data[(A->n - 1) * A->n + (A->n - 1)];

        m_copy_on(H, H_shifted);
        for (uint64_t i = 0; i < A->n; i++) {
            H_shifted->data[i * A->n + i] = c_sub(H_shifted->data[i * A->n + i], shift);
        }
        m_mgs(H_shifted, Q, R);
        m_mult(R, Q, H);
        for (uint64_t i = 0; i < A->n; i++) {
            H->data[i * A->n + i] = c_add(H->data[i * A->n + i], shift);
        }
        if (m_norm2_offdiag(H) < tol * tol)
            break;
    }
    for (uint64_t i = 0; i < A->n; i++) {
        result->data[i] = H->data[i * A->n + i];
    }
    printf("H : \n");
    m_printf(H);
    m_free(R);
    m_free(Q);
    m_free(H);
    m_free(I);
    m_free(H_shifted);
    return result;
}

vector **qr_eigenvectors(matrix *A, double tol, int max_iter, vector **result) {
    if (!A || A->m != A->n) {
        fprintf(stderr, "ERROR : Wrong input for qr_eigenvectors\n");
        return NULL;
    }
    if (result) {
        for (uint64_t i = 0; i < A->m; i++) {
            if (result[i]->m != A->m) {
                fprintf(stderr, "ERROR : incompatible dim for result eigen\n");
                return NULL;
            }
        }
    } else {
        result = malloc(A->m * sizeof(vector *));
        for (uint64_t i = 0; i < A->m; i++) {
            result[i] = v_init(A->m);
            if (!result[i]) {
                fprintf(stderr, "ERROR : Failed to allocate result eigen vector\n");
                return NULL;
            }
        }
    }
    matrix *R = m_init(A->m, A->n);
    matrix *Q = m_init(A->m, A->n);
    matrix *H = m_init(A->m, A->n);
    m_copy_on(A, H);
    matrix *V = m_init(A->m, A->n);
    m_id(V);
    matrix *V_temp = m_init(A->m, A->n);
    matrix *I = m_init(A->m, A->n);
    m_id(I);
    matrix *H_shifted = m_init(A->m, A->n);

    for (uint64_t i = 0; i < max_iter; i++) {
        complex shift = H->data[(A->n - 1) * A->n + (A->n - 1)];

        m_copy_on(H, H_shifted);
        for (uint64_t i = 0; i < A->n; i++) {
            H_shifted->data[i * A->n + i] = c_sub(H_shifted->data[i * A->n + i], shift);
        }
        m_mgs(H_shifted, Q, R);
        m_mult(R, Q, H);
        for (uint64_t i = 0; i < A->n; i++) {
            H->data[i * A->n + i] = c_add(H->data[i * A->n + i], shift);
        }
        if (m_norm2_offdiag(H) < tol * tol)
            break;
        // Update eigenvectors
        m_mult(V, Q, V_temp);
        m_copy_on(V_temp, V);
    }
    for (uint64_t j = 0; j < A->n; j++) {
        for (uint64_t i = 0; i < A->m; i++) {
            result[j]->data[i] = V->data[i * A->n + j];
        }
    }
    m_free(R);
    m_free(Q);
    m_free(H);
    m_free(V);
    m_free(V_temp);
    m_free(I);
    m_free(H_shifted);
    return result;
}

matrix *m_exp_pade(matrix *A, complex alpha, matrix *Result) {
    if (!A || A->m != A->n) {
        fprintf(stderr, "ERROR: Input matrix is NULL or not square.\n");
        return NULL;
    }
    uint64_t n = A->n;
    // Coefficients for Pade-13
    double c[] = {64764752532480000.0,
                  32382376266240000.0,
                  7771770303897600.0,
                  1187353796428800.0,
                  129060195264000.0,
                  10559470521600.0,
                  670442572800.0,
                  33522128640.0,
                  1323241920.0,
                  40840800.0,
                  960960.0,
                  16380.0,
                  182.0,
                  1.0};
    // Scale A by alpha
    matrix *A_scaled = m_init(n, n);
    for (uint64_t i = 0; i < n * n; i++)
        A_scaled->data[i] = c_mult(A->data[i], alpha);
    // Compute norm for scaling
    double norm = m_norm(A_scaled);
    int s = 0;
    if (norm > 0.5) {
        s = (int)fmax(0, ceil(log2(norm / 0.5)));
        complex inv_s = c_real(1.0 / pow(2.0, s));
        m_scale(A_scaled, inv_s, A_scaled);
    }
    // Precompute powers of A
    matrix *A2 = m_mult(A_scaled, A_scaled, NULL);
    matrix *A4 = m_mult(A2, A2, NULL);
    matrix *A6 = m_mult(A2, A4, NULL);
    // U = A * (A6 * c13 + A4 * c11 + A2 * c9 + I * c7)
    matrix *U = m_init(n, n);
    matrix *term = m_init(n, n);
    m_scale_r(A6, c[13], U);
    m_scale_r(A4, c[11], term);
    m_add(U, term, U);
    m_scale_r(A2, c[9], term);
    m_add(U, term, U);
    matrix *I = m_init(n, n);
    m_id(I);
    m_scale_r(I, c[7], term);
    m_add(U, term, U);
    m_mult(A_scaled, U, U); // U = A * U
    // V = A6 * c12 + A4 * c10 + A2 * c8 + I * c6
    matrix *V = m_init(n, n);
    m_scale_r(A6, c[12], V);
    m_scale_r(A4, c[10], term);
    m_add(V, term, V);
    m_scale_r(A2, c[8], term);
    m_add(V, term, V);
    m_scale_r(I, c[6], term);
    m_add(V, term, V);
    // P = V + U
    // Q = V - U
    matrix *P = m_add(V, U, NULL);
    matrix *Q = m_sub(V, U, NULL);

    // Compute (Q)^-1
    matrix *Qinv = m_inv_lu(Q, NULL);
    if (!Qinv) {
        fprintf(stderr, "ERROR: Matrix inversion failed.\n");
        // Free and return
        m_free(A_scaled);
        m_free(A2);
        m_free(A4);
        m_free(A6);
        m_free(U);
        m_free(V);
        m_free(I);
        m_free(term);
        m_free(P);
        m_free(Q);
        return NULL;
    }

    // Final result = Qinv * P
    Result = m_mult(Qinv, P, Result);

    // Undo scaling
    for (int i = 0; i < s; i++) {
        Result = m_mult(Result, Result, Result);
    }

    // Cleanup
    m_free(A_scaled);
    m_free(A2);
    m_free(A4);
    m_free(A6);
    m_free(U);
    m_free(V);
    m_free(I);
    m_free(term);
    m_free(P);
    m_free(Q);
    m_free(Qinv);

    return Result;
}
