#include "the_header_file.h"

// TODO: make all the matrix/vect operation possible/needed
// TODO: optimze them as much as possible in my range of skills

// Matrix Operations List

void Mtx_update_active(Mtx *m) {
    m->active_num_row = 0;
    m->active_num_col = 0;
    for (size_t i = 0; i < DIM_MAT; i++) {
        uint8_t is_row_empty = 1;
        uint8_t is_col_empty = 1;
        for (size_t j = 0; j < DIM_MAT; j++) {
            if (m->mat[j * DIM_MAT + i] != 0.0) {
                is_row_empty = 0;
            }
            if (m->mat[i * DIM_MAT + j] != 0.0) {
                is_col_empty = 0;
            }
        }
        m->active_row[i] = !is_row_empty;
        m->active_col[i] = !is_col_empty;
        if (!is_row_empty) {
            m->active_num_row++;
        }
        if (!is_col_empty) {
            m->active_num_col++;
        }
    }
}

void Mtx_init(Mtx *m, size_t dimR, size_t dimC, float *mat) {
    m->dimR = dimR;
    m->dimC = dimC;
    if (mat == NULL) {
        for (size_t i = 0; i < DIM_MAT; i++) {
            for (size_t j = 0; j < DIM_MAT; j++) {
                m->mat[i * DIM_MAT + j] = 0.0;
            }
        }
    } else {
        for (size_t i = 0; i < DIM_MAT; i++) {
            for (size_t j = 0; j < DIM_MAT; j++) {
                if (i < dimR && j < dimC) {
                    m->mat[i * DIM_MAT + j] = mat[i * DIM_MAT + j];
                } else {
                    m->mat[i * DIM_MAT + j] = 0.0;
                }
            }
        }
    }
    Mtx_update_active(m);
}

uint8_t Mtx_is_sparse(Mtx *m) {
    // Consider the matrix sparse if less than 25% of elements are non-zero
    size_t total_active = m->active_num_row * m->active_num_col;
    size_t threshold = (DIM_MAT * DIM_MAT) / 4;
    return total_active < threshold;
}

void Mtx_copy(Mtx *m1, Mtx *m2) {
    // copy m1 on m2
    if (m1->dimR != m2->dimR || m1->dimC != m2->dimC) {
        printf("Error: wrong dim for Mtx_copy : m1 dim(%zux%zu) but m2 dim(%zux%zu)\n",
               m1->dimR, m1->dimC, m2->dimR, m2->dimC);
        return;
    }
    memcpy(m2->mat, m1->mat, DIM_MAT * DIM_MAT * sizeof(float));
}

void Mtx_print(Mtx *m) {
    printf("[");
    for (size_t i = 0; i < m->dimR; i++) {
        printf("[");
        for (size_t j = 0; j < m->dimC; j++) {
            printf("%f, ", m->mat[i * DIM_MAT + j]);
        }
        printf("]\n ");
    }
    printf("]\n");
}

Mtx *Mtx_copy_list(Mtx *m, float *mat) {
    for (size_t i = 0; i < m->dimR; i++) {
        for (size_t j = 0; j < m->dimC; j++) {
            m->mat[i * DIM_MAT + j] = mat[i * m->dimC + j];
        }
    }
    return m;
}

/*
1. DONE: Basic Matrix Operations
    - Addition: C[i][j] = A[i][j] + B[i][j] // Element-wise addition of two matrices
    - Subtraction: C[i][j] = A[i][j] - B[i][j] // Element-wise subtraction of two matrices
    - Scalar Multiplication: C[i][j] = scalar * A[i][j] // Multiply each element by a
scalar
    - Matrix Multiplication: C[i][k] += A[i][j] * B[j][k] // Standard matrix
multiplication
    - Transpose: C[j][i] = A[i][j] // Swap rows and columns
    - Element-wise Multiplication: C[i][j] = A[i][j] * B[i][j] // Multiply corresponding
elements
    - Element-wise Division: C[i][j] = A[i][j] / B[i][j] // Divide corresponding elements
*/

Mtx *Mtx_add(Mtx *m1, Mtx *m2, Mtx *mr) {
    if (m1->dimR != m2->dimR || m1->dimC != m2->dimC) {
        printf("Error: wrong dim for Mtx_add : m1 dim(%zux%zu) but m2 dim(%zux%zu)\n",
               m1->dimR, m1->dimC, m2->dimR, m2->dimC);
        return NULL;
    }
    if (mr == NULL) { // Modify m1 as the result
        size_t index;
        for (size_t i = 0; i < m1->dimR; i++) {
            if (!m2->active_row[i])
                continue; // Skip inactive rows in m2
            for (size_t j = 0; j < m1->dimC; j++) {
                index = i * DIM_MAT + j;
                m1->mat[index] += m2->mat[index];
            }
        }
        return m1;

    } else { // Store result in mr
        size_t index;
        for (size_t i = 0; i < m1->dimR; i++) {
            for (size_t j = 0; j < m1->dimC; j++) {
                index = i * DIM_MAT + j;
                mr->mat[index] = m1->mat[index] + m2->mat[index];
            }
        }
        return mr;
    }
}

Mtx *Mtx_sub(Mtx *m1, Mtx *m2, Mtx *mr) {
    if (m1->dimR != m2->dimR || m1->dimC != m2->dimC) {
        printf("Error: wrong dim for Mtx_sub : m1 dim(%zux%zu) but m2 dim(%zux%zu)\n",
               m1->dimR, m1->dimC, m2->dimR, m2->dimC);
        return NULL;
    }
    if (mr == NULL) { // Modify m1 as the result
        size_t index;
        for (size_t i = 0; i < m1->dimR; i++) {
            if (!m2->active_row[i])
                continue; // Skip inactive rows in m2
            for (size_t j = 0; j < m1->dimC; j++) {
                index = i * DIM_MAT + j;
                m1->mat[index] -= m2->mat[index];
            }
        }
        return m1;
    } else { // Store result in mr
        size_t index;
        for (size_t i = 0; i < m1->dimR; i++) {
            for (size_t j = 0; j < m1->dimC; j++) {
                index = i * DIM_MAT + j;
                mr->mat[index] = m1->mat[index] - m2->mat[index];
            }
        }
        return mr;
    }
}

Mtx *Mtx_scal(Mtx *m1, float scalar, Mtx *mr) {
    if (mr == NULL) { // Modify m1 as the result
        for (size_t i = 0; i < m1->dimR; i++) {
            for (size_t j = 0; j < m1->dimC; j++) {
                m1->mat[i * DIM_MAT + j] *= scalar;
            }
        }
        return m1;
    } else { // Store result in mr
        for (size_t i = 0; i < m1->dimR; i++) {
            for (size_t j = 0; j < m1->dimC; j++) {
                mr->mat[i * DIM_MAT + j] = m1->mat[i * DIM_MAT + j] * scalar;
            }
        }
        return mr;
    }
}

Mtx *Mtx_mult(Mtx *m1, Mtx *m2, Mtx *mr) {
    // Check that the dimensions correspond for matrix multiplication
    if (m1->dimC != m2->dimR) {
        printf("Error: wrong dim for Mtx_mult: m1->dimY = %zu but m2->dimX = %zu\n",
               m1->dimC, m2->dimR);
        return NULL;
    }
    if (mr == NULL || mr == m2) { // Modify `m1` as the result
        // Temporary buffer to store a column of the result
        float temp_mat[DIM_MAT * DIM_MAT] = {0};

        // Matrix multiplication logic
        for (size_t i = 0; i < m1->dimR; i++) {     // Iterate over rows of `m1`
            for (size_t j = 0; j < m2->dimC; j++) { // Iterate over columns of `m2`
                temp_mat[i * DIM_MAT + j] = 0;
                for (size_t k = 0; k < m1->dimC;
                     k++) { // Iterate over the shared dimension
                    temp_mat[i * DIM_MAT + j] +=
                        m1->mat[i * DIM_MAT + k] * m2->mat[k * DIM_MAT + j];
                }
            }
        }

        if (mr == NULL) {
            // Copy the result from the temporary buffer back into `m1`
            for (size_t i = 0; i < m1->dimR; i++) {
                for (size_t j = 0; j < m2->dimC; j++) {
                    m1->mat[i * DIM_MAT + j] = temp_mat[i * DIM_MAT + j];
                }
            }
            // Update dimensions of `m1`
            m1->dimC = m2->dimC;
            // Update active rows and columns
            // don't use the dedicate function beacause
            // this is faster here
            for (size_t i = 0; i < m1->dimR; i++) {
                m1->active_row[i] = 1; // All rows are now active
            }
            for (size_t j = 0; j < m1->dimC; j++) {
                m1->active_col[j] = 1; // All columns are now active
            }
            m1->active_num_row = m1->dimR;
            m1->active_num_col = m1->dimC;
            return m1;
        } else {
            // Copy the result from the temporary buffer back into `m1`
            for (size_t i = 0; i < m1->dimR; i++) {
                for (size_t j = 0; j < m2->dimC; j++) {
                    m2->mat[i * DIM_MAT + j] = temp_mat[i * DIM_MAT + j];
                }
            }
            // Update dimensions of `m1`
            m2->dimR = m2->dimR;
            // Update active rows and columns
            // don't use the dedicate function beacause
            // this is faster here
            for (size_t i = 0; i < m2->dimR; i++) {
                m2->active_row[i] = 1; // All rows are now active
            }
            for (size_t j = 0; j < m2->dimC; j++) {
                m2->active_col[j] = 1; // All columns are now active
            }
            m2->active_num_row = m2->dimR;
            m2->active_num_col = m2->dimC;
            return m2;
        }
    } else { // Modify `mr`
        // Ensure `mr` has correct dimensions
        mr->dimR = m1->dimR;
        mr->dimC = m2->dimC;
        // Perform matrix multiplication
        size_t index;
        for (size_t i = 0; i < m1->dimR; i++) {
            for (size_t j = 0; j < m2->dimC; j++) {
                mr->mat[i * DIM_MAT + j] = 0; // Initialize to zero
                index = i * DIM_MAT + j;
                for (size_t k = 0; k < m1->dimC; k++) {
                    mr->mat[index] += m1->mat[i * DIM_MAT + k] * m2->mat[k * DIM_MAT + j];
                }
            }
        }
        // Update active rows and columns in `mr`
        // don't use the dedicate function beacause
        // this is faster here
        for (size_t i = 0; i < mr->dimR; i++) {
            mr->active_row[i] = 1;
        }
        for (size_t j = 0; j < mr->dimC; j++) {
            mr->active_col[j] = 1;
        }
        mr->active_num_row = mr->dimR;
        mr->active_num_col = mr->dimC;

        return mr;
    }
}

Mtx *Mtx_T(Mtx *m1, Mtx *mr) {
    if (mr == NULL) { // In-place transpose
        for (size_t i = 0; i < m1->dimC; i++) {
            for (size_t j = i + 1; j < m1->dimR; j++) {
                float temp = m1->mat[j * DIM_MAT + i];
                m1->mat[j * DIM_MAT + i] = m1->mat[i * DIM_MAT + j];
                m1->mat[i * DIM_MAT + j] = temp;
            }
        }
        // Swap dimensions
        size_t temp = m1->dimR;
        m1->dimR = m1->dimC;
        m1->dimC = temp;

        // Swap active rows and columns
        uint8_t temp_active[DIM_MAT];
        memcpy(temp_active, m1->active_row, DIM_MAT * sizeof(uint8_t));
        memcpy(m1->active_row, m1->active_col, DIM_MAT * sizeof(uint8_t));
        memcpy(m1->active_col, temp_active, DIM_MAT * sizeof(uint8_t));
        size_t temp_count = m1->active_num_row;
        m1->active_num_row = m1->active_num_col;
        m1->active_num_col = temp_count;

        return m1;

    } else { // Transpose into mr
        mr->dimR = m1->dimC;
        mr->dimC = m1->dimR;
        for (size_t i = 0; i < m1->dimC; i++) {
            for (size_t j = 0; j < m1->dimR; j++) {
                mr->mat[i * DIM_MAT + j] = m1->mat[j * DIM_MAT + i];
            }
        }
        // Update active rows and columns
        for (size_t i = 0; i < m1->dimC; i++) {
            mr->active_col[i] = m1->active_row[i];
        }
        for (size_t j = 0; j < m1->dimR; j++) {
            mr->active_row[j] = m1->active_col[j];
        }
        mr->active_num_row = m1->active_num_col;
        mr->active_num_col = m1->active_num_row;

        return mr;
    }
}

Mtx *Mtx_mult_w(Mtx *m1, Mtx *m2, Mtx *mr) {
    if (m1->dimR != m2->dimR || m1->dimC != m2->dimC) {
        printf("Error: wrong dim for Mtx_mult_w : m1 dim(%zux%zu) but m2 dim(%zux%zu)\n",
               m1->dimR, m1->dimC, m2->dimR, m2->dimC);
        return NULL;
    }
    if (mr == NULL) { // Modify m1 as the result
        size_t index;
        for (size_t i = 0; i < m1->dimR; i++) {
            if (!m2->active_row[i])
                continue; // Skip inactive rows in m2
            for (size_t j = 0; j < m1->dimC; j++) {
                index = i * DIM_MAT + j;
                m1->mat[index] *= m2->mat[index];
            }
        }
        return m1;
    } else { // Store result in mr
        size_t index;
        for (size_t i = 0; i < DIM_MAT; i++) {
            for (size_t j = 0; j < DIM_MAT; j++) {
                index = i * DIM_MAT + j;
                mr->mat[index] = m1->mat[index] * m2->mat[index];
            }
        }
        return mr;
    }
}

Mtx *Mtx_div_w(Mtx *m1, Mtx *m2, Mtx *mr) {
    if (m1->dimR != m2->dimR || m1->dimC != m2->dimC) {
        printf("Error: wrong dim for Mtx_div_w : m1 dim(%zux%zu) but m2 dim(%zux%zu)\n",
               m1->dimR, m1->dimC, m2->dimR, m2->dimC);
        return NULL;
    }
    if (mr == NULL) { // Modify m1 as the result
        size_t index;
        for (size_t i = 0; i < m1->dimR; i++) {
            if (!m2->active_row[i])
                continue; // Skip inactive rows in m2
            for (size_t j = 0; j < m1->dimC; j++) {
                if (!m2->active_col[j])
                    continue;            // Skip inactive columns in m2
                index = i * DIM_MAT + j; // 1D index calculation
                if (m2->mat[index] != 0.0) {
                    m1->mat[index] /= m2->mat[index];
                } else {
                    m1->mat[index] = INFINITY; // Division by zero results in infinity
                }
            }
        }
        return m1;
    } else { // Store result in mr
        size_t index;
        for (size_t i = 0; i < m1->dimR; i++) {
            for (size_t j = 0; j < m1->dimC; j++) {
                index = i * DIM_MAT + j; // 1D index calculation
                if (m2->mat[index] != 0.0) {
                    mr->mat[index] = m1->mat[index] / m2->mat[index];
                } else {
                    mr->mat[index] = INFINITY; // Division by zero results in infinity
                }
            }
        }
        return mr;
    }
}

Mtx *Mtx_eye(Mtx *m) {
    // add one to each value of the main diagonal
    if (m->dimR != m->dimC) {
        printf("Error: wrong dim for Mtx_eye m must be square but m->dimR = %zu != %zu = "
               "m->dimC\n",
               m->dimR, m->dimC);
        return NULL;
    }
    for (size_t i = 0; i < m->dimR; i++) {
        m->mat[i * DIM_MAT + i] += 1;
    }
    return m;
}

Mtx *Mtx_shift(Mtx *m, float scalar) {
    // add one to each value of the main diagonal
    if (m->dimR != m->dimC) {
        printf(
            "Error: wrong dim for Mtx_shift m must be square but m->dimR = %zu != %zu = "
            "m->dimC\n",
            m->dimR, m->dimC);
        return NULL;
    }
    for (size_t i = 0; i < m->dimR; i++) {
        m->mat[i * DIM_MAT + i] += scalar;
    }
    return m;
}

/*
2. TODO: Decompositions
    ok - LU Decomposition: A = L * U // Factor into lower and upper triangular matrices
    ok - QR Decomposition: A = Q * R // Factor into orthogonal (Q) and triangular (R)
    ok - Cholesky Decomposition: A = L * L^T // For symmetric positive-definite matrices
    - Eigen Decomposition: A = Q * Lambda * Q^-1 // Break into eigenvalues/vectors
    - Singular Value Decomposition (SVD): A = U * Sigma * V^T // General matrix
factorization
    - Schur Decomposition: A = Q * T * Q^T // Similar to triangular factorization
    - Jordan Decomposition: A = P * J * P^-1 // Converts to Jordan normal form
*/

void Mtx_eigen_value(Mtx *m, Vect *eigen_val) {}

void Mtx_cholesky(Mtx *m, Mtx *ml) {
    if (m->dimR != m->dimC) {
        printf("Error: wrong for Mtx_cholesky m must be square but m->dimR = %zu != %zu "
               "= m->dimC\n",
               m->dimR, m->dimC);
    }
    // check if define positive and symmetric
    // this must be done before since it's to complex in term of time

    // cholesky
    for (size_t i = 0; i < m->dimR; i++) {
        float lii;
        for (size_t j = 0; j < m->dimC; j++) {
            if (i == j) {
                lii = m->mat[i * DIM_MAT + j];
                float liik = 0;
                for (size_t k = 0; k < j; k++) {
                    liik += ml->mat[j * DIM_MAT + k] * ml->mat[j * DIM_MAT + k];
                }
                lii = sqrtf(lii - liik);
                ml->mat[i * DIM_MAT + j] = lii;
            } else if (j < i) {
                float lijk = 0;
                for (size_t k = 0; k < j; k++) {
                    lijk += ml->mat[i * DIM_MAT + k] * ml->mat[j * DIM_MAT + k];
                }
                ml->mat[i * DIM_MAT + j] = (1.0 / lii) * (m->mat[i * DIM_MAT + j] - lijk);
            }
        }
    }
}

void Mtx_cholesky_d(Mtx *m, Mtx *ml, Mtx *md) {
    if (m->dimR != m->dimC) {
        printf("Error: wrong for Mtx_cholesky_d m must be square but m->dimR = %zu != "
               "%zu = m->dimC\n",
               m->dimR, m->dimC);
    }
    // check if symmetric

    // cholesky with diag
}

void Mtx_qr(Mtx *m, Mtx *q, Mtx *r) { // good for dimR == dimC not sure for the rest
    if (m->dimR < m->dimC) {
        printf("Error: Matrix m must have rows >= columns for QR decomposition.\nbut "
               "here dimR = %zu and dimC = %zu\n",
               m->dimR, m->dimC);
        return;
    }
    if (q == NULL || r == NULL) {
        printf("the matrix q and/or r doesn't exist\n");
        return;
    }
    // by householder reflection
    Mtx_init(q, m->dimR, m->dimR, NULL);
    Mtx_eye(q);
    Mtx_init(r, m->dimR, m->dimC, m->mat);
    Mtx qi, h;
    Mtx_init(&h, m->dimR, m->dimR, NULL);
    Mtx id;

    // find Q
    for (size_t i = 0; i < m->dimR; i++) {
        Vect vect;
        Vect_init(&vect, m->dimR - i, NULL);
        for (size_t j = i; j < m->dimR; j++) {
            vect.v[j - i] = r->mat[j * DIM_MAT + i];
        }
        // qi  = I - (2 * vv^T)
        Mtx_init(&qi, m->dimR - i, m->dimR - i, NULL);
        Mtx_init(&id, m->dimR - i, m->dimR - i, NULL);
        Mtx_eye(&id);
        // alpha
        float alpha = Vect_norm(&vect, true);
        if (vect.v[0] > 0) {
            alpha = -alpha;
        }
        // printf("############################################\n");
        // printf("        %zu\n",i);
        // printf("############################################\n");
        vect.v[0] -= alpha;
        Vect_scal(&vect, 1.0 / Vect_norm(&vect, true), NULL);
        // printf("alpha : %f \n",alpha);
        // printf("vect : ");
        // Vect_print(&vect);
        Vect_mult_to_Mtx(&vect, &vect, &qi);
        // printf("q%zu prime: \n",i);
        // Mtx_print(&qi);
        Mtx_scal(&qi, 2.0, NULL);
        Mtx_sub(&id, &qi, &qi);
        // printf("q%zu prime: \n",i);
        // Mtx_print(&qi);
        //  Q: (m-i) x (m-i) ->  m x m
        Mtx_init(&h, m->dimR, m->dimC, NULL);
        for (size_t j = 0; j < m->dimR; j++) {
            for (size_t k = 0; k < m->dimC; k++) {
                if (j < i) {
                    if (k == j) {
                        h.mat[j * DIM_MAT + k] = 1.0;
                    } else {
                        h.mat[j * DIM_MAT + k] = 0.0;
                    }
                } else {
                    h.mat[j * DIM_MAT + k] = qi.mat[(j - i) * DIM_MAT + (k - i)];
                }
            }
        }
        // printf(" q%zu : \n",i);
        // Mtx_print(&h);
        Mtx_mult(&h, r, r);
        // printf(" r%zu : \n",i);
        // Mtx_print(r);
        Mtx_mult(&h, q, q);
    }
    Mtx_T(q, NULL);
}

void Mtx_lu_pp(Mtx *m, Mtx *mp, Mtx *ml, Mtx *mu, size_t *rowswap) {
    // with partial pivoting
    if (m->dimR != m->dimC) {
        printf("Error: wrong dim for Mtx_lu_pp m must be square but m->dimR = %zu != %zu "
               "= m->dimC",
               m->dimR, m->dimC);
        return;

    } else if (ml != NULL && mu != NULL && mp != NULL) {
        // return the ml and mu in there correct matrix
        Mtx_init(mp, m->dimR, m->dimC, NULL);
        Mtx_init(ml, m->dimR, m->dimC, NULL);
        Mtx_init(mu, m->dimR, m->dimC, NULL);
        Mtx_eye(mp);
        Mtx_copy(m, mu);

        for (size_t i = 0; i < m->dimR; i++) {
            // find the index of max_pivot_col
            float max_pivot_col = fabsf(mu->mat[i * DIM_MAT + i]);
            size_t index_max_col = i;
            for (size_t k = i + 1; k < m->dimR; k++) {
                float temp = fabsf(mu->mat[k * DIM_MAT + i]);
                if (temp > max_pivot_col) {
                    max_pivot_col = temp;
                    index_max_col = k;
                }
            }
            // swap the line
            if (index_max_col != i) {
                if (rowswap != NULL) {
                    *rowswap += 1;
                }
                float temp;
                for (size_t j = 0; j < m->dimC; j++) {
                    // swap the row in mu
                    temp = mu->mat[index_max_col * DIM_MAT + j];
                    mu->mat[index_max_col * DIM_MAT + j] = mu->mat[i * DIM_MAT + j];
                    mu->mat[i * DIM_MAT + j] = temp;
                    // swap the row in mp
                    temp = mp->mat[index_max_col * DIM_MAT + j];
                    mp->mat[index_max_col * DIM_MAT + j] = mp->mat[i * DIM_MAT + j];
                    mp->mat[i * DIM_MAT + j] = temp;
                    // swap the row in ml up to the current column
                    if (j < i) {
                        temp = ml->mat[index_max_col * DIM_MAT + j];
                        ml->mat[index_max_col * DIM_MAT + j] = ml->mat[i * DIM_MAT + j];
                        ml->mat[i * DIM_MAT + j] = temp;
                    }
                }
            }
            // elimination element
            for (size_t j = i + 1; j < m->dimR; j++) {
                float multiplier = mu->mat[j * DIM_MAT + i] / mu->mat[i * DIM_MAT + i];
                ml->mat[j * DIM_MAT + i] = multiplier;
                for (size_t k = 0; k < m->dimC; k++) {
                    mu->mat[j * DIM_MAT + k] -= multiplier * mu->mat[i * DIM_MAT + k];
                }
            }
        }
        // fill the diagonal of ml with 1
        Mtx_eye(ml);
        return;
    } else {
        printf("Error: wrong there is an issue with the matrix mp,ml,mu\n");
        return;
    }
}

/*
3. TODO: Intermediate Operations
    - Determinant: Compute det(A) // Scalar value for invertibility
    - Inverse: A^-1 such that A * A^-1 = I // Requires a square matrix
    - Trace: trace = sum(A[i][i]) // Sum of diagonal elements
    - Rank: Maximum number of linearly independent rows or columns
    WARNING:LATER
    - Outer Product: C[i][j] = A[i] * B[j] // Builds a matrix from two vectors
    - Norms:
      * Frobenius norm: sqrt(sum(A[i][j]^2)) // Measure of matrix size
      * p-norms: (sum(abs(A[i][j])^p))^(1/p) // Generalized size measure
*/

Mtx *Mtx_inv(Mtx *m, Mtx *mr) {
    // WARNING: doesn't really work to much error!
    Mtx mp;
    Mtx ml;
    Mtx mu;
    Mtx temp;
    Mtx_init(&temp, m->dimR, m->dimC, NULL);
    Mtx_lu_pp(m, &mp, &ml, &mu, NULL);

    // invers ml
    for (size_t i = 0; i < m->dimR; i++) {
        for (size_t j = 0; j < m->dimC; j++) {
            if (j < i) {
                ml.mat[i * DIM_MAT + j] *= -1;
            }
        }
    }
    // inverse mu
    Mtx_eye(&temp);
    for (size_t col = 0; col < m->dimR; col++) {
        for (size_t i = m->dimR; i-- > 0;) {
            float sum = 0.0f;
            for (size_t j = i + 1; j < m->dimR; j++) {
                sum += mu.mat[i * DIM_MAT + j] * temp.mat[j * DIM_MAT + col];
            }
            temp.mat[i * DIM_MAT + col] =
                (temp.mat[i * DIM_MAT + col] - sum) / mu.mat[i * DIM_MAT + i];
        }
    }
    Mtx_copy(&temp, &mu);
    // inverse m : pm = lu -> m^-1 p^-1 = u^-1 l^-1
    //         m^-1 = u^-1 l^-1 p
    if (mr == NULL) {
        Mtx_mult(Mtx_mult(&mu, &ml, m), &mp, NULL);
        return m;
    } else {
        Mtx_mult(Mtx_mult(&mu, &ml, mr), &mp, NULL);
        return mr;
    }
}

size_t Mtx_rank(Mtx *m) {
    if (m->dimR == m->dimC) {
        Mtx mp_stack;
        Mtx *mp = &mp_stack;
        Mtx ml_stack;
        Mtx *ml = &ml_stack;
        Mtx mu_stack;
        Mtx *mu = &mu_stack;
        size_t rank = m->dimR;
        Mtx_lu_pp(m, mp, ml, mu, NULL);
        for (size_t i = m->dimR - 1; i >= 0; i--) {
            if (mu->mat[i * DIM_MAT + m->dimR - 1] == 0.0) {
                rank--;
            } else {
                break;
            }
        }
        return rank;

    } else {
        // TODO:
        // call qr methode
        return m->dimR;
    }
}

float Mtx_det(Mtx *m) {
    if (m->dimC != m->dimR) {
        printf("Error: wrong dim for Mtx_det : m->dimC = %zu but m2->dimR = %zu\n",
               m->dimC, m->dimR);
        return 0.0f;
    }
    float det = 1;
    Mtx mp;
    Mtx ml;
    Mtx mu;
    size_t rowswap = 0;

    Mtx_lu_pp(m, &mp, &ml, &mu, &rowswap);

    for (size_t i = 0; i < m->dimR; i++) {
        det *= mu.mat[i * DIM_MAT + i];
    }
    if (rowswap % 2 != 0) {
        det = -det;
    }
    printf("%zu = rowswap\n", rowswap);

    return det;
}

float Mtx_trace(Mtx *m) {
    if (m->dimC != m->dimR) {
        printf("Error: wrong dim for Mtx_trace : m->dimC = %zu but m2->dimR = %zu\n",
               m->dimC, m->dimR);
        return 0.0f;
    }
    float trace = 0.0;
    for (size_t i = 0; i < m->dimR; i++) {
        trace += m->mat[i * DIM_MAT + i];
    }
    return trace;
}

/*
4. TODO: Advanced Matrix Operations
    - Matrix Exponential: Compute e^A = I + A + A^2/2! + A^3/3! ... // Solve differential
systems
    - Matrix Logarithm: Compute ln(A) // Inverse of matrix exponential
    - Matrix Power: Compute A^n // Raise matrix to integer or fractional powers
    - Diagonalization: A = P * D * P^-1 // Convert to diagonal form
    - Pseudo-Inverse: Compute A^+ // General inverse for non-square matrices
    - Kronecker Product: Builds larger block matrix from A and B
    - Block Matrix Operations: Perform arithmetic at the block level
*/

/*
5. TODO: Special Matrix Computations
    - Solving Linear Systems: Solve Ax = b // Use decomposition or iterative solvers
    - Least Squares Solution: Solve min(||Ax - b||^2) // Best fit for over-determined
systems
    - Condition Number: kappa(A) = ||A|| * ||A^-1|| // Measures sensitivity to changes
    - Projection:
      * Orthogonal projection: Projects vector onto a subspace
      * Projection matrix: P = A * (A^T * A)^-1 * A^T
    - Null Space: Find x such that Ax = 0 // Basis of null space
    - Column Space: Span of matrix columns
    - Cofactor Expansion: Use minors and cofactors for determinant calculation
    - Adjugate Matrix: Transpose of the cofactor matrix
*/

/*
6. TODO: Applications of Matrices
    - Markov Chains: Use transition matrices for probability modeling
    - Principal Component Analysis (PCA): Reduce dimensions using eigen decomposition
    - Graph Theory:
      * Adjacency matrices: Represent connections in a graph
      * Laplacian matrices: Analyze graph properties
    - Fourier Transforms: Matrix representation of the discrete Fourier transform
    - Image Processing:
      * Convolutions: Use matrix filters for images
      * Transformations: Scaling, rotation, and translation
    - Control Theory:
      * State transition matrices: Model dynamic systems
      * Lyapunov equations: Stability analysis
    - Quantum Computing:
      * Unitary matrices: Represent quantum gate operations
      * Tensor products: Model quantum entanglement
*/

/*
7. TODO: Specialized Matrix Types
    - Symmetric Matrix: A[i][j] == A[j][i]
    - Hermitian Matrix: Complex symmetric matrix with real eigenvalues
    - Positive Definite Matrix: x^T * A * x > 0 for all x != 0
    - Sparse Matrix: Most elements are zero
    - Diagonal Matrix: Non-zero elements only on the diagonal
    - Orthogonal Matrix: Q * Q^T = I
    - Permutation Matrix: Rearranges rows or columns
    - Toeplitz Matrix: Constant diagonals for efficient computation
*/

/*
8. TODO: Numerical Analysis Techniques
    - Iterative Solvers:
      * Conjugate Gradient: Solve sparse systems iteratively
      * GMRES: Solve non-symmetric sparse systems
    - Conditioning and Stability: Assess numerical robustness
    - Deflation: Reduce problem size by factoring out known solutions
*/
