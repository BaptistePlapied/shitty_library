#include "../../headers/general.h"
#include "../../headers/sequential.h"
#include "../../headers/utils.h"
#include <CUnit/Basic.h>
#include <CUnit/CUnit.h>
#include <stdint.h>
#include <stdio.h>

void test_forward_sub() {
    uint64_t dim_m = rand_dim_gen(5);
    matrix *A = m_init(dim_m, dim_m);
    rand_m_data_gen(A, 10.0);
    matrix *L = m_init(dim_m, dim_m);
    matrix *U = m_init(dim_m, dim_m);
    matrix *P = m_init(dim_m, dim_m);
    lu_pp(A, L, U, P, NULL);
    matrix *x_m = m_init(dim_m, 1); // "vector" as a matrix
    rand_m_data_gen(x_m, 10.0);
    matrix *b_m = m_init(dim_m, 1); // "vector" as a matrix
    m_mult(L, x_m, b_m);
    vector *b = v_init_l(dim_m, b_m->data);
    vector *result = forward_sub(L, b, NULL);
    /* printf("L : \n"); */
    /* m_printf(L); */
    /* printf("x : \n"); */
    /* m_printf(x_m); */
    /* printf("b : \n"); */
    /* m_printf(b_m); */
    /* printf("result : \n"); */
    /* v_printf(result); */
    for (uint64_t i = 0; i < dim_m; i++) {
        CU_ASSERT_DOUBLE_EQUAL(result->data[i].Re, x_m->data[i].Re, 1e-6);
        CU_ASSERT_DOUBLE_EQUAL(result->data[i].Im, x_m->data[i].Im, 1e-6);
    }
    v_free(result);
    m_free(A);
    m_free(L);
    m_free(U);
    m_free(P);
    m_free(b_m);
    m_free(x_m);
    v_free(b);
}

void test_back_sub() {
    uint64_t dim_m = rand_dim_gen(5);
    matrix *A = m_init(dim_m, dim_m);
    rand_m_data_gen(A, 10.0);
    matrix *L = m_init(dim_m, dim_m);
    matrix *U = m_init(dim_m, dim_m);
    matrix *P = m_init(dim_m, dim_m);
    lu_pp(A, L, U, P, NULL);
    matrix *x_m = m_init(dim_m, 1); // "vector" as a matrix
    rand_m_data_gen(x_m, 10.0);
    matrix *b_m = m_init(dim_m, 1); // "vector" as a matrix
    m_mult(U, x_m, b_m);
    vector *b = v_init_l(dim_m, b_m->data);
    vector *result = back_sub(U, b, NULL);
    /* printf("U : \n"); */
    /* m_printf(U); */
    /* printf("x : \n"); */
    /* m_printf(x_m); */
    /* printf("b : \n"); */
    /* m_printf(b_m); */
    /* printf("result : \n"); */
    /* v_printf(result); */
    for (uint64_t i = 0; i < dim_m; i++) {
        CU_ASSERT_DOUBLE_EQUAL(result->data[i].Re, x_m->data[i].Re, 1e-6);
        CU_ASSERT_DOUBLE_EQUAL(result->data[i].Im, x_m->data[i].Im, 1e-6);
    }
    v_free(result);
    m_free(A);
    m_free(L);
    m_free(U);
    m_free(P);
    m_free(b_m);
    m_free(x_m);
    v_free(b);
}
// TODO: Fix this test
// Doesn't work yet (i guess)
void test_m_echelon() {
    uint64_t dim_m = rand_dim_gen(15);
    matrix *A = m_init(dim_m, dim_m);
    rand_m_data_gen(A, 5);
    echelon_expr echelon;
    matrix *c = m_echelon(A, &echelon, NULL);
    echelon_expr_printf(echelon);
    vector *b = v_init(dim_m);
    rand_v_data_gen(b, 5);
    vector *d = v_apply_echelon(b, &echelon, NULL);
    matrix *b_m = m_init_l(b->m, 1, b->data);
    matrix *A_inv = m_inv_lu(A, NULL);
    matrix *x_m = m_mult(A_inv, b_m, NULL);
    /* printf("A : \n"); */
    /* m_printf(A); */
    /* printf("A echelon = C : \n"); */
    /* m_printf(c); */
    /* printf("backsub :\n\t"); */
    /* printf(" b : \n\t\t"); */
    /* v_printf(b); */
    /* printf(" result : \n\t\t"); */
    /* v_printf(d); */
    /* printf(" x: \n\t\t"); */
    /* m_printf(m_transpose(x_m, NULL)); */
    /* printf("I = A_inv * A: \n"); */
    /* m_printf(m_mult(A_inv, A, NULL)); */
    v_free(d);
    v_free(b);
    m_free(c);
    m_free(A);
    m_free(A_inv);
    m_free(b_m);
    m_free(x_m);
}

void test_lu_pp() {
    uint64_t dim_m = rand_dim_gen(10);
    matrix *A = m_init(dim_m, dim_m);
    rand_m_data_gen(A, 100);
    /* printf("A : \n"); */
    /* m_printf(A); */
    matrix *U = m_init(dim_m, dim_m);
    matrix *L = m_init(dim_m, dim_m);
    matrix *P = m_init(dim_m, dim_m);
    lu_pp(A, L, U, P, NULL);
    /* printf("L : \n"); */
    /* m_printf(L); */
    /* printf("U : \n"); */
    /* m_printf(U); */
    /* printf("P : \n"); */
    /* m_printf(P); */
    matrix *PA = m_init(dim_m, dim_m);
    matrix *LU = m_init(dim_m, dim_m);
    m_mult(P, A, PA);
    m_mult(L, U, LU);
    /* m_printf(PA); */
    /* printf("LU : \n"); */
    /* m_printf(LU); */
    // Check if PA == LU
    for (uint64_t i = 0; i < dim_m; i++) {
        for (uint64_t j = 0; j < dim_m; j++) {
            CU_ASSERT_DOUBLE_EQUAL(PA->data[i * dim_m + j].Re, LU->data[i * dim_m + j].Re,
                                   1e-6);
            CU_ASSERT_DOUBLE_EQUAL(PA->data[i * dim_m + j].Im, LU->data[i * dim_m + j].Im,
                                   1e-6);
        }
    }
    m_free(A);
    m_free(L);
    m_free(U);
    m_free(P);
    m_free(PA);
    m_free(LU);
}

void test_m_inv_lu() {
    uint64_t dim_m = rand_dim_gen(10);
    matrix *A = m_init(dim_m, dim_m);
    rand_m_data_gen(A, 100);
    /* printf("A : \n"); */
    /* m_printf(A); */
    matrix *Result = m_inv_lu(A, NULL);
    /* printf("Result : \n"); */
    /* m_printf(Result); */
    // Check if Result is the inverse of A
    matrix *I = m_mult(A, Result, NULL);
    /* printf("I = A * Result : \n"); */
    /* m_printf(I); */
    for (uint64_t i = 0; i < dim_m; i++) {
        for (uint64_t j = 0; j < dim_m; j++) {
            if (i == j) {
                CU_ASSERT_DOUBLE_EQUAL(I->data[i * dim_m + j].Re, 1.0, 1e-6);
                CU_ASSERT_DOUBLE_EQUAL(I->data[i * dim_m + j].Im, 0.0, 1e-6);
            } else {
                CU_ASSERT_DOUBLE_EQUAL(I->data[i * dim_m + j].Re, 0.0, 1e-6);
                CU_ASSERT_DOUBLE_EQUAL(I->data[i * dim_m + j].Im, 0.0, 1e-6);
            }
        }
    }
    m_free(A);
    m_free(Result);
    m_free(I);
}

void test_m_det() {
    uint64_t dim_m = rand_dim_gen(10);
    matrix *A = m_init(dim_m, dim_m);
    rand_m_data_gen(A, 5);
    /* printf("A : \n"); */
    /* m_printf(A); */
    complex det = m_det(A);
    /* printf("det : %f + %fi\n", det.Re, det.Im); */
    // Check if det is correct
    matrix *L = m_init(dim_m, dim_m);
    matrix *U = m_init(dim_m, dim_m);
    matrix *P = m_init(dim_m, dim_m);
    uint64_t swap;
    lu_pp(A, L, U, P, &swap);
    complex tr_det = m_tr_det(U);
    if (swap % 2 == 1) {
        tr_det = c_mult((complex){-1.0, 0.0}, tr_det);
    }
    CU_ASSERT_DOUBLE_EQUAL(det.Re, tr_det.Re, 1e-6);
    CU_ASSERT_DOUBLE_EQUAL(det.Im, tr_det.Im, 1e-6);
    m_free(A);
    m_free(L);
    m_free(U);
    m_free(P);
}

void test_m_mgs() {
    uint64_t gid = find_empty_group_id();
    uint64_t dim_m = rand_dim_gen(10);
    matrix *A = m_init(dim_m, dim_m);
    rand_m_data_gen(A, 10.0);
    matrix *Q = m_init(dim_m, dim_m);
    matrix *R = m_init(dim_m, dim_m);
    m_mgs(A, Q, R);
    /* printf("A : \n"); */
    /* m_printf(A); */
    /* printf("Q : \n"); */
    /* m_printf(Q); */
    /* printf("R : \n"); */
    /* m_printf(R); */
    /* printf("I = Q^H * Q : \n"); */
    /* m_printf(m_track_in(m_mult(m_hermitian(Q, NULL), Q, NULL), gid)); */
    /* printf("Q * R : \n"); */
    /* m_printf(m_track_in(m_mult(Q, R, NULL), gid)); */
    // Check if Q is orthogonal and R is upper triangular
    matrix *QHQ = m_track_get_g(gid, 0);
    for (uint64_t i = 0; i < dim_m; i++) {
        for (uint64_t j = 0; j < dim_m; j++) {
            if (i == j) {
                CU_ASSERT_DOUBLE_EQUAL(QHQ->data[i * dim_m + j].Re, 1.0, 1e-6);
                CU_ASSERT_DOUBLE_EQUAL(QHQ->data[i * dim_m + j].Im, 0.0, 1e-6);
            } else {
                CU_ASSERT_DOUBLE_EQUAL(QHQ->data[i * dim_m + j].Re, 0.0, 1e-6);
                CU_ASSERT_DOUBLE_EQUAL(QHQ->data[i * dim_m + j].Im, 0.0, 1e-6);
            }
        }
    }
    for (uint64_t i = 0; i < dim_m; i++) {
        for (uint64_t j = 0; j < dim_m; j++) {
            if (i < j) {
                CU_ASSERT_DOUBLE_EQUAL(R->data[j * dim_m + i].Re, 0.0, 1e-6);
                CU_ASSERT_DOUBLE_EQUAL(R->data[j * dim_m + i].Im, 0.0, 1e-6);
            }
        }
    }
    // Check if Q * R = A
    matrix *QR = m_track_get_g(gid, 1);
    for (uint64_t i = 0; i < dim_m; i++) {
        for (uint64_t j = 0; j < dim_m; j++) {
            CU_ASSERT_DOUBLE_EQUAL(QR->data[i * dim_m + j].Re, A->data[i * dim_m + j].Re,
                                   1e-6);
            CU_ASSERT_DOUBLE_EQUAL(QR->data[i * dim_m + j].Im, A->data[i * dim_m + j].Im,
                                   1e-6);
        }
    }
    m_free(A);
    m_free(Q);
    m_free(R);
    /* printf("Tracked items in group %llu: %llu\n", gid, track_group_count(gid)); */
    track_group_clear(gid);
    /* printf("Tracked items in group %llu: %llu\n", gid, track_group_count(gid)); */
}

int main() {
    rand_init_seed();
    printf("Running Algo Tests...\n");
    CU_initialize_registry();
    CU_pSuite suite = CU_add_suite("Math Library Tests", NULL, NULL);

    CU_add_test(suite, "test_forward_sub", test_forward_sub);
    CU_add_test(suite, "test_back_sub", test_back_sub);
    /* CU_add_test(suite, "test_m_echelon", test_m_echelon); */
    CU_add_test(suite, "test_lu_pp", test_lu_pp);
    CU_add_test(suite, "test_m_inv_lu", test_m_inv_lu);
    CU_add_test(suite, "test_m_det", test_m_det);
    CU_add_test(suite, "test_m_mgs", test_m_mgs);

    CU_basic_run_tests();
    CU_cleanup_registry();
    return 0;
}
