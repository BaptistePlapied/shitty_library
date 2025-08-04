#include "../../headers/general.h"
#include "../../headers/sequential.h"
#include "../../headers/utils.h"
#include <CUnit/Basic.h>
#include <CUnit/CUnit.h>
#include <stdint.h>
#include <stdio.h>

void test_m_add() {
    uint64_t m = rand_dim_gen(100);
    uint64_t n = rand_dim_gen(100);
    matrix *A = m_init(m, n);
    matrix *B = m_init(m, n);
    matrix *result = m_init(m, n);
    rand_m_data_gen(A, 10.0);
    rand_m_data_gen(B, 10.0);
    m_add(A, B, result);
    for (uint64_t i = 0; i < m; i++) {
        for (uint64_t j = 0; j < n; j++) {
            CU_ASSERT_DOUBLE_EQUAL(result->data[i * n + j].Re,
                                   A->data[i * n + j].Re + B->data[i * n + j].Re, 1e-6);
            CU_ASSERT_DOUBLE_EQUAL(result->data[i * n + j].Im,
                                   A->data[i * n + j].Im + B->data[i * n + j].Im, 1e-6);
        }
    }
    m_free(A);
    m_free(B);
    m_free(result);
}
void test_m_sub() {
    uint64_t m = rand_dim_gen(100);
    uint64_t n = rand_dim_gen(100);
    matrix *A = m_init(m, n);
    matrix *B = m_init(m, n);
    matrix *result = m_init(m, n);
    rand_m_data_gen(A, 10.0);
    rand_m_data_gen(B, 10.0);
    m_sub(A, B, result);
    for (uint64_t i = 0; i < m; i++) {
        for (uint64_t j = 0; j < n; j++) {
            CU_ASSERT_DOUBLE_EQUAL(result->data[i * n + j].Re,
                                   A->data[i * n + j].Re - B->data[i * n + j].Re, 1e-6);
            CU_ASSERT_DOUBLE_EQUAL(result->data[i * n + j].Im,
                                   A->data[i * n + j].Im - B->data[i * n + j].Im, 1e-6);
        }
    }
    m_free(A);
    m_free(B);
    m_free(result);
}
void test_m_mult_e() {
    uint64_t m = rand_dim_gen(100);
    uint64_t n = rand_dim_gen(100);
    matrix *A = m_init(m, n);
    matrix *B = m_init(m, n);
    matrix *result = m_init(m, n);
    rand_m_data_gen(A, 10.0);
    rand_m_data_gen(B, 10.0);
    m_mult_e(A, B, result);
    for (uint64_t i = 0; i < m; i++) {
        for (uint64_t j = 0; j < n; j++) {
            CU_ASSERT_DOUBLE_EQUAL(result->data[i * n + j].Re,
                                   A->data[i * n + j].Re * B->data[i * n + j].Re -
                                       A->data[i * n + j].Im * B->data[i * n + j].Im,
                                   1e-6);
            CU_ASSERT_DOUBLE_EQUAL(result->data[i * n + j].Im,
                                   A->data[i * n + j].Re * B->data[i * n + j].Im +
                                       A->data[i * n + j].Im * B->data[i * n + j].Re,
                                   1e-6);
        }
    }
    m_free(A);
    m_free(B);
    m_free(result);
}
void test_m_div_e() {
    uint64_t m = rand_dim_gen(100);
    uint64_t n = rand_dim_gen(100);
    matrix *A = m_init(m, n);
    matrix *B = m_init(m, n);
    matrix *result = m_init(m, n);
    rand_m_data_gen(A, 10.0);
    rand_m_data_gen(B, 10.0);
    m_div_e(A, B, result);
    for (uint64_t i = 0; i < m; i++) {
        for (uint64_t j = 0; j < n; j++) {
            double denominator = B->data[i * n + j].Re * B->data[i * n + j].Re +
                                 B->data[i * n + j].Im * B->data[i * n + j].Im;
            CU_ASSERT_DOUBLE_EQUAL(result->data[i * n + j].Re,
                                   (A->data[i * n + j].Re * B->data[i * n + j].Re +
                                    A->data[i * n + j].Im * B->data[i * n + j].Im) /
                                       denominator,
                                   1e-6);
            CU_ASSERT_DOUBLE_EQUAL(result->data[i * n + j].Im,
                                   (A->data[i * n + j].Im * B->data[i * n + j].Re -
                                    A->data[i * n + j].Re * B->data[i * n + j].Im) /
                                       denominator,
                                   1e-6);
        }
    }
    m_free(A);
    m_free(B);
    m_free(result);
}
void test_m_mult() {
    uint64_t m = rand_dim_gen(100);
    uint64_t n = rand_dim_gen(100);
    uint64_t p = rand_dim_gen(100);
    matrix *A = m_init(m, n);
    matrix *B = m_init(n, p);
    matrix *result = m_init(m, p);
    rand_m_data_gen(A, 10.0);
    rand_m_data_gen(B, 10.0);
    m_mult(A, B, result);
    for (uint64_t i = 0; i < m; i++) {
        for (uint64_t j = 0; j < p; j++) {
            complex sum = c_zero();
            for (uint64_t k = 0; k < n; k++) {
                sum = c_add(sum, c_mult(A->data[i * n + k], B->data[k * p + j]));
            }
            CU_ASSERT_DOUBLE_EQUAL(result->data[i * p + j].Re, sum.Re, 1e-6);
            CU_ASSERT_DOUBLE_EQUAL(result->data[i * p + j].Im, sum.Im, 1e-6);
        }
    }
    m_free(A);
    m_free(B);
    m_free(result);
}
void test_m_scale() {
    uint64_t m = rand_dim_gen(100);
    uint64_t n = rand_dim_gen(100);
    matrix *A = m_init(m, n);
    matrix *result = m_init(m, n);
    complex alpha = rand_c_gen(10.0);
    rand_m_data_gen(A, 10.0);
    m_scale(A, alpha, result);
    for (uint64_t i = 0; i < m; i++) {
        for (uint64_t j = 0; j < n; j++) {
            CU_ASSERT_DOUBLE_EQUAL(result->data[i * n + j].Re,
                                   A->data[i * n + j].Re * alpha.Re -
                                       A->data[i * n + j].Im * alpha.Im,
                                   1e-6);
            CU_ASSERT_DOUBLE_EQUAL(result->data[i * n + j].Im,
                                   A->data[i * n + j].Re * alpha.Im +
                                       A->data[i * n + j].Im * alpha.Re,
                                   1e-6);
        }
    }
    m_free(A);
    m_free(result);
}
void test_m_scale_r() {
    uint64_t m = rand_dim_gen(100);
    uint64_t n = rand_dim_gen(100);
    matrix *A = m_init(m, n);
    matrix *result = m_init(m, n);
    double alpha = rand_dim_gen(10.0);
    rand_m_data_gen(A, 10.0);
    m_scale_r(A, alpha, result);
    for (uint64_t i = 0; i < m; i++) {
        for (uint64_t j = 0; j < n; j++) {
            CU_ASSERT_DOUBLE_EQUAL(result->data[i * n + j].Re,
                                   A->data[i * n + j].Re * alpha, 1e-6);
            CU_ASSERT_DOUBLE_EQUAL(result->data[i * n + j].Im,
                                   A->data[i * n + j].Im * alpha, 1e-6);
        }
    }
    m_free(A);
    m_free(result);
}
void test_m_conj() {
    uint64_t m = rand_dim_gen(100);
    uint64_t n = rand_dim_gen(100);
    matrix *A = m_init(m, n);
    matrix *result = m_init(m, n);
    rand_m_data_gen(A, 10.0);
    m_conj(A, result);
    for (uint64_t i = 0; i < m; i++) {
        for (uint64_t j = 0; j < n; j++) {
            CU_ASSERT_DOUBLE_EQUAL(result->data[i * n + j].Re, A->data[i * n + j].Re,
                                   1e-6);
            CU_ASSERT_DOUBLE_EQUAL(result->data[i * n + j].Im, -A->data[i * n + j].Im,
                                   1e-6);
        }
    }
    m_free(A);
    m_free(result);
}
/*
void test_m_Re() {
    uint64_t m = rand_dim_gen(100);
    uint64_t n = rand_dim_gen(100);
    matrix *A = m_init(m, n);
    matrix *result = m_init(m, n);
    rand_m_data_gen(A, 10.0);
    m_Re(A, result);
    for (uint64_t i = 0; i < m; i++) {
        for (uint64_t j = 0; j < n; j++) {
            CU_ASSERT_DOUBLE_EQUAL(result->data[i * n + j].Re, A->data[i * n + j].Re,
1e-6); CU_ASSERT_DOUBLE_EQUAL(result->data[i * n + j].Im, 0.0, 1e-6);
        }
    }
    m_free(A);
    m_free(result);
}
void test_m_Im() {
    uint64_t m = rand_dim_gen(100);
    uint64_t n = rand_dim_gen(100);
    matrix *A = m_init(m, n);
    matrix *result = m_init(m, n);
    rand_m_data_gen(A, 10.0);
    m_Im(A, result);
    for (uint64_t i = 0; i < m; i++) {
        for (uint64_t j = 0; j < n; j++) {
            CU_ASSERT_DOUBLE_EQUAL(result->data[i * n + j].Re, 0.0, 1e-6);
            CU_ASSERT_DOUBLE_EQUAL(result->data[i * n + j].Im, A->data[i * n + j].Im,
1e-6);
        }
    }
    m_free(A);
    m_free(result);
}
void test_m_map() {
    uint64_t m = rand_dim_gen(100);
    uint64_t n = rand_dim_gen(100);
    matrix *A = m_init(m, n);
    matrix *result = m_init(m, n);
    rand_m_data_gen(A, 10.0);
    m_map(A, c_sqrt, result);
    for (uint64_t i = 0; i < m; i++) {
        for (uint64_t j = 0; j < n; j++) {
            complex expected = c_sqrt(A->data[i * n + j]);
            CU_ASSERT_DOUBLE_EQUAL(result->data[i * n + j].Re, expected.Re, 1e-6);
            CU_ASSERT_DOUBLE_EQUAL(result->data[i * n + j].Im, expected.Im, 1e-6);
        }
    }
    m_free(A);
    m_free(result);
} */
void test_m_transpose() {
    uint64_t m = rand_dim_gen(100);
    uint64_t n = rand_dim_gen(100);
    matrix *A = m_init(m, n);
    matrix *result = m_init(n, m);
    rand_m_data_gen(A, 10.0);
    m_transpose(A, result);
    for (uint64_t i = 0; i < m; i++) {
        for (uint64_t j = 0; j < n; j++) {
            CU_ASSERT_DOUBLE_EQUAL(result->data[j * m + i].Re, A->data[i * n + j].Re,
                                   1e-6);
            CU_ASSERT_DOUBLE_EQUAL(result->data[j * m + i].Im, A->data[i * n + j].Im,
                                   1e-6);
        }
    }
    m_free(A);
    m_free(result);
}
void test_m_hermitian() {
    uint64_t m = rand_dim_gen(100);
    matrix *A = m_init(m, m);
    matrix *result = m_init(m, m);
    rand_m_data_gen(A, 10.0);
    m_hermitian(A, result);
    for (uint64_t i = 0; i < m; i++) {
        for (uint64_t j = 0; j < m; j++) {
            CU_ASSERT_DOUBLE_EQUAL(result->data[i * m + j].Re, A->data[j * m + i].Re,
                                   1e-6);
            CU_ASSERT_DOUBLE_EQUAL(result->data[i * m + j].Im, -A->data[j * m + i].Im,
                                   1e-6);
        }
    }
    m_free(A);
    m_free(result);
}
void test_m_norm() {
    uint64_t m = rand_dim_gen(100);
    uint64_t n = rand_dim_gen(100);
    matrix *A = m_init(m, n);
    rand_m_data_gen(A, 10.0);
    double result = m_norm(A);
    double expected = 0.0;
    for (uint64_t i = 0; i < m; i++) {
        for (uint64_t j = 0; j < n; j++) {
            expected += c_norm2(A->data[i * n + j]);
        }
    }
    expected = sqrt(expected);
    CU_ASSERT_DOUBLE_EQUAL(result, expected, 1e-6);
    m_free(A);
}
void test_m_norm2() {
    uint64_t m = rand_dim_gen(100);
    uint64_t n = rand_dim_gen(100);
    matrix *A = m_init(m, n);
    rand_m_data_gen(A, 10.0);
    double result = m_norm2(A);
    double expected = 0.0;
    for (uint64_t i = 0; i < m; i++) {
        for (uint64_t j = 0; j < n; j++) {
            expected += c_norm2(A->data[i * n + j]);
        }
    }
    CU_ASSERT_DOUBLE_EQUAL(result, expected, 1e-6);
    m_free(A);
}

void test_v_outer() {
    uint64_t m = rand_dim_gen(100);
    vector *a = v_init(m);
    vector *b = v_init(m);
    matrix *result = m_init(m, m);
    rand_v_data_gen(a, 10.0);
    rand_v_data_gen(b, 10.0);
    v_outer(a, b, result);
    for (uint64_t i = 0; i < m; i++) {
        for (uint64_t j = 0; j < m; j++) {
            CU_ASSERT_DOUBLE_EQUAL(
                result->data[i * m + j].Re,
                a->data[i].Re * b->data[j].Re + a->data[i].Im * b->data[j].Im, 1e-6);
            CU_ASSERT_DOUBLE_EQUAL(
                result->data[i * m + j].Im,
                -a->data[i].Re * b->data[j].Im + a->data[i].Im * b->data[j].Re, 1e-6);
        }
    }
    v_free(a);
    v_free(b);
    m_free(result);
}
int main() {
    rand_init_seed();
    printf("Running Matrix Tests...\n");
    CU_initialize_registry();
    CU_pSuite suite = CU_add_suite("Matrix Tests", NULL, NULL);

    CU_add_test(suite, "test_m_add", test_m_add);
    CU_add_test(suite, "test_m_sub", test_m_sub);
    CU_add_test(suite, "test_m_mult_e", test_m_mult_e);
    CU_add_test(suite, "test_m_div_e", test_m_div_e);
    CU_add_test(suite, "test_m_mult", test_m_mult);
    CU_add_test(suite, "test_m_scale", test_m_scale);
    CU_add_test(suite, "test_m_scale_r", test_m_scale_r);
    CU_add_test(suite, "test_m_conj", test_m_conj);
    // CU_add_test(suite, "test_m_Re", test_m_Re);
    // CU_add_test(suite, "test_m_Im", test_m_Im);
    // CU_add_test(suite, "test_m_map", test_m_map);
    CU_add_test(suite, "test_m_transpose", test_m_transpose);
    CU_add_test(suite, "test_m_hermitian", test_m_hermitian);
    CU_add_test(suite, "test_m_norm", test_m_norm);
    CU_add_test(suite, "test_m_norm2", test_m_norm2);
    CU_add_test(suite, "test_v_outer", test_v_outer);

    CU_basic_run_tests();
    CU_cleanup_registry();

    return 0;
}
