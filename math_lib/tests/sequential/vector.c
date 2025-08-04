#include "../../headers/general.h"
#include "../../headers/sequential.h"
#include "../../headers/utils.h"
#include <CUnit/Basic.h>
#include <CUnit/CUnit.h>
#include <stdint.h>
#include <stdio.h>

void test_v_add() {
    uint64_t m = rand_dim_gen(100);
    vector *a = v_init(m);
    vector *b = v_init(m);
    vector *result = v_init(m);
    rand_v_data_gen(a, 10.0);
    rand_v_data_gen(b, 10.0);
    v_add(a, b, result);
    for (uint64_t i = 0; i < m; i++) {
        CU_ASSERT_DOUBLE_EQUAL(result->data[i].Re, a->data[i].Re + b->data[i].Re, 1e-6);
        CU_ASSERT_DOUBLE_EQUAL(result->data[i].Im, a->data[i].Im + b->data[i].Im, 1e-6);
    }
    v_free(a);
    v_free(b);
    v_free(result);
}
void test_v_sub() {
    uint64_t m = rand_dim_gen(100);
    vector *a = v_init(m);
    vector *b = v_init(m);
    vector *result = v_init(m);
    rand_v_data_gen(a, 10.0);
    rand_v_data_gen(b, 10.0);
    v_sub(a, b, result);
    for (uint64_t i = 0; i < m; i++) {
        CU_ASSERT_DOUBLE_EQUAL(result->data[i].Re, a->data[i].Re - b->data[i].Re, 1e-6);
        CU_ASSERT_DOUBLE_EQUAL(result->data[i].Im, a->data[i].Im - b->data[i].Im, 1e-6);
    }
    v_free(a);
    v_free(b);
    v_free(result);
}
void test_v_mult_e() {
    uint64_t m = rand_dim_gen(100);
    vector *a = v_init(m);
    vector *b = v_init(m);
    vector *result = v_init(m);
    rand_v_data_gen(a, 10.0);
    rand_v_data_gen(b, 10.0);
    v_mult_e(a, b, result);
    for (uint64_t i = 0; i < m; i++) {
        CU_ASSERT_DOUBLE_EQUAL(
            result->data[i].Re,
            a->data[i].Re * b->data[i].Re - a->data[i].Im * b->data[i].Im, 1e-6);
        CU_ASSERT_DOUBLE_EQUAL(
            result->data[i].Im,
            a->data[i].Re * b->data[i].Im + a->data[i].Im * b->data[i].Re, 1e-6);
    }
    v_free(a);
    v_free(b);
    v_free(result);
}
void test_v_div_e() {
    uint64_t m = rand_dim_gen(100);
    vector *a = v_init(m);
    vector *b = v_init(m);
    vector *result = v_init(m);
    rand_v_data_gen(a, 10.0);
    rand_v_data_gen(b, 10.0);
    for (uint64_t i = 0; i < m; i++) {
        if (b->data[i].Re == 0 && b->data[i].Im == 0) {
            b->data[i].Re = 1e-6; // avoid division by zero
        }
    }
    v_div_e(a, b, result);
    for (uint64_t i = 0; i < m; i++) {
        double denom = b->data[i].Re * b->data[i].Re + b->data[i].Im * b->data[i].Im;
        CU_ASSERT_DOUBLE_EQUAL(
            result->data[i].Re,
            (a->data[i].Re * b->data[i].Re + a->data[i].Im * b->data[i].Im) / denom,
            1e-6);
        CU_ASSERT_DOUBLE_EQUAL(
            result->data[i].Im,
            (a->data[i].Im * b->data[i].Re - a->data[i].Re * b->data[i].Im) / denom,
            1e-6);
    }
    v_free(a);
    v_free(b);
    v_free(result);
}
void test_v_scale() {
    uint64_t m = rand_dim_gen(100);
    vector *a = v_init(m);
    vector *result = v_init(m);
    complex alpha = rand_c_gen(10.0);
    rand_v_data_gen(a, 10.0);
    v_scale(a, alpha, result);
    for (uint64_t i = 0; i < m; i++) {
        CU_ASSERT_DOUBLE_EQUAL(result->data[i].Re,
                               a->data[i].Re * alpha.Re - a->data[i].Im * alpha.Im, 1e-6);
        CU_ASSERT_DOUBLE_EQUAL(result->data[i].Im,
                               a->data[i].Re * alpha.Im + a->data[i].Im * alpha.Re, 1e-6);
    }
    v_free(a);
    v_free(result);
}
void test_v_scale_r() {
    uint64_t m = rand_dim_gen(100);
    vector *a = v_init(m);
    vector *result = v_init(m);
    double alpha = rand_c_Int_gen(10.0).Re; // using Re part
    rand_v_data_gen(a, 10.0);
    v_scale_r(a, alpha, result);
    for (uint64_t i = 0; i < m; i++) {
        CU_ASSERT_DOUBLE_EQUAL(result->data[i].Re, a->data[i].Re * alpha, 1e-6);
        CU_ASSERT_DOUBLE_EQUAL(result->data[i].Im, a->data[i].Im * alpha, 1e-6);
    }
    v_free(a);
    v_free(result);
}
void test_v_conj() {
    uint64_t m = rand_dim_gen(100);
    vector *a = v_init(m);
    vector *result = v_init(m);
    rand_v_data_gen(a, 10.0);
    v_conj(a, result);
    for (uint64_t i = 0; i < m; i++) {
        CU_ASSERT_DOUBLE_EQUAL(result->data[i].Re, a->data[i].Re, 1e-6);
        CU_ASSERT_DOUBLE_EQUAL(result->data[i].Im, -a->data[i].Im, 1e-6);
    }
    v_free(a);
    v_free(result);
}
void test_v_Re() {
    uint64_t m = rand_dim_gen(100);
    vector *a = v_init(m);
    vector *result = v_init(m);
    rand_v_data_gen(a, 10.0);
    v_Re(a, result);
    for (uint64_t i = 0; i < m; i++) {
        CU_ASSERT_DOUBLE_EQUAL(result->data[i].Re, a->data[i].Re, 1e-6);
        CU_ASSERT_DOUBLE_EQUAL(result->data[i].Im, 0.0, 1e-6);
    }
    v_free(a);
    v_free(result);
}
void test_v_Im() {
    uint64_t m = rand_dim_gen(100);
    vector *a = v_init(m);
    vector *result = v_init(m);
    rand_v_data_gen(a, 10.0);
    v_Im(a, result);
    for (uint64_t i = 0; i < m; i++) {
        CU_ASSERT_DOUBLE_EQUAL(result->data[i].Re, 0.0, 1e-6);
        CU_ASSERT_DOUBLE_EQUAL(result->data[i].Im, a->data[i].Im, 1e-6);
    }
    v_free(a);
    v_free(result);
}
void test_v_map() {
    uint64_t m = rand_dim_gen(100);
    vector *a = v_init(m);
    vector *result = v_init(m);
    rand_v_data_gen(a, 10.0);
    v_map(a, c_sqrt, result);
    for (uint64_t i = 0; i < m; i++) {
        complex expected = c_sqrt(a->data[i]);
        CU_ASSERT_DOUBLE_EQUAL(result->data[i].Re, expected.Re, 1e-6);
        CU_ASSERT_DOUBLE_EQUAL(result->data[i].Im, expected.Im, 1e-6);
    }
    v_free(a);
    v_free(result);
}
void test_v_equal() {
    uint64_t m = rand_dim_gen(100);
    vector *a = v_init(m);
    vector *b = v_init(m);
    vector *result = v_init(m);
    rand_v_data_gen(a, 10.0);
    rand_v_data_gen(b, 10.0);
    v_equal(a, b, result);
    for (uint64_t i = 0; i < m; i++) {
        if (a->data[i].Re == b->data[i].Re && a->data[i].Im == b->data[i].Im) {
            CU_ASSERT_DOUBLE_EQUAL(result->data[i].Re, 1.0, 1e-6);
            CU_ASSERT_DOUBLE_EQUAL(result->data[i].Im, 0.0, 1e-6);
        } else {
            CU_ASSERT_DOUBLE_EQUAL(result->data[i].Re, 0.0, 1e-6);
            CU_ASSERT_DOUBLE_EQUAL(result->data[i].Im, 0.0, 1e-6);
        }
    }
    v_free(a);
    v_free(b);
    v_free(result);
}
void test_v_dot_prod() {
    uint64_t m = rand_dim_gen(100);
    vector *a = v_init(m);
    vector *b = v_init(m);
    rand_v_data_gen(a, 10.0);
    rand_v_data_gen(b, 10.0);
    complex result = v_dot_prod(a, b);
    complex expected = c_zero();
    for (uint64_t i = 0; i < m; i++) {
        expected = c_add(expected, c_mult(a->data[i], c_conj(b->data[i])));
    }
    CU_ASSERT_DOUBLE_EQUAL(result.Re, expected.Re, 1e-6);
    CU_ASSERT_DOUBLE_EQUAL(result.Im, expected.Im, 1e-6);
    v_free(a);
    v_free(b);
}
void test_v_norm() {
    uint64_t m = rand_dim_gen(100);
    vector *a = v_init(m);
    rand_v_data_gen(a, 10.0);
    double result = v_norm(a);
    double expected = 0.0;
    for (uint64_t i = 0; i < m; i++) {
        expected += c_norm2(a->data[i]);
    }
    expected = sqrt(expected);
    CU_ASSERT_DOUBLE_EQUAL(result, expected, 1e-6);
    v_free(a);
}
void test_v_norm2() {
    uint64_t m = rand_dim_gen(100);
    vector *a = v_init(m);
    rand_v_data_gen(a, 10.0);
    double result = v_norm2(a);
    double expected = 0.0;
    for (uint64_t i = 0; i < m; i++) {
        expected += c_norm2(a->data[i]);
    }
    CU_ASSERT_DOUBLE_EQUAL(result, expected, 1e-6);
    v_free(a);
}

int main() {
    rand_init_seed();
    printf("Running vector operations tests...\n");
    CU_initialize_registry();
    CU_pSuite suite = CU_add_suite("Vector Operations", NULL, NULL);

    CU_add_test(suite, "test_v_add", test_v_add);
    CU_add_test(suite, "test_v_sub", test_v_sub);
    CU_add_test(suite, "test_v_mult_e", test_v_mult_e);
    CU_add_test(suite, "test_v_div_e", test_v_div_e);
    CU_add_test(suite, "test_v_scale", test_v_scale);
    CU_add_test(suite, "test_v_scale_r", test_v_scale_r);
    CU_add_test(suite, "test_v_conj", test_v_conj);
    CU_add_test(suite, "test_v_Re", test_v_Re);
    CU_add_test(suite, "test_v_Im", test_v_Im);
    CU_add_test(suite, "test_v_map", test_v_map);
    CU_add_test(suite, "test_v_equal", test_v_equal);
    CU_add_test(suite, "test_v_dot_prod", test_v_dot_prod);
    CU_add_test(suite, "test_v_norm", test_v_norm);
    CU_add_test(suite, "test_v_norm2", test_v_norm2);

    CU_basic_run_tests();
    CU_cleanup_registry();

    return 0;
}
