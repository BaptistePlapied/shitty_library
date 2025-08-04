#include "../../headers/general.h"
#include "../../headers/sequential.h"
#include "../../headers/utils.h"
#include <CUnit/Basic.h>
#include <CUnit/CUnit.h>
#include <stdint.h>
#include <stdio.h>

void test_c_add() {
    complex a = rand_c_gen(10);
    complex b = rand_c_gen(10);
    complex result = c_add(a, b);
    complex expected = {a.Re + b.Re, a.Im + b.Im};
    CU_ASSERT_DOUBLE_EQUAL(result.Re, expected.Re, 1e-6);
    CU_ASSERT_DOUBLE_EQUAL(result.Im, expected.Im, 1e-6);
}
void test_c_sub() {
    complex a = rand_c_gen(10);
    complex b = rand_c_gen(10);
    complex result = c_sub(a, b);
    complex expected = {a.Re - b.Re, a.Im - b.Im};
    CU_ASSERT_DOUBLE_EQUAL(result.Re, expected.Re, 1e-6);
    CU_ASSERT_DOUBLE_EQUAL(result.Im, expected.Im, 1e-6);
}
void test_c_scale() {
    complex a = rand_c_gen(10);
    double scale = rand_dim_gen(10);
    complex result = c_scale(a, scale);
    complex expected = {a.Re * scale, a.Im * scale};
    CU_ASSERT_DOUBLE_EQUAL(result.Re, expected.Re, 1e-6);
    CU_ASSERT_DOUBLE_EQUAL(result.Im, expected.Im, 1e-6);
}
void test_c_zero() {
    complex result = c_zero();
    complex expected = {0.0, 0.0};
    CU_ASSERT_DOUBLE_EQUAL(result.Re, expected.Re, 1e-6);
    CU_ASSERT_DOUBLE_EQUAL(result.Im, expected.Im, 1e-6);
}
void test_c_real() {
    double a = rand_dim_gen(10);
    complex result = c_real(a);
    complex expected = {a, 0.0};
    CU_ASSERT_DOUBLE_EQUAL(result.Re, expected.Re, 1e-6);
    CU_ASSERT_DOUBLE_EQUAL(result.Im, expected.Im, 1e-6);
}
void test_c_Re() {
    complex a = rand_c_gen(10);
    complex result = c_Re(a);
    complex expected = {a.Re, 0.0};
    CU_ASSERT_DOUBLE_EQUAL(result.Re, expected.Re, 1e-6);
    CU_ASSERT_DOUBLE_EQUAL(result.Im, expected.Im, 1e-6);
}
void test_c_Im() {
    complex a = rand_c_gen(10);
    complex result = c_Im(a);
    complex expected = {0.0, a.Im};
    CU_ASSERT_DOUBLE_EQUAL(result.Re, expected.Re, 1e-6);
    CU_ASSERT_DOUBLE_EQUAL(result.Im, expected.Im, 1e-6);
}
void test_c_conj() {
    complex a = rand_c_gen(10);
    complex result = c_conj(a);
    complex expected = {a.Re, -a.Im};
    CU_ASSERT_DOUBLE_EQUAL(result.Re, expected.Re, 1e-6);
    CU_ASSERT_DOUBLE_EQUAL(result.Im, expected.Im, 1e-6);
}
void test_c_abs() {
    complex a = rand_c_gen(10);
    complex result = c_abs(a);
    complex expected = {fabs(a.Re), fabs(a.Im)};
    CU_ASSERT_DOUBLE_EQUAL(result.Re, expected.Re, 1e-6);
    CU_ASSERT_DOUBLE_EQUAL(result.Im, expected.Im, 1e-6);
}
void test_c_norm2() {
    complex a = rand_c_gen(10);
    double result = c_norm2(a);
    double expected = a.Re * a.Re + a.Im * a.Im;
    CU_ASSERT_DOUBLE_EQUAL(result, expected, 1e-6);
}
void test_c_norm() {
    complex a = rand_c_gen(10);
    double result = c_norm(a);
    double expected = sqrt(c_norm2(a));
    CU_ASSERT_DOUBLE_EQUAL(result, expected, 1e-6);
}
void test_c_arg() {
    complex a = rand_c_gen(10);
    double result = c_arg(a);
    double expected = atan2(a.Im, a.Re);
    CU_ASSERT_DOUBLE_EQUAL(result, expected, 1e-6);
}
void test_c_mult() {
    complex a = rand_c_gen(10);
    complex b = rand_c_gen(10);
    complex result = c_mult(a, b);
    complex expected = {a.Re * b.Re - a.Im * b.Im, a.Re * b.Im + a.Im * b.Re};
    CU_ASSERT_DOUBLE_EQUAL(result.Re, expected.Re, 1e-6);
    CU_ASSERT_DOUBLE_EQUAL(result.Im, expected.Im, 1e-6);
}
void test_c_div() {
    complex a = rand_c_gen(10);
    complex b = rand_c_gen(10);
    if (c_equal(b, c_zero())) {
        CU_ASSERT_TRUE(1); // Skip division by zero
        return;
    }
    complex result = c_div(a, b);
    complex expected = {(a.Re * b.Re + a.Im * b.Im) / (b.Re * b.Re + b.Im * b.Im),
                        (a.Im * b.Re - a.Re * b.Im) / (b.Re * b.Re + b.Im * b.Im)};
    CU_ASSERT_DOUBLE_EQUAL(result.Re, expected.Re, 1e-6);
    CU_ASSERT_DOUBLE_EQUAL(result.Im, expected.Im, 1e-6);
}
void test_c_sqrt() {
    complex a = rand_c_gen(10);
    complex result = c_sqrt(a);
    complex expected = {sqrt((a.Re + sqrt(c_norm2(a))) / 2),
                        (a.Im >= 0 ? 1 : -1) * sqrt((-a.Re + sqrt(c_norm2(a))) / 2)};
    CU_ASSERT_DOUBLE_EQUAL(result.Re, expected.Re, 1e-6);
    CU_ASSERT_DOUBLE_EQUAL(result.Im, expected.Im, 1e-6);
}
void test_c_exp() {
    complex a = rand_c_gen(10);
    complex result = c_exp(a);
    complex expected = {exp(a.Re) * cos(a.Im), exp(a.Re) * sin(a.Im)};
    CU_ASSERT_DOUBLE_EQUAL(result.Re, expected.Re, 1e-6);
    CU_ASSERT_DOUBLE_EQUAL(result.Im, expected.Im, 1e-6);
}
void test_c_log() {
    complex a = rand_c_gen(10);
    if (c_equal(a, c_zero())) {
        CU_ASSERT_TRUE(1); // Skip log of zero
        return;
    }
    complex result = c_log(a);
    complex expected = {log(c_norm(a)), atan2(a.Im, a.Re)};
    CU_ASSERT_DOUBLE_EQUAL(result.Re, expected.Re, 1e-6);
    CU_ASSERT_DOUBLE_EQUAL(result.Im, expected.Im, 1e-6);
}
void test_c_equal() {
    complex a = rand_c_gen(10);
    complex b = a;
    uint8_t result = c_equal(a, b);
    CU_ASSERT_TRUE(result);
    b.Re += 1e-6; // Slightly modify b
    result = c_equal(a, b);
    CU_ASSERT_FALSE(result);
}
void test_c_compare() {
    complex a = rand_c_gen(10);
    complex b = a;
    uint8_t cmp_result = c_compare(a, b);
    CU_ASSERT_EQUAL(cmp_result, 1); // Should be equal
    if (a.Re > 0) {
        b.Re += 1e-7; // Slightly modify a
    } else {
        b.Re -= 1e-7; // Slightly modify a
    }
    cmp_result = c_compare(a, b);
    CU_ASSERT_EQUAL(cmp_result, 0); // a < b
    if (a.Re > 0) {
        a.Re += 2e-7; // Further modify a
    } else {
        a.Re -= 2e-7; // Further modify a
    }
    cmp_result = c_compare(a, b);
    CU_ASSERT_EQUAL(cmp_result, 2); // a > b
}

int main() {
    rand_init_seed();
    printf("Running Complex Operations Tests\n");
    CU_initialize_registry();
    CU_pSuite suite = CU_add_suite("Complex Operations Suite", NULL, NULL);

    CU_add_test(suite, "test_c_add", test_c_add);
    CU_add_test(suite, "test_c_sub", test_c_sub);
    CU_add_test(suite, "test_c_scale", test_c_scale);
    CU_add_test(suite, "test_c_zero", test_c_zero);
    CU_add_test(suite, "test_c_real", test_c_real);
    CU_add_test(suite, "test_c_Re", test_c_Re);
    CU_add_test(suite, "test_c_Im", test_c_Im);
    CU_add_test(suite, "test_c_conj", test_c_conj);
    CU_add_test(suite, "test_c_abs", test_c_abs);
    CU_add_test(suite, "test_c_norm2", test_c_norm2);
    CU_add_test(suite, "test_c_norm", test_c_norm);
    CU_add_test(suite, "test_c_arg", test_c_arg);
    CU_add_test(suite, "test_c_mult", test_c_mult);
    CU_add_test(suite, "test_c_div", test_c_div);
    CU_add_test(suite, "test_c_sqrt", test_c_sqrt);
    CU_add_test(suite, "test_c_exp", test_c_exp);
    CU_add_test(suite, "test_c_log", test_c_log);
    CU_add_test(suite, "test_c_equal", test_c_equal);
    CU_add_test(suite, "test_c_compare", test_c_compare);

    CU_basic_run_tests();
    CU_cleanup_registry();

    return 0;
}
