
#include "headers/general.h"
#include "headers/sequential.h"
#include "headers/utils.h"
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

void test_v() {
    uint64_t dim_v = rand_dim_gen(5);
    vector *a = v_init(dim_v);
    rand_v_data_gen(a, 1);

    vector *b = v_init(dim_v);
    rand_v_data_gen(b, 1);

    printf("a : \n");
    v_printf(a);
    printf("b : \n");
    v_printf(b);
    printf("a / b = c : \n");
    vector *c = v_div_e(a, b, NULL);
    v_printf(c);
}
void test_m() {
    uint64_t dim_m = rand_dim_gen(5);
    uint64_t dim_n = rand_dim_gen(5);
    uint64_t dim_k = rand_dim_gen(5);
    matrix *A = m_init(dim_m, dim_k);
    rand_m_data_gen(A, 5);

    matrix *B = m_init(dim_k, dim_n);
    rand_m_data_gen(B, 5);

    printf("A : \n");
    m_printf(A);
    printf("B : \n");
    m_printf(B);
    printf("A * B = C : \n");
    matrix *c = m_mult(A, B, NULL);
    m_printf(c);
}

void test_m_echelon() {
    uint64_t dim_m = rand_dim_gen(15);
    // uint64_t dim_n = rand_dim_gen(15);
    matrix *A = m_init(dim_m, dim_m);
    rand_m_data_gen(A, 5);
    printf("A : \n");
    m_printf(A);
    printf("A echelon = C : \n");
    echelon_expr echelon;
    matrix *c = m_echelon(A, &echelon, NULL);
    m_printf(c);
    echelon_expr_printf(echelon);
    vector *b = v_init(dim_m);
    rand_v_data_gen(b, 5);
    vector *d = backsub(c, b, NULL);
    printf("backsub :\n\t");
    printf(" b : \n\t\t");
    v_printf(b);
    printf(" result : \n\t\t");
    v_printf(d);
}
int main() {
    rand_init_seed();
    test_m_echelon();
}
