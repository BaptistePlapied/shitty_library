
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

void test_m_LU() {
    uint64_t dim_m = rand_dim_gen(10);
    // uint64_t dim_n = rand_dim_gen(15);
    matrix *A = m_init(dim_m, dim_m);
    rand_m_data_gen(A, 100);
    printf("A : \n");
    m_printf(A);
    matrix *U = m_init(dim_m, dim_m);
    matrix *L = m_init(dim_m, dim_m);
    matrix *P = m_init(dim_m, dim_m);
    lu_pp(A, L, U, P, NULL);
    printf("L : \n");
    m_printf(L);
    printf("U : \n");
    m_printf(U);
    printf("P : \n");
    m_printf(P);

    matrix *PA = m_init(dim_m, dim_m);
    matrix *LU = m_init(dim_m, dim_m);
    m_mult(P, A, PA);
    m_mult(L, U, LU);

    printf("PA : \n");
    m_printf(PA);
    printf("LU : \n");
    m_printf(LU);
}

void test_m_QR() {
    uint64_t dim_m = rand_dim_gen(5);
    // uint64_t dim_n = rand_dim_gen(15);
    matrix *A = m_init(dim_m, dim_m);
    rand_m_data_Re_gen(A, 10);
    printf("A : \n");
    m_printf(A);
    matrix *Q = m_init(dim_m, dim_m);
    matrix *R = m_init(dim_m, dim_m);
    m_mgs_qr(A, Q, R);
    printf("Q : \n");
    m_printf(Q);
    printf("R : \n");
    m_printf(R);

    printf("Q * R : \n");
    m_printf(m_mult(Q, R, 0));
}

void test_eigenvalues() {
    uint64_t dim_m = rand_dim_gen(5);
    matrix *A = m_init(dim_m, dim_m);
    rand_m_data_Int_gen(A, 10);
    printf("A : \n");
    m_printf(A);
    vector *eigenvals = qr_eigenvalues(A, 1e-8, 50, NULL);
    printf("\n eigenvals :\n\t");
    v_printf(eigenvals);
}

void test_hessenberg_reduction() {
    uint64_t dim_m = rand_dim_gen(5);
    matrix *A = m_init(dim_m, dim_m);
    rand_m_data_Int_gen(A, 10);
    printf("A : \n");
    m_printf(A);
    matrix *A_res = m_init(dim_m, dim_m);
    matrix *Q_tot = m_init(dim_m, dim_m);
    m_hessen_reduc(A, A_res, Q_tot);
    printf("A_res :\n");
    m_printf(A_res);
    printf("Q_tot :\n");
    m_printf(Q_tot);
    matrix *I = m_mult(m_track(m_transpose(Q_tot, NULL)), Q_tot, NULL);
    printf("I : \n");
    m_printf(I);
    uint64_t max_index = track_count();
    printf("max_index : %llu\n", max_index);
    matrix *Q_Tot = m_track_get(max_index - 1);
    matrix *A_rescon = m_mult(Q_Tot, m_track(m_mult(A, Q_tot, NULL)), NULL);
    printf("A_res reconstitute: \n");
    m_printf(A_rescon);
    printf("max_index : %llu\n", track_count());
    track_clear_all();
    printf("max_index : %llu\n", track_count());
}
int main() {
    rand_init_seed();
    test_m_QR();
    test_eigenvalues();
    // test_hessenberg_reduction();
}
