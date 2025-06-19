#include "../headers/general.h"
#include "../headers/sequential.h"
#include "../headers/utils.h"
#include <stdint.h>

void rand_m_data_gen(matrix *A, double range) {
    uint64_t m = A->m;
    uint64_t n = A->n;
    for (uint64_t i = 0; i < m * n; i++) {
        A->data[i] = rand_c_gen(range);
    }
}

void rand_m_data_Re_gen(matrix *A, double range) {
    uint64_t m = A->m;
    uint64_t n = A->n;
    for (uint64_t i = 0; i < m * n; i++) {
        A->data[i] = c_Re(rand_c_gen(range));
    }
}
