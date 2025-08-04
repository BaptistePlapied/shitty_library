#include "../headers/general.h"
#include "../headers/sequential.h"
#include "../headers/utils.h"
#include <stdint.h>
#include <stdlib.h>

void rand_v_data_gen(vector *a, double range) {
    uint64_t m = a->m;
    for (uint64_t i = 0; i < m; i++) {
        a->data[i] = rand_c_gen(range);
    }
}

void rand_v_data_Re_gen(vector *a, double range) {
    uint64_t m = a->m;
    for (uint64_t i = 0; i < m; i++) {
        a->data[i] = c_Re(rand_c_gen(range));
    }
}

void rand_v_data_Int_gen(vector *a, double range) {
    uint64_t m = a->m;
    for (uint64_t i = 0; i < m; i++) {
        a->data[i] = c_Re(rand_c_Int_gen(range));
    }
}
