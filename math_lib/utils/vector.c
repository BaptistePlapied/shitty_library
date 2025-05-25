#include "../headers/general.h"
#include "../headers/utils.h"
#include <stdint.h>
#include <stdlib.h>

void rand_v_data_gen(vector *a, double range) {
    uint64_t m = a->m;
    for (uint64_t i = 0; i < m; i++) {
        a->data[i] = rand_c_gen(range);
    }
}
