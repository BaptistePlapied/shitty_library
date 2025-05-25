#include "general.h"
#include <stdint.h>

// all.c
void rand_init_seed();
uint64_t rand_dim_gen(uint64_t max_size);
complex rand_c_gen(double range);

// vector.c
void rand_v_data_gen(vector *a, double range);

// matrix.c
void rand_m_data_gen(matrix *A, double range);

// tensor.c
