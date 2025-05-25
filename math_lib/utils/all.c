#include "../headers/general.h"
#include "../headers/utils.h"
#include <stdlib.h>
#include <time.h>

void rand_init_seed() {
    struct timespec ts;
    clock_gettime(CLOCK_REALTIME, &ts); // You could also use CLOCK_REALTIME
    unsigned int seed = (unsigned int)(ts.tv_nsec ^ ts.tv_sec);
    srand(seed);
}

uint64_t rand_dim_gen(uint64_t max_size) { return ((uint64_t)rand() % max_size) + 1; }

complex rand_c_gen(double range) {
    return (complex){((double)rand() / RAND_MAX) * 2 * range - range,
                     ((double)rand() / RAND_MAX) * 2 * range - range};
}
