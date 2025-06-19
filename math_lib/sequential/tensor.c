#include "../headers/general.h"
#include "../headers/sequential.h"
#include <stdlib.h>

void t_free(tensor *T) {
    if (!T)
        return;
    if (T->data) {
        free(T->data);
        T->data = NULL;
    }
    if (T->shape) {
        free(T->shape);
        T->shape = NULL;
    }
    // Free the tensor itself if dynamically allocated
    free(T);
}
