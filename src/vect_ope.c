#include "the_header_file.h"

void Vect_init(Vect *vect, size_t dim, float *v) {
    if(dim > DIM_VEC){
        fprintf(stderr,"Warning the trying to assign a dim of %zu but DIM_VECT = %d\n",
                dim, DIM_VEC);
        vect->dim = DIM_VEC;
    }else{
        vect->dim = dim;
    }
    if (v == NULL) {
        for (size_t i = 0; i < DIM_VEC; i++) {
            vect->v[i] = 0;
        }
    } else {
        for (size_t i = 0; i < DIM_VEC; i++) {
            if (i < vect->dim) {
                vect->v[i] = v[i];
            } else {
                vect->v[i] = 0.0;
            }
        }
    }
}

Vect *Vect_copy(Vect *vect_c, Vect *vect) {
    if (vect->dim != vect_c->dim) {
        fprintf(stderr,
                "Error: wrong dim for Vect_copy : vect_c->dim = %zu but vect->dim = %zu\n",
                vect_c->dim, vect->dim);
        return NULL;
    }
    for (size_t i = 0; i < vect->dim; i++) {
        vect->v[i] = vect_c->v[i];
    }
    return vect;
}

void Vect_print(Vect *vect) {
    printf("[ ");
    for (size_t i = 0; i < vect->dim; i++) {
        printf("%f", vect->v[i]);
        if (i < vect->dim - 1) {
            printf(", ");
        }
    }
    printf("]\n");
}

Vect *Vect_add(Vect *vect1, Vect *vect2, Vect *vectR) {
    if (vect1->dim != vect2->dim) {
        printf("Error: wrong dim for Vect_add : dim v1 = %zu but dim v2 = %zu\n",
               vect1->dim, vect2->dim);
        return NULL;
    }
    if (vectR == NULL) {
        for (size_t i = 0; i < vect1->dim; i++) {
            vect1->v[i] += vect2->v[i];
        }
        return vect1;
    } else {
        for (size_t i = 0; i < vect1->dim; i++) {
            vectR->v[i] = vect1->v[i] + vect2->v[i];
        }
        return vectR;
    }
}

Vect *Vect_sub(Vect *vect1, Vect *vect2, Vect *vectR) {
    if (vect1->dim != vect2->dim) {
        printf("Error: wrong dim for Vect_sub : dim v1 = %zu but dim v2 = %zu\n",
               vect1->dim, vect2->dim);
        return NULL;
    }
    if (vectR == NULL) {
        for (size_t i = 0; i < vect1->dim; i++) {
            vect1->v[i] -= vect2->v[i];
        }
        return vect1;
    } else {
        for (size_t i = 0; i < vect1->dim; i++) {
            vectR->v[i] = vect1->v[i] - vect2->v[i];
        }
        return vectR;
    }
}

Vect *Vect_mult_w(Vect *vect1, Vect *vect2, Vect *vectR) {
    if (vect1->dim != vect2->dim) {
        printf("Error: wrong dim for Vect_add : dim v1 = %zu but dim v2 = %zu\n",
               vect1->dim, vect2->dim);
        return NULL;
    }
    if (vectR == NULL) {
        for (size_t i = 0; i < vect1->dim; i++) {
            vect1->v[i] *= vect2->v[i];
        }
        return vect1;
    } else {
        for (size_t i = 0; i < vect1->dim; i++) {
            vectR->v[i] = vect1->v[i] * vect2->v[i];
        }
        return vectR;
    }
}

Vect *Vect_div_w(Vect *vect1, Vect *vect2, Vect *vectR) {
    if (vect1->dim != vect2->dim) {
        printf("Error: wrong dim for Vect_add : dim v1 = %zu but dim v2 = %zu\n",
               vect1->dim, vect2->dim);
        return NULL;
    }
    if (vectR == NULL) {
        for (size_t i = 0; i < vect1->dim; i++) {
            if (vect2->v[i] == 0) {
                vect1->v[i] = INFINITY;
            } else {
                vect1->v[i] = vect1->v[i] / vect2->v[i];
            }
        }
        return vect1;
    } else {
        for (size_t i = 0; i < vect1->dim; i++) {
            if (vect2->v[i] == 0) {
                vectR->v[i] = INFINITY;
            } else {
                vectR->v[i] = vect1->v[i] / vect2->v[i];
            }
        }
        return vectR;
    }
}

Vect *Vect_scal(Vect *vect, float scalar, Vect *vectR) {
    if (vectR == NULL) {
        for (size_t i = 0; i < vect->dim; i++) {
            vect->v[i] *= scalar;
        }
        return vect;
    } else {
        for (size_t i = 0; i < vect->dim; i++) {
            vectR->v[i] = vect->v[i] * scalar;
        }
        return vectR;
    }
}

float Vect_dot_prod(Vect *vect1, Vect *vect2) {
    if (vect1->dim != vect2->dim) {
        printf("Error: wrong dim for Vect_dot_prod : dim v1 = %zu but dim v2 = %zu\n",
               vect1->dim, vect2->dim);
        return 0.0;
    }
    float dot_prod = 0;
    for (size_t i = 0; i < vect1->dim; i++) {
        dot_prod += vect1->v[i] * vect2->v[i];
    }
    return dot_prod;
}

float Vect_sum(Vect *vect) {
    float sum = 0;
    for (size_t i = 0; i < vect->dim; i++) {
        sum += vect->v[i];
    }
    return sum;
}

float Vect_norm(Vect *vect, uint8_t root) {
    float norm = 0.0;
    for (size_t i = 0; i < vect->dim; i++) {
        norm += vect->v[i] * vect->v[i];
    }
    if (root) {
        return sqrtf(norm);
    } else {
        return norm;
    }
}

Mtx *Vect_mult_to_Mtx(Vect *vect1, Vect *vect2, Mtx *m) {
    if (vect1->dim != vect2->dim || vect1->dim != m->dimR || m->dimR != m->dimC) {
        printf("Error: wrong dim for Vect_mult_to_Mtx : dim v1 = %zu but dim v2 = %zu\n",
               vect1->dim, vect2->dim);
        return NULL;
    }
    if (m == NULL) {
        printf("Error: m doesn't exist\n");
        return NULL;
    }
    if (m->dimR != m->dimC) {
        printf("Error: wrong dim for m must be square : m dimR = %zu but m dimC = %zu\n",
               m->dimR, m->dimC);
        return NULL;
    }
    if (vect1->dim != m->dimR) {
        printf("Error: wrong dim between Vect and Mtx : dim v1 = %zu but m dim = %zu\n",
               vect1->dim, m->dimR);
        return NULL;
    }
    Mtx_init(m, vect1->dim, vect2->dim, NULL);
    for (size_t i = 0; i < vect1->dim; i++) {
        for (size_t j = 0; j < vect1->dim; j++) {
            m->mat[i * DIM_MAT + j] = vect1->v[i] * vect2->v[j];
        }
    }
    return m;
}


Vect *Vect_mult_l_mtx(Vect *vect, Mtx *m, Vect *vect_r){
    if (vect->dim != m->dimR){
        fprintf(stderr,"Error wrong dim for Vect_mult_l_mtx :"
                " vect->dim = %zu but m->dimR = %zu\n",
                vect->dim, m->dimR);
        exit(0);
        return NULL;
    }else if(vect_r == NULL){
        Vect temp_v;
        Vect_init(&temp_v, m->dimC, NULL);
        for (size_t i = 0; i < m->dimC; i++){
            float temp_f = 0;
            for(size_t j = 0; j < m->dimR; j++){
                temp_f += vect->v[j] * m->mat[i * DIM_MAT + j];
            }
            temp_v.v[i] = temp_f;
        }
        vect->dim = m->dimC;
        Vect_copy(&temp_v, vect);
        return vect;
    }else{
        vect_r->dim = m->dimC;
        for (size_t i = 0; i < m->dimC; i++){
            float temp_f = 0;
            for(size_t j = 0; j < m->dimR; j++){
                temp_f += vect->v[j] * m->mat[i * DIM_MAT + j];
            }
            vect_r->v[i] = temp_f;
        }
        return vect_r;
    }
}


Vect *Vect_mult_g_mtx(Vect *vect, Mtx *m, Vect *vect_r){
    if (vect->dim != m->dimC){
        fprintf(stderr,"Error wrong dim for Vect_mult_g_mtx :"
                " vect->dim = %zu but m->dimC = %zu\n",
                vect->dim, m->dimC);
        exit(0);
        return NULL;
    }else if(vect_r == NULL){
        Vect temp_v;
        Vect_init(&temp_v, m->dimR, NULL);
        for(size_t i = 0; i < m->dimR; i++){
            float temp_f = 0;
            for(size_t j = 0; j < m->dimC; j++){
                temp_f += vect->v[j] * m->mat[i * DIM_MAT + j];
            }
            temp_v.v[i] = temp_f;
        }
        vect->dim = m->dimR;
        Vect_copy(&temp_v, vect);
        return vect;
    }else{
        vect_r->dim = m->dimR;
        for(size_t i = 0; i < m->dimR; i++){
            float temp_f = 0;
            for(size_t j = 0; j < m->dimC; j++){
                temp_f += vect->v[j] * m->mat[i * DIM_MAT + j];
            }
            vect_r->v[i] = temp_f;
        }
        return vect_r;

    }
}

Vect *Vect_shift_dim(Vect *vect, int shift, Vect *vect_r){
    if ( shift >= 0){
        if (shift >= DIM_VEC){
            fprintf(stderr,"Error : in Vect_shift : trying to shift by" 
                    " %d but vect->dim is only = %zu "
                    "will return a vect bigger than the limit DIM_VEC = %d\n"
                    , shift, vect->dim, DIM_VEC);
            exit(0);
            return NULL;
        }
        if (vect_r == NULL){
            size_t dim_init = vect->dim;
            size_t dim = dim_init +  (size_t)shift;
            vect->dim = (dim > DIM_VEC) ? DIM_VEC : dim;
            Vect temp_v;
            Vect_init(&temp_v, vect->dim, NULL);
            for (size_t i = 0; i < dim_init; i++){
                temp_v.v[i+shift] = vect->v[i];
            }
            Vect_copy(&temp_v, vect);
            return vect;
        }else{
            size_t dim_init = vect->dim;
            size_t dim = dim_init +  (size_t)shift;
            vect_r->dim = (dim > DIM_VEC) ? DIM_VEC : dim;
            for (size_t i = 0; i < dim_init; i++){
                vect_r->v[i+shift] = vect->v[i];
            }
            return vect_r;
        }
    }else{
        if(vect->dim <= (size_t)(-shift)){
            fprintf(stderr,"Error: in Vect_shift : trying to shift by" 
                    " %d but vect->dim is only = %zu will return a vect of dim < 0\n"
                    ,shift,vect->dim);
            exit(0);
            return NULL;
        }
        if(vect_r == NULL){
            size_t dim = vect->dim + shift;
            Vect temp_v;
            Vect_init(&temp_v, dim, NULL);
            for(size_t i = 0; i < dim; i++){
                temp_v.v[i] = vect->v[i-shift];
            }
            vect->dim = dim;
            Vect_copy(&temp_v, vect);
            return vect;
        }else{
            size_t dim = vect->dim + shift;
            vect_r->dim = dim;
            for(size_t i = 0; i < dim; i++){
                vect_r->v[i] = vect->v[i-shift];
            }
            return vect_r;
        }
    }
}


Vect *Vect_shift(Vect *vect, int shift, Vect *vect_r){
    if ( shift >= 0){
        if (vect_r == NULL){
            size_t dim_init = vect->dim;
            Vect temp_v;
            Vect_init(&temp_v, vect->dim, NULL);
            for (size_t i = 0; i+shift < dim_init; i++){
                temp_v.v[i+shift] = vect->v[i];
            }
            Vect_copy(&temp_v, vect);
            return vect;
        }else{
            size_t dim_init = vect->dim;
            vect_r->dim = dim_init;
            for (size_t i = 0; i+shift < dim_init; i++){
                vect_r->v[i+shift] = vect->v[i];
            }
            return vect_r;
        }
    }else{
        if(vect_r == NULL){
            size_t dim = vect->dim;
            Vect temp_v;
            Vect_init(&temp_v, dim, NULL);
            for(size_t i = 0; i-shift < dim; i++){
                temp_v.v[i] = vect->v[i-shift];
            }
            Vect_copy(&temp_v, vect);
            return vect;
        }else{
            size_t dim = vect->dim;
            vect_r->dim = dim;
            for(size_t i = 0; i-shift < dim; i++){
                vect_r->v[i] = vect->v[i-shift];
            }
            return vect_r;
        }
    }
}
