#include "the_header_file.h"

int main() {
    
    Vect v1;
    Vect v2;
    Vect v3;
    
    Vect v1_r;
    Vect v2_r;
    Vect v3_r;

    float v1_f[] ={1,2,13,0};
    float v2_f[] ={1,2,3,4};
    float v3_f[] ={5,6,7};

    Vect_init(&v1, 4, v1_f);
    Vect_init(&v2, 4, v2_f);
    Vect_init(&v3, 3, v3_f);

    Vect_shift_dim(&v1, 1, &v1_r);
    Vect_shift_dim(&v2, 1, &v2_r);
    Vect_shift(&v3, 0, &v3_r);
    
    v3_r.dim++;


    v1_r.dim = 4;
    printf("dim v1_r = %zu\t dim v2 = %zu\n",v1_r.dim,v2.dim);
    Mtx m1;
    Mtx_init(&m1, 4, 4, NULL);
    Vect_mult_to_Mtx(&v2, &v2, &m1);


    printf("vect : \n");
    Vect_print(&v1);
    Vect_print(&v2);
    Vect_print(&v3);

    printf("\n");
    
    printf("vect shifted : \n");
    Vect_print(&v1_r);
    Vect_print(&v2_r);
    Vect_print(&v3_r);

    printf("mtx : \n");
    Mtx_print(&m1);
    printf("vect l and g :\n");
    
    proj_init_Mtx();
    
    float near = 1.0f;
    float far = 100.0f;
    float fov = 80.0f;
    float aspect_ratio = (float)Window_Dim_Y / (float)Window_Dim_Y;
    proj_update_all(near, far, fov, aspect_ratio);
    
    printf("mtx : \n");
    Mtx_print(&mtx_projection);
    Vect_mult_l_mtx(&v1, &mtx_projection, &v2);
    Vect_mult_g_mtx(&v1, &mtx_projection, NULL);

    Vect_print(&v2);
    Vect_print(&v1);

    return 0;
}
