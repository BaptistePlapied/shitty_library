#include <SDL2/SDL.h>
#include <math.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define DIM_MAT 5
#define DIM_VEC 5

#define Window_Dim_X 700
#define Window_Dim_Y 400

// don't know if i should make another one that support complex number
typedef struct {
    float real;
    float imag;
} Complex;

typedef struct {
    size_t dim;
    float v[DIM_VEC];
} Vect;

//
//  the function of Vect;
//

void Vect_init(Vect *vect, size_t dim, float *v);
void Vect_print(Vect *vect);

Vect *Vect_copy(Vect *vect_c, Vect *vect);
Vect *Vect_add(Vect *vect1, Vect *vect2, Vect *vectR);
Vect *Vect_sub(Vect *vect1, Vect *vect2, Vect *vectR);
Vect *Vect_mult_w(Vect *vect1, Vect *vect2, Vect *vectR);
Vect *Vect_div_w(Vect *vect1, Vect *vect2, Vect *vectR);
Vect *Vect_scal(Vect *vect, float scalar, Vect *vectR);
Vect *Vect_shift_dim(Vect *vect, int shift, Vect *vect_r);
Vect *Vect_shift(Vect *vect, int shift, Vect *vect_r);

float Vect_dot_prod(Vect *vect1, Vect *vect2);
float Vect_sum(Vect *vect);
float Vect_norm(Vect *vect, uint8_t root);

typedef struct {
    size_t dimR;                  // the actual number of row
    size_t dimC;                  // the actual number of col
    size_t active_num_row;        // the number of row that are !NULL
    size_t active_num_col;        // the number of col that are !NULL
    uint8_t active_row[DIM_MAT];  // row that are active == 1
    uint8_t active_col[DIM_MAT];  // col that are active == 1
    float mat[DIM_MAT * DIM_MAT]; // the matrix
} Mtx;

//
//  the function of Mtx
//

void Mtx_update_active(Mtx *m);
void Mtx_init(Mtx *m, size_t dimR, size_t dimC, float *mat);
void Mtx_copy(Mtx *m1, Mtx *m2);
void Mtx_print(Mtx *m);
void Mtx_lu_pp(Mtx *m, Mtx *mp, Mtx *ml, Mtx *mu, size_t *rowswap);
void Mtx_qr(Mtx *m, Mtx *q, Mtx *r);
void Mtx_cholesky(Mtx *m, Mtx *ml);

Mtx *Mtx_add(Mtx *m1, Mtx *m2, Mtx *mr);
Mtx *Mtx_sub(Mtx *m1, Mtx *m2, Mtx *mr);
Mtx *Mtx_scal(Mtx *m1, float scalar, Mtx *mr);
Mtx *Mtx_mult(Mtx *m1, Mtx *m2, Mtx *mr);
Mtx *Mtx_T(Mtx *m1, Mtx *mr);
Mtx *Mtx_mult_w(Mtx *m1, Mtx *m2, Mtx *mr);
Mtx *Mtx_div_w(Mtx *m1, Mtx *m2, Mtx *mr);
Mtx *Mtx_eye(Mtx *m);
Mtx *Mtx_inv(Mtx *m, Mtx *mr);
Mtx *Vect_mult_to_Mtx(Vect *vect1, Vect *vect2, Mtx *m);
Mtx *Mtx_copy_list(Mtx *m, float *mat);

Vect *Vect_mult_l_mtx(Vect *vect, Mtx *m, Vect *vect_r); // Vect at left of the Mtx
Vect *Vect_mult_g_mtx(Vect *vect, Mtx *m, Vect *vect_r); // Vect at right of the Mtx

float Mtx_det(Mtx *m);
float Mtx_trace(Mtx *m);

size_t Mtx_rank(Mtx *m);

uint8_t Mtx_is_sparse(Mtx *m);

extern Mtx mtx_projection;
Mtx *proj_init_Mtx();
Mtx *proj_update_all(float near, float far, float fov_deg, float aspect_ratio);
Mtx *proj_update_fov(float fov, float aspect_ratio);
Mtx *proj_update_nearfar(float near, float far);

typedef struct {
    union {
        Vect V;
        struct {
            float x, y;
        };
    };
} Vect2d;

void Vect2d_init(Vect2d *vect2d, float *vect);

typedef struct {
    union {
        Vect V;
        struct {
            float x, y, z;
        };
    };
} Vect3d;

void Vect3d_init(Vect3d *vect3d, float *vect);

typedef struct {
    union {
        Vect P;
        struct {
            uint8_t r, g, b, a;
        };
    };
} Pixel;

void Pixel_init(Pixel *pixel, float *vect);

typedef struct {
    union {
        Vect Q;  // Store as a 4D vector
        struct { // Allow direct named access
            float w, x, y, z;
        };
    };
} Quaternion;

void Quaternion_init(Quaternion *quaternion, float *vect);

typedef struct {
    Vect3d p[3];       // Vertex positions
    Vect3d normal[3];  // Per-vertex normals for smooth shading
    Vect3d tangent[3]; // Per-vertex tangents for texture mapping
    Vect3d uv[3];      // Texture coordinates for each vertex
    Vect3d bezier[3];  // Bézier control points for cubic interpolation
    Pixel color[3];  // Per-vertex color (optional)
} Triangle;

typedef struct {
    Triangle *tris;
    int numtris;
} Mesh;
