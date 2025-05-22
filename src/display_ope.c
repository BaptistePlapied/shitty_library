#include "the_header_file.h"

// TODO: implement frustum culling !!!
// TODO: implement Z-buffer !!!
// TODO: implement occlusion culling

// projection matrix and update it

// float near = 1.0f;
// float far = 100.0f;
// float fov = 90.0f;
// float aspect_ratio = (float)Window_Dim_Y / Window_Dim_Y;
// float fov_rad = 1.0f / tanf(fov * 0.5f / 180.0f * M_PI);
Mtx mtx_projection;

Mtx *proj_init_Mtx() {
    Mtx_init(&mtx_projection, 4, 4, NULL);
    mtx_projection.mat[0 + (0 * DIM_MAT)] = 1; // fAspectRatio * fFovRad;
    mtx_projection.mat[1 + (1 * DIM_MAT)] = 1; // fFovRad;
    mtx_projection.mat[2 + (2 * DIM_MAT)] = 1; // fFar / (fFar - fNear);
    mtx_projection.mat[2 + (3 * DIM_MAT)] = 1; // (-fFar * fNear) / (fFar - fNear);
    mtx_projection.mat[3 + (2 * DIM_MAT)] = 1; // 1.0f;
    mtx_projection.mat[3 + (3 * DIM_MAT)] = 0; // 0.0f;
    return &mtx_projection;
}

Mtx *proj_update_all(float near, float far, float fov_deg, float aspect_ratio) {
    float fov_rad = 1.0 / tanf(fov_deg * (M_PI / (180.0 * 2)));
    mtx_projection.mat[0 + (0 * DIM_MAT)] = aspect_ratio * fov_rad;
    mtx_projection.mat[1 + (1 * DIM_MAT)] = fov_rad;
    mtx_projection.mat[2 + (2 * DIM_MAT)] = far / (far - near);
    mtx_projection.mat[2 + (3 * DIM_MAT)] = (-far * near) / (far - near);
    mtx_projection.mat[3 + (2 * DIM_MAT)] = 1.0f;
    mtx_projection.mat[3 + (3 * DIM_MAT)] = 0.0f;
    return &mtx_projection;
}

Mtx *proj_update_fov(float fov, float aspect_ratio) {
    float fov_rad = 1.0 / tanf(fov * (M_PI / (180.0 * 2)));
    mtx_projection.mat[0 + (0 * DIM_MAT)] = aspect_ratio * fov_rad;
    mtx_projection.mat[1 + (1 * DIM_MAT)] = fov_rad;
    return &mtx_projection;
}

Mtx *proj_update_nearfar(float near, float far) {
    mtx_projection.mat[2 + (2 * DIM_MAT)] = far / (far - near);
    mtx_projection.mat[2 + (3 * DIM_MAT)] = (-far * near) / (far - near);
    return &mtx_projection;
}

// display function

void Pixel_init(Pixel *pixel, float *vect) {
    if (vect == NULL) {
        Vect_init(&pixel->P, 4, NULL);
        pixel->r = 0.0;
        pixel->g = 0.0;
        pixel->b = 0.0;
        pixel->a = 0.0;

    } else {
        Vect_init(&pixel->P, 4, vect);
    }
}

void Vect3d_init(Vect3d *vect3d, float *vect) {
    if (vect == NULL) {
        Vect_init(&vect3d->V, 3, NULL);
        vect3d->x = 0.0;
        vect3d->y = 0.0;
        vect3d->z = 0.0;

    } else {
        Vect_init(&vect3d->V, 3, vect);
    }
}

void Vect2d_init(Vect2d *vect2d, float *vect) {
    if (vect == NULL) {
        Vect_init(&vect2d->V, 2, NULL);
        vect2d->x = 0.0;
        vect2d->y = 0.0;

    } else {
        Vect_init(&vect2d->V, 2, vect);
    }
}

Pixel canvas[Window_Dim_X * Window_Dim_Y];

void init_canvas(Pixel *canvas) {
    for (int i = 0; i < Window_Dim_X * Window_Dim_Y; i++) {
        Vect_init(&(&canvas[i])->P, 4, NULL);
    }
}

void display_line(Pixel *canvas, Vect2d ps1, Vect2d ps2, Pixel color) {
    ps1.x = (int)ps1.x;
    ps1.y = (int)ps1.y;
    ps2.x = (int)ps2.x;
    ps2.y = (int)ps2.y;

    int dx = fabsf(ps2.x - ps1.x), dy = fabsf(ps2.y - ps1.y);
    int sx = (ps1.x < ps2.x) ? 1 : -1;
    int sy = (ps1.y < ps2.y) ? 1 : -1;
    int err = dx - dy;

    while (1) {
        if (ps1.x >= 0 && ps1.x < Window_Dim_X && ps1.y >= 0 && ps1.y < Window_Dim_Y)
            canvas[(int)ps1.y * Window_Dim_X + (int)ps1.x] = color;
        if (ps1.x == ps2.x && ps1.y == ps2.y)
            break;
        int e2 = err * 2;
        if (e2 > -dy) {
            err -= dy;
            ps1.x += sx;
        }
        if (e2 < dx) {
            err += dx;
            ps1.y += sy;
        }
    }
}

// TODO: awfull
void display_line_bezier(Pixel *canvas, Vect2d P0, Vect2d P1, Vect2d P2, Vect2d P3,
                         Pixel c0, Pixel c1) {
    for (float t = 0; t <= 1; t += 0.01) {
        float u = 1 - t;
        float x = u * u * u * P0.V.v[0] + 3 * u * u * t * P1.V.v[0] +
                  3 * u * t * t * P2.V.v[0] + t * t * t * P3.V.v[0];
        float y = u * u * u * P0.V.v[1] + 3 * u * u * t * P1.V.v[1] +
                  3 * u * t * t * P2.V.v[1] + t * t * t * P3.V.v[1];

        // Interpolate color
        Pixel color = {.r = (uint8_t)(1 - t) * c0.r + t * c1.r,
                       .g = (uint8_t)(1 - t) * c0.g + t * c1.g,
                       .b = (uint8_t)(1 - t) * c0.b + t * c1.b,
                       .a = 255};
        canvas[(int)y * Window_Dim_X + (int)x] = color;
    }
}

/*
Vect display_line_AA(){
}

Vect display_triangle(){
}

Vect display_mesh(){

}

Vect display_light(){
}
*/
