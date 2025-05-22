#include "the_header_file.h"



// matrix of rotation 

// X axe
Mtx *Mtx_rot_X(Mtx *m, float angle){ // angle in degs
    if (m == NULL){
        fprintf(stderr,"there is no matrix for Mtx_rot_X");
        return NULL;
    }
    m->dimC = 4;
    m->dimR = 4;
    float cosA = cosf((angle / 180.0) * M_PI);
    float sinA = sinf((angle / 180.0) * M_PI);
    float mtxRot[16] = {1,    0,     0, 0,
                        0, cosA, -sinA, 0,
                        0, sinA,  cosA, 0,
                        0,    0,     0, 1 };
    Mtx_copy_list(m, mtxRot);
    return m;
}

// Y axe
Mtx *Mtx_rot_Y(Mtx *m, float angle){ // angle in degs
    if (m == NULL){
        fprintf(stderr,"there is no matrix for Mtx_rot_X");
        return NULL;
    }
    m->dimC = 4;
    m->dimR = 4;
    float cosA = cosf((angle / 180.0) * M_PI);
    float sinA = sinf((angle / 180.0) * M_PI);
    float mtxRot[16] = { cosA, 0 , sinA, 0,
                           0 , 1 ,   0 , 0,
                        -sinA, 0 , cosA, 0,
                           0 , 0 ,   0 , 1 };
    Mtx_copy_list(m, mtxRot);
    return m;
}

// Z axe
Mtx *Mtx_rot_Z(Mtx *m, float angle){ // angle in degs
    if (m == NULL){
        fprintf(stderr,"there is no matrix for Mtx_rot_X");
        return NULL;
    }
    m->dimC = 4;
    m->dimR = 4;
    float cosA = cosf((angle / 180.0) * M_PI);
    float sinA = sinf((angle / 180.0) * M_PI);
    float mtxRot[16] = {cosA, -sinA, 0 , 0,
                        sinA, cosA , 0 , 0,
                         0  ,   0  , 1 , 0,
                         0  ,   0  , 0 , 1 };
    Mtx_copy_list(m, mtxRot);
    return m;
}


Vect Triangle_center(Triangle *triangle){
    Vect center;
    Vect_init(&center, 3, NULL);
    
    Vect_add(&center, &triangle->p[0].V, NULL);
    Vect_add(&center, &triangle->p[1].V, NULL);
    Vect_add(&center, &triangle->p[2].V, NULL);
    
    Vect_scal(&center, 1.0/3.0, NULL);
    return center;
}

Vect Mesh_center(Mesh *mesh){
    Vect center;
    Vect_init(&center, 3, NULL);

    for (size_t i = 0; i < mesh->numtris; i ++){
        for (size_t j = 0; j < mesh->numtris; j++){
            Vect_add(&center, &mesh->tris[i].p[j].V, NULL);
        }
    }
    Vect_scal(&center, 1.0/ (mesh->numtris * 3), NULL);

    return center;
}


// translation :  Mesh -> triangle -> vect


Vect *Vect_tr(Vect *vect, float dx, float dy, float dz){
    vect->v[0] += dx;
    vect->v[1] += dy;
    vect->v[2] += dz;
    return vect;
}

Triangle *Triangle_tr(Triangle *triangle, float dx, float dy, float dz){
    // only the p move for now
    Vect_tr(&triangle->p[0].V, dx, dy, dz);
    Vect_tr(&triangle->p[1].V, dx, dy, dz);
    Vect_tr(&triangle->p[2].V, dx, dy, dz);
    return triangle;
}

Mesh *Mesh_tr(Mesh *mesh, float dx, float dy, float dz){
    for (size_t i = 0; i < mesh->numtris ; i++){
        Triangle_tr(&mesh->tris[i], dx, dy , dz);
    }
    return mesh;
}


// scaling : mesh -> triangle 

Triangle *Triangle_scal(Triangle *triangle, Vect *center, float ratio){
    Vect temp;
    Vect_init(&temp, 3, NULL);
    
    if (center == NULL){
        Vect center_triangle;
        Vect_init(&center_triangle, 3, NULL);
        center_triangle = Triangle_center(triangle);
       
        Vect_add(&center_triangle,Vect_scal(Vect_sub(&triangle->p[0].V, &center_triangle, &temp), ratio, NULL),&triangle->p[0].V);
        Vect_add(&center_triangle,Vect_scal(Vect_sub(&triangle->p[1].V, &center_triangle, &temp), ratio, NULL),&triangle->p[1].V);
        Vect_add(&center_triangle,Vect_scal(Vect_sub(&triangle->p[2].V, &center_triangle, &temp), ratio, NULL),&triangle->p[2].V);
    }else{
        Vect_add(center,Vect_scal(Vect_sub(&triangle->p[0].V, center, &temp), ratio, NULL),&triangle->p[0].V);
        Vect_add(center,Vect_scal(Vect_sub(&triangle->p[1].V, center, &temp), ratio, NULL),&triangle->p[1].V);
        Vect_add(center,Vect_scal(Vect_sub(&triangle->p[2].V, center, &temp), ratio, NULL),&triangle->p[2].V);

    }
    return triangle;
}

Mesh *Mesh_scal(Mesh * mesh, float ratio){
    Vect center;
    Vect_init(&center, 3, NULL);
    center = Mesh_center(mesh);

    for( size_t i = 0; i < mesh->numtris; i++){
        Triangle_scal(&mesh->tris[i], &center, ratio);
    }
    return mesh;
}

// rotation mesh -> triangle -> vect   QUATERNION (with pivot)


void Quaternion_init(Quaternion *quaternion, float *vect){
    if (vect == NULL){
        Vect_init(&quaternion->Q, 4, NULL);
        quaternion->w = 0.0;
        quaternion->x = 0.0;
        quaternion->y = 0.0;
        quaternion->z = 0.0;

    }else{
        Vect_init(&quaternion->Q, 4, vect);
    }
}

// From Axis and Angle(deg)
Quaternion *Quaternion_FAA(Quaternion *quaternion, Vect *axis, float angle){
    angle *= (M_PI / 180.0)/2.0;// covnert to rad and divide by 2
    float sinA = sinf(angle);
    quaternion->w = cosf(angle);
    quaternion->x = axis->v[0] * sinA; 
    quaternion->y = axis->v[1] * sinA; 
    quaternion->z = axis->v[2] * sinA; 

    return quaternion;
}

Quaternion *Quaternion_inv(Quaternion *quaternion, Quaternion *q_result){
    if(q_result == NULL){
        quaternion->x *= -1;
        quaternion->y *= -1;
        quaternion->z *= -1;
        
        return quaternion;
    }else{
        q_result->w = quaternion->w;
        q_result->x = quaternion->x * -1;
        q_result->y = quaternion->y * -1;
        q_result->z = quaternion->z * -1;

        return q_result;
    }
}

Quaternion* Quaternion_mult(Quaternion *q1, Quaternion *q2, Quaternion *qr){
    if(qr == NULL || qr == q2){
        Quaternion q_temp;
        Quaternion_init(&q_temp, NULL);
        q_temp.w = q1->w * q2->w - q1->x * q2->x - q1->y * q2->y - q1->z * q2->z;
        q_temp.x = q1->w * q2->x + q1->x * q2->w + q1->y * q2->z - q1->z * q2->y;
        q_temp.y = q1->w * q2->y - q1->x * q2->z + q1->y * q2->w + q1->z * q2->x;
        q_temp.z = q1->w * q2->z + q1->x * q2->y - q1->y * q2->x + q1->z * q2->w;
        
        q1->w = q_temp.w;
        q1->x = q_temp.x;
        q1->y = q_temp.y;
        q1->z = q_temp.z;
        return q1;
    }else{
        qr->w = q1->w * q2->w - q1->x * q2->x - q1->y * q2->y - q1->z * q2->z;
        qr->x = q1->w * q2->x + q1->x * q2->w + q1->y * q2->z - q1->z * q2->y;
        qr->y = q1->w * q2->y - q1->x * q2->z + q1->y * q2->w + q1->z * q2->x;
        qr->z = q1->w * q2->z + q1->x * q2->y - q1->y * q2->x + q1->z * q2->w;
        
        return qr;
    }
}

Vect *Vect_rot(Vect *vect, Vect *pivot, Quaternion *quaternion, Vect* vect_r){
    Quaternion vect_q;
    float temp[4] = {0, vect->v[0] - pivot->v[0], vect->v[1] - pivot->v[1], vect->v[2] - pivot->v[2]}; 
    Quaternion_init(&vect_q, temp);

    Quaternion quaternion_inv;
    Quaternion_init(&quaternion_inv, NULL);
    Quaternion_inv(quaternion, &quaternion_inv);

    Quaternion_mult(&quaternion_inv, &vect_q, NULL);
    Quaternion_mult(&quaternion_inv, quaternion, &vect_q);

    if(vect_r == NULL){
        vect->v[0] = vect_q.Q.v[1] + pivot->v[0];
        vect->v[1] = vect_q.Q.v[2] + pivot->v[1];
        vect->v[2] = vect_q.Q.v[3] + pivot->v[2];
        return vect;
    }else{
        vect_r->v[0] = vect_q.Q.v[1] + pivot->v[0];
        vect_r->v[1] = vect_q.Q.v[2] + pivot->v[1];
        vect_r->v[2] = vect_q.Q.v[3] + pivot->v[2];
        return vect_r;
    }
}

Triangle *Triangle_rot(Triangle *triangle, Vect *pivot, Quaternion *quaternion, Triangle *triangle_r){
    if (triangle_r == NULL){
        Vect_rot(&triangle->p[0].V, pivot, quaternion, NULL);
        Vect_rot(&triangle->p[1].V, pivot, quaternion, NULL);
        Vect_rot(&triangle->p[2].V, pivot, quaternion, NULL);
        return triangle;
    }else{
        Vect_rot(&triangle->p[0].V, pivot, quaternion, &triangle_r->p[0].V);
        Vect_rot(&triangle->p[1].V, pivot, quaternion, &triangle_r->p[1].V);
        Vect_rot(&triangle->p[2].V, pivot, quaternion, &triangle_r->p[2].V);
        return triangle_r;
    }
}

Mesh *Mesh_rot(Mesh *mesh, Vect *pivot, Quaternion *quaternion, Mesh *mesh_r){
    if (mesh == NULL){
        for(size_t i = 0; i < mesh->numtris; i++){
            Triangle_rot(&mesh->tris[i], pivot, quaternion, NULL);
        }
        return mesh;
    }else{
        for(size_t i = 0; i < mesh->numtris; i++){
            Triangle_rot(&mesh->tris[i], pivot, quaternion, &mesh_r->tris[i]);
        }
        return mesh_r;
    }
}

// shearing : mesh -> triangle -> vect

//not use though
Mtx *Mtx_shearing(Mtx *m, float sh_xy, float sh_xz, float sh_yx, float sh_yz, float sh_zx, float sh_zy){ //streching
    if (m == NULL){
        fprintf(stderr,"there is no matrix for Mtx_rot_X");
        return NULL;
    }
    m->dimC = 4;
    m->dimR = 4;
    float mtxRot[16] = {   1  , sh_xy, sh_xz,
                         sh_yx,   1  , sh_yz,
                         sh_zx, sh_zy,   1   };
                         
    Mtx_copy_list(m, mtxRot);
    return m;
}

Vect *Vect_shearing(Vect *vect, float sh_xy, float sh_xz, float sh_yx, float sh_yz, float sh_zx, float sh_zy, Vect *vect_r){
    float x = vect->v[0];
    float y = vect->v[1];
    float z = vect->v[2];
    if (vect_r == NULL){
        vect->v[0] = x + sh_xy * y + sh_xz * z;
        vect->v[1] = sh_yx * x + y + sh_yz * z;
        vect->v[2] = sh_zx * x + sh_zy * y + z;
        return vect;
    }else{
        vect_r->v[0] = x + sh_xy * y + sh_xz * z;
        vect_r->v[1] = sh_yx * x + y + sh_yz * z;
        vect_r->v[2] = sh_zx * x + sh_zy * y + z;
        return vect_r;
    }
}

Triangle *Triangle_shearing(Triangle * triangle, float sh_xy, float sh_xz, float sh_yx, float sh_yz, float sh_zx, float sh_zy, Triangle *triangle_r){
    if(triangle_r == NULL){
        Vect_shearing(&triangle->p[0].V, sh_xy, sh_xz, sh_yx, sh_yz, sh_zx, sh_zy, NULL);
        Vect_shearing(&triangle->p[1].V, sh_xy, sh_xz, sh_yx, sh_yz, sh_zx, sh_zy, NULL);
        Vect_shearing(&triangle->p[2].V, sh_xy, sh_xz, sh_yx, sh_yz, sh_zx, sh_zy, NULL);
        return triangle;
    }else{
        Vect_shearing(&triangle->p[0].V, sh_xy, sh_xz, sh_yx, sh_yz, sh_zx, sh_zy, &triangle_r->p[0].V);
        Vect_shearing(&triangle->p[1].V, sh_xy, sh_xz, sh_yx, sh_yz, sh_zx, sh_zy, &triangle_r->p[0].V);
        Vect_shearing(&triangle->p[2].V, sh_xy, sh_xz, sh_yx, sh_yz, sh_zx, sh_zy, &triangle_r->p[0].V);
    return triangle_r;
    }
}
