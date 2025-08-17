#include "../headers/general.h"
#include "../headers/sequential.h"
#include <assert.h>
#include <stddef.h>
#include <stdint.h>

// TODO: 3. Générateur d'opérateurs aux différences finies (FD)  <-- RÉFÉRENCE AJOUTÉE
//   But: construire automatiquement la matrice discrète pour un opérateur linéaire
//   spatial donné API proposée:
//     matrix *build_fd_operator_1D(uint64_t nx, double dx,
//                                 operator_descriptor *ops, size_t nops,
//                                 bc_descriptor *bcs);
//     matrix *build_fd_operator_2D(uint64_t nx, uint64_t ny, double dx, double dy,
//                                 operator_descriptor *ops, size_t nops,
//                                 bc_descriptor *bcs);
//   where operator_descriptor contains: { type: LAPLACIAN | GRADIENT_X | ADVECTION_X |
//   CUSTOM, coefficient (complex), stencil_order } and bc_descriptor describes: {
//   boundary_type: DIRICHLET | NEUMANN | PERIODIC, value or function, node indices }
//   Sous-tâches:
//     - 4.1 Définir structures:
//         typedef struct { int type; complex coeff; int order; } operator_descriptor;
//         typedef struct { int type; /* dirichlet/neumann/periodic */ double
//         (*val)(double x, double y); } bc_descriptor;
//     - 4.2 Implémenter generate_laplacian_1D(nx, dx) réutilisable par
//     build_fd_operator_1D
//     - 4.3 Implémenter operator appliqué: add_operator_to_matrix(matrix *M,
//     operator_descriptor *op, position info)
//     - 4.4 Appliquer BCs: apply_dirichlet(M, rhs, node, value), apply_neumann(...) — API
//     commune
//     - 4.5 Assembler matrice globale creuse (ou dense si prototype): support row-major
//     indexing
//     - 4.6 Tests unitaires:
//         - Poisson 1D avec Dirichlet: comparer solution analytique
//         - Diffusion explicite/implicite : stability check
//     - 4.7 Retourner matrice et vecteur source (si opérateur inclut termes constants)
//   Notes:
//     - Cette approche couvre la majorité des PDE linéaires (diffusion, advection
//     linéaire, Poisson, Helmholtz).
//     - Pour des schémas upwind ou d’ordre >2, ajouter des stencils personnalisés dans
//     operator_descriptor.

#define IDX(i, j, nx) ((i) + (j) * (nx))

// --- Operators you can combine linearly
typedef enum {
    OP_I,   // identity (reaction)
    OP_DX,  // first derivative x
    OP_DY,  // first derivative y
    OP_DXX, // second derivative x
    OP_DYY, // second derivative y
    OP_DXY  // mixed second derivative
} op_type;

typedef enum { SCHEME_CENTERED, SCHEME_UPWIND } scheme_type;

typedef struct {
    op_type type;
    complex coeff;      // scale for this operator
    int order;          // 1 or 2 supported here (accuracy order of scheme)
    scheme_type scheme; // upwind only used for DX/DY
    // optional parameters if needed later (e.g., wind for upwind direction):
    double ax, ay; // advection direction (sign decides upwind)
} operator_descriptor;

typedef enum { BC_DIRICHLET, BC_NEUMANN, BC_PERIODIC } bc_type;

typedef struct {
    bc_type left, right, bottom, top;
    // Value callbacks (can be constant-returning functions) for Dirichlet/Neumann
    complex (*g_left)(double x, double y);
    complex (*g_right)(double x, double y);
    complex (*g_bottom)(double x, double y);
    complex (*g_top)(double x, double y);
} bc_descriptor;

typedef struct {
    uint64_t nx, ny;
    double dx, dy;
    double x0, y0; // physical origin (optional, defaults 0,0)
} grid2d;

static inline double x_at(const grid2d *G, uint64_t i) { return G->x0 + i * G->dx; }
static inline double y_at(const grid2d *G, uint64_t j) { return G->y0 + j * G->dy; }

// rhs_vec = vec
static inline void m_add_ee(matrix *A, uint64_t r, uint64_t c, complex v) {
    A->data[r * A->n + c] = c_add(A->data[r * A->n + c], v);
}
static void add_entry_with_bc(matrix *A, vector *rhs, const grid2d *G,
                              const bc_descriptor *BC, uint64_t i, uint64_t j, int64_t di,
                              int64_t dj, complex w) {
    int64_t ii = (int64_t)i + di;
    int64_t jj = (int64_t)j + dj;
    // check noundaries
    if (ii >= 0 && ii < (int64_t)G->nx && jj >= 0 && jj < (int64_t)G->ny) {
        m_add_ee(A, IDX(i, j, G->nx), IDX((uint64_t)ii, (uint64_t)jj, G->nx), w);
        return;
    }
    // did we hit boundary ?
    // check left/right
    if (ii < 0) {
        if (BC->left == BC_PERIODIC) {
            int64_t wrap_i = (ii % (int64_t)G->nx + (int64_t)G->nx) % (int64_t)G->nx;
            m_add_ee(A, IDX(i, j, G->nx), IDX((uint64_t)wrap_i, (uint64_t)jj, G->nx), w);
            return;
        } else if (BC->left == BC_DIRICHLET) {
            complex val = BC->g_left ? BC->g_left(x_at(G, 0), y_at(G, j)) : c_zero();
            rhs->data[IDX(i, j, G->nx)] =
                c_sub(rhs->data[IDX(i, j, G->nx)], c_mult(w, val));
            return;
        } else if (BC->left == BC_NEUMANN) {
            complex q = BC->g_left ? BC->g_left(x_at(G, 0), y_at(G, j)) : c_zero();
            if (di == -1 && dj == 0) {
                if (i + 1 < G->nx) {
                    m_add_ee(A, IDX(i, j, G->nx), IDX(i + 1, j, G->nx), w);
                }
                rhs->data[IDX(i, j, G->nx)] =
                    c_sub(rhs->data[IDX(i, j, G->nx)],
                          c_mult(c_mult(w, q), c_real(2.0 * G->dx)));
                return;
            }
        } else {
            assert(0 && "Unhandled BC type at left boundary");
        }
    }
    if (ii >= (int64_t)G->nx) {
        if (BC->right == BC_PERIODIC) {
            int64_t wrap_i = ii % (int64_t)G->nx;
            m_add_ee(A, IDX(i, j, G->nx), IDX((uint64_t)wrap_i, (uint64_t)jj, G->nx), w);
            return;
        } else if (BC->right == BC_DIRICHLET) {
            complex val =
                BC->g_right ? BC->g_right(x_at(G, G->nx - 1), y_at(G, j)) : c_zero();
            rhs->data[IDX(i, j, G->nx)] =
                c_sub(rhs->data[IDX(i, j, G->nx)], c_mult(w, val));
            return;
        } else if (BC->right == BC_NEUMANN) {
            complex q =
                BC->g_right ? BC->g_right(x_at(G, G->nx - 1), y_at(G, j)) : c_zero();
            if (di == +1 && dj == 0) {
                if (i >= 1) {
                    m_add_ee(A, IDX(i, j, G->nx), IDX(i - 1, j, G->nx), w);
                }
                rhs->data[IDX(i, j, G->nx)] =
                    c_add(rhs->data[IDX(i, j, G->nx)],
                          c_mult(c_mult(w, q), c_real(2.0 * G->dx)));
                return;
            }
        } else {
            assert(0 && "Unhandled BC type at right boundary");
        }
    }
    // Check Bottom/Top
    if (jj < 0) {
        if (BC->bottom == BC_PERIODIC) {
            int64_t wrap_j = (jj % (int64_t)G->ny + (int64_t)G->ny) % (int64_t)G->ny;
            m_add_ee(A, IDX(i, j, G->nx), IDX((uint64_t)ii, (uint64_t)wrap_j, G->nx), w);
            return;
        } else if (BC->bottom == BC_DIRICHLET) {
            complex val = BC->g_bottom ? BC->g_bottom(x_at(G, i), y_at(G, 0)) : c_zero();
            rhs->data[IDX(i, j, G->nx)] =
                c_sub(rhs->data[IDX(i, j, G->nx)], c_mult(w, val));
            return;
        } else if (BC->bottom == BC_NEUMANN) {
            complex q = BC->g_bottom ? BC->g_bottom(x_at(G, i), y_at(G, 0)) : c_zero();
            if (di == 0 && dj == -1) {
                if (j + 1 < G->ny) {
                    m_add_ee(A, IDX(i, j, G->nx), IDX(i, j + 1, G->nx), w);
                }
                rhs->data[IDX(i, j, G->nx)] =
                    c_sub(rhs->data[IDX(i, j, G->nx)],
                          c_mult(c_mult(w, q), c_real(2.0 * G->dy)));
                return;
            }
        } else {
            assert(0 && "Unhandled BC type at bottom boundary");
        }
    }
    if (jj >= (int64_t)G->ny) {
        if (BC->top == BC_PERIODIC) {
            int64_t wrap_j = jj % (int64_t)G->ny;
            m_add_ee(A, IDX(i, j, G->nx), IDX((uint64_t)ii, (uint64_t)wrap_j, G->nx), w);
            return;
        } else if (BC->top == BC_DIRICHLET) {
            complex val =
                BC->g_top ? BC->g_top(x_at(G, i), y_at(G, G->ny - 1)) : c_zero();
            rhs->data[IDX(i, j, G->nx)] =
                c_sub(rhs->data[IDX(i, j, G->nx)], c_mult(w, val));
            return;
        } else if (BC->top == BC_NEUMANN) {
            complex q = BC->g_top ? BC->g_top(x_at(G, i), y_at(G, G->ny - 1)) : c_zero();
            if (di == 0 && dj == +1) {
                if (j >= 1) {
                    m_add_ee(A, IDX(i, j, G->nx), IDX(i, j - 1, G->nx), w);
                }
                rhs->data[IDX(i, j, G->nx)] =
                    c_add(rhs->data[IDX(i, j, G->nx)],
                          c_mult(c_mult(w, q), c_real(2.0 * G->dy)));
                return;
            }
        } else {
            assert(0 && "Unhandled BC type at top boundary");
        }
    }
}

typedef struct {
    int n;
    int di[9];
    int dj[9];
    double w[9];
} stencil;

static stencil S_identity(void) {
    stencil s = {.n = 1, .di = {0}, .dj = {0}, .w = {1.0}};
    return s;
}

// First derivatives
static stencil S_dx(double dx, scheme_type scheme, int order, double ax) {
    stencil s = {0};
    if (scheme == SCHEME_UPWIND && ax != 0.0) {
        // first-order upwind
        if (ax > 0) {
            s = (stencil){
                .n = 2, .di = {0, -1}, .dj = {0, 0}, .w = {1.0 / dx, -1.0 / dx}};
        } else {
            s = (stencil){
                .n = 2, .di = {+1, 0}, .dj = {0, 0}, .w = {1.0 / dx, -1.0 / dx}};
        }
    } else {
        // centered
        if (order == 2) {
            s = (stencil){
                .n = 2, .di = {+1, -1}, .dj = {0, 0}, .w = {+0.5 / dx, -0.5 / dx}};
        }
        if (order == 5) {
            s = (stencil){.n = 5,
                          .di = {+2, +1, 0, -1, -2},
                          .dj = {0, 0, 0, 0, 0},
                          .w = {-1.0 / (12 * dx), +8.0 / (12 * dx), 0.0, -8.0 / (12 * dx),
                                +1.0 / (12 * dx)}};
        } else { // fallback to 2^rd
            s = (stencil){
                .n = 2, .di = {+1, -1}, .dj = {0, 0}, .w = {+0.5 / dx, -0.5 / dx}};
        }
    }
    return s;
}
static stencil S_dy(double dy, scheme_type scheme, int order, double ay) {
    stencil s = {0};
    if (scheme == SCHEME_UPWIND && ay != 0.0) {
        if (ay > 0) {
            s = (stencil){
                .n = 2, .di = {0, 0}, .dj = {0, -1}, .w = {1.0 / dy, -1.0 / dy}};
        } else {
            s = (stencil){
                .n = 2, .di = {0, 0}, .dj = {+1, 0}, .w = {1.0 / dy, -1.0 / dy}};
        }
    } else {
        if (order == 2) {
            s = (stencil){
                .n = 2, .di = {0, 0}, .dj = {+1, -1}, .w = {+0.5 / dy, -0.5 / dy}};
        } else if (order == 5) {
            s = (stencil){.n = 5,
                          .di = {0, 0, 0, 0, 0},
                          .dj = {+2, +1, 0, -1, -2},
                          .w = {-1.0 / (12 * dy), +8.0 / (12 * dy), 0.0, -8.0 / (12 * dy),
                                +1.0 / (12 * dy)}};
        } else { // fallback to 2^rd
            s = (stencil){
                .n = 2, .di = {0, 0}, .dj = {+1, -1}, .w = {+0.5 / dy, -0.5 / dy}};
        }
    }
    return s;
}

// Second derivatives (2nd order centered)
static stencil S_dxx(double dx) {
    double c = 1.0 / (dx * dx);
    stencil s = {.n = 3, .di = {-1, 0, +1}, .dj = {0, 0, 0}, .w = {c, -2.0 * c, c}};
    return s;
}
static stencil S_dyy(double dy) {
    double c = 1.0 / (dy * dy);
    stencil s = {.n = 3, .di = {0, 0, 0}, .dj = {-1, 0, +1}, .w = {c, -2.0 * c, c}};
    return s;
}

// Mixed derivative d2/(dx dy) using 4-point cross (second-order accurate)
static stencil S_dxy(double dx, double dy) {
    double c = 1.0 / (4.0 * dx * dy);
    stencil s = {
        .n = 4, .di = {+1, +1, -1, -1}, .dj = {+1, -1, +1, -1}, .w = {+c, -c, -c, +c}};
    return s;
}

// ---- get stencil for an operator -------------------------------------------
static inline stencil stencil_for(const operator_descriptor *op, const grid2d *G) {
    switch (op->type) {
    case OP_I:
        return S_identity();
    case OP_DX:
        return S_dx(G->dx, op->scheme, op->order, op->ax);
    case OP_DY:
        return S_dy(G->dy, op->scheme, op->order, op->ay);
    case OP_DXX:
        return S_dxx(G->dx);
    case OP_DYY:
        return S_dyy(G->dy);
    case OP_DXY:
        return S_dxy(G->dx, G->dy);
    default:
        return S_identity(); // safe fallback
    }
}

// ---- apply ONE operator to ONE row (node i,j) -------------------------------
static inline void apply_operator_to_row(matrix *A, vector *rhs, const grid2d *G,
                                         const bc_descriptor *BC, uint64_t i, uint64_t j,
                                         const operator_descriptor *op) {
    stencil s = stencil_for(op, G);

    for (int t = 0; t < s.n; ++t) {
        // Real stencil weight -> complex, then scale by operator's complex coeff
        complex w_real = c_real(s.w[t]);            // (s.w[t] is double)
        complex w_full = c_mult(op->coeff, w_real); // complex weight

        // NOTE: di/dj MUST be signed in add_entry_with_bc! (see signature below)
        add_entry_with_bc(A, rhs, G, BC, i, j, (int)s.di[t], (int)s.dj[t], w_full);
    }
}

// ---- assemble all operators over the whole grid -----------------------------
static inline void add_all_ops(matrix *A, vector *rhs, const grid2d *G,
                               const bc_descriptor *BC, const operator_descriptor *ops,
                               size_t nops) {
    for (uint64_t j = 0; j < G->ny; ++j) {
        for (uint64_t i = 0; i < G->nx; ++i) {
            for (size_t k = 0; k < nops; ++k) {
                apply_operator_to_row(A, rhs, G, BC, i, j, &ops[k]);
            }
        }
    }
}

typedef struct {
    matrix *A;
    vector *b;
} system_lin;

static system_lin build_fd_operator_2D(const grid2d *G, const operator_descriptor *ops,
                                       size_t nops, const bc_descriptor *BC) {
    uint64_t N = G->nx * G->ny;
    system_lin sys = {.A = m_init(N, N), .b = v_init(N)};
    vector *rhs = sys.b;
    for (uint64_t j = 0; j < G->ny; ++j) {
        for (uint64_t i = 0; i < G->nx; ++i) {
            for (size_t k = 0; k < nops; ++k) {
                const operator_descriptor *op = &ops[k];
                stencil s = stencil_for(op, G);
                for (int t = 0; t < s.n; ++t) {
                    complex w = c_mult(op->coeff, c_real(s.w[t]));
                    add_entry_with_bc(sys.A, rhs, G, BC, i, j, s.di[t], s.dj[t], w);
                }
            }
        }
    }
    return sys;
}
