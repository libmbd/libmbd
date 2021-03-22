// vim: set ft=c:

extern const _Bool cmbd_with_mpi;
extern const _Bool cmbd_with_scalapack;
extern const int cmbd_version_major;
extern const int cmbd_version_minor;
extern const int cmbd_version_patch;
extern const char cmbd_version_suffix[30];

struct geom_t* cmbd_init_geom(
    int n_atoms,
    double* coords,
    double* lattice,
    int* k_grid,
    int n_kpts,
    double* custom_k_pts,
    int n_freq,
    _Bool do_rpa,
    _Bool get_spectrum,
    _Bool get_rpa_orders,
    _Bool rpa_rescale_eigs,
    int max_atoms_per_block
);

void cmbd_update_coords(struct geom_t* geom, double* coords);

void cmbd_update_lattice(struct geom_t* geom, double* lattice);

void cmbd_destroy_geom(struct geom_t* geom);

void cmbd_get_exception(
    struct geom_t* geom,
    int* code,
    char origin[50],
    char msg[150]
);

struct cmbd_damping* cmbd_init_damping(
    int n_atoms,
    char* version,
    double* R_vdw,
    double* sigma,
    double beta,
    double a
);

void cmbd_destroy_damping(struct cmbd_damping* damping);

void cmbd_print_timing(struct geom_t* geom);

struct result_t* cmbd_ts_energy(
    struct geom_t* geom,
    double* alpha_0,
    double* C6,
    struct cmbd_damping* damping,
    _Bool grad
);

struct result_t* cmbd_mbd_energy(
    struct geom_t* geom,
    double* alpha_0,
    double* C6,
    struct cmbd_damping* damping,
    _Bool grad
);

struct result_t* cmbd_mbd_scs_energy(
    struct geom_t* geom,
    char* variant,
    double* alpha_0,
    double* C6,
    struct cmbd_damping* damping,
    _Bool grad
);

void cmbd_get_results(
    struct result_t* result,
    double* energy,
    double* gradients,
    double* lattice_gradients,
    double* eigvals,
    double* eigvecs,
    double* rpa_orders,
    double* eigvals_k,  // is actually complex double
    double* eigvecs_k  // is actually complex double
);

void cmbd_destroy_result(struct result_t* result);

double cmbd_dipole_matrix(
    struct geom_t* geom,
    struct cmbd_damping* damping,
    double* q_point,
    double* dipmat
);

double cmbd_coulomb_energy(
    struct geom_t* geom,
    int n_atoms,
    double* q,
    double* m,
    double* w_t,
    char* version,
    double* r_vdw,
    double beta,
    double a,
    double* C
);

double cmbd_dipole_energy(
    struct geom_t* geom,
    int n_atoms,
    double* a0,
    double* w,
    double* w_t,
    char* version,
    double* r_vdw,
    double beta,
    double a,
    double* C
);
