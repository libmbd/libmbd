// vim: set ft=c:

extern const _Bool with_scalapack;

struct MBD_calc {
    int n_freq;
    double* omega_grid;
    double* omega_grid_w;
};

struct MBD_calc* mbd_init_calc(int n_freq);
void mbd_set_parallel(struct MBD_calc* calc);
void mbd_destroy_calc(struct MBD_calc* calc);

struct MBD_system* mbd_init_system(
    struct MBD_calc* calc,
    int n_atoms,
    double* coords,
    double* lattice,
    int* k_grid
);

void mbd_destroy_system(struct MBD_system* sys);

struct MBD_damping* mbd_init_damping(
    int n_atoms,
    char* version,
    double* R_vdw,
    double* sigma,
    double beta,
    double a
);

void mbd_destroy_damping(struct MBD_damping* damping);

double calc_ts_energy(
    struct MBD_system* sys,
    int n_atoms,
    double* alpha_0,
    double* C6,
    struct MBD_damping* damping,
    double* gradients
);

double calc_mbd_energy(
    struct MBD_system* sys,
    int n_atoms,
    double* alpha_0,
    double* C6,
    struct MBD_damping* damping,
    double* gradients
);

double calc_rpa_energy(
    struct MBD_system* sys,
    int n_atoms,
    double* alpha_0,
    double* C6,
    struct MBD_damping* damping,
    double* gradients
);

double calc_mbd_rsscs_energy(
    struct MBD_system* sys,
    int n_atoms,
    double* alpha_0,
    double* C6,
    struct MBD_damping* damping,
    double* gradients,
    double* eigvals,
    double* eigvecs
);

double calc_mbd_scs_energy(
    struct MBD_system* sys,
    int n_atoms,
    double* alpha_0,
    double* C6,
    struct MBD_damping* damping,
    double* gradients
);

double calc_dipole_matrix(
    struct MBD_system* sys,
    struct MBD_damping* damping,
    double* k_point,
    double* dipmat
);

void calc_full_coulomb(
    int n,
    double* coords,
    double* C,
    double* w,
    double* w0,
    double* a0,
    double* rvdw0,
    double alpha,
    double beta,
    char* version,
    double dampswitch,
    double* ecoul,
    double* en,
    double* ee,
    double* nn
);

double calc_coulomb_energy(
    struct MBD_system* sys,
    int n_atoms,
    double* q,
    double* m,
    double* w_t,
    double* C
);

double calc_get_dipole_energy(
    struct MBD_calc* calc,
    int n,
    char* version,
    double* R,
    double* a0,
    double* w,
    double* w_t,
    double* r0,
    double beta,
    double alpha,
    double* C
);
