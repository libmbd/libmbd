// vim: set ft=c:

struct MBD_calc {
    int n_freq;
    double* omega_grid;
    double* omega_grid_w;
};

struct MBD_calc* mbd_init_calc(void);
void mbd_set_parallel(struct MBD_calc* calc, int rank, int n_proc);
void mbd_destroy_calc(struct MBD_calc* calc);

struct MBD_system {
    double* forces;
    _Bool* do_force;
};

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

double calc_mbd_energy(
    struct MBD_system* sys,
    int n_atoms,
    double* alpha_0,
    double* C6,
    struct MBD_damping* damping
);

double calc_rpa_energy(
    struct MBD_system* sys,
    int n_atoms,
    double* alpha_0,
    double* C6,
    struct MBD_damping* damping
);

double calc_mbd_rsscs_energy(
    struct MBD_system* sys,
    int n_atoms,
    double* alpha_0,
    double* C6,
    struct MBD_damping* damping
);

double calc_mbd_scs_energy(
    struct MBD_system* sys,
    int n_atoms,
    double* alpha_0,
    double* C6,
    struct MBD_damping* damping
);

double calc_dipole_matrix(
    struct MBD_system* sys,
    struct MBD_damping* damping,
    double* k_point,
    double* dipmat
);
