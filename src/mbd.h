// vim: set ft=c:

struct MBD_calc* mbd_init_calc(void);

void mbd_destroy_calc(struct MBD_calc* calc);

struct MBD_system {
    double* forces;
    _Bool* do_force;
};

struct MBD_system* mbd_init_system(
    struct MBD_calc* calc,
    int n_atoms,
    double* coords,
    _Bool periodic,
    double* lattice,
    int* k_grid
);

void mbd_destroy_system(struct MBD_system* sys);

struct MBD_damping* mbd_init_damping(
    int n_atoms,
    double* R_vdw,
    double beta,
    double a
);

void mbd_destroy_damping(struct MBD_damping* damping);

double calc_mbd_energy(
    struct MBD_system* sys,
    int n_atoms,
    double* alpha_0,
    double* omega,
    struct MBD_damping* damping
);

double calc_rpa_energy(
    struct MBD_system* sys,
    int n_atoms,
    double* alpha_0,
    double* omega,
    struct MBD_damping* damping
);

double calc_mbd_rsscs_energy(
    struct MBD_system* sys,
    int n_atoms,
    double* alpha_0,
    double* omega,
    struct MBD_damping* damping
);

double calc_mbd_scs_energy(
    struct MBD_system* sys,
    int n_atoms,
    double* alpha_0,
    double* omega,
    struct MBD_damping* damping
);
