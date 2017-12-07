// vim: set ft=c:

struct MBD_calc* mbd_init_calc(int n_grid);

void mbd_destroy_calc(struct MBD_calc* calc);

struct MBD_damping* mbd_init_damping(
    int n_atoms,
    double* R_vdw,
    double beta,
    double a
);

void mbd_destroy_damping(struct MBD_damping* damping);

void calc_mbd_energy(
    struct MBD_calc* calc,
    int n_atoms,
    double* coords,
    double* alpha_0,
    double* omega,
    struct MBD_damping* damping,
    double* energy
);

void calc_mbd_rsscs_energy(
    struct MBD_calc* calc,
    int n_atoms,
    double* coords,
    double* alpha_0,
    double* omega,
    struct MBD_damping* damping,
    double* energy
);
