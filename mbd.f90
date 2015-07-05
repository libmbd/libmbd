module mbd

use mbd_interface, only: &
    sync_sum_array, sync_sum_number, print_error, print_warning, print_log

implicit none

! these variables need to be initialized before the module is used
integer :: n_atoms
real*8, allocatable :: coords(:, :)
logical :: is_periodic
real*8 :: lattice_vector(3, 3)
real*8 :: electric_field(3)
integer :: n_tasks, my_task

! these variables can be modified but have defaults
real*8 :: &
    pi = acos(-1.d0), &
    nan = sqrt(-sin(0.d0))
real*8, parameter :: &
    bohr = 0.529177249
real*8 :: &
    param_scs_dip_cutoff = 120.d0/bohr, &
    param_mbd_supercell_cutoff = 25.d0/bohr, &
    param_mbd_dip_cutoff = 100.d0/bohr, &
    param_energy_accuracy = 1.d-6
integer :: &
    param_mbd_nbody_max = 3, &
    param_rpa_order_max = 10
logical :: &
    param_vacuum_axis(3) = (/ .false., .false., .false. /)
integer, parameter :: &
    n_grid_omega = 20
real*8 :: omega_grid(0:n_grid_omega) = (/ &
    0.0000000, &
    0.0392901, 0.1183580, 0.1989120, 0.2820290, 0.3689190, &
    0.4610060, 0.5600270, 0.6681790, 0.7883360, 0.9243900, &
    1.0817900, 1.2684900, 1.4966100, 1.7856300, 2.1691700, &
    2.7106200, 3.5457300, 5.0273400, 8.4489600, 25.451700 /)
real*8 :: omega_grid_w(0:n_grid_omega) = (/ &
    0.0000000, &
    0.0786611, 0.0796400, 0.0816475, 0.0847872, 0.0892294, &
    0.0952317, 0.1031720, 0.1136050, 0.1273500, 0.1456520, &
    0.1704530, 0.2049170, 0.2544560, 0.3289620, 0.4480920, &
    0.6556060, 1.0659600, 2.0635700, 5.6851000, 50.955800 /)

contains

function get_ts_energy( &
        C6, alpha_0, version, R_vdw, s_R, d, damping_custom, overlap) &
        result(ene)
    implicit none

    real*8, intent(in) :: C6(:), alpha_0(size(C6))
    character(len=*), intent(in) :: version
    real*8, intent(in), optional :: R_vdw(size(C6))
    real*8, intent(in), optional :: s_R, d
    real*8, intent(in), optional :: &
        damping_custom(size(C6), size(C6)), &
        overlap(size(C6), size(C6))
    real*8 :: ene

    real*8 :: C6_ij, r(3), r_norm, f_damp, R_vdw_sum, overlap_ij
    real*8 :: ene_add, ene_contrib
    real*8 :: cell_shift(3)
    integer :: n_shell, i_cell, j_cell, k_cell, ijk_cell(3)
    integer :: i_atom, j_atom

    ene = 0.d0
    n_shell = 0
    do
        ene_add = 0.d0
        do i_cell = -n_shell, n_shell
        do j_cell = -n_shell, n_shell
        do k_cell = -n_shell, n_shell
            ijk_cell = (/ i_cell, j_cell, k_cell /)
            if (all(abs(ijk_cell) /= n_shell)) cycle
            cell_shift = matmul(lattice_vector, ijk_cell)
            do i_atom = 1, n_atoms
                if (my_task /= modulo(i_atom, n_tasks)) cycle
                do j_atom = i_atom, n_atoms
                    if (n_shell == 0 .and. i_atom == j_atom) cycle
                    r = coords(:, i_atom)-coords(:, j_atom)-cell_shift
                    r_norm = sqrt(sum(r**2))
                    C6_ij = combine_C6(C6(i_atom), C6(j_atom), &
                                       alpha_0(i_atom), alpha_0(j_atom))
                    if (present(R_vdw)) then
                        R_vdw_sum = R_vdw(i_atom)+R_vdw(j_atom)
                    end if
                    if (present(overlap)) then
                        if (i_atom < j_atom) then
                            overlap_ij = overlap(i_atom, j_atom)
                        else
                            overlap_ij = overlap(j_atom, i_atom)
                        end if
                    end if
                    select case (version)
                        case ("fermi")
                            f_damp = damping_fermi(r_norm, s_R*R_vdw_sum, d)
                        case ("fermi2")
                            f_damp = damping_fermi(r_norm, s_R*R_vdw_sum, d)**2
                        case ("erf")
                            f_damp = damping_erf(r_norm, s_R*R_vdw_sum, d)
                        case ("1mexp")
                            f_damp = damping_1mexp(r_norm, s_R*R_vdw_sum, d)
                        case ("overlap")
                            f_damp = damping_overlap( &
                                r_norm, overlap_ij, C6_ij, s_R, d)
                        case ("custom")
                            f_damp = damping_custom(i_atom, j_atom)
                    end select
                    ene_contrib = -C6_ij/r_norm**6*f_damp
                    if (i_atom == j_atom) then
                        ene_contrib = ene_contrib/2
                    endif
                    ene_add = ene_add + ene_contrib
                enddo
            enddo
        enddo
        enddo
        enddo
        call sync_sum_number (ene_add)
        ene = ene + ene_add
        if (.not. is_periodic) exit
        if (n_shell > 0 .and. abs(ene_add) < param_energy_accuracy) then
            exit
        else
            n_shell = n_shell+1
        endif
    enddo
end function get_ts_energy

function run_scs(alpha, version, R_vdw, beta, a, damping_custom) & 
        result(relay)
    implicit none

    real*8, intent(in) :: alpha(:)
    character(len=*), intent(in) :: version
    real*8, intent(in), optional :: R_vdw(size(alpha))
    real*8, intent(in), optional :: beta, a
    real*8, intent(in), optional :: damping_custom(size(alpha), size(alpha))
    real*8 :: relay(3*size(alpha), 3*size(alpha))

    integer :: i_atom, j_atom, i_xyz, j_xyz
    real*8 :: lattice_dim(3)
    integer :: i_cell, j_cell, k_cell, ijk_cell(3), supercell_dim(3)
    real*8 :: r(3), r_norm, Tpp(3, 3), cell_shift(3)
    real*8 :: sigma, R_vdw_sum

    if (.not. is_periodic) then
        supercell_dim(:) = 0
    else
        lattice_dim = sqrt(sum(lattice_vector**2, 1))
        supercell_dim = ceiling(param_scs_dip_cutoff/lattice_dim)
    endif
    relay(:, :) = 0.d0
    do i_atom = 1, n_atoms
        if (my_task /= modulo(i_atom, n_tasks)) cycle
        do i_xyz = 1, 3
            relay(3*(i_atom-1)+i_xyz, 3*(i_atom-1)+i_xyz) = 1.d0/alpha(i_atom)
        enddo
    enddo
    do i_cell = -supercell_dim(1), supercell_dim(1)
    do j_cell = -supercell_dim(2), supercell_dim(2)
    do k_cell = -supercell_dim(3), supercell_dim(3)
        ijk_cell = (/ i_cell, j_cell, k_cell /)
        cell_shift = matmul(lattice_vector, ijk_cell)
        do i_atom = 1, n_atoms
            if (my_task /= modulo(i_atom, n_tasks)) cycle
            do j_atom = i_atom, n_atoms
                if (i_atom == j_atom) then
                    if (all(ijk_cell(:) == 0)) cycle
                endif
                r = coords(:, i_atom)-coords(:, j_atom)-cell_shift
                r_norm = sqrt(sum(r**2))
                if (r_norm > param_scs_dip_cutoff) cycle
                sigma = sqrt(sum( &
                    get_sigma_selfint((/ alpha(i_atom), alpha(j_atom) /))**2))
                if (present(R_vdw)) then
                    R_vdw_sum = R_vdw(i_atom)+R_vdw(j_atom)
                end if
                select case (version)
                    case ("dip,gg")
                        Tpp = T_erf_coulomb(r, sigma, 1.d0)
                    case ("1mexp,dip,gg")
                        Tpp = (1.d0-damping_1mexp(r_norm, beta*R_vdw_sum, a)) &
                            *T_erf_coulomb(r, sigma, 1.d0)
                    case ("erf,dip,gg")
                        Tpp = (1.d0-damping_erf(r_norm, beta*R_vdw_sum, a)) & 
                            *T_erf_coulomb(r, sigma, 1.d0)
                    case ("fermi,dip,gg")
                        Tpp = (1.d0-damping_fermi(r_norm, beta*R_vdw_sum, a)) &
                            * T_erf_coulomb(r, sigma, 1.d0)
                    case ("custom,dip,gg")
                        Tpp = (1.d0-damping_custom(i_atom, j_atom)) &
                            * T_erf_coulomb(r, sigma, 1.d0)
                endselect
                do i_xyz = 1, 3
                do j_xyz = 1, 3
                    relay(3*(i_atom-1)+i_xyz, 3*(j_atom-1)+j_xyz) &
                        = relay(3*(i_atom-1)+i_xyz, 3*(j_atom-1)+j_xyz) &
                            -Tpp(i_xyz, j_xyz)
                    if (i_atom /= j_atom) then
                        relay(3*(j_atom-1)+j_xyz, 3*(i_atom-1)+i_xyz) &
                            = relay(3*(j_atom-1)+j_xyz, 3*(i_atom-1)+i_xyz) &
                                -Tpp(i_xyz, j_xyz)
                    endif
                enddo
                enddo
            enddo
        enddo
    enddo
    enddo
    enddo
    call sync_sum_array (relay, size(relay))
    relay = invert_matrix(relay)
end function run_scs

function contract_polarizability(relay) result(alpha)
    implicit none

    real*8, intent(in) :: relay(:, :)
    real*8 :: alpha(size(relay, 1)/3)

    integer :: i_atom, i_xyz, j_xyz
    real*8 :: alpha_3_3(3, 3), alpha_diag(3)

    do i_atom = 1, n_atoms
        do i_xyz = 1, 3
            do j_xyz = 1, 3
                alpha_3_3(i_xyz, j_xyz) = sum(relay(i_xyz::3, 3*(i_atom-1)+j_xyz))
            enddo
        enddo
        alpha_diag = diagonalize_matrix(alpha_3_3)
        alpha(i_atom) = sum(alpha_diag)/3
    enddo
end function

function get_mbd_energy( &
        omega, alpha_0, version, R_vdw, beta, a, &
        damping_custom, potential_custom, overlap, C6, &
        omegap, modes) &
        result(ene)
    implicit none

    real*8, intent(in) :: omega(:), alpha_0(size(omega))
    character(len=*), intent(in) :: version
    real*8, intent(in), optional :: R_vdw(size(omega))
    real*8, intent(in), optional :: beta, a
    real*8, intent(in), optional :: &
        damping_custom(size(omega), size(omega)), &
        overlap(size(omega), size(omega)), C6(size(omega))
    real*8, intent(in), optional :: &
        potential_custom(size(omega), size(omega), 3, 3)
    real*8, intent(out), optional :: &
        omegap(size(omega)), modes(3*size(omega), 3*size(omega))
    real*8 :: ene

    character(len=1000) :: info_str
    real*8, allocatable :: C(:, :), C_eigs(:), coords_hyper(:, :)
    real*8, allocatable :: xi_delta(:), rhs(:)
    real*8 :: delta_ene
    logical :: is_electric_field
    integer :: i_atom, j_atom, i_xyz, j_xyz
    real*8 :: lattice_dim(3), cell_shift(3)
    integer :: i_cell, j_cell, k_cell, i_cell_total, ijk_cell(3), supercell_dim(3)
    integer :: i_atom_hyp, j_atom_hyp, hypercell_dim(3), n_cells_hyper, &
               n_atoms_hyper
    real*8 :: lattice_vector_hyper(3, 3), lattice_dim_hyper(3)
    real*8 :: prefactor, Tpp(3, 3), R_vdw_sum, r(3), r_norm, C6_ij, overlap_ij
    integer :: n_negative_eigs

    if (.not. is_periodic) then
        hypercell_dim(:) = 1
    else
        lattice_dim = sqrt(sum(lattice_vector**2, 1))
        hypercell_dim = ceiling(param_mbd_supercell_cutoff/lattice_dim)
    endif
    is_electric_field = any(electric_field /= 0.d0)
    n_cells_hyper = product(hypercell_dim(:))
    n_atoms_hyper = n_cells_hyper*n_atoms
    allocate (C(3*n_atoms_hyper, 3*n_atoms_hyper))
    C(:, :) = 0.d0
    allocate (C_eigs(3*n_atoms_hyper))
    allocate (coords_hyper(3, n_atoms_hyper))
    if (is_electric_field) then
        allocate (xi_delta(3*n_atoms_hyper))
        allocate (rhs(3*n_atoms_hyper))
        rhs(:) = 0.d0
    endif
    i_cell_total = 0
    do i_cell = 0, hypercell_dim(1)-1
    do j_cell = 0, hypercell_dim(2)-1
    do k_cell = 0, hypercell_dim(3)-1
        ijk_cell = (/ i_cell, j_cell, k_cell /)
        cell_shift = matmul(lattice_vector, ijk_cell)
        do i_atom = 1, n_atoms
            coords_hyper(:, i_cell_total*n_atoms+i_atom) &
                = coords(:, i_atom)+cell_shift
        enddo
        i_cell_total = i_cell_total+1
    enddo
    enddo
    enddo
    do i_xyz = 1, 3
        lattice_vector_hyper(i_xyz, :) = hypercell_dim*lattice_vector(i_xyz, :)
    enddo
    if (.not. is_periodic) then
        supercell_dim(:) = 0
    else
        lattice_dim_hyper = sqrt(sum(lattice_vector_hyper**2, 1))
        supercell_dim = ceiling(param_mbd_dip_cutoff/lattice_dim_hyper)
    endif
    do i_atom_hyp = 1, n_atoms_hyper
        i_atom = ind_orig(i_atom_hyp)
        if (my_task /= modulo(i_atom, n_tasks)) cycle
        do i_xyz = 1, 3
            C(3*(i_atom_hyp-1)+i_xyz, 3*(i_atom_hyp-1)+i_xyz) = omega(i_atom)**2
            if (is_electric_field) then
                rhs(3*(i_atom_hyp-1)+i_xyz) &
                    = -sqrt(alpha_0(i_atom))*omega(i_atom)*electric_field(i_xyz)
            endif
        enddo
    enddo
    do i_cell = -supercell_dim(1), supercell_dim(1)
    do j_cell = -supercell_dim(2), supercell_dim(2)
    do k_cell = -supercell_dim(3), supercell_dim(3)
        ijk_cell = (/ i_cell, j_cell, k_cell /)
        cell_shift = matmul(lattice_vector_hyper, ijk_cell)
        do i_atom_hyp = 1, n_atoms_hyper
            i_atom = ind_orig(i_atom_hyp)
            if (my_task /= modulo(i_atom, n_tasks)) cycle
            do j_atom_hyp = i_atom_hyp, n_atoms_hyper
                j_atom = ind_orig(j_atom_hyp)
                if (i_atom_hyp == j_atom_hyp) then
                    if (all(ijk_cell(:) == 0)) cycle
                endif
                r = coords_hyper(:, i_atom_hyp) &
                    -coords_hyper(:, j_atom_hyp)-cell_shift
                r_norm = sqrt(sum(r**2))
                if (r_norm > param_mbd_dip_cutoff) cycle
                if (present(R_vdw)) then
                    R_vdw_sum = R_vdw(i_atom)+R_vdw(j_atom)
                end if
                if (present(overlap)) then
                    if (i_atom < j_atom) then
                        overlap_ij = overlap(i_atom, j_atom)
                    else
                        overlap_ij = overlap(j_atom, i_atom)
                    end if
                end if
                if (present(C6)) then
                    C6_ij = combine_C6( &
                        C6(i_atom), C6(j_atom), &
                        alpha_0(i_atom), alpha_0(j_atom))
                end if
                select case (version)
                    case ("dip,1mexp")
                        Tpp = T_1mexp_coulomb(r, beta*R_vdw_sum, a)
                    case ("dip,erf")
                        Tpp = T_erf_coulomb(r, beta*R_vdw_sum, a)
                    case ("dip,fermi")
                        Tpp = T_fermi_coulomb(r, beta*R_vdw_sum, a)
                    case ("dip,overlap")
                        Tpp = T_overlap_coulomb(r, overlap_ij, C6_ij, beta, a)
                    case ("1mexp,dip")
                        Tpp = damping_1mexp(r_norm, beta*R_vdw_sum, a)*T_bare(r)
                    case ("erf,dip")
                        Tpp = damping_erf(r_norm, beta*R_vdw_sum, a)*T_bare(r)
                    case ("fermi,dip", "fermi@TS,dip", "fermi@rsSCS,dip")
                        Tpp = damping_fermi(r_norm, beta*R_vdw_sum, a)*T_bare(r)
                    case ("overlap,dip")
                        Tpp = damping_overlap(r_norm, overlap_ij, C6_ij, beta, a) & 
                            *T_bare(r)
                    case ("custom,dip")
                        Tpp = damping_custom(i_atom, j_atom)*T_bare(r)
                    case ("dip,custom")
                        Tpp = potential_custom(i_atom, j_atom, :, :)
                endselect
                prefactor = omega(i_atom)*omega(j_atom) &
                    *sqrt(alpha_0(i_atom)*alpha_0(j_atom))
                do i_xyz = 1, 3
                do j_xyz = 1, 3
                    C(3*(i_atom_hyp-1)+i_xyz, 3*(j_atom_hyp-1)+j_xyz) &
                        = C(3*(i_atom_hyp-1)+i_xyz, 3*(j_atom_hyp-1)+j_xyz) &
                          -prefactor*Tpp(i_xyz, j_xyz)
                    if (i_atom_hyp /= j_atom_hyp) then
                        C(3*(j_atom_hyp-1)+j_xyz, 3*(i_atom_hyp-1)+i_xyz) &
                            = C(3*(j_atom_hyp-1)+j_xyz, 3*(i_atom_hyp-1)+i_xyz) &
                              -prefactor*Tpp(i_xyz, j_xyz)
                    endif
                enddo
                enddo
            enddo
        enddo
    enddo
    enddo
    enddo
    call sync_sum_array (C, size(C))

    if (present(omegap)) then
        C_eigs = diagonalize_matrix(C, modes)
        omegap = sqrt(C_eigs)
    else
        C_eigs = diagonalize_matrix(C)
    end if
    n_negative_eigs = count(C_eigs(:) < 0)
    if (n_negative_eigs > 0) then
        write (info_str, "(a,1x,i10,1x,a)") &
            "CFDM Hamiltonian has", n_negative_eigs, "negative eigenvalues"
        call print_warning (info_str)
    endif
    where (C_eigs < 0) C_eigs = 0.d0
    ene = 1.d0/2*sum(sqrt(C_eigs))/n_cells_hyper-3.d0/2*sum(omega)
    if (is_electric_field) then
        call sync_sum_array (rhs, size(rhs))
        xi_delta(:) = solve_lin_sys(C, rhs)
        delta_ene = &
            (1.d0/2*expect_value(C, xi_delta) &
            -dot_product(rhs, xi_delta))/n_cells_hyper &
            +1.d0/2*sum(alpha_0)*sum(electric_field**2)
        ene = ene+delta_ene
    endif

    deallocate (coords_hyper)
    deallocate (C)
    deallocate (C_eigs)
    if (is_electric_field) then
        deallocate (xi_delta)
        deallocate (rhs)
    endif
end function get_mbd_energy

function nbody_mbd( &
        omega, alpha_0, version, R_vdw, beta, a, &
        damping_custom, potential_custom, overlap, C6) &
        result(ene)
    implicit none

    real*8, intent(in) :: omega(:), alpha_0(size(omega))
    character(len=*), intent(in) :: version
    real*8, intent(in), optional :: R_vdw(size(omega))
    real*8, intent(in), optional :: beta, a
    real*8, intent(in), optional :: &
        damping_custom(size(omega), size(omega)), &
        overlap(size(omega), size(omega)), C6(size(omega))
    real*8, intent(in), optional :: &
        potential_custom(size(omega), size(omega), 3, 3)
    real*8 :: ene(3)

    character(len=1000) :: info_str
    real*8, allocatable :: C(:, :), C_eigs(:)
    real*8 :: delta_ene
    logical :: is_electric_field
    integer :: i_atom, j_atom, i_xyz, j_xyz
    real*8 :: lattice_dim(3), cell_shift(3)
    integer :: i_cell, j_cell, k_cell, i_cell_total, ijk_cell(3), supercell_dim(3)
    integer :: i_atom_hyp, j_atom_hyp, hypercell_dim(3), n_cells_hyper, &
               n_atoms_hyper
    real*8 :: lattice_vector_hyper(3, 3), lattice_dim_hyper(3)
    real*8 :: prefactor, Tpp(3, 3), R_vdw_sum, r(3), r_norm, C6_ij, overlap_ij
    integer :: n_negative_eigs
    integer :: multi_index(param_mbd_nbody_max), i_body, j_body, i_tuple, &
               i_atom_ind, j_atom_ind, i_index

    ene(:) = 0.d0
    do i_body = 2, param_mbd_nbody_max
        allocate (C(i_body*3, i_body*3))
        allocate (C_eigs(i_body*3))
        i_tuple = 0
        multi_index(1:i_body-1) = 1
        multi_index(i_body:param_mbd_nbody_max) = 0
        do
            multi_index(i_body) = multi_index(i_body)+1
            do i_index = i_body, 2, -1
                if (multi_index(i_index) > n_atoms) then
                    multi_index(i_index) = 1
                    multi_index(i_index-1) = multi_index(i_index-1)+1
                end if
            end do
            if (multi_index(1) > n_atoms) exit
            if (any(multi_index(1:i_body-1)-multi_index(2:i_body) >= 0)) cycle
            i_tuple = i_tuple+1
            if (my_task /= modulo(i_tuple, n_tasks)) cycle
            C(:, :) = 0.d0
            do i_atom_ind = 1, i_body
                i_atom = multi_index(i_atom_ind)
                do i_xyz = 1, 3
                    C((i_atom_ind-1)*3+i_xyz, (i_atom_ind-1)*3+i_xyz) &
                        = omega(i_atom)**2
                enddo
            end do
            do i_atom_ind = 1, i_body
            do j_atom_ind = i_atom_ind+1, i_body
                i_atom = multi_index(i_atom_ind)
                j_atom = multi_index(j_atom_ind)
                r = coords(:, i_atom)-coords(:, j_atom)
                r_norm = sqrt(sum(r**2))
                if (present(R_vdw)) then
                    R_vdw_sum = R_vdw(i_atom)+R_vdw(j_atom)
                end if
                if (present(overlap)) then
                    if (i_atom < j_atom) then
                        overlap_ij = overlap(i_atom, j_atom)
                    else
                        overlap_ij = overlap(j_atom, i_atom)
                    end if
                end if
                if (present(C6)) then
                    C6_ij = combine_C6( &
                        C6(i_atom), C6(j_atom), alpha_0(i_atom), alpha_0(j_atom))
                end if
                select case (version)
                    case ("dip,1mexp")
                        Tpp = T_1mexp_coulomb(r, beta*R_vdw_sum, a)
                    case ("dip,erf")
                        Tpp = T_erf_coulomb(r, beta*R_vdw_sum, a)
                    case ("dip,fermi")
                        Tpp = T_fermi_coulomb(r, beta*R_vdw_sum, a)
                    case ("dip,overlap")
                        Tpp = T_overlap_coulomb(r, overlap_ij, C6_ij, beta, a)
                    case ("1mexp,dip")
                        Tpp = damping_1mexp(r_norm, beta*R_vdw_sum, a)*T_bare(r)
                    case ("erf,dip")
                        Tpp = damping_erf(r_norm, beta*R_vdw_sum, a)*T_bare(r)
                    case ("fermi,dip", "fermi@TS,dip", "fermi@rsSCS,dip")
                        Tpp = damping_fermi(r_norm, beta*R_vdw_sum, a)*T_bare(r)
                    case ("overlap,dip")
                        Tpp = damping_overlap(r_norm, overlap_ij, C6_ij, beta, a) &
                                *T_bare(r)
                    case ("custom,dip")
                        Tpp = damping_custom(i_atom, j_atom)*T_bare(r)
                    case ("dip,custom")
                        Tpp = potential_custom(i_atom, j_atom, :, :)
                endselect
                prefactor = omega(i_atom)*omega(j_atom) &
                            *sqrt(alpha_0(i_atom)*alpha_0(j_atom))
                do i_xyz = 1, 3
                do j_xyz = 1, 3
                    C(3*(i_atom_ind-1)+i_xyz, 3*(j_atom_ind-1)+j_xyz) &
                        = C(3*(i_atom_ind-1)+i_xyz, 3*(j_atom_ind-1)+j_xyz) &
                        -prefactor*Tpp(i_xyz, j_xyz)
                    C(3*(j_atom_ind-1)+j_xyz, 3*(i_atom_ind-1)+i_xyz) &
                        = C(3*(j_atom_ind-1)+j_xyz, 3*(i_atom_ind-1)+i_xyz) &
                        -prefactor*Tpp(i_xyz, j_xyz)
                enddo
                enddo
            end do
            end do
            C_eigs = diagonalize_matrix(C)
            where (C_eigs < 0) C_eigs = 0.d0
            ene(i_body) = ene(i_body)+1.d0/2*sum(sqrt(C_eigs))
        enddo ! end cycling over tuples
        deallocate (C)
        deallocate (C_eigs)
    enddo ! end cycling over i_body
    call sync_sum_array(ene, size(ene))
    ene(1) = 3.d0/2*sum(omega)
    do i_body = 2, min(param_mbd_nbody_max, n_atoms)
        do j_body = 1, i_body-1
            ene(i_body) = ene(i_body) &
                -nbody_coeffs(j_body, i_body, n_atoms)*ene(j_body)
        end do
    end do
    ene(1) = sum(ene(2:param_mbd_nbody_max))
end function nbody_mbd

function nbody_coeffs(k, m, N) result(a)
    integer, intent(in) :: k, m, N
    integer :: a

    integer :: i

    a = 1
    do i = N-m+1, N-k
        a = a*i
    end do
    do i = 1, m-k
        a = a/i
    end do
end function nbody_coeffs

function pairwise_mbd( &
        omega, alpha_0, version, R_vdw, beta, a, &
        damping_custom, potential_custom, overlap, C6, &
        omegap, modes) &
        result(ene)
    implicit none

    real*8, intent(in) :: omega(:), alpha_0(size(omega))
    character(len=*), intent(in) :: version
    real*8, intent(in), optional :: R_vdw(size(omega))
    real*8, intent(in), optional :: beta, a
    real*8, intent(in), optional :: &
        damping_custom(size(omega), size(omega)), &
        overlap(size(omega), size(omega)), C6(size(omega))
    real*8, intent(in), optional :: &
        potential_custom(size(omega), size(omega), 3, 3)
    real*8, intent(out), optional :: modes(6, 6, 3*size(omega), 3*size(omega))
    real*8, intent(out), optional :: omegap(6, 3*size(omega), 3*size(omega))
    real*8 :: ene

    character(len=1000) :: info_str
    real*8 :: C(6, 6), C_eigs(6), C_modes(6, 6)
    real*8 :: delta_ene
    logical :: is_electric_field
    integer :: i_atom, j_atom, i_xyz, j_xyz
    real*8 :: lattice_dim(3), cell_shift(3)
    integer :: i_cell, j_cell, k_cell, i_cell_total, ijk_cell(3), supercell_dim(3)
    integer :: i_atom_hyp, j_atom_hyp, hypercell_dim(3), n_cells_hyper, &
               n_atoms_hyper
    real*8 :: lattice_vector_hyper(3, 3), lattice_dim_hyper(3)
    real*8 :: prefactor, Tpp(3, 3), R_vdw_sum, r(3), r_norm, C6_ij, overlap_ij
    integer :: n_negative_eigs


    if (present(omegap)) then
        modes(:, :, :, :) = 0.d0
        omegap(:, :, :) = 0.d0
    end if
    ene = 0.d0
    do i_atom = 1, n_atoms
        if (my_task /= modulo(i_atom, n_tasks)) cycle
        do j_atom = i_atom+1, n_atoms
            C(:, :) = 0.d0
            do i_xyz = 1, 3
                C(i_xyz, i_xyz) = omega(i_atom)**2
                C(3+i_xyz, 3+i_xyz) = omega(j_atom)**2
            enddo
            r = coords(:, i_atom)-coords(:, j_atom)
            r_norm = sqrt(sum(r**2))
            if (present(R_vdw)) then
                R_vdw_sum = R_vdw(i_atom)+R_vdw(j_atom)
            end if
            if (present(overlap)) then
                if (i_atom < j_atom) then
                    overlap_ij = overlap(i_atom, j_atom)
                else
                    overlap_ij = overlap(j_atom, i_atom)
                end if
            end if
            if (present(C6)) then
                C6_ij = combine_C6( &
                    C6(i_atom), C6(j_atom), alpha_0(i_atom), alpha_0(j_atom))
            end if
            select case (version)
                case ("dip,1mexp")
                    Tpp = T_1mexp_coulomb(r, beta*R_vdw_sum, a)
                case ("dip,erf")
                    Tpp = T_erf_coulomb(r, beta*R_vdw_sum, a)
                case ("dip,fermi")
                    Tpp = T_fermi_coulomb(r, beta*R_vdw_sum, a)
                case ("dip,overlap")
                    Tpp = T_overlap_coulomb(r, overlap_ij, C6_ij, beta, a)
                case ("1mexp,dip")
                    Tpp = damping_1mexp(r_norm, beta*R_vdw_sum, a)*T_bare(r)
                case ("erf,dip")
                    Tpp = damping_erf(r_norm, beta*R_vdw_sum, a)*T_bare(r)
                case ("fermi,dip", "fermi@TS,dip", "fermi@rsSCS,dip")
                    Tpp = damping_fermi(r_norm, beta*R_vdw_sum, a)*T_bare(r)
                case ("overlap,dip")
                    Tpp = damping_overlap(r_norm, overlap_ij, C6_ij, beta, a) &
                            *T_bare(r)
                case ("custom,dip")
                    Tpp = damping_custom(i_atom, j_atom)*T_bare(r)
                case ("dip,custom")
                    Tpp = potential_custom(i_atom, j_atom, :, :)
            endselect
            prefactor = omega(i_atom)*omega(j_atom) &
                        *sqrt(alpha_0(i_atom)*alpha_0(j_atom))
            do i_xyz = 1, 3
                do j_xyz = 1, 3
                    C(i_xyz, 3+j_xyz) = C(i_xyz, 3+j_xyz) &
                        -prefactor*Tpp(i_xyz, j_xyz)
                    C(3+j_xyz, i_xyz) = C(3+j_xyz, i_xyz) &
                        -prefactor*Tpp(i_xyz, j_xyz)
                enddo
            enddo
            C_eigs = diagonalize_matrix(C, C_modes)
            if (present(omegap)) then
                omegap(:, i_atom, j_atom) = sqrt(C_eigs)
                modes(:, :, i_atom, j_atom) = C_modes
            end if
            where (C_eigs < 0) C_eigs = 0.d0
            ene = ene &
                +1.d0/2*sum(sqrt(C_eigs))-3.d0/2*(omega(i_atom)+omega(j_atom))
        enddo
    enddo
    call sync_sum_number (ene)
    if (present(omegap)) then
        call sync_sum_array (omegap, size(omegap))
        call sync_sum_array (modes, size(modes))
    end if
end function pairwise_mbd

function get_qho_rpa_energy( &
        alpha, version, R_vdw, beta, a, & 
        damping_custom, potential_custom) &
        result(ene)
    implicit none

    real*8, intent(in) :: alpha(:, :) ! indexing from 1 here because of f2py
    character(len=*), intent(in) :: version
    real*8, intent(in), optional :: R_vdw(size(alpha, 2))
    real*8, intent(in), optional :: beta, a
    real*8, intent(in), optional :: &
        damping_custom(size(alpha, 2), size(alpha, 2))
    real*8, intent(in), optional :: &
        potential_custom(size(alpha, 2), size(alpha, 2), 3, 3)
    real*8 :: ene(10)

    character(len=1000) :: info_str
    real*8, allocatable :: AT(:, :), onemAT(:, :), AT_eigs(:), coords_hyper(:, :)
    real*8, allocatable :: xi_delta(:), rhs(:)
    real*8 :: delta_ene
    logical :: is_electric_field
    integer :: i_atom, j_atom, i_xyz, j_xyz
    real*8 :: lattice_dim(3), cell_shift(3)
    integer :: i_cell, j_cell, k_cell, i_cell_total, ijk_cell(3), supercell_dim(3)
    integer :: i_atom_hyp, j_atom_hyp, hypercell_dim(3), n_cells_hyper
    integer :: n_atoms_hyper
    integer :: i_grid_omega, i_order
    real*8 :: lattice_vector_hyper(3, 3), lattice_dim_hyper(3)
    real*8 :: prefactor, Tpp(3, 3), sigma, R_vdw_sum, r(3), r_norm
    integer :: n_negative_eigs
    real*8 :: AT_eigs_grid(param_rpa_order_max, n_grid_omega+1)

    if (.not. is_periodic) then
        hypercell_dim(:) = 1
    else
        lattice_dim = sqrt(sum(lattice_vector**2, 1))
        hypercell_dim = ceiling(param_mbd_supercell_cutoff/lattice_dim)
    endif
    is_electric_field = any(electric_field /= 0.d0)
    n_cells_hyper = product(hypercell_dim(:))
    n_atoms_hyper = n_cells_hyper*n_atoms
    allocate (AT(3*n_atoms_hyper, 3*n_atoms_hyper))
    allocate (onemAT(3*n_atoms_hyper, 3*n_atoms_hyper))
    AT(:, :) = 0.d0
    allocate (AT_eigs(3*n_atoms_hyper))
    allocate (coords_hyper(3, n_atoms_hyper))
    if (is_electric_field) then
        allocate (xi_delta(3*n_atoms_hyper))
        allocate (rhs(3*n_atoms_hyper))
        rhs(:) = 0.d0
    endif
    i_cell_total = 0
    do i_cell = 0, hypercell_dim(1)-1
    do j_cell = 0, hypercell_dim(2)-1
    do k_cell = 0, hypercell_dim(3)-1
        ijk_cell = (/ i_cell, j_cell, k_cell /)
        cell_shift = matmul(lattice_vector, ijk_cell)
        do i_atom = 1, n_atoms
            coords_hyper(:, i_cell_total*n_atoms+i_atom) &
                = coords(:, i_atom)+cell_shift
        enddo
        i_cell_total = i_cell_total+1
    enddo
    enddo
    enddo
    do i_xyz = 1, 3
        lattice_vector_hyper(i_xyz, :) = hypercell_dim*lattice_vector(i_xyz, :)
    enddo
    if (.not. is_periodic) then
        supercell_dim(:) = 0
    else
        lattice_dim_hyper = sqrt(sum(lattice_vector_hyper**2, 1))
        supercell_dim = ceiling(param_mbd_dip_cutoff/lattice_dim_hyper)
    endif
    AT_eigs_grid(:, :) = 0.d0
    do i_grid_omega = 1, n_grid_omega+1
        if (my_task /= modulo(i_grid_omega, n_tasks)) cycle
        AT(:, :) = 0.d0
        do i_cell = -supercell_dim(1), supercell_dim(1)
        do j_cell = -supercell_dim(2), supercell_dim(2)
        do k_cell = -supercell_dim(3), supercell_dim(3)
            ijk_cell = (/ i_cell, j_cell, k_cell /)
            cell_shift = matmul(lattice_vector_hyper, ijk_cell)
            do i_atom_hyp = 1, n_atoms_hyper
                i_atom = ind_orig(i_atom_hyp)
                do j_atom_hyp = i_atom_hyp, n_atoms_hyper
                    j_atom = ind_orig(j_atom_hyp)
                    if (i_atom_hyp == j_atom_hyp) then
                        if (all(ijk_cell(:) == 0)) cycle
                    endif
                    r = coords_hyper(:, i_atom_hyp) &
                        -coords_hyper(:, j_atom_hyp)-cell_shift
                    r_norm = sqrt(sum(r**2))
                    if (r_norm > param_mbd_dip_cutoff) cycle
                    sigma = sqrt(sum(get_sigma_selfint((/ &
                        alpha(i_grid_omega, i_atom), & 
                        alpha(i_grid_omega, j_atom) &
                        /))**2))
                    if (present(R_vdw)) then
                        R_vdw_sum = R_vdw(i_atom)+R_vdw(j_atom)
                    end if
                    select case (version)
                        case ("dip,1mexp")
                            Tpp = T_1mexp_coulomb(r, beta*R_vdw_sum, a)
                        case ("dip,erf")
                            Tpp = T_erf_coulomb(r, beta*R_vdw_sum, a)
                        case ("dip,fermi")
                            Tpp = T_fermi_coulomb(r, beta*R_vdw_sum, a)
                        case ("1mexp,dip")
                            Tpp = damping_1mexp(r_norm, beta*R_vdw_sum, a) & 
                                *T_bare(r)
                        case ("erf,dip")
                            Tpp = damping_erf(r_norm, beta*R_vdw_sum, a) &
                                *T_bare(r)
                        case ("fermi,dip", "fermi@TS,dip", "fermi@rsSCS,dip")
                            Tpp = damping_fermi(r_norm, beta*R_vdw_sum, a) &
                                *T_bare(r)
                        case ("dip")
                            Tpp = T_bare(r)
                        case ("dip,gg")
                            Tpp = T_erf_coulomb(r, sigma, 1.d0)
                        case ("1mexp,dip,gg")
                            Tpp = damping_1mexp(r_norm, beta*R_vdw_sum, a) &
                                *T_erf_coulomb(r, sigma, 1.d0)
                        case ("erf,dip,gg")
                            Tpp = damping_erf(r_norm, beta*R_vdw_sum, a) &
                                *T_erf_coulomb(r, sigma, 1.d0)
                        case ("fermi,dip,gg")
                            Tpp = damping_fermi(r_norm, beta*R_vdw_sum, a) &
                                *T_erf_coulomb(r, sigma, 1.d0)
                        case ("custom,dip")
                            Tpp = damping_custom(i_atom, j_atom)*T_bare(r)
                        case ("dip,custom")
                            Tpp = potential_custom(i_atom, j_atom, :, :)
                    endselect
                    prefactor = sqrt( &
                        alpha(i_grid_omega, i_atom)*alpha(i_grid_omega, j_atom))
                    do i_xyz = 1, 3
                    do j_xyz = 1, 3
                        AT(3*(i_atom_hyp-1)+i_xyz, 3*(j_atom_hyp-1)+j_xyz) &
                            = AT(3*(i_atom_hyp-1)+i_xyz, 3*(j_atom_hyp-1)+j_xyz) &
                              +prefactor*Tpp(i_xyz, j_xyz)
                        if (i_atom_hyp /= j_atom_hyp) then
                            AT(3*(j_atom_hyp-1)+j_xyz, 3*(i_atom_hyp-1)+i_xyz) &
                                = AT(3*(j_atom_hyp-1)+j_xyz, 3*(i_atom_hyp-1)+i_xyz) &
                                  +prefactor*Tpp(i_xyz, j_xyz)
                        endif
                    enddo
                    enddo
                enddo
            enddo
        enddo
        enddo
        enddo
        onemAT = -AT
        do i_atom_hyp = 1, n_atoms_hyper
            do i_xyz = 1, 3
                onemAT(3*(i_atom_hyp-1)+i_xyz, 3*(i_atom_hyp-1)+i_xyz) = 1.d0
            enddo
        enddo
        AT_eigs = diagonalize_matrix(onemAT)
        n_negative_eigs = count(AT_eigs(:) < 0)
        if (n_negative_eigs > 0) then
            write (info_str, "(a,1x,i10,1x,a)") &
                "RPA matrix has", n_negative_eigs, "negative eigenvalues"
            call print_warning (info_str)
        endif
        where (AT_eigs < 0) AT_eigs = 0.d0
        AT_eigs_grid(1, i_grid_omega) = sum(log(AT_eigs(:)))
        AT_eigs = diagonalize_matrix(AT)
        do i_order = 2, param_rpa_order_max
            AT_eigs_grid(i_order, i_grid_omega) = sum(AT_eigs(:)**i_order)
        end do
    enddo
    call sync_sum_array (AT_eigs_grid, size(AT_eigs_grid))
    ene(1) = 1.d0/(2*pi)*sum(AT_eigs_grid(1, :)*omega_grid_w(:))/n_cells_hyper
    do i_order = 2, param_rpa_order_max
        ene(i_order) = &
            -1.d0/(2*pi*i_order)*sum(AT_eigs_grid(i_order, :)*omega_grid_w(:))/n_cells_hyper
    end do
    deallocate (coords_hyper)
    deallocate (AT)
    deallocate (onemAT)
    deallocate (AT_eigs)
    if (is_electric_field) then
        deallocate (xi_delta)
        deallocate (rhs)
    endif
end function get_qho_rpa_energy

function ind_orig(i_atom_sup) result(i_atom)
    implicit none

    integer, intent(in) :: i_atom_sup
    integer :: i_atom

    i_atom = modulo(i_atom_sup-1, n_atoms)+1
end function ind_orig

function alpha_dynamic_ts(alpha_0, C6, omega) result(alpha)
    implicit none

    real*8, intent(in) :: alpha_0(:), C6(:), omega
    real*8 :: alpha(size(alpha_0))

    alpha(:) = alpha_osc(alpha_0, omega_eff(C6, alpha_0), omega)
end function alpha_dynamic_ts

elemental function alpha_osc(alpha_0, omega, u) result(alpha)
    implicit none

    real*8, intent(in) :: alpha_0, omega, u
    real*8 :: alpha

    alpha = alpha_0/(1+(u/omega)**2)
end function alpha_osc

elemental function combine_C6 (C6_i, C6_j, alpha_0_i, alpha_0_j) result(C6_ij)
    implicit none

    real*8, intent(in) :: C6_i, C6_j, alpha_0_i, alpha_0_j
    real*8 :: C6_ij

    C6_ij = 2*C6_i*C6_j/(alpha_0_j/alpha_0_i*C6_i+alpha_0_i/alpha_0_j*C6_j)
end function combine_C6

elemental function V_to_R(V) result(R)
    implicit none

    real*8, intent(in) :: V
    real*8 :: R

    R = (3.d0*V/(4.d0*pi))**(1.d0/3)
end function V_to_R

!> Evaluates a local polarizability model of Vydrov & van Voorhis
!> \cite vydrov_dispersion_2010 .
elemental function vv_polarizability(rho, rho_grad, omega, C) result(alpha)
    implicit none

    real*8, intent(in) :: rho, rho_grad, omega, C
    real*8 :: alpha

    alpha = rho/(4*pi/3*rho+C*(rho_grad/rho)**4+omega**2)
end function vv_polarizability

function omega_eff(C6, alpha) result(omega)
    implicit none

    real*8, intent(in) :: C6(:), alpha(size(C6))
    real*8 :: omega(size(C6))

    omega = 4.d0/3*C6/alpha**2
end function omega_eff

elemental function get_sigma_selfint(alpha) result(sigma)
    implicit none

    real*8, intent(in) :: alpha
    real*8 :: sigma

    sigma = (sqrt(2.d0/pi)*alpha/3.d0)**(1.d0/3)
end function get_sigma_selfint

function get_C6_from_alpha(alpha) result(C6)
    implicit none

    real*8, intent(in) :: alpha(:, :)
    real*8 :: C6(size(alpha, 2))
    integer :: i_atom

    do i_atom = 1, n_atoms
        C6(i_atom) = 3.d0/pi*sum((alpha(:, i_atom)**2)*omega_grid_w(:))
    enddo
end function get_C6_from_alpha

function get_total_C6_from_alpha(alpha) result(C6)
    implicit none

    real*8, intent(in) :: alpha(:, :)
    real*8 :: C6

    C6 = 3.d0/pi*sum((sum(alpha, 2)**2)*omega_grid_w(:))
end function get_total_C6_from_alpha

function T_bare(rxyz) result(T)
    implicit none

    real*8, intent(in) :: rxyz(3)
    real*8 :: T(3, 3)

    integer :: i, j
    real*8 :: r_sq, r_5
    real*8 :: rarb(3, 3)

    r_sq = sum(rxyz(:)**2)
    r_5 = sqrt(r_sq)**5
    do i = 1, 3
        T(i, i) = (3.d0*rxyz(i)**2-r_sq)/r_5
        do j = i+1, 3
            T(i, j) = 3.d0*rxyz(i)*rxyz(j)/r_5
            T(j, i) = T(i, j)
        end do
    end do
end function T_bare

function damping_fermi(r, sigma, a) result(f)
    implicit none

    real*8, intent(in) :: r, sigma, a
    real*8 :: f

    f = 1.d0/(1+exp(-a*(r/sigma-1)))
end function damping_fermi

function damping_erf(r, sigma, a) result(f)
    implicit none

    real*8, intent(in) :: r, sigma, a
    real*8 :: f

    f = erf((r/sigma)**a)
end function damping_erf

function damping_1mexp(r, sigma, a) result(f)
    implicit none

    real*8, intent(in) :: r, sigma, a
    real*8 :: f

    f = 1-exp(-(r/sigma)**a)
end function damping_1mexp

function damping_overlap(r, overlap, C6, beta, a) result(f)
    implicit none

    real*8, intent(in) :: r, overlap, C6, beta, a
    real*8 :: f

    f = 1.d0-terf(-overlap/(erf(r/6)**6*C6/r**6), beta, a)
end function damping_overlap

function T_overlap_coulomb(rxyz, overlap, C6, beta, a) result(T)
    implicit none

    real*8, intent(in) :: rxyz(3), overlap, C6, beta, a
    real*8 :: T(3, 3)

    real*8 :: zeta_1, zeta_2
    real*8 :: r, erff, exp36, qene, qenep, qenepp

    r = sqrt(sum(rxyz**2))
    erff = erf(r/6)
    exp36 = exp(r**2/36)
    qene = overlap*r**6/(C6*erff**6)
    qenep = 2.d0*overlap*r**5*(-(1.d0/exp36)*r/sqrt(pi)+3.d0*erff)/(C6*erff**7)
    qenepp = (1.d0/exp36**2)*overlap*r**4/(9*C6*pi*erff**8) &
        *(42*r**2+exp36*sqrt(pi)*r*(-216+r**2)*erff & 
        +270.d0*exp36**2*pi*erff**2)
    zeta_1 = 1.d0/2*(2.d0-erf(a*(beta-qene))+erf(a*(beta+qene)) &
        +2*a*r*qenep/sqrt(pi)*(-exp(-a**2*(beta-qene)**2) &
        -exp(-a**2*(beta+qene)**2)))
    zeta_2 = 1.d0/sqrt(pi)*a*exp(-a**2*(beta+qene)**2)*r**2 &
        *(2*a**2*qenep**2*(beta*(-1.d0+exp(4*a**2*beta*qene)) &
        -qene*(1.d0+exp(4*a**2*beta*qene))) &
        +qenepp*(1.d0+exp(4*a**2*beta*qene)))
    T = zeta_1*T_bare(rxyz)+zeta_2*cart_prod(rxyz, rxyz)/sqrt(sum(rxyz**2))**5
end function T_overlap_coulomb

function T_fermi_coulomb(rxyz, sigma, a) result(T)
    implicit none

    real*8, intent(in) :: rxyz(3), sigma, a
    real*8 :: T(3, 3)

    real*8 :: r_sigma, d_r_sigma_m_1, zeta_1, zeta_2

    r_sigma = sqrt(sum(rxyz**2))/sigma
    d_r_sigma_m_1 = a*(r_sigma-1)
    zeta_1 = 1.d0/(1.d0+exp(-d_r_sigma_m_1)) &
        -a/2.d0*r_sigma/(1.d0+cosh(-d_r_sigma_m_1))
    zeta_2 = 2.d0*a**2*r_sigma**2/sinh(-d_r_sigma_m_1)**3 &
        *sinh(-d_r_sigma_m_1/2.d0)**4
    T = zeta_1*T_bare(rxyz)+zeta_2*cart_prod(rxyz, rxyz)/sqrt(sum(rxyz**2))**5
end function T_fermi_coulomb

function T_erf_coulomb(rxyz, sigma, a) result(T)
    implicit none

    real*8, intent(in) :: rxyz(3), sigma, a
    real*8 :: T(3, 3)

    real*8 :: r_sigma, zeta_1, zeta_2

    r_sigma = (sqrt(sum(rxyz**2))/sigma)**a
    zeta_1 = erf(r_sigma)-2.d0/sqrt(pi)*a*r_sigma*exp(-r_sigma**2)
    zeta_2 = -2.d0/sqrt(pi)*a*r_sigma*exp(-r_sigma**2) &
        *(1.d0+a*(-1.d0+2.d0*r_sigma**2))
    T = zeta_1*T_bare(rxyz)+zeta_2*cart_prod(rxyz, rxyz)/sqrt(sum(rxyz**2))**5
end function T_erf_coulomb

function T_1mexp_coulomb(rxyz, sigma, a) result(T)
    implicit none

    real*8, intent(in) :: rxyz(3), sigma, a
    real*8 :: T(3, 3)

    real*8 :: r_sigma, zeta_1, zeta_2

    r_sigma = (sqrt(sum(rxyz**2))/sigma)**a
    zeta_1 = 1.d0-exp(-r_sigma)-a*r_sigma*exp(-r_sigma)
    zeta_2 = -r_sigma*a*exp(-r_sigma)*(1+a*(-1+r_sigma))
    T = zeta_1*T_bare(rxyz)+zeta_2*cart_prod(rxyz, rxyz)/sqrt(sum(rxyz**2))**5
end function T_1mexp_coulomb

subroutine get_damping_parameters ( &
        xc, ts_d, ts_s_r, mbd_scs_a, mbd_ts_a, mbd_ts_erf_beta, &
        mbd_ts_fermi_beta, mbd_rsscs_a, mbd_rsscs_beta)
    implicit none

    character(len=*), intent(in) :: xc
    real*8, intent(out) :: &
        ts_d, ts_s_r, mbd_scs_a, mbd_ts_a, mbd_ts_erf_beta, &
        mbd_ts_fermi_beta, mbd_rsscs_a, mbd_rsscs_beta

    ts_d = 20.d0
    ts_s_r = 1.d0
    mbd_scs_a = 2.d0
    mbd_ts_a = 6.d0
    mbd_ts_erf_beta = 1.d0
    mbd_ts_fermi_beta = 1.d0
    mbd_rsscs_a = 6.d0
    mbd_rsscs_beta = 1.d0
    select case (xc)
        case ("pbe")
            ts_s_r = 0.94d0
            mbd_scs_a = 2.56d0
            mbd_ts_erf_beta = 1.07d0
            mbd_ts_fermi_beta = 0.81d0
            mbd_rsscs_beta = 0.83d0
        case ("pbe0")
            ts_s_r = 0.96d0
            mbd_scs_a = 2.53d0
            mbd_ts_erf_beta = 1.08d0
            mbd_ts_fermi_beta = 0.83d0
            mbd_rsscs_beta = 0.85d0
        case ("hse")
            ts_s_r = 0.96d0
            mbd_scs_a = 2.53d0
            mbd_ts_erf_beta = 1.08d0
            mbd_ts_fermi_beta = 0.83d0
            mbd_rsscs_beta = 0.85d0
        case ("blyp")
            ts_s_r = 0.62d0
        case ("b3lyp")
            ts_s_r = 0.84d0
        case ("revpbe")
            ts_s_r = 0.60d0
        case ("am05")
            ts_s_r = 0.84d0
    endselect
end subroutine get_damping_parameters

function solve_lin_sys(A_arg, B_arg) result(X)
    implicit none

    real*8, intent(in) :: A_arg(:, :), B_arg(size(A_arg, 1))
    real*8 :: X(size(B_arg))

    real*8 :: A(size(B_arg), size(B_arg))
    real*8 :: B(size(B_arg))
    integer :: i_pivot(size(B_arg))
    integer :: n
    integer :: error_flag

    A(:, :) = A_arg(:, :)
    B(:) = B_arg(:)
    n = size(B_arg)
    call DGESV (n, 1, A, n, i_pivot, B, n, error_flag)
    X(:) = B(:)
end function solve_lin_sys

function invert_matrix(A_arg, n_work_arr_arg) result(A_inv)
    implicit none

    real*8, intent(in) :: A_arg(:, :)
    integer, intent(in), optional :: n_work_arr_arg
    real*8 :: A_inv(size(A_arg, 1), size(A_arg, 1))

    character(len=1000) :: info_str
    real*8 :: A(size(A_arg, 1), size(A_arg, 2))
    integer :: i_pivot(size(A_arg, 1))
    real*8, allocatable :: work_arr(:)
    integer :: n
    integer :: n_work_arr
    real*8 :: n_work_arr_optim
    integer :: error_flag

    A(:, :) = A_arg(:, :)
    n = size(A, 1)
    call DGETRF (n, n, A, n, i_pivot, error_flag)
    if (error_flag /= 0) then
        write (info_str, "(a,1x,i5)") &
            "Matrix inversion failed in module mbd with error code", error_flag
        call print_error (info_str)
    endif
    if (present(n_work_arr_arg)) then
        n_work_arr = n_work_arr_arg
    else
        call DGETRI (n, A, n, i_pivot, n_work_arr_optim, -1, error_flag)
        n_work_arr = nint(n_work_arr_optim)
    endif
    allocate (work_arr(n_work_arr))
    call DGETRI (n, A, n, i_pivot, work_arr, n_work_arr, error_flag)
    deallocate (work_arr)
    if (error_flag /= 0) then
        write (info_str, "(a,1x,i5)") &
            "Matrix inversion failed in module mbd with error code", error_flag
        call print_error (info_str)
    endif
    A_inv(:, :) = A(:, :)
end function invert_matrix

function diagonalize_matrix(matrix, eigvecs) result(eigs)
    implicit none

    real*8, intent(in) :: matrix(:, :)
    real*8, intent(out), optional :: eigvecs(size(matrix, 1), size(matrix, 1))
    real*8 :: eigs(size(matrix, 1))

    character(len=1000) :: info_str
    real*8 :: A(size(matrix, 1), size(matrix, 2))
    real*8, allocatable :: work_arr(:)
    integer :: n
    real*8 :: n_work_arr
    integer :: error_flag
    character(len=1) :: mode

    if (present(eigvecs)) then
        mode = 'V'
    else
        mode = 'N'
    end if
    A = matrix
    n = size(A, 1)
    call DSYEV (mode, "U", n, A, n, eigs, n_work_arr, -1, error_flag)
    allocate (work_arr(nint(n_work_arr)))
    call DSYEV (mode, "U", n, A, n, eigs, work_arr, size(work_arr), error_flag)
    deallocate (work_arr)
    if (error_flag /= 0) then
        write (info_str, "(a,1x,i5)") &
            "Matrix diagonalization failed in module mbd with error code", &
            error_flag
        call print_error (info_str)
    endif
    if (mode == 'V') then
        eigvecs = A
    end if
end function diagonalize_matrix

function eye(N)
    implicit none

    integer, intent(in) :: N
    real*8 :: eye(N, N)

    integer :: i

    eye(:, :) = 0.d0
    do i = 1, N
        eye(i, i) = 1.d0
    enddo
end function eye

function cart_prod(a, b) result(c)
    implicit none

    real*8, intent(in) :: a(:), b(:)
    real*8 :: c(size(a), size(b))

    integer :: i, j

    do i = 1, size(a)
        do j = 1, size(b)
            c(i, j) = a(i)*b(j)
        enddo
    enddo
end function cart_prod

function expect_value(O, x) result(y)
    implicit none

    real*8, intent(in) :: O(:, :), x(:)
    real*8 :: y
    integer :: i, j

    y = 0.d0
    do i = 1, size(x)
        do j = 1, size(x)
            y = y+x(i)*x(j)*O(i, j)
        enddo
    enddo
end function expect_value

elemental function terf(r, r0, a)
    implicit none

    real*8, intent(in) :: r, r0, a
    real*8 :: terf

    terf = 0.5d0*(erf(a*(r+r0))+erf(a*(r-r0)))
end function terf

end module mbd
