module mbd

use mbd_interface, only: &
    sync_sum_array, sync_sum_number, print_error, print_warning, print_log

implicit none

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
        coords, &
        C6, &
        alpha_0, &
        version, &
        R_vdw, s_R, d, &
        overlap, &
        damping_custom, &
        lattice_vector, &
        my_task, n_tasks) &
        result(ene)
    implicit none

    real*8, intent(in) :: &
       coords(:, :), &
       C6(size(coords, 1)), &
       alpha_0(size(coords, 1))
    character(len=*), intent(in) :: version
    real*8, intent(in), optional :: &
        R_vdw(size(coords, 1)), &
        s_R, d, &
        overlap(size(coords, 1), size(coords, 1)), &
        damping_custom(size(coords, 1), size(coords, 1)), &
        lattice_vector(3, 3)
    integer, intent(in), optional :: my_task, n_tasks
    real*8 :: ene

    real*8 :: C6_ij, r(3), r_norm, f_damp, R_vdw_ij, overlap_ij
    real*8 :: ene_add, ene_contrib
    real*8 :: r_cell(3)
    integer :: n_shell, i_cell, j_cell, k_cell, ijk_cell(3)
    integer :: i_atom, j_atom, i_range, range_g_cell(3), g_cell(3)
    real*8, parameter :: step = 100.d0
    logical :: &
        is_periodic = .false., &
        is_parallel = .false.

    if (present(lattice_vector)) then
        if (any(lattice_vector > 0.d0)) then
            is_periodic = .true.
        end if
    end if
    if (present(n_tasks)) then
        if (n_tasks > 0) then
            is_parallel = .true.
        end if
    end if
    ene = 0.d0
    i_range = 0
    do
        ene_add = 0.d0
        i_range = i_range+1
        if (is_periodic) then
            range_g_cell = supercell_circum(lattice_vector, i_range*step)
        else
            range_g_cell(:) = 0
        end if
        g_cell = (/ 0, 0, -1 /)
        do i_cell = 1, product(1+2*range_g_cell)
            call shift_cell (g_cell, -range_g_cell, range_g_cell)
            if (is_parallel .and. is_periodic) then
                if (my_task /= modulo(i_cell, n_tasks)) cycle
            end if
            if (is_periodic) then
                r_cell = matmul(g_cell, lattice_vector)
            else
                r_cell(:) = 0.d0
            end if
            do i_atom = 1, size(coords, 1)
                if (is_parallel .and. .not. is_periodic) then
                    if (my_task /= modulo(i_atom, n_tasks)) cycle
                end if
                do j_atom = 1, i_atom
                    if (i_cell == 1) then
                        if (i_atom == j_atom) cycle
                    end if
                    r = coords(:, i_atom)-coords(:, j_atom)-r_cell
                    r_norm = sqrt(sum(r**2))
                    if (r_norm > i_range*step .or. r_norm < (i_range-1)*step) cycle
                    C6_ij = combine_C6(C6(i_atom), C6(j_atom), &
                        alpha_0(i_atom), alpha_0(j_atom))
                    if (present(R_vdw)) then
                        R_vdw_ij = R_vdw(i_atom)+R_vdw(j_atom)
                    end if
                    if (present(overlap)) then
                        overlap_ij = overlap(i_atom, j_atom)
                    end if
                    select case (version)
                        case ("fermi")
                            f_damp = damping_fermi(r_norm, s_R*R_vdw_ij, d)
                        case ("fermi2")
                            f_damp = damping_fermi(r_norm, s_R*R_vdw_ij, d)**2
                        case ("erf")
                            f_damp = damping_erf(r_norm, s_R*R_vdw_ij, d)
                        case ("1mexp")
                            f_damp = damping_1mexp(r_norm, s_R*R_vdw_ij, d)
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
                    ene_add = ene_add + (-C6_ij/r_norm**6*f_damp)
                enddo
            enddo
        enddo
        call sync_sum_number (ene_add)
        ene = ene + ene_add
        if (.not. is_periodic) exit
        if (i_range > 1 .and. abs(ene_add) < param_energy_accuracy) then
            exit
        endif
    enddo
end function get_ts_energy

function build_dipole_matrix( &
        coords, &
        version, &
        alpha, &
        R_vdw, beta, a, &
        overlap, &
        C6, &
        damping_custom, &
        potential_custom, &
        lattice_vector, &
        dipole_cutoff, &
        k_point, &
        my_task, n_tasks) &
        result(T)
    implicit none

    real*8, intent(in) :: coords(:, :)
    character(len=*), intent(in) :: version
    real*8, intent(in), optional :: &
        alpha(size(coords, 1)), &
        R_vdw(size(coords, 1)), &
        beta, a, &
        overlap(size(coords, 1), size(coords, 1)), &
        C6(size(coords, 1)), &
        damping_custom(size(coords, 1), size(coords, 1)), &
        potential_custom(size(coords, 1), size(coords, 1), 3, 3), &
        lattice_vector(3, 3), &
        dipole_cutoff, &
        k_point(3)
    integer, intent(in), optional :: my_task, n_tasks
    real*8 :: T(3*size(coords, 1), 3*size(coords, 1))

    real*8 :: Tpp(3, 3)
    real*8 :: r_cell(3), r(3), r_norm
    real*8 :: R_vdw_ij, C6_ij, overlap_ij, sigma_ij
    integer :: i_atom, j_atom, i_cell, g_cell(3), range_g_cell(3)
    logical :: &
        is_periodic = .false., &
        is_parallel = .false.

    if (present(lattice_vector)) then
        if (any(lattice_vector > 0.d0)) then
            is_periodic = .true.
        end if
    end if
    if (present(n_tasks)) then
        if (n_tasks > 0) then
            is_parallel = .true.
        end if
    end if
    if (is_periodic) then
        range_g_cell = supercell_circum(lattice_vector, dipole_cutoff)
    else
        range_g_cell(:) = 0
    end if
    T(:, :) = 0.d0
    g_cell = (/ 0, 0, -1 /)
    do i_cell = 1, product(1+2*range_g_cell)
        call shift_cell (g_cell, -range_g_cell, range_g_cell)
        if (is_parallel .and. is_periodic) then
            if (my_task /= modulo(i_cell, n_tasks)) cycle
        end if
        if (is_periodic) then
            r_cell = matmul(g_cell, lattice_vector)
        else
            r_cell(:) = 0.d0
        end if
        do i_atom = 1, size(coords, 1)
            if (is_parallel .and. .not. is_periodic) then
                if (my_task /= modulo(i_atom, n_tasks)) cycle
            end if
            do j_atom = 1, i_atom
                if (i_cell == 1) then
                    if (i_atom == j_atom) cycle
                end if
                r = coords(:, i_atom)-coords(:, j_atom)-r_cell
                r_norm = sqrt(sum(r**2))
                if (is_periodic) then
                    if (r_norm > dipole_cutoff) cycle
                end if
                if (present(R_vdw)) then
                    R_vdw_ij = R_vdw(i_atom)+R_vdw(j_atom)
                end if
                if (present(alpha)) then
                    sigma_ij = sqrt(sum(get_sigma_selfint(alpha((/ i_atom , j_atom /)))**2))
                end if
                if (present(overlap)) then
                    overlap_ij = overlap(i_atom, j_atom)
                end if
                if (present(C6)) then
                    C6_ij = combine_C6(C6(i_atom), C6(j_atom), alpha(i_atom), alpha(j_atom))
                end if
                select case (version)
                    case ("dip,1mexp")
                        Tpp = T_1mexp_coulomb(r, beta*R_vdw_ij, a)
                    case ("dip,erf")
                        Tpp = T_erf_coulomb(r, beta*R_vdw_ij, a)
                    case ("dip,fermi")
                        Tpp = T_fermi_coulomb(r, beta*R_vdw_ij, a)
                    case ("dip,overlap")
                        Tpp = T_overlap_coulomb(r, overlap_ij, C6_ij, beta, a)
                    case ("1mexp,dip")
                        Tpp = damping_1mexp(r_norm, beta*R_vdw_ij, a)*T_bare(r)
                    case ("erf,dip")
                        Tpp = damping_erf(r_norm, beta*R_vdw_ij, a)*T_bare(r)
                    case ("fermi,dip", "fermi@TS,dip", "fermi@rsSCS,dip")
                        Tpp = damping_fermi(r_norm, beta*R_vdw_ij, a)*T_bare(r)
                    case ("overlap,dip")
                        Tpp = damping_overlap(r_norm, overlap_ij, C6_ij, beta, a)*T_bare(r)
                    case ("custom,dip")
                        Tpp = damping_custom(i_atom, j_atom)*T_bare(r)
                    case ("dip,custom")
                        Tpp = potential_custom(i_atom, j_atom, :, :)
                    case ("dip,gg")
                        Tpp = T_erf_coulomb(r, sigma_ij, 1.d0)
                    case ("1mexp,dip,gg")
                        Tpp = (1.d0-damping_1mexp(r_norm, beta*R_vdw_ij, a)) &
                            *T_erf_coulomb(r, sigma_ij, 1.d0)
                    case ("erf,dip,gg")
                        Tpp = (1.d0-damping_erf(r_norm, beta*R_vdw_ij, a)) & 
                            *T_erf_coulomb(r, sigma_ij, 1.d0)
                    case ("fermi,dip,gg")
                        Tpp = (1.d0-damping_fermi(r_norm, beta*R_vdw_ij, a)) &
                            *T_erf_coulomb(r, sigma_ij, 1.d0)
                    case ("custom,dip,gg")
                        Tpp = (1.d0-damping_custom(i_atom, j_atom)) &
                            *T_erf_coulomb(r, sigma_ij, 1.d0)
                end select
                T(3*(i_atom-1)+1:3*(i_atom-1)+3, &
                    3*(j_atom-1)+1:3*(j_atom-1)+3) = &
                    T(3*(i_atom-1)+1:3*(i_atom-1)+3, &
                    3*(j_atom-1)+1:3*(j_atom-1)+3) + Tpp
                T(3*(j_atom-1)+1:3*(j_atom-1)+3, &
                    3*(i_atom-1)+1:3*(i_atom-1)+3) = &
                    T(3*(j_atom-1)+1:3*(j_atom-1)+3, &
                    3*(i_atom-1)+1:3*(i_atom-1)+3) + transpose(Tpp)
            end do ! j_atom
        end do ! i_atom
    end do ! i_cell
    if (is_parallel) then
        call sync_sum_array (T, size(T))
    end if
end function build_dipole_matrix

function run_scs( &
        coords, &
        alpha, &
        version, &
        R_vdw, beta, a, &
        damping_custom, &
        lattice_vector, &
        my_task, n_tasks) & 
        result(alpha_scs)
    implicit none

    real*8, intent(in) :: &
        coords(:, :), &
        alpha(size(coords, 1))
    character(len=*), intent(in) :: version
    real*8, intent(in), optional :: &
        R_vdw(size(coords, 1)), &
        beta, a, &
        damping_custom(size(coords, 1), size(coords, 1)), &
        lattice_vector(3, 3)
    integer, intent(in), optional :: my_task, n_tasks
    real*8 :: alpha_scs(size(alpha, 1), size(alpha, 2))

    real*8 :: alpha_full(3*size(coords, 1), 3*size(coords, 1))
    real*8 :: T(3*size(coords, 1), 3*size(coords, 1))
    integer :: i_atom, i_xyz, i_grid_omega
    logical :: is_parallel

    is_parallel = .false.
    if (present(n_tasks)) then
        if (n_tasks > 0) then
            is_parallel = .true.
        end if
    end if
    alpha_scs(:, :) = 0.d0
    do i_grid_omega = 0, n_grid_omega
        if (is_parallel) then
            if (my_task /= modulo(i_grid_omega, n_tasks)) cycle
        end if
    T = build_dipole_matrix( &
            coords, version, alpha(i_grid_omega+1, :), R_vdw, beta, a, &
            damping_custom=damping_custom, lattice_vector=lattice_vector, &
            dipole_cutoff=param_mbd_dip_cutoff)
        alpha_full = -T
    do i_atom = 1, size(coords, 1)
        do i_xyz = 1, 3
                alpha_full(3*(i_atom-1)+i_xyz, 3*(i_atom-1)+i_xyz) = &
                    1.d0/alpha(i_grid_omega+1, i_atom)
        end do
    end do
        alpha_full = invert_matrix(alpha_full)
        alpha_scs(i_grid_omega+1, :) = contract_polarizability(alpha_full)
    end do
    if (is_parallel) then
        call sync_sum_array (alpha_scs, size(alpha_scs))
    end if
end function run_scs

function contract_polarizability(alpha_3n_3n) result(alpha)
    implicit none

    real*8, intent(in) :: alpha_3n_3n(:, :)
    real*8 :: alpha(size(alpha_3n_3n, 1)/3)

    integer :: i_atom, i_xyz, j_xyz
    real*8 :: alpha_3_3(3, 3), alpha_diag(3)

    do i_atom = 1, size(alpha)
        do i_xyz = 1, 3
            do j_xyz = 1, 3
                alpha_3_3(i_xyz, j_xyz) = sum(alpha_3n_3n(i_xyz::3, 3*(i_atom-1)+j_xyz))
            enddo
        enddo
        alpha_diag = diagonalize_matrix(alpha_3_3)
        alpha(i_atom) = sum(alpha_diag)/3
    enddo
end function

function get_mbd_energy( &
        coords, &
        alpha_0, &
        omega, &
        version, &
        R_vdw, beta, a, &
        overlap, &
        C6, &
        damping_custom, &
        potential_custom, &
        lattice_vector, &
        k_point, &
        eigenomega, modes, &
        my_task, n_tasks) &
        result(ene)
    implicit none

    real*8, intent(in) :: &
        coords(:, :), &
        alpha_0(size(coords, 1)), &
        omega(size(coords, 1))
    character(len=*), intent(in) :: version
    real*8, intent(in), optional :: &
        R_vdw(size(coords, 1)), &
        beta, a, &
        overlap(size(coords, 1), size(coords, 1)), &
        C6(size(coords, 1)), &
        damping_custom(size(coords, 1), size(coords, 1)), &
        potential_custom(size(coords, 1), size(coords, 1), 3, 3), &
        lattice_vector(3, 3), &
        k_point(3)
    real*8, intent(out), optional :: &
        eigenomega(3*size(coords, 1)), &
        modes(3*size(coords, 1), 3*size(coords, 1))
    integer, intent(in), optional :: my_task, n_tasks
    real*8 :: ene

    character(len=1000) :: info_str
    real*8, dimension(3*size(coords, 1), 3*size(coords, 1)) :: T, V
    real*8 :: eigs(3*size(coords, 1))
    integer :: i_atom, j_atom, i_xyz, j_xyz
    integer :: i_cell, j_cell, k_cell, i_cell_total, ijk_cell(3), range_g_cell(3)
    real*8 :: prefactor, Tpp(3, 3), R_vdw_sum, r(3), r_norm, C6_ij, overlap_ij
    integer :: n_negative_eigs, g_cell(3), i_order

    V(:, :) = 0.d0
    do i_atom = 1, size(coords, 1)
        do i_xyz = 1, 3
            V(3*(i_atom-1)+i_xyz, 3*(i_atom-1)+i_xyz) = omega(i_atom)**2
        end do
    end do
    T = build_dipole_matrix( &
        coords, version, alpha_0, R_vdw, beta, a, overlap, C6, damping_custom, &
        potential_custom, lattice_vector, param_mbd_dip_cutoff, k_point, &
        my_task, n_tasks)
    do i_atom = 1, size(coords, 1)
        do j_atom = 1, size(coords, 1)
            if (i_atom == j_atom) cycle
            V(3*(i_atom-1)+1:3*(i_atom-1)+3, 3*(j_atom-1)+1:3*(j_atom-1)+3) = &
                omega(i_atom)*omega(j_atom)*sqrt(alpha_0(i_atom)*alpha_0(j_atom))* &
                T(3*(i_atom-1)+1:3*(i_atom-1)+3, 3*(j_atom-1)+1:3*(j_atom-1)+3)
        end do
    end do
    if (present(eigenomega)) then
        eigs = diagonalize_matrix(V, modes)
        eigenomega = sqrt(eigs)
    else
        eigs = diagonalize_matrix(V)
    end if
    n_negative_eigs = count(eigs(:) < 0)
    if (n_negative_eigs > 0) then
        write (info_str, "(a,1x,i10,1x,a)") &
            "CFDM Hamiltonian has", n_negative_eigs, "negative eigenvalues"
        call print_warning (info_str)
    endif
    where (eigs < 0) eigs = 0.d0
    ene = 1.d0/2*sum(sqrt(eigs))-3.d0/2*sum(omega)
end function get_mbd_energy

function mbd_nbody( &
        coords, &
        alpha_0, &
        omega, &
        version, &
        R_vdw, beta, a, &
        my_task, n_tasks) &
        result(ene_orders)
    implicit none

    real*8, intent(in) :: &
        coords(:, :), &
        alpha_0(size(coords, 1)), &
        omega(size(coords, 1)), &
        R_vdw(size(coords, 1)), &
        beta, a
    character(len=*), intent(in) :: version
    integer, intent(in), optional :: my_task, n_tasks
    real*8 :: ene_orders(20)

    integer :: &
        multi_index(param_mbd_nbody_max), i_body, j_body, i_tuple, &
        i_atom_ind, j_atom_ind, i_index
    real*8 :: ene
    logical :: &
        is_parallel = .false.

    if (present(n_tasks)) then
        if (n_tasks > 0) then
            is_parallel = .true.
        end if
    end if
    ene_orders(:) = 0.d0
    do i_body = 2, param_mbd_nbody_max
        i_tuple = 0
        multi_index(1:i_body-1) = 1
        multi_index(i_body:param_mbd_nbody_max) = 0
        do
            multi_index(i_body) = multi_index(i_body)+1
            do i_index = i_body, 2, -1
                if (multi_index(i_index) > size(coords, 1)) then
                    multi_index(i_index) = 1
                    multi_index(i_index-1) = multi_index(i_index-1)+1
                end if
            end do
            if (multi_index(1) > size(coords, 1)) exit
            if (any(multi_index(1:i_body-1)-multi_index(2:i_body) >= 0)) cycle
            i_tuple = i_tuple+1
            if (is_parallel) then
                if (my_task /= modulo(i_tuple, n_tasks)) cycle
            end if
            ene = get_mbd_energy( &
                coords(multi_index(1:i_body), :), &
                alpha_0(multi_index(1:i_body)), &
                omega(multi_index(1:i_body)), &
                version, &
                R_vdw(multi_index(1:i_body)), &
                beta, a)
            ene_orders(i_body) = ene_orders(i_body) &
                +ene+3.d0/2*sum(omega(multi_index(1:i_body)))
        end do ! i_tuple
    end do ! i_body
    if (is_parallel) then
        call sync_sum_array(ene_orders, size(ene_orders))
    end if
    ene_orders(1) = 3.d0/2*sum(omega)
    do i_body = 2, min(param_mbd_nbody_max, size(coords, 1))
        do j_body = 1, i_body-1
            ene_orders(i_body) = ene_orders(i_body) &
                -nbody_coeffs(j_body, i_body, size(coords, 1))*ene_orders(j_body)
        end do
    end do
    ene_orders(1) = sum(ene_orders(2:param_mbd_nbody_max))
end function mbd_nbody

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

function get_qho_rpa_energy( &
        coords, &
        alpha, &
        version, &
        R_vdw, beta, a, &
        overlap, &
        C6, &
        damping_custom, &
        potential_custom, &
        rpa_orders, &
        my_task, n_tasks) &
        result(ene)
    implicit none

    real*8, intent(in) :: &
        coords(:, :), &
        alpha(:, :)
    character(len=*), intent(in) :: version
    real*8, intent(in), optional :: &
        R_vdw(size(coords, 1)), &
        beta, a, &
        overlap(size(coords, 1), size(coords, 1)), &
        C6(size(coords, 1)), &
        damping_custom(size(coords, 1), size(coords, 1)), &
        potential_custom(size(coords, 1), size(coords, 1), 3, 3)
    real*8, intent(out), optional :: rpa_orders(20)
    integer, intent(in), optional :: my_task, n_tasks
    real*8 :: ene

    real*8, dimension(3*size(coords, 1), 3*size(coords, 1)) :: &
        T, big_alpha, AT, one
    real*8 :: eigs(3*size(coords, 1))
    integer :: i_atom, i_xyz, i_grid_omega
    integer :: n_order
    logical :: &
        is_parallel = .false.

    if (present(n_tasks)) then
        if (n_tasks > 0) then
            is_parallel = .true.
        end if
    end if
    one = identity_matrix(3*size(coords, 1))
    do i_grid_omega = 0, n_grid_omega
        if (is_parallel) then
            if (my_task /= modulo(i_grid_omega, n_tasks)) cycle
        end if
        T = build_dipole_matrix( &
            coords, version, alpha(i_grid_omega+1, :), R_vdw, beta, a, &
            overlap, C6, damping_custom, potential_custom)
        big_alpha(:, :) = 0.d0
        do i_atom = 1, size(coords, 1)
            do i_xyz = 1, 3
                big_alpha(3*(i_atom-1)+i_xyz, 3*(i_atom-1)+i_xyz) = &
                    alpha(i_grid_omega+1, i_atom)
            end do
        end do
        AT = matmul(big_alpha, T)
        eigs = diagonalize_matrix(one-AT)
        ene = ene + 1.d0/(2*pi)*sum(log(eigs))*omega_grid_w(i_grid_omega)
        if (present(rpa_orders)) then
            eigs = diagonalize_matrix(AT)
            do n_order = 2, param_rpa_order_max
                rpa_orders(n_order) = rpa_orders(n_order) + &
                    (-1.d0/(2*pi)*sum(eigs**n_order)/n_order) &
                    *omega_grid_w(i_grid_omega)
            end do
        end if
    end do
    if (is_parallel) then
        call sync_sum_number (ene)
        if (present(rpa_orders)) then
            call sync_sum_array (rpa_orders, size(rpa_orders))
        end if
    end if
end function get_qho_rpa_energy

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

    do i_atom = 1, size(alpha, 2)
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

function invert_matrix(A_arg) result(A_inv)
    implicit none

    real*8, intent(in) :: A_arg(:, :)
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
    call DGETRI (n, A, n, i_pivot, n_work_arr_optim, -1, error_flag)
    n_work_arr = nint(n_work_arr_optim)
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

function supercell_circum(uc, radius) result(sc)
    implicit none

    real*8, intent(in) :: uc(3, 3), radius
    real*8 :: sc(3)

    real*8 :: ruc(3, 3)

    ruc = 2*pi*invert_matrix(transpose(uc))
    sc = &
        ceiling(radius/sqrt(sum((uc*(diag(1.d0/sqrt(sum(ruc**2, 2)))*ruc))**2, 2))-.5d0)
end function supercell_circum

subroutine shift_cell (ijk, first_cell, last_cell)
    implicit none

    integer, intent(inout) :: ijk(3)
    integer, intent(in) :: first_cell(3), last_cell(3)

    integer :: i_dim, i

    do i_dim = 3, 1, -1
        i = ijk(i_dim)+1
        if (i <= last_cell(i_dim)) then
            ijk(i_dim) = i
            exit
        else
            ijk(i_dim) = first_cell(i_dim)
        end if
    end do
end subroutine shift_cell

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

function identity_matrix(n) result(A)
    implicit none

    integer, intent(in) :: n
    real*8 :: A(n, n)

    integer :: i

    A(:, :) = 0.d0
    do i = 1, n
        A(i, i) = 1.d0
    enddo
end function identity_matrix

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

function diag(a) result(b)
    implicit none

    real*8, intent(in) :: a(:)
    real*8 :: b(size(a), size(a))

    integer :: i

    b = 0.d0
    do i = 1, size(a)
        b(i, i) = a(i)
    end do
end function diag

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
