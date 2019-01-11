! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
#ifndef MBD_TYPE
module mbd_rpa

use mbd_constants
use mbd_damping, only: damping_t
use mbd_dipole, only: dipole_matrix
use mbd_formulas, only: sigma_selfint
use mbd_geom, only: geom_t
use mbd_matrix, only: matrix_re_t, matrix_cplx_t
use mbd_utils, only: result_t, tostr

implicit none

private
public :: get_mbd_rpa_energy

interface get_mbd_rpa_energy
    module procedure get_mbd_rpa_energy_real
    module procedure get_mbd_rpa_energy_complex
end interface

contains

#   define MBD_TYPE 0
#endif

#if MBD_TYPE == 0
type(result_t) function get_mbd_rpa_energy_real( &
        geom, alpha, damp) result(res)
#elif MBD_TYPE == 1
type(result_t) function get_mbd_rpa_energy_complex( &
        geom, alpha, damp, q) result(res)
#endif
    type(geom_t), intent(inout) :: geom
    real(dp), intent(in) :: alpha(:, 0:)
    type(damping_t), intent(in) :: damp
#if MBD_TYPE == 1
    real(dp), intent(in) :: q(3)
#endif

#if MBD_TYPE == 0
    type(matrix_re_t) :: relay, AT
#elif MBD_TYPE == 1
    type(matrix_cplx_t) :: relay, AT
#endif
    complex(dp), allocatable :: eigs(:)
    integer :: i_freq, i, my_i_atom, n_order, n_negative_eigs
    type(damping_t) :: damp_alpha

    res%energy = 0d0
    damp_alpha = damp
    allocate (eigs(3*geom%siz()))
    do i_freq = 0, ubound(geom%freq, 1)
        damp_alpha%sigma = sigma_selfint(alpha(:, i_freq))
        ! relay = T
#if MBD_TYPE == 0
        relay = dipole_matrix(geom, damp_alpha)
#elif MBD_TYPE == 1
        relay = dipole_matrix(geom, damp_alpha, q=q)
#endif
        do my_i_atom = 1, size(relay%idx%i_atom)
            associate ( &
                    i_atom => relay%idx%i_atom(my_i_atom), &
                    relay_sub => relay%val(3*(my_i_atom-1)+1:, :) &
            )
                relay_sub(:3, :) = relay_sub(:3, :)*alpha(i_atom, i_freq)
            end associate
        end do
        ! relay = alpha*T
        if (geom%get_rpa_orders) AT = relay
        ! relay = 1+alpha*T
        call relay%add_diag_scalar(1d0)
        call geom%clock(23)
        eigs = relay%eigvals(geom%exc, destroy=.true.)
        call geom%clock(-23)
        if (geom%has_exc()) return
        ! The count construct won't work here due to a bug in Cray compiler
        ! Has to manually unroll the counting
        n_negative_eigs = 0
        do i = 1, size(eigs)
           if (dble(eigs(i)) < 0) n_negative_eigs = n_negative_eigs + 1
        end do
        if (n_negative_eigs > 0) then
            geom%exc%code = MBD_EXC_NEG_EIGVALS
            geom%exc%msg = "1+AT matrix has " // &
                trim(tostr(n_negative_eigs)) // " negative eigenvalues"
            return
        end if
        res%energy = res%energy + &
            1d0/(2*pi)*sum(log(dble(eigs)))*geom%freq(i_freq)%weight
        if (geom%get_rpa_orders) then
            call geom%clock(24)
            eigs = AT%eigvals(geom%exc, destroy=.true.)
            call geom%clock(-24)
            if (geom%has_exc()) return
            allocate (res%rpa_orders(geom%param%rpa_order_max))
            do n_order = 2, geom%param%rpa_order_max
                res%rpa_orders(n_order) = res%rpa_orders(n_order) &
                    +(-1d0/(2*pi)*(-1)**n_order &
                    *sum(dble(eigs)**n_order)/n_order) &
                    *geom%freq(i_freq)%weight
            end do
        end if
    end do
end function

#if MBD_TYPE == 0
#   undef MBD_TYPE
#   define MBD_TYPE 1
#   include "mbd_rpa.F90"

end module

#endif
