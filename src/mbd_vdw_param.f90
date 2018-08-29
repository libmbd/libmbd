! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
module mbd_vdw_param

use mbd_constants
use mbd_common, only: lower

implicit none

private
public :: ts_vdw_params, tssurf_vdw_params, species_index

real(dp), parameter :: ts_vdw_params(3, 86) = reshape([ &
    4.5d0, 6.5d0, 3.1d0, &  ! H
    1.38d0, 1.46d0, 2.65d0, &  ! He
    164.2d0, 1387.0d0, 4.16d0, &  ! Li
    38.0d0, 214.0d0, 4.17d0, &  ! Be
    21.0d0, 99.5d0, 3.89d0, &  ! B
    12.0d0, 46.6d0, 3.59d0, &  ! C
    7.4d0, 24.2d0, 3.34d0, &  ! N
    5.4d0, 15.6d0, 3.19d0, &  ! O
    3.8d0, 9.52d0, 3.04d0, &  ! F
    2.67d0, 6.38d0, 2.91d0, &  ! Ne
    162.7d0, 1556.0d0, 3.73d0, &  ! Na
    71.0d0, 627.0d0, 4.27d0, &  ! Mg
    60.0d0, 528.0d0, 4.33d0, &  ! Al
    37.0d0, 305.0d0, 4.2d0, &  ! Si
    25.0d0, 185.0d0, 4.01d0, &  ! P
    19.6d0, 134.0d0, 3.86d0, &  ! S
    15.0d0, 94.6d0, 3.71d0, &  ! Cl
    11.1d0, 64.3d0, 3.55d0, &  ! Ar
    292.9d0, 3897.0d0, 3.71d0, &  ! K
    160.0d0, 2221.0d0, 4.65d0, &  ! Ca
    120.0d0, 1383.0d0, 4.59d0, &  ! Sc
    98.0d0, 1044.0d0, 4.51d0, &  ! Ti
    84.0d0, 832.0d0, 4.44d0, &  ! V
    78.0d0, 602.0d0, 3.99d0, &  ! Cr
    63.0d0, 552.0d0, 3.97d0, &  ! Mn
    56.0d0, 482.0d0, 4.23d0, &  ! Fe
    50.0d0, 408.0d0, 4.18d0, &  ! Co
    48.0d0, 373.0d0, 3.82d0, &  ! Ni
    42.0d0, 253.0d0, 3.76d0, &  ! Cu
    40.0d0, 284.0d0, 4.02d0, &  ! Zn
    60.0d0, 498.0d0, 4.19d0, &  ! Ga
    41.0d0, 354.0d0, 4.2d0, &  ! Ge
    29.0d0, 246.0d0, 4.11d0, &  ! As
    25.0d0, 210.0d0, 4.04d0, &  ! Se
    20.0d0, 162.0d0, 3.93d0, &  ! Br
    16.8d0, 129.6d0, 3.82d0, &  ! Kr
    319.2d0, 4691.0d0, 3.72d0, &  ! Rb
    199.0d0, 3170.0d0, 4.54d0, &  ! Sr
    126.737d0, 1968.58d0, 4.8151d0, &  ! Y
    119.97d0, 1677.91d0, 4.53d0, &  ! Zr
    101.603d0, 1263.61d0, 4.2365d0, &  ! Nb
    88.4225785d0, 1028.73d0, 4.099d0, &  ! Mo
    80.083d0, 1390.87d0, 4.076d0, &  ! Tc
    65.895d0, 609.754d0, 3.9953d0, &  ! Ru
    56.1d0, 469.0d0, 3.95d0, &  ! Rh
    23.68d0, 157.5d0, 3.66d0, &  ! Pd
    50.6d0, 339.0d0, 3.82d0, &  ! Ag
    39.7d0, 452.0d0, 3.99d0, &  ! Cd
    70.22d0, 707.046d0, 4.23198d0, &  ! In
    55.95d0, 587.417d0, 4.303d0, &  ! Sn
    43.67197d0, 459.322d0, 4.276d0, &  ! Sb
    37.65d0, 396.0d0, 4.22d0, &  ! Te
    35.0d0, 385.0d0, 4.17d0, &  ! I
    27.3d0, 285.9d0, 4.08d0, &  ! Xe
    427.12d0, 6582.08d0, 3.78d0, &  ! Cs
    275.0d0, 5727.0d0, 4.77d0, &  ! Ba
    0d0, 0d0, 0d0, &  ! La
    0d0, 0d0, 0d0, &  ! Ce
    0d0, 0d0, 0d0, &  ! Pr
    0d0, 0d0, 0d0, &  ! Nd
    0d0, 0d0, 0d0, &  ! Pm
    0d0, 0d0, 0d0, &  ! Sm
    0d0, 0d0, 0d0, &  ! Eu
    0d0, 0d0, 0d0, &  ! Gd
    0d0, 0d0, 0d0, &  ! Tb
    0d0, 0d0, 0d0, &  ! Dy
    0d0, 0d0, 0d0, &  ! Ho
    0d0, 0d0, 0d0, &  ! Er
    0d0, 0d0, 0d0, &  ! Tm
    0d0, 0d0, 0d0, &  ! Yb
    0d0, 0d0, 0d0, &  ! Lu
    99.52d0, 1274.8d0, 4.21d0, &  ! Hf
    82.53d0, 1019.92d0, 4.15d0, &  ! Ta
    71.041d0, 847.93d0, 4.08d0, &  ! W
    63.04d0, 710.2d0, 4.02d0, &  ! Re
    55.055d0, 596.67d0, 3.84d0, &  ! Os
    42.51d0, 359.1d0, 4.0d0, &  ! Ir
    39.68d0, 347.1d0, 3.92d0, &  ! Pt
    36.5d0, 298.0d0, 3.86d0, &  ! Au
    33.9d0, 392.0d0, 3.98d0, &  ! Hg
    69.92d0, 717.44d0, 3.91d0, &  ! Tl
    61.8d0, 697.0d0, 4.31d0, &  ! Pb
    49.02d0, 571.0d0, 4.32d0, &  ! Bi
    45.013d0, 530.92d0, 4.097d0, &  ! Po
    38.93d0, 457.53d0, 4.07d0, &  ! At
    33.54d0, 390.63d0, 4.23d0 &  ! Rn
], [3, 86])

real(dp), parameter :: tssurf_vdw_params(3, 86) = reshape([ &
    4.5d0, 6.5d0, 3.1d0, &  ! H
    1.38d0, 1.46d0, 2.65d0, &  ! He
    164.2d0, 1387.0d0, 4.16d0, &  ! Li
    38.0d0, 214.0d0, 4.17d0, &  ! Be
    21.0d0, 99.5d0, 3.89d0, &  ! B
    12.0d0, 46.6d0, 3.59d0, &  ! C
    7.4d0, 24.2d0, 3.34d0, &  ! N
    5.4d0, 15.6d0, 3.19d0, &  ! O
    3.8d0, 9.52d0, 3.04d0, &  ! F
    2.67d0, 6.38d0, 2.91d0, &  ! Ne
    162.7d0, 1556.0d0, 3.73d0, &  ! Na
    71.0d0, 627.0d0, 4.27d0, &  ! Mg
    60.0d0, 528.0d0, 4.33d0, &  ! Al
    37.0d0, 305.0d0, 4.2d0, &  ! Si
    25.0d0, 185.0d0, 4.01d0, &  ! P
    19.6d0, 134.0d0, 3.86d0, &  ! S
    15.0d0, 94.6d0, 3.71d0, &  ! Cl
    11.1d0, 64.3d0, 3.55d0, &  ! Ar
    292.9d0, 3897.0d0, 3.71d0, &  ! K
    160.0d0, 2221.0d0, 4.65d0, &  ! Ca
    120.0d0, 1383.0d0, 4.59d0, &  ! Sc
    98.0d0, 1044.0d0, 4.51d0, &  ! Ti
    84.0d0, 832.0d0, 4.44d0, &  ! V
    78.0d0, 602.0d0, 3.99d0, &  ! Cr
    63.0d0, 552.0d0, 3.97d0, &  ! Mn
    56.0d0, 482.0d0, 4.23d0, &  ! Fe
    50.0d0, 408.0d0, 4.18d0, &  ! Co
    10.22d0, 59.2d0, 2.28d0, &  ! Ni
    10.88d0, 58.9d0, 2.4d0, &  ! Cu
    13.77d0, 46.0d0, 2.82d0, &  ! Zn
    60.0d0, 498.0d0, 4.19d0, &  ! Ga
    41.0d0, 354.0d0, 4.2d0, &  ! Ge
    29.0d0, 246.0d0, 4.11d0, &  ! As
    25.0d0, 210.0d0, 4.04d0, &  ! Se
    20.0d0, 162.0d0, 3.93d0, &  ! Br
    16.8d0, 129.6d0, 3.82d0, &  ! Kr
    319.2d0, 4691.0d0, 3.72d0, &  ! Rb
    199.0d0, 3170.0d0, 4.54d0, &  ! Sr
    126.737d0, 1968.58d0, 4.8151d0, &  ! Y
    119.97d0, 1677.91d0, 4.53d0, &  ! Zr
    101.603d0, 1263.61d0, 4.2365d0, &  ! Nb
    88.4225785d0, 1028.73d0, 4.099d0, &  ! Mo
    80.083d0, 1390.87d0, 4.076d0, &  ! Tc
    65.895d0, 609.754d0, 3.9953d0, &  ! Ru
    56.1d0, 469.0d0, 3.95d0, &  ! Rh
    13.9d0, 102.0d0, 3.06d0, &  ! Pd
    15.36d0, 122.0d0, 2.57d0, &  ! Ag
    39.7d0, 452.0d0, 3.99d0, &  ! Cd
    70.22d0, 707.046d0, 4.23198d0, &  ! In
    55.95d0, 587.417d0, 4.303d0, &  ! Sn
    43.67197d0, 459.322d0, 4.276d0, &  ! Sb
    37.65d0, 396.0d0, 4.22d0, &  ! Te
    35.0d0, 385.0d0, 4.17d0, &  ! I
    27.3d0, 285.9d0, 4.08d0, &  ! Xe
    427.12d0, 6582.08d0, 3.78d0, &  ! Cs
    275.0d0, 5727.0d0, 4.77d0, &  ! Ba
    0d0, 0d0, 0d0, &  ! La
    0d0, 0d0, 0d0, &  ! Ce
    0d0, 0d0, 0d0, &  ! Pr
    0d0, 0d0, 0d0, &  ! Nd
    0d0, 0d0, 0d0, &  ! Pm
    0d0, 0d0, 0d0, &  ! Sm
    0d0, 0d0, 0d0, &  ! Eu
    0d0, 0d0, 0d0, &  ! Gd
    0d0, 0d0, 0d0, &  ! Tb
    0d0, 0d0, 0d0, &  ! Dy
    0d0, 0d0, 0d0, &  ! Ho
    0d0, 0d0, 0d0, &  ! Er
    0d0, 0d0, 0d0, &  ! Tm
    0d0, 0d0, 0d0, &  ! Yb
    0d0, 0d0, 0d0, &  ! Lu
    99.52d0, 1274.8d0, 4.21d0, &  ! Hf
    82.53d0, 1019.92d0, 4.15d0, &  ! Ta
    71.041d0, 847.93d0, 4.08d0, &  ! W
    63.04d0, 710.2d0, 4.02d0, &  ! Re
    55.055d0, 596.67d0, 3.84d0, &  ! Os
    42.51d0, 359.1d0, 4.0d0, &  ! Ir
    14.45d0, 120.5d0, 2.8d0, &  ! Pt
    15.62d0, 133.9d0, 2.91d0, &  ! Au
    33.9d0, 392.0d0, 3.98d0, &  ! Hg
    69.92d0, 717.44d0, 3.91d0, &  ! Tl
    61.8d0, 697.0d0, 4.31d0, &  ! Pb
    49.02d0, 571.0d0, 4.32d0, &  ! Bi
    45.013d0, 530.92d0, 4.097d0, &  ! Po
    38.93d0, 457.53d0, 4.07d0, &  ! At
    33.54d0, 390.63d0, 4.23d0 &  ! Rn
], [3, 86])

contains


integer elemental function species_index(species)
    character(len=*), intent(in) :: species

    integer :: i

    select case (lower(trim(species)))
    case ('h');  i = 1;   case ('he'); i = 2;   case ('li'); i = 3
    case ('be'); i = 4;   case ('b');  i = 5;   case ('c');  i = 6
    case ('n');  i = 7;   case ('o');  i = 8;   case ('f');  i = 9
    case ('ne'); i = 10;  case ('na'); i = 11;  case ('mg'); i = 12
    case ('al'); i = 13;  case ('si'); i = 14;  case ('p');  i = 15
    case ('s');  i = 16;  case ('cl'); i = 17;  case ('ar'); i = 18
    case ('k');  i = 19;  case ('ca'); i = 20;  case ('sc'); i = 21
    case ('ti'); i = 22;  case ('v');  i = 23;  case ('cr'); i = 24
    case ('mn'); i = 25;  case ('fe'); i = 26;  case ('co'); i = 27
    case ('ni'); i = 28;  case ('cu'); i = 29;  case ('zn'); i = 30
    case ('ga'); i = 31;  case ('ge'); i = 32;  case ('as'); i = 33
    case ('se'); i = 34;  case ('br'); i = 35;  case ('kr'); i = 36
    case ('rb'); i = 37;  case ('sr'); i = 38;  case ('y');  i = 39
    case ('zr'); i = 40;  case ('nb'); i = 41;  case ('mo'); i = 42
    case ('tc'); i = 43;  case ('ru'); i = 44;  case ('rh'); i = 45
    case ('pd'); i = 46;  case ('ag'); i = 47;  case ('cd'); i = 48
    case ('in'); i = 49;  case ('sn'); i = 50;  case ('sb'); i = 51
    case ('te'); i = 52;  case ('i');  i = 53;  case ('xe'); i = 54
    case ('cs'); i = 55;  case ('ba'); i = 56;  case ('la'); i = 57
    case ('ce'); i = 58;  case ('pr'); i = 59;  case ('nd'); i = 60
    case ('pm'); i = 61;  case ('sm'); i = 62;  case ('eu'); i = 63
    case ('gd'); i = 64;  case ('tb'); i = 65;  case ('dy'); i = 66
    case ('ho'); i = 67;  case ('er'); i = 68;  case ('tm'); i = 69
    case ('yb'); i = 70;  case ('lu'); i = 71;  case ('hf'); i = 72
    case ('ta'); i = 73;  case ('w');  i = 74;  case ('re'); i = 75
    case ('os'); i = 76;  case ('ir'); i = 77;  case ('pt'); i = 78
    case ('au'); i = 79;  case ('hg'); i = 80;  case ('tl'); i = 81
    case ('pb'); i = 82;  case ('bi'); i = 83;  case ('po'); i = 84
    case ('at'); i = 85;  case ('rn'); i = 86;  case ('fr'); i = 87
    case ('ra'); i = 88;  case ('ac'); i = 89;  case ('th'); i = 90
    case ('pa'); i = 91;  case ('u');  i = 92;  case ('np'); i = 93
    case ('pu'); i = 94;  case ('am'); i = 95;  case ('cm'); i = 96
    case ('bk'); i = 97;  case ('cf'); i = 98;  case ('es'); i = 99
    case ('fm'); i = 100; case ('md'); i = 101; case ('no'); i = 102
    case ('lr'); i = 103; case ('rf'); i = 104; case ('db'); i = 105
    case ('sg'); i = 106; case ('bh'); i = 107; case ('hs'); i = 108
    case ('mt'); i = 109; case ('ds'); i = 110; case ('rg'); i = 111
    case ('cn'); i = 112
    case default
        i = -1
    end select
    species_index = i
end function species_index

end module
