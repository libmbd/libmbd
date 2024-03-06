! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
#include "version.h"

module mbd_version

implicit none

integer, parameter, public :: mbd_version_major = MBD_VERSION_MAJOR
integer, parameter, public :: mbd_version_minor = MBD_VERSION_MINOR
integer, parameter, public :: mbd_version_patch = MBD_VERSION_PATCH
character(len=30), parameter, public :: mbd_version_suffix = MBD_VERSION_SUFFIX

end module
