! -
!
! SPDX-FileCopyrightText: Copyright (c) 2017-2022 Pedro Costa and the CaNS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
module mod_common_mpi
  implicit none
  public
  integer :: myid,ierr
  integer :: halo(3),halo_s(3),halo_big(3)
  integer :: nh_p,nh_v,nh_s,nh_b
  integer :: ipencil_axis
end module mod_common_mpi
