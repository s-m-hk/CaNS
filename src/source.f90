!
! SPDX-License-Identifier: MIT
!
module mod_source
  !
  use mod_types
  use mod_common_mpi, only: myid,ierr
#if defined(_BOUSSINESQ)
  use mod_param     , only: gacc,tmp0,beta_th
#else
  use mod_param     , only: gacc
#endif
  implicit none
  private
  public :: grav_src
  contains
  subroutine grav_src(nx,ny,nz, &
                      factor12, &
#if defined(_HEAT_TRANSFER)
                      nh_s,tmp, &
#endif
                      u,v,w)
    implicit none
    integer , intent(in   )                      :: nx,ny,nz
    real(rp), intent(in   )                      :: factor12
#if defined(_HEAT_TRANSFER)
    integer , intent(in   )                      :: nh_s
    real(rp), intent(in   ), dimension(1-nh_s:,1-nh_s:,1-nh_s:) :: tmp
#endif
    real(rp), intent(inout), dimension(0:,0:,0:) :: u,v,w
    !
    real(rp) :: gacc1, gacc2, gacc3
    real(rp) :: termx,termy,termz,tmppx,tmppy,tmppz
    integer  :: i,j,k
    !
    gacc1 = gacc(1)
    gacc2 = gacc(2)
    gacc3 = gacc(3)
    !
    !$acc parallel loop collapse(3) default(present) private(termx,termy,termz,tmppx,tmppy,tmppz) async(1)
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i,j,k) &
    !$OMP PRIVATE(tmppx,tmppy,tmppz) &
    !$OMP PRIVATE(termx,termy,termz) &
    !$OMP SHARED(nx,ny,nz,dudt,dvdt,dwdt)
    do k=1,nz
      do j=1,ny
        do i=1,nx
#if defined(_HEAT_TRANSFER) && defined(_BOUSSINESQ)
          tmppx = 0.5_rp*(tmp(i+1,j,k)+tmp(i,j,k))
          tmppy = 0.5_rp*(tmp(i,j+1,k)+tmp(i,j,k))
          tmppz = 0.5_rp*(tmp(i,j,k+1)+tmp(i,j,k))
          termx = gacc1*(1.0_rp-beta_th*(tmppx-tmp0))
          termy = gacc2*(1.0_rp-beta_th*(tmppy-tmp0))
          termz = gacc3*(1.0_rp-beta_th*(tmppz-tmp0))
#else
          termx = gacc1
          termy = gacc2
          termz = gacc3
#endif
          !
          u(i,j,k) = u(i,j,k) + factor12*termx
          v(i,j,k) = v(i,j,k) + factor12*termy
          w(i,j,k) = w(i,j,k) + factor12*termz
          !
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO
  end subroutine grav_src
  !
end module mod_source
