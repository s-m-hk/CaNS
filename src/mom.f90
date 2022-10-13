! -
!
! SPDX-FileCopyrightText: Copyright (c) 2017-2022 Pedro Costa and the CaNS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
module mod_mom
  use mpi
  use mod_common_mpi, only: ierr
  use mod_types
  implicit none
  private
  public momx_a,momy_a,momz_a, &
         momx_d,momy_d,momz_d, &
         momx_p,momy_p,momz_p, cmpt_wallshear
#if defined(_IMPDIFF_1D)
  public momx_d_xy,momy_d_xy,momz_d_xy, &
         momx_d_z ,momy_d_z ,momz_d_z
#endif
  public mom_xyz_ad
  contains
  !
  subroutine momx_a(nx,ny,nz,dxi,dyi,dzfi,u,v,w, &
#if defined(_IBM)
                    psi_u,psi_v,psi_w, &
#endif
                    dudt)
    implicit none
    integer , intent(in) :: nx,ny,nz
    real(rp), intent(in) :: dxi,dyi
    real(rp), intent(in), dimension(0:) :: dzfi
    real(rp), dimension(0:,0:,0:), intent(in   ) :: u,v,w
#if defined(_IBM)
    real(rp), dimension(0:,0:,0:), intent(in   ) :: psi_u,psi_v,psi_w
#endif
    real(rp), dimension( :, :, :), intent(inout) :: dudt
    integer :: i,j,k
    real(rp) :: uuip,uuim,uvjp,uvjm,uwkp,uwkm
    !
    !$acc parallel loop collapse(3) default(present) private(uuip,uuim,uvjp,uvjm,uwkp,uwkm) async(1)
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(uuip,uuim,uvjp,uvjm,uwkp,uwkm) &
    !$OMP SHARED(nx,ny,nz,dxi,dyi,u,v,w,psi_u,psi_v,psi_w,dudt,dzfi)
    do k=1,nz
      do j=1,ny
        do i=1,nx
          uuip  = 0.25*( u(i+1,j,k)+u(i,j,k) )*( u(i+1,j  ,k  )+u(i,j  ,k  ) )
          uuim  = 0.25*( u(i-1,j,k)+u(i,j,k) )*( u(i-1,j  ,k  )+u(i,j  ,k  ) )
          uvjp  = 0.25*( u(i,j+1,k)+u(i,j,k) )*( v(i+1,j  ,k  )+v(i,j  ,k  ) )
          uvjm  = 0.25*( u(i,j-1,k)+u(i,j,k) )*( v(i+1,j-1,k  )+v(i,j-1,k  ) )
          uwkp  = 0.25*( u(i,j,k+1)+u(i,j,k) )*( w(i+1,j  ,k  )+w(i,j  ,k  ) )
          uwkm  = 0.25*( u(i,j,k-1)+u(i,j,k) )*( w(i+1,j  ,k-1)+w(i,j  ,k-1) )
#if defined(_IBM) && defined(_SIMPLE)
          if (psi_u(i,j,k) == 0.) then ! if not in solid
            if (psi_u(i,j+1,k) == 1.) uvjp  = 0._rp
            if (psi_u(i,j-1,k) == 1.) uvjm  = 0._rp
            if (psi_u(i,j,k+1) == 1.) uwkp  = 0._rp
            if (psi_u(i,j,k-1) == 1.) uwkm  = 0._rp
          endif
#endif
          !
          ! Momentum balance
          !
          dudt(i,j,k) = dudt(i,j,k) + &
                        dxi*(     -uuip + uuim ) + &
                        dyi*(     -uvjp + uvjm ) + &
                        dzfi(k)*( -uwkp + uwkm )
        end do
      end do
    end do
  end subroutine momx_a
  !
  subroutine momy_a(nx,ny,nz,dxi,dyi,dzfi,u,v,w, &
#if defined(_IBM)
                    psi_u,psi_v,psi_w, &
#endif
                    dvdt)
    implicit none
    integer , intent(in) :: nx,ny,nz
    real(rp), intent(in) :: dxi,dyi
    real(rp), intent(in), dimension(0:) :: dzfi
    real(rp), dimension(0:,0:,0:), intent(in   ) :: u,v,w
#if defined(_IBM)
    real(rp), dimension(0:,0:,0:), intent(in   ) :: psi_u,psi_v,psi_w
#endif
    real(rp), dimension( :, :, :), intent(inout) :: dvdt
    integer :: i,j,k
    real(rp) :: uvip,uvim,vvjp,vvjm,wvkp,wvkm
    !
    !$acc parallel loop collapse(3) default(present) private(uvip,uvim,vvjp,vvjm,wvkp,wvkm) async(1)
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(uvip,uvim,vvjp,vvjm,wvkp,wvkm) &
    !$OMP SHARED(nx,ny,nz,dxi,dyi,dzfi,u,v,w,psi_u,psi_v,psi_w,dvdt)
    do k=1,nz
      do j=1,ny
        do i=1,nx
          uvip  = 0.25*( u(i  ,j,k  )+u(i  ,j+1,k  ) )*( v(i,j,k)+v(i+1,j,k) )
          uvim  = 0.25*( u(i-1,j,k  )+u(i-1,j+1,k  ) )*( v(i,j,k)+v(i-1,j,k) )
          vvjp  = 0.25*( v(i  ,j,k  )+v(i  ,j+1,k  ) )*( v(i,j,k)+v(i,j+1,k) )
          vvjm  = 0.25*( v(i  ,j,k  )+v(i  ,j-1,k  ) )*( v(i,j,k)+v(i,j-1,k) )
          wvkp  = 0.25*( w(i  ,j,k  )+w(i  ,j+1,k  ) )*( v(i,j,k+1)+v(i,j,k) )
          wvkm  = 0.25*( w(i  ,j,k-1)+w(i  ,j+1,k-1) )*( v(i,j,k-1)+v(i,j,k) )
#if defined(_IBM) && defined(_SIMPLE)
          if (psi_v(i,j,k) == 0.) then ! if not in solid
            if (psi_v(i+1,j,k) == 1.) uvip = 0._rp
            if (psi_v(i-1,j,k) == 1.) uvim = 0._rp
            if (psi_v(i,j,k+1) == 1.) wvkp = 0._rp
            if (psi_v(i,j,k-1) == 1.) wvkm = 0._rp
          endif
#endif
          !
          ! Momentum balance
          !
          dvdt(i,j,k) = dvdt(i,j,k) + &
                        dxi*(     -uvip + uvim ) + &
                        dyi*(     -vvjp + vvjm ) + &
                        dzfi(k)*( -wvkp + wvkm )
        end do
      end do
    end do
  end subroutine momy_a
  !
  subroutine momz_a(nx,ny,nz,dxi,dyi,dzci,u,v,w, &
#if defined(_IBM)
                    psi_u,psi_v,psi_w, &
#endif
                    dwdt)
    implicit none
    integer , intent(in) :: nx,ny,nz
    real(rp), intent(in) :: dxi,dyi
    real(rp), intent(in), dimension(0:) :: dzci
    real(rp), dimension(0:,0:,0:), intent(in   ) :: u,v,w
#if defined(_IBM)
    real(rp), dimension(0:,0:,0:), intent(in   ) :: psi_u,psi_v,psi_w
#endif
    real(rp), dimension( :, :, :), intent(inout) :: dwdt
    integer :: i,j,k
    real(rp) :: uwip,uwim,vwjp,vwjm,wwkp,wwkm
    !
    !$acc parallel loop collapse(3) default(present) private(uwip,uwim,vwjp,vwjm,wwkp,wwkm) async(1)
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(uwip,uwim,vwjp,vwjm,wwkp,wwkm) &
    !$OMP SHARED(nx,ny,nz,dxi,dyi,dzci,u,v,w,psi_u,psi_v,psi_w,dwdt)
    do k=1,nz
      do j=1,ny
        do i=1,nx
          uwip  = 0.25*( w(i,j,k)+w(i+1,j,k) )*( u(i  ,j  ,k)+u(i  ,j  ,k+1) )
          uwim  = 0.25*( w(i,j,k)+w(i-1,j,k) )*( u(i-1,j  ,k)+u(i-1,j  ,k+1) )
          vwjp  = 0.25*( w(i,j,k)+w(i,j+1,k) )*( v(i  ,j  ,k)+v(i  ,j  ,k+1) )
          vwjm  = 0.25*( w(i,j,k)+w(i,j-1,k) )*( v(i  ,j-1,k)+v(i  ,j-1,k+1) )
          wwkp  = 0.25*( w(i,j,k)+w(i,j,k+1) )*( w(i  ,j  ,k)+w(i  ,j  ,k+1) )
          wwkm  = 0.25*( w(i,j,k)+w(i,j,k-1) )*( w(i  ,j  ,k)+w(i  ,j  ,k-1) )
#if defined(_IBM) && defined(_SIMPLE)
          if (psi_w(i,j,k) == 0.) then ! if not in solid
           if (psi_w(i+1,j,k) == 1.) uwip = 0._rp
           if (psi_w(i-1,j,k) == 1.) uwim = 0._rp
           if (psi_w(i,j+1,k) == 1.) vwjp = 0._rp
           if (psi_w(i,j-1,k) == 1.) vwjm = 0._rp
          endif
#endif
          !
          ! Momentum balance
          !
          dwdt(i,j,k) = dwdt(i,j,k) + &
                        dxi*(     -uwip + uwim ) + &
                        dyi*(     -vwjp + vwjm ) + &
                        dzci(k)*( -wwkp + wwkm )
        end do
      end do
    end do
  end subroutine momz_a
  !
  subroutine momx_d(nx,ny,nz,dxi,dyi,dzci,dzfi,visc,u, &
#if defined(_IBM)
                    psi_u,psi_v,psi_w, &
#endif
                    dudt)
    implicit none
    integer , intent(in) :: nx,ny,nz
    real(rp), intent(in) :: dxi,dyi,visc
    real(rp), intent(in), dimension(0:) :: dzci,dzfi
    real(rp), dimension(0:,0:,0:), intent(in   ) :: u
#if defined(_IBM)
    real(rp), dimension(0:,0:,0:), intent(in   ) :: psi_u,psi_v,psi_w
#endif
    real(rp), dimension( :, :, :), intent(inout) :: dudt
    real(rp) :: dudxp,dudxm,dudyp,dudym,dudzp,dudzm
    integer :: i,j,k
    !
    !$acc parallel loop collapse(3) default(present) private(dudxp,dudxm,dudyp,dudym,dudzp,dudzm) async(1)
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(dudxp,dudxm,dudyp,dudym,dudzp,dudzm) &
    !$OMP SHARED(nx,ny,nz,dxi,dyi,dzci,dzfi,u,psi_u,psi_v,psi_w,dudt,visc)
    do k=1,nz
      do j=1,ny
        do i=1,nx
          dudxp = (u(i+1,j,k)-u(i,j,k))*dxi
          dudxm = (u(i,j,k)-u(i-1,j,k))*dxi
          dudyp = (u(i,j+1,k)-u(i,j,k))*dyi
          dudym = (u(i,j,k)-u(i,j-1,k))*dyi
          dudzp = (u(i,j,k+1)-u(i,j,k))*dzci(k  )
          dudzm = (u(i,j,k)-u(i,j,k-1))*dzci(k-1)
#if defined(_IBM) && defined(_SIMPLE)
          if (psi_u(i,j,k) == 0.) then ! if not in solid
            if (psi_u(i,j+1,k) == 1.) dudyp = -2._rp*u(i,j,k)*dyi
            if (psi_u(i,j-1,k) == 1.) dudym =  2._rp*u(i,j,k)*dyi
            if (psi_u(i,j,k+1) == 1.) dudzp = -2._rp*u(i,j,k)*dzci(k)
            if (psi_u(i,j,k-1) == 1.) dudzm =  2._rp*u(i,j,k)*dzci(k)
          endif
#endif
          dudt(i,j,k) = dudt(i,j,k) + &
                        (dudxp-dudxm)*visc*dxi + &
                        (dudyp-dudym)*visc*dyi + &
                        (dudzp-dudzm)*visc*dzfi(k)
        end do
      end do
    end do
  end subroutine momx_d
  !
  subroutine momy_d(nx,ny,nz,dxi,dyi,dzci,dzfi,visc,v, &
#if defined(_IBM)
                    psi_u,psi_v,psi_w, &
#endif
                    dvdt)
    implicit none
    integer , intent(in) :: nx,ny,nz
    real(rp), intent(in) :: dxi,dyi,visc
    real(rp), intent(in), dimension(0:) :: dzci,dzfi
    real(rp), dimension(0:,0:,0:), intent(in   ) :: v
#if defined(_IBM)
    real(rp), dimension(0:,0:,0:), intent(in   ) :: psi_u,psi_v,psi_w
#endif
    real(rp), dimension( :, :, :), intent(inout) :: dvdt
    real(rp) :: dvdxp,dvdxm,dvdyp,dvdym,dvdzp,dvdzm
    integer :: i,j,k
    !
    !$acc parallel loop collapse(3) default(present) private(dvdxp,dvdxm,dvdyp,dvdym,dvdzp,dvdzm) async(1)
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(dvdxp,dvdxm,dvdyp,dvdym,dvdzp,dvdzm) &
    !$OMP SHARED(nx,ny,nz,dxi,dyi,dzci,dzfi,v,psi_u,psi_v,psi_w,dvdt,visc)
    do k=1,nz
      do j=1,ny
        do i=1,nx
          dvdxp = (v(i+1,j,k)-v(i,j,k))*dxi
          dvdxm = (v(i,j,k)-v(i-1,j,k))*dxi
          dvdyp = (v(i,j+1,k)-v(i,j,k))*dyi
          dvdym = (v(i,j,k)-v(i,j-1,k))*dyi
          dvdzp = (v(i,j,k+1)-v(i,j,k))*dzci(k  )
          dvdzm = (v(i,j,k)-v(i,j,k-1))*dzci(k-1)
#if defined(_IBM) && defined(_SIMPLE)
          if (psi_v(i,j,k) == 0.) then ! if not in solid
            if (psi_v(i+1,j,k) == 1.) dvdxp = -2._rp*v(i,j,k)*dxi
            if (psi_v(i-1,j,k) == 1.) dvdxm =  2._rp*v(i,j,k)*dxi
            if (psi_v(i,j,k+1) == 1.) dvdzp = -2._rp*v(i,j,k)*dzci(k)
            if (psi_v(i,j,k-1) == 1.) dvdzm =  2._rp*v(i,j,k)*dzci(k)
          endif
#endif
          dvdt(i,j,k) = dvdt(i,j,k) + &
                        (dvdxp-dvdxm)*visc*dxi + &
                        (dvdyp-dvdym)*visc*dyi + &
                        (dvdzp-dvdzm)*visc*dzfi(k)
        end do
      end do
    end do
  end subroutine momy_d
  !
  subroutine momz_d(nx,ny,nz,dxi,dyi,dzci,dzfi,visc,w, &
#if defined(_IBM)
                    psi_u,psi_v,psi_w, &
#endif
                    dwdt)
    implicit none
    integer , intent(in) :: nx,ny,nz
    real(rp), intent(in) :: dxi,dyi,visc
    real(rp), intent(in), dimension(0:) :: dzci,dzfi
    real(rp), dimension(0:,0:,0:), intent(in   ) :: w
#if defined(_IBM)
    real(rp), dimension(0:,0:,0:), intent(in   ) :: psi_u,psi_v,psi_w
#endif
    real(rp), dimension( :, :, :), intent(inout) :: dwdt
    integer :: i,j,k
    real(rp) :: dwdxp,dwdxm,dwdyp,dwdym,dwdzp,dwdzm
    !
    !$acc parallel loop collapse(3) default(present) private(dwdxp,dwdxm,dwdyp,dwdym,dwdzp,dwdzm) async(1)
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(dwdxp,dwdxm,dwdyp,dwdym,dwdzp,dwdzm) &
    !$OMP SHARED(nx,ny,nz,dxi,dyi,dzci,dzfi,w,psi_u,psi_v,psi_w,dwdt,visc)
    do k=1,nz
      do j=1,ny
        do i=1,nx
          dwdxp = (w(i+1,j,k)-w(i,j,k))*dxi
          dwdxm = (w(i,j,k)-w(i-1,j,k))*dxi
          dwdyp = (w(i,j+1,k)-w(i,j,k))*dyi
          dwdym = (w(i,j,k)-w(i,j-1,k))*dyi
          dwdzp = (w(i,j,k+1)-w(i,j,k))*dzfi(k+1)
          dwdzm = (w(i,j,k)-w(i,j,k-1))*dzfi(k  )
#if defined(_IBM) && defined(_SIMPLE)
          if (psi_w(i,j,k) == 0.) then ! if not in solid
           if (psi_w(i+1,j,k) == 1.) dwdxp = -2._rp*w(i,j,k)*dxi
           if (psi_w(i-1,j,k) == 1.) dwdxm =  2._rp*w(i,j,k)*dxi
           if (psi_w(i,j+1,k) == 1.) dwdyp = -2._rp*w(i,j,k)*dyi
           if (psi_w(i,j-1,k) == 1.) dwdym =  2._rp*w(i,j,k)*dyi
          endif
#endif
          dwdt(i,j,k) = dwdt(i,j,k) + &
                        (dwdxp-dwdxm)*visc*dxi + &
                        (dwdyp-dwdym)*visc*dyi + &
                        (dwdzp-dwdzm)*visc*dzci(k)
        end do
      end do
    end do
  end subroutine momz_d
  !
  subroutine momx_p(nx,ny,nz,dxi,bforce,p,dudt)
    implicit none
    integer , intent(in) :: nx,ny,nz
    real(rp), intent(in) :: dxi
    real(rp), intent(in) :: bforce
    real(rp), dimension(0:,0:,0:), intent(in ) :: p
    real(rp), dimension( :, :, :), intent(out) :: dudt
    integer :: i,j,k
    !
    !$acc parallel loop collapse(3) default(present) async(1)
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP SHARED(nx,ny,nz,dxi,bforce,p,dudt)
    do k=1,nz
      do j=1,ny
        do i=1,nx
          dudt(i,j,k) = - dxi*( p(i+1,j,k)-p(i,j,k) ) + bforce
        end do
      end do
    end do
  end subroutine momx_p
  !
  subroutine momy_p(nx,ny,nz,dyi,bforce,p,dvdt)
    implicit none
    integer , intent(in) :: nx,ny,nz
    real(rp), intent(in) :: dyi
    real(rp), intent(in) :: bforce
    real(rp), dimension(0:,0:,0:), intent(in ) :: p
    real(rp), dimension( :, :, :), intent(out) :: dvdt
    integer :: i,j,k
    !
    !$acc parallel loop collapse(3) default(present) async(1)
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP SHARED(nx,ny,nz,dyi,bforce,p,dvdt)
    do k=1,nz
      do j=1,ny
        do i=1,nx
          dvdt(i,j,k) = - dyi*( p(i,j+1,k)-p(i,j,k) ) + bforce
        end do
      end do
    end do
  end subroutine momy_p
  !
  subroutine momz_p(nx,ny,nz,dzci,bforce,p,dwdt)
    implicit none
    integer , intent(in) :: nx,ny,nz
    real(rp), intent(in), dimension(0:) :: dzci
    real(rp), intent(in) :: bforce
    real(rp), dimension(0:,0:,0:), intent(in ) :: p
    real(rp), dimension( :, :, :), intent(out) :: dwdt
    integer :: i,j,k
    !
    !$acc parallel loop collapse(3) default(present) async(1)
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP SHARED(nx,ny,nz,dzci,bforce,p,dwdt)
    do k=1,nz
      do j=1,ny
        do i=1,nx
          dwdt(i,j,k) = - dzci(k)*( p(i,j,k+1)-p(i,j,k) ) + bforce
        end do
      end do
    end do
  end subroutine momz_p
  !
#if defined(_IMPDIFF_1D)
  subroutine momx_d_z(nx,ny,nz,dzci,dzfi,visc,u,dudt)
    implicit none
    integer , intent(in) :: nx,ny,nz
    real(rp), intent(in) :: visc
    real(rp), intent(in), dimension(0:) :: dzci,dzfi
    real(rp), dimension(0:,0:,0:), intent(in   ) :: u
    real(rp), dimension( :, :, :), intent(inout) :: dudt
    real(rp) :: dudzp,dudzm
    integer :: i,j,k
    !
    !$acc parallel loop collapse(3) default(present) private(dudzp,dudzm) async(1)
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(dudzp,dudzm) &
    !$OMP SHARED(nx,ny,nz,dzci,dzfi,u,dudt,visc)
    do k=1,nz
      do j=1,ny
        do i=1,nx
          dudzp = (u(i,j,k+1)-u(i,j,k))*dzci(k  )
          dudzm = (u(i,j,k)-u(i,j,k-1))*dzci(k-1)
          dudt(i,j,k) = dudt(i,j,k) + &
                        (dudzp-dudzm)*visc*dzfi(k)
        end do
      end do
    end do
  end subroutine momx_d_z
  !
  subroutine momy_d_z(nx,ny,nz,dzci,dzfi,visc,v,dvdt)
    implicit none
    integer , intent(in) :: nx,ny,nz
    real(rp), intent(in) :: visc
    real(rp), intent(in), dimension(0:) :: dzci,dzfi
    real(rp), dimension(0:,0:,0:), intent(in   ) :: v
    real(rp), dimension( :, :, :), intent(inout) :: dvdt
    real(rp) :: dvdzp,dvdzm
    integer :: i,j,k
    !
    !$acc parallel loop collapse(3) default(present) private(dvdzp,dvdzm) async(1)
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(dvdzp,dvdzm) &
    !$OMP SHARED(nx,ny,nz,dzci,dzfi,v,dvdt,visc)
    do k=1,nz
      do j=1,ny
        do i=1,nx
          dvdzp = (v(i,j,k+1)-v(i,j,k))*dzci(k  )
          dvdzm = (v(i,j,k)-v(i,j,k-1))*dzci(k-1)
          dvdt(i,j,k) = dvdt(i,j,k) + &
                       (dvdzp-dvdzm)*visc*dzfi(k)
        end do
      end do
    end do
  end subroutine momy_d_z
  !
  subroutine momz_d_z(nx,ny,nz,dzci,dzfi,visc,w,dwdt)
    implicit none
    integer , intent(in) :: nx,ny,nz
    real(rp), intent(in) :: visc
    real(rp), intent(in), dimension(0:) :: dzci,dzfi
    real(rp), dimension(0:,0:,0:), intent(in   ) :: w
    real(rp), dimension( :, :, :), intent(inout) :: dwdt
    integer :: i,j,k
    real(rp) :: dwdzp,dwdzm
    !
    !$acc parallel loop collapse(3) default(present) private(dwdzp,dwdzm) async(1)
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(dwdzp,dwdzm) &
    !$OMP SHARED(nx,ny,nz,dzci,dzfi,w,dwdt,visc)
    do k=1,nz
      do j=1,ny
        do i=1,nx
          dwdzp = (w(i,j,k+1)-w(i,j,k))*dzfi(k+1)
          dwdzm = (w(i,j,k)-w(i,j,k-1))*dzfi(k  )
          dwdt(i,j,k) = dwdt(i,j,k) + &
                        (dwdzp-dwdzm)*visc*dzci(k)
        end do
      end do
    end do
  end subroutine momz_d_z
  !
  subroutine momx_d_xy(nx,ny,nz,dxi,dyi,visc,u,dudt)
    implicit none
    integer , intent(in) :: nx,ny,nz
    real(rp), intent(in) :: dxi,dyi,visc
    real(rp), dimension(0:,0:,0:), intent(in   ) :: u
    real(rp), dimension( :, :, :), intent(inout) :: dudt
    real(rp) :: dudxp,dudxm,dudyp,dudym
    integer :: i,j,k
    !
    !$acc parallel loop collapse(3) default(present) private(dudxp,dudxm,dudyp,dudym) async(1)
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(dudxp,dudxm,dudyp,dudym) &
    !$OMP SHARED(nx,ny,nz,dxi,dyi,u,dudt,visc)
    do k=1,nz
      do j=1,ny
        do i=1,nx
          dudxp = (u(i+1,j,k)-u(i,j,k))*dxi
          dudxm = (u(i,j,k)-u(i-1,j,k))*dxi
          dudyp = (u(i,j+1,k)-u(i,j,k))*dyi
          dudym = (u(i,j,k)-u(i,j-1,k))*dyi
          dudt(i,j,k) = dudt(i,j,k) + &
                        (dudxp-dudxm)*visc*dxi + &
                        (dudyp-dudym)*visc*dyi
        end do
      end do
    end do
  end subroutine momx_d_xy
  !
  subroutine momy_d_xy(nx,ny,nz,dxi,dyi,visc,v,dvdt)
    implicit none
    integer , intent(in) :: nx,ny,nz
    real(rp), intent(in) :: dxi,dyi,visc
    real(rp), dimension(0:,0:,0:), intent(in   ) :: v
    real(rp), dimension( :, :, :), intent(inout) :: dvdt
    real(rp) :: dvdxp,dvdxm,dvdyp,dvdym
    integer :: i,j,k
    !
    !$acc parallel loop collapse(3) default(present) private(dvdxp,dvdxm,dvdyp,dvdym) async(1)
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(dvdxp,dvdxm,dvdyp,dvdym) &
    !$OMP SHARED(nx,ny,nz,dxi,dyi,v,dvdt,visc)
    do k=1,nz
      do j=1,ny
        do i=1,nx
          dvdxp = (v(i+1,j,k)-v(i,j,k))*dxi
          dvdxm = (v(i,j,k)-v(i-1,j,k))*dxi
          dvdyp = (v(i,j+1,k)-v(i,j,k))*dyi
          dvdym = (v(i,j,k)-v(i,j-1,k))*dyi
          dvdt(i,j,k) = dvdt(i,j,k) + &
                        (dvdxp-dvdxm)*visc*dxi + &
                        (dvdyp-dvdym)*visc*dyi
        end do
      end do
    end do
  end subroutine momy_d_xy
  !
  subroutine momz_d_xy(nx,ny,nz,dxi,dyi,visc,w,dwdt)
    implicit none
    integer , intent(in) :: nx,ny,nz
    real(rp), intent(in) :: dxi,dyi,visc
    real(rp), dimension(0:,0:,0:), intent(in   ) :: w
    real(rp), dimension( :, :, :), intent(inout) :: dwdt
    integer :: i,j,k
    real(rp) :: dwdxp,dwdxm,dwdyp,dwdym
    !
    !$acc parallel loop collapse(3) default(present) private(dwdxp,dwdxm,dwdyp,dwdym) async(1)
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(dwdxp,dwdxm,dwdyp,dwdym) &
    !$OMP SHARED(nx,ny,nz,dxi,dyi,w,dwdt,visc)
    do k=1,nz
      do j=1,ny
        do i=1,nx
          dwdxp = (w(i+1,j,k)-w(i,j,k))*dxi
          dwdxm = (w(i,j,k)-w(i-1,j,k))*dxi
          dwdyp = (w(i,j+1,k)-w(i,j,k))*dyi
          dwdym = (w(i,j,k)-w(i,j-1,k))*dyi
          dwdt(i,j,k) = dwdt(i,j,k) + &
                        (dwdxp-dwdxm)*visc*dxi + &
                        (dwdyp-dwdym)*visc*dyi
        end do
      end do
    end do
  end subroutine momz_d_xy
#endif
  !
  subroutine cmpt_wallshear(n,is_bound,l,dli,dzci,dzfi,visc,u,v,w,taux,tauy,tauz)
    use mod_param, only: cbcpre
    implicit none
    integer , intent(in ), dimension(3) :: n
    logical , intent(in ), dimension(0:1,3) :: is_bound
    real(rp), intent(in ), dimension(3)     :: l,dli
    real(rp), intent(in ), dimension(0:)    :: dzci,dzfi
    real(rp), intent(in )                   :: visc
    real(rp), intent(in ), dimension(0:,0:,0:) :: u,v,w
    real(rp), intent(out), dimension(3) :: taux,tauy,tauz
    real(rp) :: dudyp,dudym,dudzp,dudzm, &
                dvdxp,dvdxm,dvdzp,dvdzm, &
                dwdxp,dwdxm,dwdyp,dwdym
    !
    ! n.b.: replace scalars with reduction of tau(1:3,1:3) once the
    !       nvfortran bug for array reductions on Pascal architectures
    !       is solved; this subroutine is not used in production anyway
    !
    real(rp) :: tau21,tau31,tau12,tau32,tau13,tau23
    integer :: i,j,k,nx,ny,nz
    real(rp) :: dxi,dyi,lx,ly,lz
    real(rp) :: tau(3,3)
    !
    nx = n(1); ny = n(2); nz = n(3)
    dxi = dli(1); dyi = dli(2)
    lx = l(1); ly = l(2); lz = l(3)
    tau21 = 0._rp
    !$acc data copy(tau21) async(1)
    if(is_bound(0,2).and.cbcpre(0,2)//cbcpre(1,2) /= 'PP') then
      !$acc parallel loop collapse(2) default(present) private(dudyp) reduction(+:tau21) async(1)
      !$OMP PARALLEL DO DEFAULT(shared) private(dudyp) reduction(+:tau21)
      do k=1,nz
        do i=1,nx
          dudyp = (u(i,1 ,k)-u(i,0   ,k))*dyi*visc
          tau21 = tau21 + dudyp/(dxi*dzfi(k)*lx*lz)
        end do
      end do
    end if
    if(is_bound(1,2).and.cbcpre(0,2)//cbcpre(1,2) /= 'PP') then
      !$acc parallel loop collapse(2) default(present) private(dudym) reduction(+:tau21) async(1)
      !$OMP PARALLEL DO DEFAULT(shared) private(dudym) reduction(+:tau21)
      do k=1,nz
        do i=1,nx
          dudym = (u(i,ny,k)-u(i,ny+1,k))*dyi*visc
          tau21 = tau21 + dudym/(dxi*dzfi(k)*lx*lz)
        end do
      end do
    end if
    !$acc end data
    tau31 = 0._rp
    !$acc data copy(tau31) async(1)
    if(is_bound(0,3).and.cbcpre(0,3)//cbcpre(1,3) /= 'PP') then
      !$acc parallel loop collapse(2) default(present) private(dudzp) reduction(+:tau31) async(1)
      !$OMP PARALLEL DO DEFAULT(shared) private(dudzp) reduction(+:tau31)
      do j=1,ny
        do i=1,nx
          dudzp = (u(i,j,1 )-u(i,j,0   ))*dzci(0)*visc
          tau31 = tau31 + dudzp/(dxi*dyi*lx*ly)
        end do
      end do
    end if
    if(is_bound(1,3).and.cbcpre(0,3)//cbcpre(1,3) /= 'PP') then
      !$acc parallel loop collapse(2) default(present) private(dudzm) reduction(+:tau31) async(1)
      !$OMP PARALLEL DO DEFAULT(shared) private(dudzm) reduction(+:tau31)
      do j=1,ny
        do i=1,nx
          dudzm = (u(i,j,nz)-u(i,j,nz+1))*dzci(nz)*visc
          tau31 = tau31 + dudzm/(dxi*dyi*lx*ly)
        end do
      end do
    end if
    !$acc end data
    !
    tau12 = 0._rp
    !$acc data copy(tau12) async(1)
    if(is_bound(0,1).and.cbcpre(0,1)//cbcpre(1,1) /= 'PP') then
      !$acc parallel loop collapse(2) default(present) private(dvdxp) reduction(+:tau12) async(1)
      !$OMP PARALLEL DO DEFAULT(shared) private(dvdxp) reduction(+:tau12)
      do k=1,nz
        do j=1,ny
          dvdxp = (v(1  ,j,k)-v(0  ,j,k))*dxi*visc
          tau12 = tau12 + dvdxp/(dyi*dzfi(k)*ly*lz)
        end do
      end do
    end if
    if(is_bound(1,1).and.cbcpre(0,1)//cbcpre(1,1) /= 'PP') then
      !$acc parallel loop collapse(2) default(present) private(dvdxm) reduction(+:tau12) async(1)
      !$OMP PARALLEL DO DEFAULT(shared) private(dvdxm) reduction(+:tau12)
      do k=1,nz
        do j=1,ny
          dvdxm = (v(nx,j,k)-v(nx+1,j,k))*dxi*visc
          tau12 = tau12 + dvdxm/(dyi*dzfi(k)*ly*lz)
        end do
      end do
    end if
    !$acc end data
    tau32 = 0._rp
    !$acc data copy(tau32) async(1)
    if(is_bound(0,3).and.cbcpre(0,3)//cbcpre(1,3) /= 'PP') then
      !$acc parallel loop collapse(2) default(present) private(dvdzp) reduction(+:tau32) async(1)
      !$OMP PARALLEL DO DEFAULT(shared) private(dvdzp) reduction(+:tau32)
      do j=1,ny
        do i=1,nx
          dvdzp = (v(i,j,1 )-v(i,j,0   ))*dzci(0)*visc
          tau32 = tau32 + dvdzp/(dxi*dyi*lx*ly)
        end do
      end do
    end if
    if(is_bound(1,3).and.cbcpre(0,3)//cbcpre(1,3) /= 'PP') then
      !$acc parallel loop collapse(2) default(present) private(dvdzm) reduction(+:tau32) async(1)
      !$OMP PARALLEL DO DEFAULT(shared) private(dvdzm) reduction(+:tau32)
      do j=1,ny
        do i=1,nx
          dvdzm = (v(i,j,nz)-v(i,j,nz+1))*dzci(nz)*visc
          tau32 = tau32 + dvdzm/(dxi*dyi*lx*ly)
        end do
      end do
    end if
    !$acc end data
    !
    tau13 = 0._rp
    !$acc data copy(tau13) async(1)
    if(is_bound(0,1).and.cbcpre(0,1)//cbcpre(1,1) /= 'PP') then
      !$acc parallel loop collapse(2) default(present) private(dwdxp) reduction(+:tau13) async(1)
      !$OMP PARALLEL DO DEFAULT(shared) private(dwdxp) reduction(+:tau13)
      do k=1,nz
        do j=1,ny
          dwdxp = (w(1 ,j,k)-w(0   ,j,k))*dxi*visc
          tau13 = tau13 + dwdxp/(dyi*dzfi(k)*ly*lz)
        end do
      end do
    end if
    if(is_bound(1,1).and.cbcpre(0,1)//cbcpre(1,1) /= 'PP') then
      !$acc parallel loop collapse(2) default(present) private(dwdxm) reduction(+:tau13) async(1)
      !$OMP PARALLEL DO DEFAULT(shared) private(dwdxm) reduction(+:tau13)
      do k=1,nz
        do j=1,ny
          dwdxm = (w(nx,j,k)-w(nx+1,j,k))*dxi*visc
          tau13 = tau13 + dwdxm/(dyi*dzfi(k)*ly*lz)
        end do
      end do
    end if
    !$acc end data
    tau23 = 0._rp
    !$acc data copy(tau23) async(1)
    if(is_bound(0,2).and.cbcpre(0,2)//cbcpre(1,2) /= 'PP') then
      !$acc parallel loop collapse(2) default(present) private(dwdyp) reduction(+:tau23) async(1)
      !$OMP PARALLEL DO DEFAULT(shared) private(dwdyp) reduction(+:tau23)
      do k=1,nz
        do i=1,nx
          dwdyp = (w(i,1,k )-w(i,0   ,k))*dyi*visc
          tau23 = tau23 + dwdyp/(dxi*dzfi(k)*lx*lz)
        end do
      end do
    end if
    if(is_bound(1,2).and.cbcpre(0,2)//cbcpre(1,2) /= 'PP') then
      !$acc parallel loop collapse(2) default(present) private(dwdym) reduction(+:tau23) async(1)
      !$OMP PARALLEL DO DEFAULT(shared) private(dwdym) reduction(+:tau23)
      do k=1,nz
        do i=1,nx
          dwdym = (w(i,ny,k)-w(i,ny+1,k))*dyi*visc
          tau23 = tau23 + dwdym/(dxi*dzfi(k)*lx*lz)
        end do
      end do
    end if
    !$acc end data
    !$acc wait(1)
    tau(:,:) = 0._rp
    tau(2,1) = tau21
    tau(3,1) = tau31
    tau(1,2) = tau12
    tau(3,2) = tau32
    tau(1,3) = tau13
    tau(2,3) = tau23
    call MPI_ALLREDUCE(MPI_IN_PLACE,tau(1,1),9,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    taux(:) = tau(:,1)
    tauy(:) = tau(:,2)
    tauz(:) = tau(:,3)
  end subroutine cmpt_wallshear
  !
  subroutine mom_xyz_ad(nx,ny,nz,dxi,dyi,dzci,dzfi,visc,u,v,w,dudt,dvdt,dwdt,dudtd,dvdtd,dwdtd,psi_u,psi_v,psi_w)
    !
    ! lump all r.h.s. of momentum terms (excluding pressure) into a single fast kernel
    !
    integer , intent(in   ) :: nx,ny,nz
    real(rp), intent(in   ) :: dxi,dyi
    real(rp), intent(in   ), dimension(0:) :: dzci,dzfi
    real(rp), intent(in   ) :: visc
    real(rp), intent(in   ), dimension(0:,0:,0:) :: u,v,w
    real(rp), intent(in   ), dimension(0:,0:,0:), optional :: psi_u,psi_v,psi_w
    real(rp), intent(inout), dimension( :, :, :) :: dudt,dvdt,dwdt
    real(rp), intent(inout), dimension( :, :, :), optional :: dudtd,dvdtd,dwdtd
    integer  :: i,j,k
    real(rp) :: u_ccm,u_pcm,u_cpm,u_cmc,u_pmc,u_mcc,u_ccc,u_pcc,u_mpc,u_cpc,u_cmp,u_mcp,u_ccp, &
                v_ccm,v_pcm,v_cpm,v_cmc,v_pmc,v_mcc,v_ccc,v_pcc,v_mpc,v_cpc,v_cmp,v_mcp,v_ccp, &
                w_ccm,w_pcm,w_cpm,w_cmc,w_pmc,w_mcc,w_ccc,w_pcc,w_mpc,w_cpc,w_cmp,w_mcp,w_ccp
    real(rp) :: uuip,uuim,uvjp,uvjm,uwkp,uwkm, &
                uvip,uvim,vvjp,vvjm,wvkp,wvkm, &
                uwip,uwim,vwjp,vwjm,wwkp,wwkm
    real(rp) :: dudxp,dudxm,dudyp,dudym,dudzp,dudzm, &
                dvdxp,dvdxm,dvdyp,dvdym,dvdzp,dvdzm, &
                dwdxp,dwdxm,dwdyp,dwdym,dwdzp,dwdzm
    real(rp) :: dudt_s ,dvdt_s ,dwdt_s , &
                dudtd_s,dvdtd_s,dwdtd_s, &
                dudtd_xy_s,dudtd_z_s   , &
                dvdtd_xy_s,dvdtd_z_s   , &
                dwdtd_xy_s,dwdtd_z_s
    !$acc parallel loop collapse(3) default(present) async(1) &
    !$acc private(u_ccm,u_pcm,u_cpm,u_cmc,u_pmc,u_mcc,u_ccc,u_pcc,u_mpc,u_cpc,u_cmp,u_mcp,u_ccp) &
    !$acc private(v_ccm,v_pcm,v_cpm,v_cmc,v_pmc,v_mcc,v_ccc,v_pcc,v_mpc,v_cpc,v_cmp,v_mcp,v_ccp) &
    !$acc private(w_ccm,w_pcm,w_cpm,w_cmc,w_pmc,w_mcc,w_ccc,w_pcc,w_mpc,w_cpc,w_cmp,w_mcp,w_ccp) &
    !$acc private(uuip,uuim,uvjp,uvjm,uwkp,uwkm) &
    !$acc private(uvip,uvim,vvjp,vvjm,wvkp,wvkm) &
    !$acc private(uwip,uwim,vwjp,vwjm,wwkp,wwkm) &
    !$acc private(dudxp,dudxm,dudyp,dudym,dudzp,dudzm) &
    !$acc private(dvdxp,dvdxm,dvdyp,dvdym,dvdzp,dvdzm) &
    !$acc private(dwdxp,dwdxm,dwdyp,dwdym,dwdzp,dwdzm) &
    !$acc private(dudt_s ,dvdt_s ,dwdt_s ) &
    !$acc private(dudtd_s,dvdtd_s,dwdtd_s) &
    !$acc private(dudtd_xy_s,dudtd_z_s) &
    !$acc private(dvdtd_xy_s,dvdtd_z_s) &
    !$acc private(dwdtd_xy_s,dwdtd_z_s)
    !$OMP PARALLEL DO DEFAULT(shared) &
    !$OMP PRIVATE(u_ccm,u_pcm,u_cpm,u_cmc,u_pmc,u_mcc,u_ccc,u_pcc,u_mpc,u_cpc,u_cmp,u_mcp,u_ccp) &
    !$OMP PRIVATE(v_ccm,v_pcm,v_cpm,v_cmc,v_pmc,v_mcc,v_ccc,v_pcc,v_mpc,v_cpc,v_cmp,v_mcp,v_ccp) &
    !$OMP PRIVATE(w_ccm,w_pcm,w_cpm,w_cmc,w_pmc,w_mcc,w_ccc,w_pcc,w_mpc,w_cpc,w_cmp,w_mcp,w_ccp) &
    !$OMP PRIVATE(uuip,uuim,uvjp,uvjm,uwkp,uwkm) &
    !$OMP PRIVATE(uvip,uvim,vvjp,vvjm,wvkp,wvkm) &
    !$OMP PRIVATE(uwip,uwim,vwjp,vwjm,wwkp,wwkm) &
    !$OMP PRIVATE(dudxp,dudxm,dudyp,dudym,dudzp,dudzm) &
    !$OMP PRIVATE(dvdxp,dvdxm,dvdyp,dvdym,dvdzp,dvdzm) &
    !$OMP PRIVATE(dwdxp,dwdxm,dwdyp,dwdym,dwdzp,dwdzm) &
    !$OMP PRIVATE(dudt_s ,dvdt_s ,dwdt_s ) &
    !$OMP PRIVATE(dudtd_s,dvdtd_s,dwdtd_s) &
    !$OMP PRIVATE(dudtd_xy_s,dudtd_z_s) &
    !$OMP PRIVATE(dvdtd_xy_s,dvdtd_z_s) &
    !$OMP PRIVATE(dwdtd_xy_s,dwdtd_z_s)
    do k=1,nz
      do j=1,ny
        do i=1,nx
          !
          ! touch u,v,w sequentially
          !
          u_ccm = u(i  ,j  ,k-1)
          u_pcm = u(i+1,j  ,k-1)
          u_cpm = u(i  ,j+1,k-1)
          u_cmc = u(i  ,j-1,k  )
          u_pmc = u(i+1,j-1,k  )
          u_mcc = u(i-1,j  ,k  )
          u_ccc = u(i  ,j  ,k  )
          u_pcc = u(i+1,j  ,k  )
          u_mpc = u(i-1,j+1,k  )
          u_cpc = u(i  ,j+1,k  )
          u_cmp = u(i  ,j-1,k+1)
          u_mcp = u(i-1,j  ,k+1)
          u_ccp = u(i  ,j  ,k+1)
          !
          v_ccm = v(i  ,j  ,k-1)
          v_pcm = v(i+1,j  ,k-1)
          v_cpm = v(i  ,j+1,k-1)
          v_cmc = v(i  ,j-1,k  )
          v_pmc = v(i+1,j-1,k  )
          v_mcc = v(i-1,j  ,k  )
          v_ccc = v(i  ,j  ,k  )
          v_pcc = v(i+1,j  ,k  )
          v_mpc = v(i-1,j+1,k  )
          v_cpc = v(i  ,j+1,k  )
          v_cmp = v(i  ,j-1,k+1)
          v_mcp = v(i-1,j  ,k+1)
          v_ccp = v(i  ,j  ,k+1)
          !
          w_ccm = w(i  ,j  ,k-1)
          w_pcm = w(i+1,j  ,k-1)
          w_cpm = w(i  ,j+1,k-1)
          w_cmc = w(i  ,j-1,k  )
          w_pmc = w(i+1,j-1,k  )
          w_mcc = w(i-1,j  ,k  )
          w_ccc = w(i  ,j  ,k  )
          w_pcc = w(i+1,j  ,k  )
          w_mpc = w(i-1,j+1,k  )
          w_cpc = w(i  ,j+1,k  )
          w_cmp = w(i  ,j-1,k+1)
          w_mcp = w(i-1,j  ,k+1)
          w_ccp = w(i  ,j  ,k+1)
          !
          ! x diffusion
          !
          dudxp = (u_pcc-u_ccc)*dxi
          dudxm = (u_ccc-u_mcc)*dxi
          dudyp = (u_cpc-u_ccc)*dyi
          dudym = (u_ccc-u_cmc)*dyi
          dudzp = (u_ccp-u_ccc)*dzci(k  )
          dudzm = (u_ccc-u_ccm)*dzci(k-1)
#if defined(_IBM) && defined(_SIMPLE)
          if (psi_u(i,j,k) == 0._rp) then ! if not in solid
            if (psi_u(i,j+1,k) == 1._rp) dudyp = -2._rp*u(i,j,k)*dyi
            if (psi_u(i,j-1,k) == 1._rp) dudym =  2._rp*u(i,j,k)*dyi
            if (psi_u(i,j,k+1) == 1._rp) dudzp = -2._rp*u(i,j,k)*dzci(k)
            if (psi_u(i,j,k-1) == 1._rp) dudzm =  2._rp*u(i,j,k)*dzci(k)
          endif
          ! dudyp = dudyp - ((1. - psi_u(i,j,k))*psi_u(i,j+1,k)*dudyp) &
                        ! + ((1. - psi_u(i,j,k))*psi_u(i,j+1,k)*(-2._rp*u(i,j,k)*dyi))
          ! dudym = dudym - ((1. - psi_u(i,j,k))*psi_u(i,j-1,k)*dudym) &
                        ! + ((1. - psi_u(i,j,k))*psi_u(i,j-1,k)*( 2._rp*u(i,j,k)*dyi))
          ! dudzp = dudzp - ((1. - psi_u(i,j,k))*psi_u(i,j,k+1)*dudzp) &
                        ! + ((1. - psi_u(i,j,k))*psi_u(i,j,k+1)*(-2._rp*u(i,j,k)*dzci(k)))
          ! dudzm = dudzm - ((1. - psi_u(i,j,k))*psi_u(i,j,k-1)*dudzm) &
                        ! + ((1. - psi_u(i,j,k))*psi_u(i,j,k-1)*( 2._rp*u(i,j,k)*dzci(k)))
#endif
          !
          ! x advection
          !
          uuip  = 0.25_rp*(u_pcc+u_ccc)*(u_pcc+u_ccc)
          uuim  = 0.25_rp*(u_mcc+u_ccc)*(u_mcc+u_ccc)
          uvjp  = 0.25_rp*(u_cpc+u_ccc)*(v_pcc+v_ccc)
          uvjm  = 0.25_rp*(u_cmc+u_ccc)*(v_pmc+v_cmc)
          uwkp  = 0.25_rp*(u_ccp+u_ccc)*(w_pcc+w_ccc)
          uwkm  = 0.25_rp*(u_ccm+u_ccc)*(w_pcm+w_ccm)
#if defined(_IBM) && defined(_SIMPLE)
          if (psi_u(i,j,k) == 0._rp) then ! if not in solid
            if (psi_u(i,j+1,k) == 1._rp) uvjp  = 0._rp
            if (psi_u(i,j-1,k) == 1._rp) uvjm  = 0._rp
            if (psi_u(i,j,k+1) == 1._rp) uwkp  = 0._rp
            if (psi_u(i,j,k-1) == 1._rp) uwkm  = 0._rp
          endif
          ! uvjp = uvjp - ((1. - psi_u(i,j,k))*psi_u(i,j+1,k)*uvjp)
          ! uvjm = uvjm - ((1. - psi_u(i,j,k))*psi_u(i,j-1,k)*uvjm)
          ! uwkp = uwkp - ((1. - psi_u(i,j,k))*psi_u(i,j,k+1)*uwkp)
          ! uwkm = uwkm - ((1. - psi_u(i,j,k))*psi_u(i,j,k-1)*uwkm)
#endif
          !
          dudtd_xy_s = visc*(dudxp-dudxm)*dxi + &
                       visc*(dudyp-dudym)*dyi
          dudtd_z_s  = visc*(dudzp-dudzm)*dzfi(k)
          dudt_s     = -(uuip -uuim )*dxi - &
                        (uvjp -uvjm )*dyi - &
                        (uwkp -uwkm )*dzfi(k)
          !
          ! y diffusion
          !
          dvdxp = (v_pcc-v_ccc)*dxi
          dvdxm = (v_ccc-v_mcc)*dxi
          dvdyp = (v_cpc-v_ccc)*dyi
          dvdym = (v_ccc-v_cmc)*dyi
          dvdzp = (v_ccp-v_ccc)*dzci(k  )
          dvdzm = (v_ccc-v_ccm)*dzci(k-1)
#if defined(_IBM) && defined(_SIMPLE)
          if (psi_v(i,j,k) == 0._rp) then ! if not in solid
            if (psi_v(i+1,j,k) == 1._rp) dvdxp = -2._rp*v(i,j,k)*dxi
            if (psi_v(i-1,j,k) == 1._rp) dvdxm =  2._rp*v(i,j,k)*dxi
            if (psi_v(i,j,k+1) == 1._rp) dvdzp = -2._rp*v(i,j,k)*dzci(k)
            if (psi_v(i,j,k-1) == 1._rp) dvdzm =  2._rp*v(i,j,k)*dzci(k)
          endif
          ! dvdxp = dvdxp - ((1. - psi_v(i,j,k))*psi_v(i+1,j,k)*dvdxp) &
                        ! + ((1. - psi_v(i,j,k))*psi_v(i+1,j,k)*(-2._rp*v(i,j,k)*dxi))
          ! dvdxm = dvdxm - ((1. - psi_v(i,j,k))*psi_v(i-1,j,k)*dvdxm) &
                        ! + ((1. - psi_v(i,j,k))*psi_v(i-1,j,k)*( 2._rp*v(i,j,k)*dxi))
          ! dvdzp = dvdzp - ((1. - psi_v(i,j,k))*psi_v(i,j,k+1)*dvdzp) &
                        ! + ((1. - psi_v(i,j,k))*psi_v(i,j,k+1)*(-2._rp*v(i,j,k)*dzci(k)))
          ! dvdzm = dvdzm - ((1. - psi_v(i,j,k))*psi_v(i,j,k-1)*dvdzm) &
                        ! + ((1. - psi_v(i,j,k))*psi_v(i,j,k-1)*( 2._rp*v(i,j,k)*dzci(k)))
#endif
          !
          ! y advection
          !
          uvip  = 0.25_rp*(u_ccc+u_cpc)*(v_ccc+v_pcc)
          uvim  = 0.25_rp*(u_mcc+u_mpc)*(v_ccc+v_mcc)
          vvjp  = 0.25_rp*(v_ccc+v_cpc)*(v_ccc+v_cpc)
          vvjm  = 0.25_rp*(v_ccc+v_cmc)*(v_ccc+v_cmc)
          wvkp  = 0.25_rp*(w_ccc+w_cpc)*(v_ccp+v_ccc)
          wvkm  = 0.25_rp*(w_ccm+w_cpm)*(v_ccm+v_ccc)
#if defined(_IBM) && defined(_SIMPLE)
          if (psi_v(i,j,k) == 0._rp) then ! if not in solid
            if (psi_v(i+1,j,k) == 1._rp) uvip = 0._rp
            if (psi_v(i-1,j,k) == 1._rp) uvim = 0._rp
            if (psi_v(i,j,k+1) == 1._rp) wvkp = 0._rp
            if (psi_v(i,j,k-1) == 1._rp) wvkm = 0._rp
          endif
          ! uvip = uvip - ((1. - psi_v(i,j,k))*psi_v(i+1,j,k)*uvip)
          ! uvim = uvim - ((1. - psi_v(i,j,k))*psi_v(i-1,j,k)*uvim)
          ! wvkp = wvkp - ((1. - psi_v(i,j,k))*psi_v(i,j,k+1)*wvkp)
          ! wvkm = wvkm - ((1. - psi_v(i,j,k))*psi_v(i,j,k-1)*wvkm)
#endif
          !
          dvdtd_xy_s = visc*(dvdxp-dvdxm)*dxi + &
                       visc*(dvdyp-dvdym)*dyi
          dvdtd_z_s  = visc*(dvdzp-dvdzm)*dzfi(k)
          dvdt_s     = -(uvip -uvim )*dxi - &
                        (vvjp -vvjm )*dyi - &
                        (wvkp -wvkm )*dzfi(k)
          !
          ! z diffusion
          !
          dwdxp = (w_pcc-w_ccc)*dxi
          dwdxm = (w_ccc-w_mcc)*dxi
          dwdyp = (w_cpc-w_ccc)*dyi
          dwdym = (w_ccc-w_cmc)*dyi
          dwdzp = (w_ccp-w_ccc)*dzfi(k+1)
          dwdzm = (w_ccc-w_ccm)*dzfi(k  )
#if defined(_IBM) && defined(_SIMPLE)
          if (psi_w(i,j,k) == 0._rp) then ! if not in solid
           if (psi_w(i+1,j,k) == 1._rp) dwdxp = -2._rp*w(i,j,k)*dxi
           if (psi_w(i-1,j,k) == 1._rp) dwdxm =  2._rp*w(i,j,k)*dxi
           if (psi_w(i,j+1,k) == 1._rp) dwdyp = -2._rp*w(i,j,k)*dyi
           if (psi_w(i,j-1,k) == 1._rp) dwdym =  2._rp*w(i,j,k)*dyi
          endif
          ! dwdxp = dwdxp - ((1. - psi_v(i,j,k))*psi_v(i+1,j,k)*dwdxp) &
                        ! + ((1. - psi_v(i,j,k))*psi_v(i+1,j,k)*(-2._rp*w(i,j,k)*dxi))
          ! dwdxm = dwdxm - ((1. - psi_v(i,j,k))*psi_v(i-1,j,k)*dwdxm) &
                        ! + ((1. - psi_v(i,j,k))*psi_v(i-1,j,k)*( 2._rp*w(i,j,k)*dxi))
          ! dwdyp = dwdyp - ((1. - psi_v(i,j,k))*psi_v(i,j,k+1)*dwdyp) &
                        ! + ((1. - psi_v(i,j,k))*psi_v(i,j,k+1)*(-2._rp*w(i,j,k)*dyi))
          ! dwdym = dwdym - ((1. - psi_v(i,j,k))*psi_v(i,j,k-1)*dwdym) &
                        ! + ((1. - psi_v(i,j,k))*psi_v(i,j,k-1)*( 2._rp*w(i,j,k)*dyi))
#endif
          !
          ! z advection
          !
          uwip  = 0.25_rp*(w_ccc+w_pcc)*(u_ccc+u_ccp)
          uwim  = 0.25_rp*(w_ccc+w_mcc)*(u_mcc+u_mcp)
          vwjp  = 0.25_rp*(w_ccc+w_cpc)*(v_ccc+v_ccp)
          vwjm  = 0.25_rp*(w_ccc+w_cmc)*(v_cmc+v_cmp)
          wwkp  = 0.25_rp*(w_ccc+w_ccp)*(w_ccc+w_ccp)
          wwkm  = 0.25_rp*(w_ccc+w_ccm)*(w_ccc+w_ccm)
#if defined(_IBM) && defined(_SIMPLE)
          if (psi_w(i,j,k) == 0._rp) then ! if not in solid
           if (psi_w(i+1,j,k) == 1._rp) uwip = 0._rp
           if (psi_w(i-1,j,k) == 1._rp) uwim = 0._rp
           if (psi_w(i,j+1,k) == 1._rp) vwjp = 0._rp
           if (psi_w(i,j-1,k) == 1._rp) vwjm = 0._rp
          endif
          ! uwip = uvip - ((1. - psi_w(i,j,k))*psi_w(i+1,j,k)*uwip)
          ! uwim = uwim - ((1. - psi_w(i,j,k))*psi_w(i-1,j,k)*uwim)
          ! vwjp = vwjp - ((1. - psi_w(i,j,k))*psi_w(i,j+1,k)*vwjp)
          ! vwjm = vwjm - ((1. - psi_w(i,j,k))*psi_w(i,j-1,k)*vwjm)
#endif
          !
          dwdtd_xy_s = visc*(dwdxp-dwdxm)*dxi + &
                       visc*(dwdyp-dwdym)*dyi
          dwdtd_z_s  = visc*(dwdzp-dwdzm)*dzci(k)
          dwdt_s     =  -(uwip -uwim )*dxi - &
                         (vwjp -vwjm )*dyi - &
                         (wwkp -wwkm )*dzci(k)
#if defined(_IMPDIFF)
#if defined(_IMPDIFF_1D)
          dudt_s = dudt_s + dudtd_xy_s
          dvdt_s = dvdt_s + dvdtd_xy_s
          dwdt_s = dwdt_s + dwdtd_xy_s
          dudtd_s = dudtd_z_s
          dvdtd_s = dvdtd_z_s
          dwdtd_s = dwdtd_z_s
#else
          dudtd_s = dudtd_xy_s + dudtd_z_s
          dvdtd_s = dvdtd_xy_s + dvdtd_z_s
          dwdtd_s = dwdtd_xy_s + dwdtd_z_s
#endif
          dudt( i,j,k) = dudt_s
          dvdt( i,j,k) = dvdt_s
          dwdt( i,j,k) = dwdt_s
          dudtd(i,j,k) = dudtd_s
          dvdtd(i,j,k) = dvdtd_s
          dwdtd(i,j,k) = dwdtd_s
#else
          dudt_s = dudt_s + dudtd_xy_s + dudtd_z_s
          dvdt_s = dvdt_s + dvdtd_xy_s + dvdtd_z_s
          dwdt_s = dwdt_s + dwdtd_xy_s + dwdtd_z_s
          dudt(i,j,k) = dudt_s
          dvdt(i,j,k) = dvdt_s
          dwdt(i,j,k) = dwdt_s
#endif
        end do
      end do
    end do
  end subroutine mom_xyz_ad
end module mod_mom
