! -
!
! SPDX-FileCopyrightText: Copyright (c) 2017-2022 Pedro Costa and the CaNS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
#define _FAST_MOM_KERNELS
module mod_rk
  use mod_mom    ,  only: momx_a,momy_a,momz_a, &
                          momx_d,momy_d,momz_d, &
                          momx_p,momy_p,momz_p, &
                          cmpt_wallshear,bulk_forcing
#if defined(_IMPDIFF_1D)
  use mod_mom    ,  only: momx_d_xy,momy_d_xy,momz_d_xy, &
                          momx_d_z ,momy_d_z ,momz_d_z
#endif
#if defined(_FAST_MOM_KERNELS)
  use mod_mom    ,  only: mom_xyz_ad
#endif
  use mod_scal   ,  only: scal,cmpt_scalflux,scal_forcing
  use mod_source ,  only: grav_src
#if defined(_IBM)
  use mod_forcing,  only: ib_force,bulk_vel, &
#if defined(_HEAT_TRANSFER)
                          force_scal, &
#endif
                          force_vel
#endif
  use mod_common_mpi, only: myid
  use mod_utils,      only: bulk_mean,swap
  use mod_types
  implicit none
  private
  public rk,rk_scal
  contains
  subroutine rk(time_scheme,rkpar,n,nh_v,nh_s,dli,l,zc,zf,dzci,dzfi,grid_vol_ratio_c,grid_vol_ratio_f,visc,dt,dto,p, &
                is_bound,is_forced,velf,bforce,tauxo,tauyo,tauzo,u,v,w, &
                dl,dzc,dzf, &
#if defined(_IBM)
                psi_s,psi_u,psi_v,psi_w, &
                fx,fy,fz, &
                fibm, &
#endif
#if defined(_HEAT_TRANSFER)
                s, &
#endif
#if defined(_LPP)
                duconv,dvconv,dwconv, &
#endif
                f)
    !
    ! low-storage 3rd-order Runge-Kutta scheme
    ! for time integration of the momentum equations.
    !
    implicit none
    character(len=3), intent(in   )              :: time_scheme
    real(rp), intent(in   ), dimension(2)        :: rkpar
    integer , intent(in   ), dimension(3)        :: n
    integer , intent(in   )                      :: nh_v,nh_s
    real(rp), intent(in   )                      :: visc,dt,dto
    real(rp), intent(in   ), dimension(3)        :: dli,l
    real(rp), intent(in   ), dimension(3)        :: dl
    real(rp), intent(in   ), dimension(0:)       :: dzc,dzf
    real(rp), intent(in   ), dimension(0:)       :: zc,zf,dzci,dzfi
    real(rp), intent(in   ), dimension(0:)       :: grid_vol_ratio_c,grid_vol_ratio_f
    real(rp), intent(in   ), dimension(0:,0:,0:) :: p
    logical , intent(in   ), dimension(0:1,3)    :: is_bound
    logical , intent(in   ), dimension(4)        :: is_forced
    real(rp), intent(in   ), dimension(3)        :: velf,bforce
    real(rp), intent(inout), dimension(3)        :: tauxo,tauyo,tauzo
    real(rp), intent(inout), dimension(0:,0:,0:) :: u,v,w
    real(rp), intent(inout), dimension(4) :: f
    real(rp), target     , allocatable, dimension(:,:,:), save :: dudtrk_t ,dvdtrk_t ,dwdtrk_t , &
                                                                  dudtrko_t,dvdtrko_t,dwdtrko_t
    real(rp), pointer    , contiguous , dimension(:,:,:), save :: dudtrk   ,dvdtrk   ,dwdtrk   , &
                                                                  dudtrko  ,dvdtrko  ,dwdtrko
    real(rp),              allocatable, dimension(:,:,:), save :: dudtrkd  ,dvdtrkd  ,dwdtrkd
#if defined(_IBM)
    real(rp), intent(in   ), dimension(0:,0:,0:) :: psi_s,psi_u,psi_v,psi_w
    real(rp), intent(inout), dimension(0:,0:,0:) :: fx,fy,fz
    real(rp), intent(out), dimension(4) :: fibm
#endif
#if defined(_HEAT_TRANSFER)
    real(rp), intent(in   ), dimension(1-nh_s:,1-nh_s:,1-nh_s:) :: s
#endif
#if defined(_LPP)
    real(rp), intent(out), dimension(:,:,:) :: duconv,dvconv,dwconv
#endif
    logical, save :: is_first = .true.
    real(rp) :: factor1,factor2,factor12
    real(rp), dimension(3) :: taux,tauy,tauz
    integer :: i,j,k
    real(rp) :: mean
    !
    if(    time_scheme.eq.'ab2') then !1st order Euler used for first time-step
     if(is_first) then
      factor1 = 1.0_rp*dt
      factor2 = 0.0_rp*dt
     else
      factor1 = (1.0_rp+0.5_rp*(dt/dto) )*dt
      factor2 = (      -0.5_rp*(dt/dto) )*dt
     endif
    elseif(time_scheme.eq.'rk3') then
      factor1 = rkpar(1)*dt
      factor2 = rkpar(2)*dt
    endif
    factor12 = factor1 + factor2
    !
    ! initialization
    !
    if(is_first) then ! leverage save attribute to allocate these arrays on the device only once
      is_first = .false.
      allocate(dudtrk_t( n(1),n(2),n(3)),dvdtrk_t( n(1),n(2),n(3)),dwdtrk_t( n(1),n(2),n(3)))
      allocate(dudtrko_t(n(1),n(2),n(3)),dvdtrko_t(n(1),n(2),n(3)),dwdtrko_t(n(1),n(2),n(3)))
      !$acc enter data create(dudtrk_t ,dvdtrk_t ,dwdtrk_t ) async(1)
      !$acc enter data create(dudtrko_t,dvdtrko_t,dwdtrko_t) async(1)
      !$acc kernels default(present) async(1) ! not really necessary
      dudtrko_t(:,:,:) = 0.0_rp
      dvdtrko_t(:,:,:) = 0.0_rp
      dwdtrko_t(:,:,:) = 0.0_rp
      !$acc end kernels
#if defined(_IMPDIFF)
      allocate(dudtrkd(n(1),n(2),n(3)),dvdtrkd(n(1),n(2),n(3)),dwdtrkd(n(1),n(2),n(3)))
      !$acc enter data create(dudtrkd,dvdtrkd,dwdtrkd) async(1)
#endif
      dudtrk  => dudtrk_t
      dvdtrk  => dvdtrk_t
      dwdtrk  => dwdtrk_t
      dudtrko => dudtrko_t
      dvdtrko => dvdtrko_t
      dwdtrko => dwdtrko_t
    end if
    !
#if defined(_FAST_MOM_KERNELS)
    call mom_xyz_ad(n(1),n(2),n(3),dli(1),dli(2),dzci,dzfi,visc,u,v,w,dudtrk,dvdtrk,dwdtrk,dudtrkd,dvdtrkd,dwdtrkd &
#if defined(_IBM) && defined(_SIMPLE)
                    ,psi_u,psi_v,psi_w &
#endif
                    )
#else
    !$acc kernels default(present) async(1)
    !$OMP PARALLEL WORKSHARE
    dudtrk(:,:,:) = 0.0_rp
    dvdtrk(:,:,:) = 0.0_rp
    dwdtrk(:,:,:) = 0.0_rp
    !$OMP END PARALLEL WORKSHARE
    !$acc end kernels
#if !defined(_IMPDIFF)
    call momx_d(n(1),n(2),n(3),dli(1),dli(2),dzci,dzfi,visc,u, &
#if defined(_IBM)
                    psi_u,psi_v,psi_w, &
#endif
                    dudtrk)
    call momy_d(n(1),n(2),n(3),dli(1),dli(2),dzci,dzfi,visc,v, &
#if defined(_IBM)
                    psi_u,psi_v,psi_w, &
#endif
                    dvdtrk)
    call momz_d(n(1),n(2),n(3),dli(1),dli(2),dzci,dzfi,visc,w, &
#if defined(_IBM)
                    psi_u,psi_v,psi_w, &
#endif
                    dwdtrk)
#else
    !$acc kernels default(present) async(1)
    !$OMP PARALLEL WORKSHARE
    dudtrkd(:,:,:) = 0.0_rp
    dvdtrkd(:,:,:) = 0.0_rp
    dwdtrkd(:,:,:) = 0.0_rp
    !$OMP END PARALLEL WORKSHARE
    !$acc end kernels
#if !defined(_IMPDIFF_1D)
    call momx_d(n(1),n(2),n(3),dli(1),dli(2),dzci,dzfi,visc,u,
     &
#if defined(_IBM)
                    psi_u,psi_v,psi_w, &
#endif
                    dudtrkd)
    call momy_d(n(1),n(2),n(3),dli(1),dli(2),dzci,dzfi,visc,v, &
#if defined(_IBM)
                    psi_u,psi_v,psi_w, &
#endif
                    dvdtrkd)
    call momz_d(n(1),n(2),n(3),dli(1),dli(2),dzci,dzfi,visc,w, &
#if defined(_IBM)
                    psi_u,psi_v,psi_w, &
#endif
                    dwdtrkd)
#else
    call momx_d_xy(n(1),n(2),n(3),dli(1),dli(2),visc,u,dudtrk )
    call momy_d_xy(n(1),n(2),n(3),dli(1),dli(2),visc,v,dvdtrk )
    call momz_d_xy(n(1),n(2),n(3),dli(1),dli(2),visc,w,dwdtrk )
    call momx_d_z( n(1),n(2),n(3),dzci  ,dzfi  ,visc,u,dudtrkd)
    call momy_d_z( n(1),n(2),n(3),dzci  ,dzfi  ,visc,v,dvdtrkd)
    call momz_d_z( n(1),n(2),n(3),dzci  ,dzfi  ,visc,w,dwdtrkd)
#endif
#endif
    call momx_a(n(1),n(2),n(3),dli(1),dli(2),dzfi,u,v,w, &
#if defined(_IBM)
                psi_u,psi_v,psi_w, &
#endif
                dudtrk)
    call momy_a(n(1),n(2),n(3),dli(1),dli(2),dzfi,u,v,w, &
#if defined(_IBM)
                psi_u,psi_v,psi_w, &
#endif
                dvdtrk)
    call momz_a(n(1),n(2),n(3),dli(1),dli(2),dzci,u,v,w, &
#if defined(_IBM)
                psi_u,psi_v,psi_w, &
#endif
                dwdtrk)
#endif
    !
    !$acc parallel loop collapse(3) default(present) async(1)
    !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
#if !defined(_FAST_MOM_KERNELS)
          u(i,j,k) = u(i,j,k) + factor1*dudtrk(i,j,k) + factor2*dudtrko(i,j,k)
          v(i,j,k) = v(i,j,k) + factor1*dvdtrk(i,j,k) + factor2*dvdtrko(i,j,k)
          w(i,j,k) = w(i,j,k) + factor1*dwdtrk(i,j,k) + factor2*dwdtrko(i,j,k)
#else
          u(i,j,k) = u(i,j,k) + factor1*dudtrk(i,j,k) + factor2*dudtrko(i,j,k) + &
                                factor12*(bforce(1) - dli(1)*( p(i+1,j,k)-p(i,j,k)))
          !
          v(i,j,k) = v(i,j,k) + factor1*dvdtrk(i,j,k) + factor2*dvdtrko(i,j,k) + &
                                factor12*(bforce(2) - dli(2)*( p(i,j+1,k)-p(i,j,k)))
          !
          w(i,j,k) = w(i,j,k) + factor1*dwdtrk(i,j,k) + factor2*dwdtrko(i,j,k) + &
                                factor12*(bforce(3) - dzci(k)*(p(i,j,k+1)-p(i,j,k)))
#endif
#if defined(_IMPDIFF)
          u(i,j,k) = u(i,j,k) + factor12*dudtrkd(i,j,k)
          v(i,j,k) = v(i,j,k) + factor12*dvdtrkd(i,j,k)
          w(i,j,k) = w(i,j,k) + factor12*dwdtrkd(i,j,k)
#endif
        end do
      end do
    end do
    !
    ! swap d?dtrk <-> d?dtrko
    !
    call swap(dudtrk,dudtrko)
    call swap(dvdtrk,dvdtrko)
    call swap(dwdtrk,dwdtrko)
!#if 0 /*pressure gradient term treated explicitly later */
!    !$acc kernels
!    !$OMP PARALLEL WORKSHARE
!    dudtrk(:,:,:) = 0._rp
!    dvdtrk(:,:,:) = 0._rp
!    dwdtrk(:,:,:) = 0._rp
!    !$OMP END PARALLEL WORKSHARE
!    !$acc end kernels
!    call momx_p(n(1),n(2),n(3),dli(1),bforce(1),p,dudtrk)
!    call momy_p(n(1),n(2),n(3),dli(2),bforce(2),p,dvdtrk)
!    call momz_p(n(1),n(2),n(3),dzci  ,bforce(3),p,dwdtrk)
!    !$acc parallel loop collapse(3)
!    !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared)
!    do k=1,n(3)
!      do j=1,n(2)
!        do i=1,n(1)
!          u(i,j,k) = u(i,j,k) + factor12*dudtrk(i,j,k)
!          v(i,j,k) = v(i,j,k) + factor12*dvdtrk(i,j,k)
!          w(i,j,k) = w(i,j,k) + factor12*dwdtrk(i,j,k)
!        end do
!      end do
!    end do
!#endif
#if !defined(_FAST_MOM_KERNELS)
    !$acc parallel loop collapse(3) default(present) async(1)
    !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          u(i,j,k) = u(i,j,k) + factor12*(bforce(1) - dli(1)*( p(i+1,j,k)-p(i,j,k)))
          v(i,j,k) = v(i,j,k) + factor12*(bforce(2) - dli(2)*( p(i,j+1,k)-p(i,j,k)))
          w(i,j,k) = w(i,j,k) + factor12*(bforce(3) - dzci(k)*(p(i,j,k+1)-p(i,j,k)))
        end do
      end do
    end do
#endif
    !
    call grav_src(n(1),n(2),n(3), &
                  factor12, &
#if defined(_HEAT_TRANSFER)
                  nh_s,s, &
#endif
                  u,v,w)
    !
    ! compute bulk velocity forcing
    !
    call cmpt_bulk_forcing(n,l,dl,dzf,dzc,zf,zc,is_forced,velf,grid_vol_ratio_c,grid_vol_ratio_f,&
#if defined(_IBM)
                           psi_u,psi_v,psi_w,&
#endif
                           u,v,w,f)
    call bulk_forcing(n,is_forced,f,u,v,w)
#if defined(_IBM)
    call ib_force(n,dl,dzc,dzf,l,psi_u,psi_v,psi_w,u,v,w,fx,fy,fz,fibm)
#endif
#if defined(_IMPDIFF)
    !
    ! compute rhs of helmholtz equation
    !
    !$acc parallel loop collapse(3) default(present) async(1)
    !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          u(i,j,k) = u(i,j,k) - 0.5_rp*factor12*dudtrkd(i,j,k)
          v(i,j,k) = v(i,j,k) - 0.5_rp*factor12*dvdtrkd(i,j,k)
          w(i,j,k) = w(i,j,k) - 0.5_rp*factor12*dwdtrkd(i,j,k)
        end do
      end do
    end do
#endif
  end subroutine rk
  subroutine rk_scal(time_scheme,rkpar,n,nh_v,nh_s,dli,zc,zf,dzci,dzfi,grid_vol_ratio_f,dt,dto,l,u,v,w,alph_f, &
                     is_bound,is_forced,is_cmpt_wallflux,tmpf,ssource, &
                     dl,dzc,dzf, &
#if defined(_IBM)
                     alph_s,al, &
                     psi_s,psi_u, &
#if defined(_VOLUME) && defined(_HEAT_TRANSFER) && defined(_ISOTHERMAL)
                     fibm, &
#endif
#endif
                     s,f,fluxo,wall_flux)
    !
    ! low-storage 3rd-order Runge-Kutta scheme
    ! for time integration of the scalar field.
    !
    implicit none
    character(len=3) , intent(in   )                            :: time_scheme
    real(rp), intent(in   ), dimension(2)                       :: rkpar
    integer , intent(in   ), dimension(3)                       :: n
    integer , intent(in   )                                     :: nh_v,nh_s
    real(rp), intent(in   ), dimension(3)                       :: dli,l
    real(rp), intent(in   ), dimension(0:)                      :: zc,zf,dzci,dzfi
    real(rp), intent(in   ), dimension(0:)                      :: grid_vol_ratio_f
    logical , intent(in   ), dimension(0:1,3)                   :: is_bound
    logical , intent(in   ), dimension(4)                       :: is_forced
    logical , intent(in   )                                     :: is_cmpt_wallflux
    real(rp), intent(in   )                                     :: tmpf,ssource
    real(rp), intent(in   ), dimension(3)                       :: dl
    real(rp), intent(in   ), dimension(0:)                      :: dzc,dzf
    real(rp), intent(in   )                                     :: alph_f,dt,dto
    real(rp), intent(in   ), dimension(0:,0:,0:)                :: u,v,w
    real(rp), intent(inout), dimension(1-nh_s:,1-nh_s:,1-nh_s:) :: s
    real(rp), intent(inout), dimension(4)                       :: f
    real(rp), intent(inout), dimension(3)                       :: fluxo
    real(rp), intent(out)                                       :: wall_flux
    real(rp), target     , allocatable, dimension(:,:,:), save  :: dsdtrk_t, dsdtrko_t
    real(rp), pointer    , contiguous , dimension(:,:,:), save  :: dsdtrk, dsdtrko
#if defined(_IMPDIFF)
    real(rp),              allocatable, dimension(:,:,:), save  :: dsdtrkd
#if defined(_CONSTANT_COEFFS_DIFF)
    real(rp),              allocatable, dimension(:,:,:), save  :: dsdtrkdc
#endif
#endif
    real(rp)                                                    :: factor1,factor2,factor12
#if defined(_IBM)
    real(rp), intent(in   )                                     :: alph_s
    real(rp), intent(in   ), dimension(0:,0:,0:)                :: al
    real(rp), intent(in   ), dimension(0:,0:,0:)                :: psi_s,psi_u
#if defined(_VOLUME) && defined(_HEAT_TRANSFER) && defined(_ISOTHERMAL)
    real(rp), intent(out  ), dimension(4)                       :: fibm
#endif
#endif
    real(rp), dimension(3)                                      :: flux
    real(rp)                                                    :: mean,vol_f
    integer                                                     :: i,j,k
    logical, save                                               :: is_first = .true.
    !
    if(    time_scheme.eq.'ab2') then
     if(is_first) then !1st order Euler used for first time-step
      factor1 = 1.0_rp*dt
      factor2 = 0.0_rp*dt
     else
      factor1 = (1.0_rp+0.5_rp*(dt/dto) )*dt
      factor2 = (      -0.5_rp*(dt/dto) )*dt
     endif
    elseif(time_scheme.eq.'rk3') then
      factor1 = rkpar(1)*dt
      factor2 = rkpar(2)*dt
    endif
    factor12 = factor1 + factor2
    !
    ! initialization
    !
    if(is_first) then ! leverage save attribute to allocate these arrays on the device only once
      is_first = .false.
      allocate(dsdtrk_t(n(1),n(2),n(3)),dsdtrko_t(n(1),n(2),n(3)))
      !$acc enter data create(dsdtrk_t, dsdtrko_t) async(1)
      !$acc kernels default(present) async(1)
      dsdtrko_t(:,:,:) = 0.0_rp
      !$acc end kernels
      dsdtrk  => dsdtrk_t
      dsdtrko => dsdtrko_t
#if defined(_IMPDIFF)
      allocate(dsdtrkd(n(1),n(2),n(3)))
      !$acc enter data create(dsdtrkd) async(1)
      !$acc kernels default(present) async(1)
      !$OMP PARALLEL WORKSHARE
      dsdtrkd(:,:,:) = 0.0_rp
      !$OMP END PARALLEL WORKSHARE
      !$acc end kernels
#if defined(_CONSTANT_COEFFS_DIFF)
      allocate(dsdtrkdc(n(1),n(2),n(3)))
      !$acc enter data create(dsdtrkdc) async(1)
      !$acc kernels default(present) async(1)
      !$OMP PARALLEL WORKSHARE
      dsdtrkdc(:,:,:) = 0.0_rp
      !$OMP END PARALLEL WORKSHARE
      !$acc end kernels
#endif
#endif
    end if
    !
    if(is_cmpt_wallflux) then
      call cmpt_scalflux(n,is_bound,l,dli,dzci,dzfi,alph_f,s,flux)
      wall_flux = (factor1*sum(flux(:)/l(:)) + factor2*sum(fluxo(:)/l(:)))
      fluxo(:) = flux(:)
    end if
    !
    call scal(n(1),n(2),n(3),dli(1),dli(2),dli(3),dzci,dzfi,alph_f, &
#if defined(_IBM)
              alph_s,al,psi_s, &
#endif
              u,v,w,s, &
#if defined(_IMPDIFF)
              dsdtrkd, &
#if defined(_CONSTANT_COEFFS_DIFF)
              dsdtrkdc, &
#endif
#endif
              dsdtrk)
    !$acc parallel loop collapse(3) default(present) async(1)
    !$OMP PARALLEL DO COLLAPSE(3) DEFAULT(shared)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          s(i,j,k) = s(i,j,k) + factor1*dsdtrk(i,j,k) + factor2*dsdtrko(i,j,k) + factor12*ssource
#if defined(_IMPDIFF)
          s(i,j,k) = s(i,j,k) + factor12*dsdtrkd(i,j,k)
#endif
        end do
      end do
    end do
    !
    ! swap dsdtrk <-> dsdtrko
    !
    call swap(dsdtrk,dsdtrko)
    !
    ! Maintain a constant bulk temperature
    !
    if(is_forced(4)) then
#if !defined(_IBM)
     call bulk_mean(n,nh_s,grid_vol_ratio_f,s,mean)
     f(4) = tmpf - mean
#else
     call bulk_vel(n,nh_s,dl,dzf,zc,l,psi_s,s,mean,vol_f)
     call force_vel(n,1,dl,dzf,l,psi_s,s,tmpf,vol_f,mean,f(4))
#endif
     call scal_forcing(n,is_forced(4),f(4),s)
    endif
#if defined(_IBM) && defined(_VOLUME) && defined(_HEAT_TRANSFER) && defined(_ISOTHERMAL)
    !
    ! Volume penalization to impose isothermal conditions
    !
    call force_scal(n,nh_s,dl,dzf,l,psi_s,s,fibm)
#endif
#if defined(_IMPDIFF)
    !
    ! compute rhs of helmholtz equation
    !
    !$acc parallel loop collapse(3) default(present) async(1)
    !$OMP PARALLEL DO COLLAPSE(3) DEFAULT(shared)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
#if defined(_CONSTANT_COEFFS_DIFF)
          s(i,j,k) = s(i,j,k) - 0.5_rp*factor12*dsdtrkdc(i,j,k)
#else
          s(i,j,k) = s(i,j,k) - 0.5_rp*factor12*dsdtrkd(i,j,k)
#endif
        end do
      end do
    end do
#endif
  end subroutine rk_scal
  !
  subroutine cmpt_bulk_forcing(n,l,dl,dzf,dzc,zf,zc,is_forced,velf,grid_vol_ratio_c,grid_vol_ratio_f,&
#if defined(_IBM)
                               psi_u,psi_v,psi_w,&
#endif
                               u,v,w,f)
    use mod_param,      only: force_fluid_only
    implicit none
    integer , intent(in   ), dimension(3)  :: n
    real(rp), intent(in   ), dimension(3)  :: l,dl
    real(rp), intent(in   ), dimension(0:) :: dzc,dzf,zf,zc
    logical , intent(in   ), dimension(3)  :: is_forced
    real(rp), intent(in   ), dimension(3)  :: velf
    real(rp), intent(in   ), dimension(0:) :: grid_vol_ratio_c,grid_vol_ratio_f
#if defined(_IBM)
    real(rp), intent(in), dimension(0:,0:,0:) :: psi_u,psi_v,psi_w
#endif
    real(rp), intent(inout), dimension(0:,0:,0:) :: u,v,w
    real(rp), intent(out  ), dimension(3) :: f
    real(rp) :: vol_f,mean
    !
    ! bulk velocity forcing
    !
    f(1:3) = 0.0_rp
#if !defined(_IBM)
    if(is_forced(1)) then
     call bulk_mean(n,1,grid_vol_ratio_c,u,mean)
     f(1) = velf(1) - mean
    end if
    if(is_forced(2)) then
     call bulk_mean(n,1,grid_vol_ratio_c,v,mean)
     f(2) = velf(2) - mean
    end if
    if(is_forced(3)) then
     call bulk_mean(n,1,grid_vol_ratio_c,w,mean)
     f(3) = velf(3) - mean
    end if
#else
    if(is_forced(1)) then
     call bulk_vel(n,1,dl,dzf,zc,l,psi_u,u,mean,vol_f)
     if(force_fluid_only) then
      call force_vel(n,1,dl,dzf,l,psi_u,u,velf(1),vol_f,mean,f(1))
     else
      f(1) = velf(1) - mean
     endif
    endif
    if(is_forced(2)) then
     call bulk_vel(n,1,dl,dzf,zc,l,psi_v,v,mean,vol_f)
     if(force_fluid_only) then
      call force_vel(n,1,dl,dzf,l,psi_v,v,velf(2),vol_f,mean,f(2))
     else
      f(2) = velf(2) - mean
     endif
    endif
    if(is_forced(3)) then
     call bulk_vel(n,1,dl,dzc,zf,l,psi_w,w,mean,vol_f)
     if(force_fluid_only) then
      call force_vel(n,1,dl,dzc,l,psi_w,w,velf(3),vol_f,mean,f(3))
     else
      f(3) = velf(3) - mean
     endif
    endif
#endif
  end subroutine cmpt_bulk_forcing
  !
  subroutine cmpt_bulk_forcing_alternative(rkpar,n,dli,l,dzci,dzfi,visc,dt,is_bound,is_forced,u,v,w,tauxo,tauyo,tauzo,f,is_first)
    !
    ! computes the pressure gradient to be added to the flow that perfectly balances the wall shear stresses
    ! this effectively prescribes zero net acceleration, which allows to sustain a constant mass flux
    !
    implicit none
    real(rp), intent(in), dimension(2) :: rkpar
    integer , intent(in), dimension(3) :: n
    real(rp), intent(in) :: visc,dt
    real(rp), intent(in   ), dimension(3) :: dli,l
    real(rp), intent(in   ), dimension(0:) :: dzci,dzfi
    logical , intent(in   ), dimension(0:1,3)    :: is_bound
    logical , intent(in   ), dimension(3) :: is_forced
    real(rp), intent(in   ), dimension(0:,0:,0:) :: u,v,w
    real(rp), intent(inout), dimension(3) :: tauxo,tauyo,tauzo
    real(rp), intent(inout), dimension(3) :: f
    real(rp), dimension(3) :: f_aux
    logical , intent(in   ) :: is_first
    real(rp), dimension(3) :: taux,tauy,tauz
    real(rp) :: factor1,factor2,factor12
    !
    factor1 = rkpar(1)*dt
    factor2 = rkpar(2)*dt
    factor12 = (factor1 + factor2)/2.0_rp
    !
    taux(:) = 0._rp
    tauy(:) = 0._rp
    tauz(:) = 0._rp
    call cmpt_wallshear(n,is_forced,is_bound,l,dli,dzci,dzfi,visc,u,v,w,taux,tauy,tauz)
#if !defined(_IMPDIFF)
    if(is_first) then
      f(1) = (factor1*sum(taux(:)/l(:)) + factor2*sum(tauxo(:)/l(:)))
      f(2) = (factor1*sum(tauy(:)/l(:)) + factor2*sum(tauyo(:)/l(:)))
      f(3) = (factor1*sum(tauz(:)/l(:)) + factor2*sum(tauzo(:)/l(:)))
      tauxo(:) = taux(:)
      tauyo(:) = tauy(:)
      tauzo(:) = tauz(:)
    end if
#else
#if defined(_IMPDIFF_1D)
    f_aux(1) = factor12*taux(3)/l(3)
    f_aux(2) = factor12*tauy(3)/l(3)
    if(is_first) then
      f(1) = factor1*taux(2)/l(2) + factor2*tauxo(2)/l(2) + f_aux(1)
      f(2) = factor1*tauy(1)/l(1) + factor2*tauyo(1)/l(1) + f_aux(2)
      f(3) = factor1*sum(tauz(1:2)/l(1:2)) + factor2*sum(tauzo(1:2)/l(1:2))
      tauxo(1:2) = taux(1:2)
      tauyo(1:2) = tauy(1:2)
      tauzo(1:2) = tauz(1:2)
    else
      f(1) = f(1) + f_aux(1)
      f(2) = f(2) + f_aux(2)
    end if
#else
    f_aux(:) = factor12*[sum(taux(:)/l(:)), &
                         sum(tauy(:)/l(:)), &
                         sum(tauz(:)/l(:))]
    if(is_first) then
       f(:) = f_aux(:)
    else
       f(:) = f(:) + f_aux(:)
    endif
#endif
#endif
  end subroutine cmpt_bulk_forcing_alternative
end module mod_rk
