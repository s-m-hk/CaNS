! -
!
! SPDX-FileCopyrightText: Copyright (c) 2017-2022 Pedro Costa and the CaNS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
module mod_scal
  use mod_weno , only: weno5
  use mod_types
  implicit none
  private
  public scal,cmpt_scalflux,scal_forcing
  contains
#if !defined(_WENO)
  subroutine scal(nx,ny,nz,dxi,dyi,dzi,dzci,dzfi,alph_f, &
#if defined(_IBM)
                  alph_s,alph,psi, &
#endif
                  u,v,w,s, &
#if defined(_IMPDIFF)
                  dsdtd, &
#if defined(_CONSTANT_COEFFS_DIFF)
                  dsdtdc, &
#endif
#endif
                  dsdt)
    !
    implicit none
    integer , intent(in) :: nx,ny,nz
    real(rp), intent(in) :: dxi,dyi,dzi,alph_f
    real(rp), intent(in), dimension(0:) :: dzci,dzfi
#if defined(_IBM)
    real(rp), intent(in) :: alph_s
    real(rp), intent(in), dimension(0:,0:,0:) :: alph
    real(rp), intent(in), dimension(0:,0:,0:) :: psi
#endif
    real(rp), dimension(0:,0:,0:), intent(in) :: u,v,w,s
#if defined(_IMPDIFF)
    real(rp), dimension(:,:,:), intent(out) :: dsdtd
#if defined(_CONSTANT_COEFFS_DIFF)
    real(rp), dimension(:,:,:), intent(out) :: dsdtdc
#endif
#endif
    real(rp), dimension(:,:,:), intent(out) :: dsdt
    integer :: i,j,k
    real(rp) :: usip,usim,vsjp,vsjm,wskp,wskm
    real(rp) :: dsdxp,dsdxm,dsdyp,dsdym,dsdzp,dsdzm
    real(rp) :: source,maxalph 
    !
#if defined(_IBM)
    maxalph = max(alph_f,alph_s)
#endif
    !
    !$acc parallel loop collapse(3) default(present) private(usip,usim,vsjp,vsjm,wskp,wskm,dsdxp,dsdxm,dsdyp,dsdym,dsdzp,dsdzm) async(1)
    !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared)  PRIVATE(usip,usim,vsjp,vsjm,wskp,wskm,dsdxp,dsdxm,dsdyp,dsdym,dsdzp,dsdzm)
    do k=1,nz
      do j=1,ny
        do i=1,nx
#if defined(_HEAT_SOURCE)
          source  = u(i,j,k)
#else
          source  = 0.0_rp
#endif
          usim  = 0.5*( s(i-1,j,k)+s(i,j,k) )*u(i-1,j,k)
          usip  = 0.5*( s(i+1,j,k)+s(i,j,k) )*u(i  ,j,k)
          vsjm  = 0.5*( s(i,j-1,k)+s(i,j,k) )*v(i,j-1,k)
          vsjp  = 0.5*( s(i,j+1,k)+s(i,j,k) )*v(i,j  ,k)
          wskm  = 0.5*( s(i,j,k-1)+s(i,j,k) )*w(i,j,k-1)
          wskp  = 0.5*( s(i,j,k+1)+s(i,j,k) )*w(i,j,k  )
#if defined(_IBM)
          dsdxp = (0.5_rp*(alph(i+1,j,k)+alph(i,j,k)))*(s(i+1,j,k)-s(i  ,j,k))*dxi
          dsdxm = (0.5_rp*(alph(i-1,j,k)+alph(i,j,k)))*(s(i  ,j,k)-s(i-1,j,k))*dxi
          dsdyp = (0.5_rp*(alph(i,j+1,k)+alph(i,j,k)))*(s(i,j+1,k)-s(i,j  ,k))*dyi
          dsdym = (0.5_rp*(alph(i,j-1,k)+alph(i,j,k)))*(s(i,j  ,k)-s(i,j-1,k))*dyi
          dsdzp = (0.5_rp*(alph(i,j,k+1)+alph(i,j,k)))*(s(i,j,k+1)-s(i,j,k  ))*dzci(k  )
          dsdzm = (0.5_rp*(alph(i,j,k-1)+alph(i,j,k)))*(s(i,j,k  )-s(i,j,k-1))*dzci(k-1)
          !
#if defined(_IMPDIFF)
          dsdtd(i,j,k) = (dsdxp-dsdxm)*dxi + &
                         (dsdyp-dsdym)*dyi + &
                         (dsdzp-dsdzm)*dzfi(k)
#if defined(_CONSTANT_COEFFS_DIFF)           
          dsdxp = maxalph*(s(i+1,j,k)-s(i  ,j,k))*dxi
          dsdxm = maxalph*(s(i  ,j,k)-s(i-1,j,k))*dxi
          dsdyp = maxalph*(s(i,j+1,k)-s(i,j  ,k))*dyi
          dsdym = maxalph*(s(i,j  ,k)-s(i,j-1,k))*dyi
          dsdzp = maxalph*(s(i,j,k+1)-s(i,j,k  ))*dzci(k  )
          dsdzm = maxalph*(s(i,j,k  )-s(i,j,k-1))*dzci(k-1)
          !
          dsdtdc(i,j,k) = (dsdxp-dsdxm)*dxi + &
                          (dsdyp-dsdym)*dyi + &
                          (dsdzp-dsdzm)*dzfi(k)
#endif
          !
          dsdt(i,j,k) = dxi*(     -usip + usim ) + &
                        dyi*(     -vsjp + vsjm ) + &
                        dzfi(k)*( -wskp + wskm ) + &
                        source
#else
          dsdt(i,j,k) = dxi*(     -usip + usim ) + (dsdxp-dsdxm)*dxi + &
                        dyi*(     -vsjp + vsjm ) + (dsdyp-dsdym)*dyi + &
                        dzfi(k)*( -wskp + wskm ) + (dsdzp-dsdzm)*dzfi(k) + &
                        source
#endif
#else
          dsdxp = (s(i+1,j,k)-s(i  ,j,k))*dxi
          dsdxm = (s(i  ,j,k)-s(i-1,j,k))*dxi
          dsdyp = (s(i,j+1,k)-s(i,j  ,k))*dyi
          dsdym = (s(i,j  ,k)-s(i,j-1,k))*dyi
          dsdzp = (s(i,j,k+1)-s(i,j,k  ))*dzci(k  )
          dsdzm = (s(i,j,k  )-s(i,j,k-1))*dzci(k-1)
          !
#if defined(_IMPDIFF)
          dsdtd(i,j,k) = (dsdxp-dsdxm)*alph_f*dxi + &
                         (dsdyp-dsdym)*alph_f*dyi + &
                         (dsdzp-dsdzm)*alph_f*dzfi(k)
          !
          dsdt(i,j,k) = dxi*(     -usip + usim ) + &
                        dyi*(     -vsjp + vsjm ) + &
                        dzfi(k)*( -wskp + wskm ) + &
                        source
#else
          dsdt(i,j,k) = dxi*(     -usip + usim ) + (dsdxp-dsdxm)*alph_f*dxi + &
                        dyi*(     -vsjp + vsjm ) + (dsdyp-dsdym)*alph_f*dyi + &
                        dzfi(k)*( -wskp + wskm ) + (dsdzp-dsdzm)*alph_f*dzfi(k) + &
                        source
#endif
#endif
        end do
      end do
    end do
  end subroutine scal
#else
  subroutine scal(nx,ny,nz,dxi,dyi,dzi,dzci,dzfi,alph_f, &
#if defined(_IBM)
                  alph,psi, &
#endif
                  u,v,w,s, &
#if defined(_IMPDIFF)
                  dsdtd, &
#endif
                  dsdt)
    !
    implicit none
    integer , intent(in) :: nx,ny,nz
    real(rp), intent(in) :: dxi,dyi,dzi,alph_f
    real(rp), intent(in), dimension(0:) :: dzci,dzfi
    real(rp), dimension(0:,0:,0:), intent(in) :: u,v,w
    real(rp), dimension(-2:,-2:,-2:), intent(in) :: s
#if defined(_IBM)
    real(rp), intent(in), dimension(0:,0:,0:) :: alph
    real(rp), intent(in), dimension(0:,0:,0:) :: psi
#endif
#if defined(_IMPDIFF)
    real(rp), dimension(:,:,:), intent(out) :: dsdtd
#endif
    real(rp), dimension(:,:,:), intent(out) :: dsdt
    integer :: i,j,k
    real(rp) :: usip,usim,vsjp,vsjm,wskp,wskm
    real(rp) :: dsdip,dsdim,dsdjp,dsdjm,dsdkp,dsdkm,dsdtd_xy,dsdtd_z,dsdta,dsdt_s
    integer  :: psis,psisip,psisim,psisjp,psisjm,psiskp,psiskm
    real(rp) :: alphip,alphim,alphjp,alphjm,alphkp,alphkm
#if defined(_WENO)
    call weno5(nx,ny,nz,nh_s,dxi,dyi,dzi,s,u,v,w,dsdt)
#endif
    !$acc parallel loop collapse(3) default(present) private(psis,psisip,psisim,psisjp,psisjm,psiskp,psiskm,alphip,alphim,alphjp,alphjm,alphkp,alphkm,usip,usim,vsjp,vsjm,wskp,wskm,dsdip,dsdim,dsdjp,dsdjm,dsdkp,dsdkm,dsdtd_xy,dsdtd_z,dsdta,dsdt_s) async(1)
    !$OMP PARALLEL DO COLLAPSE(3) DEFAULT(shared) private(psis,psisip,psisim,psisjp,psisjm,psiskp,psiskm,alphip,alphim,alphjp,alphjm,alphkp,alphkm,usip,usim,vsjp,vsjm,wskp,wskm,dsdip,dsdim,dsdjp,dsdjm,dsdkp,dsdkm,dsdtd_xy,dsdtd_z,dsdta,dsdt_s)
    do k=1,nz
      do j=1,ny
        do i=1,nx
#if !defined(_WENO) && defined(_IBM) && defined(_SIMPLE)
         alphip = 0.5_rp*(alph(i+1,j,k)+alph(i,j,k))
         alphim = 0.5_rp*(alph(i-1,j,k)+alph(i,j,k))
         alphjp = 0.5_rp*(alph(i,j+1,k)+alph(i,j,k))
         alphjm = 0.5_rp*(alph(i,j-1,k)+alph(i,j,k))
         alphkp = 0.5_rp*(alph(i,j,k+1)+alph(i,j,k))
         alphkm = 0.5_rp*(alph(i,j,k-1)+alph(i,j,k))
         !
#if defined(_ISOTHERMAL)
         psis   = psi(i,j,k)  
         psisip = psi(i+1,j,k)
         psisim = psi(i-1,j,k)
         psisjp = psi(i,j+1,k)
         psisjm = psi(i,j-1,k)
         psiskp = psi(i,j,k+1)
         psiskm = psi(i,j,k-1)
#endif
#endif
         !
         dsdip = (s(i+1,j,k)-s(i  ,j,k))*dxi
         dsdim = (s(i  ,j,k)-s(i-1,j,k))*dxi
         dsdjp = (s(i,j+1,k)-s(i,j  ,k))*dyi
         dsdjm = (s(i,j  ,k)-s(i,j-1,k))*dyi
         dsdkp = (s(i,j,k+1)-s(i,j,k  ))*dzci(k  )
         dsdkm = (s(i,j,k  )-s(i,j,k-1))*dzci(k-1)
         !
#if !defined(_WENO) && defined(_IBM) && defined(_SIMPLE) && defined(_ISOTHERMAL)
         if (psis == 0.0_rp) then ! if not in solid
          if (psisip == 1.0_rp) s(i+1,j,k) = 2.0_rp*solidtemp; s(i,j,k) = 2.0_rp*s(i,j,k)
          if (psisim == 1.0_rp) s(i-1,j,k) = 2.0_rp*solidtemp; s(i,j,k) = 2.0_rp*s(i,j,k)
          if (psisjp == 1.0_rp) s(i,j+1,k) = 2.0_rp*solidtemp; s(i,j,k) = 2.0_rp*s(i,j,k)
          if (psisjm == 1.0_rp) s(i,j-1,k) = 2.0_rp*solidtemp; s(i,j,k) = 2.0_rp*s(i,j,k)
          if (psiskp == 1.0_rp) s(i,j,k+1) = 2.0_rp*solidtemp; s(i,j,k) = 2.0_rp*s(i,j,k)
          if (psiskm == 1.0_rp) s(i,j,k-1) = 2.0_rp*solidtemp; s(i,j,k) = 2.0_rp*s(i,j,k)
         endif
#endif
         !
#if !defined(_WENO)
         usip  = 0.5_rp*( s(i+1,j,k)+s(i,j,k) )*u(i  ,j,k)
         usim  = 0.5_rp*( s(i-1,j,k)+s(i,j,k) )*u(i-1,j,k)
         vsjp  = 0.5_rp*( s(i,j+1,k)+s(i,j,k) )*v(i,j  ,k)
         vsjm  = 0.5_rp*( s(i,j-1,k)+s(i,j,k) )*v(i,j-1,k)
         wskp  = 0.5_rp*( s(i,j,k+1)+s(i,j,k) )*w(i,j,k  )
         wskm  = 0.5_rp*( s(i,j,k-1)+s(i,j,k) )*w(i,j,k-1)
#endif
         !
#if !defined(_WENO) && defined(_IBM) && defined(_SIMPLE) && defined(_ISOTHERMAL)
         if (psis == 0.0_rp) then ! if not in solid
          if (psisip == 1.0_rp) usip = solidtemp*u(i ,j,k)
          if (psisim == 1.0_rp) usim = solidtemp*u(i-1,j,k)
          if (psisjp == 1.0_rp) vsjp = solidtemp*v(i,j ,k)
          if (psisjm == 1.0_rp) vsjm = solidtemp*v(i,j-1,k)
          if (psiskp == 1.0_rp) wskp = solidtemp*w(i,j,k )
          if (psiskm == 1.0_rp) wskm = solidtemp*w(i,j,k-1)
         endif
#endif
#if !defined(_WENO) & defined(_IBM) && defined(_SIMPLE)
         dsdtd_xy = (alphip*dsdip-alphim*dsdim)*dxi + &
                    (alphjp*dsdjp-alphjm*dsdjm)*dyi
         dsdtd_z  = (alphkp*dsdkp-alphkm*dsdkm)*dzfi(k)
         dsdta    = - (usip - usim)*dxi &    
                    - (vsjp - vsjm)*dyi &    
                    - (wskp - wskm)*dzfi(k)
#elif defined(_WENO) & !defined(_IBM)
         dsdtd_xy = alph_f*( &
                             (dsdip-dsdim)*dxi + &
                             (dsdjp-dsdjm)*dyi   &
                           )
         dsdtd_z  = alph_f*(dsdkp-dsdkm)*dzfi(k)
         dsdta    = dsdt(i,j,k)
#else
         dsdtd_xy = alph_f*( &
                            (dsdip-dsdim)*dxi + &
                            (dsdjp-dsdjm)*dyi   &
                           )
         dsdtd_z  = alph_f*(dsdkp-dsdkm)*dzfi(k)
         dsdta    = - (usip - usim)*dxi &    
                    - (vsjp - vsjm)*dyi &    
                    - (wskp - wskm)*dzfi(k)
#endif
#if defined(_IMPDIFF)
         dsdtd(i,j,k) = dsdtd_xy + dsdtd_z
         dsdt_s = dsdta
#else
         dsdt_s = dsdta + dsdtd_xy + dsdtd_z
#endif
         dsdt(i,j,k) = dsdt_s
        end do
      end do
    end do
 end subroutine scal
#endif
 subroutine cmpt_scalflux(n,is_bound,l,dli,dzci,dzfi,alpha,s,flux)
  use mpi
  use mod_param, only: cbctmp
  implicit none
  integer , intent(in ), dimension(3) :: n
  logical , intent(in ), dimension(0:1,3) :: is_bound
  real(rp), intent(in ), dimension(3)     :: l,dli
  real(rp), intent(in ), dimension(0:)    :: dzci,dzfi
  real(rp), intent(in )                   :: alpha
  real(rp), intent(in ), dimension(0:,0:,0:) :: s
  real(rp), intent(out), dimension(3) :: flux
  real(rp) :: dsdxp,dsdxm,dsdyp,dsdym,dsdzp,dsdzm
  real(rp) :: flux_x,flux_y,flux_z
  integer :: ierr
  !
  integer :: i,j,k,nx,ny,nz
  real(rp) :: dxi,dyi,lx,ly,lz
  !
  nx = n(1); ny = n(2); nz = n(3)
  dxi = dli(1); dyi = dli(2)
  lx = l(1); ly = l(2); lz = l(3)
  flux_x = 0.0_rp
  !$acc data copy(flux_x) async(1)
  if(is_bound(0,1).and.cbctmp(0,1)//cbctmp(1,1) /= 'PP') then
    !$acc parallel loop collapse(2) default(present) private(dsdxp) reduction(+:flux_x) async(1)
    !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(shared) private(dsdxp) reduction(+:flux_x)
    do k=1,nz
      do j=1,ny
        dsdxp = (s(1 ,j,k)-s(0   ,j,k))*dxi*alpha
        flux_x = flux_x + dsdxp/(dyi*dzfi(k)*ly*lz)
      end do
    end do
  end if
  if(is_bound(1,1).and.cbctmp(0,1)//cbctmp(1,1) /= 'PP') then
    !$acc parallel loop collapse(2) default(present) private(dsdxm) reduction(+:flux_x) async(1)
    !$OMP PARALLEL DO   collapse(2) default(shared)  private(dsdxm) reduction(+:flux_x)
    do k=1,nz
      do j=1,ny
        dsdxm = (s(nx,j,k)-s(nx+1,j,k))*dxi*alpha
        flux_x = flux_x + dsdxm/(dyi*dzfi(k)*ly*lz)
      end do
    end do
  end if
  !$acc end data
  flux_y = 0.0_rp
  !$acc data copy(flux_y) async(1)
  if(is_bound(0,2).and.cbctmp(0,2)//cbctmp(1,2) /= 'PP') then
    !$acc parallel loop collapse(2) default(present) private(dsdyp) reduction(+:flux_y) async(1)
    !$OMP PARALLEL DO   collapse(2) default(shared)  private(dsdyp) reduction(+:flux_y)
    do k=1,nz
      do i=1,nx
        dsdyp = (s(i,1 ,k)-s(i,0   ,k))*dyi*alpha
        flux_y = flux_y + dsdyp/(dxi*dzfi(k)*lx*lz)
      end do
    end do
  end if
  if(is_bound(1,2).and.cbctmp(0,2)//cbctmp(1,2) /= 'PP') then
    !$acc parallel loop collapse(2) default(present) private(dsdym) reduction(+:flux_y) async(1)
    !$OMP PARALLEL DO   collapse(2) default(shared)  private(dsdym) reduction(+:flux_y)
    do k=1,nz
      do i=1,nx
        dsdym = (s(i,ny,k)-s(i,ny+1,k))*dyi*alpha
        flux_y = flux_y + dsdym/(dxi*dzfi(k)*lx*lz)
      end do
    end do
  end if
  !$acc end data
  flux_z = 0.0_rp
  !$acc data copy(flux_z) async(1)
  if(is_bound(0,3).and.cbctmp(0,3)//cbctmp(1,3) /= 'PP') then
    !$acc parallel loop collapse(2) default(present) private(dsdzp) reduction(+:flux_z) async(1)
    !$OMP PARALLEL DO   collapse(2) default(shared)  private(dsdzp) reduction(+:flux_z)
    do j=1,ny
      do i=1,nx
        dsdzp = (s(i,j,1 )-s(i,j,0   ))*dzci(0)*alpha
        flux_z = flux_z + dsdzp/(dxi*dyi*lx*ly)
      end do
    end do
  end if
  if(is_bound(1,3).and.cbctmp(0,3)//cbctmp(1,3) /= 'PP') then
    !$acc parallel loop collapse(2) default(present) private(dsdzm) reduction(+:flux_z) async(1)
    !$OMP PARALLEL DO   collapse(2) default(shared)  private(dsdzm) reduction(+:flux_z)
    do j=1,ny
      do i=1,nx
        dsdzm = (s(i,j,nz)-s(i,j,nz+1))*dzci(nz)*alpha
        flux_z = flux_z + dsdzm/(dxi*dyi*lx*ly)
      end do
    end do
  end if
  !$acc end data
  !$acc wait(1)
  flux(:) = [flux_x,flux_y,flux_z]
  call MPI_ALLREDUCE(MPI_IN_PLACE,flux(3),3,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
 end subroutine cmpt_scalflux
!
  subroutine scal_forcing(n,is_forced,f,s)
    integer , intent(in   ), dimension(3) :: n
    logical , intent(in   ) :: is_forced
    real(rp), intent(in   ) :: f
    real(rp), intent(inout), dimension(0:,0:,0:) :: s
    if(is_forced) then
      !$acc kernels default(present) async(1)
      s(1:n(1),1:n(2),1:n(3)) = s(1:n(1),1:n(2),1:n(3)) + f
      !$acc end kernels
    end if
  end subroutine scal_forcing
!
end module mod_scal
