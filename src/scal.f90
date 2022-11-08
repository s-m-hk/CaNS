! -
!
! SPDX-FileCopyrightText: Copyright (c) 2017-2022 Pedro Costa and the CaNS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
module mod_scal
  use mod_gradls ,only: weno5
  use mod_types, only: rp
  implicit none
  private
  public scal
  contains
  subroutine scal(nx,ny,nz,dxi,dyi,dzi,nh_s,dzci,dzfi,alph_f, &
#if defined(_IBM)
                  alph_s,alph,psi, &
#endif
                  u,v,w,s,dsdt)
    !
    implicit none
    integer , intent(in) :: nx,ny,nz,nh_s
    real(rp), intent(in) :: dxi,dyi,dzi,alph_f
    real(rp), intent(in), dimension(0:) :: dzci,dzfi
    real(rp), dimension(0:,0:,0:), intent(in) :: u,v,w
    real(rp), dimension(1-nh_s:,1-nh_s:,1-nh_s:), intent(in) :: s
#if defined(_IBM)
    real(rp), intent(in), dimension(0:,0:,0:) :: alph
    real(rp), intent(in), dimension(0:,0:,0:) :: psi
    real(rp), intent(in) :: alph_s
#endif
    real(rp), dimension(:,:,:), intent(out) :: dsdt
    integer :: i,j,k
    real(rp) :: usip,usim,vsjp,vsjm,wskp,wskm
    real(rp) :: dsdip,dsdim,dsdjp,dsdjm,dsdkp,dsdkm,dsdtd_xy,dsdtd_z,dsdta,dsdt_s
    real(rp) :: psis,psisip,psisim,psisjp,psisjm,psiskp,psiskm
    real(rp) :: alphip,alphim,alphjp,alphjm,alphkp,alphkm
#if defined(_WENO)
    call weno5(nx,ny,nz,nh_s,dxi,dyi,dzi,s,u,v,w,dsdt)
#endif
    !$acc parallel loop collapse(3) default(present) private(psis,psisip,psisim,psisjp,psisjm,psiskp,psiskm,alphip,alphim,alphjp,alphjm,alphkp,alphkm,usip,usim,vsjp,vsjm,wskp,wskm,dsdip,dsdim,dsdjp,dsdjm,dsdkp,dsdkm,dsdtd_xy,dsdtd_z,dsdta,dsdt_s) async(1)
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(usip,usim,vsjp,vsjm,wskp,wskm) &
    !$OMP PRIVATE(dsdxp,dsdxm,dsdyp,dsdym,dsdzp,dsdzm) &
    !$OMP SHARED(nx,ny,nz,dxi,dyi,dzi,visc,u,v,w,s,dsdt,dzci,dzfi)
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
         psis   = psi(i,j,k)
         psisip = psi(i+1,j,k)
         psisim = psi(i-1,j,k)
         psisjp = psi(i,j+1,k)
         psisjm = psi(i,j-1,k)
         psiskp = psi(i,j,k+1)
         psiskm = psi(i,j,k-1)
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
         if (psis == 0.) then ! if not in solid
          if (psisip == 1.) s(i+1,j,k) = 2._rp*solidtemp; s(i,j,k) = 2._rp*s(i,j,k)
          if (psisim == 1.) s(i-1,j,k) = 2._rp*solidtemp; s(i,j,k) = 2._rp*s(i,j,k)
          if (psisjp == 1.) s(i,j+1,k) = 2._rp*solidtemp; s(i,j,k) = 2._rp*s(i,j,k)
          if (psisjm == 1.) s(i,j-1,k) = 2._rp*solidtemp; s(i,j,k) = 2._rp*s(i,j,k)
          if (psiskp == 1.) s(i,j,k+1) = 2._rp*solidtemp; s(i,j,k) = 2._rp*s(i,j,k)
          if (psiskm == 1.) s(i,j,k-1) = 2._rp*solidtemp; s(i,j,k) = 2._rp*s(i,j,k)
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
         if (psis == 0.) then ! if not in solid
          if (psisip == 1.) usip = solidtemp*u(i ,j,k)
          if (psisim == 1.) usim = solidtemp*u(i-1,j,k)
          if (psisjp == 1.) vsjp = solidtemp*v(i,j ,k)
          if (psisjm == 1.) vsjm = solidtemp*v(i,j-1,k)
          if (psiskp == 1.) wskp = solidtemp*w(i,j,k )
          if (psiskm == 1.) wskm = solidtemp*w(i,j,k-1)
#endif
#if !defined(_WENO) & defined(_IBM) && defined(_SIMPLE) && !defined(_ISOTHERMAL)
         dsdtd_xy = (alphip*dsdip-alphim*dsdim)*dxi + &
                    (alphjp*dsdjp-alphjm*dsdjm)*dyi
         dsdtd_z  = (alphkp*dsdkp-alphkm*dsdkm)*dzfi(k)
         dsdta    = - ( usip - usim )*dxi &    
                    - ( vsjp - vsjm )*dyi &    
                    - ( wskp - wskm )*dzfi(k)
#elif !defined(_WENO) & !defined(_IBM)
         dsdtd_xy = alph_f*(dsdip-dsdim)*dxi + &
                    alph_f*(dsdjp-dsdjm)*dyi
         dsdtd_z  = alph_f*(dsdkp-dsdkm)*dzfi(k)
         dsdta    = - ( usip - usim )*dxi &    
                    - ( vsjp - vsjm )*dyi &    
                    - ( wskp - wskm )*dzfi(k)
#elif defined(_WENO)
         dsdtd_xy = alph_f*(dsdip-dsdim)*dxi + &
                    alph_f*(dsdjp-dsdjm)*dyi
         dsdtd_z  = alph_f*(dsdkp-dsdkm)*dzfi(k)
         dsdta    = + dsdt(i,j,k)
#endif
         dsdt_s = dsdta + dsdtd_xy + dsdtd_z
         dsdt(i,j,k) = dsdt_s
        end do
      end do
    end do
    !
  end subroutine scal
end module mod_scal
