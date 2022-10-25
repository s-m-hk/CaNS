! -
!
! SPDX-FileCopyrightText: Copyright (c) 2017-2022 Pedro Costa and the CaNS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
module mod_scal
  use mod_types
  implicit none
  private
  public scal
  contains
  subroutine scal(nx,ny,nz,dxi,dyi,dzi,dzci,dzfi,alph_f,alph_s, &
#if defined(_IBM)
                  psi, &
#endif
                  u,v,w,s,dsdt)
    !
    !
    !
    implicit none
    integer , intent(in) :: nx,ny,nz
    real(rp), intent(in) :: dxi,dyi,dzi,alph_f,alph_s
    real(rp), intent(in), dimension(0:) :: dzci,dzfi
    real(rp), dimension(0:,0:,0:), intent(in) :: u,v,w,s
#if defined(_IBM)
    real(rp), intent(in), dimension(0:,0:,0:) :: psi
#endif
    real(rp), dimension(:,:,:), intent(out) :: dsdt
    integer :: i,j,k
    real(rp) :: usip,usim,vsjp,vsjm,wskp,wskm
    real(rp) :: dsdxp,dsdxm,dsdyp,dsdym,dsdzp,dsdzm
    !
    !$acc parallel loop collapse(3) default(present) private(usip,usim,vsjp,vsjm,wskp,wskm,dsdxp,dsdxm,dsdyp,dsdym,dsdzp,dsdzm) async(1)
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(usip,usim,vsjp,vsjm,wskp,wskm) &
    !$OMP PRIVATE(dsdxp,dsdxm,dsdyp,dsdym,dsdzp,dsdzm) &
    !$OMP SHARED(nx,ny,nz,dxi,dyi,dzi,visc,u,v,w,s,dsdt,dzci,dzfi)
    do k=1,nz
      do j=1,ny
        do i=1,nx
#if defined(_IBM) && defined(_SIMPLE) && defined(_ISOTHERMAL)
          if (psi(i,j,k) == 0.) then ! if not in solid
              if (psi(i+1,j,k) == 1.) s(i+1,j,k) = 2.*solidtemp; s(i,j,k) = 2.*s(i,j,k)
              if (psi(i-1,j,k) == 1.) s(i-1,j,k) = 2.*solidtemp; s(i,j,k) = 2.*s(i,j,k)
              if (psi(i,j+1,k) == 1.) s(i,j+1,k) = 2.*solidtemp; s(i,j,k) = 2.*s(i,j,k)
              if (psi(i,j-1,k) == 1.) s(i,j-1,k) = 2.*solidtemp; s(i,j,k) = 2.*s(i,j,k)
              if (psi(i,j,k+1) == 1.) s(i,j,k+1) = 2.*solidtemp; s(i,j,k) = 2.*s(i,j,k)
              if (psi(i,j,k-1) == 1.) s(i,j,k-1) = 2.*solidtemp; s(i,j,k) = 2.*s(i,j,k)
          endif
#endif
          usim  = 0.5*( s(i-1,j,k)+s(i,j,k) )*u(i-1,j,k)
          usip  = 0.5*( s(i+1,j,k)+s(i,j,k) )*u(i  ,j,k)
          vsjm  = 0.5*( s(i,j-1,k)+s(i,j,k) )*v(i,j-1,k)
          vsjp  = 0.5*( s(i,j+1,k)+s(i,j,k) )*v(i,j  ,k)
          wskm  = 0.5*( s(i,j,k-1)+s(i,j,k) )*w(i,j,k-1)
          wskp  = 0.5*( s(i,j,k+1)+s(i,j,k) )*w(i,j,k  )
#if defined(_IBM) && defined(_SIMPLE) && defined(_ISOTHERMAL)
          if (psi(i,j,k) == 0.) then ! if not in solid
              if (psi(i+1,j,k) == 1.) usip = solidtemp*u(i ,j,k)
              if (psi(i-1,j,k) == 1.) usim = solidtemp*u(i-1,j,k)
              if (psi(i,j+1,k) == 1.) vsjp = solidtemp*v(i,j ,k)
              if (psi(i,j-1,k) == 1.) vsjm = solidtemp*v(i,j-1,k)
              if (psi(i,j,k+1) == 1.) wskp = solidtemp*w(i,j,k )
              if (psi(i,j,k-1) == 1.) wskm = solidtemp*w(i,j,k-1)
          endif
#endif
          dsdxp = (s(i+1,j,k)-s(i  ,j,k))*dxi
          dsdxm = (s(i  ,j,k)-s(i-1,j,k))*dxi
          dsdyp = (s(i,j+1,k)-s(i,j  ,k))*dyi
          dsdym = (s(i,j  ,k)-s(i,j-1,k))*dyi
          dsdzp = (s(i,j,k+1)-s(i,j,k  ))*dzci(k  )
          dsdzm = (s(i,j,k  )-s(i,j,k-1))*dzci(k-1)
          !
#if defined(_IBM) && defined(_SIMPLE)
          dsdt(i,j,k) = dxi*(     -usip + usim ) + (dsdxp-dsdxm)* &
                        ( (1. - psi(i,j,k)) * alph_f + psi(i,j,k) * alph_s )*dxi + &
                        dyi*(     -vsjp + vsjm ) + (dsdyp-dsdym)* &
                        ( (1. - psi(i,j,k)) * alph_f + psi(i,j,k) * alph_s )*dyi + &
                        dzfi(k)*( -wskp + wskm ) + (dsdzp-dsdzm)* &
                        ( (1. - psi(i,j,k)) * alph_f + psi(i,j,k) * alph_s )*dzfi(k)
#else
          dsdt(i,j,k) = dxi*(     -usip + usim ) + (dsdxp-dsdxm)*alph_f*dxi + &
                        dyi*(     -vsjp + vsjm ) + (dsdyp-dsdym)*alph_f*dyi + &
                        dzfi(k)*( -wskp + wskm ) + (dsdzp-dsdzm)*alph_f*dzfi(k)
#endif
        end do
      end do
    end do
  end subroutine scal
end module mod_scal
