!
! SPDX-License-Identifier: MIT
!
module mod_weno
  !
  use mod_types, only: rp
  !
  implicit none
  !
  private
  public  :: weno5
  !
  contains
  !
  subroutine weno5(nx,ny,nz,nh_s,dxi,dyi,dzi,tmp,ux,uy,uz,dphidt)
    !
    implicit none
    !
    real(rp), parameter :: c11    =  2._rp/6._rp, &
                           c21    = -7._rp/6._rp, & 
                           c31    = 11._rp/6._rp, & 
                           c12    = -1._rp/6._rp, & 
                           c22    =  5._rp/6._rp, & 
                           c32    =  2._rp/6._rp, & 
                           c13    =  2._rp/6._rp, & 
                           c23    =  5._rp/6._rp, & 
                           c33    = -1._rp/6._rp 
    real(rp), parameter :: sigma1 = 0.1_rp, &
                           sigma2 = 0.6_rp, & 
                           sigma3 = 0.3_rp 
    real(rp), parameter :: eps    = 10._rp**(-6)
    !
    integer , intent(in )                      :: nx,ny,nz,nh_s
    real(rp), intent(in )                      :: dxi,dyi,dzi
    real(rp), intent(in ), dimension(1-nh_s:,1-nh_s:,1-nh_s:) :: tmp
    real(rp), intent(in ), dimension(0:,0:,0:) :: ux,uy,uz
    real(rp), intent(out), dimension(1:,1:,1:) :: dphidt
    !
    real(rp) :: fm2,fm1,f0,fp1,fp2
    real(rp) :: beta1,beta2,beta3
    real(rp) :: we1,we2,we3,sum_we
    real(rp) :: dfdlh1,dfdlh2,dfdlh3
    real(rp) :: dphidx,dphidy,dphidz
    integer  :: a,i,j,k,p
    real(rp) :: uxc,uyc,uzc
    !
    !$acc parallel loop collapse(3) default(present) private(fm2,fm1,f0,fp1,fp2,beta1,beta2,beta3,we1,we2,we3,sum_we,dfdlh1,dfdlh2,dfdlh3,dphidx,dphidy,dphidz,a) async(1)
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i,j,k) &
    !$OMP PRIVATE(uxc,uyc,uzc,a) &
    !$OMP PRIVATE(fm2,fm1,f0,fp1,fp2) &
    !$OMP PRIVATE(beta1,beta2,beta3,we1,we2,we3,eps) &
    !$OMP PRIVATE(c11,c21,c31,c12,c22,c32,c13,c23,c33) &
    !$OMP PRIVATE(dfdlh1,dfdlh2,dfdlh3,dphidx,dphidy,dphidz) &
    !$OMP SHARED(nx,ny,nz,dxi,dyi,dzi,ux,uy,uz,tmp,dphidt)
    do k=1,nz
      do j=1,ny
        do i=1,nx
          !
          uxc = 0.5_rp*(ux(i-1,j,k)+ux(i,j,k))
          a   = nint(sign(1._rp,uxc))
          fm2 = a*(tmp(i-2*a,j,k) - tmp(i-3*a,j,k))*dxi
          fm1 = a*(tmp(i-1*a,j,k) - tmp(i-2*a,j,k))*dxi
          f0  = a*(tmp(i+0*a,j,k) - tmp(i-1*a,j,k))*dxi
          fp1 = a*(tmp(i+1*a,j,k) - tmp(i+0*a,j,k))*dxi
          fp2 = a*(tmp(i+2*a,j,k) - tmp(i+1*a,j,k))*dxi
          beta1 = (13._rp/12._rp)*(      fm2 - 2._rp*fm1 +       f0 )**2 + &
                  ( 1._rp/4._rp )*(      fm2 - 4._rp*fm1 + 3._rp*f0 )**2
          beta2 = (13._rp/12._rp)*(      fm1 - 2._rp*f0  +       fp1)**2 + &
                  ( 1._rp/4._rp )*(      fm1             -       fp1)**2
          beta3 = (13._rp/12._rp)*(      f0  - 2._rp*fp1 +       fp2)**2 + &
                  ( 1._rp/4._rp )*(3._rp*f0  - 4._rp*fp1 +       fp2)**2
          we1 = sigma1*(1._rp+(abs(beta1-beta3)/(eps+beta1)))
          we2 = sigma2*(1._rp+(abs(beta1-beta3)/(eps+beta2)))
          we3 = sigma3*(1._rp+(abs(beta1-beta3)/(eps+beta3)))
          sum_we = we1+we2+we3
          we1 = we1/sum_we
          we2 = we2/sum_we
          we3 = we3/sum_we
          dfdlh1 = c11*fm2+c21*fm1+c31*f0
          dfdlh2 = c12*fm1+c22*f0 +c32*fp1
          dfdlh3 = c13*f0 +c23*fp1+c33*fp2
          dphidx = we1*dfdlh1+we2*dfdlh2+we3*dfdlh3
          !
          uyc = 0.5_rp*(uy(i,j-1,k)+uy(i,j,k))
          a   = nint(sign(1._rp,uyc))
          fm2 = a*(tmp(i,j-2*a,k) - tmp(i,j-3*a,k))*dyi
          fm1 = a*(tmp(i,j-1*a,k) - tmp(i,j-2*a,k))*dyi
          f0  = a*(tmp(i,j+0*a,k) - tmp(i,j-1*a,k))*dyi
          fp1 = a*(tmp(i,j+1*a,k) - tmp(i,j+0*a,k))*dyi
          fp2 = a*(tmp(i,j+2*a,k) - tmp(i,j+1*a,k))*dyi
          beta1 = (13._rp/12._rp)*(      fm2 - 2._rp*fm1 +       f0 )**2 + &
                  ( 1._rp/4._rp )*(      fm2 - 4._rp*fm1 + 3._rp*f0 )**2
          beta2 = (13._rp/12._rp)*(      fm1 - 2._rp*f0  +       fp1)**2 + &
                  ( 1._rp/4._rp )*(      fm1             -       fp1)**2
          beta3 = (13._rp/12._rp)*(      f0  - 2._rp*fp1 +       fp2)**2 + &
                  ( 1._rp/4._rp )*(3._rp*f0  - 4._rp*fp1 +       fp2)**2
          we1 = sigma1*(1._rp+(abs(beta1-beta3)/(eps+beta1)))
          we2 = sigma2*(1._rp+(abs(beta1-beta3)/(eps+beta2)))
          we3 = sigma3*(1._rp+(abs(beta1-beta3)/(eps+beta3)))
          sum_we = we1+we2+we3
          we1 = we1/sum_we
          we2 = we2/sum_we
          we3 = we3/sum_we
          dfdlh1 = c11*fm2+c21*fm1+c31*f0
          dfdlh2 = c12*fm1+c22*f0 +c32*fp1
          dfdlh3 = c13*f0 +c23*fp1+c33*fp2
          dphidy = we1*dfdlh1+we2*dfdlh2+we3*dfdlh3
          !
          uzc = 0.5_rp*(uz(i,j,k-1)+uz(i,j,k))
          a   = nint(sign(1._rp,uzc))
          fm2 = a*(tmp(i,j,k-2*a) - tmp(i,j,k-3*a))*dzi
          fm1 = a*(tmp(i,j,k-1*a) - tmp(i,j,k-2*a))*dzi
          f0  = a*(tmp(i,j,k+0*a) - tmp(i,j,k-1*a))*dzi
          fp1 = a*(tmp(i,j,k+1*a) - tmp(i,j,k+0*a))*dzi
          fp2 = a*(tmp(i,j,k+2*a) - tmp(i,j,k+1*a))*dzi
          beta1 = (13._rp/12._rp)*(      fm2 - 2._rp*fm1 +       f0 )**2 + &
                  ( 1._rp/4._rp )*(      fm2 - 4._rp*fm1 + 3._rp*f0 )**2
          beta2 = (13._rp/12._rp)*(      fm1 - 2._rp*f0  +       fp1)**2 + &
                  ( 1._rp/4._rp )*(      fm1             -       fp1)**2
          beta3 = (13._rp/12._rp)*(      f0  - 2._rp*fp1 +       fp2)**2 + &
                  ( 1._rp/4._rp )*(3._rp*f0  - 4._rp*fp1 +       fp2)**2
          we1 = sigma1*(1._rp+(abs(beta1-beta3)/(eps+beta1)))
          we2 = sigma2*(1._rp+(abs(beta1-beta3)/(eps+beta2)))
          we3 = sigma3*(1._rp+(abs(beta1-beta3)/(eps+beta3)))
          sum_we = we1+we2+we3
          we1 = we1/sum_we
          we2 = we2/sum_we
          we3 = we3/sum_we
          dfdlh1 = c11*fm2+c21*fm1+c31*f0
          dfdlh2 = c12*fm1+c22*f0 +c32*fp1
          dfdlh3 = c13*f0 +c23*fp1+c33*fp2
          dphidz = we1*dfdlh1+we2*dfdlh2+we3*dfdlh3
          !
          dphidt(i,j,k) = - (uxc*dphidx + uyc*dphidy + uzc*dphidz)
          !
        enddo
      enddo
    enddo
    !
    return
  end subroutine weno5
end module mod_weno
