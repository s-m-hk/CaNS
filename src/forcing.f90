#if defined(_IBM)
module mod_forcing
  use mpi
  use mod_types
  use mod_common_mpi, only: ierr
  implicit none
  private
  public force_vel,force_scal,force_bulk_vel,bulk_mean_ibm
  contains
  subroutine force_vel(n,dl,dzc,dzf,l,psi_u,psi_v,psi_w,u,v,w,fibm)
    !
    ! Force velocity field using volume-penalization IBM:
    ! the force is proportional to the volume fraction of
    ! solid in a grid cell i,j,k.
    ! The volume fraction is defined in the cell centers and
    ! is interpolated to the cell faces where the velocity
    ! components are defined. This results in a smoothing of the
    ! local volume fraction field.
    ! 
    ! Note: for some systems it may be convenient to save the force
    !       distribution
    !
    implicit none
    integer , intent(in), dimension(3) :: n
    real(rp), intent(in   ), dimension(3) :: dl,l
    real(rp), intent(in   ), dimension(0:) :: dzc,dzf
    real(rp), intent(in   ), dimension(0:,0:,0:) :: psi_u,psi_v,psi_w
    real(rp), intent(inout), dimension(0:,0:,0:) :: u,v,w
    real(rp), intent(out  ), dimension(4) :: fibm
    real(rp) :: psix,psiy,psiz,fx,fy,fz,fxtot,fytot,fztot,dx,dy
    integer :: i,j,k,nx,ny,nz
    !
    nx = n(1)
    ny = n(2)
    nz = n(3)
    dx = dl(1)
    dy = dl(2)
    fxtot = 0._rp
    fytot = 0._rp
    fztot = 0._rp
    !$acc data copy(fxtot,fytot,fztot) async(1)
    !$acc parallel loop collapse(3) default(present) private(psix,psiy,psiz,fx,fy,fz) reduction(+:fxtot,fytot,fztot) async(1)
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP SHARED(n,u,v,w,psi_u,psi_v,psi_w,dl,dzc,dzf) &
    !$OMP PRIVATE(i,j,k,fx,fy,fz,psix,psiy,psiz) &
    !$OMP REDUCTION(+:fxtot,fytot,fztot)
    do k=1,nz
      do j=1,ny
        do i=1,nx
          psix  = psi_u(i,j,k)
          psiy  = psi_v(i,j,k)
          psiz  = psi_w(i,j,k)
          fx    = - u(i,j,k)*psix ! (u(i,j,k)*(1.-psix)-u(i,j,k))*dti
          fy    = - v(i,j,k)*psiy ! (v(i,j,k)*(1.-psiy)-v(i,j,k))*dti
          fz    = - w(i,j,k)*psiz ! (w(i,j,k)*(1.-psiz)-w(i,j,k))*dti
          u(i,j,k) = u(i,j,k) + fx
          v(i,j,k) = v(i,j,k) + fy
          w(i,j,k) = w(i,j,k) + fz
          fxtot = fxtot + fx*dx*dy*dzf(k)
          fytot = fytot + fy*dx*dy*dzf(k)
          fztot = fztot + fz*dx*dy*dzc(k)
        enddo
      enddo
    enddo
    !$acc end data
    !$acc wait(1)
    call MPI_ALLREDUCE(MPI_IN_PLACE,fxtot,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE,fytot,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE,fztot,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    fibm(1) = fxtot/(l(1)*l(2)*l(3))
    fibm(2) = fytot/(l(1)*l(2)*l(3))
    fibm(3) = fztot/(l(1)*l(2)*l(3))
  end subroutine force_vel
  !
  subroutine force_scal(n,nh_s,dl,dz,l,psi,s,fibm)
    !
    implicit none
    integer , intent(in), dimension(3) :: n
    integer , intent(in)               :: nh_s
    real(rp), intent(in   ) , dimension(3) :: dl,l
    real(rp), intent(in   ) , dimension(0:) :: dz
    real(rp), intent(in   ) , dimension(0:,0:,0:) :: psi
    real(rp), intent(inout) , dimension(1-nh_s:,1-nh_s:,1-nh_s:) :: s
    real(rp), intent(inout  ), dimension(4) :: fibm
    real(rp) :: psis,fs,fstot,dx,dy
    integer :: i,j,k,nx,ny,nz
    !
    nx = n(1)
    ny = n(2)
    nz = n(3)
    dx = dl(1)
    dy = dl(2)
    fstot = 0._rp
    !$acc data copy(fstot) async(1)
    !$acc parallel loop collapse(3) default(present) private(psis,fs) reduction(+:fstot) async(1)
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP SHARED(n,s,psi,dl,dzc,dzf) &
    !$OMP PRIVATE(i,j,k,fs,psis) &
    !$OMP REDUCTION(+:fstot)
    do k=1,nz
      do j=1,ny
        do i=1,nx
          psis  = psi(i,j,k)
          fs    = - s(i,j,k)*psis
          s(i,j,k) = s(i,j,k) + fs
          fstot = fstot + fs*dx*dy*dz(k)
        enddo
      enddo
    enddo
    !$acc end data
    !$acc wait(1)
    call MPI_ALLREDUCE(MPI_IN_PLACE,fstot,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    fibm(4) = fstot/(l(1)*l(2)*l(3))
  end subroutine force_scal
  !
  subroutine force_bulk_vel(n,nh,dl,dz,l,psi,p,velf,f)
    !
    ! bulk velocity forcing only in a region of the domain
    ! where psi is non-zero
    !
    implicit none
    integer , intent(in   ), dimension(3) :: n
    integer , intent(in   ) :: nh
    real(rp), intent(in   ) , dimension(3) :: dl,l
    real(rp), intent(in   ) , dimension(0:) :: dz
    real(rp), intent(in   ), dimension(0:,0:,0:) :: psi
    real(rp), intent(inout), dimension(1-nh:,1-nh:,1-nh:) :: p
    real(rp), intent(in   ) :: velf
    real(rp), intent(inout  ) :: f
    integer :: i,j,k,nx,ny,nz
    real(rp) :: psis,mean_val,mean_psi,dx,dy
    !
    nx = n(1)
    ny = n(2)
    nz = n(3)
    dx = dl(1)
    dy = dl(2)
    mean_val = 0._rp
    mean_psi = 0._rp
    !
    !$acc data copy(mean_val,mean_psi) async(1)
    !$acc parallel loop collapse(3) default(present) private(psis) reduction(+:mean_val) reduction(+:mean_psi) async(1)
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP SHARED(n,psi,p,dl,dz) &
    !$OMP PRIVATE(i,j,k,psis) &
    !$OMP REDUCTION(+:mean_val,mean_psi)
    do k=1,nz
      do j=1,ny
        do i=1,nx
          psis = 1.0_rp - psi(i,j,k)
          mean_val = mean_val + p(i,j,k)*psis*dx*dy*dz(k)
          mean_psi = mean_psi + psis*dx*dy*dz(k)
        enddo
      enddo
    enddo
    !$acc end data
    !$acc wait(1)
    call mpi_allreduce(MPI_IN_PLACE,mean_val,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    call mpi_allreduce(MPI_IN_PLACE,mean_psi,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    mean_val = mean_val/mean_psi
    !$acc parallel loop collapse(3) default(present) private(psis) async(1)
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP SHARED(n,psi,p,mean_val,velf,dl,dz) &
    !$OMP PRIVATE(i,j,k,psis) &
    !$OMP REDUCTION(+:f)
    do k=1,nz
      do j=1,ny
        do i=1,nx
#if defined(_FORCE_FLUID_ONLY)
          psis = 1._rp-psi(i,j,k) ! (if bulk velocity forced only inside the fluid)
#else
          psis = 1._rp
#endif
          p(i,j,k) = p(i,j,k) + (velf-mean_val)*psis
          f = f + (velf-mean_val)*psis*dx*dy*dz(k)
        enddo
      enddo
    enddo
    !$acc wait(1)
    call mpi_allreduce(MPI_IN_PLACE,f,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
#if defined(_FORCE_FLUID_ONLY)
    f = f/mean_psi ! (if bulk velocity forced only inside the fluid)
#else
    f = f/(l(1)*l(2)*l(3))
#endif
  end subroutine force_bulk_vel
  !
  subroutine bulk_mean_ibm(n,dl,dz,psi,p,mean)
    !
    implicit none
    integer , intent(in   ), dimension(3) :: n
    real(rp), intent(in   ) , dimension(3) :: dl
    real(rp), intent(in   ) , dimension(0:) :: dz
    real(rp), intent(in   ), dimension(0:,0:,0:) :: psi
    real(rp), intent(inout), dimension(0:,0:,0:) :: p
    real(rp), intent(out  ) :: mean
    real(rp)              :: mean_val,mean_psi
    integer :: i,j,k,nx,ny,nz
    real(rp) :: psis,dx,dy
    !
    nx = n(1)
    ny = n(2)
    nz = n(3)
    dx = dl(1)
    dy = dl(2)
    mean_val = 0._rp
    mean_psi = 0._rp
    !
    !$acc data copy(mean_val,mean_psi) async(1)
    !$acc parallel loop collapse(3) default(present) private(psis) reduction(+:mean_val) reduction(+:mean_psi) async(1)
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP SHARED(n,psi,p,dl,dz) &
    !$OMP PRIVATE(i,j,k,psis) &
    !$OMP REDUCTION(+:mean_val,mean_psi)
    do k=1,nz
      do j=1,ny
        do i=1,nx
          psis = 1._rp-psi(i,j,k)
          mean_val = mean_val + p(i,j,k)*psis*dx*dy*dz(k)
          mean_psi = mean_psi + psis*dx*dy*dz(k)
        enddo
      enddo
    enddo
    !$acc end data
    !$acc wait(1)
    call mpi_allreduce(MPI_IN_PLACE,mean_val,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    call mpi_allreduce(MPI_IN_PLACE,mean_psi,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    mean = mean_val/mean_psi
  end subroutine bulk_mean_ibm
end module mod_forcing
#endif