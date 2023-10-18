! -
!
! SPDX-FileCopyrightText: Copyright (c) 2017-2022 Pedro Costa and the CaNS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
module mod_output
  use mpi
  use decomp_2d_io
  use mod_common_mpi, only:ierr,myid
  use mod_types
  implicit none
  private
  public out0d,gen_alias,out1d,out1d_chan,out1d_chan_tmp,out2d,out3d,write_log_output,write_visu_2d,write_visu_3d
  contains
  subroutine out0d(fname,n,var)
    !
    ! appends the first n entries of an array
    ! var to a file
    ! fname -> name of the file
    ! n     -> number of entries
    ! var   -> input array of real values
    !
    implicit none
    character(len=*), intent(in) :: fname
    integer , intent(in) :: n
    real(rp), intent(in), dimension(:) :: var
    integer :: iunit
    !
    if (myid  ==  0) then
      open(newunit=iunit,file=fname,position='append')
      write(iunit,'(*(E16.7e3))') var(1:n)
      close(iunit)
    end if
  end subroutine out0d
  !
  subroutine gen_alias(myid,datadir,fname,fname_alias)
    !
    ! this subroutine generates a symlink with name `fname_alias`, pointing to
    ! file `datadir//fname` using the `execute_command_line` intrinsic;
    ! it is called by task `myid`
    !
    integer, intent(in) :: myid
    character(len=*), intent(in) :: datadir,fname,fname_alias
    if(myid == 0) call execute_command_line('ln -sf '//trim(fname)//' '//trim(datadir)//fname_alias)
  end subroutine gen_alias
  !
  subroutine out1d(fname,ng,lo,hi,idir,l,dl,z_g,dz,p)
    !
    ! writes the profile of a variable averaged
    ! over two domain directions
    !
    ! fname -> name of the file
    ! ng    -> global domain sizes
    ! lo,hi -> upper and lower extents of the input array
    ! idir  -> direction of the profile
    ! dl,dl -> uniform grid spacing and length arrays
    ! z_g   -> global z coordinate array (grid is non-uniform in z)
    ! dz    -> local z grid spacing array (should work also with the global one)
    ! p     -> 3D input scalar field
    !
    implicit none
    character(len=*), intent(in) :: fname
    integer , intent(in), dimension(3) :: ng,lo,hi
    integer , intent(in) :: idir
    real(rp), intent(in), dimension(3) :: l,dl
    real(rp), intent(in), dimension(0:       ) :: z_g
    real(rp), intent(in), dimension(lo(3)-1:) :: dz
    real(rp), intent(in), dimension(lo(1)-1:,lo(2)-1:,lo(3)-1:) :: p
    real(rp), allocatable, dimension(:) :: p1d
    integer :: i,j,k
    integer :: iunit
    real(rp) :: grid_area_ratio,p1d_s
    !
    allocate(p1d(ng(idir)))
    !$acc enter data create(p1d)
    !$acc kernels default(present)
    p1d(:) = 0._rp
    !$acc end kernels
    select case(idir)
    case(3)
      grid_area_ratio = dl(1)*dl(2)/(l(1)*l(2))
      !$acc parallel loop gang default(present) private(p1d_s)
      do k=lo(3),hi(3)
        p1d_s = 0._rp
        !$acc loop collapse(2) reduction(+:p1d_s)
        do j=lo(2),hi(2)
          do i=lo(1),hi(1)
            p1d_s = p1d_s + p(i,j,k)*grid_area_ratio
          end do
        end do
        p1d(k) = p1d_s
      end do
      !$acc exit data copyout(p1d)
      call MPI_ALLREDUCE(MPI_IN_PLACE,p1d(1),ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      if(myid == 0) then
        open(newunit=iunit,file=fname)
        do k=1,ng(3)
          write(iunit,'(2E16.7e3)') z_g(k),p1d(k)
        end do
        close(iunit)
      end if
    case(2)
      grid_area_ratio = dl(1)/(l(1)*l(3))
      !$acc parallel loop gang default(present) private(p1d_s)
      do j=lo(2),hi(2)
        p1d_s = 0._rp
        !$acc loop collapse(2) reduction(+:p1d_s)
        do k=lo(3),hi(3)
          do i=lo(1),hi(1)
            p1d_s = p1d_s + p(i,j,k)*dz(k)*grid_area_ratio
          end do
        end do
        p1d(j) = p1d_s
      end do
      !$acc exit data copyout(p1d)
      call MPI_ALLREDUCE(MPI_IN_PLACE,p1d(1),ng(2),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      if(myid == 0) then
        open(newunit=iunit,file=fname)
        do j=1,ng(2)
          write(iunit,'(2E16.7e3)') (j-.5)*dl(2),p1d(j)
        end do
        close(iunit)
      end if
    case(1)
      grid_area_ratio = dl(2)/(l(2)*l(3))
      !$acc parallel loop gang default(present) private(p1d_s)
      do i=lo(1),hi(1)
        p1d_s = 0._rp
        !$acc loop collapse(2) reduction(+:p1d_s)
        do k=lo(3),hi(3)
          do j=lo(2),hi(2)
            p1d_s = p1d_s + p(i,j,k)*dz(k)*grid_area_ratio
          end do
        end do
        p1d(i) = p1d_s
      end do
      !$acc exit data copyout(p1d)
      call MPI_ALLREDUCE(MPI_IN_PLACE,p1d(1),ng(1),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      if(myid == 0) then
        open(newunit=iunit,file=fname)
        do i=1,ng(1)
          write(iunit,'(2E16.7e3)') (i-.5)*dl(1),p1d(i)
        end do
        close(iunit)
      end if
    end select
  end subroutine out1d
  !
  subroutine out2d(fname,fname_fld,varname,lo,hi,inorm,islice,p)
    use mod_common_mpi, only: ipencil => ipencil_axis
    !
    ! saves a planar slice of a scalar field into a binary file
    !
    ! fname  -> name of the output file
    ! inorm  -> plane is perpendicular to direction
    !           inorm (1,2,3)
    ! islice -> plane is of constant index islice
    !           in direction inorm
    ! p      -> 3D input scalar field
    !
    implicit none
    character(len=*), intent(in) :: fname,fname_fld,varname
    integer, dimension(3) :: lo,hi
    integer , intent(in) :: inorm,islice
    real(rp), intent(in), dimension(:,:,:) :: p
    !
    ! masked in case of _SINGLE_PRECISION_POISSON since 2DECOMP does not yet support two precisions
    !
#if defined(_USE_HDF5)
    block
      use decomp_2d
      use mod_load, only: io_field_hdf5
      integer, dimension(3) :: ng
      !
      select case(inorm)
      case(1)
       ng(:) = [islice,ny_global,nz_global]
      case(2)
       ng(:) = [nx_global,islice,nz_global]
      case(3)
       ng(:) = [nx_global,ny_global,islice]
      end select
      !
      call io_field_hdf5('w',fname,varname,ng,[1,1,1],lo,hi,p)
    end block
#else
    select case(inorm)
    case(1) !normal to x --> yz plane
       call decomp_2d_write_plane(ipencil,p,inorm,islice,'.',fname,'dummy')
    case(2) !normal to y --> zx plane
       call decomp_2d_write_plane(ipencil,p,inorm,islice,'.',fname,'dummy')
    case(3) !normal to z --> xy plane
       call decomp_2d_write_plane(ipencil,p,inorm,islice,'.',fname,'dummy')
    end select
#endif
  end subroutine out2d
  !
  subroutine out3d(fname,fname_fld,varname,nskip,p)
    use mod_common_mpi, only: ipencil => ipencil_axis
    !
    ! saves a 3D scalar field into a binary file
    !
    ! fname  -> name of the output file
    ! nskip  -> array with the step size for which the
    !           field is written; i.e.: [1,1,1]
    !           writes the full field
    ! p      -> 3D input scalar field
    !
    implicit none
    character(len=*), intent(in) :: fname,fname_fld,varname
    integer , intent(in), dimension(3) :: nskip
    real(rp), intent(in), dimension(:,:,:) :: p
    integer :: fh
    integer(kind=MPI_OFFSET_KIND) :: filesize,disp
    !
    ! masked in case of _SINGLE_PRECISION_POISSON since 2DECOMP does not yet support two precisions
    !
    call MPI_FILE_OPEN(MPI_COMM_WORLD, fname_fld, &
         MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL,fh, ierr)
    filesize = 0_MPI_OFFSET_KIND
    call MPI_FILE_SET_SIZE(fh,filesize,ierr)
    disp = 0_MPI_OFFSET_KIND
#if !defined(_OPENACC)
    call decomp_2d_write_every(ipencil,p,nskip(1),nskip(2),nskip(3),fname_fld,.true.)
#else
    !
    ! temporary workaround due to 2DECOMP's lack of support for different real kinds in a single build
    !
    block
      use decomp_2d
#if defined(_USE_HDF5)
      use mod_load, only: io_field, io_field_hdf5
#else
      use mod_load, only: io_field
#endif
      integer, dimension(3) :: ng,lo,hi
      ng(:) = [nx_global,ny_global,nz_global]
      select case(ipencil)
      case(1)
        lo(:) = xstart(:)
        hi(:) = xend(:)
      case(2)
        lo(:) = ystart(:)
        hi(:) = yend(:)
      case(3)
        lo(:) = zstart(:)
        hi(:) = zend(:)
      end select
      if(any(nskip /= 1) .and. myid == 0) &
        print*, 'Warning: `nskip` should be `[1,1,1]` if `io_field()` is used to output 3D field data'
#if defined(_USE_HDF5)
      call io_field_hdf5('w',fname,varname,ng,[1,1,1],lo,hi,p)
#else
      call io_field('w',fh,ng,[0,0,0],lo,hi,disp,p)
#endif
    end block
#endif
    call MPI_FILE_CLOSE(fh,ierr)
  end subroutine out3d
  !
  subroutine write_log_output(fname,fname_fld,varname,nmin,nmax,nskip,time,istep)
    !
    ! appends information about a saved binary file to a log file
    ! this file is used to generate a xdmf file for visualization of field data
    !
    ! fname     -> name of the output log file
    ! fname_fld -> name of the saved binary file (excluding the directory)
    ! varname   -> name of the variable that is saved
    ! nmin      -> first element of the field that is saved in each direction, e.g. [1,1,1]
    ! nmax      -> last  element of the field that is saved in each direction, e.g. [ng(1),ng(2),ng(3)]
    ! nskip     -> step size between nmin and nmax, e.g. [1,1,1] if the whole array is saved
    ! time      -> physical time
    ! istep     -> time step number
    !
    implicit none
    character(len=*), intent(in) :: fname,fname_fld,varname
    integer , intent(in), dimension(3) :: nmin,nmax,nskip
    real(rp), intent(in)               :: time
    integer , intent(in)               :: istep
    character(len=100) :: cfmt
    integer :: iunit
    !
    write(cfmt, '(A)') '(A,A,A,9i5,E16.7e3,i7)'
    if (myid  ==  0) then
      open(newunit=iunit,file=fname,position='append')
      write(iunit,trim(cfmt)) trim(fname_fld),' ',trim(varname),nmin,nmax,nskip,time,istep
      close(iunit)
    end if
  end subroutine write_log_output
  !
  subroutine write_visu_3d(datadir,fname,fname_bin,fname_log,varname,nmin,nmax,nskip,time,istep,p)
    !
    ! wraps the calls of out3d and write_log_output into the same subroutine
    !
    implicit none
    character(len=*), intent(in)          :: datadir,fname,fname_bin,fname_log,varname
    integer , intent(in), dimension(3)    :: nmin,nmax,nskip
    real(rp), intent(in)                  :: time
    integer , intent(in)                  :: istep
    real(rp), intent(in), dimension(:,:,:) :: p
    !
    call out3d(trim(datadir)//trim(fname),trim(datadir)//trim(fname_bin),trim(varname),nskip,p)
    call write_log_output(trim(datadir)//trim(fname_log),trim(fname_bin),trim(varname),nmin,nmax,nskip,time,istep)
  end subroutine write_visu_3d
  !
  subroutine write_visu_2d(datadir,fname,fname_bin,fname_log,varname,inorm,nslice,ng,time,istep,p)
    !
    ! wraps the calls of out2d and write-log_output into the same subroutine
    !
    implicit none
    character(len=*), intent(in)          :: datadir,fname,fname_bin,fname_log,varname
    integer , intent(in)                  :: inorm,nslice
    integer , intent(in), dimension(3)    :: ng
    real(rp), intent(in)                  :: time
    integer , intent(in)                  :: istep
    real(rp), intent(in), dimension(:,:,:) :: p
    integer , dimension(3) :: nmin_2d,nmax_2d
    !
    select case(inorm)
    case(1)
      nmin_2d(:) = [nslice,1    ,1    ]
      nmax_2d(:) = [nslice,ng(2),ng(3)]
    case(2)
      nmin_2d(:) = [1    ,nslice,1    ]
      nmax_2d(:) = [ng(1),nslice,ng(3)]
    case(3)
      nmin_2d(:) = [1    ,1    ,nslice]
      nmax_2d(:) = [ng(1),ng(2),nslice]
    end select
    !
    call out2d(trim(datadir)//trim(fname),trim(datadir)//trim(fname_bin),trim(varname),nmin_2d,nmax_2d,inorm,nslice,p)
    call write_log_output(trim(datadir)//trim(fname_log),trim(fname_bin),trim(varname),nmin_2d,nmax_2d,[1,1,1],time,istep)
  end subroutine write_visu_2d
  !
  subroutine out1d_chan(fname,ng,lo,hi,idir,l,dl,dz,z_g,u,v,w,p,psi) ! e.g. for a channel with streamwise dir in x
    implicit none
    character(len=*), intent(in) :: fname
    integer , intent(in), dimension(3) :: ng,lo,hi
    integer , intent(in) :: idir
    real(rp), intent(in), dimension(3) :: l,dl
    real(rp), intent(in), dimension(0:) :: z_g
    real(rp), intent(in), dimension(lo(3)-1:) :: dz
    real(rp), intent(in), dimension(lo(1)-1:,lo(2)-1:,lo(3)-1:) :: u,v,w,p
    real(rp), intent(in), dimension(lo(1)-1:,lo(2)-1:,lo(3)-1:), optional :: psi
    real(rp), allocatable, dimension(:) :: um,vm,wm,pm,u2,v2,w2,p2,uw
#if defined(_IBM)
    real(rp), allocatable, dimension(:) :: ud,vd,wd,pd,ud2,vd2,wd2,pd2,udwd
    real(rp), allocatable, dimension(:) :: uf,vf,wf,pf,uf2,vf2,wf2,pf2,ufwf
    real(rp), allocatable, dimension(:) :: mean_psif
#endif
    integer :: i,j,k
    integer :: iunit
    integer :: q
    real(rp) :: grid_area_ratio
#if defined(_IBM)
    real(rp) :: dx,dy,psif
#endif
    !
    q = ng(idir)
    select case(idir)
    case(3)
      grid_area_ratio = dl(1)*dl(2)/(l(1)*l(2))
      allocate(um(0:q+1),vm(0:q+1),wm(0:q+1),pm(0:q+1),u2(0:q+1),v2(0:q+1),w2(0:q+1),p2(0:q+1),uw(0:q+1))
#if defined(_IBM)
      dx = dl(1)
      dy = dl(2)
      allocate(ud(0:q+1),vd(0:q+1),wd(0:q+1),pd(0:q+1),ud2(0:q+1),vd2(0:q+1),wd2(0:q+1),pd2(0:q+1),udwd(0:q+1))
      allocate(uf(0:q+1),vf(0:q+1),wf(0:q+1),pf(0:q+1),uf2(0:q+1),vf2(0:q+1),wf2(0:q+1),pf2(0:q+1),ufwf(0:q+1))
      allocate(mean_psif(0:q+1))
#endif
      um(:) = 0.0_rp
      vm(:) = 0.0_rp
      wm(:) = 0.0_rp
      pm(:) = 0.0_rp
      u2(:) = 0.0_rp
      v2(:) = 0.0_rp
      w2(:) = 0.0_rp
      p2(:) = 0.0_rp
      uw(:) = 0.0_rp
#if defined(_IBM)
      ! Phase-averaged velocities and stresses (fluid phase)
      uf(:) = 0.0_rp
      vf(:) = 0.0_rp
      wf(:) = 0.0_rp
      pf(:) = 0.0_rp
      uf2(:) = 0.0_rp
      vf2(:) = 0.0_rp
      wf2(:) = 0.0_rp
      pf2(:) = 0.0_rp
      ufwf(:) = 0.0_rp
      mean_psif(:) = 0.0_rp
#endif
      do k=lo(3),hi(3)
        do j=lo(2),hi(2)
          do i=lo(1),hi(1)
#if defined(_IBM)
            psif = 1.0_rp - psi(i,j,k)
            mean_psif(k) = mean_psif(k) + psif*dx*dy*dz(k)
#endif
            um(k) = um(k) + u(i,j,k)
            vm(k) = vm(k) + v(i,j,k)
            wm(k) = wm(k) + 0.5_rp*(w(i,j,k-1) + w(i,j,k))
            pm(k) = pm(k) + p(i,j,k)
            u2(k) = u2(k) + u(i,j,k)**2
            v2(k) = v2(k) + v(i,j,k)**2
            w2(k) = w2(k) + 0.5_rp*(w(i,j,k)**2+w(i,j,k-1)**2)
            p2(k) = p2(k) + p(i,j,k)**2
            uw(k) = uw(k) + 0.25_rp*(u(i-1,j,k) + u(i,j,k))* &
                                    (w(i,j,k-1) + w(i,j,k))
#if defined(_IBM)
            uf(k) = uf(k)     + u(i,j,k)*psif*dx*dy*dz(k)
            vf(k) = vf(k)     + v(i,j,k)*psif*dx*dy*dz(k)
            wf(k) = wf(k)     + 0.5_rp*(w(i,j,k-1) + w(i,j,k))*psif*dx*dy*dz(k)
            pf(k) = pf(k)     + p(i,j,k)*psif*dx*dy*dz(k)
            uf2(k) = uf2(k)   + (u(i,j,k)**2)*psif*dx*dy*dz(k)
            vf2(k) = vf2(k)   + (v(i,j,k)**2)*psif*dx*dy*dz(k)
            wf2(k) = wf2(k)   + 0.5_rp*(w(i,j,k)**2+w(i,j,k-1)**2)*psif*dx*dy*dz(k)
            pf2(k) = pf2(k)   + (p(i,j,k)**2)*psif*dx*dy*dz(k)
            ufwf(k) = ufwf(k) + 0.25_rp*(u(i-1,j,k) + u(i,j,k))* &
                                        (w(i,j,k-1) + w(i,j,k))*psif*dx*dy*dz(k)
#endif
          end do
        end do
      end do
      call MPI_ALLREDUCE(MPI_IN_PLACE,um(1),ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,vm(1),ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,wm(1),ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,pm(1),ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,u2(1),ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,v2(1),ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,w2(1),ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,p2(1),ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,uw(1),ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
#if defined(_IBM)
      call MPI_ALLREDUCE(MPI_IN_PLACE,uf(1),ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,vf(1),ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,wf(1),ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,pf(1),ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,uf2(1),ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,vf2(1),ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,wf2(1),ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,pf2(1),ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,ufwf(1),ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,mean_psif(1),ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
#endif
      um(:) = um(:)*grid_area_ratio
      vm(:) = vm(:)*grid_area_ratio
      wm(:) = wm(:)*grid_area_ratio
      pm(:) = pm(:)*grid_area_ratio
      u2(:) = u2(:)*grid_area_ratio - um(:)**2
      v2(:) = v2(:)*grid_area_ratio - vm(:)**2
      w2(:) = w2(:)*grid_area_ratio - wm(:)**2
      p2(:) = p2(:)*grid_area_ratio - pm(:)**2
      uw(:) = uw(:)*grid_area_ratio - um(:)*wm(:)
#if defined(_IBM)
      do k=1,ng(3)
       if (mean_psif(k)/=0.0_rp) then
        uf(k)   = uf(k)/mean_psif(k)
        vf(k)   = vf(k)/mean_psif(k)
        wf(k)   = wf(k)/mean_psif(k)
        pf(k)   = pf(k)/mean_psif(k)
        uf2(k)  = uf2(k)/mean_psif(k) - uf(k)**2
        vf2(k)  = vf2(k)/mean_psif(k) - vf(k)**2
        wf2(k)  = wf2(k)/mean_psif(k) - wf(k)**2
        pf2(k)  = pf2(k)/mean_psif(k) - pf(k)**2
        ufwf(k) = ufwf(k)/mean_psif(k) - uf(k)*wf(k)
       end if
      end do
#endif
#if defined(_IBM)
      ! Dispersive velocities and stresses
      ud(:)   = 0._rp
      vd(:)   = 0._rp
      wd(:)   = 0._rp
      pd(:)   = 0._rp
      ud2(:)  = 0._rp
      vd2(:)  = 0._rp
      wd2(:)  = 0._rp
      pd2(:)  = 0._rp
      udwd(:) = 0._rp
      do k=lo(3),hi(3)
        do j=lo(2),hi(2)
          do i=lo(1),hi(1)
            ud(k)   = ud(k)   + (u(i,j,k) - um(k))
            vd(k)   = vd(k)   + (v(i,j,k) - vm(k))
            wd(k)   = wd(k)   + (0.5_rp*(w(i,j,k-1) + w(i,j,k)) - wm(k))
            pd(k)   = pd(k)   + (p(i,j,k) - pm(k))
            ud2(k)  = ud2(k)  + (u(i,j,k) - um(k))**2
            vd2(k)  = vd2(k)  + (v(i,j,k) - vm(k))**2
            wd2(k)  = wd2(k)  + 0.5_rp*((w(i,j,k) - wm(k))**2+(w(i,j,k-1) - wm(k))**2)
            pd2(k)  = pd2(k)  + (p(i,j,k) - pm(k))**2
            udwd(k) = udwd(k) + 0.25_rp*((u(i-1,j,k) - um(k)) + (u(i,j,k) - um(k)))* &
                                        ((w(i,j,k-1) - wm(k)) + (w(i,j,k) - wm(k)))
          end do
        end do
      end do
      call MPI_ALLREDUCE(MPI_IN_PLACE,ud(1)  ,ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,vd(1)  ,ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,wd(1)  ,ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,pd(1)  ,ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,ud2(1) ,ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,vd2(1) ,ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,wd2(1) ,ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,pd2(1) ,ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,udwd(1),ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      ud(:) = ud(:)*grid_area_ratio
      vd(:) = vd(:)*grid_area_ratio
      wd(:) = wd(:)*grid_area_ratio
      pd(:) = pd(:)*grid_area_ratio
      ud2(:) = ud2(:)*grid_area_ratio
      vd2(:) = vd2(:)*grid_area_ratio
      wd2(:) = wd2(:)*grid_area_ratio
      pd2(:) = pd2(:)*grid_area_ratio
      udwd(:) = udwd(:)*grid_area_ratio
#endif
      !
#if defined(_IBM)
      if(myid == 0) then
        open(newunit=iunit,file=fname)
        do k=1,ng(3)
          write(iunit,'(28E16.7e3)') z_g(k),um(k),vm(k),wm(k),pm(k), &
                                            u2(k),v2(k),w2(k),p2(k), &
                                            uw(k), &
                                            ud(k) ,vd(k) ,wd(k) ,pd(k) , &
                                            ud2(k),vd2(k),wd2(k),pd2(k), &
                                            udwd(k), &
                                            uf(k) ,vf(k) ,wf(k) ,pf(k) , &
                                            uf2(k),vf2(k),wf2(k),pf2(k), &
                                            ufwf(k)
        end do
        close(iunit)
      end if
#else
      if(myid == 0) then
        open(newunit=iunit,file=fname)
        do k=1,ng(3)
          write(iunit,'(10E16.7e3)') z_g(k),um(k),vm(k),wm(k),pm(k), &
                                            u2(k),v2(k),w2(k),p2(k), &
                                            uw(k)
        end do
        close(iunit)
      end if
#endif
    case(2)
    case(1)
    end select
  end subroutine out1d_chan
  !
  subroutine out1d_chan_tmp(fname,ng,lo,hi,idir,l,dl,dz,z_g,u,v,w,s,psi)
    implicit none
    character(len=*), intent(in) :: fname
    integer , intent(in), dimension(3)  :: ng,lo,hi
    integer , intent(in) :: idir
    real(rp), intent(in), dimension(3) :: l,dl
    real(rp), intent(in), dimension(0:) :: z_g
    real(rp), intent(in), dimension(lo(3)-1:) :: dz
    real(rp), intent(in), dimension(lo(1)-1:,lo(2)-1:,lo(3)-1:) :: u,v,w,s
    real(rp), intent(in), dimension(lo(1)-1:,lo(2)-1:,lo(3)-1:), optional :: psi
    real(rp), allocatable, dimension(:) :: um,vm,wm,sm,s2,us,vs,ws
#if defined(_IBM)
    real(rp), allocatable, dimension(:) :: ud,vd,wd,sd,sd2,udsd,vdsd,wdsd
    real(rp), allocatable, dimension(:) :: uf,vf,wf,sf,sf2,ufsf,vfsf,wfsf
    real(rp), allocatable, dimension(:) :: ssol,ssol2
    real(rp), allocatable, dimension(:) :: mean_psif,mean_psis
#endif
    integer :: i,j,k
    integer :: iunit
    integer :: q
    real(rp) :: grid_area_ratio
#if defined(_IBM)
    real(rp) :: dx,dy,psif,psis
#endif
    !
    q = ng(idir)
    grid_area_ratio = dl(1)*dl(2)/(l(1)*l(2))
    allocate(um(0:q+1),vm(0:q+1),wm(0:q+1),sm(0:q+1),s2(0:q+1),us(0:q+1),vs(0:q+1),ws(0:q+1))
#if defined(_IBM)
    dx = dl(1)
    dy = dl(2)
    allocate(ud(0:q+1),vd(0:q+1),wd(0:q+1),sd(0:q+1),sd2(0:q+1),udsd(0:q+1),vdsd(0:q+1),wdsd(0:q+1))
    allocate(uf(0:q+1),vf(0:q+1),wf(0:q+1),sf(0:q+1),sf2(0:q+1),ufsf(0:q+1),vfsf(0:q+1),wfsf(0:q+1))
    allocate(ssol(0:q+1),ssol2(0:q+1))
    allocate(mean_psif(0:q+1),mean_psis(0:q+1))
#endif
    um(:) = 0.0_rp
    vm(:) = 0.0_rp
    wm(:) = 0.0_rp
    sm(:) = 0.0_rp
    s2(:) = 0.0_rp
    us(:) = 0.0_rp
    vs(:) = 0.0_rp
    ws(:) = 0.0_rp
#if defined(_IBM)
    uf(:)    = 0.0_rp
    vf(:)    = 0.0_rp
    wf(:)    = 0.0_rp
    sf(:)    = 0.0_rp
    sf2(:)   = 0.0_rp
    ufsf(:)  = 0.0_rp
    vfsf(:)  = 0.0_rp
    wfsf(:)  = 0.0_rp
    ssol(:)  = 0.0_rp
    ssol2(:) = 0.0_rp
    mean_psif(:) = 0.0_rp
    mean_psis(:) = 0.0_rp
#endif
    do k=lo(3),hi(3)
      do j=lo(2),hi(2)
        do i=lo(1),hi(1)
#if defined(_IBM)
          psis = psi(i,j,k)
          psif = 1.0_rp - psi(i,j,k)
          mean_psis(k) = mean_psis(k) + psis*dx*dy*dz(k)
          mean_psif(k) = mean_psif(k) + psif*dx*dy*dz(k)
#endif
          sm(k) = sm(k) + s(i,j,k)
          s2(k) = s2(k) + s(i,j,k)**2
          um(k) = um(k) + u(i,j,k)
          vm(k) = vm(k) + v(i,j,k)
          wm(k) = wm(k) + 0.5_rp*(w(i,j,k-1) + w(i,j,k))
          us(k) = us(k) + 0.5_rp*(u(i-1,j,k) + u(i,j,k))*s(i,j,k)
          vs(k) = vs(k) + 0.5_rp*(v(i,j-1,k) + v(i,j,k))*s(i,j,k)
          ws(k) = ws(k) + 0.5_rp*(w(i,j,k-1) + w(i,j,k))*s(i,j,k)
#if defined(_IBM)
          sf(k)    = sf(k)    + s(i,j,k)*psif*dx*dy*dz(k)
          sf2(k)   = sf2(k)   + (s(i,j,k)**2)*psif*dx*dy*dz(k)
          uf(k)    = uf(k)    + u(i,j,k)*psif*dx*dy*dz(k)
          vf(k)    = vf(k)    + v(i,j,k)*psif*dx*dy*dz(k)
          wf(k)    = wf(k)    + 0.5_rp*(w(i,j,k-1) + w(i,j,k))*psif*dx*dy*dz(k)
          ufsf(k)  = ufsf(k)  + 0.5_rp*(u(i-1,j,k) + u(i,j,k))*s(i,j,k)*psif*dx*dy*dz(k)
          vfsf(k)  = vfsf(k)  + 0.5_rp*(v(i,j-1,k) + v(i,j,k))*s(i,j,k)*psif*dx*dy*dz(k)
          wfsf(k)  = wfsf(k)  + 0.5_rp*(w(i,j,k-1) + w(i,j,k))*s(i,j,k)*psif*dx*dy*dz(k)
          ssol(k)  = ssol(k)  + s(i,j,k)*psis*dx*dy*dz(k)
          ssol2(k) = ssol2(k) + (s(i,j,k)**2)*psis*dx*dy*dz(k)
#endif
        enddo
      enddo
    enddo
    call MPI_ALLREDUCE(MPI_IN_PLACE,sm(1),ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE,um(1),ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE,vm(1),ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE,wm(1),ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE,s2(1),ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE,us(1),ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE,vs(1),ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE,ws(1),ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
#if defined(_IBM)
    call MPI_ALLREDUCE(MPI_IN_PLACE,sf(1),ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE,uf(1),ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE,vf(1),ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE,wf(1),ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE,sf2(1),ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE,ufsf(1),ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE,vfsf(1),ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE,wfsf(1),ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE,ssol(1),ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE,ssol2(1),ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE,mean_psis(1),ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE,mean_psif(1),ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
#endif
    um(:) = um(:)*grid_area_ratio
    vm(:) = vm(:)*grid_area_ratio
    wm(:) = wm(:)*grid_area_ratio
    sm(:) = sm(:)*grid_area_ratio  
    s2(:) = s2(:)*grid_area_ratio - sm(:)**2
    us(:) = us(:)*grid_area_ratio - um(:)*sm(:)
    vs(:) = vs(:)*grid_area_ratio - vm(:)*sm(:)
    ws(:) = ws(:)*grid_area_ratio - wm(:)*sm(:)
#if defined(_IBM)
    do k=1,ng(3)
     if (mean_psif(k)/=0.0_rp) then
      uf(k)    = uf(k)/mean_psif(k)
      vf(k)    = vf(k)/mean_psif(k)
      wf(k)    = wf(k)/mean_psif(k)
      sf(k)    = sf(k)/mean_psif(k)  
      sf2(k)   = sf2(k)/mean_psif(k)  - sf(k)**2
      ufsf(k)  = ufsf(k)/mean_psif(k) - uf(k)*sf(k)
      vfsf(k)  = vfsf(k)/mean_psif(k) - vf(k)*sf(k)
      wfsf(k)  = wfsf(k)/mean_psif(k) - wf(k)*sf(k)
     end if
     if (mean_psis(k)/=0.0_rp) then
      ssol(k)  = ssol(k)/mean_psis(k)  
      ssol2(k) = ssol2(k)/mean_psis(k) - ssol(k)**2
     end if
    end do
#endif
#if defined(_IBM)
    ! Dispersive components
    sd(k)   = 0.0_rp
    sd2(k)  = 0.0_rp
    ud(k)   = 0.0_rp
    vd(k)   = 0.0_rp
    wd(k)   = 0.0_rp
    udsd(k) = 0.0_rp
    vdsd(k) = 0.0_rp
    wdsd(k) = 0.0_rp
    do k=lo(3),hi(3)
      do j=lo(2),hi(2)
        do i=lo(1),hi(1)
          sd(k)  = sd(k)  + (s(i,j,k) - sm(k))
          sd2(k) = sd2(k) + (s(i,j,k) - sm(k))**2
          ud(k)  = ud(k)  + (u(i,j,k) - um(k))
          vd(k)  = vd(k)  + (v(i,j,k) - vm(k))
          wd(k)  = wd(k)  + (0.5_rp*(w(i,j,k-1) + w(i,j,k)) - wm(k))
          udsd(k) = udsd(k) + 0.5_rp*((u(i-1,j,k) - um(k)) + (u(i,j,k) - um(k)))*(s(i,j,k) - sm(k))
          vdsd(k) = vdsd(k) + 0.5_rp*((v(i,j-1,k) - vm(k)) + (v(i,j,k) - vm(k)))*(s(i,j,k) - sm(k))
          wdsd(k) = wdsd(k) + 0.5_rp*((w(i,j,k-1) - wm(k)) + (w(i,j,k) - wm(k)))*(s(i,j,k) - sm(k))
        enddo
      enddo
    enddo
    call MPI_ALLREDUCE(MPI_IN_PLACE,sd(1),ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE,ud(1),ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE,vd(1),ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE,wd(1),ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE,sd2(1),ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE,udsd(1),ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE,vdsd(1),ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE,wdsd(1),ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    ud(:) = ud(:)*grid_area_ratio
    vd(:) = vd(:)*grid_area_ratio
    wd(:) = wd(:)*grid_area_ratio
    sd(:) = sd(:)*grid_area_ratio  
    sd2(:) = sd2(:)*grid_area_ratio
    udsd(:) = udsd(:)*grid_area_ratio
    vdsd(:) = vdsd(:)*grid_area_ratio
    wdsd(:) = wdsd(:)*grid_area_ratio
#endif
#if defined(_IBM)
    if(myid == 0) then
       open(newunit=iunit,file=fname)
       do k=1,ng(3)
         write(iunit,'(18E16.7e3)') z_g(k),sm(k),s2(k), &
                                           us(k),vs(k),ws(k), &
                                           sd(k),sd2(k), &
                                           udsd(k),vdsd(k),wdsd(k), &
                                           sf(k),sf2(k), &
                                           ufsf(k),vfsf(k),wfsf(k), &
                                           ssol(k),ssol2(k)
       enddo
      close(iunit)
    endif
#else
    if(myid == 0) then
       open(newunit=iunit,file=fname)
       do k=1,ng(3)
         write(iunit,'(6E16.7e3)') z_g(k),sm(k),s2(k), &
                                          us(k),vs(k),ws(k)
       enddo
      close(iunit)
    endif
#endif
  end subroutine out1d_chan_tmp
  !
  subroutine out2d_duct(fname,ng,lo,hi,idir,l,dl,z_g,u,v,w) ! e.g. for a duct
    !
    implicit none
    character(len=*), intent(in) :: fname
    integer , intent(in), dimension(3) :: ng,lo,hi
    integer , intent(in) :: idir
    real(rp), intent(in), dimension(3) :: l,dl
    real(rp), intent(in), dimension(0:) :: z_g
    real(rp), intent(in), dimension(lo(1)-1:,lo(2)-1:,lo(3)-1:) :: u,v,w
    real(rp), allocatable, dimension(:,:) :: um,vm,wm,u2,v2,w2,uv,uw,vw
    integer :: i,j,k
    integer :: iunit
    integer :: p,q
    real(rp) :: x_g,y_g,grid_area_ratio
    !
    select case(idir) ! streamwise direction
    case(3)
    case(2)
      grid_area_ratio = dl(2)/l(2)
      p = ng(1)
      q = ng(3)
      allocate(um(p,q),vm(p,q),wm(p,q),u2(p,q),v2(p,q),w2(p,q),uv(p,q),vw(p,q))
      !
      um(:,:) = 0.0_rp
      vm(:,:) = 0.0_rp
      wm(:,:) = 0.0_rp
      u2(:,:) = 0.0_rp
      v2(:,:) = 0.0_rp
      w2(:,:) = 0.0_rp
      uv(:,:) = 0.0_rp
      vw(:,:) = 0.0_rp
      do k=lo(3),hi(3)
        do i=lo(1),hi(1)
          do j=lo(2),hi(2)
            um(i,k) = um(i,k) + 0.5_rp*(u(i-1,j,k)+u(i,j,k))
            vm(i,k) = vm(i,k) + v(i,j,k)
            wm(i,k) = wm(i,k) + 0.5_rp*(w(i,j,k-1)+w(i,j,k))
            u2(i,k) = u2(i,k) + 0.5_rp*(u(i-1,j,k)**2+u(i,j,k)**2)
            v2(i,k) = v2(i,k) + v(i,j,k)**2
            w2(i,k) = w2(i,k) + 0.5_rp*(w(i,j,k-1)**2+w(i,j,k)**2)
            vw(i,k) = vw(i,k) + 0.25_rp*(v(i,j-1,k) + v(i,j,k))* &
                                     (w(i,j,k-1) + w(i,j,k))
            uv(i,k) = uv(i,k) + 0.25_rp*(u(i-1,j,k) + u(i,j,k))* &
                                     (v(i,j-1,k) + v(i,j,k))
          end do
        end do
      end do
      call MPI_ALLREDUCE(MPI_IN_PLACE,um(1,1),ng(1)*ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,vm(1,1),ng(1)*ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,wm(1,1),ng(1)*ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,u2(1,1),ng(1)*ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,v2(1,1),ng(1)*ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,w2(1,1),ng(1)*ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,vw(1,1),ng(1)*ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,uv(1,1),ng(1)*ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      um(:,:) =      um(:,:)*grid_area_ratio
      vm(:,:) =      vm(:,:)*grid_area_ratio
      wm(:,:) =      wm(:,:)*grid_area_ratio
      u2(:,:) =      u2(:,:)*grid_area_ratio - um(:,:)**2
      v2(:,:) =      v2(:,:)*grid_area_ratio - vm(:,:)**2
      w2(:,:) =      w2(:,:)*grid_area_ratio - wm(:,:)**2
      vw(:,:) =      vw(:,:)*grid_area_ratio - vm(:,:)*wm(:,:)
      uv(:,:) =      uv(:,:)*grid_area_ratio - um(:,:)*vm(:,:)
      if(myid == 0) then
        open(newunit=iunit,file=fname)
        do k=1,ng(3)
          do i=1,ng(1)
            x_g = (i-.5)*dl(1)
            write(iunit,'(10E16.7e3)') x_g,z_g(k),um(i,k),vm(i,k),wm(i,k), &
                                                  u2(i,k),v2(i,k),w2(i,k), &
                                                  vw(i,k),uv(i,k)
          end do
        end do
        close(iunit)
      end if
    case(1)
      grid_area_ratio = dl(1)/l(1)
      p = ng(2)
      q = ng(3)
      allocate(um(p,q),vm(p,q),wm(p,q),u2(p,q),v2(p,q),w2(p,q),uv(p,q),uw(p,q))
      !
      um(:,:) = 0.0_rp
      vm(:,:) = 0.0_rp
      wm(:,:) = 0.0_rp
      u2(:,:) = 0.0_rp
      v2(:,:) = 0.0_rp
      w2(:,:) = 0.0_rp
      uv(:,:) = 0.0_rp
      uw(:,:) = 0.0_rp
      do k=lo(3),hi(3)
        do j=lo(2),hi(2)
          do i=lo(1),hi(1)
            um(j,k) = um(j,k) + u(i,j,k)
            vm(j,k) = vm(j,k) + 0.5_rp*(v(i,j-1,k)+v(i,j,k))
            wm(j,k) = wm(j,k) + 0.5_rp*(w(i,j,k-1)+w(i,j,k))
            u2(j,k) = u2(j,k) + u(i,j,k)**2
            v2(j,k) = v2(j,k) + 0.5_rp*(v(i,j-1,k)**2+v(i,j,k)**2)
            w2(j,k) = w2(j,k) + 0.5_rp*(w(i,j,k-1)**2+w(i,j,k)**2)
            uv(j,k) = uv(j,k) + 0.25_rp*(u(i-1,j,k) + u(i,j,k))* &
                                     (v(i,j-1,k) + v(i,j,k))
            uw(j,k) = uw(j,k) + 0.25_rp*(u(i-1,j,k) + u(i,j,k))* &
                                     (w(i,j,k-1) + w(i,j,k))
          end do
        end do
      end do
      call MPI_ALLREDUCE(MPI_IN_PLACE,um(1,1),ng(1)*ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,vm(1,1),ng(1)*ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,wm(1,1),ng(1)*ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,u2(1,1),ng(1)*ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,v2(1,1),ng(1)*ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,w2(1,1),ng(1)*ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,uv(1,1),ng(1)*ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,uw(1,1),ng(1)*ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      um(:,:) =      um(:,:)*grid_area_ratio
      vm(:,:) =      vm(:,:)*grid_area_ratio
      wm(:,:) =      wm(:,:)*grid_area_ratio
      u2(:,:) = u2(:,:)*grid_area_ratio - um(:,:)**2
      v2(:,:) = v2(:,:)*grid_area_ratio - vm(:,:)**2
      w2(:,:) = w2(:,:)*grid_area_ratio - wm(:,:)**2
      uv(:,:) =      uv(:,:)*grid_area_ratio - um(:,:)*vm(:,:)
      uw(:,:) =      uw(:,:)*grid_area_ratio - um(:,:)*wm(:,:)
      if(myid == 0) then
        open(newunit=iunit,file=fname)
        do k=1,ng(3)
          do j=1,ng(2)
            y_g = (j-.5)*dl(2)
            write(iunit,'(10E16.7e3)') y_g,z_g(k),um(j,k),vm(j,k),wm(j,k), &
                                                  u2(j,k),v2(j,k),w2(j,k), &
                                                  uv(j,k),uw(j,k)
          end do
        end do
        close(iunit)
      end if
    end select
  end subroutine out2d_duct
!
end module mod_output
