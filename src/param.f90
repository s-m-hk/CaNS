! -
!
! SPDX-FileCopyrightText: Copyright (c) 2017-2022 Pedro Costa and the CaNS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
module mod_param
use mod_types
!@acc use cudecomp
implicit none
public
!
! parameters
!
real(rp), parameter :: pi = 4._rp*atan(1._rp)
#if !defined(_EPS_EXACT_ZERO) /* recommended */
real(sp), parameter :: eps_sp = epsilon(1._sp)
real(dp), parameter :: eps_dp = epsilon(1._dp)
real(rp), parameter :: eps_rp = epsilon(1._rp)
#else
real(sp), parameter :: eps_sp = 0._sp
real(dp), parameter :: eps_dp = 0._dp
real(rp), parameter :: eps_rp = 0._rp
#endif
real(gp), parameter :: small = epsilon(1._gp)*10**(precision(1._gp)/2)
character(len=100), parameter :: datadir = 'data/'
real(rp), parameter, dimension(2,3) :: rkcoeff = reshape([32._rp/60._rp,  0._rp        , &
                                                          25._rp/60._rp, -17._rp/60._rp, &
                                                          45._rp/60._rp, -25._rp/60._rp], shape(rkcoeff))
real(rp), parameter, dimension(3)   :: rkcoeff12 = rkcoeff(1,:)+rkcoeff(2,:)
!
! variables to be determined from the input file 'dns.in'
!
integer , protected :: itot,jtot,ktot
real(rp), protected :: lx,ly,dx,dy,dz,dxi,dyi,dzi,gr
real(rp) :: lz
real(rp), protected :: cfl,dtmin
real(rp), protected :: uref,lref,rey,visc
!
character(len=100), protected :: inivel
logical, protected :: is_wallturb
!
integer , protected :: nstep
real(rp), protected :: time_max,tw_max
logical , protected, dimension(3) :: stop_type
logical , protected :: restart,is_overwrite_save,reset_time
integer , protected :: nsaves_max
integer , protected :: ioutput,icheck,iout0d,iout1d,iout2d,iout3d,ioutLPP,isave
!
integer , dimension(2) :: dims
!
integer, dimension(0:1,3) :: nb
logical, dimension(0:1,3) :: is_bound
character(len=1), protected, dimension(0:1,3,3) ::  cbcvel
real(rp)        , protected, dimension(0:1,3,3) ::  bcvel
character(len=1), protected, dimension(0:1,3)   ::  cbcpre
real(rp)        , protected, dimension(0:1,3)   ::  bcpre
!
real(rp), protected, dimension(3) :: gacc 
real(rp), protected, dimension(3) :: bforce
logical , protected, dimension(4) :: is_forced
real(rp), protected, dimension(3) :: velf
!
integer , protected, dimension(3) :: ng
real(rp), dimension(3) :: l
real(rp), protected, dimension(3) :: dl
real(rp), protected, dimension(3) :: dli
!
character(len=3)   :: itmp
real(rp)           :: alph_f,alph_s,Pr,dr,tg0
real(rp)           :: solidtemp
real(rp)           :: tmp0,beta_th
character(len=1), dimension(0:1,3) :: cbctmp
real(rp)        , dimension(0:1,3) ::  bctmp
real(rp) :: solid_height_ratio
real(rp) :: Rotation_angle
real(rp) :: sx, sy, sz, depth, rod
character(len=15) :: surface_type
logical :: height_map
!
#if defined(_OPENACC)
!
! cuDecomp input parameters
!
integer, protected :: cudecomp_t_comm_backend,cudecomp_h_comm_backend
logical, protected :: cudecomp_is_t_comm_autotune ,cudecomp_is_h_comm_autotune , &
                      cudecomp_is_t_enable_nccl   ,cudecomp_is_h_enable_nccl   , &
                      cudecomp_is_t_enable_nvshmem,cudecomp_is_h_enable_nvshmem, &
                      cudecomp_is_t_in_place
logical :: exists
#endif
contains
  subroutine read_input(myid)
  use mpi
  implicit none
  integer, intent(in) :: myid
  integer :: iunit,ierr
    open(newunit=iunit,file='dns.in',status='old',action='read',iostat=ierr)
      if( ierr == 0 ) then
        read(iunit,*,iostat=ierr) itot,jtot,ktot
        read(iunit,*,iostat=ierr) lx,ly,lz
        read(iunit,*,iostat=ierr) gr
        read(iunit,*,iostat=ierr) cfl,dtmin
        read(iunit,*,iostat=ierr) uref,lref,rey
        read(iunit,*,iostat=ierr) inivel
        read(iunit,*,iostat=ierr) is_wallturb
        read(iunit,*,iostat=ierr) nstep, time_max,tw_max
        read(iunit,*,iostat=ierr) stop_type(1),stop_type(2),stop_type(3)
        read(iunit,*,iostat=ierr) restart,is_overwrite_save,nsaves_max,reset_time; if( ierr /= 0 ) nsaves_max = 0 ! a good default, for backward compatibility
        read(iunit,*,iostat=ierr) ioutput,icheck,iout0d,iout1d,iout2d,iout3d,ioutLPP,isave
        read(iunit,*,iostat=ierr) cbcvel(0,1,1),cbcvel(1,1,1),cbcvel(0,2,1),cbcvel(1,2,1),cbcvel(0,3,1),cbcvel(1,3,1)
        read(iunit,*,iostat=ierr) cbcvel(0,1,2),cbcvel(1,1,2),cbcvel(0,2,2),cbcvel(1,2,2),cbcvel(0,3,2),cbcvel(1,3,2)
        read(iunit,*,iostat=ierr) cbcvel(0,1,3),cbcvel(1,1,3),cbcvel(0,2,3),cbcvel(1,2,3),cbcvel(0,3,3),cbcvel(1,3,3)
        read(iunit,*,iostat=ierr) cbcpre(0,1  ),cbcpre(1,1  ),cbcpre(0,2  ),cbcpre(1,2  ),cbcpre(0,3  ),cbcpre(1,3  )
        read(iunit,*,iostat=ierr) bcvel(0,1,1), bcvel(1,1,1), bcvel(0,2,1), bcvel(1,2,1), bcvel(0,3,1), bcvel(1,3,1)
        read(iunit,*,iostat=ierr) bcvel(0,1,2), bcvel(1,1,2), bcvel(0,2,2), bcvel(1,2,2), bcvel(0,3,2), bcvel(1,3,2)
        read(iunit,*,iostat=ierr) bcvel(0,1,3), bcvel(1,1,3), bcvel(0,2,3), bcvel(1,2,3), bcvel(0,3,3), bcvel(1,3,3)
        read(iunit,*,iostat=ierr) bcpre(0,1  ), bcpre(1,1  ), bcpre(0,2  ), bcpre(1,2  ), bcpre(0,3  ), bcpre(1,3  )
        read(iunit,*,iostat=ierr) gacc(1), gacc(2), gacc(3)
        read(iunit,*,iostat=ierr) bforce(1),bforce(2),bforce(3)
        read(iunit,*,iostat=ierr) is_forced(1),is_forced(2),is_forced(3)
        read(iunit,*,iostat=ierr) velf(1),velf(2),velf(3)
        read(iunit,*,iostat=ierr) dims(1),dims(2)
      else
        if(myid == 0) print*, 'Error reading main input file'
        if(myid == 0) print*, 'Aborting...'
        call MPI_FINALIZE(ierr)
        error stop
      end if
    close(iunit)
    !  
    ! load immersed boundary method input files
    !
#if defined(_IBM)
    open(newunit=iunit,file='ibm.in',status='old',action='read',iostat=ierr)
      if( ierr == 0 ) then
        read(iunit,*,iostat=ierr) surface_type
        read(iunit,*,iostat=ierr) solid_height_ratio
        read(iunit,*,iostat=ierr) Rotation_angle
        read(iunit,*,iostat=ierr) sx, sy, sz
        read(iunit,*,iostat=ierr) rod
        read(iunit,*,iostat=ierr) depth
      else
        if(myid == 0) print*, 'Error reading IBM input file' 
        if(myid == 0) print*, 'Aborting...'
        call MPI_FINALIZE(ierr)
        error stop
      endif
    close(iunit)
    if (trim(surface_type) == 'HeightMap') height_map = .true.
#endif
    !  
    ! load heat transfer input files
    !
#if defined(_HEAT_TRANSFER)
    open(newunit=iunit,file='heat.in',status='old',action='read',iostat=ierr)
    if( ierr.eq.0 ) then
      read(iunit,*,iostat=ierr) itmp
      read(iunit,*,iostat=ierr) tg0
      read(iunit,*,iostat=ierr) Pr
      read(iunit,*,iostat=ierr) dr
      read(iunit,*,iostat=ierr) cbctmp(0,1 ),cbctmp(1,1 ),cbctmp(0,2 ),cbctmp(1,2 ),cbctmp(0,3 ),cbctmp(1,3 )
      read(iunit,*,iostat=ierr) bctmp(0,1  ),bctmp(1,1  ),bctmp(0,2  ),bctmp(1,2  ),bctmp(0,3  ),bctmp(1,3  )
      read(iunit,*,iostat=ierr) tmp0,beta_th
      read(iunit,*,iostat=ierr) solidtemp
      read(iunit,*,iostat=ierr) is_forced(4)
      read(iunit,*,iostat=ierr) velf(4)
    else
      if(myid == 0) print*, 'Error reading heat input file' 
      if(myid == 0) print*, 'Aborting...'
      call MPI_FINALIZE(ierr)
      call exit
    endif
    close(iunit)
#endif
    dx = lx/(1.*itot)
    dy = ly/(1.*jtot)
    dz = lz/(1.*ktot)
    dxi = dx**(-1)
    dyi = dy**(-1)
    dzi = dz**(-1)
    !
    visc = uref*lref/rey
#if defined(_HEAT_TRANSFER)
    alph_f = visc/Pr
    alph_s = dr*alph_f
#endif
    ng  = [itot,jtot,ktot]
    l   = [lx,ly,lz]
    dl  = [dx,dy,dz]
    dli = [dxi,dyi,dzi]
#if defined(_OPENACC)
    !
    ! read cuDecomp parameter file cudecomp.in, if it exists
    !
    ! defaults
    !
    cudecomp_is_t_comm_autotune  = .true.
    cudecomp_is_h_comm_autotune  = .true.
    cudecomp_is_t_enable_nccl    = .true.
    cudecomp_is_h_enable_nccl    = .true.
    cudecomp_is_t_enable_nvshmem = .true.
    cudecomp_is_h_enable_nvshmem = .true.
    inquire(file='cudecomp.in', exist = exists)
    if(exists) then
      open(newunit=iunit,file='cudecomp.in',status='old',action='read',iostat=ierr)
        if( ierr == 0 ) then
          read(iunit,*) cudecomp_t_comm_backend,cudecomp_is_t_enable_nccl,cudecomp_is_t_enable_nvshmem
          read(iunit,*) cudecomp_h_comm_backend,cudecomp_is_h_enable_nccl,cudecomp_is_h_enable_nvshmem
          if(cudecomp_t_comm_backend >= 1 .and. cudecomp_t_comm_backend <= 7) then
            cudecomp_is_t_comm_autotune = .false. ! do not autotune if backend is prescribed
            select case(cudecomp_t_comm_backend)
            case(1)
              cudecomp_t_comm_backend = CUDECOMP_TRANSPOSE_COMM_MPI_P2P
            case(2)
              cudecomp_t_comm_backend = CUDECOMP_TRANSPOSE_COMM_MPI_P2P_PL
            case(3)
              cudecomp_t_comm_backend = CUDECOMP_TRANSPOSE_COMM_MPI_A2A
            case(4)
              cudecomp_t_comm_backend = CUDECOMP_TRANSPOSE_COMM_NCCL
            case(5)
              cudecomp_t_comm_backend = CUDECOMP_TRANSPOSE_COMM_NCCL_PL
            case(6)
              cudecomp_t_comm_backend = CUDECOMP_TRANSPOSE_COMM_NVSHMEM
            case(7)
              cudecomp_t_comm_backend = CUDECOMP_TRANSPOSE_COMM_NVSHMEM_PL
            case default
              cudecomp_t_comm_backend = CUDECOMP_TRANSPOSE_COMM_MPI_P2P
            end select
          end if
          if(cudecomp_h_comm_backend >= 1 .and. cudecomp_h_comm_backend <= 4) then
            cudecomp_is_h_comm_autotune = .false. ! do not autotune if backend is prescribed
            select case(cudecomp_h_comm_backend)
            case(1)
              cudecomp_h_comm_backend = CUDECOMP_HALO_COMM_MPI
            case(2)
              cudecomp_h_comm_backend = CUDECOMP_HALO_COMM_MPI_BLOCKING
            case(3)
              cudecomp_h_comm_backend = CUDECOMP_HALO_COMM_NCCL
            case(4)
              cudecomp_h_comm_backend = CUDECOMP_HALO_COMM_NVSHMEM
            case(5)
              cudecomp_h_comm_backend = CUDECOMP_HALO_COMM_NVSHMEM_BLOCKING
            case default
              cudecomp_h_comm_backend = CUDECOMP_HALO_COMM_MPI
            end select
          end if
        else
          if(myid == 0) print*, 'Error reading cuDecomp input file'
          if(myid == 0) print*, 'Aborting...'
          call MPI_FINALIZE(ierr)
          error stop
        end if
      close(iunit)
    end if
    !
    ! manually set cuDecomp out-of-place transposes by default
    !
    cudecomp_is_t_in_place = .false.
#endif
  end subroutine read_input
end module mod_param
