module mod_initIBM
#if defined(_IBM)
use mpi
#if !defined(_IBM_BC)
use mod_IBM, only: IBM_Mask
#else
use mod_IBM, only: IBM_Mask, normal_vectors,intersect,mirrorpoints, &
                   interpolation_2D_velocity,interpolation_dphi, &
                   mirrorpoints_ijk, interpolation_mirror,InterpolationWeights
#endif
use mod_bound, only: bounduvw,boundp
use mod_load,  only: load_all,load_one
use mod_param, only: datadir
use mod_common_mpi, only: myid
use mod_types
implicit none
private
public initIBM
!
contains
#if defined(_SIMPLE)
subroutine initIBM(cbcvel,cbcpre,bcvel,bcpre,is_bound,n,nh,halo,ng,nb,lo,hi,cell_u_tag,cell_v_tag,cell_w_tag,cell_phi_tag, &
                   ldz,zc,zf,zf_g,dzc,dzf,dl,dli)
implicit none
character(len=1), intent(in), dimension(0:1,3,3)            :: cbcvel
real(rp), intent(in), dimension(0:1,3,3)                    :: bcvel
character(len=1), intent(in), dimension(0:1,3)              :: cbcpre
real(rp), intent(in), dimension(0:1,3)                      :: bcpre
integer , intent(in), dimension(0:1,3  )                    :: nb
integer , intent(in)                                        :: nh
integer , intent(in), dimension(3)                          :: halo
logical , intent(in), dimension(0:1,3  )                    :: is_bound
integer , intent(in), dimension(3)                          :: n,ng,lo,hi
integer , intent(in )                                       :: ldz
real(rp), intent(in ), dimension(ldz:)                      :: zc,zf,zf_g,dzc,dzf,dl,dli
real(rp), intent(out),dimension(0:,0:,0:) :: cell_u_tag,cell_v_tag,cell_w_tag,cell_phi_tag
real(rp) :: dummy_time
integer  :: dummy_istep
integer  :: i,j,k,nx,ny,nz,idir
logical  :: is_data
!
nx = n(1)
ny = n(2)
nz = n(3)
!
cell_u_tag(1:nx,1:ny,1:nz)    = 0.0_rp
cell_v_tag(1:nx,1:ny,1:nz)    = 0.0_rp
cell_w_tag(1:nx,1:ny,1:nz)    = 0.0_rp
cell_phi_tag(1:nx,1:ny,1:nz)  = 0.0_rp
!
inquire(file=trim(datadir)//'IBM.bin',exist=is_data)

if (.not.is_data) then
   call IBM_mask(n,ng,lo,hi,zc,zf,zf_g,dzc,dzf,cell_phi_tag)
   !$acc enter data copyin(cell_phi_tag)
   call boundp(cbcpre,n,nh,halo,bcpre,nb,is_bound,dl,dzc,cell_phi_tag)
   !
   do k=1,nz
    do j=1,ny
     do i=1,nx
     if((cell_phi_tag(i,j,k) + cell_phi_tag(i+1,j,k)) > 0.5_rp) cell_u_tag(i,j,k) = 1.0_rp
     if((cell_phi_tag(i,j,k) + cell_phi_tag(i,j+1,k)) > 0.5_rp) cell_v_tag(i,j,k) = 1.0_rp
     if((cell_phi_tag(i,j,k) + cell_phi_tag(i,j,k+1)) > 0.5_rp) cell_w_tag(i,j,k) = 1.0_rp
     enddo
    enddo
   enddo
   !$acc enter data copyin(cell_u_tag,cell_v_tag,cell_w_tag)
   call boundp(cbcpre,n,nh,halo,bcpre,nb,is_bound,dl,dzc,cell_u_tag)
   call boundp(cbcpre,n,nh,halo,bcpre,nb,is_bound,dl,dzc,cell_v_tag)
   call boundp(cbcpre,n,nh,halo,bcpre,nb,is_bound,dl,dzc,cell_w_tag)
   !
   dummy_time = 0.0_rp; dummy_istep = 0
   !$acc update self(cell_u_tag,cell_v_tag,cell_w_tag,cell_phi_tag)
   call load_all('w',trim(datadir)//'IBM.bin',MPI_COMM_WORLD,ng,[1,1,1],lo,hi,dummy_time,dummy_istep,cell_u_tag,cell_v_tag,cell_w_tag,cell_phi_tag)
else
   call load_all('r',trim(datadir)//'IBM.bin',MPI_COMM_WORLD,ng,[1,1,1],lo,hi,dummy_time,dummy_istep,cell_u_tag,cell_v_tag,cell_w_tag,cell_phi_tag)
   !
   !$acc enter data copyin(cell_phi_tag,cell_u_tag,cell_v_tag,cell_w_tag)
    call boundp(cbcpre,n,nh,halo,bcpre,nb,is_bound,dl,dzc,cell_phi_tag)
    call bounduvw(cbcvel,n,nh,halo,bcvel,nb,is_bound,.false.,dl,dzc,dzf,cell_u_tag,cell_v_tag,cell_w_tag)
  if (myid == 0)  print*, '*** Saved IBM data loaded ***'
   !---------------------------------------------------------------------
endif

end subroutine initIBM
#endif
!
#if defined(_VOLUME)
subroutine initIBM(cbcvel,cbcpre,bcvel,bcpre,is_bound,n,nh,halo,ng,nb,lo,hi,cell_u_tag,cell_v_tag,cell_w_tag,cell_phi_tag, &
                   Level_set, &
                   surf_height, &
                   ldz,zc,zf,zf_g,dzc,dzf,dl,dli &
#if defined(_IBM_BC)
                   ,nx_surf,ny_surf,nz_surf,nabs_surf,i_mirror,j_mirror,k_mirror, &
                   i_IP1,j_IP1,k_IP1,i_IP2,j_IP2,k_IP2,WP1,WP2,deltan &
#endif
                   )
implicit none
character(len=1), intent(in), dimension(0:1,3,3)          :: cbcvel
real(rp), intent(in), dimension(0:1,3,3)                  :: bcvel
character(len=1), intent(in), dimension(0:1,3)            :: cbcpre
real(rp), intent(in), dimension(0:1,3)                    :: bcpre
integer , intent(in), dimension(0:1,3  )                  :: nb
integer , intent(in)                                      :: nh
integer , intent(in), dimension(3)                        :: halo
logical , intent(in), dimension(0:1,3  )                  :: is_bound
integer , intent(in), dimension(3)                        :: n,ng,lo,hi
integer , intent(in )                                     :: ldz
real(rp), intent(in ), dimension(ldz:)                    :: zc,zf,zf_g,dzc,dzf,dl,dli
real(rp), intent(in ), dimension(1:n(1),1:n(2))           :: surf_height
#if defined(_IBM_BC)
real(rp), intent(out),dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6) :: cell_u_tag,cell_v_tag,cell_w_tag,cell_phi_tag
integer,  intent(out),dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6) :: Level_set
#else
real(rp), intent(out),dimension(0:n(1)+1,0:n(2)+1,0:n(3)+1) :: cell_u_tag,cell_v_tag,cell_w_tag,cell_phi_tag
integer,  intent(out),dimension(0:n(1)+1,0:n(2)+1,0:n(3)+1) :: Level_set
#endif
real(rp), dimension(:,:,:), allocatable :: tmp
#if defined(_IBM_BC)
integer,optional,intent(out),dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6) :: i_mirror,j_mirror,k_mirror, &
                                                                          i_IP1,j_IP1,k_IP1, &
                                                                          i_IP2,j_IP2,k_IP2
real(rp), dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6),     intent(out), optional :: nx_surf,ny_surf,nz_surf,nabs_surf,deltan
real(rp), dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6,1:7), intent(out), optional :: WP1,WP2
real(rp), dimension(:,:,:), allocatable :: z_intersect(:,:,:),y_intersect(:,:,:), x_intersect(:,:,:)
real(rp), dimension(:,:,:), allocatable :: x_mirror(:,:,:), y_mirror(:,:,:), z_mirror(:,:,:)
real(rp), dimension(:,:,:), allocatable :: x_IP1(:,:,:), y_IP1(:,:,:), z_IP1(:,:,:) !can be deallocated later
real(rp), dimension(:,:,:), allocatable :: x_IP2(:,:,:), y_IP2(:,:,:), z_IP2(:,:,:) !can be deallocated later
#endif
real(rp) :: dummy_time
integer  :: dummy_istep
integer  :: i,j,k,nx,ny,nz,h,idir
logical  :: is_data
!
nx = n(1)
ny = n(2)
nz = n(3)
!
cell_u_tag(1:nx,1:ny,1:nz)    = 0.0_rp
cell_v_tag(1:nx,1:ny,1:nz)    = 0.0_rp
cell_w_tag(1:nx,1:ny,1:nz)    = 0.0_rp
cell_phi_tag(1:nx,1:ny,1:nz)  = 0.0_rp
Level_set(1:nx,1:ny,1:nz)     = 0
!
inquire(file=trim(datadir)//'IBM.bin',exist=is_data)

if (.not.is_data) then
#if defined(_IBM_BC)
  allocate(z_intersect(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6),y_intersect(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6),&
           x_intersect(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6))
  allocate(x_mirror(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6), y_mirror(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6), &
           z_mirror(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6))
  allocate(x_IP1(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6), y_IP1(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6), &
           z_IP1(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6))
  allocate(x_IP2(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6), y_IP2(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6), &
           z_IP2(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6))
  allocate(tmp(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6))
#else
  allocate(tmp(0:n(1)+1,0:n(2)+1,0:n(3)+1))
#endif

  !$acc enter data copyin(cell_phi_tag,cell_u_tag,cell_v_tag,cell_w_tag,Level_set)
  call IBM_mask(n,ng,lo,hi,zc,zf,zf_g,dzc,dzf,is_bound,cell_phi_tag,Level_set,surf_height)

      !$acc enter data copyin(tmp)
       tmp(1:nx,1:ny,1:nz) = 1._rp*Level_set(1:nx,1:ny,1:nz)
       call boundp(cbcpre,n,nh,halo,bcpre,nb,is_bound,dl,dzc,tmp)
       Level_set(0:nx+1,0:ny+1,0:nz+1) = int(tmp(0:nx+1,0:ny+1,0:nz+1))
      !$acc exit data copyout(tmp) ! not needed on the device
   if (myid == 0)  print*, '*** Solid marker set! ***'
  !---------------------------------------------------------------------
       call boundp(cbcpre,n,nh,halo,bcpre,nb,is_bound,dl,dzc,cell_phi_tag)
       !$acc kernels default(present) async(1)
       do k=1,nz
        do j=1,ny
         do i=1,nx
          cell_u_tag(i,j,k) = 0.5_rp*(cell_phi_tag(i+1,j,k)+cell_phi_tag(i,j,k))
          cell_v_tag(i,j,k) = 0.5_rp*(cell_phi_tag(i,j+1,k)+cell_phi_tag(i,j,k))
          cell_w_tag(i,j,k) = 0.5_rp*(cell_phi_tag(i,j,k+1)+cell_phi_tag(i,j,k))
         enddo
        enddo
       enddo
       !$acc end kernels
       call boundp(cbcpre,n,nh,halo,bcpre,nb,is_bound,dl,dzc,cell_u_tag)
       call boundp(cbcpre,n,nh,halo,bcpre,nb,is_bound,dl,dzc,cell_v_tag)
       call boundp(cbcpre,n,nh,halo,bcpre,nb,is_bound,dl,dzf,cell_w_tag)
   if (myid == 0)  print*, '*** Volume fractions calculated! ***'
  !---------------------------------------------------------------------

#if defined(_IBM_BC)
  nx_surf(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6)      = 0.0_rp
  ny_surf(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6)      = 0.0_rp
  nz_surf(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6)      = 0.0_rp
  nabs_surf(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6)    = 0.0_rp
  nx_surf_nonnorm(-1:n(1)+1,-1:n(2)+1,-1:n(3)+1) = 0.0_rp
  ny_surf_nonnorm(-1:n(1)+1,-1:n(2)+1,-1:n(3)+1) = 0.0_rp
  nz_surf_nonnorm(-1:n(1)+1,-1:n(2)+1,-1:n(3)+1) = 0.0_rp

  call normal_vectors(n,lo,hi,Level_set,cell_phi_tag,nx_surf,ny_surf,nz_surf,nabs_surf,zc,zf,dzc,dzf,dl,dli,n,surf_height)
     !$acc enter data copyin(nx_surf,ny_surf,nz_surf,nabs_surf)
     call bounduvw(cbcvel,n,nh,halo,bcvel,nb,is_bound,.false.,dl,dzc,dzf,nx_surf,ny_surf,nz_surf)
     call boundp(cbcpre,n,nh,halo,bcpre,nb,is_bound,dl,dzc,nabs_surf)
  if (myid.eq.0)  print*, '*** Normal vectors calculated! ***'
  !*********************************************************************
  x_intersect(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6) = 0.0_rp
  y_intersect(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6) = 0.0_rp
  z_intersect(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6) = 0.0_rp

  call intersect(n,lo,hi,nx_surf,ny_surf,nz_surf,nabs_surf,x_intersect,y_intersect,z_intersect,zc,dzc,surf_height)
      !$acc enter data copyin(x_intersect,y_intersect,z_intersect)
      call boundp(cbcpre,n,nh,halo,bcpre,nb,is_bound,dl,dzc,x_intersect)
      call boundp(cbcpre,n,nh,halo,bcpre,nb,is_bound,dl,dzc,y_intersect)
      call boundp(cbcpre,n,nh,halo,bcpre,nb,is_bound,dl,dzc,z_intersect)
  if (myid.eq.0)  print*, '*** Intersect points have  been calculated! ***'
  !*********************************************************************

  x_mirror(-5:i1+5,-5:j1+5,-5:k1+5)    =   -1000.0_rp
  y_mirror(-5:i1+5,-5:j1+5,-5:k1+5)    =   -1000.0_rp
  z_mirror(-5:i1+5,-5:j1+5,-5:k1+5)    =   -1000.0_rp
  x_IP1(-5:i1+5,-5:j1+5,-5:k1+5)       =   -1000.0_rp
  x_IP2(-5:i1+5,-5:j1+5,-5:k1+5)       =   -1000.0_rp
  y_IP1(-5:i1+5,-5:j1+5,-5:k1+5)       =   -1000.0_rp
  y_IP2(-5:i1+5,-5:j1+5,-5:k1+5)       =   -1000.0_rp
  z_IP1(-5:i1+5,-5:j1+5,-5:k1+5)       =   -1000.0_rp
  z_IP2(-5:i1+5,-5:j1+5,-5:k1+5)       =   -1000.0_rp
  deltan(-5:i1+5,-5:j1+5,-5:k1+5)      =   -1000.0_rp

  call mirrorpoints(n,lo,hi, &
                    nx_surf,ny_surf,nz_surf,nabs_surf, &
                    x_intersect,y_intersect,z_intersect, &
                    x_mirror,y_mirror,z_mirror, &
                    deltan, &
                    x_IP1,y_IP1,z_IP1,x_IP2,y_IP2,z_IP2,& 
                    zc,dzc)

      !$acc enter data copyin(deltan,x_mirror,y_mirror,z_mirror,x_IP1,x_IP2,y_IP1,y_IP2,z_IP1,z_IP2)
      call boundp(cbcpre,n,nh,halo,bcpre,nb,is_bound,dl,dzc,deltan)
      call boundp(cbcpre,n,nh,halo,bcpre,nb,is_bound,dl,dzc,x_mirror)
      call boundp(cbcpre,n,nh,halo,bcpre,nb,is_bound,dl,dzc,y_mirror)
      call boundp(cbcpre,n,nh,halo,bcpre,nb,is_bound,dl,dzc,z_mirror)
      call boundp(cbcpre,n,nh,halo,bcpre,nb,is_bound,dl,dzc,x_IP1)
      call boundp(cbcpre,n,nh,halo,bcpre,nb,is_bound,dl,dzc,x_IP2)
      call boundp(cbcpre,n,nh,halo,bcpre,nb,is_bound,dl,dzc,y_IP1)
      call boundp(cbcpre,n,nh,halo,bcpre,nb,is_bound,dl,dzc,y_IP2)
      call boundp(cbcpre,n,nh,halo,bcpre,nb,is_bound,dl,dzc,z_IP1)
      call boundp(cbcpre,n,nh,halo,bcpre,nb,is_bound,dl,dzc,z_IP2)
  !*********************************************************************
      i_mirror(-5:i1+5,-5:j1+5,-5:k1+5)    =   -1000
      j_mirror(-5:i1+5,-5:j1+5,-5:k1+5)    =   -1000
      k_mirror(-5:i1+5,-5:j1+5,-5:k1+5)    =   -1000
      i_IP1(-5:i1+5,-5:j1+5,-5:k1+5)       =   -1000
      j_IP1(-5:i1+5,-5:j1+5,-5:k1+5)       =   -1000
      k_IP1(-5:i1+5,-5:j1+5,-5:k1+5)       =   -1000
      i_IP2(-5:i1+5,-5:j1+5,-5:k1+5)       =   -1000
      j_IP2(-5:i1+5,-5:j1+5,-5:k1+5)       =   -1000
      k_IP2(-5:i1+5,-5:j1+5,-5:k1+5)       =   -1000

  call mirrorpoints_ijk(n,lo,hi, &
                        nabs_surf,x_mirror,y_mirror,z_mirror, &
                        x_IP1,y_IP1,z_IP1,x_IP2,y_IP2,z_IP2, &
                        i_mirror,j_mirror,k_mirror, &
                        i_IP1,j_IP1,k_IP1,i_IP2,j_IP2,k_IP2, &
                        zf,dzc)
   do h = 1,nh
      i_mirror(:,:,lo(3)-h)    = i_mirror(:,:,lo(3))
      i_mirror(:,:,hi(3)+h)    = i_mirror(:,:,hi(3))
      j_mirror(:,:,lo(3)-h)    = j_mirror(:,:,lo(3))
      j_mirror(:,:,hi(3)+h)    = j_mirror(:,:,hi(3))
      k_mirror(:,:,lo(3)-h)    = k_mirror(:,:,lo(3))
      k_mirror(:,:,hi(3)+h)    = k_mirror(:,:,hi(3))
      !
      i_IP1(:,:,lo(3)-h)       = i_IP1(:,:,lo(3))
      i_IP1(:,:,hi(3)+h)       = i_IP1(:,:,hi(3))
      j_IP1(:,:,lo(3)-h)       = j_IP1(:,:,lo(3))
      j_IP1(:,:,hi(3)+h)       = j_IP1(:,:,hi(3))
      k_IP1(:,:,lo(3)-h)       = k_IP1(:,:,lo(3))
      k_IP1(:,:,hi(3)+h)       = k_IP1(:,:,hi(3))
      !
      i_IP2(:,:,lo(3)-h)       = i_IP2(:,:,lo(3))
      i_IP2(:,:,hi(3)+h)       = i_IP2(:,:,hi(3))
      j_IP2(:,:,lo(3)-h)       = j_IP2(:,:,lo(3))
      j_IP2(:,:,hi(3)+h)       = j_IP2(:,:,hi(3))
      k_IP2(:,:,lo(3)-h)       = k_IP2(:,:,lo(3))
      k_IP2(:,:,hi(3)+h)       = k_IP2(:,:,hi(3))
   end do
  
deallocate(z_intersect,y_intersect,x_intersect)
! deallocate(i_mirror,j_mirror,k_mirror)
deallocate(x_mirror,y_mirror,z_mirror)


  if (myid == 0)  print*, '*** Mirror points have been calculated! ***'
  !------------------------------------------------------------
  WP1(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6,:) = 0.0_rp
  WP2(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6,:) = 0.0_rp
  call InterpolationWeights(n,lo,hi, &
                            nabs_surf,Level_set, &
                            x_IP1,y_IP1,z_IP1, &
                            x_IP2,y_IP2,z_IP2, &
                            i_IP1,j_IP1,k_IP1, &
                            i_IP2,j_IP2,k_IP2, &
                            WP1,WP2,zc)
  do i = 1,7
   do h = 1,nh
     WP1(:,:,lo(3)-h,i)       = WP1(:,:,lo(3),i)
     WP1(:,:,hi(3)-h,i)       = WP1(:,:,hi(3),i)
     WP2(:,:,lo(3)-h,i)       = WP2(:,:,lo(3),i)
     WP2(:,:,hi(3)-h,i)       = WP2(:,:,hi(3),i)
   end do
  end do

  if (myid.eq.0)  print*, '*** Interpolation Weights have been calculated! ***'

  deallocate(x_IP1,x_IP2,y_IP1,y_IP2,z_IP1,z_IP2,tmp)
#endif
#if defined(_IBM_BC)
  call loadIBM('w',trim(datadir)//'IBM.bin',n,cell_phi_tag,cell_u_tag,cell_v_tag,cell_w_tag,level_set,nx_surf, &
               ny_surf,nz_surf,nabs_surf,deltan,i_mirror,j_mirror,k_mirror, &
               i_IP1,j_IP1,k_IP1,i_IP2,j_IP2,k_IP2, &
               WP1,WP2)
#else
   dummy_time = 0.; dummy_istep = 0
   !$acc update self(cell_u_tag,cell_v_tag,cell_w_tag,cell_phi_tag)
   call load_all('w',trim(datadir)//'IBM.bin',MPI_COMM_WORLD,ng,[1,1,1],lo,hi,dummy_time,dummy_istep,cell_u_tag,cell_v_tag,cell_w_tag,cell_phi_tag)
#endif
else
#if defined(_IBM_BC)
  call loadIBM('r',trim(datadir)//'IBM.bin',n,cell_phi_tag,cell_u_tag,cell_v_tag,cell_w_tag,level_set,nx_surf, &
               ny_surf,nz_surf,nabs_surf,deltan,i_mirror,j_mirror,k_mirror, &
               i_IP1,j_IP1,k_IP1,i_IP2,j_IP2,k_IP2, &
               WP1,WP2)

  ! allocate(tmp(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6))
  !$acc enter data copyin(cell_u_tag,cell_v_tag,cell_w_tag,cell_phi_tag)
  ! tmp(:,:,:) = real(Level_set(:,:,:),rp)
  call bounduvw(cbcvel,n,nh,halo,bcvel,nb,is_bound,.false.,dl,dzc,dzf,cell_u_tag,cell_v_tag,cell_w_tag)
  call boundp(cbcpre,n,nh,halo,bcpre,nb,is_bound,dl,dzc,cell_phi_tag)
  ! call boundp(cbcpre,n,bcpre,nb,is_bound,dl,dzc,tmp)
  ! Level_set(:,:,:) = int(tmp(:,:,:),i8)
  !$acc enter data copyin(nx_surf,ny_surf,nz_surf,nabs_surf,deltan)
  call bounduvw(cbcvel,n,nh,halo,bcvel,nb,is_bound,.false.,dl,dzc,dzf,nx_surf,ny_surf,nz_surf)
  call boundp(cbcpre,n,nh,halo,bcpre,nb,is_bound,dl,dzc,nabs_surf)
  call boundp(cbcpre,n,nh,halo,bcpre,nb,is_bound,dl,dzc,deltan)

   do h = 1,nh
      i_mirror(:,:,lo(3)-h)    = i_mirror(:,:,lo(3))
      i_mirror(:,:,hi(3)+h)    = i_mirror(:,:,hi(3))
      j_mirror(:,:,lo(3)-h)    = j_mirror(:,:,lo(3))
      j_mirror(:,:,hi(3)+h)    = j_mirror(:,:,hi(3))
      k_mirror(:,:,lo(3)-h)    = k_mirror(:,:,lo(3))
      k_mirror(:,:,hi(3)+h)    = k_mirror(:,:,hi(3))
      !
      i_IP1(:,:,lo(3)-h)       = i_IP1(:,:,lo(3))
      i_IP1(:,:,hi(3)+h)       = i_IP1(:,:,hi(3))
      j_IP1(:,:,lo(3)-h)       = j_IP1(:,:,lo(3))
      j_IP1(:,:,hi(3)+h)       = j_IP1(:,:,hi(3))
      k_IP1(:,:,lo(3)-h)       = k_IP1(:,:,lo(3))
      k_IP1(:,:,hi(3)+h)       = k_IP1(:,:,hi(3))
      !
      i_IP2(:,:,lo(3)-h)       = i_IP2(:,:,lo(3))
      i_IP2(:,:,hi(3)+h)       = i_IP2(:,:,hi(3))
      j_IP2(:,:,lo(3)-h)       = j_IP2(:,:,lo(3))
      j_IP2(:,:,hi(3)+h)       = j_IP2(:,:,hi(3))
      k_IP2(:,:,lo(3)-h)       = k_IP2(:,:,lo(3))
      k_IP2(:,:,hi(3)+h)       = k_IP2(:,:,hi(3))
   end do
#else
   call load_all('r',trim(datadir)//'IBM.bin',MPI_COMM_WORLD,ng,[1,1,1],lo,hi,dummy_time,dummy_istep,cell_u_tag,cell_v_tag,cell_w_tag,          cell_phi_tag)
   ! allocate(tmp(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6))
   !$acc enter data copyin(cell_u_tag,cell_v_tag,cell_w_tag,cell_phi_tag)
   ! tmp(:,:,:) = real(Level_set(:,:,:),rp)
   call boundp(cbcpre,n,nh,halo,bcpre,nb,is_bound,dl,dzc,cell_phi_tag)
   call boundp(cbcpre,n,nh,halo,bcpre,nb,is_bound,dl,dzc,cell_u_tag)
   call boundp(cbcpre,n,nh,halo,bcpre,nb,is_bound,dl,dzc,cell_v_tag)
   call boundp(cbcpre,n,nh,halo,bcpre,nb,is_bound,dl,dzc,cell_w_tag)
   ! call boundp(cbcpre,n,bcpre,nb,is_bound,dl,dzc,tmp)
   ! Level_set(:,:,:) = int(tmp(:,:,:),i8)
#endif
     if (myid == 0)  print*, '*** Saved IBM data loaded ***'
  !---------------------------------------------------------------------
endif
return
end subroutine initIBM
#endif
#endif
end module mod_initIBM