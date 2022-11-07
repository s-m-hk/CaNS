module mod_IBM
#if defined(_IBM)
use mpi
use decomp_2d
use mod_param, only: lx,ly,lz,dx,dy,dxi,dyi,dims, &
                     surface_type,solid_height_ratio,Rotation_angle, &
                     sx, sy, sz, depth, rod, &
                     small, nb
use mod_common_mpi, only: myid, ipencil => ipencil_axis
use mod_types
!@acc use openacc
!@acc use cudecomp
implicit none
private
#if defined(_SIMPLE)
public IBM_Mask
#else
public IBM_Mask, normal_vectors,intersect,mirrorpoints, &
       interpolation_2D_velocity,interpolation_dphi, &
       mirrorpoints_ijk, interpolation_mirror,InterpolationWeights
#endif
contains
!
#if defined(_SIMPLE)
subroutine IBM_Mask(n,ng,lo,hi,zc,zf,zf_g,dzc,dzf,cell_phi_tag)
implicit none
integer ,intent(in), dimension(3) :: n,ng,lo,hi
real(rp),intent(in), dimension(0:) :: zc,zf,zf_g,dzc,dzf
real(rp),intent(out),dimension(0:,0:,0:) :: cell_phi_tag
real(rp)::xxx,yyy,zzz,ratio
integer i,j,k,nx,ny,nz
logical :: ghost
nx = n(1)
ny = n(2)
nz = n(3)
#if !defined(_DECOMP_Z)
ratio = 1.
#else
ratio = solid_height_ratio
#endif
if(trim(surface_type) == 'Lattice') then
do k=1,int(ratio*nz)
 do j=1,n(2)
  do i=1,n(1)
     xxx = (i+lo(1)-1-.5)*dx
     yyy = (j+lo(2)-1-.5)*dy
     zzz = zc(k)
     ghost = lattice(xxx,yyy,zzz,zf_g,dzc,lx,ly,lz)
     if (ghost) then
       cell_phi_tag(i,j,k) = 1.
     endif
  enddo
 enddo
enddo
endif

end subroutine IBM_Mask
#endif
!
#if defined(_VOLUME)
subroutine IBM_Mask(n,ng,lo,hi,zc,zf,zf_g,dzc,dzf,is_bound,cell_phi_tag,Level_set,surf_height)
implicit none
integer ,intent(in), dimension(3) :: n,ng,lo,hi
logical ,intent(in), dimension(0:1,3) :: is_bound
real(rp),intent(in), dimension(0:) :: zc,zf,zf_g,dzc,dzf
real(rp),intent(in), dimension(1:,1:) :: surf_height
#if defined(_IBM_BC)
real(rp),intent(out),dimension(-5:,-5:,-5:) :: cell_phi_tag
integer, intent(out),dimension(-5:,-5:,-5:) :: Level_set
#else
real(rp),intent(out),dimension(0:,0:,0:) :: cell_phi_tag
integer, intent(out),dimension(0:,0:,0:) :: Level_set
#endif
#if !defined(_DECOMP_Z)
real(rp), allocatable, dimension(:,:,:) :: var_u,var_v,var_w,var_phi
real(rp), allocatable, dimension(:,:,:) :: tmp_x,tmp_y,tmp_z
#endif
integer i,j,k,l,nn,m,number_of_divisions
real(rp):: length_z,xxx,yyy,zzz,dxx,dyy,dzz,dxl,dyl
real(rp):: cell_start_x,cell_end_x,cell_start_y,cell_end_y,cell_start_z,cell_end_z
real(rp):: counter,ratio
logical :: inside,ghost
integer, dimension(3) :: nh
integer :: nx,ny,nz,kk
dxl = dx
dyl = dy
length_z = lz
nx = n(1)
ny = n(2)
nz = n(3)
#if defined(_IBM_BC)
nh(1:3) = 6
#else
nh(1:3) = 1
#endif
#if !defined(_DECOMP_Z)
ratio = 1.
#else
ratio = solid_height_ratio
#endif
!$acc enter data create(xxx,yyy,zzz,dxx,dyy,dzz,cell_start_x,cell_end_x,cell_start_y,cell_end_y,cell_start_z,cell_end_z,ghost,inside) async(1)
! Wall Geometry
#if !defined(_GPU)
number_of_divisions = 50
#else
number_of_divisions = 100
#endif
if (myid == 0) print*, '*** Calculating volume fractions ***'
if(trim(surface_type) == 'HeightMap') then
if(.not.is_bound(1,3)) then
!$acc parallel loop gang collapse(3) default(present) private(xxx,yyy,zzz,dxx,dyy,dzz,cell_start_x,cell_end_x,cell_start_y,cell_end_y,cell_start_z,cell_end_z,ghost,inside) async(1)
do k=1,int(ratio*nz) ! Lower wall
  ! if (myid == 0) print*, '*** Calculating volume fractions at k = ', k
  do j=1,ny
    do i=1,nx
      xxx = (i+lo(1)-1-.5)*dxl
      yyy = (j+lo(2)-1-.5)*dyl
      zzz = zc(k)
	  ghost = height_map_ghost(xxx,yyy,zzz,i,j,dzc(k),n,surf_height)
      if (ghost) then
          cell_phi_tag(i,j,k) = 1.
          cell_u_tag(i,j,k)   = 1.
          cell_v_tag(i,j,k)   = 1.
          cell_w_tag(i,j,k)   = 1.
          Level_set(i,j,k)    = 1
      endif

! Cell Center
      inside = .false.
      cell_start_x = (i+lo(1)-1-1.)*dxl
      cell_end_x   = (i+lo(1)-1-.0)*dxl

      cell_start_y = (j+lo(2)-1-1.)*dyl
      cell_end_y   = (j+lo(2)-1-.0)*dyl

      cell_start_z = zf(k-1)
      cell_end_z   = zf(k)

      dxx = (cell_end_x-cell_start_x)/number_of_divisions
      dyy = (cell_end_y-cell_start_y)/number_of_divisions
      dzz = (cell_end_z-cell_start_z)/number_of_divisions

      counter = 0
      !$acc loop seq
      do nn= 1,number_of_divisions
         zzz = cell_start_z + (nn-1)*dzz
        !$acc loop seq
        do m = 1,number_of_divisions
           yyy = cell_start_y + (m-1 )*dyy
            !$acc loop seq
            do l = 1,number_of_divisions
               xxx = cell_start_x + (l-1 )*dxx
	           inside = height_map(xxx,yyy,zzz,dxl,dyl,i,j,dzc(k),n,surf_height,lo,hi)
               if (inside) counter = counter +1
            enddo
        enddo
      enddo
     cell_phi_tag(i,j,k) = counter/(1.*number_of_divisions**3) !Solid volume fraction

    enddo
  enddo
enddo
endif
if(.not.is_bound(0,3)) then
#if !defined(_DECOMP_Z)
!$acc parallel loop gang collapse(3) default(present) private(xxx,yyy,zzz,dxx,dyy,dzz,cell_start_x,cell_end_x,cell_start_y,cell_end_y,cell_start_z,cell_end_z,ghost,inside) async(1)
do k=nz,(nz-int(ratio*nz)),-1 ! Upper wall
  do j=1,ny
    do i=1,nx
      xxx = (i+lo(1)-1-.5)*dxl
      yyy = (j+lo(2)-1-.5)*dyl
      zzz = length_z - zc(k)
	  ghost = height_map_ghost(xxx,yyy,zzz,i,j,dzc(k),n,surf_height)
      if (ghost) then
        cell_phi_tag(i,j,k) = 1.
        cell_u_tag(i,j,k)   = 1.
        cell_v_tag(i,j,k)   = 1.
        cell_w_tag(i,j,k)   = 1.
        Level_set(i,j,k)    = 1
      endif

! Cell Center
      inside = .false.
      cell_start_x = (i+lo(1)-1-1.)*dxl
      cell_end_x   = (i+lo(1)-1-.0)*dxl

      cell_start_y = (j+lo(2)-1-1.)*dyl
      cell_end_y   = (j+lo(2)-1-.0)*dyl

      cell_start_z = length_z - zf(k-1)
      cell_end_z   = length_z - zf(k)

      dxx = (cell_end_x-cell_start_x)/number_of_divisions
      dyy = (cell_end_y-cell_start_y)/number_of_divisions
      dzz = (cell_end_z-cell_start_z)/number_of_divisions

      counter = 0
      !$acc loop seq
      do nn= 1,number_of_divisions
         zzz = cell_start_z + (nn-1)*dzz
        !$acc loop seq
        do m = 1,number_of_divisions
           yyy = cell_start_y + (m-1 )*dyy
            !$acc loop seq
            do l = 1,number_of_divisions
               xxx = cell_start_x + (l-1 )*dxx
	           inside = height_map(xxx,yyy,zzz,dxl,dyl,i,j,dzc(k),n,surf_height,lo,hi)
               if (inside) counter = counter +1
            enddo
        enddo
      enddo
     cell_phi_tag(i,j,k) = counter/(1.*number_of_divisions**3) !Solid volume fraction

    enddo
  enddo
enddo
#endif
endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
else ! Solids generated using functions [not GPU driven for now]
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if(.not.is_bound(1,3)) then
!!$acc parallel loop gang collapse(3) default(present) private(xxx,yyy,zzz,dxx,dyy,dzz,cell_start_x,cell_end_x,cell_start_y,cell_end_y,cell_start_z,cell_end_z,ghost,inside) async(1)
do k=1,int(ratio*nz) ! Lower wall
  do j=1,ny
    do i=1,nx
      xxx = (i+lo(1)-1-.5)*dx
      yyy = (j+lo(2)-1-.5)*dy
      zzz = zc(k)
	  if(trim(surface_type) == 'Wall') then
          ghost = wall_ghost(xxx,yyy,zzz,i,j,dzc(k))
      elseif(trim(surface_type) == 'Sphere') then
          ghost = sphere_ghost(xxx,yyy,zzz,i,j,k,dzc(k))
      endif
      if (ghost) then
          cell_phi_tag(i,j,k) = 1.
          cell_u_tag(i,j,k)   = 1.
          cell_v_tag(i,j,k)   = 1.
          cell_w_tag(i,j,k)   = 1.
          Level_set(i,j,k)    = 1
      endif

! Cell Center
      inside = .false.
      cell_start_x = (i+lo(1)-1-1.)*dxl
      cell_end_x   = (i+lo(1)-1-.0)*dxl

      cell_start_y = (j+lo(2)-1-1.)*dyl
      cell_end_y   = (j+lo(2)-1-.0)*dyl

      cell_start_z = zf(k-1)
      cell_end_z   = zf(k)

      dxx = (cell_end_x-cell_start_x)/number_of_divisions
      dyy = (cell_end_y-cell_start_y)/number_of_divisions
      dzz = (cell_end_z-cell_start_z)/number_of_divisions

      counter = 0
!      !$acc loop seq
      do nn= 1,number_of_divisions
          zzz = cell_start_z+(nn-1)*dzz
!        !$acc loop seq
        do m = 1,number_of_divisions
            yyy = cell_start_y + (m-1)*dyy
!            !$acc loop seq
            do l = 1,number_of_divisions
              xxx = cell_start_x + (l-1)*dxx
	          if(trim(surface_type) == 'Wall') then
                  inside = wall(xxx,yyy,zzz,i,j,dzc(k))
              elseif(trim(surface_type) == 'Sphere') then
                  inside = sphere(xxx,yyy,zzz,i,j,dzc(k),lo,hi)
              endif
              if (inside) counter = counter +1
            enddo
        enddo
      enddo
     cell_phi_tag(i,j,k) = counter/(1.*number_of_divisions**3) !Solid volume fraction

    enddo
  enddo
enddo
!
endif
if(.not.is_bound(0,3)) then
#if !defined(_DECOMP_Z)
!!$acc parallel loop gang collapse(3) default(present) private(xxx,yyy,zzz,dxx,dyy,dzz,cell_start_x,cell_end_x,cell_start_y,cell_end_y,cell_start_z,cell_end_z,ghost,inside) async(1)
do k=nz,(nz-int(ratio*nz)),-1 ! Upper wall
  do j=1,ny
    do i=1,nx
      xxx = (i+lo(1)-1-.5)*dx
      yyy = (j+lo(2)-1-.5)*dy
      zzz = length_z - zc(k)
	  if(trim(surface_type) == 'Wall') then
          ghost = wall_ghost(xxx,yyy,zzz,i,j,dzc(k))
      elseif(trim(surface_type) == 'Sphere') then
          ghost = sphere_ghost(xxx,yyy,zzz,i,j,k,dzc(k))
      endif
      if (ghost) then
          cell_phi_tag(i,j,k) = 1.
          cell_u_tag(i,j,k)   = 1.
          cell_v_tag(i,j,k)   = 1.
          cell_w_tag(i,j,k)   = 1.
          Level_set(i,j,k)    = 1
      endif

! Cell Center
      inside = .false.
      cell_start_x = (i+lo(1)-1-1.)*dxl
      cell_end_x   = (i+lo(1)-1-.0)*dxl

      cell_start_y = (j+lo(2)-1-1.)*dyl
      cell_end_y   = (j+lo(2)-1-.0)*dyl

      cell_start_z = length_z - zf(k-1)
      cell_end_z   = length_z - zf(k)

      dxx = (cell_end_x-cell_start_x)/number_of_divisions
      dyy = (cell_end_y-cell_start_y)/number_of_divisions
      dzz = (cell_end_z-cell_start_z)/number_of_divisions

      counter = 0
!      !$acc loop seq
      do nn= 1,number_of_divisions
          zzz = cell_start_z + (nn-1)*dzz
!        !$acc loop seq
        do m = 1,number_of_divisions
            yyy = cell_start_y + (m-1)*dyy
!            !$acc loop seq
            do l = 1,number_of_divisions
              xxx = cell_start_x + (l-1)*dxx
	          if(trim(surface_type) == 'Wall') then
                  inside = wall(xxx,yyy,zzz,i,j,dzc(k))
              elseif(trim(surface_type) == 'Sphere') then
                  inside = sphere(xxx,yyy,zzz,i,j,dzc(k),lo,hi)
              endif
              if (inside) counter = counter +1
            enddo
        enddo
      enddo
     cell_phi_tag(i,j,k) = counter/(1.*number_of_divisions**3) !Solid volume fraction

    enddo
  enddo
enddo
#endif
endif
!
endif

! Duplicate volume fractions on the upper part of the domain (symmetric channel)
#if defined(_DECOMP_Z)
kk = 1
!$acc parallel loop gang default(present) async(1)
do k = nz,(nz-int(ratio*nz)),-1
     cell_phi_tag(1:n(1),1:n(2),k) = cell_phi_tag(1:n(1),1:n(2),kk)
	 kk = kk + 1
enddo
#endif

end subroutine IBM_Mask
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!! Functions for geometries !!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function lattice(xIn, yIn, zIn, zf_g, dz, lengthx, lengthy, lengthz)
  implicit none
  logical :: lattice,read_z
  real(rp), dimension(0:), intent(in) :: zf_g,dz
  real(rp), intent(in) :: xIn,yIn,zIn,lengthx,lengthy,lengthz
  real(rp) :: x0,y0,x1,y1,z0,z1,z_offset
  real(rp) :: d
  integer :: i,j,k,N_x,N_y,N_z
	   
	   lattice = .false.
	   
	   d = rod  !Rod thickness
	   
	   N_z = INT(depth/sz)
	   N_y = INT(lengthy/sy)
	   N_x = INT(lengthx/sx)

       inquire(file='zf.dat',exist=read_z)
	   if (read_z) then
	    z_offset = zIn - zf_g(0)
       else
	    z_offset = zIn
       endif
	   
	  if (z_offset.le.(d+(N_z*sz))) then
	   
	  !Streamwise
	   x0 = 0.; x1 = d; y0 = 0.; y1 = d
	   do i = 0,N_x
	    if ((xIn.gt.x0) .and. (xIn.lt.x1)) then
	     do j = 0,N_y
		  if ((yIn.gt.y0).and.(yIn.lt.y1)) then
		   lattice = .true.
		  endif
		  y0 = y0 + sy
		  y1 = y1 + sy
		 enddo
	    endif
		x0 = x0 + sx
		x1 = x1 + sx
	   enddo
	   
	   x0 = 0.; x1 = d; z0 = 0.; z1 = INT(d/dz(1))*dz(1)
	   do i = 0,N_x
	    if ((xIn.gt.x0) .and. (xIn.lt.x1)) then
	     do k = 0,N_z
		  if ((z_offset.gt.z0).and.(z_offset.lt.z1)) then
		   lattice = .true.
		  endif
		  z0 = z0 + sz
		  z1 = z1 + sz
		 enddo
	    endif
		x0 = x0 + sx
		x1 = x1 + sx
	   enddo
	   
	  !Spanwise
	   y0 = 0.; y1 = d; z0 = 0.; z1 = INT(d/dz(1))*dz(1)
	   do j = 0,N_y
	    if ((yIn.gt.y0) .and. (yIn.lt.y1)) then
	     do k = 0,N_z
		  if ((z_offset.gt.z0).and.(z_offset.lt.z1)) then
		   lattice = .true.
		  endif
		  z0 = z0 + sz
		  z1 = z1 + sz
		 enddo
	    endif
		y0 = y0 + sy
		y1 = y1 + sy
	   enddo
	   
	   y0 = 0.; y1 = d; x0 = 0.; x1 = d
	   do j = 0,N_y
	    if ((yIn.gt.y0) .and. (yIn.lt.y1)) then
	     do i = 0,N_x
		  if ((xIn.gt.x0).and.(xIn.lt.x1)) then
		   lattice = .true.
		  endif
		  x0 = x0 + sx
		  x1 = x1 + sx
		 enddo
	    endif
		y0 = y0 + sy
		y1 = y1 + sy
	   enddo
	   
	  !Wall-normal
	   z0 = 0.; z1 = INT(d/dz(1))*dz(1); y0 = 0.; y1 = d
	   do k = 0,N_z
	    if ((z_offset.gt.z0) .and. (z_offset.lt.z1)) then
	     do j = 0,N_y
		  if ((yIn.gt.y0).and.(yIn.lt.y1)) then
		   lattice = .true.
		  endif
		  y0 = y0 + sy
		  y1 = y1 + sy
		 enddo
	    endif
		z0 = z0 + sz
		z1 = z1 + sz
	   enddo
	   
	   z0 = 0.; z1 = INT(d/dz(1))*dz(1); x0 = 0.; x1 = d
	   do k = 0,N_z
	    if ((z_offset.gt.z0) .and. (z_offset.lt.z1)) then
	     do i = 0,N_x
		  if ((xIn.gt.x0).and.(xIn.lt.x1)) then
		   lattice = .true.
		  endif
		  x0 = x0 + sx
		  x1 = x1 + sx
		 enddo
	    endif
		z0 = z0 + sz
		z1 = z1 + sz
	   enddo
	   
	   endif

end function lattice

function height_map(xxx,yyy,zzz,dxl,dyl,ii,jj,dz,n,surf_height,lo,hi)
implicit none
logical :: height_map,cond1,cond2,cond3
integer ,intent(in), dimension(3) :: n,lo,hi
real(rp), intent(in ), dimension(1:n(1),1:n(2)) :: surf_height
real(rp), intent(in):: xxx,yyy,zzz,dxl,dyl,dz
integer, intent(in):: ii,jj

     height_map=.false.
	 cond1 = zzz.le.surf_height(ii,jj)
	 cond2 = abs((ii+lo(1)-1-.5)*dxl-xxx).le.(0.5*dxl)
	 cond3 = abs((jj+lo(2)-1-.5)*dyl-yyy).le.(0.5*dyl)
	 if (cond1.and.cond2.and.cond3) height_map=.true.

end function height_map
!
function height_map_ghost(xxx,yyy,zzz,ii,jj,dz,n,surf_height)
implicit none
logical :: height_map_ghost,cond1,cond2,cond3
integer ,intent(in), dimension(3) :: n
real(rp), intent(in ), dimension(1:n(1),1:n(2)) :: surf_height
real(rp), intent(in):: xxx,yyy,zzz,dz
integer, intent(in):: ii,jj

     height_map_ghost=.false.
	 if (zzz.lt.surf_height(ii,jj)) height_map_ghost=.true.

end function height_map_ghost

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function wall(xxx,yyy,zzz,ii,jj,dz)
implicit none
logical :: wall,cond1,cond2,cond3
real(rp), intent(in):: xxx,yyy,zzz,dz
real(rp):: x_center,y_center,z_center
integer, intent(in):: ii,jj

     wall=.false.
     z_center = lz/40.
     if (zzz.le.z_center) wall=.true.
	 
end function wall
!
function wall_ghost(xxx,yyy,zzz,ii,jj,dz)
implicit none
logical :: wall_ghost,cond1
real(rp), intent(in):: xxx,yyy,zzz,dz
real(rp):: x_center,y_center,z_center
integer, intent(in):: ii,jj

     wall_ghost=.false.
     z_center = lz/40.

     if (zzz.lt.z_center) wall_ghost=.true.

end function wall_ghost

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function sphere(xxx,yyy,zzz,ii,jj,dz,lo,hi)
implicit none
logical :: sphere,cond1,cond2,cond3,cond4
real(rp), intent(in):: xxx,yyy,zzz,dz
real(rp):: RR,LL,x_center,y_center,z_center
integer ,intent(in), dimension(3) :: lo,hi
integer, intent(in):: ii,jj

     sphere=.false.
     RR =  lz/4.
     x_center = lx/2.
     y_center = ly/2.
     z_center = lz/8. + 0.0002
	 
     LL = (zzz-z_center)*(zzz-z_center)+(yyy-y_center)*(yyy-y_center)+(xxx-x_center)*(xxx-x_center)
	 cond1 = LL.le.(RR*RR)
	 cond2 = zzz.le.0.0002
	 cond3 = abs((ii+lo(1)-1-.5)*dx-xxx).le.(0.5*dx)
	 cond4 = abs((jj+lo(2)-1-.5)*dy-yyy).le.(0.5*dy)
     if (cond1.and.cond2.and.cond3.and.cond4) sphere=.true.

end function sphere
!
function sphere_ghost(xxx,yyy,zzz,ii,jj,kk,dz)
implicit none
logical :: sphere_ghost,cond1,cond2,cond3,cond4
real(rp), intent(in):: xxx,yyy,zzz,dz
real(rp) :: x1,x2,y1,y2,z1,z2,RR,LL,x_center,y_center,z_center, &
            threshold
integer, intent(in):: ii,jj,kk

     sphere_ghost=.false.

     RR       = lz/4.
     x_center = lx/2.
     y_center = ly/2.
     z_center = lz/8. + 0.0002
     
     LL = (zzz-z_center)*(zzz-z_center)+(yyy-y_center)*(yyy-y_center)+(xxx-x_center)*(xxx-x_center)
     if (LL.lt.(RR*RR)) sphere_ghost=.true.
	 
	 ! threshold = 1.5*sqrt(dx(ii)**2+dy(jj)**2+dz(kk)**2+eps)
     ! RR =  zLength/4.0
     ! x_center = xLength/2.0
     ! y_center = yLength/2.0
     ! z_center = (zLength/8.0) + 0.0002
     ! x1 = x_center + sqrt(RR*RR-(zzz-z_center)*(zzz-z_center)-(yyy-y_center)*(yyy-y_center)+eps)
     ! x2 = x_center - sqrt(RR*RR-(zzz-z_center)*(zzz-z_center)-(yyy-y_center)*(yyy-y_center)+eps)
     
     ! y1 = y_center + sqrt(RR*RR-(zzz-z_center)*(zzz-z_center)-(xxx-x_center)*(xxx-x_center)+eps)
     ! y2 = y_center - sqrt(RR*RR-(zzz-z_center)*(zzz-z_center)-(xxx-x_center)*(xxx-x_center)+eps)
     
     ! z1 = z_center + sqrt(RR*RR-(xxx-x_center)*(xxx-x_center)-(yyy-y_center)*(yyy-y_center)+eps)
     ! z2 = z_center - sqrt(RR*RR-(xxx-x_center)*(xxx-x_center)-(yyy-y_center)*(yyy-y_center)+eps)
     ! cond1 = (abs(x1-xxx).lt.threshold).or.(abs(x2-xxx).lt.threshold)
     ! cond2 = (abs(y1-yyy).lt.threshold).or.(abs(y2-yyy).lt.threshold)
     ! cond3 = (abs(z1-zzz).lt.threshold).or.(abs(z2-zzz).lt.threshold)
	 ! cond4 = zzz.lt.0.0002
     ! if (cond1.or.cond2.or.cond3.or.cond4) sphere_ghost=.true.

end function sphere_ghost

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function wavywall(xxx,yyy,zzz,ii,jj,dz,lo,hi) !Maaß and Schumann, ERCOFTAC classic case
implicit none
logical :: wavywall,cond1,cond2,cond3
real(rp), intent(in):: xxx,yyy,zzz,dz
real(rp):: AA
integer ,intent(in), dimension(3) :: lo,hi
integer, intent(in):: ii,jj

     wavywall=.false.
     AA = 0.05*cos(2.*(4.*atan(1.))*((ii+lo(1)-1-.5)*dx))
	 cond1 = zzz.le.AA
	 cond2 = abs((ii+lo(1)-1-.5)*dx-xxx).le.(0.5*dx)
	 cond3 = abs((jj+lo(2)-1-.5)*dy-yyy).le.(0.5*dy)
	 if (cond1.and.cond2.and.cond3) wavywall=.true.

end function wavywall
!
function wavywall_ghost(xxx,yyy,zzz,ii,jj,dz,lo,hi)
implicit none
logical :: wavywall_ghost,cond1,cond2
real(rp), intent(in):: xxx,yyy,zzz,dz
real(rp):: AA
integer ,intent(in), dimension(3) :: lo,hi
integer, intent(in):: ii,jj

     AA = 0.05*2.*cos(2.*(4.*atan(1._rp))*((ii+lo(1)-1-.5)*dx/2.))
	 wavywall_ghost=.false.
	 cond1 = (AA.gt.zzz)
	 cond2 = (abs(AA-zzz).lt.dz)
     if (cond1.and.cond2) wavywall_ghost=.true.

end function wavywall_ghost
#if defined(_VOLUME)
Subroutine normal_vectors(lo,hi,Level_set,cell_phi_tag,nx_surf,ny_surf,nz_surf,nabs_surf, &
                          nx_surf_nonnorm,ny_surf_nonnorm,nz_surf_nonnorm,zc,zf,dzc,dzf,dl,dli,n,surf_height)
implicit none
integer ,dimension(3),intent(in) :: n,lo,hi
real(rp),dimension(1:n(1),1:n(2)),intent(in),optional :: surf_height
real(rp),dimension(0:n(3)+1),intent(in) :: zc,zf,dzc,dzf,dl,dli
real(rp),dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6),intent(in)  :: cell_phi_tag
integer ,dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6),intent(in)  :: Level_set
real(rp),dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6),intent(out) :: nx_surf
real(rp),dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6),intent(out) :: ny_surf
real(rp),dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6),intent(out) :: nz_surf
real(rp),dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6),intent(out) :: nabs_surf
real(rp),dimension(0:n(1)+1,0:n(2)+1,0:n(3)+1),intent(out):: nx_surf_nonnorm,ny_surf_nonnorm,nz_surf_nonnorm
integer::i,j,k,l
real(rp):: nx, ny, nz, n_abs,n_abs_p
real(rp):: m_av_x,m_av_y,m_av_z,normal_denum
real(rp):: mx1,mx2,mx3,mx4,mx5,mx6,mx7,mx8
real(rp):: my1,my2,my3,my4,my5,my6,my7,my8
real(rp):: mz1,mz2,mz3,mz4,mz5,mz6,mz7,mz8
real(rp):: cell_tag_wall = 0.5
real(rp) :: dli1, dli2, dli3
real(rp) :: dl1, dl2, dl3
integer:: ip,jp,kp
integer:: im,jm,km
integer:: iii,jjj,kkk
logical :: inside_1,inside_2,inside_3,inside_4,inside_5,inside_6,inside_7,inside_8,inside_9
logical :: inside_10,inside_11,inside_12,inside_13,inside_14,inside_15,inside_16,inside_17,inside_18
logical :: inside_19,inside_20,inside_21,inside_22,inside_23,inside_24,inside_25,inside_26,inside_27
logical ::inside_jpkm,inside_jmkp,inside_jmkm,inside_jk,ghost_cond, inside
real(rp)::xxx,yyy,zzz,dxl,dyl
integer, dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6) :: ghost_cell_tag
integer:: Level_set_all
logical:: ghost

dli1=dli(1)
dli2=dli(2)
dli3=dli(3)
dl1=dl(1)
dl2=dl(2)
dl3=dl(3)
dxl = dx
dyl = dy

! Preparing normal vectors to the surfaces

! Identifying the ghost cells
ghost_cell_tag(:,:,:) = 0
do k=1,int(solid_height_ratio*(n(3)+1))
if (myid.eq.0) print*, 'Calculating Normal Vectors at k = ', k
  do j=1,n(2)
   do i=1,n(1)
            ghost = .false.
            xxx   =  (i+lo(1)-1-.5)*dx
            yyy   =  (j+lo(2)-1-.5)*dy
            zzz   =  zc(k)
            if(trim(surface_type) == 'HeightMap') then 
	            ghost = height_map_ghost(xxx,yyy,zzz,i,j,dzc(k),n,surf_height)
	        elseif(trim(surface_type) == 'Wall') then
                ghost = wall_ghost(xxx,yyy,zzz,i,j,dzc(k))
            elseif(trim(surface_type) == 'Sphere') then
                ghost = sphere_ghost(xxx,yyy,zzz,i,j,k,dzc(k))
            endif
            if (.not.ghost)  cycle
            ip = i+1
            jp = j+1
            kp = k+1
            im = i-1
            jm = j-1
            km = k-1
            ghost_cond = .false.
            Level_set_all = Level_set(im,jm,kp) + Level_set(im,j,kp)+Level_set(im,jp,kp)+ &
                            Level_set(i,jm,kp)  + Level_set(i,j,kp) +Level_set(i,jp,kp) + &
                            Level_set(ip,jm,kp) + Level_set(ip,j,kp)+Level_set(ip,jp,kp)+ &
                            Level_set(im,jm,k)  + Level_set(im,j,k) +Level_set(im,jp,k) + &
                            Level_set(i,jm,k)   + Level_set(i,j,k)  +Level_set(i,jp,k)  + &
                            Level_set(ip,jm,k)  + Level_set(ip,j,k) +Level_set(ip,jp,k) + &
                            Level_set(im,jm,kp) + Level_set(im,j,kp)+Level_set(im,jp,kp)+ &
                            Level_set(i,jm,km)  + Level_set(i,j,km) +Level_set(i,jp,km) + &
                            Level_set(ip,jm,km) + Level_set(ip,j,km)+Level_set(ip,jp,km)

            if ((Level_set_all.gt.0).and.(Level_set(i,j,k)).eq.0) ghost_cond =.true.
            if (ghost_cond)  ghost_cell_tag(i,j,k) = 1
   enddo
 enddo
enddo

do k=1,n(3)
 do j=1,n(2)
  do i=1,n(1)
   if (ghost_cell_tag(i,j,k) .eq. 1) then

      nx     = 0.
      ny     = 0.
      nz     = 0.
      n_abs  = 0.
      m_av_x = 0.
      m_av_y = 0.
      m_av_z = 0.
      normal_denum = 0.

          !i+1/2 j+1/2 k+1/2
          mx1 = ((cell_phi_tag(i+1,j  ,k  )+cell_phi_tag(i+1,j+1,k  )+cell_phi_tag(i+1,j  ,k+1)+cell_phi_tag(i+1,j+1,k+1)) - &
                   (cell_phi_tag(i  ,j  ,k  )+cell_phi_tag(i  ,j+1,k  )+cell_phi_tag(i  ,j  ,k+1)+cell_phi_tag(i  ,j+1,k+1)))*dli1*0.25_rp
          !i+1/2 j-1/2 k+1/2
          mx2 = ((cell_phi_tag(i+1,j  ,k  )+cell_phi_tag(i+1,j-1,k  )+cell_phi_tag(i+1,j  ,k+1)+cell_phi_tag(i+1,j-1,k+1)) - &
                   (cell_phi_tag(i  ,j  ,k  )+cell_phi_tag(i  ,j-1,k  )+cell_phi_tag(i  ,j  ,k+1)+cell_phi_tag(i  ,j-1,k+1)))*dli1*0.25_rp
          !i+1/2 j+1/2 k-1/2
          mx3 = ((cell_phi_tag(i+1,j  ,k  )+cell_phi_tag(i+1,j+1,k  )+cell_phi_tag(i+1,j  ,k-1)+cell_phi_tag(i+1,j+1,k-1)) - &
                   (cell_phi_tag(i  ,j  ,k  )+cell_phi_tag(i  ,j+1,k  )+cell_phi_tag(i  ,j  ,k-1)+cell_phi_tag(i  ,j+1,k-1)))*dli1*0.25_rp
          !i+1/2 j-1/2 k-1/2
          mx4 = ((cell_phi_tag(i+1,j  ,k  )+cell_phi_tag(i+1,j-1,k  )+cell_phi_tag(i+1,j  ,k-1)+cell_phi_tag(i+1,j-1,k-1)) - &
                   (cell_phi_tag(i  ,j  ,k  )+cell_phi_tag(i  ,j-1,k  )+cell_phi_tag(i  ,j  ,k-1)+cell_phi_tag(i  ,j-1,k-1)))*dli1*0.25_rp
          !i-1/2 j+1/2 k+1/2
          mx5 = ((cell_phi_tag(i  ,j  ,k  )+cell_phi_tag(i  ,j+1,k  )+cell_phi_tag(i  ,j  ,k+1)+cell_phi_tag(i  ,j+1,k+1)) - &
                   (cell_phi_tag(i-1,j  ,k  )+cell_phi_tag(i-1,j+1,k  )+cell_phi_tag(i-1,j  ,k+1)+cell_phi_tag(i-1,j+1,k+1)))*dli1*0.25_rp
          !i-1/2 j-1/2 k+1/2
          mx6 = ((cell_phi_tag(i  ,j  ,k  )+cell_phi_tag(i  ,j-1,k  )+cell_phi_tag(i  ,j  ,k+1)+cell_phi_tag(i  ,j-1,k+1)) - &
                   (cell_phi_tag(i-1,j  ,k  )+cell_phi_tag(i-1,j-1,k  )+cell_phi_tag(i-1,j  ,k+1)+cell_phi_tag(i-1,j-1,k+1)))*dli1*0.25_rp
          !i-1/2 j+1/2 k-1/2
          mx7 = ((cell_phi_tag(i  ,j  ,k  )+cell_phi_tag(i  ,j+1,k  )+cell_phi_tag(i  ,j  ,k-1)+cell_phi_tag(i  ,j+1,k-1)) - &
                   (cell_phi_tag(i-1,j  ,k  )+cell_phi_tag(i-1,j+1,k  )+cell_phi_tag(i-1,j  ,k-1)+cell_phi_tag(i-1,j+1,k-1)))*dli1*0.25_rp
          !i-1/2 j-1/2 k-1/2
          mx8 = ((cell_phi_tag(i  ,j  ,k  )+cell_phi_tag(i  ,j-1,k  )+cell_phi_tag(i  ,j  ,k-1)+cell_phi_tag(i  ,j-1,k-1)) - &
                   (cell_phi_tag(i-1,j  ,k  )+cell_phi_tag(i-1,j-1,k  )+cell_phi_tag(i-1,j  ,k-1)+cell_phi_tag(i-1,j-1,k-1)))*dli1*0.25_rp
          !
      m_av_x= 0.125_rp*(mx1+mx2+mx3+mx4+mx5+mx6+mx7+mx8)
          ! 
          !i+1/2 j+1/2 k+1/2
          my1 = ((cell_phi_tag(i  ,j+1,k  )+cell_phi_tag(i+1,j+1,k  )+cell_phi_tag(i  ,j+1,k+1)+cell_phi_tag(i+1,j+1,k+1)) - &
                   (cell_phi_tag(i  ,j  ,k  )+cell_phi_tag(i+1,j  ,k  )+cell_phi_tag(i  ,j  ,k+1)+cell_phi_tag(i+1,j  ,k+1)))*dli2*0.25_rp
          !i+1/2 j-1/2 k+1/2
          my2 = ((cell_phi_tag(i  ,j  ,k  )+cell_phi_tag(i+1,j  ,k  )+cell_phi_tag(i  ,j  ,k+1)+cell_phi_tag(i+1,j  ,k+1)) - &
                   (cell_phi_tag(i  ,j-1,k  )+cell_phi_tag(i+1,j-1,k  )+cell_phi_tag(i  ,j-1,k+1)+cell_phi_tag(i+1,j-1,k+1)))*dli2*0.25_rp
          !i+1/2 j+1/2 k-1/2
          my3 = ((cell_phi_tag(i  ,j+1,k  )+cell_phi_tag(i+1,j+1,k  )+cell_phi_tag(i  ,j+1,k-1)+cell_phi_tag(i+1,j+1,k-1)) - &
                   (cell_phi_tag(i  ,j  ,k  )+cell_phi_tag(i+1,j  ,k  )+cell_phi_tag(i  ,j  ,k-1)+cell_phi_tag(i+1,j  ,k-1)))*dli2*0.25_rp
          !i+1/2 j-1/2 k-1/2
          my4 = ((cell_phi_tag(i  ,j  ,k  )+cell_phi_tag(i+1,j  ,k  )+cell_phi_tag(i  ,j  ,k-1)+cell_phi_tag(i+1,j  ,k-1)) - &
                   (cell_phi_tag(i  ,j-1,k  )+cell_phi_tag(i+1,j-1,k  )+cell_phi_tag(i  ,j-1,k-1)+cell_phi_tag(i+1,j-1,k-1)))*dli2*0.25_rp
          !i-1/2 j+1/2 k+1/2
          my5 = ((cell_phi_tag(i  ,j+1,k  )+cell_phi_tag(i-1,j+1,k  )+cell_phi_tag(i  ,j+1,k+1)+cell_phi_tag(i-1,j+1,k+1)) - &
                   (cell_phi_tag(i  ,j  ,k  )+cell_phi_tag(i-1,j  ,k  )+cell_phi_tag(i  ,j  ,k+1)+cell_phi_tag(i-1,j  ,k+1)))*dli2*0.25_rp
          !i-1/2 j-1/2 k+1/2
          my6 = ((cell_phi_tag(i  ,j  ,k  )+cell_phi_tag(i-1,j  ,k  )+cell_phi_tag(i  ,j  ,k+1)+cell_phi_tag(i-1,j  ,k+1)) - &
                   (cell_phi_tag(i  ,j-1,k  )+cell_phi_tag(i-1,j-1,k  )+cell_phi_tag(i  ,j-1,k+1)+cell_phi_tag(i-1,j-1,k+1)))*dli2*0.25_rp
          !i-1/2 j+1/2 k-1/2
          my7 = ((cell_phi_tag(i  ,j+1,k  )+cell_phi_tag(i-1,j+1,k  )+cell_phi_tag(i  ,j+1,k-1)+cell_phi_tag(i-1,j+1,k-1)) - &
                   (cell_phi_tag(i  ,j  ,k  )+cell_phi_tag(i-1,j  ,k  )+cell_phi_tag(i  ,j  ,k-1)+cell_phi_tag(i-1,j  ,k-1)))*dli2*0.25_rp
          !i-1/2 j-1/2 k-1/2
          my8 = ((cell_phi_tag(i  ,j  ,k  )+cell_phi_tag(i-1,j  ,k  )+cell_phi_tag(i  ,j  ,k-1)+cell_phi_tag(i-1,j  ,k-1)) - &
                   (cell_phi_tag(i  ,j-1,k  )+cell_phi_tag(i-1,j-1,k  )+cell_phi_tag(i  ,j-1,k-1)+cell_phi_tag(i-1,j-1,k-1)))*dli2*0.25_rp
          !
      m_av_y= 0.125_rp*(my1+my2+my3+my4+my5+my6+my7+my8)
          !
          !i+1/2 j+1/2 k+1/2
          mz1 = ((cell_phi_tag(i  ,j  ,k+1)+cell_phi_tag(i+1,j  ,k+1)+cell_phi_tag(i  ,j+1,k+1)+cell_phi_tag(i+1,j+1,k+1)) - &
                   (cell_phi_tag(i  ,j  ,k  )+cell_phi_tag(i+1,j  ,k  )+cell_phi_tag(i  ,j+1,k  )+cell_phi_tag(i+1,j+1,k  )))*dli3*0.25_rp
          !i+1/2 j-1/2 k+1/2
          mz2 = ((cell_phi_tag(i  ,j  ,k+1)+cell_phi_tag(i+1,j  ,k+1)+cell_phi_tag(i  ,j-1,k+1)+cell_phi_tag(i+1,j-1,k+1)) - &
                   (cell_phi_tag(i  ,j  ,k  )+cell_phi_tag(i+1,j  ,k  )+cell_phi_tag(i  ,j-1,k  )+cell_phi_tag(i+1,j-1,k  )))*dli3*0.25_rp
          !i+1/2 j+1/2 k-1/2
          mz3 = ((cell_phi_tag(i  ,j  ,k  )+cell_phi_tag(i+1,j  ,k  )+cell_phi_tag(i  ,j+1,k  )+cell_phi_tag(i+1,j+1,k  )) - &
                   (cell_phi_tag(i  ,j  ,k-1)+cell_phi_tag(i+1,j  ,k-1)+cell_phi_tag(i  ,j+1,k-1)+cell_phi_tag(i+1,j+1,k-1)))*dli3*0.25_rp
          !i+1/2 j-1/2 k-1/2
          mz4 = ((cell_phi_tag(i  ,j  ,k  )+cell_phi_tag(i+1,j  ,k  )+cell_phi_tag(i  ,j-1,k  )+cell_phi_tag(i+1,j-1,k  )) - &
                   (cell_phi_tag(i  ,j  ,k-1)+cell_phi_tag(i+1,j  ,k-1)+cell_phi_tag(i  ,j-1,k-1)+cell_phi_tag(i+1,j-1,k-1)))*dli3*0.25_rp
          !i-1/2 j+1/2 k+1/2
          mz5 = ((cell_phi_tag(i  ,j  ,k+1)+cell_phi_tag(i-1,j  ,k+1)+cell_phi_tag(i  ,j+1,k+1)+cell_phi_tag(i-1,j+1,k+1)) - &
                   (cell_phi_tag(i  ,j  ,k  )+cell_phi_tag(i-1,j  ,k  )+cell_phi_tag(i  ,j+1,k  )+cell_phi_tag(i-1,j+1,k  )))*dli3*0.25_rp
          !i-1/2 j-1/2 k+1/2
          mz6 = ((cell_phi_tag(i  ,j  ,k+1)+cell_phi_tag(i-1,j  ,k+1)+cell_phi_tag(i  ,j-1,k+1)+cell_phi_tag(i-1,j-1,k+1)) - &
                   (cell_phi_tag(i  ,j  ,k  )+cell_phi_tag(i-1,j  ,k  )+cell_phi_tag(i  ,j-1,k  )+cell_phi_tag(i-1,j-1,k  )))*dli3*0.25_rp
          !i-1/2 j+1/2 k-1/2
          mz7 = ((cell_phi_tag(i  ,j  ,k  )+cell_phi_tag(i-1,j  ,k  )+cell_phi_tag(i  ,j+1,k  )+cell_phi_tag(i-1,j+1,k  )) - &
                   (cell_phi_tag(i  ,j  ,k-1)+cell_phi_tag(i-1,j  ,k-1)+cell_phi_tag(i  ,j+1,k-1)+cell_phi_tag(i-1,j+1,k-1)))*dli3*0.25_rp
          !i-1/2 j-1/2 k-1/2
          mz8 = ((cell_phi_tag(i  ,j  ,k  )+cell_phi_tag(i-1,j  ,k  )+cell_phi_tag(i  ,j-1,k  )+cell_phi_tag(i-1,j-1,k  )) - &
                   (cell_phi_tag(i  ,j  ,k-1)+cell_phi_tag(i-1,j  ,k-1)+cell_phi_tag(i  ,j-1,k-1)+cell_phi_tag(i-1,j-1,k-1)))*dli3*0.25_rp
          !
      m_av_z = 0.125_rp*(mz1+mz2+mz3+mz4+mz5+mz6+mz7+mz8)
      
      nx_surf_nonnorm(i,j,k) = m_av_x
      ny_surf_nonnorm(i,j,k) = m_av_y
      nz_surf_nonnorm(i,j,k) = m_av_z

      normal_denum = sqrt(m_av_x**2 + m_av_y**2 + m_av_z**2)

      nx_surf(i,j,k) = m_av_x/max(normal_denum, small)
      ny_surf(i,j,k) = m_av_y/max(normal_denum, small)
      nz_surf(i,j,k) = m_av_z/max(normal_denum, small)

      nabs_surf(i,j,k) = sqrt(nx_surf(i,j,k)**2 + ny_surf(i,j,k)**2 + nz_surf(i,j,k)**2)

     endif

    enddo
  enddo
enddo
return

end subroutine normal_vectors

Subroutine intersect(n,lo,hi,nx_surf,ny_surf,nz_surf,nabs_surf,x_intersect,y_intersect,z_intersect,zc,dzc,surf_height)

implicit none
integer,dimension(3),intent(in) :: n,lo,hi
real(rp),dimension(1:n(1),1:n(2)),intent(in),optional :: surf_height
real(rp),dimension(0:n(3)+1),intent(in) :: zc,dzc
real(rp),dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6),intent(in) ::  nx_surf
real(rp),dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6),intent(in) ::  ny_surf
real(rp),dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6),intent(in) ::  nz_surf
real(rp),dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6),intent(in) ::  nabs_surf
real(rp),dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6),intent(out):: x_intersect
real(rp),dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6),intent(out):: y_intersect
real(rp),dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6),intent(out):: z_intersect
integer::i,j,k,l
real(rp):: nx, ny, nz,n_abs
real(rp):: step
real(rp):: xxx,yyy,zzz,dxl,dyl
logical :: inside
logical :: confirmation
integer:: lmax
real(rp)::distance_ghost_intersect
real(rp):: x_ghost,y_ghost,z_ghost

dxl = dx
dyl = dy
distance_ghost_intersect= 1000._rp

!***************************************************************************************
! Part one: cell-centered
!***************************************************************************************
do k=1,int(solid_height_ratio*(n(3)+1))
 if (myid.eq.0) print*, 'Calculating Intersect Points at k = ', k
 step = 1.0e-5*dzc(k)
 lmax=100*int(2.*sqrt(3.)*dzc(k)/step)
 do j= 1,n(2)
   do i= 1,n(1)
    if (nabs_surf(i,j,k).gt.small) then
       confirmation = .false.
       nx =   nx_surf(i,j,k)
       ny =   ny_surf(i,j,k)
       nz =   nz_surf(i,j,k)
       n_abs = nabs_surf(i,j,k)
       x_ghost = (i+lo(1)-1-.5)*dx
       y_ghost = (j+lo(2)-1-.5)*dy
       z_ghost = zc(k)
       do l=0,lmax
         yyy = y_ghost+l*(ny/(n_abs+small))*step
         zzz = z_ghost+l*(nz/(n_abs+small))*step
         xxx = x_ghost+l*(nx/(n_abs+small))*step
         if(trim(surface_type) == 'HeightMap') then 
	         inside = height_map(xxx,yyy,zzz,dxl,dyl,i,j,dzc(k),n,surf_height,lo,hi)
	     elseif(trim(surface_type) == 'Wall') then
             inside = wall(xxx,yyy,zzz,i,j,dzc(k))
         elseif(trim(surface_type) == 'Sphere') then
             inside = sphere(xxx,yyy,zzz,i,j,dzc(k),lo,hi)
         endif
         if (.not.inside) then
           x_intersect(i,j,k) = xxx
           y_intersect(i,j,k) = yyy
           z_intersect(i,j,k) = zzz
           confirmation = .true.
          exit
         endif
      enddo
      if (.not.confirmation) then
       print*, '--------------------------------------------------------------------------'
       print*,'Error in detecting intersect point at  i , j ,k=' &
              ,i,j,k,'at processor ',myid, 'where the normal vector components are ', &
               nx,ny,nz,n_abs
       print*, '--------------------------------------------------------------------------'
      endif

       distance_ghost_intersect =  sqrt((x_intersect(i,j,k)-x_ghost)**2 + &
                                        (y_intersect(i,j,k)-y_ghost)**2 + &
                                        (z_intersect(i,j,k)-z_ghost)**2)

       if (distance_ghost_intersect.gt.sqrt(dx**2 + dy**2 + dzc(k)**2)) then
        print*, '--------------------------------------------------------------------------'
        print*, ' Error in detecting intersect point  processor   ', &
        myid,' : check IBM.f90 - distance_ghost_intesect is ',distance_ghost_intersect, &
       'where the normal vector components are ', nx,ny,nz,n_abs
        print*, '--------------------------------------------------------------------------'
       endif

    endif
   enddo
  enddo
enddo

return
end subroutine intersect

Subroutine mirrorpoints(n,lo,hi,&
                        nx_surf,ny_surf,nz_surf,nabs_surf, &
                        x_intersect,y_intersect,z_intersect, &
                        x_mirror,y_mirror,z_mirror, &
                        deltan, &
                        x_IP1,y_IP1,z_IP1,x_IP2,y_IP2,z_IP2, &
                        zc,dzc)

implicit none
integer,dimension(3),intent(in) :: n,lo,hi
real(rp),dimension(0:n(3)+1),intent(in) :: zc,dzc
real(rp),dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6),intent(in) :: nx_surf
real(rp),dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6),intent(in) :: ny_surf
real(rp),dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6),intent(in) :: nz_surf
real(rp),dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6),intent(in) :: nabs_surf
real(rp),dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6),intent(in) :: x_intersect
real(rp),dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6),intent(in) :: y_intersect
real(rp),dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6),intent(in) :: z_intersect
real(rp),dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6),intent(out):: x_mirror
real(rp),dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6),intent(out):: y_mirror
real(rp),dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6),intent(out):: z_mirror
real(rp),dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6),intent(out):: deltan
real(rp),dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6),intent(out):: x_IP1
real(rp),dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6),intent(out):: y_IP1
real(rp),dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6),intent(out):: z_IP1
real(rp),dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6),intent(out):: x_IP2
real(rp),dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6),intent(out):: y_IP2
real(rp),dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6),intent(out):: z_IP2
integer::i,j,k,l,m,nn
real(rp):: nx, ny, nz, n_abs
real(rp):: step,step_2
real(rp):: xxx,yyy,zzz
real(rp):: cell_start_y,cell_end_y,cell_start_z,cell_end_z
real(rp)::distance_ghost_intersect

!***************************************************************************************
! Part one: cell-centered
!***************************************************************************************
do k=1,int(solid_height_ratio*(n(3)+1))
if (myid.eq.0) print*, 'Calculating Mirror Points at k = ', k
 step   = sqrt(dx**2 + dy**2 + dzc(k)**2)
 step_2 = 0.5_rp*sqrt(dx**2 + dy**2 + dzc(k)**2)
 do j= 1,n(2)
   do i= 1,n(1)
    if (nabs_surf(i,j,k).gt.small) then

       yyy = (i+lo(1)-1-.5)*dx
       xxx = (j+lo(2)-1-.5)*dy
       zzz = zc(k)

       nx =   nx_surf(i,j,k)
       ny =   ny_surf(i,j,k)
       nz =   nz_surf(i,j,k)
       n_abs = nabs_surf(i,j,k)
       distance_ghost_intersect = sqrt((x_intersect(i,j,k)-xxx)**2 + &
                                       (y_intersect(i,j,k)-yyy)**2 + &
                                       (z_intersect(i,j,k)-zzz)**2)


       if  (distance_ghost_intersect.gt.sqrt(dx**2 + dy**2 + dzc(k)**2)) then
       print*, '--------------------------------------------------------------------------'
           print*, ' Error: in mirror point detection at cell-center at processor ', &
            myid,' : check IBM.f90 - distance_ghost_intesect is ',distance_ghost_intersect, &
            'where the normal vector components are ', nx,ny,nz,n_abs
       print*, '--------------------------------------------------------------------------'
       endif

       x_mirror(i,j,k) = x_intersect(i,j,k)+(nx/(n_abs+small))*step
       y_mirror(i,j,k) = y_intersect(i,j,k)+(ny/(n_abs+small))*step !change for 3D
       z_mirror(i,j,k) = z_intersect(i,j,k)+(nz/(n_abs+small))*step !change for 3D

       deltan(i,j,k)   = distance_ghost_intersect

       x_IP1(i,j,k)    = x_mirror(i,j,k)+(nx/(n_abs+small))*step_2
       y_IP1(i,j,k)    = y_mirror(i,j,k)+(ny/(n_abs+small))*step_2
       z_IP1(i,j,k)    = z_mirror(i,j,k)+(nz/(n_abs+small))*step_2

       x_IP2(i,j,k)    = x_IP1(i,j,k)+(nx/(n_abs+small))*step_2
       y_IP2(i,j,k)    = y_IP1(i,j,k)+(ny/(n_abs+small))*step_2
       z_IP2(i,j,k)    = z_IP1(i,j,k)+(nz/(n_abs+small))*step_2

    endif
   enddo
  enddo
enddo



return
end subroutine mirrorpoints

Subroutine mirrorpoints_ijk(n,lo,hi, &
                            nabs_surf,x_mirror,y_mirror,z_mirror, &
                            x_IP1,y_IP1,z_IP1,x_IP2,y_IP2,z_IP2, &
                            i_mirror,j_mirror,k_mirror, &
                            i_IP1,j_IP1,k_IP1,i_IP2,j_IP2,k_IP2, &
                            zf,dzc)

implicit none
integer,dimension(3),intent(in) :: n,lo,hi
real(rp),dimension(0:n(3)+1),intent(in) :: zf,dzc
real(rp),dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6),intent(in):: nabs_surf
real(rp),dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6),intent(in):: x_mirror
real(rp),dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6),intent(in):: y_mirror
real(rp),dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6),intent(in):: z_mirror
real(rp),dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6),intent(in):: x_IP1
real(rp),dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6),intent(in):: y_IP1
real(rp),dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6),intent(in):: z_IP1
real(rp),dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6),intent(in):: x_IP2
real(rp),dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6),intent(in):: y_IP2
real(rp),dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6),intent(in):: z_IP2


integer,dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6),intent(out):: i_mirror
integer,dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6),intent(out):: j_mirror
integer,dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6),intent(out):: k_mirror
integer,dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6),intent(out):: i_IP1
integer,dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6),intent(out):: j_IP1
integer,dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6),intent(out):: k_IP1
integer,dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6),intent(out):: i_IP2
integer,dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6),intent(out):: j_IP2
integer,dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6),intent(out):: k_IP2

integer::i,j,k,l,m,nn
real(rp):: nx, ny, nz, n_abs
real(rp):: xxx,yyy,zzz
real(rp)::cell_start_x,cell_end_x,cell_start_y,cell_end_y,cell_start_z,cell_end_z
real(rp)::distance_ghost_intesect

do k=1,int(solid_height_ratio*(n(3)+1))
if (myid.eq.0) print*, 'Calculating mirror points indices at k = ', k
 do j= 1,n(2)
   do i= 1,n(1)
   !***************************************************************************************
   ! Part one: cell-centered
   !***************************************************************************************
    if (nabs_surf(i,j,k).gt.small) then
    ! Mirror points
       xxx = x_mirror(i,j,k)
       yyy = y_mirror(i,j,k)
       zzz = z_mirror(i,j,k)
       do l = 1,n(3)
         do m = -5,n(2)+6
           do nn = -5,n(1)+6
             cell_start_x = (i+lo(1)-1-1.)*dx
             cell_end_x   = (i+lo(1)-1-.0)*dx 
             cell_start_y = (j+lo(2)-1-1.)*dy
             cell_end_y   = (j+lo(2)-1-.0)*dy
             cell_start_z = zf(l-1)
             cell_end_z   = zf(l)
             if ((yyy.ge.cell_start_y).and.(yyy.lt.cell_end_y).and.&
                 (zzz.ge.cell_start_z).and.(zzz.lt.cell_end_z).and.&
                 (xxx.ge.cell_start_x).and.(xxx.lt.cell_end_x)) then
                  i_mirror(i,j,k) = nn
                  j_mirror(i,j,k) = m
                  k_mirror(i,j,k) = l
              exit


             endif
            enddo
          enddo
       enddo
       if ((i_mirror(i,j,k).eq.i).and.(j_mirror(i,j,k).eq.j).and.(k_mirror(i,j,k).eq.k)) then
             print*, 'Error: Ghost and mirror point are the same'
       endif
       if ((i_mirror(i,j,k).eq.-1000).or.(j_mirror(i,j,k).eq.-1000).or.(k_mirror(i,j,k).eq.-1000)) then
       print*, '--------------------------------------------------------------------------'
         print*,'Error: no grid point detected for mirror point at center i = ', &
         i, ' j= ',j,' k= ',k, 'at processor ',myid
       print*, '--------------------------------------------------------------------------'
      endif

    ! First interpolation points
       xxx = x_IP1(i,j,k)
       yyy = y_IP1(i,j,k)
       zzz = z_IP1(i,j,k)
       do l = 1,n(3)
         do m = -5,n(2)+6
           do nn = -5,n(1)+6
             cell_start_x = (i+lo(1)-1-1.)*dx
             cell_end_x   = (i+lo(1)-1-.0)*dx 
             cell_start_y = (j+lo(2)-1-1.)*dy
             cell_end_y   = (j+lo(2)-1-.0)*dy
             cell_start_z = zf(l-1)
             cell_end_z   = zf(l)
             if ((yyy.ge.cell_start_y).and.(yyy.lt.cell_end_y).and.&
                 (zzz.ge.cell_start_z).and.(zzz.lt.cell_end_z).and.&
                 (xxx.ge.cell_start_x).and.(xxx.lt.cell_end_x)) then
                  i_IP1(i,j,k) = nn
                  j_IP1(i,j,k) = m
                  k_IP1(i,j,k) = l
              exit
             endif
            enddo
          enddo
       enddo


       if ((i_IP1(i,j,k).eq.-1000).or.(j_IP1(i,j,k).eq.-1000).or.(k_IP1(i,j,k).eq.-1000)) then
       print*, '--------------------------------------------------------------------------'
         print*,'Error: no grid point detected for the first interpolation point at center i = '&
         ,i, ' j= ',j,' k= ',k, 'at processor ',myid
       print*, '--------------------------------------------------------------------------'
      endif

    ! Second interpolation points
       xxx = x_IP2(i,j,k)
       yyy = y_IP2(i,j,k)
       zzz = z_IP2(i,j,k)
       do l = 1,n(3)
         do m = -5,n(2)+6
           do nn = -5,n(1)+6
             cell_start_x = (i+lo(1)-1-1.)*dx
             cell_end_x   = (i+lo(1)-1-.0)*dx 
             cell_start_y = (j+lo(2)-1-1.)*dy
             cell_end_y   = (j+lo(2)-1-.0)*dy
             cell_start_z = zf(l-1)
             cell_end_z   = zf(l)
             if ((yyy.ge.cell_start_y).and.(yyy.lt.cell_end_y).and.&
                 (zzz.ge.cell_start_z).and.(zzz.lt.cell_end_z).and.&
                 (xxx.ge.cell_start_x).and.(xxx.lt.cell_end_x)) then
                  i_IP2(i,j,k) = nn
                  j_IP2(i,j,k) = m
                  k_IP2(i,j,k) = l
              exit
             endif

            enddo
          enddo
       enddo


       if ((i_IP2(i,j,k).eq.-1000).or.(j_IP2(i,j,k).eq.-1000).or.(k_IP2(i,j,k).eq.-1000)) then
       print*, '--------------------------------------------------------------------------'
         print*,'Error: no grid point detected for the second interpolation  point at center i = ', &
         i, ' j= ',j,' k= ',k, 'at processor ',myid
       print*, '--------------------------------------------------------------------------'
      endif
    endif
  enddo
 enddo
enddo
return
end Subroutine mirrorpoints_ijk


Subroutine InterpolationWeights(n,lo,hi, &
                                nabs_surf,Level_set, &
                                x_IP1,y_IP1,z_IP1, &
                                x_IP2,y_IP2,z_IP2, &
                                i_IP1,j_IP1,k_IP1, &
                                i_IP2,j_IP2,k_IP2, &
                                WP1,WP2, &
                                zc)

implicit none
integer,dimension(3),intent(in) :: n,lo,hi
real(rp),dimension(0:n(3)+1),intent(in) :: zc
real(rp),dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6),intent(in):: nabs_surf
integer, dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6),intent(in):: Level_set
real(rp),dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6),intent(in):: x_IP1
real(rp),dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6),intent(in):: y_IP1
real(rp),dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6),intent(in):: z_IP1
real(rp),dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6),intent(in):: x_IP2
real(rp),dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6),intent(in):: y_IP2
real(rp),dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6),intent(in):: z_IP2
integer ,dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6),intent(in):: i_IP1
integer ,dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6),intent(in):: j_IP1
integer ,dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6),intent(in):: k_IP1
integer ,dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6),intent(in):: i_IP2
integer ,dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6),intent(in):: j_IP2
integer ,dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6),intent(in):: k_IP2

real(rp),dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6,1:7),intent(out):: WP1
real(rp),dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6,1:7),intent(out):: WP2

real(rp):: x1(7),y1(7),z1(7),x2(7),y2(7),z2(7)
real(rp):: h1(7),h2(7)
integer :: i_p_1(7),j_p_1(7),k_p_1(7),i_p_2(7),j_p_2(7),k_p_2(7)
integer :: ii,i,j,k,ii1,ii2,jj1,jj2,kk1,kk2
real(rp):: xx1,xx2,yy1,yy2,zz1,zz2
logical :: cond11,cond21,cond12,cond22
real(rp):: contribution1(7),contribution2(7)

do k=1,int(solid_height_ratio*(n(3)+1))
if (myid.eq.0) print*, 'Calculating Interpolation Weights at k = ', k
 do j= 1,n(2)
   do i= 1,n(1)
   !***************************************************************************************
   ! Part one: cell-centered
   !***************************************************************************************
     !Initialization
     xx1               = 0.
     xx2               = 0.
     yy1               = 0.
     yy2               = 0.
     zz1               = 0.
     zz2               = 0.
     cond11            = .false.
     cond21            = .false.
     cond12            = .false.
     cond22            = .false.
     contribution1(:)  =  0.
     contribution2(:)  =  0.
     x1(:)             =  0.
     y1(:)             =  0.
     z1(:)             =  0.
     x2(:)             =  0.
     y2(:)             =  0.
     z2(:)             =  0.
     h1(:)             =  0.
     h2(:)             =  0.
     i_p_1(:)          =  0
     j_p_1(:)          =  0
     k_p_1(:)          =  0
     i_p_2(:)          =  0
     j_p_2(:)          =  0
     k_p_2(:)          =  0
     if (nabs_surf(i,j,k).gt.small) then
         xx1  =  x_IP1(i,j,k)
         yy1  =  y_IP1(i,j,k)
         zz1  =  z_IP1(i,j,k)
         ii1  =  i_IP1(i,j,k)
         jj1  =  j_IP1(i,j,k)
         kk1  =  k_IP1(i,j,k)
         
         xx2  =  x_IP2(i,j,k)
         yy2  =  y_IP2(i,j,k)
         zz2  =  z_IP2(i,j,k)
         ii2  =  i_IP2(i,j,k)
         jj2  =  j_IP2(i,j,k)
         kk2  =  k_IP2(i,j,k)
         
         i_p_1(1) = ii1
         j_p_1(1) = jj1+1
         k_p_1(1) = kk1
         
         i_p_2(1) = ii2
         j_p_2(1) = jj2+1
         k_p_2(1) = kk2
         
         i_p_1(2) = ii1
         j_p_1(2) = jj1
         k_p_1(2) = kk1+1
         
         i_p_2(2) = ii2
         j_p_2(2) = jj2
         k_p_2(2) = kk2+1
         
         i_p_1(3) = ii1
         j_p_1(3) = jj1-1
         k_p_1(3) = kk1
         
         i_p_2(3) = ii2
         j_p_2(3) = jj2-1
         k_p_2(3) = kk2
         
         i_p_1(4) = ii1
         j_p_1(4) = jj1
         k_p_1(4) = kk1-1
         
         i_p_2(4) = ii2
         j_p_2(4) = jj2
         k_p_2(4) = kk2-1
         
         i_p_1(5) = ii1
         j_p_1(5) = jj1
         k_p_1(5) = kk1
         
         i_p_2(5) = ii2
         j_p_2(5) = jj2
         k_p_2(5) = kk2
         
         i_p_1(6) = ii1-1
         j_p_1(6) = jj1
         k_p_1(6) = kk1
         
         i_p_2(6) = ii2-1
         j_p_2(6) = jj2
         k_p_2(6) = kk2
         
         i_p_1(7) = ii1+1
         j_p_1(7) = jj1
         k_p_1(7) = kk1
         
         i_p_2(7) = ii2+1
         j_p_2(7) = jj2
         k_p_2(7) = kk2

        do ii = 1,7
           cond11 = (Level_set(i_p_1(ii),j_p_1(ii),k_p_1(ii)).eq.1)
           cond21 = (nabs_surf(i_p_1(ii),j_p_1(ii),k_p_1(ii)).eq.0.)
           cond12 = (Level_set(i_p_2(ii),j_p_2(ii),k_p_2(ii)).eq.1)
           cond22 = (nabs_surf(i_p_2(ii),j_p_2(ii),k_p_2(ii)).eq.0.)

           if (cond11.and.cond21) contribution1(ii) = 1.
           if (cond12.and.cond22) contribution2(ii) = 1.

           x1(ii) = (i_p_1(ii)+lo(1)-1-.5)*dx 
           y1(ii) = (j_p_1(ii)+lo(2)-1-.5)*dy
           z1(ii) = zc(k_p_1(ii)           )
           x2(ii) = (i_p_2(ii)+lo(1)-1-.5)*dx 
           y2(ii) = (j_p_2(ii)+lo(2)-1-.5)*dy
           z2(ii) = zc(k_p_2(ii)           )
        enddo

        do ii = 1,7
         h1(ii) = sqrt( (x1(ii)-xx1)**2 + (y1(ii)-yy1)**2 + (z1(ii)-zz1)**2 )
         h2(ii) = sqrt( (x2(ii)-xx2)**2 + (y2(ii)-yy2)**2 + (z2(ii)-zz2)**2 )
        enddo

        do ii = 1,7
         WP1(i,j,k,ii)  = (1./(h1(ii)*h1(ii)))*contribution1(ii)
         WP2(i,j,k,ii)  = (1./(h2(ii)*h2(ii)))*contribution2(ii)
        enddo

        !-------- Exceptional cases ---------------------
        do ii = 1,7
         if ((h1(ii).lt.(1.e-8_rp)).and.(contribution1(ii).eq.1._rp)) then
             WP1(i,j,k,1)  = 0.
             WP1(i,j,k,2)  = 0.
             WP1(i,j,k,3)  = 0.
             WP1(i,j,k,4)  = 0.
             WP1(i,j,k,5)  = 0.
             WP1(i,j,k,6)  = 0.
             WP1(i,j,k,7)  = 0.
             WP1(i,j,k,ii) = 1.
         endif
         if ((h2(ii).lt.(1.e-8_rp)).and.(contribution2(ii).eq.1._rp)) then
             WP2(i,j,k,1)  = 0.
             WP2(i,j,k,2)  = 0.
             WP2(i,j,k,3)  = 0.
             WP2(i,j,k,4)  = 0.
             WP2(i,j,k,5)  = 0.
             WP2(i,j,k,6)  = 0.
             WP2(i,j,k,7)  = 0.
             WP2(i,j,k,ii) = 1.
         endif
       enddo
       !-------------------------------------------
     endif

   enddo
  enddo
enddo

return
end subroutine InterpolationWeights

subroutine interpolation_mirror(n,AA,iii,jjj,kkk, &
                                i_IP1,j_IP1,k_IP1, &
                                i_IP2,j_IP2,k_IP2, &
                                WP1,WP2,BB)

implicit none
integer, dimension(3), intent(in)  :: n
real(rp),dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6,1:7),intent(in):: WP1,WP2
real(rp),dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6),intent(in)    :: AA
integer, dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6),intent(in)    :: i_IP1,j_IP1,k_IP1,i_IP2,j_IP2,k_IP2
integer, intent(in)  :: iii,jjj,kkk
real(rp),intent(out) :: BB
real(rp),dimension(7) ::WW1,WW2,B1,B2
real(rp) :: q1,q2,B_1,B_2
integer  :: ii1,jj1,kk1,ii2,jj2,kk2,l

!initialization
q1       = 0.
q2       = 0.
WW1(:)   = 0.
WW2(:)   = 0.
B1(:)    = 0.
B2(:)    = 0.
B_1      = 0.
B_2      = 0.
BB       = 0.

ii1    = i_IP1(iii,jjj,kkk) ! Coordinates for cell of interpolation point 1 (IP1)
jj1    = j_IP1(iii,jjj,kkk)
kk1    = k_IP1(iii,jjj,kkk)
ii2    = i_IP2(iii,jjj,kkk) ! Coordinates for cell of interpolation point 2  (IP2)
jj2    = j_IP2(iii,jjj,kkk)
kk2    = k_IP2(iii,jjj,kkk)

WW1(1) = WP1(iii,jjj,kkk,1) ! Weights for interpolating to IP1
WW1(2) = WP1(iii,jjj,kkk,2)
WW1(3) = WP1(iii,jjj,kkk,3)
WW1(4) = WP1(iii,jjj,kkk,4)
WW1(5) = WP1(iii,jjj,kkk,5)
WW1(6) = WP1(iii,jjj,kkk,6)
WW1(7) = WP1(iii,jjj,kkk,7)
WW2(1) = WP2(iii,jjj,kkk,1) ! Weights for interpolating to IP2
WW2(2) = WP2(iii,jjj,kkk,2)
WW2(3) = WP2(iii,jjj,kkk,3)
WW2(4) = WP2(iii,jjj,kkk,4)
WW2(5) = WP2(iii,jjj,kkk,5)
WW2(6) = WP2(iii,jjj,kkk,6)
WW2(7) = WP2(iii,jjj,kkk,7)
B1(1)   = AA(ii1,jj1+1,kk1)  ! Values for interpolating to IP1
B1(2)   = AA(ii1,jj1,kk1+1)
B1(3)   = AA(ii1,jj1-1,kk1)
B1(4)   = AA(ii1,jj1,kk1-1)
B1(5)   = AA(ii1,jj1,kk1)
B1(6)   = AA(ii1-1,jj1,kk1)
B1(7)   = AA(ii1+1,jj1,kk1)
B2(1)   = AA(ii2,jj2+1,kk2)  ! Values for interpolating to IP2
B2(2)   = AA(ii2,jj2,kk2+1)
B2(3)   = AA(ii2,jj2-1,kk2)
B2(4)   = AA(ii2,jj2,kk2-1)
B2(5)   = AA(ii2,jj2,kk2)
B2(6)   = AA(ii2-1,jj2,kk2)
B2(7)   = AA(ii2+1,jj2,kk2)

do l = 1,7 ! Sum weights
   q1 =  max(WW1(l) + q1, small)
   q2 =  max(WW2(l) + q2, small)
enddo

do l = 1,7 ! Compute value at IP1 and IP2
   B_1 = (1./q1)*(WW1(l)*B1(l))+B_1
   B_2 = (1./q2)*(WW2(l)*B2(l))+B_2
enddo
BB = 2.*B_1-B_2 ! Value at mirror point

end subroutine interpolation_mirror

Subroutine interpolation_2D_velocity(n,UU,VV,WW,iii,jjj,kkk, &
                                     i_IP1,j_IP1,k_IP1, &
                                     i_IP2,j_IP2,k_IP2, &
                                     WP1,WP2,U_m,V_m,W_m)

implicit none
integer, dimension(3), intent(in)  :: n
integer, intent(in) :: iii,jjj,kkk
integer, dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6),intent(in)    :: i_IP1
integer, dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6),intent(in)    :: j_IP1
integer, dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6),intent(in)    :: k_IP1
integer, dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6),intent(in)    :: i_IP2
integer, dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6),intent(in)    :: j_IP2
integer, dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6),intent(in)    :: k_IP2
real(rp),dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6,1:7),intent(in):: WP1
real(rp),dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6,1:7),intent(in):: WP2

real(rp), intent(out) ::U_m, V_m, W_m
real(rp), dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6), intent(in):: UU ,VV ,WW
real(rp) :: q1,q2,WW1(7),WW2(7),U1(7),U2(7),V1(7),V2(7),W1(7),W2(7)
real(rp)::  U_1,U_2,V_1,V_2,W_1,W_2
integer:: ii1,jj1,kk1,ii2,jj2,kk2,l

!initialization
q1      = 0.
q2      = 0.
WW1(:)  = 0.
WW2(:)  = 0.
U1(:)   = 0.
U2(:)   = 0.
V1(:)   = 0.
V2(:)   = 0.
W1(:)   = 0.
W2(:)   = 0.
U_1     = 0.
U_2     = 0.
V_1     = 0.
V_2     = 0.
W_1     = 0.
W_2     = 0.

ii1     = i_IP1(iii,jjj,kkk)
jj1     = j_IP1(iii,jjj,kkk)
kk1     = k_IP1(iii,jjj,kkk)
ii2     = i_IP2(iii,jjj,kkk)
jj2     = j_IP2(iii,jjj,kkk)
kk2     = k_IP2(iii,jjj,kkk)

WW1(1) = WP1(iii,jjj,kkk,1)
WW1(2) = WP1(iii,jjj,kkk,2)
WW1(3) = WP1(iii,jjj,kkk,3)
WW1(4) = WP1(iii,jjj,kkk,4)
WW1(5) = WP1(iii,jjj,kkk,5)
WW1(6) = WP1(iii,jjj,kkk,6)
WW1(7) = WP1(iii,jjj,kkk,7)
WW2(1) = WP2(iii,jjj,kkk,1)
WW2(2) = WP2(iii,jjj,kkk,2)
WW2(3) = WP2(iii,jjj,kkk,3)
WW2(4) = WP2(iii,jjj,kkk,4)
WW2(5) = WP2(iii,jjj,kkk,5)
WW2(6) = WP2(iii,jjj,kkk,6)
WW2(7) = WP2(iii,jjj,kkk,7)

U1(1)    = UU(ii1,jj1+1,kk1)
U1(2)    = UU(ii1,jj1,  kk1+1)
U1(3)    = UU(ii1,jj1-1,kk1)
U1(4)    = UU(ii1,jj1,  kk1-1)
U1(5)    = UU(ii1,jj1,  kk1)
U1(6)    = UU(ii1-1,jj1,  kk1)
U1(7)    = UU(ii1+1,jj1,  kk1)
U2(1)    = UU(ii2,jj2+1,kk2)
U2(2)    = UU(ii2,jj2,  kk2+1)
U2(3)    = UU(ii2,jj2-1,kk2)
U2(4)    = UU(ii2,jj2,  kk2-1)
U2(5)    = UU(ii2,jj2,  kk2)
U2(6)    = UU(ii2-1,jj2,  kk2)
U2(7)    = UU(ii2+1,jj2,  kk2)

V1(1)    = VV(ii1,jj1+1,kk1)
V1(2)    = VV(ii1,jj1,  kk1+1)
V1(3)    = VV(ii1,jj1-1,kk1)
V1(4)    = VV(ii1,jj1,  kk1-1)
V1(5)    = VV(ii1,jj1,  kk1)
V1(6)    = VV(ii1-1,jj1,  kk1)
V1(7)    = VV(ii1+1,jj1,  kk1)
V2(1)    = VV(ii2,jj2+1,kk2)
V2(2)    = VV(ii2,jj2,  kk2+1)
V2(3)    = VV(ii2,jj2-1,kk2)
V2(4)    = VV(ii2,jj2,  kk2-1)
V2(5)    = VV(ii2,jj2,  kk2)
V2(6)    = VV(ii2-1,jj2,  kk2)
V2(7)    = VV(ii2+1,jj2,  kk2)

W1(1)    = WW(ii1,jj1+1,kk1)
W1(2)    = WW(ii1,jj1,  kk1+1)
W1(3)    = WW(ii1,jj1-1,kk1)
W1(4)    = WW(ii1,jj1,  kk1-1)
W1(5)    = WW(ii1,jj1,  kk1)
W1(6)    = WW(ii1-1,jj1,  kk1)
W1(7)    = WW(ii1+1,jj1,  kk1)
W2(1)    = WW(ii2,jj2+1,kk2)
W2(2)    = WW(ii2,jj2,  kk2+1)
W2(3)    = WW(ii2,jj2-1,kk2)
W2(4)    = WW(ii2,jj2,  kk2-1)
W2(5)    = WW(ii2,jj2,  kk2)
W2(6)    = WW(ii2-1,jj2,  kk2)
W2(7)    = WW(ii2+1,jj2,  kk2)

do l = 1,7
   q1 =  max(WW1(l) + q1, small)
   q2 =  max(WW2(l) + q2, small)
enddo

do l = 1,7
 U_1 = (1./q1)*(WW1(l)*U1(l))+U_1
 U_2 = (1./q2)*(WW2(l)*U2(l))+U_2
 V_1 = (1./q1)*(WW1(l)*V1(l))+V_1
 V_2 = (1./q2)*(WW2(l)*V2(l))+V_2
 W_1 = (1./q1)*(WW1(l)*W1(l))+W_1
 W_2 = (1./q2)*(WW2(l)*W2(l))+W_2
enddo

u_m = 2.*U_1-U_2
V_m = 2.*V_1-V_2
W_m = 2.*W_1-W_2

return
end subroutine interpolation_2D_velocity

Subroutine interpolation_dphi(n,PFM_phi,dzc, &
                              iii,jjj,kkk,  &
                              i_IP1,j_IP1,k_IP1,  &
                              i_IP2,j_IP2,k_IP2,  &
                              WP1,WP2,  &
                              dphidx,dphidy,dphidz)
implicit none
integer, dimension(3), intent(in)  :: n
integer, intent(in) :: iii,jjj,kkk
real(rp),dimension(0:n(3)+1),intent(in) :: dzc
real(rp),dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6):: PFM_phi
integer, dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6),intent(in)      ::i_IP1
integer, dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6),intent(in)      ::j_IP1
integer, dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6),intent(in)      ::k_IP1
integer, dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6),intent(in)      ::i_IP2
integer, dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6),intent(in)      ::j_IP2
integer, dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6),intent(in)      ::k_IP2
real(rp),dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6,1:7),intent(in) :: WP1
real(rp),dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6,1:7),intent(in) :: WP2
real(rp),intent(out):: dphidx,dphidy,dphidz
integer  ::ii,l
real(rp) :: q1,q2,WW1(7),WW2(7)!,B1(5),B2(5),B_1,B_2
integer  :: ii1,jj1,kk1,ii2,jj2,kk2
real(rp) :: phi_ip1(7),phi_im1(7), phi_jp1(7),phi_jm1(7),phi_kp1(7),phi_km1(7)
real(rp) :: phi_ip2(7),phi_im2(7),phi_jp2(7),phi_jm2(7),phi_kp2(7),phi_km2(7)
real(rp) :: phi_x1(7),phi_y1(7),phi_z1(7),phi_x2(7),phi_y2(7),phi_z2(7)
real(rp) :: dphidx1,dphidy1,dphidz1,dphidx2,dphidy2,dphidz2,dzi

!initialization
q1          = 0.
q2          = 0.
WW1(:)      = 0.
WW2(:)      = 0.
phi_ip1(:)  = 0.
phi_im1(:)  = 0.
phi_jp1(:)  = 0.
phi_jm1(:)  = 0.
phi_kp1(:)  = 0.
phi_km1(:)  = 0.
phi_ip2(:)  = 0.
phi_im2(:)  = 0.
phi_jp2(:)  = 0.
phi_jm2(:)  = 0.
phi_kp2(:)  = 0.
phi_km2(:)  = 0.
phi_x1(:)   = 0.
phi_y1(:)   = 0.
phi_z1(:)   = 0.
phi_x2(:)   = 0.
phi_y2(:)   = 0.
phi_z2(:)   = 0.
dphidx1     = 0.
dphidy1     = 0.
dphidz1     = 0.
dphidx2     = 0.
dphidy2     = 0.
dphidz2     = 0.

  ii1     = i_IP1(iii,jjj,kkk)
  jj1     = j_IP1(iii,jjj,kkk)
  kk1     = k_IP1(iii,jjj,kkk)
  ii2     = i_IP2(iii,jjj,kkk)
  jj2     = j_IP2(iii,jjj,kkk)
  kk2     = k_IP2(iii,jjj,kkk)

  WW1(1) = WP1(iii,jjj,kkk,1) ! Weights for interpolating to IP1
  WW1(2) = WP1(iii,jjj,kkk,2)
  WW1(3) = WP1(iii,jjj,kkk,3)
  WW1(4) = WP1(iii,jjj,kkk,4)
  WW1(5) = WP1(iii,jjj,kkk,5)
  WW1(6) = WP1(iii,jjj,kkk,6)
  WW1(7) = WP1(iii,jjj,kkk,7)
  WW2(1) = WP2(iii,jjj,kkk,1) ! Weights for interpolating to IP2
  WW2(2) = WP2(iii,jjj,kkk,2)
  WW2(3) = WP2(iii,jjj,kkk,3)
  WW2(4) = WP2(iii,jjj,kkk,4)
  WW2(5) = WP2(iii,jjj,kkk,5)
  WW2(6) = WP2(iii,jjj,kkk,6)
  WW2(7) = WP2(iii,jjj,kkk,7)

  ! Cell edge values for AP for IP1
  phi_ip1(1) = 0.*(PFM_phi(ii1,jj1+1,kk1)+PFM_phi(ii1+1,jj1+1,kk1))
  phi_im1(1) = 0.*(PFM_phi(ii1,jj1+1,kk1)+PFM_phi(ii1-1,jj1+1,kk1))
  phi_jp1(1) = 0.*(PFM_phi(ii1,jj1+1,kk1)+PFM_phi(ii1,jj1+2,kk1))
  phi_jm1(1) = 0.*(PFM_phi(ii1,jj1+1,kk1)+PFM_phi(ii1,jj1,kk1))
  phi_kp1(1) = 0.*(PFM_phi(ii1,jj1+1,kk1)+PFM_phi(ii1,jj1+1,kk1+1))
  phi_km1(1) = 0.*(PFM_phi(ii1,jj1+1,kk1)+PFM_phi(ii1,jj1+1,kk1-1))


  phi_ip1(2) = 0.*(PFM_phi(ii1,jj1,kk1+1)+PFM_phi(ii1+1,jj1,kk1+1))
  phi_im1(2) = 0.*(PFM_phi(ii1,jj1,kk1+1)+PFM_phi(ii1-1,jj1-1,kk1+1))
  phi_jp1(2) = 0.*(PFM_phi(ii1,jj1,kk1+1)+PFM_phi(ii1,jj1+1,kk1+1))
  phi_jm1(2) = 0.*(PFM_phi(ii1,jj1,kk1+1)+PFM_phi(ii1,jj1-1,kk1+1))
  phi_kp1(2) = 0.*(PFM_phi(ii1,jj1,kk1+1)+PFM_phi(ii1,jj1,kk1+2))
  phi_km1(2) = 0.*(PFM_phi(ii1,jj1,kk1+1)+PFM_phi(ii1,jj1,kk1))


  phi_ip1(3) = 0.*(PFM_phi(ii1,jj1-1,kk1)+PFM_phi(ii1+1,jj1-1,kk1))
  phi_im1(3) = 0.*(PFM_phi(ii1,jj1-1,kk1)+PFM_phi(ii1-1,jj1-1,kk1))
  phi_jp1(3) = 0.*(PFM_phi(ii1,jj1-1,kk1)+PFM_phi(ii1,jj1,kk1))
  phi_jm1(3) = 0.*(PFM_phi(ii1,jj1-1,kk1)+PFM_phi(ii1,jj1-2,kk1))
  phi_kp1(3) = 0.*(PFM_phi(ii1,jj1-1,kk1)+PFM_phi(ii1,jj1-1,kk1+1))
  phi_km1(3) = 0.*(PFM_phi(ii1,jj1-1,kk1)+PFM_phi(ii1,jj1-1,kk1-1))

  phi_ip1(4) = 0.*(PFM_phi(ii1,jj1,kk1-1)+PFM_phi(ii1+1,jj1,kk1-1))
  phi_im1(4) = 0.*(PFM_phi(ii1,jj1,kk1-1)+PFM_phi(ii1-1,jj1-1,kk1-1))
  phi_jp1(4) = 0.*(PFM_phi(ii1,jj1,kk1-1)+PFM_phi(ii1,jj1+1,kk1-1))
  phi_jm1(4) = 0.*(PFM_phi(ii1,jj1,kk1-1)+PFM_phi(ii1,jj1-1,kk1-1))
  phi_kp1(4) = 0.*(PFM_phi(ii1,jj1,kk1-1)+PFM_phi(ii1,jj1,kk1))
  phi_km1(4) = 0.*(PFM_phi(ii1,jj1,kk1-1)+PFM_phi(ii1,jj1,kk1-2))


  phi_ip1(5) = 0.*(PFM_phi(ii1,jj1,kk1)+PFM_phi(ii1+1,jj1,kk1))
  phi_im1(5) = 0.*(PFM_phi(ii1,jj1,kk1)+PFM_phi(ii1-1,jj1,kk1))
  phi_jp1(5) = 0.*(PFM_phi(ii1,jj1,kk1)+PFM_phi(ii1,jj1+1,kk1))
  phi_jm1(5) = 0.*(PFM_phi(ii1,jj1,kk1)+PFM_phi(ii1,jj1-1,kk1))
  phi_kp1(5) = 0.*(PFM_phi(ii1,jj1,kk1)+PFM_phi(ii1,jj1,kk1+1))
  phi_km1(5) = 0.*(PFM_phi(ii1,jj1,kk1)+PFM_phi(ii1,jj1,kk1-1))


  phi_ip1(6) = 0.*(PFM_phi(ii1-1,jj1,kk1)+PFM_phi(ii1,jj1,kk1))
  phi_im1(6) = 0.*(PFM_phi(ii1-1,jj1,kk1)+PFM_phi(ii1-2,jj1,kk1))
  phi_jp1(6) = 0.*(PFM_phi(ii1-1,jj1,kk1)+PFM_phi(ii1-1,jj1+1,kk1))
  phi_jm1(6) = 0.*(PFM_phi(ii1-1,jj1,kk1)+PFM_phi(ii1-1,jj1-1,kk1))
  phi_kp1(6) = 0.*(PFM_phi(ii1-1,jj1,kk1)+PFM_phi(ii1-1,jj1,kk1+1))
  phi_km1(6) = 0.*(PFM_phi(ii1-1,jj1,kk1)+PFM_phi(ii1-1,jj1,kk1-1))


  phi_ip1(7) = 0.*(PFM_phi(ii1+1,jj1,kk1)+PFM_phi(ii1+2,jj1,kk1))
  phi_im1(7) = 0.*(PFM_phi(ii1+1,jj1,kk1)+PFM_phi(ii1,jj1,kk1))
  phi_jp1(7) = 0.*(PFM_phi(ii1+1,jj1,kk1)+PFM_phi(ii1+1,jj1+1,kk1))
  phi_jm1(7) = 0.*(PFM_phi(ii1+1,jj1,kk1)+PFM_phi(ii1+1,jj1-1,kk1))
  phi_kp1(7) = 0.*(PFM_phi(ii1+1,jj1,kk1)+PFM_phi(ii1+1,jj1,kk1+1))
  phi_km1(7) = 0.*(PFM_phi(ii1+1,jj1,kk1)+PFM_phi(ii1+1,jj1,kk1-1))

  ! Cell edge values for AP for IP2
  phi_ip2(1) = 0.*(PFM_phi(ii2,jj2+1,kk2)+PFM_phi(ii2+1,jj2+1,kk2))
  phi_im2(1) = 0.*(PFM_phi(ii2,jj2+1,kk2)+PFM_phi(ii2-1,jj2+1,kk2))
  phi_jp2(1) = 0.*(PFM_phi(ii2,jj2+1,kk2)+PFM_phi(ii2,jj2+2,kk2))
  phi_jm2(1) = 0.*(PFM_phi(ii2,jj2+1,kk2)+PFM_phi(ii2,jj2,kk2))
  phi_kp2(1) = 0.*(PFM_phi(ii2,jj2+1,kk2)+PFM_phi(ii2,jj2+1,kk2+1))
  phi_km2(1) = 0.*(PFM_phi(ii2,jj2+1,kk2)+PFM_phi(ii2,jj2+1,kk2-1))

  phi_ip2(2) = 0.*(PFM_phi(ii2,jj2,kk2+1)+PFM_phi(ii2+1,jj2,kk2+1))
  phi_im2(2) = 0.*(PFM_phi(ii2,jj2,kk2+1)+PFM_phi(ii2-1,jj2-1,kk2+1))
  phi_jp2(2) = 0.*(PFM_phi(ii2,jj2,kk2+1)+PFM_phi(ii2,jj2+1,kk2+1))
  phi_jm2(2) = 0.*(PFM_phi(ii2,jj2,kk2+1)+PFM_phi(ii2,jj2-1,kk2+1))
  phi_kp2(2) = 0.*(PFM_phi(ii2,jj2,kk2+1)+PFM_phi(ii2,jj2,kk2+2))
  phi_km2(2) = 0.*(PFM_phi(ii2,jj2,kk2+1)+PFM_phi(ii2,jj2,kk2))

  phi_ip2(3) = 0.*(PFM_phi(ii2,jj2-1,kk2)+PFM_phi(ii2+1,jj2-1,kk2))
  phi_im2(3) = 0.*(PFM_phi(ii2,jj2-1,kk2)+PFM_phi(ii2-1,jj2-1,kk2))
  phi_jp2(3) = 0.*(PFM_phi(ii2,jj2-1,kk2)+PFM_phi(ii2,jj2,kk2))
  phi_jm2(3) = 0.*(PFM_phi(ii2,jj2-1,kk2)+PFM_phi(ii2,jj2-2,kk2))
  phi_kp2(3) = 0.*(PFM_phi(ii2,jj2-1,kk2)+PFM_phi(ii2,jj2-1,kk2+1))
  phi_km2(3) = 0.*(PFM_phi(ii2,jj2-1,kk2)+PFM_phi(ii2,jj2-1,kk2-1))

  phi_ip2(4) = 0.*(PFM_phi(ii2,jj2,kk2-1)+PFM_phi(ii2+1,jj2,kk2-1))
  phi_im2(4) = 0.*(PFM_phi(ii2,jj2,kk2-1)+PFM_phi(ii2-1,jj2-1,kk2-1))
  phi_jp2(4) = 0.*(PFM_phi(ii2,jj2,kk2-1)+PFM_phi(ii2,jj2+1,kk2-1))
  phi_jm2(4) = 0.*(PFM_phi(ii2,jj2,kk2-1)+PFM_phi(ii2,jj2-1,kk2-1))
  phi_kp2(4) = 0.*(PFM_phi(ii2,jj2,kk2-1)+PFM_phi(ii2,jj2,kk2))
  phi_km2(4) = 0.*(PFM_phi(ii2,jj2,kk2-1)+PFM_phi(ii2,jj2,kk2-2))

  phi_ip2(5) = 0.*(PFM_phi(ii2,jj2,kk2)+PFM_phi(ii2+1,jj2,kk2))
  phi_im2(5) = 0.*(PFM_phi(ii2,jj2,kk2)+PFM_phi(ii2-1,jj2,kk2))
  phi_jp2(5) = 0.*(PFM_phi(ii2,jj2,kk2)+PFM_phi(ii2,jj2+1,kk2))
  phi_jm2(5) = 0.*(PFM_phi(ii2,jj2,kk2)+PFM_phi(ii2,jj2-1,kk2))
  phi_kp2(5) = 0.*(PFM_phi(ii2,jj2,kk2)+PFM_phi(ii2,jj2,kk2+1))
  phi_km2(5) = 0.*(PFM_phi(ii2,jj2,kk2)+PFM_phi(ii2,jj2,kk2-1))

  phi_ip2(6) = 0.*(PFM_phi(ii2-1,jj2,kk2)+PFM_phi(ii2,jj2,kk2))
  phi_im2(6) = 0.*(PFM_phi(ii2-1,jj2,kk2)+PFM_phi(ii2-2,jj2,kk2))
  phi_jp2(6) = 0.*(PFM_phi(ii2-1,jj2,kk2)+PFM_phi(ii2-1,jj2+1,kk2))
  phi_jm2(6) = 0.*(PFM_phi(ii2-1,jj2,kk2)+PFM_phi(ii2-1,jj2-1,kk2))
  phi_kp2(6) = 0.*(PFM_phi(ii2-1,jj2,kk2)+PFM_phi(ii2-1,jj2,kk2+1))
  phi_km2(6) = 0.*(PFM_phi(ii2-1,jj2,kk2)+PFM_phi(ii2-1,jj2,kk2-1))

  phi_ip2(7) = 0.*(PFM_phi(ii2+1,jj2,kk2)+PFM_phi(ii2+2,jj2,kk2))
  phi_im2(7) = 0.*(PFM_phi(ii2+1,jj2,kk2)+PFM_phi(ii2,jj2,kk2))
  phi_jp2(7) = 0.*(PFM_phi(ii2+1,jj2,kk2)+PFM_phi(ii2+1,jj2+1,kk2))
  phi_jm2(7) = 0.*(PFM_phi(ii2+1,jj2,kk2)+PFM_phi(ii2+1,jj2-1,kk2))
  phi_kp2(7) = 0.*(PFM_phi(ii2+1,jj2,kk2)+PFM_phi(ii2+1,jj2,kk2+1))
  phi_km2(7) = 0.*(PFM_phi(ii2+1,jj2,kk2)+PFM_phi(ii2+1,jj2,kk2-1))


dzi = 1./dzc(kkk)

! Derivatives at AP
  do ii=1,7
    phi_x1(ii) = (phi_ip1(ii)-phi_im1(ii))*dxi
    phi_y1(ii) = (phi_jp1(ii)-phi_jm1(ii))*dyi
    phi_z1(ii) = (phi_kp1(ii)-phi_km1(ii))*dzi
    phi_x2(ii) = (phi_ip2(ii)-phi_im2(ii))*dxi
    phi_y2(ii) = (phi_jp2(ii)-phi_jm2(ii))*dyi
    phi_z2(ii) = (phi_kp2(ii)-phi_km2(ii))*dzi
    q1 =  max(WW1(ii) + q1, small)
    q2 =  max(WW2(ii) + q2, small)
  enddo



! Derivatives at IP1 and IP2
  do l = 1,7
   dphidx1 = (1./q1)*(WW1(l)*phi_x1(l))+dphidx1
   dphidy1 = (1./q1)*(WW1(l)*phi_y1(l))+dphidy1
   dphidz1 = (1./q1)*(WW1(l)*phi_z1(l))+dphidz1
   dphidx2 = (1./q2)*(WW2(l)*phi_x2(l))+dphidx2
   dphidy2 = (1./q2)*(WW2(l)*phi_y2(l))+dphidy2
   dphidz2 = (1./q2)*(WW2(l)*phi_z2(l))+dphidz2
  enddo

! Extrapolate to mirror point
  dphidx = 2.*dphidx1-dphidx2
  dphidy = 2.*dphidy1-dphidy2
  dphidz = 2.*dphidz1-dphidz2

end subroutine interpolation_dphi
#endif
#endif
end module mod_IBM
