#if defined(_LPP_GPU)
module mod_lag_part
   !@cuf use cudafor
   use mpi
   ! use decomp_2d
   ! use decomp_2d_io
   use mod_common_mpi, only: myid, mydev, coord, comm_cart, ierr, &
                             cross_buff_id, cross_buff_id_rank, cross_buff_new, cross_buff_new_rank, &
                             xcOld_r_buf, ycOld_r_buf, zcOld_r_buf, ucOld_r_buf, vcOld_r_buf, wcOld_r_buf, uc_r_buf, vc_r_buf, wc_r_buf, uf_r_buf, vf_r_buf, wf_r_buf, xcOld_s_buf, ycOld_s_buf, zcOld_s_buf, ucOld_s_buf, vcOld_s_buf, wcOld_s_buf, uc_s_buf, vc_s_buf, wc_s_buf, uf_s_buf, vf_s_buf, wf_s_buf, ic_s_buf, jc_s_buf, kc_s_buf, ic_r_buf, jc_r_buf, kc_r_buf
   use mod_param,      only: lx,ly,lz,dx,dy,dxi,dyi,dims,pi,small,rho,visc, &
                             gacc,datadir
   use mod_types
   implicit none
   private
   public initLPP,SeedParticles,boundlpp,lppsweeps,ComputeSubDerivativeVel,outLPP, &
          StorePartOld,AveragePartSol

   integer :: istat
   logical :: LPP_initialized = .false.
   integer, dimension(:), allocatable, managed :: num_part

   type particle
      real(rp), managed :: xc,yc,zc,uc,vc,wc,duc,dvc,dwc,vol,sur, &
                                        xcOld,ycOld,zcOld,ucOld,vcOld,wcOld,uf,vf,wf
      integer, managed :: id,ic,jc,kc,dummyint
      ! Note: OPEN_MPI sometimes failed to communicate the last variable in 
      !       MPI_TYPE_STRUCT correctly, dummyint is included to go around 
      !       this bug in mpi
   end type particle 
   type (particle), dimension(:,:), allocatable, managed :: parts

   type particle_collect
      real(rp), managed :: xc,yc,zc,hfx,hfy,hfz,vol
      integer,  managed :: ic,jc,kc,dummyint
   end type particle_collect

   type (particle_collect), dimension(:,:), allocatable, managed :: parts_collect
   type (particle_collect), dimension(:),   allocatable, managed :: parts_collect_rank

   integer, dimension(:,:), allocatable, managed :: parts_cross_id
   integer, dimension(:),   allocatable, managed :: parts_cross_id_rank
   integer, dimension(:,:), allocatable, managed :: parts_cross_newrank
   integer, dimension(:),   allocatable, managed :: parts_cross_newrank_rank
   integer, dimension(:),   allocatable, managed :: num_part_cross

   ! substantial derivative of velocity
   real(rp), dimension(:,:,:), allocatable, managed :: sdu,sdv,sdw
   integer,  parameter :: CRAZY_INT = 9999999 
   real(rp), parameter :: CRAZY_REAl = 123456789.123456789e-16_rp
   real(rp), parameter :: OneThird = 0.3333333333333333_rp

   integer, parameter :: dragmodel_Stokes = 1 ! Stokes drag
   integer, parameter :: dragmodel_SN = 2     ! Schiller & Nauman
   integer, parameter :: dragmodel_CG = 3     ! Clift & Gauvin
   integer, parameter :: dragmodel_MKL = 4    ! Mei, Klausner,Lawrence 1994

   integer :: dragmodel 
   integer :: ntimesteptag
   real(rp) :: vol_cut, xlpp_min,ylpp_min,zlpp_min, & 
               xlpp_max,ylpp_max,zlpp_max
   real(rp) :: vol_debris
   character(20), managed :: lppbdry_cond(6)

   integer :: max_num_part
   integer :: max_num_part_cross

   integer  :: outputlpp_format
   real(rp) :: ConvertRegSizeToDiam
   
   integer :: nStepConverion

   real(rp), parameter :: AspRatioSphere = 1._rp
   real(rp) :: AspRatioTol

   logical :: WallEffectSettling
   logical :: output_lpp_evolution
   logical :: DoOutputLPP
   
   integer :: TwoWayCouplingFlag
   integer, parameter :: TwoWayCouplingIgnore      = 0
   integer, parameter :: TwoWayCouplingFilterForce = 1
   integer, parameter :: TwoWayCouplingFilterVel   = 2
   real(rp) :: LengthLPP2dp, LengthLPP,mfLPP_ref

   logical :: UnsteadyPartForce

   integer :: SeedParticlesFlag
   logical :: SeedWithFluidVel
   logical :: AvoidSeedAtBlockBdry,AvoidSeedPartTooClose
   integer, parameter :: SeedParticlesNone      = 0
   integer, parameter :: SeedParticlesOneTime   = 1
   integer, parameter :: SeedParticlesRegular   = 2
   real(rp) :: xmin_part_seed, ymin_part_seed, zmin_part_seed, dmin_part_seed, & 
               xmax_part_seed, ymax_part_seed, zmax_part_seed, dmax_part_seed, &
               umin_part_seed, vmin_part_seed, wmin_part_seed, & 
               umax_part_seed, vmax_part_seed, wmax_part_seed
   integer  :: num_part_seed
   real(rp) :: time_part_seed
   integer  :: tag_opened=0

contains
!=================================================================================================
   subroutine initLPP(zmin,zmax,dzc,n,ng)
      implicit none
      integer , intent(in), dimension(3) :: n,ng
      real(rp), intent(in), dimension(0:n(3)+1) :: dzc
      real(rp), intent(in) :: zmin,zmax
      integer :: i,j,k,ipart
      attributes(managed) :: dzc

      call ReadLPPParameters(zmin,zmax,dzc,n)

      allocate( parts(max_num_part,0:(dims(1)*dims(2)-1)) )
      allocate( num_part(0:(dims(1)*dims(2))-1) )
      if ( UnsteadyPartForce ) then 
       allocate( sdu(0:n(1)+1,0:n(2)+1,0:n(3)+1), & 
                 sdv(0:n(1)+1,0:n(2)+1,0:n(3)+1), &
                 sdw(0:n(1)+1,0:n(2)+1,0:n(3)+1)   )
      end if

      if ( TwoWayCouplingFlag == TwoWayCouplingFilterForce) then  
         allocate( parts_collect(max_num_part,0:(dims(1)*dims(2)-1)) )
         allocate( parts_collect_rank(max_num_part) )
      end if

      ! set default values
      num_part = 0

      LPP_initialized = .true.

      if( .not. allocated( cross_buff_id ) ) allocate( cross_buff_id( max_num_part_cross,0:dims(1)*dims(2)-1 ) )
      if( .not. allocated( cross_buff_id_rank ) ) allocate( cross_buff_id_rank( max_num_part_cross ) )
      if( .not. allocated( cross_buff_new ) ) allocate( cross_buff_new( max_num_part_cross,0:dims(1)*dims(2)-1 ) )
      if( .not. allocated( cross_buff_new_rank ) ) allocate( cross_buff_new_rank( max_num_part_cross ) )
      if( .not. allocated( xcOld_r_buf ) ) allocate( xcOld_r_buf( max_num_part_cross,0:(dims(1)*dims(2)-1)) )
      if( .not. allocated( ycOld_r_buf ) ) allocate( ycOld_r_buf( max_num_part_cross,0:(dims(1)*dims(2)-1)) )
      if( .not. allocated( zcOld_r_buf ) ) allocate( zcOld_r_buf( max_num_part_cross,0:(dims(1)*dims(2)-1)) )
      if( .not. allocated( ucOld_r_buf ) ) allocate( ucOld_r_buf( max_num_part_cross,0:(dims(1)*dims(2)-1)) )
      if( .not. allocated( vcOld_r_buf ) ) allocate( vcOld_r_buf( max_num_part_cross,0:(dims(1)*dims(2)-1)) )
      if( .not. allocated( wcOld_r_buf ) ) allocate( wcOld_r_buf( max_num_part_cross,0:(dims(1)*dims(2)-1)) )
      if( .not. allocated( uf_r_buf ) ) allocate( uf_r_buf( max_num_part_cross,0:(dims(1)*dims(2)-1)) )
      if( .not. allocated( vf_r_buf ) ) allocate( vf_r_buf( max_num_part_cross,0:(dims(1)*dims(2)-1)) )
      if( .not. allocated( wf_r_buf ) ) allocate( wf_r_buf( max_num_part_cross,0:(dims(1)*dims(2)-1)) )
      if( .not. allocated( ic_r_buf ) ) allocate( ic_r_buf( max_num_part_cross,0:(dims(1)*dims(2)-1)) )
      if( .not. allocated( jc_r_buf ) ) allocate( jc_r_buf( max_num_part_cross,0:(dims(1)*dims(2)-1)) )
      if( .not. allocated( kc_r_buf ) ) allocate( kc_r_buf( max_num_part_cross,0:(dims(1)*dims(2)-1)) )
      if( .not. allocated( xcOld_s_buf ) ) allocate( xcOld_s_buf( max_num_part_cross,0:(dims(1)*dims(2)-1)) )
      if( .not. allocated( ycOld_s_buf ) ) allocate( ycOld_s_buf( max_num_part_cross,0:(dims(1)*dims(2)-1)) )
      if( .not. allocated( zcOld_s_buf ) ) allocate( zcOld_s_buf( max_num_part_cross,0:(dims(1)*dims(2)-1)) )
      if( .not. allocated( ucOld_s_buf ) ) allocate( ucOld_s_buf( max_num_part_cross,0:(dims(1)*dims(2)-1)) )
      if( .not. allocated( vcOld_s_buf ) ) allocate( vcOld_s_buf( max_num_part_cross,0:(dims(1)*dims(2)-1)) )
      if( .not. allocated( wcOld_s_buf ) ) allocate( wcOld_s_buf( max_num_part_cross,0:(dims(1)*dims(2)-1)) )
      if( .not. allocated( uf_s_buf ) ) allocate( uf_s_buf( max_num_part_cross,0:(dims(1)*dims(2)-1)) )
      if( .not. allocated( vf_s_buf ) ) allocate( vf_s_buf( max_num_part_cross,0:(dims(1)*dims(2)-1)) )
      if( .not. allocated( wf_s_buf ) ) allocate( wf_s_buf( max_num_part_cross,0:(dims(1)*dims(2)-1)) )
      if( .not. allocated( ic_s_buf ) ) allocate( ic_s_buf( max_num_part_cross,0:(dims(1)*dims(2)-1)) )
      if( .not. allocated( jc_s_buf ) ) allocate( jc_s_buf( max_num_part_cross,0:(dims(1)*dims(2)-1)) )
      if( .not. allocated( kc_s_buf ) ) allocate( kc_s_buf( max_num_part_cross,0:(dims(1)*dims(2)-1)) )
#ifndef GPU_MPI
     istat = cudaMemAdvise(        cross_buff_id, size(cross_buff_id), cudaMemAdviseSetPreferredLocation, cudaCpuDeviceId )
     istat = cudaMemAdvise(        cross_buff_id, size(cross_buff_id), cudaMemAdviseSetAccessedBy, mydev )
     istat = cudaMemPrefetchAsync( cross_buff_id, size(cross_buff_id), cudaCpuDeviceId, 0 )
     istat = cudaMemAdvise(        cross_buff_id_rank, size(cross_buff_id_rank), cudaMemAdviseSetPreferredLocation, cudaCpuDeviceId )
     istat = cudaMemAdvise(        cross_buff_id_rank, size(cross_buff_id_rank), cudaMemAdviseSetAccessedBy, mydev )
     istat = cudaMemPrefetchAsync( cross_buff_id_rank, size(cross_buff_id_rank), cudaCpuDeviceId, 0 )
     istat = cudaMemAdvise(        cross_buff_new, size(cross_buff_new), cudaMemAdviseSetPreferredLocation, cudaCpuDeviceId )
     istat = cudaMemAdvise(        cross_buff_new, size(cross_buff_new), cudaMemAdviseSetAccessedBy, mydev )
     istat = cudaMemPrefetchAsync( cross_buff_new, size(cross_buff_new), cudaCpuDeviceId, 0 )
     istat = cudaMemAdvise(        cross_buff_new_rank, size(cross_buff_new_rank), cudaMemAdviseSetPreferredLocation, cudaCpuDeviceId )
     istat = cudaMemAdvise(        cross_buff_new_rank, size(cross_buff_new_rank), cudaMemAdviseSetAccessedBy, mydev )
     istat = cudaMemPrefetchAsync( cross_buff_new_rank, size(cross_buff_new_rank), cudaCpuDeviceId, 0 )
     istat = cudaMemAdvise(        xcOld_r_buf, size(xcOld_r_buf), cudaMemAdviseSetPreferredLocation, cudaCpuDeviceId )
     istat = cudaMemAdvise(        xcOld_r_buf, size(xcOld_r_buf), cudaMemAdviseSetAccessedBy, mydev )
     istat = cudaMemPrefetchAsync( xcOld_r_buf, size(xcOld_r_buf), cudaCpuDeviceId, 0 )
     istat = cudaMemAdvise(        ycOld_r_buf, size(ycOld_r_buf), cudaMemAdviseSetPreferredLocation, cudaCpuDeviceId )
     istat = cudaMemAdvise(        ycOld_r_buf, size(ycOld_r_buf), cudaMemAdviseSetAccessedBy, mydev )
     istat = cudaMemPrefetchAsync( ycOld_r_buf, size(ycOld_r_buf), cudaCpuDeviceId, 0 )
     istat = cudaMemAdvise(        zcOld_r_buf, size(zcOld_r_buf), cudaMemAdviseSetPreferredLocation, cudaCpuDeviceId )
     istat = cudaMemAdvise(        zcOld_r_buf, size(zcOld_r_buf), cudaMemAdviseSetAccessedBy, mydev )
     istat = cudaMemPrefetchAsync( zcOld_r_buf, size(zcOld_r_buf), cudaCpuDeviceId, 0 )
     istat = cudaMemAdvise(        ucOld_r_buf, size(ucOld_r_buf), cudaMemAdviseSetPreferredLocation, cudaCpuDeviceId )
     istat = cudaMemAdvise(        ucOld_r_buf, size(ucOld_r_buf), cudaMemAdviseSetAccessedBy, mydev )
     istat = cudaMemPrefetchAsync( ucOld_r_buf, size(ucOld_r_buf), cudaCpuDeviceId, 0 )
     istat = cudaMemAdvise(        vcOld_r_buf, size(vcOld_r_buf), cudaMemAdviseSetPreferredLocation, cudaCpuDeviceId )
     istat = cudaMemAdvise(        vcOld_r_buf, size(vcOld_r_buf), cudaMemAdviseSetAccessedBy, mydev )
     istat = cudaMemPrefetchAsync( vcOld_r_buf, size(vcOld_r_buf), cudaCpuDeviceId, 0 )
     istat = cudaMemAdvise(        wcOld_r_buf, size(wcOld_r_buf), cudaMemAdviseSetPreferredLocation, cudaCpuDeviceId )
     istat = cudaMemAdvise(        wcOld_r_buf, size(wcOld_r_buf), cudaMemAdviseSetAccessedBy, mydev )
     istat = cudaMemPrefetchAsync( wcOld_r_buf, size(wcOld_r_buf), cudaCpuDeviceId, 0 )
     istat = cudaMemAdvise(        uf_r_buf, size(uf_r_buf), cudaMemAdviseSetPreferredLocation, cudaCpuDeviceId )
     istat = cudaMemAdvise(        uf_r_buf, size(uf_r_buf), cudaMemAdviseSetAccessedBy, mydev )
     istat = cudaMemPrefetchAsync( uf_r_buf, size(uf_r_buf), cudaCpuDeviceId, 0 )
     istat = cudaMemAdvise(        vf_r_buf, size(vf_r_buf), cudaMemAdviseSetPreferredLocation, cudaCpuDeviceId )
     istat = cudaMemAdvise(        vf_r_buf, size(vf_r_buf), cudaMemAdviseSetAccessedBy, mydev )
     istat = cudaMemPrefetchAsync( vf_r_buf, size(vf_r_buf), cudaCpuDeviceId, 0 )
     istat = cudaMemAdvise(        wf_r_buf, size(wf_r_buf), cudaMemAdviseSetPreferredLocation, cudaCpuDeviceId )
     istat = cudaMemAdvise(        wf_r_buf, size(wf_r_buf), cudaMemAdviseSetAccessedBy, mydev )
     istat = cudaMemPrefetchAsync( wf_r_buf, size(wf_r_buf), cudaCpuDeviceId, 0 )
     istat = cudaMemAdvise(        xcOld_s_buf, size(xcOld_s_buf), cudaMemAdviseSetPreferredLocation, cudaCpuDeviceId )
     istat = cudaMemAdvise(        xcOld_s_buf, size(xcOld_s_buf), cudaMemAdviseSetAccessedBy, mydev )
     istat = cudaMemPrefetchAsync( xcOld_s_buf, size(xcOld_s_buf), cudaCpuDeviceId, 0 )
     istat = cudaMemAdvise(        ycOld_s_buf, size(ycOld_s_buf), cudaMemAdviseSetPreferredLocation, cudaCpuDeviceId )
     istat = cudaMemAdvise(        ycOld_s_buf, size(ycOld_s_buf), cudaMemAdviseSetAccessedBy, mydev )
     istat = cudaMemPrefetchAsync( ycOld_s_buf, size(ycOld_s_buf), cudaCpuDeviceId, 0 )
     istat = cudaMemAdvise(        zcOld_s_buf, size(zcOld_s_buf), cudaMemAdviseSetPreferredLocation, cudaCpuDeviceId )
     istat = cudaMemAdvise(        zcOld_s_buf, size(zcOld_s_buf), cudaMemAdviseSetAccessedBy, mydev )
     istat = cudaMemPrefetchAsync( zcOld_s_buf, size(zcOld_s_buf), cudaCpuDeviceId, 0 )
     istat = cudaMemAdvise(        ucOld_s_buf, size(ucOld_s_buf), cudaMemAdviseSetPreferredLocation, cudaCpuDeviceId )
     istat = cudaMemAdvise(        ucOld_s_buf, size(ucOld_s_buf), cudaMemAdviseSetAccessedBy, mydev )
     istat = cudaMemPrefetchAsync( ucOld_s_buf, size(ucOld_s_buf), cudaCpuDeviceId, 0 )
     istat = cudaMemAdvise(        vcOld_s_buf, size(vcOld_s_buf), cudaMemAdviseSetPreferredLocation, cudaCpuDeviceId )
     istat = cudaMemAdvise(        vcOld_s_buf, size(vcOld_s_buf), cudaMemAdviseSetAccessedBy, mydev )
     istat = cudaMemPrefetchAsync( vcOld_s_buf, size(vcOld_s_buf), cudaCpuDeviceId, 0 )
     istat = cudaMemAdvise(        wcOld_s_buf, size(wcOld_s_buf), cudaMemAdviseSetPreferredLocation, cudaCpuDeviceId )
     istat = cudaMemAdvise(        wcOld_s_buf, size(wcOld_s_buf), cudaMemAdviseSetAccessedBy, mydev )
     istat = cudaMemPrefetchAsync( wcOld_s_buf, size(wcOld_s_buf), cudaCpuDeviceId, 0 )
     istat = cudaMemAdvise(        uf_s_buf, size(uf_s_buf), cudaMemAdviseSetPreferredLocation, cudaCpuDeviceId )
     istat = cudaMemAdvise(        uf_s_buf, size(uf_s_buf), cudaMemAdviseSetAccessedBy, mydev )
     istat = cudaMemPrefetchAsync( uf_s_buf, size(uf_s_buf), cudaCpuDeviceId, 0 )
     istat = cudaMemAdvise(        vf_s_buf, size(vf_s_buf), cudaMemAdviseSetPreferredLocation, cudaCpuDeviceId )
     istat = cudaMemAdvise(        vf_s_buf, size(vf_s_buf), cudaMemAdviseSetAccessedBy, mydev )
     istat = cudaMemPrefetchAsync( vf_s_buf, size(vf_s_buf), cudaCpuDeviceId, 0 )
     istat = cudaMemAdvise(        wf_s_buf, size(wf_s_buf), cudaMemAdviseSetPreferredLocation, cudaCpuDeviceId )
     istat = cudaMemAdvise(        wf_s_buf, size(wf_s_buf), cudaMemAdviseSetAccessedBy, mydev )
     istat = cudaMemPrefetchAsync( wf_s_buf, size(wf_s_buf), cudaCpuDeviceId, 0 )
#endif
   if ( UnsteadyPartForce ) then 
    istat = cudaMemAdvise( sdu,     size(sdu),      cudaMemAdviseSetPreferredLocation, mydev )
    istat = cudaMemAdvise( sdv,     size(sdv),      cudaMemAdviseSetPreferredLocation, mydev )
    istat = cudaMemAdvise( sdw,     size(sdw),      cudaMemAdviseSetPreferredLocation, mydev )
   endif
   istat = cudaMemAdvise( parts%ic, size(parts%ic), cudaMemAdviseSetPreferredLocation, mydev )
   istat = cudaMemAdvise( parts%jc, size(parts%jc), cudaMemAdviseSetPreferredLocation, mydev )
   istat = cudaMemAdvise( parts%kc, size(parts%kc), cudaMemAdviseSetPreferredLocation, mydev )
   istat = cudaMemAdvise( parts%uf, size(parts%uf), cudaMemAdviseSetPreferredLocation, mydev )
   istat = cudaMemAdvise( parts%vf, size(parts%vf), cudaMemAdviseSetPreferredLocation, mydev )
   istat = cudaMemAdvise( parts%wf, size(parts%wf), cudaMemAdviseSetPreferredLocation, mydev )
   istat = cudaMemAdvise( parts%xcOld, size(parts%xcOld), cudaMemAdviseSetPreferredLocation, mydev )
   istat = cudaMemAdvise( parts%ycOld, size(parts%ycOld), cudaMemAdviseSetPreferredLocation, mydev )
   istat = cudaMemAdvise( parts%zcOld, size(parts%zcOld), cudaMemAdviseSetPreferredLocation, mydev )
   istat = cudaMemAdvise( parts%ucOld, size(parts%ucOld), cudaMemAdviseSetPreferredLocation, mydev )
   istat = cudaMemAdvise( parts%vcOld, size(parts%vcOld), cudaMemAdviseSetPreferredLocation, mydev )
   istat = cudaMemAdvise( parts%wcOld, size(parts%wcOld), cudaMemAdviseSetPreferredLocation, mydev )
   istat = cudaMemAdvise( parts%xc, size(parts%xc), cudaMemAdviseSetPreferredLocation, mydev )
   istat = cudaMemAdvise( parts%yc, size(parts%yc), cudaMemAdviseSetPreferredLocation, mydev )
   istat = cudaMemAdvise( parts%zc, size(parts%zc), cudaMemAdviseSetPreferredLocation, mydev )
   istat = cudaMemAdvise( parts%uc, size(parts%uc), cudaMemAdviseSetPreferredLocation, mydev )
   istat = cudaMemAdvise( parts%vc, size(parts%vc), cudaMemAdviseSetPreferredLocation, mydev )
   istat = cudaMemAdvise( parts%wc, size(parts%wc), cudaMemAdviseSetPreferredLocation, mydev )
   istat = cudaMemAdvise( parts%duc, size(parts%duc), cudaMemAdviseSetPreferredLocation, mydev )
   istat = cudaMemAdvise( parts%dvc, size(parts%dvc), cudaMemAdviseSetPreferredLocation, mydev )
   istat = cudaMemAdvise( parts%dwc, size(parts%dwc), cudaMemAdviseSetPreferredLocation, mydev )
   istat = cudaMemAdvise( parts%vol, size(parts%vol), cudaMemAdviseSetPreferredLocation, mydev )
   istat = cudaMemAdvise( parts%sur, size(parts%sur), cudaMemAdviseSetPreferredLocation, mydev )
   istat = cudaMemAdvise( num_part, size(num_part), cudaMemAdviseSetPreferredLocation, mydev )

  istat = cudaMemPrefetchAsync( sdu, size(sdu), mydev, 0 )
  istat = cudaMemPrefetchAsync( sdv, size(sdv), mydev, 0 )
  istat = cudaMemPrefetchAsync( sdw, size(sdw), mydev, 0 )

  if ( UnsteadyPartForce ) then
   !$cuf kernel do(3) <<<*,*>>> 
   do k=1,n(3)
     do j=1,n(2)
       do i=1,n(1)
         sdu(i,j,k) = 0._rp
         sdv(i,j,k) = 0._rp
         sdw(i,j,k) = 0._rp
       enddo
     enddo
   enddo
  end if

   end subroutine initLPP

   subroutine ReadLPPParameters(zmin,zmax,dzc,n)
      implicit none
      integer , intent(in), dimension(3) :: n
      real(rp), intent(in), dimension(0:n(3)+1) :: dzc
      real(rp), intent(in) :: zmin,zmax
      integer :: i
      logical file_is_there
      namelist /lppparameters/ dragmodel, nTimeStepTag,   &
                               ConvertRegSizeToDiam, & 
                               vol_cut,xlpp_min,xlpp_max,ylpp_min,ylpp_max,zlpp_min,zlpp_max,    &
                               max_num_part, max_num_part_cross, &
                               outputlpp_format,nStepConverion,AspRatioTol, &
                               WallEffectSettling,output_lpp_evolution,lppbdry_cond,& 
                               TwoWayCouplingFlag,LengthLPP2dp,mfLPP_ref,UnsteadyPartForce,&  
                               SeedParticlesFlag, SeedWithFluidVel, & 
                               AvoidSeedAtBlockBdry,AvoidSeedPartTooClose, & 
                               num_part_seed, time_part_seed, & 
                               xmin_part_seed, ymin_part_seed, zmin_part_seed, dmin_part_seed, & 
                               xmax_part_seed, ymax_part_seed, zmax_part_seed, dmax_part_seed, &
                               umin_part_seed, vmin_part_seed, wmin_part_seed, & 
                               umax_part_seed, vmax_part_seed, wmax_part_seed, & 
                               DoOutputLPP,vol_debris

      ! Set default values 
      dragmodel    = 1
      nTimeStepTag = 10
      vol_cut  = 1.e-9_rp
      xlpp_min = (1   +coord(1)*n(1))*dx-0.5_rp*dx
      xlpp_max = (n(1)+coord(1)*n(1))*dx-0.5_rp*dx
      ylpp_min = (1   +coord(2)*n(2))*dy-0.5_rp*dy
      ylpp_max = (n(2)+coord(2)*n(2))*dy-0.5_rp*dy
      zlpp_min = zmin
      zlpp_max = zmax
      max_num_part = 2
      max_num_part_cross = 2
      outputlpp_format = 1
      ConvertRegSizeToDiam = 2._rp
      nStepConverion = 0
      AspRatioTol = 1.5_rp
      WallEffectSettling = .false.
      output_lpp_evolution = .false.
      DoOutputLpp = .true.
      lppbdry_cond=['undefined','undefined','undefined','undefined','undefined','undefined']
      TwoWayCouplingFlag = 0 
      LengthLPP2dp = 8._rp
      mfLPP_ref = 0.5_rp
      UnsteadyPartForce = .false.
      SeedParticlesFlag = 0
      SeedWithFluidVel =.false. 
      AvoidSeedAtBlockBdry =.false. 
      AvoidSeedPartTooClose =.false. 
      num_part_seed = 0; time_part_seed = 0._rp;
      xmin_part_seed=0._rp; ymin_part_seed=0._rp; zmin_part_seed=0._rp; dmin_part_seed=0._rp 
      xmax_part_seed=0._rp; ymax_part_seed=0._rp; zmax_part_seed=0._rp; dmax_part_seed=0._rp
      umin_part_seed=0._rp; vmin_part_seed=0._rp; wmin_part_seed=0._rp
      umax_part_seed=0._rp; vmax_part_seed=0._rp; wmax_part_seed=0._rp; 
      vol_debris = dx*dy*minval(dzc)

      inquire(file='lpp.in',exist=file_is_there)
      open(unit=32, file='lpp.in', status='old', action='read', iostat=ierr)
      if (file_is_there) then
         if(ierr == 0) then
            read(UNIT=32,NML=lppparameters)
            if(myid.eq.0) print*, 'Lagrangian point-particle parameters read successfully'
         else
            if(myid.eq.0) print*, 'Error LPP input file' 
            if(myid.eq.0) print*, 'Aborting...'
            call MPI_FINALIZE(ierr)
            call exit
         endif
      else
         if(myid.eq.0) print*, 'No LPP input file' 
         if(myid.eq.0) print*, 'Aborting...'
         call MPI_FINALIZE(ierr)
         call exit
      endif
      close(32)

      do i=1,6
         if(lppbdry_cond(i) == 'undefined') then 
            lppbdry_cond(i) = 'periodic'
            if ( myid == 0 ) write(*,*) 'Undefined lpp condition in direction', i,' now is set to periodic by default!' 
         end if ! lppbdry_cond 
      end do !i
      
   istat = cudaMemPrefetchAsync( lppbdry_cond, size(lppbdry_cond), mydev, 0 )

   end subroutine ReadLPPParameters

   subroutine lppsweeps(time,dt,up,vp,wp,u,v,w,zf,dzf,zc,n,ng,istage)
      implicit none
      integer , intent(in), dimension(3) :: n,ng
      integer,  intent(in) :: istage 
      real(rp), dimension(0:,0:,0:), intent(inout) :: up,vp,wp
      real(rp), dimension(0:n(3)+1), intent(in) :: zf,dzf,zc
      real(rp), dimension(0:n(1)+1,0:n(2)+1,0:n(3)+1), intent(in) :: u,v,w
      real(rp), intent(in) :: time,dt
      attributes(managed) :: up,vp,wp,u,v,w,zf,dzf,zc

      call ComputePartPropertiesGPU(u,v,w,time,dt,zc,zf,n,ng)
      !@cuf istat=cudaDeviceSynchronize()
      if ( TwoWayCouplingFlag.ne.TwoWayCouplingIgnore ) then
       call TwoWayCouplingForce(TwoWayCouplingFlag,up,vp,wp,zc,zf,dzf,n,ng)
       !@cuf istat=cudaDeviceSynchronize()
      endif
      call UpdatePartSol(dt,zf,dzf,n,ng,istage)
      !@cuf istat=cudaDeviceSynchronize()

   end subroutine lppsweeps

   function Rep(uf,up,vf,vp,wf,wp,nuf,dp)
      implicit none
      real(rp), intent(in) :: uf,up,vf,vp,wf,wp,nuf,dp 
      real(rp) :: Rep

      Rep = sqrt((uf-up)**2 + (vf-vp)**2 + (wf-wp)**2)*dp/nuf

   end function Rep

   subroutine Surface2VolumeIntrpl(ux1,ux2,uy1,uy2,uz1,uz2,x1,x2,y1,y2,z1,z2,uf)
      implicit none

      real(rp), intent(in)  :: ux1,ux2,uy1,uy2,uz1,uz2,x1,x2,y1,y2,z1,z2
      real(rp), intent(out) :: uf
      real(rp) :: wtsum,x1p,x2p,y1p,y2p,z1p,z2p

      x1p = 1._rp/x1 
      x2p = 1._rp/x2  
      y1p = 1._rp/y1
      y2p = 1._rp/y2
      z1p = 1._rp/z1
      z2p = 1._rp/z2
      wtsum = x1p + x2p + y1p + y2p + z1p + z2p 
       
      uf =(x1p*ux1 + x2p*ux2 & 
         + y1p*uy1 + y2p*uy2 & 
         + z1p*uz1 + z2p*uz2)/wtsum

   end subroutine Surface2VolumeIntrpl

   subroutine FindPartLocCell(xp,yp,zp,ip,jp,kp,zf,n,ng)
      implicit none
      integer , intent(in), dimension(3) :: n,ng
      real(rp), dimension(0:n(3)+1),intent(in) :: zf
      real(rp), intent(in)  :: xp,yp,zp
      integer,  intent(out) :: ip,jp,kp

      if ( xp < (0+coord(1)*n(1))*dx .or. xp > ((n(1)+1)+coord(1)*n(1))*dx .or. & 
           yp < (0+coord(2)*n(2))*dy .or. yp > ((n(2)+1)+coord(2)*n(2))*dy .or. & 
           zp < zf(0) .or. zp > zf(n(3)+1) )    & 
           call lpperror("Failed to find cell index for particle outside of domain!") 

      do ip = 0+coord(1)*n(1),(n(1)+1)+coord(1)*n(1)
         if (xp <= (ip*dx)) exit
      end do ! ip

      do jp = 0+coord(2)*n(2),(n(2)+1)+coord(2)*n(2)
         if (yp <= (jp*dy)) exit
      end do ! jp

      do kp = 0,n(3)+1
         if (zp <= zf(kp)) exit
      end do ! kp

   end subroutine 

   ! subroutine FindCellIndexBdryConvertReg(xc,yc,zc,l,i1,ic,i2,j1,jc,j2,k1,kc,k2,zf)
      ! implicit none
      ! real(rp), dimension(0:n(3)),intent(in) :: zf
      ! real(rp), intent (in) :: xc,yc,zc,l
      ! integer, intent (in) :: ic,jc,kc
      ! integer, intent(out) :: i1,i2,j1,j2,k1,k2
      ! real(rp) :: l2,xl,xr,yl,yr,zl,zr

      ! l2 = 0.5*l

      ! xl = xc - l2
      ! do i1 = ic,0,-1
         ! if ( xl > (i1+coord(1)*imax)*dx ) exit 
      ! end do ! i
      
      ! xr = xc + l2
      ! do i2 = ic,n(1)
         ! if ( xr < (i2+coord(1)*imax)*dx ) exit 
      ! end do ! i

      ! yl = yc - l2
      ! do j1 = jc,0,-1
         ! if ( yl > (j1+coord(2)*jmax)*dy ) exit 
      ! end do ! j

      ! yr = yc + l2
      ! do j2 = jc,n(2)
         ! if ( yr < (j2+coord(2)*jmax)*dy ) exit 
      ! end do ! j

      ! zl = zc - l2
      ! do k1 = kc,0,-1
         ! if ( zl > zf(k1-1) ) exit 
      ! end do ! k

      ! zr = zc + l2
      ! do k2 = kc,n(3)
         ! if ( zr < zf(k2  ) ) exit 
      ! end do ! k

   ! end subroutine FindCellIndexBdryConvertReg

   ! subroutine FindCellIndexBdryConvertRegUnifMesh(l,i1,ic,i2,j1,jc,j2,k1,kc,k2)
      ! implicit none

      ! real(rp), intent (in) :: l
      ! integer, intent (in) :: ic,jc,kc
      ! integer, intent(out) :: i1,i2,j1,j2,k1,k2
      ! real(rp) :: l2
      ! integer :: ncl2

      ! l2 = 0.5_rp*l
      ! ncl2 = NINT(l2/dx)
      ! i1 = ic - ncl2
      ! i2 = ic + ncl2
      ! j1 = jc - ncl2
      ! j2 = jc + ncl2
      ! k1 = kc - ncl2
      ! k2 = kc + ncl2

      ! i1 = max(0,i1)
      ! i2 = min(n(1),i2)
      ! j1 = max(0,j1)
      ! j2 = min(n(2),j2)
      ! k1 = max(0,k1)
      ! k2 = min(n(3),k2)
   ! end subroutine FindCellIndexBdryConvertRegUnifMesh

   subroutine ComputePartPropertiesGPU(u,v,w,time,dt,z_c,z_f,n,ng)
      implicit none
      integer , intent(in), dimension(3) :: n,ng
      real(rp), dimension(0:n(1)+1,0:n(2)+1,0:n(3)+1), intent(in) :: u,v,w
      real(rp), dimension(0:n(3)+1), intent(in) :: z_c,z_f
      real(rp), intent(in) :: time,dt
      real(rp), dimension(:), allocatable :: relvel, partforce
      real(rp) :: Cm
      real(rp) :: dp, mp, Rep, muf, mup, rhof, rhop, taup
      real(rp) :: phi_c
      real(rp) :: xp,yp,zp
      real(rp) :: up,vp,wp,uf,vf,wf,DufDt,DvfDt,DwfDt
      real(rp) :: fhx,fhy,fhz,gx,gy,gz
      real(rp) :: dp2Lx
      real(rp) :: xh1,xh2,yh1,yh2,temp1,temp2,temp3,temp4,temp5,temp6, &
                  xl,xr,yl,yr,zl,zr, &
                  switch_stokes, switch_SN, &
                  switch_CG, switch_MKL
      integer  :: ip,jp,kp,si,sj,sk
      integer  :: ipart
      integer  :: MPI_errorcode, switch_wall, switch_force
      attributes(managed) :: u,v,w,z_c,z_f,relvel,partforce

      ! write(*,*)"\Rank=",myid,"Entered particle calculation subroutine"
      
      allocate( relvel(1:4), partforce(1:3) )

      rhof = rho; rhop = rho
      muf  = rho*visc; mup  = rho*visc
      
      gx = gacc(1); gy = gacc(2); gz = gacc(3)

      MPI_errorcode = 0
      switch_wall = 0._rp; switch_force = 0._rp
      switch_stokes = 0._rp; switch_SN = 0._rp
      switch_CG = 0._rp; switch_MKL = 0._rp
      
      Cm = 0.5

      if ( WallEffectSettling ) switch_wall = 1
      if ( UnsteadyPartForce )  switch_force = 1

      if ( dragmodel == dragmodel_Stokes ) switch_stokes = 1._rp
      if ( dragmodel == dragmodel_SN )     switch_SN = 1._rp
      if ( dragmodel == dragmodel_CG )     switch_CG = 1._rp  
      if ( dragmodel == dragmodel_MKL )    switch_MKL = 1._rp  

      if ( num_part(mydev) > 0 ) then
       !$cuf kernel do(1) <<<*,*>>>
       do ipart = 1,num_part(mydev)

        uf = 0._rp; vf = 0._rp; wf = 0._rp
        DufDt = 0._rp; DvfDt = 0._rp; DwfDt = 0._rp
        parts(ipart,mydev)%uf = 0._rp; parts(ipart,mydev)%vf = 0._rp; parts(ipart,mydev)%wf = 0._rp

        ip = parts(ipart,mydev)%ic; jp = parts(ipart,mydev)%jc; kp = parts(ipart,mydev)%kc
        xp = parts(ipart,mydev)%xc; yp = parts(ipart,mydev)%yc; zp = parts(ipart,mydev)%zc 
        up = parts(ipart,mydev)%uc; vp = parts(ipart,mydev)%vc; wp = parts(ipart,mydev)%wc

        ! write(*,*)"\Rank=",myid," \P-Pos-B ",xp,yp,zp

        dp = (parts(ipart,mydev)%vol*6._rp/PI)**(1._rp/3._rp)
        taup = rhop*dp*dp/18._rp/muf*(3._rp + 3._rp*muf/mup)/(3._rp + 2._rp*muf/mup)
        partforce(1) = 0._rp; partforce(2) = 0._rp; partforce(3) = 0._rp
        relvel(1) = 0._rp; relvel(2) = 0._rp; relvel(3) = 0._rp; relvel(4) = 0._rp

        !! Get Fluid Property
            ! Trilinear interpolation for u velocity
             si = -1;sj = -1;sk = -1
             
             if ( yp > (jp*dy-0.5_rp*dy) ) sj =  0
             if ( zp > z_c(kp) ) sk =  0
             
             xh1 = (ip+si)*dx; xh2 = (ip+1+si)*dx
             yh1 = (jp+sj)*dy; yh2 = (jp+1+sj)*dy
             
             xl = (xp-xh1)/(xh2-xh1)
             xr = 1._rp - xh2
             temp1  = u(ip+si, jp+sj, kp+sk)*xr + u(ip+1+si, jp+sj, kp+sk)*xl
             
             xl = (xp-xh1)/(xh2-xh1)
             xr = 1._rp - xh2
             temp2  = u(ip+si, jp+1+sj, kp+sk)*xr + u(ip+1+si, jp+1+sj, kp+sk)*xl
             
             xl = (xp-xh1)/(xh2-xh1)
             xr = 1._rp - xh2
             temp3  = u(ip+si, jp+sj, kp+1+sk)*xr + u(ip+1+si, jp+sj, kp+1+sk)*xl
             
             xl = (xp-xh1)/(xh2-xh1)
             xr = 1._rp - xh2
             temp4  = u(ip+si, jp+1+sj, kp+1+sk)*xr + u(ip+1+si, jp+1+sj, kp+1+sk)*xl
             
             yl = (yp-yh1)/(yh2-yh1)
             yr = 1._rp - yh2
             temp5  = temp1*yr + temp2*yl
             
             yl = (yp-yh1)/(yh2-yh1)
             yr = 1._rp - yh2
             temp6  = temp3*yr + temp4*yl
             
             zl = (zp-z_f(kp+sk))/(z_f(kp+1+sk)-z_f(kp+sk))
             zr = 1._rp - z_f(kp+1+sk)
             uf  = temp5*zr + temp6*zl
             if ( switch_force == 1 ) then
              xl = (xp-xh1)/(xh2-xh1)
              xr = 1._rp - xh2
              temp1  = sdu(ip+si, jp+sj, kp+sk)*xr + sdu(ip+1+si, jp+sj, kp+sk)*xl
              
              xl = (xp-xh1)/(xh2-xh1)
              xr = 1._rp - xh2
              temp2  = sdu(ip+si, jp+1+sj, kp+sk)*xr + sdu(ip+1+si, jp+1+sj, kp+sk)*xl
              
              xl = (xp-xh1)/(xh2-xh1)
              xr = 1._rp - xh2
              temp3  = sdu(ip+si, jp+sj, kp+1+sk)*xr + sdu(ip+1+si, jp+sj, kp+1+sk)*xl
              
              xl = (xp-xh1)/(xh2-xh1)
              xr = 1._rp - xh2
              temp4  = sdu(ip+si, jp+1+sj, kp+1+sk)*xr + sdu(ip+1+si, jp+1+sj, kp+1+sk)*xl
              
              yl = (yp-yh1)/(yh2-yh1)
              yr = 1._rp - yh2
              temp5  = temp1*yr + temp2*yl
              
              yl = (yp-yh1)/(yh2-yh1)
              yr = 1._rp - yh2
              temp6  = temp3*yr + temp4*yl
              
              zl = (zp-z_f(kp+sk))/(z_f(kp+1+sk)-z_f(kp+sk))
              zr = 1._rp - z_f(kp+1+sk)
              DufDt  = temp5*zr + temp6*zl
             end if
            ! Trilinear interpolation for v velocity
             si = -1;sj = -1;sk = -1
             
             if ( xp > (ip*dx-0.5_rp*dx) ) si =  0 
             if ( zp > z_c(kp) ) sk =  0
             
             xh1 = (ip+si)*dx; xh2 = (ip+1+si)*dx
             yh1 = (jp+sj)*dy; yh2 = (jp+1+sj)*dy

             xl = (xp-xh1)/(xh2-xh1)
             xr = 1._rp - xh2
             temp1  = v(ip+si, jp+sj, kp+sk)*xr + v(ip+1+si, jp+sj, kp+sk)*xl
             
             xl = (xp-xh1)/(xh2-xh1)
             xr = 1._rp - xh2
             temp2  = v(ip+si, jp+1+sj, kp+sk)*xr + v(ip+1+si, jp+1+sj, kp+sk)*xl
             
             xl = (xp-xh1)/(xh2-xh1)
             xr = 1._rp - xh2
             temp3  = v(ip+si, jp+sj, kp+1+sk)*xr + v(ip+1+si, jp+sj, kp+1+sk)*xl
             
             xl = (xp-xh1)/(xh2-xh1)
             xr = 1._rp - xh2
             temp4  = v(ip+si, jp+1+sj, kp+1+sk)*xr + v(ip+1+si, jp+1+sj, kp+1+sk)*xl
             
             yl = (yp-yh1)/(yh2-yh1)
             yr = 1._rp - yh2
             temp5  = temp1*yr + temp2*yl
             
             yl = (yp-yh1)/(yh2-yh1)
             yr = 1._rp - yh2
             temp6  = temp3*yr + temp4*yl
             
             zl = (zp-z_f(kp+sk))/(z_f(kp+1+sk)-z_f(kp+sk))
             zr = 1._rp - z_f(kp+1+sk)
             vf  = temp5*zr + temp6*zl
             if ( switch_force == 1 ) then
              xl = (xp-xh1)/(xh2-xh1)
              xr = 1._rp - xh2
              temp1  = sdv(ip+si, jp+sj, kp+sk)*xr + sdv(ip+1+si, jp+sj, kp+sk)*xl
              
              xl = (xp-xh1)/(xh2-xh1)
              xr = 1._rp - xh2
              temp2  = sdv(ip+si, jp+1+sj, kp+sk)*xr + sdv(ip+1+si, jp+1+sj, kp+sk)*xl
              
              xl = (xp-xh1)/(xh2-xh1)
              xr = 1._rp - xh2
              temp3  = sdv(ip+si, jp+sj, kp+1+sk)*xr + sdv(ip+1+si, jp+sj, kp+1+sk)*xl
              
              xl = (xp-xh1)/(xh2-xh1)
              xr = 1._rp - xh2
              temp4  = sdv(ip+si, jp+1+sj, kp+1+sk)*xr + sdv(ip+1+si, jp+1+sj, kp+1+sk)*xl
              
              yl = (yp-yh1)/(yh2-yh1)
              yr = 1._rp - yh2
              temp5  = temp1*yr + temp2*yl
              
              yl = (yp-yh1)/(yh2-yh1)
              yr = 1._rp - yh2
              temp6  = temp3*yr + temp4*yl
              
              zl = (zp-z_f(kp+sk))/(z_f(kp+1+sk)-z_f(kp+sk))
              zr = 1._rp - z_f(kp+1+sk)
              DvfDt  = temp5*zr + temp6*zl
             end if
            ! Trilinear interpolation for w
             si = -1;sj = -1;sk = -1
             
             if ( xp > (ip*dx-0.5_rp*dx) ) si =  0
             if ( yp > (jp*dy-0.5_rp*dy) ) sj =  0
             
             xh1 = (ip+si)*dx; xh2 = (ip+1+si)*dx
             yh1 = (jp+sj)*dy; yh2 = (jp+1+sj)*dy

             xl = (xp-xh1)/(xh2-xh1)
             xr = 1._rp - xh2
             temp1  = w(ip  +si,jp  +sj,kp  +sk)*xr + w(ip+1+si,jp  +sj,kp  +sk)*xl
             
             xl = (xp-xh1)/(xh2-xh1)
             xr = 1._rp - xh2
             temp2  = w(ip  +si,jp+1+sj,kp  +sk)*xr + w(ip+1+si,jp+1+sj,kp  +sk)*xl
             
             xl = (xp-xh1)/(xh2-xh1)
             xr = 1._rp - xh2
             temp3  = w(ip  +si,jp  +sj,kp+1+sk)*xr + w(ip+1+si,jp  +sj,kp+1+sk)*xl
             
             xl = (xp-xh1)/(xh2-xh1)
             xr = 1._rp - xh2
             temp4  = w(ip  +si,jp+1+sj,kp+1+sk)*xr + w(ip+1+si,jp+1+sj,kp+1+sk)*xl
             
             yl = (yp-yh1)/(yh2-yh1)
             yr = 1._rp - yh2
             temp5  = temp1*yr + temp2*yl
             
             yl = (yp-yh1)/(yh2-yh1)
             yr = 1._rp - yh2
             temp6  = temp3*yr + temp4*yl
             
             zl = (zp-z_f(kp  +sk))/(z_f(kp+1+sk)-z_f(kp  +sk))
             zr = 1._rp - z_f(kp+1+sk)
             wf  = temp5*zr + temp6*zl
             if ( switch_force == 1 ) then
              xl = (xp-xh1)/(xh2-xh1)
              xr = 1._rp - xh2
              temp1  = sdw(ip  +si,jp  +sj,kp  +sk)*xr + sdw(ip+1+si,jp  +sj,kp  +sk)*xl
              
              xl = (xp-xh1)/(xh2-xh1)
              xr = 1._rp - xh2
              temp2  = sdw(ip  +si,jp+1+sj,kp  +sk)*xr + sdw(ip+1+si,jp+1+sj,kp  +sk)*xl
              
              xl = (xp-xh1)/(xh2-xh1)
              xr = 1._rp - xh2
              temp3  = sdw(ip  +si,jp  +sj,kp+1+sk)*xr + sdw(ip+1+si,jp  +sj,kp+1+sk)*xl
              
              xl = (xp-xh1)/(xh2-xh1)
              xr = 1._rp - xh2
              temp4  = sdw(ip  +si,jp+1+sj,kp+1+sk)*xr + sdw(ip+1+si,jp+1+sj,kp+1+sk)*xl
              
              yl = (yp-yh1)/(yh2-yh1)
              yr = 1._rp - yh2
              temp5  = temp1*yr + temp2*yl
              
              yl = (yp-yh1)/(yh2-yh1)
              yr = 1._rp - yh2
              temp6  = temp3*yr + temp4*yl
              
              zl = (zp-z_f(kp  +sk))/(z_f(kp+1+sk)-z_f(kp  +sk))
              zr = 1._rp - z_f(kp+1+sk)
              DwfDt  = temp5*zr + temp6*zl
             end if

        ! write(*,*)"\Rank=",myid," \P-Pos-A ",xp,yp,zp

        if ( taup < 5._rp*dt ) then ! If taup is much smaller than dt, no need to compute force, particle velocity is set to fluid velocity directly
           partforce(1) = CRAZY_REAL
           partforce(2) = CRAZY_REAL
           partforce(3) = CRAZY_REAL
           parts(ipart,mydev)%uc = uf
           parts(ipart,mydev)%vc = vf
           parts(ipart,mydev)%wc = wf
        else
           relvel(1) = uf - up
           relvel(2) = vf - vp
           relvel(3) = wf - wp 
           relvel(4) = sqrt(relvel(1)**2 + relvel(2)**2 + relvel(3)**2)
           Rep  = rhof*relvel(4)*dp/muf
           
           if ( switch_wall == 1 ) then ! correction correlation by De Felice 1996
              dp2Lx = dp/lz   ! assuming Lx=Lz, y is settling direction
              taup  = taup*((1._rp-dp2Lx)/(1._rp-0.33_rp*dp2Lx))**2.7_rp
           end if
           
           fhx=0._rp; fhy=0._rp; fhz=0._rp ! Note: set history to be zero for now
           
           phi_c = 0._rp
            
           phi_c = phi_c + switch_stokes*1._rp
           phi_c = phi_c + switch_SN*(1._rp + 0.15_rp*Rep**0.687_rp)  
           phi_c = phi_c + switch_CG*(1._rp + 0.15_rp*Rep**0.687_rp + 1.75e-2_rp*Rep/(1.0_rp + 4.25e4_rp/Rep**1.16_rp))
           phi_c = phi_c + switch_MKL*(1._rp + 1.0_rp/(8._rp/Rep + 0.5_rp *(1._rp + 3.315_rp/Rep**0.5)))
           
           if ( switch_force == 1 ) then ! [BUG] Causes illegal memory access error
              partforce(1) =((uf - up)/taup*phi_c + (1._rp-rhof/rhop)*gx + (1._rp+Cm)*rhof/rhop*DufDt + fhx )/(1._rp+Cm*rhof/rhop)
              partforce(2) =((vf - vp)/taup*phi_c + (1._rp-rhof/rhop)*gy + (1._rp+Cm)*rhof/rhop*DvfDt + fhy )/(1._rp+Cm*rhof/rhop)
              partforce(3) =((wf - wp)/taup*phi_c + (1._rp-rhof/rhop)*gz + (1._rp+Cm)*rhof/rhop*DwfDt + fhz )/(1._rp+Cm*rhof/rhop)
           else
              partforce(1) = (uf - up)/taup*phi_c + (1._rp-rhof/rhop)*gx
              partforce(2) = (vf - vp)/taup*phi_c + (1._rp-rhof/rhop)*gy
              partforce(3) = (wf - wp)/taup*phi_c + (1._rp-rhof/rhop)*gz
           end if ! UnsteadyPartForce
        end if  ! taup 

        ! if ( partforce(1) /= partforce(1) ) then 
           ! write(*,*) ipart,mydev,num_part(mydev),UnsteadyPartForce,partforce(1:3),relvel(1:3),taup
           ! write(*,*) parts(ipart,mydev)%ic,parts(ipart,mydev)%jc,parts(ipart,mydev)%kc,parts(ipart,mydev)%vol
           ! write(*,*) xp,yp,zp,up,vp,wp
           ! MPI_errorcode=1
           ! print*, "LPP ERROR *** ","particle force in x direction is NaN!"," *** STOP at processor: ", mydev
           ! exit
        ! else if ( partforce(2) /= partforce(2) ) then 
           ! write(*,*) ipart,mydev,num_part(mydev),UnsteadyPartForce,partforce(1:3),relvel(1:3),taup
           ! write(*,*) parts(ipart,mydev)%ic,parts(ipart,mydev)%jc,parts(ipart,mydev)%kc,parts(ipart,mydev)%vol
           ! write(*,*) xp,yp,zp,up,vp,wp
           ! MPI_errorcode=1
           ! print*, "LPP ERROR *** ","particle force in y direction is NaN!"," *** STOP at processor: ", mydev
           ! exit
        ! else if ( partforce(3) /= partforce(3) ) then 
           ! write(*,*) ipart,mydev,num_part(mydev),UnsteadyPartForce,partforce(1:3),relvel(1:3),taup
           ! write(*,*) parts(ipart,mydev)%ic,parts(ipart,mydev)%jc,parts(ipart,mydev)%kc,parts(ipart,mydev)%kc,parts(ipart,mydev)%vol
           ! write(*,*) xp,yp,zp,up,vp,wp
           ! MPI_errorcode=1
           ! print*, "LPP ERROR *** ","particle force in z direction is NaN!"," *** STOP at processor: ", mydev
           ! exit
        ! else
           ! parts(ipart,mydev)%duc = partforce(1) 
           ! parts(ipart,mydev)%dvc = partforce(2)
           ! parts(ipart,mydev)%dwc = partforce(3)
        ! end if ! partforce

        parts(ipart,mydev)%uf = uf
        parts(ipart,mydev)%vf = vf
        parts(ipart,mydev)%wf = wf

        if ( TwoWayCouplingFlag == TwoWayCouplingFilterForce) then 
           parts_collect(ipart,mydev)%xc = xp
           parts_collect(ipart,mydev)%yc = yp
           parts_collect(ipart,mydev)%zc = zp
           parts_collect(ipart,mydev)%ic = parts(ipart,mydev)%ic
           parts_collect(ipart,mydev)%jc = parts(ipart,mydev)%jc
           parts_collect(ipart,mydev)%kc = parts(ipart,mydev)%kc
           mp =  rho*parts(ipart,mydev)%vol
           parts_collect(ipart,mydev)%hfx = &
              mp*(partforce(1) - (1._rp-rho/rho)*gx )  
           parts_collect(ipart,mydev)%hfy = &
              mp*(partforce(2) - (1._rp-rho/rho)*gy )  
           parts_collect(ipart,mydev)%hfz = &
              mp*(partforce(3) - (1._rp-rho/rho)*gz ) 
           parts_collect(ipart,mydev)%vol = parts(ipart,mydev)%vol
        end if ! TwoWayCouplingFlag
       end do ! ipart

       deallocate( relvel(1:4), partforce(1:3) )

       if(MPI_errorcode.eq.1) then
          call MPI_ABORT(MPI_COMM_WORLD, MPI_errorcode, ierr)
          call MPI_finalize(ierr)
          error stop
       endif

      end if ! num_part(mydev) 
   end subroutine ComputePartPropertiesGPU

   subroutine GetFluidProp(ip,jp,kp,xp,yp,zp,uf,vf,wf,DufDt,DvfDt,DwfDt,volp,TwoWayCouplingFlag,u,v,w,z_c,z_f,n)
      implicit none
      integer, dimension(3), intent(in) :: n
      real(rp), dimension(0:n(1)+1,0:n(2)+1,0:n(3)+1), intent(in) :: u,v,w
      real(rp), dimension(0:n(3)+1), intent(in) :: z_c,z_f
      integer, intent(in)   :: ip,jp,kp
      real(rp), intent(in)  :: xp,yp,zp 
      real(rp), intent(in)  :: volp 
      real(rp), intent(inout) :: uf,vf,wf,DufDt,DvfDt,DwfDt
      integer, intent(in)  :: TwoWayCouplingFlag

      real(rp) :: dp, dxi,dyj,dzk,max_gridsize,min_gridsize

      real(rp) :: ConvertRegSize
      real(rp) :: xh1,xh2,x1,x2,yh1,yh2,y1,y2,zh1,zh2,z1,z2
      integer :: i1,i2,j1,j2,k1,k2

      if ( xp > (n(1)+coord(1)*n(1))*dx+dx .or. xp < ((0+coord(1)*n(1))*dx-dx) ) then 
         write(*,*) ip,jp,kp,xp,yp,zp,((0+coord(1)*n(1))*dx-dx),(n(1)+coord(1)*n(1))*dx
         call lpperror("Failed to get fluid properties as particle outside of x range")
      else if ( yp > (n(2)+coord(2)*n(2))*dy+dy .or. yp < ((0+coord(2)*n(2))*dy-dy) ) then
         write(*,*) ip,jp,kp,xp,yp,zp,((0+coord(2)*n(2))*dy),(n(2)+coord(2)*n(2))*dy
         call lpperror("Failed to get fluid properties as particle outside of y range")
      else if ( zp > z_f(n(3)+1) .or. zp < z_f(0) ) then 
         write(*,*) ip,jp,kp,xp,yp,zp,z_f(0),z_f(n(3))
         call lpperror("Failed to get fluid properties as particle outside of z range")
      end if ! ip,jp,kp

      if ( TwoWayCouplingFlag == TwoWayCouplingIgnore .or. TwoWayCouplingFlag == TwoWayCouplingFilterForce )  then 
           call TrilinearIntrplFluidVel(xp,yp,zp,ip,jp,kp,uf,DufDt,1,u,v,w,z_c,z_f,n)
           call TrilinearIntrplFluidVel(xp,yp,zp,ip,jp,kp,vf,DvfDt,2,u,v,w,z_c,z_f,n)
           call TrilinearIntrplFluidVel(xp,yp,zp,ip,jp,kp,wf,DwfDt,3,u,v,w,z_c,z_f,n)
      else if ( TwoWayCouplingFlag == TwoWayCouplingFilterVel ) then
           call lpperror("Two way coupling with filtered fluid velocity is not ready yet!")
      end if ! dp

   end subroutine GetFluidProp

   subroutine TwoWayCouplingForce(TwoWayCouplingFlag,up,vp,wp,z_c,z_f,dzf,n,ng)
	  implicit none
      integer , intent(in), dimension(3) :: n,ng
      real(rp), dimension(0:n(3)+1), intent(in) :: z_c,z_f,dzf
      real(rp), dimension(0:,0:,0:), intent(inout) :: up,vp,wp
      integer,  intent(in) :: TwoWayCouplingFlag
      integer :: irank
      integer :: MPI_part_collect_type,MPI_part_collect_row, & 
                 oldtypes(0:3), blockcounts(0:3), & 
                 intextent, rpextent
      integer(KIND=MPI_ADDRESS_KIND) :: offsets(0:3)

      integer  :: i,j,k,ipart,ip,jp,kp,nx,ny,nz,ntx,nty,ntz,npx,npy,npz
      real(rp) :: xp,yp,zp,fx,fy,fz,x2,y2,z2,s0,s1,dp
      real(rp) :: massLPP, massVOF, mfLPP
      real(rp) :: mfLPPsmall = 1.e-40_rp
      attributes(managed) :: z_f,dzf,z_c,up,vp,wp

      nx = n(1); ntx = ng(1); npx = coord(1)
      ny = n(2); nty = ng(2); npy = coord(2)
      nz = n(3); ntz = ng(3)

         !  Setup MPI derived type for particle 
         call MPI_TYPE_EXTENT(MPI_REAL_RP, rpextent,  ierr) 
         call MPI_TYPE_EXTENT(MPI_INTEGER, intextent, ierr)

         offsets    (0) = 0 
         oldtypes   (0) = MPI_REAL_RP 
         blockcounts(0) = 7 
         offsets    (1) = offsets(0) + blockcounts(0)*rpextent 
         oldtypes   (1) = MPI_INTEGER  
         blockcounts(1) = 4  
         
         call MPI_TYPE_CREATE_STRUCT(2, blockcounts, offsets, oldtypes, & 
                                     MPI_part_collect_type, ierr) 
         call MPI_TYPE_COMMIT(MPI_part_collect_type, ierr)
         call MPI_TYPE_CONTIGUOUS(max_num_part, MPI_part_collect_type, MPI_part_collect_row, ierr)
         call MPI_TYPE_COMMIT(MPI_part_collect_row, ierr)

         parts_collect_rank(1:max_num_part) = parts_collect(1:max_num_part,mydev)
         call MPI_ALLGATHER(parts_collect_rank(1:max_num_part), 1, MPI_part_collect_row, &
                            parts_collect(1:max_num_part,:),    1, MPI_part_collect_row, & 
                            MPI_COMM_WORLD, ierr)

         ! Compute LPP mass fraction
         massLPP = 0._rp
         if ( num_part(mydev) > 0 ) then
            !$cuf kernel do(1) <<<*,*>>>
            do ipart = 1,num_part(mydev)
               massLPP = massLPP + parts(ipart,mydev)%vol
            end do ! i
            massLPP = massLPP*rho
         end if ! num_part(mydev)
         massVOF = 0._rp
         !$cuf kernel do(1) <<<*,*>>>
         do k=1,nz
             massVOF = massVOF + rho*dx*dy*dzf(k)
         end do

            do irank = 0,(dims(1)*dims(2)-1)
             if ( num_part(irank) > 0 ) then

              !$cuf kernel do(4) <<<*,*>>>
              do ipart = 1,num_part(irank)
               do k=1,nz
                do j=1,ny
                 do i=1,nx
                  mfLPP = (massLPP-parts_collect(ipart,irank)%vol*rho)/massVOF
                  mfLPP = max(mfLPP,small)
                  dp = (6._rp*parts_collect(ipart,irank)%vol/PI)**OneThird 
                  LengthLPP = LengthLPP2dp*dp/erf(mfLPP/mfLPP_ref)
                  s1 = max(LengthLPP,small)
                  if ( i > 0 .and. i < ntx .and. j > 0 .and. j < nty .and. k > 0 .and. k < ntz ) then
                   xp = parts_collect(ipart,irank)%xc
                   yp = parts_collect(ipart,irank)%yc
                   zp = parts_collect(ipart,irank)%zc
                   fx = parts_collect(ipart,irank)%hfx
                   fy = parts_collect(ipart,irank)%hfy
                   fz = parts_collect(ipart,irank)%hfz
                  
                   x2 = ((i+npx*nx)*dx-xp)**2 + (((j+npy*ny)*dy-0.5_rp*dy)-yp)**2 + (z_c(k)-zp)**2
                   y2 = (((i+npx*nx)*dx-0.5_rp*dx)-xp)**2 + ((j+npy*ny)*dy-yp)**2 + (z_c(k)-zp)**2
                   z2 = (((i+npx*nx)*dx-0.5_rp*dx)-xp)**2 + (((j+npy*ny)*dy-0.5_rp*dy)-yp)**2 + (z_f(k)-zp)**2
                   
                   s0 = 1/max((2._rp*PI*s1*s1)**(1.5_rp),small)
                   up(i,j,k) = up(i,j,k) - fx*s0*exp(-x2*0.5_rp/s1/s1)
                   vp(i,j,k) = vp(i,j,k) - fy*s0*exp(-y2*0.5_rp/s1/s1)
                   wp(i,j,k) = wp(i,j,k) - fz*s0*exp(-z2*0.5_rp/s1/s1)
                  end if ! i,j,k
                 end do
                end do
               end do
              end do ! ipart
             endif  ! num_part(irank)
            end do ! irank

   end subroutine TwoWayCouplingForce

   subroutine UpdatePartSol(dt,z_f,dzf,n,ng,istage)
      implicit none
      integer , intent(in), dimension(3) :: n,ng
      real(rp), dimension(0:n(3)+1), intent(in) :: z_f,dzf
      real(rp), intent(in) :: dt
      integer, intent(in) :: istage
      real(rp) :: volcell
      integer :: ipart
      integer :: MPI_errorcode
      attributes(managed) :: z_f,dzf


      ! write(*,*)"\Rank=",myid,"Entered particle location subroutine"

      MPI_errorcode = 0

      if ( num_part(mydev) > 0 ) then
         !$cuf kernel do(1) <<<*,*>>>
         do ipart = 1,num_part(mydev)
            parts(ipart,mydev)%xc = parts(ipart,mydev)%xc + parts(ipart,mydev)%uc*dt 
            parts(ipart,mydev)%yc = parts(ipart,mydev)%yc + parts(ipart,mydev)%vc*dt 
            parts(ipart,mydev)%zc = parts(ipart,mydev)%zc + parts(ipart,mydev)%wc*dt

! TEMPORARY - Check particle location inside domain or not  
            ! if ( parts(ipart,mydev)%xc < -dx .or. parts(ipart,mydev)%xc > (ng(1)+1)*dx ) then
               ! write(*,*) ipart,mydev,parts(ipart,mydev)%xc,& 
                          ! parts(ipart,mydev)%uc,parts(ipart,mydev)%uf, & 
                          ! parts(ipart,mydev)%duc, parts(ipart,mydev)%vol
               ! MPI_errorcode=1
               ! print*, "LPP ERROR *** ","Particle location is out of bound in x direction!"," *** STOP at processor: ", mydev
               ! exit
            ! else if ( parts(ipart,mydev)%yc < -dy .or. parts(ipart,mydev)%yc > (ng(2)+1)*dy ) then
               ! write(*,*) ipart,mydev,parts(ipart,mydev)%yc,& 
                          ! parts(ipart,mydev)%vc,parts(ipart,mydev)%dvc,&
                          ! parts(ipart,mydev)%vol
               ! MPI_errorcode=1
               ! print*, "LPP ERROR *** ","Particle location is out of bound in y direction!"," *** STOP at processor: ", mydev
               ! exit
            ! else if ( parts(ipart,mydev)%zc < -dzf(0) .or. parts(ipart,mydev)%zc > z_f(ng(3)+1) ) then
               ! write(*,*) ipart,mydev,parts(ipart,mydev)%zc,& 
                          ! parts(ipart,mydev)%wc,parts(ipart,mydev)%dwc,&
                          ! parts(ipart,mydev)%vol
               ! MPI_errorcode=1
               ! print*, "LPP ERROR *** ","Particle location is out of bound in z direction!"," *** STOP at processor: ", mydev
               ! exit 
            ! end if ! i1
! END TEMPORARY

            if ( parts(ipart,mydev)%duc /= CRAZY_REAL ) then 
                 parts(ipart,mydev)%uc = parts(ipart,mydev)%uc + &
                                                parts(ipart,mydev)%duc*dt 
                 parts(ipart,mydev)%vc = parts(ipart,mydev)%vc + &
                                                parts(ipart,mydev)%dvc*dt 
                 parts(ipart,mydev)%wc = parts(ipart,mydev)%wc + &
                                                parts(ipart,mydev)%dwc*dt
            end if ! parts(ipart.mydev)%duc
         end do ! ipart
         if(MPI_errorcode.eq.1) then
            call MPI_ABORT(MPI_COMM_WORLD, MPI_errorcode, ierr)
            call MPI_finalize(ierr)
            error stop
         endif
         call UpdatePartLocCell(z_f,dzf,n,ng)
      end if ! num_part(mydev)
   end subroutine UpdatePartSol

   subroutine StorePartOld
      implicit none
      integer :: i

      if ( num_part(mydev) > 0 ) then
      !$cuf kernel do(1) <<<*,*>>>
       do i = 1,num_part(mydev)
           parts(i,mydev)%xcOld = parts(i,mydev)%xc
           parts(i,mydev)%ycOld = parts(i,mydev)%yc
           parts(i,mydev)%zcOld = parts(i,mydev)%zc
           !
           parts(i,mydev)%ucOld = parts(i,mydev)%uc
           parts(i,mydev)%vcOld = parts(i,mydev)%vc
           parts(i,mydev)%wcOld = parts(i,mydev)%wc
       enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! old code !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! parts(1:num_part(mydev),mydev)%xcOld = parts(1:num_part(mydev),mydev)%xc 
         ! parts(1:num_part(mydev),mydev)%ycOld = parts(1:num_part(mydev),mydev)%yc 
         ! parts(1:num_part(mydev),mydev)%zcOld = parts(1:num_part(mydev),mydev)%zc 
   
         ! parts(1:num_part(mydev),mydev)%ucOld = parts(1:num_part(mydev),mydev)%uc 
         ! parts(1:num_part(mydev),mydev)%vcOld = parts(1:num_part(mydev),mydev)%vc 
         ! parts(1:num_part(mydev),mydev)%wcOld = parts(1:num_part(mydev),mydev)%wc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! old code !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end if
      !@cuf istat=cudaDeviceSynchronize()
   end subroutine StorePartOld

   subroutine AveragePartSol(z_f,dzf,n,ng)  ! Update particle solution for 2nd order time integration
      implicit none
      integer :: i
      integer , intent(in), dimension(3) :: n,ng
      real(rp), dimension(0:n(3)+1), intent(in) :: z_f,dzf
      attributes(managed) :: z_f,dzf

       !$cuf kernel do(1) <<<*,*>>>
       do i = 1,num_part(mydev)
           parts(i,mydev)%xc = 0.5_rp*( parts(i,mydev)%xc + parts(i,mydev)%xcOld )
           parts(i,mydev)%yc = 0.5_rp*( parts(i,mydev)%yc + parts(i,mydev)%ycOld )
           parts(i,mydev)%zc = 0.5_rp*( parts(i,mydev)%zc + parts(i,mydev)%zcOld )
           !
           parts(i,mydev)%uc = 0.5_rp*( parts(i,mydev)%uc + parts(i,mydev)%ucOld )
           parts(i,mydev)%vc = 0.5_rp*( parts(i,mydev)%vc + parts(i,mydev)%vcOld )
           parts(i,mydev)%wc = 0.5_rp*( parts(i,mydev)%wc + parts(i,mydev)%wcOld )
       enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! old code !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! parts(1:num_part(mydev),mydev)%xc = 0.5*( parts(1:num_part(mydev),mydev)%xc + parts(1:num_part(mydev),mydev)%xcOld )
     ! parts(1:num_part(mydev),mydev)%yc = 0.5*( parts(1:num_part(mydev),mydev)%yc + parts(1:num_part(mydev),mydev)%ycOld )
     ! parts(1:num_part(mydev),mydev)%zc = 0.5*( parts(1:num_part(mydev),mydev)%zc + parts(1:num_part(mydev),mydev)%zcOld )
     
     ! parts(1:num_part(mydev),mydev)%uc = 0.5*( parts(1:num_part(mydev),mydev)%uc + parts(1:num_part(mydev),mydev)%ucOld )
     ! parts(1:num_part(mydev),mydev)%vc = 0.5*( parts(1:num_part(mydev),mydev)%vc + parts(1:num_part(mydev),mydev)%vcOld )
     ! parts(1:num_part(mydev),mydev)%wc = 0.5*( parts(1:num_part(mydev),mydev)%wc + parts(1:num_part(mydev),mydev)%wcOld )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! old code !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !@cuf istat=cudaDeviceSynchronize()
      call UpdatePartLocCell(z_f,dzf,n,ng)
   end subroutine AveragePartSol

   subroutine boundlpp(n,ng) ! Transfer particles crossing blocks
      implicit none
      integer , intent(in), dimension(3) :: n,ng
      if ( dims(1)*dims(2) > 1 ) then  
         call CollectPartCrossBlocks(n,ng)
         !@cuf istat=cudaDeviceSynchronize()
         call TransferPartCrossBlocks(n,ng)
         !@cuf istat=cudaDeviceSynchronize()
      end if

      call SetPartBC(n,ng)
   end subroutine boundlpp

   subroutine UpdatePartLocCell(z_f,dzf,n,ng)
      implicit none
      integer , intent(in), dimension(3) :: n,ng
      real(rp), dimension(0:n(3)+1), intent(in) :: z_f,dzf
      real(rp) :: xp,yp,zp
      integer :: i,j,k,ipart,i1,j1,k1,kk
      integer :: MPI_errorcode
      logical :: found
      attributes(managed) :: z_f,dzf

      found=.false.
      MPI_errorcode = 0

!     !$cuf kernel do(1) <<<*,*>>>
      ! A fast version for uniform mesh
      do ipart = 1,num_part(mydev)
         ! x direction 
         i  = parts(ipart,mydev)%ic
         xp = parts(ipart,mydev)%xc
         i1 = INT(xp/dx) + 1       ! [BUG] Using intrinsic INT conversion inside CUDA kernel causes numerical errors
         if ( i1 < 0 .or. i1 > ng(1)+1 ) then
          MPI_errorcode=1
          print*, "LPP ERROR *** ","Cell index is out of bound in x direction!"," *** STOP at processor: ", mydev
          write(*,*) i,i1,ipart,xp,dx
          exit
         end if
!         if      ( i <= Ng+1 ) then
!            i1 = INT((xp + real(Ng,rp)*dx(2))/dx(2)) + 1
!         else if ( i > Nx+Ng ) then 
!            i1 = Nx+2*Ng - INT((xh(Nx+2*Ng)-xp)/dx(Nx+Ng+1)) 
!         else if ( xp > xh(i-1) .and. xp <= xh(i  ) ) then
!            i1 = i
!         else if ( xp > xh(i  ) .and. xp <= xh(i+1) ) then 
!            i1 = i+1
!         else if ( xp > xh(i+1) .and. xp <= xh(i+2) ) then 
!            i1 = i+2
!         else if ( xp > xh(i-2) .and. xp <= xh(i-1) ) then 
!            i1 = i-1
!         else if ( xp > xh(i-3) .and. xp <= xh(i-2) ) then 
!            i1 = i-2
!         else 
!            call lpperror("Particle moves out of tracking range in x direction !")
!         end if !xp
        
         ! y direction 
         j  = parts(ipart,mydev)%jc
         yp = parts(ipart,mydev)%yc 
         j1 = INT(yp/dy) + 1
         if ( j1 < 0 .or. j1 > ng(2)+1 ) then
          MPI_errorcode=1
          print*, "LPP ERROR *** ","Cell index is out of bound in y direction!"," *** STOP at processor: ", mydev
          exit
         end if
!         if      ( j <= Ng+1 ) then
!            j1 = INT((yp + real(Ng,rp)*dy(2))/dy(2)) + 1
!         else if ( j > Ny+Ng ) then 
!            j1 = Ny+2*Ng - INT((yh(Ny+2*Ng)-yp)/dy(Ny+Ng+1)) 
!         else if ( yp > yh(j-1) .and. yp <= yh(j  ) ) then
!            j1 = j
!         else if ( yp > yh(j  ) .and. yp <= yh(j+1) ) then 
!            j1 = j+1
!         else if ( yp > yh(j+1) .and. yp <= yh(j+2) ) then 
!            j1 = j+2
!         else if ( yp > yh(j-2) .and. yp <= yh(j-1) ) then 
!            j1 = j-1
!         else if ( yp > yh(j-3) .and. yp <= yh(j-2) ) then 
!            j1 = j-2
!         else
!            call lpperror("Particle moves out of tracking range in y direction !")
!         end if !yp

         ! z direction 
         k  = parts(ipart,mydev)%kc
         zp = parts(ipart,mydev)%zc
         ! k1 = INT((zp/dzf(k)) + 1
         ! if ( k1 < 1 .or. k1 > ng(3) ) call lpperror("Cell index is out of bound in z direction!") 
        if      ( k <= 1 ) then
           k1 = INT(zp/dzf(k)) + 1
           found=.true.
        else if ( k > ng(3) ) then 
           k1 = ng(3)+1 - INT((z_f(ng(3)+1)-zp)/dzf(ng(3)))
           found=.true.
        else if ( zp > z_f(k-1) .and. zp <= z_f(k  ) ) then
           k1 = k
           found=.true.
        else if ( zp > z_f(k  ) .and. zp <= z_f(k+1) ) then 
           k1 = k+1
           found=.true.
        else if ( zp > z_f(k+1) .and. zp <= z_f(k+2) ) then 
           k1 = k+2
           found=.true.
        else if ( zp > z_f(k-2) .and. zp <= z_f(k-1) ) then 
           k1 = k-1
           found=.true.
        else if ( zp > z_f(k-3) .and. zp <= z_f(k-2) ) then 
           k1 = k-2
           found=.true.
        else
         do kk = 1, ng(3)+1
          if ( zp > z_f(kk-1) .and. zp <= z_f(kk) ) then
            k1 = kk
            found=.true.
            exit
          endif
         enddo
        end if !zp
        
        if (.not.found) then
           MPI_errorcode=1
           print*, "LPP ERROR *** ","Particle moves out of tracking range in z direction !"," *** STOP at processor: ", mydev
           write(*,*) k,zp,parts(ipart,mydev)%wc
           write(*,*) z_f(k)
           ! exit
        end if

         parts(ipart,mydev)%ic = i1 
         parts(ipart,mydev)%jc = j1
         parts(ipart,mydev)%kc = k1 

      end do ! ipart
      !@cuf istat=cudaDeviceSynchronize()
      if(MPI_errorcode.eq.1) then
         call MPI_ABORT(MPI_COMM_WORLD, MPI_errorcode, ierr)
         call MPI_finalize(ierr)
         error stop
      endif
   end subroutine UpdatePartLocCell
   
   subroutine CollectPartCrossBlocks(n,ng)
      implicit none
      integer , intent(in), dimension(3) :: n,ng
      integer :: ipart,i,j,k
      integer :: ranknew
      integer :: c1,c2,c3
      integer :: MPI_errorcode
      logical :: PartNeedTransfer, PartExitDomain
      allocate( num_part_cross(0:dims(1)*dims(2)-1) )
      allocate( parts_cross_id(max_num_part_cross,0:dims(1)*dims(2)-1) )
      allocate( parts_cross_id_rank(max_num_part_cross) )
      allocate( parts_cross_newrank(max_num_part_cross,0:dims(1)*dims(2)-1) )
      allocate( parts_cross_newrank_rank(max_num_part_cross) )

      num_part_cross(:) = 0
      parts_cross_id     (:,:) = CRAZY_INT 
      parts_cross_newrank(:,:) = CRAZY_INT
      MPI_errorcode = 0

      if ( num_part(mydev) > 0 ) then
         do ipart = 1,num_part(mydev)
            i = parts(ipart,mydev)%ic
            j = parts(ipart,mydev)%jc
            k = parts(ipart,mydev)%kc

            ! Check if the particle need to be transferred across blocks
            PartNeedTransfer = .false.
            PartExitDomain = .false.
            if ( i < (1+coord(1)*n(1)) ) then
               if ( i > 0 ) then
                  PartNeedTransfer = .true.
               else 
                  if ( lppbdry_cond(1) == 'periodic' ) then 
                     i = i + ng(1)
                     PartNeedTransfer = .true.
                  else if ( lppbdry_cond(1) == 'exit' ) then 
                     PartExitDomain = .true. 
                  end if ! lppbdry
               end if ! i
            else if ( i > (n(1)+coord(1)*n(1)) ) then 
               if ( i <= ng(1) ) then 
                  PartNeedTransfer = .true.
               else 
                  if ( lppbdry_cond(4) == 'periodic' ) then
                     i = i - ng(1)
                     PartNeedTransfer = .true.
                  else if ( lppbdry_cond(4) == 'exit' ) then 
                     PartExitDomain = .true. 
                  end if ! lppbdry
               end if ! i
            end if ! i

            if ( j < (1+coord(2)*n(2)) ) then 
               if ( j > 0 ) then
                  PartNeedTransfer = .true.
               else 
                  if ( lppbdry_cond(2) == 'periodic' ) then 
                     j = j + ng(2)
                     PartNeedTransfer = .true.
                  else if ( lppbdry_cond(2) == 'exit' ) then
                     PartExitDomain = .true. 
                  end if ! lppbdry
               end if ! j
            else if ( j > (n(2)+coord(2)*n(2)) ) then 
               if ( j <= ng(2) ) then 
                  PartNeedTransfer = .true.
               else 
                  if ( lppbdry_cond(5) == 'periodic' ) then 
                     j = j - ng(2)
                     PartNeedTransfer = .true.
                  else if ( lppbdry_cond(5) == 'exit' ) then
                     PartExitDomain = .true. 
                  end if ! lppbdry
               end if ! j
            end if ! j

            if ( k < 1 ) then 
               if ( k > 0 ) then 
                  PartNeedTransfer = .true.
               else 
                  if ( lppbdry_cond(3) == 'periodic' ) then 
                     k = k + ng(3)
                     PartNeedTransfer = .true.
                  else if ( lppbdry_cond(3) == 'exit' ) then
                     PartExitDomain = .true. 
                  end if ! lppbdry
               end if ! k
            else if ( k > n(3) ) then 
               if ( k <= ng(3) ) then
                  PartNeedTransfer = .true.
               else 
                  if ( lppbdry_cond(6) == 'periodic' ) then 
                     k = k - ng(3)
                     PartNeedTransfer = .true.
                  else if ( lppbdry_cond(6) == 'exit' ) then
                     PartExitDomain = .true. 
                  end if ! lppbdry
               end if ! k 
            end if ! k

            ! if particle locates at boundary cell in any one direction and the
            ! BC is not periodic, then no need to transfer particle 
            if ( PartExitDomain ) PartNeedTransfer = .false.

            ! Note: here only collect and transfer particles which cross blocks 
            !       due to periodic BC, the location and cell index stored in 
            !       parts(ipart,mydev) will not be changed until SetPartBC is called

            ! Compute new rank for the particle if transfer is needed
            if ( PartNeedTransfer ) then  
               c1 = (i-1)/n(1)  
               c2 = (j-1)/n(2)
               c3 = (k-1)/n(3)
               ranknew = c1*dims(2)*1 + c2*1 + c3
               if ( ranknew > ((dims(1)*dims(2))-1) .or. ranknew < 0 ) then
                  MPI_errorcode=1
                  print*, "LPP ERROR *** ","new rank of particle out of range!"," *** STOP at processor: ", mydev
                  exit
               else if ( ranknew /= mydev ) then 
                  num_part_cross(mydev)  = num_part_cross(mydev) + 1
                  parts_cross_id     (num_part_cross(mydev),mydev) = ipart 
                  parts_cross_newrank(num_part_cross(mydev),mydev) = ranknew 
               end if ! ranknew
            end if ! PartNeedTransfer

         end do ! ipart
         if(MPI_errorcode.eq.1) then
            call MPI_ABORT(MPI_COMM_WORLD, MPI_errorcode, ierr)
            call MPI_finalize(ierr)
            error stop
         endif
      end if ! num_part(mydev)
   end subroutine CollectPartCrossBlocks

   subroutine TransferPartCrossBlocks(n,ng)
	  implicit none
      integer , intent(in), dimension(3) :: n,ng
      integer :: ipart,ipart_cross,ipart1,i,j
      integer :: irank
      integer :: ranknew
      integer :: req(24),sta(MPI_STATUS_SIZE,24)
      integer :: MPI_particle_type, oldtypes(0:3), blockcounts(0:3), & 
                 intextent, rpextent
      integer(KIND=MPI_ADDRESS_KIND) :: offsets(0:3)
      integer :: max_num_part_cross_use, MPI_int_row
      type(particle), managed :: element_NULL
      integer :: num_part_cross_rank,num_part_rank
   
      element_NULL%xc = 0._rp;element_NULL%yc = 0._rp;element_NULL%zc = 0._rp
      element_NULL%uc = 0._rp;element_NULL%vc = 0._rp;element_NULL%wc = 0._rp
      element_NULL%duc = 0._rp;element_NULL%dvc = 0._rp;element_NULL%dwc = 0._rp
      element_NULL%vol = 0._rp;element_NULL%id = CRAZY_INT

      num_part_cross_rank = num_part_cross(mydev)
      call MPI_ALLGATHER(num_part_cross_rank, 1, MPI_INTEGER, &
                         num_part_cross,      1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
      max_num_part_cross_use = maxval(num_part_cross)
      if ( max_num_part_cross_use > max_num_part_cross ) call lpperror("Increase max_num_part_cross")

      if ( max_num_part_cross_use  > 0 ) then

         call MPI_TYPE_CONTIGUOUS(max_num_part_cross_use, MPI_INTEGER, MPI_int_row, ierr)
         call MPI_TYPE_COMMIT(MPI_int_row, ierr)

       !$cuf kernel do(1) <<<*,*>>>
       do i = 1,max_num_part_cross_use
        cross_buff_id_rank(i) = parts_cross_id(i,mydev)
        cross_buff_new_rank(i) = parts_cross_newrank(i,mydev)
       enddo
       !@cuf istat = cudaDeviceSynchronize()
       
            call MPI_ALLGATHER(cross_buff_id_rank(1:max_num_part_cross_use), 1, MPI_int_row, &
                               cross_buff_id(1:max_num_part_cross_use,:),    1, MPI_int_row, & 
                               MPI_COMM_WORLD, ierr)
            call MPI_ALLGATHER(cross_buff_new_rank(1:max_num_part_cross_use), 1, MPI_int_row, &
                               cross_buff_new(1:max_num_part_cross_use,:),    1, MPI_int_row, & 
                               MPI_COMM_WORLD, ierr)
                               
       !$cuf kernel do(2) <<<*,*>>>
       do i = 1,max_num_part_cross_use
        do j = 0,(dims(1)*dims(2))-1
         parts_cross_id(i,j) = cross_buff_id(i,j)
         parts_cross_newrank(i,j) = cross_buff_new(i,j)
        enddo
       enddo
    !@cuf istat = cudaDeviceSynchronize()

       !  Setup MPI derived type for particle 
       call MPI_TYPE_EXTENT(MPI_REAL_RP, rpextent,  ierr) 
       call MPI_TYPE_EXTENT(MPI_INTEGER, intextent, ierr) 
       offsets    (0) = 0 
       oldtypes   (0) = MPI_REAL_RP 
       blockcounts(0) = 11
       offsets    (1) = offsets(0) + blockcounts(0)*rpextent 
       oldtypes   (1) = MPI_INTEGER  
       blockcounts(1) = 1  
       offsets    (2) = offsets(1) + blockcounts(1)*intextent 
       oldtypes   (2) = MPI_REAL_RP  
       blockcounts(2) = 9
       offsets    (3) = offsets(2) + blockcounts(2)*rpextent
       oldtypes   (3) = MPI_INTEGER  
       blockcounts(3) = 5  
       
       call MPI_TYPE_CREATE_STRUCT(4, blockcounts, offsets, oldtypes, & 
                            MPI_particle_type, ierr) 
       call MPI_TYPE_COMMIT(MPI_particle_type, ierr)

       ! Store particle information to be exchanged in GPU buffer arrays
       do irank = 0,(dims(1)*dims(2)-1)
        if ( num_part_cross(irank) > 0 ) then
         !$cuf kernel do(1) <<<*,*>>>
         do ipart_cross = 1,num_part_cross(irank)
            ipart   = parts_cross_id     (ipart_cross,irank)
            ranknew = parts_cross_newrank(ipart_cross,irank)
            if ( mydev == irank ) then
             xcOld_s_buf(ipart,irank) = parts(ipart,irank)%xcOld
             ycOld_s_buf(ipart,irank) = parts(ipart,irank)%ycOld
             zcOld_s_buf(ipart,irank) = parts(ipart,irank)%zcOld
             ucOld_s_buf(ipart,irank) = parts(ipart,irank)%ucOld
             vcOld_s_buf(ipart,irank) = parts(ipart,irank)%vcOld
             wcOld_s_buf(ipart,irank) = parts(ipart,irank)%wcOld
             uf_s_buf(ipart,irank) = parts(ipart,irank)%uf
             vf_s_buf(ipart,irank) = parts(ipart,irank)%vf
             wf_s_buf(ipart,irank) = parts(ipart,irank)%wf
             uf_s_buf(ipart,irank) = parts(ipart,irank)%uc
             vf_s_buf(ipart,irank) = parts(ipart,irank)%vc
             wf_s_buf(ipart,irank) = parts(ipart,irank)%wc
             ic_s_buf(ipart,irank) = parts(ipart,irank)%ic
             jc_s_buf(ipart,irank) = parts(ipart,irank)%jc
             kc_s_buf(ipart,irank) = parts(ipart,irank)%kc
            else if ( mydev == ranknew ) then
             num_part(ranknew) = num_part(ranknew) + 1
            endif
         enddo
        endif
       enddo
    
       ! Transfer particle information to buffers in destination block
       do irank = 0,(dims(1)*dims(2)-1)
         ranknew = parts_cross_newrank(ipart_cross,irank)
          if ( num_part_cross(irank) > 0 ) then
                if ( mydev == irank ) then
                   call MPI_ISEND(xcOld_s_buf(:,irank),size(xcOld_s_buf(:,irank)), MPI_REAL_RP, & 
                                  ranknew, 15, MPI_COMM_WORLD, req(1), ierr)
                   call MPI_WAIT(req(1),sta(:,1),ierr)
                   !
                   call MPI_ISEND(ycOld_s_buf(:,irank),size(ycOld_s_buf(:,irank)), MPI_REAL_RP, & 
                                  ranknew, 16, MPI_COMM_WORLD, req(2), ierr)
                   call MPI_WAIT(req(2),sta(:,2),ierr)
                   !
                   call MPI_ISEND(zcOld_s_buf(:,irank),size(zcOld_s_buf(:,irank)), MPI_REAL_RP, & 
                                  ranknew, 17, MPI_COMM_WORLD, req(3), ierr)
                   call MPI_WAIT(req(3),sta(:,3),ierr)
                   !
                   call MPI_ISEND(ucOld_s_buf(:,irank),size(ucOld_s_buf(:,irank)), MPI_REAL_RP, & 
                                  ranknew, 18, MPI_COMM_WORLD, req(4), ierr)
                   call MPI_WAIT(req(4),sta(:,4),ierr)
                   !
                   call MPI_ISEND(vcOld_s_buf(:,irank),size(vcOld_s_buf(:,irank)), MPI_REAL_RP, & 
                                  ranknew, 19, MPI_COMM_WORLD, req(5), ierr)
                   call MPI_WAIT(req(5),sta(:,5),ierr)
                   !
                   call MPI_ISEND(wcOld_s_buf(:,irank),size(wcOld_s_buf(:,irank)), MPI_REAL_RP, & 
                                  ranknew, 20, MPI_COMM_WORLD, req(6), ierr)
                   call MPI_WAIT(req(6),sta(:,6),ierr)
                   !
                   call MPI_ISEND(uf_s_buf(:,irank),size(uf_s_buf(:,irank)), MPI_REAL_RP, & 
                                  ranknew, 21, MPI_COMM_WORLD, req(7), ierr)
                   call MPI_WAIT(req(7),sta(:,7),ierr)
                   !
                   call MPI_ISEND(vf_s_buf(:,irank),size(vf_s_buf(:,irank)), MPI_REAL_RP, & 
                                  ranknew, 22, MPI_COMM_WORLD, req(8), ierr)
                   call MPI_WAIT(req(8),sta(:,8),ierr)
                   !
                   call MPI_ISEND(wf_s_buf(:,irank),size(wf_s_buf(:,irank)), MPI_REAL_RP, & 
                                  ranknew, 23, MPI_COMM_WORLD, req(9), ierr)
                   call MPI_WAIT(req(9),sta(:,9),ierr)
                   !
                   call MPI_ISEND(ic_s_buf(:,irank),size(ic_s_buf(:,irank)), MPI_INTEGER, & 
                                  ranknew, 24, MPI_COMM_WORLD, req(10), ierr)
                   call MPI_WAIT(req(10),sta(:,10),ierr)
                   !
                   call MPI_ISEND(jc_s_buf(:,irank),size(jc_s_buf(:,irank)), MPI_INTEGER, & 
                                  ranknew, 25, MPI_COMM_WORLD, req(11), ierr)
                   call MPI_WAIT(req(11),sta(:,11),ierr)
                   !
                   call MPI_ISEND(kc_s_buf(:,irank),size(kc_s_buf(:,irank)), MPI_INTEGER, & 
                                  ranknew, 26, MPI_COMM_WORLD, req(12), ierr)
                   call MPI_WAIT(req(12),sta(:,12),ierr)
                else if ( mydev == ranknew ) then 
                   call MPI_IRECV(xcOld_r_buf(:,ranknew),size(xcOld_r_buf(:,ranknew)),MPI_REAL_RP, & 
                                  irank, 15, MPI_COMM_WORLD, req(13), ierr)
                   call MPI_WAIT(req(13),sta(:,13),ierr)
                   !
                   call MPI_IRECV(ycOld_r_buf(:,ranknew),size(ycOld_r_buf(:,ranknew)),MPI_REAL_RP, & 
                                  irank, 16, MPI_COMM_WORLD, req(14), ierr)
                   call MPI_WAIT(req(14),sta(:,14),ierr)
                   !
                   call MPI_IRECV(zcOld_r_buf(:,ranknew),size(zcOld_r_buf(:,ranknew)),MPI_REAL_RP, & 
                                  irank, 17, MPI_COMM_WORLD, req(15), ierr)
                   call MPI_WAIT(req(15),sta(:,15),ierr)
                   !
                   call MPI_IRECV(ucOld_r_buf(:,ranknew),size(ucOld_r_buf(:,ranknew)),MPI_REAL_RP, & 
                                  irank, 18, MPI_COMM_WORLD, req(16), ierr)
                   call MPI_WAIT(req(16),sta(:,16),ierr)
                   !
                   call MPI_IRECV(vcOld_r_buf(:,ranknew),size(vcOld_r_buf(:,ranknew)),MPI_REAL_RP, & 
                                  irank, 19, MPI_COMM_WORLD, req(17), ierr)
                   call MPI_WAIT(req(17),sta(:,17),ierr)
                   !
                   call MPI_IRECV(wcOld_r_buf(:,ranknew),size(wcOld_r_buf(:,ranknew)),MPI_REAL_RP, & 
                                  irank, 20, MPI_COMM_WORLD, req(18), ierr)
                   call MPI_WAIT(req(18),sta(:,18),ierr)
                   !
                   call MPI_IRECV(uf_r_buf(:,ranknew),size(uf_r_buf(:,ranknew)),MPI_REAL_RP, & 
                                  irank, 21, MPI_COMM_WORLD, req(19), ierr)
                   call MPI_WAIT(req(19),sta(:,19),ierr)
                   !
                   call MPI_IRECV(vf_r_buf(:,ranknew),size(vf_r_buf(:,ranknew)),MPI_REAL_RP, & 
                                  irank, 22, MPI_COMM_WORLD, req(20), ierr)
                   call MPI_WAIT(req(20),sta(:,20),ierr)
                   !
                   call MPI_IRECV(wf_r_buf(:,ranknew),size(wf_r_buf(:,ranknew)),MPI_REAL_RP, & 
                                  irank, 23, MPI_COMM_WORLD, req(21), ierr)
                   call MPI_WAIT(req(21),sta(:,21),ierr)
                   !
                   call MPI_IRECV(ic_r_buf(:,ranknew),size(ic_r_buf(:,ranknew)),MPI_INTEGER, & 
                                  irank, 24, MPI_COMM_WORLD, req(22), ierr)
                   call MPI_WAIT(req(22),sta(:,22),ierr)
                   !
                   call MPI_IRECV(jc_r_buf(:,ranknew),size(jc_r_buf(:,ranknew)),MPI_INTEGER, & 
                                  irank, 25, MPI_COMM_WORLD, req(23), ierr)
                   call MPI_WAIT(req(23),sta(:,23),ierr)
                   !
                   call MPI_IRECV(kc_r_buf(:,ranknew),size(kc_r_buf(:,ranknew)),MPI_INTEGER, & 
                                  irank, 26, MPI_COMM_WORLD, req(24), ierr)
                   call MPI_WAIT(req(24),sta(:,24),ierr)
                end if ! mydev
          end if ! num_part_cross(irank)
       end do ! irank
      !@cuf istat = cudaDeviceSynchronize()

       ! Update new particle info in destination block
       do irank = 0,(dims(1)*dims(2)-1)
        if ( num_part_cross(irank) > 0 ) then
         !$cuf kernel do(1) <<<*,*>>>
         do ipart_cross = 1,num_part_cross(irank)
            ipart   = parts_cross_id     (ipart_cross,irank)
            ranknew = parts_cross_newrank(ipart_cross,irank)
            if ( mydev == ranknew ) then
             parts(num_part(ranknew),ranknew)%xcOld = xcOld_r_buf(ipart,ranknew)
             parts(num_part(ranknew),ranknew)%ycOld = ycOld_r_buf(ipart,ranknew) 
             parts(num_part(ranknew),ranknew)%zcOld = zcOld_r_buf(ipart,ranknew) 
             parts(num_part(ranknew),ranknew)%ucOld = ucOld_r_buf(ipart,ranknew) 
             parts(num_part(ranknew),ranknew)%vcOld = vcOld_r_buf(ipart,ranknew) 
             parts(num_part(ranknew),ranknew)%wcOld = wcOld_r_buf(ipart,ranknew) 
             parts(num_part(ranknew),ranknew)%uf = uf_r_buf(ipart,ranknew)
             parts(num_part(ranknew),ranknew)%vf = vf_r_buf(ipart,ranknew)
             parts(num_part(ranknew),ranknew)%wf = wf_r_buf(ipart,ranknew)
             parts(num_part(ranknew),ranknew)%uc = uf_r_buf(ipart,ranknew)
             parts(num_part(ranknew),ranknew)%vc = vf_r_buf(ipart,ranknew)
             parts(num_part(ranknew),ranknew)%wc = wf_r_buf(ipart,ranknew)
             parts(num_part(ranknew),ranknew)%ic = ic_r_buf(ipart,ranknew)
             parts(num_part(ranknew),ranknew)%jc = jc_r_buf(ipart,ranknew)
             parts(num_part(ranknew),ranknew)%kc = kc_r_buf(ipart,ranknew)
            endif
         enddo
        endif
       enddo

      ! Remove particle in the origin block
      if ( num_part_cross(mydev) > 0 ) then
         do ipart_cross = 1,num_part_cross(mydev)
            ipart = parts_cross_id(ipart_cross,mydev)
            do ipart1 = ipart,num_part(mydev)-1
               parts(ipart1,mydev) = parts(ipart1+1,mydev)
            end do ! ipart1
            parts(num_part(mydev),mydev) = element_NULL
            num_part(mydev) = num_part(mydev) - 1
            if ( ipart_cross < num_part_cross(mydev) ) & 
               parts_cross_id(ipart_cross+1:num_part_cross(mydev),mydev) = &
               parts_cross_id(ipart_cross+1:num_part_cross(mydev),mydev) - 1
         end do ! ipart_cross 
      end if ! num_part_cross(irank)

      ! Allgather updated number of particles in every block
      num_part_rank = num_part(mydev)
      call MPI_ALLGATHER(num_part_rank, 1, MPI_INTEGER, &
                         num_part(:)  , 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
      end if ! max_num_part_cross_use

      ! final
      deallocate(num_part_cross)
      deallocate(parts_cross_id)
      deallocate(parts_cross_id_rank)
      deallocate(parts_cross_newrank)
      deallocate(parts_cross_newrank_rank)

   end subroutine TransferPartCrossBlocks

   subroutine SetPartBC(n,ng)
	  implicit none
      integer , intent(in), dimension(3) :: n,ng
      integer :: ipart,ipart1
      type(particle), managed :: element_NULL
      integer :: num_part_rank
   
      element_NULL%xc = 0._rp;element_NULL%yc = 0._rp;element_NULL%zc = 0._rp
      element_NULL%uc = 0._rp;element_NULL%vc = 0._rp;element_NULL%wc = 0._rp
      element_NULL%duc = 0._rp;element_NULL%dvc = 0._rp;element_NULL%dwc = 0._rp
      element_NULL%vol = 0._rp;element_NULL%id = CRAZY_INT

      if ( num_part(mydev) > 0 ) then
         do ipart = 1,num_part(mydev)
            if ( parts(ipart,mydev)%ic <= 0 ) then
               if ( lppbdry_cond(1) == 'exit') then 
                  call PartBC_exit(ipart,mydev)
                  cycle 
               else 
                  call ImposePartBC(ipart,mydev,1,n,ng)
               end if ! lppbdry_cond
            else if ( parts(ipart,mydev)%ic > ng(1) ) then 
               if ( lppbdry_cond(4) == 'exit') then 
                  parts(ipart,mydev)%vol = -1._rp
                  cycle 
               else 
                  call ImposePartBC(ipart,mydev,4,n,ng)
               end if ! lppbdry_cond
            end if ! ic

            if ( parts(ipart,mydev)%jc <= 0 ) then 
               if ( lppbdry_cond(2) == 'exit') then 
                  parts(ipart,mydev)%vol = -1._rp
                  cycle 
               else 
                  call ImposePartBC(ipart,mydev,2,n,ng)
               end if ! lppbdry_cond
            else if ( parts(ipart,mydev)%jc > ng(2) ) then 
               if ( lppbdry_cond(5) == 'exit') then 
                  parts(ipart,mydev)%vol = -1._rp
                  cycle 
               else 
                  call ImposePartBC(ipart,mydev,5,n,ng)
               end if ! lppbdry_cond
            end if ! jc 

            if ( parts(ipart,mydev)%kc <= 0 ) then  
               if ( lppbdry_cond(3) == 'exit') then 
                  parts(ipart,mydev)%vol = -1._rp
                  cycle 
               else 
                  call ImposePartBC(ipart,mydev,3,n,ng)
               end if ! lppbdry_cond
            else if ( parts(ipart,mydev)%kc > ng(3) ) then 
               if ( lppbdry_cond(6) == 'exit') then 
                  parts(ipart,mydev)%vol = -1._rp
                  cycle 
               else 
                  call ImposePartBC(ipart,mydev,6,n,ng)
               end if ! lppbdry_cond
            end if ! kc
         
         end do ! ipart

         ! Remove particles that exit the domain from the list (vol=-1)
         ipart = 0 
         do
            ipart = ipart + 1
            if ( parts(ipart,mydev)%vol < 0._rp ) then
               if ( ipart < num_part(mydev) ) then
                  do ipart1 = ipart+1,num_part(mydev)
                     parts(ipart1-1,mydev) = parts(ipart1,mydev)
                  end do ! ipart1
               end if ! ipart
               parts(num_part(mydev),mydev) = element_NULL
               num_part(mydev) = num_part(mydev) - 1
               ipart          = ipart          - 1 ! Note: the index also needs to be shifted
            end if ! parts(ipart,mydev)%vol
            if ( ipart >= num_part(mydev) ) exit
         end do ! ipart

      end if ! num_part(mydev)

      ! Allgather updated number of particles after imposing BC
      num_part_rank = num_part(mydev)
      call MPI_ALLGATHER(num_part_rank, 1, MPI_INTEGER, &
                         num_part(:)  , 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)

   end subroutine SetPartBC

   subroutine ImposePartBC(ipart,mydev,d,n,ng)
      implicit none
      integer , intent(in), dimension(3) :: n,ng
      integer, intent (in) :: ipart,mydev,d

      if ( lppbdry_cond(d) == 'periodic' ) then
         call PartBC_periodic(ipart,mydev,d,n,ng)
      else if ( lppbdry_cond(d) == 'reflect' ) then
         call PartBC_reflect(ipart,mydev,d,n,ng)
      else 
         call lpperror("Unknown particle boundary condition!")
      end if ! lppbdry_cond
   end subroutine ImposePartBC

   subroutine PartBC_periodic(ipart,mydev,d,n,ng)
      implicit none
      integer , intent(in), dimension(3) :: n,ng
      integer, intent (in) :: ipart,mydev,d
      
      if ( d == 1 ) then 
         parts(ipart,mydev)%ic = parts(ipart,mydev)%ic + ng(1)
         parts(ipart,mydev)%xc = parts(ipart,mydev)%xc + lx
      else if ( d == 4 ) then 
         parts(ipart,mydev)%ic = parts(ipart,mydev)%ic - ng(1)
         parts(ipart,mydev)%xc = parts(ipart,mydev)%xc - lx
      else if ( d == 2 ) then 
         parts(ipart,mydev)%jc = parts(ipart,mydev)%jc + ng(2)
         parts(ipart,mydev)%yc = parts(ipart,mydev)%yc + ly
      else if ( d == 5 ) then 
         parts(ipart,mydev)%jc = parts(ipart,mydev)%jc - ng(2)
         parts(ipart,mydev)%yc = parts(ipart,mydev)%yc - ly
      else if ( d == 3 ) then 
         parts(ipart,mydev)%kc = parts(ipart,mydev)%kc + ng(3)
         parts(ipart,mydev)%zc = parts(ipart,mydev)%zc + lz
      else if ( d == 6 ) then 
         parts(ipart,mydev)%kc = parts(ipart,mydev)%kc - ng(3)
         parts(ipart,mydev)%zc = parts(ipart,mydev)%zc - lz
      end if ! d
   end subroutine PartBC_periodic

   subroutine PartBC_exit(ipart,mydev)
      implicit none

      integer, intent (in) :: ipart,mydev
      
      ! For particle that will exit the domain, mark volume=-1 
      parts(ipart,mydev)%vol = -1._rp
   end subroutine PartBC_exit

   subroutine PartBC_reflect(ipart,mydev,d,n,ng)
      implicit none
      integer , intent(in), dimension(3) :: n,ng
      integer, intent (in) :: ipart,mydev,d
      
      if ( d == 1 ) then 
         parts(ipart,mydev)%ic = 2* 0    +1-parts(ipart,mydev)%ic 
         parts(ipart,mydev)%xc = -parts(ipart,mydev)%xc 
         parts(ipart,mydev)%uc = -parts(ipart,mydev)%uc 
      else if ( d == 4 ) then 
         parts(ipart,mydev)%ic = 2*(ng(1)+0)+1-parts(ipart,mydev)%ic 
         parts(ipart,mydev)%xc = 2._rp*lx - parts(ipart,mydev)%xc 
         parts(ipart,mydev)%uc = -parts(ipart,mydev)%uc 
      else if ( d == 2 ) then 
         parts(ipart,mydev)%jc = 2* 0    +1-parts(ipart,mydev)%jc 
         parts(ipart,mydev)%yc = -parts(ipart,mydev)%yc 
         parts(ipart,mydev)%vc = -parts(ipart,mydev)%vc 
      else if ( d == 5 ) then 
         parts(ipart,mydev)%jc = 2*(ng(2)+0)+1-parts(ipart,mydev)%jc 
         parts(ipart,mydev)%yc = 2._rp*ly - parts(ipart,mydev)%yc 
         parts(ipart,mydev)%vc = -parts(ipart,mydev)%vc 
      else if ( d == 3 ) then 
         parts(ipart,mydev)%kc = 2* 0    +1-parts(ipart,mydev)%kc 
         parts(ipart,mydev)%zc = -parts(ipart,mydev)%zc 
         parts(ipart,mydev)%wc = -parts(ipart,mydev)%wc 
      else if ( d == 6 ) then 
         parts(ipart,mydev)%kc = 2*(ng(3)+0)+1-parts(ipart,mydev)%kc 
         parts(ipart,mydev)%zc = 2._rp*lz - parts(ipart,mydev)%zc 
         parts(ipart,mydev)%wc = -parts(ipart,mydev)%wc 
      end if ! d
   end subroutine PartBC_reflect
   
   subroutine LinearIntrpl(x,x0,x1,f0,f1,f)
      implicit none
      real(rp), intent (in) :: x,x0,x1,f0,f1
      real(rp), intent(inout) :: f      
      real(rp) :: xl,xr

      xl = (x-x0)/(x1-x0)
      xr = 1._rp - xl
      f  = f0*xr + f1*xl
   end subroutine LinearIntrpl

   subroutine BilinearIntrpl(x,y,x0,y0,x1,y1,f00,f01,f10,f11,f)
      implicit none
      real(rp), intent (in) :: x,y,x0,y0,x1,y1,f00,f01,f10,f11
      real(rp), intent(inout) :: f      
      real(rp) :: f0,f1

      call LinearIntrpl(x,x0,x1,f00,f10,f0)
      call LinearIntrpl(x,x0,x1,f01,f11,f1)
      call LinearIntrpl(y,y0,y1,f0 ,f1 ,f)
   end subroutine BilinearIntrpl

   subroutine TrilinearIntrpl(x,y,z,x0,y0,z0,x1,y1,z1,f000,f001,f010,f011,f100,f101,f110,f111,f)
      implicit none
      real(rp), intent (in) :: x,y,z,x0,y0,z0,x1,y1,z1, & 
                               f000,f001,f010,f011,f100,f101,f110,f111
      real(rp), intent(inout) :: f
      real(rp) :: f0,f1,f00,f01,f10,f11

      call LinearIntrpl(x,x0,x1,f000,f100,f00)
      call LinearIntrpl(x,x0,x1,f010,f110,f10)
      call LinearIntrpl(x,x0,x1,f001,f101,f01)
      call LinearIntrpl(x,x0,x1,f011,f111,f11)
      call LinearIntrpl(y,y0,y1,f00,f10,f0)
      call LinearIntrpl(y,y0,y1,f01,f11,f1)
      call LinearIntrpl(z,z0,z1,f0 ,f1 ,f)
   end subroutine TrilinearIntrpl

   subroutine TrilinearIntrplFluidVel(xp,yp,zp,ip,jp,kp,vel,sdvel,dir,u,v,w,zc,zf,n)
      implicit none
      integer , intent(in), dimension(3) :: n
      real(rp), dimension(0:n(1)+1,0:n(2)+1,0:n(3)+1), intent(in) :: u,v,w
      real(rp), dimension(0:n(3)+1), intent(in) :: zc,zf
      real(rp), intent (in) :: xp,yp,zp
      integer , intent (in) :: ip,jp,kp,dir
      real(rp), intent(inout) :: vel,sdvel
      real(rp) :: xh1,xh2,yh1,yh2
      integer :: si,sj,sk

      if ( dir == 1 ) then ! Trilinear interpolation for u
         si = -1;sj = -1;sk = -1

         if ( yp > (jp*dy-0.5_rp*dy) ) sj =  0
         if ( zp > zc(kp) ) sk =  0

         xh1 = (ip+si)*dx; xh2 = (ip+1+si)*dx
         yh1 = (jp+sj)*dy; yh2 = (jp+1+sj)*dy

         call TrilinearIntrpl(xp,yp,zp,xh1,yh1,zf(kp  +sk),        &
                                       xh2,yh2,zf(kp+1+sk),        &
                                       u(ip  +si,jp  +sj,kp  +sk), &
                                       u(ip  +si,jp  +sj,kp+1+sk), &
                                       u(ip  +si,jp+1+sj,kp  +sk), &
                                       u(ip  +si,jp+1+sj,kp+1+sk), &
                                       u(ip+1+si,jp  +sj,kp  +sk), &
                                       u(ip+1+si,jp  +sj,kp+1+sk), &
                                       u(ip+1+si,jp+1+sj,kp  +sk), &
                                       u(ip+1+si,jp+1+sj,kp+1+sk), vel)

         if ( UnsteadyPartForce ) &
         call TrilinearIntrpl(xp,yp,zp,xh1,yh1,zf(kp  +sk),          &
                                       xh2,yh2,zf(kp+1+sk),          &
                                       sdu(ip  +si,jp  +sj,kp  +sk), &
                                       sdu(ip  +si,jp  +sj,kp+1+sk), &
                                       sdu(ip  +si,jp+1+sj,kp  +sk), &
                                       sdu(ip  +si,jp+1+sj,kp+1+sk), &
                                       sdu(ip+1+si,jp  +sj,kp  +sk), &
                                       sdu(ip+1+si,jp  +sj,kp+1+sk), &
                                       sdu(ip+1+si,jp+1+sj,kp  +sk), &
                                       sdu(ip+1+si,jp+1+sj,kp+1+sk), sdvel)

      else if ( dir == 2 ) then ! Trilinear interpolation for v
         si = -1;sj = -1;sk = -1

         if ( xp > (ip*dx-0.5_rp*dx) ) si =  0 
         if ( zp > zc(kp) ) sk =  0

         xh1 = (ip+si)*dx; xh2 = (ip+1+si)*dx
         yh1 = (jp+sj)*dy; yh2 = (jp+1+sj)*dy

         call TrilinearIntrpl(xp,yp,zp,xh1,yh1,zf(kp  +sk),        &
                                       xh2,yh2,zf(kp+1+sk),        &
                                       v(ip  +si,jp  +sj,kp  +sk), &
                                       v(ip  +si,jp  +sj,kp+1+sk), &
                                       v(ip  +si,jp+1+sj,kp  +sk), &
                                       v(ip  +si,jp+1+sj,kp+1+sk), &
                                       v(ip+1+si,jp  +sj,kp  +sk), &
                                       v(ip+1+si,jp  +sj,kp+1+sk), &
                                       v(ip+1+si,jp+1+sj,kp  +sk), &
                                       v(ip+1+si,jp+1+sj,kp+1+sk), vel)

         if ( UnsteadyPartForce ) & 
         call TrilinearIntrpl(xp,yp,zp,xh1,yh1,zf(kp  +sk),          &
                                       xh2,yh2,zf(kp+1+sk),          &
                                       sdv(ip  +si,jp  +sj,kp  +sk), &
                                       sdv(ip  +si,jp  +sj,kp+1+sk), &
                                       sdv(ip  +si,jp+1+sj,kp  +sk), &
                                       sdv(ip  +si,jp+1+sj,kp+1+sk), &
                                       sdv(ip+1+si,jp  +sj,kp  +sk), &
                                       sdv(ip+1+si,jp  +sj,kp+1+sk), &
                                       sdv(ip+1+si,jp+1+sj,kp  +sk), &
                                       sdv(ip+1+si,jp+1+sj,kp+1+sk), sdvel)

      else if ( dir == 3 ) then ! Trilinear interpolation for w
         si = -1;sj = -1;sk = -1

         if ( xp > (ip*dx-0.5_rp*dx) ) si =  0
         if ( yp > (jp*dy-0.5_rp*dy) ) sj =  0

         xh1 = (ip+si)*dx; xh2 = (ip+1+si)*dx
         yh1 = (jp+sj)*dy; yh2 = (jp+1+sj)*dy

         call TrilinearIntrpl(xp,yp,zp,xh1,yh1,zf(kp  +sk),        & 
                                       xh2,yh2,zf(kp+1+sk),        &
                                       w(ip  +si,jp  +sj,kp  +sk), &
                                       w(ip  +si,jp  +sj,kp+1+sk), &
                                       w(ip  +si,jp+1+sj,kp  +sk), &
                                       w(ip  +si,jp+1+sj,kp+1+sk), &
                                       w(ip+1+si,jp  +sj,kp  +sk), &
                                       w(ip+1+si,jp  +sj,kp+1+sk), &
                                       w(ip+1+si,jp+1+sj,kp  +sk), &
                                       w(ip+1+si,jp+1+sj,kp+1+sk), vel)

         if ( UnsteadyPartForce ) & 
         call TrilinearIntrpl(xp,yp,zp,xh1,yh1,zf(kp  +sk),          &
                                       xh2,yh2,zf(kp+1+sk),          &
                                       sdw(ip  +si,jp  +sj,kp  +sk), &
                                       sdw(ip  +si,jp  +sj,kp+1+sk), &
                                       sdw(ip  +si,jp+1+sj,kp  +sk), &
                                       sdw(ip  +si,jp+1+sj,kp+1+sk), &
                                       sdw(ip+1+si,jp  +sj,kp  +sk), &
                                       sdw(ip+1+si,jp  +sj,kp+1+sk), &
                                       sdw(ip+1+si,jp+1+sj,kp  +sk), &
                                       sdw(ip+1+si,jp+1+sj,kp+1+sk), sdvel)
      else 
         call lpperror("Wrong direction in velocity interpolation!")
      end if ! dir
   end subroutine TrilinearIntrplFluidVel

   subroutine ComputeSubDerivativeVel(n,dli,dzci,dt,p,duconv,dvconv,dwconv,du,dv,dw)
      implicit none
      integer , intent(in), dimension(3)  :: n
      real(rp), intent(in), dimension(3)  :: dli
      real(rp), intent(in), dimension(0:) :: dzci
      real(rp), intent(in) :: dt
      real(rp), intent(in), dimension(0:,0:,0:) :: p,du,dv,dw
      real(rp), intent(in), dimension(:,:,:) :: duconv,dvconv,dwconv
      integer  :: i,j,k,ip,jp,kp
      attributes(managed) :: p,duconv,dvconv,dwconv,du,dv,dw,dzci

      ! Note: The current way to compute substantial derivatives of velocity
      ! only works when diffusion terms are solved EXPLICITLY

      if ( UnsteadyPartForce ) then !Correct velocity with convective terms substracted

         !$cuf kernel do(3) <<<*,*>>>
         do k=1,n(3)
           do j=1,n(2)
             do i=1,n(1)
               ip = i+1
               sdu(i,j,k) = (du(i,j,k) - duconv(i,j,k)) - dxi*(p(ip,j,k)-p(i,j,k))
             enddo
           enddo
         enddo


         !$cuf kernel do(3) <<<*,*>>>
         do k=1,n(3)
           do j=1,n(2)
              jp = j+1
             do i=1,n(1)
               sdv(i,j,k) = (dv(i,j,k) - dvconv(i,j,k)) - dyi*(p(i,jp,k)-p(i,j,k))
             enddo
           enddo
         enddo


         !$cuf kernel do(3) <<<*,*>>>
         do k=1,n(3)
            kp = k+1
           do j=1,n(2)
             do i=1,n(1)
               sdw(i,j,k) = (dw(i,j,k) - dwconv(i,j,k)) - dzci(k)*(p(i,j,kp)-p(i,j,k))
             enddo
           enddo
         enddo
      end if
    !@cuf istat=cudaDeviceSynchronize()
   end subroutine ComputeSubDerivativeVel

   subroutine SeedParticles(u,v,w,time,dt,z_c,z_f,n,ng)
	  implicit none
      integer , intent(in), dimension(3) :: n,ng
      real(rp), dimension(0:n(3)+1), intent(in) :: z_c,z_f
      real(rp), dimension(0:n(1)+1,0:n(2)+1,0:n(3)+1), intent(in) :: u,v,w
      real(rp), intent(in) :: time,dt
      logical :: SeedParticlesNow
      integer :: ipart_seed,itime,itime_part_seed
      real(rp) :: xp,yp,zp,up,vp,wp,dp,volp
      real(rp) :: rand
      integer :: ip,jp,kp
!      integer :: c1,c2,c3,rankpart
      integer :: num_part_rank
      real(rp) :: dist,dp1
      logical :: PartSeedTooClose
      integer :: ipart,ipart1,iter,niters
      real(rp) :: xb1,xb2,yb1,yb2,zb1,zb2
      integer :: num_part_seeded(0:(dims(1)*dims(2))-1),rankseed,c1,c2,c3
      real(rp) :: dummyreal,uf,vf,wf
      real(rp) :: xh1,xh2,yh1,yh2,temp1,temp2,temp3,temp4,temp5,temp6, &
                  xl,xr,yl,yr,zl,zr
      integer :: si,sj,sk
      integer :: MPI_errorcode
      integer, parameter :: largeint = 1234567890
      attributes(managed) :: u,v,w,z_c,z_f

      istat = cudaMemPrefetchAsync( parts%ic,  size(parts%ic), cudaCpuDeviceId, 0 )
      istat = cudaMemPrefetchAsync( parts%jc,  size(parts%jc), cudaCpuDeviceId, 0 )
      istat = cudaMemPrefetchAsync( parts%kc,  size(parts%kc), cudaCpuDeviceId, 0 )
      istat = cudaMemPrefetchAsync( parts%xc,  size(parts%xc), cudaCpuDeviceId, 0 )
      istat = cudaMemPrefetchAsync( parts%yc,  size(parts%yc), cudaCpuDeviceId, 0 )
      istat = cudaMemPrefetchAsync( parts%zc,  size(parts%zc), cudaCpuDeviceId, 0 )
      istat = cudaMemPrefetchAsync( parts%uc,  size(parts%uc), cudaCpuDeviceId, 0 )
      istat = cudaMemPrefetchAsync( parts%vc,  size(parts%vc), cudaCpuDeviceId, 0 )
      istat = cudaMemPrefetchAsync( parts%wc,  size(parts%wc), cudaCpuDeviceId, 0 )
      istat = cudaMemPrefetchAsync( parts%vol, size(parts%vol), cudaCpuDeviceId, 0 )
      istat = cudaMemPrefetchAsync( num_part,  size(num_part), cudaCpuDeviceId, 0 )

      MPI_errorcode = 0
      up = 0._rp; vp = 0._rp; wp = 0._rp
      uf = 0._rp; vf = 0._rp; wf = 0._rp

      SeedParticlesNow = .false.
      if ( SeedParticlesFlag == SeedParticlesNone ) then 
         return 
      else if ( SeedParticlesFlag == SeedParticlesOneTime ) then
         itime           = nint(time/dt)
         itime_part_seed = nint(time_part_seed/dt)
         if ( itime == itime_part_seed ) then 
            SeedParticlesNow = .true.
         end if ! t
      else if ( SeedParticlesFlag == SeedParticlesRegular ) then
         itime           = nint(time/dt)
         itime_part_seed = nint(time_part_seed/dt)
         if ( mod(itime,itime_part_seed) == 0 ) then 
            SeedParticlesNow = .true.
         end if ! itime_part_seed
      else 
         call lpperror("Unknown seeding particle flag!")
      end if ! SeedParticlesFlag

      num_part_seeded = 0 
      if ( SeedParticlesNow ) then
         call random_seed()
         if ( myid == 0 ) write(*,*) num_part_seed,'particles to seed at t =',time
         ipart_seed = 0
         do
            ipart_seed = ipart_seed + 1
            call generate_random_particle(dmin_part_seed,dmax_part_seed,     & 
                 xmin_part_seed,xmax_part_seed,ymin_part_seed,ymax_part_seed,&
                 zmin_part_seed,zmax_part_seed,umin_part_seed,umax_part_seed,&
                 vmin_part_seed,vmax_part_seed,wmin_part_seed,wmax_part_seed,& 
                 dp,volp,xp,yp,zp,up,vp,wp)

            if ( xp > (0+coord(1)*n(1))*dx .and. xp <= (n(1)+coord(1)*n(1))*dx .and. & 
                 yp > (0+coord(2)*n(2))*dy .and. yp <= (n(2)+coord(2)*n(2))*dy .and. & 
                 zp > z_f(0) .and. zp <= z_f(n(3)) ) then   
               num_part(myid) = num_part(myid) + 1
               parts(num_part(myid),myid)%xc = xp
               parts(num_part(myid),myid)%yc = yp
               parts(num_part(myid),myid)%zc = zp
               parts(num_part(myid),myid)%vol = volp
               call FindPartLocCell(xp,yp,zp,ip,jp,kp,z_f,n,ng)
               parts(num_part(myid),myid)%ic = ip
               parts(num_part(myid),myid)%jc = jp
               parts(num_part(myid),myid)%kc = kp

               if ( SeedWithFluidVel ) then
                   call GetFluidProp(ip,jp,kp, &
                                     xp,yp,zp, &
                                     uf,vf,wf,dummyreal,dummyreal,dummyreal,   & 
                                     volp,&
                                     TwoWayCouplingFlag, &
                                     u,v,w,z_c,z_f,n)
                  !! GetFluidProp
                  ! Check particle location
                  ! if( xp > (n(1)+coord(1)*n(1))*dx .or. xp < ((0+coord(1)*n(1))*dx-dx) ) then 
                     ! write(*,*) ip,jp,kp,xp,yp,zp
                     ! MPI_errorcode=1
                     ! print*, "LPP ERROR *** ","Failed to get fluid properties as particle out x range", " *** STOP at processor: ", mydev
                     ! exit
                  ! else if( yp > (n(2)+coord(2)*n(2))*dy .or. yp < ((0+coord(2)*n(2))*dy-dy) ) then
                     ! write(*,*) ip,jp,kp,xp,yp,zp
                     ! MPI_errorcode=1
                     ! print*, "LPP ERROR *** ","Failed to get fluid properties as particle out y range", " *** STOP at processor: ", mydev
                     ! exit
                  ! else if( zp > z_f(n(3)+1) .or. zp < z_f(0) ) then 
                     ! write(*,*) ip,jp,kp,xp,yp,zp
                     ! MPI_errorcode=1
                     ! print*, "LPP ERROR *** ","Failed to get fluid properties as particle out z range", " *** STOP at processor: ", mydev
                     ! exit
                  ! end if ! ip,jp,kp
                  
                  ! if (TwoWayCouplingFlag == TwoWayCouplingIgnore .or. TwoWayCouplingFlag == TwoWayCouplingFilterForce ) then        
                   
                   
                      ! si = -1;sj = -1;sk = -1 ! Trilinear interpolation for u velocity
                      ! if ( yp > (jp+coord(2)*n(2))*dy ) sj =  0 
                      ! if ( zp > z_c(kp) ) sk =  0

                       ! xh1 = ((ip+si)+coord(1)*n(1))*dx; xh2 = ((ip+1+si)+coord(1)*n(1))*dx
                       ! yh1 = ((jp+sj)+coord(2)*n(2))*dy; yh2 = ((jp+1+sj)+coord(2)*n(2))*dy
                       
                       ! xl = (xp-xh1)/(xh2-xh1)
                       ! xr = 1._rp - xh2
                       ! temp1  = u(ip  +si,jp  +sj,kp  +sk)*xr + u(ip+1+si,jp  +sj,kp  +sk)*xl
                       
                       ! xl = (xp-xh1)/(xh2-xh1)
                       ! xr = 1._rp - xh2
                       ! temp2  = u(ip  +si,jp+1+sj,kp  +sk)*xr + u(ip+1+si,jp+1+sj,kp  +sk)*xl
                       
                       ! xl = (xp-xh1)/(xh2-xh1)
                       ! xr = 1._rp - xh2
                       ! temp3  = u(ip  +si,jp  +sj,kp+1+sk)*xr + u(ip+1+si,jp  +sj,kp+1+sk)*xl
                       
                       ! xl = (xp-xh1)/(xh2-xh1)
                       ! xr = 1._rp - xh2
                       ! temp4  = u(ip  +si,jp+1+sj,kp+1+sk)*xr + u(ip+1+si,jp+1+sj,kp+1+sk)*xl
                       
                       ! yl = (yp-yh1)/(yh2-yh1)
                       ! yr = 1._rp - yh2
                       ! temp5  = temp1*yr + temp2*yl
                       
                       ! yl = (yp-yh1)/(yh2-yh1)
                       ! yr = 1._rp - yh2
                       ! temp6  = temp3*yr + temp4*yl
                       
                       ! zl = (zp-z_f(kp  +sk))/(z_f(kp+1+sk)-z_f(kp  +sk))
                       ! zr = 1._rp - z_f(kp+1+sk)
                       ! uf  = temp5*zr + temp6*zl
                  
                  
                      ! si = -1;sj = -1;sk = -1 ! Trilinear interpolation for v velocity
                      ! if ( xp > (ip+coord(1)*n(1))*dx ) si =  0 
                      ! if ( zp > z_c(kp) ) sk =  0

                       ! xh1 = ((ip+si)+coord(1)*n(1))*dx; xh2 = ((ip+1+si)+coord(1)*n(1))*dx
                       ! yh1 = ((jp+sj)+coord(2)*n(2))*dy; yh2 = ((jp+1+sj)+coord(2)*n(2))*dy
                  
                       ! xl = (xp-xh1)/(xh2-xh1)
                       ! xr = 1._rp - xh2
                       ! temp1  = v(ip  +si,jp  +sj,kp  +sk)*xr + v(ip+1+si,jp  +sj,kp  +sk)*xl
                       
                       ! xl = (xp-xh1)/(xh2-xh1)
                       ! xr = 1._rp - xh2
                       ! temp2  = v(ip  +si,jp+1+sj,kp  +sk)*xr + v(ip+1+si,jp+1+sj,kp  +sk)*xl
                       
                       ! xl = (xp-xh1)/(xh2-xh1)
                       ! xr = 1._rp - xh2
                       ! temp3  = v(ip  +si,jp  +sj,kp+1+sk)*xr + v(ip+1+si,jp  +sj,kp+1+sk)*xl
                       
                       ! xl = (xp-xh1)/(xh2-xh1)
                       ! xr = 1._rp - xh2
                       ! temp4  = v(ip  +si,jp+1+sj,kp+1+sk)*xr + v(ip+1+si,jp+1+sj,kp+1+sk)*xl
                       
                       ! yl = (yp-yh1)/(yh2-yh1)
                       ! yr = 1._rp - yh2
                       ! temp5  = temp1*yr + temp2*yl
                       
                       ! yl = (yp-yh1)/(yh2-yh1)
                       ! yr = 1._rp - yh2
                       ! temp6  = temp3*yr + temp4*yl
                       
                       ! zl = (zp-z_f(kp  +sk))/(z_f(kp+1+sk)-z_f(kp  +sk))
                       ! zr = 1._rp - z_f(kp+1+sk)
                       ! vf  = temp5*zr + temp6*zl
                  
                   
                      ! si = -1;sj = -1;sk = -1 ! Trilinear interpolation for w
                      ! if ( xp > (ip+coord(1)*n(1))*dx ) si =  0 
                      ! if ( yp > (jp+coord(2)*n(2))*dy ) sj =  0

                       ! xh1 = ((ip+si)+coord(1)*n(1))*dx; xh2 = ((ip+1+si)+coord(1)*n(1))*dx
                       ! yh1 = ((jp+sj)+coord(2)*n(2))*dy; yh2 = ((jp+1+sj)+coord(2)*n(2))*dy
                  
                       ! xl = (xp-xh1)/(xh2-xh1)
                       ! xr = 1._rp - xh2
                       ! temp1  = w(ip  +si,jp  +sj,kp  +sk)*xr + w(ip+1+si,jp  +sj,kp  +sk)*xl
                       
                       ! xl = (xp-xh1)/(xh2-xh1)
                       ! xr = 1._rp - xh2
                       ! temp2  = w(ip  +si,jp+1+sj,kp  +sk)*xr + w(ip+1+si,jp+1+sj,kp  +sk)*xl
                       
                       ! xl = (xp-xh1)/(xh2-xh1)
                       ! xr = 1._rp - xh2
                       ! temp3  = w(ip  +si,jp  +sj,kp+1+sk)*xr + w(ip+1+si,jp  +sj,kp+1+sk)*xl
                       
                       ! xl = (xp-xh1)/(xh2-xh1)
                       ! xr = 1._rp - xh2
                       ! temp4  = w(ip  +si,jp+1+sj,kp+1+sk)*xr + w(ip+1+si,jp+1+sj,kp+1+sk)*xl
                       
                       ! yl = (yp-yh1)/(yh2-yh1)
                       ! yr = 1._rp - yh2
                       ! temp5  = temp1*yr + temp2*yl
                       
                       ! yl = (yp-yh1)/(yh2-yh1)
                       ! yr = 1._rp - yh2
                       ! temp6  = temp3*yr + temp4*yl
                       
                       ! zl = (zp-z_f(kp  +sk))/(z_f(kp+1+sk)-z_f(kp  +sk))
                       ! zr = 1._rp - z_f(kp+1+sk)
                       ! wf  = temp5*zr + temp6*zl
                  ! else if ( TwoWayCouplingFlag == TwoWayCouplingFilterVel ) then
                   ! MPI_errorcode=1
                   ! print*, "LPP ERROR *** ","Two way coupling with filtered fluid velocity is not ready yet!", " *** STOP at processor: ", myid
                   ! exit
                  ! end if
                  !! GetFluidProp
                  up = uf; vp = vf; wp = wf
               end if ! SeedWithFluidVel
               parts(num_part(myid),myid)%uc = up
               parts(num_part(myid),myid)%vc = vp
               parts(num_part(myid),myid)%wc = wp
            end if ! xp,yp,zp 
            if ( ipart_seed >= num_part_seed ) exit
         end do ! ipart_seed
         if(MPI_errorcode.eq.1) then
            call MPI_ABORT(MPI_COMM_WORLD, MPI_errorcode, ierr)
            call MPI_finalize(ierr)
            error stop
         endif
         ! Adjust particle positions to avoid locating at block bdry 
         if ( AvoidSeedAtBlockBdry .and. num_part(myid) > 0 ) then
            do ipart = 1,num_part(myid)
               dp = (parts(ipart,myid)%vol*6._rp/PI)**OneThird
               xp = parts(ipart,myid)%xc
               yp = parts(ipart,myid)%yc
               zp = parts(ipart,myid)%zc
               if ( xp - 0.65_rp*dp < (0+coord(1)*n(1))*dx )    xp = (0+coord(1)*n(1))*dx    + 0.65_rp*dp
               if ( xp + 0.65_rp*dp > (n(1)+coord(1)*n(1))*dx ) xp = (n(1)+coord(1)*n(1))*dx - 0.65_rp*dp
               if ( yp - 0.65_rp*dp < (0+coord(2)*n(2))*dy )    yp = (0+coord(2)*n(2))*dy    + 0.65_rp*dp
               if ( yp + 0.65_rp*dp > (n(2)+coord(2)*n(2))*dy ) yp = (n(2)+coord(2)*n(2))*dy - 0.65_rp*dp
               if ( zp - 0.65_rp*dp < z_f(0) )                  zp = z_f(0)                   + 0.65_rp*dp
               if ( zp + 0.65_rp*dp > z_f(n(3))   )             zp = z_f(n(3))                - 0.65_rp*dp
               parts(ipart,myid)%xc = xp 
               parts(ipart,myid)%yc = yp
               parts(ipart,myid)%zc = zp
            end do !ipart
         end if !AvoidSeedAtBlockBdry
           
         ! Adjust particle positions to avoid particle too close to each other 
         if ( AvoidSeedPartTooClose .and. num_part(myid) > 1 ) then
            niters = 5000

            do iter = 1,niters
               do ipart = 1,num_part(myid)
                  ! Check if a particle is too close with others
                  PartSeedTooClose = .false.
                  dp = (parts(ipart,myid)%vol*6._rp/PI)**OneThird
                  do ipart1 = 1,num_part(myid)
                     if ( ipart1 /= ipart ) then 
                        dist = (parts(ipart ,myid)%xc - & 
                                parts(ipart1,myid)%xc)**2 & 
                             + (parts(ipart ,myid)%yc - &
                                parts(ipart1,myid)%yc)**2 &
                             + (parts(ipart ,myid)%zc - &
                                parts(ipart1,myid)%zc)**2
                        dp1  = (parts(ipart1,myid)%vol*6._rp/PI)**OneThird
                        if ( sqrt(dist) < (dp+dp1)*1.3_rp ) then 
                           PartSeedTooClose = .true.
                           exit
                        end if ! sqrt(dist)
                     end if ! ipart1 
                  end do ! ipart1
                  ! relocate the particle within the block
                  if ( PartSeedTooClose ) then
                     xb1 = max((0+coord(1)*n(1))*dx   +0.65_rp*dp,xmin_part_seed)
                     xb2 = min((n(1)+coord(1)*n(1))*dx-0.65_rp*dp,xmax_part_seed)
                     call random_number(rand) 
                     parts(ipart,myid)%xc = xb1 + rand*(xb2-xb1)
                                 
                     yb1 = max((0+coord(2)*n(2))*dy   +0.65_rp*dp,ymin_part_seed) 
                     yb2 = min((n(2)+coord(2)*n(2))*dy-0.65_rp*dp,ymax_part_seed)
                     call random_number(rand) 
                     parts(ipart,myid)%yc = yb1 + rand*(yb2-yb1)

                     zb1 = max(z_f(0)   +0.65_rp*dp,zmin_part_seed) 
                     zb2 = min(z_f(n(3))-0.65_rp*dp,zmax_part_seed) 
                     call random_number(rand) 
                     parts(ipart,myid)%zc = zb1 + rand*(zb2-zb1)
                     exit
                  end if ! PartSeedTooClose
               end do ! ipart
               if ( .NOT.PartSeedTooClose ) exit
            end do ! iter 
            write(*,*) 'iteration on adjusting particle locations',myid,iter,num_part(myid)
         end if !AvoidSeedPartTooClose

         ! recompute particle velocity if seed with fluid vel and part positions
         ! have been adjusted
         if ( SeedWithFluidVel .and. (AvoidSeedAtBlockBdry .or. AvoidSeedPartTooClose) ) then
            do ipart = 1, num_part(myid)
               xp = parts(ipart,myid)%xc
               yp = parts(ipart,myid)%yc
               zp = parts(ipart,myid)%zc
               call FindPartLocCell(xp,yp,zp,ip,jp,kp,z_f,n,ng)
               parts(ipart,myid)%ic = ip
               parts(ipart,myid)%jc = jp
               parts(ipart,myid)%kc = kp
               !! Get Fluid Property         
               ! Check particle location
               if( xp > (n(1)+coord(1)*n(1))*dx .or. xp < ((0+coord(1)*n(1))*dx-dx) ) then 
                  write(*,*) ip,jp,kp,xp,yp,zp
                  MPI_errorcode=1
                  print*, "LPP ERROR *** ","Failed to get fluid properties as particle out x range", " *** STOP at processor: ", myid
                  exit
               else if( yp > (n(2)+coord(2)*n(2))*dy .or. yp < ((0+coord(2)*n(2))*dy-dy) ) then
                  write(*,*) ip,jp,kp,xp,yp,zp
                  MPI_errorcode=1
                  print*, "LPP ERROR *** ","Failed to get fluid properties as particle out y range", " *** STOP at processor: ", myid
                  exit
               else if( zp > z_f(n(3)+1) .or. zp < z_f(0) ) then 
                  write(*,*) ip,jp,kp,xp,yp,zp
                  MPI_errorcode=1
                  print*, "LPP ERROR *** ","Failed to get fluid properties as particle out z range", " *** STOP at processor: ", myid
                  exit
               end if ! ip,jp,kp
               
               if (TwoWayCouplingFlag == TwoWayCouplingIgnore .or. TwoWayCouplingFlag == TwoWayCouplingFilterForce ) then        
                
                ! Trilinear interpolation for u velocity
                   si = -1;sj = -1;sk = -1 
                   if ( yp > (jp+coord(2)*n(2))*dy ) sj =  0 
                   if ( zp > z_c(kp) ) sk =  0

                    xh1 = ((ip+si)+coord(1)*n(1))*dx; xh2 = ((ip+1+si)+coord(1)*n(1))*dx
                    yh1 = ((jp+sj)+coord(2)*n(2))*dy; yh2 = ((jp+1+sj)+coord(2)*n(2))*dy
                    
                    xl = (xp-xh1)/(xh2-xh1)
                    xr = 1._rp - xh2
                    temp1  = u(ip  +si,jp  +sj,kp  +sk)*xr + u(ip+1+si,jp  +sj,kp  +sk)*xl
                    !
                    xl = (xp-xh1)/(xh2-xh1)
                    xr = 1._rp - xh2
                    temp2  = u(ip  +si,jp+1+sj,kp  +sk)*xr + u(ip+1+si,jp+1+sj,kp  +sk)*xl
                    !
                    xl = (xp-xh1)/(xh2-xh1)
                    xr = 1._rp - xh2
                    temp3  = u(ip  +si,jp  +sj,kp+1+sk)*xr + u(ip+1+si,jp  +sj,kp+1+sk)*xl
                    !
                    xl = (xp-xh1)/(xh2-xh1)
                    xr = 1._rp - xh2
                    temp4  = u(ip  +si,jp+1+sj,kp+1+sk)*xr + u(ip+1+si,jp+1+sj,kp+1+sk)*xl
                    !
                    yl = (yp-yh1)/(yh2-yh1)
                    yr = 1._rp - yh2
                    temp5  = temp1*yr + temp2*yl
                    !
                    yl = (yp-yh1)/(yh2-yh1)
                    yr = 1._rp - yh2
                    temp6  = temp3*yr + temp4*yl
                    !
                    zl = (zp-z_f(kp  +sk))/(z_f(kp+1+sk)-z_f(kp  +sk))
                    zr = 1._rp - z_f(kp+1+sk)
                    uf  = temp5*zr + temp6*zl
               
               ! Trilinear interpolation for v velocity
                   si = -1;sj = -1;sk = -1 
                   if ( xp > (ip+coord(1)*n(1))*dx ) si =  0 
                   if ( zp > z_c(kp) ) sk =  0

                    xh1 = ((ip+si)+coord(1)*n(1))*dx; xh2 = ((ip+1+si)+coord(1)*n(1))*dx
                    yh1 = ((jp+sj)+coord(2)*n(2))*dy; yh2 = ((jp+1+sj)+coord(2)*n(2))*dy
               
                    xl = (xp-xh1)/(xh2-xh1)
                    xr = 1._rp - xh2
                    temp1  = v(ip  +si,jp  +sj,kp  +sk)*xr + v(ip+1+si,jp  +sj,kp  +sk)*xl
                    !
                    xl = (xp-xh1)/(xh2-xh1)
                    xr = 1._rp - xh2
                    temp2  = v(ip  +si,jp+1+sj,kp  +sk)*xr + v(ip+1+si,jp+1+sj,kp  +sk)*xl
                    !
                    xl = (xp-xh1)/(xh2-xh1)
                    xr = 1._rp - xh2
                    temp3  = v(ip  +si,jp  +sj,kp+1+sk)*xr + v(ip+1+si,jp  +sj,kp+1+sk)*xl
                    !
                    xl = (xp-xh1)/(xh2-xh1)
                    xr = 1._rp - xh2
                    temp4  = v(ip  +si,jp+1+sj,kp+1+sk)*xr + v(ip+1+si,jp+1+sj,kp+1+sk)*xl
                    !
                    yl = (yp-yh1)/(yh2-yh1)
                    yr = 1._rp - yh2
                    temp5  = temp1*yr + temp2*yl
                    !
                    yl = (yp-yh1)/(yh2-yh1)
                    yr = 1._rp - yh2
                    temp6  = temp3*yr + temp4*yl
                    !
                    zl = (zp-z_f(kp  +sk))/(z_f(kp+1+sk)-z_f(kp  +sk))
                    zr = 1._rp - z_f(kp+1+sk)
                    vf  = temp5*zr + temp6*zl
               
                ! Trilinear interpolation for w
                   si = -1;sj = -1;sk = -1 
                   if ( xp > (ip+coord(1)*n(1))*dx ) si =  0 
                   if ( yp > (jp+coord(2)*n(2))*dy ) sj =  0

                    xh1 = ((ip+si)+coord(1)*n(1))*dx; xh2 = ((ip+1+si)+coord(1)*n(1))*dx
                    yh1 = ((jp+sj)+coord(2)*n(2))*dy; yh2 = ((jp+1+sj)+coord(2)*n(2))*dy
               
                    xl = (xp-xh1)/(xh2-xh1)
                    xr = 1._rp - xh2
                    temp1  = w(ip  +si,jp  +sj,kp  +sk)*xr + w(ip+1+si,jp  +sj,kp  +sk)*xl
                    !
                    xl = (xp-xh1)/(xh2-xh1)
                    xr = 1._rp - xh2
                    temp2  = w(ip  +si,jp+1+sj,kp  +sk)*xr + w(ip+1+si,jp+1+sj,kp  +sk)*xl
                    !
                    xl = (xp-xh1)/(xh2-xh1)
                    xr = 1._rp - xh2
                    temp3  = w(ip  +si,jp  +sj,kp+1+sk)*xr + w(ip+1+si,jp  +sj,kp+1+sk)*xl
                    !
                    xl = (xp-xh1)/(xh2-xh1)
                    xr = 1._rp - xh2
                    temp4  = w(ip  +si,jp+1+sj,kp+1+sk)*xr + w(ip+1+si,jp+1+sj,kp+1+sk)*xl
                    !
                    yl = (yp-yh1)/(yh2-yh1)
                    yr = 1._rp - yh2
                    temp5  = temp1*yr + temp2*yl
                    !
                    yl = (yp-yh1)/(yh2-yh1)
                    yr = 1._rp - yh2
                    temp6  = temp3*yr + temp4*yl
                    !
                    zl = (zp-z_f(kp  +sk))/(z_f(kp+1+sk)-z_f(kp  +sk))
                    zr = 1._rp - z_f(kp+1+sk)
                    wf  = temp5*zr + temp6*zl
               else if ( TwoWayCouplingFlag == TwoWayCouplingFilterVel ) then
                MPI_errorcode=1
                print*, "LPP ERROR *** ","Two way coupling with filtered fluid velocity is not ready yet!", " *** STOP at processor: ", myid
                exit
               end if
               !! GetFluidProp
               up = uf; vp = vf; wp = wf
               parts(ipart,myid)%uc = up
               parts(ipart,myid)%vc = vp
               parts(ipart,myid)%wc = wp
            end do ! ipart
            if(MPI_errorcode.eq.1) then
               call MPI_ABORT(MPI_COMM_WORLD, MPI_errorcode, ierr)
               call MPI_finalize(ierr)
               error stop
            endif
         end if ! SeedWithFluidVel

         ! Allgather updated number of particles after seeding
         num_part_rank = num_part(myid)
         call MPI_ALLGATHER(num_part_rank, 1, MPI_INTEGER, &
                            num_part(:)  , 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
         if ( myid == 0 ) write(*,*) sum(num_part), "particles seeded!"
      end if ! SeedParticlesNow

  istat = cudaMemPrefetchAsync( parts%ic, size(parts%ic), mydev, 0 )
  istat = cudaMemPrefetchAsync( parts%jc, size(parts%jc), mydev, 0 )
  istat = cudaMemPrefetchAsync( parts%kc, size(parts%kc), mydev, 0 )
  istat = cudaMemPrefetchAsync( parts%xc, size(parts%xc), mydev, 0 )
  istat = cudaMemPrefetchAsync( parts%yc, size(parts%yc), mydev, 0 )
  istat = cudaMemPrefetchAsync( parts%zc, size(parts%zc), mydev, 0 )
  istat = cudaMemPrefetchAsync( parts%uc, size(parts%uc), mydev, 0 )
  istat = cudaMemPrefetchAsync( parts%vc, size(parts%vc), mydev, 0 )
  istat = cudaMemPrefetchAsync( parts%wc, size(parts%wc), mydev, 0 )
  istat = cudaMemPrefetchAsync( parts%vol, size(parts%vol), mydev, 0 )
  istat = cudaMemPrefetchAsync( num_part,  size(num_part), mydev, 0 )

   end subroutine SeedParticles

   subroutine generate_random_particle(dmin_part_seed,dmax_part_seed,& 
                                       xmin_part_seed,xmax_part_seed,ymin_part_seed,ymax_part_seed,&
                                       zmin_part_seed,zmax_part_seed,umin_part_seed,umax_part_seed,&
                                       vmin_part_seed,vmax_part_seed,wmin_part_seed,wmax_part_seed,& 
                                       dp,volp,xp,yp,zp,up,vp,wp)
      implicit none
      real(rp),intent(in) :: dmin_part_seed,dmax_part_seed,   &
                             xmin_part_seed,xmax_part_seed,ymin_part_seed,ymax_part_seed,&
                             zmin_part_seed,zmax_part_seed,umin_part_seed,umax_part_seed,&
                             vmin_part_seed,vmax_part_seed,wmin_part_seed,wmax_part_seed
      real(rp),intent(out) :: dp,volp,xp,yp,zp,up,vp,wp
      real(rp) :: rand

      call random_number(rand) 
      dp = dmin_part_seed + rand*(dmax_part_seed-dmin_part_seed)
      volp = PI*dp*dp*dp/6._rp
      call random_number(rand) 
      xp = xmin_part_seed + rand*(xmax_part_seed-xmin_part_seed)
      call random_number(rand) 
      yp = ymin_part_seed + rand*(ymax_part_seed-ymin_part_seed)
      call random_number(rand) 
      zp = zmin_part_seed + rand*(zmax_part_seed-zmin_part_seed)
      call random_number(rand) 
      up = umin_part_seed + rand*(umax_part_seed-umin_part_seed)
      call random_number(rand) 
      vp = vmin_part_seed + rand*(vmax_part_seed-vmin_part_seed)
      call random_number(rand) 
      wp = wmin_part_seed + rand*(wmax_part_seed-wmin_part_seed)
   end subroutine generate_random_particle

   subroutine lpperror(message) 
	  implicit none
      integer :: MPI_errorcode
      character(*) :: message
      MPI_errorcode=1
      print*, "LPP ERROR *** ",message, " *** STOP at processor ranked: ", mydev
      call MPI_ABORT(MPI_COMM_WORLD, MPI_errorcode, ierr)
      call MPI_finalize(ierr)
      error stop
   end subroutine lpperror

! ==================================================================================================
! Output routines for Lagrangian particles
! ==================================================================================================
   subroutine outLPP(nf,time,i1,i2,j1,j2,k1,k2)
      implicit none

      real(rp), intent(in) :: time
      integer,  intent(in) :: nf,i1,i2,j1,j2,k1,k2
      integer,  parameter  :: Plot3D = 1
      integer,  parameter  :: Decomp = 2

      if ( outputlpp_format == Plot3D ) then 
         call output_LPP_Plot3D(nf,time,i1,i2,j1,j2,k1,k2)
      else if ( outputlpp_format == Decomp ) then
         call output_LPP_2DECOMP(nf,time)
      end if

   end subroutine outLPP
!-------------------------------------------------------------------------------------------------
   subroutine append_LPP_visit_file(rootname)
      implicit none

      character(*) :: rootname
      character(len=30) :: fldnum
      character(len=100) :: filename
      integer :: prank, lpp_opened=0

      if(mydev.ne.0) call lpperror('mydev.ne.0 in append_LPP')
      if(lpp_opened==0) then
         OPEN(UNIT=88,FILE=TRIM(datadir)//'Particle/'//'lpp.visit')
         write(88,10) (dims(1)*dims(2))
10       format('!NBLOCKS ',I4)
         lpp_opened=1
      else
         OPEN(UNIT=88,FILE=TRIM(datadir)//'Particle/'//'lpp.visit',position='append')
      endif
      do prank=0,(dims(1)*dims(2))-1
          write(fldnum,'(i4.4)') prank
          filename = TRIM(rootname)//TRIM(fldnum)//'.3D'
          write(88,11) TRIM(filename)
11        format(A)
      enddo
      close(88)
   end subroutine  append_LPP_visit_file
!-------------------------------------------------------------------------------------------------
   subroutine output_LPP_Plot3D(istep,time,i1,i2,j1,j2,k1,k2)
      implicit none
      real(rp),intent(in)   :: time
      integer,intent(in):: istep,i1,i2,j1,j2,k1,k2
      character(len=30) :: rootname, fldnum1, fldnum2
      character(len=100) :: filename
      integer :: ipart
      
      write(fldnum1,'(i7.7)') istep
      write(fldnum2,'(i4.4)') mydev

      ! output lpp data in plot3d format 
      rootname = 'LPP'//TRIM(fldnum2)//'_'
      if(mydev==0) call append_LPP_visit_file(rootname)
      filename = TRIM(rootname)//TRIM(fldnum1)//'.3D'
      OPEN(UNIT=8,FILE=TRIM(datadir)//'Particle/'//TRIM(filename))
!      write(8,10)
      write(8,11)
!10    format('# plot3D data file')
11    format('x y z vol')

      if ( num_part(mydev) > 0 ) then 
         do ipart = 1,num_part(mydev) 
            write(8,320) &
            parts(ipart,mydev)%xc,& 
            parts(ipart,mydev)%yc, & 
            parts(ipart,mydev)%zc, &  
            parts(ipart,mydev)%vol
         enddo
      end if ! num_part(mydev)

      ! add virtual particles at the 8 corners of the block, so that
      ! molecular plot in VisIt would work in the same scale of vof
      ! write(8,320) x(i1),y(j1),z(k1),volsmall 
      ! write(8,320) x(i1),y(j1),z(k2),volsmall 
      ! write(8,320) x(i1),y(j2),z(k1),volsmall 
      ! write(8,320) x(i1),y(j2),z(k2),volsmall 
      ! write(8,320) x(i2),y(j1),z(k1),volsmall 
      ! write(8,320) x(i2),y(j1),z(k2),volsmall 
      ! write(8,320) x(i2),y(j2),z(k1),volsmall 
      ! write(8,320) x(i2),y(j2),z(k2),volsmall 
320   format(e14.5,e14.5,e14.5,e14.5)
      close(8)

      ! write complete lpp data for other post-processing treatment
      ! rootname = 'LPPDATA'//TRIM(fldnum2)//'_'
      ! filename = TRIM(rootname)//TRIM(fldnum1)//'.dat'
      ! OPEN(UNIT=8,FILE=TRIM(datadir)//'Particle/'//TRIM(filename))
      ! write(8,*) time,num_part(mydev) 
      ! if ( num_part(mydev) > 0 ) then 
         ! do ipart = 1,num_part(mydev) 
            ! write(8,'(11(E15.8,1X))') & 
            ! parts(ipart,mydev)%xc,& 
            ! parts(ipart,mydev)%yc, & 
            ! parts(ipart,mydev)%zc, &  
            ! parts(ipart,mydev)%uc, &  
            ! parts(ipart,mydev)%vc, &  
            ! parts(ipart,mydev)%wc, &  
            ! parts(ipart,mydev)%uf, &  
            ! parts(ipart,mydev)%vf, &  
            ! parts(ipart,mydev)%wf, &  
            ! parts(ipart,mydev)%vol, &
            ! parts(ipart,mydev)%sur
         ! enddo
      ! end if ! num_part(mydev)
      ! close(8)

   end subroutine output_LPP_Plot3D
!-------------------------------------------------------------------------------------------------
   subroutine output_LPP_2DECOMP(istep,time)

      implicit none
      real(rp), intent(in) :: time
      real(rp), dimension(3) :: fldinfo
      integer,  intent(in) :: istep
      real(rp), dimension(num_part(mydev)):: axp,ayp,azp,aup,avp,awp
      integer :: fh, ipart
      integer(kind=MPI_OFFSET_KIND) :: filesize,disp,good
      character(len=7) :: fldnum
      character(len=100) :: filename
      
      ! write(fldnum,'(i7.7)') istep
      
      ! if ( num_part(mydev) > 0 ) then 
         ! do ipart = 1,num_part(mydev)
          ! axp(ipart) = parts(ipart,mydev)%xc
          ! ayp(ipart) = parts(ipart,mydev)%yc
          ! azp(ipart) = parts(ipart,mydev)%zc
          ! aup(ipart) = parts(ipart,mydev)%uc
          ! avp(ipart) = parts(ipart,mydev)%vc
          ! awp(ipart) = parts(ipart,mydev)%wc
         ! enddo
      ! endif

      ! filename = trim(datadir)//'Particle/'//'pos_'//trim(fldnum)//'.bin'
      ! call MPI_FILE_OPEN(MPI_COMM_WORLD, filename, &
                         ! MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL,fh, ierr)
      ! filesize = 0_MPI_OFFSET_KIND
      ! call MPI_FILE_SET_SIZE(fh,filesize,ierr)  ! guarantee overwriting
      ! disp = 0_MPI_OFFSET_KIND
      ! call decomp_2d_write_var(fh,disp,3,axp)
      ! call decomp_2d_write_var(fh,disp,3,ayp)
      ! call decomp_2d_write_var(fh,disp,3,azp)
      ! call MPI_FILE_CLOSE(fh,ierr)
      
      ! filename = trim(datadir)//'Particle/'//'vel_'//trim(fldnum)//'.bin'
      ! call MPI_FILE_OPEN(MPI_COMM_WORLD, filename, &
                         ! MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL,fh, ierr)
      ! filesize = 0_MPI_OFFSET_KIND
      ! call MPI_FILE_SET_SIZE(fh,filesize,ierr)  ! guarantee overwriting
      ! disp = 0_MPI_OFFSET_KIND
      ! call decomp_2d_write_var(fh,disp,3,aup)
      ! call decomp_2d_write_var(fh,disp,3,avp)
      ! call decomp_2d_write_var(fh,disp,3,awp)
      ! call MPI_FILE_CLOSE(fh,ierr)
      
      ! filename = trim(datadir)//'Particle/'//'part_time'//fldnum//'.bin'
      ! call MPI_FILE_OPEN(MPI_COMM_WORLD, filename, &
                         ! MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL,fh, ierr)
      ! filesize = 0_MPI_OFFSET_KIND
      ! call MPI_FILE_SET_SIZE(fh,filesize,ierr)  ! guarantee overwriting
      ! disp = 0_MPI_OFFSET_KIND
      ! call decomp_2d_write_scalar(fh,disp,1,time)
      ! call MPI_FILE_CLOSE(fh,ierr)

   end subroutine output_LPP_2DECOMP
!-------------------------------------------------------------------------------------------------
   subroutine backup_LPP_write(time,itimestep)
      implicit none
      real(rp), intent(in) :: time, itimestep
      integer :: ipart,padding
      character(len=100) :: filename

      filename = trim(datadir)//'Particles/'//'backuplpp_'//int2text(mydev,padding)
      call system('touch '//trim(filename)//'; mv '//trim(filename)//' '//trim(filename)//'.old')
      OPEN(UNIT=7,FILE=trim(filename),status='REPLACE',ACTION='write')
      write(7,1100)time,itimestep,num_part(mydev)
      if ( num_part(mydev) > 0 ) then 
         do ipart=1,num_part(mydev)
            write(7,1200) parts(ipart,mydev)%xc, & 
                          parts(ipart,mydev)%yc, & 
                          parts(ipart,mydev)%zc, & 
                          parts(ipart,mydev)%uc, & 
                          parts(ipart,mydev)%vc, & 
                          parts(ipart,mydev)%wc, & 
                          parts(ipart,mydev)%duc, & 
                          parts(ipart,mydev)%dvc, & 
                          parts(ipart,mydev)%dwc, & 
                          parts(ipart,mydev)%vol, &  
                          parts(ipart,mydev)%sur, &  
                          parts(ipart,mydev)%ic, & 
                          parts(ipart,mydev)%jc, & 
                          parts(ipart,mydev)%kc 
         end do! ipart
      end if ! num_part(mydev)
      if(mydev==0)print*,'Backup LPP written at t=',time
      1100 FORMAT(es17.8e3,2I10)
      1200 FORMAT(11es17.8e3,3I5)
      CLOSE(7)

   end subroutine backup_LPP_write

!-------------------------------------------------------------------------------------------------
   subroutine backup_LPP_read
      implicit none
      real(rp) :: time, itimestep
      integer  :: ipart,padding
      character(len=100) :: filename

      filename = trim(datadir)//'Particles/'//'backuplpp_'//int2text(mydev,padding)
      OPEN(UNIT=7,FILE=trim(filename),status='old',action='read')
      read(7,*)time,itimestep,num_part(mydev)
      if ( num_part(mydev) < 0 ) &
         call lpperror("Error: backuplpp_read")
      if ( num_part(mydev) > 0 ) then 
         do ipart=1,num_part(mydev)
            read(7,*    ) parts(ipart,mydev)%xc, & 
                          parts(ipart,mydev)%yc, & 
                          parts(ipart,mydev)%zc, & 
                          parts(ipart,mydev)%uc, & 
                          parts(ipart,mydev)%vc, & 
                          parts(ipart,mydev)%wc, & 
                          parts(ipart,mydev)%duc, & 
                          parts(ipart,mydev)%dvc, & 
                          parts(ipart,mydev)%dwc, & 
                          parts(ipart,mydev)%vol, &  
                          parts(ipart,mydev)%sur, &  
                          parts(ipart,mydev)%ic, & 
                          parts(ipart,mydev)%jc, & 
                          parts(ipart,mydev)%kc 
         end do !ipart
      end if ! num_part(mydev)
      CLOSE(7)

   end subroutine backup_LPP_read
!=================================================================================================
! function text
!   Returns 'number' as a string with length of 'length'
!   called in:    function output
!-------------------------------------------------------------------------------------------------
   function int2text(number,length)
      implicit none
      integer :: number, length, i, MPI_errorcode
      character(len=length) :: int2text
      character, dimension(0:9) :: num = (/'0', '1', '2', '3', '4', '5', '6', '7', '8', '9'/)
      
      MPI_errorcode = 1
      if(number>=10**length) then
       print*, "int2text error: string is not large enough"
       call MPI_ABORT(MPI_COMM_WORLD, MPI_errorcode, ierr)
       call MPI_finalize(ierr)
       error stop
      endif
      do i=1,length
       int2text(length+1-i:length+1-i) = num(mod(number/(10**(i-1)),10))
      enddo
   end function
end module mod_Lag_part
#endif