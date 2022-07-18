module mod_sanity
  use iso_c_binding, only: C_PTR
  use mpi
  use decomp_2d
  use mod_bound     , only: boundp,bounduvw,updt_rhs_b
  use mod_chkdiv    , only: chkdiv
  use mod_common_mpi, only: myid,ierr
  use mod_correc    , only: correc
  use mod_debug     , only: chk_helmholtz
  use mod_fft       , only: fftend
  use mod_fillps    , only: fillps
  use mod_initflow  , only: add_noise
  use mod_initmpi   , only: initmpi
  use mod_initsolver, only: initsolver
  use mod_param     , only: small
  use mod_solver    , only: solver
  use mod_types
  implicit none
  private
  public test_sanity
  contains
  subroutine test_sanity(ng,dims,n,n_z,lo,hi,stop_type,cbcvel,cbcpre,bcvel,bcpre,is_forced, &
                         nb,is_bound,dli,dzci_g,dzfi_g,dzci,dzfi)
    !
    ! performs some a priori checks of the input files before the calculation starts
    !
    implicit none
    integer , intent(in), dimension(3) :: ng
    integer , intent(in), dimension(2) :: dims
    integer , intent(in), dimension(3) :: n,n_z,lo,hi
    logical , intent(in), dimension(3) :: stop_type
    character(len=1), intent(in), dimension(0:1,3,3) :: cbcvel
    character(len=1), intent(in), dimension(0:1,3)   :: cbcpre
    real(rp)        , intent(in), dimension(0:1,3,3) :: bcvel
    real(rp)        , intent(in), dimension(0:1,3)   :: bcpre
    logical         , intent(in), dimension(3)       :: is_forced
    integer         , intent(in), dimension(0:1,3)   :: nb
    logical         , intent(in), dimension(0:1,3)   :: is_bound
    real(rp)        , intent(in), dimension(3)       :: dli
    real(rp)        , intent(in), dimension(0:)      :: dzci,dzfi,dzci_g,dzfi_g
    logical :: passed
    !
    call chk_dims(ng,dims,passed);                 if(.not.passed) call abortit
    call chk_stop_type(stop_type,passed);          if(.not.passed) call abortit
    call chk_bc(cbcvel,cbcpre,bcvel,bcpre,passed); if(.not.passed) call abortit
    call chk_forcing(cbcpre,is_forced,passed);     if(.not.passed) call abortit
#if defined(_IMPDIFF_1D) && !defined(_IMPDIFF)
    if(myid == 0)  print*, 'ERROR: `_IMPDIFF_1D` cpp macro requires building with `_IMPDIFF` too.'; call abortit
#endif
#if defined(_IMPDIFF_1D) && !defined(_DECOMP_Z)
    if(myid == 0)  print*, 'WARNING: a run with implicit Z diffusion (`_IMPDIFF_1D`) is much more efficient &
                                     when combined with a Z-pencils parallelization (`_DECOMP_Z`).'
#endif
#if defined(_SINGLE_PRECISION_POISSON) && defined(IMPDIFF) && !(defined(IMPDIFF_1D) || defined(_DECOMP_Z))
    ! note: the code won't even compile as is now
    if(myid == 0)  print*, 'ERROR: a run with `_SINGLE_PRECISION_POISSON` can only accomodate implicit diffusion along Z,
                                   and requires building with building with `_IMPDIFF_1D` and `_DECOMP_Z`.'; call abortit
#endif
#if defined(_DEBUG) && !defined(_SINGLE_PRECISION) && !defined(_SINGLE_PRECISION_POISSON)
    call chk_solvers(ng,n,n_z,lo,hi,dli,dzci_g,dzfi_g,dzci,dzfi,nb,is_bound,cbcvel,cbcpre,bcvel,bcpre,passed)
    if(.not.passed) call abortit
#endif
  end subroutine test_sanity
  !
  subroutine chk_stop_type(stop_type,passed)
  implicit none
  logical, intent(in), dimension(3) :: stop_type
  logical, intent(out) :: passed
  passed = .true.
  if(.not.any(stop_type(:))) then
    if(myid == 0) print*, 'ERROR: stopping criterion not chosen.'
    passed = .false.
  end if
  end subroutine chk_stop_type
  !
  subroutine chk_dims(ng,dims,passed)
    implicit none
    integer, intent(in), dimension(3) :: ng
    integer, intent(in), dimension(2) :: dims
    logical, intent(out) :: passed
    integer, dimension(2) :: ii
    logical :: passed_loc
    passed = .true.
    !passed_loc = all(mod(ng(1:2),2) == 0)
    !if(myid == 0.and.(.not.passed_loc)) &
    !  print*, 'ERROR: itot and jtot should be even.'
    !passed = passed.and.passed_loc
#if !defined(_DECOMP_Y) && !defined(_DECOMP_Z)
    ii = [2,3]
#elif defined(_DECOMP_Y)
    ii = [1,3]
#elif defined(_DECOMP_Z)
    ii = [1,2]
#endif
    passed_loc = all(dims(:)<=ng(ii)).and.all(dims(:)>=1)
    if(myid == 0.and.(.not.passed_loc)) &
      print*, 'ERROR: 1 <= dims(:) <= [itot,jtot], or [itot,ktot], or [jtot ktot] depending on the decomposition.'
    passed = passed.and.passed_loc
  end subroutine chk_dims
  !
  subroutine chk_bc(cbcvel,cbcpre,bcvel,bcpre,passed)
  implicit none
  character(len=1), intent(in), dimension(0:1,3,3) :: cbcvel
  character(len=1), intent(in), dimension(0:1,3  ) :: cbcpre
  real(rp)        , intent(in), dimension(0:1,3,3) :: bcvel
  real(rp)        , intent(in), dimension(0:1,3  ) :: bcpre
  logical         , intent(out) :: passed
  character(len=2) :: bc01v,bc01p
  integer :: ivel,idir
  logical :: passed_loc
  passed = .true.
  !
  ! check validity of pressure and velocity BCs
  !
  passed_loc = .true.
  do ivel = 1,3
    do idir=1,3
      bc01v = cbcvel(0,idir,ivel)//cbcvel(1,idir,ivel)
      passed_loc = passed_loc.and.( (bc01v == 'PP').or. &
                                    (bc01v == 'ND').or. &
                                    (bc01v == 'DN').or. &
                                    (bc01v == 'NN').or. &
                                    (bc01v == 'DD') )
    end do
  end do
  if(myid == 0.and.(.not.passed_loc)) print*, 'ERROR: velocity BCs not valid.'
  passed = passed.and.passed_loc
  !
  passed_loc = .true.
  do idir=1,3
    bc01p = cbcpre(0,idir)//cbcpre(1,idir)
    passed_loc = passed_loc.and.( (bc01p == 'PP').or. &
                                  (bc01p == 'ND').or. &
                                  (bc01p == 'DN').or. &
                                  (bc01p == 'NN').or. &
                                  (bc01p == 'DD') )
  end do
  if(myid == 0.and.(.not.passed_loc)) print*, 'ERROR: pressure BCs not valid.'
  passed = passed.and.passed_loc
  !
  passed_loc = .true.
  do idir=1,3
    ivel = idir
    bc01v = cbcvel(0,idir,ivel)//cbcvel(1,idir,ivel)
    bc01p = cbcpre(0,idir)//cbcpre(1,idir)
    passed_loc = passed_loc.and.( (bc01v == 'PP'.and.bc01p == 'PP').or. &
                                  (bc01v == 'ND'.and.bc01p == 'DN').or. &
                                  (bc01v == 'DN'.and.bc01p == 'ND').or. &
                                  (bc01v == 'DD'.and.bc01p == 'NN').or. &
                                  (bc01v == 'NN'.and.bc01p == 'DD') )
  end do
  if(myid == 0.and.(.not.passed_loc)) print*, 'ERROR: velocity and pressure BCs not compatible.'
  passed = passed.and.passed_loc
  !
  passed_loc = .true.
  do idir=1,2
    passed_loc = passed_loc.and.((bcpre(0,idir) == 0.).and.(bcpre(1,idir) == 0.))
  end do
  if(myid == 0.and.(.not.passed_loc)) &
    print*, 'ERROR: pressure BCs in directions x and y must be homogeneous (value = 0.).'
  passed = passed.and.passed_loc
#if defined(_IMPDIFF)
  passed_loc = .true.
  do ivel = 1,3
    do idir=1,2
      bc01v = cbcvel(0,idir,ivel)//cbcvel(1,idir,ivel)
      passed_loc = passed_loc.and.(bc01v /= 'NN')
    end do
  end do
  if(myid == 0.and.(.not.passed_loc)) &
    print*, 'ERROR: Neumann-Neumann velocity BCs with implicit diffusion currently not supported in x and y; only in z.'
  passed = passed.and.passed_loc
  !
  passed_loc = .true.
  do ivel = 1,3
    do idir=1,2
      passed_loc = passed_loc.and.((bcvel(0,idir,ivel) == 0.).and.(bcvel(1,idir,ivel) == 0.))
    end do
  end do
  if(myid == 0.and.(.not.passed_loc)) &
    print*, 'ERROR: velocity BCs with implicit diffusion in directions x and y must be homogeneous (value = 0.).'
  passed = passed.and.passed_loc
#endif
  end subroutine chk_bc
  !
  subroutine chk_forcing(cbcpre,is_forced,passed)
  implicit none
  character(len=1), intent(in), dimension(0:1,3) :: cbcpre
  logical         , intent(in), dimension(3) :: is_forced
  logical         , intent(out) :: passed
  integer :: idir
  passed = .true.
  !
  ! 1) check for compatibility between pressure BCs and flow forcing
  !
  do idir=1,3
    if(is_forced(idir)) then
      passed = passed.and.(cbcpre(0,idir)//cbcpre(1,idir) == 'PP')
    end if
  end do
  if(myid == 0.and.(.not.passed)) &
  print*, 'ERROR: Flow cannot be forced in a non-periodic direction; check the BCs and is_forced in dns.in.'
  end subroutine chk_forcing
  !
  subroutine chk_solvers(ng,n,n_z,lo,hi,dli,dzci_g,dzfi_g,dzci,dzfi,nb,is_bound,cbcvel,cbcpre,bcvel,bcpre,passed)
    implicit none
    integer , intent(in), dimension(3) :: ng,n,n_z,lo,hi
    real(rp), intent(in), dimension(3) :: dli
    real(rp), intent(in), dimension(0:) :: dzci_g,dzfi_g,dzci,dzfi
    integer , intent(in), dimension(0:1,3) :: nb
    logical , intent(in), dimension(0:1,3) :: is_bound
    character(len=1), intent(in), dimension(0:1,3,3) :: cbcvel
    character(len=1), intent(in), dimension(0:1,3)   :: cbcpre
    real(rp), intent(in), dimension(0:1,3,3)          :: bcvel
    real(rp), intent(in), dimension(0:1,3)            :: bcpre
    logical , intent(out) :: passed
    real(rp), dimension(0:n(1)+1,0:n(2)+1,0:n(3)+1) :: u,v,w,p,up,vp,wp
    type(C_PTR), dimension(2,2) :: arrplan
    real(gp), dimension(n_z(1),n_z(2)) :: lambdaxy
    real(rp) :: normfft
    real(gp), dimension(n_z(3)) :: a,b,c,bb
    real(rp), dimension(n(2),n(3),0:1) :: rhsbx
    real(rp), dimension(n(1),n(3),0:1) :: rhsby
    real(rp), dimension(n(1),n(2),0:1) :: rhsbz
    real(rp), dimension(3) :: dl
    real(rp), dimension(0:n(3)+1) :: dzc,dzf
    real(rp) :: dt,dti,alpha
    real(rp) :: divtot,divmax,resmax
    logical :: passed_loc
    passed = .true.
    !
    ! initialize velocity below with some random noise
    !
    u(:,:,:) = 0.
    v(:,:,:) = 0.
    w(:,:,:) = 0.
    call add_noise(ng,lo,123,.5_rp,u(1:n(1),1:n(2),1:n(3)))
    call add_noise(ng,lo,456,.5_rp,v(1:n(1),1:n(2),1:n(3)))
    call add_noise(ng,lo,789,.5_rp,w(1:n(1),1:n(2),1:n(3)))
    !
    ! test pressure correction
    !
    call initsolver(ng,zstart,zend,dli,dzci_g,dzfi_g,cbcpre,bcpre(:,:),lambdaxy,['c','c','c'],a,b,c,arrplan,normfft, &
                    rhsbx,rhsby,rhsbz)
    dl  = dli**(-1)
    dzc = dzci**(-1)
    dzf = dzfi**(-1)
    dt  = acos(-1.) ! value is irrelevant
    dti = dt**(-1)
    call bounduvw(cbcvel,n,bcvel,nb,is_bound,.false.,dl,dzc,dzf,u,v,w)
    call fillps(n,dli,dzfi,dti,u,v,w,p)
    call updt_rhs_b(['c','c','c'],cbcpre,n,is_bound,rhsbx,rhsby,rhsbz,p)
    call solver(n,arrplan,normfft,lambdaxy,a,b,c,cbcpre(:,3),['c','c','c'],p)
    call boundp(cbcpre,n,bcpre,nb,is_bound,dl,dzc,p)
    call correc(n,dli,dzci,dt,p,u,v,w)
    call bounduvw(cbcvel,n,bcvel,nb,is_bound,.true.,dl,dzc,dzf,u,v,w)
    call chkdiv(lo,hi,dli,dzfi,u,v,w,divtot,divmax)
    passed_loc = divmax < small
    if(myid == 0.and.(.not.passed_loc)) &
    print*, 'ERROR: Pressure correction: Divergence is too large, with maximum = ', divmax
    passed = passed.and.passed_loc
    call fftend(arrplan)
#if defined(_IMPDIFF) && !defined(_IMPDIFF_1D)
    alpha = acos(-1.) ! irrelevant
    up(:,:,:) = 0.
    vp(:,:,:) = 0.
    wp(:,:,:) = 0.
    call add_noise(ng,lo,123,.5_rp,up(1:n(1),1:n(2),1:n(3)))
    call add_noise(ng,lo,456,.5_rp,vp(1:n(1),1:n(2),1:n(3)))
    call add_noise(ng,lo,789,.5_rp,wp(1:n(1),1:n(2),1:n(3)))
    call initsolver(ng,zstart,zend,dli,dzci_g,dzfi_g,cbcvel(:,:,1),bcvel(:,:,1),lambdaxy,['f','c','c'],a,b,c,arrplan,normfft, &
                    rhsbx,rhsby,rhsbz)
    call bounduvw(cbcvel,n,bcvel,nb,is_bound,.false.,dl,dzc,dzf,up,vp,wp)
    up(:,:,:) = up(:,:,:)*alpha
    u( :,:,:) = up(:,:,:)
    bb(:) = b(:) + alpha
    call updt_rhs_b(['f','c','c'],cbcvel(:,:,1),n,is_bound,rhsbx,rhsby,rhsbz,up)
    call solver(n,arrplan,normfft,lambdaxy,a,bb,c,cbcvel(:,3,1),['f','c','c'],up)
    call fftend(arrplan)
    call bounduvw(cbcvel,n,bcvel,nb,is_bound,.false.,dl,dzc,dzf,up,vp,wp) ! actually we are only interested in boundary condition in up
    call chk_helmholtz(lo,hi,dli,dzci,dzfi,alpha,u,up,cbcvel(:,:,1),is_bound,['f','c','c'],resmax)
    passed_loc = resmax < small
    if(myid == 0.and.(.not.passed_loc)) &
    print*, 'ERROR: wrong solution of Helmholtz equation in x direction.'
    passed = passed.and.passed_loc
    !
    call initsolver(ng,zstart,zend,dli,dzci_g,dzfi_g,cbcvel(:,:,2),bcvel(:,:,2),lambdaxy,['c','f','c'],a,b,c,arrplan,normfft, &
                    rhsbx,rhsby,rhsbz)
    call bounduvw(cbcvel,n,bcvel,nb,is_bound,.false.,dl,dzc,dzf,up,vp,wp)
    vp(:,:,:) = vp(:,:,:)*alpha
    v( :,:,:) = vp(:,:,:)
    bb(:) = b(:) + alpha
    call updt_rhs_b(['c','f','c'],cbcvel(:,:,2),n,is_bound,rhsbx,rhsby,rhsbz,vp)
    call solver(n,arrplan,normfft,lambdaxy,a,bb,c,cbcvel(:,3,2),['c','f','c'],vp)
    call fftend(arrplan)
    call bounduvw(cbcvel,n,bcvel,nb,is_bound,.false.,dl,dzc,dzf,up,vp,wp) ! actually we are only interested in boundary condition in vp
    call chk_helmholtz(lo,hi,dli,dzci,dzfi,alpha,v,vp,cbcvel(:,:,2),is_bound,['c','f','c'],resmax)
    passed_loc = resmax < small
    if(myid == 0.and.(.not.passed_loc)) &
    print*, 'ERROR: wrong solution of Helmholtz equation in y direction.'
    passed = passed.and.passed_loc
    !
    call initsolver(ng,zstart,zend,dli,dzci_g,dzfi_g,cbcvel(:,:,3),bcvel(:,:,3),lambdaxy,['c','c','f'],a,b,c,arrplan,normfft, &
                    rhsbx,rhsby,rhsbz)
    call bounduvw(cbcvel,n,bcvel,nb,is_bound,.false.,dl,dzc,dzf,up,vp,wp)
    wp(:,:,:) = wp(:,:,:)*alpha
    w( :,:,:) = wp(:,:,:)
    bb(:) = b(:) + alpha
    call updt_rhs_b(['c','c','f'],cbcvel(:,:,3),n,is_bound,rhsbx,rhsby,rhsbz,wp)
    call solver(n,arrplan,normfft,lambdaxy,a,bb,c,cbcvel(:,3,3),['c','c','f'],wp)
    call fftend(arrplan)
    call bounduvw(cbcvel,n,bcvel,nb,is_bound,.false.,dl,dzc,dzf,up,vp,wp) ! actually we are only interested in boundary condition in wp
    call chk_helmholtz(lo,hi,dli,dzci,dzfi,alpha,w,wp,cbcvel(:,:,3),is_bound,['c','c','f'],resmax)
    passed_loc = resmax < small
    if(myid == 0.and.(.not.passed_loc)) &
    print*, 'ERROR: wrong solution of Helmholtz equation in z direction.'
    passed = passed.and.passed_loc
#endif
  end subroutine chk_solvers
  !
  subroutine abortit
    implicit none
    if(myid == 0) print*, ''
    if(myid == 0) print*, '*** Simulation aborted due to errors in the input file ***'
    if(myid == 0) print*, '    check dns.in'
    call decomp_2d_finalize
    call MPI_FINALIZE(ierr)
    error stop
  end subroutine abortit
end module mod_sanity
