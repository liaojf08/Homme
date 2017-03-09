#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#define constLev 128
!#define _DBG_ print *,"File:",__FILE__," at ",__LINE__
!#define _DBG_ !DBG
!
!
module prim_advance_mod
  use edge_mod, only : EdgeBuffer_t
  use kinds, only : real_kind, iulog
  use perf_mod, only: t_startf, t_stopf, t_barrierf ! _EXTERNAL
  use parallel_mod, only : abortmp
use control_mod, only : north, south, east, west, neast, nwest, seast, swest
  implicit none
  private
  save
  public :: prim_advance_exp, prim_advance_si, prim_advance_init, preq_robert3,&
       applyCAMforcing_dynamics, applyCAMforcing, smooth_phis

  type (EdgeBuffer_t) :: edge1
  type (EdgeBuffer_t) :: edge2
  type (EdgeBuffer_t) :: edge3p1

  real (kind=real_kind) :: initialized_for_dt   = 0

  real (kind=real_kind), allocatable :: ur_weights(:)


contains

  subroutine prim_advance_init(integration)
    use edge_mod, only : initEdgeBuffer
    use dimensions_mod, only : nlev
    use control_mod, only : qsplit,rsplit
    character(len=*)    , intent(in) :: integration
    integer :: i

    if (rsplit==0) then
       call initEdgeBuffer(edge3p1,3*nlev+1)
    else
       ! need extra buffer space for dp3d
       call initEdgeBuffer(edge3p1,4*nlev+1)
    endif

    if(integration == 'semi_imp') then
       call initEdgeBuffer(edge1,nlev)
       call initEdgeBuffer(edge2,2*nlev)
    end if

    ! compute averaging weights for RK+LF (tstep_type=1) timestepping:
    allocate(ur_weights(qsplit))
    ur_weights(:)=0.0d0

    if(mod(qsplit,2).NE.0)then
       ur_weights(1)=1.0d0/qsplit
       do i=3,qsplit,2
         ur_weights(i)=2.0d0/qsplit
       enddo
    else
       do i=2,qsplit,2
         ur_weights(i)=2.0d0/qsplit
       enddo
    endif

  end subroutine prim_advance_init




  subroutine prim_advance_exp(elem, deriv, hvcoord, hybrid,&
       dt, tl,  nets, nete, compute_diagnostics)
    use bndry_mod, only : bndry_exchangev
    use control_mod, only : prescribed_wind, qsplit, tstep_type, rsplit, qsplit, moisture
    use derivative_mod, only : derivative_t, vorticity, divergence, gradient, gradient_wk
    use dimensions_mod, only : np, nlev, nlevp
    use edge_mod, only : EdgeBuffer_t, edgevpack, edgevunpack, initEdgeBuffer
    use element_mod, only : element_t
    use filter_mod, only : filter_t
    use hybvcoord_mod, only : hvcoord_t
    use hybrid_mod, only : hybrid_t
    use reduction_mod, only : reductionbuffer_ordered_1d_t
    use time_mod, only : TimeLevel_t,  timelevel_qdp
    use diffusion_mod, only :  prim_diffusion

#ifndef CAM
    use asp_tests, only : asp_advection_vertical
#endif

    implicit none

    type (element_t), intent(inout), target   :: elem(:)
    type (derivative_t), intent(in)   :: deriv
    type (hvcoord_t)                  :: hvcoord

    type (hybrid_t)    , intent(in):: hybrid

    real (kind=real_kind), intent(in) :: dt
    type (TimeLevel_t)   , intent(in) :: tl
    integer              , intent(in) :: nets
    integer              , intent(in) :: nete
    logical, intent(in)               :: compute_diagnostics

    ! =================
    ! Local
    ! =================
    real (kind=real_kind) ::  dt2, time, dt_vis, eta_ave_w
    real (kind=real_kind) ::  eta_dot_dpdn(np,np,nlevp)
    real (kind=real_kind) ::  dp(np,np)
    integer :: ie,nm1,n0,np1,nstep,method,qsplit_stage,k, qn0


    call t_barrierf('sync_prim_advance_exp', hybrid%par%comm)
    call t_startf('prim_advance_exp')
    nm1   = tl%nm1
    n0    = tl%n0
    np1   = tl%np1
    nstep = tl%nstep

    ! timelevel to use for accessing Qdp() to compute virtual temperature
    qn0 = -1    ! -1 = disabled (assume dry dynamics)
    if ( moisture /= "dry") then
       call TimeLevel_Qdp(tl, qsplit, qn0)  ! compute current Qdp() timelevel
    endif



!
!   tstep_type=0  pure leapfrog except for very first timestep   CFL=1
!                    typically requires qsplit=4 or 5
!   tstep_type=1  RK2 followed by qsplit-1 leapfrog steps        CFL=close to qsplit
!                    typically requires qsplit=4 or 5
!   tstep_type=2  RK2-SSP 3 stage (as used by tracers)           CFL=.58
!                    optimal in terms of SSP CFL, but not        CFLSSP=2
!                    optimal in terms of CFL
!                    typically requires qsplit=3
!                    but if windspeed > 340m/s, could use this
!                    with qsplit=1
!   tstep_type=3  classic RK3                                    CFL=1.73 (sqrt(3))
!
!   tstep_type=4  Kinnmark&Gray RK4 4 stage                      CFL=sqrt(8)
!                 should we replace by standard RK4 (CFL=sqrt(8))?
!   tstep_type=5  Kinnmark&Gray RK3 5 stage 3rd order            CFL=3.87  (sqrt(15))
!                    optimal: for windspeeds ~120m/s,gravity: 340m/2
!                    run with qsplit=1
!

    if(tstep_type==0)then
       method=0                ! pure leapfrog
       if (nstep==0) method=1  ! always use RK2 for first timestep
    else if (tstep_type==1) then
       method=0                           ! LF
       qsplit_stage = mod(nstep,qsplit)
       if (qsplit_stage==0) method=1      ! RK2 on first of qsplit steps
    else if (tstep_type>1) then
       method = tstep_type                ! other RK variants
    endif


    ! weights for computing mean dynamics fluxes
    eta_ave_w = 1d0/qsplit
    if(tstep_type==1)then
       ! RK2 + LF scheme has tricky weights:
       eta_ave_w=ur_weights(qsplit_stage+1)
    endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
    ! fix dynamical variables, skip dynamics
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
    if (1==prescribed_wind) then
       time=tl%nstep*dt
       do ie=nets,nete
#ifdef CAM
          ! read in omega met data and store it in derived%omega_prescribed
          ! --- compute eta_dot_dpdn from omega_prescribed ...
          eta_dot_dpdn(:,:,nlev+1) = 0.0d0

          do k = nlev,2,-1
             eta_dot_dpdn(:,:,k) = eta_dot_dpdn(:,:,k+1) - 2.d0*elem(ie)%derived%omega_prescribed(:,:,k)
          enddo

          eta_dot_dpdn(:,:,1) = 0.0d0
#else
          ! assume most fields are constant in time
          elem(ie)%state%ps_v(:,:,np1) = elem(ie)%state%ps_v(:,:,n0)
          elem(ie)%state%lnps(:,:,np1) = elem(ie)%state%lnps(:,:,n0)
          elem(ie)%state%T(:,:,:,np1)   = elem(ie)%state%T(:,:,:,n0)
          elem(ie)%state%v(:,:,:,:,np1) = elem(ie)%state%v(:,:,:,:,n0)
          elem(ie)%derived%div = 0
	  !get eta_dot_dpdn at n0, ps_v is not used or calculated in this routine
          call asp_advection_vertical(time,hvcoord,elem(ie)%state%ps_v(:,:,n0),&
              eta_dot_dpdn)
#endif
          ! accumulate mean fluxes for advection
          elem(ie)%derived%eta_dot_dpdn(:,:,:) = &
               elem(ie)%derived%eta_dot_dpdn(:,:,:) + eta_dot_dpdn(:,:,:)*eta_ave_w

          ! subcycling code uses a mean flux to advect tracers
#if (defined ELEMENT_OPENMP)
          !$omp parallel do private(k,dp)
#endif
          do k=1,nlev
             dp(:,:) =&
                  ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
                  ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(:,:,tl%n0)

             elem(ie)%derived%vn0(:,:,1,k)=elem(ie)%derived%vn0(:,:,1,k)+&
                  eta_ave_w*elem(ie)%state%v(:,:,1,k,n0)*dp(:,:)
             elem(ie)%derived%vn0(:,:,2,k)=elem(ie)%derived%vn0(:,:,2,k)+&
                  eta_ave_w*elem(ie)%state%v(:,:,2,k,n0)*dp(:,:)
          enddo
       end do
       call t_stopf('prim_advance_exp')
       return
    endif




    ! ==================================
    ! Take timestep
    ! ==================================
    dt_vis = dt
    if (method==0) then
       ! regular LF step
       dt2 = 2*dt
       call compute_and_apply_rhs(np1,nm1,n0,qn0,dt2,elem,hvcoord,hybrid,&
            deriv,nets,nete,compute_diagnostics,eta_ave_w)
       dt_vis = dt2  ! dt to use for time-split dissipation
    else if (method==1) then
       ! RK2
       ! forward euler to u(dt/2) = u(0) + (dt/2) RHS(0)  (store in u(np1))
       call compute_and_apply_rhs(np1,n0,n0,qn0,dt/2,elem,hvcoord,hybrid,&
            deriv,nets,nete,compute_diagnostics,0d0)
       ! leapfrog:  u(dt) = u(0) + dt RHS(dt/2)     (store in u(np1))
       call compute_and_apply_rhs(np1,n0,np1,qn0,dt,elem,hvcoord,hybrid,&
            deriv,nets,nete,.false.,eta_ave_w)
    else if (method==2) then
       ! RK2-SSP 3 stage.  matches tracer scheme. optimal SSP CFL, but
       ! not optimal for regular CFL
       ! u1 = u0 + dt/2 RHS(u0)
       call compute_and_apply_rhs(np1,n0,n0,qn0,dt/2,elem,hvcoord,hybrid,&
            deriv,nets,nete,compute_diagnostics,eta_ave_w/3)
       ! u2 = u1 + dt/2 RHS(u1)
       call compute_and_apply_rhs(np1,np1,np1,qn0,dt/2,elem,hvcoord,hybrid,&
            deriv,nets,nete,.false.,eta_ave_w/3)
       ! u3 = u2 + dt/2 RHS(u2)
       call compute_and_apply_rhs(np1,np1,np1,qn0,dt/2,elem,hvcoord,hybrid,&
            deriv,nets,nete,.false.,eta_ave_w/3)
       ! unew = u/3 +2*u3/3  = u + 1/3 (RHS(u) + RHS(u1) + RHS(u2))
       do ie=nets,nete
          elem(ie)%state%v(:,:,:,:,np1)= elem(ie)%state%v(:,:,:,:,n0)/3 &
               + 2*elem(ie)%state%v(:,:,:,:,np1)/3
          elem(ie)%state%T(:,:,:,np1)= elem(ie)%state%T(:,:,:,n0)/3 &
               + 2*elem(ie)%state%T(:,:,:,np1)/3
          if (rsplit==0) then
             elem(ie)%state%ps_v(:,:,np1)= elem(ie)%state%ps_v(:,:,n0)/3 &
                  + 2*elem(ie)%state%ps_v(:,:,np1)/3
          else
             elem(ie)%state%dp3d(:,:,:,np1)= elem(ie)%state%dp3d(:,:,:,n0)/3 &
                  + 2*elem(ie)%state%dp3d(:,:,:,np1)/3
          endif
       enddo
    else if (method==3) then
       ! classic RK3  CFL=sqrt(3)
       ! u1 = u0 + dt/3 RHS(u0)
       call compute_and_apply_rhs(np1,n0,n0,qn0,dt/3,elem,hvcoord,hybrid,&
            deriv,nets,nete,compute_diagnostics,0d0)
       ! u2 = u0 + dt/2 RHS(u1)
       call compute_and_apply_rhs(np1,n0,np1,qn0,dt/2,elem,hvcoord,hybrid,&
            deriv,nets,nete,.false.,0d0)
       ! u3 = u0 + dt RHS(u2)
       call compute_and_apply_rhs(np1,n0,np1,qn0,dt,elem,hvcoord,hybrid,&
            deriv,nets,nete,.false.,eta_ave_w)
    else if (method==4) then
       ! KG 4th order 4 stage:   CFL=sqrt(8)
       ! low storage version of classic RK4
       ! u1 = u0 + dt/4 RHS(u0)
       call compute_and_apply_rhs(np1,n0,n0,qn0,dt/4,elem,hvcoord,hybrid,&
            deriv,nets,nete,compute_diagnostics,0d0)
       ! u2 = u0 + dt/3 RHS(u1)
       call compute_and_apply_rhs(np1,n0,np1,qn0,dt/3,elem,hvcoord,hybrid,&
            deriv,nets,nete,.false.,0d0)
       ! u3 = u0 + dt/2 RHS(u2)
       call compute_and_apply_rhs(np1,n0,np1,qn0,dt/2,elem,hvcoord,hybrid,&
            deriv,nets,nete,.false.,0d0)
       ! u4 = u0 + dt RHS(u3)
       call compute_and_apply_rhs(np1,n0,np1,qn0,dt,elem,hvcoord,hybrid,&
            deriv,nets,nete,.false.,eta_ave_w)
    else if (method==5) then
       ! KG 3nd order 5 stage:   CFL=sqrt( 4^2 -1) = 3.87
       ! u1 = u0 + dt/5 RHS(u0)
       call compute_and_apply_rhs(np1,n0,n0,qn0,dt/5,elem,hvcoord,hybrid,&
            deriv,nets,nete,compute_diagnostics,0d0)
       ! u2 = u0 + dt/5 RHS(u1)
       call compute_and_apply_rhs(np1,n0,np1,qn0,dt/5,elem,hvcoord,hybrid,&
            deriv,nets,nete,.false.,0d0)
       ! u3 = u0 + dt/3 RHS(u2)
       call compute_and_apply_rhs(np1,n0,np1,qn0,dt/3,elem,hvcoord,hybrid,&
            deriv,nets,nete,.false.,0d0)
       ! u4 = u0 + dt/2 RHS(u3)
       call compute_and_apply_rhs(np1,n0,np1,qn0,dt/2,elem,hvcoord,hybrid,&
            deriv,nets,nete,.false.,0d0)
       ! u5 = u0 + dt RHS(u4)
       call compute_and_apply_rhs(np1,n0,np1,qn0,dt,elem,hvcoord,hybrid,&
            deriv,nets,nete,.false.,eta_ave_w)
    else
       call abortmp('ERROR: bad choice of tstep_type')
    endif



    ! ==============================================
    ! Time-split Horizontal diffusion: nu.del^2 or nu.del^4
    ! U(*) = U(t+1)  + dt2 * HYPER_DIFF_TERM(t+1)
    ! ==============================================
#ifdef ENERGY_DIAGNOSTICS
    if (compute_diagnostics) then
       do ie = nets,nete
          elem(ie)%accum%DIFF(:,:,:,:)=elem(ie)%state%v(:,:,:,:,np1)
          elem(ie)%accum%DIFFT(:,:,:)=elem(ie)%state%T(:,:,:,np1)
       enddo
    endif
#endif

    ! note:time step computes u(t+1)= u(t*) + RHS.
    ! for consistency, dt_vis = t-1 - t*, so this is timestep method dependent
    if (tstep_type==0) then
       ! leapfrog special case
       call advance_hypervis_lf(edge3p1,elem,hvcoord,hybrid,deriv,nm1,n0,np1,nets,nete,dt_vis)
    else
       if (rsplit==0) then
          ! forward-in-time, maybe hypervis applied to PS
          call advance_hypervis(edge3p1,elem,hvcoord,hybrid,deriv,np1,nets,nete,dt_vis,eta_ave_w)
       else
          ! forward-in-time, hypervis applied to dp3d
          call advance_hypervis_dp(edge3p1,elem,hvcoord,hybrid,deriv,np1,nets,nete,dt_vis,eta_ave_w)
       endif
    endif


#ifdef ENERGY_DIAGNOSTICS
    if (compute_diagnostics) then
       do ie = nets,nete
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k)
#endif
         do k=1,nlev  !  Loop index added (AAM)
          elem(ie)%accum%DIFF(:,:,:,k)=( elem(ie)%state%v(:,:,:,k,np1) -&
               elem(ie)%accum%DIFF(:,:,:,k) ) / dt_vis
          elem(ie)%accum%DIFFT(:,:,k)=( elem(ie)%state%T(:,:,k,np1) -&
               elem(ie)%accum%DIFFT(:,:,k) ) / dt_vis
         enddo
       enddo
    endif
#endif
    call t_stopf('prim_advance_exp')
    end subroutine prim_advance_exp



subroutine prim_advance_si(elem, nets, nete, cg, blkjac, red, &
          refstate, hvcoord, deriv, flt, hybrid, tl, dt)
       use bndry_mod, only : bndry_exchangev
       use cg_mod, only : cg_t, cg_create
       use control_mod, only : filter_freq,debug_level, precon_method
       use derivative_mod, only : derivative_t, vorticity, divergence, gradient, gradient_wk
       use dimensions_mod, only : np, nlev, nlevp
       use edge_mod, only : EdgeBuffer_t, edgevpack, edgevunpack, initEdgeBuffer
       use element_mod, only : element_t
       use filter_mod, only : filter_t, preq_filter
       use hybvcoord_mod, only : hvcoord_t
       use hybrid_mod, only : hybrid_t
       use prim_si_ref_mod, only : ref_state_t, set_vert_struct_mat
       use reduction_mod, only : reductionbuffer_ordered_1d_t
       use solver_mod, only : pcg_solver, blkjac_t, blkjac_init
       use time_mod, only : TimeLevel_t
       use prim_si_mod, only : preq_vertadv, preq_omegap, preq_pressure
       use diffusion_mod, only :  prim_diffusion
       use physical_constants, only : kappa, rrearth, rgas, cp, rwater_vapor
       use physics_mod, only : virtual_temperature, virtual_specific_heat
       implicit none

       integer, intent(in)               :: nets,nete
       type (element_t), intent(inout), target :: elem(:)
       type (blkjac_t)                   :: blkjac(nets:nete)

       type (cg_t)                       :: cg

       type (ReductionBuffer_ordered_1d_t), intent(inout) :: red

       type (ref_state_t), intent(in), target :: refstate
       type (hvcoord_t), intent(in)      :: hvcoord
       type (derivative_t), intent(in)   :: deriv
       type (filter_t), intent(in)       :: flt
       type (hybrid_t), intent(in)       :: hybrid
       type (TimeLevel_t), intent(in)    :: tl
       real(kind=real_kind), intent(in)  :: dt
       real(kind=real_kind)              :: time_adv
       ! ==========================
       ! Local variables...
       ! ==========================

       real(kind=real_kind)                           :: ps0
       real(kind=real_kind)                           :: psref

       real(kind=real_kind), dimension(np,np)         :: ps
       real(kind=real_kind), dimension(np,np)         :: rps
       real(kind=real_kind), dimension(np,np,nlev)    :: rpmid
       real(kind=real_kind), dimension(np,np,nlev)    :: omegap
       real(kind=real_kind), dimension(np,np,nlev)    :: rpdel

       real(kind=real_kind) :: pintref(nlevp)
       real(kind=real_kind) :: pdelref(nlev)
       real(kind=real_kind) :: pmidref(nlev)
       real(kind=real_kind) :: rpdelref(nlev)
       real(kind=real_kind) :: rpmidref(nlev)

       real(kind=real_kind) :: pint(np,np,nlevp)
       real(kind=real_kind) :: pdel(np,np,nlev)
       real(kind=real_kind) :: pmid(np,np,nlev)

       real(kind=real_kind), dimension(np,np,nlevp) :: eta_dot_dp_deta
       real(kind=real_kind), dimension(np,np,nlev)  :: vgrad_ps

       real(kind=real_kind), dimension(np,np,nlev)   :: T_vadv
       real(kind=real_kind), dimension(np,np,2,nlev) :: v_vadv

       real(kind=real_kind), dimension(np,np)      :: HT
       real(kind=real_kind), dimension(np,np)      :: HrefT
       real(kind=real_kind), dimension(np,np)      :: HrefTm1

       real(kind=real_kind), dimension(np,np)      :: Gref0
       real(kind=real_kind), dimension(np,np)      :: Grefm1
       real(kind=real_kind), dimension(np,np)      :: E
       real(kind=real_kind), dimension(np,np)      :: Phi
       real(kind=real_kind), dimension(np,np)      :: dGref

       real(kind=real_kind), dimension(np,np,2)    :: vco
       real(kind=real_kind), dimension(np,np,2)    :: gradT
       real(kind=real_kind), dimension(np,np,2)    :: grad_Phi

       real(kind=real_kind), dimension(:,:), pointer  :: Emat
       real(kind=real_kind), dimension(:,:), pointer  :: Emat_inv
       real(kind=real_kind), dimension(:,:), pointer  :: Amat
       real(kind=real_kind), dimension(:,:), pointer  :: Amat_inv
       real(kind=real_kind), dimension(:), pointer    :: Lambda

       real(kind=real_kind), dimension(:), pointer    :: Tref
       real(kind=real_kind), dimension(:), pointer    :: RTref
       real(kind=real_kind), dimension(:), pointer    :: Pvec
       real(kind=real_kind), dimension(:,:), pointer  :: Href
       real(kind=real_kind), dimension(:,:), pointer  :: Tmat

       real(kind=real_kind) :: Vscript(np,np,2,nlev,nets:nete)
       real(kind=real_kind) :: Tscript(np,np,nlev,nets:nete)
       real(kind=real_kind) :: Pscript(np,np,nets:nete)

       real(kind=real_kind), dimension(np,np)      :: HrefTscript
       real(kind=real_kind), dimension(np,np)      :: suml
       real(kind=real_kind), dimension(np,np,2)    :: gVscript
       real(kind=real_kind), dimension(np,np,nlev) :: div_Vscript

       real(kind=real_kind) :: B(np,np,nlev,nets:nete)
       real(kind=real_kind) :: C(np,np,nlev,nets:nete)
       real(kind=real_kind) :: D(np,np,nlev,nets:nete)

       real(kind=real_kind) :: Gamma_ref(np,np,nlev,nets:nete)

       real(kind=real_kind) :: Gref(np,np,nlev,nets:nete)
       real(kind=real_kind) :: grad_dGref(np,np,2,nlev)
       real(kind=real_kind) :: grad_Gref(np,np,2,nlev)

       real(kind=real_kind) :: div(np,np)
       real(kind=real_kind) :: gv(np,np,2)

       real(kind=real_kind) :: dt2
       real(kind=real_kind) :: rpsref
       real(kind=real_kind) :: rdt
       real(kind=real_kind) :: hkk, hkl
       real(kind=real_kind) :: ddiv

       real(kind=real_kind) :: vgradT
       real(kind=real_kind) :: hybfac
       real(kind=real_kind) :: Crkk
       real(kind=real_kind) :: v1,v2
       real(kind=real_kind) :: term

       real(kind=real_kind) :: Vs1,Vs2
       real(kind=real_kind) :: glnps1, glnps2
       real(kind=real_kind) :: gGr1,gGr2

       real (kind=real_kind),allocatable :: solver_wts(:,:)  ! solver weights array for nonstaggered grid

       integer              :: nm1,n0,np1,nfilt
       integer              :: nstep
       integer              :: i,j,k,l,ie,kptr

       call t_barrierf('sync_prim_advance_si', hybrid%par%comm)
       call t_startf('prim_advance_si')

       nm1   = tl%nm1
       n0    = tl%n0
       np1   = tl%np1
       nstep = tl%nstep


       if ( dt /= initialized_for_dt ) then
          if(hybrid%par%masterproc) print *,'Initializing semi-implicit matricies for dt=',dt

#if (! defined ELEMENT_OPENMP)
          !$OMP MASTER
#endif
          call set_vert_struct_mat(dt, refstate, hvcoord, hybrid%masterthread)
#if (! defined ELEMENT_OPENMP)
          !$OMP END MASTER
#endif

          allocate(solver_wts(np*np,nete-nets+1))
          do ie=nets,nete
             kptr=1
             do j=1,np
                do i=1,np

                   ! so this code is BFB  with old code.  should change to simpler formula below
                   solver_wts(kptr,ie-nets+1) = 1d0/nint(1d0/(elem(ie)%mp(i,j)*elem(ie)%rmp(i,j)))
                   !solver_wts(kptr,ie-nets+1) = elem(ie)%mp(i,j)*elem(ie)%rmp(i,j)

                   kptr=kptr+1
                end do
             end do
          end do
          call cg_create(cg, np*np, nlev, nete-nets+1, hybrid, debug_level, solver_wts)
          deallocate(solver_wts)
          if (precon_method == "block_jacobi") then
             call blkjac_init(elem, deriv,refstate%Lambda,nets,nete,blkjac)
          end if
          initialized_for_dt = dt
       endif


       nfilt = tl%nm1     ! time level at which filter is applied (time level n)
       dt2   = 2.0_real_kind*dt
       rdt   = 1.0_real_kind/dt

       ps0      = hvcoord%ps0
       psref    = refstate%psr

       Emat     => refstate%Emat
       Emat_inv => refstate%Emat_inv
       Amat     => refstate%Amat
       Amat_inv => refstate%Amat_inv
       Lambda   => refstate%Lambda

       RTref    => refstate%RTref
       Tref     => refstate%Tref
       Href     => refstate%Href
       Tmat     => refstate%Tmat
       Pvec     => refstate%Pvec

       ! ============================================================
       ! If the time is right, apply a filter to the state variables
       ! ============================================================

       if (nstep > 0 .and. filter_freq > 0 .and. MODULO(nstep,filter_freq) == 0 ) then
          call preq_filter(elem, edge3p1, flt, cg%hybrid, nfilt, nets, nete)
       end if

       do ie = nets, nete

          elem(ie)%derived%grad_lnps(:,:,:) = gradient(elem(ie)%state%lnps(:,:,n0),deriv)*rrearth

       end do

       ! ================================================
       ! boundary exchange grad_lnps
       ! ================================================

       do ie = nets, nete

          do k=1,nlevp
             pintref(k)  = hvcoord%hyai(k)*ps0 + hvcoord%hybi(k)*psref
          end do

          do k=1,nlev
             pmidref(k)  = hvcoord%hyam(k)*ps0 + hvcoord%hybm(k)*psref
             pdelref(k)  = pintref(k+1) - pintref(k)
             rpmidref(k) = 1.0_real_kind/pmidref(k)
             rpdelref(k) = 1.0_real_kind/pdelref(k)
          end do

          rpsref   = 1.0_real_kind/psref

          ps(:,:) = EXP(elem(ie)%state%lnps(:,:,n0))
          rps(:,:) = 1.0_real_kind/ps(:,:)

          call preq_pressure(ps0,ps,hvcoord%hyai,hvcoord%hybi,hvcoord%hyam,hvcoord%hybm,pint,pmid,pdel)

          rpmid = 1.0_real_kind/pmid
          rpdel = 1.0_real_kind/pdel

#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k,i,j,v1,v2,vco)
#endif
          do k=1,nlev
             do j=1,np
                do i=1,np
                   v1     = elem(ie)%state%v(i,j,1,k,n0)
                   v2     = elem(ie)%state%v(i,j,2,k,n0)

                   ! Contravariant velocities

                   vco(i,j,1) = elem(ie)%Dinv(1,1,i,j)*v1 + elem(ie)%Dinv(1,2,i,j)*v2
                   vco(i,j,2) = elem(ie)%Dinv(2,1,i,j)*v1 + elem(ie)%Dinv(2,2,i,j)*v2

                   vgrad_ps(i,j,k) = ps(i,j)*(vco(i,j,1)*elem(ie)%derived%grad_lnps(i,j,1) + &
                        vco(i,j,2)*elem(ie)%derived%grad_lnps(i,j,2))

                end do
             end do
          end do

          call preq_omegap(elem(ie)%derived%div(:,:,:,n0),vgrad_ps,pdel,rpmid, &
               hvcoord%hybm,hvcoord%hybd,elem(ie)%derived%omega_p)

          Pscript(:,:,ie)        = 0.0_real_kind
          eta_dot_dp_deta(:,:,1) = 0.0_real_kind

          do k=1,nlev
             do j=1,np
                do i=1,np
                   eta_dot_dp_deta(i,j,k+1) = eta_dot_dp_deta(i,j,k) + &
                        vgrad_ps(i,j,k)*hvcoord%hybd(k) + elem(ie)%derived%div(i,j,k,n0)*pdel(i,j,k)
                   ddiv = elem(ie)%derived%div(i,j,k,n0) - 0.5_real_kind*elem(ie)%derived%div(i,j,k,nm1)
                   Pscript(i,j,ie) = Pscript(i,j,ie) + ddiv*pdelref(k)
                end do
             end do
          end do

          do j=1,np
             do i=1,np
                Pscript(i,j,ie) = elem(ie)%state%lnps(i,j,nm1) + &
                     dt2*( rpsref*Pscript(i,j,ie) - rps(i,j)*eta_dot_dp_deta(i,j,nlev+1) )
             end do
          end do

#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k,i,j)
#endif
          do k=1,nlev-1
             do j=1,np
                do i=1,np
                   eta_dot_dp_deta(i,j,k+1) = hvcoord%hybi(k+1)*eta_dot_dp_deta(i,j,nlev+1) - &
                        eta_dot_dp_deta(i,j,k+1)
                end do
             end do
          end do

          eta_dot_dp_deta(:,:,nlev+1) = 0.0_real_kind

          call preq_vertadv(elem(ie)%state%T(:,:,:,n0),elem(ie)%state%v(:,:,:,:,n0), &
               eta_dot_dp_deta,rpdel,T_vadv,v_vadv)

          suml(:,:) = 0.0_real_kind

          do k=1,nlev

             gradT(:,:,:) = gradient(elem(ie)%state%T(:,:,k,n0),deriv)*rrearth
             Crkk       = 0.5_real_kind

             do j=1,np
                do i=1,np
                   term = Crkk*(elem(ie)%derived%div(i,j,k,n0) - &
                        0.5_real_kind*elem(ie)%derived%div(i,j,k,nm1))*pdelref(k)
                   suml(i,j)  = suml(i,j) + term

                   v1     = elem(ie)%state%v(i,j,1,k,n0)
                   v2     = elem(ie)%state%v(i,j,2,k,n0)

                   ! Contravariant velocities

                   vco(i,j,1) = elem(ie)%Dinv(1,1,i,j)*v1 + elem(ie)%Dinv(1,2,i,j)*v2
                   vco(i,j,2) = elem(ie)%Dinv(2,1,i,j)*v1 + elem(ie)%Dinv(2,2,i,j)*v2

                   vgradT = vco(i,j,1)*gradT(i,j,1) + vco(i,j,2)*gradT(i,j,2)

                   Tscript(i,j,k,ie) = elem(ie)%state%T(i,j,k,nm1) &
                        + dt2*(- vgradT - T_vadv(i,j,k)           &
                        + kappa*(elem(ie)%state%T(i,j,k,n0)*elem(ie)%derived%omega_p(i,j,k) &
                        + Tref(k)*rpmidref(k)*suml(i,j)))
                   suml(i,j)  = suml(i,j) + term
                end do
             end do
          end do

          HrefT(:,:)   = 0.0_real_kind
          HrefTm1(:,:) = 0.0_real_kind
          HT(:,:)      = 0.0_real_kind

          do k=nlev,1,-1

             do j=1,np
                do i=1,np
                   hkl = rpmidref(k)*pdelref(k)
                   hkk = hkl*0.5_real_kind
                   Gref0(i,j)   = HrefT(i,j) + Rgas*hkk*elem(ie)%state%T(i,j,k,n0)
                   HrefT(i,j)   = HrefT(i,j) + Rgas*hkl*elem(ie)%state%T(i,j,k,n0)
                   Grefm1(i,j)  = HrefTm1(i,j) + Rgas*hkk*elem(ie)%state%T(i,j,k,nm1)
                   HrefTm1(i,j) = HrefTm1(i,j) + Rgas*hkl*elem(ie)%state%T(i,j,k,nm1)
                   hkl = rpmid(i,j,k)*pdel(i,j,k)
                   hkk = hkl*0.5_real_kind
                   Phi(i,j) = HT(i,j) + Rgas*hkk*elem(ie)%state%T(i,j,k,n0)
                   HT(i,j)  = HT(i,j) + Rgas*hkl*elem(ie)%state%T(i,j,k,n0)
                end do
             end do

             do j=1,np
                do i=1,np
                   v1     = elem(ie)%state%v(i,j,1,k,n0)
                   v2     = elem(ie)%state%v(i,j,2,k,n0)

                   ! covariant velocity

                   vco(i,j,1) = elem(ie)%D(1,1,i,j)*v1 + elem(ie)%D(2,1,i,j)*v2
                   vco(i,j,2) = elem(ie)%D(1,2,i,j)*v1 + elem(ie)%D(2,2,i,j)*v2

                   E(i,j) = 0.5_real_kind*( v1*v1 + v2*v2 )

                   Gref0(i,j)  =  Gref0(i,j)  + elem(ie)%state%phis(i,j) + RTref(k)*elem(ie)%state%lnps(i,j,n0)
                   Grefm1(i,j) =  Grefm1(i,j) + elem(ie)%state%phis(i,j) + RTref(k)*elem(ie)%state%lnps(i,j,nm1)

                   Phi(i,j)    =  Phi(i,j) + E(i,j) + elem(ie)%state%phis(i,j)
                   dGref(i,j)  =  -(Gref0(i,j)  - 0.5_real_kind*Grefm1(i,j))
                end do
             end do

             elem(ie)%derived%zeta(:,:,k) = vorticity(vco,deriv)*rrearth
             grad_Phi(:,:,:)     = gradient(Phi,deriv)*rrearth
             grad_dGref(:,:,:,k) = gradient_wk(dGref,deriv)*rrearth

             do j=1,np
                do i=1,np

                   elem(ie)%derived%zeta(i,j,k) = elem(ie)%rmetdet(i,j)*elem(ie)%derived%zeta(i,j,k)
                   hybfac =  hvcoord%hybm(k)*(ps(i,j)*rpmid(i,j,k))

                   glnps1 = elem(ie)%Dinv(1,1,i,j)*elem(ie)%derived%grad_lnps(i,j,1) + &
                        elem(ie)%Dinv(2,1,i,j)*elem(ie)%derived%grad_lnps(i,j,2)
                   glnps2 = elem(ie)%Dinv(1,2,i,j)*elem(ie)%derived%grad_lnps(i,j,1) + &
                        elem(ie)%Dinv(2,2,i,j)*elem(ie)%derived%grad_lnps(i,j,2)

                   v1 = elem(ie)%Dinv(1,1,i,j)*grad_Phi(i,j,1) + elem(ie)%Dinv(2,1,i,j)*grad_Phi(i,j,2)
                   v2 = elem(ie)%Dinv(1,2,i,j)*grad_Phi(i,j,1) + elem(ie)%Dinv(2,2,i,j)*grad_Phi(i,j,2)

                   Vscript(i,j,1,k,ie) = - v_vadv(i,j,1,k) &
                        + elem(ie)%state%v(i,j,2,k,n0) * (elem(ie)%fcor(i,j) + elem(ie)%derived%zeta(i,j,k)) &
                        - v1 - Rgas*hybfac*elem(ie)%state%T(i,j,k,n0)*glnps1

                   Vscript(i,j,2,k,ie) = - v_vadv(i,j,2,k) &
                        - elem(ie)%state%v(i,j,1,k,n0) * (elem(ie)%fcor(i,j) + elem(ie)%derived%zeta(i,j,k)) &
                        - v2 - Rgas*hybfac*elem(ie)%state%T(i,j,k,n0)*glnps2

                end do
             end do

          end do

#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k,i,j,Vs1,Vs2)
#endif
          do k=1,nlev

             do j=1,np
                do i=1,np
                   Vs1 = elem(ie)%Dinv(1,1,i,j)*grad_dGref(i,j,1,k) + elem(ie)%Dinv(2,1,i,j)*grad_dGref(i,j,2,k)
                   Vs2 = elem(ie)%Dinv(1,2,i,j)*grad_dGref(i,j,1,k) + elem(ie)%Dinv(2,2,i,j)*grad_dGref(i,j,2,k)

                   Vscript(i,j,1,k,ie) = elem(ie)%mp(i,j)*Vscript(i,j,1,k,ie) + Vs1
                   Vscript(i,j,2,k,ie) = elem(ie)%mp(i,j)*Vscript(i,j,2,k,ie) + Vs2

                   Vscript(i,j,1,k,ie) = elem(ie)%mp(i,j)*elem(ie)%state%v(i,j,1,k,nm1) + dt2*Vscript(i,j,1,k,ie)
                   Vscript(i,j,2,k,ie) = elem(ie)%mp(i,j)*elem(ie)%state%v(i,j,2,k,nm1) + dt2*Vscript(i,j,2,k,ie)
                end do
             end do

          end do

          HrefTscript(:,:) = 0.0_real_kind

          do k=nlev,1,-1

             do j=1,np
                do i=1,np
                   hkl = rpmidref(k)*pdelref(k)
                   hkk = hkl*0.5_real_kind
                   B(i,j,k,ie)      = HrefTscript(i,j) + Rgas*hkk*Tscript(i,j,k,ie)
                   B(i,j,k,ie)      = B(i,j,k,ie) +  elem(ie)%state%phis(i,j) + RTref(k)*Pscript(i,j,ie)
                   HrefTscript(i,j) = HrefTscript(i,j) + Rgas*hkl*Tscript(i,j,k,ie)
                end do
             end do

          end do

          kptr=0
          call edgeVpack(edge2, Vscript(1,1,1,1,ie),2*nlev,kptr,elem(ie)%desc)

       end do

       call bndry_exchangeV(cg%hybrid,edge2)

       do ie = nets, nete

          kptr=0
          call edgeVunpack(edge2, Vscript(1,1,1,1,ie), 2*nlev, kptr, elem(ie)%desc)
#ifdef DEBUGOMP
#if (! defined ELEMENT_OPENMP)
!$OMP BARRIER
#endif
#endif
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k,i,j,gVscript,deriv)
#endif
          do k=1,nlev

             do j=1,np
                do i=1,np
                   Vscript(i,j,1,k,ie) = elem(ie)%rmp(i,j)*Vscript(i,j,1,k,ie)
                   Vscript(i,j,2,k,ie) = elem(ie)%rmp(i,j)*Vscript(i,j,2,k,ie)
                end do
             end do

             do j=1,np
                do i=1,np

                   ! Contravariant Vscript

                   gVscript(i,j,1) = elem(ie)%Dinv(1,1,i,j)*Vscript(i,j,1,k,ie) + &
                        elem(ie)%Dinv(1,2,i,j)*Vscript(i,j,2,k,ie)
                   gVscript(i,j,2) = elem(ie)%Dinv(2,1,i,j)*Vscript(i,j,1,k,ie) + &
                        elem(ie)%Dinv(2,2,i,j)*Vscript(i,j,2,k,ie)

                   gVscript(i,j,1) = elem(ie)%metdet(i,j)*gVscript(i,j,1)
                   gVscript(i,j,2) = elem(ie)%metdet(i,j)*gVscript(i,j,2)

                end do
             end do

             div_Vscript(:,:,k) = divergence(gVscript,deriv)*rrearth

          end do

#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k,i,j,l)
#endif
          do k=1,nlev

             do j=1,np
                do i=1,np
                   C(i,j,k,ie) = elem(ie)%metdet(i,j)*B(i,j,k,ie)
                end do
             end do

             do l=1,nlev
                do j=1,np
                   do i=1,np
                      C(i,j,k,ie) = C(i,j,k,ie) - dt*Amat(l,k)*div_Vscript(i,j,l)
                   end do
                end do
             end do

          end do

          ! ===============================================================
          !  Weight C (the RHS of the helmholtz problem) by the mass matrix
          ! ===============================================================

#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k,i,j)
#endif
          do k=1,nlev
             do j=1,np
                do i=1,np
                   C(i,j,k,ie) = elem(ie)%mp(i,j)*C(i,j,k,ie)
                end do
             end do
          end do

          ! ===================================
          ! Pack C into the edge1 buffer
          ! ===================================

          kptr=0
          call edgeVpack(edge1,C(1,1,1,ie),nlev,kptr,elem(ie)%desc)

       end do

       ! ==================================
       ! boundary exchange C
       ! ==================================

       call bndry_exchangeV(cg%hybrid,edge1)

       do ie=nets,nete

          ! ===================================
          ! Unpack C from the edge1 buffer
          ! ===================================

          kptr=0
          call edgeVunpack(edge1, C(1,1,1,ie), nlev, kptr, elem(ie)%desc)

          ! ===============================================
          ! Complete global assembly by normalizing by rmp
          ! ===============================================

#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k,i,j)
#endif
          do k=1,nlev

             do j=1,np
                do i=1,np
                   D(i,j,k,ie) = elem(ie)%rmp(i,j)*C(i,j,k,ie)
                end do
             end do

          end do

#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k,i,j,l)
#endif
          do k=1,nlev

             do j=1,np
                do i=1,np
                   C(i,j,k,ie) = 0.0_real_kind
                end do
             end do

             do l=1,nlev
                do j=1,np
                   do i=1,np
                      C(i,j,k,ie) = C(i,j,k,ie) + Emat_inv(l,k)*D(i,j,l,ie)
                   end do
                end do
             end do

          end do

       end do
#ifdef DEBUGOMP
#if (! defined ELEMENT_OPENMP)
!$OMP BARRIER
#endif
#endif
       ! ==========================================
       ! solve for Gamma_ref, given C as RHS input
       ! ==========================================

       Gamma_ref = pcg_solver(elem, &
            C,          &
            cg,         &
            red,        &
            edge1,      &
            edge2,      &
            Lambda,     &
            deriv,      &
            nets,       &
            nete,       &
            blkjac)


       ! ================================================================
       ! Backsubstitute Gamma_ref into semi-implicit system of equations
       ! to find prognostic variables at time level n+1
       ! ================================================================

       do ie = nets, nete

#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k,i,j,l)
#endif
          do k=1,nlev

             do j=1,np
                do i=1,np
                   Gref(i,j,k,ie) = 0.0_real_kind
                end do
             end do

             do l=1,nlev
                do j=1,np
                   do i=1,np
                      Gref(i,j,k,ie) = Gref(i,j,k,ie) + Emat(l,k)*Gamma_ref(i,j,l,ie)
                   end do
                end do
             end do

             do j=1,np
                do i=1,np
                   B(i,j,k,ie) = elem(ie)%mp(i,j) * dt * (B(i,j,k,ie) - Gref(i,j,k,ie))
                end do
             end do

          end do

          kptr=0
          call edgeVpack(edge1,B(:,:,:,ie),nlev,kptr,elem(ie)%desc)

       end do

       call bndry_exchangeV(cg%hybrid,edge1)

       do ie = nets, nete

          kptr=0
          call edgeVunpack(edge1, B(:,:,:,ie), nlev, kptr, elem(ie)%desc)
#ifdef DEBUGOMP
#if (! defined ELEMENT_OPENMP)
!$OMP BARRIER
#endif
#endif
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k,i,j)
#endif
          do k=1,nlev

             do j=1,np
                do i=1,np
                   B(i,j,k,ie) = elem(ie)%rmp(i,j)*B(i,j,k,ie)
                end do
             end do

          end do

#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k,i,j,l)
#endif
          do k=1,nlev

             do j=1,np
                do i=1,np
                   D(i,j,k,ie) = 0.0_real_kind
                end do
             end do

             do l=1,nlev
                do j=1,np
                   do i=1,np
                      D(i,j,k,ie) = D(i,j,k,ie) + Emat_inv(l,k)*B(i,j,l,ie)
                   end do
                end do
             end do

          end do

#if 1
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k,i,j,l)
#endif
          do k=1,nlev

             do j=1,np
                do i=1,np
                   elem(ie)%derived%div(i,j,k,np1) = 0.0_real_kind
                end do
             end do

             do l=1,nlev
                do j=1,np
                   do i=1,np
                      elem(ie)%derived%div(i,j,k,np1) = elem(ie)%derived%div(i,j,k,np1) + Emat(l,k)*D(i,j,l,ie)/Lambda(l)
                   end do
                end do
             end do

          end do
#endif

          do k=1,nlev

             grad_Gref(:,:,:,k)=gradient_wk(Gref(:,:,k,ie),deriv)*rrearth

             do j=1,np
                do i=1,np
                   gGr1 = grad_Gref(i,j,1,k)
                   gGr2 = grad_Gref(i,j,2,k)
                   elem(ie)%state%v(i,j,1,k,np1) = elem(ie)%Dinv(1,1,i,j)*gGr1 + elem(ie)%Dinv(2,1,i,j)*gGr2
                   elem(ie)%state%v(i,j,2,k,np1) = elem(ie)%Dinv(1,2,i,j)*gGr1 + elem(ie)%Dinv(2,2,i,j)*gGr2
                end do
             end do

             do j=1,np
                do i=1,np
                   Pscript(i,j,ie) = Pscript(i,j,ie) - dt*Pvec(k)*elem(ie)%derived%div(i,j,k,np1)
                end do
             end do


             do l=1,nlev
                do j=1,np
                   do i=1,np
                      Tscript(i,j,k,ie) = Tscript(i,j,k,ie) - dt*Tmat(l,k)*elem(ie)%derived%div(i,j,l,np1)
                   end do
                end do
             end do

          end do

          do j=1,np
             do i=1,np
                Pscript(i,j,ie) = elem(ie)%mp(i,j)*Pscript(i,j,ie)
             end do
          end do
          do k=1,nlev
             do j=1,np
                do i=1,np
                   Tscript(i,j,k,ie) = elem(ie)%mp(i,j)*Tscript(i,j,k,ie)
                end do
             end do
          end do

          ! ===============================================
          ! Pack v at time level n+1 into the edge3p1 buffer
          ! ===============================================

          kptr=0
          call edgeVpack(edge3p1, elem(ie)%state%v(:,:,:,:,np1),2*nlev,kptr,elem(ie)%desc)

          kptr=2*nlev
          call edgeVpack(edge3p1, Tscript(:,:,:,ie),nlev,kptr,elem(ie)%desc)

          kptr=3*nlev
          call edgeVpack(edge3p1, Pscript(:,:,ie),1,kptr,elem(ie)%desc)

       end do

       ! ======================================
       ! boundary exchange v at time level n+1
       ! ======================================

       call bndry_exchangeV(cg%hybrid,edge3p1)

       do ie=nets,nete

          ! ===================================
          ! Unpack v from the edge2 buffer
          ! ===================================

          kptr=0
          call edgeVunpack(edge3p1, elem(ie)%state%v(:,:,:,:,np1), 2*nlev, kptr, elem(ie)%desc)

          kptr=2*nlev
          call edgeVunpack(edge3p1, Tscript(:,:,:,ie), nlev, kptr, elem(ie)%desc)

          kptr=3*nlev
          call edgeVunpack(edge3p1, Pscript(:,:,ie), 1, kptr, elem(ie)%desc)

          ! ==========================================================
          ! Complete global assembly by normalizing velocity by rmp
          ! Vscript = Vscript - dt*grad(Gref)
          ! ==========================================================

#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k,i,j)
#endif
          do k=1,nlev

             do j=1,np
                do i=1,np
                   elem(ie)%state%v(i,j,1,k,np1) = Vscript(i,j,1,k,ie) + dt*elem(ie)%rmp(i,j)*elem(ie)%state%v(i,j,1,k,np1)
                   elem(ie)%state%v(i,j,2,k,np1) = Vscript(i,j,2,k,ie) + dt*elem(ie)%rmp(i,j)*elem(ie)%state%v(i,j,2,k,np1)
                end do
             end do

          end do

          do j=1,np
             do i=1,np
                elem(ie)%state%lnps(i,j,np1) = elem(ie)%rmp(i,j)*Pscript(i,j,ie)
             end do
          end do

#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k,i,j)
#endif
          do k=1,nlev
             do j=1,np
                do i=1,np
                   elem(ie)%state%T(i,j,k,np1) = elem(ie)%rmp(i,j)*Tscript(i,j,k,ie)
                end do
             end do
          end do

       end do
#ifdef DEBUGOMP
#if (! defined ELEMENT_OPENMP)
!$OMP BARRIER
#endif
#endif
#if 1
       call prim_diffusion(elem, nets,nete,np1,deriv,dt2,cg%hybrid)
#endif
       call t_stopf('prim_advance_si')
  end subroutine prim_advance_si


  subroutine preq_robert3(nm1,n0,np1,elem,hvcoord,nets,nete)
  use dimensions_mod, only : np, nlev, qsize
  use hybvcoord_mod, only : hvcoord_t
  use element_mod, only : element_t
  use time_mod, only: smooth
  use control_mod, only : integration

  implicit none
  integer              , intent(in) :: nm1,n0,np1,nets,nete
  type (hvcoord_t), intent(in)      :: hvcoord
  type (element_t)     , intent(inout) :: elem(:)


  integer :: i,j,k,ie,q
  real (kind=real_kind) :: dp
  logical :: filter_ps = .false.
  if (integration == "explicit") filter_ps = .true.

  call t_startf('preq_robert')
  do ie=nets,nete
     if (filter_ps) then
        elem(ie)%state%ps_v(:,:,n0) = elem(ie)%state%ps_v(:,:,n0) + smooth*(elem(ie)%state%ps_v(:,:,nm1) &
             - 2.0D0*elem(ie)%state%ps_v(:,:,n0)   + elem(ie)%state%ps_v(:,:,np1))
        elem(ie)%state%lnps(:,:,n0) = LOG(elem(ie)%state%ps_v(:,:,n0))
     else
        elem(ie)%state%lnps(:,:,n0) = elem(ie)%state%lnps(:,:,n0) + smooth*(elem(ie)%state%lnps(:,:,nm1) &
             - 2.0D0*elem(ie)%state%lnps(:,:,n0)   + elem(ie)%state%lnps(:,:,np1))
        elem(ie)%state%ps_v(:,:,n0) = EXP(elem(ie)%state%lnps(:,:,n0))
     endif

     elem(ie)%state%T(:,:,:,n0) = elem(ie)%state%T(:,:,:,n0) + smooth*(elem(ie)%state%T(:,:,:,nm1) &
          - 2.0D0*elem(ie)%state%T(:,:,:,n0)   + elem(ie)%state%T(:,:,:,np1))
     elem(ie)%state%v(:,:,:,:,n0) = elem(ie)%state%v(:,:,:,:,n0) + smooth*(elem(ie)%state%v(:,:,:,:,nm1) &
          - 2.0D0*elem(ie)%state%v(:,:,:,:,n0) + elem(ie)%state%v(:,:,:,:,np1))

  end do
  call t_stopf('preq_robert')

  end subroutine preq_robert3




  subroutine applyCAMforcing(elem,hvcoord,np1,np1_qdp,dt_q,nets,nete)
  use dimensions_mod, only : np, nlev, qsize
  use element_mod, only : element_t
  use hybvcoord_mod, only : hvcoord_t
  use control_mod, only : moisture

  implicit none
  type (element_t)     , intent(inout) :: elem(:)
  real (kind=real_kind), intent(in) :: dt_q
  type (hvcoord_t), intent(in)      :: hvcoord
  integer,  intent(in) :: np1,nets,nete,np1_qdp

  ! local
  integer :: i,j,k,ie,q
  real (kind=real_kind) :: v1,dp
  logical :: wet

  wet = (moisture /= "dry")

  do ie=nets,nete
     ! apply forcing to Qdp
     elem(ie)%derived%FQps(:,:,1)=0
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(q,k,i,j,v1)
#endif
     do q=1,qsize
        do k=1,nlev
           do j=1,np
              do i=1,np
                 v1 = dt_q*elem(ie)%derived%FQ(i,j,k,q,1)
                 !if (elem(ie)%state%Qdp(i,j,k,q,np1) + v1 < 0 .and. v1<0) then
                 if (elem(ie)%state%Qdp(i,j,k,q,np1_qdp) + v1 < 0 .and. v1<0) then
                    !if (elem(ie)%state%Qdp(i,j,k,q,np1) < 0 ) then
                    if (elem(ie)%state%Qdp(i,j,k,q,np1_qdp) < 0 ) then
                       v1=0  ! Q already negative, dont make it more so
                    else
                       !v1 = -elem(ie)%state%Qdp(i,j,k,q,np1)
                       v1 = -elem(ie)%state%Qdp(i,j,k,q,np1_qdp)
                    endif
                 endif
                 !elem(ie)%state%Qdp(i,j,k,q,np1) = elem(ie)%state%Qdp(i,j,k,q,np1)+v1
                 elem(ie)%state%Qdp(i,j,k,q,np1_qdp) = elem(ie)%state%Qdp(i,j,k,q,np1_qdp)+v1
                 if (q==1) then
                    elem(ie)%derived%FQps(i,j,1)=elem(ie)%derived%FQps(i,j,1)+v1/dt_q
                 endif
              enddo
           enddo
        enddo
     enddo

     if (wet .and. qsize>0) then
        ! to conserve dry mass in the precese of Q1 forcing:
        elem(ie)%state%ps_v(:,:,np1) = elem(ie)%state%ps_v(:,:,np1) + &
             dt_q*elem(ie)%derived%FQps(:,:,1)
     endif

     ! Qdp(np1) and ps_v(np1) were updated by forcing - update Q(np1)
     ! ps_v(n0) may also have been changed if using Robert,
     ! but Q(n0) will be updated after robert filter
     ! so no need to do that now
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(q,k,i,j,dp)
#endif
     do q=1,qsize
        do k=1,nlev
           do j=1,np
              do i=1,np
                 dp = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
                      ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(i,j,np1)
                 elem(ie)%state%Q(i,j,k,q) = elem(ie)%state%Qdp(i,j,k,q,np1_qdp)/dp
              enddo
           enddo
        enddo
     enddo

     elem(ie)%state%T(:,:,:,np1) = elem(ie)%state%T(:,:,:,np1) + &
          dt_q*elem(ie)%derived%FT(:,:,:,1)
     elem(ie)%state%v(:,:,:,:,np1) = elem(ie)%state%v(:,:,:,:,np1) + &
          dt_q*elem(ie)%derived%FM(:,:,:,:,1)
  enddo
  end subroutine applyCAMforcing



  subroutine applyCAMforcing_dynamics(elem,hvcoord,np1,dt_q,nets,nete)
  use dimensions_mod, only : np, nlev, qsize
  use element_mod, only : element_t
  use hybvcoord_mod, only : hvcoord_t

  implicit none
  type (element_t)     , intent(inout) :: elem(:)
  real (kind=real_kind), intent(in) :: dt_q
  type (hvcoord_t), intent(in)      :: hvcoord
  integer,  intent(in) :: np1,nets,nete

  ! local
  integer :: i,j,k,ie,q
  real (kind=real_kind) :: v1,dp
  logical :: wet

  do ie=nets,nete
     elem(ie)%state%T(:,:,:,np1) = elem(ie)%state%T(:,:,:,np1) + &
          dt_q*elem(ie)%derived%FT(:,:,:,1)
     elem(ie)%state%v(:,:,:,:,np1) = elem(ie)%state%v(:,:,:,:,np1) + &
          dt_q*elem(ie)%derived%FM(:,:,:,:,1)
  enddo
  end subroutine applyCAMforcing_dynamics



  subroutine advance_hypervis(edge3,elem,hvcoord,hybrid,deriv,nt,nets,nete,dt2,eta_ave_w)
  !
  !  take one timestep of:
  !          u(:,:,:,np) = u(:,:,:,np) +  dt2*nu*laplacian**order ( u )
  !          T(:,:,:,np) = T(:,:,:,np) +  dt2*nu_s*laplacian**order ( T )
  !
  !
  !  For correct scaling, dt2 should be the same 'dt2' used in the leapfrog advace
  !
  !
  use dimensions_mod, only : np, np, nlev
  use control_mod, only : nu, nu_div, nu_s, hypervis_order, hypervis_subcycle, nu_p, nu_top, psurf_vis
  use hybrid_mod, only : hybrid_t
  use hybvcoord_mod, only : hvcoord_t
  use element_mod, only : element_t
  use derivative_mod, only : derivative_t, laplace_sphere_wk, vlaplace_sphere_wk
  use edge_mod, only : EdgeBuffer_t, edgevpack, edgevunpack
  use bndry_mod, only : bndry_exchangev
  use viscosity_mod, only : biharmonic_wk
  use physical_constants, only: Cp
!  use time_mod, only : TimeLevel_t
  implicit none

  type (hybrid_t)      , intent(in) :: hybrid
  type (element_t)     , intent(inout), target :: elem(:)
  type (EdgeBuffer_t)  , intent(inout) :: edge3
  type (derivative_t)  , intent(in) :: deriv
  type (hvcoord_t), intent(in)      :: hvcoord
!  type (TimeLevel_t)   , intent(in) :: tl

  real (kind=real_kind) :: dt2
  integer :: nets,nete

  ! local
  real (kind=real_kind) :: eta_ave_w  ! weighting for mean flux terms
  real (kind=real_kind) :: nu_scale, dpdn,dpdn0, nu_scale_top,nu_ratio
  integer :: k,kptr,i,j,ie,ic,nt
  real (kind=real_kind), dimension(np,np,2,nlev,nets:nete)      :: vtens
  real (kind=real_kind), dimension(np,np,nlev,nets:nete)        :: ptens
  real (kind=real_kind), dimension(np,np,nets:nete) :: pstens
  real (kind=real_kind), dimension(np,np,nlev) :: p
  real (kind=real_kind), dimension(np,np) :: dptemp1,dptemp2


! NOTE: PGI compiler bug: when using spheremp, rspheremp and ps as pointers to elem(ie)% members,
  !       data is incorrect (offset by a few numbers actually)
  !       removed for now.
  !       real (kind=real_kind), dimension(:,:), pointer :: spheremp,rspheremp
  !       real (kind=real_kind), dimension(:,:,:), pointer   :: ps

  real (kind=real_kind), dimension(np,np) :: lap_p
  real (kind=real_kind), dimension(np,np,2) :: lap_v
  real (kind=real_kind) :: v1,v2,dt,heating,utens_tmp,vtens_tmp,ptens_tmp


  if (nu_s == 0 .and. nu == 0 .and. nu_p==0 ) return;
  call t_barrierf('sync_advance_hypervis', hybrid%par%comm)
  call t_startf('advance_hypervis')


  dt=dt2/hypervis_subcycle
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  regular viscosity
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (hypervis_order == 1) then
     if (nu_p>0) call abortmp( 'ERROR: hypervis_order == 1 not coded for nu_p>0')
     do ic=1,hypervis_subcycle
        do ie=nets,nete

#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k,lap_p,lap_v,deriv,i,j)
#endif
           do k=1,nlev
              lap_p=laplace_sphere_wk(elem(ie)%state%T(:,:,k,nt),deriv,elem(ie),var_coef=.false.)
              lap_v=vlaplace_sphere_wk(elem(ie)%state%v(:,:,:,k,nt),deriv,elem(ie),var_coef=.false.)
              ! advace in time.  (note: DSS commutes with time stepping, so we
              ! can time advance and then DSS.  this has the advantage of
              ! not letting any discontinuties accumulate in p,v via roundoff
              do j=1,np
                 do i=1,np
                    elem(ie)%state%T(i,j,k,nt)=elem(ie)%state%T(i,j,k,nt)*elem(ie)%spheremp(i,j)  +  dt*nu_s*lap_p(i,j)
                    elem(ie)%state%v(i,j,1,k,nt)=elem(ie)%state%v(i,j,1,k,nt)*elem(ie)%spheremp(i,j) + dt*nu*lap_v(i,j,1)
                    elem(ie)%state%v(i,j,2,k,nt)=elem(ie)%state%v(i,j,2,k,nt)*elem(ie)%spheremp(i,j) + dt*nu*lap_v(i,j,2)
                 enddo
              enddo
           enddo

           kptr=0
           call edgeVpack(edge3, elem(ie)%state%T(:,:,:,nt),nlev,kptr,elem(ie)%desc)
           kptr=nlev
           call edgeVpack(edge3,elem(ie)%state%v(:,:,:,:,nt),2*nlev,kptr,elem(ie)%desc)
        enddo

        call bndry_exchangeV(hybrid,edge3)

        do ie=nets,nete

           kptr=0
           call edgeVunpack(edge3, elem(ie)%state%T(:,:,:,nt), nlev, kptr, elem(ie)%desc)
           kptr=nlev
           call edgeVunpack(edge3, elem(ie)%state%v(:,:,:,:,nt), 2*nlev, kptr, elem(ie)%desc)

           ! apply inverse mass matrix
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k,i,j)
#endif
           do k=1,nlev
              do j=1,np
                 do i=1,np
                    elem(ie)%state%T(i,j,k,nt)=elem(ie)%rspheremp(i,j)*elem(ie)%state%T(i,j,k,nt)
                    elem(ie)%state%v(i,j,1,k,nt)=elem(ie)%rspheremp(i,j)*elem(ie)%state%v(i,j,1,k,nt)
                    elem(ie)%state%v(i,j,2,k,nt)=elem(ie)%rspheremp(i,j)*elem(ie)%state%v(i,j,2,k,nt)
                 enddo
              enddo
           enddo
        enddo
#ifdef DEBUGOMP
#if (! defined ELEMENT_OPENMP)
!$OMP BARRIER
#endif
#endif
     enddo  ! subcycle
  endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  hyper viscosity
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! nu_p=0:
!   scale T dissipaton by dp  (conserve IE, dissipate T^2)
! nu_p>0
!   dont scale:  T equation IE dissipation matches (to truncation error)
!                IE dissipation from continuity equation
!                (1 deg: to about 0.1 W/m^2)
!
  if (hypervis_order == 2) then
     nu_ratio = nu_div/nu ! possibly weight div component more inside biharmonc_wk
     do ic=1,hypervis_subcycle
        call biharmonic_wk(elem,pstens,ptens,vtens,deriv,edge3,hybrid,nt,nets,nete,nu_ratio)
        do ie=nets,nete

           ! comptue mean flux
           if (nu_p>0) then
#if 0
              elem(ie)%derived%psdiss_ave(:,:)=&
                   elem(ie)%derived%psdiss_ave(:,:)+eta_ave_w*elem(ie)%state%ps_v(:,:,nt)/hypervis_subcycle
              elem(ie)%derived%psdiss_biharmonic(:,:)=&
                   elem(ie)%derived%psdiss_biharmonic(:,:)+eta_ave_w*pstens(:,:,ie)/hypervis_subcycle
#else
              do k=1,nlev
                 dptemp1(:,:) = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
                      ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(:,:,nt)
                 elem(ie)%derived%dpdiss_ave(:,:,k)=elem(ie)%derived%dpdiss_ave(:,:,k)+eta_ave_w*dptemp1(:,:)/hypervis_subcycle

                 dptemp2(:,:) = (hvcoord%hybi(k+1)-hvcoord%hybi(k))*pstens(:,:,ie)
                 elem(ie)%derived%dpdiss_biharmonic(:,:,k)=&
                      elem(ie)%derived%dpdiss_biharmonic(:,:,k)+eta_ave_w*dptemp2(:,:)/hypervis_subcycle
              enddo
#endif
           endif
           nu_scale=1
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k,i,j,lap_p,lap_v,nu_scale_top,dpdn,dpdn0,nu_scale,utens_tmp,vtens_tmp,ptens_tmp)
#endif
           do k=1,nlev
              ! advace in time.
              ! note: DSS commutes with time stepping, so we can time advance and then DSS.
              ! note: weak operators alreayd have mass matrix "included"

              ! add regular diffusion in top 3 layers:
              if (nu_top>0 .and. k<=3) then
                 lap_p=laplace_sphere_wk(elem(ie)%state%T(:,:,k,nt),deriv,elem(ie),var_coef=.false.)
                 lap_v=vlaplace_sphere_wk(elem(ie)%state%v(:,:,:,k,nt),deriv,elem(ie),var_coef=.false.)
              endif
              nu_scale_top = 1
              if (k==1) nu_scale_top=4
              if (k==2) nu_scale_top=2

              do j=1,np
                 do i=1,np
                    if (nu_p==0) then
                       ! normalize so as to conserve IE
                       ! scale by 1/rho (normalized to be O(1))
                       ! dp/dn = O(ps0)*O(delta_eta) = O(ps0)/O(nlev)
                       dpdn = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
                            ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(i,j,nt)
                       dpdn0 = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
                            ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*hvcoord%ps0
                       nu_scale = dpdn0/dpdn
                    endif

                    ! biharmonic terms need a negative sign:
                    if (nu_top>0 .and. k<=3) then
                       utens_tmp=(-nu*vtens(i,j,1,k,ie) + nu_scale_top*nu_top*lap_v(i,j,1))
                       vtens_tmp=(-nu*vtens(i,j,2,k,ie) + nu_scale_top*nu_top*lap_v(i,j,2))
                       ptens_tmp=nu_scale*(-nu_s*ptens(i,j,k,ie) + nu_scale_top*nu_top*lap_p(i,j) )
                    else
                       utens_tmp=-nu*vtens(i,j,1,k,ie)
                       vtens_tmp=-nu*vtens(i,j,2,k,ie)
                       ptens_tmp=-nu_scale*nu_s*ptens(i,j,k,ie)
                    endif

                    ptens(i,j,k,ie) = ptens_tmp
                    vtens(i,j,1,k,ie)=utens_tmp
                    vtens(i,j,2,k,ie)=vtens_tmp
                 enddo
              enddo
           enddo

           pstens(:,:,ie)  =  -nu_p*pstens(:,:,ie)
           kptr=0
           call edgeVpack(edge3, ptens(:,:,:,ie),nlev,kptr,elem(ie)%desc)
           kptr=nlev
           call edgeVpack(edge3,vtens(:,:,:,:,ie),2*nlev,kptr,elem(ie)%desc)
           kptr=3*nlev
           call edgeVpack(edge3,pstens(:,:,ie),1,kptr,elem(ie)%desc)
        enddo


        call bndry_exchangeV(hybrid,edge3)

        do ie=nets,nete

           kptr=0
           call edgeVunpack(edge3, ptens(:,:,:,ie), nlev, kptr, elem(ie)%desc)
           kptr=nlev
           call edgeVunpack(edge3, vtens(:,:,:,:,ie), 2*nlev, kptr, elem(ie)%desc)


           ! apply inverse mass matrix, accumulate tendencies
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k,i,j,v1,v2,heating)
#endif
           do k=1,nlev
              vtens(:,:,1,k,ie)=dt*vtens(:,:,1,k,ie)*elem(ie)%rspheremp(:,:)
              vtens(:,:,2,k,ie)=dt*vtens(:,:,2,k,ie)*elem(ie)%rspheremp(:,:)
              ptens(:,:,k,ie)=dt*ptens(:,:,k,ie)*elem(ie)%rspheremp(:,:)
           enddo

           ! apply hypervis to u -> u+utens:
           ! E0 = dpdn * .5*u dot u + dpdn * T  + dpdn*PHIS
           ! E1 = dpdn * .5*(u+utens) dot (u+utens) + dpdn * (T-X) + dpdn*PHIS
           ! E1-E0:   dpdn (u dot utens) + dpdn .5 utens dot utens   - dpdn X
           !      X = (u dot utens) + .5 utens dot utens
           !  alt:  (u+utens) dot utens
           do k=1,nlev
              do j=1,np
                 do i=1,np
                    ! update v first (gives better results than updating v after heating)
                    elem(ie)%state%v(i,j,:,k,nt)=elem(ie)%state%v(i,j,:,k,nt) + &
                         vtens(i,j,:,k,ie)

                    v1=elem(ie)%state%v(i,j,1,k,nt)
                    v2=elem(ie)%state%v(i,j,2,k,nt)
                    heating = (vtens(i,j,1,k,ie)*v1  + vtens(i,j,2,k,ie)*v2 )
                    elem(ie)%state%T(i,j,k,nt)=elem(ie)%state%T(i,j,k,nt) &
                         +ptens(i,j,k,ie)-heating/cp

                 enddo
              enddo
           enddo

           if (nu_p>0) then
              kptr=3*nlev
              call edgeVunpack(edge3, pstens(:,:,ie), 1, kptr, elem(ie)%desc)
              pstens(:,:,ie)=dt*pstens(:,:,ie)*elem(ie)%rspheremp(:,:)
              elem(ie)%state%ps_v(:,:,nt)=elem(ie)%state%ps_v(:,:,nt) + pstens(:,:,ie)
           endif

        enddo
#ifdef DEBUGOMP
#if (! defined ELEMENT_OPENMP)
!$OMP BARRIER
#endif
#endif
     enddo
  endif

  call t_stopf('advance_hypervis')

  end subroutine advance_hypervis




  subroutine advance_hypervis_dp(edge3,elem,hvcoord,hybrid,deriv,nt,nets,nete,dt2,eta_ave_w)
  !
  !  take one timestep of:
  !          u(:,:,:,np) = u(:,:,:,np) +  dt2*nu*laplacian**order ( u )
  !          T(:,:,:,np) = T(:,:,:,np) +  dt2*nu_s*laplacian**order ( T )
  !
  !
  !  For correct scaling, dt2 should be the same 'dt2' used in the leapfrog advace
  !
  !
  use dimensions_mod, only : np, np, nlev
  use control_mod, only : nu, nu_div, nu_s, hypervis_order, hypervis_subcycle, nu_p, nu_top, psurf_vis
  use hybrid_mod, only : hybrid_t
  use hybvcoord_mod, only : hvcoord_t
  use element_mod, only : element_t
  use derivative_mod, only : derivative_t, laplace_sphere_wk, vlaplace_sphere_wk
  use edge_mod, only : EdgeBuffer_t, edgevpack, edgevunpack
  use bndry_mod, only : bndry_exchangev
  use viscosity_mod, only : biharmonic_wk_dp3d
  use physical_constants, only: Cp, rrearth
!  use time_mod, only : TimeLevel_t
  implicit none

  type (hybrid_t)      , intent(in) :: hybrid
  type (element_t)     , intent(inout), target :: elem(:)
  type (EdgeBuffer_t)  , intent(inout) :: edge3
  type (derivative_t)  , intent(in) :: deriv
  type (hvcoord_t), intent(in)      :: hvcoord
!  type (TimeLevel_t)   , intent(in) :: tl

  real (kind=real_kind) :: dt2
  integer :: nets,nete

  ! local
  real (kind=real_kind) :: eta_ave_w  ! weighting for mean flux terms
  real (kind=real_kind) :: dpdn,dpdn0, nu_scale_top,nu_ratio
  integer :: k,kptr,i,j,ie,ic,nt
  real (kind=real_kind), dimension(np,np,2,nlev,nets:nete)      :: vtens
  real (kind=real_kind), dimension(np,np,nlev,nets:nete)        :: ttens
  real (kind=real_kind), dimension(np,np,nlev,nets:nete)        :: dptens
  real (kind=real_kind), dimension(np,np,nlev) :: p
  real (kind=real_kind), dimension(np,np) :: dptemp1,dptemp2


! NOTE: PGI compiler bug: when using spheremp, rspheremp and ps as pointers to elem(ie)% members,
  !       data is incorrect (offset by a few numbers actually)
  !       removed for now.
  !       real (kind=real_kind), dimension(:,:), pointer :: spheremp,rspheremp
  !       real (kind=real_kind), dimension(:,:,:), pointer   :: ps

  real (kind=real_kind), dimension(np,np) :: lap_t,lap_dp
  real (kind=real_kind), dimension(np,np,2) :: lap_v
  real (kind=real_kind) :: v1,v2,dt,heating,utens_tmp,vtens_tmp,ttens_tmp,dptens_tmp

  real (kind=real_kind), dimension(np,np,nlev) :: elem_derived_dpdiss_ave
  pointer(elem_derived_dpdiss_ave_ptr, elem_derived_dpdiss_ave)

  real (kind=real_kind), dimension(np,np,nlev) :: elem_state_dp3d
  pointer(elem_state_dp3d_ptr, elem_state_dp3d)
  
  real (kind=real_kind), dimension(np,np,nlev) :: elem_derived_dpdiss_biharmonic
  pointer(elem_derived_dpdiss_biharmonic_ptr, elem_derived_dpdiss_biharmonic)

  real (kind=real_kind), dimension(np,np,nlev) :: elem_state_T
  pointer(elem_state_T_ptr, elem_state_T)

  real (kind=real_kind), dimension(np,np,2,nlev) :: elem_state_v
  pointer(elem_state_v_ptr, elem_state_v)
 
  real (kind=real_kind), dimension(np,np) :: elem_metdet
  pointer(elem_metdet_ptr, elem_metdet)
 
  real (kind=real_kind), dimension(np,np) :: elem_rmetdet
  pointer(elem_rmetdet_ptr, elem_rmetdet)

  real (kind=real_kind), dimension(2,2,np,np) :: elem_Dinv
  pointer(elem_Dinv_ptr, elem_Dinv)

  real (kind=real_kind), dimension(2,2,np,np) :: elem_D
  pointer(elem_D_ptr, elem_D)

  real (kind=real_kind), dimension(2,2,np,np) :: elem_metinv
  pointer(elem_metinv_ptr, elem_metinv)

  real (kind=real_kind), dimension(np,np) :: deriv_dvv
  pointer(deriv_dvv_ptr, deriv_dvv)

  real (kind=real_kind), dimension(np,np) :: elem_mp
  pointer(elem_mp_ptr, elem_mp)

  real (kind=real_kind), dimension(np,np) :: elem_spheremp
  pointer(elem_spheremp_ptr, elem_spheremp)

  integer(kind=8), dimension(12,nets:nete) :: elem_array

  if (nu_s == 0 .and. nu == 0 .and. nu_p==0 ) return;
  call t_barrierf('sync_advance_hypervis', hybrid%par%comm)
  call t_startf('advance_hypervis_dp')


  dt=dt2/hypervis_subcycle


  do ie = nets,nete
      elem_array(1,ie) = loc(elem(ie)%derived%dpdiss_ave)
      elem_array(2,ie) = loc(elem(ie)%state%dp3d(:,:,:,nt))
      elem_array(3,ie) = loc(elem(ie)%derived%dpdiss_biharmonic)
      elem_array(4,ie) = loc(elem(ie)%state%T(:,:,:,nt))
      elem_array(5,ie) = loc(elem(ie)%state%v(:,:,:,:,nt))
      elem_array(6,ie) = loc(elem(ie)%metdet)
      elem_array(7,ie) = loc(elem(ie)%rmetdet)
      elem_array(8,ie) = loc(elem(ie)%D)
      elem_array(9,ie) = loc(elem(ie)%Dinv)
      elem_array(10,ie)= loc(elem(ie)%mp)
      elem_array(11,ie)= loc(elem(ie)%metinv)
      elem_array(12,ie)= loc(elem(ie)%spheremp)
  enddo
  deriv_dvv_ptr = loc(deriv%dvv)
     nu_ratio = nu_div/nu ! possibly weight div component more inside biharmonc_wk
     do ic=1,hypervis_subcycle
        call biharmonic_wk_dp3d(elem,dptens,ttens,vtens,deriv,edge3,hybrid,nt,nets,nete,nu_ratio)

        ! conghui: this may be optimized
        call t_startf("hypervis_dp before bndry")
        do ie=nets,nete
           do k=1,nlev
              elem_derived_dpdiss_ave_ptr        = elem_array(1,ie)
              elem_state_dp3d_ptr                = elem_array(2,ie)
              elem_derived_dpdiss_biharmonic_ptr = elem_array(3,ie)
              elem_state_T_ptr                   = elem_array(4,ie)
              elem_state_v_ptr                   = elem_array(5,ie)
              elem_metdet_ptr                    = elem_array(6,ie)
              elem_rmetdet_ptr                   = elem_array(7,ie)
              elem_D_ptr                         = elem_array(8,ie)
              elem_Dinv_ptr                      = elem_array(9,ie)
              elem_mp_ptr                        = elem_array(10,ie)  
              elem_metinv_ptr                    = elem_array(11,ie)
              elem_spheremp_ptr                  = elem_array(12,ie)

              elem_derived_dpdiss_ave(:,:,k)=elem_derived_dpdiss_ave(:,:,k)+&
                   eta_ave_w*elem_state_dp3d(:,:,k)/hypervis_subcycle
              elem_derived_dpdiss_biharmonic(:,:,k)=elem_derived_dpdiss_biharmonic(:,:,k)+&
                   eta_ave_w*dptens(:,:,k,ie)/hypervis_subcycle
              if (nu_top>0 .and. k<=3) then
                 lap_t=my_laplace_sphere_wk(elem_state_T(:,:,k),deriv,elem(ie))
                 lap_dp=my_laplace_sphere_wk2(elem_state_dp3d(:,:,k),deriv,elem(ie))
                 lap_v=my_laplace_sphere_wk_new(elem_state_v(:,:,:,k),elem(ie),rrearth, elem_metdet, &
                 elem_rmetdet, elem_D, elem_Dinv, deriv_dvv, elem_mp,elem_metinv,elem_spheremp)
              endif
              nu_scale_top = 1
              if (k==1) nu_scale_top=4
              if (k==2) nu_scale_top=2

              do j=1,np
                 do i=1,np
                    ! biharmonic terms need a negative sign:
                    if (nu_top>0 .and. k<=3) then
                       utens_tmp=(-nu*vtens(i,j,1,k,ie) + nu_scale_top*nu_top*lap_v(i,j,1))
                       vtens_tmp=(-nu*vtens(i,j,2,k,ie) + nu_scale_top*nu_top*lap_v(i,j,2))
                       ttens_tmp=(-nu_s*ttens(i,j,k,ie) + nu_scale_top*nu_top*lap_t(i,j) )
                       dptens_tmp=(-nu_p*dptens(i,j,k,ie) + nu_scale_top*nu_top*lap_dp(i,j) )
                    else
                       utens_tmp=-nu*vtens(i,j,1,k,ie)
                       vtens_tmp=-nu*vtens(i,j,2,k,ie)
                       ttens_tmp=-nu_s*ttens(i,j,k,ie)
                       dptens_tmp=-nu_p*dptens(i,j,k,ie)
                    endif
                    ttens(i,j,k,ie) = ttens_tmp
                    dptens(i,j,k,ie) =dptens_tmp
                    vtens(i,j,1,k,ie)=utens_tmp
                    vtens(i,j,2,k,ie)=vtens_tmp
                 enddo
              enddo
           enddo
        enddo
        do ie=nets,nete

           kptr=0
           call edgeVpack(edge3, ttens(:,:,:,ie),nlev,kptr,elem(ie)%desc)
           kptr=nlev
           call edgeVpack(edge3,vtens(:,:,:,:,ie),2*nlev,kptr,elem(ie)%desc)
           kptr=3*nlev
           call edgeVpack(edge3,dptens(:,:,:,ie),nlev,kptr,elem(ie)%desc)
        enddo

        call t_stopf("hypervis_dp before bndry")

        call bndry_exchangeV(hybrid,edge3)

        call t_startf("hypervis_dp after bndry")
        do ie=nets,nete

           kptr=0
           call edgeVunpack(edge3, ttens(:,:,:,ie), nlev, kptr, elem(ie)%desc)
           kptr=nlev
           call edgeVunpack(edge3, vtens(:,:,:,:,ie), 2*nlev, kptr, elem(ie)%desc)
           kptr=3*nlev
           call edgeVunpack(edge3, dptens(:,:,:,ie), nlev, kptr, elem(ie)%desc)
        !enddo
        !do ie=nets,nete
           do k=1,nlev
              vtens(:,:,1,k,ie)=dt*vtens(:,:,1,k,ie)*elem(ie)%rspheremp(:,:)
              vtens(:,:,2,k,ie)=dt*vtens(:,:,2,k,ie)*elem(ie)%rspheremp(:,:)
              ttens(:,:,k,ie)=dt*ttens(:,:,k,ie)*elem(ie)%rspheremp(:,:)
              dptens(:,:,k,ie)=dt*dptens(:,:,k,ie)*elem(ie)%rspheremp(:,:)
           enddo

           ! apply hypervis to u -> u+utens:
           ! E0 = dpdn * .5*u dot u + dpdn * T  + dpdn*PHIS
           ! E1 = dpdn * .5*(u+utens) dot (u+utens) + dpdn * (T-X) + dpdn*PHIS
           ! E1-E0:   dpdn (u dot utens) + dpdn .5 utens dot utens   - dpdn X
           !      X = (u dot utens) + .5 utens dot utens
           !  alt:  (u+utens) dot utens
           do k=1,nlev
              do j=1,np
                 do i=1,np
                    ! update v first (gives better results than updating v after heating)
                    elem(ie)%state%v(i,j,:,k,nt)=elem(ie)%state%v(i,j,:,k,nt) + &
                         vtens(i,j,:,k,ie)

                    v1=elem(ie)%state%v(i,j,1,k,nt)
                    v2=elem(ie)%state%v(i,j,2,k,nt)
                    heating = (vtens(i,j,1,k,ie)*v1  + vtens(i,j,2,k,ie)*v2 )
                    elem(ie)%state%T(i,j,k,nt)=elem(ie)%state%T(i,j,k,nt) &
                         +ttens(i,j,k,ie)-heating/cp
                    elem(ie)%state%dp3d(i,j,k,nt)=elem(ie)%state%dp3d(i,j,k,nt) + &
                         dptens(i,j,k,ie)
                 enddo
              enddo
           enddo
         enddo ! ie
         call t_stopf("hypervis_dp after bndry")
       enddo ! cycle

  call t_stopf('advance_hypervis_dp')

  end subroutine advance_hypervis_dp

  function my_laplace_sphere_wk2(s,deriv,elem) result(laplace)

!
!   input:  s = scalar
!   ouput:  -< grad(PHI), grad(s) >   = weak divergence of grad(s)
!     note: for this form of the operator, grad(s) does not need to be made C0
!            
    use derivative_mod, only : derivative_t, divergence_sphere_wk, gradient_sphere
    use element_mod, only: element_t
    real(kind=real_kind), intent(in) :: s(4,4) 
    type (derivative_t),intent(in)              :: deriv

    type (element_t),intent(in)                 :: elem
    real(kind=real_kind)             :: laplace(4,4)
    integer i,j

    ! Local
    real(kind=real_kind) :: grads(4,4,2), oldgrads(4,4,2)

    grads=gradient_sphere(s,deriv,elem%Dinv)
 
    laplace=divergence_sphere_wk(grads,deriv,elem)

  end function my_laplace_sphere_wk2

  function my_laplace_sphere_wk(s,deriv,elem) result(laplace)

!
!   input:  s = scalar
!   ouput:  -< grad(PHI), grad(s) >   = weak divergence of grad(s)
!     note: for this form of the operator, grad(s) does not need to be made C0
!            
    use derivative_mod, only : derivative_t, divergence_sphere_wk, gradient_sphere
    use element_mod, only: element_t
    real(kind=real_kind), intent(in) :: s(4,4) 
    type (derivative_t),intent(in)              :: deriv

    type (element_t),intent(in)                 :: elem
    real(kind=real_kind)             :: laplace(4,4)
    integer i,j

    ! Local
    real(kind=real_kind) :: grads(4,4,2), oldgrads(4,4,2)

    grads=gradient_sphere(s,deriv,elem%Dinv)
 
    laplace=divergence_sphere_wk(grads,deriv,elem)

  end function my_laplace_sphere_wk

  function my_laplace_sphere_wk_new(v,elem,rrearth, elem_metdet,  elem_rmetdet, elem_D, elem_Dinv, deriv_dvv, elem_mp, elem_metinv, elem_spheremp) result(laplace)
!
!   input:  v = vector in lat-lon coordinates
!   ouput:  weak laplacian of v, in lat-lon coordinates
!
!   du/dt = laplace(u) = grad(div) - curl(vor)
!   < PHI du/dt > = < PHI laplace(u) >        PHI = covariant, u = contravariant
!                 = < PHI grad(div) >  - < PHI curl(vor) >
!                 = grad_wk(div) - curl_wk(vor)               
    use derivative_mod, only : derivative_t, divergence_sphere, gradient_sphere, vorticity_sphere, gradient_sphere_wk_testcov, curl_sphere_wk_testcov
    use element_mod, only: element_t
    real(kind=real_kind), intent(in) :: v(4,4,2) 
    !type (derivative_t)              :: deriv
    type (element_t)                 :: elem
    real(kind=real_kind) :: laplace(4,4,2)
    real(kind=real_kind) :: rrearth
    real(kind=real_kind), dimension(4,4) :: elem_metdet, elem_rmetdet,deriv_dvv, elem_mp, elem_spheremp
    real(kind=real_kind), dimension(2,2,4,4) :: elem_D, elem_Dinv, elem_metinv
    ! Local

    integer i,j,l,m,n
    real(kind=real_kind) :: vor(4,4),div(4,4)
    real(kind=real_kind) :: v1,v2,div1,div2,vor1,vor2,phi_x,phi_y


     call my_divergence_sphere(v,deriv_Dvv(:,:),elem_metdet(:,:),elem_rmetdet(:,:),elem_Dinv(:,:,:,:),rrearth, div)
     call my_vorticity_sphere(v,deriv_Dvv(:,:),elem_D(:,:,:,:), elem_rmetdet(:,:),rrearth,vor)

    laplace = my_gradient_sphere_wk_testcov(div,deriv_dvv,elem, rrearth,elem_metdet, elem_D, elem_mp, elem_metinv) - &
         my_curl_sphere_wk_testcov(vor,deriv_dvv, rrearth, elem_mp, elem_D)

    do n=1, 4
       do m=1, 4
          ! add in correction so we dont damp rigid rotation
          laplace(m,n,1)=laplace(m,n,1) + 2*elem_spheremp(m,n)*v(m,n,1)*(rrearth**2)
          laplace(m,n,2)=laplace(m,n,2) + 2*elem_spheremp(m,n)*v(m,n,2)*(rrearth**2)
       enddo
    enddo
  end function my_laplace_sphere_wk_new

  
  function my_curl_sphere_wk_testcov(s,dvv,rrearth, elem_mp, elem_D) result(ds)
    !type (derivative_t)              :: deriv
    real(kind=real_kind)             :: dvv(4,4)
    real(kind=real_kind), intent(in) :: s(4,4)

    real(kind=real_kind) :: ds(4,4,2)
    real(kind=real_kind) :: rrearth

    integer i,j,l,m,n
    real(kind=real_kind) ::  dscontra(4,4,2)

    real(kind=real_kind), dimension(4,4) :: elem_mp
    real(kind=real_kind), dimension(2,2,4,4) :: elem_D
    dscontra=0
    do n=1,4
       do m=1,4
          do j=1,4
             dscontra(m,n,1)=dscontra(m,n,1)-(elem_mp(m,j)*s(m,j)*Dvv(n,j) )*rrearth
             dscontra(m,n,2)=dscontra(m,n,2)+(elem_mp(j,n)*s(j,n)*Dvv(m,j) )*rrearth
          enddo
       enddo
    enddo

    do j=1,4
       do i=1,4
          ds(i,j,1)=(elem_D(1,1,i,j)*dscontra(i,j,1) + elem_D(1,2,i,j)*dscontra(i,j,2))
          ds(i,j,2)=(elem_D(2,1,i,j)*dscontra(i,j,1) + elem_D(2,2,i,j)*dscontra(i,j,2))
       enddo
    enddo
    end function my_curl_sphere_wk_testcov

  function my_gradient_sphere_wk_testcov(s,dvv,elem,rrearth, elem_metdet, elem_D, elem_mp, elem_metinv) result(ds)
    !type (derivative_t)              :: deriv
   use element_mod, only: element_t
    real(kind=real_kind) :: dvv(4,4)
    type (element_t)                 :: elem
    real(kind=real_kind), intent(in) :: s(4,4)

    real(kind=real_kind) :: ds(4,4,2)
    real(kind=real_kind) :: rrearth

    integer i,j,l,m,n
    real(kind=real_kind) ::  dscontra(4,4,2)

    real(kind=real_kind), dimension(4,4) :: elem_mp,  elem_metdet
    real(kind=real_kind), dimension(2,2,4,4) :: elem_D, elem_metinv
    
    dscontra=0
    do n=1,4
       do m=1,4
          do j=1,4
             dscontra(m,n,1)=dscontra(m,n,1)-(&
                  (elem_mp(j,n)*elem_metinv(1,1,m,n)*elem_metdet(m,n)*s(j,n)*Dvv(m,j) ) +&
                  (elem_mp(m,j)*elem_metinv(2,1,m,n)*elem_metdet(m,n)*s(m,j)*Dvv(n,j) ) &
                  ) *rrearth

             dscontra(m,n,2)=dscontra(m,n,2)-(&
                  (elem_mp(j,n)*elem_metinv(1,2,m,n)*elem_metdet(m,n)*s(j,n)*Dvv(m,j) ) +&
                  (elem_mp(m,j)*elem_metinv(2,2,m,n)*elem_metdet(m,n)*s(m,j)*Dvv(n,j) ) &
                  ) *rrearth
          enddo
       enddo
    enddo
    do j=1,4
       do i=1,4
          ds(i,j,1)=(elem_D(1,1,i,j)*dscontra(i,j,1) + elem_D(1,2,i,j)*dscontra(i,j,2))
          ds(i,j,2)=(elem_D(2,1,i,j)*dscontra(i,j,1) + elem_D(2,2,i,j)*dscontra(i,j,2))
       enddo
    enddo

    end function my_gradient_sphere_wk_testcov


  subroutine advance_hypervis_lf(edge3,elem,hvcoord,hybrid,deriv,nm1,n0,nt,nets,nete,dt2)
  !
  !  take one timestep of:
  !          u(:,:,:,np) = u(:,:,:,np) +  dt2*nu*laplacian**order ( u )
  !          T(:,:,:,np) = T(:,:,:,np) +  dt2*nu_s*laplacian**order ( T )
  !
  !
  !  For correct scaling, dt2 should be the same 'dt2' used in the leapfrog advace
  !
  !
  use dimensions_mod, only : np, np, nlev
  use control_mod, only : nu, nu_div, nu_s, hypervis_order, hypervis_subcycle, nu_p, nu_top, psurf_vis
  use hybrid_mod, only : hybrid_t
  use hybvcoord_mod, only : hvcoord_t
  use element_mod, only : element_t
  use derivative_mod, only : derivative_t, laplace_sphere_wk, vlaplace_sphere_wk
  use edge_mod, only : EdgeBuffer_t, edgevpack, edgevunpack
  use bndry_mod, only : bndry_exchangev
  use viscosity_mod, only : biharmonic_wk
  use physical_constants, only: Cp
!  use time_mod, only : TimeLevel_t
  implicit none

  type (hybrid_t)      , intent(in) :: hybrid
  type (element_t)     , intent(inout), target :: elem(:)
  type (EdgeBuffer_t)  , intent(inout) :: edge3
  type (derivative_t)  , intent(in) :: deriv
  type (hvcoord_t), intent(in)      :: hvcoord
!  type (TimeLevel_t)   , intent(in) :: tl

  real (kind=real_kind) :: dt2
  integer :: nets,nete

  ! local
  real (kind=real_kind) :: nu_scale, dpdn,dpdn0, nu_scale_top,nu_ratio
  integer :: k,kptr,i,j,ie,ic,n0,nt,nm1
  real (kind=real_kind), dimension(np,np,2,nlev,nets:nete)      :: vtens
  real (kind=real_kind), dimension(np,np,nlev,nets:nete)        :: ptens
  real (kind=real_kind), dimension(np,np,nets:nete) :: pstens
  real (kind=real_kind), dimension(np,np,nlev) :: p
  real (kind=real_kind), dimension(np,np) :: dXdp


! NOTE: PGI compiler bug: when using spheremp, rspheremp and ps as pointers to elem(ie)% members,
  !       data is incorrect (offset by a few numbers actually)
  !       removed for now.
  !       real (kind=real_kind), dimension(:,:), pointer :: spheremp,rspheremp
  !       real (kind=real_kind), dimension(:,:,:), pointer   :: ps

  real (kind=real_kind), dimension(np,np) :: lap_p
  real (kind=real_kind), dimension(np,np,2) :: lap_v
  real (kind=real_kind) :: v1,v2,dt,heating,utens_tmp,vtens_tmp,ptens_tmp


  if (nu_s == 0 .and. nu == 0 .and. nu_p==0 ) return;
  call t_barrierf('sync_advance_hypervis_lf', hybrid%par%comm)
  call t_startf('advance_hypervis_lf')

! for non-leapfrog,nt=n0=nmt
!
!  nm1 = tl%nm1   ! heating term uses U,V at average of nt and nm1 levels
!  n0 = tl%n0     ! timelevel used for ps scaling.  use n0 for leapfrog.
!  nt = tl%np1    ! apply viscosity to this timelevel  (np1)


  dt=dt2/hypervis_subcycle
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  regular viscosity
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (hypervis_order == 1) then
     if (nu_p>0) stop 'ERROR: hypervis_order == 1 not coded for nu_p>0'
     do ic=1,hypervis_subcycle
        do ie=nets,nete

#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k,lap_p,lap_v,deriv,i,j)
#endif
           do k=1,nlev
              lap_p=laplace_sphere_wk(elem(ie)%state%T(:,:,k,nt),deriv,elem(ie),var_coef=.false.)
              lap_v=vlaplace_sphere_wk(elem(ie)%state%v(:,:,:,k,nt),deriv,elem(ie),var_coef=.false.)
              ! advace in time.  (note: DSS commutes with time stepping, so we
              ! can time advance and then DSS.  this has the advantage of
              ! not letting any discontinuties accumulate in p,v via roundoff
              do j=1,np
                 do i=1,np
                    elem(ie)%state%T(i,j,k,nt)=elem(ie)%state%T(i,j,k,nt)*elem(ie)%spheremp(i,j)  +  dt*nu_s*lap_p(i,j)
                    elem(ie)%state%v(i,j,1,k,nt)=elem(ie)%state%v(i,j,1,k,nt)*elem(ie)%spheremp(i,j) + dt*nu*lap_v(i,j,1)
                    elem(ie)%state%v(i,j,2,k,nt)=elem(ie)%state%v(i,j,2,k,nt)*elem(ie)%spheremp(i,j) + dt*nu*lap_v(i,j,2)
                 enddo
              enddo
           enddo

           kptr=0
           call edgeVpack(edge3, elem(ie)%state%T(:,:,:,nt),nlev,kptr,elem(ie)%desc)
           kptr=nlev
           call edgeVpack(edge3,elem(ie)%state%v(:,:,:,:,nt),2*nlev,kptr,elem(ie)%desc)
        enddo

        call bndry_exchangeV(hybrid,edge3)

        do ie=nets,nete

           kptr=0
           call edgeVunpack(edge3, elem(ie)%state%T(:,:,:,nt), nlev, kptr, elem(ie)%desc)
           kptr=nlev
           call edgeVunpack(edge3, elem(ie)%state%v(:,:,:,:,nt), 2*nlev, kptr, elem(ie)%desc)

           ! apply inverse mass matrix
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k,i,j)
#endif
           do k=1,nlev
              do j=1,np
                 do i=1,np
                    elem(ie)%state%T(i,j,k,nt)=elem(ie)%rspheremp(i,j)*elem(ie)%state%T(i,j,k,nt)
                    elem(ie)%state%v(i,j,1,k,nt)=elem(ie)%rspheremp(i,j)*elem(ie)%state%v(i,j,1,k,nt)
                    elem(ie)%state%v(i,j,2,k,nt)=elem(ie)%rspheremp(i,j)*elem(ie)%state%v(i,j,2,k,nt)
                 enddo
              enddo
           enddo
        enddo
#ifdef DEBUGOMP
#if (! defined ELEMENT_OPENMP)
!$OMP BARRIER
#endif
#endif
     enddo  ! subcycle
  endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  hyper viscosity
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (hypervis_order == 2) then
     nu_ratio = nu_div/nu ! possibly weight div component more inside biharmonc_wk
     do ic=1,hypervis_subcycle
        call biharmonic_wk(elem,pstens,ptens,vtens,deriv,edge3,hybrid,nt,nets,nete,nu_ratio)
        do ie=nets,nete

           nu_scale=1
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k,i,j,lap_p,lap_v,nu_scale_top,dpdn,dpdn0,nu_scale,utens_tmp,vtens_tmp,ptens_tmp)
#endif
           do k=1,nlev
              ! advace in time.
              ! note: DSS commutes with time stepping, so we can time advance and then DSS.
              ! note: weak operators alreayd have mass matrix "included"

              ! add regular diffusion in top 3 layers:
              if (nu_top>0 .and. k<=3) then
                 lap_p=laplace_sphere_wk(elem(ie)%state%T(:,:,k,nt),deriv,elem(ie),var_coef=.false.)
                 lap_v=vlaplace_sphere_wk(elem(ie)%state%v(:,:,:,k,nt),deriv,elem(ie),var_coef=.false.)
              endif
              nu_scale_top = 1
              if (k==1) nu_scale_top=4
              if (k==2) nu_scale_top=2

              do j=1,np
                 do i=1,np
                    if (psurf_vis==0) then
                       ! normalize so as to conserve IE  (not needed when using p-surface viscosity)
                       ! scale velosity by 1/rho (normalized to be O(1))
                       ! dp/dn = O(ps0)*O(delta_eta) = O(ps0)/O(nlev)
                       dpdn = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
                            ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(i,j,n0)  ! nt ?
                       dpdn0 = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
                            ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*hvcoord%ps0
                       nu_scale = dpdn0/dpdn
                    endif

                    ! biharmonic terms need a negative sign:
                    if (nu_top>0 .and. k<=3) then
                       utens_tmp=(-nu*vtens(i,j,1,k,ie) + nu_scale_top*nu_top*lap_v(i,j,1))
                       vtens_tmp=(-nu*vtens(i,j,2,k,ie) + nu_scale_top*nu_top*lap_v(i,j,2))
                       ptens_tmp=nu_scale*(-nu_s*ptens(i,j,k,ie) + nu_scale_top*nu_top*lap_p(i,j) )
                    else
                       utens_tmp=-nu*vtens(i,j,1,k,ie)
                       vtens_tmp=-nu*vtens(i,j,2,k,ie)
                       ptens_tmp=-nu_scale*nu_s*ptens(i,j,k,ie)
                    endif

                    ptens(i,j,k,ie) = ptens_tmp
                    vtens(i,j,1,k,ie)=utens_tmp
                    vtens(i,j,2,k,ie)=vtens_tmp
                 enddo
              enddo
           enddo

           pstens(:,:,ie)  =  -nu_p*pstens(:,:,ie)
           kptr=0
           call edgeVpack(edge3, ptens(:,:,:,ie),nlev,kptr,elem(ie)%desc)
           kptr=nlev
           call edgeVpack(edge3,vtens(:,:,:,:,ie),2*nlev,kptr,elem(ie)%desc)
           kptr=3*nlev
           call edgeVpack(edge3,pstens(:,:,ie),1,kptr,elem(ie)%desc)
        enddo


        call bndry_exchangeV(hybrid,edge3)

        do ie=nets,nete

           kptr=0
           call edgeVunpack(edge3, ptens(:,:,:,ie), nlev, kptr, elem(ie)%desc)
           kptr=nlev
           call edgeVunpack(edge3, vtens(:,:,:,:,ie), 2*nlev, kptr, elem(ie)%desc)
           kptr=3*nlev
           call edgeVunpack(edge3, pstens(:,:,ie), 1, kptr, elem(ie)%desc)

           if (psurf_vis == 1 ) then
              ! apply p-surface correction
              do k=1,nlev
                 p(:,:,k)   = hvcoord%hyam(k)*hvcoord%ps0 + hvcoord%hybm(k)*elem(ie)%state%ps_v(:,:,nt)
              enddo
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k,dXdp)
#endif
              do k=1,nlev
                 if (k.eq.1) then
                    ! no correction needed
                 else if (k.eq.nlev) then
                    ! one-sided difference
                    dXdp = (elem(ie)%state%T(:,:,k,nt) - elem(ie)%state%T(:,:,k-1,nt)) / &
                        (p(:,:,k)-p(:,:,k-1))
                    ptens(:,:,k,ie) = ptens(:,:,k,ie) - dXdp(:,:)*hvcoord%hybm(k)*pstens(:,:,ie)
                 else
                    dXdp = (elem(ie)%state%T(:,:,k+1,nt) - elem(ie)%state%T(:,:,k-1,nt)) / &
                         (p(:,:,k+1)-p(:,:,k-1))
                    ptens(:,:,k,ie) = ptens(:,:,k,ie) - dXdp(:,:)*hvcoord%hybm(k)*pstens(:,:,ie)
                 endif
              enddo
           endif


           ! apply inverse mass matrix, accumulate tendencies
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k,i,j,v1,v2,heating)
#endif
           do k=1,nlev
              do j=1,np
                 do i=1,np

                    elem(ie)%state%v(i,j,1,k,nt)=elem(ie)%state%v(i,j,1,k,nt) + &
                         dt*elem(ie)%rspheremp(i,j)*vtens(i,j,1,k,ie)
                    elem(ie)%state%v(i,j,2,k,nt)=elem(ie)%state%v(i,j,2,k,nt) +  &
                         dt*elem(ie)%rspheremp(i,j)*vtens(i,j,2,k,ie)

                    ! better E conservation if we use v after adding in vtens:
                    v1=.5*(elem(ie)%state%v(i,j,1,k,nt)+elem(ie)%state%v(i,j,1,k,nm1))
                    v2=.5*(elem(ie)%state%v(i,j,2,k,nt)+elem(ie)%state%v(i,j,2,k,nm1))
                    heating = (vtens(i,j,1,k,ie)*v1  + vtens(i,j,2,k,ie)*v2 )

                    elem(ie)%state%T(i,j,k,nt)=elem(ie)%state%T(i,j,k,nt)     + &
                         dt*elem(ie)%rspheremp(i,j)*(cp*ptens(i,j,k,ie) - heating)/cp

                 enddo
              enddo
           enddo
           elem(ie)%state%ps_v(:,:,nt)=elem(ie)%state%ps_v(:,:,nt) + dt*elem(ie)%rspheremp(:,:)*pstens(:,:,ie)
        enddo
#ifdef DEBUGOMP
#if (! defined ELEMENT_OPENMP)
!$OMP BARRIER
#endif
#endif
     enddo
  endif

  call t_stopf('advance_hypervis_lf')

  end subroutine advance_hypervis_lf

  ! added by conghui
  subroutine my_Virtual_Temperature1d(Tin,rin, Rwater_vapor, Rgas, Tv)
    real (kind=8),intent(in) :: Tin
    real (kind=8),intent(in) :: rin
    real (kind=8),intent(in) :: Rwater_vapor
    real (kind=8),intent(in) :: Rgas
    real (kind=8),intent(out) :: Tv

    Tv = Tin*(1 + (Rwater_vapor/Rgas - 1.0)*rin)
  end subroutine

  subroutine my_Virtual_Specific_Heat(rin, Cp, Cpwater_vapor, Cp_star)
    real (kind=8),intent(in) :: rin
    real (kind=8),intent(in) :: Cp
    real (kind=8),intent(in) :: Cpwater_vapor
    real (kind=8),intent(out) :: Cp_star

    Cp_star = Cp*(1.0 + (Cpwater_vapor/Cp - 1.0)*rin)
  end subroutine

  subroutine my_preq_hydrostatic(phi,phis,T_v,p,dp, rgas)
    implicit none

    !------------------------------Arguments---------------------------------------------------------------
    real(kind=8), intent(out) :: phi(4,4,constLev)
    real(kind=8), intent(in) :: phis(4,4)
    real(kind=8), intent(in) :: T_v(4,4,constLev)
    real(kind=8), intent(in) :: p(4,4,constLev)
    real(kind=8), intent(in) :: dp(4,4,constLev)
    real(kind=8), intent(in) :: rgas

    !---------------------------Local workspace-----------------------------
    integer i,j,k                         ! longitude, level indices
    real(kind=8) Hkk,Hkl          ! diagonal term of energy conversion matrix
    real(kind=8), dimension(4,4,constLev+1) :: phii       ! Geopotential at interfaces

      phii(:,:,constLev+1) = 0
      do k=constLev,1,-1
        do j=1,4   !   Loop inversion (AAM)
          do i=1,4
            ! hkk = dp*ckk
            hkk = dp(i,j,k)*0.5d0/p(i,j,k)
            hkl = 2*hkk
            phii(i,j,k) = phii(i,j,k+1) + Rgas*T_v(i,j,k)*hkl
            phi(i,j,k) = phis(i,j) + phii(i,j,k+1) + Rgas*T_v(i,j,k)*hkk
          enddo
        enddo
      enddo
  end subroutine

  subroutine my_preq_omega_ps(omega_p,p,vgrad_p,divdp,suml)
    implicit none

    real(kind=8), intent(in) :: divdp(4,4,constLev)      ! divergence
    real(kind=8), intent(in) :: vgrad_p(4,4,constLev) ! v.grad(p)
    real(kind=8), intent(in) :: p(4,4,constLev)     ! layer thicknesses (pressure)
    real(kind=8), intent(out):: omega_p(4,4,constLev)   ! vertical pressure velocity
    real(kind=8), intent(in) :: suml(4,4,constLev+1)      ! partial sum over l = (1, k-1)

    integer i,j,k                         ! longitude, level indices
    real(kind=8) term             ! one half of basic term in omega/p summation
    real(kind=8) Ckk,Ckl          ! diagonal term of energy conversion matrix

    do k=1,constLev
      do j=1,4   !   Loop inversion (AAM)
        do i=1,4
          ckk = 0.5d0/p(i,j,k)
          ckl = 2*ckk
          omega_p(i,j,k) = vgrad_p(i,j,k)/p(i,j,k) - ckl*suml(i,j,k) - ckk*divdp(i,j,k)
        end do
      end do
    enddo
  end subroutine

    subroutine my_gradient_sphere(s,deriv_Dvv,Dinv,ds,my_rrearth)
!
!   i4ut s:  scalar
!   output  ds: spherical gradient of s, lat-lon coordinates
!

    real(kind=8), intent(in), dimension(4,4)          :: deriv_Dvv
    real(kind=8), intent(in), dimension(2,2,4,4) :: Dinv
    real(kind=8), intent(in) :: s(4,4)
    real(kind=8), intent(in) :: my_rrearth

    real(kind=8), intent(out) :: ds(4,4,2)

    integer i
    integer j
    integer l

    real(kind=real_kind) ::  dsdx00
    real(kind=real_kind) ::  dsdy00
    real(kind=real_kind) ::  v1(4,4),v2(4,4)

    do j=1,4
       do l=1,4
          dsdx00=0.0d0
          dsdy00=0.0d0
          do i=1,4
             dsdx00 = dsdx00 + deriv_Dvv(i,l  )*s(i,j  )
             dsdy00 = dsdy00 + deriv_Dvv(i,l  )*s(j  ,i)
          end do
          v1(l  ,j  ) = dsdx00*my_rrearth
          v2(j  ,l  ) = dsdy00*my_rrearth
       end do
    end do
    ! convert covarient to latlon
    do j=1,4
       do i=1,4
          ds(i,j,1)=Dinv(1,1,i,j)*v1(i,j) + Dinv(2,1,i,j)*v2(i,j)
          ds(i,j,2)=Dinv(1,2,i,j)*v1(i,j) + Dinv(2,2,i,j)*v2(i,j)
       enddo
    enddo

  end subroutine

  subroutine my_divergence_sphere(v,deriv_Dvv,elem_metdet,elem_rmetdet,elem_Dinv,my_rrearth, div)
    real(kind=8), intent(in) :: v(4,4,2)  ! in lat-lon coordinates
    real(kind=8), intent(in):: deriv_Dvv(4,4)
    real(kind=8), intent(in):: elem_metdet(4,4)
    real(kind=8), intent(in):: elem_rmetdet(4,4)
    real(kind=8), intent(in):: elem_Dinv(2,2,4,4)
    real(kind=8), intent(in):: my_rrearth
    real(kind=8), intent(out):: div(4,4)

    ! Local
    integer i
    integer j
    integer l

    real(kind=8) ::  dudx00
    real(kind=8) ::  dvdy00
    real(kind=8) ::  gv(4,4,2),vvtemp(4,4)

    ! convert to contra variant form and multiply by g
    do j=1,4
       do i=1,4
          gv(i,j,1)=elem_metdet(i,j)*(elem_Dinv(1,1,i,j)*v(i,j,1) + elem_Dinv(1,2,i,j)*v(i,j,2))
          gv(i,j,2)=elem_metdet(i,j)*(elem_Dinv(2,1,i,j)*v(i,j,1) + elem_Dinv(2,2,i,j)*v(i,j,2))
       enddo
    enddo

    ! compute d/dx and d/dy
    do j=1,4
       do l=1,4
          dudx00=0.0d0
          dvdy00=0.0d0
          do i=1,4
             dudx00 = dudx00 + deriv_Dvv(i,l  )*gv(i,j  ,1)
             dvdy00 = dvdy00 + deriv_Dvv(i,l  )*gv(j  ,i,2)
          end do
          div(l  ,j  ) = dudx00
          vvtemp(j  ,l  ) = dvdy00
       end do
    end do

    do j=1,4
       do i=1,4
          div(i,j)=(div(i,j)+vvtemp(i,j))*(elem_rmetdet(i,j)*my_rrearth)
       end do
    end do
  end subroutine

  subroutine my_vorticity_sphere(v,deriv_Dvv,elem_D, elem_rmetdet,my_rrearth,vort)
    real(kind=8), intent(in) :: deriv_Dvv(4,4)
    real(kind=8), intent(in) :: elem_D(2,2,4,4)
    real(kind=8), intent(in) :: elem_rmetdet(4,4)
    real(kind=8), intent(in) :: v(4,4,2)
    real(kind=8), intent(in) :: my_rrearth
    real(kind=8), intent(out) :: vort(4,4)

    integer i
    integer j
    integer l

    real(kind=8) ::  dvdx00
    real(kind=8) ::  dudy00
    real(kind=8) ::  vco(4,4,2)
    real(kind=8) ::  vtemp(4,4)

    ! convert to covariant form
    do j=1,4
      do i=1,4
        vco(i,j,1)=(elem_D(1,1,i,j)*v(i,j,1) + elem_D(2,1,i,j)*v(i,j,2))
        vco(i,j,2)=(elem_D(1,2,i,j)*v(i,j,1) + elem_D(2,2,i,j)*v(i,j,2))
      enddo
    enddo

    do j=1,4
      do l=1,4

        dudy00=0.0d0
        dvdx00=0.0d0

        do i=1,4
          dvdx00 = dvdx00 + deriv_Dvv(i,l  )*vco(i,j  ,2)
          dudy00 = dudy00 + deriv_Dvv(i,l  )*vco(j  ,i,1)
        enddo

        vort(l  ,j  ) = dvdx00
        vtemp(j  ,l  ) = dudy00
      enddo
    enddo

    do j=1,4
      do i=1,4
        vort(i,j)=(vort(i,j)-vtemp(i,j))*(elem_rmetdet(i,j)*my_rrearth)
      end do
    end do
  end subroutine

  subroutine my_advance_acc(np1,nm1,n0,qn0,dt2,my_rrearth,elem,hvcoord,&
       deriv,nets,nete,compute_diagnostics,eta_ave_w,use_cpstar,rsplit, &
       cp, cpwater_vapor, Rgas, kappa, Rwater_vapor,elem_array)

  use element_mod, only : element_t
  use derivative_mod, only : derivative_t
  use hybvcoord_mod, only : hvcoord_t



  implicit none
  !!$ACC routine memstack(size=409600; name=private_liaojf) reuseldm(size=49152)

  integer, parameter :: np = 4
  integer, parameter :: nlev = constLev 
  integer :: np1,nm1,n0,qn0,nets,nete
  integer :: use_cpstar,rsplit
  real (kind=8) :: cp, cpwater_vapor, Rgas, kappa, Rwater_vapor
  real*8 :: dt2
  logical  :: compute_diagnostics

  real(kind=8), intent(in) :: my_rrearth
  type (hvcoord_t)     , intent(in) :: hvcoord
  type (element_t),dimension(nets:nete), intent(inout) :: elem
  integer(kind=8), dimension(23,nets:nete), intent(inout) :: elem_array
  type (derivative_t)  , intent(in) :: deriv
  real (kind=8) :: eta_ave_w  ! weighting for eta_dot_dpdn mean flux

  ! local
  real (kind=8), dimension(np,np)   :: omega_p
  !real (kind=8), dimension(np,np,nlev)   :: T_v
  real (kind=8), dimension(np,np)   :: divdp
  real (kind=8), dimension(np,np,2)      :: vtemp
  real (kind=8), dimension(np,np)        :: vgrad_T
  real (kind=8), dimension(np,np)        :: Ephi
  real (kind=8), dimension(np,np,2)      :: grad_ps
  real (kind=8), dimension(np,np,2) :: grad_p
  real (kind=8), dimension(np,np)   :: vort
  real (kind=8), dimension(np,np,nlev)   :: p
  !real (kind=8), dimension(np,np,nlev)   :: dp
  real (kind=8), dimension(np,np)   :: vgrad_p
  real (kind=8), dimension(np,np)   :: suml 
    real(kind=8), dimension(np,np,constLev+1) :: phii       ! Geopotential at interfaces

  real (kind=8) ::  cp2,cp_ratio,E,de,Qt,v1,v2, hkk
  real (kind=8) ::  glnps1,glnps2,gpterm
  integer :: i,j,k,kptr,ie

  real (kind=8), dimension(np,np,nlev)    :: elem_derived_phi
  pointer(elem_derived_phi_ptr, elem_derived_phi)

  real (kind=8), dimension(nlev)          :: hvcoord_hyai
  pointer(hvcoord_hyai_ptr, hvcoord_hyai)

  real (kind=8) :: hvcoord_ps0

  real (kind=8), dimension(np,np,nlev) :: elem_state_dp3d_n0
  pointer(elem_state_dp3d_n0_ptr, elem_state_dp3d_n0)

  real (kind=8), dimension(np,np,nlev) :: elem_state_dp3d_np1
  pointer(elem_state_dp3d_np1_ptr, elem_state_dp3d_np1)

  real (kind=8), dimension(np,np,nlev) :: elem_state_dp3d_nm1
  pointer(elem_state_dp3d_nm1_ptr, elem_state_dp3d_nm1)

  real (kind=8), dimension(np,np,2,2) :: elem_Dinv
  pointer(elem_Dinv_ptr, elem_Dinv)

  real (kind=8), dimension(np,np,2,nlev) :: elem_state_v_n0
  pointer(elem_state_v_n0_ptr, elem_state_v_n0)

  real (kind=8), dimension(np,np,2,nlev) :: elem_state_v_np1
  pointer(elem_state_v_np1_ptr, elem_state_v_np1)

  real (kind=8), dimension(np,np,2,nlev) :: elem_state_v_nm1
  pointer(elem_state_v_nm1_ptr, elem_state_v_nm1)

  real (kind=8), dimension(np,np,2,nlev) :: elem_derived_vn0
  pointer(elem_derived_vn0_ptr, elem_derived_vn0)

  real (kind=8), dimension(np,np,nlev)   :: elem_state_Qdp_1_qn0
  pointer(elem_state_Qdp_1_qn0_ptr, elem_state_Qdp_1_qn0)

  real (kind=8), dimension(np,np,nlev)   :: elem_state_T_n0
  pointer(elem_state_T_n0_ptr, elem_state_T_n0)

  real (kind=8), dimension(np,np,nlev)   :: elem_state_T_np1
  pointer(elem_state_T_np1_ptr, elem_state_T_np1)

  real (kind=8), dimension(np,np,nlev)   :: elem_state_T_nm1
  pointer(elem_state_T_nm1_ptr, elem_state_T_nm1)

  real (kind=8), dimension(np,np) :: elem_state_phis
  pointer(elem_state_phis_ptr, elem_state_phis)

  real (kind=8), dimension(np,np,nlev)   :: elem_derived_omega_p
  pointer(elem_derived_omega_p_ptr, elem_derived_omega_p)

  real (kind=8), dimension(np,np,nlev)   :: elem_derived_pecnd
  pointer(elem_derived_pecnd_ptr, elem_derived_pecnd)

  real (kind=8), dimension(np,np)   :: elem_fcor
  pointer(elem_fcor_ptr, elem_fcor)

  real (kind=8), dimension(np,np)   :: elem_spheremp
  pointer(elem_spheremp_ptr, elem_spheremp)

  real (kind=8), dimension(np,np)   :: elem_state_ps_v_np1
  pointer(elem_state_ps_v_np1_ptr, elem_state_ps_v_np1)

  real (kind=8), dimension(np,np)   :: elem_state_ps_v_nm1
  pointer(elem_state_ps_v_nm1_ptr, elem_state_ps_v_nm1)

  real (kind=8), dimension(np,np)   :: deriv_Dvv
  pointer(deriv_Dvv_ptr, deriv_Dvv)

  real (kind=8), dimension(np,np)   :: elem_metdet
  pointer(elem_metdet_ptr, elem_metdet)

  real (kind=8), dimension(np,np)   :: elem_rmetdet
  pointer(elem_rmetdet_ptr, elem_rmetdet)

  real (kind=8), dimension(2,2,np,np) :: elem_D
  pointer(elem_D_ptr, elem_D)

  real (kind=8) :: conghui_tmp_var

  !write(*,*) "we built this city"

  hvcoord_ps0      = hvcoord%ps0
  hvcoord_hyai_ptr = loc(hvcoord%hyai)
  deriv_Dvv_ptr    = loc(deriv%Dvv)

  !$ACC PARALLEL LOOP local(omega_p,divdp,vtemp,vgrad_T,Ephi,grad_ps,grad_p,vort,p,vgrad_p, suml, phii) copyin(elem_array,hvcoord_hyai,deriv_dvv)  annotate(entire(deriv_dvv))
  do ie=nets,nete

     elem_derived_phi_ptr          = elem_array(1,ie)
     elem_state_dp3d_n0_ptr        = elem_array(2,ie)
     elem_state_dp3d_np1_ptr       = elem_array(3,ie)
     elem_state_dp3d_nm1_ptr       = elem_array(4,ie)
     elem_Dinv_ptr                 = elem_array(5,ie)
     elem_state_v_n0_ptr           = elem_array(6,ie)
     elem_state_v_np1_ptr          = elem_array(7,ie)
     elem_state_v_nm1_ptr          = elem_array(8,ie)
     elem_derived_vn0_ptr          = elem_array(9,ie)
     !elem_state_Qdp_1_qn0_ptr      = elem_array(10,ie)
     elem_state_T_n0_ptr           = elem_array(11,ie)
     elem_state_T_np1_ptr          = elem_array(12,ie)
     elem_state_T_nm1_ptr          = elem_array(13,ie)
     elem_state_phis_ptr           = elem_array(14,ie)
     elem_derived_omega_p_ptr      = elem_array(15,ie)
     elem_derived_pecnd_ptr        = elem_array(16,ie)
     elem_fcor_ptr                 = elem_array(17,ie)
     elem_spheremp_ptr             = elem_array(18,ie)
     elem_state_ps_v_np1_ptr       = elem_array(19,ie)
     elem_state_ps_v_nm1_ptr       = elem_array(20,ie)
     elem_metdet_ptr               = elem_array(21,ie)
     elem_rmetdet_ptr              = elem_array(22,ie)
     elem_D_ptr                    = elem_array(23,ie)

     !$ACC DATA COPYIN(elem_Dinv, elem_state_phis, elem_fcor, elem_spheremp, elem_state_ps_v_nm1, elem_metdet, elem_rmetdet, elem_D) copyout(elem_state_ps_v_np1) annotate(entire(deriv_dvv, elem_state_ps_v_nm1, elem_spheremp, elem_state_ps_v_np1, elem_fcor, elem_Dinv, elem_state_phis, elem_rmetdet, elem_D, elem_metdet))
     !phi => elem_derived_phi
     suml(:,:)= 0
     phii(:,:,nlev+1) = 0
     do k=1,nlev
           !$ACC DATA COPYIN(elem_state_dp3d_n0(*,*,k),elem_state_dp3d_n0(*,*,k-1))
           if (k /= 1) then
             p(:,:,k)=p(:,:,k-1) + elem_state_dp3d_n0(:,:,k-1)*0.5d0 + elem_state_dp3d_n0(:,:,k)*0.5d0
           else
             p(:,:,1)=hvcoord_hyai(1)*hvcoord_ps0 + elem_state_dp3d_n0(:,:,1)*0.5d0
           endif
           !$ACC END DATA
     enddo
     do k=nlev,1,-1
        !$ACC DATA COPYIN(elem_state_dp3d_n0(*,*,k),elem_state_T_n0(*,*,k))
        phii(:,:,k) = phii(:,:,k+1) + Rgas*elem_state_T_n0(:,:,k)*elem_state_dp3d_n0(:,:,k)/p(:,:,k)
        !$ACC END DATA
     enddo
     do k = 1, nlev
       !$ACC DATA COPYIN(elem_state_dp3d_nm1(*,*,k), elem_state_v_n0(*,*,*,k), elem_state_v_nm1(*,*,*,k), elem_state_T_nm1(*,*,k), elem_derived_pecnd(*,*,k)) copyout(elem_derived_phi(*,*,k), elem_state_dp3d_np1(*,*,k), elem_state_v_np1(*,*,*,k), elem_state_T_np1(*,*,k)) copy(elem_derived_vn0(*,*,*,k), elem_derived_omega_p(*,*,k),elem_state_dp3d_n0(*,*,k),elem_state_T_n0(*,*,k))
       call my_gradient_sphere(p(:,:,k),deriv_Dvv(:,:),elem_Dinv(:,:,:,:),grad_p(:,:,:),my_rrearth)
       do j=1,np
         do i=1,np
           v1 = elem_state_v_n0(i,j,1,k)
           v2 = elem_state_v_n0(i,j,2,k)
           vgrad_p(i,j) = (v1*grad_p(i,j,1) + v2*grad_p(i,j,2))
           vtemp(i,j,1) = v1*elem_state_dp3d_n0(i,j,k)
           vtemp(i,j,2) = v2*elem_state_dp3d_n0(i,j,k)
         end do
       end do
       elem_derived_vn0(:,:,:,k)=elem_derived_vn0(:,:,:,k)+eta_ave_w*vtemp(:,:,:)
       ! edit by conghui
       call my_divergence_sphere(vtemp(:,:,:),deriv_Dvv(:,:),elem_metdet(:,:),elem_rmetdet(:,:),elem_Dinv(:,:,:,:),my_rrearth, divdp(:,:))
       call my_vorticity_sphere(elem_state_v_n0(:,:,:,k),deriv_Dvv(:,:),elem_D(:,:,:,:), elem_rmetdet(:,:),my_rrearth,vort(:,:))
           elem_derived_phi(:,:,k) = elem_state_phis(:,:) + phii(:,:,k+1) + Rgas*elem_state_T_n0(:,:,k)*elem_state_dp3d_n0(:,:,k)*0.5d0/p(:,:,k)
       do j=1,4   !   Loop inversion (AAM)
         do i=1,4
           hkk = 0.5d0/p(i,j,k)
           omega_p(i,j) = vgrad_p(i,j)/p(i,j,k) - 2*hkk*suml(i,j) - hkk*divdp(i,j)
         end do
       end do
       suml(:,:) = suml(:,:) + divdp(:,:)
       elem_derived_omega_p(:,:,k)      = elem_derived_omega_p(:,:,k) + eta_ave_w*omega_p(:,:)
       do j=1,np
         do i=1,np
           v1     = elem_state_v_n0(i,j,1,k)
           v2     = elem_state_v_n0(i,j,2,k)
           E = 0.5D0*( v1*v1 + v2*v2 )
           Ephi(i,j)=E+elem_derived_phi(i,j,k)+elem_derived_pecnd(i,j,k)
          end do
        end do
        ! ================================================
        ! compute gradp term (ps/p)*(dp/dps)*T
        ! ================================================
        !vtemp(:,:,:)   = gradient_sphere(elem_state_T_n0(:,:,k),deriv,elem_Dinv)
        call my_gradient_sphere(elem_state_T_n0(:,:,k),deriv_Dvv(:,:),elem_Dinv(:,:,:,:),vtemp(:,:,:),my_rrearth)
        do j=1,np
          do i=1,np
            v1     = elem_state_v_n0(i,j,1,k)
            v2     = elem_state_v_n0(i,j,2,k)
            vgrad_T(i,j) =  v1*vtemp(i,j,1) + v2*vtemp(i,j,2)
          end do
        end do


        ! vtemp = grad ( E + PHI )
        !vtemp = gradient_sphere(Ephi(:,:),deriv,elem_Dinv)
        call my_gradient_sphere(Ephi(:,:),deriv_Dvv(:,:),elem_Dinv(:,:,:,:),vtemp(:,:,:),my_rrearth)

        do j=1,np
          do i=1,np
            gpterm = elem_state_T_n0(i,j,k)/p(i,j,k)
            glnps1 = Rgas*gpterm*grad_p(i,j,1)
            glnps2 = Rgas*gpterm*grad_p(i,j,2)
            v1     = elem_state_v_n0(i,j,1,k)
            v2     = elem_state_v_n0(i,j,2,k)

            grad_p(i,j,1) = + v2*(elem_fcor(i,j) + vort(i,j)) - vtemp(i,j,1) - glnps1
            grad_p(i,j,2) = - v1*(elem_fcor(i,j) + vort(i,j)) - vtemp(i,j,2) - glnps2
            p(i,j,k)  = -vgrad_T(i,j) + kappa*elem_state_T_n0(i,j,k)*omega_p(i,j)
          end do
        end do
        elem_state_v_np1(:,:,1,k) = elem_spheremp(:,:)*( elem_state_v_nm1(:,:,1,k) + dt2*grad_p(:,:,1) )
        elem_state_v_np1(:,:,2,k) = elem_spheremp(:,:)*( elem_state_v_nm1(:,:,2,k) + dt2*grad_p(:,:,2) )
        elem_state_T_np1(:,:,k) = elem_spheremp(:,:)*(elem_state_T_nm1(:,:,k) + dt2*p(:,:,k))
        elem_state_dp3d_np1(:,:,k) = elem_spheremp(:,:)* (elem_state_dp3d_nm1(:,:,k)-dt2*divdp(:,:) )
        !$ACC END DATA
      enddo
      elem_state_ps_v_np1(:,:) = elem_spheremp(:,:)*( elem_state_ps_v_nm1(:,:))
      !$ACC END DATA
    enddo
    !$ACC END PARALLEL LOOP
  end subroutine my_advance_acc

  subroutine compute_and_apply_rhs(np1,nm1,n0,qn0,dt2,elem,hvcoord,hybrid,&
       deriv,nets,nete,compute_diagnostics,eta_ave_w)
  ! ===================================
  ! compute the RHS, accumulate into u(np1) and apply DSS
  !
  !           u(np1) = u(nm1) + dt2*DSS[ RHS(u(n0)) ]
  !
  ! This subroutine is normally called to compute a leapfrog timestep
  ! but by adjusting np1,nm1,n0 and dt2, many other timesteps can be
  ! accomodated.  For example, setting nm1=np1=n0 this routine will
  ! take a forward euler step, overwriting the input with the output.
  !
  !    qn0 = timelevel used to access Qdp() in order to compute virtual Temperature
  !          qn0=-1 for the dry case
  !
  ! if  dt2<0, then the DSS'd RHS is returned in timelevel np1
  !
  ! Combining the RHS and DSS pack operation in one routine
  ! allows us to fuse these two loops for more cache reuse
  !
  ! Combining the dt advance and DSS unpack operation in one routine
  ! allows us to fuse these two loops for more cache reuse
  !
  ! note: for prescribed velocity case, velocity will be computed at
  ! "real_time", which should be the time of timelevel n0.
  !
  !
  ! ===================================
  use kinds, only : real_kind
  use dimensions_mod, only : np, np, nlev
  use hybrid_mod, only : hybrid_t
  use element_mod, only : element_t
  use derivative_mod, only : derivative_t, divergence_sphere, gradient_sphere, vorticity_sphere
  use edge_mod, only : edgevpack, edgevunpack
  use bndry_mod, only : bndry_exchangev
  use control_mod, only : moisture, qsplit, use_cpstar, rsplit
  use hybvcoord_mod, only : hvcoord_t

  use physical_constants, only : cp, cpwater_vapor, Rgas, kappa, rrearth, Rwater_vapor
  use physics_mod, only : virtual_specific_heat, virtual_temperature
  use prim_si_mod, only : preq_vertadv, preq_omega_ps, preq_hydrostatic


  implicit none
  integer :: np1,nm1,n0,qn0,nets,nete
  real*8 :: dt2
  logical  :: compute_diagnostics

  type (hvcoord_t)     , intent(in) :: hvcoord
  type (hybrid_t)      , intent(in) :: hybrid
  type (element_t)     , intent(inout), target :: elem(:)
  type (derivative_t)  , intent(in) :: deriv
  real (kind=real_kind) :: eta_ave_w  ! weighting for eta_dot_dpdn mean flux

  ! local
  real (kind=real_kind), pointer, dimension(:,:)      :: ps         ! surface pressure for current tiime level
  real (kind=real_kind), pointer, dimension(:,:,:)   :: phi

  real (kind=real_kind), dimension(np,np,nlev)   :: omega_p
  real (kind=real_kind), dimension(np,np,nlev)   :: T_v
  real (kind=real_kind), dimension(np,np,nlev)   :: divdp
  real (kind=real_kind), dimension(np,np,nlev+1)   :: eta_dot_dpdn  ! half level vertical velocity on p-grid
  real (kind=real_kind), dimension(np,np)      :: sdot_sum   ! temporary field
  real (kind=real_kind), dimension(np,np,2)    :: vtemp     ! generic gradient storage
  real (kind=real_kind), dimension(np,np)      :: vgrad_T    ! v.grad(T)
  real (kind=real_kind), dimension(np,np)      :: Ephi       ! kinetic energy + PHI term
  real (kind=real_kind), dimension(np,np,2)      :: grad_ps    ! lat-lon coord version
  real (kind=real_kind), dimension(np,np,2,nlev) :: grad_p
  real (kind=real_kind), dimension(np,np,nlev)   :: vort       ! vorticity
  real (kind=real_kind), dimension(np,np,nlev)   :: p          ! pressure
  real (kind=real_kind), dimension(np,np,nlev)   :: dp         ! delta pressure
  real (kind=real_kind), dimension(np,np,nlev)   :: rdp        ! inverse of delta pressure
  real (kind=real_kind), dimension(np,np,nlev)   :: T_vadv     ! temperature vertical advection
  real (kind=real_kind), dimension(np,np,nlev)   :: vgrad_p    ! v.grad(p)
  real (kind=real_kind), dimension(np,np,nlev+1) :: ph               ! half level pressures on p-grid
  real (kind=real_kind), dimension(np,np,2,nlev) :: v_vadv   ! velocity vertical advection
  real (kind=real_kind) ::  Kappa_star(np,np,nlev)
  real (kind=real_kind) ::  vtens1(np,np,nlev)
  real (kind=real_kind) ::  vtens2(np,np,nlev)
  real (kind=real_kind) ::  ttens(np,np,nlev)

  real (kind=real_kind) ::  cp2,cp_ratio,E,de,Qt,v1,v2
  real (kind=real_kind) ::  glnps1,glnps2,gpterm
  integer :: i,j,k,kptr,ie
  integer(kind=8) :: count_start, count_stop, count_rate, count_max
  integer(kind=8), dimension(23,nets:nete) :: elem_array

  real (kind=8), dimension(np,np)   :: elem_state_ps_v_np1
  pointer(elem_state_ps_v_np1_ptr, elem_state_ps_v_np1)
  
  real(kind=8), dimension(nlev) :: edge_buf_5, edge_buf_6, edge_buf_7, edge_buf_8
  pointer(edge_buf_5_ptr, edge_buf_5)
  pointer(edge_buf_6_ptr, edge_buf_6)
  pointer(edge_buf_7_ptr, edge_buf_7)
  pointer(edge_buf_8_ptr, edge_buf_8)

  real(kind=8), dimension(nlev) :: edge_buf_in_1, edge_buf_in_2, edge_buf_in_3, edge_buf_in_4
  pointer(edge_buf_in_1_ptr, edge_buf_in_1)
  pointer(edge_buf_in_2_ptr, edge_buf_in_2)
  pointer(edge_buf_in_3_ptr, edge_buf_in_3)
  pointer(edge_buf_in_4_ptr, edge_buf_in_4)


  real(kind=8), dimension(nlev) :: edge_buf_is_1, edge_buf_is_2, edge_buf_is_3, edge_buf_is_4
  pointer(edge_buf_is_1_ptr, edge_buf_is_1)
  pointer(edge_buf_is_2_ptr, edge_buf_is_2)
  pointer(edge_buf_is_3_ptr, edge_buf_is_3)
  pointer(edge_buf_is_4_ptr, edge_buf_is_4)

  real(kind=8), dimension(nlev) :: edge_buf_ie_1, edge_buf_ie_2, edge_buf_ie_3, edge_buf_ie_4
  pointer(edge_buf_ie_1_ptr, edge_buf_ie_1)
  pointer(edge_buf_ie_2_ptr, edge_buf_ie_2)
  pointer(edge_buf_ie_3_ptr, edge_buf_ie_3)
  pointer(edge_buf_ie_4_ptr, edge_buf_ie_4)

  real(kind=8), dimension(nlev) :: edge_buf_iw_1, edge_buf_iw_2, edge_buf_iw_3, edge_buf_iw_4
  pointer(edge_buf_iw_1_ptr, edge_buf_iw_1)
  pointer(edge_buf_iw_2_ptr, edge_buf_iw_2)
  pointer(edge_buf_iw_3_ptr, edge_buf_iw_3)
  pointer(edge_buf_iw_4_ptr, edge_buf_iw_4)
  
  logical, dimension(8) :: elem_desc_reverse
  pointer(elem_desc_reverse_ptr, elem_desc_reverse)

  integer, dimension(8) :: elem_desc_putmapP
  pointer(elem_desc_putmapP_ptr, elem_desc_putmapP)
  
  call t_barrierf('sync_compute_and_apply_rhs', hybrid%par%comm)
  call t_startf('compute_and_apply_rhs')

  !call system_clock(count_start, count_rate, count_max)
  do ie = nets,nete

     elem_array(1,ie)      = loc(elem(ie)%derived%phi)
     elem_array(2,ie)      = loc(elem(ie)%state%dp3d(:,:,:,n0))
     elem_array(3,ie)      = loc(elem(ie)%state%dp3d(:,:,:,np1))
     elem_array(4,ie)      = loc(elem(ie)%state%dp3d(:,:,:,nm1))
     elem_array(5,ie)      = loc(elem(ie)%Dinv)
     elem_array(6,ie)      = loc(elem(ie)%state%v(:,:,:,:,n0))
     elem_array(7,ie)      = loc(elem(ie)%state%v(:,:,:,:,np1))
     elem_array(8,ie)      = loc(elem(ie)%state%v(:,:,:,:,nm1))
     elem_array(9,ie)      = loc(elem(ie)%derived%vn0)
     elem_array(10,ie)      = loc(elem(ie)%state%Qdp(:,:,:,1,qn0))
     elem_array(11,ie)      = loc(elem(ie)%state%T(:,:,:,n0))
     elem_array(12,ie)      = loc(elem(ie)%state%T(:,:,:,np1))
     elem_array(13,ie)      = loc(elem(ie)%state%T(:,:,:,nm1))
     elem_array(14,ie)      = loc(elem(ie)%state%phis)
     elem_array(15,ie)      = loc(elem(ie)%derived%omega_p)
     elem_array(16,ie)      = loc(elem(ie)%derived%pecnd)
     elem_array(17,ie)      = loc(elem(ie)%fcor)
     elem_array(18,ie)      = loc(elem(ie)%spheremp)
     elem_array(19,ie)      = loc(elem(ie)%state%ps_v(:,:,np1))
     elem_array(20,ie)      = loc(elem(ie)%state%ps_v(:,:,nm1))
     elem_array(21,ie)      = loc(elem(ie)%metdet)
     elem_array(22,ie)      = loc(elem(ie)%rmetdet)
     elem_array(23,ie)      = loc(elem(ie)%D)

  enddo
  call t_startf('my_advance')
   call my_advance_acc(np1,nm1,n0,qn0,dt2,rrearth,elem,hvcoord,&
     deriv,nets,nete,compute_diagnostics,eta_ave_w,use_cpstar,rsplit, &
     cp, cpwater_vapor, Rgas, kappa, Rwater_vapor,elem_array)
  call t_stopf('my_advance')
  !call system_clock(count_stop, count_rate, count_max)
   !write(*,*) "My_advance_acc time=",count_stop-count_start
   do ie=nets,nete
      elem_desc_reverse_ptr            = loc(elem(ie)%desc%reverse)
      elem_desc_putmapP_ptr            = loc(elem(ie)%desc%putmapP)
      edge_buf_5_ptr    = loc(edge3p1%buf(:, elem_desc_putmapP(5)+1))
      edge_buf_6_ptr    = loc(edge3p1%buf(:, elem_desc_putmapP(6)+1))
      edge_buf_7_ptr    = loc(edge3p1%buf(:, elem_desc_putmapP(7)+1))
      edge_buf_8_ptr    = loc(edge3p1%buf(:, elem_desc_putmapP(8)+1))
      edge_buf_in_1_ptr = loc(edge3p1%buf(:, elem_desc_putmapP(north)+1))
      edge_buf_in_2_ptr = loc(edge3p1%buf(:, elem_desc_putmapP(north)+2))
      edge_buf_in_3_ptr = loc(edge3p1%buf(:, elem_desc_putmapP(north)+3))
      edge_buf_in_4_ptr = loc(edge3p1%buf(:, elem_desc_putmapP(north)+4))
      edge_buf_is_1_ptr = loc(edge3p1%buf(:, elem_desc_putmapP(south)+1))
      edge_buf_is_2_ptr = loc(edge3p1%buf(:, elem_desc_putmapP(south)+2))
      edge_buf_is_3_ptr = loc(edge3p1%buf(:, elem_desc_putmapP(south)+3))
      edge_buf_is_4_ptr = loc(edge3p1%buf(:, elem_desc_putmapP(south)+4))
      edge_buf_ie_1_ptr = loc(edge3p1%buf(:, elem_desc_putmapP(east)+1))
      edge_buf_ie_2_ptr = loc(edge3p1%buf(:, elem_desc_putmapP(east)+2))
      edge_buf_ie_3_ptr = loc(edge3p1%buf(:, elem_desc_putmapP(east)+3))
      edge_buf_ie_4_ptr = loc(edge3p1%buf(:, elem_desc_putmapP(east)+4))
      edge_buf_iw_1_ptr = loc(edge3p1%buf(:, elem_desc_putmapP(west)+1))
      edge_buf_iw_2_ptr = loc(edge3p1%buf(:, elem_desc_putmapP(west)+2))
      edge_buf_iw_3_ptr = loc(edge3p1%buf(:, elem_desc_putmapP(west)+3))
      edge_buf_iw_4_ptr = loc(edge3p1%buf(:, elem_desc_putmapP(west)+4))
      
      !elem_state_ps_v_np1 = elem_array(19,ie)

     kptr=0
     call edgeVpack(edge3p1, elem(ie)%state%ps_v(:,:,np1),1,kptr,elem(ie)%desc)
      
     !call my_edgeVpack_acc(elem(ie)%state%T(:,:,:,np1) , nlev , 1 , elem_desc_putmapP(:), &
     !    elem_desc_reverse, &
     !    edge_buf_5, edge_buf_6, edge_buf_7, edge_buf_8, &
     !    edge_buf_in_1, edge_buf_in_2, edge_buf_in_3, edge_buf_in_4, &
     !    edge_buf_is_1, edge_buf_is_2, edge_buf_is_3, edge_buf_is_4, &
     !    edge_buf_ie_1, edge_buf_ie_2, edge_buf_ie_3, edge_buf_ie_4, &
     !    edge_buf_iw_1, edge_buf_iw_2, edge_buf_iw_3, edge_buf_iw_4)

     kptr=1
     call edgeVpack(edge3p1, elem(ie)%state%T(:,:,:,np1),nlev,kptr,elem(ie)%desc)

     kptr=nlev+1
     call edgeVpack(edge3p1, elem(ie)%state%v(:,:,:,:,np1),2*nlev,kptr,elem(ie)%desc)

     kptr=kptr+2*nlev
     call edgeVpack(edge3p1, elem(ie)%state%dp3d(:,:,:,np1),nlev,kptr,elem(ie)%desc)
  end do

  ! =============================================================
    ! Insert communications here: for shared memory, just a single
  ! sync is required
  ! =============================================================

  call bndry_exchangeV(hybrid,edge3p1)

  do ie=nets,nete
     ! ===========================================================
     ! Unpack the edges for vgrad_T and v tendencies...
     ! ===========================================================
     kptr=0
     call edgeVunpack(edge3p1, elem(ie)%state%ps_v(:,:,np1), 1, kptr, elem(ie)%desc)

     kptr=1
     call edgeVunpack(edge3p1, elem(ie)%state%T(:,:,:,np1), nlev, kptr, elem(ie)%desc)

     kptr=nlev+1
     call edgeVunpack(edge3p1, elem(ie)%state%v(:,:,:,:,np1), 2*nlev, kptr, elem(ie)%desc)

        kptr=kptr+2*nlev
        call edgeVunpack(edge3p1, elem(ie)%state%dp3d(:,:,:,np1),nlev,kptr,elem(ie)%desc)

     ! ====================================================
     ! Scale tendencies by inverse mass matrix
     ! ====================================================

     do k=1,nlev
        elem(ie)%state%T(:,:,k,np1)   = elem(ie)%rspheremp(:,:)*elem(ie)%state%T(:,:,k,np1)
        elem(ie)%state%v(:,:,1,k,np1) = elem(ie)%rspheremp(:,:)*elem(ie)%state%v(:,:,1,k,np1)
        elem(ie)%state%v(:,:,2,k,np1) = elem(ie)%rspheremp(:,:)*elem(ie)%state%v(:,:,2,k,np1)
     end do

     if (rsplit>0) then
        ! vertically lagrangian: complete dp3d timestep:
        do k=1,nlev
           elem(ie)%state%dp3d(:,:,k,np1)= elem(ie)%rspheremp(:,:)*elem(ie)%state%dp3d(:,:,k,np1)
        enddo
        ! when debugging: also update ps_v
        !elem(ie)%state%ps_v(:,:,np1) = elem(ie)%rspheremp(:,:)*elem(ie)%state%ps_v(:,:,np1)
     else
        ! vertically eulerian: complete ps_v timestep:
        elem(ie)%state%ps_v(:,:,np1) = elem(ie)%rspheremp(:,:)*elem(ie)%state%ps_v(:,:,np1)
     endif


  end do
#ifdef DEBUGOMP
#if (! defined ELEMENT_OPENMP)
!$OMP BARRIER
#endif
#endif
  call t_stopf('compute_and_apply_rhs')
  end subroutine compute_and_apply_rhs



  subroutine my_edgeVpack_acc(v,vlyr,kptr,desc_putmapP, &
      desc_reverse, &
      edge_buf_5,edge_buf_6,edge_buf_7,edge_buf_8, &
      edge_buf_in_1,edge_buf_in_2,edge_buf_in_3,edge_buf_in_4, &
      edge_buf_is_1,edge_buf_is_2,edge_buf_is_3,edge_buf_is_4, &
      edge_buf_ie_1,edge_buf_ie_2,edge_buf_ie_3,edge_buf_ie_4, &
      edge_buf_iw_1,edge_buf_iw_2,edge_buf_iw_3,edge_buf_iw_4)


   !integer, public, parameter :: west  = 1
   !integer, public, parameter :: east  = 2
   !integer, public, parameter :: south = 3
   !integer, public, parameter :: north = 4

   !integer, public, parameter :: swest = 5
   !integer, public, parameter :: seast = 6
   !integer, public, parameter :: nwest = 7
   !integer, public, parameter :: neast = 8
   !integer, public, parameter :: max_corner_elem = 1
   implicit none
    integer,                                         intent(in)     :: vlyr
    real (kind=8),dimension(4,4,constLev),                 intent(in)     :: v
    integer,                                         intent(in)     :: kptr
    integer, dimension(8),                           intent(in)     :: desc_putmapP
    logical,dimension(8),                            intent(in)     :: desc_reverse
    real(kind=8), dimension(constLev),                     intent(inout)  :: edge_buf_5, edge_buf_6, edge_buf_7, edge_buf_8
    real(kind=8), dimension(constLev),                     intent(inout)  :: edge_buf_in_1, edge_buf_in_2, edge_buf_in_3, edge_buf_in_4
    real(kind=8), dimension(constLev),                     intent(inout)  :: edge_buf_is_1, edge_buf_is_2, edge_buf_is_3, edge_buf_is_4
    real(kind=8), dimension(constLev),                     intent(inout)  :: edge_buf_ie_1, edge_buf_ie_2, edge_buf_ie_3, edge_buf_ie_4
    real(kind=8), dimension(constLev),                     intent(inout)  :: edge_buf_iw_1, edge_buf_iw_2, edge_buf_iw_3, edge_buf_iw_4


    ! Local variables
    logical, parameter :: UseUnroll = .TRUE.
    integer :: i,k,ir,ll

    integer :: is,ie,in,iw

  !integer :: rank,ierr
  !include "mpif.h"
  !call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

    is = desc_putmapP(3)
    ie = desc_putmapP(2)
    in = desc_putmapP(4)
    iw = desc_putmapP(1)
       do k=1,vlyr
          !do i=1,np,2
             !edge%buf(kptr+k,is+i)   = v(i  ,1 ,k)
             !edge%buf(kptr+k,is+i+1) = v(i+1,1 ,k)
             !edge%buf(kptr+k,ie+i)   = v(np ,i ,k)
             !edge%buf(kptr+k,ie+i+1) = v(np ,i+1 ,k)
             !edge%buf(kptr+k,in+i)   = v(i  ,np,k)
             !edge%buf(kptr+k,in+i+1) = v(i+1  ,np,k)
             !edge%buf(kptr+k,iw+i)   = v(1  ,i ,k)
             !edge%buf(kptr+k,iw+i+1) = v(1  ,i+1 ,k)
          !enddo
          edge_buf_in_1(k) = v(1,4,k)
          edge_buf_in_2(k) = v(2,4,k)
          edge_buf_in_3(k) = v(3,4,k)
          edge_buf_in_4(k) = v(4,4,k)
          edge_buf_is_1(k) = v(1,1,k)
          edge_buf_is_2(k) = v(2,1,k)
          edge_buf_is_3(k) = v(3,1,k)
          edge_buf_is_4(k) = v(4,1,k)
          edge_buf_ie_1(k) = v(4,1,k)
          edge_buf_ie_2(k) = v(4,2,k)
          edge_buf_ie_3(k) = v(4,3,k)
          edge_buf_ie_4(k) = v(4,4,k)
          edge_buf_iw_1(k) = v(1,1,k)
          edge_buf_iw_2(k) = v(1,2,k)
          edge_buf_iw_3(k) = v(1,3,k)
          edge_buf_iw_4(k) = v(1,4,k)
       end do



    if(desc_reverse(3)) then
       do k=1,vlyr
             edge_buf_is_4(k)=v(1,1,k)
             edge_buf_is_3(k)=v(2,1,k)
             edge_buf_is_2(k)=v(3,1,k)
             edge_buf_is_1(k)=v(4,1,k)
        enddo
    endif

    if(desc_reverse(2)) then
       do k=1,vlyr
             edge_buf_ie_4(k)=v(4,1,k)
             edge_buf_ie_3(k)=v(4,2,k)
             edge_buf_ie_2(k)=v(4,3,k)
             edge_buf_ie_1(k)=v(4,4,k)
       enddo
    endif

    if(desc_reverse(4)) then
       do k=1,vlyr
             edge_buf_in_4(k)=v(1,4,k)
             edge_buf_in_3(k)=v(2,4,k)
             edge_buf_in_2(k)=v(3,4,k)
             edge_buf_in_1(k)=v(4,4,k)
       enddo
    endif

    if(desc_reverse(1)) then
       do k=1,vlyr
             edge_buf_iw_4(k)=v(1,1,k)
             edge_buf_iw_3(k)=v(1,2,k)
             edge_buf_iw_2(k)=v(1,3,k)
             edge_buf_iw_1(k)=v(1,4,k)
       enddo
    endif


! SWEST
    ll = 5
        if (desc_putmapP(ll) /= -1) then
            do k=1,vlyr
                edge_buf_5(k)=v(1  ,1 ,k)
            end do
        end if

! SEAST
    ll = 6
        if (desc_putmapP(ll) /= -1) then
            do k=1,vlyr
                edge_buf_6(k)=v(4 ,1 ,k)
            end do
        end if

! NEAST
    ll = 8
        if (desc_putmapP(ll) /= -1) then
            do k=1,vlyr
                edge_buf_8(k)=v(4 ,4,k)
            end do
        end if

! NWEST
     ll = 7
        if (desc_putmapP(ll) /= -1) then
            do k=1,vlyr
                edge_buf_7(k)=v(1  ,4,k)
            end do
        end if


  end subroutine my_edgeVpack_acc


  subroutine compute_frontogenesis(frontga,frontgf,psurf_ref,hvcoord,tl,&
       elem,deriv,hybrid,nets,nete)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! compute frontogenis function F
  !   F =  -gradth dot C
  ! with:
  !   theta  = potential temperature
  !   gradth = grad(theta)
  !   C = ( gradth dot grad ) U
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  use physical_constants, only : kappa
  use dimensions_mod, only : np,nlev
  use hybvcoord_mod, only : hvcoord_t
  use element_mod, only : element_t
  use hybrid_mod, only : hybrid_t
  use edge_mod, only : EdgeBuffer_t, edgevpack, edgevunpack
  use bndry_mod, only : bndry_exchangev
  use derivative_mod, only : derivative_t, gradient_sphere, ugradv_sphere
  type (hybrid_t)      , intent(in) :: hybrid
  type (hvcoord_t)      , intent(in) :: hvcoord
  type (element_t)     , intent(inout), target :: elem(:)
  type (derivative_t)  , intent(in) :: deriv
  integer :: nets,nete
  integer :: tl ! timelevel to use
  real(kind=real_kind), intent(in) ::  psurf_ref  ! from CAM's ref_pres module
  real(kind=real_kind), intent(out) ::  frontga(np,np,nlev,nets:nete)
  real(kind=real_kind), intent(out) ::  frontgf(np,np,nlev,nets:nete)

  ! local
  integer :: k,kptr,i,j,ie,component
  real(kind=real_kind)  ::  gradth(np,np,2,nlev,nets:nete)  ! grad(theta)
  real(kind=real_kind)  :: p(np,np)        ! pressure at mid points
  real(kind=real_kind)  :: theta(np,np)    ! potential temperature at mid points
  real(kind=real_kind)  ::  C(np,np,2)

  do ie=nets,nete
     do k=1,nlev
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
        ! potential temperature: theta = T (p/p0)^kappa
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

        ! pressure at mid points
        p(:,:)   = hvcoord%hyam(k)*hvcoord%ps0 + hvcoord%hybm(k)*elem(ie)%state%ps_v(:,:,tl)
        theta(:,:) = elem(ie)%state%T(:,:,k,tl)*(psurf_ref / p(:,:))**kappa
        gradth(:,:,:,k,ie) = gradient_sphere(theta,deriv,elem(ie)%Dinv)

        ! compute C = (grad(theta) dot grad ) u
        C(:,:,:) = ugradv_sphere(gradth(:,:,:,k,ie), elem(ie)%state%v(:,:,:,k,tl),deriv,elem(ie))

        ! gradth dot C
        frontgf(:,:,k,ie) = -( C(:,:,1)*gradth(:,:,1,k,ie) +  C(:,:,2)*gradth(:,:,2,k,ie)  )

        ! apply mass matrix
        gradth(:,:,1,k,ie)=gradth(:,:,1,k,ie)*elem(ie)%spheremp(:,:)
        gradth(:,:,2,k,ie)=gradth(:,:,2,k,ie)*elem(ie)%spheremp(:,:)
        frontgf(:,:,k,ie)=frontgf(:,:,k,ie)*elem(ie)%spheremp(:,:)

     enddo
     ! pack
     call edgeVpack(edge3p1, frontgf(:,:,:,ie),nlev,0,elem(ie)%desc)
     call edgeVpack(edge3p1, gradth(:,:,:,:,ie),2*nlev,nlev,elem(ie)%desc)
  enddo
  call bndry_exchangeV(hybrid,edge3p1)
  do ie=nets,nete
     call edgeVunpack(edge3p1, frontgf(:,:,:,ie),nlev,0,elem(ie)%desc)
     call edgeVunpack(edge3p1, gradth(:,:,:,:,ie),2*nlev,nlev,elem(ie)%desc)
     ! apply inverse mass matrix,
     do k=1,nlev
        gradth(:,:,1,k,ie)=gradth(:,:,1,k,ie)*elem(ie)%rspheremp(:,:)
        gradth(:,:,2,k,ie)=gradth(:,:,2,k,ie)*elem(ie)%rspheremp(:,:)
        frontgf(:,:,k,ie)=frontgf(:,:,k,ie)*elem(ie)%rspheremp(:,:)

        frontga(:,:,k,ie) = atan2 ( gradth(:,:,2,k,ie) , gradth(:,:,1,k,ie) + 1.e-10_real_kind )
     enddo
  enddo
  end subroutine compute_frontogenesis








  subroutine smooth_phis(phis,elem,hybrid,deriv,nets,nete,minf,numcycle)
  use dimensions_mod, only : np, np, nlev
  use control_mod, only : smooth_phis_nudt
  use hybrid_mod, only : hybrid_t
  use edge_mod, only : EdgeBuffer_t, edgevpack, edgevunpack, edgevunpackmax, edgevunpackmin
  use bndry_mod, only : bndry_exchangev
  use element_mod, only : element_t
  use derivative_mod, only : derivative_t , laplace_sphere_wk
  use time_mod, only : TimeLevel_t
  implicit none

  integer :: nets,nete
  real (kind=real_kind), dimension(np,np,nets:nete), intent(inout)   :: phis
  type (hybrid_t)      , intent(in) :: hybrid
  type (element_t)     , intent(inout), target :: elem(:)
  type (derivative_t)  , intent(in) :: deriv
  real (kind=real_kind), intent(in)   :: minf
  integer,               intent(in) :: numcycle

  ! local
  real (kind=real_kind), dimension(np,np,nets:nete) :: pstens
  real (kind=real_kind), dimension(nets:nete) :: pmin,pmax
  real (kind=real_kind) :: mx,mn
  integer :: nt,ie,ic,i,j,order,order_max, iuse


  ! compute local element neighbor min/max
  do ie=nets,nete
     pstens(:,:,ie)=minval(phis(:,:,ie))
     call edgeVpack(edge3p1,pstens(:,:,ie),1,0,elem(ie)%desc)
  enddo
  call bndry_exchangeV(hybrid,edge3p1)
  do ie=nets,nete
     call edgeVunpackMin(edge3p1, pstens(:,:,ie), 1, 0, elem(ie)%desc)
     pmin(ie)=minval(pstens(:,:,ie))
  enddo
  do ie=nets,nete
     pstens(:,:,ie)=maxval(phis(:,:,ie))
     call edgeVpack(edge3p1,pstens(:,:,ie),1,0,elem(ie)%desc)
  enddo
  call bndry_exchangeV(hybrid,edge3p1)
  do ie=nets,nete
     call edgeVunpackMax(edge3p1, pstens(:,:,ie), 1, 0, elem(ie)%desc)
     pmax(ie)=maxval(pstens(:,:,ie))
  enddo


  do ic=1,numcycle

     pstens=phis

     ! order = 1   laplacian
     ! order = 2   grad^4 (need to add a negative sign)
     ! order = 3   grad^6
     ! order = 4   grad^8 (need to add a negative sign)
     order_max = 1

     do order=1,order_max-1
        do ie=nets,nete
           pstens(:,:,ie)=laplace_sphere_wk(pstens(:,:,ie),deriv,elem(ie),var_coef=.true.)
           call edgeVpack(edge3p1,pstens(:,:,ie),1,0,elem(ie)%desc)
        enddo
        call bndry_exchangeV(hybrid,edge3p1)
        do ie=nets,nete
           call edgeVunpack(edge3p1, pstens(:,:,ie), 1, 0, elem(ie)%desc)
           pstens(:,:,ie)=pstens(:,:,ie)*elem(ie)%rspheremp(:,:)
        enddo
#ifdef DEBUGOMP
#if (! defined ELEMENT_OPENMP)
!$OMP BARRIER
#endif
#endif
     enddo
     do ie=nets,nete
        pstens(:,:,ie)=laplace_sphere_wk(pstens(:,:,ie),deriv,elem(ie),var_coef=.true.)
     enddo
     if (mod(order_max,2)==0) pstens=-pstens

     do ie=nets,nete
        !  ps(t+1) = ps(t) + Minv * DSS * M * RHS
        !  ps(t+1) = Minv * DSS * M [ ps(t) +  RHS ]
        ! but output of biharminc_wk is of the form M*RHS.  rewrite as:
        !  ps(t+1) = Minv * DSS * M [ ps(t) +  M*RHS/M ]
        ! so we can apply limiter to ps(t) +  (M*RHS)/M
#if 0
        mn=minval(phis(:,:,ie))
        mx=maxval(phis(:,:,ie))
        iuse = numcycle/2  ! apply first half of iterations
#else
        mn=pmin(ie)
        mx=pmax(ie)
        iuse = numcycle+1  ! always apply
#endif
        phis(:,:,ie)=phis(:,:,ie) + &
           smooth_phis_nudt*pstens(:,:,ie)/elem(ie)%spheremp(:,:)


        ! remove new extrema.  could use conservative reconstruction from advection
        ! but no reason to conserve mean PHI.
        if (ic < iuse) then
        do i=1,np
        do j=1,np
           if (phis(i,j,ie)>mx) phis(i,j,ie)=mx
           if (phis(i,j,ie)<mn) phis(i,j,ie)=mn
        enddo
        enddo
        endif


        ! user specified minimum
        do i=1,np
        do j=1,np
           if (phis(i,j,ie)<minf) phis(i,j,ie)=minf
        enddo
        enddo

        phis(:,:,ie)=phis(:,:,ie)*elem(ie)%spheremp(:,:)
        call edgeVpack(edge3p1,phis(:,:,ie),1,0,elem(ie)%desc)
     enddo
     call bndry_exchangeV(hybrid,edge3p1)
     do ie=nets,nete
        call edgeVunpack(edge3p1, phis(:,:,ie), 1, 0, elem(ie)%desc)
        phis(:,:,ie)=phis(:,:,ie)*elem(ie)%rspheremp(:,:)
     enddo
#ifdef DEBUGOMP
#if (! defined ELEMENT_OPENMP)
!$OMP BARRIER
#endif
#endif
  enddo
  end subroutine smooth_phis

  subroutine overwrite_SEdensity(elem, fvm, hybrid,nets,nete, np1)
    use fvm_reconstruction_mod, only: reconstruction
    use fvm_filter_mod, only: monotonic_gradient_cart, recons_val_cart
    use dimensions_mod, only : np, nlev, nc,nhe
    use control_mod, only : smooth_phis_nudt
    use hybrid_mod, only : hybrid_t
    use edge_mod, only : EdgeBuffer_t, edgevpack, edgevunpack, edgevunpackmax, edgevunpackmin
    use bndry_mod, only : bndry_exchangev
    use element_mod, only : element_t
    use derivative_mod, only : derivative_t , laplace_sphere_wk
    use time_mod, only : TimeLevel_t
    use fvm_control_volume_mod, only : fvm_struct
    use spelt_mod, only : spelt_struct


    type (element_t) , intent(inout)        :: elem(:)

#if defined(_SPELT)
      type(spelt_struct), intent(inout) :: fvm(:)
#else
      type(fvm_struct), intent(inout) :: fvm(:)
#endif
    type (hybrid_t), intent(in)           :: hybrid  ! distributed parallel structure (shared)

    integer, intent(in)                     :: nets  ! starting thread element number (private)
    integer, intent(in)                     :: nete  ! ending thread element number   (private)
    integer, intent(in)                     :: np1
    integer :: ie, k

    real (kind=real_kind)             :: xp,yp, tmpval
    integer                           :: i, j,ix, jy, starti,endi,tmpi

    real (kind=real_kind), dimension(5,1-nhe:nc+nhe,1-nhe:nc+nhe)      :: recons

    if ((nc .ne. 4) .or. (np .ne. 4)) then
      if(hybrid%masterthread) then
        print *,"You are in OVERWRITE SE AIR DENSITY MODE"
        print *,"This only works for nc=4 and np=4"
        print *,"Write a new search algorithm or pay $10000!"
      endif
      stop
    endif
#if defined(_FVM)
    do ie=nets,nete
      call reconstruction(fvm(ie)%psc, fvm(ie),recons)
      call monotonic_gradient_cart(fvm(ie)%psc, fvm(ie),recons, elem(ie)%desc)
      do j=1,np
        do i=1,np
          xp=tan(elem(ie)%cartp(i,j)%x)
          yp=tan(elem(ie)%cartp(i,j)%y)
          ix=i
          jy=j
          ! Search index along "x"  (bisection method)
!           starti = 1
!           endi = nc+1
!           do
!              if  ((endi-starti) <=  1)  exit
!              tmpi = (endi + starti)/2
!              if (xp  >  fvm%acartx(tmpi)) then
!                 starti = tmpi
!              else
!                 endi = tmpi
!              endif
!           enddo
!           ix = starti
!
!         ! Search index along "y"
!           starti = 1
!           endi = nc+1
!           do
!              if  ((endi-starti) <=  1)  exit
!              tmpi = (endi + starti)/2
!              if (yp  >  fvm%acarty(tmpi)) then
!                 starti = tmpi
!              else
!                 endi = tmpi
!              endif
!           enddo
!           jy = starti

          call recons_val_cart(fvm(ie)%psc, xp,yp,fvm(ie)%spherecentroid,recons,ix,jy,tmpval)
          elem(ie)%state%ps_v(i,j,np1)=tmpval
        end do
      end do
      elem(ie)%state%ps_v(:,:,np1)=elem(ie)%state%ps_v(:,:,np1)*elem(ie)%spheremp(:,:)
     call edgeVpack(edge3p1,elem(ie)%state%ps_v(:,:,np1),1,0,elem(ie)%desc)
  enddo
  call bndry_exchangeV(hybrid,edge3p1)
  do ie=nets,nete
     call edgeVunpack(edge3p1, elem(ie)%state%ps_v(:,:,np1), 1, 0, elem(ie)%desc)
     elem(ie)%state%ps_v(:,:,np1)=elem(ie)%state%ps_v(:,:,np1)*elem(ie)%rspheremp(:,:)
  enddo
#endif
  end subroutine overwrite_SEdensity


end module prim_advance_mod

! TPP Processed 06af35b4b6e1d5ae0884d64cc00c796d

! TPP Processed 06af35b4b6e1d5ae0884d64cc00c796d

! TPP Processed 06af35b4b6e1d5ae0884d64cc00c796d

! TPP Processed 06af35b4b6e1d5ae0884d64cc00c796d
