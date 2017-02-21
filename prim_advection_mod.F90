#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#define constLev 128
#define partLev 32 

#if 0
SUBROUTINES:
   prim_advec_tracers_remap_rk2()
      SEM 2D RK2 + monotone remap + hyper viscosity
      SEM 2D RK2 can use sign-preserving or monotone reconstruction

Notes on Lagrange+REMAP advection
dynamics will compute mean fluxes, so that (i.e. for qsplit=3)

    dp(t+3)-dp(t) = -3dt div(Udp_sum/3) - 3dt d(eta_dot_dpdn_sum/3)  + 3dt D(dpdiss_sum/3)

Where the floating lagrangian component:
    dp_star(t+3) = dp(t)  -3dt div(Udp_sum/3)  + 3dt D(dpdiss_sum/3)
OR:
    dp_star(t+3) = dp(t+1) + 3dt d( eta_dot_dpdn_ave(t) )


For RK2 advection of Q:  (example of 2 stage RK for tracers):   dtq = qsplit*dt
For consistency, if Q=1
  dp1  = dp(t)- dtq div[ U1 dp(t)]
  dp2  = dp1  - dtq div[ U2 dp1  ]  + 2*dtq D( dpdiss_ave )
  dp*  = (dp(t) + dp2 )/2
       =  dp(t) - dtq  div[ U1 dp(t) + U2 dp1 ]/2   + dtq D( dpdiss_ave )

so we require:
  U1 = Udp_ave / dp(t)
  U2 = Udp_ave / dp1

For tracer advection:
  Qdp1  = Qdp(t)- dtq div[ U1 Qdp(t)]
  Qdp2  = Qdp1  - dtq div[ U2 Qdp1  ]  + 2*dtq D( Q dpdiss_ave )
  Qdp*  = (Qdp(t) + Qdp2 )/2
       =  Qdp(t) - dtq  div[ U1 Qdp(t) + U2 Qdp1 ]   + dtq D( Q dpdiss_ave )

Qdp1:  limit Q, with Q = Qdp1-before-DSS/(dp1-before-DSS)      with dp1 as computed above
Qdp2:  limit Q, with Q = Qdp2-before-DSS/(dp2-before-DSS)      with dp2 as computed above

For dissipation: Q = Qdp1-after-DSS / dp1-after-DSS


last step:
  remap Qdp* to Qdp(t+1)   [ dp_star(t+1) -> dp(t+1) ]



#endif




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Begin GPU remap module  !!
!! by Rick Archibald, 2010  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module vertremap_mod

  !**************************************************************************************
  !
  !  Purpose:
  !        Construct sub-grid-scale polynomials using piecewise spline method with
  !        monotone filters.
  !
  !  References: PCM - Zerroukat et al., Q.J.R. Meteorol. Soc., 2005. (ZWS2005QJR)
  !              PSM - Zerroukat et al., Int. J. Numer. Meth. Fluids, 2005. (ZWS2005IJMF)
  !
  !**************************************************************************************

  use kinds, only                  : real_kind,int_kind
  use dimensions_mod, only         : np,nlev,qsize,nlevp,npsq,ntrac,nc
  use hybvcoord_mod, only          : hvcoord_t
  use element_mod, only            : element_t
  use fvm_control_volume_mod, only : fvm_struct
  use spelt_mod, only              : spelt_struct
  use perf_mod, only               : t_startf, t_stopf  ! _EXTERNAL
  use parallel_mod, only           : abortmp
  use control_mod, only : vert_remap_q_alg

  public my_vertical_remap_acc
  public remap1                  ! remap any field, splines, monotone
  public remap1_nofilter         ! remap any field, splines, no filter
! todo: tweak interface to match remap1 above, rename remap1_ppm:
  public remap_q_ppm             ! remap state%Q, PPM, monotone

  contains

!=======================================================================================================!

!remap_calc_grids computes the vertical pressures and pressure differences for one vertical column for the reference grid
!and for the deformed Lagrangian grid. This was pulled out of each routine since it was a repeated task.
  subroutine my_vertical_remap_acc(elem, hvcoord, ps0, my_nets, my_nete, my_nlev, my_qsize, my_np)
  
  type (element_t), intent(inout)   :: elem(:)
  type (hvcoord_t), intent(in)      :: hvcoord
  !integer(kind=8), intent(in), dimension(7, my_nets:my_nete) :: elem_array
  integer, intent(in) :: my_nets, my_nete, my_nlev, my_qsize, my_np
  real(kind=8) :: ps0


  integer :: ie,i,j,k
  real (kind=real_kind), dimension(my_np,my_np,my_nlev)  :: dp,dp_star
  real (kind=real_kind), dimension(my_np,my_np,my_nlev,2)  :: ttmp

  real(kind=8) :: elem_state_ps_v(my_np,my_np)
  pointer(elem_state_ps_v_ptr, elem_state_ps_v)

  real(kind=8) :: elem_state_dp3d(my_np,my_np,my_nlev)
  pointer(elem_state_dp3d_ptr, elem_state_dp3d)

  real(kind=8) :: elem_state_t(my_np,my_np,my_nlev)
  pointer(elem_state_t_ptr, elem_state_t)

  real(kind=8) :: elem_state_v(my_np,my_np,2,my_nlev)
  pointer(elem_state_v_ptr, elem_state_v)

  real(kind=8) :: elem_state_qdp(my_np,my_np,my_nlev,my_qsize)
  pointer(elem_state_qdp_ptr, elem_state_qdp)

  real(kind=8) :: hvcoord_hyai(nlev+1)
  pointer(hvcoord_hyai_ptr, hvcoord_hyai)

  real(kind=8) :: hvcoord_hybi(nlev+1)
  pointer(hvcoord_hybi_ptr, hvcoord_hybi)

  do ie=nets,nete
        elem(ie)%state%ps_v(:,:,np1) = hvcoord%hyai(1)*hvcoord%ps0 + &
             sum(elem(ie)%state%dp3d(:,:,:,np1),3)
        do k=1,nlev
           dp(:,:,k) = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
                ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(:,:,np1)
           dp_star(:,:,k) = elem(ie)%state%dp3d(:,:,k,np1)
        enddo

        ttmp(:,:,:,1)=elem(ie)%state%t(:,:,:,np1)
        ttmp(:,:,:,1)=ttmp(:,:,:,1)*dp_star
        call remap1(ttmp,np,1,dp_star,dp)
        elem(ie)%state%t(:,:,:,np1)=ttmp(:,:,:,1)/dp

        ttmp(:,:,:,1)=elem(ie)%state%v(:,:,1,:,np1)*dp_star
        ttmp(:,:,:,2)=elem(ie)%state%v(:,:,2,:,np1)*dp_star
        call remap1(ttmp,np,2,dp_star,dp)
        elem(ie)%state%v(:,:,1,:,np1)=ttmp(:,:,:,1)/dp
        elem(ie)%state%v(:,:,2,:,np1)=ttmp(:,:,:,2)/dp

        call remap1(elem(ie)%state%Qdp(:,:,:,:,np1_qdp),np,qsize,dp_star,dp)



  enddo
  end subroutine
subroutine remap_calc_grids( hvcoord , ps , dt , eta_dot_dpdn , p_lag , p_ref , dp_lag , dp_ref )
  implicit none
  type(hvcoord_t)      , intent(in   ) :: hvcoord               !Derived type to hold vertical sigma grid parameters
  real(kind=real_kind) , intent(in   ) :: ps                    !Surface pressure for this column
  real(kind=real_kind) , intent(in   ) :: dt                    !Time step
  real(kind=real_kind) , intent(in   ) :: eta_dot_dpdn(nlev+1)  !Looks like a vertical pressure flux
                                                                !to compute deformed grid spacing
  real(kind=real_kind) , intent(  out) :: p_lag(nlev+1)         !Pressures at interfaces of the Lagrangian deformed grid
  real(kind=real_kind) , intent(  out) :: p_ref(nlev+1)         !Pressures at interfaces of the reference grid
  real(kind=real_kind) , intent(  out) :: dp_lag(nlev)          !Pressure differences on Lagrangian deformed grid
  real(kind=real_kind) , intent(  out) :: dp_ref(nlev)          !Pressure differences on reference grid
  integer :: k                                                  !Iterator
  p_ref(1) = 0  !Both grids have a model top pressure of zero
  p_lag(1) = 0  !Both grids have a model top pressure of zero
  do k = 1 , nlev
    dp_ref(k) = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) ) * hvcoord%ps0 + &
         ( hvcoord%hybi(k+1) - hvcoord%hybi(k) ) * ps  !Reference pressure difference
    ! Lagrangian pressure difference (flux in - flux out over the time step)
    dp_lag(k) = dp_ref(k) + dt * ( eta_dot_dpdn(k+1) - eta_dot_dpdn(k) )
    p_ref(k+1) = p_ref(k) + dp_ref(k) !Pressure at interfaces accumulated using difference over each cell
    p_lag(k+1) = p_lag(k) + dp_lag(k) !Pressure at interfaces accumulated using difference over each cell
  enddo
end subroutine remap_calc_grids

!=======================================================================================================!



subroutine remap1(Qdp,nx,qsize,dp1,dp2)
  ! remap 1 field
  ! input:  Qdp   field to be remapped (NOTE: MASS, not MIXING RATIO)
  !         dp1   layer thickness (source)
  !         dp2   layer thickness (target)
  !
  ! output: remaped Qdp, conserving mass, monotone on Q=Qdp/dp
  !
  implicit none
  integer, intent(in) :: nx,qsize
  real (kind=real_kind), intent(inout) :: Qdp(nx,nx,nlev,qsize)
  real (kind=real_kind), intent(in) :: dp1(nx,nx,nlev),dp2(nx,nx,nlev)
  ! ========================
  ! Local Variables
  ! ========================

  real (kind=real_kind), dimension(nlev+1)    :: rhs,lower_diag,diag,upper_diag,q_diag,zgam,z1c,z2c,zv
  real (kind=real_kind), dimension(nlev)      :: h,Qcol,dy,za0,za1,za2,zarg,zhdp,dp_star,dp_np1
  real (kind=real_kind)  :: f_xm,level1,level2,level3,level4,level5, &
                            peaks_min,peaks_max,tmp_cal,xm,xm_d,zv1,zv2, &
                            zero = 0,one = 1,tiny = 1e-12,qmax = 1d50
  integer(kind=int_kind) :: zkr(nlev+1),filter_code(nlev),peaks,im1,im2,im3,ip1,ip2, &
                            lt1,lt2,lt3,t0,t1,t2,t3,t4,tm,tp,ie,i,ilev,j,jk,k,q
  logical :: abort=.false.
  call t_startf('remap1')

  if (vert_remap_q_alg == 1 .or. vert_remap_q_alg == 2) then
     call remap_Q_ppm(qdp,nx,qsize,dp1,dp2)
     call t_stopf('remap1')
     return
  endif

#if (defined ELEMENT_OPENMP)
!$omp parallel do private(qsize,i,j,z1c,z2c,zv,k,dp_np1,dp_star,Qcol,zkr,ilev) &
!$omp    private(jk,zgam,zhdp,h,zarg,rhs,lower_diag,diag,upper_diag,q_diag,tmp_cal,filter_code) &
!$omp    private(dy,im1,im2,im3,ip1,t1,t2,t3,za0,za1,za2,xm_d,xm,f_xm,t4,tm,tp,peaks,peaks_min) &
!$omp    private(peaks_max,ip2,level1,level2,level3,level4,level5,lt1,lt2,lt3,zv1,zv2)
#endif
  do q=1,qsize
  do i=1,nx
    do j=1,nx

      z1c(1)=0 ! source grid
      z2c(1)=0 ! target grid
      do k=1,nlev
         z1c(k+1)=z1c(k)+dp1(i,j,k)
         z2c(k+1)=z2c(k)+dp2(i,j,k)
      enddo

      zv(1)=0
      do k=1,nlev
        Qcol(k)=Qdp(i,j,k,q)!  *(z1c(k+1)-z1c(k)) input is mass
        zv(k+1) = zv(k)+Qcol(k)
      enddo

      if (ABS(z2c(nlev+1)-z1c(nlev+1)).GE.0.000001) then
        write(6,*) 'SURFACE PRESSURE IMPLIED BY ADVECTION SCHEME'
        write(6,*) 'NOT CORRESPONDING TO SURFACE PRESSURE IN    '
        write(6,*) 'DATA FOR MODEL LEVELS'
        write(6,*) 'PLEVMODEL=',z2c(nlev+1)
        write(6,*) 'PLEV     =',z1c(nlev+1)
        write(6,*) 'DIFF     =',z2c(nlev+1)-z1c(nlev+1)
        abort=.true.
      endif

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! quadratic splies with UK met office monotonicity constraints  !!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      zkr  = 99
      ilev = 2
      zkr(1) = 1
      zkr(nlev+1) = nlev
      kloop: do k = 2,nlev
        do jk = ilev,nlev+1
          if (z1c(jk).ge.z2c(k)) then
            ilev      = jk
            zkr(k)   = jk-1
            cycle kloop
          endif
        enddo
      enddo kloop

      zgam  = (z2c(1:nlev+1)-z1c(zkr)) / (z1c(zkr+1)-z1c(zkr))
      zgam(1)      = 0.0
      zgam(nlev+1) = 1.0
      zhdp = z1c(2:nlev+1)-z1c(1:nlev)


      h = 1/zhdp
      zarg = Qcol * h
      rhs = 0
      lower_diag = 0
      diag = 0
      upper_diag = 0

      rhs(1)=3*zarg(1)
      rhs(2:nlev) = 3*(zarg(2:nlev)*h(2:nlev) + zarg(1:nlev-1)*h(1:nlev-1))
      rhs(nlev+1)=3*zarg(nlev)

      lower_diag(1)=1
      lower_diag(2:nlev) = h(1:nlev-1)
      lower_diag(nlev+1)=1

      diag(1)=2
      diag(2:nlev) = 2*(h(2:nlev) + h(1:nlev-1))
      diag(nlev+1)=2

      upper_diag(1)=1
      upper_diag(2:nlev) = h(2:nlev)
      upper_diag(nlev+1)=0

      q_diag(1)=-upper_diag(1)/diag(1)
      rhs(1)= rhs(1)/diag(1)

      do k=2,nlev+1
        tmp_cal    =  1/(diag(k)+lower_diag(k)*q_diag(k-1))
        q_diag(k) = -upper_diag(k)*tmp_cal
        rhs(k) =  (rhs(k)-lower_diag(k)*rhs(k-1))*tmp_cal
      enddo
      do k=nlev,1,-1
        rhs(k)=rhs(k)+q_diag(k)*rhs(k+1)
      enddo

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!  monotonicity modifications  !!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      filter_code = 0
      dy(1:nlev-1) = zarg(2:nlev)-zarg(1:nlev-1)
      dy(nlev) = dy(nlev-1)

      dy = merge(zero, dy, abs(dy) < tiny )

      do k=1,nlev
        im1=MAX(1,k-1)
        im2=MAX(1,k-2)
        im3=MAX(1,k-3)
        ip1=MIN(nlev,k+1)
        t1 = merge(1,0,(zarg(k)-rhs(k))*(rhs(k)-zarg(im1)) >= 0)
        t2 = merge(1,0,dy(im2)*(rhs(k)-zarg(im1)) > 0 .AND. dy(im2)*dy(im3) > 0 &
             .AND. dy(k)*dy(ip1) > 0 .AND. dy(im2)*dy(k) < 0 )
        t3 = merge(1,0,ABS(rhs(k)-zarg(im1)) > ABS(rhs(k)-zarg(k)))

        filter_code(k) = merge(0,1,t1+t2 > 0)
        rhs(k) = (1-filter_code(k))*rhs(k)+filter_code(k)*(t3*zarg(k)+(1-t3)*zarg(im1))
        filter_code(im1) = MAX(filter_code(im1),filter_code(k))
      enddo

      rhs = merge(qmax,rhs,rhs > qmax)
      rhs = merge(zero,rhs,rhs < zero)

      za0 = rhs(1:nlev)
      za1 = -4*rhs(1:nlev) - 2*rhs(2:nlev+1) + 6*zarg
      za2 =  3*rhs(1:nlev) + 3*rhs(2:nlev+1) - 6*zarg

      dy(1:nlev) = rhs(2:nlev+1)-rhs(1:nlev)
      dy = merge(zero, dy, abs(dy) < tiny )

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! Compute the 3 quadratic spline coeffients {za0, za1, za2}				   !!
      !! knowing the quadratic spline parameters {rho_left,rho_right,zarg}		   !!
      !! Zerroukat et.al., Q.J.R. Meteorol. Soc., Vol. 128, pp. 2801-2820 (2002).   !!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      h = rhs(2:nlev+1)

      do k=1,nlev
        xm_d = merge(one,2*za2(k),abs(za2(k)) < tiny)
        xm = merge(zero,-za1(k)/xm_d, abs(za2(k)) < tiny)
        f_xm = za0(k) + za1(k)*xm + za2(k)*xm**2

        t1 = merge(1,0,ABS(za2(k)) > tiny)
        t2 = merge(1,0,xm <= zero .OR. xm >= 1)
        t3 = merge(1,0,za2(k) > zero)
        t4 = merge(1,0,za2(k) < zero)
        tm = merge(1,0,t1*((1-t2)+t3) .EQ. 2)
        tp = merge(1,0,t1*((1-t2)+(1-t3)+t4) .EQ. 3)

        peaks=0
        peaks = merge(-1,peaks,tm .EQ. 1)
        peaks = merge(+1,peaks,tp .EQ. 1)
        peaks_min = merge(f_xm,MIN(za0(k),za0(k)+za1(k)+za2(k)),tm .EQ. 1)
        peaks_max = merge(f_xm,MAX(za0(k),za0(k)+za1(k)+za2(k)),tp .EQ. 1)

        im1=MAX(1,k-1)
        im2=MAX(1,k-2)
        ip1=MIN(nlev,k+1)
        ip2=MIN(nlev,k+2)

        t1 = merge(abs(peaks),0,(dy(im2)*dy(im1) <= tiny) .OR. &
             (dy(ip1)*dy(ip2) <= tiny) .OR. (dy(im1)*dy(ip1) >= tiny) .OR. &
             (dy(im1)*float(peaks) <= tiny))

        filter_code(k) = merge(1,t1+(1-t1)*filter_code(k),(rhs(k) >= qmax) .OR. &
             (rhs(k) <= zero) .OR. (peaks_max > qmax) .OR. (peaks_min < tiny))

        if (filter_code(k) > 0) then
          level1 = rhs(k)
          level2 = (2*rhs(k)+h(k))/3
          level3 = 0.5*(rhs(k)+h(k))
          level4 = (1/3d0)*rhs(k)+2*(1/3d0)*h(k)
          level5 = h(k)

          t1 = merge(1,0,h(k) >= rhs(k))
          t2 = merge(1,0,zarg(k) <= level1 .OR.  zarg(k) >= level5)
          t3 = merge(1,0,zarg(k) >  level1 .AND. zarg(k) <  level2)
          t4 = merge(1,0,zarg(k) >  level4 .AND. zarg(k) <  level5)

          lt1 = t1*t2
          lt2 = t1*(1-t2+t3)
          lt3 = t1*(1-t2+1-t3+t4)

          za0(k) = merge(zarg(k),za0(k),lt1 .EQ. 1)
          za1(k) = merge(zero,za1(k),lt1 .EQ. 1)
          za2(k) = merge(zero,za2(k),lt1 .EQ. 1)

          za0(k) = merge(rhs(k),za0(k),lt2 .EQ. 2)
          za1(k) = merge(zero,za1(k),lt2 .EQ. 2)
          za2(k) = merge(3*(zarg(k)-rhs(k)),za2(k),lt2 .EQ. 2)

          za0(k) = merge(-2*h(k)+3*zarg(k),za0(k),lt3 .EQ. 3)
          za1(k) = merge(+6*h(k)-6*zarg(k),za1(k),lt3 .EQ. 3)
          za2(k) = merge(-3*h(k)+3*zarg(k),za2(k),lt3 .EQ. 3)

          t2 = merge(1,0,zarg(k) >= level1 .OR.  zarg(k) <= level5)
          t3 = merge(1,0,zarg(k) <  level1 .AND. zarg(k) >  level2)
          t4 = merge(1,0,zarg(k) <  level4 .AND. zarg(k) >  level5)

          lt1 = (1-t1)*t2
          lt2 = (1-t1)*(1-t2+t3)
          lt3 = (1-t1)*(1-t2+1-t3+t4)

          za0(k) = merge(zarg(k),za0(k),lt1 .EQ. 1)
          za1(k) = merge(zero,za1(k),lt1 .EQ. 1)
          za2(k) = merge(zero,za2(k),lt1 .EQ. 1)

          za0(k) = merge(rhs(k),za0(k),lt2 .EQ. 2)
          za1(k) = merge(zero,za1(k),lt2 .EQ. 2)
          za2(k) = merge(3*(zarg(k)-rhs(k)),za2(k),lt2 .EQ. 2)

          za0(k) = merge(-2*h(k)+3*zarg(k),za0(k),lt3 .EQ. 3)
          za1(k) = merge(+6*h(k)-6*zarg(k),za1(k),lt3 .EQ. 3)
          za2(k) = merge(-3*h(k)+3*zarg(k),za2(k),lt3 .EQ. 3)
        endif
      enddo

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! start iteration from top to bottom of atmosphere !!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      zv1 = 0
      do k=1,nlev
        if (zgam(k+1)>1d0) then
          WRITE(*,*) 'r not in [0:1]', zgam(k+1)
          abort=.true.
        endif
        zv2 = zv(zkr(k+1))+(za0(zkr(k+1))*zgam(k+1)+(za1(zkr(k+1))/2)*(zgam(k+1)**2)+ &
             (za2(zkr(k+1))/3)*(zgam(k+1)**3))*zhdp(zkr(k+1))
        Qdp(i,j,k,q) = (zv2 - zv1) ! / (z2c(k+1)-z2c(k) ) dont convert back to mixing ratio
        zv1 = zv2
      enddo
    enddo
  enddo
  enddo ! q loop
  if (abort) call abortmp('Bad levels in remap1.  usually CFL violatioin')
  call t_stopf('remap1')
end subroutine remap1

subroutine remap1_nofilter(Qdp,nx,qsize,dp1,dp2)
  ! remap 1 field
  ! input:  Qdp   field to be remapped (NOTE: MASS, not MIXING RATIO)
  !         dp1   layer thickness (source)
  !         dp2   layer thickness (target)
  !
  ! output: remaped Qdp, conserving mass
  !
  implicit none
  integer, intent(in) :: nx,qsize
  real (kind=real_kind), intent(inout) :: Qdp(nx,nx,nlev,qsize)
  real (kind=real_kind), intent(in) :: dp1(nx,nx,nlev),dp2(nx,nx,nlev)
  ! ========================
  ! Local Variables
  ! ========================

  real (kind=real_kind), dimension(nlev+1)    :: rhs,lower_diag,diag,upper_diag,q_diag,zgam,z1c,z2c,zv
  real (kind=real_kind), dimension(nlev)      :: h,Qcol,dy,za0,za1,za2,zarg,zhdp,dp_star,dp_np1
  real (kind=real_kind)  :: f_xm,level1,level2,level3,level4,level5, &
                            peaks_min,peaks_max,tmp_cal,xm,xm_d,zv1,zv2, &
                            zero = 0,one = 1,tiny = 1e-12,qmax = 1d50
  integer(kind=int_kind) :: zkr(nlev+1),filter_code(nlev),peaks,im1,im2,im3,ip1,ip2, &
                            lt1,lt2,lt3,t0,t1,t2,t3,t4,tm,tp,ie,i,ilev,j,jk,k,q
  logical :: abort=.false.
  call t_startf('remap1_nofilter')

#if (defined ELEMENT_OPENMP)
!$omp parallel do private(qsize,i,j,z1c,z2c,zv,k,dp_np1,dp_star,Qcol,zkr,ilev) &
!$omp    private(jk,zgam,zhdp,h,zarg,rhs,lower_diag,diag,upper_diag,q_diag,tmp_cal,filter_code) &
!$omp    private(dy,im1,im2,im3,ip1,t1,t2,t3,za0,za1,za2,xm_d,xm,f_xm,t4,tm,tp,peaks,peaks_min) &
!$omp    private(peaks_max,ip2,level1,level2,level3,level4,level5,lt1,lt2,lt3,zv1,zv2)
#endif
  do q=1,qsize
  do i=1,nx
    do j=1,nx

      z1c(1)=0 ! source grid
      z2c(1)=0 ! target grid
      do k=1,nlev
         z1c(k+1)=z1c(k)+dp1(i,j,k)
         z2c(k+1)=z2c(k)+dp2(i,j,k)
      enddo

      zv(1)=0
      do k=1,nlev
        Qcol(k)=Qdp(i,j,k,q)!  *(z1c(k+1)-z1c(k)) input is mass
        zv(k+1) = zv(k)+Qcol(k)
      enddo

      if (ABS(z2c(nlev+1)-z1c(nlev+1)).GE.0.000001) then
        write(6,*) 'SURFACE PRESSURE IMPLIED BY ADVECTION SCHEME'
        write(6,*) 'NOT CORRESPONDING TO SURFACE PRESSURE IN    '
        write(6,*) 'DATA FOR MODEL LEVELS'
        write(6,*) 'PLEVMODEL=',z2c(nlev+1)
        write(6,*) 'PLEV     =',z1c(nlev+1)
        write(6,*) 'DIFF     =',z2c(nlev+1)-z1c(nlev+1)
        abort=.true.
      endif

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! quadratic splies with UK met office monotonicity constraints  !!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      zkr  = 99
      ilev = 2
      zkr(1) = 1
      zkr(nlev+1) = nlev
      kloop: do k = 2,nlev
        do jk = ilev,nlev+1
          if (z1c(jk).ge.z2c(k)) then
            ilev      = jk
            zkr(k)   = jk-1
            cycle kloop
          endif
        enddo
      enddo kloop

      zgam  = (z2c(1:nlev+1)-z1c(zkr)) / (z1c(zkr+1)-z1c(zkr))
      zgam(1)      = 0.0
      zgam(nlev+1) = 1.0
      zhdp = z1c(2:nlev+1)-z1c(1:nlev)


      h = 1/zhdp
      zarg = Qcol * h
      rhs = 0
      lower_diag = 0
      diag = 0
      upper_diag = 0

      rhs(1)=3*zarg(1)
      rhs(2:nlev) = 3*(zarg(2:nlev)*h(2:nlev) + zarg(1:nlev-1)*h(1:nlev-1))
      rhs(nlev+1)=3*zarg(nlev)

      lower_diag(1)=1
      lower_diag(2:nlev) = h(1:nlev-1)
      lower_diag(nlev+1)=1

      diag(1)=2
      diag(2:nlev) = 2*(h(2:nlev) + h(1:nlev-1))
      diag(nlev+1)=2

      upper_diag(1)=1
      upper_diag(2:nlev) = h(2:nlev)
      upper_diag(nlev+1)=0

      q_diag(1)=-upper_diag(1)/diag(1)
      rhs(1)= rhs(1)/diag(1)

      do k=2,nlev+1
        tmp_cal    =  1/(diag(k)+lower_diag(k)*q_diag(k-1))
        q_diag(k) = -upper_diag(k)*tmp_cal
        rhs(k) =  (rhs(k)-lower_diag(k)*rhs(k-1))*tmp_cal
      enddo
      do k=nlev,1,-1
        rhs(k)=rhs(k)+q_diag(k)*rhs(k+1)
      enddo

      za0 = rhs(1:nlev)
      za1 = -4*rhs(1:nlev) - 2*rhs(2:nlev+1) + 6*zarg
      za2 =  3*rhs(1:nlev) + 3*rhs(2:nlev+1) - 6*zarg


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! start iteration from top to bottom of atmosphere !!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      zv1 = 0
      do k=1,nlev
        if (zgam(k+1)>1d0) then
          WRITE(*,*) 'r not in [0:1]', zgam(k+1)
          abort=.true.
        endif
        zv2 = zv(zkr(k+1))+(za0(zkr(k+1))*zgam(k+1)+(za1(zkr(k+1))/2)*(zgam(k+1)**2)+ &
             (za2(zkr(k+1))/3)*(zgam(k+1)**3))*zhdp(zkr(k+1))
        Qdp(i,j,k,q) = (zv2 - zv1) ! / (z2c(k+1)-z2c(k) ) dont convert back to mixing ratio
        zv1 = zv2
      enddo
    enddo
  enddo
  enddo ! q loop
  if (abort) call abortmp('Bad levels in remap1_nofilter.  usually CFL violatioin')
  call t_stopf('remap1_nofilter')
end subroutine remap1_nofilter

!=======================================================================================================!


!This uses the exact same model and reference grids and data as remap_Q, but it interpolates
!using PPM instead of splines.
subroutine remap_Q_ppm(Qdp,nx,qsize,dp1,dp2)
  ! remap 1 field
  ! input:  Qdp   field to be remapped (NOTE: MASS, not MIXING RATIO)
  !         dp1   layer thickness (source)
  !         dp2   layer thickness (target)
  !
  ! output: remaped Qdp, conserving mass
  !
  use control_mod, only        : prescribed_wind, vert_remap_q_alg
  implicit none
  integer,intent(in) :: nx,qsize
  real (kind=real_kind), intent(inout) :: Qdp(nx,nx,nlev,qsize)
  real (kind=real_kind), intent(in) :: dp1(nx,nx,nlev),dp2(nx,nx,nlev)
  ! Local Variables
  integer, parameter :: gs = 2                              !Number of cells to place in the ghost region
  real(kind=real_kind), dimension(       nlev+2 ) :: pio    !Pressure at interfaces for old grid
  real(kind=real_kind), dimension(       nlev+1 ) :: pin    !Pressure at interfaces for new grid
  real(kind=real_kind), dimension(       nlev+1 ) :: masso  !Accumulate mass up to each interface
  real(kind=real_kind), dimension(  1-gs:nlev+gs) :: ao     !Tracer value on old grid
  real(kind=real_kind), dimension(  1-gs:nlev+gs) :: dpo    !change in pressure over a cell for old grid
  real(kind=real_kind), dimension(  1-gs:nlev+gs) :: dpn    !change in pressure over a cell for old grid
  real(kind=real_kind), dimension(3,     nlev   ) :: coefs  !PPM coefficients within each cell
  real(kind=real_kind), dimension(       nlev   ) :: z1, z2
  real(kind=real_kind) :: ppmdx(10,0:nlev+1)  !grid spacings
  real(kind=real_kind) :: mymass, massn1, massn2
  integer :: i, j, k, q, kk, kid(nlev)

  call t_startf('remap_Q_ppm')
  do j = 1 , nx
    do i = 1 , nx

      pin(1)=0
      pio(1)=0
      do k=1,nlev
         dpn(k)=dp2(i,j,k)
         dpo(k)=dp1(i,j,k)
         pin(k+1)=pin(k)+dpn(k)
         pio(k+1)=pio(k)+dpo(k)
      enddo



      pio(nlev+2) = pio(nlev+1) + 1.  !This is here to allow an entire block of k threads to run in the remapping phase.
                                      !It makes sure there's an old interface value below the domain that is larger.
      pin(nlev+1) = pio(nlev+1)       !The total mass in a column does not change.
                                      !Therefore, the pressure of that mass cannot either.
      !Fill in the ghost regions with mirrored values. if vert_remap_q_alg is defined, this is of no consequence.
      do k = 1 , gs
        dpo(1   -k) = dpo(       k)
        dpo(nlev+k) = dpo(nlev+1-k)
      enddo

      !Compute remapping intervals once for all tracers. Find the old grid cell index in which the
      !k-th new cell interface resides. Then integrate from the bottom of that old cell to the new
      !interface location. In practice, the grid never deforms past one cell, so the search can be
      !simplified by this. Also, the interval of integration is usually of magnitude close to zero
      !or close to dpo because of minimial deformation.
      !Numerous tests confirmed that the bottom and top of the grids match to machine precision, so
      !I set them equal to each other.
      do k = 1 , nlev
        kk = k  !Keep from an order n^2 search operation by assuming the old cell index is close.
        !Find the index of the old grid cell in which this new cell's bottom interface resides.
        do while ( pio(kk) <= pin(k+1) )
          kk = kk + 1
        enddo
        kk = kk - 1                   !kk is now the cell index we're integrating over.
        if (kk == nlev+1) kk = nlev   !This is to keep the indices in bounds.
                                      !Top bounds match anyway, so doesn't matter what coefficients are used
        kid(k) = kk                   !Save for reuse
        z1(k) = -0.5D0                !This remapping assumes we're starting from the left interface of an old grid cell
                                      !In fact, we're usually integrating very little or almost all of the cell in question
        z2(k) = ( pin(k+1) - ( pio(kk) + pio(kk+1) ) * 0.5 ) / dpo(kk)  !PPM interpolants are normalized to an independent
                                                                        !coordinate domain [-0.5,0.5].
      enddo

      !This turned out a big optimization, remembering that only parts of the PPM algorithm depends on the data, namely the
      !limiting. So anything that depends only on the grid is pre-computed outside the tracer loop.
      ppmdx(:,:) = compute_ppm_grids( dpo )

      !From here, we loop over tracers for only those portions which depend on tracer data, which includes PPM limiting and
      !mass accumulation
      do q = 1 , qsize
        !Accumulate the old mass up to old grid cell interface locations to simplify integration
        !during remapping. Also, divide out the grid spacing so we're working with actual tracer
        !values and can conserve mass. The option for ifndef ZEROHORZ I believe is there to ensure
        !tracer consistency for an initially uniform field. I copied it from the old remap routine.
        masso(1) = 0.
        do k = 1 , nlev
          ao(k) = Qdp(i,j,k,q)
          masso(k+1) = masso(k) + ao(k) !Accumulate the old mass. This will simplify the remapping
          ao(k) = ao(k) / dpo(k)        !Divide out the old grid spacing because we want the tracer mixing ratio, not mass.
        enddo
        !Fill in ghost values. Ignored if vert_remap_q_alg == 2
        do k = 1 , gs
          ao(1   -k) = ao(       k)
          ao(nlev+k) = ao(nlev+1-k)
        enddo
        !Compute monotonic and conservative PPM reconstruction over every cell
        coefs(:,:) = compute_ppm( ao , ppmdx )
        !Compute tracer values on the new grid by integrating from the old cell bottom to the new
        !cell interface to form a new grid mass accumulation. Taking the difference between
        !accumulation at successive interfaces gives the mass inside each cell. Since Qdp is
        !supposed to hold the full mass this needs no normalization.
        massn1 = 0.
        do k = 1 , nlev
          kk = kid(k)
          massn2 = masso(kk) + integrate_parabola( coefs(:,kk) , z1(k) , z2(k) ) * dpo(kk)
          Qdp(i,j,k,q) = massn2 - massn1
          massn1 = massn2
        enddo
      enddo
    enddo
  enddo
  call t_stopf('remap_Q_ppm')
end subroutine remap_Q_ppm


!=======================================================================================================!


!THis compute grid-based coefficients from Collela & Woodward 1984.
function compute_ppm_grids( dx )   result(rslt)
  use control_mod, only: vert_remap_q_alg
  implicit none
  real(kind=real_kind), intent(in) :: dx(-1:nlev+2)  !grid spacings
  real(kind=real_kind)             :: rslt(10,0:nlev+1)  !grid spacings
  integer :: j
  integer :: indB, indE

  !Calculate grid-based coefficients for stage 1 of compute_ppm
  if (vert_remap_q_alg == 2) then
    indB = 2
    indE = nlev-1
  else
    indB = 0
    indE = nlev+1
  endif
  do j = indB , indE
    rslt( 1,j) = dx(j) / ( dx(j-1) + dx(j) + dx(j+1) )
    rslt( 2,j) = ( 2.*dx(j-1) + dx(j) ) / ( dx(j+1) + dx(j) )
    rslt( 3,j) = ( dx(j) + 2.*dx(j+1) ) / ( dx(j-1) + dx(j) )
  enddo

  !Caculate grid-based coefficients for stage 2 of compute_ppm
  if (vert_remap_q_alg == 2) then
    indB = 2
    indE = nlev-2
  else
    indB = 0
    indE = nlev
  endif
  do j = indB , indE
    rslt( 4,j) = dx(j) / ( dx(j) + dx(j+1) )
    rslt( 5,j) = 1. / sum( dx(j-1:j+2) )
    rslt( 6,j) = ( 2. * dx(j+1) * dx(j) ) / ( dx(j) + dx(j+1 ) )
    rslt( 7,j) = ( dx(j-1) + dx(j  ) ) / ( 2. * dx(j  ) + dx(j+1) )
    rslt( 8,j) = ( dx(j+2) + dx(j+1) ) / ( 2. * dx(j+1) + dx(j  ) )
    rslt( 9,j) = dx(j  ) * ( dx(j-1) + dx(j  ) ) / ( 2.*dx(j  ) +    dx(j+1) )
    rslt(10,j) = dx(j+1) * ( dx(j+1) + dx(j+2) ) / (    dx(j  ) + 2.*dx(j+1) )
  enddo
end function compute_ppm_grids

!=======================================================================================================!



!This computes a limited parabolic interpolant using a net 5-cell stencil, but the stages of computation are broken up into 3 stages
function compute_ppm( a , dx )    result(coefs)
  use control_mod, only: vert_remap_q_alg
  implicit none
  real(kind=real_kind), intent(in) :: a    (    -1:nlev+2)  !Cell-mean values
  real(kind=real_kind), intent(in) :: dx   (10,  0:nlev+1)  !grid spacings
  real(kind=real_kind) ::             coefs(0:2,   nlev  )  !PPM coefficients (for parabola)
  real(kind=real_kind) :: ai (0:nlev  )                     !fourth-order accurate, then limited interface values
  real(kind=real_kind) :: dma(0:nlev+1)                     !An expression from Collela's '84 publication
  real(kind=real_kind) :: da                                !Ditto
  ! Hold expressions based on the grid (which are cumbersome).
  real(kind=real_kind) :: dx1, dx2, dx3, dx4, dx5, dx6, dx7, dx8, dx9, dx10
  real(kind=real_kind) :: al, ar                            !Left and right interface values for cell-local limiting
  integer :: j
  integer :: indB, indE

  ! Stage 1: Compute dma for each cell, allowing a 1-cell ghost stencil below and above the domain
  if (vert_remap_q_alg == 2) then
    indB = 2
    indE = nlev-1
  else
    indB = 0
    indE = nlev+1
  endif
  do j = indB , indE
    da = dx(1,j) * ( dx(2,j) * ( a(j+1) - a(j) ) + dx(3,j) * ( a(j) - a(j-1) ) )
    dma(j) = minval( (/ abs(da) , 2. * abs( a(j) - a(j-1) ) , 2. * abs( a(j+1) - a(j) ) /) ) * sign(1.D0,da)
    if ( ( a(j+1) - a(j) ) * ( a(j) - a(j-1) ) <= 0. ) dma(j) = 0.
  enddo

  ! Stage 2: Compute ai for each cell interface in the physical domain (dimension nlev+1)
  if (vert_remap_q_alg == 2) then
    indB = 2
    indE = nlev-2
  else
    indB = 0
    indE = nlev
  endif
  do j = indB , indE
    ai(j) = a(j) + dx(4,j) * ( a(j+1) - a(j) ) + dx(5,j) * ( dx(6,j) * ( dx(7,j) - dx(8,j) ) &
         * ( a(j+1) - a(j) ) - dx(9,j) * dma(j+1) + dx(10,j) * dma(j) )
  enddo

  ! Stage 3: Compute limited PPM interpolant over each cell in the physical domain
  ! (dimension nlev) using ai on either side and ao within the cell.
  if (vert_remap_q_alg == 2) then
    indB = 3
    indE = nlev-2
  else
    indB = 1
    indE = nlev
  endif
  do j = indB , indE
    al = ai(j-1)
    ar = ai(j  )
    if ( (ar - a(j)) * (a(j) - al) <= 0. ) then
      al = a(j)
      ar = a(j)
    endif
    if ( (ar - al) * (a(j) - (al + ar)/2.) >  (ar - al)**2/6. ) al = 3.*a(j) - 2. * ar
    if ( (ar - al) * (a(j) - (al + ar)/2.) < -(ar - al)**2/6. ) ar = 3.*a(j) - 2. * al
    !Computed these coefficients from the edge values and cell mean in Maple. Assumes normalized coordinates: xi=(x-x0)/dx
    coefs(0,j) = 1.5 * a(j) - ( al + ar ) / 4.
    coefs(1,j) = ar - al
    coefs(2,j) = -6. * a(j) + 3. * ( al + ar )
  enddo

  !If we're not using a mirrored boundary condition, then make the two cells bordering the top and bottom
  !material boundaries piecewise constant. Zeroing out the first and second moments, and setting the zeroth
  !moment to the cell mean is sufficient to maintain conservation.
  if (vert_remap_q_alg == 2) then
    coefs(0,1:2) = a(1:2)
    coefs(1:2,1:2) = 0.
    coefs(0,nlev-1:nlev) = a(nlev-1:nlev)
    coefs(1:2,nlev-1:nlev) = 0.D0
  endif
end function compute_ppm

!=======================================================================================================!


!Simple function computes the definite integral of a parabola in normalized coordinates, xi=(x-x0)/dx,
!given two bounds. Make sure this gets inlined during compilation.
function integrate_parabola( a , x1 , x2 )    result(mass)
  implicit none
  real(kind=real_kind), intent(in) :: a(0:2)  !Coefficients of the parabola
  real(kind=real_kind), intent(in) :: x1      !lower domain bound for integration
  real(kind=real_kind), intent(in) :: x2      !upper domain bound for integration
  real(kind=real_kind)             :: mass
  mass = a(0) * (x2 - x1) + a(1) * (x2 ** 2 - x1 ** 2) / 0.2D1 + a(2) * (x2 ** 3 - x1 ** 3) / 0.3D1
end function integrate_parabola


!=============================================================================================!



end module vertremap_mod

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! End GPU remap module    !!
!! by Rick Archibald, 2010  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!=======================================================================================================!




module prim_advection_mod
!
! two formulations.  both are conservative
! u grad Q formulation:
!
!    d/dt[ Q] +  U grad Q  +  eta_dot dp/dn dQ/dp  = 0
!                            ( eta_dot dQ/dn )
!
!    d/dt[ dp/dn ] = div( dp/dn U ) + d/dn ( eta_dot dp/dn )
!
! total divergence formulation:
!    d/dt[dp/dn Q] +  div( U dp/dn Q ) + d/dn ( eta_dot dp/dn Q ) = 0
!
! for convience, rewrite this as dp Q:  (since dn does not depend on time or the horizonal):
! equation is now:
!    d/dt[dp Q] +  div( U dp Q ) + d( eta_dot_dpdn Q ) = 0
!
!
  use kinds, only              : real_kind
  use dimensions_mod, only     : nlev, nlevp, np, qsize, ntrac, nc, nep
  use physical_constants, only : rgas, Rwater_vapor, kappa, g, rearth, rrearth, cp
  use derivative_mod, only     : gradient, vorticity, gradient_wk, derivative_t, divergence, &
                                 gradient_sphere, divergence_sphere
  use element_mod, only        : element_t
  use fvm_control_volume_mod, only        : fvm_struct
  use spelt_mod, only          : spelt_struct
  use filter_mod, only         : filter_t, filter_P
  use hybvcoord_mod, only      : hvcoord_t
  use time_mod, only           : TimeLevel_t, smooth, TimeLevel_Qdp
  use prim_si_mod, only        : preq_pressure
  use diffusion_mod, only      : scalar_diffusion, diffusion_init
  use control_mod, only        : integration, test_case, filter_freq_advection,  hypervis_order, &
        statefreq, moisture, TRACERADV_TOTAL_DIVERGENCE, TRACERADV_UGRADQ, &
        prescribed_wind, nu_q, nu_p, limiter_option, hypervis_subcycle_q, rsplit
  use edge_mod, only           : EdgeBuffer_t, edgevpack, edgerotate, edgevunpack, initedgebuffer, edgevunpackmin
  use hybrid_mod, only         : hybrid_t
  use bndry_mod, only          : bndry_exchangev
  use viscosity_mod, only      : biharmonic_wk_scalar, biharmonic_wk_scalar_minmax, neighbor_minmax
  use perf_mod, only           : t_startf, t_stopf, t_barrierf ! _EXTERNAL
  use parallel_mod, only : parallel_t, abortmp, status, srequest, rrequest,  mpireal_t
  implicit none

  private
  save

  public :: Prim_Advec_Init, Prim_Advec_Tracers_remap_rk2
  public :: prim_advec_tracers_fvm, prim_advec_tracers_spelt
  public :: vertical_remap

  type (EdgeBuffer_t) :: edgeAdv, edgeAdvQ3, edgeAdv_p1, edgeAdvQ2, edgeAdv1

  integer,parameter :: DSSeta = 1
  integer,parameter :: DSSomega = 2
  integer,parameter :: DSSdiv_vdp_ave = 3
  integer,parameter :: DSSno_var = -1

  real(kind=real_kind), allocatable :: qmin(:,:,:), qmax(:,:,:)


contains

  subroutine Prim_Advec_Init()
    use dimensions_mod, only : nlev, qsize, nelemd

    ! Shared buffer pointers.
    ! Using "=> null()" in a subroutine is usually bad, because it makes
    ! the variable have an implicit "save", and therefore shared between
    ! threads. But in this case we want shared pointers.
    real(kind=real_kind), pointer :: buf_ptr(:) => null()
    real(kind=real_kind), pointer :: receive_ptr(:) => null()

    ! this might be called with qsize=0
    ! allocate largest one first
    ! Currently this is never freed. If it was, only this first one should
    ! be freed, as only it knows the true size of the buffer.
    call initEdgeBuffer(edgeAdvQ3,max(nlev,qsize*nlev*3), buf_ptr, receive_ptr)  ! Qtens,Qmin, Qmax

    ! remaining edge buffers can share %buf and %receive with edgeAdvQ3
    ! (This is done through the optional 1D pointer arguments.)
    call initEdgeBuffer(edgeAdv1,nlev,buf_ptr,receive_ptr)
    call initEdgeBuffer(edgeAdv,qsize*nlev,buf_ptr,receive_ptr)
    call initEdgeBuffer(edgeAdv_p1,qsize*nlev + nlev,buf_ptr,receive_ptr)
    call initEdgeBuffer(edgeAdvQ2,qsize*nlev*2,buf_ptr,receive_ptr)  ! Qtens,Qmin, Qmax

    ! Don't actually want these saved, if this is ever called twice.
    nullify(buf_ptr)
    nullify(receive_ptr)

    ! this static array is shared by all threads, so dimension for all threads (nelemd), not nets:nete:
    allocate (qmin(nlev,qsize,nelemd))
    allocate (qmax(nlev,qsize,nelemd))

  end subroutine Prim_Advec_Init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SPELT driver
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine Prim_Advec_Tracers_spelt(elem, spelt, deriv,hvcoord,hybrid,&
        dt,tl,nets,nete)
    use perf_mod, only : t_startf, t_stopf            ! _EXTERNAL
    use spelt_mod, only : spelt_run, spelt_runair, edgeveloc, spelt_mcgregordss, spelt_rkdss
    use derivative_mod, only : interpolate_gll2spelt_points
    use vertremap_mod, only: remap1_nofilter ! _EXTERNAL (actually INTERNAL)

    implicit none
    type (element_t), intent(inout)   :: elem(:)
    type (spelt_struct), intent(inout)  :: spelt(:)
    type (derivative_t), intent(in)   :: deriv
    type (hvcoord_t)                  :: hvcoord
    type (hybrid_t),     intent(in):: hybrid
    type (TimeLevel_t)                :: tl

    real(kind=real_kind) , intent(in) :: dt
    integer,intent(in)                :: nets,nete

    real (kind=real_kind), dimension(np,np,nlev)    :: dp_star
    real (kind=real_kind), dimension(np,np,nlev)    :: dp

    integer :: np1,ie,k, i, j


    call t_barrierf('sync_prim_advec_tracers_spelt', hybrid%par%comm)
    call t_startf('prim_advec_tracers_spelt')
    np1 = tl%np1


    ! interpolate t+1 velocity from reference levels to lagrangian levels
    ! For rsplit=0, we need to first compute lagrangian levels based on vertical velocity
    ! which requires we first DSS mean vertical velocity from dynamics
    !
    if (rsplit==0) then
       do ie=nets,nete
          do k=1,nlev
             elem(ie)%derived%eta_dot_dpdn(:,:,k) = elem(ie)%spheremp(:,:)*elem(ie)%derived%eta_dot_dpdn(:,:,k)
          enddo
          call edgeVpack(edgeAdv1,elem(ie)%derived%eta_dot_dpdn(:,:,1:nlev),nlev,0,elem(ie)%desc)
       enddo
       call bndry_exchangeV(hybrid,edgeAdv1)
       do ie=nets,nete
          call edgeVunpack(edgeAdv1,elem(ie)%derived%eta_dot_dpdn(:,:,1:nlev),nlev,0,elem(ie)%desc)
          do k=1,nlev
             elem(ie)%derived%eta_dot_dpdn(:,:,k)=elem(ie)%derived%eta_dot_dpdn(:,:,k)*elem(ie)%rspheremp(:,:)
          enddo
          ! SET VERTICAL VELOCITY TO ZERO FOR DEBUGGING
          elem(ie)%derived%eta_dot_dpdn(:,:,:)=0

          ! elem%state%u(np1)  = velocity at time t+1 on reference levels
          ! elem%derived%vstar = velocity at t+1 on floating levels (computed below) using eta_dot_dpdn
!           call remap_UV_ref2lagrange(np1,dt,elem,hvcoord,ie)
          do k=1,nlev
             dp(:,:,k) = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
                  ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(:,:,np1)
             dp_star(:,:,k) = dp(:,:,k) + dt*(elem(ie)%derived%eta_dot_dpdn(:,:,k+1) -&
                  elem(ie)%derived%eta_dot_dpdn(:,:,k))
          enddo
          elem(ie)%derived%vstar=elem(ie)%state%v(:,:,:,:,np1)
          call remap1_nofilter(elem(ie)%derived%vstar,np,1,dp,dp_star)
          !take the average on level, should be improved later, because we know the SE velocity at t+1/2
          spelt(ie)%vn12=(spelt(ie)%vn0+elem(ie)%derived%vstar)/2.0D0
       end do
    else
       ! for rsplit>0:  dynamics is also vertically lagrangian, so we do not need to
       ! remap the velocities
       stop 'FVM need to use lagrangian winds here'
    endif


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! 2D advection step
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!------------------------------------------------------------------------------------

    call t_startf('spelt_depalg')
!     call spelt_mcgregordss(elem,spelt,nets,nete, hybrid, deriv, dt, 3)
    call spelt_rkdss(elem,spelt,nets,nete, hybrid, deriv, dt, 3)
    call t_stopf('spelt_depalg')

    ! ! end mcgregordss
    ! spelt departure calcluation should use vstar.
    ! from c(n0) compute c(np1):
!     call spelt_run(elem,spelt,hybrid,deriv,dt,tl,nets,nete)
    call spelt_runair(elem,spelt,hybrid,deriv,dt,tl,nets,nete)

    call t_stopf('prim_advec_tracers_spelt')
  end subroutine Prim_Advec_Tracers_spelt

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! fvm driver
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine Prim_Advec_Tracers_fvm(elem, fvm, deriv,hvcoord,hybrid,&
        dt,tl,nets,nete)
    use perf_mod, only : t_startf, t_stopf            ! _EXTERNAL
    use vertremap_mod, only: remap1_nofilter  ! _EXTERNAL (actually INTERNAL)
!    use fvm_mod, only : cslam_run, cslam_runairdensity, edgeveloc, fvm_mcgregor, fvm_mcgregordss
    use fvm_mod, only : cslam_runairdensity, edgeveloc, fvm_mcgregor, fvm_mcgregordss, fvm_rkdss

    implicit none
    type (element_t), intent(inout)   :: elem(:)
    type (fvm_struct), intent(inout)   :: fvm(:)
    type (derivative_t), intent(in)   :: deriv
    type (hvcoord_t)                  :: hvcoord
    type (hybrid_t),     intent(in):: hybrid
    type (TimeLevel_t)                :: tl

    real(kind=real_kind) , intent(in) :: dt
    integer,intent(in)                :: nets,nete


    real (kind=real_kind), dimension(np,np,nlev)    :: dp_star
    real (kind=real_kind), dimension(np,np,nlev)    :: dp

    integer :: np1,ie,k

    real (kind=real_kind)  :: vstar(np,np,2)
    real (kind=real_kind)  :: vhat(np,np,2)
    real (kind=real_kind), dimension(np, np) :: v1, v2


    call t_barrierf('sync_prim_advec_tracers_fvm', hybrid%par%comm)
    call t_startf('prim_advec_tracers_fvm')
    np1 = tl%np1

    ! interpolate t+1 velocity from reference levels to lagrangian levels
    ! For rsplit=0, we need to first compute lagrangian levels based on vertical velocity
    ! which requires we first DSS mean vertical velocity from dynamics
    !
    if (rsplit==0) then
       do ie=nets,nete
          do k=1,nlev
             elem(ie)%derived%eta_dot_dpdn(:,:,k) = elem(ie)%spheremp(:,:)*elem(ie)%derived%eta_dot_dpdn(:,:,k)
          enddo
          call edgeVpack(edgeAdv1,elem(ie)%derived%eta_dot_dpdn(:,:,1:nlev),nlev,0,elem(ie)%desc)
       enddo
       call bndry_exchangeV(hybrid,edgeAdv1)
       do ie=nets,nete
          call edgeVunpack(edgeAdv1,elem(ie)%derived%eta_dot_dpdn(:,:,1:nlev),nlev,0,elem(ie)%desc)
          do k=1,nlev
             elem(ie)%derived%eta_dot_dpdn(:,:,k)=elem(ie)%derived%eta_dot_dpdn(:,:,k)*elem(ie)%rspheremp(:,:)
          enddo

          ! SET VERTICAL VELOCITY TO ZERO FOR DEBUGGING
          !        elem(ie)%derived%eta_dot_dpdn(:,:,:)=0
          ! elem%state%u(np1)  = velocity at time t+1 on reference levels
          ! elem%derived%vstar = velocity at t+1 on floating levels (computed below)
!           call remap_UV_ref2lagrange(np1,dt,elem,hvcoord,ie)
          do k=1,nlev
             dp(:,:,k) = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
                  ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(:,:,np1)
             dp_star(:,:,k) = dp(:,:,k) + dt*(elem(ie)%derived%eta_dot_dpdn(:,:,k+1) -&
                  elem(ie)%derived%eta_dot_dpdn(:,:,k))
          enddo
          elem(ie)%derived%vstar=elem(ie)%state%v(:,:,:,:,np1)
          call remap1_nofilter(elem(ie)%derived%vstar,np,1,dp,dp_star)
       end do
    else
       ! for rsplit>0:  dynamics is also vertically lagrangian, so we do not need to
       ! remap the velocities
       stop 'FVM need to use lagrangian winds here'
    endif


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! 2D advection step
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!------------------------------------------------------------------------------------
    call t_startf('fvm_depalg')

!     call fvm_mcgregordss(elem,fvm,nets,nete, hybrid, deriv, dt, 3)
    call fvm_rkdss(elem,fvm,nets,nete, hybrid, deriv, dt, 3)
    call t_stopf('fvm_depalg')

!------------------------------------------------------------------------------------

    ! fvm departure calcluation should use vstar.
    ! from c(n0) compute c(np1):
    call cslam_runairdensity(elem,fvm,hybrid,deriv,dt,tl,nets,nete)

    call t_stopf('prim_advec_tracers_fvm')
  end subroutine Prim_Advec_Tracers_fvm



!=================================================================================================!

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
! forward-in-time 2 level vertically lagrangian step
!  this code takes a lagrangian step in the horizontal
! (complete with DSS), and then applies a vertical remap
!
! This routine may use dynamics fields at timelevel np1
! In addition, other fields are required, which have to be
! explicitly saved by the dynamics:  (in elem(ie)%derived struct)
!
! Fields required from dynamics: (in
!    omega_p   it will be DSS'd here, for later use by CAM physics
!              we DSS omega here because it can be done for "free"
!    Consistent mass/tracer-mass advection (used if subcycling turned on)
!       dp()   dp at timelevel n0
!       vn0()  mean flux  < U dp > going from n0 to np1
!
! 3 stage
!    Euler step from t     -> t+.5
!    Euler step from t+.5  -> t+1.0
!    Euler step from t+1.0 -> t+1.5
!    u(t) = u(t)/3 + u(t+2)*2/3
!
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
  subroutine Prim_Advec_Tracers_remap_rk2( elem , deriv , hvcoord , flt , hybrid , dt , tl , nets , nete )
    use perf_mod      , only : t_startf, t_stopf            ! _EXTERNAL
    use derivative_mod, only : divergence_sphere
    use control_mod   , only : vert_remap_q_alg, qsplit
    implicit none
    type (element_t)     , intent(inout) :: elem(:)
    type (derivative_t)  , intent(in   ) :: deriv
    type (hvcoord_t)     , intent(in   ) :: hvcoord
    type (filter_t)      , intent(in   ) :: flt
    type (hybrid_t)      , intent(in   ) :: hybrid
    real(kind=real_kind) , intent(in   ) :: dt
    type (TimeLevel_t)   , intent(inout) :: tl
    integer              , intent(in   ) :: nets
    integer              , intent(in   ) :: nete

    real (kind=real_kind), dimension(np,np,2     ) :: gradQ
    real (kind=real_kind), dimension(np,np  ,nlev) :: dp_star
    real (kind=real_kind), dimension(np,np  ,nlev) :: dp_np1
    integer :: i,j,k,l,ie,q,nmin
    integer :: nfilt,rkstage,rhs_multiplier
    integer :: n0_qdp, np1_qdp

    call t_barrierf('sync_prim_advec_tracers_remap_k2', hybrid%par%comm)
    call t_startf('prim_advec_tracers_remap_rk2')
    call TimeLevel_Qdp( tl, qsplit, n0_qdp, np1_qdp) !time levels for qdp are not the same
    rkstage = 3 !   3 stage RKSSP scheme, with optimal SSP CFL

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! RK2 2D advection step
    ! note: stage 3 we take the oppertunity to DSS omega
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! use these for consistent advection (preserve Q=1)
    ! derived%vdp_ave        =  mean horiz. flux:   U*dp
    ! derived%eta_dot_dpdn    =  mean vertical velocity (used for remap)
    ! derived%omega_p         =  advection code will DSS this for the physics, but otherwise
    !                            it is not needed
    ! Also: save a copy of div(U dp) in derived%div(:,:,:,1), which will be DSS'd
    !       and a DSS'ed version stored in derived%div(:,:,:,2)
    do ie=nets,nete
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k, gradQ)
#endif
      do k=1,nlev
        ! div( U dp Q),
        gradQ(:,:,1)=elem(ie)%derived%vn0(:,:,1,k)
        gradQ(:,:,2)=elem(ie)%derived%vn0(:,:,2,k)
        elem(ie)%derived%divdp(:,:,k) = divergence_sphere(gradQ,deriv,elem(ie))
      enddo
      elem(ie)%derived%divdp_proj(:,:,:) = elem(ie)%derived%divdp(:,:,:)
    enddo

    !rhs_multiplier is for obtaining dp_tracers at each stage:
    !dp_tracers(stage) = dp - rhs_multiplier*dt*divdp_proj
    rhs_multiplier = 0
    call euler_step( np1_qdp , n0_qdp  , dt/2 , elem , hvcoord , hybrid , deriv , nets , nete , DSSdiv_vdp_ave , rhs_multiplier )

    rhs_multiplier = 1
    call euler_step( np1_qdp , np1_qdp , dt/2 , elem , hvcoord , hybrid , deriv , nets , nete , DSSeta         , rhs_multiplier )

    rhs_multiplier = 2
    call euler_step( np1_qdp , np1_qdp , dt/2 , elem , hvcoord , hybrid , deriv , nets , nete , DSSomega       , rhs_multiplier )

    !to finish the 2D advection step, we need to average the t and t+2 results to get a second order estimate for t+1.
    call qdp_time_avg( elem , rkstage , n0_qdp , np1_qdp , limiter_option , nu_p , nets , nete )

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  Dissipation
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if ( limiter_option == 8 .or. nu_p > 0 ) then
      ! dissipation was applied in RHS.
    else
      call advance_hypervis_scalar(edgeadv,elem,hvcoord,hybrid,deriv,tl%np1,np1_qdp,nets,nete,dt)
    endif

    call t_stopf('prim_advec_tracers_remap_rk2')
  end subroutine prim_advec_tracers_remap_rk2

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

  subroutine qdp_time_avg( elem , rkstage , n0_qdp , np1_qdp , limiter_option , nu_p , nets , nete )
#ifdef _ACCEL
    use cuda_mod, only: qdp_time_avg_cuda
#endif
    implicit none
    type(element_t)     , intent(inout) :: elem(:)
    integer             , intent(in   ) :: rkstage , n0_qdp , np1_qdp , nets , nete , limiter_option
    real(kind=real_kind), intent(in   ) :: nu_p
    integer :: ie
#ifdef _ACCEL
    call qdp_time_avg_cuda( elem , rkstage , n0_qdp , np1_qdp , limiter_option , nu_p , nets , nete )
    return
#endif
    do ie=nets,nete
      elem(ie)%state%Qdp(:,:,:,1:qsize,np1_qdp) =               &
                   ( elem(ie)%state%Qdp(:,:,:,1:qsize,n0_qdp) + &
                     (rkstage-1)*elem(ie)%state%Qdp(:,:,:,1:qsize,np1_qdp) ) / rkstage
    enddo
  end subroutine qdp_time_avg

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
subroutine my_edgeVunpack_euler(edge_nlyr, edge_nbuf, &
    desc_getmapP, v, my_vlyr,kptr, my_swest, my_max_corner_elem, &
    swest_buf, seast_buf, neast_buf, nwest_buf, &
    edge_buf_is_1, edge_buf_is_2, edge_buf_is_3, edge_buf_is_4, &
    edge_buf_iw_1, edge_buf_iw_2, edge_buf_iw_3, edge_buf_iw_4, &
    edge_buf_ie_1, edge_buf_ie_2, edge_buf_ie_3, edge_buf_ie_4, &
    edge_buf_in_1, edge_buf_in_2, edge_buf_in_3, edge_buf_in_4)

  integer, intent(in) :: edge_nlyr, edge_nbuf, my_swest, my_max_corner_elem
  real(kind=8), dimension(constLev), intent(in) :: swest_buf
  real(kind=8), dimension(constLev), intent(in) :: seast_buf
  real(kind=8), dimension(constLev), intent(in) :: neast_buf
  real(kind=8), dimension(constLev), intent(in) :: nwest_buf
  real(kind=8), dimension(constLev), intent(in) :: edge_buf_is_1, edge_buf_is_2, edge_buf_is_3, edge_buf_is_4
  real(kind=8), dimension(constLev), intent(in) :: edge_buf_iw_1, edge_buf_iw_2, edge_buf_iw_3, edge_buf_iw_4
  real(kind=8), dimension(constLev), intent(in) :: edge_buf_ie_1, edge_buf_ie_2, edge_buf_ie_3, edge_buf_ie_4
  real(kind=8), dimension(constLev), intent(in) :: edge_buf_in_1, edge_buf_in_2, edge_buf_in_3, edge_buf_in_4
  integer, dimension(:), intent(in) :: desc_getmapP

  integer,               intent(in)  :: my_vlyr
  real (kind=8), intent(inout) :: v(4,4,constLev)
  integer,               intent(in)  :: kptr

  ! Local
  logical, parameter :: UseUnroll = .TRUE.
  integer :: i,k,ll

  do k=1,my_vlyr
    v(1  ,1  ,k) = v(1  ,1  ,k)+edge_buf_is_1(k)
    v(1+1,1  ,k) = v(1+1,1  ,k)+edge_buf_is_2(k)
    v(3  ,1  ,k) = v(3  ,1  ,k)+edge_buf_is_3(k)
    v(3+1,1  ,k) = v(3+1,1  ,k)+edge_buf_is_4(k)

    v(1  ,1  ,k) = v(1  ,1  ,k)+edge_buf_iw_1(k)
    v(1  ,1+1,k) = v(1  ,1+1,k)+edge_buf_iw_2(k)
    v(1  ,3  ,k) = v(1  ,3  ,k)+edge_buf_iw_3(k)
    v(1  ,3+1,k) = v(1  ,3+1,k)+edge_buf_iw_4(k)

    v(np ,1  ,k) = v(np ,1  ,k)+edge_buf_ie_1(k)
    v(np ,1+1,k) = v(np ,1+1,k)+edge_buf_ie_2(k)
    v(np ,3  ,k) = v(np ,3  ,k)+edge_buf_ie_3(k)
    v(np ,3+1,k) = v(np ,3+1,k)+edge_buf_ie_4(k)

    v(1  ,np ,k) = v(1  ,np ,k)+edge_buf_in_1(k)
    v(1+1,np ,k) = v(1+1,np ,k)+edge_buf_in_2(k)
    v(3  ,np ,k) = v(3  ,np ,k)+edge_buf_in_3(k)
    v(3+1,np ,k) = v(3+1,np ,k)+edge_buf_in_4(k)
  end do

! swest
  do ll=my_swest,my_swest+my_max_corner_elem-1
    !write(*,*), 'swest loop, ll = ', ll, 'desc_getmapP(ll) = ', desc_getmapP(ll) ! 5
      if(desc_getmapP(ll) /= -1) then
          do k=1,my_vlyr
              v(1  ,1 ,k)=v(1 ,1 ,k)+swest_buf(k)
          enddo
      endif
  end do

! SEAST
  do ll=my_swest+my_max_corner_elem,my_swest+2*my_max_corner_elem-1
    !write(*,*), 'seast loop, ll = ', ll, 'desc_getmapP(ll) = ', desc_getmapP(ll) !6
      if(desc_getmapP(ll) /= -1) then
          do k=1,my_vlyr
              v(np ,1 ,k)=v(np,1 ,k)+seast_buf(k)
          enddo
      endif
  end do

! NEAST
  do ll=my_swest+3*my_max_corner_elem,my_swest+4*my_max_corner_elem-1
    !write(*,*), 'neast loop, ll = ', ll, 'desc_getmapP(ll) = ', desc_getmapP(ll) !8
      if(desc_getmapP(ll) /= -1) then
          do k=1,my_vlyr
              v(np ,np,k)=v(np,np,k)+ neast_buf(k)
          enddo
      endif
  end do

! NWEST
  do ll=my_swest+2*my_max_corner_elem,my_swest+3*my_max_corner_elem-1
    !write(*,*), 'nwest loop, ll = ', ll, 'desc_getmapP(ll) = ', desc_getmapP(ll) !7
      if(desc_getmapP(ll) /= -1) then
          do k=1,my_vlyr
              v(1  ,np,k)=v(1 ,np,k)+nwest_buf(k)
          enddo
      endif
  end do

end subroutine

subroutine my_unpack_acc(nets, nete, edge_nlyr, edge_nbuf, &
                        edge_buf, my_south, my_east, my_north, my_west, my_nlev, &
                        my_swest, my_max_corner_elem, my_elem, DSSopt, np1_qdp, calLen, calList)
    implicit none

    integer, intent(in) :: calLen
    integer, intent(in), dimension(calLen) :: calList
    integer, intent(in) ::  edge_nlyr, edge_nbuf, my_nlev, my_south, my_east, my_north, my_west, my_swest, my_max_corner_elem
    integer, intent(in) :: nets, nete, DSSopt, np1_qdp
    real(kind=8), dimension(edge_nlyr, edge_nbuf), intent(in) :: edge_buf

    type (element_t), intent(inout) :: my_elem(nets:nete)

    real(kind=8), dimension(4,4,constLev) :: elem_state_Qdp
    pointer(elem_state_Qdp_ptr, elem_state_Qdp)

    real(kind=8), dimension(4,4) :: elem_rspheremp
    pointer(elem_rspheremp_ptr, elem_rspheremp)

    real(kind=8), dimension(4,4,constLev) :: elem_derived_eta_dot_dpdn
    pointer(elem_derived_eta_dot_dpdn_ptr, elem_derived_eta_dot_dpdn)

    real(kind=8), dimension(4,4,constLev) :: elem_derived_omega_p
    pointer(elem_derived_omega_p_ptr, elem_derived_omega_p)

    real(kind=8), dimension(4,4,constLev) :: elem_derived_divdp_proj
    pointer(elem_derived_divdp_proj_ptr, elem_derived_divdp_proj)

    integer, dimension(my_swest+4*my_max_corner_elem-1) :: getmapP
    pointer(getmapP_ptr, getmapP)

    integer :: is, ie, iw, in
    pointer(is_ptr, is)
    pointer(ie_ptr, ie)
    pointer(iw_ptr, iw)
    pointer(in_ptr, in)

    real(kind=8), dimension(constLev) :: swest_buf, seast_buf, neast_buf, nwest_buf
    pointer(swest_buf_ptr, swest_buf)
    pointer(seast_buf_ptr, seast_buf)
    pointer(neast_buf_ptr, neast_buf)
    pointer(nwest_buf_ptr, nwest_buf)

    real(kind=8), dimension(constLev) :: edge_buf_is_1, edge_buf_is_2, edge_buf_is_3, edge_buf_is_4
    real(kind=8), dimension(constLev) :: edge_buf_iw_1, edge_buf_iw_2, edge_buf_iw_3, edge_buf_iw_4
    real(kind=8), dimension(constLev) :: edge_buf_ie_1, edge_buf_ie_2, edge_buf_ie_3, edge_buf_ie_4
    real(kind=8), dimension(constLev) :: edge_buf_in_1, edge_buf_in_2, edge_buf_in_3, edge_buf_in_4
    pointer(edge_buf_is_1_ptr, edge_buf_is_1)
    pointer(edge_buf_is_2_ptr, edge_buf_is_2)
    pointer(edge_buf_is_3_ptr, edge_buf_is_3)
    pointer(edge_buf_is_4_ptr, edge_buf_is_4)
    pointer(edge_buf_iw_1_ptr, edge_buf_iw_1)
    pointer(edge_buf_iw_2_ptr, edge_buf_iw_2)
    pointer(edge_buf_iw_3_ptr, edge_buf_iw_3)
    pointer(edge_buf_iw_4_ptr, edge_buf_iw_4)
    pointer(edge_buf_ie_1_ptr, edge_buf_ie_1)
    pointer(edge_buf_ie_2_ptr, edge_buf_ie_2)
    pointer(edge_buf_ie_3_ptr, edge_buf_ie_3)
    pointer(edge_buf_ie_4_ptr, edge_buf_ie_4)
    pointer(edge_buf_in_1_ptr, edge_buf_in_1)
    pointer(edge_buf_in_2_ptr, edge_buf_in_2)
    pointer(edge_buf_in_3_ptr, edge_buf_in_3)
    pointer(edge_buf_in_4_ptr, edge_buf_in_4)

    !local
    integer :: iee, k, q, i
    integer(kind=8) :: count_start, count_stop, count_max, count_rate

    !$ACC PARALLEL LOOP collapse(2)
    do iee = nets , nete
      do q = 1, qsize
        !iee = calList(i)
        elem_state_Qdp_ptr = loc(my_elem(iee)%state%Qdp(:,:,:,q,np1_qdp))
        elem_rspheremp_ptr = loc(my_elem(iee)%rspheremp)
        getmapP_ptr        = loc(my_elem(iee)%desc%getmapP)
        is_ptr             = loc(my_elem(iee)%desc%getmapP(my_south))
        ie_ptr             = loc(my_elem(iee)%desc%getmapP(my_east))
        iw_ptr             = loc(my_elem(iee)%desc%getmapP(my_west))
        in_ptr             = loc(my_elem(iee)%desc%getmapP(my_north))
        edge_buf_is_1_ptr  = loc(edge_buf((q-1)*nlev+1, is + 1))
        edge_buf_is_2_ptr  = loc(edge_buf((q-1)*nlev+1, is + 2))
        edge_buf_is_3_ptr  = loc(edge_buf((q-1)*nlev+1, is + 3))
        edge_buf_is_4_ptr  = loc(edge_buf((q-1)*nlev+1, is + 4))
        edge_buf_iw_1_ptr  = loc(edge_buf((q-1)*nlev+1, iw + 1))
        edge_buf_iw_2_ptr  = loc(edge_buf((q-1)*nlev+1, iw + 2))
        edge_buf_iw_3_ptr  = loc(edge_buf((q-1)*nlev+1, iw + 3))
        edge_buf_iw_4_ptr  = loc(edge_buf((q-1)*nlev+1, iw + 4))
        edge_buf_ie_1_ptr  = loc(edge_buf((q-1)*nlev+1, ie + 1))
        edge_buf_ie_2_ptr  = loc(edge_buf((q-1)*nlev+1, ie + 2))
        edge_buf_ie_3_ptr  = loc(edge_buf((q-1)*nlev+1, ie + 3))
        edge_buf_ie_4_ptr  = loc(edge_buf((q-1)*nlev+1, ie + 4))
        edge_buf_in_1_ptr  = loc(edge_buf((q-1)*nlev+1, in + 1))
        edge_buf_in_2_ptr  = loc(edge_buf((q-1)*nlev+1, in + 2))
        edge_buf_in_3_ptr  = loc(edge_buf((q-1)*nlev+1, in + 3))
        edge_buf_in_4_ptr  = loc(edge_buf((q-1)*nlev+1, in + 4))
        swest_buf_ptr      = loc(edge_buf((q-1)*nlev+1, getmapP(5) + 1))
        seast_buf_ptr      = loc(edge_buf((q-1)*nlev+1, getmapP(6) + 1))
        neast_buf_ptr      = loc(edge_buf((q-1)*nlev+1, getmapP(8) + 1))
        nwest_buf_ptr      = loc(edge_buf((q-1)*nlev+1, getmapP(7) + 1))

        !$ACC DATA copy(elem_state_Qdp) copyin(getmapP,is,ie,iw,in,edge_buf_is_1,edge_buf_is_2,edge_buf_is_3,edge_buf_is_4,edge_buf_iw_1,edge_buf_iw_2,edge_buf_iw_3,edge_buf_iw_4,edge_buf_ie_1,edge_buf_ie_2,edge_buf_ie_3,edge_buf_ie_4,edge_buf_in_1,edge_buf_in_2,edge_buf_in_3,edge_buf_in_4,swest_buf,seast_buf,neast_buf,nwest_buf, elem_rspheremp)
        call my_edgeVunpack_euler(edge_nlyr, edge_nbuf, getmapP(:), &
          elem_state_Qdp(:,:,:), nlev, (q-1)*nlev, my_swest, my_max_corner_elem, &
          swest_buf, seast_buf, neast_buf, nwest_buf, &
          edge_buf_is_1, edge_buf_is_2, edge_buf_is_3, edge_buf_is_4, &
          edge_buf_iw_1, edge_buf_iw_2, edge_buf_iw_3, edge_buf_iw_4, &
          edge_buf_ie_1, edge_buf_ie_2, edge_buf_ie_3, edge_buf_ie_4, &
          edge_buf_in_1, edge_buf_in_2, edge_buf_in_3, edge_buf_in_4)

        do k = 1 , nlev
          elem_state_Qdp(:,:,k) = elem_rspheremp(:,:) * elem_state_Qdp(:,:,k)
        enddo

        !$ACC END DATA
      end do
    end do
    !$ACC END PARALLEL LOOP

    if (DSSopt == DSSeta) then
      !$ACC PARALLEL LOOP
      do iee = nets, nete
        !iee = calList(i) 
        elem_derived_eta_dot_dpdn_ptr = loc(my_elem(iee)%derived%eta_dot_dpdn)
        elem_rspheremp_ptr            = loc(my_elem(iee)%rspheremp)
        getmapP_ptr                   = loc(my_elem(iee)%desc%getmapP)
        is_ptr                        = loc(my_elem(iee)%desc%getmapP(my_south))
        ie_ptr                        = loc(my_elem(iee)%desc%getmapP(my_east))
        iw_ptr                        = loc(my_elem(iee)%desc%getmapP(my_west))
        in_ptr                        = loc(my_elem(iee)%desc%getmapP(my_north))
        edge_buf_is_1_ptr             = loc(edge_buf(qsize*nlev+1, is + 1))
        edge_buf_is_2_ptr             = loc(edge_buf(qsize*nlev+1, is + 2))
        edge_buf_is_3_ptr             = loc(edge_buf(qsize*nlev+1, is + 3))
        edge_buf_is_4_ptr             = loc(edge_buf(qsize*nlev+1, is + 4))
        edge_buf_iw_1_ptr             = loc(edge_buf(qsize*nlev+1, iw + 1))
        edge_buf_iw_2_ptr             = loc(edge_buf(qsize*nlev+1, iw + 2))
        edge_buf_iw_3_ptr             = loc(edge_buf(qsize*nlev+1, iw + 3))
        edge_buf_iw_4_ptr             = loc(edge_buf(qsize*nlev+1, iw + 4))
        edge_buf_ie_1_ptr             = loc(edge_buf(qsize*nlev+1, ie + 1))
        edge_buf_ie_2_ptr             = loc(edge_buf(qsize*nlev+1, ie + 2))
        edge_buf_ie_3_ptr             = loc(edge_buf(qsize*nlev+1, ie + 3))
        edge_buf_ie_4_ptr             = loc(edge_buf(qsize*nlev+1, ie + 4))
        edge_buf_in_1_ptr             = loc(edge_buf(qsize*nlev+1, in + 1))
        edge_buf_in_2_ptr             = loc(edge_buf(qsize*nlev+1, in + 2))
        edge_buf_in_3_ptr             = loc(edge_buf(qsize*nlev+1, in + 3))
        edge_buf_in_4_ptr             = loc(edge_buf(qsize*nlev+1, in + 4))
        swest_buf_ptr                 = loc(edge_buf(qsize*nlev+1, getmapP(5) + 1))
        seast_buf_ptr                 = loc(edge_buf(qsize*nlev+1, getmapP(6) + 1))
        neast_buf_ptr                 = loc(edge_buf(qsize*nlev+1, getmapP(8) + 1))
        nwest_buf_ptr                 = loc(edge_buf(qsize*nlev+1, getmapP(7) + 1))

        !$ACC DATA copy(elem_derived_eta_dot_dpdn) copyin(elem_rspheremp,getmapP,is,ie,iw,in,edge_buf_is_1,edge_buf_is_2,edge_buf_is_3,edge_buf_is_4,edge_buf_iw_1,edge_buf_iw_2,edge_buf_iw_3,edge_buf_iw_4,edge_buf_ie_1,edge_buf_ie_2,edge_buf_ie_3,edge_buf_ie_4,edge_buf_in_1,edge_buf_in_2,edge_buf_in_3,edge_buf_in_4,swest_buf,seast_buf,neast_buf,nwest_buf)
        call my_edgeVunpack_euler(edge_nlyr, edge_nbuf, getmapP(:), &
          elem_derived_eta_dot_dpdn(:,:,:), nlev, qsize*nlev, my_swest, my_max_corner_elem, &
          swest_buf(:), seast_buf(:), neast_buf(:), nwest_buf(:), &
          edge_buf_is_1(:), edge_buf_is_2(:), edge_buf_is_3(:), edge_buf_is_4(:), &
          edge_buf_iw_1(:), edge_buf_iw_2(:), edge_buf_iw_3(:), edge_buf_iw_4(:), &
          edge_buf_ie_1(:), edge_buf_ie_2(:), edge_buf_ie_3(:), edge_buf_ie_4(:), &
          edge_buf_in_1(:), edge_buf_in_2(:), edge_buf_in_3(:), edge_buf_in_4(:))

        do k = 1, nlev
          elem_derived_eta_dot_dpdn(:,:,k) = elem_derived_eta_dot_dpdn(:,:,k) * elem_rspheremp(:,:)
        end do
        !$ACC END DATA
      end do !end ie
      !$ACC END PARALLEL LOOP

    else if (DSSopt == DSSomega) then
      !$ACC PARALLEL LOOP
      do iee = nets, nete
        !iee = calList(i) 
        elem_derived_omega_p_ptr = loc(my_elem(iee)%derived%omega_p)
        elem_rspheremp_ptr       = loc(my_elem(iee)%rspheremp)
        getmapP_ptr              = loc(my_elem(iee)%desc%getmapP)
        is_ptr                   = loc(my_elem(iee)%desc%getmapP(my_south))
        ie_ptr                   = loc(my_elem(iee)%desc%getmapP(my_east))
        iw_ptr                   = loc(my_elem(iee)%desc%getmapP(my_west))
        in_ptr                   = loc(my_elem(iee)%desc%getmapP(my_north))
        edge_buf_is_1_ptr        = loc(edge_buf(qsize*nlev+1, is + 1))
        edge_buf_is_2_ptr        = loc(edge_buf(qsize*nlev+1, is + 2))
        edge_buf_is_3_ptr        = loc(edge_buf(qsize*nlev+1, is + 3))
        edge_buf_is_4_ptr        = loc(edge_buf(qsize*nlev+1, is + 4))
        edge_buf_iw_1_ptr        = loc(edge_buf(qsize*nlev+1, iw + 1))
        edge_buf_iw_2_ptr        = loc(edge_buf(qsize*nlev+1, iw + 2))
        edge_buf_iw_3_ptr        = loc(edge_buf(qsize*nlev+1, iw + 3))
        edge_buf_iw_4_ptr        = loc(edge_buf(qsize*nlev+1, iw + 4))
        edge_buf_ie_1_ptr        = loc(edge_buf(qsize*nlev+1, ie + 1))
        edge_buf_ie_2_ptr        = loc(edge_buf(qsize*nlev+1, ie + 2))
        edge_buf_ie_3_ptr        = loc(edge_buf(qsize*nlev+1, ie + 3))
        edge_buf_ie_4_ptr        = loc(edge_buf(qsize*nlev+1, ie + 4))
        edge_buf_in_1_ptr        = loc(edge_buf(qsize*nlev+1, in + 1))
        edge_buf_in_2_ptr        = loc(edge_buf(qsize*nlev+1, in + 2))
        edge_buf_in_3_ptr        = loc(edge_buf(qsize*nlev+1, in + 3))
        edge_buf_in_4_ptr        = loc(edge_buf(qsize*nlev+1, in + 4))
        swest_buf_ptr            = loc(edge_buf(qsize*nlev+1, getmapP(5) + 1))
        seast_buf_ptr            = loc(edge_buf(qsize*nlev+1, getmapP(6) + 1))
        neast_buf_ptr            = loc(edge_buf(qsize*nlev+1, getmapP(8) + 1))
        nwest_buf_ptr            = loc(edge_buf(qsize*nlev+1, getmapP(7) + 1))

        !$ACC DATA copy(elem_derived_omega_p) copyin(elem_rspheremp,getmapP,is,ie,iw,in,edge_buf_is_1,edge_buf_is_2,edge_buf_is_3,edge_buf_is_4,edge_buf_iw_1,edge_buf_iw_2,edge_buf_iw_3,edge_buf_iw_4,edge_buf_ie_1,edge_buf_ie_2,edge_buf_ie_3,edge_buf_ie_4,edge_buf_in_1,edge_buf_in_2,edge_buf_in_3,edge_buf_in_4,swest_buf,seast_buf,neast_buf,nwest_buf)
        call my_edgeVunpack_euler(edge_nlyr, edge_nbuf, getmapP(:), &
          elem_derived_omega_p(:,:,:), nlev, qsize*nlev, my_swest, my_max_corner_elem, &
          swest_buf(:), seast_buf(:), neast_buf(:), nwest_buf(:), &
          edge_buf_is_1(:), edge_buf_is_2(:), edge_buf_is_3(:), edge_buf_is_4(:), &
          edge_buf_iw_1(:), edge_buf_iw_2(:), edge_buf_iw_3(:), edge_buf_iw_4(:), &
          edge_buf_ie_1(:), edge_buf_ie_2(:), edge_buf_ie_3(:), edge_buf_ie_4(:), &
          edge_buf_in_1(:), edge_buf_in_2(:), edge_buf_in_3(:), edge_buf_in_4(:))

        do k = 1, nlev
          elem_derived_omega_p(:,:,k) = elem_derived_omega_p(:,:,k) * elem_rspheremp(:,:)
        end do
        !$ACC END DATA
      end do

    else if (DSSopt == DSSdiv_vdp_ave) then
      !call system_clock(count_start, count_rate, count_max)

      !$ACC PARALLEL LOOP
      do iee = nets, nete
        !iee = calList(i) 
        elem_derived_divdp_proj_ptr = loc(my_elem(iee)%derived%divdp_proj)
        elem_rspheremp_ptr          = loc(my_elem(iee)%rspheremp)
        getmapP_ptr                 = loc(my_elem(iee)%desc%getmapP)
        is_ptr                      = loc(my_elem(iee)%desc%getmapP(my_south))
        ie_ptr                      = loc(my_elem(iee)%desc%getmapP(my_east))
        iw_ptr                      = loc(my_elem(iee)%desc%getmapP(my_west))
        in_ptr                      = loc(my_elem(iee)%desc%getmapP(my_north))
        edge_buf_is_1_ptr           = loc(edge_buf(qsize*nlev+1, is + 1))
        edge_buf_is_2_ptr           = loc(edge_buf(qsize*nlev+1, is + 2))
        edge_buf_is_3_ptr           = loc(edge_buf(qsize*nlev+1, is + 3))
        edge_buf_is_4_ptr           = loc(edge_buf(qsize*nlev+1, is + 4))
        edge_buf_iw_1_ptr           = loc(edge_buf(qsize*nlev+1, iw + 1))
        edge_buf_iw_2_ptr           = loc(edge_buf(qsize*nlev+1, iw + 2))
        edge_buf_iw_3_ptr           = loc(edge_buf(qsize*nlev+1, iw + 3))
        edge_buf_iw_4_ptr           = loc(edge_buf(qsize*nlev+1, iw + 4))
        edge_buf_ie_1_ptr           = loc(edge_buf(qsize*nlev+1, ie + 1))
        edge_buf_ie_2_ptr           = loc(edge_buf(qsize*nlev+1, ie + 2))
        edge_buf_ie_3_ptr           = loc(edge_buf(qsize*nlev+1, ie + 3))
        edge_buf_ie_4_ptr           = loc(edge_buf(qsize*nlev+1, ie + 4))
        edge_buf_in_1_ptr           = loc(edge_buf(qsize*nlev+1, in + 1))
        edge_buf_in_2_ptr           = loc(edge_buf(qsize*nlev+1, in + 2))
        edge_buf_in_3_ptr           = loc(edge_buf(qsize*nlev+1, in + 3))
        edge_buf_in_4_ptr           = loc(edge_buf(qsize*nlev+1, in + 4))
        swest_buf_ptr               = loc(edge_buf(qsize*nlev+1, getmapP(5) + 1))
        seast_buf_ptr               = loc(edge_buf(qsize*nlev+1, getmapP(6) + 1))
        neast_buf_ptr               = loc(edge_buf(qsize*nlev+1, getmapP(8) + 1))
        nwest_buf_ptr               = loc(edge_buf(qsize*nlev+1, getmapP(7) + 1))

        !$ACC DATA copy(elem_derived_divdp_proj) copyin(elem_rspheremp,getmapP,is,ie,iw,in,edge_buf_is_1,edge_buf_is_2,edge_buf_is_3,edge_buf_is_4,edge_buf_iw_1,edge_buf_iw_2,edge_buf_iw_3,edge_buf_iw_4,edge_buf_ie_1,edge_buf_ie_2,edge_buf_ie_3,edge_buf_ie_4,edge_buf_in_1,edge_buf_in_2,edge_buf_in_3,edge_buf_in_4,swest_buf,seast_buf,neast_buf,nwest_buf)
        call my_edgeVunpack_euler(edge_nlyr, edge_nbuf, getmapP(:), &
          elem_derived_divdp_proj(:,:,:), nlev, qsize*nlev, my_swest, my_max_corner_elem, &
          swest_buf(:), seast_buf(:), neast_buf(:), nwest_buf(:), &
          edge_buf_is_1(:), edge_buf_is_2(:), edge_buf_is_3(:), edge_buf_is_4(:), &
          edge_buf_iw_1(:), edge_buf_iw_2(:), edge_buf_iw_3(:), edge_buf_iw_4(:), &
          edge_buf_ie_1(:), edge_buf_ie_2(:), edge_buf_ie_3(:), edge_buf_ie_4(:), &
          edge_buf_in_1(:), edge_buf_in_2(:), edge_buf_in_3(:), edge_buf_in_4(:))

        do k = 1, nlev
          elem_derived_divdp_proj(:,:,k) = elem_derived_divdp_proj(:,:,k) * elem_rspheremp(:,:)
        end do
        !$ACC END DATA

      end do
      !$ACC END PARALLEL LOOP

      !call system_clock(count_stop, count_rate, count_max)
      !write(*,*), 'acc count = ', count_stop - count_start
    end if

  end subroutine


  subroutine cal_qmin_qmax_rhs_multiplier_0(nets, nete, my_qsize, my_nlev, &
                                            my_np, Qtens_biharmonic, my_qmin, my_qmax)
    integer, intent(in) :: nets
    integer, intent(in) :: nete
    integer, intent(in) :: my_qsize
    integer, intent(in) :: my_nlev
    integer, intent(in) :: my_np
    real(kind=8), dimension(my_np, my_np, my_nlev, my_qsize,nets:nete), intent(in) :: Qtens_biharmonic
    real(kind=8), dimension(my_nlev, my_qsize, nets:nete), intent(out) :: my_qmin
    real(kind=8), dimension(my_nlev, my_qsize, nets:nete), intent(out) :: my_qmax

    !local
    integer :: ie, q, k

    !!$ACC PARALLEL LOOP collapse(2) copyout(my_qmin, my_qmax) copyin(Qtens_biharmonic)
    do ie = nets, nete
      do q = 1, my_qsize
        do k = 1, my_nlev
          my_qmin(k,q,ie)=minval(Qtens_biharmonic(:,:,k,q,ie))
          my_qmax(k,q,ie)=maxval(Qtens_biharmonic(:,:,k,q,ie))
          my_qmin(k,q,ie)=max(my_qmin(k,q,ie),0d0)
        enddo
      enddo
    enddo
    !!$ACC END PARALLEL LOOP
  end subroutine

  subroutine cal_qmin_qmax_rhs_multiplier_1(nets, nete, my_qsize, my_nlev, &
                                            my_np, Qtens_biharmonic, my_qmin, my_qmax)
    integer, intent(in) :: nets
    integer, intent(in) :: nete
    integer, intent(in) :: my_qsize
    integer, intent(in) :: my_nlev
    integer, intent(in) :: my_np
    real(kind=8), dimension(my_np, my_np, my_nlev, my_qsize,nets:nete), intent(in) :: Qtens_biharmonic
    real(kind=8), dimension(my_nlev, my_qsize, nets:nete), intent(inout) :: my_qmin
    real(kind=8), dimension(my_nlev, my_qsize, nets:nete), intent(inout) :: my_qmax

    !local
    integer :: ie, q, k

    !!$ACC PARALLEL LOOP collapse(2) copy(my_qmin, my_qmax) copyin(Qtens_biharmonic)
    do ie = nets, nete
      do q = 1, my_qsize
        do k = 1 , my_nlev    !  Loop index added with implicit inversion (AAM)
          my_qmin(k,q,ie)=min(my_qmin(k,q,ie),minval(Qtens_biharmonic(:,:,k,q,ie)))
          my_qmin(k,q,ie)=max(my_qmin(k,q,ie),0d0)
          my_qmax(k,q,ie)=max(my_qmax(k,q,ie),maxval(Qtens_biharmonic(:,:,k,q,ie)))
        enddo
      enddo
    enddo
    !!$ACC END PARALLEL LOOP
  end subroutine

  subroutine cal_qtens_biharmonic_nu_p(nets, nete, my_qsize, my_nlev, &
                                            my_np, Qtens_biharmonic, &
                                          elem_derived_dpdiss_ave, &
                                        hvcoord_hyai, hvcoord_hybi, hvcoord_ps0)
    integer, intent(in) :: nets
    integer, intent(in) :: nete
    integer, intent(in) :: my_qsize
    integer, intent(in) :: my_nlev
    integer, intent(in) :: my_np
    real(kind=8), dimension(my_np, my_np, my_nlev, nets:nete), intent(in)   :: elem_derived_dpdiss_ave
    real(kind=8), dimension(my_nlev+1), intent(in)         :: hvcoord_hyai
    real(kind=8), dimension(my_nlev+1), intent(in)         :: hvcoord_hybi
    real(kind=8), intent(in)                               :: hvcoord_ps0
    real(kind=8), dimension(my_np, my_np, my_nlev, my_qsize, nets:nete), intent(inout) :: Qtens_biharmonic

    !local
    integer :: ie, q, k

    !!$ACC PARALLEL LOOP copy(Qtens_biharmonic) copyin(elem_derived_dpdiss_ave, hvcoord_hyai, hvcoord_hybi) annotate(entire(hvcoord_hyai, hvcoord_hybi)) collapse(2) tile(q:1)
    do ie = nets , nete
      do q = 1 , my_qsize
        do k = 1 , my_nlev
          Qtens_biharmonic(:,:,k,q,ie) = Qtens_biharmonic(:,:,k,q,ie) * elem_derived_dpdiss_ave(:,:,k, ie) / &
              ( (hvcoord_hyai(k+1) - hvcoord_hyai(k)) * hvcoord_ps0 + &
                ( hvcoord_hybi(k+1) - hvcoord_hybi(k) )*hvcoord_ps0 )
        enddo
      enddo
    enddo
    !!$ACC END PARALLEL LOOP
  end subroutine

  subroutine cal_qtens_biharmonic(nets, nete, my_qsize, my_nlev, &
                                            my_np, Qtens_biharmonic, &
                                            rhs_viss, dt, my_nu_q, hvcoord_hyai, hvcoord_hybi, &
                                            hvcoord_ps0, elem_spheremp_array)
    integer, intent(in) :: nets
    integer, intent(in) :: nete
    integer, intent(in) :: my_qsize
    integer, intent(in) :: my_nlev
    integer, intent(in) :: my_np
    real(kind=8), dimension(my_np, my_np, my_nlev, my_qsize, nets:nete), intent(inout) :: Qtens_biharmonic
    integer, intent(in) :: rhs_viss
    real(kind=8), intent(in) :: dt
    real(kind=8), intent(in) :: my_nu_q
    real(kind=8), dimension(my_nlev+1), intent(in)         :: hvcoord_hyai
    real(kind=8), dimension(my_nlev+1), intent(in)         :: hvcoord_hybi
    real(kind=8), intent(in)                               :: hvcoord_ps0
    real(kind=8), dimension(my_np, my_np, nets:nete), intent(in) :: elem_spheremp_array

    !local
    integer :: ie, q, k

   !!$ACC PARALLEL LOOP copy(Qtens_biharmonic) copyin(hvcoord_hyai, hvcoord_hybi, elem_spheremp_array) annotate(entire(hvcoord_hyai, hvcoord_hybi)) collapse(2) tile(q:1)
   do ie = nets , nete
    do k = 1 , my_nlev    !  Loop inversion (AAM)
      do q = 1 , my_qsize
        qtens_biharmonic(:,:,k,q,ie) = &
                 -rhs_viss*dt*my_nu_q*( ( hvcoord_hyai(k+1) - hvcoord_hyai(k))*hvcoord_ps0 + ( hvcoord_hybi(k+1) - hvcoord_hybi(k))*hvcoord_ps0) *Qtens_biharmonic(:,:,k,q,ie) / elem_spheremp_array(:,:,ie)
      enddo
    enddo
  enddo
  !!$ACC END PARALLEL LOOP
  end subroutine
  subroutine my_divergence_sphere(gradQ, derivDvvi, elemmetdevt, elemDinv, elemrmetdet, my_rrearth, div)
    ! ===========================
    ! edited by Conghui
    ! used as a sub function in acc
    !==============================
    implicit none
    real(8), intent(in) :: gradQ(4,4,2)  ! in lat-lon coordinates
    real(8), intent(in) :: derivDvvi(4, 4)
    real(8), intent(in) :: elemmetdevt(4, 4)
    real(8), intent(in) :: elemDinv(2,2,4,4)     ! Map vector field on the sphere to covariant gradQ on cube
    real(8), intent(in) :: elemrmetdet(4,4)      ! 1/metdet on velocity pressure grid
    real(8), intent(out) :: div(4,4)
    real(8), intent(in) :: my_rrearth

    ! Local
    integer :: i
    integer :: j
    integer :: l

    real(8) ::  dudx00
    real(8) ::  dvdy00
    real(8) ::  gv(4,4,2),vvtemp(4,4)

    ! convert to contra variant form and multiply by g
    do j=1,4
     do i=1,4
        gv(i,j,1)=elemmetdevt(i,j)*(elemDinv(1,1,i,j)*gradQ(i,j,1) + elemDinv(1,2,i,j)*gradQ(i,j,2))
        gv(i,j,2)=elemmetdevt(i,j)*(elemDinv(2,1,i,j)*gradQ(i,j,1) + elemDinv(2,2,i,j)*gradQ(i,j,2))
     enddo
    enddo
    ! compute d/dx and d/dy
    do j=1,4
     do l=1,4
        dudx00=0.0d0
        dvdy00=0.0d0
        do i=1,4
           dudx00 = dudx00 + derivDvvi(i,l  )*gv(i,j  ,1)
           dvdy00 = dvdy00 + derivDvvi(i,l  )*gv(j  ,i,2)
        end do
        div(l  ,j  ) = dudx00
        vvtemp(j  ,l  ) = dvdy00
     end do
    end do

    do j=1,4
     do i=1,4
        div(i,j)=(div(i,j)+vvtemp(i,j))*(elemrmetdet(i,j)*my_rrearth)
     end do
    end do
  end subroutine

  subroutine my_euler_step_acc(nets, nete, rhs_multiplier, rhs_viss, &
               dt, my_rrearth, my_nu_p, my_nu_q, &
               deriv_dvv, Qtens_biharmonic, &
               my_qmin,my_qmax,elem_array,  &
               my_n0_qdp, my_np1_qdp, const_qsize, calList, calLen)

    implicit none
    !integer, parameter :: const_qsize = 25
    integer, parameter :: const_np    = 4
    integer, parameter :: const_nlev  = constLev
    integer, parameter :: tile_size   = 2

    integer, intent(in)       :: nets, nete, rhs_multiplier, rhs_viss, my_n0_qdp, my_np1_qdp, const_qsize, calLen
    integer, dimension(calLen), intent(in) :: calList
    real(kind=8), intent(in)  :: dt, my_rrearth, my_nu_p, my_nu_q
    real(kind=8), dimension(const_np,const_np), intent(in)                                   :: deriv_dvv
    real(kind=8), dimension(const_np,const_np,const_nlev,const_qsize,nets:nete), intent(in)  :: Qtens_biharmonic
    real(kind=8), dimension(const_nlev,const_qsize,nets:nete), intent(inout)                 :: my_qmin
    real(kind=8), dimension(const_nlev,const_qsize,nets:nete), intent(inout)                 :: my_qmax
    integer(kind=8), dimension(11,nets:nete), intent(inout) :: elem_array

    ! cray pointers
    real(kind=8), dimension(const_np,const_np,const_nlev) :: elem_derived_dp
    pointer(elem_derived_dp_ptr, elem_derived_dp)

    real(kind=8), dimension(const_np,const_np,const_nlev) :: elem_derived_divdp_proj
    pointer(elem_derived_divdp_proj_ptr, elem_derived_divdp_proj)

    real(kind=8), dimension(const_np,const_np,const_nlev)        :: elem_derived_divdp
    pointer(elem_derived_divdp_ptr, elem_derived_divdp)

    real(kind=8), dimension(const_np,const_np,const_nlev)        :: elem_derived_dpdiss_biharmonic
    pointer(elem_derived_dpdiss_biharmonic_ptr, elem_derived_dpdiss_biharmonic)

    real(kind=8), dimension(const_np,const_np,2,const_nlev)      :: elem_derived_vn0
    pointer(elem_derived_vn0_ptr, elem_derived_vn0)

    real(kind=8), dimension(2,2,const_np,const_np)            :: elem_Dinv
    pointer(elem_Dinv_ptr, elem_Dinv)

    real(kind=8), dimension(const_np,const_np)                :: elem_metdet
    pointer(elem_metdet_ptr, elem_metdet)

    real(kind=8), dimension(const_np,const_np)                :: elem_rmetdet
    pointer(elem_rmetdet_ptr, elem_rmetdet)

    real(kind=8), dimension(const_np,const_np)                :: elem_spheremp
    pointer(elem_spheremp_ptr, elem_spheremp)

    real(kind=8), dimension(const_np,const_np,const_nlev,const_qsize) :: elem_state_Qdp_in
    pointer(elem_state_Qdp_in_ptr, elem_state_Qdp_in)

    real(kind=8), dimension(const_np,const_np,const_nlev,const_qsize) :: elem_state_Qdp_out
    pointer(elem_state_Qdp_out_ptr, elem_state_Qdp_out)

    ! local
    real(kind=8), dimension(const_np,const_np)   :: qdp_val
    real(kind=8), dimension(const_np,const_np,2) :: gradQ
    real(kind=8), dimension(const_np,const_np)   :: Qtens
    real(kind=8), dimension(const_np,const_np)   :: dp_star
    integer :: ie, k, q, q1,  qs, qe, iee
    integer :: startLev, endLev
    integer :: time1, time2
    !----for divergence----
    
    integer :: i
    integer :: j
    integer :: l
    
    real(8) ::  dudx00
    real(8) ::  dvdy00
    real(8) ::  gv(4,4,2),vvtemp(4,4)
    !call penv_cg0_gld_count_init()
    !call penv_cg0_gld_count(time1)
    !$ACC PARALLEL LOOP tile(q:2,ie:1)  collapse(2) local(qdp_val, gradQ, Qtens, dp_star, gv, vvtemp) copyin(deriv_dvv) annotate(entire(deriv_dvv); readonly(deriv_dvv))
    do ie = nets , nete
      do q = 1 , const_qsize

        
        !$ACC DATA COPYIN(elem_array(*,ie))        
        elem_derived_dp_ptr                = elem_array(1,ie)
        elem_derived_divdp_proj_ptr        = elem_array(2,ie)
        elem_derived_divdp_ptr             = elem_array(3,ie)
        elem_derived_dpdiss_biharmonic_ptr = elem_array(4,ie)
        elem_derived_vn0_ptr               = elem_array(5,ie)
        elem_Dinv_ptr                      = elem_array(6,ie)
        elem_metdet_ptr                    = elem_array(7,ie)
        elem_rmetdet_ptr                   = elem_array(8,ie)
        elem_spheremp_ptr                  = elem_array(9,ie)
        elem_state_Qdp_in_ptr              = elem_array(10,ie)
        elem_state_Qdp_out_ptr             = elem_array(11,ie)


        !$ACC DATA COPYIN(elem_Dinv, elem_metdet, elem_rmetdet, elem_spheremp, elem_derived_dp(*,*,1:partLev), elem_derived_divdp_proj(*,*,1:partLev), elem_derived_divdp(*,*,1:partLev), elem_derived_dpdiss_biharmonic(*,*,1:partLev), elem_derived_vn0(*,*,*,1:partLev))
        !$ACC DATA COPYOUT(elem_state_Qdp_out(*,*,1:partLev,q)) COPY(my_qmin(1:partLev,q,ie), my_qmax(1:partLev,q,ie)) COPYIN(Qtens_biharmonic(*,*,1:partLev,q,ie), elem_state_Qdp_in(*,*,1:partLev,q))
        do k = 1 , partLev
          qdp_val = elem_state_Qdp_in(:,:,k,q)
          gradQ(:,:,1) = elem_derived_vn0(:,:,1,k) / (elem_derived_dp(:,:,k) - rhs_multiplier * dt * elem_derived_divdp_proj(:,:,k)) * qdp_val
          gradQ(:,:,2) = elem_derived_vn0(:,:,2,k) / (elem_derived_dp(:,:,k) - rhs_multiplier * dt * elem_derived_divdp_proj(:,:,k)) * qdp_val

          call my_divergence_sphere(gradQ, deriv_dvv, elem_metdet(:,:), elem_Dinv(:,:,:,:), elem_rmetdet(:,:), my_rrearth, dp_star(:,:))

          !----------------my_divergence_sphere-------------------

         ! ! convert to contra variant form and multiply by g
         ! do j=1,4
         !  do i=1,4
         !     gv(i,j,1)=elem_metdet(i,j)*(elem_Dinv(1,1,i,j)*gradQ(i,j,1) + elem_Dinv(1,2,i,j)*gradQ(i,j,2))
         !     gv(i,j,2)=elem_metdet(i,j)*(elem_Dinv(2,1,i,j)*gradQ(i,j,1) + elem_Dinv(2,2,i,j)*gradQ(i,j,2))
         !  enddo
         ! enddo
         ! ! compute d/dx and d/dy
         ! do j=1,4
         !  do l=1,4
         !     dudx00=0.0d0
         !     dvdy00=0.0d0
         !     do i=1,4
         !        dudx00 = dudx00 + deriv_Dvv(i,l  )*gv(i,j  ,1)
         !        dvdy00 = dvdy00 + deriv_Dvv(i,l  )*gv(j  ,i,2)
         !     end do
         !     dp_star(l  ,j  ) = dudx00
         !     vvtemp(j  ,l  ) = dvdy00
         !  end do
         ! end do
      
         ! do j=1,4
         !  do i=1,4
         !     dp_star(i,j)=(dp_star(i,j)+vvtemp(i,j))*(elem_rmetdet(i,j)*my_rrearth)
         !  end do
         ! end do

          !----------------end my_divergence_sphere----------------


          Qtens(:,:) = qdp_val - dt * dp_star(:,:)
          if ( rhs_viss /= 0 ) Qtens(:,:) = Qtens(:,:) + Qtens_biharmonic(:,:,k,q,ie)
          dp_star(:,:) = (elem_derived_dp(:,:,k) - rhs_multiplier * dt * elem_derived_divdp_proj(:,:,k)) - dt * elem_derived_divdp(:,:,k)
          if ( my_nu_p > 0 .and. rhs_viss /= 0 ) then
            dp_star(:,:) = dp_star(:,:) - rhs_viss * dt * my_nu_q * elem_derived_dpdiss_biharmonic(:,:,k) / elem_spheremp(:,:)
          endif

          call limiter_optim_iter_full( Qtens(:,:) , elem_spheremp(:,:) , my_qmin(k,q,ie), my_qmax(k,q,ie) , dp_star(:,:))

          elem_state_Qdp_out(:,:,k,q) = elem_spheremp(:,:) * Qtens(:,:)
        enddo ! k
        !$ACC END DATA
        !$ACC END DATA
        !$ACC END DATA
      enddo ! end qsize
    end do
    !$ACC END PARALLEL LOOP
    !call penv_cg0_gld_count(time2)
    !write(*,*) "Asher gld",time2-time1

    
    !call penv_cg0_gld_count(time1)
    !$ACC PARALLEL LOOP  collapse(2) local(qdp_val, gradQ, Qtens, dp_star) copyin(deriv_dvv, calList) annotate(entire(deriv_dvv); readonly(deriv_dvv))
    do ie =  nets, nete
      do q = 1 , const_qsize

        !ie = calList(iee) 
        
        !$ACC DATA COPYIN(elem_array(*,ie))        
        elem_derived_dp_ptr                = elem_array(1,ie)
        elem_derived_divdp_proj_ptr        = elem_array(2,ie)
        elem_derived_divdp_ptr             = elem_array(3,ie)
        elem_derived_dpdiss_biharmonic_ptr = elem_array(4,ie)
        elem_derived_vn0_ptr               = elem_array(5,ie)
        elem_Dinv_ptr                      = elem_array(6,ie)
        elem_metdet_ptr                    = elem_array(7,ie)
        elem_rmetdet_ptr                   = elem_array(8,ie)
        elem_spheremp_ptr                  = elem_array(9,ie)
        elem_state_Qdp_in_ptr              = elem_array(10,ie)
        elem_state_Qdp_out_ptr             = elem_array(11,ie)

        !$ACC DATA COPYIN(elem_Dinv, elem_metdet, elem_rmetdet, elem_spheremp, elem_derived_dp(*,*,33:64), elem_derived_divdp_proj(*,*,33:64), elem_derived_divdp(*,*,33:64), elem_derived_dpdiss_biharmonic(*,*,33:64), elem_derived_vn0(*,*,*,33:64))
        !$ACC DATA COPYOUT(elem_state_Qdp_out(*,*,33:64,q)) COPY(my_qmin(33:64,q,ie), my_qmax(33:64,q,ie)) COPYIN(Qtens_biharmonic(*,*,33:64,q,ie), elem_state_Qdp_in(*,*,33:64,q))
        do k = 33, 64 
          qdp_val = elem_state_Qdp_in(:,:,k,q)
          gradQ(:,:,1) = elem_derived_vn0(:,:,1,k) / (elem_derived_dp(:,:,k) - rhs_multiplier * dt * elem_derived_divdp_proj(:,:,k)) * qdp_val
          gradQ(:,:,2) = elem_derived_vn0(:,:,2,k) / (elem_derived_dp(:,:,k) - rhs_multiplier * dt * elem_derived_divdp_proj(:,:,k)) * qdp_val

          call my_divergence_sphere(gradQ, deriv_dvv, elem_metdet(:,:), elem_Dinv(:,:,:,:), elem_rmetdet(:,:), my_rrearth, dp_star(:,:))

          Qtens(:,:) = qdp_val - dt * dp_star(:,:)
          if ( rhs_viss /= 0 ) Qtens(:,:) = Qtens(:,:) + Qtens_biharmonic(:,:,k,q,ie)
          dp_star(:,:) = (elem_derived_dp(:,:,k) - rhs_multiplier * dt * elem_derived_divdp_proj(:,:,k)) - dt * elem_derived_divdp(:,:,k)
          if ( my_nu_p > 0 .and. rhs_viss /= 0 ) then
            dp_star(:,:) = dp_star(:,:) - rhs_viss * dt * my_nu_q * elem_derived_dpdiss_biharmonic(:,:,k) / elem_spheremp(:,:)
          endif

          call limiter_optim_iter_full( Qtens(:,:) , elem_spheremp(:,:) , my_qmin(k,q,ie), my_qmax(k,q,ie) , dp_star(:,:))

          elem_state_Qdp_out(:,:,k,q) = elem_spheremp(:,:) * Qtens(:,:)
        enddo ! k
        !$ACC END DATA
        !$ACC END DATA
        !$ACC END DATA
      enddo ! end qsize
    end do
    !$ACC END PARALLEL LOOP
    !call penv_cg0_gld_count(time2)
    !write(*,*) "Asher gld2",time2-time1

    !$ACC PARALLEL LOOP  collapse(2) local(qdp_val, gradQ, Qtens, dp_star) copyin(deriv_dvv, calList) annotate(entire(deriv_dvv); readonly(deriv_dvv))
    do ie = nets , nete
      do q = 1 , const_qsize

        !ie = calList(iee) 
        
        !$ACC DATA COPYIN(elem_array(*,ie))        
        elem_derived_dp_ptr                = elem_array(1,ie)
        elem_derived_divdp_proj_ptr        = elem_array(2,ie)
        elem_derived_divdp_ptr             = elem_array(3,ie)
        elem_derived_dpdiss_biharmonic_ptr = elem_array(4,ie)
        elem_derived_vn0_ptr               = elem_array(5,ie)
        elem_Dinv_ptr                      = elem_array(6,ie)
        elem_metdet_ptr                    = elem_array(7,ie)
        elem_rmetdet_ptr                   = elem_array(8,ie)
        elem_spheremp_ptr                  = elem_array(9,ie)
        elem_state_Qdp_in_ptr              = elem_array(10,ie)
        elem_state_Qdp_out_ptr             = elem_array(11,ie)

        !$ACC DATA COPYIN(elem_Dinv, elem_metdet, elem_rmetdet, elem_spheremp, elem_derived_dp(*,*,65:96), elem_derived_divdp_proj(*,*,65:96), elem_derived_divdp(*,*,65:96), elem_derived_dpdiss_biharmonic(*,*,65:96), elem_derived_vn0(*,*,*,65:96))
        !$ACC DATA COPYOUT(elem_state_Qdp_out(*,*,65:96)) COPY(my_qmin(65:96,q,ie), my_qmax(65:96,q,ie)) COPYIN(Qtens_biharmonic(*,*,65:96,q,ie), elem_state_Qdp_in(*,*,65:96,q))
        do k = 65, 96 
          qdp_val = elem_state_Qdp_in(:,:,k,q)
          gradQ(:,:,1) = elem_derived_vn0(:,:,1,k) / (elem_derived_dp(:,:,k) - rhs_multiplier * dt * elem_derived_divdp_proj(:,:,k)) * qdp_val
          gradQ(:,:,2) = elem_derived_vn0(:,:,2,k) / (elem_derived_dp(:,:,k) - rhs_multiplier * dt * elem_derived_divdp_proj(:,:,k)) * qdp_val

          call my_divergence_sphere(gradQ, deriv_dvv, elem_metdet(:,:), elem_Dinv(:,:,:,:), elem_rmetdet(:,:), my_rrearth, dp_star(:,:))

          Qtens(:,:) = qdp_val - dt * dp_star(:,:)
          if ( rhs_viss /= 0 ) Qtens(:,:) = Qtens(:,:) + Qtens_biharmonic(:,:,k,q,ie)
          dp_star(:,:) = (elem_derived_dp(:,:,k) - rhs_multiplier * dt * elem_derived_divdp_proj(:,:,k)) - dt * elem_derived_divdp(:,:,k)
          if ( my_nu_p > 0 .and. rhs_viss /= 0 ) then
            dp_star(:,:) = dp_star(:,:) - rhs_viss * dt * my_nu_q * elem_derived_dpdiss_biharmonic(:,:,k) / elem_spheremp(:,:)
          endif

          call limiter_optim_iter_full( Qtens(:,:) , elem_spheremp(:,:) , my_qmin(k,q,ie), my_qmax(k,q,ie) , dp_star(:,:))

          elem_state_Qdp_out(:,:,k,q) = elem_spheremp(:,:) * Qtens(:,:)
        enddo ! k
        !$ACC END DATA
        !$ACC END DATA
        !$ACC END DATA
      enddo ! end qsize
    end do
    !$ACC END PARALLEL LOOP


    !$ACC PARALLEL LOOP  collapse(2) local(qdp_val, gradQ, Qtens, dp_star) copyin(deriv_dvv, calList) annotate(entire(deriv_dvv); readonly(deriv_dvv))
    do ie = nets , nete
      do q = 1 , const_qsize

        !ie = calList(iee) 
        
        !$ACC DATA COPYIN(elem_array(*,ie))        
        elem_derived_dp_ptr                = elem_array(1,ie)
        elem_derived_divdp_proj_ptr        = elem_array(2,ie)
        elem_derived_divdp_ptr             = elem_array(3,ie)
        elem_derived_dpdiss_biharmonic_ptr = elem_array(4,ie)
        elem_derived_vn0_ptr               = elem_array(5,ie)
        elem_Dinv_ptr                      = elem_array(6,ie)
        elem_metdet_ptr                    = elem_array(7,ie)
        elem_rmetdet_ptr                   = elem_array(8,ie)
        elem_spheremp_ptr                  = elem_array(9,ie)
        elem_state_Qdp_in_ptr              = elem_array(10,ie)
        elem_state_Qdp_out_ptr             = elem_array(11,ie)

        !$ACC DATA COPYIN(elem_Dinv, elem_metdet, elem_rmetdet, elem_spheremp, elem_derived_dp(*,*,97:128), elem_derived_divdp_proj(*,*,97:128), elem_derived_divdp(*,*,97:128), elem_derived_dpdiss_biharmonic(*,*,97:128), elem_derived_vn0(*,*,*,97:128))
        !$ACC DATA COPYOUT(elem_state_Qdp_out(*,*,97:128,q)) COPY(my_qmin(97:128,q,ie), my_qmax(97:128,q,ie)) COPYIN(Qtens_biharmonic(*,*,97:128,q,ie), elem_state_Qdp_in(*,*,97:128,q))
        do k = 97, 128 
          qdp_val = elem_state_Qdp_in(:,:,k,q)
          gradQ(:,:,1) = elem_derived_vn0(:,:,1,k) / (elem_derived_dp(:,:,k) - rhs_multiplier * dt * elem_derived_divdp_proj(:,:,k)) * qdp_val
          gradQ(:,:,2) = elem_derived_vn0(:,:,2,k) / (elem_derived_dp(:,:,k) - rhs_multiplier * dt * elem_derived_divdp_proj(:,:,k)) * qdp_val

          call my_divergence_sphere(gradQ, deriv_dvv, elem_metdet(:,:), elem_Dinv(:,:,:,:), elem_rmetdet(:,:), my_rrearth, dp_star(:,:))

          Qtens(:,:) = qdp_val - dt * dp_star(:,:)
          if ( rhs_viss /= 0 ) Qtens(:,:) = Qtens(:,:) + Qtens_biharmonic(:,:,k,q,ie)
          dp_star(:,:) = (elem_derived_dp(:,:,k) - rhs_multiplier * dt * elem_derived_divdp_proj(:,:,k)) - dt * elem_derived_divdp(:,:,k)
          if ( my_nu_p > 0 .and. rhs_viss /= 0 ) then
            dp_star(:,:) = dp_star(:,:) - rhs_viss * dt * my_nu_q * elem_derived_dpdiss_biharmonic(:,:,k) / elem_spheremp(:,:)
          endif

          call limiter_optim_iter_full( Qtens(:,:) , elem_spheremp(:,:) , my_qmin(k,q,ie), my_qmax(k,q,ie) , dp_star(:,:))

          elem_state_Qdp_out(:,:,k,q) = elem_spheremp(:,:) * Qtens(:,:)
        enddo ! k
        !$ACC END DATA
        !$ACC END DATA
        !$ACC END DATA
      enddo ! end qsize
    end do
    !$ACC END PARALLEL LOOP

  end subroutine

  subroutine my_edgeVpack_acc(v,vlyr,kptr,desc_putmapP, &
      desc_reverse_north,desc_reverse_south,desc_reverse_east,desc_reverse_west, &
      edge_buf_5,edge_buf_6,edge_buf_7,edge_buf_8, &
      edge_buf_in_1,edge_buf_in_2,edge_buf_in_3,edge_buf_in_4, &
      edge_buf_is_1,edge_buf_is_2,edge_buf_is_3,edge_buf_is_4, &
      edge_buf_ie_1,edge_buf_ie_2,edge_buf_ie_3,edge_buf_ie_4, &
      edge_buf_iw_1,edge_buf_iw_2,edge_buf_iw_3,edge_buf_iw_4)


    integer,  parameter :: west  = 1
    integer,  parameter :: east  = 2
    integer,  parameter :: south = 3
    integer,  parameter :: north = 4

    integer,  parameter :: swest = 5
    integer,  parameter :: seast = 6
    integer,  parameter :: nwest = 7
    integer,  parameter :: neast = 8
    integer,  parameter :: max_corner_elem = 1
    integer,  parameter :: my_np = 4
    integer,                                         intent(in)   :: vlyr
    real (kind=8),dimension(4,4,128),                 intent(in)   :: v
    integer,                                         intent(in)   :: kptr
    integer, dimension(8),                           intent(in)   :: desc_putmapP
    logical,                                         intent(in)   :: desc_reverse_north, desc_reverse_south, desc_reverse_east, desc_reverse_west
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

    is = desc_putmapP(south)
    ie = desc_putmapP(east)
    in = desc_putmapP(north)
    iw = desc_putmapP(west)
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
          edge_buf_in_1(k) = v(1,my_np,k)
          edge_buf_in_2(k) = v(2,my_np,k)
          edge_buf_in_3(k) = v(3,my_np,k)
          edge_buf_in_4(k) = v(4,my_np,k)
          edge_buf_is_1(k) = v(1,1,k)
          edge_buf_is_2(k) = v(2,1,k)
          edge_buf_is_3(k) = v(3,1,k)
          edge_buf_is_4(k) = v(4,1,k)
          edge_buf_ie_1(k) = v(my_np,1,k)
          edge_buf_ie_2(k) = v(my_np,2,k)
          edge_buf_ie_3(k) = v(my_np,3,k)
          edge_buf_ie_4(k) = v(my_np,4,k)
          edge_buf_iw_1(k) = v(1,1,k)
          edge_buf_iw_2(k) = v(1,2,k)
          edge_buf_iw_3(k) = v(1,3,k)
          edge_buf_iw_4(k) = v(1,4,k)
       end do



    if(desc_reverse_south) then
       do k=1,vlyr
             edge_buf_is_4(k)=v(1,1,k)
             edge_buf_is_3(k)=v(2,1,k)
             edge_buf_is_2(k)=v(3,1,k)
             edge_buf_is_1(k)=v(4,1,k)
        enddo
    endif

    if(desc_reverse_east) then
       do k=1,vlyr
             edge_buf_ie_4(k)=v(my_np,1,k)
             edge_buf_ie_3(k)=v(my_np,2,k)
             edge_buf_ie_2(k)=v(my_np,3,k)
             edge_buf_ie_1(k)=v(my_np,4,k)
       enddo
    endif

    if(desc_reverse_north) then
       do k=1,vlyr
             edge_buf_in_4(k)=v(1,my_np,k)
             edge_buf_in_3(k)=v(2,my_np,k)
             edge_buf_in_2(k)=v(3,my_np,k)
             edge_buf_in_1(k)=v(4,my_np,k)
       enddo
    endif

    if(desc_reverse_west) then
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
                edge_buf_6(k)=v(my_np ,1 ,k)
            end do
        end if

! NEAST
    ll = 8
        if (desc_putmapP(ll) /= -1) then
            do k=1,vlyr
                edge_buf_8(k)=v(my_np ,my_np,k)
            end do
        end if

! NWEST
     ll = 7
        if (desc_putmapP(ll) /= -1) then
            do k=1,vlyr
                edge_buf_7(k)=v(1  ,my_np,k)
            end do
        end if


  end subroutine my_edgeVpack_acc
  subroutine my_packCode_acc(np1_qdp,nets,nete,DSSopt,pack_elem_array,pack_buf_array,my_qsize,calLen, calList, calFlag)
  implicit none

  integer,  parameter :: west  = 1
  integer,  parameter :: east  = 2
  integer,  parameter :: south = 3
  integer,  parameter :: north = 4

  integer,  parameter :: swest = 5
  integer,  parameter :: seast = 6
  integer,  parameter :: nwest = 7
  integer,  parameter :: neast = 8
  integer,  parameter :: max_corner_elem = 1

  integer              , intent(in   )         :: np1_qdp, my_qsize, calLen
  integer              , intent(in   )         :: nets
  integer              , intent(in   )         :: nete
  integer              , intent(in   )         :: DSSopt
  integer(kind=8)      , intent(inout), dimension(18,nets:nete)      :: pack_elem_array
  integer(kind=8)      , intent(inout), dimension(20,qsize+1,nets:nete) :: pack_buf_array
  integer, dimension(calLen), intent(in) :: calList
  logical, intent(in):: calFlag

  !------------local------------
  integer :: ie,k,q, iee
  real(kind=8), pointer, dimension(:,:,:)               :: DSSvar
  integer, parameter :: my_nlev=constLev
  !integer, parameter :: my_qsize=25
  integer, parameter :: my_np=4

  !-------cray pointee------
  real(kind=8), dimension(my_np,my_np,my_nlev) :: elem_derived_eta_dot_dpdn
  pointer(elem_derived_eta_dot_dpdn_ptr, elem_derived_eta_dot_dpdn)

  real(kind=8), dimension(my_np,my_np,my_nlev) :: elem_derived_omega_p
  pointer(elem_derived_omega_p_ptr, elem_derived_omega_p)

  real(kind=8), dimension(my_np,my_np,my_nlev) :: elem_derived_divdp_proj
  pointer(elem_derived_divdp_proj_ptr, elem_derived_divdp_proj)

  real(kind=8), dimension(my_np,my_np,my_nlev,my_qsize) :: elem_state_Qdp
  pointer(elem_state_Qdp_ptr, elem_state_Qdp)

  real(kind=8), dimension(my_np,my_np) :: elem_spheremp
  pointer(elem_spheremp_ptr, elem_spheremp)

  logical :: elem_desc_reverse_south
  pointer(elem_desc_reverse_south_ptr, elem_desc_reverse_south)

  logical :: elem_desc_reverse_north
  pointer(elem_desc_reverse_north_ptr, elem_desc_reverse_north)

  logical :: elem_desc_reverse_east
  pointer(elem_desc_reverse_east_ptr, elem_desc_reverse_east)

  logical :: elem_desc_reverse_west
  pointer(elem_desc_reverse_west_ptr, elem_desc_reverse_west)

  integer, dimension(swest+4*max_corner_elem-1) :: elem_desc_putmapP
  pointer(elem_desc_putmapP_ptr, elem_desc_putmapP)

  real(kind=8), dimension(my_nlev) :: edge_buf_5, edge_buf_6, edge_buf_7, edge_buf_8
  pointer(edge_buf_5_ptr, edge_buf_5)
  pointer(edge_buf_6_ptr, edge_buf_6)
  pointer(edge_buf_7_ptr, edge_buf_7)
  pointer(edge_buf_8_ptr, edge_buf_8)

  real(kind=8), dimension(my_nlev) :: edge_buf_in_1, edge_buf_in_2, edge_buf_in_3, edge_buf_in_4
  pointer(edge_buf_in_1_ptr, edge_buf_in_1)
  pointer(edge_buf_in_2_ptr, edge_buf_in_2)
  pointer(edge_buf_in_3_ptr, edge_buf_in_3)
  pointer(edge_buf_in_4_ptr, edge_buf_in_4)


  real(kind=8), dimension(my_nlev) :: edge_buf_is_1, edge_buf_is_2, edge_buf_is_3, edge_buf_is_4
  pointer(edge_buf_is_1_ptr, edge_buf_is_1)
  pointer(edge_buf_is_2_ptr, edge_buf_is_2)
  pointer(edge_buf_is_3_ptr, edge_buf_is_3)
  pointer(edge_buf_is_4_ptr, edge_buf_is_4)

  real(kind=8), dimension(my_nlev) :: edge_buf_ie_1, edge_buf_ie_2, edge_buf_ie_3, edge_buf_ie_4
  pointer(edge_buf_ie_1_ptr, edge_buf_ie_1)
  pointer(edge_buf_ie_2_ptr, edge_buf_ie_2)
  pointer(edge_buf_ie_3_ptr, edge_buf_ie_3)
  pointer(edge_buf_ie_4_ptr, edge_buf_ie_4)

  real(kind=8), dimension(my_nlev) :: edge_buf_iw_1, edge_buf_iw_2, edge_buf_iw_3, edge_buf_iw_4
  pointer(edge_buf_iw_1_ptr, edge_buf_iw_1)
  pointer(edge_buf_iw_2_ptr, edge_buf_iw_2)
  pointer(edge_buf_iw_3_ptr, edge_buf_iw_3)
  pointer(edge_buf_iw_4_ptr, edge_buf_iw_4)


      !$ACC PARALLEL LOOP collapse(2) 
      do ie = nets , nete
          do q = 1,my_qsize
            !!ie = calList(iee)
            !$ACC DATA COPYIN(pack_elem_array(*,ie),pack_buf_array(*,q,ie))
            elem_state_Qdp_ptr            = pack_elem_array(4,ie)
            elem_desc_reverse_south_ptr   = pack_elem_array(6,ie)
            elem_desc_reverse_north_ptr   = pack_elem_array(7,ie)
            elem_desc_reverse_east_ptr    = pack_elem_array(8,ie)
            elem_desc_reverse_west_ptr    = pack_elem_array(9,ie)
            elem_desc_putmapP_ptr         = pack_elem_array(10,ie)
            edge_buf_in_1_ptr = pack_buf_array(1,q,ie)
            edge_buf_in_2_ptr = pack_buf_array(2,q,ie)
            edge_buf_in_3_ptr = pack_buf_array(3,q,ie)
            edge_buf_in_4_ptr = pack_buf_array(4,q,ie)
            edge_buf_5_ptr  = pack_buf_array(5,q,ie)
            edge_buf_6_ptr  = pack_buf_array(6,q,ie)
            edge_buf_7_ptr  = pack_buf_array(7,q,ie)
            edge_buf_8_ptr  = pack_buf_array(8,q,ie)
            edge_buf_is_1_ptr = pack_buf_array(9,q,ie)
            edge_buf_is_2_ptr = pack_buf_array(10,q,ie)
            edge_buf_is_3_ptr = pack_buf_array(11,q,ie)
            edge_buf_is_4_ptr = pack_buf_array(12,q,ie)
            edge_buf_ie_1_ptr = pack_buf_array(13,q,ie)
            edge_buf_ie_2_ptr = pack_buf_array(14,q,ie)
            edge_buf_ie_3_ptr = pack_buf_array(15,q,ie)
            edge_buf_ie_4_ptr = pack_buf_array(16,q,ie)
            edge_buf_iw_1_ptr = pack_buf_array(17,q,ie)
            edge_buf_iw_2_ptr = pack_buf_array(18,q,ie)
            edge_buf_iw_3_ptr = pack_buf_array(19,q,ie)
            edge_buf_iw_4_ptr = pack_buf_array(20,q,ie)
            !$ACC DATA COPYIN(elem_state_Qdp(*,*,*,q), elem_desc_reverse_north, elem_desc_reverse_south, elem_desc_reverse_east, elem_desc_reverse_west, elem_desc_reverse_south,elem_desc_reverse_north,elem_desc_reverse_east,elem_desc_reverse_west,elem_desc_putmapP) COPY(edge_buf_5, edge_buf_6, edge_buf_7, edge_buf_8) COPYOUT(edge_buf_in_1, edge_buf_in_2, edge_buf_in_3, edge_buf_in_4,  edge_buf_is_1, edge_buf_is_2, edge_buf_is_3, edge_buf_is_4, edge_buf_ie_1, edge_buf_ie_2, edge_buf_ie_3, edge_buf_ie_4, edge_buf_iw_1, edge_buf_iw_2, edge_buf_iw_3, edge_buf_iw_4)
            call my_edgeVpack_acc(elem_state_Qdp(:,:,:,q) , my_nlev , my_nlev*(q-1) , elem_desc_putmapP(:), &
             elem_desc_reverse_north, elem_desc_reverse_south, elem_desc_reverse_east, elem_desc_reverse_west, &
             edge_buf_5, edge_buf_6, edge_buf_7, edge_buf_8, &
             edge_buf_in_1, edge_buf_in_2, edge_buf_in_3, edge_buf_in_4, &
             edge_buf_is_1, edge_buf_is_2, edge_buf_is_3, edge_buf_is_4, &
             edge_buf_ie_1, edge_buf_ie_2, edge_buf_ie_3, edge_buf_ie_4, &
             edge_buf_iw_1, edge_buf_iw_2, edge_buf_iw_3, edge_buf_iw_4)
           !$ACC END DATA
           !$ACC END DATA
        enddo
     enddo
     !$ACC END PARALLEL LOOP

     if (DSSopt == 1) then
       !$ACC PARALLEL LOOP  
       do ie = nets, nete 
          !ie = calList(iee)
          q=my_qsize+1
          !$ACC DATA COPYIN(pack_elem_array(*,ie), pack_buf_array(*,q,ie))
          elem_derived_eta_dot_dpdn_ptr = pack_elem_array(1,ie)
          elem_state_Qdp_ptr            = pack_elem_array(4,ie)
          elem_spheremp_ptr             = pack_elem_array(5,ie)
          elem_desc_reverse_south_ptr   = pack_elem_array(6,ie)
          elem_desc_reverse_north_ptr   = pack_elem_array(7,ie)
          elem_desc_reverse_east_ptr    = pack_elem_array(8,ie)
          elem_desc_reverse_west_ptr    = pack_elem_array(9,ie)
          elem_desc_putmapP_ptr         = pack_elem_array(10,ie)
          edge_buf_in_1_ptr = pack_buf_array(1,q,ie)
          edge_buf_in_2_ptr = pack_buf_array(2,q,ie)
          edge_buf_in_3_ptr = pack_buf_array(3,q,ie)
          edge_buf_in_4_ptr = pack_buf_array(4,q,ie)
          edge_buf_5_ptr    = pack_buf_array(5,q,ie)
          edge_buf_6_ptr    = pack_buf_array(6,q,ie)
          edge_buf_7_ptr    = pack_buf_array(7,q,ie)
          edge_buf_8_ptr    = pack_buf_array(8,q,ie)
          edge_buf_is_1_ptr = pack_buf_array(9,q,ie)
          edge_buf_is_2_ptr = pack_buf_array(10,q,ie)
          edge_buf_is_3_ptr = pack_buf_array(11,q,ie)
          edge_buf_is_4_ptr = pack_buf_array(12,q,ie)
          edge_buf_ie_1_ptr = pack_buf_array(13,q,ie)
          edge_buf_ie_2_ptr = pack_buf_array(14,q,ie)
          edge_buf_ie_3_ptr = pack_buf_array(15,q,ie)
          edge_buf_ie_4_ptr = pack_buf_array(16,q,ie)
          edge_buf_iw_1_ptr = pack_buf_array(17,q,ie)
          edge_buf_iw_2_ptr = pack_buf_array(18,q,ie)
          edge_buf_iw_3_ptr = pack_buf_array(19,q,ie)
          edge_buf_iw_4_ptr = pack_buf_array(20,q,ie)
          !$ACC DATA COPYIN(elem_desc_reverse_south, elem_desc_reverse_north, elem_desc_reverse_east, elem_desc_reverse_west,  elem_spheremp, elem_desc_putmapP) COPY(elem_derived_eta_dot_dpdn, edge_buf_5, edge_buf_6, edge_buf_7, edge_buf_8) COPYOUT(edge_buf_in_1, edge_buf_in_2, edge_buf_in_3, edge_buf_in_4, edge_buf_is_1, edge_buf_is_2, edge_buf_is_3, edge_buf_is_4, edge_buf_ie_1, edge_buf_ie_2, edge_buf_ie_3, edge_buf_ie_4, edge_buf_iw_1, edge_buf_iw_2, edge_buf_iw_3, edge_buf_iw_4)
          do k = 1 , my_nlev
            elem_derived_eta_dot_dpdn(:,:,k) = elem_spheremp(:,:) * elem_derived_eta_dot_dpdn(:,:,k)
          enddo
          call my_edgeVpack_acc(  elem_derived_eta_dot_dpdn(:,:,:) , my_nlev , my_nlev*my_qsize , elem_desc_putmapP(:), &
           elem_desc_reverse_north, elem_desc_reverse_south, elem_desc_reverse_east, elem_desc_reverse_west, &
           edge_buf_5, edge_buf_6, edge_buf_7, edge_buf_8, &
           edge_buf_in_1, edge_buf_in_2, edge_buf_in_3, edge_buf_in_4, &
           edge_buf_is_1, edge_buf_is_2, edge_buf_is_3, edge_buf_is_4, &
           edge_buf_ie_1, edge_buf_ie_2, edge_buf_ie_3, edge_buf_ie_4, &
           edge_buf_iw_1, edge_buf_iw_2, edge_buf_iw_3, edge_buf_iw_4)
        !$ACC END DATA
        !$ACC END DATA
      enddo
      !$ACC END PARALLEL LOOP
    endif


     if (DSSopt == 2) then
       !$ACC PARALLEL LOOP  
       do ie = nets, nete 
          !ie = calList(iee)
          q=my_qsize+1
          !$ACC DATA COPYIN(pack_elem_array(*,ie), pack_buf_array(*,q,ie))
          elem_derived_omega_p_ptr      = pack_elem_array(2,ie)
          elem_state_Qdp_ptr            = pack_elem_array(4,ie)
          elem_spheremp_ptr             = pack_elem_array(5,ie)
          elem_desc_reverse_south_ptr   = pack_elem_array(6,ie)
          elem_desc_reverse_north_ptr   = pack_elem_array(7,ie)
          elem_desc_reverse_east_ptr    = pack_elem_array(8,ie)
          elem_desc_reverse_west_ptr    = pack_elem_array(9,ie)
          elem_desc_putmapP_ptr         = pack_elem_array(10,ie)
          edge_buf_in_1_ptr = pack_buf_array(1,q,ie)
          edge_buf_in_2_ptr = pack_buf_array(2,q,ie)
          edge_buf_in_3_ptr = pack_buf_array(3,q,ie)
          edge_buf_in_4_ptr = pack_buf_array(4,q,ie)
          edge_buf_5_ptr    = pack_buf_array(5,q,ie)
          edge_buf_6_ptr    = pack_buf_array(6,q,ie)
          edge_buf_7_ptr    = pack_buf_array(7,q,ie)
          edge_buf_8_ptr    = pack_buf_array(8,q,ie)
          edge_buf_is_1_ptr = pack_buf_array(9,q,ie)
          edge_buf_is_2_ptr = pack_buf_array(10,q,ie)
          edge_buf_is_3_ptr = pack_buf_array(11,q,ie)
          edge_buf_is_4_ptr = pack_buf_array(12,q,ie)
          edge_buf_ie_1_ptr = pack_buf_array(13,q,ie)
          edge_buf_ie_2_ptr = pack_buf_array(14,q,ie)
          edge_buf_ie_3_ptr = pack_buf_array(15,q,ie)
          edge_buf_ie_4_ptr = pack_buf_array(16,q,ie)
          edge_buf_iw_1_ptr = pack_buf_array(17,q,ie)
          edge_buf_iw_2_ptr = pack_buf_array(18,q,ie)
          edge_buf_iw_3_ptr = pack_buf_array(19,q,ie)
          edge_buf_iw_4_ptr = pack_buf_array(20,q,ie)
          !$ACC DATA COPYIN(elem_desc_reverse_south, elem_desc_reverse_north, elem_desc_reverse_east, elem_desc_reverse_west,  elem_spheremp, elem_desc_putmapP) COPY(elem_derived_omega_p, edge_buf_5, edge_buf_6, edge_buf_7, edge_buf_8) COPYOUT(edge_buf_in_1, edge_buf_in_2, edge_buf_in_3, edge_buf_in_4, edge_buf_is_1, edge_buf_is_2, edge_buf_is_3, edge_buf_is_4, edge_buf_ie_1, edge_buf_ie_2, edge_buf_ie_3, edge_buf_ie_4, edge_buf_iw_1, edge_buf_iw_2, edge_buf_iw_3, edge_buf_iw_4)

          do k = 1 , my_nlev
            elem_derived_omega_p(:,:,k) = elem_spheremp(:,:) * elem_derived_omega_p(:,:,k)
          enddo
          call my_edgeVpack_acc(  elem_derived_omega_p(:,:,:) , my_nlev , my_nlev*my_qsize , elem_desc_putmapP, &
           elem_desc_reverse_north, elem_desc_reverse_south, elem_desc_reverse_east, elem_desc_reverse_west, &
           edge_buf_5, edge_buf_6, edge_buf_7, edge_buf_8, &
           edge_buf_in_1, edge_buf_in_2, edge_buf_in_3, edge_buf_in_4, &
           edge_buf_is_1, edge_buf_is_2, edge_buf_is_3, edge_buf_is_4, &
           edge_buf_ie_1, edge_buf_ie_2, edge_buf_ie_3, edge_buf_ie_4, &
           edge_buf_iw_1, edge_buf_iw_2, edge_buf_iw_3, edge_buf_iw_4)

        !$ACC END DATA
        !$ACC END DATA
      enddo
      !$ACC END PARALLEL LOOP
    endif


     if (DSSopt == 3) then
       !$ACC PARALLEL LOOP  
       do ie = nets, nete 
          !ie = calList(iee)
          q=my_qsize+1
          !$ACC DATA COPYIN(pack_elem_array(*,ie), pack_buf_array(*,q,ie))
          elem_derived_divdp_proj_ptr   = pack_elem_array(3,ie)
          elem_state_Qdp_ptr            = pack_elem_array(4,ie)
          elem_spheremp_ptr             = pack_elem_array(5,ie)
          elem_desc_reverse_south_ptr   = pack_elem_array(6,ie)
          elem_desc_reverse_north_ptr   = pack_elem_array(7,ie)
          elem_desc_reverse_east_ptr    = pack_elem_array(8,ie)
          elem_desc_reverse_west_ptr    = pack_elem_array(9,ie)
          elem_desc_putmapP_ptr         = pack_elem_array(10,ie)
          edge_buf_in_1_ptr = pack_buf_array(1,q,ie)
          edge_buf_in_2_ptr = pack_buf_array(2,q,ie)
          edge_buf_in_3_ptr = pack_buf_array(3,q,ie)
          edge_buf_in_4_ptr = pack_buf_array(4,q,ie)
          edge_buf_5_ptr    = pack_buf_array(5,q,ie)
          edge_buf_6_ptr    = pack_buf_array(6,q,ie)
          edge_buf_7_ptr    = pack_buf_array(7,q,ie)
          edge_buf_8_ptr    = pack_buf_array(8,q,ie)
          edge_buf_is_1_ptr = pack_buf_array(9,q,ie)
          edge_buf_is_2_ptr = pack_buf_array(10,q,ie)
          edge_buf_is_3_ptr = pack_buf_array(11,q,ie)
          edge_buf_is_4_ptr = pack_buf_array(12,q,ie)
          edge_buf_ie_1_ptr = pack_buf_array(13,q,ie)
          edge_buf_ie_2_ptr = pack_buf_array(14,q,ie)
          edge_buf_ie_3_ptr = pack_buf_array(15,q,ie)
          edge_buf_ie_4_ptr = pack_buf_array(16,q,ie)
          edge_buf_iw_1_ptr = pack_buf_array(17,q,ie)
          edge_buf_iw_2_ptr = pack_buf_array(18,q,ie)
          edge_buf_iw_3_ptr = pack_buf_array(19,q,ie)
          edge_buf_iw_4_ptr = pack_buf_array(20,q,ie)
          !$ACC DATA COPYIN(elem_desc_reverse_south, elem_desc_reverse_north, elem_desc_reverse_east, elem_desc_reverse_west,  elem_spheremp, elem_desc_putmapP) COPY(elem_derived_divdp_proj, edge_buf_5, edge_buf_6, edge_buf_7, edge_buf_8) COPYOUT(edge_buf_in_1, edge_buf_in_2, edge_buf_in_3, edge_buf_in_4, edge_buf_is_1, edge_buf_is_2, edge_buf_is_3, edge_buf_is_4, edge_buf_ie_1, edge_buf_ie_2, edge_buf_ie_3, edge_buf_ie_4, edge_buf_iw_1, edge_buf_iw_2, edge_buf_iw_3, edge_buf_iw_4)

          do k = 1 , my_nlev
            elem_derived_divdp_proj(:,:,k) = elem_spheremp(:,:) * elem_derived_divdp_proj(:,:,k)
          enddo
          call my_edgeVpack_acc(  elem_derived_divdp_proj(:,:,:) , my_nlev , my_nlev*my_qsize , elem_desc_putmapP, &
           elem_desc_reverse_north, elem_desc_reverse_south, elem_desc_reverse_east, elem_desc_reverse_west, &
           edge_buf_5, edge_buf_6, edge_buf_7, edge_buf_8, &
           edge_buf_in_1, edge_buf_in_2, edge_buf_in_3, edge_buf_in_4, &
           edge_buf_is_1, edge_buf_is_2, edge_buf_is_3, edge_buf_is_4, &
           edge_buf_ie_1, edge_buf_ie_2, edge_buf_ie_3, edge_buf_ie_4, &
           edge_buf_iw_1, edge_buf_iw_2, edge_buf_iw_3, edge_buf_iw_4)
        !$ACC END DATA
        !$ACC END DATA
      enddo
      !$ACC END PARALLEL LOOP
    endif

  end subroutine my_packCode_acc

  subroutine my_euler_initCal(my_nets,my_nete,my_qsize,my_nlev,my_hvcoord,my_qmin,my_qmax,Qtens_biharmonic,my_rhs_multiplier,my_dt,my_np,n0_qdp,initCal_array)

    integer              , intent(in   )         :: my_nets
    integer              , intent(in   )         :: my_nete
    integer              , intent(in   )         :: my_qsize
    integer              , intent(in   )         :: my_nlev
    integer              , intent(in   )         :: my_np
    integer              , intent(in   )         :: my_rhs_multiplier
    integer              , intent(in   )         :: n0_qdp
    real (kind=real_kind), intent(in   )         :: my_dt
    type (hvcoord_t)     , intent(in   )         :: my_hvcoord
    integer(kind=8)      , intent(in   ),    dimension(4,my_nets:my_nete) :: initCal_array
    real(kind=8)         , intent(inout),    dimension(my_nlev,my_qsize,my_nets:my_nete)        :: my_qmin
    real(kind=8)         , intent(inout),    dimension(my_nlev,my_qsize,my_nets:my_nete)        :: my_qmax
    real(kind=real_kind) , intent(inout),    dimension(my_np,my_np  ,my_nlev,my_qsize,my_nets:my_nete) :: Qtens_biharmonic
    
    real(kind=8)         , dimension(my_np,my_np,my_nlev) :: elem_derived_dp
    pointer(elem_derived_dp_ptr, elem_derived_dp)
 
    real(kind=8)         , dimension(my_np,my_np,my_nlev) :: elem_derived_divdp_proj
    pointer(elem_derived_divdp_proj_ptr, elem_derived_divdp_proj)

    real(kind=8)         , dimension(my_np,my_np,my_nlev,my_qsize) :: elem_state_Qdp
    pointer(elem_state_Qdp_ptr, elem_state_Qdp)

    real(kind=real_kind) , dimension(my_np,my_np,my_nlev)   :: elem_derived_dpdiss_ave
    pointer(elem_derived_dpdiss_ave_ptr, elem_derived_dpdiss_ave)

    real(kind=8)         , dimension(my_nlev+1) :: hvcoord_hyai
    pointer(hvcoord_hyai_ptr, hvcoord_hyai)

    real(kind=8)         , dimension(my_nlev+1) :: hvcoord_hybi
    pointer(hvcoord_hybi_ptr, hvcoord_hybi)

    real(kind=8) :: hvcoord_ps0

    real(kind=real_kind), dimension(my_np,my_np  ,my_nlev                ) :: dp
    integer :: ie,q,k


    hvcoord_ps0      = my_hvcoord%ps0
    hvcoord_hyai_ptr = loc(my_hvcoord%hyai)
    hvcoord_hybi_ptr = loc(my_hvcoord%hybi)

    !$ACC PARALLEL LOOP  collapse(2) local(dp) copyin(initCal_array, hvcoord_hyai, hvcoord_hybi) copyout(my_qmin,my_qmax) annotate(entire(hvcoord_hyai, hvcoord_hybi))
    do ie = my_nets , my_nete 
      do q = 1 , my_qsize
        elem_derived_dp_ptr         = initCal_array(1,ie)
        elem_derived_divdp_proj_ptr = initCal_array(2,ie)
        elem_state_Qdp_ptr          = initCal_array(3,ie)
        elem_derived_dpdiss_ave_ptr = initCal_array(4,ie)
        !$ACC DATA copyin(elem_derived_dp(*,*,1:64), elem_derived_divdp_proj(*,*,1:64), elem_state_Qdp(*,*,1:64,q), elem_derived_dpdiss_ave(*,*,1:64)) copyout(Qtens_biharmonic(*,*,1:64,q,ie)) 
        do k = 1 , 64   
          dp(:,:,k) = elem_derived_dp(:,:,k) - my_rhs_multiplier*my_dt*elem_derived_divdp_proj(:,:,k)
          Qtens_biharmonic(:,:,k,q,ie) = elem_state_Qdp(:,:,k,q)/dp(:,:,k)
          my_qmin(k,q,ie)=minval(Qtens_biharmonic(:,:,k,q,ie))
          my_qmax(k,q,ie)=maxval(Qtens_biharmonic(:,:,k,q,ie))
          my_qmin(k,q,ie)=max(my_qmin(k,q,ie),0d0)
          if ( my_rhs_multiplier == 2) then
            Qtens_biharmonic(:,:,k,q,ie) = Qtens_biharmonic(:,:,k,q,ie) * elem_derived_dpdiss_ave(:,:,k) / &
              ( (hvcoord_hyai(k+1) - hvcoord_hyai(k)) * hvcoord_ps0 + &
                (hvcoord_hybi(k+1) - hvcoord_hybi(k) )* hvcoord_ps0 )
          endif
        enddo
        !$ACC END DATA
      enddo
    enddo
    !$ACC END PARALLEL LOOP
    
    !$ACC PARALLEL LOOP  collapse(2) local(dp) copyin(initCal_array, hvcoord_hyai, hvcoord_hybi) copyout(my_qmin,my_qmax) annotate(entire(hvcoord_hyai, hvcoord_hybi))
    do ie = my_nets , my_nete 
      do q = 1 , my_qsize
        elem_derived_dp_ptr         = initCal_array(1,ie)
        elem_derived_divdp_proj_ptr = initCal_array(2,ie)
        elem_state_Qdp_ptr          = initCal_array(3,ie)
        elem_derived_dpdiss_ave_ptr = initCal_array(4,ie)
        !$ACC DATA copyin(elem_derived_dp(*,*,65:128), elem_derived_divdp_proj(*,*,65:128), elem_state_Qdp(*,*,65:128,q), elem_derived_dpdiss_ave(*,*,65:128)) copyout(Qtens_biharmonic(*,*,65:128,q,ie)) 
        do k = 65 , 128   
          dp(:,:,k) = elem_derived_dp(:,:,k) - my_rhs_multiplier*my_dt*elem_derived_divdp_proj(:,:,k)
          Qtens_biharmonic(:,:,k,q,ie) = elem_state_Qdp(:,:,k,q)/dp(:,:,k)
          my_qmin(k,q,ie)=minval(Qtens_biharmonic(:,:,k,q,ie))
          my_qmax(k,q,ie)=maxval(Qtens_biharmonic(:,:,k,q,ie))
          my_qmin(k,q,ie)=max(my_qmin(k,q,ie),0d0)
          if ( my_rhs_multiplier == 2) then
            Qtens_biharmonic(:,:,k,q,ie) = Qtens_biharmonic(:,:,k,q,ie) * elem_derived_dpdiss_ave(:,:,k) / &
              ( (hvcoord_hyai(k+1) - hvcoord_hyai(k)) * hvcoord_ps0 + &
                (hvcoord_hybi(k+1) - hvcoord_hybi(k) )* hvcoord_ps0 )
          endif
        enddo
        !$ACC END DATA
      enddo
    enddo
    !$ACC END PARALLEL LOOP
  end subroutine my_euler_initCal

  subroutine euler_step( np1_qdp , n0_qdp , dt , elem , hvcoord , hybrid , deriv , nets , nete , DSSopt , rhs_multiplier )
  ! ===================================
  ! This routine is the basic foward
  ! euler component used to construct RK SSP methods
  !
  !           u(np1) = u(n0) + dt2*DSS[ RHS(u(n0)) ]
  !
  ! n0 can be the same as np1.
  !
  ! DSSopt = DSSeta or DSSomega:   also DSS eta_dot_dpdn or omega
  !
  ! ===================================
  use kinds          , only : real_kind
  use dimensions_mod , only : np, npdg, nlev, max_corner_elem
  use hybrid_mod     , only : hybrid_t
  use element_mod    , only : element_t
  use derivative_mod , only : derivative_t, divergence_sphere, gradient_sphere, vorticity_sphere
  use edge_mod       , only : edgevpack, edgevunpack
  use bndry_mod      , only : bndry_exchangev
  use hybvcoord_mod  , only : hvcoord_t
  use control_mod, only : hypervis_scaling, north, south, east, west, neast, nwest, seast, swest, which_vlaplace
  use schedule_mod, only : schedule_t, cycle_t, schedule
#ifdef _ACCEL
  use cuda_mod, only: euler_step_cuda
#endif
  implicit none
  integer              , intent(in   )         :: np1_qdp, n0_qdp
  real (kind=real_kind), intent(in   )         :: dt
  type (element_t)     , intent(inout), target :: elem(:)
  type (hvcoord_t)     , intent(in   )         :: hvcoord
  type (hybrid_t)      , intent(in   )         :: hybrid
  type (derivative_t)  , intent(in   )         :: deriv
  integer              , intent(in   )         :: nets
  integer              , intent(in   )         :: nete
  integer              , intent(in   )         :: DSSopt
  integer              , intent(in   )         :: rhs_multiplier

  ! local
  real(kind=real_kind), dimension(np,np                       ) :: divdp, dpdiss
  real(kind=real_kind), dimension(np,np,2                     ) :: gradQ
  real(kind=real_kind), dimension(np,np,2,nlev                ) :: Vstar
  real(kind=real_kind), dimension(np,np  ,nlev                ) :: Qtens
  real(kind=real_kind), dimension(np,np  ,nlev                ) :: dp,dp_star
  real(kind=real_kind), dimension(np,np  ,nlev,qsize,nets:nete) :: Qtens_biharmonic
  real(kind=real_kind), pointer, dimension(:,:,:)               :: DSSvar
  real(kind=real_kind) :: dp0
  integer :: ie,q,i,j,k, ieInner, ieOuter
  integer :: rhs_viss = 0
  integer(kind=8) :: count_start, count_stop, count_max, count_rate

  ! local added by conghui
  real(kind=real_kind), dimension(np,np)                      ::  qdp_val
  real(kind=real_kind), dimension(np,np,nets:nete)            ::  elem_spheremp
  real(kind=real_kind), dimension(np,np,nlev,nets:nete)   :: elem_derived_dpdiss_ave
  integer(kind=8), dimension(11,nets:nete) :: elem_array
  integer(kind=8), dimension(18,nets:nete)      :: pack_elem_array
  integer(kind=8), dimension(20,qsize+1,nets:nete) :: pack_buf_array
  integer(kind=8), dimension(4,nets:nete) :: initCal_array

  real,    dimension(nlev*(qsize+1),(nete-nets+1)*20):: buf1,buf2
  integer, dimension((nete-nets+1)*20) :: buffPosToIe
  integer, dimension(nete-nets+1) :: innerList, outerList
  logical, dimension(nete-nets+1) :: outFlag
  integer :: outerLen, innerLen
  type (Schedule_t),pointer                     :: pSchedule
  type (Cycle_t),pointer                        :: pCycle
  integer :: nSendCycles, nRecvCycles, iCycle, iptr, buffLen, dest, source, length, nlyr
  integer :: rank, ierr, tag
  include "mpif.h"
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)



  if ( npdg > 0 ) then
    call euler_step_dg( np1_qdp , n0_qdp , dt , elem , hvcoord , hybrid , deriv , nets , nete , DSSopt , rhs_multiplier )
    return
  endif
#ifdef _ACCEL
  call euler_step_cuda( np1_qdp , n0_qdp , dt , elem , hvcoord , hybrid , deriv , nets , nete , DSSopt , rhs_multiplier )
  return
#endif
! call t_barrierf('sync_euler_step', hybrid%par%comm)
  call t_startf('euler_step')

  ! added by conghui
  ! AOS TO SOA
  do ie = nets, nete
    elem_derived_dpdiss_ave(:,:,:,ie) = elem(ie)%derived%dpdiss_ave(:,:,:)
    elem_spheremp(:,:,ie)                    = elem(ie)%spheremp(:,:)
  end do

  rhs_viss = 0
    call t_startf('euler init cal')
    do ie = nets, nete
        initCal_array(1,ie) = loc(elem(ie)%derived%dp) 
        initCal_array(2,ie) = loc(elem(ie)%derived%divdp_proj)
        initCal_array(3,ie) = loc(elem(ie)%state%Qdp(:,:,:,:,n0_qdp))
        initCal_array(4,ie) = loc(elem(ie)%derived%dpdiss_ave)
    enddo
    call my_euler_initCal(nets,nete,qsize,nlev,hvcoord,qmin,qmax,Qtens_biharmonic,rhs_multiplier,dt,np,n0_qdp,initCal_array)
    call t_stopf('euler init cal')

    
    call t_startf('euler_init cal2')
    ! compute element qmin/qmax
    if ( rhs_multiplier == 0 ) then
           call neighbor_minmax(elem,hybrid,edgeAdvQ2,nets,nete,qmin(:,:,nets:nete),qmax(:,:,nets:nete))
    endif

    if ( rhs_multiplier == 1 ) then
      call cal_qmin_qmax_rhs_multiplier_1(nets, nete, qsize, nlev, &
                                            np, Qtens_biharmonic, qmin, qmax)
    endif

    if ( rhs_multiplier == 2 ) then
      rhs_viss = 3
        call biharmonic_wk_scalar_minmax( elem , qtens_biharmonic , deriv , edgeAdvQ3 , hybrid , &
                                          nets , nete , qmin(:,:,nets:nete) , qmax(:,:,nets:nete) )
      call cal_qtens_biharmonic(nets, nete, qsize, nlev, &
                                            np, Qtens_biharmonic, &
                                            rhs_viss, dt, nu_q, hvcoord%hyai, hvcoord%hybi, &
                                            hvcoord%ps0, elem_spheremp)
    endif
 

  call t_stopf('euler_init cal2')
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !   2D Advection step
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  do ie = nets, nete
     buffPosToIe(elem(ie)%desc%putmapP(north)+1) = ie
     buffPosToIe(elem(ie)%desc%putmapP(north)+2) = ie  
     buffPosToIe(elem(ie)%desc%putmapP(north)+3) = ie 
     buffPosToIe(elem(ie)%desc%putmapP(north)+4) = ie 
     buffPosToIe(elem(ie)%desc%putmapP(5)+1    ) = ie  
     buffPosToIe(elem(ie)%desc%putmapP(6)+1    ) = ie
     buffPosToIe(elem(ie)%desc%putmapP(7)+1    ) = ie
     buffPosToIe(elem(ie)%desc%putmapP(8)+1    ) = ie
     buffPosToIe(elem(ie)%desc%putmapP(south)+1) = ie 
     buffPosToIe(elem(ie)%desc%putmapP(south)+2) = ie 
     buffPosToIe(elem(ie)%desc%putmapP(south)+3) = ie 
     buffPosToIe(elem(ie)%desc%putmapP(south)+4) = ie 
     buffPosToIe(elem(ie)%desc%putmapP(east)+1 ) = ie
     buffPosToIe(elem(ie)%desc%putmapP(east)+2 ) = ie
     buffPosToIe(elem(ie)%desc%putmapP(east)+3 ) = ie
     buffPosToIe(elem(ie)%desc%putmapP(east)+4 ) = ie
     buffPosToIe(elem(ie)%desc%putmapP(west)+1 ) = ie
     buffPosToIe(elem(ie)%desc%putmapP(west)+2 ) = ie
     buffPosToIe(elem(ie)%desc%putmapP(west)+3 ) = ie
     buffPosToIe(elem(ie)%desc%putmapP(west)+4 ) = ie
  end do
  !if (rank==149) write(*,*) "buffToIe:",buffPosToIe

  buffLen = 0
  outFlag(:) = .false.
  pSchedule => Schedule(1)
  nSendCycles = pSchedule%nSendCycles
  !if (rank==194) write(*,*) "ljf cc", nSendCycles
  do iCycle = 1, nSendCycles
      pCycle      => pSchedule%SendCycle(icycle)
      iptr            = pCycle%ptrP
      buffLen = buffLen + pCycle%lengthP
      do i = 1, pCycle%lengthP
          outFlag(buffPosToIe(iptr+i-1)) = .true.
      enddo
  end do
  !if (rank==149) then
  !   do iCycle = 1, nSendCycles
  !        pCycle      => pSchedule%SendCycle(icycle)
  !        write(*,*) "iCycle",iCycle,"dest",pCycle%dest-1,"source",pCycle%source-1
  !   enddo
  !endif
  innerLen = 0
  outerLen = 0
  do ie = nets, nete
     if (outFlag(ie) == .true.) then
         outerLen = outerLen + 1
         outerList(outerLen) = ie
     else 
         innerLen = innerLen + 1
         innerList(innerLen) = ie
     end if
  end do
 !if (innerLen /=0) then
 !   write(*,*) "not empty rank",rank,"innerLen",innerLen
 !endif
! if (rank==1) then
!    write(*,*) "rank =",rank,"outerLen",outerLen
!    write(*,*) "outerList",outerList
! endif
! if (rank /= 300) then
     outerLen = nete - nets + 1
     innerLen = 0
     do i = nets, nete
         outerList(i-nets+1) = i
     enddo
! endif


  !call system_clock(count_start, count_rate, count_max)
  call t_startf('my_euler step')
  do ie = nets , nete
    elem_array(1,ie)  = loc(elem(ie)%derived%dp)
    elem_array(2,ie)  = loc(elem(ie)%derived%divdp_proj)
    elem_array(3,ie)  = loc(elem(ie)%derived%divdp)
    elem_array(4,ie)  = loc(elem(ie)%derived%dpdiss_biharmonic)
    elem_array(5,ie)  = loc(elem(ie)%derived%vn0)
    elem_array(6,ie)  = loc(elem(ie)%Dinv)
    elem_array(7,ie)  = loc(elem(ie)%metdet)
    elem_array(8,ie)  = loc(elem(ie)%rmetdet)
    elem_array(9,ie)  = loc(elem(ie)%spheremp)
    elem_array(10,ie) = loc(elem(ie)%state%Qdp(:,:,:,:,n0_qdp))
    elem_array(11,ie) = loc(elem(ie)%state%Qdp(:,:,:,:,np1_qdp))
  end do
  if (outerLen > 0) then
      call my_euler_step_acc(nets, nete, rhs_multiplier, rhs_viss, &
               dt, rrearth, nu_p, nu_q, &
               deriv%dvv, Qtens_biharmonic, &
               qmin,qmax,elem_array, n0_qdp, np1_qdp, qsize, outerList, outerLen)
  endif
  
  !write(*,*) "rank =",rank,"buffLen =",buffLen,"totalLen =", (nete-nets+1)*20, "propotion =", buffLen/((nete-nets+1)*20)
  !call system_clock(count_stop, count_rate, count_max)
  !write(*,*) 'eluer computing count = ', (count_stop - count_start), 'rate=', count_rate

  !call system_clock(count_start, count_rate, count_max)
  call t_stopf('my_euler step')
  
  call t_startf('my_euler pack')
  do ie=nets,nete
    pack_elem_array(1,ie) = loc(elem(ie)%derived%eta_dot_dpdn)
    pack_elem_array(2,ie) = loc(elem(ie)%derived%omega_p)
    pack_elem_array(3,ie) = loc(elem(ie)%derived%divdp_proj)
    pack_elem_array(4,ie) = loc(elem(ie)%state%Qdp(:,:,:,:,np1_qdp))
    pack_elem_array(5,ie) = loc(elem(ie)%spheremp)
    pack_elem_array(6,ie) = loc(elem(ie)%desc%reverse(south))
    pack_elem_array(7,ie) = loc(elem(ie)%desc%reverse(north))
    pack_elem_array(8,ie) = loc(elem(ie)%desc%reverse(east))
    pack_elem_array(9,ie) = loc(elem(ie)%desc%reverse(west))
    pack_elem_array(10,ie)= loc(elem(ie)%desc%putmapP)
  enddo
  do ie = nets, nete
    do q = 1, qsize+1
      pack_buf_array(1,q,ie) = loc(edgeAdv_p1%buf((q-1)*nlev+1,elem(ie)%desc%putmapP(north)+1))
      pack_buf_array(2,q,ie) = loc(edgeAdv_p1%buf((q-1)*nlev+1,elem(ie)%desc%putmapP(north)+2))
      pack_buf_array(3,q,ie) = loc(edgeAdv_p1%buf((q-1)*nlev+1,elem(ie)%desc%putmapP(north)+3))
      pack_buf_array(4,q,ie) = loc(edgeAdv_p1%buf((q-1)*nlev+1,elem(ie)%desc%putmapP(north)+4))
      pack_buf_array(5,q,ie) = loc(edgeAdv_p1%buf((q-1)*nlev+1,elem(ie)%desc%putmapP(5)+1))
      pack_buf_array(6,q,ie) = loc(edgeAdv_p1%buf((q-1)*nlev+1,elem(ie)%desc%putmapP(6)+1))
      pack_buf_array(7,q,ie) = loc(edgeAdv_p1%buf((q-1)*nlev+1,elem(ie)%desc%putmapP(7)+1))
      pack_buf_array(8,q,ie) = loc(edgeAdv_p1%buf((q-1)*nlev+1,elem(ie)%desc%putmapP(8)+1))
      pack_buf_array(9,q,ie) = loc(edgeAdv_p1%buf((q-1)*nlev+1,elem(ie)%desc%putmapP(south)+1))
      pack_buf_array(10,q,ie) = loc(edgeAdv_p1%buf((q-1)*nlev+1,elem(ie)%desc%putmapP(south)+2))
      pack_buf_array(11,q,ie) = loc(edgeAdv_p1%buf((q-1)*nlev+1,elem(ie)%desc%putmapP(south)+3))
      pack_buf_array(12,q,ie) = loc(edgeAdv_p1%buf((q-1)*nlev+1,elem(ie)%desc%putmapP(south)+4))
      pack_buf_array(13,q,ie) = loc(edgeAdv_p1%buf((q-1)*nlev+1,elem(ie)%desc%putmapP(east)+1))
      pack_buf_array(14,q,ie) = loc(edgeAdv_p1%buf((q-1)*nlev+1,elem(ie)%desc%putmapP(east)+2))
      pack_buf_array(15,q,ie) = loc(edgeAdv_p1%buf((q-1)*nlev+1,elem(ie)%desc%putmapP(east)+3))
      pack_buf_array(16,q,ie) = loc(edgeAdv_p1%buf((q-1)*nlev+1,elem(ie)%desc%putmapP(east)+4))
      pack_buf_array(17,q,ie) = loc(edgeAdv_p1%buf((q-1)*nlev+1,elem(ie)%desc%putmapP(west)+1))
      pack_buf_array(18,q,ie) = loc(edgeAdv_p1%buf((q-1)*nlev+1,elem(ie)%desc%putmapP(west)+2))
      pack_buf_array(19,q,ie) = loc(edgeAdv_p1%buf((q-1)*nlev+1,elem(ie)%desc%putmapP(west)+3))
      pack_buf_array(20,q,ie) = loc(edgeAdv_p1%buf((q-1)*nlev+1,elem(ie)%desc%putmapP(west)+4))
    enddo
  enddo

  if (outerLen /=0) call my_packCode_acc(np1_qdp,nets,nete,DSSopt,pack_elem_array,pack_buf_array, qsize, outerLen, outerList, .true.)
  
  
  


  !call system_clock(count_stop, count_rate, count_max)
  !write(*,*) "edge pack count=",count_stop-count_start

  call t_stopf('my_euler pack')

  call t_startf('my_euler bndry')
  !call system_clock(count_start, count_rate, count_max)
    !call bndry_exchangeV( hybrid , edgeAdv_p1 )
  !--------------bndry_exchangeV---------------------

       nlyr = edgeAdv_p1%nlyr

       nSendCycles = pSchedule%nSendCycles
       nRecvCycles = pSchedule%nRecvCycles

       !==================================================
       !  Fire off the sends
       !==================================================
       do icycle=1,nSendCycles
          pCycle      => pSchedule%SendCycle(icycle)
          dest            = pCycle%dest - 1
          length      = nlyr * pCycle%lengthP
          tag             = pCycle%tag
          iptr            = pCycle%ptrP
          !DBG if(Debug) print *,'bndry_exchangeV: MPI_Isend: DEST:',dest,'LENGTH:',length,'TAG: ',tag
          call MPI_Isend(edgeAdv_p1%buf(1,iptr),length,MPIreal_t,dest,tag,hybrid%par%comm,Srequest(icycle),ierr)
          !if(ierr .ne. MPI_SUCCESS) then
          !   errorcode=ierr
          !   call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
          !   print *,'bndry_exchangeV: Error after call to MPI_Isend: ',errorstring
          !endif
       end do    ! icycle

       !==================================================
       !  Post the Receives 
       !==================================================
       do icycle=1,nRecvCycles
          pCycle         => pSchedule%RecvCycle(icycle)
          source          = pCycle%source - 1
          length      = nlyr * pCycle%lengthP
          tag             = pCycle%tag
          iptr            = pCycle%ptrP
          !DBG if(Debug) print *,'bndry_exchangeV: MPI_Irecv: SRC:',source,'LENGTH:',length,'TAG: ',tag
          call MPI_Irecv(edgeAdv_p1%receive(1,iptr),length,MPIreal_t, &
               source,tag,hybrid%par%comm,Rrequest(icycle),ierr)
          !if(ierr .ne. MPI_SUCCESS) then
          !   errorcode=ierr
          !   call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
          !   print *,'bndry_exchangeV: Error after call to MPI_Irecv: ',errorstring
          !endif
       end do    ! icycle

       if (innerLen > 0) then
          !write(*,*) "inner is not empty", rank, innerLen             
          call my_euler_step_acc(nets, nete, rhs_multiplier, rhs_viss, &
               dt, rrearth, nu_p, nu_q, &
               deriv%dvv, Qtens_biharmonic, &
               qmin,qmax,elem_array, n0_qdp, np1_qdp, qsize, innerList, innerLen)
       endif    

       !==================================================
       !  Wait for all the receives to complete
       !==================================================

       call MPI_Waitall(nSendCycles,Srequest,status,ierr)
       call MPI_Waitall(nRecvCycles,Rrequest,status,ierr)

       do icycle=1,nRecvCycles
          pCycle         => pSchedule%RecvCycle(icycle)
          length             = pCycle%lengthP
          iptr            = pCycle%ptrP
          do i=0,length-1
             edgeAdv_p1%buf(1:nlyr,iptr+i) = edgeAdv_p1%receive(1:nlyr,iptr+i)
          enddo
       end do   ! icycle
       if (innerLen > 0) call my_packCode_acc(np1_qdp,nets,nete,DSSopt,pack_elem_array,pack_buf_array, qsize, innerLen, innerList, .true.)
  !-----------------------------------
  call t_stopf('my_euler bndry')
  
  
  !call system_clock(count_stop, count_rate, count_max)
  !write(*,*) 'bndry exchange count = ', count_stop - count_start

  !call system_clock(count_start, count_rate, count_max)
  
  call t_startf('my_euler unpack')
  !outerLen = nete - nets + 1
  !do i = 1, outerLen
  !    outerList(i) = nets + i - 1
  !enddo
  if (outerLen > 0) then 
      call  my_unpack_acc(nets, nete, edgeAdv_p1%nlyr, edgeAdv_p1%nbuf, &
                        edgeAdv_p1%buf, south, east, north, west, nlev, &
                        swest, max_corner_elem, elem, DSSopt, np1_qdp, outerLen, outerList)
  endif
  if (innerLen > 0) then
      call  my_unpack_acc(nets, nete, edgeAdv_p1%nlyr, edgeAdv_p1%nbuf, &
                        edgeAdv_p1%buf, south, east, north, west, nlev, &
                        swest, max_corner_elem, elem, DSSopt, np1_qdp, innerLen, innerList)

  endif
  !call system_clock(count_stop, count_rate, count_max)
  !write(*,*) 'edge unpack count = ', count_stop - count_start
  call t_stopf('my_euler unpack')
  
  call t_stopf('euler_step')
  end subroutine euler_step

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

  subroutine euler_step_dg(np1_qdp, n0_qdp, dt,elem,hvcoord,hybrid,deriv,nets,nete,&
      DSSopt,rhs_multiplier)
  ! ===================================
  ! This routine is the basic foward
  ! euler component used to construct RK SSP methods
  !
  !           u(np1) = u(n0) + dt2*DSS[ RHS(u(n0)) ]
  !
  ! n0 can be the same as np1.
  !
  ! DSSopt = DSSeta or DSSomega:   also DSS eta_dot_dpdn or omega
  !
  ! ===================================
  use kinds, only : real_kind
  use dimensions_mod, only : np, npdg, nlev
  use hybrid_mod, only : hybrid_t
  use element_mod, only : element_t
  use derivative_mod, only : derivative_t, divergence_sphere_wk, edge_flux_u_cg, gll_to_dgmodal, dgmodal_to_gll
  use edge_mod, only : edgevpack, edgevunpack, edgedgvunpack
  use bndry_mod, only : bndry_exchangev
  use hybvcoord_mod, only : hvcoord_t

  implicit none
  integer :: np1_qdp, n0_qdp, nets, nete, DSSopt, rhs_multiplier
  real (kind=real_kind), intent(in)  :: dt

  type (hvcoord_t)     , intent(in) :: hvcoord
  type (hybrid_t)      , intent(in) :: hybrid
  type (element_t)     , intent(inout), target :: elem(:)
  type (derivative_t)  , intent(in) :: deriv

  ! local
  real (kind=real_kind), dimension(np,np)    :: divdp
  real (kind=real_kind), dimension(npdg,npdg)    :: pshat
  real (kind=real_kind), dimension(0:np+1,0:np+1,nlev,qsize)    :: qedges
  real (kind=real_kind), dimension(np,np,2)    :: vtemp
  real(kind=real_kind), dimension(np,np,nlev) :: dp,dp_star
  real(kind=real_kind), dimension(np,np,2,nlev) :: Vstar
  real (kind=real_kind), pointer, dimension(:,:,:)   :: DSSvar
! nelemd

  real(kind=real_kind) :: dp0
  integer :: ie,q,i,j,k
  integer :: rhs_viss=0

  call t_barrierf('sync_euler_step_dg', hybrid%par%comm)
  call t_startf('euler_step_dg')


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !   compute Q min/max values for lim8
  !   compute biharmonic mixing term f
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  rhs_viss=0
  if (limiter_option == 8) then
     call abortmp('limiter_opiton=8 not supported for dg advection')
     ! todo:  we need to track a 'dg' mass, and use that to back out Q
     ! then compute Qmin/Qmax here
  endif



  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !   2D Advection step
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  write(*,*) "columns ljf:", nete-nets
  do ie=nets,nete

     ! note: eta_dot_dpdn is actually dimension nlev+1, but nlev+1 data is
     ! all zero so we only have to DSS 1:nlev
     if ( DSSopt == DSSeta) DSSvar => elem(ie)%derived%eta_dot_dpdn(:,:,:)
     if ( DSSopt == DSSomega) DSSvar => elem(ie)%derived%omega_p(:,:,:)
     if ( DSSopt == DSSdiv_vdp_ave) DSSvar => elem(ie)%derived%divdp_proj(:,:,:)

     if(DSSopt==DSSno_var)then
	call edgeVpack(edgeAdv,elem(ie)%state%Qdp(:,:,:,:,n0_qdp),nlev*qsize,0,elem(ie)%desc)
     else
	call edgeVpack(edgeAdv_p1,elem(ie)%state%Qdp(:,:,:,:,n0_qdp),nlev*qsize,0,elem(ie)%desc)
	! also DSS extra field
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k)
#endif
	do k=1,nlev
	    DSSvar(:,:,k) = elem(ie)%spheremp(:,:)*DSSvar(:,:,k)
	enddo
	call edgeVpack(edgeAdv_p1,DSSvar(:,:,1:nlev),nlev,nlev*qsize,elem(ie)%desc)
     endif

  end do

  if(DSSopt==DSSno_var)then
     call bndry_exchangeV(hybrid,edgeAdv)
  else
     call bndry_exchangeV(hybrid,edgeAdv_p1)
  endif

  do ie=nets,nete

     if ( DSSopt == DSSeta) DSSvar => elem(ie)%derived%eta_dot_dpdn(:,:,:)
     if ( DSSopt == DSSomega) DSSvar => elem(ie)%derived%omega_p(:,:,:)
     if ( DSSopt == DSSdiv_vdp_ave) DSSvar => elem(ie)%derived%divdp_proj(:,:,:)

     if(DSSopt==DSSno_var)then
	call edgeDGVunpack(edgeAdv,qedges,nlev*qsize,0,elem(ie)%desc)
     else
	call edgeDGVunpack(edgeAdv_p1,qedges,nlev*qsize,0,elem(ie)%desc)
	call edgeVunpack(edgeAdv_p1,DSSvar(:,:,1:nlev),nlev,qsize*nlev,elem(ie)%desc)
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k,q)
#endif
	do k=1,nlev
	  DSSvar(:,:,k)=DSSvar(:,:,k)*elem(ie)%rspheremp(:,:)
	enddo
     endif

     ! compute flux and advection term
     do k=1,nlev
        dp(:,:,k) = elem(ie)%derived%dp(:,:,k) - &
             rhs_multiplier*dt*elem(ie)%derived%divdp_proj(:,:,k)
        Vstar(:,:,1,k) = elem(ie)%derived%vn0(:,:,1,k)/dp(:,:,k)
        Vstar(:,:,2,k) = elem(ie)%derived%vn0(:,:,2,k)/dp(:,:,k)
     enddo

     do q=1,qsize
        do k=1,nlev
           vtemp(:,:,1)=elem(ie)%state%Qdp(:,:,k,q,n0_qdp)*Vstar(:,:,1,k)
           vtemp(:,:,2)=elem(ie)%state%Qdp(:,:,k,q,n0_qdp)*Vstar(:,:,2,k)

           divdp = divergence_sphere_wk(vtemp,deriv,elem(ie)) + &
                edge_flux_u_cg( Vstar(:,:,:,k), elem(ie)%state%Qdp(:,:,k,q,n0_qdp),qedges(:,:,k,q),&
                deriv, elem(ie), u_is_contra=.false.)

           ! advance in time. GLL quadrature, cardinal function basis, under-integrated.
           ! local mass matrix is diagonal, with entries elem(ie)%spheremp(),
           ! so we divide through by elem(ie)%spheremp().
           elem(ie)%state%Qdp(:,:,k,q,np1_qdp)=elem(ie)%state%Qdp(:,:,k,q,n0_qdp) - dt*divdp/elem(ie)%spheremp

           if (npdg<np) then
              ! modal timestep, with exact integration.  using prognostic variable: p*metdet
              ! local mass matrix is diagonal assuming npdg<np so that GLL quadrature is exact)
              ! (note: GLL/modal conversion comutes with time-stepping)

              ! compute modal coefficients of p*metdet
              ! (spherical inner-product of Legendre polynomial and p)
              pshat = gll_to_dgmodal(elem(ie)%state%Qdp(:,:,k,q,np1_qdp)*elem(ie)%metdet(:,:),deriv)

              ! modal based limiter goes here
              ! apply a little dissipation to last mode:
              do j=1,npdg
              do i=1,npdg
                 !if ( (i-1)+(j-1) == 4) pshat(i,j)=pshat(i,j)*.75
                 !if ( (i-1)+(j-1) == 3) pshat(i,j)=pshat(i,j)*.90
                 if ( i==npdg) pshat(i,j)=pshat(i,j)*.90
                 if ( j==npdg) pshat(i,j)=pshat(i,j)*.90
              enddo
              enddo


              ! evalute modal expanion of p*metdet on GLL points
              divdp=dgmodal_to_gll(pshat,deriv)

              ! convert from p*metdet back to p:
              elem(ie)%state%Qdp(:,:,k,q,np1_qdp)=divdp/elem(ie)%metdet(:,:)
           endif
        enddo
        if(limiter_option == 4)then
           ! reuse CG limiter, which wants Qdp*spheremp:
           do k=1,nlev
              elem(ie)%state%Qdp(:,:,k,q,np1_qdp)=elem(ie)%state%Qdp(:,:,k,q,np1_qdp)*elem(ie)%spheremp(:,:)
           enddo
           call limiter2d_zero(elem(ie)%state%Qdp(:,:,:,q,np1_qdp),hvcoord)
           do k=1,nlev
              elem(ie)%state%Qdp(:,:,k,q,np1_qdp)=elem(ie)%state%Qdp(:,:,k,q,np1_qdp)/elem(ie)%spheremp(:,:)
           enddo
        endif
     end do
  end do
#ifdef DEBUGOMP
#if (! defined ELEMENT_OPENMP)
!$OMP BARRIER
#endif
#endif
  call t_stopf('euler_step_dg')

  end subroutine euler_step_dg

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

  subroutine limiter_optim_iter_full(ptens,sphweights,minp,maxp,dpmass)
    !THIS IS A NEW VERSION OF LIM8, POTENTIALLY FASTER BECAUSE INCORPORATES KNOWLEDGE FROM
    !PREVIOUS ITERATIONS

    !The idea here is the following: We need to find a grid field which is closest
    !to the initial field (in terms of weighted sum), but satisfies the min/max constraints.
    !So, first we find values which do not satisfy constraints and bring these values
    !to a closest constraint. This way we introduce some mass change (addmass),
    !so, we redistribute addmass in the way that l2 error is smallest.
    !This redistribution might violate constraints thus, we do a few iterations.
    real (kind=8), dimension(4*4), intent(inout)            :: ptens
    real (kind=8), dimension(4*4), intent(in   )            :: sphweights
    real (kind=8),  intent(inout)            :: minp
    real (kind=8),  intent(inout)            :: maxp
    real (kind=8), dimension(4*4), intent(in   ), optional  :: dpmass

    !local
    real (kind=8), dimension(4*4) :: weights
    integer  k1, i, j, iter, i1, i2
    integer :: whois_neg(4*4), whois_pos(4*4), neg_counter, pos_counter
    real (kind=8) :: addmass, weightssum, mass
    real (kind=8) :: x(4*4),c(4*4)
    real (kind=8) :: al_neg(4*4), al_pos(4*4), howmuch
    real (kind=8) :: tol_limiter = 1e-15
    integer, parameter :: maxiter = 5

      weights(:) = sphweights(:) * dpmass(:)
      ptens(:) = ptens(:) / dpmass(:)

      c = weights(:)
      x = ptens(:)

      mass = sum(c*x)

      ! relax constraints to ensure limiter has a solution:
      ! This is only needed if runnign with the SSP CFL>1 or
      ! due to roundoff errors
      if( (mass / sum(c)) < minp ) then
        minp = mass / sum(c)
      endif
      if( (mass / sum(c)) > maxp ) then
        maxp = mass / sum(c)
      endif

      addmass = 0.0d0
      pos_counter = 0;
      neg_counter = 0;

      ! apply constraints, compute change in mass caused by constraints
      do k1 = 1 , 4*4
        if ( ( x(k1) >= maxp ) ) then
          addmass = addmass + ( x(k1) - maxp ) * c(k1)
          x(k1) = maxp
          whois_pos(k1) = -1
        else
          pos_counter = pos_counter+1;
          whois_pos(pos_counter) = k1;
        endif
        if ( ( x(k1) <= minp ) ) then
          addmass = addmass - ( minp - x(k1) ) * c(k1)
          x(k1) = minp
          whois_neg(k1) = -1
        else
          neg_counter = neg_counter+1;
          whois_neg(neg_counter) = k1;
        endif
      enddo

      ! iterate to find field that satifies constraints and is l2-norm closest to original
      weightssum = 0.0d0
      if ( addmass > 0 ) then
        do i2 = 1 , maxIter
          weightssum = 0.0
          do k1 = 1 , pos_counter
            i1 = whois_pos(k1)
            weightssum = weightssum + c(i1)
            al_pos(i1) = maxp - x(i1)
          enddo

          if( ( pos_counter > 0 ) .and. ( addmass > tol_limiter * abs(mass) ) ) then
            do k1 = 1 , pos_counter
              i1 = whois_pos(k1)
              howmuch = addmass / weightssum
              if ( howmuch > al_pos(i1) ) then
                howmuch = al_pos(i1)
                whois_pos(k1) = -1
              endif
              addmass = addmass - howmuch * c(i1)
              weightssum = weightssum - c(i1)
              x(i1) = x(i1) + howmuch
            enddo
            !now sort whois_pos and get a new number for pos_counter
            !here neg_counter and whois_neg serve as temp vars
            neg_counter = pos_counter
            whois_neg = whois_pos
            whois_pos = -1
            pos_counter = 0
            do k1 = 1 , neg_counter
              if ( whois_neg(k1) .ne. -1 ) then
                pos_counter = pos_counter+1
                whois_pos(pos_counter) = whois_neg(k1)
              endif
            enddo
          else
            exit
          endif
        enddo
      else
         do i2 = 1 , maxIter
           weightssum = 0.0
           do k1 = 1 , neg_counter
             i1 = whois_neg(k1)
             weightssum = weightssum + c(i1)
             al_neg(i1) = x(i1) - minp
           enddo

           if ( ( neg_counter > 0 ) .and. ( (-addmass) > tol_limiter * abs(mass) ) ) then
             do k1 = 1 , neg_counter
               i1 = whois_neg(k1)
               howmuch = -addmass / weightssum
               if ( howmuch > al_neg(i1) ) then
                 howmuch = al_neg(i1)
                 whois_neg(k1) = -1
               endif
               addmass = addmass + howmuch * c(i1)
               weightssum = weightssum - c(i1)
               x(i1) = x(i1) - howmuch
             enddo
             !now sort whois_pos and get a new number for pos_counter
             !here pos_counter and whois_pos serve as temp vars
             pos_counter = neg_counter
             whois_pos = whois_neg
             whois_neg = -1
             neg_counter = 0
             do k1 = 1 , pos_counter
               if ( whois_pos(k1) .ne. -1 ) then
                 neg_counter = neg_counter+1
                 whois_neg(neg_counter) = whois_pos(k1)
               endif
             enddo
           else
             exit
           endif
         enddo
      endif

      ptens(:) = x

      ptens(:) = ptens(:) * dpmass(:)
  end subroutine limiter_optim_iter_full


!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

  subroutine limiter2d_minmax(Q,dp,hvcoord,spheremp,qmin,qmax)
!
! mass conserving limiter (2D only).  to be called just before DSS
!
! in pure 2D advection, the element mass will not be negative before DSS
! this routine will redistribute to remove negative values (conservative)
!
! if used in 3D, should be applied with 2D/vertical split advection
!
! call with Qdp and assocated dp
!
!
  implicit none
  real (kind=real_kind), intent(inout) :: Q(np,np,nlev)
  real (kind=real_kind), intent(in   ) :: spheremp(np,np)
  real (kind=real_kind), intent(in   ) ::  dp(np,np,nlev)
  type (hvcoord_t)     , intent(in   ) :: hvcoord

  ! local
  integer i,j,k
  real (kind=real_kind) :: mass,mass_new,area,qmin(nlev),qmax(nlev),mass2

#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k,mass,area,mass2,mass_new,i,j)
#endif
  do k = 1 , nlev
    mass = sum( Q(:,:,k)*spheremp(:,:) )
    area = sum( dp(:,:,k)*spheremp(:,:) )

    Q(:,:,k) = Q(:,:,k) / dp(:,:,k)  ! % convert to concentration

!   if (mass>0) print *,k,mass/area,qmin(k),qmax(k)
    ! max limiter
    if ( maxval(Q(:,:,k)) > qmax(k) ) then
      Q(:,:,k) = qmax(k) - Q(:,:,k)      ! some of these will be negative
      mass2 = area * qmax(k) - mass

      if (mass2 < 0) Q(:,:,k) = -Q(:,:,k)
      mass_new = 0
      do j = 1 , np
        do i = 1 , np
          if ( Q(i,j,k) < 0 ) then
            Q(i,j,k) = 0
          else
            mass_new = mass_new + Q(i,j,k) * dp(i,j,k) * spheremp(i,j)
          endif
        enddo
      enddo

      ! now scale the all positive values to restore mass
      if ( mass_new > 0 ) Q(:,:,k) = Q(:,:,k) * abs(mass2) / mass_new
      if ( mass2    < 0 ) Q(:,:,k) = -Q(:,:,k)

      Q(:,:,k) = qmax(k) - Q(:,:,k)
    endif

    ! min limiter
    if ( minval(Q(:,:,k)) < qmin(k) ) then
      Q(:,:,k) = Q(:,:,k) - qmin(k)
      mass2 = mass - area * qmin(k)
      ! negative mass.  so reduce all postive values to zero
      ! then increase negative values as much as possible
      if ( mass2 < 0 ) Q(:,:,k) = -Q(:,:,k)
      mass_new = 0
      do j = 1 , np
        do i = 1 , np
          if ( Q(i,j,k) < 0 ) then
            Q(i,j,k) = 0
          else
            mass_new = mass_new + Q(i,j,k) * dp(i,j,k) * spheremp(i,j)
          endif
        enddo
      enddo

      ! now scale the all positive values to restore mass
      if ( mass_new > 0 ) Q(:,:,k) = Q(:,:,k) * abs(mass2) / mass_new
      if ( mass2    < 0 ) Q(:,:,k) = -Q(:,:,k)

      Q(:,:,k) = Q(:,:,k) + qmin(k)
    endif

    Q(:,:,k) = Q(:,:,k) * dp(:,:,k)
  enddo
  end subroutine limiter2d_minmax

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

  subroutine limiter2d_zero(Q,hvcoord)
  ! mass conserving zero limiter (2D only).  to be called just before DSS
  !
  ! this routine is called inside a DSS loop, and so Q had already
  ! been multiplied by the mass matrix.  Thus dont include the mass
  ! matrix when computing the mass = integral of Q over the element
  !
  ! ps is only used when advecting Q instead of Qdp
  ! so ps should be at one timelevel behind Q
  implicit none
  real (kind=real_kind), intent(inout) :: Q(np,np,nlev)
  type (hvcoord_t)     , intent(in   ) :: hvcoord

  ! local
  real (kind=real_kind) :: dp(np,np)
  real (kind=real_kind) :: mass,mass_new,ml
  integer i,j,k

  do k = nlev , 1 , -1
    mass = 0
    do j = 1 , np
      do i = 1 , np
        !ml = Q(i,j,k)*dp(i,j)*spheremp(i,j)  ! see above
        ml = Q(i,j,k)
        mass = mass + ml
      enddo
    enddo

    ! negative mass.  so reduce all postive values to zero
    ! then increase negative values as much as possible
    if ( mass < 0 ) Q(:,:,k) = -Q(:,:,k)
    mass_new = 0
    do j = 1 , np
      do i = 1 , np
        if ( Q(i,j,k) < 0 ) then
          Q(i,j,k) = 0
        else
          ml = Q(i,j,k)
          mass_new = mass_new + ml
        endif
      enddo
    enddo

    ! now scale the all positive values to restore mass
    if ( mass_new > 0 ) Q(:,:,k) = Q(:,:,k) * abs(mass) / mass_new
    if ( mass     < 0 ) Q(:,:,k) = -Q(:,:,k)
  enddo
  end subroutine limiter2d_zero

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

  subroutine advance_hypervis_scalar( edgeAdv , elem , hvcoord , hybrid , deriv , nt , nt_qdp , nets , nete , dt2 )
  !  hyperviscsoity operator for foward-in-time scheme
  !  take one timestep of:
  !          Q(:,:,:,np) = Q(:,:,:,np) +  dt2*nu*laplacian**order ( Q )
  !
  !  For correct scaling, dt2 should be the same 'dt2' used in the leapfrog advace
#ifdef _ACCEL
  use cuda_mod       , only : advance_hypervis_scalar_cuda
#endif
  use kinds          , only : real_kind
  use dimensions_mod , only : np, nlev
  use hybrid_mod     , only : hybrid_t
  use element_mod    , only : element_t
  use derivative_mod , only : derivative_t
  use edge_mod       , only : EdgeBuffer_t, edgevpack, edgevunpack
  use bndry_mod      , only : bndry_exchangev
  use perf_mod       , only : t_startf, t_stopf                          ! _EXTERNAL
  implicit none
  type (EdgeBuffer_t)  , intent(inout)         :: edgeAdv
  type (element_t)     , intent(inout), target :: elem(:)
  type (hvcoord_t)     , intent(in   )         :: hvcoord
  type (hybrid_t)      , intent(in   )         :: hybrid
  type (derivative_t)  , intent(in   )         :: deriv
  integer              , intent(in   )         :: nt
  integer              , intent(in   )         :: nt_qdp
  integer              , intent(in   )         :: nets
  integer              , intent(in   )         :: nete
  real (kind=real_kind), intent(in   )         :: dt2

  ! local
  real (kind=real_kind), dimension(np,np,nlev,qsize,nets:nete) :: Qtens
  real (kind=real_kind), dimension(np,np,nlev                ) :: dp
  real (kind=real_kind), dimension(      nlev,qsize,nets:nete) :: min_neigh
  real (kind=real_kind), dimension(      nlev,qsize,nets:nete) :: max_neigh
  integer :: k , kptr , i , j , ie , ic , q

! NOTE: PGI compiler bug: when using spheremp, rspheremp and ps as pointers to elem(ie)% members,
!       data is incorrect (offset by a few numbers actually)
!       removed for now.
!  real (kind=real_kind), dimension(:,:), pointer :: spheremp,rspheremp
  real (kind=real_kind), dimension(np,np) :: lap_p
  real (kind=real_kind) :: v1,v2,dt,dp0
  integer :: density_scaling = 0
  if ( nu_q           == 0 ) return
  if ( hypervis_order /= 2 ) return
#ifdef _ACCEL
  call advance_hypervis_scalar_cuda( edgeAdv , elem , hvcoord , hybrid , deriv , nt , nt_qdp , nets , nete , dt2 )
  return
#endif
  call t_barrierf('sync_advance_hypervis_scalar', hybrid%par%comm)
  call t_startf('advance_hypervis_scalar')

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  hyper viscosity
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  dt = dt2 / hypervis_subcycle_q

  do ic = 1 , hypervis_subcycle_q
    do ie = nets , nete
      ! Qtens = Q/dp   (apply hyperviscsoity to dp0 * Q, not Qdp)
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k,q)
#endif
      do k = 1 , nlev
         ! apply dissipation to Q, not Qdp, for tracer/mass consistency
        dp(:,:,k) = elem(ie)%derived%dp(:,:,k) - dt2*elem(ie)%derived%divdp_proj(:,:,k)
        do q = 1 , qsize
          Qtens(:,:,k,q,ie) = elem(ie)%state%Qdp(:,:,k,q,nt_qdp) / dp(:,:,k)
        enddo
      enddo
    enddo

    ! compute biharmonic operator. Qtens = input and output
    call biharmonic_wk_scalar( elem , Qtens , deriv , edgeAdv , hybrid , nets , nete )
    do ie = nets , nete
      !spheremp     => elem(ie)%spheremp
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(q,k,i,j)
#endif
      do q = 1 , qsize
        do k = 1 , nlev
          dp0 = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) ) * hvcoord%ps0 + &
                ( hvcoord%hybi(k+1) - hvcoord%hybi(k) ) * hvcoord%ps0
          do j = 1 , np
            do i = 1 , np
              ! advection Qdp.  For mass advection consistency:
              ! DIFF( Qdp) ~   dp0 DIFF (Q)  =  dp0 DIFF ( Qdp/dp )
              elem(ie)%state%Qdp(i,j,k,q,nt_qdp) = elem(ie)%state%Qdp(i,j,k,q,nt_qdp) * elem(ie)%spheremp(i,j) &
                                                   - dt * nu_q * dp0 * Qtens(i,j,k,q,ie)
            enddo
          enddo
        enddo

        ! smooth some of the negativities introduced by diffusion:
        call limiter2d_zero( elem(ie)%state%Qdp(:,:,:,q,nt_qdp) , hvcoord )
      enddo
      call edgeVpack  ( edgeAdv , elem(ie)%state%Qdp(:,:,:,:,nt_qdp) , qsize*nlev , 0 , elem(ie)%desc )
    enddo

    call bndry_exchangeV( hybrid , edgeAdv )

    do ie = nets , nete
      call edgeVunpack( edgeAdv , elem(ie)%state%Qdp(:,:,:,:,nt_qdp) , qsize*nlev , 0 , elem(ie)%desc )
      !rspheremp     => elem(ie)%rspheremp
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(q,k)
#endif
      do q = 1 , qsize
        ! apply inverse mass matrix
        do k = 1 , nlev
          elem(ie)%state%Qdp(:,:,k,q,nt_qdp) = elem(ie)%rspheremp(:,:) * elem(ie)%state%Qdp(:,:,k,q,nt_qdp)
        enddo
      enddo
    enddo
#ifdef DEBUGOMP
#if (! defined ELEMENT_OPENMP)
!$OMP BARRIER
#endif
#endif
  enddo
  call t_stopf('advance_hypervis_scalar')
  end subroutine advance_hypervis_scalar





  subroutine vertical_remap(elem,fvm,hvcoord,dt,np1,np1_qdp,nets,nete)
  ! This routine is called at the end of the vertically Lagrangian
  ! dynamics step to compute the vertical flux needed to get back
  ! to reference eta levels
  !
  ! input:
  !     derived%dp()  delta p on levels at beginning of timestep
  !     state%dp3d(np1)  delta p on levels at end of timestep
  ! output:
  !     state%ps_v(np1)          surface pressure at time np1
  !     derived%eta_dot_dpdn()   vertical flux from final Lagrangian
  !                              levels to reference eta levels
  !
  use kinds, only : real_kind
  use hybvcoord_mod, only : hvcoord_t
  use vertremap_mod, only : my_vertical_remap_acc, remap1, remap1_nofilter, remap_q_ppm ! _EXTERNAL (actually INTERNAL)
  use control_mod, only :  rsplit
  use parallel_mod, only : abortmp
  use fvm_control_volume_mod, only : fvm_struct

  type(fvm_struct), intent(inout) :: fvm(:)

  !    type (hybrid_t), intent(in)       :: hybrid  ! distributed parallel structure (shared)
  type (element_t), intent(inout)   :: elem(:)
  type (hvcoord_t)                  :: hvcoord
  real (kind=real_kind)             :: dt,sga

  integer :: ie,i,j,k,np1,nets,nete,np1_qdp
  real (kind=real_kind), dimension(np,np,nlev)  :: dp,dp_star
  real (kind=real_kind), dimension(np,np,nlev,2)  :: ttmp

  real(kind=real_kind) :: ps0
  integer(kind=8), dimension(7, nets:nete) :: elem_array

  call t_startf('vertical_remap')
!  do ie = nets, nete
!     elem_array(1,ie) = loc(elem(ie)%state%ps_v(:,:,np1))
!     elem_array(2,ie) = loc(elem(ie)%state%dp3d(:,:,:,np1))
!     elem_array(3,ie) = loc(elem(ie)%state%t(:,:,:,np1))
!     elem_array(4,ie) = loc(elem(ie)%state%v(:,:,:,:,np1))
!     elem_array(5,ie) = loc(elem(ie)%state%Qdp(:,:,:,:,np1))
!     elem_array(6,ie) = loc(hvcoord%hyai(:)) 
!     elem_array(7,ie) = loc(hvcoord%hybi(:))
!  enddo
!  ps0 = hvcoord%ps0
!  call my_vertical_remap_acc(elem, hvcoord, ps0, nets, nete, nlev, qsize, np)
      

  do ie=nets,nete
        elem(ie)%state%ps_v(:,:,np1) = hvcoord%hyai(1)*hvcoord%ps0 + &
             sum(elem(ie)%state%dp3d(:,:,:,np1),3)
        do k=1,nlev
           dp(:,:,k) = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
                ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(:,:,np1)
           dp_star(:,:,k) = elem(ie)%state%dp3d(:,:,k,np1)
        enddo

        ttmp(:,:,:,1)=elem(ie)%state%t(:,:,:,np1)
        ttmp(:,:,:,1)=ttmp(:,:,:,1)*dp_star
        call remap1(ttmp,np,1,dp_star,dp)
        elem(ie)%state%t(:,:,:,np1)=ttmp(:,:,:,1)/dp

        ttmp(:,:,:,1)=elem(ie)%state%v(:,:,1,:,np1)*dp_star
        ttmp(:,:,:,2)=elem(ie)%state%v(:,:,2,:,np1)*dp_star
        call remap1(ttmp,np,2,dp_star,dp)
        elem(ie)%state%v(:,:,1,:,np1)=ttmp(:,:,:,1)/dp
        elem(ie)%state%v(:,:,2,:,np1)=ttmp(:,:,:,2)/dp

        call remap1(elem(ie)%state%Qdp(:,:,:,:,np1_qdp),np,qsize,dp_star,dp)



  enddo
  
  call t_stopf('vertical_remap')
  end subroutine vertical_remap




end module prim_advection_mod











! TPP Processed 06af35b4b6e1d5ae0884d64cc00c796d

! TPP Processed 06af35b4b6e1d5ae0884d64cc00c796d

! TPP Processed 06af35b4b6e1d5ae0884d64cc00c796d

! TPP Processed 06af35b4b6e1d5ae0884d64cc00c796d
