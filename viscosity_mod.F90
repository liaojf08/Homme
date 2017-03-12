#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#define constLev 128

module viscosity_mod
!
!  This module should be renamed "global_deriv_mod.F90"
!
!  It is a collection of derivative operators that must be applied to the field
!  over the sphere (as opposed to derivative operators that can be applied element
!  by element)
!
!
use kinds, only : real_kind, iulog
use dimensions_mod, only : np, nlev,qsize, max_corner_elem
use hybrid_mod, only : hybrid_t
use element_mod, only : element_t
use derivative_mod, only : derivative_t, laplace_sphere_wk, vlaplace_sphere_wk, vorticity_sphere, derivinit, divergence_sphere
use edge_mod, only : EdgeBuffer_t, edgevpack, edgerotate, edgevunpack, edgevunpackmin, &
    edgevunpackmax, initEdgeBuffer, FreeEdgeBuffer, EdgeDescriptor_t
use bndry_mod, only : bndry_exchangev
use control_mod, only : hypervis_scaling, north, south, east, west, neast, nwest, seast, swest, which_vlaplace

implicit none
save

public :: biharmonic_wk
#ifdef _PRIM
public :: biharmonic_wk_scalar
public :: biharmonic_wk_scalar_minmax
#endif
public :: compute_zeta_C0
public :: compute_div_C0
public :: compute_zeta_C0_2d
public :: compute_div_C0_2d
public :: test_ibyp

interface compute_zeta_C0_2d
    module procedure compute_zeta_C0_2d_sphere
    module procedure compute_zeta_C0_2d_contra
end interface

interface compute_div_C0_2d
    module procedure compute_div_C0_2d_sphere
    module procedure compute_div_C0_2d_contra
end interface


type (EdgeBuffer_t)          :: edge1

contains

#ifdef _PRIM
subroutine biharmonic_wk(elem,pstens,ptens,vtens,deriv,edge3,hybrid,nt,nets,nete,nu_ratio)
#else
subroutine biharmonic_wk(elem,ptens,vtens,deriv,edge3,hybrid,nt,nets,nete,nu_ratio)
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! compute weak biharmonic operator
!    input:  h,v (stored in elem()%, in lat-lon coordinates
!    output: ptens,vtens  overwritten with weak biharmonic of h,v (output in lat-lon coordinates)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
type (hybrid_t)      , intent(in) :: hybrid
type (element_t)     , intent(inout), target :: elem(:)
integer :: nt,nets,nete
real (kind=real_kind), dimension(np,np,2,nlev,nets:nete)  :: vtens
real (kind=real_kind), dimension(np,np,nlev,nets:nete) :: ptens
type (EdgeBuffer_t)  , intent(inout) :: edge3
type (derivative_t)  , intent(in) :: deriv
real (kind=real_kind) ::  nu_ratio
#ifdef _PRIM
real (kind=real_kind), dimension(np,np,nets:nete) :: pstens
#endif

! local
integer :: k,kptr,i,j,ie,ic
real (kind=real_kind), dimension(:,:), pointer :: rspheremv
real (kind=real_kind), dimension(np,np) :: lap_ps
real (kind=real_kind), dimension(np,np,nlev) :: T
real (kind=real_kind), dimension(np,np,2) :: v

   do ie=nets,nete

#ifdef _PRIM
      ! should filter lnps + PHI_s/RT?

      !if tensor hyperviscosity with tensor V is used, then biharmonic operator is (\grad\cdot V\grad) (\grad \cdot \grad)
      if(hypervis_scaling > 0)then
	pstens(:,:,ie)=laplace_sphere_wk(elem(ie)%state%ps_v(:,:,nt),deriv,elem(ie),var_coef=.false.)
      else
	pstens(:,:,ie)=laplace_sphere_wk(elem(ie)%state%ps_v(:,:,nt),deriv,elem(ie),var_coef=.true.)
      endif
#endif

#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k, j, i)
#endif
      do k=1,nlev
         do j=1,np
            do i=1,np
#ifdef _PRIM
               T(i,j,k)=elem(ie)%state%T(i,j,k,nt)
#elif defined _PRIMDG
            T(i,j,k)=elem(ie)%state%p(i,j,k,nt) + elem(ie)%state%phis(i,j)
#else
               ! filter surface height, not thickness
               T(i,j,k)=elem(ie)%state%p(i,j,k,nt) + elem(ie)%state%ps(i,j)
#endif
            enddo
         enddo

         if(hypervis_scaling > 0)then
           ptens(:,:,k,ie)=laplace_sphere_wk(T(:,:,k),deriv,elem(ie),var_coef=.false.)
           vtens(:,:,:,k,ie)=vlaplace_sphere_wk(elem(ie)%state%v(:,:,:,k,nt),deriv,&
                elem(ie),var_coef=.false.,nu_ratio=nu_ratio)
         else
           ptens(:,:,k,ie)=laplace_sphere_wk(T(:,:,k),deriv,elem(ie),var_coef=.true.)
           vtens(:,:,:,k,ie)=vlaplace_sphere_wk(elem(ie)%state%v(:,:,:,k,nt),deriv,&
                elem(ie),var_coef=.true.,nu_ratio=nu_ratio)
         endif

      enddo
      kptr=0
      call edgeVpack(edge3, ptens(1,1,1,ie),nlev,kptr,elem(ie)%desc)
      kptr=nlev
      call edgeVpack(edge3, vtens(1,1,1,1,ie),2*nlev,kptr,elem(ie)%desc)

#ifdef _PRIM
      kptr=3*nlev
      call edgeVpack(edge3, pstens(1,1,ie),1,kptr,elem(ie)%desc)
#endif
   enddo

   call bndry_exchangeV(hybrid,edge3)

   do ie=nets,nete
      rspheremv     => elem(ie)%rspheremp(:,:)

      kptr=0
      call edgeVunpack(edge3, ptens(1,1,1,ie), nlev, kptr, elem(ie)%desc)
      kptr=nlev
      call edgeVunpack(edge3, vtens(1,1,1,1,ie), 2*nlev, kptr, elem(ie)%desc)

      ! apply inverse mass matrix, then apply laplace again
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k, j, i, v)
#endif
      do k=1,nlev
         do j=1,np
            do i=1,np
               T(i,j,k)=rspheremv(i,j)*ptens(i,j,k,ie)
               v(i,j,1)=rspheremv(i,j)*vtens(i,j,1,k,ie)
               v(i,j,2)=rspheremv(i,j)*vtens(i,j,2,k,ie)
            enddo
         enddo
         ptens(:,:,k,ie)=laplace_sphere_wk(T(:,:,k),deriv,elem(ie),var_coef=.true.)
         vtens(:,:,:,k,ie)=vlaplace_sphere_wk(v(:,:,:),deriv,elem(ie),var_coef=.true.,&
              nu_ratio=nu_ratio)
      enddo

#ifdef _PRIM
      kptr=3*nlev
      call edgeVunpack(edge3, pstens(1,1,ie), 1, kptr, elem(ie)%desc)
      ! apply inverse mass matrix, then apply laplace again
      lap_ps(:,:)=rspheremv(:,:)*pstens(:,:,ie)
      pstens(:,:,ie)=laplace_sphere_wk(lap_ps,deriv,elem(ie),var_coef=.true.)
#endif

   enddo
#ifdef DEBUGOMP
#if (! defined ELEMENT_OPENMP)
!$OMP BARRIER
#endif
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end subroutine

subroutine my_edgeVunpack_viscosity_test(edge_nlyr, edge_nbuf, edge_buf, is, ie, in, iw, &
    desc_getmapP, v, my_vlyr,kptr, my_swest, my_max_corner_elem, &
    swest_buf, seast_buf, neast_buf, nwest_buf)

  integer, intent(in) :: edge_nlyr, edge_nbuf, is,ie,in,iw, my_swest, my_max_corner_elem
  real(kind=8), dimension(edge_nlyr, edge_nbuf), intent(in) :: edge_buf
  real(kind=8), dimension(constLev), intent(in) :: swest_buf
  real(kind=8), dimension(constLev), intent(in) :: seast_buf
  real(kind=8), dimension(constLev), intent(in) :: neast_buf
  real(kind=8), dimension(constLev), intent(in) :: nwest_buf
  integer, dimension(:), intent(in) :: desc_getmapP

  integer,               intent(in)  :: my_vlyr
  real (kind=8), intent(inout) :: v(4,4,constLev)
  integer,               intent(in)  :: kptr

  ! Local
  logical, parameter :: UseUnroll = .TRUE.
  integer :: i,k,ll

  if(MODULO(np,2) == 0 .and. UseUnroll) then
     do k=1,my_vlyr
        do i=1,np,2
           v(i  ,1  ,k) = v(i  ,1  ,k)+edge_buf(kptr+k,is+i  )
           v(i+1,1  ,k) = v(i+1,1  ,k)+edge_buf(kptr+k,is+i+1)
           v(np ,i  ,k) = v(np ,i  ,k)+edge_buf(kptr+k,ie+i  )
           v(np ,i+1,k) = v(np ,i+1,k)+edge_buf(kptr+k,ie+i+1)
           v(i  ,np ,k) = v(i  ,np ,k)+edge_buf(kptr+k,in+i  )
           v(i+1,np ,k) = v(i+1,np ,k)+edge_buf(kptr+k,in+i+1)
           v(1  ,i  ,k) = v(1  ,i  ,k)+edge_buf(kptr+k,iw+i  )
           v(1  ,i+1,k) = v(1  ,i+1,k)+edge_buf(kptr+k,iw+i+1)
        end do
     end do
  else
     do k=1,my_vlyr
        do i=1,np
           v(i  ,1  ,k) = v(i  ,1  ,k)+edge_buf(kptr+k,is+i  )
           v(np ,i  ,k) = v(np ,i  ,k)+edge_buf(kptr+k,ie+i  )
           v(i  ,np ,k) = v(i  ,np ,k)+edge_buf(kptr+k,in+i  )
           v(1  ,i  ,k) = v(1  ,i  ,k)+edge_buf(kptr+k,iw+i  )
        end do
     end do
  endif

! swest
  do ll=my_swest,my_swest+my_max_corner_elem-1
    !write(*,*), 'swest loop, ll = ', ll, 'desc_getmapP(ll) = ', desc_getmapP(ll) ! 5
      if(desc_getmapP(ll) /= -1) then
          do k=1,my_vlyr
              !v(1  ,1 ,k)=v(1 ,1 ,k)+edge_buf(kptr+k,desc_getmapP(ll)+1)
              v(1  ,1 ,k)=v(1 ,1 ,k)+swest_buf(kptr + k)
          enddo
      endif
  end do

! SEAST
  do ll=my_swest+my_max_corner_elem,my_swest+2*my_max_corner_elem-1
    !write(*,*), 'seast loop, ll = ', ll, 'desc_getmapP(ll) = ', desc_getmapP(ll) !6
      if(desc_getmapP(ll) /= -1) then
          do k=1,my_vlyr
              !v(np ,1 ,k)=v(np,1 ,k)+edge_buf(kptr+k,desc_getmapP(ll)+1)
              v(np ,1 ,k)=v(np,1 ,k)+seast_buf(kptr + k)
          enddo
      endif
  end do

! NEAST
  do ll=my_swest+3*my_max_corner_elem,my_swest+4*my_max_corner_elem-1
    !write(*,*), 'neast loop, ll = ', ll, 'desc_getmapP(ll) = ', desc_getmapP(ll) !8
      if(desc_getmapP(ll) /= -1) then
          do k=1,my_vlyr
              !v(np ,np,k)=v(np,np,k)+edge_buf(kptr+k,desc_getmapP(ll)+1)
              v(np ,np,k)=v(np,np,k)+ neast_buf(kptr + k)
          enddo
      endif
  end do

! NWEST
  do ll=my_swest+2*my_max_corner_elem,my_swest+3*my_max_corner_elem-1
    !write(*,*), 'nwest loop, ll = ', ll, 'desc_getmapP(ll) = ', desc_getmapP(ll) !7
      if(desc_getmapP(ll) /= -1) then
          do k=1,my_vlyr
              !v(1  ,np,k)=v(1 ,np,k)+edge_buf(kptr+k,desc_getmapP(ll)+1)
              v(1  ,np,k)=v(1 ,np,k)+nwest_buf(kptr+k)
          enddo
      endif
  end do

end subroutine

subroutine my_edgeVunpack_viscosity(edge_nlyr, edge_nbuf, &
    desc_getmapP, v, my_vlyr,kptr, my_swest, my_max_corner_elem, &
    swest_buf, seast_buf, neast_buf, nwest_buf, &
    edge_buf_is, edge_buf_iw, edge_buf_ie, edge_buf_in)

  integer, intent(in) :: edge_nlyr, edge_nbuf, my_swest, my_max_corner_elem
  real(kind=8), dimension(edge_nlyr), intent(in) :: swest_buf
  real(kind=8), dimension(edge_nlyr), intent(in) :: seast_buf
  real(kind=8), dimension(edge_nlyr), intent(in) :: neast_buf
  real(kind=8), dimension(edge_nlyr), intent(in) :: nwest_buf
  real(kind=8), dimension(edge_nlyr, 4), intent(in) :: edge_buf_is
  real(kind=8), dimension(edge_nlyr, 4), intent(in) :: edge_buf_iw
  real(kind=8), dimension(edge_nlyr, 4), intent(in) :: edge_buf_ie
  real(kind=8), dimension(edge_nlyr, 4), intent(in) :: edge_buf_in
  integer, dimension(:), intent(in) :: desc_getmapP

  integer,               intent(in)  :: my_vlyr
  real (kind=8), intent(inout) :: v(4,4,constLev)
  integer,               intent(in)  :: kptr

  ! Local
  logical, parameter :: UseUnroll = .TRUE.
  integer :: i,k,ll

  if(MODULO(np,2) == 0 .and. UseUnroll) then
     do k=1,my_vlyr
        do i=1,np,2
           v(i  ,1  ,k) = v(i  ,1  ,k)+edge_buf_is(kptr+k,i  )
           v(i+1,1  ,k) = v(i+1,1  ,k)+edge_buf_is(kptr+k,i+1)
           v(1  ,i  ,k) = v(1  ,i  ,k)+edge_buf_iw(kptr+k,i  )
           v(1  ,i+1,k) = v(1  ,i+1,k)+edge_buf_iw(kptr+k,i+1)
           v(np ,i  ,k) = v(np ,i  ,k)+edge_buf_ie(kptr+k,i  )
           v(np ,i+1,k) = v(np ,i+1,k)+edge_buf_ie(kptr+k,i+1)
           v(i  ,np ,k) = v(i  ,np ,k)+edge_buf_in(kptr+k,i  )
           v(i+1,np ,k) = v(i+1,np ,k)+edge_buf_in(kptr+k,i+1)
        end do
     end do
  else
     do k=1,my_vlyr
        do i=1,np
           v(i  ,1  ,k) = v(i  ,1  ,k)+edge_buf_is(kptr+k,i  )
           v(1  ,i  ,k) = v(1  ,i  ,k)+edge_buf_iw(kptr+k,i  )
           v(np ,i  ,k) = v(np ,i  ,k)+edge_buf_ie(kptr+k,i  )
           v(i  ,np ,k) = v(i  ,np ,k)+edge_buf_in(kptr+k,i  )
        end do
     end do
  endif

! swest
  do ll=my_swest,my_swest+my_max_corner_elem-1
    !write(*,*), 'swest loop, ll = ', ll, 'desc_getmapP(ll) = ', desc_getmapP(ll) ! 5
      if(desc_getmapP(ll) /= -1) then
          do k=1,my_vlyr
              v(1  ,1 ,k)=v(1 ,1 ,k)+swest_buf(kptr + k)
          enddo
      endif
  end do

! SEAST
  do ll=my_swest+my_max_corner_elem,my_swest+2*my_max_corner_elem-1
    !write(*,*), 'seast loop, ll = ', ll, 'desc_getmapP(ll) = ', desc_getmapP(ll) !6
      if(desc_getmapP(ll) /= -1) then
          do k=1,my_vlyr
              v(np ,1 ,k)=v(np,1 ,k)+seast_buf(kptr + k)
          enddo
      endif
  end do

! NEAST
  do ll=my_swest+3*my_max_corner_elem,my_swest+4*my_max_corner_elem-1
    !write(*,*), 'neast loop, ll = ', ll, 'desc_getmapP(ll) = ', desc_getmapP(ll) !8
      if(desc_getmapP(ll) /= -1) then
          do k=1,my_vlyr
              v(np ,np,k)=v(np,np,k)+ neast_buf(kptr + k)
          enddo
      endif
  end do

! NWEST
  do ll=my_swest+2*my_max_corner_elem,my_swest+3*my_max_corner_elem-1
    !write(*,*), 'nwest loop, ll = ', ll, 'desc_getmapP(ll) = ', desc_getmapP(ll) !7
      if(desc_getmapP(ll) /= -1) then
          do k=1,my_vlyr
              v(1  ,np,k)=v(1 ,np,k)+nwest_buf(kptr+k)
          enddo
      endif
  end do

end subroutine

subroutine my_edgeVunpack_viscosity_all(nets, nete, edge_nlyr, edge_nbuf, &
  edge_buf, my_south, my_east, my_north, my_west, my_elem, ptens, vtens, dptens, my_nlev, &
  my_swest, my_max_corner_elem)

  implicit none
  integer, intent(in) :: nets, nete, edge_nlyr, edge_nbuf, my_nlev, my_south, my_east, my_north, my_west, my_swest, my_max_corner_elem

  real(kind=8), dimension(edge_nlyr, edge_nbuf), intent(in) :: edge_buf
  type (element_t), intent(in) :: my_elem(nets:nete)

  real (kind=8), intent(inout) :: ptens(4,4,constLev,nets:nete)
  real (kind=8), intent(inout) :: vtens(4,4,2,constLev,nets:nete)
  real (kind=8), intent(inout) :: dptens(4,4,constLev,nets:nete)

  integer, dimension(my_swest+4*my_max_corner_elem-1) :: getmapP
  pointer(getmapP_ptr, getmapP)

  integer :: is, ie, iw, in
  pointer(is_ptr, is)
  pointer(ie_ptr, ie)
  pointer(iw_ptr, iw)
  pointer(in_ptr, in)

  real(kind=8), dimension(edge_nlyr) :: swest_buf, seast_buf, neast_buf, nwest_buf
  pointer(swest_buf_ptr, swest_buf)
  pointer(seast_buf_ptr, seast_buf)
  pointer(neast_buf_ptr, neast_buf)
  pointer(nwest_buf_ptr, nwest_buf)

  real(kind=8), dimension(edge_nlyr, 4) :: edge_buf_is, edge_buf_iw, edge_buf_ie, edge_buf_in
  pointer(edge_buf_is_ptr, edge_buf_is)
  pointer(edge_buf_iw_ptr, edge_buf_iw)
  pointer(edge_buf_ie_ptr, edge_buf_ie)
  pointer(edge_buf_in_ptr, edge_buf_in)

  ! local
  integer :: iee, kptr
  integer(kind=8) :: count_start, count_stop, count_rate, count_max

  !call system_clock(count_start)

  !$ACC PARALLEL LOOP copy(ptens,dptens) 
  do iee=nets,nete

    getmapP_ptr = loc(my_elem(iee)%desc%getmapP)
    is_ptr = loc(my_elem(iee)%desc%getmapP(my_south))
    ie_ptr = loc(my_elem(iee)%desc%getmapP(my_east))
    iw_ptr = loc(my_elem(iee)%desc%getmapP(my_west))
    in_ptr = loc(my_elem(iee)%desc%getmapP(my_north))

    swest_buf_ptr = loc(edge_buf(:, getmapP(5) + 1))
    seast_buf_ptr = loc(edge_buf(:, getmapP(6) + 1))
    neast_buf_ptr = loc(edge_buf(:, getmapP(8) + 1))
    nwest_buf_ptr = loc(edge_buf(:, getmapP(7) + 1))

    edge_buf_is_ptr = loc(edge_buf(:,is+1:is+4))
    edge_buf_iw_ptr = loc(edge_buf(:,iw+1:iw+4))
    edge_buf_ie_ptr = loc(edge_buf(:,ie+1:ie+4))
    edge_buf_in_ptr = loc(edge_buf(:,in+1:in+4))

    !!$ACC DATA copyin(getmapP, is, ie, iw, in, swest_buf, seast_buf, neast_buf, nwest_buf, edge_buf_is, edge_buf_in, edge_buf_ie, edge_buf_iw)
    !$ACC DATA copyin(getmapP, is, ie, iw, in, swest_buf, seast_buf, neast_buf, nwest_buf)
    kptr=0
    call my_edgeVunpack_viscosity(edge_nlyr, edge_nbuf, &
      getmapP(:), &
      ptens(:,:,:,iee), my_nlev, kptr, my_swest, my_max_corner_elem, &
      swest_buf, seast_buf, neast_buf, nwest_buf, &
      edge_buf_is, edge_buf_iw, edge_buf_ie, edge_buf_in)

    kptr=3*my_nlev
    call my_edgeVunpack_viscosity(edge_nlyr, edge_nbuf, &
      getmapP(:), &
      dptens(:,:,:,iee), my_nlev, kptr, my_swest, my_max_corner_elem, &
      swest_buf, seast_buf, neast_buf, nwest_buf, &
      edge_buf_is, edge_buf_iw, edge_buf_ie, edge_buf_in)
    
    !$ACC END DATA
  end do
  !$ACC END PARALLEL LOOP
    
  !$ACC PARALLEL LOOP copy(vtens) 
  do iee=nets,nete

    getmapP_ptr = loc(my_elem(iee)%desc%getmapP)
    is_ptr = loc(my_elem(iee)%desc%getmapP(my_south))
    ie_ptr = loc(my_elem(iee)%desc%getmapP(my_east))
    iw_ptr = loc(my_elem(iee)%desc%getmapP(my_west))
    in_ptr = loc(my_elem(iee)%desc%getmapP(my_north))

    swest_buf_ptr = loc(edge_buf(:, getmapP(5) + 1))
    seast_buf_ptr = loc(edge_buf(:, getmapP(6) + 1))
    neast_buf_ptr = loc(edge_buf(:, getmapP(8) + 1))
    nwest_buf_ptr = loc(edge_buf(:, getmapP(7) + 1))

    edge_buf_is_ptr = loc(edge_buf(:,is+1:is+4))
    edge_buf_iw_ptr = loc(edge_buf(:,iw+1:iw+4))
    edge_buf_ie_ptr = loc(edge_buf(:,ie+1:ie+4))
    edge_buf_in_ptr = loc(edge_buf(:,in+1:in+4))
  
    !!$ACC DATA copyin(getmapP, is, ie, iw, in, swest_buf, seast_buf, neast_buf, nwest_buf, edge_buf_is, edge_buf_in, edge_buf_ie, edge_buf_iw)
    !$ACC DATA copyin(getmapP, is, ie, iw, in, swest_buf, seast_buf, neast_buf, nwest_buf)
    kptr=my_nlev
    call my_edgeVunpack_viscosity(edge_nlyr, edge_nbuf, &
      getmapP(:), &
      vtens(:,:,:,:,iee), 2*my_nlev, kptr, my_swest, my_max_corner_elem, &
      swest_buf, seast_buf, neast_buf, nwest_buf, &
      edge_buf_is, edge_buf_iw, edge_buf_ie, edge_buf_in)

    !$ACC END DATA
  end do
  !$ACC END PARALLEL LOOP

  !call system_clock(count_stop, count_rate, count_max)
  !write(*,*) 'normal count = ', (count_stop - count_start)

end subroutine

    subroutine my_gradient_sphere_vis(s,deriv_Dvv,Dinv,ds,my_rrearth)

    real(kind=8), intent(in), dimension(4,4)          :: deriv_Dvv
    real(kind=8), intent(in), dimension(2,2,4,4) :: Dinv
    real(kind=8), intent(in) :: s(4,4)
    real(kind=8), intent(in) :: my_rrearth

    real(kind=8), intent(out) :: ds(4,4,2)

    integer i
    integer j
    integer l

    real(kind=8) ::  dsdx00
    real(kind=8) ::  dsdy00
    real(kind=8) ::  v1(4,4),v2(4,4)

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

  function my_divergence_sphere_wk(v,elem_Dinv, elem_spheremp, deriv_Dvv, my_rrearth) result(div)
    real(kind=8), intent(in) :: v(4,4,2)  ! in lat-lon coordinates
    real(kind=8), dimension(2,2,4,4), intent(in) :: elem_Dinv
    real(kind=8), dimension(4,4), intent(in) :: elem_spheremp
    real(kind=8), dimension(4,4), intent(in) :: deriv_Dvv
    real(kind=8), intent(in) :: my_rrearth

    real(kind=real_kind) :: div(4,4)

    ! Local
    integer i,j,m,n

    real(kind=real_kind) ::  vtemp(4,4,2)
    real(kind=real_kind) ::  ggtemp(4,4,2)
    real(kind=real_kind) ::  gtemp(4,4,2)
    real(kind=real_kind) ::  psi(4,4)
    real(kind=real_kind) :: xtmp

    ! latlon- > contra
    do j=1,4
       do i=1,4
          vtemp(i,j,1)=(elem_Dinv(1,1,i,j)*v(i,j,1) + elem_Dinv(1,2,i,j)*v(i,j,2))
          vtemp(i,j,2)=(elem_Dinv(2,1,i,j)*v(i,j,1) + elem_Dinv(2,2,i,j)*v(i,j,2))
       enddo
    enddo

    do n=1,4
       do m=1,4
          div(m,n)=0
          do j=1,4
             div(m,n)=div(m,n)-(elem_spheremp(j,n)*vtemp(j,n,1)*deriv_Dvv(m,j) &
                               +elem_spheremp(m,j)*vtemp(m,j,2)*deriv_Dvv(n,j)) &
                              * my_rrearth
          enddo
       end do
    end do
  end function

subroutine my_laplace_sphere_wk(s, deriv_Dvv, elem_spheremp, elem_Dinv, elem_variable_hyperviscosity, elem_tensorVisc, laplace, my_rrearth, var_coef)

  real(kind=8), intent(in) :: s(4,4)
  real(kind=8), dimension(4,4), intent(in) :: deriv_Dvv
  real(kind=8), dimension(4,4), intent(in) :: elem_spheremp
  real(kind=8), dimension(2,2,4,4), intent(in) :: elem_Dinv
  real(kind=8), dimension(4,4), intent(in) :: elem_variable_hyperviscosity
  real(kind=8), dimension(2,2,4,4), intent(in) :: elem_tensorVisc
  real(kind=8), intent(in) :: my_rrearth
  real(kind=8), intent(out) :: laplace(4,4)
  logical :: var_coef

  ! Local
  real(kind=8) :: grads(4,4,2), oldgrads(4,4,2)
  integer i,j

  call my_gradient_sphere_vis(s,deriv_Dvv,elem_Dinv,grads,my_rrearth)

  if (var_coef) then
     if (hypervis_scaling==0 ) then
        ! const or variable viscosity, (1) or (2)
        grads(:,:,1) = grads(:,:,1)*elem_variable_hyperviscosity(:,:)
        grads(:,:,2) = grads(:,:,2)*elem_variable_hyperviscosity(:,:)
     else
        ! tensor hv, (3)
        oldgrads=grads
        do j=1,4
           do i=1,4
              grads(i,j,1) = sum(oldgrads(i,j,:)*elem_tensorVisc(1,:,i,j))
              grads(i,j,2) = sum(oldgrads(i,j,:)*elem_tensorVisc(2,:,i,j))
           end do
        end do
     endif
  endif

  laplace = my_divergence_sphere_wk(grads, elem_Dinv, elem_spheremp, deriv_Dvv, my_rrearth)
end subroutine

  function my_divergence_sphere(v, deriv_Dvv, elem_Dinv, elem_metdet, elem_rmetdet, my_rrearth) result(div)
    real(kind=8), intent(in) :: v(4,4,2)  ! in lat-lon coordinates
    real(kind=8), dimension(4, 4), intent(in) :: deriv_Dvv
    real(kind=8), dimension(2,2,4, 4), intent(in) :: elem_Dinv
    real(kind=8), dimension(4, 4), intent(in) :: elem_rmetdet
    real(kind=8), dimension(4, 4), intent(in) :: elem_metdet
    real(kind=8), intent(in) :: my_rrearth
    real(kind=8) :: div(4,4)

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

  end function

  function my_vorticity_sphere(v,deriv_Dvv, elem_D, elem_rmetdet, my_rrearth) result(vort)
    real(kind=8), intent(in) :: v(4,4,2)
    real(kind=8), dimension(4,4), intent(in) :: deriv_Dvv
    real(kind=8), dimension(2,2,4,4), intent(in) :: elem_D
    real(kind=8), dimension(4,4), intent(in) :: elem_rmetdet
    real(kind=8), intent(in) :: my_rrearth

    real(kind=8) :: vort(4,4)

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

  end function

  function my_gradient_sphere_wk_testcov(my_rrearth, s,deriv_Dvv, elem_D, elem_metdet, elem_mp, elem_metinv) result(ds)
    real(kind=8), intent(in) :: s(4,4)

    real(kind=8), dimension(4,4), intent(in) :: deriv_Dvv
    real(kind=8), dimension(2,2,4,4), intent(in) :: elem_D

    real(kind=8), dimension(4,4), intent(in) :: elem_metdet
    real(kind=8), dimension(4,4), intent(in) :: elem_mp
    real(kind=8), dimension(2,2,4,4), intent(in) :: elem_metinv
    real(kind=8), intent(in) :: my_rrearth

    real(kind=8) :: ds(4,4,2)

    integer i,j,l,m,n
    real(kind=8) ::  dscontra(4,4,2)

    dscontra=0
    do n=1,4
       do m=1,4
          do j=1,4
             dscontra(m,n,1)=dscontra(m,n,1)-(&
                  (elem_mp(j,n)*elem_metinv(1,1,m,n)*elem_metdet(m,n)*s(j,n)*deriv_Dvv(m,j) ) +&
                  (elem_mp(m,j)*elem_metinv(2,1,m,n)*elem_metdet(m,n)*s(m,j)*deriv_Dvv(n,j) ) &
                  ) *my_rrearth

             dscontra(m,n,2)=dscontra(m,n,2)-(&
                  (elem_mp(j,n)*elem_metinv(1,2,m,n)*elem_metdet(m,n)*s(j,n)*deriv_Dvv(m,j) ) +&
                  (elem_mp(m,j)*elem_metinv(2,2,m,n)*elem_metdet(m,n)*s(m,j)*deriv_Dvv(n,j) ) &
                  ) *my_rrearth
          enddo
       enddo
    enddo
    ! convert contra -> latlon
    do j=1,4
       do i=1,4
          ds(i,j,1)=(elem_D(1,1,i,j)*dscontra(i,j,1) + elem_D(1,2,i,j)*dscontra(i,j,2))
          ds(i,j,2)=(elem_D(2,1,i,j)*dscontra(i,j,1) + elem_D(2,2,i,j)*dscontra(i,j,2))
       enddo
    enddo

  end function

  function my_curl_sphere_wk_testcov(my_rrearth, s, deriv_Dvv, elem_D, elem_mp) result(ds)
    real(kind=8), intent(in) :: s(4,4)
    real(kind=8), dimension(4,4), intent(in) :: deriv_Dvv
    real(kind=8), dimension(2,2,4,4), intent(in) :: elem_D
    real(kind=8), dimension(4,4), intent(in) :: elem_mp
    real(kind=8), intent(in) :: my_rrearth

    real(kind=8) :: ds(4,4,2)

    integer i,j,l,m,n
    real(kind=8) ::  dscontra(4,4,2)

    dscontra=0
    do n=1,4
       do m=1,4
          do j=1,4
             dscontra(m,n,1)=dscontra(m,n,1)-(elem_mp(m,j)*s(m,j)*deriv_Dvv(n,j) )*my_rrearth
             dscontra(m,n,2)=dscontra(m,n,2)+(elem_mp(j,n)*s(j,n)*deriv_Dvv(m,j) )*my_rrearth
          enddo
       enddo
    enddo

    ! convert contra -> latlon
    do j=1,4
       do i=1,4
          ds(i,j,1)=(elem_D(1,1,i,j)*dscontra(i,j,1) + elem_D(1,2,i,j)*dscontra(i,j,2))
          ds(i,j,2)=(elem_D(2,1,i,j)*dscontra(i,j,1) + elem_D(2,2,i,j)*dscontra(i,j,2))
       enddo
    enddo
    end function

  function my_laplace_sphere_wk_new(my_rrearth, deriv_Dvv,elem_D,elem_metdet,elem_mp,elem_metinv, &
      elem_Dinv,elem_rmetdet,elem_spheremp,elem_variable_hyperviscosity, &
      v,var_coef,nu_ratio) result(laplace)

    real(kind=8), intent(in) :: v(4,4,2)
    logical :: var_coef
    real(kind=8) :: laplace(4,4,2)
    real(kind=8), optional :: nu_ratio
    real(kind=8), intent(in) :: my_rrearth

    real(kind=8), dimension(4,4), intent(in)      :: deriv_Dvv
    real(kind=8), dimension(2,2,4,4), intent(in)  :: elem_D
    real(kind=8), dimension(4,4), intent(in)      :: elem_metdet
    real(kind=8), dimension(4,4), intent(in)      :: elem_mp
    real(kind=8), dimension(2,2,4,4), intent(in)  :: elem_metinv
    real(kind=8), dimension(2,2,4, 4), intent(in) :: elem_Dinv
    real(kind=8), dimension(4, 4), intent(in)     :: elem_rmetdet
    real(kind=8), dimension(4, 4), intent(in)     :: elem_spheremp
    real(kind=8), dimension(4, 4), intent(in)     :: elem_variable_hyperviscosity

    ! Local
    integer i,j,l,m,n
    real(kind=8) :: vor(4,4),div(4,4)
    real(kind=8) :: v1,v2,div1,div2,vor1,vor2,phi_x,phi_y

    div=my_divergence_sphere(v,deriv_Dvv,elem_Dinv, elem_metdet, elem_rmetdet, my_rrearth)
    vor= my_vorticity_sphere(v,deriv_Dvv, elem_D, elem_rmetdet, my_rrearth)

    if (var_coef) then
       div = div*elem_variable_hyperviscosity(:,:)
       vor = vor*elem_variable_hyperviscosity(:,:)
    end if
    if (present(nu_ratio)) div = nu_ratio*div

    laplace = my_gradient_sphere_wk_testcov(my_rrearth, div,deriv_Dvv,elem_D, elem_metdet, elem_mp, elem_metinv) - &
         my_curl_sphere_wk_testcov(my_rrearth, vor,deriv_Dvv, elem_D, elem_mp)

    do n=1,4
       do m=1,4
          ! add in correction so we dont damp rigid rotation
          laplace(m,n,1)=laplace(m,n,1) + 2*elem_spheremp(m,n)*v(m,n,1)*(my_rrearth**2)
          laplace(m,n,2)=laplace(m,n,2) + 2*elem_spheremp(m,n)*v(m,n,2)*(my_rrearth**2)
       enddo
    enddo
  end function

subroutine my_unpack_after(nets, nete, my_elem, deriv, ptens, vtens, dptens,  nu_ratio, my_rrearth)
  implicit none

  integer, parameter :: my_np = 4
  integer, parameter :: my_nlev = constLev 

  integer, intent(in) :: nets, nete
  real (kind=8), intent(in) ::  nu_ratio
  real(kind=8), intent(in) :: my_rrearth
  type (element_t), intent(in) :: my_elem(nets:nete)
  real (kind=8), intent(inout) :: ptens(4,4,constLev,nets:nete)
  real (kind=8), intent(inout) :: vtens(4,4,2,constLev,nets:nete)
  real (kind=8), intent(inout) :: dptens(4,4,constLev,nets:nete)
  type (derivative_t)  , intent(in) :: deriv

  !pointer
  real(kind=8), dimension(my_np, my_np) :: elem_rspheremp
  pointer(elem_rspheremp_ptr, elem_rspheremp)

  real(kind=8), dimension(np,np) :: deriv_Dvv
  pointer(deriv_Dvv_ptr, deriv_Dvv)

  real(kind=8), dimension(np,np) :: elem_spheremp
  pointer(elem_spheremp_ptr, elem_spheremp)

  real(kind=8), dimension(2,2,np,np) :: elem_Dinv
  pointer(elem_Dinv_ptr, elem_Dinv)

  real(kind=8), dimension(np,np) :: elem_variable_hyperviscosity
  pointer(elem_variable_hyperviscosity_ptr, elem_variable_hyperviscosity)

  real(kind=8), dimension(2,2,np,np) :: elem_tensorVisc
  pointer(elem_tensorVisc_ptr, elem_tensorVisc)

  real(kind=8), dimension(2,2,4,4)  :: elem_D
  pointer(elem_D_ptr, elem_D)

  real(kind=8), dimension(4,4)     :: elem_metdet
  pointer(elem_metdet_ptr, elem_metdet)

  real(kind=8), dimension(4,4)      :: elem_mp
  pointer(elem_mp_ptr, elem_mp)

  real(kind=8), dimension(2,2,4,4)  :: elem_metinv
  pointer(elem_metinv_ptr, elem_metinv)

  real(kind=8), dimension(4, 4)    :: elem_rmetdet
  pointer(elem_rmetdet_ptr, elem_rmetdet)

  !local
  real(kind=8), dimension(my_np, my_np) :: tmp
  real(kind=8), dimension(my_np, my_np, 2) :: v
  integer :: ie, k

  deriv_Dvv_ptr = loc(deriv%Dvv)

  !$ACC PARALLEL LOOP collapse(2) local(tmp, v) copy(ptens, dptens,vtens) copyin(deriv_Dvv) annotate(entire(deriv_Dvv))
  do ie = nets, nete
    do k=1,my_nlev
      elem_rspheremp_ptr               = loc(my_elem(ie)%rspheremp(:,:))
      elem_spheremp_ptr                = loc(my_elem(ie)%spheremp)
      elem_Dinv_ptr                    = loc(my_elem(ie)%Dinv)
      elem_variable_hyperviscosity_ptr = loc(my_elem(ie)%variable_hyperviscosity)
      elem_tensorVisc_ptr              = loc(my_elem(ie)%tensorVisc)
      elem_D_ptr                       = loc(my_elem(ie)%D)
      elem_metdet_ptr                  = loc(my_elem(ie)%metdet)
      elem_mp_ptr                      = loc(my_elem(ie)%mp)
      elem_metinv_ptr                  = loc(my_elem(ie)%metinv)
      elem_rmetdet_ptr                 = loc(my_elem(ie)%rmetdet)

      !$ACC DATA copyin(elem_rspheremp, elem_Dinv, elem_spheremp, elem_variable_hyperviscosity, elem_tensorVisc,elem_D,elem_metdet,elem_mp,elem_metinv,elem_rmetdet)
      tmp(:,:)=elem_rspheremp(:,:)*ptens(:,:,k,ie)

      call my_laplace_sphere_wk(tmp, deriv_Dvv, elem_spheremp, elem_Dinv, elem_variable_hyperviscosity, &
                                elem_tensorVisc, ptens(:,:,k,ie), my_rrearth, var_coef=.true. )

      tmp(:,:)=elem_rspheremp(:,:)*dptens(:,:,k,ie)
      call my_laplace_sphere_wk(tmp, deriv_Dvv, elem_spheremp, elem_Dinv, elem_variable_hyperviscosity, &
                                elem_tensorVisc, dptens(:,:,k,ie), my_rrearth, var_coef=.true. )

      v(:,:,1)=elem_rspheremp(:,:)*vtens(:,:,1,k,ie)
      v(:,:,2)=elem_rspheremp(:,:)*vtens(:,:,2,k,ie)
      vtens(:,:,:,k,ie) =  my_laplace_sphere_wk_new(my_rrearth, deriv_Dvv, elem_D,elem_metdet,elem_mp,elem_metinv, &
                                              elem_Dinv,elem_rmetdet,elem_spheremp,elem_variable_hyperviscosity, &
                                              v(:,:,:),var_coef=.true.,nu_ratio=nu_ratio)
      !$ACC END DATA
    end do
  end do
  !$ACC END PARALLEL LOOP

end subroutine

subroutine my_pack_before(nt, nets, nete, my_elem, deriv, ptens, vtens, dptens,  nu_ratio, my_rrearth)
  integer, parameter :: my_nlev = constLev 
  integer, parameter :: my_np   = 4

  integer, intent(in) :: nt, nets, nete
  type (element_t), intent(in) :: my_elem(nets:nete)
  real (kind=8), intent(in) ::  nu_ratio
  real(kind=8), intent(in) :: my_rrearth
  real (kind=8), intent(inout) :: ptens(4,4,constLev,nets:nete)
  real (kind=8), intent(inout) :: vtens(4,4,2,constLev,nets:nete)
  real (kind=8), intent(inout) :: dptens(4,4,constLev,nets:nete)
  type (derivative_t)  , intent(in) :: deriv

  !pointer
  real(kind = 8), dimension(my_np, my_np) :: elem_state_T
  pointer(elem_state_T_ptr, elem_state_T)

  real(kind = 8), dimension(my_np, my_np) :: elem_state_dp3d
  pointer(elem_state_dp3d_ptr, elem_state_dp3d)

  real(kind = 8), dimension(my_np, my_np, 2) :: elem_state_v
  pointer(elem_state_v_ptr, elem_state_v)

  real(kind=8), dimension(np,np) :: deriv_Dvv
  pointer(deriv_Dvv_ptr, deriv_Dvv)

  real(kind=8), dimension(np,np) :: elem_spheremp
  pointer(elem_spheremp_ptr, elem_spheremp)

  real(kind=8), dimension(2,2,np,np) :: elem_Dinv
  pointer(elem_Dinv_ptr, elem_Dinv)

  real(kind=8), dimension(np,np) :: elem_variable_hyperviscosity
  pointer(elem_variable_hyperviscosity_ptr, elem_variable_hyperviscosity)

  real(kind=8), dimension(2,2,np,np) :: elem_tensorVisc
  pointer(elem_tensorVisc_ptr, elem_tensorVisc)

  real(kind=8), dimension(2,2,4,4)  :: elem_D
  pointer(elem_D_ptr, elem_D)

  real(kind=8), dimension(4,4)     :: elem_metdet
  pointer(elem_metdet_ptr, elem_metdet)

  real(kind=8), dimension(4,4)      :: elem_mp
  pointer(elem_mp_ptr, elem_mp)

  real(kind=8), dimension(2,2,4,4)  :: elem_metinv
  pointer(elem_metinv_ptr, elem_metinv)

  real(kind=8), dimension(4, 4)    :: elem_rmetdet
  pointer(elem_rmetdet_ptr, elem_rmetdet)

  !local
  integer :: ie, k

  deriv_Dvv_ptr = loc(deriv%Dvv)

  !$ACC PARALLEL LOOP collapse(2)  copy(ptens, dptens,vtens) copyin(deriv_Dvv) annotate(entire(deriv_Dvv))
  do ie=nets,nete
    do k=1,my_nlev
      elem_state_T_ptr = loc(my_elem(ie)%state%T(:,:,k,nt))
      elem_state_dp3d_ptr = loc(my_elem(ie)%state%dp3d(:,:,k,nt))
      elem_state_v_ptr = loc(my_elem(ie)%state%v(:,:,:,k,nt))
      elem_spheremp_ptr                = loc(my_elem(ie)%spheremp)
      elem_Dinv_ptr                    = loc(my_elem(ie)%Dinv)
      elem_variable_hyperviscosity_ptr = loc(my_elem(ie)%variable_hyperviscosity)
      elem_tensorVisc_ptr              = loc(my_elem(ie)%tensorVisc)
      elem_D_ptr                       = loc(my_elem(ie)%D)
      elem_metdet_ptr                  = loc(my_elem(ie)%metdet)
      elem_mp_ptr                      = loc(my_elem(ie)%mp)
      elem_metinv_ptr                  = loc(my_elem(ie)%metinv)
      elem_rmetdet_ptr                 = loc(my_elem(ie)%rmetdet)

      !$ACC DATA copyin(elem_state_T, elem_state_dp3d, elem_state_v, elem_Dinv, elem_spheremp, elem_variable_hyperviscosity, elem_tensorVisc,elem_D,elem_metdet,elem_mp,elem_metinv,elem_rmetdet)
      call my_laplace_sphere_wk(elem_state_T(:,:), deriv_Dvv(:,:), elem_spheremp(:,:), elem_Dinv, elem_variable_hyperviscosity, &
                                elem_tensorVisc, ptens(:,:,k,ie), my_rrearth, var_coef=.true. )

      call my_laplace_sphere_wk(elem_state_dp3d(:,:), deriv_Dvv, elem_spheremp, elem_Dinv, elem_variable_hyperviscosity, &
                                elem_tensorVisc, dptens(:,:,k,ie), my_rrearth, var_coef=.true. )

      vtens(:,:,:,k,ie) =  my_laplace_sphere_wk_new(my_rrearth, deriv_Dvv(:,:), elem_D(:,:,:,:),elem_metdet(:,:),elem_mp(:,:),elem_metinv(:,:,:,:), &
                                              elem_Dinv(:,:,:,:),elem_rmetdet(:,:),elem_spheremp(:,:),elem_variable_hyperviscosity(:,:), &
                                              elem_state_v(:,:,:),var_coef=.true.,nu_ratio=nu_ratio)
      !$ACC END DATA
    enddo
  end do
  !$ACC END PARALLEL LOOP
end subroutine


#ifdef _PRIM
subroutine biharmonic_wk_dp3d(elem,dptens,ptens,vtens,deriv,edge3,hybrid,nt,nets,nete,nu_ratio)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! compute weak biharmonic operator
!    input:  h,v (stored in elem()%, in lat-lon coordinates
!    output: ptens,vtens  overwritten with weak biharmonic of h,v (output in lat-lon coordinates)
! conghui
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use physical_constants, only : rrearth
use perf_mod, only: t_startf, t_stopf, t_barrierf ! _EXTERNAL
type (hybrid_t)      , intent(in) :: hybrid
type (element_t)     , intent(inout), target :: elem(:)
integer :: nt,nets,nete
real (kind=real_kind), dimension(np,np,2,nlev,nets:nete)  :: vtens
real (kind=real_kind), dimension(np,np,nlev,nets:nete) :: ptens,dptens
type (EdgeBuffer_t)  , intent(inout) :: edge3
type (derivative_t)  , intent(in) :: deriv
real (kind=real_kind) ::  nu_ratio

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
! local
integer :: k,kptr,ie
real (kind=real_kind), dimension(:,:), pointer :: rspheremv
real (kind=real_kind), dimension(np,np) :: tmp
real (kind=real_kind), dimension(np,np,2) :: v

integer (kind=8) :: count_start, count_stop, count_rate, count_max
integer(kind=8), dimension(3,nets:nete)      :: pack_elem_array

integer, dimension(8) :: elem_desc_putmapP
pointer(elem_desc_putmapP_ptr, elem_desc_putmapP)
  
logical, dimension(8) :: elem_desc_reverse
pointer(elem_desc_reverse_ptr, elem_desc_reverse)

integer(kind=8), dimension(20,4,nets:nete) :: pack_buf_array

    !call system_clock(count_start, count_rate, count_max)
    call my_pack_before(nt, nets, nete, elem, deriv, ptens, vtens, dptens,  nu_ratio, rrearth)
    !call system_clock(count_stop, count_rate, count_max)
    !write(*,*) 'normal count = ', (count_stop - count_start)
    call t_startf("biharmonic pack")
    do ie =nets, nete
            pack_buf_array(1,1,ie) = loc(edge3%buf(1,elem(ie)%desc%putmapP(north)+1))
            pack_buf_array(2,1,ie) = loc(edge3%buf(1,elem(ie)%desc%putmapP(north)+2))
            pack_buf_array(3,1,ie) = loc(edge3%buf(1,elem(ie)%desc%putmapP(north)+3))
            pack_buf_array(4,1,ie) = loc(edge3%buf(1,elem(ie)%desc%putmapP(north)+4))
            pack_buf_array(5,1,ie) = loc(edge3%buf(1,elem(ie)%desc%putmapP(5)+1))
            pack_buf_array(6,1,ie) = loc(edge3%buf(1,elem(ie)%desc%putmapP(6)+1))
            pack_buf_array(7,1,ie) = loc(edge3%buf(1,elem(ie)%desc%putmapP(7)+1))
            pack_buf_array(8,1,ie) = loc(edge3%buf(1,elem(ie)%desc%putmapP(8)+1))
            pack_buf_array(9,1,ie) = loc(edge3%buf(1,elem(ie)%desc%putmapP(south)+1))
            pack_buf_array(10,1,ie) = loc(edge3%buf(1,elem(ie)%desc%putmapP(south)+2))
            pack_buf_array(11,1,ie) = loc(edge3%buf(1,elem(ie)%desc%putmapP(south)+3))
            pack_buf_array(12,1,ie) = loc(edge3%buf(1,elem(ie)%desc%putmapP(south)+4))
            pack_buf_array(13,1,ie) = loc(edge3%buf(1,elem(ie)%desc%putmapP(east)+1))
            pack_buf_array(14,1,ie) = loc(edge3%buf(1,elem(ie)%desc%putmapP(east)+2))
            pack_buf_array(15,1,ie) = loc(edge3%buf(1,elem(ie)%desc%putmapP(east)+3))
            pack_buf_array(16,1,ie) = loc(edge3%buf(1,elem(ie)%desc%putmapP(east)+4))
            pack_buf_array(17,1,ie) = loc(edge3%buf(1,elem(ie)%desc%putmapP(west)+1))
            pack_buf_array(18,1,ie) = loc(edge3%buf(1,elem(ie)%desc%putmapP(west)+2))
            pack_buf_array(19,1,ie) = loc(edge3%buf(1,elem(ie)%desc%putmapP(west)+3))
            pack_buf_array(20,1,ie) = loc(edge3%buf(1,elem(ie)%desc%putmapP(west)+4))
            
            pack_buf_array(1,2,ie) = loc(edge3%buf(nlev+1,elem(ie)%desc%putmapP(north)+1))
            pack_buf_array(2,2,ie) = loc(edge3%buf(nlev+1,elem(ie)%desc%putmapP(north)+2))
            pack_buf_array(3,2,ie) = loc(edge3%buf(nlev+1,elem(ie)%desc%putmapP(north)+3))
            pack_buf_array(4,2,ie) = loc(edge3%buf(nlev+1,elem(ie)%desc%putmapP(north)+4))
            pack_buf_array(5,2,ie) = loc(edge3%buf(nlev+1,elem(ie)%desc%putmapP(5)+1))
            pack_buf_array(6,2,ie) = loc(edge3%buf(nlev+1,elem(ie)%desc%putmapP(6)+1))
            pack_buf_array(7,2,ie) = loc(edge3%buf(nlev+1,elem(ie)%desc%putmapP(7)+1))
            pack_buf_array(8,2,ie) = loc(edge3%buf(nlev+1,elem(ie)%desc%putmapP(8)+1))
            pack_buf_array(9,2,ie) = loc(edge3%buf(nlev+1,elem(ie)%desc%putmapP(south)+1))
            pack_buf_array(10,2,ie) = loc(edge3%buf(nlev+1,elem(ie)%desc%putmapP(south)+2))
            pack_buf_array(11,2,ie) = loc(edge3%buf(nlev+1,elem(ie)%desc%putmapP(south)+3))
            pack_buf_array(12,2,ie) = loc(edge3%buf(nlev+1,elem(ie)%desc%putmapP(south)+4))
            pack_buf_array(13,2,ie) = loc(edge3%buf(nlev+1,elem(ie)%desc%putmapP(east)+1))
            pack_buf_array(14,2,ie) = loc(edge3%buf(nlev+1,elem(ie)%desc%putmapP(east)+2))
            pack_buf_array(15,2,ie) = loc(edge3%buf(nlev+1,elem(ie)%desc%putmapP(east)+3))
            pack_buf_array(16,2,ie) = loc(edge3%buf(nlev+1,elem(ie)%desc%putmapP(east)+4))
            pack_buf_array(17,2,ie) = loc(edge3%buf(nlev+1,elem(ie)%desc%putmapP(west)+1))
            pack_buf_array(18,2,ie) = loc(edge3%buf(nlev+1,elem(ie)%desc%putmapP(west)+2))
            pack_buf_array(19,2,ie) = loc(edge3%buf(nlev+1,elem(ie)%desc%putmapP(west)+3))
            pack_buf_array(20,2,ie) = loc(edge3%buf(nlev+1,elem(ie)%desc%putmapP(west)+4))
           
            pack_buf_array(1,3,ie) = loc(edge3%buf(2*nlev+1,elem(ie)%desc%putmapP(north)+1))
            pack_buf_array(2,3,ie) = loc(edge3%buf(2*nlev+1,elem(ie)%desc%putmapP(north)+2))
            pack_buf_array(3,3,ie) = loc(edge3%buf(2*nlev+1,elem(ie)%desc%putmapP(north)+3))
            pack_buf_array(4,3,ie) = loc(edge3%buf(2*nlev+1,elem(ie)%desc%putmapP(north)+4))
            pack_buf_array(5,3,ie) = loc(edge3%buf(2*nlev+1,elem(ie)%desc%putmapP(5)+1))
            pack_buf_array(6,3,ie) = loc(edge3%buf(2*nlev+1,elem(ie)%desc%putmapP(6)+1))
            pack_buf_array(7,3,ie) = loc(edge3%buf(2*nlev+1,elem(ie)%desc%putmapP(7)+1))
            pack_buf_array(8,3,ie) = loc(edge3%buf(2*nlev+1,elem(ie)%desc%putmapP(8)+1))
            pack_buf_array(9,3,ie) = loc(edge3%buf(2*nlev+1,elem(ie)%desc%putmapP(south)+1))
            pack_buf_array(10,3,ie) = loc(edge3%buf(2*nlev+1,elem(ie)%desc%putmapP(south)+2))
            pack_buf_array(11,3,ie) = loc(edge3%buf(2*nlev+1,elem(ie)%desc%putmapP(south)+3))
            pack_buf_array(12,3,ie) = loc(edge3%buf(2*nlev+1,elem(ie)%desc%putmapP(south)+4))
            pack_buf_array(13,3,ie) = loc(edge3%buf(2*nlev+1,elem(ie)%desc%putmapP(east)+1))
            pack_buf_array(14,3,ie) = loc(edge3%buf(2*nlev+1,elem(ie)%desc%putmapP(east)+2))
            pack_buf_array(15,3,ie) = loc(edge3%buf(2*nlev+1,elem(ie)%desc%putmapP(east)+3))
            pack_buf_array(16,3,ie) = loc(edge3%buf(2*nlev+1,elem(ie)%desc%putmapP(east)+4))
            pack_buf_array(17,3,ie) = loc(edge3%buf(2*nlev+1,elem(ie)%desc%putmapP(west)+1))
            pack_buf_array(18,3,ie) = loc(edge3%buf(2*nlev+1,elem(ie)%desc%putmapP(west)+2))
            pack_buf_array(19,3,ie) = loc(edge3%buf(2*nlev+1,elem(ie)%desc%putmapP(west)+3))
            pack_buf_array(20,3,ie) = loc(edge3%buf(2*nlev+1,elem(ie)%desc%putmapP(west)+4))

            pack_buf_array(1,4,ie) = loc(edge3%buf(3*nlev+1,elem(ie)%desc%putmapP(north)+1))
            pack_buf_array(2,4,ie) = loc(edge3%buf(3*nlev+1,elem(ie)%desc%putmapP(north)+2))
            pack_buf_array(3,4,ie) = loc(edge3%buf(3*nlev+1,elem(ie)%desc%putmapP(north)+3))
            pack_buf_array(4,4,ie) = loc(edge3%buf(3*nlev+1,elem(ie)%desc%putmapP(north)+4))
            pack_buf_array(5,4,ie) = loc(edge3%buf(3*nlev+1,elem(ie)%desc%putmapP(5)+1))
            pack_buf_array(6,4,ie) = loc(edge3%buf(3*nlev+1,elem(ie)%desc%putmapP(6)+1))
            pack_buf_array(7,4,ie) = loc(edge3%buf(3*nlev+1,elem(ie)%desc%putmapP(7)+1))
            pack_buf_array(8,4,ie) = loc(edge3%buf(3*nlev+1,elem(ie)%desc%putmapP(8)+1))
            pack_buf_array(9,4,ie) = loc(edge3%buf(3*nlev+1,elem(ie)%desc%putmapP(south)+1))
            pack_buf_array(10,4,ie) = loc(edge3%buf(3*nlev+1,elem(ie)%desc%putmapP(south)+2))
            pack_buf_array(11,4,ie) = loc(edge3%buf(3*nlev+1,elem(ie)%desc%putmapP(south)+3))
            pack_buf_array(12,4,ie) = loc(edge3%buf(3*nlev+1,elem(ie)%desc%putmapP(south)+4))
            pack_buf_array(13,4,ie) = loc(edge3%buf(3*nlev+1,elem(ie)%desc%putmapP(east)+1))
            pack_buf_array(14,4,ie) = loc(edge3%buf(3*nlev+1,elem(ie)%desc%putmapP(east)+2))
            pack_buf_array(15,4,ie) = loc(edge3%buf(3*nlev+1,elem(ie)%desc%putmapP(east)+3))
            pack_buf_array(16,4,ie) = loc(edge3%buf(3*nlev+1,elem(ie)%desc%putmapP(east)+4))
            pack_buf_array(17,4,ie) = loc(edge3%buf(3*nlev+1,elem(ie)%desc%putmapP(west)+1))
            pack_buf_array(18,4,ie) = loc(edge3%buf(3*nlev+1,elem(ie)%desc%putmapP(west)+2))
            pack_buf_array(19,4,ie) = loc(edge3%buf(3*nlev+1,elem(ie)%desc%putmapP(west)+3))
            pack_buf_array(20,4,ie) = loc(edge3%buf(3*nlev+1,elem(ie)%desc%putmapP(west)+4))
    
            pack_elem_array(1,ie) = loc(elem(ie)%desc%reverse)
            pack_elem_array(2,ie) = loc(elem(ie)%desc%putmapP)
        enddo
         !$ACC parallel loop copyin(pack_elem_array, pack_buf_array, ptens)
    do ie = nets, nete
            elem_desc_reverse_ptr   = pack_elem_array(1,ie)
            elem_desc_putmapP_ptr         = pack_elem_array(2,ie)
            edge_buf_in_1_ptr = pack_buf_array(1,1,ie)
            edge_buf_in_2_ptr = pack_buf_array(2,1,ie)
            edge_buf_in_3_ptr = pack_buf_array(3,1,ie)
            edge_buf_in_4_ptr = pack_buf_array(4,1,ie)
            edge_buf_5_ptr  = pack_buf_array(5,1,ie)
            edge_buf_6_ptr  = pack_buf_array(6,1,ie)
            edge_buf_7_ptr  = pack_buf_array(7,1,ie)
            edge_buf_8_ptr  = pack_buf_array(8,1,ie)
            edge_buf_is_1_ptr = pack_buf_array(9,1,ie)
            edge_buf_is_2_ptr = pack_buf_array(10,1,ie)
            edge_buf_is_3_ptr = pack_buf_array(11,1,ie)
            edge_buf_is_4_ptr = pack_buf_array(12,1,ie)
            edge_buf_ie_1_ptr = pack_buf_array(13,1,ie)
            edge_buf_ie_2_ptr = pack_buf_array(14,1,ie)
            edge_buf_ie_3_ptr = pack_buf_array(15,1,ie)
            edge_buf_ie_4_ptr = pack_buf_array(16,1,ie)
            edge_buf_iw_1_ptr = pack_buf_array(17,1,ie)
            edge_buf_iw_2_ptr = pack_buf_array(18,1,ie)
            edge_buf_iw_3_ptr = pack_buf_array(19,1,ie)
            edge_buf_iw_4_ptr = pack_buf_array(20,1,ie)
           !$ACC data  copyin(elem_desc_putmapP, elem_desc_reverse) copyout(edge_buf_in_1, edge_buf_in_2, edge_buf_in_3, edge_buf_in_4, edge_buf_5, edge_buf_6, edge_buf_7, edge_buf_8, edge_buf_is_1, edge_buf_is_2, edge_buf_is_3, edge_buf_is_4, edge_buf_ie_1, edge_buf_ie_2, edge_buf_ie_3, edge_buf_ie_4, edge_buf_iw_1, edge_buf_iw_2, edge_buf_iw_3, edge_buf_iw_4)
         call my_edgeVpack_acc(ptens(:,:,:,ie) , nlev , 0 , elem_desc_putmapP(:), &
             elem_desc_reverse, &
             edge_buf_5, edge_buf_6, edge_buf_7, edge_buf_8, &
             edge_buf_in_1, edge_buf_in_2, edge_buf_in_3, edge_buf_in_4, &
             edge_buf_is_1, edge_buf_is_2, edge_buf_is_3, edge_buf_is_4, &
             edge_buf_ie_1, edge_buf_ie_2, edge_buf_ie_3, edge_buf_ie_4, &
             edge_buf_iw_1, edge_buf_iw_2, edge_buf_iw_3, edge_buf_iw_4)
             !$ACC END DATA
      !kptr=0
      !call edgeVpack(edge3, ptens(1,1,1,ie),nlev,kptr,elem(ie)%desc)
    enddo
   !$ACC END parallel loop
         !$ACC parallel loop copyin(pack_elem_array, pack_buf_array, vtens)
         do ie=nets, nete  

            elem_desc_reverse_ptr   = pack_elem_array(1,ie)
            elem_desc_putmapP_ptr         = pack_elem_array(2,ie)
            edge_buf_in_1_ptr = pack_buf_array(1,2,ie)
            edge_buf_in_2_ptr = pack_buf_array(2,2,ie)
            edge_buf_in_3_ptr = pack_buf_array(3,2,ie)
            edge_buf_in_4_ptr = pack_buf_array(4,2,ie)
            edge_buf_5_ptr  = pack_buf_array(5,2,ie)
            edge_buf_6_ptr  = pack_buf_array(6,2,ie)
            edge_buf_7_ptr  = pack_buf_array(7,2,ie)
            edge_buf_8_ptr  = pack_buf_array(8,2,ie)
            edge_buf_is_1_ptr = pack_buf_array(9,2,ie)
            edge_buf_is_2_ptr = pack_buf_array(10,2,ie)
            edge_buf_is_3_ptr = pack_buf_array(11,2,ie)
            edge_buf_is_4_ptr = pack_buf_array(12,2,ie)
            edge_buf_ie_1_ptr = pack_buf_array(13,2,ie)
            edge_buf_ie_2_ptr = pack_buf_array(14,2,ie)
            edge_buf_ie_3_ptr = pack_buf_array(15,2,ie)
            edge_buf_ie_4_ptr = pack_buf_array(16,2,ie)
            edge_buf_iw_1_ptr = pack_buf_array(17,2,ie)
            edge_buf_iw_2_ptr = pack_buf_array(18,2,ie)
            edge_buf_iw_3_ptr = pack_buf_array(19,2,ie)
            edge_buf_iw_4_ptr = pack_buf_array(20,2,ie)
           !$ACC data  copyin(elem_desc_putmapP, elem_desc_reverse) copyout(edge_buf_in_1, edge_buf_in_2, edge_buf_in_3, edge_buf_in_4, edge_buf_5, edge_buf_6, edge_buf_7, edge_buf_8,edge_buf_is_1, edge_buf_is_2, edge_buf_is_3, edge_buf_is_4, edge_buf_ie_1, edge_buf_ie_2, edge_buf_ie_3, edge_buf_ie_4, edge_buf_iw_1, edge_buf_iw_2, edge_buf_iw_3, edge_buf_iw_4)
             call my_edgeVpack_acc(vtens(:,:,:,1:64,ie) , nlev , nlev+1 , elem_desc_putmapP(:), &
                elem_desc_reverse, &
                edge_buf_5, edge_buf_6, edge_buf_7, edge_buf_8, &
                edge_buf_in_1, edge_buf_in_2, edge_buf_in_3, edge_buf_in_4, &
                edge_buf_is_1, edge_buf_is_2, edge_buf_is_3, edge_buf_is_4, &
                edge_buf_ie_1, edge_buf_ie_2, edge_buf_ie_3, edge_buf_ie_4, &
                edge_buf_iw_1, edge_buf_iw_2, edge_buf_iw_3, edge_buf_iw_4)
                !$ACC end data 
            enddo
            !$ACC end parallel loop
         
         !$ACC parallel loop copyin(pack_elem_array, pack_buf_array, vtens)
        do ie=nets,nete

            elem_desc_reverse_ptr   = pack_elem_array(1,ie)
            elem_desc_putmapP_ptr         = pack_elem_array(2,ie)
            edge_buf_in_1_ptr = pack_buf_array(1,3,ie)
            edge_buf_in_2_ptr = pack_buf_array(2,3,ie)
            edge_buf_in_3_ptr = pack_buf_array(3,3,ie)
            edge_buf_in_4_ptr = pack_buf_array(4,3,ie)
            edge_buf_5_ptr  = pack_buf_array(5,3,ie)
            edge_buf_6_ptr  = pack_buf_array(6,3,ie)
            edge_buf_7_ptr  = pack_buf_array(7,3,ie)
            edge_buf_8_ptr  = pack_buf_array(8,3,ie)
            edge_buf_is_1_ptr = pack_buf_array(9,3,ie)
            edge_buf_is_2_ptr = pack_buf_array(10,3,ie)
            edge_buf_is_3_ptr = pack_buf_array(11,3,ie)
            edge_buf_is_4_ptr = pack_buf_array(12,3,ie)
            edge_buf_ie_1_ptr = pack_buf_array(13,3,ie)
            edge_buf_ie_2_ptr = pack_buf_array(14,3,ie)
            edge_buf_ie_3_ptr = pack_buf_array(15,3,ie)
            edge_buf_ie_4_ptr = pack_buf_array(16,3,ie)
            edge_buf_iw_1_ptr = pack_buf_array(17,3,ie)
            edge_buf_iw_2_ptr = pack_buf_array(18,3,ie)
            edge_buf_iw_3_ptr = pack_buf_array(19,3,ie)
            edge_buf_iw_4_ptr = pack_buf_array(20,3,ie)
           !!$ACC data  copyin(elem_desc_putmapP, elem_desc_reverse) copyout(edge_buf_in_1, edge_buf_in_2, edge_buf_in_3, edge_buf_in_4, edge_buf_5, edge_buf_6, edge_buf_7, edge_buf_8,edge_buf_is_1, edge_buf_is_2, edge_buf_is_3, edge_buf_is_4, edge_buf_ie_1, edge_buf_ie_2, edge_buf_ie_3, edge_buf_ie_4, edge_buf_iw_1, edge_buf_iw_2, edge_buf_iw_3, edge_buf_iw_4)
           !$ACC data  copyin(elem_desc_putmapP, elem_desc_reverse) copyout(edge_buf_in_1, edge_buf_in_2, edge_buf_in_3, edge_buf_in_4, edge_buf_5, edge_buf_6, edge_buf_7, edge_buf_8,edge_buf_is_1, edge_buf_is_2, edge_buf_is_3, edge_buf_is_4, edge_buf_ie_1, edge_buf_ie_2, edge_buf_ie_3, edge_buf_ie_4, edge_buf_iw_1, edge_buf_iw_2, edge_buf_iw_3, edge_buf_iw_4)
             call my_edgeVpack_acc(vtens(:,:,:,65:128,ie) , nlev , 2*nlev , elem_desc_putmapP(:), &
                elem_desc_reverse, &
                edge_buf_5, edge_buf_6, edge_buf_7, edge_buf_8, &
                edge_buf_in_1, edge_buf_in_2, edge_buf_in_3, edge_buf_in_4, &
                edge_buf_is_1, edge_buf_is_2, edge_buf_is_3, edge_buf_is_4, &
                edge_buf_ie_1, edge_buf_ie_2, edge_buf_ie_3, edge_buf_ie_4, &
                edge_buf_iw_1, edge_buf_iw_2, edge_buf_iw_3, edge_buf_iw_4)
                !$ACC END DATA
                !kptr=nlev
           !call edgeVpack(edge3,vtens(:,:,:,:,ie),2*nlev,kptr,elem(ie)%desc)
       enddo
       !$ACC end parallel loop

       !$ACC parallel loop copyin(pack_elem_array, pack_buf_array, dptens)
       do ie=nets, nete 

            elem_desc_reverse_ptr   = pack_elem_array(1,ie)
            elem_desc_putmapP_ptr         = pack_elem_array(2,ie)
            edge_buf_in_1_ptr = pack_buf_array(1,4,ie)
            edge_buf_in_2_ptr = pack_buf_array(2,4,ie)
            edge_buf_in_3_ptr = pack_buf_array(3,4,ie)
            edge_buf_in_4_ptr = pack_buf_array(4,4,ie)
            edge_buf_5_ptr  = pack_buf_array(5,4,ie)
            edge_buf_6_ptr  = pack_buf_array(6,4,ie)
            edge_buf_7_ptr  = pack_buf_array(7,4,ie)
            edge_buf_8_ptr  = pack_buf_array(8,4,ie)
            edge_buf_is_1_ptr = pack_buf_array(9,4,ie)
            edge_buf_is_2_ptr = pack_buf_array(10,4,ie)
            edge_buf_is_3_ptr = pack_buf_array(11,4,ie)
            edge_buf_is_4_ptr = pack_buf_array(12,4,ie)
            edge_buf_ie_1_ptr = pack_buf_array(13,4,ie)
            edge_buf_ie_2_ptr = pack_buf_array(14,4,ie)
            edge_buf_ie_3_ptr = pack_buf_array(15,4,ie)
            edge_buf_ie_4_ptr = pack_buf_array(16,4,ie)
            edge_buf_iw_1_ptr = pack_buf_array(17,4,ie)
            edge_buf_iw_2_ptr = pack_buf_array(18,4,ie)
            edge_buf_iw_3_ptr = pack_buf_array(19,4,ie)
            edge_buf_iw_4_ptr = pack_buf_array(20,4,ie)
           !!$ACC data  copyin(elem_desc_putmapP, elem_desc_reverse) copyout(edge_buf_in_1, edge_buf_in_2, edge_buf_in_3, edge_buf_in_4, edge_buf_5, edge_buf_6, edge_buf_7, edge_buf_8,edge_buf_is_1, edge_buf_is_2, edge_buf_is_3, edge_buf_is_4, edge_buf_ie_1, edge_buf_ie_2, edge_buf_ie_3, edge_buf_ie_4, edge_buf_iw_1, edge_buf_iw_2, edge_buf_iw_3, edge_buf_iw_4)
           !$ACC data  copyin(elem_desc_putmapP, elem_desc_reverse) copyout(edge_buf_in_1, edge_buf_in_2, edge_buf_in_3, edge_buf_in_4, edge_buf_5, edge_buf_6, edge_buf_7, edge_buf_8,edge_buf_is_1, edge_buf_is_2, edge_buf_is_3, edge_buf_is_4, edge_buf_ie_1, edge_buf_ie_2, edge_buf_ie_3, edge_buf_ie_4, edge_buf_iw_1, edge_buf_iw_2, edge_buf_iw_3, edge_buf_iw_4)
            call my_edgeVpack_acc(dptens(:,:,:,ie) , nlev , 3*nlev , elem_desc_putmapP(:), &
              elem_desc_reverse, &
              edge_buf_5, edge_buf_6, edge_buf_7, edge_buf_8, &
              edge_buf_in_1, edge_buf_in_2, edge_buf_in_3, edge_buf_in_4, &
              edge_buf_is_1, edge_buf_is_2, edge_buf_is_3, edge_buf_is_4, &
              edge_buf_ie_1, edge_buf_ie_2, edge_buf_ie_3, edge_buf_ie_4, &
              edge_buf_iw_1, edge_buf_iw_2, edge_buf_iw_3, edge_buf_iw_4)
              !$ACC END DATA
              !kptr=3*nlev  
            !call edgeVpack(edge3,dptens(:,:,:,ie),nlev,kptr,elem(ie)%desc)
        enddo
        !$ACC end parallel loop

    !do ie=nets, nete
    !  kptr=0
    !  call edgeVpack(edge3, ptens(1,1,1,ie),nlev,kptr,elem(ie)%desc)
    !  kptr=nlev
    !  call edgeVpack(edge3, vtens(1,1,1,1,ie),2*nlev,kptr,elem(ie)%desc)
    !  kptr=3*nlev
    !  call edgeVpack(edge3, dptens(1,1,1,ie),nlev,kptr,elem(ie)%desc)
    !enddo
    call t_stopf("biharmonic pack")

    call bndry_exchangeV(hybrid,edge3)

    call my_edgeVunpack_viscosity_all(nets, nete, edge3%nlyr, edge3%nbuf, &
            edge3%buf, south, east, north, west, elem, ptens, vtens, dptens, nlev, &
            swest, max_corner_elem)

    call my_unpack_after(nets, nete, elem, deriv, ptens, vtens, dptens, nu_ratio, rrearth)

end subroutine


subroutine biharmonic_wk_scalar(elem,qtens,deriv,edgeq,hybrid,nets,nete)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! compute weak biharmonic operator
!    input:  qtens = Q
!    output: qtens = weak biharmonic of Q
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
type (hybrid_t)      , intent(in) :: hybrid
type (element_t)     , intent(inout), target :: elem(:)
integer :: nets,nete
real (kind=real_kind), dimension(np,np,nlev,qsize,nets:nete) :: qtens
type (EdgeBuffer_t)  , intent(inout) :: edgeq
type (derivative_t)  , intent(in) :: deriv

! local
integer :: k,kptr,i,j,ie,ic,q
real (kind=real_kind), dimension(np,np) :: lap_p

   do ie=nets,nete
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k, q, lap_p)
#endif
      do q=1,qsize
         do k=1,nlev    !  Potential loop inversion (AAM)
            lap_p(:,:)=qtens(:,:,k,q,ie)
! Original use of qtens on left and right hand sides caused OpenMP errors (AAM)
           qtens(:,:,k,q,ie)=laplace_sphere_wk(lap_p,deriv,elem(ie),var_coef=.true.)
         enddo
      enddo
      call edgeVpack(edgeq, qtens(:,:,:,:,ie),qsize*nlev,0,elem(ie)%desc)
   enddo

   call bndry_exchangeV(hybrid,edgeq)

   do ie=nets,nete
      call edgeVunpack(edgeq, qtens(:,:,:,:,ie),qsize*nlev,0,elem(ie)%desc)

      ! apply inverse mass matrix, then apply laplace again
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k, q, lap_p)
#endif
      do q=1,qsize
      do k=1,nlev    !  Potential loop inversion (AAM)
         lap_p(:,:)=elem(ie)%rspheremp(:,:)*qtens(:,:,k,q,ie)
         qtens(:,:,k,q,ie)=laplace_sphere_wk(lap_p,deriv,elem(ie),var_coef=.true.)
      enddo
      enddo
   enddo
#ifdef DEBUGOMP
#if (! defined ELEMENT_OPENMP)
!$OMP BARRIER
#endif
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end subroutine

subroutine my_edgeVunpack_bihamomic(edge_nlyr, edge_nbuf, &
    desc_getmapP, v, my_vlyr, my_swest, my_max_corner_elem, &
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

subroutine my_edgeVunpack_min_bihamomic(edge_nlyr, edge_nbuf, &
    desc_getmapP, v, my_vlyr, my_swest, my_max_corner_elem, &
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

  ! Local
  logical, parameter :: UseUnroll = .TRUE.
  integer :: i,k,ll

  do k=1,my_vlyr
    v(1  ,1  ,k) = min(v(1  ,1  ,k),edge_buf_is_1(k))
    v(1+1,1  ,k) = min(v(1+1,1  ,k),edge_buf_is_2(k))
    v(3  ,1  ,k) = min(v(3  ,1  ,k),edge_buf_is_3(k))
    v(3+1,1  ,k) = min(v(3+1,1  ,k),edge_buf_is_4(k))

    v(1  ,1  ,k) = min(v(1  ,1  ,k),edge_buf_iw_1(k))
    v(1  ,1+1,k) = min(v(1  ,1+1,k),edge_buf_iw_2(k))
    v(1  ,3  ,k) = min(v(1  ,3  ,k),edge_buf_iw_3(k))
    v(1  ,3+1,k) = min(v(1  ,3+1,k),edge_buf_iw_4(k))

    v(np ,1  ,k) = min(v(np ,1  ,k),edge_buf_ie_1(k))
    v(np ,1+1,k) = min(v(np ,1+1,k),edge_buf_ie_2(k))
    v(np ,3  ,k) = min(v(np ,3  ,k),edge_buf_ie_3(k))
    v(np ,3+1,k) = min(v(np ,3+1,k),edge_buf_ie_4(k))

    v(1  ,np ,k) = min(v(1  ,np ,k),edge_buf_in_1(k))
    v(1+1,np ,k) = min(v(1+1,np ,k),edge_buf_in_2(k))
    v(3  ,np ,k) = min(v(3  ,np ,k),edge_buf_in_3(k))
    v(3+1,np ,k) = min(v(3+1,np ,k),edge_buf_in_4(k))
  end do

! swest
  do ll=my_swest,my_swest+my_max_corner_elem-1
    !write(*,*), 'swest loop, ll = ', ll, 'desc_getmapP(ll) = ', desc_getmapP(ll) ! 5
      if(desc_getmapP(ll) /= -1) then
          do k=1,my_vlyr
              v(1  ,1 ,k)=min(v(1 ,1 ,k),swest_buf(k))
          enddo
      endif
  end do

! SEAST
  do ll=my_swest+my_max_corner_elem,my_swest+2*my_max_corner_elem-1
    !write(*,*), 'seast loop, ll = ', ll, 'desc_getmapP(ll) = ', desc_getmapP(ll) !6
      if(desc_getmapP(ll) /= -1) then
          do k=1,my_vlyr
              v(np ,1 ,k)=min(v(np,1 ,k),seast_buf(k))
          enddo
      endif
  end do

! NEAST
  do ll=my_swest+3*my_max_corner_elem,my_swest+4*my_max_corner_elem-1
    !write(*,*), 'neast loop, ll = ', ll, 'desc_getmapP(ll) = ', desc_getmapP(ll) !8
      if(desc_getmapP(ll) /= -1) then
          do k=1,my_vlyr
              v(np ,np,k)=min(v(np,np,k), neast_buf(k))
          enddo
      endif
  end do

! NWEST
  do ll=my_swest+2*my_max_corner_elem,my_swest+3*my_max_corner_elem-1
    !write(*,*), 'nwest loop, ll = ', ll, 'desc_getmapP(ll) = ', desc_getmapP(ll) !7
      if(desc_getmapP(ll) /= -1) then
          do k=1,my_vlyr
              v(1  ,np,k)=min(v(1 ,np,k),nwest_buf(k))
          enddo
      endif
  end do

end subroutine

subroutine my_edgeVunpack_max_bihamomic(edge_nlyr, edge_nbuf, &
    desc_getmapP, v, my_vlyr, my_swest, my_max_corner_elem, &
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

  ! Local
  logical, parameter :: UseUnroll = .TRUE.
  integer :: i,k,ll

  do k=1,my_vlyr
    v(1  ,1  ,k) = max(v(1  ,1  ,k),edge_buf_is_1(k))
    v(1+1,1  ,k) = max(v(1+1,1  ,k),edge_buf_is_2(k))
    v(3  ,1  ,k) = max(v(3  ,1  ,k),edge_buf_is_3(k))
    v(3+1,1  ,k) = max(v(3+1,1  ,k),edge_buf_is_4(k))

    v(1  ,1  ,k) = max(v(1  ,1  ,k),edge_buf_iw_1(k))
    v(1  ,1+1,k) = max(v(1  ,1+1,k),edge_buf_iw_2(k))
    v(1  ,3  ,k) = max(v(1  ,3  ,k),edge_buf_iw_3(k))
    v(1  ,3+1,k) = max(v(1  ,3+1,k),edge_buf_iw_4(k))

    v(np ,1  ,k) = max(v(np ,1  ,k),edge_buf_ie_1(k))
    v(np ,1+1,k) = max(v(np ,1+1,k),edge_buf_ie_2(k))
    v(np ,3  ,k) = max(v(np ,3  ,k),edge_buf_ie_3(k))
    v(np ,3+1,k) = max(v(np ,3+1,k),edge_buf_ie_4(k))

    v(1  ,np ,k) = max(v(1  ,np ,k),edge_buf_in_1(k))
    v(1+1,np ,k) = max(v(1+1,np ,k),edge_buf_in_2(k))
    v(3  ,np ,k) = max(v(3  ,np ,k),edge_buf_in_3(k))
    v(3+1,np ,k) = max(v(3+1,np ,k),edge_buf_in_4(k))
  end do

! swest
  do ll=my_swest,my_swest+my_max_corner_elem-1
    !write(*,*), 'swest loop, ll = ', ll, 'desc_getmapP(ll) = ', desc_getmapP(ll) ! 5
      if(desc_getmapP(ll) /= -1) then
          do k=1,my_vlyr
              v(1  ,1 ,k)=max(v(1 ,1 ,k),swest_buf(k))
          enddo
      endif
  end do

! SEAST
  do ll=my_swest+my_max_corner_elem,my_swest+2*my_max_corner_elem-1
    !write(*,*), 'seast loop, ll = ', ll, 'desc_getmapP(ll) = ', desc_getmapP(ll) !6
      if(desc_getmapP(ll) /= -1) then
          do k=1,my_vlyr
              v(np ,1 ,k)=max(v(np,1 ,k),seast_buf(k))
          enddo
      endif
  end do

! NEAST
  do ll=my_swest+3*my_max_corner_elem,my_swest+4*my_max_corner_elem-1
    !write(*,*), 'neast loop, ll = ', ll, 'desc_getmapP(ll) = ', desc_getmapP(ll) !8
      if(desc_getmapP(ll) /= -1) then
          do k=1,my_vlyr
              v(np ,np,k)=max(v(np,np,k), neast_buf(k))
          enddo
      endif
  end do

! NWEST
  do ll=my_swest+2*my_max_corner_elem,my_swest+3*my_max_corner_elem-1
    !write(*,*), 'nwest loop, ll = ', ll, 'desc_getmapP(ll) = ', desc_getmapP(ll) !7
      if(desc_getmapP(ll) /= -1) then
          do k=1,my_vlyr
              v(1  ,np,k)=max(v(1 ,np,k),nwest_buf(k))
          enddo
      endif
  end do

end subroutine

subroutine my_biharmonic_wk_scalar_minmax_after_bndry(nets, nete, edge_nlyr, edge_nbuf, &
                        edge_buf, my_south, my_east, my_north, my_west, my_qsize, &
                        my_swest, my_max_corner_elem, my_elem, emin, emax, qtens, deriv, my_rrearth)
  implicit none

  integer, intent(in) ::  edge_nlyr, edge_nbuf, my_south, my_east, my_north, my_west, my_qsize, my_swest, my_max_corner_elem
  integer, intent(in) :: nets, nete
  real(kind=8), dimension(edge_nlyr, edge_nbuf), intent(in) :: edge_buf
  type (element_t), intent(inout) :: my_elem(nets:nete)
  real (kind=8), intent(out), dimension(constLev,my_qsize,nets:nete) :: emin,emax
  real (kind=8), dimension(4,4,constLev,my_qsize,nets:nete), intent(inout) :: qtens
  type (derivative_t)  , intent(in) :: deriv
  real(kind=8), intent(in) :: my_rrearth

  !local
  real (kind=8), dimension(4,4,constLev) :: Qmin
  real (kind=8), dimension(4,4,constLev) :: Qmax
  real (kind=8), dimension(4,4) :: lap_p
  integer :: iee, q, k

  ! pointers
  real(kind=8), dimension(4,4) :: elem_rspheremp
  pointer(elem_rspheremp_ptr, elem_rspheremp)

  real(kind=8), dimension(4,4) :: deriv_Dvv
  pointer(deriv_Dvv_ptr, deriv_Dvv)

  real(kind=8), dimension(4,4) :: elem_spheremp
  pointer(elem_spheremp_ptr, elem_spheremp)

  real(kind=8), dimension(2,2,4,4) :: elem_Dinv
  pointer(elem_Dinv_ptr, elem_Dinv)

  real(kind=8), dimension(4,4) :: elem_variable_hyperviscosity
  pointer(elem_variable_hyperviscosity_ptr, elem_variable_hyperviscosity)

  real(kind=8), dimension(2,2,4,4) :: elem_tensorVisc
  pointer(elem_tensorVisc_ptr, elem_tensorVisc)

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

  deriv_Dvv_ptr = loc(deriv%Dvv)

  !$ACC PARALLEL LOOP collapse(2)  local(lap_p) copy(qtens) copyin(deriv_Dvv) annotate(entire(deriv_Dvv))
  do iee = nets , nete
    do q = 1, my_qsize

      getmapP_ptr        = loc(my_elem(iee)%desc%getmapP)
      is_ptr             = loc(my_elem(iee)%desc%getmapP(my_south))
      ie_ptr             = loc(my_elem(iee)%desc%getmapP(my_east))
      iw_ptr             = loc(my_elem(iee)%desc%getmapP(my_west))
      in_ptr             = loc(my_elem(iee)%desc%getmapP(my_north))
      !$ACC DATA copyin(getmapP,is,ie,iw,in)

      edge_buf_is_1_ptr  = loc(edge_buf((q-1)*constLev+1, is + 1))
      edge_buf_is_2_ptr  = loc(edge_buf((q-1)*constLev+1, is + 2))
      edge_buf_is_3_ptr  = loc(edge_buf((q-1)*constLev+1, is + 3))
      edge_buf_is_4_ptr  = loc(edge_buf((q-1)*constLev+1, is + 4))
      edge_buf_iw_1_ptr  = loc(edge_buf((q-1)*constLev+1, iw + 1))
      edge_buf_iw_2_ptr  = loc(edge_buf((q-1)*constLev+1, iw + 2))
      edge_buf_iw_3_ptr  = loc(edge_buf((q-1)*constLev+1, iw + 3))
      edge_buf_iw_4_ptr  = loc(edge_buf((q-1)*constLev+1, iw + 4))
      edge_buf_ie_1_ptr  = loc(edge_buf((q-1)*constLev+1, ie + 1))
      edge_buf_ie_2_ptr  = loc(edge_buf((q-1)*constLev+1, ie + 2))
      edge_buf_ie_3_ptr  = loc(edge_buf((q-1)*constLev+1, ie + 3))
      edge_buf_ie_4_ptr  = loc(edge_buf((q-1)*constLev+1, ie + 4))
      edge_buf_in_1_ptr  = loc(edge_buf((q-1)*constLev+1, in + 1))
      edge_buf_in_2_ptr  = loc(edge_buf((q-1)*constLev+1, in + 2))
      edge_buf_in_3_ptr  = loc(edge_buf((q-1)*constLev+1, in + 3))
      edge_buf_in_4_ptr  = loc(edge_buf((q-1)*constLev+1, in + 4))
      swest_buf_ptr      = loc(edge_buf((q-1)*constLev+1, getmapP(5) + 1))
      seast_buf_ptr      = loc(edge_buf((q-1)*constLev+1, getmapP(6) + 1))
      neast_buf_ptr      = loc(edge_buf((q-1)*constLev+1, getmapP(8) + 1))
      nwest_buf_ptr      = loc(edge_buf((q-1)*constLev+1, getmapP(7) + 1))

      !$ACC DATA copyin(edge_buf_is_1,edge_buf_is_2,edge_buf_is_3,edge_buf_is_4,edge_buf_iw_1,edge_buf_iw_2,edge_buf_iw_3,edge_buf_iw_4,edge_buf_ie_1,edge_buf_ie_2,edge_buf_ie_3,edge_buf_ie_4,edge_buf_in_1,edge_buf_in_2,edge_buf_in_3,edge_buf_in_4,swest_buf,seast_buf,neast_buf,nwest_buf)
      call my_edgeVunpack_bihamomic(edge_nlyr, edge_nbuf, getmapP(:), &
        qtens(:,:,:,q,iee), constLev, my_swest, my_max_corner_elem, &
        swest_buf, seast_buf, neast_buf, nwest_buf, &
        edge_buf_is_1, edge_buf_is_2, edge_buf_is_3, edge_buf_is_4, &
        edge_buf_iw_1, edge_buf_iw_2, edge_buf_iw_3, edge_buf_iw_4, &
        edge_buf_ie_1, edge_buf_ie_2, edge_buf_ie_3, edge_buf_ie_4, &
        edge_buf_in_1, edge_buf_in_2, edge_buf_in_3, edge_buf_in_4)
      !$ACC END DATA

      elem_rspheremp_ptr               = loc(my_elem(iee)%rspheremp)
      elem_spheremp_ptr                = loc(my_elem(iee)%spheremp)
      elem_Dinv_ptr                    = loc(my_elem(iee)%Dinv)
      elem_variable_hyperviscosity_ptr = loc(my_elem(iee)%variable_hyperviscosity)
      elem_tensorVisc_ptr              = loc(my_elem(iee)%tensorVisc)

      !$ACC DATA copyin(elem_rspheremp,elem_spheremp,elem_Dinv,elem_variable_hyperviscosity,elem_tensorVisc)
      do k=1,constLev
        lap_p(:,:)=elem_rspheremp(:,:)*qtens(:,:,k,q,iee)
        call my_laplace_sphere_wk(lap_p(:,:), deriv_Dvv(:,:), elem_spheremp(:,:), elem_Dinv, elem_variable_hyperviscosity, &
                                elem_tensorVisc, qtens(:,:,k,q,iee), my_rrearth, var_coef=.true. )
      enddo !k
      !$ACC END DATA

      ! end for first part variables
      !$ACC END DATA
    end do !q
  end do !ie
  !$ACC END PARALLEL LOOP

  !$ACC PARALLEL LOOP collapse(2) tile(q:3) local(qmin) copy(emin) 
  do iee = nets , nete
    do q = 1, my_qsize
      do k=1,constLev
        Qmin(:,:,k)=emin(k,q,iee)
      enddo

      getmapP_ptr        = loc(my_elem(iee)%desc%getmapP)
      is_ptr             = loc(my_elem(iee)%desc%getmapP(my_south))
      ie_ptr             = loc(my_elem(iee)%desc%getmapP(my_east))
      iw_ptr             = loc(my_elem(iee)%desc%getmapP(my_west))
      in_ptr             = loc(my_elem(iee)%desc%getmapP(my_north))
      !$ACC DATA copyin(getmapP,is,ie,iw,in)

      edge_buf_is_1_ptr  = loc(edge_buf(my_qsize*constLev+(q-1)*constLev+1, is + 1))
      edge_buf_is_2_ptr  = loc(edge_buf(my_qsize*constLev+(q-1)*constLev+1, is + 2))
      edge_buf_is_3_ptr  = loc(edge_buf(my_qsize*constLev+(q-1)*constLev+1, is + 3))
      edge_buf_is_4_ptr  = loc(edge_buf(my_qsize*constLev+(q-1)*constLev+1, is + 4))
      edge_buf_iw_1_ptr  = loc(edge_buf(my_qsize*constLev+(q-1)*constLev+1, iw + 1))
      edge_buf_iw_2_ptr  = loc(edge_buf(my_qsize*constLev+(q-1)*constLev+1, iw + 2))
      edge_buf_iw_3_ptr  = loc(edge_buf(my_qsize*constLev+(q-1)*constLev+1, iw + 3))
      edge_buf_iw_4_ptr  = loc(edge_buf(my_qsize*constLev+(q-1)*constLev+1, iw + 4))
      edge_buf_ie_1_ptr  = loc(edge_buf(my_qsize*constLev+(q-1)*constLev+1, ie + 1))
      edge_buf_ie_2_ptr  = loc(edge_buf(my_qsize*constLev+(q-1)*constLev+1, ie + 2))
      edge_buf_ie_3_ptr  = loc(edge_buf(my_qsize*constLev+(q-1)*constLev+1, ie + 3))
      edge_buf_ie_4_ptr  = loc(edge_buf(my_qsize*constLev+(q-1)*constLev+1, ie + 4))
      edge_buf_in_1_ptr  = loc(edge_buf(my_qsize*constLev+(q-1)*constLev+1, in + 1))
      edge_buf_in_2_ptr  = loc(edge_buf(my_qsize*constLev+(q-1)*constLev+1, in + 2))
      edge_buf_in_3_ptr  = loc(edge_buf(my_qsize*constLev+(q-1)*constLev+1, in + 3))
      edge_buf_in_4_ptr  = loc(edge_buf(my_qsize*constLev+(q-1)*constLev+1, in + 4))
      swest_buf_ptr      = loc(edge_buf(my_qsize*constLev+(q-1)*constLev+1, getmapP(5) + 1))
      seast_buf_ptr      = loc(edge_buf(my_qsize*constLev+(q-1)*constLev+1, getmapP(6) + 1))
      neast_buf_ptr      = loc(edge_buf(my_qsize*constLev+(q-1)*constLev+1, getmapP(8) + 1))
      nwest_buf_ptr      = loc(edge_buf(my_qsize*constLev+(q-1)*constLev+1, getmapP(7) + 1))

      !$ACC DATA copyin(edge_buf_is_1,edge_buf_is_2,edge_buf_is_3,edge_buf_is_4,edge_buf_iw_1,edge_buf_iw_2,edge_buf_iw_3,edge_buf_iw_4,edge_buf_ie_1,edge_buf_ie_2,edge_buf_ie_3,edge_buf_ie_4,edge_buf_in_1,edge_buf_in_2,edge_buf_in_3,edge_buf_in_4,swest_buf,seast_buf,neast_buf,nwest_buf)
      call my_edgeVunpack_min_bihamomic(edge_nlyr, edge_nbuf, getmapP(:), &
        Qmin(:,:,:), constLev, my_swest, my_max_corner_elem, &
        swest_buf, seast_buf, neast_buf, nwest_buf, &
        edge_buf_is_1, edge_buf_is_2, edge_buf_is_3, edge_buf_is_4, &
        edge_buf_iw_1, edge_buf_iw_2, edge_buf_iw_3, edge_buf_iw_4, &
        edge_buf_ie_1, edge_buf_ie_2, edge_buf_ie_3, edge_buf_ie_4, &
        edge_buf_in_1, edge_buf_in_2, edge_buf_in_3, edge_buf_in_4)
      !$ACC END DATA

      do k=1,constLev
        emin(k,q,iee)=min(qmin(1,1,k),qmin(1,4,k),qmin(4,1,k),qmin(4,4,k))
        emin(k,q,iee)=max(emin(k,q,iee),0d0)
      enddo !k

      ! end for first part variables
      !$ACC END DATA
    end do !q
  end do !ie
  !$ACC END PARALLEL LOOP


  !$ACC PARALLEL LOOP collapse(2) tile(q:3) local(qmax) copy(emax) 
  do iee = nets , nete
    do q = 1, my_qsize
      do k=1,constLev
        Qmax(:,:,k)=emax(k,q,iee)
      enddo

      getmapP_ptr        = loc(my_elem(iee)%desc%getmapP)
      is_ptr             = loc(my_elem(iee)%desc%getmapP(my_south))
      ie_ptr             = loc(my_elem(iee)%desc%getmapP(my_east))
      iw_ptr             = loc(my_elem(iee)%desc%getmapP(my_west))
      in_ptr             = loc(my_elem(iee)%desc%getmapP(my_north))
      !$ACC DATA copyin(getmapP,is,ie,iw,in)

      edge_buf_is_1_ptr  = loc(edge_buf(2*my_qsize*constLev+(q-1)*constLev+1, is + 1))
      edge_buf_is_2_ptr  = loc(edge_buf(2*my_qsize*constLev+(q-1)*constLev+1, is + 2))
      edge_buf_is_3_ptr  = loc(edge_buf(2*my_qsize*constLev+(q-1)*constLev+1, is + 3))
      edge_buf_is_4_ptr  = loc(edge_buf(2*my_qsize*constLev+(q-1)*constLev+1, is + 4))
      edge_buf_iw_1_ptr  = loc(edge_buf(2*my_qsize*constLev+(q-1)*constLev+1, iw + 1))
      edge_buf_iw_2_ptr  = loc(edge_buf(2*my_qsize*constLev+(q-1)*constLev+1, iw + 2))
      edge_buf_iw_3_ptr  = loc(edge_buf(2*my_qsize*constLev+(q-1)*constLev+1, iw + 3))
      edge_buf_iw_4_ptr  = loc(edge_buf(2*my_qsize*constLev+(q-1)*constLev+1, iw + 4))
      edge_buf_ie_1_ptr  = loc(edge_buf(2*my_qsize*constLev+(q-1)*constLev+1, ie + 1))
      edge_buf_ie_2_ptr  = loc(edge_buf(2*my_qsize*constLev+(q-1)*constLev+1, ie + 2))
      edge_buf_ie_3_ptr  = loc(edge_buf(2*my_qsize*constLev+(q-1)*constLev+1, ie + 3))
      edge_buf_ie_4_ptr  = loc(edge_buf(2*my_qsize*constLev+(q-1)*constLev+1, ie + 4))
      edge_buf_in_1_ptr  = loc(edge_buf(2*my_qsize*constLev+(q-1)*constLev+1, in + 1))
      edge_buf_in_2_ptr  = loc(edge_buf(2*my_qsize*constLev+(q-1)*constLev+1, in + 2))
      edge_buf_in_3_ptr  = loc(edge_buf(2*my_qsize*constLev+(q-1)*constLev+1, in + 3))
      edge_buf_in_4_ptr  = loc(edge_buf(2*my_qsize*constLev+(q-1)*constLev+1, in + 4))
      swest_buf_ptr      = loc(edge_buf(2*my_qsize*constLev+(q-1)*constLev+1, getmapP(5) + 1))
      seast_buf_ptr      = loc(edge_buf(2*my_qsize*constLev+(q-1)*constLev+1, getmapP(6) + 1))
      neast_buf_ptr      = loc(edge_buf(2*my_qsize*constLev+(q-1)*constLev+1, getmapP(8) + 1))
      nwest_buf_ptr      = loc(edge_buf(2*my_qsize*constLev+(q-1)*constLev+1, getmapP(7) + 1))

      !$ACC DATA copyin(edge_buf_is_1,edge_buf_is_2,edge_buf_is_3,edge_buf_is_4,edge_buf_iw_1,edge_buf_iw_2,edge_buf_iw_3,edge_buf_iw_4,edge_buf_ie_1,edge_buf_ie_2,edge_buf_ie_3,edge_buf_ie_4,edge_buf_in_1,edge_buf_in_2,edge_buf_in_3,edge_buf_in_4,swest_buf,seast_buf,neast_buf,nwest_buf)
      call my_edgeVunpack_max_bihamomic(edge_nlyr, edge_nbuf, getmapP(:), &
        Qmax(:,:,:), constLev, my_swest, my_max_corner_elem, &
        swest_buf, seast_buf, neast_buf, nwest_buf, &
        edge_buf_is_1, edge_buf_is_2, edge_buf_is_3, edge_buf_is_4, &
        edge_buf_iw_1, edge_buf_iw_2, edge_buf_iw_3, edge_buf_iw_4, &
        edge_buf_ie_1, edge_buf_ie_2, edge_buf_ie_3, edge_buf_ie_4, &
        edge_buf_in_1, edge_buf_in_2, edge_buf_in_3, edge_buf_in_4)
      !$ACC END DATA

      do k=1,constLev
        emax(k,q,iee)=max(qmax(1,1,k),qmax(1,4,k),qmax(4,1,k),qmax(4,4,k))
      enddo !k

      ! end for first part variables
      !$ACC END DATA
    end do !q
  end do !ie
  !$ACC END PARALLEL LOOP

end subroutine

subroutine my_biharmonic_wk_scalar_minmax_before_bndry(nets, nete, my_elem, emin, emax, qtens, deriv, edgeq_buf, edge_nbuf, edge_nlyr, my_rrearth, my_qsize)
  implicit none
  !integer, parameter :: west  = 1
  !integer, parameter :: east  = 2
  !integer, parameter :: south = 3
  !integer, parameter :: north = 4

  !integer, parameter :: swest = 5
  !integer, parameter :: seast = 6
  !integer, parameter :: nwest = 7
  !integer, parameter :: neast = 8
  integer, parameter :: max_corner_elem = 1
  integer, parameter :: my_nlev = constLev 
  !integer, parameter :: my_qsize = 25
  integer, parameter :: my_np   = 4
  integer,                                                 intent(in)    :: nets, nete, my_qsize
  type (element_t), dimension(nets:nete),                  intent(inout) :: my_elem
  real (kind=8),    dimension(nlev,qsize,nets:nete) ,      intent(in)    :: emin,emax
  real (kind=8),    dimension(np,np,nlev,qsize,nets:nete), intent(inout) :: qtens
  type (derivative_t),                                     intent(in)    :: deriv

  integer,                                                 intent(in)    :: edge_nbuf, edge_nlyr
  real(kind=8),dimension(edge_nlyr, edge_nbuf),            intent(inout) :: edgeq_buf

  real (kind=8),                                           intent(in)    :: my_rrearth
  !local
  real (kind=8), dimension(np,np,nlev) :: Qmin
  real (kind=8), dimension(np,np,nlev) :: Qmax
  real (kind=8), dimension(np,np) :: lap_p
  logical :: var_coef
  integer :: ie, q, k

  real(kind=8), dimension(my_np, my_np) :: deriv_dvv
  pointer(deriv_dvv_ptr, deriv_dvv)

  real(kind=8), dimension(my_np, my_np) :: elem_spheremp
  pointer(elem_spheremp_ptr, elem_spheremp)

  real(kind=8), dimension(2, 2, my_np, my_np) :: elem_Dinv
  pointer(elem_Dinv_ptr, elem_Dinv)

  real(kind=8), dimension(my_np, my_np) :: elem_variable_hyperviscosity
  pointer(elem_variable_hyperviscosity_ptr, elem_variable_hyperviscosity)

  real(kind=8), dimension(2, 2, my_np, my_np) :: elem_tensorVisc
  pointer(elem_tensorVisc_ptr, elem_tensorVisc)


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

  logical, dimension(8) :: elem_desc_reverse
  pointer(elem_desc_reverse_ptr, elem_desc_reverse)

  integer, dimension(8) :: elem_desc_putmapP
  pointer(elem_desc_putmapP_ptr, elem_desc_putmapP)

  deriv_Dvv_ptr = loc(deriv%dvv)

  !$ACC PARALLEL LOOP collapse(2)  local(lap_p)  copy(qtens) copyin(deriv_Dvv) annotate(entire(deriv_Dvv))
  do ie=nets,nete
    do q=1,my_qsize
      elem_spheremp_ptr                = loc(my_elem(ie)%spheremp)
      elem_Dinv_ptr                    = loc(my_elem(ie)%Dinv)
      elem_variable_hyperviscosity_ptr = loc(my_elem(ie)%variable_hyperviscosity)
      elem_tensorVisc_ptr              = loc(my_elem(ie)%tensorVisc)

      !$ACC DATA COPYIN(elem_spheremp,elem_Dinv,elem_variable_hyperviscosity,elem_tensorVisc)
      do k=1,my_nlev    !  Potential loop inversion (AAM)
        lap_p(:,:)   =qtens(:,:,k,q,ie)
        call  my_laplace_sphere_wk(lap_p, deriv_Dvv, elem_spheremp, elem_Dinv, elem_variable_hyperviscosity, elem_tensorVisc, qtens(:,:,k,q,ie), my_rrearth, var_coef=.true.)
      enddo
      !$ACC END DATA

      elem_desc_reverse_ptr            = loc(my_elem(ie)%desc%reverse)
      elem_desc_putmapP_ptr            = loc(my_elem(ie)%desc%putmapP)

      !$ACC DATA COPYIN(elem_desc_putmapP, elem_desc_reverse)

      edge_buf_5_ptr    = loc(edgeq_buf(1+my_nlev*(q-1), elem_desc_putmapP(5)+1))
      edge_buf_6_ptr    = loc(edgeq_buf(1+my_nlev*(q-1), elem_desc_putmapP(6)+1))
      edge_buf_7_ptr    = loc(edgeq_buf(1+my_nlev*(q-1), elem_desc_putmapP(7)+1))
      edge_buf_8_ptr    = loc(edgeq_buf(1+my_nlev*(q-1), elem_desc_putmapP(8)+1))
      edge_buf_in_1_ptr = loc(edgeq_buf(1+my_nlev*(q-1), elem_desc_putmapP(north)+1))
      edge_buf_in_2_ptr = loc(edgeq_buf(1+my_nlev*(q-1), elem_desc_putmapP(north)+2))
      edge_buf_in_3_ptr = loc(edgeq_buf(1+my_nlev*(q-1), elem_desc_putmapP(north)+3))
      edge_buf_in_4_ptr = loc(edgeq_buf(1+my_nlev*(q-1), elem_desc_putmapP(north)+4))
      edge_buf_is_1_ptr = loc(edgeq_buf(1+my_nlev*(q-1), elem_desc_putmapP(south)+1))
      edge_buf_is_2_ptr = loc(edgeq_buf(1+my_nlev*(q-1), elem_desc_putmapP(south)+2))
      edge_buf_is_3_ptr = loc(edgeq_buf(1+my_nlev*(q-1), elem_desc_putmapP(south)+3))
      edge_buf_is_4_ptr = loc(edgeq_buf(1+my_nlev*(q-1), elem_desc_putmapP(south)+4))
      edge_buf_ie_1_ptr = loc(edgeq_buf(1+my_nlev*(q-1), elem_desc_putmapP(east)+1))
      edge_buf_ie_2_ptr = loc(edgeq_buf(1+my_nlev*(q-1), elem_desc_putmapP(east)+2))
      edge_buf_ie_3_ptr = loc(edgeq_buf(1+my_nlev*(q-1), elem_desc_putmapP(east)+3))
      edge_buf_ie_4_ptr = loc(edgeq_buf(1+my_nlev*(q-1), elem_desc_putmapP(east)+4))
      edge_buf_iw_1_ptr = loc(edgeq_buf(1+my_nlev*(q-1), elem_desc_putmapP(west)+1))
      edge_buf_iw_2_ptr = loc(edgeq_buf(1+my_nlev*(q-1), elem_desc_putmapP(west)+2))
      edge_buf_iw_3_ptr = loc(edgeq_buf(1+my_nlev*(q-1), elem_desc_putmapP(west)+3))
      edge_buf_iw_4_ptr = loc(edgeq_buf(1+my_nlev*(q-1), elem_desc_putmapP(west)+4))

      !$ACC DATA COPY(edge_buf_5,edge_buf_6,edge_buf_7,edge_buf_8) COPYOUT(edge_buf_in_1,edge_buf_in_2,edge_buf_in_3,edge_buf_in_4,edge_buf_is_1,edge_buf_is_2,edge_buf_is_3,edge_buf_is_4,edge_buf_ie_1,edge_buf_ie_2,edge_buf_ie_3,edge_buf_ie_4,edge_buf_iw_1,edge_buf_iw_2,edge_buf_iw_3,edge_buf_iw_4)
      call my_edgeVpack_acc(qtens(:,:,:,q,ie) , my_nlev , 0+my_nlev*(q-1) , elem_desc_putmapP(:), &
         elem_desc_reverse, &
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


  !$ACC PARALLEL LOOP collapse(2) tile(q:3) local(Qmin) copyin(emin) 
  do ie=nets,nete
    do q=1,my_qsize

      do k=1,my_nlev    !  Potential loop inversion (AAM)
        Qmin(:,:,k)=emin(k,q,ie)  ! need to set all values in element for
      enddo

      elem_desc_reverse_ptr            = loc(my_elem(ie)%desc%reverse)
      elem_desc_putmapP_ptr            = loc(my_elem(ie)%desc%putmapP)

      !$ACC DATA COPYIN(elem_desc_putmapP, elem_desc_reverse)

      edge_buf_5_ptr    = loc(edgeq_buf(1+my_nlev*(q-1)+my_nlev*my_qsize, elem_desc_putmapP(5)+1))
      edge_buf_6_ptr    = loc(edgeq_buf(1+my_nlev*(q-1)+my_nlev*my_qsize, elem_desc_putmapP(6)+1))
      edge_buf_7_ptr    = loc(edgeq_buf(1+my_nlev*(q-1)+my_nlev*my_qsize, elem_desc_putmapP(7)+1))
      edge_buf_8_ptr    = loc(edgeq_buf(1+my_nlev*(q-1)+my_nlev*my_qsize, elem_desc_putmapP(8)+1))
      edge_buf_in_1_ptr = loc(edgeq_buf(1+my_nlev*(q-1)+my_nlev*my_qsize, elem_desc_putmapP(north)+1))
      edge_buf_in_2_ptr = loc(edgeq_buf(1+my_nlev*(q-1)+my_nlev*my_qsize, elem_desc_putmapP(north)+2))
      edge_buf_in_3_ptr = loc(edgeq_buf(1+my_nlev*(q-1)+my_nlev*my_qsize, elem_desc_putmapP(north)+3))
      edge_buf_in_4_ptr = loc(edgeq_buf(1+my_nlev*(q-1)+my_nlev*my_qsize, elem_desc_putmapP(north)+4))
      edge_buf_is_1_ptr = loc(edgeq_buf(1+my_nlev*(q-1)+my_nlev*my_qsize, elem_desc_putmapP(south)+1))
      edge_buf_is_2_ptr = loc(edgeq_buf(1+my_nlev*(q-1)+my_nlev*my_qsize, elem_desc_putmapP(south)+2))
      edge_buf_is_3_ptr = loc(edgeq_buf(1+my_nlev*(q-1)+my_nlev*my_qsize, elem_desc_putmapP(south)+3))
      edge_buf_is_4_ptr = loc(edgeq_buf(1+my_nlev*(q-1)+my_nlev*my_qsize, elem_desc_putmapP(south)+4))
      edge_buf_ie_1_ptr = loc(edgeq_buf(1+my_nlev*(q-1)+my_nlev*my_qsize, elem_desc_putmapP(east)+1))
      edge_buf_ie_2_ptr = loc(edgeq_buf(1+my_nlev*(q-1)+my_nlev*my_qsize, elem_desc_putmapP(east)+2))
      edge_buf_ie_3_ptr = loc(edgeq_buf(1+my_nlev*(q-1)+my_nlev*my_qsize, elem_desc_putmapP(east)+3))
      edge_buf_ie_4_ptr = loc(edgeq_buf(1+my_nlev*(q-1)+my_nlev*my_qsize, elem_desc_putmapP(east)+4))
      edge_buf_iw_1_ptr = loc(edgeq_buf(1+my_nlev*(q-1)+my_nlev*my_qsize, elem_desc_putmapP(west)+1))
      edge_buf_iw_2_ptr = loc(edgeq_buf(1+my_nlev*(q-1)+my_nlev*my_qsize, elem_desc_putmapP(west)+2))
      edge_buf_iw_3_ptr = loc(edgeq_buf(1+my_nlev*(q-1)+my_nlev*my_qsize, elem_desc_putmapP(west)+3))
      edge_buf_iw_4_ptr = loc(edgeq_buf(1+my_nlev*(q-1)+my_nlev*my_qsize, elem_desc_putmapP(west)+4))

      !$ACC DATA COPY(edge_buf_5,edge_buf_6,edge_buf_7,edge_buf_8) COPYOUT(edge_buf_in_1,edge_buf_in_2,edge_buf_in_3,edge_buf_in_4,edge_buf_is_1,edge_buf_is_2,edge_buf_is_3,edge_buf_is_4,edge_buf_ie_1,edge_buf_ie_2,edge_buf_ie_3,edge_buf_ie_4,edge_buf_iw_1,edge_buf_iw_2,edge_buf_iw_3,edge_buf_iw_4)
     call my_edgeVpack_acc(Qmin(:,:,:) , my_nlev , my_nlev*qsize+my_nlev*(q-1) , elem_desc_putmapP(:), &
         elem_desc_reverse, &
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


  !$ACC PARALLEL LOOP collapse(2) tile(q:3) local(Qmax) copyin(emax) 
  do ie=nets,nete
    do q=1,my_qsize

      do k=1,my_nlev    !  Potential loop inversion (AAM)
        Qmax(:,:,k)=emax(k,q,ie)  ! edgeVpack routine below
      enddo

      elem_desc_reverse_ptr            = loc(my_elem(ie)%desc%reverse)
      elem_desc_putmapP_ptr            = loc(my_elem(ie)%desc%putmapP)

      !$ACC DATA COPYIN(elem_desc_putmapP, elem_desc_reverse)
      edge_buf_5_ptr = loc(edgeq_buf(1+my_nlev*(q-1)+2*my_nlev*my_qsize, elem_desc_putmapP(5)+1))
      edge_buf_6_ptr = loc(edgeq_buf(1+my_nlev*(q-1)+2*my_nlev*my_qsize, elem_desc_putmapP(6)+1))
      edge_buf_7_ptr = loc(edgeq_buf(1+my_nlev*(q-1)+2*my_nlev*my_qsize, elem_desc_putmapP(7)+1))
      edge_buf_8_ptr = loc(edgeq_buf(1+my_nlev*(q-1)+2*my_nlev*my_qsize, elem_desc_putmapP(8)+1))

      edge_buf_in_1_ptr = loc(edgeq_buf(1+my_nlev*(q-1)+2*my_nlev*my_qsize, elem_desc_putmapP(north)+1))
      edge_buf_in_2_ptr = loc(edgeq_buf(1+my_nlev*(q-1)+2*my_nlev*my_qsize, elem_desc_putmapP(north)+2))
      edge_buf_in_3_ptr = loc(edgeq_buf(1+my_nlev*(q-1)+2*my_nlev*my_qsize, elem_desc_putmapP(north)+3))
      edge_buf_in_4_ptr = loc(edgeq_buf(1+my_nlev*(q-1)+2*my_nlev*my_qsize, elem_desc_putmapP(north)+4))

      edge_buf_is_1_ptr = loc(edgeq_buf(1+my_nlev*(q-1)+2*my_nlev*my_qsize, elem_desc_putmapP(south)+1))
      edge_buf_is_2_ptr = loc(edgeq_buf(1+my_nlev*(q-1)+2*my_nlev*my_qsize, elem_desc_putmapP(south)+2))
      edge_buf_is_3_ptr = loc(edgeq_buf(1+my_nlev*(q-1)+2*my_nlev*my_qsize, elem_desc_putmapP(south)+3))
      edge_buf_is_4_ptr = loc(edgeq_buf(1+my_nlev*(q-1)+2*my_nlev*my_qsize, elem_desc_putmapP(south)+4))

      edge_buf_ie_1_ptr = loc(edgeq_buf(1+my_nlev*(q-1)+2*my_nlev*my_qsize, elem_desc_putmapP(east)+1))
      edge_buf_ie_2_ptr = loc(edgeq_buf(1+my_nlev*(q-1)+2*my_nlev*my_qsize, elem_desc_putmapP(east)+2))
      edge_buf_ie_3_ptr = loc(edgeq_buf(1+my_nlev*(q-1)+2*my_nlev*my_qsize, elem_desc_putmapP(east)+3))
      edge_buf_ie_4_ptr = loc(edgeq_buf(1+my_nlev*(q-1)+2*my_nlev*my_qsize, elem_desc_putmapP(east)+4))

      edge_buf_iw_1_ptr = loc(edgeq_buf(1+my_nlev*(q-1)+2*my_nlev*my_qsize, elem_desc_putmapP(west)+1))
      edge_buf_iw_2_ptr = loc(edgeq_buf(1+my_nlev*(q-1)+2*my_nlev*my_qsize, elem_desc_putmapP(west)+2))
      edge_buf_iw_3_ptr = loc(edgeq_buf(1+my_nlev*(q-1)+2*my_nlev*my_qsize, elem_desc_putmapP(west)+3))
      edge_buf_iw_4_ptr = loc(edgeq_buf(1+my_nlev*(q-1)+2*my_nlev*my_qsize, elem_desc_putmapP(west)+4))

      !$ACC DATA COPY(edge_buf_5,edge_buf_6,edge_buf_7,edge_buf_8) COPYOUT(edge_buf_in_1,edge_buf_in_2,edge_buf_in_3,edge_buf_in_4,edge_buf_is_1,edge_buf_is_2,edge_buf_is_3,edge_buf_is_4,edge_buf_ie_1,edge_buf_ie_2,edge_buf_ie_3,edge_buf_ie_4,edge_buf_iw_1,edge_buf_iw_2,edge_buf_iw_3,edge_buf_iw_4)
      call my_edgeVpack_acc(Qmax(:,:,:) , my_nlev , 2*my_nlev*qsize+my_nlev*(q-1) , elem_desc_putmapP(:), &
         elem_desc_reverse, &
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
end subroutine

subroutine biharmonic_wk_scalar_minmax(elem,qtens,deriv,edgeq,hybrid,nets,nete,emin,emax)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! compute weak biharmonic operator
!    input:  qtens = Q
!    output: qtens = weak biharmonic of Q and Q element min/max
!
!    note: emin/emax must be initialized with Q element min/max.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use physical_constants, only : rrearth
type (hybrid_t)      , intent(in) :: hybrid
type (element_t)     , intent(inout), target :: elem(:)
integer :: nets,nete
real (kind=real_kind), dimension(np,np,nlev,qsize,nets:nete) :: qtens
type (EdgeBuffer_t)  , intent(inout) :: edgeq
type (derivative_t)  , intent(in) :: deriv
real (kind=real_kind), intent(out), dimension(nlev,qsize,nets:nete) :: emin,emax

! local
integer :: k,kptr,i,j,ie,ic,q
real (kind=real_kind), dimension(np,np) :: lap_p
real (kind=real_kind) :: Qmin(np,np,nlev,qsize)
real (kind=real_kind) :: Qmax(np,np,nlev,qsize)
integer(kind=8) :: count_start, count_stop, count_rate, count_max

   !do ie=nets,nete
      !do q=1,qsize
      !do k=1,nlev    !  Potential loop inversion (AAM)
         !Qmin(:,:,k,q)=emin(k,q,ie)  ! need to set all values in element for
         !Qmax(:,:,k,q)=emax(k,q,ie)  ! edgeVpack routine below
         !lap_p(:,:) = qtens(:,:,k,q,ie)
!! Original use of qtens on left and right hand sides caused OpenMP errors (AAM)
         !qtens(:,:,k,q,ie)=laplace_sphere_wk(lap_p,deriv,elem(ie),var_coef=.true.)
      !enddo
      !enddo
      !call edgeVpack(edgeq, qtens(:,:,:,:,ie),qsize*nlev,0,elem(ie)%desc)
      !call edgeVpack(edgeq,Qmin,nlev*qsize,nlev*qsize,elem(ie)%desc)
      !call edgeVpack(edgeq,Qmax,nlev*qsize,2*nlev*qsize,elem(ie)%desc)
   !enddo

   !call system_clock(count_start, count_rate, count_max)
   call my_biharmonic_wk_scalar_minmax_before_bndry(nets, nete, elem, emin, emax, qtens, deriv, edgeq%buf, edgeq%nbuf, edgeq%nlyr, rrearth, qsize)
   !call system_clock(count_stop, count_rate, count_max)
   !write(*,*), 'acc my_biharmonic_wk_scalar_minmax_before_bndry count = ',   count_stop - count_start

   call bndry_exchangeV(hybrid,edgeq)

   !call system_clock(count_start, count_rate, count_max)
   call my_biharmonic_wk_scalar_minmax_after_bndry(nets, nete, edgeq%nlyr, edgeq%nbuf, &
                        edgeq%buf, south, east, north, west, qsize,&
                        swest, max_corner_elem, elem, emin, emax, qtens, deriv,rrearth)
   !call system_clock(count_stop, count_rate, count_max)
   !write(*,*), 'acc my_biharmonic_wk_scalar_minmax_after_bndry count = ',   count_stop - count_start

end subroutine

#endif





subroutine make_C0_2d(zeta,elem,hybrid,nets,nete)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! apply DSS (aka assembly procedure) to zeta.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

type (hybrid_t)      , intent(in) :: hybrid
type (element_t)     , intent(in), target :: elem(:)
integer :: nets,nete
real (kind=real_kind), dimension(np,np,nets:nete) :: zeta

! local
integer :: k,i,j,ie,ic,kptr


call initEdgeBuffer(edge1,1)

do ie=nets,nete
   zeta(:,:,ie)=zeta(:,:,ie)*elem(ie)%spheremp(:,:)
   kptr=0
   call edgeVpack(edge1, zeta(1,1,ie),1,kptr,elem(ie)%desc)
enddo
call bndry_exchangeV(hybrid,edge1)
do ie=nets,nete
   kptr=0
   call edgeVunpack(edge1, zeta(1,1,ie),1,kptr,elem(ie)%desc)
   zeta(:,:,ie)=zeta(:,:,ie)*elem(ie)%rspheremp(:,:)
enddo

call FreeEdgeBuffer(edge1)
end subroutine


subroutine make_C0(zeta,elem,hybrid,nets,nete)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! apply DSS (aka assembly procedure) to zeta.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

type (hybrid_t)      , intent(in) :: hybrid
type (element_t)     , intent(in), target :: elem(:)
integer :: nets,nete
real (kind=real_kind), dimension(np,np,nlev,nets:nete) :: zeta

! local
integer :: k,i,j,ie,ic,kptr


call initEdgeBuffer(edge1,nlev)

do ie=nets,nete
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k)
#endif
   do k=1,nlev
      zeta(:,:,k,ie)=zeta(:,:,k,ie)*elem(ie)%spheremp(:,:)
   enddo
   kptr=0
   call edgeVpack(edge1, zeta(1,1,1,ie),nlev,kptr,elem(ie)%desc)
enddo
call bndry_exchangeV(hybrid,edge1)
do ie=nets,nete
   kptr=0
   call edgeVunpack(edge1, zeta(1,1,1,ie),nlev,kptr,elem(ie)%desc)
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k)
#endif
   do k=1,nlev
      zeta(:,:,k,ie)=zeta(:,:,k,ie)*elem(ie)%rspheremp(:,:)
   enddo
enddo
#ifdef DEBUGOMP
#if (! defined ELEMENT_OPENMP)
!$OMP BARRIER
#endif
#endif

call FreeEdgeBuffer(edge1)
end subroutine


subroutine make_C0_vector(v,elem,hybrid,nets,nete)
type (hybrid_t)      , intent(in) :: hybrid
type (element_t)     , intent(in), target :: elem(:)
integer :: nets,nete
real (kind=real_kind), dimension(np,np,2,nlev,nets:nete) :: v

! local
integer :: k,i,j,ie,ic,kptr
type (EdgeBuffer_t)          :: edge2


call initEdgeBuffer(edge2,2*nlev)

do ie=nets,nete
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k)
#endif
   do k=1,nlev
      v(:,:,1,k,ie)=v(:,:,1,k,ie)*elem(ie)%spheremp(:,:)
      v(:,:,2,k,ie)=v(:,:,2,k,ie)*elem(ie)%spheremp(:,:)
   enddo
   kptr=0
   call edgeVpack(edge2, v(1,1,1,1,ie),2*nlev,kptr,elem(ie)%desc)
enddo
call bndry_exchangeV(hybrid,edge2)
do ie=nets,nete
   kptr=0
   call edgeVunpack(edge2, v(1,1,1,1,ie),2*nlev,kptr,elem(ie)%desc)
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k)
#endif
   do k=1,nlev
      v(:,:,1,k,ie)=v(:,:,1,k,ie)*elem(ie)%rspheremp(:,:)
      v(:,:,2,k,ie)=v(:,:,2,k,ie)*elem(ie)%rspheremp(:,:)
   enddo
enddo
#ifdef DEBUGOMP
#if (! defined ELEMENT_OPENMP)
!$OMP BARRIER
#endif
#endif

call FreeEdgeBuffer(edge2)
end subroutine





subroutine compute_zeta_C0_2d_sphere(zeta,elem,hybrid,nets,nete,nt,k)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! compute C0 vorticity.  That is, solve:
!     < PHI, zeta > = <PHI, curl(elem%state%v >
!
!    input:  v (stored in elem()%, in lat-lon coordinates)
!    output: zeta(:,:,:,:)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

type (hybrid_t)      , intent(in) :: hybrid
type (element_t)     , intent(in), target :: elem(:)
integer :: nt,nets,nete,k
real (kind=real_kind), dimension(np,np,nets:nete) :: zeta

! local
integer :: i,j,ie,ic
type (derivative_t)          :: deriv

call derivinit(deriv)

do ie=nets,nete
   !    zeta(:,:,k,ie)=elem(ie)%state%zeta(:,:,k)
   zeta(:,:,ie)=vorticity_sphere(elem(ie)%state%v(:,:,:,k,nt),deriv,elem(ie))
enddo

call make_C0_2d(zeta,elem,hybrid,nets,nete)

end subroutine

subroutine compute_zeta_C0_2d_contra(zeta,elem,hybrid,nets,nete,nt)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! compute C0 vorticity.  That is, solve:
!     < PHI, zeta > = <PHI, curl(elem%state%v >
!
!    input:  v (stored in elem()%, in contra-variant coordinates)
!    output: zeta(:,:,:,:)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

type (hybrid_t)      , intent(in) :: hybrid
type (element_t)     , intent(in), target :: elem(:)
integer :: nt,nets,nete
real (kind=real_kind), dimension(np,np,nlev,nets:nete) :: zeta
real (kind=real_kind), dimension(np,np,2) :: ulatlon
real (kind=real_kind), dimension(np,np) :: v1,v2

! local
integer :: k,ie
type (derivative_t)          :: deriv

call derivinit(deriv)

do k=1,nlev
do ie=nets,nete
    v1 = elem(ie)%state%v(:,:,1,k,nt)
    v2 = elem(ie)%state%v(:,:,2,k,nt)
    ulatlon(:,:,1) = elem(ie)%D(1,1,:,:)*v1 + elem(ie)%D(1,2,:,:)*v2
    ulatlon(:,:,2) = elem(ie)%D(2,1,:,:)*v1 + elem(ie)%D(2,2,:,:)*v2
   !    zeta(:,:,k,ie)=elem(ie)%state%zeta(:,:,k)
   zeta(:,:,k,ie)=vorticity_sphere(ulatlon,deriv,elem(ie))
enddo
enddo

call make_C0(zeta,elem,hybrid,nets,nete)

end subroutine


subroutine compute_div_C0_2d_sphere(zeta,elem,hybrid,nets,nete,nt,k)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! compute C0 divergence. That is, solve:
!     < PHI, zeta > = <PHI, div(elem%state%v >
!
!    input:  v (stored in elem()%, in lat-lon coordinates)
!    output: zeta(:,:,:,:)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

type (hybrid_t)      , intent(in) :: hybrid
type (element_t)     , intent(in), target :: elem(:)
integer :: nt,nets,nete,k
real (kind=real_kind), dimension(np,np,nets:nete) :: zeta

! local
integer :: i,j,ie,ic
type (derivative_t)          :: deriv

call derivinit(deriv)

do ie=nets,nete
   !    zeta(:,:,k,ie)=elem(ie)%state%zeta(:,:,k)
   zeta(:,:,ie)=divergence_sphere(elem(ie)%state%v(:,:,:,k,nt),deriv,elem(ie))
enddo

call make_C0_2d(zeta,elem,hybrid,nets,nete)

end subroutine

subroutine compute_div_C0_2d_contra(zeta,elem,hybrid,nets,nete,nt)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! compute C0 divergence. That is, solve:
!     < PHI, zeta > = <PHI, div(elem%state%v >
!
!    input:  v (stored in elem()%, in contra-variant coordinates)
!    output: zeta(:,:,:,:)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

type (hybrid_t)      , intent(in) :: hybrid
type (element_t)     , intent(in), target :: elem(:)
integer :: nt,nets,nete
real (kind=real_kind), dimension(np,np,nlev,nets:nete) :: zeta
real (kind=real_kind), dimension(np,np,2) :: ulatlon
real (kind=real_kind), dimension(np,np) :: v1,v2

! local
integer :: k,ie
type (derivative_t)          :: deriv

call derivinit(deriv)

do k=1,nlev
do ie=nets,nete
    v1 = elem(ie)%state%v(:,:,1,k,nt)
    v2 = elem(ie)%state%v(:,:,2,k,nt)
    ulatlon(:,:,1) = elem(ie)%D(1,1,:,:)*v1 + elem(ie)%D(1,2,:,:)*v2
    ulatlon(:,:,2) = elem(ie)%D(2,1,:,:)*v1 + elem(ie)%D(2,2,:,:)*v2
   !    zeta(:,:,k,ie)=elem(ie)%state%zeta(:,:,k)
   zeta(:,:,k,ie)=divergence_sphere(ulatlon,deriv,elem(ie))
enddo
enddo

call make_C0(zeta,elem,hybrid,nets,nete)

end subroutine

subroutine compute_zeta_C0(zeta,elem,hybrid,nets,nete,nt)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! compute C0 vorticity.  That is, solve:
!     < PHI, zeta > = <PHI, curl(elem%state%v >
!
!    input:  v (stored in elem()%, in lat-lon coordinates)
!    output: zeta(:,:,:,:)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

type (hybrid_t)      , intent(in) :: hybrid
type (element_t)     , intent(in), target :: elem(:)
integer :: nt,nets,nete
real (kind=real_kind), dimension(np,np,nlev,nets:nete) :: zeta

! local
integer :: k,i,j,ie,ic
type (derivative_t)          :: deriv

call derivinit(deriv)

do ie=nets,nete
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k)
#endif
do k=1,nlev
   !    zeta(:,:,k,ie)=elem(ie)%state%zeta(:,:,k)
   zeta(:,:,k,ie)=vorticity_sphere(elem(ie)%state%v(:,:,:,k,nt),deriv,elem(ie))
enddo
enddo

call make_C0(zeta,elem,hybrid,nets,nete)

end subroutine


subroutine compute_div_C0(zeta,elem,hybrid,nets,nete,nt)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! compute C0 divergence. That is, solve:
!     < PHI, zeta > = <PHI, div(elem%state%v >
!
!    input:  v (stored in elem()%, in lat-lon coordinates)
!    output: zeta(:,:,:,:)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

type (hybrid_t)      , intent(in) :: hybrid
type (element_t)     , intent(in), target :: elem(:)
integer :: nt,nets,nete
real (kind=real_kind), dimension(np,np,nlev,nets:nete) :: zeta

! local
integer :: k,i,j,ie,ic
type (derivative_t)          :: deriv

call derivinit(deriv)

do ie=nets,nete
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k)
#endif
do k=1,nlev
   !    zeta(:,:,k,ie)=elem(ie)%state%zeta(:,:,k)
   zeta(:,:,k,ie)=divergence_sphere(elem(ie)%state%v(:,:,:,k,nt),deriv,elem(ie))
enddo
enddo

call make_C0(zeta,elem,hybrid,nets,nete)

end subroutine








#ifdef _PRIM


subroutine my_neighbor_minmax_after_bndry(nets, nete, edge_nlyr, edge_nbuf, &
                        edge_buf, my_south, my_east, my_north, my_west, my_qsize, &
                        my_swest, my_max_corner_elem, my_elem, &
                        emin, emax, my_rrearth)
  implicit none

  integer, intent(in) ::  edge_nlyr, edge_nbuf, my_south, my_east, my_north, my_west, my_qsize, my_swest, my_max_corner_elem
  integer, intent(in) :: nets, nete
  real(kind=8), dimension(edge_nlyr, edge_nbuf), intent(in) :: edge_buf
  type (element_t), intent(in) :: my_elem(nets:nete)
  real (kind=8), intent(out), dimension(constLev,my_qsize,nets:nete) :: emin,emax
  real(kind=8), intent(in) :: my_rrearth

  !local
  real (kind=8), dimension(4,4,constLev) :: Qmin
  real (kind=8), dimension(4,4,constLev) :: Qmax
  integer :: iee, q, k

  ! pointers
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

  !$ACC PARALLEL LOOP collapse(2) tile(q:3) local(qmin) copy(emin)
  do iee = nets , nete
    do q = 1, my_qsize

      do k=1,constLev
        Qmin(:,:,k)=emin(k,q,iee) ! restore element data.  we could avoid
      enddo

      getmapP_ptr        = loc(my_elem(iee)%desc%getmapP)
      is_ptr             = loc(my_elem(iee)%desc%getmapP(my_south))
      ie_ptr             = loc(my_elem(iee)%desc%getmapP(my_east))
      iw_ptr             = loc(my_elem(iee)%desc%getmapP(my_west))
      in_ptr             = loc(my_elem(iee)%desc%getmapP(my_north))
      !$ACC DATA copyin(getmapP,is,ie,iw,in)

      edge_buf_is_1_ptr  = loc(edge_buf((q-1)*constLev+1, is + 1))
      edge_buf_is_2_ptr  = loc(edge_buf((q-1)*constLev+1, is + 2))
      edge_buf_is_3_ptr  = loc(edge_buf((q-1)*constLev+1, is + 3))
      edge_buf_is_4_ptr  = loc(edge_buf((q-1)*constLev+1, is + 4))
      edge_buf_iw_1_ptr  = loc(edge_buf((q-1)*constLev+1, iw + 1))
      edge_buf_iw_2_ptr  = loc(edge_buf((q-1)*constLev+1, iw + 2))
      edge_buf_iw_3_ptr  = loc(edge_buf((q-1)*constLev+1, iw + 3))
      edge_buf_iw_4_ptr  = loc(edge_buf((q-1)*constLev+1, iw + 4))
      edge_buf_ie_1_ptr  = loc(edge_buf((q-1)*constLev+1, ie + 1))
      edge_buf_ie_2_ptr  = loc(edge_buf((q-1)*constLev+1, ie + 2))
      edge_buf_ie_3_ptr  = loc(edge_buf((q-1)*constLev+1, ie + 3))
      edge_buf_ie_4_ptr  = loc(edge_buf((q-1)*constLev+1, ie + 4))
      edge_buf_in_1_ptr  = loc(edge_buf((q-1)*constLev+1, in + 1))
      edge_buf_in_2_ptr  = loc(edge_buf((q-1)*constLev+1, in + 2))
      edge_buf_in_3_ptr  = loc(edge_buf((q-1)*constLev+1, in + 3))
      edge_buf_in_4_ptr  = loc(edge_buf((q-1)*constLev+1, in + 4))
      swest_buf_ptr      = loc(edge_buf((q-1)*constLev+1, getmapP(5) + 1))
      seast_buf_ptr      = loc(edge_buf((q-1)*constLev+1, getmapP(6) + 1))
      neast_buf_ptr      = loc(edge_buf((q-1)*constLev+1, getmapP(8) + 1))
      nwest_buf_ptr      = loc(edge_buf((q-1)*constLev+1, getmapP(7) + 1))

      !$ACC DATA copyin(edge_buf_is_1,edge_buf_is_2,edge_buf_is_3,edge_buf_is_4,edge_buf_iw_1,edge_buf_iw_2,edge_buf_iw_3,edge_buf_iw_4,edge_buf_ie_1,edge_buf_ie_2,edge_buf_ie_3,edge_buf_ie_4,edge_buf_in_1,edge_buf_in_2,edge_buf_in_3,edge_buf_in_4,swest_buf,seast_buf,neast_buf,nwest_buf)
      call my_edgeVunpack_min_bihamomic(edge_nlyr, edge_nbuf, getmapP(:), &
        Qmin(:,:,:), constLev, my_swest, my_max_corner_elem, &
        swest_buf, seast_buf, neast_buf, nwest_buf, &
        edge_buf_is_1, edge_buf_is_2, edge_buf_is_3, edge_buf_is_4, &
        edge_buf_iw_1, edge_buf_iw_2, edge_buf_iw_3, edge_buf_iw_4, &
        edge_buf_ie_1, edge_buf_ie_2, edge_buf_ie_3, edge_buf_ie_4, &
        edge_buf_in_1, edge_buf_in_2, edge_buf_in_3, edge_buf_in_4)
      !$ACC END DATA
        do k=1,nlev
           emin(k,q,iee)=min(qmin(1,1,k),qmin(1,np,k),qmin(np,1,k),qmin(np,np,k))
           emin(k,q,iee)=max(emin(k,q,iee),0d0)
        enddo

      ! end for first part variables
      !$ACC END DATA
    end do !q
  end do !ie
  !$ACC END PARALLEL LOOP

  !$ACC PARALLEL LOOP collapse(2) tile(q:3) local(qmax) copy(emax)
  do iee = nets , nete
    do q = 1, my_qsize

      do k=1,constLev
        Qmax(:,:,k)=emax(k,q,iee) ! this by adding a "iee" index to Qmin/max
      enddo

      getmapP_ptr        = loc(my_elem(iee)%desc%getmapP)
      is_ptr             = loc(my_elem(iee)%desc%getmapP(my_south))
      ie_ptr             = loc(my_elem(iee)%desc%getmapP(my_east))
      iw_ptr             = loc(my_elem(iee)%desc%getmapP(my_west))
      in_ptr             = loc(my_elem(iee)%desc%getmapP(my_north))
      !$ACC DATA copyin(getmapP,is,ie,iw,in)
      edge_buf_is_1_ptr  = loc(edge_buf(my_qsize*constLev+(q-1)*constLev+1, is + 1))
      edge_buf_is_2_ptr  = loc(edge_buf(my_qsize*constLev+(q-1)*constLev+1, is + 2))
      edge_buf_is_3_ptr  = loc(edge_buf(my_qsize*constLev+(q-1)*constLev+1, is + 3))
      edge_buf_is_4_ptr  = loc(edge_buf(my_qsize*constLev+(q-1)*constLev+1, is + 4))
      edge_buf_iw_1_ptr  = loc(edge_buf(my_qsize*constLev+(q-1)*constLev+1, iw + 1))
      edge_buf_iw_2_ptr  = loc(edge_buf(my_qsize*constLev+(q-1)*constLev+1, iw + 2))
      edge_buf_iw_3_ptr  = loc(edge_buf(my_qsize*constLev+(q-1)*constLev+1, iw + 3))
      edge_buf_iw_4_ptr  = loc(edge_buf(my_qsize*constLev+(q-1)*constLev+1, iw + 4))
      edge_buf_ie_1_ptr  = loc(edge_buf(my_qsize*constLev+(q-1)*constLev+1, ie + 1))
      edge_buf_ie_2_ptr  = loc(edge_buf(my_qsize*constLev+(q-1)*constLev+1, ie + 2))
      edge_buf_ie_3_ptr  = loc(edge_buf(my_qsize*constLev+(q-1)*constLev+1, ie + 3))
      edge_buf_ie_4_ptr  = loc(edge_buf(my_qsize*constLev+(q-1)*constLev+1, ie + 4))
      edge_buf_in_1_ptr  = loc(edge_buf(my_qsize*constLev+(q-1)*constLev+1, in + 1))
      edge_buf_in_2_ptr  = loc(edge_buf(my_qsize*constLev+(q-1)*constLev+1, in + 2))
      edge_buf_in_3_ptr  = loc(edge_buf(my_qsize*constLev+(q-1)*constLev+1, in + 3))
      edge_buf_in_4_ptr  = loc(edge_buf(my_qsize*constLev+(q-1)*constLev+1, in + 4))
      swest_buf_ptr      = loc(edge_buf(my_qsize*constLev+(q-1)*constLev+1, getmapP(5) + 1))
      seast_buf_ptr      = loc(edge_buf(my_qsize*constLev+(q-1)*constLev+1, getmapP(6) + 1))
      neast_buf_ptr      = loc(edge_buf(my_qsize*constLev+(q-1)*constLev+1, getmapP(8) + 1))
      nwest_buf_ptr      = loc(edge_buf(my_qsize*constLev+(q-1)*constLev+1, getmapP(7) + 1))

      !$ACC DATA copyin(edge_buf_is_1,edge_buf_is_2,edge_buf_is_3,edge_buf_is_4,edge_buf_iw_1,edge_buf_iw_2,edge_buf_iw_3,edge_buf_iw_4,edge_buf_ie_1,edge_buf_ie_2,edge_buf_ie_3,edge_buf_ie_4,edge_buf_in_1,edge_buf_in_2,edge_buf_in_3,edge_buf_in_4,swest_buf,seast_buf,neast_buf,nwest_buf)
      call my_edgeVunpack_max_bihamomic(edge_nlyr, edge_nbuf, getmapP(:), &
        Qmax(:,:,:), constLev, my_swest, my_max_corner_elem, &
        swest_buf, seast_buf, neast_buf, nwest_buf, &
        edge_buf_is_1, edge_buf_is_2, edge_buf_is_3, edge_buf_is_4, &
        edge_buf_iw_1, edge_buf_iw_2, edge_buf_iw_3, edge_buf_iw_4, &
        edge_buf_ie_1, edge_buf_ie_2, edge_buf_ie_3, edge_buf_ie_4, &
        edge_buf_in_1, edge_buf_in_2, edge_buf_in_3, edge_buf_in_4)
      !$ACC END DATA

        do k=1,nlev
           emax(k,q,iee)=max(qmax(1,1,k),qmax(1,np,k),qmax(np,1,k),qmax(np,np,k))
        enddo

      ! end for first part variables
      !$ACC END DATA
    end do !q
  end do !ie
  !$ACC END PARALLEL LOOP

end subroutine

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

subroutine my_neighbor_minmax_acc_before(nets,nete,min_neigh,max_neigh,my_elem,edgeMinMax,edgeMinMax_nlyr,edgeMinMax_nbuf, edgeMinMax_buf, my_qsize)
  implicit none
  integer, parameter :: my_nlev  = constLev
  !integer, parameter :: my_qsize = 25
  integer, parameter :: my_np    = 4
  integer, intent(in) :: nets, nete, my_qsize
  real(kind=8), dimension(my_nlev,my_qsize,nets:nete), intent(in) :: min_neigh, max_neigh
  type (element_t), dimension(nets:nete), intent(in) :: my_elem
  type (EdgeBuffer_t)  , intent(inout) :: edgeMinMax
  integer, intent(in) :: edgeMinMax_nlyr, edgeMinMax_nbuf
  real(kind=8), dimension(edgeMinMax_nlyr, edgeMinMax_nbuf), intent(inout) :: edgeMinMax_buf

  !-----local-----
  integer :: ie, q, k
  real(kind=8), dimension(my_np,my_np,my_nlev) :: Qmin, Qmax

  logical, dimension(8) :: elem_desc_reverse
  pointer(elem_desc_reverse_ptr, elem_desc_reverse)

  integer, dimension(8) :: elem_desc_putmapP
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
  !$ACC PARALLEL LOOP collapse(2) tile(ie:1,q:3) local(Qmin) copyin(min_neigh,max_neigh)
  do ie=nets,nete
    do q=1,my_qsize
      do k=1,my_nlev
        Qmin(:,:,k)=min_neigh(k,q,ie)
      enddo

      elem_desc_reverse_ptr            = loc(my_elem(ie)%desc%reverse)
      elem_desc_putmapP_ptr            = loc(my_elem(ie)%desc%putmapP)

      !$ACC DATA COPYIN(elem_desc_putmapP, elem_desc_reverse)
      edge_buf_5_ptr    = loc(edgeMinMax_buf(1+my_nlev*(q-1), elem_desc_putmapP(5)+1))
      edge_buf_6_ptr    = loc(edgeMinMax_buf(1+my_nlev*(q-1), elem_desc_putmapP(6)+1))
      edge_buf_7_ptr    = loc(edgeMinMax_buf(1+my_nlev*(q-1), elem_desc_putmapP(7)+1))
      edge_buf_8_ptr    = loc(edgeMinMax_buf(1+my_nlev*(q-1), elem_desc_putmapP(8)+1))
      edge_buf_in_1_ptr = loc(edgeMinMax_buf(1+my_nlev*(q-1), elem_desc_putmapP(north)+1))
      edge_buf_in_2_ptr = loc(edgeMinMax_buf(1+my_nlev*(q-1), elem_desc_putmapP(north)+2))
      edge_buf_in_3_ptr = loc(edgeMinMax_buf(1+my_nlev*(q-1), elem_desc_putmapP(north)+3))
      edge_buf_in_4_ptr = loc(edgeMinMax_buf(1+my_nlev*(q-1), elem_desc_putmapP(north)+4))
      edge_buf_is_1_ptr = loc(edgeMinMax_buf(1+my_nlev*(q-1), elem_desc_putmapP(south)+1))
      edge_buf_is_2_ptr = loc(edgeMinMax_buf(1+my_nlev*(q-1), elem_desc_putmapP(south)+2))
      edge_buf_is_3_ptr = loc(edgeMinMax_buf(1+my_nlev*(q-1), elem_desc_putmapP(south)+3))
      edge_buf_is_4_ptr = loc(edgeMinMax_buf(1+my_nlev*(q-1), elem_desc_putmapP(south)+4))
      edge_buf_ie_1_ptr = loc(edgeMinMax_buf(1+my_nlev*(q-1), elem_desc_putmapP(east)+1))
      edge_buf_ie_2_ptr = loc(edgeMinMax_buf(1+my_nlev*(q-1), elem_desc_putmapP(east)+2))
      edge_buf_ie_3_ptr = loc(edgeMinMax_buf(1+my_nlev*(q-1), elem_desc_putmapP(east)+3))
      edge_buf_ie_4_ptr = loc(edgeMinMax_buf(1+my_nlev*(q-1), elem_desc_putmapP(east)+4))
      edge_buf_iw_1_ptr = loc(edgeMinMax_buf(1+my_nlev*(q-1), elem_desc_putmapP(west)+1))
      edge_buf_iw_2_ptr = loc(edgeMinMax_buf(1+my_nlev*(q-1), elem_desc_putmapP(west)+2))
      edge_buf_iw_3_ptr = loc(edgeMinMax_buf(1+my_nlev*(q-1), elem_desc_putmapP(west)+3))
      edge_buf_iw_4_ptr = loc(edgeMinMax_buf(1+my_nlev*(q-1), elem_desc_putmapP(west)+4))

      !$ACC DATA COPY(edge_buf_5,edge_buf_6,edge_buf_7,edge_buf_8) COPYOUT(edge_buf_in_1,edge_buf_in_2,edge_buf_in_3,edge_buf_in_4,edge_buf_is_1,edge_buf_is_2,edge_buf_is_3,edge_buf_is_4,edge_buf_ie_1,edge_buf_ie_2,edge_buf_ie_3,edge_buf_ie_4,edge_buf_iw_1,edge_buf_iw_2,edge_buf_iw_3,edge_buf_iw_4)
      call my_edgeVpack_acc(Qmin(:,:,:) , my_nlev , 0+my_nlev*(q-1) , elem_desc_putmapP(:), &
         elem_desc_reverse, &
         edge_buf_5, edge_buf_6, edge_buf_7, edge_buf_8, &
         edge_buf_in_1, edge_buf_in_2, edge_buf_in_3, edge_buf_in_4, &
         edge_buf_is_1, edge_buf_is_2, edge_buf_is_3, edge_buf_is_4, &
         edge_buf_ie_1, edge_buf_ie_2, edge_buf_ie_3, edge_buf_ie_4, &
         edge_buf_iw_1, edge_buf_iw_2, edge_buf_iw_3, edge_buf_iw_4)
      !$ACC END DATA

      !$ACC END DATA
      end do
    enddo
    !$ACC END PARALLEL LOOP

   !$ACC PARALLEL LOOP collapse(2) tile(ie:1,q:3) local(Qmax) copyin(min_neigh,max_neigh)
    do ie=nets,nete
      do q=1,my_qsize
        do k=1,my_nlev
          Qmax(:,:,k)=max_neigh(k,q,ie)
        enddo

      elem_desc_reverse_ptr            = loc(my_elem(ie)%desc%reverse)
      elem_desc_putmapP_ptr            = loc(my_elem(ie)%desc%putmapP)

      !$ACC DATA COPYIN(elem_desc_putmapP, elem_desc_reverse)
      edge_buf_5_ptr    = loc(edgeMinMax_buf(1+my_nlev*(q-1)+my_nlev*my_qsize, elem_desc_putmapP(5)+1))
      edge_buf_6_ptr    = loc(edgeMinMax_buf(1+my_nlev*(q-1)+my_nlev*my_qsize, elem_desc_putmapP(6)+1))
      edge_buf_7_ptr    = loc(edgeMinMax_buf(1+my_nlev*(q-1)+my_nlev*my_qsize, elem_desc_putmapP(7)+1))
      edge_buf_8_ptr    = loc(edgeMinMax_buf(1+my_nlev*(q-1)+my_nlev*my_qsize, elem_desc_putmapP(8)+1))
      edge_buf_in_1_ptr = loc(edgeMinMax_buf(1+my_nlev*(q-1)+my_nlev*my_qsize, elem_desc_putmapP(north)+1))
      edge_buf_in_2_ptr = loc(edgeMinMax_buf(1+my_nlev*(q-1)+my_nlev*my_qsize, elem_desc_putmapP(north)+2))
      edge_buf_in_3_ptr = loc(edgeMinMax_buf(1+my_nlev*(q-1)+my_nlev*my_qsize, elem_desc_putmapP(north)+3))
      edge_buf_in_4_ptr = loc(edgeMinMax_buf(1+my_nlev*(q-1)+my_nlev*my_qsize, elem_desc_putmapP(north)+4))
      edge_buf_is_1_ptr = loc(edgeMinMax_buf(1+my_nlev*(q-1)+my_nlev*my_qsize, elem_desc_putmapP(south)+1))
      edge_buf_is_2_ptr = loc(edgeMinMax_buf(1+my_nlev*(q-1)+my_nlev*my_qsize, elem_desc_putmapP(south)+2))
      edge_buf_is_3_ptr = loc(edgeMinMax_buf(1+my_nlev*(q-1)+my_nlev*my_qsize, elem_desc_putmapP(south)+3))
      edge_buf_is_4_ptr = loc(edgeMinMax_buf(1+my_nlev*(q-1)+my_nlev*my_qsize, elem_desc_putmapP(south)+4))
      edge_buf_ie_1_ptr = loc(edgeMinMax_buf(1+my_nlev*(q-1)+my_nlev*my_qsize, elem_desc_putmapP(east)+1))
      edge_buf_ie_2_ptr = loc(edgeMinMax_buf(1+my_nlev*(q-1)+my_nlev*my_qsize, elem_desc_putmapP(east)+2))
      edge_buf_ie_3_ptr = loc(edgeMinMax_buf(1+my_nlev*(q-1)+my_nlev*my_qsize, elem_desc_putmapP(east)+3))
      edge_buf_ie_4_ptr = loc(edgeMinMax_buf(1+my_nlev*(q-1)+my_nlev*my_qsize, elem_desc_putmapP(east)+4))
      edge_buf_iw_1_ptr = loc(edgeMinMax_buf(1+my_nlev*(q-1)+my_nlev*my_qsize, elem_desc_putmapP(west)+1))
      edge_buf_iw_2_ptr = loc(edgeMinMax_buf(1+my_nlev*(q-1)+my_nlev*my_qsize, elem_desc_putmapP(west)+2))
      edge_buf_iw_3_ptr = loc(edgeMinMax_buf(1+my_nlev*(q-1)+my_nlev*my_qsize, elem_desc_putmapP(west)+3))
      edge_buf_iw_4_ptr = loc(edgeMinMax_buf(1+my_nlev*(q-1)+my_nlev*my_qsize, elem_desc_putmapP(west)+4))

      !$ACC DATA COPY(edge_buf_5,edge_buf_6,edge_buf_7,edge_buf_8) COPYOUT(edge_buf_in_1,edge_buf_in_2,edge_buf_in_3,edge_buf_in_4,edge_buf_is_1,edge_buf_is_2,edge_buf_is_3,edge_buf_is_4,edge_buf_ie_1,edge_buf_ie_2,edge_buf_ie_3,edge_buf_ie_4,edge_buf_iw_1,edge_buf_iw_2,edge_buf_iw_3,edge_buf_iw_4)
      call my_edgeVpack_acc(Qmax(:,:,:) , my_nlev , 0+my_nlev*(q-1) , elem_desc_putmapP(:), &
         elem_desc_reverse, &
         edge_buf_5, edge_buf_6, edge_buf_7, edge_buf_8, &
         edge_buf_in_1, edge_buf_in_2, edge_buf_in_3, edge_buf_in_4, &
         edge_buf_is_1, edge_buf_is_2, edge_buf_is_3, edge_buf_is_4, &
         edge_buf_ie_1, edge_buf_ie_2, edge_buf_ie_3, edge_buf_ie_4, &
         edge_buf_iw_1, edge_buf_iw_2, edge_buf_iw_3, edge_buf_iw_4)
      !$ACC END DATA

      !$ACC END DATA
    end do
  enddo
  !$ACC END PARALLEL LOOP
end subroutine

subroutine neighbor_minmax(elem,hybrid,edgeMinMax,nets,nete,min_neigh,max_neigh)
!
! compute Q min&max over the element and all its neighbors
!
!
use physical_constants, only : rrearth
integer :: nets,nete
type (element_t)     , intent(in) :: elem(:)
type (hybrid_t)      , intent(in) :: hybrid
type (EdgeBuffer_t)  , intent(inout) :: edgeMinMax
real (kind=real_kind) :: min_neigh(nlev,qsize,nets:nete)
real (kind=real_kind) :: max_neigh(nlev,qsize,nets:nete)

! local
integer :: ie,k,q
real (kind=real_kind) :: Qmin(np,np,nlev,qsize)
real (kind=real_kind) :: Qmax(np,np,nlev,qsize)
integer(kind=8) :: count_start, count_stop, count_rate, count_max

    call my_neighbor_minmax_acc_before(nets,nete,min_neigh,max_neigh,elem,edgeMinMax,edgeMinMax%nlyr, edgeMinMax%nbuf, edgeMinMax%buf, qsize)

    call bndry_exchangeV(hybrid,edgeMinMax)

    call my_neighbor_minmax_after_bndry(nets, nete, edgeMinMax%nlyr, edgeMinMax%nbuf, &
                        edgeMinMax%buf, south, east, north, west, qsize, &
                        swest, max_corner_elem, elem, &
                        min_neigh, max_neigh, rrearth)

end subroutine


#else


subroutine neighbor_minmax(elem,hybrid,edgeMinMax,nets,nete,nt,min_neigh,max_neigh,min_var,max_var,kmass)
!
! compute Q min&max over the element and all its neighbors
!
!
integer :: nets,nete,nt
type (hybrid_t)      , intent(in) :: hybrid
type (element_t)     , intent(inout) :: elem(:)
type (EdgeBuffer_t)  , intent(in) :: edgeMinMax
real (kind=real_kind) :: min_neigh(nlev,nets:nete)
real (kind=real_kind) :: max_neigh(nlev,nets:nete)
real (kind=real_kind),optional :: min_var(nlev,nets:nete)
real (kind=real_kind),optional :: max_var(nlev,nets:nete)
real (kind=real_kind) :: Qmin(np,np,nlev)
real (kind=real_kind) :: Qmax(np,np,nlev)
real (kind=real_kind) :: Qvar(np,np,nlev)
type (EdgeBuffer_t)          :: edgebuf
integer, optional :: kmass

! local
integer :: ie,k,q

  if(present(kmass))then
!the check if kmass is a valid number is done in sweq_mod
    do k=1,nlev
      if(k.ne.kmass)then
         do ie=nets,nete
            elem(ie)%state%p(:,:,k,nt)=elem(ie)%state%p(:,:,k,nt)/&
            elem(ie)%state%p(:,:,kmass,nt)
         enddo
      endif
    enddo
  endif



    ! create edge buffer for 3 fields
    call initEdgeBuffer(edgebuf,3*nlev)


    ! compute p min, max
    do ie=nets,nete
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k)
#endif
       do k=1,nlev
          Qmin(:,:,k)=minval(elem(ie)%state%p(:,:,k,nt))
          Qmax(:,:,k)=maxval(elem(ie)%state%p(:,:,k,nt))
          ! max - min - crude approximation to TV within the element:
          Qvar(:,:,k)=Qmax(1,1,k)-Qmin(1,1,k)
       enddo
       call edgeVpack(edgebuf,Qmax,nlev,0,elem(ie)%desc)
       call edgeVpack(edgebuf,Qmin,nlev,nlev,elem(ie)%desc)
       call edgeVpack(edgebuf,Qvar,nlev,2*nlev,elem(ie)%desc)
    enddo

    call bndry_exchangeV(hybrid,edgebuf)

    do ie=nets,nete
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k)
#endif
       do k=1,nlev
          Qmin(:,:,k)=minval(elem(ie)%state%p(:,:,k,nt))
          Qmax(:,:,k)=maxval(elem(ie)%state%p(:,:,k,nt))
       enddo

       ! now unpack the min
       if (present(min_var)) then
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k)
#endif
          do k=1,nlev
             Qvar(:,:,k)=Qmax(1,1,k)-Qmin(1,1,k)
          enddo
! WARNING - edgeVunpackMin/Max take second argument as input/ouput
          call edgeVunpackMin(edgebuf,Qvar,nlev,2*nlev,elem(ie)%desc)
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k)
#endif
          do k=1,nlev
             min_var(k,ie)=minval(Qvar(:,:,k))
          enddo
       endif

       ! now unpack the max
       if (present(max_var)) then
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k)
#endif
          do k=1,nlev
             Qvar(:,:,k)=Qmax(1,1,k)-Qmin(1,1,k)
          enddo
! WARNING - edgeVunpackMin/Max take second argument as input/ouput
          call edgeVunpackMax(edgebuf,Qvar,nlev,2*nlev,elem(ie)%desc)
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k)
#endif
          do k=1,nlev
             max_var(k,ie)=maxval(Qvar(:,:,k))
          enddo
       endif


! WARNING - edgeVunpackMin/Max take second argument as input/ouput
       call edgeVunpackMax(edgebuf,Qmax,nlev,0,elem(ie)%desc)
       call edgeVunpackMin(edgebuf,Qmin,nlev,nlev,elem(ie)%desc)
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k)
#endif
       do k=1,nlev
          max_neigh(k,ie)=maxval(Qmax(:,:,k))
          min_neigh(k,ie)=minval(Qmin(:,:,k))
       enddo

    end do

    call FreeEdgeBuffer(edgebuf)
#ifdef DEBUGOMP
#if (! defined ELEMENT_OPENMP)
!$OMP BARRIER
#endif
#endif

  if(present(kmass))then
    do k=1,nlev
       if(k.ne.kmass)then
          do ie=nets,nete
             elem(ie)%state%p(:,:,k,nt)=elem(ie)%state%p(:,:,k,nt)*&
             elem(ie)%state%p(:,:,kmass,nt)
          enddo
       endif
    enddo
  endif
end subroutine
#endif





  subroutine test_ibyp(elem, hybrid,  nets,   nete)
!
! Note: vector test functions should be co-variant since u is contra-variant
!  PHIvec = PHIcov  (test function)
!  PHIcon = DtD PHIcov
!
! weak grad:
!  < PHIcov du/dt > = < PHIcon grad(p) >    (output of grad is covariant)
!  < PHIcov du/dt > = -< div(PHIcon) p >    (input of div is contra)
!  verify:
!    gradient_sphere_wk(p) = - <div(PHIcon) p >
!    gradient_sphere_wk(p) + MASS*grad(p) = b.c. (b.c. are covariant)
!
! weak div:
!   < PHI div(u) > = -< grad(PHI) dot u >     u=contra, output of grad is covariant
! verify:
!   divergence_sphere_wk(u) = -<grad(PHI) dot u>
!   divergence_sphere_wk(u) + MASS*div(u) = b.c.  (b.c. are scalars)
!
! weak curl:
!  < PHIcov du/dt > = < PHIcov curl( a ) >    (output of curl is contra)
!  < PHIcov du/dt > = < vor(PHIcov) a >       (input to vor is covariant)
! verify:
!    curl_sphere_wk(a) = < vor(PHIcov) a >
!    curl_sphere_wk(a) - MASS*curl(a) = b.c. (b.c. are contra)
!
    ! ---------------------
    use kinds, only : real_kind
    ! ---------------------
    use physical_constants, only : rearth
    ! ---------------------
    use dimensions_mod, only : np, nlev
    ! ---------------------
    use element_mod, only : element_t
    ! ---------------------
    use hybrid_mod, only : hybrid_t
    ! ---------------------
    use derivative_mod, only : derivative_t, gradient_sphere, divergence_sphere,vorticity_sphere,&
                               divergence_sphere_wk, curl_sphere
    use global_norms_mod

    implicit none

    type (element_t)     , intent(inout), target :: elem(:)

    type (hybrid_t)      , intent(in) :: hybrid

    integer              , intent(in) :: nets
    integer              , intent(in) :: nete

#if 0
#undef CURLGRAD_TEST
#define IBYP_TEST
    ! =================
    ! Local
    ! =================
    ! pointer ...
    real (kind=real_kind), dimension(:,:), pointer :: rspheremv,spheremv

    ! Thread private working set ...

    real (kind=real_kind), dimension(np,np,nets:nete) :: ptens
    real (kind=real_kind), dimension(np,np,nets:nete) :: ptens2
    real (kind=real_kind), dimension(np,np,nets:nete) :: ptens3

    real (kind=real_kind), dimension(np,np,2,nets:nete)    :: pv      ! p*v lat-lon
    real (kind=real_kind), dimension(np,np,nets:nete)            :: E          ! kinetic energy term
    real (kind=real_kind), dimension(np,np,nets:nete)  :: divbig
    real (kind=real_kind), dimension(np,np,nets:nete)  :: wdivbig
    real (kind=real_kind), dimension(np,np,2,nets:nete)    :: gradbig

    real (kind=real_kind), dimension(np,np,2)      :: ulatlon
    real (kind=real_kind), dimension(np,np,2)      :: grade
    real (kind=real_kind), dimension(np,np,2)      :: grade2
    real (kind=real_kind), dimension(np,np)      :: vor
    real (kind=real_kind), dimension(np,np)      :: div
    real (kind=real_kind), dimension(np,np)      :: wdiv

    real (kind=real_kind) ::  v1,v2,v3

    real*8                :: st,et, time_adv
    integer    :: i,j,k,ie,iie,ii,jj
    integer    :: kptr
    integer    :: nm1,n0,np1
    integer    :: nstep

    type (derivative_t)          :: deriv
    call derivinit(deriv)



    ! ===================================
    ! construct test function
    ! ===================================

    do iie=nets,nete
       do jj=1,np
          do ii=1,np
             ! test for cardinal function at iie,jj,ii

    write(iulog,'(a,3i4,2e20.10)') 'carinal function:  ',iie,ii,jj

    ! construct two global cardinal functions  pv and E
    do ie=nets,nete
       do j=1,np
          do i=1,np


             E(i,j,ie)=0
#ifdef CURLGRAD_TEST
             if (ie==iie .and. i==ii .and. j==jj) E(i,j,ie)=1
#else
             if (ie==1 .and. i==1 .and. j==1) E(i,j,ie)=1
             if (ie==1 .and. i==1 .and. j==2) E(i,j,ie)=1
             if (ie==1 .and. i==2 .and. j==1) E(i,j,ie)=1
             if (ie==1 .and. i==3 .and. j==3) E(i,j,ie)=1
#endif

             ! delta function in contra coordinates
             v1     = 0
             v2     = 0
             if (ie==iie .and. i==ii .and. j==jj) then
                !v1=rearth
                v2=rearth
             endif

             ulatlon(i,j,1)=elem(ie)%D(1,1,i,j)*v1 + elem(ie)%D(1,2,i,j)*v2   ! contra->latlon
             ulatlon(i,j,2)=elem(ie)%D(2,1,i,j)*v1 + elem(ie)%D(2,2,i,j)*v2   ! contra->latlon
             pv(i,j,1,ie) = ulatlon(i,j,1)
             pv(i,j,2,ie) = ulatlon(i,j,2)

          end do
       end do
    enddo
    call make_C0(E,elem,hybrid,nets,nete)
    call make_C0_vector(pv,elem,hybrid,nets,nete)


#ifdef CURLGRAD_TEST
    ! check curl(grad(E))
    do ie=nets,nete
       if ( maxval(abs(E(:,:,ie))) > 0 ) then
       !write(iulog,'(a,i4,2e20.10)') 'maxval: E =',ie,maxval(E(:,:,ie))
       grade=curl_sphere(E(:,:,ie),deriv,elem(ie))
       div=divergence_sphere(grade,deriv,elem(ie))
       vor=vorticity_sphere(grade,deriv,elem(ie))
       if (maxval(abs(div))*rearth**2 > .2e-11) then
          write(iulog,'(a,i4,2e20.10)') 'maxval: div(curl),  vor(curl)=',ie,maxval(abs(div))*rearth**2,maxval(abs(vor))*rearth**2
       endif

       grade=gradient_sphere(E(:,:,ie),deriv,elem(ie)%Dinv)
       vor=vorticity_sphere(grade,deriv,elem(ie))
       div=divergence_sphere(grade,deriv,elem(ie))
       if (maxval(abs(vor))*rearth**2 > .2e-11) then
          write(iulog,'(a,i4,2e20.10)') 'maxval: curl(grad), div(grad)=',ie,maxval(abs(vor))*rearth**2,maxval(abs(div))*rearth**2
       endif
       endif
    enddo

    ! check div(curl(E)) with DSS
    do ie=nets,nete
       gradbig(:,:,:,ie)=curl_sphere(E(:,:,ie),deriv,elem(ie))
    enddo
    call make_C0_vector(gradbig,elem,hybrid,nets,nete)
    do ie=nets,nete
       divbig(:,:,ie)=divergence_sphere(gradbig(:,:,:,ie),deriv,elem(ie))
    enddo
    call make_C0(divbig,elem,hybrid,nets,nete)
    do ie=nets,nete
       if (maxval(abs(divbig(:,:,ie)))*rearth**2 > .8e-12) then
          write(iulog,'(a,i4,2e20.10)') 'maxval: [div([curl])]=',ie,maxval(abs(divbig(:,:,ie)))*rearth**2
       endif
    enddo


    ! check curl(grad(E)) with DSS
    do ie=nets,nete
       gradbig(:,:,:,ie)=gradient_sphere(E(:,:,ie),deriv,elem(ie)%Dinv)
    enddo
    call make_C0_vector(gradbig,elem,hybrid,nets,nete)
    do ie=nets,nete
       divbig(:,:,ie)=vorticity_sphere(gradbig(:,:,:,ie),deriv,elem(ie))
    enddo
    call make_C0(divbig,elem,hybrid,nets,nete)

    do ie=nets,nete
       if (maxval(abs(divbig(:,:,ie)))*rearth**2 > .8e-12) then
          write(iulog,'(a,i4,2e20.10)') 'maxval: [curl([gradl])]=',ie,maxval(abs(divbig(:,:,ie)))*rearth**2
       endif
    enddo
#endif


#ifdef IBYP_TEST
    ! compare <grad(E) dot pv> and <E div(pv)>  < E weak_div(pv) >
    v2=0
    do ie=nets,nete
       spheremv     => elem(ie)%spheremp(:,:)

       div = divergence_sphere(pv(1,1,1,ie),deriv,elem(ie))      ! latlon vector -> scalar
       grade = gradient_sphere(E(1,1,ie),deriv,elem(ie)%Dinv)
       wdiv = divergence_sphere_wk(pv(1,1,1,ie),deriv,elem(ie))



       do j=1,np
          do i=1,np
!             write(iulog,'(3i3,3e22.14)') ie,i,j,pv(i,j,1,ie),pv(i,j,2,ie),div(i,j)

             ! (grad(E) dot pv )
             ptens3(i,j,ie) = grade(i,j,1)*pv(i,j,1,ie) + grade(i,j,2)*pv(i,j,2,ie)
             v2=v2+wdiv(i,j)*E(i,j,ie)
          end do
       end do
       ptens(:,:,ie)=div(:,:)*E(:,:,ie)   ! < E div(pv) >
       ptens2(:,:,ie)=wdiv(:,:)*E(:,:,ie)/spheremv(:,:)   ! < wdiv E >
       ! ===================================================
       ! Pack cube edges of tendencies, rotate velocities
       ! ===================================================
       divbig(:,:,ie)=div(:,:)*spheremv(:,:)
       wdivbig(:,:,ie)=wdiv(:,:)
    end do
    call make_C0(divbig,elem,hybrid,nets,nete)
    call make_C0(wdivbig,elem,hybrid,nets,nete)


!!$    v1=global_integral(elem,ptens,hybrid,np,nets,nete)
!!$    v3=global_integral(elem,ptens3,hybrid,np,nets,nete)
!!$    print *,'< E div(pv) >   =',v1
!!$    print *,'< E div_wk(pv) >=',v2/(4*4*atan(1d0))
!!$    v2=global_integral(elem,ptens2,hybrid,np,nets,nete)
!!$    print *,'< E div_wk(pv) >=',v2
!!$    print *,'-<grad(E),pv >  =',-v3
!!$!    print *,'sum1-sum2/max',(v1-v2)/max(abs(v1),abs(v2))
!!$


    do ie=nets,nete
       div(:,:)=divbig(:,:,ie)
       wdiv(:,:)=wdivbig(:,:,ie)
       ! ===========================================================
       ! Compute velocity and pressure tendencies for all levels
       ! ===========================================================
          do j=1,np
             do i=1,np
                ! < div(pv) >   vs < div_wk(pv) >
                if ( abs(div(i,j)-wdiv(i,j)) > .15e-17) then
!                   write(iulog,'(3i3,4e22.14)') ie,i,j,div(i,j),wdiv(i,j),div(i,j)-wdiv(i,j),E(i,j,ie)
                endif
                if ( E(i,j,ie)/=0 ) then
!                   write(iulog,'(3i3,4e22.14)') ie,i,j,div(i,j),wdiv(i,j),div(i,j)-wdiv(i,j),E(i,j,ie)
                endif

             end do
          end do


    end do
#ifdef DEBUGOMP
#if (! defined ELEMENT_OPENMP)
!$OMP BARRIER
#endif
#endif
    write(iulog,'(a,3i4,2e20.10)') 'max diff div-wdiv: ',iie,ii,jj,maxval(abs(divbig(:,:,:)-wdivbig(:,:,:))),maxval(divbig(:,:,:))
#endif


    enddo
    enddo
    enddo
    stop
#endif
  end subroutine test_ibyp









  subroutine check_edge_flux(elem,deriv,nets,nete)
!
!  check local element vector dentities:
!*****
!  1. div and weak div are adjoints: (for all scalar test functions)
!     integral[  p div(u) ] + integral[ grad(p) dot u ] = boundary_integral[ p u dot n]
!       PHI div(u) spheremp - div_wk(u)(i,j) = boundary_integral[ u PHI]
!       where PHI = the delta function at (i,j)
!
!*****
!  2. grad and weak grad are adjoints:
!     weak gradient is defined with CONTRA vector test functions
!     i.e. it returns vector:   [  integral[ p div(PHIcontra_1) ]
!                               [  integral[ p div(PHIcontra_2) ]
!
!   integral[  p div(u) ] + integral[ grad(p) dot u ] = boundary_integral[ p u dot n]
! take u = PHIcontra_1 = (1,0) vector delta funciton at (i,j):
!  -grad_wk(p)_1(i,j) + spheremp PHIcontra_1 dot grad(p) = boundary_integral[ PHIcontra_1 p]
! and then take u = PHIcontra_2 = (0,1) vector delta function at (i,j):
!  -grad_wk(p)_2(i,j) + spheremp PHIcontra_2 dot grad(p) = boundary_integral[ PHIcontra_2 p]
!
! which is an equation for each covariant component:
! -grad_wk(p)_cov1 + spheremp grad(p)_cov1 = boundary_integral[ PHIcontra_1 p dot n]
! -grad_wk(p)_cov2 + spheremp grad(p)_cov2 = boundary_integral[ PHIcontra_2 p dot n]
!
! HOMME-SE works in latlon, so convert cov->lat/lon:
!
! -grad_wk(p) + spheremp grad(p) = D^-t * B
!
! with
!    B1 = boundary_integral[ PHIcontra_1 p]
!    B2 = boundary_integral[ PHIcontra_2 p]
!
!*****
! 3.  weak grid with COVARIANT test functions!
!   integral[  p div(u) ] + integral[ grad(p) dot u ] = boundary_integral[ p u dot n]
! take u = PHIcov_1 = (1,0) vector delta funciton at (i,j):
!  -grad_wk(p)_1(i,j) + spheremp PHIcov_1 dot grad(p) = boundary_integral[ PHIcov_1 p]
! and then take u = PHIcov_2 = (0,1) vector delta function at (i,j):
!  -grad_wk(p)_2(i,j) + spheremp PHIcov_2 dot grad(p) = boundary_integral[ PHIcov_2 p]
!
! which is an equation for each CONTRA component:
! -grad_wk(p)_contra1 + spheremp grad(p)_contra1 = B1
! -grad_wk(p)_contra2 + spheremp grad(p)_contra2 = B2
!
! HOMME-SE works in latlon, so convert contra ->lat/lon:
!
! -grad_wk(p) + spheremp grad(p) = D * B
!
! with
!    B1 = boundary_integral[ PHIcov_1 p]
!    B2 = boundary_integral[ PHIcov_2 p]
!
!*****
! 4.  weak curl with COVARIANT test functions!
!  integral[ u dot curl(v)] - integral[v dot curl(u)] = boundary_integral[ v cross u dot n]
!  curl(p) = curl(p*khat) = horizontal vector
!  vor(U) =  s*khat       = (which we treat as a scalar)
!   integral[ p * vor(u)  ] - integral[ u dot curl(p) ] = boundary_integral[ u cross p*khat  dot n]
!
! take u = PHIcov_1 = (1,0) vector delta funciton at (i,j):
!   curl_wk(p)_1(i,j) - spheremp PHIcov_1 dot curl(p) = boundary_integral[ perp(PHIcov_1) p]
! and then take u = PHIcov_2 = (0,1) vector delta function at (i,j):
!   curl_wk(p)_2(i,j) - spheremp PHIcov_2 dot curl(p) = boundary_integral[ perp(PHIcov_2) p]
!
! which is an equation for each CONTRA component:
! curl_wk(p)_contra1 - spheremp curl(p)_contra1 = B1
! curl_wk(p)_contra2 - spheremp curl(p)_contra2 = B2
!
! HOMME-SE works in latlon, so convert contra ->lat/lon:
!
! curl_wk(p) + spheremp curl(p) = D * B
!
! with
!    B1 = boundary_integral[ PHIcov_1 p]
!    B2 = boundary_integral[ PHIcov_2 p]
!
  use dimensions_mod, only : np, np, nlev
  use element_mod, only    : element_t
  use derivative_mod, only  : derivative_t, divergence_sphere, divergence_sphere_wk, &
                             element_boundary_integral, gradient_sphere, &
                             gradient_sphere_wk_testcontra,gradient_sphere_wk_testcov, &
                             curl_sphere, curl_sphere_wk_testcov
  use physical_constants, only : rrearth

  implicit none

  type (element_t)     , intent(inout), target :: elem(:)
  type (derivative_t)  , intent(in) :: deriv
  integer :: nets,nete
  ! local
  real (kind=real_kind), dimension(np,np,2) :: ucontra,ulatlon,gradp,gradp_wk,ucov
  real (kind=real_kind), dimension(np,np) :: phidivu,ugradphi,rhs,lhs,p
  real (kind=real_kind), dimension(np,np) :: rhs2,lhs2
  integer :: i,j,ie

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  print *,'integration by parts identity: check div/weak div:'
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! test integration by parts identity for each Cardinal function PHI:
  ! div(u)*spheremp - div_wk(u) = boundary integral phi u dot n
  do ie=nets,nete
     call random_number(ucontra)
     ! contra->latlon
     ulatlon(:,:,1)=(elem(ie)%D(1,1,:,:)*ucontra(:,:,1) + elem(ie)%D(1,2,:,:)*ucontra(:,:,2))
     ulatlon(:,:,2)=(elem(ie)%D(2,1,:,:)*ucontra(:,:,1) + elem(ie)%D(2,2,:,:)*ucontra(:,:,2))
     phidivu = elem(ie)%spheremp(:,:)*divergence_sphere(ulatlon,deriv,elem(ie))
     ugradphi = divergence_sphere_wk(ulatlon,deriv,elem(ie))
     lhs = phidivu - ugradphi

     rhs = element_boundary_integral(ulatlon,deriv,elem(ie))


     do j=1,np
        do i=1,np
           if ( abs(lhs(i,j)-rhs(i,j)) .gt. 1d-20) then
              write(*,'(a)') 'ERROR: div/div_wk integration by parts failure!'
              write(*,'(a,2i3,a,3e12.5)') 'for test function (i,j)=',i,j,' lhs,rhs=',lhs(i,j),rhs(i,j),lhs(i,j)-rhs(i,j)
           endif
        enddo
     enddo
  enddo

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  print *,'check grad/weak grad (gradient_sphere_wk_testcontra)'
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! PHIVEC = contra cardinal function
  !          check each contra component seperately

  do ie=nets,nete
     call random_number(p)

     ! grad(p)  (lat/lon vector)
     gradp = gradient_sphere(p,deriv,elem(ie)%Dinv)
     gradp(:,:,1)=gradp(:,:,1)*elem(ie)%spheremp(:,:)
     gradp(:,:,2)=gradp(:,:,2)*elem(ie)%spheremp(:,:)
     gradp_wk = gradient_sphere_wk_testcontra(p,deriv,elem(ie))

     ucontra(:,:,1)=p  ! PHIvec_1 * p
     ucontra(:,:,2)=0
     ! contra->latlon
     ulatlon(:,:,1)=(elem(ie)%D(1,1,:,:)*ucontra(:,:,1) + elem(ie)%D(1,2,:,:)*ucontra(:,:,2))
     ulatlon(:,:,2)=(elem(ie)%D(2,1,:,:)*ucontra(:,:,1) + elem(ie)%D(2,2,:,:)*ucontra(:,:,2))

     rhs = element_boundary_integral(ulatlon,deriv,elem(ie))
     lhs = gradp(:,:,1)-gradp_wk(:,:,1)

     ucontra(:,:,1)=0  ! PHIvec_2 * p
     ucontra(:,:,2)=p
     ! contra->latlon
     ulatlon(:,:,1)=(elem(ie)%D(1,1,:,:)*ucontra(:,:,1) + elem(ie)%D(1,2,:,:)*ucontra(:,:,2))
     ulatlon(:,:,2)=(elem(ie)%D(2,1,:,:)*ucontra(:,:,1) + elem(ie)%D(2,2,:,:)*ucontra(:,:,2))
     rhs2 = element_boundary_integral(ulatlon,deriv,elem(ie))
     lhs2 = gradp(:,:,2)-gradp_wk(:,:,2)


     ! boundary integral gives covariant components. (see above) convert to latlon:
     ! cov -> latlon
     gradp(:,:,1)=rhs
     gradp(:,:,2)=rhs2
     rhs(:,:)=elem(ie)%Dinv(1,1,:,:)*gradp(:,:,1) + elem(ie)%Dinv(2,1,:,:)*gradp(:,:,2)
     rhs2(:,:)=elem(ie)%Dinv(1,2,:,:)*gradp(:,:,1) + elem(ie)%Dinv(2,2,:,:)*gradp(:,:,2)


     do j=1,np
        do i=1,np
           if ( abs(lhs(i,j)-rhs(i,j)) .gt. 1d-20) then
              write(*,'(a)') 'ERROR: grad/grad_wk CONTRA (1) integration by parts failure!'
              write(*,'(a,2i3,a,4e12.4)') '(i,j)=',i,j,' lhs,rhs=',lhs(i,j),rhs(i,j),&
                   lhs(i,j)-rhs(i,j),lhs(i,j)/rhs(i,j)
           endif
        enddo
     enddo


     do j=1,np
        do i=1,np
           if ( abs(lhs2(i,j)-rhs2(i,j)) .gt. 1d-20) then
              write(*,'(a)') 'ERROR: grad/grad_wk CONTRA (2) integration by parts failure!'
              write(*,'(a,2i2,a,3e12.4)') '(i,j)=',i,j,' lhs,rhs=',lhs2(i,j),rhs2(i,j),lhs2(i,j)-rhs2(i,j)
           endif
        enddo
     enddo
  enddo

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  print *,'check grad/weak grad (gradient_sphere_wk_testcov)'
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  do ie=nets,nete
     call random_number(p)


     ! grad(p)  (lat/lon vector)
     gradp = gradient_sphere(p,deriv,elem(ie)%Dinv)
     gradp(:,:,1)=gradp(:,:,1)*elem(ie)%spheremp(:,:)
     gradp(:,:,2)=gradp(:,:,2)*elem(ie)%spheremp(:,:)
     gradp_wk = gradient_sphere_wk_testcov(p,deriv,elem(ie))
     lhs = gradp(:,:,1)-gradp_wk(:,:,1)
     lhs2 = gradp(:,:,2)-gradp_wk(:,:,2)

     ucov(:,:,1)=p  ! PHIvec_1 * p
     ucov(:,:,2)=0
     ! cov->latlon
     ulatlon(:,:,1)=(elem(ie)%Dinv(1,1,:,:)*ucov(:,:,1) + elem(ie)%Dinv(2,1,:,:)*ucov(:,:,2))
     ulatlon(:,:,2)=(elem(ie)%Dinv(1,2,:,:)*ucov(:,:,1) + elem(ie)%Dinv(2,2,:,:)*ucov(:,:,2))
     rhs = element_boundary_integral(ulatlon,deriv,elem(ie))

     ucov(:,:,1)=0  ! PHIvec_2 * p
     ucov(:,:,2)=p
     ! cov->latlon
     ulatlon(:,:,1)=(elem(ie)%Dinv(1,1,:,:)*ucov(:,:,1) + elem(ie)%Dinv(2,1,:,:)*ucov(:,:,2))
     ulatlon(:,:,2)=(elem(ie)%Dinv(1,2,:,:)*ucov(:,:,1) + elem(ie)%Dinv(2,2,:,:)*ucov(:,:,2))
     rhs2 = element_boundary_integral(ulatlon,deriv,elem(ie))


     ! boundary integral gives contra components. (see above) convert to latlon:
     ! contra -> latlon
     gradp(:,:,1)=rhs
     gradp(:,:,2)=rhs2
     rhs(:,:) =elem(ie)%D(1,1,:,:)*gradp(:,:,1) + elem(ie)%D(1,2,:,:)*gradp(:,:,2)
     rhs2(:,:)=elem(ie)%D(2,1,:,:)*gradp(:,:,1) + elem(ie)%D(2,2,:,:)*gradp(:,:,2)


     do j=1,np
        do i=1,np
           if ( abs(lhs(i,j)-rhs(i,j)) .gt. 1d-20) then
              write(*,'(a)') 'ERROR: grad/grad_wk COV (1) integration by parts failure!'
              write(*,'(a,2i2,a,4e12.4)') '(i,j)=',i,j,' lhs,rhs=',lhs(i,j),rhs(i,j),lhs(i,j)-rhs(i,j),lhs(i,j)/rhs(i,j)
           endif
        enddo
     enddo

     do j=1,np
        do i=1,np
           if ( abs(lhs2(i,j)-rhs2(i,j)) .gt. 1d-20) then
              write(*,'(a)') 'ERROR: grad/grad_wk COV (2) integration by parts failure!'
              write(*,'(a,2i2,a,3e12.4)') '(i,j)=',i,j,' lhs,rhs=',lhs2(i,j),rhs2(i,j),lhs2(i,j)-rhs2(i,j)
           endif
        enddo
     enddo
  enddo

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  print *,'check curl/weak curl (curl_sphere_wk_testcov)'
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  do ie=nets,nete
     call random_number(p)

     ! grad(p)  (lat/lon vector)
     gradp = curl_sphere(p,deriv,elem(ie))
     gradp(:,:,1)=gradp(:,:,1)*elem(ie)%spheremp(:,:)
     gradp(:,:,2)=gradp(:,:,2)*elem(ie)%spheremp(:,:)
     gradp_wk = curl_sphere_wk_testcov(p,deriv,elem(ie))
     lhs =  gradp_wk(:,:,1)-gradp(:,:,1)
     lhs2 = gradp_wk(:,:,2)-gradp(:,:,2)

     ucov(:,:,1)=p  ! PHIvec_1 * p
     ucov(:,:,2)=0
     ! cov->latlon, and then u cross khat:
     ulatlon(:,:,2)=-(elem(ie)%Dinv(1,1,:,:)*ucov(:,:,1) + elem(ie)%Dinv(2,1,:,:)*ucov(:,:,2))
     ulatlon(:,:,1)= (elem(ie)%Dinv(1,2,:,:)*ucov(:,:,1) + elem(ie)%Dinv(2,2,:,:)*ucov(:,:,2))
     rhs = element_boundary_integral(ulatlon,deriv,elem(ie))

     ucov(:,:,1)=0  ! PHIvec_2 * p
     ucov(:,:,2)=p
     ! cov->latlon, and u cross khat:
     ulatlon(:,:,2)=-(elem(ie)%Dinv(1,1,:,:)*ucov(:,:,1) + elem(ie)%Dinv(2,1,:,:)*ucov(:,:,2))
     ulatlon(:,:,1)= (elem(ie)%Dinv(1,2,:,:)*ucov(:,:,1) + elem(ie)%Dinv(2,2,:,:)*ucov(:,:,2))
     rhs2 = element_boundary_integral(ulatlon,deriv,elem(ie))


     ! boundary integral gives contra components. (see above) convert to latlon:
     ! contra -> latlon
     gradp(:,:,1)=rhs
     gradp(:,:,2)=rhs2
     rhs(:,:) =elem(ie)%D(1,1,:,:)*gradp(:,:,1) + elem(ie)%D(1,2,:,:)*gradp(:,:,2)
     rhs2(:,:)=elem(ie)%D(2,1,:,:)*gradp(:,:,1) + elem(ie)%D(2,2,:,:)*gradp(:,:,2)


     do j=1,np
        do i=1,np
           if ( abs(lhs(i,j)-rhs(i,j)) .gt. 1d-20) then
              write(*,'(a)') 'ERROR: curl/curl_wk COV (1) integration by parts failure!'
              write(*,'(a,2i2,a,4e12.4)') '(i,j)=',i,j,' lhs,rhs=',lhs(i,j),rhs(i,j),lhs(i,j)-rhs(i,j),lhs(i,j)/rhs(i,j)
           endif
        enddo
     enddo
     stop

     do j=1,np
        do i=1,np
           if ( abs(lhs2(i,j)-rhs2(i,j)) .gt. 1d-20) then
              write(*,'(a)') 'ERROR: curl/curl_wk COV (2) integration by parts failure!'
              write(*,'(a,2i2,a,3e12.4)') '(i,j)=',i,j,' lhs,rhs=',lhs2(i,j),rhs2(i,j),lhs2(i,j)-rhs2(i,j)
           endif
        enddo
     enddo
  enddo

  print *,'done. integration by parts identity check:'
  end subroutine
end module

! TPP Processed 06af35b4b6e1d5ae0884d64cc00c796d

! TPP Processed 06af35b4b6e1d5ae0884d64cc00c796d

! TPP Processed 06af35b4b6e1d5ae0884d64cc00c796d

! TPP Processed 06af35b4b6e1d5ae0884d64cc00c796d
