#include "REAL.H"
#include "SPACE.H"
#include "CONSTANTS.H"

#include "CONSTANTS.H"
      subroutine APPLYLAP_4D(
     &           lofphi
     &           ,ilofphilo0,ilofphilo1,ilofphilo2,ilofphilo3
     &           ,ilofphihi0,ilofphihi1,ilofphihi2,ilofphihi3
     &           ,nlofphicomp
     &           ,flux
     &           ,ifluxlo0,ifluxlo1,ifluxlo2,ifluxlo3
     &           ,ifluxhi0,ifluxhi1,ifluxhi2,ifluxhi3
     &           ,nfluxcomp
     &           ,phi
     &           ,iphilo0,iphilo1,iphilo2,iphilo3
     &           ,iphihi0,iphihi1,iphihi2,iphihi3
     &           ,nphicomp
     &           ,iregionlo0,iregionlo1,iregionlo2,iregionlo3
     &           ,iregionhi0,iregionhi1,iregionhi2,iregionhi3
     &           ,ifluxregionlo0,ifluxregionlo1,ifluxregionlo2,ifluxregionlo3
     &           ,ifluxregionhi0,ifluxregionhi1,ifluxregionhi2,ifluxregionhi3
     &           ,dx
     &           ,idir
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer nlofphicomp
      integer ilofphilo0,ilofphilo1,ilofphilo2,ilofphilo3
      integer ilofphihi0,ilofphihi1,ilofphihi2,ilofphihi3
      REAL_T lofphi(
     &           ilofphilo0:ilofphihi0,
     &           ilofphilo1:ilofphihi1,
     &           ilofphilo2:ilofphihi2,
     &           ilofphilo3:ilofphihi3,
     &           0:nlofphicomp-1)
      integer nfluxcomp
      integer ifluxlo0,ifluxlo1,ifluxlo2,ifluxlo3
      integer ifluxhi0,ifluxhi1,ifluxhi2,ifluxhi3
      REAL_T flux(
     &           ifluxlo0:ifluxhi0,
     &           ifluxlo1:ifluxhi1,
     &           ifluxlo2:ifluxhi2,
     &           ifluxlo3:ifluxhi3,
     &           0:nfluxcomp-1)
      integer nphicomp
      integer iphilo0,iphilo1,iphilo2,iphilo3
      integer iphihi0,iphihi1,iphihi2,iphihi3
      REAL_T phi(
     &           iphilo0:iphihi0,
     &           iphilo1:iphihi1,
     &           iphilo2:iphihi2,
     &           iphilo3:iphihi3,
     &           0:nphicomp-1)
      integer iregionlo0,iregionlo1,iregionlo2,iregionlo3
      integer iregionhi0,iregionhi1,iregionhi2,iregionhi3
      integer ifluxregionlo0,ifluxregionlo1,ifluxregionlo2,ifluxregionlo3
      integer ifluxregionhi0,ifluxregionhi1,ifluxregionhi2,ifluxregionhi3
      REAL_T dx
      integer idir
      REAL_T dxinv
      integer n,ncomp
      integer ii,i,jj,j,kk,k,ll,l
      ncomp = nphicomp
      if ((ncomp .ne. nlofphicomp).or.(ncomp.ne.nfluxcomp)) then
      endif
      dxinv = one/dx
      
      ii = CHF_ID(idir, 0)
      jj = CHF_ID(idir, 1)
      kk = CHF_ID(idir, 2)
      ll = CHF_ID(idir, 3)
      do n = 0, ncomp-1
         
      do l = ifluxregionlo3,ifluxregionhi3
      do k = ifluxregionlo2,ifluxregionhi2
      do j = ifluxregionlo1,ifluxregionhi1
      do i = ifluxregionlo0,ifluxregionhi0

             flux(i,j,k,l,n) =
     &           ( -phi(i   ,j   ,k   ,l   ,n)
     &           + (phi(i-ii,j-jj,k-kk,l-ll,n))
     &           )*dxinv
         
      enddo
      enddo
      enddo
      enddo
         
      do l = iregionlo3,iregionhi3
      do k = iregionlo2,iregionhi2
      do j = iregionlo1,iregionhi1
      do i = iregionlo0,iregionhi0

            lofphi(i,j,k,l,n) =
     &           lofphi(i,j,k,l,n) -
     &           ( (flux(i+ii,j+jj,k+kk,l+ll,n)
     &            - flux(i   ,j   ,k   ,l   ,n))
     &           )*dxinv
         
      enddo
      enddo
      enddo
      enddo
       enddo
      return
      end
      subroutine HOAVGDOWN_4D(
     &           phi
     &           ,iphilo0,iphilo1,iphilo2,iphilo3
     &           ,iphihi0,iphihi1,iphihi2,iphihi3
     &           ,nphicomp
     &           ,phicoarse
     &           ,iphicoarselo0,iphicoarselo1,iphicoarselo2,iphicoarselo3
     &           ,iphicoarsehi0,iphicoarsehi1,iphicoarsehi2,iphicoarsehi3
     &           ,nphicoarsecomp
     &           ,nrefine
     &           ,iregionlo0,iregionlo1,iregionlo2,iregionlo3
     &           ,iregionhi0,iregionhi1,iregionhi2,iregionhi3
     &           ,iavstencillo0,iavstencillo1,iavstencillo2,iavstencillo3
     &           ,iavstencilhi0,iavstencilhi1,iavstencilhi2,iavstencilhi3
     &           ,navstencil
     &           ,ilapstencillo0,ilapstencillo1,ilapstencillo2,ilapstencillo3
     &           ,ilapstencilhi0,ilapstencilhi1,ilapstencilhi2,ilapstencilhi3
     &           ,nlapstencil
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nphicomp
      integer iphilo0,iphilo1,iphilo2,iphilo3
      integer iphihi0,iphihi1,iphihi2,iphihi3
      REAL_T phi(
     &           iphilo0:iphihi0,
     &           iphilo1:iphihi1,
     &           iphilo2:iphihi2,
     &           iphilo3:iphihi3,
     &           0:nphicomp-1)
      integer nphicoarsecomp
      integer iphicoarselo0,iphicoarselo1,iphicoarselo2,iphicoarselo3
      integer iphicoarsehi0,iphicoarsehi1,iphicoarsehi2,iphicoarsehi3
      REAL_T phicoarse(
     &           iphicoarselo0:iphicoarsehi0,
     &           iphicoarselo1:iphicoarsehi1,
     &           iphicoarselo2:iphicoarsehi2,
     &           iphicoarselo3:iphicoarsehi3,
     &           0:nphicoarsecomp-1)
      integer nrefine
      integer iregionlo0,iregionlo1,iregionlo2,iregionlo3
      integer iregionhi0,iregionhi1,iregionhi2,iregionhi3
      integer iavstencillo0,iavstencillo1,iavstencillo2,iavstencillo3
      integer iavstencilhi0,iavstencilhi1,iavstencilhi2,iavstencilhi3
      integer navstencil
      integer ilapstencillo0,ilapstencillo1,ilapstencillo2,ilapstencillo3
      integer ilapstencilhi0,ilapstencilhi1,ilapstencilhi2,ilapstencilhi3
      integer nlapstencil
      REAL_T lofphi,avlphi
      integer n,ncomp,spacedim
      integer i,ic,is,j,jc,js,k,kc,ks,l,lc,ls
      ncomp = nphicomp
      if (ncomp .ne. nphicoarsecomp) then
      endif
      ncomp = nphicomp
      spacedim = CH_SPACEDIM
      avlphi = 0.0
      do n = 0, ncomp-1
         
      do lc = iregionlo3,iregionhi3
      do kc = iregionlo2,iregionhi2
      do jc = iregionlo1,iregionhi1
      do ic = iregionlo0,iregionhi0

         
      do ls = ilapstencillo3,ilapstencilhi3
      do ks = ilapstencillo2,ilapstencilhi2
      do js = ilapstencillo1,ilapstencilhi1
      do is = ilapstencillo0,ilapstencilhi0

             
             i = ic*nrefine + is
             j = jc*nrefine + js
             k = kc*nrefine + ks
             l = lc*nrefine + ls
             lofphi =
     &              -2*spacedim*phi(i,j,k,l,n) 
     &           +  phi(i+1 ,j   ,k   ,l  ,n)
     &           +  phi(i-1 ,j   ,k   ,l  ,n) 
     &           +  phi(i   ,j+1 ,k   ,l  ,n)
     &           +  phi(i   ,j-1 ,k   ,l  ,n) 
     &           +  phi(i   ,j   ,k+1 ,l  ,n)
     &           +  phi(i   ,j   ,k-1 ,l  ,n) 
     &           +  phi(i   ,j   ,k   ,l+1,n)
     &           +  phi(i   ,j   ,k   ,l-1,n) 
             avlphi = avlphi + lofphi
         
      enddo
      enddo
      enddo
      enddo
         avlphi = avlphi / nlapstencil
         phicoarse(ic,jc,kc,lc,n) = zero
         
      do ls = iavstencillo3,iavstencilhi3
      do ks = iavstencillo2,iavstencilhi2
      do js = iavstencillo1,iavstencilhi1
      do is = iavstencillo0,iavstencilhi0

             
             i = ic*nrefine + is
             j = jc*nrefine + js
             k = kc*nrefine + ks
             l = lc*nrefine + ls
             phicoarse(ic,jc,kc,lc,n) =
     &                      phicoarse(ic,jc,kc,lc,n) +
     &                      phi(i,j,k,l,n)
         
      enddo
      enddo
      enddo
      enddo
             phicoarse(ic,jc,kc,lc,n) =
     &                      phicoarse(ic,jc,kc,lc,n)/navstencil
     &                      -avlphi/24
         
      enddo
      enddo
      enddo
      enddo
       enddo
      return
      end
      subroutine ADDHYPERVISCOUSFLUX_4D(
     &           flux
     &           ,ifluxlo0,ifluxlo1,ifluxlo2,ifluxlo3
     &           ,ifluxhi0,ifluxhi1,ifluxhi2,ifluxhi3
     &           ,nfluxcomp
     &           ,u
     &           ,iulo0,iulo1,iulo2,iulo3
     &           ,iuhi0,iuhi1,iuhi2,iuhi3
     &           ,nucomp
     &           ,iboxlo0,iboxlo1,iboxlo2,iboxlo3
     &           ,iboxhi0,iboxhi1,iboxhi2,iboxhi3
     &           ,mu
     &           ,dx
     &           ,dir
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer nfluxcomp
      integer ifluxlo0,ifluxlo1,ifluxlo2,ifluxlo3
      integer ifluxhi0,ifluxhi1,ifluxhi2,ifluxhi3
      REAL_T flux(
     &           ifluxlo0:ifluxhi0,
     &           ifluxlo1:ifluxhi1,
     &           ifluxlo2:ifluxhi2,
     &           ifluxlo3:ifluxhi3,
     &           0:nfluxcomp-1)
      integer nucomp
      integer iulo0,iulo1,iulo2,iulo3
      integer iuhi0,iuhi1,iuhi2,iuhi3
      REAL_T u(
     &           iulo0:iuhi0,
     &           iulo1:iuhi1,
     &           iulo2:iuhi2,
     &           iulo3:iuhi3,
     &           0:nucomp-1)
      integer iboxlo0,iboxlo1,iboxlo2,iboxlo3
      integer iboxhi0,iboxhi1,iboxhi2,iboxhi3
      REAL_T mu
      REAL_T dx
      integer dir
      integer i,j,k,l, n, s, d
      integer in,jn,kn,ln
      integer it,jt,kt,lt
      REAL_T hvflux, coeff
      coeff = mu
      in = CHF_ID(0,dir)
      jn = CHF_ID(1,dir)
      kn = CHF_ID(2,dir)
      ln = CHF_ID(3,dir)
      do n=0, nucomp-1
         
      do l = iboxlo3,iboxhi3
      do k = iboxlo2,iboxhi2
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0

#if 0
          if (dir.eq.0) then
             hvflux =
     &                       ( u(i+1,j  ,n) - u(i-2,j  ,n) )
     &              -  5.0 * ( u(i  ,j  ,n) - u(i-1,j  ,n) )
     &              +        ( u(i  ,j+1,n) - u(i-1,j+1,n) )
     &              +        ( u(i  ,j-1,n) - u(i-1,j-1,n) )
          else
             hvflux =
     &                       ( u(i  ,j+1,n) - u(i  ,j-2,n) )
     &              -  5.0 * ( u(i  ,j  ,n) - u(i  ,j-1,n) )
     &              +        ( u(i+1,j  ,n) - u(i+1,j-1,n) )
     &              +        ( u(i-1,j  ,n) - u(i-1,j-1,n) )
          endif
          if (dir.eq.0) then
             hvflux =
     &              -        ( u(i+2,j  ,n) - u(i-3,j  ,n) )
     &              +  9.0 * ( u(i+1,j  ,n) - u(i-2,j  ,n) )
     &              - 28.0 * ( u(i  ,j  ,n) - u(i-1,j  ,n) )
     &              -  2.0 * ( u(i+1,j+1,n) - u(i-2,j+1,n) )
     &              + 10.0 * ( u(i  ,j+1,n) - u(i-1,j+1,n) )
     &              -  2.0 * ( u(i+1,j-1,n) - u(i-2,j-1,n) )
     &              + 10.0 * ( u(i  ,j-1,n) - u(i-1,j-1,n) )
     &              -        ( u(i  ,j+2,n) - u(i-1,j+2,n) )
     &              -        ( u(i  ,j-2,n) - u(i-1,j-2,n) )
          else
             hvflux =
     &              -        ( u(i  ,j+2,n) - u(i  ,j-3,n) )
     &              +  9.0 * ( u(i  ,j+1,n) - u(i  ,j-2,n) )
     &              - 28.0 * ( u(i  ,j  ,n) - u(i  ,j-1,n) )
     &              -  2.0 * ( u(i+1,j+1,n) - u(i+1,j-2,n) )
     &              + 10.0 * ( u(i+1,j  ,n) - u(i+1,j-1,n) )
     &              -  2.0 * ( u(i-1,j+1,n) - u(i-1,j-2,n) )
     &              + 10.0 * ( u(i-1,j  ,n) - u(i-1,j-1,n) )
     &              -        ( u(i+2,j  ,n) - u(i+2,j-1,n) )
     &              -        ( u(i-2,j  ,n) - u(i-2,j-1,n) )
          endif
#endif
          hvflux = coeff * hvflux
          flux(i,j,k,l,n) = flux(i,j,k,l,n) + hvflux
          
      enddo
      enddo
      enddo
      enddo
      enddo
      return
      end
      subroutine SECONDORDERGRADIENT_4D(
     &           grad_u
     &           ,igrad_ulo0,igrad_ulo1,igrad_ulo2,igrad_ulo3
     &           ,igrad_uhi0,igrad_uhi1,igrad_uhi2,igrad_uhi3
     &           ,ngrad_ucomp
     &           ,u
     &           ,iulo0,iulo1,iulo2,iulo3
     &           ,iuhi0,iuhi1,iuhi2,iuhi3
     &           ,nucomp
     &           ,iboxlo0,iboxlo1,iboxlo2,iboxlo3
     &           ,iboxhi0,iboxhi1,iboxhi2,iboxhi3
     &           ,dx
     &           ,facedir
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer ngrad_ucomp
      integer igrad_ulo0,igrad_ulo1,igrad_ulo2,igrad_ulo3
      integer igrad_uhi0,igrad_uhi1,igrad_uhi2,igrad_uhi3
      REAL_T grad_u(
     &           igrad_ulo0:igrad_uhi0,
     &           igrad_ulo1:igrad_uhi1,
     &           igrad_ulo2:igrad_uhi2,
     &           igrad_ulo3:igrad_uhi3,
     &           0:ngrad_ucomp-1)
      integer nucomp
      integer iulo0,iulo1,iulo2,iulo3
      integer iuhi0,iuhi1,iuhi2,iuhi3
      REAL_T u(
     &           iulo0:iuhi0,
     &           iulo1:iuhi1,
     &           iulo2:iuhi2,
     &           iulo3:iuhi3,
     &           0:nucomp-1)
      integer iboxlo0,iboxlo1,iboxlo2,iboxlo3
      integer iboxhi0,iboxhi1,iboxhi2,iboxhi3
      REAL_T dx(0:3)
      integer facedir
      integer i,j,k,l, n, nn, d
      integer in,jn,kn,ln, it,jt,kt,lt
      REAL_T avghi, avglo, factor(0:CH_SPACEDIM-1), result
      
      in = CHF_ID(0,facedir)
      jn = CHF_ID(1,facedir)
      kn = CHF_ID(2,facedir)
      ln = CHF_ID(3,facedir)
      do d=0, CH_SPACEDIM-1
        factor(d) = 1.0 / dx(d)
        if (facedir.ne.d) then
          factor(d) = factor(d) * 0.25
        endif
      enddo
      do n = 0, nucomp-1
        do d=0, CH_SPACEDIM-1
          nn = d + n * CH_SPACEDIM
          
      do l = iboxlo3,iboxhi3
      do k = iboxlo2,iboxhi2
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0

            if (facedir.eq.d) then
              result = ( u(i,j,k,l,n) - u(i-in,j-jn,k-kn,l-ln,n) )
     &               * factor(d)
           else
              
              it = CHF_ID(0,d)
              jt = CHF_ID(1,d)
              kt = CHF_ID(2,d)
              lt = CHF_ID(3,d)
              avghi = u(i+it,j+jt,k+kt,l+lt,n)
     &              + u(i+it-in,j+jt-jn,k+kt-kn,l+lt-ln,n)
              avglo = u(i-it,j-jt,k-kt,l-lt,n)
     &              + u(i-it-in,j-jt-jn,k-kt-kn,l-lt-ln,n)
              result = factor(d) * ( avghi - avglo );
            endif
            grad_u(i,j,k,l,nn) = result;
          
      enddo
      enddo
      enddo
      enddo
        enddo
      enddo
      return
      end
      subroutine SCALARFIELDMULTIPLY_4D(
     &           vector_field
     &           ,ivector_fieldlo0,ivector_fieldlo1,ivector_fieldlo2,ivector_fieldlo3
     &           ,ivector_fieldhi0,ivector_fieldhi1,ivector_fieldhi2,ivector_fieldhi3
     &           ,nvector_fieldcomp
     &           ,scalar_field
     &           ,iscalar_fieldlo0,iscalar_fieldlo1,iscalar_fieldlo2,iscalar_fieldlo3
     &           ,iscalar_fieldhi0,iscalar_fieldhi1,iscalar_fieldhi2,iscalar_fieldhi3
     &           ,nscalar_fieldcomp
     &           ,iboxlo0,iboxlo1,iboxlo2,iboxlo3
     &           ,iboxhi0,iboxhi1,iboxhi2,iboxhi3
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nvector_fieldcomp
      integer ivector_fieldlo0,ivector_fieldlo1,ivector_fieldlo2,ivector_fieldlo3
      integer ivector_fieldhi0,ivector_fieldhi1,ivector_fieldhi2,ivector_fieldhi3
      REAL_T vector_field(
     &           ivector_fieldlo0:ivector_fieldhi0,
     &           ivector_fieldlo1:ivector_fieldhi1,
     &           ivector_fieldlo2:ivector_fieldhi2,
     &           ivector_fieldlo3:ivector_fieldhi3,
     &           0:nvector_fieldcomp-1)
      integer nscalar_fieldcomp
      integer iscalar_fieldlo0,iscalar_fieldlo1,iscalar_fieldlo2,iscalar_fieldlo3
      integer iscalar_fieldhi0,iscalar_fieldhi1,iscalar_fieldhi2,iscalar_fieldhi3
      REAL_T scalar_field(
     &           iscalar_fieldlo0:iscalar_fieldhi0,
     &           iscalar_fieldlo1:iscalar_fieldhi1,
     &           iscalar_fieldlo2:iscalar_fieldhi2,
     &           iscalar_fieldlo3:iscalar_fieldhi3,
     &           0:nscalar_fieldcomp-1)
      integer iboxlo0,iboxlo1,iboxlo2,iboxlo3
      integer iboxhi0,iboxhi1,iboxhi2,iboxhi3
      integer i,j,k,l, n
      do n = 0, nvector_fieldcomp-1
        
      do l = iboxlo3,iboxhi3
      do k = iboxlo2,iboxhi2
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0

              vector_field(i,j,k,l,n) = vector_field(i,j,k,l,n)
     &                                      * scalar_field(i,j,k,l,0)
        
      enddo
      enddo
      enddo
      enddo
      enddo
      return
      end
