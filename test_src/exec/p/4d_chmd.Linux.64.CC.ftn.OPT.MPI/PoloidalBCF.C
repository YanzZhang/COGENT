#include "REAL.H"
#include "SPACE.H"
#include "CONSTANTS.H"

#include "CONSTANTS.H"
      subroutine SET_POLOIDAL_BC_4D(
     &           iboundaryboxlo0,iboundaryboxlo1,iboundaryboxlo2,iboundaryboxlo3
     &           ,iboundaryboxhi0,iboundaryboxhi1,iboundaryboxhi2,iboundaryboxhi3
     &           ,flipdist
     &           ,iflipdistlo0,iflipdistlo1,iflipdistlo2,iflipdistlo3
     &           ,iflipdisthi0,iflipdisthi1,iflipdisthi2,iflipdisthi3
     &           ,phi
     &           ,iphilo0,iphilo1,iphilo2,iphilo3
     &           ,iphihi0,iphihi1,iphihi2,iphihi3
     &           ,vp
     &           ,ivphi0
     &           ,vp_lo
     &           ,psi
     &           ,ipsihi0
     &           ,psilo
     &           ,psi_lim
     &           ,f
     &           ,iflo0,iflo1,iflo2,iflo3
     &           ,ifhi0,ifhi1,ifhi2,ifhi3
     &           ,charge
     &           ,ilohi
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iboundaryboxlo0,iboundaryboxlo1,iboundaryboxlo2,iboundaryboxlo3
      integer iboundaryboxhi0,iboundaryboxhi1,iboundaryboxhi2,iboundaryboxhi3
      integer iflipdistlo0,iflipdistlo1,iflipdistlo2,iflipdistlo3
      integer iflipdisthi0,iflipdisthi1,iflipdisthi2,iflipdisthi3
      REAL_T flipdist(
     &           iflipdistlo0:iflipdisthi0,
     &           iflipdistlo1:iflipdisthi1,
     &           iflipdistlo2:iflipdisthi2,
     &           iflipdistlo3:iflipdisthi3)
      integer iphilo0,iphilo1,iphilo2,iphilo3
      integer iphihi0,iphihi1,iphihi2,iphihi3
      REAL_T phi(
     &           iphilo0:iphihi0,
     &           iphilo1:iphihi1,
     &           iphilo2:iphihi2,
     &           iphilo3:iphihi3)
      integer ivphi0
      REAL_T vp(
     &           0:ivphi0)
      integer vp_lo
      integer ipsihi0
      REAL_T psi(
     &           0:ipsihi0)
      integer psilo
      REAL_T psi_lim
      integer iflo0,iflo1,iflo2,iflo3
      integer ifhi0,ifhi1,ifhi2,ifhi3
      REAL_T f(
     &           iflo0:ifhi0,
     &           iflo1:ifhi1,
     &           iflo2:ifhi2,
     &           iflo3:ifhi3)
      REAL_T charge
      integer ilohi
      double precision vploc,vpsq,df,phibndy
      integer i,j,k,l
      integer jbndy,jbndy1,vp_hi,nvp,twojbndy1,jsrc,jbndy2,jbndygtr,jbndyless,iparsign
      if (ilohi == 1) then
           jbndy = iboundaryboxhi1
           jbndy1 = jbndy+1
           jbndy2 = jbndy+2
           jbndygtr = jbndy2
           jbndyless = jbndy1
           iparsign = 1
      else
           jbndy = iboundaryboxlo1
           jbndy1 = jbndy-1
           jbndy2 = jbndy-2
           jbndygtr = jbndy1
           jbndyless = jbndy2
           iparsign = -1
      endif
      twojbndy1 = jbndy+jbndy1
      
      do l = iboundaryboxlo3,iboundaryboxhi3
      do k = iboundaryboxlo2,iboundaryboxhi2
      do j = iboundaryboxlo1,iboundaryboxhi1
      do i = iboundaryboxlo0,iboundaryboxhi0

#if CH_SPACEDIM==5
         phibndy = phi(i,jbndy,k,iphilo3,iphilo4)
#else
         phibndy = phi(i,jbndy,iphilo2,iphilo3)
#endif
         if (psi(i-psilo) .ge. psi_lim) then
#if CH_SPACEDIM==5
            df = f(i,jbndyless,k,l,m)-f(i,jbndygtr,k,l,m)
#else
            df = f(i,jbndyless,k,l)-f(i,jbndygtr,k,l)
#endif
            vploc = iparsign*vp(k-vp_lo)
            vpsq=vploc**2
            if (1 .lt. 0) then
#if CH_SPACEDIM==5
                f(i,j,k,l,m) = f(i,jbndy1,k,l,m)+(jbndy1-j)*df
#else
                f(i,j,k,l) = f(i,jbndy1,k,l)+(jbndy1-j)*df
#endif
            else
                jsrc = twojbndy1-j
                vpsq = vploc**2
                if (vpsq .gt.  (-charge*phibndy)) then
#if CH_SPACEDIM==5
                     f(i,j,k,l,m) = 0.
#else
                     f(i,j,k,l) = 0.
#endif
                else
#if CH_SPACEDIM==5
                    f(i,j,k,l,m) = flipdist(i,jsrc,k,-l,m)
#else
                    f(i,j,k,l) = flipdist(i,jsrc,-k,l)
#endif
                endif
            endif
         endif
      
      enddo
      enddo
      enddo
      enddo
      return
      end
