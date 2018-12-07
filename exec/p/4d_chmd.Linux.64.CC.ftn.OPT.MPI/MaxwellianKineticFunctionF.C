#include "REAL.H"
#include "SPACE.H"
#include "CONSTANTS.H"

#include "CONSTANTS.H"
      subroutine SET_MAXWELL4D_4D(
     &           f
     &           ,iflo0,iflo1,iflo2,iflo3
     &           ,ifhi0,ifhi1,ifhi2,ifhi3
     &           ,nfcomp
     &           ,iboxlo0,iboxlo1,iboxlo2,iboxlo3
     &           ,iboxhi0,iboxhi1,iboxhi2,iboxhi3
     &           ,coords
     &           ,icoordslo0,icoordslo1,icoordslo2,icoordslo3
     &           ,icoordshi0,icoordshi1,icoordshi2,icoordshi3
     &           ,ncoordscomp
     &           ,dens
     &           ,idenslo0,idenslo1,idenslo2,idenslo3
     &           ,idenshi0,idenshi1,idenshi2,idenshi3
     &           ,temp
     &           ,itemplo0,itemplo1,itemplo2,itemplo3
     &           ,itemphi0,itemphi1,itemphi2,itemphi3
     &           ,vshift
     &           ,ivshiftlo0,ivshiftlo1,ivshiftlo2,ivshiftlo3
     &           ,ivshifthi0,ivshifthi1,ivshifthi2,ivshifthi3
     &           ,b
     &           ,iblo0,iblo1,iblo2,iblo3
     &           ,ibhi0,ibhi1,ibhi2,ibhi3
     &           ,mass
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nfcomp
      integer iflo0,iflo1,iflo2,iflo3
      integer ifhi0,ifhi1,ifhi2,ifhi3
      REAL_T f(
     &           iflo0:ifhi0,
     &           iflo1:ifhi1,
     &           iflo2:ifhi2,
     &           iflo3:ifhi3,
     &           0:nfcomp-1)
      integer iboxlo0,iboxlo1,iboxlo2,iboxlo3
      integer iboxhi0,iboxhi1,iboxhi2,iboxhi3
      integer ncoordscomp
      integer icoordslo0,icoordslo1,icoordslo2,icoordslo3
      integer icoordshi0,icoordshi1,icoordshi2,icoordshi3
      REAL_T coords(
     &           icoordslo0:icoordshi0,
     &           icoordslo1:icoordshi1,
     &           icoordslo2:icoordshi2,
     &           icoordslo3:icoordshi3,
     &           0:ncoordscomp-1)
      integer idenslo0,idenslo1,idenslo2,idenslo3
      integer idenshi0,idenshi1,idenshi2,idenshi3
      REAL_T dens(
     &           idenslo0:idenshi0,
     &           idenslo1:idenshi1,
     &           idenslo2:idenshi2,
     &           idenslo3:idenshi3)
      integer itemplo0,itemplo1,itemplo2,itemplo3
      integer itemphi0,itemphi1,itemphi2,itemphi3
      REAL_T temp(
     &           itemplo0:itemphi0,
     &           itemplo1:itemphi1,
     &           itemplo2:itemphi2,
     &           itemplo3:itemphi3)
      integer ivshiftlo0,ivshiftlo1,ivshiftlo2,ivshiftlo3
      integer ivshifthi0,ivshifthi1,ivshifthi2,ivshifthi3
      REAL_T vshift(
     &           ivshiftlo0:ivshifthi0,
     &           ivshiftlo1:ivshifthi1,
     &           ivshiftlo2:ivshifthi2,
     &           ivshiftlo3:ivshifthi3)
      integer iblo0,iblo1,iblo2,iblo3
      integer ibhi0,ibhi1,ibhi2,ibhi3
      REAL_T b(
     &           iblo0:ibhi0,
     &           iblo1:ibhi1,
     &           iblo2:ibhi2,
     &           iblo3:ibhi3)
      REAL_T mass
      integer i,j,k,l,n
      integer ivp,imu
      REAL_T denloc, temploc, vshiftloc, bloc
      REAL_T vpar, mu
      REAL_T eparnorm, munorm
      REAL_T factor, val
      REAL_T minf, maxf
      minf=1.0d30
      maxf=zero
      factor = dsqrt(PI*(two/mass)**3)
      ivp = ncoordscomp-2
      imu = ncoordscomp-1
      do n=0,nfcomp-1
      
      do l = iboxlo3,iboxhi3
      do k = iboxlo2,iboxhi2
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0

#if CH_SPACEDIM==4
         vshiftloc = vshift(i,j,itemplo2,itemplo3)
         temploc   = temp(i,j,itemplo2,itemplo3)
         denloc    = dens(i,j,idenslo2,idenslo3)
         bloc      = b(i,j,iblo2,iblo3)
#else
         vshiftloc = vshift(i,j,k,itemplo3,itemplo4)
         temploc   = temp(i,j,k,itemplo3,itemplo4)
         denloc    = dens(i,j,k,idenslo3,idenslo4)
         bloc      = b(i,j,k,iblo3,iblo4)
#endif
         vpar     = coords(i,j,k,l,ivp)
         mu       = coords(i,j,k,l,imu)
         eparnorm = half * mass * (vpar-vshiftloc)**2 / temploc
         munorm   = half * bloc * mu / temploc
         val      = dexp( -( eparnorm + munorm ) )
         val    = val * denloc / ( factor * dsqrt( temploc) * temploc )
         minf = min(minf,val)
         maxf = max(maxf,val)
         f(i,j,k,l,n) = val
      
      enddo
      enddo
      enddo
      enddo
      enddo
      return
      end
      subroutine SET_CANONICAL_MAXWELL_4D(
     &           f
     &           ,iflo0,iflo1,iflo2,iflo3
     &           ,ifhi0,ifhi1,ifhi2,ifhi3
     &           ,nfcomp
     &           ,iboxlo0,iboxlo1,iboxlo2,iboxlo3
     &           ,iboxhi0,iboxhi1,iboxhi2,iboxhi3
     &           ,coords
     &           ,icoordslo0,icoordslo1,icoordslo2,icoordslo3
     &           ,icoordshi0,icoordshi1,icoordshi2,icoordshi3
     &           ,ncoordscomp
     &           ,toroidal_coords
     &           ,itoroidal_coordslo0,itoroidal_coordslo1,itoroidal_coordslo2,itoroidal_coordslo3
     &           ,itoroidal_coordshi0,itoroidal_coordshi1,itoroidal_coordshi2,itoroidal_coordshi3
     &           ,ntoroidal_coordscomp
     &           ,B_inj
     &           ,iB_injlo0,iB_injlo1,iB_injlo2,iB_injlo3
     &           ,iB_injhi0,iB_injhi1,iB_injhi2,iB_injhi3
     &           ,psi_inj
     &           ,ipsi_injlo0,ipsi_injlo1,ipsi_injlo2,ipsi_injlo3
     &           ,ipsi_injhi0,ipsi_injhi1,ipsi_injhi2,ipsi_injhi3
     &           ,RBtor
     &           ,psi_p
     &           ,dpsidr_p
     &           ,n_val
     &           ,n_kappa
     &           ,n_width
     &           ,T_val
     &           ,T_kappa
     &           ,T_width
     &           ,mode_coeff
     &           ,mass
     &           ,charge
     &           ,larmor
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nfcomp
      integer iflo0,iflo1,iflo2,iflo3
      integer ifhi0,ifhi1,ifhi2,ifhi3
      REAL_T f(
     &           iflo0:ifhi0,
     &           iflo1:ifhi1,
     &           iflo2:ifhi2,
     &           iflo3:ifhi3,
     &           0:nfcomp-1)
      integer iboxlo0,iboxlo1,iboxlo2,iboxlo3
      integer iboxhi0,iboxhi1,iboxhi2,iboxhi3
      integer ncoordscomp
      integer icoordslo0,icoordslo1,icoordslo2,icoordslo3
      integer icoordshi0,icoordshi1,icoordshi2,icoordshi3
      REAL_T coords(
     &           icoordslo0:icoordshi0,
     &           icoordslo1:icoordshi1,
     &           icoordslo2:icoordshi2,
     &           icoordslo3:icoordshi3,
     &           0:ncoordscomp-1)
      integer ntoroidal_coordscomp
      integer itoroidal_coordslo0,itoroidal_coordslo1,itoroidal_coordslo2,itoroidal_coordslo3
      integer itoroidal_coordshi0,itoroidal_coordshi1,itoroidal_coordshi2,itoroidal_coordshi3
      REAL_T toroidal_coords(
     &           itoroidal_coordslo0:itoroidal_coordshi0,
     &           itoroidal_coordslo1:itoroidal_coordshi1,
     &           itoroidal_coordslo2:itoroidal_coordshi2,
     &           itoroidal_coordslo3:itoroidal_coordshi3,
     &           0:ntoroidal_coordscomp-1)
      integer iB_injlo0,iB_injlo1,iB_injlo2,iB_injlo3
      integer iB_injhi0,iB_injhi1,iB_injhi2,iB_injhi3
      REAL_T B_inj(
     &           iB_injlo0:iB_injhi0,
     &           iB_injlo1:iB_injhi1,
     &           iB_injlo2:iB_injhi2,
     &           iB_injlo3:iB_injhi3)
      integer ipsi_injlo0,ipsi_injlo1,ipsi_injlo2,ipsi_injlo3
      integer ipsi_injhi0,ipsi_injhi1,ipsi_injhi2,ipsi_injhi3
      REAL_T psi_inj(
     &           ipsi_injlo0:ipsi_injhi0,
     &           ipsi_injlo1:ipsi_injhi1,
     &           ipsi_injlo2:ipsi_injhi2,
     &           ipsi_injlo3:ipsi_injhi3)
      REAL_T RBtor
      REAL_T psi_p
      REAL_T dpsidr_p
      REAL_T n_val
      REAL_T n_kappa
      REAL_T n_width
      REAL_T T_val
      REAL_T T_kappa
      REAL_T T_width
      REAL_T mode_coeff(0:3)
      REAL_T mass
      REAL_T charge
      REAL_T larmor
      integer i,j,k,l,n
      integer ivp,imu
      REAL_T den, temp, B, psi
      REAL_T psi_inv, vpar, mu, phi, theta
      REAL_T eparnorm, munorm
      REAL_T factor, norm, val, arg, pi
      pi = four*datan(one)
      norm      = larmor * mass / charge
      factor    = dsqrt(PI*(two/mass)**3)
      ivp = ncoordscomp-2
      imu = ncoordscomp-1
      do n=0,nfcomp-1
      
      do l = iboxlo3,iboxhi3
      do k = iboxlo2,iboxhi2
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0

#if CH_SPACEDIM==4
         B      = B_inj(i,j,iB_injlo2,iB_injlo3)
         psi    = psi_inj(i,j,ipsi_injlo2,ipsi_injlo3)
#else
         B      = B_inj(i,j,k,iB_injlo3,iB_injlo4)
         psi    = psi_inj(i,j,k,ipsi_injlo3,ipsi_injlo4)
         phi    = toroidal_coords(i,j,k,itoroidal_coordslo3,itoroidal_coordslo4,1)
         theta  = toroidal_coords(i,j,k,itoroidal_coordslo3,itoroidal_coordslo4,2)
#endif
         vpar     = coords(i,j,k,l,ivp)
         mu       = coords(i,j,k,l,imu)
         psi_inv = psi + norm*RBtor/B*vpar  
         den   = n_val * dexp(-n_kappa * n_width * dtanh((psi_inv-psi_p)/n_width))  
         temp  = T_val * dexp(-T_kappa * T_width * dtanh((psi_inv-psi_p)/T_width))  
         eparnorm = half * mass * vpar**2 / temp
         munorm   = half * B * mu / temp
         val      = dexp( -( eparnorm + munorm ) )
         val      = val * den / ( factor * dsqrt( temp) * temp )
         arg      = mode_coeff(1)*phi + mode_coeff(2)*theta
         val      = val * (one + mode_coeff(0)*cos(arg))
         f(i,j,k,l,n) = val
      
      enddo
      enddo
      enddo
      enddo
      enddo
      return
      end
