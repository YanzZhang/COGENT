#include "REAL.H"
#include "SPACE.H"
#include "CONSTANTS.H"

#include "CONSTANTS.H"
      subroutine COMPUTE_GK_VELOCITY_4D(
     &           dir
     &           ,igridboxlo0,igridboxlo1,igridboxlo2,igridboxlo3
     &           ,igridboxhi0,igridboxhi1,igridboxhi2,igridboxhi3
     &           ,h
     &           ,Z
     &           ,mass
     &           ,larmor
     &           ,include_drifts
     &           ,include_par_streaming
     &           ,include_gradb
     &           ,mag_drifts_only
     &           ,use_field_alignment
     &           ,E
     &           ,iElo0,iElo1,iElo2,iElo3
     &           ,iEhi0,iEhi1,iEhi2,iEhi3
     &           ,nEcomp
     &           ,Bvec
     &           ,iBveclo0,iBveclo1,iBveclo2,iBveclo3
     &           ,iBvechi0,iBvechi1,iBvechi2,iBvechi3
     &           ,nBveccomp
     &           ,gradB
     &           ,igradBlo0,igradBlo1,igradBlo2,igradBlo3
     &           ,igradBhi0,igradBhi1,igradBhi2,igradBhi3
     &           ,ngradBcomp
     &           ,curlb
     &           ,icurlblo0,icurlblo1,icurlblo2,icurlblo3
     &           ,icurlbhi0,icurlbhi1,icurlbhi2,icurlbhi3
     &           ,ncurlbcomp
     &           ,velocity
     &           ,ivelocitylo0,ivelocitylo1,ivelocitylo2,ivelocitylo3
     &           ,ivelocityhi0,ivelocityhi1,ivelocityhi2,ivelocityhi3
     &           ,nvelocitycomp
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer dir
      integer igridboxlo0,igridboxlo1,igridboxlo2,igridboxlo3
      integer igridboxhi0,igridboxhi1,igridboxhi2,igridboxhi3
      REAL_T h(0:3)
      REAL_T Z
      REAL_T mass
      REAL_T larmor
      integer include_drifts
      integer include_par_streaming
      integer include_gradb
      integer mag_drifts_only
      integer use_field_alignment
      integer nEcomp
      integer iElo0,iElo1,iElo2,iElo3
      integer iEhi0,iEhi1,iEhi2,iEhi3
      REAL_T E(
     &           iElo0:iEhi0,
     &           iElo1:iEhi1,
     &           iElo2:iEhi2,
     &           iElo3:iEhi3,
     &           0:nEcomp-1)
      integer nBveccomp
      integer iBveclo0,iBveclo1,iBveclo2,iBveclo3
      integer iBvechi0,iBvechi1,iBvechi2,iBvechi3
      REAL_T Bvec(
     &           iBveclo0:iBvechi0,
     &           iBveclo1:iBvechi1,
     &           iBveclo2:iBvechi2,
     &           iBveclo3:iBvechi3,
     &           0:nBveccomp-1)
      integer ngradBcomp
      integer igradBlo0,igradBlo1,igradBlo2,igradBlo3
      integer igradBhi0,igradBhi1,igradBhi2,igradBhi3
      REAL_T gradB(
     &           igradBlo0:igradBhi0,
     &           igradBlo1:igradBhi1,
     &           igradBlo2:igradBhi2,
     &           igradBlo3:igradBhi3,
     &           0:ngradBcomp-1)
      integer ncurlbcomp
      integer icurlblo0,icurlblo1,icurlblo2,icurlblo3
      integer icurlbhi0,icurlbhi1,icurlbhi2,icurlbhi3
      REAL_T curlb(
     &           icurlblo0:icurlbhi0,
     &           icurlblo1:icurlbhi1,
     &           icurlblo2:icurlbhi2,
     &           icurlblo3:icurlbhi3,
     &           0:ncurlbcomp-1)
      integer nvelocitycomp
      integer ivelocitylo0,ivelocitylo1,ivelocitylo2,ivelocitylo3
      integer ivelocityhi0,ivelocityhi1,ivelocityhi2,ivelocityhi3
      REAL_T velocity(
     &           ivelocitylo0:ivelocityhi0,
     &           ivelocitylo1:ivelocityhi1,
     &           ivelocitylo2:ivelocityhi2,
     &           ivelocitylo3:ivelocityhi3,
     &           0:nvelocitycomp-1)
      double precision vpar, mu, Bmag, Bstar_par, G(0:2), b(0:2), Bstar(0:2),
     &     bxG(0:2), mag_drift_fac(0:2), Bstar_dot_G, fac, this_Bvec(0:2),
     &     this_curlb(0:2), this_gradB(0:2), this_E(0:2)
      integer i0,i1,i2,i3, n, phys_comp(0:2)
      if (include_drifts .ne. 0) then
         fac = larmor / Z
      else
         fac = zero
      endif
      phys_comp(0) = 0
      phys_comp(1) = 6 - CH_SPACEDIM
      phys_comp(2) = 2
      
      do i3 = igridboxlo3,igridboxhi3
      do i2 = igridboxlo2,igridboxhi2
      do i1 = igridboxlo1,igridboxhi1
      do i0 = igridboxlo0,igridboxhi0

#if CH_SPACEDIM==4
         vpar  = ( i2 + half*(1-CHF_ID(CH_SPACEDIM-2,dir)) ) * h(CH_SPACEDIM-2)
         mu    = ( i3 + half*(1-CHF_ID(CH_SPACEDIM-1,dir)) ) * h(CH_SPACEDIM-1)
#endif
#if CH_SPACEDIM==5
         vpar  = ( i3 + half*(1-CHF_ID(CH_SPACEDIM-2,dir)) ) * h(CH_SPACEDIM-2)
         mu    = ( i4 + half*(1-CHF_ID(CH_SPACEDIM-1,dir)) ) * h(CH_SPACEDIM-1)
#endif
         Bmag = zero
         do n = 0, 2
#if CH_SPACEDIM==4
            this_Bvec(n)  =  Bvec(i0,i1,iBveclo2 ,iBveclo3 ,n)
            this_curlb(n) = curlb(i0,i1,icurlblo2,icurlblo3,n)
            this_gradB(n) = gradB(i0,i1,igradBlo2,igradBlo3,n)
            this_E(n)     =     E(i0,i1,iElo2    ,iElo3    ,n)
#endif
#if CH_SPACEDIM==5
            this_Bvec(n)  =  Bvec(i0,i1,i2,iBveclo3 ,iBveclo4 ,n)
            this_curlb(n) = curlb(i0,i1,i2,icurlblo3,icurlblo4,n)
            this_gradB(n) = gradB(i0,i1,i2,igradBlo3,igradBlo4,n)
            this_E(n)     =     E(i0,i1,i2,iElo3    ,iElo4    ,n)
#endif
            Bmag = Bmag + this_Bvec(n)**2
         enddo
         Bmag = dsqrt(Bmag)
         Bstar_par = zero
         do n = 0, 2
            G(n) = -Z * this_E(n)
            if (mag_drifts_only.ne.0) then
               G(n) = 0.0
            endif 
            b(n) = this_Bvec(n) / Bmag
            mag_drift_fac(n) = fac * mass * vpar * this_curlb(n)
            Bstar(n) = this_Bvec(n) + mag_drift_fac(n)
            Bstar_par = Bstar_par + b(n) * Bstar(n)
         enddo
         if (include_gradb .ne. 0) then
            do n = 0, 2
               G(n) = G(n) + half * mu * this_gradB(n)
            enddo
         endif
         bxG(0) = b(1)*G(2) - b(2)*G(1)
         bxG(1) = b(2)*G(0) - b(0)*G(2)
         bxG(2) = b(0)*G(1) - b(1)*G(0)
         do n = 0, CH_SPACEDIM-3
            velocity(i0,i1,i2,i3,n)  = fac*bxG(phys_comp(n)) + vpar * mag_drift_fac(phys_comp(n))
            if (include_par_streaming .ne. 0) then
               velocity(i0,i1,i2,i3,n) = velocity(i0,i1,i2,i3,n) 
     &              + vpar * this_Bvec(phys_comp(n))
            endif
            if ((use_field_alignment .ne. 0). and. (dir.eq.0)) then
                velocity(i0,i1,i2,i3,n)  = fac * bxG(phys_comp(n)) + vpar * mag_drift_fac(phys_comp(n))
            endif
            velocity(i0,i1,i2,i3,n)  = velocity(i0,i1,i2,i3,n) / Bstar_par
         enddo
         Bstar_dot_G = zero
         do n = 0, 2
            Bstar_dot_G = Bstar_dot_G + Bstar(n) * G(n)
         enddo
         velocity(i0,i1,i2,i3,CH_SPACEDIM-2) = - Bstar_dot_G / (mass * Bstar_par)
         velocity(i0,i1,i2,i3,CH_SPACEDIM-1) = zero
      
      enddo
      enddo
      enddo
      enddo
      return
      end
      subroutine COMPUTE_BFIELD_VELOCITY_4D(
     &           igridboxlo0,igridboxlo1,igridboxlo2,igridboxlo3
     &           ,igridboxhi0,igridboxhi1,igridboxhi2,igridboxhi3
     &           ,Bvec
     &           ,iBveclo0,iBveclo1,iBveclo2,iBveclo3
     &           ,iBvechi0,iBvechi1,iBvechi2,iBvechi3
     &           ,nBveccomp
     &           ,velocity
     &           ,ivelocitylo0,ivelocitylo1,ivelocitylo2,ivelocitylo3
     &           ,ivelocityhi0,ivelocityhi1,ivelocityhi2,ivelocityhi3
     &           ,nvelocitycomp
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer igridboxlo0,igridboxlo1,igridboxlo2,igridboxlo3
      integer igridboxhi0,igridboxhi1,igridboxhi2,igridboxhi3
      integer nBveccomp
      integer iBveclo0,iBveclo1,iBveclo2,iBveclo3
      integer iBvechi0,iBvechi1,iBvechi2,iBvechi3
      REAL_T Bvec(
     &           iBveclo0:iBvechi0,
     &           iBveclo1:iBvechi1,
     &           iBveclo2:iBvechi2,
     &           iBveclo3:iBvechi3,
     &           0:nBveccomp-1)
      integer nvelocitycomp
      integer ivelocitylo0,ivelocitylo1,ivelocitylo2,ivelocitylo3
      integer ivelocityhi0,ivelocityhi1,ivelocityhi2,ivelocityhi3
      REAL_T velocity(
     &           ivelocitylo0:ivelocityhi0,
     &           ivelocitylo1:ivelocityhi1,
     &           ivelocitylo2:ivelocityhi2,
     &           ivelocitylo3:ivelocityhi3,
     &           0:nvelocitycomp-1)
      integer i0,i1,i2,i3, kk, itor4vel, phys_comp(0:2), dimless1, dimless2, dimless3, n
      kk = iBveclo2
      phys_comp(0) = 0
      phys_comp(1) = 6-CH_SPACEDIM
      phys_comp(2) = 2
      dimless1 = CH_SPACEDIM-1
      dimless2 = CH_SPACEDIM-2
      dimless3 = CH_SPACEDIM-3
      
      do i3 = igridboxlo3,igridboxhi3
      do i2 = igridboxlo2,igridboxhi2
      do i1 = igridboxlo1,igridboxhi1
      do i0 = igridboxlo0,igridboxhi0

#if CH_SPACEDIM==5
         kk = i2
#endif
        do n = 0, dimless3
           velocity(i0,i1,i2,i3,n) = Bvec(i0,i1,kk,iBveclo3,phys_comp(n))
        enddo
        do n=dimless2, dimless1
           velocity(i0,i1,i2,i3,n) = zero
        enddo
      
      enddo
      enddo
      enddo
      enddo
      return
      end
      subroutine ANNULUS_POLVEL_TEST_4D(
     &           dir
     &           ,igridboxlo0,igridboxlo1,igridboxlo2,igridboxlo3
     &           ,igridboxhi0,igridboxhi1,igridboxhi2,igridboxhi3
     &           ,h
     &           ,rmin
     &           ,rbar
     &           ,R0
     &           ,velocity
     &           ,ivelocitylo0,ivelocitylo1,ivelocitylo2,ivelocitylo3
     &           ,ivelocityhi0,ivelocityhi1,ivelocityhi2,ivelocityhi3
     &           ,nvelocitycomp
     &           ,l_const_minorrad
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer dir
      integer igridboxlo0,igridboxlo1,igridboxlo2,igridboxlo3
      integer igridboxhi0,igridboxhi1,igridboxhi2,igridboxhi3
      REAL_T h(0:3)
      REAL_T rmin
      REAL_T rbar
      REAL_T R0
      integer nvelocitycomp
      integer ivelocitylo0,ivelocitylo1,ivelocitylo2,ivelocitylo3
      integer ivelocityhi0,ivelocityhi1,ivelocityhi2,ivelocityhi3
      REAL_T velocity(
     &           ivelocitylo0:ivelocityhi0,
     &           ivelocitylo1:ivelocityhi1,
     &           ivelocitylo2:ivelocityhi2,
     &           ivelocitylo3:ivelocityhi3,
     &           0:nvelocitycomp-1)
      integer l_const_minorrad
      double precision r, theta,Rmaj,costheta,radfac
      integer i0,i1,i2,i3, m
      
      do i3 = igridboxlo3,igridboxhi3
      do i2 = igridboxlo2,igridboxhi2
      do i1 = igridboxlo1,igridboxhi1
      do i0 = igridboxlo0,igridboxhi0

         theta = ( i1 + half*(1-CHF_ID(1,dir)) )*h(1)
         costheta = dcos(theta)
         if (l_const_minorrad .eq. 1) then
             r = rbar
             Rmaj = R0
         else
             r     = ( i0 + half*(1-CHF_ID(0,dir)) )*h(0) + rmin
            Rmaj = R0 + r*costheta
         endif
         radfac = r*R0/(rbar*Rmaj)
         velocity(i0,i1,i2,i3,0) = -radfac * dsin(theta)
         velocity(i0,i1,i2,i3,1) =  radfac * costheta
         velocity(i0,i1,i2,i3,2) = zero
         velocity(i0,i1,i2,i3,3) = zero
#if CH_SPACEDIM==5
         velocity(i0,i1,i2,i3,4) = zero
#endif
      
      enddo
      enddo
      enddo
      enddo
      return
      end
      subroutine ANNULUS_RADVEL_TEST_4D(
     &           dir
     &           ,igridboxlo0,igridboxlo1,igridboxlo2,igridboxlo3
     &           ,igridboxhi0,igridboxhi1,igridboxhi2,igridboxhi3
     &           ,h
     &           ,rmin
     &           ,rbar
     &           ,R0
     &           ,velocity
     &           ,ivelocitylo0,ivelocitylo1,ivelocitylo2,ivelocitylo3
     &           ,ivelocityhi0,ivelocityhi1,ivelocityhi2,ivelocityhi3
     &           ,nvelocitycomp
     &           ,l_const_minorrad
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer dir
      integer igridboxlo0,igridboxlo1,igridboxlo2,igridboxlo3
      integer igridboxhi0,igridboxhi1,igridboxhi2,igridboxhi3
      REAL_T h(0:3)
      REAL_T rmin
      REAL_T rbar
      REAL_T R0
      integer nvelocitycomp
      integer ivelocitylo0,ivelocitylo1,ivelocitylo2,ivelocitylo3
      integer ivelocityhi0,ivelocityhi1,ivelocityhi2,ivelocityhi3
      REAL_T velocity(
     &           ivelocitylo0:ivelocityhi0,
     &           ivelocitylo1:ivelocityhi1,
     &           ivelocitylo2:ivelocityhi2,
     &           ivelocitylo3:ivelocityhi3,
     &           0:nvelocitycomp-1)
      integer l_const_minorrad
      double precision r, theta,Rmaj,costheta,radfac
      integer i0,i1,i2,i3, m
      
      do i3 = igridboxlo3,igridboxhi3
      do i2 = igridboxlo2,igridboxhi2
      do i1 = igridboxlo1,igridboxhi1
      do i0 = igridboxlo0,igridboxhi0

         theta = ( i1 + half*(1-CHF_ID(1,dir)) )*h(1)
         costheta = dcos(theta)
         if (l_const_minorrad .eq. 1) then
             r = rbar
             Rmaj = R0
         else
             r     = ( i0 + half*(1-CHF_ID(0,dir)) )*h(0) + rmin
            Rmaj = R0 + r*costheta
         endif
         radfac = (rbar*R0)/r*Rmaj
         velocity(i0,i1,i2,i3,0) = costheta*radfac
         velocity(i0,i1,i2,i3,1) = dsin(theta)*radfac
         velocity(i0,i1,i2,i3,2) = zero
         velocity(i0,i1,i2,i3,3) = zero
#if CH_SPACEDIM==5
         velocity(i0,i1,i2,i3,4) = zero
#endif
      
      enddo
      enddo
      enddo
      enddo
      return
      end
      subroutine ANNULUS_RADPOLVEL_TEST_4D(
     &           dir
     &           ,igridboxlo0,igridboxlo1,igridboxlo2,igridboxlo3
     &           ,igridboxhi0,igridboxhi1,igridboxhi2,igridboxhi3
     &           ,h
     &           ,rmin
     &           ,rbar
     &           ,R0
     &           ,velocity
     &           ,ivelocitylo0,ivelocitylo1,ivelocitylo2,ivelocitylo3
     &           ,ivelocityhi0,ivelocityhi1,ivelocityhi2,ivelocityhi3
     &           ,nvelocitycomp
     &           ,l_const_minorrad
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer dir
      integer igridboxlo0,igridboxlo1,igridboxlo2,igridboxlo3
      integer igridboxhi0,igridboxhi1,igridboxhi2,igridboxhi3
      REAL_T h(0:3)
      REAL_T rmin
      REAL_T rbar
      REAL_T R0
      integer nvelocitycomp
      integer ivelocitylo0,ivelocitylo1,ivelocitylo2,ivelocitylo3
      integer ivelocityhi0,ivelocityhi1,ivelocityhi2,ivelocityhi3
      REAL_T velocity(
     &           ivelocitylo0:ivelocityhi0,
     &           ivelocitylo1:ivelocityhi1,
     &           ivelocitylo2:ivelocityhi2,
     &           ivelocitylo3:ivelocityhi3,
     &           0:nvelocitycomp-1)
      integer l_const_minorrad
      double precision r, theta,delta_r,polvelmult,radfac,sintheta,costheta
      double precision Rmaj,radfac1
      integer i0,i1,i2,i3, m, imax
      imax = igridboxhi0
      delta_r = ( imax+half*(1-CHF_ID(0,dir)))*h(0)
      polvelmult = 2.*3.14159/delta_r
      
      do i3 = igridboxlo3,igridboxhi3
      do i2 = igridboxlo2,igridboxhi2
      do i1 = igridboxlo1,igridboxhi1
      do i0 = igridboxlo0,igridboxhi0

         r     = ( i0 + half*(1-CHF_ID(0,dir)) )*h(0) + rmin
         theta = ( i1 + half*(1-CHF_ID(1,dir)) )*h(1)
         costheta = dcos(theta)
         sintheta = dsin(theta)
         if (l_const_minorrad .eq. 1) then
             r = rbar
             Rmaj = R0
         else
             r     = ( i0 + half*(1-CHF_ID(0,dir)) )*h(0) + rmin
            Rmaj = R0 + r*costheta
         endif
         radfac = (rbar*R0)/r*Rmaj
         radfac1 = r*R0/(rbar*Rmaj)
         velocity(i0,i1,i2,i3,0) = -radfac1*polvelmult*sintheta + costheta*radfac
         velocity(i0,i1,i2,i3,1) =  radfac1*polvelmult*costheta + sintheta*radfac
         velocity(i0,i1,i2,i3,2) = zero
         velocity(i0,i1,i2,i3,3) = zero
#if CH_SPACEDIM==5
         velocity(i0,i1,i2,i3,4) = zero
#endif
      
      enddo
      enddo
      enddo
      enddo
      return
      end
      subroutine MAJOR_RADIAL_VEL_TEST_4D(
     &           dir
     &           ,igridboxlo0,igridboxlo1,igridboxlo2,igridboxlo3
     &           ,igridboxhi0,igridboxhi1,igridboxhi2,igridboxhi3
     &           ,velocity
     &           ,ivelocitylo0,ivelocitylo1,ivelocitylo2,ivelocitylo3
     &           ,ivelocityhi0,ivelocityhi1,ivelocityhi2,ivelocityhi3
     &           ,nvelocitycomp
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer dir
      integer igridboxlo0,igridboxlo1,igridboxlo2,igridboxlo3
      integer igridboxhi0,igridboxhi1,igridboxhi2,igridboxhi3
      integer nvelocitycomp
      integer ivelocitylo0,ivelocitylo1,ivelocitylo2,ivelocitylo3
      integer ivelocityhi0,ivelocityhi1,ivelocityhi2,ivelocityhi3
      REAL_T velocity(
     &           ivelocitylo0:ivelocityhi0,
     &           ivelocitylo1:ivelocityhi1,
     &           ivelocitylo2:ivelocityhi2,
     &           ivelocitylo3:ivelocityhi3,
     &           0:nvelocitycomp-1)
      integer i0,i1,i2,i3
      
      do i3 = igridboxlo3,igridboxhi3
      do i2 = igridboxlo2,igridboxhi2
      do i1 = igridboxlo1,igridboxhi1
      do i0 = igridboxlo0,igridboxhi0

         velocity(i0,i1,i2,i3,0) = one
         velocity(i0,i1,i2,i3,1) = zero
         velocity(i0,i1,i2,i3,2) = zero
         velocity(i0,i1,i2,i3,3) = zero
#if CH_SPACEDIM==5
         velocity(i0,i1,i2,i3,4) = zero
#endif
      
      enddo
      enddo
      enddo
      enddo
      return
      end
      subroutine FREE_STREAM_VEL_TEST_4D(
     &           igridboxlo0,igridboxlo1,igridboxlo2,igridboxlo3
     &           ,igridboxhi0,igridboxhi1,igridboxhi2,igridboxhi3
     &           ,component
     &           ,major_radius
     &           ,imajor_radiuslo0,imajor_radiuslo1,imajor_radiuslo2,imajor_radiuslo3
     &           ,imajor_radiushi0,imajor_radiushi1,imajor_radiushi2,imajor_radiushi3
     &           ,axisymmetric
     &           ,velocity
     &           ,ivelocitylo0,ivelocitylo1,ivelocitylo2,ivelocitylo3
     &           ,ivelocityhi0,ivelocityhi1,ivelocityhi2,ivelocityhi3
     &           ,nvelocitycomp
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer igridboxlo0,igridboxlo1,igridboxlo2,igridboxlo3
      integer igridboxhi0,igridboxhi1,igridboxhi2,igridboxhi3
      REAL_T component(0:3)
      integer imajor_radiuslo0,imajor_radiuslo1,imajor_radiuslo2,imajor_radiuslo3
      integer imajor_radiushi0,imajor_radiushi1,imajor_radiushi2,imajor_radiushi3
      REAL_T major_radius(
     &           imajor_radiuslo0:imajor_radiushi0,
     &           imajor_radiuslo1:imajor_radiushi1,
     &           imajor_radiuslo2:imajor_radiushi2,
     &           imajor_radiuslo3:imajor_radiushi3)
      integer axisymmetric
      integer nvelocitycomp
      integer ivelocitylo0,ivelocitylo1,ivelocitylo2,ivelocitylo3
      integer ivelocityhi0,ivelocityhi1,ivelocityhi2,ivelocityhi3
      REAL_T velocity(
     &           ivelocitylo0:ivelocityhi0,
     &           ivelocitylo1:ivelocityhi1,
     &           ivelocitylo2:ivelocityhi2,
     &           ivelocitylo3:ivelocityhi3,
     &           0:nvelocitycomp-1)
      integer i0,i1,i2,i3,n,kk
      kk = imajor_radiuslo2
      
      do i3 = igridboxlo3,igridboxhi3
      do i2 = igridboxlo2,igridboxhi2
      do i1 = igridboxlo1,igridboxhi1
      do i0 = igridboxlo0,igridboxhi0

         do n = 0, nvelocitycomp-1
           velocity(i0,i1,i2,i3,n) = component(n)
         enddo
      
      enddo
      enddo
      enddo
      enddo
      if (axisymmetric .ne. 0) then
         
      do i3 = igridboxlo3,igridboxhi3
      do i2 = igridboxlo2,igridboxhi2
      do i1 = igridboxlo1,igridboxhi1
      do i0 = igridboxlo0,igridboxhi0

#if CH_SPACEDIM==5
            kk = i2
#endif
            velocity(i0,i1,i2,i3,0) = velocity(i0,i1,i2,i3,0)
     &         / major_radius(i0,i1,kk,imajor_radiuslo3)
         
      enddo
      enddo
      enddo
      enddo
      endif
      return
      end
      subroutine COMPUTE_BSTAR_4D(
     &           igridboxlo0,igridboxlo1,igridboxlo2,igridboxlo3
     &           ,igridboxhi0,igridboxhi1,igridboxhi2,igridboxhi3
     &           ,h
     &           ,prefactor
     &           ,B
     &           ,iBlo0,iBlo1,iBlo2,iBlo3
     &           ,iBhi0,iBhi1,iBhi2,iBhi3
     &           ,nBcomp
     &           ,B_magnitude
     &           ,iB_magnitudelo0,iB_magnitudelo1,iB_magnitudelo2,iB_magnitudelo3
     &           ,iB_magnitudehi0,iB_magnitudehi1,iB_magnitudehi2,iB_magnitudehi3
     &           ,curlb
     &           ,icurlblo0,icurlblo1,icurlblo2,icurlblo3
     &           ,icurlbhi0,icurlbhi1,icurlbhi2,icurlbhi3
     &           ,ncurlbcomp
     &           ,bdotcurlb
     &           ,ibdotcurlblo0,ibdotcurlblo1,ibdotcurlblo2,ibdotcurlblo3
     &           ,ibdotcurlbhi0,ibdotcurlbhi1,ibdotcurlbhi2,ibdotcurlbhi3
     &           ,BStar
     &           ,iBStarlo0,iBStarlo1,iBStarlo2,iBStarlo3
     &           ,iBStarhi0,iBStarhi1,iBStarhi2,iBStarhi3
     &           ,nBStarcomp
     &           ,BStarPar
     &           ,iBStarParlo0,iBStarParlo1,iBStarParlo2,iBStarParlo3
     &           ,iBStarParhi0,iBStarParhi1,iBStarParhi2,iBStarParhi3
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer igridboxlo0,igridboxlo1,igridboxlo2,igridboxlo3
      integer igridboxhi0,igridboxhi1,igridboxhi2,igridboxhi3
      REAL_T h(0:3)
      REAL_T prefactor
      integer nBcomp
      integer iBlo0,iBlo1,iBlo2,iBlo3
      integer iBhi0,iBhi1,iBhi2,iBhi3
      REAL_T B(
     &           iBlo0:iBhi0,
     &           iBlo1:iBhi1,
     &           iBlo2:iBhi2,
     &           iBlo3:iBhi3,
     &           0:nBcomp-1)
      integer iB_magnitudelo0,iB_magnitudelo1,iB_magnitudelo2,iB_magnitudelo3
      integer iB_magnitudehi0,iB_magnitudehi1,iB_magnitudehi2,iB_magnitudehi3
      REAL_T B_magnitude(
     &           iB_magnitudelo0:iB_magnitudehi0,
     &           iB_magnitudelo1:iB_magnitudehi1,
     &           iB_magnitudelo2:iB_magnitudehi2,
     &           iB_magnitudelo3:iB_magnitudehi3)
      integer ncurlbcomp
      integer icurlblo0,icurlblo1,icurlblo2,icurlblo3
      integer icurlbhi0,icurlbhi1,icurlbhi2,icurlbhi3
      REAL_T curlb(
     &           icurlblo0:icurlbhi0,
     &           icurlblo1:icurlbhi1,
     &           icurlblo2:icurlbhi2,
     &           icurlblo3:icurlbhi3,
     &           0:ncurlbcomp-1)
      integer ibdotcurlblo0,ibdotcurlblo1,ibdotcurlblo2,ibdotcurlblo3
      integer ibdotcurlbhi0,ibdotcurlbhi1,ibdotcurlbhi2,ibdotcurlbhi3
      REAL_T bdotcurlb(
     &           ibdotcurlblo0:ibdotcurlbhi0,
     &           ibdotcurlblo1:ibdotcurlbhi1,
     &           ibdotcurlblo2:ibdotcurlbhi2,
     &           ibdotcurlblo3:ibdotcurlbhi3)
      integer nBStarcomp
      integer iBStarlo0,iBStarlo1,iBStarlo2,iBStarlo3
      integer iBStarhi0,iBStarhi1,iBStarhi2,iBStarhi3
      REAL_T BStar(
     &           iBStarlo0:iBStarhi0,
     &           iBStarlo1:iBStarhi1,
     &           iBStarlo2:iBStarhi2,
     &           iBStarlo3:iBStarhi3,
     &           0:nBStarcomp-1)
      integer iBStarParlo0,iBStarParlo1,iBStarParlo2,iBStarParlo3
      integer iBStarParhi0,iBStarParhi1,iBStarParhi2,iBStarParhi3
      REAL_T BStarPar(
     &           iBStarParlo0:iBStarParhi0,
     &           iBStarParlo1:iBStarParhi1,
     &           iBStarParlo2:iBStarParhi2,
     &           iBStarParlo3:iBStarParhi3)
      integer i0,i1,i2,i3, m, kk, indvpar
      double precision v_par_avg, dv_parallel
      kk = iBlo2
      dv_parallel = h(2)
      
      do i3 = igridboxlo3,igridboxhi3
      do i2 = igridboxlo2,igridboxhi2
      do i1 = igridboxlo1,igridboxhi1
      do i0 = igridboxlo0,igridboxhi0

         indvpar = i2
#if CH_SPACEDIM==5
         kk = i2
         indvpar = i3
#endif
         v_par_avg = (indvpar + half)*dv_parallel
         do m = 0, 2
            BStar(i0,i1,i2,i3,m) = B(i0,i1,kk,iBlo3,m)
     &           + prefactor * v_par_avg * curlb(i0,i1,kk,icurlblo3,m)
         enddo
         BStarPar(i0,i1,i2,i3) = B_magnitude(i0,i1,kk,iB_magnitudelo3)
     &           + prefactor * v_par_avg * bdotcurlb(i0,i1,kk,ibdotcurlblo3)
      
      enddo
      enddo
      enddo
      enddo
      return
      end
      subroutine PSCSPOINTDOTPRODCFG_4D(
     &           c
     &           ,iclo0,iclo1,iclo2,iclo3
     &           ,ichi0,ichi1,ichi2,ichi3
     &           ,a
     &           ,ialo0,ialo1,ialo2,ialo3
     &           ,iahi0,iahi1,iahi2,iahi3
     &           ,nacomp
     &           ,astart
     &           ,b
     &           ,iblo0,iblo1,iblo2,iblo3
     &           ,ibhi0,ibhi1,ibhi2,ibhi3
     &           ,nbcomp
     &           ,bstart
     &           ,ncomp
     &           ,igridboxlo0,igridboxlo1,igridboxlo2,igridboxlo3
     &           ,igridboxhi0,igridboxhi1,igridboxhi2,igridboxhi3
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iclo0,iclo1,iclo2,iclo3
      integer ichi0,ichi1,ichi2,ichi3
      REAL_T c(
     &           iclo0:ichi0,
     &           iclo1:ichi1,
     &           iclo2:ichi2,
     &           iclo3:ichi3)
      integer nacomp
      integer ialo0,ialo1,ialo2,ialo3
      integer iahi0,iahi1,iahi2,iahi3
      REAL_T a(
     &           ialo0:iahi0,
     &           ialo1:iahi1,
     &           ialo2:iahi2,
     &           ialo3:iahi3,
     &           0:nacomp-1)
      integer astart
      integer nbcomp
      integer iblo0,iblo1,iblo2,iblo3
      integer ibhi0,ibhi1,ibhi2,ibhi3
      REAL_T b(
     &           iblo0:ibhi0,
     &           iblo1:ibhi1,
     &           iblo2:ibhi2,
     &           iblo3:ibhi3,
     &           0:nbcomp-1)
      integer bstart
      integer ncomp
      integer igridboxlo0,igridboxlo1,igridboxlo2,igridboxlo3
      integer igridboxhi0,igridboxhi1,igridboxhi2,igridboxhi3
      integer i,j,k,l, kk
      integer comp, na, nb
      
      do l = igridboxlo3,igridboxhi3
      do k = igridboxlo2,igridboxhi2
      do j = igridboxlo1,igridboxhi1
      do i = igridboxlo0,igridboxhi0

         c(i,j,k,l) = zero
      
      enddo
      enddo
      enddo
      enddo
      kk = iblo2
      do comp=0, ncomp-1
         na = astart+comp
         nb = bstart+comp
         
      do l = igridboxlo3,igridboxhi3
      do k = igridboxlo2,igridboxhi2
      do j = igridboxlo1,igridboxhi1
      do i = igridboxlo0,igridboxhi0

#if CH_SPACEDIM==5
           kk = k
#endif
           c(i,j,k,l) = c(i,j,k,l)
     &     +a(i,j,k,l,na)*b(i,j,kk,iblo3,nb)
         
      enddo
      enddo
      enddo
      enddo
      enddo
      return
      end
      subroutine PSCSPOINTDOTPRODVEL_4D(
     &           c
     &           ,iclo0,iclo1,iclo2,iclo3
     &           ,ichi0,ichi1,ichi2,ichi3
     &           ,a
     &           ,ialo0,ialo1,ialo2,ialo3
     &           ,iahi0,iahi1,iahi2,iahi3
     &           ,nacomp
     &           ,astart
     &           ,b
     &           ,iblo0,iblo1,iblo2,iblo3
     &           ,ibhi0,ibhi1,ibhi2,ibhi3
     &           ,nbcomp
     &           ,bstart
     &           ,ncomp
     &           ,igridboxlo0,igridboxlo1,igridboxlo2,igridboxlo3
     &           ,igridboxhi0,igridboxhi1,igridboxhi2,igridboxhi3
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iclo0,iclo1,iclo2,iclo3
      integer ichi0,ichi1,ichi2,ichi3
      REAL_T c(
     &           iclo0:ichi0,
     &           iclo1:ichi1,
     &           iclo2:ichi2,
     &           iclo3:ichi3)
      integer nacomp
      integer ialo0,ialo1,ialo2,ialo3
      integer iahi0,iahi1,iahi2,iahi3
      REAL_T a(
     &           ialo0:iahi0,
     &           ialo1:iahi1,
     &           ialo2:iahi2,
     &           ialo3:iahi3,
     &           0:nacomp-1)
      integer astart
      integer nbcomp
      integer iblo0,iblo1,iblo2,iblo3
      integer ibhi0,ibhi1,ibhi2,ibhi3
      REAL_T b(
     &           iblo0:ibhi0,
     &           iblo1:ibhi1,
     &           iblo2:ibhi2,
     &           iblo3:ibhi3,
     &           0:nbcomp-1)
      integer bstart
      integer ncomp
      integer igridboxlo0,igridboxlo1,igridboxlo2,igridboxlo3
      integer igridboxhi0,igridboxhi1,igridboxhi2,igridboxhi3
      integer i,j,k,l, kk
      integer comp, na, nb
      
      do l = igridboxlo3,igridboxhi3
      do k = igridboxlo2,igridboxhi2
      do j = igridboxlo1,igridboxhi1
      do i = igridboxlo0,igridboxhi0

         c(i,j,k,l) = zero
      
      enddo
      enddo
      enddo
      enddo
      do comp=0, ncomp-1
         na = astart+comp
         nb = bstart+comp
         
      do l = igridboxlo3,igridboxhi3
      do k = igridboxlo2,igridboxhi2
      do j = igridboxlo1,igridboxhi1
      do i = igridboxlo0,igridboxhi0

          kk = k
#if CH_SPACEDIM==5
          kk = iblo2
#endif
          c(i,j,k,l) = c(i,j,k,l)
     &     +a(i,j,k,l,na)*b(iblo0,iblo1,kk,l,nb)
         
      enddo
      enddo
      enddo
      enddo
      enddo
      return
      end
      subroutine PSCSPOINTWISEREDUCEDDOTPRODCFG_4D(
     &           c
     &           ,iclo0,iclo1,iclo2,iclo3
     &           ,ichi0,ichi1,ichi2,ichi3
     &           ,a
     &           ,ialo0,ialo1,ialo2,ialo3
     &           ,iahi0,iahi1,iahi2,iahi3
     &           ,nacomp
     &           ,astart
     &           ,b
     &           ,iblo0,iblo1,iblo2,iblo3
     &           ,ibhi0,ibhi1,ibhi2,ibhi3
     &           ,nbcomp
     &           ,bstart
     &           ,numcomp
     &           ,igridboxlo0,igridboxlo1,igridboxlo2,igridboxlo3
     &           ,igridboxhi0,igridboxhi1,igridboxhi2,igridboxhi3
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iclo0,iclo1,iclo2,iclo3
      integer ichi0,ichi1,ichi2,ichi3
      REAL_T c(
     &           iclo0:ichi0,
     &           iclo1:ichi1,
     &           iclo2:ichi2,
     &           iclo3:ichi3)
      integer nacomp
      integer ialo0,ialo1,ialo2,ialo3
      integer iahi0,iahi1,iahi2,iahi3
      REAL_T a(
     &           ialo0:iahi0,
     &           ialo1:iahi1,
     &           ialo2:iahi2,
     &           ialo3:iahi3,
     &           0:nacomp-1)
      integer astart
      integer nbcomp
      integer iblo0,iblo1,iblo2,iblo3
      integer ibhi0,ibhi1,ibhi2,ibhi3
      REAL_T b(
     &           iblo0:ibhi0,
     &           iblo1:ibhi1,
     &           iblo2:ibhi2,
     &           iblo3:ibhi3,
     &           0:nbcomp-1)
      integer bstart
      integer numcomp
      integer igridboxlo0,igridboxlo1,igridboxlo2,igridboxlo3
      integer igridboxhi0,igridboxhi1,igridboxhi2,igridboxhi3
      integer i,j,k,l, kk
      integer comp, na, nb
      kk = ialo2
      
      do l = igridboxlo3,igridboxhi3
      do k = igridboxlo2,igridboxhi2
      do j = igridboxlo1,igridboxhi1
      do i = igridboxlo0,igridboxhi0

         c(i,j,k,l) = zero
      
      enddo
      enddo
      enddo
      enddo
      do comp=0, numcomp-1
         na = comp + astart;
         nb = comp + bstart;
         
      do l = igridboxlo3,igridboxhi3
      do k = igridboxlo2,igridboxhi2
      do j = igridboxlo1,igridboxhi1
      do i = igridboxlo0,igridboxhi0

#if CH_SPACEDIM==5
           kk = k
#endif
           c(i,j,k,l) = c(i,j,k,l)
     &  +a(i,j,kk,ialo3,na)*b(i,j,k,l,nb)
         
      enddo
      enddo
      enddo
      enddo
      enddo
      return
      end
      subroutine PSCSPOINTWISEREDUCEDDOTPRODVEL_4D(
     &           c
     &           ,iclo0,iclo1,iclo2,iclo3
     &           ,ichi0,ichi1,ichi2,ichi3
     &           ,a
     &           ,ialo0,ialo1,ialo2,ialo3
     &           ,iahi0,iahi1,iahi2,iahi3
     &           ,nacomp
     &           ,astart
     &           ,b
     &           ,iblo0,iblo1,iblo2,iblo3
     &           ,ibhi0,ibhi1,ibhi2,ibhi3
     &           ,nbcomp
     &           ,bstart
     &           ,numcomp
     &           ,igridboxlo0,igridboxlo1,igridboxlo2,igridboxlo3
     &           ,igridboxhi0,igridboxhi1,igridboxhi2,igridboxhi3
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iclo0,iclo1,iclo2,iclo3
      integer ichi0,ichi1,ichi2,ichi3
      REAL_T c(
     &           iclo0:ichi0,
     &           iclo1:ichi1,
     &           iclo2:ichi2,
     &           iclo3:ichi3)
      integer nacomp
      integer ialo0,ialo1,ialo2,ialo3
      integer iahi0,iahi1,iahi2,iahi3
      REAL_T a(
     &           ialo0:iahi0,
     &           ialo1:iahi1,
     &           ialo2:iahi2,
     &           ialo3:iahi3,
     &           0:nacomp-1)
      integer astart
      integer nbcomp
      integer iblo0,iblo1,iblo2,iblo3
      integer ibhi0,ibhi1,ibhi2,ibhi3
      REAL_T b(
     &           iblo0:ibhi0,
     &           iblo1:ibhi1,
     &           iblo2:ibhi2,
     &           iblo3:ibhi3,
     &           0:nbcomp-1)
      integer bstart
      integer numcomp
      integer igridboxlo0,igridboxlo1,igridboxlo2,igridboxlo3
      integer igridboxhi0,igridboxhi1,igridboxhi2,igridboxhi3
      integer i,j,k,l, kk
      integer comp, na, nb
      
      do l = igridboxlo3,igridboxhi3
      do k = igridboxlo2,igridboxhi2
      do j = igridboxlo1,igridboxhi1
      do i = igridboxlo0,igridboxhi0

         c(i,j,k,l) = zero
      
      enddo
      enddo
      enddo
      enddo
      do comp=0, numcomp-1
         na = comp + astart;
         nb = comp + bstart;
         
      do l = igridboxlo3,igridboxhi3
      do k = igridboxlo2,igridboxhi2
      do j = igridboxlo1,igridboxhi1
      do i = igridboxlo0,igridboxhi0

          kk = k
#if CH_SPACEDIM==5
          kk = ialo2
#endif
           c(i,j,k,l) = c(i,j,k,l)
     &  +a(ialo0,ialo1,kk,l,na)*b(i,j,k,l,nb)
         
      enddo
      enddo
      enddo
      enddo
      enddo
      return
      end
      subroutine PSCSTANGRADFACE_4D(
     &           gradPhi
     &           ,igradPhilo0,igradPhilo1,igradPhilo2,igradPhilo3
     &           ,igradPhihi0,igradPhihi1,igradPhihi2,igradPhihi3
     &           ,phi
     &           ,iphilo0,iphilo1,iphilo2,iphilo3
     &           ,iphihi0,iphihi1,iphihi2,iphihi3
     &           ,igradBoxlo0,igradBoxlo1,igradBoxlo2,igradBoxlo3
     &           ,igradBoxhi0,igradBoxhi1,igradBoxhi2,igradBoxhi3
     &           ,gradDir
     &           ,dx
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer igradPhilo0,igradPhilo1,igradPhilo2,igradPhilo3
      integer igradPhihi0,igradPhihi1,igradPhihi2,igradPhihi3
      REAL_T gradPhi(
     &           igradPhilo0:igradPhihi0,
     &           igradPhilo1:igradPhihi1,
     &           igradPhilo2:igradPhihi2,
     &           igradPhilo3:igradPhihi3)
      integer iphilo0,iphilo1,iphilo2,iphilo3
      integer iphihi0,iphihi1,iphihi2,iphihi3
      REAL_T phi(
     &           iphilo0:iphihi0,
     &           iphilo1:iphihi1,
     &           iphilo2:iphihi2,
     &           iphilo3:iphihi3)
      integer igradBoxlo0,igradBoxlo1,igradBoxlo2,igradBoxlo3
      integer igradBoxhi0,igradBoxhi1,igradBoxhi2,igradBoxhi3
      integer gradDir
      REAL_T dx
      integer i,j,k,l
      integer ii,jj,kk,ll
      REAL_T halfOnDx
      halfOnDx = half/dx
      ii = CHF_ID(0,gradDir)
                jj = CHF_ID(1,gradDir)
                kk = CHF_ID(2,gradDir)
                ll = CHF_ID(3,gradDir)
      
      do l = igradBoxlo3,igradBoxhi3
      do k = igradBoxlo2,igradBoxhi2
      do j = igradBoxlo1,igradBoxhi1
      do i = igradBoxlo0,igradBoxhi0

      gradPhi(i,j,k,l) = halfOnDx*(phi(i+ii,j+jj,k+kk,l+ll)
     &                                    -phi(i-ii,j-jj,k-kk,l-ll))
      
      enddo
      enddo
      enddo
      enddo
      return
      end
      subroutine FLUXNORMALDIVERGENCE_4D(
     &           Flux
     &           ,iFluxlo0,iFluxlo1,iFluxlo2,iFluxlo3
     &           ,iFluxhi0,iFluxhi1,iFluxhi2,iFluxhi3
     &           ,nFluxcomp
     &           ,div
     &           ,idivlo0,idivlo1,idivlo2,idivlo3
     &           ,idivhi0,idivhi1,idivhi2,idivhi3
     &           ,ndivcomp
     &           ,igridIntlo0,igridIntlo1,igridIntlo2,igridIntlo3
     &           ,igridInthi0,igridInthi1,igridInthi2,igridInthi3
     &           ,dx
     &           ,dir
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer nFluxcomp
      integer iFluxlo0,iFluxlo1,iFluxlo2,iFluxlo3
      integer iFluxhi0,iFluxhi1,iFluxhi2,iFluxhi3
      REAL_T Flux(
     &           iFluxlo0:iFluxhi0,
     &           iFluxlo1:iFluxhi1,
     &           iFluxlo2:iFluxhi2,
     &           iFluxlo3:iFluxhi3,
     &           0:nFluxcomp-1)
      integer ndivcomp
      integer idivlo0,idivlo1,idivlo2,idivlo3
      integer idivhi0,idivhi1,idivhi2,idivhi3
      REAL_T div(
     &           idivlo0:idivhi0,
     &           idivlo1:idivhi1,
     &           idivlo2:idivhi2,
     &           idivlo3:idivhi3,
     &           0:ndivcomp-1)
      integer igridIntlo0,igridIntlo1,igridIntlo2,igridIntlo3
      integer igridInthi0,igridInthi1,igridInthi2,igridInthi3
      REAL_T dx(0:3)
      integer dir
      integer i,j,k,l
      integer ii,jj,kk,ll
      integer comp
      do comp=0,ndivcomp-1
      
      do l = igridIntlo3,igridInthi3
      do k = igridIntlo2,igridInthi2
      do j = igridIntlo1,igridInthi1
      do i = igridIntlo0,igridInthi0

      
      ii = i+CHF_ID(0,dir)
      jj = j+CHF_ID(1,dir)
      kk = k+CHF_ID(2,dir)
      ll = l+CHF_ID(3,dir)
      div(i,j,k,l,comp) = div(i,j,k,l,comp)
     &                              + ( Flux(ii,jj,kk,ll,comp)
     &                              -   Flux(i,j,k,l,comp) )/dx(dir)
      
      enddo
      enddo
      enddo
      enddo
      enddo
      return
      end
      subroutine PSCSDIVERGENCE_4D(
     &           uEdge
     &           ,iuEdgelo0,iuEdgelo1,iuEdgelo2,iuEdgelo3
     &           ,iuEdgehi0,iuEdgehi1,iuEdgehi2,iuEdgehi3
     &           ,nuEdgecomp
     &           ,div
     &           ,idivlo0,idivlo1,idivlo2,idivlo3
     &           ,idivhi0,idivhi1,idivhi2,idivhi3
     &           ,ndivcomp
     &           ,igridIntlo0,igridIntlo1,igridIntlo2,igridIntlo3
     &           ,igridInthi0,igridInthi1,igridInthi2,igridInthi3
     &           ,dx
     &           ,dir
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer nuEdgecomp
      integer iuEdgelo0,iuEdgelo1,iuEdgelo2,iuEdgelo3
      integer iuEdgehi0,iuEdgehi1,iuEdgehi2,iuEdgehi3
      REAL_T uEdge(
     &           iuEdgelo0:iuEdgehi0,
     &           iuEdgelo1:iuEdgehi1,
     &           iuEdgelo2:iuEdgehi2,
     &           iuEdgelo3:iuEdgehi3,
     &           0:nuEdgecomp-1)
      integer ndivcomp
      integer idivlo0,idivlo1,idivlo2,idivlo3
      integer idivhi0,idivhi1,idivhi2,idivhi3
      REAL_T div(
     &           idivlo0:idivhi0,
     &           idivlo1:idivhi1,
     &           idivlo2:idivhi2,
     &           idivlo3:idivhi3,
     &           0:ndivcomp-1)
      integer igridIntlo0,igridIntlo1,igridIntlo2,igridIntlo3
      integer igridInthi0,igridInthi1,igridInthi2,igridInthi3
      REAL_T dx
      integer dir
      integer i,j,k,l
      integer ii,jj,kk,ll
      integer comp
      REAL_T one_on_dx
      one_on_dx = one/dx
      do comp=0,ndivcomp-1
      
      do l = igridIntlo3,igridInthi3
      do k = igridIntlo2,igridInthi2
      do j = igridIntlo1,igridInthi1
      do i = igridIntlo0,igridInthi0

      
      ii = i+CHF_ID(0,dir)
      jj = j+CHF_ID(1,dir)
      kk = k+CHF_ID(2,dir)
      ll = l+CHF_ID(3,dir)
      div(i,j,k,l,comp) = div(i,j,k,l,comp)
     &                      +one_on_dx*(uEdge(ii,jj,kk,ll,comp)
     &                                 -uEdge(i,j,k,l,comp))
      
      enddo
      enddo
      enddo
      enddo
      enddo
      return
      end
      subroutine MULT_CFG_NT_4D(
     &           igridboxlo0,igridboxlo1,igridboxlo2,igridboxlo3
     &           ,igridboxhi0,igridboxhi1,igridboxhi2,igridboxhi3
     &           ,N
     &           ,iNlo0,iNlo1,iNlo2,iNlo3
     &           ,iNhi0,iNhi1,iNhi2,iNhi3
     &           ,nNcomp
     &           ,velocity
     &           ,ivelocitylo0,ivelocitylo1,ivelocitylo2,ivelocitylo3
     &           ,ivelocityhi0,ivelocityhi1,ivelocityhi2,ivelocityhi3
     &           ,nvelocitycomp
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer igridboxlo0,igridboxlo1,igridboxlo2,igridboxlo3
      integer igridboxhi0,igridboxhi1,igridboxhi2,igridboxhi3
      integer nNcomp
      integer iNlo0,iNlo1,iNlo2,iNlo3
      integer iNhi0,iNhi1,iNhi2,iNhi3
      REAL_T N(
     &           iNlo0:iNhi0,
     &           iNlo1:iNhi1,
     &           iNlo2:iNhi2,
     &           iNlo3:iNhi3,
     &           0:nNcomp-1)
      integer nvelocitycomp
      integer ivelocitylo0,ivelocitylo1,ivelocitylo2,ivelocitylo3
      integer ivelocityhi0,ivelocityhi1,ivelocityhi2,ivelocityhi3
      REAL_T velocity(
     &           ivelocitylo0:ivelocityhi0,
     &           ivelocitylo1:ivelocityhi1,
     &           ivelocitylo2:ivelocityhi2,
     &           ivelocitylo3:ivelocityhi3,
     &           0:nvelocitycomp-1)
      double precision NT11, NT12, NT13, NT21, NT22, NT23, NT31, NT32, NT33, v0, v1, v2
      integer i0,i1,i2,i3, kk, k12, k21, k22
      kk = iNlo2
      k12 = 2
      k22 = 3
#if CH_SPACEDIM==5
       k12 = 3
       k22 = 4
#endif
      
      do i3 = igridboxlo3,igridboxhi3
      do i2 = igridboxlo2,igridboxhi2
      do i1 = igridboxlo1,igridboxhi1
      do i0 = igridboxlo0,igridboxhi0

#if CH_SPACEDIM==5
         kk = i2
#endif
         NT11 = N(i0,i1,kk,iNlo3,0)
         NT12 = N(i0,i1,kk,iNlo3,k12)
         NT21 = N(i0,i1,kk,iNlo3,1)
         NT22 = N(i0,i1,kk,iNlo3,k22)
#if CH_SPACEDIM==5
         NT13 = N(i0,i1,kk,iNlo3,6)
         NT23 = N(i0,i1,kk,iNlo3,7)
         NT31 = N(i0,i1,kk,iNlo3,2)
         NT32 = N(i0,i1,kk,iNlo3,5)
         NT33 = N(i0,i1,kk,iNlo3,8)
#endif
         v0 = NT11 * velocity(i0,i1,i2,i3,0) + NT12 * velocity(i0,i1,i2,i3,1)
         v1 = NT21 * velocity(i0,i1,i2,i3,0) + NT22 * velocity(i0,i1,i2,i3,1)
#if CH_SPACEDIM==5
         v0 = v0 + NT13*velocity(i0,i1,i2,i3,2)
         v1 = v1 + NT23*velocity(i0,i1,i2,i3,2)
         v2 = NT31 * velocity(i0,i1,i2,i3,0) + NT32 * velocity(i0,i1,i2,i3,1)
     &            + NT33 * velocity(i0,i1,i2,i3,2)
#endif
         velocity(i0,i1,i2,i3,0) = v0
         velocity(i0,i1,i2,i3,1) = v1
#if CH_SPACEDIM==5
         velocity(i0,i1,i2,i3,2) = v2
#endif
      
      enddo
      enddo
      enddo
      enddo
      return
      end
      subroutine MULT_VEL_NT_4D(
     &           igridboxlo0,igridboxlo1,igridboxlo2,igridboxlo3
     &           ,igridboxhi0,igridboxhi1,igridboxhi2,igridboxhi3
     &           ,TwoPirRmaj
     &           ,iTwoPirRmajlo0,iTwoPirRmajlo1,iTwoPirRmajlo2,iTwoPirRmajlo3
     &           ,iTwoPirRmajhi0,iTwoPirRmajhi1,iTwoPirRmajhi2,iTwoPirRmajhi3
     &           ,velocity
     &           ,ivelocitylo0,ivelocitylo1,ivelocitylo2,ivelocitylo3
     &           ,ivelocityhi0,ivelocityhi1,ivelocityhi2,ivelocityhi3
     &           ,nvelocitycomp
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer igridboxlo0,igridboxlo1,igridboxlo2,igridboxlo3
      integer igridboxhi0,igridboxhi1,igridboxhi2,igridboxhi3
      integer iTwoPirRmajlo0,iTwoPirRmajlo1,iTwoPirRmajlo2,iTwoPirRmajlo3
      integer iTwoPirRmajhi0,iTwoPirRmajhi1,iTwoPirRmajhi2,iTwoPirRmajhi3
      REAL_T TwoPirRmaj(
     &           iTwoPirRmajlo0:iTwoPirRmajhi0,
     &           iTwoPirRmajlo1:iTwoPirRmajhi1,
     &           iTwoPirRmajlo2:iTwoPirRmajhi2,
     &           iTwoPirRmajlo3:iTwoPirRmajhi3)
      integer nvelocitycomp
      integer ivelocitylo0,ivelocitylo1,ivelocitylo2,ivelocitylo3
      integer ivelocityhi0,ivelocityhi1,ivelocityhi2,ivelocityhi3
      REAL_T velocity(
     &           ivelocitylo0:ivelocityhi0,
     &           ivelocitylo1:ivelocityhi1,
     &           ivelocitylo2:ivelocityhi2,
     &           ivelocitylo3:ivelocityhi3,
     &           0:nvelocitycomp-1)
      integer i0,i1,i2,i3, dir, kk
      kk = iTwoPirRmajlo2
      
      do i3 = igridboxlo3,igridboxhi3
      do i2 = igridboxlo2,igridboxhi2
      do i1 = igridboxlo1,igridboxhi1
      do i0 = igridboxlo0,igridboxhi0

#if CH_SPACEDIM==5
         kk = i2
#endif
         do dir = CH_SPACEDIM-2, CH_SPACEDIM-1
            velocity(i0,i1,i2,i3,dir) = velocity(i0,i1,i2,i3,dir)
     &           * TwoPirRmaj(i0,i1,kk,iTwoPirRmajlo3)
         enddo
      
      enddo
      enddo
      enddo
      enddo
      return
      end
      subroutine MULT_CFG_4D(
     &           igridboxlo0,igridboxlo1,igridboxlo2,igridboxlo3
     &           ,igridboxhi0,igridboxhi1,igridboxhi2,igridboxhi3
     &           ,factor
     &           ,ifactorlo0,ifactorlo1,ifactorlo2,ifactorlo3
     &           ,ifactorhi0,ifactorhi1,ifactorhi2,ifactorhi3
     &           ,data
     &           ,idatalo0,idatalo1,idatalo2,idatalo3
     &           ,idatahi0,idatahi1,idatahi2,idatahi3
     &           ,ndatacomp
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer igridboxlo0,igridboxlo1,igridboxlo2,igridboxlo3
      integer igridboxhi0,igridboxhi1,igridboxhi2,igridboxhi3
      integer ifactorlo0,ifactorlo1,ifactorlo2,ifactorlo3
      integer ifactorhi0,ifactorhi1,ifactorhi2,ifactorhi3
      REAL_T factor(
     &           ifactorlo0:ifactorhi0,
     &           ifactorlo1:ifactorhi1,
     &           ifactorlo2:ifactorhi2,
     &           ifactorlo3:ifactorhi3)
      integer ndatacomp
      integer idatalo0,idatalo1,idatalo2,idatalo3
      integer idatahi0,idatahi1,idatahi2,idatahi3
      REAL_T data(
     &           idatalo0:idatahi0,
     &           idatalo1:idatahi1,
     &           idatalo2:idatahi2,
     &           idatalo3:idatahi3,
     &           0:ndatacomp-1)
      integer i0,i1,i2,i3, comp, kk
      kk = ifactorlo2
      do comp = 0, ndatacomp-1
         
      do i3 = igridBoxlo3,igridBoxhi3
      do i2 = igridBoxlo2,igridBoxhi2
      do i1 = igridBoxlo1,igridBoxhi1
      do i0 = igridBoxlo0,igridBoxhi0

#if CH_SPACEDIM==5
            kk = i2
#endif
            data(i0,i1,i2,i3,comp) = data(i0,i1,i2,i3,comp)
     &        * factor(i0,i1,kk,ifactorlo3)
         
      enddo
      enddo
      enddo
      enddo
      enddo
      return
      end
      subroutine MULT_VEL_4D(
     &           igridboxlo0,igridboxlo1,igridboxlo2,igridboxlo3
     &           ,igridboxhi0,igridboxhi1,igridboxhi2,igridboxhi3
     &           ,factor
     &           ,ifactorlo0,ifactorlo1,ifactorlo2,ifactorlo3
     &           ,ifactorhi0,ifactorhi1,ifactorhi2,ifactorhi3
     &           ,data
     &           ,idatalo0,idatalo1,idatalo2,idatalo3
     &           ,idatahi0,idatahi1,idatahi2,idatahi3
     &           ,ndatacomp
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer igridboxlo0,igridboxlo1,igridboxlo2,igridboxlo3
      integer igridboxhi0,igridboxhi1,igridboxhi2,igridboxhi3
      integer ifactorlo0,ifactorlo1,ifactorlo2,ifactorlo3
      integer ifactorhi0,ifactorhi1,ifactorhi2,ifactorhi3
      REAL_T factor(
     &           ifactorlo0:ifactorhi0,
     &           ifactorlo1:ifactorhi1,
     &           ifactorlo2:ifactorhi2,
     &           ifactorlo3:ifactorhi3)
      integer ndatacomp
      integer idatalo0,idatalo1,idatalo2,idatalo3
      integer idatahi0,idatahi1,idatahi2,idatahi3
      REAL_T data(
     &           idatalo0:idatahi0,
     &           idatalo1:idatahi1,
     &           idatalo2:idatahi2,
     &           idatalo3:idatahi3,
     &           0:ndatacomp-1)
      integer i0,i1,i2,i3, comp, kk
      do comp = 0, ndatacomp-1
         
      do i3 = igridBoxlo3,igridBoxhi3
      do i2 = igridBoxlo2,igridBoxhi2
      do i1 = igridBoxlo1,igridBoxhi1
      do i0 = igridBoxlo0,igridBoxhi0

            kk = i2
#if CH_SPACEDIM==5
            kk = ifactorlo2
#endif
            data(i0,i1,i2,i3,comp) = data(i0,i1,i2,i3,comp)
     &        * factor(ifactorlo0,ifactorlo1,kk,i3)
         
      enddo
      enddo
      enddo
      enddo
      enddo
      return
      end
      subroutine MULT_RADIAL_COMPONENT_4D(
     &           igridboxlo0,igridboxlo1,igridboxlo2,igridboxlo3
     &           ,igridboxhi0,igridboxhi1,igridboxhi2,igridboxhi3
     &           ,factor
     &           ,ifactorlo0,ifactorlo1,ifactorlo2,ifactorlo3
     &           ,ifactorhi0,ifactorhi1,ifactorhi2,ifactorhi3
     &           ,data
     &           ,idatalo0,idatalo1,idatalo2,idatalo3
     &           ,idatahi0,idatahi1,idatahi2,idatahi3
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer igridboxlo0,igridboxlo1,igridboxlo2,igridboxlo3
      integer igridboxhi0,igridboxhi1,igridboxhi2,igridboxhi3
      integer ifactorlo0,ifactorlo1,ifactorlo2,ifactorlo3
      integer ifactorhi0,ifactorhi1,ifactorhi2,ifactorhi3
      REAL_T factor(
     &           ifactorlo0:ifactorhi0,
     &           ifactorlo1:ifactorhi1,
     &           ifactorlo2:ifactorhi2,
     &           ifactorlo3:ifactorhi3)
      integer idatalo0,idatalo1,idatalo2,idatalo3
      integer idatahi0,idatahi1,idatahi2,idatahi3
      REAL_T data(
     &           idatalo0:idatahi0,
     &           idatalo1:idatahi1,
     &           idatalo2:idatahi2,
     &           idatalo3:idatahi3)
      integer i0,i1,i2,i3,kk
      kk = ifactorlo2
      
      do i3 = igridBoxlo3,igridBoxhi3
      do i2 = igridBoxlo2,igridBoxhi2
      do i1 = igridBoxlo1,igridBoxhi1
      do i0 = igridBoxlo0,igridBoxhi0

#if CH_SPACEDIM==5
         kk = i2
#endif
         data(i0,i1,i2,i3) = data(i0,i1,i2,i3)
     &        * factor(i0,i1,kk,ifactorlo3)
      
      enddo
      enddo
      enddo
      enddo
      return
      end
      subroutine COMPUTE_Q0X_4D(
     &           igridboxlo0,igridboxlo1,igridboxlo2,igridboxlo3
     &           ,igridboxhi0,igridboxhi1,igridboxhi2,igridboxhi3
     &           ,integrals
     &           ,iintegralslo0,iintegralslo1,iintegralslo2,iintegralslo3
     &           ,iintegralshi0,iintegralshi1,iintegralshi2,iintegralshi3
     &           ,nintegralscomp
     &           ,dx
     &           ,charge
     &           ,larmor
     &           ,include_ExB_drift
     &           ,include_magnetic_drift
     &           ,Q
     &           ,iQlo0,iQlo1,iQlo2,iQlo3
     &           ,iQhi0,iQhi1,iQhi2,iQhi3
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer igridboxlo0,igridboxlo1,igridboxlo2,igridboxlo3
      integer igridboxhi0,igridboxhi1,igridboxhi2,igridboxhi3
      integer nintegralscomp
      integer iintegralslo0,iintegralslo1,iintegralslo2,iintegralslo3
      integer iintegralshi0,iintegralshi1,iintegralshi2,iintegralshi3
      REAL_T integrals(
     &           iintegralslo0:iintegralshi0,
     &           iintegralslo1:iintegralshi1,
     &           iintegralslo2:iintegralshi2,
     &           iintegralslo3:iintegralshi3,
     &           0:nintegralscomp-1)
      REAL_T dx(0:3)
      REAL_T charge
      REAL_T larmor
      integer include_ExB_drift
      integer include_magnetic_drift
      integer iQlo0,iQlo1,iQlo2,iQlo3
      integer iQhi0,iQhi1,iQhi2,iQhi3
      REAL_T Q(
     &           iQlo0:iQhi0,
     &           iQlo1:iQhi1,
     &           iQlo2:iQhi2,
     &           iQlo3:iQhi3)
      integer i,j,k,l, comp, kk
      double precision eta(0:1), vpar, mu_lo, mu_hi, mu_integral0, mu_integral1, sum
      kk = iintegralslo2
      
      do l = igridBoxlo3,igridBoxhi3
      do k = igridBoxlo2,igridBoxhi2
      do j = igridBoxlo1,igridBoxhi1
      do i = igridBoxlo0,igridBoxhi0

#if CH_SPACEDIM==5
         kk = k
         vpar = l * dx(3)
         mu_lo = m * dx(4)
         mu_hi = mu_lo + dx(4)
         mu_integral0 = dx(4)
         mu_integral1 = half * (mu_hi**2 - mu_lo**2)
#else
         vpar = k * dx(2)
         mu_lo = l * dx(3)
         mu_hi = mu_lo + dx(3)
         mu_integral0 = dx(3)
         mu_integral1 = half * (mu_hi**2 - mu_lo**2)
#endif
         if (include_ExB_drift .ne. 0) then
            eta(0) = mu_integral0
         else
            eta(0) = zero
         endif            
         if (include_magnetic_drift .ne. 0) then
            eta(1) = half * mu_integral1 / charge
         else
            eta(1) = zero
         endif
         sum = zero
         do comp = 0, 1
            sum = sum + eta(comp) * integrals(i,j,kk,iintegralslo3, comp)
         enddo
         Q(i,j,k,l) = vpar * larmor * sum
      
      enddo
      enddo
      enddo
      enddo
      return
      end
      subroutine COMPUTE_Q12_4D(
     &           igridboxlo0,igridboxlo1,igridboxlo2,igridboxlo3
     &           ,igridboxhi0,igridboxhi1,igridboxhi2,igridboxhi3
     &           ,dir
     &           ,side
     &           ,nodal_integrals
     &           ,inodal_integralslo0,inodal_integralslo1,inodal_integralslo2,inodal_integralslo3
     &           ,inodal_integralshi0,inodal_integralshi1,inodal_integralshi2,inodal_integralshi3
     &           ,nnodal_integralscomp
     &           ,dx
     &           ,mass
     &           ,charge
     &           ,larmor
     &           ,no_par_stream
     &           ,include_magnetic_drift
     &           ,Q
     &           ,iQlo0,iQlo1,iQlo2,iQlo3
     &           ,iQhi0,iQhi1,iQhi2,iQhi3
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer igridboxlo0,igridboxlo1,igridboxlo2,igridboxlo3
      integer igridboxhi0,igridboxhi1,igridboxhi2,igridboxhi3
      integer dir
      integer side
      integer nnodal_integralscomp
      integer inodal_integralslo0,inodal_integralslo1,inodal_integralslo2,inodal_integralslo3
      integer inodal_integralshi0,inodal_integralshi1,inodal_integralshi2,inodal_integralshi3
      REAL_T nodal_integrals(
     &           inodal_integralslo0:inodal_integralshi0,
     &           inodal_integralslo1:inodal_integralshi1,
     &           inodal_integralslo2:inodal_integralshi2,
     &           inodal_integralslo3:inodal_integralshi3,
     &           0:nnodal_integralscomp-1)
      REAL_T dx(0:3)
      REAL_T mass
      REAL_T charge
      REAL_T larmor
      integer no_par_stream
      integer include_magnetic_drift
      integer iQlo0,iQlo1,iQlo2,iQlo3
      integer iQhi0,iQhi1,iQhi2,iQhi3
      REAL_T Q(
     &           iQlo0:iQhi0,
     &           iQlo1:iQhi1,
     &           iQlo2:iQhi2,
     &           iQlo3:iQhi3)
      integer i,j,k,l, ii,jj,kk,ll, tdir, comp, k_fix, k_kk_fix
      double precision eta(0:1), vpar_lo, vpar_hi, sum
      tdir = 1 - dir
      ii = CHF_ID(0,tdir)
                jj = CHF_ID(1,tdir)
                kk = CHF_ID(2,tdir)
                ll = CHF_ID(3,tdir)
      k_fix = inodal_integralslo2
      k_kk_fix = inodal_integralslo2
      
      do l = igridBoxlo3,igridBoxhi3
      do k = igridBoxlo2,igridBoxhi2
      do j = igridBoxlo1,igridBoxhi1
      do i = igridBoxlo0,igridBoxhi0

#if CH_SPACEDIM==5
         k_fix = k
         k_kk_fix = k + kk
         vpar_lo = l * dx(3)
         vpar_hi = vpar_lo + dx(3)
#else
         vpar_lo = k * dx(2)
         vpar_hi = vpar_lo + dx(2)
#endif
         if (no_par_stream .ne. 0) then
            eta(0) = zero
         else
            eta(0) = half * (vpar_hi**2 - vpar_lo**2)
         endif
         if (include_magnetic_drift .ne. 0) then
            eta(1) = third * (vpar_hi**3 - vpar_lo**3) * mass * larmor / charge
         else
            eta(1) = zero
         endif
         if (side .eq. 0) then
            sum = zero
            do comp = 0, 1
               sum = sum + eta(comp) * nodal_integrals(i,j,k_fix,inodal_integralslo3, comp)
            enddo
         else
            sum = zero
            do comp = 0, 1
               sum = sum + eta(comp) * nodal_integrals(i+ii,j+jj,k_kk_fix,inodal_integralslo3, comp)
            enddo
         endif
#if CH_SPACEDIM==5
         Q(i,j,k,l) = dx(4) * sum
#else
         Q(i,j,k,l) = dx(3) * sum
#endif
      
      enddo
      enddo
      enddo
      enddo
      return
      end
      subroutine COMPUTE_UHAT_4D(
     &           igridboxlo0,igridboxlo1,igridboxlo2,igridboxlo3
     &           ,igridboxhi0,igridboxhi1,igridboxhi2,igridboxhi3
     &           ,integrals
     &           ,iintegralslo0,iintegralslo1,iintegralslo2,iintegralslo3
     &           ,iintegralshi0,iintegralshi1,iintegralshi2,iintegralshi3
     &           ,nintegralscomp
     &           ,dx
     &           ,mass
     &           ,charge
     &           ,no_zero_order_par
     &           ,uhat
     &           ,iuhatlo0,iuhatlo1,iuhatlo2,iuhatlo3
     &           ,iuhathi0,iuhathi1,iuhathi2,iuhathi3
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer igridboxlo0,igridboxlo1,igridboxlo2,igridboxlo3
      integer igridboxhi0,igridboxhi1,igridboxhi2,igridboxhi3
      integer nintegralscomp
      integer iintegralslo0,iintegralslo1,iintegralslo2,iintegralslo3
      integer iintegralshi0,iintegralshi1,iintegralshi2,iintegralshi3
      REAL_T integrals(
     &           iintegralslo0:iintegralshi0,
     &           iintegralslo1:iintegralshi1,
     &           iintegralslo2:iintegralshi2,
     &           iintegralslo3:iintegralshi3,
     &           0:nintegralscomp-1)
      REAL_T dx(0:3)
      REAL_T mass
      REAL_T charge
      integer no_zero_order_par
      integer iuhatlo0,iuhatlo1,iuhatlo2,iuhatlo3
      integer iuhathi0,iuhathi1,iuhathi2,iuhathi3
      REAL_T uhat(
     &           iuhatlo0:iuhathi0,
     &           iuhatlo1:iuhathi1,
     &           iuhatlo2:iuhathi2,
     &           iuhatlo3:iuhathi3)
      integer i,j,k,l, comp, kk
      double precision mu_lo, mu_hi, mu_integral0, mu_integral1, cfg_area
      kk = iintegralslo2
      
      do l = igridBoxlo3,igridBoxhi3
      do k = igridBoxlo2,igridBoxhi2
      do j = igridBoxlo1,igridBoxhi1
      do i = igridBoxlo0,igridBoxhi0

#if CH_SPACEDIM==5
         kk = k
         mu_lo = m * dx(4)
         mu_hi = mu_lo + dx(4)
         mu_integral0 = dx(4)
#else
         mu_lo = l * dx(3)
         mu_hi = mu_lo + dx(3)
         mu_integral0 = dx(3)
#endif
         if (no_zero_order_par .ne. 0) then
           mu_integral1 = 0.0
         else 
           mu_integral1 = half * (mu_hi**2 - mu_lo**2)
         endif
         uhat(i,j,k,l) 
     &        = (mu_integral0 * charge * integrals(i,j,kk,iintegralslo3,0) 
     &        +  half * mu_integral1 * integrals(i,j,kk,iintegralslo3,1) ) / mass
      
      enddo
      enddo
      enddo
      enddo
      return
      end
      subroutine COMPUTE_RADIAL_PROJECTION_4D(
     &           iboxlo0,iboxlo1,iboxlo2,iboxlo3
     &           ,iboxhi0,iboxhi1,iboxhi2,iboxhi3
     &           ,bunit
     &           ,ibunitlo0,ibunitlo1,ibunitlo2,ibunitlo3
     &           ,ibunithi0,ibunithi1,ibunithi2,ibunithi3
     &           ,nbunitcomp
     &           ,vector
     &           ,ivectorlo0,ivectorlo1,ivectorlo2,ivectorlo3
     &           ,ivectorhi0,ivectorhi1,ivectorhi2,ivectorhi3
     &           ,nvectorcomp
     &           ,vector_r
     &           ,ivector_rlo0,ivector_rlo1,ivector_rlo2,ivector_rlo3
     &           ,ivector_rhi0,ivector_rhi1,ivector_rhi2,ivector_rhi3
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iboxlo0,iboxlo1,iboxlo2,iboxlo3
      integer iboxhi0,iboxhi1,iboxhi2,iboxhi3
      integer nbunitcomp
      integer ibunitlo0,ibunitlo1,ibunitlo2,ibunitlo3
      integer ibunithi0,ibunithi1,ibunithi2,ibunithi3
      REAL_T bunit(
     &           ibunitlo0:ibunithi0,
     &           ibunitlo1:ibunithi1,
     &           ibunitlo2:ibunithi2,
     &           ibunitlo3:ibunithi3,
     &           0:nbunitcomp-1)
      integer nvectorcomp
      integer ivectorlo0,ivectorlo1,ivectorlo2,ivectorlo3
      integer ivectorhi0,ivectorhi1,ivectorhi2,ivectorhi3
      REAL_T vector(
     &           ivectorlo0:ivectorhi0,
     &           ivectorlo1:ivectorhi1,
     &           ivectorlo2:ivectorhi2,
     &           ivectorlo3:ivectorhi3,
     &           0:nvectorcomp-1)
      integer ivector_rlo0,ivector_rlo1,ivector_rlo2,ivector_rlo3
      integer ivector_rhi0,ivector_rhi1,ivector_rhi2,ivector_rhi3
      REAL_T vector_r(
     &           ivector_rlo0:ivector_rhi0,
     &           ivector_rlo1:ivector_rhi1,
     &           ivector_rlo2:ivector_rhi2,
     &           ivector_rlo3:ivector_rhi3)
      integer i,j,k,l 
      double precision fac
      
      do l = iboxlo3,iboxhi3
      do k = iboxlo2,iboxhi2
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0

#if CH_SPACEDIM==4
         fac = sqrt(bunit(i,j,ibunitlo2,ibunitlo3,0) ** 2
     &             +bunit(i,j,ibunitlo2,ibunitlo3,2) ** 2)                
         vector_r(i,j,k,l) = vector(i,j,k,l,0) * bunit(i,j,ibunitlo2,ibunitlo3,2)/fac
     &                               - vector(i,j,k,l,1) * bunit(i,j,ibunitlo2,ibunitlo3,0)/fac 
#else
         fac = sqrt(bunit(i,j,k,ibunitlo3,ibunitlo4,0) ** 2
     &             +bunit(i,j,k,ibunitlo3,ibunitlo4,2) ** 2)                
         vector_r(i,j,k,l) = vector(i,j,k,l,0) * bunit(i,j,k,ibunitlo3,ibunitlo4,2)/fac
     &                               - vector(i,j,k,l,2) * bunit(i,j,k,ibunitlo3,ibunitlo4,0)/fac 
#endif
      
      enddo
      enddo
      enddo
      enddo
      return
      end
      subroutine COMPUTE_POLOIDAL_PROJECTION_4D(
     &           iboxlo0,iboxlo1,iboxlo2,iboxlo3
     &           ,iboxhi0,iboxhi1,iboxhi2,iboxhi3
     &           ,bunit
     &           ,ibunitlo0,ibunitlo1,ibunitlo2,ibunitlo3
     &           ,ibunithi0,ibunithi1,ibunithi2,ibunithi3
     &           ,nbunitcomp
     &           ,vector
     &           ,ivectorlo0,ivectorlo1,ivectorlo2,ivectorlo3
     &           ,ivectorhi0,ivectorhi1,ivectorhi2,ivectorhi3
     &           ,nvectorcomp
     &           ,vector_pol
     &           ,ivector_pollo0,ivector_pollo1,ivector_pollo2,ivector_pollo3
     &           ,ivector_polhi0,ivector_polhi1,ivector_polhi2,ivector_polhi3
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iboxlo0,iboxlo1,iboxlo2,iboxlo3
      integer iboxhi0,iboxhi1,iboxhi2,iboxhi3
      integer nbunitcomp
      integer ibunitlo0,ibunitlo1,ibunitlo2,ibunitlo3
      integer ibunithi0,ibunithi1,ibunithi2,ibunithi3
      REAL_T bunit(
     &           ibunitlo0:ibunithi0,
     &           ibunitlo1:ibunithi1,
     &           ibunitlo2:ibunithi2,
     &           ibunitlo3:ibunithi3,
     &           0:nbunitcomp-1)
      integer nvectorcomp
      integer ivectorlo0,ivectorlo1,ivectorlo2,ivectorlo3
      integer ivectorhi0,ivectorhi1,ivectorhi2,ivectorhi3
      REAL_T vector(
     &           ivectorlo0:ivectorhi0,
     &           ivectorlo1:ivectorhi1,
     &           ivectorlo2:ivectorhi2,
     &           ivectorlo3:ivectorhi3,
     &           0:nvectorcomp-1)
      integer ivector_pollo0,ivector_pollo1,ivector_pollo2,ivector_pollo3
      integer ivector_polhi0,ivector_polhi1,ivector_polhi2,ivector_polhi3
      REAL_T vector_pol(
     &           ivector_pollo0:ivector_polhi0,
     &           ivector_pollo1:ivector_polhi1,
     &           ivector_pollo2:ivector_polhi2,
     &           ivector_pollo3:ivector_polhi3)
      integer i,j,k,l 
      double precision fac
      
      do l = iboxlo3,iboxhi3
      do k = iboxlo2,iboxhi2
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0

#if CH_SPACEDIM==4
         fac = sqrt(bunit(i,j,ibunitlo2,ibunitlo3,0) ** 2
     &             +bunit(i,j,ibunitlo2,ibunitlo3,2) ** 2)                
         vector_pol(i,j,k,l) = vector(i,j,k,l,0) * bunit(i,j,ibunitlo2,ibunitlo3,0)/fac
     &                                 + vector(i,j,k,l,1) * bunit(i,j,ibunitlo2,ibunitlo3,2)/fac 
#else
         fac = sqrt(bunit(i,j,k,ibunitlo3,ibunitlo4,0) ** 2
     &             +bunit(i,j,k,ibunitlo3,ibunitlo4,2) ** 2)                
         vector_pol(i,j,k,l) = vector(i,j,k,l,0) * bunit(i,j,k,ibunitlo3,ibunitlo4,0)/fac
     &                                 + vector(i,j,k,l,2) * bunit(i,j,k,ibunitlo3,ibunitlo4,2)/fac 
#endif
      
      enddo
      enddo
      enddo
      enddo
      return
      end
