#include "REAL.H"
#include "SPACE.H"
#include "CONSTANTS.H"

#include "CONSTANTS.H"
      subroutine GET_SLAB_FIELD_DATA_2D(
     &           igridboxlo0,igridboxlo1
     &           ,igridboxhi0,igridboxhi1
     &           ,xix
     &           ,ixixlo0,ixixlo1
     &           ,ixixhi0,ixixhi1
     &           ,ByInner
     &           ,ByOuter
     &           ,BzInner
     &           ,BzOuter
     &           ,xmax
     &           ,ximax
     &           ,b_pt
     &           ,ib_ptlo0,ib_ptlo1
     &           ,ib_pthi0,ib_pthi1
     &           ,nb_ptcomp
     &           ,Bmag_pt
     &           ,iBmag_ptlo0,iBmag_ptlo1
     &           ,iBmag_pthi0,iBmag_pthi1
     &           ,bunit_pt
     &           ,ibunit_ptlo0,ibunit_ptlo1
     &           ,ibunit_pthi0,ibunit_pthi1
     &           ,nbunit_ptcomp
     &           ,gradb_pt
     &           ,igradb_ptlo0,igradb_ptlo1
     &           ,igradb_pthi0,igradb_pthi1
     &           ,ngradb_ptcomp
     &           ,curlb_pt
     &           ,icurlb_ptlo0,icurlb_ptlo1
     &           ,icurlb_pthi0,icurlb_pthi1
     &           ,ncurlb_ptcomp
     &           ,bdotcurlb_pt
     &           ,ibdotcurlb_ptlo0,ibdotcurlb_ptlo1
     &           ,ibdotcurlb_pthi0,ibdotcurlb_pthi1
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer igridboxlo0,igridboxlo1
      integer igridboxhi0,igridboxhi1
      integer ixixlo0,ixixlo1
      integer ixixhi0,ixixhi1
      REAL_T xix(
     &           ixixlo0:ixixhi0,
     &           ixixlo1:ixixhi1)
      REAL_T ByInner
      REAL_T ByOuter
      REAL_T BzInner
      REAL_T BzOuter
      REAL_T xmax
      REAL_T ximax
      integer nb_ptcomp
      integer ib_ptlo0,ib_ptlo1
      integer ib_pthi0,ib_pthi1
      REAL_T b_pt(
     &           ib_ptlo0:ib_pthi0,
     &           ib_ptlo1:ib_pthi1,
     &           0:nb_ptcomp-1)
      integer iBmag_ptlo0,iBmag_ptlo1
      integer iBmag_pthi0,iBmag_pthi1
      REAL_T Bmag_pt(
     &           iBmag_ptlo0:iBmag_pthi0,
     &           iBmag_ptlo1:iBmag_pthi1)
      integer nbunit_ptcomp
      integer ibunit_ptlo0,ibunit_ptlo1
      integer ibunit_pthi0,ibunit_pthi1
      REAL_T bunit_pt(
     &           ibunit_ptlo0:ibunit_pthi0,
     &           ibunit_ptlo1:ibunit_pthi1,
     &           0:nbunit_ptcomp-1)
      integer ngradb_ptcomp
      integer igradb_ptlo0,igradb_ptlo1
      integer igradb_pthi0,igradb_pthi1
      REAL_T gradb_pt(
     &           igradb_ptlo0:igradb_pthi0,
     &           igradb_ptlo1:igradb_pthi1,
     &           0:ngradb_ptcomp-1)
      integer ncurlb_ptcomp
      integer icurlb_ptlo0,icurlb_ptlo1
      integer icurlb_pthi0,icurlb_pthi1
      REAL_T curlb_pt(
     &           icurlb_ptlo0:icurlb_pthi0,
     &           icurlb_ptlo1:icurlb_pthi1,
     &           0:ncurlb_ptcomp-1)
      integer ibdotcurlb_ptlo0,ibdotcurlb_ptlo1
      integer ibdotcurlb_pthi0,ibdotcurlb_pthi1
      REAL_T bdotcurlb_pt(
     &           ibdotcurlb_ptlo0:ibdotcurlb_pthi0,
     &           ibdotcurlb_ptlo1:ibdotcurlb_pthi1)
      integer i0,i1, l
      double precision
     &     xi, Bz, By, BzBy, dBydx, dBzdx, Bmag, bunitDotcurlb,
     &     bunit(0:2), Bvec(0:2), gradB(0:2), curlb(0:2)
      dBydx = (ByOuter - ByInner)/xmax
      dBzdx = (BzOuter - BzInner)/xmax
      
      do i1 = igridboxlo1,igridboxhi1
      do i0 = igridboxlo0,igridboxhi0

         xi = xix(i0,i1)
         By = ByInner + dBydx * xi * (xmax/ximax) 
         Bz = BzInner + dBzdx * xi * (xmax/ximax)
         Bvec(0) = 0.
         Bvec(1) = By
         Bvec(2) = Bz
         Bmag = dsqrt(By*By + Bz*Bz)
         bunit(0) = 0.
         bunit(1) = By/Bmag
         bunit(2) = Bz/Bmag
         gradB(0) = (Bz*dBzdx + By*dBydx)/Bmag
         gradB(1) = zero
         gradB(2) = zero 
         curlb(0) = zero
         curlb(1) = -(By**2 * dBzdx - By*Bz*dBydx)/(Bmag**3)         
         curlb(2) =  (Bz**2 * dBydx - Bz*By*dBzdx)/(Bmag**3)         
         if (curlb(1)**2 + curlb(2)**2 .gt. 1.0e-10) then
           bunitDotcurlb = bunit(1) * curlb(1) + bunit(2) * curlb(2) 
           bunitDotcurlb = bunitDotcurlb / dsqrt(curlb(1)**2 + curlb(2)**2)
	 else
           bunitDotcurlb = zero
         endif
         Bmag_pt(i0,i1) = Bmag
         bdotcurlb_pt(i0,i1) = bunitDotcurlb
         do l = 0, 2
            bunit_pt(i0,i1,l) = bunit(l)
            b_pt(i0,i1,l)     = Bvec(l)
            gradb_pt(i0,i1,l) = gradB(l)
            curlb_pt(i0,i1,l) = curlb(l)
         end do
      
      enddo
      enddo
      return
      end
