#include "REAL.H"
#include "SPACE.H"
#include "CONSTANTS.H"

      subroutine CELL_CENTERED_FIELD_COMPONENT_2D(
     &           iboxlo0,iboxlo1
     &           ,iboxhi0,iboxhi1
     &           ,dir
     &           ,phi
     &           ,iphilo0,iphilo1
     &           ,iphihi0,iphihi1
     &           ,h
     &           ,order
     &           ,Efield
     &           ,iEfieldlo0,iEfieldlo1
     &           ,iEfieldhi0,iEfieldhi1
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer dir
      integer iphilo0,iphilo1
      integer iphihi0,iphihi1
      REAL_T phi(
     &           iphilo0:iphihi0,
     &           iphilo1:iphihi1)
      REAL_T h(0:1)
      integer order
      integer iEfieldlo0,iEfieldlo1
      integer iEfieldhi0,iEfieldhi1
      REAL_T Efield(
     &           iEfieldlo0:iEfieldhi0,
     &           iEfieldlo1:iEfieldhi1)
      integer i,j, ii,jj
      
      ii = CHF_ID(0,dir)
      jj = CHF_ID(1,dir)
      if (order .eq. 4) then
         
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0

            Efield(i,j) = -(
     &        phi(i-2*ii,j-2*jj)
     &        - eight * phi(i-  ii,j-  jj)
     &        + eight * phi(i+  ii,j+  jj)
     &        - phi(i+2*ii,j+2*jj)
     &        ) / (twelve * h(dir))
         
      enddo
      enddo
      else if (order .eq. 2) then
         
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0

            Efield(i,j) = -(
     &        - phi(i-  ii,j-  jj)
     &        + phi(i+  ii,j+  jj)
     &        ) / (two * h(dir))
         
      enddo
      enddo
      endif
      return
      end
      subroutine FACE_CENTERED_FIELD_COMPONENT_2D(
     &           iboxlo0,iboxlo1
     &           ,iboxhi0,iboxhi1
     &           ,dir
     &           ,phi
     &           ,iphilo0,iphilo1
     &           ,iphihi0,iphihi1
     &           ,h
     &           ,order
     &           ,Efield
     &           ,iEfieldlo0,iEfieldlo1
     &           ,iEfieldhi0,iEfieldhi1
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer dir
      integer iphilo0,iphilo1
      integer iphihi0,iphihi1
      REAL_T phi(
     &           iphilo0:iphihi0,
     &           iphilo1:iphihi1)
      REAL_T h(0:1)
      integer order
      integer iEfieldlo0,iEfieldlo1
      integer iEfieldhi0,iEfieldhi1
      REAL_T Efield(
     &           iEfieldlo0:iEfieldhi0,
     &           iEfieldlo1:iEfieldhi1)
      integer i,j, ii,jj
      
      ii = CHF_ID(0,dir)
      jj = CHF_ID(1,dir)
      if (order .eq. 4) then
         
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0

            Efield(i,j) = -(
     &        27.d0 * (phi(i     ,j     )
     &               - phi(i-  ii,j-  jj))
     &               - phi(i+  ii,j+  jj)
     &               + phi(i-2*ii,j-2*jj)
     &        ) / (24.d0 * h(dir))
         
      enddo
      enddo
      else if (order .eq. 2) then
         
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0

            Efield(i,j) = -(
     &        - phi(i-  ii,j-  jj)
     &        + phi(i     ,j     )
     &        ) / h(dir)
         
      enddo
      enddo
      endif
      return
      end
      subroutine FACE_INTERPOLATE_2D(
     &           dir
     &           ,iboxlo0,iboxlo1
     &           ,iboxhi0,iboxhi1
     &           ,order
     &           ,phi
     &           ,iphilo0,iphilo1
     &           ,iphihi0,iphihi1
     &           ,data
     &           ,idatalo0,idatalo1
     &           ,idatahi0,idatahi1
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer dir
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer order
      integer iphilo0,iphilo1
      integer iphihi0,iphihi1
      REAL_T phi(
     &           iphilo0:iphihi0,
     &           iphilo1:iphihi1)
      integer idatalo0,idatalo1
      integer idatahi0,idatahi1
      REAL_T data(
     &           idatalo0:idatahi0,
     &           idatalo1:idatahi1)
      integer i,j, ii,jj
      
      ii = CHF_ID(0,dir)
      jj = CHF_ID(1,dir)
      if (order .eq. 4) then
         
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0

            data(i,j) = (
     &        9.d0 * (phi(i     ,j     )
     &              + phi(i-  ii,j-  jj))
     &             - (phi(i+  ii,j+  jj)
     &              + phi(i-2*ii,j-2*jj))
     &        ) / 16.d0
         
      enddo
      enddo
      else if (order .eq. 2) then
         
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0

         data(i,j) = (
     &            phi(i-  ii,j-  jj)
     &          + phi(i     ,j     )
     &        ) / 2.d0
         
      enddo
      enddo
      endif
      return
      end
      subroutine EXTRAP_FOR_CC_OPS_2D(
     &           dir
     &           ,side
     &           ,order
     &           ,ifaceboxlo0,ifaceboxlo1
     &           ,ifaceboxhi0,ifaceboxhi1
     &           ,iinteriorboxlo0,iinteriorboxlo1
     &           ,iinteriorboxhi0,iinteriorboxhi1
     &           ,array
     &           ,iarraylo0,iarraylo1
     &           ,iarrayhi0,iarrayhi1
     &           ,narraycomp
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer dir
      integer side
      integer order
      integer ifaceboxlo0,ifaceboxlo1
      integer ifaceboxhi0,ifaceboxhi1
      integer iinteriorboxlo0,iinteriorboxlo1
      integer iinteriorboxhi0,iinteriorboxhi1
      integer narraycomp
      integer iarraylo0,iarraylo1
      integer iarrayhi0,iarrayhi1
      REAL_T array(
     &           iarraylo0:iarrayhi0,
     &           iarraylo1:iarrayhi1,
     &           0:narraycomp-1)
      integer i,id,ni,j,jd,nj, m, n, comp, ncomp
      double precision sum, coef2(0:2,2), coef4(0:4,3)
      data coef2 /  3.d0, -3.d0, 1.d0,
     &             6.d0, -8.d0, 3.d0 /
      data coef4 /  5.d0,  -10.d0,  10.d0,  -5.d0,  1.d0,
     &            15.d0,  -40.d0,  45.d0, -24.d0,  5.d0,
     &            35.d0, -105.d0, 126.d0, -70.d0, 15.d0 /
      ncomp = narraycomp
      
      id = CHF_ID(0,dir)*side
      jd = CHF_ID(1,dir)*side
      
      do j = ifaceboxlo1,ifaceboxhi1
      do i = ifaceboxlo0,ifaceboxhi0

        if (side .eq. -1) then
           
           ni = id*(i-iinteriorboxlo0)
           n = ni
           nj = jd*(j-iinteriorboxlo1)
           n = n + nj
         else if (side .eq. 1) then
            
            ni = id*(i-iinteriorboxhi0)
            n = ni
            nj = jd*(j-iinteriorboxhi1)
            n = n + nj
          endif
          do comp = 0, ncomp-1
             sum = zero
             if (order .eq. 4) then
                do m = 0, 4
                   sum = sum + coef4(m,n)*array(i-id*(ni+m),j-jd*(nj+m),comp)
                enddo
                array(i,j,comp) = sum
             else if (order .eq. 2) then
                do m = 0, 2
                   sum = sum + coef2(m,n)*array(i-id*(ni+m),j-jd*(nj+m),comp)
                enddo
                array(i,j,comp) = sum
             endif
          enddo
      
      enddo
      enddo
      return
      end
      subroutine EXTRAP_FOR_FC_OPS_2D(
     &           dir
     &           ,side
     &           ,order
     &           ,ifaceboxlo0,ifaceboxlo1
     &           ,ifaceboxhi0,ifaceboxhi1
     &           ,iinteriorboxlo0,iinteriorboxlo1
     &           ,iinteriorboxhi0,iinteriorboxhi1
     &           ,array
     &           ,iarraylo0,iarraylo1
     &           ,iarrayhi0,iarrayhi1
     &           ,narraycomp
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer dir
      integer side
      integer order
      integer ifaceboxlo0,ifaceboxlo1
      integer ifaceboxhi0,ifaceboxhi1
      integer iinteriorboxlo0,iinteriorboxlo1
      integer iinteriorboxhi0,iinteriorboxhi1
      integer narraycomp
      integer iarraylo0,iarraylo1
      integer iarrayhi0,iarrayhi1
      REAL_T array(
     &           iarraylo0:iarrayhi0,
     &           iarraylo1:iarrayhi1,
     &           0:narraycomp-1)
      integer i,id,ni,j,jd,nj, m, n, comp, ncomp
      double precision sum, coef2(0:2,2), coef4(0:4,2)
      data coef2 /  3.d0, -3.d0, 1.d0,
     &             6.d0, -8.d0, 3.d0 /
      data coef4 /  4.625d0, -8.5d0, 7.75d0, -3.5d0, 0.625d0,
     &             11.25d0, -25.d0, 22.5d0,  -9.d0, 1.25d0 /
      ncomp = narraycomp
      
      id = CHF_ID(0,dir)*side
      jd = CHF_ID(1,dir)*side
      
      do j = ifaceboxlo1,ifaceboxhi1
      do i = ifaceboxlo0,ifaceboxhi0

        if (side .eq. -1) then
           
           ni = id*(i-iinteriorboxlo0)
           n = ni
           nj = jd*(j-iinteriorboxlo1)
           n = n + nj
         else if (side .eq. 1) then
            
            ni = id*(i-iinteriorboxhi0)
            n = ni
            nj = jd*(j-iinteriorboxhi1)
            n = n + nj
          endif
          do comp = 0, ncomp-1
             sum = zero
             if (order .eq. 4) then
                do m = 0, 4
                   sum = sum + coef4(m,n)*array(i-id*(ni+m),j-jd*(nj+m),comp)
                enddo
                array(i,j,comp) = sum
             else if (order .eq. 2) then
                do m = 0, 2
                   sum = sum + coef2(m,n)*array(i-id*(ni+m),j-jd*(nj+m),comp)
                enddo
                array(i,j,comp) = sum
             endif
          enddo
      
      enddo
      enddo
      return
      end
