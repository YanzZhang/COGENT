#include "REAL.H"
#include "SPACE.H"
#include "CONSTANTS.H"

      DOUBLE PRECISION FUNCTION minmod(a,b)
        DOUBLE PRECISION, INTENT(IN) :: a,b
        IF ((a .LT. 0) .AND. (b .LT. 0)) THEN
          minmod = MAX(a,b)
        ELSEIF ((a .GT. 0) .AND. (b .GT. 0)) THEN
          minmod = MIN(a,b)
        ELSE                                  
          minmod = 0.0_8
        ENDIF
        RETURN
      END FUNCTION 
      DOUBLE PRECISION FUNCTION mc(a,b)
        DOUBLE PRECISION, INTENT(IN)  :: a,b
        DOUBLE PRECISION              :: minmod
        mc = minmod(2.0_8*minmod(a,b),0.5_8*(a+b))
        RETURN
      END FUNCTION 
      DOUBLE PRECISION FUNCTION med(a,b,c)
        DOUBLE PRECISION, INTENT(IN)  :: a,b,c
        med = a + b + c - MAX(a,b,c) - MIN(a,b,c)
        RETURN
      END FUNCTION 
      SUBROUTINE RECONSTRUCT_DFN_FACE_VPAR_4D(
     &           f
     &           ,iflo0,iflo1,iflo2,iflo3
     &           ,ifhi0,ifhi1,ifhi2,ifhi3
     &           ,a
     &           ,ialo0,ialo1,ialo2,ialo3
     &           ,iahi0,iahi1,iahi2,iahi3
     &           ,u
     &           ,iulo0,iulo1,iulo2,iulo3
     &           ,iuhi0,iuhi1,iuhi2,iuhi3
     &           ,igridlo0,igridlo1,igridlo2,igridlo3
     &           ,igridhi0,igridhi1,igridhi2,igridhi3
     &           ,flag
     &           ,Nvpar
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iflo0,iflo1,iflo2,iflo3
      integer ifhi0,ifhi1,ifhi2,ifhi3
      REAL_T f(
     &           iflo0:ifhi0,
     &           iflo1:ifhi1,
     &           iflo2:ifhi2,
     &           iflo3:ifhi3)
      integer ialo0,ialo1,ialo2,ialo3
      integer iahi0,iahi1,iahi2,iahi3
      REAL_T a(
     &           ialo0:iahi0,
     &           ialo1:iahi1,
     &           ialo2:iahi2,
     &           ialo3:iahi3)
      integer iulo0,iulo1,iulo2,iulo3
      integer iuhi0,iuhi1,iuhi2,iuhi3
      REAL_T u(
     &           iulo0:iuhi0,
     &           iulo1:iuhi1,
     &           iulo2:iuhi2,
     &           iulo3:iuhi3)
      integer igridlo0,igridlo1,igridlo2,igridlo3
      integer igridhi0,igridhi1,igridhi2,igridhi3
      integer flag
      integer Nvpar
        INTEGER                     :: i,j,k,l
        DOUBLE PRECISION            :: fc, fu, fd, f3, fq, med
        DOUBLE PRECISION, PARAMETER :: flc_fac = 0.25, one_by_eight = 0.125,
     &                                 one_by_thirty        = 1.0_8/30.0_8,
     &                                 thirteen_by_sixty    = 13.0_8/60.0_8,
     &                                 fortyseven_by_sixty  = 47.0_8/60.0_8,
     &                                 twentyseven_by_sixty = 27.0_8/60.0_8,
     &                                 one_by_twenty        = 1.0_8/20.0_8
        IF (flag .EQ. 1) THEN
#if CH_SPACEDIM==5
          
      do k = igridlo3,igridhi3
      do m = igridlo2,igridhi2
      do j = igridlo1,igridhi1
      do i = igridlo0,igridhi0

            fc = 0.5_8 * (u(i,j,m,k-1) + u(i,j,m,k))
            IF (a(i,j,m,k) .LT. 0.0_8) THEN
              fu = u(i,j,m,k-1)
              fd = u(i,j,m,k)
              fq = fc - one_by_eight * (u(i,j,m,k-2)+u(i,j,m,k)-2.0_8*u(i,j,m,k-1))
            ELSE
              fu = u(i,j,m,k)
              fd = u(i,j,m,k-1)
              fq = fc - one_by_eight * (u(i,j,m,k+1)+u(i,j,m,k-1)-2.0_8*u(i,j,m,k))
            ENDIF
            f3 = fu + flc_fac * (fc - fu)
            f(i,j,m,k) = med(fc,med(fc,fd,f3),fq)
          
      enddo
      enddo
      enddo
      enddo
#else
          
      do l = igridlo3,igridhi3
      do k = igridlo2,igridhi2
      do j = igridlo1,igridhi1
      do i = igridlo0,igridhi0

            fc = 0.5_8 * (u(i,j,k-1,l) + u(i,j,k,l))
            IF (a(i,j,k,l) .LT. 0.0_8) THEN
              fu = u(i,j,k-1,l)
              fd = u(i,j,k,l)
              fq = fc - one_by_eight * (u(i,j,k-2,l)+u(i,j,k,l)-2.0_8*u(i,j,k-1,l))
            ELSE
              fu = u(i,j,k,l)
              fd = u(i,j,k-1,l)
              fq = fc - one_by_eight * (u(i,j,k+1,l)+u(i,j,k-1,l)-2.0_8*u(i,j,k,l))
            ENDIF
            f3 = fu + flc_fac * (fc - fu)
            f(i,j,k,l) = med(fc,med(fc,fd,f3),fq)
          
      enddo
      enddo
      enddo
      enddo
#endif
        ELSE
#if CH_SPACEDIM==5
          
      do k = igridlo3,igridhi3
      do m = igridlo2,igridhi2
      do j = igridlo1,igridhi1
      do i = igridlo0,igridhi0

            IF (a(i,j,m,k) .LT. 0.0_8) THEN
              f(i,j,m,k) =   one_by_thirty        * u(i,j,m,k-3)
     &                               - thirteen_by_sixty    * u(i,j,m,k-2)
     &                               + fortyseven_by_sixty  * u(i,j,m,k-1)
     &                               + twentyseven_by_sixty * u(i,j,m,k  )
     &                               - one_by_twenty        * u(i,j,m,k+1)
            ELSE
              f(i,j,m,k) =   one_by_thirty        * u(i,j,m,k+2)
     &                               - thirteen_by_sixty    * u(i,j,m,k+1)
     &                               + fortyseven_by_sixty  * u(i,j,m,k  )
     &                               + twentyseven_by_sixty * u(i,j,m,k-1)
     &                               - one_by_twenty        * u(i,j,m,k-2)
            ENDIF
          
      enddo
      enddo
      enddo
      enddo
#else
          
      do l = igridlo3,igridhi3
      do k = igridlo2,igridhi2
      do j = igridlo1,igridhi1
      do i = igridlo0,igridhi0

            IF (a(i,j,k,l) .LT. 0.0_8) THEN
              f(i,j,k,l) =   one_by_thirty        * u(i,j,k-3,l)
     &                             - thirteen_by_sixty    * u(i,j,k-2,l)
     &                             + fortyseven_by_sixty  * u(i,j,k-1,l)
     &                             + twentyseven_by_sixty * u(i,j,k  ,l)
     &                             - one_by_twenty        * u(i,j,k+1,l)
            ELSE
              f(i,j,k,l) =   one_by_thirty        * u(i,j,k+2,l)
     &                             - thirteen_by_sixty    * u(i,j,k+1,l)
     &                             + fortyseven_by_sixty  * u(i,j,k  ,l)
     &                             + twentyseven_by_sixty * u(i,j,k-1,l)
     &                             - one_by_twenty        * u(i,j,k-2,l)
            ENDIF
          
      enddo
      enddo
      enddo
      enddo
#endif
        ENDIF
      END SUBROUTINE 
      SUBROUTINE RECONSTRUCT_DFN_FACE_MU_4D(
     &           f
     &           ,iflo0,iflo1,iflo2,iflo3
     &           ,ifhi0,ifhi1,ifhi2,ifhi3
     &           ,a
     &           ,ialo0,ialo1,ialo2,ialo3
     &           ,iahi0,iahi1,iahi2,iahi3
     &           ,u
     &           ,iulo0,iulo1,iulo2,iulo3
     &           ,iuhi0,iuhi1,iuhi2,iuhi3
     &           ,igridlo0,igridlo1,igridlo2,igridlo3
     &           ,igridhi0,igridhi1,igridhi2,igridhi3
     &           ,flag
     &           ,Nmu
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iflo0,iflo1,iflo2,iflo3
      integer ifhi0,ifhi1,ifhi2,ifhi3
      REAL_T f(
     &           iflo0:ifhi0,
     &           iflo1:ifhi1,
     &           iflo2:ifhi2,
     &           iflo3:ifhi3)
      integer ialo0,ialo1,ialo2,ialo3
      integer iahi0,iahi1,iahi2,iahi3
      REAL_T a(
     &           ialo0:iahi0,
     &           ialo1:iahi1,
     &           ialo2:iahi2,
     &           ialo3:iahi3)
      integer iulo0,iulo1,iulo2,iulo3
      integer iuhi0,iuhi1,iuhi2,iuhi3
      REAL_T u(
     &           iulo0:iuhi0,
     &           iulo1:iuhi1,
     &           iulo2:iuhi2,
     &           iulo3:iuhi3)
      integer igridlo0,igridlo1,igridlo2,igridlo3
      integer igridhi0,igridhi1,igridhi2,igridhi3
      integer flag
      integer Nmu
        INTEGER                     :: i,j,k,l
        DOUBLE PRECISION            :: fc, fu, fd, f3, fq, med
        DOUBLE PRECISION, PARAMETER :: flc_fac = 0.25, one_by_eight = 0.125,
     &                                 one_by_thirty        = 1.0_8/30.0_8,
     &                                 thirteen_by_sixty    = 13.0_8/60.0_8,
     &                                 fortyseven_by_sixty  = 47.0_8/60.0_8,
     &                                 twentyseven_by_sixty = 27.0_8/60.0_8,
     &                                 one_by_twenty        = 1.0_8/20.0_8
        IF (flag .EQ. 1) THEN
#if CH_SPACEDIM==5
          
      do k = igridlo3,igridhi3
      do m = igridlo2,igridhi2
      do j = igridlo1,igridhi1
      do i = igridlo0,igridhi0

            fc = 0.5_8 * (u(i,j,m,k) + u(i,j,m,k))
            IF (a(i,j,m,k) .LT. 0.0_8) THEN
              fu = u(i,j,m,k)
              fd = u(i,j,m,k)
              fq = fc - one_by_eight * (u(i,j,m,k)+u(i,j,m,k)-2.0_8*u(i,j,m,k))
            ELSE
              fu = u(i,j,m,k)
              fd = u(i,j,m,k)
              fq = fc - one_by_eight * (u(i,j,m,k)+u(i,j,m,k)-2.0_8*u(i,j,m,k))
            ENDIF
            f3 = fu + flc_fac * (fc - fu)
            f(i,j,m,k) = med(fc,med(fc,fd,f3),fq)
          
      enddo
      enddo
      enddo
      enddo
#else
          
      do l = igridlo3,igridhi3
      do k = igridlo2,igridhi2
      do j = igridlo1,igridhi1
      do i = igridlo0,igridhi0

            fc = 0.5_8 * (u(i,j,k,l-1) + u(i,j,k,l))
            IF (a(i,j,k,l) .LT. 0.0_8) THEN
              fu = u(i,j,k,l-1)
              fd = u(i,j,k,l)
              fq = fc - one_by_eight * (u(i,j,k,l-2)+u(i,j,k,l)-2.0_8*u(i,j,k,l-1))
            ELSE
              fu = u(i,j,k,l)
              fd = u(i,j,k,l-1)
              fq = fc - one_by_eight * (u(i,j,k,l+1)+u(i,j,k,l-1)-2.0_8*u(i,j,k,l))
            ENDIF
            f3 = fu + flc_fac * (fc - fu)
            f(i,j,k,l) = med(fc,med(fc,fd,f3),fq)
          
      enddo
      enddo
      enddo
      enddo
#endif
        ELSE
#if CH_SPACEDIM==5
          
      do k = igridlo3,igridhi3
      do m = igridlo2,igridhi2
      do j = igridlo1,igridhi1
      do i = igridlo0,igridhi0

            IF (a(i,j,m,k) .LT. 0.0_8) THEN
              f(i,j,m,k) =   one_by_thirty        * u(i,j,m,k)
     &                               - thirteen_by_sixty    * u(i,j,m,k)
     &                               + fortyseven_by_sixty  * u(i,j,m,k)
     &                               + twentyseven_by_sixty * u(i,j,m,k)
     &                               - one_by_twenty        * u(i,j,m,k)
            ELSE
              f(i,j,m,k) =   one_by_thirty        * u(i,j,m,k)
     &                               - thirteen_by_sixty    * u(i,j,m,k)
     &                               + fortyseven_by_sixty  * u(i,j,m,k)
     &                               + twentyseven_by_sixty * u(i,j,m,k)
     &                               - one_by_twenty        * u(i,j,m,k)
            ENDIF
          
      enddo
      enddo
      enddo
      enddo
#else
          
      do l = igridlo3,igridhi3
      do k = igridlo2,igridhi2
      do j = igridlo1,igridhi1
      do i = igridlo0,igridhi0

            IF (a(i,j,k,l) .LT. 0.0_8) THEN
              f(i,j,k,l) =   one_by_thirty        * u(i,j,k,l-3)
     &                             - thirteen_by_sixty    * u(i,j,k,l-2)
     &                             + fortyseven_by_sixty  * u(i,j,k,l-1)
     &                             + twentyseven_by_sixty * u(i,j,k,l  )
     &                             - one_by_twenty        * u(i,j,k,l+1)
            ELSE
              f(i,j,k,l) =   one_by_thirty        * u(i,j,k,l+2)
     &                             - thirteen_by_sixty    * u(i,j,k,l+1)
     &                             + fortyseven_by_sixty  * u(i,j,k,l  )
     &                             + twentyseven_by_sixty * u(i,j,k,l-1)
     &                             - one_by_twenty        * u(i,j,k,l-2)
            ENDIF
          
      enddo
      enddo
      enddo
      enddo
#endif
        ENDIF
      END SUBROUTINE 
