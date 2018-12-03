#include "REAL.H"
#include "SPACE.H"
#include "CONSTANTS.H"

#include "CONSTANTS.H"
      subroutine COMPUTE_LIMITS_3D(
     &           ibeg
     &           ,iend
     &           ,istride
     &           ,ibdryboxlo0,ibdryboxlo1,ibdryboxlo2
     &           ,ibdryboxhi0,ibdryboxhi1,ibdryboxhi2
     &           ,idir
     &           ,iidirhi0
     &           ,iside
     &           ,iisidehi0
     &           ,icodim
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ibeg(0:2)
      integer iend(0:2)
      integer istride(0:2)
      integer ibdryboxlo0,ibdryboxlo1,ibdryboxlo2
      integer ibdryboxhi0,ibdryboxhi1,ibdryboxhi2
      integer iidirhi0
      integer idir(
     &           0:iidirhi0)
      integer iisidehi0
      integer iside(
     &           0:iisidehi0)
      integer icodim
      integer ic,m,itmp
      
        ibeg(0) = ibdryboxlo0
        ibeg(1) = ibdryboxlo1
        ibeg(2) = ibdryboxlo2
      
        iend(0) = ibdryboxhi0
        iend(1) = ibdryboxhi1
        iend(2) = ibdryboxhi2
      do m=0,CH_SPACEDIM-1
        istride(m) = 1
      enddo
      do ic=0,iidirhi0
        if ((iside(ic).eq.0)) then
          itmp = ibeg(idir(ic))
          ibeg(idir(ic))    = iend(idir(ic))
          iend(idir(ic))    = itmp
          istride(idir(ic)) = -1
        endif
      enddo
      return
      end
      subroutine FILL_CODIM_GHOST_CELLS_3D(
     &           f
     &           ,iflo0,iflo1,iflo2
     &           ,ifhi0,ifhi1,ifhi2
     &           ,nfcomp
     &           ,ibdryboxlo0,ibdryboxlo1,ibdryboxlo2
     &           ,ibdryboxhi0,ibdryboxhi1,ibdryboxhi2
     &           ,idir
     &           ,iidirhi0
     &           ,iside
     &           ,iisidehi0
     &           ,icodim
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nfcomp
      integer iflo0,iflo1,iflo2
      integer ifhi0,ifhi1,ifhi2
      REAL_T f(
     &           iflo0:ifhi0,
     &           iflo1:ifhi1,
     &           iflo2:ifhi2,
     &           0:nfcomp-1)
      integer ibdryboxlo0,ibdryboxlo1,ibdryboxlo2
      integer ibdryboxhi0,ibdryboxhi1,ibdryboxhi2
      integer iidirhi0
      integer idir(
     &           0:iidirhi0)
      integer iisidehi0
      integer iside(
     &           0:iisidehi0)
      integer icodim
      if (icodim.eq.2) then
      call FILL_CODIM2_GHOST_CELLS_3D(
     &           f
     &           ,iflo0,iflo1,iflo2
     &           ,ifhi0,ifhi1,ifhi2
     &           ,nfcomp
     &           ,ibdryboxlo0,ibdryboxlo1,ibdryboxlo2
     &           ,ibdryboxhi0,ibdryboxhi1,ibdryboxhi2
     &           ,idir
     &           ,iidirhi0
     &           ,iside
     &           ,iisidehi0
     &           )
      else if (icodim.eq.3) then
      call FILL_CODIM3_GHOST_CELLS_3D(
     &           f
     &           ,iflo0,iflo1,iflo2
     &           ,ifhi0,ifhi1,ifhi2
     &           ,nfcomp
     &           ,ibdryboxlo0,ibdryboxlo1,ibdryboxlo2
     &           ,ibdryboxhi0,ibdryboxhi1,ibdryboxhi2
     &           ,idir
     &           ,iidirhi0
     &           ,iside
     &           ,iisidehi0
     &           )
      else if (icodim.eq.4) then
      call FILL_CODIM4_GHOST_CELLS_3D(
     &           f
     &           ,iflo0,iflo1,iflo2
     &           ,ifhi0,ifhi1,ifhi2
     &           ,nfcomp
     &           ,ibdryboxlo0,ibdryboxlo1,ibdryboxlo2
     &           ,ibdryboxhi0,ibdryboxhi1,ibdryboxhi2
     &           ,idir
     &           ,iidirhi0
     &           ,iside
     &           ,iisidehi0
     &           )
      endif
      return
      end
      subroutine FILL_CODIM2_GHOST_CELLS_3D(
     &           f
     &           ,iflo0,iflo1,iflo2
     &           ,ifhi0,ifhi1,ifhi2
     &           ,nfcomp
     &           ,ibdryboxlo0,ibdryboxlo1,ibdryboxlo2
     &           ,ibdryboxhi0,ibdryboxhi1,ibdryboxhi2
     &           ,idir
     &           ,iidirhi0
     &           ,iside
     &           ,iisidehi0
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer nfcomp
      integer iflo0,iflo1,iflo2
      integer ifhi0,ifhi1,ifhi2
      REAL_T f(
     &           iflo0:ifhi0,
     &           iflo1:ifhi1,
     &           iflo2:ifhi2,
     &           0:nfcomp-1)
      integer ibdryboxlo0,ibdryboxlo1,ibdryboxlo2
      integer ibdryboxhi0,ibdryboxhi1,ibdryboxhi2
      integer iidirhi0
      integer idir(
     &           0:iidirhi0)
      integer iisidehi0
      integer iside(
     &           0:iisidehi0)
      integer i,j,k
      integer i10,j10,k10
      integer i01,j01,k01
      integer i11,j11,k11
      integer ibeg(0:CH_SPACEDIM-1)
      integer iend(0:CH_SPACEDIM-1)
      integer istride(0:CH_SPACEDIM-1)
      integer isign(0:1)
      integer ic,icodim,n
      icodim = 2
      do ic=0,icodim-1
        isign(ic) = 1-2*iside(ic)
      enddo
      call COMPUTE_LIMITS_3D(
     &           ibeg
     &           ,iend
     &           ,istride
     &           ,ibdryboxlo0,ibdryboxlo1,ibdryboxlo2
     &           ,ibdryboxhi0,ibdryboxhi1,ibdryboxhi2
     &           ,idir
     &           ,iidirhi0
     &           ,iside
     &           ,iisidehi0
     &           ,icodim
     &           )
      do n=0,nfcomp-1
      
      do k=ibeg(2),iend(2),istride(2)
      do j=ibeg(1),iend(1),istride(1)
      do i=ibeg(0),iend(0),istride(0)
        
        i10 = i+CHF_ID(idir(0),0)*isign(0)
        j10 = j+CHF_ID(idir(0),1)*isign(0)
        k10 = k+CHF_ID(idir(0),2)*isign(0)
        
        i01 = i+CHF_ID(idir(1),0)*isign(1)
        j01 = j+CHF_ID(idir(1),1)*isign(1)
        k01 = k+CHF_ID(idir(1),2)*isign(1)
        
        i11 = i10+i01-i
        j11 = j10+j01-j
        k11 = k10+k01-k
        f(i,j,k,n) = f(i10,j10,k10,n)
     &                       + f(i01,j01,k01,n)
     &                       - f(i11,j11,k11,n)
      
      enddo
      enddo
      enddo
      enddo
      return
      end
      subroutine FILL_CODIM3_GHOST_CELLS_3D(
     &           f
     &           ,iflo0,iflo1,iflo2
     &           ,ifhi0,ifhi1,ifhi2
     &           ,nfcomp
     &           ,ibdryboxlo0,ibdryboxlo1,ibdryboxlo2
     &           ,ibdryboxhi0,ibdryboxhi1,ibdryboxhi2
     &           ,idir
     &           ,iidirhi0
     &           ,iside
     &           ,iisidehi0
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer nfcomp
      integer iflo0,iflo1,iflo2
      integer ifhi0,ifhi1,ifhi2
      REAL_T f(
     &           iflo0:ifhi0,
     &           iflo1:ifhi1,
     &           iflo2:ifhi2,
     &           0:nfcomp-1)
      integer ibdryboxlo0,ibdryboxlo1,ibdryboxlo2
      integer ibdryboxhi0,ibdryboxhi1,ibdryboxhi2
      integer iidirhi0
      integer idir(
     &           0:iidirhi0)
      integer iisidehi0
      integer iside(
     &           0:iisidehi0)
      integer i,j,k
      integer i100,j100,k100
      integer i010,j010,k010
      integer i001,j001,k001
      integer i111,j111,k111
      integer ibeg(0:CH_SPACEDIM-1)
      integer iend(0:CH_SPACEDIM-1)
      integer istride(0:CH_SPACEDIM-1)
      integer isign(0:2)
      integer ic,icodim,n
      icodim = 3
      do ic=0,icodim-1
        isign(ic) = 1-2*iside(ic)
      enddo
      call COMPUTE_LIMITS_3D(
     &           ibeg
     &           ,iend
     &           ,istride
     &           ,ibdryboxlo0,ibdryboxlo1,ibdryboxlo2
     &           ,ibdryboxhi0,ibdryboxhi1,ibdryboxhi2
     &           ,idir
     &           ,iidirhi0
     &           ,iside
     &           ,iisidehi0
     &           ,icodim
     &           )
      do n=0,nfcomp-1
      
      do k=ibeg(2),iend(2),istride(2)
      do j=ibeg(1),iend(1),istride(1)
      do i=ibeg(0),iend(0),istride(0)
        
        i100 = i+CHF_ID(idir(0),0)*isign(0)
        j100 = j+CHF_ID(idir(0),1)*isign(0)
        k100 = k+CHF_ID(idir(0),2)*isign(0)
        
        i010 = i+CHF_ID(idir(1),0)*isign(1)
        j010 = j+CHF_ID(idir(1),1)*isign(1)
        k010 = k+CHF_ID(idir(1),2)*isign(1)
        
        i001 = i+CHF_ID(idir(2),0)*isign(2)
        j001 = j+CHF_ID(idir(2),1)*isign(2)
        k001 = k+CHF_ID(idir(2),2)*isign(2)
        
        i111 = i100+i010+i001-2*i
        j111 = j100+j010+j001-2*j
        k111 = k100+k010+k001-2*k
        f(i,j,k,n) = f(i100,j100,k100,n)
     &                       + f(i010,j010,k010,n)
     &                       + f(i001,j001,k001,n)
     &                 - two * f(i111,j111,k111,n)
      
      enddo
      enddo
      enddo
      enddo
      return
      end
      subroutine FILL_CODIM4_GHOST_CELLS_3D(
     &           f
     &           ,iflo0,iflo1,iflo2
     &           ,ifhi0,ifhi1,ifhi2
     &           ,nfcomp
     &           ,ibdryboxlo0,ibdryboxlo1,ibdryboxlo2
     &           ,ibdryboxhi0,ibdryboxhi1,ibdryboxhi2
     &           ,idir
     &           ,iidirhi0
     &           ,iside
     &           ,iisidehi0
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer nfcomp
      integer iflo0,iflo1,iflo2
      integer ifhi0,ifhi1,ifhi2
      REAL_T f(
     &           iflo0:ifhi0,
     &           iflo1:ifhi1,
     &           iflo2:ifhi2,
     &           0:nfcomp-1)
      integer ibdryboxlo0,ibdryboxlo1,ibdryboxlo2
      integer ibdryboxhi0,ibdryboxhi1,ibdryboxhi2
      integer iidirhi0
      integer idir(
     &           0:iidirhi0)
      integer iisidehi0
      integer iside(
     &           0:iisidehi0)
      integer i,j,k
      integer i1000,j1000,k1000
      integer i0100,j0100,k0100
      integer i0010,j0010,k0010
      integer i0001,j0001,k0001
      integer i1111,j1111,k1111
      integer ibeg(0:CH_SPACEDIM-1)
      integer iend(0:CH_SPACEDIM-1)
      integer istride(0:CH_SPACEDIM-1)
      integer isign(0:3)
      integer ic,icodim,n
      icodim = 4
      do ic=0,icodim-1
        isign(ic) = 1-2*iside(ic)
      enddo
      call COMPUTE_LIMITS_3D(
     &           ibeg
     &           ,iend
     &           ,istride
     &           ,ibdryboxlo0,ibdryboxlo1,ibdryboxlo2
     &           ,ibdryboxhi0,ibdryboxhi1,ibdryboxhi2
     &           ,idir
     &           ,iidirhi0
     &           ,iside
     &           ,iisidehi0
     &           ,icodim
     &           )
      do n=0,nfcomp-1
      
      do k=ibeg(2),iend(2),istride(2)
      do j=ibeg(1),iend(1),istride(1)
      do i=ibeg(0),iend(0),istride(0)
        
        i1000 = i+CHF_ID(idir(0),0)*isign(0)
        j1000 = j+CHF_ID(idir(0),1)*isign(0)
        k1000 = k+CHF_ID(idir(0),2)*isign(0)
        
        i0100 = i+CHF_ID(idir(1),0)*isign(1)
        j0100 = j+CHF_ID(idir(1),1)*isign(1)
        k0100 = k+CHF_ID(idir(1),2)*isign(1)
        
        i0010 = i+CHF_ID(idir(2),0)*isign(2)
        j0010 = j+CHF_ID(idir(2),1)*isign(2)
        k0010 = k+CHF_ID(idir(2),2)*isign(2)
        
        i0001 = i+CHF_ID(idir(3),0)*isign(3)
        j0001 = j+CHF_ID(idir(3),1)*isign(3)
        k0001 = k+CHF_ID(idir(3),2)*isign(3)
        
        i1111 = i1000+i0100+i0010+i0001-3*i
        j1111 = j1000+j0100+j0010+j0001-3*j
        k1111 = k1000+k0100+k0010+k0001-3*k
        f(i,j,k,n) = f(i1000,j1000,k1000,n)
     &                       + f(i0100,j0100,k0100,n)
     &                       + f(i0010,j0010,k0010,n)
     &                       + f(i0001,j0001,k0001,n)
     &                - four * f(i1111,j1111,k1111,n)
      
      enddo
      enddo
      enddo
      enddo
      return
      end
