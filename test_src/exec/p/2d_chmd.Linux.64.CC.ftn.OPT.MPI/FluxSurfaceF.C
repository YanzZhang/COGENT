#include "REAL.H"
#include "SPACE.H"
#include "CONSTANTS.H"

      subroutine ADD_FLUX_SURFACE_ARRAY_2D(
     &           igridboxlo0,igridboxlo1
     &           ,igridboxhi0,igridboxhi1
     &           ,radial
     &           ,iradiallo0,iradiallo1
     &           ,iradialhi0,iradialhi1
     &           ,data
     &           ,idatalo0,idatalo1
     &           ,idatahi0,idatahi1
     &           ,ndatacomp
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer igridboxlo0,igridboxlo1
      integer igridboxhi0,igridboxhi1
      integer iradiallo0,iradiallo1
      integer iradialhi0,iradialhi1
      REAL_T radial(
     &           iradiallo0:iradialhi0,
     &           iradiallo1:iradialhi1)
      integer ndatacomp
      integer idatalo0,idatalo1
      integer idatahi0,idatahi1
      REAL_T data(
     &           idatalo0:idatahi0,
     &           idatalo1:idatahi1,
     &           0:ndatacomp-1)
      integer i0,i1, m, ncomp
      ncomp = ndatacomp
      do m = 0, ncomp-1
         
      do i1 = igridboxlo1,igridboxhi1
      do i0 = igridboxlo0,igridboxhi0

#if CH_SPACEDIM == 3
            data(i0,i1,m) = data(i0,i1,m)
     &           + radial(i0,iradiallo1,iradiallo2)
#else
            data(i0,i1,m) = data(i0,i1,m)
     &           + radial(i0,iradiallo1)
#endif
         
      enddo
      enddo
      enddo
      return
      end
      subroutine SUBTRACT_FLUX_SURFACE_ARRAY_2D(
     &           igridboxlo0,igridboxlo1
     &           ,igridboxhi0,igridboxhi1
     &           ,radial
     &           ,iradiallo0,iradiallo1
     &           ,iradialhi0,iradialhi1
     &           ,data
     &           ,idatalo0,idatalo1
     &           ,idatahi0,idatahi1
     &           ,ndatacomp
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer igridboxlo0,igridboxlo1
      integer igridboxhi0,igridboxhi1
      integer iradiallo0,iradiallo1
      integer iradialhi0,iradialhi1
      REAL_T radial(
     &           iradiallo0:iradialhi0,
     &           iradiallo1:iradialhi1)
      integer ndatacomp
      integer idatalo0,idatalo1
      integer idatahi0,idatahi1
      REAL_T data(
     &           idatalo0:idatahi0,
     &           idatalo1:idatahi1,
     &           0:ndatacomp-1)
      integer i0,i1, m, ncomp
      ncomp = ndatacomp
      do m = 0, ncomp-1
         
      do i1 = igridboxlo1,igridboxhi1
      do i0 = igridboxlo0,igridboxhi0

#if CH_SPACEDIM == 3
            data(i0,i1,m) = data(i0,i1,m)
     &           - radial(i0,iradiallo1,iradiallo2)
#else
            data(i0,i1,m) = data(i0,i1,m)
     &           - radial(i0,iradiallo1)
#endif
         
      enddo
      enddo
      enddo
      return
      end
      subroutine MULTIPLY_FLUX_SURFACE_ARRAY_2D(
     &           igridboxlo0,igridboxlo1
     &           ,igridboxhi0,igridboxhi1
     &           ,radial
     &           ,iradiallo0,iradiallo1
     &           ,iradialhi0,iradialhi1
     &           ,data
     &           ,idatalo0,idatalo1
     &           ,idatahi0,idatahi1
     &           ,ndatacomp
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer igridboxlo0,igridboxlo1
      integer igridboxhi0,igridboxhi1
      integer iradiallo0,iradiallo1
      integer iradialhi0,iradialhi1
      REAL_T radial(
     &           iradiallo0:iradialhi0,
     &           iradiallo1:iradialhi1)
      integer ndatacomp
      integer idatalo0,idatalo1
      integer idatahi0,idatahi1
      REAL_T data(
     &           idatalo0:idatahi0,
     &           idatalo1:idatahi1,
     &           0:ndatacomp-1)
      integer i0,i1, m, ncomp
      ncomp = ndatacomp
      do m = 0, ncomp-1
         
      do i1 = igridboxlo1,igridboxhi1
      do i0 = igridboxlo0,igridboxhi0

#if CH_SPACEDIM == 3
            data(i0,i1,m) = data(i0,i1,m)
     &           * radial(i0,iradiallo1,iradiallo2)
#else
            data(i0,i1,m) = data(i0,i1,m)
     &           * radial(i0,iradiallo1)
#endif
         
      enddo
      enddo
      enddo
      return
      end
      subroutine DIVIDE_FLUX_SURFACE_ARRAY_2D(
     &           igridboxlo0,igridboxlo1
     &           ,igridboxhi0,igridboxhi1
     &           ,radial
     &           ,iradiallo0,iradiallo1
     &           ,iradialhi0,iradialhi1
     &           ,data
     &           ,idatalo0,idatalo1
     &           ,idatahi0,idatahi1
     &           ,ndatacomp
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer igridboxlo0,igridboxlo1
      integer igridboxhi0,igridboxhi1
      integer iradiallo0,iradiallo1
      integer iradialhi0,iradialhi1
      REAL_T radial(
     &           iradiallo0:iradialhi0,
     &           iradiallo1:iradialhi1)
      integer ndatacomp
      integer idatalo0,idatalo1
      integer idatahi0,idatahi1
      REAL_T data(
     &           idatalo0:idatahi0,
     &           idatalo1:idatahi1,
     &           0:ndatacomp-1)
      integer i0,i1, m, ncomp
      ncomp = ndatacomp
      do m = 0, ncomp-1
         
      do i1 = igridboxlo1,igridboxhi1
      do i0 = igridboxlo0,igridboxhi0

#if CH_SPACEDIM == 3
            data(i0,i1,m) = data(i0,i1,m)
     &           / radial(i0,iradiallo1,iradiallo2)
#else
            data(i0,i1,m) = data(i0,i1,m)
     &           / radial(i0,iradiallo1)
#endif
         
      enddo
      enddo
      enddo
      return
      end
