#include "REAL.H"
#include "SPACE.H"
#include "CONSTANTS.H"

      subroutine COMPUTE_CONDUCTIVITY_2D(
     &           iboxlo0,iboxlo1
     &           ,iboxhi0,iboxhi1
     &           ,sigma
     &           ,isigmalo0,isigmalo1
     &           ,isigmahi0,isigmahi1
     &           ,T
     &           ,iTlo0,iTlo1
     &           ,iThi0,iThi1
     &           ,coeff
     &           ,sigma_max
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer isigmalo0,isigmalo1
      integer isigmahi0,isigmahi1
      REAL_T sigma(
     &           isigmalo0:isigmahi0,
     &           isigmalo1:isigmahi1)
      integer iTlo0,iTlo1
      integer iThi0,iThi1
      REAL_T T(
     &           iTlo0:iThi0,
     &           iTlo1:iThi1)
      REAL_T coeff
      REAL_T sigma_max
      integer i,j
      
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0

       sigma(i,j) = coeff * T(i,j)**(three/two)
       if (sigma(i,j).ge.sigma_max) then
          sigma(i,j) = sigma_max
       endif
      
      enddo
      enddo
      return
      end
