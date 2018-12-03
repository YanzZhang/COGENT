#include "REAL.H"
#include "SPACE.H"
#include "CONSTANTS.H"

      subroutine PHASE_BLOCK_PROJECT_4D(
     &           iboxlo0,iboxlo1,iboxlo2,iboxlo3
     &           ,iboxhi0,iboxhi1,iboxhi2,iboxhi3
     &           ,vec_src
     &           ,ivec_srclo0,ivec_srclo1,ivec_srclo2,ivec_srclo3
     &           ,ivec_srchi0,ivec_srchi1,ivec_srchi2,ivec_srchi3
     &           ,nvec_srccomp
     &           ,vec_dst
     &           ,ivec_dstlo0,ivec_dstlo1,ivec_dstlo2,ivec_dstlo3
     &           ,ivec_dsthi0,ivec_dsthi1,ivec_dsthi2,ivec_dsthi3
     &           ,nvec_dstcomp
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iboxlo0,iboxlo1,iboxlo2,iboxlo3
      integer iboxhi0,iboxhi1,iboxhi2,iboxhi3
      integer nvec_srccomp
      integer ivec_srclo0,ivec_srclo1,ivec_srclo2,ivec_srclo3
      integer ivec_srchi0,ivec_srchi1,ivec_srchi2,ivec_srchi3
      REAL_T vec_src(
     &           ivec_srclo0:ivec_srchi0,
     &           ivec_srclo1:ivec_srchi1,
     &           ivec_srclo2:ivec_srchi2,
     &           ivec_srclo3:ivec_srchi3,
     &           0:nvec_srccomp-1)
      integer nvec_dstcomp
      integer ivec_dstlo0,ivec_dstlo1,ivec_dstlo2,ivec_dstlo3
      integer ivec_dsthi0,ivec_dsthi1,ivec_dsthi2,ivec_dsthi3
      REAL_T vec_dst(
     &           ivec_dstlo0:ivec_dsthi0,
     &           ivec_dstlo1:ivec_dsthi1,
     &           ivec_dstlo2:ivec_dsthi2,
     &           ivec_dstlo3:ivec_dsthi3,
     &           0:nvec_dstcomp-1)
      integer i0,i1,i2,i3, comp, ncomp, kk
      double precision dotprod
      ncomp = nvec_srccomp
      kk = ivec_srclo2
      
      do i3 = iboxlo3,iboxhi3
      do i2 = iboxlo2,iboxhi2
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0

#if CH_SPACEDIM==5
        kk = i2
#endif
        dotprod = zero
        do comp = 0, ncomp - 1
          dotprod = dotprod + vec_src(i0,i1,kk,ivec_srclo3,comp)
     &                      * vec_dst(i0,i1,i2,i3,comp)
        enddo     
        do comp = 0, ncomp - 1
          vec_dst(i0,i1,i2,i3,comp) = dotprod * vec_src(i0,i1,kk,ivec_srclo3,comp)
        enddo     
      
      enddo
      enddo
      enddo
      enddo
      return
      end
      subroutine PHASE_BLOCK_PSITHETA_PROJECTIONS_4D(
     &           iboxlo0,iboxlo1,iboxlo2,iboxlo3
     &           ,iboxhi0,iboxhi1,iboxhi2,iboxhi3
     &           ,vec_psi
     &           ,ivec_psilo0,ivec_psilo1,ivec_psilo2,ivec_psilo3
     &           ,ivec_psihi0,ivec_psihi1,ivec_psihi2,ivec_psihi3
     &           ,nvec_psicomp
     &           ,vec_theta
     &           ,ivec_thetalo0,ivec_thetalo1,ivec_thetalo2,ivec_thetalo3
     &           ,ivec_thetahi0,ivec_thetahi1,ivec_thetahi2,ivec_thetahi3
     &           ,nvec_thetacomp
     &           ,vec_dst
     &           ,ivec_dstlo0,ivec_dstlo1,ivec_dstlo2,ivec_dstlo3
     &           ,ivec_dsthi0,ivec_dsthi1,ivec_dsthi2,ivec_dsthi3
     &           ,nvec_dstcomp
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iboxlo0,iboxlo1,iboxlo2,iboxlo3
      integer iboxhi0,iboxhi1,iboxhi2,iboxhi3
      integer nvec_psicomp
      integer ivec_psilo0,ivec_psilo1,ivec_psilo2,ivec_psilo3
      integer ivec_psihi0,ivec_psihi1,ivec_psihi2,ivec_psihi3
      REAL_T vec_psi(
     &           ivec_psilo0:ivec_psihi0,
     &           ivec_psilo1:ivec_psihi1,
     &           ivec_psilo2:ivec_psihi2,
     &           ivec_psilo3:ivec_psihi3,
     &           0:nvec_psicomp-1)
      integer nvec_thetacomp
      integer ivec_thetalo0,ivec_thetalo1,ivec_thetalo2,ivec_thetalo3
      integer ivec_thetahi0,ivec_thetahi1,ivec_thetahi2,ivec_thetahi3
      REAL_T vec_theta(
     &           ivec_thetalo0:ivec_thetahi0,
     &           ivec_thetalo1:ivec_thetahi1,
     &           ivec_thetalo2:ivec_thetahi2,
     &           ivec_thetalo3:ivec_thetahi3,
     &           0:nvec_thetacomp-1)
      integer nvec_dstcomp
      integer ivec_dstlo0,ivec_dstlo1,ivec_dstlo2,ivec_dstlo3
      integer ivec_dsthi0,ivec_dsthi1,ivec_dsthi2,ivec_dsthi3
      REAL_T vec_dst(
     &           ivec_dstlo0:ivec_dsthi0,
     &           ivec_dstlo1:ivec_dsthi1,
     &           ivec_dstlo2:ivec_dsthi2,
     &           ivec_dstlo3:ivec_dsthi3,
     &           0:nvec_dstcomp-1)
      integer i0,i1,i2,i3, comp, kk
      double precision dotprod_psi, dotprod_theta 
      kk = ivec_psilo2
      
      do i3 = iboxlo3,iboxhi3
      do i2 = iboxlo2,iboxhi2
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0

#if CH_SPACEDIM==5
        kk = i2
#endif
        dotprod_psi = zero
        do comp = 0, 1
          dotprod_psi = dotprod_psi + vec_psi(i0,i1,kk,ivec_psilo3,comp)
     &                         * vec_dst(i0,i1,i2,i3,comp)
        enddo     
        dotprod_theta = zero
        do comp = 0, 1
          dotprod_theta = dotprod_theta + vec_theta(i0,i1,kk,ivec_thetalo3,comp)
     &                         * vec_dst(i0,i1,i2,i3,comp)
        enddo     
        vec_dst(i0,i1,i2,i3,0) = dotprod_psi
        vec_dst(i0,i1,i2,i3,1) = dotprod_theta
      
      enddo
      enddo
      enddo
      enddo
      return
      end
      subroutine PHASE_BLOCK_GRADF_FACTOR_4D(
     &           iboxlo0,iboxlo1,iboxlo2,iboxlo3
     &           ,iboxhi0,iboxhi1,iboxhi2,iboxhi3
     &           ,vec_src
     &           ,ivec_srclo0,ivec_srclo1,ivec_srclo2,ivec_srclo3
     &           ,ivec_srchi0,ivec_srchi1,ivec_srchi2,ivec_srchi3
     &           ,vec_dst
     &           ,ivec_dstlo0,ivec_dstlo1,ivec_dstlo2,ivec_dstlo3
     &           ,ivec_dsthi0,ivec_dsthi1,ivec_dsthi2,ivec_dsthi3
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iboxlo0,iboxlo1,iboxlo2,iboxlo3
      integer iboxhi0,iboxhi1,iboxhi2,iboxhi3
      integer ivec_srclo0,ivec_srclo1,ivec_srclo2,ivec_srclo3
      integer ivec_srchi0,ivec_srchi1,ivec_srchi2,ivec_srchi3
      REAL_T vec_src(
     &           ivec_srclo0:ivec_srchi0,
     &           ivec_srclo1:ivec_srchi1,
     &           ivec_srclo2:ivec_srchi2,
     &           ivec_srclo3:ivec_srchi3)
      integer ivec_dstlo0,ivec_dstlo1,ivec_dstlo2,ivec_dstlo3
      integer ivec_dsthi0,ivec_dsthi1,ivec_dsthi2,ivec_dsthi3
      REAL_T vec_dst(
     &           ivec_dstlo0:ivec_dsthi0,
     &           ivec_dstlo1:ivec_dsthi1,
     &           ivec_dstlo2:ivec_dsthi2,
     &           ivec_dstlo3:ivec_dsthi3)
      integer i0,i1,i2,i3, kk
      
      do i3 = iboxlo3,iboxhi3
      do i2 = iboxlo2,iboxhi2
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0

#if CH_SPACEDIM==5
        kk = i2
#endif
        vec_dst(i0,i1,i2,i3) = vec_src(i0,i1,kk,ivec_srclo3)
      
      enddo
      enddo
      enddo
      enddo
      return
      end
