/* Copyright (C) 1991-2015 Free Software Foundation, Inc.
   This file is part of the GNU C Library.
   The GNU C Library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.
   The GNU C Library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.
   You should have received a copy of the GNU Lesser General Public
   License along with the GNU C Library; if not, see
   <http://www.gnu.org/licenses/>.  */
/* This header is separate from features.h so that the compiler can
   include it implicitly at the start of every compilation.  It must
   not itself include <features.h> or any other header that includes
   <features.h> because the implicit include comes before any feature
   test macros that may be defined in a source file before it first
   explicitly includes a system header.  GCC knows the name of this
   header in order to preinclude it.  */
/* glibc's intent is to support the IEC 559 math functionality, real
   and complex.  If the GCC (4.9 and later) predefined macros
   specifying compiler intent are available, use them to determine
   whether the overall intent is to support these features; otherwise,
   presume an older compiler has intent to support these features and
   define these macros by default.  */
/* wchar_t uses Unicode 7.0.0.  Version 7.0 of the Unicode Standard is
   synchronized with ISO/IEC 10646:2012, plus Amendments 1 (published
   on April, 2013) and 2 (not yet published as of February, 2015).
   Additionally, it includes the accelerated publication of U+20BD
   RUBLE SIGN.  Therefore Unicode 7.0.0 is between 10646:2012 and
   10646:2014, and so we use the date ISO/IEC 10646:2012 Amd.1 was
   published.  */
/* We do not support C11 <threads.h>.  */
      subroutine GET_NC_MAPPED_COORDS_2D(
     & iboxlo0,iboxlo1
     & ,iboxhi0,iboxhi1
     & ,dx
     & ,xi
     & ,ixilo0,ixilo1
     & ,ixihi0,ixihi1
     & ,nxicomp
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      REAL*8 dx(0:1)
      integer nxicomp
      integer ixilo0,ixilo1
      integer ixihi0,ixihi1
      REAL*8 xi(
     & ixilo0:ixihi0,
     & ixilo1:ixihi1,
     & 0:nxicomp-1)
      integer i,j
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0
        xi(i,j,0) = i*dx(0)
                  xi(i,j,1) = j*dx(1)
      enddo
      enddo
      return
      end
      subroutine GET_CC_MAPPED_COORDS_2D(
     & iboxlo0,iboxlo1
     & ,iboxhi0,iboxhi1
     & ,dx
     & ,xi
     & ,ixilo0,ixilo1
     & ,ixihi0,ixihi1
     & ,nxicomp
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      REAL*8 dx(0:1)
      integer nxicomp
      integer ixilo0,ixilo1
      integer ixihi0,ixihi1
      REAL*8 xi(
     & ixilo0:ixihi0,
     & ixilo1:ixihi1,
     & 0:nxicomp-1)
      integer i,j
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0
        xi(i,j,0) = (i + (0.500d0))*dx(0)
                  xi(i,j,1) = (j + (0.500d0))*dx(1)
      enddo
      enddo
      return
      end
      subroutine GET_FC_MAPPED_COORDS_2D(
     & iboxlo0,iboxlo1
     & ,iboxhi0,iboxhi1
     & ,dir
     & ,dx
     & ,xi
     & ,ixilo0,ixilo1
     & ,ixihi0,ixihi1
     & ,nxicomp
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer dir
      REAL*8 dx(0:1)
      integer nxicomp
      integer ixilo0,ixilo1
      integer ixihi0,ixihi1
      REAL*8 xi(
     & ixilo0:ixihi0,
     & ixilo1:ixihi1,
     & 0:nxicomp-1)
      integer i,j
      double precision offset(0:2 -1)
      offset(0) = (0.500d0)
                offset(1) = (0.500d0)
      offset(dir) = (0.0d0)
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0
        xi(i,j,0) = (i + offset(0))*dx(0)
                  xi(i,j,1) = (j + offset(1))*dx(1)
      enddo
      enddo
      return
      end
      subroutine INCREMENTLAPLACIAN2_2D(
     & lapPhi
     & ,ilapPhilo0,ilapPhilo1
     & ,ilapPhihi0,ilapPhihi1
     & ,nlapPhicomp
     & ,phi
     & ,iphilo0,iphilo1
     & ,iphihi0,iphihi1
     & ,nphicomp
     & ,igridBoxlo0,igridBoxlo1
     & ,igridBoxhi0,igridBoxhi1
     & ,dir
     & ,factor
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer nlapPhicomp
      integer ilapPhilo0,ilapPhilo1
      integer ilapPhihi0,ilapPhihi1
      REAL*8 lapPhi(
     & ilapPhilo0:ilapPhihi0,
     & ilapPhilo1:ilapPhihi1,
     & 0:nlapPhicomp-1)
      integer nphicomp
      integer iphilo0,iphilo1
      integer iphihi0,iphihi1
      REAL*8 phi(
     & iphilo0:iphihi0,
     & iphilo1:iphihi1,
     & 0:nphicomp-1)
      integer igridBoxlo0,igridBoxlo1
      integer igridBoxhi0,igridBoxhi1
      integer dir
      REAL*8 factor
      integer i0,i1
      integer ii0,ii1
      integer comp
      REAL*8 thisLap
      ii0=CHF_ID(0,dir)
      ii1=CHF_ID(1,dir)
      do comp=0, nphicomp-1
      do i1 = igridBoxlo1,igridBoxhi1
      do i0 = igridBoxlo0,igridBoxhi0
            thisLap = phi(i0 +ii0,i1 +ii1, comp)
     & + phi(i0 -ii0,i1 -ii1, comp)
     & -(2.0d0)*phi(i0,i1, comp)
            lapPhi(i0,i1, comp) =
     & lapPhi(i0,i1, comp) + factor*thisLap
      enddo
      enddo
      enddo
      return
      end
      subroutine UNIT_FS_TANGENT_2D(
     & iboxlo0,iboxlo1
     & ,iboxhi0,iboxhi1
     & ,dXdxi
     & ,idXdxilo0,idXdxilo1
     & ,idXdxihi0,idXdxihi1
     & ,ndXdxicomp
     & ,data
     & ,idatalo0,idatalo1
     & ,idatahi0,idatahi1
     & ,ndatacomp
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer ndXdxicomp
      integer idXdxilo0,idXdxilo1
      integer idXdxihi0,idXdxihi1
      REAL*8 dXdxi(
     & idXdxilo0:idXdxihi0,
     & idXdxilo1:idXdxihi1,
     & 0:ndXdxicomp-1)
      integer ndatacomp
      integer idatalo0,idatalo1
      integer idatahi0,idatahi1
      REAL*8 data(
     & idatalo0:idatahi0,
     & idatalo1:idatahi1,
     & 0:ndatacomp-1)
      integer i,j
      double precision fac
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0
        fac = (1.0d0) / dsqrt( dXdxi(i,j,1)**2 + dXdxi(i,j,3)**2 )
        data(i,j,0) = dXdxi(i,j,1) * fac
        data(i,j,1) = dXdxi(i,j,3) * fac
      enddo
      enddo
      return
      end
      subroutine UNIT_FS_NORMAL_2D(
     & iboxlo0,iboxlo1
     & ,iboxhi0,iboxhi1
     & ,dXdxi
     & ,idXdxilo0,idXdxilo1
     & ,idXdxihi0,idXdxihi1
     & ,ndXdxicomp
     & ,data
     & ,idatalo0,idatalo1
     & ,idatahi0,idatahi1
     & ,ndatacomp
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer ndXdxicomp
      integer idXdxilo0,idXdxilo1
      integer idXdxihi0,idXdxihi1
      REAL*8 dXdxi(
     & idXdxilo0:idXdxihi0,
     & idXdxilo1:idXdxihi1,
     & 0:ndXdxicomp-1)
      integer ndatacomp
      integer idatalo0,idatalo1
      integer idatahi0,idatahi1
      REAL*8 data(
     & idatalo0:idatahi0,
     & idatalo1:idatahi1,
     & 0:ndatacomp-1)
      integer i,j
      double precision fac
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0
        fac = (1.0d0) / dsqrt( dXdxi(i,j,1)**2 + dXdxi(i,j,3)**2 )
        data(i,j,0) = dXdxi(i,j,3) * fac
        data(i,j,1) = -dXdxi(i,j,1) * fac
      enddo
      enddo
      return
      end
      subroutine UNIT_RADIAL_2D(
     & iboxlo0,iboxlo1
     & ,iboxhi0,iboxhi1
     & ,dXdxi
     & ,idXdxilo0,idXdxilo1
     & ,idXdxihi0,idXdxihi1
     & ,ndXdxicomp
     & ,data
     & ,idatalo0,idatalo1
     & ,idatahi0,idatahi1
     & ,ndatacomp
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer ndXdxicomp
      integer idXdxilo0,idXdxilo1
      integer idXdxihi0,idXdxihi1
      REAL*8 dXdxi(
     & idXdxilo0:idXdxihi0,
     & idXdxilo1:idXdxihi1,
     & 0:ndXdxicomp-1)
      integer ndatacomp
      integer idatalo0,idatalo1
      integer idatahi0,idatahi1
      REAL*8 data(
     & idatalo0:idatahi0,
     & idatalo1:idatahi1,
     & 0:ndatacomp-1)
      integer i,j
      double precision fac
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0
        fac = (1.0d0) / dsqrt( dXdxi(i,j,0)**2 + dXdxi(i,j,2)**2 )
        data(i,j,0) = dXdxi(i,j,0) * fac
        data(i,j,1) = dXdxi(i,j,2) * fac
      enddo
      enddo
      return
      end
      subroutine GRADF_FACTOR_2D(
     & iboxlo0,iboxlo1
     & ,iboxhi0,iboxhi1
     & ,dXdxi
     & ,idXdxilo0,idXdxilo1
     & ,idXdxihi0,idXdxihi1
     & ,ndXdxicomp
     & ,data
     & ,idatalo0,idatalo1
     & ,idatahi0,idatahi1
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer ndXdxicomp
      integer idXdxilo0,idXdxilo1
      integer idXdxihi0,idXdxihi1
      REAL*8 dXdxi(
     & idXdxilo0:idXdxihi0,
     & idXdxilo1:idXdxihi1,
     & 0:ndXdxicomp-1)
      integer idatalo0,idatalo1
      integer idatahi0,idatahi1
      REAL*8 data(
     & idatalo0:idatahi0,
     & idatalo1:idatahi1)
      integer i,j
      double precision fac1, fac2, v1(0:1), v2(0:1)
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0
        fac1 = (1.0d0) / dsqrt( dXdxi(i,j,0)**2 + dXdxi(i,j,2)**2 )
        fac2 = (1.0d0) / dsqrt( dXdxi(i,j,1)**2 + dXdxi(i,j,3)**2 )
        v1(0) = dXdxi(i,j,0) * fac1
        v1(1) = dXdxi(i,j,2) * fac1
        v2(0) = dXdxi(i,j,3) * fac2
        v2(1) = -dXdxi(i,j,1) * fac2
        data(i,j) = fac1 * ( v1(0)*v2(0) + v1(1)*v2(1) )
      enddo
      enddo
      return
      end
      subroutine MAG_BLOCK_PROJECT_2D(
     & iboxlo0,iboxlo1
     & ,iboxhi0,iboxhi1
     & ,vec_src
     & ,ivec_srclo0,ivec_srclo1
     & ,ivec_srchi0,ivec_srchi1
     & ,nvec_srccomp
     & ,vec_dst
     & ,ivec_dstlo0,ivec_dstlo1
     & ,ivec_dsthi0,ivec_dsthi1
     & ,nvec_dstcomp
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer nvec_srccomp
      integer ivec_srclo0,ivec_srclo1
      integer ivec_srchi0,ivec_srchi1
      REAL*8 vec_src(
     & ivec_srclo0:ivec_srchi0,
     & ivec_srclo1:ivec_srchi1,
     & 0:nvec_srccomp-1)
      integer nvec_dstcomp
      integer ivec_dstlo0,ivec_dstlo1
      integer ivec_dsthi0,ivec_dsthi1
      REAL*8 vec_dst(
     & ivec_dstlo0:ivec_dsthi0,
     & ivec_dstlo1:ivec_dsthi1,
     & 0:nvec_dstcomp-1)
      integer i,j, comp, ncomp
      double precision dotprod
      ncomp = nvec_srccomp
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0
        dotprod = (0.0d0)
        do comp = 0, ncomp - 1
          dotprod = dotprod + vec_src(i,j,comp) * vec_dst(i,j,comp)
        enddo
        do comp = 0, ncomp - 1
          vec_dst(i,j,comp) = dotprod * vec_src(i,j,comp)
        enddo
      enddo
      enddo
      return
      end
      subroutine MAG_BLOCK_PSITHETA_PROJECTIONS_2D(
     & iboxlo0,iboxlo1
     & ,iboxhi0,iboxhi1
     & ,vec_psi
     & ,ivec_psilo0,ivec_psilo1
     & ,ivec_psihi0,ivec_psihi1
     & ,nvec_psicomp
     & ,vec_theta
     & ,ivec_thetalo0,ivec_thetalo1
     & ,ivec_thetahi0,ivec_thetahi1
     & ,nvec_thetacomp
     & ,vec_dst
     & ,ivec_dstlo0,ivec_dstlo1
     & ,ivec_dsthi0,ivec_dsthi1
     & ,nvec_dstcomp
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer nvec_psicomp
      integer ivec_psilo0,ivec_psilo1
      integer ivec_psihi0,ivec_psihi1
      REAL*8 vec_psi(
     & ivec_psilo0:ivec_psihi0,
     & ivec_psilo1:ivec_psihi1,
     & 0:nvec_psicomp-1)
      integer nvec_thetacomp
      integer ivec_thetalo0,ivec_thetalo1
      integer ivec_thetahi0,ivec_thetahi1
      REAL*8 vec_theta(
     & ivec_thetalo0:ivec_thetahi0,
     & ivec_thetalo1:ivec_thetahi1,
     & 0:nvec_thetacomp-1)
      integer nvec_dstcomp
      integer ivec_dstlo0,ivec_dstlo1
      integer ivec_dsthi0,ivec_dsthi1
      REAL*8 vec_dst(
     & ivec_dstlo0:ivec_dsthi0,
     & ivec_dstlo1:ivec_dsthi1,
     & 0:nvec_dstcomp-1)
      integer i,j, comp, ncomp
      double precision dotprod_psi, dotprod_theta
      ncomp = nvec_dstcomp
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0
        dotprod_psi = (0.0d0)
        do comp = 0, ncomp - 1
          dotprod_psi = dotprod_psi + vec_psi(i,j,comp) * vec_dst(i,j,co
     &mp)
        enddo
        dotprod_theta = (0.0d0)
        do comp = 0, ncomp - 1
          dotprod_theta = dotprod_theta + vec_theta(i,j,comp) * vec_dst(
     &i,j,comp)
        enddo
        vec_dst(i,j,0) = dotprod_psi
        vec_dst(i,j,1) = dotprod_theta
      enddo
      enddo
      return
      end
      subroutine GET_NODAL_FIELD_DATA_2D(
     & igridboxlo0,igridboxlo1
     & ,igridboxhi0,igridboxhi1
     & ,RZ
     & ,iRZlo0,iRZlo1
     & ,iRZhi0,iRZhi1
     & ,nRZcomp
     & ,psi
     & ,ipsilo0,ipsilo1
     & ,ipsihi0,ipsihi1
     & ,RB
     & ,iRBlo0,iRBlo1
     & ,iRBhi0,iRBhi1
     & ,nRBcomp
     & ,RBtor
     & ,A
     & ,iAlo0,iAlo1
     & ,iAhi0,iAhi1
     & ,nAcomp
     & ,bunit
     & ,ibunitlo0,ibunitlo1
     & ,ibunithi0,ibunithi1
     & ,nbunitcomp
     & ,Bmag
     & ,iBmaglo0,iBmaglo1
     & ,iBmaghi0,iBmaghi1
     & ,nBmagcomp
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer igridboxlo0,igridboxlo1
      integer igridboxhi0,igridboxhi1
      integer nRZcomp
      integer iRZlo0,iRZlo1
      integer iRZhi0,iRZhi1
      REAL*8 RZ(
     & iRZlo0:iRZhi0,
     & iRZlo1:iRZhi1,
     & 0:nRZcomp-1)
      integer ipsilo0,ipsilo1
      integer ipsihi0,ipsihi1
      REAL*8 psi(
     & ipsilo0:ipsihi0,
     & ipsilo1:ipsihi1)
      integer nRBcomp
      integer iRBlo0,iRBlo1
      integer iRBhi0,iRBhi1
      REAL*8 RB(
     & iRBlo0:iRBhi0,
     & iRBlo1:iRBhi1,
     & 0:nRBcomp-1)
      REAL*8 RBtor
      integer nAcomp
      integer iAlo0,iAlo1
      integer iAhi0,iAhi1
      REAL*8 A(
     & iAlo0:iAhi0,
     & iAlo1:iAhi1,
     & 0:nAcomp-1)
      integer nbunitcomp
      integer ibunitlo0,ibunitlo1
      integer ibunithi0,ibunithi1
      REAL*8 bunit(
     & ibunitlo0:ibunithi0,
     & ibunitlo1:ibunithi1,
     & 0:nbunitcomp-1)
      integer nBmagcomp
      integer iBmaglo0,iBmaglo1
      integer iBmaghi0,iBmaghi1
      REAL*8 Bmag(
     & iBmaglo0:iBmaghi0,
     & iBmaglo1:iBmaghi1,
     & 0:nBmagcomp-1)
      integer i,j, n
      double precision R, Z, Bmagnitude, B(0:2)
      do j = igridboxlo1,igridboxhi1
      do i = igridboxlo0,igridboxhi0
         R = RZ(i,j,0);
         Z = RZ(i,j,1);
         B(0) = RB(i,j,0) / R
         B(1) = RBtor / R
         B(2) = RB(i,j,1) / R
         Bmagnitude = (0.0d0)
         do n = 0, 2
            Bmagnitude = Bmagnitude + B(n)**2
         enddo
         Bmagnitude = sqrt(Bmagnitude)
         do n = 0, 2
            bunit(i,j,n) = B(n) / Bmagnitude
         enddo
         Bmag(i,j,0) = Bmagnitude
         A(i,j,0) = Z * RBtor / R
         A(i,j,1) = psi(i,j) / R
         A(i,j,2) = (0.0d0)
      enddo
      enddo
      return
      end
      subroutine GET_FIELD_DATA_2D(
     & igridboxlo0,igridboxlo1
     & ,igridboxhi0,igridboxhi1
     & ,RZ
     & ,iRZlo0,iRZlo1
     & ,iRZhi0,iRZhi1
     & ,nRZcomp
     & ,RB
     & ,iRBlo0,iRBlo1
     & ,iRBhi0,iRBhi1
     & ,nRBcomp
     & ,dRBdR
     & ,idRBdRlo0,idRBdRlo1
     & ,idRBdRhi0,idRBdRhi1
     & ,ndRBdRcomp
     & ,dRBdZ
     & ,idRBdZlo0,idRBdZlo1
     & ,idRBdZhi0,idRBdZhi1
     & ,ndRBdZcomp
     & ,RBtor
     & ,B
     & ,iBlo0,iBlo1
     & ,iBhi0,iBhi1
     & ,nBcomp
     & ,Bmag
     & ,iBmaglo0,iBmaglo1
     & ,iBmaghi0,iBmaghi1
     & ,bunit
     & ,ibunitlo0,ibunitlo1
     & ,ibunithi0,ibunithi1
     & ,nbunitcomp
     & ,gradB
     & ,igradBlo0,igradBlo1
     & ,igradBhi0,igradBhi1
     & ,ngradBcomp
     & ,curlb
     & ,icurlblo0,icurlblo1
     & ,icurlbhi0,icurlbhi1
     & ,ncurlbcomp
     & ,bdotcurlb
     & ,ibdotcurlblo0,ibdotcurlblo1
     & ,ibdotcurlbhi0,ibdotcurlbhi1
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer igridboxlo0,igridboxlo1
      integer igridboxhi0,igridboxhi1
      integer nRZcomp
      integer iRZlo0,iRZlo1
      integer iRZhi0,iRZhi1
      REAL*8 RZ(
     & iRZlo0:iRZhi0,
     & iRZlo1:iRZhi1,
     & 0:nRZcomp-1)
      integer nRBcomp
      integer iRBlo0,iRBlo1
      integer iRBhi0,iRBhi1
      REAL*8 RB(
     & iRBlo0:iRBhi0,
     & iRBlo1:iRBhi1,
     & 0:nRBcomp-1)
      integer ndRBdRcomp
      integer idRBdRlo0,idRBdRlo1
      integer idRBdRhi0,idRBdRhi1
      REAL*8 dRBdR(
     & idRBdRlo0:idRBdRhi0,
     & idRBdRlo1:idRBdRhi1,
     & 0:ndRBdRcomp-1)
      integer ndRBdZcomp
      integer idRBdZlo0,idRBdZlo1
      integer idRBdZhi0,idRBdZhi1
      REAL*8 dRBdZ(
     & idRBdZlo0:idRBdZhi0,
     & idRBdZlo1:idRBdZhi1,
     & 0:ndRBdZcomp-1)
      REAL*8 RBtor
      integer nBcomp
      integer iBlo0,iBlo1
      integer iBhi0,iBhi1
      REAL*8 B(
     & iBlo0:iBhi0,
     & iBlo1:iBhi1,
     & 0:nBcomp-1)
      integer iBmaglo0,iBmaglo1
      integer iBmaghi0,iBmaghi1
      REAL*8 Bmag(
     & iBmaglo0:iBmaghi0,
     & iBmaglo1:iBmaghi1)
      integer nbunitcomp
      integer ibunitlo0,ibunitlo1
      integer ibunithi0,ibunithi1
      REAL*8 bunit(
     & ibunitlo0:ibunithi0,
     & ibunitlo1:ibunithi1,
     & 0:nbunitcomp-1)
      integer ngradBcomp
      integer igradBlo0,igradBlo1
      integer igradBhi0,igradBhi1
      REAL*8 gradB(
     & igradBlo0:igradBhi0,
     & igradBlo1:igradBhi1,
     & 0:ngradBcomp-1)
      integer ncurlbcomp
      integer icurlblo0,icurlblo1
      integer icurlbhi0,icurlbhi1
      REAL*8 curlb(
     & icurlblo0:icurlbhi0,
     & icurlblo1:icurlbhi1,
     & 0:ncurlbcomp-1)
      integer ibdotcurlblo0,ibdotcurlblo1
      integer ibdotcurlbhi0,ibdotcurlbhi1
      REAL*8 bdotcurlb(
     & ibdotcurlblo0:ibdotcurlbhi0,
     & ibdotcurlblo1:ibdotcurlbhi1)
      integer i,j, n
      double precision R, Z, Bmag2
      do j = igridboxlo1,igridboxhi1
      do i = igridboxlo0,igridboxhi0
         R = RZ(i,j,0);
         Z = RZ(i,j,1);
         B(i,j,0) = RB(i,j,0) / R
         B(i,j,1) = RBtor / R
         B(i,j,2) = RB(i,j,1) / R
         Bmag2 = (0.0d0)
         do n = 0, 2
            Bmag2 = Bmag2 + B(i,j,n)**2
         enddo
         Bmag(i,j) = sqrt(Bmag2)
         do n = 0, 2
            bunit(i,j,n) = B(i,j,n) / Bmag(i,j)
         enddo
         gradB(i,j,0) = (B(i,j,0) * dRBdR(i,j,0) + B(i,j,2) * dRBdR(i,j,
     &1))
     & / (R * Bmag(i,j)) - Bmag(i,j) / R
         gradB(i,j,1) = (0.0d0)
         gradB(i,j,2) = (B(i,j,0) * dRBdZ(i,j,0) + B(i,j,2) * dRBdZ(i,j,
     &1))
     & / (R * Bmag(i,j))
         curlb(i,j,0) = bunit(i,j,1) * gradB(i,j,2) / Bmag(i,j)
         curlb(i,j,1) = (dRBdZ(i,j,0) - dRBdR(i,j,1)) / (R * Bmag(i,j))
     & + (bunit(i,j,2) * gradB(i,j,0)
     & - bunit(i,j,0) * gradB(i,j,2)) / Bmag(i,j)
     & + bunit(i,j,2) / R
         curlb(i,j,2) = - bunit(i,j,1) * gradB(i,j,0) / Bmag(i,j)
         bdotcurlb(i,j) = (0.0d0)
         do n = 0, 2
            bdotcurlb(i,j) = bdotcurlb(i,j) + bunit(i,j,n) * curlb(i,j,n
     &)
         enddo
      enddo
      enddo
      return
      end
