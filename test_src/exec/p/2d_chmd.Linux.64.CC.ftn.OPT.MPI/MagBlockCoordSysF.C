#include "REAL.H"
#include "SPACE.H"
#include "CONSTANTS.H"

      subroutine GET_NC_MAPPED_COORDS_2D(
     &           iboxlo0,iboxlo1
     &           ,iboxhi0,iboxhi1
     &           ,dx
     &           ,xi
     &           ,ixilo0,ixilo1
     &           ,ixihi0,ixihi1
     &           ,nxicomp
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      REAL_T dx(0:1)
      integer nxicomp
      integer ixilo0,ixilo1
      integer ixihi0,ixihi1
      REAL_T xi(
     &           ixilo0:ixihi0,
     &           ixilo1:ixihi1,
     &           0:nxicomp-1)
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
     &           iboxlo0,iboxlo1
     &           ,iboxhi0,iboxhi1
     &           ,dx
     &           ,xi
     &           ,ixilo0,ixilo1
     &           ,ixihi0,ixihi1
     &           ,nxicomp
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      REAL_T dx(0:1)
      integer nxicomp
      integer ixilo0,ixilo1
      integer ixihi0,ixihi1
      REAL_T xi(
     &           ixilo0:ixihi0,
     &           ixilo1:ixihi1,
     &           0:nxicomp-1)
      integer i,j
      
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0

        xi(i,j,0) = (i + half)*dx(0)
                  xi(i,j,1) = (j + half)*dx(1)
      
      enddo
      enddo
      return
      end
      subroutine GET_FC_MAPPED_COORDS_2D(
     &           iboxlo0,iboxlo1
     &           ,iboxhi0,iboxhi1
     &           ,dir
     &           ,dx
     &           ,xi
     &           ,ixilo0,ixilo1
     &           ,ixihi0,ixihi1
     &           ,nxicomp
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer dir
      REAL_T dx(0:1)
      integer nxicomp
      integer ixilo0,ixilo1
      integer ixihi0,ixihi1
      REAL_T xi(
     &           ixilo0:ixihi0,
     &           ixilo1:ixihi1,
     &           0:nxicomp-1)
      integer i,j
      double precision offset(0:CH_SPACEDIM-1)
      offset(0) = half
                offset(1) = half
      offset(dir) = zero
      
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0

        xi(i,j,0) = (i + offset(0))*dx(0)
                  xi(i,j,1) = (j + offset(1))*dx(1)
      
      enddo
      enddo
      return
      end
      subroutine INCREMENTLAPLACIAN2_2D(
     &           lapPhi
     &           ,ilapPhilo0,ilapPhilo1
     &           ,ilapPhihi0,ilapPhihi1
     &           ,nlapPhicomp
     &           ,phi
     &           ,iphilo0,iphilo1
     &           ,iphihi0,iphihi1
     &           ,nphicomp
     &           ,igridBoxlo0,igridBoxlo1
     &           ,igridBoxhi0,igridBoxhi1
     &           ,dir
     &           ,factor
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer nlapPhicomp
      integer ilapPhilo0,ilapPhilo1
      integer ilapPhihi0,ilapPhihi1
      REAL_T lapPhi(
     &           ilapPhilo0:ilapPhihi0,
     &           ilapPhilo1:ilapPhihi1,
     &           0:nlapPhicomp-1)
      integer nphicomp
      integer iphilo0,iphilo1
      integer iphihi0,iphihi1
      REAL_T phi(
     &           iphilo0:iphihi0,
     &           iphilo1:iphihi1,
     &           0:nphicomp-1)
      integer igridBoxlo0,igridBoxlo1
      integer igridBoxhi0,igridBoxhi1
      integer dir
      REAL_T factor
      integer i0,i1
      integer ii0,ii1
      integer comp
      REAL_T thisLap
      
      ii0=CHF_ID(0,dir)

      ii1=CHF_ID(1,dir)

      do comp=0, nphicomp-1
         
      do i1 = igridBoxlo1,igridBoxhi1
      do i0 = igridBoxlo0,igridBoxhi0

            thisLap = phi(i0 +ii0,i1 +ii1, comp)
     &              + phi(i0 -ii0,i1 -ii1, comp)
     &           -two*phi(i0,i1, comp)
            lapPhi(i0,i1, comp) =
     &           lapPhi(i0,i1, comp) + factor*thisLap
         
      enddo
      enddo
      enddo
      return
      end
      subroutine UNIT_FS_TANGENT_2D(
     &           iboxlo0,iboxlo1
     &           ,iboxhi0,iboxhi1
     &           ,dXdxi
     &           ,idXdxilo0,idXdxilo1
     &           ,idXdxihi0,idXdxihi1
     &           ,ndXdxicomp
     &           ,data
     &           ,idatalo0,idatalo1
     &           ,idatahi0,idatahi1
     &           ,ndatacomp
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer ndXdxicomp
      integer idXdxilo0,idXdxilo1
      integer idXdxihi0,idXdxihi1
      REAL_T dXdxi(
     &           idXdxilo0:idXdxihi0,
     &           idXdxilo1:idXdxihi1,
     &           0:ndXdxicomp-1)
      integer ndatacomp
      integer idatalo0,idatalo1
      integer idatahi0,idatahi1
      REAL_T data(
     &           idatalo0:idatahi0,
     &           idatalo1:idatahi1,
     &           0:ndatacomp-1)
      integer i,j
      double precision fac
      
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0

        fac = one / dsqrt( dXdxi(i,j,1)**2 + dXdxi(i,j,3)**2 )
        data(i,j,0) = dXdxi(i,j,1) * fac
        data(i,j,1) = dXdxi(i,j,3) * fac
      
      enddo
      enddo
      return
      end
      subroutine UNIT_FS_NORMAL_2D(
     &           iboxlo0,iboxlo1
     &           ,iboxhi0,iboxhi1
     &           ,dXdxi
     &           ,idXdxilo0,idXdxilo1
     &           ,idXdxihi0,idXdxihi1
     &           ,ndXdxicomp
     &           ,data
     &           ,idatalo0,idatalo1
     &           ,idatahi0,idatahi1
     &           ,ndatacomp
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer ndXdxicomp
      integer idXdxilo0,idXdxilo1
      integer idXdxihi0,idXdxihi1
      REAL_T dXdxi(
     &           idXdxilo0:idXdxihi0,
     &           idXdxilo1:idXdxihi1,
     &           0:ndXdxicomp-1)
      integer ndatacomp
      integer idatalo0,idatalo1
      integer idatahi0,idatahi1
      REAL_T data(
     &           idatalo0:idatahi0,
     &           idatalo1:idatahi1,
     &           0:ndatacomp-1)
      integer i,j
      double precision fac
      
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0

        fac = one / dsqrt( dXdxi(i,j,1)**2 + dXdxi(i,j,3)**2 )
        data(i,j,0) =  dXdxi(i,j,3) * fac
        data(i,j,1) = -dXdxi(i,j,1) * fac
      
      enddo
      enddo
      return
      end
      subroutine UNIT_RADIAL_2D(
     &           iboxlo0,iboxlo1
     &           ,iboxhi0,iboxhi1
     &           ,dXdxi
     &           ,idXdxilo0,idXdxilo1
     &           ,idXdxihi0,idXdxihi1
     &           ,ndXdxicomp
     &           ,data
     &           ,idatalo0,idatalo1
     &           ,idatahi0,idatahi1
     &           ,ndatacomp
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer ndXdxicomp
      integer idXdxilo0,idXdxilo1
      integer idXdxihi0,idXdxihi1
      REAL_T dXdxi(
     &           idXdxilo0:idXdxihi0,
     &           idXdxilo1:idXdxihi1,
     &           0:ndXdxicomp-1)
      integer ndatacomp
      integer idatalo0,idatalo1
      integer idatahi0,idatahi1
      REAL_T data(
     &           idatalo0:idatahi0,
     &           idatalo1:idatahi1,
     &           0:ndatacomp-1)
      integer i,j
      double precision fac
      
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0

        fac = one / dsqrt( dXdxi(i,j,0)**2 + dXdxi(i,j,2)**2 )
        data(i,j,0) = dXdxi(i,j,0) * fac
        data(i,j,1) = dXdxi(i,j,2) * fac
      
      enddo
      enddo
      return
      end
      subroutine GRADF_FACTOR_2D(
     &           iboxlo0,iboxlo1
     &           ,iboxhi0,iboxhi1
     &           ,dXdxi
     &           ,idXdxilo0,idXdxilo1
     &           ,idXdxihi0,idXdxihi1
     &           ,ndXdxicomp
     &           ,data
     &           ,idatalo0,idatalo1
     &           ,idatahi0,idatahi1
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer ndXdxicomp
      integer idXdxilo0,idXdxilo1
      integer idXdxihi0,idXdxihi1
      REAL_T dXdxi(
     &           idXdxilo0:idXdxihi0,
     &           idXdxilo1:idXdxihi1,
     &           0:ndXdxicomp-1)
      integer idatalo0,idatalo1
      integer idatahi0,idatahi1
      REAL_T data(
     &           idatalo0:idatahi0,
     &           idatalo1:idatahi1)
      integer i,j
      double precision fac1, fac2, v1(0:1), v2(0:1)
      
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0

        fac1 = one / dsqrt( dXdxi(i,j,0)**2 + dXdxi(i,j,2)**2 )
        fac2 = one / dsqrt( dXdxi(i,j,1)**2 + dXdxi(i,j,3)**2 )
        v1(0) = dXdxi(i,j,0) * fac1
        v1(1) = dXdxi(i,j,2) * fac1
        v2(0) =  dXdxi(i,j,3) * fac2
        v2(1) = -dXdxi(i,j,1) * fac2
        data(i,j) = fac1 * ( v1(0)*v2(0) + v1(1)*v2(1) )
      
      enddo
      enddo
      return
      end
      subroutine MAG_BLOCK_PROJECT_2D(
     &           iboxlo0,iboxlo1
     &           ,iboxhi0,iboxhi1
     &           ,vec_src
     &           ,ivec_srclo0,ivec_srclo1
     &           ,ivec_srchi0,ivec_srchi1
     &           ,nvec_srccomp
     &           ,vec_dst
     &           ,ivec_dstlo0,ivec_dstlo1
     &           ,ivec_dsthi0,ivec_dsthi1
     &           ,nvec_dstcomp
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer nvec_srccomp
      integer ivec_srclo0,ivec_srclo1
      integer ivec_srchi0,ivec_srchi1
      REAL_T vec_src(
     &           ivec_srclo0:ivec_srchi0,
     &           ivec_srclo1:ivec_srchi1,
     &           0:nvec_srccomp-1)
      integer nvec_dstcomp
      integer ivec_dstlo0,ivec_dstlo1
      integer ivec_dsthi0,ivec_dsthi1
      REAL_T vec_dst(
     &           ivec_dstlo0:ivec_dsthi0,
     &           ivec_dstlo1:ivec_dsthi1,
     &           0:nvec_dstcomp-1)
      integer i,j, comp, ncomp
      double precision dotprod
      ncomp = nvec_srccomp
      
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0

        dotprod = zero
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
     &           iboxlo0,iboxlo1
     &           ,iboxhi0,iboxhi1
     &           ,vec_psi
     &           ,ivec_psilo0,ivec_psilo1
     &           ,ivec_psihi0,ivec_psihi1
     &           ,nvec_psicomp
     &           ,vec_theta
     &           ,ivec_thetalo0,ivec_thetalo1
     &           ,ivec_thetahi0,ivec_thetahi1
     &           ,nvec_thetacomp
     &           ,vec_dst
     &           ,ivec_dstlo0,ivec_dstlo1
     &           ,ivec_dsthi0,ivec_dsthi1
     &           ,nvec_dstcomp
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer nvec_psicomp
      integer ivec_psilo0,ivec_psilo1
      integer ivec_psihi0,ivec_psihi1
      REAL_T vec_psi(
     &           ivec_psilo0:ivec_psihi0,
     &           ivec_psilo1:ivec_psihi1,
     &           0:nvec_psicomp-1)
      integer nvec_thetacomp
      integer ivec_thetalo0,ivec_thetalo1
      integer ivec_thetahi0,ivec_thetahi1
      REAL_T vec_theta(
     &           ivec_thetalo0:ivec_thetahi0,
     &           ivec_thetalo1:ivec_thetahi1,
     &           0:nvec_thetacomp-1)
      integer nvec_dstcomp
      integer ivec_dstlo0,ivec_dstlo1
      integer ivec_dsthi0,ivec_dsthi1
      REAL_T vec_dst(
     &           ivec_dstlo0:ivec_dsthi0,
     &           ivec_dstlo1:ivec_dsthi1,
     &           0:nvec_dstcomp-1)
      integer i,j, comp, ncomp
      double precision dotprod_psi, dotprod_theta
      ncomp = nvec_dstcomp
      
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0

        dotprod_psi = zero
        do comp = 0, ncomp - 1
          dotprod_psi = dotprod_psi + vec_psi(i,j,comp) * vec_dst(i,j,comp)
        enddo     
        dotprod_theta = zero
        do comp = 0, ncomp - 1
          dotprod_theta = dotprod_theta + vec_theta(i,j,comp) * vec_dst(i,j,comp)
        enddo     
        vec_dst(i,j,0) = dotprod_psi
        vec_dst(i,j,1) = dotprod_theta
      
      enddo
      enddo
      return
      end
      subroutine GET_NODAL_FIELD_DATA_2D(
     &           igridboxlo0,igridboxlo1
     &           ,igridboxhi0,igridboxhi1
     &           ,RZ
     &           ,iRZlo0,iRZlo1
     &           ,iRZhi0,iRZhi1
     &           ,nRZcomp
     &           ,psi
     &           ,ipsilo0,ipsilo1
     &           ,ipsihi0,ipsihi1
     &           ,RB
     &           ,iRBlo0,iRBlo1
     &           ,iRBhi0,iRBhi1
     &           ,nRBcomp
     &           ,RBtor
     &           ,A
     &           ,iAlo0,iAlo1
     &           ,iAhi0,iAhi1
     &           ,nAcomp
     &           ,bunit
     &           ,ibunitlo0,ibunitlo1
     &           ,ibunithi0,ibunithi1
     &           ,nbunitcomp
     &           ,Bmag
     &           ,iBmaglo0,iBmaglo1
     &           ,iBmaghi0,iBmaghi1
     &           ,nBmagcomp
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer igridboxlo0,igridboxlo1
      integer igridboxhi0,igridboxhi1
      integer nRZcomp
      integer iRZlo0,iRZlo1
      integer iRZhi0,iRZhi1
      REAL_T RZ(
     &           iRZlo0:iRZhi0,
     &           iRZlo1:iRZhi1,
     &           0:nRZcomp-1)
      integer ipsilo0,ipsilo1
      integer ipsihi0,ipsihi1
      REAL_T psi(
     &           ipsilo0:ipsihi0,
     &           ipsilo1:ipsihi1)
      integer nRBcomp
      integer iRBlo0,iRBlo1
      integer iRBhi0,iRBhi1
      REAL_T RB(
     &           iRBlo0:iRBhi0,
     &           iRBlo1:iRBhi1,
     &           0:nRBcomp-1)
      REAL_T RBtor
      integer nAcomp
      integer iAlo0,iAlo1
      integer iAhi0,iAhi1
      REAL_T A(
     &           iAlo0:iAhi0,
     &           iAlo1:iAhi1,
     &           0:nAcomp-1)
      integer nbunitcomp
      integer ibunitlo0,ibunitlo1
      integer ibunithi0,ibunithi1
      REAL_T bunit(
     &           ibunitlo0:ibunithi0,
     &           ibunitlo1:ibunithi1,
     &           0:nbunitcomp-1)
      integer nBmagcomp
      integer iBmaglo0,iBmaglo1
      integer iBmaghi0,iBmaghi1
      REAL_T Bmag(
     &           iBmaglo0:iBmaghi0,
     &           iBmaglo1:iBmaghi1,
     &           0:nBmagcomp-1)
      integer i,j, n
      double precision R, Z, Bmagnitude, B(0:2)
      
      do j = igridboxlo1,igridboxhi1
      do i = igridboxlo0,igridboxhi0

         R = RZ(i,j,0);
         Z = RZ(i,j,1);
         B(0) = RB(i,j,0) / R
         B(1) = RBtor / R
         B(2) = RB(i,j,1) / R
         Bmagnitude = zero
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
         A(i,j,2) = zero
      
      enddo
      enddo
      return
      end
      subroutine GET_FIELD_DATA_2D(
     &           igridboxlo0,igridboxlo1
     &           ,igridboxhi0,igridboxhi1
     &           ,RZ
     &           ,iRZlo0,iRZlo1
     &           ,iRZhi0,iRZhi1
     &           ,nRZcomp
     &           ,RB
     &           ,iRBlo0,iRBlo1
     &           ,iRBhi0,iRBhi1
     &           ,nRBcomp
     &           ,dRBdR
     &           ,idRBdRlo0,idRBdRlo1
     &           ,idRBdRhi0,idRBdRhi1
     &           ,ndRBdRcomp
     &           ,dRBdZ
     &           ,idRBdZlo0,idRBdZlo1
     &           ,idRBdZhi0,idRBdZhi1
     &           ,ndRBdZcomp
     &           ,RBtor
     &           ,B
     &           ,iBlo0,iBlo1
     &           ,iBhi0,iBhi1
     &           ,nBcomp
     &           ,Bmag
     &           ,iBmaglo0,iBmaglo1
     &           ,iBmaghi0,iBmaghi1
     &           ,bunit
     &           ,ibunitlo0,ibunitlo1
     &           ,ibunithi0,ibunithi1
     &           ,nbunitcomp
     &           ,gradB
     &           ,igradBlo0,igradBlo1
     &           ,igradBhi0,igradBhi1
     &           ,ngradBcomp
     &           ,curlb
     &           ,icurlblo0,icurlblo1
     &           ,icurlbhi0,icurlbhi1
     &           ,ncurlbcomp
     &           ,bdotcurlb
     &           ,ibdotcurlblo0,ibdotcurlblo1
     &           ,ibdotcurlbhi0,ibdotcurlbhi1
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer igridboxlo0,igridboxlo1
      integer igridboxhi0,igridboxhi1
      integer nRZcomp
      integer iRZlo0,iRZlo1
      integer iRZhi0,iRZhi1
      REAL_T RZ(
     &           iRZlo0:iRZhi0,
     &           iRZlo1:iRZhi1,
     &           0:nRZcomp-1)
      integer nRBcomp
      integer iRBlo0,iRBlo1
      integer iRBhi0,iRBhi1
      REAL_T RB(
     &           iRBlo0:iRBhi0,
     &           iRBlo1:iRBhi1,
     &           0:nRBcomp-1)
      integer ndRBdRcomp
      integer idRBdRlo0,idRBdRlo1
      integer idRBdRhi0,idRBdRhi1
      REAL_T dRBdR(
     &           idRBdRlo0:idRBdRhi0,
     &           idRBdRlo1:idRBdRhi1,
     &           0:ndRBdRcomp-1)
      integer ndRBdZcomp
      integer idRBdZlo0,idRBdZlo1
      integer idRBdZhi0,idRBdZhi1
      REAL_T dRBdZ(
     &           idRBdZlo0:idRBdZhi0,
     &           idRBdZlo1:idRBdZhi1,
     &           0:ndRBdZcomp-1)
      REAL_T RBtor
      integer nBcomp
      integer iBlo0,iBlo1
      integer iBhi0,iBhi1
      REAL_T B(
     &           iBlo0:iBhi0,
     &           iBlo1:iBhi1,
     &           0:nBcomp-1)
      integer iBmaglo0,iBmaglo1
      integer iBmaghi0,iBmaghi1
      REAL_T Bmag(
     &           iBmaglo0:iBmaghi0,
     &           iBmaglo1:iBmaghi1)
      integer nbunitcomp
      integer ibunitlo0,ibunitlo1
      integer ibunithi0,ibunithi1
      REAL_T bunit(
     &           ibunitlo0:ibunithi0,
     &           ibunitlo1:ibunithi1,
     &           0:nbunitcomp-1)
      integer ngradBcomp
      integer igradBlo0,igradBlo1
      integer igradBhi0,igradBhi1
      REAL_T gradB(
     &           igradBlo0:igradBhi0,
     &           igradBlo1:igradBhi1,
     &           0:ngradBcomp-1)
      integer ncurlbcomp
      integer icurlblo0,icurlblo1
      integer icurlbhi0,icurlbhi1
      REAL_T curlb(
     &           icurlblo0:icurlbhi0,
     &           icurlblo1:icurlbhi1,
     &           0:ncurlbcomp-1)
      integer ibdotcurlblo0,ibdotcurlblo1
      integer ibdotcurlbhi0,ibdotcurlbhi1
      REAL_T bdotcurlb(
     &           ibdotcurlblo0:ibdotcurlbhi0,
     &           ibdotcurlblo1:ibdotcurlbhi1)
      integer i,j, n
      double precision R, Z, Bmag2
      
      do j = igridboxlo1,igridboxhi1
      do i = igridboxlo0,igridboxhi0

         R = RZ(i,j,0);
         Z = RZ(i,j,1);
         B(i,j,0) = RB(i,j,0) / R
         B(i,j,1) = RBtor / R
         B(i,j,2) = RB(i,j,1) / R
         Bmag2 = zero
         do n = 0, 2
            Bmag2 = Bmag2 + B(i,j,n)**2
         enddo
         Bmag(i,j) = sqrt(Bmag2)
         do n = 0, 2
            bunit(i,j,n) = B(i,j,n) / Bmag(i,j)
         enddo
         gradB(i,j,0) = (B(i,j,0) * dRBdR(i,j,0) + B(i,j,2) * dRBdR(i,j,1))
     &                            / (R * Bmag(i,j)) - Bmag(i,j) / R
         gradB(i,j,1) = zero
         gradB(i,j,2) = (B(i,j,0) * dRBdZ(i,j,0) + B(i,j,2) * dRBdZ(i,j,1))
     &                            / (R * Bmag(i,j))
         curlb(i,j,0) =  bunit(i,j,1) * gradB(i,j,2) / Bmag(i,j)
         curlb(i,j,1) = (dRBdZ(i,j,0) - dRBdR(i,j,1)) / (R * Bmag(i,j))
     &                          + (bunit(i,j,2) * gradB(i,j,0) 
     &                           - bunit(i,j,0) * gradB(i,j,2)) / Bmag(i,j)
     &                          + bunit(i,j,2) / R
         curlb(i,j,2) = - bunit(i,j,1) * gradB(i,j,0) / Bmag(i,j)
         bdotcurlb(i,j) = zero
         do n = 0, 2
            bdotcurlb(i,j) = bdotcurlb(i,j) + bunit(i,j,n) * curlb(i,j,n) 
         enddo
      
      enddo
      enddo
      return
      end
