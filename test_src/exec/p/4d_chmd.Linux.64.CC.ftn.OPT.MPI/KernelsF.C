#include "REAL.H"
#include "SPACE.H"
#include "CONSTANTS.H"

#include "CONSTANTS.H"
      subroutine COMPUTE_VEL_CELL_4D(
     &           velCell
     &           ,ivelCelllo0,ivelCelllo1,ivelCelllo2,ivelCelllo3
     &           ,ivelCellhi0,ivelCellhi1,ivelCellhi2,ivelCellhi3
     &           ,nvelCellcomp
     &           ,velFace
     &           ,ivelFacelo0,ivelFacelo1,ivelFacelo2,ivelFacelo3
     &           ,ivelFacehi0,ivelFacehi1,ivelFacehi2,ivelFacehi3
     &           ,nvelFacecomp
     &           ,dir
     &           ,igridboxlo0,igridboxlo1,igridboxlo2,igridboxlo3
     &           ,igridboxhi0,igridboxhi1,igridboxhi2,igridboxhi3
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer nvelCellcomp
      integer ivelCelllo0,ivelCelllo1,ivelCelllo2,ivelCelllo3
      integer ivelCellhi0,ivelCellhi1,ivelCellhi2,ivelCellhi3
      REAL_T velCell(
     &           ivelCelllo0:ivelCellhi0,
     &           ivelCelllo1:ivelCellhi1,
     &           ivelCelllo2:ivelCellhi2,
     &           ivelCelllo3:ivelCellhi3,
     &           0:nvelCellcomp-1)
      integer nvelFacecomp
      integer ivelFacelo0,ivelFacelo1,ivelFacelo2,ivelFacelo3
      integer ivelFacehi0,ivelFacehi1,ivelFacehi2,ivelFacehi3
      REAL_T velFace(
     &           ivelFacelo0:ivelFacehi0,
     &           ivelFacelo1:ivelFacehi1,
     &           ivelFacelo2:ivelFacehi2,
     &           ivelFacelo3:ivelFacehi3,
     &           0:nvelFacecomp-1)
      integer dir
      integer igridboxlo0,igridboxlo1,igridboxlo2,igridboxlo3
      integer igridboxhi0,igridboxhi1,igridboxhi2,igridboxhi3
      integer i,j,k,l, ii,jj,kk,ll, cfg_dim, comp
      double precision fac
#if CH_SPACEDIM==4
      cfg_dim = 2
      fac = 1.0/4.0
#else
      cfg_dim = 3
      fac = 1.0/6.0	
#endif
      
      ii = CHF_ID(0,dir)
      jj = CHF_ID(1,dir)
      kk = CHF_ID(2,dir)
      ll = CHF_ID(3,dir)
      
      do l = igridboxlo3,igridboxhi3
      do k = igridboxlo2,igridboxhi2
      do j = igridboxlo1,igridboxhi1
      do i = igridboxlo0,igridboxhi0

         do comp = 0, cfg_dim - 1
           velCell(i,j,k,l,comp) = velCell(i,j,k,l,comp) 
     &                                     + fac * velFace(i,j,k,l,comp)
     &                                     + fac * velFace(i+ ii,j+ jj,k+ kk,l+ ll,comp) 
         enddo
      
      enddo
      enddo
      enddo
      enddo
      return
      end
      subroutine COMPUTE_PERP_VEL_4D(
     &           iboxlo0,iboxlo1,iboxlo2,iboxlo3
     &           ,iboxhi0,iboxhi1,iboxhi2,iboxhi3
     &           ,result
     &           ,iresultlo0,iresultlo1,iresultlo2,iresultlo3
     &           ,iresulthi0,iresulthi1,iresulthi2,iresulthi3
     &           ,nresultcomp
     &           ,dfn
     &           ,idfnlo0,idfnlo1,idfnlo2,idfnlo3
     &           ,idfnhi0,idfnhi1,idfnhi2,idfnhi3
     &           ,gkVel
     &           ,igkVello0,igkVello1,igkVello2,igkVello3
     &           ,igkVelhi0,igkVelhi1,igkVelhi2,igkVelhi3
     &           ,ngkVelcomp
     &           ,velCoords
     &           ,ivelCoordslo0,ivelCoordslo1,ivelCoordslo2,ivelCoordslo3
     &           ,ivelCoordshi0,ivelCoordshi1,ivelCoordshi2,ivelCoordshi3
     &           ,nvelCoordscomp
     &           ,b
     &           ,iblo0,iblo1,iblo2,iblo3
     &           ,ibhi0,ibhi1,ibhi2,ibhi3
     &           ,nbcomp
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iboxlo0,iboxlo1,iboxlo2,iboxlo3
      integer iboxhi0,iboxhi1,iboxhi2,iboxhi3
      integer nresultcomp
      integer iresultlo0,iresultlo1,iresultlo2,iresultlo3
      integer iresulthi0,iresulthi1,iresulthi2,iresulthi3
      REAL_T result(
     &           iresultlo0:iresulthi0,
     &           iresultlo1:iresulthi1,
     &           iresultlo2:iresulthi2,
     &           iresultlo3:iresulthi3,
     &           0:nresultcomp-1)
      integer idfnlo0,idfnlo1,idfnlo2,idfnlo3
      integer idfnhi0,idfnhi1,idfnhi2,idfnhi3
      REAL_T dfn(
     &           idfnlo0:idfnhi0,
     &           idfnlo1:idfnhi1,
     &           idfnlo2:idfnhi2,
     &           idfnlo3:idfnhi3)
      integer ngkVelcomp
      integer igkVello0,igkVello1,igkVello2,igkVello3
      integer igkVelhi0,igkVelhi1,igkVelhi2,igkVelhi3
      REAL_T gkVel(
     &           igkVello0:igkVelhi0,
     &           igkVello1:igkVelhi1,
     &           igkVello2:igkVelhi2,
     &           igkVello3:igkVelhi3,
     &           0:ngkVelcomp-1)
      integer nvelCoordscomp
      integer ivelCoordslo0,ivelCoordslo1,ivelCoordslo2,ivelCoordslo3
      integer ivelCoordshi0,ivelCoordshi1,ivelCoordshi2,ivelCoordshi3
      REAL_T velCoords(
     &           ivelCoordslo0:ivelCoordshi0,
     &           ivelCoordslo1:ivelCoordshi1,
     &           ivelCoordslo2:ivelCoordshi2,
     &           ivelCoordslo3:ivelCoordshi3,
     &           0:nvelCoordscomp-1)
      integer nbcomp
      integer iblo0,iblo1,iblo2,iblo3
      integer ibhi0,ibhi1,ibhi2,ibhi3
      REAL_T b(
     &           iblo0:ibhi0,
     &           iblo1:ibhi1,
     &           iblo2:ibhi2,
     &           iblo3:ibhi3,
     &           0:nbcomp-1)
      integer i,j,k,l 
      double precision bR, bphi, bZ, vpar_R, vpar_Phi, vpar_Z, vperp_R, vperp_Phi, vperp_Z 
      
      do l = iboxlo3,iboxhi3
      do k = iboxlo2,iboxhi2
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0

#if CH_SPACEDIM==4
          bR = b(i,j,iblo2,iblo3,0) 
          bZ = b(i,j,iblo2,iblo3,2) 
          vpar_R = velCoords(i,j,k,l,0) * bR
          vpar_Z = velCoords(i,j,k,l,0) * bZ
          vperp_R =  gkVel(i,j,k,l,0) - vpar_R    
          vperp_Z =  gkVel(i,j,k,l,1) - vpar_Z
          result(i,j,k,l,0) = vperp_R * dfn(i,j,k,l)
          result(i,j,k,l,1) = vperp_Z * dfn(i,j,k,l)
#else        
          bR = b(i,j,k,iblo2,iblo3,0) 
          bZ = b(i,j,k,iblo2,iblo3,2)
          bphi =  b(i,j,k,iblo2,iblo3,1)
          vpar_R   = velCoords(i,j,k,l,0) * bR
          vpar_Phi = velCoords(i,j,k,l,0) * bPhi
          vpar_Z   = velCoords(i,j,k,l,0) * bZ
          vperp_R   =  gkVel(i,j,k,l,0) - vpar_R    
          vperp_Phi =  gkVel(i,j,k,l,1) - vpar_Phi
          vperp_Z   =  gkVel(i,j,k,l,2) - vpar_Z
          result(i,j,k,l,0) = vperp_R   * dfn(i,j,k,l)
          result(i,j,k,l,1) = vperp_Phi * dfn(i,j,k,l)
          result(i,j,k,l,2) = vperp_Z   * dfn(i,j,k,l)
#endif
      
      enddo
      enddo
      enddo
      enddo
      return
      end
      subroutine COMPUTE_PAR_MOM_4D(
     &           result
     &           ,iresultlo0,iresultlo1,iresultlo2,iresultlo3
     &           ,iresulthi0,iresulthi1,iresulthi2,iresulthi3
     &           ,nresultcomp
     &           ,parVel
     &           ,iparVello0,iparVello1,iparVello2,iparVello3
     &           ,iparVelhi0,iparVelhi1,iparVelhi2,iparVelhi3
     &           ,igridboxlo0,igridboxlo1,igridboxlo2,igridboxlo3
     &           ,igridboxhi0,igridboxhi1,igridboxhi2,igridboxhi3
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nresultcomp
      integer iresultlo0,iresultlo1,iresultlo2,iresultlo3
      integer iresulthi0,iresulthi1,iresulthi2,iresulthi3
      REAL_T result(
     &           iresultlo0:iresulthi0,
     &           iresultlo1:iresulthi1,
     &           iresultlo2:iresulthi2,
     &           iresultlo3:iresulthi3,
     &           0:nresultcomp-1)
      integer iparVello0,iparVello1,iparVello2,iparVello3
      integer iparVelhi0,iparVelhi1,iparVelhi2,iparVelhi3
      REAL_T parVel(
     &           iparVello0:iparVelhi0,
     &           iparVello1:iparVelhi1,
     &           iparVello2:iparVelhi2,
     &           iparVello3:iparVelhi3)
      integer igridboxlo0,igridboxlo1,igridboxlo2,igridboxlo3
      integer igridboxhi0,igridboxhi1,igridboxhi2,igridboxhi3
      integer i,j,k,l, comp
      
      do l = igridboxlo3,igridboxhi3
      do k = igridboxlo2,igridboxhi2
      do j = igridboxlo1,igridboxhi1
      do i = igridboxlo0,igridboxhi0

         do comp = 0, nresultcomp-1
           result(i,j,k,l,comp) = result(i,j,k,l,comp) * parVel(i,j,k,l)
         enddo
      
      enddo
      enddo
      enddo
      enddo
      return
      end
      subroutine COMPUTE_PRESSURE_4D(
     &           result
     &           ,iresultlo0,iresultlo1,iresultlo2,iresultlo3
     &           ,iresulthi0,iresulthi1,iresulthi2,iresulthi3
     &           ,nresultcomp
     &           ,velCoords
     &           ,ivelCoordslo0,ivelCoordslo1,ivelCoordslo2,ivelCoordslo3
     &           ,ivelCoordshi0,ivelCoordshi1,ivelCoordshi2,ivelCoordshi3
     &           ,nvelCoordscomp
     &           ,vparShift
     &           ,ivparShiftlo0,ivparShiftlo1,ivparShiftlo2,ivparShiftlo3
     &           ,ivparShifthi0,ivparShifthi1,ivparShifthi2,ivparShifthi3
     &           ,B
     &           ,iBlo0,iBlo1,iBlo2,iBlo3
     &           ,iBhi0,iBhi1,iBhi2,iBhi3
     &           ,mass
     &           ,igridboxlo0,igridboxlo1,igridboxlo2,igridboxlo3
     &           ,igridboxhi0,igridboxhi1,igridboxhi2,igridboxhi3
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nresultcomp
      integer iresultlo0,iresultlo1,iresultlo2,iresultlo3
      integer iresulthi0,iresulthi1,iresulthi2,iresulthi3
      REAL_T result(
     &           iresultlo0:iresulthi0,
     &           iresultlo1:iresulthi1,
     &           iresultlo2:iresulthi2,
     &           iresultlo3:iresulthi3,
     &           0:nresultcomp-1)
      integer nvelCoordscomp
      integer ivelCoordslo0,ivelCoordslo1,ivelCoordslo2,ivelCoordslo3
      integer ivelCoordshi0,ivelCoordshi1,ivelCoordshi2,ivelCoordshi3
      REAL_T velCoords(
     &           ivelCoordslo0:ivelCoordshi0,
     &           ivelCoordslo1:ivelCoordshi1,
     &           ivelCoordslo2:ivelCoordshi2,
     &           ivelCoordslo3:ivelCoordshi3,
     &           0:nvelCoordscomp-1)
      integer ivparShiftlo0,ivparShiftlo1,ivparShiftlo2,ivparShiftlo3
      integer ivparShifthi0,ivparShifthi1,ivparShifthi2,ivparShifthi3
      REAL_T vparShift(
     &           ivparShiftlo0:ivparShifthi0,
     &           ivparShiftlo1:ivparShifthi1,
     &           ivparShiftlo2:ivparShifthi2,
     &           ivparShiftlo3:ivparShifthi3)
      integer iBlo0,iBlo1,iBlo2,iBlo3
      integer iBhi0,iBhi1,iBhi2,iBhi3
      REAL_T B(
     &           iBlo0:iBhi0,
     &           iBlo1:iBhi1,
     &           iBlo2:iBhi2,
     &           iBlo3:iBhi3)
      REAL_T mass
      integer igridboxlo0,igridboxlo1,igridboxlo2,igridboxlo3
      integer igridboxhi0,igridboxhi1,igridboxhi2,igridboxhi3
      integer i,j,k,l, comp
      double precision vperp2, vparLoc, v2
      
      do l = igridboxlo3,igridboxhi3
      do k = igridboxlo2,igridboxhi2
      do j = igridboxlo1,igridboxhi1
      do i = igridboxlo0,igridboxhi0

#if CH_SPACEDIM==4
         vperp2 = velCoords(i,j,k,l,1) * B(i,j,iblo2,iblo3)
         vparLoc = velCoords(i,j,k,l,0) - vparShift(i,j,ivparShiftlo2,ivparShiftlo3)
#else
         vperp2 = velCoords(i,j,k,l,1) * B(i,j,k,iblo3,iblo4)
         vparLoc = velCoords(i,j,k,l,0) - vparShift(i,j,k,ivparShiftlo3,ivparShiftlo4)
#endif
         v2 = mass * vparLoc**2 + vperp2	 
         do comp = 0, nresultcomp-1
           result(i,j,k,l,comp) = result(i,j,k,l,comp) * v2 / 3.0
         enddo
      
      enddo
      enddo
      enddo
      enddo
      return
      end
      subroutine COMPUTE_MAXWELLIAN_4D(
     &           result
     &           ,iresultlo0,iresultlo1,iresultlo2,iresultlo3
     &           ,iresulthi0,iresulthi1,iresulthi2,iresulthi3
     &           ,nresultcomp
     &           ,velCoords
     &           ,ivelCoordslo0,ivelCoordslo1,ivelCoordslo2,ivelCoordslo3
     &           ,ivelCoordshi0,ivelCoordshi1,ivelCoordshi2,ivelCoordshi3
     &           ,nvelCoordscomp
     &           ,vparShift
     &           ,ivparShiftlo0,ivparShiftlo1,ivparShiftlo2,ivparShiftlo3
     &           ,ivparShifthi0,ivparShifthi1,ivparShifthi2,ivparShifthi3
     &           ,n
     &           ,inlo0,inlo1,inlo2,inlo3
     &           ,inhi0,inhi1,inhi2,inhi3
     &           ,T
     &           ,iTlo0,iTlo1,iTlo2,iTlo3
     &           ,iThi0,iThi1,iThi2,iThi3
     &           ,B
     &           ,iBlo0,iBlo1,iBlo2,iBlo3
     &           ,iBhi0,iBhi1,iBhi2,iBhi3
     &           ,mass
     &           ,igridboxlo0,igridboxlo1,igridboxlo2,igridboxlo3
     &           ,igridboxhi0,igridboxhi1,igridboxhi2,igridboxhi3
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nresultcomp
      integer iresultlo0,iresultlo1,iresultlo2,iresultlo3
      integer iresulthi0,iresulthi1,iresulthi2,iresulthi3
      REAL_T result(
     &           iresultlo0:iresulthi0,
     &           iresultlo1:iresulthi1,
     &           iresultlo2:iresulthi2,
     &           iresultlo3:iresulthi3,
     &           0:nresultcomp-1)
      integer nvelCoordscomp
      integer ivelCoordslo0,ivelCoordslo1,ivelCoordslo2,ivelCoordslo3
      integer ivelCoordshi0,ivelCoordshi1,ivelCoordshi2,ivelCoordshi3
      REAL_T velCoords(
     &           ivelCoordslo0:ivelCoordshi0,
     &           ivelCoordslo1:ivelCoordshi1,
     &           ivelCoordslo2:ivelCoordshi2,
     &           ivelCoordslo3:ivelCoordshi3,
     &           0:nvelCoordscomp-1)
      integer ivparShiftlo0,ivparShiftlo1,ivparShiftlo2,ivparShiftlo3
      integer ivparShifthi0,ivparShifthi1,ivparShifthi2,ivparShifthi3
      REAL_T vparShift(
     &           ivparShiftlo0:ivparShifthi0,
     &           ivparShiftlo1:ivparShifthi1,
     &           ivparShiftlo2:ivparShifthi2,
     &           ivparShiftlo3:ivparShifthi3)
      integer inlo0,inlo1,inlo2,inlo3
      integer inhi0,inhi1,inhi2,inhi3
      REAL_T n(
     &           inlo0:inhi0,
     &           inlo1:inhi1,
     &           inlo2:inhi2,
     &           inlo3:inhi3)
      integer iTlo0,iTlo1,iTlo2,iTlo3
      integer iThi0,iThi1,iThi2,iThi3
      REAL_T T(
     &           iTlo0:iThi0,
     &           iTlo1:iThi1,
     &           iTlo2:iThi2,
     &           iTlo3:iThi3)
      integer iBlo0,iBlo1,iBlo2,iBlo3
      integer iBhi0,iBhi1,iBhi2,iBhi3
      REAL_T B(
     &           iBlo0:iBhi0,
     &           iBlo1:iBhi1,
     &           iBlo2:iBhi2,
     &           iBlo3:iBhi3)
      REAL_T mass
      integer igridboxlo0,igridboxlo1,igridboxlo2,igridboxlo3
      integer igridboxhi0,igridboxhi1,igridboxhi2,igridboxhi3
      integer i,j,k,l, comp
      double precision v_parallel, mu, eparnorm, munorm, val, factor
      factor = sqrt(PI * (two/mass)**three);
      
      do l = igridboxlo3,igridboxhi3
      do k = igridboxlo2,igridboxhi2
      do j = igridboxlo1,igridboxhi1
      do i = igridboxlo0,igridboxhi0

         v_parallel = velCoords(i,j,k,l,0)
         mu = velCoords(i,j,k,l,1)
#if CH_SPACEDIM==4
         eparnorm = (one/two) * mass * (v_parallel-vparShift(i,j,ivparShiftlo2,ivparShiftlo3))**two / T(i,j,iTlo2,iTlo3)
         munorm   = (one/two) * B(i,j,iblo2,iblo3) * mu / T(i,j,iTlo2,iTlo3)
         val      = exp( -( eparnorm + munorm ) )
         val    = val * n(i,j,inlo2,inlo3) / ( factor * T(i,j,iTlo2,iTlo3)**(three/two) )
#else
         eparnorm = (one/two) * mass * (v_parallel-vparShift(i,j,k,ivparShiftlo3,ivparShiftlo4))**two / T(i,j,k,iTlo3,iTlo4)
         munorm   = (one/two) * B(i,j,k,iblo3,iblo4) * mu / T(i,j,k,iTlo3,iTlo4)
         val      = exp( -( eparnorm + munorm ) )
         val    = val * n(i,j,k,inlo3,inlo4) / ( factor * T(i,j,k,iTlo3,iTlo4)**(three/two) )
#endif
         do comp = 0, nresultcomp-1
           result(i,j,k,l,comp) = val
         enddo
      
      enddo
      enddo
      enddo
      enddo
      return
      end
