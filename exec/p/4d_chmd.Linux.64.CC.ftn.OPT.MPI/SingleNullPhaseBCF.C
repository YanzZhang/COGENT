#include "REAL.H"
#include "SPACE.H"
#include "CONSTANTS.H"

#include "CONSTANTS.H"
      subroutine SET_LOGICAL_SHEATH_BC_4D(
     &           f
     &           ,iflo0,iflo1,iflo2,iflo3
     &           ,ifhi0,ifhi1,ifhi2,ifhi3
     &           ,ibdryboxlo0,ibdryboxlo1,ibdryboxlo2,ibdryboxlo3
     &           ,ibdryboxhi0,ibdryboxhi1,ibdryboxhi2,ibdryboxhi3
     &           ,f_rflct
     &           ,if_rflctlo0,if_rflctlo1,if_rflctlo2,if_rflctlo3
     &           ,if_rflcthi0,if_rflcthi1,if_rflcthi2,if_rflcthi3
     &           ,vel
     &           ,ivello0,ivello1,ivello2,ivello3
     &           ,ivelhi0,ivelhi1,ivelhi2,ivelhi3
     &           ,nvelcomp
     &           ,vn
     &           ,ivnlo0,ivnlo1,ivnlo2,ivnlo3
     &           ,ivnhi0,ivnhi1,ivnhi2,ivnhi3
     &           ,phi
     &           ,iphilo0,iphilo1,iphilo2,iphilo3
     &           ,iphihi0,iphihi1,iphihi2,iphihi3
     &           ,mass
     &           ,charge
     &           ,iside
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
      integer ibdryboxlo0,ibdryboxlo1,ibdryboxlo2,ibdryboxlo3
      integer ibdryboxhi0,ibdryboxhi1,ibdryboxhi2,ibdryboxhi3
      integer if_rflctlo0,if_rflctlo1,if_rflctlo2,if_rflctlo3
      integer if_rflcthi0,if_rflcthi1,if_rflcthi2,if_rflcthi3
      REAL_T f_rflct(
     &           if_rflctlo0:if_rflcthi0,
     &           if_rflctlo1:if_rflcthi1,
     &           if_rflctlo2:if_rflcthi2,
     &           if_rflctlo3:if_rflcthi3)
      integer nvelcomp
      integer ivello0,ivello1,ivello2,ivello3
      integer ivelhi0,ivelhi1,ivelhi2,ivelhi3
      REAL_T vel(
     &           ivello0:ivelhi0,
     &           ivello1:ivelhi1,
     &           ivello2:ivelhi2,
     &           ivello3:ivelhi3,
     &           0:nvelcomp-1)
      integer ivnlo0,ivnlo1,ivnlo2,ivnlo3
      integer ivnhi0,ivnhi1,ivnhi2,ivnhi3
      REAL_T vn(
     &           ivnlo0:ivnhi0,
     &           ivnlo1:ivnhi1,
     &           ivnlo2:ivnhi2,
     &           ivnlo3:ivnhi3)
      integer iphilo0,iphilo1,iphilo2,iphilo3
      integer iphihi0,iphihi1,iphihi2,iphihi3
      REAL_T phi(
     &           iphilo0:iphihi0,
     &           iphilo1:iphihi1,
     &           iphilo2:iphihi2,
     &           iphilo3:iphihi3)
      REAL_T mass
      REAL_T charge
      integer iside
      integer i,j,k,l,m
      integer isign
      integer vn_jbdry, phi_jbdry
      integer jbdry,jsrc,jsrc_offset
      real  pot_energy_min
      if (iside.eq.0) then
        jbdry = ibdryboxhi1
      else
        jbdry = ibdryboxlo1
      endif
      isign = 2*iside-1
      jsrc_offset = 2 * jbdry - isign
      
      do l = ibdryboxlo3,ibdryboxhi3
      do k = ibdryboxlo2,ibdryboxhi2
      do j = ibdryboxlo1,ibdryboxhi1
      do i = ibdryboxlo0,ibdryboxhi0

          vn_jbdry = jbdry - isign + iside
#if CH_SPACEDIM==4
          if (isign*vn(i,vn_jbdry,k,l).le.zero) then
             phi_jbdry = jbdry - isign
             pot_energy_min = -charge * phi(i,phi_jbdry,iphilo2,iphilo3)
             if (mass * (vel(i,jbdry,k,l,0)**2).gt.pot_energy_min) then
               f(i,j,k,l) = zero
             else
               jsrc = jsrc_offset - j
               f(i,j,k,l) = f_rflct(i,jsrc,-k-1,l)
             endif
          endif
#else
          if (isign*vn(i,vn_jbdry,k,l,m).le.zero) then
             phi_jbdry = jbdry - isign
             pot_energy_min = -charge * phi(i,phi_jbdry,k,iphilo2,iphilo3)
             if (mass * (vel(i,jbdry,k,l,m,0)**2).gt.pot_energy_min) then
               f(i,j,k,l,m) = zero
             else
               jsrc = jsrc_offset - j
               f(i,j,k,l,m) = f_rflct(i,jsrc,k,-l-1,m)
             endif
          endif
#endif
      
      enddo
      enddo
      enddo
      enddo
      return
      end
      subroutine SET_POLOIDAL_DIVERTER_BC_4D(
     &           f
     &           ,iflo0,iflo1,iflo2,iflo3
     &           ,ifhi0,ifhi1,ifhi2,ifhi3
     &           ,nfcomp
     &           ,ibdryboxlo0,ibdryboxlo1,ibdryboxlo2,ibdryboxlo3
     &           ,ibdryboxhi0,ibdryboxhi1,ibdryboxhi2,ibdryboxhi3
     &           ,f_rflct
     &           ,if_rflctlo0,if_rflctlo1,if_rflctlo2,if_rflctlo3
     &           ,if_rflcthi0,if_rflcthi1,if_rflcthi2,if_rflcthi3
     &           ,nf_rflctcomp
     &           ,vn
     &           ,ivnlo0,ivnlo1,ivnlo2,ivnlo3
     &           ,ivnhi0,ivnhi1,ivnhi2,ivnhi3
     &           ,nvncomp
     &           ,qphi
     &           ,iqphilo0,iqphilo1,iqphilo2,iqphilo3
     &           ,iqphihi0,iqphihi1,iqphihi2,iqphihi3
     &           ,iside
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nfcomp
      integer iflo0,iflo1,iflo2,iflo3
      integer ifhi0,ifhi1,ifhi2,ifhi3
      REAL_T f(
     &           iflo0:ifhi0,
     &           iflo1:ifhi1,
     &           iflo2:ifhi2,
     &           iflo3:ifhi3,
     &           0:nfcomp-1)
      integer ibdryboxlo0,ibdryboxlo1,ibdryboxlo2,ibdryboxlo3
      integer ibdryboxhi0,ibdryboxhi1,ibdryboxhi2,ibdryboxhi3
      integer nf_rflctcomp
      integer if_rflctlo0,if_rflctlo1,if_rflctlo2,if_rflctlo3
      integer if_rflcthi0,if_rflcthi1,if_rflcthi2,if_rflcthi3
      REAL_T f_rflct(
     &           if_rflctlo0:if_rflcthi0,
     &           if_rflctlo1:if_rflcthi1,
     &           if_rflctlo2:if_rflcthi2,
     &           if_rflctlo3:if_rflcthi3,
     &           0:nf_rflctcomp-1)
      integer nvncomp
      integer ivnlo0,ivnlo1,ivnlo2,ivnlo3
      integer ivnhi0,ivnhi1,ivnhi2,ivnhi3
      REAL_T vn(
     &           ivnlo0:ivnhi0,
     &           ivnlo1:ivnhi1,
     &           ivnlo2:ivnhi2,
     &           ivnlo3:ivnhi3,
     &           0:nvncomp-1)
      integer iqphilo0,iqphilo1,iqphilo2,iqphilo3
      integer iqphihi0,iqphihi1,iqphihi2,iqphihi3
      REAL_T qphi(
     &           iqphilo0:iqphihi0,
     &           iqphilo1:iqphihi1,
     &           iqphilo2:iqphihi2,
     &           iqphilo3:iqphihi3)
      integer iside
      integer i,j,k,l,m
      integer isign
      integer n
      integer jbdry,jsrc,jsrc_offset
      REAL_T  pot_energy_min
      if (iside.eq.0) then
        jbdry = ibdryboxhi1
      else
        jbdry = ibdryboxlo1
      endif
      isign = 2*iside-1
      jsrc_offset = 2 * jbdry + isign
      do n=0,nfcomp-1
      
      do l = ibdryboxlo3,ibdryboxhi3
      do k = ibdryboxlo2,ibdryboxhi2
      do j = ibdryboxlo1,ibdryboxhi1
      do i = ibdryboxlo0,ibdryboxhi0

#if CH_SPACEDIM==5
          pot_energy_min = (-qphi(i,jbdry,k,iqphilo3,iqphilo4))
          if (isign*vn(i,jbdry,k,l,m,0).le.zero) then
            if ((vn(i,jbdry,k,l,m,0)**2).gt.pot_energy_min) then
              f(i,j,k,l,n,m) = zero
            else
              jsrc = jsrc_offset - j
              f(i,j,k,l,n,m) = f_rflct(i,jsrc,k,-l,m,n)
            endif
          endif
#else
          pot_energy_min = (-qphi(i,jbdry,iqphilo2,iqphilo3))
          if (isign*vn(i,jbdry,k,l,0).le.zero) then
            if ((vn(i,jbdry,k,l,0)**2).gt.pot_energy_min) then
              f(i,j,k,l,n) = zero
            else
              jsrc = jsrc_offset - j
              f(i,j,k,l,n) = f_rflct(i,jsrc,-k,l,n)
            endif
          endif
#endif
      
      enddo
      enddo
      enddo
      enddo
      enddo
      return
      end
