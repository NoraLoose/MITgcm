#include "ECCO_OPTIONS.h"

      subroutine ecco_phys( mythid )

c     ==================================================================
c     SUBROUTINE ecco_phys
c     ==================================================================
c
c     ==================================================================
c     SUBROUTINE ecco_phys
c     ==================================================================

      implicit none

c     == global variables ==

#include "EEPARAMS.h"
#include "SIZE.h"
#include "PARAMS.h"
#include "DYNVARS.h"
#include "FFIELDS.h"
#include "GRID.h"
#ifdef ALLOW_ECCO
# include "ecco.h"
#endif
#ifdef ALLOW_PTRACERS
# include "PTRACERS_SIZE.h"
# include "PTRACERS_FIELDS.h"
#endif

c     == routine arguments ==

      integer mythid

c     == local variables ==
C     localT           :: local copy of tracer field to compute gmredi transports
C     tmp2dmskW        :: copy of cell area mask to compute gmredi transports
C     tmp2dmskS        :: copy of cell area mask to compute gmredi transports
C     tmp3dmskW        :: copy of section mask
C     tmp3dmskS        :: copy of section mask
C     trIdentity       :: =1 for theta, =2 for salt
C     diffKh           :: diffKhT or diffKhS
C     df	           :: field for diffusive fluxes

      integer bi,bj
      integer i,j,k
      integer jmin,jmax
      integer imin,imax
#ifdef ALLOW_GENCOST_CONTRIBUTION
      integer kgen, kgen3d, itr
      integer trIdentity
      _RL diffKh
      _RL areavolTile(nSx,nSy), areavolGlob
      _RL tmpfld, tmpvol, tmpmsk, tmpmsk2, tmpmskW, tmpmskS
      _RL localT(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL tmp2dmskW(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL tmp2dmskS(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL tmp3dmskW(1:sNx,1:sNy,Nr)
      _RL tmp3dmskS(1:sNx,1:sNy,Nr)
      _RL df(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
#endif

c- note defined with overlap here, not needed, but more efficient
      _RL trVolW(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL trVolS(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL trHeatW(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL trHeatS(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL trSaltW(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL trSaltS(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)

#ifdef ATMOSPHERIC_LOADING
      _RL sIceLoadFac
#endif
#ifdef ALLOW_PSBAR_STERIC
      _RL RHOInSituLoc(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL VOLsumTile(nSx,nSy),RHOsumTile(nSx,nSy)
#endif

c need to include halos for find_rho_2d
      iMin = 1-OLx
      iMax = sNx+OLx
      jMin = 1-OLy
      jMax = sNy+OLy

#ifdef ALLOW_PSBAR_STERIC

      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
          do k = 1,nr
            CALL FIND_RHO_2D(
     I                iMin, iMax, jMin, jMax, k,
     I                theta(1-OLx,1-OLy,k,bi,bj),
     I                salt (1-OLx,1-OLy,k,bi,bj),
     O                RHOInSituLoc(1-OLx,1-OLy,k,bi,bj),
     I                k, bi, bj, myThid )
          enddo
        enddo
      enddo

      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
          RHOsumTile(bi,bj)=0. _d 0
          VOLsumTile(bi,bj)=0. _d 0
          VOLsumGlob=0. _d 0
          RHOsumGlob=0. _d 0
          do k = 1,nr
            do j = 1,sNy
              do i =  1,sNx
                RHOsumTile(bi,bj)=RHOsumTile(bi,bj)+
     &            (rhoConst+RHOInSituLoc(i,j,k,bi,bj))*
     &            hFacC(i,j,k,bi,bj)*drF(k)*rA(i,j,bi,bj)
                VOLsumTile(bi,bj)=VOLsumTile(bi,bj)+
     &            hFacC(i,j,k,bi,bj)*drF(k)*rA(i,j,bi,bj)
              enddo
            enddo
          enddo
        enddo
      enddo
      CALL GLOBAL_SUM_TILE_RL( VOLsumTile, VOLsumGlob, myThid )
      CALL GLOBAL_SUM_TILE_RL( RHOsumTile, RHOsumGlob, myThid )
      RHOsumGlob=RHOsumGlob/VOLsumGlob

      if (RHOsumGlob_0.GT.0. _d 0) then
        sterGloH=VOLsumGlob_0/globalArea
     &        *(1. _d 0 - RHOsumGlob/RHOsumGlob_0)
      else
        sterGloH=0. _d 0
      endif

c     WRITE(msgBuf,'(A,1PE21.14)') ' sterGloH= ', sterGloH
c        CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
c    &                       SQUEEZE_RIGHT, myThid )

#endif

#ifdef ATMOSPHERIC_LOADING
      sIceLoadFac=zeroRL
      IF ( useRealFreshWaterFlux ) sIceLoadFac=recip_rhoConst
#endif

      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
            do j = jmin,jmax
              do i =  imin,imax
                m_eta(i,j,bi,bj)=
     &                 etan(i,j,bi,bj)
#ifdef ATMOSPHERIC_LOADING
     &                +sIceLoad(i,j,bi,bj)*sIceLoadFac
#endif
#ifdef ALLOW_PSBAR_STERIC
     &                +sterGloH
#endif
              enddo
            enddo
        enddo
      enddo

      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
          do k = 1,nr
            do j = 1,sNy
              do i =  1,sNx
                m_UE(i,j,k,bi,bj)=0. _d 0
                m_VN(i,j,k,bi,bj)=0. _d 0
              enddo
            enddo
          enddo
        enddo
      enddo

      CALL ROTATE_UV2EN_RL(
     U          uVel, vVel, m_UE, m_VN,
     I          .TRUE., .TRUE., .FALSE., Nr, mythid )


c--   trVol : volume flux    --- [m^3/sec] (order of 10^6 = 1 Sv)
c--   trHeat: heat transport --- [Watt] (order of 1.E15 = PW)
c--   trSalt: salt transport --- [kg/sec] (order 1.E9 equiv. 1 Sv in vol.)
c--       convert from [ppt*m^3/sec] via rhoConst/1000.
c--       ( 1ppt = 1000*[mass(salt)]/[mass(seawater)] )

c-- init
      call ecco_zero(trVol,Nr,zeroRL,myThid)
      call ecco_zero(trHeat,Nr,zeroRL,myThid)
      call ecco_zero(trSalt,Nr,zeroRL,myThid)

#ifdef ALLOW_GENCOST_CONTRIBUTION

      do kgen=1,NGENCOST

      itr = gencost_itracer(kgen)

      call ecco_zero(gencost_storefld(1-OLx,1-OLy,1,1,kgen),
     &     1,zeroRL,myThid)

      do bj=myByLo(myThid),myByHi(myThid)
       do bi=myBxLo(myThid),myBxHi(myThid)
         areavolTile(bi,bj)=0. _d 0
       enddo
      enddo
      areavolGlob=0. _d 0

      do bj=myByLo(myThid),myByHi(myThid)
       do bi=myBxLo(myThid),myBxHi(myThid)
        do j = 1,sNy
         do i =  1,sNx
c---------
          do k = 1,nr
            tmpvol=hFacC(i,j,k,bi,bj)*drF(k)*rA(i,j,bi,bj)
c
            tmpmsk=0. _d 0
            if (.NOT.gencost_msk_is3d(kgen)) then
              tmpmsk=gencost_mskCsurf(i,j,bi,bj,kgen)*
     &               gencost_mskVertical(k,kgen)
#ifdef ALLOW_GENCOST3D
            else
              kgen3d=gencost_msk_pointer3d(kgen)
              tmpmsk=gencost_mskC(i,j,k,bi,bj,kgen3d)
#endif
            endif
c
            tmpfld=0. _d 0
            tmpmsk2=0. _d 0
            if (gencost_barfile(kgen)(1:15).EQ.'m_boxmean_theta') then
              tmpfld=theta(i,j,k,bi,bj)
              if (tmpmsk.NE.0. _d 0) tmpmsk2=1. _d 0
            elseif (gencost_barfile(kgen)(1:14).EQ.'m_boxmean_salt')
     &        then
              tmpfld=salt(i,j,k,bi,bj)
              if (tmpmsk.NE.0. _d 0) tmpmsk2=1. _d 0
#ifdef ALLOW_PTRACERS
            elseif (gencost_barfile(kgen)(1:17).EQ.'m_boxmean_ptracer')
     &        then
              tmpfld=pTracer(i,j,k,bi,bj,itr)
              if (tmpmsk.NE.0. _d 0) tmpmsk2=1. _d 0
#endif
            endif
c
            gencost_storefld(i,j,bi,bj,kgen) =
     &          gencost_storefld(i,j,bi,bj,kgen)
     &          +tmpmsk*tmpfld*tmpvol
            areavolTile(bi,bj)=areavolTile(bi,bj)
     &          +tmpmsk2*eccoVol_0(i,j,k,bi,bj)
c
          enddo
c---------
          tmpmsk=maskC(i,j,1,bi,bj)*gencost_mskCsurf(i,j,bi,bj,kgen)
          tmpfld=0. _d 0
          tmpmsk2=0. _d 0
          if (gencost_barfile(kgen)(1:13).EQ.'m_boxmean_eta') then
            tmpfld=m_eta(i,j,bi,bj)
            if (tmpmsk.NE.0. _d 0) tmpmsk2=1. _d 0
          endif
c
          gencost_storefld(i,j,bi,bj,kgen) =
     &        gencost_storefld(i,j,bi,bj,kgen)
     &        +tmpmsk*tmpfld*rA(i,j,bi,bj)
          areavolTile(bi,bj)=areavolTile(bi,bj)
     &        +tmpmsk2*rA(i,j,bi,bj)
c---------
          do k = 1,nr
c
            tmpmskW=0. _d 0
            tmpmskS=0. _d 0
            if (.NOT.gencost_msk_is3d(kgen)) then
              tmpmskW=gencost_mskWsurf(i,j,bi,bj,kgen)
     &          *gencost_mskVertical(k,kgen)
              tmpmskS=gencost_mskSsurf(i,j,bi,bj,kgen)
     &          *gencost_mskVertical(k,kgen)
#ifdef ALLOW_GENCOST3D
            else
              kgen3d=gencost_msk_pointer3d(kgen)
              tmpmskW=gencost_mskW(i,j,k,bi,bj,kgen3d)
              tmpmskS=gencost_mskS(i,j,k,bi,bj,kgen3d)
#endif
            endif

            if (
     &        (gencost_barfile(kgen)(1:15).eq.'m_horflux_theta').or.
     &        (gencost_barfile(kgen)(1:14).eq.'m_horflux_salt')
     &      ) then
c---------- Store local copy of section masks to add diffusive fluxes further down
              tmp3dmskW(i,j,k)=tmpmskW
              tmp3dmskS(i,j,k)=tmpmskS
            endif

            tmpmskW=tmpmskW*hFacW(i,j,k,bi,bj)*dyG(i,j,bi,bj)*drF(k)
            tmpmskS=tmpmskS*hFacS(i,j,k,bi,bj)*dxG(i,j,bi,bj)*drF(k)
c
            if (gencost_barfile(kgen)(1:13).EQ.'m_horflux_vol') then
              gencost_storefld(i,j,bi,bj,kgen) =
     &          gencost_storefld(i,j,bi,bj,kgen)
     &          +uVel(i,j,k,bi,bj)*tmpmskW
     &          +vVel(i,j,k,bi,bj)*tmpmskS
c
            elseif (gencost_barfile(kgen)(1:15).eq.'m_horflux_theta')
     & then
c------------ advective flux; exact if using enum_centered_2nd, cf. cad_c2_adv_x/y
              gencost_storefld(i,j,bi,bj,kgen) =
     &          gencost_storefld(i,j,bi,bj,kgen)
     &          +uVel(i,j,k,bi,bj)
     &            *(theta(i,j,k,bi,bj)+theta(i-1,j,k,bi,bj))*halfRL
     &            *HeatCapacity_Cp*rhoConst
     &            *tmpmskW
     &          +vVel(i,j,k,bi,bj)
     &            *(theta(i,j,k,bi,bj)+theta(i,j-1,k,bi,bj))*halfRL
     &            *HeatCapacity_Cp*rhoConst
     &            *tmpmskS
c
            elseif (gencost_barfile(kgen)(1:14).eq.'m_horflux_salt')
     & then
c------------ advective flux; exact if using enum_centered_2nd, cf. cad_c2_adv_x/y
              gencost_storefld(i,j,bi,bj,kgen) =
     &          gencost_storefld(i,j,bi,bj,kgen)
     &          +uVel(i,j,k,bi,bj)
     &            *(salt(i,j,k,bi,bj)+salt(i-1,j,k,bi,bj))*halfRL
     &            *rhoConst/1000.
     &            *tmpmskW
     &          +vVel(i,j,k,bi,bj)
     &            *(salt(i,j,k,bi,bj)+salt(i,j-1,k,bi,bj))*halfRL
     &            *rhoConst/1000.
     &            *tmpmskS
c------------
            ! Only compute tr[Vol,Heat,Salt] if necessary, use
            ! gencost_mask[W/S] rather than old msktrVol  
            elseif ( gencost_barfile(kgen)(1:7).eq.'m_trVol' .or. 
     &               gencost_barfile(kgen)(1:8).eq.'m_trHeat'.or.
     &               gencost_barfile(kgen)(1:8).eq.'m_trSalt'    ) then

                trVolW(i,j,k) =
     &                 uVel(i,j,k,bi,bj)*tmpmskW
     &                *maskInW(i,j,bi,bj)
                trVolS(i,j,k) =
     &                 vVel(i,j,k,bi,bj)*tmpmskS
     &                *maskInS(i,j,bi,bj)

                trHeatW(i,j,k) = trVolW(i,j,k)
     &                *(theta(i,j,k,bi,bj)+theta(i-1,j,k,bi,bj))*halfRL
     &                *HeatCapacity_Cp*rhoConst
                trHeatS(i,j,k) = trVolS(i,j,k)
     &                *(theta(i,j,k,bi,bj)+theta(i,j-1,k,bi,bj))*halfRL
     &                *HeatCapacity_Cp*rhoConst

                trSaltW(i,j,k) = trVolW(i,j,k)
     &                *(salt(i,j,k,bi,bj)+salt(i-1,j,k,bi,bj))*halfRL
     &                *rhoConst/1000.
                trSaltS(i,j,k) = trVolS(i,j,k)
     &                *(salt(i,j,k,bi,bj)+salt(i,j-1,k,bi,bj))*halfRL
     &                *rhoConst/1000.
c now summing
                trVol(i,j,k,bi,bj)=trVolW(i,j,k)+trVolS(i,j,k)
                trHeat(i,j,k,bi,bj)=trHeatW(i,j,k)+trHeatS(i,j,k)
                trSalt(i,j,k,bi,bj)=trSaltW(i,j,k)+trSaltS(i,j,k)

            endif
          enddo !k
c---------
         enddo !i
        enddo !j


c------ Add diffusive and GM/Redi fluxes to horflux_[theta/salt]
        if (
     &      (gencost_barfile(kgen)(1:15).eq.'m_horflux_theta').or.
     &      (gencost_barfile(kgen)(1:14).eq.'m_horflux_salt')
     &    ) then
c-------- Make local copy of tracer array (for this tile)
          if (gencost_barfile(kgen)(1:15).eq.'m_horflux_theta') then
            do k=1,nr
              do j=1-OLy,sNy+OLy
                do i=1-OLx,sNx+OLx
                  localT(i,j,k)=theta(i,j,k,bi,bj)
                enddo
              enddo
            enddo
            trIdentity = 1
            diffKh = diffKhT

          elseif (gencost_barfile(kgen)(1:14).eq.'m_horflux_salt') then
            do k=1,nr
              do j=1-OLy,sNy+OLy
                do i=1-OLx,sNx+OLx
                  localT(i,j,k)=salt(i,j,k,bi,bj)
                enddo
              enddo
            enddo
            trIdentity = 2
            diffKh = diffKhS
          endif
c----------
          do k=1,nr
            do j=1-OLy,sNy+OLy
              do i=1-OLx,sNx+OLx
                tmp2dmskW(i,j)=hFacW(i,j,k,bi,bj)*dyG(i,j,bi,bj)*drF(k)
                tmp2dmskS(i,j)=hFacS(i,j,k,bi,bj)*dxG(i,j,bi,bj)*drF(k)
              enddo
            enddo
c---------- Diffusive flux in X; assume deepFacC=1; ignore biharmonic
            if (diffKh.ne.0.) then
              call gad_diff_x(
     I          bi,bj,k,tmp2dmskW,diffKh,localT,
     O          df,
     I          myThid)
            else
              do j=1-OLy,sNy+OLy
                do i=1-OLx,sNx+OLx
                  df(i,j) = 0. _d 0
                enddo
              enddo
            endif
c---------- Add GM/Redi flux in X to df
#ifdef ALLOW_GMREDI
            if ( useGMredi ) then
              iMin = 1-OLx+1
              iMax = sNx+OLx
              jMin = 1-OLy
              jMax = sNy+Oly

              call gmredi_xtransport(
     I          trIdentity, bi, bj, k, iMin, iMax, jMin, jMax,
     I          tmp2dmskW,localT,
     U          df,
     I          myThid )
            endif
#endif
c---------- Add fluxes in X to gencost; assume rhoFacC = 1
            if (gencost_barfile(kgen)(1:15).eq.'m_horflux_theta')
     &      then
              do j=1,sNy
                do i=1,sNx
                  tmpmskW = tmp3dmskW(i,j,k)

                  gencost_storefld(i,j,bi,bj,kgen) =
     &              gencost_storefld(i,j,bi,bj,kgen)
     &              +df(i,j)*HeatCapacity_Cp*rhoConst
     &                *tmpmskW
                enddo
              enddo
c
            elseif (gencost_barfile(kgen)(1:14).eq.'m_horflux_salt')
     &      then
              do j=1,sNy
                do i=1,sNx
                  tmpmskW = tmp3dmskW(i,j,k)

                  gencost_storefld(i,j,bi,bj,kgen) =
     &              gencost_storefld(i,j,bi,bj,kgen)
     &              +df(i,j)*rhoConst/1000.
     &                *tmpmskW
                enddo
              enddo
            endif
c
c---------- Diffusive flux in Y; assume deepFacC=1; ignore biharmonic
            if (diffKh.ne.0.) then
              call gad_diff_y(
     I          bi,bj,k,tmp2dmskS,diffKh,localT,
     O          df,
     I          myThid)
            else
              do j=1-OLy,sNy+OLy
                do i=1-OLx,sNx+OLx
                  df(i,j) = 0. _d 0
                enddo
              enddo
            endif
c
c---------- Add GM/Redi flux in Y to df
#ifdef ALLOW_GMREDI
            if ( useGMredi ) then
              iMin = 1-OLx
              iMax = sNx+OLx
              jMin = 1-OLy+1
              jMax = sNy+Oly
              call gmredi_ytransport(
     I          trIdentity, bi, bj, k, iMin, iMax, jMin, jMax,
     I          tmp2dmskS,localT,
     U          df,
     I          myThid )
            endif
#endif
c
c---------- Add fluxes in Y to gencost; assume rhoFacC = 1
            if (gencost_barfile(kgen)(1:15).eq.'m_horflux_theta')
     &      then
              do j=1,sNy
                do i=1,sNx
                  tmpmskS = tmp3dmskS(i,j,k)

                  gencost_storefld(i,j,bi,bj,kgen) =
     &              gencost_storefld(i,j,bi,bj,kgen)
     &              +df(i,j)*HeatCapacity_Cp*rhoConst
     &                *tmpmskS
                enddo
              enddo
c
            elseif (gencost_barfile(kgen)(1:14).eq.'m_horflux_salt')
     &      then
              do j=1,sNy
                do i=1,sNx
                  tmpmskS = tmp3dmskS(i,j,k)

                  gencost_storefld(i,j,bi,bj,kgen) =
     &              gencost_storefld(i,j,bi,bj,kgen)
     &              +df(i,j)*rhoConst/1000.
     &                *tmpmskS
                enddo
              enddo
            endif
c
          enddo !k
        endif !gencost_barfile = m_horflux_[theta/salt]
c-----------------------------------------------------------------
       enddo
      enddo

      if (gencost_barfile(kgen)(1:9).EQ.'m_boxmean') then
        CALL GLOBAL_SUM_TILE_RL( areavolTile, areavolGlob, myThid )
        CALL ecco_div( gencost_storefld(1-OLx,1-OLy,1,1,kgen),
     &                 1, areavolGlob, myThid )
      endif

      enddo

#endif /* ALLOW_GENCOST_CONTRIBUTION */


      return
      end
