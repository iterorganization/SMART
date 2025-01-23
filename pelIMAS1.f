C===============================
      Subroutine pelIMAS1
     >(YAM,YVP,YVOL,YCOS0,YEFF,YDL,
     >  YDABL,YDDEP,YPELSRS,yswitch,
     >  ne,ni,Te,Ti,NHYDR,NDEUT,NTRIT,FP,
     >  ametr,shif,vr,mu,
     >  HRO,ROC,BTOR,RTOR,NA1,NRD)
! no common blocks, all inputs are explicit Polevoi 22-11-2022
!     Subroutine pelIMAS1 version 05-OCT-2023
!     >	(YAM,YVP,YVOL,YCOS0,YEFF,YDL,YDABL,YDDEP,YPELSRS,yswitch)
!       Subroutine pelnew(YAM,YVP,YVOL,YCOS0,YEFF,YDL,YDABL,YDDEP,yswitch)
!       Subroutine pelinj(YAM,YVP,YRP,YCOS0,YEFF,YSHAPE)
!       Subroutine psmar2(YAM,YVP,YRP,YCOS0,YEFF,YSHAPE)
! Simplified Mass Ablation and Relocation Treatment Polevoy 22-09-2000
! corrected 17-01-2001 converted to double precision 6-DEC-2012
!
! reduced output of the pellet source due to pellet ablation Polevoi 6-OCT-15
c=========================================================================
C Pellet ablation model by B.Kuteev, NF, 35 (1995) 431
c Cloud size by P.B.Parks, Phys. of Plasmas, 7 (2000) 1968
c Mass relocation model by H.R.Strauss, 5 (1998) 2676
c==========================================================================
c pelnew: c*19-NOV-2013
c YRP is replaced by pellet volume YVOL
c   YSHAPE is replaced by factor for effective length along the trajectory:
c  YDL = dRtraj/dRmid-plane
c  YCOS0 fror relacement is back
c
c=========================================================================
c Sperica/Cylindrical/Cubic pellet
c Input:
c   YAM 1-3 pellet atomic mass H,D,T or HD,DT pellets only
c  YVP > 0 [km/s] pellet velocity
c  YRP > 0 [cm] Rp initial radius of a pellet
c  YCOS = (1 , -1) for High Field Side (LFS) injection
c  YEFF > 0 fuelling efficiency
c  YSHAPE = 1(Spherical d=2Rp)
c    2(Cylinrical d=2Rp h=2Rp)
c    3(Cubic h=2Rp)
c  YVOL [mm3] pellet volume
c  yswitch (< 0 = does not affect plasma parameters Ne,Ni,Te,Ti,etc)
c Output:
c  YDABL,YDDEP depth of ablation and deposition (RHO/ROC)
c Ablation rate is represented in the form:
c dN/dt[atoms/s] =
c coeff * Ne**an * Te**at *rp**ap * Mi**am,
c where  rp is the pellet radius in cm,
c  Mi is the mass of the pellet material in atomic units
c  Te is plasma electron temperature in eV
c  Ne is the plasma density in cm-3
c-----------------------------------------------------------------
C Warnings  F1,F2,F3 densities must be prescribed in ASTRA
c   model explicily if proper pellet is used
c Subroutine SMART should be called after NE,TE,TI transport
c
c equations
c=======================================================================
c The values of the coefficients:
c Parks' model:
c an=0.333 at=1.64 ap=1.333 am=-0.333
c Hydrogen: coeff = 1.12e16
c Krypton:  coeff = 2.1e15
c
c Kuteev model:
c Deuterium: an=0.453 at=1.72 ap=1.443 am=-0.283
C coeff=3.46e14 (2e14 New)
c------------------------------------------------------------------
c=================================================================
c Input parameters:
c  Np  is the relative mass of the pellet in g/cm+3
c  Vb  is the pellet velosity in km/s
c  rp0 is the pellet radius in cm
c  Mi  is the mass of the pellet material in atomic units
c  cf = coeff/10e15
c  Zp  is the charge number of the pellet material
c  Eion is energy required for full ionization
c  of a pellet atom 13ev +molecule dissosiation 2.2 eV
c=====================================================================
         use mathematical_constants, only: M_PI
         implicit none
! include 'for/parameter.inc'
! include 'for/status.inc'
! include 'for/const.inc'
         integer NA,NA1,NRD !05-12-2022
         double precision
     >     ne(*),ni(*),Te(*),Ti(*),NHYDR(*),NDEUT(*),NTRIT(*),FP(*),
     >     ametr(*),shif(*),vr(*),mu(*),YPELSRS(*),
     >     HRO,ROC,BTOR,SHIFT,ABC,RTOR
         integer JS1,J0,JJ,JABS,JDEL,JBEG,JEND,JS,J,JMIN
         double precision YAM,YCN,YCE,YCI,YSTNE,YSTNI,YSTNE1,YSTNI1,YTST
         double precision YSDNE,YSDNI,YX,YDX,YF,YDF,YNE,YNI,YDNE,YDNI,YDL
         double precision YR,YDA,YDV,YLC,YFI,YPSI,YBETB,YR1,YR2,YF1,YF2
         double precision YRP1,YRP2,YA1,YA2,YA3,Y16,YPP,YP23,YCOS
         double precision YRP,YVP,YKC1,YEFF,YCOS0,YNp,YSIG,YHRO
         double precision Zp,an,Eion,at,ap,am,coeff,ycoef,yshape,ALFA
         double precision YDABL,YDDEP,YVOL,yswitch
c>26-APR-2023 M.H
         double precision, parameter :: eps = 1.0d-5
         integer n
         integer, parameter :: nmax = 1000
c<26-APR-2023 M.H
         double precision, allocatable ::
     >   YTE(:),YTI(:),YNTE(:),YNTI(:),YXJ(:),YX12(:),
     >   DNI(:),DNE(:),YKCR(:),YSFT(:),JJFP(:)

c======================================================================
C write(*,*) 'YAM,YVP,YRP,YCOS0,YEFF,YSHAPE',
C     . YAM,YVP,YRP,YCOS0,YEFF,YSHAPE
!====== 22-11-2022
         allocate (YTE(NRD))
         allocate (YTI(NRD))
         allocate (YNTE(NRD))
         allocate (YNTI(NRD))
         allocate (YXJ(NRD))
         allocate (YX12(NRD))
         allocate (DNI(NRD))
         allocate (DNE(NRD))
         allocate (YKCR(NRD))
         allocate (YSFT(NRD))
         allocate (JJFP(NRD))
         NA=NA1-1
         SHIFT=shif(NA1)
         ABC=AMETR(NA1)
c
c======================================================================
C	write(*,*) 'YAM,YVP,YRP,YCOS0,YEFF,YSHAPE',
C     .	YAM,YVP,YRP,YCOS0,YEFF,YSHAPE
c*19-NOV-2013 vvvvvvvvvvvvvvvvvv
         yshape=1.d0
         yrp=1.d-1*(.75d0/M_PI*yvol)**.333333 !YRP [cm]
         JMIN=NA1

c*19-NOV-2013 ^^^^^^^^^^^^^^^^^^^
c*6-OCT-15
!   yswitch=-1.d0 !!!!! No changes of plasma parameterts
!   yswitch=1.d0 !!!!! changes of plasma parameterts

         do j=1,na1
            YPELSRS(j) = NE(j)
         enddo
c*6-OCT-15
         if(YAM.lt.1.d0.or.YAM.gt.3.d0) then
            write(*,*) 'Warning from SMART: wrong input data'
            write(*,*) ' pellet mass must be H/D/T/HD/DT only'
            return
         endif
         if((YAM*YVP*YRP*(yshape-.999)*YEFF).le.0..or.dabs(YCOS0).gt.1.d0)
     .   then
            write(*,*) 'Warning from SMART: wrong input parameters'
            return
         endif
Constants for Kuteev's model
         Eion=15.2d-3   !(13.+2.2)*1.d-3
         Zp=1.d0
         an=0.453d0
         at=1.72d0
         ap=1.443d0
         am=-0.283d0
         YNp=2.98d0 !e22 molecules/cm3 of H/D/T
         coeff=2.d0 !3.46 !e14 coeff for dN/dt
cShape
         if(yshape.le.1.d0) then
ccc write(*,*) 'spherical pellet'
            ycoef=1.d0
         endif
         if(yshape.le.2.d0.and.yshape.gt.1.d0) then
ccc write(*,*) 'cylindrical pellet'
            ycoef=1.5d0
         endif
         if(yshape.gt.2.d0) then
ccc write(*,*) 'cubic pellet'
            ycoef=6.d0/M_PI
         endif

c=====================================================================
Constants
C effective opacity for inclined pellet
C*NEW vvvvvvv
c*19-NOV-2013 YDL=dSQRT(ELONG**2+YCOS0**2*(1.d0-ELONG**2))
c*19-NOV-2013 YCOSA=dabs(YCOS0)
c*19-NOV-2013 if(YCOSA.gt.1.d-6) then
         if(YCOS0.gt.0.d0) then
            YSIG=1.d0
         else
            YSIG=-1.d0
         endif
c*19-NOV-2013  YCOS=YSIG
         YCOS=YCOS0
c temporarily YCOS0 - inclination *to take into account in absorbtion
c YCOS - normal to the magn. surf. is parallel to the midplain
C*NEW ^^^^^^^^

         YA1=8.d0*M_PI*YNP  !4.d0*M_PI*2.d0*YNP
C 10**11=10**14*100[m->cm]/10**5[km/s->cm/s]/10**22
C /2 Simpson integration
         YA2=coeff*(3.d0-ap)*1.d-11*(1.d13)**an*1000.**at*YAM**am
         YA2=YA2/YVP/YA1/2.d0*YDL
C 1000=10**22/10**19 * shape correction (1 for sphere)
         YA3=1.d3/3.d0*YA1*YEFF *ycoef
         YP23=2.d0/3.d0/(3.d0-ap)
         YPP=3.d0/(3.d0-ap)
C for Kc = Rperp/Rp (Parks) cm->m
         Y16 =1.d0/6.d0
C*NEW vvvvvvv
c TSTAR =2 eV
         YTST =2.d0
         YKC1 =1.54d0*(10.d0/YAM)**Y16*(.01d0)**.667
     .   *dsqrt(13.33d0*YTST+31.6d0) /1.4d0
!     . *sqrt((4*5/3*YTST+13.6+2.2)/(1.-0.5))  /1.4
C*NEW /2.
C*NEW ^^^^^^^^

         do J=1,NA1
            DNI(J)=0.d0
            DNE(J)=0.d0
            YTE(J)=0.d0
            YTI(J)=0.d0
            YNTE(J)=0.d0
            YNTI(J)=0.d0

            YKCR(J)=0.d0
            YSFT(J)=0.d0
            YXJ(J)=(j-1)*HRO/ROC
            YX12(J)=(j-.5)*HRO/ROC

         enddo
         YXJ(NA1)=1.d0
         YX12(NA1)=1.d0

         JS=1
         J=NA1
         YR1=SHIFT-YSIG*JS*ABC
         YF1=NE(NA1)**an*TE(NA1)**at
         YRP1=YRP**(3.d0-ap)

c write(*,*) YA1,YA2,YA3
c>26-APR-2023 M.H
         n = 0
c<26-APR-2023 M.H
    1    J=J-JS
         YR2=YR1
         YF2=YF1
         YRP2=YRP1
         YR1=SHIF(J)-YSIG*JS*AMETR(J)
         YF1=NE(J)**an*TE(J)**at
         YRP1=YRP2 - YA2*dabs(YR2-YR1)*(YF1+YF2)

! YKCR(J)=YKC1*TE(J)**Y16/ (NE(J)*ALOG(2000.*TE(J)/7.5) )**.3333
         YKCR(J)=YKC1*TE(J)**Y16/(NE(J)*DLOG(2.667d2*TE(J)))**.3333
     .   *YRP2**YP23

c write(*,*) YF1,YF2,YR1,YR2,YRP1,YRP2,J
         if(YRP1.gt.0.d0) then
            DNI(J)= YA3*(YRP2**YPP-YRP1**YPP)+DNI(J)

         else
            DNI(J)= YA3*YRP2**YPP+DNI(J)

         endif
C*NEW
         JABS=J

         if(YRP1.le.0.) goto 2
! write(*,*) j
c>26-APR-2023 M.H
         n = n + 1
         if (n .gt. nmax) then
            write(*,*) 'Warning from SMART: Exited due to given profile'
            return
         endif
c<26-APR-2023 M.H
         if(J.eq.1) then
            J=0
            JS=-1
            goto 1
         endif
         if(J.lt.NA1) goto 1

    2    continue
c*NEW-1 vvvvvvvvvvvvv
         YKCR(NA1)=YKCR(NA)
         do J=NA1,JABS,-1
            JDEL=idint(2.0d0*YKCR(J)/HRO)+1
            JBEG=JDEL/2
            JBEG=J+JBEG
            YHRO=2.d0*YKCR(J)/dble(JDEL)
            if(JBEG.gt.NA1) JBEG=NA1
            JEND=JBEG-JDEL
            if(JEND.lt.1) JEND=1
            YDV=0.d0
            YDNE=0.d0
            do JJ=JBEG,JEND,-1
               if(JJ.eq.1) then
! YDV=VOLUM(1)
                  YDV=VR(1)*HRO
               else
                  YDV=VR(JJ)*YHRO+YDV
CVOLUM(J)-VOLUM(J-1)
               endif
               YDNE=YDNE+DNI(JJ)
            enddo
            YDNI = YDNE/YDV
            YNE =NE(J)+YDNI*Zp
            YNI =NI(J)+YDNI
C*NEW vvvvvvv
Cshift by Strauss dPSI=-q*B*bet*dN/dFI*cosT/(a*R*(n+dn)(1+da/dFI/a))
C dFI=Lc/R/q - toroidal angle length of a cloud
C Channel size by Parks: da=kc*rp
C    Lc = SQRT(kc*rp*R)
C
            YR=RTOR+SHIF(J)
C YLC=SQRT(2.*YKCR(J)*YR)
            YDA=YKCR(J)
            YLC=dSQRT(YKCR(J)*YR)
            YFI=YLC/(YR/MU(J))
c YFI=YLC/YR
            YBETB=4.d-3*(TE(J)*NE(J)+TI(J)*NI(J))
C 2*GP for ASTRA units
            YPSI=-YR/MU(J)*YBETB/(BTOR)*YDNE/YLC/
     .        AMETR(J)/YNE/(1.d0+YDA/AMETR(J)/YFI)*YCOS

c*YSIG
C*NEW ^^^^^^^^^^^^^
            YSFT(J)=YPSI

         enddo

c*NEW-1 ^^^^^^^^^^^^^^
         do J=1,NA1 !============================================
            if(J.eq.1) then
               YDV=VR(2)*HRO
!VOLUM(1)
            else
               YDV=VR(J)*HRO
            endif
            YDNI = DNI(J)/YDV
            YNE =NE(J)+YDNI*Zp
            YNI =NI(J)+YDNI

            YTE(J)=(TE(J)-0.667d0*YDNI/NE(J)*Eion)*NE(J)/YNE
            YTI(J)=TI(J)*NI(J)/YNI

         enddo
ccc
c index for density shift

         YDX=(1.d0/NRD)
         J=2
         YDF=FP(NA1)-FP(1)
         YF=(FP(J)-FP(1))/YDF

         do JJ=1,NRD
            YX=YDX*JJ
            if(YX.gt.YF) then
               J=J+1
               J= min(J,NA1)
               YF=(FP(J)-FP(1))/YDF
            endif
            JJFP(JJ)=J-1
         enddo
         JJFP(NRD)=NA1
C density, energy shift
         do j=NA1,1,-1

            if(DNI(J).ne.0.d0) then
               JJ=idint(dble(NRD)*(FP(J)+YSFT(J)-FP(1))/YDF) + 1
               if(JJ.lt.0) JJ=-JJ
               if(JJ.lt.NRD) then

                  J0=idint(JJFP(JJ))
                  JDEL=idint(2.0d0*YKCR(J)/HRO)+1
                  Jdel=min(JDEL,1)
                  JBEG=JDEL/2
                  JBEG=J0-JBEG
c*19-NOV-2013 vvvvvvvvvvvvvvvvvv

                  JMIN=min(jbeg,jmin)

c*19-NOV-2013 ^^^^^^^^^^^^^^^^^^^
c>03-MAY-2023 M.H
                  if(JMIN .le. 0) then
                     write(*,*)
     &                  'Warning from SMART: Exited due to n/E shift'
                     return
                  endif
c<03-MAY-2023 M.H
                  do JS1=1,JDEL
                     JS=JBEG+JS1-1
                     if(JS.lt.0)JS=-JS
                     if(JS.eq.0) JS=1
                     if(JS.le.NA1) then
                        DNE(JS)=DNI(J)/JDEL+DNE(JS)
                        YNTE(JS)=YTE(J)*DNI(J)/JDEL+YNTE(JS)
                        YNTI(JS)=YTI(J)*DNI(J)/JDEL+YNTI(JS)
                     endif
                  enddo
               endif
            endif
         enddo
C smoothing with energy/particle conservation
         ALFA = 0.001d0
         call SMTH(ALFA,NA1,DNE,YXJ,NA1,DNI,YX12,NRD)
         call SMTH(ALFA,NA1,YNTE,YXJ,NA1,YKCR,YX12,NRD)
         call SMTH(ALFA,NA1,YNTI,YXJ,NA1,YSFT,YX12,NRD)
         YSDNE=0.d0
         YSDNI=0.d0
         YSTNE=0.d0
         YSTNI=0.d0
         YSTNE1=0.d0
         YSTNI1=0.d0
         do J=1,NA1
            YSDNE=YSDNE+DNE(J)
            YSDNI=YSDNE+DNI(J)
            YSTNE=YSTNE+YNTE(J)
            YSTNI=YSTNI+YNTI(J)
            YSTNE1=YSTNE1+YKCR(J)
            YSTNI1=YSTNI1+YSFT(J)
            YNTE(J)=YKCR(J)
            YNTI(J)=YSFT(J)
         enddo

         YCN=YSDNE/YSDNI
         YCE=YSTNE/YSTNE1
         YCI=YSTNI/YSTNI1
         if(yci.le.0.d0.or.yce.le.0.d0.or.ycn.le.0.d0) then
            write(*,*) 'j,yni,yne,ydni'
            write(*,*) j,yci,yce,ycn
         endif
!===============================vvvv 26-03-2018
! extra smoothing for edge
!
!        ALFA = 0.1d0
!        do j=jmin,na1
!           YXJ(j-jmin+1)=YXJ(j)
!           YX12(j-jmin+1)=YX12(j)
! write(*,*) YXJ(j),YX12(j),jmin,j,na1
!        enddo
!        j=na1-jmin+1
!        call SMTH(ALFA,j,DNI,YXJ,j,DNE,YX12,NRD)
!        call SMTH(ALFA,j,DNE,YX12,j,DNI,YXj,NRD)

!===============================^^^^ 26-03-2018
         DO J=1,NA1

cc if(DNE(J).gt.0.d0) then
            if(J.eq.1) then
               YDV=VR(2)*HRO
            else
               YDV=VR(J)*HRO

            endif
            if(ydv.le.0.d0) write(*,*) 'j,dV',j,ydv
            YDNI = DNI(J)/YDV*YCN
            YNE =NE(J)+YDNI*Zp
            YNI =NI(J)+YDNI
            if(yni.le.0.d0.or.yne.le.0.d0) then
               write(*,*) 'j,yni,yne'
               write(*,*) j,yni,yne
            endif
c*19-NOV-2013 vvvvvvvvvvvvvvvvvv
            if(yswitch.ge.0.d0) then !=================== start switch

c*19-NOV-2013 ^^^^^^^^^^^^^^^^^^^

Ctemporary output
c Te
               TE(J) =(YTE(J)*NE(J)+YNTE(J)/YDV*YCE)/YNE
c Ti
               TI(J) =(YTI(J)*NI(J)+YNTI(J)/YDV*YCI)/YNI
c Ne
               NE(J) =NE(J)+YDNI*Zp

c>11-APR-2023 M.H
c       H
               if (abs(YAM-1.d0).lt.eps) then
                  NI(J)=NI(J)-NHYDR(J)
                  NHYDR(J)=NHYDR(J)+YDNI
                  NI(J)=NI(J)+NHYDR(J)
c       D
               else if (abs(YAM-2.d0).lt.eps) then
                  NI(J)=NI(J)-NDEUT(J)
                  NDEUT(J)=NDEUT(J)+YDNI
                  NI(J)=NI(J)+NDEUT(J)
c       T
               else if (abs(YAM-3.d0).lt.eps) then
                  NI(J)=NI(J)-NTRIT(J)
                  NTRIT(J)=NTRIT(J)+YDNI
                  NI(J)=NI(J)+NTRIT(J)
c       H/D
               else if (YAM.gt.1.d0.and.YAM.lt.2.d0) then
                  NI(J)=NI(J)-NHYDR(J)-NDEUT(J)
                  NHYDR(J)=NHYDR(J)+(2.d0-YAM)*YDNI
                  NDEUT(J)=NDEUT(J)+(YAM-1.d0)*YDNI
                  NI(J)=NI(J)+NHYDR(J)+NDEUT(J)
c       D/T
               else if (YAM.gt.2.d0.and.YAM.lt.3.d0) then
                  NI(J)=NI(J)-NTRIT(J)-NDEUT(J)
                  NTRIT(J)=NTRIT(J)+(YAM-2.d0)*YDNI
                  NDEUT(J)=NDEUT(J)+(3.d0-YAM)*YDNI
                  NI(J)=NI(J)+NTRIT(J)+NDEUT(J)
               endif
c*6-OCT-15 vvvv
            endif   !!!!!!!!!!!!!!!!!!!! end of yswitch

            YPELSRS(j) = YPELSRS(j)+YDNI*Zp

         enddo
c*19-NOV-2013 vvvvvvvvvvvvvvvvvv

         YDABL=1.d0 - dble(JABS)/dble(NA1)
         YDDEP=1.d0 - dble(JMIN)/dble(NA1)
c write(*,*) 'JABS,JMIN,rp= ',JABS,JMIN,rp
         YDABL=dmax1(YDABL,1.d-8)
         YDDEP=dmax1(YDDEP,1.d-8)

c*19-NOV-2013 ^^^^^^^^^^^^^^^^^^^


c*6-OCT-15 vvvv
         do j=1,na1
            YPELSRS(j) = YPELSRS(j)-NE(j)
         enddo
         deallocate (YTE)
         deallocate (YTI)
         deallocate (YNTE)
         deallocate (YNTI)
         deallocate (YXJ)
         deallocate (YX12)
         deallocate (DNI)
         deallocate (DNE)
         deallocate (YKCR)
         deallocate (YSFT)
         deallocate (JJFP)
         return
      end
!========================
