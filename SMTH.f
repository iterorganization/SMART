      subroutine MIXF19(RmDRs,YESCU,NA1,RHO,VNEWR,AMAIN,ZEF,
     > G33,G22,IPOL,
     > F1,F2,F3,F4,F5,F6,F7,F8,F9,TE,TI,NE,NI,MU,PFAST,PBLON,PBPER,
     > CU,CUTOR,EQPF,EQFF,NHYDR,NDEUT,NTRIT,NALF,NHE3,FP,
     > RTOR,BTOR,SHIFT,IMIX)
!==========================================Saw tooth mixing: 20-09-25
!	Polevoi A.R. 07-DEC-93, 18-MAY-09, 13-MAY-2011, 20-SEP-25
! Mix the arrays:
!	Ti,e, Ne, Ni, ion species F1-F9, fast ions over Rmix, if q=<1, rq'cr < rq'
!	with particle, energy and current conservation within the mixing area
! Input:
!	RmDRs = Rmix/Rs(q=1) the ratio of mixing radius to q=1 location (1.4d0 is recommended)
!	YESCU = 1.d0 (0.d0) switch to suppress current density CU, q and psi mixing (YESCU = 0.)
!	NA1 	! number of radial points (core_profiles)
!	RHO	! [m] dimentional radiall grid (Srt(ToroidalFlux/BTOR/3.1415))
!	VNEWR	! [m2] derivativ of the volume V'
!	AMAIN	! [a.u.] ion mass
!	Zef 	! [a.u.] effective ion charge
!	G33,IPOL ! Metric coefficients
!	BTOR(RTOR)! [T] |BTOR| module of vacuum magnetic field
!	RTOR	! [m] major radius
!	SHIFT	! [m] shift SHIFT= (RMAX+RMIN)/2-RTOR
! I/O:
! 	profiles before/after ST mix (if happened)
! Output:
!	IMIX (core profile grid) index of the surface within which ST is happened
!	IMIX =0 should be the last line of the outer time loop
!========================================================================
! ST stabilization by pressure/density gradient according to
! F.M.Levinton et al, Phys. Rev. Lett. Vol. 72, N 18, 2 May (1994)
!	beta_e=8\pi*press/B_0^2;
!	rq'cr=1.4(Ai/2Zef)^1/6*R(bet|n'|/n))^2/3(|p'|/p)^1/3
!	rq1'=rq': 	rq'cr > rq' = stable
!
! Flag 1 is written to the file dat/mix in the case of current mixing
! 	open(1,file='dat/mix')
!	write(1,*) IMIX
!	close(1)
!==================================================================
	implicit none
!	include	'for/parameter.inc'
!	include	'for/const.inc'
!	include	'for/status.inc'
!	include	'for/outcmn.inc'
	integer I,J,IRS,I1,IMIX,j1,NA,NA1
	double precision  RmDRs,YESCU,YNE,YTE,YTI,YV,YNI
     .         ,ysurf1,ytok,YCU,rqst,rqst1,YNA,YS,y
	double precision	YCB,YG,YTH2
	double precision	HH,YAJ,YCJ,GP,HRO,HROA,ROC,RTOR,BTOR,SHIFT
	double precision YF0,YF1,YF2,YF3,YF4,YF5,YF6,YF7,YF8,YF9
	double precision YNIZ1,YNIZ2,YNIZ3,YALF,YDEUT,YHYDR,YTRIT,YHE3
	double precision YI(NA1),YB(NA1),YA(NA1)
	double precision RHO(*),VNEWR(*),AMAIN(*),ZEF(*),G33(*),IPOL(*),
     > 	F1(*),F2(*),F3(*),F4(*),F5(*),F6(*),F7(*),F8(*),F9(*),
     > 	TE(*),TI(*),NE(*),NI(*),MU(*),PFAST(*),PBLON(*),PBPER(*),
     > 	CU(*),CUTOR(*),EQPF(*),EQFF(*),G22(*),FP(*),
     >  NHYDR(*),NDEUT(*),NTRIT(*),NALF(*),NHE3(*)
	NA=NA-1
	HRO	=RHO(3)-RHO(2)
	HROA	=RHO(NA1)-RHO(NA)
	ROC	=RHO(NA1)
	GP	=3.1414926d0

!	real YARR(*),yCUA(NRD)
	IMIX	=0
	IRS	=1
	I	=1
1	J	=NA1-I
     	I	=I+1
	IF(I.GE.NA1) RETURN
	IF(MU(J).GE.1.)	IRS=J
	IF(IRS.LE.1) GOTO 1
		I1	=IRS*RmDRs+1
	IF(I1.GE.NA) I1=NA
	IF(I1.LT.3) RETURN
c	write(*,*) 'Jrs=',Jrs
!	include 'fml/rqst'
C rqst:	ST stability limit
C F.M.Levinton et al, Phys. Rev. Lett. Vol. 72, N 18, 2 May (1994)
C		Polevoi (18-05-09)
C	beta_e=8\pi*press/B_0^2;
C	rq'cr=1.4(Ai/2Zef)^1/6*R(bet|n'|/n))^2/3(|p'|/p)^1/3
C
	if(j.lt.na1) then
	j1=j+1
	y= NI(j)*TI(j)+NE(j)*TE(j)+.5*(PBLON(j)+PBPER(j))+PFAST(j)
	rqst=
     .	NI(j1)*TI(j1)+NE(j1)*TE(j1)+.5*(PBLON(j1)+PBPER(j1))+PFAST(j1)
	rqst =1.4*(0.5*AMAIN(j)/Zef(j))**0.1667*(RTOR+SHIFT)/HRO*2.*
     .          (2.01E-3*(rqst+y)/(RTOR*BTOR/(RTOR+SHIFT))**2*
     .          abs( (NE(j1)-NE(j))/(NE(j1)+NE(j)) ))**.667*
     .  	abs((rqst-y)/(rqst+y))**0.333
	endif

!	include 'fml/rqst1'
C rqst1:	ST stability limit
C F.M.Levinton et al, Phys. Rev. Lett. Vol. 72, N 18, 2 May (1994)
C		Polevoi (18-05-09)
C	beta_e=8\pi*press/B_0^2;
C	rq'cr=1.4(Ai/2Zef)^1/6*R(bet|n'|/n))^2/3(|p'|/p)^1/3
C	rq1'=rq': 	rq'cr > rq' = stable
C
	if(j.lt.na1) then
	rqst1 = j*(1./MU(j+1)-1./MU(j))
	endif

c	write(*,*) 'rqst>rqst1=',rqst,rqst1
	if(rqst.gt.rqst1) return
		IMIX	=I1-1
		YNE	=0.d0
		YTE	=0.d0
		YTI	=0.d0
		YV	=0.d0
		YNI	=0.d0
		YNA	=0.d0
		YS  =0.d0
		YF0	=0.d0
		YF1	=0.d0
		YF2	=0.d0
		YF3	=0.d0
		YF4	=0.d0
		YF5	=0.d0
		YF6	=0.d0
		YF7	=0.d0
		YF8	=0.d0
		YF9	=0.d0
		YNIZ1	=0.d0
		YNIZ2	=0.d0
		YNIZ3	=0.d0
		YALF	=0.d0
		YDEUT	=0.d0
		YHYDR	=0.d0
		YTRIT	=0.d0
		YHE3	=0.d0
	IF(YESCU.NE.0.)	THEN
	   ysurf1=0.d0
	   ytok=0.d0
	       YCU=1.25/(GP*GP*RTOR)
	   do j=1,IMIX
	      y=RHO(j)*HRO/G33(J)/IPOL(J)**3
	      ysurf1=ysurf1 + y
	      ytok=ytok+CU(J)*y
       		YI(j)=HRO**3*(j-0.5)/(YCU*G33(J)*IPOL(J)**3)
	      enddo
C ==== average current density
	      ytok=ytok/ysurf1
	   endif

	DO J	=1,IMIX
            	YS    	=YS+MU(J)*J
		YV	=YV+VNEWR(J)
		YTE	=YTE+NE(J)*TE(J)*VNEWR(J)
		YNE	=YNE+NE(J)*VNEWR(J)
		YNI	=YNI+NI(J)*VNEWR(J)
		YTI	=YTI+NI(J)*TI(J)*VNEWR(J)
!		YF0	=YF0+F0(J)*VNEWR(J)
		YF1	=YF1+F1(J)*VNEWR(J)
		YF2	=YF2+F2(J)*VNEWR(J)
		YF3	=YF3+F3(J)*VNEWR(J)
		YF4	=YF4+F4(J)*VNEWR(J)
		YF5	=YF5+F5(J)*VNEWR(J)
		YF6	=YF6+F6(J)*VNEWR(J)
		YF7	=YF7+F7(J)*VNEWR(J)
		YF8	=YF8+F8(J)*VNEWR(J)
		YF9	=YF9+F9(J)*VNEWR(J)
!		YNIZ1	=YNIZ1+NIZ1(J)*VNEWR(J)
!		YNIZ2	=YNIZ2+NIZ2(J)*VNEWR(J)
!		YNIZ3	=YNIZ3+NIZ3(J)*VNEWR(J)
		YNIZ1	=YNIZ1+PFAST(J)*VNEWR(J)
		YNIZ2	=YNIZ2+PBLON(J)*VNEWR(J)
		YNIZ3	=YNIZ3+PBPER(J)*VNEWR(J)
		YALF	=YALF+NALF(J)*VNEWR(J)
		YDEUT	=YDEUT+NDEUT(J)*VNEWR(J)
		YHYDR	=YHYDR+NHYDR(J)*VNEWR(J)
		YTRIT	=YTRIT+NTRIT(J)*VNEWR(J)
		YHE3	=YHE3+NHE3(J)*VNEWR(J)

	enddo
            	YS    =(YS+1.E-6)/MU(I1)
c=27-MAY-09
	IF(YESCU.NE.0.d0) then
cc          CU(IMIX)= YTOK
!====================================nporOHKA==========for fp
		YB(1) = -YTOK*YI(1)/G22(1)
		YA(1) = 1.d0
	do j=2,IMIX
		j1=j-1
	y=G22(j1)*(YA(j1)-1.d0)-G22(j)
		YA(j)=-G22(j)/y
		YB(j)=(-YB(j1)*G22(j1)+YTOK*YI(j))/y
	enddo

      do J=IMIX-1,1,-1
	 j1=j+1
cc        CU(J)= YTOK
	fp(j) = YA(j)*fp(j1)+YB(j)
	enddo

	endif
       DO 	J	=1,IMIX
C...If you don't want to mix Ne,i, Te,i then comment the proper lines
		NE(J)	=YNE/YV
		NI(J)	=YNI/YV
		TE(J)	=YTE/YNE
		TI(J)	=YTI/YNI
!		F0(J)	=YF0/YV
		F1(J)	=YF1/YV
		F2(J)	=YF2/YV
		F3(J)	=YF3/YV
		F4(J)	=YF4/YV
		F5(J)	=YF5/YV
		F6(J)	=YF6/YV
		F7(J)	=YF7/YV
		F8(J)	=YF8/YV
		F9(J)	=YF9/YV
	!	NIZ1(J)	=YNIZ1/YV
	!	NIZ2(J)	=YNIZ2/YV
	!	NIZ3(J)	=YNIZ3/YV
		PFAST(J)=YNIZ1/YV
		PBLON(J)=YNIZ2/YV
		PBPER(J)=YNIZ3/YV
		NALF(J)	=YALF/YV
		NDEUT(J)=YDEUT/YV
		NTRIT(J)=YTRIT/YV
		NHE3(J)	=YHE3/YV
	enddo

call CUOFP !======================================================================!

	HH = HRO*HRO
	YAJ = 0.
C	write(*,*)"CU"
C	write(*,100)(CU(j),j=NA1-5,NA1)
	do	J=1,IMIX-1
	   YCJ = YAJ
	   if (j .lt. NA)	then
	      YAJ = (FP(j+1)-FP(j))/HH
	      MU(j) = YAJ/j
	      YAJ = G22(j)*YAJ
	      CU(j) = (YAJ-YCJ)/HRO
	   else
	      YAJ = (FP(j+1)-FP(j))/HRO/HROA
	      MU(j) = YAJ/j
	      YAJ = G22(j)*YAJ
C30-11	      CU(j) = 2.*(YAJ-YCJ)/(HRO+HROA)
	      CU(j) = (YAJ-YCJ)/HRO
	   endif
	   CU(j) = CU(j)/(j-0.5)
 	enddo

	YCJ = 1.25/(GP*GP*RTOR)
	YAJ = 0.5/(GP*BTOR)
	do	J=1,IMIX-1
	   CU(j) = YCJ*CU(j)*G33(J)*IPOL(J)**3
           MU(J) = YAJ*MU(j)
 	enddo
C Preparing input for equilibrium solver:
	YCB = 1.6E-3*RTOR/(BTOR*HRO*HRO)
	do	J=1,IMIX-1
	   CUTOR(J) = RHO(J)/ROC
	   if (j.eq.NA) YCB = YCB*HRO/HROA

	   EQFF(J) = (NE(J+1)*TE(J+1)-NE(J)*TE(J))
     +		    +(NI(J+1)*TI(J+1)-NI(J)*TI(J))
     +			+(pfast(J+1)-pfast(J))
     +			+0.5*(PBLON(J+1)-PBLON(J)+PBPER(J+1)-PBPER(J))

!     +			+0.5*NB2EQL*(PBLON(J+1)-PBLON(J)+PBPER(J+1)-PBPER(J))

	   EQFF(J) = -YCB*EQFF(J)/(MU(J)*J)
	enddo
	do	J=1,IMIX-1

	   EQPF(j) = EQFF(j)
	   YTH2	= RHO(j)*G22(J)*(MU(J)/RTOR)**2
	   YG	= (1.+YTH2)*G33(J)
	   EQFF(J) = (CU(J)/IPOL(J)-EQPF(J))/YG
	   CUTOR(J) = (CU(J)/IPOL(J)+YTH2*EQPF(J))/(1.+YTH2)

	enddo
c============================================

	open(1,file='mix')
	write(1,*) IMIX
	close(1)

	RETURN
	END
C===================================================================
C======================================================================|
	subroutine	CUOFPy(NA1,RHO,RTOR,BTOR,FP,MU,CU,G22,G33,IPOL)
C----------------------------------------------------------------------|
C Compute CU(rho) and MU(rho) from FP(rho)
C----------------------------------------------------------------------|
C	Input:	HRO	- radial step (m)
C		HROA	- edge radial step (m)
C		NA1	- number of grid points
C		G22(1:NA)	- <g22/g>*............
C		G33(1:NA)	-
C		IPOL(1:NA)	-
C		CD(1:NA)	- external (+bootstrap) current
C		FP(1:NA1)	- poloidal flux
C	Output:	CU(1:NA1) - (1/rho)d{K*dF/d(rho)}/d(rho) current density
C		MU(1:NA1) - (1/rho)dF/d(rho)	    rotational transform
! made from CUOFP Polevoi 22-JUL-2025
C----------------------------------------------------------------------|
	implicit none
!	include	'for/parameter.inc'
!	include	'for/const.inc'
!	include	'for/status.inc'
	integer	j,NA,NA1
	double precision	HH,YAJ,YCJ,ARRNA1
	double precision FP(*),MU(*),CU(*),G22(*),G33(*),IPOL(*),RHO(*)
	double precision HRO,HROA,RTOR,BTOR,GP,ROC
      external ARRNA1
		NA=NA1-1
        HRO=RHO(4)-RHO(3)
        ROC=RHO(NA1)
        NA=NA1-1
        HROA=ROC-RHO(NA)
        GP=3.14159263359d0
        HROA=HRO

        HH = HRO*HRO
        YAJ = 0.
C	write(*,*)"CU"
C	write(*,100)(CU(j),j=NA1-5,NA1)
	do	1	J=1,NA
	   YCJ = YAJ
	if (j .lt. NA)	then
	      YAJ = (FP(j+1)-FP(j))/HH
	      MU(j) = YAJ/j
	      YAJ = G22(j)*YAJ
	      CU(j) = (YAJ-YCJ)/HRO
	else
	      YAJ = (FP(j+1)-FP(j))/HRO/HROA
	      MU(j) = YAJ/j
	      YAJ = G22(j)*YAJ
C30-11	      CU(j) = 2.*(YAJ-YCJ)/(HRO+HROA)
	      CU(j) = (YAJ-YCJ)/HRO
	endif

		CU(j) = CU(j)/(j-0.5)
 1	continue
      MU(1) =MU(2)
C Suppress the next line to avoid CU(NA1) =/= 0
		CU(NA1) = ARRNA1(CU(NA-2),HROA/HRO)
		HH = HROA/(HRO+HROA)		! Attributes MU(NA1)
		MU(NA1) = ARRNA1(MU(NA-2),HH)	! to the point ROC
C	MU(NA1) = MU(NA)*NA/(NA-0.5+HROA/HRO)
		YCJ = 1.25/(GP*GP*RTOR)
		YAJ = 0.5/(GP*BTOR)
	do	2	J=1,NA1
		CU(j) = YCJ*CU(j)*G33(J)*IPOL(J)**3
		MU(J) = YAJ*MU(j)
 2	continue
C	CU(NA)=CUBS(NA)+CD(NA)+CC(NA)*ULON(NA)/(RTOR*GP2)
C	CU(NA1)=CUBS(NA1)+CD(NA1)+CC(NA1)*ULON(NA1)/(RTOR*GP2)
C	write(*,100)Time
C	write(*,100)(CU(j),j=NA1-5,NA1)
C	write(*,100)(MU(j),j=NA1-5,NA1)
 100	format(3(2F10.5,2X))
	end
C======================================================================|
	double precision function ARRNA1(ARR,H)
C----------------------------------------------------------------------|
C Compute edge value of the array ARR(NA1) assuming that d^3(ARR)/dr^3=0
C----------------------------------------------------------------------|
	double precision ARR(3),H
C	ARRNA1 = ARR(1)+H*(2.*ARR(1)+ARR(-2)-3.*ARR(-1))
		ARRNA1 = ARR(3)*(1+(1.5+0.5*H)*H)+ARR(1)*H*0.5*(1+H)-
     >		ARR(2)*H*(2+H)
	end
C======================================================================|
      Subroutine CUBSy(NA1,RTOR,BTOR,IPL,
     > FP,MU,ZEF,TE,TI,NE,NI,AMAIN,ZMAIN,
     > BMINT,BMAXT,BDB0,BDB02,FOFB,SQEPS,RHO,
     > CUBS,CC) !output bootstrap [MA/m2] current conductivity [1/mkOM*m]
! Bootsrap current and current conductvity by Sauter
!                       Sauter, Angioni, Lin-Liu
!			Physics of Plasmas, 6, 2834 (1999)
!	(Polevoi 22-JUL-2025)
    	implicit none
	integer	j,NA,NA1
	double precision
     > YC,YD,YA,ZZ,ZDF,ZFT,ZFTE,ZFTE1,ZFTE2,ZFTE3,ZFTE4,
     > ZFTI,ZFTI1,ZFTI2,ZFTI3,ZFTI4,A0,A1,ALP,
     > BETPL,COULG,NUEE,NUES,NUI,NUIS,HCEE,HCEI,
     > DCSA,HCSA,XCSA,CCSP,CNSA
	double precision
     > FP(*),MU(*),CUBS(*),CC(*),ZEF(*),TE(*),TI(*),NE(*),NI(*),
     > BMINT(*),BMAXT(*),BDB0(*),BDB02(*),FOFB(*),
     > RHO(*),SQEPS(*),AMAIN(*),ZMAIN(*)
	double precision HRO,HROA,RTOR,BTOR,GP,FTLLMRy,HC,XC,DC,ROC,IPL
	external FTLLMRy
        HRO=RHO(3)-RHO(2)
        ROC=RHO(NA1)
        NA=NA1-1
        HROA=ROC-RHO(NA)
        GP=3.14159263359d0
!    call	markloc("Current equation"//char(0))
      YC=.4*GP*RTOR/BTOR
      YD=-0.8*GP*GP*RTOR
      YA=2./(HRO**2*YD)
!         write(*,*) 'cu',
!     > FTLLMRy(BTOR,BMAXT(1),BDB0(1),BDB02(1),FOFB(1))
      do 2224 J=1,NA1
	if(J.ne.NA1)	then
	   BETPL = 1.6E-4*GP*(NE(J+1)+NE(J))*(TE(J+1)+TE(J))
	   BETPL = BETPL*(RTOR/(BTOR*J*HRO*ABS(MU(J)+.0001)))**2
	else
	   BETPL = 6.4E-4*GP*NE(J)*TE(J)
	   BETPL = BETPL*(RTOR/(BTOR*ROC*ABS(MU(J)+.0001)))**2
	endif
	ZZ=ZEF(J)
        ZFT=FTLLMRy(BTOR,BMAXT(J),BDB0(J),BDB02(J),FOFB(J))
!                 write(*,*) 'j',j, ZFT

	COULG=15.9d0-0.5d0*LOG(NE(J))+log(TE(J))
	NUEE=670.*COULG*NE(J)/(TE(J)*SQRT(TE(J)))
	NUES=NUEE*1.4*ZEF(J)*RTOR
     *	/ABS(MU(J)*(SQEPS(J)**3)*SQRT(TE(J))*1.875E7)
	ZDF = 1.+(1.-0.1*ZFT)*SQRT(NUES)
	ZDF = ZDF + 0.5*(1.-ZFT)*NUES/ZZ
	ZFT = ZFT/ZDF
!	               write(*,*) 'j', j, ZFT
       	DCSA = (1+1.4/(ZZ+1))*ZFT-1.9/(ZZ+1)*ZFT*ZFT
        DCSA = DCSA + (0.3*ZFT*ZFT+0.2*ZFT*ZFT*ZFT)*ZFT/(ZZ+1)
	DCSA=DCSA*BETPL*(1.+TI(J)/(ZZ*TE(J)))
!      DC(J)=DCSA
!      write(*,*) 'j, DCSA ', j, DCSA
        DC=DCSA
!	if(J.ne.NA1)	then
!	   BETPL = 1.6E-4*GP*(NE(J+1)+NE(J))*(TE(J+1)+TE(J))
!	   BETPL = BETPL*(RTOR/(BTOR*J*HRO*ABS(MU(J)+.0001)))**2
!	else
!	   BETPL = 6.4E-4*GP*NE(J)*TE(J)
!	   BETPL = BETPL*(RTOR/(BTOR*ROC*ABS(MU(J)+.0001)))**2
!	endif
	ZZ=ZEF(J)
        ZFT=FTLLMRy(BTOR,BMAXT(J),BDB0(J),BDB02(J),FOFB(J))
!	COULG=15.9d0-0.5d0*LOG(NE(J))+log(TE(J))
	NUEE=670.*COULG*NE(J)/(TE(J)*SQRT(TE(J)))
	NUES=NUEE*1.4*ZEF(J)*RTOR
     *	/ABS(MU(J)*(SQEPS(J)**3)*SQRT(TE(J))*1.875E7)
	ZDF = 1.+0.26*(1.-ZFT)*SQRT(NUES)
	ZDF = ZDF + 0.18*(1.-0.37*ZFT)*NUES/SQRT(ZZ)
        ZFTE = ZFT/ZDF
        ZFTE2=ZFTE*ZFTE
        ZFTE3=ZFTE*ZFTE2
	ZFTE4=ZFTE2*ZFTE2
        ZDF = 1.+(1.+0.6*ZFT)*SQRT(NUES)
	ZDF = ZDF + 0.85*(1.-0.37*ZFT)*NUES*(1+ZZ)
	ZFTI = ZFT/ZDF
        ZFTI2=ZFTI*ZFTI
        ZFTI3=ZFTI*ZFTI2
	ZFTI4=ZFTI2*ZFTI2
	HCEE = (0.05+0.62*ZZ)/ZZ/(1.+0.44*ZZ)*(ZFTE-ZFTE4)
	HCEE = HCEE +(ZFTE2-ZFTE4-1.2*(ZFTE3-ZFTE4))/(1+0.22*ZZ)
	HCEE = HCEE + 1.2/(1+0.5*ZZ)*ZFTE4
        HCEI = -(0.56+1.93*ZZ)/ZZ/(1.+0.44*ZZ)*(ZFTI-ZFTI4)
	HCEI = HCEI+4.95/(1+2.48*ZZ)*(ZFTI2-ZFTI4-0.55*(ZFTI3-ZFTI4))
	HCEI = HCEI-1.2/(1+0.5*ZZ)*ZFTI4
!	if(J.ne.NA1)	then
!	   BETPL = 1.6E-4*GP*(NE(J+1)+NE(J))*(TE(J+1)+TE(J))
!	   BETPL = BETPL*(RTOR/(BTOR*J*HRO*ABS(MU(J)+.0001)))**2
!	else
!	   BETPL = 6.4E-4*GP*NE(J)*TE(J)
!	   BETPL = BETPL*(RTOR/(BTOR*ROC*ABS(MU(J)+.0001)))**2
!	endif
	ZZ=ZEF(J)
        ZFT=FTLLMRy(BTOR,BMAXT(J),BDB0(J),BDB02(J),FOFB(J))
!	COULG=15.9d0-0.5d0*LOG(NE(J))+log(TE(J))
	NUEE=670.*COULG*NE(J)/(TE(J)*SQRT(TE(J)))
	NUES=NUEE*1.4*ZEF(J)*RTOR
     *	/ABS(MU(J)*(SQEPS(J)**3)*SQRT(TE(J))*1.875E7)
	ZDF = 1.+(1.-0.1*ZFT)*SQRT(NUES)
	ZDF = ZDF + 0.5*(1.-ZFT)*NUES/ZZ
	ZFT = ZFT/ZDF
       	DCSA = (1.+1.4/(ZZ+1))*ZFT-1.9/(ZZ+1)*ZFT*ZFT
        DCSA = DCSA + (0.3*ZFT*ZFT+0.2*ZFT*ZFT*ZFT)*ZFT/(ZZ+1.)
	DCSA=DCSA*BETPL*(1.+TI(J)/(ZZ*TE(J)))
        HCSA = BETPL*(HCEE+HCEI)+DCSA/(1.+TI(J)/(ZZ*TE(J)))
!      HC(J)=HCSA
!            write(*,*) 'j, HCSA ', j, HCSA
      HC=HCSA
!	if(J.ne.NA1)	then
!	   BETPL = 1.6E-4*GP*(NE(J+1)+NE(J))*(TE(J+1)+TE(J))
!	   BETPL = BETPL*(RTOR/(BTOR*J*HRO*ABS(MU(J)+.0001)))**2
!	else
!	   BETPL = 6.4E-4*GP*NE(J)*TE(J)
!	   BETPL = BETPL*(RTOR/(BTOR*ROC*ABS(MU(J)+.0001)))**2
!	endif
	ZZ=ZEF(J)
        ZFT=FTLLMRy(BTOR,BMAXT(J),BDB0(J),BDB02(J),FOFB(J))
!	COULG=15.9d0-0.5d0*LOG(NE(J))+log(TE(J))
	NUEE=670.*COULG*NE(J)/(TE(J)*SQRT(TE(J)))
	NUES=NUEE*1.4*ZEF(J)*RTOR
     *	/ABS(MU(J)*(SQEPS(J)**3)*SQRT(TE(J))*1.875E7)
	NUI=ZMAIN(J)**4*NI(J)*322./(TI(J)*SQRT(TI(J)*AMAIN(J)))
	NUIS=3.2E-6*NUI*RTOR/(ABS(MU(J)+.0001)*
     *		SQEPS(J)**3*SQRT(TI(J)/AMAIN(J)))
	ZDF = 1.+(1.-0.1*ZFT)*SQRT(NUES)
	ZDF = ZDF + 0.5*(1.-0.5*ZFT)*NUES/ZZ
        ZFTE = ZFT/ZDF
	XCSA = (1+1.4/(ZZ+1))*ZFTE-1.9/(ZZ+1.)*ZFTE*ZFTE
        XCSA = XCSA + (0.3*ZFTE*ZFTE+0.2*ZFTE*ZFTE*ZFTE)*ZFTE/(ZZ+1)
	A0 = -(1.17)*(1.-ZFT)
        A0 =A0/(1.-0.22*ZFT-0.19*ZFT*ZFT)
	ALP =(A0+0.25*(1-ZFT*ZFT)*SQRT(NUIS))/(1.+0.5*SQRT(NUIS))
        A1 = NUIS*NUIS*ZFT*ZFT*ZFT*ZFT*ZFT*ZFT
        ALP = (ALP + 0.315*A1)/(1.+0.15*A1)
!	if(J.ne.NA1)	then
!	   BETPL = 1.6E-4*GP*(NE(J+1)+NE(J))*(TE(J+1)+TE(J))
!	   BETPL = BETPL*(RTOR/(BTOR*J*HRO*ABS(MU(J)+.0001)))**2
!	else
!	   BETPL = 6.4E-4*GP*NE(J)*TE(J)
!	   BETPL = BETPL*(RTOR/(BTOR*ROC*ABS(MU(J)+.0001)))**2
!	endif
	ZZ=ZEF(J)
        ZFT=FTLLMRy(BTOR,BMAXT(J),BDB0(J),BDB02(J),FOFB(J))
!	COULG=15.9d0-0.5d0*LOG(NE(J))+log(TE(J))
	NUEE=670.*COULG*NE(J)/(TE(J)*SQRT(TE(J)))
	NUES=NUEE*1.4*ZEF(J)*RTOR
     *	/ABS(MU(J)*(SQEPS(J)**3)*SQRT(TE(J))*1.875E7)
	ZDF = 1.+(1.-0.1*ZFT)*SQRT(NUES)
	ZDF = ZDF + 0.5*(1.-ZFT)*NUES/ZZ
	ZFT = ZFT/ZDF
       	DCSA = (1.+1.4/(ZZ+1.))*ZFT-1.9/(ZZ+1.)*ZFT*ZFT
        DCSA = DCSA + (0.3*ZFT*ZFT+0.2*ZFT*ZFT*ZFT)*ZFT/(ZZ+1)
	DCSA=DCSA*BETPL*(1.+TI(J)/(ZZ*TE(J)))
        XCSA = BETPL*(XCSA*ALP)*TI(J)/ZZ/TE(J)
        XCSA = XCSA+DCSA/(1.+ZZ*TE(J)/TI(J))
        XC=XCSA
!      XC(J)=XCSA
!      CD(J)=CUBM(J)+CUECR(J)+CUICR(J)+CULH(J)+CUFW(J)+CAR8(J)
!     + +CAR18(J)+CAR28(J)
	ZZ=ZEF(J)
        ZFT=FTLLMRy(BTOR,BMAXT(J),BDB0(J),BDB02(J),FOFB(J))
!	COULG=15.9d0-0.5d0*LOG(NE(J))+log(TE(J))
	NUEE=670.*COULG*NE(J)/(TE(J)*SQRT(TE(J)))
	NUES=NUEE*1.4*ZEF(J)*RTOR
     *	/ABS(MU(J)*(SQEPS(J)**3)*SQRT(TE(J))*1.875E7)
        ZDF = 1.+(0.55-0.1*ZFT)*SQRT(NUES)
	ZDF = ZDF + 0.45*(1.-ZFT)*NUES/ZZ/SQRT(ZZ)
        ZFT = ZFT/ZDF
	CNSA=1.-(1.+0.36/ZZ)*ZFT
        CNSA=CNSA+0.59/ZZ*ZFT*ZFT-0.23/ZZ*ZFT*ZFT*ZFT
!	COULG=15.9d0-0.5d0*LOG(NE(J))+log(TE(J))
	CCSP = 1041.d0/COULG*TE(J)*SQRT(TE(J))*
     *		(1.13d0+ZEF(J))/(2.67d0+ZEF(J))/ZEF(J)
	CNSA=CCSP*CNSA
      CC(J)=CNSA
!      write(*,*) 'j, CC', j , CC(j)
      if (j .eq. NA1) goto 2222
      if (j .eq. NA) YA = 2./(YD*HROA**2)
      CUBS(J)=YA*(FP(J+1)-FP(J))*(0.
     > +HC*(TE(J+1)-TE(J))/(TE(J+1)+TE(J))
     > +XC*(TI(J+1)-TI(J))/(TI(J+1)+TI(J))
     > +DC*(NE(J+1)-NE(J))/(NE(J+1)+NE(J))
     > )
      if(j.eq.NA1) CUBS(J)= max(0.,(2*CUBS(NA)-CUBS(NA-1)))
      goto 2223
 2222    continue
 2223    continue
 2224    continue
         return
         end
!=============================
C======================================================================|
	subroutine	RHSEQy(NA1,RTOR,BTOR,RHO,
     > NE,NI,TE,TI,PBLON,PBPER,PFAST,MU,CU,
     > G22,G33,IPOL,AMETR,
     > CUTOR,EQFF,EQPF)
C-----------------------------------------------------------------------
C Input: RTOR,BTOR,NA,NA1,HRO,HROA,   NB2EQL,
C	 NE,NI,TE,TI,MU,CU,AMETR,RHO, PBLON,PBPER,G22,G33,IPOL
C Output:
C       EQPF
C       EQFF
C       CUTOR
C Both quantities EQPF (~p') and EQFF (~II') are given in [MA/m^2]
C	EQPF = -1.6E-3*(2*\pi*R_0)\prti{n_13*T_keV}{\psi[Vs=T*m^2]}
C	     = -1.E-6*/(2*\pi*R_0)\prti{p[J/m^3=Pascal]}{\psi[Vs]}
C	EQFF = -1.E-6*2*\pi/(R_0*\mu_0)*I*\prti{I}{\psi}
C	     = -5./R_0*I*\prti{I}{\psi}
C Local toroidal current density j[MA/m^2] is EQPF*r/R_0+EQFF*R_0/r, i.e.
C       j(r,z) = r*(\vec j\cdot\nabla\zeta) = EQPF*r/R_0+EQFF*R_0/r ,
C ASTRA average toroidal current density is
C	  R_0*<\vec j\cdot\nabla\zeta> = EQPF+EQFF*<R_0^2/r^2>
C made from RHSEQ (Polevoi 24-JUL-2025)
C-----------------------------------------------------------------------
	implicit none

	integer	j,NA,NA1
	double precision HRO,ROC,RTOR,BTOR,HROA,GP
      double precision
     >  NE(*),NI(*),TE(*),TI(*),PBLON(*),PBPER(*),PFAST(*),
     > G22(*),G33(*),IPOL(*),AMETR(*),
     > CUTOR(*),RHO(*),EQFF(*),EQPF(*),MU(*),CU(*)
	double precision	YCB,YG,YTH2,NB2EQL
        NA=NA1-1
        HRO=RHO(3)-RHO(2)
        ROC=RHO(NA1)
        NA=NA1-1
        HROA=ROC-RHO(NA)
        GP=3.14159263359d0
        HROA=HRO
        NB2EQL=1.
C Preparing input for the 3M equilibrium solver:
		YCB = 1.6E-3*RTOR/(BTOR*HRO*HRO)
	 do	J=2,NA
	   if (j.eq.NA) YCB = YCB*HRO/HROA
C	   EQFF(J) = -YCB/(MU(J))*(
C     +		.5*(PBLON(J+1)-PBLON(J)+PBPER(J+1)-PBPER(J))/(J+1)
C     +	        +( (NE(J+1)*TE(J+1)-NE(J)*TE(J))+
C     +		   (NI(J+1)*TI(J+1)-NI(J)*TI(J)) )/J )
	   EQFF(J) = ( (NE(J+1)*TE(J+1)-NE(J)*TE(J))+
     +		       (NI(J+1)*TI(J+1)-NI(J)*TI(J)) )/J
	   EQFF(J) = EQFF(J)+0.5*NB2EQL*
     +			(PBLON(J+1)-PBLON(J)+PBPER(J+1)-PBPER(J))/J
	   EQFF(J) = EQFF(J)+(PFAST(J+1)-PFAST(J))/J
	   EQFF(J) = -YCB*EQFF(J)/(MU(J))
	  enddo
		EQFF(1)	=EQFF(2)
		EQFF(NA1)=EQFF(NA)+(EQFF(NA)-EQFF(NA-1))
     .		*(AMETR(NA)-AMETR(NA-1))/(AMETR(NA1)-AMETR(NA))
		EQFF(NA1)=EQFF(NA)+(EQFF(NA)-EQFF(NA-1))/HRO*HROA
	  do	J=1,NA1
	   EQPF(j) = EQFF(j)
	   YTH2	= RHO(j)*G22(J)*(MU(J)/RTOR)**2
	   YG	= (1.+YTH2)*G33(J)
	   EQFF(J) = (CU(J)/IPOL(J)-EQPF(J))/YG
C	   if (J.eq.NA1) EQFF(J) = -EQPF(J) /YG
	   CUTOR(J) = (CU(J)/IPOL(J)+YTH2*EQPF(J))/(1.+YTH2)
	  enddo
	  return
	  end
C======================================================================
C FTLLMR [%]:
C Effective trapped particle fraction Lin-Liu and Miller GA-A21820
C                                                        (Oct 94)
C Y.R.Lin-Liu and R.L.Miller, Phys.Plasmas 2(5), May 1995, pp.1666-1668
C made from fnc/ftllm.f 	Polevoi 21-JUL-2025
	double precision function FTLLMRy(BTOR,BMAXT,BDB0,BDB02,FOFB)
	implicit none
	double precision YR,YYR,YH,YFTUP,YFTLO,BTOR,BMAXT,BDB0,BDB02,FOFB
        FTLLMRy=0.1
        IF (ABS(BMAXT).LT.0.0001 .OR.
     +      ABS(BDB0).LT.0.0001) THEN
           FTLLMRy=0.1
           GOTO 99
        ENDIF
        YH=min(.999999d0,BDB0/BMAXT*BTOR)
        YFTUP=1.-BDB02/BDB0**2*
     +        (1.-SQRT(1.-YH)*(1.+.5*YH))
        YFTLO=1.-BDB02*FOFB
		FTLLMRy=.75*YFTUP+.25*YFTLO
 99	  RETURN
	  END
!========================
C======================================================================|
	SUBROUTINE	ECH2a(YR0,YDR,YQ,YEFF1,YP,YC,NA1,RHO,VR,G33,IPOL,TE,NE)
!			RHOEC,RHODR,QEC,YEFFec,PEECR,CUECR,NA1,RHO,VR,G33,IPOL,TE,NE
C-----------------------------------------30-JUL-2025
C	Artificial parabolic Heating + CD		POLEVOY
C	YR0 	distance from plasma centre [m]
C	YDR	half width [m]
C	Q	ECH Power [MW]
C	YP(r)	ECR power density [MW/m3] ~ Q*exp(-((r-r0)/Dr)**2)
C	YC(r)	ECR current density [MW/m3]
C	YEFF	ECR current drive efficiency [A/W]
	implicit none
	integer j,NA1,jres,NA
	double precision YP(*),YC(*),RHO(*),VR(*),G33(*),IPOL(*),TE(*),NE(*)
!	double precision IINTa,VINT
	double precision IINTa,VINTa,YS,YM,YR0,YDR,YQ,YEFF,ROC,YEFF1
	external IINTa,VINTa
		YS=0.d0
		YM=0.d0
		NA=NA1-1
		ROC=RHO(NA1)
	if(ydr.le.1.d-7) then
	write(*,*) 'sbr/ech2a: too small iput power width'
!	 	YDR=ABC
		YDR=RHO(NA1)
	endif
	DO J	=1,NA
		YP(J)=-((RHO(J)-YR0)/YDR)**2
	if((RHO(J)/ROC-YR0)*(RHO(J+1)/ROC-YR0).le.0.) jres=j
	if(YP(J).gt.-18.d0) then
		YP(J)=dexp(YP(J))
	else
		YP(J)=0.d0
	endif
	enddo
         if(jres.lt.1) jres=1
         if(jres.ge.na1) then
            YP(1:na1)	=0.
            YC(1:na1)	=0.
            write(*,*)
     >  'Warning (ECH2a in SMTH.f): ECR location is out of plasma. Rec =',YR0,'> ROC=',RHO(NA1)
         return
         endif
		YM=IINTa(YP,ROC,RHO,G33,IPOL,NA1)
		YS=VINTa(YP,ROC,RHO,VR,NA1)
!	write(*,*) 'YS=',YS
!	if(ym.le.1.d-7.or.ys.le.1.d-7) then
		YEFF=YEFF1*TE(jres)/NE(jres)
	if(ys.le.0.) then
		YP(1:na1)	=0.
	else
		YP(1:na1)	=YQ/(YS+1.d-7)*YP(1:na1)
	endif
	if(ym.le.0.) then
		YC(1:na1)	=0.
	else
		YC(1:na1)	=YP(1:na1)*YEFF/(YM+1.d-7)
	endif
	return
		RETURN
	END

C=======================================================================
C VINT:	Volume integral {0,R} of any array
C Only a radially dependent array may be the 1st parameter of the function
C Examples:
C    out\Vint(CAR3)	!Radial profile of CAR3 volume integral
C    out_Vint(CAR3,Ro); !Volume integral {0,Ro} of CAR3
C    out_Vint(CAR3B)    !Total volume integral of CAR3 (0,ROC)
C			(Yushmanov 26-DEC-90)
	double precision function VINTa(ARR,YR,RHO,VR,NA1)
	implicit none
	double precision ARR(*),VR(*),RHO(*),YR,YDR,YR1,HRO,HROA
	integer JK,J,NA1,NA
	if(YR.le.0.)	return
		NA=NA1-1
		HRO=RHO(3)-RHO(2)
		HROA=RHO(NA1)-RHO(NA)
	if(YR.le. RHO(NA)) 	then
		JK = YR/HRO+1
		YR1 =  YR
	else
		JK = NA
		YR1 = min(YR,RHO(NA)+.5d0*HROA)
	endif
		YDR=(JK-YR1/HRO)*VR(JK)
		VINTa=0.
	do 1 J=1,JK
 1		VINTa=VINTa+ARR(J)*VR(J)
		VINTa=HRO*(VINTa-ARR(JK)*YDR)
	end
C======================================================================|
C======================================================================|
C IINTa:	Integral {0,R} of current density
C 	Iint=integral {0,R} (ARR/IPOL**2)dV*IPOL/(GP2*Ro)
C 	Only radial dependent array may be a parameter of the function
C 	Examples:
C    out\Iint(CU)	!Radial profile of toroidal current
C    out_Iint(CD,Ro)    !Toroidal driven current inside {0,Ro}
C    out_Iint(CUB)      !Total toroidal current =Iint(CU,ROC); (=IPL)
C			(Pereverzev 23-OCT-99)
	double precision function IINTa(ARR,YR,RHO,G33,IPOL,NA1)
	implicit none
	integer J,JK,NA1,NA
	double precision
     > ARR(*),RHO(*),G33(*),IPOL(*),YR,YDR,YA,GP2,GP,HRO,HROA,ROC
        NA=NA1-1
        HRO=RHO(3)-RHO(2)
        ROC=RHO(NA1)
        NA=NA1-1
        HROA=ROC-RHO(NA)
        GP=3.14159263359d0
        GP2=2.*GP
		IINTa = 0.
	if (YR .le. 0.)	return
		JK = YR/HRO+1.-1.E-4
	if (JK .gt. NA)	JK = NA
		YA = 0.
	do   1	J=1,JK
		IINTa = IINTa+YA
		YA   = ARR(J)*RHO(J)/(G33(J)*IPOL(J)**3)
		YDR = YR-JK*HRO+HRO
 1	continue
	if (JK .ge. NA)	then
	   YDR = min(YDR,0.5*(HRO+HROA))
	endif
		IINTa = GP2*IPOL(JK)*(HRO*IINTa+YDR*YA)
	end
C=======================================================================
C SVD1 [10^19m^3/s]:	The formula is a fit to D-D reaction rate
C	according to Putvinskiy
C	D+D=n(2.452MeV)+3He(0.817MeV)
C	Use:	P3He=Nd*Nd*SVD1*817./625.	[MW/m^3]
C	Use:	Neutron_Source=Nd*Nd*SVD1	[10^19/m^3/s](Neutrons)
C			(Yushmanov 11-JUN-87)
	double precision function SVD1y(TI)
	implicit none
	double precision YSVD,Ti
		YSVD=TI**.333333
		SVD1y=0.16247+0.001741*TI-0.029*EXP(-0.3843*SQRT(TI))
		SVD1y=SVD1y*EXP(-18.8085/YSVD)/YSVD**2
	return
	end
C SVD2 [10^19m^3/s]:	The formula is a fit to D-D reaction rate
C	according to Putvinskiy
C	d+d=t(1.008MeV)+p(3.025MeV)
C	Use:	P(p+t)=Nd*Nd*SVD2*4033./625.	[MW/m^3]
C	Use:	Neutron_Source=Nd*Nd*SVD1	[10^19/m^3/s](Neutrons)
C			(Yushmanov 11-JUN-87)
	double precision function SVD2y(TI)
	implicit none
	double precision YSVD,Ti
		YSVD=TI**.333333
		SVD2y=0.16052+0.001176*TI-0.01877*EXP(-0.3807*SQRT(TI))
		SVD2y=SVD2y*EXP(-18.8085/YSVD)/YSVD**2
	return
	end
C SVDBH [10#-19m#3/s]:	The formula is a fit to D-D reaction rate
C	according to Bosch and Hale, NUC. FUS.,32,p. 611, (1992)
C	D+D=n(2.452MeV)+3He(0.817MeV)
C	Use:	P3He=Nd*Nd*SVD1*817./625. [MW/m#3]
C	Use:	Neutron_Source=Nd*Nd*SVD1	[10^19/m^3/s](Neutrons)
C			(Stober 15-APR-99)
	double precision function SVDBHy(TI)
	implicit none
	double precision YSVD,Ti
        SVDBHy=1.0+TI*(7.68222E-3-TI*2.96400E-6)
        SVDBHy=TI/(1.0-TI*5.85778E-3/SVDBHy)
        YSVD =(31.3970**2/4.0/SVDBHy)**.333333
        SVDBHy=5.43360E-12*SVDBHy*sqrt(YSVD/937814.0/TI**3)
        SVDBHy=SVDBHy*EXP(-3.0*YSVD)*1.0E13
	return
	end
!==============================================================
C	stbrn[a.u.]
C	Fprb, Probability for 1.008 MeV T NBI to burn out on D target
C	approximation of cold Deuterium (Ti = 0) DT (50:50 plasmas)
C	crossection by S.V.Putvinskij//VANT,v.2,1988,p.3
C	Corrected by S.V.Putvinskij 22.11.89
C	Use:		Fprb = stbrn*Nd	[1/m#3]
C  		 D-T	Fusion power thermal D  + 1.008 MeV T fast
C
C				Polevoy		06-JUL-99,08-JUN-07
c	ABEAM = 3.
c	EBEAM = 1008.

	double precision function STBRNy(TE,NE,Z2NdA)
	implicit none
	integer j,jk
	double precision TE,NE,Z2NdA
	double precision YECM,YSQ,YECDEB,YXC3,YX,YX2,YX3,YS2,YE,YSIG
		YECM	=1008.d0*.4d0
		YSQ	=DSQRT(YECM)
		stbrny	=0.d0
		YECDEB  = 14.6d0*TE*3.d0/1008.d0*(Z2NdA/NE)**0.667
		YXC3 =YECDEB*DSQRT(YECDEB)
	DO 	Jk = 1,1000
		YX	=(1.d-3*Jk)
		YX2	=YX*YX
		YX3	=YX2*YX
		YE	=YECM*YX2
		YS2	=(YE-48.7878d0)**2
		YSIG=
     >		DEXP(-34.3812d0/YX/YSQ)*(1.d0+1.1177d-5*YS2)/
     >		(1.d0+6.433d-4*YS2)/YE
		stbrny	=stbrny+ YSIG/(1.d0+YXC3/YX3)
	 enddo
	stbrny	=stbrny*26798.d0*1.d-3*
     * 	4.38d-4*DSQRT(1008.d0/3.d0)*
     *	2.d0*3.d0/(NE+1.d-9)*TE*DSQRT(TE)/
     /	(15.85d0+DLOG(TE/DSQRT(NE+1.d-9)))
		return
		end
!!=========================================================
      double precision function svcxy(TI,AMAIN)
      implicit none
      double precision TI,AMAIN
		svcxy=0.
      if (ti.gt.0. and. amain.ge.1.) then
		SVCXy	=10.**(5.9+0.3*LOG10(TI/AMAIN))
      endif
      return
      end
!!=========================================================
      double precision function SVRCy(TE)
         implicit none
         double precision TE
      IF(TE.LT..0001) THEN
         SVRCy=0.
      ELSE
		SVRCy	=13.6E-3/TE
		SVRCy	=1.27*SVRCy*sqrt(SVRCy)/(SVRCy+.59)
      ENDIF
         return
         end
!!=========================================================
      double precision function SVIEy(TE)
      implicit none
		double precision	TE
		SVIEy=0.
      IF(TE.GT..01)	THEN
		SVIEy	=.0136/TE
      IF(TE.GT..01)	THEN
		SVIEy=
     >	9.7E5*EXP(-SVIEy)*SQRT(SVIEy/(1.+SVIEy))/(SVIEy+.73)
				ELSE
		SVIEy=2.958E5*EXP(-SVIEy)*SQRT(SVIEy)
				ENDIF
				endif
      return
      end
!=========================================================
C PENLI [MW/m#3]: Electron heat losses due to cold neutrals Radiation
C	PENLI=<sigma*v>*Ne*Nn*(0.0102keV)
C	PE=...-PENLI = ...-	PENLIy*NE*f0
C			(Pereverzev 09-JULY-90)
	double precision function PENLIy(TE)
	implicit none
	double precision TE,YY
		PENLIy=0.
	if(TE.gt.0.) then
		YY	=.0102/TE
		PENLIy	=0.48/(0.28+YY)*SQRT(YY*(1.+YY))*EXP(-YY)
		PENLIy	=16.3*PENLIy
	endif
		return
		end
!========================================================
C PAIONy [MW/m#3]   D-T Fraction of fusion alpha power deposited to ions
C     P.Pavlo  22.06.89/A.Polevoi 20-MAY-94
	double precision function PAIONy(TE,EBDECTE)
	implicit none
	double precision TE,EBDECTE,y,y2
		PAIONy=1.
	if(TE.gt.0.) then
        y2 = EBDECTE/TE
        y = sqrt(y2)
        PAIONy    = 2.* (  0.166666667*LOG( (1.-y+y2)/(1.+2.*y+y2) ) +
     + 0.57735026*(ATAN( 0.57735026*(2.*y-1.) )+0.52359874) ) / y2
	endif
	return
	end
C======================================================================|
	subroutine	CUOFPy(NA1,RHO,RTOR,BTOR,FP,MU,CU,G22,G33,IPOL)
C----------------------------------------------------------------------|
C Compute CU(rho) and MU(rho) from FP(rho)
C----------------------------------------------------------------------|
C	Input:	HRO	- radial step (m)
C		HROA	- edge radial step (m)
C		NA1	- number of grid points
C		G22(1:NA)	- <g22/g>*............
C		G33(1:NA)	-
C		IPOL(1:NA)	-
C		CD(1:NA)	- external (+bootstrap) current
C		FP(1:NA1)	- poloidal flux
C	Output:	CU(1:NA1) - (1/rho)d{K*dF/d(rho)}/d(rho) current density
C		MU(1:NA1) - (1/rho)dF/d(rho)	    rotational transform
! made from CUOFP Polevoi 22-JUL-2025
C----------------------------------------------------------------------|
	implicit none
!	include	'for/parameter.inc'
!	include	'for/const.inc'
!	include	'for/status.inc'
	integer	j,NA,NA1
	double precision	HH,YAJ,YCJ,ARRNA1
	double precision FP(*),MU(*),CU(*),G22(*),G33(*),IPOL(*),RHO(*)
	double precision HRO,HROA,RTOR,BTOR,GP,ROC
      external ARRNA1
		NA=NA1-1
        HRO=RHO(3)-RHO(2)
        ROC=RHO(NA1)
        NA=NA1-1
        HROA=ROC-RHO(NA)
        GP=3.14159263359d0
        HROA=HRO

        HH = HRO*HRO
        YAJ = 0.
C	write(*,*)"CU"
C	write(*,100)(CU(j),j=NA1-5,NA1)
	do	1	J=1,NA
	   YCJ = YAJ
	if (j .lt. NA)	then
	      YAJ = (FP(j+1)-FP(j))/HH
	      MU(j) = YAJ/j
	      YAJ = G22(j)*YAJ
	      CU(j) = (YAJ-YCJ)/HRO
	else
	      YAJ = (FP(j+1)-FP(j))/HRO/HROA
	      MU(j) = YAJ/j
	      YAJ = G22(j)*YAJ
C30-11	      CU(j) = 2.*(YAJ-YCJ)/(HRO+HROA)
	      CU(j) = (YAJ-YCJ)/HRO
	endif

		CU(j) = CU(j)/(j-0.5)
 1	continue
      MU(1) =MU(2)
C Suppress the next line to avoid CU(NA1) =/= 0
		CU(NA1) = ARRNA1(CU(NA-2),HROA/HRO)
		HH = HROA/(HRO+HROA)		! Attributes MU(NA1)
		MU(NA1) = ARRNA1(MU(NA-2),HH)	! to the point ROC
C	MU(NA1) = MU(NA)*NA/(NA-0.5+HROA/HRO)
		YCJ = 1.25/(GP*GP*RTOR)
		YAJ = 0.5/(GP*BTOR)
	do	2	J=1,NA1
		CU(j) = YCJ*CU(j)*G33(J)*IPOL(J)**3
		MU(J) = YAJ*MU(j)
 2	continue
C	CU(NA)=CUBS(NA)+CD(NA)+CC(NA)*ULON(NA)/(RTOR*GP2)
C	CU(NA1)=CUBS(NA1)+CD(NA1)+CC(NA1)*ULON(NA1)/(RTOR*GP2)
C	write(*,100)Time
C	write(*,100)(CU(j),j=NA1-5,NA1)
C	write(*,100)(MU(j),j=NA1-5,NA1)
 100	format(3(2F10.5,2X))
	end
C======================================================================|
	double precision function ARRNA1(ARR,H)
C----------------------------------------------------------------------|
C Compute edge value of the array ARR(NA1) assuming that d^3(ARR)/dr^3=0
C----------------------------------------------------------------------|
	double precision ARR(3),H
C	ARRNA1 = ARR(1)+H*(2.*ARR(1)+ARR(-2)-3.*ARR(-1))
		ARRNA1 = ARR(3)*(1+(1.5+0.5*H)*H)+ARR(1)*H*0.5*(1+H)-
     >		ARR(2)*H*(2+H)
	end
C======================================================================|
      Subroutine CUBSy(NA1,RTOR,BTOR,IPL,
     > FP,MU,ZEF,TE,TI,NE,NI,AMAIN,ZMAIN,
     > BMINT,BMAXT,BDB0,BDB02,FOFB,SQEPS,RHO,
     > CUBS,CC) !output bootstrap [MA/m2] current conductivity [1/mkOM*m]
! Bootsrap current and current conductvity by Sauter
!                       Sauter, Angioni, Lin-Liu
!			Physics of Plasmas, 6, 2834 (1999)
!	(Polevoi 22-JUL-2025)
    	implicit none
	integer	j,NA,NA1
	double precision
     > YC,YD,YA,ZZ,ZDF,ZFT,ZFTE,ZFTE1,ZFTE2,ZFTE3,ZFTE4,
     > ZFTI,ZFTI1,ZFTI2,ZFTI3,ZFTI4,A0,A1,ALP,
     > BETPL,COULG,NUEE,NUES,NUI,NUIS,HCEE,HCEI,
     > DCSA,HCSA,XCSA,CCSP,CNSA
	double precision
     > FP(*),MU(*),CUBS(*),CC(*),ZEF(*),TE(*),TI(*),NE(*),NI(*),
     > BMINT(*),BMAXT(*),BDB0(*),BDB02(*),FOFB(*),
     > RHO(*),SQEPS(*),AMAIN(*),ZMAIN(*)
	double precision HRO,HROA,RTOR,BTOR,GP,FTLLMRy,HC,XC,DC,ROC,IPL
	external FTLLMRy
        HRO=RHO(3)-RHO(2)
        ROC=RHO(NA1)
        NA=NA1-1
        HROA=ROC-RHO(NA)
        GP=3.14159263359d0
!    call	markloc("Current equation"//char(0))
      YC=.4*GP*RTOR/BTOR
      YD=-0.8*GP*GP*RTOR
      YA=2./(HRO**2*YD)
!         write(*,*) 'cu',
!     > FTLLMRy(BTOR,BMAXT(1),BDB0(1),BDB02(1),FOFB(1))
      do 2224 J=1,NA1
	if(J.ne.NA1)	then
	   BETPL = 1.6E-4*GP*(NE(J+1)+NE(J))*(TE(J+1)+TE(J))
	   BETPL = BETPL*(RTOR/(BTOR*J*HRO*ABS(MU(J)+.0001)))**2
	else
	   BETPL = 6.4E-4*GP*NE(J)*TE(J)
	   BETPL = BETPL*(RTOR/(BTOR*ROC*ABS(MU(J)+.0001)))**2
	endif
	ZZ=ZEF(J)
        ZFT=FTLLMRy(BTOR,BMAXT(J),BDB0(J),BDB02(J),FOFB(J))
!                 write(*,*) 'j',j, ZFT

	COULG=15.9d0-0.5d0*LOG(NE(J))+log(TE(J))
	NUEE=670.*COULG*NE(J)/(TE(J)*SQRT(TE(J)))
	NUES=NUEE*1.4*ZEF(J)*RTOR
     *	/ABS(MU(J)*(SQEPS(J)**3)*SQRT(TE(J))*1.875E7)
	ZDF = 1.+(1.-0.1*ZFT)*SQRT(NUES)
	ZDF = ZDF + 0.5*(1.-ZFT)*NUES/ZZ
	ZFT = ZFT/ZDF
!	               write(*,*) 'j', j, ZFT
       	DCSA = (1+1.4/(ZZ+1))*ZFT-1.9/(ZZ+1)*ZFT*ZFT
        DCSA = DCSA + (0.3*ZFT*ZFT+0.2*ZFT*ZFT*ZFT)*ZFT/(ZZ+1)
	DCSA=DCSA*BETPL*(1.+TI(J)/(ZZ*TE(J)))
!      DC(J)=DCSA
!      write(*,*) 'j, DCSA ', j, DCSA
        DC=DCSA
!	if(J.ne.NA1)	then
!	   BETPL = 1.6E-4*GP*(NE(J+1)+NE(J))*(TE(J+1)+TE(J))
!	   BETPL = BETPL*(RTOR/(BTOR*J*HRO*ABS(MU(J)+.0001)))**2
!	else
!	   BETPL = 6.4E-4*GP*NE(J)*TE(J)
!	   BETPL = BETPL*(RTOR/(BTOR*ROC*ABS(MU(J)+.0001)))**2
!	endif
	ZZ=ZEF(J)
        ZFT=FTLLMRy(BTOR,BMAXT(J),BDB0(J),BDB02(J),FOFB(J))
!	COULG=15.9d0-0.5d0*LOG(NE(J))+log(TE(J))
	NUEE=670.*COULG*NE(J)/(TE(J)*SQRT(TE(J)))
	NUES=NUEE*1.4*ZEF(J)*RTOR
     *	/ABS(MU(J)*(SQEPS(J)**3)*SQRT(TE(J))*1.875E7)
	ZDF = 1.+0.26*(1.-ZFT)*SQRT(NUES)
	ZDF = ZDF + 0.18*(1.-0.37*ZFT)*NUES/SQRT(ZZ)
        ZFTE = ZFT/ZDF
        ZFTE2=ZFTE*ZFTE
        ZFTE3=ZFTE*ZFTE2
	ZFTE4=ZFTE2*ZFTE2
        ZDF = 1.+(1.+0.6*ZFT)*SQRT(NUES)
	ZDF = ZDF + 0.85*(1.-0.37*ZFT)*NUES*(1+ZZ)
	ZFTI = ZFT/ZDF
        ZFTI2=ZFTI*ZFTI
        ZFTI3=ZFTI*ZFTI2
	ZFTI4=ZFTI2*ZFTI2
	HCEE = (0.05+0.62*ZZ)/ZZ/(1.+0.44*ZZ)*(ZFTE-ZFTE4)
	HCEE = HCEE +(ZFTE2-ZFTE4-1.2*(ZFTE3-ZFTE4))/(1+0.22*ZZ)
	HCEE = HCEE + 1.2/(1+0.5*ZZ)*ZFTE4
        HCEI = -(0.56+1.93*ZZ)/ZZ/(1.+0.44*ZZ)*(ZFTI-ZFTI4)
	HCEI = HCEI+4.95/(1+2.48*ZZ)*(ZFTI2-ZFTI4-0.55*(ZFTI3-ZFTI4))
	HCEI = HCEI-1.2/(1+0.5*ZZ)*ZFTI4
!	if(J.ne.NA1)	then
!	   BETPL = 1.6E-4*GP*(NE(J+1)+NE(J))*(TE(J+1)+TE(J))
!	   BETPL = BETPL*(RTOR/(BTOR*J*HRO*ABS(MU(J)+.0001)))**2
!	else
!	   BETPL = 6.4E-4*GP*NE(J)*TE(J)
!	   BETPL = BETPL*(RTOR/(BTOR*ROC*ABS(MU(J)+.0001)))**2
!	endif
	ZZ=ZEF(J)
        ZFT=FTLLMRy(BTOR,BMAXT(J),BDB0(J),BDB02(J),FOFB(J))
!	COULG=15.9d0-0.5d0*LOG(NE(J))+log(TE(J))
	NUEE=670.*COULG*NE(J)/(TE(J)*SQRT(TE(J)))
	NUES=NUEE*1.4*ZEF(J)*RTOR
     *	/ABS(MU(J)*(SQEPS(J)**3)*SQRT(TE(J))*1.875E7)
	ZDF = 1.+(1.-0.1*ZFT)*SQRT(NUES)
	ZDF = ZDF + 0.5*(1.-ZFT)*NUES/ZZ
	ZFT = ZFT/ZDF
       	DCSA = (1.+1.4/(ZZ+1))*ZFT-1.9/(ZZ+1)*ZFT*ZFT
        DCSA = DCSA + (0.3*ZFT*ZFT+0.2*ZFT*ZFT*ZFT)*ZFT/(ZZ+1.)
	DCSA=DCSA*BETPL*(1.+TI(J)/(ZZ*TE(J)))
        HCSA = BETPL*(HCEE+HCEI)+DCSA/(1.+TI(J)/(ZZ*TE(J)))
!      HC(J)=HCSA
!            write(*,*) 'j, HCSA ', j, HCSA
      HC=HCSA
!	if(J.ne.NA1)	then
!	   BETPL = 1.6E-4*GP*(NE(J+1)+NE(J))*(TE(J+1)+TE(J))
!	   BETPL = BETPL*(RTOR/(BTOR*J*HRO*ABS(MU(J)+.0001)))**2
!	else
!	   BETPL = 6.4E-4*GP*NE(J)*TE(J)
!	   BETPL = BETPL*(RTOR/(BTOR*ROC*ABS(MU(J)+.0001)))**2
!	endif
	ZZ=ZEF(J)
        ZFT=FTLLMRy(BTOR,BMAXT(J),BDB0(J),BDB02(J),FOFB(J))
!	COULG=15.9d0-0.5d0*LOG(NE(J))+log(TE(J))
	NUEE=670.*COULG*NE(J)/(TE(J)*SQRT(TE(J)))
	NUES=NUEE*1.4*ZEF(J)*RTOR
     *	/ABS(MU(J)*(SQEPS(J)**3)*SQRT(TE(J))*1.875E7)
	NUI=ZMAIN(J)**4*NI(J)*322./(TI(J)*SQRT(TI(J)*AMAIN(J)))
	NUIS=3.2E-6*NUI*RTOR/(ABS(MU(J)+.0001)*
     *		SQEPS(J)**3*SQRT(TI(J)/AMAIN(J)))
	ZDF = 1.+(1.-0.1*ZFT)*SQRT(NUES)
	ZDF = ZDF + 0.5*(1.-0.5*ZFT)*NUES/ZZ
        ZFTE = ZFT/ZDF
	XCSA = (1+1.4/(ZZ+1))*ZFTE-1.9/(ZZ+1.)*ZFTE*ZFTE
        XCSA = XCSA + (0.3*ZFTE*ZFTE+0.2*ZFTE*ZFTE*ZFTE)*ZFTE/(ZZ+1)
	A0 = -(1.17)*(1.-ZFT)
        A0 =A0/(1.-0.22*ZFT-0.19*ZFT*ZFT)
	ALP =(A0+0.25*(1-ZFT*ZFT)*SQRT(NUIS))/(1.+0.5*SQRT(NUIS))
        A1 = NUIS*NUIS*ZFT*ZFT*ZFT*ZFT*ZFT*ZFT
        ALP = (ALP + 0.315*A1)/(1.+0.15*A1)
!	if(J.ne.NA1)	then
!	   BETPL = 1.6E-4*GP*(NE(J+1)+NE(J))*(TE(J+1)+TE(J))
!	   BETPL = BETPL*(RTOR/(BTOR*J*HRO*ABS(MU(J)+.0001)))**2
!	else
!	   BETPL = 6.4E-4*GP*NE(J)*TE(J)
!	   BETPL = BETPL*(RTOR/(BTOR*ROC*ABS(MU(J)+.0001)))**2
!	endif
	ZZ=ZEF(J)
        ZFT=FTLLMRy(BTOR,BMAXT(J),BDB0(J),BDB02(J),FOFB(J))
!	COULG=15.9d0-0.5d0*LOG(NE(J))+log(TE(J))
	NUEE=670.*COULG*NE(J)/(TE(J)*SQRT(TE(J)))
	NUES=NUEE*1.4*ZEF(J)*RTOR
     *	/ABS(MU(J)*(SQEPS(J)**3)*SQRT(TE(J))*1.875E7)
	ZDF = 1.+(1.-0.1*ZFT)*SQRT(NUES)
	ZDF = ZDF + 0.5*(1.-ZFT)*NUES/ZZ
	ZFT = ZFT/ZDF
       	DCSA = (1.+1.4/(ZZ+1.))*ZFT-1.9/(ZZ+1.)*ZFT*ZFT
        DCSA = DCSA + (0.3*ZFT*ZFT+0.2*ZFT*ZFT*ZFT)*ZFT/(ZZ+1)
	DCSA=DCSA*BETPL*(1.+TI(J)/(ZZ*TE(J)))
        XCSA = BETPL*(XCSA*ALP)*TI(J)/ZZ/TE(J)
        XCSA = XCSA+DCSA/(1.+ZZ*TE(J)/TI(J))
        XC=XCSA
!      XC(J)=XCSA
!      CD(J)=CUBM(J)+CUECR(J)+CUICR(J)+CULH(J)+CUFW(J)+CAR8(J)
!     + +CAR18(J)+CAR28(J)
	ZZ=ZEF(J)
        ZFT=FTLLMRy(BTOR,BMAXT(J),BDB0(J),BDB02(J),FOFB(J))
!	COULG=15.9d0-0.5d0*LOG(NE(J))+log(TE(J))
	NUEE=670.*COULG*NE(J)/(TE(J)*SQRT(TE(J)))
	NUES=NUEE*1.4*ZEF(J)*RTOR
     *	/ABS(MU(J)*(SQEPS(J)**3)*SQRT(TE(J))*1.875E7)
        ZDF = 1.+(0.55-0.1*ZFT)*SQRT(NUES)
	ZDF = ZDF + 0.45*(1.-ZFT)*NUES/ZZ/SQRT(ZZ)
        ZFT = ZFT/ZDF
	CNSA=1.-(1.+0.36/ZZ)*ZFT
        CNSA=CNSA+0.59/ZZ*ZFT*ZFT-0.23/ZZ*ZFT*ZFT*ZFT
!	COULG=15.9d0-0.5d0*LOG(NE(J))+log(TE(J))
	CCSP = 1041.d0/COULG*TE(J)*SQRT(TE(J))*
     *		(1.13d0+ZEF(J))/(2.67d0+ZEF(J))/ZEF(J)
	CNSA=CCSP*CNSA
      CC(J)=CNSA
!      write(*,*) 'j, CC', j , CC(j)
      if (j .eq. NA1) goto 2222
      if (j .eq. NA) YA = 2./(YD*HROA**2)
      CUBS(J)=YA*(FP(J+1)-FP(J))*(0.
     > +HC*(TE(J+1)-TE(J))/(TE(J+1)+TE(J))
     > +XC*(TI(J+1)-TI(J))/(TI(J+1)+TI(J))
     > +DC*(NE(J+1)-NE(J))/(NE(J+1)+NE(J))
     > )
      if(j.eq.NA1) CUBS(J)= max(0.,(2*CUBS(NA)-CUBS(NA-1)))
      goto 2223
 2222    continue
 2223    continue
 2224    continue
         return
         end
!=============================
C======================================================================|
	subroutine	RHSEQy(NA1,RTOR,BTOR,RHO,
     > NE,NI,TE,TI,PBLON,PBPER,PFAST,MU,CU,
     > G22,G33,IPOL,AMETR,
     > CUTOR,EQFF,EQPF)
C-----------------------------------------------------------------------
C Input: RTOR,BTOR,NA,NA1,HRO,HROA,   NB2EQL,
C	 NE,NI,TE,TI,MU,CU,AMETR,RHO, PBLON,PBPER,G22,G33,IPOL
C Output:
C       EQPF
C       EQFF
C       CUTOR
C Both quantities EQPF (~p') and EQFF (~II') are given in [MA/m^2]
C	EQPF = -1.6E-3*(2*\pi*R_0)\prti{n_13*T_keV}{\psi[Vs=T*m^2]}
C	     = -1.E-6*/(2*\pi*R_0)\prti{p[J/m^3=Pascal]}{\psi[Vs]}
C	EQFF = -1.E-6*2*\pi/(R_0*\mu_0)*I*\prti{I}{\psi}
C	     = -5./R_0*I*\prti{I}{\psi}
C Local toroidal current density j[MA/m^2] is EQPF*r/R_0+EQFF*R_0/r, i.e.
C       j(r,z) = r*(\vec j\cdot\nabla\zeta) = EQPF*r/R_0+EQFF*R_0/r ,
C ASTRA average toroidal current density is
C	  R_0*<\vec j\cdot\nabla\zeta> = EQPF+EQFF*<R_0^2/r^2>
C made from RHSEQ (Polevoi 24-JUL-2025)
C-----------------------------------------------------------------------
	implicit none

	integer	j,NA,NA1
	double precision HRO,ROC,RTOR,BTOR,HROA,GP
      double precision
     >  NE(*),NI(*),TE(*),TI(*),PBLON(*),PBPER(*),PFAST(*),
     > G22(*),G33(*),IPOL(*),AMETR(*),
     > CUTOR(*),RHO(*),EQFF(*),EQPF(*),MU(*),CU(*)
	double precision	YCB,YG,YTH2,NB2EQL
        NA=NA1-1
        HRO=RHO(3)-RHO(2)
        ROC=RHO(NA1)
        NA=NA1-1
        HROA=ROC-RHO(NA)
        GP=3.14159263359d0
        HROA=HRO
        NB2EQL=1.
C Preparing input for the 3M equilibrium solver:
		YCB = 1.6E-3*RTOR/(BTOR*HRO*HRO)
	 do	J=2,NA
	   if (j.eq.NA) YCB = YCB*HRO/HROA
C	   EQFF(J) = -YCB/(MU(J))*(
C     +		.5*(PBLON(J+1)-PBLON(J)+PBPER(J+1)-PBPER(J))/(J+1)
C     +	        +( (NE(J+1)*TE(J+1)-NE(J)*TE(J))+
C     +		   (NI(J+1)*TI(J+1)-NI(J)*TI(J)) )/J )
	   EQFF(J) = ( (NE(J+1)*TE(J+1)-NE(J)*TE(J))+
     +		       (NI(J+1)*TI(J+1)-NI(J)*TI(J)) )/J
	   EQFF(J) = EQFF(J)+0.5*NB2EQL*
     +			(PBLON(J+1)-PBLON(J)+PBPER(J+1)-PBPER(J))/J
	   EQFF(J) = EQFF(J)+(PFAST(J+1)-PFAST(J))/J
	   EQFF(J) = -YCB*EQFF(J)/(MU(J))
	  enddo
		EQFF(1)	=EQFF(2)
		EQFF(NA1)=EQFF(NA)+(EQFF(NA)-EQFF(NA-1))
     .		*(AMETR(NA)-AMETR(NA-1))/(AMETR(NA1)-AMETR(NA))
		EQFF(NA1)=EQFF(NA)+(EQFF(NA)-EQFF(NA-1))/HRO*HROA
	  do	J=1,NA1
	   EQPF(j) = EQFF(j)
	   YTH2	= RHO(j)*G22(J)*(MU(J)/RTOR)**2
	   YG	= (1.+YTH2)*G33(J)
	   EQFF(J) = (CU(J)/IPOL(J)-EQPF(J))/YG
C	   if (J.eq.NA1) EQFF(J) = -EQPF(J) /YG
	   CUTOR(J) = (CU(J)/IPOL(J)+YTH2*EQPF(J))/(1.+YTH2)
	  enddo
	  return
	  end
C======================================================================
C FTLLMR [%]:
C Effective trapped particle fraction Lin-Liu and Miller GA-A21820
C                                                        (Oct 94)
C Y.R.Lin-Liu and R.L.Miller, Phys.Plasmas 2(5), May 1995, pp.1666-1668
C made from fnc/ftllm.f 	Polevoi 21-JUL-2025
	double precision function FTLLMRy(BTOR,BMAXT,BDB0,BDB02,FOFB)
	implicit none
	double precision YR,YYR,YH,YFTUP,YFTLO,BTOR,BMAXT,BDB0,BDB02,FOFB
        FTLLMRy=0.1
        IF (ABS(BMAXT).LT.0.0001 .OR.
     +      ABS(BDB0).LT.0.0001) THEN
           FTLLMRy=0.1
           GOTO 99
        ENDIF
        YH=min(.999999d0,BDB0/BMAXT*BTOR)
        YFTUP=1.-BDB02/BDB0**2*
     +        (1.-SQRT(1.-YH)*(1.+.5*YH))
        YFTLO=1.-BDB02*FOFB
		FTLLMRy=.75*YFTUP+.25*YFTLO
 99	  RETURN
	  END
!========================
C======================================================================|
	SUBROUTINE	ECH2a(YR0,YDR,YQ,YEFF1,YP,YC,NA1,RHO,VR,G33,IPOL,TE,NE)
!			RHOEC,RHODR,QEC,YEFFec,PEC,CUECR,NA1,RHO,VR,G33,IPOL,TE,NE
C-----------------------------------------30-JUL-2025
C	Artificial parabolic Heating + CD		POLEVOY
C	YR0 	distance from plasma centre [m]
C	YDR	half width [m]
C	Q	ECH Power [MW]
C	YP(r)	ECR power density [MW/m3] ~ Q*exp(-((r-r0)/Dr)**2)
C	YC(r)	ECR current density [MW/m3]
C	YEFF	ECR current drive efficiency [A/W]
	implicit none
	integer j,NA1,jres,NA
	double precision YP(*),YC(*),RHO(*),VR(*),G33(*),IPOL(*),TE(*),NE(*)
!	double precision IINTa,VINT
	double precision IINTa,VINTa,YS,YM,YR0,YDR,YQ,YEFF,ROC,YEFF1
	external IINTa,VINTa
		YS=0.d0
		YM=0.d0
		NA=NA1-1
		ROC=RHO(NA1)
	if(ydr.le.1.d-7) then
	write(*,*) 'sbr/ech2a: too small iput power width'
!	 	YDR=ABC
		YDR=RHO(NA1)
	endif
	DO J	=1,NA
		YP(J)=-((RHO(J)-YR0)/YDR)**2
	if((RHO(J)/ROC-YR0)*(RHO(J+1)/ROC-YR0).le.0.) jres=j
	if(YP(J).gt.-18.d0) then
		YP(J)=dexp(YP(J))
	else
		YP(J)=0.d0
	endif
	enddo
		YM=IINTa(YP,ROC,RHO,G33,IPOL,NA1)
		YS=VINTa(YP,ROC,RHO,VR,NA1)
!	write(*,*) 'YS=',YS
!	if(ym.le.1.d-7.or.ys.le.1.d-7) then
		YEFF=YEFF1*TE(jres)/NE(jres)
	if(ys.le.0.) then
		YP(1:na1)	=0.
	else
		YP(1:na1)	=YQ/(YS+1.d-7)*YP(1:na1)
	endif
	if(ym.le.0.) then
		YC(1:na1)	=0.
	else
		YC(1:na1)	=YP(1:na1)*YEFF/(YM+1.d-7)
	endif
	return
		RETURN
	END

C=======================================================================
C VINT:	Volume integral {0,R} of any array
C Only a radially dependent array may be the 1st parameter of the function
C Examples:
C    out\Vint(CAR3)	!Radial profile of CAR3 volume integral
C    out_Vint(CAR3,Ro); !Volume integral {0,Ro} of CAR3
C    out_Vint(CAR3B)    !Total volume integral of CAR3 (0,ROC)
C			(Yushmanov 26-DEC-90)
	double precision function VINTa(ARR,YR,RHO,VR,NA1)
	implicit none
	double precision ARR(*),VR(*),RHO(*),YR,YDR,YR1,HRO,HROA
	integer JK,J,NA1,NA
	if(YR.le.0.)	return
		NA=NA1-1
		HRO=RHO(3)-RHO(2)
		HROA=RHO(NA1)-RHO(NA)
	if(YR.le. RHO(NA)) 	then
		JK = YR/HRO+1
		YR1 =  YR
	else
		JK = NA
		YR1 = min(YR,RHO(NA)+.5d0*HROA)
	endif
		YDR=(JK-YR1/HRO)*VR(JK)
		VINTa=0.
	do 1 J=1,JK
 1		VINTa=VINTa+ARR(J)*VR(J)
		VINTa=HRO*(VINTa-ARR(JK)*YDR)
	end
C======================================================================|
C======================================================================|
C IINTa:	Integral {0,R} of current density
C 	Iint=integral {0,R} (ARR/IPOL**2)dV*IPOL/(GP2*Ro)
C 	Only radial dependent array may be a parameter of the function
C 	Examples:
C    out\Iint(CU)	!Radial profile of toroidal current
C    out_Iint(CD,Ro)    !Toroidal driven current inside {0,Ro}
C    out_Iint(CUB)      !Total toroidal current =Iint(CU,ROC); (=IPL)
C			(Pereverzev 23-OCT-99)
	double precision function IINTa(ARR,YR,RHO,G33,IPOL,NA1)
	implicit none
	integer J,JK,NA1,NA
	double precision
     > ARR(*),RHO(*),G33(*),IPOL(*),YR,YDR,YA,GP2,GP,HRO,HROA,ROC
        NA=NA1-1
        HRO=RHO(3)-RHO(2)
        ROC=RHO(NA1)
        NA=NA1-1
        HROA=ROC-RHO(NA)
        GP=3.14159263359d0
        GP2=2.*GP
		IINTa = 0.
	if (YR .le. 0.)	return
		JK = YR/HRO+1.-1.E-4
	if (JK .gt. NA)	JK = NA
		YA = 0.
	do   1	J=1,JK
		IINTa = IINTa+YA
		YA   = ARR(J)*RHO(J)/(G33(J)*IPOL(J)**3)
		YDR = YR-JK*HRO+HRO
 1	continue
	if (JK .ge. NA)	then
	   YDR = min(YDR,0.5*(HRO+HROA))
	endif
		IINTa = GP2*IPOL(JK)*(HRO*IINTa+YDR*YA)
	end
C=======================================================================
C SVD1 [10^19m^3/s]:	The formula is a fit to D-D reaction rate
C	according to Putvinskiy
C	D+D=n(2.452MeV)+3He(0.817MeV)
C	Use:	P3He=Nd*Nd*SVD1*817./625.	[MW/m^3]
C	Use:	Neutron_Source=Nd*Nd*SVD1	[10^19/m^3/s](Neutrons)
C			(Yushmanov 11-JUN-87)
	double precision function SVD1y(TI)
	implicit none
	double precision YSVD,Ti
		YSVD=TI**.333333
		SVD1y=0.16247+0.001741*TI-0.029*EXP(-0.3843*SQRT(TI))
		SVD1y=SVD1y*EXP(-18.8085/YSVD)/YSVD**2
	return
	end
C SVD2 [10^19m^3/s]:	The formula is a fit to D-D reaction rate
C	according to Putvinskiy
C	d+d=t(1.008MeV)+p(3.025MeV)
C	Use:	P(p+t)=Nd*Nd*SVD2*4033./625.	[MW/m^3]
C	Use:	Neutron_Source=Nd*Nd*SVD1	[10^19/m^3/s](Neutrons)
C			(Yushmanov 11-JUN-87)
	double precision function SVD2y(TI)
	implicit none
	double precision YSVD,Ti
		YSVD=TI**.333333
		SVD2y=0.16052+0.001176*TI-0.01877*EXP(-0.3807*SQRT(TI))
		SVD2y=SVD2y*EXP(-18.8085/YSVD)/YSVD**2
	return
	end
C SVDBH [10#-19m#3/s]:	The formula is a fit to D-D reaction rate
C	according to Bosch and Hale, NUC. FUS.,32,p. 611, (1992)
C	D+D=n(2.452MeV)+3He(0.817MeV)
C	Use:	P3He=Nd*Nd*SVD1*817./625. [MW/m#3]
C	Use:	Neutron_Source=Nd*Nd*SVD1	[10^19/m^3/s](Neutrons)
C			(Stober 15-APR-99)
	double precision function SVDBHy(TI)
	implicit none
	double precision YSVD,Ti
        SVDBHy=1.0+TI*(7.68222E-3-TI*2.96400E-6)
        SVDBHy=TI/(1.0-TI*5.85778E-3/SVDBHy)
        YSVD =(31.3970**2/4.0/SVDBHy)**.333333
        SVDBHy=5.43360E-12*SVDBHy*sqrt(YSVD/937814.0/TI**3)
        SVDBHy=SVDBHy*EXP(-3.0*YSVD)*1.0E13
	return
	end
!==============================================================
C	stbrn[a.u.]
C	Fprb, Probability for 1.008 MeV T NBI to burn out on D target
C	approximation of cold Deuterium (Ti = 0) DT (50:50 plasmas)
C	crossection by S.V.Putvinskij//VANT,v.2,1988,p.3
C	Corrected by S.V.Putvinskij 22.11.89
C	Use:		Fprb = stbrn*Nd	[1/m#3]
C  		 D-T	Fusion power thermal D  + 1.008 MeV T fast
C
C				Polevoy		06-JUL-99,08-JUN-07
c	ABEAM = 3.
c	EBEAM = 1008.

	double precision function STBRNy(TE,NE,Z2NdA)
	implicit none
	integer j,jk
	double precision TE,NE,Z2NdA
	double precision YECM,YSQ,YECDEB,YXC3,YX,YX2,YX3,YS2,YE,YSIG
		YECM	=1008.d0*.4d0
		YSQ	=DSQRT(YECM)
		stbrny	=0.d0
		YECDEB  = 14.6d0*TE*3.d0/1008.d0*(Z2NdA/NE)**0.667
		YXC3 =YECDEB*DSQRT(YECDEB)
	DO 	Jk = 1,1000
		YX	=(1.d-3*Jk)
		YX2	=YX*YX
		YX3	=YX2*YX
		YE	=YECM*YX2
		YS2	=(YE-48.7878d0)**2
		YSIG=
     >		DEXP(-34.3812d0/YX/YSQ)*(1.d0+1.1177d-5*YS2)/
     >		(1.d0+6.433d-4*YS2)/YE
		stbrny	=stbrny+ YSIG/(1.d0+YXC3/YX3)
	 enddo
	stbrny	=stbrny*26798.d0*1.d-3*
     * 	4.38d-4*DSQRT(1008.d0/3.d0)*
     *	2.d0*3.d0/(NE+1.d-9)*TE*DSQRT(TE)/
     /	(15.85d0+DLOG(TE/DSQRT(NE+1.d-9)))
		return
		end
!!=========================================================
      double precision function svcxy(TI,AMAIN)
      implicit none
      double precision TI,AMAIN
		svcxy=0.
      if (ti.gt.0. and. amain.ge.1.) then
		SVCXy	=10.**(5.9+0.3*LOG10(TI/AMAIN))
      endif
      return
      end
!!=========================================================
      double precision function SVRCy(TE)
         implicit none
         double precision TE
      IF(TE.LT..0001) THEN
         SVRCy=0.
      ELSE
		SVRCy	=13.6E-3/TE
		SVRCy	=1.27*SVRCy*sqrt(SVRCy)/(SVRCy+.59)
      ENDIF
         return
         end
!!=========================================================
      double precision function SVIEy(TE)
      implicit none
		double precision	TE
		SVIEy=0.
      IF(TE.GT..01)	THEN
		SVIEy	=.0136/TE
      IF(TE.GT..01)	THEN
		SVIEy=
     >	9.7E5*EXP(-SVIEy)*SQRT(SVIEy/(1.+SVIEy))/(SVIEy+.73)
				ELSE
		SVIEy=2.958E5*EXP(-SVIEy)*SQRT(SVIEy)
				ENDIF
				endif
      return
      end
!=========================================================
C PENLI [MW/m#3]: Electron heat losses due to cold neutrals Radiation
C	PENLI=<sigma*v>*Ne*Nn*(0.0102keV)
C	PE=...-PENLI = ...-	PENLIy*NE*f0
C			(Pereverzev 09-JULY-90)
	double precision function PENLIy(TE)
	implicit none
	double precision TE,YY
		PENLIy=0.
	if(TE.gt.0.) then
		YY	=.0102/TE
		PENLIy	=0.48/(0.28+YY)*SQRT(YY*(1.+YY))*EXP(-YY)
		PENLIy	=16.3*PENLIy
	endif
		return
		end
!========================================================
C PAIONy [MW/m#3]   D-T Fraction of fusion alpha power deposited to ions
C     P.Pavlo  22.06.89/A.Polevoi 20-MAY-94
	double precision function PAIONy(TE,EBDECTE)
	implicit none
	double precision TE,EBDECTE,y,y2
		PAIONy=1.
	if(TE.gt.0.) then
        y2 = EBDECTE/TE
        y = sqrt(y2)
        PAIONy    = 2.* (  0.166666667*LOG( (1.-y+y2)/(1.+2.*y+y2) ) +
     + 0.57735026*(ATAN( 0.57735026*(2.*y-1.) )+0.52359874) ) / y2
	endif
	return
	end

C======================================================================|
		subroutine	SMTH(ALFA,NO,FO,XO,N,FN,XN,NRD)
C----------------------------------------------------------------------|
C copy of SMTH() (see in for/surv.f), but in/out names can be the same
	implicit none
!	include	'for/parameter.inc'
	integer	NO,N,J,I,NRD
	double precision ALFA,XO(*),FO(*),XN(*),FN(*),P(NRD),YFN(NRD)
	double precision YF,YX,YP,YQ,YD,FJ
	if (N .gt. NRD .or. NO .le. 0)	then
		write(*,*)' >>> SMTH: array is out of limits'
!		call	a_stop
	endif
	if (NO .eq. 1)	then
	   do	j=1,N
		FN(j) = FO(1)
	   enddo
	   return
	endif
	if (NO .eq. 2)	then
	  do	j=1,N
	   FN(j)=(FO(2)*(XN(j)-XO(1))-FO(1)*(XN(j)-XO(2)))/(XO(2)-XO(1))
	  enddo
	  return
	endif
	if (N .lt. 2)	then
		write(*,*)' >>> SMTH: no output grid is provided'
!		call	a_stop
	endif
	if (abs(XO(NO)-XN(N)) .gt. XN(N)/N)	then
	    write(*,*)'>>> SMTH: grids are not aligned'
	    write(*,'(1A23,I4,F8.4)')'     Old grid size/edge',NO,XO(NO)
	    write(*,'(1A23,I4,F8.4)')'     New grid size/edge',N,XN(N)
!	    call	a_stop
	endif
	do	1	j=2,N
	   YP = (XN(j)-XN(j-1))
	   if (YP .le. 0.d0)	then
	write(*,*)'>>> SMTH: new grid is not increasing monotonically'
	      write(*,'(A,I4,A,F8.4)')'Node ',j-1,'   Value',XN(j-1)
	      write(*,'(A,I4,A,F8.4)')'Node ',j,  '   Value',XN(j)
!	      call	a_stop
	   endif
	   P(j)	=ALFA/YP/XO(NO)**2
 1	continue
	P(1)	=0.
	YFN(1)	=0.
	I	=1
	YF	=(FO(2)-FO(1))/(XO(2)-XO(1))
	YX	=2./(XN(2)+XN(1))
	YP	=0.
	YQ	=0.
	do	5	j=1,N-1
		if(XO(I) .gt. XN(j))	GO TO 4
 3		I	=I+1
		if(I .gt. NO)	I=NO
		if(I .ne. NO .and. XO(I) .lt. XN(j))	GOTO	3
		YF	=(FO(I)-FO(I-1))/(XO(I)-XO(I-1))
 4		FJ	=FO(I)+YF*(XN(j)-XO(I))
		YD=1.+YX*(YP+P(j+1))
		P(j)	=YX*P(j+1)/YD
		YFN(j)	=(FJ+YX*YQ)/YD
		if (j .eq. N-1)	goto	5
		YX	=2./(XN(j+2)-XN(j))
		YP	=(1.-P(j))*P(j+1)
		YQ	=YFN(j)*P(j+1)
 5	continue
	FN(N)	=FO(NO)
	do	6	j=N-1,2,-1
		FN(j)	=P(j)*FN(j+1)+YFN(j)
 6	continue
	FN(1)	=FO(1)
	end
C======================================================================|
