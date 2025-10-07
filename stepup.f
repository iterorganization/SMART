      Subroutine STEPUPN0(
     >   NA1,TAU,HRO,VRo,VR,G11,SLAT,RHO,
     >   CN,DN,NEo,NE,NEX,QNB,SN,SNN,SNTOT,QN,GN,GNX)

! diffusive model for hydrogen neutrals QNB[10^19at/s] neutral influx > 0.
! sign is introduced in line 47
         implicit none
         integer NA,NA1,ND,ND1,NA1N,J
         double precision
     >    TAU,HRO,VRo(*),VR(*),G11(*),SLAT(*),RHo(*),
     >    CN(*),DN(*),NEo(*),NE(*),NEX(*),
     >    SN(*),SNN(*),SNTOT(*),QN(*),GN(*),GNX(*),
     >    YHRO,HROA,QNB
         double precision, allocatable ::
     >     YWA(:),YWB(:),YWC(:),DSN(:)
         allocate(
     >     YWA(NA1),YWB(NA1),YWC(NA1),DSN(NA1)
     >      )
         YWA = 0.d0
         YWB = 0.d0
         YWC = 0.d0
         DSN = 0.d0

!      callmarkloc("tmp/eqns.inc"//char(0))
!2100    continue
C **** Density equation
!      callmarkloc("NE equation"//char(0))
         YHRO = HRO
         NA=NA1-1
         HROA=RHO(na1)-RHO(NA)
         do 2212 J=1,NA1
!           DN(J)=0.d0
!           CN(J)=0.d0
!           SN(J)=SNEBM(J)
!           SNN(J)=0.d0
            DSN(J)=0.d0
            SNTOT(J)=SN(J)
            YWA(J)=DN(J)
            if (j.gt.NA) goto 2212
            if (j.eq.NA) YHRO = HROA
            YWB(J)=-CN(J)
 2212    continue
         ND1 = NA1
         NA1N = ND1
         ND = ND1-1
         YWC(1) = RHO(ND1)-RHO(ND)
         if (ND1.lt.NA1) then
            do j=ND1,NA1
               NE(J)=NEX(J)
            enddo
         endif
!      NE(ND1)=NEB
         NEO(ND1)=NE(ND1)
         YWC(2)=-QNB
         YWC(3)=0.d0
         YWC(4)=-1.d0
         call RUNNa(YWA,YWB,DSN,SNN,SN,NEO,ND,TAU,HRO,YWC,NE,VRO,VR,G11)
         do 2214 J=1,NA
            QN(J)=-G11(J)*(YWA(J)*NE(J+1)-YWB(J)*NE(J))
            QN(J)=QN(J)+SLAT(J)*GNX(J)
            GN(J)=QN(J)/SLAT(J)
            SNTOT(J)=SNTOT(J)+SNN(J)*NE(J)
 2214    continue
         QN(NA1)=QN(NA)
         GN(NA1)=QN(NA1)/SLAT(NA1)
         SNTOT(NA1)=SNTOT(NA)
         deallocate(
     >     YWA,YWB,YWC,DSN
     >      )

         return
      end
!===================================================
      Subroutine STEPUPN(
     >  NA1,TAU,HRO,VRo,VR,G11,SLAT,RHO,
     > CN,DN,DSN,NEo,NE,NEX,NEB,SN,SNN,SNTOT,QN,GN,GNX
     > )
         !  use physics_module_level1
         implicit none
         integer NA,NA1,ND,ND1,NA1N,J
         double precision
     >    TAU,HRO,VRo(*),VR(*),G11(*),SLAT(*),RHo(*),
     >    CN(*),DN(*),DSN(*),NEo(*),NE(*),NEX(*),
     >    SN(*),SNN(*),SNTOT(*),QN(*),GN(*),GNX(*),
     >    YHRO,HROA,NEB
         double precision, allocatable ::
     >     YWA(:),YWB(:),YWC(:)
         allocate(
     >     YWA(NA1),YWB(NA1),YWC(NA1)
     >      )
         YWA = 0.d0
         YWB = 0.d0
         YWC = 0.d0

!      callmarkloc("tmp/eqns.inc"//char(0))
!2100    continue
C **** Density equation
!      callmarkloc("NE equation"//char(0))
         YHRO = HRO
         NA=NA1-1
         HROA=RHO(na1)-RHO(NA)
         do 2212 J=1,NA1
!           DN(J)=0.d0
!           CN(J)=0.d0
!           SN(J)=SNEBM(J)
            SNTOT(J)=SN(J)
            YWA(J)=DN(J)
            if (j.gt.NA) goto 2212
            if (j.eq.NA) YHRO = HROA
            YWB(J)=-CN(J)
 2212    continue
         ND1 = NA1
         NA1N = ND1
         ND = ND1-1
         YWC(1) = RHO(ND1)-RHO(ND)
         if (ND1.lt.NA1) then
            do j=ND1,NA1
               NE(J)=NEX(J)
            enddo
         endif
         NE(ND1)=NEB
         NEO(ND1)=NE(ND1)
         YWC(4)=1.d0
         call RUNNa(YWA,YWB,DSN,SNN,SN,NEO,ND,TAU,HRO,YWC,NE,VRO,VR,G11)
         do 2214 J=1,NA
            QN(J)=-G11(J)*(YWA(J)*NE(J+1)-YWB(J)*NE(J))
            QN(J)=QN(J)+SLAT(J)*GNX(J)
            GN(J)=QN(J)/SLAT(J)
            SNTOT(J)=SNTOT(J)+SNN(J)*NE(J)
 2214    continue
         QN(NA1)=QN(NA)
         GN(NA1)=QN(NA1)/SLAT(NA1)
         SNTOT(NA1)=SNTOT(NA)
         deallocate(
     >     YWA,YWB,YWC
     >      )
         return
      end
!=====
            subroutine STEPUPT(
     >  NA1,NB1,TAU,HRO,VRo,VR,G11,SLAT,RHO,
     >  XI,HE,DSI,DSE,
     >  PE,PET,PETOT,PI,PIT,PITOT,Z2NdA,
     >  TEX,TE,TEo,TEB,TIX,TI,TIo,TIB,
     >  NEX,NEo,NE,NIX,NIo,NI,Qe,Qi,GNX,
     >  GN2E,GN2I)

         !  use physics_module_level1

         implicit none
         integer NA,NA1,ND,ND1,NA1E,NA1I,NB1,J
         double precision
     >    TAU,HRO,TEB,TIB,YHRO,HROA,GN2E,GN2I
         double precision VRo(*),VR(*),G11(*),SLAT(*),RHo(*),
     >    XI(*),HE(*),GNX(*),QE(*),QI(*),DSE(*),DSI(*),
     >    PE(*),PET(*),PETOT(*),
     >    PI(*),PIT(*),PITOT(*)
         double precision
     >    TEX(*),TE(*),TEo(*),TIX(*),TI(*),TIo(*),
     >    NEX(*),NEo(*),NE(*),NIX(*),NIo(*),NI(*),Z2NdA(*)
         double precision, allocatable ::
     >     YWA(:),YWB(:),YWC(:),YWD(:),
     >     PDI(:),PDE(:),WORK1(:,:)
         external RUNTTa,NURTTa
         allocate(
     >     YWA(NA1),YWB(NA1),YWC(NA1),YWD(NA1),
     >     PDI(NA1),PDE(NA1),WORK1(NA1,24)
     >      )
         YWA = 0.d0
         YWB = 0.d0
         YWC = 0.d0
         YWD = 0.d0
         PDI = 0.d0
         PDE = 0.d0
         WORK1 = 0.d0
         GN2E=0.d0   !tmp
         GN2I=0.d0   !tmp
C **** Electron temperature equation
!      callmarkloc("TE equation"//char(0))
         NA=NA1-1
         HROA=RHO(NA1)-RHO(NA)

         do 2221 J=1,NA1
!      HE(J)=0.d0
!      PE(J)=PEBM(J)
            if (j.gt.NA) goto 2221
            YWA(J)=HE(J)
            YWA(J)=YWA(J)*(NE(J+1)+NE(J))*0.5
            YWD(J)=0.
            YWD(J)=YWD(J)*(NE(J+1)+NE(J))*0.5+GN2E*GNX(J)*SLAT(J)/G11(J)
 2221    continue
!2222    YHRO = HRO
         YHRO = HRO
         do 2223 J=1,NA1
!      PET(J)=0.
            PET(J)=0.d0
            PETOT(J)=PE(J)
            if (j.gt.NA) goto 2223
            YWB(J)=-YWD(J)
            YWC(J)=PE(j)
 2223    continue
         ND1 = NA1
         NA1E = ND1
         ND = ND1-1
         QE(1) = RHO(ND1)-RHO(ND)
         if (ND1.lt.NA1) then
            do j=ND1,NA1
               TE(J)=TEX(J)
            enddo
         endif
         TE(ND1)=TEB
         TEO(ND1)=TE(ND1)
         QE(4)=1.d0
         YWD(ND1)=0.d0
         DSE(ND1)=0.d0
         call RUNTTa(YWA,YWB,PET,YWC,NEO,NE,TEO,
     >    ND,TAU,HRO,QE(1),YWD,DSE,VRO,VR,G11,WORK1,Z2NdA,TE,NE)
!     >    ND,TAU,HRO,QE(1),YWD,DSE,VRO,VR,G11,WORK1,PEI)
         do J=ND1,NB1
            PDE(j) = 0.d0
         enddo
         do J=1,ND
            PDE(j) = 0.d0
         enddo
C **** Ion temperature equation
C      callmarkloc("TI equation"//char(0))
         do 2231 J=1,NA1
!      XI(J)=0.d0
!      PI(J)=PIBM(J)
            PIT(J)=0.d0
            if (j.gt.NA) goto 2231
            YWA(J)=XI(J)
            YWA(J)=YWA(J)*(NI(J+1)+NI(J))*0.5
            YWD(J)=2.*GN2I*GNX(J)*SLAT(J)/G11(J)/(NE(J+1)+NE(J))
            YWD(J)=0.5*YWD(J)*(NI(J+1)+NI(J))
 2231    continue
!2232    YHRO = HRO
         YHRO = HRO
         do 2233 J=1,NA1
!      PIT(J)=0.
            PITOT(J)=PI(J)
            if (j.gt.NA) goto 2233
            YWB(J)=-YWD(J)
            YWC(J)=PI(j)
 2233    continue
         ND1 = NA1
         NA1I = ND1
         ND = ND1-1
         QI(1) = RHO(ND1)-RHO(ND)
         if (ND1.lt.NA1) then
            do j=ND1,NA1
               TI(J)=TIX(J)
            enddo
         endif
         TI(ND1)=TIB
         TIO(ND1)=TI(ND1)
         QI(4)=1.d0
         YWD(ND1)=0.d0
         DSI(ND1)=0.d0
         call RUNTTa(YWA,YWB,PIT,YWC,NIO,NI,TIO,
     >   ND,TAU,HRO,QI(1),YWD,DSI,VRO,VR,G11,WORK1,Z2NdA,TE,NE)
!     >   ND,TAU,HRO,QI(1),YWD,DSI,VRO,VR,G11,WORK1,PEI)
         do J=ND1,NB1
            PDI(j) = 0.d0
         enddo
         do J=1,ND
            PDI(j) = 0.d0
         enddo
         call NURTTa(TE,TI,QE,QI,PETOT,PITOT,ND,WORK1)
         if (ND1.lt.NA1) then
            do j=ND1+1,NA1
               QE(j)=QE(ND1)
               QI(j)=QI(ND1)
            enddo
         endif
         deallocate(
     >     YWA,YWB,YWC,YWD,PDE,PDI,WORK1
     >      )
         return
      end

!==================================================
      subroutine stepupf(NA1,RHO,TAU,RTOR,BTOR,IPL,
     > CUBS,CD,CC,G22,G33,IPOL,
     > FP,FPo,MU,CU,UPL,ULON,FV) ! output
      implicit none
      integer j,NA1,NA
      double precision HRO,ROC,TAU,RTOR,BTOR,IPL,HROA,
     > CUBS(*),CD(*),CC(*),G22(*),G33(*),IPOL(*),
     > FP(*),FPo(*),MU(*),CU(*),UPL(*),ULON(*),RHO(*),FV(*)
        double precision YD,YC,YYD,GP,ARRNA1
        double precision, allocatable ::
     >  YWA(:),YWB(:),YWC(:),YWD(:)
      external RUNF,CUOFPy,ARRNA1
	allocate(
     >  YWA(NA1),YWB(NA1),YWC(NA1),YWD(NA1)
     >   )
        NA=NA1-1
        HRO=RHO(3)-RHO(2)
        ROC=RHO(NA1)
        NA=NA1-1
        HROA=ROC-RHO(NA)
        GP=3.14159263359d0
      YD=-0.8*GP*GP*RTOR
      YC=.4*GP
      do 2225 J=1,NA1
            YYD=CUBS(J)+CD(J)
        YWD(J)=YYD*YD/(IPOL(J)**3*G33(J))
        YWB(J)=CC(J)*YC/IPOL(J)**2
 2225 continue

C Prescribed plasma current:
      FP(NA)=.4*GP*HROA*RTOR/(G22(NA)*IPOL(NA1))
      FP(NA1)=(HROA-.5*HRO)*ROC*CC(NA1)/RTOR/IPOL(NA1)/TAU
      FP(NA1)=0.
      FP(NA-1)=1.+FP(NA)*FP(NA1)
      FP(NA1)=FP(NA)*(IPL+FPO(NA1)*FP(NA1))
      FP(NA)=-1.
      call RUNF(G22,YWB,YWC,YWD,FPO,NA,TAU,HRO,HROA,FP,FV)
!                AK,B,C,D,FO,N,GT,H,HB,F,FV
      do 2226 J=1,NA1
      UPL(J)=YWD(J)
      ULON(J)=IPOL(J)*G33(J)*UPL(J)
 2226 continue
      UPL(NA1)=ARRNA1(UPL(NA-2),HROA/HRO)
      ULON(NA1)=IPOL(NA1)*G33(NA1)*UPL(NA1)
      call CUOFPy(NA1,RHO,RTOR,BTOR,FP,MU,CU,G22,G33,IPOL)
      deallocate(
     >  YWA,YWB,YWC,YWD
     >   )
      return
      end
!===================================================
