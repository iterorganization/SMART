      Subroutine STEPUPN0(
     >   NA1,NB1,TAU,HRO,VRo,VR,G11,SLAT,RHO,
     >   CN,DN,NEo,NE,NEX,QNB,SN,SNN,SNTOT,QN,GN,GNX)

! diffusive model for hydrogen neutrals QNB[10^19at/s] neutral influx > 0.
! sign is introduced in line 47
         implicit none
         integer NA,NA1,ND,ND1,NA1N,NA1E,NB1,J,JIT,JEX,Jcall
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
!      call	markloc("tmp/eqns.inc"//char(0))
 2100    continue
C **** Density equation
!      call	markloc("NE equation"//char(0))
         YHRO = HRO
         NA=NA1-1
         do 2212 J=1,NA1
!           DN(J)=0.d0
!           CN(J)=0.d0
!           SN(J)=SNEBM(J)
            SNN(J)=0.d0
            DSN(J)=0.d0
            GNX(J)=0.d0
            SNTOT(J)=SN(J)
            YWA(J)=+DN(J)
            if (j.gt.NA)	goto 2212
            if (j.eq.NA)	YHRO = HROA
            YWB(J)=0.
     >       -CN(J)
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
         YWC(3)=0.
         YWC(4)=-1.
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
     >  NA1,NB1,TAU,HRO,VRo,VR,G11,SLAT,RHO,
     > CN,DN,NEo,NE,NEX,NEB,SN,SNN,SNTOT,QN,GN,GNX
     > )
         !  use physics_module_level1
         implicit none
         integer NA,NA1,ND,ND1,NA1N,NA1E,NB1,J,JIT,JEX,Jcall
         double precision
     >    TAU,HRO,VRo(*),VR(*),G11(*),SLAT(*),RHo(*),
     >    CN(*),DN(*),NEo(*),NE(*),NEX(*),
     >    SN(*),SNN(*),SNTOT(*),QN(*),GN(*),GNX(*),
     >    YHRO,HROA,NEB
         double precision, allocatable ::
     >     YWA(:),YWB(:),YWC(:),DSN(:)
         allocate(
     >     YWA(NA1),YWB(NA1),YWC(NA1),DSN(NA1)
     >      )
!      call	markloc("tmp/eqns.inc"//char(0))
 2100    continue
C **** Density equation
!      call	markloc("NE equation"//char(0))
         YHRO = HRO
         NA=NA1-1
         HROA=RHO(na1)-RHO(NA)
         do 2212 J=1,NA1
!           DN(J)=0.d0
!           CN(J)=0.d0
!           SN(J)=SNEBM(J)
            SNN(J)=0.d0
            DSN(J)=0.d0
            YWA(J)=0.d0
            GNX(J)=0.d0
            SNTOT(J)=SN(J)
            YWA(J)=+DN(J)
            if (j.gt.NA)	goto 2212
            if (j.eq.NA)	YHRO = HROA
            YWB(J)=0.
     >       -CN(J)
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
         YWC(4)=1.
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
      subroutine STEPUPT(
     > NA1,NB1,TAU,HRO,VRo,VR,G11,SLAT,RHO,
     > XI,HE,
     > PE,PET,PETOT,PI,PIT,PITOT,PEI,
     > TEX,TE,TEo,TEB,TIX,TI,TIo,TIB,
     > NEX,NEo,NE,NIX,NIo,NI,Qe,Qi,GNX,
     > GN2E,GN2I
     > )
         !  use physics_module_level1

         implicit none
         integer NA,NA1,ND,ND1,NA1N,NA1E,NA1I,NB1,J,JIT,JEX,Jcall
         double precision
     >    TAU,HRO,VRo(*),VR(*),G11(*),SLAT(*),RHo(*),
     >    XI(*),HE(*),GNX(*),QE(*),QI(*),
     >    PE(*),PET(*),PETOT(*),
     >    PI(*),PIT(*),PITOT(*),PEI(*),
     >    TEX(*),TE(*),TEo(*),TEB,TIX(*),TI(*),TIo(*),TIB,
     >    NEo(*),NE(*),NEX(*),NIo(*),NI(*),NIX(*),
     >    YHRO,HROA,GN2E,GN2I
         double precision, allocatable ::
     >     YWA(:),YWB(:),YWC(:),YWD(:),
     >     DSE(:),DSI(:),PDI(:),PDE(:),WORK1(:,:)
         allocate(
     >     YWA(NA1),YWB(NA1),YWC(NA1),YWD(NA1),
     >     DSE(NA1),DSI(NA1),PDI(NA1),PDE(NA1),WORK1(NA1,24)
     >      )

C **** Electron temperature equation
!      call	markloc("TE equation"//char(0))
         NA=NA1-1
         HROA=RHO(NA1)-RHO(NA)

         do 2221 J=1,NA1
!      HE(J)=0.d0
!      PE(J)=PEBM(J)
            if (j.gt.NA)	goto 2221
            YWA(J)=0.
     >       +HE(J)
            YWA(J)=YWA(J)*(NE(J+1)+NE(J))*0.5
            YWD(J)=0.
            YWD(J)=YWD(J)*(NE(J+1)+NE(J))*0.5+GN2E*GNX(J)*SLAT(J)/G11(J)
 2221    continue
 2222    YHRO = HRO
         do 2223 J=1,NA1
!      PET(J)=0.
            PET(J)=0.
            PETOT(J)=PE(J)
            if (j.gt.NA)	goto 2223
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
         QE(4)=1.
         YWD(ND1)=0.
         DSE(ND1)=0.
         call RUNTTa(YWA,YWB,PET,YWC,NEO,NE,TEO,
     >    ND,TAU,HRO,QE(1),YWD,DSE,VRO,VR,G11,WORK1,PEI)
         do	J=ND1,NB1
            PDE(j) = 0.
         enddo
         do	J=1,ND
            PDE(j) = 0.
         enddo
C **** Ion temperature equation
C      call	markloc("TI equation"//char(0))
         do 2231 J=1,NA1
!      XI(J)=0.d0
!      PI(J)=PIBM(J)
            PIT(J)=0.
            if (j.gt.NA)	goto 2231
            YWA(J)=0.
     >       +XI(J)
            YWA(J)=YWA(J)*(NI(J+1)+NI(J))*0.5
            YWD(J)=2.*GN2I*GNX(J)*SLAT(J)/G11(J)/(NE(J+1)+NE(J))
            YWD(J)=0.5*YWD(J)*(NI(J+1)+NI(J))
 2231    continue
 2232    YHRO = HRO
         do 2233 J=1,NA1
!      PIT(J)=0.
            PITOT(J)=PI(J)
            if (j.gt.NA)	goto 2233
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
         QI(4)=1.
         YWD(ND1)=0.
         DSI(ND1)=0.
         call RUNTTa(YWA,YWB,PIT,YWC,NIO,NI,TIO,
     >   ND,TAU,HRO,QI(1),YWD,DSI,VRO,VR,G11,WORK1,PEI)
         do	J=ND1,NB1
            PDI(j) = 0.
         enddo
         do	J=1,ND
            PDI(j) = 0.
         enddo
         call NURTTa(TE,TI,QE,QI,PETOT,PITOT,ND,WORK1)
         if (ND1.lt.NA1) then
            do j=ND1+1,NA1
               QE(j)=QE(ND1)
               QI(j)=QI(ND1)
            enddo
         endif
         deallocate(
     >     YWA,YWB,YWC,YWD,DSE,DSI,PDE,PDI,WORK1
     >      )
         return
      end
C **** Current profile adjustment
