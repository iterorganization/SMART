C======================================================================|
      subroutine RUNNA(A,B,DS,C,D,NO,N,GT,H,E,NN,VO,VN,G11)
C     call       RUNNa(YWA,YWB,DSN,SNN,SN,NEO,ND,TAU,HRO,YWC,NE,VRO,VR,G11)
C----------------------------------------------------------------------|
C Example call:
C      SNN(J)=SNNEU
C      YWA(J)=DN(J)
C      YWB(J)=-CN(J)   ! CN=-VP*VRHH
C      NE(ND1)=...
C      NEO(ND1)=NE(ND1)
C      YWC(4)=1.
C      call RUNN(YWA,YWB,SNN,SN,YWD,NEO,ND,TAU,HRO,YWC,NE,VRO,VR,G11)
C----------------------------------------------------------------------|
C Exponential scheme
C----------------------------------------------------------------------|
C       The subroutine provides run of the equation
C       dn/dt=1/V'*d[V'<(\nabla\rho)^2>*(A*dn/dr+B*n)]/dr+C*n+D
C       Scheme: h_j/GT*V'(j)*(NN(j)-NO(j))=
C          =GT*G11(j+1/2)*A(j+1/2)/h(j+1/2)
C               *(NN(j+1)*f(j+1/2)-NN(j)*g(j+1/2))-
C          -GT*G11(j-1/2)*A(j-1/2)/h(j-1/2)
C               *(NN(j)*f(j-1/2)-NN(j-1)*g(j-1/2))-
C          +h_j*(V'(j)*C(j)*NN(j)+D(j))
C       Here    A,B,C,D         - arrays(1:N)
C               NOld,NNew       - arrays(1:N+1)
C               VOld,VNew       - arrays(1:N+1)
C       Input:  H       - radial step (m)
C               HB      - edge grid cell size (m)
C               N+1     - number of mesh points
C               GT      - time step (sec)
C               A(N),B(N),C(N),D(N),VO(N),VN(N),NO(N),NN(N+1),E(1:4)
C               E(1) = HROA
C               E(4) < 0 if (NEB isn't set) .and. (QNB .or. QNNB is set)
C                       then E(2)=QNB or E(3)=QNNB
C       Not used: SN(N+1)
C                 SNN(N+1)
C       Output:        NN(N) - new quantity
C               A(N)  - coefficient at (j+1) for flux calculation
C               B(N)  - coefficient at  (j)  for flux calculation
C       Internal use:
C               E(*)
C----------------------------------------------------------------------|
         implicit none
         integer N,j
         double precision
     1   A(*),B(*),DS(*),C(*),D(*),NO(*),NN(*),VO(*),VN(*),E(*),
     2   G11(*),H,HB,Q2,Q3,BC,HJ,GT,G0,G1,P0,P1,Q0,Q1,
     3   Y0,Y1,YA,YB,YS,YJ,AJ,BJ,CJ,DJ
         HB = E(1)
         Q2 = E(2)
         Q3 = E(3)
         BC = E(4)
         G1 = 0.d0
         P1 = 0.d0
         Q1 = 0.d0
         Y1 = 0.d0
         YS = 0.d0
         HJ = H
         do 10 J=1,N
            G0 = G1
            P0 = P1
            Q0 = Q1
            G1 = GT*G11(j)
            YA = A(j)
            YB = B(j)
            if (YA .lt. 0.d0) then
               Y1 = j*H
               goto 99
            endif
C Power-law scheme: (1 line)
            P1 = 0.5*(dabs(YB)+YB)               ! 0.5(|B|+B)
            if (YA .eq. 0.d0) goto 1
            if (j .eq. N) HJ = HB
C Exponential scheme: (7 lines)
            YJ = HJ*YB/YA                       ! |\xi|
            if (dabs(YJ) .ge. 4.d1) goto 1       ! Use (A/h)*f(\xi) = .5*(|B|+B)
            if (dabs(YJ) .ge. 1.d-5) then
               P1 = YB/(1.-dexp(-YJ))            ! Use (A/h)*f(\xi) = B/(1-exp{-\xi})
            else
               P1 = YA/HJ*(1.+0.5*YJ)           ! Use (A/h)*f(\xi) = A/h/(1-\xi/2)
            endif
C Power-law scheme: (6 lines)
C          YJ = abs(HJ*YB/YA)                   ! \xi = h B/A
C          if (YJ .ge. 1.d1) goto 1
C          Y1 = 1.d0-1.d-1*YJ                   ! (1 - 0.1|\xi|)
C          Y2 = Y1*Y1
C          YJ = Y2*Y2*Y1*YA/HJ                  ! (A/h)*(1 - 0.1|\xi|)^5
C          P1 = P1+YJ                           ! (A/h)*f(\xi) is done
    1       continue
            Q1 = P1-YB
            Y0 = Y1
            Y1 = DS(j)/HJ
            AJ = G1*(P1+Y1)
            BJ = G1*(Q1+Y1)+G0*(P0+Y0)+H*(VN(J)-C(J)*VN(J)*GT)
            CJ = G0*(Q0+Y0)
            YJ = YS
            YS = G1*Y1*(NO(j+1)-NO(j))
            DJ = H*(D(J)*VN(J)*GT+NO(J)*VO(J))-YS+YJ
            if(J .ne. 1) then
               BJ = BJ-CJ*E(J-1)
               DJ = DJ+CJ*NN(J-1)
            endif
            E(J) = AJ/BJ
            NN(J) = DJ/BJ
            A(j) = P1
            B(j) = Q1
   10    continue
         if (BC .lt. 0.) NN(N+1) = (G11(N)*Q1*NN(N)-Q2)
     .   /(G11(N)*(P1-E(N)*Q1)+Q3)
         do J=N,1,-1
            NN(J) = E(J)*NN(J+1)+NN(J)
         enddo
         return
   99    write(*,'(2A,F10.6,A,I4,A,F10.6)')
     &   " >>> ERROR >>> Diffusion coefficient shouldn't be negative.",
     &   "  RHO =",Y1,"  node =",j,"   D =",YA
!        call IFKEY(ichar(' '))
      end
C======================================================================|
      subroutine RUNTTA(
     >  A,B,C,D,NO,NN,TO,N,GT,H,HB,DV,DS,VO,VN,G11,W,Z2NdA,TE,NE
     >  )
C----------------------------------------------------------------------|
C     call RUNTT (YWA,YWB,PET,YWC,NEO,NE,TEO,
C                 ND,TAU,HRO,QE(1),YWD,DSE,VRO,VR,G11,WORK1,PEI)
C----------------------------------------------------------------------|
C Exponential scheme
C----------------------------------------------------------------------|
C       The subroutine makes time step in the matrix equation:
C               d(N*T)/dt=1/V'*d[V'*(A*dT/dr+B*T)]/dr+625.*(C*T+D)
C
C       Here    A,B,C,D   - arrays(1:N)
C               NOld,NNew - arrays(1:N+1)       (densities)
C               VOld,VNew - arrays(1:N+1)       (dV/drho)
C               TO      - array (1:N+1) (temperature)
C               W(*)    - work space (e.g., work1(*))
C       Input:  H       - radial step (m)
C               N       - number of mesh points
C               GT      - time step (sec)
C               A(1:N)  - diffusivity n_e*\chi_e (or n_i*\chi_i)
C               B(1:N)  - convective velocity e.g. 5/2*GNX(J)*SLAT(J)/G11(J)
C               C(1:N)  - PET or PIT without equipartition
C               D(1:N)  - PE  or PI
C               NO(1:N) - old density (previous time step)
C               NN(1:N) - new density (next time step)
C               VO(1:N) - old V'
C               VN(1:N) - new V'
C               DV(1:N) - Aux. heat conductiviy compensated by advection
C               DS(1:N) - Aux. heat conductiviy compensated by source
C               DV(N+1) - Enable DV treatment if DV(N+1) is nonzero
C               DS(N+1) - Enable DS treatment if DS(N+1) is nonzero
C               TO(1:N+1) - old T_e (or T_i)
C               HB      - edge cell size
C       Output: DV(1:N) - contribution to rhs due to DV
C----------------------------------------------------------------------|
         implicit none
         double precision
     1   A(*),B(*),C(*),D(*),NO(*),NN(*),TO(*),DV(*),DS(*),VO(*),
     2   VN(*),G11(*),H,HB,GT,GT23,Y625,AJ,BJ,CJ,DJ,W1,P0,P1,
     3   Q0,Q1,G0,G1,H1,HJ,YA,YB,Y0,Y1,Y2,YS,YJ,Y11,Y12,Y21,Y22
!        integer N,j,j1,jn,icall,N0,N1
         integer N,j,N0,N1
!       double precision GETPEI,W(N,*)
!       external GETPEI
         double precision Z2NdA(*),TE(*),NE(*),W(N,*)

!        save icall,N0
!        data icall/0/
         integer, save :: icall = 0
         save N0
!        call add2loc("Subroutine RUNTT"//char(0))
         GT23 = GT/1.5
         Y625 = 625.*GT23
         HJ = H
         G1 = 0.d0
         P1 = 0.d0
         Q1 = 0.d0
         H1 = H
         CJ = 0.d0
         W1 = 0.d0
         Y1 = 0.d0
         YS = 0.d0
         N1 = N+1
         do 2 J=1,N
            G0 = G1
            P0 = P1
            Q0 = Q1
            Y0 = Y1
            if (j .eq. N) H1 = HB
            G1 = GT23*G11(j)
            YA = A(j)
            YB = B(j)
            if (DV(N1) .gt. 0.) then
               Y11 = DV(j)*0.5*(NN(j)+NN(j+1))
               Y12 = Y11*dlog(TO(j)/TO(j+1))/H1
               YA = YA+Y11
               YB = YB+Y12
               Y2 = TO(j)-TO(j+1)
               if (abs(Y2) .lt. 1.d-6) then
                  DV(j) = Y11*2./(TO(j)+TO(j+1))
               else
                  DV(j) = Y12/Y2       ! Contribution to the source
               endif
            endif
            if (YA .lt. 0.d0) write(*,*) 'j,a(j)',j,a(j)
            if (YA .lt. 0.d0) goto 99
            P1 = 0.5*(dabs(YB)+YB)          ! 0.5(|B|+B) (Used if vh/D is big)
            if (YA .eq. 0.d0) goto 1
            YJ = H1*YB/YA                  ! |\xi| Peclet number
            if (dabs(YJ) .ge. 4.d1) goto 1  ! Use (A/h)*f(\xi) = .5*(|B|+B)
            if (dabs(YJ) .ge. 1.d-5) then
               P1 = YB/(1.-dexp(-YJ))       ! Use (A/h)*f(\xi) = B/(1-exp{-\xi})
            else
               P1 = YA/H1*(1.+0.5*YJ)      ! Use (A/h)*f(\xi) = A/h/(1-\xi/2)
            endif
    1       continue
            Q1 = P1-YB                     ! (A/h)*g(\xi) = (A/h)*f(\xi)-B
            if (DS(N1) .gt. 5.d-1) Y1 = DS(j)*0.5*(NN(j)+NN(j+1))/H1
            if ( icall .eq. 0) then
               N0 = N
               W(j,17) = P1*G11(j)
               W(j,18) = Q1*G11(j)
!              W(j,10) = -GETPEI(J)
               W(j,10) = -0.00246*(15.9 - .5*dlog(NE(j)) + dlog(TE(j)))
     >          *NE(j)*Z2NdA(j)/TE(j)/dsqrt(TE(j))
            else
               if (N0.ne.N) goto 98
               W(j,19) = P1*G11(j)
               W(j,20) = Q1*G11(j)
            endif
C          AJ = G1*P1
            AJ = G1*(P1+Y1)
C          BJ = G1*Q1+G0*P0
            BJ = G1*(Q1+Y1)+G0*(P0+Y0)
C          CJ = G0*Q0
            CJ = G0*(Q0+Y0)
            BJ = BJ+HJ*VN(J)*(NN(J)-Y625*(W(j,10)+C(J)))
            YJ = (VO(J)/VN(J))**0.666667
            DJ = HJ*(Y625*D(J)*VN(J)+NO(J)*TO(J)*VO(J)*YJ)
            YJ = YS
            YS = G1*Y1*(TO(j+1)-TO(j))
            DJ = DJ-YS+YJ
            W(j,1+icall) = AJ   ! A_e or A_i | P_k
            W(j,3+icall) = BJ   ! B_e or B_i | Q_k+P_{k-1}+C_k
            W(j,5+icall) = CJ   ! C_e or C_i | Q_{k-1}
            W(j,7+icall) = DJ   ! D_e or D_i |
            W(j,23+icall) = C(j)
            W(j,21+icall) = D(j)+B(j)/(Y625*HJ*VN(j))   ! PDE or PDI
    2    continue
         if (icall .eq. 0) then
            icall = 1
            return
         else
            icall = 0
         endif
C  2nd call:
C     Work array usage:         W(1:N,1<->24)
C        Exchange with RUNTT:   W(1:N,1<->8)
C        Exchange with NURTTa:  W(1:N,10<->24)
C  W(1:N,1) - A_e       W(1:N,2) - A_i
C  W(1:N,3) - B_e       W(1:N,4) - B_i
C  W(1:N,5) - C_e       W(1:N,6) - C_i
C  W(1:N,7) - D_e       W(1:N,8) - D_i
C  W(1:N,9) - not used
C  W(1:N,10)  - Pe->i
C  W(1:N,11) - E_11     W(1:N,12) - E_12
C  W(1:N,13) - E_21     W(1:N,14) - E_22
C  W(1:N,15) - G_1      W(1:N,16) - G_2
C  W(1:N,17),   W(1:N,18) - Coefficients for flux evaluation
C  W(1:N,19),   W(1:N,20) -   (used in NURTTa to compute QE, QI)
C  W(1:N,21),   W(1:N,22) - Coefficients for RHS evaluation
C  W(1:N,23),   W(1:N,24) -   (used in NURTTa to compute PETOT, PITOT)
         do 3 j=1,N
            Y1 = W(j,3)
            Y2 = W(j,4)
            YA = HJ*VN(j)*Y625*W(j,10)
            YB = YA
            if (j .gt. 1) then
               Y1 = Y1-W(j,5)*W(j-1,11) ! Direct = (B_k - Q_{k-1}*E_{k-1})
               Y2 = Y2-W(j,6)*W(j-1,14) !          (Y1 YA)
               YA = YA-W(j,5)*W(j-1,12) !          (YB Y2)
               YB = YB-W(j,6)*W(j-1,13)
            endif
            YJ = 1./(Y1*Y2-YA*YB)
            Y11 = Y2*YJ                 !           ( Y2 -YA)      (Y11 Y12)
            Y12 =-YA*YJ                 ! Inversed =         /det =
            Y21 =-YB*YJ                 !           (-YB  Y1)      (Y21 Y22)
            Y22 = Y1*YJ
            W(j,11) = Y11*W(j,1)        !
            W(j,12) = Y12*W(j,2)        !            (-YA  Y1) / det
            W(j,13) = Y21*W(j,1)        !      (W) = ( 11 12 )
            W(j,14) = Y22*W(j,2)        !          = ( 13 22 )
            Y1 = W(j,7)
            Y2 = W(j,8)
            if (j .gt. 1) then
               Y1 = Y1+W(j,5)*W(j-1,15)
               Y2 = Y2+W(j,6)*W(j-1,16)
            endif
            W(j,15) = Y1*Y11+Y2*Y12
            W(j,16) = Y1*Y21+Y2*Y22
    3    continue
         return
   98    write(*,'(2A)')
     &   ">>> ERROR >>> The same boundary is required",
     &        " for both TE and TI"
!       call IFKEY(ichar(' '))
   99    write(*,'(2A,1F10.6)')
     &   " >>> ERROR >>> Heat conductivity shouldn't be negative.",
     &   "  RHO =",j*H
!       call IFKEY(ichar(' '))
      end
C======================================================================|
      subroutine NURTTa(T1,T2,Q1,Q2,P1,P2,N,W)
C----------------------------------------------------------------------|
C TE,TI or TI,TE in the same order as by calling RUNTT
C
C In:  Boundary conditions -
C        Q1(4) > 0   ->   given T1(N+1)
C        Q2(4) > 0   ->   given T2(N+1)
C        Q1(4) < 0   ->   given Q1(N+1) = Q1(2)+T{e,i}(N+1)*Q1(3)
C        Q2(4) < 0   ->   given Q2(N+1) = Q2(2)+T{e,i}(N+1)*Q2(3)
C
C Out: Fluxes
C        Q_j(1:N+1) = -0.0016*(p_j*T_{j,k+1}-q_j*T_{j,N})
C      RHSs
C        P_j(1:N+1) =
C----------------------------------------------------------------------|
         implicit none
         double precision T1(*),T2(*),Q1(*),Q2(*),P1(*),P2(*)
         integer N,j
         double precision W(N,*)
         double precision Y1,Y2,Y11,Y12,Y21,Y22,YD
         if (Q1(4).lt.0. .and. Q2(4).lt.0.) then
C Both eqns use fluxes as boundary conditions:
            Y11 = W(N,19)*W(N,11)-W(N,17)-625.*Q1(3)
            Y22 = W(N,20)*W(N,14)-W(N,18)-625.*Q2(3)
            Y12 = W(N,19)*W(N,12)
            Y21 = W(N,20)*W(N,13)
            YD = Y11*Y22-Y12*Y21
            Y1 = 625.*Q1(2)-W(N,19)*W(N,15)
            Y2 = 625.*Q2(2)-W(N,20)*W(N,16)
            T1(N+1) = (Y11*Y1-Y12*Y2)/YD
            T2(N+1) = (Y22*Y2-Y21*Y1)/YD
         elseif (Q1(4) .lt. 0.) then
C Mixed boundary conditions: 1st eqn flux, 2nd eqn temperature
            Y11 = W(N,19)*W(N,11)-W(N,17)-625.*Q1(3)
            Y1  = 625.*Q1(2)-W(N,19)*(W(N,15)+W(N,12)*T2(N+1))
            T1(N+1) = Y1/Y11
         elseif (Q2(4) .lt. 0.) then
C Mixed boundary conditions: 2nd eqn flux, 1st eqn temperature
            Y22 = W(N,20)*W(N,14)-W(N,18)-625.*Q2(3)
            Y2  = 625.*Q2(2)-W(N,20)*(W(N,16)+W(N,13)*T1(N+1))
            T2(N+1) = Y2/Y22
         endif
C Define new temperatures:
         do J=N,1,-1
            T1(J) = W(j,11)*T1(J+1)+W(j,12)*T2(J+1)+W(j,15)
            T2(J) = W(j,13)*T1(J+1)+W(j,14)*T2(J+1)+W(j,16)
         enddo
C Define fluxes:
         do J=1,N
            Q1(j) = -0.0016*(W(j,17)*T1(j+1)-W(j,18)*T1(j))
            Q2(j) = -0.0016*(W(j,19)*T2(j+1)-W(j,20)*T2(j))
         enddo
         Q1(N+1) = Q1(N)
         Q2(N+1) = Q2(N)
C Define RHSs:
         do j=1,N
            P1(j) = W(j,21)+W(j,23)*T1(j)+W(j,10)*(T1(j)-T2(j))
            P2(j) = W(j,22)+W(j,24)*T2(j)+W(j,10)*(T2(j)-T1(j))
         enddo
         P1(N+1) = P1(N)
         P2(N+1) = P2(N)
         do J=1,N+1
            T1(j) = dmax1(T1(j),1.d-4)
            T2(j) = dmax1(T2(j),1.d-4)
         enddo
      end
!==================================
C======================================================================|
        subroutine
     >          RUNF(AK,B,C,D,FO,N,GT,H,HB,F,FV)
C----------------------------------------------------------------------|
C       The subroutine provides inversion of the matrix equation for FP
C       The boundary condition is supplied in the following form:
C               F(N-1)*Psi(N+1)+F(N)*Psi(N)=F(N+1)
C           where F(N-1), F(N) and F(N+1) are input parameters.
C       Scheme:
C       Input:  H       - radial step (m)
C               HB      - edge radial step (m)
C               N+1     - number of grid points
C               GT      - time step (sec)
C               AK(N)   - G22
C               B(N)    - conductivity
C               D(N)    - external (+bootstrap) current
C               FO(N)   - old poloidal flux
C               F(N-1),F(N),F(N+1)      - edge conditions
C       Output: F(1:N+1)- new poloidal flux
C               C(1:N+1)- (1/rho)d{K*dF/d(rho)}/d(rho) ~ current density
C               B(1:N+1)- (1/rho)dF/d(rho)      ~ rotational transform
C               D(1:N+1)- dF/dt toroidal loop voltage
C----------------------------------------------------------------------|
        implicit none
        integer N,j
        double precision
     1          AK(N+1),B(N+1),C(N+1),D(N+1),FV(N+1),F(N+1),FO(N+1),
     2          H,HB,GT,HH,AJ,BJ,CJ,DJ,RJ,RJHH
        HH = H*H
        AJ = 0.
        RJ = -0.5*H
        do      1       J=1,N
           CJ = AJ
           RJ = RJ+H
           RJHH = RJ*HH
           if (j .eq. N)        then
              CJ = CJ*HB/H
              RJHH = RJ*HB*H
           endif
           AJ = AK(J)
           DJ = RJHH*B(J)/GT
           BJ = AJ+CJ+DJ
           DJ = DJ*FO(J)+RJHH*D(J)-AJ*(FV(J+1)-FV(J))
           if(J .ne. 1) then
              BJ = BJ-CJ*C(J-1)
              DJ = DJ+CJ*(D(J-1)-FV(J-1)+FV(J))
           endif
           C(J) = AJ/BJ
           D(J) = DJ/BJ
 1      continue
        F(N+1) = (F(N+1)-F(N)*D(N))/(F(N-1)+F(N)*C(N))
        do      2       J=N,1,-1
        F(J) = C(J)*F(J+1)+D(J)
 2      continue
        do      j=1,N+1
           D(j) = (F(j)-FO(j))/GT
C          FO(j) = F(j)
        enddo
        end
C======================================================================|

