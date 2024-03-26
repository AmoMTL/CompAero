C     MAIN PROGRAM
      COMMON /BLC0/ RL,XCTR
      COMMON /BLC7/ ETA(201),DETA(201),A(201)
      COMMON /BLCS/ DLS(200),VW(200),CF(200),THT(200)
      COMMON /BLC2/ NX,NXT,NP,NPT,NTR,NVRS,IT,INVRS,ISWPT
      COMMON /BLC6/ DELF(201),DELU(201),DELV(201),DELW(201)
      COMMON /BLC8/ F(201,2),U(201,2),V(201,2),W(201,2),B(201,2)
      COMMON /BLC3/ X(201),UE(200),UE0(200),P1(200),P2(200),GMTR(200)
      COMMON /BLCG/ XC(200), YC(200)
	CHARACTER*80 input_name, output_name
	DIMENSION RS(200), RTHETA(200)
C
	WRITE(6,*) "Enter input file name (include extension name)"
	READ(5,*) input_name
	OPEN(unit=55,file=input_name,STATUS="OLD")
	WRITE(6,*) "Enter output file name"
	READ(5,*) output_name
	OPEN(unit=66,file=output_name)                                                  

      CALL INPUT
      CALL IVPL
C
      NX     = 0
 100  NX     = NX+1
      WRITE (66,3000) NX,X(NX)
      IT     = 0
      IGROW  = 0
 120  IT     = IT+1
      IF (IT.LE.5) GO TO 140
      WRITE (66,1000) P2(NX)
      CALL OUTPUT
      GO TO 160
 140  IF (NX .GE. NTR) CALL EDDY
C
      CALL COEF3
      CALL SOLV3
      WRITE (66,4000) IT,V(1,2),DELV(1)
C
      IF ( V(1,2) .LE. 0.0) THEN
C  FLOW SEPARATION
      WRITE(66,1100)
      GO TO 160
      ENDIF
C  CHECK CONVERGENCE
C  LAMINAR FLOW
      IF (NX .LT. NTR) THEN
         IF(ABS(DELV(1)) .GT. 0.0001) GO TO 120
      ELSE
C  TURBULENT FLOW
         IF(ABS(DELV(1)/V(1,2)) .GT. 0.01) GO TO 120
      ENDIF
C  CHECK BOUNDARY LAYER GROWTH
C
      IF (NP.LT.NPT) THEN
         IF (ABS(V(NP,2)) .GT. 0.0005 ) THEN
            IF (IGROW .LE. 3) THEN
               IGROW  = IGROW+1
               CALL GROWTH(1)
               IT     = 0
               GO TO 120
            ENDIF
         ENDIF
      ENDIF
      CALL OUTPUT
C
      IF (NX .LT. NXT) GO TO 100
  160 CONTINUE
C
      DO I = 1, NX
	   RS(I) = UE(I)*X(I)*RL
	   RTHETA(I)= UE(I)*THT(I)*RL
	ENDDO
C      WRITE (66,9100) (I,XC(I),X(I),
C     *                UE(I),VW(I),CF(I),DLS(I),THT(I),I=1,NX)
      WRITE (66,9200) (I,X(I),RS(I),UE(I),DLS(I),
     *                 THT(I),VW(I),CF(I),RTHETA(I),I=1,NX)
 9200 FORMAT(///2X,'*** BOUNDARY LAYER PARAMETERS  '/
     +  /2X,'I',6X,'S ',10X,'RS',
     +   10X,'UE',8X,'DELS',8X,'THETA',8X,'VW',10X,'CF',9X,'RTHETA'/
     +   (I3,8E12.4))
C
      STOP
 1000 FORMAT(//10X,'ITERATIONS EXCEEDED MAX P2 = ',E12.5)
 1100 FORMAT(//10X,'FLOW SEPARATED')
 3000 FORMAT(/2X,'NX =',I3,3X,'X=',F8.4/3X,'IT',5X,'VW',8X,'DELV'/)
 4000 FORMAT (I5,2E11.4)
 9100 FORMAT(///2X,'*** BOUNDARY LAYER PARAMETERS  '/
     +  /2X,'NX',6X,'XC',7X,'S',
     +   8X,'UE',8X,'VW',8X,'CF',8X,'DLS',8X,'THT'/
     +   (I5,2F8.4,5E11.4))
      END


      SUBROUTINE INPUT
C----------------------------------------------------
      COMMON /BLC0/ RL,XCTR
      COMMON /BLC7/ ETA(201),DETA(201),A(201)
      COMMON /BLC2/ NX,NXT,NP,NPT,NTR,NVRS,IT,INVRS,ISWPT
      COMMON /BLC3/ X(201),UE(200),UE0(200),P1(200),P2(200),GMTR(200)
      COMMON /BLCG/ XC(200), YC(200)
      DIMENSION     D1(200)

C
      GRANG1(X1,X2,X3,Y1,Y2,Y3,X0)= (2.*X0-X2-X3)/(X1-X2)/(X1-X3)*Y1
     +       +(2.*X0-X1-X3)/(X2-X1)/(X2-X3)*Y2+(2.*X0-X1-X2)
     +       /(X3-X1)/(X3-X2)*Y3
C
      READ ( 55,* ) NXT,NPT
      READ ( 55,* ) (XC(I),I=1,NXT)
      READ ( 55,* ) (YC(I),I=1,NXT)
      READ ( 55,* ) (UE(I),I=1,NXT)
      READ ( 55,* ) RL,XCTR,ETAE,VGP,DETA(1),P2(1)
C
      DO I=1,NXT
      WRITE(50,*) XC(I),YC(I),UE(I)
      END DO 
 
      WRITE(66,6000) RL,XCTR,ETAE,VGP,DETA(1),P2(1)
C  CALCULATION OF SURFACE DISTANCE
      X(1)=0.0
      DO I=2,NXT
        X(I) = X(I-1)+SQRT((XC(I)-XC(I-1))**2+(YC(I)-YC(I-1))**2)
      ENDDO
C
      IF((VGP-1.0) .GT. 0.001) THEN
      NP     = ALOG((ETAE/DETA(1))*(VGP-1.0)+1.0)/ALOG(VGP)+1.0001
      ELSE
      NP = ETAE/DETA(1)+1.0001
      END IF
C
      IF(NP .GT. NPT) THEN
      WRITE (66,9000)
      STOP
      END IF
C
      ETA(1) = 0.0
      DO 20 J=2,NPT
        DETA(J)  = VGP*DETA(J-1)
        ETA(J)   = ETA(J-1)+DETA(J-1)
        A(J)     = 0.5*DETA(J-1)
   20 CONTINUE
C
C  DETERMINE TRANSITION LOCATION
      DO 30 I=1,NXT
         IF (XC(I) .GE. XCTR) GO TO 40
   30 CONTINUE
   40 NTR = I
C
C  CALCULATE GAMTR
      PGAMTR  = 1200.
      RXNTR   = X(NTR-1)* UE(NTR-1) * RL
      GGFT    = RL**2/RXNTR**1.34*UE(NTR-1)**3
      UEINTG  = 0.0
      U1      = 0.5/UE(NTR-1)/ PGAMTR
      DO 60 I = NTR,NXT
         U2      = 0.5/UE(I)/PGAMTR
         UEINTG  = UEINTG+(U1+U2)*(X(I)-X(I-1))
         U1      = U2
         GG      = GGFT*UEINTG*(X(I)-X(NTR-1))
         IF(GG .GT. 10.0) THEN
            GMTR(I)=1.0
         ELSE
            GMTR(I)= 1.0-EXP(-GG)
         ENDIF
   60 CONTINUE
C
C  CALCULATION OF PRESSURE GRADIENT PARAMETERS
C
      DO  80 I = 2,NXT
      IF (I.LT.NXT) THEN
      D1(I) = GRANG1(X(I-1),X(I),X(I+1),UE(I-1),UE(I),UE(I+1),X(I))
      ELSE
      D1(I) = GRANG1(X(I-2),X(I-1),X(I),UE(I-2),UE(I-1),UE(I),X(I))
      ENDIF
C
      P2(I)   = X(I) * D1(I) /UE(I)
      P1(I)   = 0.5 * (1.0 + P2(I))
   80 CONTINUE
C
      P1(1)   = 0.5 * (1.0 + P2(1))
C
      RETURN
 6000 FORMAT(//5X,'RL =',E10.3,5X,'XCTR =',F7.3/'ETAE=',F8.3,3X,
     &      'VGP=',F7.3,3X,'DETA(1)=',F7.3,5X,'P2(1)=',F7.3//)
 9000 FORMAT('NP EXCEEDED NPT - PROGRAM TERMINATED')
      END


      SUBROUTINE IVPL
      COMMON /BLC7/ ETA(201),DETA(201),A(201)
      COMMON /BLC2/ NX,NXT,NP,NPT,NTR,NVRS,IT,INVRS,ISWPT
      COMMON /BLC8/ F(201,2),U(201,2),V(201,2),W(201,2),B(201,2)
C
      DO 20  J=1,NP
         ETAB = ETA(J)/ETA(NP)
         ETAB2 = ETAB**2
         F(J,2) = 0.25*ETA(NP)*ETAB2*(3.0 - 0.5*ETAB2)
         U(J,2) = 0.5*ETAB*(3.0 - ETAB2)
         V(J,2) = 1.5*(1.0 - ETAB2)/ETA(NP)
         B(J,2) = 1.0
  20  CONTINUE
      RETURN
      END

      SUBROUTINE GROWTH(INDEX)
      COMMON /BLC7/ ETA(201),DETA(201),A(201)
      COMMON /BLC2/ NX,NXT,NP,NPT,NTR,NVRS,IT,INVRS,ISWPT
      COMMON /BLC8/ F(201,2),U(201,2),V(201,2),W(201,2),B(201,2)
C
      NP1 = NP + 1
      IF (INDEX.EQ.1) THEN
      NP = MIN0(NPT,NP+2)
      NPEND = NP
      ELSE
      NPEND = NPT
      ENDIF
      DO 30 J=NP1,NPEND
      F(J,2) = F(J-1,2) + DETA(J-1)*U(J-1,2)
      U(J,2) = U(J-1,2)
      V(J,2) = 0.0
      B(J,2) = B(J-1,2)
  30  CONTINUE
      RETURN
      END


      SUBROUTINE COEF3
      COMMON /BLC7/ ETA(201),DETA(201),A(201)
      COMMON /BLC2/ NX,NXT,NP,NPT,NTR,NVRS,IT,INVRS,ISWPT
      COMMON /BLC8/ F(201,2),U(201,2),V(201,2),W(201,2),B(201,2)
      COMMON /BLC9/ S1(201),S2(201),S3(201),S4(201),S5(201),S6(201),
     +              S7(201),S8(201),R1(201),R2(201),R3(201),R4(201)
      COMMON /BLC3/ X(201),UE(200),UE0(200),P1(200),P2(200),GMTR(200)
C
      IF (NX .EQ. 1) THEN
         CEL = 0.0
      ELSE
         CEL = 0.5 * (X(NX)+X(NX-1))/(X(NX)-X(NX-1))
      ENDIF
C
      CELH= 0.5 * CEL
      P1P   = P1(NX)+CEL
      P2P   = P2(NX)+CEL
      P1PH   = 0.5 * P1P
C
      DO 100 J=2,NP
C  PRESENT STATION
         USB   = 0.5*(U(J,2)**2+U(J-1,2)**2)
         FVB   = 0.5*(F(J,2)*V(J,2)+F(J-1,2)*V(J-1,2))
         UB    = 0.5*(U(J,2)+U(J-1,2))
         VB    = 0.5*(V(J,2)+V(J-1,2))
         FB    = 0.5*(F(J,2)+F(J-1,2))
         DERBV = (B(J,2)*V(J,2)-B(J-1,2)*V(J-1,2))/DETA(J-1)
C
         IF(NX.EQ.1) THEN
            R2(J)  = - (DERBV +P1P*FVB- P2P*USB + P2(NX) )
         ELSE
C  PREVIOUS STATION
            CFB   = 0.5*(F(J,1)+F(J-1,1))
            CVB   = 0.5*(V(J,1)+V(J-1,1))
            CFVB  = 0.5*(F(J,1)*V(J,1)+F(J-1,1)*V(J-1,1))
            CUSB  = 0.5*(U(J,1)**2+U(J-1,1)**2)
            CDERBV= (B(J,1)*V(J,1)-B(J-1,1)*V(J-1,1))/DETA(J-1)
C
            CLB   = CDERBV + P1(NX-1)*CFVB + P2(NX-1) *(1.0-CUSB)
            CRB   = -CLB - P2(NX) - CEL * CUSB + CEL*CFVB
            R2(J) = CRB-(DERBV +P1P*FVB- P2P*USB + CEL*(FB*CVB-VB*CFB))
         ENDIF
C
         S1(J) = B(J,2)/DETA(J-1) + P1PH*F(J,2) - CELH*CFB
         S2(J) = -B(J-1,2)/DETA(J-1) + P1PH*F(J-1,2) - CELH*CFB
         S3(J) = P1PH*V(J,2) + CELH*CVB
         S4(J) = P1PH*V(J-1,2) + CELH*CVB
         S5(J) = -P2P*U(J,2)
         S6(J) = -P2P*U(J-1,2)
C
C  DEFINITIONS OF RJ
C
      R1(J)  = F(J-1,2)-F(J,2) + UB*DETA(J-1)
      R3(J-1)= U(J-1,2)-U(J,2) + VB*DETA(J-1)
C
  100 CONTINUE
C
      R1(1) = 0.0
      R2(1) = 0.0
      R3(NP)= 0.0
C
      RETURN
      END

      SUBROUTINE SOLV3
C
      COMMON /BLC2/ NX,NXT,NP,NPT,NTR,NVRS,IT,INVRS,ISWPT
      COMMON /BLC7/ ETA(201),DETA(201),A(201)
      COMMON /BLC8/ F(201,2),U(201,2),V(201,2),W(201,2),B(201,2)
      COMMON /BLC9/ S1(201),S2(201),S3(201),S4(201),S5(201),S6(201),
     +              S7(201),S8(201),R1(201),R2(201),R3(201),R4(201)
      COMMON /BLC6/ DELF(201),DELU(201),DELV(201),DELW(201)
      DIMENSION     A11(201),A12(201),A13(201),A14(201),
     +              A21(201),A22(201),A23(201),A24(201)
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      A11(1)= 1.0
      A12(1)= 0.0
      A13(1)= 0.0
      A21(1)= 0.0
      A22(1)= 1.0
      A23(1)= 0.0
      G11   =-1.0
      G12   =-A(2)
      G13   = 0.0
      G21   = S4(2)
      G23   =-S2(2)/A(2)
      G22   = G23 + S6(2)
      A11(2)= 1.0
      A12(2)=-A(2)-G13
      A13(2)= A(2)*G13
      A21(2)= S3(2)
      A22(2)= S5(2)-G23
      A23(2)= S1(2)+A(2)*G23
      R1(2) = R1(2)-G11*R1(1)-G12*R2(1)-G13*R3(1)
      R2(2) = R2(2)-G21*R1(1)-G22*R2(1)-G23*R3(1)
C
C  FORWARD SWEEP
C
      DO 500 J=2,NP
         DEN   = (A13(J-1)*A21(J-1)-A23(J-1)*A11(J-1)-A(J)*
     1           (A12(J-1)*A21(J-1)-A22(J-1)*A11(J-1)))
         DEN1  = A22(J-1)*A(J)-A23(J-1)
         G11   = (A23(J-1)+A(J)*(A(J)*A21(J-1)-A22(J-1)))/DEN
         G12   = -(A(J)*A(J)+G11*(A12(J-1)*A(J)-A13(J-1)))/DEN1
         G13   = (G11*A13(J-1)+G12*A23(J-1))/A(J)
         G21   = (S2(J)*A21(J-1)-S4(J)*A23(J-1)+A(J)*(S4(J)*
     1           A22(J-1)-S6(J)*A21(J-1)))/DEN
         G22   = (-S2(J)+S6(J)*A(J)-G21*(A(J)*A12(J-1)-A13(J-1)))/DEN1
         G23   = G21*A12(J-1)+G22*A22(J-1)-S6(J)
         A11(J)= 1.0
         A12(J)=-A(J)-G13
         A13(J)= A(J)*G13
         A21(J)= S3(J)
         A22(J)= S5(J)-G23
         A23(J)= S1(J)+A(J)*G23
         R1(J) = R1(J)-G11*R1(J-1)-G12*R2(J-1)-G13*R3(J-1)
         R2(J) = R2(J)-G21*R1(J-1)-G22*R2(J-1)-G23*R3(J-1)
  500 CONTINUE
C
C  BACKWARD SWEEP
C
      DELU(NP) = R3(NP)
      E1       = R1(NP)-A12(NP)*DELU(NP)
      E2       = R2(NP)-A22(NP)*DELU(NP)
      DELV(NP) = (E2*A11(NP)-E1*A21(NP))/(A23(NP)*A11(NP)-A13(NP)*
     1           A21(NP))
      DELF(NP) = (E1-A13(NP)*DELV(NP))/A11(NP)
      J     = NP
  600 J     = J-1
      E3    = R3(J)-DELU(J+1)+A(J+1)*DELV(J+1)
      DEN2  = A21(J)*A12(J)*A(J+1)-A21(J)*A13(J)-A(J+1)*A22(J)*A11(J)+
     1        A23(J)*A11(J)
      DELV(J) = (A11(J)*(R2(J)+E3*A22(J))-A21(J)*R1(J)-E3*A21(J)*A12(J)
     1           )/DEN2
      DELU(J)  =-A(J+1)*DELV(J)-E3
      DELF(J)  = (R1(J)-A12(J)*DELU(J)-A13(J)*DELV(J))/A11(J)
      IF(J .GT. 1) GO TO 600
C
C
      DO 700 J=1,NP
         F(J,2)= F(J,2)+DELF(J)
         U(J,2)= U(J,2)+DELU(J)
         V(J,2)= V(J,2)+DELV(J)
  700 CONTINUE
      U(1,2)= 0.0
      RETURN
      END

C------------------------------------------------------------
      SUBROUTINE OUTPUT
      COMMON /BLC0/ RL,XCTR
      COMMON /BLC7/ ETA(201),DETA(201),A(201)
      COMMON /BLCS/ DLS(200),VW(200),CF(200),THT(200)
      COMMON /BLC2/ NX,NXT,NP,NPT,NTR,NVRS,IT,INVRS,ISWPT
      COMMON /BLC8/ F(201,2),U(201,2),V(201,2),W(201,2),B(201,2)
      COMMON /BLC3/ X(201),UE(200),UE0(200),P1(200),P2(200),GMTR(200)
C
      IF(NX.EQ.1 ) THEN
         DLS(NX)= 0.0
         THT(NX)= 0.0
         CF(NX) = 0.0
         VW(NX) = V(1,2)
      ELSE
C
         CF(NX) = (2*V(1,2)*B(1,2))/SQRT(UE(NX)*X(NX)*RL)
         VW(NX) = 0
         DLS(NX)= (X(NX)*(ETA(NP)-F(NP,2)))/(SQRT(UE(NX)*X(NX)*RL))
         U1     = U(1,2) * (1.0 -U(1,2))
         SUM    = 0.0
         DO 20 J=2,NP
            U2     = U(J,2) * (1.0 -U(J,2))
            SUM    = SUM + A(J) * (U1 + U2)
            U1     = U2
  20     CONTINUE
         THT(NX)= X(NX)/SQRX * SUM
      ENDIF
      WRITE (66,220) (J,ETA(J),F(J,2),U(J,2),V(J,2),B(J,2),J=1,NP)

C write out data to ftn01 which will be used as 
C input data of stability code
      WRITE (1,666) NX, NP
      WRITE (1,888) (ETA(J), J=1, NP)
      WRITE (1,888) (U(J,2), J=1, NP)
      WRITE (1,888) (V(J,2), J=1, NP)
 666  FORMAT(2I5)
 888  FORMAT(6E14.6)
C
C     SHIFT PROFILES FOR THE NEXT STATION
C
      CALL GROWTH(2)
      DO 175 J=1,NPT
            F(J,1)  = F(J,2)
            U(J,1)  = U(J,2)
            V(J,1)  = V(J,2)
            B(J,1)  = B(J,2)
 175  CONTINUE
C
      RETURN
 220  FORMAT('  J ',4X,'ETA',9X,'F',13X,'U',13X,'V',13X,'B'/
     &       (I3,F8.3,4E14.5))
      END

      SUBROUTINE EDDY
      COMMON /BLC0/ RL,XCTR
      COMMON /BLC2/ NX,NXT,NP,NPT,NTR,NVRS,IT,INVRS,ISWPT
      COMMON /BLC3/ X(201),UE(200),UE0(200),P1(200),P2(200),GMTR(200)
      COMMON /BLC7/ ETA(201),DETA(201),A(201)
      COMMON /BLC8/ F(201,2),U(201,2),V(201,2),W(201,2),B(201,2)
      DIMENSION     EDVI(201),FINT(201)
C
      RL2    = SQRT(RL*UE(NX) * X(NX))
      RL4    = SQRT(RL2)
      RL216  = 0.16 * RL2
      UDEL = 0.995*U(NP,2)
C
      DO 10 J=2,NP
      JJ = J
      IF (U(J,2).GT.UDEL) GO TO 15
   10 CONTINUE
      DEL = ETA(NP)
      GO TO 20
   15 DEL = ETA(JJ-1)+(ETA(JJ)-ETA(JJ-1))/(U(JJ,2)-U(JJ-1,2))*(UDEL-
     1      U(JJ-1,2))
   20 DO 25 J=1,NP
      FINT(J) = 1./(1.0+5.5*(ETA(J)/DEL)**6)
   25 CONTINUE
C
      EDVO   = 0.0168*RL2*GMTR(NX)*(U(NP,2)*ETA(NP)-F(NP,2))
      EDVI(1)= 0.0
C
      PLUS = 1./(RL4*ABS(V(1,2))**1.5)*P2(NX)
      CN = SQRT(1.-11.8*PLUS)
      YBAJ   = CN*RL4*SQRT(ABS(V(1,2)))/26.0
C
      DO 70 J=2,NP
         JJ     = J
         YBA    = YBAJ*ETA(J)
         EL     = 1.0
         IF(YBA .LT. 10.0) EL = 1.0 - EXP(-YBA)
C
         EDVI(J) = RL216*GMTR(NX)*(EL*ETA(J))**2 * ABS(V(J,2))
         IF(EDVI(J) .GT. EDVO) GO TO 90
         IF (EDVI(J) .LE. EDVI(J-1)) EDVI(J)= EDVI(J-1)
         B(J,2) = 1.0 + EDVI(J)*FINT(J)
   70 CONTINUE
   90 DO 100 J=JJ,NPT
  100 B(J,2) = 1.0 + EDVO*FINT(J)
      B(1,2) = 1.0
C
      RETURN
      END
