      PROGRAM ISOMC2
C  **** modified to do a ZAMS and use Kurucz, Lejeune or Worthy
C  **** Phoenix or V&C  colours
C  **** YALEISO.F   (BC. 1994)   ---from
C  **** ISOAUTO.FOR (B.C. 1990)  ---from
C  **** COLORISO.FOR (P.D. 1989) ---from
C  **** ISO.FOR --- (REV.E.M.GREEN,1984) --- ****
C  **** ISOCHRONE  CONSTRUCTION  W/ PARROT INTERPOLATION ****
C  **** THIS PROGRAM READS EVOLUTIONARY TRACKS AND CONSTRUCTS
C         ISOCHRONES IN THE (LOGTEFF-LOGL)-PLANE  ****
C  **** AND CALCULATES COLOR ISOCHRONES  ****
C
C modification by B.C. to allow automatic processing of YREC3B iso files
C tailored for use with the Yale isochrone project files
C basic COLORISO.FOR program heavily modified
C PROGRAM ASSUMES NO MORE THAN THE FIRST 99 LINES IN THE ISO FILES ARE
C RESCALEING LINES
c
c SRB 7/03 
c Modified the way the main sequence turnoff EEP is defined. Originally
c the EEP is defined in terms of helium abundance: Y = 1.0 - Z. Since some stars
c fail to reach this value at the MSTO (perhaps because of heavy element 
c diffusion?) we define the MSTO EEP by a small but nonzero helium core mass: 
c Mc = 0.01Msun. 
c
c Also modified the program so that it skips over any input track files
c that are empty or contain less than 1000 models.

C*************BC
C  modified to read in correct [alfa/Fe], etc from the headers.
C modified to automatically process Monte Carlo runs
C
C ON INPUT
C  --- UNITS 11, 12, 13,... etc. are the YREC3B iso files in order of 
C      increasing mass.  SUBROUTINE TRACKIN reads unit 95 (Namelist file) 
C      for information
C **** modified 10/23/91 by B.C. to use more than 12 mass tracks
C
C  **** modified by B.C. March 18/90 to include B. Green's new bolometric
C       corrections
C
C ON OUTPUT
C  --- EACH ISOCHRONE IS PRECEDED BY ONE LINE CONTAINING
C         NPTS,ALFA, OV, AGE, Y, Z, ZEFF, [Fe/H], [alpha/Fe]
C         IN FORMAT(1x,I3,F9.6,F8.4,F7.0,F8.4,2(E12.4),2(F6.2))
C
C  NPTS       ACTUAL NUMBER OF POINTS ON ISOCHRONE                         
C  ALFA       RATIO OF MIXING LENGTH TO PRESSURE SCALE HEIGHT IN           
C             CONVECTIVE ENVELOPE                                          
C  OV         AMOUNT OF OVERSHOOT (IN UNITS OF PRESSURE SCALE HEIGHT)      
C             BELOW THE BASE OF THE CONVECTIVE ENVELOPE                    
C  AGE        AGE IN MYRS                                                  
C  Y          HELIUM ABUNDANCE (MASS FRACTION)                             
C  Z          TRUE Z USED IN THE STELLAR MODEL CONSTRUCTION                
C  ZEFF       EFFECTIVE Z USED TO CONVERT TO [Fe/H] -- TAKES INTO ACCOUNT  
C             ALPHA ELEMENT ENHANCEMENT                                    
C  [Fe/H]     METALLICITY                                                  
C  [alpha/Fe] alpha ELEMENT ENHACEMENT OVER SOLAR RATIO                    
C
C  --- EACH LINE OF THE TABLE INCLUDES
C       J,AX,SMT(J),LOGG,XXT(K),YYT(J),ANUM(J,1),ANUM(J,2)
C          ANUM(J,3) IN FORMAT(I3,A1,F9.6,3F10.7,3E12.6E1)
C       AX='*' IF EEP IS EXTRAPOLATED FROM INPUT TRACKS OR INTERPOLATED
C          BETWEEN WIDELY SPACED TRACKS, AND ' ' OTHERWISE
C       SMT(J)=MASS AT J (SOLAR MASS UNITS)
C       LOGG == Log G
C       XXT(J)=LOGTEFF, YYT(J)=LOG(L/LSUN) AT J
C       ANUM(J,K)=RELATIVE NUMBER OF STARS IN MASS INTERVAL (J,J+1)
C          FOR IMF INDEX S(K).  N.B. S USED HERE=-(X+1) WHERE X IS
C          TINSLEY'S INDEX.
C
C
C introduce a parameter for the maixmum number of models for a single mass
      PARAMETER(NMODMAX=11000,NMAXISO=800)
      PARAMETER (MAXIN=30, NMAXAGE=60)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 ISOAGE(NMAXAGE)
      INTEGER QT(12),QSUM,QTMS,QTB,QBMP
      CHARACTER*1 AX(NMAXISO)
      LOGICAL LEND,LLARGE,LINEAR,LSMALL
      COMMON/PARMIN/REFZ,SOLBOL,ALFA
      COMMON/PAROT/IS
      COMMON/QTCOM/QT
      COMMON/ISOIN/ISOAGE,NUMAGE,NTRK
      COMMON/EEPTRK/XX(NMODMAX,MAXIN),YY(NMODMAX,MAXIN),
     $              TL(NMODMAX,MAXIN)
      DIMENSION SML(MAXIN),NEEP(MAXIN,13),SP1(3),
     1 LASTPT(MAXIN),XXT(NMODMAX),YYT(NMODMAX),
     3 BOL(NMODMAX),SMT(NMODMAX),ANUM(NMAXISO,3),TQ(MAXIN),SQ(MAXIN),
     4 XQ(MAXIN),YQ(MAXIN),IQ(MAXIN),NDIVI(MAXIN), TEMP(NMODMAX)
C SP1 = S + 1, WHERE S IS IMF EXPONENT: N(M) = M**S
      DATA SP1/1.,-1.35,-3./
      DATA EPS/1.D-8/
C
      IS=1
C
      CALL TRACKIN(SML,Y,FEH,AFE,Z,ZEFF,OV,NEEP,NDIVI,LASTPT,NTAUP)
c      WRITE(*,111) Y,Z,ZEFF,FEH,AFE
 111  FORMAT(' Y =',F4.2,', Z =',F7.5,
     *     ', ZEFF =',F7.5,', [Fe/H] = ',F6.3,', [alpha/Fe] = ',F5.2)
C QTMS = EEP AT END OF MS; QTB = EEP AT BASE OF RGB;
C QBMP = EEP AT RGB BUMP
      QTMS=1+QT(1)+QT(2)+QT(3)
      QTB=QTMS+QT(4)
      NQB=4
      QBMP=QTB+QT(5)+1
C
C     MAIN INTERPOLATION LOOP BETWEEN TRACKS
C
      SMT(1)=-1.
C  READ IN ISOCHRONE AGE
      DO 301 IAGE=1,NUMAGE
         TAU=ISOAGE(IAGE)
         TAUL=DLOG10(TAU+1.D0)
         IGB=0
         QSUM=2
         L=0
         IF(TAU.EQ.0.) THEN
            CALL ZAMSISO(SML,SP1,Y,FEH,AFE,Z,ZEFF,OV)
            GOTO 301
         ENDIF
         DO 81 J=2,NTAUP
C CHECK FIRST SUBEEP IN EACH MAJOR INTERVAL TO SEE WHICH TRACKS ARE DEFINED
            IF(J.EQ.QSUM) THEN
               CALL TRKDEF(SML,TAUL,J,NQB,QTMS,LASTPT,L,QSUM,NTAUP,
     $                  IT,IQ,SQ,LLARGE,LEND)
               IF(LEND) GOTO 301
               IF(LLARGE) THEN
C MASS TOO LARGE OR END OF DEFINED TRACKS
                  SMT(J)=0.
                  NT=J-1
                  GO TO 88
               ENDIF
            ENDIF
C FOR EACH SUBEEP J, PUT TEFF, LOGL, AND TIME FOR DEFINED TRACKS INTO
C CONTINUOUS ARRAYS FOR INTERPOLATION
            DO 82 I=1,IT
               TQ(I)=TL(J,IQ(I))
               XQ(I)=XX(J,IQ(I))
               YQ(I)=YY(J,IQ(I))
 82         CONTINUE
C CHECK FOR LOW MASS EXTRAPOLATION, AND DETERMINE K SUCH THAT
C TQ(K) < TAUL < TQ(K-1)
            CALL MASSK(TQ,TAUL,J,IT,SMT,K,LINEAR,LSMALL,LEND)
            IF(LEND) GOTO 301
            IF(LSMALL) GOTO 81
            IF(.NOT.LINEAR) THEN
C
C INTERPOLATE IN TIME TO GET MASS FOR EEP J
               CALL PARROT(TQ,SQ,K,IT,TAUL,GAMMA,1)

C USE MASS TO INTERPOLATE FOR TEFF AND LOG L SINCE THESE ARE SLIGHTLY
C MORE LINEAR WITH RESP TO LOG M THAN TO LOG TIME
               CALL PARROT(SQ,XQ,K,IT,GAMMA,XXT(J),1)
               CALL PARROT(SQ,YQ,K,IT,GAMMA,YYT(J),1)
            ELSE
C LINEAR EXTRAPOLATION OR THERE ARE ONLY TWO TRACKS FOR INTERPOLATION
               FRAC=(TQ(K-1)-TAUL)/(TQ(K-1)-TQ(K))
               GAMMA=SQ(K-1)-FRAC*(SQ(K-1)-SQ(K))
               XXT(J)=XQ(K-1)-FRAC*(XQ(K-1)-XQ(K))
               YYT(J)=YQ(K-1)-FRAC*(YQ(K-1)-YQ(K))
            ENDIF

            BOL(J)=SOLBOL-(2.5*YYT(J))
            SMT(J)=10.**GAMMA
C LABEL EXTRAPOLATED POINTS AND LONG INTERPOLATIONS WITH '*'
            AX(J)='*'
            IF(GAMMA.GE.SQ(1)) THEN
               DO 52 I=2,IT
                  IF (GAMMA.LE.SQ(I)) THEN
                     IF (IQ(I)-IQ(I-1).EQ.1) AX(J)=' '
                     GO TO 56
c                     IF (IQ(I)-IQ(I-1).EQ.1) then
c                        AX(J) = ' ' 
c                        go to 56
c                    else
c                        write(*,183) IT, J, (IQ(III), III=1, IT)
c183                     format('IT', I2, ' EEP', I4, 'IQ =', 9I2)
c                     end if             
                  ENDIF
 52            CONTINUE
            ENDIF

 56         IF(SMT(J).LT.SMT(J-1)) THEN
C IF MASS IS NONMONOTONIC, FIX POSSIBLE NUMERICAL ROUNDOFF OR RGB ERRORS
C DUE TO TOO FEW POINTS ON PUBLISHED SG GIANT TRACKS
               SMT(J)=SMT(J)+EPS
               IF(J.EQ.QBMP.OR.J.EQ.QBMP+1.OR.J.EQ.NTAUP) 
     $              SMT(J)=SMT(J-1)
            ENDIF
 81      CONTINUE
         NT=NTAUP
         IF(NT.GT.QTB) IGB=-1
C
C END OF INTERPOLATION LOOP
C
C ELIMINATE MASS,LOGL DISCONTINUITIES ACROSS EEP BOUNDARY AT BASE OF RGB
C DUE TO INTERPOLATION CHANGING FROM QUADRATIC TO LINEAR.
         IF(NT.GT.QTB) THEN
            DELM=SMT(QTB)+3.*(SMT(QTB+2)-SMT(QTB+1))-SMT(QTB+3)
            DELT=XXT(QTB)+3.*(XXT(QTB+2)-XXT(QTB+1))-XXT(QTB+3)
            DELL=YYT(QTB)+3.*(YYT(QTB+2)-YYT(QTB+1))-YYT(QTB+3)
            DO 45 J=QTB+1,NT
               SMT(J)=SMT(J)+DELM
               XXT(J)=XXT(J)+DELT
               YYT(J)=YYT(J)+DELL
 45         CONTINUE
         ENDIF
C
C FIND ENDPOINTS OF CALCULATED ISOCHRONE
C NT=HIGHEST DEFINED ISOCHRONE POINT; N0=LOWEST DEFINED POINT
 88      DO 84 J=2,NT
            IF(SMT(J).GT.0.D0) GO TO 83
 84      CONTINUE
         J = NT
 83      N0=J
         N2=NT-1
         NPTS=NT-N0+1
C
C COMPUTE THE NUMBER OF STARS (ANUM(J,M)) IN THE INTERVAL BETWEEN
C EEP J AND EEP #J+1
         CALL NUMSTAR(SP1,SMT,ANUM,NT,N0,N2)
         NPTS=N2-N0+2
C
C ISOCHRONES ARE WRITTEN TO OUTPUT FILE
         CALL OUTISO(XXT,YYT,SMT,BOL,ANUM,Y,FEH,AFE,Z,ZEFF,OV,TAU,
     $               IGB,NPTS,N0,NTAUP,AX)
C
C NEXT AGE
 301  CONTINUE
      write(*,*) done
      END
C
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE PARROT (X,Y,J,N,X0,Y0,K)
C
C PARROT INTERPOLATES BETWEEN A SET OF POINTS IN ARRAYS X AND Y
C (MONOTONIC IN X) TO PRODUCE A VERY SMOOTH CURVE.
C THE FIRST AND SECOND DERIVATIVES ARE EVERYWHERE CONTINUOUS, AND
C POLYNOMIAL "WIGGLES" BETWEEN POINTS ARE ALMOST ENTIRELY ELIMINATED.
C PARROT FINDS ONE ROTATED PARABOLA PASSING THROUGH THE POINTS J-2, J-1,
C J, WITH VERTEX EXACTLY AT J-1, AND A SECOND PARABOLA PASSING THROUGH J
C J, AND J+1, WITH VERTEX AT J, THEN DETERMINES THE Y VALUES CORRESPONDI
C TO X0 FOR EACH CURVE.  Y0 IS THE WEIGHTED SUM OF THE TWO Y VALUES, WIT
C THE WEIGHTS DETERMINED BY THE DISTANCE OF THE INTERPOLATED POINT FROM
C POINTS J AND J-1, RESPECTIVELY.
C (IN THE FIRST AND LAST INTERVALS, THE RESULT IS A WEIGHTED AVERAGE OF
C A PARABOLA AND A STRAIGHT LINE.)
C
C  INPUT PARAMETERS (UNCHANGED BY PARROT) :
C
C     X(I),Y(I) - ARRAYS OF KNOWN POINTS TO BE INTERPOLATED; X MONOTONIC
C             N - LARGEST DEFINED INDEX OF X(I) AND Y(I)
C            X0 - POINT FOR WHICH INTERPOLATED VALUE IS TO BE RETURNED
C             J - X(J) AND X(J-1) MUST BRACKET X0
C             K - SET K.EQ.1 IF CALLING PARROT FOR THE FIRST TIME IN THI
C                 INTERVAL; SET K.NE.1 IF X,Y, AND J ARE THE SAME AS THE
C                 PREVIOUS CALL TO PARROT (FOR THIS VALUE OF IS).
C                 AS PARROT IS RELATIVELY SLOW, MOST EFFICIENT USE REQUI
C                 ORDERING X0 AND SETTING K.NE.1 WHENEVER POSSIBLE.
C           [IS - SETTING IS=1, 2, OR 3 IN THE CALLING PROGRAM (WITH
C                 COMMON/PAROT/IS) ALLOWS THREE DIFFERENT SETS OF PARAME
C                 TO BE REMEMBERED SIMULTANEOUSLY]
C
C  OUTPUT PARAMETER :
C            Y0 - RETURNED VALUE CORRESPONDING TO X0
C
C
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/PAROT/IS
      COMMON/SAVE/CW(2,3),SW(2,3),A(2,3),B(2,3),FX(2,3),FY(2,3),
     1     ILINE(2)
      DIMENSION X(N),Y(N),RES(2),RT(3)
C
      DO 300 I=1,2
         J0=I+J-2
         XO=X0-X(J0)
         IF(K.EQ.1) GO TO 100
C SKIP PARABOLA CALCULATIONS IF X AND Y POINTS ARE THE SAME AS LAST TIME
         XO=XO/FX(I,IS)
         Y1=(Y(J+1-I)-Y(J0))/FY(I,IS)
         GO TO 200
C
C TRANSLATE COORDINATES TO VERTEX OF EACH FUTURE PARABOLA, WITH X1
C ALWAYS CHOSEN TO BE CLOSEST TO X0 (+ LINEAR EXTRAPOLATION PAST ENDS)
100      J1=J-I+1
         X1=X(J1)-X(J0)
         Y1=Y(J1)-Y(J0)
         IF(J0.GT.1) GO TO 110
         X3=X(1)-X(2)
         Y3=Y(1)-Y(2)
         GO TO 125
110      IF(J0.LT.N) GO TO 120
         X3=X(N)-X(N-1)
         Y3=Y(N)-Y(N-1)
         IF(X3.EQ.0.0.AND.Y3.EQ.0.0) PRINT *,'+++',N,X(N),X(N-1)
         IF(X3.EQ.0.0.AND.Y3.EQ.0.0) PRINT *,'+++',N,Y(N),Y(N-1)
         GO TO 125
120      J3=J+3*I-5
         X3=X(J3)-X(J0)
         Y3=Y(J3)-Y(J0)
C SCALE X COORDINATES TO APPROXIMATELY SAME SIZE AS Y COORDINATES
125      FX(I,IS)=DMAX1(DABS(X3),DABS(X1))
         FY(I,IS)=DMAX1(DABS(Y3),DABS(Y1))
         IF(FY(I,IS).EQ.0.) FY(I,IS)=1.
         X1=X1/FX(I,IS)
         X3=X3/FX(I,IS)
         XO=XO/FX(I,IS)
         Y1=Y1/FY(I,IS)
         Y3=Y3/FY(I,IS)
C DETERMINE ANGLE OF ROTATION W SUCH THAT A PARABOLA OF THE FORM Y=A*X**
C PASSES THRU THE THREE POINTS (IF THETA IS THE ANGLE FROM THE ORIGIN
C AND (X1,Y1) TO THE NEW (ROTATED) Y AXIS, NEED TO SOLVE A CUBIC EQUATIO
C IN TAN(THETA).)
         D1=X1*X1+Y1*Y1
         D3=X3*X3+Y3*Y3
         R=DSQRT(D3/D1)
         D=DSQRT(D1*D3)
         CA=(X1*X3+Y1*Y3)/D
         SA=(X1*Y3-X3*Y1)/D
         IF(DABS(SA).LT..001) GO TO 150
         RCA=R*CA
         P=CA*(1.-RCA)/SA
         P2=P*P
         Q=2.*RCA-P2/3.
         R=2.*(P2*P/27.-P*RCA/3.)-R*SA
         Q3=Q*Q*Q/27.
         ROOT=.25*R*R+Q3
         IF(ROOT) 130,129,128
128      ROOT=DSQRT(ROOT)
         THETA=DATAN(QR(-0.5*R+ROOT)+QR(-0.5*R-ROOT)-P/3.)
         GO TO 180
C EVALUATE 2 SOLUTIONS; FIND ANGLE BETWEEN 0 AND A
129      C=QR(-R/2.)
         TA=SA/CA
         LL=1
         RT(1)=2.*C-P/3.
         RT(2)=-C-P/3.
         IF(RT(2)*(TA-RT(2)).GT.RT(1)*(TA-RT(1))) LL=2
         GO TO 132
C EVALUATE 3 SOLUTIONS; FIND ANGLE BETWEEN 0 AND A
130      PHI=DACOS(-R/2./DSQRT(-Q3))/3.
         C=2.*DSQRT(-Q/3.)
         TA=SA/CA
         LAST=-1.
         DO 131 L=1,3
            RT(L)=C*DCOS(PHI+(L-1)*2.0943951) - P/3.
            TST=RT(L)*(TA-RT(L))
            IF(TST.LT.LAST) GO TO 131
            LL=L
            LAST=TST
131      CONTINUE
132      THETA=DATAN(RT(LL))
         GO TO 180
150      THETA=SA
         IF(CA.LT.0.) THETA=3.1415927-SA
         THETA=THETA/2.
180      W=THETA-DATAN2(X1,Y1)
C SAVE COS AND SIN OF W, PLUS A AND B
         CW(I,IS)=DCOS(W)
         SW(I,IS)=DSIN(W)
         XP=X1*CW(I,IS)+Y1*SW(I,IS)
         YP=Y1*CW(I,IS)-X1*SW(I,IS)
         CW2=CW(I,IS)*CW(I,IS)
         XP2=XP*XP
         C=4.*SW(I,IS)*YP
C NOTE C AND XP2 CAN NEVER BE SIMULTANEOUSLY ZERO, NOR XP2 AND YP,
C NOR YP AND CW2, NOR C AND CW2 (UNLESS X1=0, WHICH IS NOT LEGAL
C IN PARROT)
         IF(DABS(C).GT.1.D8*XP2*CW2) GO TO 190
         A(I,IS)=YP/XP2
         B(I,IS)=C/XP2
         ILINE(I)=0
         GO TO 200
190      A(I,IS)=XP2/YP
         B(I,IS)=XP2/C
         ILINE(I)=1
C
C CALCULATE RESULTING Y FOR THIS PARABOLA
200      IF(X0.NE.X(J-1)) GO TO 201
           IF(I*K.EQ.1) GO TO 300
           Y0=Y(J-1)
           RETURN
201      IF(X0.NE.X(J)) GO TO 202
           IF(I*K.EQ.1) GO TO 300
           Y0=Y(J)
           RETURN
C
202      CW2=CW(I,IS)*CW(I,IS)
         IF(ILINE(I).EQ.1) GO TO 250
         C=B(I,IS)*XO
         AA=A(I,IS)
         IF(DABS(C).GT.1.D10*CW2) GO TO 240
210      IF(DABS(C).LT..0001*CW2) GO TO 230
         IF(DABS(SW(I,IS)).LT..0001) GO TO 220
C FIND INTERSECTION OF PARABOLA AND X0 LINE IN ROTATED COORDINATES;
C CALCULATE RESULT BACK IN UNROTATED COORDINATES.
         C=CW2-C
         IF(C.LT.0.) GO TO 400
         ROOT=DSQRT(C)
         C=XO*CW(I,IS)/SW(I,IS)
         D=2.*AA*SW(I,IS)*SW(I,IS)
         RES(I)=(CW(I,IS)-ROOT)/D - C
         ALT=(CW(I,IS)+ROOT)/D - C
C CHOOSE ROOT IN INTERVAL [Y(J-1),Y(J)] (OR NEAREST TO THAT INTERVAL)
         IF((Y1-ALT)*ALT.GT.(Y1-RES(I))*RES(I)) RES(I)=ALT
         GO TO 280
C
C SPECIAL CASES 220,230,240,250,260 BELOW:
C ROTATION ANGLE CLOSE TO ZERO, BUT A NOT SMALL
220      RES(I)=XO*(AA*XO/(1.-C/2.)+SW(I,IS))/CW(I,IS)
         GO TO 280
C PARABOLA DEGENERATES TOWARDS A STRAIGHT HORIZONTAL LINE
230      C=XO/CW(I,IS)
         RES(I)=C*(AA*C/CW(I,IS)+SW(I,IS))
         GO TO 280
C XO VERY LARGE
240      XA=XO/A(I,IS)
         GO TO 260
C ILINE = 1
250      XA=XO*A(I,IS)
         IF(DABS(XO).GT.1.D8*DABS(B(I,IS))*CW2) GO TO 260
C XO IS ALWAYS NONZERO, THEREFORE, B, XP2, AND A MUST BE NONZERO;
C RETURN TO NORMAL CASE
         C=XO/B(I,IS)
         AA=1./A(I,IS)
         GO TO 210
C PARABOLA DEGENERATES TOWARDS A FOLDED VERTICAL LINE
260      ROOT=DSQRT(-XA/SW(I,IS))
         C=XO*CW(I,IS)
         RES(I)=(ROOT-C)/SW(I,IS)
         ALT=(-ROOT-C)/SW(I,IS)
         IF((Y1-ALT)*ALT.GT.(Y1-RES(I))*RES(I)) RES(I)=ALT
C
280      RES(I)=RES(I)*FY(I,IS)+Y(J0)
300   CONTINUE
C
C WEIGHT TWO PARABOLAS ACCORDING TO DISTANCE FROM (J-1) AND (J)
      DX=X(J)-X0
      DY=Y(J)-RES(2)
      WT1=DSQRT(DX*DX+DY*DY)
      DX=X0-X(J-1)
      DY=RES(1)-Y(J-1)
      WT2=DSQRT(DX*DX+DY*DY)
      Y0=(RES(1)*WT1+RES(2)*WT2)/(WT1+WT2)
      RETURN
C
400   write(*,*)' X0 NOT BRACKETED BY X(J-1) AND X(J)'
      RETURN
      END
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      FUNCTION QR(X)
      IMPLICIT REAL*8(A-H,O-Z)
      IF(X.EQ.0.) THEN 
         QR = 0.0D0
      ELSE
         QR = DSIGN(DEXP(DLOG(DABS(X))/3.D0),X)
      ENDIF
      RETURN
      END
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE COLOR(FEH,TEFF,LOGG,COLR)
C
C INTERPOLATES IN COLOR.TBL
C    INPUT: FEH, TEFF, LOG G  (UNCHANGED BY COLOR)
C    OUTPUT: COLR(5), WHERE COLR(1) = BC
C                           COLR(2) = U-B
C                           COLR(3) = B-V
C                           COLR(4) = V-R
C                           COLR(5) = R-I
C
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 FEH,TEFF,LOGG,COLR(8)
      DIMENSION COL(5),TABL(15,50,10,5),GLOG(15,50,10),JMAX(50,10),
     1  IMAX(10),CZ(10),TEF(50,10),COLT(50),COLZ(10)
      COMMON/PARMIN/REFZ,SOLBOL,ALFA
      DATA NCOLOR/0/
      SAVE 
C
C SET UP AND READ IN COLOR TABLE ON FIRST CALL
C
      IF(NCOLOR.EQ.0) THEN
C
        OPEN(UNIT=61,STATUS='OLD')
C
C N.B. REFZ IS THE REFERENCE (SOLAR) METALLICITY
        SOLARZ=DLOG10(REFZ)
        NCOLOR=5
C
C READ IN COLOR TABLE
        K=0
        TLAST=0.
        FLAST=10.
        READ(61,*)
1       READ(61,*,END=20) FE,TE,G,(COL(M),M=1,NCOLOR)
        IF(TE.EQ.TLAST) GOTO 15
        IF(FE.EQ.FLAST) GOTO 10
        K=K+1
        KMAX=K
        CZ(K)=FE+SOLARZ
        FLAST=FE
        I=0
10      I=I+1
        IMAX(K)=I
        TEF(I,K)=DLOG10(TE)
        TLAST=TE
        J=0
15      J=J+1
        JMAX(I,K)=J
        GLOG(J,I,K)=G
        DO M=1,NCOLOR
          TABL(J,I,K,M)=COL(M)
        ENDDO
        GOTO 1
C
      ENDIF
C
C INTERPOLATE IN LOG G, LOG TEFF, AND [FE/H]
C
20    ZZ=FEH+SOLARZ
      IF((CZ(KMAX)-ZZ)*(ZZ-CZ(1)).LT.0.) GO TO 980
      KK=1
30    KK=KK+1
      IF(CZ(KK).GT.ZZ) GOTO 30
      TE=DLOG10(TEFF)
      G=LOGG
      DO M=1,NCOLOR
        K0=MAX0(1,KK-2)
        KEND=MIN0(KMAX,KK+1)
        NK=KEND-K0+1
        DO K=K0,KEND
           IF((TEF(IMAX(K),K)-TE)*(TE-TEF(1,K)).LT.0.) THEN
C put in a hack -- instead of aborting when outside colour table, put
C in -9.999
              write(*,*)'Teff=',teff,' OUTSIDE LIMITS OF COLOR FILE'
              DO 33 III=1,5
                 COLR(III)=-9.999
 33           CONTINUE
              RETURN
          ENDIF
          II=1
40        II=II+1
          IF(TE.GT.TEF(II,K).AND.II.LT.IMAX(K)) GO TO 40
          I0=MAX0(1,II-2)
          IEND=MIN0(IMAX(K),II+1)
          NI=IEND-I0+1
          DO I=I0,IEND
            IF((GLOG(JMAX(I,K),I,K)-G)*(G-GLOG(1,I,K)).LT.0.) GO TO 990
            JJ=1
45          JJ=JJ+1
            IF(G.GT.GLOG(JJ,I,K)) GO TO 45
            J0=MAX0(1,JJ-2)
            JEND=MIN0(JMAX(I,K),JJ+1)
            NJ=JEND-J0+1
            CALL PARROT(GLOG(J0,I,K),TABL(J0,I,K,M),JJ-J0+1,NJ,G,
     1                                                 COLT(I),1)
          ENDDO
          CALL PARROT(TEF(I0,K),COLT(I0),II-I0+1,NI,TE,COLZ(K),1)
        ENDDO
        CALL PARROT(CZ(K0),COLZ(K0),KK-K0+1,NK,ZZ,COL(M),1)
        COLR(M)=COL(M)
      ENDDO
C
      RETURN
C
 980  write(*,*)'[Fe/H] =',FEH,'OUTSIDE LIMITS OF COLORFILE'
      DO 333 III=1,5
         COLR(III)=-9.999
 333  CONTINUE
      RETURN
 985  write(*,*)'Teff=',teff,' OUTSIDE LIMITS OF COLOR CALIBRATION FILE'
      STOP
 990  write(*,*)G,', ',XXT,', Z =',ZZ,
     1  ' OUTSIDE LIMITS OF COLOR CALIBRATION FILE'
      STOP
      END
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE KURCOL(FEH,TEFF,LOGG,COLR)
C
C INTERPOLATES IN the KURUCZ COLOR.TBL
C    INPUT: FEH, TEFF, LOG G  (UNCHANGED BY COLOR)
C    OUTPUT: COLR(8), WHERE COLR(1) = BC
C                           COLR(2) = U-B
C                           COLR(3) = B-V
C                           COLR(4) = V-R
C                           COLR(5) = V-I
C                           COLR(6) = V-J 
C                           COLR(7) = V-H
C                           COLR(8) = V-K 
C
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 FEH,TEFF,LOGG,COLR(8)   ! use this line for all colours
C      REAL*8 FEH,TEFF,LOGG,COLR(5)
      INTEGER ONE
      DIMENSION COL(8),TABL(15,70,10,8),GLOG(15,70,10),JMAX(70,10),
     1  IMAX(10),CZ(10),TEF(70,10),COLT(70),COLZ(10)
C TABL(J,I,K,M):   M==colors
C                  K==[Fe/H]
C                  I==Teff
C                  J=Log G
C ***BC 5/92 dimension statement for RYI colour table
c      DIMENSION COL(5),TABL(15,50,10,5),GLOG(15,50,10),JMAX(50,10),
c     1  IMAX(10),CZ(10),TEF(50,10),COLT(50),COLZ(10)
      COMMON/PARMIN/REFZ,SOLBOL,ALFA
      DATA NCOLOR/0/
      DATA ONE/60/
      SAVE
C
C SET UP AND READ IN COLOR TABLE ON FIRST CALL
C
      IF(NCOLOR.EQ.0) THEN
C
        OPEN(UNIT=ONE,STATUS='OLD') 
C
C N.B. REFZ IS THE REFERENCE (SOLAR) METALLICITY
        SOLARZ=DLOG10(REFZ)
        NCOLOR=8           !use this line for all colours
C
C READ IN COLOR TABLE
        K=0
        TLAST=0.
        FLAST=10.
        READ(ONE,*)
1       READ(ONE,*,END=20) FE,TE,G,(COL(M),M=1,NCOLOR)
c        write(*,*)fe,te,g
        IF(TE.EQ.TLAST) GOTO 15
        IF(FE.EQ.FLAST) GOTO 10
        K=K+1
        KMAX=K
        CZ(K)=FE+SOLARZ
        FLAST=FE
        I=0
10      I=I+1
        IMAX(K)=I
        TEF(I,K)=DLOG10(TE)
        TLAST=TE
        J=0
15      J=J+1
        JMAX(I,K)=J
        GLOG(J,I,K)=G
c***BC put in test
        if(K.GT.10 .OR. I.GT.70 .OR. J.GT.15) THEN
           write(*,*) 'need to resize arrays in KURCOL'
           write(*,*)'i,j,k=',i,j,k
           stop
        endif
        DO M=1,NCOLOR
          TABL(J,I,K,M)=COL(M)
        ENDDO
        GOTO 1
C
      ENDIF
C
C INTERPOLATE IN LOG G, LOG TEFF, AND [FE/H]
C
20    ZZ=FEH+SOLARZ
C ****BC 5/92
C if temp < 3500 K, just put -9.999 for the colours and exit
      IF(TEFF.LT.3500.0) THEN
         DO 666 III=1,NCOLOR
            COLR(III)=-9.999
 666     CONTINUE
         RETURN
      ENDIF
C ***BC 5/92 end of modification
      IF((CZ(KMAX)-ZZ)*(ZZ-CZ(1)).LT.0.) GO TO 980
      KK=1
30    KK=KK+1
      IF(CZ(KK).GT.ZZ) GOTO 30
      TE=DLOG10(TEFF)
      G=LOGG
      DO M=1,NCOLOR
        K0=MAX0(1,KK-2)
        KEND=MIN0(KMAX,KK+1)
        NK=KEND-K0+1
        DO K=K0,KEND
c          write(*,*)  10**tef(imax(k),k), 10**tef(1,k)
c          write(*,*) 10**te
          IF((TEF(IMAX(K),K)-TE)*(TE-TEF(1,K)).LT.0.) GO TO 985
          II=1
40        II=II+1
          IF(TE.GT.TEF(II,K).AND.II.LT.IMAX(K)) GO TO 40
          I0=MAX0(1,II-2)
          IEND=MIN0(IMAX(K),II+1)
          NI=IEND-I0+1
          DO I=I0,IEND
c            write(*,*) jmax(i,k), glog(jmax(i,k),i,k), glog(1,i,k)
c            write(*,*) g
C ***BC 5/92 
C Kurucz does not have high G low Teff tables, so just put -9.999 for colour
C and continue
c            IF((GLOG(JMAX(I,K),I,K)-G)*(G-GLOG(1,I,K)).LT.0.) GO TO 990
            IF((GLOG(JMAX(I,K),I,K)-G)*(G-GLOG(1,I,K)).LT.0.) THEN
               DO 555 III=1,NCOLOR
                  COLR(III)=-9.999
 555           CONTINUE
               RETURN
            ENDIF
C ***BC end of modification
            JJ=1
45          JJ=JJ+1
            IF(G.GT.GLOG(JJ,I,K)) GO TO 45
            J0=MAX0(1,JJ-2)
            JEND=MIN0(JMAX(I,K),JJ+1)
            NJ=JEND-J0+1
            CALL PARROT(GLOG(J0,I,K),TABL(J0,I,K,M),JJ-J0+1,NJ,G,
     1                                                 COLT(I),1)
          ENDDO
          CALL PARROT(TEF(I0,K),COLT(I0),II-I0+1,NI,TE,COLZ(K),1)
        ENDDO
        CALL PARROT(CZ(K0),COLZ(K0),KK-K0+1,NK,ZZ,COL(M),1)
        COLR(M)=COL(M)
      ENDDO
C
      RETURN
C
 980  write(*,*)'[Fe/H] =',FEH,
     1  ' OUTSIDE LIMITS OF COLOR CALIBRATION FILE'
      write(*,*)'feh problem'
      STOP
 985  write(*,*)TEFF,', [Fe/H] =',FEH, 
     1  ' OUTSIDE LIMITS OF COLOR CALIBRATION FILE'
      write(*,*)'teff problem'
      STOP
c990  write(*,*)LOGG,', ',TEFF,', [Fe/H] =',FEH,
c    1  ' OUTSIDE LIMITS OF COLOR CALIBRATION FILE'
c     STOP
      END
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE LEJCOL(FEH,TEFF,LOGG,COLR)
C
C INTERPOLATES IN the Lejeune COLOR.TBL
C    INPUT: FEH, TEFF, LOG G  (UNCHANGED BY COLOR)
C    OUTPUT: COLR(8), WHERE COLR(1) = BC
C                           COLR(2) = U-B
C                           COLR(3) = B-V
C                           COLR(4) = V-R
C                           COLR(5) = V-I
C                           COLR(6) = V-K 
C                           COLR(7) = R-I
C                           COLR(8) = I-K
c  there are lots more colors in the table that we could read in...
C
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 FEH,TEFF,LOGG,COLR(8)   ! use this line for all colours
      INTEGER ONE
      DIMENSION COL(8),TABL(15,70,20,8),GLOG(15,70,20),JMAX(70,20),
     1  IMAX(20),CZ(20),TEF(70,20),COLT(70),COLZ(20)
C TABL(J,I,K,M):   M==colors
C                  K==[Fe/H]
C                  I==Teff
C                  J=Log G
C ***BC 5/92 dimension statement for RYI colour table
c      DIMENSION COL(5),TABL(15,50,10,5),GLOG(15,50,10),JMAX(50,10),
c     1  IMAX(10),CZ(10),TEF(50,10),COLT(50),COLZ(10)
      COMMON/PARMIN/REFZ,SOLBOL,ALFA
      DATA NCOLOR/0/
      DATA ONE/62/
      SAVE
C
C SET UP AND READ IN COLOR TABLE ON FIRST CALL
C
      IF(NCOLOR.EQ.0) THEN
C
        OPEN(UNIT=ONE,STATUS='OLD') 
C
C N.B. REFZ IS THE REFERENCE (SOLAR) METALLICITY
        SOLARZ=DLOG10(REFZ)
        NCOLOR=8           !use this line for all colours
C
C READ IN COLOR TABLE
        K=0
        TLAST=0.
        FLAST=10.
        READ(ONE,*)
        READ(ONE,*)
1       READ(ONE,*,END=20) TE,G,FE,xxx,(COL(M),M=1,NCOLOR)
c        write(*,*)fe,te,g
c        write(*,*)fe,te,k
        IF(TE.EQ.TLAST) GOTO 15
        IF(FE.EQ.FLAST) GOTO 10
        K=K+1
        KMAX=K
        CZ(K)=FE+SOLARZ
        FLAST=FE
        I=0
10      I=I+1
        IMAX(K)=I
        TEF(I,K)=DLOG10(TE)
        TLAST=TE
        J=0
15      J=J+1
        JMAX(I,K)=J
        GLOG(J,I,K)=G
c***BC put in test
        if(K.GT.20 .OR. I.GT.70 .OR. J.GT.15) THEN
           write(*,*) 'need to resize arrays in KURCOL'
           write(*,*)'i,j,k=',i,j,k
           stop
        endif
        DO M=1,NCOLOR
          TABL(J,I,K,M)=COL(M)
        ENDDO
        GOTO 1
C
      ENDIF
C
C INTERPOLATE IN LOG G, LOG TEFF, AND [FE/H]
C
20    ZZ=FEH+SOLARZ
C ****BC 5/92
C if temp < 2000 K, just put -9.999 for the colours and exit
      IF(TEFF.LT.2000.0) THEN
         DO 666 III=1,NCOLOR
            COLR(III)=-9.999
 666     CONTINUE
         RETURN
      ENDIF
C ***BC 5/92 end of modification
      IF((CZ(KMAX)-ZZ)*(ZZ-CZ(1)).LT.0.) GO TO 980
      KK=1
30    KK=KK+1
      IF(CZ(KK).GT.ZZ) GOTO 30
      TE=DLOG10(TEFF)
      G=LOGG
      DO M=1,NCOLOR
        K0=MAX0(1,KK-2)
        KEND=MIN0(KMAX,KK+1)
        NK=KEND-K0+1
        DO K=K0,KEND
c          write(*,*)  10**tef(imax(k),k), 10**tef(1,k)
c          write(*,*) 10**te
          IF((TEF(IMAX(K),K)-TE)*(TE-TEF(1,K)).LT.0.) GO TO 985
          II=1
40        II=II+1
          IF(TE.GT.TEF(II,K).AND.II.LT.IMAX(K)) GO TO 40
c  SRB 2/03 
c  To keep from going outside the lejeune color table as often, use only
c  two temperature points in the interpolation.
c         I0=MAX0(1,II-2)
c         IEND=MIN0(IMAX(K),II+1)
          I0=MAX0(1,II-1)
          IEND=MIN0(IMAX(K),II)
c  SRB end changes
          NI=IEND-I0+1
          DO I=I0,IEND
c            write(*,*) jmax(i,k), glog(jmax(i,k),i,k), glog(1,i,k)
c            write(*,*) g
C ***BC 5/92 
C Kurucz does not have high G low Teff tables, so just put -9.999 for colour
C and continue
C CHECK THIS!!!!
c            IF((GLOG(JMAX(I,K),I,K)-G)*(G-GLOG(1,I,K)).LT.0.) GO TO 990
            IF((GLOG(JMAX(I,K),I,K)-G)*(G-GLOG(1,I,K)).LT.0.) THEN
               DO 555 III=1,NCOLOR
                  COLR(III)=-9.999
 555           CONTINUE
               RETURN
            ENDIF
C ***BC end of modification
            JJ=1
45          JJ=JJ+1
            IF(G.GT.GLOG(JJ,I,K)) GO TO 45
            J0=MAX0(1,JJ-2)
            JEND=MIN0(JMAX(I,K),JJ+1)
            NJ=JEND-J0+1
            CALL PARROT(GLOG(J0,I,K),TABL(J0,I,K,M),JJ-J0+1,NJ,G,
     1                                                 COLT(I),1)
          ENDDO
          CALL PARROT(TEF(I0,K),COLT(I0),II-I0+1,NI,TE,COLZ(K),1)
        ENDDO
        CALL PARROT(CZ(K0),COLZ(K0),KK-K0+1,NK,ZZ,COL(M),1)
        COLR(M)=COL(M)
      ENDDO
C
      RETURN
C
 980  write(*,*)'[Fe/H] =',FEH,
     1  ' OUTSIDE LIMITS OF COLOR CALIBRATION FILE'
      write(*,*)'feh problem'
      STOP
 985  write(*,*)TEFF,', [Fe/H] =',FEH, 
     1  ' OUTSIDE LIMITS OF COLOR CALIBRATION FILE'
      write(*,*)'teff problem'
      STOP
c990  write(*,*)LOGG,', ',TEFF,', [Fe/H] =',FEH,
c    1  ' OUTSIDE LIMITS OF COLOR CALIBRATION FILE'
c     STOP
      END
C
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      FUNCTION DBC(TEFF)
      IMPLICIT REAL*8(A-H,O-Z)
C
C Formula for BC correction for old (published) RYI tables, for Teff > 5600; 
C with double precision arithmetic, correction is accurate to < .01 mag for 
C all Teff up to 20000.
C    Correction is -.01 at 5600 K
C                  -.05 at 7550
C                  -.10 at 9000
C                  -.90 at 20000
      REAL*8 C(4)
      T = TEFF
      IF((4.32-T)*(T-4.14).GE.0.) THEN
        C(1) = 2425.712
        C(2) = -1669.781
        C(3) = 383.1679
        C(4) = -29.32235
       ELSE IF((4.14-T)*(T-3.724).GE.0.) THEN
        C(1) = 528.387
        C(2) = -416.602
        C(3) = 109.5516
        C(4) = -9.60872
       ELSE
        C(1) = 0.
        C(2) = 0.
        C(3) = 0.
        C(4) = 0.
      ENDIF
      DBC = C(1) + C(2)*T + C(3)*T*T + C(4)*T*T*T
      RETURN
      END
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE TRACKIN(SML,Y,FEH,AFE,Z,ZEFF,OV,
     $                  NEEP,NDIVI,LASTPT,NTAUP)
      PARAMETER(NMODMAX=11000,NMAXISO=800)
      PARAMETER (MAXIN=30, NMAXAGE=60, NEEPMAX=12)
C
C This subroutine processs isochrone files from YREC3B and puts things into
C the appropriate variables for use by the isochrone construction routines
C It also determines the models numbers for the primary eep's.
C Requires the user to assign the input *.ISO files to units 11, 12, ....
C where the track files go from the lowest to highest mass.  
C If two successive files have the same mass, the program assumes that 
C the second file is a contination of the first.
C Many parameters are read in from the namelist file, unit 95
C
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8  EEP(12),Y,Z,PHE,AGE,ISOAGE(NMAXAGE),SM(MAXIN),
     $        TEFF,LUM,G,YC(MAXIN,NMODMAX),RADIUS,MCORE,YCC,
     $        MHECORE(MAXIN,NMODMAX),SML(MAXIN),TRUEAGE,OV,OVSH
      INTEGER EEPNUM(13),KEEP,MASS(MAXIN),IEEP,
     *        NTRK,TRACK,JJ,MODEL,JJI,MOD,KJ,MODNEG,NUMEEP,NUMAGE,
     $        NEEP(MAXIN,13),NDIVI(MAXIN),LASTPT(MAXIN),QT(12),
     $        NMOD(MAXIN),NEND,endf,i,modlast,NTRKLOW
      LOGICAL LRGB,LAUTOEEP,LKURCOL,LLEJCOL,LWORCOL,LBASCOL,
     $        LPHXCOL,LVCSCOL
      COMMON/PARMIN/REFZ,SOLBOL,ALFA
      COMMON/QTCOM/QT
      COMMON/ISOIN/ISOAGE,NUMAGE,NTRK
      COMMON/ATRACK/ATEFF(MAXIN,NMODMAX),ALOGL(MAXIN,NMODMAX),
     $        AAGE(MAXIN,NMODMAX),ADIST(MAXIN,NMODMAX)
      COMMON/ZAMS/XQQ(MAXIN),YQQ(MAXIN)
      COMMON/COLTAB/LKURCOL,LLEJCOL,LWORCOL,LBASCOL,LPHXCOL,LVCSCOL
      DATA  LAUTOEEP,LKURCOL,LLEJCOL,LWORCOL,LBASCOL,LPHXCOL,LVCSCOL
     $     /.TRUE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE./
      DATA  ALPFE,OVSH/0.0D0,0.0D0/
      DATA  REFZ,SOLBOL/0.0188,4.790D0/
      DATA  REFX/0.7183/
C QT'S ARE NUMBER OF SUBEEPS TO BE PLACED IN EACH PRIMARY EEP INTERVAL
C     DATA QT/47,9,18,35,50,2,40,1,0,0,0,0/
C     DATA QT/50,15,20,20,50,35,10,8,8,8,8,8/
      DATA QT/50,15,20,20,50,35,10,8,8,8,8,0/

c SRB 6/03 read in alpha/fe and convestive overshoot parameters
c from the first line of the *.iso files, rather than from the namelist
c      NAMELIST /DATA/ALPFE,OVSH,NUMEEP,EEP,NTRK,NUMAGE,ISOAGE,
c     *               REFZ,REFX,QT,LAUTOEEP,LKURCOL,LLEJCOL
      NAMELIST /DATA/NUMEEP,EEP,NTRK,NUMAGE,ISOAGE,
     *               REFZ,REFX,QT,LAUTOEEP,LKURCOL,LLEJCOL,
     $               LWORCOL,LBASCOL,LPHXCOL,LVCSCOL,NTRKLOW


C 
C-------------------- NAMELIST VARIABLES  -----------------------------
C ALPFE    :[alpha/Fe] of the isochrones
C OVSH     :overshoot at the base of the convective envelope
C NUMEEP   :NUMBER OF EEPS
C EEP(I)   :THE PRIMARY EEPS -- FIRST 3 ARE IN TERMS OF CENTRAL Y (FOR MS)
C           THE REST ARE IN TERMS OF MASS OF THE He CORE (RGB)
C LAUTOEEP :TURN ON AUTOMATIC DETERMINATION OF THE EEPS
C NTRK     :TOTAL NUMBER OF TRACK FILES 
C NUMAGE   :NUMBER OF AGES TO CALCULATE ISOCHRONES FOR
C ISOAGE(I):AGE OF ISOCHRONES, IN Myr
C REFZ     :Z OF THE SUN
C SOLBOL   :BOLOMETRIC LUMINOSITY OF THE SUN
C QT(I)    :NUMBER OF SECONDARY EEPs IN THE ITH EEP INTERVAL
C LKURCOL  :use the Kurucz color tables (unit=32), default is the RYI table
C LLEJCOL  : use the Lejeune color table
C LWORCOL  :use the Worthy & Lee 2006 color table
C LBASCOL  :use the BaSeL v3.1 color table
C
C NTRKLOW  : not used by this program, used to lowmass.f which replaces
C            low mass isochrone points
C-------------------- VARIABLE DICTIONARY -----------------------------
C AGE     :AGE (YRS)
C MAGE    :AGE (MYRS) 
C TEFF    : TEFF, IN K
C LUM     :LUMINOSITY (ERG/S)
C YC      :CENTRAL VALUE OF Y  (HE COMP.)
C MHECORE :MASS OF HE CORE  (IN SOLAR MASS UNITS)
C
C output variables used by the isochrone construction part of the code
C FEH     :[Fe/H] = LOG10(Z/REFZ)
C ALPHA   :MIXING LENGTH PARAMETER USED IN THE EVOLUTIONARY CALCULATIONS
C NEEP(I,J): MODEL NUMBER ON THE ITH INPUT TRACK CORRECSPONDING TO 
C            PRIMARY EEP J
C NDIVI(I):NUMBER OF PRIMARY EEP INTERVALS FOR EACH MASS I
C NTAUP   :LARGEST NUMBER OF EEPs  
C LASTPT(I): LAST POINT DEFINDED IN THE ITH MASS TRACK ONCE CONVERTED TO
C            EEPS
C SML(I)  :LOG10(MASS) (in solar units) for THE ITH TRACK
C following variables are defined in TRANEEP
C XX(I,J) :LOG Teff FOR THE JTH EEP POINT ON THE ITH TRACK
C YY(I,J) :LOG L/Lsun FOR THE JTH EEP POINT ON THE ITH TRACK
C TL(I,J) :LOG AGE (Myr) FOR THE JTH EEP POINT ON THE ITH TRACK
C XQQ(I)  :LOG Teff for ZAMSpoint on the Ith track
C YQQ(I)  :Log L/Lsun for ZAMS point on the Ith track
C----------------------------------------------------------------------
C
      OPEN(UNIT=95,STATUS='OLD')
      READ(UNIT=95, NML=DATA)
      CLOSE(95)
c
c BCC 1/2007  -- get colour table selection by reading it in from 
c                a MC var file
      open(unit=75,status='old')
      do i = 1,20
         read(75,*)
      enddo
      read(75,337)itable
 337  format(i2)
      close(75)
      if(itable.eq.0) then
         LPHXCOL=.TRUE.
      elseif(itable.eq.1) then
         LVCSCOL=.TRUE.
      else
         write(*,*)'incorrect colour table in var file'
         stop
      endif

c      AFE = ALPFE
c      OV  = OVSH
C set up luminosity and mass of the sun, in units of erg/s & grams
      CLSUN = 3.8515D33
      CMSUN = 1.9891D33
C set up some constants
      JJ = 1
      LRGB=.FALSE.
      NEND = NTRK
      DO 300 JJJ = 1,NEND
         TRACK=10+JJJ

c SRB 6/03
c check to make sure the *.iso file contains a good evolutionary track
c (many tracks fail during rescaling or around the MS turnoff)
c If file has less than 1000 models, we consider it bad. Decrement NTRK
c and proceed to the next file.

         i = 0
 10      read (TRACK, fmt=*, iostat=endf)
         if (endf .eq. 0) then
            i = i+1
            go to 10
         else if (endf.gt.0 .or. (endf.lt.0 .and. i.lt.1000)) then
c            write(*,*) 'input file', TRACK, 'is bad.'
            close(TRACK)
            NTRK = NTRK-1
            go to 290
         else if (endf.lt.0 .and. i.gt.1000) then
            rewind(TRACK)
         end if

C read in the header, and determine the mass
         CALL HEADIN(SM,CMSUN,JJ,TRACK,YMOD,ZMOD,IMASS,ALPFE,
     $        OVSH,KTTAU,FEH)
         AFE = ALPFE
         OV = OVSH
c         write(*,*)'overshoot = ',ov
c         write(*,*)'[a/Fe] = ',afe
c         IF (JJJ.EQ.1) THEN 
         IF (JJ.EQ.1) THEN
            Z = ZMOD
            Y = YMOD
            XX = 1.0 - Y - Z
            CALL SETUPS(AFE,Z,ZEFF,NEEPMAX,IEEP,EEPNUM)
c            FEH = DLOG10(ZEFF/REFZ) - DLOG10(XX/REFX)
C PUT IN THE SOLAR BOLOMETRIC LUMINOSITY FOR THE RYI, 
C KURUCZ & Lejeune COLOR TABLES
            SOLBOL= 4.79D0
            IF(LKURCOL)SOLBOL = 4.750D0
            IF(LLEJCOL)SOLBOL = 4.720D0
            IF(LWORCOL)SOLBOL = 4.750D0
            IF(LBASCOL)SOLBOL = 4.720D0
            IF(LPHXCOL)SOLBOL = 4.750D0               
            IF(LVCSCOL)SOLBOL = 4.750D0
         ENDIF
         MASS(JJ) = IMASS
c         write(*,874)imass,afe,feh,z,zeff
c 874     format(I4,2x,' [a/Fe]=',F5.2,' [Fe/H]=',F5.2,' Z=',F7.5,
c     $        ' Zeff=',F7.5)
C check for non-increasing mass, or RGB continuation
         IF (JJ.GT.1) CALL MASSCHK(MASS,JJ,LRGB)
c         WRITE(*,32)JJJ,MASS(JJ)
 32      FORMAT(' ISOFILE NUMBER',I3, ' MASS IS ',I3)
 54      READ(TRACK,*) AGE,TEFF,G,LUM,RADIUS,YCC,Zcore,ZXenv,xLh,xLhe,
     $                 MCORE,xMcoCore
         MODEL=1
         IF(MODEL.LT.1) GOTO 54

c***BC change for new format
c         MCORE = MCORE/CMSUN
C set up model number, phe, age, etc.
         IF(LRGB) THEN
            JJ = JJ-1
            NTRK = NTRK-1
            MOD = MODEL + 1 - MODNEG
         ELSE
            MOD    = 1
            MODNEG = MODEL
            TRUEAGE = AGE/1D6
         ENDIF
         MHECORE(JJ,MOD) = MCORE
         YC(JJ,MOD) = YCC
         AAGE(JJ,MOD)  = AGE/1.D6 - TRUEAGE
         ATEFF(JJ,MOD) = TEFF
         ALOGL(JJ,MOD) = LUM
         ADIST(JJ,MOD) = 0.0D0
C finish reading in the track
         DO 50 JJI = 1,10000
            READ(TRACK,*,END=250) AGE,TEFF,G,LUM,RADIUS,YCC,Zcore,
     $                 ZXenv,xLh,xLhe,MCORE,xMcoCore
            MODEL= MODEL + 1
            MOD  = MODEL+1-MODNEG
            IF(MOD.GT.NMODMAX) THEN
               WRITE(*,*)'TOO MANY POINTS IN THIS TRACK'
               STOP
            ENDIF
            YC(JJ,MOD) = YCC
            MHECORE(JJ,MOD) = MCORE
            AAGE(JJ,MOD)  = AGE/1.D6 - TRUEAGE
            ATEFF(JJ,MOD) = TEFF
            ALOGL(JJ,MOD) = LUM
            DX = ATEFF(JJ,MOD) - ATEFF(JJ,MOD-1)
            DY = ALOGL(JJ,MOD) - ALOGL(JJ,MOD-1)
            ADIST(JJ,MOD) = ADIST(JJ,MOD-1) + DSQRT(DX*DX + DY*DY)
            IF(MOD.EQ. 2) THEN
               XQQ(JJ) = ATEFF(JJ,MOD)
               YQQ(JJ) = ALOGL(JJ,MOD)
            ENDIF
 50      CONTINUE
 250     NMOD(JJ) = MOD
         JJ = JJ + 1
 290  continue
 300  CONTINUE


C now, determine the eep's
      IF(LAUTOEEP) CALL GETEEP(MHECORE,NUMEEP,NTRK,NMOD,YMOD,Z,EEP)
      DO 400 JJ=1,NTRK
C set things up
         PHE  = YC(JJ,1)
         KEEP = 1
         DO 500 JJI = 1,NMOD(JJ)

            CALL EEPCHEK(PHE,KEEP,YC(JJ,JJI),MHECORE(JJ,JJI),
     $                   EEPNUM,JJI-1,EEP)
            IF(KEEP.GE.IEEP) THEN
C finished with the eep's --end this loop
               GOTO 230
            ENDIF

c SRB 6/03 Changed the third EEP to be in terms of helium core mass
c rather than central helium abundance
c            IF(KEEP.LT.4) THEN
            IF(KEEP.LT.3) THEN
               PHE = YC(JJ,JJI)
            ELSE
               PHE = MHECORE(JJ,JJI)
            ENDIF
 500     CONTINUE
         DO 45 KJ = KEEP,IEEP
            EEPNUM(KJ)=-1
 45      CONTINUE
 230     DO 47 KJ = 1,12
            NEEP(JJ,KJ+1) = EEPNUM(KJ)
 47      CONTINUE

c SRB   
c         write(*,'(A5,I1,A4,13I5)') 'NEEP(',JJ,',J)=',
c     *        (NEEP(JJ,III), III=1,13)

 400  CONTINUE
c      DO JJJ=1,NTRK
c         write(33,88)(neep(jjj,kj),kj=2,13)
c 88      format(12(I4,1X))
c      ENDDO
C determine the number of primary eep intervals for each mass, and 
C ensure that number of defined tracks is never equal to 1 for any interval
      CALL INTCHK(NEEP,NDIVI,NTRK)
C find the largest number of EEP's in any mass track
      CALL MAXEEP(NDIVI,NTRK,NTAUP)
C transform track from input points to eeps
      CALL TRANEEP(SM,SML,LASTPT,NDIVI,NEEP,NTRK)
      RETURN
      END 
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE HEADIN(SM,CMSUN,JJ,TRACK,YMOD,ZMOD,IMASS,
     $     ALPHAFE,OVERSH,KTTAU,FEH)
      PARAMETER (MAXIN=30)
C this subroutine reads in the header info, and determines the mass of the
C track as well as Z and alfa, and [Fe/H]
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 SM(MAXIN),SSTANDARD(7)
      INTEGER TRACK,IMASS
      CHARACTER*1 INLINE(8)
      COMMON/PARMIN/REFZ,SOLBOL,ALFA
C
      OPEN(UNIT=TRACK,STATUS='OLD')
c****BC  -- change for the MC input files
c      read in the MC header line
      REWIND(TRACK)
c     READ(TRACK,*)
c     READ(TRACK, 1494) xx, feh, Yprim,
c    *     alphafe, cmixlen, FGRY, FGRZ, KTTAU, oversh,
c    *     (SSTANDARD(I), I=1,7), alexcoef, opalcoef2, talphacoef,
c    *     plascoef, cocoef
c1494 format(7F12.8, I6, 14F12.8)
c      write(*,*)'feh=',feh,'[a/Fe]=',alphafe,'MixL=',cmixlen

C read in the standard track header line
c     READ(TRACK,55)SMASS, XENV,ZENV,AFE,CMIXL
c55   FORMAT(4X,F5.3,3x,E10.4,3x,E10.4,18x,F5.2,6x,F6.4)
c     READ(TRACK,*)
c      write(*,*)'Mass=',SMASS, 'X=',XENV,'Z=',ZENV
c      write(*,*)'AFE=',AFE,'CMIXL=',CMIXL
c
C*******BC -- some more MC changes
c      write(*,*)

c****MY read MC parameters from Varfiles
      OPEN(UNIT=75,STATUS='OLD')
      REWIND(UNIT=75)
      read(unit=75, fmt=*) feh
      read(unit=75, fmt=*) Yprim
      read(unit=75, fmt=*) alphafe
      read(unit=75, fmt=*) cmixlen
      read(unit=75, fmt=*) FGRY
      read(unit=75, fmt=*) FGRZ
      read(unit=75, fmt=*) KTTAU
      read(unit=75, fmt=*) ALPHAE
      read(unit=75, fmt=*) 
      read(unit=75, fmt=*) SSTANDARD(1)
      read(unit=75, fmt=*) SSTANDARD(2)
      read(unit=75, fmt=*) SSTANDARD(3)
      read(unit=75, fmt=*) SSTANDARD(4)
      read(unit=75, fmt=*) SSTANDARD(5)
      read(unit=75, fmt=*) SSTANDARD(6)
      read(unit=75, fmt=*) SSTANDARD(7)
      read(unit=75, fmt=*) alexcoef
      read(unit=75, fmt=*) opalcoef2
      read(unit=75, fmt=*) talphacoef
      read(unit=75, fmt=*) plascoef
      read(unit=75, fmt=*) cocoef
      close(unit=75)
C read in the standard track header line
      READ(TRACK,55)SMASS, XENV,ZENV,AFE,CMIXL
 55   FORMAT(4X,F5.3,3x,E10.4,3x,E10.4,18x,F5.2,6x,F6.4)
      READ(TRACK,*)
C
      YMOD1 = Yprim
      ZMOD1 = ZENV
      XENV1 = XENV
      
C
      IF(JJ.EQ.1) THEN
         YMOD = Yprim
         ALFA = CMIXL
         ZMOD = ZENV
c         write(*,*)'Y:',ymod,ymod1
c         write(*,*)'X:',XENV,XENV1
c         write(*,*)'Z:',ZMOD,ZMOD1
c         write(*,*)'MIX:',ALFA,ALFA1
         zmod = zmod1
         zenv = zmod1
         ymod = ymod1
         xenv = xenv1
      ELSE
C perform some consitency checks
         zmod = zmod1
         zenv = zmod1
         ymod = ymod1
         xenv = xenv1
         Y = 1.0D0-XENV-ZENV
         IF(DABS(Y - YMOD).GT.1.0D-3 .AND. Y.NE.1.0D0) THEN
            WRITE(*,*)'ISOFILE NUMBER ',TRACK-10, ' HAS AN INCORRECT Y'
            WRITE(*,*)'TRACK Y = ',Y,' PREVIOUS Y = ',YMOD
            STOP
         ENDIF
         Z = ZENV
         IF(DABS(ZMOD-Z).GT.1.0D-7 .AND. Z.NE.0.0D0) THEN
            WRITE(*,*)'ISOFILE NUMBER ',TRACK-10, ' HAS AN INCORRECT Z'
            WRITE(*,*)'TRACK Z = ',Z,' FIRST TRACK Z = ',ZMOD
            STOP
         ENDIF
         IF(ABS(CMIXL-ALFA).GT.1.0D-5) THEN
            WRITE(*,*)'ISOFILE NUMBER ',TRACK-10, 
     $                ' HAS AN INCORRECT MIXING LENGTH'
            WRITE(*,*)'TRACK ALFA = ',CMIXL,' PREVIOUS ALFA = ',ALFA
            STOP
         ENDIF
      ENDIF
C convert mass from grams to solar units*100 (integer)
      IMASS = NINT(SMASS*100.0D0)
c      write(*,*)'Mass:',smass,imass
      SM(JJ) = SMASS

c skip over the pre-main sequence models
c***ASSUME no more than 1200 pre-main sequence models
      read(TRACK,*)age,tlog,glog,blog,rlog,ycore1,zcore,zx,
     $       hl,hel,hecore,cocore
      do iii=2,1200
         read(TRACK,*)age,tlog,glog,blog,rlog,ycore,zcore,zx,
     $        hl,hel,hecore,cocore
c         write(*,200)age,ycore,hl
c 200     format('Age=',E10.4,' Ycore=',E12.6,' M_He_core=',E12.6)
         if(ycore.gt.(ycore1 + 0.0025)) then 
c            write(*,*)'end of pre-ms',iii,yprim
            RETURN
         endif
      enddo

      RETURN
      END
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE MASSCHK(MASS,JJ,LRGB)
C This subroutine ensures that input iso files are in order of increasing
C mass.  If the present track has the same mass as the previous track
C we have a rgb continuation track and lrbg=.true. 
      PARAMETER (MAXIN=30)
C
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER MASS(MAXIN),JJ
      LOGICAL LRGB
      LRGB=.FALSE.
      IF(MASS(JJ).LT.MASS(JJ-1)) THEN
         WRITE(*,*)'INPUT ISO FILES NOT IN ORDER OF INCREASEING MASS'
         WRITE(*,100)JJ,MASS(JJ)
 100     FORMAT('ISOFILE NUMBER',I3,' HAS A MASS OF ',I3)
         STOP
      ENDIF
      IF(MASS(JJ).EQ.MASS(JJ-1)) THEN
         LRGB = .TRUE.
         WRITE(*,*)'RGB CONTINUATION FILE'
      ENDIF
      RETURN
      END
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE ISOREAD(CLSUN,CMSUN,TRUEAGE,TRACK,MODNEG,JJ)
      PARAMETER(NMODMAX=11000,NMAXISO=800)
      PARAMETER (MAXIN=30, NMAXAGE=60)
C finish reading a track file, without checking for eep's
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 MHECORE,TEFF,LUM,G,RADIUS,CMSUN,CLSUN
      INTEGER TRACK
      COMMON/ATRACK/ATEFF(MAXIN,NMODMAX),ALOGL(MAXIN,NMODMAX),
     $        AAGE(MAXIN,NMODMAX),ADIST(MAXIN,NMODMAX)
      DO 60 JJI = 1,10000
         READ(TRACK,*,END=999) MODEL,AGE,LUM,RADIUS,TEFF,G,YC,
     *        MHECORE
         MOD  = MODEL+1-MODNEG
         IF(MOD.GT.NMODMAX) THEN
            WRITE(*,*)'TOO MANY POINTS IN THIS TRACK'
            STOP
         ENDIF
         AAGE(JJ,MOD)  = AGE/1.D6 - TRUEAGE
         ATEFF(JJ,MOD) = DLOG10(TEFF)
         ALOGL(JJ,MOD) = DLOG10(LUM/CLSUN)
         DX = ATEFF(JJ,MOD) - ATEFF(JJ,MOD-1)
         DY = ALOGL(JJ,MOD) - ALOGL(JJ,MOD-1)
         ADIST(JJ,MOD) = ADIST(JJ,MOD-1) + DSQRT(DX*DX + DY*DY)
 60   CONTINUE
 999  RETURN
      END
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE SETUPS(AFE,Z,ZEFF,NEEPMAX,IEEP,EEPNUM)
C
      REAL*8 ZEFF,Z,AFE
      INTEGER NEEPMAX,IEEP,EEPNUM(13)
C
      IEEP = NEEPMAX + 1
 200  DO 110 JJ=IEEP,12
         EEPNUM(JJ)=-1
 110  CONTINUE
c calculate ZEFF using relation of Salris, Chieffi & Straniero
C (1993, ApJ, 414, 580
c      WRITE(*,*)'IN SETUPS, afe=',afe
      ZEFF = Z/(0.638D0*(10.0D0**AFE) + 0.362D0)
      RETURN
      END
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE EEPCHEK(PHE,KEEP,YC,MHECORE,EEPNUM,MODEL,EEP)
C
C This subroutine checks to see if the previous model was an EEP
C If so, it stuffs the model number into EEPNUM, and increases KEEP by 1

c SRB 6/03 changed the third EEP from being in terms of central
c helium abundance to being in terms of helum core mass. GETEEP
c also incorporates the change

C
      REAL*8 PHE,YC,MHECORE,HE,EEP(12),X,Y
      INTEGER EEPNUM(12),KEEP,MODEL
C
c      IF (KEEP.LT.4) THEN
      IF (KEEP.LT.3) THEN
         HE=YC
      ELSE
         HE=MHECORE
      ENDIF
      X=ABS(PHE - EEP(KEEP))
      Y=ABS(HE - EEP(KEEP))
      IF(X.LT.Y.AND.X.LT.0.04) THEN
         EEPNUM(KEEP) = MODEL
         KEEP=KEEP+1
      ENDIF
      RETURN
      END
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE INTCHK(NEEP,NDIVI,NTRK)
      PARAMETER (MAXIN=30)
C This subroutine determines the number of primary eep intervals for each 
C mass, and ensures that the number of defined tracks is never equal
C to 1 for any interval
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER NEEP(MAXIN,13),NDIVI(MAXIN)
C
      DO 6 I=1,NTRK
         DO 5 J=2,13
            IF (NEEP(I,J).LT.0) THEN
               NDIVI(I) = J-2
               GOTO 6
            ENDIF
 5       CONTINUE
         NDIVI(I) = 12
 6    CONTINUE
C CHECK THAT NUMBER OF DEFINED TRACKS IS NEVER EQUAL TO 1 FOR ANY INTERV
      DO 415 J=2,13
         NSUM=0
         DO 405 I=1,NTRK
            IF(NEEP(I,J).GT.0) THEN
               NSUM=NSUM+NEEP(I,J)
               IM=I
            ENDIF
 405     CONTINUE
         IF(NSUM.EQ.0) GOTO 425
         IF(NEEP(IM,J).EQ.NSUM) GO TO 420
 415  CONTINUE
      J=NDIVI(IM)+2
 420  NDIVI(IM)=J-2
 425  CONTINUE


C SRB
c      write(*,'(A9,9I4)') 'NDIVI=', (NDIVI(I), I=1, NTRK)

      RETURN
      END
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE MAXEEP(NDIVI,NTRK,NTAUP)
      PARAMETER (MAXIN=30)
C This subroutine determines the maximum number of eeps in the mass tracks
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER NDIVI(MAXIN)
      INTEGER QT(12),QSUM
      COMMON/QTCOM/QT
C
      NTAUP = 1
      DO 100 I=1,NTRK
C QSUM IS TOTAL NUMBER OF EEPS IN THIS MASS TRACK; NTAUP IS LARGEST QSUM
         QSUM=1
         DO 200 J=1,NDIVI(I)
            QSUM=QSUM+QT(J)
 200     CONTINUE
         IF (QSUM.GT.NTAUP) NTAUP=QSUM
 100  CONTINUE
      RETURN
      END
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE TRANEEP(SM,SML,LASTPT,NDIVI,NEEP,NTRK)
      PARAMETER(NMODMAX=11000,NMAXISO=800)
      PARAMETER (MAXIN=30)
C This subroutine transforms the track files from the input points to EEPs
C
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 SM(MAXIN),SML(MAXIN),DIST(NMODMAX),X(NMODMAX),
     $       Y(NMODMAX),Z(NMODMAX)
      INTEGER NDIVI(MAXIN),NEEP(MAXIN,13),QT(12),LASTPT(MAXIN),
     $        FIRST,LAST
      COMMON/EEPTRK/XX(NMODMAX,MAXIN),YY(NMODMAX,MAXIN),
     $              TL(NMODMAX,MAXIN)
      COMMON/ATRACK/ATEFF(MAXIN,NMODMAX),ALOGL(MAXIN,NMODMAX),
     $        AAGE(MAXIN,NMODMAX),ADIST(MAXIN,NMODMAX)
      COMMON/QTCOM/QT
      
C
C TRANSFORM TRACK FROM INPUT POINTS TO EEPS; KAPPA=1,NTAUP IS INDEX
C NUMBER OF EEPS AFTER DIVIDING UP PRIMARY INTERVALS ACC. TO QT'S
      DO 100 I=1,NTRK
         DO 200 IJK=1,NMODMAX
            X(IJK) = ATEFF(I,IJK)
            Y(IJK) = ALOGL(I,IJK)
            Z(IJK) = DLOG10(AAGE(I,IJK)+1.D0)
 200     CONTINUE
         KAPPA=1
C FIRST EEP IS SET EQUAL TO ZERO AGE MAIN SEQUENCE
         NEEP(I,1)=1
         XX(1,I)=X(1)
         YY(1,I)=Y(1)
         TL(1,I)=Z(1)
C
         DO 61 J=1,NDIVI(I)
            NPTA=NEEP(I,J)
            NPTB=NEEP(I,J+1)
            IF(NPTA.EQ.NPTB) THEN
               WRITE(*,*)'we have equal eeps!!'
               STOP
            ENDIF
            IF(QT(J).GT.1) THEN
               FIRST=MAX0(NPTA-1,1)
               LAST=MIN0(NPTB+1,NEEP(I,NDIVI(I)+1))
               IF(LAST.LE.0) PRINT *,NEEP(I,J+1),NEEP(I,NDIVI(I)+1)
               IF(LAST.LE.0) PRINT *,'LAST=',LAST
C COMPUTE ARC LENGTH BETWEEN PRIMARY EEPS ON A EACH TRACK;
C CHOOSE METRIC SO THAT ARC LENGTH AND TIME ALONG TRACK ARE
C EQUALLY WEIGHTED OVER EACH SEGMENT
               DMETR=ADIST(I,NPTB)-ADIST(I,NPTA)
               TMETR=Z(NPTB)-Z(NPTA)
               DIST(FIRST)=0.0D0
               DO 25 IN=FIRST+1,LAST
                  IF(IN.GT.NMODMAX) THEN
                     WRITE(*,*)'uh?'
                     STOP
                  ENDIF
                  DD=(ADIST(I,IN)-ADIST(I,IN-1))/DMETR
                  DT=(Z(IN)-Z(IN-1))/TMETR
                  DIST(IN)=DIST(IN-1)+DSQRT(DD*DD + DT*DT)
                  IF(DIST(IN).EQ.DIST(IN-1)) DIST(IN-1)=
     1                           DIST(IN-1)-.01*(DIST(IN)-DIST(IN-2))
 25            CONTINUE
C DIVIDE METRIC 'DISTANCE' BETWEEN EEPS INTO EQUAL SEGMENTS BY
C INTERPOLATION BETWEEN ORIGINAL TRACK POINTS
c               NPTA=NEEP(I,J)
               NPTC=NPTA+1
               ALPHA=(DIST(NPTB)-DIST(NPTA))/QT(J)
               KAPPA=KAPPA+1
               DIST0=DIST(NPTA)
               K=1
               DO 501 IK=KAPPA,KAPPA+QT(J)-2
                  DIST0=DIST0+ALPHA
 502              IF(DIST(NPTC).LE.DIST0) THEN
                     NPTC=NPTC+1
                     K=1
                     GO TO 502
                  ENDIF
                  IS=1
                  CALL PARROT(DIST,X,NPTC,LAST,DIST0,XX(IK,I),K)
                  IS=2
                  CALL PARROT(DIST,Y,NPTC,LAST,DIST0,YY(IK,I),K)
                  IS=3
                  CALL PARROT(DIST,Z,NPTC,LAST,DIST0,TL(IK,I),K)
                  K=2
 501           CONTINUE
               KAPPA=KAPPA+QT(J)-2
            ENDIF
C SET LAST EEP IN INTERVAL TO LAST TRACK POINT IN INTERVAL
            KAPPA=KAPPA+1
            XX(KAPPA,I)=X(NPTB)
            YY(KAPPA,I)=Y(NPTB)
            TL(KAPPA,I)=Z(NPTB)
 61      CONTINUE
C SAVE INDEX OF LAST DEFINED TRACK POINT
         LASTPT(I)=KAPPA
         SML(I)=DLOG10(SM(I))
 100  CONTINUE

C SRB
      open(unit=60, file='eeptracks', status='unknown')

      do I = 1, NTRK
         write(60,'(A,I3,3A15)') '#', I,'EEP#','LogTeff','LogL/Lsun'
         do J = 1, LASTPT(I)
            if (XX(J,I) .gt. 0) then
            write(60,'(4X,I15,2F15.4)') J, XX(J,I), YY(J,I)
            end if
         end do
         write(60,*)
         write(60,*)
      end do

      close(60)

      RETURN
      END
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE TRKDEF(SML,TAUL,J,NQB,QTMS,LASTPT,L,QSUM,NTAUP,
     $                  IT,IQ,SQ,LLARGE,LEND)
      PARAMETER(NMODMAX=11000,NMAXISO=800)
      PARAMETER (MAXIN=30, NMAXAGE=60)
C This subroutine checks to which tracks are defined in a given eep
C
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 ISOAGE(NMAXAGE),SQ(MAXIN),SML(MAXIN)
      INTEGER QT(12),QTMS,QSUM,IQ(MAXIN),LASTPT(MAXIN)
      LOGICAL LEND,LLARGE
      COMMON/EEPTRK/XX(NMODMAX,MAXIN),YY(NMODMAX,MAXIN),
     $              TL(NMODMAX,MAXIN)
      COMMON/QTCOM/QT
      COMMON/ISOIN/ISOAGE,NUMAGE,NTRK
C
      LEND=.FALSE.
      LLARGE=.FALSE.
      L=L+1
      QSUM=QSUM+QT(L)
c      IF(L.GT.NQB) QSUM=NTAUP+1
      IT=0
      DO 72 I=1,NTRK
         IF(J.LE.LASTPT(I)) THEN
            IT=IT+1
            IQ(IT)=I
            SQ(IT)=SML(I)
         ENDIF
 72   CONTINUE
C IF NO TRACKS ARE DEFINED FOR THIS EEP, END LOOP
      IF(IT.EQ.0) THEN
         LLARGE=.TRUE.
      ELSE
C
C CHECK LAST SUBEEP IN THIS INTERVAL TO SEE IF HIGH MASS EXTRAPOLATION
C WILL BE NEEDED.  ISOCHRONE IS COMPUTED IN WHOLE INTERVAL SEGMENTS,
C FOR CONSISTENCY.  IF LAST SUBEEP (#IN) IS OUTSIDE AGE LIMIT, END LOOP
         IN=QSUM-1
         IF(TAUL.LT.TL(IN,IQ(IT))) THEN
C HIGH MASS EXTRAPOLATION BEFORE END OF MS NOT ALLOWED; INPUT AGE
C TOO LOW FOR GOOD ISOCHRONE.
            IF(J.LE.QTMS) THEN
               WRITE(*,*)'ISOCHRONE FAILED TURNOFF MASS TOO HIGH'
               LEND=.TRUE.
            ELSE
C ALLOWING HIGH MASS EXTRAPOLATION ON RGB SHOULD ONLY BE DONE WITH
C EXTREME CAUTION: IN GENERAL, IT PRODUCES POOR RESULTS AT EEP QTB+1
C (WHERE INTERPOLATION GOES FROM AVERAGED QUADRATIC TO LINEAR)
C AND OFTEN NEAR HOOK.  ALSO, ALLOWED HERE ONLY IF NEXT LOWER MASS
C TRACK EXISTS.
               IF(IQ(IT)-1.NE.IQ(IT-1).OR.TAUL.LT.TL(IN,IQ(IT))-
     $              .2*(TL(IN,IQ(IT)-1)-TL(IN,IQ(IT)))) LLARGE=.TRUE.

C SRB 3/03
c               IF(IQ(IT)-1.NE.IQ(IT-1).OR.TAUL.LT.TL(IN,IQ(IT))-
c     $              .2*(TL(IN,IQ(IT)-1)-TL(IN,IQ(IT)))) THEN
c                  LLARGE=.TRUE.
c                  write(*,*) 'high mass extrapolation disallowed?'
c               end if


            ENDIF
         ENDIF
      ENDIF
      RETURN
      END
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE MASSK(TQ,TAUL,J,IT,SMT,K,LINEAR,LSMALL,LEND)
      PARAMETER(NMODMAX=11000,NMAXISO=800)
      PARAMETER (MAXIN=30, NMAXAGE=60)
C This subroutine checks low mass extrapolation, and determines K such
C that TQ(K) < TAUL < TQ(K-1)
C
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 TQ(MAXIN),SMT(NMODMAX)
      INTEGER QT(12)
      LOGICAL LINEAR,LSMALL,LEND
      COMMON/QTCOM/QT
      LEND=.FALSE.
      LSMALL=LEND
      LINEAR=LEND
C SET K=2 BEFORE BRANCH
      K=2
C CHECK EVERY J TO SEE IF LOW MASS EXTRAPOLATION NEEDED
      IF(TAUL.GT.TQ(1)) THEN
C LOW MASS EXTRAPOLATION ALLOWED ONLY IN FIRST INTERVAL (LOWER MS)
C AND ON RGB; OTHERWISE INPUT AGE TOO LARGE FOR GOOD ISOCHRONE
c         IF(J.EQ.QT(1)+1) THEN
c            WRITE(*,*) 'ISOCHRONE FAILED: MS MASS TOO LOW'
c            LEND=.TRUE.
c            RETURN
c         ENDIF
         IF(TAUL.LE.TQ(1)+2.*(TQ(1)-TQ(2))) THEN
            LINEAR=.TRUE.
            RETURN
         ENDIF
C MASS TOO SMALL
         SMT(J)=-1.
         LSMALL=.TRUE.
         RETURN
      ENDIF
      IF(IT.EQ.2) THEN
         LINEAR=.TRUE.
      ELSE
C FIND K SUCH THAT TQ(K) < TAUL < TQ(K-1)
 80      IF (TAUL.LE.TQ(K)) THEN
            IF(K.EQ.IT) THEN
               LINEAR=.TRUE.
               RETURN
            ENDIF
            K=K+1
            GO TO 80
         ENDIF
c         DO 86 I=2,IT
c            IF(TQ(I).GE.TQ(I-1)) THEN
c               WRITE(*,880) J
c 880           FORMAT(//1X,'THERE IS A PROBLEM WITH CHOICE OF ',
c     $              'EEPS OR METRIC:',/,6X,'AT EEP ',I3,', TIME ',
c     $              'IS NOT A MONOTONIC FUNCTION OF MASS')
c            ENDIF
c 86      CONTINUE
      ENDIF
      RETURN
      END
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE NUMSTAR(SP1,SMT,ANUM,NT,N0,N2)
      PARAMETER(NMODMAX=11000,NMAXISO=800)
      PARAMETER (MAXIN=30, NMAXAGE=60)
C This subroutine computes the number of stars (ANUM(J,M)) in the 
C interval between EEP J and EEP J+1
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 SP1(3),SMT(NMODMAX),ANUM(NMAXISO,3)
      LOGICAL LLOW,LHIGH
C
C Do a check, and make sure that all points in the range N0 to N2 really
C are defined
      LLOW=.FALSE.
      LHIGH=.FALSE.
      DO 10 J=N0+1,N2+1
         IF(SMT(J).LE.0.0) THEN
C something went wrong earlier, and we have to change n0 or n2
            IF(J.LT.( (N0+N2)/2 ) ) THEN
               ILOW = J +1
               LLOW=.TRUE.
            ELSE
               IF(.NOT.LHIGH) THEN
                  IHIGH=J - 2
                  LHIGH=.TRUE.
               ENDIF
            ENDIF
         ENDIF
 10   CONTINUE
C change N0,N2 if necessary
      IF(LLOW)  N0=ILOW
      IF(LHIGH) N2=IHIGH
      DO 99 M=1,3
C SFACTR NORMALIZES TO 1000 STARS IN THE MS MASS RANGE 0.5 TO 1.0
         SFACTR=1000./(1.-0.5**SP1(M))
         SMJ=SMT(N0)**SP1(M)
         DO 97 J=N0,N2
            SMNXT=SMT(J+1)**SP1(M)
            ANUM(J,M)=(SMNXT-SMJ)*SFACTR
            SMJ=SMNXT
 97      CONTINUE
 99   CONTINUE
      DO 98 M=1,3
         ANUM(NT,M)=0.
 98   CONTINUE
      RETURN
      END
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE OUTISO(XXT,YYT,SMT,BOL,ANUM,Y,FEH,AFE,Z,ZEFF,OV,TAU,
     $               IGB,NPTS,N0,NTAUP,AX)
      IMPLICIT NONE
      INTEGER NMODMAX,NMAXISO,MAXIN,NMAXAGE,itable
      PARAMETER(NMODMAX=11000,NMAXISO=800)
      PARAMETER (MAXIN=30, NMAXAGE=60)
C This subroutine outputs the isochrones into file number 32
      REAL*8 XXT(NMODMAX),YYT(NMODMAX),SMT(NMODMAX),BOL(NMODMAX),
     $    ANUM(NMAXISO,3),ALFA,FEH,TAU,AFE,Z,ZEFF,OV,Y,afe1,filter(13)
      REAL*8 ASUM1,ASUM2,ASUM3,R2RL,LOGG,V,COLR(8),CLN,TEMP,REFZ,
     $       SOLBOL,DBC
      INTEGER IGB,NPTS,N0,NTAUP,J,ICOL,NCOLOR,rnpts
      LOGICAL LKURCOL,LLEJCOL,LWORCOL,LBASCOL,LPHXCOL,LVCSCOL
      logical lfirst
      CHARACTER*1 AX(NMAXISO)
      COMMON/PARMIN/REFZ,SOLBOL,ALFA
      COMMON/COLTAB/LKURCOL,LLEJCOL,LWORCOL,LBASCOL,LPHXCOL,LVCSCOL
      save lfirst
      data lfirst/.true./
C THESE 3 LINES PRECEDE EACH ISOCHRONE IN OUTPUT FILE
      WRITE(92,119)

c SRB 6/03 only output EEPs number 30 or above. The number of points,
c then, may be reduced
c      WRITE(92,120) NPTS,ALFA,OV,TAU,Y,Z,ZEF,FEH,AFE
      rnpts = min(NPTS,NPTS-(30-N0))
      WRITE(92,120) rnpts,ALFA,OV,TAU,Y,Z,ZEFF,FEH,AFE
 119  FORMAT('#NPTS MIX-LEN  OVERSH AGE(MYR)   Y        Z     ',
     $       '    ZEFF    [Fe/H] [alpha/Fe]')
 120  FORMAT('#',I3,F9.6,F8.4,F7.0,F8.4,2(E12.4),2(F6.2))
c SRB 6/03 end changes

C********************BC
C big hack here -- NEED TO change the LEJCOL routine so that it
C returns the correct IR colors....  Right now, V-J, V-H and V-K are not 
C being outputed....
c      WRITE(*,111) TAU
 111  FORMAT(' AGE = ',F6.0,' Myr')
      ASUM1=0.0D0
      ASUM2=ASUM1
      ASUM3=ASUM1
      NCOLOR=5
      IF(LKURCOL) NCOLOR=8
      IF(LLEJCOL) NCOLOR=8
      IF(LWORCOL) NCOLOR=8
      IF(LBASCOL) NCOLOR=8
      IF(LPHXCOL) THEN 
         NCOLOR=4
         itable = 4
c [a/Fe] in 0.2 dex steps in the tables...
         if(afe.lt.-0.10d0) afe = -0.2d0
         if(afe.ge.-0.10d0 .and. afe.lt.0.10) afe1 = 0.0d0
         if(afe.ge.0.10d0 .and. afe.lt.0.30) afe1 = 0.2d0
         if(afe.ge.0.30d0 .and. afe.lt.0.50) afe1 = 0.4d0
         if(afe.ge.0.50d0 .and. afe.lt.0.70) afe1 = 0.6d0
         if(afe.ge.0.70) afe1 = 0.8d0
         if(lfirst) call color_init(itable,feh,afe1)
         itable = 5
         if(lfirst) call color_init(itable,feh,afe1)
      endif
      IF(LVCSCOL) THEN
         NCOLOR=4
         itable = 8
         if(lfirst) call color_init(itable,feh,afe)      
      endif
      lfirst=.false.
      IF(NCOLOR.GT.4) WRITE(92,121)
      IF(NCOLOR.EQ.4) WRITE(92,144)
 121  FORMAT('#EEP  MASS     Log G    Log Teff   Log(L/Ls)    NUM0',
     $       '           NUM-2.35        NUM-4.0         V      U-B',
     $       '    B-V    V-R    V-I    V-J    V-H    V-K')
 144  FORMAT('#EEP  MASS     Log G    Log Teff   Log(L/Ls)    NUM0',
     $       '           NUM-2.35        NUM-4.0         V      V-I',
     $       '    F606W  F606W-F814W')

      DO 90 J=N0-1,NTAUP
         IF (SMT(J).GT.0.0D0) THEN
C COLOR ISOCHRONES ARE CONSTRUCTED
            R2RL=-4.D0*XXT(J)+YYT(J)+2.D0*1.08426D1+4.D0*(3.763D0)
            LOGG=-7.17600D0+DLOG10(SMT(J))+3.329863D1-R2RL
            CLN=DLOG(10.0D0)
            TEMP=DEXP(XXT(J)*CLN)
c            write(*,*)temp,xxt(j)
C     
            IF(LKURCOL) THEN
               CALL KURCOL(FEH,TEMP,LOGG,COLR)
            ELSEIF(LLEJCOL) THEN
               CALL LEJCOL(FEH,TEMP,LOGG,COLR)
            ELSEIF(LWORCOL) THEN
               CALL WORCOL(FEH,TEMP,LOGG,COLR)
            ELSEIF(LBASCOL) THEN
               CALL BASCOL(FEH,TEMP,LOGG,COLR)
            elseif(lphxcol) then
               itable=4
               call get_mags(itable,feh,afe1,bol(j),logg,xxt(j),filter)
               v = filter(3)
               colr(2) = filter(3) - filter(5)
               itable = 5
               call get_mags(itable,feh,afe1,bol(j),logg,xxt(j),filter)
               colr(3) = filter(4)
               colr(4) = filter(4) - filter(7)
            elseif(lvcscol)then
               itable = 8
               call get_mags(itable,feh,afe,bol(j),logg,xxt(j),filter)
               v = filter(2)
               colr(2) = filter(2) - filter(3)
               colr(3) = filter(4)
               colr(4) = filter(4) - filter(5)
            ELSE
               CALL COLOR(FEH,TEMP,LOGG,COLR)
            ENDIF
            IF(COLR(2).GT.-9.999D0 .AND. 
     $           .NOT.LPHXCOL .AND. .NOT. LVCSCOL) THEN
               V=BOL(J)-COLR(1)
C correct the bolometric correction
               IF(.NOT.LKURCOL .AND. .NOT. LLEJCOL .AND. 
     $            .NOT.LWORCOL .AND. .NOT. LBASCOL) V=V - DBC(XXT(J))
            ELSEIF(.NOT.LPHXCOL .AND. .NOT. LVCSCOL) THEN
               V = -9.999D0
            ENDIF
C
            ASUM1=ANUM(J,1)+ASUM1
            ASUM2=ANUM(J,2)+ASUM2
            ASUM3=ANUM(J,3)+ASUM3
c            write(*,*)'ncolor=',ncolor
c            write(*,*)'ir cols:',colr(6),colr(7),colr(8)

c SRB 6/03 luminosity function program fails if there are more than 205
c points in the input isochrone. So don't write out EEPS less than 30.
            if (J .ge. 30) then

            WRITE(92,151) J,AX(J),SMT(J),LOGG,
     $           XXT(J),YYT(J),ANUM(J,1),ANUM(J,2),
     $           ANUM(J,3),V,(COLR(ICOL),ICOL=2,NCOLOR)
 151        FORMAT (I3,A1,F9.6,2F10.6,F11.7,3E16.9E1,8F7.3)

            end if
c SRB 6/03 end changes

         ELSE
            IF(SMT(J).EQ.0.0D0) THEN
 91            WRITE(*,*)'MASS IS TOO LARGE OR TRACKS NOT DEFINED '
     $              ,'IN NEXT INTERVAL'
               write(92,*)
               write(92,*)
               RETURN
            ENDIF
         ENDIF
 90   CONTINUE
c  SRB 2/03 added a couple of blank lines between ages so gnuplot will interpret
c  them as separate data sets
      write(92,*)
      write(92,*)
      RETURN
      END
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE GETEEP(MHECORE,NUMEEP,NTRK,NMOD,YMOD,Z,EEP)
C This routine automatically determines the EEPs.  The QT (number of 
C points in each EEP interval) are also rearranged for optimal isochrone
C construction.  This routine only works for input tracks which have good
C sized giant brances
C 
      IMPLICIT NONE
      INTEGER NMODMAX,NMAXISO,MAXIN,NMAXAGE
      REAL*8 EPS
      PARAMETER(NMODMAX=11000,NMAXISO=800)
      PARAMETER (MAXIN=30, NMAXAGE=60)
      PARAMETER (EPS = 1.0D-4)
C
      REAL*8 MHECORE(MAXIN,NMODMAX),EEP(12),YMOD,Z
      INTEGER NUMEEP,NTRK,NMOD(MAXIN),QT(12)
      COMMON/QTCOM/QT
C
      REAL*8 YMIN,YMAX,YADD,MAXCORE(MAXIN),ENDCORE,COREADD
      INTEGER I
C default the number of EEPs to 12
      NUMEEP = 12
C first three eep's are on the main sequence in terms of central He abundance
C Put EEP(3) at the main sequence turn-off, the other two are equally spaced.
      YMIN = YMOD
      YMAX = 1.0D0 - Z - EPS

c SRB 6/03 
c With heavy element diffusion, some stars never reach a central helium
c abundance equal to YMAX = 1 - Z - EPS. Therefore we changed the third EEP
c to be in terms of the helium core mass (a very small value will be near
c the main sequence turnoff) rather than central helium abundance. Subroutine
c EEPCHK also incorporates the change.

c      EEP(3) = YMAX
      EEP(3) = 0.01D0

      YADD = (YMAX - YMIN)/3.0D0
      IF(YADD.LT.0.0) THEN
         WRITE(*,*)'something wrong in geteep'
         STOP
      ENDIF
      EEP(1) = YMIN + YADD
      EEP(2) = EEP(1) + YADD
C now do the red giant branch.  EEP's are in terms of the mass of the
C He core. for the rgb, use default values for the 3.  Next N values, are 
C equally spaced, to give full coverage of the input RGB tracks
      EEP(4) = 0.08D0
      EEP(5) = 0.16D0
c      EEP(6) = 0.24D0
      EEP(6) = 0.235D0
C find the 2nd highest He core max attained -- this sets the upper
C limit to the EEPs
C first, put the maxiumum He core max attained into their own array
      DO 100 I=1,NTRK
         MAXCORE(I) = MHECORE(I,NMOD(I))
 100  CONTINUE
C do a shell sort, to put things into ascending order
      CALL SHELL(NTRK,MAXCORE)
      ENDCORE = MAXCORE(NTRK-1) - EPS
      IF(ENDCORE.LT.0.26) THEN
         WRITE(*,*)'INPUT FILES HAVE SMALL GIANT BRANCHES'
         IF(ENDCORE.LT.0.02) THEN
            WRITE(*,*)'NO GIANT BRANCH WILL BE CALCULATED'
            NUMEEP=3
         ELSE
            NUMEEP=6
            EEP(6) = ENDCORE
            EEP(4) = ENDCORE/3.0D0
            EEP(5) = 2.0D0*EEP(4)
         ENDIF
      ELSE
         EEP(NUMEEP) = ENDCORE
c         COREADD = (ENDCORE - 0.24)/(NUMEEP - 6) 
         COREADD = (ENDCORE - 0.235)/(NUMEEP - 6)
         DO 200 I=7,NUMEEP
            EEP(I) = EEP(I-1) + COREADD
 200     CONTINUE
      ENDIF
      DO 300 i=1,numeep
c         write(*,*)'i=',i,' eep = ',eep(i)
 300  continue
      RETURN
      END
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE ZAMSISO(SQ,SP1,Y,FEH,AFE,Z,ZEFF,OV)
C This routine calculates a ZAMS isochrone
      IMPLICIT NONE
      INTEGER NMODMAX,NMAXISO,MAXIN,NMAXAGE
      PARAMETER(NMODMAX=11000,NMAXISO=800)
      PARAMETER (MAXIN=30, NMAXAGE=60)
C
      REAL*8  ISOAGE(NMAXAGE),SQ(MAXIN),SP1(3),REFZ,SOLBOL,ALFA,
     $     Y,FEH,AFE,Z,ZEFF,OV,XQQ,YQQ
      INTEGER QT(12),NUMAGE,NTRK
      COMMON/QTCOM/QT
      COMMON/ISOIN/ISOAGE,NUMAGE,NTRK
      COMMON/PARMIN/REFZ,SOLBOL,ALFA
      COMMON/ZAMS/XQQ(MAXIN),YQQ(MAXIN)
C
      REAL*8 SM(MAXIN),XM(100),XML(100),XQ(MAXIN),YQ(MAXIN),
     $     TEFFL(100),LOGL(100),BOL(100),ANUM(NMAXISO,3),DIFFM
      INTEGER I,K,NMAX,N0,N2,NPTS,NT
      CHARACTER*1 AX(100)
C
      NMAX=QT(1)
      DO 82 I=1,NTRK
         XQ(I)=XQQ(I)
         YQ(I)=YQQ(I)
         SM(I)=10.0D0**SQ(I)
 82   CONTINUE
      DIFFM = (SM(NTRK) - SM(1))/(NMAX-1)
c      write(*,*)nmax,diffm
      XM(1) = SM(1)
c      write(*,*)'xm(1)=',xm(1)
      AX(1) = ' '
      DO 100 I=2,NMAX
         XM(I) = XM(I-1) + DIFFM
         XML(I) = DLOG10(XM(I))
         AX(I) = AX(1)
c         write(*,*)i,xm(i)
 100  CONTINUE
      TEFFL(1)   = XQ(1)
      LOGL(1)    = YQ(1)
      BOL(1)     = SOLBOL-(2.5D0*LOGL(1))
      TEFFL(NMAX)= XQ(NTRK)
      LOGL(NMAX) = YQ(NTRK)
      BOL(NMAX)  = SOLBOL-(2.5D0*LOGL(NMAX))
      K = 2
      DO 200 I=2,NMAX-1
C find K such that SM(K) < XM(I) < SM(K-1).  THIS ASSUMES A FINER MASS SPACING
C FOR THE ISOCHRONE POINTS THAN THE INPUT *.ISO FILES!
         IF(XM(I).GT.SM(K)) K = K+1
         CALL PARROT(SQ,XQ,K,NTRK,XML(I),TEFFL(I),1)
         CALL PARROT(SQ,YQ,K,NTRK,XML(I),LOGL(I),1)
         BOL(I) = SOLBOL-(2.5D0*LOGL(I))
 200  CONTINUE
      N0 = 1
      NT = NMAX
      NPTS=NMAX
      N2 = NMAX -1
C COMPUTE THE NUMBER OF STARS (ANUM(J,M)) IN THE INTERVAL BETWEEN
C EEP J AND EEP #J+1
      CALL NUMSTAR(SP1,XM,ANUM,NT,N0,N2)
C ISOCHRONES ARE WRITTEN TO OUTPUT FILE
      CALL OUTISO(TEFFL,LOGL,XM,BOL,ANUM,Y,FEH,AFE,Z,ZEFF,OV,0.0d0,
     $               0,NPTS,2,NPTS,AX)
      RETURN
      END
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE shell(n,a)
C sorts an array a(1:n) into ascending numerical order by Shell's method
C N is input; A is replaced on output by its sorted rearrangement
      IMPLICIT NONE
      INTEGER n
      REAL*8 a(n)
      INTEGER i,j,inc
      REAL*8 v
      inc=1
1     inc=3*inc+1
      if(inc.le.n)goto 1
2     continue
        inc=inc/3
        do 11 i=inc+1,n
          v=a(i)
          j=i
3         if(a(j-inc).gt.v)then
            a(j)=a(j-inc)
            j=j-inc
            if(j.le.inc)goto 4
          goto 3
          endif
4         a(j)=v
11      continue
      if(inc.gt.1)goto 2
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software 0@.1Y..


C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE WORCOL(FE,TEFF,GL,WCOL)
      IMPLICIT REAL*8(A-H,O-Z)
c for a given value of [Fe/H], Teff and log G, return the following colors
c using the Worthey & Lee 2006 colour tables:
c           WCOL(8) where   WCOL(1) = BC
c                           WCOL(2) = U-B
c                           WCOL(3) = B-V
C                           WCOL(4) = V-R
C                           WCOL(5) = V-I
C                           WCOL(6) = V-J
C                           WCOL(7) = V-H
C                           WCOL(8) = V-K
      integer nind,nvk,nfe,ng,i,isave
      parameter    (nind=9,nvk=86,nfe=7,ng=13)
      real*8        aaa(nvk,nfe,ng,nind),wcol(8),clrs(9)
      real*8        fe,teff,gl,theta
      data isave/0/
      save aaa,isave
      
c
      if(isave.eq.0) then
         call readtable(aaa)
c         write(*,*)'readtable W'
         isave = 1
      endif
      theta = 5040.0d0/teff
      call teffinterp(aaa,gl,fe,theta,clrs,iflag)
c      write(*,*)(clrs(i),i=1,8)
      if(iflag.ne.0) then
c outside color table
         do i = 1, 8
            wcol(i) = -9.999
         enddo
      else
         wcol(1) = clrs(9)
         wcol(2) = clrs(1)
         wcol(3) = clrs(2)
         wcol(4) = clrs(3)
         wcol(5) = clrs(4)
         wcol(6) = clrs(7) - clrs(5)
         wcol(7) = clrs(7) - clrs(6)
         wcol(8) = clrs(7)
      endif
      return
      end

C--------------------------------------------------------------------------
      subroutine readtable(aaa)

      implicit none

      integer      nind,nvk,nfe,ng,i
      parameter    (nind=9,nvk=86,nfe=7,ng=13)
C     the 8th nind is the temperature
C     the 9th nind is BC_v
      real*8         aaa(nvk,nfe,ng,nind),g(ng),fe(nfe)
      integer       ife,ig,ivk,k
      real*8        x1,x2,xteff,x3


      data g / -0.5,0.0,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5 /
      data fe /  -2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5  /

      open(unit=63,status='old')
      do ife = 1,nfe
         do ig = 1,ng
            do ivk=1,75
               read(63,*) x1,x2,xteff,x3,
     z                 (aaa(ivk,ife,ig,k),k=1,7),aaa(ivk,ife,ig,9)
               if ( abs(x1 - fe(ife)).gt. 0.01 ) print*, 'bad fe'
               if ( abs(x2 - g(ig)).gt. 0.01 ) print*, 'bad g'
C              transform to THETA rather than Teff . . . 
               aaa(ivk,ife,ig,8) = 5040./xteff
            end do
         end do
      end do

      return
      end
C----------------------------------------------------------------------------

C     subroutine teffinterp
      subroutine teffinterp(aaa,grav,feh,theta,clrs,iflag)
C     given a THETA (=5040/Teff), return colors
      implicit none
      integer      nind,nvk,nfe,ng,iflag
      parameter    (nind=9,nvk=86,nfe=7,ng=13)
      real*8         grav,feh,theta,clrs(nind)
      real*8         aaa(nvk,nfe,ng,nind),g(ng),fe(nfe),vk(nvk)
C     the 8th nind is the temperature, 9th BC_v
C     local variables
      integer      jg,jfe,ivk,ind,jt,jt0,i
      real*8       c(nvk,nind),ffe,gg,dy,frac

      data g / -0.5,0.0,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5 /
      data fe /  -2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5  /
C     clrs are  1 U-B  2 B-V  3 V-R  4 V-I  5 J-K  6 H-K  7 V-K  (8 Teff)
C     a() is organised: a(nvk,nfe,ng,nind)                       (9 BCv)

C     find jfe and jg interpolation corners
      call locate(g,ng,grav,jg)
      if (jg.eq.0) jg = 1
      if (jg.eq.ng) jg = ng-1
      call locate(fe,nfe,feh,jfe)
      if (jfe.eq.0) jfe=1
      if (jfe.eq.nfe) jfe=nfe-1
C     fill c array with bilinear-interp results
      ffe = ( feh - fe(jfe) )/( fe(jfe+1)-fe(jfe)  )
      gg  = ( grav - g(jg) )/( g(jg+1) - g(jg)  )
      do ivk=1,75
         do ind=1,nind
            c(ivk,ind) = (1.0-ffe)*(1.0-gg)*aaa(ivk,jfe,jg,ind)
     z                 +  ffe     *(1.0-gg)*aaa(ivk,jfe+1,jg,ind)
     z                 +  ffe     *gg      *aaa(ivk,jfe+1,jg+1,ind)
     z                 + (1.0-ffe)*gg      *aaa(ivk,jfe,jg+1,ind)
         end do
      end do

C     find temperature (it's really THETA) and interpolate colors
      call locate(c(1,8),75,theta,jt)
      jt0 = jt
      if (jt.eq.0) jt = 1
      if (jt.gt.70) jt = 71
      do i=1,7
         call polint(c(jt,8),c(jt,i),5,theta,clrs(i),dy)
      end do
      call polint(c(jt,8),c(jt,9),5,theta,clrs(9),dy)
C     if requested temperature is out-of-bounds, use linear interpolation to
C     extrapolate. Return iflag = int(number of segments beyond the tabulated)
      if ( jt0 .eq. 0 ) then
         frac = (theta - c(1,8))/(c(2,8)-c(1,8))
         iflag = -1 + int(frac)
         do i=1,7
            clrs(i) = (1.0-frac)*c(1,i) + frac*c(2,i)
         end do
            clrs(9) = (1.0-frac)*c(1,9) + frac*c(2,9)
      end if
      if ( jt0 .eq. 75 ) then
         frac = (theta - c(74,8))/(c(75,8)-c(74,8))
         iflag = 1 + int(frac)
         do i=1,7
            clrs(i) = (1.0-frac)*c(74,i) + frac*c(75,i)
         end do
            clrs(9) = (1.0-frac)*c(74,9) + frac*c(75,9)
      end if

      return
      end
C-----------------------------------------------------------------------

C     subroutine vkinterp
      subroutine vkinterp(aaa,grav,feh,theta,clrs)
C     given a V-K (in clrs(7)), return colors and THETA = 5040/Teff
      implicit none
      integer      nind,nvk,nfe,ng
      parameter    (nind=9,nvk=86,nfe=7,ng=13)
      real*8         grav,feh,theta,clrs(nind)
      real*8         aaa(nvk,nfe,ng,nind),g(ng),fe(nfe),vk(nvk)
C     the 8th nind is the temperature, 9th BC_v
C     local variables
      integer      jg,jfe,ivk,ind,jt,i
      real*8       c(nvk,nind),ffe,gg,dy

      data g / -0.5,0.0,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5 /
      data fe /  -2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5  /
C     clrs are  1 U-B  2 B-V  3 V-R  4 V-I  5 J-K  6 H-K  7 V-K  (8 Teff)
C     a() is organised: a(nvk,nfe,ng,nind)                       (9 BCv)

C     find jfe and jg interpolation corners
      call locate(g,ng,grav,jg)
      if (jg.eq.0) jg = 1
      if (jg.eq.ng) jg = ng-1
      call locate(fe,nfe,feh,jfe)
      if (jfe.eq.0) jfe=1
      if (jfe.eq.nfe) jfe=nfe-1
C     fill c array with bilinear-interp results
      ffe = ( feh - fe(jfe) )/( fe(jfe+1)-fe(jfe)  )
      gg  = ( grav - g(jg) )/( g(jg+1) - g(jg)  )
      do ivk=1,75
         do ind=1,nind
            c(ivk,ind) = (1.0-ffe)*(1.0-gg)*aaa(ivk,jfe,jg,ind)
     z                 +  ffe     *(1.0-gg)*aaa(ivk,jfe+1,jg,ind)
     z                 +  ffe     *gg      *aaa(ivk,jfe+1,jg+1,ind)
     z                 + (1.0-ffe)*gg      *aaa(ivk,jfe,jg+1,ind)
         end do
      end do

C     find V-K and interpolate THETA and colors
      call locate(c(1,7),75,clrs(7),jt)
      if (jt.eq.0) jt = 1
      if (jt.gt.70) jt = 71
      do i=1,6
         call polint(c(jt,7),c(jt,i),5,clrs(7),clrs(i),dy)
      end do
      call polint(c(jt,7),c(jt,8),5,clrs(7),clrs(8),dy)
      call polint(c(jt,7),c(jt,9),5,clrs(7),clrs(9),dy)
      theta = clrs(8)

      return
      end
C-----------------------------------------------------------------------
C     subroutine viinterp
      subroutine viinterp(aaa,grav,feh,theta,clrs)
C     given a V-I (in clrs(4)), return colors and THETA = 5040/Teff
      implicit none
      integer      nind,nvk,nfe,ng
      parameter    (nind=9,nvk=86,nfe=7,ng=13)
      real*8         grav,feh,theta,clrs(nind)
      real*8         aaa(nvk,nfe,ng,nind),g(ng),fe(nfe),vk(nvk)
C     the 8th nind is the temperature, 9th BCv
C     local variables
      integer      jg,jfe,ivk,ind,jt,i
      real*8       c(nvk,nind),ffe,gg,dy

      data g / -0.5,0.0,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5 /
      data fe /  -2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5  /
C     clrs are  1 U-B  2 B-V  3 V-R  4 V-I  5 J-K  6 H-K  7 V-K  (8 Teff)
C     a() is organised: a(nvk,nfe,ng,nind)

C     find jfe and jg interpolation corners
      call locate(g,ng,grav,jg)
      if (jg.eq.0) jg = 1
      if (jg.eq.ng) jg = ng-1
      call locate(fe,nfe,feh,jfe)
      if (jfe.eq.0) jfe=1
      if (jfe.eq.nfe) jfe=nfe-1
C     fill c array with bilinear-interp results
      ffe = ( feh - fe(jfe) )/( fe(jfe+1)-fe(jfe)  )
      gg  = ( grav - g(jg) )/( g(jg+1) - g(jg)  )
      do ivk=1,75
         do ind=1,nind
            c(ivk,ind) = (1.0-ffe)*(1.0-gg)*aaa(ivk,jfe,jg,ind)
     z                 +  ffe     *(1.0-gg)*aaa(ivk,jfe+1,jg,ind)
     z                 +  ffe     *gg      *aaa(ivk,jfe+1,jg+1,ind)
     z                 + (1.0-ffe)*gg      *aaa(ivk,jfe,jg+1,ind)
         end do
      end do

C     find V-I and interpolate THETA and colors
      call locate(c(1,4),75,clrs(4),jt)
      if (jt.eq.0) jt = 1
      if (jt.gt.70) jt = 71
      do i=1,3
         call polint(c(jt,4),c(jt,i),5,clrs(4),clrs(i),dy)
      end do
      do i=5,9
         call polint(c(jt,4),c(jt,i),5,clrs(4),clrs(i),dy)
      end do
      theta = clrs(8)

      return
      end
C-----------------------------------------------------------------------
C------------------------------------------------------------------------
C     subroutine linear
C     quick linear interpolation
      subroutine linear(x,y,nxy,xin,yout,iOK)
      implicit none

      integer    nxy,iOK
      real*8       x(nxy),y(nxy),xin,yout
C     x and y are nxy long
C     xin is the nontabulated input x value
C     yout is the interpolated y guess.
C     iOK is -3 if X is not in ascending order, -2 if input x is
C     out-of-bounds by more than 1 xpoint spacing, -1 if out-of-bounds
C     by less than 1 xpoint spacing, 0 if all is OK

C     local variables
      integer    j,jl,ju,jm
      real*8       frac

      iOK = 0
      if ( x(2).lt.x(1) ) then
         print*, 
     z   'Error. Sub LINEAR. Input X array must be in ascending order.'
         yout = 0.0
         iOK = -3
         return
      end if

C     locate correct array element by bisection
      jl=0
      ju=nxy+1
10    if (ju-jl.gt.1) then
         jm = (ju+jl)/2
         if ((x(nxy).gt.x(1)).eqv.(xin.gt.x(jm))) then
            jl=jm
         else
            ju=jm
         end if
         goto 10
      end if
      j = jl
C     j is 0 or nxy if xin is off the grid

C     if off-grid, reset j and set output flag iOK
      if ( j.eq.0) then
         if ( xin.lt.(x(1)-(x(2)-x(1))) ) then
            iOK = -2
         else
            iOK = -1
         end if
         j=1
      endif
      if ( j.eq.nxy) then
         if ( xin.gt.(x(nxy)+(x(nxy)-x(nxy-1))) ) then
            iOK = -2
         else
            iOK = -1
         end if
         j = nxy-1
      end if

C     now interpolate/extrapolate
      frac = (xin - x(j))/(x(j+1)-x(j))
      yout = (1.0-frac)*y(j) + frac*y(j+1)


      return
      end
C     end subroutine linear ---------------------------------------------

C     -----NUMERICAL RECIPES routines: locate and polint
      SUBROUTINE LOCATE(XX,N,X,J)
      implicit none
      integer n,j
      real*8 xx(n),x
      integer jl,ju,jm


      JL=0
      JU=N+1
10    IF(JU-JL.GT.1)THEN
        JM=(JU+JL)/2
        IF((XX(N).GT.XX(1)).EQV.(X.GT.XX(JM)))THEN
          JL=JM
        ELSE
          JU=JM
        ENDIF
      GO TO 10
      ENDIF
      J=JL
      RETURN
      END

      SUBROUTINE POLINT(XA,YA,N,X,Y,DY)
      implicit none

      integer nmax,n
      PARAMETER (NMAX=10)
      real*8 XA(N),YA(N),C(NMAX),D(NMAX),x,y,dy

      integer ns,m,i
      real*8 dif,dift,ho,hp,w,den


      NS=1
      DIF=ABS(X-XA(1))
      DO 11 I=1,N 
        DIFT=ABS(X-XA(I))
        IF (DIFT.LT.DIF) THEN
          NS=I
          DIF=DIFT
        ENDIF
        C(I)=YA(I)
        D(I)=YA(I)
11    CONTINUE
      Y=YA(NS)
      NS=NS-1
      DO 13 M=1,N-1
        DO 12 I=1,N-M
          HO=XA(I)-X
          HP=XA(I+M)-X
          W=C(I+1)-D(I)
          DEN=HO-HP
          IF(DEN.EQ.0.) then
             write(*,*)'problem with DEN=',DEN
             stop
          ENDIF
          DEN=W/DEN
          D(I)=HP*DEN
          C(I)=HO*DEN
12      CONTINUE
        IF (2*NS.LT.N-M)THEN
          DY=C(NS+1)
        ELSE
          DY=D(NS)
          NS=NS-1
        ENDIF
        Y=Y+DY
13    CONTINUE
      RETURN
      END

C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE BASCOL(FEH,TEFF,LOGG,COLR)
C
C INTERPOLATES IN the Basel v3.1 (Padova 2000 verson) COLOR.TBL
C    INPUT: FEH, TEFF, LOG G  (UNCHANGED BY COLOR)
C    OUTPUT: COLR(8), WHERE COLR(1) = BC
C                           COLR(2) = U-B
C                           COLR(3) = B-V
C                           COLR(4) = V-R
C                           COLR(5) = V-I
C                           COLR(6) = V-J 
C                           COLR(7) = V-H
C                           COLR(8) = V-K 
c  there are lots more colors in the table that we could read in...
C
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 FEH,TEFF,LOGG,COLR(8)   ! use this line for all colours
      INTEGER ONE
      DIMENSION COL(8),TABL(15,70,20,8),GLOG(15,70,20),JMAX(70,20),
     1  IMAX(20),CZ(20),TEF(70,20),COLT(70),COLZ(20),colraw(15)
C TABL(J,I,K,M):   M==colors
C                  K==[Fe/H]
C                  I==Teff
C                  J=Log G
C ***BC 5/92 dimension statement for RYI colour table
c      DIMENSION COL(5),TABL(15,50,10,5),GLOG(15,50,10),JMAX(50,10),
c     1  IMAX(10),CZ(10),TEF(50,10),COLT(50),COLZ(10)
      COMMON/PARMIN/REFZ,SOLBOL,ALFA
      DATA NCOLOR/0/
      DATA ONE/64/
      SAVE
C
C SET UP AND READ IN COLOR TABLE ON FIRST CALL
C
      IF(NCOLOR.EQ.0) THEN
C
        OPEN(UNIT=ONE,STATUS='OLD') 
C
C N.B. REFZ IS THE REFERENCE (SOLAR) METALLICITY
	SOLARZ=DLOG10(REFZ)
	NCOLOR=8           !use this line for all colours
C
C READ IN COLOR TABLE
	K=0
	TLAST=0.
	FLAST=10.
	READ(ONE,*)
        READ(ONE,*)
1       READ(ONE,*,END=20) TE,G,FE,xxx,(colraw(M),M=1,15)
c        write(*,*)fe,te,g
c        write(*,*)fe,te,k
c convert the colours to V-J, V-H and V-K
        do m = 1, 5
           col(m) = colraw(m)
        enddo
        col(6) = colraw(6) - colraw(9) - colraw(10)
        col(7) = colraw(6) - colraw(10)
        col(8) = colraw(6)
c

	IF(TE.EQ.TLAST) GOTO 15
	IF(FE.EQ.FLAST) GOTO 10
	K=K+1
	KMAX=K
	CZ(K)=FE+SOLARZ
	FLAST=FE
	I=0
10      I=I+1
	IMAX(K)=I
	TEF(I,K)=DLOG10(TE)
	TLAST=TE
	J=0
15      J=J+1
	JMAX(I,K)=J
	GLOG(J,I,K)=G
c***BC put in test
        if(K.GT.20 .OR. I.GT.70 .OR. J.GT.15) THEN
           write(*,*) 'need to resize arrays in BASCOL'
           write(*,*)'i,j,k=',i,j,k
           stop
        endif
	DO M=1,NCOLOR
	  TABL(J,I,K,M)=COL(M)
	ENDDO
	GOTO 1
C
      ENDIF
C
C INTERPOLATE IN LOG G, LOG TEFF, AND [FE/H]
C
20    ZZ=FEH+SOLARZ
C ****BC 5/92
C if temp < 2000 K, or [Fe/H] outside limits 
C just put -9.999 for the colours and exit
      IF(TEFF.LT.2000.0 .or. feh.lt.-2.0001 .or. feh.gt.0.5001) THEN
         DO 666 III=1,NCOLOR
            COLR(III)=-9.999
 666     CONTINUE
         RETURN
      ENDIF
C ***BC 5/92 end of modification
      IF((CZ(KMAX)-ZZ)*(ZZ-CZ(1)).LT.0.) GO TO 980
      KK=1
30    KK=KK+1
      IF(CZ(KK).GT.ZZ) GOTO 30
      TE=DLOG10(TEFF)
      G=LOGG
      DO M=1,NCOLOR
	K0=MAX0(1,KK-2)
	KEND=MIN0(KMAX,KK+1)
	NK=KEND-K0+1
	DO K=K0,KEND
c          write(*,*)  10**tef(imax(k),k), 10**tef(1,k)
c          write(*,*) 10**te
	  IF((TEF(IMAX(K),K)-TE)*(TE-TEF(1,K)).LT.0.) GO TO 985
	  II=1
40        II=II+1
	  IF(TE.GT.TEF(II,K).AND.II.LT.IMAX(K)) GO TO 40
	  I0=MAX0(1,II-2)
	  IEND=MIN0(IMAX(K),II+1)
	  NI=IEND-I0+1
	  DO I=I0,IEND
c            write(*,*) jmax(i,k), glog(jmax(i,k),i,k), glog(1,i,k)
c            write(*,*) g
C ***BC 5/92 
C Kurucz does not have high G low Teff tables, so just put -9.999 for colour
C and continue
C CHECK THIS!!!!
c	    IF((GLOG(JMAX(I,K),I,K)-G)*(G-GLOG(1,I,K)).LT.0.) GO TO 990
	    IF((GLOG(JMAX(I,K),I,K)-G)*(G-GLOG(1,I,K)).LT.0.) THEN
               DO 555 III=1,NCOLOR
                  COLR(III)=-9.999
 555           CONTINUE
               RETURN
            ENDIF
C ***BC end of modification
	    JJ=1
45          JJ=JJ+1
	    IF(G.GT.GLOG(JJ,I,K)) GO TO 45
	    J0=MAX0(1,JJ-2)
	    JEND=MIN0(JMAX(I,K),JJ+1)
	    NJ=JEND-J0+1
	    CALL PARROT(GLOG(J0,I,K),TABL(J0,I,K,M),JJ-J0+1,NJ,G,
     1                                                 COLT(I),1)
	  ENDDO
	  CALL PARROT(TEF(I0,K),COLT(I0),II-I0+1,NI,TE,COLZ(K),1)
	ENDDO
	CALL PARROT(CZ(K0),COLZ(K0),KK-K0+1,NK,ZZ,COL(M),1)
	COLR(M)=COL(M)
      ENDDO
C
      RETURN
C
 980  write(*,*)'[Fe/H] =',FEH,
     1  ' OUTSIDE LIMITS OF COLOR CALIBRATION FILE'
      write(*,*)'feh problem'
      STOP
 985  write(*,*)TEFF,', [Fe/H] =',FEH, 
     1  ' OUTSIDE LIMITS OF COLOR CALIBRATION FILE'
      write(*,*)'teff problem'
      STOP
c990  TYPE *,LOGG,', ',TEFF,', [Fe/H] =',FEH,
c    1  ' OUTSIDE LIMITS OF COLOR CALIBRATION FILE'
c     STOP
      END


c----------------------------------------------------------------------------
c color-teff transforms: 
c NCT=1 Vandenberg & Clem - BVRI
c     4 PHOENIX           - UBVRI & 2MASS JHK
c     5    "              - HST ACS
c     6    "              - HST WFPC2
c     7    "              - SDSS ugriz
c     8 GCTreasury        - V & C/PHOENIX: V, I and ACS F606W, F814W
c Indices for these filters are:
c (Index) | 1      2      3      4      5      6      7      8
c ----------------------------------------------------------------
c         | U      B      V      R      I      J      H      K
c (ACS)   | F435W  F475W  F555W  F606W  F625W  F775W  F814W  F850LP
c (WFPC2) | F336W  F439W  F450W  F555W  F606W  F791W  F814W  F850LW  
c         | u      g      r      i      z
c GCTrsry | B      V      I      F606W  F814W  F606W  F814W (ACS/WFPC2)
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine color_init(nct,feh,afe)
      
      real*8 feh,feha,afe,t0,g0,bmv,vmr,vmi,bc,sby,sm1,sc1 
      integer nct

      feha=feh+5.0d-1*afe   !V&C, GW do not account for [a/Fe] variation
      if(nct.eq.1.or.nct.eq.8) CALL BVRI(0,FEHA,G0,T0,BMV,VMR,VMI,BC) !BVRI
      if(nct.eq.4) call phx_color_init(feh,afe,1) !UBVRIJHKs-PHX
      if(nct.eq.5) call phx_color_init(feh,afe,2) !HST-ACSWF 
      if(nct.eq.6) call phx_color_init(feh,afe,3) !HST-WFPC2
      if(nct.eq.7) call phx_color_init(feh,afe,4) !SDSS ugriz
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine get_mags(nct,feh,afe,mbol,grav,teff,filter)

      real*8 filter(13),color(13),feh,afe,grav,teff,feha,mv
      real*8 bmv,vmr,vmi,bc,mbol
      real*8 f606w_acs,f814w_acs,f606w_wfpc2,f814w_wfpc2
      real clrs(9)
      integer ncol,icol(9),nct,j

      data icol/4,4,8,8,8,8,5,7,8/

!adjust [Fe/H] for semi-empirical color-Teff transformations
      feha=feh+5.0d-1*afe

!***************************
!for Vandenberg & Clem BVRI
!***************************
      if(nct.eq.1.or.nct.eq.8)then
         call bvri(1,feha,grav,teff,bmv,vmr,vmi,bc)
         mv=mbol-bc
         filter(1)=bmv+mv       !B from B-V
         filter(2)=mv           !V
         filter(3)=mv-vmr       !R from V-R
         filter(4)=mv-vmi       !I from V-I
         if(nct.eq.8)then
!***************************
!     for the GC Treasury project
!***************************
!     V&C versions of HST mags: filters 1, 2, 3 = B, V, I
            filter(3)=filter(4) !don't keep R, change to I
! filters 4, 5, 6, 7 = ACS F606W + F814W, WFPC2 F606W + F814W
            filter(4)=f606w_acs(filter(2),filter(3)) !F606W (ACS)
            filter(5)=f814w_acs(filter(2),filter(3)) !F814W (ACS)
            filter(6)=f606w_wfpc2(filter(4),filter(5)) !F606W (WFPC2)
            filter(7)=f814w_wfpc2(filter(4),filter(5)) !F814W (WFPC2)
         endif
!***************************
!for PHOENIX synthetic
!***************************
      elseif(nct.gt.3.and.nct.lt.8)then !PHOENIX
         call phx_color_interp(teff,grav,color,nct)
         mv=mbol+color(1)
         if(nct.eq.4)then       !UBVRIJHK
            filter(1)=color(2)+color(3)+mv !U from U-B and B-V
            filter(2)=color(3)+mv !B from B-V
            filter(3)=mv        !V
            do j=4,icol(nct)
               filter(j)=mv-color(j) ![filter] from V-[filter]
            enddo
         elseif(nct.eq.5.or.nct.eq.6)then !HST ACS/WFPC2
            do j=1,icol(nct)
               filter(j)=mv-color(j+1) ![filter] from V-[filter]
            enddo
         elseif(nct.eq.7)then   !SDSS
            filter(1)=color(2)+mv !g: here "mv" is really "mg"
            do j=2,icol(nct)
               filter(j)=filter(j-1)-color(j)
            enddo
         endif
!anything else...
      else
         stop "INVALID COLOR-TEFF TRANSFORMATION SPECIFIED"
      endif

      return
      end
 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real*8 function f606w_acs(v,i)
c     Sirianni et al 2005 HST-ACS transformation equation from V and I to
c     F606W (Vega mag system)
      real*8 v,i
      if(v-i.ge.4.0d-1) then
         f606w_acs=v-2.6331d1+2.6398d1-3.4d-1*(v-i)+3.8d-2*(v-i)*(v-i)
      else
         f606w_acs=v-2.6394d1+2.6398d1-1.53d-1*(v-i)-9.6d-2*(v-i)*(v-i)
      endif
       
      return
      end
 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real*8 function f606w_wfpc2(v,i)
c     Sirianni et al 2005 HST-ACS transformation equation from ACS to WFPC2
c     Here V=F606W and I=F814W
      real*8 v,i
 
      f606w_wfpc2= v + 1.6d-2*(v-i) - 3.0d-2*(v-i)*(v-i) !- 4.168d0
      return
      end
  
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real*8 function f814w_acs(v,i)
c     Sirianni et al 2005 HST-ACS transformation equation from V and I to
c     F814W (Vega mag system)
      real*8 v,i
      if(v-i.ge.1.0d-1)then
         f814w_acs=i-2.5496d1+2.5501d1+1.4d-2*(v-i)-1.5d-2*(v-i)*(v-i)
      else
         f814w_acs=i-2.5489d1+2.5501d1-4.1d-2*(v-i)+9.3d-2*(v-i)*(v-i)
      endif
       
      return
      end
 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real*8 function f814w_wfpc2(v,i)
c     Sirianni et al 2005 HST-ACS transformation equation from ACS to WFPC2
c     Here V=F606W and I=F814W
      real*8 v,i
      f814w_wfpc2= i + 2.3d-2*(v-i) - 2.0d-3*(v-i)*(v-i) !- 4.601d0
      return
      end
 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                       END OF ALL_COLOR FILE                             c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                        C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C       PHOENIX COLOR INTERPOLATION IN TEFF, LOG G, [Fe/H],& [a/Fe]      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                        C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

ccccccccccccccccccccccc subroutine phx_color_init cccccccccccccccccccccccc
c this subroutine determines which filter set will be used, reads in     c
c the appropriate tables, and stores them in arrays in a common block    c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine phx_color_init(feh,afe,filter)

      implicit none
      integer filter, iz, iafe, nz, nzafe
      parameter(nz=10,nzafe=7)
c      parameter(nz=7,nzafe=7)
      real*8 feh, afe, z_0(nz), z_afe(nzafe)
      data z_0/-4.0,-3.5,-3.0,-2.5,-2.0,-1.5,-1.0,-0.5,0.0,0.5/
c      data   z_0/-2.5,-2.0,-1.5,-1.0,-0.5,0.0,0.5/
c      data z_afe/-2.5,-2.0,-1.5,-1.0,-0.5,0.0,0.5/
c filter is an integer representing the choice of filter set
c options include:
c     1. J-C UBVRI and 2MASS JHK
c     2. HST-ACS WFC
c     3. HST-WFPC2
c     4. SDSS ugriz

c set [a/Fe] index parameter
      iafe=0
      if(afe.eq.-0.2d0) iafe=1
      if(afe.eq. 0.0d0) iafe=2
      if(afe.eq. 0.2d0) iafe=3
      if(afe.eq. 0.4d0) iafe=4
      if(afe.eq. 0.6d0) iafe=5
      if(afe.eq. 0.8d0) iafe=6
      if(iafe.eq.0) stop "PROBLEM WITH [a/Fe] IN ISOCHRONE FILE"

c locate [Fe/H] in the array of available values
      call hunt(z_0,nz,feh,iz)

c check to make sure [Fe/H] contained within available values
c note that the alpha-enhanced colors begin at -2.5 (not -4)
c$$$      if(iz.eq.0) stop '[Fe/H] out of bounds!'
c$$$      if(iz.eq.nz) stop '[Fe/H] out of bounds!'
      if(iafe.eq.2.and.iz.lt.2) iz=2
      if(iafe.ne.2.and.iz.lt.5) iz=5
      if(iz.gt.nz-2) iz=nz-2
      if(iafe.eq.6.and.iz.gt.nz-3) iz=nz-3

c read in the required tables for a range of [Fe/H] at fixed [a/Fe]    
      call phx_color_read(iz,iafe,filter)

c interpolate in [Fe/H], only need once
      call z_interp(feh)

c all set for interpolation in T_eff and log(g)      
      return
      end


ccccccccccccccccccccccccc subroutine phx_color_read ccccccccccccccccccccccc
c reads in color tables based on specified Z, [a/Fe], filter set          c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine phx_color_read(iz,iafe,filter)

c put in a hack here, so that ACS tables (filter=2) are 
c saved to a separate common  block.... this means we only will have 
c to read in things once for V,I F606W, F814W colours

      implicit none
      integer iz,iafe,filter,nmax,ifile,nt,it,nz,nfltr
      parameter(nt=66,nz=10,nfltr=4,nmax=1000,ifile=21)
c      parameter(nt=66,nz=7,nfltr=4,nmax=1000,ifile=21)
      integer icol(nfltr),i,j,k
      integer ncol,itnum,isize,itemp
      integer ncolA,itnumA,isizeA,itempA
      character filename(4)*23,zfile(nz)*5,afile(6)*3,suffix(nfltr)*9
      real*8 zl,fz,coltbl,teff,ggl
      real*8 coltblA,teffA,gglA
      data zfile/'Zm4d0','Zm3d5','Zm3d0','Zm2d5','Zm2d0','Zm1d5',
     *     'Zm1d0','Zm0d5','Zp0d0','Zp0d5'/
      data afile/'am2','ap0','ap2','ap4','ap6','ap8'/
      data suffix/'STD_2MASS','HST_ACSWF','HST_WFPC2','SDSSugriz'/
      data icol/8,8,8,5/
      common/phx/coltbl(4,nmax,13),teff(4,nmax),ggl(4,nmax),
     *     itemp(4,nt),itnum(4),isize(4),ncol
      common/phxA/coltblA(4,nmax,13),teffA(4,nmax),gglA(4,nmax),
     *     itempA(4,nt),itnumA(4),isizeA(4),ncolA
      common/z/zl(4),fz(4)
c      data zl/-4.0,-3.5,-3.0,-2.5,-2.0,-1.5,-1.0,-0.5,0.0,0.5/

      ncol=icol(filter)
      ncolA=icol(filter)   !NOT USIN THIS RIGHT NOW.....
      do i=1,4
c set up filenames based on the input variables
         filename(i)(1:4)="phx/"
         filename(i)(5:9)=zfile(iz+i-2)
         filename(i)(10:12)=afile(iafe)
         filename(i)(13:13)="."
         do j=1,nfltr
            if(filter.eq.j) filename(i)(14:23)=suffix(j)
         enddo
c open table for reading
         write(*,*)filename(i)
         open(ifile,file=filename(i),status='old')
         read(ifile,*)
c read data in table
         do j=1,nmax
c***BCC 1/2007 my  hack -- ASSUMES Z is the same in each of the files
            if(filter.ne.2) then
               read(ifile,*,end=2) 
     *           teff(i,j),ggl(i,j),zl(i),(coltbl(i,j,k),k=1,ncol)
               if(j.ge.2.and.teff(i,j).eq.teff(i,j-1).and.
     *              ggl(i,j).eq.ggl(i,j-1)) then
c print warning message if duplicate lines exist in table
                  write(*,*) "DUPLICATE LINE IN COLOR TABLE"
                  write(*,*)zl(i),teff(i,j),teff(i,j-1),
     $                 ggl(i,j),ggl(i,j-1)
            endif
            elseif(filter.eq.2) then
               read(ifile,*,end=2) 
     *           teffA(i,j),gglA(i,j),zl(i),(coltblA(i,j,k),k=1,ncol)
               if(j.ge.2.and.teffA(i,j).eq.teffA(i,j-1).and.
     *              gglA(i,j).eq.gglA(i,j-1)) then
c print warning message if duplicate lines exist in table
                  write(*,*) "DUPLICATE LINE IN COLOR TABLE"
                  write(*,*)zl(i),teffa(i,j),teffa(i,j-1),
     $                 ggla(i,j),ggla(i,j-1)
               endif
            else
               write(*,*)'something wrong here'
               stop
            endif


         enddo
c close data file
 2       close(ifile)
c isize counts number of data points, itnum counts number of distinct T's,
c itemp gives location of each T value in the j-array
         if(filter.ne.2) then 
            isize(i)=j-1
            itemp(i,1)=1
            itnum(i)=1
            do j=2,isize(i)
               if(teff(i,j).gt.teff(i,j-1)) then
                  itnum(i)=itnum(i)+1
                  itemp(i,itnum(i))=j
               endif
            enddo
         else

C***BCC more of the hack

            isizeA(i)=j-1
            itempA(i,1)=1
            itnumA(i)=1
            do j=2,isizeA(i)
               if(teffA(i,j).gt.teffA(i,j-1)) then
                  itnumA(i)=itnumA(i)+1
                  itempA(i,itnumA(i))=j
               endif
            enddo
         endif



      enddo
      
      return
      end

ccccccccccccccccccccccc subroutine phx_color_interp ccccccccccccccccccccccc
c     interpolates the mags and colors based on T, g, [Fe/H], [a/Fe]      c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine phx_color_interp(tl,gl,color,nct)
      implicit none
      integer nmax,nt,i,j,k,iii,jj,inewt,inewg,tint,gint,nct
      integer ncol,itnum,isize,itemp
      integer ncolA,itnumA,isizeA,itempA
      parameter(nmax=1000,nt=66)
      real*8 zl,fz,fg(4),ft(4),coltbl,tl,gl,teff,ggl
      real*8 coltblA,teffA,gglA
      real*8 cln,t,qt(4),qg(4),colr(4,13),col(4,13,4),color(13)
      common/phx/coltbl(4,nmax,13),teff(4,nmax),ggl(4,nmax),itemp(4,nt),
     *     itnum(4),isize(4),ncol
      common/phxA/coltblA(4,nmax,13),teffA(4,nmax),gglA(4,nmax),
     *     itempA(4,nt),itnumA(4),isizeA(4),ncolA
      common/z/zl(4),fz(4)
      
      cln=dlog(1.0d1)



      if(nct.ne.5) then
      do iii=1,4
c     locate T in the Teff array
         t=dexp(cln*tl)
         do i=1,itnum(iii)-1
            if(t.ge.teff(iii,itemp(iii,i)).and.
     *           t.lt.teff(iii,itemp(iii,i+1))) inewt=i
         enddo
         if(inewt.lt.2) inewt=2
         if(inewt.gt.itnum(iii)-2) inewt=itnum(iii)-2

c     find interpolation coeff.'s in T
         do i=1,4
            qt(i)=teff(iii,itemp(iii,inewt+i-2))
         enddo
         call interp(qt,ft,t,4)

c     locate Log G in the Log G array for each T
         do j=1,4
            jj=inewt+j-2
            do k=itemp(iii,jj),itemp(iii,jj+1)-1
               if(gl.ge.ggl(iii,k).and.gl.lt.ggl(iii,k+1)) inewg=k
            enddo
            if(inewg.lt.2) inewg=2
            if(inewg.gt.itemp(iii,jj+1)-2) inewg=itemp(iii,jj+1)-2
            do k=1,4
               qg(k)=ggl(iii,inewg+k-2)
            enddo

c     find interpolation coefficients in Log G
            call interp(qg,fg,gl,4)

c     g-interpolation
            do k=1,ncol
               col(iii,k,j) = fg(1)*coltbl(iii,inewg-1,k)
     *                      + fg(2)*coltbl(iii,inewg  ,k)
     *                      + fg(3)*coltbl(iii,inewg+1,k)
     *                      + fg(4)*coltbl(iii,inewg+2,k)
            enddo               !k-loop
         enddo                  !j-loop

c     T-interpolation
         do i=1,ncol
            colr(iii,i)=ft(1)*col(iii,i,1)+ft(2)*col(iii,i,2)
     *           +ft(3)*col(iii,i,3)+ft(4)*col(iii,i,4)
         enddo                  !i-loop
        
      enddo                     !iii-loop

      else
c do the interpolation for the ACS filters, which are in different common block
      do iii=1,4
c     locate T in the Teff array
         t=dexp(cln*tl)
         do i=1,itnumA(iii)-1
            if(t.ge.teffA(iii,itempA(iii,i)).and.
     *           t.lt.teffA(iii,itempA(iii,i+1))) inewt=i
         enddo
         if(inewt.lt.2) inewt=2
         if(inewt.gt.itnumA(iii)-2) inewt=itnumA(iii)-2

c     find interpolation coeff.'s in T
         do i=1,4
            qt(i)=teffA(iii,itempA(iii,inewt+i-2))
         enddo
         call interp(qt,ft,t,4)

c     locate Log G in the Log G array for each T
         do j=1,4
            jj=inewt+j-2
            do k=itempA(iii,jj),itempA(iii,jj+1)-1
               if(gl.ge.gglA(iii,k).and.gl.lt.gglA(iii,k+1)) inewg=k
            enddo
            if(inewg.lt.2) inewg=2
            if(inewg.gt.itempA(iii,jj+1)-2) inewg=itempA(iii,jj+1)-2
            do k=1,4
               qg(k)=gglA(iii,inewg+k-2)
            enddo

c     find interpolation coefficients in Log G
            call interp(qg,fg,gl,4)

c     g-interpolation
            do k=1,ncol
               col(iii,k,j) = fg(1)*coltblA(iii,inewg-1,k)
     *                      + fg(2)*coltblA(iii,inewg  ,k)
     *                      + fg(3)*coltblA(iii,inewg+1,k)
     *                      + fg(4)*coltblA(iii,inewg+2,k)
            enddo               !k-loop
         enddo                  !j-loop

c     T-interpolation
         do i=1,ncol
            colr(iii,i)=ft(1)*col(iii,i,1)+ft(2)*col(iii,i,2)
     *           +ft(3)*col(iii,i,3)+ft(4)*col(iii,i,4)
         enddo                  !i-loop
        
      enddo                     !iii-loop
      endif



c     Z-interpolation
      do i=1,ncol
         color(i)=fz(1)*colr(1,i) + fz(2)*colr(2,i)
     *        + fz(3)*colr(3,i) + fz(4)*colr(4,i)
      enddo

      return
      end

ccccccccccccccccccccccccc subroutine z_interp cccccccccccccccccccccccccccc
c      Get interpolation coefficients in [Fe/H], only needed once.       c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine z_interp(feh)
      implicit none
      real*8 zl,fz,feh
      common/z/zl(4),fz(4)

      call interp(zl,fz,feh,4)

      return
      end

**********************************************************************
*                                HUNT                                *
**********************************************************************
      subroutine hunt(xx,n,x,jlo)
C Subroutine HUNT searches an array of length N to return a value
C of JLO such that X is located between array elements JLO and JLO+1.
C This search is based on the last value of JLO returned. JLO = 0 or
C JLO = N is returned to indicate that X is out of range of the array.
C This routine is taken from Numerical Recipes, pp.91-92
      implicit none
      integer n, jlo, jhi, inc, jm
      logical lascnd
      real*8 xx(n), x
      
      lascnd=xx(n).gt.xx(1)
      if(jlo.le.0.or.jlo.gt.n)then
C Input guess not useful; go immediately to bisection.
         jlo=0
         jhi=n+1
         goto 3
      endif
C Set up the hunting increment. 
      inc=1
C Hunt up.
      if(x.ge.xx(jlo).eqv.lascnd)then
 1       jhi=jlo+inc
         if(jhi.gt.n)then
C Done hunting, since off end of the table.
            jhi=n+1
         else if(x.ge.xx(jhi).eqv.lascnd)then
C Not done hunting.
            jlo=jhi
            inc=inc+inc
            goto 1
C Done hunting.
         endif
      else
C Hunt down.
         jhi=jlo
 2       jlo=jhi-inc
         if(jlo.lt.1) then
C Done hunting, since off end of table.
            jlo=0
         else if(x.lt.xx(jlo).eqv.lascnd)then
C Not done hunting.
            jhi=jlo
            inc=inc+inc
            goto 2
C Done hunting.
         endif
      endif
C Hunt is done; begin the final bisection phase.
 3    if(jhi-jlo.eq.1) return
      jm=(jhi+jlo)/2
      if(x.gt.xx(jm).eqv.lascnd)then
         jlo=jm
      else
         jhi=jm
      endif
      goto 3
      
      return
      end
      
C******************************************************************************
C               N-PT. Lagrangian Interpolation Coefficients
C******************************************************************************
      subroutine interp(a,b,x,n)
c {a} are the tabulated values for use in interpolation
c {b} are coefficients of the interpolating polynomial
c  x  is the abscissa to be interpolated
c  n  is the number of points to be used, interpolating polynomial
c     has order n-1 
      implicit none
      integer i,j,n
      real*8 a(n),b(n),x
      do i=1,n
         b(i)=1.0d0
         do j=1,n
            if(j.ne.i) b(i)=b(i)*(x-a(j))/(a(i)-a(j))
         enddo
      enddo
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                           END of PHX_COLOR.F                           c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
c BEGIN VRSUB SUBROUTINES: comprises all the useful SR's in the V-R code c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c BEGIN VANDENBERG & CLEM COLOR TRANSFORM SUBROUTINES                    c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE BVRI(MODE,FE,GV,TEFF,BMV,VMR,VMI,BCV)
C ----------------------------------------------------------------------
C *** THIS READS, AND INTERPOLATES IN, bvrilo.data AND bvrihi.data (SEE
C     VANDENBERG & CLEM 2003, AJ, 126, IN PRESS).            
C *** SET MODE=0 TO READ THE TABLES AND TO INTERPOLATE THEM TO THE
C                DESIRED VALUE OF [Fe/H]  (= FE IN THE ARGUMENT LIST)   
C *** SET MODE=1 TO INTERPOLATE FOR THE COLORS AT THE DESIRED VALUES OF
C                [Fe/H], log g, and log Teff (GV AND TEFF REPRESENT THE
C                LAST TWO OF THESE QUANTITIES IN THE ARGUMENT LIST).
C                (NOTE THAT THE [Fe/H] INTERPOLATION IS CARRIED OUT ONLY
C                WHEN THE SUBROUTINE IS CALLED WITH MODE=0 OR MODE=-1.)
C *** SET MODE=-1 IF THE INPUT TABLES, WHICH HAVE ALREADY BEEN READ
C                USING MODE=0, ARE TO BE RE-INTERPOLATED TO BE
C                CONSISTENT WITH A NEW VALUE OF [Fe/H].        
C *** THE OUTPUT CONSISTS OF THE B-V, V-R, AND V-I COLORS, ALONG WITH
C     THE BOLOMETRIC CORRECTION TO V (ON THE SCALE WHERE THE SUN HAS
C     M_bol = 4.75 and M_V = 4.84).  THESE QUANTITIES ARE REPRESENTED,
C     IN TURN, BY BMV, VMR, VMI, AND BCV IN THE ARGUMENT LIST. (B AND V
C     ARE ON THE JOHNSON SYSTEM, R AND I ON THE COUSINS SYSTEM.)
C *** MODE MUST BE SET TO ZERO THE FIRST TIME THAT THIS SUBROUTINE IS
C     CALLED, AS IT IS ONLY WHEN MODE=0 THAT BVRILO.DATA AND BVRIHI.DATA
C     ARE READ.  ONCE THE TABLES HAVE BEEN INTERPOLATED TO A PARTICULAR
C     VALUE OF [Fe/H] (USING MODE=0 OR MODE=-1), INTERPOLATIONS FOR THE
C     COLORS APPROPRIATE TO INPUT VALUES OF [Fe/H], log g, and log Teff
C     MAY BE CARRIED OUT *ANY NUMBER OF TIMES* USING MODE=1 (THE ONLY
C     VALUE OF MODE FOR WHICH COLOR INFORMATION IS OBTAINED).
C *** NOTE THAT THE BVRILO.DATA AND BVRIHI.DATA FILES ARE "ATTACHED" TO
C     UNITS 10 AND 11, RESPECTIVELY.  WHENEVER AN [Fe/H] INTERPOLATION 
C     IS CARRIED OUT, A MESSAGE IS SENT TO UNIT 6 ("ATTACHED" TO THE  
C     TERMINAL) TO INDICATE THE VALUE OF [FE/H] THAT APPLIES TO
C     SUBSEQUENT INTERPOLATIONS FOR log g and log Teff. 
C *** CHKBVRI.FOR (ALSO ON DISK) MAY BE USED TO TEST THE INTERPOLATION
C     SUBROUTINE
C *** NOTE, IN THE OPEN STATEMENTS BELOW, THAT THE NAMES OF THE INPUT
C *** DATA FILES ARE WRITTEN IN LOWER-CASE TEXT
C ----------------------------------------------------------------------
*23456789012345678901234567890123456789012345678901234567890123456789012
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(34,7,6,4),B(34,7,4),BV(11,12,6),VR(11,12,6), 
     1   VI(11,12,6),BC(11,12,6),C(11,12),D(11,12),E(11,12),F(11,12), 
     2   H(4,4),TB(34),TT(34),GG(7),TA(11),T(11),G(12),FEH(6),P(4),Q(4),    
     3   R(4),S(4),AG(4),AT(4),X(11),Y(11)
      CHARACTER*19 NAMCOL
      CHARACTER*10 NAMFE                                       
      INTEGER ILO, IHI
      PARAMETER(ILO=17,IHI=18)
      SAVE
 1000 FORMAT(I3,13X,I3,13X,I3,13X,I3)                                   
 1001 FORMAT(13F6.0)                                                    
 1002 FORMAT(20F4.1)                                                    
 1003 FORMAT(13F6.3)                                                    
 1004 FORMAT('     Color grid interpolated to [FE/H] =',F6.2)           
 1005 FORMAT(' log g =',F6.3,'  log Teff =',F7.4,' ARE OUTSIDE THE',          
     1   1X,'RANGE OF THE COLOR TABLES')                          
 1006 FORMAT(A10,F5.2,A19)
 1007 FORMAT(1X,' **** INPUT DATA FILE DOES NOT EXIST **** ')      
      IF(MODE) 11,3,16                                                  
C *** WHEN MODE=0 THE COLOR TRANSFORMATION TABLES ARE READ
    3 TBND=LOG10(5000.)                                                
      GBND=2.
      OPEN(UNIT=ILO,FILE='vdb/bvrilo.data',ERR=49,STATUS='OLD')
      OPEN(UNIT=IHI,FILE='vdb/bvrihi.data',ERR=49,STATUS='OLD')      
      READ(ILO,1000) NT,NG,NFE,NDX                                       
      READ(ILO,1001) (TA(I),I=1,NT)                                      
      READ(ILO,1002) (G(I),I=1,NG)                                       
      DO 4 I=1,NT                                                       
    4 T(I)=LOG10(TA(I))                                                
      DO 5 K=1,NFE                                                      
      READ(ILO,1006) NAMFE,FEH(K),NAMCOL                                    
      DO 5 J=1,NG                                                       
      READ(ILO,1003) (BV(I,J,K),I=1,NT)                                  
    5 CONTINUE                                                          
      DO 6 K=1,NFE                                                      
      READ(ILO,1006) NAMFE,FEH(K),NAMCOL                                 
      DO 6 J=1,NG                                                       
      READ(ILO,1003) (VR(I,J,K),I=1,NT)                                  
    6 CONTINUE                                                          
      DO 7 K=1,NFE                                                      
      READ(ILO,1006) NAMFE,FEH(K),NAMCOL                                   
      DO 7 J=1,NG                                                       
      READ(ILO,1003) (VI(I,J,K),I=1,NT)                                  
    7 CONTINUE                                                          
      DO 8 K=1,NFE                                                      
      READ(ILO,1006) NAMFE,FEH(K),NAMCOL                                   
      DO 8 J=1,NG
      READ(ILO,1003) (BC(I,J,K),I=1,NT)
    8 CONTINUE                                             
      READ(IHI,1000) NTT,NGG,NFE,NDXX                                    
      READ(IHI,1001) (TB(I),I=1,NTT)                                     
      READ(IHI,1002) (GG(I),I=1,NGG)                                     
      DO 9 I=1,NTT                                                      
    9 TT(I)=LOG10(TB(I))                                               
      DO 10 L=1,NDXX                                                    
      DO 10 K=1,NFE                                                     
      READ(IHI,1006) NAMFE,FEH(K),NAMCOL
      DO 10 J=1,NGG                                                     
      READ(IHI,1003) (A(I,J,K,L),I=1,NTT)
   10 CONTINUE
      CLOSE(UNIT=ILO,STATUS='KEEP')
      CLOSE(UNIT=IHI,STATUS='KEEP')         
C *** WHEN MODE=0 OR MODE=-1, COLOR TRANSFORMATION TABLES ARE
C *** CREATED FOR THE INPUT [Fe/H] VALUE USING LINEAR INTERPOLATION
   11 DO 12 M=2,NFE                                                     
      K=M-1                                                             
      IF(FE.LE.FEH(M)) GO TO 13                                         
   12 CONTINUE                                                          
   13 M=K+1                                                             
      SLOPE=(FE-FEH(K))/(FEH(M)-FEH(K))                                 
      DO 14 J=1,NG                                                      
      DO 14 I=1,NT                                                      
      C(I,J)=BV(I,J,K)+SLOPE*(BV(I,J,M)-BV(I,J,K))                      
      D(I,J)=VR(I,J,K)+SLOPE*(VR(I,J,M)-VR(I,J,K))                      
      E(I,J)=VI(I,J,K)+SLOPE*(VI(I,J,M)-VI(I,J,K))                      
      F(I,J)=BC(I,J,K)+SLOPE*(BC(I,J,M)-BC(I,J,K))                      
   14 CONTINUE                                                          
      DO 15 L=1,NDXX                                                    
      DO 15 J=1,NGG                                                     
      DO 15 I=1,NTT                                                     
      B(I,J,L)=A(I,J,K,L)+SLOPE*(A(I,J,M,L)-A(I,J,K,L))                 
   15 CONTINUE                                                          
      !WRITE(0,1004) FE                                                  
C *** WHEN MODE=0 OR MODE=-1, CONTROL RETURNS TO THE CALLING
C *** PROGRAM ONCE THE [Fe/H] INTERPOLATION IS CARRIED OUT (I.E., NO
C *** INTERPOLATIONS ARE PERFORMED FOR INPUT VALUES OF GV AND TEFF)
      GO TO 50
C *** WHEN MODE=1, INTERPOLATIONS ARE CARRIED OUT FOR THE INPUT VALUES
C *** OF [Fe/H], log g, AND log Teff (= FE, GV, AND TEFF IN THE ARGUMENT
C *** LIST)                                    
   16 IF(TEFF.LE.TBND) GO TO 18                                         
      IF(GV.GE.GBND) GO TO 40                                           
      IF(TEFF.LE.T(NT)) GO TO 18                                        
C      WRITE(0,1005) GV,TEFF                                             
C *** EXECUTION HALTS WITH A "STOP30" CODE IF THE INPUT TEMPERATURE IS
C *** OUTSIDE THE RANGE OF THE TABLES ON THE HIGH SIDE
C      STOP30                                                            
      BCV=-99.999d0
      BMV=-99.999d0
      VMR=-99.999d0
      VMI=-99.999d0
      RETURN
C *** THE NEXT SECTION ASSUMES THAT THE LOW-TEMPERATURE TABLES ARE THE
C *** RELEVANT ONES TO USE
   18 NGM=NG-1                                                          
      DO 19 I=3,NGM                                                     
      MG=I-2                                                            
      IF(GV.LE.G(I)) GO TO 20                                           
   19 CONTINUE                                                          
      IF(GV.LE.G(NG)) GO TO 20
C *** SOME EXTRAPOLATION OUTSIDE THE RANGE OF log g CONSIDERED IN THE
C *** TABLES IS PERMITTED.  TO AVOID ANY EXTRAPOLATIONS, THE NEXT TWO
C *** LINES SHOULD BE ACTIVATED                                          
C     WRITE(6,1005) GV,TEFF                                             
C     STOP31                                                            
   20 R(1)=G(MG)                                                        
      R(2)=G(MG+1)                                                      
      R(3)=G(MG+2)                                                      
      R(4)=G(MG+3)                                                      
      CALL LGRAN4(R,AG,GV)                                              
      NTM=NT-1                                                          
      DO 25 I=3,NTM                                                     
      MT=I-2                                                            
      IF(TEFF.LE.T(I)) GO TO 30                                         
   25 CONTINUE                                                          
      IF(TEFF.LE.T(NT)) GO TO 30                                        
C *** SOME EXTRAPOLATION OUTSIDE THE RANGE OF log Teff CONSIDERED IN
C *** THE TABLES IS PERMITTED.  TO AVOID ANY EXTRAPOLATIONS, THE NEXT
C *** TWO LINES SHOULD BE ACTIVATED.
C     WRITE(6,1005) GV,TEFF                                             
C     STOP32                                                            
   30 R(1)=T(MT)                                                        
      R(2)=T(MT+1)                                                      
      R(3)=T(MT+2)                                                      
      R(4)=T(MT+3)                                                      
      CALL LGRAN4(R,AT,TEFF)                                            
      L=MG-1                                                            
      DO 35 I=1,4                                                       
      L=L+1                                                             
      P(I)=AT(1)*C(MT,L)+AT(2)*C(MT+1,L)+AT(3)*C(MT+2,L)+AT(4)*C(MT+3,L)
      Q(I)=AT(1)*D(MT,L)+AT(2)*D(MT+1,L)+AT(3)*D(MT+2,L)+AT(4)*D(MT+3,L)
      R(I)=AT(1)*E(MT,L)+AT(2)*E(MT+1,L)+AT(3)*E(MT+2,L)+AT(4)*E(MT+3,L)
   35 S(I)=AT(1)*F(MT,L)+AT(2)*F(MT+1,L)+AT(3)*F(MT+2,L)+AT(4)*F(MT+3,L)
C *** 4-POINT LAGRANGIAN INTERPOLATION IS USED TO DERIVE THE VALUES OF
C *** B-V, V-R, V-I, AND BC_V CORRESPONDING TO THE INPUT VALUES OF
C *** [Fe/H], log g, and log Teff.
      BMV=AG(1)*P(1)+AG(2)*P(2)+AG(3)*P(3)+AG(4)*P(4)                   
      VMR=AG(1)*Q(1)+AG(2)*Q(2)+AG(3)*Q(3)+AG(4)*Q(4)                   
      VMI=AG(1)*R(1)+AG(2)*R(2)+AG(3)*R(3)+AG(4)*R(4)                   
      BCV=AG(1)*S(1)+AG(2)*S(2)+AG(3)*S(3)+AG(4)*S(4)                  
C *** SPLINE INTERPOLATION IS USED TO FIND THE VALUE OF BC_V AT LOW
C *** TEMPERATURES.
      DO 38 I=1,NT                                                      
      X(I)=T(I)                                                         
   38 Y(I)=AG(1)*F(I,MG)+AG(2)*F(I,MG+1)+AG(3)*F(I,MG+2)+AG(4)*F(I,MG+3)
      CALL INTEP(TEFF,BCV,X,Y,NT,IER)                                  
      GO TO 50                                                          
C *** THE NEXT SECTION ASSUMES THAT THE HIGH-TEMPERATURE TABLES ARE THE
C *** RELEVANT ONES TO USE.
   40 IF(TEFF.LE.TT(NTT)) GO TO 42                                      
      WRITE(0,1005) GV,TEFF                                             
      STOP33                                                            
   42 NGM=NGG-1                                                         
      DO 43 I=3,NGM                                                     
      MG=I-2                                                            
      IF(GV.LE.GG(I)) GO TO 44                                          
   43 CONTINUE                                                          
      IF(GV.LE.GG(NGG)) GO TO 44                                         
C *** SOME EXTRAPOLATION OUTSIDE THE RANGE OF log g CONSIDERED IN THE
C *** TABLES IS PERMITTED.  TO AVOID ANY EXTRAPOLATIONS, THE NEXT TWO
C *** LINES SHOULD BE ACTIVATED.
C     WRITE(6,1005) GV,TEFF                                             
C     STOP34                                                            
   44 R(1)=GG(MG)                                                       
      R(2)=GG(MG+1)                                                     
      R(3)=GG(MG+2)                                                     
      R(4)=GG(MG+3)                                                     
      CALL LGRAN4(R,AG,GV)                                              
      NTM=NTT-1                                                         
      DO 45 I=3,NTM                                                     
      MT=I-2                                                            
      IF(TEFF.LE.TT(I)) GO TO 46                                        
   45 CONTINUE                                                          
      IF(TEFF.LE.TT(NTT)) GO TO 46                                      
C *** SOME EXTRAPOLATION OUTSIDE THE RANGE OF log Teff CONSIDERED IN
C *** THE TABLES IS PERMITTED.  TO AVOID ANY EXTRAPOLATIONS, THE NEXT
C *** TWO LINES SHOULD SHOULD BE ACTIVATED.
C     WRITE(6,1005) GV,TEFF                                             
C     STOP35                                                            
   46 R(1)=TT(MT)                                                       
      R(2)=TT(MT+1)                                                     
      R(3)=TT(MT+2)                                                     
      R(4)=TT(MT+3)                                                     
      CALL LGRAN4(R,AT,TEFF)                                            
      L=MG-1                                                            
      DO 48 I=1,4                                                       
      L=L+1                                                             
      DO 48 K=1,NDXX                                                    
      H(I,K)=AT(1)*B(MT,L,K)+AT(2)*B(MT+1,L,K)+AT(3)*B(MT+2,L,K)+       
     1   AT(4)*B(MT+3,L,K)                                              
   48 CONTINUE                                                          
C *** 4-POINT LAGRANGIAN INTERPOLATION IS USED TO DERIVE THE VALUES OF
C *** B-V, V-R, V-I, AND BC_V CORRESPONDING TO THE INPUT VALUES OF
C *** [Fe/H], log g, and log Teff.
      BMV=AG(1)*H(1,1)+AG(2)*H(2,1)+AG(3)*H(3,1)+AG(4)*H(4,1)           
      VMR=AG(1)*H(1,2)+AG(2)*H(2,2)+AG(3)*H(3,2)+AG(4)*H(4,2)           
      VMI=AG(1)*H(1,3)+AG(2)*H(2,3)+AG(3)*H(3,3)+AG(4)*H(4,3)           
      BCV=AG(1)*H(1,4)+AG(2)*H(2,4)+AG(3)*H(3,4)+AG(4)*H(4,4)
      GO TO 50
   49 WRITE(0,1007)
      STOP          
   50 RETURN                                                            
      END

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
             
      SUBROUTINE UVBY(MODE,FE,GV,TEFF,SBY,SM1,SC1,BCV)           
C ----------------------------------------------------------------------
C *** THIS READS, AND INTERPOLATES IN, UVBYLO.DATA AND UVBYHI.DATA (SEE
C     CLEM, VANDENBERG, GRUNDAHL, & BELL, AJ, IN PRESS).            
C *** SET MODE=0 TO READ THE TABLES AND TO INTERPOLATE THEM TO THE
C                DESIRED VALUE OF [Fe/H]  (= FE IN THE ARGUMENT LIST)   
C *** SET MODE=1 TO INTERPOLATE FOR THE COLORS AT THE DESIRED VALUES OF
C                [Fe/H], log g, and log Teff (GV AND TEFF REPRESENT THE
C                LAST TWO OF THESE QUANTITIES IN THE ARGUMENT LIST).
C                (NOTE THAT THE [Fe/H] INTERPOLATION IS CARRIED OUT ONLY
C                WHEN THE SUBROUTINE IS CALLED WITH MODE=0 OR MODE=-1.)
C *** SET MODE=-1 IF THE INPUT TABLES, WHICH HAVE ALREADY BEEN READ
C                USING MODE=0, ARE TO BE RE-INTERPOLATED TO BE
C                CONSISTENT WITH A NEW VALUE OF [Fe/H].        
C *** THE OUTPUT CONSISTS OF THE b-y, m1, AND c1 INDICES, ALONG 
C     WITH THE BOLOMETRIC CORRECTION TO V (ON THE SCALE WHERE THE SUN 
C     HAS M_bol = 4.75 and M_V = 4.84).  THESE QUANTITIES ARE 
C     REPRESENTED, IN TURN, BY SBY, SM1, SC1, AND BCV IN THE 
C     ARGUMENT LIST.   NOTE THAT IT IS ASSUMED THE BOLOMETRIC 
C     CORRECTIONS TO STROMGREN y ARE IDENTICAL TO THOSE IN JOHNSON V.
C *** MODE MUST BE SET TO ZERO THE FIRST TIME THAT THIS SUBROUTINE IS
C     CALLED, AS IT IS ONLY WHEN MODE=0 THAT UVBYLO.DATA AND UVBYHI.DATA
C     ARE READ.  ONCE THE TABLES HAVE BEEN INTERPOLATED TO A PARTICULAR
C     VALUE OF [Fe/H] (USING MODE=0 OR MODE=-1), INTERPOLATIONS FOR THE
C     COLORS APPROPRIATE TO INPUT VALUES OF [Fe/H], log g, and log Teff
C     MAY BE CARRIED OUT *ANY NUMBER OF TIMES* USING MODE=1 (THE ONLY
C     VALUE OF MODE FOR WHICH COLOR INFORMATION IS OBTAINED).
C *** NOTE THAT THE UVBYLO.DATA AND UVBYHI.DATA FILES ARE "ATTACHED" TO
C     UNITS 12 AND 13, RESPECTIVELY.  WHENEVER AN [Fe/H] INTERPOLATION 
C     IS CARRIED OUT, A MESSAGE IS SENT TO UNIT 6 ("ATTACHED" TO THE  
C     TERMINAL) TO INDICATE THE VALUE OF [FE/H] THAT APPLIES TO
C     SUBSEQUENT INTERPOLATIONS FOR log g and log Teff. 
C *** CHKUVBY.FOR (ALSO ON DISK) MAY BE USED TO TEST THE INTERPOLATION
C     SUBROUTINE
C ----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8 M1(13,12,8)
      DIMENSION A(51,7,8,4),B(51,7,4),BY(13,12,8),C1(13,12,8),
     1  BC(13,12,8),C(13,12),D(13,12),E(13,12),F(13,12),H(4,4),TB(51),
     2  TT(51),GG(7),TA(13),T(13),G(12),FEH(8),O(4),P(4),Q(4),R(4),
     3  AG(4),AT(4),X(13),Y(13)
      CHARACTER*19 NAMCOL
      CHARACTER*10 NAMFE
      INTEGER ILO,IHI
      PARAMETER(ILO=19,IHI=20)
      SAVE
 1000 FORMAT(I3,13X,I3,13X,I3,13X,I3)                                   
 1001 FORMAT(13F6.0)                                                    
 1002 FORMAT(20F4.1)                                                    
 1003 FORMAT(13F6.3)                                                    
 1004 FORMAT('     Color grid interpolated to [FE/H] =',F6.2)           
 1005 FORMAT(1X,7HLOG G =,F6.3,11H LOG TEFF =,F7.4,8H OUTSIDE,          
     1   1X,11HCOLOR TABLE)                                             
 1006 FORMAT(A10,F5.2,A19)
 1007 FORMAT(1X,' **** INPUT DATA FILE DOES NOT EXIST **** ')      
      IF(MODE) 11,3,16                                                  
C *** WHEN MODE=0 THE COLOR TRANSFORMATION TABLES ARE READ
    3 TBND=LOG10(5500.)                                                
      GBND=2.
      OPEN(UNIT=ILO,ERR=49,FILE='uvbylo.data',STATUS='OLD')
      OPEN(UNIT=IHI,ERR=49,FILE='uvbyhi.data',STATUS='OLD')      
      READ(ILO,1000) NT,NG,NFE,NDX                                       
      READ(ILO,1001) (TA(I),I=1,NT)                                      
      READ(ILO,1002) (G(I),I=1,NG)                                       
      DO 4 I=1,NT                                                       
    4 T(I)=LOG10(TA(I))                                                
      DO 5 K=1,NFE                                                      
      READ(ILO,1006) NAMFE,FEH(K),NAMCOL                                    
      DO 5 J=1,NG                                                       
      READ(ILO,1003) (BY(I,J,K),I=1,NT)                                  
    5 CONTINUE                                                          
      DO 6 K=1,NFE                                                      
      READ(ILO,1006) NAMFE,FEH(K),NAMCOL                                 
      DO 6 J=1,NG                                                       
      READ(ILO,1003) (M1(I,J,K),I=1,NT)                                  
    6 CONTINUE                                                          
      DO 7 K=1,NFE                                                      
      READ(ILO,1006) NAMFE,FEH(K),NAMCOL                                   
      DO 7 J=1,NG                                                       
      READ(ILO,1003) (C1(I,J,K),I=1,NT)                                  
    7 CONTINUE                                                          
      DO 8 K=1,NFE                                                      
      READ(ILO,1006) NAMFE,FEH(K),NAMCOL                                   
      DO 8 J=1,NG
      READ(ILO,1003) (BC(I,J,K),I=1,NT)
    8 CONTINUE                                             
      READ(IHI,1000) NTT,NGG,NFE,NDXX                                    
      READ(IHI,1001) (TB(I),I=1,NTT)                                     
      READ(IHI,1002) (GG(I),I=1,NGG)                                     
      DO 9 I=1,NTT                                                      
    9 TT(I)=LOG10(TB(I))                                               
      DO 10 L=1,NDXX                                                    
      DO 10 K=1,NFE                                                     
      READ(IHI,1006) NAMFE,FEH(K),NAMCOL
      DO 10 J=1,NGG                                                     
      READ(IHI,1003) (A(I,J,K,L),I=1,NTT)
   10 CONTINUE
      CLOSE(UNIT=ILO,STATUS='KEEP')
      CLOSE(UNIT=IHI,STATUS='KEEP')         
C *** WHEN MODE=0 OR MODE=-1, COLOR TRANSFORMATION TABLES ARE
C *** CREATED FOR THE INPUT [Fe/H] VALUE USING LINEAR INTERPOLATION
   11 DO 12 M=2,NFE                                                     
      K=M-1                                                             
      IF(FE.LE.FEH(M)) GO TO 13                                         
   12 CONTINUE                                                          
   13 M=K+1                                                             
      SLOPE=(FE-FEH(K))/(FEH(M)-FEH(K))                                 
      DO 14 J=1,NG                                                      
      DO 14 I=1,NT                                                      
      C(I,J)=BY(I,J,K)+SLOPE*(BY(I,J,M)-BY(I,J,K))
      D(I,J)=M1(I,J,K)+SLOPE*(M1(I,J,M)-M1(I,J,K))
      E(I,J)=C1(I,J,K)+SLOPE*(C1(I,J,M)-C1(I,J,K))
      F(I,J)=BC(I,J,K)+SLOPE*(BC(I,J,M)-BC(I,J,K))
   14 CONTINUE                                                          
      DO 15 L=1,NDXX                                                    
      DO 15 J=1,NGG                                                     
      DO 15 I=1,NTT                                                     
      B(I,J,L)=A(I,J,K,L)+SLOPE*(A(I,J,M,L)-A(I,J,K,L))                 
   15 CONTINUE                                                          
      !WRITE(0,1004) FE  
C *** WHEN MODE=0 OR MODE=-1, CONTROL RETURNS TO THE CALLING
C *** PROGRAM ONCE THE [Fe/H] INTERPOLATION IS CARRIED OUT (I.E., NO
C *** INTERPOLATIONS ARE PERFORMED FOR INPUT VALUES OF GV AND TEFF)
      GO TO 50
C *** WHEN MODE=1, INTERPOLATIONS ARE CARRIED OUT FOR THE INPUT VALUES
C *** OF [Fe/H], log g, AND log Teff (= FE, GV, AND TEFF IN THE ARGUMENT
C *** LIST)
   16 IF(TEFF.LE.TBND) GO TO 18                                         
      IF(GV.GE.GBND) GO TO 40                                           
      IF(TEFF.LE.T(NT)) GO TO 18                                        
C *** EXECUTION HALTS WITH A "STOP30" CODE IF THE INPUT TEMPERATURE IS
C *** OUTSIDE THE RANGE OF THE TABLES ON THE HIGH SIDE
C      WRITE(0,1005) GV,TEFF                                             
C      STOP30                                                            
      BCV=-99.999d0
      BMV=-99.999d0
      VMR=-99.999d0
      VMI=-99.999d0
      RETURN
C *** THE NEXT SECTION ASSUMES THAT THE LOW-TEMPERATURE TABLES ARE THE
C *** RELEVANT ONES TO USE
   18 NGM=NG-1                                                          
      DO 19 I=3,NGM                                                     
      MG=I-2                                                            
      IF(GV.LE.G(I)) GO TO 20                                           
   19 CONTINUE                                                          
      IF(GV.LE.G(NG)) GO TO 20                                          
C *** SOME EXTRAPOLATION OUTSIDE THE RANGE OF log g CONSIDERED IN THE 
C *** TABLES IS PERMITTED.  TO AVOID ANY EXTRAPOLATIONS, THE NEXT TWO
C *** LINES SHOULD BE ACTIVATED   
C     WRITE(6,1005) GV,TEFF
C     STOP31
   20 R(1)=G(MG)                                                        
      R(2)=G(MG+1)                                                      
      R(3)=G(MG+2)                                                      
      R(4)=G(MG+3)                                                      
      CALL LGRAN4(R,AG,GV)                                              
      NTM=NT-1                                                          
      DO 25 I=3,NTM                                                     
      MT=I-2                                                            
      IF(TEFF.LE.T(I)) GO TO 30                                         
   25 CONTINUE                                                          
      IF(TEFF.LE.T(NT)) GO TO 30                                        
C *** SOME EXTRAPOLATION OUTSIDE THE RANGE OF log Teff CONSIDERED IN
C *** THE TABLES IS PERMITTED.  TO AVOID ANY EXTRAPOLATIONS, THE NEXT
C *** TWO LINES SHOULD BE ACTIVATED.
C     WRITE(6,1005) GV,TEFF
C     STOP32
   30 R(1)=T(MT)                                                        
      R(2)=T(MT+1)                                                      
      R(3)=T(MT+2)                                                      
      R(4)=T(MT+3)                                                      
      CALL LGRAN4(R,AT,TEFF)                                            
      L=MG-1                                                            
      DO 35 I=1,4                                                       
      L=L+1                                                             
      O(I)=AT(1)*C(MT,L)+AT(2)*C(MT+1,L)+AT(3)*C(MT+2,L)+AT(4)*C(MT+3,L)
      P(I)=AT(1)*D(MT,L)+AT(2)*D(MT+1,L)+AT(3)*D(MT+2,L)+AT(4)*D(MT+3,L)
      Q(I)=AT(1)*E(MT,L)+AT(2)*E(MT+1,L)+AT(3)*E(MT+2,L)+AT(4)*E(MT+3,L)
   35 R(I)=AT(1)*F(MT,L)+AT(2)*F(MT+1,L)+AT(3)*F(MT+2,L)+AT(4)*F(MT+3,L)
C *** 4-POINT LAGRANGIAN INTERPOLATION IS USED TO DERIVE THE VALUES OF
C *** b-y, m1, c1, AND BC_V CORRESPONDING TO THE INPUT VALUES OF
C *** [Fe/H], log g, and log Teff.
      SBY=AG(1)*O(1)+AG(2)*O(2)+AG(3)*O(3)+AG(4)*O(4)
      SM1=AG(1)*P(1)+AG(2)*P(2)+AG(3)*P(3)+AG(4)*P(4)                   
      SC1=AG(1)*Q(1)+AG(2)*Q(2)+AG(3)*Q(3)+AG(4)*Q(4)                   
      BCV=AG(1)*R(1)+AG(2)*R(2)+AG(3)*R(3)+AG(4)*R(4)                  
C *** SPLINE INTERPOLATION IS USED TO FIND THE VALUE OF BC_V AT LOW
C *** TEMPERATURES.
      DO 38 I=1,NT                                                      
      X(I)=T(I)                                                         
   38 Y(I)=AG(1)*F(I,MG)+AG(2)*F(I,MG+1)+AG(3)*F(I,MG+2)+AG(4)*F(I,MG+3)
      CALL INTEP(TEFF,BCV,X,Y,NT,IER)                                  
      GO TO 50                                                          
C *** THE NEXT SECTION ASSUMES THAT THE HIGH-TEMPERATURE TABLES ARE THE
C *** RELEVANT ONES TO USE.
   40 IF(TEFF.LE.TT(NTT)) GO TO 42                                      
      WRITE(0,1005) GV,TEFF                                             
      STOP33                                                            
   42 NGM=NGG-1                                                         
      DO 43 I=3,NGM                                                     
      MG=I-2                                                            
      IF(GV.LE.GG(I)) GO TO 44                                          
   43 CONTINUE                                                          
      IF(GV.LE.GG(NGG)) GO TO 44                                        
C *** SOME EXTRAPOLATION OUTSIDE THE RANGE OF log g CONSIDERED IN THE
C *** TABLES IS PERMITTED.  TO AVOID ANY EXTRAPOLATIONS, THE NEXT TWO
C *** LINES SHOULD BE ACTIVATED.
C     WRITE(6,1005) GV,TEFF   
C     STOP34
   44 R(1)=GG(MG)                                                       
      R(2)=GG(MG+1)                                                     
      R(3)=GG(MG+2)                                                     
      R(4)=GG(MG+3)                                                     
      CALL LGRAN4(R,AG,GV)                                              
      NTM=NTT-1                                                         
      DO 45 I=3,NTM                                                     
      MT=I-2                                                            
      IF(TEFF.LE.TT(I)) GO TO 46                                        
   45 CONTINUE                                                          
      IF(TEFF.LE.TT(NTT)) GO TO 46                                      
C *** SOME EXTRAPOLATION OUTSIDE THE RANGE OF log Teff CONSIDERED IN
C *** THE TABLES IS PERMITTED.  TO AVOID ANY EXTRAPOLATIONS, THE NEXT
C *** TWO LINES SHOULD SHOULD BE ACTIVATED.
C     WRITE(6,1005) GV,TEFF
C     STOP35
   46 R(1)=TT(MT)                                                       
      R(2)=TT(MT+1)                                                     
      R(3)=TT(MT+2)                                                     
      R(4)=TT(MT+3)                                                     
      CALL LGRAN4(R,AT,TEFF)                                            
      L=MG-1                                                            
      DO 48 I=1,4                                                       
      L=L+1                                                             
      DO 48 K=1,NDXX                                                    
      H(I,K)=AT(1)*B(MT,L,K)+AT(2)*B(MT+1,L,K)+AT(3)*B(MT+2,L,K)+       
     1   AT(4)*B(MT+3,L,K)                                              
   48 CONTINUE                                                          
C *** 4-POINT LAGRANGIAN INTERPOLATION IS USED TO DERIVE THE VALUES OF
C *** b-y, m1, c1, AND BC_V CORRESPONDING TO THE INPUT VALUES OF
C *** [Fe/H], log g, and log Teff.
      SBY=AG(1)*H(1,1)+AG(2)*H(2,1)+AG(3)*H(3,1)+AG(4)*H(4,1)           
      SM1=AG(1)*H(1,2)+AG(2)*H(2,2)+AG(3)*H(3,2)+AG(4)*H(4,2)           
      SC1=AG(1)*H(1,3)+AG(2)*H(2,3)+AG(3)*H(3,3)+AG(4)*H(4,3)           
      BCV=AG(1)*H(1,4)+AG(2)*H(2,4)+AG(3)*H(3,4)+AG(4)*H(4,4)
      GO TO 50
   49 WRITE(0,1007)
      STOP          
   50 RETURN                                                            
      END

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE LGRAN4(X,A,XX)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      SAVE
      DIMENSION X(4),A(4)
      R1=(X(1)-X(2))*(X(1)-X(3))*(X(1)-X(4))
      R2=(X(2)-X(1))*(X(2)-X(3))*(X(2)-X(4))
      R3=(X(3)-X(1))*(X(3)-X(2))*(X(3)-X(4))
      R4=(X(4)-X(1))*(X(4)-X(2))*(X(4)-X(3))
      A(1)=((XX-X(2))*(XX-X(3))*(XX-X(4)))/R1
      A(2)=((XX-X(1))*(XX-X(3))*(XX-X(4)))/R2
      A(3)=((XX-X(1))*(XX-X(2))*(XX-X(4)))/R3
      A(4)=((XX-X(1))*(XX-X(2))*(XX-X(3)))/R4
      RETURN
      END

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE INTEP(XP,P,X,F,N,IER)
C *** PURPOSE:  To interpolate a function value P for a given argument
C     XP from a table of N values (X,F).  This is a spline interpolation
C     scheme based on Hermite polynomials.  The source is U.S. Airforce
C     Surveys in Geophysics No 272.
C *** USAGE:  For random values of XP
C               CALL INTEP(XP,P,X,F,N,IER)
C     or after the first call to INTEP with monotonically increasing or
C     decreasing values of XP consistent with the X vector
C               CALL EINTEP(XP,P,X,F,N,IER)
C     DESCRIPTION OF PARAMETERS:
C     XP  - the chosen argument value
C     P   - the resultant interpolated value
C     X   - the vector of independent values
C     F   - the vector of dependent values
C     N   - the number of points in the (X,F) vectors
C     IER - the resultant error parameter (set = 2 if XP is beyond
C           either extreme of X; in which case P is set to the value of
C           F at that extremum)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8 LP1,LP2,L1,L2
      DIMENSION F(*),X(*)
      IER=1
      IO=1
      IUP=0
      IF(X(2).LT.X(1)) IUP=1
      N1=N-1
      IF((XP.GE.X(N).AND.IUP.EQ.0).OR.(XP.LE.X(N).AND.IUP.EQ.1)) THEN
*    5  P=F(N)
       P=F(N)
*       GO TO 6
       IER=2
       RETURN
      ELSE IF((XP.LE.X(1).AND.IUP.EQ.0).OR.
     1   (XP.GE.X(1).AND.IUP.EQ.1)) THEN
       P=F(1)
*    6  IER=2
       IER=2
       RETURN
      ENDIF
      ENTRY EINTEP(XP,P,X,F,N,IER)
    8 DO 1 I=IO,N
      IF(XP.LT.X(I).AND.IUP.EQ.0) GO TO 2
      IF(XP.GT.X(I).AND.IUP.EQ.1) GO TO 2
    1 CONTINUE
      P=F(N)
      IER=2
      RETURN
*      GO TO 5
    2 I=I-1
      IF(I.EQ.IO-1) GO TO 4
      IO=I+1
      LP1=1./(X(I)-X(I+1))
      LP2=1./(X(I+1)-X(I))
      IF(I.EQ.1) FP1=(F(2)-F(1))/(X(2)-X(1))
      IF(I.EQ.1) GO TO 3
      FP1=(F(I+1)-F(I-1))/(X(I+1)-X(I-1))
    3 IF(I.GE.N1) FP2=(F(N)-F(N-1))/(X(N)-X(N-1))
      IF(I.EQ.N1) GO TO 4
      FP2=(F(I+2)-F(I))/(X(I+2)-X(I))
    4 XPI1=XP-X(I+1)
      XPI=XP-X(I)
      L1=XPI1*LP1
      L2=XPI*LP2
      P=F(I)*(1.-2.*LP1*XPI)*L1*L1+F(I+1)*(1.-2.*LP2*XPI1)*L2*L2+
     1   FP2*XPI1*L2*L2+FP1*XPI*L1*L1
      RETURN
      END      

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c END OF VANDENBERG & CLEM COLOR TRANSFORM SUBROUTINES
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c BEGIN AKIMA SPLINE SUBROUTINES
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE AKM_CSPL (n,x,f,c)

*  DESCRIPTION
*  Subroutine akm_cspl is based on the method of interpolation described
*  by Hiroshi Akima (1970 Journal of the Association for Computing Machinery,
*  Vol. 17, pp. 580-602) and also by Carl de Boor (1978 Applied Mathematical
*  Sciences, Vol. 27, "A Practical Guide to Splines").

*  This routine computes the Akima spline interpolant with the "not-a-knot
*  end conditions which are useful when there is no information about the
*  derivatives at the endpoints. In this case P(1)=P(2) and P(n-2)=P(n-1);
*  in other words, x(2) and x(n-1) are not knots even though they define
*  the endpoints of a polynomial piece --- the spline must pass through
*  these two points. Thus, it is a requirement that f''' be continuous
*  across x(2) and x(n-1). 

*  DIMENSIONS
*  The internal work arrays dd, w, and s are restricted to 2000 elements 
*  by the parameter nmx; consequently the data arrays x and f are also
*  limited to 2000 elements. The spline coefficients are contained in
*  the array c, a 4x2000 2 dimensional array.

*  CALL PARAMETERS
*     n      - number of data pairs
*     x(n)   - array, independent variable
*     f(n)   - array, dependent variable
*     c(4,n) - array, spline coefficients

      implicit double precision (a-h, o-z)
      parameter (nmx = 2000)
      parameter (zero = 0.0d0, two = 2.0d0, three = 3.0d0)
      dimension x(n),f(n),c(4,n) 
      dimension dd(nmx),w(nmx),s(nmx)

      nm1 = n - 1
      nm2 = n - 2
      nm3 = n - 3

*  Create two additional points outside the region of the spline by fitting
*  a quadratic to the three adjacent data points. This is a special feature
*  of the Akima spline.

      call akm_ends (x(1),x(2),x(3),f(1),f(2),f(3),f0,fm1)
      call akm_ends (x(n),x(nm1),x(nm2),f(n),f(nm1),f(nm2),fnp1,fnp2)

*  Compute the divided differences

      ddm1 = (f0 - fm1)/(x(2) - x(1))
      dd0 = (f(1) - f0)/(x(3) - x(2))

      do i = 1, nm1
        ip1 = i + 1
        dd(i) = (f(ip1) - f(i))/(x(ip1) - x(i))
      end do

      dd(n) = (f(n) - fnp1)/(x(nm2) - x(nm1))
      ddnp1 = (fnp1 - fnp2)/(x(nm1) - x(n))

*  Compute the Akima weights

      w0 = abs (dd0 - ddm1)
      w(1) = abs (dd(1) - dd0)

      do i = 2, nm1
        im1 = i - 1
        w(i) = abs(dd(i) - dd(im1))
      end do

      w(n) = abs (dd(n) - dd(nm1))
      wnp1 = abs (ddnp1 - dd(n))

*  Compute Akima slopes at interior knots

      if (w(2).eq.zero.and.w0.eq.zero) then
        s(1) = 5.0d-1*(dd(1) + dd0)
      else
        s(1) = (w(2)*dd0 + w0*dd(1))/(w0 + w(2))
      end if

      do i = 2, nm1
        im1 = i - 1
        ip1 = i + 1
        if (w(ip1).eq.zero.and.w(im1).eq.zero) then
          s(i) = 5.0d-1*(dd(i) + dd(im1))
        else
          s(i) = (w(ip1)*dd(im1) + w(im1)*dd(i))/(w(im1) + w(ip1))
        end if
      end do

      if (wnp1.eq.zero.and.w(nm1).eq.zero) then
        s(n) = 5.0d-1*(dd(n) + dd(nm1))
      else
        s(n) = (wnp1*dd(nm1) + w(nm1)*dd(n))/(w(nm1) + wnp1)
      end if

*  Spline coefficients for all the polynomial pieces

      do i = 1, nm1
        ip1 = i + 1
        dx = x(ip1) - x(i)
        c(1,i) = f(i)
        c(2,i) = s(i)
        c(3,i) = (three*dd(i) - two*s(i) - s(ip1))/dx
        c(4,i) = (s(ip1) + s(i) - two*dd(i))/dx/dx
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE AKM_ENDS(x1,x2,x3,f1,f2,f3,f00,f01)
      
*  DESCRIPTION
*  Subroutine akm_ends is based on the method of interpolation described
*  by Hiroshi Akima (1970 Journal of the Association for Computing Machinery,
*  Vol. 17, pp. 580-602).

*  This routine is required by the routine AKM_CSPL. It sets up the two 
*  additional external points required by the Akima method of constructing
*  a spline.

*  CALL PARAMETERS
*  (x1,f1), (x2,f2) and (x3,f3) are the known data pairs at the ends of the
*  region of interpolation. f00 and f01 are the computed values of the
*  interpolating polynomial external to the region of interpolation.

      implicit double precision (a-h,o-z)

      g0 = f1
      df21 = f2 - f1
      df31 = f3 - f1
      dx21 = x2 - x1
      dx31 = x3 - x1
      dx32 = x3 - x2
      den = dx21*dx32*dx31

      g1 = (df21*dx31*dx31 - df31*dx21*dx21)/den
      g2 = (df31*dx21 - df21*dx31)/den

      f00 = g0 - g1*dx32 + g2*dx32*dx32
      f01 = g0 - g1*dx31 + g2*dx31*dx31

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE AKM_EVAL (n,x,c,xp,fp,dfp)
      
*  DESCRIPTION
*  This routine is used to evaluate the spline defined by the
*  routine AKM_CSPL and its derivative at an interpolating point.

*  CALL PARAMETERS
*     n      - number of data pairs
*     x(n)   - array, independent variable
*     c(4,n) - array, spline coefficients
*     xp     - interpolating point
*     fp     - value of the spline at xp
*     dfp    - first derivative of the spline at xp   
      
      implicit double precision (a-h, o-z)
      dimension x(n), c(4,n)
      logical more

*  Statement function defines cubic spline interpolation

      csval(dx,a,b,e,d) = a +dx*(b + dx*(e + dx*d))

*  Statement function defines cubic spline derivative

      csder(dx,b,e,d) = b + dx*(2.0d0*e + dx*3.0d0*d)

*  Check direction of spline

      if (x(n).gt.x(1)) then

*  Check to see that x(1) < xp < x(n)

      if (xp.lt.x(1)) then
        write(*,'('' ***Test point is less than x(1)***'')')
        stop
      else if (xp.gt.x(n)) then
        write(*,'('' ***Test point is greater than x(n)***'')')
        stop
      end if

*  Find the polynomial piece containing the test point

      i = 2
      more = .true.
      do while (more)
        if (xp.lt.x(i).or.i.eq.n) then
          more = .false.
          nl = i - 1
        else 
          i = i + 1
        end if
      end do

      else

*  The spline is backwards
*  Check to see that x(n) < xp < x(1)

      if (xp.gt.x(1)) then
        write(*,'('' ***Test point is greater than x(1)***'')')
        stop
      else if (xp.lt.x(n)) then
        write(*,'('' ***Test point is less than x(n)***'')')
        stop
      end if

*  Find the polynomial piece containing the test point

      i = 2
      more = .true.
      do while (more)
        if (xp.gt.x(i).or.i.eq.n) then
          more = .false.
          nl = i - 1
        else 
          i = i + 1
        end if
      end do

      end if

*  Evaluate spline at the test point xp

      dx = xp - x(nl)
      fp = csval(dx,c(1,nl),c(2,nl),c(3,nl),c(4,nl))

*  Evaluate the derivative at the test point

      dfp = csder(dx,c(2,nl),c(3,nl),c(4,nl))
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   
      SUBROUTINE AKM_INT (a,b,n,x,c,sum)

*  DESCRIPTION
*  This routine is used to integrate the spline defined by the
*  routine AKM_CSPL.

*  CALL PARAMETERS
*     a      - lower bound
*     b      - upper bound
*     n      - number of data pairs
*     x(n)   - array, independent variable
*     c(4,n) - array, spline coefficients
*     sum     - value of the integral over the interval [a,b]

      implicit double precision (a-h, o-z)
      dimension x(n),c(4,n)

*  Statement function defined by cubic spline integration

      csint(dx,a,b,g,d)=abs(dx*(a+dx*(b/2.0d0+dx*(g/3.0d0+dx*d/4.0d0))))

*  Check direction of integration

      if (x(1).lt.x(n)) then

*  Forward integration
*  Check limits of integration

      if (a.lt.x(1)) then
        write(*,'('' ***Lower bound is less than x(1)***'')')
        stop
      else if (b.gt.x(n)) then
        write(*,'('' ***Upper bound is greater than x(n)***'')')
        stop
      end if

*  Find the polynomial piece containing the lower bound

      i = 2
      do while (a.gt.x(i))
        i = i + 1
      end do
      nl = i - 1

*  Find the polynomial piece containing the upper bound

      i = n - 1
      do while (b.lt.x(i))
        i = i - 1
      end do
      nu = i

      else

*  Backward integration
*  Check limits of integration

      if (a.gt.x(1)) then
        write(*,'('' ***Upper bound is greater than x(1)***'')')
        stop
      else if (b.lt.x(n)) then
        write(*,'('' ***Lower bound is less than x(n)***'')')
        stop
      end if

*  Find the polynomial piece containing the upper bound

      i = 2
      do while (a.lt.x(i))
        i = i + 1
      end do
      nl = i - 1

*  Find the polynomial piece containing the lower bound

      i = n - 1
      do while (b.gt.x(i))
        i = i - 1
      end do
      nu = i

      end if

*  Initialize the sum

      sum = 0.0d0

*  Subtract the portion of the 1st polynomial piece outside the lower bound

      dx = a - x(nl)
      if (dx.ne.0.0d0) then
        sum = sum - csint(dx,c(1,nl),c(2,nl),c(3,nl),c(4,nl))
      end if

*  Integrate from x(nl) to x(nu)

      if (nu.gt.nl) then
        nt = nu - 1
        do i = nl, nt
          dx = x(i+1) - x(i)
          sum = sum + csint(dx,c(1,i),c(2,i),c(3,i),c(4,i))
        end do
      end if

*  Add the portion of the last polynomial piece within the upper bound

      dx = b - x(nu)
      if (dx.ne.0.0d0) then
        sum = sum + csint(dx,c(1,nu),c(2,nu),c(3,nu),c(4,nu))
      end if

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE AKM_BISECTION(n,c,x,x1,x2,froot,xroot)

*  DESCRIPTION
*  This routine is used to invert the spline defined by the
*  routine AKM_CSPL in order to find a root by the bisection method.
*  The routine AKM_EVAL is called within this routine.

*  CALL PARAMETERS
*     n      - number of data pairs
*     c(4,n) - array, spline coefficients
*     x(n)   - array, independent variable
*     x1     - estimated lower limit of region to be searched
*     x2     - estimated upper limit of region to be searched
*     froot  - value of the spline at the location of the root
*     xroot  - calculated value of the root

      implicit double precision (a-h, o-z)
      dimension x(*), c(4,n)
      logical more

*  This routine must be preceded by a call to akm_cspl(n,x,f,c) to set up
*  the spline which is used to find the root.

      xa = x1
      xb = x2

      call akm_eval(n,x,c,xa,fa,dfdx)
      fa = froot - fa
      call akm_eval(n,x,c,xb,fb,dfdx)
      fb = froot - fb
      kl = 0

      do while (fa*fb.ge.0.0d0.and.kl.lt.20) 

        kl = kl + 1

        if (abs(fa).lt.abs(fb)) then
          xa = xa + 1.6d0*(xa - xb)
          xa = max(xa,x(1))
          call akm_eval(n,x,c,xa,fa,dfdx)
          fa = froot - fa
        else
          xb = xb + 1.6d0*(xb - xa)
          xb = min(xb,x(n))
          call akm_eval(n,x,c,xb,fb,dfdx)
          fb = froot - fb
        end if

       write(*,'(i4,1pd14.7,0pf7.3,1pd14.7,0pf7.3)') kl,xa,fa,xb,fb
      end do

      if (kl.eq.20) then
        write(*,'('' *** akm_bisection failed to find the region'')')
        write(*,'(''     containing the root after 20 iterations.'')')
        stop
      end if

      diff = abs(xb - xa)
      more = .true.
      kl = 0

      do while (more)

        kl = kl + 1
        xab = 5.0d-1*(xa + xb)
        call akm_eval(n,x,c,xab,fb,dfdx)
        fb = froot - fb

        if (fa*fb.le.0.0d0) then
          xb = xab
        else
          xa = xab
          fa = fb
        end if

      diff = abs(xb - xa)/abs(5.0d-1*(xa + xb))
      more = (diff.gt.1.0d-8.and.kl.lt.40)

      end do

      if (kl.eq.40) then
        write(*,'('' ***akm_bisection root did not converge within'')')
        write(*,'(''    one in 1.0d-8 after 40 iterations.'')')
        write(*,'(''    Converged within '',1pd12.5)') diff
      end if

      xroot = 5.0d-1*(xa + xb)

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c END AKIMA SPLINE SUBROUTINES
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c BEGIN UTILITY SUBROUTINES
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE PARSE (FILESPEC,FIELD,OUTSPEC,LENGTH)
      CHARACTER FILESPEC*(*),FIELD*(*),OUTSPEC*(*)

      IF (FIELD.EQ.'DEVICE') THEN

        NC1 = 1
        NC2 = INDEX(filespec,':')
        OUTSPEC = FILESPEC(NC1:NC2)
        LENGTH = NC2 - NC1 + 1

      ELSE IF (FIELD.EQ.'DIRECTORY') THEN

        NC1 = index(filespec,'[')
        NC2 = index(filespec,']')
        OUTSPEC = FILESPEC(NC1:NC2)
        LENGTH = NC2 - NC1 + 1

      ELSE IF (FIELD.EQ.'NAME') THEN

        NC1 = index(filespec,']')
        NC2 = index(filespec,'.')
        OUTSPEC = FILESPEC(NC1+1:NC2-1)
        LENGTH = NC2 - NC1 - 1

      ELSE IF (FIELD.EQ.'TYPE') THEN

        NC1 = index(filespec,'.')
        NC2 = index(filespec,';')
        OUTSPEC = FILESPEC(NC1:NC2-1)
        LENGTH = NC2 - NC1

      ELSE IF (FIELD.EQ.'VERSION') THEN

        NC1 = index(filespec,';')
        NC2 = index(filespec,' ')
        OUTSPEC = FILESPEC(NC1+1:NC2-1)
        LENGTH = NC2 - NC1 - 1

      ELSE

        OUTSPEC = 'INVALID FILE SPECIFICATION'
        LENGTH = 26

      END IF

      RETURN
      END
      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine str_trim(nl,str,nch)
      character*(*) str
      
      nch = nl
      do while (str(nch:nch).eq.' '.and.nch.gt.0)
        nch = nch - 1
      end do
      
      return
      end    
        
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        SUBROUTINE OPEN_FILE(NUNIT,FSTAT,FNME,ASK,PROMPT)
        INTEGER NUNIT
        CHARACTER*(*) FSTAT,FNME,PROMPT
        CHARACTER*1 ANS
        LOGICAL ASK, MORE

        MORE = .TRUE.

        DO WHILE (MORE)
          IF (ASK) THEN
           WRITE(*,'(T25,A,$)') PROMPT
            READ(*,'(A)') FNME
          END IF
          
          NCH=1
          DO WHILE (FNME(NCH:NCH).NE.' ')
            NCH = NCH + 1
          END DO
          NCH = NCH - 1

          IF (FSTAT.EQ.'IN'.OR.FSTAT.EQ.'in') THEN
            OPEN (NUNIT,IOSTAT=IERR,FILE=FNME,STATUS='OLD')
            IF (IERR.EQ.0) THEN
              MORE=.FALSE.
            ELSE
              WRITE(*,*)'*** ',FNME(1:NCH),' DOES NOT EXIST***'
            END IF

          ELSE IF(FSTAT.EQ.'OUT'.OR.FSTAT.EQ.'out') THEN
            OPEN (NUNIT,IOSTAT=IERR,FILE=FNME,STATUS='NEW')
            IF (IERR.EQ.0) THEN
              MORE=.FALSE.
            ELSE
              WRITE(*,*)'***',FNME(1:NCH),' ALREADY EXISTS***'
              WRITE(*,'(T25,''DO YOU WANT TO OVERWRITE? '',$)')
              READ(*,'(A)') ANS
              IF (ANS.EQ.'Y'.OR.ANS.EQ.'y') THEN
                OPEN(NUNIT,FILE=FNME,STATUS='OLD')
                CLOSE(NUNIT,STATUS='DELETE')
                OPEN(NUNIT,FILE=FNME,STATUS='NEW')
                MORE = .FALSE.
              ELSE
                WRITE(*,'(T25,''New file name: '',$)')
                READ(*,'(A)') FNME
              END IF
            END IF
          END IF

       END DO

       RETURN
       END

*23456789012345678901234567890123456789012345678901234567890123456789012
        SUBROUTINE IOFILE(NUNIT,FSTAT,DNME,DASK,DPROM,FNME,ASK,PROM)
        INTEGER NUNIT
        CHARACTER*(*) FSTAT,FNME,PROM,DPROM,DNME
        CHARACTER*40 TNME
        CHARACTER*1 ANS
        LOGICAL ASK, MORE, DASK, TASK

        IF (DASK) THEN
          LP = LEN(DPROM)
          LD = LEN(DNME)
          WRITE(*,'(T25,A,A,A,A,$)') DPROM(1:LP),'[',DNME(1:LD),']: '
          READ(*,'(A)') TNME
          IF (TNME.EQ.' ') THEN
            FNME = DNME
          ELSE
            FNME = TNME
          END IF
        END IF
        
        MORE = .TRUE.
        TASK = ASK

        DO WHILE (MORE)
          IF (TASK) THEN
           WRITE(*,'(T25,A,$)') PROM
            READ(*,'(A)') FNME
          END IF
          
          NCH=1
          DO WHILE (FNME(NCH:NCH).NE.' ')
            NCH = NCH + 1
          END DO
          NCH = NCH - 1

          IF (FSTAT.EQ.'IN'.or.fstat.eq.'in') THEN
            OPEN (NUNIT,iostat=ierr,FILE=FNME,STATUS='OLD')
            if (ierr.eq.0) then
              more=.false.
            else
              WRITE(*,*)'***',FNME(1:NCH),' DOES NOT EXIST***'
              if (.not.ask) task = .true.
            end if

          ELSE IF(FSTAT.EQ.'OUT'.or.fstat.eq.'out') THEN
            OPEN (NUNIT,iostat=ierr,FILE=FNME,STATUS='NEW')
            if (ierr.eq.0) then
              more=.false.
            else
              WRITE(*,*)'***',FNME(1:NCH),' ALREADY EXISTS***'
              write(*,'(''       Do you want to OVERWRITE? '',$)')
              read(*,'(a)') ans
              if (ans.eq.'y'.or.ans.eq.'Y') then
                open(nunit,file=fnme,status='old')
                close(nunit,status='delete')
                open(nunit,file=fnme,status='new')
                more = .false.
              else         
                WRITE(*,'(T25,''New file name: '',$)')
                READ(*,'(A)') FNME
              end if
            end if

          ELSE
            OPEN(NUNIT,iostat=ierr,STATUS='SCRATCH')
            if (ierr.eq.0) then
              more=.false.
            else
              WRITE(*,*)'***SCRATCH FILE ALREADY EXISTS***'
              write(*,'(''       Do you want to OVERWRITE? '',$)')
              read(*,'(a)') ans
              if (ans.eq.'y'.or.ans.eq.'Y') then
                open(nunit,file=fnme,status='old')
                close(nunit,status='delete')
                open(nunit,file=fnme,status='new')
                more = .false.
              end if
            end if

          END IF

       end do

       return
       END

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      LOGICAL FUNCTION YES(QUESTION)

*  Prompts question on 'stdout', checks for valid answer on 'stdin'
*  Returns true if answer valid and affirmative. Else false.

      CHARACTER QUESTION*(*),LETTER*1
      LOGICAL MORE

      MORE = .TRUE.
      DO WHILE (MORE)
         WRITE(*,'(T25,A,'' [Y/N] '',$)') QUESTION
         READ (*,'(A)') LETTER
         IF (LETTER.EQ.'Y'.OR.LETTER.EQ.'y') THEN
            YES = .TRUE.
            MORE = .FALSE.
         ELSE IF (LETTER.EQ.'N'.OR.LETTER.EQ.'n') THEN
            YES = .FALSE.
            MORE = .FALSE.
         ELSE
            WRITE(*,'('' **** INVALID ANSWER **** '')')
         END IF
      END DO

      RETURN
      END

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c END UTILITY SUBROUTINES
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c $$$$$$$$$$$$$$$$$$$$$$ END VRSUB SUBROUTINES $$$$$$$$$$$$$$$$$$$$$$$$$$
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
