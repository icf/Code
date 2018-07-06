* *******************************************************************
* COPYRIGHT (c) 1969 Hyprotech UK
* All rights reserved.
*
* None of the comments in this Copyright notice between the lines
* of asterisks shall be removed or altered in any way.
*
* This Package is intended for compilation without modification,
* so most of the embedded comments have been removed.
*
* ALL USE IS SUBJECT TO LICENCE. For full details of an HSL ARCHIVE
* Licence, see http://hsl.rl.ac.uk/archive/cou.html
*
* Please note that for an HSL ARCHIVE Licence:
*
* 1. The Package must not be copied for use by any other person.
*    Supply of any part of the library by the Licensee to a third party
*    shall be subject to prior written agreement between AEA
*    Hyprotech UK Limited and the Licensee on suitable terms and
*    conditions, which will include financial conditions.
* 2. All information on the Package is provided to the Licensee on the
*    understanding that the details thereof are confidential.
* 3. All publications issued by the Licensee that include results obtained
*    with the help of one or more of the Packages shall acknowledge the
*    use of the Packages. The Licensee will notify the Numerical Analysis
*    Group at Rutherford Appleton Laboratory of any such publication.
* 4. The Packages may be modified by or on behalf of the Licensee
*    for such use in research applications but at no time shall such
*    Packages or modifications thereof become the property of the
*    Licensee. The Licensee shall make available free of charge to the
*    copyright holder for any purpose all information relating to
*    any modification.
* 5. Neither CCLRC nor Hyprotech UK Limited shall be liable for any
*    direct or consequential loss or damage whatsoever arising out of
*    the use of Packages by the Licensee.
* *******************************************************************
*
*######DATE 3 February 1994
C       Toolpack tool decs employed.
C	Arg dimensions made *.
C	NWU removed.
C	DFLOAT -> DBLE.
C
      SUBROUTINE VA05AD(M,N,F,X,DSTEP,DMAX,ACC,MAXFUN,IPRINT,W)
      DOUBLE PRECISION ACC,DMAX,DSTEP
      INTEGER IPRINT,M,MAXFUN,N
      DOUBLE PRECISION F(*),W(*),X(*)
      DOUBLE PRECISION AD,ANMULT,AP,DD,DM,DMULT,DN,DPAR,DS,DSS,DTEST,DW,
     +                 FMIN,FNP,FSQ,PAR,PARM,PJ,PPAR,PRED,PTM,SP,SPP,SS,
     +                 ST,TINC
      INTEGER I,IC,IPC,IS,J,K,KK,KS,MAXC,MPN,NT,NTEST,NTPAR,NWC,NWD,NWF,
     +        NWI,NWT,NWV,NWW,NWX
      EXTERNAL CALFUN,MB11AD
      INTRINSIC DABS,DMAX1,DMIN1,DSQRT,FLOAT,IABS
      MAXC = 0
      MPN = M + N
      NT = N + 2
      NTEST = 0
      DTEST = FLOAT(N+N) - 0.5
      NWI = M*N
      NWX = NWI + MPN*N
      NWF = NWX + N
      NWC = NWF + M
      NWD = NWC + N
      NWW = NWD + N*N
      NWV = NWW + N
      NWT = NWV + M
      FMIN = -1.0D0
      DD = 0.0D0
      DSS = DSTEP*DSTEP
      DM = DMAX*DMAX
      PARM = DSQRT(ACC)/DMAX
      DPAR = 10.0D0*DM
      IS = 4
      IC = 0
      TINC = 1.0D0
      IF (IPRINT) 1,3,1
    1 WRITE (6,FMT=2)
    2 FORMAT (1H1,4X,'THE FOLLOWING OUTPUT IS PROVIDED BY SUBROUTINE',
     +       ' VA05AD',/,/)
      IPC = 0
      GO TO 3
    4 IF (MAXFUN-MAXC) 5,5,3
    5 IF (IPRINT) 139,140,139
  140 IPRINT = 2
      GO TO 19
  139 WRITE (6,FMT=6) MAXC
    6 FORMAT (/,/,/,5X,'ERROR RETURN FROM VA05A BECAUSE THERE HAVE BEEN'
     +       ,I5,' CALLS OF CALFUN')
      GO TO 7
    3 MAXC = MAXC + 1
      CALL CALFUN(M,N,F,X)
      FSQ = 0.0D0
      DO 8 I = 1,M
        FSQ = FSQ + F(I)*F(I)
    8 CONTINUE
      GO TO (9,10,9,10) IS
    9 IF (FSQ-FMIN) 11,12,12
   12 IF (DD-DSS) 13,13,10
   13 NTEST = NTEST - 1
      IF (NTEST) 14,14,10
   14 IF (IPRINT) 15,17,15
   17 IPRINT = 1
      GO TO 19
   15 WRITE (6,FMT=16)
   16 FORMAT (/,/,/,5X,'ERROR RETURN FROM VA05A BECAUSE F(X) NO LONGER',
     +       ' DECREASES',/,/,5X,
     +       'THIS MAY BE DUE TO THE VALUES OF DSTEP',
     +       ' AND ACC, OR TO LOSS OF RANK IN THE JACOBIAN MATRIX')
    7 IF (IPRINT) 18,19,18
   18 WRITE (6,FMT=20) MAXC
   20 FORMAT (/,/,/,5X,'THE FINAL SOLUTION CALCULATED BY VA05A REQUIRED'
     +       ,I5,' CALLS OF CALFUN, AND IS')
      WRITE (6,FMT=21) (I,W(NWX+I),I=1,N)
   21 FORMAT (/,/,4X,'I',7X,'X(I)',10X,'I',7X,'X(I)',10X,'I',7X,'X(I)',
     +       10X,'I',7X,'X(I)',10X,'I',7X,'X(I)',/,/,5 (I5,D17.8))
      WRITE (6,FMT=22) (I,W(NWF+I),I=1,M)
   22 FORMAT (/,/,4X,'I',7X,'F(I)',10X,'I',7X,'F(I)',10X,'I',7X,'F(I)',
     +       10X,'I',7X,'F(I)',10X,'I',7X,'F(I)',/,/,5 (I5,D17.8))
      WRITE (6,FMT=23) FMIN
   23 FORMAT (/,5X,'THE SUM OF SQUARES IS',D17.8)
   19 DO 135 I = 1,N
        X(I) = W(NWX+I)
  135 CONTINUE
      DO 136 I = 1,M
        F(I) = W(NWF+I)
  136 CONTINUE
      RETURN
   11 NTEST = NT
   10 IF (IABS(IPRINT)-1) 39,38,40
   38 WRITE (6,FMT=41) MAXC
   41 FORMAT (/,/,/,5X,'AT THE',I5,' TH CALL OF CALFUN WE HAVE')
   42 WRITE (6,FMT=21) (I,X(I),I=1,N)
      WRITE (6,FMT=23) FSQ
      IF (IPRINT) 39,39,142
  142 WRITE (6,FMT=22) (I,F(I),I=1,M)
      GO TO 39
   40 IPC = IPC - 1
      IF (IPC) 43,43,39
   43 WRITE (6,FMT=44) MAXC
   44 FORMAT (/,/,/,5X,'THE BEST ESTIMATE AFTER',I5,
     +       ' CALLS OF CALFUN IS')
      IPC = IABS(IPRINT)
      IF (FSQ-FMIN) 42,45,45
   45 IF (FMIN) 42,46,46
   46 WRITE (6,FMT=21) (I,W(NWX+I),I=1,N)
      WRITE (6,FMT=23) FMIN
      IF (IPRINT) 39,39,143
  143 WRITE (6,FMT=22) (I,W(NWF+I),I=1,M)
   39 GO TO (49,47,47,48) IS
   48 IF (IC) 50,50,51
   50 DO 52 I = 1,N
        W(NWX+I) = X(I)
   52 CONTINUE
      GO TO 54
   51 K = IC
      DO 55 I = 1,M
        W(K) = (F(I)-W(NWF+I))/DSTEP
        K = K + N
   55 CONTINUE
      IF (FMIN-FSQ) 56,56,57
   56 X(IC) = W(NWX+IC)
      GO TO 58
   57 W(NWX+IC) = X(IC)
   54 DO 53 I = 1,M
        W(NWF+I) = F(I)
   53 CONTINUE
      FMIN = FSQ
   58 IC = IC + 1
      IF (IC-N) 59,59,60
   59 X(IC) = W(NWX+IC) + DSTEP
      GO TO 3
   60 K = NWD
      DO 61 I = 1,N
        DO 62 J = 1,N
          K = K + 1
          W(K) = 0.0D0
   62   CONTINUE
        W(K+I-N) = 1.D0
        W(NWC+I) = 1.0D0 + DBLE(N-I)
   61 CONTINUE
   24 PAR = PARM
   25 PPAR = PAR*PAR
      NTPAR = 0
   63 KK = 0
      K = NWI + NWI
      DO 26 I = 1,N
        DO 141 J = 1,M
          KK = KK + 1
          W(KK+NWI) = W(KK)
  141   CONTINUE
        DO 27 J = 1,N
          K = K + 1
          W(K) = 0.0D0
   27   CONTINUE
        W(K+I-N) = PAR
   26 CONTINUE
      CALL MB11AD(N,MPN,W(NWI+1),N,W(NWW+1))
   64 IF (FMIN-ACC) 7,7,65
   65 DS = 0.0D0
      DN = 0.0D0
      SP = 0.0D0
      DO 66 I = 1,N
        X(I) = 0.0D0
        F(I) = 0.0D0
        K = I
        DO 67 J = 1,M
          X(I) = X(I) - W(K)*W(NWF+J)
          F(I) = F(I) - W(NWI+K)*W(NWF+J)
          K = K + N
   67   CONTINUE
        DS = DS + X(I)*X(I)
        DN = DN + F(I)*F(I)
        SP = SP + X(I)*F(I)
   66 CONTINUE
      PRED = SP + SP
      DMULT = 0.0D0
      K = 0
      DO 68 I = 1,M
        AP = 0.0D0
        AD = 0.0D0
        DO 69 J = 1,N
          K = K + 1
          AP = AP + W(K)*F(J)
          AD = AD + W(K)*X(J)
   69   CONTINUE
        PRED = PRED - AP*AP
        DMULT = DMULT + AD*AD
   68 CONTINUE
      IF (DN-DM) 28,28,29
   28 AP = DSQRT(DN)
      IF (PRED+2.0D0*PPAR*AP* (DMAX-AP)-ACC) 7,7,70
   29 IF (PRED+PPAR* (DM-DN)-ACC) 7,7,70
   70 DMULT = DS/DMULT
      DS = DS*DMULT*DMULT
   71 IS = 2
      IF (DN-DD) 72,72,73
   72 IF (PAR-PARM) 30,30,24
   30 DD = DMAX1(DN,DSS)
      DS = 0.25D0*DN
      TINC = 1.0D0
      IF (DN-DSS) 74,132,132
   74 IS = 3
      GO TO 103
   73 IF (DN-DPAR) 31,31,32
   31 NTPAR = 0
      GO TO 33
   32 IF (NTPAR) 34,34,35
   34 NTPAR = 1
      PTM = DN
      GO TO 33
   35 NTPAR = NTPAR + 1
      PTM = DMIN1(PTM,DN)
      IF (NTPAR-NT) 33,36,36
   36 PAR = PAR* (PTM/DM)**0.25D0
      IF (6.0D0*DD-DM) 137,25,25
  137 AP = DSQRT(PRED/DN)
      IF (AP-PAR) 25,25,138
  138 PAR = DMIN1(AP,PAR* (DM/ (6.0D0*DD))**0.25D0)
      GO TO 25
   33 IF (DS-DD) 75,76,76
   76 IF (DD) 77,77,78
   77 DD = DMIN1(DM,DS)
      IF (DD-DSS) 79,78,78
   79 DD = DSS
      GO TO 71
   78 ANMULT = 0.D0
      DMULT = DMULT*DSQRT(DD/DS)
      GO TO 80
   75 SP = SP*DMULT
      ANMULT = (DD-DS)/ ((SP-DS)+DSQRT((SP-DD)**2+ (DN-DD)* (DD-DS)))
      DMULT = DMULT* (1.0D0-ANMULT)
   80 DN = 0.0D0
      SP = 0.0D0
      DO 81 I = 1,N
        F(I) = DMULT*X(I) + ANMULT*F(I)
        DN = DN + F(I)*F(I)
        SP = SP + F(I)*W(NWD+I)
   81 CONTINUE
      DS = 0.25D0*DN
      IF (W(NWC+1)-DTEST) 132,132,82
   82 IF (SP*SP-DS) 83,132,132
   83 DO 84 I = 1,N
        X(I) = W(NWX+I) + DSTEP*W(NWD+I)
        W(NWC+I) = W(NWC+I+1) + 1.D0
   84 CONTINUE
      W(NWD) = 1.D0
      IF (N.LE.1) GO TO 4
      DO 85 I = 1,N
        K = NWD + I
        SP = W(K)
        DO 86 J = 2,N
          W(K) = W(K+N)
          K = K + N
   86   CONTINUE
        W(K) = SP
   85 CONTINUE
      GO TO 4
  132 IF (N.GE.2) GO TO 153
      IS = 1
      GO TO 152
  153 SP = 0D0
      K = NWD
      DW = 0.0D0
      DO 87 I = 1,N
        X(I) = DW
        DW = 0.0D0
        DO 88 J = 1,N
          K = K + 1
          DW = DW + F(J)*W(K)
   88   CONTINUE
        GO TO (89,90) IS
   90   W(NWC+I) = W(NWC+I) + 1.D0
        SP = SP + DW*DW
        IF (SP-DS) 87,87,91
   91   IS = 1
        KK = I
        X(1) = DW
        GO TO 92
   89   X(I) = DW
   92   W(NWC+I) = W(NWC+I+1) + 1.D0
   87 CONTINUE
      W(NWD) = 1.D0
      IF (KK-1) 93,93,94
   94 KS = NWC + KK*N
      DO 95 I = 1,N
        K = KS + I
        SP = W(K)
        DO 96 J = 2,KK
          W(K) = W(K-N)
          K = K - N
   96   CONTINUE
        W(K) = SP
   95 CONTINUE
   93 DO 97 I = 1,N
        W(NWW+I) = 0.D0
   97 CONTINUE
      SP = X(1)*X(1)
      K = NWD
      DO 98 I = 2,N
        DS = DSQRT(SP* (SP+X(I)*X(I)))
        DW = SP/DS
        DS = X(I)/DS
        SP = SP + X(I)*X(I)
        DO 99 J = 1,N
          K = K + 1
          W(NWW+J) = W(NWW+J) + X(I-1)*W(K)
          W(K) = DW*W(K+N) - DS*W(NWW+J)
   99   CONTINUE
   98 CONTINUE
      SP = 1.0D0/DSQRT(DN)
      DO 100 I = 1,N
        K = K + 1
        W(K) = SP*F(I)
  100 CONTINUE
  152 FNP = 0.0D0
      K = 0
      DO 101 I = 1,M
        W(NWW+I) = W(NWF+I)
        DO 102 J = 1,N
          K = K + 1
          W(NWW+I) = W(NWW+I) + W(K)*F(J)
  102   CONTINUE
        FNP = FNP + W(NWW+I)**2
  101 CONTINUE
  103 DO 104 I = 1,N
        X(I) = W(NWX+I) + F(I)
  104 CONTINUE
      GO TO 4
   49 DMULT = 0.9D0*FMIN + 0.1D0*FNP - FSQ
      IF (DMULT) 105,108,108
  105 DD = DMAX1(DSS,0.25D0*DD)
      TINC = 1.0D0
      IF (FSQ-FMIN) 106,107,107
  108 SP = 0.0D0
      SS = 0.0D0
      DO 109 I = 1,M
        SP = SP + DABS(F(I)* (F(I)-W(NWW+I)))
        SS = SS + (F(I)-W(NWW+I))**2
  109 CONTINUE
      PJ = 1.0D0 + DMULT/ (SP+DSQRT(SP*SP+DMULT*SS))
      SP = DMIN1(4.0D0,TINC,PJ)
      TINC = PJ/SP
      DD = DMIN1(DM,SP*DD)
      GO TO 106
   47 IF (FSQ-FMIN) 106,110,110
  106 FMIN = FSQ
      DO 111 I = 1,N
        SP = X(I)
        X(I) = W(NWX+I)
        W(NWX+I) = SP
  111 CONTINUE
      DO 112 I = 1,M
        SP = F(I)
        F(I) = W(NWF+I)
        W(NWF+I) = SP
  112 CONTINUE
  110 GO TO (107,107,113) IS
  113 IS = 2
      IF (FMIN-ACC) 7,7,83
  107 DS = 0.0D0
      DO 114 I = 1,N
        X(I) = X(I) - W(NWX+I)
        DS = DS + X(I)*X(I)
  114 CONTINUE
      DO 115 I = 1,M
        F(I) = F(I) - W(NWF+I)
  115 CONTINUE
      K = NWI
      SS = 0.0D0
      DO 116 I = 1,MPN
        SP = 0.0D0
        DO 117 J = 1,N
          K = K + 1
          SP = SP + W(K)*X(J)
  117   CONTINUE
        W(NWV+I) = SP
        SS = SS + SP*SP
  116 CONTINUE
      DO 118 I = 1,N
        ST = 0.0D0
        K = NWI + I
        DO 119 J = 1,MPN
          ST = ST + W(K)*W(J+NWV)
          K = K + N
  119   CONTINUE
        ST = ST/SS
        K = NWI + I
        DO 120 J = 1,MPN
          W(K) = W(K) - ST*W(J+NWV)
          K = K + N
  120   CONTINUE
        ST = PPAR*X(I)
        K = I
        DO 121 J = 1,M
          ST = ST + W(K)*F(J)
          K = K + N
  121   CONTINUE
        W(NWW+I) = ST
  118 CONTINUE
      IC = 0
      K = 0
      KK = NWI
      SP = 0.0D0
      SPP = 0.0D0
      DO 122 I = 1,M
        SS = F(I)
        ST = F(I)
        DO 123 J = 1,N
          IC = IC + 1
          KK = KK + 1
          SS = SS - W(IC)*X(J)
          ST = ST - W(KK)*W(NWW+J)
  123   CONTINUE
        SS = SS/DS
        W(NWV+I) = ST
        SP = SP + F(I)*ST
        SPP = SPP + ST*ST
        DO 124 J = 1,N
          K = K + 1
          W(K) = W(K) + SS*X(J)
  124   CONTINUE
  122 CONTINUE
      DO 125 I = 1,N
        ST = PAR*X(I)
        DO 126 J = 1,N
          KK = KK + 1
          ST = ST - W(KK)*W(NWW+J)
  126   CONTINUE
        W(NWT+I) = ST
        SP = SP + PAR*X(I)*ST
        SPP = SPP + ST*ST
  125 CONTINUE
      IF (0.01D0*SPP-DABS(SP-SPP)) 63,63,127
  127 DO 128 I = 1,N
        K = NWI + I
        ST = X(I)
        DO 129 J = 1,M
          ST = ST - W(K)*F(J)
          K = K + N
  129   CONTINUE
        SS = 0.0D0
        DO 130 J = 1,N
          SS = SS + W(K)*X(J)
          K = K + N
  130   CONTINUE
        ST = (ST-PAR*SS)/SP
        K = NWI + I
        DO 131 J = 1,MPN
          W(K) = W(K) + ST*W(NWV+J)
          K = K + N
  131   CONTINUE
  128 CONTINUE
      GO TO 64
      END



* *******************************************************************
* COPYRIGHT (c) 1969 Hyprotech UK
* All rights reserved.
*
* None of the comments in this Copyright notice between the lines
* of asterisks shall be removed or altered in any way.
*
* This Package is intended for compilation without modification,
* so most of the embedded comments have been removed.
*
* ALL USE IS SUBJECT TO LICENCE. For full details of an HSL ARCHIVE
* Licence, see http://hsl.rl.ac.uk/archive/cou.html
*
* Please note that for an HSL ARCHIVE Licence:
*
* 1. The Package must not be copied for use by any other person.
*    Supply of any part of the library by the Licensee to a third party
*    shall be subject to prior written agreement between AEA
*    Hyprotech UK Limited and the Licensee on suitable terms and
*    conditions, which will include financial conditions.
* 2. All information on the Package is provided to the Licensee on the
*    understanding that the details thereof are confidential.
* 3. All publications issued by the Licensee that include results obtained
*    with the help of one or more of the Packages shall acknowledge the
*    use of the Packages. The Licensee will notify the Numerical Analysis
*    Group at Rutherford Appleton Laboratory of any such publication.
* 4. The Packages may be modified by or on behalf of the Licensee
*    for such use in research applications but at no time shall such
*    Packages or modifications thereof become the property of the
*    Licensee. The Licensee shall make available free of charge to the
*    copyright holder for any purpose all information relating to
*    any modification.
* 5. Neither CCLRC nor Hyprotech UK Limited shall be liable for any
*    direct or consequential loss or damage whatsoever arising out of
*    the use of Packages by the Licensee.
* *******************************************************************
*
*######DATE 7 Dec 1992
C       Toolpack tool decs employed.
C       Make ZERO and ONE PARAMETER.
C       Change DFLOAT to DBLE.
C       Change arg dimensions to *.
C       Remove MB11CD reference from MB11AD
C       SAVE statements added.
C
C  EAT 21/6/93 EXTERNAL statement put in for block data so will work on VAXs.
C
C
      SUBROUTINE MB11AD(M,N,A,IA,W)
      DOUBLE PRECISION ONE,ZERO
      PARAMETER (ONE=1.0D0,ZERO=0.0D0)
      INTEGER IA,M,N
      DOUBLE PRECISION A(IA,*),W(*)
      DOUBLE PRECISION AKK,BSQ,RMAX,SIGMA,SUM,WKK
      INTEGER I,IR,J,KK,KP,MMK,N1,N2,N3,N4,N5,N6,NCW,NRW
      EXTERNAL MB11DD,MB11ED,MB11FD
      EXTERNAL MB11CD
      INTRINSIC DABS,DBLE,DSIGN,DSQRT,IDINT
      COMMON /MB11BD/LP
      INTEGER LP
      SAVE /MB11BD/
      NRW = M
      NCW = M + M
      DO 1 I = 1,M
        N1 = NRW + I
        W(N1) = 0.5D0 + DBLE(I)
    1 CONTINUE
      DO 2 I = 1,N
        N1 = NCW + I
        W(N1) = 0.5D0 + DBLE(I)
    2 CONTINUE
      KK = 1
    3 RMAX = 0.0D0
      DO 4 I = KK,M
        SUM = 0.0D0
        DO 5 J = KK,N
          SUM = SUM + A(I,J)**2
    5   CONTINUE
        IF (RMAX-SUM) 6,4,4
    6   RMAX = SUM
        IR = I
    4 CONTINUE
      IF (RMAX.EQ.0.0D0) GO TO 81
      IF (IR-KK) 7,7,8
    8 N3 = NRW + KK
      SUM = W(N3)
      N4 = NRW + IR
      W(N3) = W(N4)
      W(N4) = SUM
      DO 9 J = 1,N
        SUM = A(KK,J)
        A(KK,J) = A(IR,J)
        A(IR,J) = SUM
    9 CONTINUE
    7 RMAX = 0.0D0
      SUM = 0.0D0
      DO 10 J = KK,N
        SUM = SUM + A(KK,J)**2
        IF (RMAX-DABS(A(KK,J))) 11,10,10
   11   RMAX = DABS(A(KK,J))
        IR = J
   10 CONTINUE
      IF (IR-KK) 12,12,13
   13 N5 = NCW + KK
      RMAX = W(N5)
      N6 = NCW + IR
      W(N5) = W(N6)
      W(N6) = RMAX
      DO 14 I = 1,M
        RMAX = A(I,KK)
        A(I,KK) = A(I,IR)
        A(I,IR) = RMAX
   14 CONTINUE
   12 SIGMA = DSQRT(SUM)
      BSQ = DSQRT(SUM+SIGMA*DABS(A(KK,KK)))
      W(KK) = DSIGN(SIGMA+DABS(A(KK,KK)),A(KK,KK))/BSQ
      A(KK,KK) = -DSIGN(SIGMA,A(KK,KK))
      KP = KK + 1
      IF (KP-N) 15,15,16
   15 DO 17 J = KP,N
        A(KK,J) = A(KK,J)/BSQ
   17 CONTINUE
      IF (KP-M) 18,18,16
   18 WKK = W(KK)
      CALL MB11DD(A(KK+1,KK+1),A(KK+1,KK),A(KK,KK+1),IA,WKK,M-KK,N-KK)
      KK = KP
      GO TO 3
   16 KK = M
      KP = M + 1
      SUM = -W(M)/A(M,M)
      IF (N-M) 33,33,34
   34 DO 35 J = KP,N
        A(M,J) = SUM*A(M,J)
   35 CONTINUE
   33 A(M,M) = 1.0D0/A(M,M) + SUM*W(M)
   36 KP = KK
      KK = KP - 1
      IF (KK) 37,37,38
   38 WKK = W(KK)
      CALL MB11ED(A(KK+1,KK+1),A(KK,KK+1),IA,W(KK+1),WKK,M-KK,N-KK)
      AKK = ONE/A(KK,KK)
      CALL MB11FD(A(KK+1,KK+1),A(KK,KK+1),A(KK+1,KK),IA,WKK,AKK,M-KK,
     +            N-KK)
      SUM = 1.0D0 - WKK**2
      DO 44 I = KP,M
        SUM = SUM - A(I,KK)*W(I)
        A(I,KK) = W(I)
   44 CONTINUE
      A(KK,KK) = SUM/A(KK,KK)
      GO TO 36
   37 DO 45 I = 1,M
   46   N1 = NRW + I
        IR = IDINT(W(N1))
        IF (I-IR) 47,45,45
   47   SUM = W(N1)
        N2 = NRW + IR
        W(N1) = W(N2)
        W(N2) = SUM
        DO 48 J = 1,N
          SUM = A(I,J)
          A(I,J) = A(IR,J)
          A(IR,J) = SUM
   48   CONTINUE
        GO TO 46
   45 CONTINUE
      DO 49 J = 1,N
   50   N1 = NCW + J
        IR = IDINT(W(N1))
        IF (J-IR) 51,49,49
   51   SUM = W(N1)
        N2 = NCW + IR
        W(N1) = W(N2)
        W(N2) = SUM
        DO 52 I = 1,M
          SUM = A(I,J)
          A(I,J) = A(I,IR)
          A(I,IR) = SUM
   52   CONTINUE
        GO TO 50
   49 CONTINUE
   80 RETURN
   81 IF (LP.LE.0) GO TO 80
      MMK = M - KK
      WRITE (LP,FMT=82) MMK
   82 FORMAT (1H0,22H *** MB11AD ERROR *** ,I3,8H REDUCED,
     +       22H ROWS FOUND TO BE ZERO)
      STOP
      END
      BLOCK DATA MB11CD
      COMMON /MB11BD/LP
      INTEGER LP
      SAVE /MB11BD/
      DATA LP/6/
      END
C@PROCESS DIRECTIVE('IBMD')
      SUBROUTINE MB11DD(A,B,C,IA,WKK,MKK,NKK)
      DOUBLE PRECISION WKK
      INTEGER IA,MKK,NKK
      DOUBLE PRECISION A(IA,NKK),B(*),C(IA,NKK)
      DOUBLE PRECISION SUM
      INTEGER I,J
      DO 19 I = 1,MKK
        SUM = WKK*B(I)
CIBMD PREFER SCALAR
        DO 20 J = 1,NKK
          SUM = SUM + C(1,J)*A(I,J)
   20   CONTINUE
        SUM = -SUM
        B(I) = B(I) + SUM*WKK
CIBMD PREFER SCALAR
        DO 21 J = 1,NKK
          A(I,J) = A(I,J) + SUM*C(1,J)
   21   CONTINUE
   19 CONTINUE
      RETURN
      END
C@PROCESS DIRECTIVE('IBMD')
      SUBROUTINE MB11ED(A,B,IA,W,WKK,MKK,NKK)
      DOUBLE PRECISION WKK
      INTEGER IA,MKK,NKK
      DOUBLE PRECISION A(IA,NKK),B(IA,NKK),W(*)
      DOUBLE PRECISION SUM
      INTEGER I,J
      DO 39 I = 1,MKK
        SUM = 0.0D0
CIBMD PREFER SCALAR
        DO 40 J = 1,NKK
          SUM = SUM + B(1,J)*A(I,J)
   40   CONTINUE
        SUM = -SUM
CIBMD PREFER SCALAR
        DO 41 J = 1,NKK
          A(I,J) = A(I,J) + SUM*B(1,J)
   41   CONTINUE
        W(I) = SUM*WKK
   39 CONTINUE
      RETURN
      END
C@PROCESS DIRECTIVE('IBMD')
      SUBROUTINE MB11FD(A,B,C,IA,WKK,AKK,MKK,NKK)
      DOUBLE PRECISION AKK,WKK
      INTEGER IA,MKK,NKK
      DOUBLE PRECISION A(IA,*),B(IA,*),C(*)
      DOUBLE PRECISION SUM
      INTEGER I,J
      DO 42 J = 1,NKK
        SUM = WKK*B(1,J)
        DO 43 I = 1,MKK
          SUM = SUM + C(I)*A(I,J)
   43   CONTINUE
        SUM = -SUM
        B(1,J) = SUM*AKK
   42 CONTINUE
      RETURN
      END

