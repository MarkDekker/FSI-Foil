
      SUBROUTINE READI(N,IVAR,ERROR)
      DIMENSION IVAR(N)
      LOGICAL ERROR
C--------------------------------------------------
C     Reads N integer variables, leaving unchanged 
C     if only <return> is entered.
C--------------------------------------------------
      DIMENSION IVTMP(40)
      CHARACTER*80 LINE
C
      READ(*,1000) LINE
 1000 FORMAT(A80)
C
      DO 10 I=1, N
        IVTMP(I) = IVAR(I)
 10   CONTINUE
C
      NTMP = 40
      CALL GETINT(LINE,IVTMP,NTMP,ERROR)
C
      IF(ERROR) RETURN
C
      DO 20 I=1, N
        IVAR(I) = IVTMP(I)
 20   CONTINUE
C
      RETURN
      END ! READI



      SUBROUTINE READR(N,VAR,ERROR)
      DIMENSION VAR(N)
      LOGICAL ERROR
C-------------------------------------------------
C     Reads N real variables, leaving unchanged 
C     if only <return> is entered.
C-------------------------------------------------
      DIMENSION VTMP(40)
      CHARACTER*80 LINE
C
      READ(*,1000) LINE
 1000 FORMAT(A80)
C
      DO 10 I=1, N
        VTMP(I) = VAR(I)
 10   CONTINUE
C
      NTMP = 40
      CALL GETFLT(LINE,VTMP,NTMP,ERROR)
C
      IF(ERROR) RETURN
C
      DO 20 I=1, N
        VAR(I) = VTMP(I)
 20   CONTINUE
C
      RETURN
      END ! READR



      SUBROUTINE GETINT(INPUT,INUM,NI,ERROR)
      CHARACTER*(*) INPUT
      INTEGER INUM(*)
      LOGICAL ERROR
C----------------------------------------------------------------
C     Parses character string INPUT into an array
C     of integer numbers returned in INUM(1..NI).
C
C     Will attempt to extract no more than NI numbers,
C     unless NI = 0, in which case all numbers present 
C     in INPUT will be extracted.
C
C     NI returns how many numbers were actually extracted.
C----------------------------------------------------------------
C
C---- number of characters to be examined
      ILEN = LEN(INPUT)
C
C---- ignore everything after a "!" character
      K = INDEX(INPUT,'!')
      IF(K.GT.0) ILEN = K-1
C
C---- set limit on numbers to be read
      NINP = NI
      IF(NINP.EQ.0) NINP = ILEN/2 + 1
C
      NI = 0
C
      IF(ILEN.EQ.0) RETURN
C
C---- extract numbers
      N = 0
      K = 1
      DO 10 IPASS=1, ILEN
C------ find next space (pretend there's one after the end of the string)
        KSPACE = INDEX(INPUT(K:ILEN),' ') + K - 1
        IF(KSPACE.EQ.K-1) KSPACE = ILEN + 1
C
        IF(KSPACE.EQ.K) THEN
C------- just skip this space
         K = K+1
         GO TO 9
        ENDIF
C
C------ also find next comma
        KCOMMA = INDEX(INPUT(K:ILEN),',') + K - 1
        IF(KCOMMA.EQ.K-1) KCOMMA = ILEN + 1
C
C------ space is farther down, so we ran into something...
        N = N+1
C
C------ bug out early if no more numbers are to be read
        IF(N.GT.NINP) GO TO 11
C
C------ set ending delimiter position for this number
        KDELIM = MIN(KSPACE,KCOMMA)
C
        IF(K.EQ.KDELIM) THEN
C------- nothing but a comma... just keep looking
         K = K+1
         GO TO 9
        ENDIF
C
C------ whatever we have, it is in substring K:KEND
        KEND = KDELIM - 1
C
C------ search for floating-point number indicator in substring
        KFLOAT = MAX( INDEX(INPUT(K:KEND),'.'),
     &                INDEX(INPUT(K:KEND),'E'),
     &                INDEX(INPUT(K:KEND),'e'),
     &                INDEX(INPUT(K:KEND),'D'),
     &                INDEX(INPUT(K:KEND),'d')  ) + K - 1
C
        IF(KFLOAT.EQ.K) THEN
C------- only ".eEdD" was input (ugh!)... pretend it's zero anyway
         INUM(N) = 0
        ELSE
C------- read normally, ignoring any stuff after floating-point indicator
         IF(KFLOAT.GT.K) KEND = KFLOAT - 1
         READ(INPUT(K:KEND),*,ERR=20) INUM(N)
        ENDIF
C
        NI = N
C
C------ keep looking after delimiter
        K = KDELIM + 1
C
  9     IF(K.GE.ILEN) GO TO 11
 10   CONTINUE
C
C---- normal return
 11   CONTINUE
      ERROR = .FALSE.
      RETURN
C
C---- bzzzt !!!
 20   CONTINUE
ccc   WRITE(*,*) 'GETINT: List-directed read error.'
      ERROR = .TRUE.
      RETURN
      END ! GETINT



      SUBROUTINE GETFLT(INPUT,RNUM,NR,ERROR)
      CHARACTER*(*) INPUT
      REAL RNUM(*)
      LOGICAL ERROR
C----------------------------------------------------------------
C     Parses character string INPUT into an array
C     of real numbers returned in RNUM(1..NR).
C
C     Will attempt to extract no more than NR numbers,
C     unless NR = 0, in which case all numbers present 
C     in INPUT will be extracted.
C
C     NR returns how many numbers were actually extracted.
C----------------------------------------------------------------
C
C---- number of characters to be examined
      ILEN = LEN(INPUT)
C
C---- ignore everything after a "!" character
      K = INDEX(INPUT,'!')
      IF(K.GT.0) ILEN = K-1
C
C---- set limit on numbers to be read
      NINP = NR
      IF(NINP.EQ.0) NINP = ILEN/2 + 1
C
      NR = 0
C
      IF(ILEN.EQ.0) RETURN
C
C---- extract numbers
      N = 0
      K = 1
      DO 10 IPASS=1, ILEN
C------ find next space (pretend there's one after the end of the string)
        KSPACE = INDEX(INPUT(K:ILEN),' ') + K - 1
        IF(KSPACE.EQ.K-1) KSPACE = ILEN + 1
C
        IF(KSPACE.EQ.K) THEN
C------- just skip this space
         K = K+1
         GO TO 9
        ENDIF
C
C------ also find next comma
        KCOMMA = INDEX(INPUT(K:ILEN),',') + K - 1
        IF(KCOMMA.EQ.K-1) KCOMMA = ILEN + 1
C
C------ space is farther down, so we ran into something...
        N = N+1
C
C------ bug out early if no more numbers are to be read
        IF(N.GT.NINP) GO TO 11
C
C------ set ending delimiter position for this number
        KDELIM = MIN(KSPACE,KCOMMA)
C
        IF(K.EQ.KDELIM) THEN
C------- nothing but a comma... just keep looking
         K = K+1
         GO TO 9
        ENDIF
C
C------ whatever we have, it is in substring K:KEND
        KEND = KDELIM - 1
        READ(INPUT(K:KEND),*,ERR=20) RNUM(N)
        NR = N
C
C------ keep looking after delimiter
        K = KDELIM + 1
C
  9     IF(K.GE.ILEN) GO TO 11
 10   CONTINUE
C
C---- normal return
 11   CONTINUE
      ERROR = .FALSE.
      RETURN
C
C---- bzzzt !!!
 20   CONTINUE
ccc   WRITE(*,*) 'GETFLT: List-directed read error.'
      ERROR = .TRUE.
      RETURN
      END ! GETFLT



      SUBROUTINE GETNUM(INPUT,INUM,RNUM,NI,NR,NUMTYP,ERROR)
      CHARACTER*(*) INPUT, NUMTYP
      INTEGER INUM(*)
      REAL RNUM(*)
      LOGICAL ERROR
C----------------------------------------------------------------
C     Parses character string INPUT into separate arrays
C     of integer and real numbers returned in 
C     INUM(1..NI), RNUM(1..NR).
C
C     Will attempt to extract no more than NI,NR numbers
C     of each type, unless NI,NR = 0, in which case all 
C     numbers present in INPUT will be extracted.
C
C     NI,NR return how many numbers were actually extracted.
C
C     String NUMTYP indicates into which array each number went...
C
C     NUMTYP(N:N) = 'i'   N'th number in INPUT went into INUM(N)
C                   'r'   N'th number in INPUT went into RNUM(N)
C                   'n'   N'th number in INPUT was blank (just a comma)
C----------------------------------------------------------------
C
C---- number of characters to be examined
      ILEN = LEN(INPUT)
C
C---- ignore everything after a "!" character
      K = INDEX(INPUT,'!')
      IF(K.GT.0) ILEN = K-1
C
C---- set limit on numbers to be read
      NIINP = NI
      NRINP = NR
      IF(NIINP.EQ.0) NIINP = ILEN/2 + 1
      IF(NRINP.EQ.0) NRINP = ILEN/2 + 1
      NINP = MAX( NIINP , NRINP )
C
      NI = 0
      NR = 0
      NUMTYP = ' '
C
      IF(ILEN.EQ.0) RETURN
C
C---- extract numbers
      N = 0
      K = 1
      DO 10 IPASS=1, ILEN
C------ find next space (pretend there's one after the end of the string)
        KSPACE = INDEX(INPUT(K:ILEN),' ') + K - 1
        IF(KSPACE.EQ.K-1) KSPACE = ILEN + 1
C
        IF(KSPACE.EQ.K) THEN
C------- just skip this space
         K = K+1
         GO TO 9
        ENDIF
C
C------ also find next comma
        KCOMMA = INDEX(INPUT(K:ILEN),',') + K - 1
        IF(KCOMMA.EQ.K-1) KCOMMA = ILEN + 1
C
C------ space is farther down, so we ran into something...
        N = N+1
C
C------ bug out early if no more numbers are to be read
        IF(N.GT.NINP) GO TO 11
C
C------ set ending delimiter position for this number
        KDELIM = MIN(KSPACE,KCOMMA)
C
        IF(K.EQ.KDELIM) THEN
C------- nothing but a comma... just set null type indicator and keep looking
         NUMTYP(N:N) = 'n'
         K = K+1
         GO TO 9
        ENDIF
C
C------ whatever we have, it is in substring K:KEND
        KEND = KDELIM - 1
C
C------ search for floating-point number indicator in substring
        KFLOAT = MAX( INDEX(INPUT(K:KEND),'.'),
     &                INDEX(INPUT(K:KEND),'E'),
     &                INDEX(INPUT(K:KEND),'e'),
     &                INDEX(INPUT(K:KEND),'D'),
     &                INDEX(INPUT(K:KEND),'d')  ) + K - 1
C
        IF(KFLOAT.GE.K .AND. KFLOAT.LE.KEND) THEN
C------- real number... read it only if max has not been reached
         IF(N.LE.NRINP) THEN
          READ(INPUT(K:KEND),*,ERR=20) RNUM(N)
          NUMTYP(N:N) = 'r'
          NR = N
         ENDIF
        ELSE
C------- integer number...
         IF(N.LE.NIINP) THEN
          READ(INPUT(K:KEND),*,ERR=20) INUM(N)
          NUMTYP(N:N) = 'i'
          NI = N
         ENDIF
        ENDIF
C
C------ keep looking after delimiter
        K = KDELIM + 1
C
  9     IF(K.GE.ILEN) GO TO 11
 10   CONTINUE
C
C---- normal return
 11   CONTINUE
      ERROR = .FALSE.
      RETURN
C
C---- bzzzt !!!
 20   CONTINUE
ccc   WRITE(*,*) 'GETNUM: List-directed read error.'
      ERROR = .TRUE.
      RETURN
      END



      SUBROUTINE GETARG0(IARG,ARG)
C------------------------------------------------
C     Same as GETARG, but...
C
C     ...in the case of Intel Fortran, this one
C     doesn't barf if there's no Unix argument 
C      (just returns blank string instead)
C------------------------------------------------
      CHARACTER*(*) ARG
C
      NARG = IARGC()
      IF(NARG.GE.IARG) THEN
       CALL GETARG(IARG,ARG)
      ELSE
       ARG = ' '
      ENDIF
C
      RETURN
      END ! GETARG0
