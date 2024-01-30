module svd_lib
  use kind_parameters
  implicit none

!  private

contains


!! ------------------------------------------------------------------------------------------------
  subroutine svd_solve(Amat,NSIZE,bvec)
!!  This routine calculates the singular value decomposition in the form A=U.W.V^T
!!  and uses it to solve the system Ax=b
!!
!!  Takes as input:
!!             Amat(nsize,nsize)   !! LHS
!!             nsize               !! Size of problem
!!             bvec                !! RHS
!!
!!     OUTPUT:
!!       MATRIX U(NSIZE,NSIZE) STORED IN Amat
!!       MATRIX W = DIAG(NSIZE) STORED AS Wvec OF PHYSICAL DIMENSION nsize
!!       MATRIX V(NSIZE,NSIZE) (NOT THE TRANSPOSE) STORED IN Vmat(nsize,nsize)
!!
!!     BASED ON NUMERICAL RECIPES SUBROUTINE SVDCMP
!!
!!     *************************************************************************

!C     PARAMETERS
!C     ==========
      real(rkind) ZERO,ONE,TWO
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      integer(ikind) NMAX
      PARAMETER(NMAX=500)

!C     ARGUMENTS
!C     =========
      real(rkind), intent(in) :: Amat(nsize,nsize)
      real(rkind) Vmat(nsize,nsize),umat(nsize,nsize)
      real(rkind) Wvec(nsize)
      integer(ikind) NSIZE


!C     ARGUMENTS
!C     =========
      real(rkind), intent(inout) :: bvec(nsize)



!C     EXTERNAL FUNCTION
!C     =================
!      real(rkind) PYTHAG
!      EXTERNAL PYTHAG


!C     LOCAL DATA
!C     ==========
      real(rkind) TMP(NMAX)
      real(rkind) RV1(NMAX)
      real(rkind) ANORM,CVAR,FVAR,GVAR,HVAR,SVAR,SSCALE
      real(rkind) XVAR,YVAR,ZVAR
      integer(ikind) IC,ITS,JC,JJ,KC,LC,NM


!!    Set umat to Amat
     umat = Amat

!C     BEGIN
!C     =====
!C     HOUSEHOLDER REDUCTION TO BIDIAGONAL FORM
      GVAR = ZERO
      SSCALE = ZERO
      ANORM = ZERO
      DO IC = 1,NSIZE
        LC = IC + 1
        RV1(IC) = SSCALE*GVAR
        GVAR = ZERO
        SVAR = ZERO
        SSCALE = ZERO
        IF(IC.LE.NSIZE)THEN
          DO KC = IC,NSIZE
            SSCALE = SSCALE + ABS(umat(KC,IC))
          ENDDO
          IF(SSCALE.NE.ZERO)THEN
            DO KC = IC,NSIZE
              umat(KC,IC) = umat(KC,IC)/SSCALE
              SVAR = SVAR + umat(KC,IC)*umat(KC,IC)
            ENDDO
            FVAR = umat(IC,IC)
            GVAR = -SIGN(SQRT(SVAR),FVAR)
            HVAR = FVAR*GVAR - SVAR
            umat(IC,IC) = FVAR - GVAR
            DO JC = LC,NSIZE
              SVAR = ZERO
              DO KC = IC,NSIZE
                SVAR = SVAR + umat(KC,IC)*umat(KC,JC)
              ENDDO
              FVAR = SVAR/HVAR
              DO KC = IC,NSIZE
                umat(KC,JC) = umat(KC,JC) + FVAR*umat(KC,IC)
              ENDDO
            ENDDO
            DO KC = IC,NSIZE
              umat(KC,IC) = SSCALE*umat(KC,IC)
            ENDDO
          ENDIF
        ENDIF
        Wvec(IC) = SSCALE*GVAR
        GVAR = ZERO
        SVAR = ZERO
        SSCALE = ZERO
        IF(IC.LE.NSIZE)THEN
          DO KC = LC,NSIZE
            SSCALE = SSCALE + ABS(umat(IC,KC))
          ENDDO
          IF(SSCALE.NE.ZERO)THEN
            DO KC = LC,NSIZE
              umat(IC,KC) = umat(IC,KC)/SSCALE
              SVAR = SVAR + umat(IC,KC)*umat(IC,KC)
            ENDDO
            FVAR = umat(IC,LC)
            GVAR = -SIGN(SQRT(SVAR),FVAR)
            HVAR = FVAR*GVAR - SVAR
            umat(IC,LC) = FVAR - GVAR
            DO KC = LC,NSIZE
              RV1(KC)=umat(IC,KC)/HVAR
            ENDDO
            DO JC = LC,NSIZE
              SVAR = ZERO
              DO KC = LC,NSIZE
                SVAR = SVAR + umat(JC,KC)*umat(IC,KC)
              ENDDO
              DO KC = LC,NSIZE
                umat(JC,KC) = umat(JC,KC) + SVAR*RV1(KC)
              ENDDO
            ENDDO
            DO KC = LC,NSIZE
              umat(IC,KC)=SSCALE*umat(IC,KC)
            ENDDO
          ENDIF
        ENDIF
        ANORM = MAX(ANORM,(ABS(Wvec(IC))+ABS(RV1(IC))))
      ENDDO

!C     ACCUMULATION OF RIGHT-HAND TRANSFORMATIONS
      DO IC = NSIZE,1,-1
        IF(IC.LT.NSIZE)THEN
          IF(GVAR.NE.ZERO)THEN
!C           DOUBLE DIVISION TO AVOID POSSIBLE UNDERFLOW
            DO JC = LC,NSIZE
              Vmat(JC,IC) = (umat(IC,JC)/umat(IC,LC))/GVAR
            ENDDO
            DO JC = LC,NSIZE
              SVAR = ZERO
              DO KC = LC,NSIZE
                SVAR = SVAR + umat(IC,KC)*Vmat(KC,JC)
              ENDDO
              DO KC = LC,NSIZE
                Vmat(KC,JC) = Vmat(KC,JC) + SVAR*Vmat(KC,IC)
              ENDDO
            ENDDO
          ENDIF
          DO JC = LC,NSIZE
            Vmat(IC,JC) = ZERO
            Vmat(JC,IC) = ZERO
          ENDDO
        ENDIF
        Vmat(IC,IC) = ONE
        GVAR = RV1(IC)
        LC = IC
      ENDDO

!C     ACCUMULATION OF LEFT-HAND TRANSFORMATIONS
      DO IC = NSIZE,1,-1
        LC = IC + 1
        GVAR = Wvec(IC)
        DO JC = LC,NSIZE
          umat(IC,JC) = ZERO
        ENDDO
        IF(GVAR.NE.ZERO)THEN
          GVAR = ONE/GVAR
          DO JC=LC,NSIZE
            SVAR = ZERO
            DO KC = LC,NSIZE
              SVAR = SVAR + umat(KC,IC)*umat(KC,JC)
            ENDDO
            FVAR = (SVAR/umat(IC,IC))*GVAR
            DO KC = IC,NSIZE
              umat(KC,JC) = umat(KC,JC) + FVAR*umat(KC,IC)
            ENDDO
          ENDDO
          DO JC = IC,NSIZE
            umat(JC,IC) = umat(JC,IC)*GVAR
          ENDDO
        ELSE
          DO JC= IC,NSIZE
            umat(JC,IC) = ZERO
          ENDDO
        ENDIF
        umat(IC,IC) = umat(IC,IC) + ONE
      ENDDO

!C     DIAGONALIZATION OF THE BIDIAGONAL FORM
!C     LOOP OVER SINGULAR VALUES
      NM = 1
      DO KC = NSIZE,1,-1
!C       LOOP OVER ALLOWED ITERATIONS
        DO ITS = 1,30
          DO LC = KC,1,-1
!C           TEST FOR SPLITTING
            NM = LC-1
!C           NOTE THAT RV1(1) IS ALWAYS ZERO
            IF((ABS(RV1(LC))+ANORM).EQ.ANORM) GOTO 2000
            IF((ABS(Wvec(NM))+ANORM).EQ.ANORM) GOTO 1000
          ENDDO

!C         CANCELLATION OF RV1(LC) IF LC > 1
1000      CVAR = ZERO
          SVAR = ONE
          DO IC = LC,KC
            FVAR = SVAR*RV1(IC)
            RV1(IC) = CVAR*RV1(IC)
            IF((ABS(FVAR)+ANORM).EQ.ANORM) GOTO 2000
            GVAR = Wvec(IC)
            HVAR = PYTHAG(FVAR,GVAR)
            Wvec(IC) = HVAR
            HVAR = ONE/HVAR
            CVAR =  (GVAR*HVAR)
            SVAR = -(FVAR*HVAR)
            DO JC = 1,NSIZE
              YVAR = umat(JC,NM)
              ZVAR = umat(JC,IC)
              umat(JC,NM) =  (YVAR*CVAR)+(ZVAR*SVAR)
              umat(JC,IC) = -(YVAR*SVAR)+(ZVAR*CVAR)
            ENDDO
          ENDDO

2000      ZVAR = Wvec(KC)
          IF(LC.EQ.KC)THEN
!C           CONVERGENCE
            IF(ZVAR.LT.ZERO)THEN
!C             SINGULAR VALUE IS MADE NONNEGATIVE
              Wvec(KC) = -ZVAR
              DO JC = 1,NSIZE
                Vmat(JC,KC) = -Vmat(JC,KC)
              ENDDO
            ENDIF
            GOTO 3000
          ENDIF
          IF(ITS.EQ.30)THEN
            WRITE(6,*)'SVDCMP: no convergence'
            WRITE(6,'(4I7)')NSIZE
            DO IC = 1,NSIZE
              WRITE(6,'(19(1PE15.7))')(umat(IC,JC),JC=1,NSIZE)
            ENDDO
            WRITE(6,*)
            DO IC = 1,NSIZE
              WRITE(6,'(19(1PE15.7))')Wvec(IC)
            ENDDO
            WRITE(6,*)
            DO IC = 1,NSIZE
              WRITE(6,'(19(1PE15.7))')(Vmat(IC,JC),JC=1,NSIZE)
            ENDDO
            WRITE(6,*)
            STOP
          ENDIF

!C         SHIFT FROM BOTTOM 2-BY-2 MINOR
          XVAR = Wvec(LC)
          NM = KC-1
          YVAR = Wvec(NM)
          GVAR = RV1(NM)
          HVAR = RV1(KC)
          FVAR = ((YVAR-ZVAR)*(YVAR+ZVAR) +  (GVAR-HVAR)*(GVAR+HVAR))/(TWO*HVAR*YVAR)
          GVAR = PYTHAG(FVAR,ONE)
          FVAR = ((XVAR-ZVAR)*(XVAR+ZVAR) +  HVAR*((YVAR/(FVAR+SIGN(GVAR,FVAR)))-HVAR))/XVAR
!C         NEXT QR TRANSFORMATION
          CVAR = ONE
          SVAR = ONE
          DO JC = LC,NM
            IC = JC+1
            GVAR = RV1(IC)
            YVAR = Wvec(IC)
            HVAR = SVAR*GVAR
            GVAR = CVAR*GVAR
            ZVAR = PYTHAG(FVAR,HVAR)
            RV1(JC) = ZVAR
            CVAR = FVAR/ZVAR
            SVAR = HVAR/ZVAR
            FVAR =  (XVAR*CVAR)+(GVAR*SVAR)
            GVAR = -(XVAR*SVAR)+(GVAR*CVAR)
            HVAR = YVAR*SVAR
            YVAR = YVAR*CVAR
            DO JJ = 1,NSIZE
              XVAR = Vmat(JJ,JC)
              ZVAR = Vmat(JJ,IC)
              Vmat(JJ,JC) =  (XVAR*CVAR)+(ZVAR*SVAR)
              Vmat(JJ,IC) = -(XVAR*SVAR)+(ZVAR*CVAR)
            ENDDO
            ZVAR = PYTHAG(FVAR,HVAR)
            Wvec(JC) = ZVAR
!C           ROTATION CAN BE ARBITRARY IF ZVAR = 0
            IF(ZVAR.NE.ZERO)THEN
              ZVAR = ONE/ZVAR
              CVAR = FVAR*ZVAR
              SVAR = HVAR*ZVAR
            ENDIF
            FVAR =  (CVAR*GVAR)+(SVAR*YVAR)
            XVAR = -(SVAR*GVAR)+(CVAR*YVAR)
            DO JJ = 1,NSIZE
              YVAR = umat(JJ,JC)
              ZVAR = umat(JJ,IC)
              umat(JJ,JC) =  (YVAR*CVAR)+(ZVAR*SVAR)
              umat(JJ,IC) = -(YVAR*SVAR)+(ZVAR*CVAR)
            ENDDO
          ENDDO
          RV1(LC) = ZERO
          RV1(KC) = FVAR
          Wvec(KC) = XVAR
        ENDDO

3000    CONTINUE

      ENDDO


!     *************************************************************************
      !! Now use the SVD to solve the linear system

!     BEGIN
!     =====
!     CALCULATE U^T.B
      DO JC = 1,NSIZE
        SVAR = ZERO
!       NONZERO RESULT ONLY IF W(J) IS NONZERO
        IF(Wvec(JC).NE.ZERO)THEN
          DO IC = 1,NSIZE
            SVAR = SVAR + umat(IC,JC)*bvec(IC)
          ENDDO
!         THIS IS THE DIVIDE BY W(J) .
          SVAR = SVAR/Wvec(JC)
        ENDIF
        TMP(JC) = SVAR
      ENDDO

!     MATRIX MULTIPLY BY V TO GET THE ANSWER
      DO JC = 1,NSIZE
        SVAR = ZERO
        DO JJ = 1,NSIZE
          SVAR = SVAR + Vmat(JC,JJ)*TMP(JJ)
        ENDDO
        bvec(JC) = SVAR
      ENDDO

      RETURN
   end subroutine svd_solve

!! ------------------------------------------------------------------------------------------------
   function pythag(AVALUE,BVALUE)
   !! This function calculates sqrt(A^2 + B^2) without destructive under/over -flow
      real(rkind) ZERO,ONE
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)

      real(rkind) PYTHAG
      real(rkind) AVALUE,BVALUE
      real(rkind) ABSA,ABSB   !! Local

      ABSA = ABS(AVALUE)
      ABSB = ABS(BVALUE)
      IF(ABSA.GT.ABSB)THEN
        PYTHAG = ABSA*SQRT(ONE+(ABSB/ABSA)**2)
      ELSE
        IF(ABSB.EQ.ZERO)THEN
          PYTHAG = ZERO
        ELSE
          PYTHAG = ABSB*SQRT(ONE+(ABSA/ABSB)**2)
        ENDIF
      ENDIF


      return
   end function pythag


!! ------------------------------------------------------------------------------------------------
end module svd_lib
