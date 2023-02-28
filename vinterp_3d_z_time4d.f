C NCLFORTSTART
      SUBROUTINE DINTERP3DZ(V4D,Z,LOC,V3D,NX,NY,NZ,NT)
      IMPLICIT NONE
      INTEGER NX,NY,NZ,NT
      DOUBLE PRECISION V4D(NT,NZ,NX,NY),V3D(NT,NX,NY)
      DOUBLE PRECISION Z(NT,NZ,NX,NY)
      DOUBLE PRECISION LOC
C NCLEND
      DOUBLE PRECISION VMSG
      INTEGER I,J,KP,IP,IM,IT
      LOGICAL INTERP
      DOUBLE PRECISION HEIGHT,W1,W2

      HEIGHT = LOC
      VMSG = -999999.      
c does vertical coordinate increase or decrease with increasing k?
c set offset appropriately

      IP = 0
      IM = 1
c      IF (Z(1,1,1).GT.Z(NZ,1,1)) THEN
c          IP = 1
c          IM = 0
c      END IF
      DO IT = 1,NT
        IF (Z(IT,1,1,1).GT.Z(IT,NZ,1,1)) THEN
          IP = 1
          IM = 0
        END IF
        DO I = 1,NX
          DO J = 1,NY
C Initialize to missing.  Was initially hard-coded to -999999.
              V3D(IT,I,J) = VMSG
              INTERP = .false.
              KP = NZ

              DO WHILE ((.NOT.INTERP) .AND. (KP.GE.2))

                  IF (((Z(IT,KP-IM,I,J).LE.HEIGHT).AND. (Z(IT,KP-IP,
     +                I,J).GT.HEIGHT))) THEN
                      W2 = (HEIGHT-Z(IT,KP-IM,I,J))/
     +                     (Z(IT,KP-IP,I,J)-Z(IT,KP-IM,I,J))
                      W1 = 1.D0 - W2
                      V3D(IT,I,J) = W1*V4D(IT,KP-IM,I,J)+W2*
     +                           V4D(IT,KP-IP,I,J)
                      INTERP = .true.
                  END IF
                  KP = KP - 1

              END DO

          END DO
        END DO
      END DO

      RETURN
      END
