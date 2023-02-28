!! USAGE:    CALHEL(UST,VST,HELI)
!!
!!   INPUT ARGUMENT LIST:
!!     DPTH      - DEPTH IN METERS OVER WHICH HELICITY SHOULD BE COMPUTED;
!!   ALLOWS ONE TO DISTINGUISH 0-3 KM, 0-1 KM, and 0-500m  VALUES
!!     U : 3D ARRAY OF U COMPONENT OF WIND
!!     V : 3D ARRAY OF V COMPONENT OF WIND
!!     Z : 3D ARRAY of AGL HEIGHTS
!!   OUTPUT ARGUMENT LIST: 
!!     HELI     - STORM-RELATIVE HELICITY (M**2/S**2)
!!
!!     LANGUAGE: FORTRAN 90
      subroutine CALHELPLEV(U,V,U10,V10,Z2,DEPTH,NDX,IX,JY,KZ,HELI)
      implicit none
!     
      real,parameter :: P150=15000.0,P300=30000.0,S15=15.0
      real,parameter :: D3000=3000.0,PI6=0.5235987756,PI9=0.34906585
      real,parameter :: D5500=5500.0,D6000=6000.0,D7000=7000.0
      real,parameter :: D500=500.0
      real,parameter :: D1000=1000.0
      real,parameter :: D1500=1500.0
!     
!     DECLARE VARIABLES
!     
      real,dimension(NDX),intent(in)        :: DEPTH
      real,dimension(IX,JY,KZ),intent(in)   :: U,V,Z2
      real,dimension(IX,JY),intent(in)   :: U10,V10
      integer,intent(in):: NDX,IX,JY,KZ
      real,dimension(IX,JY,NDX),intent(out) ::HELI
      real,dimension(IX,JY) :: UST6
      real,dimension(IX,JY) :: VST6
      real,dimension(IX,JY) :: UST5
      real,dimension(IX,JY) :: UST
      real,dimension(IX,JY) :: VST
      real,dimension(IX,JY) :: VST5
      real,dimension(IX,JY) :: UST1
      real,dimension(IX,JY) :: VST1
      real,dimension(IX,JY) :: USHR6
      real,dimension(IX,JY) :: VSHR6
      real,dimension(IX,JY) :: U1
      real,dimension(IX,JY) :: V1
      real,dimension(IX,JY) :: U2
      real,dimension(IX,JY) :: V2
      real,dimension(IX,JY) :: HGT1
      real,dimension(IX,JY) :: HGT2
      real,dimension(IX,JY) :: UMEAN
      real,dimension(IX,JY) :: VMEAN
      integer, dimension(IX,JY) :: COUNT6, COUNT5,COUNT1,L1, L2
      integer :: I,J,L,N
      real :: DZABV,UMEAN5,VMEAN5,UMEAN1,VMEAN1,UMEAN6,VMEAN6
      real :: DENOM,DZ,DZ1,DZ2,DU1,DU2,DV1,DV2

      
!$omp  parallel do private(i,j)
      DO J=1,JY
        DO I=1,IX
          UST(I,J)    = 0.0
          VST(I,J)    = 0.0
          HELI(I,J,:)   = 0.0
          UST1(I,J)   = U10(I,J)!0.0
          VST1(I,J)   = V10(I,J)!0.0
          UST5(I,J)   = 0.0
          VST5(I,J)   = 0.0
          UST6(I,J)   = U10(I,J)!0.0
          VST6(I,J)   = V10(I,J)!0.0
          COUNT6(I,J) = 1
          COUNT5(I,J) = 0
          COUNT1(I,J) = 1
          USHR6(I,J)  = 0.0
          VSHR6(I,J)  = 0.0
          U1(I,J)     = 0.0
          U2(I,J)     = 0.0
          V1(I,J)     = 0.0
          V2(I,J)     = 0.0
          UMEAN(I,J)  = 0.0
          VMEAN(I,J)  = 0.0
          HGT1(I,J)   = 0.0
          HGT2(I,J)   = 0.0
          L1(I,J)     = 0
          L2(I,J)     = 0

        ENDDO
      ENDDO
!     LOOP OVER HORIZONTAL GRID.
      DO L = 1,KZ
        DO J=1,JY
          DO I=1,IX
            DZABV = Z2(I,J,L)
            !0-6km mean wind
            IF (DZABV <= D6000 .AND. DZABV > 10.0) THEN
               UST6 (I,J)   = UST6(I,J) + U(I,J,L)
               VST6 (I,J)   = VST6(I,J) + V(I,J,L)
               COUNT6 (I,J) = COUNT6(I,J) + 1
            ENDIF
!            UST6 (I,J)   = UST6(I,J) + U10(I,J)
!            VST6 (I,J)   = VST6(I,J) + V10(I,J)
!            COUNT6 (I,J) = COUNT6(I,J) + 1
            IF (DZABV < D6000 .AND. DZABV >= D5500) THEN
               UST5(I,J)   = UST5(I,J) + U(I,J,L)
               VST5(I,J)   = VST5(I,J) + V(I,J,L)
               COUNT5(I,J) = COUNT5(I,J) + 1
            ENDIF
            !0-500 m wind
            IF (DZABV < D500 .AND. DZABV > 10.0) THEN
               UST1(I,J)   = UST1(I,J) + U(I,J,L)
               VST1(I,J)   = VST1(I,J) + V(I,J,L)
               COUNT1(I,J) = COUNT1(I,J) + 1
            ENDIF
!            UST1 (I,J)   = UST1(I,J) + U10(I,J)
!            VST1 (I,J)   = VST1(I,J) + V10(I,J)
!            COUNT1 (I,J) = COUNT1(I,J) + 1

            IF (DZABV >= D1000 .AND. DZABV <= D1500) THEN 
               U2(I,J)   = U(I,J,L)
               V2(I,J)   = V(I,J,L)
               HGT2(I,J) = DZABV
               L2(I,J)   = L
            ENDIF

            IF (DZABV >= D500 .AND. DZABV < D1000 .AND.               &
               L1(I,J) <= L2(I,J)) THEN
               U1(I,J)   = U(I,J,L)
               V1(I,J)   = V(I,J,L)
               HGT1(I,J) = DZABV
               L1(I,J)   = L
            ENDIF
          ENDDO
        ENDDO
      ENDDO

      DO J=1,JY
        DO I=1,IX
          IF (COUNT6(I,J) > 0 .AND. COUNT1(I,J) > 0 .AND.    &
             COUNT5(I,J) >0) THEN
            UMEAN5 = UST5(I,J) / COUNT5(I,J)
            VMEAN5 = VST5(I,J) / COUNT5(I,J)
            UMEAN1 = UST1(I,J) / COUNT1(I,J)
            VMEAN1 = VST1(I,J) / COUNT1(I,J)
            UMEAN6 = UST6(I,J) / COUNT6(I,J)
            VMEAN6 = VST6(I,J) / COUNT6(I,J)


!           COMPUTE STORM MOTION VECTOR
!           IT IS DEFINED AS 7.5 M/S TO THE RIGHT OF THE 0-6 KM MEAN
!           WIND CONSTRAINED ALONG A LINE WHICH IS BOTH PERPENDICULAR
!           TO THE 0-6 KM MEAN VERTICAL WIND SHEAR VECTOR AND PASSES
!           THROUGH THE 0-6 KM MEAN WIND.  THE WIND SHEAR VECTOR IS
!           SET AS THE DIFFERENCE BETWEEN THE 5.5-6 KM WIND (THE HEAD
!           OF THE SHEAR VECTOR) AND THE 0-0.5 KM WIND (THE TAIL).
!           THIS IS FOR THE RIGHT-MOVING CASE;  WE IGNORE THE LEFT MOVER.

            USHR6(I,J) = UMEAN5 - UMEAN1
            VSHR6(I,J) = VMEAN5 - VMEAN1

            DENOM = USHR6(I,J)*USHR6(I,J)+VSHR6(I,J)*VSHR6(I,J)
            IF (DENOM /= 0.0) THEN
              UST(I,J) = UMEAN6 + (7.5*VSHR6(I,J)/SQRT(DENOM))
              VST(I,J) = VMEAN6 - (7.5*USHR6(I,J)/SQRT(DENOM))
            ELSE
              UST(I,J) = 0
              VST(I,J) = 0
            ENDIF
          ELSE
            UST(I,J) = 0.0
            VST(I,J) = 0.0
            USHR6(I,J) = 0.0
            VSHR6(I,J) = 0.0
          ENDIF

          IF(L1(I,J) > 0 .AND. L2(I,J) > 0) THEN
            UMEAN(I,J) = U1(I,J) + (D1000 - HGT1(I,J))*(U2(I,J)- &
                                    U1(I,J))/(HGT2(I,J) - HGT1(I,J))
            VMEAN(I,J) = V1(I,J) + (D1000 - HGT1(I,J))*(V2(I,J)- &
                                    V1(I,J))/(HGT2(I,J) - HGT1(I,J))
          ELSE IF(L1(I,J) > 0 .AND. L2(I,J) == 0) THEN
            UMEAN(I,J) = U1(I,J)
            VMEAN(I,J) = V1(I,J)
          ELSE IF(L1(I,J) == 0 .AND. L2(I,J) > 0) THEN
            UMEAN(I,J) = U2(I,J)
            VMEAN(I,J) = V2(I,J)
          ELSE
            UMEAN(I,J) = 0.0
            VMEAN(I,J) = 0.0
          ENDIF

        ENDDO
      ENDDO
!
!       COMPUTE STORM-RELATIVE HELICITY
!
      DO N=1,NDX ! for different helicity depths
       DO L = 2,KZ-1
          DO J=1,JY
            DO I=1,IX
              DZABV=Z2(I,J,L)
              IF(DZABV <= DEPTH(N) .AND. DZABV > 10.0) THEN
                IF (L .EQ. 2) THEN
                    !if L=2 (first loop), compute diffs against 10m U/V
                    DZ = Z2(I,J,L)-10.0
                    DZ1 = Z2(I,J,L+1)-Z2(I,J,L)
                    DZ2 = Z2(I,J,L)-10.0
                    DU1 = U(I,J,L+1)-U(I,J,L)
                    DU2 = U(I,J,L)-U10(I,J)
                    DV1 = V(I,J,L+1)-V(I,J,L)
                    DV2 = V(I,J,L)-V10(I,J)
                
                ELSE
                    DZ = Z2(I,J,L)-Z2(I,J,L-1)
                    DZ1 = Z2(I,J,L+1)-Z2(I,J,L)
                    DZ2 = Z2(I,J,L)-Z2(I,J,L-1)
                    DU1 = U(I,J,L+1)-U(I,J,L)
                    DU2 = U(I,J,L)-U(I,J,L-1)
                    DV1 = V(I,J,L+1)-V(I,J,L)
                    DV2 = V(I,J,L)-V(I,J,L-1)
                END IF


                HELI(I,J,N) = ((V(I,J,L)-VST(I,J))* &
                               (DZ2*(DU1/DZ1)+DZ1*(DU2/DZ2)) &
                            -  (U(I,J,L)-UST(I,J))* &
                               (DZ2*(DV1/DZ1)+DZ1*(DV2/DZ2))) &
                               *DZ/(DZ1+DZ2)+HELI(I,J,N)
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      END DO  
      end subroutine CALHELPLEV
   
