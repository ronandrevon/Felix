!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! felixrefine
!
! Richard Beanland, Keith Evans, Rudolf A Roemer and Alexander Hubert
!
! (C) 2013/14, all right reserved
!
! Version: :VERSION:
! Date:    :DATE:
! Time:    :TIME:
! Status:  :RLSTATUS:
! Build:   :BUILD:
! Author:  :AUTHOR:
! 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!  This file is part of felixrefine.
!
!  felixrefine is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
!  
!  felixrefine is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!  
!  You should have received a copy of the GNU General Public License
!  along with felixrefine.  If not, see <http://www.gnu.org/licenses/>.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! $Id: Felixrefine.f90,v 1.89 2014/04/28 12:26:19 phslaz Exp $
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PROGRAM Felixrefine
 
  USE MyNumbers
  
  USE CConst; USE IConst; USE RConst
  USE IPara; USE RPara; USE SPara; USE CPara
  USE BlochPara

  USE IChannels

  USE MPI
  USE MyMPI

  !--------------------------------------------------------------------
  ! local variable definitions
  IMPLICIT NONE

  INTEGER(IKIND) :: IHours,IMinutes,ISeconds,IErr,IMilliSeconds,IIterationFLAG,&
       ind,jnd,knd,ICalls,Iter,ICutOff,IHOLZgPoolMag,IBSMaxLocGVecAmp,&
	   ILaueLevel,INumTotalReflections,ITotalLaueZoneLevel,INhkl,IExitFLAG,&
	   INumInitReflections,IZerothLaueZoneLevel,INumFinalReflections
  INTEGER(IKIND) :: IStartTime,ICurrentTime,IRate
  INTEGER(IKIND), DIMENSION(:),ALLOCATABLE :: IOriginGVecIdentifier
  REAL(RKIND) :: StartTime, CurrentTime, Duration, TotalDurationEstimate,&
       RFigureOfMerit,RHOLZAcceptanceAngle,RLaueZoneGz
  REAL(RKIND), DIMENSION(:,:), ALLOCATABLE :: RSimplexVariable
  REAL(RKIND),DIMENSION(:),ALLOCATABLE :: RSimplexFoM,RIndependentVariable
  REAL(RKIND), DIMENSION(:,:), ALLOCATABLE :: RgDummyVecMat,RgPoolMagLaue
  REAL(RKIND) :: RBCASTREAL,RStandardDeviation,RMean,RGzUnitVec,RMinLaueZoneValue,&
       RMaxLaueZoneValue,RMaxAcceptanceGVecMag,RLaueZoneElectronWaveVectorMag
  REAL(RKIND) :: RdeltaUg,Rtol,RpointA,RpointB,RpointC,RfitA,RfitB,RfitC,RbestFit
  LOGICAL :: LInitialSimulationFLAG = .TRUE.
  CHARACTER*40 :: my_rank_string
  CHARACTER*20 :: Sind
  CHARACTER*200 :: SPrintString

  !-------------------------------------------------------------------
  ! constants
  CALL Init_Numbers
  
  !-------------------------------------------------------------------
  ! set the error value to zero, will change upon error
  IErr=0

  !--------------------------------------------------------------------
  ! MPI initialization
  CALL MPI_Init(IErr)  
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixrefine(", my_rank, ") error in MPI_Init()"
     GOTO 9999
  END IF

  ! Get the rank of the current process
  CALL MPI_Comm_rank(MPI_COMM_WORLD,my_rank,IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixrefine(", my_rank, ") error in MPI_Comm_rank()"
     GOTO 9999
  END IF

  ! Get the size of the current communicator
  CALL MPI_Comm_size(MPI_COMM_WORLD,p,IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixrefine(", my_rank, ") error in MPI_Comm_size()"
     GOTO 9999
  END IF

  !--------------------------------------------------------------------
  ! startup
  IF((IWriteFLAG.GE.0.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"--------------------------------------------------------------"
     PRINT*,"felixrefine: ", RStr
     PRINT*,"             ", DStr
     PRINT*,"             ", AStr
     PRINT*,"    on rank= ", my_rank, " of ", p, " in total."
     PRINT*,"--------------------------------------------------------------"
  END IF
  ISoftwareMode =2 ! felixrefinemode

  !--------------------------------------------------------------------
  ! timing startup
  CALL SYSTEM_CLOCK(count_rate=IRate)
  CALL SYSTEM_CLOCK(IStarttime)

  !--------------------------------------------------------------------
  ! INPUT section 
   CALL ReadInput (IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"felixrefine(",my_rank,") error in ReadInput"
     GOTO 9999
  END IF  

  ALLOCATE(RImageExpi(2*IPixelCount,2*IPixelCount,INoOfLacbedPatterns),STAT=IErr)  
  IF( IErr.NE.0 ) THEN
     PRINT*,"felixrefine(",my_rank,")error allocating RImageExpi"
     GOTO 9999
  END IF
  
  CALL ReadExperimentalImages(IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"felixrefine(",my_rank,") error in ReadExperimentalImages"
     GOTO 9999
  END IF

  !--------------------------------------------------------------------
  ! Initial simulation and variable setup
  CALL MicroscopySettings( IErr )
  IF( IErr.NE.0 ) THEN
    PRINT*,"felixrefine(",my_rank,") error in MicroscopySettings"
    GOTO 9999
  ENDIF  
 
  CALL ReciprocalLattice(IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"felixrefine(",my_rank,")error in ReciprocalLattice"
     GOTO 9999
  ENDIF

  ALLOCATE(RAtomCoordinate(ITotalAtoms,THREEDIM),STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"felixrefine(",my_rank,")error allocating RAtomCoordinate"
     GOTO 9999
  ENDIF
  CALL AllAtomPositions(IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"felixrefine(",my_rank,")error in AllAtomPositions"
     GOTO 9999
  ENDIF

!deallocations
DEALLOCATE(RFullPartialOccupancy,SMNP,MNP,RFullAtomicFracCoordVec,SFullAtomicNameVec, &
RFullIsotropicDebyeWallerFactor,IFullAtomicNumber,IFullAnisotropicDWFTensor)

! set up reflection pool
!-----------------------------------------
!zz from diffractionpatterninitialisation/reflectiondetermination
  RHOLZAcceptanceAngle=TWODEG2RADIAN!RB seems way too low?
  IHKLMAXValue = 5!RB starting value, increments in loop below
! Count the reflections that make up the pool of g-vectors
!Note the application of acceptance angle is incorrect since it uses hkl;
!it should use the reciprocal lattice vectors as calculated in RgPool
  CALL HKLCount(IHKLMAXValue,RZDirC,INhkl,RHOLZAcceptanceAngle,IErr)
  DO WHILE (INhkl.LT.IMinReflectionPool) 
     IHKLMAXValue = IHKLMAXValue*2
     CALL HKLCount(IHKLMAXValue,RZDirC,INhkl,RHOLZAcceptanceAngle,IErr)
  END DO
! Fill the list of reflections Rhkl as indices h,k,l
  ALLOCATE(Rhkl(INhkl,THREEDIM),STAT=IErr)
  CALL HKLMake(IHKLMAXValue,RZDirC,RHOLZAcceptanceAngle,IErr)
  
  ! sort them in descending order of magnitude
  ! may result in an error when the reflection pool does not reach the highest hkl of the experimental data? 
  CALL SortHKL(Rhkl,INhkl,IErr) 
  IF( IErr.NE.0 ) THEN
     PRINT*,"felixrefine(",my_rank,")error in SortHKL"
     GOTO 9999
  END IF

  !Assign numbers to the different reflections in IOutputReflections
  CALL SpecificReflectionDetermination (IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"felixrefine(",my_rank,")error in SpecificReflectionDetermination"
     GOTO 9999
  END IF

!allocations-----------------------------------  
  ALLOCATE(RgPoolT(INhkl,THREEDIM),STAT=IErr)
  IF(IErr.NE.0) THEN
     PRINT*,"felixrefine(",my_rank,")error allocating RgPoolT"
     GOTO 9999
  END IF
  ALLOCATE(RgDummyVecMat(INhkl,THREEDIM),STAT=IErr)
  IF(IErr.NE.0) THEN
     PRINT*,"felixrefine(",my_rank,")error allocating RgDummyVecMat"
     GOTO 9999
  END IF
 
!Calculate the g vector list RgPoolT in reciprocal angstrom units (in the microscope reference frame?)
  ICutOff = 1
  DO ind=1,INhkl
     WRITE(Sind,'(I10.1)')ind
     DO jnd=1,THREEDIM
        RgPoolT(ind,jnd)= &
             Rhkl(ind,1)*RarVecM(jnd) + &
             Rhkl(ind,2)*RbrVecM(jnd) + &
             Rhkl(ind,3)*RcrVecM(jnd)
        RgDummyVecMat(ind,jnd)=RgPoolT(ind,jnd)
     ENDDO
	 !If a g-vector has a non-zero z-component it is not in the ZOLZ
     IF((RgPoolT(ind,3).GT.TINY.OR.RgPoolT(ind,3).LT.-TINY).AND.ICutOff.NE.0) THEN
        RGzUnitVec=ABS(RgPoolT(ind,3))
        ICutOff=0
     END IF
  ENDDO
  !No higher order Laue Zones with ZOLZFlag switched on
  IF(ICutOff.EQ.1) THEN
     RGzUnitVec=ZERO
  END IF
  
  !sort into Laue Zones 
  WHERE(RgDummyVecMat(:,3).GT.TINY.OR.RgDummyVecMat(:,3).LT.-TINY)
     RgDummyVecMat(:,3)=RgDummyVecMat(:,3)/RGzUnitVec!possible divide by zero from line 237?
  END WHERE
  !min&max Laue Zones 
  RMaxLaueZoneValue=MAXVAL(RgDummyVecMat(:,3),DIM=1)
  RMinLaueZoneValue=MINVAL(RgDummyVecMat(:,3),DIM=1)
  ITotalLaueZoneLevel=NINT(RMaxLaueZoneValue+ABS(RMinLaueZoneValue)+1,IKIND)
  !doing what, here? More sorting inlo Laue zones? 
  !deallocation
  DEALLOCATE(RgDummyVecMat)
  ALLOCATE(RgPoolMagLaue(INhkl,ITotalLaueZoneLevel),STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"felixrefine(",my_rank,")error allocating RgPoolMagLaue"
     GOTO 9999
  END IF
  IF(RAcceptanceAngle.NE.ZERO.AND.IZOLZFLAG.EQ.0) THEN
     INumtotalReflections=0
     DO ind=1,ITotalLaueZoneLevel
        ILaueLevel=ind-IZerothLaueZoneLevel
        RLaueZoneGz=RGzUnitVec*ILaueLevel
        DO jnd=1,INhkl
           IF(RgPoolT(jnd,3).GE.(RLaueZoneGz-TINY).AND. &
                RgPoolT(jnd,3).LE.(RLaueZoneGz+TINY)) THEN
              RgPoolMagLaue(jnd,ind)=SQRT((RgPoolT(jnd,1)**2)+(RgPoolT(jnd,2)**2))              
           ELSE
              RgPoolMagLaue(jnd,ind)=NEGHUGE
           END IF
        END DO
        INumInitReflections=COUNT(RgPoolMagLaue(:,ind).NE.NEGHUGE)
        RLaueZoneElectronWaveVectorMag=RElectronWaveVectorMagnitude-ABS(RLaueZoneGz)           
        RMaxAcceptanceGVecMag=(RLaueZoneElectronWaveVectorMag*TAN(RAcceptanceAngle*DEG2RADIAN))
        WHERE(ABS(RgPoolMagLaue(:,ind)).GT.RMaxAcceptanceGVecMag)
           RgPoolMagLaue(:,ind)=NEGHUGE
        END WHERE
        INumFinalReflections=COUNT(RgPoolMagLaue(:,ind).NE.NEGHUGE)
        INumTotalReflections=INumTotalReflections+INumInitReflections
     END DO
     knd=0
     DO ind=1,INhkl
        IF(SUM(RgPoolMagLaue(ind,:))/REAL(ITotalLaueZoneLevel,RKIND).GT.NEGHUGE) THEN
           knd=knd+1
        END IF
     END DO
     IHOLZgPoolMag=knd

     ALLOCATE(IOriginGVecIdentifier(IHOLZgPoolMag),STAT=IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"felixrefine(",my_rank,")error allocating IOriginGVecIdentifier"
        GOTO 9999
     END IF

     IOriginGVecIdentifier=0
     knd=1
     DO ind=1,INhkl
        IF((SUM(RgPoolMagLaue(ind,:))/REAL(ITotalLaueZoneLevel,RKIND)).GT.NEGHUGE) THEN
           IOriginGVecIdentifier(knd)=ind
           knd=knd+1
        END IF
     END DO
  END IF

  !calculate g vector magnitudes for the reflection pool RgPoolMag in reciprocal Angstrom units
  ALLOCATE(RgPoolMag(INhkl),STAT=IErr)
  IF(IErr.NE.0) THEN
     PRINT*,"felixrefine(",my_rank,")error allocating RgPoolMag"
     GOTO 9999
  END IF
  DO ind=1,INhkl
     RgPoolMag(ind)= SQRT(DOT_PRODUCT(RgPoolT(ind,:),RgPoolT(ind,:)))
  ENDDO

  !some other basic numbers
  ALLOCATE(RgVecVec(INhkl),STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"felixrefine(",my_rank,") error allocating RgVecVec"
     GOTO 9999
  END IF
  CALL DiffractionPatternCalculation(IErr)
    IF( IErr.NE.0 ) THEN
     PRINT*,"felixrefine(",my_rank,")error in DiffractionPatternCalculation"
     GOTO 9999
  END IF

  !acceptance angle
  IF(RAcceptanceAngle.NE.ZERO.AND.IZOLZFLAG.EQ.1) THEN
     RMaxAcceptanceGVecMag=(RElectronWaveVectorMagnitude*TAN(RAcceptanceAngle*DEG2RADIAN))
     IF(RgPoolMag(IMinReflectionPool).GT.RMaxAcceptanceGVecMag) THEN
        RBSMaxGVecAmp = RMaxAcceptanceGVecMag 
     ELSE
        RBSMaxGVecAmp = RgPoolMag(IMinReflectionPool)
     END IF
  ELSEIF(RAcceptanceAngle.NE.ZERO.AND.IZOLZFLAG.EQ.0) THEN
     IBSMaxLocGVecAmp=MAXVAL(IOriginGVecIdentifier)
     RBSMaxGVecAmp=RgPoolMag(IBSMaxLocGVecAmp)
     IF(RgPoolMag(IBSMaxLocGVecAmp).LT.RgPoolMag(IMinreflectionPool)) THEN

     ELSE
        RBSMaxGVecAmp = RgPoolMag(IMinReflectionPool)
     END IF
  ELSE
     RBSMaxGVecAmp = RgPoolMag(IMinReflectionPool)
  END IF
  
  IThicknessCount= (RFinalThickness-RInitialThickness)/RDeltaThickness + 1

  !count reflections up to cutoff magnitude
  nReflections = 0
  nStrongBeams = 0
  nWeakBeams = 0
  DO ind=1,INhkl
     IF (ABS(RgPoolMag(ind)).LE.RBSMaxGVecAmp) THEN
        nReflections = nReflections + 1
     END IF
  ENDDO

!deallocation
  DEALLOCATE(RgPoolMagLaue)!
  IF (RAcceptanceAngle.NE.ZERO.AND.IZOLZFLAG.EQ.0) THEN
    DEALLOCATE(IOriginGVecIdentifier)
    PRINT*,"felixrefine deallocating IOriginGVecIdentifier"
  END IF

  !Calculate Ug matrix--------------------------------------------------------
  ALLOCATE(CUgMatNoAbs(nReflections,nReflections),STAT=IErr)  !RB Matrix without absorption
  ALLOCATE(CUgMatPrime(nReflections,nReflections),STAT=IErr)  !RB Matrix of absorption  
  ALLOCATE(CUgMat(nReflections,nReflections),STAT=IErr)  !RB Matrix including absorption

  !RB Matrix of sums of indices - for symmetry equivalence  in the Ug matrix, only for Ug refinement
  IF(IRefineMode(1).EQ.1 .OR. IRefineMode(12).EQ.1) THEN
    ALLOCATE(RgSumMat(nReflections,nReflections),STAT=IErr)  
  END IF
  IF( IErr.NE.0 ) THEN
     PRINT*,"felixrefine(",my_rank,") error allocating CUgMat or its components"
     GOTO 9999
  END IF
  
  CALL StructureFactorSetup(IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"felixrefine(",my_rank,")error in StructureFactorSetup"
     GOTO 9999
  END IF

  IF(IRefineMode(1).EQ.1 .OR. IRefineMode(12).EQ.1) THEN !It's a Ug refinement
    DEALLOCATE(RAtomCoordinate,STAT=IErr)!Don't need this any more
    !Identify unique Ug's and count the number of independent variables INoOfVariables
	!using the Hermitian matrix CUgMatNoAbs
    !We count over INoofUgs, specified in felix.inp
    !The count excludes Ug components that are zero and U(000), the inner potential
	IUgOffset=1!choose how many Ug's to skip in the refinement, 1 is the inner potential...

    !Matrix with numbers marking equivalent Ug's
    ALLOCATE(ISymmetryRelations(nReflections,nReflections),STAT=IErr)
    IF( IErr.NE.0 ) THEN
      PRINT*,"felixrefine(",my_rank,")error allocating ISymmetryRelations"
      GOTO 9999
    ENDIF
    CALL SetupUgsToRefine(IErr)!NB IEquivalentUgKey and CUgToRefine allocated in here
    IF( IErr.NE.0 ) THEN
      PRINT*,"felixrefine(",my_rank,")error in SetupUgsToRefine"
      GOTO 9999
    END IF
	
    IF(my_rank.EQ.0) THEN
      IF ( INoOfVariables.EQ.1 ) THEN 
        PRINT*,"Only one independent variable"
	  ELSE
        WRITE(SPrintString,FMT='(I3,1X,A21)') INoOfVariables,"independent variables"
        PRINT*,TRIM(ADJUSTL(SPrintString))
      END IF
    END IF
	
    ALLOCATE(RIndependentVariable(INoOfVariables),STAT=IErr)  
    IF( IErr.NE.0 ) THEN
      PRINT*,"felixrefine(",my_rank,")error allocating RIndependentVariable"
      GOTO 9999
    END IF
    !Fill up the IndependentVariable list with CUgMatNoAbs components
    jnd=1
    DO ind = 1+IUgOffset,INoofUgs+IUgOffset !=== temp changes so real part only***
      IF ( ABS(REAL(CUgToRefine(ind),RKIND)).GE.RTolerance ) THEN
        RIndependentVariable(jnd) = REAL(CUgToRefine(ind),RKIND)
        jnd=jnd+1
	  END IF
      IF ( ABS(AIMAG(CUgToRefine(ind))).GE.RTolerance ) THEN!===
        RIndependentVariable(jnd) = AIMAG(CUgToRefine(ind))!===
        jnd=jnd+1!===
      END IF!===
    END DO
    RIndependentVariable(jnd) = RAbsorptionPercentage![[[!===RB absorption always included in structure factor refinement as last variable
	
	IF( IErr.NE.0 ) THEN
      PRINT*,"felixrefine(",my_rank,")error in deallocation CUgMatNoAbs,CUgMatPrime"
      GOTO 9999
    END IF
  END IF
 
  !--------------------------------------------------------------------
  ! Setup Simplex Variables
  !--------------------------------------------------------------------!RB restore later for other types of refinement
!  IF(IRefineMode(2).EQ.1) THEN !It's an atom coordinate refinement
!     CALL SetupAtomicVectorMovements(IErr)
!     IF( IErr.NE.0 ) THEN
!        PRINT*,"felixrefine (", my_rank, ") error in SetupAtomicVectorMovements"
!        GOTO 9999
!     END IF
!  END IF
!  INoofElementsForEachRefinementType(2)=IRefineMode(2)*IAllowedVectors
!  INoofElementsForEachRefinementType(3)=IRefineMode(3)*SIZE(IAtomicSitesToRefine)
!  INoofElementsForEachRefinementType(4)=IRefineMode(4)*SIZE(IAtomicSitesToRefine)
!  INoofElementsForEachRefinementType(5)=IRefineMode(5)*SIZE(IAtomicSitesToRefine)*6
!  INoofElementsForEachRefinementType(6)=IRefineMode(6)*3
!  INoofElementsForEachRefinementType(7)=IRefineMode(7)*3
!  INoofElementsForEachRefinementType(8)=IRefineMode(8)
!  INoofElementsForEachRefinementType(9)=IRefineMode(9)
!  INoofElementsForEachRefinementType(10)=IRefineMode(10)
!  INoofElementsForEachRefinementType(11)=IRefineMode(11)
  !Number of independent variables
!  INoOfVariables = SUM(INoofElementsForEachRefinementType)
!This has been calculated in SetupUgsToRefine
!  IF(my_rank.EQ.0) THEN
!    IF ( INoOfVariables.EQ.1 ) THEN 
!      PRINT*,"Only one independent variable"!
!	ELSE
!      WRITE(SPrintString,FMT='(I3,1X,A21))') INoOfVariables,"independent variables"
!      PRINT*,TRIM(ADJUSTL(SPrintString))
!    END IF
!  END IF


  !--------------------------------------------------------------------
  !  Assign IDs - is this needed for a Ug refinement?
!  ALLOCATE(IIterativeVariableUniqueIDs(INoOfVariables,5),STAT=IErr)
!  IF( IErr.NE.0 ) THEN
!     PRINT*,"felixrefine(",my_rank,")error allocating IIterativeVariableUniqueIDs"
!     GOTO 9999
!  ENDIF
!  IIterativeVariableUniqueIDs = 0
!  ICalls = 0
!  DO ind = 1,IRefinementVariableTypes !Loop over all possible iterative variables
!     IF(IRefineMode(ind).EQ.1) THEN
!        DO jnd = 1,INoofElementsForEachRefinementType(ind)
!           ICalls = ICalls + 1
!           IIterativeVariableUniqueIDs(ICalls,1) = ICalls
!           CALL AssignArrayLocationsToIterationVariables(ind,jnd,IIterativeVariableUniqueIDs,IErr)
!        END DO
!     END IF
!  END DO 

  !--------------------------------------------------------------------
  ! Setup Images for output
  ALLOCATE(RhklPositions(nReflections,2),STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"felixrefine(",my_rank,") error allocating RhklPositions"
     GOTO 9999
  END IF
  CALL ImageSetup( IErr )
  IF( IErr.NE.0 ) THEN
     PRINT*,"felixrefine(",my_rank,") error in ImageSetup"
     GOTO 9999
  END IF 
  !All the individual calculations go into RSimulatedPatterns later with MPI_GATHERV
  !(Note that RSimulatedPatterns is a vector with respect to pixels, not a 2D image)
  ALLOCATE(RSimulatedPatterns(INoOfLacbedPatterns,IThicknessCount,IPixelTotal),STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"felixrefine(",my_rank,") error allocating Root Reflections"
     GOTO 9999
  END IF
  RSimulatedPatterns = ZERO
  !Allocations for the pixels to be calculated by this core  
  ILocalPixelCountMin= (IPixelTotal*(my_rank)/p)+1
  ILocalPixelCountMax= (IPixelTotal*(my_rank+1)/p) 
  ALLOCATE(RIndividualReflections(INoOfLacbedPatterns,IThicknessCount,&
         (ILocalPixelCountMax-ILocalPixelCountMin)+1),STAT=IErr)
  IF( IErr.NE.0 ) THEN
    PRINT*,"felixrefine(",my_rank,")error allocating RIndividualReflections"
    GOTO 9999
  END IF
  !position of pixels calculated by this core
  ALLOCATE(IDisplacements(p),ICount(p),STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"felixrefine(",my_rank,")error allocating IDisplacements and/or ICount"
     GOTO 9999
  END IF
  DO ind = 1,p
     IDisplacements(ind) = (IPixelTotal*(ind-1)/p)*INoOfLacbedPatterns*IThicknessCount
     ICount(ind) = (((IPixelTotal*(ind)/p) - (IPixelTotal*(ind-1)/p)))*INoOfLacbedPatterns*IThicknessCount    
  END DO

  !--------------------------------------------------------------------
  ! Weighting parameter
  ALLOCATE(RWeightingCoefficients(INoOfLacbedPatterns),STAT=IErr) 
  SELECT CASE (IWeightingFLAG)
  CASE(0)!uniform weighting
     RWeightingCoefficients = ONE
  CASE(1)!smaller g's more important
     DO ind = 1,INoOfLacbedPatterns!NB untested, does RgPoolMag(ind)match output reflection (ind)???
        RWeightingCoefficients(ind) = RgPoolMag(ind)/MAXVAL(RgPoolMag)
     END DO
  CASE(2)!larger g's more important
     DO ind = 1,INoOfLacbedPatterns
        RWeightingCoefficients(ind) = MAXVAL(RgPoolMag)/RgPoolMag(ind)
     END DO
  END SELECT

  !--------------------------------------------------------------------
  ! Allocate memory for deviation parameter and bloch calc in main loop
  ALLOCATE(RDevPara(nReflections),STAT=IErr)
  ALLOCATE(IStrongBeamList(nReflections),STAT=IErr)
  ALLOCATE(IWeakBeamList(nReflections),STAT=IErr)
  ALLOCATE(CFullWaveFunctions(nReflections),STAT=IErr)
  ALLOCATE(RFullWaveIntensity(nReflections),STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"felixrefine(",my_rank,") error in allocations for Bloch calculation"
     GOTO 9999
  END IF
  !--------------------------------------------------------------------
  !baseline simulation
  RFigureofMerit=9.999!Inital value
  Iter = 0
  CALL FelixFunction(LInitialSimulationFLAG,IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"felixrefine(",my_rank,")error in FelixFunction"
     GOTO 9999
  END IF
  !Baseline simulation output, core 0 only
  IExitFLAG = 0 !Do not exit
  IPreviousPrintedIteration = -100!RB ensuring baseline simulation is printed
  IF(my_rank.EQ.0) THEN   
    CALL CreateImagesAndWriteOutput(Iter,IExitFLAG,IErr) 
    IF( IErr.NE.0 ) THEN
      PRINT*,"felixrefine(",my_rank,")error in CreateImagesAndWriteOutput"
      GOTO 9999
    END IF
  END IF

  IF (IRefineMode(12).EQ.1) THEN
    !bisection
	RdeltaUg=0.01!RSimplexLengthScale/100.0!use simplex length scale
	Rtol=0.002! precision 0.01 (in what units, eh? Do find out...)
	DO jnd=1,10!10 cycles to see how it converges
	 DO ind = 1,INoOfVariables!work through Ug components one at a time
	  Iter=Iter+IPrint
	  RpointA=RIndependentVariable(ind)
	  RpointB=RIndependentVariable(ind)+ABS(RdeltaUg*RIndependentVariable(ind))!b must be > a
	  RpointC=ZERO!doesn't matter since this will be the intermediate point returned by mnbrak
	  !bracket the minimum between point A and B
      IF(my_rank.EQ.0) THEN
        PRINT*,"--------------------------------"
        WRITE(SPrintString,FMT='(A20,I2,A4,I3)') "Optimising variable ",ind," of ",INoOfVariables
        PRINT*,TRIM(ADJUSTL(SPrintString))
	    WRITE(SPrintString,FMT='(A14,F8.6,A18,F8.6)') "Initial value ",RpointA,": figure of merit ",RFigureofMerit
        PRINT*,TRIM(ADJUSTL(SPrintString))
		!PRINT*,"Bracketing..."
      END IF	  
	  CALL mnbrak(RIndependentVariable,RpointA,RpointB,RpointC,RfitA,RfitB,RfitC,ind,IErr)
      IF(my_rank.EQ.0) THEN
        !PRINT*,"--------------------------------"
        WRITE(SPrintString,FMT='(A19,F8.6,A1,F8.6,A6,F8.6,A1,F8.6,A1)')&
		"Minimum is between ",RpointA,"(",RfitA,") and ",RpointC,"(",RfitC,")"
        !PRINT*,TRIM(ADJUSTL(SPrintString))
	    WRITE(SPrintString,FMT='(A14,F8.6,A18,F8.6)') "Current value ",RpointB,": figure of merit ",RfitB
        !PRINT*,TRIM(ADJUSTL(SPrintString))
		!PRINT*,"Finding best fit..."
      END IF	  
	  !find the minimum using Brent's method, pass best figure of merit in
	  RFigureofMerit=RfitB
	  CALL BRENT(RFigureofMerit,RIndependentVariable,RpointA,RpointB,RpointC,Rtol,RbestFit,ind,IErr)
	  RIndependentVariable(ind)=RbestFit
      IF(my_rank.EQ.0) THEN
        CALL CreateImagesAndWriteOutput(Iter,IExitFLAG,IErr) 
        IF( IErr.NE.0 ) THEN
          PRINT*,"felixrefine(",my_rank,")error in CreateImagesAndWriteOutput"
          GOTO 9999
        END IF
	    WRITE(SPrintString,FMT='(A12,F8.6,A18,F8.6)') "Final value ",RbestFit,": figure of merit ",RFigureofMerit
        PRINT*,TRIM(ADJUSTL(SPrintString))
      END IF
     END DO	  
	END DO
	
  ELSE
  
    ! Initialise Simplex
    ALLOCATE(RSimplexVariable(INoOfVariables+1,INoOfVariables), STAT=IErr)  
    ALLOCATE(RSimplexFoM(INoOfVariables+1),STAT=IErr)  
    IF( IErr.NE.0 ) THEN
      PRINT*,"felixrefine(",my_rank,")error allocating Simplex variables"
      GOTO 9999
    END IF
    IF(my_rank.EQ.0) THEN
      CALL CreateRandomisedSimplex(RSimplexVariable,RIndependentVariable,IErr)
    END IF
    !=====================================send RSimplexVariable out to all cores
    CALL MPI_BCAST(RSimplexVariable,(INoOfVariables+1)*(INoOfVariables),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IErr)
    !=====================================
    ! Perform initial simplex simulations
    DO ind = 1,(INoOfVariables+1)
      IF((IWriteFLAG.GE.0.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
        PRINT*,"--------------------------------"
        WRITE(SPrintString,FMT='(A8,I2,A4,I3)') "Simplex ",ind," of ",INoOfVariables+1
        PRINT*,TRIM(ADJUSTL(SPrintString))
        PRINT*,"--------------------------------"
      END IF
      CALL SimulateAndFit(RSimplexFoM(ind),RSimplexVariable(ind,:),1,0,IErr)
      IF( IErr.NE.0 ) THEN
        PRINT*,"SimplexInitialisation(",my_rank,") error in SimulateAndFit"
        GOTO 9999
      ENDIF
      IF((IWriteFLAG.GE.0.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
        WRITE(SPrintString,FMT='(A16,F7.5)') "Figure of merit ",RSimplexFoM(ind)
        PRINT*,TRIM(ADJUSTL(SPrintString))
      END IF
    END DO
    ! Apply Simplex Method and iterate
    Iter = 1  
    CALL NDimensionalDownhillSimplex(RSimplexVariable,RSimplexFoM,&
       INoOfVariables+1,INoOfVariables,INoOfVariables,&
       RExitCriteria,Iter,RStandardDeviation,RMean,IErr)
    IF( IErr.NE.0 ) THEN
      PRINT*,"felixrefine(",my_rank,")error in NDimensionalDownhillSimplex"
      GOTO 9999
    ENDIF

  END IF
  !--------------------------------------------------------------------
  ! Deallocate Memory
  !--------------------------------------------------------------------
  DEALLOCATE(CUgMatNoAbs,STAT=IErr)
  DEALLOCATE(CUgMatPrime,STAT=IErr)
  DEALLOCATE(RWeightingCoefficients,STAT=IErr)
  DEALLOCATE(RImageExpi,STAT=IErr)  
  DEALLOCATE(ISymmetryRelations,STAT=IErr)
  DEALLOCATE(IEquivalentUgKey,STAT=IErr)
  DEALLOCATE(CUgToRefine,STAT=IErr)
  DEALLOCATE(RIndividualReflections,STAT=IErr)
  DEALLOCATE(IDisplacements,STAT=IErr)
  DEALLOCATE(ICount,STAT=IErr)
  DEALLOCATE(Rhkl,STAT=IErr)
  DEALLOCATE(RgPoolMag,STAT=IErr)
  DEALLOCATE(RgPoolT,STAT=IErr)
  DEALLOCATE(CUgMat,STAT=IErr)
  DEALLOCATE(RSimulatedPatterns,STAT=IErr)
  IF (IRefineMode(1).EQ.1 .OR. IRefineMode(12).EQ.1) THEN
    DEALLOCATE(RgSumMat,STAT=IErr)
  ELSE
	DEALLOCATE(RAtomCoordinate,STAT=IErr)
  END IF  
  IF( IErr.NE.0 ) THEN
     PRINT*,"felixrefine(",my_rank,")error in final deallocations"
     GOTO 9999
  END IF
  
  !--------------------------------------------------------------------
  ! finish off
  WRITE(my_rank_string,*) my_rank
  CALL SYSTEM_CLOCK(ICurrentTime)
  Duration=REAL(ICurrentTime-IStartTime)/REAL(IRate)
  IHours = FLOOR(Duration/3600.0D0)
  IMinutes = FLOOR(MOD(Duration,3600.0D0)/60.0D0)
  ISeconds = INT(MOD(Duration,3600.0D0)-IMinutes*60)
  IMilliSeconds = INT((Duration-(IHours*3600+IMinutes*60+ISeconds))*1000,IKIND)

  IF(my_rank.EQ.0) THEN
    PRINT*,"--------------------------------"
    WRITE(SPrintString,FMT='(A24,I3,A5,I2,A6,I2,A4)')&
    "Refinement completed in ",IHours," hrs ",IMinutes," mins ",ISeconds," sec"
    PRINT*,TRIM(ADJUSTL(SPrintString))
    PRINT*,"--------------------------------"
    PRINT*,"||||||||||||||||||||||||||||||||"
  END IF
  !--------------------------------------------------------------------
  ! Shut down MPI
  !--------------------------------------------------------------------

9999 &
  CALL MPI_Finalize(IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixrefine(", my_rank, ") error ", IErr, " in MPI_Finalize()"
     STOP
  ENDIF
  
  ! clean shutdown
  STOP
  
END PROGRAM Felixrefine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE mnbrak(RIndependentVariable,Rax,Rbx,Rcx,Rfa,Rfb,Rfc,ind,IErr)
!From Numerical recipes in Fortran section 10.1, adapted for use here

  USE MyNumbers
  
  USE CConst; USE IConst; USE RConst
  USE IPara; USE RPara; USE SPara; USE CPara
  USE BlochPara

  USE IChannels

  USE MPI
  USE MyMPI

  IMPLICIT NONE
  
  REAL(RKIND) :: Rax,Rbx,Rcx,Rfa,Rfb,Rfc,Rgold,RGlimit
  REAL(RKIND) :: Rdum,Rfu,Rq,Rr,Ru,Rulim
  REAL(RKIND),DIMENSION(INoOfVariables) :: RIndependentVariable
  INTEGER(IKIND) :: IErr,Iiter,IExitFLAG,ind
  PARAMETER (Rgold=1.618034, RGlimit=100.0)
  !Given distinct initial points Rax and Rbx this routine searches in the downhill direction
  !and returns new points Rax, Rbx and Rcx that bracket a minimum, with their values Rfa, Rfb, Rfc.
  !Rgold is the (golden) ratio by which intervals are magnified
  !RGlimit is the maximum magnification allowed
  
  Iiter=0!we don't write out while bracketing
  IExitFLAG=0!we never exit from this subroutine
  !Rfa=F(Rax)
  RIndependentVariable(ind)=Rax
  CALL SimulateAndFit(Rfa,RIndependentVariable,Iiter,IExitFLAG,IErr)
  !PRINT*,"a ",Rax,": ",Rfa
  !Rfb=F(Rbx)
  RIndependentVariable(ind)=Rbx
  CALL SimulateAndFit(Rfb,RIndependentVariable,Iiter,IExitFLAG,IErr)
  !PRINT*,"b ",Rbx,": ",Rfb
  IF (Rfb.GT.Rfa) THEN!switch a and b so that a->b is downhill
  !PRINT*,"a->b is uphill,swap"
    Rdum=Rax
	Rax=Rbx
	Rbx=Rdum
	Rdum=Rfb
	Rfb=Rfa
	Rfa=Rdum
  END IF
  Rcx=Rbx+Rgold*(Rbx-Rax)!first guess for c
  !Rfc=F(Rcx)
  RIndependentVariable(ind)=Rcx
  CALL SimulateAndFit(Rfc,RIndependentVariable,Iiter,IExitFLAG,IErr)
  !PRINT*,"c ",Rcx,": ",Rfc
1 IF (Rfb.GE.Rfc) THEN
    Rr=(Rbx-Rax)*(Rfb-Rfc)
    Rq=(Rbx-Rcx)*(Rfb-Rfa)
	Ru=Rbx-((Rbx-Rcx)*Rq-(Rbx-Rax)*Rr)/(TWO*SIGN(MAX(ABS(Rq-Rr),TINY),Rq-Rr))
	Rulim=Rbx+RGlimit*(Rcx-Rbx)
	IF ((Rbx-Ru)*(Ru-Rcx).GT.ZERO) THEN
	  !Rfu=F(Ru)
      RIndependentVariable(ind)=Ru
      CALL SimulateAndFit(Rfu,RIndependentVariable,Iiter,IExitFLAG,IErr)
      !PRINT*,"u ",Ru,": ",Rfu
      IF (Rfu.LT.Rfc) THEN!min is between b and c
        Rax=Rbx
        Rfa=Rfb
        Rbx=Ru
        Rfb=Rfu
        GOTO 2
      ELSE IF(Rfu.GT.Rfb) THEN!min is between a and u
        Rcx=Ru
        Rfc=Rfu
        GOTO 2
      END IF
      Ru=Rcx+Rgold*(Rcx-Rbx)!parabolic fit was no good, default to golden ratio
      !Rfu=F(Ru)
      RIndependentVariable(ind)=Ru
      CALL SimulateAndFit(Rfu,RIndependentVariable,Iiter,IExitFLAG,IErr)
      !PRINT*,"u ",Ru,": ",Rfu
    ELSE IF((Rcx-Ru)*(Ru-Rulim).GT.0) THEN
      !Rfu=F(Ru)
      RIndependentVariable(ind)=Ru
      CALL SimulateAndFit(Rfu,RIndependentVariable,Iiter,IExitFLAG,IErr)
	  !PRINT*,"u ",Ru,": ",Rfu
      IF(Rfu.LT.Rfc) THEN
        Rbx=Rcx
        Rcx=Ru
        Ru=Rcx+Rgold*(Rcx-Rbx)
        Rfb=Rfc
        Rfc=Rfu
        !Rfu=F(Ru)
        RIndependentVariable(ind)=Ru
        CALL SimulateAndFit(Rfu,RIndependentVariable,Iiter,IExitFLAG,IErr)
		!PRINT*,"u ",Ru,": ",Rfu
      ENDIF
    ELSE IF((Ru-Rulim)*(Rulim-Rcx).GE.0) THEN
      Ru=Rulim
      !Rfu=F(Ru)
      RIndependentVariable(ind)=Ru
      CALL SimulateAndFit(Rfu,RIndependentVariable,Iiter,IExitFLAG,IErr)
	  !PRINT*,"u ",Ru,": ",Rfu
    ELSE
      Ru=Rcx+Rgold*(Rcx-Rbx)
      !Rfu=F(Ru)
      RIndependentVariable(ind)=Ru
      CALL SimulateAndFit(Rfu,RIndependentVariable,Iiter,IExitFLAG,IErr)
	  !PRINT*,"u ",Ru,": ",Rfu
    END IF
    Rax=Rbx!shimmy along a bit, do
    Rbx=Rcx
    Rcx=Ru
    Rfa=Rfb
    Rfb=Rfc
    Rfc=Rfu
    GOTO 1
  END IF
  
  !put in order
2  IF (Rax.NE.MIN(Rax,Rbx,Rcx)) THEN
    IF (Rbx.EQ.MIN(Rax,Rbx,Rcx)) THEN!swap a and b
      Rdum=Rax
	  Rax=Rbx
	  Rbx=Rdum
	  Rdum=Rfa
	  Rfa=Rfb
	  Rfb=Rdum
	ELSE!swap a and c
      Rdum=Rax
	  Rax=Rcx
	  Rcx=Rdum
	  Rdum=Rfa
	  Rfa=Rfc
	  Rfc=Rdum
	END IF
  END IF
  IF (Rcx.NE.MAX(Rax,Rbx,Rcx)) THEN!swap b and c
    Rdum=Rbx
	Rbx=Rcx
	Rcx=Rdum
	Rdum=Rfb
	Rfb=Rfc
	Rfc=Rdum
  END IF
  
  RETURN

END SUBROUTINE mnbrak

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE BRENT(Rbrent,RIndependentVariable,Rax,Rbx,Rcx,Rtol,RbestFit,ind,IErr)
  !given a function F(Rx) and a bracketing triplet of abscissas Rax,Rbx,Rcx 
  !where bx is between ax and cx, and F(bx) is less than F(ax) or F(cx)
  !this routine isolates the minimum to a precision of tol using Brent's method
  !position of minimum is RbestFit and its value there is Rbrent
  USE MyNumbers
  USE WriteToScreen

  USE CConst; USE IConst
  USE IPara; USE RPara; USE CPara
  USE BlochPara

  USE IChannels

  USE MPI
  USE MyMPI

  IMPLICIT NONE
  
  INTEGER(IKIND) :: Iiter,Iitmax,IExitFLAG,ind,IErr
  REAL(RKIND) :: Ra,Rb,Rax,Rbx,Rcx,Rd,Re,Rfx,Rfu,Rfv,Rfw,Rp,Rq,Rr,Ru,Rv,Rw,Rx,Rxm
  REAL(RKIND) :: Rtol,Rtol1,Rtol2,RzEPS,ReTemp,RcGold,Rbrent,RbestFit
  REAL(RKIND),DIMENSION(INoOfVariables) :: RIndependentVariable
  PARAMETER (Iitmax=100,RcGold=0.3819660,RzEPS=1.0E-10)

  Iiter=IPreviousPrintedIteration-IPrint!write out while minimising  
  IExitFLAG=0!We never exit felixrefine from this subroutine
  Ra=MIN(Rax,Rcx) 
  Rb=MAX(Rax,Rcx) 
  Rv=Rbx 
  Rw=Rbx 
  Rx=Rbx 
  Re=ZERO
  Rfx=Rbrent
  Rfv=Rbrent 
  Rfw=Rbrent
  DO Iiter=1,Iitmax !main loop
    Rxm=0.5*(Ra+Rb) 
    Rtol1=Rtol*ABS(Rx)+RzEPS 
    Rtol2=2.*Rtol1 
    IF (ABS(Rx-Rxm).LE.(Rtol2-.5*(Rb-Ra))) THEN!Test for done
      GOTO 3
	END IF
    IF( ABS(Re).GT.Rtol1) THEN 
      Rr=(Rx-Rw)*(Rfx-Rfv) 
      Rq=(Rx-Rv)*(Rfx-Rfw) 
      Rp=(Rx-Rv)*Rq-(Rx-Rw)*Rr 
      Rq=2.*(Rq-Rr) 
      IF (Rq.GT.ZERO) THEN
	    Rp=-Rp 
      END IF
      Rq=ABS(Rq) 
      ReTemp=Re 
      Re=Rd 
      IF (ABS(Rp).GE.ABS(.5*Rq*ReTemp).OR.Rp.LE.Rq*(Ra-Rx).OR. Rp.GE.Rq*(Rb-Rx)) GOTO 1 
      Rd=Rp/Rq !parabolic fit
      Ru=Rx+Rd 
      IF (Ru-Ra.LT.Rtol2 .OR. Rb-Ru.LT.Rtol2) THEN 
        Rd=SIGN(Rtol1,Rxm-Rx)
      END IF
      GOTO 2!skip golden section
    END IF 
1   IF (Rx.GE.Rxm) THEN!golden section 
      Re=Ra-Rx 
    ELSE 
      Re=Rb-Rx 
    END IF 
    Rd=RcGold*Re 
2   IF (ABS(Rd).GE.Rtol1) THEN!end of branching, Rd is either golden section or parabolic
      Ru=Rx+Rd 
    ELSE 
      Ru=Rx+SIGN(Rtol1,Rd) 
    END IF 
    RIndependentVariable(ind)=Ru
	CALL SimulateAndFit(Rfu,RIndependentVariable,Iiter,IExitFLAG,IErr)
    !IF(my_rank.EQ.0) THEN
	!  PRINT*,Ru,": ",Rfu
	!END IF
    IF (Rfu.LE.Rfx) THEN 
      IF (Ru.GE.Rx) THEN 
        Ra=Rx 
      ELSE 
        Rb=Rx 
      ENDIF 
      Rv=Rw 
      Rfv=Rfw 
      Rw=Rx 
      Rfw=Rfx 
      Rx=Ru 
      Rfx=Rfu 
    ELSE
      IF(Ru.LT.Rx) THEN
	    Ra=Ru
      ELSE
	    Rb=Ru
      END IF
      IF (Rfu.LE.Rfw .OR. Rw.EQ.Rx) THEN
	    Rv=Rw
		Rfv=Rfw
		Rw=Ru
		Rfw=Rfu
      ELSE IF (Rfu.LE.Rfv .OR. Rv.EQ.Rx .OR. Rv.EQ.Rw) THEN
	    Rv=Ru
	    Rfv=Rfu
      END IF
    END IF
    IF(Ru.LT.Rx) THEN 
      Ra=Ru 
    ELSE 
      Rb=Ru 
    END IF 
    IF(Rfu.LE.Rfw .OR. Rw.EQ.Rx) THEN 
      Rv=Rw 
      Rfv=Rfw 
      Rw=Ru 
      Rfw=Rfu 
    ELSE IF(Rfu.LE.Rfv .OR. Rv.EQ.Rx .OR. Rv.EQ.Rw) THEN 
      Rv=Ru 
      Rfv=Rfu 
    END IF 
  END DO 

  IF (my_rank.EQ.0) THEN
    PRINT*,"Brent exceeded maximum iterations."
    IErr=1
  END IF
  
3 IF (Rfx.LT.Rbrent) THEN!only accept if outgoing fit is better than incoming one
    RbestFit=Rx 
    Rbrent=Rfx 
  ELSE
    RbestFit=Rbx
  END IF
  
  RETURN
  
END SUBROUTINE BRENT

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE AssignArrayLocationsToIterationVariables(IIterativeVariableType,IVariableNo,IArrayToFill,IErr)
!NB IArrayToFill here is equivalent to IIterativeVariableUniqueIDs outside this subroutine
  USE MyNumbers
  
  USE CConst; USE IConst; USE RConst
  USE IPara; USE RPara; USE SPara; USE CPara
  USE BlochPara

  USE IChannels

  USE MPI
  USE MyMPI
  
  IMPLICIT NONE

  INTEGER(IKIND) :: IIterativeVariableType,IVariableNo,IErr,IArrayIndex,&
       IAnisotropicDebyeWallerFactorElementNo
  INTEGER(IKIND),DIMENSION(INoOfVariables,5),INTENT(OUT) :: IArrayToFill  

!!$  Calculate How Many of Each Variable Type There are
!  CALL DetermineNumberofRefinementVariablesPerType(INoofElementsForEachRefinementType,IErr)
  
!!$  Where am I in the Array Right Now?
  IArrayIndex = SUM(INoofElementsForEachRefinementType(:(IIterativeVariableType-1)))+IVariableNo

  SELECT CASE(IIterativeVariableType)

  CASE(1) ! Ugs
     IArrayToFill(IArrayIndex,2) = IIterativeVariableType
     IArrayToFill(IArrayIndex,3) = &
          NINT(REAL(INoofUgs,RKIND)*(REAL(IVariableNo/REAL(INoofUgs,RKIND),RKIND)-&
          CEILING(REAL(IVariableNo/REAL(INoofUgs,RKIND),RKIND)))+REAL(INoofUgs,RKIND))

  CASE(2) ! Coordinates (x,y,z)
     IArrayToFill(IArrayIndex,2) = IIterativeVariableType
     IArrayToFill(IArrayIndex,3) = IVariableNo

  CASE(3) ! Atomic Site Occupancies
     IArrayToFill(IArrayIndex,2) = IIterativeVariableType
     IArrayToFill(IArrayIndex,3) = IAtomicSitesToRefine(IVariableNo)

  CASE(4) ! Isotropic Debye Waller Factors 
     IArrayToFill(IArrayIndex,2) = IIterativeVariableType
     IArrayToFill(IArrayIndex,3) = IAtomicSitesToRefine(IVariableNo)

  CASE(5) ! Anisotropic Debye Waller Factors (a11-a33)
     IArrayToFill(IArrayIndex,2) = IIterativeVariableType
     IArrayToFill(IArrayIndex,3) = IAtomicSitesToRefine(INT(CEILING(REAL(IVariableNo/6.0D0,RKIND))))
     IAnisotropicDebyeWallerFactorElementNo = &
          NINT(6.D0*(REAL(IVariableNo/6.0D0,RKIND)-CEILING(REAL(IVariableNo/6.0D0,RKIND)))+6.0D0)

     SELECT CASE(IAnisotropicDebyeWallerFactorElementNo)

        CASE(1)
           IArrayToFill(IArrayIndex,4:5) = [1,1]
        CASE(2)
           IArrayToFill(IArrayIndex,4:5) = [2,1]
        CASE(3)
           IArrayToFill(IArrayIndex,4:5) = [2,2]
        CASE(4)
           IArrayToFill(IArrayIndex,4:5) = [3,1]
        CASE(5)
           IArrayToFill(IArrayIndex,4:5) = [3,2]
        CASE(6)
           IArrayToFill(IArrayIndex,4:5) = [3,3]

        END SELECT

  CASE(6) ! Lattice Parameters
     IArrayToFill(IArrayIndex,2) = IIterativeVariableType
     IArrayToFill(IArrayIndex,3) = &
          NINT(3.D0*(REAL(IVariableNo/3.0D0,RKIND)-CEILING(REAL(IVariableNo/3.0D0,RKIND)))+3.0D0)
     
  CASE(7) ! Lattice Angles
     IArrayToFill(IArrayIndex,2) = IIterativeVariableType
     IArrayToFill(IArrayIndex,3) = &
          NINT(3.D0*(REAL(IVariableNo/3.0D0,RKIND)-CEILING(REAL(IVariableNo/3.0D0,RKIND)))+3.0D0)

  CASE(8) 
     IArrayToFill(IArrayIndex,2) = IIterativeVariableType

  CASE(9)  
     IArrayToFill(IArrayIndex,2) = IIterativeVariableType

  CASE(10)
     IArrayToFill(IArrayIndex,2) = IIterativeVariableType

  CASE(11)
     IArrayToFill(IArrayIndex,2) = IIterativeVariableType
     
  END SELECT
  
END SUBROUTINE AssignArrayLocationsToIterationVariables

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE RefinementVariableSetup(RIndependentVariable,IErr)
  
  USE MyNumbers
  
  USE CConst; USE IConst; USE RConst
  USE IPara; USE RPara; USE SPara; USE CPara
  USE BlochPara
  
  USE IChannels
  
  USE MPI
  USE MyMPI
  
  IMPLICIT NONE
  
  INTEGER(IKIND) :: IErr,ind,IVariableType
  REAL(RKIND),DIMENSION(INoOfVariables),INTENT(OUT) :: RIndependentVariable
  
  IF((IWriteFLAG.GE.10.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"RefinementVariableSetup(",my_rank,")"
  END IF
  
!!$  Fill the Independent Value array with values

  DO ind = 1,INoOfVariables
     IVariableType = IIterativeVariableUniqueIDs(ind,2)
     SELECT CASE (IVariableType)
     CASE(1)
	    !Structure factor refinement, define in SymmetryRelatedStructureFactorDetermination
    CASE(2)
        RIndependentVariable(ind) = &
             RAllowedVectorMagnitudes(IIterativeVariableUniqueIDs(ind,3))
     CASE(3)
        RIndependentVariable(ind) = &
             RAtomicSitePartialOccupancy(IIterativeVariableUniqueIDs(ind,3))
     CASE(4)
        RIndependentVariable(ind) = &
             RIsotropicDebyeWallerFactors(IIterativeVariableUniqueIDs(ind,3))
     CASE(5)
        RIndependentVariable(ind) = &
             RAnisotropicDebyeWallerFactorTensor(&
             IIterativeVariableUniqueIDs(ind,3),&
             IIterativeVariableUniqueIDs(ind,4),&
             IIterativeVariableUniqueIDs(ind,5))
     CASE(6)
        SELECT CASE(IIterativeVariableUniqueIDs(ind,3))
        CASE(1)
           RIndependentVariable(ind) = RLengthX
        CASE(2)
           RIndependentVariable(ind) = RLengthY
        CASE(3)
           RIndependentVariable(ind) = RLengthZ
        END SELECT
     CASE(7)
        SELECT CASE(IIterativeVariableUniqueIDs(ind,3))
        CASE(1)
           RIndependentVariable(ind) = RAlpha
        CASE(2)
           RIndependentVariable(ind) = RBeta
        CASE(3)
           RIndependentVariable(ind) = RGamma
        END SELECT
     CASE(8)
        RIndependentVariable(ind) = &
             RConvergenceAngle
     CASE(9)
        RIndependentVariable(ind) = &
             RAbsorptionPercentage
     CASE(10)
        RIndependentVariable(ind) = &
             RAcceleratingVoltage
     CASE(11)
        RIndependentVariable(ind) = &
             RRSoSScalingFactor
     CASE(12)
	    !Structure factor refinement, define in SymmetryRelatedStructureFactorDetermination
     END SELECT
  END DO

END SUBROUTINE RefinementVariableSetup
 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE SetupUgsToRefine(IErr)
!Identify unique Ug's and count the number of independent variables INoOfVariables
!using the Hermitian matrix CUgMatNoAbs
!We count over INoofUgs, specified in felix.inp
!The count excludes Ug components that are zero and starts at 1+IUgOffset
!IUgOffset should be at least 1 to exclude U(000), the inner potential
  
  USE MyNumbers
  
  USE CConst; USE IConst; USE RConst
  USE IPara; USE RPara; USE CPara
  USE BlochPara

  USE IChannels

  USE MPI
  USE MyMPI
  
  IMPLICIT NONE 
  
  INTEGER(IKIND) :: IErr,ind,jnd,Iuid,knd,lnd,IUgID
  INTEGER(IKIND),DIMENSION(2) :: ILoc
  CHARACTER*200 :: SPrintString

  IF((IWriteFLAG.GE.5.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"SetupUgsToRefine(",my_rank,")"
  END IF

  !Count equivalent Ugs
  !Equivalent Ug's are identified by the sum of their abs(indices)plus the sum of abs(Ug)'s with no absorption
  ISymmetryRelations = 0_IKIND 
  Iuid=0_IKIND 
  DO ind = 1,nReflections
     DO jnd = 1,ind
        IF(ISymmetryRelations(ind,jnd).NE.0_IKIND) THEN
           CYCLE
        ELSE
          Iuid=Iuid+1_IKIND
          !Fill the symmetry relation matrix with incrementing numbers that have the sign of the imaginary part
          WHERE (ABS(ABS(RgSumMat)-ABS(RgSumMat(ind,jnd))).LE.RTolerance)
            ISymmetryRelations = Iuid*SIGN(1_IKIND,NINT(AIMAG(CUgMatNoAbs)/TINY))
          END WHERE
        END IF
     END DO
  END DO

  IF((IWriteFLAG.GE.0.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     WRITE(SPrintString,FMT='(I5,A25)') Iuid," unique structure factors"
     PRINT*,TRIM(ADJUSTL(SPrintString))
  END IF
  IF(IWriteFLAG.EQ.3.AND.my_rank.EQ.0) THEN
    PRINT*,"Ug matrix:"
    DO ind =1,99
     WRITE(SPrintString,FMT='(2(2X,F5.2,1X,F5.2))') CUgMatNoAbs(ind,1:2)
     PRINT*,TRIM(SPrintString)
    END DO
    !PRINT*,"RgSum matrix:"
    !DO ind =1,50
    ! WRITE(SPrintString,FMT='(8(2X,F5.2))') RgSumMat(ind,1:8)
    ! PRINT*,TRIM(ADJUSTL(SPrintString))
    !END DO
	PRINT*,"hkl: symmetry matrix"
    DO ind =1,200
     WRITE(SPrintString,FMT='(I4,3(1X,I4),A1,8(2X,I4))') ind,NINT(Rhkl(ind,:)),",",ISymmetryRelations(ind,1:2)
     PRINT*,TRIM(SPrintString)
    END DO
  END IF

!Link each key with its Ug, from 1 to the number of unique Ug's Iuid
  ALLOCATE(IEquivalentUgKey(Iuid),STAT=IErr)
  ALLOCATE(CUgToRefine(Iuid),STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"SetupUgsToRefine(",my_rank,")error allocating IEquivalentUgKey or CUgToRefine"
     RETURN
  END IF
  DO ind = 1,Iuid
     ILoc = MINLOC(ABS(ISymmetryRelations-ind))
     IEquivalentUgKey(ind) = ind
     CUgToRefine(ind) = CUgMatNoAbs(ILoc(1),ILoc(2))
    IF (IWriteFLAG.EQ.3.AND.my_rank.EQ.0) THEN
	  WRITE(SPrintString,FMT='(A11,I5,A3,E15.7,1X,E15.7)') "CUgToRefine",ind," : ",CUgToRefine(ind)
      PRINT*,TRIM(ADJUSTL(SPrintString))
	  WRITE(SPrintString,FMT='(A11,I5,A1,I5,A3,E15.7,1X,E15.7)')&
        "CUgMatNoAbs",ILoc(1),",",ILoc(2)," : ",CUgMatNoAbs(ILoc(1),ILoc(2))
      PRINT*,TRIM(ADJUSTL(SPrintString))
    END IF
  END DO

  IF (IWriteFLAG.EQ.6.AND.my_rank.EQ.0) THEN
    PRINT*,"Before sorting:"
    DO ind = 1,50
	  WRITE(SPrintString,FMT='(A11,I3,A3,E15.7,1X,E15.7)') "CUgToRefine",ind," : ",CUgToRefine(ind)
      PRINT*,TRIM(ADJUSTL(SPrintString))
    END DO
  END IF 
  
!Put them in descending order of magnitude  
  CALL ReSortUgs(IEquivalentUgKey,CUgToRefine,Iuid)

  IF (IWriteFLAG.EQ.3.AND.my_rank.EQ.0) THEN
    PRINT*,"After sorting:"
    DO ind = 1,50
	  WRITE(SPrintString,FMT='(A11,I3,A3,I3,A3,E15.7,1X,E15.7)')&
        "CUgToRefine",ind," : ",IEquivalentUgKey(ind)," : ",CUgToRefine(ind)
      PRINT*,TRIM(ADJUSTL(SPrintString))
    END DO
  END IF
  
  !Refinement continuation [hack]
  IF (IContinueFLAG.EQ.1) THEN
    !111
    CUgToRefine(2)=CMPLX(1.687904897922,1.82240814132662)
	IUgID=ISymmetryRelations(2,1)
	WHERE (ABS(ISymmetryRelations-IUgID).EQ.0)
      CUgMatNoAbs = CUgToRefine(2)
    END WHERE
	WHERE (ABS(ISymmetryRelations+IUgID).EQ.0)
      CUgMatNoAbs = CONJG(CUgToRefine(2))
    END WHERE
    !!-2-20
    CUgToRefine(3)=CMPLX(2.3745415053797,0.0)
	IUgID=ISymmetryRelations(8,1)
	WHERE (ABS(ISymmetryRelations-IUgID).EQ.0)
      CUgMatNoAbs = CUgToRefine(3)
    END WHERE
	WHERE (ABS(ISymmetryRelations+IUgID).EQ.0)
      CUgMatNoAbs = CONJG(CUgToRefine(3))
    END WHERE
    !004
    CUgToRefine(4)=CMPLX(1.72720261804284,0.0)
	IUgID=ISymmetryRelations(18,1)
	WHERE (ABS(ISymmetryRelations-IUgID).EQ.0)
      CUgMatNoAbs = CUgToRefine(4)
    END WHERE
	WHERE (ABS(ISymmetryRelations+IUgID).EQ.0)
      CUgMatNoAbs = CONJG(CUgToRefine(4))
    END WHERE
    !11-3
    CUgToRefine(5)=CMPLX(0.94836494246152,1.04943761482216)
	IUgID=ISymmetryRelations(11,1)
	WHERE (ABS(ISymmetryRelations-IUgID).EQ.0)
      CUgMatNoAbs = CUgToRefine(5)
    END WHERE
	WHERE (ABS(ISymmetryRelations+IUgID).EQ.0)
      CUgMatNoAbs = CONJG(CUgToRefine(5))
    END WHERE
    !22-4
    CUgToRefine(6)=CMPLX(1.25676764796887,0.0)
	IUgID=ISymmetryRelations(24,1)
	WHERE (ABS(ISymmetryRelations-IUgID).EQ.0)
      CUgMatNoAbs = CUgToRefine(6)
    END WHERE
	WHERE (ABS(ISymmetryRelations+IUgID).EQ.0)
      CUgMatNoAbs = CONJG(CUgToRefine(6))
    END WHERE
    !440
    CUgToRefine(7)=CMPLX(1.07087394971792,0.0)
	IUgID=ISymmetryRelations(36,1)
	WHERE (ABS(ISymmetryRelations-IUgID).EQ.0)
      CUgMatNoAbs = CUgToRefine(7)
    END WHERE
	WHERE (ABS(ISymmetryRelations+IUgID).EQ.0)
      CUgMatNoAbs = CONJG(CUgToRefine(7))
    END WHERE
    !-3-31
    CUgToRefine(8)=CMPLX(0.689044546419982,0.763792118843367)
	IUgID=ISymmetryRelations(21,1)
	WHERE (ABS(ISymmetryRelations-IUgID).EQ.0)
      CUgMatNoAbs = CUgToRefine(8)
    END WHERE
	WHERE (ABS(ISymmetryRelations+IUgID).EQ.0)
      CUgMatNoAbs = CONJG(CUgToRefine(8))
    END WHERE
    !33-3
    CUgToRefine(9)=CMPLX(0.570891083894945,0.601661598489738)
	IUgID=ISymmetryRelations(33,1)
	WHERE (ABS(ISymmetryRelations-IUgID).EQ.0)
      CUgMatNoAbs = CUgToRefine(9)
    END WHERE
	WHERE (ABS(ISymmetryRelations+IUgID).EQ.0)
      CUgMatNoAbs = CONJG(CUgToRefine(9))
    END WHERE
    !115
    CUgToRefine(10)=CMPLX(0.561361147893234,0.626843766384188)
	IUgID=ISymmetryRelations(28,1)
	WHERE (ABS(ISymmetryRelations-IUgID).EQ.0)
      CUgMatNoAbs = CUgToRefine(10)
    END WHERE
	WHERE (ABS(ISymmetryRelations+IUgID).EQ.0)
      CUgMatNoAbs = CONJG(CUgToRefine(10))
    END WHERE
    !-4-44
    CUgToRefine(11)=CMPLX(0.80107282828933,0.0)
	IUgID=ISymmetryRelations(52,1)
	WHERE (ABS(ISymmetryRelations-IUgID).EQ.0)
      CUgMatNoAbs = CUgToRefine(11)
    END WHERE
	WHERE (ABS(ISymmetryRelations+IUgID).EQ.0)
      CUgMatNoAbs = CONJG(CUgToRefine(11))
    END WHERE
    !008
    CUgToRefine(12)=CMPLX(0.630494960029028,0.0)
	IUgID=ISymmetryRelations(68,1)
	WHERE (ABS(ISymmetryRelations-IUgID).EQ.0)
      CUgMatNoAbs = CUgToRefine(12)
    END WHERE
	WHERE (ABS(ISymmetryRelations+IUgID).EQ.0)
      CUgMatNoAbs = CONJG(CUgToRefine(12))
    END WHERE
    !-3-35
    CUgToRefine(13)=CMPLX(0.412261267514878,0.467054305188535)
	IUgID=ISymmetryRelations(44,1)
	WHERE (ABS(ISymmetryRelations-IUgID).EQ.0)
      CUgMatNoAbs = CUgToRefine(13)
    END WHERE
	WHERE (ABS(ISymmetryRelations+IUgID).EQ.0)
      CUgMatNoAbs = CONJG(CUgToRefine(13))
    END WHERE
    !22-8
    CUgToRefine(14)=CMPLX(0.534308704133266,0.0)
	IUgID=ISymmetryRelations(78,1)
	WHERE (ABS(ISymmetryRelations-IUgID).EQ.0)
      CUgMatNoAbs = CUgToRefine(14)
    END WHERE
	WHERE (ABS(ISymmetryRelations+IUgID).EQ.0)
      CUgMatNoAbs = CONJG(CUgToRefine(14))
    END WHERE
    !-6-60
    CUgToRefine(15)=CMPLX(0.553765683779452,0.0)
	IUgID=ISymmetryRelations(82,1)
	WHERE (ABS(ISymmetryRelations-IUgID).EQ.0)
      CUgMatNoAbs = CUgToRefine(15)
    END WHERE
	WHERE (ABS(ISymmetryRelations+IUgID).EQ.0)
      CUgMatNoAbs = CONJG(CUgToRefine(15))
    END WHERE
    !-1-1-7
    CUgToRefine(16)=CMPLX(0.380276235359774,0.397870473611188)
	IUgID=ISymmetryRelations(56,1)
	WHERE (ABS(ISymmetryRelations-IUgID).EQ.0)
      CUgMatNoAbs = CUgToRefine(16)
    END WHERE
	WHERE (ABS(ISymmetryRelations+IUgID).EQ.0)
      CUgMatNoAbs = CONJG(CUgToRefine(16))
    END WHERE
    !-5-51
    CUgToRefine(17)=CMPLX(0.380109613273699,0.398029660861434)
	IUgID=ISymmetryRelations(58,1)
	WHERE (ABS(ISymmetryRelations-IUgID).EQ.0)
      CUgMatNoAbs = CUgToRefine(17)
    END WHERE
	WHERE (ABS(ISymmetryRelations+IUgID).EQ.0)
      CUgMatNoAbs = CONJG(CUgToRefine(17))
    END WHERE
    !-5-53
    CUgToRefine(18)=CMPLX(0.321160169037219,0.358283860049556)
	IUgID=ISymmetryRelations(64,1)
	WHERE (ABS(ISymmetryRelations-IUgID).EQ.0)
      CUgMatNoAbs = CUgToRefine(18)
    END WHERE
	WHERE (ABS(ISymmetryRelations+IUgID).EQ.0)
      CUgMatNoAbs = CONJG(CUgToRefine(18))
    END WHERE
    !-6-6-4
    CUgToRefine(19)=CMPLX(0.458801755210435,0.0)
	IUgID=ISymmetryRelations(96,1)
	WHERE (ABS(ISymmetryRelations-IUgID).EQ.0)
      CUgMatNoAbs = CUgToRefine(19)
    END WHERE
	WHERE (ABS(ISymmetryRelations+IUgID).EQ.0)
      CUgMatNoAbs = CONJG(CUgToRefine(19))
    END WHERE
    !-3-37
    CUgToRefine(20)=CMPLX(0.283391669499492,0.314015644984719)
	IUgID=ISymmetryRelations(70,1)
	WHERE (ABS(ISymmetryRelations-IUgID).EQ.0)
      CUgMatNoAbs = CUgToRefine(20)
    END WHERE
	WHERE (ABS(ISymmetryRelations+IUgID).EQ.0)
      CUgMatNoAbs = CONJG(CUgToRefine(20))
    END WHERE
    !-4-4-8
    CUgToRefine(21)=CMPLX(0.423656396467032,0.0)
	IUgID=ISymmetryRelations(100,1)
	WHERE (ABS(ISymmetryRelations-IUgID).EQ.0)
      CUgMatNoAbs = CUgToRefine(21)
    END WHERE
	WHERE (ABS(ISymmetryRelations+IUgID).EQ.0)
      CUgMatNoAbs = CONJG(CUgToRefine(21))
    END WHERE
    !-5-55
    CUgToRefine(22)=CMPLX(0.2598195569662342,0.284488486779033)
	IUgID=ISymmetryRelations(84,1)
	WHERE (ABS(ISymmetryRelations-IUgID).EQ.0)
      CUgMatNoAbs = CUgToRefine(22)
    END WHERE
	WHERE (ABS(ISymmetryRelations+IUgID).EQ.0)
      CUgMatNoAbs = CONJG(CUgToRefine(22))
    END WHERE
    !119
    CUgToRefine(23)=CMPLX(0.234458160087768,0.264990601930051)
	IUgID=ISymmetryRelations(92,1)
	WHERE (ABS(ISymmetryRelations-IUgID).EQ.0)
      CUgMatNoAbs = CUgToRefine(23)
    END WHERE
	WHERE (ABS(ISymmetryRelations+IUgID).EQ.0)
      CUgMatNoAbs = CONJG(CUgToRefine(23))
    END WHERE
    !-8-80
    CUgToRefine(24)=CMPLX(0.299705500267044,0.0)
	IUgID=ISymmetryRelations(138,1)
	WHERE (ABS(ISymmetryRelations-IUgID).EQ.0)
      CUgMatNoAbs = CUgToRefine(24)
    END WHERE
	WHERE (ABS(ISymmetryRelations+IUgID).EQ.0)
      CUgMatNoAbs = CONJG(CUgToRefine(24))
    END WHERE
    !771
    CUgToRefine(25)=CMPLX(0.199809232775208,0.224462634318766)
	IUgID=ISymmetryRelations(105,1)
	WHERE (ABS(ISymmetryRelations-IUgID).EQ.0)
      CUgMatNoAbs = CUgToRefine(25)
    END WHERE
	WHERE (ABS(ISymmetryRelations+IUgID).EQ.0)
      CUgMatNoAbs = CONJG(CUgToRefine(25))
    END WHERE
    !339
    CUgToRefine(26)=CMPLX(0.199881022519487,0.224373244214161)
	IUgID=ISymmetryRelations(104,1)
	WHERE (ABS(ISymmetryRelations-IUgID).EQ.0)
      CUgMatNoAbs = CUgToRefine(26)
    END WHERE
	WHERE (ABS(ISymmetryRelations+IUgID).EQ.0)
      CUgMatNoAbs = CONJG(CUgToRefine(26))
    END WHERE
    !55-7
    CUgToRefine(27)=CMPLX(0.199881022519487,0.224373244214161)
	IUgID=ISymmetryRelations(106,1)
	WHERE (ABS(ISymmetryRelations-IUgID).EQ.0)
      CUgMatNoAbs = CUgToRefine(27)
    END WHERE
	WHERE (ABS(ISymmetryRelations+IUgID).EQ.0)
      CUgMatNoAbs = CONJG(CUgToRefine(27))
    END WHERE
    !668
    CUgToRefine(28)=CMPLX(0.300807268960559,0.0)
	IUgID=ISymmetryRelations(152,1)
	WHERE (ABS(ISymmetryRelations-IUgID).EQ.0)
      CUgMatNoAbs = CUgToRefine(28)
    END WHERE
	WHERE (ABS(ISymmetryRelations+IUgID).EQ.0)
      CUgMatNoAbs = CONJG(CUgToRefine(28))
    END WHERE
    !-7-7-3
    CUgToRefine(29)=CMPLX(0.176769844421415,0.195817605460189)
	IUgID=ISymmetryRelations(118,1)
	WHERE (ABS(ISymmetryRelations-IUgID).EQ.0)
      CUgMatNoAbs = CUgToRefine(29)
    END WHERE
	WHERE (ABS(ISymmetryRelations+IUgID).EQ.0)
      CUgMatNoAbs = CONJG(CUgToRefine(29))
    END WHERE
    !0012
    CUgToRefine(30)=CMPLX(0.28147196375376,0.0)
	IUgID=ISymmetryRelations(164,1)
	WHERE (ABS(ISymmetryRelations-IUgID).EQ.0)
      CUgMatNoAbs = CUgToRefine(30)
    END WHERE
	WHERE (ABS(ISymmetryRelations+IUgID).EQ.0)
      CUgMatNoAbs = CONJG(CUgToRefine(30))
    END WHERE
    !884
    CUgToRefine(31)=CMPLX(0.266581374365749,0.0)
	IUgID=ISymmetryRelations(160,1)
	WHERE (ABS(ISymmetryRelations-IUgID).EQ.0)
      CUgMatNoAbs = CUgToRefine(31)
    END WHERE
	WHERE (ABS(ISymmetryRelations+IUgID).EQ.0)
      CUgMatNoAbs = CONJG(CUgToRefine(31))
    END WHERE
    PRINT*,"After update:"
    DO ind = 1,50
	  WRITE(SPrintString,FMT='(A11,I3,A3,I3,A3,E15.7,1X,E15.7)')&
        "CUgToRefine",ind," : ",IEquivalentUgKey(ind)," : ",CUgToRefine(ind)
      PRINT*,TRIM(ADJUSTL(SPrintString))
    END DO
  END IF
  
  !Count the number of Independent Variables
  jnd=1
  DO ind = 1+IUgOffset,INoofUgs+IUgOffset !=== temp comment out !=== for real part only***
    IF ( ABS(REAL(CUgToRefine(ind),RKIND)).GE.RTolerance ) THEN!===
      jnd=jnd+1
	END IF!===
    IF ( ABS(AIMAG(CUgToRefine(ind))).GE.RTolerance ) THEN!===
      jnd=jnd+1!===
	END IF!===
  END DO
  INoOfVariables = jnd![[[-1 !===the last increment is for absorption ![[[ delete the -1 to include absorption***
  
END SUBROUTINE SetupUgsToRefine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE CreateRandomisedSimplex(RSimplexVariable,RIndependentVariable,IErr)

USE MyNumbers
  
  USE CConst; USE IConst; USE RConst
  USE IPara; USE RPara; USE SPara; USE CPara
  USE BlochPara 
  USE IChannels
  USE MPI
  USE MyMPI
  
  IMPLICIT NONE
  
  INTEGER(IKIND) :: IErr,ind
  REAL(RKIND),DIMENSION(:),ALLOCATABLE :: RRandomSigns,RRandomNumbers
  REAL(RKIND),DIMENSION(INoOfVariables+1,INoOfVariables),INTENT(OUT) :: RSimplexVariable
  REAL(RKIND),DIMENSION(INoOfVariables),INTENT(INOUT) :: RIndependentVariable
  
  IF(IRefineMode(2).EQ.1) THEN!it's an atomic coordinate refinement
     DO ind = 1,(INoOfVariables+1)
        ALLOCATE(RRandomSigns(IAllowedVectors),RRandomNumbers(IAllowedVectors),&
             STAT=IErr)       
        
!!$           Randomise Atomic Displacements
        CALL RandomSequence(RRandomNumbers,IAllowedVectors,ind,IErr)
        CALL RandomSequence(RRandomSigns,IAllowedVectors,2*ind,IErr)
        WHERE (RRandomSigns.LT.HALF)
           RRandomSigns = ONE
        ELSEWHERE
           RRandomSigns = -ONE
        END WHERE
        RSimplexVariable(ind,:IAllowedVectors) = &
             RRandomNumbers*RRandomSigns*RSimplexLengthScale
        DEALLOCATE(RRandomSigns,RRandomNumbers) 
        ALLOCATE(RRandomSigns(INoOfVariables-IAllowedVectors),&
             RRandomNumbers(INoOfVariables-IAllowedVectors),STAT=IErr)
        
!!$           Randomise Everything else
        CALL RandomSequence(RRandomNumbers,&
             INoOfVariables-IAllowedVectors,ind,IErr)
        CALL RandomSequence(RRandomSigns,&
             INoOfVariables-IAllowedVectors,2*ind,IErr)
        WHERE (RRandomSigns.LT.HALF)
           RRandomSigns = ONE
        ELSEWHERE
           RRandomSigns = -ONE
        END WHERE
        RSimplexVariable(ind,(IAllowedVectors+1):) = &
             RIndependentVariable((IAllowedVectors+1):)*&
             (1+(RRandomNumbers*RRandomSigns*RSimplexLengthScale))
        DEALLOCATE(RRandomSigns,RRandomNumbers)
        
     END DO
     
  ELSE
     ALLOCATE(RRandomSigns(INoOfVariables),&
          RRandomNumbers(INoOfVariables),STAT=IErr)
     
     DO ind = 1,(INoOfVariables+1)
        CALL RandomSequence(RRandomNumbers,INoOfVariables,ind,IErr)
        CALL RandomSequence(RRandomSigns,INoOfVariables,2*ind,IErr)
        WHERE (RRandomSigns.LT.HALF)
           RRandomSigns = ONE
        ELSEWHERE
           RRandomSigns = -ONE
        END WHERE
        RSimplexVariable(ind,:) = &
             RIndependentVariable(:)*&
             (1+(RRandomNumbers*RRandomSigns*RSimplexLengthScale))
     END DO
        DEALLOCATE(RRandomSigns,RRandomNumbers)
     
  END IF
  

END SUBROUTINE CreateRandomisedSimplex

!!$  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE InitialiseAtomicVectorMagnitudes(IVariableID,RCorrectedMovement,IErr)
  
!!$  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!$  % Creates pseudo random movements of atoms using allowed vectors
!!$  % to initialise the simplex, proposed movements which exit the unit
!!$  $ cell are corrected to bring the atom back in on the opposite side
!!$  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  USE MyNumbers
  
  USE CConst; USE IConst; USE RConst
  USE IPara; USE RPara; USE SPara; USE CPara
  USE BlochPara
  
  USE IChannels
  
  USE MPI
  USE MyMPI
  
  IMPLICIT NONE

  INTEGER(IKIND) :: IErr,ind,IVariableID
  REAL(RKIND) :: RNegativeMovement,RPositiveMovement,RCorrectedMovement,RANDOMNUMBER

  RNegativeMovement = RSimplexLengthScale*(-1.0_RKIND)
  RPositiveMovement = RSimplexLengthScale
!RB this check can be done in less lines than it takes to call the subroutine
  IF(RANDOMNUMBER(IVariableID,IErr).LT.0.5_RKIND) THEN
     CALL OutofUnitCellCheck(IVariableID,RNegativeMovement,RCorrectedMovement,IErr)
  ELSE
     CALL OutofUnitCellCheck(IVariableID,RPositiveMovement,RCorrectedMovement,IErr)
  END IF

END SUBROUTINE InitialiseAtomicVectorMagnitudes

!!$  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE RandomSequence(RRandomSequence,IRandomSequenceLength,ISeedModifier,IErr)
!!$  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!$  % Sets up a pseudo random sequence and selects a number

  USE MyNumbers
  
  USE CConst; USE IConst; USE RConst
  USE IPara; USE RPara; USE SPara; USE CPara
  USE BlochPara

  USE IChannels

  USE MPI
  USE MyMPI

  IMPLICIT NONE

  INTEGER(IKIND) :: IErr,Ivalues(1:8), k,IRandomSequenceLength,ISeedModifier
  INTEGER(IKIND), DIMENSION(:), ALLOCATABLE :: seed
  REAL(RKIND),DIMENSION(IRandomSequenceLength) :: RRandomSequence
  
  CALL DATE_AND_TIME(VALUES=Ivalues)

  IValues = IValues*ISeedModifier
!!$  CALL SYSTEM_CLOCK(
  IF (IRandomFLAG.EQ.0) THEN
     CALL RANDOM_SEED(size=k)
     allocate(seed(1:k))
     seed(:) = Ivalues(8)
     CALL RANDOM_SEED(put=seed)
  ELSE
     CALL RANDOM_SEED(size=k)
     ALLOCATE(seed(1:k))
     seed(:) = IFixedSeed*ISeedModifier
     CALL RANDOM_SEED(put=seed)
  END IF
   
  DEALLOCATE(seed)

  CALL RANDOM_NUMBER(RRandomSequence)
  
!!$  RANDOMSEQUENCE = RRandomNumberSequence(IRequestedNumber)
  
END SUBROUTINE  RANDOMSEQUENCE

!!$  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

REAL(RKIND) FUNCTION RANDOMNUMBER(IRequestedNumber,IErr)
!!$  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!$  % Sets up a pseudo random sequence and selects a number
!!$  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  USE MyNumbers
  
  USE CConst; USE IConst; USE RConst
  USE IPara; USE RPara; USE SPara; USE CPara
  USE BlochPara

  USE IChannels

  USE MPI
  USE MyMPI

  IMPLICIT NONE

  INTEGER(IKIND) :: IErr,values(1:8), k,IRequestedNumber
  INTEGER(IKIND), DIMENSION(:), ALLOCATABLE :: seed
  REAL(RKIND),DIMENSION(IRequestedNumber) :: RRandomNumberSequence
  
  CALL DATE_AND_TIME(values=values)
  
  IF (IRandomFLAG.EQ.0) THEN
     CALL RANDOM_SEED(size=k)
     allocate(seed(1:k))
     seed(:) = values(8)
     CALL RANDOM_SEED(put=seed)
  ELSE
     CALL RANDOM_SEED(size=k)
     ALLOCATE(seed(1:k))
     seed(:) = IFixedSeed*IRequestedNumber
     CALL RANDOM_SEED(put=seed)
  END IF
   
  CALL RANDOM_NUMBER(RRandomNumberSequence)
  
  RANDOMNUMBER = RRandomNumberSequence(IRequestedNumber)
  
END FUNCTION RANDOMNUMBER

!!$  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE OutofUnitCellCheck(IVariableID,RProposedMovement,RCorrectedMovement,IErr)
!!$  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!$  % Checks that vector movement applied by the simplex initialisation
!!$  % does not move an atom out fo the unit cell, and if it does
!!$  % the atom is moved back into the unit cell on the opposite side
!!$  % as if the atom had moved from one unit cell into the neighbouring one
!!$  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  USE MyNumbers
  
  USE CConst; USE IConst; USE RConst
  USE IPara; USE RPara; USE SPara; USE CPara
  USE BlochPara

  USE IChannels

  USE MPI
  USE MyMPI

  IMPLICIT NONE
  
  INTEGER(IKIND) :: ind,IErr,IVariableID,IAtomID,IVectorID
  REAL(RKIND),DIMENSION(THREEDIM) :: RProposedAtomicCoordinate,RDummyMovement
  REAL(RKIND),INTENT(IN) :: RProposedMovement
  REAL(RKIND),INTENT(OUT) :: RCorrectedMovement
  
  IVectorID = IIterativeVariableUniqueIDs(IVariableID,3)
  
  IAtomID = IAllowedVectorIDs(IVectorID)
 
  RProposedAtomicCoordinate(:) = RAtomSiteFracCoordVec(IAtomID,:) + &
       RProposedMovement*RAllowedVectors(IVectorID,:)

  RDummyMovement = RProposedMovement

  IF(ANY(RProposedAtomicCoordinate.GT.ONE).OR.ANY(RProposedAtomicCoordinate.LT.ZERO)) THEN
     DO ind = 1,THREEDIM
        IF (RProposedAtomicCoordinate(ind).GT.ONE) THEN
           RDummyMovement(ind) = (ONE-RAtomSiteFracCoordVec(IAtomID,ind))/RAllowedVectors(IVectorID,ind)
        ELSEIF(RProposedAtomicCoordinate(ind).LT.ZERO) THEN
           RDummyMovement(ind) = (-RAtomSiteFracCoordVec(IAtomID,ind))/RAllowedVectors(IVectorID,ind)
        ELSE
           RDummyMovement(ind) = RProposedMovement
        END IF
     END DO
  END IF

  IF(RProposedMovement.LT.ZERO) THEN
     RCorrectedMovement = MAXVAL(RDummyMovement)
  ELSE
     RCorrectedMovement = MINVAL(RDummyMovement)
  END IF

END SUBROUTINE OutofUnitCellCheck

!!$  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE CreateIdentityMatrix(IIdentityMatrix,ISize,IErr)

!!$  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!$  % This Subroutine creates an identity matrix of size
!!$  % ISize * ISize
!!$  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!why do we have this????
  USE MyNumbers
  
  USE CConst; USE IConst; USE RConst
  USE IPara; USE RPara; USE SPara; USE CPara
  USE BlochPara

  USE IChannels

  USE MPI
  USE MyMPI

  IMPLICIT NONE
  
  INTEGER(IKIND) :: &
       IErr,ISize,ind
  INTEGER(IKIND),DIMENSION(ISize,ISize) :: &
       IIdentityMatrix

  IIdentityMatrix = 0

  DO ind = 1,ISize
     IIdentityMatrix(ind,ind) = 1
  END DO

END SUBROUTINE CreateIdentityMatrix

!!$  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE RecoverSavedSimplex(RSimplexVariable,RSimplexFoM,RStandardDeviation,RMean,Iter,IErr)

!!$  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!$  % This Subroutine reads the fr-simplex.txt file from a previous
!!$  % refinement run, and recreates the simplex volume and tolerances
!!$  % allowing for the continuation of a previous refinement
!!$  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  USE MyNumbers
  
  USE CConst; USE IConst; USE RConst
  USE IPara; USE RPara; USE SPara; USE CPara
  USE BlochPara

  USE IChannels

  USE MPI
  USE MyMPI

  IMPLICIT NONE

  INTEGER(IKIND) :: &
       IErr,ind,Iter
  REAL(RKIND),DIMENSION(INoOfVariables+1,INoOfVariables) :: RSimplexVariable
  REAL(RKIND),DIMENSION(INoOfVariables+1) :: RSimplexFoM
  REAL(RKIND) :: RStandardDeviation,RMean
  CHARACTER*200 :: CSizeofData,SFormatString,filename

  WRITE(filename,*) "fr-Simplex.txt"

  OPEN(UNIT=IChOutSimplex,STATUS='UNKNOWN',&
        FILE=TRIM(ADJUSTL(filename)))
  
  WRITE(CSizeofData,*) INoOfVariables+1
  WRITE(SFormatString,*) "("//TRIM(ADJUSTL(CSizeofData))//"(1F6.3,1X),A1)"

  DO ind = 1,(INoOfVariables+1)
     READ(IChOutSimplex,FMT=SFormatString) RSimplexVariable(ind,:),RSimplexFoM(ind)
  END DO
    
  READ(IChOutSimplex,FMT="(2(1F6.3,1X),I5.1,I5.1,A1)") RStandardDeviation,RMean,IStandardDeviationCalls,Iter

  CLOSE(IChOutSimplex)

END SUBROUTINE RecoverSavedSimplex

!!$  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE SetupAtomicVectorMovements(IErr)

  USE MyNumbers
  
  USE CConst; USE IConst; USE RConst
  USE IPara; USE RPara; USE SPara; USE CPara
  USE BlochPara

  USE IChannels

  USE MPI
  USE MyMPI
  
  IMPLICIT NONE

  INTEGER(IKIND) :: IErr,knd,jnd,ind,ISpaceGrp
  INTEGER(IKIND),DIMENSION(:),ALLOCATABLE :: IVectors
  
  CALL ConvertSpaceGroupToNumber(ISpaceGrp,IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"SetupAtomicVectorMovements(",my_rank,")error in ConvertSpaceGroupToNumber"
     RETURN
  ENDIF

  ALLOCATE(IVectors(SIZE(SWyckoffSymbols)),STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"felixrefine (", my_rank, ") error in Allocation() of IVectors"
     RETURN
  ENDIF
  
  DO ind = 1,SIZE(SWyckoffSymbols)!NB SIZE(SWyckoffSymbols)=IAtomicSitesToRefine
     CALL CountAllowedMovements(ISpaceGrp,SWyckoffSymbols(ind),IVectors(ind),IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"SetupAtomicVectorMovements(",my_rank,")error in CountAllowedMovements "
        RETURN
     ENDIF    
  END DO
  
  IAllowedVectors = SUM(IVectors)
  
  ALLOCATE(IAllowedVectorIDs(IAllowedVectors),STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"SetupAtomicVectorMovements(", my_rank, ") error in Allocation() of IAllowedVectorIDs"
     RETURN
  ENDIF
  
  knd = 0
  
  DO ind = 1,SIZE(SWyckoffSymbols)
     DO jnd = 1,IVectors(ind)
        knd = knd + 1
        IAllowedVectorIDs(knd) = IAtomicSitesToRefine(ind)
     END DO
  END DO
  
  ALLOCATE(RAllowedVectors(IAllowedVectors,THREEDIM),STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"SetupAtomicVectorMovements(",my_rank,")error Allocation  RAllowedVectors"
     RETURN
  ENDIF
  ALLOCATE(RAllowedVectorMagnitudes(IAllowedVectors),STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"SetupAtomicVectorMovements(",my_rank,")error Allocation  RAllowedVectorMagnitudes"
     RETURN
  ENDIF
  
  RAllowedVectorMagnitudes = ZERO
  
  DO ind = 1,SIZE(SWyckoffSymbols)
     CALL DetermineAllowedMovements(ISpaceGrp,SWyckoffSymbols(ind),&
          RAllowedVectors(SUM(IVectors(:(ind-1)))+1:SUM(IVectors(:(ind))),:),&
          IVectors(ind),IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"SetupAtomicVectorMovements(",my_rank,")error in DetermineAllowedMovements"
        RETURN
     ENDIF
     
  END DO
  
  !--------------------------------------------------------------------
  ! Save Atomic Coordinates  
  !--------------------------------------------------------------------
  
  ALLOCATE(RInitialAtomSiteFracCoordVec(&
       SIZE(RAtomSiteFracCoordVec,DIM=1),SIZE(RAtomSiteFracCoordVec,DIM=2)),STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"SetupAtomicVectorMovements(",my_rank,")error ALLOCATE RInitialAtomSiteFracCoordVec "
     RETURN
  ENDIF
  
  RInitialAtomSiteFracCoordVec = RAtomSiteFracCoordVec
  
END SUBROUTINE SetupAtomicVectorMovements
