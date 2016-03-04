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
       ind,jnd,knd,ICalls,IIterationCount,ICutOff,IHOLZgPoolMag,IBSMaxLocGVecAmp,&
	   ILaueLevel,INumTotalReflections,ITotalLaueZoneLevel,INhkl,Iuid,IFinishFLAG,&
	   INumInitReflections,IZerothLaueZoneLevel,INumFinalReflections,IThicknessIndex
  INTEGER(IKIND) :: IStartTime,ICurrentTime,IRate
  INTEGER(IKIND),DIMENSION(2) :: ILoc
  INTEGER(IKIND), DIMENSION(:),ALLOCATABLE :: IOriginGVecIdentifier
  REAL(RKIND) :: StartTime,CurrentTime,Duration,&
       RFigureOfMerit,RHOLZAcceptanceAngle,RLaueZoneGz
  REAL(RKIND), DIMENSION(:,:), ALLOCATABLE :: RSimplexVariable
  REAL(RKIND),DIMENSION(:),ALLOCATABLE :: RSimplexFoM,RIndependentVariable
  REAL(RKIND), DIMENSION(:,:), ALLOCATABLE :: RgDummyVecMat,RgPoolMagLaue
  REAL(RKIND) :: RBCASTREAL,RStandardDeviation,RMean,RGzUnitVec,RMinLaueZoneValue,&
       RMaxLaueZoneValue,RMaxAcceptanceGVecMag,RLaueZoneElectronWaveVectorMag
  CHARACTER*40 :: my_rank_string
  CHARACTER*20 :: Sind
  CHARACTER*200 :: SPrintString
  LOGICAL :: LInitialSimulationFLAG = .TRUE.

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
  ! on screen output
  IF((IWriteFLAG.GE.0.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"--------------------------------------------------------------"
     PRINT*,"Felixrefine: ", RStr
     PRINT*,"          ", DStr
     PRINT*,"          ", AStr
     PRINT*,"          on rank= ", my_rank, " of ", p, " in total."
     PRINT*,"--------------------------------------------------------------"
  END IF
  ISoftwareMode =2 ! felixrefinemode

  !--------------------------------------------------------------------
  ! timing startup
  CALL SYSTEM_CLOCK(count_rate=IRate)
  CALL SYSTEM_CLOCK(IStarttime)
  
  !--------------------------------------------------------------------
  ! read input files
  CALL ReadInpFile(IErr)
  CALL ReadHklFile(IErr)
  CALL ReadScaFile(IErr)
  CALL ReadCif(IErr)
  !--allocate--
  ALLOCATE(RImageExpi(2*IPixelCount,2*IPixelCount,INoOfLacbedPatterns),STAT=IErr) 
  CALL ReadExperimentalImages(IErr)
  
  !--------------------------------------------------------------------
  ! some basic calculations  (was MicroscopySettings)
  RElectronVelocity= RSpeedOfLight * SQRT( ONE - ( (RElectronMass*RSpeedOfLight**2) / &
       (RElectronCharge*RAcceleratingVoltage*1000.0_RKIND + &
        RElectronMass*RSpeedOfLight**2) )**2 )
  RElectronWaveLength= RPlanckConstant / &
       ( SQRT(TWO*RElectronMass*RElectronCharge*RAcceleratingVoltage*1000.0_RKIND) * &
         SQRT(ONE + (RElectronCharge*RAcceleratingVoltage*1000.0_RKIND) / &
          (TWO*RElectronMass*RSpeedOfLight**2) )) * RAngstromConversion
  RElectronWaveVectorMagnitude=TWOPI/RElectronWaveLength
  RRelativisticCorrection= ONE/SQRT(ONE-(RElectronVelocity/RSpeedOfLight)**2 )
  RRelativisticMass= RRelativisticCorrection*RElectronMass  
  
  CALL ReciprocalLattice(IErr) 

  !**********************************
  INAtomsUnitCell=8
  ITotalAtoms=8
  ALLOCATE(IAtoms(ITotalAtoms),STAT=IErr)
  IAtoms(1)=31
  IAtoms(2)=31
  IAtoms(3)=31
  IAtoms(4)=31
  IAtoms(5)=33
  IAtoms(6)=33
  IAtoms(7)=33
  IAtoms(8)=33
  ALLOCATE(RAtomCoordinate(ITotalAtoms,THREEDIM),STAT=IErr)
  RAtomCoordinate(1,1) = ZERO
  RAtomCoordinate(1,2) = ZERO
  RAtomCoordinate(1,3) = ZERO
  RAtomCoordinate(2,1) = 0.5D0
  RAtomCoordinate(2,2) = 0.5D0
  RAtomCoordinate(2,3) = ZERO
  RAtomCoordinate(3,1) = 0.5D0
  RAtomCoordinate(3,2) = ZERO
  RAtomCoordinate(3,3) = 0.5D0
  RAtomCoordinate(4,1) = ZERO
  RAtomCoordinate(4,2) = 0.5D0
  RAtomCoordinate(4,3) = 0.5D0
  RAtomCoordinate(5,1) = 0.25D0
  RAtomCoordinate(5,2) = 0.25D0
  RAtomCoordinate(5,3) = 0.25D0
  RAtomCoordinate(6,1) = 0.75D0
  RAtomCoordinate(6,2) = 0.75D0
  RAtomCoordinate(6,3) = 0.25D0
  RAtomCoordinate(7,1) = 0.75D0
  RAtomCoordinate(7,2) = 0.25D0
  RAtomCoordinate(7,3) = 0.75D0
  RAtomCoordinate(8,1) = 0.25D0
  RAtomCoordinate(8,2) = 0.75D0
  RAtomCoordinate(8,3) = 0.75D0
  ALLOCATE(ROcc(ITotalAtoms),STAT=IErr)
  ROcc(1)=1
  ROcc(2)=1
  ROcc(3)=1
  ROcc(4)=1
  ROcc(5)=1
  ROcc(6)=1
  ROcc(7)=1
  ROcc(8)=1
  ALLOCATE(RDWF(ITotalAtoms),STAT=IErr)
  RDWF(1)=0.4
  RDWF(2)=0.4
  RDWF(3)=0.4
  RDWF(4)=0.4
  RDWF(5)=0.3
  RDWF(6)=0.3
  RDWF(7)=0.3
  RDWF(8)=0.3
  !**********************************  
 
! Count the reflections that make up the pool of g-vectors 
  RHOLZAcceptanceAngle=TWODEG2RADIAN  
  IHKLMAXValue = 5!RB starting value, increments in loop below
  CALL HKLCount(IHKLMAXValue,RZDirC,INhkl,RHOLZAcceptanceAngle,IErr)
  DO WHILE (INhkl.LT.IMinReflectionPool) 
     IHKLMAXValue = IHKLMAXValue*2
     CALL HKLCount(IHKLMAXValue,RZDirC,INhkl,RHOLZAcceptanceAngle,IErr)
  END DO  

  ! Fill the list of reflections Rhkl as indices h,k,l
  ALLOCATE(Rhkl(INhkl,THREEDIM),STAT=IErr)
  CALL HKLMake(IHKLMAXValue,RZDirC,RHOLZAcceptanceAngle,IErr)
  CALL SortHKL(Rhkl,INhkl,IErr)   
  CALL SpecificReflectionDetermination (IErr)
  ALLOCATE(RgPoolT(INhkl,THREEDIM),STAT=IErr)

  !Calculate the g vector list RgPoolT in reciprocal angstrom units (in the microscope reference frame?)
  ICutOff = 1
  DO ind=1,INhkl
    DO jnd=1,THREEDIM
      RgPoolT(ind,jnd) = Rhkl(ind,1)*RarVecM(jnd) + &
      Rhkl(ind,2)*RbrVecM(jnd) + Rhkl(ind,3)*RcrVecM(jnd)
    ENDDO
	!If a g-vector has a non-zero z-component it is not in the ZOLZ
    IF((RgPoolT(ind,3).GT.TINY.OR.RgPoolT(ind,3).LT.-TINY).AND.ICutOff.NE.0) THEN
        RGzUnitVec=ABS(RgPoolT(ind,3))
        ICutOff=0
    END IF
  ENDDO
  
  !g-vector magnitudes
  ALLOCATE(RgPoolMag(INhkl),STAT=IErr)
  DO ind=1,INhkl
     RgPoolMag(ind)= SQRT(DOT_PRODUCT(RgPoolT(ind,:),RgPoolT(ind,:)))
  ENDDO
  RBSMaxGVecAmp = RgPoolMag(IMinReflectionPool)
  
  ALLOCATE(RgVecVec(INhkl),STAT=IErr)
  CALL DiffractionPatternCalculation(IErr)  
  
  IThicknessCount= (RFinalThickness-RInitialThickness)/RDeltaThickness + 1

  !count reflections up to cutoff magnitude
  nReflections = 0
  DO ind=1,INhkl
     IF (ABS(RgPoolMag(ind)).LE.RBSMaxGVecAmp) THEN
        nReflections = nReflections + 1
     END IF
  ENDDO
  nStrongBeams = nReflections
  nWeakBeams = 0

  ALLOCATE(CUgMatNoAbs(nReflections,nReflections),STAT=IErr)
  ALLOCATE(CUgMat(nReflections,nReflections),STAT=IErr)
  ALLOCATE(RgSumMat(nReflections,nReflections),STAT=IErr) 
  ALLOCATE(RgMatMat(nReflections,nReflections,THREEDIM),STAT=IErr)
  ALLOCATE(RgMatMag(nReflections,nReflections),STAT=IErr)
  
  CALL GMatrixInitialisation (IErr)
  CALL StructureFactorInitialisation (IErr)
 !   DO ind =1,6
 !    WRITE(SPrintString,FMT='(10(1X,F5.2))') CUgMatNoAbs(ind,1:5)
 !    PRINT*,TRIM(ADJUSTL(SPrintString))
 !   END DO
  DEALLOCATE(RgMatMag,STAT=IErr)
  DEALLOCATE(RgMatMat,STAT=IErr)
  
  ALLOCATE(ISymmetryRelations(nReflections,nReflections),STAT=IErr)  
  !Count equivalent Ugs
  !Equivalent Ug's are identified by the sum of their abs(indices)plus the sum of abs(Ug)'s with no absorption
  RgSumMat = RgSumMat+ABS(REAL(CUgMatNoAbs))+ABS(AIMAG(CUgMatNoAbs))
  ISymmetryRelations = 0_IKIND 
  Iuid = 0_IKIND 
  DO ind = 1,nReflections
     DO jnd = 1,ind
        IF(ISymmetryRelations(ind,jnd).NE.0) THEN
           CYCLE
        ELSE
           Iuid = Iuid + 1_IKIND
           !Ug Fill the symmetry relation matrix with incrementing numbers that have the sign of the imaginary part
		   WHERE (ABS(RgSumMat-ABS(RgSumMat(ind,jnd))).LE.RTolerance)
              ISymmetryRelations = Iuid*SIGN(1_IKIND,NINT(AIMAG(CUgMatNoAbs)/TINY**2))
           END WHERE
        END IF
     END DO
  END DO
  WRITE(SPrintString,FMT='(I5,A25)') Iuid," unique structure factors"
  PRINT*,TRIM(ADJUSTL(SPrintString))
  
  ALLOCATE(IEquivalentUgKey(Iuid),STAT=IErr)
  ALLOCATE(CUgToRefine(Iuid),STAT=IErr)  
  
  DO ind = 1,Iuid
     ILoc = MINLOC(ABS(ISymmetryRelations-ind))
     IEquivalentUgKey(ind) = ind
     CUgToRefine(ind) = CUgMatNoAbs(ILoc(1),ILoc(2))
  END DO

  CALL ReSortUgs(IEquivalentUgKey,CUgToRefine,Iuid)
  
  jnd=1
  DO ind = 2,INoofUgs+1!RB ignore the first one as it is the internal potential
    IF ( ABS(REAL(CUgToRefine(ind),RKIND)).GE.RTolerance ) THEN
      jnd=jnd+1
	END IF
    IF ( ABS(AIMAG(CUgToRefine(ind))).GE.RTolerance ) THEN
      jnd=jnd+1
	END IF
  END DO
  INoOfVariables = jnd
  
  IF(my_rank.EQ.0) THEN
    IF ( INoOfVariables.EQ.1 ) THEN 
      PRINT*,"Only one independent variable"
	ELSE
      WRITE(SPrintString,FMT='(I3,1X,A21))') INoOfVariables,"independent variables"
      PRINT*,TRIM(ADJUSTL(SPrintString))
    END IF
  END IF
  
  ALLOCATE(RIndependentVariable(INoOfVariables),STAT=IErr)
    !Fill up the IndependentVariable list with CUgMatNoAbs components
    jnd=1
    DO ind = 2,INoofUgs+1!start from 2 since the first one is the inner potential
      IF ( ABS(REAL(CUgToRefine(ind),RKIND)).GE.RTolerance ) THEN
        RIndependentVariable(jnd) = REAL(CUgToRefine(ind),RKIND)
        jnd=jnd+1
	  END IF
      IF ( ABS(AIMAG(CUgToRefine(ind))).GE.RTolerance ) THEN
        RIndependentVariable(jnd) = AIMAG(CUgToRefine(ind))
        jnd=jnd+1
      END IF
    END DO
    RIndependentVariable(jnd) = RAbsorptionPercentage!RB absorption always included in structure factor refinement as last variable

  ALLOCATE(RSimplexVariable(INoOfVariables+1,INoOfVariables), STAT=IErr)  
  ALLOCATE(RSimplexFoM(INoOfVariables+1),STAT=IErr) 
  
  ! Setup Images for output
  ALLOCATE(RhklPositions(nReflections,2),STAT=IErr)
  CALL ImageSetup(IErr)
  ALLOCATE(RSimulatedPatterns(INoOfLacbedPatterns,IThicknessCount,IPixelTotal),STAT=IErr)
  RSimulatedPatterns = ZERO
  
  !For Bloch wave calculation
  ALLOCATE(RDevPara(nReflections),STAT=IErr)
  ALLOCATE(IStrongBeamList(nReflections),STAT=IErr)
  ALLOCATE(IWeakBeamList(nReflections),STAT=IErr)
  ALLOCATE(CFullWaveFunctions(nReflections),STAT=IErr)
  ALLOCATE(RFullWaveIntensity(nReflections),STAT=IErr)
  
  !===========================================CORE SPECIFIC
  !Allocations for the pixels to be calculated by this core  
  ILocalPixelCountMin= (IPixelTotal*(my_rank)/p)+1
  ILocalPixelCountMax= (IPixelTotal*(my_rank+1)/p) 
  ALLOCATE(RIndividualReflections(INoOfLacbedPatterns,IThicknessCount,&
         (ILocalPixelCountMax-ILocalPixelCountMin)+1),STAT=IErr)
  ALLOCATE(IDisplacements(p),ICount(p),STAT=IErr)
  DO ind = 1,p
     IDisplacements(ind) = (IPixelTotal*(ind-1)/p)*INoOfLacbedPatterns*IThicknessCount
     ICount(ind) = (((IPixelTotal*(ind)/p) - (IPixelTotal*(ind-1)/p)))*INoOfLacbedPatterns*IThicknessCount    
  END DO
  !===========================================CORE SPECIFIC  
  
  IIterationCount = 0  
  
  CALL FelixFunction(LInitialSimulationFLAG,IErr)
  
  ALLOCATE(RWeightingCoefficients(INoOfLacbedPatterns),STAT=IErr)
  RWeightingCoefficients = ONE
   
  !>>>>>>>>>>Core 0 only
  IF(my_rank.EQ.0) THEN   
    IFinishFLAG=0
    IPreviousPrintedIteration = -IPrint!RB ensuring baseline simulation is printed
    CALL CalculateFigureofMeritandDetermineThickness(IThicknessIndex,IErr)
	CALL WriteIterationOutput(IIterationCount,IThicknessIndex,IFinishFLAG,IErr)
	CALL CreateRandomisedSimplex(RSimplexVariable,RIndependentVariable,IErr)
  END IF
  !>>>>>>>>>>>

  !+++++++++++++++++++++++++++++++++++++++++++MPI_BCAST
  CALL MPI_BCAST(RSimplexVariable,(INoOfVariables+1)*(INoOfVariables),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IErr)
  
  DO ind = 1,(INoOfVariables+1)
    IF(my_rank.EQ.0) THEN
      PRINT*,"--------------------------------"
      WRITE(SPrintString,FMT='(A8,I2,A4,I3)') "Simplex ",ind," of ",INoOfVariables+1
      PRINT*,TRIM(ADJUSTL(SPrintString))
      PRINT*,"--------------------------------"
    END IF
	CALL UpdateStructureFactors(RSimplexVariable(ind,:),IErr)
    IF(my_rank.EQ.0) THEN	
	  CALL PrintVariables(IErr)
    END IF	
    CALL FelixFunction(LInitialSimulationFLAG,IErr)	
    CALL CalculateFigureofMeritandDetermineThickness(IThicknessIndex,IErr)
	RSimplexFoM(ind)=RCrossCorrelation
	CALL WriteIterationOutput(IIterationCount,IThicknessIndex,IFinishFLAG,IErr)
	
  END DO
  
  
  
  
  
  
  IIterationCount = 1   
  
  
  
  
  
  
  
  
  
  
		 
  !**********************************  
   
  !--------------------------------------------------------------------
  ! tidy up
  DEALLOCATE(RFullWaveIntensity,STAT=IErr)
  DEALLOCATE(CFullWaveFunctions,STAT=IErr)
  DEALLOCATE(IWeakBeamList,STAT=IErr)
  DEALLOCATE(IStrongBeamList,STAT=IErr)
  DEALLOCATE(RDevPara,STAT=IErr)
  DEALLOCATE(RSimulatedPatterns,STAT=IErr)
  DEALLOCATE(RhklPositions,STAT=IErr)
  DEALLOCATE(RSimplexFoM,STAT=IErr) 
  DEALLOCATE(RSimplexVariable,STAT=IErr)  
  DEALLOCATE(RIndependentVariable,STAT=IErr)
  DEALLOCATE(CUgToRefine,STAT=IErr)  
  DEALLOCATE(IEquivalentUgKey,STAT=IErr)
  DEALLOCATE(ISymmetryRelations,STAT=IErr)  
  DEALLOCATE(RgSumMat,STAT=IErr)
  DEALLOCATE(CUgMat,STAT=IErr)
  DEALLOCATE(CUgMatNoAbs,STAT=IErr) 
  DEALLOCATE(RgVecVec,STAT=IErr)
  DEALLOCATE(RgPoolMag,STAT=IErr)
  DEALLOCATE(RgPoolT,STAT=IErr)
  DEALLOCATE(Rhkl,STAT=IErr)
  DEALLOCATE(RDWF,STAT=IErr)
  DEALLOCATE(ROcc,STAT=IErr)
  DEALLOCATE(RAtomCoordinate,STAT=IErr)
  DEALLOCATE(IAtoms,STAT=IErr)
  DEALLOCATE(RImageExpi,STAT=IErr)    
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
!  IF(IRefineModeSelectionArray(1).EQ.1) THEN
!     DO ind = 1,INoofUgs !RB ignore the first one as it is the internal potential
!        RIndependentVariable((ind-1)*2+1) = &!yy
!             REAL(CUgToRefine(ind+1),RKIND)!yy ind+1 instead of ind
!        RIndependentVariable((ind-1)*2+2) = &
!             AIMAG(CUgToRefine(ind+1))!yy ind+1 instead of ind
!     END DO
!  END IF
!  RIndependentVariable(2*INoofUgs+1) = RAbsorptionPercentage!RB absorption always included in structure factor refinement as last variable

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
     END SELECT
  END DO

END SUBROUTINE RefinementVariableSetup
 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!now redundant?
SUBROUTINE StructureFactorRefinementSetup(RIndependentVariable,IIterationCount,IErr)
  
  USE MyNumbers
  
  USE CConst; USE IConst; USE RConst
  USE IPara; USE RPara; USE SPara; USE CPara
  USE BlochPara
  
  USE IChannels
  
  USE MPI
  USE MyMPI
  
  IMPLICIT NONE
  
  INTEGER(IKIND) :: IErr,ind,jnd
  REAL(RKIND),DIMENSION(INoOfVariables),INTENT(OUT) :: RIndependentVariable
  INTEGER(IKIND),INTENT(IN) :: IIterationCount
  CHARACTER*200 :: SPrintString

  IF((IWriteFLAG.GE.5.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"StructureFactorRefinementSetup(",my_rank,")"
  END IF

  !jnd =0
  DO ind = 1,INoofUgs !RB ignore the first one as it is the internal potential
     RIndependentVariable((ind-1)*2+1) = &!yy
          REAL(CUgToRefine(ind+1),RKIND)!yy ind+1 instead of ind
     RIndependentVariable((ind-1)*2+2) = &
          AIMAG(CUgToRefine(ind+1))!yy ind+1 instead of ind
  END DO
  RIndependentVariable(2*INoofUgs+1) = RAbsorptionPercentage!RB absorption always included in structure factor refinement as last variable

END SUBROUTINE StructureFactorRefinementSetup

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE SetupUgsToRefine(IErr)
!Identify unique Ug's and count the number of independent variables INoOfVariables
!using the Hermitian matrix CUgMatNoAbs
!We count over INoofUgs, specified in felix.inp
!The count excludes Ug components that are zero and U(000), the inner potential
  
  USE MyNumbers
  
  USE CConst; USE IConst; USE RConst
  USE IPara; USE RPara; USE CPara
  USE BlochPara

  USE IChannels

  USE MPI
  USE MyMPI
  
  IMPLICIT NONE 
  
  INTEGER(IKIND) :: IErr,ind,jnd,Iuid
  INTEGER(IKIND),DIMENSION(2) :: ILoc
  CHARACTER*200 :: SPrintString

  IF((IWriteFLAG.GE.5.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"SetupUgsToRefine(",my_rank,")"
  END IF

  !Matrix with numbers marking equivalent Ug's
  ALLOCATE(ISymmetryRelations(nReflections,nReflections),STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"SetupUgsToRefine(",my_rank,")error allocating ISymmetryRelations"
     RETURN
  ENDIF

!Count equivalent Ugs
!Equivalent Ug's are identified by the sum of their abs(indices)plus the sum of abs(Ug)'s with no absorption
  RgSumMat = RgSumMat+ABS(REAL(CUgMatNoAbs))+ABS(AIMAG(CUgMatNoAbs))
  ISymmetryRelations = 0_IKIND 
  Iuid = 0_IKIND 
  DO ind = 1,nReflections
     DO jnd = 1,ind
        IF(ISymmetryRelations(ind,jnd).NE.0) THEN
           CYCLE
        ELSE
           Iuid = Iuid + 1_IKIND
           !Ug Fill the symmetry relation matrix with incrementing numbers that have the sign of the imaginary part
		   WHERE (ABS(RgSumMat-ABS(RgSumMat(ind,jnd))).LE.RTolerance)
              ISymmetryRelations = Iuid*SIGN(1_IKIND,NINT(AIMAG(CUgMatNoAbs)/TINY**2))
           END WHERE
        END IF
     END DO
  END DO

  IF((IWriteFLAG.GE.0.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     WRITE(SPrintString,FMT='(I5,A25)') Iuid," unique structure factors"
     PRINT*,TRIM(ADJUSTL(SPrintString))
  END IF
  IF(IWriteFLAG.EQ.3.AND.my_rank.EQ.0) THEN
    DO ind =1,8
     WRITE(SPrintString,FMT='(8(2X,F5.2,1X,F5.2))') CUgMatNoAbs(ind,1:8)
     PRINT*,TRIM(ADJUSTL(SPrintString))
    END DO
    DO ind =1,8
     WRITE(SPrintString,FMT='(3(1X,I3),8(1X,I3))') NINT(Rhkl(ind,:)),ISymmetryRelations(ind,1:8)
     PRINT*,TRIM(ADJUSTL(SPrintString))
    END DO
  END IF

!Link each number with its Ug, up to the number of unique Ug's Iuid
  ALLOCATE(IEquivalentUgKey(Iuid),STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"SetupUgsToRefine(",my_rank,")error allocating IEquivalentUgKey"
     RETURN
  END IF
  ALLOCATE(CUgToRefine(Iuid),STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"SetupUgsToRefine(",my_rank,")error allocating CUgToRefine"
     RETURN
  END IF
  
  DO ind = 1,Iuid
     ILoc = MINLOC(ABS(ISymmetryRelations-ind))
     IEquivalentUgKey(ind) = ind
     CUgToRefine(ind) = CUgMatNoAbs(ILoc(1),ILoc(2))
  END DO

   IF(IWriteFLAG.EQ.3.AND.my_rank.EQ.0) THEN
   PRINT*,"before sorting"
    DO ind =1,11
     WRITE(SPrintString,FMT='(I3,A1,F5.2,A1,2(1X,F5.2))')&
	 IEquivalentUgKey(ind),":",ABS(CUgToRefine(ind)),",",CUgToRefine(ind)
     PRINT*,TRIM(ADJUSTL(SPrintString))
    END DO
  END IF
  
!Put them in descending order of magnitude  
  CALL ReSortUgs(IEquivalentUgKey,CUgToRefine,Iuid)

   IF(IWriteFLAG.EQ.3.AND.my_rank.EQ.0) THEN
   PRINT*,"after sorting"
    DO ind =1,11
     WRITE(SPrintString,FMT='(I3,A1,F5.2,A1,2(1X,F5.2))')&
	 IEquivalentUgKey(ind),":",ABS(CUgToRefine(ind)),",",CUgToRefine(ind)
     PRINT*,TRIM(ADJUSTL(SPrintString))
    END DO
  END IF

!Count the number of Independent Variables
 
  jnd=1
  DO ind = 2,INoofUgs+1!RB ignore the first one as it is the internal potential
    IF ( ABS(REAL(CUgToRefine(ind),RKIND)).GE.RTolerance ) THEN
!      WRITE(SPrintString,FMT='(3(I3,A1),F5.2)') ind,",",jnd,":",IEquivalentUgKey(ind),":",ABS(REAL(CUgToRefine(ind),RKIND))
!      PRINT*,TRIM(ADJUSTL(SPrintString))
      jnd=jnd+1
	END IF
    IF ( ABS(AIMAG(CUgToRefine(ind))).GE.RTolerance ) THEN
!      WRITE(SPrintString,FMT='(3(I3,A1),F5.2)') ind,",",jnd,":",IEquivalentUgKey(ind),":",ABS(AIMAG(CUgToRefine(ind)))
!      PRINT*,TRIM(ADJUSTL(SPrintString))
      jnd=jnd+1
	END IF
  END DO
  
  INoOfVariables = jnd !RB last +1 is for absorption
  
END SUBROUTINE SetupUgsToRefine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE SimplexInitialisation(RSimplexVariable,RSimplexFoM,RIndependentVariable,&
     IIterationCount,RStandardDeviation,RMean,IErr)
  
  USE MyNumbers
  
  USE CConst; USE IConst; USE RConst
  USE IPara; USE RPara; USE SPara; USE CPara
  USE BlochPara
  
  USE IChannels
  
  USE MPI
  USE MyMPI
  
  IMPLICIT NONE

  INTEGER(IKIND) :: IErr,ind,jnd,IExitFLAG
  LOGICAL :: LInitialSimulationFLAG = .TRUE.
  REAL(RKIND),DIMENSION(INoOfVariables+1,INoOfVariables),INTENT(OUT) :: RSimplexVariable
  REAL(RKIND),DIMENSION(INoOfVariables+1),INTENT(OUT) :: RSimplexFoM
  REAL(RKIND) :: RSimplexDummy
  REAL(RKIND),DIMENSION(INoOfVariables),INTENT(INOUT) :: RIndependentVariable
  INTEGER(IKIND),INTENT(INOUT) :: IIterationCount
  REAL(RKIND),INTENT(OUT) :: RStandardDeviation,RMean
  REAL(RKIND) :: RStandardError,RStandardTolerance
  CHARACTER*200 :: SPrintString

  IF(IWriteFLAG.GE.10.AND.my_rank.EQ.0) THEN
     PRINT*,"SimplexInitialisation(",my_rank,")"
  END IF

  CALL FelixFunction(LInitialSimulationFLAG,IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"SimplexInitialisation(",my_rank,")error in FelixFunction"
     RETURN
  ENDIF

  CALL InitialiseWeightingCoefficients(IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"SimplexInitialisation(",my_rank,")error in InitialiseWeightingCoefficients"
     RETURN
  ENDIF

   IExitFLAG = 0; ! Do not exit
   IPreviousPrintedIteration = -IPrint!RB ensuring baseline simulation is printed
   IF(my_rank.EQ.0) THEN   
     CALL CreateImagesAndWriteOutput(IIterationCount,IExitFLAG,IErr) 
     IF( IErr.NE.0 ) THEN
       PRINT*,"SimplexInitialisation(",my_rank,")error in CreateImagesAndWriteOutput"
       RETURN
     ENDIF
  END IF

!  CALL RefinementVariableSetup(RIndependentVariable,IErr)
!  IF( IErr.NE.0 ) THEN
!     PRINT*,"SimplexInitialisation(", my_rank, ") error in RefinementVariableSetup()"
!     RETURN
!  ENDIF
 
!  IF(IRefineModeSelectionArray(1).EQ.1) THEN
!     CALL StructureFactorRefinementSetup(RIndependentVariable,IIterationCount,IErr)
!     IF( IErr.NE.0 ) THEN
!        PRINT*,"SimplexInitialisation(", my_rank, ") error in StructureFactorRefinementSetup()"
!        RETURN
!     ENDIF
!  ENDIF

!!$ RandomSequence
  !IF(IContinueFLAG.EQ.0) THEN
     IF(my_rank.EQ.0) THEN
        CALL CreateRandomisedSimplex(RSimplexVariable,RIndependentVariable,IErr)
     END IF
     CALL MPI_BCAST(RSimplexVariable,(INoOfVariables+1)*(INoOfVariables),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IErr)

     DO ind = 1,(INoOfVariables+1)
        IF((IWriteFLAG.GE.0.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
           PRINT*,"--------------------------------"
           WRITE(SPrintString,FMT='(A8,I2,A4,I3)') "Simplex ",ind," of ",INoOfVariables+1
           PRINT*,TRIM(ADJUSTL(SPrintString))
           PRINT*,"--------------------------------"
        END IF
        CALL SimplexFunction(RSimplexDummy,RSimplexVariable(ind,:),1,0,IErr)
        IF( IErr.NE.0 ) THEN
           PRINT*,"SimplexInitialisation(",my_rank,") error in SimplexFunction"
           RETURN
        ENDIF
        RStandardTolerance = RStandardError(RStandardDeviation,RMean,RSimplexDummy,IErr)
        RSimplexFoM(ind) =  RSimplexDummy
        IF((IWriteFLAG.GE.0.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
          WRITE(SPrintString,FMT='(A16,F7.5))') "Figure of merit ",RSimplexFoM(ind)
          PRINT*,TRIM(ADJUSTL(SPrintString))
        END IF
     END DO

  !ELSE
    
   !  CALL RecoverSavedSimplex(RSimplexVariable,RSimplexFoM,RStandardDeviation,RMean,IIterationCount,IErr)
   !  IF( IErr.NE.0 ) THEN
   !     PRINT*,"SimplexInitialisation (", my_rank, ") error in RecoverSavedSimplex()"
   !     RETURN
   !  ENDIF
     
  !END IF
     
END SUBROUTINE SimplexInitialisation

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
  
  IF(IRefineModeSelectionArray(2).EQ.1) THEN
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
             RRandomNumbers(INoOfVariables-IAllowedVectors),&
             STAT=IErr)
        
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

SUBROUTINE ApplyNewStructureFactors(IErr)

!!$  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!$  % Subroutine to place iteratively determined Structure factors
!!$  % to Ug Matrix
!!$  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  USE MyNumbers
  
  USE CConst; USE IConst; USE RConst
  USE IPara; USE RPara; USE SPara; USE CPara
  USE BlochPara

  USE IChannels

  USE MPI
  USE MyMPI

  IMPLICIT NONE

  INTEGER(IKIND) :: IErr,ind
  COMPLEX(CKIND),DIMENSION(nReflections,nReflections) :: CUgMatDummy

!!$  Dummy Matrix to contain new iterative values
  
   CUgMatDummy = CZERO

!!$  Populate Ug Matrix with new iterative elements, include proportional absorption here for now !RB not good
  DO ind = 1,INoofUgs
     WHERE(ISymmetryRelations.EQ.IEquivalentUgKey(ind))
        CUgMatDummy = CUgToRefine(ind)+&
		CUgToRefine(ind)*EXP(CIMAGONE*PI/2)*(RAbsorptionPercentage/100_RKIND)
     END WHERE
     WHERE(ISymmetryRelations.EQ.-IEquivalentUgKey(ind))!RB 
        CUgMatDummy = CONJG(CUgToRefine(ind))+&
		CONJG(CUgToRefine(ind))*EXP(CIMAGONE*PI/2)*(RAbsorptionPercentage/100_RKIND)
     END WHERE
  END DO

  WHERE(ABS(CUgMatDummy).GT.TINY)
     CUgMat = CUgMatDummy!RB
  END WHERE
  
END SUBROUTINE ApplyNewStructureFactors

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

SUBROUTINE RecoverSavedSimplex(RSimplexVariable,RSimplexFoM,RStandardDeviation,RMean,IIterationCount,IErr)

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
       IErr,ind,IIterationCount
  REAL(RKIND),DIMENSION(INoOfVariables+1,INoOfVariables) :: &
       RSimplexVariable
  REAL(RKIND),DIMENSION(INoOfVariables+1) :: &
       RSimplexFoM
  REAL(RKIND) :: &
       RStandardDeviation,RMean
  CHARACTER*200 :: &
       CSizeofData,SFormatString,filename

  WRITE(filename,*) "fr-Simplex.txt"

  OPEN(UNIT=IChOutSimplex,STATUS='UNKNOWN',&
        FILE=TRIM(ADJUSTL(filename)))
  
  WRITE(CSizeofData,*) INoOfVariables+1
  WRITE(SFormatString,*) "("//TRIM(ADJUSTL(CSizeofData))//"(1F6.3,1X),A1)"

  DO ind = 1,(INoOfVariables+1)
     READ(IChOutSimplex,FMT=SFormatString) RSimplexVariable(ind,:),RSimplexFoM(ind)
  END DO
    
  READ(IChOutSimplex,FMT="(2(1F6.3,1X),I5.1,I5.1,A1)") RStandardDeviation,RMean,IStandardDeviationCalls,IIterationCount

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
