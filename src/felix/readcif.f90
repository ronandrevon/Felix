!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Felix
!
! Richard Beanland, Keith Evans & Rudolf A Roemer
!
! (C) 2013-17, all rights reserved
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
!  Felix is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
!  
!  Felix is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!  
!  You should have received a copy of the GNU General Public License
!  along with Felix.  If not, see <http://www.gnu.org/licenses/>.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


!>
!! Procedure-description: Reads from felix.cif input file using toolbox
!!
!! Major-Authors: Keith Evans (2014), Richard Beanland (2016)
!!
SUBROUTINE ReadCif(IErr)

  ! -----------------------------------------------------------------------
  ! ReadCif: Read the input file
  !
  !	IErr	error code
  ! ----------------------------------------------------------------------
  !
  ! based on 
  ! 
  !                     CIF Tool Box Application 'tbx_ex.f'
  !                     ------------------------
  !                     Version: June 1998
  !
  ! ----------------------------------------------------------------------
  
  !
  ! The tool box common variable file 'ciftbx.cmn' must be present 
  ! at the start of EACH routine using the CIFtbx functions.

  USE MyNumbers
  
  USE CConst; USE IConst
  USE IPara; USE RPara; USE SPara
  USE CPara
  USE IChannels
  USE message_mod
  USE MyMPI
  
  IMPLICIT NONE

  INCLUDE       'ciftbx-f90.cmn'

  LOGICAL       f1,f2,f3
  CHARACTER*32  name
  CHARACTER*80  line,SPrintString
  CHARACTER*4   label(6)
  CHARACTER*1   SAlphabetarray(52)
  CHARACTER*52  alphabet
  CHARACTER*2   rs
  CHARACTER*1   slash
  CHARACTER string*(30)
  REAL(RKIND),DIMENSION(:),ALLOCATABLE :: RPrint
  REAL          cela,celb,celc,siga,sigb,sigc
  REAL          x,y,z,u,su,sx,sy,sz,B,sB,sOcc,Uso,suso,Occ
  REAL          numb,sdev,dum
  REAL          xf(6),yf(6),zf(6),uij(6,6)
  INTEGER       i,j,nsite, iset, imark
  DATA SAlphabetarray /"A","B","C","D","E","F","G","H","I","J","K","L","M","N", &
       "O","P","Q","R","S","T","U","V","W","X","Y","Z","a","b","c","d","e", &
       "f","g","h","i","j","k","l","m","n","o","p","q","r","s","t","u","v", &
       "w","x","y","z"/
  DATA alphabet /"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"/
  DATA          cela,celb,celc,siga,sigb,sigc/6*0./
  DATA          x,y,z,u,sx,sy,sz,su/8*0./
  DATA          xf,yf,zf,uij/54*0./
  DATA          rs/'\\'/

  INTEGER IAtomCount, ICommaPosLeft, ICommaPosRight, &
       Ipos,Idpos, IXYZminus,IFRACminus, Inum,Idenom,IAtomID
  CHARACTER*32 Csym(ITHREE)
  INTEGER IErr,ind,jnd

  ! fudge to deal with gfortran vs. g77
  slash = rs(1:1)

  ! Assign the CIFtbx files 
  f1 = init_( 1, 2, 3, 6 )

  ! Request dictionary validation check
  IF(.NOT.dict_('cif_core.dic','valid dtype')) THEN
    ! Restore the old clipping action for text fields
    clipt_ = .TRUE.
    pclipt_ = .TRUE.
  END IF
  ! Open the CIF to be accessed

  name='felix.cif'
  IF(.NOT.ocif_(name)) THEN
    IF(my_rank.EQ.0) THEN
      PRINT*,"Error:Cannot find .cif, exiting"
	  END IF
    IErr=1
    RETURN
  END IF

  ! Assign the data block to be accessed
  IF(.NOT.data_(' ')) THEN
    IF(my_rank.EQ.0) THEN
      PRINT*,"Error:No cif data_ statement found"
    END IF
    IErr=1
  END IF
 
  ! Extracts crystal forumla
  f1 = char_('_chemical_formula_structural', name)
  IF(.NOT.f1) THEN
     IF(my_rank.EQ.0) THEN
        PRINT*,"Error:ReadCif(", my_rank, ") Chemical formula missing!"
     END IF
     IErr=1
  END IF
  ChemicalFormula = StripChars(name) ! Strips spaces, brackets etc. from formula string

  ! Extract some cell dimensions; test all is OK
  ! NEED TO PUT IN A CHECK FOR LENGTH UNITS
  siga = 0.
  sigb = 0.
  sigc = 0.
  f1 = numb_('_cell_length_a', cela, siga)
  f2 = numb_('_cell_length_b', celb, sigb)
  f3 = numb_('_cell_length_c', celc, sigc)
  !error call
  IF(.NOT.(f1.AND.f2.AND.f3)) THEN
     IF(my_rank.EQ.0) THEN
        PRINT*,"Error:ReadCif(", my_rank, ") Cell dimension(s) missing!"
     END IF
     IErr=1
  END IF
  RLengthX=cela; RLengthY=celb; RLengthZ=celc

  siga = 0.
  sigb = 0.
  sigc = 0.
  f1 = numb_('_cell_angle_alpha', cela, siga)
  f2 = numb_('_cell_angle_beta', celb, sigb)
  f3 = numb_('_cell_angle_gamma', celc, sigc)
   IF(.NOT.(f1.AND.f2.AND.f3)) THEN
    IF(my_rank.EQ.0) THEN
      PRINT*,"Error:ReadCif(", my_rank, ") Cell angle(s) missing!"
    END IF
    IErr=1
  END IF

  ! convert angles from degrees to radians
  IF (cela.GT.TWOPI) THEN!assume this angle is expressed in degrees
    RAlpha=cela*DEG2RADIAN;
  END IF
  IF (celb.GT.TWOPI) THEN!assume this angle is expressed in degrees
    RBeta=celb*DEG2RADIAN;
  END IF
  IF (celc.GT.TWOPI) THEN!assume this angle is expressed in degrees
    RGamma=celc*DEG2RADIAN;
  END IF

  CALL message( LL, dbg14, "alpha,beta,gamma",(/ RAlpha,RBeta,RGamma /) )

  CALL message( LL, dbg14, "siga,sigb,sigc", Real( (/ siga,sigb, sigc /) ,RKIND) )


  f1 = numb_('_cell_volume', cela, siga)
  !Cell volume
  IF((f1) .EQV. .FALSE.) THEN
    IVolumeFLAG= 0
    RVolume= RLengthX*RLengthY*RLengthZ* &
          SQRT(1.0D0-COS(RAlpha)*COS(RAlpha)-COS(RBeta)*COS(RBeta)-COS(RGamma)*COS(RGamma) + &
          2.0D0*COS(RAlpha)*COS(RBeta)*COS(RGamma))
  ELSE 
    RVolume= cela
    IVolumeFLAG= 1
  END IF

  CALL message ( LL, "Unit cell volume", RVolume )

  DO      
    f1 = char_('_atom_type_symbol', name)
    CALL message (LL, dbg14, "_atom_type_symbol",f1 )
    IF(loop_.NEQV. .TRUE.) EXIT
  END DO

  ! Extract space group notation (expected char string)
  f1 = char_('_symmetry_space_group_name_H-M', name)
	
  !different types of space groups as well as different phrasing of Hall space groups
  IF (SCAN(name,'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz').EQ.0) THEN
    f1 = char_('_symmetry_space_group_name_Hall',name)
    IF (SCAN(name,'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz').EQ.0) THEN
      f1 = numb_('_symmetry_Int_tables_number',numb,sx)
      IF (numb.LT.TINY) THEN
        f1 = numb_('_space_group_IT_number',numb,sx)
        IF (numb.LT.TINY) THEN
          !Error message
          CALL message(LS,"No Space Group")
          IErr = 1
        ELSE
          name = CSpaceGrp(NINT(numb))
        END IF
      ELSE 
        name = CSpaceGrp(NINT(numb))
      END IF
    END IF
  END IF

  SSpaceGroupName=TRIM(name(1:1))
  SSpaceGrp = TRIM(ADJUSTL(name))
  
  !sometimes space group is input in lowercase letters - below changes the first letter to uppercase
  IF (SCAN(alphabet,SSpaceGroupName).GT.26) THEN
     SSpaceGroupName=SAlphabetarray(SCAN(alphabet,SSpaceGroupName)-26)
  END IF

  CALL message( LL, dbg14, "_symmetry_space_group_name_H-M ", SSpaceGrp)
  
  ! ----------------------------------------------------------
  ! Extract atom site data
  ! count how many atoms
  IAtomCount=0
  DO 
    f1 = char_('_atom_site_label', name)
    IAtomCount= IAtomCount+1
    IF(loop_ .NEQV. .TRUE.) EXIT
  END DO
  !check for consistency with IAtomicSitesToRefine
  IF (SIZE(IAtomicSitesToRefine,DIM=1).GT.IAtomCount) THEN
    CALL message("Number of atomic sites to refine is larger than the number of atoms")
    CALL message("Please correct in felix.inp")
    IErr=1
    RETURN
  END IF
  !allocate variables
  ALLOCATE(RBasisAtomPosition(IAtomCount,ITHREE),STAT=IErr)
  ALLOCATE(SBasisAtomLabel(IAtomCount),STAT=IErr)
  ALLOCATE(SBasisAtomName(IAtomCount),STAT=IErr)
  ALLOCATE(IBasisAtomicNumber(IAtomCount),STAT=IErr)
  ALLOCATE(RBasisIsoDW(IAtomCount),STAT=IErr)
  ALLOCATE(RBasisOccupancy(IAtomCount),STAT=IErr)
  ALLOCATE(RAnisotropicDebyeWallerFactorTensor(IAtomCount,ITHREE,ITHREE),STAT=IErr)
  ALLOCATE(IBasisAnisoDW(IAtomCount),STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Error:ReadCif(", my_rank, ") error ", IErr, " in ALLOCATE()"
     RETURN
  END IF

  !initialise variables
  IBasisAtomicNumber = 0
  RAnisotropicDebyeWallerFactorTensor = ZERO
  IMaxPossibleNAtomsUnitCell = 0
  ! input data loop
  DO ind=1,IAtomCount
    B = ZERO
    Uso = ZERO
    f1 = char_('_atom_site_label', name)
    SBasisAtomLabel(ind)=name
    f1 = char_('_atom_site_type_symbol', name)
    SBasisAtomName(ind)=name(1:2)
    ! remove the oxidation state numbers
    Ipos=SCAN(SBasisAtomName(ind),"1234567890")
    IF (Ipos.GT.0) WRITE(SBasisAtomName(ind),'(A1,A1)') name(1:1)," "
    !get atomic number
    DO jnd=1,INElements!NB must match SElementSymbolMatrix defined in smodules line 73
      IF(TRIM(SBasisAtomName(ind)).EQ.TRIM(SElementSymbolMatrix(jnd))) THEN
        IBasisAtomicNumber(ind)=jnd
      END IF
    END DO
    IF (IBasisAtomicNumber(ind).EQ.0.AND.my_rank.EQ.0) THEN
      WRITE(SPrintString,FMT='(A35,I3,A9,A5,1X,A2,A4,I3,A3,3F7.4,A1)')&
      "ReadCif: Could not find Z for atom ",ind,", symbol ",SBasisAtomName(ind)
      PRINT*,"Error:",TRIM(ADJUSTL(SPrintString))
      IErr=1
    END IF
    f2 = numb_('_atom_site_fract_x', x, sx)
    RBasisAtomPosition(ind,1)= x
    f2 = numb_('_atom_site_fract_y', y, sy)
    RBasisAtomPosition(ind,2)= y
    f2 = numb_('_atom_site_fract_z', z, sz)
    RBasisAtomPosition(ind,3)= z
    f2 = numb_('_atom_site_B_iso_or_equiv',B,sB)
    f2 = numb_('_atom_site_U_iso_or_equiv',Uso,suso)	
    IF(ABS(B).GT.TINY) THEN
      RBasisIsoDW(ind) = B
    ELSE
	  IF(ABS(Uso).GT.TINY) THEN
        RBasisIsoDW(ind) = Uso*(8*PI**2)
      ELSE
        RBasisIsoDW(ind) = RDebyeWallerConstant
      END IF
    END IF
    f2 = numb_('_atom_site_occupancy',Occ, sOcc)
    RBasisOccupancy(ind) = Occ

    !todo - a case where message subroutine struggles against classical print
    !IF(IWriteFLAG.EQ.7.AND.my_rank.EQ.0) THEN!optional output
    !  WRITE(SPrintString,FMT='(A4,I3,A2,A5,1X,A2,A3,I2,A3,3F7.4,A7,F5.3,A12,F7.4)')&
    !  "Atom",ind,": ",SBasisAtomLabel(ind),SBasisAtomName(ind)," Z=",IBasisAtomicNumber(ind),&
    !  ", [",RBasisAtomPosition(ind,:),"], DWF=",RBasisIsoDW(ind),", occupancy=",RBasisOccupancy(ind)
    !  PRINT*,TRIM(ADJUSTL(SPrintString))
    !END IF

    CALL message( LL, dbg7, "For Atom ",ind)
    CALL message( LL, dbg7, SBasisAtomLabel(ind)//SBasisAtomName(ind)//" Z=",IBasisAtomicNumber(ind) )
    CALL message( LL, dbg7, "RBasisAtomPosition", RBasisAtomPosition(ind,:) )
    CALL message( LL, dbg7, "(DWF, occupancy) respectively = ", (/ RBasisIsoDW(ind), RBasisOccupancy(ind) /) )
    
    IF(loop_ .NEQV. .TRUE.) EXIT
  END DO

  !Branch here depending on felixsim or felixrefine
  IF (ISoftwareMode.NE.0) THEN !felixrefine
    ALLOCATE(SWyckoffSymbols(SIZE(IAtomicSitesToRefine)),STAT=IErr)
    IF( IErr.NE.0 ) THEN
      PRINT*,"Error:Read Cif (refine)(", my_rank, ") error ", IErr, " in ALLOCATE()"
     RETURN
    END IF
    DO ind=1,IAtomCount
      f2 = char_('_atom_site_Wyckoff_symbol',name)
      WHERE (IAtomicSitesToRefine.EQ.ind)
	    SWyckoffSymbols = name
	  END WHERE
    END DO
  ELSE !felixsim
    ALLOCATE(SWyckoffSymbols(IAtomCount),STAT=IErr)
    IF( IErr.NE.0 ) THEN
      PRINT*,"Error:Read Cif (sim)(", my_rank, ") error ", IErr, " in ALLOCATE()"
      RETURN
    END IF
    DO ind=1,IAtomCount
      f2 = char_('_atom_site_Wyckoff_symbol',name)
	  SWyckoffSymbols(ind) = name
    END DO
  END IF
  
  DO ind=1,IAtomCount
    f2 = numb_('_atom_site_aniso_U_11',u,su) 
    RAnisotropicDebyeWallerFactorTensor(ind,1,1) = u
    f2 = numb_('_atom_site_aniso_U_22',u,su) 
    RAnisotropicDebyeWallerFactorTensor(ind,2,2) = u
    f2 = numb_('_atom_site_aniso_U_33',u,su) 
    RAnisotropicDebyeWallerFactorTensor(ind,3,3) = u
    f2 = numb_('_atom_site_aniso_U_23',u,su) 
    RAnisotropicDebyeWallerFactorTensor(ind,2,3) = u
    f2 = numb_('_atom_site_aniso_U_12',u,su) 
    RAnisotropicDebyeWallerFactorTensor(ind,1,2) = u
    f2 = numb_('_atom_site_aniso_U_13',u,su) 
    RAnisotropicDebyeWallerFactorTensor(ind,1,3) = u
    IBasisAnisoDW(ind) = ind
    IF(IAnisoDebyeWallerFactorFlag.EQ.1) THEN
      CALL message( LM,"RAnisotropicDebyeWallerFactorTensor, index = ",ind)
      CALL message( LM, "..",RAnisotropicDebyeWallerFactorTensor(ind,:,:) )
    END IF
  END DO

  ! counting loop
  ISymCount=0
  DO 
    f1 = char_('_symmetry_equiv_pos_as_xyz', name)
    DO 
      f2 = char_(name, line)
      ISymCount=ISymCount+1
      IF(text_ .NEQV. .TRUE.) EXIT
    END DO
    IF(loop_ .NEQV. .TRUE.) EXIT
  END DO

  ALLOCATE(RSymVec(ISymCount,ITHREE),STAT=IErr)
  ALLOCATE(RSymMat(ISymCount,ITHREE,ITHREE),STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Error:ReadCif(",my_rank,")error",IErr,"in allocations"
     RETURN
  END IF
  
  RSymVec=ZERO
  RSymMat=ZERO
  
  !data loop
  ISymCount=0
  DO 
    f1 = char_('_symmetry_equiv_pos_as_xyz', name)
    DO 
      f2 = char_(name, line)
      ISymCount=ISymCount+1
      ICommaPosLeft = SCAN(name, ",")
      ICommaPosRight= SCAN(name, ",",.TRUE.)
      Csym(1)= name(1:ICommaPosLeft-1)
      Csym(2)= name(ICommaPosLeft+1:ICommaPosRight-1)
      Csym(3)= name(ICommaPosRight+1:LEN_TRIM(name))
      DO ind=1,ITHREE
        IXYZminus=1
        IFRACminus=1
        name= Csym(ind)
        Ipos= SCAN(name, "xX")
        IF(Ipos > 0) THEN ! there is an X
          IF(Ipos>1) THEN
            IF(name(Ipos-1:Ipos-1)=="-") IXYZminus=-1
          END IF
          RSymMat(ISymCount,ind,1)=IXYZminus
        END IF
           
        Ipos= SCAN(name, "yY")
        IF(Ipos > 0) THEN ! there is a Y
          IF(Ipos>1) THEN
            IF(name(Ipos-1:Ipos-1)=="-") IXYZminus=-1
          END IF
          RSymMat(ISymCount, ind,2)=IXYZminus
        END IF
           
        Ipos= SCAN(name, "zZ")
        IF(Ipos > 0) THEN ! there is a Z
          IF(Ipos>1) THEN
            IF(name(Ipos-1:Ipos-1)=="-") IXYZminus=-1
          END IF
          RSymMat(ISymCount, ind,3)=IXYZminus
        END IF
           
        Ipos= SCAN(name, "/")
        IF(Ipos > 1) THEN
          IF(Ipos < LEN_TRIM(NAME) ) THEN ! there is an /
            Inum  = IACHAR(name(Ipos-1:Ipos-1))-48
            Idenom= IACHAR(name(Ipos+1:Ipos+1))-48
            IF(Ipos>2) THEN
              IF(name(Ipos-2:Ipos-2)=="-") IFRACminus=-1
            END IF
          END IF
          RSymVec(ISymCount,ind)=IFRACminus*REAL(Inum)/REAL(Idenom)
        END IF
      END DO
      IF(text_ .NEQV. .TRUE.) EXIT
    END DO

    IF(loop_ .NEQV. .TRUE.) EXIT
  END DO

  ! closes the cif file
  CALL close_

  ! Reset Message Counter
  IMessageCounter =0  
  
  RETURN

CONTAINS

  !>
  !! Procedure-description: Strips any character from string which is not simply
  !! numeric or alphabetic
  !!
  !! Major-Authors: Keith Evans (2014), Richard Beanland (2016)
  !! 
  CHARACTER(len(string)) FUNCTION StripChars(string)
    CHARACTER(*) :: string
    INTEGER :: i, n = 0
    DO i = 1,LEN_TRIM(string) 
      IF ( VERIFY(string(i:i), &
      "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890") &
      == 1 ) THEN
        CYCLE
      END IF 
      n = n+1 
      StripChars(n:n) = string(i:i) 
    END DO 
  END FUNCTION

END SUBROUTINE ReadCif