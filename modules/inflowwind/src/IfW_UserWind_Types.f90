!STARTOFREGISTRYGENERATEDFILE 'IfW_UserWind_Types.f90'
!
! WARNING This file is generated automatically by the FAST registry.
! Do not edit.  Your changes to this file will be lost.
!
! FAST Registry
!*********************************************************************************************************************************
! IfW_UserWind_Types
!.................................................................................................................................
! This file is part of IfW_UserWind.
!
! Copyright (C) 2012-2016 National Renewable Energy Laboratory
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.
!
!
! W A R N I N G : This file was automatically generated from the FAST registry.  Changes made to this file may be lost.
!
!*********************************************************************************************************************************
!> This module contains the user-defined types needed in IfW_UserWind. It also contains copy, destroy, pack, and
!! unpack routines associated with each defined data type. This code is automatically generated by the FAST Registry.
MODULE IfW_UserWind_Types
!---------------------------------------------------------------------------------------------------------------------------------
USE NWTC_Library
IMPLICIT NONE
! =========  IfW_UserWind_InitInputType  =======
  TYPE, PUBLIC :: IfW_UserWind_InitInputType
    CHARACTER(1024)  :: WindFileName      !< Name of the wind file to use [-]
  END TYPE IfW_UserWind_InitInputType
! =======================
! =========  IfW_UserWind_InitOutputType  =======
  TYPE, PUBLIC :: IfW_UserWind_InitOutputType
    TYPE(ProgDesc)  :: Ver      !< Version information off HHWind submodule [-]
  END TYPE IfW_UserWind_InitOutputType
! =======================
! =========  IfW_UserWind_MiscVarType  =======
  TYPE, PUBLIC :: IfW_UserWind_MiscVarType
    REAL(ReKi)  :: DummyMiscVar      !< Remove this variable if you have misc variables [-]
  END TYPE IfW_UserWind_MiscVarType
! =======================
! =========  IfW_UserWind_ParameterType  =======
  TYPE, PUBLIC :: IfW_UserWind_ParameterType
    REAL(SiKi)  :: dummy      !< remove if you have parameters [-]
  END TYPE IfW_UserWind_ParameterType
! =======================
CONTAINS
 SUBROUTINE IfW_UserWind_CopyInitInput( SrcInitInputData, DstInitInputData, CtrlCode, ErrStat, ErrMsg )
   TYPE(IfW_UserWind_InitInputType), INTENT(IN) :: SrcInitInputData
   TYPE(IfW_UserWind_InitInputType), INTENT(INOUT) :: DstInitInputData
   INTEGER(IntKi),  INTENT(IN   ) :: CtrlCode
   INTEGER(IntKi),  INTENT(  OUT) :: ErrStat
   CHARACTER(*),    INTENT(  OUT) :: ErrMsg
! Local 
   INTEGER(IntKi)                 :: i,j,k
   INTEGER(IntKi)                 :: ErrStat2
   CHARACTER(ErrMsgLen)           :: ErrMsg2
   CHARACTER(*), PARAMETER        :: RoutineName = 'IfW_UserWind_CopyInitInput'
! 
   ErrStat = ErrID_None
   ErrMsg  = ""
    DstInitInputData%WindFileName = SrcInitInputData%WindFileName
 END SUBROUTINE IfW_UserWind_CopyInitInput

 SUBROUTINE IfW_UserWind_DestroyInitInput( InitInputData, ErrStat, ErrMsg, DEALLOCATEpointers )
  TYPE(IfW_UserWind_InitInputType), INTENT(INOUT) :: InitInputData
  INTEGER(IntKi),  INTENT(  OUT) :: ErrStat
  CHARACTER(*),    INTENT(  OUT) :: ErrMsg
  LOGICAL,OPTIONAL,INTENT(IN   ) :: DEALLOCATEpointers
  
  INTEGER(IntKi)                 :: i, i1, i2, i3, i4, i5 
  LOGICAL                        :: DEALLOCATEpointers_local
  INTEGER(IntKi)                 :: ErrStat2
  CHARACTER(ErrMsgLen)           :: ErrMsg2
  CHARACTER(*),    PARAMETER :: RoutineName = 'IfW_UserWind_DestroyInitInput'

  ErrStat = ErrID_None
  ErrMsg  = ""

  IF (PRESENT(DEALLOCATEpointers)) THEN
     DEALLOCATEpointers_local = DEALLOCATEpointers
  ELSE
     DEALLOCATEpointers_local = .true.
  END IF
  
 END SUBROUTINE IfW_UserWind_DestroyInitInput

 SUBROUTINE IfW_UserWind_PackInitInput( ReKiBuf, DbKiBuf, IntKiBuf, Indata, ErrStat, ErrMsg, SizeOnly )
  REAL(ReKi),       ALLOCATABLE, INTENT(  OUT) :: ReKiBuf(:)
  REAL(DbKi),       ALLOCATABLE, INTENT(  OUT) :: DbKiBuf(:)
  INTEGER(IntKi),   ALLOCATABLE, INTENT(  OUT) :: IntKiBuf(:)
  TYPE(IfW_UserWind_InitInputType),  INTENT(IN) :: InData
  INTEGER(IntKi),   INTENT(  OUT) :: ErrStat
  CHARACTER(*),     INTENT(  OUT) :: ErrMsg
  LOGICAL,OPTIONAL, INTENT(IN   ) :: SizeOnly
    ! Local variables
  INTEGER(IntKi)                 :: Re_BufSz
  INTEGER(IntKi)                 :: Re_Xferred
  INTEGER(IntKi)                 :: Db_BufSz
  INTEGER(IntKi)                 :: Db_Xferred
  INTEGER(IntKi)                 :: Int_BufSz
  INTEGER(IntKi)                 :: Int_Xferred
  INTEGER(IntKi)                 :: i,i1,i2,i3,i4,i5
  LOGICAL                        :: OnlySize ! if present and true, do not pack, just allocate buffers
  INTEGER(IntKi)                 :: ErrStat2
  CHARACTER(ErrMsgLen)           :: ErrMsg2
  CHARACTER(*), PARAMETER        :: RoutineName = 'IfW_UserWind_PackInitInput'
 ! buffers to store subtypes, if any
  REAL(ReKi),      ALLOCATABLE   :: Re_Buf(:)
  REAL(DbKi),      ALLOCATABLE   :: Db_Buf(:)
  INTEGER(IntKi),  ALLOCATABLE   :: Int_Buf(:)

  OnlySize = .FALSE.
  IF ( PRESENT(SizeOnly) ) THEN
    OnlySize = SizeOnly
  ENDIF
    !
  ErrStat = ErrID_None
  ErrMsg  = ""
  Re_BufSz  = 0
  Db_BufSz  = 0
  Int_BufSz  = 0
      Int_BufSz  = Int_BufSz  + 1*LEN(InData%WindFileName)  ! WindFileName
  IF ( Re_BufSz  .GT. 0 ) THEN 
     ALLOCATE( ReKiBuf(  Re_BufSz  ), STAT=ErrStat2 )
     IF (ErrStat2 /= 0) THEN 
       CALL SetErrStat(ErrID_Fatal, 'Error allocating ReKiBuf.', ErrStat, ErrMsg,RoutineName)
       RETURN
     END IF
  END IF
  IF ( Db_BufSz  .GT. 0 ) THEN 
     ALLOCATE( DbKiBuf(  Db_BufSz  ), STAT=ErrStat2 )
     IF (ErrStat2 /= 0) THEN 
       CALL SetErrStat(ErrID_Fatal, 'Error allocating DbKiBuf.', ErrStat, ErrMsg,RoutineName)
       RETURN
     END IF
  END IF
  IF ( Int_BufSz  .GT. 0 ) THEN 
     ALLOCATE( IntKiBuf(  Int_BufSz  ), STAT=ErrStat2 )
     IF (ErrStat2 /= 0) THEN 
       CALL SetErrStat(ErrID_Fatal, 'Error allocating IntKiBuf.', ErrStat, ErrMsg,RoutineName)
       RETURN
     END IF
  END IF
  IF(OnlySize) RETURN ! return early if only trying to allocate buffers (not pack them)

  Re_Xferred  = 1
  Db_Xferred  = 1
  Int_Xferred = 1

    DO I = 1, LEN(InData%WindFileName)
      IntKiBuf(Int_Xferred) = ICHAR(InData%WindFileName(I:I), IntKi)
      Int_Xferred = Int_Xferred + 1
    END DO ! I
 END SUBROUTINE IfW_UserWind_PackInitInput

 SUBROUTINE IfW_UserWind_UnPackInitInput( ReKiBuf, DbKiBuf, IntKiBuf, Outdata, ErrStat, ErrMsg )
  REAL(ReKi),      ALLOCATABLE, INTENT(IN   ) :: ReKiBuf(:)
  REAL(DbKi),      ALLOCATABLE, INTENT(IN   ) :: DbKiBuf(:)
  INTEGER(IntKi),  ALLOCATABLE, INTENT(IN   ) :: IntKiBuf(:)
  TYPE(IfW_UserWind_InitInputType), INTENT(INOUT) :: OutData
  INTEGER(IntKi),  INTENT(  OUT) :: ErrStat
  CHARACTER(*),    INTENT(  OUT) :: ErrMsg
    ! Local variables
  INTEGER(IntKi)                 :: Buf_size
  INTEGER(IntKi)                 :: Re_Xferred
  INTEGER(IntKi)                 :: Db_Xferred
  INTEGER(IntKi)                 :: Int_Xferred
  INTEGER(IntKi)                 :: i
  INTEGER(IntKi)                 :: ErrStat2
  CHARACTER(ErrMsgLen)           :: ErrMsg2
  CHARACTER(*), PARAMETER        :: RoutineName = 'IfW_UserWind_UnPackInitInput'
 ! buffers to store meshes, if any
  REAL(ReKi),      ALLOCATABLE   :: Re_Buf(:)
  REAL(DbKi),      ALLOCATABLE   :: Db_Buf(:)
  INTEGER(IntKi),  ALLOCATABLE   :: Int_Buf(:)
    !
  ErrStat = ErrID_None
  ErrMsg  = ""
  Re_Xferred  = 1
  Db_Xferred  = 1
  Int_Xferred  = 1
    DO I = 1, LEN(OutData%WindFileName)
      OutData%WindFileName(I:I) = CHAR(IntKiBuf(Int_Xferred))
      Int_Xferred = Int_Xferred + 1
    END DO ! I
 END SUBROUTINE IfW_UserWind_UnPackInitInput

 SUBROUTINE IfW_UserWind_CopyInitOutput( SrcInitOutputData, DstInitOutputData, CtrlCode, ErrStat, ErrMsg )
   TYPE(IfW_UserWind_InitOutputType), INTENT(IN) :: SrcInitOutputData
   TYPE(IfW_UserWind_InitOutputType), INTENT(INOUT) :: DstInitOutputData
   INTEGER(IntKi),  INTENT(IN   ) :: CtrlCode
   INTEGER(IntKi),  INTENT(  OUT) :: ErrStat
   CHARACTER(*),    INTENT(  OUT) :: ErrMsg
! Local 
   INTEGER(IntKi)                 :: i,j,k
   INTEGER(IntKi)                 :: ErrStat2
   CHARACTER(ErrMsgLen)           :: ErrMsg2
   CHARACTER(*), PARAMETER        :: RoutineName = 'IfW_UserWind_CopyInitOutput'
! 
   ErrStat = ErrID_None
   ErrMsg  = ""
      CALL NWTC_Library_Copyprogdesc( SrcInitOutputData%Ver, DstInitOutputData%Ver, CtrlCode, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg,RoutineName)
         IF (ErrStat>=AbortErrLev) RETURN
 END SUBROUTINE IfW_UserWind_CopyInitOutput

 SUBROUTINE IfW_UserWind_DestroyInitOutput( InitOutputData, ErrStat, ErrMsg, DEALLOCATEpointers )
  TYPE(IfW_UserWind_InitOutputType), INTENT(INOUT) :: InitOutputData
  INTEGER(IntKi),  INTENT(  OUT) :: ErrStat
  CHARACTER(*),    INTENT(  OUT) :: ErrMsg
  LOGICAL,OPTIONAL,INTENT(IN   ) :: DEALLOCATEpointers
  
  INTEGER(IntKi)                 :: i, i1, i2, i3, i4, i5 
  LOGICAL                        :: DEALLOCATEpointers_local
  INTEGER(IntKi)                 :: ErrStat2
  CHARACTER(ErrMsgLen)           :: ErrMsg2
  CHARACTER(*),    PARAMETER :: RoutineName = 'IfW_UserWind_DestroyInitOutput'

  ErrStat = ErrID_None
  ErrMsg  = ""

  IF (PRESENT(DEALLOCATEpointers)) THEN
     DEALLOCATEpointers_local = DEALLOCATEpointers
  ELSE
     DEALLOCATEpointers_local = .true.
  END IF
  
  CALL NWTC_Library_Destroyprogdesc( InitOutputData%Ver, ErrStat2, ErrMsg2 )
     CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
 END SUBROUTINE IfW_UserWind_DestroyInitOutput

 SUBROUTINE IfW_UserWind_PackInitOutput( ReKiBuf, DbKiBuf, IntKiBuf, Indata, ErrStat, ErrMsg, SizeOnly )
  REAL(ReKi),       ALLOCATABLE, INTENT(  OUT) :: ReKiBuf(:)
  REAL(DbKi),       ALLOCATABLE, INTENT(  OUT) :: DbKiBuf(:)
  INTEGER(IntKi),   ALLOCATABLE, INTENT(  OUT) :: IntKiBuf(:)
  TYPE(IfW_UserWind_InitOutputType),  INTENT(IN) :: InData
  INTEGER(IntKi),   INTENT(  OUT) :: ErrStat
  CHARACTER(*),     INTENT(  OUT) :: ErrMsg
  LOGICAL,OPTIONAL, INTENT(IN   ) :: SizeOnly
    ! Local variables
  INTEGER(IntKi)                 :: Re_BufSz
  INTEGER(IntKi)                 :: Re_Xferred
  INTEGER(IntKi)                 :: Db_BufSz
  INTEGER(IntKi)                 :: Db_Xferred
  INTEGER(IntKi)                 :: Int_BufSz
  INTEGER(IntKi)                 :: Int_Xferred
  INTEGER(IntKi)                 :: i,i1,i2,i3,i4,i5
  LOGICAL                        :: OnlySize ! if present and true, do not pack, just allocate buffers
  INTEGER(IntKi)                 :: ErrStat2
  CHARACTER(ErrMsgLen)           :: ErrMsg2
  CHARACTER(*), PARAMETER        :: RoutineName = 'IfW_UserWind_PackInitOutput'
 ! buffers to store subtypes, if any
  REAL(ReKi),      ALLOCATABLE   :: Re_Buf(:)
  REAL(DbKi),      ALLOCATABLE   :: Db_Buf(:)
  INTEGER(IntKi),  ALLOCATABLE   :: Int_Buf(:)

  OnlySize = .FALSE.
  IF ( PRESENT(SizeOnly) ) THEN
    OnlySize = SizeOnly
  ENDIF
    !
  ErrStat = ErrID_None
  ErrMsg  = ""
  Re_BufSz  = 0
  Db_BufSz  = 0
  Int_BufSz  = 0
   ! Allocate buffers for subtypes, if any (we'll get sizes from these) 
      Int_BufSz   = Int_BufSz + 3  ! Ver: size of buffers for each call to pack subtype
      CALL NWTC_Library_Packprogdesc( Re_Buf, Db_Buf, Int_Buf, InData%Ver, ErrStat2, ErrMsg2, .TRUE. ) ! Ver 
        CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
        IF (ErrStat >= AbortErrLev) RETURN

      IF(ALLOCATED(Re_Buf)) THEN ! Ver
         Re_BufSz  = Re_BufSz  + SIZE( Re_Buf  )
         DEALLOCATE(Re_Buf)
      END IF
      IF(ALLOCATED(Db_Buf)) THEN ! Ver
         Db_BufSz  = Db_BufSz  + SIZE( Db_Buf  )
         DEALLOCATE(Db_Buf)
      END IF
      IF(ALLOCATED(Int_Buf)) THEN ! Ver
         Int_BufSz = Int_BufSz + SIZE( Int_Buf )
         DEALLOCATE(Int_Buf)
      END IF
  IF ( Re_BufSz  .GT. 0 ) THEN 
     ALLOCATE( ReKiBuf(  Re_BufSz  ), STAT=ErrStat2 )
     IF (ErrStat2 /= 0) THEN 
       CALL SetErrStat(ErrID_Fatal, 'Error allocating ReKiBuf.', ErrStat, ErrMsg,RoutineName)
       RETURN
     END IF
  END IF
  IF ( Db_BufSz  .GT. 0 ) THEN 
     ALLOCATE( DbKiBuf(  Db_BufSz  ), STAT=ErrStat2 )
     IF (ErrStat2 /= 0) THEN 
       CALL SetErrStat(ErrID_Fatal, 'Error allocating DbKiBuf.', ErrStat, ErrMsg,RoutineName)
       RETURN
     END IF
  END IF
  IF ( Int_BufSz  .GT. 0 ) THEN 
     ALLOCATE( IntKiBuf(  Int_BufSz  ), STAT=ErrStat2 )
     IF (ErrStat2 /= 0) THEN 
       CALL SetErrStat(ErrID_Fatal, 'Error allocating IntKiBuf.', ErrStat, ErrMsg,RoutineName)
       RETURN
     END IF
  END IF
  IF(OnlySize) RETURN ! return early if only trying to allocate buffers (not pack them)

  Re_Xferred  = 1
  Db_Xferred  = 1
  Int_Xferred = 1

      CALL NWTC_Library_Packprogdesc( Re_Buf, Db_Buf, Int_Buf, InData%Ver, ErrStat2, ErrMsg2, OnlySize ) ! Ver 
        CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
        IF (ErrStat >= AbortErrLev) RETURN

      IF(ALLOCATED(Re_Buf)) THEN
        IntKiBuf( Int_Xferred ) = SIZE(Re_Buf); Int_Xferred = Int_Xferred + 1
        IF (SIZE(Re_Buf) > 0) ReKiBuf( Re_Xferred:Re_Xferred+SIZE(Re_Buf)-1 ) = Re_Buf
        Re_Xferred = Re_Xferred + SIZE(Re_Buf)
        DEALLOCATE(Re_Buf)
      ELSE
        IntKiBuf( Int_Xferred ) = 0; Int_Xferred = Int_Xferred + 1
      ENDIF
      IF(ALLOCATED(Db_Buf)) THEN
        IntKiBuf( Int_Xferred ) = SIZE(Db_Buf); Int_Xferred = Int_Xferred + 1
        IF (SIZE(Db_Buf) > 0) DbKiBuf( Db_Xferred:Db_Xferred+SIZE(Db_Buf)-1 ) = Db_Buf
        Db_Xferred = Db_Xferred + SIZE(Db_Buf)
        DEALLOCATE(Db_Buf)
      ELSE
        IntKiBuf( Int_Xferred ) = 0; Int_Xferred = Int_Xferred + 1
      ENDIF
      IF(ALLOCATED(Int_Buf)) THEN
        IntKiBuf( Int_Xferred ) = SIZE(Int_Buf); Int_Xferred = Int_Xferred + 1
        IF (SIZE(Int_Buf) > 0) IntKiBuf( Int_Xferred:Int_Xferred+SIZE(Int_Buf)-1 ) = Int_Buf
        Int_Xferred = Int_Xferred + SIZE(Int_Buf)
        DEALLOCATE(Int_Buf)
      ELSE
        IntKiBuf( Int_Xferred ) = 0; Int_Xferred = Int_Xferred + 1
      ENDIF
 END SUBROUTINE IfW_UserWind_PackInitOutput

 SUBROUTINE IfW_UserWind_UnPackInitOutput( ReKiBuf, DbKiBuf, IntKiBuf, Outdata, ErrStat, ErrMsg )
  REAL(ReKi),      ALLOCATABLE, INTENT(IN   ) :: ReKiBuf(:)
  REAL(DbKi),      ALLOCATABLE, INTENT(IN   ) :: DbKiBuf(:)
  INTEGER(IntKi),  ALLOCATABLE, INTENT(IN   ) :: IntKiBuf(:)
  TYPE(IfW_UserWind_InitOutputType), INTENT(INOUT) :: OutData
  INTEGER(IntKi),  INTENT(  OUT) :: ErrStat
  CHARACTER(*),    INTENT(  OUT) :: ErrMsg
    ! Local variables
  INTEGER(IntKi)                 :: Buf_size
  INTEGER(IntKi)                 :: Re_Xferred
  INTEGER(IntKi)                 :: Db_Xferred
  INTEGER(IntKi)                 :: Int_Xferred
  INTEGER(IntKi)                 :: i
  INTEGER(IntKi)                 :: ErrStat2
  CHARACTER(ErrMsgLen)           :: ErrMsg2
  CHARACTER(*), PARAMETER        :: RoutineName = 'IfW_UserWind_UnPackInitOutput'
 ! buffers to store meshes, if any
  REAL(ReKi),      ALLOCATABLE   :: Re_Buf(:)
  REAL(DbKi),      ALLOCATABLE   :: Db_Buf(:)
  INTEGER(IntKi),  ALLOCATABLE   :: Int_Buf(:)
    !
  ErrStat = ErrID_None
  ErrMsg  = ""
  Re_Xferred  = 1
  Db_Xferred  = 1
  Int_Xferred  = 1
      Buf_size=IntKiBuf( Int_Xferred )
      Int_Xferred = Int_Xferred + 1
      IF(Buf_size > 0) THEN
        ALLOCATE(Re_Buf(Buf_size),STAT=ErrStat2)
        IF (ErrStat2 /= 0) THEN 
           CALL SetErrStat(ErrID_Fatal, 'Error allocating Re_Buf.', ErrStat, ErrMsg,RoutineName)
           RETURN
        END IF
        Re_Buf = ReKiBuf( Re_Xferred:Re_Xferred+Buf_size-1 )
        Re_Xferred = Re_Xferred + Buf_size
      END IF
      Buf_size=IntKiBuf( Int_Xferred )
      Int_Xferred = Int_Xferred + 1
      IF(Buf_size > 0) THEN
        ALLOCATE(Db_Buf(Buf_size),STAT=ErrStat2)
        IF (ErrStat2 /= 0) THEN 
           CALL SetErrStat(ErrID_Fatal, 'Error allocating Db_Buf.', ErrStat, ErrMsg,RoutineName)
           RETURN
        END IF
        Db_Buf = DbKiBuf( Db_Xferred:Db_Xferred+Buf_size-1 )
        Db_Xferred = Db_Xferred + Buf_size
      END IF
      Buf_size=IntKiBuf( Int_Xferred )
      Int_Xferred = Int_Xferred + 1
      IF(Buf_size > 0) THEN
        ALLOCATE(Int_Buf(Buf_size),STAT=ErrStat2)
        IF (ErrStat2 /= 0) THEN 
           CALL SetErrStat(ErrID_Fatal, 'Error allocating Int_Buf.', ErrStat, ErrMsg,RoutineName)
           RETURN
        END IF
        Int_Buf = IntKiBuf( Int_Xferred:Int_Xferred+Buf_size-1 )
        Int_Xferred = Int_Xferred + Buf_size
      END IF
      CALL NWTC_Library_Unpackprogdesc( Re_Buf, Db_Buf, Int_Buf, OutData%Ver, ErrStat2, ErrMsg2 ) ! Ver 
        CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
        IF (ErrStat >= AbortErrLev) RETURN

      IF(ALLOCATED(Re_Buf )) DEALLOCATE(Re_Buf )
      IF(ALLOCATED(Db_Buf )) DEALLOCATE(Db_Buf )
      IF(ALLOCATED(Int_Buf)) DEALLOCATE(Int_Buf)
 END SUBROUTINE IfW_UserWind_UnPackInitOutput

 SUBROUTINE IfW_UserWind_CopyMisc( SrcMiscData, DstMiscData, CtrlCode, ErrStat, ErrMsg )
   TYPE(IfW_UserWind_MiscVarType), INTENT(IN) :: SrcMiscData
   TYPE(IfW_UserWind_MiscVarType), INTENT(INOUT) :: DstMiscData
   INTEGER(IntKi),  INTENT(IN   ) :: CtrlCode
   INTEGER(IntKi),  INTENT(  OUT) :: ErrStat
   CHARACTER(*),    INTENT(  OUT) :: ErrMsg
! Local 
   INTEGER(IntKi)                 :: i,j,k
   INTEGER(IntKi)                 :: ErrStat2
   CHARACTER(ErrMsgLen)           :: ErrMsg2
   CHARACTER(*), PARAMETER        :: RoutineName = 'IfW_UserWind_CopyMisc'
! 
   ErrStat = ErrID_None
   ErrMsg  = ""
    DstMiscData%DummyMiscVar = SrcMiscData%DummyMiscVar
 END SUBROUTINE IfW_UserWind_CopyMisc

 SUBROUTINE IfW_UserWind_DestroyMisc( MiscData, ErrStat, ErrMsg, DEALLOCATEpointers )
  TYPE(IfW_UserWind_MiscVarType), INTENT(INOUT) :: MiscData
  INTEGER(IntKi),  INTENT(  OUT) :: ErrStat
  CHARACTER(*),    INTENT(  OUT) :: ErrMsg
  LOGICAL,OPTIONAL,INTENT(IN   ) :: DEALLOCATEpointers
  
  INTEGER(IntKi)                 :: i, i1, i2, i3, i4, i5 
  LOGICAL                        :: DEALLOCATEpointers_local
  INTEGER(IntKi)                 :: ErrStat2
  CHARACTER(ErrMsgLen)           :: ErrMsg2
  CHARACTER(*),    PARAMETER :: RoutineName = 'IfW_UserWind_DestroyMisc'

  ErrStat = ErrID_None
  ErrMsg  = ""

  IF (PRESENT(DEALLOCATEpointers)) THEN
     DEALLOCATEpointers_local = DEALLOCATEpointers
  ELSE
     DEALLOCATEpointers_local = .true.
  END IF
  
 END SUBROUTINE IfW_UserWind_DestroyMisc

 SUBROUTINE IfW_UserWind_PackMisc( ReKiBuf, DbKiBuf, IntKiBuf, Indata, ErrStat, ErrMsg, SizeOnly )
  REAL(ReKi),       ALLOCATABLE, INTENT(  OUT) :: ReKiBuf(:)
  REAL(DbKi),       ALLOCATABLE, INTENT(  OUT) :: DbKiBuf(:)
  INTEGER(IntKi),   ALLOCATABLE, INTENT(  OUT) :: IntKiBuf(:)
  TYPE(IfW_UserWind_MiscVarType),  INTENT(IN) :: InData
  INTEGER(IntKi),   INTENT(  OUT) :: ErrStat
  CHARACTER(*),     INTENT(  OUT) :: ErrMsg
  LOGICAL,OPTIONAL, INTENT(IN   ) :: SizeOnly
    ! Local variables
  INTEGER(IntKi)                 :: Re_BufSz
  INTEGER(IntKi)                 :: Re_Xferred
  INTEGER(IntKi)                 :: Db_BufSz
  INTEGER(IntKi)                 :: Db_Xferred
  INTEGER(IntKi)                 :: Int_BufSz
  INTEGER(IntKi)                 :: Int_Xferred
  INTEGER(IntKi)                 :: i,i1,i2,i3,i4,i5
  LOGICAL                        :: OnlySize ! if present and true, do not pack, just allocate buffers
  INTEGER(IntKi)                 :: ErrStat2
  CHARACTER(ErrMsgLen)           :: ErrMsg2
  CHARACTER(*), PARAMETER        :: RoutineName = 'IfW_UserWind_PackMisc'
 ! buffers to store subtypes, if any
  REAL(ReKi),      ALLOCATABLE   :: Re_Buf(:)
  REAL(DbKi),      ALLOCATABLE   :: Db_Buf(:)
  INTEGER(IntKi),  ALLOCATABLE   :: Int_Buf(:)

  OnlySize = .FALSE.
  IF ( PRESENT(SizeOnly) ) THEN
    OnlySize = SizeOnly
  ENDIF
    !
  ErrStat = ErrID_None
  ErrMsg  = ""
  Re_BufSz  = 0
  Db_BufSz  = 0
  Int_BufSz  = 0
      Re_BufSz   = Re_BufSz   + 1  ! DummyMiscVar
  IF ( Re_BufSz  .GT. 0 ) THEN 
     ALLOCATE( ReKiBuf(  Re_BufSz  ), STAT=ErrStat2 )
     IF (ErrStat2 /= 0) THEN 
       CALL SetErrStat(ErrID_Fatal, 'Error allocating ReKiBuf.', ErrStat, ErrMsg,RoutineName)
       RETURN
     END IF
  END IF
  IF ( Db_BufSz  .GT. 0 ) THEN 
     ALLOCATE( DbKiBuf(  Db_BufSz  ), STAT=ErrStat2 )
     IF (ErrStat2 /= 0) THEN 
       CALL SetErrStat(ErrID_Fatal, 'Error allocating DbKiBuf.', ErrStat, ErrMsg,RoutineName)
       RETURN
     END IF
  END IF
  IF ( Int_BufSz  .GT. 0 ) THEN 
     ALLOCATE( IntKiBuf(  Int_BufSz  ), STAT=ErrStat2 )
     IF (ErrStat2 /= 0) THEN 
       CALL SetErrStat(ErrID_Fatal, 'Error allocating IntKiBuf.', ErrStat, ErrMsg,RoutineName)
       RETURN
     END IF
  END IF
  IF(OnlySize) RETURN ! return early if only trying to allocate buffers (not pack them)

  Re_Xferred  = 1
  Db_Xferred  = 1
  Int_Xferred = 1

    ReKiBuf(Re_Xferred) = InData%DummyMiscVar
    Re_Xferred = Re_Xferred + 1
 END SUBROUTINE IfW_UserWind_PackMisc

 SUBROUTINE IfW_UserWind_UnPackMisc( ReKiBuf, DbKiBuf, IntKiBuf, Outdata, ErrStat, ErrMsg )
  REAL(ReKi),      ALLOCATABLE, INTENT(IN   ) :: ReKiBuf(:)
  REAL(DbKi),      ALLOCATABLE, INTENT(IN   ) :: DbKiBuf(:)
  INTEGER(IntKi),  ALLOCATABLE, INTENT(IN   ) :: IntKiBuf(:)
  TYPE(IfW_UserWind_MiscVarType), INTENT(INOUT) :: OutData
  INTEGER(IntKi),  INTENT(  OUT) :: ErrStat
  CHARACTER(*),    INTENT(  OUT) :: ErrMsg
    ! Local variables
  INTEGER(IntKi)                 :: Buf_size
  INTEGER(IntKi)                 :: Re_Xferred
  INTEGER(IntKi)                 :: Db_Xferred
  INTEGER(IntKi)                 :: Int_Xferred
  INTEGER(IntKi)                 :: i
  INTEGER(IntKi)                 :: ErrStat2
  CHARACTER(ErrMsgLen)           :: ErrMsg2
  CHARACTER(*), PARAMETER        :: RoutineName = 'IfW_UserWind_UnPackMisc'
 ! buffers to store meshes, if any
  REAL(ReKi),      ALLOCATABLE   :: Re_Buf(:)
  REAL(DbKi),      ALLOCATABLE   :: Db_Buf(:)
  INTEGER(IntKi),  ALLOCATABLE   :: Int_Buf(:)
    !
  ErrStat = ErrID_None
  ErrMsg  = ""
  Re_Xferred  = 1
  Db_Xferred  = 1
  Int_Xferred  = 1
    OutData%DummyMiscVar = ReKiBuf(Re_Xferred)
    Re_Xferred = Re_Xferred + 1
 END SUBROUTINE IfW_UserWind_UnPackMisc

 SUBROUTINE IfW_UserWind_CopyParam( SrcParamData, DstParamData, CtrlCode, ErrStat, ErrMsg )
   TYPE(IfW_UserWind_ParameterType), INTENT(IN) :: SrcParamData
   TYPE(IfW_UserWind_ParameterType), INTENT(INOUT) :: DstParamData
   INTEGER(IntKi),  INTENT(IN   ) :: CtrlCode
   INTEGER(IntKi),  INTENT(  OUT) :: ErrStat
   CHARACTER(*),    INTENT(  OUT) :: ErrMsg
! Local 
   INTEGER(IntKi)                 :: i,j,k
   INTEGER(IntKi)                 :: ErrStat2
   CHARACTER(ErrMsgLen)           :: ErrMsg2
   CHARACTER(*), PARAMETER        :: RoutineName = 'IfW_UserWind_CopyParam'
! 
   ErrStat = ErrID_None
   ErrMsg  = ""
    DstParamData%dummy = SrcParamData%dummy
 END SUBROUTINE IfW_UserWind_CopyParam

 SUBROUTINE IfW_UserWind_DestroyParam( ParamData, ErrStat, ErrMsg, DEALLOCATEpointers )
  TYPE(IfW_UserWind_ParameterType), INTENT(INOUT) :: ParamData
  INTEGER(IntKi),  INTENT(  OUT) :: ErrStat
  CHARACTER(*),    INTENT(  OUT) :: ErrMsg
  LOGICAL,OPTIONAL,INTENT(IN   ) :: DEALLOCATEpointers
  
  INTEGER(IntKi)                 :: i, i1, i2, i3, i4, i5 
  LOGICAL                        :: DEALLOCATEpointers_local
  INTEGER(IntKi)                 :: ErrStat2
  CHARACTER(ErrMsgLen)           :: ErrMsg2
  CHARACTER(*),    PARAMETER :: RoutineName = 'IfW_UserWind_DestroyParam'

  ErrStat = ErrID_None
  ErrMsg  = ""

  IF (PRESENT(DEALLOCATEpointers)) THEN
     DEALLOCATEpointers_local = DEALLOCATEpointers
  ELSE
     DEALLOCATEpointers_local = .true.
  END IF
  
 END SUBROUTINE IfW_UserWind_DestroyParam

 SUBROUTINE IfW_UserWind_PackParam( ReKiBuf, DbKiBuf, IntKiBuf, Indata, ErrStat, ErrMsg, SizeOnly )
  REAL(ReKi),       ALLOCATABLE, INTENT(  OUT) :: ReKiBuf(:)
  REAL(DbKi),       ALLOCATABLE, INTENT(  OUT) :: DbKiBuf(:)
  INTEGER(IntKi),   ALLOCATABLE, INTENT(  OUT) :: IntKiBuf(:)
  TYPE(IfW_UserWind_ParameterType),  INTENT(IN) :: InData
  INTEGER(IntKi),   INTENT(  OUT) :: ErrStat
  CHARACTER(*),     INTENT(  OUT) :: ErrMsg
  LOGICAL,OPTIONAL, INTENT(IN   ) :: SizeOnly
    ! Local variables
  INTEGER(IntKi)                 :: Re_BufSz
  INTEGER(IntKi)                 :: Re_Xferred
  INTEGER(IntKi)                 :: Db_BufSz
  INTEGER(IntKi)                 :: Db_Xferred
  INTEGER(IntKi)                 :: Int_BufSz
  INTEGER(IntKi)                 :: Int_Xferred
  INTEGER(IntKi)                 :: i,i1,i2,i3,i4,i5
  LOGICAL                        :: OnlySize ! if present and true, do not pack, just allocate buffers
  INTEGER(IntKi)                 :: ErrStat2
  CHARACTER(ErrMsgLen)           :: ErrMsg2
  CHARACTER(*), PARAMETER        :: RoutineName = 'IfW_UserWind_PackParam'
 ! buffers to store subtypes, if any
  REAL(ReKi),      ALLOCATABLE   :: Re_Buf(:)
  REAL(DbKi),      ALLOCATABLE   :: Db_Buf(:)
  INTEGER(IntKi),  ALLOCATABLE   :: Int_Buf(:)

  OnlySize = .FALSE.
  IF ( PRESENT(SizeOnly) ) THEN
    OnlySize = SizeOnly
  ENDIF
    !
  ErrStat = ErrID_None
  ErrMsg  = ""
  Re_BufSz  = 0
  Db_BufSz  = 0
  Int_BufSz  = 0
      Re_BufSz   = Re_BufSz   + 1  ! dummy
  IF ( Re_BufSz  .GT. 0 ) THEN 
     ALLOCATE( ReKiBuf(  Re_BufSz  ), STAT=ErrStat2 )
     IF (ErrStat2 /= 0) THEN 
       CALL SetErrStat(ErrID_Fatal, 'Error allocating ReKiBuf.', ErrStat, ErrMsg,RoutineName)
       RETURN
     END IF
  END IF
  IF ( Db_BufSz  .GT. 0 ) THEN 
     ALLOCATE( DbKiBuf(  Db_BufSz  ), STAT=ErrStat2 )
     IF (ErrStat2 /= 0) THEN 
       CALL SetErrStat(ErrID_Fatal, 'Error allocating DbKiBuf.', ErrStat, ErrMsg,RoutineName)
       RETURN
     END IF
  END IF
  IF ( Int_BufSz  .GT. 0 ) THEN 
     ALLOCATE( IntKiBuf(  Int_BufSz  ), STAT=ErrStat2 )
     IF (ErrStat2 /= 0) THEN 
       CALL SetErrStat(ErrID_Fatal, 'Error allocating IntKiBuf.', ErrStat, ErrMsg,RoutineName)
       RETURN
     END IF
  END IF
  IF(OnlySize) RETURN ! return early if only trying to allocate buffers (not pack them)

  Re_Xferred  = 1
  Db_Xferred  = 1
  Int_Xferred = 1

    ReKiBuf(Re_Xferred) = InData%dummy
    Re_Xferred = Re_Xferred + 1
 END SUBROUTINE IfW_UserWind_PackParam

 SUBROUTINE IfW_UserWind_UnPackParam( ReKiBuf, DbKiBuf, IntKiBuf, Outdata, ErrStat, ErrMsg )
  REAL(ReKi),      ALLOCATABLE, INTENT(IN   ) :: ReKiBuf(:)
  REAL(DbKi),      ALLOCATABLE, INTENT(IN   ) :: DbKiBuf(:)
  INTEGER(IntKi),  ALLOCATABLE, INTENT(IN   ) :: IntKiBuf(:)
  TYPE(IfW_UserWind_ParameterType), INTENT(INOUT) :: OutData
  INTEGER(IntKi),  INTENT(  OUT) :: ErrStat
  CHARACTER(*),    INTENT(  OUT) :: ErrMsg
    ! Local variables
  INTEGER(IntKi)                 :: Buf_size
  INTEGER(IntKi)                 :: Re_Xferred
  INTEGER(IntKi)                 :: Db_Xferred
  INTEGER(IntKi)                 :: Int_Xferred
  INTEGER(IntKi)                 :: i
  INTEGER(IntKi)                 :: ErrStat2
  CHARACTER(ErrMsgLen)           :: ErrMsg2
  CHARACTER(*), PARAMETER        :: RoutineName = 'IfW_UserWind_UnPackParam'
 ! buffers to store meshes, if any
  REAL(ReKi),      ALLOCATABLE   :: Re_Buf(:)
  REAL(DbKi),      ALLOCATABLE   :: Db_Buf(:)
  INTEGER(IntKi),  ALLOCATABLE   :: Int_Buf(:)
    !
  ErrStat = ErrID_None
  ErrMsg  = ""
  Re_Xferred  = 1
  Db_Xferred  = 1
  Int_Xferred  = 1
    OutData%dummy = REAL(ReKiBuf(Re_Xferred), SiKi)
    Re_Xferred = Re_Xferred + 1
 END SUBROUTINE IfW_UserWind_UnPackParam

END MODULE IfW_UserWind_Types
!ENDOFREGISTRYGENERATEDFILE
