!STARTOFREGISTRYGENERATEDFILE 'SCDataEx_Types.f90'
!
! WARNING This file is generated automatically by the FAST registry.
! Do not edit.  Your changes to this file will be lost.
!
! FAST Registry
!*********************************************************************************************************************************
! SCDataEx_Types
!.................................................................................................................................
! This file is part of SCDataEx.
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
!> This module contains the user-defined types needed in SCDataEx. It also contains copy, destroy, pack, and
!! unpack routines associated with each defined data type. This code is automatically generated by the FAST Registry.
MODULE SCDataEx_Types
!---------------------------------------------------------------------------------------------------------------------------------
USE NWTC_Library
IMPLICIT NONE
! =========  SC_DX_InitInputType_C  =======
  TYPE, BIND(C) :: SC_DX_InitInputType_C
   TYPE(C_PTR) :: object = C_NULL_PTR
    INTEGER(KIND=C_INT) :: NumSC2Ctrl 
    INTEGER(KIND=C_INT) :: NumSC2CtrlGlob 
    INTEGER(KIND=C_INT) :: NumCtrl2SC 
  END TYPE SC_DX_InitInputType_C
  TYPE, PUBLIC :: SC_DX_InitInputType
    TYPE( SC_DX_InitInputType_C ) :: C_obj
    INTEGER(IntKi)  :: NumSC2Ctrl = 0_IntKi      !< number of turbine specific controller inputs [from supercontroller] [-]
    INTEGER(IntKi)  :: NumSC2CtrlGlob = 0_IntKi      !< number of global controller inputs [from supercontroller] [-]
    INTEGER(IntKi)  :: NumCtrl2SC = 0_IntKi      !< number of controller outputs [to supercontroller] [-]
  END TYPE SC_DX_InitInputType
! =======================
! =========  SC_DX_InitOutputType_C  =======
  TYPE, BIND(C) :: SC_DX_InitOutputType_C
   TYPE(C_PTR) :: object = C_NULL_PTR
  END TYPE SC_DX_InitOutputType_C
  TYPE, PUBLIC :: SC_DX_InitOutputType
    TYPE( SC_DX_InitOutputType_C ) :: C_obj
    TYPE(ProgDesc)  :: Ver      !< This module's name, version, and date [-]
  END TYPE SC_DX_InitOutputType
! =======================
! =========  SC_DX_ParameterType_C  =======
  TYPE, BIND(C) :: SC_DX_ParameterType_C
   TYPE(C_PTR) :: object = C_NULL_PTR
    LOGICAL(KIND=C_BOOL) :: useSC 
  END TYPE SC_DX_ParameterType_C
  TYPE, PUBLIC :: SC_DX_ParameterType
    TYPE( SC_DX_ParameterType_C ) :: C_obj
    LOGICAL  :: useSC = .FALSE.      !< Flag that tells this module if supercontroller is on. [-]
  END TYPE SC_DX_ParameterType
! =======================
! =========  SC_DX_InputType_C  =======
  TYPE, BIND(C) :: SC_DX_InputType_C
   TYPE(C_PTR) :: object = C_NULL_PTR
    TYPE(C_ptr) :: toSC = C_NULL_PTR 
    INTEGER(C_int) :: toSC_Len = 0 
  END TYPE SC_DX_InputType_C
  TYPE, PUBLIC :: SC_DX_InputType
    TYPE( SC_DX_InputType_C ) :: C_obj
    REAL(KIND=C_FLOAT) , DIMENSION(:), POINTER  :: toSC => NULL()      !< inputs to the super controller (from the turbine controller) [-]
  END TYPE SC_DX_InputType
! =======================
! =========  SC_DX_OutputType_C  =======
  TYPE, BIND(C) :: SC_DX_OutputType_C
   TYPE(C_PTR) :: object = C_NULL_PTR
    TYPE(C_ptr) :: fromSC = C_NULL_PTR 
    INTEGER(C_int) :: fromSC_Len = 0 
    TYPE(C_ptr) :: fromSCglob = C_NULL_PTR 
    INTEGER(C_int) :: fromSCglob_Len = 0 
  END TYPE SC_DX_OutputType_C
  TYPE, PUBLIC :: SC_DX_OutputType
    TYPE( SC_DX_OutputType_C ) :: C_obj
    REAL(KIND=C_FLOAT) , DIMENSION(:), POINTER  :: fromSC => NULL()      !< global outputs of the super controller (to the turbine controller) [-]
    REAL(KIND=C_FLOAT) , DIMENSION(:), POINTER  :: fromSCglob => NULL()      !< turbine specific outputs of the super controller (to the turbine controller) [-]
  END TYPE SC_DX_OutputType
! =======================
CONTAINS

subroutine SC_DX_CopyInitInput(SrcInitInputData, DstInitInputData, CtrlCode, ErrStat, ErrMsg)
   type(SC_DX_InitInputType), intent(in) :: SrcInitInputData
   type(SC_DX_InitInputType), intent(inout) :: DstInitInputData
   integer(IntKi),  intent(in   ) :: CtrlCode
   integer(IntKi),  intent(  out) :: ErrStat
   character(*),    intent(  out) :: ErrMsg
   character(*), parameter        :: RoutineName = 'SC_DX_CopyInitInput'
   ErrStat = ErrID_None
   ErrMsg  = ''
   DstInitInputData%NumSC2Ctrl = SrcInitInputData%NumSC2Ctrl
   DstInitInputData%C_obj%NumSC2Ctrl = SrcInitInputData%C_obj%NumSC2Ctrl
   DstInitInputData%NumSC2CtrlGlob = SrcInitInputData%NumSC2CtrlGlob
   DstInitInputData%C_obj%NumSC2CtrlGlob = SrcInitInputData%C_obj%NumSC2CtrlGlob
   DstInitInputData%NumCtrl2SC = SrcInitInputData%NumCtrl2SC
   DstInitInputData%C_obj%NumCtrl2SC = SrcInitInputData%C_obj%NumCtrl2SC
end subroutine

subroutine SC_DX_DestroyInitInput(InitInputData, ErrStat, ErrMsg)
   type(SC_DX_InitInputType), intent(inout) :: InitInputData
   integer(IntKi),  intent(  out) :: ErrStat
   character(*),    intent(  out) :: ErrMsg
   character(*), parameter        :: RoutineName = 'SC_DX_DestroyInitInput'
   ErrStat = ErrID_None
   ErrMsg  = ''
end subroutine

subroutine SC_DX_PackInitInput(RF, Indata)
   type(RegFile), intent(inout) :: RF
   type(SC_DX_InitInputType), intent(in) :: InData
   character(*), parameter         :: RoutineName = 'SC_DX_PackInitInput'
   if (RF%ErrStat >= AbortErrLev) return
   if (c_associated(InData%C_obj%object)) then
      call SetErrStat(ErrID_Severe,'C_obj%object cannot be packed.', RF%ErrStat, RF%ErrMsg, RoutineName)
      return
   end if
   call RegPack(RF, InData%NumSC2Ctrl)
   call RegPack(RF, InData%NumSC2CtrlGlob)
   call RegPack(RF, InData%NumCtrl2SC)
   if (RegCheckErr(RF, RoutineName)) return
end subroutine

subroutine SC_DX_UnPackInitInput(RF, OutData)
   type(RegFile), intent(inout)    :: RF
   type(SC_DX_InitInputType), intent(inout) :: OutData
   character(*), parameter            :: RoutineName = 'SC_DX_UnPackInitInput'
   if (RF%ErrStat /= ErrID_None) return
   call RegUnpack(RF, OutData%NumSC2Ctrl); if (RegCheckErr(RF, RoutineName)) return
   OutData%C_obj%NumSC2Ctrl = OutData%NumSC2Ctrl
   call RegUnpack(RF, OutData%NumSC2CtrlGlob); if (RegCheckErr(RF, RoutineName)) return
   OutData%C_obj%NumSC2CtrlGlob = OutData%NumSC2CtrlGlob
   call RegUnpack(RF, OutData%NumCtrl2SC); if (RegCheckErr(RF, RoutineName)) return
   OutData%C_obj%NumCtrl2SC = OutData%NumCtrl2SC
end subroutine

SUBROUTINE SC_DX_C2Fary_CopyInitInput(InitInputData, ErrStat, ErrMsg, SkipPointers)
   TYPE(SC_DX_InitInputType), INTENT(INOUT) :: InitInputData
   INTEGER(IntKi),  INTENT(  OUT) :: ErrStat
   CHARACTER(*),    INTENT(  OUT) :: ErrMsg
   LOGICAL,OPTIONAL,INTENT(IN   ) :: SkipPointers
   ! 
   LOGICAL                        :: SkipPointers_local
   ErrStat = ErrID_None
   ErrMsg  = ""
   
   IF (PRESENT(SkipPointers)) THEN
      SkipPointers_local = SkipPointers
   ELSE
      SkipPointers_local = .false.
   END IF
   InitInputData%NumSC2Ctrl = InitInputData%C_obj%NumSC2Ctrl
   InitInputData%NumSC2CtrlGlob = InitInputData%C_obj%NumSC2CtrlGlob
   InitInputData%NumCtrl2SC = InitInputData%C_obj%NumCtrl2SC
END SUBROUTINE

SUBROUTINE SC_DX_F2C_CopyInitInput( InitInputData, ErrStat, ErrMsg, SkipPointers  )
   TYPE(SC_DX_InitInputType), INTENT(INOUT) :: InitInputData
   INTEGER(IntKi),  INTENT(  OUT) :: ErrStat
   CHARACTER(*),    INTENT(  OUT) :: ErrMsg
   LOGICAL,OPTIONAL,INTENT(IN   ) :: SkipPointers
   ! 
   LOGICAL                        :: SkipPointers_local
   ErrStat = ErrID_None
   ErrMsg  = ''
   
   IF (PRESENT(SkipPointers)) THEN
      SkipPointers_local = SkipPointers
   ELSE
      SkipPointers_local = .false.
   END IF
   InitInputData%C_obj%NumSC2Ctrl = InitInputData%NumSC2Ctrl
   InitInputData%C_obj%NumSC2CtrlGlob = InitInputData%NumSC2CtrlGlob
   InitInputData%C_obj%NumCtrl2SC = InitInputData%NumCtrl2SC
END SUBROUTINE

subroutine SC_DX_CopyInitOutput(SrcInitOutputData, DstInitOutputData, CtrlCode, ErrStat, ErrMsg)
   type(SC_DX_InitOutputType), intent(in) :: SrcInitOutputData
   type(SC_DX_InitOutputType), intent(inout) :: DstInitOutputData
   integer(IntKi),  intent(in   ) :: CtrlCode
   integer(IntKi),  intent(  out) :: ErrStat
   character(*),    intent(  out) :: ErrMsg
   integer(IntKi)                 :: ErrStat2
   character(ErrMsgLen)           :: ErrMsg2
   character(*), parameter        :: RoutineName = 'SC_DX_CopyInitOutput'
   ErrStat = ErrID_None
   ErrMsg  = ''
   call NWTC_Library_CopyProgDesc(SrcInitOutputData%Ver, DstInitOutputData%Ver, CtrlCode, ErrStat2, ErrMsg2)
   call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   if (ErrStat >= AbortErrLev) return
end subroutine

subroutine SC_DX_DestroyInitOutput(InitOutputData, ErrStat, ErrMsg)
   type(SC_DX_InitOutputType), intent(inout) :: InitOutputData
   integer(IntKi),  intent(  out) :: ErrStat
   character(*),    intent(  out) :: ErrMsg
   integer(IntKi)                 :: ErrStat2
   character(ErrMsgLen)           :: ErrMsg2
   character(*), parameter        :: RoutineName = 'SC_DX_DestroyInitOutput'
   ErrStat = ErrID_None
   ErrMsg  = ''
   call NWTC_Library_DestroyProgDesc(InitOutputData%Ver, ErrStat2, ErrMsg2)
   call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
end subroutine

subroutine SC_DX_PackInitOutput(RF, Indata)
   type(RegFile), intent(inout) :: RF
   type(SC_DX_InitOutputType), intent(in) :: InData
   character(*), parameter         :: RoutineName = 'SC_DX_PackInitOutput'
   if (RF%ErrStat >= AbortErrLev) return
   if (c_associated(InData%C_obj%object)) then
      call SetErrStat(ErrID_Severe,'C_obj%object cannot be packed.', RF%ErrStat, RF%ErrMsg, RoutineName)
      return
   end if
   call NWTC_Library_PackProgDesc(RF, InData%Ver) 
   if (RegCheckErr(RF, RoutineName)) return
end subroutine

subroutine SC_DX_UnPackInitOutput(RF, OutData)
   type(RegFile), intent(inout)    :: RF
   type(SC_DX_InitOutputType), intent(inout) :: OutData
   character(*), parameter            :: RoutineName = 'SC_DX_UnPackInitOutput'
   if (RF%ErrStat /= ErrID_None) return
   call NWTC_Library_UnpackProgDesc(RF, OutData%Ver) ! Ver 
end subroutine

SUBROUTINE SC_DX_C2Fary_CopyInitOutput(InitOutputData, ErrStat, ErrMsg, SkipPointers)
   TYPE(SC_DX_InitOutputType), INTENT(INOUT) :: InitOutputData
   INTEGER(IntKi),  INTENT(  OUT) :: ErrStat
   CHARACTER(*),    INTENT(  OUT) :: ErrMsg
   LOGICAL,OPTIONAL,INTENT(IN   ) :: SkipPointers
   ! 
   LOGICAL                        :: SkipPointers_local
   ErrStat = ErrID_None
   ErrMsg  = ""
   
   IF (PRESENT(SkipPointers)) THEN
      SkipPointers_local = SkipPointers
   ELSE
      SkipPointers_local = .false.
   END IF
END SUBROUTINE

SUBROUTINE SC_DX_F2C_CopyInitOutput( InitOutputData, ErrStat, ErrMsg, SkipPointers  )
   TYPE(SC_DX_InitOutputType), INTENT(INOUT) :: InitOutputData
   INTEGER(IntKi),  INTENT(  OUT) :: ErrStat
   CHARACTER(*),    INTENT(  OUT) :: ErrMsg
   LOGICAL,OPTIONAL,INTENT(IN   ) :: SkipPointers
   ! 
   LOGICAL                        :: SkipPointers_local
   ErrStat = ErrID_None
   ErrMsg  = ''
   
   IF (PRESENT(SkipPointers)) THEN
      SkipPointers_local = SkipPointers
   ELSE
      SkipPointers_local = .false.
   END IF
END SUBROUTINE

subroutine SC_DX_CopyParam(SrcParamData, DstParamData, CtrlCode, ErrStat, ErrMsg)
   type(SC_DX_ParameterType), intent(in) :: SrcParamData
   type(SC_DX_ParameterType), intent(inout) :: DstParamData
   integer(IntKi),  intent(in   ) :: CtrlCode
   integer(IntKi),  intent(  out) :: ErrStat
   character(*),    intent(  out) :: ErrMsg
   character(*), parameter        :: RoutineName = 'SC_DX_CopyParam'
   ErrStat = ErrID_None
   ErrMsg  = ''
   DstParamData%useSC = SrcParamData%useSC
   DstParamData%C_obj%useSC = SrcParamData%C_obj%useSC
end subroutine

subroutine SC_DX_DestroyParam(ParamData, ErrStat, ErrMsg)
   type(SC_DX_ParameterType), intent(inout) :: ParamData
   integer(IntKi),  intent(  out) :: ErrStat
   character(*),    intent(  out) :: ErrMsg
   character(*), parameter        :: RoutineName = 'SC_DX_DestroyParam'
   ErrStat = ErrID_None
   ErrMsg  = ''
end subroutine

subroutine SC_DX_PackParam(RF, Indata)
   type(RegFile), intent(inout) :: RF
   type(SC_DX_ParameterType), intent(in) :: InData
   character(*), parameter         :: RoutineName = 'SC_DX_PackParam'
   if (RF%ErrStat >= AbortErrLev) return
   if (c_associated(InData%C_obj%object)) then
      call SetErrStat(ErrID_Severe,'C_obj%object cannot be packed.', RF%ErrStat, RF%ErrMsg, RoutineName)
      return
   end if
   call RegPack(RF, InData%useSC)
   if (RegCheckErr(RF, RoutineName)) return
end subroutine

subroutine SC_DX_UnPackParam(RF, OutData)
   type(RegFile), intent(inout)    :: RF
   type(SC_DX_ParameterType), intent(inout) :: OutData
   character(*), parameter            :: RoutineName = 'SC_DX_UnPackParam'
   if (RF%ErrStat /= ErrID_None) return
   call RegUnpack(RF, OutData%useSC); if (RegCheckErr(RF, RoutineName)) return
   OutData%C_obj%useSC = OutData%useSC
end subroutine

SUBROUTINE SC_DX_C2Fary_CopyParam(ParamData, ErrStat, ErrMsg, SkipPointers)
   TYPE(SC_DX_ParameterType), INTENT(INOUT) :: ParamData
   INTEGER(IntKi),  INTENT(  OUT) :: ErrStat
   CHARACTER(*),    INTENT(  OUT) :: ErrMsg
   LOGICAL,OPTIONAL,INTENT(IN   ) :: SkipPointers
   ! 
   LOGICAL                        :: SkipPointers_local
   ErrStat = ErrID_None
   ErrMsg  = ""
   
   IF (PRESENT(SkipPointers)) THEN
      SkipPointers_local = SkipPointers
   ELSE
      SkipPointers_local = .false.
   END IF
   ParamData%useSC = ParamData%C_obj%useSC
END SUBROUTINE

SUBROUTINE SC_DX_F2C_CopyParam( ParamData, ErrStat, ErrMsg, SkipPointers  )
   TYPE(SC_DX_ParameterType), INTENT(INOUT) :: ParamData
   INTEGER(IntKi),  INTENT(  OUT) :: ErrStat
   CHARACTER(*),    INTENT(  OUT) :: ErrMsg
   LOGICAL,OPTIONAL,INTENT(IN   ) :: SkipPointers
   ! 
   LOGICAL                        :: SkipPointers_local
   ErrStat = ErrID_None
   ErrMsg  = ''
   
   IF (PRESENT(SkipPointers)) THEN
      SkipPointers_local = SkipPointers
   ELSE
      SkipPointers_local = .false.
   END IF
   ParamData%C_obj%useSC = ParamData%useSC
END SUBROUTINE

subroutine SC_DX_CopyInput(SrcInputData, DstInputData, CtrlCode, ErrStat, ErrMsg)
   type(SC_DX_InputType), intent(in) :: SrcInputData
   type(SC_DX_InputType), intent(inout) :: DstInputData
   integer(IntKi),  intent(in   ) :: CtrlCode
   integer(IntKi),  intent(  out) :: ErrStat
   character(*),    intent(  out) :: ErrMsg
   integer(B8Ki)                  :: LB(1), UB(1)
   integer(IntKi)                 :: ErrStat2
   character(*), parameter        :: RoutineName = 'SC_DX_CopyInput'
   ErrStat = ErrID_None
   ErrMsg  = ''
   if (associated(SrcInputData%toSC)) then
      LB(1:1) = lbound(SrcInputData%toSC, kind=B8Ki)
      UB(1:1) = ubound(SrcInputData%toSC, kind=B8Ki)
      if (.not. associated(DstInputData%toSC)) then
         allocate(DstInputData%toSC(LB(1):UB(1)), stat=ErrStat2)
         if (ErrStat2 /= 0) then
            call SetErrStat(ErrID_Fatal, 'Error allocating DstInputData%toSC.', ErrStat, ErrMsg, RoutineName)
            return
         end if
         DstInputData%C_obj%toSC_Len = size(DstInputData%toSC)
         if (DstInputData%C_obj%toSC_Len > 0) &
            DstInputData%C_obj%toSC = c_loc(DstInputData%toSC(LB(1)))
      end if
      DstInputData%toSC = SrcInputData%toSC
   end if
end subroutine

subroutine SC_DX_DestroyInput(InputData, ErrStat, ErrMsg)
   type(SC_DX_InputType), intent(inout) :: InputData
   integer(IntKi),  intent(  out) :: ErrStat
   character(*),    intent(  out) :: ErrMsg
   character(*), parameter        :: RoutineName = 'SC_DX_DestroyInput'
   ErrStat = ErrID_None
   ErrMsg  = ''
   if (associated(InputData%toSC)) then
      deallocate(InputData%toSC)
      InputData%toSC => null()
      InputData%C_obj%toSC = c_null_ptr
      InputData%C_obj%toSC_Len = 0
   end if
end subroutine

subroutine SC_DX_PackInput(RF, Indata)
   type(RegFile), intent(inout) :: RF
   type(SC_DX_InputType), intent(in) :: InData
   character(*), parameter         :: RoutineName = 'SC_DX_PackInput'
   logical         :: PtrInIndex
   if (RF%ErrStat >= AbortErrLev) return
   if (c_associated(InData%C_obj%object)) then
      call SetErrStat(ErrID_Severe,'C_obj%object cannot be packed.', RF%ErrStat, RF%ErrMsg, RoutineName)
      return
   end if
   call RegPackPtr(RF, InData%toSC)
   if (RegCheckErr(RF, RoutineName)) return
end subroutine

subroutine SC_DX_UnPackInput(RF, OutData)
   type(RegFile), intent(inout)    :: RF
   type(SC_DX_InputType), intent(inout) :: OutData
   character(*), parameter            :: RoutineName = 'SC_DX_UnPackInput'
   integer(B8Ki)   :: LB(1), UB(1)
   integer(IntKi)  :: stat
   logical         :: IsAllocAssoc
   integer(B8Ki)   :: PtrIdx
   type(c_ptr)     :: Ptr
   if (RF%ErrStat /= ErrID_None) return
   call RegUnpackPtr(RF, OutData%toSC); if (RegCheckErr(RF, RoutineName)) return
   if (associated(OutData%toSC)) then
      OutData%C_obj%toSC_Len = size(OutData%toSC)
      if (OutData%C_obj%toSC_Len > 0) OutData%C_obj%toSC = c_loc(OutData%toSC(LB(1)))
   end if
end subroutine

SUBROUTINE SC_DX_C2Fary_CopyInput(InputData, ErrStat, ErrMsg, SkipPointers)
   TYPE(SC_DX_InputType), INTENT(INOUT) :: InputData
   INTEGER(IntKi),  INTENT(  OUT) :: ErrStat
   CHARACTER(*),    INTENT(  OUT) :: ErrMsg
   LOGICAL,OPTIONAL,INTENT(IN   ) :: SkipPointers
   ! 
   LOGICAL                        :: SkipPointers_local
   ErrStat = ErrID_None
   ErrMsg  = ""
   
   IF (PRESENT(SkipPointers)) THEN
      SkipPointers_local = SkipPointers
   ELSE
      SkipPointers_local = .false.
   END IF
   
   ! -- toSC Input Data fields
   IF ( .NOT. SkipPointers_local ) THEN
      IF ( .NOT. C_ASSOCIATED( InputData%C_obj%toSC ) ) THEN
         NULLIFY( InputData%toSC )
      ELSE
         CALL C_F_POINTER(InputData%C_obj%toSC, InputData%toSC, [InputData%C_obj%toSC_Len])
      END IF
   END IF
END SUBROUTINE

SUBROUTINE SC_DX_F2C_CopyInput( InputData, ErrStat, ErrMsg, SkipPointers  )
   TYPE(SC_DX_InputType), INTENT(INOUT) :: InputData
   INTEGER(IntKi),  INTENT(  OUT) :: ErrStat
   CHARACTER(*),    INTENT(  OUT) :: ErrMsg
   LOGICAL,OPTIONAL,INTENT(IN   ) :: SkipPointers
   ! 
   LOGICAL                        :: SkipPointers_local
   ErrStat = ErrID_None
   ErrMsg  = ''
   
   IF (PRESENT(SkipPointers)) THEN
      SkipPointers_local = SkipPointers
   ELSE
      SkipPointers_local = .false.
   END IF
   
   ! -- toSC Input Data fields
   IF (.NOT. SkipPointers_local ) THEN
      IF (.NOT. ASSOCIATED(InputData%toSC)) THEN 
         InputData%C_obj%toSC_Len = 0
         InputData%C_obj%toSC = C_NULL_PTR
      ELSE
         InputData%C_obj%toSC_Len = SIZE(InputData%toSC)
         IF (InputData%C_obj%toSC_Len > 0) &
            InputData%C_obj%toSC = C_LOC(InputData%toSC(LBOUND(InputData%toSC,1, kind=B8Ki)))
      END IF
   END IF
END SUBROUTINE

subroutine SC_DX_CopyOutput(SrcOutputData, DstOutputData, CtrlCode, ErrStat, ErrMsg)
   type(SC_DX_OutputType), intent(in) :: SrcOutputData
   type(SC_DX_OutputType), intent(inout) :: DstOutputData
   integer(IntKi),  intent(in   ) :: CtrlCode
   integer(IntKi),  intent(  out) :: ErrStat
   character(*),    intent(  out) :: ErrMsg
   integer(B8Ki)                  :: LB(1), UB(1)
   integer(IntKi)                 :: ErrStat2
   character(*), parameter        :: RoutineName = 'SC_DX_CopyOutput'
   ErrStat = ErrID_None
   ErrMsg  = ''
   if (associated(SrcOutputData%fromSC)) then
      LB(1:1) = lbound(SrcOutputData%fromSC, kind=B8Ki)
      UB(1:1) = ubound(SrcOutputData%fromSC, kind=B8Ki)
      if (.not. associated(DstOutputData%fromSC)) then
         allocate(DstOutputData%fromSC(LB(1):UB(1)), stat=ErrStat2)
         if (ErrStat2 /= 0) then
            call SetErrStat(ErrID_Fatal, 'Error allocating DstOutputData%fromSC.', ErrStat, ErrMsg, RoutineName)
            return
         end if
         DstOutputData%C_obj%fromSC_Len = size(DstOutputData%fromSC)
         if (DstOutputData%C_obj%fromSC_Len > 0) &
            DstOutputData%C_obj%fromSC = c_loc(DstOutputData%fromSC(LB(1)))
      end if
      DstOutputData%fromSC = SrcOutputData%fromSC
   end if
   if (associated(SrcOutputData%fromSCglob)) then
      LB(1:1) = lbound(SrcOutputData%fromSCglob, kind=B8Ki)
      UB(1:1) = ubound(SrcOutputData%fromSCglob, kind=B8Ki)
      if (.not. associated(DstOutputData%fromSCglob)) then
         allocate(DstOutputData%fromSCglob(LB(1):UB(1)), stat=ErrStat2)
         if (ErrStat2 /= 0) then
            call SetErrStat(ErrID_Fatal, 'Error allocating DstOutputData%fromSCglob.', ErrStat, ErrMsg, RoutineName)
            return
         end if
         DstOutputData%C_obj%fromSCglob_Len = size(DstOutputData%fromSCglob)
         if (DstOutputData%C_obj%fromSCglob_Len > 0) &
            DstOutputData%C_obj%fromSCglob = c_loc(DstOutputData%fromSCglob(LB(1)))
      end if
      DstOutputData%fromSCglob = SrcOutputData%fromSCglob
   end if
end subroutine

subroutine SC_DX_DestroyOutput(OutputData, ErrStat, ErrMsg)
   type(SC_DX_OutputType), intent(inout) :: OutputData
   integer(IntKi),  intent(  out) :: ErrStat
   character(*),    intent(  out) :: ErrMsg
   character(*), parameter        :: RoutineName = 'SC_DX_DestroyOutput'
   ErrStat = ErrID_None
   ErrMsg  = ''
   if (associated(OutputData%fromSC)) then
      deallocate(OutputData%fromSC)
      OutputData%fromSC => null()
      OutputData%C_obj%fromSC = c_null_ptr
      OutputData%C_obj%fromSC_Len = 0
   end if
   if (associated(OutputData%fromSCglob)) then
      deallocate(OutputData%fromSCglob)
      OutputData%fromSCglob => null()
      OutputData%C_obj%fromSCglob = c_null_ptr
      OutputData%C_obj%fromSCglob_Len = 0
   end if
end subroutine

subroutine SC_DX_PackOutput(RF, Indata)
   type(RegFile), intent(inout) :: RF
   type(SC_DX_OutputType), intent(in) :: InData
   character(*), parameter         :: RoutineName = 'SC_DX_PackOutput'
   logical         :: PtrInIndex
   if (RF%ErrStat >= AbortErrLev) return
   if (c_associated(InData%C_obj%object)) then
      call SetErrStat(ErrID_Severe,'C_obj%object cannot be packed.', RF%ErrStat, RF%ErrMsg, RoutineName)
      return
   end if
   call RegPackPtr(RF, InData%fromSC)
   call RegPackPtr(RF, InData%fromSCglob)
   if (RegCheckErr(RF, RoutineName)) return
end subroutine

subroutine SC_DX_UnPackOutput(RF, OutData)
   type(RegFile), intent(inout)    :: RF
   type(SC_DX_OutputType), intent(inout) :: OutData
   character(*), parameter            :: RoutineName = 'SC_DX_UnPackOutput'
   integer(B8Ki)   :: LB(1), UB(1)
   integer(IntKi)  :: stat
   logical         :: IsAllocAssoc
   integer(B8Ki)   :: PtrIdx
   type(c_ptr)     :: Ptr
   if (RF%ErrStat /= ErrID_None) return
   call RegUnpackPtr(RF, OutData%fromSC); if (RegCheckErr(RF, RoutineName)) return
   if (associated(OutData%fromSC)) then
      OutData%C_obj%fromSC_Len = size(OutData%fromSC)
      if (OutData%C_obj%fromSC_Len > 0) OutData%C_obj%fromSC = c_loc(OutData%fromSC(LB(1)))
   end if
   call RegUnpackPtr(RF, OutData%fromSCglob); if (RegCheckErr(RF, RoutineName)) return
   if (associated(OutData%fromSCglob)) then
      OutData%C_obj%fromSCglob_Len = size(OutData%fromSCglob)
      if (OutData%C_obj%fromSCglob_Len > 0) OutData%C_obj%fromSCglob = c_loc(OutData%fromSCglob(LB(1)))
   end if
end subroutine

SUBROUTINE SC_DX_C2Fary_CopyOutput(OutputData, ErrStat, ErrMsg, SkipPointers)
   TYPE(SC_DX_OutputType), INTENT(INOUT) :: OutputData
   INTEGER(IntKi),  INTENT(  OUT) :: ErrStat
   CHARACTER(*),    INTENT(  OUT) :: ErrMsg
   LOGICAL,OPTIONAL,INTENT(IN   ) :: SkipPointers
   ! 
   LOGICAL                        :: SkipPointers_local
   ErrStat = ErrID_None
   ErrMsg  = ""
   
   IF (PRESENT(SkipPointers)) THEN
      SkipPointers_local = SkipPointers
   ELSE
      SkipPointers_local = .false.
   END IF
   
   ! -- fromSC Output Data fields
   IF ( .NOT. SkipPointers_local ) THEN
      IF ( .NOT. C_ASSOCIATED( OutputData%C_obj%fromSC ) ) THEN
         NULLIFY( OutputData%fromSC )
      ELSE
         CALL C_F_POINTER(OutputData%C_obj%fromSC, OutputData%fromSC, [OutputData%C_obj%fromSC_Len])
      END IF
   END IF
   
   ! -- fromSCglob Output Data fields
   IF ( .NOT. SkipPointers_local ) THEN
      IF ( .NOT. C_ASSOCIATED( OutputData%C_obj%fromSCglob ) ) THEN
         NULLIFY( OutputData%fromSCglob )
      ELSE
         CALL C_F_POINTER(OutputData%C_obj%fromSCglob, OutputData%fromSCglob, [OutputData%C_obj%fromSCglob_Len])
      END IF
   END IF
END SUBROUTINE

SUBROUTINE SC_DX_F2C_CopyOutput( OutputData, ErrStat, ErrMsg, SkipPointers  )
   TYPE(SC_DX_OutputType), INTENT(INOUT) :: OutputData
   INTEGER(IntKi),  INTENT(  OUT) :: ErrStat
   CHARACTER(*),    INTENT(  OUT) :: ErrMsg
   LOGICAL,OPTIONAL,INTENT(IN   ) :: SkipPointers
   ! 
   LOGICAL                        :: SkipPointers_local
   ErrStat = ErrID_None
   ErrMsg  = ''
   
   IF (PRESENT(SkipPointers)) THEN
      SkipPointers_local = SkipPointers
   ELSE
      SkipPointers_local = .false.
   END IF
   
   ! -- fromSC Output Data fields
   IF (.NOT. SkipPointers_local ) THEN
      IF (.NOT. ASSOCIATED(OutputData%fromSC)) THEN 
         OutputData%C_obj%fromSC_Len = 0
         OutputData%C_obj%fromSC = C_NULL_PTR
      ELSE
         OutputData%C_obj%fromSC_Len = SIZE(OutputData%fromSC)
         IF (OutputData%C_obj%fromSC_Len > 0) &
            OutputData%C_obj%fromSC = C_LOC(OutputData%fromSC(LBOUND(OutputData%fromSC,1, kind=B8Ki)))
      END IF
   END IF
   
   ! -- fromSCglob Output Data fields
   IF (.NOT. SkipPointers_local ) THEN
      IF (.NOT. ASSOCIATED(OutputData%fromSCglob)) THEN 
         OutputData%C_obj%fromSCglob_Len = 0
         OutputData%C_obj%fromSCglob = C_NULL_PTR
      ELSE
         OutputData%C_obj%fromSCglob_Len = SIZE(OutputData%fromSCglob)
         IF (OutputData%C_obj%fromSCglob_Len > 0) &
            OutputData%C_obj%fromSCglob = C_LOC(OutputData%fromSCglob(LBOUND(OutputData%fromSCglob,1, kind=B8Ki)))
      END IF
   END IF
END SUBROUTINE
END MODULE SCDataEx_Types
!ENDOFREGISTRYGENERATEDFILE
