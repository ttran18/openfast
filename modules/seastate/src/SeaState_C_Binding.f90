!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2025 National Renewable Energy Lab
!
! This file is part of SeaState.
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
!**********************************************************************************************************************************
MODULE SeaState_C_Binding

   USE ISO_C_BINDING
   USE SeaSt_WaveField
   USE SeaState
   USE SeaState_Types
   USE SeaState_Output
   USE NWTC_Library
   USE NWTC_C_Binding, ONLY: ErrMsgLen_C, IntfStrLen, SetErrStat_F2C, FileNameFromCString
   USE VersionInfo

   implicit none
   save

   PUBLIC :: SeaSt_C_PreInit
   PUBLIC :: SeaSt_C_Init
   PUBLIC :: SeaSt_C_CalcOutput
   PUBLIC :: SeaSt_C_End
   PUBLIC :: SeaSt_C_GetWaveFieldPointer
   PUBLIC :: SeaSt_C_SetWaveFieldPointer
   PUBLIC :: SeaSt_C_GetFluidVelAcc
   PUBLIC :: SeaSt_C_GetSurfElev
   PUBLIC :: SeaSt_C_GetSurfNorm

   !------------------------------------------------------------------------------------
   !  Debugging: DebugLevel -- passed at PreInit
   !     0  - none
   !     1  - some summary info
   !     2  - above + all position/orientation info
   !     3  - above + input files (if direct passed)
   !     4  - above + meshes
   integer(IntKi)                         :: DebugLevel = 4
   logical                                :: PreInitDone = .false.

   !------------------------------------------------------------------------------------
   !  Visualization
   type VTKvis
      character(1024)                     :: outdir
      character(1024)                     :: OutRootName       ! includes directory
      integer(IntKi)                      :: write             ! 0 off, 1 init, 2 animate
      real(DbKi)                          :: dt
      integer(IntKi)                      :: NWaveElevPts(2)   ! number of points in x/y directions
      real(SiKi), allocatable             :: WaveElevVisX(:),WaveElevVisY(:)    ! x, y locations of points
      real(SiKi), allocatable             :: WaveElevVisGrid(:,:,:)             ! the actual surface data for full time series
      integer(IntKi)                      :: tWidth = 5        ! Should calculate this, but not going to
      integer(IntKi)                      :: LastWaveIndx
      integer(IntKi)                      :: lastCount = -1
   end type VTKvis
   type(VTKvis)                           :: vtk

   !------------------------------
   !  Primary derived types
   type(SeaSt_InputType)                  :: u           !< inputs to SS
   type(SeaSt_InitInputType)              :: InitInp     !< initialization input
   type(SeaSt_InitOutputType)             :: InitOutData !< Initial output data
   type(SeaSt_ParameterType), target      :: p           !< Parameters
   type(SeaSt_OutputType)                 :: y           !< Initial output (outputs are not calculated; only the output mesh is initialized)
   type(SeaSt_MiscVarType)                :: m           !< Misc variables for optimization (not copied in glue code)
   type(SeaSt_ContinuousStateType)        :: x           !< Initial continuous states
   type(SeaSt_DiscreteStateType)          :: xd          !< Initial discrete states
   type(SeaSt_ConstraintStateType)        :: z           !< Initial guess of the constraint states
   type(SeaSt_OtherStateType)             :: OtherState  !< Initial other states            

contains


!> Set environment variables
subroutine SeaSt_C_PreInit(Gravity_C, WtrDens_C, WtrDpth_C, MSL2SWL_C, DebugLevel_In, OutVTKDir_C, WrVTK_in, WrVTK_inDT, ErrStat_C, ErrMsg_C) BIND (C, NAME='SeaSt_C_PreInit')
#ifndef IMPLICIT_DLLEXPORT
!DEC$ ATTRIBUTES DLLEXPORT :: SeaSt_C_PreInit
!GCC$ ATTRIBUTES DLLEXPORT :: SeaSt_C_PreInit
#endif
   real(c_float),              intent(in   ) :: Gravity_C
   real(c_float),              intent(in   ) :: WtrDens_C
   real(c_float),              intent(in   ) :: WtrDpth_C
   real(c_float),              intent(in   ) :: MSL2SWL_C
   integer(c_int),             intent(in   ) :: DebugLevel_In
   character(kind=c_char),     intent(in   ) :: OutVTKDir_C(IntfStrLen)       !< Directory to put all vtk output
   integer(c_int),             intent(in   ) :: WrVTK_in                      !< Write VTK outputs [0: none, 1: init only, 2: animation]
   real(c_double),             intent(in   ) :: WrVTK_inDT                    !< Timestep between VTK writes
   integer(c_int),             intent(  out) :: ErrStat_C
   character(kind=c_char),     intent(  out) :: ErrMsg_C(ErrMsgLen_C)

   character(kind=C_CHAR, len=IntfStrLen), pointer :: InputString             !< Input string as a single string with NULL chracter separating lines
   integer                          :: ErrStat, ErrStat2
   character(ErrMsgLen)             :: ErrMsg,  ErrMsg2
   integer                          :: i,j,k
   character(*), parameter          :: RoutineName = 'SeaSt_C_PreInit'

   ! Initialize error handling
   ErrStat =  ErrID_None
   ErrMsg  =  ""

   call NWTC_Init( ProgNameIn=  SeaSt_ProgDesc%Name )
   call DispCopyrightLicense(   SeaSt_ProgDesc%Name )
   call DispCompileRuntimeInfo( SeaSt_ProgDesc%Name )

   ! interface debugging
   DebugLevel = int(DebugLevel_in,IntKi)

   ! check valid debug level, show passed data if >0
   if (DebugLevel < 0_IntKi) then
      ErrStat2 = ErrID_Fatal
      ErrMsg2  = "Interface debug level must be 0 or greater"//NewLine// &
      "  0  - none"//NewLine// &
      "  1  - some summary info and variables passed through interface (init only)"//NewLine// &
      "  2  - above + all position info on all calls"
      call ShowPassedData()
      if (Failed()) return;
   elseif (DebugLevel > 0_IntKi) THEN
      call WrScr("   Interface debugging level "//trim(Num2Lstr(DebugLevel))//" requested.")
      call ShowPassedData()
   endif

   ! store environment values
   InitInp%Gravity      = Gravity_C
   InitInp%defWtrDens   = WtrDens_C
   InitInp%defWtrDpth   = WtrDpth_C
   InitInp%defMSL2SWL   = MSL2SWL_C

   !----------------------
   ! store VTK output info
   vtk%write  = int(WrVTK_in, IntKi)
   vtk%dt     = real(WrVTK_inDT, DbKi)

   if (vtk%write < 0_IntKi .or. vtk%write > 2_IntKi) then
      ErrStat2 = ErrID_Warn
      ErrMSg2  = "WrVTK_in must be 0 (off), 1 (init), 2 (animation), but "//trim(Num2LStr(vtk%write))//" was passed.  Turning off VTK surface export."
      vtk%write = 0_IntKi
      if (Failed()) return
   endif

   if (vtk%write > 0_IntKi) then
      ! Store the out root dir
      vtk%outdir = TRANSFER( OutVTKDir_C, vtk%outdir )
      i = INDEX(vtk%outdir,C_NULL_CHAR) - 1               ! if this has a c null character at the end...
      if ( i > 0 ) vtk%outdir = vtk%outdir(1:I)            ! remove it
      ! Tell SeaState to generate the visualization using default grid
      InitInp%SurfaceVis     = .true.
      InitInp%SurfaceVisNx   = 0    ! use the WaveField grid resolution
      InitInp%SurfaceVisNy   = 0    ! use the WaveField grid resolution
   endif


   ! If we got this far, we are initialized
   PreInitDone = .true.

   call Cleanup()
   return
contains
   logical function Failed()
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      Failed = ErrStat >= AbortErrLev
      if (Failed) call Cleanup()
   end function Failed
   subroutine Cleanup()    ! NOTE: we are ignoring any error reporting from here
      call SetErrStat_F2C(ErrStat,ErrMsg,ErrStat_C,ErrMsg_C)
   end subroutine Cleanup
   subroutine ShowPassedData()
      ! character(1) :: TmpFlag
      ! integer      :: i,j
      call WrScr("-----------------------------------------------------------")
      call WrScr("Interface debugging:  SeaSt_C_PreInit")
      call WrScr("   --------------------------------------------------------")
      call WrScr("   Gravity_C              -> "//trim(Num2LStr(Gravity_C)))
      call WrScr("   WtrDens_C              -> "//trim(Num2LStr(WtrDens_C)))
      call WrScr("   WtrDpth_C              -> "//trim(Num2LStr(WtrDpth_C)))
      call WrScr("   MSL2SWL_C              -> "//trim(Num2LStr(MSL2SWL_C)))
      call WrScr("   DebugLevel_In          -> "//trim(Num2LStr(DebugLevel_In)))
      call WrScr("   OutVTKDir_C (ptr addr) -> "//trim(Num2LStr(LOC(OutVTKDir_C))))
      call WrScr("   WrVTK_in               -> "//trim(Num2LStr(WrVTK_in)))
      call WrScr("   WrVTK_inDT             -> "//trim(Num2LStr(WrVTK_inDT)))
      call WrScr("-----------------------------------------------------------")
   end subroutine ShowPassedData
end subroutine SeaSt_C_PreInit


!> Initialize the library (PreInit must be called first)
subroutine SeaSt_C_Init(InputFile_C, OutRootName_C, NSteps_C, TimeInterval_C, WaveElevSeriesFlag_C, WrWvKinMod_C, NumChannels_C, OutputChannelNames_C, OutputChannelUnits_C, ErrStat_C, ErrMsg_C) BIND (C, NAME='SeaSt_C_Init')
#ifndef IMPLICIT_DLLEXPORT
!DEC$ ATTRIBUTES DLLEXPORT :: SeaSt_C_Init
!GCC$ ATTRIBUTES DLLEXPORT :: SeaSt_C_Init
#endif
   type(c_ptr),                intent(in   ) :: InputFile_C
   type(c_ptr),                intent(in   ) :: OutRootName_C
   integer(c_int),             intent(in   ) :: NSteps_C
   real(c_float),              intent(in   ) :: TimeInterval_C
   integer(c_int),             intent(in   ) :: WaveElevSeriesFlag_C
   integer(c_int),             intent(in   ) :: WrWvKinMod_C
   integer(c_int),             intent(  out) :: NumChannels_C
   character(kind=c_char),     intent(  out) :: OutputChannelNames_C(ChanLen*MaxOutPts+1)
   character(kind=c_char),     intent(  out) :: OutputChannelUnits_C(ChanLen*MaxOutPts+1)
   integer(c_int),             intent(  out) :: ErrStat_C
   character(kind=c_char),     intent(  out) :: ErrMsg_C(ErrMsgLen_C)

   ! Local variables
   character(kind=C_CHAR, len=IntfStrLen), pointer :: InputFileString          !< Input file as a single string with NULL chracter separating lines
   character(kind=C_CHAR, len=IntfStrLen), pointer :: OutputFileString          !< Input file as a single string with NULL chracter separating lines
   character(IntfStrLen)            :: InputFileName
   character(IntfStrLen)            :: OutRootName
   character(1024)                  :: vtkroot
   real(DbKi)                       :: Interval        !< DT for calling
   integer                          :: ErrStat, ErrStat2
   character(ErrMsgLen)             :: ErrMsg,  ErrMsg2
   integer                          :: i,j,k
   character(*), parameter          :: RoutineName = 'SeaSt_C_Init'  !< for error handling

   if (.not. PreInitDone) then
      ErrStat = ErrID_Fatal
      ErrMSg  = "SeaSt_C_PreInit must be called before SeaSt_C_Init"
      call Cleanup()
   endif

   ! Initialize error handling
   ErrStat =  ErrID_None
   ErrMsg  =  ""

   call NWTC_Init( ProgNameIn=  SeaSt_ProgDesc%Name )
   call DispCopyrightLicense(   SeaSt_ProgDesc%Name )
   call DispCompileRuntimeInfo( SeaSt_ProgDesc%Name )

   if (DebugLevel > 0_IntKi) call ShowPassedData()

   ! Input files
   call C_F_POINTER(InputFile_C, InputFileString)  ! Get a pointer to the input file string
   InputFileName = FileNameFromCString(InputFileString, IntfStrLen)  ! convert the input file name from c_char to fortran character

   call C_F_POINTER(OutRootName_C, OutputFileString)  ! Get a pointer to the input file string
   OutRootName = FileNameFromCString(OutputFileString, IntfStrLen)  ! convert the input file name from c_char to fortran character

   ! Set other inputs for calling SeaSt_Init
   InitInp%InputFile    = InputFileName
   InitInp%UseInputFile = .TRUE. 
   InitInp%OutRootName  = OutRootName
   InitInp%TMax         = (NSteps_C - 1) * TimeInterval_C   ! Using this to match the SeaState driver; could otherwise get TMax directly
   InitInp%WaveFieldMod = WaveElevSeriesFlag_C
   ! REAL(ReKi)  :: PtfmLocationX = 0.0_ReKi      !< Supplied by Driver:  X coordinate of platform location in the wave field [m]
   ! REAL(ReKi)  :: PtfmLocationY = 0.0_ReKi      !< Supplied by Driver:  Y coordinate of platform location in the wave field [m]
   InitInp%WrWvKinMod = WrWvKinMod_C
   ! LOGICAL  :: HasIce = .false.      !< Supplied by Driver:  Whether this simulation has ice loading (flag) [-]
   ! LOGICAL  :: Linearize = .FALSE.      !< Flag that tells this module if the glue code wants to linearize. [-]
   ! LOGICAL  :: SurfaceVis = .FALSE.      !< Turn on grid surface visualization outputs [-]
   ! INTEGER(IntKi)  :: SurfaceVisNx = 0      !< Number of points in X direction to output for visualization grid.  Use 0 or negative to set to SeaState resolution. [-]
   ! INTEGER(IntKi)  :: SurfaceVisNy = 0      !< Number of points in Y direction to output for visualization grid.  Use 0 or negative to set to SeaState resolution. [-]

   call SeaSt_Init( InitInp, u, p, x, xd, z, OtherState, y, m, Interval, InitOutData, ErrStat2, ErrMsg2 )
      if (Failed()) return

   ! Number of channels
   NumChannels_C = size(InitOutData%WriteOutputHdr)

   ! transfer the output channel names and units to c_char arrays for returning
   k=1
   do i=1,NumChannels_C
      do j=1,ChanLen    ! max length of channel name.  Same for units
         OutputChannelNames_C(k)=InitOutData%WriteOutputHdr(i)(j:j)
         OutputChannelUnits_C(k)=InitOutData%WriteOutputUnt(i)(j:j)
         k=k+1
      enddo
   enddo

   ! null terminate the string
   OutputChannelNames_C(k) = C_NULL_CHAR
   OutputChannelUnits_C(k) = C_NULL_CHAR

   if (vtk%write > 0_IntKi) then
      ! check dt (can't check against Interval since that is never set, so just make sure it is positive)
      if (vtk%dt <= 0.0) vtk%dt = 0.25
      if (allocated(InitOutData%WaveElevVisGrid)) then
         vtk%NWaveElevPts(1) = size(InitOutData%WaveElevVisX)
         vtk%NWaveElevPts(2) = size(InitOutData%WaveElevVisY)
         call move_alloc(InitOutData%WaveElevVisX, vtk%WaveElevVisX)
         call move_alloc(InitOutData%WaveElevVisY, vtk%WaveElevVisY)
         call move_alloc(InitOutData%WaveElevVisGrid,vtk%WaveElevVisGrid )
      else
         vtk%NWaveElevPts = 0
         vtk%write = 0     ! FIXME throw warning if we do this
      endif
      ! get the name of the output directory for vtk files (in a subdirectory called "vtk" of the output directory), and
      ! create the VTK directory if it does not exist
      call GetPath ( OutRootName, vtk%OutRootName, vtkroot ) ! the returned vtk%OutRootName includes a file separator character at the end
      if (PathIsRelative(trim(vtk%OutRootName))) then
         vtk%OutRootName = trim(vtk%OutRootName) // trim(vtk%outdir)
      else
         vtk%OutRootName = trim(vtk%outdir)
      endif
      call MKDIR( trim(vtk%OutRootName) )
      vtk%OutRootName = trim( vtk%OutRootName ) // PathSep // trim(vtkroot)
      call WrVTK_WaveElevVisGrid  (0.0_DbKi, vtk, ErrStat2, ErrMsg2)
      if (Failed()) return
   endif



   call Cleanup()

contains
   logical function Failed()
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      Failed = ErrStat >= AbortErrLev
      if (Failed) call Cleanup()
   end function Failed
   subroutine Cleanup()    ! NOTE: we are ignoring any error reporting from here
      call SetErrStat_F2C(ErrStat,ErrMsg,ErrStat_C,ErrMsg_C)
   end subroutine Cleanup
   subroutine ShowPassedData()
      ! character(1) :: TmpFlag
      ! integer      :: i,j
      call WrScr("-----------------------------------------------------------")
      call WrScr("Interface debugging:  SeaSt_C_Init")
      call WrScr("   --------------------------------------------------------")
      call WrScr("   FIXME: THIS SECTION IS MISSING!!!!!!!")
      call WrScr("-----------------------------------------------------------")
   end subroutine ShowPassedData
end subroutine SeaSt_C_Init

subroutine SeaSt_C_CalcOutput(Time_C, OutputChannelValues_C, ErrStat_C, ErrMsg_C) BIND (C, NAME='SeaSt_C_CalcOutput')
#ifndef IMPLICIT_DLLEXPORT
!DEC$ ATTRIBUTES DLLEXPORT :: SeaSt_C_CalcOutput
!GCC$ ATTRIBUTES DLLEXPORT :: SeaSt_C_CalcOutput
#endif

   real(c_double),             intent(in   ) :: Time_C
   real(c_float),              intent(  out) :: OutputChannelValues_C(p%NumOuts)
   integer(c_int),             intent(  out) :: ErrStat_C
   character(kind=c_char),     intent(  out) :: ErrMsg_C(ErrMsgLen_C)

   ! Local variables
   type(SeaSt_InputType)           :: u           !< An initial guess for the input; input mesh must be defined
   type(SeaSt_ContinuousStateType) :: x           !< Initial continuous states
   type(SeaSt_DiscreteStateType)   :: xd          !< Initial discrete states
   type(SeaSt_ConstraintStateType) :: z           !< Initial guess of the constraint states
   type(SeaSt_OtherStateType)      :: OtherState  !< Initial other states            

   real(DbKi)                 :: Time
   integer                    :: ErrStat, ErrStat2
   character(ErrMsgLen)       :: ErrMsg,  ErrMsg2
   character(*), parameter    :: RoutineName = 'SeaSt_C_CalcOutput'  !< for error handling

   ! Initialize error handling
   ErrStat =  ErrID_None
   ErrMsg  =  ""

   ! Debugging
   if (DebugLevel > 1) call ShowPassedData()

   ! Convert the inputs from C to Fortran
   Time = REAL(Time_C,DbKi)

   call SeaSt_CalcOutput( Time, u, p, x, xd, z, OtherState, y, m, ErrStat2, ErrMsg2 )
   if (Failed()) return

   ! Get the output channel info out of y
   OutputChannelValues_C = REAL(y%WriteOutput, C_FLOAT)

   if (vtk%write > 1_IntKi) then
      call WrVTK_WaveElevVisGrid  (Time, vtk, ErrStat2, ErrMsg2)
      if (Failed()) return
   endif

   call Cleanup()

   ! Debugging
   if (DebugLevel > 1) call ShowReturnData()

contains
   logical function Failed()
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      Failed = ErrStat >= AbortErrLev
      if (Failed) call Cleanup()
   end function Failed
   subroutine Cleanup()    ! NOTE: we are ignoring any error reporting from here
      CALL SetErrStat_F2C(ErrStat,ErrMsg,ErrStat_C,ErrMsg_C)
   end subroutine Cleanup
   subroutine ShowPassedData()
      call WrScr("-----------------------------------------------------------")
      call WrScr("Interface debugging:  SeaSt_C_CalcOutput")
      call WrScr("   --------------------------------------------------------")
      call WrScr("   Time_C                 -> "//trim(Num2LStr(Time_C)))
   end subroutine ShowPassedData
   subroutine ShowReturnData()
      call WrScr("   OutputChannelValues_C  <-")
      call WrMatrix(OutputChannelValues_C,CU,'g15.6')
      call WrScr("-----------------------------------------------------------")
   end subroutine ShowReturnData
end subroutine

subroutine SeaSt_C_End(ErrStat_C,ErrMsg_C) BIND (C, NAME='SeaSt_C_End')
#ifndef IMPLICIT_DLLEXPORT
!DEC$ ATTRIBUTES DLLEXPORT :: SeaSt_C_End
!GCC$ ATTRIBUTES DLLEXPORT :: SeaSt_C_End
#endif
   integer(C_INT),             intent(  out) :: ErrStat_C
   character(kind=C_CHAR),     intent(  out) :: ErrMsg_C(ErrMsgLen_C)

   integer                    :: ErrStat                          !< aggregated error status
   character(ErrMsgLen)       :: ErrMsg                           !< aggregated error message
   integer                    :: ErrStat2                         !< temporary error status  from a call
   character(ErrMsgLen)       :: ErrMsg2                          !< temporary error message from a call
   character(*), parameter    :: RoutineName = 'SeaSt_C_End'  !< for error handling
   ErrStat  =  ErrID_None
   ErrMsg   =  ""
   call SeaSt_End(u, p, x, xd, z, OtherState, y, m, ErrStat2, ErrMsg2)
   call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   call SetErrStat_F2C( ErrStat, ErrMsg, ErrStat_C, ErrMsg_C )
end subroutine


!> return the pointer to the WaveField data
subroutine SeaSt_C_GetWaveFieldPointer(WaveFieldPointer_C,ErrStat_C,ErrMsg_C) BIND (C, NAME='SeaSt_C_GetWaveFieldPointer')
#ifndef IMPLICIT_DLLEXPORT
!DEC$ ATTRIBUTES DLLEXPORT :: SeaSt_C_GetWaveFieldPointer
!GCC$ ATTRIBUTES DLLEXPORT :: SeaSt_C_GetWaveFieldPointer
#endif
   type(c_ptr),               intent(  out)  :: WaveFieldPointer_C
   integer(c_int),            intent(  out)  :: ErrStat_C
   character(kind=c_char),    intent(  out)  :: ErrMsg_C(ErrMsgLen_C)
   integer                                   :: ErrStat
   character(ErrMsgLen)                      :: ErrMsg
   character(*),              parameter      :: RoutineName = 'SeaSt_C_GetWaveFieldPointer'
   ErrStat = ErrID_None
   ErrMSg = ""
   if (associated(p%WaveField)) then
      WaveFieldPointer_C = C_LOC(p%WaveField)
   else
      WaveFieldPointer_C = C_NULL_PTR
      call SetErrStat(ErrID_Fatal,"Pointer to WaveField data not valid: data not initialized",ErrStat,ErrMsg,RoutineName)
   endif
   call SetErrStat_F2C( ErrStat, ErrMsg, ErrStat_C, ErrMsg_C )
   if (DebugLevel > 1) call ShowPassedData()
   return
contains
   subroutine ShowPassedData()
      call WrScr("-----------------------------------------------------------")
      call WrScr("Interface debugging:  SeaSt_C_GetWaveFieldPointer")
      call WrScr("   --------------------------------------------------------")
      call WrScr("   WaveFieldPointer_C     -> "//trim(Num2LStr(loc(p%WaveField))))
      call WrScr("-----------------------------------------------------------")
   end subroutine ShowPassedData
end subroutine


!> set the pointer to the WaveField data
subroutine SeaSt_C_SetWaveFieldPointer(WaveFieldPointer_C,ErrStat_C,ErrMsg_C) BIND (C, NAME='SeaSt_C_SetWaveFieldPointer')
#ifndef IMPLICIT_DLLEXPORT
!DEC$ ATTRIBUTES DLLEXPORT :: SeaSt_C_SetWaveFieldPointer
!GCC$ ATTRIBUTES DLLEXPORT :: SeaSt_C_SetWaveFieldPointer
#endif
   type(c_ptr),               intent(in   )  :: WaveFieldPointer_C
   integer(c_int),            intent(  out)  :: ErrStat_C
   character(kind=c_char),    intent(  out)  :: ErrMsg_C(ErrMsgLen_C)
   integer                                   :: ErrStat
   character(ErrMsgLen)                      :: ErrMsg
   character(*),              parameter      :: RoutineName = 'SeaSt_C_SetWaveFieldPointer'
   ErrStat = ErrID_None
   ErrMSg = ""
   call C_F_POINTER(WaveFieldPointer_C, p%WaveField)
   if (associated(p%WaveField)) then
      ! basic sanity check
      if (.not. allocated(p%WaveField%WaveTime)) then
         call SetErrStat(ErrID_Fatal,"Invalid pointer passed in, or WaveField not initialized",ErrStat,ErrMsg,RoutineName)
      endif
   else
      call SetErrStat(ErrID_Fatal,"Invalid pointer passed in, or WaveField not initialized",ErrStat,ErrMsg,RoutineName)
   endif
   call SetErrStat_F2C( ErrStat, ErrMsg, ErrStat_C, ErrMsg_C )
   if (DebugLevel > 1) call ShowPassedData()
   return
contains
   subroutine ShowPassedData()
      call WrScr("-----------------------------------------------------------")
      call WrScr("Interface debugging:  SeaSt_C_SetWaveFieldPointer")
      call WrScr("   --------------------------------------------------------")
      call WrScr("   WaveFieldPointer_C     <- "//trim(Num2LStr(loc(p%WaveField))))
      call WrScr("-----------------------------------------------------------")
   end subroutine ShowPassedData
end subroutine


!> Get the fluid velocity, acceleration, and node-in-water status at time+position coordinate
subroutine SeaSt_C_GetFluidVelAcc(Time_C, Pos_C, Vel_C, Acc_C, NodeInWater_C, ErrStat_C,ErrMsg_C) BIND (C, NAME='SeaSt_C_GetFluidVelAcc')
#ifndef IMPLICIT_DLLEXPORT
!DEC$ ATTRIBUTES DLLEXPORT :: SeaSt_C_GetFluidVelAcc
!GCC$ ATTRIBUTES DLLEXPORT :: SeaSt_C_GetFluidVelAcc
#endif
   real(c_double),            intent(in   ) :: Time_C
   real(c_float),             intent(in   ) :: Pos_c(3)
   real(c_float),             intent(  out) :: Vel_c(3)
   real(c_float),             intent(  out) :: Acc_c(3)
   integer(c_int),            intent(  out) :: NodeInWater_C
   integer(c_int),            intent(  out) :: ErrStat_C
   character(kind=c_char),    intent(  out) :: ErrMsg_C(ErrMsgLen_C)
   real(DbKi)                 :: Time
   real(ReKi)                 :: Pos(3)
   real(SiKi)                 :: Vel(3)
   real(SiKi)                 :: Acc(3)
   logical                    :: forceNodeInWater
!FIXME:dev-tc uncomment next line
!   logical                    :: fetchDynCurrent
   integer(IntKi)             :: nodeInWater
   integer                    :: ErrStat, ErrStat2
   character(ErrMsgLen)       :: ErrMsg,  ErrMsg2
   character(*), parameter    :: RoutineName = 'SeaSt_C_GetFluidVelAcc'

   ! Initialize error handling
   ErrStat  =  ErrID_None
   ErrMsg   =  ""
   
   forceNodeInWater = .false.

   if (DebugLevel > 1) call ShowPassedData()

   ! convert position and time to fortran types
   Time = real(Time_C, DbKi)
   Pos = real(Pos_C, ReKi)

   ! get wave field velocity and acceleration (current is included in this)
   ! Notes:
   !     - if node is out of water, velocity and acceleration are zero
   !     - if position is outside the wave field boundary, it will simply return boundary edge value
   !     - time must be positive or a fatal error occurs
   call WaveField_GetNodeWaveVelAcc( p%WaveField, m%WaveField_m, Time, pos, forceNodeInWater, nodeInWater, Vel, Acc, ErrStat, ErrMsg )
!FIXME:dev-tc use next line instead of above
!   call WaveField_GetNodeWaveVelAcc( p%WaveField, m%WaveField_m, Time, pos, forceNodeInWater, fetchDynCurrent, nodeInWater, Vel, Acc, ErrStat, ErrMsg )

   ! Store resulting velocity and acceleration as C type
   Vel_c = real(Vel,c_float)
   Acc_c = real(Acc,c_float)

   ! Density value and node status to return
   if (nodeInWater == 1_IntKi) then
      NodeInWater_C = 1_c_int
   else
      NodeInWater_C = 0_c_int
   endif

   call SetErrStat_F2C( ErrStat, ErrMsg, ErrStat_C, ErrMsg_C )    ! convert error from fortran to C for return
   if (DebugLevel > 1) call ShowReturnData()
   return
contains
   subroutine ShowPassedData()
      call WrScr("-----------------------------------------------------------")
      call WrScr("Interface debugging:  SeaSt_C_GetFluidVelAccDens")
      call WrScr("   --------------------------------------------------------")
      call WrScr("   Time_C                 -> "//trim(Num2LStr(Time_C)))
      call WrScr("   Pos_C                  -> ("//trim(Num2LStr(Pos_C(1)))//","//trim(Num2LStr(Pos_C(2)))//","//trim(Num2LStr(Pos_C(3)))//")")
   end subroutine ShowPassedData
   subroutine ShowReturnData()
      call WrScr("   Vel_C                  <- ("//trim(Num2LStr(Vel_C(1)))//","//trim(Num2LStr(Vel_C(2)))//","//trim(Num2LStr(Vel_C(3)))//")")
      call WrScr("   Acc_C                  <- ("//trim(Num2LStr(Acc_C(1)))//","//trim(Num2LStr(Acc_C(2)))//","//trim(Num2LStr(Acc_C(3)))//")")
      call WrScr("   NodeInWater_C          <- "//trim(Num2LStr(NodeInWater_C)))
      call WrScr("-----------------------------------------------------------")
   end subroutine ShowReturnData
end subroutine SeaSt_C_GetFluidVelAcc



!> return the surface elevation and normal vector at a point.
subroutine SeaSt_C_GetSurfElev(Time_C, Pos_C, Elev_C, ErrStat_C,ErrMsg_C) BIND (C, NAME='SeaSt_C_GetSurfElev')
#ifndef IMPLICIT_DLLEXPORT
!DEC$ ATTRIBUTES DLLEXPORT :: SeaSt_C_GetSurfElev
!GCC$ ATTRIBUTES DLLEXPORT :: SeaSt_C_GetSurfElev
#endif
   real(c_double),            intent(in   ) :: Time_C
   real(c_float),             intent(in   ) :: Pos_c(2)
   real(c_float),             intent(  out) :: Elev_C
   integer(c_int),            intent(  out) :: ErrStat_C
   character(kind=c_char),    intent(  out) :: ErrMsg_C(ErrMsgLen_C)

   real(DbKi)                 :: Time
   real(ReKi)                 :: Pos(2)
   real(SiKi)                 :: Elev
   integer                    :: ErrStat
   character(ErrMsgLen)       :: ErrMsg
   character(*), parameter    :: RoutineName = 'SeaSt_C_GetSurfElev'

   if (DebugLevel > 1) call ShowPassedData()

   ! convert position and time to fortran types
   Time = real(Time_C, DbKi)
   Pos = real(Pos_C(1:2), ReKi)

   ! get wave elevation (total combined first and second order)
   ! Notes:
   !     - if position is outside the wave field boundary, it will simply return boundary edge value
   !     - time must be positive or a fatal error occurs
   Elev = WaveField_GetNodeTotalWaveElev( p%WaveField, m%WaveField_m, Time, pos, ErrStat, ErrMsg )

   ! Store resulting elevation as C type
   Elev_C = real(Elev,c_float)

   if (DebugLevel > 1) call ShowReturnData()
   call Cleanup()
   return
contains
   subroutine Cleanup()    ! NOTE: we are ignoring any error reporting from here
      CALL SetErrStat_F2C(ErrStat,ErrMsg,ErrStat_C,ErrMsg_C)
   end subroutine Cleanup
   subroutine ShowPassedData()
      call WrScr("-----------------------------------------------------------")
      call WrScr("Interface debugging:  SeaSt_C_GetSurfElev")
      call WrScr("   --------------------------------------------------------")
      call WrScr("   Time_C                 -> "//trim(Num2LStr(Time_C)))
      call WrScr("   Pos_C                  -> ("//trim(Num2LStr(Pos_C(1)))//","//trim(Num2LStr(Pos_C(2)))//")")
   end subroutine ShowPassedData
   subroutine ShowReturnData()
      call WrScr("   Elev_C                 <- "//trim(Num2LStr(Elev_C)))
      call WrScr("-----------------------------------------------------------")
   end subroutine ShowReturnData
end subroutine SeaSt_C_GetSurfElev



!> return the surface normal vector at a point.
subroutine SeaSt_C_GetSurfNorm(Time_C, Pos_C, NormVec_C, ErrStat_C,ErrMsg_C) BIND (C, NAME='SeaSt_C_GetSurfNorm')
#ifndef IMPLICIT_DLLEXPORT
!DEC$ ATTRIBUTES DLLEXPORT :: SeaSt_C_GetSurfNorm
!GCC$ ATTRIBUTES DLLEXPORT :: SeaSt_C_GetSurfNorm
#endif
   real(c_double),            intent(in   ) :: Time_C
   real(c_float),             intent(in   ) :: Pos_c(2)
   real(c_float),             intent(  out) :: NormVec_C(3)
   integer(c_int),            intent(  out) :: ErrStat_C
   character(kind=c_char),    intent(  out) :: ErrMsg_C(ErrMsgLen_C)

   real(DbKi)                 :: Time
   real(ReKi)                 :: Pos(2), r
   real(SiKi)                 :: NormVec(3)
   integer                    :: ErrStat
   character(ErrMsgLen)       :: ErrMsg
   character(*), parameter    :: RoutineName = 'SeaSt_C_GetSurfNorm'

   if (DebugLevel > 1) call ShowPassedData()

   ! convert position and time to fortran types
   Time = real(Time_C, DbKi)
   Pos = real(Pos_C(1:2), ReKi)

   ! get the normal vector at the point (set to vertical if outside region)
   r = 0.0_ReKi   ! distance for central differencing - use default in code
   call WaveField_GetNodeWaveNormal( p%WaveField, m%WaveField_m, Time, pos, r, NormVec, ErrStat, ErrMsg )

   ! Store resulting normal vector as C type
   NormVec_C = real(NormVec,c_float)

   if (DebugLevel > 1) call ShowReturnData()
   call Cleanup()
   return
contains
   subroutine Cleanup()    ! NOTE: we are ignoring any error reporting from here
      CALL SetErrStat_F2C(ErrStat,ErrMsg,ErrStat_C,ErrMsg_C)
   end subroutine Cleanup
   subroutine ShowPassedData()
      call WrScr("-----------------------------------------------------------")
      call WrScr("Interface debugging:  SeaSt_C_GetSurfNorm")
      call WrScr("   --------------------------------------------------------")
      call WrScr("   Time_C                 -> "//trim(Num2LStr(Time_C)))
      call WrScr("   Pos_C                  -> ("//trim(Num2LStr(Pos_C(1)))//","//trim(Num2LStr(Pos_C(2)))//")")
   end subroutine ShowPassedData
   subroutine ShowReturnData()
      call WrScr("   NormVec_C              <- ("//trim(Num2LStr(NormVec_C(1)))//","//trim(Num2LStr(NormVec_C(2)))//","//trim(Num2LStr(NormVec_C(3)))//")")
      call WrScr("-----------------------------------------------------------")
   end subroutine ShowReturnData
end subroutine SeaSt_C_GetSurfNorm




!FIXME: the following visualization writer should be merged into the library vtk.f90
!        this is a modified duplicate of the routine from FAST_Subs by the same name.
!----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine writes the wave elevation data for a given time step
subroutine WrVTK_WaveElevVisGrid(t_global, vtk, ErrStat, ErrMsg)
   real(DbKi),                intent(in   )  :: t_global       !< Current global time
   type(VTKvis),              intent(inout)  :: vtk
   integer(IntKi),            intent(  out)  :: ErrStat
   character(ErrMsgLen),      intent(  out)  :: ErrMsg
   integer(IntKi)                            :: Un             ! fortran unit number
   integer(IntKi)                            :: n, iy, ix      ! loop counters
   real(SiKi)                                :: t
   integer(IntKi)                            :: count          ! which file is this?
   character(1024)                           :: FileName
   integer(IntKi)                            :: NumberOfPoints
   integer(IntKi), parameter                 :: NumberOfLines = 0
   integer(IntKi)                            :: NumberOfPolys
   character(1024)                           :: Tstr
   integer(IntKi)                            :: ErrStat2
   character(ErrMsgLen)                      :: ErrMsg2
   character(*),parameter                    :: RoutineName = 'WrVTK_WaveElevVisGrid'

   ! Initialize error handling
   ErrStat = ErrID_None
   ErrMsg  = ""

   NumberOfPoints = vtk%NWaveElevPts(1) * vtk%NWaveElevPts(2)
   ! I'm going to make triangles for now. we should probably just make this a structured file at some point
   NumberOfPolys  = ( vtk%NWaveElevPts(1) - 1 ) * &
                    ( vtk%NWaveElevPts(2) - 1 ) * 2

   count = nint(t_global / vtk%dt)
   if (count == vtk%lastCount) return     ! already wrote this one
   vtk%lastCount = count                  ! store the current number to make sure we don't write it twice

   !.................................................................
   ! write the data that potentially changes each time step:
   !.................................................................
   ! construct the string for the zero-padded VTK write-out step
   write(Tstr, '(i' // trim(Num2LStr(vtk%tWidth)) //'.'// trim(Num2LStr(vtk%tWidth)) // ')') count

   ! PolyData (.vtp) - Serial vtkPolyData (unstructured) file
   FileName = TRIM(vtk%OutRootName)//'.WaveSurface.'//TRIM(Tstr)//'.vtp'

   call WrVTK_header( FileName, NumberOfPoints, NumberOfLines, NumberOfPolys, Un, ErrStat, ErrMsg )
      if (ErrStat >= AbortErrLev) return

! points (nodes, augmented with NumSegments):
   write(Un,'(A)')         '      <Points>'
   write(Un,'(A)')         '        <DataArray type="Float32" NumberOfComponents="3" format="ascii">'

   ! I'm not going to interpolate in time; I'm just going to get the index of the closest wave time value
   t = REAL(t_global,SiKi)
   call GetWaveElevIndx( t, p%WaveField%WaveTime, vtk%LastWaveIndx )

   do ix=1,vtk%NWaveElevPts(1)
      do iy=1,vtk%NWaveElevPts(2)
         write(Un,VTK_AryFmt) vtk%WaveElevVisX(ix), vtk%WaveElevVisY(iy), vtk%WaveElevVisGrid(vtk%LastWaveIndx,ix,iy)
      end do
   end do

   write(Un,'(A)')         '        </DataArray>'
   write(Un,'(A)')         '      </Points>'
   write(Un,'(A)')         '      <Polys>'
   write(Un,'(A)')         '        <DataArray type="Int32" Name="connectivity" format="ascii">'

   do ix=1,vtk%NWaveElevPts(1)-1
      do iy=1,vtk%NWaveElevPts(2)-1
         n = vtk%NWaveElevPts(2)*(ix-1)+iy - 1 ! points start at 0
         write(Un,'(3(i7))') n,   n+1,                                    n+vtk%NWaveElevPts(2)
         write(Un,'(3(i7))') n+1, n+1+vtk%NWaveElevPts(2), n+vtk%NWaveElevPts(2)
      end do
   end do
   write(Un,'(A)')         '        </DataArray>'
   write(Un,'(A)')         '        <DataArray type="Int32" Name="offsets" format="ascii">'
   do n=1,NumberOfPolys
      WRITE(Un,'(i7)') 3*n
   end do
   write(Un,'(A)')         '        </DataArray>'
   write(Un,'(A)')         '      </Polys>'
   call WrVTK_footer( Un )
contains
   !----------------------------------------------------------------------------------------------------------------------------------
   !> This function returns the index, Ind, of the XAry closest to XValIn, where XAry is assumed to be periodic. It starts
   !! searching at the value of Ind from a previous step.
   subroutine GetWaveElevIndx( XValIn, XAry, Ind )
      integer, intent(inout)       :: Ind                ! Initial and final index into the arrays.
      real(SiKi), intent(in)       :: XAry    (:)        !< Array of X values to be interpolated.
      real(SiKi), intent(in)       :: XValIn             !< X value to be found
      integer                      :: AryLen             ! Length of the arrays.
      real(SiKi)                   :: XVal               !< X to be found (wrapped/periodic)

      AryLen = size(XAry)

      ! Wrap XValIn into the range XAry(1) to XAry(AryLen)
      XVal = MOD(XValIn, XAry(AryLen))

      ! Let's check the limits first.
      if ( XVal <= XAry(1) )  then
         Ind = 1
         return
      else if ( XVal >= XAry(AryLen) )  then
         Ind = AryLen
         return
      else
         ! Set the Ind to the first index if we are at the beginning of XAry
         if ( XVal <= XAry(2) )  then
            Ind = 1
         end if
      end if

      ! Let's interpolate!
      Ind = MAX( MIN( Ind, AryLen-1 ), 1 )
      do
         if ( XVal < XAry(Ind) )  then
            Ind = Ind - 1
         else if ( XVal >= XAry(Ind+1) )  then
            Ind = Ind + 1
         else
            ! XAry(Ind) <= XVal < XAry(Ind+1)
            ! this would make it the "closest" node, but I'm not going to worry about that for visualization purposes
            !if ( XVal > (XAry(Ind+1) + XAry(Ind))/2.0_SiKi ) Ind = Ind + 1
            return
         end if
      end do
      return
   end subroutine GetWaveElevIndx
end subroutine WrVTK_WaveElevVisGrid
end module SeaState_C_Binding
