!**********************************************************************************************************************************
! UnsteadyAero_DriverCode: This code tests a stand-alone version of the UnsteadyAero module
!..................................................................................................................................
! LICENSING
! Copyright (C) 2015  National Renewable Energy Laboratory
!
!    This file is part of UnsteadyAero.
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
program UnsteadyAero_Driver

   use NWTC_Library
   use AirfoilInfo
   use AirfoilInfo_Types
   use UnsteadyAero_Types
   use UnsteadyAero
   use UA_Dvr_Subs
   use VersionInfo

   use LinDyn

   implicit none
    ! Variables
   
   real(DbKi)  :: t
   integer     :: i, j, n, iu

   ! --- All Data
   type(Dvr_Data)         :: dvr

   integer(IntKi)                                :: ErrStat              ! Status of error message
   character(ErrMsgLen)                          :: ErrMsg               ! Error message if ErrStat /= ErrID_None
   
   CHARACTER(1024)                               :: dvrFilename          ! Filename and path for the driver input file.  This is passed in as a command line argument when running the Driver exe.
   character(*), parameter                       :: RoutineName = 'UnsteadyAero_Driver'
   
   CHARACTER(200)                                :: git_commit
   TYPE(ProgDesc), PARAMETER   :: version   = ProgDesc( 'UnsteadyAero Driver', '', '' )  ! The version number of this program.
      ! Initialize the NWTC library
   call NWTC_Init()
   
      ! Initialize error handling variables
   ErrMsg  = ''
   ErrStat = ErrID_None
   
      ! Display the copyright notice
   CALL DispCopyrightLicense( version%Name )
      ! Obtain OpenFAST git commit hash
   git_commit = QueryGitVersion()
      ! Tell our users what they're running
   CALL WrScr( ' Running '//TRIM( version%Name )//' a part of OpenFAST - '//TRIM(git_Commit)//NewLine//' linked with '//TRIM( NWTC_Ver%Name )//NewLine )
   
   
   ! --- Parse the driver file if one
   if ( command_argument_count() > 1 ) then
      call print_help()
      call NormStop()
   endif
   call get_command_argument(1, dvrFilename)
   call ReadDriverInputFile( dvrFilename, dvr%p, errStat, errMsg ); call checkError()

   ! --- Driver Parameters
   call Dvr_SetParameters(dvr%p, errStat, errMsg); call checkError()

   ! --- Initialize Elastic Section
   if ( dvr%p%SimMod == 3 ) then
      call LD_InitInputData(3, dvr%LD_InitInData, errStat, errMsg); call checkError()
      dvr%LD_InitInData%dt        = dvr%p%dt
      dvr%LD_InitInData%IntMethod = 1  ! 1=RK4, TODO expose to user
      dvr%LD_InitInData%prefix    = '' ! for output channel names
      dvr%LD_InitInData%MM         = dvr%p%MM
      dvr%LD_InitInData%CC         = dvr%p%CC
      dvr%LD_InitInData%KK         = dvr%p%KK
      dvr%LD_InitInData%x0         = dvr%p%initPos
      dvr%LD_InitInData%xd0        = dvr%p%initVel
      dvr%LD_InitInData%activeDOFs = dvr%p%activeDOFs
      dvr%LD_InitInData%DOFsNames = (/'x  ','y  ','th '/)
      dvr%LD_InitInData%DOFsUnits = (/'m  ','m  ','rad'/)
      call LD_Init(dvr%LD_InitInData, dvr%LD_u(1), dvr%LD_p, dvr%LD_x, dvr%LD_xd, dvr%LD_z, dvr%LD_OtherState, dvr%LD_y, dvr%LD_m, dvr%LD_InitOutData, errStat, errMsg); call checkError()
      ! Allocate other inputs of LD
      do iu = 2,NumInp
         call AllocAry(dvr%LD_u(iu)%Fext, dvr%LD_p%nx, 'Fext', errStat, errMsg); call checkError()
      enddo
   end if

   ! --- Init UA input data based on driver inputs
   call driverInputsToUAInitData(dvr%p, dvr%UA_InitInData, dvr%AFI_Params, dvr%AFIndx, errStat, errMsg); call checkError()

   ! --- Initialize UnsteadyAero (need AFI)
   if ( dvr%p%SimMod == 3 ) then
      ! TODO
      dvr%UA_InitInData%UA_OUTS = 1  ! 0=None, 1=Write Outputs, 2=Separate File
   endif
   call UA_Init( dvr%UA_InitInData, dvr%UA_u(1), dvr%UA_p, dvr%UA_x, dvr%UA_xd, dvr%UA_OtherState, dvr%UA_y, dvr%UA_m, dvr%p%dt, dvr%AFI_Params, dvr%AFIndx, dvr%UA_InitOutData, errStat, errMsg ); call checkError()
   if (dvr%UA_p%NumOuts <= 0) then
      ErrStat = ErrID_Warn
      ErrMsg = "No outputs from UA are generated."
      call checkError()
   end if

   ! --- Driver Outputs
   dvr%out%Root = dvr%p%OutRootName
   if ( dvr%p%SimMod == 3 ) then
      call Dvr_InitializeDriverOutputs(dvr, dvr%out, errStat, errMsg); call checkError()
   endif

   ! --- Initialize Inputs
   !u(1) = time at n=1  (t=   0)
   !u(2) = time at n=0  (t= -dt)
   !u(3) = time at n=-1 (t= -2dt) if NumInp > 2
   if ( dvr%p%SimMod == 3 ) then
      ! General inputs
      do iu = 1, NumInp !u(NumInp) is overwritten in time-sim loop, so no need to init here 
         dvr%uTimes(iu) = (2-iu-1)*dvr%p%dt
      enddo
      ! Inflow "inputs"
      do iu = 1,NumInp
         call setInflow(t=dvr%uTimes(iu), p=dvr%p, U0=dvr%U0(iu,:))
      enddo
      ! UA inputs at t=0, stored in u(1)
      do iu = 1, NumInp-1 !u(NumInp) is overwritten in time-sim loop, so no need to init here 
         call setUAinputs(dvr%U0(iu,:), dvr%LD_x, dvr%p, dvr%m, dvr%UA_u(iu))
      enddo
   else
      ! UA inputs at t=0, stored in u(1)
      do iu = 1, NumInp-1 !u(NumInp) is overwritten in time-sim loop, so no need to init here 
         call setUAinputsAlphaSim(2-iu,  dvr%UA_u(iu), dvr%uTimes(iu), dvr%p, errStat, errMsg); call checkError()
      end do
   endif

   i = 1 ! nodes per blade
   j = 1 ! number of blades
   ! --- Time marching loop

   if ( dvr%p%SimMod == 3 ) then

      call Dvr_InitializeOutputs(dvr%out, dvr%p%numSteps, errStat, errMsg)

      dvr%LD_u(1)%Fext=0.0_ReKi ! TODO TODO
      dvr%LD_u(2)%Fext=0.0_ReKi ! TODO TODO

      ! --- time marching loop
      call WrScr(' Time simulation - TMax = '//trim(num2lstr(dvr%p%numSteps*dvr%p%dt)))
      do n = 1, dvr%p%numSteps

         ! --- Set inputs at t by storing in u(2) what was in u(1) at previous time step
         !u(1) = time at n=n+1  (t=t+dt)
         !u(2) = time at n=n    (t=t   )
         do iu = NumInp-1, 1, -1
            dvr%uTimes(iu+1)  = dvr%uTimes(iu)
            dvr%U0(    iu+1,:)= dvr%U0(iu,:)
            dvr%UA_u(  iu+1)  = dvr%UA_u(  iu)
            dvr%LD_u(  iu+1)  = dvr%LD_u(  iu)
         end do

         ! --- Calc Outputs at t 
         iu = 2 ! Index 2 is t
         t = dvr%uTimes(iu)
         ! Use existing states to compute the outputs
         call LD_CalcOutput(t, dvr%LD_u(iu), dvr%LD_p, dvr%LD_x, dvr%LD_xd, dvr%LD_z, dvr%LD_OtherState, dvr%LD_y, dvr%LD_m, errStat, errMsg); call checkError()
         !! Use existing states to compute the outputs
         call UA_CalcOutput(i, j, t, dvr%UA_u(iu),  dvr%UA_p, dvr%UA_x, dvr%UA_xd, dvr%UA_OtherState, dvr%AFI_Params(dvr%AFIndx(i,j)), dvr%UA_y, dvr%UA_m, errStat, errMsg ); call checkError()
         ! Generate file outputs
         call UA_WriteOutputToFile(t, dvr%UA_p, dvr%UA_y)
         ! Driver outputs
         call AeroKinetics(dvr%U0(iu,:), dvr%LD_x%q(1:3), dvr%LD_x%q(4:6), (/dvr%UA_y%Cl, dvr%UA_y%Cd, dvr%UA_y%Cm/), dvr%p, dvr%m)
         ! Write/Store outputs 
         call Dvr_WriteOutputs(n, t, dvr, dvr%out, errStat, errMsg); call checkError()

         ! --- Set inputs at t+dt in u(1)
         iu = 1 ! Index 1 is t+dt
         ! Basic inputs
         dvr%uTimes(iu) = (n+1-1)*dvr%p%dt
         ! Inflow inputs
         call setInflow(t=dvr%uTimes(iu), p=dvr%p, U0=dvr%U0(iu,:))
         ! LinDyn inputs at t+dt
         ! Using everything at t!!!! Dynamics are deterministic (only a function of values at t)
         ! NOTE: should use extrap interp maybe
         call setLDinputs(dvr%U0(2,:), dvr%LD_x, dvr%UA_y, dvr%p, dvr%m, dvr%LD_u(iu))

         ! --- Integrate LinDyn from t to t+dt
         call LD_UpdateStates(t, n, dvr%LD_u, dvr%uTimes, dvr%LD_p, dvr%LD_x, dvr%LD_xd, dvr%LD_z, dvr%LD_OtherState, dvr%LD_m, errStat, errMsg); call checkError()

         ! Calc LinDyn outputs at t+dt
         iu = 1 ! Index 1 is t+dt
         call LD_CalcOutput(t, dvr%LD_u(iu), dvr%LD_p, dvr%LD_x, dvr%LD_xd, dvr%LD_z, dvr%LD_OtherState, dvr%LD_y, dvr%LD_m, errStat, errMsg); call checkError()

         ! --- Set UA Inputs at t+dt
         call setUAinputs(dvr%U0(iu,:), dvr%LD_x, dvr%p, dvr%m, dvr%UA_u(iu))

         ! --- Integrate LinDyn from t to t+dt
         call UA_UpdateStates(i, j, t, n, dvr%UA_u, dvr%uTimes, dvr%UA_p, dvr%UA_x, dvr%UA_xd, dvr%UA_OtherState, dvr%AFI_Params(dvr%AFIndx(i,j)), dvr%UA_m, errStat, errMsg ); call checkError()
      end do

      call Dvr_EndSim(dvr, errStat, errMsg)

   else
      ! --- time marching loop
      do n = 1, dvr%p%numSteps
        
         ! set inputs:
         DO iu = NumInp-1, 1, -1
            dvr%UA_u(  iu+1) = dvr%UA_u(     iu)
            dvr%uTimes(iu+1) = dvr%uTimes(iu)
         END DO
     
         ! first value of uTimes/u contain inputs at t+dt
         call setUAinputsAlphaSim(n+1,  dvr%UA_u(1), dvr%uTimes(1), dvr%p, errStat, errMsg); call checkError()
           
         t = dvr%uTimes(2)

         ! Use existing states to compute the outputs
         call UA_CalcOutput(i, j, t, dvr%UA_u(2), dvr%UA_p, dvr%UA_x, dvr%UA_xd, dvr%UA_OtherState, dvr%AFI_Params(dvr%AFIndx(i,j)), dvr%UA_y, dvr%UA_m, errStat, errMsg ); call checkError()
               
         ! Generate file outputs
         call UA_WriteOutputToFile(t, dvr%UA_p, dvr%UA_y)
         
         ! Prepare states for next time step
         call UA_UpdateStates(i, j, t, n, dvr%UA_u, dvr%uTimes, dvr%UA_p, dvr%UA_x, dvr%UA_xd, dvr%UA_OtherState, dvr%AFI_Params(dvr%AFIndx(i,j)), dvr%UA_m, errStat, errMsg ); call checkError()
         
      end do
   endif
   
   ! --- Exit
   call Cleanup()
   call NormStop()

contains
   
   !====================================================================================================
   subroutine Cleanup()
      call UA_End(dvr%UA_p)
      ! probably should also deallocate driver variables here...
      
   end subroutine Cleanup

   !----------------------------------------------------------------------------------------------------  
   subroutine checkError()
      
      if (ErrStat >= AbortErrLev) then
         
         call Cleanup()
         call ProgAbort(ErrMsg)
            
      elseif ( ErrStat /= ErrID_None ) then
         
         call WrScr( trim(ErrMsg) )
            
      end if
      
   end subroutine checkError
   !----------------------------------------------------------------------------------------------------  
   
   subroutine print_help()
    print '(a)', 'usage: '
    print '(a)', ''
    print '(a)', 'UnsteadyAero_Driver.exe [driverfilename]'
    print '(a)', ''
    print '(a)', 'Where the optional argument, driverfilename, is the name of the UnsteadyAero driver input file.'
    print '(a)', ''

   end subroutine print_help
   
end program UnsteadyAero_Driver

