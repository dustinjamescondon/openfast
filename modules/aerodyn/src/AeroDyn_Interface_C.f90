!> This file of is used to export a DLL of Aerodyn_Interface. It aliases each of the exported subroutines so that the 
!! calling C/C++ subroutine just has to use the alias, without dealing with how the subroutine symbol is mangled by the 
!! compiling process. Also it deals with accepting a C string as an argument by including the length of the string as 
!! another argument, so FORTRAN knows its size. This has to be done because FORTRAN does not follow the null-terminated 
!! character procedure that C/C++ does.
      
   module Aerodyn_Interface_C

   USE AeroDyn_Interface_Types
   USE AeroDyn_Interface_Subs  ! this module holds the data needed for the simulation and loads the lower-level modules

   implicit none
   
   type, private :: SimInstance
      type(AeroDyn_Data)                 :: AD
      type(AeroDyn_Data)                 :: AD_saved ! Calling Interface_SaveCurrentStates does this: "AD_saved = AD"
      type(Dvr_SimData)                  :: DvrData
      logical(kind=C_BOOL)               :: ADInstance_Initialized

   end type SimInstance
      
   CONTAINS

   !---------------------------------------------------------------------------------------------------------------------------------------------------------------
   !> Initialize Aerodyn, which includes reading input files for number of blades and nodes.
   !! However, unlike Aerodyn_Driver, we don't set the initial inflows here. This is because
   !! inflows are set from an outside source, and that source needs to know the number of blades
   !! and nodes before sending inflows. Therefore, inflows must be set after this via initInflows
   !! once we know the number of nodes and blades.
   subroutine Interface_InitAeroDyn_C(AD_inputFile,inputFile_len,AD_outputFile,outputFile_len,timestep,useAddedMass,coeffAddedMass,fluidDensity,kinematicFluidVisc,&
      hubPos,hubOri,hubVel,hubAcc,hubRotVel,hubRotAcc,nBlades,bladePitch,hubRadius,precone,nNodes,turbineDiameter,instancePtr,&
      errStat_out, errMsg_out) !BIND(C, NAME="INTERFACE_INITAERODYN")
      !DEC$ ATTRIBUTES DLLEXPORT::Interface_InitAerodyn_C
      !GCC$ ATTRIBUTES DLLEXPORT::Interface_InitAerodyn_C
      !DEC$ ATTRIBUTES DECORATE, ALIAS: "INTERFACE_INITAERODYN":: Interface_InitAerodyn_C
      use, intrinsic :: ISO_C_BINDING, ONLY: C_CHAR, C_NULL_CHAR, C_PTR, C_LOC
      
      integer(C_INT),                      intent(in   )  :: inputFile_len
      character(inputFile_len, kind=C_CHAR),   intent(in   )  :: AD_inputFile
      integer(C_INT),                      intent(in   )  :: outputFile_len
      character(outputFile_len, kind=C_CHAR), intent(in   )  :: AD_outputFile
      real(C_DOUBLE),                      intent(in   )  :: timestep
      logical(kind=C_BOOL),                intent(in   )  :: useAddedMass
      real(C_DOUBLE),                      intent(in   )  :: coeffAddedMass ! coefficient of added mass (only used when useAddedMass=true)
      real(C_DOUBLE),                      intent(in   )  :: fluidDensity
      real(C_DOUBLE),                      intent(in   )  :: kinematicFluidVisc
      real(C_DOUBLE), dimension(1:3),      intent(in   )  :: hubPos
      real(C_DOUBLE), dimension(1:3,1:3),  intent(in   )  :: hubOri ! The hub's global to local orientation matrix
      real(C_DOUBLE), dimension(1:3),      intent(in   )  :: hubVel
      real(C_DOUBLE), dimension(1:3),      intent(in   )  :: hubAcc
      real(C_DOUBLE), dimension(1:3),      intent(in   )  :: hubRotVel
      real(C_DOUBLE), dimension(1:3),      intent(in   )  :: hubRotAcc
      integer(C_INT),                      intent(in   )  :: nBlades
      real(C_DOUBLE),                      intent(in   )  :: bladePitch
      real(C_DOUBLE),                      intent(in   )  :: hubRadius
      real(C_DOUBLE),                      intent(in   )  :: precone
      integer(C_INT),                      intent(  out)  :: nNodes
      real(C_DOUBLE),                      intent(  out)  :: turbineDiameter
      type(C_PTR),                         intent(  out)  :: instancePtr
      integer(C_INT),                      intent(  out)  :: ErrStat_out
      character(ErrMsgLen + 1,kind=C_CHAR),intent(  out)  :: ErrMsg_out

      character(*), parameter                     :: RoutineName = 'Interface_InitAeroDyn_C'
      real(ReKi)    :: rmax    ! storing the maximum radius of a blade node relative to hub center
      integer       :: j       ! for iterating through number of nodes on one blade

      integer(IntKi)             :: ErrStat  ! Status of error message
      character(ErrMsgLen)	      :: ErrMsg	! Error message if ErrStat /= ErrID_None
      type(SimInstance), pointer :: simIns   ! Pointer to dynamically allocated Simulation Instance 
      
      allocate(simIns)
      
      ! assigns instance pointer the memory address of simIns
      instancePtr = C_LOC(simIns)
      
      simIns%ADInstance_Initialized = .false.

      ! Call this function to validate the given parameters   
      call Interface_Init_Parameters(AD_InputFile,AD_outputFile,timestep,nBlades,hubRadius,precone,simIns%DvrData,ErrStat,ErrMsg)
      !CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )   ! update error list if applicable
      if ( NeedToAbort(ErrStat) ) then
         return
      end if
   
      ! Initialize AeroDyn
      ! Set the Initialization input data for AeroDyn based on the Driver input file data, and initialize AD
      ! (this also initializes inputs to AD for first time step)
      call Init_AeroDyn(simIns%DvrData,simIns%AD,hubPos,hubOri,hubVel,hubAcc,hubRotVel,hubRotAcc,bladePitch,simIns%DvrData%dt,useAddedMass,&
         real(coeffAddedMass, kind=DbKi),errStat, errMsg)
      CALL SetErrStat( ErrStat, ErrMsg, ErrStat_out, ErrMsg_out, RoutineName )   ! update error list if applicable
      if( NeedToAbort(ErrStat) ) then
         return
      end if
   
      ! Set these parameters from the arguments, which will overwrite what was specified in the input file
      !  Note: Setting these from subroutine arguments makes sense with the coupling because PDS can just set these according to
      !        what its own simulation is using, rather than having to update the AD input file to match them.
      ! TODO: NOTE I COMMENTED THESE OUT BECAUSE RIGHT NOW I'M ASSUMING THESE ARE SET FROM THE AD INPUT FILE
      !simIns%AD%p%AirDens = fluidDensity
      !simIns%AD%p%KinVisc = kinematicFluidVisc

      simIns%ADInstance_Initialized = .true.

      nNodes = simIns%AD%p%NumBlNds   

      ! Initialize the main output file
      call Dvr_InitializeOutputFile( simIns%DvrData%OutFileData, errStat, errMsg)
      CALL SetErrStat( ErrStat, ErrMsg, ErrStat_out, ErrMsg_out, RoutineName )   ! update error list if applicable
      if ( NeedToAbort(ErrStat) ) then
         return
      endif

      turbineDiameter = CalcTurbineDiameter(simIns%AD, nNodes)
   
      if ( NeedToAbort(errStat) ) then
         errStat_out = errStat
         ! set the c string reference so the calling function will have
         errMsg_out = TRIM(errMsg) // C_NULL_CHAR
         return 
      endif
      
      ! Allocate and initialize our saved AeroDyn states
      call AD_Dvr_CopyAeroDyn_Data(simIns%AD, simIns%AD_saved, MESH_NEWCOPY, ErrStat, ErrMsg)

      contains
      function CalcTurbineDiameter(AD,nNodes) result(turbineDiameter)
         type(AeroDyn_Data), intent(in) :: AD
         integer(IntKi),     intent(in) :: nNodes
      
         real(ReKi)                     :: turbineDiameter
         real(ReKi)                     :: rmax
         integer(IntKi)                 :: i
      
            ! I know the first blade is always pointing up and already "pitched" to the appropriate precone angle, 
            ! so just find the highest z value on that blade
            rmax = 0.0_ReKi
            do i=1,nNodes
               rmax = max(rmax, AD%u(1)%BladeMotion(1)%Position(3,i) )
            end do !i=nodes
            turbineDiameter = rmax * 2.0_ReKi
      end function CalcTurbineDiameter

   !------------------------------------------------------------------
   END SUBROUTINE Interface_InitAeroDyn_C


   !---------------------------------------------------------------------------------------------------------------------------------------------------------------
   SUBROUTINE Interface_InitInputs_Inflow_C(simInsAddr, nBlades, nNodes, bladeNodeInflowVel, bladeNodeInflowAcc) BIND(C, NAME='INTERFACE_INITINPUTS_INFLOW')
      !DEC$ ATTRIBUTES DLLEXPORT :: Interface_InitInputs_Inflow_C
      !GCC$ ATTRIBUTES DLLEXPORT :: Interface_InitInputs_Inflow_C

      type(C_PTR), intent(in), value:: simInsAddr
      integer(C_INT), intent(in) :: nBlades, nNodes
      real(C_DOUBLE), dimension(3,nNodes,nBlades), intent(in) :: bladeNodeInflowVel
      real(C_DOUBLE), dimension(3,nNodes,nBlades), intent(in) :: bladeNodeInflowAcc
      type(SimInstance), pointer :: v
      
      call C_F_POINTER(simInsAddr, v)

      call Init_AD_Inputs_Inflow(v%AD, nBlades, nNodes, bladeNodeInflowVel, bladeNodeInflowAcc)

   END SUBROUTINE Interface_InitInputs_Inflow_C
   
   !SUBROUTINE Interface_SaveCurrentStates(simInsAddr, statePtr_out) BIND(C, NAME="INTERFACE_SAVECURRENTADSTATES")
   
   
   !---------------------------------------------------------------------------------------------------------------------------------------------------------------
   SUBROUTINE Interface_SaveCurrentStates(simInsAddr) BIND(C, NAME="INTERFACE_SAVECURRENTSTATES")
      !DEC$ ATTRIBUTES DLLEXPORT::Interface_SaveCurrentStates
      !GCC$ ATTRIBUTES DLLEXPORT::Interface_SaveCurrentStates
      
      use, intrinsic :: ISO_C_BINDING, ONLY: C_PTR, C_DOUBLE
      
      type(C_PTR),   value,               intent(in) :: simInsAddr
      type(SimInstance), pointer                 :: v
      
      integer(IntKi)             :: ErrStat  ! Status of error message
      character(ErrMsgLen)	      :: ErrMsg	! Error message if ErrStat /= ErrID_None
      
      call C_F_POINTER(simInsAddr, v)
      call AD_Dvr_CopyAeroDyn_Data(v%AD, v%AD_saved, MESH_UPDATECOPY, ErrStat, ErrMsg)

   END SUBROUTINE Interface_SaveCurrentStates

   !---------------------------------------------------------------------------------------------------------------------------------------------------------------
   SUBROUTINE Interface_RestoreSavedStates(simInsAddr) BIND(C, NAME="INTERFACE_RESTORESAVEDSTATES")
      !DEC$ ATTRIBUTES DLLEXPORT::Interface_RestoreSavedStates
      !GCC$ ATTRIBUTES DLLEXPORT::Interface_RestoreSavedStates
      
      use, intrinsic :: ISO_C_BINDING, ONLY: C_PTR, C_DOUBLE
      
      type(C_PTR),   value,               intent(in) :: simInsAddr
      type(SimInstance), pointer                 :: v
      
      integer(IntKi)             :: ErrStat  ! Status of error message
      character(ErrMsgLen)	      :: ErrMsg	! Error message if ErrStat /= ErrID_None
      
      call C_F_POINTER(simInsAddr, v)
      call AD_Dvr_CopyAeroDyn_Data(v%AD_saved, v%AD, MESH_UPDATECOPY, ErrStat, ErrMsg)

   END SUBROUTINE Interface_RestoreSavedStates
   
   !---------------------------------------------------------------------------------------------------------------------------------------------------------------
   SUBROUTINE Interface_SetInputs_Hub_C(simInsAddr, time, hubPos, hubOri, hubVel, hubAcc, hubRotVel, hubRotAcc, bladePitch) BIND(C, NAME='INTERFACE_SETINPUTS_HUB')
      !DEC$ ATTRIBUTES DLLEXPORT::Interface_SetInputs_Hub_C
      !GCC$ ATTRIBUTES DLLEXPORT::Interface_SetInputs_Hub_C

      use, intrinsic :: ISO_C_BINDING, ONLY: C_PTR, C_DOUBLE


      type(C_PTR),   value,               intent(in) :: simInsAddr
      type(SimInstance), pointer                 :: v
      real(C_DOUBLE),	                  intent(in) :: time              ! time of these inputs
      real(C_DOUBLE), dimension(1:3),     intent(in) :: hubPos
      real(C_DOUBLE), dimension(1:3,1:3), intent(in) :: hubOri
      real(C_DOUBLE), dimension(1:3),     intent(in) :: hubVel
      real(C_DOUBLE), dimension(1:3),     intent(in) :: hubAcc
      real(C_DOUBLE), dimension(1:3),     intent(in) :: hubRotVel
      real(C_DOUBLE), dimension(1:3),     intent(in) :: hubRotAcc
      real(C_DOUBLE),                     intent(in) :: bladePitch
      integer(C_INT)                        :: ErrStat
      character(ErrMsgLen + 1,kind=C_CHAR)  :: ErrMsg
      
      call C_F_POINTER(simInsAddr, v)
      
      call Set_AD_Inputs_Hub(real(time, kind=DbKi),v%DvrData,v%AD%u(1),hubPos,hubOri,hubVel,hubAcc,hubRotVel,hubRotAcc,bladePitch)
      v%AD%InputTime(1) = real(time, kind=DbKi)
         
   END SUBROUTINE Interface_SetInputs_Hub_C
   
   !----------------------------------------------------------------------------------
   !> 
   SUBROUTINE Interface_SetInputs_Inflow_C(simInsAddr, nBlades, nNodes, bladeNodeInflowVel, bladeNodeInflowAcc) BIND(C, NAME='INTERFACE_SETINPUTS_INFLOW')
      !DEC$ ATTRIBUTES DLLEXPORT::Interface_SetInputs_Inflow_C
      !GCC$ ATTRIBUTES DLLEXPORT::Interface_SetInputs_Inflow_C
      use, intrinsic :: ISO_C_BINDING, ONLY: C_PTR, C_DOUBLE

      type(C_PTR),   value,                     intent(in)    :: simInsAddr
      type(SimInstance), pointer                              :: v
      integer(C_INT),                           intent(in)    :: nBlades, nNodes
      real(C_DOUBLE), dimension(3,nNodes,nBlades),  intent(in)    :: bladeNodeInflowVel
      real(C_DOUBLE), dimension(3,nNodes,nBlades),  intent(in)    :: bladeNodeInflowAcc

      call C_F_POINTER(simInsAddr, v)
            
      call Set_AD_Inputs_Inflow(v%AD, nBlades, nNodes, bladeNodeInflowVel, bladeNodeInflowAcc)

   END SUBROUTINE Interface_SetInputs_Inflow_C
      
   !----------------------------------------------------------------------------------
   !> 
   SUBROUTINE Interface_Advance_InputWindow_C(simInsAddr) BIND(C, NAME='INTERFACE_ADVANCE_INPUTWINDOW')
      !DEC$ ATTRIBUTES DLLEXPORT:: Interface_Advance_InputWindow_C
      !GCC$ ATTRIBUTES DLLEXPORT:: Interface_Advance_InputWindow_C
      use, intrinsic :: ISO_C_BINDING, ONLY: C_PTR, C_DOUBLE

      type(C_PTR),   value,                     intent(in   ) :: simInsAddr
      type(SimInstance), pointer                              :: v
      
      integer(IntKi)             :: ErrStat  ! Status of error message
      character(ErrMsgLen)	      :: ErrMsg	! Error message if ErrStat /= ErrID_None

      call C_F_POINTER(simInsAddr, v)

      call Advance_AD_InputWindow(v%AD, errStat, errMsg)
      
      ! TODO do error checking
      
   END SUBROUTINE Interface_Advance_InputWindow_C
   
   !----------------------------------------------------------------------------------------------------------------------------------
   !< AeroDyn (with added mass) has two inputs that are accelerations and are therefore part of the direct feedthrough problem: linear acceleration
   !< and angular acceleration, and both have 3 components; therefore, there are 6 inputs to perturb in total
   subroutine Interface_SetInputs_HubAcceleration( simInsAddr, linearAcc, rotationAcc ) BIND(C, NAME='INTERFACE_SETINPUTS_HUBACCELERATION')
      !DEC$ ATTRIBUTES DLLEXPORT::Interface_SetInputs_HubAcceleration
      !GCC$ ATTRIBUTES DLLEXPORT::Interface_SetInputs_HubAcceleration
      use, intrinsic :: ISO_C_BINDING, ONLY: C_PTR, C_DOUBLE
   !..................................................................................................................................
      type(C_PTR), value,             intent(in   )   :: simInsAddr
      type(SimInstance), pointer                      :: v
      real(C_DOUBLE),         intent(in   ) :: linearAcc(3)    !< Hub's linear acceleration
      real(C_DOUBLE),         intent(in   ) :: rotationAcc(3)  !< Hub's rotational acceleration (axis-angle vector)
      
      call C_F_POINTER(simInsAddr, v)
      
      call Calc_AD_NodeKinematics(linearAcc, rotationAcc, v%AD%u(1), v%DvrData)

   end subroutine Interface_SetInputs_HubAcceleration
   
   !----------------------------------------------------------------------------------
   !< 
   SUBROUTINE Interface_CalcOutput_C(simInsAddr, force, moment, power, tsr) BIND(C, NAME='INTERFACE_CALCOUTPUT')
      !DEC$ ATTRIBUTES DLLEXPORT::Interface_CalcOutput_C
      !GCC$ ATTRIBUTES DLLEXPORT::Interface_CalcOutput_C
      use, intrinsic :: ISO_C_BINDING, ONLY: C_PTR, C_DOUBLE
      use AeroDyn_IO, only: RtAeroPwr, RtAeroFxh, RtAeroFyh, RtAeroFzh, RtAeroMxh, RtAeroMyh, RtAeroMzh, RtTSR

      
      type(C_PTR), value,             intent(in   )   :: simInsAddr
      type(SimInstance), pointer                      :: v
      real(C_DOUBLE), dimension(1:3), intent(  out)	:: force     ! 3D force (Fx, Fy, Fz)
      real(C_DOUBLE), dimension(1:3), intent(  out)   :: moment    ! moment (Mx, My, Mz) at rotor hub
      real(C_DOUBLE),                 intent(  out)   :: power
      real(C_DOUBLE),                 intent(  out)   :: tsr

      
      real(DbKi)                                  :: time_ad
      integer(IntKi)             :: ErrStat  ! Status of error message 
      character(ErrMsgLen)	      :: ErrMsg	! Error message if ErrStat /= ErrID_None
      
      call C_F_POINTER(simInsAddr, v)
            

      time_ad = v%AD%stepNum * v%DvrData%dt

      ! This will update the AllOuts in the main AeroDyn output states
      call AD_CalcOutput( time_ad, v%AD%u(1), v%AD%p, v%AD%x(STATE_PRED), v%AD%xd(STATE_PRED),&
         v%AD%z(STATE_PRED), v%AD%OtherState(STATE_PRED), v%AD%y, v%AD%m, errStat, errMsg )

      ! So then set the output arguments accordingly
      force(1) =  v%AD%m%AllOuts( RtAeroFxh )
      force(2) =  v%AD%m%AllOuts( RtAeroFyh )
      force(3) =  v%AD%m%AllOuts( RtAeroFzh )

      moment(1) = v%AD%m%AllOuts( RtAeroMxh )
      moment(2) = v%AD%m%AllOuts( RtAeroMyh )
      moment(3) = v%AD%m%AllOuts( RtAeroMzh )
      
      power = v%AD%m%AllOuts( RtAeroPwr )
      tsr = v%AD%m%AllOuts (RtTSR)
      
   END SUBROUTINE Interface_CalcOutput_C
   
   !-------------------------------------------------------------------------------------------------------
   !> 
   SUBROUTINE Interface_CopyStates_Pred_to_Curr(simInsAddr) BIND(C, NAME='INTERFACE_COPYSTATES_PRED_TO_CURR')
      !DEC$ ATTRIBUTES DLLEXPORT::Interface_CopyStates_Pred_to_Curr
      !GCC$ ATTRIBUTES DLLEXPORT::Interface_CopyStates_Pred_to_Curr
      type(C_PTR), value,             intent(in   )   :: simInsAddr
      type(SimInstance), pointer                      :: v
      integer(IntKi)             :: ErrStat2  ! Status of error message
      character(ErrMsgLen)	      :: ErrMsg2	! Error message if ErrStat /= ErrID_None

      call C_F_POINTER(simInsAddr, v)
      
      ! Copy current states to pred states
      CALL AD_CopyContState   (v%AD%x( STATE_PRED), v%AD%x( STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
      !   CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL AD_CopyDiscState   (v%AD%xd(STATE_PRED), v%AD%xd(STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)  
      !   CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL AD_CopyConstrState (v%AD%z( STATE_PRED), v%AD%z( STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
      !   CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL AD_CopyOtherState( v%AD%OtherState(STATE_PRED), v%AD%OtherState(STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
      !   CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   END SUBROUTINE Interface_CopyStates_Pred_to_Curr
   
   !---------------------------------------------------------------------------------------------------------------------------------------------------------------
   !> TODO remove the "contains" subroutine indirection. Only added this because of "isRealStep" but now I've removed it
   SUBROUTINE Interface_UpdateStates_C(simInsAddr) BIND(C, NAME='INTERFACE_UPDATESTATES')
      !DEC$ ATTRIBUTES DLLEXPORT::Interface_UpdateStates_C
      !GCC$ ATTRIBUTES DLLEXPORT::Interface_UpdateStates_C

      type(C_PTR), value,             intent(in   )      :: simInsAddr
      type(SimInstance), pointer                      :: v

      integer(IntKi)		        :: errStat  ! Status of error message
      character(ErrMsgLen)	     :: errMsg   ! Error message if ErrStat /= ErrID_None

      call C_F_POINTER(simInsAddr, v)
      
      call Interface_UpdateStates(v%AD, errStat, errMsg)
      
   contains    
   SUBROUTINE Interface_UpdateStates(AD, errStat, errMsg)

   type(AeroDyn_Data),             intent(inout) :: AD

   integer(IntKi),                          intent(  out)   :: errStat              ! Status of error message
   character(ErrMsgLen),                    intent(  out)   :: errMsg               ! Error message if ErrStat /= ErrID_None

   ! loop counters
   integer(IntKi) 	::	i
   integer(IntKi) 	::	j
   integer(IntKi) 	::	k

   real(ReKi)                                  :: theta(3)
   real(ReKi)                                  :: position(3)
   real(ReKi)                                  :: orientation(3,3)
   real(ReKi)                                  :: rotateMat(3,3)

   real(DbKi), parameter                       :: epsilon = 1e-6

   INTEGER(IntKi)                         :: ErrStat2                     ! error status from any called subroutines
   CHARACTER(ErrMsgLen)                   :: ErrMsg2                      ! error message from any called subroutines
   character(*), parameter                :: RoutineName = 'Interface_Simulate'   ! name of this subroutine for error printing

   errStat     = ErrID_None
   errMsg      = ''
      
   ! SEE ED_AllocOutput in ElastoDyn.f90 for examples of setting up meshes for rotor
   
   ! Copy current states to pred states
   CALL AD_CopyContState   (AD%x( STATE_CURR), AD%x( STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2)
      CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AD_CopyDiscState   (AD%xd(STATE_CURR), AD%xd(STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2)  
      CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AD_CopyConstrState (AD%z( STATE_CURR), AD%z( STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2)
      CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AD_CopyOtherState( AD%OtherState(STATE_CURR), AD%OtherState(STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2)
      CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   !.............................
   ! Update states from t_i to t_(i+1)

   ! Then call UpdateStates to bring the states from "time" to "time + dt"
   call AD_UpdateStates(AD%InputTime(2), AD%stepNum, AD%u, AD%InputTime, AD%p, AD%x(STATE_PRED), &
      AD%xd(STATE_PRED), AD%z(STATE_PRED), AD%OtherState(STATE_PRED), AD%m, errStat2, errMsg2 )
   call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
   
   END SUBROUTINE Interface_UpdateStates

   END SUBROUTINE Interface_UpdateStates_C

   !---------------------------------------------------------------------------------------------------------------------------------------------------------------
   !> 
   SUBROUTINE Interface_WriteOutputLine(simInsAddr) BIND(C, NAME='INTERFACE_WRITEOUTPUTLINE')
   !DEC$ ATTRIBUTES DLLEXPORT::Interface_WriteOutputLine
   !GCC$ ATTRIBUTES DLLEXPORT::Interface_WriteOutputLine
   type(C_PTR), value,             intent(in   )      :: simInsAddr
   type(SimInstance), pointer                      :: v
   
   integer(IntKi) :: ErrStat
   character(ErrMsgLen) :: ErrMsg
   
      call C_F_POINTER(simInsAddr, v)
   
      call Dvr_WriteOutputLine(v%DvrData%OutFileData, v%AD%InputTime(1), v%AD%y%WriteOutput, errStat, errMsg)
      !call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
   
   END SUBROUTINE Interface_WriteOutputLine
   
   !---------------------------------------------------------------------------------------------------------------------------------------------------------------
   !> TODO remove the "contains" subroutine indirection. Only added this because of "isRealStep" but now I've removed it   
   SUBROUTINE Interface_GetBladeNodePos_C(simInsAddr, nodePos) BIND(C, NAME='INTERFACE_GETBLADENODEPOS')
   !DEC$ ATTRIBUTES DLLEXPORT::Interface_GetBladeNodePos_C
   !GCC$ ATTRIBUTES DLLEXPORT::Interface_GetBladeNodePos_C
      use, intrinsic :: ISO_C_BINDING, ONLY: C_PTR, C_F_POINTER

      type(C_PTR),     intent(in ), value           :: simInsAddr
      type(SimInstance), pointer                    :: v
      real(C_DOUBLE), dimension(1:*), intent(  out) :: nodePos	! variable for outgoing node positions

      integer(C_INT)                              :: errStat   ! Status of error message
      character(ErrMsgLen)		                    :: errMsg    ! Error message if ErrStat /= ErrID_None

      call C_F_POINTER(simInsAddr, v)
      
      ! Access AD internals and calculate positions, passing array reference backwards
      call Interface_GetBladeNodePos(v%AD, v%DvrData, nodePos)
   
   contains
      ! gets the node coordinates of all the blade nodes from AeroDyn
   SUBROUTINE Interface_GetBladeNodePos(AD, DvrData, bladeNodePos)

   type(AeroDyn_Data),        intent(inout) :: AD
   type(Dvr_SimData),         intent(in   ) :: DvrData
   real(R8Ki), dimension(1:*), intent(out)  :: bladeNodePos	!output variable for blade node positions

   integer(intKi)		:: j             ! loop counter for nodes
   integer(intKi)		:: k             ! loop counter for blades

   real(ReKi)      :: nodePos(3)
   integer(8)		:: curIndex

   ! Blade motions:
   do k=1,DvrData%numBlades

      do j=1,AD%u(1)%BladeMotion(k)%nnodes

         curIndex = ((k-1) * AD%u(1)%BladeMotion(k)%nnodes * 3)  + ((j-1) * 3)

         nodePos = AD%u(1)%BladeMotion(k)%Position(:,j) + AD%u(1)%BladeMotion(k)%TranslationDisp(:,j)

         bladeNodePos(curIndex+1) = nodePos(1)
         bladeNodePos(curIndex+2) = nodePos(2)
         bladeNodePos(curIndex+3) = nodePos(3)

      end do !j=nnodes

   end do !k=numBlades

   END SUBROUTINE Interface_GetBladeNodePos

   END SUBROUTINE Interface_GetBladeNodePos_C

   !---------------------------------------------------------------------------------------------------------------------------------------------------------------
   SUBROUTINE Interface_End_C(simInsAddr) BIND(C, NAME='INTERFACE_END')
      !DEC$ ATTRIBUTES DLLEXPORT::Interface_End_C
      !GCC$ ATTRIBUTES DLLEXPORT::Interface_End_C

      use, intrinsic :: ISO_C_BINDING, ONLY:  C_PTR, C_F_POINTER

      type(C_PTR), intent(in), value           :: simInsAddr
      type(SimInstance), pointer               :: v ! the memory pointed to by p
      integer(IntKi)                           :: ErrStat	! Status of error message
      character(ErrMsgLen)	                    :: ErrMsg	! Error message if ErrStat /= ErrID_None
      character(*), parameter                  :: RoutineName = 'Interface_End_C'

      call C_F_POINTER(simInsAddr, v)

      if ( ASSOCIATED(v) ) then

          ! Close the output file
         if (v%DvrData%OutFileData%unOutFile > 0) close(v%DvrData%OutFileData%unOutFile)

         if ( v%ADInstance_Initialized ) then	! (?) Do we need to call AD_End? Doesn't AD_Dvr_DestroyDvr_SimData do all we need?
            ! TODO confirm that we don't need to call this (comment out for now)
            ! AD_End( v%AD%u(1), v%AD%p, V%AD%x(STATE_CURR), v%AD%xd(STATE_CURR), v%AD%z(STATE_CURR), v%AD%OtherState(STATE_CURR), v%AD%y, v%AD%m, errStat, errMsg )
               !call SetErrStat( errStat, errMsg, errStat_out, errMsg_out, RoutineName )
         end if

         call AD_Dvr_DestroyDvr_SimData( v%DvrData, ErrStat, ErrMsg )
            !call SetErrStat( errStat, errMsg, errStat, errMsg, RoutineName )

         call AD_Dvr_DestroyAeroDyn_Data( v%AD, ErrStat, ErrMsg )
            !call SetErrStat( errStat, errMsg, errStat_out, errMsg_out, RoutineName )
   
         call AD_Dvr_DestroyAeroDyn_Data( v%AD_saved, ErrStat, ErrMsg )
            !call SetErrStat( errStat, errMsg, errStat_out, errMsg_out, RoutineName )


         if (ErrStat >= AbortErrLev) then
            !CALL ProgAbort( 'AeroDyn Driver encountered simulation error level: '&
            !//TRIM(GetErrStr(ErrStat)), TrapErrors=.FALSE., TimeWait=3._ReKi )  ! wait 3 seconds (in case they double-clicked and got an error)
         else
            !call NormStop() ! don't call this because it will end the entire process that calls it
         end if
         
         deallocate(v)
      end if

      END SUBROUTINE Interface_End_C

   end module AeroDyn_Interface_C
