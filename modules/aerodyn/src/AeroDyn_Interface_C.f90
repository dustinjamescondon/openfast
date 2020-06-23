!> This file of is used to export a DLL of Aerodyn_Interface. It aliases each of the exported subroutines so that the 
!! calling C/C++ subroutine just has to use the alias, without dealing with how the subroutine symbol is mangled by the 
!! compiling process. Also it deals with accepting a C string as an argument by including the length of the string as 
!! another argument, so FORTRAN knows its size. This has to be done because FORTRAN does not follow the null-terminated 
!! character procedure that C/C++ does.
      
   module Aerodyn_Interface_C

   USE AeroDyn_Interface_Types
   USE AeroDyn_Interface  ! this module holds the data needed for the simulation and loads the lower-level modules

   implicit none
   
   PUBLIC :: Interface_InitAeroDyn_C		! Wrap initialization function
   PUBLIC :: Interface_UpdateStates_C		! Steps AD ahead one internal step
   PUBLIC :: Interface_GetBladeNodePos_C		! Pass blade node positions up to wrapper to go to PDS
   PUBLIC :: Interface_End_C 		! Destroy data structures used by AD and PDS, also generate output files?
   
   type, private :: SimInstance
      type(AeroDyn_Data)                 :: AD
      type(AeroDyn_Data)                 :: AD_fake  ! Instance of AeorDyn's states for a fake update which won't affect AD
      type(Dvr_SimData)                  :: DvrData
      logical                            :: ADInstance_Initialized

   end type SimInstance
   
   CONTAINS

   !---------------------------------------------------------------------------------------------------------------------------------------------------------------
   !> Initialize Aerodyn, which includes reading input files for number of blades and nodes.
   !! However, unlike Aerodyn_Driver, we don't set the initial inflows here. This is because
   !! inflows are set from an outside source, and that source needs to know the number of blades
   !! and nodes before sending inflows. Therefore, inflows must be set after this via initInflows
   !! once we know the number of nodes and blades.
   subroutine Interface_InitAeroDyn_C(AD_inputFile,inputFile_len,AD_outputFile,outputFile_len,timestep,useAddedMass,fluidDensity,kinematicFluidVisc,&
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
      logical,                             intent(in   )  :: useAddedMass
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

      integer(IntKi)             :: ErrStat  ! Status of error message
      character(ErrMsgLen)	      :: ErrMsg	! Error message if ErrStat /= ErrID_None
      type(SimInstance), pointer :: simIns   ! Pointer to dynamically allocated Simulation Instance 
      
      allocate(simIns)
      
      ! assigns instance pointer the memory address of simIns
      instancePtr = C_LOC(simIns)
      
      simIns%ADInstance_Initialized = .false.

      CALL Interface_Init(simIns%AD, simIns%DvrData, &
         simIns%ADInstance_Initialized,AD_inputFile,AD_outputFile,timestep,useAddedMass, &
         fluidDensity, kinematicFluidVisc, hubPos,hubOri,hubAcc,hubVel,         &
         hubRotVel,hubRotAcc,nBlades,bladePitch,hubRadius,precone,nNodes,turbineDiameter, errStat, errMsg)
      
      if ( NeedToAbort(errStat) ) then
         errStat_out = errStat
         ! set the c string reference so the calling function will have
         errMsg_out = TRIM(errMsg) // C_NULL_CHAR
         return 
      endif
      
   !------------------------------------------------------------------
   END SUBROUTINE Interface_InitAeroDyn_C


   !---------------------------------------------------------------------------------------------------------------------------------------------------------------
   SUBROUTINE Interface_InitInflows_C(simInsAddr, nBlades, nNodes, bladeNodeInflowVel, bladeNodeInflowAcc) BIND(C, NAME='INTERFACE_INITINFLOWS')
      !DEC$ ATTRIBUTES DLLEXPORT :: Interface_InitInflows_C
      !GCC$ ATTRIBUTES DLLEXPORT :: Interface_InitInflows_C

      type(C_PTR), intent(in), value:: simInsAddr
      integer(C_INT), intent(in) :: nBlades, nNodes
      real(C_DOUBLE), dimension(3,nNodes,nBlades), intent(in) :: bladeNodeInflowVel
      real(C_DOUBLE), dimension(3,nNodes,nBlades), intent(in) :: bladeNodeInflowAcc
      type(SimInstance), pointer :: v
      
      call C_F_POINTER(simInsAddr, v)

      call Init_AD_Inflows(v%AD, nBlades, nNodes, bladeNodeInflowVel, bladeNodeInflowAcc)

   END SUBROUTINE Interface_InitInflows_C
   
   
   !---------------------------------------------------------------------------------------------------------------------------------------------------------------
   SUBROUTINE Interface_SetHubMotion_C(simInsAddr, time, hubPos, hubOri, hubVel, hubAcc, hubRotVel, hubRotAcc, bladePitch) BIND(C, NAME='INTERFACE_SETHUBMOTION')
      !DEC$ ATTRIBUTES DLLEXPORT::Interface_SetHubMotion_C
      !GCC$ ATTRIBUTES DLLEXPORT::Interface_SetHubMotion_C

      use, intrinsic :: ISO_C_BINDING, ONLY: C_PTR, C_DOUBLE


      type(C_PTR),   value,    intent(in) :: simInsAddr
      type(SimInstance), pointer                 :: v
      real(C_DOUBLE),	                  intent(in) :: time              ! time of these inputs
      real(C_DOUBLE), dimension(1:3),     intent(in) :: hubPos
      real(C_DOUBLE), dimension(1:3,1:3), intent(in) :: hubOri
      real(C_DOUBLE), dimension(1:3),     intent(in) :: hubVel
      real(C_DOUBLE), dimension(1:3),     intent(in) :: hubAcc
      real(C_DOUBLE), dimension(1:3),     intent(in) :: hubRotVel
      real(C_DOUBLE), dimension(1:3),     intent(in) :: hubRotAcc
      real(C_DOUBLE),                     intent(in) :: bladePitch

      call C_F_POINTER(simInsAddr, v)
      
      call Interface_SetHubMotion(v%AD,v%DvrData,time,hubPos,hubOri,hubVel,hubAcc,hubRotVel,hubRotAcc,bladePitch)

   END SUBROUTINE Interface_SetHubMotion_C
   
   !---------------------------------------------------------------------------------------------------------------------------------------------------------------
   !> Sets the hub kinematics for a fake step, which will not affect the main AD states
   SUBROUTINE Interface_SetHubMotion_Fake_C(simInsAddr, time, hubPos, hubOri, hubVel, hubAcc, hubRotVel, hubRotAcc, bladePitch) BIND(C, NAME='INTERFACE_SETHUBMOTION_FAKE')
      !DEC$ ATTRIBUTES DLLEXPORT::Interface_SetHubMotion_Fake_C
      !GCC$ ATTRIBUTES DLLEXPORT::Interface_SetHubMotion_Fake_C

      use, intrinsic :: ISO_C_BINDING, ONLY: C_PTR, C_DOUBLE


      type(C_PTR),   value,    intent(in)            :: simInsAddr
      type(SimInstance), pointer                     :: v
      real(C_DOUBLE),	                  intent(in) :: time       ! time of these inputs
      real(C_DOUBLE), dimension(1:3),     intent(in) :: hubPos
      real(C_DOUBLE), dimension(1:3,1:3), intent(in) :: hubOri
      real(C_DOUBLE), dimension(1:3),     intent(in) :: hubVel
      real(C_DOUBLE), dimension(1:3),     intent(in) :: hubAcc
      real(C_DOUBLE), dimension(1:3),     intent(in) :: hubRotVel
      real(C_DOUBLE), dimension(1:3),     intent(in) :: hubRotAcc
      real(C_DOUBLE),                     intent(in) :: bladePitch
      
      integer(IntKi)                             :: ErrStat
      character(ErrMsgLen)                       :: ErrMsg

      call C_F_POINTER(simInsAddr, v)
      
      ! Copy all of our real AD states into the fake one
      call AD_Dvr_CopyAeroDyn_Data(v%AD, v%AD_fake, MESH_NEWCOPY, ErrStat, ErrMsg)
      
      call Interface_SetHubMotion(v%AD_fake,v%DvrData,time,hubPos,hubOri,hubVel,hubAcc,hubRotVel,hubRotAcc,bladePitch)

   END SUBROUTINE Interface_SetHubMotion_Fake_C
   
   !----------------------------------------------------------------------------------
   !> 
   SUBROUTINE Interface_SetInflows_C(simInsAddr, nBlades, nNodes, bladeNodeInflowVel, bladeNodeInflowAcc) BIND(C, NAME='INTERFACE_SETINFLOWS')
      !DEC$ ATTRIBUTES DLLEXPORT::Interface_SetInflows_C
      !GCC$ ATTRIBUTES DLLEXPORT::Interface_SetInflows_C
      use, intrinsic :: ISO_C_BINDING, ONLY: C_PTR, C_DOUBLE

      type(C_PTR),   value,                     intent(in)    :: simInsAddr
      type(SimInstance), pointer                              :: v
      integer(IntKi),                           intent(in)    :: nBlades, nNodes
      real(ReKi), dimension(3,nNodes,nBlades),  intent(in)    :: bladeNodeInflowVel
      real(ReKi), dimension(3,nNodes,nBlades),  intent(in)    :: bladeNodeInflowAcc

      ! Use to reshape the one-dimensional array of the blade node inflows

      call C_F_POINTER(simInsAddr, v)
            
      call Set_AD_Inflows(v%AD, nBlades, nNodes, bladeNodeInflowVel, bladeNodeInflowAcc)

   END SUBROUTINE Interface_SetInflows_C
   
      !----------------------------------------------------------------------------------
   !> 
   SUBROUTINE Interface_SetInflows_Fake_C(simInsAddr, nBlades, nNodes, bladeNodeInflowVel, bladeNodeInflowAcc) BIND(C, NAME='INTERFACE_SETINFLOWS_FAKE')
      !DEC$ ATTRIBUTES DLLEXPORT::Interface_SetInflows_Fake_C
      !GCC$ ATTRIBUTES DLLEXPORT::Interface_SetInflows_Fake_C
      use, intrinsic :: ISO_C_BINDING, ONLY: C_PTR, C_DOUBLE

      type(C_PTR),   value,                     intent(in)    :: simInsAddr
      type(SimInstance), pointer                              :: v
      integer(IntKi),                           intent(in)    :: nBlades, nNodes
      real(ReKi), dimension(3,nNodes,nBlades),  intent(in)    :: bladeNodeInflowVel
      real(ReKi), dimension(3,nNodes,nBlades),  intent(in)    :: bladeNodeInflowAcc

      call C_F_POINTER(simInsAddr, v)
      
      call Set_AD_Inflows(v%AD, nBlades, nNodes, bladeNodeInflowVel, bladeNodeInflowAcc)

   END SUBROUTINE Interface_SetInflows_Fake_C


   
   !---------------------------------------------------------------------------------------------------------------------------------------------------------------
   SUBROUTINE Interface_UpdateStates_C(simInsAddr, force, moment, power, tsr,massMatrix, addedMassMatrix) BIND(C, NAME='INTERFACE_UPDATESTATES')
      !DEC$ ATTRIBUTES DLLEXPORT::Interface_UpdateStates_C
      !GCC$ ATTRIBUTES DLLEXPORT::Interface_UpdateStates_C

      type(C_PTR), value, intent(in) :: simInsAddr
      type(SimInstance), pointer     :: v
      real(C_DOUBLE), dimension(1:3), intent(  out)	:: force     ! 3D force (Fx, Fy, Fz)
      real(C_DOUBLE), dimension(1:3), intent(  out)   :: moment    ! moment (Mx, My, Mz) at rotor hub
      real(C_DOUBLE),                 intent(  out)   :: power     ! power
      real(C_DOUBLE),                 intent(  out)   :: tsr       ! tip-speed ratio
      real(C_DOUBLE), dimension(1:6, 1:6), intent(  out)	:: massMatrix        ! 6x6 mass matrix with mass, inertia, and added mass
      real(C_DOUBLE), dimension(1:6, 1:6), intent(  out)	:: addedMassMatrix   ! 6x6 mass matrix with mass, inertia, and added mass

      integer(IntKi)		        :: errStat  ! Status of error message
      character(ErrMsgLen)	     :: errMsg   ! Error message if ErrStat /= ErrID_None

      call C_F_POINTER(simInsAddr, v)
      
      call Interface_UpdateStates(v%AD, v%DvrData, .TRUE., force, moment, power, tsr, massMatrix, addedMassMatrix, errStat, errMsg)

   END SUBROUTINE Interface_UpdateStates_C
   
   !---------------------------------------------------------------------------------------------------------------------------------------------------------------
   !> Performs a fake advance of states, returning the resulting loads. However, the main states remain unchanged.
   SUBROUTINE Interface_UpdateStates_Fake_C(simInsAddr, force, moment, power, tsr,massMatrix, addedMassMatrix) BIND(C, NAME='INTERFACE_UPDATESTATES_FAKE')
      !DEC$ ATTRIBUTES DLLEXPORT::Interface_UpdateStates_Fake_C
      !GCC$ ATTRIBUTES DLLEXPORT::Interface_UpdateStates_Fake_C

      type(C_PTR), value, intent(in) :: simInsAddr
      type(SimInstance), pointer     :: v
      real(C_DOUBLE), dimension(1:3), intent(  out)	:: force     ! 3D force (Fx, Fy, Fz)
      real(C_DOUBLE), dimension(1:3), intent(  out)   :: moment    ! moment (Mx, My, Mz) at rotor hub
      real(C_DOUBLE),                 intent(  out)   :: power     ! power
      real(C_DOUBLE),                 intent(  out)   :: tsr       ! tip-speed ratio
      real(C_DOUBLE), dimension(1:6, 1:6), intent(  out)	:: massMatrix        ! 6x6 mass matrix with mass, inertia, and added mass
      real(C_DOUBLE), dimension(1:6, 1:6), intent(  out)	:: addedMassMatrix   ! 6x6 mass matrix with mass, inertia, and added mass

      integer(IntKi)		        :: errStat  ! Status of error message
      character(ErrMsgLen)	     :: errMsg   ! Error message if ErrStat /= ErrID_None

      call C_F_POINTER(simInsAddr, v)
      
      call Interface_UpdateStates(v%AD_fake, v%DvrData, .FALSE., force, moment, power, tsr, massMatrix, addedMassMatrix, errStat, errMsg)

   END SUBROUTINE Interface_UpdateStates_Fake_C
   

   !---------------------------------------------------------------------------------------------------------------------------------------------------------------
   SUBROUTINE Interface_GetBladeNodePos_C(simInsAddr, nodePos) BIND(C, NAME='INTERFACE_GETBLADENODEPOS')
      !DEC$ ATTRIBUTES DLLEXPORT::Interface_GetBladeNodePos_C
      !GCC$ ATTRIBUTES DLLEXPORT::Interface_GetBladeNodePos_C
      use, intrinsic :: ISO_C_BINDING, ONLY: C_PTR, C_F_POINTER

      type(C_PTR),     intent(in ), value :: simInsAddr
      type(SimInstance), pointer                    :: v
      real(C_DOUBLE), dimension(1:*), intent(  out) :: nodePos	! variable for outgoing node positions

      integer(C_INT)                              :: errStat   ! Status of error message
      character(ErrMsgLen)		                    :: errMsg    ! Error message if ErrStat /= ErrID_None

      call C_F_POINTER(simInsAddr, v)
      
      ! Access AD internals and calculate positions, passing array reference backwards
      call Interface_GetBladeNodePos(v%AD, v%DvrData, nodePos)

   END SUBROUTINE Interface_GetBladeNodePos_C
   
   !---------------------------------------------------------------------------------------------------------------------------------------------------------------
   SUBROUTINE Interface_GetBladeNodePos_Fake_C(simInsAddr, nodePos) BIND(C, NAME='INTERFACE_GETBLADENODEPOS_FAKE')
      !DEC$ ATTRIBUTES DLLEXPORT::Interface_GetBladeNodePos_Fake_C
      !GCC$ ATTRIBUTES DLLEXPORT::Interface_GetBladeNodePos_Fake_C
      use, intrinsic :: ISO_C_BINDING, ONLY: C_PTR, C_F_POINTER

      type(C_PTR),     intent(in ), value           :: simInsAddr
      type(SimInstance), pointer                    :: v
      real(C_DOUBLE), dimension(1:*), intent(  out) :: nodePos	! variable for outgoing node positions

      integer(C_INT)                                :: errStat   ! Status of error message
      character(ErrMsgLen)		                      :: errMsg    ! Error message if ErrStat /= ErrID_None

      call C_F_POINTER(simInsAddr, v)
      
      ! Access AD internals and calculate positions, passing array reference backwards
      call Interface_GetBladeNodePos(v%AD_fake, v%DvrData, nodePos)

   END SUBROUTINE Interface_GetBladeNodePos_Fake_C


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

         call Interface_End(v%AD, v%DvrData, v%ADInstance_Initialized, ErrStat, ErrMsg)

         call AD_Dvr_DestroyAeroDyn_Data( v%AD_fake, ErrStat, ErrMsg )

         deallocate(v)
      end if

      END SUBROUTINE Interface_End_C

   end module AeroDyn_Interface_C
