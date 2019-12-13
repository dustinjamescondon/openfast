   !
   ! Dustin Condon
   ! 30 July 2019
   !--------------------------------
   ! Built upon the code produced by:
   ! Patrick Connolly
   ! June 26, 2018



   MODULE AeroDyn_Interface

   USE AeroDyn_Interface_Types
   USE AeroDyn_Interface_Subs

   implicit none

   public :: Interface_Init
   public :: Interface_InitInflows
   public :: Interface_SetInflows
   public :: Interface_SetHubMotion
   public :: Interface_GetBladeNodePos
   public :: Interface_UpdateStates
   public :: Interface_End

   Contains


   !-----------------------------------------------------------------------------------------------------------------------------------------------------------
   !>
   SUBROUTINE Interface_Init(AD,DvrData,AD_Initialized,driverFileName,useAddedMass,fluidDensity,kinematicFluidVisc, &
      hubPos,hubOri,hubVel,hubRotVel,bladePitch,nBlades,nNodes,turbineDiameter,errStat,errMsg)
   type(AeroDyn_Data),                 intent(inout) :: AD
   type(Dvr_SimData),                  intent(inout) :: DvrData
   logical,                            intent(inout) :: AD_Initialized
   character(*),                       intent(in   ) :: driverFileName
   logical,                            intent(in   ) :: useAddedMass
   real(C_DOUBLE),                     intent(in   ) :: fluidDensity
   real(C_DOUBLE),                     intent(in   ) :: kinematicFluidVisc
   real(C_DOUBLE), dimension(1:3),     intent(in   ) :: hubPos
   real(C_DOUBLE), dimension(1:3,1:3), intent(in   ) :: hubOri ! The hub's global to local orientation matrix
   real(C_DOUBLE), dimension(1:3),     intent(in   ) :: hubVel
   real(C_DOUBLE), dimension(1:3),     intent(in   ) :: hubRotVel
   real(C_DOUBLE),                     intent(in   ) :: bladePitch
   real(C_DOUBLE),                     intent(out  ) :: turbineDiameter
   integer(IntKi),                     intent(out  ) :: nBlades ! number of blades
   integer(IntKi),                     intent(out  ) :: nNodes  ! number of nodes per blade, as read from the input file
   integer(IntKi),                     intent(out  ) :: errStat                      ! Status of error message
   character(ErrMsgLen),               intent(out  ) :: errMsg                       ! Error message if ErrStat /= ErrID_None


   INTEGER(IntKi)                         :: ErrStat2                     ! error status from any called subroutines
   CHARACTER(ErrMsgLen)                   :: ErrMsg2                      ! error message from any called subroutines
   character(*), parameter                :: RoutineName = 'Interface_Init' ! name of this subroutine for error printing

   character(50) :: atm	! testing variables
   character(50) :: s_y_couple
   character(50) :: p_moment
   real(ReKi)    :: rmax    ! storing the maximum radius of a blade node relative to hub center
   integer       :: j       ! for iterating through number of nodes on one blade

   ! call function to initialize driver data, load in driver files, etc
   ! @mth: note: this function needs to be updated to do all driver things we want
   call Dvr_Init( driverFileName, DvrData, ErrStat, ErrMsg)
   !CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )   ! update error list if applicable
   if ( NeedToAbort(ErrStat) ) then
      return
   end if

   nBlades = DvrData%numBlades
   ! set the these to parameters from the arguments, which will overwrite what was specified in the input file
   AD%p%AirDens = fluidDensity
   AD%p%KinVisc = kinematicFluidVisc

   ! Initialize AeroDyn
   ! Set the Initialization input data for AeroDyn based on the Driver input file data, and initialize AD
   ! (this also initializes inputs to AD for first time step)
   call Init_AeroDyn(DvrData,AD,hubPos,hubOri,hubVel,hubRotVel,bladePitch,DvrData%dt,useAddedMass,errStat2, errMsg2)
   CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )   ! update error list if applicable
   if( NeedToAbort(ErrStat) ) then
      return
   end if

   AD_Initialized = .true.


   nNodes = AD%p%NumBlNds
   
   AD%p%IncludeAddedMass = useAddedMass


   ! Initialize the main output file
   call Dvr_InitializeOutputFile( DvrData%OutFileData, errStat2, errMsg2)
   CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )   ! update error list if applicable
   if ( NeedToAbort(ErrStat) ) then
      return
   endif

   DvrData%iADstep = -1  ! start off the AeroDyn time step counter at -1

   turbineDiameter = CalcTurbineDiameter(AD, nNodes)
   
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

   END SUBROUTINE Interface_Init

   !------------------------------------------------------------------------------------------------------------------------------
   !> Needed to seperate Set_AD_Input from couple subroutine, because this sets the hubstate for the next timestep,
   !! and once we do that, we can get the node positions of that time, which PDS can use to get new inflows, and
   !! finally we can call the coupling subroutine to set the inflows and do the timestep from
   SUBROUTINE Interface_SetHubMotion(AD, DvrData, time, hubPos, hubOri, hubVel, hubRotVel, bladePitch)
      type(AeroDyn_Data),                 intent(inout) :: AD
      type(Dvr_SimData),                  intent(inout) :: DvrData
      real(C_DOUBLE),                     intent(in)    :: time              ! time of these inputs
      real(C_DOUBLE), dimension(1:3),     intent(in)    :: hubPos            !
      real(C_DOUBLE), dimension(1:3,1:3), intent(in)    :: hubOri            ! The global to local hub orientation matrix
      real(C_DOUBLE), dimension(1:3),     intent(in)    :: hubVel            !
      real(C_DOUBLE), dimension(1:3),     intent(in)    :: hubRotVel         !
      real(C_DOUBLE),                     intent(in)    :: bladePitch        !

      ! variables
      real(DbKi)                                   :: dbTime            ! used to cast our C_DOUBLE to DbKi
      integer(IntKi)                               :: errStat              ! Status of error message
      character(ErrMsgLen)                         :: errMsg               ! Error message if ErrStat /= ErrID_None

      dbTime = time
      ! Moves current input information in u(1) to u(2) and then updates u(1) with these hub states.
      call Set_AD_Inputs(dbTime, DvrData, AD, hubPos, hubOri, hubVel, hubRotVel, bladePitch, errStat, errMsg)

   END SUBROUTINE Interface_SetHubMotion

   
   !-------------------------------------------------------------------------------------------
   SUBROUTINE Interface_InitInflows(AD, nBlades, nNodes, bladeNodeInflows)
   type(AeroDyn_Data),                       intent(inout) :: AD
   integer(IntKi),                           intent(in) :: nBlades, nNodes
   real(ReKi), dimension(3,nNodes,nBlades), intent(in) :: bladeNodeInflows

   call Init_AD_Inflows(AD, nBlades, nNodes, bladeNodeInflows)

   END SUBROUTINE Interface_InitInflows

   !----------------------------------------------------------------------------------------------
   SUBROUTINE Interface_SetInflows(AD, nBlades, nNodes, bladeNodeInflows)
   type(AeroDyn_Data),                       intent(inout) :: AD
   integer(IntKi),                              intent(in) :: nBlades, nNodes
   real(ReKi), dimension(3,nNodes,nBlades),  intent(in)    :: bladeNodeInflows

   call Set_AD_Inflows(AD, nBlades, nNodes, bladeNodeInflows)

   END SUBROUTINE Interface_SetInflows


   !-----------------------------------------------------------------------------------------------------------------------------------------
   !> Brings simulation from t_i to t_(i+1), where t_i is the end of the last time-step, and t_(i+1) is the time passed to
   !! Interface_SetHubMotion. Returns the force and moment on the hub at t_(i+1). Note, AD_Step_SetHubstate is
   !! expected to be called before this subroutine.
   SUBROUTINE Interface_UpdateStates(AD, DvrData, isPermanent, force, moment, power, tsr, massMatrix, addedMassMatrix, errStat, errMsg)
   ! include these specific indices used to index into the output array so we can return these values to the calling subroutine
   use AeroDyn_IO, only: RtAeroPwr, RtAeroFxh, RtAeroFyh, RtAeroFzh, RtAeroMxh, RtAeroMyh, RtAeroMzh, RtTSR

   type(AeroDyn_Data),             intent(inout) :: AD
   type(Dvr_SimData),              intent(inout) :: DvrData
   logical, intent(in)                                     :: isPermanent       ! If true, results will be written to output file, and step number will be incremented. If false, neither will happen.

   real(C_DOUBLE), dimension(1:3),          intent(  out)	:: force     ! 3D force (Fx, Fy, Fz)
   real(C_DOUBLE), dimension(1:3),          intent(  out)   :: moment    ! moment (Mx, My, Mz) at rotor hub
   real(C_DOUBLE),                          intent(  out)   :: power     ! power
   real(C_DOUBLE),                          intent(  out)   :: tsr       ! tip-speed ratio
   real(C_DOUBLE), dimension(1:6, 1:6),     intent(  out)	:: massMatrix            ! 6x6 mass matrix with mass, inertia, and added mass
   real(C_DOUBLE), dimension(1:6, 1:6),     intent(  out)	:: addedMassMatrix       ! 6x6 mass matrix with mass, inertia, and added mass

   integer(IntKi),                          intent(  out)   :: errStat              ! Status of error message
   character(ErrMsgLen),                    intent(  out)   :: errMsg               ! Error message if ErrStat /= ErrID_None


   real(DbKi)                 :: dt_total ! variable to hold the dt since last input time
   real(DbKi)                 :: time_ad ! The current time of the states

   ! loop counters
   integer(IntKi) 	::	i
   integer(IntKi) 	::	j
   integer(IntKi) 	::	k

   integer(IntKi)  :: nSubSteps

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

   !.............................
   ! Update states from t_i to t_(i+1)

   ! calculate the time that Aerodyn's internal states are at
   time_ad = DvrData%iADStep * DvrData%dt

   ! calculate the timestep taken by the calling function, which is t_(i+1) - t_i
   dt_total = AD%InputTime(1) - time_ad


   ! if the time-step is larger than AeroDyn's dt, then we break it down into substeps of Aerodyn's dt
   if(dt_total >= (DvrData%dt - epsilon)) then
      ! Adding epsilon so that floating point truncation error won't prevent a substep from happening
      nSubSteps = (dt_total + epsilon) / DvrData%dt

      do i = 1, nSubSteps

         ! Then call UpdateStates with to bring the states from "time" to "time + dt"
         call AD_UpdateStates(time_ad, DvrData%iADstep, AD%u, AD%InputTime, AD%p, AD%x, AD%xd, AD%z, AD%OtherState, AD%m, errStat2, errMsg2 )
         call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )

         ! Update the time to reflect what time the current states represent
         time_ad = time_ad + DvrData%dt

         ! Advance step counter (only if this is a real step)
         if ( isPermanent ) then
            DvrData%iADStep = DvrData%iADStep + 1
         endif

         ! Use the input and states at t_(i+1) to calculate the output at t_(i+1)
         call AD_CalcOutput( time_ad, AD%u(1), AD%p, AD%x, AD%xd, AD%z, AD%OtherState, AD%y, AD%m, errStat2, errMsg2 )
         call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )

         ! If this is a permanent call to simulate then we want to write the results to the output file
         if ( isPermanent ) then
            ! Write the output for t_(i+1) to the output file
            call Dvr_WriteOutputLine(DvrData%OutFileData, time_ad, AD%y%WriteOutput, errStat2, errMsg2)
            call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
         endif
      end do

      ! if there's a remainder of time left, then boohoo
   else 

      ! donothing and then AD_CalcOuput will return the outputs at the last timestep 
   endif 
   
   ! @djc Currently the added mass isn't working because we still need to get the acceleration as the nodes to calculate the added mass
   ! and then we have to figure out how to do the iterative solve to actually get the added mass force right
   !call CalcAddedMassMatrix(addedMassMatrix, AD%m, AD%p)

   !..................................
   ! Calculate outputs at t_(i+1)


   ! @djc: usually would be in this order, with these inputs
   !
   !  calculate the output given by inputs at prev time
   !call AD_CalcOutput( time, AD%u(2), AD%p, AD%x, AD%xd, AD%z, AD%OtherState, AD%y, AD%m, errStat2, errMsg2 )
   !call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
   !
   !  write the ouputs at prev time to the output file
   !call Dvr_WriteOutputLine(DvrData%OutFileData, time, AD%y%WriteOutput, errStat2, errMsg2)
   !call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
   !
   !  update states from prev time to current time
   !call AD_UpdateStates( time , DvrData%iADstep, AD%u, AD%InputTime, AD%p, AD%x, AD%xd, AD%z, AD%OtherState, AD%m, errStat2, errMsg2 )
   !call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )

   ! ------------------------- Mass and Inertia matrix ------------------------------
   ! cross-coupling terms exist since coupling point is at the hub

   massMatrix(1,1)  = DvrData%MassAndInertia(1)	! xx
   massMatrix(2,2)  = DvrData%MassAndInertia(1)	! yy

   massMatrix(6,2) = DvrData%MassAndInertia(2)	! sway-yaw coupling
   massMatrix(2,6) = DvrData%MassAndInertia(2)	! yaw-sway coupling

   massMatrix(3,3) = DvrData%MassAndInertia(1)	! zz

   massMatrix(5,3) = DvrData%MassAndInertia(3)	! heave-pitch coupling
   massMatrix(3,5) = DvrData%MassAndInertia(3)	! pitch-heave coupling

   massMatrix(4,4) = DvrData%MassAndInertia(4)	! Roll moment of inertia, about rotor axis
   massMatrix(5,5) = DvrData%MassAndInertia(5)	! Pitch moment of inerita, about axis perpendicular to rotor shaft
   massMatrix(6,6) = DvrData%MassAndInertia(6)	! Yaw moment of inerita, about axis perpendicular to rotor shaft

   ! return the hub force and moment in the hub coordinate system back to the calling function
   force(1) = AD%m%AllOuts( RtAeroFxh )
   force(2) = AD%m%AllOuts( RtAeroFyh )
   force(3) = AD%m%AllOuts( RtAeroFzh )

   moment(1) = AD%m%AllOuts( RtAeroMxh )
   moment(2) = AD%m%AllOuts( RtAeroMyh )
   moment(3) = AD%m%AllOuts( RtAeroMzh )

   power = AD%m%AllOuts( RtAeroPwr )
   tsr = AD%m%AllOuts (RtTSR)

   END SUBROUTINE Interface_UpdateStates

   ! gets the node coordinates of all the blade nodes from AeroDyn
   SUBROUTINE Interface_GetBladeNodePos(AD, DvrData, bladeNodePos)

   type(AeroDyn_Data),        intent(inout) :: AD
   type(Dvr_SimData),         intent(in   ) :: DvrData
   real(ReKi), dimension(1:*), intent(out)  :: bladeNodePos	!output variable for blade node positions

   integer(intKi)		:: j             ! loop counter for nodes
   integer(intKi)		:: k             ! loop counter for blades

   real(ReKi)      :: nodePos(3)
   integer(8)		:: curIndex

   ! Blade motions:
   do k=1,DvrData%numBlades

      do j=1,AD%u(1)%BladeMotion(k)%nnodes

         curIndex = ((k-1) * AD%u(1)%BladeMotion(k)%nnodes * 3)  + ((j-1) * 3)

         nodePos = AD%u(1)%BladeMotion(k)%Position(:,j) + AD%u(1)%BladeMotion(k)%TranslationDisp(:,j)

         bladeNodePos(curIndex+1) = nodePos(1)  ! and these are the even more absolute locations?
         bladeNodePos(curIndex+2) = nodePos(2)
         bladeNodePos(curIndex+3) = nodePos(3)

      end do !j=nnodes

   end do !k=numBlades

   END SUBROUTINE Interface_GetBladeNodePos


   !SUBROUTINE Force_Moment_Torque_Mass(time, hubForceAndMoment, genTorque, massMatrix)
   !
   !   real(DbKi),      	intent(inout)	:: time	! current simulation time
   !   real(DbKi),	Dimension(:),	intent(  out)	:: hubForceAndMoment	! 3D force (Fx, Fy, Fz) and moment (Mx, My, Mz) at rotor hub
   !   real(DbKi),		intent(  out)	:: genTorque	! generator torque 1D
   !   real(DbKi),	Dimension(:),	intent(  out)	:: massMatrix	! 6x6 mass matrix with mass, inertia, and added mass
   !
   !   !integer(IntKi)			:: errStat	! Status of error message
   !   !character(ErrMsgLen)		:: errMsg	! Error message if ErrStat /= ErrID_None
   !
   !
   !   real(DbKi) 		:: total_AM_vector(6)	! holds the sum of the added mass force vectors
   !   real(DbKi)		:: AM_Dir_UnitVec(6)	! unit vector in the direction of the added mass force (flapwise)
   !   real(DbKi)		:: square_mag		! magnitude of FAddedMass_vec(:,j,k)
   !
   !   ! loop counters
   !   integer(IntKi) 	::	i
   !   integer(IntKi) 	::	j
   !   integer(IntKi) 	::	k
   !
   !   ! Loops and logic to get added mass matrix  - todo
   !   ! AM
   !   do k=1, DvrData%NumBlades
   !
   !	do j=1, AD%p%numBlNds
   !
   !		!summing added mass in vector form node by node for each blade
   !
   !		! reset force magnitude
   !		square_mag = 0.0
   !
   !		do i=1, 6
   !			! getting magnitude of added_mass force vector so that the unit vector can be used
   !			square_mag = square_mag + ((AD%m%FAddedMass_vec(i,j,k))**2)
   !
   !		! Define force direction unit vector to sort into flapwise direction
   !		AM_Dir_UnitVec = AD%m%FAddedMass_vec(:,j,k) / (square_mag**(1/2))
   !
   !		! Add the added mass vector from this node to the total for this blade
   !		total_AM_vector = total_AM_vector + (AM_Dir_UnitVec * AD%m%addedMass(j,k))
   !
   !		end do ! dofs
   !
   !	end do ! blade nodes
   !
   !	! now add the total from this blade to the mass matrix
   !	massMatrix(1) = massMatrix(1) + total_AM_vector(1)	! x,x -> 1,1 -> (1-1)*6 + (1*1) -> 1
   !	massMatrix(8) = massMatrix(8) + total_AM_vector(2)	! y,y -> 2,2 -> (2-1)*6 + (2*1) -> 8
   !	massMatrix(15) = massMatrix(15) + total_AM_vector(3)	! z,z -> 3,3 -> (3-1)*6 + (3*1) -> 15
   !	massMatrix(22) = massMatrix(22) + total_AM_vector(4)	! r,r -> 1,1 -> (4-1)*6 + (4*1) -> 22
   !	massMatrix(29) = massMatrix(29) + total_AM_vector(5)	! p,p -> 1,1 -> (5-1)*6 + (5*1) -> 29
   !	massMatrix(36) = massMatrix(36) + total_AM_vector(6)	! w,w -> 1,1 -> (6-1)*6 + (6*1) -> 36
   !
   !	!reset the local variables
   !	total_AM_vector(1:6) = 0.0
   !	AM_Dir_UnitVec(1:6) = 0.0
   !
   !   end do ! num blades
   !
   !
   !   ! @PC: assuming DvrData is accessible at this level
   !
   !   ! Mass and Inertia
   !   ! cross-coupling terms exist since coupling point is at the hub
   !
   !   massMatrix(1) = massMatrix(1) + DvrData%MassAndInertia(1)	! xx
   !   massMatrix(8) = massMatrix(8) + DvrData%MassAndInertia(1)	! yy
   !
   !   massMatrix(12) = massMatrix(12) + DvrData%MassAndInertia(2)	! sway-yaw coupling
   !   massMatrix(32) = massMatrix(32) + DvrData%MassAndInertia(2)	! yaw-sway coupling
   !
   !   massMatrix(15) = massMatrix(15) + DvrData%MassAndInertia(1)	! zz
   !
   !   massMatrix(17) = massMatrix(17) + DvrData%MassAndInertia(3)	! heave-pitch coupling
   !   massMatrix(27) = massMatrix(27) + DvrData%MassAndInertia(3)	! pitch-heave coupling
   !
   !   massMatrix(22) = massMatrix(22) + DvrData%MassAndInertia(4)	! Roll moment of inertia, about rotor axis
   !   massMatrix(29) = massMatrix(29) + DvrData%MassAndInertia(5)	! Pitch moment of inerita, about axis perpendicular to rotor shaft
   !   massMatrix(36) = massMatrix(36) + DvrData%MassAndInertia(6)	! Yaw moment of inerita, about axis perpendicular to rotor shaft
   !
   !   ! logic to get generator torque - todo
   !   ! - store a single value and return it?
   !   genTorque = 0.0_DbKi
   !
   !   ! logic to get Force and Moment at the hub
   !   hubForceAndMoment(1) = AD%m%HubLoad%Force(1,1)
   !   hubForceAndMoment(2) = AD%m%HubLoad%Force(2,1)
   !   hubForceAndMoment(3) = AD%m%HubLoad%Force(3,1)
   !
   !   hubForceAndMoment(1) = AD%m%HubLoad%Moment(1,1)
   !   hubForceAndMoment(2) = AD%m%HubLoad%Moment(2,1)
   !   hubForceAndMoment(3) = AD%m%HubLoad%Moment(3,1)
   !
   !
   !
   !END SUBROUTINE ! Force_Moment_Torque_Mass


   ! Returns true if we need to abort, and false if otherwise
   FUNCTION NeedToAbort(ErrStat) result(abort)
   integer(IntKi),         intent(in) :: ErrStat              ! Status of error message

   logical                            :: abort

   abort = .false.

   if (ErrStat >= AbortErrLev) then
      abort = .true.
   end if

   END FUNCTION NeedToAbort


   SUBROUTINE Interface_End(AD, DvrData, AD_Initialized, errStat, errMsg)

   ! Local variables
   type(AeroDyn_Data),  intent(inout)  :: AD
   type(Dvr_SimData),   intent(inout)  :: DvrData
   logical,	               intent(in)  :: AD_Initialized
   integer(IntKi),         intent(out) :: errStat	! Status of error message
   character(ErrMsgLen),   intent(out) :: errMsg	! Error message if ErrStat /= ErrID_None

   character(ErrMsgLen)                :: errMsg2                 ! temporary Error message if ErrStat /= ErrID_None
   integer(IntKi)		                  :: errStat2                ! temporary Error status of the operation

   character(*), parameter	:: RoutineName = 'Interface_End'
   ! Close the output file
   if (DvrData%OutFileData%unOutFile > 0) close(DvrData%OutFileData%unOutFile)

   if ( AD_Initialized ) then	! (?)
      call AD_End( AD%u(1), AD%p, AD%x, AD%xd, AD%z, AD%OtherState, AD%y, AD%m, errStat2, errMsg2 )
      call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
   end if

   call AD_Dvr_DestroyDvr_SimData( DvrData, ErrStat2, ErrMsg2 )
   call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )

   call AD_Dvr_DestroyAeroDyn_Data( AD, ErrStat2, ErrMsg2 )
   call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )

   if (ErrStat >= AbortErrLev) then
      !CALL ProgAbort( 'AeroDyn Driver encountered simulation error level: '&
      !//TRIM(GetErrStr(ErrStat)), TrapErrors=.FALSE., TimeWait=3._ReKi )  ! wait 3 seconds (in case they double-clicked and got an error)
   else
      !call NormStop() ! don't call this because it will end the entire process that calls it
   end if

   END SUBROUTINE Interface_End

   END Module AeroDyn_Interface
