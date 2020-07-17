   !**********************************************************************************************************************************
   ! LICENSING
   ! Copyright (C) 2015-2016  National Renewable Energy Laboratory
   ! Copyright (C) 2016-2018  Envision Energy USA, LTD
   !
   !    This file is part of AeroDyn.
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
   module AeroDyn_Interface_Subs

   use AeroDyn_Interface_Types
   use AeroDyn
   use AeroDyn_Types
   use VersionInfo

   implicit none
   
   ! state array indexes
   INTEGER(IntKi), PARAMETER :: STATE_CURR              = 1          !< index for "current" (t_global) states
   INTEGER(IntKi), PARAMETER :: STATE_PRED              = 2          !< index for "predicted" (t_global_next) states

   TYPE(ProgDesc), PARAMETER   :: version   = ProgDesc( 'AeroDyn_Interface_Subs', '', '' )  ! The version number of this program.
   
   contains
   !---------------------------------------------------------------
   !> Returns true if we need to abort, and false if otherwise
   FUNCTION NeedToAbort(ErrStat) result(abort)
   integer(IntKi),         intent(in) :: ErrStat              ! Status of error message

   logical                          :: abort

   abort = .false.

   if (ErrStat >= AbortErrLev) then
      abort = .true.
   end if

   END FUNCTION NeedToAbort

   !----------------------------------------------------------------------------------------------------------------------------------
   !< TODO comment this
   subroutine Interface_Init_Parameters(inputFile,outputfile,timestep,numBlades,hubRad,precone,DvrData,errStat,errMsg )

   CHARACTER(*),              intent(in   ) :: InputFile     ! file name for AeroDyn input file
   CHARACTER(*),              intent(in   ) :: OutputFile    ! file name for the output of the simulation
   real(ReKi),                intent(in   ) :: timestep
   integer(IntKi),            intent(in   ) :: numBlades 
   real(ReKi),                intent(in   ) :: hubRad
   real(ReKi),                intent(in   ) :: precone       ! radians
   type(Dvr_SimData),         intent(  out) :: DvrData       ! driver data
   integer(IntKi),            intent(  out) :: errStat       ! Status of error message
   character(*),              intent(  out) :: errMsg        ! Error message if ErrStat /= ErrID_None

   ! local variables
   integer(IntKi)                              :: errStat2      ! local status of error message
   character(ErrMsgLen)                        :: errMsg2       ! local error message if ErrStat /= ErrID_None
   character(*), parameter                     :: RoutineName = 'Interface_Init_Parameters'

   CHARACTER(200)                              :: git_commit    ! String containing the current git commit hash

   TYPE(ProgDesc), PARAMETER                   :: version   = ProgDesc( 'AeroDyn Interface', '', '' )  ! The version number of this program.

   ErrStat = ErrID_None
   ErrMsg  = ""


   DvrData%OutFileData%unOutFile   = -1

   CALL NWTC_Init()
   ! Display the copyright notice
   CALL DispCopyrightLicense( version ) ! djc: commented this out before because it kept hanging on a Write statement
   ! Obtain OpenFAST git commit hash
   git_commit = QueryGitVersion()
   ! Tell our users what they're running
   CALL WrScr( ' Running '//GetNVD( version )//' a part of OpenFAST - '//TRIM(git_Commit)//NewLine//' linked with '//TRIM( GetNVD( NWTC_Ver ))//NewLine )

   ! set the interface parameters
   DvrData%dt     = timestep
   DvrData%hubRad = hubRad
   DvrData%numBlades = numBlades 
   DvrData%AD_InputFile = inputFile ! Note: may need to do something special with this
   DvrData%OutFileData%Root = outputFile ! hard-code for now 
   DvrData%precone = precone
   DvrData%OutFileData%delim = TAB ! hard-code for now
   DvrData%OutFileData%OutFmt = "ES10.3E2" 
   
   ! validate the inputs
   call WrScr('validating input')
   call ValidateInputs(DvrData, errStat2, errMsg2)
   call SetErrStat(errStat2, errMsg2, ErrStat, ErrMsg, RoutineName)

   end subroutine Interface_Init_Parameters


   !----------------------------------------------------------------------------------------------------------------------------------
   subroutine Init_AeroDyn(DvrData, AD, hubPos, hubOri, hubVel, hubAcc, hubRotVel, hubRotAcc, bladePitch, dt, useAddedMass, coeffAddedMass, errStat, errMsg)

   type(Dvr_SimData),              intent(inout) :: DvrData       ! Input data for initialization
   type(AeroDyn_Data),             intent(inout) :: AD            ! AeroDyn data
   real(R8Ki), dimension(1:3),     intent(in   ) :: hubPos
   real(R8Ki), dimension(1:3,1:3), intent(in   ) :: hubOri        ! The hub's initial orientation matrix (local to global)
   real(R8Ki), dimension(1:3),     intent(in   ) :: hubVel
   real(R8Ki), dimension(1:3),     intent(in   ) :: hubAcc
   real(R8Ki), dimension(1:3),     intent(in   ) :: hubRotVel
   real(R8Ki), dimension(1:3),     intent(in   ) :: hubRotAcc
   real(R8Ki),                     intent(in   ) :: bladePitch
   real(DbKi),                     intent(inout) :: dt
   logical(kind=C_BOOL),           intent(in   ) :: useAddedMass
   real(DbKi),                     intent(in   ) :: coeffAddedMass
   integer(IntKi),                 intent(  out) :: errStat       ! Status of error message
   character(*),                   intent(  out) :: errMsg        ! Error message if ErrStat /= ErrID_None

   ! locals
   real(DbKi)                                  :: time
   real(reKi)                                  :: theta(3)
   integer(IntKi)                              :: j, k
   integer(IntKi)                              :: errStat2      ! local status of error message
   character(ErrMsgLen)                        :: errMsg2       ! local error message if ErrStat /= ErrID_None
   character(*), parameter                     :: RoutineName = 'Init_AeroDyn'

   real(R8Ki), dimension(1:3,1:3)              :: extrapOri  ! holds extrapolated orientation matrix
   ! local data

   type(AD_InitInputType)                      :: InitInData     ! Input data for initialization
   type(AD_InitOutputType)                     :: InitOutData    ! Output data from initialization
   integer, save :: turbNum = 0                                 ! Keeps track of how many turbines have been initialized
   character(len=3) :: turbNumStr                                ! String conversion of turbNum

   errStat = ErrID_None
   errMsg  = ''
      
   ! @djc: added this code to make sure each turbine instance outputs to a different file, so 
   ! multiple simulation instances can exist without crashing because of trying to access the same
   ! file. Eventually, though, we will probably silence AeroDyn so it doesn't actually ouput to a file,  
   ! so ProteusDS can be in control of what files are created and output to.
   turbNum = turbNum + 1
   write(turbNumStr, "(I3)") turbNum

   InitInData%InputFile      = DvrData%AD_InputFile
   InitInData%NumBlades      = DvrData%numBlades
   DvrData%OutFileData%Root  = trim(DvrData%outFileData%Root) // trim(turbNumStr)

   InitInData%RootName       = DvrData%outFileData%Root
   InitInData%Gravity        = 9.80665_ReKi


   ! set initialization data:
   call AllocAry( InitInData%BladeRootPosition, 3, InitInData%NumBlades, 'BladeRootPosition', errStat2, ErrMsg2 )
   call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
   call AllocAry( InitInData%BladeRootOrientation, 3, 3, InitInData%NumBlades, 'BladeRootOrientation', errStat2, ErrMsg2 )
   call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )

   if (ErrStat >= AbortErrLev) then
      call Cleanup()
      return
   end if

   ! @djc: we will always have the hub position be 0 vector in Aerodyn's view.
   InitInData%HubPosition = 0.0_ReKi
 
   ! @djc: This sets ref orientation to the identity matrix: we want the ref orientation to be the identity matrix
   theta(1) = 0.0
   theta(2) = 0.0
   theta(3) = 0.0
   InitInData%HubOrientation = EulerConstruct( theta )

   do k=1,InitInData%numBlades

      theta(1) = (k-1)*TwoPi/real(InitInData%numBlades,ReKi)
      theta(2) = DvrData%precone
      theta(3) = 0.0_ReKi 
      InitInData%BladeRootOrientation(:,:,k) = matmul( EulerConstruct( theta ), InitInData%HubOrientation )

      InitInData%BladeRootPosition(:,k)   = InitInData%HubPosition + DvrData%hubRad * InitInData%BladeRootOrientation(3,:,k)

   end do

   ! Initialize AeroDyn
   call AD_Init(InitInData, AD%u(1), AD%p, AD%x(STATE_CURR), AD%xd(STATE_CURR), AD%z(STATE_CURR), AD%OtherState(STATE_CURR), AD%y, AD%m, dt, InitOutData, ErrStat2, ErrMsg2 )
   call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
   
   CALL AD_CopyContState   (AD%x( STATE_CURR), AD%x( STATE_PRED), MESH_NEWCOPY, Errstat2, ErrMsg2)
      CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AD_CopyDiscState   (AD%xd(STATE_CURR), AD%xd(STATE_PRED), MESH_NEWCOPY, Errstat2, ErrMsg2)  
      CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AD_CopyConstrState (AD%z( STATE_CURR), AD%z( STATE_PRED), MESH_NEWCOPY, Errstat2, ErrMsg2)
      CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AD_CopyOtherState( AD%OtherState(STATE_CURR), AD%OtherState(STATE_PRED), MESH_NEWCOPY, Errstat2, ErrMsg2)
      CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   do j = 2, numInp
      call AD_CopyInput (AD%u(1),  AD%u(j),  MESH_NEWCOPY, errStat2, errMsg2)
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   end do

   if (ErrStat >= AbortErrLev) then
      call Cleanup()
      return
   end if

   ! we know exact values, so we're going to initialize inputs this way (instead of using the input guesses from AD_Init)
   AD%InputTime = -999
   AD%stepNum = -2 ! TODO confirm this is right: We want the step number to be on -1 once initialization is done?
   !  an orientation matrix for roll according to rotvel and time
   ! and multiply this on LHS by orientation at t=0
   DO j = -numInp + 1, -1
      extrapOri(:,:) = ExtrapOrientationFromRotVel(hubOri, hubRotVel, dt * j)
      time = dt * j
      call Advance_AD_InputWindow(AD, errStat2, errMsg2)
         call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
         
      AD%InputTime(1) = time

      call Set_AD_Inputs_Hub(time,DvrData,AD%u(1),hubPos,extrapOri,hubVel,hubAcc,hubRotVel,hubRotAcc,bladePitch)
         call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

   END DO

   ! We move the input window so that nstep 0 is ready for Set_Inputs
   call Advance_AD_InputWindow(AD, errstat2, errmsg2)
   
   AD%p%IncludeAddedMass = useAddedMass
   AD%p%CaBlade = coeffAddedMass
   
   ! move AD initOut data to AD Driver
   call move_alloc( InitOutData%WriteOutputHdr, DvrData%OutFileData%WriteOutputHdr )
   call move_alloc( InitOutData%WriteOutputUnt, DvrData%OutFileData%WriteOutputUnt )

   DvrData%OutFileData%AD_ver = InitOutData%ver

   contains

   function ExtrapOrientationFromRotVel(orientation, rotVel, time) result(orientation_prime)
   use NWTC_Num, only : EulerExtract
      real(R8Ki), dimension(1:3,1:3), intent(in) :: orientation ! the global to local orientation matrix
      real(R8Ki), dimension(1:3),     intent(in) :: rotVel ! rotational velocity in axis-angle form
      real(DbKi),                     intent(in) :: time   ! time to extrapolate to

      real(ReKi), dimension(1:3)             :: theta
      real(ReKi), dimension(1:3)             :: theta_prime
      real(R8Ki), dimension(1:3,1:3)         :: orientation_prime
      real(ReKi), dimension(1:3,1:3)         :: active_roll_matrix
      real(R8Ki), dimension(1:3)             :: xBase, yBase, zBase, xBase_prime, yBase_prime, zBase_prime
      real(DbKi), dimension(1:3)             :: rotVel_x_time

      ! Just assume that its going to rotate about the hub's x axis
      theta(1) = sqrt(dot_product(rotVel, rotVel)) * time
      theta(2) = 0_ReKi
      theta(3) = 0_ReKi
      active_roll_matrix(:,:) = transpose(EulerConstruct(theta))
      orientation_prime = orientation * active_roll_matrix
      
   end function ExtrapOrientationFromRotVel

   subroutine cleanup()
   call AD_DestroyInitInput( InitInData, ErrStat2, ErrMsg2 )
   call AD_DestroyInitOutput( InitOutData, ErrStat2, ErrMsg2 )
   end subroutine cleanup

   end subroutine Init_AeroDyn

 
   !-----------------------------
   !< Updates the blade node accelerations, positions, and velocities based on the hub's current rotational vel, acc, and linear acc
   subroutine Calc_AD_NodeKinematics(hub_linear_acc, hub_angular_acc, u_AD, DvrData)
   !.............................
   real(C_DOUBLE),               dimension(1:3), intent(in   ) :: hub_linear_acc
   real(C_DOUBLE),               dimension(1:3), intent(in   ) :: hub_angular_acc
   type(AD_InputType),                      intent(inout) :: u_AD
   type(Dvr_SimData),                            intent(inout) :: DvrData       ! Driver data
   
   real(ReKi)                                  :: rotVel_cross_offset(3) ! used in calculating blade node accelerations and velocities
   real(ReKi)                                  :: rotAcc_cross_offset(3) ! used in calculating blade node accelerations
   real(ReKi)                                  :: orientation(3,3)
   real(ReKi)                                  :: rotateMat(3,3)
   real(ReKi)                                  :: position(3)
   integer(IntKi) :: k ! blade index
   integer(IntKi) :: j ! blade node index

   ! First set the input hub acceleration fields directly
   u_AD%HubMotion%TranslationAcc(:,1) = hub_linear_acc
   u_AD%HubMotion%RotationAcc(:,1)    = hub_angular_acc
   
   ! Then calculate the blade node accelerations based on those previous accelerations
   do k=1,DvrData%numBlades
      rotateMat = transpose( u_AD%BladeRootMotion(k)%Orientation(  :,:,1) )
      rotateMat = matmul( rotateMat, u_AD%BladeRootMotion(k)%RefOrientation(  :,:,1) )
      orientation = transpose(rotateMat)

      rotateMat(1,1) = rotateMat(1,1) - 1.0_ReKi
      rotateMat(2,2) = rotateMat(2,2) - 1.0_ReKi
      rotateMat(3,3) = rotateMat(3,3) - 1.0_ReKi

      do j=1,u_AD%BladeMotion(k)%nnodes
         position = u_AD%BladeMotion(k)%Position(:,j)
         u_AD%BladeMotion(k)%TranslationDisp(:,j) = matmul( rotateMat, position ) + u_AD%HubMotion%Position(:,1)

         u_AD%BladeMotion(k)%Orientation(  :,:,j) = matmul( u_AD%BladeMotion(k)%RefOrientation(:,:,j), orientation )

         position =  u_AD%BladeMotion(k)%Position(:,j) + u_AD%BladeMotion(k)%TranslationDisp(:,j) &
            - u_AD%HubMotion%Position(:,1) - u_AD%HubMotion%TranslationDisp(:,1)
            
         rotVel_cross_offset = cross_product( u_AD%HubMotion%RotationVel(:,1), position )
         rotAcc_cross_offset = cross_product( u_AD%HubMotion%RotationAcc(:,1), position )
         u_AD%BladeMotion(k)%TranslationVel( :,j) = rotVel_cross_offset + u_AD%HubMotion%TranslationVel( :,1)
         u_AD%BladeMotion(k)%TranslationAcc( :,j) = rotAcc_cross_offset + rotVel_cross_offset + u_AD%HubMotion%TranslationAcc( :,1)

      end do !j=nnodes

   end do !k=numBlades   
   end subroutine Calc_AD_NodeKinematics
   
   !----------------------------------------------------------------------------------------------------------------------------------
   !> TODO use Calc_AD_NodeKinematics here instead of repeating code
   subroutine Set_AD_Inputs_Hub(time,DvrData,u,hubPos,hubOri,hubVel,hubAcc,hubRotVel,hubRotAcc,bladePitch)

   real(DbKi),                         intent(in   ) :: time          ! time of the inputs (changed from timestep number because now using variable timestep)
   type(Dvr_SimData),                  intent(inout) :: DvrData       ! Driver data
   type(AD_InputType),                 intent(inout) :: u             ! AeroDyn input
   real(C_DOUBLE), dimension(1:3),     intent(in   ) :: hubPos        ! x,y,z in meters
   real(C_DOUBLE), dimension(1:3,1:3), intent(in   ) :: hubOri        ! the hub's orientation matrix (local to global)
   real(C_DOUBLE), dimension(1:3),     intent(in   ) :: hubVel        ! x,y,z in meters/sec
   real(C_DOUBLE), dimension(1:3),     intent(in   ) :: hubAcc        ! x,y,z in meters/sec^2
   real(C_DOUBLE), dimension(1:3),     intent(in   ) :: hubRotVel     ! axis-angle in global coordinate system, (magnitute in radians/sec)
   real(C_DOUBLE), dimension(1:3),     intent(in   ) :: hubRotAcc     ! axis-angle in global coordinate system, (magnitute in radians/sec^2)
   real(C_DOUBLE),                     intent(in   ) :: bladePitch

   ! local variables
   integer(intKi)		                          :: i 	          ! loop counter for dimensions
   integer(intKi)                              :: j             ! loop counter for nodes
   integer(intKi)                              :: k             ! loop counter for blades

   real(ReKi)                                  :: theta(3)
   real(ReKi)                                  :: position(3)
   real(ReKi)                                  :: rotVel_cross_offset(3) ! used in calculating blade node accelerations and velocities
   real(ReKi)                                  :: rotAcc_cross_offset(3) ! used in calculating blade node accelerations
   real(ReKi)                                  :: orientation(3,3)
   real(ReKi)                                  :: rotateMat(3,3)
   integer(8)                                  :: curIndex

   !................
   ! calculate new values
   !................
      
   ! Tower motions:
   do j=1,u%TowerMotion%nnodes
      u%TowerMotion%Orientation(  :,:,j) = u%TowerMotion%RefOrientation(:,:,j) ! identity
      u%TowerMotion%TranslationDisp(:,j) = 0.0_ReKi
      u%TowerMotion%TranslationVel( :,j) = 0.0_ReKi
   end do !j=nnodes

   ! Hub motions:
   u%HubMotion%Position(:,1) = hubPos
   
   ! This would usually express where the hub is rotated to around the tower. Since we want hub to rotate around its center -
   ! i.e. not change position upon rotation - we just set this to zero.
   u%HubMotion%TranslationDisp(:,1) = 0.0

   ! Rather than pass the Euler angles we pass the orientation matrix itself
   ! Note that HubMotion's orientation needs to be global to local, and our parameter
   ! is local to global, so we do the transpose here
   u%HubMotion%Orientation(  :,:,1) = transpose(hubOri(:,:))
 
   ! set the rotational velocity directly from the argument (expected to be global axis-angle rotational vel)
   u%HubMotion%RotationVel(    :,1) = hubRotVel
   u%HubMotion%RotationAcc(    :,1) = hubRotAcc

   u%HubMotion%TranslationVel(:,1) = hubVel
   u%HubMotion%TranslationAcc(:,1) = hubAcc

   ! Blade root motions:

   do k=1,DvrData%numBlades
      theta(1) = (k-1)*TwoPi/real(DvrData%numBlades,ReKi)
      theta(2) =  DvrData%precone
      theta(3) = -bladePitch
      orientation = EulerConstruct(theta)

      u%BladeRootMotion(k)%Orientation(  :,:,1) = matmul( orientation, u%HubMotion%Orientation(  :,:,1) )

   end do !k=numBlades

   ! Blade motions:
   do k=1,DvrData%numBlades
      rotateMat = transpose( u%BladeRootMotion(k)%Orientation(  :,:,1) )
      rotateMat = matmul( rotateMat, u%BladeRootMotion(k)%RefOrientation(  :,:,1) )
      orientation = transpose(rotateMat)

      rotateMat(1,1) = rotateMat(1,1) - 1.0_ReKi
      rotateMat(2,2) = rotateMat(2,2) - 1.0_ReKi
      rotateMat(3,3) = rotateMat(3,3) - 1.0_ReKi

      do j=1,u%BladeMotion(k)%nnodes
         position = u%BladeMotion(k)%Position(:,j)
         u%BladeMotion(k)%TranslationDisp(:,j) = matmul( rotateMat, position ) + u%HubMotion%Position(:,1)

         u%BladeMotion(k)%Orientation(  :,:,j) = matmul( u%BladeMotion(k)%RefOrientation(:,:,j), orientation )

         position =  u%BladeMotion(k)%Position(:,j) + u%BladeMotion(k)%TranslationDisp(:,j) &
            - u%HubMotion%Position(:,1) - u%HubMotion%TranslationDisp(:,1) 
         
         rotVel_cross_offset = cross_product( u%HubMotion%RotationVel(:,1), position )
         rotAcc_cross_offset = cross_product( u%HubMotion%RotationAcc(:,1), position )
         u%BladeMotion(k)%TranslationVel( :,j) = rotVel_cross_offset + hubVel ! add hub vel because hub can have a linear velocity
         u%BladeMotion(k)%TranslationAcc( :,j) = rotAcc_cross_offset + rotVel_cross_offset + hubAcc

      end do !j=nnodes

   end do !k=numBlades

   end subroutine Set_AD_Inputs_Hub
   
   !----------------------------------------------------------------------------------------------------------------------------------
   !> this routine cycles values in the input array AD%InputTime and AD%u.
   subroutine Advance_AD_InputWindow(AD,errStat,errMsg)

   type(AeroDyn_Data),                 intent(inout) :: AD            ! AeroDyn data
   integer(IntKi),                     intent(  out) :: errStat       ! Status of error message
   character(*),                       intent(  out) :: errMsg        ! Error message if ErrStat /= ErrID_None

   ! local variables
   integer(IntKi)                              :: errStat2      ! local status of error message
   character(ErrMsgLen)                        :: errMsg2       ! local error message if ErrStat /= ErrID_None
   character(*), parameter                     :: RoutineName = 'Advance_AD_InputWindow'

   integer(intKi)                              :: j             ! loop counter inputs


   errStat = ErrID_None
   errMsg  = ""
   
   !................
   ! shift previous calculations:
   !................
   do j = numInp-1,1,-1
      call AD_CopyInput (AD%u(j),  AD%u(j+1),  MESH_UPDATECOPY, ErrStat2, ErrMsg2)
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

      AD%InputTime(j+1) = AD%InputTime(j)
   end do
   
   AD%stepNum = AD%stepNum + 1

   end subroutine Advance_AD_InputWindow
   
   !---------------------------------------------------------------------------------------------------------------------------------
   !> Refer to normal Aerodyn_Driver's Set_AD_Inputs_and_Advance_Window to see how I'm trying to keep the end result the same after taking
   !! set inflows out of that subroutine
   SUBROUTINE Init_AD_Inputs_Inflow(AD, nBlades, nNodes, bladeNodeInflowVel, bladeNodeInflowAcc)
      type(AeroDyn_Data),                       intent(inout) :: AD
      integer(IntKi),                              intent(in) :: nBlades, nNodes
      real(R8Ki), dimension(3,nNodes, nBlades), intent(in) :: bladeNodeInflowVel
      real(R8Ki), dimension(3,nNodes, nBlades), intent(in) :: bladeNodeInflowAcc


      integer :: i, j, k
      do i = 1, numInp ! This is declared globally in AeroDyn
         do k = 1, nBlades
            do j = 1, nNodes  

               AD%u(i)%InflowOnBlade(1:3,j,k) = bladeNodeInflowVel(1:3,j,k)
               AD%u(i)%InflowAccOnBlade(1:3,j,k) = bladeNodeInflowAcc(1:3,j,k)

            enddo !j=nnodes
         enddo !k=nblades
      enddo !i=numInp

   END SUBROUTINE Init_AD_Inputs_Inflow

   !------------------------------------------------------------------------------------------------------------------------------
   !> TODO comment this
   SUBROUTINE Set_AD_Inputs_Inflow(AD, nBlades, nNodes, bladeNodeInflowVel, bladeNodeInflowAcc)
      type(AeroDyn_Data),                       intent(inout)  :: AD
      integer(IntKi),                           intent(in   )  :: nBlades, nNodes
      real(ReKi), dimension(3,nNodes, nBlades), intent(in   )  :: bladeNodeInflowVel
      real(ReKi), dimension(3,nNodes, nBlades), intent(in   )  :: bladeNodeInflowAcc

      integer(IntKi) :: j, k

      do k = 1, nBlades
         do j = 1, nNodes
            AD%u(1)%InflowOnBlade(1:3,j,k) = bladeNodeInflowVel(1:3,j,k)
            AD%u(1)%InflowAccOnBlade(1:3,j,k) = bladeNodeInflowAcc(1:3,j,k)
         enddo !j=nnodes
      enddo !k=nblades

   END SUBROUTINE Set_AD_Inputs_Inflow

   !----------------------------------------------------------------------------------------------------------------------------
   !> TODO comment this
   subroutine ValidateInputs(DvrData, errStat, errMsg)

   type(Dvr_SimData),             intent(in)    :: DvrData
   integer,                       intent(  out) :: errStat           ! returns a non-zero value when an error occurs
   character(*),                  intent(  out) :: errMsg            ! Error message if errStat /= ErrID_None

   ! local variables:
   integer(intKi)                               :: i
   integer(intKi)                               :: FmtWidth          ! number of characters in string produced by DvrData%OutFmt
   integer(intKi)                               :: ErrStat2          ! temporary Error status
   character(ErrMsgLen)                         :: ErrMsg2           ! temporary Error message
   character(*), parameter                      :: RoutineName = 'ValidateInputs'



   ErrStat = ErrID_None
   ErrMsg  = ""


   ! Turbine Data:
   if ( DvrData%numBlades < 1 ) call SetErrStat( ErrID_Fatal, "There must be at least 1 blade (numBlades).", ErrStat, ErrMsg, RoutineName)
   if ( DvrData%numBlades > 3 ) call SetErrStat( ErrID_Fatal, "There can be no more than 3 blades (numBlades).", ErrStat, ErrMsg, RoutineName)
   if ( DvrData%HubRad < 0.0_ReKi .or. EqualRealNos(DvrData%HubRad, 0.0_ReKi) ) call SetErrStat( ErrID_Fatal, "HubRad must be a positive number.", ErrStat, ErrMsg, RoutineName)


   ! I-O Settings:
   ! Check that DvrData%OutFileData%OutFmt is a valid format specifier and will fit over the column headings
   call ChkRealFmtStr( DvrData%OutFileData%OutFmt, 'OutFmt', FmtWidth, ErrStat2, ErrMsg2 )
   call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   if ( FmtWidth /= ChanLen ) call SetErrStat( ErrID_Warn, 'OutFmt produces a column width of '// &
      TRIM(Num2LStr(FmtWidth))//' instead of '//TRIM(Num2LStr(ChanLen))//' characters.', ErrStat, ErrMsg, RoutineName )




   end subroutine ValidateInputs
   
   !----------------------------------------------------------------------------------------------------------------------------------
   subroutine Dvr_WriteOutputLine(OutFileData, t, output, errStat, errMsg)

   real(DbKi)             ,  intent(in   )   :: t                    ! simulation time (s)
   type(Dvr_OutputFile)   ,  intent(in   )   :: OutFileData
   real(ReKi)             ,  intent(in   )   :: output(:)            ! output data array
   integer(IntKi)         ,  intent(inout)   :: errStat              ! Status of error message
   character(*)           ,  intent(inout)   :: errMsg               ! Error message if ErrStat /= ErrID_None

   ! Local variables.

   character(200)                   :: frmt                                      ! A string to hold a format specifier
   character(15)                    :: tmpStr                                    ! temporary string to print the time output as text
   integer :: numOuts

   errStat = ErrID_None
   errMsg  = ''
   numOuts = size(output,1)
   frmt = '"'//OutFileData%delim//'"'//trim(OutFileData%outFmt)      ! format for array elements from individual modules

   ! time
   write( tmpStr, '(F15.4)' ) t
   call WrFileNR( OutFileData%unOutFile, tmpStr )
   call WrNumAryFileNR ( OutFileData%unOutFile, output,  frmt, errStat, errMsg )  ! write the outputs
   if ( errStat >= AbortErrLev ) return

   ! write a new line (advance to the next line)
   write (OutFileData%unOutFile,'()')

   end subroutine Dvr_WriteOutputLine
   !----------------------------------------------------------------------------------------------------------------------------------
   subroutine Dvr_InitializeOutputFile( OutFileData, errStat, errMsg)
   type(Dvr_OutputFile),     intent(inout)   :: OutFileData

   integer(IntKi)         ,  intent(inout)   :: errStat              ! Status of error message
   character(*)           ,  intent(inout)   :: errMsg               ! Error message if ErrStat /= ErrID_None

   ! locals
   integer(IntKi)                            ::  i
   integer(IntKi)                            :: numOuts
   character(50) 	        :: num_of_outs


   call GetNewUnit( OutFileData%unOutFile, ErrStat, ErrMsg )
   if ( ErrStat >= AbortErrLev ) then
      OutFileData%unOutFile = -1
      return
   end if

   call OpenFOutFile ( OutFileData%unOutFile, trim(outFileData%Root)//'.out', ErrStat, ErrMsg )
   if ( ErrStat >= AbortErrLev ) return

   write (OutFileData%unOutFile,'(/,A)')  'Predictions were generated on '//CurDate()//' at '//CurTime()//' using '//trim(GetNVD(version))
   write (OutFileData%unOutFile,'(1X,A)') trim(GetNVD(OutFileData%AD_ver))
   write (OutFileData%unOutFile,'()' )    !print a blank line

   write (OutFileData%unOutFile,'()' )    !print a blank line

   numOuts = size(OutFileData%WriteOutputHdr)
   !......................................................
   ! Write the names of the output parameters on one line:
   !......................................................

   call WrFileNR ( OutFileData%unOutFile, '     Time           ' )
   Write( num_of_outs, '(i10)' )  numOuts
   do i=1,NumOuts
      call WrFileNR ( OutFileData%unOutFile, OutFileData%delim//OutFileData%WriteOutputHdr(i) )
   end do ! i
   write (OutFileData%unOutFile,'()')

   !......................................................
   ! Write the units of the output parameters on one line:
   !......................................................
   call WrFileNR ( OutFileData%unOutFile, '      (s)           ' )
   do i=1,NumOuts
      call WrFileNR ( OutFileData%unOutFile, OutFileData%delim//OutFileData%WriteOutputUnt(i) )
   end do ! i
   write (OutFileData%unOutFile,'()')

   end subroutine Dvr_InitializeOutputFile
   !----------------------------------------------------------------------------------------------------------------------------------
   end module AeroDyn_Interface_Subs
