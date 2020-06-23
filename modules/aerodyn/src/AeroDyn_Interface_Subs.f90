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

   TYPE(ProgDesc), PARAMETER   :: version   = ProgDesc( 'AeroDyn_Interface_Subs', '', '' )  ! The version number of this program.
   
   contains

   !----------------------------------------------------------------------------------------------------------------------------------
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
   CALL DispCopyrightLicense( version ) ! djc: commented this out because it kept hanging on a Write statement
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
   subroutine Init_AeroDyn(DvrData, AD, hubPos, hubOri, hubVel, hubAcc, hubRotVel, hubRotAcc, bladePitch, dt, useAddedMass, errStat, errMsg)

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
   logical,                        intent(in   ) :: useAddedMass
   integer(IntKi),                 intent(  out) :: errStat       ! Status of error message
   character(*),                   intent(  out) :: errMsg        ! Error message if ErrStat /= ErrID_None

   ! locals
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
   call AD_Init(InitInData, AD%u(1), AD%p, AD%x, AD%xd, AD%z, AD%OtherState, AD%y, AD%m, dt, InitOutData, ErrStat2, ErrMsg2 )
   call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )

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
   DO j = -numInp, -1
      !extrapOri = hubOri - (hubRotVel * (dt * j) ) ! this won't work with all orientations
      extrapOri(:,:) = ExtrapOrientationFromRotVel(hubOri, hubRotVel, dt * j)
      call Set_AD_Inputs(dt * j,DvrData,AD,hubPos,extrapOri,hubVel,hubAcc,hubRotVel,hubRotAcc,bladePitch,errStat2,errMsg2)
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   END DO
   
   DvrData%iADStep = -1

   ! @djc TODO where should we actually set these?
   AD%p%IncludeAddedMass = useAddedMass
   AD%p%CaBlade = 1.0_ReKi
   
   ! move AD initOut data to AD Driver
   call move_alloc( InitOutData%WriteOutputHdr, DvrData%OutFileData%WriteOutputHdr )
   call move_alloc( InitOutData%WriteOutputUnt, DvrData%OutFileData%WriteOutputUnt )

   DvrData%OutFileData%AD_ver = InitOutData%ver

   contains

   function ExtrapOrientationFromRotVel(orientation, rotVel, time) result(orientation_prime)
      real(R8Ki), dimension(1:3,1:3), intent(in) :: orientation ! the global to local orientation matrix
      real(R8Ki), dimension(1:3),     intent(in) :: rotVel ! rotational velocity in axis-angle form
      real(DbKi),                     intent(in) :: time   ! time to extrapolate to

      real(ReKi), dimension(1:3)             :: theta_prime
      real(R8Ki), dimension(1:3,1:3)         :: orientation_prime
      real(R8Ki), dimension(1:3)             :: xBase, yBase, zBase, xBase_prime, yBase_prime, zBase_prime
      real(DbKi), dimension(1:3)             :: rotVel_x_time

      ! if time == 0, then just return the passed euler angles
      if( EqualRealNos(time, 0.0_DbKi) ) then
         orientation_prime = orientation
         return
      end if

      ! calculate and save this because we use it more than once
      rotVel_x_time = rotVel * time

      ! extract the basis vectors
      ! The matrix is global to local, so the rows are the body's basis vectors in global frame of reference
      xBase(:) = orientation(:,1)
      yBase(:) = orientation(:,2)
      zBase(:) = orientation(:,3)

      ! rotation each basis vector using rotation velocity
      xBase_prime = RotateVector(xBase, rotVel_x_time)
      yBase_prime = RotateVector(yBase, rotVel_x_time)
      zBase_prime = RotateVector(zBase, rotvel_x_time)

      ! create new orientation matrix from the rotated bases
      orientation_prime(:,1) = xBase_prime(:)
      orientation_prime(:,2) = yBase_prime(:)
      orientation_prime(:,3) = zBase_prime(:)
      
   end function ExtrapOrientationFromRotVel

   !> Implements Rodrigues' rotation formula
   function RotateVector(v, e) result(vprime)
   real(R8Ki), dimension(1:3), intent(in) :: v ! vector to rotate
   real(DbKi), dimension(1:3), intent(in) :: e ! scaled axis vector

   real(DbKi)                             :: theta
   real(ReKi), dimension(1:3)             :: e_hat
   real(ReKi), dimension(1:3)             :: vprime

   theta = norm2(e) ! magnatude of e
   
   ! to avoid division by zero
   if( EqualRealNos(theta, 0.0_DbKi) ) then
      vprime = v
      return
   end if
   
   e_hat = e / theta ! normalized axis

   vprime = cos(theta) * v + sin(theta) * cross_product(e_hat, v) * (1 - cos(theta)) * dot_product(e_hat,v) * e_hat
   end function RotateVector

   subroutine cleanup()
   call AD_DestroyInitInput( InitInData, ErrStat2, ErrMsg2 )
   call AD_DestroyInitOutput( InitOutData, ErrStat2, ErrMsg2 )
   end subroutine cleanup

   end subroutine Init_AeroDyn

   !----------------------------------------------------------------------
   !!< Hasn't been tested
   !subroutine CalcAddedMassMatrix(addedMassMatrix, m, p)
   !   real(ReKi), dimension(1:6,1:6), intent(inout) :: addedMassMatrix !< 
   !   type(AD_MiscVarType), intent(in)            :: m               !< Contains added mass data
   !   type(AD_ParameterType), intent(in)          :: p      
   !   ! local variables
   !   integer(IntKi) :: i, j, k ! iterators 
   !   real(ReKi)     :: magnitude
   !   real(ReKi), dimension(1:6) :: e ! Will hold a unit vector
   !   real(ReKi), dimension(1:6) :: totalAddedMassVec 
   !
   !   ! Make sure all entries are zero before we fill them in
   !   addedMassMatrix(:,:) = 0.0_ReKi
   !   
   !   do k=1, p%NumBlades
   !      do j=1, p%NumBlNds
   !         ! get the squared magnitude by calculating the dot product of the added mass vector with itself
   !         magnitude = dot_product(m%FAddedMass_vec(:,j,k), m%FAddedMass_vec(:,j,k))
   !         ! then get the magnitude by calculating the square root
   !         
   !         if ( .not. EqualRealNos(magnitude, 0.0_ReKi) ) then
   !            magnitude = sqrt(magnitude)
   !         
   !            ! Get the normalized added mass vector
   !            e = m%FAddedMass_vec(:,j,k) / magnitude
   !         
   !            totalAddedMassVec = totalAddedMassVec + m%AddedMass(j,k) * e
   !         
   !         end if
   !      end do ! numNodes
   !   end do ! numBlades
   !   
   !   do i=1, 6
   !      addedMassMatrix(i,i) = addedMassMatrix(i,i) + totalAddedMassVec(i)
   !   end do ! i=6
   !   
	  ! ! @mh: TODO off diagonal terms?
   !
   !   
   !end subroutine CalcAddedMassMatrix
   !

   !----------------------------------------------------------------------------------------------------------------------------------
   !> this routine returns time=(nt-1) * DvrData%Cases(iCase)%dT, and cycles values in the input array AD%InputTime and AD%u.
   !! it then sets the inputs for nt * DvrData%Cases(iCase)%dT, which are index values 1 in the arrays.
   subroutine Set_AD_Inputs(time,DvrData,AD,hubPos,hubOri,hubVel,hubAcc,hubRotVel,hubRotAcc,bladePitch,errStat,errMsg)

   real(DbKi),                         intent(in   ) :: time          ! time of the inputs (changed from timestep number because now using variable timestep)
   type(Dvr_SimData),                  intent(inout) :: DvrData       ! Driver data
   type(AeroDyn_Data),                 intent(inout) :: AD            ! AeroDyn data
   real(C_DOUBLE), dimension(1:3),     intent(in   ) :: hubPos        ! x,y,z in meters
   real(C_DOUBLE), dimension(1:3,1:3), intent(in   ) :: hubOri        ! the hub's orientation matrix (local to global)
   real(C_DOUBLE), dimension(1:3),     intent(in   ) :: hubVel        ! x,y,z in meters/sec
   real(C_DOUBLE), dimension(1:3),     intent(in   ) :: hubAcc        ! x,y,z in meters/sec^2
   real(C_DOUBLE), dimension(1:3),     intent(in   ) :: hubRotVel     ! axis-angle in global coordinate system, (magnitute in radians/sec)
   real(C_DOUBLE), dimension(1:3),     intent(in   ) :: hubRotAcc     ! axis-angle in global coordinate system, (magnitute in radians/sec^2)
   real(C_DOUBLE),                     intent(in   ) :: bladePitch
   integer(IntKi),                     intent(  out) :: errStat       ! Status of error message
   character(*),                       intent(  out) :: errMsg        ! Error message if ErrStat /= ErrID_None

   ! local variables
   integer(IntKi)                              :: errStat2      ! local status of error message
   character(ErrMsgLen)                        :: errMsg2       ! local error message if ErrStat /= ErrID_None
   character(*), parameter                     :: RoutineName = 'Set_AD_Inputs'

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
   real(C_DOUBLE), dimension(3)                :: vec_hold
   real(C_DOUBLE), dimension(3,3)              :: mat_hold
   real(DbKi)                                  :: ad_time
   real(DbKi), parameter                       :: epsilon = 1.0e-6

   errStat = ErrID_None
   errMsg  = ""

   ! note that this initialization is a little different than the general algorithm in FAST because here
   ! we can get exact values, so we are going to ignore initial guesses and not extrapolate

   
   !if ( time <= AD%InputTime(1) .OR. EqualRealNos(time, AD%InputTime(1) )) then
      !return 
   !end if
   
   !................
   ! shift previous calculations:
   !................
   do j = numInp-1,1,-1
      call AD_CopyInput (AD%u(j),  AD%u(j+1),  MESH_UPDATECOPY, ErrStat2, ErrMsg2)
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

      AD%InputTime(j+1) = AD%InputTime(j)
   end do
   AD%inputTime(1) = time

   !................
   ! calculate new values
   !................

   ! Tower motions:
   do j=1,AD%u(1)%TowerMotion%nnodes
      AD%u(1)%TowerMotion%Orientation(  :,:,j) = AD%u(1)%TowerMotion%RefOrientation(:,:,j) ! identity
      AD%u(1)%TowerMotion%TranslationDisp(:,j) = 0.0_ReKi
      AD%u(1)%TowerMotion%TranslationVel( :,j) = 0.0_ReKi
   end do !j=nnodes

   !.................
   ! Debug stuff (Visual studio's debugger doesn't display the values of non-local variables sometimes, but assigning them
   !              to a local variable allows us to see them)
   vec_hold = AD%u(1)%HubMotion%RotationVel(:,1)
   mat_hold = AD%u(1)%HubMotion%RefOrientation(:,:,1)
   !.................

   ! Hub motions:

   AD%u(1)%HubMotion%Position(:,1) = HubPos
   
   ! This would usually express where the hub is rotated to around the tower. Since we want hub to rotate around its center -
   ! i.e. not change position upon rotation - we just set this to zero.
   AD%u(1)%HubMotion%TranslationDisp(:,1) = 0.0

   !AD%u(1)%HubMotion%TranslationDisp(:,1) = 0.0 !matmul( AD%u(1)%HubMotion%Position(:,1), orientation ) - AD%u(1)%HubMotion%Position(:,1) ! = matmul( transpose(orientation) - eye(3), AD%u(1)%HubMotion%Position(:,1) )

   ! Rather than pass the Euler angles we pass the orientation matrix itself
   ! Note that HubMotion's orientation needs to be global to local, and our parameter
   ! is local to global, so we do the transpose here
   AD%u(1)%HubMotion%Orientation(  :,:,1) = transpose(hubOri(:,:))
 
   ! set the rotational velocity directly from the argument (expected to be global axis-angle rotational vel)
   AD%u(1)%HubMotion%RotationVel(    :,1) = hubRotVel
   AD%u(1)%HubMotion%RotationAcc(    :,1) = hubRotAcc

   AD%u(1)%HubMotion%TranslationVel(:,1) = hubVel
   AD%u(1)%HubMotion%TranslationAcc(:,1) = hubAcc


   ! Blade root motions:

   do k=1,DvrData%numBlades
      theta(1) = (k-1)*TwoPi/real(DvrData%numBlades,ReKi)
      theta(2) =  DvrData%precone
      theta(3) = -bladePitch
      orientation = EulerConstruct(theta)

      AD%u(1)%BladeRootMotion(k)%Orientation(  :,:,1) = matmul( orientation, AD%u(1)%HubMotion%Orientation(  :,:,1) )

   end do !k=numBlades

   ! Blade motions:
   do k=1,DvrData%numBlades
      rotateMat = transpose( AD%u(1)%BladeRootMotion(k)%Orientation(  :,:,1) )
      rotateMat = matmul( rotateMat, AD%u(1)%BladeRootMotion(k)%RefOrientation(  :,:,1) )
      orientation = transpose(rotateMat)

      rotateMat(1,1) = rotateMat(1,1) - 1.0_ReKi
      rotateMat(2,2) = rotateMat(2,2) - 1.0_ReKi
      rotateMat(3,3) = rotateMat(3,3) - 1.0_ReKi

      do j=1,AD%u(1)%BladeMotion(k)%nnodes
         position = AD%u(1)%BladeMotion(k)%Position(:,j)
         AD%u(1)%BladeMotion(k)%TranslationDisp(:,j) = matmul( rotateMat, position ) + HubPos

         AD%u(1)%BladeMotion(k)%Orientation(  :,:,j) = matmul( AD%u(1)%BladeMotion(k)%RefOrientation(:,:,j), orientation )

         position =  AD%u(1)%BladeMotion(k)%Position(:,j) + AD%u(1)%BladeMotion(k)%TranslationDisp(:,j) &
            - AD%u(1)%HubMotion%Position(:,1) - AD%u(1)%HubMotion%TranslationDisp(:,1) ! BM%Position will always be from origin, so we don't need to subtract HM%Position
            
         rotVel_cross_offset = cross_product( AD%u(1)%HubMotion%RotationVel(:,1), position )
         rotAcc_cross_offset = cross_product( AD%u(1)%HubMotion%RotationAcc(:,1), position )
         AD%u(1)%BladeMotion(k)%TranslationVel( :,j) = rotVel_cross_offset + hubVel ! add hub vel because hub can a linear velocity
         AD%u(1)%BladeMotion(k)%TranslationAcc( :,j) = rotAcc_cross_offset + rotVel_cross_offset + hubAcc

      end do !j=nnodes

   end do !k=numBlades

   ! @dustin: note, I took out code that sets the blade node inflows which usually goes here. I did this because we're setting it
   ! through an outside source (ProteusDS) and the outside source can't know where the nodes are until
   ! this subroutine ends. So I introduced a new subroutine to be called during initialization called Set_AD_Inflows
   end subroutine Set_AD_Inputs
   
   
   !---------------------------------------------------------------------------------------------------------------------------------
   !> djc: I took the code that sets the inflows out of Set_AD_Inputs because Set_AD_Inputs is called from the
   !! Init_Aerodyn subroutine, and we don't know the number of blades and nodes until that subroutine is done.
   !! So we wait for it to be done, allow PDS to get the number of nodes and blades, allocate inflows accordingly,
   !! get the node positions, and then initialize the inflows using this subroutine
   !! Must be called directly after all calls to Set_AD_Inputs.
   !! Refer to normal Aerodyn_Driver's Set_AD_Inputs to see how I'm trying to keep the end result the same after taking
   !! set inflows out of that subroutine
   SUBROUTINE Init_AD_Inflows(AD, nBlades, nNodes, bladeNodeInflowVel, bladeNodeInflowAcc)
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

   END SUBROUTINE Init_AD_Inflows

   !------------------------------------------------------------------------------------------------------------------------------
   !> djc: I took this out of Set_AD_Inputs because Set_AD_Inputs is called from the Init_Aerodyn
   !! subroutine, and we don't know the number of blades and nodes until that subroutine is done.
   !! So we wait for it to be done, allow PDS to get the number of nodes and blades, allocate inflows accordingly,
   !! get the node positions, and then initialize the inflows using this subroutine
   !! Must be called directly after all calls to Set_AD_Inputs.
   SUBROUTINE Set_AD_Inflows(AD, nBlades, nNodes, bladeNodeInflowVel, bladeNodeInflowAcc)
      type(AeroDyn_Data),                       intent(inout)  :: AD
      integer(IntKi),                              intent(in)  :: nBlades, nNodes
      real(ReKi), dimension(3,nNodes, nBlades), intent(in) :: bladeNodeInflowVel
      real(ReKi), dimension(3,nNodes, nBlades), intent(in) :: bladeNodeInflowAcc

      integer(IntKi) :: j, k

      do k = 1, nBlades
         do j = 1, nNodes
            AD%u(1)%InflowOnBlade(1:3,j,k) = bladeNodeInflowVel(1:3,j,k)
            AD%u(1)%InflowAccOnBlade(1:3,j,k) = bladeNodeInflowAcc(1:3,j,k)
         enddo !j=nnodes
      enddo !k=nblades

   END SUBROUTINE Set_AD_Inflows

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
