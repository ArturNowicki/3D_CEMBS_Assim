!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module assim_s_interior

!BOP
! !MODULE: assim_s_interior
!
! !DESCRIPTION:
!  Contains routines and variables used for determining the
!  interior potential temperature restoring.
!
! !REVISION HISTORY:
!  SVN:$Id: assim_s_interior.F90 12674 2008-10-31 22:21:32Z njn01 $
!
! !USES

   use kinds_mod
   use domain
   use constants
   use broadcast
   use io
   use forcing_tools
   use time_management
   use prognostic
   use grid
   use exit_mod

   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS:

   public :: init_assim_s_interior,      &
              get_assim_s_interior_data, &
              set_assim_s_interior

! !PUBLIC DATA MEMBERS:

   real (r8), public ::       &! public for use in restart
      assim_s_interior_interp_last , & ! time last interpolation was done
      assim_s_mask_interp_last  ! copy for mask

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  internal module variables
!
!-----------------------------------------------------------------------

   real (r8), dimension(:,:,:,:,:), allocatable :: &
      ASSIM_S_INTERIOR_DATA  ! data to restore interior pot temp towards

   real (r8), dimension(:,:,:), allocatable :: &
      S_RESTORE_RTAU  ! inverse restoring timescale for variable
                       ! interior restoring

   real (r8), dimension(:,:,:,:,:), allocatable :: &
      ASSIM_S_MASK ! mask of assimilated data
                           
   integer (int_kind), dimension(:,:,:), allocatable :: &
      S_RESTORE_MAX_LEVEL ! maximum level for applying variable
                           ! interior restoring

   real (r8), dimension(12) :: &
      assim_s_interior_data_time, &    !
      assim_s_mask_data_time    !copy for mask

   real (r8), dimension(20) :: &
      assim_s_interior_data_renorm  ! factors to convert to model units

   real (r8) ::                &
      assim_s_interior_data_inc,    &! time increment between values of forcing data
      assim_s_interior_data_next,   &! time to be used for the next value of forcing data that is needed
      assim_s_mask_data_next,   &! copy for mask
      assim_s_interior_data_update, &! time when new forcing value to be added to interpolation set
      assim_s_mask_data_update, &! copy for mask
      assim_s_interior_interp_inc,  &! time increment between interpolation
      assim_s_interior_interp_next, &! time next interpolation will be done
      assim_s_interior_restore_tau, &! restoring timescale (non-variable)
      assim_s_interior_restore_rtau  ! reciprocal of restoring timescale

   integer (int_kind) ::             &
      assim_s_interior_interp_order,      &! order of temporal interpolation
      assim_s_interior_data_time_min_loc, &! index of the third dimension of assim_s_interior_data_time containing the minimum forcing time
      assim_s_mask_data_time_min_loc, &!copy for mask
      assim_s_interior_restore_max_level



   character (char_len) ::     &
      assim_s_interior_data_type,   &! keyword for period of forcing data
      assim_s_interior_filename,    &! name of file conainting forcing data
      assim_s_interior_file_fmt,    &! format (bin or netcdf) of forcing file
      assim_s_interior_interp_freq, &! keyword for period of temporal interpolation
      assim_s_interior_interp_type, &!
      assim_s_interior_data_label,  &!
      assim_s_interior_formulation, &!
      assim_s_interior_restore_filename, &!
      assim_s_interior_restore_file_fmt, &
      assim_s_mask_filename,        &
      assim_s_mask_file_fmt

   character (char_len), dimension(:), allocatable :: &
      assim_s_interior_data_names    ! names for required input data fields

   integer (int_kind), dimension(:), allocatable :: &
      assim_s_interior_bndy_loc,    &! location and field type for
      assim_s_interior_bndy_type     !   ghost cell updates

   logical (log_kind) :: &
      assim_s_interior_variable_restore, &
      assim_s_interior_on

!EOC
!***********************************************************************

 contains

!***********************************************************************
!BOP
! !IROUTINE: init_assim_s_interior
! !INTERFACE:

 subroutine init_assim_s_interior

! !DESCRIPTION:
!  Initializes potential temperature interior forcing by either
!  calculating or reading in the 3D temperature.  Also performs
!  initial book-keeping concerning when new data is needed for
!  the temporal interpolation and when the forcing will need
!  to be updated.
!
! !REVISION HISTORY:
!  same as module

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      n,                 &! dummy loop index
      nml_error           ! namelist i/o error flag

   character (char_len) :: &
      forcing_filename,    &! full filename of forcing data file
      mask_filename,       &! full filename of mask file
      long_name             ! long name for input data field

   type (datafile) :: &
      assim_s_int_data_file, &  ! data file descriptor for interior pot temp data
      assim_s_mask_file      ! data file descriptor for assim mask data

   type (io_field_desc) :: &
      pt_data_in,   &       ! io field descriptor for input pot temp data
      mask_data_in          ! io field descriptor for input mask data
   type (io_dim) :: &
      i_dim, j_dim, &! dimension descriptors for horiz dims
      k_dim          ! dimension descriptor  for depth

   namelist /assim_s_interior_nml/ assim_s_interior_on, assim_s_interior_data_type, &
        assim_s_interior_data_inc,         assim_s_interior_interp_type,         &
        assim_s_interior_interp_freq,      assim_s_interior_interp_inc,          &
        assim_s_interior_restore_tau,      assim_s_interior_filename,            &
        assim_s_interior_file_fmt,         assim_s_interior_restore_max_level,   &
        assim_s_interior_data_renorm,      assim_s_interior_formulation,         &
        assim_s_interior_variable_restore, assim_s_interior_restore_filename,    &
        assim_s_interior_restore_file_fmt, &
        assim_s_mask_filename,             assim_s_mask_file_fmt

!-----------------------------------------------------------------------
!
!  read interior potential temperature restoring namelist input
!    after setting default values.
!
!-----------------------------------------------------------------------

   assim_s_interior_on             = .false.
   assim_s_interior_formulation    = 'restoring'
   assim_s_interior_data_type      = 'none'
   assim_s_interior_data_inc       = 1.e20_r8
   assim_s_interior_interp_type    = 'nearest'
   assim_s_interior_interp_freq    = 'never'
   assim_s_interior_interp_inc     = 1.e20_r8
   assim_s_interior_restore_tau    = 1.e20_r8
   assim_s_interior_filename       = 'unknown-assim_s_interior'
   assim_s_interior_file_fmt       = 'bin'
   assim_s_interior_restore_max_level = 0
   assim_s_interior_data_renorm       = c1
   assim_s_interior_variable_restore  = .false.
   assim_s_interior_restore_filename  = 'unknown-assim_s_interior_restore'
   assim_s_interior_restore_filename  = 'bin'
   assim_s_mask_filename       = 'unknown-assim_s_interior'
   assim_s_mask_file_fmt       = 'bin'

   if (my_task == master_task) then
      open (nml_in, file=nml_filename, status='old', iostat=nml_error)
      if (nml_error /= 0) then
         nml_error = -1
      else
         nml_error =  1
      endif
      !*** keep reading until find right namelist
      do while (nml_error > 0)
         read(nml_in, nml=assim_s_interior_nml,iostat=nml_error)
      end do
      if (nml_error == 0) close(nml_in)
   end if

   call broadcast_scalar(nml_error, master_task)
   if (nml_error /= 0) then
     call exit_POP(sigAbort,'ERROR: reading assim_s_interior_nml')
   endif

   call broadcast_scalar(assim_s_interior_on,                master_task)
   call broadcast_scalar(assim_s_interior_formulation,       master_task)
   call broadcast_scalar(assim_s_interior_data_type,         master_task)
   call broadcast_scalar(assim_s_interior_data_inc,          master_task)
   call broadcast_scalar(assim_s_interior_interp_type,       master_task)
   call broadcast_scalar(assim_s_interior_interp_freq,       master_task)
   call broadcast_scalar(assim_s_interior_interp_inc,        master_task)
   call broadcast_scalar(assim_s_interior_restore_tau,       master_task)
   call broadcast_scalar(assim_s_interior_filename,          master_task)
   call broadcast_scalar(assim_s_interior_file_fmt,          master_task)
   call broadcast_scalar(assim_s_interior_restore_max_level, master_task)
   call broadcast_scalar(assim_s_interior_variable_restore,  master_task)
   call broadcast_scalar(assim_s_interior_restore_filename,  master_task)
   call broadcast_scalar(assim_s_interior_restore_file_fmt,  master_task)
   call broadcast_array (assim_s_interior_data_renorm,       master_task)
   call broadcast_scalar (assim_s_mask_filename,              master_task)
   call broadcast_scalar (assim_s_mask_file_fmt,              master_task)

!-----------------------------------------------------------------------
!
!  convert data_type to 'monthly-calendar' if input is 'monthly'
!
!-----------------------------------------------------------------------

   if (assim_s_interior_data_type == 'monthly') &
       assim_s_interior_data_type = 'monthly-calendar'

!-----------------------------------------------------------------------
!
!  calculate inverse of restoring time scale and convert to seconds.
!
!-----------------------------------------------------------------------

   assim_s_interior_restore_rtau =c1/(seconds_in_day*assim_s_interior_restore_tau)

!-----------------------------------------------------------------------
!
!  convert interp_type to corresponding integer value.
!
!-----------------------------------------------------------------------

   select case (assim_s_interior_interp_type)
   case ('nearest')
     assim_s_interior_interp_order = 1

   case ('linear')
     assim_s_interior_interp_order = 2

   case ('4point')
     assim_s_interior_interp_order = 4

   case default
     call exit_POP(sigAbort, &
       'init_assim_s_interior: Unknown value for assim_s_interior_interp_type')

   end select

!-----------------------------------------------------------------------
!
!  set values of the interior PT array (ASSIM_S_INTERIOR_DATA)
!    depending on the type of the PT data.
!
!-----------------------------------------------------------------------

   select case (assim_s_interior_data_type)

   case ('none')

      !*** no interior forcing, therefore no interpolation in time
      !*** needed, nor are any new values to be used

      assim_s_interior_data_next = never
      assim_s_interior_data_update = never
      assim_s_interior_interp_freq = 'never'

   case ('annual')

      !*** annual mean climatological interior potential temperature
      !*** (read in from a file) that is constant in time, therefore
      !*** no new values will be needed.

      allocate(ASSIM_S_INTERIOR_DATA(nx_block,ny_block,km, &
                                max_blocks_clinic,1))

      allocate(assim_s_interior_data_names(1), &
               assim_s_interior_bndy_loc  (1), &
               assim_s_interior_bndy_type (1))

      ASSIM_S_INTERIOR_DATA = c0
      assim_s_interior_data_names(1) = 'TEMPERATURE'
      assim_s_interior_bndy_loc  (1) = field_loc_center
      assim_s_interior_bndy_type (1) = field_type_scalar

      forcing_filename = assim_s_interior_filename

      assim_s_int_data_file = construct_file(assim_s_interior_file_fmt,         &
                                   full_name=trim(forcing_filename),  &
                                   record_length=rec_type_dbl,        &
                                   recl_words=nx_global*ny_global)

      call data_set(assim_s_int_data_file, 'open_read')

      i_dim = construct_io_dim('i',nx_global)
      j_dim = construct_io_dim('j',ny_global)
      k_dim = construct_io_dim('k',km)

      pt_data_in = construct_io_field(trim(assim_s_interior_data_names(1)),  &
                             dim1=i_dim, dim2=j_dim, dim3=k_dim,        &
                             field_loc = assim_s_interior_bndy_loc(1),       &
                             field_type = assim_s_interior_bndy_type(1),     &
                             d3d_array = ASSIM_S_INTERIOR_DATA(:,:,:,:,1))

      call data_set (assim_s_int_data_file, 'define', pt_data_in)
      call data_set (assim_s_int_data_file, 'read',   pt_data_in)
      call data_set (assim_s_int_data_file, 'close')
      call destroy_io_field(pt_data_in)
      call destroy_file(assim_s_int_data_file)

      if (assim_s_interior_data_renorm(1) /= c1) &
         ASSIM_S_INTERIOR_DATA = ASSIM_S_INTERIOR_DATA*assim_s_interior_data_renorm(1)

      assim_s_interior_data_next = never
      assim_s_interior_data_update = never
      assim_s_interior_interp_freq = 'never'

      if (my_task == master_task) then
         write(stdout,blank_fmt)
         write(stdout,'(a31,a)') ' Interior PT Annual file read: ', &
                                 trim(forcing_filename)
     endif

   case ('monthly-equal','monthly-calendar')

      !*** monthly mean climatological interior potential temperature.
      !*** All 12 months are read in from a file. interpolation order
      !*** may be specified with namelist input.

      allocate(ASSIM_S_INTERIOR_DATA(nx_block,ny_block,km,   &
                                max_blocks_clinic,0:12))

      allocate(assim_s_interior_data_names(12), &
               assim_s_interior_bndy_loc  (12), &
               assim_s_interior_bndy_type (12))

      ASSIM_S_INTERIOR_DATA = c0
      call find_forcing_times(          assim_s_interior_data_time,         &
               assim_s_interior_data_inc,    assim_s_interior_interp_type,       &
               assim_s_interior_data_next,   assim_s_interior_data_time_min_loc, &
               assim_s_interior_data_update, assim_s_interior_data_type)

      forcing_filename = assim_s_interior_filename
      assim_s_int_data_file = construct_file(assim_s_interior_file_fmt,          &
                                    full_name=trim(forcing_filename),  &
                                    record_length=rec_type_dbl,        &
                                    recl_words=nx_global*ny_global)

      call data_set(assim_s_int_data_file, 'open_read')

      i_dim = construct_io_dim('i',nx_global)
      j_dim = construct_io_dim('j',ny_global)
      k_dim = construct_io_dim('k',km)

      do n=1,12
         write(assim_s_interior_data_names(n),'(a11,i2)') 'TEMPERATURE',n
         assim_s_interior_bndy_loc (n) = field_loc_center
         assim_s_interior_bndy_type(n) = field_type_scalar

         pt_data_in = construct_io_field(                               &
                             trim(assim_s_interior_data_names(n)),           &
                             dim1=i_dim, dim2=j_dim, dim3=k_dim,        &
                             field_loc = assim_s_interior_bndy_loc(n),       &
                             field_type = assim_s_interior_bndy_type(n),     &
                             d3d_array = ASSIM_S_INTERIOR_DATA(:,:,:,:,n))

         call data_set (assim_s_int_data_file, 'define', pt_data_in)
         call data_set (assim_s_int_data_file, 'read',   pt_data_in)
         call destroy_io_field(pt_data_in)
      enddo

      call data_set (assim_s_int_data_file, 'close')
      call destroy_file(assim_s_int_data_file)

      if (assim_s_interior_data_renorm(1) /= c1) &
         ASSIM_S_INTERIOR_DATA = ASSIM_S_INTERIOR_DATA*assim_s_interior_data_renorm(1)

      if (my_task == master_task) then
         write(stdout,blank_fmt)
         write(stdout,'(a32,a)') ' Interior PT Monthly file read: ', &
                                 trim(forcing_filename)
      endif

   case ('n-hour')

      !*** interior potential temperature specified every n-hours,
      !*** where the n-hour increment (assim_s_interior_data_inc) should
      !*** be specified with namelist input. only as many times as
      !*** are necessary based on the order of the temporal
      !*** interpolation scheme reside in memory at any given time.

      allocate(ASSIM_S_INTERIOR_DATA(nx_block,ny_block,km,max_blocks_clinic,&
                                          0:assim_s_interior_interp_order))

      allocate(ASSIM_S_MASK(nx_block,ny_block,km,max_blocks_clinic,&
                                          0:assim_s_interior_interp_order))

      allocate(assim_s_interior_data_names(1), &
               assim_s_interior_bndy_loc  (1), &
               assim_s_interior_bndy_type (1))

      ASSIM_S_INTERIOR_DATA = c0
      ASSIM_S_MASK = c0
      assim_s_interior_data_names(1) = 'TEMPERATURE'
      assim_s_interior_bndy_loc  (1) = field_loc_center
      assim_s_interior_bndy_type (1) = field_type_scalar

      call find_forcing_times(           assim_s_interior_data_time,        &
               assim_s_interior_data_inc,    assim_s_interior_interp_type,       &
               assim_s_interior_data_next,   assim_s_interior_data_time_min_loc, &
               assim_s_interior_data_update, assim_s_interior_data_type)

      do n = 1, assim_s_interior_interp_order

         call get_forcing_filename(mask_filename,         &
                                   assim_s_mask_filename,     &
                                   assim_s_interior_data_time(n), &
                                   assim_s_interior_data_inc)

         assim_s_mask_file = construct_file(assim_s_mask_file_fmt,       &
                                    full_name=trim(mask_filename),  &
                                    record_length=rec_type_dbl,        &
                                    recl_words=nx_global*ny_global)

         call data_set(assim_s_mask_file, 'open_read')

         i_dim = construct_io_dim('i',nx_global)
         j_dim = construct_io_dim('j',ny_global)
         k_dim = construct_io_dim('k',km)

         mask_data_in = construct_io_field(                               &
                             trim(assim_s_interior_data_names(1)),           &
                             dim1=i_dim, dim2=j_dim, dim3=k_dim,        &
                             field_loc = assim_s_interior_bndy_loc(1),       &
                             field_type = assim_s_interior_bndy_type(1),     &
                             d3d_array = ASSIM_S_MASK(:,:,:,:,n))

         call data_set (assim_s_mask_file, 'define', mask_data_in)
         call data_set (assim_s_mask_file, 'read',   mask_data_in)
         call data_set (assim_s_mask_file, 'close')
         call destroy_io_field(mask_data_in)
         call destroy_file(assim_s_mask_file)

         call get_forcing_filename(forcing_filename,         &
                                   assim_s_interior_filename,     &
                                   assim_s_interior_data_time(n), &
                                   assim_s_interior_data_inc)

         assim_s_int_data_file = construct_file(assim_s_interior_file_fmt,       &
                                    full_name=trim(forcing_filename),  &
                                    record_length=rec_type_dbl,        &
                                    recl_words=nx_global*ny_global)

         call data_set(assim_s_int_data_file, 'open_read')

         i_dim = construct_io_dim('i',nx_global)
         j_dim = construct_io_dim('j',ny_global)
         k_dim = construct_io_dim('k',km)

         pt_data_in = construct_io_field(                               &
                             trim(assim_s_interior_data_names(1)),           &
                             dim1=i_dim, dim2=j_dim, dim3=k_dim,        &
                             field_loc = assim_s_interior_bndy_loc(1),       &
                             field_type = assim_s_interior_bndy_type(1),     &
                             d3d_array = ASSIM_S_INTERIOR_DATA(:,:,:,:,n))

         call data_set (assim_s_int_data_file, 'define', pt_data_in)
         call data_set (assim_s_int_data_file, 'read',   pt_data_in)
         call data_set (assim_s_int_data_file, 'close')
         call destroy_io_field(pt_data_in)
         call destroy_file(assim_s_int_data_file)

         if (my_task == master_task) then
            write(stdout,blank_fmt)
            write(stdout,'(a31,a)') ' Interior PT n-hour file read: ', &
                                    trim(forcing_filename)
         endif
      enddo

      if (assim_s_interior_data_renorm(1) /= c1) &
         ASSIM_S_INTERIOR_DATA = ASSIM_S_INTERIOR_DATA*assim_s_interior_data_renorm(1)

   case default

     call exit_POP(sigAbort, &
       'init_assim_s_interior: Unknown value for assim_s_interior_data_type')

   end select

!-----------------------------------------------------------------------
!
!  now check interpolation period (assim_s_interior_interp_freq) to set
!    the time for the next temporal interpolation
!    (assim_s_interior_interp_next).
!
!  if no interpolation is to be done, set next interpolation time
!    to a large number so the interior PT update test in routine
!    set_surface_forcing will always be false.
!
!  if interpolation is to be done every n-hours, find the first
!    interpolation time greater than the current time.
!
!  if interpolation is to be done every timestep, set next interpolation
!    time to a large negative number so the interior PT update
!    test in routine set_surface_forcing will always be true.
!
!-----------------------------------------------------------------------

   select case (assim_s_interior_interp_freq)

   case ('never')

     assim_s_interior_interp_next = never
     assim_s_interior_interp_last = never
     assim_s_interior_interp_inc  = c0

   case ('n-hour')
     call find_interp_time(assim_s_interior_interp_inc, &
                           assim_s_interior_interp_next)

   case ('every-timestep')

     assim_s_interior_interp_next = always
     assim_s_interior_interp_inc  = c0

   case default

     call exit_POP(sigAbort, &
        'init_assim_s_interior: Unknown value for assim_s_interior_interp_freq')

   end select

   if(nsteps_total == 0) assim_s_interior_interp_last = thour00

!-----------------------------------------------------------------------
!
!  allocate and read in arrays used for variable interior restoring
!    if necessary.
!
!-----------------------------------------------------------------------

   if (assim_s_interior_variable_restore) then

      allocate(S_RESTORE_MAX_LEVEL(nx_block,ny_block,max_blocks_clinic), &
                    S_RESTORE_RTAU(nx_block,ny_block,max_blocks_clinic))

      forcing_filename = assim_s_interior_restore_filename

      assim_s_int_data_file = construct_file(assim_s_interior_restore_file_fmt,  &
                                    full_name=trim(forcing_filename),  &
                                    record_length=rec_type_dbl,        &
                                    recl_words=nx_global*ny_global)

      call data_set(assim_s_int_data_file, 'open_read')

      i_dim = construct_io_dim('i',nx_global)
      j_dim = construct_io_dim('j',ny_global)

      pt_data_in = construct_io_field('S_RESTORE_MAX_LEVEL',          &
                             dim1=i_dim, dim2=j_dim,                   &
                             field_loc = field_loc_center,             &
                             field_type = field_type_scalar,           &
                             d2d_array = S_RESTORE_RTAU)

      call data_set (assim_s_int_data_file, 'define', pt_data_in)
      call data_set (assim_s_int_data_file, 'read',   pt_data_in)
      S_RESTORE_MAX_LEVEL = nint(S_RESTORE_RTAU)
      call destroy_io_field(pt_data_in)

      pt_data_in = construct_io_field('S_RESTORE_RTAU', dim1=i_dim, dim2=j_dim, &
                             field_loc = field_loc_center,             &
                             field_type = field_type_scalar,           &
                             d2d_array = S_RESTORE_RTAU)

      call data_set (assim_s_int_data_file, 'define', pt_data_in)
      call data_set (assim_s_int_data_file, 'read',   pt_data_in)
      S_RESTORE_RTAU = S_RESTORE_RTAU/seconds_in_day ! convert days to secs
      call destroy_io_field(pt_data_in)
      call data_set (assim_s_int_data_file, 'close')
      call destroy_file(assim_s_int_data_file)

   endif

!-----------------------------------------------------------------------
!
!  echo forcing options to stdout.
!
!-----------------------------------------------------------------------

   assim_s_interior_data_label = 'Interior Potential Temperature Forcing'
   if (assim_s_interior_variable_restore .and. my_task == master_task) &
      write(stdout,'(a57)') &
      'Variable Interior Potential Temperature Restoring enabled'
   call echo_forcing_options(                 assim_s_interior_data_type,   &
                     assim_s_interior_formulation, assim_s_interior_data_inc,    &
                     assim_s_interior_interp_freq, assim_s_interior_interp_type, &
                     assim_s_interior_interp_inc,  assim_s_interior_data_label)

!-----------------------------------------------------------------------
!EOC

 end subroutine init_assim_s_interior

!***********************************************************************
!BOP
! !IROUTINE: get_assim_s_interior_data
! !INTERFACE:

 subroutine get_assim_s_interior_data

! !DESCRIPTION:
!  Determines whether new interior temperature forcing data is required
!  and reads the data if necessary.  Also interpolates data to current
!  time if required.
!
! !REVISION HISTORY:
!  same as module

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!  check if new data is necessary for interpolation.  if yes, then
!    shuffle time indices in assim_s_interior_data_time arrays
!    and read in new data if necessary ('n-hour' case).  also
!    increment values of assim_s_interior_data_time_min_loc,
!    assim_s_interior_data_next and assim_s_interior_data_update. note that no new
!    data is necessary for 'analytic' and 'annual' cases.
!
!-----------------------------------------------------------------------

   select case(assim_s_interior_data_type)

   case ('monthly-equal','monthly-calendar')

      assim_s_interior_data_label = 'assim_s_interior Monthly'
      if (thour00 >= assim_s_interior_data_update) then
         call update_forcing_data(            assim_s_interior_data_time,   &
               assim_s_interior_data_time_min_loc, assim_s_interior_interp_type, &
               assim_s_interior_data_next,         assim_s_interior_data_update, &
               assim_s_interior_data_type,         assim_s_interior_data_inc,    &
               ASSIM_S_INTERIOR_DATA(:,:,:,:,1:12),assim_s_interior_data_renorm, &
               assim_s_interior_data_label,        assim_s_interior_data_names,  &
               assim_s_interior_bndy_loc,          assim_s_interior_bndy_type,   &
               assim_s_interior_filename,          assim_s_interior_file_fmt)
      endif

      if (thour00 >= assim_s_interior_interp_next .or. nsteps_run==0) then
         call interpolate_forcing(ASSIM_S_INTERIOR_DATA(:,:,:,:,0),       &
                                  ASSIM_S_INTERIOR_DATA(:,:,:,:,1:12),    &
                   assim_s_interior_data_time, assim_s_interior_interp_type,   &
                   assim_s_interior_data_time_min_loc,                    &
                   assim_s_interior_interp_freq, assim_s_interior_interp_inc,  &
                   assim_s_interior_interp_next, assim_s_interior_interp_last, &
                   nsteps_run)

         if (nsteps_run /= 0) assim_s_interior_interp_next = &
                              assim_s_interior_interp_next + &
                              assim_s_interior_interp_inc
      endif

   case('n-hour')

      assim_s_interior_data_label = 'assim_s_interior n-hour'
      if (thour00 >= assim_s_interior_data_update) then
        assim_s_mask_data_time = assim_s_interior_data_time
        assim_s_mask_data_time_min_loc = assim_s_interior_data_time_min_loc
        assim_s_mask_data_next = assim_s_interior_data_next
        assim_s_mask_data_update = assim_s_interior_data_update
        assim_s_mask_interp_last = assim_s_interior_interp_last

         call update_forcing_data(            assim_s_interior_data_time,         &
               assim_s_interior_data_time_min_loc, assim_s_interior_interp_type, &
               assim_s_interior_data_next,         assim_s_interior_data_update, &
               assim_s_interior_data_type,         assim_s_interior_data_inc,    &
               ASSIM_S_INTERIOR_DATA(:,:,:,:,1:assim_s_interior_interp_order),   &
               assim_s_interior_data_renorm,                                      &
               assim_s_interior_data_label,        assim_s_interior_data_names,  &
               assim_s_interior_bndy_loc,          assim_s_interior_bndy_type,   &
               assim_s_interior_filename,          assim_s_interior_file_fmt)

         call update_forcing_data(            assim_s_mask_data_time,         &
               assim_s_mask_data_time_min_loc, assim_s_interior_interp_type, &
               assim_s_mask_data_next,         assim_s_mask_data_update, &
               assim_s_interior_data_type,         assim_s_interior_data_inc,    &
               ASSIM_S_MASK(:,:,:,:,1:assim_s_interior_interp_order),            &
               assim_s_interior_data_renorm,                                      &
               assim_s_interior_data_label,        assim_s_interior_data_names,                  &
               assim_s_interior_bndy_loc,          assim_s_interior_bndy_type,   &
               assim_s_mask_filename,              assim_s_mask_file_fmt)
      endif
      ! write(*,*) "AN2: forcing data",  &
      !   minval(ASSIM_S_INTERIOR_DATA(:,:,:,:,1:assim_s_interior_interp_order)), &
      !   maxval(ASSIM_S_INTERIOR_DATA(:,:,:,:,1:assim_s_interior_interp_order))
      ! write(*,*) "AN3: forcing mask",  &
      !   minval(ASSIM_S_MASK(:,:,:,:,1:assim_s_interior_interp_order)), &
      !   maxval(ASSIM_S_MASK(:,:,:,:,1:assim_s_interior_interp_order))
      if (thour00 >= assim_s_interior_interp_next .or. nsteps_run==0) then
         call interpolate_forcing(ASSIM_S_INTERIOR_DATA(:,:,:,:,0),         &
               ASSIM_S_INTERIOR_DATA(:,:,:,:,1:assim_s_interior_interp_order),   &
               assim_s_interior_data_time,         assim_s_interior_interp_type, &
               assim_s_interior_data_time_min_loc, assim_s_interior_interp_freq, &
               assim_s_interior_interp_inc,        assim_s_interior_interp_next, &
               assim_s_interior_interp_last,       nsteps_run)

         call interpolate_forcing(ASSIM_S_MASK(:,:,:,:,0),         &
               ASSIM_S_MASK(:,:,:,:,1:assim_s_interior_interp_order),   &
               assim_s_interior_data_time,         assim_s_interior_interp_type, &
               assim_s_interior_data_time_min_loc, assim_s_interior_interp_freq, &
               assim_s_interior_interp_inc,        assim_s_interior_interp_next, &
               assim_s_mask_interp_last,       nsteps_run)

         if (nsteps_run /= 0) assim_s_interior_interp_next = &
                              assim_s_interior_interp_next + &
                              assim_s_interior_interp_inc
      endif

   end select

!-----------------------------------------------------------------------
!EOC

 end subroutine get_assim_s_interior_data

!***********************************************************************
!BOP
! !IROUTINE: set_assim_s_interior
! !INTERFACE:

 subroutine set_assim_s_interior(k,this_block,S_SOURCE)

! !DESCRIPTION:
!  Computes the potential temperature restoring term using updated data.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      k     ! vertical level index

   type (block), intent(in) :: &
      this_block   ! block information for this block

! !OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block), intent(inout) :: &
      S_SOURCE    ! interior potential source term for this
                   ! block and level - accumulate restoring term
                   ! into this array

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) ::  &
      bid,                &! local block address for this block
      now                  ! index for interpolated data

   real (r8), dimension(nx_block,ny_block) :: &
      DASSIM_S_INTERIOR    ! interior potential restoring for this
                      ! block and level

!-----------------------------------------------------------------------
!
!  do interior restoring if required (no surface restoring for any)
!
!-----------------------------------------------------------------------

!   if (assim_s_interior_data_type /= 'none' .and. k > 1) then
   if (assim_s_interior_data_type /= 'none' .and. assim_s_interior_on) then

      bid = this_block%local_id

!-----------------------------------------------------------------------
!
!     set index where interpolated data located
!
!-----------------------------------------------------------------------

      select case(assim_s_interior_data_type)

      case('annual')
         now = 1

      case ('monthly-equal','monthly-calendar')
         now = 0

      case('n-hour')
         now = 0

      end select

!-----------------------------------------------------------------------
!
!     now compute restoring
!
!-----------------------------------------------------------------------

      if (assim_s_interior_variable_restore) then
         DASSIM_S_INTERIOR = S_RESTORE_RTAU(:,:,bid)*                &
                        merge((ASSIM_S_INTERIOR_DATA(:,:,k,bid,now) - &
                               TRACER(:,:,k,1,curtime,bid)),     &
                               c0, k <= S_RESTORE_MAX_LEVEL(:,:,bid))
      else
         if (k <= assim_s_interior_restore_max_level) then
            DASSIM_S_INTERIOR = assim_s_interior_restore_rtau*         &
                          (ASSIM_S_INTERIOR_DATA(:,:,k,bid,now) - &
                           TRACER(:,:,k,1,curtime,bid))
         else
            DASSIM_S_INTERIOR = c0
         endif
      endif
      ! *** add restoring to any other source terms
      where (ASSIM_S_MASK(:,:,k,bid,now) .gt. c0)
        S_SOURCE = S_SOURCE + 0.05
      endwhere

   endif ! k=1

!-----------------------------------------------------------------------
!EOC

 end subroutine set_assim_s_interior

!***********************************************************************

 end module assim_s_interior

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
