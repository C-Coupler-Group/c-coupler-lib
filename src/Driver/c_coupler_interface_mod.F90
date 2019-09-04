!*************************************************************
!  Copyright (c) 2017, Tsinghua University.
!  This is a source file of C-Coupler.
!  If you have any problem, 
!  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
!*************************************************************


 MODULE CCPL_interface_mod
   

   implicit none
   include 'mpif.h'
   private


   public :: CCPL_register_field_instance
   public :: CCPL_register_IO_field_from_data_buffer
   public :: CCPL_get_current_calendar_time
   public :: CCPL_register_V1D_grid_without_data
   public :: CCPL_register_V1D_Z_grid_via_model_data
   public :: CCPL_register_V1D_SIGMA_grid_via_model_data 
   public :: CCPL_register_V1D_HYBRID_grid_via_model_data
   public :: CCPL_register_H2D_grid_via_global_data 
   public :: CCPL_register_H2D_grid_via_local_data
   public :: CCPL_get_H2D_grid_data
   public :: CCPL_get_H2D_grid_area_in_remapping_wgts
   public :: CCPL_register_frac_based_remap_interface
   public :: CCPL_register_IO_field_from_field_instance 
   public :: CCPL_register_IO_fields_from_field_instances 
   public :: CCPL_get_number_of_current_step 
   public :: CCPL_get_number_of_total_steps 
   public :: CCPL_get_normal_time_step
   public :: CCPL_is_first_step
   public :: CCPL_is_first_restart_step
   public :: CCPL_get_current_num_days_in_year
   public :: CCPL_get_current_year 
   public :: CCPL_get_current_date 
   public :: CCPL_get_current_second 
   public :: CCPL_get_start_time
   public :: CCPL_get_stop_time
   public :: CCPL_get_previous_time
   public :: CCPL_get_current_time
   public :: CCPL_get_num_elapsed_days_from_reference 
   public :: CCPL_get_num_elapsed_days_from_start
   public :: CCPL_is_end_current_day 
   public :: CCPL_is_end_current_month
   public :: CCPL_allreduce_real16
   public :: CCPL_register_component
   public :: CCPL_get_comp_log_file_name
   public :: CCPL_get_comp_log_file_device
   public :: CCPL_get_component_id
   public :: CCPL_get_current_process_id_in_component
   public :: CCPL_is_comp_type_coupled
   public :: CCPL_get_component_process_global_id
   public :: CCPL_is_current_process_in_component
   public :: CCPL_get_num_process_in_component
   public :: CCPL_end_coupling_configuration
   public :: CCPL_register_CoR_defined_grid
   public :: CCPL_register_H2D_grid_via_file
   public :: CCPL_register_H2D_grid_from_another_component
   public :: CCPL_set_3D_grid_3D_vertical_coord_field
   public :: CCPL_set_3D_grid_variable_surface_field
   public :: CCPL_set_3D_grid_constant_surface_field
   public :: CCPL_set_3D_grid_external_surface_field
   public :: CCPL_register_MD_grid_via_multi_grids
   public :: CCPL_register_mid_point_grid
   public :: CCPL_get_grid_size
   public :: CCPL_get_grid_id
   public :: CCPL_register_normal_parallel_decomp 
   public :: CCPL_define_single_timer
   public :: CCPL_define_complex_timer 
   public :: CCPL_set_normal_time_step
   public :: CCPL_reset_current_time_to_start_time
   public :: CCPL_advance_time
   public :: CCPL_is_timer_on
   public :: CCPL_check_current_time
   public :: CCPL_finalize
   public :: CCPL_is_last_step_of_model_run
   public :: CCPL_is_model_run_ended
   public :: CCPL_register_normal_remap_interface
   public :: CCPL_register_import_interface 
   public :: CCPL_register_export_interface 
   public :: CCPL_execute_interface_using_id 
   public :: CCPL_execute_interface_using_name
   public :: CCPL_check_is_import_field_connected
   public :: CCPL_get_import_fields_sender_time
   public :: CCPL_get_local_comp_full_name 
   public :: CCPL_report_log 
   public :: CCPL_report_progress 
   public :: CCPL_report_error 
   public :: CCPL_do_restart_write_IO
   public :: CCPL_start_restart_read_IO
   public :: CCPL_is_restart_timer_on
   public :: CCPL_restart_read_fields_all
   public :: CCPL_restart_read_fields_interface
   public :: CCPL_get_restart_setting
   public :: CCPL_get_configurable_comps_full_names
   public :: CCPL_do_external_coupling_generation
   public :: CCPL_do_individual_coupling_generation
   public :: CCPL_do_family_coupling_generation
   public :: CCPL_abort


   interface CCPL_register_field_instance ; module procedure &
        CCPL_register_model_double_0D_data, &
        CCPL_register_model_double_1D_data, &
        CCPL_register_model_double_2D_data, &
        CCPL_register_model_double_3D_data, &
        CCPL_register_model_double_4D_data, &
        CCPL_register_model_float_0D_data, &
        CCPL_register_model_float_1D_data, &
        CCPL_register_model_float_2D_data, &
        CCPL_register_model_float_3D_data, &
        CCPL_register_model_float_4D_data, &
        CCPL_register_model_integer_0D_data, &
        CCPL_register_model_integer_1D_data, &
        CCPL_register_model_integer_2D_data, &
        CCPL_register_model_integer_3D_data, &
        CCPL_register_model_integer_4D_data
   end interface



   interface CCPL_register_IO_field_from_data_buffer ; module procedure &
        CCPL_register_new_IO_field_double_0D_data, &
        CCPL_register_new_IO_field_double_1D_data, &
        CCPL_register_new_IO_field_double_2D_data, &
        CCPL_register_new_IO_field_double_3D_data, &
        CCPL_register_new_IO_field_double_4D_data, &
        CCPL_register_new_IO_field_float_0D_data, &
        CCPL_register_new_IO_field_float_1D_data, &
        CCPL_register_new_IO_field_float_2D_data, &
        CCPL_register_new_IO_field_float_3D_data, &
        CCPL_register_new_IO_field_float_4D_data, &
        CCPL_register_new_IO_field_integer_0D_data, &
        CCPL_register_new_IO_field_integer_1D_data, &
        CCPL_register_new_IO_field_integer_2D_data, &
        CCPL_register_new_IO_field_integer_3D_data, &
        CCPL_register_new_IO_field_integer_4D_data
   end interface



   interface CCPL_get_current_calendar_time ; module procedure &
        CCPL_get_double_current_calendar_time, &
        CCPL_get_float_current_calendar_time
   end interface



   interface CCPL_register_V1D_Z_grid_via_model_data; module procedure &
        CCPL_register_V1D_Z_grid_via_double_data, &
        CCPL_register_V1D_Z_grid_via_float_data
   end interface



   interface CCPL_register_V1D_SIGMA_grid_via_model_data; module procedure &
        CCPL_register_V1D_SIGMA_grid_via_double_data, &
        CCPL_register_V1D_SIGMA_grid_via_float_data
   end interface



   interface CCPL_register_V1D_HYBRID_grid_via_model_data; module procedure &
        CCPL_register_V1D_HYBRID_grid_via_double_data, &
        CCPL_register_V1D_HYBRID_grid_via_float_data
   end interface



   interface CCPL_register_H2D_grid_via_global_data; module procedure &
        CCPL_register_H2D_grid_global_online_C1D_M1D_float, &
        CCPL_register_H2D_grid_global_online_C1D_M1D_double , &
        CCPL_register_H2D_grid_global_online_C2D_M2D_float, &
        CCPL_register_H2D_grid_global_online_C2D_M2D_double
   end interface



   interface CCPL_register_H2D_grid_via_local_data; module procedure &
        CCPL_register_H2D_grid_local_online_float, &
        CCPL_register_H2D_grid_local_online_double 
   end interface



   interface CCPL_get_H2D_grid_data; module procedure &
       CCPL_get_H2D_grid_integer_data, &
       CCPL_get_H2D_grid_float_data, &
       CCPL_get_H2D_grid_double_data
   end interface



   interface CCPL_get_H2D_grid_area_in_remapping_wgts; module procedure &
        CCPL_get_H2D_grid_float_area_in_remapping_wgts, &
        CCPL_get_H2D_grid_double_area_in_remapping_wgts
   end interface





   interface CCPL_register_frac_based_remap_interface; module procedure &
       CCPL_register_remap_interface_with_float_frac, &
       CCPL_register_remap_interface_with_double_frac
   end interface

   REAL,    parameter, public  :: coupling_fill_value = 1.0e30 
   integer, parameter, public  :: CCPL_NULL_INT       = 2147483647  
   integer, parameter, public  :: CCPL_NULL_COMM      = MPI_COMM_NULL  
   integer, parameter, public  :: CCPL_TAG_CPL        = 1
   integer, parameter, public  :: CCPL_TAG_REST       = 2
   integer, parameter, public  :: CCPL_TAG_CPL_REST   = 3


   integer, parameter, private :: R16 = selected_real_kind(24) ! 16-byte real
   integer, parameter, private :: R8  = selected_real_kind(12) ! 8-byte real
   integer, parameter, private :: R4  = selected_real_kind( 6) ! 4-byte real
   integer, parameter, private :: I2  = selected_int_kind(18)  ! 8-byte integer
   integer, parameter, private :: I4  = selected_int_kind(9)   ! 4-byte integer
   integer, parameter, private :: I8  = selected_int_kind(4)   ! 2-byte integer

   REAL(R16), allocatable, private :: reduce_buf_real16(:)
   integer,                private :: reduce_buf_real16_size



   interface check_CCPL_Fortran_API_int_type; module procedure &
       check_CCPL_Fortran_API_int_type_I2, &
       check_CCPL_Fortran_API_int_type_I4, &
       check_CCPL_Fortran_API_int_type_I8
   end interface



   CONTAINS 

!   
!  SUBROUTINE below
! 


   SUBROUTINE check_CCPL_Fortran_API_int_type_I2(API_int)
   implicit none
   integer(I2), INTENT(IN) ::  API_int

   call check_Fortran_API_int_type(2)
   
   END SUBROUTINE check_CCPL_Fortran_API_int_type_I2
  


   SUBROUTINE check_CCPL_Fortran_API_int_type_I4(API_int)
   implicit none
   integer(I4), INTENT(IN) ::  API_int

   call check_Fortran_API_int_type(4)
   
   END SUBROUTINE check_CCPL_Fortran_API_int_type_I4



   SUBROUTINE check_CCPL_Fortran_API_int_type_I8(API_int)
   implicit none
   integer(I8), INTENT(IN) ::  API_int

   call check_Fortran_API_int_type(8)
   
   END SUBROUTINE check_CCPL_Fortran_API_int_type_I8



   integer FUNCTION CCPL_register_model_double_0D_data(data_buf, field_name, decomp_id, comp_or_grid_id, buf_mark, usage_tag, field_unit, annotation)
   implicit none
   real(R8), INTENT(IN)                    :: data_buf
   character(len=*), intent(in)            :: field_name
   character(len=*), intent(in), optional  :: field_unit
   character(len=*), intent(in), optional  :: annotation
   integer,          intent(in)            :: decomp_id, comp_or_grid_id, buf_mark
   character *2048                         :: local_field_unit, local_annotation
   integer                                 :: field_instance_id
   integer,          intent(in), optional  :: usage_tag
   integer                                 :: local_usage_tag

   local_usage_tag = CCPL_TAG_CPL
   if (present(usage_tag)) local_usage_tag = usage_tag

   local_field_unit = "default_unit"
   if (present(field_unit)) then
       local_field_unit = field_unit
   endif
   local_annotation = ""
   if (present(annotation)) then
       local_annotation = annotation
   endif

   call register_external_field_instance(field_instance_id, trim(field_name)//char(0), loc(data_buf), 1, decomp_id, comp_or_grid_id, buf_mark, local_usage_tag, trim(local_field_unit)//char(0), trim("real8")//char(0), trim(local_annotation)//char(0))
   CCPL_register_model_double_0D_data = field_instance_id

   END FUNCTION CCPL_register_model_double_0D_data



   integer FUNCTION CCPL_register_model_double_1D_data(data_buf, field_name, decomp_id, comp_or_grid_id, buf_mark, usage_tag, field_unit, annotation)
   implicit none
   real(R8), INTENT(IN), DIMENSION(:)      :: data_buf
   character(len=*), intent(in)            :: field_name
   character(len=*), intent(in), optional  :: field_unit
   character(len=*), intent(in), optional  :: annotation
   integer,          intent(in)            :: decomp_id, comp_or_grid_id, buf_mark
   character *2048                         :: local_field_unit, local_annotation
   integer                                 :: field_instance_id
   integer,          intent(in), optional  :: usage_tag
   integer                                 :: local_usage_tag

   local_usage_tag = CCPL_TAG_CPL
   if (present(usage_tag)) local_usage_tag = usage_tag

   local_field_unit = "default_unit"
   if (present(field_unit)) then
       local_field_unit = field_unit
   endif
   local_annotation = ""
   if (present(annotation)) then
       local_annotation = annotation
   endif

   call register_external_field_instance(field_instance_id, trim(field_name)//char(0), loc(data_buf), size(data_buf), decomp_id, comp_or_grid_id, buf_mark, local_usage_tag, trim(local_field_unit)//char(0), trim("real8")//char(0), trim(local_annotation)//char(0))
   CCPL_register_model_double_1D_data = field_instance_id

   END FUNCTION CCPL_register_model_double_1D_data



   integer FUNCTION CCPL_register_model_double_2D_data(data_buf, field_name, decomp_id, comp_or_grid_id, buf_mark, usage_tag, field_unit, annotation)
   implicit none
   real(R8), INTENT(IN), DIMENSION(:,:)    :: data_buf
   character(len=*), intent(in)            :: field_name
   character(len=*), intent(in), optional  :: field_unit
   character(len=*), intent(in), optional  :: annotation
   integer,          intent(in)            :: decomp_id, comp_or_grid_id, buf_mark
   character *2048                         :: local_field_unit, local_annotation
   integer                                 :: field_instance_id
   integer,          intent(in), optional  :: usage_tag
   integer                                 :: local_usage_tag

   local_usage_tag = CCPL_TAG_CPL
   if (present(usage_tag)) local_usage_tag = usage_tag

   local_field_unit = "default_unit"
   if (present(field_unit)) then
       local_field_unit = field_unit
   endif
   local_annotation = ""
   if (present(annotation)) then
       local_annotation = annotation
   endif

   call register_external_field_instance(field_instance_id, trim(field_name)//char(0), loc(data_buf), size(data_buf), decomp_id, comp_or_grid_id, buf_mark, local_usage_tag, trim(local_field_unit)//char(0), trim("real8")//char(0), trim(local_annotation)//char(0))
   CCPL_register_model_double_2D_data = field_instance_id

   END FUNCTION CCPL_register_model_double_2D_data



   integer FUNCTION CCPL_register_model_double_3D_data(data_buf, field_name, decomp_id, comp_or_grid_id, buf_mark, usage_tag, field_unit, annotation)
   implicit none
   real(R8), INTENT(IN), DIMENSION(:,:,:)  :: data_buf
   character(len=*), intent(in)            :: field_name
   character(len=*), intent(in), optional  :: field_unit
   character(len=*), intent(in), optional  :: annotation
   integer,          intent(in)            :: decomp_id, comp_or_grid_id, buf_mark
   character *2048                         :: local_field_unit, local_annotation
   integer                                 :: field_instance_id
   integer,          intent(in), optional  :: usage_tag
   integer                                 :: local_usage_tag

   local_usage_tag = CCPL_TAG_CPL
   if (present(usage_tag)) local_usage_tag = usage_tag

   local_field_unit = "default_unit"
   if (present(field_unit)) then
       local_field_unit = field_unit
   endif
   local_annotation = ""
   if (present(annotation)) then
       local_annotation = annotation
   endif

   call register_external_field_instance(field_instance_id, trim(field_name)//char(0), loc(data_buf), size(data_buf), decomp_id, comp_or_grid_id, buf_mark, local_usage_tag, trim(local_field_unit)//char(0), trim("real8")//char(0), trim(local_annotation)//char(0))
   CCPL_register_model_double_3D_data = field_instance_id

   END FUNCTION CCPL_register_model_double_3D_data



   integer FUNCTION CCPL_register_model_double_4D_data(data_buf, field_name, decomp_id, comp_or_grid_id, buf_mark, usage_tag, field_unit, annotation)
   implicit none
   real(R8), INTENT(IN), DIMENSION(:,:,:,:)  :: data_buf
   character(len=*), intent(in)              :: field_name
   character(len=*), intent(in), optional    :: field_unit
   character(len=*), intent(in), optional    :: annotation
   integer,          intent(in)              :: decomp_id, comp_or_grid_id, buf_mark
   character *2048                           :: local_field_unit, local_annotation
   integer                                   :: field_instance_id
   integer,          intent(in), optional    :: usage_tag
   integer                                   :: local_usage_tag

   local_usage_tag = CCPL_TAG_CPL
   if (present(usage_tag)) local_usage_tag = usage_tag

   local_field_unit = "default_unit"
   if (present(field_unit)) then
       local_field_unit = field_unit
   endif
   local_annotation = ""
   if (present(annotation)) then
       local_annotation = annotation
   endif

   call register_external_field_instance(field_instance_id, trim(field_name)//char(0), loc(data_buf), size(data_buf), decomp_id, comp_or_grid_id, buf_mark, local_usage_tag, trim(local_field_unit)//char(0), trim("real8")//char(0), trim(local_annotation)//char(0))
   CCPL_register_model_double_4D_data = field_instance_id

   END FUNCTION CCPL_register_model_double_4D_data



   integer FUNCTION CCPL_register_model_float_0D_data(data_buf, field_name, decomp_id, comp_or_grid_id, buf_mark, usage_tag, field_unit, annotation)
   implicit none
   real(R4), INTENT(IN)                    :: data_buf
   character(len=*), intent(in)            :: field_name
   character(len=*), intent(in), optional  :: field_unit
   character(len=*), intent(in), optional  :: annotation
   integer,          intent(in)            :: decomp_id, comp_or_grid_id, buf_mark
   character *2048                         :: local_field_unit, local_annotation
   integer                                 :: field_instance_id
   integer,          intent(in), optional  :: usage_tag
   integer                                 :: local_usage_tag

   local_usage_tag = CCPL_TAG_CPL
   if (present(usage_tag)) local_usage_tag = usage_tag

   local_field_unit = "default_unit"
   if (present(field_unit)) then
       local_field_unit = field_unit
   endif
   local_annotation = ""
   if (present(annotation)) then
       local_annotation = annotation
   endif

   call register_external_field_instance(field_instance_id, trim(field_name)//char(0), loc(data_buf), 1, decomp_id, comp_or_grid_id, buf_mark, local_usage_tag, trim(local_field_unit)//char(0), trim("real4")//char(0), trim(local_annotation)//char(0))
   CCPL_register_model_float_0D_data = field_instance_id

   END FUNCTION CCPL_register_model_float_0D_data



   integer FUNCTION CCPL_register_model_float_1D_data(data_buf, field_name, decomp_id, comp_or_grid_id, buf_mark, usage_tag, field_unit, annotation)
   implicit none
   real(R4), INTENT(IN), DIMENSION(:)      :: data_buf
   character(len=*), intent(in)            :: field_name
   character(len=*), intent(in), optional  :: field_unit
   character(len=*), intent(in), optional  :: annotation
   integer,          intent(in)            :: decomp_id, comp_or_grid_id, buf_mark
   character *2048                         :: local_field_unit, local_annotation
   integer                                 :: field_instance_id
   integer,          intent(in), optional  :: usage_tag
   integer                                 :: local_usage_tag

   local_usage_tag = CCPL_TAG_CPL
   if (present(usage_tag)) local_usage_tag = usage_tag

   local_field_unit = "default_unit"
   if (present(field_unit)) then
       local_field_unit = field_unit
   endif
   local_annotation = ""
   if (present(annotation)) then
       local_annotation = annotation
   endif

   call register_external_field_instance(field_instance_id, trim(field_name)//char(0), loc(data_buf), size(data_buf), decomp_id, comp_or_grid_id, buf_mark, local_usage_tag, trim(local_field_unit)//char(0), trim("real4")//char(0), trim(local_annotation)//char(0))
   CCPL_register_model_float_1D_data = field_instance_id

   END FUNCTION CCPL_register_model_float_1D_data



   integer FUNCTION CCPL_register_model_float_2D_data(data_buf, field_name, decomp_id, comp_or_grid_id, buf_mark, usage_tag, field_unit, annotation)
   implicit none
   real(R4), INTENT(IN), DIMENSION(:,:)    :: data_buf
   character(len=*), intent(in)            :: field_name
   character(len=*), intent(in), optional  :: field_unit
   character(len=*), intent(in), optional  :: annotation
   integer,          intent(in)            :: decomp_id, comp_or_grid_id, buf_mark
   character *2048                         :: local_field_unit, local_annotation
   integer                                 :: field_instance_id
   integer,          intent(in), optional  :: usage_tag
   integer                                 :: local_usage_tag

   local_usage_tag = CCPL_TAG_CPL
   if (present(usage_tag)) local_usage_tag = usage_tag

   local_field_unit = "default_unit"
   if (present(field_unit)) then
       local_field_unit = field_unit
   endif
   local_annotation = ""
   if (present(annotation)) then
       local_annotation = annotation
   endif

   call register_external_field_instance(field_instance_id, trim(field_name)//char(0), loc(data_buf), size(data_buf), decomp_id, comp_or_grid_id, buf_mark, local_usage_tag, trim(local_field_unit)//char(0), trim("real4")//char(0), trim(local_annotation)//char(0))
   CCPL_register_model_float_2D_data = field_instance_id

   END FUNCTION CCPL_register_model_float_2D_data



   integer FUNCTION CCPL_register_model_float_3D_data(data_buf, field_name, decomp_id, comp_or_grid_id, buf_mark, usage_tag, field_unit, annotation)
   implicit none
   real(R4), INTENT(IN), DIMENSION(:,:,:)  :: data_buf
   character(len=*), intent(in)            :: field_name
   character(len=*), intent(in), optional  :: field_unit
   character(len=*), intent(in), optional  :: annotation
   integer,          intent(in)            :: decomp_id, comp_or_grid_id, buf_mark
   character *2048                         :: local_field_unit, local_annotation
   integer                                 :: field_instance_id
   integer,          intent(in), optional  :: usage_tag
   integer                                 :: local_usage_tag

   local_usage_tag = CCPL_TAG_CPL
   if (present(usage_tag)) local_usage_tag = usage_tag

   local_field_unit = "default_unit"
   if (present(field_unit)) then
       local_field_unit = field_unit
   endif
   local_annotation = ""
   if (present(annotation)) then
       local_annotation = annotation
   endif

   call register_external_field_instance(field_instance_id, trim(field_name)//char(0), loc(data_buf), size(data_buf), decomp_id, comp_or_grid_id, buf_mark, local_usage_tag, trim(local_field_unit)//char(0), trim("real4")//char(0), trim(local_annotation)//char(0))
   CCPL_register_model_float_3D_data = field_instance_id

   END FUNCTION CCPL_register_model_float_3D_data



   integer FUNCTION CCPL_register_model_float_4D_data(data_buf, field_name, decomp_id, comp_or_grid_id, buf_mark, usage_tag, field_unit, annotation)
   implicit none
   real(R4), INTENT(IN), DIMENSION(:,:,:,:)  :: data_buf
   character(len=*), intent(in)              :: field_name
   character(len=*), intent(in), optional    :: field_unit
   character(len=*), intent(in), optional    :: annotation
   integer,          intent(in)              :: decomp_id, comp_or_grid_id, buf_mark
   character *2048                           :: local_field_unit, local_annotation
   integer                                   :: field_instance_id
   integer,          intent(in), optional    :: usage_tag
   integer                                   :: local_usage_tag

   local_usage_tag = CCPL_TAG_CPL
   if (present(usage_tag)) local_usage_tag = usage_tag

   local_field_unit = "default_unit"
   if (present(field_unit)) then
       local_field_unit = field_unit
   endif
   local_annotation = ""
   if (present(annotation)) then
       local_annotation = annotation
   endif

   call register_external_field_instance(field_instance_id, trim(field_name)//char(0), loc(data_buf), size(data_buf), decomp_id, comp_or_grid_id, buf_mark, local_usage_tag, trim(local_field_unit)//char(0), trim("real4")//char(0), trim(local_annotation)//char(0))
   CCPL_register_model_float_4D_data = field_instance_id

   END FUNCTION CCPL_register_model_float_4D_data



   integer FUNCTION CCPL_register_model_integer_0D_data(data_buf, field_name, decomp_id, comp_or_grid_id, buf_mark, usage_tag, field_unit, annotation)
   implicit none
   integer, INTENT(IN)                     :: data_buf
   character(len=*), intent(in)            :: field_name
   character(len=*), intent(in), optional  :: field_unit
   character(len=*), intent(in), optional  :: annotation
   integer,          intent(in)            :: decomp_id, comp_or_grid_id, buf_mark
   character *2048                         :: local_field_unit, local_annotation
   integer                                 :: field_instance_id
   integer,          intent(in), optional  :: usage_tag
   integer                                 :: local_usage_tag

   local_usage_tag = CCPL_TAG_CPL
   if (present(usage_tag)) local_usage_tag = usage_tag

   local_field_unit = "default_unit"
   if (present(field_unit)) then
       local_field_unit = field_unit
   endif
   local_annotation = ""
   if (present(annotation)) then
       local_annotation = annotation
   endif

   call register_external_field_instance(field_instance_id, trim(field_name)//char(0), loc(data_buf), 1, decomp_id, comp_or_grid_id, buf_mark, local_usage_tag, trim(local_field_unit)//char(0), trim("integer")//char(0), trim(local_annotation)//char(0))
   CCPL_register_model_integer_0D_data = field_instance_id

   END FUNCTION CCPL_register_model_integer_0D_data



   integer FUNCTION CCPL_register_model_integer_1D_data(data_buf, field_name, decomp_id, comp_or_grid_id, buf_mark, usage_tag, field_unit, annotation)
   implicit none
   integer, INTENT(IN), DIMENSION(:)       :: data_buf
   character(len=*), intent(in)            :: field_name
   character(len=*), intent(in), optional  :: field_unit
   character(len=*), intent(in), optional  :: annotation
   integer,          intent(in)            :: decomp_id, comp_or_grid_id, buf_mark
   character *2048                         :: local_field_unit, local_annotation
   integer                                 :: field_instance_id
   integer,          intent(in), optional  :: usage_tag
   integer                                 :: local_usage_tag

   local_usage_tag = CCPL_TAG_CPL
   if (present(usage_tag)) local_usage_tag = usage_tag

   local_field_unit = "default_unit"
   if (present(field_unit)) then
       local_field_unit = field_unit
   endif
   local_annotation = ""
   if (present(annotation)) then
       local_annotation = annotation
   endif

   call register_external_field_instance(field_instance_id, trim(field_name)//char(0), loc(data_buf), size(data_buf), decomp_id, comp_or_grid_id, buf_mark, local_usage_tag, trim(local_field_unit)//char(0), trim("integer")//char(0), trim(local_annotation)//char(0))
   CCPL_register_model_integer_1D_data = field_instance_id

   END FUNCTION CCPL_register_model_integer_1D_data



   integer FUNCTION CCPL_register_model_integer_2D_data(data_buf, field_name, decomp_id, comp_or_grid_id, buf_mark, usage_tag, field_unit, annotation)
   implicit none
   integer, INTENT(IN), DIMENSION(:,:)     :: data_buf
   character(len=*), intent(in)            :: field_name
   character(len=*), intent(in), optional  :: field_unit
   character(len=*), intent(in), optional  :: annotation
   integer,          intent(in)            :: decomp_id, comp_or_grid_id, buf_mark
   character *2048                         :: local_field_unit, local_annotation
   integer                                 :: field_instance_id
   integer,          intent(in), optional  :: usage_tag
   integer                                 :: local_usage_tag

   local_usage_tag = CCPL_TAG_CPL
   if (present(usage_tag)) local_usage_tag = usage_tag

   local_field_unit = "default_unit"
   if (present(field_unit)) then
       local_field_unit = field_unit
   endif
   local_annotation = ""
   if (present(annotation)) then
       local_annotation = annotation
   endif

   call register_external_field_instance(field_instance_id, trim(field_name)//char(0), loc(data_buf), size(data_buf), decomp_id, comp_or_grid_id, buf_mark, local_usage_tag, trim(local_field_unit)//char(0), trim("integer")//char(0), trim(local_annotation)//char(0))
   CCPL_register_model_integer_2D_data = field_instance_id

   END FUNCTION CCPL_register_model_integer_2D_data



   integer FUNCTION CCPL_register_model_integer_3D_data(data_buf, field_name, decomp_id, comp_or_grid_id, buf_mark, usage_tag, field_unit, annotation)
   implicit none
   integer, INTENT(IN), DIMENSION(:,:,:)   :: data_buf
   character(len=*), intent(in)            :: field_name
   character(len=*), intent(in), optional  :: field_unit
   character(len=*), intent(in), optional  :: annotation
   integer,          intent(in)            :: decomp_id, comp_or_grid_id, buf_mark
   character *2048                         :: local_field_unit, local_annotation
   integer                                 :: field_instance_id
   integer,          intent(in), optional  :: usage_tag
   integer                                 :: local_usage_tag

   local_usage_tag = CCPL_TAG_CPL
   if (present(usage_tag)) local_usage_tag = usage_tag

   local_field_unit = "default_unit"
   if (present(field_unit)) then
       local_field_unit = field_unit
   endif
   local_annotation = ""
   if (present(annotation)) then
       local_annotation = annotation
   endif

   call register_external_field_instance(field_instance_id, trim(field_name)//char(0), loc(data_buf), size(data_buf), decomp_id, comp_or_grid_id, buf_mark, local_usage_tag, trim(local_field_unit)//char(0), trim("integer")//char(0), trim(local_annotation)//char(0))
   CCPL_register_model_integer_3D_data = field_instance_id

   END FUNCTION CCPL_register_model_integer_3D_data



   integer FUNCTION CCPL_register_model_integer_4D_data(data_buf, field_name, decomp_id, comp_or_grid_id, buf_mark, usage_tag, field_unit, annotation)
   implicit none
   integer, INTENT(IN), DIMENSION(:,:,:,:) :: data_buf
   character(len=*), intent(in)            :: field_name
   character(len=*), intent(in), optional  :: field_unit
   character(len=*), intent(in), optional  :: annotation
   integer,          intent(in)            :: decomp_id, comp_or_grid_id, buf_mark
   character *2048                         :: local_field_unit, local_annotation
   integer                                 :: field_instance_id
   integer,          intent(in), optional  :: usage_tag
   integer                                 :: local_usage_tag

   local_usage_tag = CCPL_TAG_CPL
   if (present(usage_tag)) local_usage_tag = usage_tag

   local_field_unit = "default_unit"
   if (present(field_unit)) then
       local_field_unit = field_unit
   endif
   local_annotation = ""
   if (present(annotation)) then
       local_annotation = annotation
   endif

   call register_external_field_instance(field_instance_id, trim(field_name)//char(0), loc(data_buf), size(data_buf), decomp_id, comp_or_grid_id, buf_mark, local_usage_tag, trim(local_field_unit)//char(0), trim("integer")//char(0), trim(local_annotation)//char(0))
   CCPL_register_model_integer_4D_data = field_instance_id

   END FUNCTION CCPL_register_model_integer_4D_data



   SUBROUTINE CCPL_register_IO_field_from_field_instance(field_inst_id, field_IO_name, annotation)
   implicit none
   integer,          intent(in)            :: field_inst_id
   character(len=*), intent(in), optional  :: field_IO_name
   character(len=*), intent(in), optional  :: annotation

   if (present(annotation)) then
       call register_an_io_field_from_field_instance(field_inst_id, trim(field_IO_name)//char(0), trim(annotation)//char(0))
   else
       call register_an_io_field_from_field_instance(field_inst_id, trim(field_IO_name)//char(0), trim("")//char(0))
   endif

   END SUBROUTINE CCPL_register_IO_field_from_field_instance



   SUBROUTINE CCPL_register_IO_fields_from_field_instances(num_field_inst, field_inst_ids, annotation)
   implicit none
   integer,          intent(in)                :: num_field_inst
   integer,          intent(in), dimension(:)  :: field_inst_ids
   character(len=*), intent(in), optional      :: annotation

   if (present(annotation)) then
       call register_io_fields_from_field_instances(num_field_inst, size(field_inst_ids), field_inst_ids, trim(annotation)//char(0))
   else
       call register_io_fields_from_field_instances(num_field_inst, size(field_inst_ids), field_inst_ids, trim("")//char(0))
   endif

   END SUBROUTINE CCPL_register_IO_fields_from_field_instances



   SUBROUTINE CCPL_register_new_IO_field_double_0D_data(data_buf, field_IO_name, &
              field_long_name, field_unit, decomp_id, comp_or_grid_id, annotation)
   implicit none
   real(R8), INTENT(IN)                    :: data_buf
   character(len=*), intent(in)            :: field_IO_name
   character(len=*), intent(in)            :: field_long_name
   character(len=*), intent(in)            :: field_unit
   integer,          intent(in)            :: comp_or_grid_id
   integer,          intent(in)            :: decomp_id
   character(len=*), intent(in), optional  :: annotation

   if (present(annotation)) then
       call register_a_new_io_field(comp_or_grid_id, decomp_id, 1, loc(data_buf), trim(field_IO_name)//char(0), trim(field_long_name)//char(0), &
                                     trim(field_unit)//char(0), trim("real8")//char(0), trim(annotation)//char(0))
   else
       call register_a_new_io_field(comp_or_grid_id, decomp_id, 1, loc(data_buf), trim(field_IO_name)//char(0), trim(field_long_name)//char(0), &
                                     trim(field_unit)//char(0), trim("real8")//char(0), trim("")//char(0))
   endif

   END SUBROUTINE CCPL_register_new_IO_field_double_0D_data 



   SUBROUTINE CCPL_register_new_IO_field_double_1D_data(data_buf, field_IO_name, &
              field_long_name, field_unit, decomp_id, comp_or_grid_id, annotation)
   implicit none
   real(R8), INTENT(IN), DIMENSION(:)      :: data_buf
   character(len=*), intent(in)            :: field_IO_name
   character(len=*), intent(in)            :: field_long_name
   character(len=*), intent(in)            :: field_unit
   integer,          intent(in)            :: comp_or_grid_id
   integer,          intent(in)            :: decomp_id
   character(len=*), intent(in), optional  :: annotation

   if (present(annotation)) then
       call register_a_new_io_field(comp_or_grid_id, decomp_id, size(data_buf), loc(data_buf), trim(field_IO_name)//char(0), trim(field_long_name)//char(0), &
                                     trim(field_unit)//char(0), trim("real8")//char(0), trim(annotation)//char(0))
   else
       call register_a_new_io_field(comp_or_grid_id, decomp_id, size(data_buf), loc(data_buf), trim(field_IO_name)//char(0), trim(field_long_name)//char(0), &
                                     trim(field_unit)//char(0), trim("real8")//char(0), trim("")//char(0))
   endif

   END SUBROUTINE CCPL_register_new_IO_field_double_1D_data 



   SUBROUTINE CCPL_register_new_IO_field_double_2D_data(data_buf, field_IO_name, &
              field_long_name, field_unit, decomp_id, comp_or_grid_id, annotation)
   implicit none
   real(R8), INTENT(IN), DIMENSION(:,:)    :: data_buf
   character(len=*), intent(in)            :: field_IO_name
   character(len=*), intent(in)            :: field_long_name
   character(len=*), intent(in)            :: field_unit
   integer,          intent(in)            :: comp_or_grid_id
   integer,          intent(in)            :: decomp_id
   character(len=*), intent(in), optional  :: annotation

   if (present(annotation)) then
       call register_a_new_io_field(comp_or_grid_id, decomp_id, size(data_buf), loc(data_buf), trim(field_IO_name)//char(0), trim(field_long_name)//char(0), &
                                     trim(field_unit)//char(0), trim("real8")//char(0), trim(annotation)//char(0))
   else
       call register_a_new_io_field(comp_or_grid_id, decomp_id, size(data_buf), loc(data_buf), trim(field_IO_name)//char(0), trim(field_long_name)//char(0), &
                                     trim(field_unit)//char(0), trim("real8")//char(0), trim("")//char(0))
   endif

   END SUBROUTINE CCPL_register_new_IO_field_double_2D_data 



   SUBROUTINE CCPL_register_new_IO_field_double_3D_data(data_buf, field_IO_name, &
              field_long_name, field_unit, decomp_id, comp_or_grid_id, annotation)
   implicit none
   real(R8), INTENT(IN), DIMENSION(:,:,:)  :: data_buf
   character(len=*), intent(in)            :: field_IO_name
   character(len=*), intent(in)            :: field_long_name
   character(len=*), intent(in)            :: field_unit
   integer,          intent(in)            :: comp_or_grid_id
   integer,          intent(in)            :: decomp_id
   character(len=*), intent(in), optional  :: annotation

   if (present(annotation)) then
       call register_a_new_io_field(comp_or_grid_id, decomp_id, size(data_buf), loc(data_buf), trim(field_IO_name)//char(0), trim(field_long_name)//char(0), &
                                     trim(field_unit)//char(0), trim("real8")//char(0), trim(annotation)//char(0))
   else
       call register_a_new_io_field(comp_or_grid_id, decomp_id, size(data_buf), loc(data_buf), trim(field_IO_name)//char(0), trim(field_long_name)//char(0), &
                                     trim(field_unit)//char(0), trim("real8")//char(0), trim("")//char(0))
   endif

   END SUBROUTINE CCPL_register_new_IO_field_double_3D_data 



   SUBROUTINE CCPL_register_new_IO_field_double_4D_data(data_buf, field_IO_name, &
              field_long_name, field_unit, decomp_id, comp_or_grid_id, annotation)
   implicit none
   real(R8), INTENT(IN), DIMENSION(:,:,:,:):: data_buf
   character(len=*), intent(in)            :: field_IO_name
   character(len=*), intent(in)            :: field_long_name
   character(len=*), intent(in)            :: field_unit
   integer,          intent(in)            :: comp_or_grid_id
   integer,          intent(in)            :: decomp_id
   character(len=*), intent(in), optional  :: annotation

   if (present(annotation)) then
       call register_a_new_io_field(comp_or_grid_id, decomp_id, size(data_buf), loc(data_buf), trim(field_IO_name)//char(0), trim(field_long_name)//char(0), &
                                     trim(field_unit)//char(0), trim("real8")//char(0), trim(annotation)//char(0))
   else
       call register_a_new_io_field(comp_or_grid_id, decomp_id, size(data_buf), loc(data_buf), trim(field_IO_name)//char(0), trim(field_long_name)//char(0), &
                                     trim(field_unit)//char(0), trim("real8")//char(0), trim("")//char(0))
   endif

   END SUBROUTINE CCPL_register_new_IO_field_double_4D_data 



   subroutine CCPL_register_new_IO_field_float_0D_data(data_buf, field_io_name, &
              field_long_name, field_unit, decomp_id, comp_or_grid_id, annotation)
   implicit none
   real(r4), intent(in)                    :: data_buf
   character(len=*), intent(in)            :: field_io_name
   character(len=*), intent(in)            :: field_long_name
   character(len=*), intent(in)            :: field_unit
   integer,          intent(in)            :: comp_or_grid_id
   integer,          intent(in)            :: decomp_id
   character(len=*), intent(in), optional  :: annotation

   if (present(annotation)) then
       call register_a_new_io_field(comp_or_grid_id, decomp_id, 1, loc(data_buf), trim(field_io_name)//char(0), trim(field_long_name)//char(0), &
                                     trim(field_unit)//char(0), trim("real4")//char(0), trim(annotation)//char(0))
   else
       call register_a_new_io_field(comp_or_grid_id, decomp_id, 1, loc(data_buf), trim(field_io_name)//char(0), trim(field_long_name)//char(0), &
                                     trim(field_unit)//char(0), trim("real4")//char(0), trim("")//char(0))
   endif

   end subroutine CCPL_register_new_IO_field_float_0D_data 



   SUBROUTINE CCPL_register_new_IO_field_float_1D_data(data_buf, field_IO_name, &
              field_long_name, field_unit, decomp_id, comp_or_grid_id, annotation)
   implicit none
   real(R4), INTENT(IN), DIMENSION(:)      :: data_buf
   character(len=*), intent(in)            :: field_IO_name
   character(len=*), intent(in)            :: field_long_name
   character(len=*), intent(in)            :: field_unit
   integer,          intent(in)            :: comp_or_grid_id
   integer,          intent(in)            :: decomp_id
   character(len=*), intent(in), optional  :: annotation

   if (present(annotation)) then
       call register_a_new_io_field(comp_or_grid_id, decomp_id, size(data_buf), loc(data_buf), trim(field_IO_name)//char(0), trim(field_long_name)//char(0), &
                                     trim(field_unit)//char(0), trim("real4")//char(0), trim(annotation)//char(0))
   else
       call register_a_new_io_field(comp_or_grid_id, decomp_id, size(data_buf), loc(data_buf), trim(field_IO_name)//char(0), trim(field_long_name)//char(0), &
                                     trim(field_unit)//char(0), trim("real4")//char(0), trim("")//char(0))
   endif

   END SUBROUTINE CCPL_register_new_IO_field_float_1D_data 



   SUBROUTINE CCPL_register_new_IO_field_float_2D_data(data_buf, field_IO_name, &
              field_long_name, field_unit, decomp_id, comp_or_grid_id, annotation)
   implicit none
   real(R4), INTENT(IN), DIMENSION(:,:)    :: data_buf
   character(len=*), intent(in)            :: field_IO_name
   character(len=*), intent(in)            :: field_long_name
   character(len=*), intent(in)            :: field_unit
   integer,          intent(in)            :: comp_or_grid_id
   integer,          intent(in)            :: decomp_id
   character(len=*), intent(in), optional  :: annotation

   if (present(annotation)) then
       call register_a_new_io_field(comp_or_grid_id, decomp_id, size(data_buf), loc(data_buf), trim(field_IO_name)//char(0), trim(field_long_name)//char(0), &
                                     trim(field_unit)//char(0), trim("real4")//char(0), trim(annotation)//char(0))
   else
       call register_a_new_io_field(comp_or_grid_id, decomp_id, size(data_buf), loc(data_buf), trim(field_IO_name)//char(0), trim(field_long_name)//char(0), &
                                     trim(field_unit)//char(0), trim("real4")//char(0), trim("")//char(0))
   endif

   END SUBROUTINE CCPL_register_new_IO_field_float_2D_data 



   SUBROUTINE CCPL_register_new_IO_field_float_3D_data(data_buf, field_IO_name, &
              field_long_name, field_unit, decomp_id, comp_or_grid_id, annotation)
   implicit none
   real(R4), INTENT(IN), DIMENSION(:,:,:)  :: data_buf
   character(len=*), intent(in)            :: field_IO_name
   character(len=*), intent(in)            :: field_long_name
   character(len=*), intent(in)            :: field_unit
   integer,          intent(in)            :: comp_or_grid_id
   integer,          intent(in)            :: decomp_id
   character(len=*), intent(in), optional  :: annotation

   if (present(annotation)) then
       call register_a_new_io_field(comp_or_grid_id, decomp_id, size(data_buf), loc(data_buf), trim(field_IO_name)//char(0), trim(field_long_name)//char(0), &
                                     trim(field_unit)//char(0), trim("real4")//char(0), trim(annotation)//char(0))
   else
       call register_a_new_io_field(comp_or_grid_id, decomp_id, size(data_buf), loc(data_buf), trim(field_IO_name)//char(0), trim(field_long_name)//char(0), &
                                     trim(field_unit)//char(0), trim("real4")//char(0), trim("")//char(0))
   endif

   END SUBROUTINE CCPL_register_new_IO_field_float_3D_data 



   SUBROUTINE CCPL_register_new_IO_field_float_4D_data(data_buf, field_IO_name, &
              field_long_name, field_unit, decomp_id, comp_or_grid_id, annotation)
   implicit none
   real(R4), INTENT(IN), DIMENSION(:,:,:,:):: data_buf
   character(len=*), intent(in)            :: field_IO_name
   character(len=*), intent(in)            :: field_long_name
   character(len=*), intent(in)            :: field_unit
   integer,          intent(in)            :: comp_or_grid_id
   integer,          intent(in)            :: decomp_id
   character(len=*), intent(in), optional  :: annotation

   if (present(annotation)) then
       call register_a_new_io_field(comp_or_grid_id, decomp_id, size(data_buf), loc(data_buf), trim(field_IO_name)//char(0), trim(field_long_name)//char(0), &
                                     trim(field_unit)//char(0), trim("real4")//char(0), trim(annotation)//char(0))
   else
       call register_a_new_io_field(comp_or_grid_id, decomp_id, size(data_buf), loc(data_buf), trim(field_IO_name)//char(0), trim(field_long_name)//char(0), &
                                     trim(field_unit)//char(0), trim("real4")//char(0), trim("")//char(0))
   endif

   END SUBROUTINE CCPL_register_new_IO_field_float_4D_data 



   subroutine CCPL_register_new_IO_field_integer_0D_data(data_buf, field_io_name, &
              field_long_name, field_unit, decomp_id, comp_or_grid_id, annotation)
   implicit none
   integer, intent(in)                     :: data_buf
   character(len=*), intent(in)            :: field_io_name
   character(len=*), intent(in)            :: field_long_name
   character(len=*), intent(in)            :: field_unit
   integer,          intent(in)            :: comp_or_grid_id
   integer,          intent(in)            :: decomp_id
   character(len=*), intent(in), optional  :: annotation

   if (present(annotation)) then
       call register_a_new_io_field(comp_or_grid_id, decomp_id, 1, loc(data_buf), trim(field_io_name)//char(0), trim(field_long_name)//char(0), &
                                     trim(field_unit)//char(0), trim("integer")//char(0), trim(annotation)//char(0))
   else
       call register_a_new_io_field(comp_or_grid_id, decomp_id, 1, loc(data_buf), trim(field_io_name)//char(0), trim(field_long_name)//char(0), &
                                     trim(field_unit)//char(0), trim("integer")//char(0), trim("")//char(0))
   endif

   end subroutine CCPL_register_new_IO_field_integer_0D_data 



   SUBROUTINE CCPL_register_new_IO_field_integer_1D_data(data_buf, field_IO_name, &
              field_long_name, field_unit, decomp_id, comp_or_grid_id, annotation)
   implicit none
   integer, INTENT(IN), DIMENSION(:)       :: data_buf
   character(len=*), intent(in)            :: field_IO_name
   character(len=*), intent(in)            :: field_long_name
   character(len=*), intent(in)            :: field_unit
   integer,          intent(in)            :: comp_or_grid_id
   integer,          intent(in)            :: decomp_id
   character(len=*), intent(in), optional  :: annotation

   if (present(annotation)) then
       call register_a_new_io_field(comp_or_grid_id, decomp_id, size(data_buf), loc(data_buf), trim(field_IO_name)//char(0), trim(field_long_name)//char(0), &
                                     trim(field_unit)//char(0), trim("integer")//char(0), trim(annotation)//char(0))
   else
       call register_a_new_io_field(comp_or_grid_id, decomp_id, size(data_buf), loc(data_buf), trim(field_IO_name)//char(0), trim(field_long_name)//char(0), &
                                     trim(field_unit)//char(0), trim("integer")//char(0), trim("")//char(0))
   endif

   END SUBROUTINE CCPL_register_new_IO_field_integer_1D_data 



   SUBROUTINE CCPL_register_new_IO_field_integer_2D_data(data_buf, field_IO_name, &
              field_long_name, field_unit, decomp_id, comp_or_grid_id, annotation)
   implicit none
   integer, INTENT(IN), DIMENSION(:,:)     :: data_buf
   character(len=*), intent(in)            :: field_IO_name
   character(len=*), intent(in)            :: field_long_name
   character(len=*), intent(in)            :: field_unit
   integer,          intent(in)            :: comp_or_grid_id
   integer,          intent(in)            :: decomp_id
   character(len=*), intent(in), optional  :: annotation

   if (present(annotation)) then
       call register_a_new_io_field(comp_or_grid_id, decomp_id, size(data_buf), loc(data_buf), trim(field_IO_name)//char(0), trim(field_long_name)//char(0), &
                                     trim(field_unit)//char(0), trim("integer")//char(0), trim(annotation)//char(0))
   else
       call register_a_new_io_field(comp_or_grid_id, decomp_id, size(data_buf), loc(data_buf), trim(field_IO_name)//char(0), trim(field_long_name)//char(0), &
                                     trim(field_unit)//char(0), trim("integer")//char(0), trim("")//char(0))
   endif

   END SUBROUTINE CCPL_register_new_IO_field_integer_2D_data 



   SUBROUTINE CCPL_register_new_IO_field_integer_3D_data(data_buf, field_IO_name, &
              field_long_name, field_unit, decomp_id, comp_or_grid_id, annotation)
   implicit none
   integer, INTENT(IN), DIMENSION(:,:,:)   :: data_buf
   character(len=*), intent(in)            :: field_IO_name
   character(len=*), intent(in)            :: field_long_name
   character(len=*), intent(in)            :: field_unit
   integer,          intent(in)            :: comp_or_grid_id
   integer,          intent(in)            :: decomp_id
   character(len=*), intent(in), optional  :: annotation

   if (present(annotation)) then
       call register_a_new_io_field(comp_or_grid_id, decomp_id, size(data_buf), loc(data_buf), trim(field_IO_name)//char(0), trim(field_long_name)//char(0), &
                                     trim(field_unit)//char(0), trim("integer")//char(0), trim(annotation)//char(0))
   else
       call register_a_new_io_field(comp_or_grid_id, decomp_id, size(data_buf), loc(data_buf), trim(field_IO_name)//char(0), trim(field_long_name)//char(0), &
                                     trim(field_unit)//char(0), trim("integer")//char(0), trim("")//char(0))
   endif

   END SUBROUTINE CCPL_register_new_IO_field_integer_3D_data 



   SUBROUTINE CCPL_register_new_IO_field_integer_4D_data(data_buf, field_IO_name, &
              field_long_name, field_unit, decomp_id, comp_or_grid_id, annotation)
   implicit none
   integer, INTENT(IN), DIMENSION(:,:,:,:) :: data_buf
   character(len=*), intent(in)            :: field_IO_name
   character(len=*), intent(in)            :: field_long_name
   character(len=*), intent(in)            :: field_unit
   integer,          intent(in)            :: comp_or_grid_id
   integer,          intent(in)            :: decomp_id
   character(len=*), intent(in), optional  :: annotation

   if (present(annotation)) then
       call register_a_new_io_field(comp_or_grid_id, decomp_id, size(data_buf), loc(data_buf), trim(field_IO_name)//char(0), trim(field_long_name)//char(0), &
                                     trim(field_unit)//char(0), trim("integer")//char(0), trim(annotation)//char(0))
   else
       call register_a_new_io_field(comp_or_grid_id, decomp_id, size(data_buf), loc(data_buf), trim(field_IO_name)//char(0), trim(field_long_name)//char(0), &
                                     trim(field_unit)//char(0), trim("integer")//char(0), trim("")//char(0))
   endif

   END SUBROUTINE CCPL_register_new_IO_field_integer_4D_data 



 integer FUNCTION CCPL_get_number_of_current_step(comp_id, annotation)
   implicit none  
   integer,          intent(in) :: comp_id
   character(len=*), intent(in), optional :: annotation
   integer                      :: nstep

   if (present(annotation)) then
      call get_ccpl_current_number_of_step(comp_id, nstep, trim(annotation)//char(0))
   else
      call get_ccpl_current_number_of_step(comp_id, nstep, trim("")//char(0))
   endif
   CCPL_get_number_of_current_step = nstep

 END FUNCTION CCPL_get_number_of_current_step



 integer FUNCTION CCPL_get_number_of_total_steps(comp_id, annotation)
   implicit none  
   integer,          intent(in) :: comp_id
   character(len=*), intent(in), optional :: annotation
   integer                      :: nstep

   if (present(annotation)) then
      call get_ccpl_num_total_step(comp_id, nstep, trim(annotation)//char(0))
   else
      call get_ccpl_num_total_step(comp_id, nstep, trim("")//char(0))
   endif
   CCPL_get_number_of_total_steps = nstep

 END FUNCTION CCPL_get_number_of_total_steps



 integer FUNCTION CCPL_get_normal_time_step(comp_id, annotation)
   implicit none  
   integer, intent(in) :: comp_id
   character(len=*), intent(in), optional :: annotation
   integer             :: step_size

   if (present(annotation)) then
      call get_ccpl_time_step(comp_id, step_size, trim(annotation)//char(0))
   else 
      call get_ccpl_time_step(comp_id, step_size, trim("")//char(0))
   endif
   CCPL_get_normal_time_step = step_size

 END FUNCTION CCPL_get_normal_time_step



 logical FUNCTION CCPL_is_first_step(comp_id, annotation)
   implicit none
   integer                      :: is_first_step
   integer,          intent(in) :: comp_id
   character(len=*), intent(in), optional :: annotation

   if (present(annotation)) then
      call is_comp_first_step(comp_id, is_first_step, trim(annotation)//char(0))
   else
      call is_comp_first_step(comp_id, is_first_step, trim("")//char(0))
   endif
   if (is_first_step .eq. 1) then
      CCPL_is_first_step = .true.
   else
      CCPL_is_first_step = .false.
   endif

 END FUNCTION CCPL_is_first_step



 logical FUNCTION CCPL_is_first_restart_step(comp_id, annotation)
   implicit none
   integer                      :: is_first_restart_step
   integer,          intent(in) :: comp_id
   character(len=*), intent(in), optional :: annotation

   if (present(annotation)) then
      call is_comp_first_restart_step(comp_id, is_first_restart_step, trim(annotation)//char(0))
   else
      call is_comp_first_restart_step(comp_id, is_first_restart_step, trim("")//char(0))
   endif
   if (is_first_restart_step .eq. 1) then
      CCPL_is_first_restart_step = .true.
   else
      CCPL_is_first_restart_step = .false.
   endif

 END FUNCTION CCPL_is_first_restart_step



 integer FUNCTION CCPL_get_current_num_days_in_year(comp_id, annotation)
   implicit none
   integer,          intent(in) :: comp_id
   character(len=*), intent(in), optional :: annotation
   integer                      :: days

   if (present(annotation)) then
      call get_ccpl_current_num_days_in_year(comp_id, days, trim(annotation)//char(0))
   else
      call get_ccpl_current_num_days_in_year(comp_id, days, trim("")//char(0))
   endif
   CCPL_get_current_num_days_in_year = days
 END FUNCTION CCPL_get_current_num_days_in_year



 integer FUNCTION CCPL_get_current_year(comp_id, annotation)
   implicit none
   character(len=*), intent(in), optional :: annotation
   integer,          intent(in) :: comp_id
   integer                      :: year

   if (present(annotation)) then
      call get_ccpl_current_year(comp_id, year, trim(annotation)//char(0))
   else
      call get_ccpl_current_year(comp_id, year, trim("")//char(0))
   endif
   CCPL_get_current_year = year
 END FUNCTION CCPL_get_current_year



 integer FUNCTION CCPL_get_current_date(comp_id, annotation)
   implicit none
   character(len=*), intent(in), optional :: annotation
   integer,          intent(in) :: comp_id
   integer                      :: date

   if (present(annotation)) then
      call get_ccpl_current_date(comp_id, date, trim(annotation)//char(0))
   else
      call get_ccpl_current_date(comp_id, date, trim("")//char(0))
   endif
   CCPL_get_current_date = date
 END FUNCTION CCPL_get_current_date



 integer FUNCTION CCPL_get_current_second(comp_id, annotation)
   implicit none
   character(len=*), intent(in), optional :: annotation
   integer,          intent(in) :: comp_id
   integer                      :: second

   if (present(annotation)) then
      call get_ccpl_current_second(comp_id, second, trim(annotation)//char(0))
   else
      call get_ccpl_current_second(comp_id, second, trim("")//char(0))
   endif
   CCPL_get_current_second = second
 END FUNCTION CCPL_get_current_second



 SUBROUTINE CCPL_get_start_time(comp_id, year, month, day, second, annotation)
   implicit none
   character(len=*), intent(in), optional :: annotation
   integer,          intent(in) :: comp_id
   integer                      :: year, month, day, second

   if (present(annotation)) then
      call get_ccpl_start_time(comp_id, year, month, day, second, trim(annotation)//char(0))
   else
      call get_ccpl_start_time(comp_id, year, month, day, second, trim("")//char(0))
   endif
 END SUBROUTINE CCPL_get_start_time



 SUBROUTINE CCPL_get_stop_time(comp_id, year, month, day, second, annotation)
   implicit none
   character(len=*), intent(in), optional :: annotation
   integer,          intent(in) :: comp_id
   integer                      :: year, month, day, second

   if (present(annotation)) then
      call get_ccpl_stop_time(comp_id, year, month, day, second, trim(annotation)//char(0))
   else
      call get_ccpl_stop_time(comp_id, year, month, day, second, trim("")//char(0))
   endif

 END SUBROUTINE CCPL_get_stop_time



 SUBROUTINE CCPL_get_previous_time(comp_id, year, month, day, second, annotation)
   implicit none
   character(len=*), intent(in), optional :: annotation
   integer,          intent(in) :: comp_id
   integer                      :: year, month, day, second

   if (present(annotation)) then
      call get_ccpl_previous_time(comp_id, year, month, day, second, trim(annotation)//char(0))
   else
      call get_ccpl_previous_time(comp_id, year, month, day, second, trim("")//char(0))
   endif

 END SUBROUTINE CCPL_get_previous_time



 SUBROUTINE CCPL_get_current_time(comp_id, year, month, day, second, shift_second, annotation)
   implicit none 
   character(len=*), intent(in), optional :: annotation
   integer,          intent(in) :: comp_id
   integer                      :: year, month, day, second
   integer,          optional   :: shift_second
   integer                      :: local_shift

   local_shift = 0
   if (present(shift_second)) local_shift = shift_second
   if (present(annotation)) then
      call get_ccpl_current_time(comp_id, year, month, day, second, local_shift, trim(annotation)//char(0))
   else 
      call get_ccpl_current_time(comp_id, year, month, day, second, local_shift, trim("")//char(0))
   endif

 END SUBROUTINE CCPL_get_current_time



 SUBROUTINE CCPL_get_num_elapsed_days_from_reference(comp_id, num_days, current_second, annotation)
   implicit none
   character(len=*), intent(in), optional :: annotation
   integer,          intent(in)           :: comp_id
   integer,          intent(out)          :: num_days, current_second

   if (present(annotation)) then
      call get_ccpl_num_elapsed_days_from_reference_date(comp_id, num_days, current_second, trim(annotation)//char(0))
   else
      call get_ccpl_num_elapsed_days_from_reference_date(comp_id, num_days, current_second, trim("")//char(0))
   endif

 END SUBROUTINE CCPL_get_num_elapsed_days_from_reference



 SUBROUTINE CCPL_get_num_elapsed_days_from_start(comp_id, num_days, current_second, annotation)
   implicit none
   character(len=*), intent(in), optional :: annotation
   integer,          intent(in)           :: comp_id
   integer,          intent(out)          :: num_days, current_second


   if (present(annotation)) then
      call get_ccpl_num_elapsed_days_from_start_date(comp_id, num_days, current_second, trim(annotation)//char(0))
   else
      call get_ccpl_num_elapsed_days_from_start_date(comp_id, num_days, current_second, trim("")//char(0))
   endif

 END SUBROUTINE CCPL_get_num_elapsed_days_from_start



 logical FUNCTION CCPL_is_end_current_day(comp_id, annotation)
   implicit none
   character(len=*), intent(in), optional :: annotation
   integer,          intent(in) :: comp_id
   integer                      :: year, month, day, second

   if (present(annotation)) then
      call get_ccpl_current_time(comp_id, year, month, day, second, 0, trim(annotation)//char(0))
   else
      call get_ccpl_current_time(comp_id, year, month, day, second, 0, trim("")//char(0))
   endif
   CCPL_is_end_current_day = (second == 0)

 END FUNCTION CCPL_is_end_current_day



 logical FUNCTION CCPL_is_end_current_month(comp_id, annotation)
   implicit none
   character(len=*), intent(in), optional :: annotation
   integer,          intent(in) :: comp_id
   integer                      :: year, month, day, second

   if (present(annotation)) then
      call get_ccpl_current_time(comp_id, year, month, day, second, 0, trim(annotation)//char(0))
   else
      call get_ccpl_current_time(comp_id, year, month, day, second, 0, trim("")//char(0))
   endif
   if (second .eq. 0 .and. day .eq. 1) then
      CCPL_is_end_current_month = .true.
   else 
      CCPL_is_end_current_month = .false.
   end if

 END FUNCTION CCPL_is_end_current_month



 SUBROUTINE CCPL_get_double_current_calendar_time(comp_id, cal_time, shift_second, annotation)
   implicit none
   real(R8)  cal_time
   character(len=*), intent(in), optional :: annotation
   integer,          intent(in) :: comp_id
   integer,          optional   ::  shift_second
   integer                      ::  local_shift

   local_shift = 0
   if (present(shift_second)) local_shift = shift_second
   if (present(annotation)) then
      call get_ccpl_double_current_calendar_time(comp_id, cal_time, local_shift, trim(annotation)//char(0))
   else
      call get_ccpl_double_current_calendar_time(comp_id, cal_time, local_shift, trim("")//char(0))
   endif
 END SUBROUTINE CCPL_get_double_current_calendar_time



 SUBROUTINE CCPL_get_float_current_calendar_time(comp_id, cal_time, shift_second, annotation)
   implicit none
   real(R4)  cal_time
   character(len=*), intent(in), optional :: annotation
   integer,          intent(in) :: comp_id
   integer,          optional   ::  shift_second
   integer                      ::  local_shift

   local_shift = 0
   if (present(shift_second)) local_shift = shift_second
   if (present(annotation)) then
      call get_ccpl_float_current_calendar_time(comp_id, cal_time, local_shift, trim(annotation)//char(0))
   else
      call get_ccpl_float_current_calendar_time(comp_id, cal_time, local_shift, trim("")//char(0))
   endif

 END SUBROUTINE CCPL_get_float_current_calendar_time



 SUBROUTINE CCPL_allreduce_real16(input_data, output_data, num_data, comm, num_proc)
   implicit none
   real(R16)         :: input_data(:), output_data(:)
   integer           :: num_data, comm
   integer           :: ierr
   integer,optional  :: num_proc
   integer           :: local_num_proc, i, k


#if ( defined  NO_MPI_REAL16 )
   if (present(num_proc)) then
       local_num_proc = num_proc
   else
       call MPI_COMM_SIZE(comm, local_num_proc, ierr)
   end if
   if (reduce_buf_real16_size < num_data*local_num_proc) then
      reduce_buf_real16_size = num_data*local_num_proc*2
      deallocate(reduce_buf_real16)
      allocate(reduce_buf_real16(reduce_buf_real16_size))
   end if
   call mpi_allgather(input_data,num_data*2,MPI_REAL8,reduce_buf_real16,num_data*2,MPI_REAL8,comm,ierr)
   do k=1,num_data
       output_data(k) = 0.0
       do i=1,local_num_proc
          output_data(k) = output_data(k) + reduce_buf_real16((i-1)*num_data+k)
       end do
    enddo
#else
   call mpi_allreduce(input_data,output_data,num_data,MPI_REAL16,MPI_SUM,comm,ierr)
#endif
 END SUBROUTINE CCPL_allreduce_real16
 



   integer FUNCTION CCPL_register_component(parent_id, comp_name, comp_type, comp_comm, considered_in_parent_coupling_gen, change_dir, annotation)
   implicit none
   integer                                 :: rcode
   integer, external                       :: getcwd   ! LINUX system call
   integer, intent(in)                     :: parent_id
   integer                                 :: comp_id
   integer                                 :: ierr
   integer, intent(inout)                  :: comp_comm
   character(len=*), intent(in)            :: comp_type
   character(len=*), intent(in)            :: comp_name
   logical, intent(in), optional           :: considered_in_parent_coupling_gen
   logical, intent(in), optional           :: change_dir
   character(len=*), intent(in), optional  :: annotation
   integer                                 :: local_change_dir
   integer                                 :: local_considered_in_parent_coupling_gen
   character *1024                         :: local_annotation
   character *1024                         :: exe_name


   local_considered_in_parent_coupling_gen = 1
   if (present(considered_in_parent_coupling_gen) .and. (.not. considered_in_parent_coupling_gen)) local_considered_in_parent_coupling_gen = 0

   local_change_dir = 0
   if (present(change_dir) .and. change_dir) local_change_dir = 1

   local_annotation = ""
   if (present(annotation)) then
       local_annotation = annotation
   endif

   if (parent_id .eq. -1) then
      call getarg(0, exe_name)
      call initialize_CCPL_mgrs
      call check_CCPL_Fortran_API_int_type(parent_id)
      call register_root_component(comp_comm, trim(comp_name)//char(0), trim(comp_type)//char(0), trim(local_annotation)//char(0), comp_id, local_considered_in_parent_coupling_gen, local_change_dir, trim(exe_name)//char(0))
   else
      call register_component(parent_id, trim(comp_name)//char(0), trim(comp_type)//char(0), comp_comm, trim(local_annotation)//char(0), local_considered_in_parent_coupling_gen, local_change_dir, comp_id)
   endif
   CCPL_register_component = comp_id

   END FUNCTION CCPL_register_component



   integer FUNCTION CCPL_get_component_id(comp_name, annotation)
   implicit none
   character(len=*), intent(in)            :: comp_name
   character(len=*), intent(in), optional  :: annotation
   character *1024                         :: local_annotation
   integer                     :: comp_id

   local_annotation = ""
   if (present(annotation)) then
       local_annotation = annotation
   endif
 
   call get_id_of_component(trim(comp_name)//char(0), trim(local_annotation)//char(0), comp_id)
   CCPL_get_component_id = comp_id

   end FUNCTION CCPL_get_component_id



   logical FUNCTION CCPL_get_comp_log_file_name(comp_id, file_name, annotation)
   integer, intent(in)                     :: comp_id
   character(len=*), intent(out)           :: file_name
   character(len=*), intent(in), optional  :: annotation
   integer                                 :: log_file_opened 
   character *1024                         :: local_annotation

   local_annotation = ""
   if (present(annotation)) local_annotation = annotation

   call get_ccpl_comp_log_file_name(comp_id, file_name, len(file_name), log_file_opened, trim(local_annotation)//char(0))
   CCPL_get_comp_log_file_name = .false.
   if (log_file_opened .eq. 1) CCPL_get_comp_log_file_name = .true.

   END FUNCTION CCPL_get_comp_log_file_name
   


   integer FUNCTION CCPL_get_comp_log_file_device(comp_id, annotation)
   integer, intent(in)                     :: comp_id
   character(len=*), intent(in), optional  :: annotation
   integer                                 :: log_file_device, log_file_opened 
   character *1024                         :: local_annotation
   character *4096                         :: file_name

   local_annotation = ""
   if (present(annotation)) local_annotation = annotation
   call get_ccpl_comp_log_file_device(comp_id, log_file_device, log_file_opened, trim(local_annotation)//char(0))
   if (log_file_opened .eq. 0) then
       call get_ccpl_comp_log_file_name(comp_id, file_name, len(file_name), log_file_opened, trim("C-Coupler CCPL_get_comp_log_file_device internal")//char(0))
       write(*,*) 'open the file to device ', log_file_device, file_name
       open(unit=log_file_device,file=file_name,position='APPEND')
   endif 
   CCPL_get_comp_log_file_device = log_file_device

   END FUNCTION CCPL_get_comp_log_file_device



   logical FUNCTION CCPL_is_current_process_in_component(comp_full_name, annotation)
   implicit none
   character(len=*), intent(in)            :: comp_full_name
   integer                                 :: is_in_comp
   character(len=*), intent(in), optional  :: annotation
   character *1024                         :: local_annotation

   local_annotation = ""
   if (present(annotation)) local_annotation = annotation

   call is_current_process_in_component(trim(comp_full_name)//char(0), is_in_comp, trim(local_annotation)//char(0))
   CCPL_is_current_process_in_component = (is_in_comp .eq. 1)

   END FUNCTION CCPL_is_current_process_in_component



   integer FUNCTION CCPL_get_current_process_id_in_component(comp_id, annotation)
   implicit none
   integer, intent(in)                     :: comp_id
   character(len=*), intent(in), optional  :: annotation
   character *1024                         :: local_annotation
   integer                                 :: proc_id
   
   local_annotation = ""
   if (present(annotation)) local_annotation = annotation

   call get_current_proc_id_in_comp(comp_id, proc_id, trim(local_annotation)//char(0))
   CCPL_get_current_process_id_in_component = proc_id

   END FUNCTION CCPL_get_current_process_id_in_component



   LOGICAL FUNCTION CCPL_is_comp_type_coupled(comp_id, comp_type, annotation)
   integer, intent(in)                     :: comp_id
   character(len=*), intent(in)            :: comp_type
   character(len=*), intent(in), optional  :: annotation
   integer                                 :: is_coupled
   character *1024                         :: local_annotation

   local_annotation = ""
   if (present(annotation)) local_annotation = annotation

   call check_is_comp_type_coupled(comp_id, trim(comp_type)//char(0), is_coupled, trim(local_annotation)//char(0))

   CCPL_is_comp_type_coupled = .false.
   if (is_coupled .eq. 1) CCPL_is_comp_type_coupled = .true.

   END FUNCTION CCPL_is_comp_type_coupled



   integer FUNCTION CCPL_get_component_process_global_id(comp_id, local_proc_id, annotation)
   implicit none
   integer, intent(in)                     :: comp_id
   integer, intent(in)                     :: local_proc_id
   character(len=*), intent(in), optional  :: annotation
   integer                                 :: global_proc_id
   character *1024                         :: local_annotation
   integer                                 :: proc_id
   
   local_annotation = ""
   if (present(annotation)) local_annotation = annotation

   call get_comp_proc_global_id(comp_id, local_proc_id, global_proc_id, annotation)    
   CCPL_get_component_process_global_id = global_proc_id

   END FUNCTION CCPL_get_component_process_global_id
  


   integer FUNCTION CCPL_get_num_process_in_component(comp_id, annotation)
   implicit none
   integer, intent(in)                     :: comp_id
   character(len=*), intent(in), optional  :: annotation
   character *1024                         :: local_annotation
   integer                                 :: num_proc
   
   local_annotation = ""
   if (present(annotation)) then
       local_annotation = annotation
   endif
   call get_num_proc_in_comp(comp_id, num_proc, annotation)
   CCPL_get_num_process_in_component = num_proc

   END FUNCTION CCPL_get_num_process_in_component



   SUBROUTINE CCPL_end_coupling_configuration(comp_id, annotation)
   implicit none
   integer                                 :: comp_id
   character(len=*), optional              :: annotation
   character *1024                         :: local_annotation
   
   local_annotation = ""
   if (present(annotation)) then
       local_annotation = annotation
   endif

   call ccpl_end_registration(comp_id, trim(local_annotation)//char(0))

   END SUBROUTINE CCPL_end_coupling_configuration



   integer FUNCTION CCPL_register_CoR_defined_grid(comp_id, CCPL_grid_name, CoR_grid_name, annotation) 
   implicit none
   integer, intent(in)                     :: comp_id
   character(len=*), intent(in)            :: CCPL_grid_name
   character(len=*), intent(in)            :: CoR_grid_name
   character(len=*), intent(in), optional  :: annotation
   integer                                 :: grid_id

   if (present(annotation)) then
       call register_cor_defined_grid(comp_id, trim(CCPL_grid_name)//char(0), trim(CoR_grid_name)//char(0), trim(annotation)//char(0), grid_id)
   else 
       call register_cor_defined_grid(comp_id, trim(CCPL_grid_name)//char(0), trim(CoR_grid_name)//char(0), trim("")//char(0), grid_id)
   endif

   CCPL_register_CoR_defined_grid = grid_id

   END FUNCTION CCPL_register_CoR_defined_grid



   integer FUNCTION CCPL_register_H2D_grid_via_file(comp_id, grid_name, data_file_name, annotation)
   implicit none
   integer, intent(in)                     :: comp_id
   character(len=*), intent(in)            :: grid_name
   character(len=*), intent(in)            :: data_file_name
   character(len=*), intent(in), optional  :: annotation
   integer                                 :: grid_id

   if (present(annotation)) then
      call register_H2D_grid_with_file(comp_id, grid_id, trim(grid_name)//char(0), trim(data_file_name)//char(0), trim(annotation)//char(0)) 
   else 
      call register_H2D_grid_with_file(comp_id, grid_id, trim(grid_name)//char(0), trim(data_file_name)//char(0), trim("")//char(0)) 
   endif
   CCPL_register_H2D_grid_via_file = grid_id

   END FUNCTION CCPL_register_H2D_grid_via_file



   integer FUNCTION CCPL_register_H2D_grid_from_another_component(comp_id, grid_name, annotation)
   implicit none
   integer, intent(in)                     :: comp_id
   character(len=*), intent(in)            :: grid_name
   character(len=*), intent(in), optional  :: annotation
   integer                                 :: grid_id

   if (present(annotation)) then
      call register_H2D_grid_from_another_component(comp_id, grid_id, trim(grid_name)//char(0), trim(annotation)//char(0)) 
   else 
      call register_H2D_grid_from_another_component(comp_id, grid_id, trim(grid_name)//char(0), trim("")//char(0)) 
   endif
   CCPL_register_H2D_grid_from_another_component = grid_id

   END FUNCTION CCPL_register_H2D_grid_from_another_component



   integer FUNCTION CCPL_register_H2D_grid_global_online_C1D_M1D_float(comp_id, grid_name, edge_type, coord_unit, cyclic_or_acyclic, dim_size1, dim_size2, &
                                                           min_lon, max_lon, min_lat, max_lat, center_lon, center_lat, mask, area, vertex_lon, vertex_lat, annotation)
   implicit none
   integer, intent(in)                                              :: comp_id
   character(len=*), intent(in)                                     :: grid_name, edge_type, coord_unit, cyclic_or_acyclic
   integer,          intent(in)                                     :: dim_size1, dim_size2
   real(R4),         intent(in)                                     :: min_lon, max_lon, min_lat, max_lat
   real(R4),         intent(in), dimension(:)                       :: center_lon, center_lat
   integer,          intent(in), dimension(:),   target, optional   :: mask
   real(R4),         intent(in), dimension(:),   target, optional   :: area
   real(R4),         intent(in), dimension(:,:), target, optional   :: vertex_lon, vertex_lat
   character(len=*), intent(in),                         optional   :: annotation
   integer,                      dimension(:),   pointer            :: temp_mask, temp_int
   real(R4),                     dimension(:),   pointer            :: temp_area, temp_float_1d
   real(R4),                     dimension(:,:), pointer            :: temp_vertex_lon, temp_vertex_lat, temp_float_2d
   integer                                                          :: grid_id
   integer                                                          :: size_center_lon, size_center_lat, size_mask, size_area, size_vertex_lon, size_vertex_lat
   
   size_center_lon = size(center_lon)
   size_center_lat = size(center_lat)
   size_mask = -1
   size_area = -1
   size_vertex_lon = -1
   size_vertex_lat = -1
   allocate(temp_int(1), temp_float_1d(1), temp_float_2d(1,1))
   temp_mask => temp_int
   temp_area => temp_float_1d
   temp_vertex_lon => temp_float_2d
   temp_vertex_lat => temp_float_2d
   if (present(mask)) then
      size_mask = size(mask)
      temp_mask => mask
   endif
   if (present(area)) then
      size_area = size(area)
      temp_area => area
   endif
   if (present(vertex_lon)) then
      size_vertex_lon = size(vertex_lon)
      temp_vertex_lon => vertex_lon
   endif
   if (present(vertex_lat)) then
      size_vertex_lat = size(vertex_lat)
      temp_vertex_lat => vertex_lat
   endif

   if (present(annotation)) then
      call register_H2D_grid_with_global_data(comp_id, grid_id, trim(grid_name)//char(0), trim(edge_type)//char(0), trim(coord_unit)//char(0), trim(cyclic_or_acyclic)//char(0), &
                                              trim("real4")//char(0), dim_size1, dim_size2, size_center_lon, size_center_lat, size_mask, size_area, size_vertex_lon, size_vertex_lat,  &
                                              min_lon, max_lon, min_lat, max_lat, center_lon, center_lat, temp_mask, temp_area, temp_vertex_lon, temp_vertex_lat, trim(annotation)//char(0))
   else
      call register_H2D_grid_with_global_data(comp_id, grid_id, trim(grid_name)//char(0), trim(edge_type)//char(0), trim(coord_unit)//char(0), trim(cyclic_or_acyclic)//char(0), &
                                              trim("real4")//char(0), dim_size1, dim_size2, size_center_lon, size_center_lat, size_mask, size_area, size_vertex_lon, size_vertex_lat,  &
                                              min_lon, max_lon, min_lat, max_lat, center_lon, center_lat, temp_mask, temp_area, temp_vertex_lon, temp_vertex_lat, trim("")//char(0))
   endif

   CCPL_register_H2D_grid_global_online_C1D_M1D_float = grid_id

   deallocate(temp_int, temp_float_1d, temp_float_2d)

   end FUNCTION CCPL_register_H2D_grid_global_online_C1D_M1D_float



   integer FUNCTION CCPL_register_H2D_grid_global_online_C1D_M1D_double(comp_id, grid_name, edge_type, coord_unit, cyclic_or_acyclic, dim_size1, dim_size2, &
                                                           min_lon, max_lon, min_lat, max_lat, center_lon, center_lat, mask, area, vertex_lon, vertex_lat, annotation)
   implicit none
   integer, intent(in)                                              :: comp_id
   character(len=*), intent(in)                                     :: grid_name, edge_type, coord_unit, cyclic_or_acyclic
   integer,          intent(in)                                     :: dim_size1, dim_size2
   real(R8),         intent(in)                                     :: min_lon, max_lon, min_lat, max_lat
   real(R8),         intent(in), dimension(:)                       :: center_lon, center_lat
   integer,          intent(in), dimension(:),   target, optional   :: mask
   real(R8),         intent(in), dimension(:),   target, optional   :: area
   real(R8),         intent(in), dimension(:,:), target, optional   :: vertex_lon, vertex_lat
   character(len=*), intent(in),                         optional   :: annotation
   integer,                      dimension(:),   pointer            :: temp_mask, temp_int
   real(R8),                     dimension(:),   pointer            :: temp_area, temp_double_1d
   real(R8),                     dimension(:,:), pointer            :: temp_vertex_lon, temp_vertex_lat, temp_double_2d
   integer                                                          :: grid_id
   integer                                                          :: size_center_lon, size_center_lat, size_mask, size_area, size_vertex_lon, size_vertex_lat
  

   size_center_lon = size(center_lon)
   size_center_lat = size(center_lat)
   size_mask = -1
   size_area = -1
   size_vertex_lon = -1
   size_vertex_lat = -1
   allocate(temp_int(1), temp_double_1d(1), temp_double_2d(1,1))
   temp_mask => temp_int
   temp_area => temp_double_1d
   temp_vertex_lon => temp_double_2d
   temp_vertex_lat => temp_double_2d
   if (present(mask)) then
      size_mask = size(mask)
      temp_mask => mask
   endif
   if (present(area)) then
      size_area = size(area)
      temp_area => area
   endif
   if (present(vertex_lon)) then
      size_vertex_lon = size(vertex_lon)
      temp_vertex_lon => vertex_lon
   endif
   if (present(vertex_lat)) then
      size_vertex_lat = size(vertex_lat)
      temp_vertex_lat => vertex_lat
   endif

   if (present(annotation)) then
      call register_H2D_grid_with_global_data(comp_id, grid_id, trim(grid_name)//char(0), trim(edge_type)//char(0), trim(coord_unit)//char(0), trim(cyclic_or_acyclic)//char(0), &
                                        trim("real8")//char(0), dim_size1, dim_size2, size_center_lon, size_center_lat, size_mask, size_area, size_vertex_lon, size_vertex_lat,  &
                                       min_lon, max_lon, min_lat, max_lat, center_lon, center_lat, temp_mask, temp_area, temp_vertex_lon, temp_vertex_lat, trim(annotation)//char(0))
   else
      call register_H2D_grid_with_global_data(comp_id, grid_id, trim(grid_name)//char(0), trim(edge_type)//char(0), trim(coord_unit)//char(0), trim(cyclic_or_acyclic)//char(0), &
                                       trim("real8")//char(0), dim_size1, dim_size2, size_center_lon, size_center_lat, size_mask, size_area, size_vertex_lon, size_vertex_lat,  &
                                       min_lon, max_lon, min_lat, max_lat, center_lon, center_lat, temp_mask, temp_area, temp_vertex_lon, temp_vertex_lat, trim("")//char(0))
   endif

   CCPL_register_H2D_grid_global_online_C1D_M1D_double = grid_id

   deallocate(temp_int, temp_double_1d, temp_double_2d)

   end FUNCTION CCPL_register_H2D_grid_global_online_C1D_M1D_double



   integer FUNCTION CCPL_register_H2D_grid_global_online_C2D_M2D_float(comp_id, grid_name, edge_type, coord_unit, cyclic_or_acyclic, dim_size1, dim_size2, &
                                                           min_lon, max_lon, min_lat, max_lat, center_lon, center_lat, mask, area, vertex_lon, vertex_lat, annotation)
   implicit none
   integer, intent(in)                                              :: comp_id
   character(len=*), intent(in)                                     :: grid_name, edge_type, coord_unit, cyclic_or_acyclic
   integer,          intent(in)                                     :: dim_size1, dim_size2
   real(R4),         intent(in)                                     :: min_lon, max_lon, min_lat, max_lat
   real(R4),         intent(in), dimension(:,:)                     :: center_lon, center_lat
   integer,          intent(in), dimension(:,:),   target, optional :: mask
   real(R4),         intent(in), dimension(:,:),   target, optional :: area
   real(R4),         intent(in), dimension(:,:,:), target, optional :: vertex_lon, vertex_lat
   character(len=*), intent(in),                           optional :: annotation
   integer,                      dimension(:,:),   pointer          :: temp_mask, temp_int
   real(R4),                     dimension(:,:),   pointer          :: temp_area, temp_float_2d
   real(R4),                     dimension(:,:,:), pointer          :: temp_vertex_lon, temp_vertex_lat, temp_float_3d
   integer                                                          :: grid_id
   integer                                                          :: size_center_lon, size_center_lat, size_mask, size_area, size_vertex_lon, size_vertex_lat
   
   size_center_lon = size(center_lon)
   size_center_lat = size(center_lat)
   size_mask = -1
   size_area = -1
   size_vertex_lon = -1
   size_vertex_lat = -1
   allocate(temp_int(1,1), temp_float_2d(1,1), temp_float_3d(1,1,1))
   temp_mask => temp_int
   temp_area => temp_float_2d
   temp_vertex_lon => temp_float_3d
   temp_vertex_lat => temp_float_3d
   if (present(mask)) then
      size_mask = size(mask)
      temp_mask => mask
   endif
   if (present(area)) then
      size_area = size(area)
      temp_area => area
   endif
   if (present(vertex_lon)) then
      size_vertex_lon = size(vertex_lon)
      temp_vertex_lon => vertex_lon
   endif
   if (present(vertex_lat)) then
      size_vertex_lat = size(vertex_lat)
      temp_vertex_lat => vertex_lat
   endif

   if (present(annotation)) then
      call register_H2D_grid_with_global_data(comp_id, grid_id, trim(grid_name)//char(0), trim(edge_type)//char(0), trim(coord_unit)//char(0), trim(cyclic_or_acyclic)//char(0), &
                                       trim("real4")//char(0), dim_size1, dim_size2, size_center_lon, size_center_lat, size_mask, size_area, size_vertex_lon, size_vertex_lat,  &
                                       min_lon, max_lon, min_lat, max_lat, center_lon, center_lat, temp_mask, temp_area, temp_vertex_lon, temp_vertex_lat, trim(annotation)//char(0))
   else
      call register_H2D_grid_with_global_data(comp_id, grid_id, trim(grid_name)//char(0), trim(edge_type)//char(0), trim(coord_unit)//char(0), trim(cyclic_or_acyclic)//char(0), &
                                       trim("real4")//char(0), dim_size1, dim_size2, size_center_lon, size_center_lat, size_mask, size_area, size_vertex_lon, size_vertex_lat,  &
                                       min_lon, max_lon, min_lat, max_lat, center_lon, center_lat, temp_mask, temp_area, temp_vertex_lon, temp_vertex_lat, trim("")//char(0))
   endif

   CCPL_register_H2D_grid_global_online_C2D_M2D_float = grid_id

   deallocate(temp_int, temp_float_2d, temp_float_3d)

   end FUNCTION CCPL_register_H2D_grid_global_online_C2D_M2D_float



   integer FUNCTION CCPL_register_H2D_grid_global_online_C2D_M2D_double(comp_id, grid_name, edge_type, coord_unit, cyclic_or_acyclic, dim_size1, dim_size2, &
                                                           min_lon, max_lon, min_lat, max_lat, center_lon, center_lat, mask, area, vertex_lon, vertex_lat, annotation)
   implicit none
   integer, intent(in)                                              :: comp_id
   character(len=*), intent(in)                                     :: grid_name, edge_type, coord_unit, cyclic_or_acyclic
   integer,          intent(in)                                     :: dim_size1, dim_size2
   real(R8),         intent(in)                                     :: min_lon, max_lon, min_lat, max_lat
   real(R8),         intent(in), dimension(:,:)                     :: center_lon, center_lat
   integer,          intent(in), dimension(:,:),   target, optional :: mask
   real(R8),         intent(in), dimension(:,:),   target, optional :: area
   real(R8),         intent(in), dimension(:,:,:), target, optional :: vertex_lon, vertex_lat
   character(len=*), intent(in),                           optional :: annotation
   integer,                      dimension(:,:),   pointer          :: temp_mask, temp_int
   real(R8),                     dimension(:,:),   pointer          :: temp_area, temp_double_2d
   real(R8),                     dimension(:,:,:), pointer          :: temp_vertex_lon, temp_vertex_lat, temp_double_3d
   integer                                                          :: grid_id
   integer                                                          :: size_center_lon, size_center_lat, size_mask, size_area, size_vertex_lon, size_vertex_lat
  

   size_center_lon = size(center_lon)
   size_center_lat = size(center_lat)
   size_mask = -1
   size_area = -1
   size_vertex_lon = -1
   size_vertex_lat = -1
   allocate(temp_int(1,1), temp_double_2d(1,1), temp_double_3d(1,1,1))
   temp_mask => temp_int
   temp_area => temp_double_2d
   temp_vertex_lon => temp_double_3d
   temp_vertex_lat => temp_double_3d
   if (present(mask)) then
      size_mask = size(mask)
      temp_mask => mask
   endif
   if (present(area)) then
      size_area = size(area)
      temp_area => area
   endif
   if (present(vertex_lon)) then
      size_vertex_lon = size(vertex_lon)
      temp_vertex_lon => vertex_lon
   endif
   if (present(vertex_lat)) then
      size_vertex_lat = size(vertex_lat)
      temp_vertex_lat => vertex_lat
   endif

   if (present(annotation)) then
      call register_H2D_grid_with_global_data(comp_id, grid_id, trim(grid_name)//char(0), trim(edge_type)//char(0), trim(coord_unit)//char(0), trim(cyclic_or_acyclic)//char(0), &
                                        trim("real8")//char(0), dim_size1, dim_size2, size_center_lon, size_center_lat, size_mask, size_area, size_vertex_lon, size_vertex_lat,  &
                                       min_lon, max_lon, min_lat, max_lat, center_lon, center_lat, temp_mask, temp_area, temp_vertex_lon, temp_vertex_lat, trim(annotation)//char(0))
   else
      call register_H2D_grid_with_global_data(comp_id, grid_id, trim(grid_name)//char(0), trim(edge_type)//char(0), trim(coord_unit)//char(0), trim(cyclic_or_acyclic)//char(0), &
                                       trim("real8")//char(0), dim_size1, dim_size2, size_center_lon, size_center_lat, size_mask, size_area, size_vertex_lon, size_vertex_lat,  &
                                       min_lon, max_lon, min_lat, max_lat, center_lon, center_lat, temp_mask, temp_area, temp_vertex_lon, temp_vertex_lat, trim("")//char(0))
   endif

   CCPL_register_H2D_grid_global_online_C2D_M2D_double = grid_id

   deallocate(temp_int, temp_double_2d, temp_double_3d)

   end FUNCTION CCPL_register_H2D_grid_global_online_C2D_M2D_double




   integer FUNCTION CCPL_register_H2D_grid_local_online_float(comp_id, grid_name, edge_type, coord_unit, cyclic_or_acyclic, grid_size, num_local_cells, local_cells_global_index, &
                                                           min_lon, max_lon, min_lat, max_lat, center_lon, center_lat, mask, area, vertex_lon, vertex_lat, decomp_name, decomp_id, annotation)
   implicit none
   integer, intent(in)                                              :: comp_id
   character(len=*), intent(in)                                     :: grid_name, edge_type, coord_unit, cyclic_or_acyclic
   integer,          intent(in)                                     :: grid_size
   integer,          intent(in)                                     :: num_local_cells
   integer,          intent(in), dimension(:)                       :: local_cells_global_index(:)
   real(R4),         intent(in)                                     :: min_lon, max_lon, min_lat, max_lat 
   real(R4),         intent(in), dimension(:)                       :: center_lon, center_lat
   integer,          intent(in), dimension(:),   target, optional   :: mask
   real(R4),         intent(in), dimension(:),   target, optional   :: area
   real(R4),         intent(in), dimension(:,:), target, optional   :: vertex_lon, vertex_lat
   character(len=*), intent(in),                         optional   :: decomp_name
   integer,          intent(out),                        optional   :: decomp_id
   character(len=*), intent(in),                         optional   :: annotation
   integer,                      dimension(:),   pointer            :: temp_mask, temp_int
   real(R4),                     dimension(:),   pointer            :: temp_area, temp_float_1d
   real(R4),                     dimension(:,:), pointer            :: temp_vertex_lon, temp_vertex_lat, temp_float_2d
   integer                                                          :: grid_id, local_decomp_id
   integer                                                          :: size_center_lon, size_center_lat, size_mask, size_area, size_vertex_lon, size_vertex_lat
   character *2048                                                   :: local_decomp_name, local_annotation
   
   size_center_lon = size(center_lon)
   size_center_lat = size(center_lat)
   size_mask = -1
   size_area = -1
   size_vertex_lon = -1
   size_vertex_lat = -1
   allocate(temp_int(1), temp_float_1d(1), temp_float_2d(1,1))
   temp_mask => temp_int
   temp_area => temp_float_1d
   temp_vertex_lon => temp_float_2d
   temp_vertex_lat => temp_float_2d

   local_decomp_name = ""
   if (present(decomp_name)) local_decomp_name = decomp_name
   local_annotation = ""
   if (present(annotation)) local_annotation = annotation
   local_decomp_id = -1
   if (present(decomp_id)) then
       local_decomp_id = 0
   endif

   if (present(mask)) then
      size_mask = size(mask)
      temp_mask => mask
   endif
   if (present(area)) then
      size_area = size(area)
      temp_area => area
   endif
   if (present(vertex_lon)) then
      size_vertex_lon = size(vertex_lon)
      temp_vertex_lon => vertex_lon
   endif
   if (present(vertex_lat)) then
      size_vertex_lat = size(vertex_lat)
      temp_vertex_lat => vertex_lat
   endif

   call register_H2D_grid_with_local_data(comp_id, grid_id, trim(grid_name)//char(0), trim(edge_type)//char(0), trim(coord_unit)//char(0), trim(cyclic_or_acyclic)//char(0), &
                                          trim("real4")//char(0), grid_size, num_local_cells, size(local_cells_global_index), size_center_lon, size_center_lat, size_mask, size_area, &
                                          size_vertex_lon, size_vertex_lat, local_cells_global_index, min_lon, max_lon, min_lat, max_lat, center_lon, center_lat, temp_mask, temp_area, &
                                          temp_vertex_lon, temp_vertex_lat, trim(local_decomp_name)//char(0), local_decomp_id, trim(local_annotation)//char(0))

   CCPL_register_H2D_grid_local_online_float = grid_id

   if (present(decomp_id)) decomp_id = local_decomp_id

   deallocate(temp_int, temp_float_1d, temp_float_2d)

   end FUNCTION CCPL_register_H2D_grid_local_online_float



   integer FUNCTION CCPL_register_H2D_grid_local_online_double(comp_id, grid_name, edge_type, coord_unit, cyclic_or_acyclic, grid_size, num_local_cells, local_cells_global_index, &
                                                           min_lon, max_lon, min_lat, max_lat, center_lon, center_lat, mask, area, vertex_lon, vertex_lat, decomp_name, decomp_id, annotation)
   implicit none
   integer, intent(in)                                              :: comp_id
   character(len=*), intent(in)                                     :: grid_name, edge_type, coord_unit, cyclic_or_acyclic
   integer,          intent(in)                                     :: grid_size
   integer,          intent(in)                                     :: num_local_cells
   integer,          intent(in), dimension(:)                       :: local_cells_global_index(:)
   real(R8),         intent(in)                                     :: min_lon, max_lon, min_lat, max_lat 
   real(R8),         intent(in), dimension(:)                       :: center_lon, center_lat
   integer,          intent(in), dimension(:),   target, optional   :: mask
   real(R8),         intent(in), dimension(:),   target, optional   :: area
   real(R8),         intent(in), dimension(:,:), target, optional   :: vertex_lon, vertex_lat
   character(len=*), intent(in),                         optional   :: decomp_name
   integer,          intent(out),                        optional   :: decomp_id
   character(len=*), intent(in),                         optional   :: annotation
   integer,                      dimension(:),   pointer            :: temp_mask, temp_int
   real(R8),                     dimension(:),   pointer            :: temp_area, temp_double_1d
   real(R8),                     dimension(:,:), pointer            :: temp_vertex_lon, temp_vertex_lat, temp_double_2d
   integer                                                          :: grid_id, local_decomp_id
   integer                                                          :: size_center_lon, size_center_lat, size_mask, size_area, size_vertex_lon, size_vertex_lat
   character *2048                                                   :: local_decomp_name, local_annotation
   
   size_center_lon = size(center_lon)
   size_center_lat = size(center_lat)
   size_mask = -1
   size_area = -1
   size_vertex_lon = -1
   size_vertex_lat = -1
   allocate(temp_int(1), temp_double_1d(1), temp_double_2d(1,1))
   temp_mask => temp_int
   temp_area => temp_double_1d
   temp_vertex_lon => temp_double_2d
   temp_vertex_lat => temp_double_2d

   local_decomp_name = ""
   if (present(decomp_name)) local_decomp_name = decomp_name
   local_annotation = ""
   if (present(annotation)) local_annotation = annotation
   local_decomp_id = -1
   if (present(decomp_id)) then
       local_decomp_id = 0
   endif

   if (present(mask)) then
      size_mask = size(mask)
      temp_mask => mask
   endif
   if (present(area)) then
      size_area = size(area)
      temp_area => area
   endif
   if (present(vertex_lon)) then
      size_vertex_lon = size(vertex_lon)
      temp_vertex_lon => vertex_lon
   endif
   if (present(vertex_lat)) then
      size_vertex_lat = size(vertex_lat)
      temp_vertex_lat => vertex_lat
   endif

   call register_H2D_grid_with_local_data(comp_id, grid_id, trim(grid_name)//char(0), trim(edge_type)//char(0), trim(coord_unit)//char(0), trim(cyclic_or_acyclic)//char(0), &
                                          trim("real8")//char(0), grid_size, num_local_cells, size(local_cells_global_index), size_center_lon, size_center_lat, size_mask, size_area, &
                                          size_vertex_lon, size_vertex_lat, local_cells_global_index, min_lon, max_lon, min_lat, max_lat, center_lon, center_lat, temp_mask, temp_area, &
                                          temp_vertex_lon, temp_vertex_lat, trim(local_decomp_name)//char(0), local_decomp_id, trim(local_annotation)//char(0))

   CCPL_register_H2D_grid_local_online_double = grid_id

   if (present(decomp_id)) decomp_id = local_decomp_id

   deallocate(temp_int, temp_double_1d, temp_double_2d)

   end FUNCTION CCPL_register_H2D_grid_local_online_double



   integer FUNCTION CCPL_register_V1D_grid_without_data(comp_id, grid_name, coord_unit, grid_size, annotation)
   implicit none
   integer, intent(in)                                     :: comp_id
   integer, intent(in)                                     :: grid_size
   character(len=*), intent(in)                            :: grid_name
   character(len=*), intent(in)                            :: coord_unit
   character(len=*), intent(in),               optional    :: annotation
   integer                                                 :: grid_id

   if (present(annotation)) then
       call register_V1D_grid_without_data(comp_id, grid_id, trim(grid_name)//char(0), trim(coord_unit)//char(0), grid_size, trim(annotation)//char(0))
   else
       call register_V1D_grid_without_data(comp_id, grid_id, trim(grid_name)//char(0), trim(coord_unit)//char(0), grid_size, trim("")//char(0))
   endif

   CCPL_register_V1D_grid_without_data = grid_id

   END FUNCTION CCPL_register_V1D_grid_without_data



   integer FUNCTION CCPL_register_V1D_Z_grid_via_double_data(comp_id, grid_name, coord_unit, coord_values, annotation)
   implicit none
   integer, intent(in)                                     :: comp_id
   character(len=*), intent(in)                            :: grid_name
   character(len=*), intent(in)                            :: coord_unit
   real(R8),         intent(in), dimension(:)              :: coord_values
   character(len=*), intent(in),               optional    :: annotation
   integer                                                 :: grid_id

   if (present(annotation)) then
       call register_V1D_grid_with_data(comp_id, grid_id, trim(grid_name)//char(0), 1, trim(coord_unit)//char(0), size(coord_values), size(coord_values), trim("real8")//char(0), coord_values, coord_values, coord_values, trim(annotation)//char(0))
   else
       call register_V1D_grid_with_data(comp_id, grid_id, trim(grid_name)//char(0), 1, trim(coord_unit)//char(0), size(coord_values), size(coord_values), trim("real8")//char(0), coord_values, coord_values, coord_values, trim("")//char(0))
   endif

   CCPL_register_V1D_Z_grid_via_double_data = grid_id

   END FUNCTION CCPL_register_V1D_Z_grid_via_double_data



   integer FUNCTION CCPL_register_V1D_Z_grid_via_float_data(comp_id, grid_name, coord_unit, coord_values, annotation)
   implicit none
   integer, intent(in)                                     :: comp_id
   character(len=*), intent(in)                            :: grid_name
   character(len=*), intent(in)                            :: coord_unit
   real(R4),         intent(in), dimension(:)              :: coord_values
   character(len=*), intent(in),               optional    :: annotation
   integer                                                 :: grid_id

   if (present(annotation)) then
       call register_V1D_grid_with_data(comp_id, grid_id, trim(grid_name)//char(0), 1, trim(coord_unit)//char(0), size(coord_values), size(coord_values), trim("real4")//char(0), coord_values, coord_values, coord_values, trim(annotation)//char(0))
   else
       call register_V1D_grid_with_data(comp_id, grid_id, trim(grid_name)//char(0), 1, trim(coord_unit)//char(0), size(coord_values), size(coord_values), trim("real4")//char(0), coord_values, coord_values, coord_values, trim("")//char(0))
   endif

   CCPL_register_V1D_Z_grid_via_float_data = grid_id

   END FUNCTION CCPL_register_V1D_Z_grid_via_float_data



   integer FUNCTION CCPL_register_V1D_SIGMA_grid_via_double_data(comp_id, grid_name, coord_unit, P0, sigma_values, annotation)
   implicit none
   integer, intent(in)                                     :: comp_id
   character(len=*), intent(in)                            :: grid_name
   character(len=*), intent(in)                            :: coord_unit
   real(R8),         intent(in)                            :: P0
   real(R8),         intent(in), dimension(:)              :: sigma_values
   character(len=*), intent(in),               optional    :: annotation
   integer                                                 :: grid_id

   if (present(annotation)) then
       call register_V1D_grid_with_data(comp_id, grid_id, trim(grid_name)//char(0), 2, trim(coord_unit)//char(0), size(sigma_values), size(sigma_values), trim("real8")//char(0), P0, sigma_values, sigma_values, trim(annotation)//char(0))
   else
       call register_V1D_grid_with_data(comp_id, grid_id, trim(grid_name)//char(0), 2, trim(coord_unit)//char(0), size(sigma_values), size(sigma_values), trim("real8")//char(0), P0, sigma_values, sigma_values, trim("")//char(0))
   endif

   CCPL_register_V1D_SIGMA_grid_via_double_data = grid_id

   END FUNCTION CCPL_register_V1D_SIGMA_grid_via_double_data



   integer FUNCTION CCPL_register_V1D_SIGMA_grid_via_float_data(comp_id, grid_name, coord_unit, P0, sigma_values, annotation)
   implicit none
   integer, intent(in)                                     :: comp_id
   character(len=*), intent(in)                            :: grid_name
   character(len=*), intent(in)                            :: coord_unit
   real(R4),         intent(in)                            :: P0
   real(R4),         intent(in), dimension(:)              :: sigma_values
   character(len=*), intent(in),               optional    :: annotation
   integer                                                 :: grid_id

   if (present(annotation)) then
       call register_V1D_grid_with_data(comp_id, grid_id, trim(grid_name)//char(0), 2, trim(coord_unit)//char(0), size(sigma_values), size(sigma_values), trim("real4")//char(0), P0, sigma_values, sigma_values, trim(annotation)//char(0))
   else
       call register_V1D_grid_with_data(comp_id, grid_id, trim(grid_name)//char(0), 2, trim(coord_unit)//char(0), size(sigma_values), size(sigma_values), trim("real4")//char(0), P0, sigma_values, sigma_values, trim("")//char(0))
   endif

   CCPL_register_V1D_SIGMA_grid_via_float_data = grid_id

   END FUNCTION CCPL_register_V1D_SIGMA_grid_via_float_data



   integer FUNCTION CCPL_register_V1D_HYBRID_grid_via_double_data(comp_id, grid_name, coord_unit, P0, coef_A, coef_B, annotation)
   implicit none
   integer, intent(in)                                     :: comp_id
   character(len=*), intent(in)                            :: grid_name
   character(len=*), intent(in)                            :: coord_unit
   real(R8),         intent(in), dimension(:)              :: coef_A
   real(R8),         intent(in), dimension(:)              :: coef_B
   character(len=*), intent(in),               optional    :: annotation
   real(R8),         intent(in)                            :: P0
   integer                                                 :: grid_id

   if (present(annotation)) then
       call register_V1D_grid_with_data(comp_id, grid_id, trim(grid_name)//char(0), 3, trim(coord_unit)//char(0), size(coef_B), size(coef_A), trim("real8")//char(0), P0, coef_B, coef_A, trim(annotation)//char(0))
   else
       call register_V1D_grid_with_data(comp_id, grid_id, trim(grid_name)//char(0), 3, trim(coord_unit)//char(0), size(coef_B), size(coef_A), trim("real8")//char(0), P0, coef_B, coef_A, trim("")//char(0))
   endif

   CCPL_register_V1D_HYBRID_grid_via_double_data = grid_id

   END FUNCTION CCPL_register_V1D_HYBRID_grid_via_double_data



   integer FUNCTION CCPL_register_V1D_HYBRID_grid_via_float_data(comp_id, grid_name, coord_unit, P0, coef_A, coef_B, annotation)
   implicit none
   integer, intent(in)                                     :: comp_id
   character(len=*), intent(in)                            :: grid_name
   character(len=*), intent(in)                            :: coord_unit
   real(R4),         intent(in), dimension(:)              :: coef_A
   real(R4),         intent(in), dimension(:)              :: coef_B
   character(len=*), intent(in),               optional    :: annotation
   real(R4),         intent(in)                            :: P0
   integer                                                 :: grid_id

   if (present(annotation)) then
       call register_V1D_grid_with_data(comp_id, grid_id, trim(grid_name)//char(0), 3, trim(coord_unit)//char(0), size(coef_B), size(coef_A), trim("real4")//char(0), P0, coef_B, coef_A, trim(annotation)//char(0))
   else
       call register_V1D_grid_with_data(comp_id, grid_id, trim(grid_name)//char(0), 3, trim(coord_unit)//char(0), size(coef_B), size(coef_A), trim("real4")//char(0), P0, coef_B, coef_A, trim("")//char(0))
   endif

   CCPL_register_V1D_HYBRID_grid_via_float_data = grid_id

   END FUNCTION CCPL_register_V1D_HYBRID_grid_via_float_data



   integer FUNCTION CCPL_register_MD_grid_via_multi_grids(comp_id, grid_name, sub_grid1_id, sub_grid2_id, sub_grid3_id, mask, annotation)
   implicit none
   integer, intent(in)                                               :: comp_id
   character(len=*), intent(in)                                      :: grid_name
   integer, intent(in)                                               :: sub_grid1_id
   integer, intent(in)                                               :: sub_grid2_id
   integer, intent(in),                                optional      :: sub_grid3_id
   character(len=*), intent(in),                       optional      :: annotation
   integer,          intent(in), dimension(:), target, optional      :: mask
   integer,                      dimension(:), pointer               :: temp_mask, temp_int
   integer                                                           :: grid_id, local_sub_grid3_id, size_mask
   
   
   size_mask = -1
   allocate(temp_int(1))
   temp_mask => temp_int
   if (present(mask)) then
       size_mask = size(mask)
       temp_mask => mask
   endif
   local_sub_grid3_id = -1
   if (present(sub_grid3_id)) local_sub_grid3_id = sub_grid3_id 

   if (present(annotation)) then
      call register_MD_grid_via_multi_grids(comp_id, grid_id, trim(grid_name)//char(0), sub_grid1_id, sub_grid2_id, local_sub_grid3_id, size_mask, temp_mask, trim(annotation)//char(0))
   else 
      call register_MD_grid_via_multi_grids(comp_id, grid_id, trim(grid_name)//char(0), sub_grid1_id, sub_grid2_id, local_sub_grid3_id, size_mask, temp_mask, trim("")//char(0))
   endif

   CCPL_register_MD_grid_via_multi_grids = grid_id

   deallocate(temp_int)

   END FUNCTION CCPL_register_MD_grid_via_multi_grids 



   SUBROUTINE CCPL_set_3D_grid_3D_vertical_coord_field(grid_id, field_id, label, annotation)
   implicit none
   integer, intent(in)                                     :: grid_id
   integer, intent(in)                                     :: field_id
   character(len=*), intent(in)                            :: label
   character(len=*), intent(in),               optional    :: annotation

   if (present(annotation)) then
       call set_3D_grid_3D_vertical_coord_field(grid_id, field_id, trim(label)//char(0), trim(annotation)//char(0))
   else
       call set_3D_grid_3D_vertical_coord_field(grid_id, field_id, trim(label)//char(0), trim("")//char(0))
   endif

   END SUBROUTINE CCPL_set_3D_grid_3D_vertical_coord_field



   SUBROUTINE CCPL_set_3D_grid_variable_surface_field(grid_id, field_id, annotation)
   implicit none
   integer, intent(in)                                     :: grid_id
   integer, intent(in)                                     :: field_id
   character(len=*), intent(in),               optional    :: annotation
 
   
   if (present(annotation)) then
       call set_3D_grid_surface_field(grid_id, field_id, 1, trim(annotation)//char(0))
   else
       call set_3D_grid_surface_field(grid_id, field_id, 1, trim("")//char(0))
   endif

   END SUBROUTINE CCPL_set_3D_grid_variable_surface_field
   


   SUBROUTINE CCPL_set_3D_grid_constant_surface_field(grid_id, field_id, annotation)
   implicit none
   integer, intent(in)                                     :: grid_id
   integer, intent(in)                                     :: field_id
   character(len=*), intent(in),               optional    :: annotation
 
   
   if (present(annotation)) then
       call set_3D_grid_surface_field(grid_id, field_id, 0, trim(annotation)//char(0))
   else
       call set_3D_grid_surface_field(grid_id, field_id, 0, trim("")//char(0))
   endif

   END SUBROUTINE CCPL_set_3D_grid_constant_surface_field



   SUBROUTINE CCPL_set_3D_grid_external_surface_field(grid_id, annotation)
   implicit none
   integer, intent(in)                                     :: grid_id
   character(len=*), intent(in),               optional    :: annotation


   if (present(annotation)) then
       call set_3D_grid_surface_field(grid_id, -1, 2, trim(annotation)//char(0))
   else
       call set_3D_grid_surface_field(grid_id, -1, 2, trim("")//char(0))
   endif

   END SUBROUTINE CCPL_set_3D_grid_external_surface_field



   SUBROUTINE CCPL_register_mid_point_grid(level_3D_grid_id, mid_3D_grid_id, mid_V1D_grid_id, mask, annotation)
   implicit none
   integer, intent(in)                                     :: level_3D_grid_id
   integer, intent(out)                                    :: mid_3D_grid_id
   integer, intent(out)                                    :: mid_V1D_grid_id
   integer, intent(in),  dimension(:), target, optional    :: mask
   integer,              dimension(:),         pointer     :: temp_mask, temp_int
   character(len=*), intent(in),               optional    :: annotation
   integer                                                 :: size_mask


   allocate(temp_int(1))
   size_mask = -1
   temp_mask => temp_int
   if (present(mask)) then
      temp_mask=>mask
      size_mask = size(mask)
   endif
   if (present(annotation)) then
      call register_mid_point_grid(level_3D_grid_id, mid_3D_grid_id, mid_V1D_grid_id, size_mask, temp_mask, trim(annotation)//char(0))
   else
      call register_mid_point_grid(level_3D_grid_id, mid_3D_grid_id, mid_V1D_grid_id, size_mask, temp_mask, trim("")//char(0))
   endif
   deallocate(temp_int)
   END SUBROUTINE CCPL_register_mid_point_grid



   integer FUNCTION CCPL_get_grid_size(grid_id, annotation) 
   implicit none
   integer, intent(in)                     :: grid_id
   integer                                 :: grid_size
   character(len=*), intent(in), optional  :: annotation

   if (present(annotation)) then
       call get_grid_size(grid_id, grid_size, trim(annotation)//char(0))
   else 
       call get_grid_size(grid_id, grid_size, trim("")//char(0))
   endif
   
   CCPL_get_grid_size = grid_size
    
   END FUNCTION CCPL_get_grid_size



   integer FUNCTION CCPL_get_grid_id(comp_id, grid_name, annotation) 
   implicit none
   integer, intent(in)                     :: comp_id
   integer                                 :: grid_id
   character(len=*), intent(in)            :: grid_name
   character(len=*), intent(in), optional  :: annotation

   if (present(annotation)) then
       call get_grid_id(comp_id, grid_name, grid_id, trim(annotation)//char(0))
   else 
       call get_grid_id(comp_id, grid_name, grid_id, trim("")//char(0))
   endif
   
   CCPL_get_grid_id = grid_id
    
   END FUNCTION CCPL_get_grid_id



   logical FUNCTION CCPL_get_H2D_grid_float_area_in_remapping_wgts(interface_id, field_index, area_array, annotation)
   implicit none
   integer, intent(in)                     :: interface_id
   integer, intent(in)                     :: field_index
   real(R4), intent(out), dimension(:)     :: area_array
   character(len=*), intent(in), optional  :: annotation
   integer                                 :: have_area

   if (present(annotation)) then
      call get_h2d_grid_area_in_remapping_weights(interface_id, field_index, area_array, size(area_array), trim("real4")//char(0), have_area, trim(annotation)//char(0))
   else
      call get_h2d_grid_area_in_remapping_weights(interface_id, field_index, area_array, size(area_array), trim("real4")//char(0), have_area, trim("")//char(0))
   endif

   CCPL_get_H2D_grid_float_area_in_remapping_wgts = .false.
   if (have_area .eq. 1) CCPL_get_H2D_grid_float_area_in_remapping_wgts = .true.

   END FUNCTION CCPL_get_H2D_grid_float_area_in_remapping_wgts



   logical FUNCTION CCPL_get_H2D_grid_double_area_in_remapping_wgts(interface_id, field_index, area_array, annotation)
   implicit none
   integer, intent(in)                     :: interface_id
   integer, intent(in)                     :: field_index
   real(R8), intent(out), dimension(:)     :: area_array
   character(len=*), intent(in), optional  :: annotation
   integer                                 :: have_area

   if (present(annotation)) then
      call get_h2d_grid_area_in_remapping_weights(interface_id, field_index, area_array, size(area_array), trim("real8")//char(0), have_area, trim(annotation)//char(0))
   else
      call get_h2d_grid_area_in_remapping_weights(interface_id, field_index, area_array, size(area_array), trim("real8")//char(0), have_area, trim("")//char(0))
   endif

   CCPL_get_H2D_grid_double_area_in_remapping_wgts = .false.
   if (have_area .eq. 1) CCPL_get_H2D_grid_double_area_in_remapping_wgts = .true.

   END FUNCTION CCPL_get_H2D_grid_double_area_in_remapping_wgts



   SUBROUTINE CCPL_get_H2D_grid_integer_data(grid_id, decomp_id, label, grid_data, annotation)
   implicit none
   integer, intent(in)                     :: grid_id
   integer, intent(in)                     :: decomp_id
   character(len=*), intent(in)            :: label
   integer, intent(inout), dimension(:)    :: grid_data
   character(len=*), intent(in), optional  :: annotation

   if (present(annotation)) then
      call get_H2D_grid_data(grid_id, decomp_id, trim(label)//char(0), trim("integer")//char(0), size(grid_data), grid_data, trim(annotation)//char(0))
   else
      call get_H2D_grid_data(grid_id, decomp_id, trim(label)//char(0), trim("integer")//char(0), size(grid_data), grid_data, trim("")//char(0))
   endif

   END SUBROUTINE CCPL_get_H2D_grid_integer_data



   SUBROUTINE CCPL_get_H2D_grid_float_data(grid_id, decomp_id, label, grid_data, annotation)
   implicit none
   integer, intent(in)                     :: grid_id
   integer, intent(in)                     :: decomp_id
   character(len=*), intent(in)            :: label
   real(R4), intent(inout), dimension(:)   :: grid_data
   character(len=*), intent(in), optional  :: annotation

   if (present(annotation)) then
      call get_H2D_grid_data(grid_id, decomp_id, trim(label)//char(0), trim("real4")//char(0), size(grid_data), grid_data, trim(annotation)//char(0))
   else
      call get_H2D_grid_data(grid_id, decomp_id, trim(label)//char(0), trim("real4")//char(0), size(grid_data), grid_data, trim("")//char(0))
   endif

   END SUBROUTINE CCPL_get_H2D_grid_float_data



   SUBROUTINE CCPL_get_H2D_grid_double_data(grid_id, decomp_id, label, grid_data, annotation)
   implicit none
   integer, intent(in)                     :: grid_id
   integer, intent(in)                     :: decomp_id
   character(len=*), intent(in)            :: label
   real(R8), intent(inout), dimension(:)   :: grid_data
   character(len=*), intent(in), optional  :: annotation

   if (present(annotation)) then
      call get_H2D_grid_data(grid_id, decomp_id, trim(label)//char(0), trim("real8")//char(0), size(grid_data), grid_data, trim(annotation)//char(0))
   else
      call get_H2D_grid_data(grid_id, decomp_id, trim(label)//char(0), trim("real8")//char(0), size(grid_data), grid_data, trim("")//char(0))
   endif

   END SUBROUTINE CCPL_get_H2D_grid_double_data



   integer FUNCTION CCPL_register_normal_parallel_decomp(decomp_name, grid_id, num_local_cells, local_cells_global_index, annotation) 
   implicit none
   character(len=*), intent(in)                :: decomp_name
   character(len=*), intent(in), optional      :: annotation
   integer,          intent(in)                :: grid_id
   integer,          intent(in)                :: num_local_cells
   integer,          intent(in), dimension(:)  :: local_cells_global_index(:)
   integer                                     :: decomp_id


   if (present(annotation)) then
       call register_parallel_decomposition(decomp_id, grid_id, num_local_cells, size(local_cells_global_index), local_cells_global_index, trim(decomp_name)//char(0), trim(annotation)//char(0))
   else 
       call register_parallel_decomposition(decomp_id, grid_id, num_local_cells, size(local_cells_global_index), local_cells_global_index, trim(decomp_name)//char(0), trim("")//char(0))
   endif

   CCPL_register_normal_parallel_decomp = decomp_id

   end FUNCTION CCPL_register_normal_parallel_decomp



   integer FUNCTION CCPL_define_single_timer(comp_id, period_unit, period_count, local_lag_count, remote_lag_count, annotation) 
   implicit none
   integer,          intent(in)                :: comp_id
   character(len=*), intent(in)                :: period_unit
   integer,          intent(in)                :: period_count
   integer,          intent(in)                :: local_lag_count
   integer,          intent(in), optional      :: remote_lag_count
   character(len=*), intent(in), optional      :: annotation
   integer                                     :: temp_remote_lag_count
   integer                                     :: timer_id
 
   temp_remote_lag_count = 0
   if (present(remote_lag_count)) temp_remote_lag_count = remote_lag_count
   if (present(annotation)) then
        call define_single_timer(comp_id, timer_id, trim(period_unit)//char(0), period_count, local_lag_count, temp_remote_lag_count, trim(annotation)//char(0))
   else
        call define_single_timer(comp_id, timer_id, trim(period_unit)//char(0), period_count, local_lag_count, temp_remote_lag_count, trim("")//char(0))
   endif
   CCPL_define_single_timer = timer_id

   end FUNCTION CCPL_define_single_timer 



   integer FUNCTION CCPL_define_complex_timer(comp_id, num_children_timers, children_timers_id, OR_or_AND, annotation)
   implicit none
   integer,          intent(in)                :: comp_id
   character(len=*), intent(in), optional      :: annotation
   integer,          intent(in)                :: OR_or_AND
   integer,          intent(in)                :: num_children_timers
   integer                                     :: timer_id
   integer,          intent(in), dimension(:)  :: children_timers_id


   if (present(annotation)) then
       call define_complex_timer(comp_id, timer_id, children_timers_id, num_children_timers, size(children_timers_id), OR_or_AND, trim(annotation)//char(0))
   else
       call define_complex_timer(comp_id, timer_id, children_timers_id, num_children_timers, size(children_timers_id), OR_or_AND, trim("")//char(0))
   endif

   CCPL_define_complex_timer = timer_id

   end FUNCTION CCPL_define_complex_timer



   SUBROUTINE CCPL_set_normal_time_step(comp_id, time_step_in_second, annotation)
   implicit none
   integer,          intent(in)                :: comp_id
   character(len=*), intent(in), optional      :: annotation
   integer,          intent(in)                :: time_step_in_second

   if (present(annotation)) then
       call set_component_time_step(comp_id, time_step_in_second, trim(annotation)//char(0))
   else
       call set_component_time_step(comp_id, time_step_in_second, trim("")//char(0))
   endif

   end SUBROUTINE CCPL_set_normal_time_step



   SUBROUTINE CCPL_reset_current_time_to_start_time(comp_id, annotation)
   implicit none
   integer,          intent(in)                :: comp_id
   character(len=*), intent(in), optional      :: annotation

   if (present(annotation)) then
       call reset_component_current_time_to_start_time(comp_id, trim(annotation)//char(0))
   else
       call reset_component_current_time_to_start_time(comp_id, trim("")//char(0))
   endif

   END SUBROUTINE CCPL_reset_current_time_to_start_time



   SUBROUTINE CCPL_advance_time(comp_id, annotation)
   implicit none
   integer,          intent(in)                :: comp_id
   character(len=*), intent(in), optional      :: annotation

   if (present(annotation)) then
        call advance_component_time(comp_id, trim(annotation)//char(0))
   else
        call advance_component_time(comp_id, trim("")//char(0))
   endif

   END SUBROUTINE CCPL_advance_time



   logical FUNCTION CCPL_is_timer_on(timer_id, annotation)
   implicit none
   integer,          intent(in)                :: timer_id
   character(len=*), intent(in), optional      :: annotation
   integer                                     :: is_on

   if ((present(annotation))) then
       call is_ccpl_timer_on(timer_id, is_on, trim(annotation)//char(0))
   else
       call is_ccpl_timer_on(timer_id, is_on, trim("")//char(0))
   endif 
   
   if (is_on .eq. 1) then
       CCPL_is_timer_on = .true.
   else
       CCPL_is_timer_on = .false.
   endif

   end FUNCTION CCPL_is_timer_on



   SUBROUTINE CCPL_check_current_time(comp_id, date, second, annotation)
   implicit none
   integer,          intent(in)                :: comp_id
   integer,          intent(in)                :: date
   integer,          intent(in)                :: second
   character(len=*), intent(in), optional      :: annotation

   if (present(annotation)) then
        call check_ccpl_component_current_time(comp_id, date, second, trim(annotation)//char(0))
   else
        call check_ccpl_component_current_time(comp_id, date, second, trim("")//char(0))
   endif

   end SUBROUTINE CCPL_check_current_time


   SUBROUTINE CCPL_finalize(to_finalize_MPI, annotation)
   implicit none
   logical,          intent(in)                :: to_finalize_MPI
   integer                                     :: local_to_finalize_MPI
   character(len=*), intent(in), optional      :: annotation
   integer ierr
   
   if (to_finalize_MPI) then
       local_to_finalize_MPI = 1
   else
       local_to_finalize_MPI = 0
   endif
   if (present(annotation)) then
       call finalize_CCPL(local_to_finalize_MPI, trim(annotation)//char(0))
   else 
       call finalize_CCPL(local_to_finalize_MPI, trim("")//char(0))
   endif

   END SUBROUTINE  CCPL_finalize


   
   logical FUNCTION CCPL_is_last_step_of_model_run(comp_id, annotation)
   implicit none
   integer,          intent(in)                :: comp_id
   character(len=*), intent(in), optional      :: annotation
   integer                                     :: is_last_step

   if (present(annotation)) then
       call check_is_ccpl_model_last_step(comp_id, is_last_step, trim(annotation)//char(0))
   else
       call check_is_ccpl_model_last_step(comp_id, is_last_step, trim("")//char(0))
   endif
   CCPL_is_last_step_of_model_run=.false.
   if (is_last_step .eq. 1) CCPL_is_last_step_of_model_run = .false.

   END FUNCTION CCPL_is_last_step_of_model_run


   logical FUNCTION CCPL_is_model_run_ended(comp_id, annotation)
   implicit none
   integer,          intent(in)                :: comp_id
   character(len=*), intent(in), optional      :: annotation
   integer                                     :: is_ended

   if (present(annotation)) then
       call check_is_ccpl_model_run_ended(comp_id, is_ended, trim(annotation)//char(0))
   else
       call check_is_ccpl_model_run_ended(comp_id, is_ended, trim("")//char(0))
   endif

   if (is_ended .eq. 1) then
       CCPL_is_model_run_ended = .true.
   else 
       CCPL_is_model_run_ended = .false.
   endif

   end FUNCTION CCPL_is_model_run_ended



   integer FUNCTION CCPL_register_normal_remap_interface(interface_name, num_field_instances, field_instance_IDs_source, field_instance_IDs_target, timer_ID, inst_or_aver, annotation)
   implicit none
   character(len=*), intent(in)                :: interface_name
   character(len=*), intent(in), optional      :: annotation
   integer,          intent(in)                :: timer_ID
   integer,          intent(in)                :: inst_or_aver
   integer,          intent(in), dimension(:)  :: field_instance_IDs_source
   integer,          intent(in), dimension(:)  :: field_instance_IDs_target
   integer,          intent(in)                :: num_field_instances
   integer                                     :: interface_id

   if (present(annotation)) then
       call register_normal_remap_interface(trim(interface_name)//char(0), interface_id, num_field_instances, field_instance_IDs_source, field_instance_IDs_target, timer_ID, inst_or_aver, size(field_instance_IDs_source), size(field_instance_IDs_target), trim(annotation)//char(0))
   else
       call register_normal_remap_interface(trim(interface_name)//char(0), interface_id, num_field_instances, field_instance_IDs_source, field_instance_IDs_target, timer_ID, inst_or_aver, size(field_instance_IDs_source), size(field_instance_IDs_target), trim("")//char(0))
   endif
   CCPL_register_normal_remap_interface = interface_id;

   END FUNCTION CCPL_register_normal_remap_interface



   integer FUNCTION CCPL_register_remap_interface_with_float_frac(interface_name, num_field_instances, field_instance_IDs_source, field_instance_IDs_target, timer_ID, inst_or_aver, frac_src, frac_dst, annotation)
   implicit none
   character(len=*), intent(in)                                 :: interface_name
   character(len=*), intent(in), optional                       :: annotation
   integer,          intent(in)                                 :: timer_ID
   integer,          intent(in)                                 :: inst_or_aver
   integer,          intent(in), dimension(:)                   :: field_instance_IDs_source
   integer,          intent(in), dimension(:)                   :: field_instance_IDs_target
   integer,          intent(in)                                 :: num_field_instances
   integer                                                      :: interface_id, size_frac_dst
   real(R4),         INTENT(IN), dimension(:)                   :: frac_src
   real(R4),         intent(in), dimension(:), target, optional :: frac_dst
   real(R4),                     dimension(:), pointer          :: temp_frac_dst, temp_float_1d
   character *2048                                               :: local_annotation
   
   allocate(temp_float_1d(1))
   temp_frac_dst => temp_float_1d
   size_frac_dst = -1
   if (present(frac_dst)) then
       temp_frac_dst => frac_dst
       size_frac_dst = size(frac_dst)
   endif
   
   local_annotation = ""
   if (present(annotation)) local_annotation = annotation

   call register_frac_based_remap_interface(trim(interface_name)//char(0), interface_id, num_field_instances, field_instance_IDs_source, field_instance_IDs_target, timer_ID, &
                                            inst_or_aver, size(field_instance_IDs_source), size(field_instance_IDs_target), loc(frac_src), loc(temp_frac_dst), size(frac_src), &
                                            size_frac_dst, trim("real4")//char(0), trim(local_annotation)//char(0))

   deallocate(temp_float_1d)

   CCPL_register_remap_interface_with_float_frac = interface_id

   END FUNCTION CCPL_register_remap_interface_with_float_frac



   integer FUNCTION CCPL_register_remap_interface_with_double_frac(interface_name, num_field_instances, field_instance_IDs_source, field_instance_IDs_target, timer_ID, inst_or_aver, frac_src, frac_dst, annotation)
   implicit none
   character(len=*), intent(in)                                 :: interface_name
   character(len=*), intent(in), optional                       :: annotation
   integer,          intent(in)                                 :: timer_ID
   integer,          intent(in)                                 :: inst_or_aver
   integer,          intent(in), dimension(:)                   :: field_instance_IDs_source
   integer,          intent(in), dimension(:)                   :: field_instance_IDs_target
   integer,          intent(in)                                 :: num_field_instances
   integer                                                      :: interface_id, size_frac_dst
   real(R8),         INTENT(IN), dimension(:)                   :: frac_src
   real(R8),         intent(in), dimension(:), target, optional :: frac_dst
   real(R8),                     dimension(:), pointer          :: temp_frac_dst, temp_double_1d
   character *2048                                               :: local_annotation
   
   allocate(temp_double_1d(1))
   temp_frac_dst => temp_double_1d
   size_frac_dst = -1
   if (present(frac_dst)) then
       temp_frac_dst => frac_dst
       size_frac_dst = size(frac_dst)
   endif
   
   local_annotation = ""
   if (present(annotation)) local_annotation = annotation

   call register_frac_based_remap_interface(trim(interface_name)//char(0), interface_id, num_field_instances, field_instance_IDs_source, field_instance_IDs_target, timer_ID, &
                                            inst_or_aver, size(field_instance_IDs_source), size(field_instance_IDs_target), loc(frac_src), loc(temp_frac_dst), size(frac_src), &
                                            size_frac_dst, trim("real8")//char(0), trim(local_annotation)//char(0))

   deallocate(temp_double_1d)

   CCPL_register_remap_interface_with_double_frac = interface_id

   END FUNCTION CCPL_register_remap_interface_with_double_frac



   logical FUNCTION CCPL_check_is_import_field_connected(import_interface_id, field_instance_id, annotation)
   implicit none
   integer,          intent(in)                         :: import_interface_id
   integer,          intent(in)                         :: field_instance_id
   character(len=*), intent(in), optional               :: annotation
   character *2048                                      :: local_annotation
   integer                                              :: check_result

   local_annotation = ""
   if (present(annotation)) local_annotation = annotation
   call check_is_ccpl_import_field_connected(import_interface_id, field_instance_id, check_result, trim(local_annotation)//char(0))
   CCPL_check_is_import_field_connected = .true.
   if (check_result .eq. 0) CCPL_check_is_import_field_connected = .false.
   
   END FUNCTION CCPL_check_is_import_field_connected
   


   SUBROUTINE CCPL_get_import_fields_sender_time(import_interface_id, sender_date, sender_second, sender_elapased_days, annotation)
   implicit none
   integer,          intent(in)                          :: import_interface_id
   integer,          intent(out), dimension(:)           :: sender_date
   integer,          intent(out), dimension(:)           :: sender_second
   integer,          intent(out), dimension(:),optional  :: sender_elapased_days ! number of elapased days since 0000-01-01
   character(len=*), intent(in), optional                :: annotation
   character *2048                                       :: local_annotation
   integer                                               :: temp_sender_elapased_days(4096) 

   local_annotation = ""
   if (present(annotation)) then
       local_annotation = annotation
   endif
   if (present(sender_elapased_days)) then
      call get_ccpl_import_fields_sender_time(import_interface_id, size(sender_date), size(sender_elapased_days), size(sender_second), sender_date, sender_elapased_days, sender_second, trim(local_annotation)//char(0))
   else
      call get_ccpl_import_fields_sender_time(import_interface_id, size(sender_date), size(temp_sender_elapased_days), size(sender_second), sender_date, temp_sender_elapased_days, sender_second, trim(local_annotation)//char(0))
   endif
   
   END SUBROUTINE CCPL_get_import_fields_sender_time



   integer FUNCTION CCPL_register_import_interface(interface_name, num_field_instances, field_instance_IDs, timer_ID, inst_or_aver, necessity, annotation)
   implicit none
   character(len=*), intent(in)                         :: interface_name
   character(len=*), intent(in), optional               :: annotation
   integer,          intent(in)                         :: timer_ID
   integer,          intent(in)                         :: inst_or_aver
   integer,          intent(in), dimension(:)           :: field_instance_IDs
   integer,          intent(in), dimension(:), optional :: necessity
   integer,          intent(in)                         :: num_field_instances
   integer                                              :: interface_id
   integer                                              :: local_necessity(2)
   character *2048                                      :: local_annotation
   
   local_annotation = ""
   if (present(annotation)) local_annotation = annotation
   
   call register_inout_interface(trim(interface_name)//char(0), interface_id, 0, num_field_instances, field_instance_IDs, timer_ID, inst_or_aver, trim(local_annotation)//char(0), size(field_instance_IDs))
   if (present(necessity)) then
      call set_import_interface_fields_necessity(interface_id, necessity, size(necessity), trim(local_annotation)//char(0))
   endif
   CCPL_register_import_interface = interface_id;

   end FUNCTION CCPL_register_import_interface



   integer FUNCTION CCPL_register_export_interface(interface_name, num_field_instances, field_instance_IDs, timer_ID, annotation)
   implicit none
   character(len=*), intent(in)                :: interface_name
   character(len=*), intent(in), optional      :: annotation
   integer,          intent(in)                :: timer_ID
   integer,          intent(in), dimension(:)  :: field_instance_IDs
   integer,          intent(in)                :: num_field_instances
   integer                                     :: interface_id

   
   if (present(annotation)) then
       call register_inout_interface(trim(interface_name)//char(0), interface_id, 1, num_field_instances, field_instance_IDs, timer_ID, field_instance_IDs, trim(annotation)//char(0), size(field_instance_IDs))
   else
       call register_inout_interface(trim(interface_name)//char(0), interface_id, 1, num_field_instances, field_instance_IDs, timer_ID, field_instance_IDs, trim("")//char(0), size(field_instance_IDs))
   endif
   CCPL_register_export_interface = interface_id;

   end FUNCTION CCPL_register_export_interface



   logical FUNCTION CCPL_execute_interface_using_id(interface_id, bypass_timer, field_update_status, annotation)
   implicit none
   integer,          intent(in)                          :: interface_id
   logical,          intent(in)                          :: bypass_timer
   character(len=*), intent(in), optional                :: annotation
   integer                                               :: local_bypass_timer
   integer,          intent(out), dimension(:), optional :: field_update_status
   integer                                               :: temp_field_update_status(4096), i, num_dst_fields
   character *2048                                       :: local_annotation


   if (bypass_timer) then
       local_bypass_timer = 1
   else
       local_bypass_timer = 0
   endif
   local_annotation = ""
   if (present(annotation)) then
       local_annotation = annotation
   endif

   if (present(field_update_status)) then
       call execute_inout_interface_with_id(interface_id, local_bypass_timer, temp_field_update_status, size(field_update_status), num_dst_fields, trim(local_annotation)//char(0))
       field_update_status(1:num_dst_fields) = temp_field_update_status(1:num_dst_fields)
   else 
       call execute_inout_interface_with_id(interface_id, local_bypass_timer, temp_field_update_status, 4000, num_dst_fields, trim(local_annotation)//char(0))
   endif

   CCPL_execute_interface_using_id = .true.

   END FUNCTION CCPL_execute_interface_using_id



   logical FUNCTION CCPL_execute_interface_using_name(component_id, interface_name, bypass_timer, field_update_status, annotation)
   implicit none
   integer,          intent(in)                          :: component_id
   logical,          intent(in)                          :: bypass_timer
   character(len=*), intent(in), optional                :: annotation
   character(len=*), intent(in)                          :: interface_name
   integer                                               :: local_bypass_timer
   integer,          intent(out), dimension(:), optional :: field_update_status
   integer                                               :: temp_field_update_status(4096), i, num_dst_fields
   character *2048                                        :: local_annotation


   if (bypass_timer) then
       local_bypass_timer = 1
   else
       local_bypass_timer = 0
   endif
   local_annotation = ""
   if (present(annotation)) then
       local_annotation = annotation
   endif

   if (present(field_update_status)) then
       call execute_inout_interface_with_name(component_id, trim(interface_name)//char(0), local_bypass_timer, temp_field_update_status, size(field_update_status), num_dst_fields, trim(local_annotation)//char(0))
       field_update_status(1:num_dst_fields) = temp_field_update_status(1:num_dst_fields)
   else 
       call execute_inout_interface_with_name(component_id, trim(interface_name)//char(0), local_bypass_timer, temp_field_update_status, 4000, num_dst_fields, trim(local_annotation)//char(0))
   endif

   CCPL_execute_interface_using_name = .true.

   END FUNCTION CCPL_execute_interface_using_name



   SUBROUTINE CCPL_get_local_comp_full_name(comp_id, comp_full_name, annotation)
   implicit none
   integer,          intent(in)                :: comp_id
   character(len=*), intent(out)               :: comp_full_name
   character(len=*), intent(in), optional      :: annotation

   if (present(annotation)) then
       call get_local_comp_full_name(comp_id, comp_full_name, len(comp_full_name), trim(annotation)//char(0))
   else
       call get_local_comp_full_name(comp_id, comp_full_name, len(comp_full_name), trim("")//char(0))
   endif

   END SUBROUTINE CCPL_get_local_comp_full_name



   SUBROUTINE CCPL_report_log(comp_id, condition, report_string, annotation)
   implicit none
   integer,          intent(in)                :: comp_id
   logical,          intent(in)                :: condition
   character(len=*), intent(in)                :: report_string
   character(len=*), intent(in), optional      :: annotation
   integer                                     :: local_condition

   local_condition = 0
   if (condition) local_condition = 1

   if (present(annotation)) then
       call CCPL_report(5, comp_id, local_condition, trim(report_string)//char(0), trim(annotation)//char(0))
   else 
       call CCPL_report(5, comp_id, local_condition, trim(report_string)//char(0), trim("")//char(0))
   endif

   END SUBROUTINE CCPL_report_log



   SUBROUTINE CCPL_report_progress(comp_id, condition, report_string, annotation)
   implicit none
   integer,          intent(in)                :: comp_id
   logical,          intent(in)                :: condition
   character(len=*), intent(in)                :: report_string
   character(len=*), intent(in), optional      :: annotation
   integer                                     :: local_condition

   local_condition = 0
   if (condition) local_condition = 1

   if (present(annotation)) then
       call CCPL_report(4, comp_id, local_condition, trim(report_string)//char(0), trim(annotation)//char(0))
   else 
       call CCPL_report(4, comp_id, local_condition, trim(report_string)//char(0), trim("")//char(0))
   endif

   END SUBROUTINE CCPL_report_progress



   SUBROUTINE CCPL_report_error(comp_id, condition, report_string, annotation)
   implicit none
   integer,          intent(in)                :: comp_id
   logical,          intent(in)                :: condition
   character(len=*), intent(in)                :: report_string
   character(len=*), intent(in), optional      :: annotation
   integer                                     :: local_condition

   local_condition = 0
   if (condition) local_condition = 1

   if (present(annotation)) then
       call CCPL_report(1, comp_id, local_condition, trim(report_string)//char(0), trim(annotation)//char(0))
   else 
       call CCPL_report(1, comp_id, local_condition, trim(report_string)//char(0), trim("")//char(0))
   endif

   END SUBROUTINE CCPL_report_error



   SUBROUTINE CCPL_do_restart_write_IO(comp_id, bypass_timer, bypass_import_fields, annotation)
   implicit none
   integer,          intent(in)                :: comp_id
   logical,          intent(in)                :: bypass_timer
   logical,          intent(in), optional      :: bypass_import_fields
   character(len=*), intent(in), optional      :: annotation
   integer                                     :: local_bypass_timer
   integer                                     :: local_bypass_import_fields

   local_bypass_timer = 0
   if (bypass_timer) local_bypass_timer = 1
   local_bypass_import_fields = 0
   if (present(bypass_import_fields) .and. bypass_import_fields) local_bypass_import_fields = 1

   if (present(annotation)) then
       call CCPL_write_restart(comp_id, local_bypass_timer, local_bypass_import_fields, trim(annotation)//char(0))
   else 
       call CCPL_write_restart(comp_id, local_bypass_timer, local_bypass_import_fields, trim("")//char(0))
   endif
   
   END SUBROUTINE CCPL_do_restart_write_IO
 


   SUBROUTINE CCPL_do_family_coupling_generation(comp_id, annotation)
   implicit none
   integer,          intent(in)                :: comp_id
   character(len=*), intent(in), optional      :: annotation

   if (present(annotation)) then
       call ccpl_family_coupling_generation(comp_id, trim(annotation)//char(0))
   else
       call ccpl_family_coupling_generation(comp_id, trim("")//char(0))
   endif
   END SUBROUTINE CCPL_do_family_coupling_generation



   SUBROUTINE CCPL_do_individual_coupling_generation(comp_id, annotation)
   implicit none
   integer,          intent(in)                :: comp_id
   character(len=*), intent(in), optional      :: annotation

   if (present(annotation)) then
       call ccpl_individual_coupling_generation(comp_id, trim(annotation)//char(0))
   else
       call ccpl_individual_coupling_generation(comp_id, trim("")//char(0))
   endif
   END SUBROUTINE CCPL_do_individual_coupling_generation


   
   SUBROUTINE CCPL_get_configurable_comps_full_names(comp_id, keyword, num_comps, comps_full_names, individual_or_family, annotation)
   integer,          intent(in)                 :: comp_id
   character(len=*), intent(in)                 :: keyword
   integer,          intent(out)                :: num_comps
   character(len=*), intent(out)                :: comps_full_names(:)
   integer,          intent(out), optional      :: individual_or_family(:)
   character(len=*), intent(in), optional       :: annotation
   character *2048                              :: local_annotation
   integer                                      :: size_comps_full_names, size_individual_or_family, i, local_individual_or_family
   

   local_annotation = ""
   if (present(annotation)) local_annotation = annotation
   size_comps_full_names = size(comps_full_names)
   size_individual_or_family = -1
   if (present(individual_or_family)) size_individual_or_family=size(individual_or_family)

   call ccpl_load_comps_full_names_from_config_file(comp_id, trim(keyword)//char(0), size_comps_full_names, size_individual_or_family, num_comps, trim(local_annotation)//char(0))
   do i = 1, num_comps
      call ccpl_get_one_comp_full_name(comp_id, keyword, len(comps_full_names(i)), i, comps_full_names(i), local_individual_or_family, trim(local_annotation)//char(0)) 
      if (present(individual_or_family)) individual_or_family(i) = local_individual_or_family
   enddo
   call ccpl_finish_getting_configurable_comps_full_names(comp_id, trim(local_annotation)//char(0)) 

   END SUBROUTINE CCPL_get_configurable_comps_full_names
   


   SUBROUTINE CCPL_do_external_coupling_generation(num_comps, comps_full_names, individual_or_family, annotation)
   implicit none
   integer,          intent(in)                :: num_comps
   character(len=*), intent(in)                :: comps_full_names(:)
   integer,          intent(in), optional      :: individual_or_family(:)
   character(len=*), intent(in), optional      :: annotation
   integer                                     :: size_comps_full_names, i, j
   character *2048                             :: local_annotation
   integer,          allocatable               :: local_individual_or_family(:)
   integer                                     :: size_individual_or_family
   

   if (num_comps > 0) then
       allocate(local_individual_or_family(num_comps))
   else
       allocate(local_individual_or_family(100))
   endif
   local_individual_or_family(:) = 1 
   size_individual_or_family = num_comps
   if (present(individual_or_family)) then
      size_individual_or_family = size(individual_or_family)
      j = num_comps
      if (size_individual_or_family .lt. j) j = size_individual_or_family
      do i = 1, j
         local_individual_or_family(i) = individual_or_family(i)
      enddo
   endif

   local_annotation = ""
   if (present(annotation)) local_annotation = annotation
   
   size_comps_full_names = size(comps_full_names)
   call ccpl_begin_external_coupling_generation(num_comps, size_comps_full_names, size_individual_or_family, trim(local_annotation)//char(0))
   do i = 1, num_comps
       call ccpl_add_comp_for_external_coupling_generation(trim(comps_full_names(i))//char(0), local_individual_or_family(i), trim(local_annotation)//char(0)) 
   enddo 
   call ccpl_end_external_coupling_generation(trim(local_annotation)//char(0))
   write(*,*) "size of comps_full_names is ", size_comps_full_names

   deallocate(local_individual_or_family)
   
   END SUBROUTINE CCPL_do_external_coupling_generation



   SUBROUTINE CCPL_start_restart_read_IO(comp_id, specified_restart_file, annotation)
   implicit none
   integer,          intent(in)                :: comp_id
   character(len=*), intent(in), optional      :: specified_restart_file
   character(len=*), intent(in), optional      :: annotation
   character *1024                             :: local_specified_restart_file


   local_specified_restart_file = ""
   if (present(specified_restart_file)) local_specified_restart_file = specified_restart_file

   if (present(annotation)) then
       call CCPL_read_restart(comp_id, trim(local_specified_restart_file)//char(0), trim(annotation)//char(0))
   else 
       call CCPL_read_restart(comp_id, trim(local_specified_restart_file)//char(0), trim("")//char(0))
   endif
   
   END SUBROUTINE CCPL_start_restart_read_IO
 


   SUBROUTINE CCPL_restart_read_fields_all(comp_id, annotation)
   implicit none
   integer,          intent(in)                :: comp_id
   character(len=*), intent(in), optional      :: annotation


   if (present(annotation)) then
       call CCPL_read_all_restart_fields(comp_id, trim(annotation)//char(0))
   else
       call CCPL_read_all_restart_fields(comp_id, trim("")//char(0))
   endif

   END SUBROUTINE CCPL_restart_read_fields_all



   SUBROUTINE CCPL_restart_read_fields_interface(interface_id, annotation)
   implicit none
   integer,          intent(in)                :: interface_id
   character(len=*), intent(in), optional      :: annotation


   if (present(annotation)) then
       call CCPL_read_import_interface_restart_fields(interface_id, trim(annotation)//char(0))
   else
       call CCPL_read_import_interface_restart_fields(interface_id, trim("")//char(0))
   endif

   END SUBROUTINE CCPL_restart_read_fields_interface



   SUBROUTINE CCPL_get_restart_setting(comp_id, restart_date, restart_second, original_case_name, run_type, annotation)
   implicit none
   integer,          intent(in)                 :: comp_id
   integer,          intent(out)                :: restart_date
   integer,          intent(out)                :: restart_second
   character(len=*), intent(out), optional      :: original_case_name
   character(len=*), intent(out), optional      :: run_type
   character(len=*), intent(in),  optional      :: annotation
   character *1024                              :: local_annotation


   local_annotation = ""
   if (present(annotation)) local_annotation = annotation
   call get_ccpl_restart_time(comp_id, restart_date, restart_second, trim(local_annotation)//char(0))
   if (present(original_case_name)) call get_ccpl_original_case_name(comp_id, len(original_case_name), original_case_name, trim(local_annotation)//char(0))
   if (present(run_type)) call get_ccpl_run_type(comp_id, len(run_type), run_type, trim(local_annotation)//char(0))

   END SUBROUTINE CCPL_get_restart_setting



   logical FUNCTION CCPL_is_restart_timer_on(comp_id, annotation)
   implicit none
   integer,          intent(in)                :: comp_id
   character(len=*), intent(in), optional      :: annotation
   integer                                     :: check_result

   if (present(annotation)) then
       call is_restart_timer_on(comp_id, check_result, trim(annotation)//char(0))
   else
       call is_restart_timer_on(comp_id, check_result, trim("")//char(0))
   endif   

   CCPL_is_restart_timer_on = .false.
   if (check_result .eq. 1) CCPL_is_restart_timer_on = .true.

   END FUNCTION CCPL_is_restart_timer_on



   SUBROUTINE CCPL_abort(error_string)
   implicit none
   character(len=*),     intent(in)    ::  error_string

   call coupling_abort(trim(error_string)//char(0))

   END SUBROUTINE CCPL_abort


 END MODULE CCPL_interface_mod
