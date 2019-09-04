/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#include <mpi.h>
#include "global_data.h"
#include "cor_global_data.h"
#include "remap_mgt.h"
#include "quick_sort.h"
#include <stdio.h>
#include <string.h>
#include <sys/time.h>
#include "coupling_interface.h"


int coupling_process_control_counter = 0;


void check_for_component_registered(int comp_id, int API_ID, const char *annotation, bool enable_minus_1)
{
    char API_label[NAME_STR_SIZE];
    

    get_API_hint(-1, API_ID, API_label);
    check_for_ccpl_managers_allocated(API_ID, annotation);

    if (comp_id == -1)
        EXECUTION_REPORT(REPORT_ERROR, comp_id, enable_minus_1, "Error happens when calling the API \"%s\": the given component model ID (-1) is not valid. Please check the model code with the annotation \"%s\"", API_label, annotation);
    else EXECUTION_REPORT(REPORT_ERROR, -1, comp_comm_group_mgt_mgr->is_legal_local_comp_id(comp_id, true) && comp_comm_group_mgt_mgr->get_current_proc_id_in_comp(comp_id, "in check_for_component_registered") >= 0, "Error happens when calling the API \"%s\": the given component ID (0x%d) is wrong. Please check the model code with the annotation \"%s\"", API_label, comp_id, annotation);
}


void copy_out_string_to_Fortran_API(int comp_id, int size_API_string, char *API_string, const char *CCPL_string, int API_id, const char *parameter_name, const char *annotation)
{
    char API_label[NAME_STR_SIZE];


    get_API_hint(comp_id, API_id,API_label);
    
    EXECUTION_REPORT(REPORT_ERROR, comp_id, size_API_string >= strlen(CCPL_string), "Error happens when calling the API \"%s\": the parameter string \"%s\" is too short: only %d while the required size is %d. Please verify the model code with the annotation \"%s\"", API_label, parameter_name, size_API_string, strlen(CCPL_string), annotation);
    strncpy(API_string, CCPL_string, strlen(CCPL_string));
    for (int i = strlen(CCPL_string); i < size_API_string; i ++)
        API_string[i] = ' ';    
}


#ifdef LINK_WITHOUT_UNDERLINE
extern "C" void finalize_ccpl
#else
extern "C" void finalize_ccpl_
#endif
(int *to_finalize_MPI, const char *annotation)
{
    if (comp_comm_group_mgt_mgr->get_current_proc_global_id() == 0)
        EXECUTION_REPORT(REPORT_PROGRESS, -1, true, "Start to finalize C-Coupler at the model code with the annotation \"%s\"", annotation);

    comp_comm_group_mgt_mgr->output_performance_timing();
    inout_interface_mgr->free_all_MPI_wins();

    delete annotation_mgr;
    delete decomps_info_mgr;
    delete decomp_grids_mgr;
    delete components_time_mgrs;
    delete timer_mgr;
    delete inout_interface_mgr;
    delete routing_info_mgr;
    delete IO_fields_mgr;
    delete components_IO_output_procedures_mgr;
    delete fields_gather_scatter_mgr;
    delete remapping_configuration_mgr;
    delete runtime_remapping_weights_mgr;
    delete fields_info;
    delete original_grid_mgr;
    delete memory_manager;
    delete coupling_generator;
    delete comp_comm_group_mgt_mgr;
    comp_comm_group_mgt_mgr = NULL;

    if (*to_finalize_MPI == 0)
        return;
    int flag;
    MPI_Finalized(&flag);
    if (!flag)
        MPI_Finalize();
}


#ifdef LINK_WITHOUT_UNDERLINE
extern "C" void get_ccpl_double_current_calendar_time
#else
extern "C" void get_ccpl_double_current_calendar_time_
#endif
(int *comp_id, double *cal_time, int *shift_seconds, const char *annotation)
{
    check_for_component_registered(*comp_id, API_ID_TIME_MGT_GET_CURRENT_CAL_TIME, annotation, false);
    *cal_time = components_time_mgrs->get_time_mgr(*comp_id)->get_double_current_calendar_time(*shift_seconds, annotation);
}


#ifdef LINK_WITHOUT_UNDERLINE
extern "C" void get_ccpl_float_current_calendar_time
#else
extern "C" void get_ccpl_float_current_calendar_time_
#endif
(int *comp_id, float *cal_time, int *shift_seconds, const char *annotation)
{
    check_for_component_registered(*comp_id, API_ID_TIME_MGT_GET_CURRENT_CAL_TIME, annotation, false);
    *cal_time = components_time_mgrs->get_time_mgr(*comp_id)->get_float_current_calendar_time(*shift_seconds, annotation);
}


#ifdef LINK_WITHOUT_UNDERLINE
extern "C" void get_ccpl_current_date
#else
extern "C" void get_ccpl_current_date_
#endif
(int *comp_id, int *date, const char *annotation)
{
    check_for_component_registered(*comp_id, API_ID_TIME_MGT_GET_CURRENT_DATE, annotation, false);
    *date = components_time_mgrs->get_time_mgr(*comp_id)->get_current_date();
}


#ifdef LINK_WITHOUT_UNDERLINE
extern "C" void get_ccpl_current_second
#else
extern "C" void get_ccpl_current_second_
#endif
(int *comp_id, int *second, const char *annotation)
{
    check_for_component_registered(*comp_id, API_ID_TIME_MGT_GET_CURRENT_SECOND, annotation, false);
    *second = components_time_mgrs->get_time_mgr(*comp_id)->get_current_second();
}


#ifdef LINK_WITHOUT_UNDERLINE
extern "C" void is_comp_first_step
#else
extern "C" void is_comp_first_step_
#endif
(int *comp_id, int *result, const char *annotation)
{
    check_for_component_registered(*comp_id, API_ID_TIME_MGT_IS_FIRST_STEP, annotation, false);
    *result = components_time_mgrs->get_time_mgr(*comp_id)->get_current_num_time_step() == 0? 1 : 0;
}


#ifdef LINK_WITHOUT_UNDERLINE
extern "C" void is_comp_first_restart_step
#else
extern "C" void is_comp_first_restart_step_
#endif
(int *comp_id, int *result, const char *annotation)
{
    check_for_component_registered(*comp_id, API_ID_TIME_MGT_IS_FIRST_RESTART_STEP, annotation, false);
    *result = components_time_mgrs->get_time_mgr(*comp_id)->is_first_restart_step()? 1 : 0;
}


#ifdef LINK_WITHOUT_UNDERLINE
extern "C" void get_ccpl_current_number_of_step
#else
extern "C" void get_ccpl_current_number_of_step_
#endif
(int *comp_id, int *nstep, const char *annotation)
{
    check_for_component_registered(*comp_id, API_ID_TIME_MGT_GET_NUM_CURRENT_STEP, annotation, false);
    *nstep = components_time_mgrs->get_time_mgr(*comp_id)->get_current_num_time_step();
}


#ifdef LINK_WITHOUT_UNDERLINE
extern "C" void get_ccpl_num_total_step
#else
extern "C" void get_ccpl_num_total_step_
#endif
(int *comp_id, int *nstep, const char *annotation)
{
    check_for_component_registered(*comp_id, API_ID_TIME_MGT_GET_NUM_TOTAL_STEPS, annotation, false);
    *nstep = (int) components_time_mgrs->get_time_mgr(*comp_id)->get_num_total_step();
}


#ifdef LINK_WITHOUT_UNDERLINE
extern "C" void get_ccpl_time_step
#else
extern "C" void get_ccpl_time_step_
#endif
(int *comp_id, int *time_step, const char *annotation)
{
    check_for_component_registered(*comp_id, API_ID_TIME_MGT_GET_NORMAL_TIME_STEP, annotation, false);
    *time_step = components_time_mgrs->get_time_mgr(*comp_id)->get_time_step_in_second();
}


#ifdef LINK_WITHOUT_UNDERLINE
extern "C" void get_ccpl_start_time
#else
extern "C" void get_ccpl_start_time_
#endif
(int *comp_id, int *year, int *month, int *day, int *seconds, const char *annotation)
{
    check_for_component_registered(*comp_id, API_ID_TIME_MGT_GET_START_TIME, annotation, false);

    *year = components_time_mgrs->get_time_mgr(*comp_id)->get_start_full_time() / 1000000000;
    *month = (components_time_mgrs->get_time_mgr(*comp_id)->get_start_full_time() / 10000000)%100;
    *day = (components_time_mgrs->get_time_mgr(*comp_id)->get_start_full_time() / 100000)%100;
    *seconds = components_time_mgrs->get_time_mgr(*comp_id)->get_start_full_time() % 100000;
}


#ifdef LINK_WITHOUT_UNDERLINE
extern "C" void get_ccpl_stop_time
#else
extern "C" void get_ccpl_stop_time_
#endif
(int *comp_id, int *year, int *month, int *day, int *second, const char *annotation)
{
    check_for_component_registered(*comp_id, API_ID_TIME_MGT_GET_STOP_TIME, annotation, false);

    *year = components_time_mgrs->get_time_mgr(*comp_id)->get_stop_year();
    *month = components_time_mgrs->get_time_mgr(*comp_id)->get_stop_month();
    *day = components_time_mgrs->get_time_mgr(*comp_id)->get_stop_day();
    *second = components_time_mgrs->get_time_mgr(*comp_id)->get_stop_second();
}


#ifdef LINK_WITHOUT_UNDERLINE
extern "C" void get_ccpl_previous_time
#else
extern "C" void get_ccpl_previous_time_
#endif
(int *comp_id, int *year, int *month, int *day, int *seconds, const char *annotation)
{
    check_for_component_registered(*comp_id, API_ID_TIME_MGT_GET_PREVIOUS_TIME, annotation, false);

    *year = components_time_mgrs->get_time_mgr(*comp_id)->get_previous_full_time() / 1000000000;
    *month = (components_time_mgrs->get_time_mgr(*comp_id)->get_previous_full_time() / 10000000)%100;
    *day = (components_time_mgrs->get_time_mgr(*comp_id)->get_previous_full_time() / 100000)%100;
    *seconds = components_time_mgrs->get_time_mgr(*comp_id)->get_previous_full_time() % 100000;
}


#ifdef LINK_WITHOUT_UNDERLINE
extern "C" void get_ccpl_current_time
#else
extern "C" void get_ccpl_current_time_
#endif
(int *comp_id, int *year, int *month, int *day, int *second, int *shift_second, const char *annotation)
{
    check_for_component_registered(*comp_id, API_ID_TIME_MGT_GET_CURRENT_TIME, annotation, false);
    components_time_mgrs->get_time_mgr(*comp_id)->get_current_time(*year, *month, *day, *second, *shift_second, annotation);
}


#ifdef LINK_WITHOUT_UNDERLINE
extern "C" void get_ccpl_current_num_days_in_year
#else
extern "C" void get_ccpl_current_num_days_in_year_
#endif
(int *comp_id, int *days, const char *annotation)
{
    check_for_component_registered(*comp_id, API_ID_TIME_MGT_GET_CURRENT_NUM_DAYS_IN_YEAR, annotation, false);
    *days = components_time_mgrs->get_time_mgr(*comp_id)->get_current_num_days_in_year();
}


#ifdef LINK_WITHOUT_UNDERLINE
extern "C" void get_ccpl_current_year
#else
extern "C" void get_ccpl_current_year_
#endif
(int *comp_id, int *year, const char *annotation)
{
    check_for_component_registered(*comp_id, API_ID_TIME_MGT_GET_CURRENT_YEAR, annotation, false);
    *year = components_time_mgrs->get_time_mgr(*comp_id)->get_current_year();
}


#ifdef LINK_WITHOUT_UNDERLINE
extern "C" void get_ccpl_num_elapsed_days_from_start_date
#else
extern "C" void get_ccpl_num_elapsed_days_from_start_date_
#endif
(int *comp_id, int *days, int *seconds, const char *annotation)
{
    check_for_component_registered(*comp_id, API_ID_TIME_MGT_GET_ELAPSED_DAYS_FROM_START, annotation, false);
    components_time_mgrs->get_time_mgr(*comp_id)->get_elapsed_days_from_start_date(days, seconds);
}


#ifdef LINK_WITHOUT_UNDERLINE
extern "C" void get_ccpl_num_elapsed_days_from_reference_date
#else
extern "C" void get_ccpl_num_elapsed_days_from_reference_date_
#endif
(int *comp_id, int *days, int *seconds, const char *annotation)
{
    check_for_component_registered(*comp_id, API_ID_TIME_MGT_GET_ELAPSED_DAYS_FROM_REF, annotation, false);
    components_time_mgrs->get_time_mgr(*comp_id)->get_elapsed_days_from_reference_date(days, seconds);
}


#ifdef LINK_WITHOUT_UNDERLINE
extern "C" void coupling_abort
#else
extern "C" void coupling_abort_
#endif
(const char *error_string)
{
        EXECUTION_REPORT(REPORT_ERROR,-1, false, error_string);
}


#ifdef LINK_WITHOUT_UNDERLINE
extern "C" void initialize_ccpl_mgrs
#else
extern "C" void initialize_ccpl_mgrs_
#endif
()
{
    execution_phase_number = 1;
    annotation_mgr = new Annotation_mgt();
    decomps_info_mgr = new Decomp_info_mgt();
    decomp_grids_mgr = new Decomp_grid_mgt();
    memory_manager = new Memory_mgt();
    components_time_mgrs = new Components_time_mgt();
    timer_mgr = new Timer_mgt();
    execution_phase_number = 2;
    inout_interface_mgr = new Inout_interface_mgt();
    IO_fields_mgr = new IO_field_mgt();
    components_IO_output_procedures_mgr = new Components_IO_output_procedures_mgt();
    fields_gather_scatter_mgr = new Fields_gather_scatter_mgt();
    remapping_configuration_mgr = new Remapping_configuration_mgt();
    routing_info_mgr = new Routing_info_mgt();
    runtime_remapping_weights_mgr = new Runtime_remapping_weights_mgt();
    all_H2D_remapping_wgt_files_info = new H2D_remapping_wgt_file_container();
    coupling_generator = new Coupling_generator();
}


#ifdef LINK_WITHOUT_UNDERLINE
extern "C" void check_fortran_api_int_type
#else
extern "C" void check_fortran_api_int_type_
#endif
(int *fortran_int_size)
{
    EXECUTION_REPORT(REPORT_ERROR, -1, *fortran_int_size == 4, "Error happens when using C-Coupler for model coupling: the size of an integer value in C-Coupler FORTRAN APIs is not 4 types. Please verify the compiler flag for C-Coupler and then recompile C-Coupler, to force the usage of 4-byte integer.");
}


#ifdef LINK_WITHOUT_UNDERLINE
extern "C" void register_root_component
#else
extern "C" void register_root_component_
#endif
(MPI_Fint *f_comm, const char *comp_name, const char *local_comp_type, const char *annotation, int *comp_id, 
int *enabled_in_parent_coupling_gen, int *change_dir, const char *executable_name)
{
    int flag;
    MPI_Comm local_comm = MPI_COMM_NULL;
    int root_comp_id;
    int current_proc_global_id;
    char file_name[NAME_STR_SIZE];
    MPI_Comm cpp_comm;


    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "start to register the root component model");

    check_and_verify_name_format_of_string_for_API(-1, comp_name, API_ID_COMP_MGT_REG_COMP, "the root component", annotation);
    check_API_parameter_string_length(-1, API_ID_COMP_MGT_REG_COMP, CCPL_NAME_STR_LEN, comp_name, "comp_name", annotation);

    if (comp_comm_group_mgt_mgr != NULL) 
        EXECUTION_REPORT(REPORT_ERROR, -1, false, "Error happens when registering the root component (\"%s\") at the model code with the annotation \"%s\": the root compnent (\"%s\") has been registered before, at the model code with the annotation \"%s\". Please note that there must be only one root component model at each MPI process", comp_name, annotation, comp_comm_group_mgt_mgr->get_root_component_model()->get_comp_name(), comp_comm_group_mgt_mgr->get_annotation_start());
    MPI_Initialized(&flag);
    if (flag == 0) {
        EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Initialize MPI when registering the root component \"%s\"", comp_name);
        MPI_Init(NULL, NULL);
    }

    synchronize_comp_processes_for_API(-1, API_ID_COMP_MGT_REG_COMP, MPI_COMM_WORLD, "registering root component", annotation);

    comp_comm_group_mgt_mgr = new Comp_comm_group_mgt_mgr(executable_name);
    import_report_setting();
	
	cpp_comm = MPI_Comm_f2c(*f_comm);
    if (cpp_comm != MPI_COMM_NULL) {
        EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Before MPI_barrier at root component \"%s\" for synchronizing the processes of the component (the corresponding model code annotation is \"%s\").", comp_name, annotation);
        EXECUTION_REPORT(REPORT_ERROR,-1, MPI_Barrier(cpp_comm) == MPI_SUCCESS);
        EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "After MPI_barrier at root component \"%s\" for synchronizing the processes of the component (the corresponding model code annotation is \"%s\").", comp_name, annotation);
        
    }

    original_grid_mgr = new Original_grid_mgt();
    root_comp_id = comp_comm_group_mgt_mgr->register_component(comp_name, local_comp_type, local_comm, -1, (*enabled_in_parent_coupling_gen) == 1, *change_dir, annotation);

    if (cpp_comm != MPI_COMM_NULL) {
        int input_comm_size, new_comm_size;
        int *input_comm_process_ids, *new_comm_process_ids, *temp_array;
        int current_proc_global_id, current_proc_local_id;
        MPI_Comm new_comm;
        EXECUTION_REPORT(REPORT_ERROR,-1, MPI_Comm_size(cpp_comm, &input_comm_size) == MPI_SUCCESS);
        new_comm = comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(root_comp_id, "C-Coupler code in register_root_component for getting component management node");
        EXECUTION_REPORT(REPORT_ERROR,-1, MPI_Comm_size(new_comm, &new_comm_size) == MPI_SUCCESS);
        EXECUTION_REPORT(REPORT_ERROR,-1, input_comm_size == new_comm_size, "Error happens when calling the API \"CCPL_register_component\" to register the root component model \"%s\": its input communicator does not match the communicator generated (the number of processes of the two communicators are different (%d VS %d)). Please check the model code with the annotation \"%s\"", comp_name, input_comm_size, new_comm_size, annotation);
        EXECUTION_REPORT(REPORT_ERROR,-1, MPI_Comm_rank(MPI_COMM_WORLD, &current_proc_global_id) == MPI_SUCCESS);
        input_comm_process_ids = new int [input_comm_size];
        new_comm_process_ids = new int [new_comm_size];
        temp_array = new int [new_comm_size];
        EXECUTION_REPORT(REPORT_ERROR,-1, MPI_Allgather(&current_proc_global_id, 1, MPI_INT, input_comm_process_ids, 1, MPI_INT, cpp_comm) == MPI_SUCCESS);
        EXECUTION_REPORT(REPORT_ERROR,-1, MPI_Allgather(&current_proc_global_id, 1, MPI_INT, new_comm_process_ids, 1, MPI_INT, new_comm) == MPI_SUCCESS);
        do_quick_sort(input_comm_process_ids, temp_array, 0, input_comm_size-1);
        do_quick_sort(new_comm_process_ids, temp_array, 0, new_comm_size-1);
        for (int i = 0; i < input_comm_size; i ++)
            EXECUTION_REPORT(REPORT_ERROR,-1, input_comm_process_ids[i] == new_comm_process_ids[i], 
                             "Error happens when calling the API \"CCPL_register_component\" to register the root component model \"%s\": its input communicator does not match the communicator generated (processes of the two communicators are not the same). Please check the model code with the annotation \"%s\"",
                             comp_name, annotation);
        delete [] input_comm_process_ids;
        delete [] new_comm_process_ids;
        delete [] temp_array;
    }
    else *f_comm = MPI_Comm_c2f(local_comm);

    *comp_id = root_comp_id;

    sprintf(file_name, "%s/all/env_run.xml", comp_comm_group_mgt_mgr->get_config_root_dir());
    components_time_mgrs->define_root_comp_time_mgr(root_comp_id, file_name);
    fields_info = new Field_info_mgt();
    remapping_configuration_mgr->add_remapping_configuration(comp_comm_group_mgt_mgr->get_global_node_root()->get_comp_id());
    if (comp_comm_group_mgt_mgr->get_global_node_of_local_comp(root_comp_id, true, "")->is_real_component_model())
        remapping_configuration_mgr->add_remapping_configuration(root_comp_id);

    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Finish registering the root component model");
}


#ifdef LINK_WITHOUT_UNDERLINE
extern "C" void register_component
#else
extern "C" void register_component_
#endif
(int *parent_comp_id, const char *comp_name, const char *local_comp_type, MPI_Fint *f_comm, const char *annotation, int *enabled_in_parent_coupling_gen, int *change_dir, int *comp_id)
{
    MPI_Comm cpp_comm = MPI_Comm_f2c(*f_comm);


    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "start to register component model \%s\"", comp_name);

    check_and_verify_name_format_of_string_for_API(-1, comp_name, API_ID_COMP_MGT_REG_COMP, "the new component", annotation);
    EXECUTION_REPORT(REPORT_ERROR, -1, comp_comm_group_mgt_mgr->is_legal_local_comp_id(*parent_comp_id, true), "Error happens when calling the API \"CCPL_register_component\" to register a component model \"%s\": The parameter of \"parent_id\" (currently is 0x%x) is wrong (not the legal ID of a component). Please check the model code with the annotation \"%s\"", comp_name, *parent_comp_id, annotation);
    check_for_coupling_registration_stage(*parent_comp_id, API_ID_COMP_MGT_REG_COMP, false, annotation);
    check_API_parameter_string_length(*parent_comp_id, API_ID_COMP_MGT_REG_COMP, CCPL_NAME_STR_LEN, comp_name, "comp_name", annotation);

    if (cpp_comm != MPI_COMM_NULL) {
        synchronize_comp_processes_for_API(*parent_comp_id, API_ID_COMP_MGT_REG_COMP, cpp_comm, "registering a component based on the parent component", annotation);
        check_API_parameter_string(*parent_comp_id, API_ID_COMP_MGT_REG_COMP, cpp_comm, "registering a component based on an available communicator", comp_name, "comp_name", annotation);
    }
    else synchronize_comp_processes_for_API(*parent_comp_id, API_ID_COMP_MGT_REG_COMP, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(*parent_comp_id, "C-Coupler code for get comm group in register_component interface"), "registering component based on the parent component", annotation);

    *comp_id = comp_comm_group_mgt_mgr->register_component(comp_name, local_comp_type, cpp_comm, *parent_comp_id, (*enabled_in_parent_coupling_gen) == 1, *change_dir, annotation);
    if (comp_comm_group_mgt_mgr->get_global_node_of_local_comp(*comp_id,true,"")->is_real_component_model())
        remapping_configuration_mgr->add_remapping_configuration(*comp_id);
    components_time_mgrs->clone_parent_comp_time_mgr(*comp_id, *parent_comp_id, annotation);

    *f_comm = MPI_Comm_c2f(cpp_comm);

    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Finish registering component model \%s\"", comp_comm_group_mgt_mgr->get_global_node_of_local_comp(*comp_id,true,"")->get_full_name());
}


#ifdef LINK_WITHOUT_UNDERLINE
extern "C" void get_ccpl_comp_log_file_name
#else
extern "C" void get_ccpl_comp_log_file_name_
#endif
(int *comp_id, char *file_name, int *size_file_name, int *log_file_opened, const char *annotation)
{
    int log_file_device_id;

    
    check_for_component_registered(*comp_id, API_ID_COMP_MGT_GET_COMP_LOG_FILE_NAME, annotation, true);
    const char *log_file_name = comp_comm_group_mgt_mgr->get_comp_model_log_file(*comp_id, log_file_device_id);
    if (log_file_device_id == -1)
        *log_file_opened = 0;
    else *log_file_opened = 1;
    copy_out_string_to_Fortran_API(*comp_id, *size_file_name, file_name, log_file_name, API_ID_COMP_MGT_GET_COMP_LOG_FILE_NAME, "file_name", annotation);
}


#ifdef LINK_WITHOUT_UNDERLINE
extern "C" void get_ccpl_comp_log_file_device
#else
extern "C" void get_ccpl_comp_log_file_device_
#endif
(int *comp_id, int *log_file_device_id, int *log_file_opened, const char *annotation)
{
    check_for_component_registered(*comp_id, API_ID_COMP_MGT_GET_COMP_LOG_FILE_DEVICE, annotation, true);
    *log_file_opened = comp_comm_group_mgt_mgr->get_global_node_of_local_comp(*comp_id, true, "get_ccpl_comp_log_file_device_")->open_comp_model_log_file(log_file_device_id);
}


#ifdef LINK_WITHOUT_UNDERLINE
extern "C" void get_id_of_component
#else
extern "C" void get_id_of_component_
#endif
(const char *comp_name, const char *annotation, int *comp_id)
{
    check_and_verify_name_format_of_string_for_API(-1, comp_name, API_ID_COMP_MGT_GET_COMP_ID, "the component", annotation);
    check_for_component_registered(-1, API_ID_COMP_MGT_GET_COMP_ID, annotation, true);
    check_API_parameter_string_length(-1, API_ID_COMP_MGT_GET_COMP_ID, CCPL_NAME_STR_LEN, comp_name, "comp_name", annotation);

    Comp_comm_group_mgt_node *node = comp_comm_group_mgt_mgr->search_comp_with_comp_name(comp_name);
    if (comp_comm_group_mgt_mgr->does_comp_name_include_reserved_prefix(comp_name))
        node = NULL;

    if (node == NULL) {
        EXECUTION_REPORT(REPORT_ERROR, -1, false, "Error happens when calling the API \"CCPL_get_component_id\" to get the ID of a component model: no component model with the name of \"%s\" has been registered. Please check the model code at the annotation \"%s\"", comp_name, annotation);
        *comp_id = -1;
    }
    else *comp_id = node->get_local_node_id();    
}


#ifdef LINK_WITHOUT_UNDERLINE
extern "C" void check_is_comp_type_coupled
#else
extern "C" void check_is_comp_type_coupled_
#endif
(int *comp_id, const char *comp_type, int *is_coupled, const char *annotation)
{
    check_for_component_registered(*comp_id, API_ID_COMP_MGT_IS_COMP_TYPE_COUPLED, annotation, true);
    synchronize_comp_processes_for_API(*comp_id, API_ID_COMP_MGT_IS_COMP_TYPE_COUPLED, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(*comp_id, "check_is_comp_type_coupled_"), "synchorization for checking whether a type of component models is coupled", annotation);
    check_API_parameter_string(*comp_id, API_ID_COMP_MGT_IS_COMP_TYPE_COUPLED, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(*comp_id, "check_is_comp_type_coupled_"), "checking whether a type of component models have been registered to C-Coupler", comp_type, "comp_type", annotation);
    if (comp_comm_group_mgt_mgr->is_comp_type_coupled(*comp_id, comp_type, annotation))
        *is_coupled = 1; 
    else *is_coupled = 0;
}


#ifdef LINK_WITHOUT_UNDERLINE
extern "C" void is_current_process_in_component
#else
extern "C" void is_current_process_in_component_
#endif
(const char *comp_full_name, int *is_in_comp, const char *annotation)
{
    check_for_component_registered(-1, API_ID_COMP_MGT_IS_CURRENT_PROC_IN_COMP, annotation, true);
    check_API_parameter_string_length(-1, API_ID_COMP_MGT_IS_CURRENT_PROC_IN_COMP, 512, comp_full_name, "comp_full_name", annotation);
    Comp_comm_group_mgt_node *comp_node = comp_comm_group_mgt_mgr->search_global_node(comp_full_name);
    if (comp_node != NULL && comp_comm_group_mgt_mgr->does_comp_name_include_reserved_prefix(comp_node->get_comp_name()))
        comp_node = NULL;
    *is_in_comp = comp_node != NULL && comp_node->get_current_proc_local_id() != -1? 1 : 0;
}


#ifdef LINK_WITHOUT_UNDERLINE
extern "C" void get_current_proc_id_in_comp
#else
extern "C" void get_current_proc_id_in_comp_
#endif
(int *comp_id, int *proc_id, const char * annotation)
{
    check_for_component_registered(*comp_id, API_ID_COMP_MGT_GET_CURRENT_PROC_ID_IN_COMP, annotation, false);
    *proc_id = comp_comm_group_mgt_mgr->get_current_proc_id_in_comp(*comp_id, annotation);
}


#ifdef LINK_WITHOUT_UNDERLINE
extern "C" void get_num_proc_in_comp
#else
extern "C" void get_num_proc_in_comp_
#endif
(int *comp_id, int *num_proc, const char * annotation)
{
    check_for_component_registered(*comp_id, API_ID_COMP_MGT_GET_NUM_PROC_IN_COMP, annotation, false);
    *num_proc = comp_comm_group_mgt_mgr->get_num_proc_in_comp(*comp_id, annotation);
}


#ifdef LINK_WITHOUT_UNDERLINE
extern "C" void get_comp_proc_global_id
#else
extern "C" void get_comp_proc_global_id_
#endif
(int *comp_id, int *local_proc_id, int *global_proc_id, const char *annotation)
{
    check_for_component_registered(*comp_id, API_ID_COMP_MGT_GET_COMP_PROC_GLOBAL_ID, annotation, false);
    int num_proc = comp_comm_group_mgt_mgr->get_num_proc_in_comp(*comp_id, annotation);
    EXECUTION_REPORT(REPORT_ERROR, *comp_id, *local_proc_id >= 0 && *local_proc_id < num_proc, "Error happens when calling the API \"CCPL_get_component_process_global_id\": the parameter \"local_proc_id\" is wrong (its value (%d) is not between 0 and %d (maximum local process ID)). Please verify the model code corresponding to the annotation \"%s\"", *local_proc_id, num_proc-1, annotation);
    *global_proc_id = comp_comm_group_mgt_mgr->search_global_node(*comp_id)->get_local_proc_global_id(*local_proc_id);
}


#ifdef LINK_WITHOUT_UNDERLINE
extern "C" void ccpl_load_comps_full_names_from_config_file
#else
extern "C" void ccpl_load_comps_full_names_from_config_file_
#endif
(int *comp_id, const char *keyword, int *size_comps_full_names, int *size_individual_or_family, int *num_comps, const char *annotation)
{
    check_for_component_registered(*comp_id, API_ID_COUPLING_GEN_GET_COMPS, annotation, false);
    EXECUTION_REPORT_LOG(REPORT_LOG, *comp_id, true, "start to the full names of a set of component models from the corresponding configuration file.");
    synchronize_comp_processes_for_API(*comp_id, API_ID_COUPLING_GEN_GET_COMPS, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(*comp_id, "ccpl_load_omps_full_names_from_config_file_"), "first synchorization for getting the full names of a set of component models", annotation);
    check_API_parameter_string_length(*comp_id, API_ID_COUPLING_GEN_GET_COMPS, CCPL_NAME_STR_LEN, keyword, "keyword", annotation);
    coupling_generator->load_comps_full_names_from_config_file(*comp_id, keyword, *size_comps_full_names, *size_individual_or_family, num_comps, annotation);
    EXECUTION_REPORT_LOG(REPORT_LOG, *comp_id, true, "finish ccpl_load_comps_full_names_from_config_file");
}


#ifdef LINK_WITHOUT_UNDERLINE
extern "C" void ccpl_get_one_comp_full_name
#else
extern "C" void ccpl_get_one_comp_full_name_
#endif
(int *comp_id, const char *keyword, int *str_size, int *index, char *comp_full_name, int *local_individual_or_family, const char *annotation)
{
    EXECUTION_REPORT_LOG(REPORT_LOG, *comp_id, true, "At the beginning of cpl_get_one_comp_full_name_");
    coupling_generator->get_one_comp_full_name(*comp_id, keyword, *str_size, *index-1, comp_full_name, local_individual_or_family , annotation);
    EXECUTION_REPORT_LOG(REPORT_LOG, *comp_id, true, "At the end of cpl_get_one_comp_full_name_");
}


#ifdef LINK_WITHOUT_UNDERLINE
extern "C" void ccpl_finish_getting_configurable_comps_full_names
#else
extern "C" void ccpl_finish_getting_configurable_comps_full_names_
#endif
(int *comp_id, const char *annotation)
{
    EXECUTION_REPORT_LOG(REPORT_LOG, *comp_id, true, "At the beginning of ccpl_finish_getting_configurable_comps_full_names_");
    coupling_generator->clear();
    EXECUTION_REPORT_LOG(REPORT_LOG, *comp_id, true, "At the end of ccpl_finish_getting_configurable_comps_full_names_");
}


#ifdef LINK_WITHOUT_UNDERLINE
extern "C" void ccpl_family_coupling_generation
#else
extern "C" void ccpl_family_coupling_generation_
#endif
(int *comp_id, const char * annotation)
{
    check_for_component_registered(*comp_id, API_ID_COUPLING_GEN_FAMILY, annotation, false);
    EXECUTION_REPORT_LOG(REPORT_LOG, *comp_id, true, "start to generate coupling procedures for the component model \"%s\" and its descendants", comp_comm_group_mgt_mgr->get_global_node_of_local_comp(*comp_id,true,"")->get_full_name());
    synchronize_comp_processes_for_API(*comp_id, API_ID_COUPLING_GEN_FAMILY, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(*comp_id, "C-Coupler code in ccpl_family_coupling_generation_"), "first synchorization for coupling generation of a component", annotation);    
    coupling_generator->generate_coupling_procedures_internal(*comp_id, true, false, annotation);
    EXECUTION_REPORT_LOG(REPORT_LOG, *comp_id, true, "Finish generating coupling procedures for the component model \"%s\" and its descendants", comp_comm_group_mgt_mgr->get_global_node_of_local_comp(*comp_id,true,"")->get_full_name());
}


#ifdef LINK_WITHOUT_UNDERLINE
extern "C" void ccpl_individual_coupling_generation
#else
extern "C" void ccpl_individual_coupling_generation_
#endif
(int *comp_id, const char * annotation)
{
    check_for_component_registered(*comp_id, API_ID_COUPLING_GEN_INDIVIDUAL, annotation, false);
    EXECUTION_REPORT_LOG(REPORT_LOG, *comp_id, true, "start to generate coupling procedures for the component model \"%s\"", comp_comm_group_mgt_mgr->get_global_node_of_local_comp(*comp_id,true,"")->get_full_name());
    synchronize_comp_processes_for_API(*comp_id, API_ID_COUPLING_GEN_INDIVIDUAL, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(*comp_id, "C-Coupler code in ccpl_individual_coupling_generation_"), "first synchorization for coupling generation of a component", annotation);    
    coupling_generator->generate_coupling_procedures_internal(*comp_id, false, true, annotation);
    EXECUTION_REPORT_LOG(REPORT_LOG, *comp_id, true, "Finish generating coupling procedures for the component model \"%s\" and its descendants", comp_comm_group_mgt_mgr->get_global_node_of_local_comp(*comp_id,true,"")->get_full_name());
}


#ifdef LINK_WITHOUT_UNDERLINE
extern "C" void ccpl_begin_external_coupling_generation
#else
extern "C" void ccpl_begin_external_coupling_generation_
#endif
(int *num_comps, int *size_comps_full_names, int *size_individual_or_family, const char *annotation)
{
    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "start to generate coupling procedures for a set of specified component model");
    EXECUTION_REPORT(REPORT_ERROR, -1, *num_comps >= 1, "ERROR happens when calling the API \"CCPL_do_external_coupling_generation\": the value of \"num_comps\" cannot be smaller than 1. Please verify the model code with the annotation \"%s\"", annotation);
    EXECUTION_REPORT(REPORT_ERROR, -1, *num_comps <= *size_comps_full_names, "ERROR happens when calling the API \"CCPL_do_external_coupling_generation\": the value of \"num_comps\" cannot be bigger than the array size of \"comps_full_names\". Please verify the model code with the annotation \"%s\"", annotation); 
    EXECUTION_REPORT(REPORT_ERROR, -1, *num_comps <= *size_individual_or_family, "ERROR happens when calling the API \"CCPL_do_external_coupling_generation\": the value of \"num_comps\" cannot be bigger than the array size of \"individual_or_family\". Please verify the model code with the annotation \"%s\"", annotation); 
    coupling_generator->begin_external_coupling_generation();    
    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Finish generating coupling procedures for a set of specified component model");
}


#ifdef LINK_WITHOUT_UNDERLINE
extern "C" void ccpl_add_comp_for_external_coupling_generation
#else
extern "C" void ccpl_add_comp_for_external_coupling_generation_
#endif
(const char *comp_full_name, int *individual_or_family, const char *annotation)
{
    coupling_generator->add_comp_for_external_coupling_generation(comp_full_name, *individual_or_family, annotation);
}


#ifdef LINK_WITHOUT_UNDERLINE
extern "C" void ccpl_end_external_coupling_generation
#else
extern "C" void ccpl_end_external_coupling_generation_
#endif
(const char *annotation)
{
    coupling_generator->do_external_coupling_generation(API_ID_COUPLING_GEN_EXTERNAL, annotation);
}


#ifdef LINK_WITHOUT_UNDERLINE
extern "C" void ccpl_end_registration
#else
extern "C" void ccpl_end_registration_
#endif
(int *comp_id, const char * annotation)
{
    check_for_component_registered(*comp_id, API_ID_COMP_MGT_END_COMP_REG, annotation, false);
    
    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "start to end the coupling registration for the component model \"%s\"", comp_comm_group_mgt_mgr->get_global_node_of_local_comp(*comp_id,true,"")->get_full_name());    
    synchronize_comp_processes_for_API(*comp_id, API_ID_COMP_MGT_END_COMP_REG, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(*comp_id, "ccpl_end_registration_"), "first synchorization for ending the registration of a component", annotation);

    if (((*comp_id) & TYPE_ID_SUFFIX_MASK) == 1) {
        coupling_generator->do_overall_coupling_generation(comp_comm_group_mgt_mgr->get_global_node_of_local_comp(*comp_id, true, "in ccpl_end_registration_")->get_comp_full_name(), annotation);
         coupling_generator->generate_IO_procedures();
        delete all_H2D_remapping_wgt_files_info;
    }
    synchronize_comp_processes_for_API(*comp_id, API_ID_COMP_MGT_END_COMP_REG, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(*comp_id, "C-Coupler code in register_component for getting component management node"), "second synchorization for ending the registration of a component", annotation);

    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Finish ending the coupling registration for the component model \"%s\"", comp_comm_group_mgt_mgr->get_global_node_of_local_comp(*comp_id, true, "")->get_full_name());
    EXECUTION_REPORT(REPORT_PROGRESS, *comp_id, true, "The coupling registration stage of the component model \"%s\" is successfully ended at the model code with the annotation \"%s\"", comp_comm_group_mgt_mgr->get_global_node_of_local_comp(*comp_id,true,"")->get_full_name(), annotation);
}


#ifdef LINK_WITHOUT_UNDERLINE
extern "C" void register_v1d_grid_without_data
#else
extern "C" void register_v1d_grid_without_data_
#endif
(int *comp_id, int *grid_id, const char *grid_name, const char *coord_unit, int *grid_size, const char *annotation)
{
	char API_label[NAME_STR_SIZE];
	
	EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "start to register V1D grid %s", grid_name);
	common_checking_for_grid_registration(*comp_id, grid_name, coord_unit, API_ID_GRID_MGT_REG_V1D_GRID_NO_DATA, annotation);
    get_API_hint(*comp_id, API_ID_GRID_MGT_REG_V1D_GRID_NO_DATA, API_label);
    EXECUTION_REPORT(REPORT_ERROR, *comp_id, *grid_size > 0, "Error happens when calling the API \"%s\" to register a vertical grid \"%s\": the input paramter \"grid_size\" is not larger than 0. Please verify the model code with the annotation \"%s\"", API_label, grid_name, annotation);
    *grid_id = original_grid_mgr->register_V1D_grid_via_data(API_ID_GRID_MGT_REG_V1D_GRID_NO_DATA, *comp_id, grid_name, 4, coord_unit, *grid_size, 0.0, NULL, NULL, annotation);
	
	EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Finish registering V1D grid %s", grid_name);
}


#ifdef LINK_WITHOUT_UNDERLINE
extern "C" void register_v1d_grid_with_data
#else
extern "C" void register_v1d_grid_with_data_
#endif
(int *comp_id, int *grid_id, const char *grid_name, int *grid_type, const char *coord_unit, int *dim_size2,  
 int *dim_size3, const char *data_type, void *value1, void *value2, void *value3, const char *annotation)
{
    double temp_value1, *temp_value2, *temp_value3;
    int API_id;
    char API_label[NAME_STR_SIZE];


    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "start to register V1D grid %s", grid_name);
        
    switch (*grid_type) {
        case 1:
            API_id = API_ID_GRID_MGT_REG_V1D_Z_GRID_VIA_MODEL;
            break;    
        case 2:            
            API_id = API_ID_GRID_MGT_REG_V1D_SIGMA_GRID_VIA_MODEL;
            break;
        case 3:            
            API_id = API_ID_GRID_MGT_REG_V1D_HYBRID_GRID_VIA_MODEL;
            break;            
        default:
            EXECUTION_REPORT(REPORT_ERROR, -1, "Software error in register_V1D_grid_with_data: wrong caller_label");
            break;
    }

    common_checking_for_grid_registration(*comp_id, grid_name, coord_unit, API_id, annotation);

    get_API_hint(*comp_id, API_id, API_label);
    EXECUTION_REPORT(REPORT_ERROR, *comp_id, *dim_size2 != 0 && *dim_size3 != 0, "Error happens when calling the API \"%s\" to register a V1D grid \"%s\": some parameters of array have not be allocated. Please verify the model code with the annotation \"%s\" (please make sure all the arrays of grid data have been allocated)", API_label, grid_name, annotation);
    EXECUTION_REPORT(REPORT_ERROR, *comp_id, *dim_size2 > 1 && *dim_size3 == *dim_size2, "Error happens when calling the API \"%s\" to register a V1D grid \"%s\": the implicit grid size that is determined by the parameter arrays is wrong: the grid size (currently is %d) is smaller than 2 or the sizes of two paramenter arrays are different. Please verify the model code with the annotation \"%s\".", API_label, grid_name, dim_size2, annotation);
    EXECUTION_REPORT(REPORT_ERROR, -1, words_are_the_same(data_type, DATA_TYPE_FLOAT) || words_are_the_same(data_type, DATA_TYPE_DOUBLE), "Software error in register_V1D_grid_with_data: wrong data type");
    temp_value2 = new double [*dim_size2];
    temp_value3 = new double [*dim_size3];
    if (words_are_the_same(data_type, DATA_TYPE_FLOAT)) {
        transform_datatype_of_arrays((float*)value1, &temp_value1, 1);
        transform_datatype_of_arrays((float*)value2, temp_value2, *dim_size2);
        transform_datatype_of_arrays((float*)value3, temp_value3, *dim_size3);
    }
    else {
        transform_datatype_of_arrays((double*)value1, &temp_value1, 1);
        transform_datatype_of_arrays((double*)value2, temp_value2, *dim_size2);
        transform_datatype_of_arrays((double*)value3, temp_value3, *dim_size3);
    }

    EXECUTION_REPORT(REPORT_ERROR, *comp_id, is_array_in_sorting_order(temp_value2,*dim_size2) != 0, "Error happens when calling the API \"%s\" to register a V1D grid \"%s\": some arrays of parameters are not in a descending/ascending order. Please check the model code with the annotation \"%s\".", API_label, grid_name, annotation);
    *grid_id = original_grid_mgr->register_V1D_grid_via_data(API_id, *comp_id, grid_name, *grid_type, coord_unit, *dim_size2, temp_value1, temp_value2, temp_value3, annotation);

    delete [] temp_value2;
    delete [] temp_value3;

    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Finish registering V1D grid %s", grid_name);
}


#ifdef LINK_WITHOUT_UNDERLINE
extern "C" void set_3d_grid_3d_vertical_coord_field
#else
extern "C" void set_3d_grid_3d_vertical_coord_field_
#endif
(int *grid_id, int *field_id, const char *static_or_dynamic, const char *annotation)
{
    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Start to execute CCPL_set_3D_grid_3D_vertical_coord_field");
	check_for_component_registered(-1, API_ID_GRID_MGT_SET_3D_GRID_3D_VERT_FLD, annotation, true);
	original_grid_mgr->set_3D_grid_3D_vertical_coord_field_inst(*grid_id, *field_id, static_or_dynamic, annotation);
	EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Finish executing CCPL_set_3D_grid_3D_vertical_coord_field");
}


#ifdef LINK_WITHOUT_UNDERLINE
extern "C" void set_3d_grid_surface_field
#else
extern "C" void set_3d_grid_surface_field_
#endif
(int *grid_id, int *field_id, int *static_or_dynamic_or_external, const char *annotation)
{
    char API_label[NAME_STR_SIZE];
    int comp_id, API_id;


    if (*static_or_dynamic_or_external == BOTTOM_FIELD_VARIATION_STATIC)
        API_id = API_ID_GRID_MGT_SET_3D_GRID_STATIC_BOT_FLD;
    else if (*static_or_dynamic_or_external == BOTTOM_FIELD_VARIATION_DYNAMIC) 
        API_id = API_ID_GRID_MGT_SET_3D_GRID_DYN_BOT_FLD;
    else if (*static_or_dynamic_or_external == BOTTOM_FIELD_VARIATION_EXTERNAL) 
        API_id = API_ID_GRID_MGT_SET_3D_GRID_EXTERNAL_BOT_FLD;
    else EXECUTION_REPORT(REPORT_ERROR, -1, false, "software error in set_3d_grid_surface_field_: wrong value of static_or_dynamic_or_external");
    get_API_hint(-1, API_id, API_label);    
    check_for_component_registered(-1, API_id, annotation, true);
    EXECUTION_REPORT(REPORT_ERROR, -1, original_grid_mgr->is_grid_id_legal(*grid_id), "Error happens when calling the API \"%s\" to set the surface field of a 3-D grid: the parameter of \"grid_id\" is wrong. Please verify the model code with the annotation \"%s\".", API_label, annotation);
    
    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Start to set surface field for the 3D grid %s", original_grid_mgr->get_name_of_grid(*grid_id));
    
    comp_id = original_grid_mgr->get_comp_id_of_grid(*grid_id);
    if (*static_or_dynamic_or_external != BOTTOM_FIELD_VARIATION_EXTERNAL) {
        EXECUTION_REPORT(REPORT_ERROR, comp_id, memory_manager->check_is_legal_field_instance_id(*field_id), "Error happens when calling the API \"%s\" to set the surface field of a 3-D grid: the parameter of \"field_id\" is wrong. Please verify the model code with the annotation \"%s\".", API_label, annotation);
        EXECUTION_REPORT(REPORT_ERROR, comp_id, comp_id == memory_manager->get_field_instance(*field_id)->get_comp_id(), "Error happens when calling the API \"%s\" to set the surface field of a 3-D grid: the components corresponding to the parameters of \"grid_id\" and \"field_id\" are different. Please verify the model code with the annotation \"%s\".", API_label, annotation);
    }    
    check_for_coupling_registration_stage(comp_id, API_id, true, annotation);
    original_grid_mgr->set_3d_grid_bottom_field(comp_id, *grid_id, *field_id, *static_or_dynamic_or_external, API_id, API_label, annotation);
    if (*static_or_dynamic_or_external != BOTTOM_FIELD_VARIATION_EXTERNAL)
        comp_comm_group_mgt_mgr->search_global_node(comp_id)->get_restart_mgr()->add_restarted_field_instance(memory_manager->get_field_instance(*field_id), true);

    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Finish to setting surface field for the 3D grid %s", original_grid_mgr->get_name_of_grid(*grid_id));
}


#ifdef LINK_WITHOUT_UNDERLINE
extern "C" void register_md_grid_via_multi_grids
#else
extern "C" void register_md_grid_via_multi_grids_
#endif
(int *comp_id, int *grid_id, const char *grid_name, int *sub_grid1_id, int *sub_grid2_id, int *sub_grid3_id, int *size_mask, int *mask, const char *annotation)
{
    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Start to register an MD grid %s", grid_name);

    common_checking_for_grid_registration(*comp_id, grid_name, NULL, API_ID_GRID_MGT_REG_MD_GRID_VIA_MULTI_GRIDS, annotation);
    *grid_id = original_grid_mgr->register_md_grid_via_multi_grids(*comp_id, grid_name, *sub_grid1_id, *sub_grid2_id, *sub_grid3_id, *size_mask, mask, annotation);

    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Finish registering an MD grid %s", grid_name);
}


#ifdef LINK_WITHOUT_UNDERLINE
extern "C" void register_h2d_grid_with_global_data
#else
extern "C" void register_h2d_grid_with_global_data_
#endif
(int *comp_id, int *grid_id, const char *grid_name, const char *edge_type, const char *coord_unit, const char *cyclic_or_acyclic, const char *data_type, int *dim_size1, int *dim_size2, int *size_center_lon, int *size_center_lat, 
 int *size_mask, int *size_area, int *size_vertex_lon, int *size_vertex_lat, char *min_lon, char *max_lon, char *min_lat, char *max_lat, char *center_lon, char *center_lat, int *mask, char *area, char *vertex_lon, char *vertex_lat, const char *annotation)
{
    char tmp_min_lon[8], tmp_max_lon[8], tmp_min_lat[8], tmp_max_lat[8], tmp_cyclic_or_acyclic[NAME_STR_SIZE];
    int data_type_size = 4;


    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Start to register an H2D grid \"%s\" with global data", grid_name);

    common_checking_for_grid_registration(*comp_id, grid_name, coord_unit, API_ID_GRID_MGT_REG_H2D_GRID_VIA_GLOBAL_DATA, annotation);

    if (words_are_the_same(data_type, DATA_TYPE_DOUBLE))
        data_type_size = 8;
    memcpy(tmp_min_lon, min_lon, data_type_size);
    memcpy(tmp_max_lon, max_lon, data_type_size);
    memcpy(tmp_min_lat, min_lat, data_type_size);
    memcpy(tmp_max_lat, max_lat, data_type_size);
    strcpy(tmp_cyclic_or_acyclic, cyclic_or_acyclic);

    *grid_id = original_grid_mgr->register_H2D_grid_via_global_data(*comp_id, grid_name, edge_type, coord_unit, tmp_cyclic_or_acyclic, data_type, *dim_size1, *dim_size2, *size_center_lon, *size_center_lat, 
                                                                    *size_mask, *size_area, *size_vertex_lon, *size_vertex_lat, tmp_min_lon, tmp_max_lon, tmp_min_lat, tmp_max_lat, center_lon, center_lat, mask, area, vertex_lon, vertex_lat, annotation,
                                                                    API_ID_GRID_MGT_REG_H2D_GRID_VIA_GLOBAL_DATA);
    if (report_error_enabled) {
        char nc_file_name[NAME_STR_SIZE];
        sprintf(nc_file_name, "%s/%s@%s.nc", comp_comm_group_mgt_mgr->get_internal_H2D_grids_dir(), grid_name, comp_comm_group_mgt_mgr->get_global_node_of_local_comp(*comp_id, true, annotation)->get_full_name());
        char temp_grid_name[NAME_STR_SIZE];
        sprintf(temp_grid_name, "%s_temp", grid_name);
        original_grid_mgr->register_H2D_grid_via_file(*comp_id, temp_grid_name, nc_file_name, annotation);
    }

    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Finish registering an H2D grid %s", grid_name);
}


#ifdef LINK_WITHOUT_UNDERLINE
extern "C" void register_h2d_grid_with_local_data
#else
extern "C" void register_h2d_grid_with_local_data_
#endif
(int *comp_id, int *grid_id, const char *grid_name, const char *edge_type, const char *coord_unit, const char *cyclic_or_acyclic, const char *data_type, int *grid_size, int *num_local_cells, int *size_local_cells_global_index, int *size_center_lon, int *size_center_lat, 
 int *size_mask, int *size_area, int *size_vertex_lon, int *size_vertex_lat, int *local_cells_global_index, char *min_lon, char *max_lon, char *min_lat, char *max_lat, char *center_lon, char *center_lat, int *mask, char *area, char *vertex_lon, char *vertex_lat, const char *decomp_name, int *decomp_id, const char *annotation)
{
    char tmp_min_lon[8], tmp_max_lon[8], tmp_min_lat[8], tmp_max_lat[8], tmp_cyclic_or_acyclic[NAME_STR_SIZE];
    int data_type_size = 4;


    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Start to register an H2D grid %s", grid_name);

    common_checking_for_grid_registration(*comp_id, grid_name, coord_unit, API_ID_GRID_MGT_REG_H2D_GRID_VIA_LOCAL_DATA, annotation);
    check_API_parameter_string_length(-1, API_ID_GRID_MGT_REG_H2D_GRID_VIA_LOCAL_DATA, CCPL_NAME_STR_LEN, decomp_name, "decomp_name", annotation);

    if (words_are_the_same(data_type, DATA_TYPE_DOUBLE))
        data_type_size = 8;
    memcpy(tmp_min_lon, min_lon, data_type_size);
    memcpy(tmp_max_lon, max_lon, data_type_size);
    memcpy(tmp_min_lat, min_lat, data_type_size);
    memcpy(tmp_max_lat, max_lat, data_type_size);
    strcpy(tmp_cyclic_or_acyclic, cyclic_or_acyclic);

    *grid_id = original_grid_mgr->register_H2D_grid_via_local_data(*comp_id, grid_name, edge_type, coord_unit, tmp_cyclic_or_acyclic, data_type, *grid_size, *num_local_cells, *size_local_cells_global_index, *size_center_lon, *size_center_lat, *size_mask, *size_area, 
                                                                   *size_vertex_lon, *size_vertex_lat, local_cells_global_index, tmp_min_lon, tmp_max_lon, tmp_min_lat, tmp_max_lat, center_lon, center_lat, mask, area, vertex_lon, vertex_lat, decomp_name, decomp_id, annotation, API_ID_GRID_MGT_REG_H2D_GRID_VIA_LOCAL_DATA);

    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Finish registering an H2D grid %s", grid_name);
}


#ifdef LINK_WITHOUT_UNDERLINE
extern "C" void register_h2d_grid_with_file
#else
extern "C" void register_h2d_grid_with_file_
#endif
(int *comp_id, int *grid_id, const char *grid_name, const char *data_file_name, const char *annotation)
{
    char full_data_file_name[NAME_STR_SIZE];


    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Start to register an H2D grid %s", grid_name);

    common_checking_for_grid_registration(*comp_id, grid_name, NULL, API_ID_GRID_MGT_REG_H2D_GRID_VIA_FILE, annotation);
    check_API_parameter_string_length(*comp_id, API_ID_GRID_MGT_REG_H2D_GRID_VIA_FILE, 1000, data_file_name, "data_file_name", annotation);
    sprintf(full_data_file_name, "%s/grids_weights/%s", comp_comm_group_mgt_mgr->get_root_comp_config_dir(), data_file_name);
    *grid_id = original_grid_mgr->register_H2D_grid_via_file(*comp_id, grid_name, full_data_file_name, annotation);

    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Finish registering an H2D grid %s", grid_name);
}


#ifdef LINK_WITHOUT_UNDERLINE
extern "C" void register_h2d_grid_from_another_component
#else
extern "C" void register_h2d_grid_from_another_component_
#endif
(int *comp_id, int *grid_id, const char *grid_name, const char *annotation)
{
    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Start to register an H2D grid %s", grid_name);

    common_checking_for_grid_registration(*comp_id, grid_name, NULL, API_ID_GRID_MGT_REG_H2D_GRID_VIA_COMP, annotation);
    check_and_verify_name_format_of_string_for_API(*comp_id, grid_name, API_ID_GRID_MGT_REG_H2D_GRID_VIA_COMP, "the C-Coupler grid", annotation);
    *grid_id = original_grid_mgr->register_H2D_grid_via_comp(*comp_id, grid_name, annotation);

    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Finish registering an H2D grid %s", grid_name);
}


#ifdef LINK_WITHOUT_UNDERLINE
extern "C" void register_cor_defined_grid
#else
extern "C" void register_cor_defined_grid_
#endif
(int *comp_id, const char *CCPL_grid_name, const char *CoR_grid_name, const char *annotation, int *grid_id)
{
    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Start to register a CoR grid %s", CCPL_grid_name);

    common_checking_for_grid_registration(*comp_id, CCPL_grid_name, NULL, API_ID_GRID_MGT_REG_GRID_VIA_COR, annotation);
    check_and_verify_name_format_of_string_for_API(*comp_id, CoR_grid_name, API_ID_GRID_MGT_REG_GRID_VIA_COR, "the CoR grid", annotation);
    check_API_parameter_string(*comp_id, API_ID_GRID_MGT_REG_GRID_VIA_COR, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(*comp_id, "C-Coupler code in register_cor_defined_grid for getting component management node"), "registering a grid", CCPL_grid_name, "CCPL_grid_name", annotation);
    check_API_parameter_string(*comp_id, API_ID_GRID_MGT_REG_GRID_VIA_COR, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(*comp_id, "C-Coupler code in register_cor_defined_grid for getting component management node"), "registering a grid", CoR_grid_name, "CoR_grid_name", annotation);
    *grid_id = original_grid_mgr->get_CoR_defined_grid(*comp_id, CCPL_grid_name, CoR_grid_name, annotation);

    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Finish registering a CoR grid %s", CCPL_grid_name);
}


#ifdef LINK_WITHOUT_UNDERLINE
extern "C" void register_mid_point_grid
#else
extern "C" void register_mid_point_grid_
#endif
(int *level_3D_grid_id, int *mid_3D_grid_id, int *mid_1D_grid_id, int *size_mask, int *mask, const char *annotation)
{
    char API_label[NAME_STR_SIZE];


    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Start to register a middle level grid");

    get_API_hint(-1, API_ID_GRID_MGT_REG_MID_POINT_GRID, API_label);    
    check_for_ccpl_managers_allocated(API_ID_GRID_MGT_REG_MID_POINT_GRID, annotation);
    EXECUTION_REPORT(REPORT_ERROR, -1, original_grid_mgr->is_grid_id_legal(*level_3D_grid_id), "Error happens when calling the API \"%s\" to register the mid-point grid of a grid: the grid ID of the interface-level grid (level_3D_grid_id) is wrong. Please verify the model code with the annotation \"%s\".", API_label, annotation);
    check_for_coupling_registration_stage(original_grid_mgr->get_comp_id_of_grid(*level_3D_grid_id), API_ID_GRID_MGT_REG_MID_POINT_GRID, true, annotation);
    original_grid_mgr->register_mid_point_grid(*level_3D_grid_id, mid_3D_grid_id, mid_1D_grid_id, *size_mask, mask, annotation, API_label);

    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Start to register a middle level grid");
}


#ifdef LINK_WITHOUT_UNDERLINE
extern "C" void get_grid_size
#else
extern "C" void get_grid_size_
#endif
(int *grid_id, int *grid_size, const char *annotation)
{
    check_for_ccpl_managers_allocated(API_ID_GRID_MGT_GET_GRID_SIZE, annotation);

    *grid_size = original_grid_mgr->get_grid_size(*grid_id, annotation);
}


#ifdef LINK_WITHOUT_UNDERLINE
extern "C" void get_grid_id
#else
extern "C" void get_grid_id_
#endif
(int *comp_id, const char *grid_name, int *grid_id, const char *annotation)
{
    check_for_ccpl_managers_allocated(API_ID_GRID_MGT_GET_GRID_ID, annotation);
    check_API_parameter_string_length(*comp_id, API_ID_GRID_MGT_GET_GRID_ID, CCPL_NAME_STR_LEN, grid_name, "grid_name", annotation);
    *grid_id = original_grid_mgr->get_grid_id(*comp_id, grid_name, annotation);
}


#ifdef LINK_WITHOUT_UNDERLINE
extern "C" void get_h2d_grid_area_in_remapping_weights
#else
extern "C" void get_h2d_grid_area_in_remapping_weights_
#endif
(int *interface_id, int *field_index, void *output_area_data, int *area_array_size, const char *data_type, int *have_area, const char *annotation)
{
    check_for_ccpl_managers_allocated(API_ID_GRID_MGT_GET_H2D_GRID_AREA_FROM_WGTS, annotation);
    *have_area = inout_interface_mgr->get_h2d_grid_area_in_remapping_weights(*interface_id, (*field_index)-1, output_area_data, *area_array_size, data_type, annotation);
}


#ifdef LINK_WITHOUT_UNDERLINE
extern "C" void get_h2d_grid_data
#else
extern "C" void get_h2d_grid_data_
#endif
(int *grid_id, int *decomp_id, const char *label, const char *data_type, int *array_size, char *grid_data, const char *annotation)
{
    char API_label[NAME_STR_SIZE];


    check_for_ccpl_managers_allocated(API_ID_GRID_MGT_GET_H2D_GRID_DATA, annotation);
    get_API_hint(-1, API_ID_GRID_MGT_GET_H2D_GRID_DATA, API_label);
    EXECUTION_REPORT(REPORT_ERROR, -1, original_grid_mgr->is_grid_id_legal(*grid_id), "Error happens when calling the API \"%s\" to get the grid data of an H2D grid: the parameter of \"grid_id\" is wrong. Please verify the model code with the annotation \"%s\".", API_label, annotation);
    EXECUTION_REPORT(REPORT_ERROR, original_grid_mgr->get_comp_id_of_grid(*grid_id), original_grid_mgr->get_original_grid(*grid_id)->is_H2D_grid(), "Error happens when calling the API \"%s\" to get the grid data of an H2D grid: the grid \"%s\" is not an H2D grid. Please verify the model code with the annotation \"%s\".", API_label, original_grid_mgr->get_original_grid(*grid_id)->get_grid_name(), annotation);
    EXECUTION_REPORT(REPORT_ERROR, original_grid_mgr->get_comp_id_of_grid(*grid_id), *decomp_id == -1 || decomps_info_mgr->is_decomp_id_legal(*decomp_id), "Error happens when calling the API \"%s\" to get the grid data of an H2D grid: the decomp_id is wrong (must be -1 or a legal decomp_id). Please verify the model code with the annotation \"%s\".", API_label, annotation);
    if (*decomp_id != -1)
        EXECUTION_REPORT(REPORT_ERROR, original_grid_mgr->get_comp_id_of_grid(*grid_id), original_grid_mgr->get_comp_id_of_grid(*grid_id) == decomps_info_mgr->get_decomp_info(*decomp_id)->get_comp_id(), "Error happens when calling the API \"%s\" to get the grid data of an H2D grid: the grid_id and decomp_id do not correspond to the same component model. Please verify the model code with the annotation \"%s\".", API_label, annotation);
    original_grid_mgr->get_original_grid(*grid_id)->get_grid_data(*decomp_id, label, data_type, *array_size, grid_data, annotation, API_label);
}


#ifdef LINK_WITHOUT_UNDERLINE
extern "C" void register_parallel_decomposition
#else
extern "C" void register_parallel_decomposition_
#endif
(int *decomp_id, int *grid_id, int *num_local_cells, int *array_size, const int *local_cells_global_indx, const char *decomp_name, const char *annotation)
{
    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Start to register a parallel decomp %s", decomp_name);
    
    check_for_ccpl_managers_allocated(API_ID_DECOMP_MGT_REG_DECOMP, annotation);
    EXECUTION_REPORT(REPORT_ERROR, -1, original_grid_mgr->is_grid_id_legal(*grid_id), "Error happens when calling the API \"CCPL_register_normal_parallel_decomp\" to register a parallel decomposition \"%s\": the parameter \"grid_id\" is wrong. Please check the model code with the annotation \"%s\"", decomp_name, annotation);
    int comp_id = original_grid_mgr->get_comp_id_of_grid(*grid_id);
    check_API_parameter_string_length(comp_id, API_ID_DECOMP_MGT_REG_DECOMP, CCPL_NAME_STR_LEN, decomp_name, "decomp_name", annotation);
    check_for_coupling_registration_stage(comp_id, API_ID_DECOMP_MGT_REG_DECOMP, true, annotation);
    check_and_verify_name_format_of_string_for_API(comp_id, decomp_name, API_ID_DECOMP_MGT_REG_DECOMP, "the parallel decomposition", annotation);

    EXECUTION_REPORT(REPORT_ERROR, comp_id, original_grid_mgr->get_original_grid(*grid_id)->is_H2D_grid(), "Error happens when calling the API \"CCPL_register_normal_parallel_decomp\" to register a parallel decomposition \"%s\": the grid \"%s\" corresponding to the parameter \"grid_id\" is not a horizontal grid. Please check the model code with the annotation \"%s\"", decomp_name, original_grid_mgr->get_original_grid(*grid_id)->get_grid_name(), annotation);
    EXECUTION_REPORT(REPORT_ERROR, comp_id, *num_local_cells >= 0, "Error happens when calling the API \"CCPL_register_normal_parallel_decomp\" to register a parallel decomposition \"%s\": the parameter \"num_local_cells\" (currently is %d) cannot be smaller than 0. Please check the model code with the annotation \"%s\"", decomp_name, *num_local_cells, annotation);
    EXECUTION_REPORT(REPORT_ERROR, comp_id, *num_local_cells <= *array_size, "Error happens when calling the API \"CCPL_register_normal_parallel_decomp\" to register a parallel decomposition \"%s\": the array size (currently is %d) of the parameter \"local_cells_global_index\" cannot be smaller than the parameter \"num_local_cells\" (currently is %d). Please check the model code with the annotation \"%s\"", decomp_name, *num_local_cells, *array_size, annotation);
    int grid_size = original_grid_mgr->get_original_grid(*grid_id)->get_original_CoR_grid()->get_grid_size();
    for (int i = 0; i < *num_local_cells; i ++)
        if (local_cells_global_indx[i] != CCPL_NULL_INT)
            EXECUTION_REPORT(REPORT_ERROR, comp_id, local_cells_global_indx[i] > 0 && local_cells_global_indx[i] <= grid_size, "Error happens when calling the API \"CCPL_register_parallel_decomp\" to register a parallel decomposition \"%s\": some values (for example %d) in parameter \"local_cells_global_indx\" are not between 1 and the size of the grid (currently is %d). Please check the model code with the annotation \"%s\"", decomp_name, local_cells_global_indx[i], grid_size, annotation);
        
    *decomp_id = decomps_info_mgr->register_H2D_parallel_decomposition(decomp_name, *grid_id, *num_local_cells, local_cells_global_indx, annotation);

    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Finish registering a parallel decomp %s", decomp_name);
}


#ifdef LINK_WITHOUT_UNDERLINE
extern "C" void register_external_field_instance
#else
extern "C" void register_external_field_instance_
#endif
(int *field_instance_id, const char *field_name, long *data_buffer_ptr, int *field_size, int *decomp_id, int *comp_or_grid_id, 
 int *buf_mark, int *usage_tag, const char *unit, const char *data_type, const char *annotation)
{
    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Start to register a field instance %s", field_name);

    check_for_ccpl_managers_allocated(API_ID_FIELD_MGT_REG_FIELD_INST, annotation);
    *field_instance_id = memory_manager->register_external_field_instance(field_name, (void*)(*data_buffer_ptr), *field_size, *decomp_id, *comp_or_grid_id, *buf_mark, *usage_tag, unit, data_type, annotation);
    Field_mem_info *field_instance = memory_manager->get_field_instance(*field_instance_id);
    if (field_instance->is_REST_field_inst())
        comp_comm_group_mgt_mgr->search_global_node(field_instance->get_host_comp_id())->get_restart_mgr()->add_restarted_field_instance(field_instance, false);
    
    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Finish registering a field instance %s", field_name);
}


#ifdef LINK_WITHOUT_UNDERLINE
extern "C" void register_an_io_field_from_field_instance
#else
extern "C" void register_an_io_field_from_field_instance_
#endif
(int *field_inst_id, const char *field_IO_name, const char *annotation)
{
    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Start to register an I/O field %s", field_IO_name);

    check_for_ccpl_managers_allocated(API_ID_FIELD_MGT_REG_IO_FIELD_from_INST, annotation);
    IO_fields_mgr->register_IO_field(*field_inst_id, field_IO_name, annotation);

    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Finish registering an I/O field %s", field_IO_name);
}


#ifdef LINK_WITHOUT_UNDERLINE
extern "C" void register_io_fields_from_field_instances
#else
extern "C" void register_io_fields_from_field_instances_
#endif
(int *num_field_inst, int *size_field_inst_ids, int *field_inst_ids, const char *annotation)
{
    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Start to register I/O fields");

    check_for_ccpl_managers_allocated(API_ID_FIELD_MGT_REG_IO_FIELDs_from_INSTs, annotation);
    EXECUTION_REPORT(REPORT_ERROR, -1, *num_field_inst > 0, "Error happers when calling the API \"CCPL_register_IO_fields_from_field_instances\": the parameter \"num_field_inst\" must be larger than 0. Please check the model code with the annotation \"%s\".", annotation);
    EXECUTION_REPORT(REPORT_ERROR, -1, *size_field_inst_ids > 0, "Error happers when calling the API \"CCPL_register_IO_fields_from_field_instances\": the parameter \"field_inst_ids\" seems to be an empty array, which means it has not been allocated. Please check the model code with the annotation \"%s\".", annotation);
    IO_fields_mgr->register_IO_fields(*num_field_inst, *size_field_inst_ids, field_inst_ids, annotation);

    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Finish registering I/O fields");
}


#ifdef LINK_WITHOUT_UNDERLINE
extern "C" void register_a_new_io_field
#else
extern "C" void register_a_new_io_field_
#endif
(int *comp_or_grid_id, int *decomp_id, int *field_size, long *data_buffer_ptr, const char *field_IO_name, 
 const char *long_name, const char *unit, const char *data_type, const char * annotation)
{
    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Start to register an I/O field %s", field_IO_name);

    check_for_ccpl_managers_allocated(API_ID_FIELD_MGT_REG_IO_FIELD_from_BUFFER, annotation);
    IO_fields_mgr->register_IO_field(*comp_or_grid_id, *decomp_id, *field_size, (void*)(*data_buffer_ptr), field_IO_name, long_name, unit, data_type, annotation);
    
    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Finish registering an I/O field %s", field_IO_name);
}


#ifdef LINK_WITHOUT_UNDERLINE
extern "C" void define_single_timer
#else
extern "C" void define_single_timer_
#endif
(int *comp_id, int *timer_id, const char *freq_unit, int *freq_count, int *local_lag_count, int *remote_lag_count, const char *annotation)
{    
    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Start to define a timer");
    
    check_for_coupling_registration_stage(*comp_id, API_ID_TIME_MGT_DEFINE_SINGLE_TIMER, true, annotation);
    EXECUTION_REPORT(REPORT_ERROR, *comp_id, components_time_mgrs->get_time_mgr(*comp_id)->get_time_step_in_second() > 0, "Error happers when calling the API \"CCPL_define_single_timer\": the time step of the corresponding component model \"%s\" has not been set yet. Please specify the time step before defining a timer at the model code with the annotation \"%s\"", 
                     comp_comm_group_mgt_mgr->get_global_node_of_local_comp(*comp_id, true, annotation)->get_comp_full_name(), annotation);
    *timer_id = timer_mgr->define_timer(*comp_id, freq_unit, *freq_count, *local_lag_count, *remote_lag_count, annotation);
    
    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Finish defining a timer");
}


#ifdef LINK_WITHOUT_UNDERLINE
extern "C" void define_complex_timer
#else
extern "C" void define_complex_timer_
#endif
(int *comp_id, int *timer_id, int *children_timers_id, int *num_children_timers, int *array_size, int *or_or_and, const char *annotation)
{
    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Start to define a timer");

    check_for_coupling_registration_stage(*comp_id, API_ID_TIME_MGT_DEFINE_COMPLEX_TIMER, true, annotation);
    EXECUTION_REPORT(REPORT_ERROR, *comp_id, components_time_mgrs->get_time_mgr(*comp_id)->get_time_step_in_second() > 0, "The time step of the component \%s\" has not been set yet. Please specify the time step before defining a timer at the model code with the annotation \"%s\"", 
                     comp_comm_group_mgt_mgr->get_global_node_of_local_comp(*comp_id, true, annotation)->get_comp_name(), annotation);
    *timer_id = timer_mgr->define_timer(*comp_id, children_timers_id, *num_children_timers, *array_size, *or_or_and, annotation);

    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Finish defining a timer");
}


#ifdef LINK_WITHOUT_UNDERLINE
extern "C" void set_component_time_step
#else
extern "C" void set_component_time_step_
#endif
(int *comp_id, int *time_step_in_second, const char *annotation)
{
    check_for_coupling_registration_stage(*comp_id, API_ID_TIME_MGT_SET_NORMAL_TIME_STEP, true, annotation);

    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Start to set the time step of component model \%s\"", comp_comm_group_mgt_mgr->get_global_node_of_local_comp(*comp_id,true,"")->get_full_name());
        
    synchronize_comp_processes_for_API(*comp_id, API_ID_TIME_MGT_SET_NORMAL_TIME_STEP, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(*comp_id, "C-Coupler code in set_component_time_step_"), "setting the time step of a component model", annotation);
    check_API_parameter_int(*comp_id, API_ID_TIME_MGT_SET_NORMAL_TIME_STEP, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(*comp_id,"C-Coupler code in set_component_time_step_"), NULL, *time_step_in_second, "time step (the unit is seconds)", annotation);
    components_time_mgrs->set_component_time_step(*comp_id, *time_step_in_second, annotation);
    
    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Finish setting the time step of component model \%s\"", comp_comm_group_mgt_mgr->get_global_node_of_local_comp(*comp_id,true,"")->get_full_name());
}


#ifdef LINK_WITHOUT_UNDERLINE
extern "C" void reset_component_current_time_to_start_time
#else
extern "C" void reset_component_current_time_to_start_time_
#endif
(int *comp_id, const char *annotation)
{
    check_for_coupling_registration_stage(*comp_id, API_ID_TIME_MGT_RESET_TIME_TO_START, true, annotation);

    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Start to reset the current time of the component model \%s\" to start time", comp_comm_group_mgt_mgr->get_global_node_of_local_comp(*comp_id,true,"")->get_full_name());
    synchronize_comp_processes_for_API(*comp_id, API_ID_TIME_MGT_SET_NORMAL_TIME_STEP, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(*comp_id, "C-Coupler code in CCPL_reset_current_time_to_start_time_"), "resetting the current time of a component model to the initial time", annotation);
    components_time_mgrs->get_time_mgr(*comp_id)->reset_current_time_to_start_time(annotation);
    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Finish resetting the current time of the component model \%s\" to start time", comp_comm_group_mgt_mgr->get_global_node_of_local_comp(*comp_id,true,"")->get_full_name());    
}


#ifdef LINK_WITHOUT_UNDERLINE
extern "C" void advance_component_time
#else
extern "C" void advance_component_time_
#endif
(int *comp_id, const char *annotation)
{
    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Start to advance time");
    
    check_for_component_registered(*comp_id, API_ID_TIME_MGT_ADVANCE_TIME, annotation, false);
//    components_IO_output_procedures_mgr->get_component_IO_output_procedures(*comp_id)->execute();
    components_time_mgrs->advance_component_time(*comp_id, annotation);
    EXECUTION_REPORT(REPORT_PROGRESS, *comp_id, true, "Component model \"%s\" advance time at the model code with the annotation \"%s\"", comp_comm_group_mgt_mgr->get_global_node_of_local_comp(*comp_id,true,"")->get_full_name(), annotation);
    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Finish advancing time");
}


#ifdef LINK_WITHOUT_UNDERLINE
extern "C" void ccpl_write_restart
#else
extern "C" void ccpl_write_restart_
#endif
(int *comp_id, int *bypass_timer, int *bypass_imported_fields, const char *annotation)
{
    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Start to do restart write");
    check_for_component_registered(*comp_id, API_ID_RESTART_MGT_WRITE_IO, annotation, false);
    if (comp_comm_group_mgt_mgr->get_global_node_of_local_comp(*comp_id,true,annotation)->is_real_component_model())
        comp_comm_group_mgt_mgr->get_global_node_of_local_comp(*comp_id,true,annotation)->get_restart_mgr()->do_restart_write(annotation, *bypass_timer == 1, *bypass_imported_fields == 1);
    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Finish doing restart write");
}


#ifdef LINK_WITHOUT_UNDERLINE
extern "C" void ccpl_read_restart
#else
extern "C" void ccpl_read_restart_
#endif
(int *comp_id, const char *specified_file_name, const char *annotation)
{
    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Start to do restart read");
    check_for_component_registered(*comp_id, API_ID_RESTART_MGT_START_READ_IO, annotation, false);
    if (comp_comm_group_mgt_mgr->get_global_node_of_local_comp(*comp_id,true,annotation)->is_real_component_model())
        comp_comm_group_mgt_mgr->get_global_node_of_local_comp(*comp_id,true,annotation)->get_restart_mgr()->read_restart_mgt_info(specified_file_name, annotation);
    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Finish doing restart read");
}


#ifdef LINK_WITHOUT_UNDERLINE
extern "C" void ccpl_read_all_restart_fields
#else
extern "C" void ccpl_read_all_restart_fields_
#endif
(int *comp_id, const char *annotation)
{
    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Start to read all restart fields from restart data file");
    check_for_component_registered(*comp_id, API_ID_RESTART_MGT_READ_ALL, annotation, false);
    comp_comm_group_mgt_mgr->search_global_node(*comp_id)->get_restart_mgr()->read_all_restarted_fields(annotation);
    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Finish reading all restart fields from restart data file");
}


#ifdef LINK_WITHOUT_UNDERLINE
extern "C" void ccpl_read_import_interface_restart_fields
#else
extern "C" void ccpl_read_import_interface_restart_fields_
#endif
(int *interface_id, const char *annotation)
{
    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Start to all fields of an import interface from restart data file");
    Inout_interface *inout_interface = inout_interface_mgr->get_interface(*interface_id);
    EXECUTION_REPORT(REPORT_ERROR, -1, inout_interface != NULL, "Error happens when calling the API \"CCPL_restart_read_fields_interface\": the parameter of \"interface_id\" is wrong");
    inout_interface->read_restart_fields(API_ID_RESTART_MGT_READ_INTERFACE, annotation);
    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Finish reading all fields of an import interface from restart data file");
}


#ifdef LINK_WITHOUT_UNDERLINE
extern "C" void get_ccpl_restart_time
#else
extern "C" void get_ccpl_restart_time_
#endif
(int *comp_id, int *restart_date, int *restart_second, const char *annotation)
{
    check_for_component_registered(*comp_id, API_ID_RESTART_MGT_GET_SETTING, annotation, false);
    Time_mgt *time_mgr = components_time_mgrs->get_time_mgr(*comp_id);
    if (time_mgr->get_runtype_mark() == RUNTYPE_MARK_CONTINUE) {
        *restart_date = time_mgr->get_common_restart_full_time()/100000;
        *restart_second = time_mgr->get_common_restart_full_time()%100000;
    }
    else if (time_mgr->get_runtype_mark() == RUNTYPE_MARK_INITIAL) {
        *restart_date = -1;
        *restart_second = -1;
    }
    else {
        *restart_date = time_mgr->get_rest_refdate();
        *restart_second = time_mgr->get_rest_refsecond();
    }
}


#ifdef LINK_WITHOUT_UNDERLINE
    extern "C" void get_ccpl_original_case_name
#else
    extern "C" void get_ccpl_original_case_name_
#endif
(int *comp_id, int *size_original_case_name, char *original_case_name, const char *annotation)
{
    check_for_component_registered(*comp_id, API_ID_RESTART_MGT_GET_SETTING, annotation, false);
    Time_mgt *time_mgr = components_time_mgrs->get_time_mgr(*comp_id);
    if (time_mgr->get_runtype_mark() == RUNTYPE_MARK_CONTINUE || time_mgr->get_runtype_mark() == RUNTYPE_MARK_INITIAL)
        copy_out_string_to_Fortran_API(*comp_id, *size_original_case_name, original_case_name, time_mgr->get_case_name(), API_ID_RESTART_MGT_GET_SETTING, "original_case_name", annotation);
    else copy_out_string_to_Fortran_API(*comp_id, *size_original_case_name, original_case_name, time_mgr->get_rest_refcase(), API_ID_RESTART_MGT_GET_SETTING, "original_case_name", annotation);
}


#ifdef LINK_WITHOUT_UNDERLINE
    extern "C" void get_ccpl_run_type
#else
    extern "C" void get_ccpl_run_type_
#endif
(int *comp_id, int *size_run_type, char *run_type, const char *annotation)
{
    check_for_component_registered(*comp_id, API_ID_RESTART_MGT_GET_SETTING, annotation, false);
    Time_mgt *time_mgr = components_time_mgrs->get_time_mgr(*comp_id);
    copy_out_string_to_Fortran_API(*comp_id, *size_run_type, run_type, time_mgr->get_run_type(), API_ID_RESTART_MGT_GET_SETTING, "run_type", annotation);
}


#ifdef LINK_WITHOUT_UNDERLINE
extern "C" void is_restart_timer_on
#else
extern "C" void is_restart_timer_on_
#endif
(int *comp_id, int *check_result, const char *annotation)
{
    check_for_component_registered(*comp_id, API_ID_RESTART_MGT_IS_TIMER_ON, annotation, false);
    EXECUTION_REPORT(REPORT_ERROR, *comp_id, comp_comm_group_mgt_mgr->get_global_node_of_local_comp(*comp_id,true,"in is_restart_timer_on_")->is_real_component_model(), "Error happens when calling the API CCPL_is_restart_timer_on: the given component model \"%s\" is not a real model. Please verify the model code related to the annotation \"%s\"",
                     comp_comm_group_mgt_mgr->get_global_node_of_local_comp(*comp_id,true,"in is_restart_timer_on_")->get_comp_full_name(), annotation);
    if (components_time_mgrs->get_time_mgr(*comp_id)->is_restart_timer_on())
        *check_result = 1;
    else *check_result = 0;
}


#ifdef LINK_WITHOUT_UNDERLINE
extern "C" void check_ccpl_component_current_time
#else
extern "C" void check_ccpl_component_current_time_
#endif
(int *comp_id, int *date, int *second, const char *annotation)
{
    check_for_component_registered(*comp_id, API_ID_TIME_MGT_CHECK_CURRENT_TIME, annotation, false);
    components_time_mgrs->check_component_current_time(*comp_id, *date, *second, annotation);
}


#ifdef LINK_WITHOUT_UNDERLINE
extern "C" void is_ccpl_timer_on
#else
extern "C" void is_ccpl_timer_on_
#endif
(int *timer_id, int *is_on, const char *annotation)
{
    check_for_ccpl_managers_allocated(API_ID_TIME_MGT_IS_TIMER_ON, annotation);
    if (timer_mgr->is_timer_on(*timer_id, annotation))
        *is_on = 1;
    else *is_on = 0;
}


#ifdef LINK_WITHOUT_UNDERLINE
extern "C" void check_is_ccpl_model_last_step
#else
extern "C" void check_is_ccpl_model_last_step_
#endif
(int *comp_id, int *is_last_step, const char *annotation)
{
    check_for_component_registered(*comp_id, API_ID_TIME_MGT_IS_MODEL_LAST_STEP, annotation, false);
    Time_mgt *time_mgr = components_time_mgrs->get_time_mgr(*comp_id);
    if (time_mgr->get_current_num_time_step() == time_mgr->get_num_total_step() - 1)
        *is_last_step = 1;
    else *is_last_step = 0;
}


#ifdef LINK_WITHOUT_UNDERLINE
extern "C" void check_is_ccpl_model_run_ended
#else
extern "C" void check_is_ccpl_model_run_ended_
#endif
(int *comp_id, int *is_ended, const char *annotation)
{
    check_for_component_registered(*comp_id, API_ID_TIME_MGT_IS_MODEL_RUN_ENDED, annotation, false);
    if (components_time_mgrs->is_model_run_ended(*comp_id, annotation))
        *is_ended = 1;
    else *is_ended = 0;
}


#ifdef LINK_WITHOUT_UNDERLINE
extern "C" void register_normal_remap_interface
#else
extern "C" void register_normal_remap_interface_
#endif
(const char *interface_name, int *interface_id, int *num_fields, int *field_ids_src, int *field_ids_dst, int *timer_id, int *inst_or_aver, int *array_size1, int *array_size2, const char *annotation)
{
    char API_label[NAME_STR_SIZE];

    
    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Start to register remap interface");
    check_for_ccpl_managers_allocated(API_ID_INTERFACE_REG_NORMAL_REMAP, annotation);    
    check_API_parameter_string_length(-1, API_ID_INTERFACE_REG_NORMAL_REMAP, CCPL_NAME_STR_LEN, interface_name, "interface_name", annotation);
    get_API_hint(-1, API_ID_INTERFACE_REG_NORMAL_REMAP, API_label);
    *interface_id = inout_interface_mgr->register_normal_remap_interface(interface_name, *num_fields, field_ids_src, field_ids_dst, *timer_id, *inst_or_aver, *array_size1, *array_size2, API_label, annotation);

    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Finish registering remap interface");
}


#ifdef LINK_WITHOUT_UNDERLINE
extern "C" void register_frac_based_remap_interface
#else
extern "C" void register_frac_based_remap_interface_
#endif
(const char *interface_name, int *interface_id, int *num_fields, int *field_ids_src, int *field_ids_dst, int *timer_id, int *inst_or_aver, int *array_size1, int *array_size2, long *frac_src, long *frac_dst, int *size_frac_src, int *size_frac_dst, const char *frac_data_type, const char *annotation)
{
    char API_label[NAME_STR_SIZE];
    

    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Start to register fraction based remap interface");

    check_for_ccpl_managers_allocated(API_ID_INTERFACE_REG_FRAC_REMAP, annotation);        
    check_API_parameter_string_length(-1, API_ID_INTERFACE_REG_FRAC_REMAP, CCPL_NAME_STR_LEN, interface_name, "interface_name", annotation);
    get_API_hint(-1, API_ID_INTERFACE_REG_FRAC_REMAP, API_label);
    *interface_id = inout_interface_mgr->register_frac_based_remap_interface(interface_name, *num_fields, field_ids_src, field_ids_dst, *timer_id, *inst_or_aver, *array_size1, *array_size2, (void*)(*frac_src), (void*)(*frac_dst), *size_frac_src, *size_frac_dst, frac_data_type, API_label, annotation);

    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Finish registering fraction based remap interface");
}


#ifdef LINK_WITHOUT_UNDERLINE
extern "C" void register_inout_interface
#else
extern "C" void register_inout_interface_
#endif
(const char *interface_name, int *interface_id, int *import_or_export, int *num_fields, int *field_ids, int *timer_id, int *inst_or_aver, const char *annotation, int *array_size1)
{
    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Start to register an import/export interface");

    if (*import_or_export == 0) {
        check_for_ccpl_managers_allocated(API_ID_INTERFACE_REG_IMPORT, annotation);
        check_API_parameter_string_length(-1, API_ID_INTERFACE_REG_IMPORT, CCPL_NAME_STR_LEN, interface_name, "interface_name", annotation);
        *interface_id = inout_interface_mgr->register_inout_interface(interface_name, *import_or_export, *num_fields, field_ids, *array_size1, *timer_id, *inst_or_aver, annotation, INTERFACE_SOURCE_REGISTER);
    }
    else {
        check_for_ccpl_managers_allocated(API_ID_INTERFACE_REG_EXPORT, annotation);
        check_API_parameter_string_length(-1, API_ID_INTERFACE_REG_EXPORT, CCPL_NAME_STR_LEN, interface_name, "interface_name", annotation);
        *interface_id = inout_interface_mgr->register_inout_interface(interface_name, *import_or_export, *num_fields, field_ids, *array_size1, *timer_id, 0, annotation, INTERFACE_SOURCE_REGISTER);
    }    

    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Finish registering an import/export interface");
}


#ifdef LINK_WITHOUT_UNDERLINE
extern "C" void check_is_ccpl_import_field_connected
#else
extern "C" void check_is_ccpl_import_field_connected_
#endif
(int *interface_id, int *field_instance_id, int *check_result, const char *annotation)
{
    check_for_ccpl_managers_allocated(API_ID_INTERFACE_CHECK_IMPORT_FIELD_CONNECTED, annotation);
    Inout_interface *import_interface = inout_interface_mgr->get_interface(*interface_id);
    EXECUTION_REPORT(REPORT_ERROR, -1, import_interface != NULL, "ERROR happens when calling the API \"CCPL_check_is_import_field_connected\": the parameter \"interface_id\" is not a legal ID of a coupling interface. Please verify the model code with the annotation \"%s\".", annotation);
    *check_result = import_interface->check_is_import_field_connected(*field_instance_id, annotation);
}


#ifdef LINK_WITHOUT_UNDERLINE
extern "C" void set_import_interface_fields_necessity
#else
extern "C" void set_import_interface_fields_necessity_
#endif
(int *import_interface_id, int *necessity, int *size_necessity, const char *annotation)
{
    inout_interface_mgr->get_interface(*import_interface_id)->set_fields_necessity(necessity, *size_necessity, annotation);
}


#ifdef LINK_WITHOUT_UNDERLINE
extern "C" void get_ccpl_import_fields_sender_time
#else
extern "C" void get_ccpl_import_fields_sender_time_
#endif
(int *import_interface_id, int *size_sender_date, int *size_sender_elapsed_days, int *size_sender_second, int *sender_date, int *sender_elapsed_days, int *sender_second, const char *annotation)
{
	EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Start to get the sender time of import fields");
	check_for_ccpl_managers_allocated(API_ID_INTERFACE_GET_SENDER_TIME, annotation);
    Inout_interface *import_interface = inout_interface_mgr->get_interface(*import_interface_id);
    EXECUTION_REPORT(REPORT_ERROR, -1, import_interface != NULL, "ERROR happens when calling the API \"CCPL_get_import_fields_sender_time\": the parameter \"interface_id\" is not a legal ID of a coupling interface. Please verify the model code with the annotation \"%s\".", annotation);
	import_interface->get_sender_time(*size_sender_date, *size_sender_elapsed_days, *size_sender_second, sender_date, sender_elapsed_days, sender_second, annotation);
	EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Finish getting the sender time of import fields");
}


#ifdef LINK_WITHOUT_UNDERLINE
extern "C" void execute_inout_interface_with_id
#else
extern "C" void execute_inout_interface_with_id_
#endif
(int *interface_id, int *bypass_timer, int *field_update_status, int *size_field_update_status, int *num_dst_fields, const char *annotation)
{
    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Start to execute an interface");

    check_for_ccpl_managers_allocated(API_ID_INTERFACE_EXECUTE_WITH_ID, annotation);
    inout_interface_mgr->execute_interface(*interface_id, API_ID_INTERFACE_EXECUTE_WITH_ID, *bypass_timer == 1, field_update_status, *size_field_update_status, num_dst_fields, annotation);

    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Finish executing an interface");
}


#ifdef LINK_WITHOUT_UNDERLINE
extern "C" void execute_inout_interface_with_name
#else
extern "C" void execute_inout_interface_with_name_
#endif
(int *comp_id, const char *interface_name, int *bypass_timer, int *field_update_status, int *size_field_update_status, int *num_dst_fields, const char *annotation)
{
    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Start to execute an interface");

    check_for_ccpl_managers_allocated(API_ID_INTERFACE_EXECUTE_WITH_NAME, annotation);
    check_API_parameter_string_length(-1, API_ID_INTERFACE_EXECUTE_WITH_NAME, CCPL_NAME_STR_LEN, interface_name, "interface_name", annotation);
    inout_interface_mgr->execute_interface(*comp_id, API_ID_INTERFACE_EXECUTE_WITH_NAME, interface_name, *bypass_timer == 1, field_update_status, *size_field_update_status, num_dst_fields, annotation);

    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Finish executing an interface");
}


#ifdef LINK_WITHOUT_UNDERLINE
extern "C" void get_local_comp_full_name
#else
extern "C" void get_local_comp_full_name_
#endif
(int *comp_id, char *comp_full_name, int *comp_full_name_size, const char *annotation)
{
    check_for_component_registered(*comp_id, API_ID_COMP_MGT_GET_LOCAL_COMP_FULL_NAME, annotation, false);
    const char *full_name = comp_comm_group_mgt_mgr->get_global_node_of_local_comp(*comp_id, true, "in get_local_comp_full_name_")->get_full_name();
    copy_out_string_to_Fortran_API(*comp_id, *comp_full_name_size, comp_full_name, full_name, API_ID_COMP_MGT_GET_LOCAL_COMP_FULL_NAME, "comp_full_name", annotation);
}


#ifdef LINK_WITHOUT_UNDERLINE
extern "C" void ccpl_report
#else
extern "C" void ccpl_report_
#endif
(int *report_type, int *comp_id, int *condition, const char *report_content, const char *annotation)
{
    int API_id;
    bool local_condition = *condition == 1? true : false;

    
    if (*report_type == REPORT_ERROR)
        API_id = API_ID_REPORT_ERROR;
    else if (*report_type == REPORT_EXTERNAL_LOG)
        API_id = API_ID_REPORT_LOG;
    else API_id = API_ID_REPORT_PROGRESS;

    check_for_ccpl_managers_allocated(API_id, annotation);
    check_API_parameter_string_length(*comp_id, API_id, 512, report_content, "report_string", annotation);
    EXECUTION_REPORT(*report_type, *comp_id, local_condition, report_content);
}

