/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#include "CCPL_api_mgt.h"
#include "global_data.h"
#include <unistd.h>
#include "quick_sort.h"



long calculate_checksum_of_array(const char *data_array, int array_size, int data_type_size, const char *original_data_type, const char *new_data_type)
{
    const char *local_buffer = data_array;
    long total_checksum = 0, temp_checksum = 0;


    for (int i = 0; i < array_size*data_type_size/sizeof(long); i ++)
        total_checksum += ((const long*) local_buffer)[i] * (i+1);
    for (int i = (array_size*data_type_size/sizeof(long))*sizeof(long); i < array_size*data_type_size; i ++) 
        temp_checksum = (temp_checksum << 8) | local_buffer[i];
    total_checksum += temp_checksum * ((array_size*data_type_size/sizeof(long))+1);

    if (local_buffer != data_array)
        delete [] local_buffer;

    return total_checksum;
}


void get_API_hint(int comp_id, int API_id, char *API_label)
{
    switch(API_id) {
        case API_ID_FINALIZE:
            sprintf(API_label, "CCPL_finalize");
            break;
        case API_ID_COMP_MGT_REG_COMP:
            sprintf(API_label, "CCPL_register_component");
            break;
        case API_ID_COMP_MGT_GET_COMP_LOG_FILE_NAME:
            sprintf(API_label, "CCPL_get_comp_log_file_name");
            break;            
        case API_ID_COMP_MGT_GET_COMP_LOG_FILE_DEVICE:
            sprintf(API_label, "CCPL_get_comp_log_file_device");
            break;
        case API_ID_COMP_MGT_END_COMP_REG:
            sprintf(API_label, "CCPL_end_coupling_configuration");
            break;
        case API_ID_COMP_MGT_GET_COMP_ID:
            sprintf(API_label, "CCPL_get_component_id");
            break;
        case API_ID_COMP_MGT_IS_COMP_TYPE_COUPLED:
            sprintf(API_label, "CCPL_is_comp_type_coupled");
            break;            
        case API_ID_COMP_MGT_IS_CURRENT_PROC_IN_COMP:
            sprintf(API_label, "CCPL_is_current_process_in_component");
            break;
        case API_ID_COMP_MGT_GET_CURRENT_PROC_ID_IN_COMP:
            sprintf(API_label, "CCPL_get_current_process_id_in_component");
            break;
        case API_ID_COMP_MGT_GET_NUM_PROC_IN_COMP:
            sprintf(API_label, "CCPL_get_num_process_in_component");
            break;    
        case API_ID_COMP_MGT_GET_COMP_PROC_GLOBAL_ID:
            sprintf(API_label, "CCPL_get_component_process_global_id");
            break;    
        case API_ID_GRID_MGT_REG_H2D_GRID_VIA_LOCAL_DATA:
            sprintf(API_label, "CCPL_register_H2D_grid_via_local_data");
            break;
        case API_ID_GRID_MGT_REG_H2D_GRID_VIA_GLOBAL_DATA:
            sprintf(API_label, "CCPL_register_H2D_grid_via_global_data");
            break;
        case API_ID_GRID_MGT_REG_H2D_GRID_VIA_FILE:
            sprintf(API_label, "CCPL_register_H2D_grid_via_data_file");
            break;
        case API_ID_GRID_MGT_REG_1D_GRID_ONLINE:
            sprintf(API_label, "CCPL_register_1D_grid");
            break;
        case API_ID_GRID_MGT_REG_GRID_VIA_COR:
            sprintf(API_label, "CCPL_register_CoR_defined_grid");
            break;
        case API_ID_GRID_MGT_GET_GRID_SIZE:
            sprintf(API_label, "CCPL_get_grid_size");
            break;            
        case API_ID_GRID_MGT_REG_GRID_VIA_LOCAL:
            sprintf(API_label, "CCPL_get_local_grid");
            break;
        case API_ID_GRID_MGT_REG_H2D_GRID_VIA_COMP:
            sprintf(API_label, "CCPL_register_H2D_grid_from_another_component");
            break;
        case API_ID_GRID_MGT_CMP_GRID_VIA_REMOTE:
            sprintf(API_label, "CCPL_compare_to_remote_grid");
            break;
        case API_ID_GRID_MGT_GET_GRID_ID:
            sprintf(API_label, "CCPL_get_grid_id");
            break;
        case API_ID_GRID_MGT_SET_GRID_DATA:
            sprintf(API_label, "CCPL_set_grid_data");
            break;
        case API_ID_GRID_MGT_SET_3D_GRID_3D_VERT_FLD:
            sprintf(API_label, "CCPL_set_3D_grid_3D_vertical_coord_field");
            break;
        case API_ID_GRID_MGT_SET_3D_GRID_DYN_BOT_FLD:
            sprintf(API_label, "CCPL_set_3D_grid_dynamic_surface_field");
            break;
        case API_ID_GRID_MGT_SET_3D_GRID_STATIC_BOT_FLD:
            sprintf(API_label, "CCPL_set_3D_grid_static_surface_field");
            break;
        case API_ID_GRID_MGT_SET_3D_GRID_EXTERNAL_BOT_FLD:
            sprintf(API_label, "CCPL_set_3D_grid_external_surface_field");
            break;            
        case API_ID_GRID_MGT_GET_H2D_GRID_DATA:
            sprintf(API_label, "CCPL_get_H2D_grid_data");
            break;
        case API_ID_GRID_MGT_GET_H2D_GRID_AREA_FROM_WGTS:
            sprintf(API_label, "CCPL_get_H2D_grid_area_in_remapping_wgts");
            break;
        case API_ID_GRID_MGT_REG_MID_POINT_GRID:
            sprintf(API_label, "CCPL_register_mid_point_grid");
            break;
        case API_ID_GRID_MGT_REG_V1D_GRID_NO_DATA:
            sprintf(API_label, "CCPL_register_V1D_grid_without_data");
            break; 
        case API_ID_GRID_MGT_REG_V1D_Z_GRID_VIA_MODEL:
            sprintf(API_label, "CCPL_register_V1D_Z_grid_via_model_data");
            break;    
        case API_ID_GRID_MGT_REG_V1D_SIGMA_GRID_VIA_MODEL:
            sprintf(API_label, "CCPL_register_V1D_SIGMA_grid_via_model_data");
            break;    
        case API_ID_GRID_MGT_REG_V1D_HYBRID_GRID_VIA_MODEL:
            sprintf(API_label, "CCPL_register_V1D_HYBRID_grid_via_model_data");
            break;    
        case API_ID_GRID_MGT_REG_MD_GRID_VIA_MULTI_GRIDS:
            sprintf(API_label, "CCPL_register_MD_grid_via_multi_grids");
            break;            
        case API_ID_DECOMP_MGT_REG_DECOMP:
            sprintf(API_label, "CCPL_register_normal_parallel_decomp");
            break;
        case API_ID_FIELD_MGT_REG_FIELD_INST:
            sprintf(API_label, "CCPL_register_field_instance");
            break;
        case API_ID_TIME_MGT_SET_NORMAL_TIME_STEP:
            sprintf(API_label, "CCPL_set_normal_time_step");
            break;
        case API_ID_TIME_MGT_RESET_TIME_TO_START:
            sprintf(API_label, "CCPL_reset_current_time_to_start_time");
            break;            
        case API_ID_TIME_MGT_ADVANCE_TIME:
            sprintf(API_label, "CCPL_advance_time");
            break;
        case API_ID_TIME_MGT_GET_CURRENT_NUM_DAYS_IN_YEAR:
            sprintf(API_label, "CCPL_get_current_num_days_in_year");
            break;            
        case API_ID_TIME_MGT_GET_CURRENT_YEAR:
            sprintf(API_label, "CCPL_get_current_year");
            break;                        
        case API_ID_TIME_MGT_GET_CURRENT_DATE:
            sprintf(API_label, "CCPL_get_current_date");
            break;                        
        case API_ID_TIME_MGT_GET_CURRENT_SECOND:
            sprintf(API_label, "CCPL_get_current_second");
            break;                        
        case API_ID_TIME_MGT_GET_START_TIME:
            sprintf(API_label, "CCPL_get_start_time");
            break;                        
        case API_ID_TIME_MGT_GET_STOP_TIME:
            sprintf(API_label, "CCPL_get_stop_time");
            break;            
        case API_ID_TIME_MGT_GET_PREVIOUS_TIME:
            sprintf(API_label, "CCPL_get_previous_time");
            break;                        
        case API_ID_TIME_MGT_GET_CURRENT_TIME:
            sprintf(API_label, "CCPL_get_current_time");
            break;            
        case API_ID_TIME_MGT_GET_ELAPSED_DAYS_FROM_REF:
            sprintf(API_label, "CCPL_get_num_elapsed_days_from_reference");
            break;            
        case API_ID_TIME_MGT_GET_ELAPSED_DAYS_FROM_START:
            sprintf(API_label, "CCPL_get_num_elapsed_days_from_start");
            break;            
        case API_ID_TIME_MGT_IS_END_CURRENT_DAY:
            sprintf(API_label, "CCPL_is_end_current_day");
            break;            
        case API_ID_TIME_MGT_IS_END_CURRENT_MONTH:
            sprintf(API_label, "CCPL_is_end_current_month");
            break;            
        case API_ID_TIME_MGT_GET_CURRENT_CAL_TIME:
            sprintf(API_label, "CCPL_get_current_calendar_time");
            break;                                            
        case API_ID_TIME_MGT_IS_FIRST_STEP:
            sprintf(API_label, "CCPL_is_first_step");
            break;
        case API_ID_TIME_MGT_IS_FIRST_RESTART_STEP:
            sprintf(API_label, "CCPL_is_first_restart_step");
            break;            
        case API_ID_TIME_MGT_GET_NUM_CURRENT_STEP:
            sprintf(API_label, "CCPL_get_number_of_current_step");
            break;
        case API_ID_TIME_MGT_GET_NUM_TOTAL_STEPS:
            sprintf(API_label, "CCPL_get_number_of_total_steps");
            break;
        case API_ID_TIME_MGT_GET_NORMAL_TIME_STEP:
            sprintf(API_label, "CCPL_get_normal_time_step");
            break;
        case API_ID_TIME_MGT_CHECK_CURRENT_TIME:
            sprintf(API_label, "CCPL_check_current_time");
            break;
        case API_ID_TIME_MGT_IS_TIMER_ON:
            sprintf(API_label, "CCPL_is_timer_on");
            break;
        case API_ID_TIME_MGT_IS_MODEL_RUN_ENDED:
            sprintf(API_label, "CCPL_is_model_run_ended");
            break;
        case API_ID_TIME_MGT_IS_MODEL_LAST_STEP:
            sprintf(API_label, "CCPL_is_last_step_of_model_run");
            break;            
        case API_ID_INTERFACE_REG_IMPORT:
            sprintf(API_label, "CCPL_register_import_interface");
            break;
        case API_ID_INTERFACE_REG_EXPORT:
            sprintf(API_label, "CCPL_register_export_interface");
            break;
        case API_ID_INTERFACE_REG_NORMAL_REMAP: 
            sprintf(API_label, "CCPL_register_normal_remap_interface");
            break;
        case API_ID_INTERFACE_REG_FRAC_REMAP:
            sprintf(API_label, "CCPL_register_frac_based_remap_interface");
            break;
        case API_ID_INTERFACE_EXECUTE_WITH_ID:
            sprintf(API_label, "CCPL_execute_interface_using_id");
            break;
        case API_ID_INTERFACE_EXECUTE_WITH_NAME:
            sprintf(API_label, "CCPL_execute_interface_using_name");
            break;
        case API_ID_INTERFACE_CHECK_IMPORT_FIELD_CONNECTED:
            sprintf(API_label, "CCPL_check_is_import_field_connected");
            break;
		case API_ID_INTERFACE_GET_SENDER_TIME:
			sprintf(API_label, "CCPL_get_import_fields_sender_time");
			break;
        case API_ID_COMP_MGT_GET_LOCAL_COMP_FULL_NAME:
            sprintf(API_label, "CCPL_get_local_comp_full_name");
            break;
        case API_ID_TIME_MGT_DEFINE_SINGLE_TIMER:
            sprintf(API_label, "CCPL_define_single_timer");
            break;
        case API_ID_TIME_MGT_DEFINE_COMPLEX_TIMER:
            sprintf(API_label, "CCPL_define_complex_timer");
            break;
        case API_ID_FIELD_MGT_REG_IO_FIELD_from_INST:
            sprintf(API_label, "CCPL_register_IO_field_from_field_instance");
            break;
        case API_ID_FIELD_MGT_REG_IO_FIELDs_from_INSTs:
            sprintf(API_label, "CCPL_register_IO_fields_from_field_instances");
            break;
        case API_ID_FIELD_MGT_REG_IO_FIELD_from_BUFFER:
            sprintf(API_label, "CCPL_register_IO_field_from_data_buffer");
            break;
        case API_ID_REPORT_LOG:
            sprintf(API_label, "CCPL_report_log");
            break;            
        case API_ID_REPORT_PROGRESS:
            sprintf(API_label, "CCPL_report_progress");
            break;    
        case API_ID_REPORT_ERROR:
            sprintf(API_label, "CCPL_report_error");
            break;    
        case API_ID_RESTART_MGT_START_READ_IO:
            sprintf(API_label, "CCPL_start_restart_read_IO");
            break;
        case API_ID_RESTART_MGT_READ_ALL:
            sprintf(API_label, "CCPL_restart_read_fields_all");
            break;            
        case API_ID_RESTART_MGT_READ_INTERFACE:
            sprintf(API_label, "CCPL_restart_read_fields_interface");
            break;    
        case API_ID_RESTART_MGT_GET_SETTING:
            sprintf(API_label, "CCPL_get_restart_setting");
            break;                
        case API_ID_RESTART_MGT_IS_TIMER_ON:
            sprintf(API_label, "CCPL_is_restart_timer_on");
            break;
        case API_ID_RESTART_MGT_WRITE_IO:
            sprintf(API_label, "CCPL_do_restart_write_IO");
            break;
        case API_ID_COUPLING_GEN_FAMILY:
            sprintf(API_label, "CCPL_do_family_coupling_generation");
            break;
        case API_ID_COUPLING_GEN_INDIVIDUAL:
            sprintf(API_label, "CCPL_do_individual_coupling_generation");
            break;
        case API_ID_COUPLING_GEN_EXTERNAL:
            sprintf(API_label, "CCPL_do_external_coupling_generation");
            break;
        case API_ID_COUPLING_GEN_GET_COMPS:
            sprintf(API_label, "CCPL_get_configurable_comps_full_names");
            break;
        default:
            EXECUTION_REPORT(REPORT_ERROR, comp_id, false, "software error1 in get_API_hint %x", API_id);
            break;
    }
}


void synchronize_comp_processes_for_API(int comp_id, int API_id, MPI_Comm comm, const char *hint, const char *annotation)
{
    char API_label_local[NAME_STR_SIZE], API_label_another[NAME_STR_SIZE];
    int local_process_id, num_processes;
    int *API_ids;
    char *annotations, *comp_names, local_annotation[NAME_STR_SIZE];


    if (!report_error_enabled)
        return;
    
    get_API_hint(-1, API_id, API_label_local);

    if (comp_id != -1)
        EXECUTION_REPORT(REPORT_ERROR, -1, comp_comm_group_mgt_mgr->is_legal_local_comp_id(comp_id, false), "Error happens when calling the interface \"%s\" for %s: the given component model ID (0x%x). Please check the model code with the annotation \"%s\"", API_label_local, hint, comp_id, annotation);

    if (comm == MPI_COMM_NULL)
        comm = comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id, "in synchronize_comp_processes_for_API");

    if (hint != NULL) {
        EXECUTION_REPORT_LOG(REPORT_LOG, comp_id, true, "Before the MPI_barrier for synchronizing all processes of a communicator for %s at C-Coupler API \"%s\" with model code annotation \"%s\"", hint, API_label_local, annotation);    
    }
    else EXECUTION_REPORT_LOG(REPORT_LOG, comp_id, true, "Before the MPI_barrier for synchronizing all processes of a communicator at C-Coupler API \"%s\" with model code annotation \"%s\"", API_label_local, annotation);
    MPI_Barrier(comm);
    if (hint != NULL) {
        EXECUTION_REPORT_LOG(REPORT_LOG, comp_id, true, "After the MPI_barrier for synchronizing all processes of a communicator for %s at C-Coupler API \"%s\" with model code annotation \"%s\"", hint, API_label_local, annotation);    
    }    
    else EXECUTION_REPORT_LOG(REPORT_LOG, comp_id, true, "After the MPI_barrier for synchronizing all processes of a communicator at C-Coupler API \"%s\" with model code annotation \"%s\"", API_label_local, annotation);
    EXECUTION_REPORT(REPORT_ERROR, -1, MPI_Comm_rank(comm, &local_process_id) == MPI_SUCCESS);
    EXECUTION_REPORT(REPORT_ERROR,-1, MPI_Comm_size(comm, &num_processes) == MPI_SUCCESS);
    API_ids = new int [num_processes];
    annotations = new char [num_processes*NAME_STR_SIZE];
    EXECUTION_REPORT(REPORT_ERROR,-1, MPI_Gather(&API_id, 1, MPI_INT, API_ids, 1, MPI_INT, 0, comm) == MPI_SUCCESS);
    EXECUTION_REPORT_LOG(REPORT_LOG, comp_id, true, "annotation is \"%s\"", annotation);
    EXECUTION_REPORT(REPORT_ERROR, -1, strlen(annotation) < NAME_STR_SIZE, "Error happens when calling the API \"%s\": the annotation is too long (%d characters, larger than %d). Please verify", annotation, strlen(annotation), NAME_STR_SIZE);
    strcpy(local_annotation, annotation);
    EXECUTION_REPORT(REPORT_ERROR,-1, MPI_Gather((void*)local_annotation, NAME_STR_SIZE, MPI_CHAR, annotations, NAME_STR_SIZE, MPI_CHAR, 0, comm) == MPI_SUCCESS);
    if (local_process_id == 0) {
        for (int i = 1; i < num_processes; i ++) {
            get_API_hint(comp_id, API_ids[i], API_label_another);
            EXECUTION_REPORT(REPORT_ERROR, comp_id, API_id == API_ids[i], "Different kinds of C-Coupler API calls (\"%s\" and \"%s\") are mapped to the same synchronization. Please check the model code related to the annotations \"%s\" and \"%s\".",
                             API_label_local, API_label_another, annotation, annotations+NAME_STR_SIZE*i);            
        }    
    }
    if (comp_id != -1) {
        comp_names = new char [num_processes*NAME_STR_SIZE];
        EXECUTION_REPORT(REPORT_ERROR, -1, MPI_Gather((void*)comp_comm_group_mgt_mgr->get_global_node_of_local_comp(comp_id,false,"C-Coupler gets component node in synchronize_comp_processes_for_API")->get_comp_full_name(), NAME_STR_SIZE, MPI_CHAR, comp_names, NAME_STR_SIZE, MPI_CHAR, 0, comm) == MPI_SUCCESS);
        if (local_process_id == 0) {
            for (int i = 1; i < num_processes; i ++)
                EXECUTION_REPORT(REPORT_ERROR, comp_id, words_are_the_same(comp_names, comp_names+NAME_STR_SIZE*i), "It is wrong that two different component models (\"%s\" and \"%s\") take part in the same API (\"%s\"). Please check the model code related to the annotations \"%s\" and \"%s\".",
                             comp_names, comp_names+NAME_STR_SIZE*i, API_label_local, annotation, annotations+NAME_STR_SIZE*i);        
        }
        delete [] comp_names;
    }
    
    delete [] API_ids;
    delete [] annotations;
}


template <class T> void check_API_parameter_scalar(int comp_id, int API_id, MPI_Comm comm, const char *hint, T value, const char *parameter_name, const char *annotation)
{
    int i, local_process_id, num_processes;
    T *values;
    char API_label[NAME_STR_SIZE];
    

    if (!report_error_enabled)
        return;

    EXECUTION_REPORT(REPORT_ERROR, -1, MPI_Comm_rank(comm, &local_process_id) == MPI_SUCCESS);
    EXECUTION_REPORT(REPORT_ERROR, -1, MPI_Comm_size(comm, &num_processes) == MPI_SUCCESS);    

    values = new T [num_processes];
    if (sizeof(T) == 1)
        EXECUTION_REPORT(REPORT_ERROR, -1, MPI_Gather(&value, 1, MPI_CHAR, values, 1, MPI_CHAR, 0, comm) == MPI_SUCCESS);
    else if (sizeof(T) == 2)
        EXECUTION_REPORT(REPORT_ERROR, -1, MPI_Gather(&value, 1, MPI_SHORT, values, 1, MPI_SHORT, 0, comm) == MPI_SUCCESS);
    else if (sizeof(T) == 4)
        EXECUTION_REPORT(REPORT_ERROR, -1, MPI_Gather(&value, 1, MPI_INT, values, 1, MPI_INT, 0, comm) == MPI_SUCCESS);
    else if (sizeof(T) == 8)
        EXECUTION_REPORT(REPORT_ERROR, -1, MPI_Gather(&value, 1, MPI_DOUBLE, values, 1, MPI_DOUBLE, 0, comm) == MPI_SUCCESS);
    else EXECUTION_REPORT(REPORT_ERROR, comp_id, true, "software error in check_API_parameter_scalar");
    if (local_process_id == 0) {
        get_API_hint(comp_id, API_id, API_label);
        for (i = 1; i < num_processes; i ++) {
            if (hint != NULL)
                EXECUTION_REPORT(REPORT_ERROR, comp_id, values[0] == values[i], "Error happens when calling the API \"%s\" for %s: parameter %s is not consistent among processes of component \"%s\". Please check the model code related to the annotation \"%s\"",
                                 API_label, hint, parameter_name, comp_comm_group_mgt_mgr->search_global_node(comp_id)->get_comp_name(), annotation);
            else EXECUTION_REPORT(REPORT_ERROR, comp_id, values[0] == values[i], "Error happens when calling the API \"%s\": parameter %s is not consistent among processes of component \"%s\". Please check the model code related to the annotation \"%s\"",
                                  API_label, parameter_name, comp_comm_group_mgt_mgr->search_global_node(comp_id)->get_comp_name(), annotation);
        }
    }

    delete [] values;
}


char *check_and_aggregate_local_grid_data(int comp_id, int API_id, MPI_Comm comm, const char *hint, int grid_size, int array_size, int data_type_size, char *array_value, 
                                          const char *parameter_name, int num_local_cells, int *local_cells_global_index, int &grid_data_size, const char *annotation)
{
    char API_label[NAME_STR_SIZE];
    int local_process_id, num_processes, *counts_for_cell_index, *displs_for_cell_index, *counts_for_array, *displs_for_array;
    int num_total_cells, total_array_size, *all_local_cells_global_index, num_point;
    char *all_array_values, *grid_data;


    grid_data_size = 0;
    get_API_hint(comp_id, API_id, API_label);
    
    int parameter_specified = array_size >= 0 ? 1 : 0;
    check_API_parameter_int(comp_id, API_id, comm, "specification (or not)", parameter_specified, parameter_name, annotation);
    if (parameter_specified == 0)
        return NULL;

    check_API_parameter_int(comp_id, API_id, comm, "data type", data_type_size, parameter_name, annotation);
    if (words_are_the_same(parameter_name, "vertex_lon") || words_are_the_same(parameter_name, "vertex_lat"))
        EXECUTION_REPORT(REPORT_ERROR, comp_id, array_size == 0 && num_local_cells == 0 || (array_size % num_local_cells) == 0, "Error happens when calling the API \"%s\" for %s: the array size (currently is %d) of the parameter \"%s\" is not an integer multiple of parameter \"num_local_cells\" (currently is %d). Please check the model code related to the annotation \"%s\"", API_label, hint, array_size, parameter_name, num_local_cells, annotation);
    else EXECUTION_REPORT(REPORT_ERROR, comp_id, array_size == num_local_cells, "Error happens when calling the API \"%s\" for %s: the array size (currently is %d) of the parameter \"%s\" must be the same as \"num_local_cells\" (currently is %d). Please check the model code related to the annotation \"%s\"", API_label, hint, array_size, parameter_name, num_local_cells, annotation);

    EXECUTION_REPORT(REPORT_ERROR, -1, MPI_Comm_rank(comm, &local_process_id) == MPI_SUCCESS);
    EXECUTION_REPORT(REPORT_ERROR, -1, MPI_Comm_size(comm, &num_processes) == MPI_SUCCESS);    

    if (local_process_id == 0) {
        counts_for_cell_index = new int [num_processes];
        displs_for_cell_index = new int [num_processes];
        counts_for_array = new int [num_processes];
        displs_for_array = new int [num_processes];
    }
    EXECUTION_REPORT(REPORT_ERROR, -1, MPI_Gather(&num_local_cells, 1, MPI_INT, counts_for_cell_index, 1, MPI_INT, 0, comm) == MPI_SUCCESS);
    EXECUTION_REPORT(REPORT_ERROR, -1, MPI_Gather(&array_size, 1, MPI_INT, counts_for_array, 1, MPI_INT, 0, comm) == MPI_SUCCESS);
    if (local_process_id == 0) {
        num_total_cells = 0;
        total_array_size = 0;
        for (int i = 0; i < num_processes; i ++) {
            displs_for_cell_index[i] = num_total_cells;
            displs_for_array[i] = total_array_size;
            num_total_cells += counts_for_cell_index[i];
            total_array_size += counts_for_array[i];
        }
        if (num_total_cells != 0)
            all_local_cells_global_index = new int [num_total_cells];
        else all_local_cells_global_index = new int [1];
        if (total_array_size != 0)
            all_array_values = new char [total_array_size*data_type_size];
        else all_array_values = new char [8];
    }
    EXECUTION_REPORT(REPORT_ERROR,-1, MPI_Gatherv(local_cells_global_index, num_local_cells, MPI_INT, all_local_cells_global_index, counts_for_cell_index, displs_for_cell_index, MPI_INT, 0, comm) == MPI_SUCCESS);
    if (data_type_size == 4)
        EXECUTION_REPORT(REPORT_ERROR,-1, MPI_Gatherv(array_value, array_size, MPI_INT, all_array_values, counts_for_array, displs_for_array, MPI_INT, 0, comm) == MPI_SUCCESS);
    else EXECUTION_REPORT(REPORT_ERROR,-1, MPI_Gatherv(array_value, array_size, MPI_DOUBLE, all_array_values, counts_for_array, displs_for_array, MPI_DOUBLE, 0, comm) == MPI_SUCCESS);
        
    if (local_process_id == 0) {
        if (words_are_the_same(parameter_name, "vertex_lon") || words_are_the_same(parameter_name, "vertex_lat")) {
            num_point = 0;
            for (int i = 0; i < num_processes; i ++) {
                if (counts_for_cell_index[i] == 0)
                    continue;
                if (num_point == 0)
                    num_point = counts_for_array[i] / counts_for_cell_index[i];
                EXECUTION_REPORT(REPORT_ERROR, comp_id, num_point == counts_for_array[i] / counts_for_cell_index[i], "Error happens when calling the API \"%s\" for %s: the number of vertexes corresponding to parameter \"%s\" does not keep the same among the processes. Please check the model code related to the annotation \"%s\"", API_label, hint, parameter_name, annotation);
            }

        }
        else num_point = 1;
        if (num_point == 0) {
            grid_data_size = 0;
            grid_data = NULL;
        }
        else {
            grid_data_size = grid_size * num_point;
            grid_data = new char [grid_data_size*data_type_size];
        }
        if (grid_data_size > 0) {
            memset(grid_data, 0, grid_data_size*data_type_size);
            int *grid_data_mark = new int [grid_size];
            memset(grid_data_mark, 0, grid_size*sizeof(int));
            for (int i = 0; i < num_processes; i ++)
                for (int j = 0; j < counts_for_cell_index[i]; j ++) {
                    int global_index = all_local_cells_global_index[displs_for_cell_index[i]+j]-1;
                    if (grid_data_mark[global_index] == 0) {
                        grid_data_mark[global_index] = 1;
                        memcpy(grid_data+global_index*num_point*data_type_size, all_array_values+(displs_for_array[i]+j*num_point)*data_type_size, num_point*data_type_size);
                    }
                    else {
                        bool is_the_same = memcmp(grid_data+global_index*num_point*data_type_size, all_array_values+(displs_for_array[i]+j*num_point)*data_type_size, num_point*data_type_size) == 0;
                        EXECUTION_REPORT(REPORT_ERROR, comp_id, is_the_same, "Error happens when calling the API \"%s\" for %s: the grid data (\"%s\") of some common cells is not the same among the processes. Please check the model code related to the annotation \"%s\"", API_label, hint, parameter_name, annotation);
                    }
                }
            delete [] grid_data_mark;
            if (report_error_enabled) {
                do_quick_sort(all_local_cells_global_index, (int*)NULL, 0, num_total_cells-1);
                EXECUTION_REPORT(REPORT_ERROR, comp_id, all_local_cells_global_index[0] == 1, "Error happens when calling the API \"%s\" for %s: no process provide grid data (\"%s\") for the first grid cell. Please check the model code related to the annotation \"%s\"", API_label, hint, parameter_name, annotation);
                EXECUTION_REPORT(REPORT_ERROR, comp_id, all_local_cells_global_index[num_total_cells-1] == grid_size, "Error happens when calling the API \"%s\" for %s: no process provide grid data (\"%s\") for the last (%dth) grid cell. Please check the model code related to the annotation \"%s\"", API_label, hint, parameter_name, grid_size, annotation);
                for (int i = 1; i < num_total_cells; i ++)                    
                    EXECUTION_REPORT(REPORT_ERROR, comp_id, all_local_cells_global_index[i] == all_local_cells_global_index[i-1] || all_local_cells_global_index[i] == all_local_cells_global_index[i-1]+1, "Error happens when calling the API \"%s\" for %s: no process provide grid data (\"%s\") for the %dth grid cell. Please check the model code related to the annotation \"%s\"", API_label, hint, parameter_name, all_local_cells_global_index[i-1]+1, annotation);
            }
        }
    }

    EXECUTION_REPORT(REPORT_ERROR, -1, MPI_Bcast(&grid_data_size, 1, MPI_INT, 0, comm) == MPI_SUCCESS);

    if (grid_data_size == 0)
        return NULL;

    if (local_process_id != 0)
        grid_data = new char [grid_data_size*data_type_size];
    EXECUTION_REPORT(REPORT_ERROR,-1, MPI_Bcast(grid_data, grid_data_size*data_type_size, MPI_CHAR, 0, comm) == MPI_SUCCESS);

    if (local_process_id == 0) {
        delete [] counts_for_cell_index;
        delete [] displs_for_cell_index;
        delete [] counts_for_array;
        delete [] displs_for_array;
        delete [] all_array_values;
        delete [] all_local_cells_global_index;
    }

    return grid_data;
}


void check_API_parameter_data_array(int comp_id, int API_id, MPI_Comm comm, const char *hint, int array_size, int data_type_size, const char *array_value, const char *parameter_name, const char *annotation)
{
    long total_checksum = 0;
    char API_label[NAME_STR_SIZE];


    get_API_hint(comp_id, API_id, API_label);
    EXECUTION_REPORT(REPORT_ERROR, comp_id, array_size != 0, "Error happens when calling the API \"%s\" for %s: the parameter array of \"%s\" may have not been allocated. Please check the model code related to the annotation \"%s\"", API_label, hint, parameter_name, annotation);

    int parameter_specified = array_size >= 0 ? 1 : 0;
    check_API_parameter_int(comp_id, API_id, comm, "specification (or not)", parameter_specified, parameter_name, annotation);
    check_API_parameter_int(comp_id, API_id, comm, "array size", array_size, parameter_name, annotation);
    check_API_parameter_int(comp_id, API_id, comm, "data type", data_type_size, parameter_name, annotation);

    if (array_size <= 0)
        return;

    total_checksum = calculate_checksum_of_array(array_value, array_size, data_type_size, NULL, NULL);
    check_API_parameter_long(comp_id, API_id, comm, "array value", total_checksum, parameter_name, annotation);
}


void check_API_parameter_float(int comp_id, int API_id, MPI_Comm comm, const char *hint, float value, const char *parameter_name, const char *annotation)
{
    check_API_parameter_scalar(comp_id, API_id, comm, hint, value, parameter_name, annotation);
}


void check_API_parameter_double(int comp_id, int API_id, MPI_Comm comm, const char *hint, double value, const char *parameter_name, const char *annotation)
{
    check_API_parameter_scalar(comp_id, API_id, comm, hint, value, parameter_name, annotation);
}


void check_API_parameter_bool(int comp_id, int API_id, MPI_Comm comm, const char *hint, bool value, const char *parameter_name, const char *annotation)
{
    check_API_parameter_scalar(comp_id, API_id, comm, hint, value, parameter_name, annotation);
}


void check_API_parameter_int(int comp_id, int API_id, MPI_Comm comm, const char *hint, int value, const char *parameter_name, const char *annotation)
{
    check_API_parameter_scalar(comp_id, API_id, comm, hint, value, parameter_name, annotation);
}


void check_API_parameter_long(int comp_id, int API_id, MPI_Comm comm, const char *hint, long value, const char *parameter_name, const char *annotation)
{
    check_API_parameter_scalar(comp_id, API_id, comm, hint, value, parameter_name, annotation);
}


void check_API_parameter_timer(int comp_id, int API_id, MPI_Comm comm, const char *hint, int timer_id, const char *parameter_name, const char *annotation)
{
    Coupling_timer *timer;

    
    EXECUTION_REPORT(REPORT_ERROR, comp_id, timer_mgr->check_is_legal_timer_id(timer_id), "Software error in check_API_parameter_timer");
    timer = timer_mgr->get_timer(timer_id);
    check_API_parameter_int(comp_id, API_id, comm, hint, timer->get_frequency_count(), parameter_name, annotation);
    check_API_parameter_int(comp_id, API_id, comm, hint, timer->get_local_lag_count(), parameter_name, annotation);
    check_API_parameter_int(comp_id, API_id, comm, hint, timer->get_remote_lag_count(), parameter_name, annotation);
    check_API_parameter_string(comp_id, API_id, comm, hint, timer->get_frequency_unit(), parameter_name, annotation);
}


void check_API_parameter_field_instance(int comp_id, int API_id, MPI_Comm comm, const char *hint, int field_id, const char *parameter_name, const char *annotation)
{
    Field_mem_info *field_instance;
    int decomp_class, is_registered;


    EXECUTION_REPORT(REPORT_ERROR, comp_id, memory_manager->check_is_legal_field_instance_id(field_id), "Software error in check_API_parameter_field_instance");
    field_instance = memory_manager->get_field_instance(field_id);
    if (field_instance->get_decomp_id() == -1)
        decomp_class = -1;
    else decomp_class = 0;
    check_API_parameter_int(comp_id, API_id, comm, hint, decomp_class, parameter_name, annotation);    
    check_API_parameter_string(comp_id, API_id, comm, hint, field_instance->get_field_name(), parameter_name, annotation);
    if (field_instance->get_decomp_id() != -1) {
        check_API_parameter_string(comp_id, API_id, comm, hint, decomps_info_mgr->get_decomp_info(field_instance->get_decomp_id())->get_decomp_name(), parameter_name, annotation);
        check_API_parameter_string(comp_id, API_id, comm, hint, original_grid_mgr->get_name_of_grid(field_instance->get_grid_id()), parameter_name, annotation);
    }
    if (field_instance->get_is_registered_model_buf())
        is_registered = 1;
    else is_registered = 0;
    check_API_parameter_int(comp_id, API_id, comm, hint, is_registered, parameter_name, annotation);    
    if (field_instance->get_is_registered_model_buf())
        check_API_parameter_int(comp_id, API_id, comm, hint, field_instance->get_buf_mark(), parameter_name, annotation);
}


void check_API_parameter_string(int comp_id, int API_id, MPI_Comm comm, const char *hint, const char *string, const char *parameter_name, const char *annotation)
{
    int local_process_id, num_processes, local_string_size, *all_string_size;
    char API_label[NAME_STR_SIZE], *all_string_para;


    if (!report_error_enabled)
        return;

    local_string_size = strlen(string);
    EXECUTION_REPORT(REPORT_ERROR, comp_id, local_string_size > 0, "Error happens when calling the API \"%s\" for %s: parameter %s is an empty string. Please check the model code related to the annotation \"%s\"", API_label, hint, parameter_name, annotation);
    EXECUTION_REPORT(REPORT_ERROR, -1, MPI_Comm_rank(comm, &local_process_id) == MPI_SUCCESS);
    EXECUTION_REPORT(REPORT_ERROR, -1, MPI_Comm_size(comm, &num_processes) == MPI_SUCCESS);
    all_string_size = new int [num_processes];
    
    EXECUTION_REPORT(REPORT_ERROR, -1, MPI_Gather(&local_string_size, 1, MPI_INT, all_string_size, 1, MPI_INT, 0, comm) == MPI_SUCCESS);
    if (local_process_id == 0) {
        get_API_hint(comp_id, API_id, API_label);
        for (int i = 1; i < num_processes; i ++)
            if (comp_id != -1)
                EXECUTION_REPORT(REPORT_ERROR, comp_id, all_string_size[0] == all_string_size[i], "Error happens when calling the API \"%s\" for %s: parameter %s is not consistent among processes of component \"%s\". Please check the model code related to the annotation \"%s\"",
                                 API_label, hint, parameter_name, comp_comm_group_mgt_mgr->search_global_node(comp_id)->get_comp_name(), annotation);
            else EXECUTION_REPORT(REPORT_ERROR, comp_id, all_string_size[0] == all_string_size[i], "Error happens when calling the API \"%s\" for %s: parameter %s is not consistent among processes. Please check the model code related to the annotation \"%s\"",
                                  API_label, hint, parameter_name, annotation);
    }
    all_string_para = new char [local_string_size*num_processes];
    EXECUTION_REPORT(REPORT_ERROR, -1, MPI_Gather((void*)string, local_string_size, MPI_CHAR, all_string_para, local_string_size, MPI_CHAR, 0, comm) == MPI_SUCCESS);
            
    if (local_process_id == 0) {
        for (int i = 1; i < num_processes; i ++)
            if (comp_id != -1)
                EXECUTION_REPORT(REPORT_ERROR, comp_id, strncmp(all_string_para, all_string_para+local_string_size*i, local_string_size) == 0, 
                                 "Error happens when calling the API \"%s\" for %s: parameter %s is not consistent among processes of component \"%s\". Please check the model code related to the annotation \"%s\"",
                                 API_label, hint, parameter_name, comp_comm_group_mgt_mgr->search_global_node(comp_id)->get_comp_name(), annotation);
            else EXECUTION_REPORT(REPORT_ERROR, comp_id, strncmp(all_string_para, all_string_para+local_string_size*i, local_string_size) == 0, 
                                  "Error happens when calling the API \"%s\" for %s: parameter %s is not consistent among processes. Please check the model code related to the annotation \"%s\"",
                                  API_label, hint, parameter_name, annotation);            
    }
    
    delete [] all_string_size;
    delete [] all_string_para;
}


bool check_and_verify_name_format_of_string(const char *string)
{
    for (int i = 0; i < strlen(string); i ++)
        if (!((string[i] >= 'a' && string[i] <= 'z') || (string[i] >= 'A' && string[i] <= 'Z') || (string[i] >= '0' && string[i] <= '9') || string[i] == '_' || string[i] == '-' || string[i] == '.'))
            return false;

    return true;
}


void check_and_verify_name_format_of_string_for_API(int comp_id, const char *string, int API_id, const char *name_owner, const char *annotation)
{
    char API_label[NAME_STR_SIZE];
    int i;


    get_API_hint(comp_id, API_id, API_label);
    
    EXECUTION_REPORT(REPORT_ERROR, comp_id, check_and_verify_name_format_of_string(string),
                     "Error happens when calling the API \"%s\": the name of %s (currently is \"%s\") is in a wrong format. Each character in the name can only be '-', '_', 'a-z', 'A-Z', '0-9', or '.'. Please verify the model code with the annotation \"%s\".",
                     API_label, name_owner, string, annotation);
}


void check_and_verify_name_format_of_string_for_XML(int comp_id, const char *string, const char *name_owner, const char *XML_file_name, int line_number)
{
    EXECUTION_REPORT(REPORT_ERROR, comp_id, check_and_verify_name_format_of_string(string),
                     "When reading the XML file \"%s\", the format of the name of %s (currently is \"%s\") is wrong. Each character in the name can only be '-', '_', 'a-z', 'A-Z', '0-9', or '.'. Please check the XML file arround the line number %d",
                     XML_file_name, name_owner, string, line_number);
}


TiXmlNode *get_XML_first_child_of_unique_root(int comp_id, const char *XML_file_name, TiXmlDocument *XML_file)
{
    TiXmlNode *root_node = XML_file->FirstChildElement();
    if (root_node == NULL)
        return NULL;

    EXECUTION_REPORT(REPORT_ERROR, comp_id, root_node->NextSibling() == NULL && words_are_the_same(root_node->Value(), "root"), "ERROR happens when reading the XML configuration file \"%s\": it now has multiple root nodes while it is allowed to have at most one root node that is named \"root\". Please verify.", XML_file_name);
    
    return root_node->FirstChild();
}


bool is_XML_setting_on(int comp_id, TiXmlElement *XML_element, const char *XML_file_name, const char *attribute_annotation, const char *XML_file_annotation)
{
    int line_number;
    const char *status = get_XML_attribute(comp_id, -1, XML_element, "status", XML_file_name, line_number, attribute_annotation, XML_file_annotation, true);
    EXECUTION_REPORT(REPORT_ERROR, comp_id, words_are_the_same(status, "on") || words_are_the_same(status, "off"), "In the XML file \"%s\" that is for %s, the value of %s is wrong (must be \"on\" or \"off\"). Please verify the XML file arround the line number %d.",
                     XML_file_name, XML_file_annotation, attribute_annotation, line_number);
    return words_are_the_same(status, "on");
}


const char *get_XML_attribute(int comp_id, int max_string_length, TiXmlElement *XML_element, const char *attribute_keyword, const char *XML_file_name, int &line_number, const char *attribute_annotation, const char *XML_file_annotation, bool check_existence)
{
    const char *attribute_value = XML_element->Attribute(attribute_keyword, &line_number);
    if (!check_existence && attribute_value == NULL)
        return NULL;
    EXECUTION_REPORT(REPORT_ERROR, comp_id, attribute_value != NULL, "In the XML file \"%s\" that is for %s, %s (the keyword is \"%s\") has not been specified. Please verify the XML file arround the line number %d.", 
                     XML_file_name, XML_file_annotation, attribute_annotation, attribute_keyword, XML_element->Row());
    EXECUTION_REPORT(REPORT_ERROR, comp_id, strlen(attribute_value) > 0, "In the XML file \"%s\" that is for %s, %s (the keyword is \"%s\") has been specified but with an empty string. Please verify the XML file arround the line number %d.", 
                     XML_file_name, XML_file_annotation, attribute_annotation, attribute_keyword, line_number);
    if (max_string_length > 0)    
        EXECUTION_REPORT(REPORT_ERROR, comp_id, strlen(attribute_value) <= max_string_length, "Error happens when using the XML configuration file \"%s\": the string size (currently is %d) of the value (\"%s\") the XML attribute \"%s\" is larger than the limit (%d). Please verify XML file arround the line %d.", XML_file_name, strlen(attribute_value), attribute_value, attribute_keyword, max_string_length, line_number);
    return attribute_value;
}


void transfer_array_from_one_comp_to_another(int current_proc_local_id_src_comp, int root_proc_global_id_src_comp, int current_proc_local_id_dst_comp, 
                                             int root_proc_global_id_dst_comp, MPI_Comm comm_dst_comp, char **array, long &array_size)
{
    MPI_Status status;

    
    if (current_proc_local_id_src_comp == 0 && current_proc_local_id_dst_comp != 0) {
        MPI_Send(&array_size, 1, MPI_LONG, root_proc_global_id_dst_comp, 0, MPI_COMM_WORLD);
        if (array_size > 0) {
            EXECUTION_REPORT(REPORT_ERROR, -1, *array != NULL, "software error in transfer_array_from_one_comp_to_another");
            MPI_Send(*array, array_size, MPI_CHAR, root_proc_global_id_dst_comp, 0, MPI_COMM_WORLD);
        }
    }
    if (current_proc_local_id_src_comp != 0 && current_proc_local_id_dst_comp == 0) {
        MPI_Recv(&array_size, 1, MPI_LONG, root_proc_global_id_src_comp, 0, MPI_COMM_WORLD, &status);
        if (array_size > 0) {
            if (*array != NULL)
                delete [] *array;
            *array = new char [array_size];
            MPI_Recv(*array, array_size, MPI_CHAR, root_proc_global_id_src_comp, 0, MPI_COMM_WORLD, &status);
        }
    }

    if (current_proc_local_id_dst_comp != -1)
        bcast_array_in_one_comp(current_proc_local_id_dst_comp, array, array_size, comm_dst_comp);
}


void gather_array_in_one_comp(int num_total_local_proc, int current_proc_local_id, void *local_array, int local_array_size, 
                              int data_type_size, int *all_array_size, void **global_array, long &global_size, MPI_Comm comm)
{
    int *displs = new int [num_total_local_proc];
    int *counts = new int [num_total_local_proc];
    bool all_array_size_empty = all_array_size == NULL;


    if (all_array_size_empty)
        all_array_size = new int [num_total_local_proc];
    MPI_Gather(&local_array_size, 1, MPI_INT, all_array_size, 1, MPI_INT, 0, comm);
    global_size = 0;
    if (current_proc_local_id == 0) {
        displs[0] = 0;
        counts[0] = all_array_size[0] * data_type_size;
        global_size = all_array_size[0];
        for (int i = 1; i < num_total_local_proc; i ++) {
            global_size += all_array_size[i];
            counts[i] = all_array_size[i] * data_type_size;
            displs[i] = displs[i-1] + counts[i-1];
        }
        *global_array = new char [displs[num_total_local_proc-1]+counts[num_total_local_proc-1]];
    }
    MPI_Gatherv(local_array, local_array_size*data_type_size, MPI_CHAR, *global_array, counts, displs, MPI_CHAR, 0, comm);

    if (all_array_size_empty)
        delete [] all_array_size;
    delete [] displs;
    delete [] counts;
}


void bcast_array_in_one_comp(int current_proc_local_id, char **array, long &array_size, MPI_Comm comm)
{
    MPI_Bcast(&array_size, 1, MPI_LONG, 0, comm);
    if (array_size > 0) {
        if (current_proc_local_id != 0 && *array != NULL)
            delete [] *array;
        if (current_proc_local_id != 0)
            *array = new char [array_size];
        MPI_Bcast(*array, array_size, MPI_CHAR, 0, comm);
    }
}


bool does_file_exist(const char *file_name)
{
    FILE *tmp_file = fopen(file_name, "r");
    if (tmp_file != NULL)
        fclose(tmp_file);
    return tmp_file != NULL;
}


TiXmlDocument *open_XML_file_to_read(int comp_id, const char *XML_file_name, MPI_Comm comm, bool wait_file)
{
    int local_process_id = 0, file_existing = 0;
    TiXmlDocument *XML_file;
    bool successful;


    EXECUTION_REPORT_LOG(REPORT_LOG, comp_id, true, "Try to load the XML configuration file \"%s\"", XML_file_name);
    if (comm != MPI_COMM_NULL)
        EXECUTION_REPORT(REPORT_ERROR, -1, MPI_Comm_rank(comm, &local_process_id) == MPI_SUCCESS);
    if (local_process_id == 0) {
        do {
            file_existing = does_file_exist(XML_file_name)? 1 : 0;
            if (file_existing == 1)
                break;
            if (wait_file)
                sleep(1);
        } while (wait_file);
    }

    if (comm != MPI_COMM_NULL)
        EXECUTION_REPORT(REPORT_ERROR, -1, MPI_Bcast(&file_existing, 1, MPI_INT, 0, comm) == MPI_SUCCESS);
    
    if (file_existing == 0)
        return NULL;

    for (int i = 0; i < 10; i ++) {
        XML_file = new TiXmlDocument(XML_file_name);
        if (comm != MPI_COMM_NULL)
            successful = XML_file->LoadFile(comm);
        else successful = XML_file->LoadFile();
        if (successful || !wait_file)
            break;
        delete XML_file;
    } 

    if (!successful) {
        EXECUTION_REPORT(REPORT_ERROR, comp_id, false, "Fail to load the XML configuration file \"%s\": the file exists while the format of the content is not legal", XML_file_name);
        delete XML_file;
        return NULL;
    }
    
    EXECUTION_REPORT_LOG(REPORT_LOG, comp_id, true, "Successfully load the XML configuration file \"%s\"", XML_file_name);

    return XML_file;
}

