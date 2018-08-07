/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#include "global_data.h"
#include "restart_mgt.h"



Restart_buffer_container::Restart_buffer_container(const char *comp_full_name, const char *buf_type, const char *keyword, Restart_mgt *restart_mgr)
{
    strcpy(this->comp_full_name, comp_full_name);
    strcpy(this->buf_type, buf_type);
    strcpy(this->keyword, keyword);
    buffer_max_size = 1000;
    buffer_content = new char [buffer_max_size];
    buffer_content_iter = 0;
    buffer_content_size = 0;
    this->restart_mgr = restart_mgr;
}


Restart_buffer_container::Restart_buffer_container(const char *array_buffer, long &buffer_content_iter, const char *file_name, Restart_mgt *restart_mgr)
{
    long total_size, str_size;


    EXECUTION_REPORT(REPORT_ERROR, -1, read_data_from_array_buffer(&total_size, sizeof(long), array_buffer, buffer_content_iter, file_name == NULL), "Fail to load the restart data file \"%s\": its format is wrong, or the information it includes is not complete. Please try a different restart time with complete restart data files.", file_name);
    
    buffer_content = load_string(NULL, str_size, -1, array_buffer, buffer_content_iter, file_name);
    this->buffer_content_iter = str_size;
    this->buffer_content_size = str_size;
    this->restart_mgr = restart_mgr;
    load_string(keyword, str_size, NAME_STR_SIZE, array_buffer, buffer_content_iter, file_name);
    load_string(buf_type, str_size, NAME_STR_SIZE, array_buffer, buffer_content_iter, file_name);
    load_string(comp_full_name, str_size, NAME_STR_SIZE, array_buffer, buffer_content_iter, file_name);
    EXECUTION_REPORT(REPORT_ERROR, -1, total_size == strlen(comp_full_name) + strlen(buf_type) + strlen(keyword) + sizeof(long)*3 + this->buffer_content_iter, "Restart_buffer_container::Restart_buffer_container: wrong format of restart data file");
}


void Restart_buffer_container::dump_in_string(const char *str, long str_size)
{
    dump_string(str, str_size, &buffer_content, buffer_max_size, buffer_content_iter);
    buffer_content_size = buffer_content_iter;
}


void Restart_buffer_container::dump_in_data(const void *data, long size)
{
    write_data_into_array_buffer(data, size, &buffer_content, buffer_max_size, buffer_content_iter);
    buffer_content_size = buffer_content_iter;
}


void Restart_buffer_container::load_restart_data(void *data, long data_size)
{
    EXECUTION_REPORT(REPORT_ERROR, -1, read_data_from_array_buffer(data, data_size, buffer_content, buffer_content_iter, false), "Fail to load the restart data file \"%s\": its format is wrong, or the information it includes is not complete. Please try a different restart time with complete restart data files.", restart_mgr->get_input_restart_mgt_info_file());
}


char *Restart_buffer_container::load_restart_string(char *str, long &str_size, long max_size)
{
    return load_string(str, str_size, max_size, buffer_content, buffer_content_iter, restart_mgr->get_input_restart_mgt_info_file());
}


void Restart_buffer_container::dump_out(char **array_buffer, long &buffer_max_size, long &buffer_content_size)
{
    dump_string(comp_full_name, -1, array_buffer, buffer_max_size, buffer_content_size);
    dump_string(buf_type, -1, array_buffer, buffer_max_size, buffer_content_size);
    dump_string(keyword, -1, array_buffer, buffer_max_size, buffer_content_size);
    dump_string(buffer_content, buffer_content_iter, array_buffer, buffer_max_size, buffer_content_size);
    long total_size = strlen(comp_full_name) + strlen(buf_type) + strlen(keyword) + sizeof(long)*3 + buffer_content_iter;
    write_data_into_array_buffer(&total_size, sizeof(long), array_buffer, buffer_max_size, buffer_content_size);
}


bool Restart_buffer_container::match(const char *buf_type, const char *keyword)
{
    return words_are_the_same(this->buf_type, buf_type) && words_are_the_same(this->keyword, keyword);
}


const char *Restart_buffer_container::get_input_restart_mgt_info_file()
{
    return restart_mgr->get_input_restart_mgt_info_file();
}


Restart_mgt::Restart_mgt(Comp_comm_group_mgt_node *comp_node)
{ 
    this->comp_node = comp_node; 
    last_restart_write_full_time = -1;
    last_restart_write_elapsed_time = -1;
    input_restart_mgt_info_file = NULL;
    restart_read_annotation = NULL;
    restart_mgt_info_written = true;
    time_mgr = NULL;
    restart_write_data_file = NULL;
    backup_restart_write_data_file = NULL;
    restart_read_data_file_name = NULL;
    are_all_restarted_fields_read = false;
    restart_normal_fields_enabled = false;
}


Restart_mgt::~Restart_mgt()
{
    clean(true);
    clean(false);
    if (input_restart_mgt_info_file != NULL)
        delete [] input_restart_mgt_info_file;
    if (restart_read_annotation != NULL)
        delete [] restart_read_annotation;
    if (restart_read_data_file_name != NULL)
        delete restart_read_data_file_name;
    if (backup_restart_write_data_file != NULL)
        delete backup_restart_write_data_file;
}


int Restart_mgt::get_comp_id() 
{ 
    return comp_node->get_comp_id(); 
}



void Restart_mgt::clean(bool is_write_buffers)
{
    if (is_write_buffers) {
        for (int i = 0; i < restart_write_buffer_containers.size(); i ++)
            delete restart_write_buffer_containers[i];
        restart_write_buffer_containers.clear();
    }
    else {
        for (int i = 0; i < restart_read_buffer_containers.size(); i ++)
            delete restart_read_buffer_containers[i];
        restart_read_buffer_containers.clear();        
    }
}


Restart_buffer_container *Restart_mgt::search_restart_buffer(const char *buf_type, const char *keyword)
{
    for (int i = 0; i < restart_read_buffer_containers.size(); i ++)
        if (restart_read_buffer_containers[i]->match(buf_type, keyword))
            return restart_read_buffer_containers[i];

    return NULL;
}


void Restart_mgt::read_restart_mgt_info(const char *specified_file_name, const char *annotation)
{
    char restart_file_full_name[NAME_STR_SIZE], restart_file_short_name[NAME_STR_SIZE];
            

    if (time_mgr == NULL)
        time_mgr = components_time_mgrs->get_time_mgr(comp_node->get_comp_id());
    
    if ((time_mgr->get_runtype_mark() == RUNTYPE_MARK_INITIAL)) {
        EXECUTION_REPORT(REPORT_PROGRESS, comp_node->get_comp_id(), true, "C-Coupler does not read the restart data file because it is a initial run (the run_type is initial)");
        return;
    }
    
    if (strlen(specified_file_name) == 0) {
        if ((time_mgr->get_runtype_mark() == RUNTYPE_MARK_CONTINUE)) {
            get_file_name_in_rpointer_file(restart_file_short_name);
            int special_pos;
            for (special_pos = strlen(restart_file_short_name)-1; special_pos >= 0; special_pos --)
                if (restart_file_short_name[special_pos] == '-')
                    break;
            sprintf(restart_file_short_name+special_pos-8, "%08d-%05d", time_mgr->get_common_restart_full_time()/100000, time_mgr->get_common_restart_full_time()%100000);
            EXECUTION_REPORT_LOG(REPORT_LOG, comp_node->get_comp_id(), true, "The active restart data file is %s corresponding to the restart full time %ld", restart_file_short_name, time_mgr->get_common_restart_full_time());
        }
        else sprintf(restart_file_short_name, "%s.%s.r.%08d-%05d", time_mgr->get_rest_refcase(), comp_node->get_comp_full_name(), time_mgr->get_rest_refdate(), time_mgr->get_rest_refsecond());
    }
    else strcpy(restart_file_short_name, specified_file_name);
    sprintf(restart_file_full_name, "%s/%s", comp_node->get_working_dir(), restart_file_short_name);
    if (input_restart_mgt_info_file != NULL)
        EXECUTION_REPORT(REPORT_ERROR, comp_node->get_comp_id(), false, "Error happens when reading the restart data file \"%s\" at the model code with the annotation \"%s\": a restart data file \"%s\" has already been read before at the model code with the annotation \"%s\", which indicates that the current component model calls the API \"CCPL_start_restart_read_IO\" more than once. Please note that a component model can call this API only once. Please verify.", restart_file_full_name, annotation, input_restart_mgt_info_file, restart_read_annotation);
    input_restart_mgt_info_file = strdup(restart_file_full_name);
    restart_read_annotation = strdup(annotation);

    if (time_mgr->get_runtype_mark() == RUNTYPE_MARK_CONTINUE || time_mgr->get_runtype_mark() == RUNTYPE_MARK_BRANCH)
        EXECUTION_REPORT(REPORT_ERROR, comp_node->get_comp_id(), !time_mgr->get_time_has_been_advanced(), "Error happens when reading the restart file \"%s\": the model time has been advanced, while a component model cannot read any restart file after its model time has been advanced. Please check the model code with the annotation \"%s\"", restart_file_full_name, annotation);
    read_restart_mgt_info((time_mgr->get_runtype_mark() == RUNTYPE_MARK_CONTINUE) || (time_mgr->get_runtype_mark() == RUNTYPE_MARK_BRANCH), restart_file_full_name, annotation);
}


void Restart_mgt::read_restart_mgt_info(bool check_existing_data, const char *file_name, const char *annotation)
{
    int local_proc_id = comp_comm_group_mgt_mgr->get_current_proc_id_in_comp(comp_node->get_comp_id(), "in Restart_mgt::do_restart");
    int temp_int;
    char *array_buffer = NULL;
    char temp_restart_read_data_file_name[NAME_STR_SIZE*2];
    long buffer_content_iter;


    restart_normal_fields_enabled = true;

    if (time_mgr == NULL)
        time_mgr = components_time_mgrs->get_time_mgr(comp_node->get_comp_id());

    if (local_proc_id == 0) {
        FILE *restart_fp = fopen(file_name, "r");
        EXECUTION_REPORT(REPORT_ERROR, comp_node->get_comp_id(), restart_fp != NULL, "Error happens when trying to read restart data at the model code with the annotation \"%s\": the data file \"%s\" does not exist. Please verify.", annotation, file_name);
        fseek(restart_fp, 0, SEEK_END);
        buffer_content_iter = ftell(restart_fp);
        array_buffer = new char [buffer_content_iter];
        fseek(restart_fp, 0, SEEK_SET);
        fread(array_buffer, buffer_content_iter, 1, restart_fp);
    }
    bcast_array_in_one_comp(local_proc_id, &array_buffer, buffer_content_iter, comp_node->get_comm_group());

    int num_restart_buffer_containers;
    EXECUTION_REPORT(REPORT_ERROR, -1, read_data_from_array_buffer(&num_restart_buffer_containers, sizeof(int), array_buffer, buffer_content_iter, false), "Fail to load the restart data file \"%s\": its format is wrong, or the information it includes is not complete. Please try a different restart time with complete restart data files.", file_name);
    for (int i = 0; i < num_restart_buffer_containers; i ++) {
        EXECUTION_REPORT(REPORT_ERROR, comp_node->get_comp_id(), buffer_content_iter > 0, "Software error in Restart_mgt::read_restart_mgt_info: wrong organization of restart data file");
        restart_read_buffer_containers.push_back(new Restart_buffer_container(array_buffer, buffer_content_iter, file_name, this));
    }
    EXECUTION_REPORT(REPORT_ERROR, -1, read_data_from_array_buffer(&temp_int, sizeof(int), array_buffer, buffer_content_iter, false), "Fail to load the restart data file \"%s\": its format is wrong, or the information it includes is not complete. Please try a different restart time with complete restart data files.", file_name);
    bypass_import_fields_at_read = (temp_int == 1);
    EXECUTION_REPORT(REPORT_ERROR, comp_node->get_comp_id(), buffer_content_iter == 0, "Software error in Restart_mgt::read_restart_mgt_info: wrong organization of restart data file");
    delete [] array_buffer;

    if ((time_mgr->get_runtype_mark() == RUNTYPE_MARK_CONTINUE) || (time_mgr->get_runtype_mark() == RUNTYPE_MARK_BRANCH)) {
        Restart_buffer_container *time_mgr_restart_buffer = search_restart_buffer(RESTART_BUF_TYPE_TIME, "local time manager");
        EXECUTION_REPORT(REPORT_ERROR, comp_node->get_comp_id(), time_mgr_restart_buffer != NULL, "Error happens when loading the restart file \"%s\" at the model code with the annotation \"%s\": this file does not include the data for restarting the time information", file_name, annotation);
        long buffer_size = time_mgr_restart_buffer->get_buffer_content_iter();
        time_mgr->import_restart_data(time_mgr_restart_buffer->get_buffer_content(), buffer_size, file_name, check_existing_data);
    }

    sprintf(temp_restart_read_data_file_name, "%s.nc", file_name);
    EXECUTION_REPORT(REPORT_ERROR, comp_node->get_comp_id(), does_file_exist(temp_restart_read_data_file_name), "Error happens when loading the restart data file \"%s\" at the model code with the annotation \"%s\": the file does not exist", temp_restart_read_data_file_name);
    restart_read_data_file_name = strdup(temp_restart_read_data_file_name);
}


void Restart_mgt::do_restart_write(const char *annotation, bool bypass_timer, bool bypass_imported_fields)
{
    int local_proc_id = comp_node->get_current_proc_local_id();
    const char *comp_full_name = comp_node->get_full_name();
    long current_full_time;
    Restart_buffer_container *time_mgr_restart_buffer;


    bypass_import_fields_at_write = bypass_imported_fields;

    if (time_mgr == NULL)
        time_mgr = components_time_mgrs->get_time_mgr(comp_node->get_comp_id());
    current_full_time = time_mgr->get_current_full_time();

    check_API_parameter_bool(comp_node->get_comp_id(), API_ID_RESTART_MGT_WRITE_IO, comp_node->get_comm_group(), "generating restart data files", bypass_timer, "bypass_timer", annotation);
    check_API_parameter_bool(comp_node->get_comp_id(), API_ID_RESTART_MGT_WRITE_IO, comp_node->get_comm_group(), "generating restart data files", bypass_imported_fields, "bypass_timer", annotation);

    if (bypass_timer || time_mgr->is_restart_timer_on()) {
        EXECUTION_REPORT(REPORT_ERROR, comp_node->get_comp_id(), current_full_time != last_restart_write_full_time, "Error happens when the component model tries to write restart data files: the corresponding API \"CCPL_do_restart_write_IO\" has been called more than once at the same time step. Please verify the model code with the annotation \"%s\"", annotation);
        if (time_mgr->get_runtype_mark() == RUNTYPE_MARK_CONTINUE || time_mgr->get_runtype_mark() == RUNTYPE_MARK_BRANCH) {
            EXECUTION_REPORT(REPORT_ERROR, comp_node->get_comp_id(), time_mgr->get_restart_full_time() != ((long)time_mgr->get_current_num_elapsed_day()*((long)100000))+time_mgr->get_current_second(), "Error happens when the component model calls the API \"CCPL_do_restart_write_IO\" to write restart data at the model time %ld: the current model run is a %s run restarted at the same model time, while the model time to write restart data cannot be the same as the restarted model time. Please verify the model code with the annotation \"%s\".", current_full_time, time_mgr->get_run_type(), annotation);
        }    
        last_restart_write_full_time = current_full_time;
        last_restart_write_elapsed_time = time_mgr->get_current_num_elapsed_day()*((long)100000)+time_mgr->get_current_second();
        int date = last_restart_write_full_time/(long)100000;
        int second = last_restart_write_full_time%(long)100000;
        if (local_proc_id == 0) {
            EXECUTION_REPORT(REPORT_ERROR, comp_node->get_comp_id(), restart_write_data_file == NULL, "Error happens when the component model tries to write restart data files: restart writing is too frequent so that a new restart writing starts before the previous restart writing does not finish. Please verify the model code with the annotation \"%s\"", annotation);
            time_mgr_restart_buffer = apply_restart_buffer(comp_full_name, RESTART_BUF_TYPE_TIME, "local time manager");
            time_mgr->write_time_mgt_into_array(time_mgr_restart_buffer->get_buffer_content_ptr(), *(time_mgr_restart_buffer->get_buffer_max_size_ptr()), *(time_mgr_restart_buffer->get_buffer_content_iter_ptr()));
            char restart_data_file_name[NAME_STR_SIZE];
            sprintf(restart_data_file_name, "%s/%s.%s.r.%08d-%05d.nc", comp_node->get_working_dir(), time_mgr->get_case_name(), comp_node->get_comp_full_name(), date, second);
            if (backup_restart_write_data_file != NULL) {
                delete backup_restart_write_data_file;
                backup_restart_write_data_file = NULL;
            }
            restart_write_data_file = new IO_netcdf(restart_data_file_name, restart_data_file_name, "w", false);
            sprintf(restart_data_file_name, "%s/%s.%s.r.%08d-%05d", comp_node->get_working_dir(), time_mgr->get_case_name(), comp_node->get_comp_full_name(), date, second);
            FILE *restart_mgt_info_file = fopen(restart_data_file_name, "w+");
            fclose(restart_mgt_info_file);
        }
        inout_interface_mgr->write_into_restart_buffers(comp_node->get_comp_id());
        restart_mgt_info_written = false;
        for (int i = 0; i < restarted_field_instances.size(); i ++) {
            if (!bypass_imported_fields)
                write_restart_field_data(restarted_field_instances[i].first, NULL, NULL, false);
            else if (!restarted_field_instances[i].second)
                write_restart_field_data(restarted_field_instances[i].first, NULL, NULL, false);
        }
    }
}


void Restart_mgt::get_field_IO_name(char *field_IO_name, Field_mem_info *field_instance, const char *interface_name, const char*label, bool use_time_info)
{
    Field_mem_info *global_field = fields_gather_scatter_mgr->gather_field(field_instance);


    if (interface_name != NULL) {
        if (use_time_info)
            sprintf(field_IO_name, "%s.%s.%s.%13ld", field_instance->get_field_name(), interface_name, label, time_mgr->get_current_full_time());
        else sprintf(field_IO_name, "%s.%s.%s", field_instance->get_field_name(), interface_name, label);    
    }
    else {
        char grid_name[NAME_STR_SIZE], decomp_name[NAME_STR_SIZE];
        sprintf(grid_name, "NULL");
        sprintf(decomp_name, "NULL");
        if (field_instance->get_grid_name() != NULL)
            strcpy(grid_name, field_instance->get_grid_name());
        if (field_instance->get_decomp_name() != NULL)
            strcpy(decomp_name, field_instance->get_decomp_name());
        sprintf(field_IO_name, "%s.%s.%s.%d", field_instance->get_field_name(), grid_name, decomp_name, field_instance->get_buf_mark());
    }
}


void Restart_mgt::write_restart_field_data(Field_mem_info *field_instance, const char *interface_name, const char*label, bool use_time_info)
{
    Field_mem_info *global_field = fields_gather_scatter_mgr->gather_field(field_instance);
    char field_IO_name[NAME_STR_SIZE*2], hint[NAME_STR_SIZE*2];


    get_field_IO_name(field_IO_name, field_instance, interface_name, label, use_time_info);
    if (comp_node->get_current_proc_local_id() == 0) {
        strcpy(global_field->get_field_data()->get_grid_data_field()->field_name_in_IO_file, field_IO_name);
        EXECUTION_REPORT(REPORT_ERROR, -1, restart_write_data_file != NULL && backup_restart_write_data_file == NULL || restart_write_data_file == NULL && backup_restart_write_data_file != NULL, "Software error in Restart_mgt::write_restart_field_data");
        IO_netcdf *active_restart_write_data_file = restart_write_data_file != NULL? restart_write_data_file : backup_restart_write_data_file;
        EXECUTION_REPORT_LOG(REPORT_LOG, comp_node->get_comp_id(), true, "Write variable \"%s\" into restart data file \"%s\"", global_field->get_field_data()->get_grid_data_field()->field_name_in_IO_file, active_restart_write_data_file->get_file_name());
        active_restart_write_data_file->write_grided_data(global_field->get_field_data(), true, -1, -1, true);
        sprintf(hint, "restart writing field \"%s\" to the file \"%s\"", field_IO_name, active_restart_write_data_file->get_file_name());
    }
    else sprintf(hint, "restart writing field \"%s\" to the file", field_IO_name);
    field_instance->check_field_sum(hint);
}


void Restart_mgt::read_all_restarted_fields(const char *annotation)
{
    if (time_mgr == NULL)
        time_mgr = components_time_mgrs->get_time_mgr(comp_node->get_comp_id());        
    if (time_mgr->get_runtype_mark() == RUNTYPE_MARK_INITIAL)
        return;

    are_all_restarted_fields_read = true;
    EXECUTION_REPORT(REPORT_ERROR, comp_node->get_comp_id(), input_restart_mgt_info_file != NULL,  "Error happens when calling the API \"CCPL_restart_read_fields_all\" to read restart fields: the API \"CCPL_start_restart_read_IO\" has not been called before. Please verify the model code corresponding to the annotation \"%s\"", annotation);
    EXECUTION_REPORT(REPORT_ERROR, comp_node->get_comp_id(), !time_mgr->get_time_has_been_advanced(), "Error happens when calling the API \"CCPL_restart_read_fields_all\" to read restart fields: the model time has already been advanced before. Please verify the model code corresponding to the annotation \"%s\"", annotation);
    EXECUTION_REPORT(REPORT_ERROR, comp_node->get_comp_id(), restart_normal_fields_enabled, "Error happens when calling the API \"CCPL_restart_read_fields_all\" to read restart fields: some import interfaces have been executed without bypassing the timer, which is not allowed. Please verify the model code corresponding to the annotation \"%s\"", annotation);

    for (int i = 0; i < restarted_field_instances.size(); i ++)
        read_restart_field_data(restarted_field_instances[i].first, NULL, NULL, false, NULL, bypass_import_fields_at_read&&restarted_field_instances[i].second, annotation);
}


void Restart_mgt::read_restart_field_data(Field_mem_info *field_instance, const char *interface_name, const char *label, bool use_time_info, const char *API_label, bool optional, const char *annotation)
{
    char field_IO_name[NAME_STR_SIZE*2], hint[NAME_STR_SIZE*2];

        
    get_field_IO_name(field_IO_name, field_instance, interface_name, label, use_time_info);

    if (interface_name == NULL && !field_instance->is_checksum_changed()) {
        EXECUTION_REPORT_LOG(REPORT_LOG, comp_node->get_comp_id(), true, "Does not read restart field \"%s\" from the file \"%s\" again at the model code with the annotation \"%s\".", field_IO_name, restart_read_data_file_name, annotation);
        return;
    }

    if (interface_name != NULL)
        restart_normal_fields_enabled = false;

    EXECUTION_REPORT(REPORT_ERROR, comp_node->get_comp_id(), does_file_exist(restart_read_data_file_name), "Error happens when loading the restart data file \"%s\" at the model code with the annotation \"%s\": the file does not exist", restart_read_data_file_name, annotation);
    IO_netcdf *restart_read_data_file = new IO_netcdf(restart_read_data_file_name, restart_read_data_file_name, "r", false);
    bool has_data_in_file = fields_gather_scatter_mgr->read_scatter_field(restart_read_data_file, field_instance, field_IO_name, -1, false);
    delete restart_read_data_file;
    if (!optional && (time_mgr->get_runtype_mark() == RUNTYPE_MARK_CONTINUE || time_mgr->get_runtype_mark() == RUNTYPE_MARK_BRANCH))
        if (interface_name != NULL)
            EXECUTION_REPORT(REPORT_ERROR, comp_node->get_comp_id(), has_data_in_file, "Error happens when loading the restart data file \"%s\" at the model code with the annotation \"%s\": the data file does not contain the variable \"%s\" for the field \"%s\" of the coupling interface \"%s\"", restart_read_data_file_name, annotation, field_IO_name, field_instance->get_field_name(), interface_name);
        else EXECUTION_REPORT(REPORT_ERROR, comp_node->get_comp_id(), has_data_in_file, "Error happens when loading the restart data file \"%s\" at the model code with the annotation \"%s\": the data file does not contain the variable \"%s\" for the field \"%s\"", restart_read_data_file_name, annotation, field_IO_name, field_instance->get_field_name());
    sprintf(hint, "restart reading field \"%s\" from the file \"%s\"", field_IO_name, restart_read_data_file_name);
    field_instance->check_field_sum(hint);
    field_instance->define_field_values(false);
    field_instance->reset_checksum();
}


void Restart_mgt::get_file_name_in_rpointer_file(char *restart_file_name)
{
    char rpointer_file_name[NAME_STR_SIZE], line[NAME_STR_SIZE*16], *line_p;
    FILE *rpointer_file;


    sprintf(rpointer_file_name, "%s/rpointer.%s", comp_comm_group_mgt_mgr->get_restart_common_dir(), comp_node->get_full_name());
    EXECUTION_REPORT(REPORT_ERROR, comp_node->get_comp_id(), does_file_exist(rpointer_file_name), "Error happens when try to restart a continue/branch run: file \"%s\" does not exist", rpointer_file_name);
    rpointer_file = fopen(rpointer_file_name, "r");
    get_next_line(line, rpointer_file);
    line_p = line;
    get_next_attr(restart_file_name, &line_p);
    fclose(rpointer_file);
}


void Restart_mgt::write_restart_mgt_into_file()
{
    char *array_buffer = NULL;
    long buffer_max_size, buffer_content_size;
    int temp_int;
    char restart_file_name[NAME_STR_SIZE], prev_rpointer_file_name[NAME_STR_SIZE], rpointer_file_name[NAME_STR_SIZE], line[NAME_STR_SIZE*16];
    FILE *restart_file, *rpointer_file;
    

    if (comp_node->get_current_proc_local_id() != 0) {
        clean(true);
        return;
    }

    if (restart_mgt_info_written)
        return;

    if (inout_interface_mgr->is_comp_in_restart_write_window(comp_node->get_comp_id()))
        return;

    restart_mgt_info_written = true;

    temp_int = bypass_import_fields_at_write? 1 : 0;
    write_data_into_array_buffer(&temp_int, sizeof(int), &array_buffer, buffer_max_size, buffer_content_size);
    
    for (int i = restart_write_buffer_containers.size()-1; i >= 0; i --)
        restart_write_buffer_containers[i]->dump_out(&array_buffer, buffer_max_size, buffer_content_size);
    temp_int = restart_write_buffer_containers.size();
    write_data_into_array_buffer(&temp_int, sizeof(int), &array_buffer, buffer_max_size, buffer_content_size);

    int date = last_restart_write_full_time/(long)100000;
    int second = last_restart_write_full_time%(long)100000;
    sprintf(restart_file_name, "%s/%s.%s.r.%08d-%05d", comp_node->get_working_dir(), time_mgr->get_case_name(), comp_node->get_comp_full_name(), date, second);
    restart_file = fopen(restart_file_name, "w+");
    EXECUTION_REPORT(REPORT_ERROR, comp_node->get_comp_id(), restart_file != NULL, "Failed to open the file \"%s\" for writing restart data", restart_file_name);
    fwrite(array_buffer, buffer_content_size, 1, restart_file);
    fclose(restart_file);
    delete [] array_buffer;    
    sprintf(rpointer_file_name, "%s/rpointer.%s", comp_comm_group_mgt_mgr->get_restart_common_dir(), comp_node->get_full_name());
    if (does_file_exist(rpointer_file_name)) {
        rpointer_file = fopen(rpointer_file_name, "r");
        get_next_line(line, rpointer_file);
        fclose(rpointer_file);        
        sprintf(prev_rpointer_file_name, "%s/prev.rpointer.%s", comp_comm_group_mgt_mgr->get_restart_common_dir(), comp_node->get_full_name());
        FILE *prev_rpointer_file = fopen(prev_rpointer_file_name, "w+");
        fprintf(prev_rpointer_file, "%s\n", line);
        fclose(prev_rpointer_file);
    }
    rpointer_file = fopen(rpointer_file_name, "w+");    
    fprintf(rpointer_file, "%s.%s.r.%08d-%05d\n", time_mgr->get_case_name(), comp_node->get_comp_full_name(), date, second);
    fclose(rpointer_file);
    EXECUTION_REPORT_LOG(REPORT_LOG, comp_node->get_comp_id(), true, "Write restart mgt information into the file \"%s\"", restart_file_name);
    EXECUTION_REPORT(REPORT_ERROR, -1, restart_write_data_file != NULL, "Software error in Restart_mgt::write_restart_mgt_into_file");

    backup_restart_write_data_file = restart_write_data_file;
    restart_write_data_file = NULL;

    clean(true);
}


const char *Restart_mgt::get_input_restart_mgt_info_file()
{
    EXECUTION_REPORT(REPORT_ERROR, comp_node->get_comp_id(), input_restart_mgt_info_file != NULL, "Failed to run the model in a continue or branch run: the restart file has not been read in. Please make sure that the API \"CCPL_start_restart_read_IO\" has been called");
    return input_restart_mgt_info_file;
}


const char *Restart_mgt::get_restart_read_annotation()
{
    EXECUTION_REPORT(REPORT_ERROR, comp_node->get_comp_id(), restart_read_annotation != NULL, "Failed to run the model in a continue or branch run: the restart file has not been read in. Please make sure that the API \"CCPL_start_restart_read_IO\" has been called");
    return restart_read_annotation;
}


Restart_buffer_container *Restart_mgt::apply_restart_buffer(const char *comp_full_name, const char *buf_type, const char *keyword)
{
    Restart_buffer_container *new_restart_buffer = new Restart_buffer_container(comp_full_name, buf_type, keyword, this);
    restart_write_buffer_containers.push_back(new_restart_buffer);

    return new_restart_buffer;
}


bool Restart_mgt::is_in_restart_write_window(long full_time, bool is_in_interface_execution)
{
    if (last_restart_write_elapsed_time == -1)
        return false;

    if (full_time == -1)
        return time_mgr->get_current_num_elapsed_day()*((long)100000)+time_mgr->get_current_second() <= last_restart_write_elapsed_time;
    else if (is_in_interface_execution) 
        return full_time <= last_restart_write_elapsed_time || time_mgr->get_current_num_elapsed_day()*((long)100000)+time_mgr->get_current_second() <= last_restart_write_elapsed_time;
    else return full_time < last_restart_write_elapsed_time || time_mgr->get_current_num_elapsed_day()*((long)100000)+time_mgr->get_current_second() <= last_restart_write_elapsed_time;
}


bool Restart_mgt::is_in_restart_read_window(long full_time)
{
    if (time_mgr == NULL)
        time_mgr = components_time_mgrs->get_time_mgr(comp_node->get_comp_id());

    if (time_mgr->get_runtype_mark() == RUNTYPE_MARK_INITIAL || time_mgr->get_runtype_mark() == RUNTYPE_MARK_HYBRID)
        return false;

    if (full_time == -1)
        return false;

    return full_time <= time_mgr->get_restart_full_time() || time_mgr->get_current_num_elapsed_day()*((long)100000)+time_mgr->get_current_second() <= time_mgr->get_restart_full_time();
}


void Restart_mgt::add_restarted_field_instance(Field_mem_info *field_instance, bool is_imported_field)
{
    for (int i = 0; i < restarted_field_instances.size(); i ++)
        if (restarted_field_instances[i].first == field_instance) {
            if (!is_imported_field)
                restarted_field_instances[i].second = is_imported_field;
            return;
        }

    restarted_field_instances.push_back(std::make_pair(field_instance, is_imported_field));
}


bool Restart_mgt::check_restart_read_started()
{
    return input_restart_mgt_info_file != NULL;
}

