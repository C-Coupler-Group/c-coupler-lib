/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#include <mpi.h>
#include "memory_mgt.h"
#include "global_data.h"
#include "cor_global_data.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>


Field_mem_info::Field_mem_info(const char *field_name, int decomp_id, int comp_or_grid_id, 
                               int buf_mark, const char *unit, const char *data_type, const char *annotation, bool check_field_name)
{
    int mem_size;
    Remap_grid_class *remap_grid_grid = NULL, *remap_grid_decomp = NULL;
    Remap_data_field *remap_data_field;


    if (decomp_id == -1) {
        if (comp_comm_group_mgt_mgr->is_legal_local_comp_id(comp_or_grid_id,false)) {
            comp_id = comp_or_grid_id;
            grid_id = -1;
            mem_size = get_data_type_size(data_type);
        }
        else {
            comp_id = original_grid_mgr->get_comp_id_of_grid(comp_or_grid_id);
            grid_id = comp_or_grid_id;
            EXECUTION_REPORT(REPORT_ERROR, comp_id, original_grid_mgr->get_original_CoR_grid(grid_id)->get_num_dimensions() == 1 && original_grid_mgr->get_original_CoR_grid(grid_id)->has_grid_coord_label(COORD_LABEL_LEV), 
                             "Error happens when calling the API \"CCPL_register_field_instance\" to register a field instance of \"%s\": when the given parallel decomposition ID (the parameter \"decomp_id\") is -1, the corresponding grid \"%s\" must be an one-dimension vertical grid. Please check the model code with the annotation \"%s\"", field_name, original_grid_mgr->get_original_grid(comp_or_grid_id)->get_grid_name(), annotation);
            mem_size = original_grid_mgr->get_grid_size(grid_id, "in Field_mem_info::Field_mem_info") * get_data_type_size(data_type);
        }
        host_comp_id = comp_id;
    }
    else {
        EXECUTION_REPORT(REPORT_ERROR, -1, original_grid_mgr->is_grid_id_legal(comp_or_grid_id), "Software error2 in new Field_mem_info");
        EXECUTION_REPORT(REPORT_ERROR, -1, decomps_info_mgr->is_decomp_id_legal(decomp_id), "Software error3 in new Field_mem_info");
        grid_id = comp_or_grid_id;
        comp_id = original_grid_mgr->get_comp_id_of_grid(comp_or_grid_id);
        host_comp_id = decomps_info_mgr->get_decomp_info(decomp_id)->get_host_comp_id();
        EXECUTION_REPORT(REPORT_ERROR, -1, comp_id == decomps_info_mgr->get_comp_id_of_decomp(decomp_id), 
                         "Software error4 in new Field_mem_info");
        remap_grid_decomp = decomps_info_mgr->get_CoR_grid_of_decomp(decomp_id);
        remap_grid_grid = original_grid_mgr->get_original_CoR_grid(comp_or_grid_id);
        mem_size = decomps_info_mgr->get_decomp_info(decomp_id)->get_num_local_cells() * get_data_type_size(data_type) * remap_grid_grid->get_grid_size()/remap_grid_decomp->get_grid_size();
        EXECUTION_REPORT(REPORT_ERROR, host_comp_id, remap_grid_decomp->is_subset_of_grid(remap_grid_grid), "Error happens when calling the API \"CCPL_register_field_instance\" to register a field instance of \"%s\": the parameters of grid ID and decomposition ID do not match each other: the grid corresponding to the decomposition (grid \"%s\") is not a subset of the grid corresponding to the grid ID (grid \"%s\"). Please check the model code with the annotation \"%s\"", field_name, decomps_info_mgr->get_decomp_info(decomp_id)->get_grid_name(), original_grid_mgr->get_original_grid(comp_or_grid_id)->get_grid_name(), annotation);
    }
    
    host_comp_time_mgr = components_time_mgrs->get_time_mgr(host_comp_id);

    const field_attr *field_attributes = fields_info->search_field_info(field_name);

    if (check_field_name) {
        EXECUTION_REPORT(REPORT_ERROR, host_comp_id, field_attributes != NULL, "Error happens when calling the API \"CCPL_register_field_instance\" to register a field instance of \"%s\": the field name \"%s\" is unknown (has not been registered through the configuration XML file public_field_attribute.xml). Please check the model code with the annotation \"%s\"", field_name, field_name, annotation);
        bool dimensions_match_grid;
        if (words_are_the_same(field_attributes->field_dim, FIELD_0_DIM))
            dimensions_match_grid = decomp_id == -1 && remap_grid_grid == NULL;
        if (words_are_the_same(field_attributes->field_dim, FIELD_V1_DIM))
            dimensions_match_grid = decomp_id == -1 && remap_grid_grid != NULL && remap_grid_grid->get_num_dimensions() == 1 && remap_grid_grid->has_grid_coord_label(COORD_LABEL_LEV);
        else if (words_are_the_same(field_attributes->field_dim, FIELD_2_DIM))
            dimensions_match_grid = decomp_id != -1 && remap_grid_grid != NULL && remap_grid_grid->get_is_sphere_grid();
        else if (words_are_the_same(field_attributes->field_dim, FIELD_3_DIM))
            dimensions_match_grid = decomp_id != -1 && remap_grid_grid != NULL && remap_grid_grid->get_num_dimensions() == 3;
        if (!dimensions_match_grid) {
            if (grid_id != -1)
                EXECUTION_REPORT(REPORT_ERROR, comp_id, false, "Error happens when calling the API \"CCPL_register_field_instance\" to register a field instance of \"%s\" at the model code with the annotation \"%s\": the dimension information (\"%s\") of the field that is specified in a configuration file does not match the dimensions of the corresponding grid \"%s\"", field_name, annotation, field_attributes->field_dim, original_grid_mgr->get_name_of_grid(grid_id));
            else EXECUTION_REPORT(REPORT_ERROR, comp_id, false, "Error happens when calling the API \"CCPL_register_field_instance\" to register a field instance of \"%s\": the dimension information (\"%s\") of the field that is specified in a configuration file does not match the dimensions of the corresponding empty grid", field_name, annotation, field_attributes->field_dim);
        }
    }    

    if (strlen(unit) > 0)
        strcpy(field_unit, unit);
    else if (field_attributes != NULL) 
        strcpy(field_unit, fields_info->search_field_info(field_name)->field_unit);
    // check the field unit

    strcpy(this->field_name, field_name);
    this->decomp_id = decomp_id;
    this->comp_or_grid_id = comp_or_grid_id;
    this->buf_mark = buf_mark;
    this->usage_tag = -1;
    is_registered_model_buf = false;
    is_field_active = false;
    define_order_count = -1;
    last_define_time = 0x7fffffffffffffff;

    remap_data_field = new Remap_data_field;
    strcpy(remap_data_field->field_name_in_application, field_name);
    strcpy(remap_data_field->field_name_in_IO_file, field_name);
    strcpy(remap_data_field->data_type_in_application, data_type);
    strcpy(remap_data_field->data_type_in_IO_file, data_type);
    remap_data_field->required_data_size = mem_size / get_data_type_size(data_type);
    remap_data_field->read_data_size = remap_data_field->required_data_size;
    if (remap_data_field->required_data_size > 0) {
        remap_data_field->data_buf = (char*) (new long [(mem_size+sizeof(long)-1)/sizeof(long)]);
        memset(remap_data_field->data_buf, 0, mem_size);
    }
    else remap_data_field->data_buf = NULL;
    if (check_field_name)
        remap_data_field->set_field_long_name(fields_info->get_field_long_name(field_name));
    remap_data_field->set_field_unit(unit);   // to complete: when strlen(unit) is 0, use default unit of the field

    if (decomp_id == -1)
        grided_field_data = new Remap_grid_data_class(NULL, remap_data_field);
    else {
        Remap_grid_class *decomp_grid = decomp_grids_mgr->search_decomp_grid_info(decomp_id, remap_grid_grid, false)->get_decomp_grid();
        grided_field_data = new Remap_grid_data_class(decomp_grid, remap_data_field);
        remap_data_field->set_fill_value(NULL);
    }
    
    last_checksum = -1;
}


Field_mem_info::~Field_mem_info()
{
    if (is_registered_model_buf)
        grided_field_data->get_grid_data_field()->data_buf = NULL;
    delete grided_field_data;
}


void Field_mem_info::reset_mem_buf(void * buf, bool is_external_field, int usage_tag)
{
	if (get_size_of_field() > 0)
	    EXECUTION_REPORT(REPORT_ERROR, host_comp_id, buf != NULL, "The data buffer corresponding to the field instance of \"%s\" is not allocated. Please verify the model code corresponding to the annotation \"%s\"", field_name, annotation_mgr->get_annotation(field_instance_id, "allocate field instance"));
    EXECUTION_REPORT(REPORT_ERROR, -1, !is_registered_model_buf, "Software error to release a registered buffer");

    if (grided_field_data->get_grid_data_field()->data_buf != NULL)
        delete [] grided_field_data->get_grid_data_field()->data_buf;

    grided_field_data->get_grid_data_field()->data_buf = buf;

    if (is_external_field) {
        this->usage_tag = usage_tag;
		is_registered_model_buf = true;
    }
}


void Field_mem_info::change_datatype_to_double()
{
    grided_field_data->change_datatype_in_application(DATA_TYPE_DOUBLE);
}


void Field_mem_info::define_field_values(bool is_restarting)
{
    if (!is_restarting)
        is_field_active = true;
    last_define_time = host_comp_time_mgr->get_current_full_time();
}


void Field_mem_info::use_field_values(const char *annotation)
{    
    if (is_registered_model_buf) 
        return;
    
    if (last_define_time == host_comp_time_mgr->get_current_full_time())
        return;

    if (last_define_time == 0x7fffffffffffffff) {
        if (is_registered_model_buf)
            EXECUTION_REPORT(REPORT_ERROR, host_comp_id, false, "field instance (field_name=\"%s\", decomp_name=\"%s\", grid_name=\"%s\", bufmark=%x) is used before defining it. Please modify the model code with the annotation \"%s\"", field_name, get_decomp_name(), get_grid_name(), buf_mark, annotation);
        else EXECUTION_REPORT(REPORT_ERROR, -1, false, "Software error in Field_mem_info::use_field_values: field instance (field_name=\"%s\", decomp_name=\"%s\", grid_name=\"%s\", bufmark=%x) is used before defining it. Please modify the model code with the annotation \"%s\"", field_name, get_decomp_name(), get_grid_name(), buf_mark, annotation);        
    }
    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, host_comp_id, last_define_time <= host_comp_time_mgr->get_current_full_time(), "C-Coupler error in Field_mem_info::use_field_values: wrong time order");
}


bool Field_mem_info::match_field_instance(const char *field_name, int decomp_id, int comp_or_grid_id, int buf_mark)
{
    return words_are_the_same(this->field_name, field_name) && this->decomp_id == decomp_id && this->comp_or_grid_id == comp_or_grid_id && this->buf_mark == buf_mark;
}


bool Field_mem_info::match_field_mem(void *data_buffer)
{
    return this->get_data_buf() == data_buffer;
}


void Field_mem_info::reset_field_name(const char *new_name)
{
    strcpy(field_name, new_name);
}


void Field_mem_info::calculate_field_conservative_sum(Field_mem_info *area_field)
{
    double partial_sum, total_sum;
    long size;

    if (report_internal_log_enabled) {
        EXECUTION_REPORT(REPORT_ERROR,-1, words_are_the_same(get_field_data()->get_grid_data_field()->data_type_in_application, DATA_TYPE_DOUBLE), "C-Coupler error in calculate_field_sum");
        size = get_field_data()->get_grid_data_field()->required_data_size;
        partial_sum = 0;
        for (long j = 0; j < size; j ++)
            partial_sum += (((double*) get_data_buf())[j])*(((double*) area_field->get_data_buf())[j]);
        MPI_Allreduce(&partial_sum, &total_sum, 1, MPI_DOUBLE, MPI_SUM, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(area_field->get_comp_id(),"in Field_mem_info::calculate_field_conservative_sum"));
    }
}


void Field_mem_info::check_field_sum(const char *hint)
{
    int partial_sum, total_sum;
    long size;


    if (report_error_enabled && is_registered_model_buf) {
        MPI_Barrier(comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(host_comp_id, "Field_mem_info::check_field_sum"));
        EXECUTION_REPORT_LOG(REPORT_LOG, host_comp_id, true, "Try to check the model data buffer of the field \"%s\" registered corresponding to the code annotation \"%s\". If it fails to pass the check (the model run is stopped), please make sure corresponding model data buffer is a global variable and has not been released", field_name, annotation_mgr->get_annotation(field_instance_id, "allocate field instance"));
        char *temp_array = new char [get_size_of_field()*get_data_type_size(get_field_data()->get_grid_data_field()->data_type_in_application)];
        memcpy(temp_array, get_data_buf(), get_size_of_field()*get_data_type_size(get_field_data()->get_grid_data_field()->data_type_in_application));
        memcpy(get_data_buf(), temp_array, get_size_of_field()*get_data_type_size(get_field_data()->get_grid_data_field()->data_type_in_application));
        MPI_Barrier(comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(host_comp_id, "Field_mem_info::check_field_sum"));        
        EXECUTION_REPORT_LOG(REPORT_LOG, host_comp_id, true, "Pass the check of the model data buffer of the field \"%s\" registered corresponding to the code annotation \"%s\". If it fails to pass the check (the model run is stopped), please make sure corresponding model data buffer is a global variable and has not been released", field_name, annotation_mgr->get_annotation(field_instance_id, "allocate field instance"));
        delete [] temp_array;
    }

    if (report_internal_log_enabled) {
        size = get_data_type_size(get_field_data()->get_grid_data_field()->data_type_in_application)*get_field_data()->get_grid_data_field()->required_data_size/4;
        partial_sum = 0;
        for (long j = 0; j < size; j ++)
            partial_sum += (((int*) get_data_buf())[j]);

        if (decomp_id != -1) {
            MPI_Allreduce(&partial_sum, &total_sum, 1, MPI_INT, MPI_SUM, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(host_comp_id, "Field_mem_info::check_field_sum"));
            EXECUTION_REPORT_LOG(REPORT_LOG, host_comp_id, true, "check sum of field \"%s\" %s is %x", get_field_name(), hint, total_sum);
//            EXECUTION_REPORT_LOG(REPORT_LOG, host_comp_id, true, "check sum of field \"%s\" %s is %x vs %x", get_field_name(), hint, total_sum, partial_sum);
        }
        else {
            total_sum = partial_sum;
            MPI_Bcast(&total_sum, 1, MPI_INT, 0, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(host_comp_id, "Field_mem_info::check_field_sum"));
            EXECUTION_REPORT(REPORT_ERROR, host_comp_id, partial_sum == total_sum, "As an instance of the field \"%s\" is not on a horizontal grid, all its values should be the same but currently are not the same across all processes of the corresponding component model. Please check the model code related to the annotation \"%s\"", field_name, annotation_mgr->get_annotation(field_instance_id, "allocate field instance"));
        }
    }
}


bool Field_mem_info::field_has_been_defined()
{
    return last_define_time != 0x7fffffffffffffff;
}


long Field_mem_info::get_size_of_field()
{
    return grided_field_data->get_grid_data_field()->required_data_size;
}


const char *Field_mem_info::get_grid_name() 
{
    if (grid_id == -1)
        return NULL;

    return original_grid_mgr->search_grid_info(grid_id)->get_grid_name();
}


const char *Field_mem_info::get_decomp_name()
{
    if (decomp_id == -1)
        return NULL;

    return decomps_info_mgr->get_decomp_info(decomp_id)->get_decomp_name();
}


const char *Field_mem_info::get_data_type()
{
    return get_field_data()->get_grid_data_field()->data_type_in_application;
}


void Field_mem_info::set_field_instance_id(int field_instance_id, const char *annotation)
{
    this->field_instance_id = field_instance_id;
    annotation_mgr->add_annotation(field_instance_id, "allocate field instance", annotation);
}


bool Field_mem_info::is_checksum_changed()
{
    if (last_checksum == -1)
        return true;

    long current_checksum = calculate_checksum_of_array((const char*)get_data_buf(), get_size_of_field(), get_data_type_size(get_data_type()), NULL, NULL);
    
    return current_checksum != last_checksum;
}


void Field_mem_info::reset_checksum()
{
    last_checksum = calculate_checksum_of_array((const char*)get_data_buf(), get_size_of_field(), get_data_type_size(get_data_type()), NULL, NULL);
}


Field_mem_info *Memory_mgt::alloc_mem(Field_mem_info *original_field_instance, int special_buf_mark, int object_id, const char *unit_or_datatype, bool check_field_name)
{
    EXECUTION_REPORT(REPORT_ERROR, -1, special_buf_mark == BUF_MARK_DATATYPE_TRANS || special_buf_mark == BUF_MARK_AVERAGED_INNER || special_buf_mark == BUF_MARK_AVERAGED_INTER || special_buf_mark == BUF_MARK_UNIT_TRANS || special_buf_mark == BUF_MARK_DATA_TRANSFER || 
                     special_buf_mark == BUF_MARK_IO_FIELD_MIRROR || special_buf_mark == BUF_MARK_REMAP_NORMAL || special_buf_mark == BUF_MARK_REMAP_DATATYPE_TRANS_SRC || special_buf_mark == BUF_MARK_REMAP_DATATYPE_TRANS_DST || special_buf_mark == BUF_MARK_REMAP_FRAC, "Software error in Field_mem_info *alloc_mem: wrong special_buf_mark");
    int new_buf_mark = (special_buf_mark ^ object_id);
    Field_mem_info *existing_field_instance = search_field_instance(original_field_instance->get_field_name(), original_field_instance->get_decomp_id(), original_field_instance->get_comp_or_grid_id(), new_buf_mark);
    if (existing_field_instance != NULL) {
        if (special_buf_mark == BUF_MARK_UNIT_TRANS)
            EXECUTION_REPORT(REPORT_ERROR, -1, words_are_the_same(existing_field_instance->get_data_type(), original_field_instance->get_data_type()) && words_are_the_same(existing_field_instance->get_unit(), unit_or_datatype), "Software error in Field_mem_info *alloc_mem: special field instance exists %lx with different data type or wrong unit", new_buf_mark);
        else if (unit_or_datatype != NULL)
            EXECUTION_REPORT(REPORT_ERROR, -1, words_are_the_same(existing_field_instance->get_data_type(), unit_or_datatype), "Software error in Field_mem_info *alloc_mem: special field instance exists %lx with wrong data type", new_buf_mark);
        else EXECUTION_REPORT(REPORT_ERROR, -1, words_are_the_same(existing_field_instance->get_data_type(), original_field_instance->get_data_type()), "Software error in Field_mem_info *alloc_mem: special field instance exists %lx with wrong data type", new_buf_mark);
        return existing_field_instance;
    }
    if (special_buf_mark == BUF_MARK_AVERAGED_INNER || special_buf_mark == BUF_MARK_AVERAGED_INTER)
        fields_mem.push_back(new Field_mem_info(original_field_instance->get_field_name(), original_field_instance->get_decomp_id(), original_field_instance->get_comp_or_grid_id(), new_buf_mark, original_field_instance->get_unit(), original_field_instance->get_data_type(), "new field instance for averaging", check_field_name));    
    else if (special_buf_mark == BUF_MARK_REMAP_FRAC)
        fields_mem.push_back(new Field_mem_info(original_field_instance->get_field_name(), original_field_instance->get_decomp_id(), original_field_instance->get_comp_or_grid_id(), new_buf_mark, original_field_instance->get_unit(), original_field_instance->get_data_type(), "new field instance for the remapping with fraction", check_field_name));
    else if (special_buf_mark == BUF_MARK_DATATYPE_TRANS || special_buf_mark == BUF_MARK_DATA_TRANSFER || special_buf_mark == BUF_MARK_REMAP_DATATYPE_TRANS_SRC || special_buf_mark == BUF_MARK_REMAP_DATATYPE_TRANS_DST) {
        get_data_type_size(unit_or_datatype);
        fields_mem.push_back(new Field_mem_info(original_field_instance->get_field_name(), original_field_instance->get_decomp_id(), original_field_instance->get_comp_or_grid_id(), new_buf_mark, original_field_instance->get_unit(), unit_or_datatype, "new field instance for data type transformation", check_field_name));
    }
    else if (special_buf_mark == BUF_MARK_IO_FIELD_MIRROR) {
        get_data_type_size(unit_or_datatype);
        fields_mem.push_back(new Field_mem_info(original_field_instance->get_field_name(), original_field_instance->get_decomp_id(), original_field_instance->get_comp_or_grid_id(), new_buf_mark, original_field_instance->get_unit(), unit_or_datatype, "new field instance for data type transformation", check_field_name));        
    }
    else if (special_buf_mark == BUF_MARK_UNIT_TRANS) {
        // check unit
        fields_mem.push_back(new Field_mem_info(original_field_instance->get_field_name(), original_field_instance->get_decomp_id(), original_field_instance->get_comp_or_grid_id(), new_buf_mark, unit_or_datatype, original_field_instance->get_data_type(), "new field instance for unit transformation", check_field_name));
    }
    else if (special_buf_mark == BUF_MARK_REMAP_NORMAL) {
        // check unit
        fields_mem.push_back(new Field_mem_info(original_field_instance->get_field_name(), original_field_instance->get_decomp_id(), original_field_instance->get_comp_or_grid_id(), new_buf_mark, original_field_instance->get_unit(), unit_or_datatype, "new field instance for remapping", check_field_name));
    }
    else if (special_buf_mark == BUF_MARK_REMAP_DATATYPE_TRANS_SRC) {
        // check unit
        fields_mem.push_back(new Field_mem_info(original_field_instance->get_field_name(), original_field_instance->get_decomp_id(), original_field_instance->get_comp_or_grid_id(), new_buf_mark, original_field_instance->get_unit(), unit_or_datatype, "new field instance for data type transformation in remapping", check_field_name));
    }
    else if (special_buf_mark == BUF_MARK_REMAP_DATATYPE_TRANS_DST) {
        // check unit
        fields_mem.push_back(new Field_mem_info(original_field_instance->get_field_name(), original_field_instance->get_decomp_id(), original_field_instance->get_comp_or_grid_id(), new_buf_mark, original_field_instance->get_unit(), unit_or_datatype, "new field instance for data type transformation in remapping", check_field_name));
    }
    else EXECUTION_REPORT(REPORT_ERROR, -1, false, "Software error in Field_mem_info *alloc_mem");

    fields_mem[fields_mem.size()-1]->set_field_instance_id(TYPE_FIELD_INST_ID_PREFIX|(fields_mem.size()-1), "in Memory_mgt::alloc_mem");

    return fields_mem[fields_mem.size()-1];
}


Field_mem_info *Memory_mgt::alloc_mem(const char *field_name, int decomp_id, int comp_or_grid_id, int buf_mark, const char *data_type, const char *field_unit, const char *annotation, bool check_field_name)
{
    Field_mem_info *field_mem, *pair_field;
    int i, comp_id;
    bool find_field_in_cfg;


    EXECUTION_REPORT(REPORT_ERROR, -1, buf_mark < 0, "Software error in Memory_mgt::alloc_mem: wrong value of buffer mark");
    EXECUTION_REPORT(REPORT_ERROR, -1, data_type != NULL, "Software error in Memory_mgt::alloc_mem: data type is NULL");
    get_data_type_size(data_type);
    if (decomp_id != -1) {
        EXECUTION_REPORT(REPORT_ERROR, -1, decomps_info_mgr->is_decomp_id_legal(decomp_id), "Software error in Memory_mgt::alloc_mem: wrong decomposition id");
        EXECUTION_REPORT(REPORT_ERROR, -1, original_grid_mgr->is_grid_id_legal(comp_or_grid_id), "Software error in Memory_mgt::alloc_mem: wrong grid id");
        comp_id = original_grid_mgr->search_grid_info(comp_or_grid_id)->get_comp_id();
    }
    else {
        EXECUTION_REPORT(REPORT_ERROR, -1, comp_comm_group_mgt_mgr->is_legal_local_comp_id(comp_or_grid_id,false) || original_grid_mgr->is_grid_id_legal(comp_or_grid_id), "Software error in Memory_mgt::alloc_mem: wrong component id or grid_id");
		if (original_grid_mgr->is_grid_id_legal(comp_or_grid_id))
			comp_id = original_grid_mgr->search_grid_info(comp_or_grid_id)->get_comp_id();
		else comp_id = comp_or_grid_id;
    }

    
    /* If memory buffer has been allocated, return it */
    for (i = 0; i < fields_mem.size(); i ++)
        if (fields_mem[i]->match_field_instance(field_name, decomp_id, comp_or_grid_id, buf_mark)) {
            // EXECUTION_REPORT(REPORT_ERROR, comp_id, field_unitÒ»ÖÂ, ...);
            EXECUTION_REPORT(REPORT_ERROR, -1, words_are_the_same(data_type, fields_mem[i]->get_field_data()->get_grid_data_field()->data_type_in_application),
                             "Software error in Memory_mgt::alloc_mem: data types conflict");
            return fields_mem[i];
        }

    /* Compute the size of the memory buffer and then allocate and return it */
    field_mem = new Field_mem_info(field_name, decomp_id, comp_or_grid_id, buf_mark, field_unit, data_type, annotation, check_field_name);
    field_mem->set_field_instance_id(TYPE_FIELD_INST_ID_PREFIX|fields_mem.size(), annotation);
    fields_mem.push_back(field_mem);

    return field_mem;
}


Memory_mgt::~Memory_mgt()
{
    for (int i = 0; i < fields_mem.size(); i ++)
        delete fields_mem[i];
}


Field_mem_info *Memory_mgt::search_field_via_data_buf(const void *data_buf, bool diag)
{
    for (int i = 0; i < fields_mem.size(); i ++)
        if (fields_mem[i]->get_data_buf() == data_buf)
            return fields_mem[i];

    if (diag)
        EXECUTION_REPORT(REPORT_ERROR,-1, false, "C-Coupler error in search_field_via_data_buf\n");
    return NULL;
}


void Memory_mgt::check_sum_of_all_fields()
{
    for (int i = 0; i < fields_mem.size(); i ++)
        fields_mem[i]->check_field_sum("when checking all fields");
}


int Memory_mgt::get_field_size(void *data_buf, const char *annotation)
{
    Field_mem_info *field = search_field_via_data_buf(data_buf, false);

    EXECUTION_REPORT(REPORT_ERROR,-1, field != NULL, "Detect a memory buffer that is not managed by C-Coupler. Please verify the model code according to annotation \"%s\"", annotation);

    return field->get_size_of_field();
}


Field_mem_info *Memory_mgt::search_field_instance(const char *field_name, int decomp_id, int comp_or_grid_id, int buf_mark)
{
    for (int i = 0; i < fields_mem.size(); i ++)
        if (fields_mem[i]->match_field_instance(field_name, decomp_id, comp_or_grid_id, buf_mark))
            return fields_mem[i];

    return NULL;
}


int Memory_mgt::register_external_field_instance(const char *field_name, void *data_buffer, int field_size, int decomp_id, int comp_or_grid_id, 
                                                 int buf_mark, int usage_tag, const char *unit, const char *data_type, const char *annotation)
{
    int comp_id, API_id;
    Field_mem_info *existing_field_instance_instance, *new_field_instance;


    if (buf_mark == BUF_MARK_IO_FIELD_REG)
        API_id = API_ID_FIELD_MGT_REG_IO_FIELD_from_BUFFER;
    else API_id = API_ID_FIELD_MGT_REG_FIELD_INST;

    EXECUTION_REPORT(REPORT_ERROR, -1, comp_comm_group_mgt_mgr->is_legal_local_comp_id(comp_or_grid_id,true) || original_grid_mgr->is_grid_id_legal(comp_or_grid_id), "Error happens when calling the API \"CCPL_register_field_instance\" to register a field instance of \"%s\": the parameter of \"comp_or_grid_id\" is not a grid id or a component id. Please check the model code with the annotation \"%s\"", field_name, annotation);
    
    if (comp_comm_group_mgt_mgr->is_legal_local_comp_id(comp_or_grid_id,true))
        comp_id = comp_or_grid_id;
    else comp_id = original_grid_mgr->get_comp_id_of_grid(comp_or_grid_id);

    check_API_parameter_string_length(comp_id, API_id, CCPL_NAME_STR_LEN, field_name, "field_name", annotation);
    check_and_verify_name_format_of_string_for_API(comp_id, field_name, API_id, "the field instance", annotation);

    if (decomp_id != -1) {
        EXECUTION_REPORT(REPORT_ERROR, -1, decomps_info_mgr->is_decomp_id_legal(decomp_id), "Error happens when calling the API \"CCPL_register_field_instance\" to register a field instance of \"%s\": the parameter of decomposition ID is wrong. Please check the model code with the annotation \"%s\"", field_name, annotation);
        EXECUTION_REPORT(REPORT_ERROR, -1, original_grid_mgr->is_grid_id_legal(comp_or_grid_id), "Error happens when calling the API \"CCPL_register_field_instance\" to register a field instance of \"%s\": the parameter of grid ID is wrong. Please check the model code with the annotation \"%s\"", field_name, annotation);        
        EXECUTION_REPORT(REPORT_ERROR, -1, comp_id == decomps_info_mgr->get_comp_id_of_decomp(decomp_id), 
                         "Error happens when calling the API \"CCPL_register_field_instance\" to register a field instance of \"%s\": the parameters of grid ID and decomposition ID do not match each other: they belong to different component models. Please check the model code with the annotation \"%s\"",
                         field_name, annotation);
        EXECUTION_REPORT(REPORT_ERROR, -1, comp_comm_group_mgt_mgr->is_legal_local_comp_id(comp_id,true), "Software error in Memory_mgt::register_external_field_instance: illegal component id from grid id");
    }
    synchronize_comp_processes_for_API(comp_id, API_id, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id,"C-Coupler code in register_external_field_instance for getting component management node"), "registering a field instance or a I/O field", annotation);
    comp_comm_group_mgt_mgr->confirm_coupling_configuration_active(comp_id, API_id, true, annotation);
    check_API_parameter_string(comp_id, API_id, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id,"C-Coupler code in register_external_field_instance for getting component management node"), "registering a field instance or a I/O field", field_name, "field name", annotation);
    check_API_parameter_string(comp_id, API_id, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id,"C-Coupler code in register_external_field_instance for getting component management node"), "registering a field instance or a I/O field", comp_comm_group_mgt_mgr->get_global_node_of_local_comp(comp_id,true,"C-Coupler code in register_external_field_instance for getting component management node")->get_comp_name(), "the component name specified by the corresponding ID", annotation);
    check_API_parameter_int(comp_id, API_id, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id,"C-Coupler code in register_external_field_instance for getting component management node"), NULL, buf_mark, "buf_mark", annotation);
    check_API_parameter_int(comp_id, API_id, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id,"C-Coupler code in register_external_field_instance for getting component management node"), NULL, usage_tag, "usage_tag", annotation);
    check_API_parameter_string(comp_id, API_id, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id,"C-Coupler code in register_external_field_instance for getting component management node"), "registering a field instance or a I/O field", data_type, "the data type (such as integer, float, and double) of the field instance", annotation);
    check_API_parameter_string(comp_id, API_id, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id,"C-Coupler code in register_external_field_instance for getting component management node"), "registering a field instance or a I/O field", unit, "unit", annotation);
    if (decomp_id == -1) {
        check_API_parameter_int(comp_id, API_id, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id,"C-Coupler code in register_external_field_instance for getting component management node"), NULL, decomp_id, "decomp_id (the ID of the parallel decomposition)", annotation);
        int temp_int = (comp_id == comp_or_grid_id)? 1 : 0;
        check_API_parameter_int(comp_id, API_id, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id,"C-Coupler code in register_external_field_instance for getting component management node"), NULL, decomp_id, "comp_or_grid_id (a grid id or a component id)", annotation);
        if (comp_id != comp_or_grid_id)    {
            check_API_parameter_string(comp_id, API_id, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id,"C-Coupler code in register_external_field_instance for getting component management node"), "registering a field instance or a I/O field", original_grid_mgr->get_name_of_grid(comp_or_grid_id), "the grid name specified by the corresponding ID", annotation);
            EXECUTION_REPORT(REPORT_ERROR, comp_id, original_grid_mgr->get_original_grid(comp_or_grid_id)->get_H2D_sub_CoR_grid() == NULL && original_grid_mgr->get_original_grid(comp_or_grid_id)->get_V1D_sub_CoR_grid() != NULL, "Error happens when calling the API \"CCPL_register_field_instance\" to register a field instance of \"%s\": the grid corresponding to the parameter of \"comp_or_grid_id\" (the grid is \"%s\") should be but not a vertical grid when the given \"decomp_id\" is -1. Please check the model code with the annotation \"%s\"", field_name, original_grid_mgr->get_original_grid(comp_or_grid_id)->get_grid_name(), annotation);
        }
        EXECUTION_REPORT(REPORT_ERROR, -1, comp_id == comp_or_grid_id, "Error happens when calling the API \"CCPL_register_field_instance\" to register a field instance of \"%s\" at the model code with the annotation \"%s\". We are sorry that C-Coupler now only supports the coupling of a scalar field or a field on a grid related to a horizontal grid that is decomposed in parallelization of a model. If you want to couple more kinds of fields, please contact us (liuli-cess@tsinghua.edu.cn)", field_name, annotation);
    }
    else {
        check_API_parameter_string(comp_id, API_id, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id,"C-Coupler code in register_external_field_instance for getting component management node"), "registering a field instance or a I/O field", decomps_info_mgr->get_decomp_info(decomp_id)->get_decomp_name(), "the parallel decomposition name specified by the corresponding ID", annotation);
        check_API_parameter_string(comp_id, API_id, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id,"C-Coupler code in register_external_field_instance for getting component management node"), "registering a field instance or a I/O field", original_grid_mgr->get_name_of_grid(comp_or_grid_id), "the grid name specified by the corresponding ID", annotation);
    }

    if (buf_mark != BUF_MARK_IO_FIELD_REG)
        EXECUTION_REPORT(REPORT_ERROR, comp_id, buf_mark >= 0, "Error happens when calling the API \"CCPL_register_field_instance\" to register a field instance of \"%s\": the parameter of the mark (\"buf_mark\") of the field instance cannot be a negative integer (currently is %d). Please check the model code with the annotation \"%s\"",
                         field_name, buf_mark, annotation);

    existing_field_instance_instance = search_field_instance(field_name, decomp_id, comp_or_grid_id, buf_mark);
    if (existing_field_instance_instance != NULL)
        EXECUTION_REPORT(REPORT_ERROR, comp_id, false, "Error happens when calling the API \"CCPL_register_field_instance\" to register a field instance: cannot register an instance of the field of \"%s\" again (the corresponding annotation is \"%s\") because this field instance has been registered before (the corresponding annotation is \"%s\")", 
                         field_name, annotation, annotation_mgr->get_annotation(existing_field_instance_instance->get_field_instance_id(), "allocate field instance"));

    new_field_instance = new Field_mem_info(field_name, decomp_id, comp_or_grid_id, buf_mark, unit, data_type, annotation, (buf_mark!=BUF_MARK_IO_FIELD_REG) && (usage_tag&REG_FIELD_TAG_CPL) == REG_FIELD_TAG_CPL);
    if (new_field_instance->get_size_of_field() > 0)
        EXECUTION_REPORT(REPORT_ERROR, comp_id, field_size == new_field_instance->get_size_of_field(), "Error happens when calling the API \"CCPL_register_field_instance\" to register a field instance of \"%s\": the size of the model data buffer (currently is %d) is different from the size determined by the parallel decomposition and grid (currently is %ld). Please check the model code with the annotation \"%s\"",
    	                 field_name, field_size, new_field_instance->get_size_of_field(), annotation);
    new_field_instance->set_field_instance_id(TYPE_FIELD_INST_ID_PREFIX|fields_mem.size(), annotation);
    new_field_instance->reset_mem_buf(data_buffer, true, usage_tag);
    EXECUTION_REPORT(REPORT_ERROR, comp_id, usage_tag >= 0 && usage_tag <= 3, "Error happens when calling the API \"CCPL_register_field_instance\" to register a field instance of \"%s\": the value of the parameter \"usage_tag\" (%d) is wrong. The right value should be between 1 and 3. Please check the model code with the annotation \"%s\"", field_name, usage_tag, annotation);
    fields_mem.push_back(new_field_instance);

    return new_field_instance->get_field_instance_id();
}


bool Memory_mgt::check_is_legal_field_instance_id(int field_instance_id)
{
    if ((field_instance_id&TYPE_ID_PREFIX_MASK) != TYPE_FIELD_INST_ID_PREFIX)
        return false;

    return (field_instance_id&TYPE_ID_SUFFIX_MASK) < fields_mem.size();
}


Field_mem_info *Memory_mgt::get_field_instance(int field_instance_id)
{
    if (!check_is_legal_field_instance_id(field_instance_id))
        return NULL;

    return fields_mem[field_instance_id&TYPE_ID_SUFFIX_MASK];
}


void Memory_mgt::copy_field_data_values(Field_mem_info *dst_field_inst, Field_mem_info *src_field_inst)
{
    EXECUTION_REPORT(REPORT_ERROR, -1, dst_field_inst->get_decomp_id() == src_field_inst->get_decomp_id() && dst_field_inst->get_comp_id() == src_field_inst->get_comp_id() && dst_field_inst->get_grid_id() == src_field_inst->get_grid_id() && words_are_the_same(src_field_inst->get_data_type(), dst_field_inst->get_data_type()),
                     "Software erorr in Memory_mgt::copy_field_data_values");
    memcpy(dst_field_inst->get_data_buf(), src_field_inst->get_data_buf(), dst_field_inst->get_size_of_field()*get_data_type_size(dst_field_inst->get_data_type()));
}

