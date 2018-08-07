/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#include "global_data.h"
#include "common_utils.h"
#include <string.h>
#include <stdlib.h>
#include <math.h>


void write_string_into_array_buffer(const char *full_string, long array_size, char **temp_array_buffer, long &buffer_max_size, long &buffer_content_size)
{
    write_data_into_array_buffer(full_string, array_size, temp_array_buffer, buffer_max_size, buffer_content_size);
    memset((*temp_array_buffer)+buffer_content_size-array_size+strlen(full_string), 0, array_size-strlen(full_string));
}


void write_data_into_array_buffer(const void *data, long data_size, char **temp_array_buffer, long &buffer_max_size, long &buffer_content_size)
{
    if (*temp_array_buffer == NULL) {
        buffer_max_size = 2 * data_size;
        buffer_content_size = 0;
        *temp_array_buffer = new char [buffer_max_size];
    }

    if (buffer_max_size < buffer_content_size+data_size) {
        buffer_max_size = (buffer_content_size+data_size) * 2;
        char *temp_buffer = new char [buffer_max_size];
        for (int i = 0; i < buffer_content_size; i ++)
            temp_buffer[i] = (*temp_array_buffer)[i];
        delete [] *temp_array_buffer;
        *temp_array_buffer = temp_buffer;
    }

    for (int i = 0; i < data_size; i ++)
        (*temp_array_buffer)[buffer_content_size++] = ((char*)data)[i];
}


bool read_data_from_array_buffer(void *data, long data_size, const char *temp_array_buffer, long &buffer_content_iter, bool report_error)
{
    if (data_size > buffer_content_iter)
        if (report_error)
            EXECUTION_REPORT(REPORT_ERROR,-1, false, "Software error in read_data_from_array_buffer");
        else return false;
    
    for (int i = 0; i < data_size; i ++)
        ((char*) data)[i] = temp_array_buffer[buffer_content_iter-data_size+i];
    
    buffer_content_iter -= data_size;

    return true;
}


bool get_next_line(char *line, FILE *fp)
{
    char c;
    int iter = 0;
    

    while (!feof(fp) && (c = getc(fp)) != -1) {
        if (c == '\n') 
            break;
        line[iter ++] = c;
    }
    line[iter ++] = '\0';
    if (iter == 1)
        return false;

    return true;
}


bool get_next_attr(char *attr, char **line)
{
    EXECUTION_REPORT(REPORT_ERROR, -1, *line != NULL, "Can not get next attribute from the configuration file. There may be problem in the configuration file");
    
    if ((*line)[0] == '\0') {
        (*line) = NULL;
        return false;
    }

    while ((*line)[0] == ' ' || (*line)[0] == '\t')
        (*line) ++;
    
    while ((*line)[0] != '\0' && (*line)[0] != '\t') {
        *attr = (*line)[0];
        attr ++;
        (*line) ++;
    }
    if ((*line)[0] == '\t')
        (*line) ++;

    if (*attr == ' ' || *attr == '\t') {
        while (*attr == ' ' || *attr == '\t')
            attr --;
        attr ++;
    }
    *attr = '\0';
    
    return true;
}


bool is_end_of_file(FILE *fp)
{
    long offset;

    offset = ftell(fp);
    getc(fp);
    if (feof(fp))
        return true;

    fseek(fp, offset, SEEK_SET);
    return false;
}


bool get_next_integer_attr(char **line, int &value)
{
    char attr[NAME_STR_SIZE];

    
    value = -1;

    if (!get_next_attr(attr, line))
        return false;
    for (int i = 0; i < strlen(attr); i ++) {
        if (i == 0 && attr[0] == '-')
            continue;
        if (attr[i]-'0'< 0 || attr[i]-'9' > 0) 
            return false;        
    }
    value = atoi(attr);

    return true;
}


bool get_next_double_attr(char **line, double &value)
{
    char attr[NAME_STR_SIZE];
    int dot_number = 0;
    int figure_number = 0;


    if (!get_next_attr(attr, line))
        return false;

    for (int i = 0; i < strlen(attr); i ++) {
        if (i == 0 && attr[0] == '-')
            continue;
        if (attr[i] == '.') {
            dot_number ++;
            if (figure_number == 0)
                return false;
            if (dot_number > 1)
                return false;
            continue;
        }
        if (attr[i]-'0' < 0 || attr[i]-'9' > 0) 
            return false;
        figure_number ++;
    }
    value = atof(attr);

    return true;
}


void check_for_ccpl_managers_allocated(int API_ID, const char *annotation)
{
    char API_label[NAME_STR_SIZE];
    

    get_API_hint(-1, API_ID, API_label);
    EXECUTION_REPORT(REPORT_ERROR, -1, comp_comm_group_mgt_mgr != NULL, "Error happens when calling the API \"%s\": the stage of registering coupling configurations has not been started or the C-Coupler has been finalized. Please call the API \"CCPL_register_component\" for the registration of the root component model (parameter \"parent_id\" should be -1) to start the configuration stage. Please check the model code related to the annotation \"%s\".", API_label, annotation);
}


void check_for_coupling_registration_stage(int comp_id, int API_ID, bool require_real_model, const char *annotation)
{
    char API_label[NAME_STR_SIZE];
    

    get_API_hint(-1, API_ID, API_label);
    check_for_ccpl_managers_allocated(API_ID, annotation);
    EXECUTION_REPORT(REPORT_ERROR, -1, comp_comm_group_mgt_mgr->is_legal_local_comp_id(comp_id,true), "Error happens when calling the API \"%s\": The ID of the given component model (currently is 0x%x) is wrong (not the legal ID of a component). Please check the model code with the annotation \"%s\"", API_label, comp_id, annotation);
    comp_comm_group_mgt_mgr->confirm_coupling_configuration_active(comp_id, API_ID, require_real_model, annotation);        
}


void common_checking_for_grid_registration(int comp_id, const char *grid_name, const char *coord_unit, int API_id, const char *annotation)
{
    char API_label[NAME_STR_SIZE];
    Original_grid_info *existing_grid;

    
    get_API_hint(comp_id, API_id, API_label);
    check_for_coupling_registration_stage(comp_id, API_id, true, annotation);
    check_API_parameter_string_length(comp_id, API_id, CCPL_NAME_STR_LEN, grid_name, "grid_name", annotation);
    if (coord_unit != NULL)
        check_API_parameter_string_length(comp_id, API_id, CCPL_NAME_STR_LEN, coord_unit, "coord_unit", annotation);
    existing_grid = original_grid_mgr->search_grid_info(grid_name, comp_id);
    if (existing_grid != NULL)
        EXECUTION_REPORT(REPORT_ERROR, comp_id, false, "Error happens when calling the API \"%s\" to register a grid \"%s\" at the model code with the annotation \"%s\": another grid with the same name has already been registered before (at the model code with the annotation \"%s\"). Please verify.", API_label, grid_name, annotation, annotation_mgr->get_annotation(existing_grid->get_grid_id(), "grid_registration"));
    check_and_verify_name_format_of_string_for_API(comp_id, grid_name, API_id, "the C-Coupler grid", annotation);
    MPI_Comm comm = comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id, "in grid registration");
    synchronize_comp_processes_for_API(comp_id, API_id, comm, "registering a grid", annotation);
    check_API_parameter_string(comp_id, API_id, comm, "registering a grid", grid_name, "grid_name", annotation);    
    if (coord_unit != NULL)
        check_API_parameter_string(comp_id, API_id, comm, "registering a grid", coord_unit, "coord_unit", annotation);
}


bool are_two_coord_arrays_same(double *array1, double *array2, int array_size1, int array_size2)
{
    double eps = 0.00001;

    
    if (array_size1 != array_size2)
        return false;

    for (int i = 0; i < array_size1; i ++)
        if (fabs(array1[i]-array2[i]) >= eps)
            return false;

    return true;
}


void transform_datatype_of_arrays(const char *src_array, char *dst_array, const char *src_data_type, const char *dst_data_type, long num_local_cells)
{
    if (words_are_the_same(src_data_type,DATA_TYPE_FLOAT) && words_are_the_same(dst_data_type,DATA_TYPE_FLOAT))
        transform_datatype_of_arrays((const float*)src_array, (float*) dst_array, num_local_cells);
    else if (words_are_the_same(src_data_type,DATA_TYPE_DOUBLE) && words_are_the_same(dst_data_type,DATA_TYPE_DOUBLE))
        transform_datatype_of_arrays((const double*)src_array, (double*) dst_array, num_local_cells);
    else if (words_are_the_same(src_data_type,DATA_TYPE_FLOAT) && words_are_the_same(dst_data_type,DATA_TYPE_DOUBLE))
        transform_datatype_of_arrays((const float*)src_array, (double*) dst_array, num_local_cells);
    else if (words_are_the_same(src_data_type,DATA_TYPE_DOUBLE) && words_are_the_same(dst_data_type,DATA_TYPE_FLOAT))
        transform_datatype_of_arrays((const double*)src_array, (float*) dst_array, num_local_cells);
    else EXECUTION_REPORT(REPORT_ERROR, -1, false, "Software error in transform_datatype_of_arrays: data type transformation from %s to %s is not supported", src_data_type, dst_data_type);
}


void check_API_parameter_string_length(int comp_id, int API_ID, int str_max_size, const char *str, const char *parameter_name, const char *annotation)
{
    char API_label[NAME_STR_SIZE];


    get_API_hint(-1, API_ID, API_label);
    
    EXECUTION_REPORT(REPORT_ERROR, comp_id, strlen(str) <= str_max_size, "Error happens when calling the API \"%s\": the string size (currently is %d) of the parameter \"%s\" (the string is \"%s\") is larger than the limit (%d). Please verify the model code with the annotation \"%s\".", API_label, strlen(str), parameter_name, str, str_max_size, annotation);
    EXECUTION_REPORT(REPORT_ERROR, comp_id, !words_are_the_same(str, "NULL"), "Error happens when calling the API \"%s\": the parameter \"%s\" cannot be \"NULL\". Please verify the model code with the annotation \"%s\".", API_label, parameter_name, str, annotation);
}


void check_XML_attribute_value_string_length(int comp_id, int str_max_size, const char *XML_attribute, const char *XML_value, const char *XML_file_name, int line_number)
{
    EXECUTION_REPORT(REPORT_ERROR, comp_id, strlen(XML_value) <= str_max_size, "Error happens when using the XML configuration file \"%s\": the string size (currently is %d) of the value (\"%s\") the XML attribute \"%s\" is larger than the limit (%d). Please verify XML file arround the line %d.", XML_file_name, strlen(XML_value), XML_value, XML_attribute, str_max_size, line_number);
}


bool is_string_decimal_number(const char *string)
{
    int length = strlen(string);

    if (length == 0)
        return false;
    
    for (int i = 0; i < length; i ++) {
        if (i == 0 && string[i] == '-' && length > 1)
            continue;
        if (!(string[i] >= '0' && string[i] <= '9'))
            return false;
    }
    
    return true;
}


char *load_string(char *str, long &str_size, long max_size, const char *array_buffer, long &buffer_content_iter, const char *file_name)
{
    char *local_str = str;
    

    if (!read_data_from_array_buffer(&str_size, sizeof(long), array_buffer, buffer_content_iter, file_name == NULL))
        if (file_name != NULL)
            EXECUTION_REPORT(REPORT_ERROR, -1, false, "Fail to load the restart data file \"%s\": its format is wrong, or the information it includes is not complete. Please try a different restart time with complete restart data files.", file_name);    
        else EXECUTION_REPORT(REPORT_ERROR, -1, false, "Software error in load_string");
    EXECUTION_REPORT(REPORT_ERROR, -1, str_size > 0, "Fail to load the restart data file \"%s\": its format is wrong, or the information it includes is not complete. Please try a different restart time with complete restart data files.", file_name);
    if (local_str != NULL)
        EXECUTION_REPORT(REPORT_ERROR, -1, str_size < max_size, "Fail to load the restart data file \"%s\": its format is wrong, or the information it includes is not complete. Please try a different restart time with complete restart data files.", file_name);
    if (local_str == NULL)
        local_str = new char [str_size+1];
    EXECUTION_REPORT(REPORT_ERROR, -1, read_data_from_array_buffer(local_str, str_size, array_buffer, buffer_content_iter, file_name == NULL), "Fail to load the restart data file \"%s\": its format is wrong, or the information it includes is not complete. Please try a different restart time with complete restart data files.", file_name);
    local_str[str_size] = '\0';

    return local_str;
}


void dump_string(const char *str, long str_size, char **array_buffer, long &buffer_max_size, long &buffer_content_size)
{
    if (str_size == -1)
        if (str != NULL)
            str_size = strlen(str);
        else str_size = 0;
    if (str != NULL)
        write_data_into_array_buffer(str, str_size, array_buffer, buffer_max_size, buffer_content_size);
    write_data_into_array_buffer(&str_size, sizeof(long), array_buffer, buffer_max_size, buffer_content_size);
}


long get_restart_time_in_rpointer_file(const char *file_name)
{
    char line[NAME_STR_SIZE*16], date_str[NAME_STR_SIZE], second_str[NAME_STR_SIZE];
    int date, second, special_pos;
    FILE *rpointer_file;


    rpointer_file = fopen(file_name, "r");
    get_next_line(line, rpointer_file);
    fclose(rpointer_file);
    for (special_pos = strlen(line)-1; special_pos >= 0; special_pos --)
        if (line[special_pos] == '-')
            break;
    EXECUTION_REPORT(REPORT_ERROR, -1, special_pos > 10, "Error happens in a continue run: the restart file name \"%s\" in the rpointer file \"%s\" is not in a right format", file_name);
    EXECUTION_REPORT(REPORT_ERROR, -1, line[special_pos-9] == '.' && line[special_pos-10] == 'r', "Error happens in a continue run: the restart file name \"%s\" in the rpointer file \"%s\" is not in a right format", file_name);
    strncpy(date_str, line+special_pos-8, 8);
    strncpy(second_str, line+special_pos+1, 5);
	date_str[8] = '\0';
	second_str[5] = '\0';
    EXECUTION_REPORT(REPORT_ERROR, -1, sscanf(date_str, "%d", &date) == 1, "Error happens in a continue run: the restart file name \"%s\" in the rpointer file \"%s\" is not in a right format", file_name);
    EXECUTION_REPORT(REPORT_ERROR, -1, sscanf(second_str, "%d", &second) == 1, "Error happens in a continue run: the restart file name \"%s\" in the rpointer file \"%s\" is not in a right format", file_name);

    return ((long)date)*((long)100000) + second;
}

