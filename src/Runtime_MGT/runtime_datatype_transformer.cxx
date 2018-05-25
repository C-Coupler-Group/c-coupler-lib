/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/

#include "global_data.h"
#include "runtime_datatype_transformer.h"
#include "cor_global_data.h"


Runtime_datatype_transformer::Runtime_datatype_transformer(Field_mem_info *src_field, Field_mem_info *dst_field)
{
    src_fields.push_back(src_field);
    dst_fields.push_back(dst_field);
}


bool Runtime_datatype_transformer::run(bool bypass_timer)
{
    transform_fields_datatype();
    
    return true;
}


void Runtime_datatype_transformer::transform_fields_datatype()
{
    char *data_type_src, *data_type_dst;
    long num_local_cells;


    for (int i = 0; i < src_fields.size(); i ++) {
        src_fields[i]->use_field_values("");
        dst_fields[i]->define_field_values(false);
        data_type_src = src_fields[i]->get_field_data()->get_grid_data_field()->data_type_in_application;
        data_type_dst = dst_fields[i]->get_field_data()->get_grid_data_field()->data_type_in_application;
        num_local_cells = src_fields[i]->get_field_data()->get_grid_data_field()->required_data_size;
        if (words_are_the_same(data_type_src, DATA_TYPE_DOUBLE)) {
            if (words_are_the_same(data_type_dst, DATA_TYPE_FLOAT))
                transform_datatype_of_arrays((double*)src_fields[i]->get_data_buf(), (float*) dst_fields[i]->get_data_buf(), num_local_cells);
            else EXECUTION_REPORT(REPORT_ERROR,-1, "C-Coupler software error1 in transform_fields_datatype of Runtime_datatype_transformer");
        }
        else if (words_are_the_same(data_type_src, DATA_TYPE_FLOAT)) {
            if (words_are_the_same(data_type_dst, DATA_TYPE_DOUBLE))
                transform_datatype_of_arrays((float*)src_fields[i]->get_data_buf(), (double*) dst_fields[i]->get_data_buf(), num_local_cells);
            else EXECUTION_REPORT(REPORT_ERROR,-1, "C-Coupler software error2 in transform_fields_datatype of Runtime_datatype_transformer");
        }
        else if (words_are_the_same(data_type_src, DATA_TYPE_LONG)) {
            if (words_are_the_same(data_type_dst, DATA_TYPE_INT))
                transform_datatype_of_arrays((long*)src_fields[i]->get_data_buf(), (int*) dst_fields[i]->get_data_buf(), num_local_cells);
            else if (words_are_the_same(data_type_dst, DATA_TYPE_SHORT))
                transform_datatype_of_arrays((long*)src_fields[i]->get_data_buf(), (short*) dst_fields[i]->get_data_buf(), num_local_cells);
            else if (words_are_the_same(data_type_dst, DATA_TYPE_BOOL))
                transform_datatype_of_arrays((long*)src_fields[i]->get_data_buf(), (bool*) dst_fields[i]->get_data_buf(), num_local_cells);
            else EXECUTION_REPORT(REPORT_ERROR,-1, "C-Coupler software error3 in transform_fields_datatype of Runtime_datatype_transformer");
        }
        else if (words_are_the_same(data_type_src, DATA_TYPE_INT)) {
            if (words_are_the_same(data_type_dst, DATA_TYPE_LONG))
                transform_datatype_of_arrays((int*)src_fields[i]->get_data_buf(), (long*) dst_fields[i]->get_data_buf(), num_local_cells);
            else if (words_are_the_same(data_type_dst, DATA_TYPE_SHORT))
                transform_datatype_of_arrays((int*)src_fields[i]->get_data_buf(), (short*) dst_fields[i]->get_data_buf(), num_local_cells);
            else if (words_are_the_same(data_type_dst, DATA_TYPE_BOOL))
                transform_datatype_of_arrays((int*)src_fields[i]->get_data_buf(), (bool*) dst_fields[i]->get_data_buf(), num_local_cells);
            else EXECUTION_REPORT(REPORT_ERROR,-1, "C-Coupler software error4 in transform_fields_datatype of Runtime_datatype_transformer");
        }
        else if (words_are_the_same(data_type_src, DATA_TYPE_SHORT)) {
            if (words_are_the_same(data_type_dst, DATA_TYPE_LONG))
                transform_datatype_of_arrays((short*)src_fields[i]->get_data_buf(), (long*) dst_fields[i]->get_data_buf(), num_local_cells);
            else if (words_are_the_same(data_type_dst, DATA_TYPE_INT))
                transform_datatype_of_arrays((short*)src_fields[i]->get_data_buf(), (int*) dst_fields[i]->get_data_buf(), num_local_cells);
            else if (words_are_the_same(data_type_dst, DATA_TYPE_BOOL))
                transform_datatype_of_arrays((short*)src_fields[i]->get_data_buf(), (bool*) dst_fields[i]->get_data_buf(), num_local_cells);
            else EXECUTION_REPORT(REPORT_ERROR,-1, "C-Coupler software error5 in transform_fields_datatype of Runtime_datatype_transformer");
        }
        else if (words_are_the_same(data_type_src, DATA_TYPE_BOOL)) {
            if (words_are_the_same(data_type_dst, DATA_TYPE_LONG))
                transform_datatype_of_arrays((bool*)src_fields[i]->get_data_buf(), (long*) dst_fields[i]->get_data_buf(), num_local_cells);
            else if (words_are_the_same(data_type_dst, DATA_TYPE_INT))
                transform_datatype_of_arrays((bool*)src_fields[i]->get_data_buf(), (int*) dst_fields[i]->get_data_buf(), num_local_cells);
            else if (words_are_the_same(data_type_dst, DATA_TYPE_SHORT))
                transform_datatype_of_arrays((bool*)src_fields[i]->get_data_buf(), (short*) dst_fields[i]->get_data_buf(), num_local_cells);
            else EXECUTION_REPORT(REPORT_ERROR,-1, "C-Coupler software error6 in transform_fields_datatype of Runtime_datatype_transformer");
        }
        else EXECUTION_REPORT(REPORT_ERROR,-1, "C-Coupler software error7 in transform_fields_datatype of Runtime_datatype_transformer");
    }
}

