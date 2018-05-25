/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#include "cor_global_data.h"
#include <string.h>

IO_mgt *io_manager;
Remap_strategy_mgt *remap_strategy_manager;
Remap_grid_mgt *remap_grid_manager;
Remap_operator_mgt *remap_operator_manager;
Remap_field_data_mgt *remap_field_data_manager;
Remap_weight_of_strategy_mgt *remap_weights_of_strategy_manager;
Remap_weight_of_operator_mgt *sequential_remap_weight_of_operator_manager;
Remap_weight_of_operator_mgt *parallel_remap_weight_of_operator_manager;
Runtime_remap_function *current_runtime_remap_function;
Remap_operator_grid *current_runtime_remap_operator_grid_src;
Remap_operator_grid *current_runtime_remap_operator_grid_dst;
Remap_operator_basis *current_runtime_remap_operator;
bool is_coord_unit_degree[256];
bool is_master_process_in_computing_node = true;



int line_number = -1;
int execution_phase_number;


int get_data_type_size(const char *data_type)
{
    if (words_are_the_same(data_type, DATA_TYPE_DOUBLE))
        return sizeof(double);
    else if (words_are_the_same(data_type, DATA_TYPE_FLOAT))
        return sizeof(float);
    else if (words_are_the_same(data_type, DATA_TYPE_LONG))
        return sizeof(long);
    else if (words_are_the_same(data_type, DATA_TYPE_INT))
        return sizeof(int);
    else if (words_are_the_same(data_type, DATA_TYPE_BOOL))
        return sizeof(bool);
    else if (words_are_the_same(data_type, DATA_TYPE_CHAR))
        return sizeof(char);
    else if (words_are_the_same(data_type, DATA_TYPE_SHORT))
        return sizeof(short);
    else if (words_are_the_same(data_type, DATA_TYPE_STRING))
        return NAME_STR_SIZE;
    else EXECUTION_REPORT(REPORT_ERROR, -1, false, 
                      "implicit data type %s is disabled in this software\n", 
                      data_type);

    return 0;
}


void check_application_io_datatype_consistency(const char *field_name, const char *datatype_application, const char *datatype_io)
{
    if (words_are_the_same(datatype_application, DATA_TYPE_DOUBLE) || words_are_the_same(datatype_application, DATA_TYPE_FLOAT)) 
        EXECUTION_REPORT(REPORT_ERROR, -1, words_are_the_same(datatype_io, DATA_TYPE_DOUBLE) || words_are_the_same(datatype_io, DATA_TYPE_FLOAT) || words_are_the_same(datatype_io, DATA_TYPE_SHORT) ,
                     "the data type of field %s in IO file must be one of short, real4 and real8, as its data type in application is %s\n", field_name, datatype_application);
    else if (words_are_the_same(datatype_application, DATA_TYPE_LONG) || words_are_the_same(datatype_application, DATA_TYPE_INT) || words_are_the_same(datatype_application, DATA_TYPE_BOOL)) 
        EXECUTION_REPORT(REPORT_ERROR, -1, words_are_the_same(datatype_io, DATA_TYPE_INT),
                     "the data type of field %s in IO file must be integer, as its data type in application is %s\n", field_name, datatype_application);
    else if (words_are_the_same(datatype_application, DATA_TYPE_SHORT))
        EXECUTION_REPORT(REPORT_ERROR, -1, words_are_the_same(datatype_io, DATA_TYPE_SHORT),
                     "the data type of field %s in IO file must be short, as its data type in application is %s\n", field_name, datatype_application);
    
    if (words_are_the_same(datatype_io, DATA_TYPE_INT))
        EXECUTION_REPORT(REPORT_ERROR, -1, words_are_the_same(datatype_application, DATA_TYPE_INT) || words_are_the_same(datatype_application, DATA_TYPE_BOOL) || words_are_the_same(datatype_application, DATA_TYPE_LONG),
                     "the data type of field %s in application must be integer, long or logical, as its data type in IO file is %s\n", field_name, datatype_io);        
    else if (words_are_the_same(datatype_io, DATA_TYPE_SHORT))
        EXECUTION_REPORT(REPORT_ERROR, -1, words_are_the_same(datatype_application, DATA_TYPE_SHORT) || words_are_the_same(datatype_application, DATA_TYPE_FLOAT) || words_are_the_same(datatype_application, DATA_TYPE_DOUBLE),
                     "the data type of field %s in application must be short, real4 or real8, as its data type in IO file is %s\n", field_name, datatype_io);   
    else if (words_are_the_same(datatype_io, DATA_TYPE_FLOAT) || words_are_the_same(datatype_io, DATA_TYPE_DOUBLE))
        EXECUTION_REPORT(REPORT_ERROR, -1, words_are_the_same(datatype_application, DATA_TYPE_FLOAT) || words_are_the_same(datatype_application, DATA_TYPE_DOUBLE),
                     "the data type of field %s in application must be real4 or real8, as its data type in IO file is %s\n", field_name, datatype_io);   
}

