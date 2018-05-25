/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#include "global_data.h"
#include "runtime_cumulate_average_algorithm.h"
#include "Runtime_Algorithm_Basis.h"
#include <string.h>


int global_algorithm_id;


Runtime_algorithm_basis::Runtime_algorithm_basis()
{
    num_src_fields = 0;
    num_dst_fields = 0;
    cumulate_average_algorithm_before_run = NULL;
    algorithm_id = global_algorithm_id;
    global_algorithm_id ++;

    comp_names = NULL;
    field_names = NULL;
    field_local_decomp_names = NULL;
    field_grid_names = NULL;
    buf_marks = NULL;
    average_mark = NULL;
}


Runtime_algorithm_basis::~Runtime_algorithm_basis()
{
//    if (num_src_fields + num_dst_fields > 0)
//        EXECUTION_REPORT(REPORT_ERROR,-1, comp_names == NULL && field_names == NULL && field_local_decomp_names == NULL && field_grid_names == NULL && buf_marks == NULL && average_mark == NULL, "C-Coupler software error when deleting Runtime_algorithm_basis");

    if (comp_names == NULL)
        return;

    for (int i = 0; i < num_src_fields+num_dst_fields; i ++) {
        delete [] comp_names[i];
        delete [] field_names[i];
        delete [] field_local_decomp_names[i];
        delete [] field_grid_names[i];
    }
    delete [] comp_names;
    delete [] field_names;
    delete [] field_local_decomp_names;
    delete [] field_grid_names;
    delete [] buf_marks;
    delete [] average_mark;
}


void Runtime_algorithm_basis::runtime_algorithm_common_initialize(const int num_src_fields, const int num_dst_fields)
{
    this->num_src_fields = num_src_fields;
    this->num_dst_fields = num_dst_fields;
}


void Runtime_algorithm_basis::cumulate_average_before_run(bool is_algorithm_in_kernel_stage)
{
    if (cumulate_average_algorithm_before_run != NULL) {
        EXECUTION_REPORT_LOG(REPORT_LOG,-1, true, "before implicit cumulating and averaging");
        cumulate_average_algorithm_before_run->run(is_algorithm_in_kernel_stage);
        EXECUTION_REPORT_LOG(REPORT_LOG,-1, true, "after implicit cumulating and averaging");
    }
}


void Runtime_algorithm_basis::transfer_fields_data_type_before_run() 
{ 
}


void Runtime_algorithm_basis::transfer_fields_data_type_after_run() 
{
}


void Runtime_algorithm_basis::allocate_basic_data_structure(int num_src_fields, int num_dst_fields)
{    
    this->num_src_fields = num_src_fields;
    this->num_dst_fields = num_dst_fields;

    EXECUTION_REPORT(REPORT_ERROR,-1, num_src_fields >= 0 && num_dst_fields >= 0,
                     "C-Coupler software error in allocate_basic_data_structure for Runtime_algorithm_basis");

    if (num_src_fields + num_dst_fields == 0)
        return;

    comp_names = new char* [num_src_fields+num_dst_fields];
    field_names = new char* [num_src_fields+num_dst_fields];
    field_local_decomp_names = new char* [num_src_fields+num_dst_fields];
    field_grid_names = new char* [num_src_fields+num_dst_fields];
    buf_marks = new int [num_src_fields+num_dst_fields];
    average_mark = new bool [num_src_fields+num_dst_fields];    

    for (int i = 0; i < num_src_fields+num_dst_fields; i ++) {
        comp_names[i] = new char [NAME_STR_SIZE];
        field_names[i] = new char [NAME_STR_SIZE];
        field_local_decomp_names[i] = new char [NAME_STR_SIZE];
        field_grid_names[i] = new char [NAME_STR_SIZE];
        buf_marks[i] = -1;
        average_mark[i] = false;
    }
}

