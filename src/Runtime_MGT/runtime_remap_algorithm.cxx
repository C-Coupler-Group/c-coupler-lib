/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#include <mpi.h>
#include "global_data.h"
#include "runtime_remap_algorithm.h"
#include "cor_cpl_interface.h"
#include "cor_global_data.h"
#include <stdio.h>


Runtime_remap_algorithm::Runtime_remap_algorithm(Runtime_remapping_weights *runtime_remapping_weights, Field_mem_info *src_field_instance, Field_mem_info *dst_field_instance, int connection_id)
{    
    comp_id = dst_field_instance->get_comp_id();
    specified_src_field_instance = src_field_instance;
    specified_dst_field_instance = dst_field_instance;
    this->runtime_remapping_weights = runtime_remapping_weights;
    
    if (words_are_the_same(src_field_instance->get_field_data()->get_grid_data_field()->data_type_in_application, DATA_TYPE_FLOAT)) {
        true_src_field_instance = memory_manager->alloc_mem(specified_src_field_instance, BUF_MARK_REMAP_DATATYPE_TRANS_SRC, connection_id, DATA_TYPE_DOUBLE, false);
        true_dst_field_instance = memory_manager->alloc_mem(specified_dst_field_instance, BUF_MARK_REMAP_DATATYPE_TRANS_DST, connection_id, DATA_TYPE_DOUBLE, false);
        transform_data_type = true;
    }
    else if (words_are_the_same(src_field_instance->get_field_data()->get_grid_data_field()->data_type_in_application, DATA_TYPE_DOUBLE)) {
        true_src_field_instance = specified_src_field_instance;
        true_dst_field_instance = specified_dst_field_instance;
        transform_data_type = false;
    }
    else EXECUTION_REPORT(REPORT_ERROR, -1, false, "Software error in Runtime_remap_algorithm::Runtime_remap_algorithm: data type is wrong: data type %s of src field %s vs data type %s of dst field %s", src_field_instance->get_field_data()->get_grid_data_field()->data_type_in_application, src_field_instance->get_field_name(), dst_field_instance->get_field_data()->get_grid_data_field()->data_type_in_application, dst_field_instance->get_field_name());
}


Runtime_remap_algorithm::~Runtime_remap_algorithm()
{
}


void Runtime_remap_algorithm::do_remap(bool is_algorithm_in_kernel_stage)
{
    long i, j, field_size_src_before_rearrange, field_size_src_after_rearrange;
    long field_size_dst, decomp_size_src_before_rearrange, decomp_size_src_after_rearrange, num_levels;
    double *src_frac_values, *dst_frac_values, *src_field_values, *dst_field_values, frac_1x;
    double *temp_src_values;
    Remap_grid_class *decomp_grid, *original_grid;
    Decomp_info *decomp;


	EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Data interpolation for field \"%s\"", specified_src_field_instance->get_field_name());

    if (runtime_remapping_weights->get_parallel_remapping_weights() == NULL)
        return;
    
    specified_src_field_instance->use_field_values("");
//    specified_src_field_instance->check_field_sum("before data interpolation");
    if (transform_data_type)
        for (int i = 0; i < specified_src_field_instance->get_size_of_field(); i ++)
            ((double*)true_src_field_instance->get_data_buf())[i] = ((float*)specified_src_field_instance->get_data_buf())[i];
	if (!words_are_the_same(specified_src_field_instance->get_field_name(),V3D_GRID_3D_LEVEL_FIELD_NAME))
	    runtime_remapping_weights->renew_dynamic_V1D_remapping_weights();
    comp_comm_group_mgt_mgr->get_global_node_of_local_comp(runtime_remapping_weights->get_dst_original_grid()->get_comp_id(),false,"")->get_performance_timing_mgr()->performance_timing_start(TIMING_TYPE_COMPUTATION, -1, -1, "remapping cal");
    runtime_remapping_weights->get_parallel_remapping_weights()->do_remap(runtime_remapping_weights->get_dst_original_grid()->get_comp_id(), true_src_field_instance->get_field_data(), true_dst_field_instance->get_field_data());
    comp_comm_group_mgt_mgr->get_global_node_of_local_comp(runtime_remapping_weights->get_dst_original_grid()->get_comp_id(),false,"")->get_performance_timing_mgr()->performance_timing_stop(TIMING_TYPE_COMPUTATION, -1, -1, "remapping cal");
    if (transform_data_type)
        for (int i = 0; i < specified_dst_field_instance->get_size_of_field(); i ++)
            ((float*)specified_dst_field_instance->get_data_buf())[i] = ((double*)true_dst_field_instance->get_data_buf())[i];
    specified_dst_field_instance->define_field_values(true);
//    specified_dst_field_instance->check_field_sum("after data interpolation");
}


bool Runtime_remap_algorithm::run(bool is_algorithm_in_kernel_stage)
{
    do_remap(is_algorithm_in_kernel_stage);

    return true;
}

