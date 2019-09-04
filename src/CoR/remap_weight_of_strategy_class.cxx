/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#include "cor_global_data.h"
#include "remap_weight_of_strategy_class.h"
#include "remap_strategy_class.h"
#include "remap_operator_basis.h"
#include "io_binary.h"
#include "io_netcdf.h"
#include "performance_timing_mgt.h"
#include "global_data.h"
#include <string.h>


Remap_weight_of_operator_instance_class::Remap_weight_of_operator_instance_class(Remap_grid_class *field_data_grid_src, Remap_grid_class *field_data_grid_dst, 
                                                               long remap_beg_iter, Remap_operator_basis *remap_operator)
{
	this->remap_weight_of_operator = NULL;
    this->remap_beg_iter = remap_beg_iter;
    this->remap_end_iter = remap_beg_iter + 1;
	this->duplicated_remap_operator = NULL;
    if (remap_operator->get_src_grid()->get_is_sphere_grid()) {
		this->duplicated_remap_operator = remap_operator->duplicate_remap_operator(true);
    }
}


Remap_weight_of_operator_instance_class::Remap_weight_of_operator_instance_class(Remap_grid_class *field_data_grid_src, Remap_grid_class *field_data_grid_dst, 
                                                               long remap_beg_iter, Remap_operator_basis *remap_operator, Remap_operator_basis *duplicated_remap_operator)
{
	this->remap_weight_of_operator = NULL;
    this->remap_beg_iter = remap_beg_iter;
    this->remap_end_iter = remap_beg_iter + 1;
    this->duplicated_remap_operator = duplicated_remap_operator;
}


Remap_weight_of_operator_instance_class *Remap_weight_of_operator_instance_class::generate_parallel_remap_weights(Remap_grid_class **decomp_original_grids, int **global_cells_local_indexes_in_decomps)
{
    Remap_weight_of_operator_instance_class *parallel_remap_weights_of_operator_instance;
    int overlap_with_decomp_counter = 0;


    parallel_remap_weights_of_operator_instance = new Remap_weight_of_operator_instance_class();
	parallel_remap_weights_of_operator_instance->remap_weight_of_operator = this->remap_weight_of_operator;

    for (int i = 0; i < 2; i ++)
        if (this->get_original_remap_operator()->get_src_grid()->have_overlap_with_grid(decomp_original_grids[i])) {
            EXECUTION_REPORT(REPORT_ERROR, -1, decomp_original_grids[i]->is_subset_of_grid(this->get_original_remap_operator()->get_src_grid()),
                         "C-Coupler error1 in generate_parallel_remap_weights of Remap_weight_of_operator_instance_class\n");
            overlap_with_decomp_counter ++;
        }
    for (int i = 0; i < 2; i ++)
        if (this->get_original_remap_operator()->get_dst_grid()->have_overlap_with_grid(decomp_original_grids[i])) {
            EXECUTION_REPORT(REPORT_ERROR, -1, decomp_original_grids[i]->is_subset_of_grid(this->get_original_remap_operator()->get_dst_grid()),
                         "C-Coupler error2 in generate_parallel_remap_weights of Remap_weight_of_operator_instance_class\n");
            overlap_with_decomp_counter ++;
        }

    if (overlap_with_decomp_counter > 0) {
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, this->duplicated_remap_operator != NULL, "C-Coupler error4 in generate_parallel_remap_weights of Remap_weight_of_operator_instance_class\n");
        parallel_remap_weights_of_operator_instance->duplicated_remap_operator = this->duplicated_remap_operator->generate_parallel_remap_operator(decomp_original_grids, global_cells_local_indexes_in_decomps);
    }
    else {
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, this->get_original_remap_operator() != NULL && parallel_remap_weights_of_operator_instance->duplicated_remap_operator == NULL, "C-Coupler error4 in generate_parallel_remap_weights of Remap_weight_of_operator_instance_class\n");
		parallel_remap_weights_of_operator_instance->duplicated_remap_operator = this->get_original_remap_operator()->duplicate_remap_operator(true);
    }

    return parallel_remap_weights_of_operator_instance;
}


Remap_weight_of_operator_instance_class::~Remap_weight_of_operator_instance_class()
{
    if (duplicated_remap_operator != NULL)
        delete duplicated_remap_operator;
}


void Remap_weight_of_operator_instance_class::renew_remapping_time_end_iter(long time_end_iter)
{
    EXECUTION_REPORT(REPORT_ERROR, -1, remap_end_iter == time_end_iter, "C-Coupler error in Remap_weight_of_operator_instance_class::renew_remapping_time_end_iter");
    remap_end_iter = time_end_iter + 1;
}


Remap_grid_class *Remap_weight_of_operator_instance_class::get_field_data_grid_src() 
{
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, remap_weight_of_operator != NULL, "Software error in Remap_weight_of_operator_instance_class::get_field_data_grid_src");
	return remap_weight_of_operator->get_field_data_grid_src(); 
}


Remap_grid_class *Remap_weight_of_operator_instance_class::get_field_data_grid_dst() 
{
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, remap_weight_of_operator != NULL, "Software error in Remap_weight_of_operator_instance_class::get_field_data_grid_dst");
	return remap_weight_of_operator->get_field_data_grid_dst(); 
}


Remap_operator_basis *Remap_weight_of_operator_instance_class::get_original_remap_operator()
{
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, remap_weight_of_operator != NULL, "Software error in Remap_weight_of_operator_instance_class::get_original_remap_operator");
	return remap_weight_of_operator->get_original_remap_operator(); 
}


Remap_grid_class *Remap_weight_of_operator_instance_class::get_operator_grid_src()
{
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, remap_weight_of_operator != NULL, "Software error in Remap_weight_of_operator_instance_class::get_operator_grid_src");
	return remap_weight_of_operator->get_operator_grid_src();
}


Remap_grid_class *Remap_weight_of_operator_instance_class::get_operator_grid_dst()
{
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, remap_weight_of_operator != NULL, "Software error in Remap_weight_of_operator_instance_class::get_operator_grid_dst");
	return remap_weight_of_operator->get_operator_grid_dst();
}


Remap_weight_of_operator_class::Remap_weight_of_operator_class(Remap_grid_class *field_data_grid_src, Remap_grid_class *field_data_grid_dst, Remap_operator_basis *remap_operator,
                                                                      Remap_grid_class *operator_grid_src, Remap_grid_class *operator_grid_dst)
{
    this->field_data_grid_src = field_data_grid_src;
    this->field_data_grid_dst = field_data_grid_dst;
    this->original_remap_operator = remap_operator;
    this->operator_grid_src = operator_grid_src;
    this->operator_grid_dst = operator_grid_dst;
    empty_remap_weight = false;
}


void Remap_weight_of_operator_class::calculate_src_decomp(long *decomp_map_src, const long *decomp_map_dst)
{
    long remap_beg_iter, remap_end_iter;
    

    for (int i = 0; i < remap_weights_of_operator_instances.size(); i ++) {
        remap_beg_iter = remap_weights_of_operator_instances[i]->remap_beg_iter;
        if (remap_weights_of_operator_instances[i]->remap_end_iter != -1)
            remap_end_iter = remap_weights_of_operator_instances[i]->remap_end_iter;
        else if (i+1 < remap_weights_of_operator_instances.size())
            remap_end_iter = remap_weights_of_operator_instances[i+1]->remap_beg_iter;
        else remap_end_iter = field_data_grid_src->get_grid_size()/operator_grid_src->get_grid_size();
        for (int j = remap_beg_iter; j < remap_end_iter; j ++) {
            EXECUTION_REPORT(REPORT_ERROR, -1, remap_weights_of_operator_instances[i]->duplicated_remap_operator != NULL, "C-Coupler error3 in do_remap of Remap_weight_of_operator_class");
            remap_weights_of_operator_instances[i]->duplicated_remap_operator->do_src_decomp_caculation(decomp_map_src, decomp_map_dst);
        }
    }
}


Remap_weight_of_operator_class *Remap_weight_of_operator_class::generate_parallel_remap_weights(Remap_grid_class **remap_related_decomp_grids, 
                                                                                                 Remap_grid_class **decomp_original_grids, 
                                                                                                 int **global_cells_local_indexes_in_decomps,
                                                                                                 int & field_data_grids_iter,
                                                                                                 Remap_weight_of_strategy_class *parallel_remap_weights_of_strategy)
{
    int i, j, k;
    Remap_weight_of_operator_instance_class *parallel_remap_weights_of_operator_instance;
    long remap_beg_iter, remap_end_iter, global_field_array_offset, local_field_array_offset;
	Remap_grid_class *field_data_grid_src, *field_data_grid_dst, *operator_grid_src, *operator_grid_dst;


    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Remap_weight_of_operator has %ld instances", remap_weights_of_operator_instances.size());

    for (i = 0; i < remap_weights_of_operator_instances.size(); i ++) {
        remap_beg_iter = remap_weights_of_operator_instances[i]->remap_beg_iter;
        if (remap_weights_of_operator_instances[i]->remap_end_iter != -1)
            remap_end_iter = remap_weights_of_operator_instances[i]->remap_end_iter;
        else {
            if (i+1 < remap_weights_of_operator_instances.size())
                remap_end_iter = remap_weights_of_operator_instances[i+1]->remap_beg_iter;
            else remap_end_iter = remap_weights_of_operator_instances[i]->get_field_data_grid_src()->get_grid_size()/remap_weights_of_operator_instances[i]->get_operator_grid_src()->get_grid_size();
        }

        if (remap_weights_of_operator_instances[i]->get_operator_grid_src()->get_is_sphere_grid()) {
            parallel_remap_weights_of_operator_instance = remap_weights_of_operator_instances[i]->generate_parallel_remap_weights(decomp_original_grids, global_cells_local_indexes_in_decomps);
            field_data_grid_src = remap_related_decomp_grids[field_data_grids_iter+0];
            field_data_grid_dst = remap_related_decomp_grids[field_data_grids_iter+1];
            operator_grid_src = remap_related_decomp_grids[field_data_grids_iter+2];
            operator_grid_dst = remap_related_decomp_grids[field_data_grids_iter+3];
            parallel_remap_weights_of_operator_instance->remap_beg_iter = remap_weights_of_operator_instances[i]->remap_beg_iter; 
            parallel_remap_weights_of_operator_instance->remap_end_iter = remap_weights_of_operator_instances[i]->remap_end_iter;
            parallel_remap_weights_of_strategy->add_remap_weight_of_operator_instance(parallel_remap_weights_of_operator_instance, field_data_grid_src, field_data_grid_dst, parallel_remap_weights_of_operator_instance->get_original_remap_operator(), operator_grid_src, operator_grid_dst);
            EXECUTION_REPORT(REPORT_ERROR, -1, remap_weights_of_operator_instances[i]->remap_beg_iter == 0, "C-Coupler error3 in generate_parallel_remap_weights of Remap_weight_of_strategy_class\n");
        }
        else {
            for (j = remap_beg_iter; j < remap_end_iter; j ++) {
                global_field_array_offset = j;
				if (remap_weights_of_operator_instances[i]->get_field_data_grid_src()->have_overlap_with_grid(decomp_original_grids[1])) // horizontal remapping is at the first constantly
                    local_field_array_offset = global_cells_local_indexes_in_decomps[1][global_field_array_offset];
                else if (remap_weights_of_operator_instances[i]->get_field_data_grid_src()->have_overlap_with_grid(decomp_original_grids[0])) {
                    EXECUTION_REPORT(REPORT_ERROR, -1, false);
                    local_field_array_offset = global_cells_local_indexes_in_decomps[0][global_field_array_offset];
                }
                else EXECUTION_REPORT(REPORT_ERROR, -1, false, "C-Coupler error4 in generate_parallel_remap_weights of Remap_weight_of_strategy_class\n");
                if (local_field_array_offset == -1)
                    continue;
                parallel_remap_weights_of_operator_instance = remap_weights_of_operator_instances[i]->generate_parallel_remap_weights(decomp_original_grids, global_cells_local_indexes_in_decomps);
                field_data_grid_src = remap_related_decomp_grids[field_data_grids_iter+0];
                field_data_grid_dst = remap_related_decomp_grids[field_data_grids_iter+1];
                operator_grid_src = remap_related_decomp_grids[field_data_grids_iter+2];
                operator_grid_dst = remap_related_decomp_grids[field_data_grids_iter+3];
                parallel_remap_weights_of_operator_instance->remap_beg_iter = local_field_array_offset; 
                parallel_remap_weights_of_operator_instance->remap_end_iter = local_field_array_offset+1;                
                parallel_remap_weights_of_strategy->add_remap_weight_of_operator_instance(parallel_remap_weights_of_operator_instance, field_data_grid_src, field_data_grid_dst, parallel_remap_weights_of_operator_instance->get_original_remap_operator(), operator_grid_src, operator_grid_dst);
            }
        }
    }
	field_data_grids_iter += 4;
}


Remap_weight_of_operator_class::~Remap_weight_of_operator_class()
{
    for (int i = 0; i < remap_weights_of_operator_instances.size(); i ++)
        delete remap_weights_of_operator_instances[i];
}


void Remap_weight_of_operator_class::do_remap(int comp_id, Remap_grid_data_class *field_data_src, Remap_grid_data_class *field_data_dst)
{

    double *data_value_src, *data_value_dst;
    int i, j, k;
    long remap_beg_iter, remap_end_iter;
    long field_array_offset;
    long field_data_size_src, field_data_size_dst;

    
    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, field_data_src->get_coord_value_grid()->is_similar_grid_with(field_data_grid_src), "C-Coupler error1 in do_remap of Remap_weight_of_operator_class");
    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, field_data_dst->get_coord_value_grid()->is_similar_grid_with(field_data_grid_dst), "C-Coupler error2 in do_remap of Remap_weight_of_operator_class");

    if (comp_id != -1)          
        comp_comm_group_mgt_mgr->get_global_node_of_local_comp(comp_id,false,"")->get_performance_timing_mgr()->performance_timing_start(TIMING_TYPE_COMPUTATION, -1, -1, "interchange data");
    field_data_src->interchange_grid_data(field_data_grid_src);
    field_data_dst->interchange_grid_data(field_data_grid_dst);
    if (comp_id != -1)          
        comp_comm_group_mgt_mgr->get_global_node_of_local_comp(comp_id,false,"")->get_performance_timing_mgr()->performance_timing_stop(TIMING_TYPE_COMPUTATION, -1, -1, "interchange data");

    field_data_size_src = field_data_src->get_grid_data_field()->read_data_size;
    field_data_size_dst = field_data_dst->get_grid_data_field()->read_data_size;

    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, !is_remap_weight_empty(), "Software error in Remap_weight_of_operator_class::do_remap: empty remap weights");
    
    for (i = 0; i < remap_weights_of_operator_instances.size(); i ++) {
        remap_beg_iter = remap_weights_of_operator_instances[i]->remap_beg_iter;
        if (remap_weights_of_operator_instances[i]->remap_end_iter != -1)
            remap_end_iter = remap_weights_of_operator_instances[i]->remap_end_iter;
        else if (i+1 < remap_weights_of_operator_instances.size())
                remap_end_iter = remap_weights_of_operator_instances[i+1]->remap_beg_iter;
        else remap_end_iter = field_data_grid_src->get_grid_size()/operator_grid_src->get_grid_size();
        for (j = remap_beg_iter; j < remap_end_iter; j ++) {
            field_array_offset = j;
            if (report_error_enabled) {
                EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, field_array_offset >= 0 && (field_array_offset+1)*remap_weights_of_operator_instances[i]->get_operator_grid_src()->get_grid_size() <= field_data_size_src,
                                 "remap software error4 in do_remap of Remap_weight_of_strategy_class");
                EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, field_array_offset >= 0 && (field_array_offset+1)*remap_weights_of_operator_instances[i]->get_operator_grid_dst()->get_grid_size() <= field_data_size_dst,
                                  "remap software error5 in do_remap of Remap_weight_of_strategy_class");
            }    
            data_value_src = ((double*) field_data_src->get_grid_data_field()->data_buf) + field_array_offset*remap_weights_of_operator_instances[i]->get_operator_grid_src()->get_grid_size();
            data_value_dst = ((double*) field_data_dst->get_grid_data_field()->data_buf) + field_array_offset*remap_weights_of_operator_instances[i]->get_operator_grid_dst()->get_grid_size();
            EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, remap_weights_of_operator_instances[i]->duplicated_remap_operator != NULL, "C-Coupler error3 in do_remap of Remap_weight_of_operator_class %s", remap_weights_of_operator_instances[i]->get_operator_grid_src()->get_grid_name());
            remap_weights_of_operator_instances[i]->duplicated_remap_operator->do_remap_values_caculation(data_value_src, data_value_dst, field_data_dst->get_grid_data_field()->required_data_size);
        }
    }
}


void Remap_weight_of_operator_class::add_remap_weight_of_operator_instance(Remap_weight_of_operator_instance_class *operator_instance)
{
    remap_weights_of_operator_instances.push_back(operator_instance);
	remap_weights_of_operator_instances[remap_weights_of_operator_instances.size()-1]->set_remap_weight_of_operator(this);
}


void Remap_weight_of_operator_class::renew_vertical_remap_weights(Remap_grid_class *runtime_remap_grid_src, Remap_grid_class *runtime_remap_grid_dst)
{
    long i;
    Remap_grid_data_class *lev_center_field_in_3D_src_grid = NULL, *lev_center_field_in_3D_dst_grid = NULL;
    Remap_operator_grid *runtime_remap_operator_grid_src = NULL, *runtime_remap_operator_grid_dst = NULL;
    Remap_operator_basis *new_remap_operator;
    double *lev_center_values_in_3D_src_grid = NULL, *lev_center_values_in_3D_dst_grid = NULL;
    long lev_grid_size_src, lev_grid_size_dst, offset;

    
    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, runtime_remap_grid_src->get_num_dimensions() == 1 && runtime_remap_grid_src->has_grid_coord_label(COORD_LABEL_LEV) && runtime_remap_grid_dst->get_num_dimensions() == 1 && runtime_remap_grid_dst->has_grid_coord_label(COORD_LABEL_LEV),
                     "C-Coupler error1 in renew_vertical_remap_weights of Remap_weight_of_operator_class");
    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, runtime_remap_grid_src->get_num_dimensions() == 1 && runtime_remap_grid_src->has_grid_coord_label(COORD_LABEL_LEV) && runtime_remap_grid_dst->get_num_dimensions() == 1 && runtime_remap_grid_dst->has_grid_coord_label(COORD_LABEL_LEV),
                     "C-Coupler error2 in renew_vertical_remap_weights of Remap_weight_of_operator_class");
    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, original_remap_operator->get_src_grid()->get_num_dimensions() == 1 && original_remap_operator->get_src_grid()->has_grid_coord_label(COORD_LABEL_LEV) && original_remap_operator->get_dst_grid()->get_num_dimensions() == 1 && original_remap_operator->get_dst_grid()->has_grid_coord_label(COORD_LABEL_LEV),
                     "C-Coupler error3 in renew_vertical_remap_weights of Remap_weight_of_operator_class");
    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, original_remap_operator->get_src_grid()->get_num_dimensions() == 1 && original_remap_operator->get_src_grid()->has_grid_coord_label(COORD_LABEL_LEV) && original_remap_operator->get_dst_grid()->get_num_dimensions() == 1 && original_remap_operator->get_dst_grid()->has_grid_coord_label(COORD_LABEL_LEV),
                     "C-Coupler error4 in renew_vertical_remap_weights of Remap_weight_of_operator_class");
    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, field_data_grid_src->is_sigma_grid() || field_data_grid_dst->is_sigma_grid() || field_data_grid_src->does_use_V3D_level_coord() || field_data_grid_dst->does_use_V3D_level_coord(), "C-Coupler error4 in renew_vertical_remap_weights of Remap_weight_of_operator_class");
    
    if (field_data_grid_src->is_sigma_grid() || field_data_grid_src->does_use_V3D_level_coord()) {
        EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, runtime_remap_grid_src->is_subset_of_grid(field_data_grid_src), "C-Coupler error5 in renew_vertical_remap_weights of Remap_weight_of_operator_class"); 
        lev_center_field_in_3D_src_grid = field_data_grid_src->get_unique_center_field();
        lev_center_field_in_3D_src_grid->interchange_grid_data(field_data_grid_src);
        lev_center_values_in_3D_src_grid = (double*) lev_center_field_in_3D_src_grid->get_grid_data_field()->data_buf;
    }
    if (field_data_grid_dst->is_sigma_grid() || field_data_grid_dst->does_use_V3D_level_coord()) {
        EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, runtime_remap_grid_dst->is_subset_of_grid(field_data_grid_dst), "C-Coupler error6 in renew_vertical_remap_weights of Remap_weight_of_operator_class"); 
        lev_center_field_in_3D_dst_grid = field_data_grid_dst->get_unique_center_field();
        lev_center_field_in_3D_dst_grid->interchange_grid_data(field_data_grid_dst);
        lev_center_values_in_3D_dst_grid = (double*) lev_center_field_in_3D_dst_grid->get_grid_data_field()->data_buf;
    }
    lev_grid_size_src = runtime_remap_grid_src->get_grid_size();
    lev_grid_size_dst = runtime_remap_grid_dst->get_grid_size();

    for (i = 0; i < remap_weights_of_operator_instances.size(); i ++) {
        EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, remap_weights_of_operator_instances[i]->get_original_remap_operator() != NULL, "C-Coupler error7 in renew_vertical_remap_weights of Remap_weight_of_operator_class"); 
        new_remap_operator = remap_weights_of_operator_instances[i]->get_original_remap_operator()->duplicate_remap_operator(true);
        new_remap_operator->set_src_grid(runtime_remap_grid_src);
        new_remap_operator->set_dst_grid(runtime_remap_grid_dst);
        offset = remap_weights_of_operator_instances[i]->remap_beg_iter;
        if (lev_center_field_in_3D_src_grid != NULL) {
            EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, offset >= 0 && offset*lev_grid_size_src+lev_grid_size_src<= lev_center_field_in_3D_src_grid->get_grid_data_field()->required_data_size, "C-Coupler error7 in renew_vertical_remap_weights of Remap_weight_of_operator_class"); 
        }    
        if (lev_center_values_in_3D_src_grid != NULL)
            runtime_remap_grid_src->renew_lev_grid_coord_values(lev_center_values_in_3D_src_grid+offset*lev_grid_size_src, NULL);
        if (lev_center_values_in_3D_dst_grid != NULL) {
            runtime_remap_grid_dst->renew_lev_grid_coord_values(lev_center_values_in_3D_dst_grid+offset*lev_grid_size_dst, NULL);
        }
        if (runtime_remap_operator_grid_src == NULL) {
            runtime_remap_operator_grid_src = new Remap_operator_grid(runtime_remap_grid_src, new_remap_operator, true, false);
            runtime_remap_operator_grid_dst = new Remap_operator_grid(runtime_remap_grid_dst, new_remap_operator, false, false);
            current_runtime_remap_operator_grid_src = runtime_remap_operator_grid_src;
            current_runtime_remap_operator_grid_dst = runtime_remap_operator_grid_dst;
        }
        if (lev_center_values_in_3D_src_grid != NULL)
            runtime_remap_operator_grid_src->update_operator_grid_data();
        if (lev_center_values_in_3D_dst_grid != NULL)
            runtime_remap_operator_grid_dst->update_operator_grid_data();
        current_runtime_remap_operator = new_remap_operator;
        new_remap_operator->calculate_remap_weights();
//        new_remap_operator->get_remap_weights_group(0)->compare_to_another_sparse_matrix(remap_weights_of_operator_instances[i]->duplicated_remap_operator->get_remap_weights_group(0));
//        new_remap_operator->get_remap_weights_group(1)->compare_to_another_sparse_matrix(remap_weights_of_operator_instances[i]->duplicated_remap_operator->get_remap_weights_group(1));
//        new_remap_operator->get_remap_weights_group(2)->compare_to_another_sparse_matrix(remap_weights_of_operator_instances[i]->duplicated_remap_operator->get_remap_weights_group(2));
//        new_remap_operator->get_remap_weights_group(3)->compare_to_another_sparse_matrix(remap_weights_of_operator_instances[i]->duplicated_remap_operator->get_remap_weights_group(3));
		if (remap_weights_of_operator_instances[i]->duplicated_remap_operator != NULL)
	        delete remap_weights_of_operator_instances[i]->duplicated_remap_operator;
        remap_weights_of_operator_instances[i]->duplicated_remap_operator = new_remap_operator;
    }

    if (runtime_remap_operator_grid_src != NULL) {
        delete runtime_remap_operator_grid_src;
        delete runtime_remap_operator_grid_dst;
    }
    
    empty_remap_weight = false;
}


void Remap_weight_of_operator_class::write_overall_remapping_weights(int comp_id)
{
	char default_wgt_file_name[NAME_STR_SIZE], full_default_wgt_file_name[NAME_STR_SIZE*2];
	Remap_operator_basis *overall_remap_operator;
	Comp_comm_group_mgt_node *comp_node = comp_comm_group_mgt_mgr->search_global_node(comp_id);

	
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, remap_weights_of_operator_instances.size() == 1, "Software error in Remap_weight_of_operator_class::write_overall_remapping_weights");
	sprintf(default_wgt_file_name, "DEFAULT_WGT_of___%s___FROM___%s___TO___%s___AT___%s.nc", remap_weights_of_operator_instances[0]->get_original_remap_operator()->get_operator_name(), remap_weights_of_operator_instances[0]->get_operator_grid_src()->get_grid_name(), remap_weights_of_operator_instances[0]->get_operator_grid_dst()->get_grid_name(), comp_node->get_full_name());
	sprintf(full_default_wgt_file_name, "%s/%s", comp_comm_group_mgt_mgr->get_internal_remapping_weights_dir(), default_wgt_file_name);
	EXECUTION_REPORT_LOG(REPORT_LOG, comp_id, true, "The default H2D weight file name is \"%s\"", default_wgt_file_name);
	overall_remap_operator = remap_weights_of_operator_instances[0]->duplicated_remap_operator->gather(comp_id);
	if (comp_node->get_current_proc_local_id() == 0) {
		Remap_weight_of_operator_instance_class *overall_remap_weight_of_operator_instance = new Remap_weight_of_operator_instance_class(operator_grid_src, operator_grid_dst, 0, remap_weights_of_operator_instances[0]->get_original_remap_operator(), overall_remap_operator);
		Remap_weight_of_operator_class *overall_remap_weight_of_operator = new Remap_weight_of_operator_class(operator_grid_src, operator_grid_dst, original_remap_operator, operator_grid_src, operator_grid_dst);
		overall_remap_weight_of_operator->remap_weights_of_operator_instances.push_back(overall_remap_weight_of_operator_instance);
		Remap_weight_of_strategy_class *overall_remap_weights = new Remap_weight_of_strategy_class("overall_remapping_weights", NULL, operator_grid_src, operator_grid_dst, NULL, false, comp_id);
		overall_remap_weights->add_remap_weights_of_operator(overall_remap_weight_of_operator);
		IO_netcdf *io_netcdf = new IO_netcdf(default_wgt_file_name, full_default_wgt_file_name, "w", true);
		int last_execution_phase_number = execution_phase_number;
		execution_phase_number = 1;
		io_netcdf->write_remap_weights(overall_remap_weights);
		execution_phase_number = last_execution_phase_number;
		delete io_netcdf;
		delete overall_remap_weights;
	}
}


void Remap_weight_of_strategy_class::initialize_object()
{
    dynamic_vertical_remapping_weights_src = false;
    dynamic_vertical_remapping_weights_dst = false;
    public_remap_weights_of_operators = false;
    num_field_data_grids_in_remapping_process = 0;
}


Remap_weight_of_strategy_class::Remap_weight_of_strategy_class(const char *object_name, const char *remap_strategy_name, 
                                                               const char *data_grid_name_src, const char *data_grid_name_dst,
                                                               const char *input_IO_file_name, const char *weight_IO_format,
                                                               bool read_from_io)
{
    initialize_object();
    strcpy(this->object_name, object_name);
    remap_strategy = remap_strategy_manager->search_remap_strategy(remap_strategy_name);
    data_grid_src = remap_grid_manager->search_remap_grid_with_grid_name(data_grid_name_src);
    data_grid_dst = remap_grid_manager->search_remap_grid_with_grid_name(data_grid_name_dst);

    EXECUTION_REPORT(REPORT_ERROR, -1, remap_strategy != NULL && data_grid_src != NULL && data_grid_dst != NULL, "C-Coupler error in Remap_weight_of_strategy_class::Remap_weight_of_strategy_class");

    generate_remapping_related_grids();

    if (!read_from_io)
        remap_strategy->calculate_remapping_weights(this, NULL, -1);
    else {
        if (words_are_the_same(weight_IO_format, "SCRIP")) 
            ((IO_netcdf*) (io_manager->search_IO_object(input_IO_file_name)))->read_remap_weights(this, remap_strategy, is_master_process_in_computing_node);
        else ((IO_binary*) (io_manager->search_IO_object(input_IO_file_name)))->read_remap_weights(this, remap_strategy, is_master_process_in_computing_node);
    }

    EXECUTION_REPORT(REPORT_ERROR, -1, data_grid_src->get_num_dimensions() == data_grid_dst->get_num_dimensions(), 
                     "grid %s and %s must have the same number of dimensions\n", data_grid_name_src, data_grid_name_dst);
}


Remap_weight_of_strategy_class::Remap_weight_of_strategy_class(const char *object_name, Remap_strategy_class *remap_strategy, 
                                                               Remap_grid_class *data_grid_src, Remap_grid_class *data_grid_dst, const char *H2D_remapping_wgt_file, bool calculate_wgts, int wgt_cal_comp_id)
{
    initialize_object();
    strcpy(this->object_name, object_name);
    this->remap_strategy = remap_strategy;
    this->data_grid_src = data_grid_src;
    this->data_grid_dst = data_grid_dst;

    EXECUTION_REPORT(REPORT_ERROR, -1, (!calculate_wgts || remap_strategy != NULL) && data_grid_src != NULL && data_grid_dst != NULL, "C-Coupler error in Remap_weight_of_strategy_class::Remap_weight_of_strategy_class");

	if (calculate_wgts) {
		EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "before generate_remapping_related_grids");
	    generate_remapping_related_grids();
		EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "After generate_remapping_related_grids, before generate_remapping_related_grids");
	    remap_strategy->calculate_remapping_weights(this, H2D_remapping_wgt_file, wgt_cal_comp_id);
	}
}


int Remap_weight_of_strategy_class::generate_remapping_related_grids()
{
    int i, j;
    Remap_grid_class *current_remap_src_data_grid, *current_remap_src_data_grid_interchanged, *current_remap_dst_data_grid, *existing_grid;
    Remap_grid_class *runtime_mask_sub_grids_src[256], *runtime_mask_sub_grids_dst[256];
    int num_runtime_mask_sub_grids_src, num_runtime_mask_sub_grids_dst;
    Remap_grid_data_class *runtime_mask_src, *runtime_mask_dst;
    Remap_grid_class *leaf_grids[256];
    int num_leaf_grids;

    
    num_field_data_grids_in_remapping_process = 0;
    field_data_grids_in_remapping_process[num_field_data_grids_in_remapping_process] = data_grid_src;
    runtime_mask_fields_in_remapping_process[num_field_data_grids_in_remapping_process++] = NULL;
    current_remap_src_data_grid = data_grid_src;
    for (i = 0; i < remap_strategy->get_num_remap_operator(); i ++) {
        if (i == 0)
            current_remap_src_data_grid->compute_remap_field_data_runtime_mask(data_grid_src,
                                                                               runtime_mask_sub_grids_src,
                                                                               &num_runtime_mask_sub_grids_src,
                                                                               &runtime_mask_src);
        else current_remap_src_data_grid->compute_remap_field_data_runtime_mask(NULL,
                                                                                runtime_mask_sub_grids_src,
                                                                                &num_runtime_mask_sub_grids_src,
                                                                                &runtime_mask_src);
        current_remap_src_data_grid->generate_interchange_grids(remap_strategy->get_remap_operator(i)->get_src_grid(), 
                                                                &current_remap_src_data_grid_interchanged, 
                                                                runtime_mask_sub_grids_src,
                                                                num_runtime_mask_sub_grids_src);
        current_remap_dst_data_grid = new Remap_grid_class(current_remap_src_data_grid_interchanged, 
                                                           remap_strategy->get_remap_operator(i)->get_src_grid(), 
                                                           remap_strategy->get_remap_operator(i)->get_dst_grid(),
                                                           remap_strategy->get_remap_operator(i)->get_is_operator_regridding()); 
        existing_grid = remap_grid_manager->search_same_remap_grid(current_remap_dst_data_grid);
        if (existing_grid != NULL) {
            delete current_remap_dst_data_grid;
            current_remap_dst_data_grid = existing_grid;
        }

        if (i == remap_strategy->get_num_remap_operator()-1 && !current_remap_dst_data_grid->is_similar_grid_with(data_grid_dst)) {
            data_grid_dst->get_leaf_grids(&num_leaf_grids, leaf_grids, data_grid_dst);
            for (j = 0; j < num_leaf_grids; j ++) 
                EXECUTION_REPORT(REPORT_ERROR, -1, leaf_grids[j]->is_subset_of_grid(current_remap_dst_data_grid), 
                                    "this remap caculation can not get 1D grid %s, which is a sub grid of the grid of destination field data. Please check the %dth remapping operator of the remapping strategy", 
                                    leaf_grids[j]->get_grid_name(), i+1);
        }
        if (i == remap_strategy->get_num_remap_operator()-1)
            current_remap_dst_data_grid->compute_remap_field_data_runtime_mask(data_grid_dst,
                                                                               runtime_mask_sub_grids_dst,
                                                                               &num_runtime_mask_sub_grids_dst,
                                                                               &runtime_mask_dst);
        else current_remap_dst_data_grid->compute_remap_field_data_runtime_mask(NULL,runtime_mask_sub_grids_dst,
                                                                                &num_runtime_mask_sub_grids_dst,
                                                                                &runtime_mask_dst);
        field_data_grids_in_remapping_process[num_field_data_grids_in_remapping_process] = current_remap_src_data_grid_interchanged;
        runtime_mask_fields_in_remapping_process[num_field_data_grids_in_remapping_process++] = runtime_mask_src;
        field_data_grids_in_remapping_process[num_field_data_grids_in_remapping_process] = current_remap_dst_data_grid;
        runtime_mask_fields_in_remapping_process[num_field_data_grids_in_remapping_process++] = runtime_mask_dst;        
        current_remap_src_data_grid = current_remap_dst_data_grid;
    }

    field_data_grids_in_remapping_process[num_field_data_grids_in_remapping_process] = data_grid_dst;
    runtime_mask_fields_in_remapping_process[num_field_data_grids_in_remapping_process++] = NULL;

    return num_field_data_grids_in_remapping_process;
}


void Remap_weight_of_strategy_class::set_basic_fields(const char *object_name, Remap_strategy_class *remap_strategy, Remap_grid_class *data_grid_src, Remap_grid_class *data_grid_dst)
{
        strcpy(this->object_name, object_name);
        this->remap_strategy = remap_strategy;
        this->data_grid_src = data_grid_src;
        this->data_grid_dst = data_grid_dst;
}


bool Remap_weight_of_strategy_class::match_object_name(const char*object_name)
{
    return words_are_the_same(this->object_name, object_name);
}


Remap_weight_of_strategy_class::~Remap_weight_of_strategy_class()
{
    if (!public_remap_weights_of_operators)
        for (int i = 0; i < remap_weights_of_operators.size(); i ++)
            delete remap_weights_of_operators[i];
}


Remap_weight_of_operator_instance_class *Remap_weight_of_strategy_class::add_remap_weight_of_operator_instance(Remap_grid_class *field_data_grid_src, Remap_grid_class *field_data_grid_dst,
                                                                  long remap_beg_iter, Remap_operator_basis *remap_operator)
{
    Remap_weight_of_operator_instance_class *remap_weight_of_operator_instance = new Remap_weight_of_operator_instance_class(field_data_grid_src, field_data_grid_dst, remap_beg_iter, remap_operator);
    add_remap_weight_of_operator_instance(remap_weight_of_operator_instance, field_data_grid_src, field_data_grid_dst, remap_operator, remap_operator->get_src_grid(), remap_operator->get_dst_grid());
    return remap_weight_of_operator_instance;
}


void Remap_weight_of_strategy_class::add_remap_weight_of_operator_instance(Remap_weight_of_operator_instance_class *weight_of_operator_instance, Remap_grid_class *field_data_grid_src, Remap_grid_class *field_data_grid_dst, Remap_operator_basis *original_remap_operator,
                                                                           Remap_grid_class *operator_grid_src, Remap_grid_class *operator_grid_dst) 
{ 
    if (remap_weights_of_operators.size() == 0 || remap_weights_of_operators[remap_weights_of_operators.size()-1]->field_data_grid_src != field_data_grid_src) {
        Remap_weight_of_operator_class *remap_weight_of_operator = new Remap_weight_of_operator_class(field_data_grid_src, field_data_grid_dst, original_remap_operator,
                                                                                                      operator_grid_src, operator_grid_dst);
        remap_weights_of_operators.push_back(remap_weight_of_operator);
    }
    remap_weights_of_operators[remap_weights_of_operators.size()-1]->add_remap_weight_of_operator_instance(weight_of_operator_instance);
}


void Remap_weight_of_strategy_class::do_remap(int comp_id, Remap_grid_data_class *field_data_src, Remap_grid_data_class *field_data_dst)
{
    Remap_grid_class *sized_sub_grids[256];
    Remap_grid_class *field_data_grid_src, *field_data_grid_dst;
    Remap_grid_data_class *tmp_field_data_src, *tmp_field_data_dst;
    Remap_operator_basis *current_remap_operator;
    double *data_value_src, *data_value_dst;
    bool is_last_remap_operator;
    int i, j, k;
    long remap_beg_iter, remap_end_iter;
    long field_array_offset;
    Remap_grid_class *original_lev_grid_src = NULL, *original_lev_grid_dst = NULL;

    
    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, field_data_src->get_coord_value_grid()->is_similar_grid_with(data_grid_src),
                 "the grid of field data \"%s\" can not match the src grid of remap weight object \"%s\"",
                 field_data_src->get_grid_data_field()->field_name_in_application, object_name);
	if (!words_are_the_same(field_data_dst->get_grid_data_field()->field_name_in_application, V3D_GRID_3D_LEVEL_FIELD_NAME))
    	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, field_data_dst->get_coord_value_grid()->is_similar_grid_with(data_grid_dst),
        	         "the grid of field data \"%s\" can not match the dst grid of remap weight object \"%s\"",
            	     field_data_dst->get_grid_data_field()->field_name_in_application, object_name);

    field_data_src->transfer_field_attributes_to_another(field_data_dst);
    if (!field_data_dst->have_data_content())
        field_data_dst->get_grid_data_field()->initialize_to_fill_value();

    tmp_field_data_dst = field_data_src;
    tmp_field_data_src = NULL;
    for (i = 0; i < remap_weights_of_operators.size(); i ++) {
        if (comp_id != -1)        
            comp_comm_group_mgt_mgr->get_global_node_of_local_comp(comp_id,false,"")->get_performance_timing_mgr()->performance_timing_start(TIMING_TYPE_COMPUTATION, -1, -1, remap_weights_of_operators[i]->get_original_remap_operator()->get_operator_name());
        if (tmp_field_data_src != NULL && tmp_field_data_src != field_data_src)
            delete tmp_field_data_src;
        tmp_field_data_src = tmp_field_data_dst;
        if (i == remap_weights_of_operators.size()-1 || remap_weights_of_operators[i]->field_data_grid_dst->is_similar_grid_with(field_data_dst->get_coord_value_grid())) {
            tmp_field_data_dst = field_data_dst;
            EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, remap_weights_of_operators[i]->field_data_grid_dst->is_similar_grid_with(tmp_field_data_dst->get_coord_value_grid()), 
                         "remap software error1 in do_remap of Remap_weight_of_strategy_class\n");
			if (i != remap_weights_of_operators.size()-1 && remap_weights_of_operators[i]->field_data_grid_dst->is_similar_grid_with(field_data_dst->get_coord_value_grid()))
				EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, words_are_the_same(field_data_dst->get_grid_data_field()->field_name_in_application, V3D_GRID_3D_LEVEL_FIELD_NAME), "Software error in Remap_weight_of_strategy_class::do_remap");
        }    
        else tmp_field_data_dst = field_data_src->duplicate_grid_data_field(remap_weights_of_operators[i]->field_data_grid_dst, 1, false, false);
        remap_weights_of_operators[i]->do_remap(comp_id, tmp_field_data_src, tmp_field_data_dst);
        if (comp_id != -1)        
            comp_comm_group_mgt_mgr->get_global_node_of_local_comp(comp_id,false,"")->get_performance_timing_mgr()->performance_timing_stop(TIMING_TYPE_COMPUTATION, -1, -1, remap_weights_of_operators[i]->get_original_remap_operator()->get_operator_name());
		if (i != remap_weights_of_operators.size()-1 && tmp_field_data_dst == field_data_dst) {
			EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, words_are_the_same(field_data_dst->get_grid_data_field()->field_name_in_application, V3D_GRID_3D_LEVEL_FIELD_NAME), "Software error in Remap_weight_of_strategy_class::do_remap");
			break;
		}
    }

    if (tmp_field_data_src != NULL && tmp_field_data_src != field_data_src)
        delete tmp_field_data_src;
    if (comp_id != -1)          
        comp_comm_group_mgt_mgr->get_global_node_of_local_comp(comp_id,false,"")->get_performance_timing_mgr()->performance_timing_start(TIMING_TYPE_COMPUTATION, -1, -1, "interchange data");
    field_data_src->interchange_grid_data(field_data_src->get_coord_value_grid());
    field_data_dst->interchange_grid_data(field_data_dst->get_coord_value_grid());
    field_data_dst->get_grid_data_field()->read_data_size = field_data_dst->get_grid_data_field()->required_data_size;
    if (comp_id != -1)          
        comp_comm_group_mgt_mgr->get_global_node_of_local_comp(comp_id,false,"")->get_performance_timing_mgr()->performance_timing_stop(TIMING_TYPE_COMPUTATION, -1, -1, "interchange data");
}


void Remap_weight_of_strategy_class::calculate_src_decomp(Remap_grid_class *grid_src, Remap_grid_class *grid_dst, long *decomp_map_src, const long *decomp_map_dst)
{
    long i, j;
    Remap_weight_of_operator_class *H2D_operator_weights = NULL;


    for (i = 0; i < grid_src->get_grid_size(); i ++)
        decomp_map_src[i] = 0;
    
     for (i = 0; i < remap_weights_of_operators.size(); i ++) {
        if (remap_weights_of_operators[i]->get_original_remap_operator()->get_src_grid()->has_grid_coord_label(COORD_LABEL_LON) || remap_weights_of_operators[i]->get_original_remap_operator()->get_src_grid()->has_grid_coord_label(COORD_LABEL_LON)) {
            EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, remap_weights_of_operators[i]->get_original_remap_operator()->get_src_grid()->get_is_sphere_grid(), "Software error in Remap_weight_of_strategy_class::calculate_src_decomp: 1-D remapping operator for lon or lat");
        }    
        if (remap_weights_of_operators[i]->get_original_remap_operator()->get_src_grid()->get_is_sphere_grid()) {
            EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, H2D_operator_weights == NULL, "Software error in Remap_weight_of_strategy_class::calculate_src_decomp: multiple H2D remapping operator");
            H2D_operator_weights = remap_weights_of_operators[i];
        }
     }
    if (H2D_operator_weights == NULL) {
        EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, grid_src == grid_dst, "Software error in Remap_weight_of_strategy_class::calculate_src_decomp: src and dst grids are different when there are no H2D operator");
        for (i = 0; i < grid_src->get_grid_size(); i ++)
            decomp_map_src[i] = decomp_map_dst[i];
        return;
    }    

    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, grid_src == H2D_operator_weights->get_original_remap_operator()->get_src_grid() && grid_dst == H2D_operator_weights->get_original_remap_operator()->get_dst_grid(), "Software error in Remap_weight_of_strategy_class::calculate_src_decomp: H2D operator grids do not match decomp grids");

    H2D_operator_weights->calculate_src_decomp(decomp_map_src, decomp_map_dst);
}


void Remap_weight_of_strategy_class::add_remap_related_grid(std::vector<std::pair<Remap_grid_class *, bool> > &all_remap_related_grids_info, Remap_grid_class *remap_grid, bool on_dst_decomp)
{
	std::pair<Remap_grid_class *, bool> remap_related_grid_info;
	
	remap_related_grid_info.first = remap_grid;
	remap_related_grid_info.second = on_dst_decomp;
	all_remap_related_grids_info.push_back(remap_related_grid_info);
}


void Remap_weight_of_strategy_class::get_remap_related_grids(std::vector<std::pair<Remap_grid_class *, bool> > &all_remap_related_grids_info)
{   
	
	bool after_H2D_remapping = false;

	add_remap_related_grid(all_remap_related_grids_info, data_grid_src, false);
	add_remap_related_grid(all_remap_related_grids_info, data_grid_dst, true);
    for (int i = 0; i < remap_weights_of_operators.size(); i ++) {
		add_remap_related_grid(all_remap_related_grids_info, remap_weights_of_operators[i]->get_field_data_grid_src(), after_H2D_remapping);
		add_remap_related_grid(all_remap_related_grids_info, remap_weights_of_operators[i]->get_field_data_grid_dst(), remap_weights_of_operators[i]->get_operator_grid_src()->get_is_sphere_grid() || after_H2D_remapping);
		add_remap_related_grid(all_remap_related_grids_info, remap_weights_of_operators[i]->get_operator_grid_src(), after_H2D_remapping);
		if (remap_weights_of_operators[i]->get_operator_grid_src()->get_is_sphere_grid())
			after_H2D_remapping = true;
		add_remap_related_grid(all_remap_related_grids_info, remap_weights_of_operators[i]->get_operator_grid_dst(), after_H2D_remapping);	
    }
}


Remap_weight_of_strategy_class *Remap_weight_of_strategy_class::generate_parallel_remap_weights(Remap_grid_class **remap_related_decomp_grids,
                                                                                                 Remap_grid_class **decomp_original_grids, 
                                                                                                 int **global_cells_local_indexes_in_decomps)
{
    int i, j, k, field_data_grids_iter = 0;
    Remap_operator_basis *current_remap_operator;
    Remap_weight_of_strategy_class *parallel_remap_weights_of_strategy = new Remap_weight_of_strategy_class;
    Remap_weight_of_operator_instance_class *parallel_remap_weights_of_operator_instance;
    Remap_grid_class *sized_sub_grids[256];
    long remap_beg_iter, remap_end_iter, global_field_array_offset, local_field_array_offset;


    EXECUTION_REPORT(REPORT_ERROR, -1, decomp_original_grids[0]->is_subset_of_grid(this->data_grid_src) && decomp_original_grids[1]->is_subset_of_grid(this->data_grid_dst),
                 "C-Coupler error1 in generate_parallel_remap_weights of Remap_weight_of_strategy_class\n");
    EXECUTION_REPORT(REPORT_ERROR, -1, this->data_grid_src->get_num_dimensions() >= 2 && this->data_grid_src->get_num_dimensions() <= 3,
                 "C-Coupler error2 in generate_parallel_remap_weights of Remap_weight_of_strategy_class\n");

    strcpy(parallel_remap_weights_of_strategy->object_name, this->object_name);
    parallel_remap_weights_of_strategy->remap_strategy = this->remap_strategy;
    parallel_remap_weights_of_strategy->data_grid_src = remap_related_decomp_grids[field_data_grids_iter++];
    parallel_remap_weights_of_strategy->data_grid_dst = remap_related_decomp_grids[field_data_grids_iter++];

    for (i = 0; i < remap_weights_of_operators.size(); i ++)
        remap_weights_of_operators[i]->generate_parallel_remap_weights(remap_related_decomp_grids, decomp_original_grids, global_cells_local_indexes_in_decomps, field_data_grids_iter, parallel_remap_weights_of_strategy);

    for (i = 0; i < remap_weights_of_operators.size(); i ++)
        if (remap_weights_of_operators[i]->is_remap_weight_empty())
            parallel_remap_weights_of_strategy->remap_weights_of_operators[i]->mark_empty_remap_weight();

    return parallel_remap_weights_of_strategy;
}


void Remap_weight_of_strategy_class::write_data_into_array(void *data, int data_size, char **array, long &current_array_size, long &max_array_size)
{
    char *new_array;

    
    if (data_size + current_array_size > max_array_size) {
        max_array_size = (data_size+current_array_size)*2;
        new_array = new char [max_array_size];
        for (long i = 0; i < current_array_size; i ++)
            new_array[i] = (*array)[i];
        delete [] (*array);
        (*array) = new_array;
    }

    for (long i = 0; i < data_size; i ++)
        (*array)[current_array_size+i] = ((char*)data)[i];
    current_array_size += data_size;
}


void Remap_weight_of_strategy_class::write_grid_info_into_array(Remap_grid_class *grid, bool consider_area_or_volumn, char **array, long &current_array_size, long &max_array_size)
{
    long grid_size;
    int grid_num_dimensions, i, id, num_leaf_grids, tmp_int_value;
    Remap_grid_class *leaf_grids[256];
    
    
    grid_size = grid->get_grid_size();
    write_data_into_array(&grid_size, sizeof(long), array, current_array_size, max_array_size);
    grid_num_dimensions = grid->get_num_dimensions();
    write_data_into_array(&grid_num_dimensions, sizeof(int), array, current_array_size, max_array_size);
    grid->get_leaf_grids(&num_leaf_grids, leaf_grids, grid);
    for (i = 0; i < num_leaf_grids; i ++) {
        if (words_are_the_same(leaf_grids[i]->get_coord_label(), COORD_LABEL_LON))
            id = 1;
        else if (words_are_the_same(leaf_grids[i]->get_coord_label(), COORD_LABEL_LAT))
            id = 2;
        else if (words_are_the_same(leaf_grids[i]->get_coord_label(), COORD_LABEL_LEV))
            id = 3;
        else if (words_are_the_same(leaf_grids[i]->get_coord_label(), COORD_LABEL_TIME))
            id = 4;
        write_data_into_array(&id, sizeof(int), array, current_array_size, max_array_size);
    }

    if (consider_area_or_volumn) {
        if (grid->get_area_or_volumn() != NULL) {
            tmp_int_value = 1;
            write_data_into_array(&tmp_int_value, sizeof(int), array, current_array_size, max_array_size);
            write_data_into_array(grid->get_area_or_volumn(), sizeof(double)*grid->get_grid_size(), array, current_array_size, max_array_size);
        }
        else {
            tmp_int_value = 0;
            write_data_into_array(&tmp_int_value, sizeof(int), array, current_array_size, max_array_size);
        }
    }
}


void Remap_weight_of_strategy_class::write_remap_weights_into_array(char **array, long &array_size, bool write_grid)
{
    Remap_grid_class *remap_grid_src, *remap_grid_dst;
    Remap_grid_class *leaf_grids[256];
    long grid_size, tmp_long_value;
    int num_leaf_grids, i, j, k, id, grid_num_dimensions, tmp_int_value;
    int num_remap_operator_instances;
    Remap_operator_basis *remap_operator_of_one_instance;
    Remap_weight_of_operator_instance_class *remap_weight_of_operator_instance;
    Remap_weight_sparse_matrix *remap_weights_group;
    char operator_name[256];
    char *output_array;
    long max_array_size;

    
    array_size = 0;
    max_array_size = 1024*1024;
    output_array = new char [max_array_size];

    remap_grid_src = get_data_grid_src();
    remap_grid_dst = get_data_grid_dst();

    if (write_grid) {
        write_grid_info_into_array(remap_grid_src, true, &output_array, array_size, max_array_size);
        write_grid_info_into_array(remap_grid_dst, true, &output_array, array_size, max_array_size);
    }
    for (i = 0, num_remap_operator_instances = 0; i < remap_weights_of_operators.size(); i ++)
        num_remap_operator_instances += remap_weights_of_operators[i]->remap_weights_of_operator_instances.size();
    write_data_into_array(&num_remap_operator_instances, sizeof(int), &output_array, array_size, max_array_size);
    for (k = 0; k < remap_weights_of_operators.size(); k ++)
        for (i = 0; i < remap_weights_of_operators[k]->remap_weights_of_operator_instances.size(); i ++) {
            remap_weight_of_operator_instance = remap_weights_of_operators[k]->remap_weights_of_operator_instances[i];
            tmp_long_value = remap_weight_of_operator_instance->get_remap_begin_iter();
            write_data_into_array(&tmp_long_value, sizeof(long), &output_array, array_size, max_array_size);
            tmp_long_value = remap_weight_of_operator_instance->get_remap_end_iter();
            write_data_into_array(&tmp_long_value, sizeof(long), &output_array, array_size, max_array_size);
            remap_operator_of_one_instance = remap_weights_of_operators[k]->remap_weights_of_operator_instances[i]->duplicated_remap_operator;
            EXECUTION_REPORT(REPORT_ERROR, -1, remap_operator_of_one_instance != NULL, "C-Coupler software error1 in write_remap_weights_into_array");
            memset(operator_name, 0, 256);
            if (write_grid) {
                strcpy(operator_name, remap_operator_of_one_instance->get_operator_name());
                write_data_into_array(operator_name, sizeof(char)*256, &output_array, array_size, max_array_size);
                write_grid_info_into_array(remap_weight_of_operator_instance->get_field_data_grid_src(), false, &output_array, array_size, max_array_size);
                write_grid_info_into_array(remap_weight_of_operator_instance->get_field_data_grid_dst(), false, &output_array, array_size, max_array_size);
                write_grid_info_into_array(remap_operator_of_one_instance->get_src_grid(), false, &output_array, array_size, max_array_size);
                write_grid_info_into_array(remap_operator_of_one_instance->get_dst_grid(), false, &output_array, array_size, max_array_size);
            }
            else {
                strcpy(operator_name, remap_operator_of_one_instance->get_object_name());
                write_data_into_array(operator_name, sizeof(char)*256, &output_array, array_size, max_array_size);
            }
            tmp_int_value = remap_operator_of_one_instance->get_num_remap_weights_groups();
            write_data_into_array(&tmp_int_value, sizeof(int), &output_array, array_size, max_array_size);
            for (j = 0; j < remap_operator_of_one_instance->get_num_remap_weights_groups(); j ++) {
                remap_weights_group = remap_operator_of_one_instance->get_remap_weights_group(j);
                tmp_long_value = remap_weights_group->get_num_weights();
                write_data_into_array(&tmp_long_value, sizeof(long), &output_array, array_size, max_array_size);
                write_data_into_array(remap_weights_group->get_indexes_src_grid(), sizeof(long)*tmp_long_value, &output_array, array_size, max_array_size);
                write_data_into_array(remap_weights_group->get_indexes_dst_grid(), sizeof(long)*tmp_long_value, &output_array, array_size, max_array_size);
                write_data_into_array(remap_weights_group->get_weight_values(), sizeof(double)*tmp_long_value, &output_array, array_size, max_array_size);
                tmp_long_value = remap_weights_group->get_num_remaped_dst_cells_indexes();
                write_data_into_array(&tmp_long_value, sizeof(long), &output_array, array_size, max_array_size);
                write_data_into_array(remap_weights_group->get_remaped_dst_cells_indexes(), sizeof(long)*tmp_long_value, &output_array, array_size, max_array_size);
            }
        }

    *array = output_array;
}


void Remap_weight_of_strategy_class::read_data_from_array(void *data, int data_size, const char *input_array, FILE *fp_binary, long &current_array_pos, long array_size, bool read_weight_values)
{
    EXECUTION_REPORT(REPORT_ERROR, -1, current_array_pos+data_size <= array_size, "the access of array is out-of-bound when reading for remapping weights %s", object_name);

    if (read_weight_values) {
        if (input_array != NULL) {
            for (long i = 0; i < data_size; i ++)
                ((char*)data)[i] = input_array[current_array_pos+i];
        }
        else fread((char*)data, 1, data_size, fp_binary);
    }
    else if (fp_binary != NULL)
        fseek(fp_binary, data_size, SEEK_CUR);
    current_array_pos += data_size;
}


void Remap_weight_of_strategy_class::read_grid_info_from_array(Remap_grid_class *grid, bool consider_area_or_volumn, const char *input_array, FILE *fp_binary, long &current_array_pos, long array_size)
{
    long grid_size;
    int grid_num_dimensions, i, gid, rid, num_leaf_grids, tmp_int_value;
    Remap_grid_class *leaf_grids[256];
    double *area_or_volumn;
    

    read_data_from_array(&grid_size, sizeof(long), input_array, fp_binary, current_array_pos, array_size, true);
    EXECUTION_REPORT(REPORT_ERROR, -1, grid_size == grid->get_grid_size(), "the grid size of %s does not match the binary file\n", grid->get_grid_name());
    read_data_from_array(&grid_num_dimensions, sizeof(int), input_array, fp_binary, current_array_pos, array_size, true);
    EXECUTION_REPORT(REPORT_ERROR, -1, grid_num_dimensions == grid->get_num_dimensions(), "the number of dimensions of grid %s does not match the binary file\n", grid->get_grid_name());
    grid->get_leaf_grids(&num_leaf_grids, leaf_grids, grid);
    for (i = 0; i < num_leaf_grids; i ++) {
        if (words_are_the_same(leaf_grids[i]->get_coord_label(), COORD_LABEL_LON))
            gid = 1;
        else if (words_are_the_same(leaf_grids[i]->get_coord_label(), COORD_LABEL_LAT))
            gid = 2;
        else if (words_are_the_same(leaf_grids[i]->get_coord_label(), COORD_LABEL_LEV))
            gid = 3;
        else if (words_are_the_same(leaf_grids[i]->get_coord_label(), COORD_LABEL_TIME))
            gid = 4;
        read_data_from_array(&rid, sizeof(int), input_array, fp_binary, current_array_pos, array_size, true);
        EXECUTION_REPORT(REPORT_ERROR, -1, gid == rid, "the arrange of coordinate systems of grid %s does not match the binary file\n", grid->get_grid_name());
    }

    if (consider_area_or_volumn) {
        read_data_from_array(&tmp_int_value, sizeof(int), input_array, fp_binary, current_array_pos, array_size, true);
        if (tmp_int_value == 1) {
            EXECUTION_REPORT(REPORT_ERROR, -1, grid->get_area_or_volumn() != NULL, "the area or volumn of grid %s does not match the binary file\n", grid->get_grid_name());
            area_or_volumn = new double [grid->get_grid_size()];
            read_data_from_array(area_or_volumn, sizeof(double)*grid->get_grid_size(), input_array, fp_binary, current_array_pos, array_size, true);
            for (long i = 0; i < grid->get_grid_size(); i ++)
                EXECUTION_REPORT(REPORT_ERROR, -1, grid->get_area_or_volumn()[i] == area_or_volumn[i], "the area or volumn of grid %s does not match the binary file\n", grid->get_grid_name());
            delete [] area_or_volumn;            
        }
        else EXECUTION_REPORT(REPORT_ERROR, -1, grid->get_area_or_volumn() == NULL, "the area or volumn of grid %s does not match the binary file\n", grid->get_grid_name());
    }
}


void Remap_weight_of_strategy_class::read_remap_operator_instance_from_array(Remap_grid_class *field_data_grid_src, Remap_grid_class *field_data_grid_dst,
                                                              Remap_grid_class *operator_grid_src, Remap_grid_class *operator_grid_dst,
                                                              Remap_operator_basis *remap_operator, long remap_begin_iter, long remap_end_iter,
                                                              const char *input_array, FILE *fp_binary, long &current_array_pos, long array_size,
                                                              bool read_weight_values)    
{
    Remap_operator_basis *duplicated_remap_operator;
    Remap_weight_of_operator_instance_class *remap_operator_instance;
    int i, num_remap_weights_groups;
    long num_weights, num_remaped_dst_cells_indexes, *indexes_src_grid, *indexes_dst_grid, *remaped_dst_cells_indexes;
    Remap_weight_sparse_matrix *weight_sparse_matrix;
    double *weight_values;


    if (read_weight_values)
        duplicated_remap_operator = remap_operator->duplicate_remap_operator(false);
    else duplicated_remap_operator = NULL;
    read_data_from_array(&num_remap_weights_groups, sizeof(int), input_array, fp_binary, current_array_pos, array_size, true);
    for (i = 0; i < num_remap_weights_groups; i ++) {
        read_data_from_array(&num_weights, sizeof(long), input_array, fp_binary, current_array_pos, array_size, true);
        if (read_weight_values) {
            indexes_src_grid = new long [num_weights];
            indexes_dst_grid = new long [num_weights];
            weight_values = new double [num_weights];
        }
        read_data_from_array(indexes_src_grid, sizeof(long)*num_weights, input_array, fp_binary, current_array_pos, array_size, read_weight_values);
        read_data_from_array(indexes_dst_grid, sizeof(long)*num_weights, input_array, fp_binary, current_array_pos, array_size, read_weight_values);
        read_data_from_array(weight_values, sizeof(double)*num_weights, input_array, fp_binary, current_array_pos, array_size, read_weight_values);
        read_data_from_array(&num_remaped_dst_cells_indexes, sizeof(long), input_array, fp_binary, current_array_pos, array_size, true);
        if (read_weight_values)
            remaped_dst_cells_indexes = new long [num_remaped_dst_cells_indexes];
        read_data_from_array(remaped_dst_cells_indexes, sizeof(long)*num_remaped_dst_cells_indexes, input_array, fp_binary, current_array_pos, array_size, read_weight_values);
        if (read_weight_values) {
            weight_sparse_matrix = new Remap_weight_sparse_matrix(remap_operator, num_weights, indexes_src_grid, indexes_dst_grid, weight_values, num_remaped_dst_cells_indexes, remaped_dst_cells_indexes);
            duplicated_remap_operator->add_weight_sparse_matrix(weight_sparse_matrix);
        }
    }
    
    remap_operator_instance = new Remap_weight_of_operator_instance_class(field_data_grid_src, field_data_grid_dst, remap_begin_iter, remap_operator, duplicated_remap_operator);
    remap_operator_instance->remap_end_iter = remap_end_iter;
    add_remap_weight_of_operator_instance(remap_operator_instance, field_data_grid_src, field_data_grid_dst, remap_operator, operator_grid_src, operator_grid_dst);
}


void Remap_weight_of_strategy_class::read_remap_weights_from_array(const char *input_array, FILE *fp_binary, long array_size, bool read_grid, Remap_grid_class **remap_related_decomp_grids, bool read_weight_values)
{
    Remap_grid_class *field_grid_src, *field_grid_dst, *current_field_grid_src, *current_field_grid_dst;
    Remap_grid_class *operator_grid_src, *operator_grid_dst;
    int i, j, k, m, num_remap_operator_instances, num_remap_operator, num_leaf_grids_all, num_leaf_grids_remap_operator;
    Remap_operator_basis *remap_operator;
    int coord_system_ids[256], tmp_grid_num_dimensions, num_sized_sub_grids;
    long tmp_grid_size, current_remap_iter, remap_end_iter;
    char operator_name[256], last_operator_name[256], tmp_grid_name[256];
    long current_array_pos = 0;
    int field_data_grids_iter = 0;
    

    if (read_grid) {
        field_grid_src = get_data_grid_src();
        field_grid_dst = get_data_grid_dst();
        m = 1;
        EXECUTION_REPORT(REPORT_ERROR, -1, remap_related_decomp_grids == NULL, "software error in read_remap_weights_from_array");
        read_grid_info_from_array(field_grid_src, true, input_array, fp_binary, current_array_pos, array_size);
        read_grid_info_from_array(field_grid_dst, true, input_array, fp_binary, current_array_pos, array_size);
    }
    else {
        EXECUTION_REPORT(REPORT_ERROR, -1, remap_related_decomp_grids != NULL, "software error in read_remap_weights_from_array");
        field_grid_src = remap_related_decomp_grids[field_data_grids_iter++];
        field_grid_dst = remap_related_decomp_grids[field_data_grids_iter++];
    }

    read_data_from_array(&num_remap_operator_instances, sizeof(int), input_array, fp_binary, current_array_pos, array_size, true);
    num_remap_operator = 0;
    current_field_grid_src = field_grid_src;
    for (i = 0; i < num_remap_operator_instances; i ++) {
        read_data_from_array(&current_remap_iter, sizeof(long), input_array, fp_binary, current_array_pos, array_size, true);
        read_data_from_array(&remap_end_iter, sizeof(long), input_array, fp_binary, current_array_pos, array_size, true);
        read_data_from_array(operator_name, sizeof(char)*256, input_array, fp_binary, current_array_pos, array_size, true);
        if (read_grid) {
            if (i == 0)
                strcpy(last_operator_name, operator_name);
            if (!words_are_the_same(last_operator_name, operator_name)) {
                num_remap_operator ++;
                strcpy(last_operator_name, operator_name);
            }
            remap_operator = remap_strategy->get_remap_operator(num_remap_operator);
            EXECUTION_REPORT(REPORT_ERROR, -1, words_are_the_same(operator_name, remap_operator->get_operator_name()),
                             "the remap operator %s does not match the binary file, which should be %s\n", 
                             operator_name, remap_operator->get_operator_name());
            read_data_from_array(&tmp_grid_size, sizeof(long), input_array, fp_binary, current_array_pos, array_size, true);
            read_data_from_array(&tmp_grid_num_dimensions, sizeof(int), input_array, fp_binary, current_array_pos, array_size, true);
            EXECUTION_REPORT(REPORT_ERROR, -1, tmp_grid_num_dimensions == field_grid_src->get_num_dimensions(), "remap software error2 in read_remap_weights_from_array binary\n");
            for (j = 0; j < tmp_grid_num_dimensions; j ++)
                read_data_from_array(&coord_system_ids[j], sizeof(int), input_array, fp_binary, current_array_pos, array_size, true);
            num_sized_sub_grids = 0;
            current_field_grid_src = field_data_grids_in_remapping_process[1+num_remap_operator*2];
            current_field_grid_dst = field_data_grids_in_remapping_process[2+num_remap_operator*2];
            read_grid_info_from_array(current_field_grid_dst, false, input_array, fp_binary, current_array_pos, array_size);
            read_grid_info_from_array(remap_operator->get_src_grid(), false, input_array, fp_binary, current_array_pos, array_size);
            read_grid_info_from_array(remap_operator->get_dst_grid(), false, input_array, fp_binary, current_array_pos, array_size);
            operator_grid_src = remap_operator->get_src_grid();
            operator_grid_dst = remap_operator->get_dst_grid();

        }
        else {
            remap_operator = remap_operator_manager->search_remap_operator(operator_name);
            EXECUTION_REPORT(REPORT_ERROR, -1, remap_operator != NULL, "software error when searching remap operator %s in read_remap_weights_from_array", operator_name);
            current_field_grid_src = remap_related_decomp_grids[field_data_grids_iter+0];
            current_field_grid_dst = remap_related_decomp_grids[field_data_grids_iter+1];
            operator_grid_src = remap_related_decomp_grids[field_data_grids_iter+2];
            operator_grid_dst = remap_related_decomp_grids[field_data_grids_iter+3];
            field_data_grids_iter += 4;
        }
        read_remap_operator_instance_from_array(current_field_grid_src, current_field_grid_dst, operator_grid_src, operator_grid_dst, remap_operator, current_remap_iter, remap_end_iter, input_array, fp_binary, current_array_pos, array_size, read_weight_values);
    }
    EXECUTION_REPORT(REPORT_ERROR, -1, current_field_grid_dst->is_similar_grid_with(field_grid_dst), "remap software error4 in read_remap_weights_from_array\n");
    EXECUTION_REPORT(REPORT_ERROR, -1, current_array_pos == array_size, "the input array does not match the remapping weights %s when reading", object_name);
}


void Remap_weight_of_strategy_class::check_remap_weights_format()
{
    Remap_grid_class *grid_src;
    int i, j, k;
    bool have_sphere_grid_remapping = false;
    

    for (k = 0; k < remap_weights_of_operators.size(); k ++)
        for (i = 0; i < remap_weights_of_operators[k]->remap_weights_of_operator_instances.size(); i ++) {
            if (remap_weights_of_operators[k]->remap_weights_of_operator_instances[i]->duplicated_remap_operator == NULL)
                continue;
            grid_src = remap_weights_of_operators[k]->remap_weights_of_operator_instances[i]->duplicated_remap_operator->get_src_grid();
            j = 0;
            if (grid_src->has_grid_coord_label(COORD_LABEL_LON))
                j ++;
            if (grid_src->has_grid_coord_label(COORD_LABEL_LAT))
                j ++;
            EXECUTION_REPORT(REPORT_ERROR, -1, j == 0 || j == 2, "the remap operator %s for coupling must remap on both longitude and latitude\n", 
                             remap_weights_of_operators[k]->remap_weights_of_operator_instances[i]->duplicated_remap_operator->get_operator_name());
            if (grid_src->has_grid_coord_label(COORD_LABEL_LON)) {
                EXECUTION_REPORT(REPORT_ERROR, -1, !have_sphere_grid_remapping, "the remap weights %s must have only one remap operator remapping on only one grid\n", 
                                 remap_weights_of_operators[k]->remap_weights_of_operator_instances[i]->duplicated_remap_operator->get_object_name());
                have_sphere_grid_remapping = true;
            }
        }
}


Remap_operator_basis *Remap_weight_of_strategy_class::get_unique_remap_operator_of_weights() 
{ 
    EXECUTION_REPORT(REPORT_ERROR, -1, remap_weights_of_operators.size() > 0 && remap_weights_of_operators[0]->remap_weights_of_operator_instances.size() > 0 &&
                     remap_weights_of_operators[0]->remap_weights_of_operator_instances[0]->duplicated_remap_operator != NULL, 
                     "C-Coupler error in get_unique_remap_operator_of_weights");

    if (remap_weights_of_operators.size() == 1 && remap_weights_of_operators[0]->remap_weights_of_operator_instances.size() == 1)
        return remap_weights_of_operators[0]->remap_weights_of_operator_instances[0]->duplicated_remap_operator;
    else return NULL;
}


void Remap_weight_of_strategy_class::add_remap_weight_of_operators_to_manager(bool are_parallel_remap_weights)
{
    public_remap_weights_of_operators = true;
    for (int i = 0; i < remap_weights_of_operators.size(); i ++)
        if (are_parallel_remap_weights)
            parallel_remap_weight_of_operator_manager->add_remap_weights_of_operator(remap_weights_of_operators[i]);
        else sequential_remap_weight_of_operator_manager->add_remap_weights_of_operator(remap_weights_of_operators[i]);
}


Remap_grid_class *Remap_weight_of_strategy_class::get_field_data_grid_in_remapping_process(int i) 
{ 
    EXECUTION_REPORT(REPORT_ERROR, -1, i < num_field_data_grids_in_remapping_process, "C-Coupler error in get_field_data_grid_in_remapping_process of Remap_weight_of_strategy_class");
    return field_data_grids_in_remapping_process[i]; 
}


Remap_grid_data_class *Remap_weight_of_strategy_class::get_runtime_mask_field_in_remapping_process(int i) 
{ 
    EXECUTION_REPORT(REPORT_ERROR, -1, i < num_field_data_grids_in_remapping_process, "C-Coupler error in get_runtime_mask_field_in_remapping_process of Remap_weight_of_strategy_class");
    return runtime_mask_fields_in_remapping_process[i]; 
}


void Remap_weight_of_strategy_class::renew_object_name(const char*new_object_name)
{
    if (words_are_the_same(object_name, new_object_name))
        return;
    
    EXECUTION_REPORT(REPORT_ERROR, -1, strncmp(object_name, "TEMP_WEIGHT", strlen("TEMP_WEIGHT")) == 0, "Remap weights %s is the same as %s. Please do not calculate the same remap weights more than once", object_name, new_object_name);
}


Remap_weight_of_operator_class *Remap_weight_of_strategy_class::get_dynamic_V1D_remap_weight_of_operator()
{
    for (int i = 0; i < remap_weights_of_operators.size(); i ++)
        if (remap_weights_of_operators[i]->operator_grid_src->has_grid_coord_label(COORD_LABEL_LEV) && (remap_weights_of_operators[i]->field_data_grid_src->is_sigma_grid() || remap_weights_of_operators[i]->field_data_grid_src->does_use_V3D_level_coord() || remap_weights_of_operators[i]->field_data_grid_dst->is_sigma_grid() || remap_weights_of_operators[i]->field_data_grid_dst->does_use_V3D_level_coord()))
            return remap_weights_of_operators[i];
    
    return NULL;
}


void Remap_weight_of_strategy_class::write_overall_H2D_remapping_weights(int comp_id)
{
	for (int i = 0; i < remap_weights_of_operators.size(); i ++)
		if (remap_weights_of_operators[i]->operator_grid_src->get_is_sphere_grid())
			remap_weights_of_operators[i]->write_overall_remapping_weights(comp_id);
}


