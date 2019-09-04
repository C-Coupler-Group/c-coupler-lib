/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#include "global_data.h"
#include "runtime_remap_function.h"
#include "cor_global_data.h"


Runtime_remap_function::Runtime_remap_function(Remap_grid_class *interchanged_grid_src,
                                          Remap_grid_class *interchanged_grid_dst,
                                          Remap_grid_class *remap_operator_runtime_grid_src,
                                          Remap_grid_class *remap_operator_runtime_grid_dst,
                                          Remap_operator_basis *runtime_remap_operator,
                                          Remap_grid_data_class *remap_field_data_src,
                                          Remap_grid_data_class *remap_field_data_dst,
                                          Remap_weight_of_strategy_class *remap_weight_of_strategy,
                                          const char *H2D_remapping_wgt_file)
{
    int num_sized_grids_of_remapping_src, num_sized_grids_of_remapping_dst, num_leaf_grids;
    Remap_grid_class *sized_grids_of_remapping_src[256], *sized_grids_of_remapping_dst[256], *leaf_grids[256];
    int num_sized_grids_of_interchanged_grid_src, num_sized_grids_of_interchanged_grid_dst;
    Remap_grid_class *sized_grids_of_interchanged_grid_src[256], *sized_grids_of_interchanged_grid_dst[256], *super_grid;
    Remap_grid_data_class *partial_redundant_mark_field;
    long i, j;


    /* Set the member variables */
    this->interchanged_grid_src = interchanged_grid_src;
    this->interchanged_grid_dst = interchanged_grid_dst;
    this->remap_operator_runtime_grid_src = remap_operator_runtime_grid_src;
    this->remap_operator_runtime_grid_dst = remap_operator_runtime_grid_dst;
    this->remap_field_data_src = remap_field_data_src;
    this->remap_field_data_dst = remap_field_data_dst;
    this->runtime_remap_operator = runtime_remap_operator;
    this->num_remapping_times = interchanged_grid_src->get_grid_size()/remap_operator_runtime_grid_src->get_grid_size();
    this->remap_weight_of_strategy = remap_weight_of_strategy;
    this->last_remapping_time_iter = -1;
    this->last_remap_weight_of_operator_instance = NULL;

    /* Check the remap software and then set num_sized_grids_of_interchanged_grid,
         sized_grids_of_interchanged_grid, index_size_array, etc */
    interchanged_grid_src->get_sized_sub_grids(&num_sized_grids_of_interchanged_grid_src, sized_grids_of_interchanged_grid_src);
    interchanged_grid_dst->get_sized_sub_grids(&num_sized_grids_of_interchanged_grid_dst, sized_grids_of_interchanged_grid_dst);
    remap_operator_runtime_grid_src->get_sized_sub_grids(&num_sized_grids_of_remapping_src, sized_grids_of_remapping_src);
    remap_operator_runtime_grid_dst->get_sized_sub_grids(&num_sized_grids_of_remapping_dst, sized_grids_of_remapping_dst);
    remap_operator_runtime_grid_src->get_leaf_grids(&this->num_leaf_grids_of_remap_operator_grid_src, 
                                                    this->leaf_grids_of_remap_operator_grid_src,
                                                    remap_operator_runtime_grid_src);
    remap_operator_runtime_grid_dst->get_leaf_grids(&this->num_leaf_grids_of_remap_operator_grid_dst, 
                                                    this->leaf_grids_of_remap_operator_grid_dst,
                                                    remap_operator_runtime_grid_dst);
    EXECUTION_REPORT(REPORT_ERROR, -1, interchanged_grid_src->get_grid_size()%remap_operator_runtime_grid_src->get_grid_size() == 0 &&
                 interchanged_grid_dst->get_grid_size()%remap_operator_runtime_grid_dst->get_grid_size() == 0, 
                 "remap software error1 in new Runtime_remap_function\n");
    for (i = 0; i < num_sized_grids_of_remapping_src; i ++)
        EXECUTION_REPORT(REPORT_ERROR, -1, sized_grids_of_remapping_src[i] == sized_grids_of_interchanged_grid_src[i], "remap software error2 in new Runtime_remap_function\n");
    for (i = 0; i < num_sized_grids_of_remapping_dst; i ++)
        EXECUTION_REPORT(REPORT_ERROR, -1, sized_grids_of_remapping_dst[i] == sized_grids_of_interchanged_grid_dst[i], "remap software error2 in new Runtime_remap_function\n");
    EXECUTION_REPORT(REPORT_ERROR, -1, num_sized_grids_of_interchanged_grid_src-num_sized_grids_of_remapping_src == num_sized_grids_of_interchanged_grid_dst-num_sized_grids_of_remapping_dst,
                 "remap software error3 in new Runtime_remap_function\n");
    for (i = num_sized_grids_of_remapping_src, j = num_sized_grids_of_remapping_dst; i < num_sized_grids_of_interchanged_grid_src; i ++)
        EXECUTION_REPORT(REPORT_ERROR, -1, sized_grids_of_interchanged_grid_src[i] == sized_grids_of_interchanged_grid_dst[j++], "remap software error5 in new Runtime_remap_function\n");
    this->num_sized_grids_of_interchanged_grid = num_sized_grids_of_interchanged_grid_src - num_sized_grids_of_remapping_src;
    for (i = 0; i < this->num_sized_grids_of_interchanged_grid; i ++) {
        this->sized_grids_of_interchanged_grid[i] = sized_grids_of_interchanged_grid_src[num_sized_grids_of_remapping_src+i];
    }

    if (runtime_remap_operator->does_require_grid_vertex_values()) {
        for (i = 0; i < num_leaf_grids_of_remap_operator_grid_src; i ++) {
            EXECUTION_REPORT(REPORT_ERROR, -1, leaf_grids_of_remap_operator_grid_src[i]->super_grid_of_setting_coord_values != NULL, "remap software error6 in new Runtime_remap_function\n");
            EXECUTION_REPORT(REPORT_ERROR, -1, leaf_grids_of_remap_operator_grid_src[i]->grid_vertex_fields.size() == 1, 
                             "remap operator %s (%s) requires users to specify vertex coordinate values in source grid %s", 
                             runtime_remap_operator->get_object_name(), runtime_remap_operator->get_operator_name(), remap_operator_runtime_grid_src->get_grid_name());
            if (leaf_grids_of_remap_operator_grid_src[i]->super_grid_of_setting_coord_values->num_dimensions > 1)
                EXECUTION_REPORT(REPORT_ERROR, -1, !leaf_grids_of_remap_operator_grid_src[i]->super_grid_of_setting_coord_values->are_vertex_values_set_in_default, 
                             "remap operator \"%s\" requires vertex values. The vertex values of \"%s\" in source grid \"%s\" are not given by users\n",
                             runtime_remap_operator->get_object_name(),
                             leaf_grids_of_remap_operator_grid_src[i]->coord_label,
                             leaf_grids_of_remap_operator_grid_src[i]->super_grid_of_setting_coord_values->grid_name);
            EXECUTION_REPORT(REPORT_ERROR, -1, leaf_grids_of_remap_operator_grid_dst[i]->super_grid_of_setting_coord_values != NULL, "remap software error8 in new Runtime_remap_function\n");
            EXECUTION_REPORT(REPORT_ERROR, -1, leaf_grids_of_remap_operator_grid_dst[i]->grid_vertex_fields.size() == 1, 
                             "remap operator %s (%s) requires users to specify vertex coordinate values in target grid %s", 
                             runtime_remap_operator->get_object_name(), runtime_remap_operator->get_operator_name(), remap_operator_runtime_grid_dst->get_grid_name());
            if (leaf_grids_of_remap_operator_grid_dst[i]->super_grid_of_setting_coord_values->num_dimensions > 1)
                EXECUTION_REPORT(REPORT_ERROR, -1, !leaf_grids_of_remap_operator_grid_dst[i]->super_grid_of_setting_coord_values->are_vertex_values_set_in_default, 
                             "remap operator \"%s\" requires vertex values. The vertex values of \"%s\" in destination grid \"%s\" are not given by users\n",
                             runtime_remap_operator->get_object_name(),
                             leaf_grids_of_remap_operator_grid_dst[i]->coord_label,
                             leaf_grids_of_remap_operator_grid_dst[i]->super_grid_of_setting_coord_values->grid_name);
        }
    }

	if (runtime_remap_operator->get_src_grid()->get_is_sphere_grid() && H2D_remapping_wgt_file != NULL) {
		runtime_remap_operator_grid_src = NULL;
		runtime_remap_operator_grid_dst = NULL;
	}
	else {
	    runtime_remap_operator_grid_src = new Remap_operator_grid(remap_operator_runtime_grid_src, runtime_remap_operator, true, false);
    	runtime_remap_operator_grid_dst = new Remap_operator_grid(remap_operator_runtime_grid_dst, runtime_remap_operator, false, false);
	}

    if (remap_field_data_dst != NULL && !remap_field_data_dst->have_data_content()) 
        remap_field_data_dst->grid_data_field->initialize_to_fill_value();

    current_mask_values_src = NULL;
    current_mask_values_dst = NULL;
    last_mask_values_src = NULL;
    last_mask_values_dst = NULL;
    if (remap_operator_runtime_grid_src->grid_mask_field != NULL) {
        current_mask_values_src = (bool*) remap_operator_runtime_grid_src->grid_mask_field->grid_data_field->data_buf;
        last_mask_values_src = new bool [remap_operator_runtime_grid_src->grid_size];
        for (i = 0; i < remap_operator_runtime_grid_src->grid_size; i ++)
            last_mask_values_src[i] = false;
    }
    if (remap_operator_runtime_grid_dst->grid_mask_field != NULL) {
        current_mask_values_dst = (bool*) remap_operator_runtime_grid_dst->grid_mask_field->grid_data_field->data_buf;
        last_mask_values_dst = new bool [remap_operator_runtime_grid_dst->grid_size];
        for (i = 0; i < remap_operator_runtime_grid_dst->grid_size; i ++)
            last_mask_values_dst[i] = false;        
    }

    /* compute overall redundant mask fields */    
    remap_field_data_redundant_mark_field_src = NULL;
    current_redundant_mark_src = NULL;
    last_redundant_mark_src = NULL;
    if (remap_operator_runtime_grid_src->redundant_cell_mark_field != NULL) {
        remap_field_data_redundant_mark_field_src = remap_operator_runtime_grid_src->redundant_cell_mark_field->duplicate_grid_data_field(interchanged_grid_src, 1, false, false);
        for (i = 0; i < remap_field_data_redundant_mark_field_src->get_coord_value_grid()->get_grid_size(); i ++)
            ((bool*)remap_field_data_redundant_mark_field_src->grid_data_field->data_buf)[i] = false;
        remap_operator_runtime_grid_src->get_leaf_grids(&num_leaf_grids, leaf_grids, remap_operator_runtime_grid_src);
        for (i = 0; i < num_leaf_grids; i ++) {
            super_grid = leaf_grids[i]->super_grid_of_setting_coord_values;
            if (super_grid == NULL)
                continue;
            if (super_grid->is_sigma_grid())
                continue;
            partial_redundant_mark_field = remap_field_data_redundant_mark_field_src->get_coord_value_grid()->expand_to_generate_full_coord_value(super_grid->redundant_cell_mark_field);
            for (j = 0; j < remap_field_data_redundant_mark_field_src->get_coord_value_grid()->get_grid_size(); j ++) 
                ((bool*)remap_field_data_redundant_mark_field_src->grid_data_field->data_buf)[j] |= ((bool*)partial_redundant_mark_field->grid_data_field->data_buf)[j];
            delete partial_redundant_mark_field;
        }
        last_redundant_mark_src = new bool [remap_operator_runtime_grid_src->grid_size];
        current_redundant_mark_src = (bool*) remap_operator_runtime_grid_src->redundant_cell_mark_field->grid_data_field->data_buf;
        for (i = 0; i < remap_operator_runtime_grid_src->grid_size; i ++) {
            last_redundant_mark_src[i] = false;
            current_redundant_mark_src[i] = ((bool*)remap_field_data_redundant_mark_field_src->grid_data_field->data_buf)[i];
        }
    }
}


Runtime_remap_function::~Runtime_remap_function()
{
	if (runtime_remap_operator_grid_src != NULL)
	    delete runtime_remap_operator_grid_src;

	if (runtime_remap_operator_grid_dst != NULL)
	    delete runtime_remap_operator_grid_dst;

    if (last_mask_values_src)
        delete [] last_mask_values_src;
    if (last_mask_values_dst)
        delete [] last_mask_values_dst;
    if (remap_field_data_redundant_mark_field_src != NULL) {
        delete [] last_redundant_mark_src;
        delete remap_field_data_redundant_mark_field_src;
    }
}


void Runtime_remap_function::calculate_static_remapping_weights(long current_remapping_time_iter, const char *H2D_remapping_wgt_file, int wgt_cal_comp_id, bool outer_mask)
/*  Calculate static remapping weights and allocate entries for dynamic remapping weights
 */
{
    long i;
    bool mask_values_have_been_changed, coord_values_have_been_changed_src, coord_values_have_been_changed_dst;
    int mask_values_status_src, mask_values_status_dst, redundant_mark_status_src;
    double *current_data_values_src, *current_data_values_dst;
    bool src_grid_changed = false, dst_grid_changed = false;


    EXECUTION_REPORT(REPORT_ERROR, -1, remap_weight_of_strategy != NULL, "Software error in Runtime_remap_function::calculate_static_remapping_weights: empty remap_weight_of_strategy");
    EXECUTION_REPORT(REPORT_ERROR, -1, current_remapping_time_iter < num_remapping_times, "Software error in Runtime_remap_function::calculate_static_remapping_weights: wrong current_remapping_time_iter. \n");
    
    if (runtime_remap_operator->get_src_grid()->has_grid_coord_label(COORD_LABEL_LEV)) {
        EXECUTION_REPORT(REPORT_ERROR, -1, runtime_remap_operator->get_src_grid()->get_num_dimensions() == 1, "Software error in Runtime_remap_function::calculate_static_remapping_weights: wrong dimension number of a remapping operator with vertical interpolation");
        if (runtime_remap_operator->get_src_grid()->get_a_leaf_grid_of_sigma_or_hybrid() || runtime_remap_operator->get_dst_grid()->get_a_leaf_grid_of_sigma_or_hybrid() || runtime_remap_operator->get_src_grid()->does_use_V3D_level_coord() || runtime_remap_operator->get_dst_grid()->does_use_V3D_level_coord()) {
			//if (outer_mask)
	            last_remap_weight_of_operator_instance = remap_weight_of_strategy->add_remap_weight_of_operator_instance(interchanged_grid_src, interchanged_grid_dst, current_remapping_time_iter, runtime_remap_operator);
            return;
        }
    }
    
    current_runtime_remap_operator_grid_src = runtime_remap_operator_grid_src;
    current_runtime_remap_operator_grid_dst = runtime_remap_operator_grid_dst;
    current_runtime_remap_operator = runtime_remap_operator;

    if (current_remapping_time_iter == 0) {
        src_grid_changed = true;
        dst_grid_changed = true;
    }

    if (remap_field_data_redundant_mark_field_src != NULL) {
        extract_runtime_field(remap_operator_runtime_grid_src, remap_field_data_redundant_mark_field_src, remap_operator_runtime_grid_src->redundant_cell_mark_field, current_remapping_time_iter);
        if (!check_mask_values_status(last_redundant_mark_src, current_redundant_mark_src, remap_operator_runtime_grid_src->grid_size))
            src_grid_changed = true;
    }
    /* Extract the runtime grid field data for runtime remapping */
    extract_runtime_field(remap_operator_runtime_grid_src, remap_operator_runtime_grid_src->original_grid_mask_field, remap_operator_runtime_grid_src->grid_mask_field, current_remapping_time_iter);
    extract_runtime_field(remap_operator_runtime_grid_dst, remap_operator_runtime_grid_dst->original_grid_mask_field, remap_operator_runtime_grid_dst->grid_mask_field, current_remapping_time_iter);
    if (!check_mask_values_status(last_mask_values_src, current_mask_values_src, remap_operator_runtime_grid_src->grid_size)) 
        src_grid_changed = true;
    if (!check_mask_values_status(last_mask_values_dst, current_mask_values_dst, remap_operator_runtime_grid_dst->grid_size))
        dst_grid_changed = true;
    
    if (src_grid_changed || dst_grid_changed) {
        if (runtime_remap_operator->get_src_grid()->get_is_sphere_grid() && H2D_remapping_wgt_file != NULL) {
            H2D_remapping_wgt_file_info *wgt_file_info = all_H2D_remapping_wgt_files_info->search_wgt_file_info(H2D_remapping_wgt_file);
			EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, wgt_cal_comp_id != -1, "Software error in Runtime_remap_function::calculate_static_remapping_weights");
			wgt_file_info->read_remapping_weights(wgt_cal_comp_id);
            EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, wgt_file_info != NULL && (wgt_file_info->get_num_wgts() == 0 || wgt_file_info->get_wgts_src_indexes() != NULL && wgt_file_info->get_wgts_dst_indexes() != NULL), "Software error in Runtime_remap_function::calculate_static_remapping_weights: empty wgt_matrix");
            Remap_weight_sparse_matrix *wgt_matrix = new Remap_weight_sparse_matrix(runtime_remap_operator, wgt_file_info->get_num_wgts(), wgt_file_info->get_wgts_src_indexes(), wgt_file_info->get_wgts_dst_indexes(), wgt_file_info->get_wgts_values(), 0, NULL);
            runtime_remap_operator->update_unique_weight_sparse_matrix(wgt_matrix);
        }
        else {
            if (src_grid_changed)
                runtime_remap_operator_grid_src->update_operator_grid_data();
            if (dst_grid_changed)
                runtime_remap_operator_grid_dst->update_operator_grid_data();
            runtime_remap_operator->calculate_remap_weights();
        }
        last_remapping_time_iter = current_remapping_time_iter;
        last_remap_weight_of_operator_instance = remap_weight_of_strategy->add_remap_weight_of_operator_instance(interchanged_grid_src, interchanged_grid_dst, current_remapping_time_iter, runtime_remap_operator);
    }
    else last_remap_weight_of_operator_instance->renew_remapping_time_end_iter(current_remapping_time_iter);
}


void Runtime_remap_function::extract_runtime_field(Remap_grid_class *remap_operator_runtime_grid, Remap_grid_data_class *global_field, Remap_grid_data_class *operator_field, long current_remapping_time_iter)
{
    if (global_field != NULL) {
        EXECUTION_REPORT(REPORT_ERROR, -1, operator_field != NULL, "remap software error8 in new extract_runtime_field\n");
        check_dimension_order_of_grid_field(global_field, remap_operator_runtime_grid);
        long grid_field_size = operator_field->grid_data_field->required_data_size;
        char *data_buf_global = (char*) global_field->grid_data_field->data_buf;
        char *data_buf_runtime = (char*) operator_field->grid_data_field->data_buf;
        memcpy(data_buf_runtime, 
               data_buf_global+current_remapping_time_iter*grid_field_size*get_data_type_size(operator_field->grid_data_field->data_type_in_application),
               grid_field_size*get_data_type_size(operator_field->grid_data_field->data_type_in_application));
    }
}


bool Runtime_remap_function::check_mask_values_status(bool *last_mask_values, bool *current_mask_values, long grid_size)
{
    bool result = true;
    long i;


    if (last_mask_values == NULL)
        return true;
    
    for (i = 0; i < grid_size; i ++) {
        if (last_mask_values[i] != current_mask_values[i]) {
            result = false;
        }
        last_mask_values[i] = current_mask_values[i];
    }

    return result;
}


void Runtime_remap_function::check_dimension_order_of_grid_field(Remap_grid_data_class *grid_data, Remap_grid_class *remap_grid)
{
    Remap_grid_class *sized_grids_of_remapping[256];
    int i, j, last_order_indx, num_sized_grids_of_remapping;


    remap_grid->get_sized_sub_grids(&num_sized_grids_of_remapping, sized_grids_of_remapping);
    for (i = 0, j = num_sized_grids_of_remapping; i < num_sized_grids_of_interchanged_grid; i ++)
        sized_grids_of_remapping[j++] = sized_grids_of_interchanged_grid[i];
    num_sized_grids_of_remapping += num_sized_grids_of_interchanged_grid;

    last_order_indx = -1;
    for (i = 0; i < grid_data->sized_grids.size(); i ++) {
        for (j = 0; j < num_sized_grids_of_remapping; j ++) {
            if (grid_data->sized_grids[i] == sized_grids_of_remapping[j]) 
                break;
        }
        EXECUTION_REPORT(REPORT_ERROR, -1, j < num_sized_grids_of_remapping, "remap software error1 in check_dimension_order_of_grid_field\n");
        EXECUTION_REPORT(REPORT_ERROR, -1, j > last_order_indx, "remap software error2 in check_dimension_order_of_grid_field\n");
        last_order_indx = j;
    }
}

