/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#include "cor_global_data.h"
#include "remap_strategy_class.h"
#include "runtime_remap_function.h"
#include "remap_operator_grid.h"
#include "remap_common_utils.h"
#include "execution_report.h"
#include <string.h>


Remap_strategy_class::Remap_strategy_class(const char *strategy_name, int num_remap_operators, Remap_operator_basis **remap_operators)
{
    int i, j, k;
    int num_temp_leaf_grids, num_all_leaf_grids;
    Remap_grid_class *temp_leaf_grids[128], *all_leaf_grids[128];
    bool src_leaf_grids_mask[128], dst_leaf_grids_mask[128];
    Remap_grid_class *remap_operator_src_grid, *remap_operator_dst_grid;
    char tmp_str[128];
    

    strcpy(this->strategy_name, strategy_name);

/*
    remap_grid_manager->get_all_leaf_remap_grids(&num_all_leaf_grids, all_leaf_grids);
    for (i = 0; i < num_all_leaf_grids; i ++) {
        src_leaf_grids_mask[i] = true;
        dst_leaf_grids_mask[i] = true;
    }

    for (i = 0; i < num_remap_operators; i ++) {
        remap_operator_src_grid = remap_operators[i]->get_src_grid();
        remap_operator_dst_grid = remap_operators[i]->get_dst_grid();
        remap_operator_src_grid->get_leaf_grids(&num_temp_leaf_grids, temp_leaf_grids, remap_operator_src_grid);
        sprintf(tmp_str, "%dth", i+1);
        for (j = 0; j < num_temp_leaf_grids; j ++) {
            for (k = 0; k < num_all_leaf_grids; k ++)
                if (all_leaf_grids[k] == temp_leaf_grids[j] && src_leaf_grids_mask[k])
                    break;
            EXECUTION_REPORT(REPORT_ERROR, -1, k < num_all_leaf_grids,
                         "the source leaf grid \"%s\" of the %s remap operator \"%s\" does not satisfy the constraint of cascading remap operators\n",
                         temp_leaf_grids[j]->get_grid_name(), tmp_str, remap_operators[i]->get_object_name());
            src_leaf_grids_mask[k] = false;
            dst_leaf_grids_mask[k] = true;
        }
        remap_operator_dst_grid->get_leaf_grids(&num_temp_leaf_grids, temp_leaf_grids, remap_operator_dst_grid);        
        for (j = 0; j < num_temp_leaf_grids; j ++) {
            for (k = 0; k < num_all_leaf_grids; k ++)
                if (all_leaf_grids[k] == temp_leaf_grids[j] && dst_leaf_grids_mask[k])
                    break;
            EXECUTION_REPORT(REPORT_ERROR, -1, k < num_all_leaf_grids, 
                         "the destination leaf grid object \"%s\" of %s remap operator object \"%s\" does not satisfy the constraint of cascading remap operators\n",
                         temp_leaf_grids[j]->get_grid_name(), tmp_str, remap_operators[i]->get_object_name());
            src_leaf_grids_mask[k] = true;
            dst_leaf_grids_mask[k] = false;
        }
    }
*/

    for (i = 0; i < num_remap_operators; i ++) {
        this->remap_operators.push_back(remap_operators[i]);
        remap_operators[i]->disable_to_set_parameters();
    }
}


bool Remap_strategy_class::match_remap_strategy(const char *strategy_name)
{
    return words_are_the_same(this->strategy_name, strategy_name);
}


void Remap_strategy_class::check_field_data_grid_center_values_for_remapping(Remap_grid_class *field_data_grid, Remap_grid_class *remap_operator_grid, bool is_remap_operator_regriding)
{
    Remap_grid_class *leaf_grids_field_data[256], *leaf_grids_remap_operator[256];
    int num_leaf_grids_field_data, num_leaf_grids_remap_operator, i;


    EXECUTION_REPORT(REPORT_ERROR, -1, remap_operator_grid->is_subset_of_grid(field_data_grid), "remap operator grid \"%s\" must be a sub grid of field data\n", remap_operator_grid->get_grid_name());

    field_data_grid->get_leaf_grids(&num_leaf_grids_field_data, leaf_grids_field_data, field_data_grid);
    remap_operator_grid->get_leaf_grids(&num_leaf_grids_remap_operator, leaf_grids_remap_operator, remap_operator_grid);

    for (i = 0; i < num_leaf_grids_field_data; i ++)
        if (leaf_grids_field_data[i]->get_super_grid_of_setting_coord_values() != NULL) {
            if (!leaf_grids_field_data[i]->get_super_grid_of_setting_coord_values()->is_sigma_grid())
                EXECUTION_REPORT(REPORT_ERROR, -1, leaf_grids_field_data[i]->get_super_grid_of_setting_coord_values()->is_subset_of_grid(field_data_grid), 
                             "The grid of setting coordinate values of 1D grid \"%s\" is \"%s\", it must be a sub grid of field data grid\n",
                             leaf_grids_field_data[i]->get_coord_label(), leaf_grids_field_data[i]->get_super_grid_of_setting_coord_values()->get_grid_name());
        }
    if (is_remap_operator_regriding)
        for (i = 0; i < num_leaf_grids_remap_operator; i ++) {
			if (leaf_grids_remap_operator[i]->does_use_V3D_level_coord()) {
			}
            else if (!(leaf_grids_remap_operator[i]->has_grid_coord_label(COORD_LABEL_LEV) && field_data_grid->is_sigma_grid())) {    
                EXECUTION_REPORT(REPORT_ERROR, -1, leaf_grids_remap_operator[i]->get_super_grid_of_setting_coord_values() != NULL,
                             "The coordinate values of \"%s\" defined in grid \"%s\" must be set for regriding\n",
                             leaf_grids_remap_operator[i]->get_coord_label(),
                             leaf_grids_remap_operator[i]->get_grid_name());
                if (!leaf_grids_remap_operator[i]->get_super_grid_of_setting_coord_values()->is_subset_of_grid(remap_operator_grid))
                    EXECUTION_REPORT(REPORT_ERROR, -1, remap_operator_grid->get_num_dimensions() == 1,
                                 "the coordinate values of %s are set in grid %s, it must be used in 1D remapping\n",
                                 leaf_grids_remap_operator[i]->get_coord_label(),
                                 leaf_grids_remap_operator[i]->get_super_grid_of_setting_coord_values()->get_grid_name());
            }
        }
}


void Remap_strategy_class::calculate_remapping_weights(Remap_weight_of_strategy_class *remap_weight_of_strategy, const char *H2D_remapping_wgt_file, int wgt_cal_wgt_id)
{
    int i, j;
    Remap_grid_class *remap_src_data_grid, *remap_dst_data_grid, *current_remap_src_data_grid;
    Remap_grid_class *current_remap_dst_data_grid, *existing_grid;
    Remap_grid_class *current_remap_src_data_grid_interchanged;
    int num_leaf_grids, num_sized_grids;
    Remap_grid_class *leaf_grids[256], *sized_grids[256];
    Remap_grid_class *runtime_remap_grid_src, *runtime_remap_grid_dst;
    Remap_grid_data_class *runtime_mask_src, *runtime_mask_dst;
    Remap_grid_class *runtime_mask_sub_grids_src[256], *runtime_mask_sub_grids_dst[256];
	bool *outer_mask;
    int num_runtime_mask_sub_grids_src, num_runtime_mask_sub_grids_dst;
    long runtime_remap_times_iter;
    double last_time, current_time;


    remap_src_data_grid = remap_weight_of_strategy->get_data_grid_src();
    remap_dst_data_grid = remap_weight_of_strategy->get_data_grid_dst();

    remap_src_data_grid->end_grid_definition_stage(NULL);
    remap_dst_data_grid->end_grid_definition_stage(NULL);

    if (remap_src_data_grid->is_sigma_grid() || remap_dst_data_grid->is_sigma_grid()) {
        for (i = 0; i < remap_operators.size(); i ++)
            if (remap_operators[i]->get_src_grid()->has_grid_coord_label(COORD_LABEL_LEV))
                break;
        EXECUTION_REPORT(REPORT_ERROR, -1, i < remap_operators.size(), "One grid in %s or %s is a sigma grid. The interpolation between them must have vertical interpolation", remap_src_data_grid->get_grid_name(), remap_dst_data_grid->get_grid_name());
    }

     j = 1;
    
    current_remap_src_data_grid = remap_src_data_grid;
    for (i = 0; i < remap_operators.size(); i ++) {
        EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "execute remap operator %s  %s %s %lx %s", remap_operators[i]->get_object_name(), remap_operators[i]->get_operator_name(), remap_operators[i]->get_dst_grid()->get_grid_name(), remap_operators[i]->get_dst_grid(), remap_dst_data_grid->get_grid_name());
        EXECUTION_REPORT(REPORT_ERROR, -1, remap_operators[i]->get_dst_grid()->is_sub_grid_of_grid(remap_dst_data_grid), "Software error: the dst grid \"%s\" of an remapping operator is not a sub grid of \"%s\"", remap_operators[i]->get_dst_grid()->get_grid_name(), remap_dst_data_grid->get_grid_name());
        current_remap_src_data_grid_interchanged = remap_weight_of_strategy->get_field_data_grid_in_remapping_process(j);
        runtime_mask_src = remap_weight_of_strategy->get_runtime_mask_field_in_remapping_process(j++);
        current_remap_dst_data_grid = remap_weight_of_strategy->get_field_data_grid_in_remapping_process(j);
        check_field_data_grid_center_values_for_remapping(current_remap_src_data_grid, get_remap_operator(i)->get_src_grid(), get_remap_operator(i)->get_is_operator_regridding());
        check_field_data_grid_center_values_for_remapping(current_remap_dst_data_grid, get_remap_operator(i)->get_dst_grid(), get_remap_operator(i)->get_is_operator_regridding());    
        runtime_mask_dst = remap_weight_of_strategy->get_runtime_mask_field_in_remapping_process(j++);
        current_remap_src_data_grid->interchange_grid_fields_for_remapping(current_remap_src_data_grid_interchanged,
                                                                           remap_operators[i]->get_src_grid(),
                                                                           runtime_mask_src);
        current_remap_dst_data_grid->interchange_grid_fields_for_remapping(current_remap_dst_data_grid,
                                                                           remap_operators[i]->get_dst_grid(),
                                                                           runtime_mask_dst);
        runtime_remap_grid_src = current_remap_src_data_grid_interchanged->generate_remap_operator_runtime_grid(remap_operators[i]->get_src_grid(), 
                                                                                                                remap_operators[i], 
                                                                                                                runtime_mask_src);
        runtime_remap_grid_dst = current_remap_dst_data_grid->generate_remap_operator_runtime_grid(remap_operators[i]->get_dst_grid(), 
                                                                                                   remap_operators[i], 
                                                                                                   runtime_mask_dst);
        current_runtime_remap_function = new Runtime_remap_function(current_remap_src_data_grid_interchanged,
                                                                    current_remap_dst_data_grid,
                                                                    runtime_remap_grid_src,
                                                                    runtime_remap_grid_dst,
                                                                    remap_operators[i],
                                                                    NULL,
                                                                    NULL,
                                                                    remap_weight_of_strategy,
                                                                    H2D_remapping_wgt_file);
        if (execution_phase_number == 1) {
			outer_mask = NULL;
			if (runtime_remap_grid_src->get_grid_mask_field() == NULL && (runtime_mask_src != NULL || runtime_mask_dst != NULL)) {
				EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, runtime_mask_src != NULL && runtime_mask_dst != NULL, "Software error in Remap_strategy_class::calculate_remapping_weights");
				EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, runtime_mask_src->get_coord_value_grid() == runtime_mask_dst->get_coord_value_grid(), "Software error in Remap_strategy_class::calculate_remapping_weights");				
				EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, runtime_mask_src->get_grid_data_field()->required_data_size == current_remap_src_data_grid_interchanged->get_grid_size()/runtime_remap_grid_src->get_grid_size(), "Software error in Remap_strategy_class::calculate_remapping_weights");
				outer_mask = (bool*) runtime_mask_src->get_grid_data_field()->data_buf;
			}
            for (runtime_remap_times_iter = 0; runtime_remap_times_iter < current_remap_src_data_grid_interchanged->get_grid_size()/runtime_remap_grid_src->get_grid_size(); runtime_remap_times_iter ++) {
                current_runtime_remap_function->calculate_static_remapping_weights(runtime_remap_times_iter, H2D_remapping_wgt_file, wgt_cal_wgt_id, outer_mask == NULL? true:outer_mask[runtime_remap_times_iter]);
            }
        }
        if (remap_operators[i]->get_src_grid()->has_grid_coord_label(COORD_LABEL_LEV))
            if (remap_operators[i]->get_src_grid()->get_a_leaf_grid_of_sigma_or_hybrid() || remap_operators[i]->get_dst_grid()->get_a_leaf_grid_of_sigma_or_hybrid())
                remap_weight_of_strategy->mark_empty_remap_weight();

        delete runtime_remap_grid_src;
        delete runtime_remap_grid_dst;
        delete current_runtime_remap_function;
        current_remap_src_data_grid = current_remap_dst_data_grid;

        if (runtime_mask_src != NULL)
            delete runtime_mask_src;

        if (runtime_mask_dst != NULL)
            delete runtime_mask_dst;
    }

    EXECUTION_REPORT(REPORT_ERROR, -1, j+1 == remap_weight_of_strategy->get_num_field_data_grids_in_remapping_process(), "C-Coupler error in calculate_remapping_weights\n");
}


void Remap_strategy_class::remap_fields(const char *field_data_name_src, const char *field_data_name_dst)
{
    int i;
    Remap_grid_class *remap_src_data_grid, *remap_dst_data_grid;
    Remap_grid_data_class *field_data_src, *field_data_dst;
    Remap_weight_of_strategy_class *remap_weight_of_strategy;


    field_data_src = remap_field_data_manager->search_remap_field_data(field_data_name_src);
    field_data_dst = remap_field_data_manager->search_remap_field_data(field_data_name_dst);
    field_data_src->transfer_field_attributes_to_another(field_data_dst);
    EXECUTION_REPORT(REPORT_ERROR, -1, field_data_src != NULL && field_data_dst != NULL, "remap software error1 in Remap_strategy_class::remap_fields\n");
    EXECUTION_REPORT(REPORT_WARNING, -1, field_data_src->have_data_content(), 
                     "source field data \"%s\" does not have essential data; before the current remapping calculation, its data should have been read from IO or calculated by remapping\n",
                     field_data_name_src);
    EXECUTION_REPORT(REPORT_WARNING, -1, !field_data_dst->have_data_content(), 
                     "destination field data \"%s\" should not have essentia data, which have been read from IO or calculated by remapping\n",
                     field_data_name_dst);
    remap_src_data_grid = field_data_src->get_coord_value_grid();
    remap_dst_data_grid = field_data_dst->get_coord_value_grid();
    EXECUTION_REPORT(REPORT_ERROR, -1, remap_src_data_grid != NULL, 
                     "\"%s\" is not a field data but a normal array data. It can not be the source data of remapping process\n",
                     field_data_src->get_grid_data_field()->field_name_in_application);
    EXECUTION_REPORT(REPORT_ERROR, -1, remap_dst_data_grid != NULL, 
                     "\"%s\" is not a field data but a normal array data. It can not be the destination data of remapping process\n",
                     field_data_dst->get_grid_data_field()->field_name_in_application);

    remap_src_data_grid->end_grid_definition_stage(NULL);
    remap_dst_data_grid->end_grid_definition_stage(NULL);

    if (remap_src_data_grid->is_sigma_grid() || remap_dst_data_grid->is_sigma_grid()) {
        for (i = 0; i < remap_operators.size(); i ++)
            if (remap_operators[i]->get_src_grid()->has_grid_coord_label(COORD_LABEL_LEV))
                break;
        EXECUTION_REPORT(REPORT_ERROR, -1, i < remap_operators.size(), "One grid in %s or %s is a sigma grid. The interpolation between them must have vertical interpolation", remap_src_data_grid->get_grid_name(), remap_dst_data_grid->get_grid_name());
    }

    remap_weight_of_strategy = remap_weights_of_strategy_manager->search_or_add_remap_weight_of_strategy(remap_src_data_grid, remap_dst_data_grid, this, NULL, NULL, NULL, false);
    remap_weight_of_strategy->do_remap(-1, field_data_src, field_data_dst);
}


