/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#include "cor_global_data.h"
#include "remap_operator_distwgt.h"
#include "remap_utils_nearest_points.h"
#include <string.h>


void Remap_operator_distwgt::set_parameter(const char *parameter_name, const char *parameter_value)
{
    EXECUTION_REPORT(REPORT_ERROR, -1, enable_to_set_parameters, 
                 "the parameter of remap operator object \"%s\" must be set before using it to build remap strategy\n",
                 object_name);

    if (words_are_the_same(parameter_name, "num_power"))
        sscanf(parameter_value, "%lf", &num_power);
    else if (words_are_the_same(parameter_name, "num_nearest_points")) {
        sscanf(parameter_value, "%d", &num_nearest_points);
        if (weigt_values_of_one_dst_cell != NULL)
            delete [] weigt_values_of_one_dst_cell;
        weigt_values_of_one_dst_cell = new double [num_nearest_points];
    }
    if (words_are_the_same(parameter_name, "enable_extrapolate")) {
        if (words_are_the_same(parameter_value, "true"))
            enable_extrapolate = true;
        else if (words_are_the_same(parameter_value, "false"))
            enable_extrapolate = false;
        else EXECUTION_REPORT(REPORT_ERROR, -1, false, "value of the parameter \"enable_extrapolate\" must be \"true\"\n");
    }
    else EXECUTION_REPORT(REPORT_ERROR, -1, false, 
                      "\"%s\" is a illegal parameter of remap operator \"%s\"\n",
                      parameter_name,
                      operator_name);
}


int Remap_operator_distwgt::check_parameter(const char *parameter_name, const char *parameter_value, char *error_string)
{
    int check_result = 0;
    if (words_are_the_same(parameter_name, "num_power")) {
        check_result = 1;
        sscanf(parameter_value, "%lf", &num_power);
        if (num_power > 0) 
            check_result = 3;
        else sprintf(error_string, "The parameter value must be larger than 0");
    }
    else if (words_are_the_same(parameter_name, "num_nearest_points")) {
        check_result = 1;
        sscanf(parameter_value, "%d", &num_nearest_points);
        if (num_nearest_points> 0) 
            check_result = 3;
        else sprintf(error_string, "The parameter value must be larger than 0");
    }
    if (words_are_the_same(parameter_name, "enable_extrapolate")) {
        check_result = 1;
        if (words_are_the_same(parameter_value, "true") || words_are_the_same(parameter_value, "false"))
            check_result = 3;
        else sprintf(error_string, "The parameter value must be \"true\" or \"false\"");
    }

    return check_result;
}


Remap_operator_distwgt::Remap_operator_distwgt()
{
    found_nearest_points_distance = NULL;
    found_nearest_points_src_indexes = NULL;
    weigt_values_of_one_dst_cell = NULL;
    enable_extrapolate = false;
}


Remap_operator_distwgt::Remap_operator_distwgt(const char *object_name, int num_remap_grids, Remap_grid_class **remap_grids)
                                       : Remap_operator_basis(object_name, 
                                                              REMAP_OPERATOR_NAME_DISTWGT, 
                                                              2, 
                                                              true, 
                                                              false, 
                                                              true,
                                                              num_remap_grids, 
                                                              remap_grids)
{
    num_nearest_points = 4;
    num_power = 1;
    enable_extrapolate = false;
    found_nearest_points_distance = new double [src_grid->get_grid_size()];
    found_nearest_points_src_indexes = new long [src_grid->get_grid_size()];
    weigt_values_of_one_dst_cell = new double [num_nearest_points];
    remap_weights_groups.push_back(new Remap_weight_sparse_matrix(this));
}


void Remap_operator_distwgt::compute_remap_weights_of_one_dst_cell(long dst_cell_index)
{
    initialize_computing_remap_weights_of_one_cell();
    compute_dist_remap_weights_of_one_dst_cell(dst_cell_index, 
                                               num_nearest_points, 
                                               num_power,
                                               &threshold_distance,
                                               found_nearest_points_distance,
                                               found_nearest_points_src_indexes,
                                               weigt_values_of_one_dst_cell,
                                               get_is_sphere_grid(),
                                               enable_extrapolate);
    finalize_computing_remap_weights_of_one_cell();    
}


void Remap_operator_distwgt::calculate_remap_weights()
{    
    threshold_distance = 1.0/6000.0;
    clear_remap_weight_info_in_sparse_matrix();
    
    for (long dst_cell_index = 0; dst_cell_index < dst_grid->get_grid_size(); dst_cell_index ++) {
        if (H2D_grid_decomp_mask != NULL && !H2D_grid_decomp_mask[dst_cell_index])
            continue;
        compute_remap_weights_of_one_dst_cell(dst_cell_index);
    }
}


Remap_operator_distwgt::~Remap_operator_distwgt()
{
    if (found_nearest_points_distance != NULL)
        delete [] found_nearest_points_distance;
    if (found_nearest_points_src_indexes != NULL)
        delete [] found_nearest_points_src_indexes;
    if (weigt_values_of_one_dst_cell != NULL)
        delete [] weigt_values_of_one_dst_cell;
}


void Remap_operator_distwgt::do_remap_values_caculation(double *data_values_src, double *data_values_dst, int dst_array_size)
{
    remap_weights_groups[0]->remap_values(data_values_src, data_values_dst, dst_array_size);
}


void Remap_operator_distwgt::do_src_decomp_caculation(long *decomp_map_src, const long *decomp_map_dst)
{
    remap_weights_groups[0]->calc_src_decomp(decomp_map_src, decomp_map_dst);
}


Remap_operator_basis *Remap_operator_distwgt::duplicate_remap_operator(bool fully_copy)
{
    Remap_operator_distwgt *duplicated_remap_operator = new Remap_operator_distwgt();
    copy_remap_operator_basic_data(duplicated_remap_operator, fully_copy);
    duplicated_remap_operator->num_power = num_power;
    duplicated_remap_operator->num_nearest_points = num_nearest_points;
    duplicated_remap_operator->threshold_distance = threshold_distance;

    return duplicated_remap_operator;
}


Remap_operator_basis *Remap_operator_distwgt::generate_parallel_remap_operator(Remap_grid_class **decomp_original_grids, int **global_cells_local_indexes_in_decomps)
{
    Remap_operator_basis *parallel_remap_operator = this->duplicate_remap_operator(false);
    this->generate_parallel_remap_weights(parallel_remap_operator, decomp_original_grids, global_cells_local_indexes_in_decomps);
    return parallel_remap_operator;
}
