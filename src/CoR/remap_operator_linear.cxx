/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#include "cor_global_data.h"
#include "remap_operator_linear.h"
#include <string.h>
#include <math.h>


void Remap_operator_linear::set_parameter(const char *parameter_name, const char *parameter_value)
{
    EXECUTION_REPORT(REPORT_ERROR, -1, enable_to_set_parameters, 
                 "the parameter of remap operator object \"%s\" must be set before using it to build remap strategy\n",
                 object_name);

    set_common_parameter(parameter_name, parameter_value);
}


int Remap_operator_linear::check_parameter(const char *parameter_name, const char *parameter_value, char *error_string)
{
    return check_common_parameter(parameter_name, parameter_value, error_string);
}


void Remap_operator_linear::compute_remap_weights_of_one_dst_cell(long dst_cell_index)
{
}


void Remap_operator_linear::calculate_remap_weights()
{
    int i;
    double remap_weight_values[2], coord_differences[2];
    long weight_src_indexes[2];
    long temp_long_value = 0;
    double temp_double_value = 0.0;

    
    clear_remap_weight_info_in_sparse_matrix();

    allocate_1D_remap_operator_common_arrays_space();
    allocate_local_arrays();

    calculate_dst_src_mapping_info();

    if (array_size_src == 0)
        return;

    EXECUTION_REPORT(REPORT_ERROR, -1, array_size_src > 1, "Less than three source cells for linear interpolation are not enough");

    for (i = 0; i < dst_grid->get_grid_size(); i ++) {
        if (src_cell_index_left[i] == -1 || src_cell_index_right[i] == -1)
            continue;
        weight_src_indexes[0] = src_cell_index_left[i];
        weight_src_indexes[1] = src_cell_index_right[i];
        if (coord_values_src[src_cell_index_left[i]] == coord_values_src[src_cell_index_right[i]]) {
            remap_weight_values[0] = 0.5;
            remap_weight_values[1] = 0.5;
        }
        else if ((coord_values_dst[i] >= coord_values_src[src_cell_index_left[i]]) == (coord_values_dst[i] <= coord_values_src[src_cell_index_right[i]])) {
            remap_weight_values[1] = (coord_values_dst[i]-coord_values_src[src_cell_index_left[i]]) / (coord_values_src[src_cell_index_right[i]]-coord_values_src[src_cell_index_left[i]]);
            remap_weight_values[0] = 1 - remap_weight_values[1];
        }
        else {
            coord_differences[0] = coord_values_dst[i] - coord_values_src[src_cell_index_left[i]];
            coord_differences[1] = coord_values_dst[i] - coord_values_src[src_cell_index_right[i]];
            if (fabs(coord_differences[0]) > fabs(coord_differences[1])) {
                remap_weight_values[0] = -fabs(coord_differences[1])/fabs(coord_differences[0]-coord_differences[1]);
                remap_weight_values[1] = fabs(coord_differences[0])/fabs(coord_differences[0]-coord_differences[1]);
            }
            else {
                remap_weight_values[0] = fabs(coord_differences[1])/fabs(coord_differences[0]-coord_differences[1]);
                remap_weight_values[1] = -fabs(coord_differences[0])/fabs(coord_differences[0]-coord_differences[1]);             
            }
        }
        add_remap_weights_to_sparse_matrix(weight_src_indexes, i, remap_weight_values, 2, 1, true);        
    }
    
    if (remap_weights_groups[1]->get_num_weights() == 0) {
        EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Encounter an empty linear remapping operator instance"); 
        return;
    }

    temp_long_value = 0; 
    temp_double_value = 0.0;
    for (i = 0; i < array_size_src; i ++)
        add_remap_weights_to_sparse_matrix(&temp_long_value, useful_src_cells_global_index[i], &temp_double_value, 1, 0, false);
}


Remap_operator_linear::Remap_operator_linear(const char *object_name, int num_remap_grids, Remap_grid_class **remap_grids)
                                       : Remap_operator_1D_basis(object_name, REMAP_OPERATOR_NAME_LINEAR, num_remap_grids, remap_grids)
{
    remap_weights_groups.push_back(new Remap_weight_sparse_matrix(this));
    remap_weights_groups.push_back(new Remap_weight_sparse_matrix(this));
    
}


void Remap_operator_linear::allocate_local_arrays()
{
    temp_decomp_map_src = (long*) (common_buffer_for_1D_remap_operator + 3*(src_grid->get_grid_size()+2));
}


Remap_operator_linear::~Remap_operator_linear()
{
}


void Remap_operator_linear::do_remap_values_caculation(double *data_values_src, double *data_values_dst, int dst_array_size)
{
    int i;
    long temp_long_value1, temp_long_value2;
    double temp_double_value, base_value;
    double eps = 1.0e-8;
    
    
    array_size_src = remap_weights_groups[0]->get_num_weights();
    if (array_size_src == 0)
        return;

    allocate_1D_remap_operator_common_arrays_space();
    allocate_local_arrays();

    for (i = 0; i < array_size_src; i ++) {
        remap_weights_groups[0]->get_weight(&temp_long_value1, &temp_long_value2, &temp_double_value, i);
        useful_src_cells_global_index[i] = temp_long_value2;
    }

    preprocess_field_value(data_values_src);

    remap_weights_groups[1]->remap_values(packed_data_values_src, data_values_dst, dst_array_size);

    postprocess_field_value(data_values_dst);
}


void Remap_operator_linear::do_src_decomp_caculation(long *decomp_map_src, const long *decomp_map_dst)
{
    EXECUTION_REPORT(REPORT_ERROR, -1, false, "Software error in Remap_operator_linear::do_src_decomp_caculation: 1-D remapping algorithm should not be used to calculate src decomp");
}


Remap_operator_basis *Remap_operator_linear::duplicate_remap_operator(bool fully_copy)
{
    Remap_operator_basis *duplicated_remap_operator = new Remap_operator_linear();
    copy_remap_operator_basic_data(duplicated_remap_operator, fully_copy);
    ((Remap_operator_linear *) duplicated_remap_operator)->initialize_1D_remap_operator();
    ((Remap_operator_linear *) duplicated_remap_operator)->copy_1D_remap_operator_info(this);

    return duplicated_remap_operator;
}


Remap_operator_basis *Remap_operator_linear::generate_parallel_remap_operator(Remap_grid_class **decomp_original_grids, int **global_cells_local_indexes_in_decomps)
{
    Remap_operator_basis *parallel_remap_operator = this->duplicate_remap_operator(false);
    this->generate_parallel_remap_weights(parallel_remap_operator, decomp_original_grids, global_cells_local_indexes_in_decomps);
    return parallel_remap_operator;
}

