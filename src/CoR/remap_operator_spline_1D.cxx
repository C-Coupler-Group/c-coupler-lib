/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file is initially finished by Mr. Yufeng Zhou,
  *  and then upgraded and merged into CoR by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#include "cor_global_data.h"
#include "remap_operator_spline_1D.h"
#include "remap_common_utils.h"
#include <string.h>
#include <math.h>
#include <assert.h>


void Remap_operator_spline_1D::solve_aperiodic_tridiagonal_system(double *a, double *b, double *c, double *f, int n)
{
    double y;


    for (int i = 1; i < n; i++) {
        y = a[i]/b[i-1];
        b[i] -= c[i-1]*y;
        f[i] -= f[i-1]*y;
    }

    f[n-1] /= b[n-1];
    for (int i = n-2; i >=0; i--)
        f[i] = (f[i]-f[i+1]*c[i])/b[i];
}


void Remap_operator_spline_1D::solve_periodic_tridiagonal_system(double *a, double *b, double *c, double *f, int n)
{
    double y;
    int i;


    temp_array_column[0] = a[0];
    temp_array_row[0] = c[n-1];
    for (i = 1; i < n-2; i++) {
        y = a[i]/b[i-1];
        b[i] -= c[i-1]*y;
        temp_array_column[i] = -temp_array_column[i-1]*y;
        f[i] -= f[i-1]*y;
        y = temp_array_row[i-1]/b[i-1];
        temp_array_row[i] = -c[i-1]*y;
        b[n-1] -= temp_array_column[i-1]*y;
        f[n-1] -= f[i-1]*y;
    }    

    y = a[n-2]/b[n-3];
    b[n-2] -= c[n-3]*y;
    c[n-2] -= temp_array_column[n-3]*y;
    f[n-2] -= f[n-3]*y;
    y = temp_array_row[n-3]/b[n-3];
    a[n-1] -= c[n-3]*y;
    b[n-1] -= temp_array_column[n-3]*y;
    f[n-1] -= f[n-3]*y;

    y = a[n-1]/b[n-2];
    b[n-1] -= c[n-2]*y;
    f[n-1] -= f[n-2]*y;
    
    // back substitution
    f[n-1] /= b[n-1];
    b[n-1] = 1.0;
    f[n-2] -= f[n-1]*c[n-2]/b[n-1];
    f[n-2] /= b[n-2];
    b[n-2] = 1.0;
    for (i = n-3; i >= 0; i--) {
        f[i] = f[i]-c[i]*f[i+1]-temp_array_column[i]*f[n-1];
        f[i] /= b[i];
        b[i] = 1.0;
    }
}


void Remap_operator_spline_1D::set_parameter(const char *parameter_name, const char *parameter_value)
{
    EXECUTION_REPORT(REPORT_ERROR, -1, enable_to_set_parameters, 
                 "the parameter of remap operator object \"%s\" must be set before using it to build remap strategy\n",
                 object_name);

    if (words_are_the_same(parameter_name, "keep_monotonicity")) {
        EXECUTION_REPORT(REPORT_ERROR, -1, !set_keep_monotonicity,
                         "The parameter \"%s\" of the 1D spline remapping operator \"%s\" has been set before. It can not been set more than once",
                         parameter_name, operator_name);
        if (words_are_the_same(parameter_value, "true")) 
            keep_monotonicity = true;
        else if (words_are_the_same(parameter_value, "false"))
            keep_monotonicity = false;
        else EXECUTION_REPORT(REPORT_ERROR, -1, false, 
                      "The value of parameter \"%s\" of the 1D spline remapping operator \"%s\" must be \"none\", \"overall\" or \"fragment\"",
                      parameter_name, operator_name);
        set_keep_monotonicity = true;
    }
    else set_common_parameter(parameter_name, parameter_value);
}


int Remap_operator_spline_1D::check_parameter(const char *parameter_name, const char *parameter_value, char *error_string)
{
    int check_result = 0;

    if (words_are_the_same(parameter_name, "keep_monotonicity")) {
        if (words_are_the_same(parameter_value, "true") || words_are_the_same(parameter_value, "false"))
            check_result = 3;
        else sprintf(error_string, "The parameter value must be \"true\" or \"false\"");
        return check_result;
    }
    else return check_common_parameter(parameter_name, parameter_value, error_string);
}


void Remap_operator_spline_1D::allocate_local_arrays()
{
    array_alpha = common_buffer_for_1D_remap_operator + 3*(src_grid->get_grid_size()+2);
    array_mu = common_buffer_for_1D_remap_operator + 4*(src_grid->get_grid_size()+2);;
    array_lambda = common_buffer_for_1D_remap_operator + 5*(src_grid->get_grid_size()+2);;
    array_h = common_buffer_for_1D_remap_operator + 6*(src_grid->get_grid_size()+2);;
    array_d = common_buffer_for_1D_remap_operator + 7*(src_grid->get_grid_size()+2);;
    temp_array_column = common_buffer_for_1D_remap_operator + 8*(src_grid->get_grid_size()+2);;
    temp_array_row = common_buffer_for_1D_remap_operator + 9*(src_grid->get_grid_size()+2);

    final_factor1 = common_buffer_for_1D_remap_operator + 12*(src_grid->get_grid_size()+2) + 3*(dst_grid->get_grid_size()+2);
    final_factor2 = common_buffer_for_1D_remap_operator + 12*(src_grid->get_grid_size()+2) + 4*(dst_grid->get_grid_size()+2);
    final_factor3 = common_buffer_for_1D_remap_operator + 12*(src_grid->get_grid_size()+2) + 5*(dst_grid->get_grid_size()+2);
    final_factor4 = common_buffer_for_1D_remap_operator + 12*(src_grid->get_grid_size()+2) + 6*(dst_grid->get_grid_size()+2);
    final_factor5 = common_buffer_for_1D_remap_operator + 12*(src_grid->get_grid_size()+2) + 7*(dst_grid->get_grid_size()+2);
    data_in_monotonicity_range = common_buffer_for_1D_remap_operator + 12*(src_grid->get_grid_size()+2) + 8*(dst_grid->get_grid_size()+2);
    dst_cell_indexes_in_monotonicity_ranges = (int*) (common_buffer_for_1D_remap_operator + 12*(src_grid->get_grid_size()+2) + 9*(dst_grid->get_grid_size()+2));
}


Remap_operator_spline_1D::Remap_operator_spline_1D(const char *object_name, int num_remap_grids, Remap_grid_class **remap_grids)
                                       : Remap_operator_1D_basis(object_name, REMAP_OPERATOR_NAME_SPLINE_1D, num_remap_grids, remap_grids)
{
    set_periodic = false;
    set_period = false;
    periodic = false;
    enable_extrapolate = false;
    set_enable_extrapolation = false;
    keep_monotonicity = false;
    set_keep_monotonicity = false;
    allocate_local_arrays();
    remap_weights_groups.push_back(new Remap_weight_sparse_matrix(this));
    remap_weights_groups.push_back(new Remap_weight_sparse_matrix(this));
    remap_weights_groups.push_back(new Remap_weight_sparse_matrix(this));
    remap_weights_groups.push_back(new Remap_weight_sparse_matrix(this));
    remap_weights_groups.push_back(new Remap_weight_sparse_matrix(this));
}


Remap_operator_spline_1D::~Remap_operator_spline_1D()
{
}


void Remap_operator_spline_1D::compute_remap_weights_of_one_dst_cell(long cell_index_dst)
{
}


void Remap_operator_spline_1D::calculate_remap_weights()
{
    int i;
    long temp_long_value = 0.0;



    allocate_1D_remap_operator_common_arrays_space();
    allocate_local_arrays();

    clear_remap_weight_info_in_sparse_matrix();
    calculate_dst_src_mapping_info();

    if (num_useful_src_cells == 0)
        return;
    
    EXECUTION_REPORT(REPORT_ERROR, -1, array_size_src > 2, "Less than three source cells for 1D spline interpolation are not enough");

    for (i = 0; i < array_size_src-1; i ++)
        array_h[i] = coord_values_src[i+1]-coord_values_src[i];

    array_mu[0] = 0.0;
    for (i = 1; i < array_size_src-1; i ++) {
        array_mu[i] = array_h[i-1]/(array_h[i-1]+array_h[i]);
        array_lambda[i] = 1.0-array_mu[i];
    }

    if (!periodic) {
        array_lambda[0] = 0.0;
        array_mu[array_size_src-1] = 0.0;
    } else {
        array_lambda[array_size_src-1] = array_h[0]/(array_h[array_size_src-2]+array_h[0]);
        array_mu[array_size_src-1] = 1.0-array_lambda[array_size_src-1];
    }

    for (i = 0; i < dst_grid->get_grid_size(); i ++) {
        final_factor1[i] = 0;
        final_factor2[i] = 0;
        final_factor3[i] = 0;
        final_factor4[i] = 0;
        final_factor5[i] = 0;
        if (src_cell_index_left[i] == -1 || src_cell_index_right[i] == -1)
            continue;
        final_factor1[i] = pow(coord_values_src[src_cell_index_right[i]]-coord_values_dst[i], 3.0)/(6.0*array_h[src_cell_index_left[i]]);
        final_factor2[i] = pow(coord_values_dst[i]-coord_values_src[src_cell_index_left[i]], 3.0)/(6.0*array_h[src_cell_index_left[i]]);
        final_factor3[i] = array_h[src_cell_index_left[i]]*array_h[src_cell_index_left[i]]/6.0;
        final_factor4[i] = (coord_values_src[src_cell_index_right[i]]-coord_values_dst[i])/array_h[src_cell_index_left[i]];
        final_factor5[i] = (coord_values_dst[i]-coord_values_src[src_cell_index_left[i]])/array_h[src_cell_index_left[i]];
    }


    for (i = 0; i < array_size_src; i ++) {
        temp_long_value = 0;
        add_remap_weights_to_sparse_matrix((long*)(&array_mu[i]), useful_src_cells_global_index[i], array_lambda+i, 1, 0, false);
        add_remap_weights_to_sparse_matrix((long*)(&coord_values_src[i]), temp_long_value, array_h+i, 1, 1, false);
    }
    for (i = 0; i < dst_grid->get_grid_size(); i ++) {
        add_remap_weights_to_sparse_matrix((long*)(&final_factor1[i]), *((long*)(&final_factor2[i])), final_factor3+i, 1, 2, false);
        add_remap_weights_to_sparse_matrix((long*)(&final_factor4[i]), *((long*)(&final_factor5[i])), coord_values_dst+i, 1, 3, false);
        temp_long_value = src_cell_index_left[i];
        add_remap_weights_to_sparse_matrix(&temp_long_value, src_cell_index_right[i], final_factor5+i, 1, 4, false);
    }
}


void Remap_operator_spline_1D::do_remap_values_caculation(double *data_values_src, double *data_values_dst, int dst_array_size)
{
    int i, j, k, m, start_index_monotonicity_range, end_index_monotonicity_range;
    int original_index1, original_index2, original_index3;
    int num_dst_cells_in_monotonicity_ranges;
    long temp_long_value1, temp_long_value2;
    double temp_double_value;
    double ratio;
    bool check_monotonicity, next_in_same_monotonicity_range;


    allocate_1D_remap_operator_common_arrays_space();
    allocate_local_arrays();

    array_size_src = remap_weights_groups[0]->get_num_weights();
    if (array_size_src == 0)
        return;
    
    for (i = 0; i < array_size_src; i ++) {
        remap_weights_groups[0]->get_weight((long*)(&array_mu[i]), &temp_long_value1, array_lambda+i, i);
        useful_src_cells_global_index[i] = temp_long_value1;
        remap_weights_groups[1]->get_weight((long*)(&coord_values_src[i]), &temp_long_value2, array_h+i, i);
    }
    for (i = 0; i < dst_grid->get_grid_size(); i ++) {
        remap_weights_groups[2]->get_weight((long*)(&final_factor1[i]), ((long*)(&final_factor2[i])), final_factor3+i, i);
        remap_weights_groups[3]->get_weight((long*)(&final_factor4[i]), ((long*)(&final_factor5[i])), coord_values_dst+i, i);
        remap_weights_groups[4]->get_weight(&temp_long_value1, &temp_long_value2, &temp_double_value, i);
        src_cell_index_left[i] = temp_long_value1;
        src_cell_index_right[i] = temp_long_value2;
    }

    for (i = 0; i < array_size_src; i++)
        array_alpha[i] = 2.0;
    
    preprocess_field_value(data_values_src);

    for (i = 1; i < array_size_src-1; i++)
        array_d[i] = 6.0*((packed_data_values_src[i+1]-packed_data_values_src[i])/array_h[i]-(packed_data_values_src[i]-packed_data_values_src[i-1])/array_h[i-1])/(array_h[i-1]+array_h[i]);

    if (!periodic) {
        array_d[0] = 0.0;
        array_d[array_size_src-1] = 0.0;
        solve_aperiodic_tridiagonal_system(array_mu, array_alpha, array_lambda, array_d, array_size_src);
    }
    else {
        array_lambda[array_size_src-1] = array_h[0]/(array_h[array_size_src-2]+array_h[0]);
        array_mu[array_size_src-1] = 1.0-array_lambda[array_size_src-1];
        array_d[array_size_src-1] = 6.0*((packed_data_values_src[1]-packed_data_values_src[0])/array_h[0]-(packed_data_values_src[array_size_src-1]-packed_data_values_src[array_size_src-2])/array_h[array_size_src-2])/(array_h[0]+array_h[array_size_src-2]);
        solve_periodic_tridiagonal_system(array_mu+1, array_alpha+1, array_lambda+1, array_d+1, array_size_src-1);
        array_d[0]=array_d[array_size_src-1];
    }

    for (i = 0; i < dst_grid->get_grid_size(); i ++) {
        if (src_cell_index_left[i] == -1 || src_cell_index_right[i] == -1)
            continue;
        if (src_cell_index_left[i] == src_cell_index_right[i])
            data_values_dst[i] = data_values_src[src_cell_index_left[i]];
        else {
            data_values_dst[i] = array_d[src_cell_index_left[i]]*final_factor1[i];
            data_values_dst[i] += array_d[src_cell_index_right[i]]*final_factor2[i];
            data_values_dst[i] += (data_values_src[src_cell_index_left[i]]-array_d[src_cell_index_left[i]]*final_factor3[i])*final_factor4[i];
            data_values_dst[i] += (data_values_src[src_cell_index_right[i]]-array_d[src_cell_index_right[i]]*final_factor3[i])*final_factor5[i];
        }
    }

    if (keep_monotonicity) {
        for (i = 0, num_dst_cells_in_monotonicity_ranges = 0; i < dst_grid->get_grid_size(); i ++) {
            if (src_cell_index_left[i] == -1 || src_cell_index_right[i] == -1)
                continue;
            if (!(coord_values_dst[i] > coord_values_src[src_cell_index_left[i]] && coord_values_dst[i] < coord_values_src[src_cell_index_right[i]]))
                continue;
            dst_cell_indexes_in_monotonicity_ranges[num_dst_cells_in_monotonicity_ranges ++] = i;
        }
        start_index_monotonicity_range = -1;
        for (i = 0; i < num_dst_cells_in_monotonicity_ranges; i ++) {
            if (start_index_monotonicity_range == -1) {
                j = 0;
                start_index_monotonicity_range = i;
                data_in_monotonicity_range[j++] = data_values_src[src_cell_index_left[dst_cell_indexes_in_monotonicity_ranges[i]]];
            }
            end_index_monotonicity_range = i;
            data_in_monotonicity_range[j++] = data_values_dst[dst_cell_indexes_in_monotonicity_ranges[i]];
            if (i == num_dst_cells_in_monotonicity_ranges - 1)
                next_in_same_monotonicity_range = false;
            else {
                original_index1 = dst_cell_indexes_in_monotonicity_ranges[start_index_monotonicity_range];
                original_index2 = dst_cell_indexes_in_monotonicity_ranges[i+1];
                next_in_same_monotonicity_range = (src_cell_index_left[original_index1] == src_cell_index_left[original_index2] && src_cell_index_right[original_index1] == src_cell_index_right[original_index2]);
            }
            if (!next_in_same_monotonicity_range) {
                original_index1 = dst_cell_indexes_in_monotonicity_ranges[start_index_monotonicity_range];
                original_index2 = dst_cell_indexes_in_monotonicity_ranges[end_index_monotonicity_range];
                EXECUTION_REPORT(REPORT_ERROR, -1, src_cell_index_left[original_index1] == src_cell_index_left[original_index2] && src_cell_index_right[original_index1] == src_cell_index_right[original_index2], 
                                 "software error: in keep monotonicity");             
                data_in_monotonicity_range[j++] = data_values_src[src_cell_index_right[original_index1]];
                check_monotonicity = true;
                for (k = 0; k < j - 1; k ++)
                    if ((data_in_monotonicity_range[k] >= data_in_monotonicity_range[k+1]) != (data_in_monotonicity_range[0] >= data_in_monotonicity_range[j-1])) {
                        check_monotonicity = false;
                        break;
                    }
                if (!check_monotonicity) {
                    original_index1 = src_cell_index_left[dst_cell_indexes_in_monotonicity_ranges[start_index_monotonicity_range]];
                    original_index2 = src_cell_index_right[dst_cell_indexes_in_monotonicity_ranges[start_index_monotonicity_range]];
                    for (k = start_index_monotonicity_range; k <= end_index_monotonicity_range; k ++) {
                        original_index3 = dst_cell_indexes_in_monotonicity_ranges[k];
                        ratio = (coord_values_dst[original_index3]-coord_values_src[original_index1]) / (coord_values_src[original_index2]-coord_values_src[original_index1]);
                        data_values_dst[original_index3] = data_values_src[original_index1]*(1-ratio) + data_values_src[original_index2]*ratio;
                    }
                }
                start_index_monotonicity_range = -1;
            }
        }
    }

    postprocess_field_value(data_values_dst);
}


void Remap_operator_spline_1D::do_src_decomp_caculation(long *decomp_map_src, const long *decomp_map_dst)
{
    EXECUTION_REPORT(REPORT_ERROR, -1, false, "Software error in Remap_operator_spline_1D::do_src_decomp_caculation: 1-D remapping algorithm should not be used to calculate src decomp");
}


Remap_operator_basis *Remap_operator_spline_1D::duplicate_remap_operator(bool fully_copy)
{
    Remap_operator_basis *duplicated_remap_operator = new Remap_operator_spline_1D();

    copy_remap_operator_basic_data(duplicated_remap_operator, fully_copy);
    ((Remap_operator_spline_1D*) duplicated_remap_operator)->initialize_1D_remap_operator();
    ((Remap_operator_spline_1D*) duplicated_remap_operator)->copy_1D_remap_operator_info(this);
    ((Remap_operator_spline_1D*) duplicated_remap_operator)->keep_monotonicity = this->keep_monotonicity;
    ((Remap_operator_spline_1D*) duplicated_remap_operator)->set_keep_monotonicity = this->set_keep_monotonicity;
    
    return duplicated_remap_operator;
}

Remap_operator_basis *Remap_operator_spline_1D::generate_parallel_remap_operator(Remap_grid_class **decomp_original_grids, int **global_cells_local_indexes_in_decomps)
{
    EXECUTION_REPORT(REPORT_ERROR, -1, false, 
                     "software error: can not generate the parallel remapping operator of the 1D spline remapping algorithm which is only used for vertical grid or time frame in the C-Coupler");    
    return NULL;
}

