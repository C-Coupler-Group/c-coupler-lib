/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#include "cor_global_data.h"
#include "remap_operator_smooth.h"
#include <string.h>


void Remap_operator_smooth::set_parameter(const char *parameter_name, const char *parameter_value)
{
    EXECUTION_REPORT(REPORT_ERROR, -1, enable_to_set_parameters, 
                 "the parameter of remap operator object \"%s\" must be set before using it to build remap strategy\n",
                 object_name);
    
    if (words_are_the_same(parameter_name, "period"))
        sscanf(parameter_value, "%d", &num_period);
    else EXECUTION_REPORT(REPORT_ERROR, -1, false, 
                      "\"%s\" is a illegal parameter of remap operator \"%s\"\n",
                      parameter_name,
                      operator_name);
}


int Remap_operator_smooth::check_parameter(const char *parameter_name, const char *parameter_value, char *error_string)
{
    EXECUTION_REPORT(REPORT_ERROR, -1, false, "Software error in Remap_operator_smooth::check_parameter");
    return 0;
}


void Remap_operator_smooth::compute_remap_weights_of_one_dst_cell(long index_dst_cell)
{
}


void Remap_operator_smooth::calculate_remap_weights()
{
}


Remap_operator_smooth::Remap_operator_smooth(const char *object_name, int num_remap_grids, Remap_grid_class **remap_grids)
                                       : Remap_operator_basis(object_name, 
                                                              REMAP_OPERATOR_NAME_SMOOTH, 
                                                              1, 
                                                              false, 
                                                              false,
                                                              true,
                                                              num_remap_grids, 
                                                              remap_grids)
{
    remap_weights_groups.push_back(new Remap_weight_sparse_matrix(this));
}


void Remap_operator_smooth::do_remap_values_caculation(double *data_values_src, double *data_values_dst, int dst_array_size)
{
    remap_weights_groups[0]->remap_values(data_values_src, data_values_dst, dst_array_size);
}


void Remap_operator_smooth::do_src_decomp_caculation(long *decomp_map_src, const long *decomp_map_dst)
{
    remap_weights_groups[0]->calc_src_decomp(decomp_map_src, decomp_map_dst);
}


Remap_operator_basis *Remap_operator_smooth::duplicate_remap_operator(bool fully_copy)
{
    Remap_operator_basis *duplicated_remap_operator = new Remap_operator_smooth();
    copy_remap_operator_basic_data(duplicated_remap_operator, fully_copy);
    return duplicated_remap_operator;
}


Remap_operator_basis *Remap_operator_smooth::generate_parallel_remap_operator(Remap_grid_class **decomp_original_grids, int **global_cells_local_indexes_in_decomps)
{
    Remap_operator_basis *parallel_remap_operator = this->duplicate_remap_operator(false);
    this->generate_parallel_remap_weights(parallel_remap_operator, decomp_original_grids, global_cells_local_indexes_in_decomps);
    return parallel_remap_operator;
}
