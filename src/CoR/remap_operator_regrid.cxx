/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#include "cor_global_data.h"
#include "remap_operator_regrid.h"
#include <string.h>
#include <math.h>


void Remap_operator_regrid::set_parameter(const char *parameter_name, const char *parameter_value)
{
    EXECUTION_REPORT(REPORT_ERROR, -1, false, "remap software error in set_parameter of Remap_operator_regrid\n");
}


int Remap_operator_regrid::check_parameter(const char *parameter_name, const char *parameter_value, char *error_string)
{
    EXECUTION_REPORT(REPORT_ERROR, -1, false, "remap software error in set_parameter of Remap_operator_regrid\n");
    return 0;
}


void Remap_operator_regrid::calculate_remap_weights()
{
    EXECUTION_REPORT(REPORT_ERROR, -1, false, "remap software error in calculate_remap_weights of Remap_operator_regrid\n");
}


Remap_operator_regrid::Remap_operator_regrid(const char *object_name, int num_remap_grids, Remap_grid_class **remap_grids, int num_dimensions)
                                       : Remap_operator_basis(object_name, 
                                                              REMAP_OPERATOR_NAME_REGRID, 
                                                              num_dimensions, 
                                                              true, 
                                                              false, 
                                                              true,
                                                              num_remap_grids, 
                                                              remap_grids)
{
}


void Remap_operator_regrid::do_remap_values_caculation(double *data_values_src, double *data_values_dst, int dst_array_size)
{
    EXECUTION_REPORT(REPORT_ERROR, -1, false, "remap software error in do_remap_values_caculation of Remap_operator_regrid\n");
}


void Remap_operator_regrid::do_src_decomp_caculation(long *decomp_map_src, const long *decomp_map_dst)
{
    EXECUTION_REPORT(REPORT_ERROR, -1, false, "remap software error in do_src_decomp_caculation of Remap_operator_regrid\n");
}


Remap_operator_basis *Remap_operator_regrid::duplicate_remap_operator(bool fully_copy)
{
    EXECUTION_REPORT(REPORT_ERROR, -1, false, "remap software error in duplicate_remap_operator of Remap_operator_regrid\n");
    return NULL;
}


Remap_operator_basis *Remap_operator_regrid::generate_parallel_remap_operator(Remap_grid_class **decomp_original_grids, int **global_cells_local_indexes_in_decomps)
{
    EXECUTION_REPORT(REPORT_ERROR, -1, false, "remap software error in generate_parallel_remap_operator of Remap_operator_regrid\n");

    return NULL;
}


void Remap_operator_regrid::compute_remap_weights_of_one_dst_cell(long index_dst_cell)
{
}

