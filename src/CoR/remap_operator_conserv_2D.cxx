/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#include "cor_global_data.h"
#include "remap_operator_conserv_2D.h"
#include <string.h>


void Remap_operator_conserv_2D::set_parameter(const char *parameter_name, const char *parameter_value)
{
    EXECUTION_REPORT(REPORT_ERROR, -1, enable_to_set_parameters, 
                 "the parameter of remap operator object \"%s\" must be set before using it to build remap strategy\n",
                 object_name);
    
    EXECUTION_REPORT(REPORT_ERROR, -1, false, 
                      "\"%s\" is a illegal parameter of remap operator \"%s\"\n",
                      parameter_name,
                      operator_name);
}


int Remap_operator_conserv_2D::check_parameter(const char *parameter_name, const char *parameter_value, char *error_string)
{
    return 0;
}


void Remap_operator_conserv_2D::compute_remap_weights_of_one_dst_cell(long cell_index_dst)
{
    double center_coord_values_dst[2], vertex_coord_values_dst[65536];
    int num_vertexes_dst, num_grid_dimensions_dst, i;
    long cell_index_src, *overlapping_src_cells_indexes;
    double common_sub_cell_vertexes_lons[65536], common_sub_cell_vertexes_lats[65536];
    double common_sub_cell_area[65536], weight_values[65536], sum_area;
    int num_overlapping_src_cells, num_common_sub_cell_vertexes, num_weights;


    get_cell_center_coord_values_of_dst_grid(cell_index_dst, center_coord_values_dst);
    get_cell_vertex_coord_values_of_dst_grid(cell_index_dst, &num_vertexes_dst, vertex_coord_values_dst, false);    
    EXECUTION_REPORT(REPORT_ERROR, -1, num_vertexes_dst <= 65536/2, "Software error in Remap_operator_conserv_2D::compute_remap_weights_of_one_dst_cell: too big num_vertexes_dst: %d", num_vertexes_dst);
    num_grid_dimensions_dst = current_runtime_remap_operator_grid_src->get_num_grid_dimensions();

    for (i = 0; i < num_vertexes_dst; i ++) {
        if (vertex_coord_values_dst[i*num_grid_dimensions_dst] == NULL_COORD_VALUE)
            continue;
        search_cell_in_src_grid(vertex_coord_values_dst+i*num_grid_dimensions_dst, &cell_index_src, true);
        if (cell_index_src != -1)
            break;
    }    
    if (cell_index_src == -1)
        return;

    num_overlapping_src_cells = displ_src_cells_overlap_with_dst_cells[cell_index_dst+1] - displ_src_cells_overlap_with_dst_cells[cell_index_dst];
    overlapping_src_cells_indexes = index_src_cells_overlap_with_dst_cells + displ_src_cells_overlap_with_dst_cells[cell_index_dst];

    if (num_overlapping_src_cells == 0)
        return;

    for (i = 0, sum_area = 0, num_weights = 0; i < num_overlapping_src_cells; i ++) {
        compute_common_sub_cell_of_src_cell_and_dst_cell_2D(overlapping_src_cells_indexes[i], cell_index_dst, num_common_sub_cell_vertexes, 
                                                            common_sub_cell_vertexes_lons, common_sub_cell_vertexes_lats);
        EXECUTION_REPORT(REPORT_ERROR, -1, num_common_sub_cell_vertexes <= 65536/2, "Software error in Remap_operator_conserv_2D::compute_remap_weights_of_one_dst_cell: too big num_common_sub_cell_vertexes: %d", num_common_sub_cell_vertexes);
        if (num_common_sub_cell_vertexes > 0) {
            common_sub_cell_area[num_weights] = compute_area_of_sphere_cell(num_common_sub_cell_vertexes, common_sub_cell_vertexes_lons, common_sub_cell_vertexes_lats);
            overlapping_src_cells_indexes[num_weights] = overlapping_src_cells_indexes[i];
            sum_area += common_sub_cell_area[num_weights];
            num_weights ++;
        }
    }

    if (num_weights > 0)
        EXECUTION_REPORT(REPORT_ERROR, -1, sum_area > 0, "remap software error in conserv_2D compute_remap_weights_of_one_dst_cell\n");
    for (i = 0; i < num_weights; i ++) {
        weight_values[i] = common_sub_cell_area[i]/sum_area;
    }

    add_remap_weights_to_sparse_matrix(overlapping_src_cells_indexes, cell_index_dst, weight_values, num_weights, 0, true);
}


void Remap_operator_conserv_2D::calculate_remap_weights()
{
    long cell_index_dst;
    bool dst_cell_mask;


    calculate_grids_overlaping();
    clear_remap_weight_info_in_sparse_matrix();

    for (cell_index_dst = 0; cell_index_dst < dst_grid->get_grid_size(); cell_index_dst ++) {
        if (H2D_grid_decomp_mask != NULL && !H2D_grid_decomp_mask[cell_index_dst])
            continue;
        get_cell_mask_of_dst_grid(cell_index_dst, &dst_cell_mask);
        if (!dst_cell_mask)
            continue;
        initialize_computing_remap_weights_of_one_cell();
        compute_remap_weights_of_one_dst_cell(cell_index_dst);
        clear_src_grid_cell_visiting_info();
        finalize_computing_remap_weights_of_one_cell();
    }
}


Remap_operator_conserv_2D::Remap_operator_conserv_2D(const char *object_name, int num_remap_grids, Remap_grid_class **remap_grids)
                                       : Remap_operator_basis(object_name, 
                                                              REMAP_OPERATOR_NAME_CONSERV_2D, 
                                                              2, 
                                                              true, 
                                                              true, 
                                                              true, 
                                                              num_remap_grids, 
                                                              remap_grids)
{
    num_order = 1;
    remap_weights_groups.push_back(new Remap_weight_sparse_matrix(this));
}


void Remap_operator_conserv_2D::do_remap_values_caculation(double *data_values_src, double *data_values_dst, int dst_array_size)
{
    remap_weights_groups[0]->remap_values(data_values_src, data_values_dst, dst_array_size);
}


void Remap_operator_conserv_2D::do_src_decomp_caculation(long *decomp_map_src, const long *decomp_map_dst)
{
    remap_weights_groups[0]->calc_src_decomp(decomp_map_src, decomp_map_dst);
}


Remap_operator_basis *Remap_operator_conserv_2D::duplicate_remap_operator(bool fully_copy)
{
    Remap_operator_conserv_2D *duplicated_remap_operator = new Remap_operator_conserv_2D();
    copy_remap_operator_basic_data(duplicated_remap_operator, fully_copy);
    duplicated_remap_operator->num_order = num_order;
    return duplicated_remap_operator;
}

Remap_operator_basis *Remap_operator_conserv_2D::generate_parallel_remap_operator(Remap_grid_class **decomp_original_grids, int **global_cells_local_indexes_in_decomps)
{
    Remap_operator_basis *parallel_remap_operator = this->duplicate_remap_operator(false);
    this->generate_parallel_remap_weights(parallel_remap_operator, decomp_original_grids, global_cells_local_indexes_in_decomps);
    return parallel_remap_operator;
}
