/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#include "cor_global_data.h"
#include "remap_utils_nearest_points.h"
#include "remap_grid_class.h"
#include "remap_operator_c_interface.h"
#include "quick_sort.h"
#include "remap_common_utils.h"
#include <math.h>


double calculate_distance_of_two_points_2D(double point1_coord1_value,
                                           double point1_coord2_value,
                                           double point2_coord1_value,
                                           double point2_coord2_value,
                                           bool is_sphere_grid)
{
    double distance, temp_value, diff1, diff2;
    double rad_lon1, rad_lon2, rad_lat1, rad_lat2;


    if (is_sphere_grid && (point2_coord2_value == 90 || point2_coord2_value == -90))
        if (point1_coord2_value == point2_coord2_value)
            return 0.0;

    if (point1_coord1_value == point2_coord1_value && point1_coord2_value == point2_coord2_value)
        return 0.0;
    
    if (is_sphere_grid) {
        rad_lon1 = DEGREE_TO_RADIAN(point1_coord1_value);
        rad_lon2 = DEGREE_TO_RADIAN(point2_coord1_value);
        rad_lat1 = DEGREE_TO_RADIAN(point1_coord2_value);
        rad_lat2 = DEGREE_TO_RADIAN(point2_coord2_value);        
        temp_value = cos(rad_lat1)*cos(rad_lat2)*cos(rad_lon1-rad_lon2) + sin(rad_lat1)*sin(rad_lat2);
        if (temp_value > 1)
            temp_value = 1;
        if (temp_value < -1)
            temp_value = -1.0;
        distance = acos(temp_value);
    }
    else {
        diff1 = compute_difference_of_two_coord_values(point1_coord1_value, point2_coord1_value, 0);
        diff2 = compute_difference_of_two_coord_values(point1_coord2_value, point2_coord2_value, 1);
        distance = sqrt(diff1*diff1 + diff2*diff2);
    }

    return distance;
}


void compute_dist_remap_weights_of_one_dst_cell(long dst_cell_index,
                                                int num_nearest_points,
                                                double num_power,
                                                double *threshold_distance,
                                                double *found_nearest_points_distance,
                                                long *found_nearest_points_src_indexes,
                                                double *weigt_values_of_one_dst_cell,
                                                bool is_sphere_grid,
                                                bool enable_extrapolate)
{
    bool successful = false, dst_cell_mask;
    long src_cell_index;
    double sum_wgt_values, src_cell_center_values[256], dst_cell_center_values[256], dst_cell_vertex_values[65536];
    int i, num_points_within_threshold_dist, num_vertexes_dst;
    double current_dist;

    
    get_cell_mask_of_dst_grid(dst_cell_index, &dst_cell_mask);
    if (!dst_cell_mask)
        return;

    get_cell_center_coord_values_of_dst_grid(dst_cell_index, dst_cell_center_values);
    get_cell_vertex_coord_values_of_dst_grid(dst_cell_index, &num_vertexes_dst, dst_cell_vertex_values, true);
    EXECUTION_REPORT(REPORT_ERROR, -1, num_vertexes_dst <= 65536/2, "Software error in compute_dist_remap_weights_of_one_dst_cell: too big number of dst vertexes: %d", num_vertexes_dst);
    
    if (num_vertexes_dst > 0 && (!enable_extrapolate && !have_overlapped_src_cells_for_dst_cell(dst_cell_index)))
        return;

    if (num_vertexes_dst == 0) {
        search_cell_in_src_grid(dst_cell_center_values, &src_cell_index, false);
        if (src_cell_index == -1 && (!enable_extrapolate))
            return;
    }

    get_current_grid2D_search_engine(true)->search_nearest_points_var_number(num_nearest_points, dst_cell_center_values[0], dst_cell_center_values[1], 
                                                                                                      num_points_within_threshold_dist, found_nearest_points_src_indexes, found_nearest_points_distance, true);

    if (num_nearest_points > num_points_within_threshold_dist)
        num_nearest_points = num_points_within_threshold_dist;
    
    if (found_nearest_points_distance[0] == 0.0) {
        weigt_values_of_one_dst_cell[0] = 1.0;
        add_remap_weights_to_sparse_matrix(found_nearest_points_src_indexes, dst_cell_index, weigt_values_of_one_dst_cell, 1, 0, true);
    }
    else {
        sum_wgt_values = 0.0;
        for (i = 0; i < num_nearest_points; i ++) {
            weigt_values_of_one_dst_cell[i] = 1.0 / pow(found_nearest_points_distance[i], num_power);
            sum_wgt_values += weigt_values_of_one_dst_cell[i];
        }
        weigt_values_of_one_dst_cell[num_nearest_points-1] = 1.0;
        for (i = 0; i < num_nearest_points-1; i ++) {
            weigt_values_of_one_dst_cell[i] = weigt_values_of_one_dst_cell[i] / sum_wgt_values;        
            weigt_values_of_one_dst_cell[num_nearest_points-1] -= weigt_values_of_one_dst_cell[i];
        }
        add_remap_weights_to_sparse_matrix(found_nearest_points_src_indexes, dst_cell_index, weigt_values_of_one_dst_cell, num_nearest_points, 0, true);
    }
}

