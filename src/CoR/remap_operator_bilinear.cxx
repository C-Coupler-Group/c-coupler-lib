/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#include "cor_global_data.h"
#include "remap_operator_bilinear.h"
#include "remap_utils_nearest_points.h"
#include "remap_common_utils.h"
#include "quick_sort.h"
#include <string.h>
#include <math.h>


void Remap_operator_bilinear::set_parameter(const char *parameter_name, const char *parameter_value)
{
    EXECUTION_REPORT(REPORT_ERROR, -1, enable_to_set_parameters, 
                 "the parameter of remap operator object \"%s\" must be set before using it to build remap strategy\n",
                 object_name);
    if (words_are_the_same(parameter_name, "enable_extrapolate")) {
        if (words_are_the_same(parameter_value, "true"))
            enable_extrapolate = true;
        else if (words_are_the_same(parameter_value, "false"))
            enable_extrapolate = false;
        else EXECUTION_REPORT(REPORT_ERROR, -1, false, "The parameter value must be \"true\" or \"false\"\n");
    }
    else EXECUTION_REPORT(REPORT_ERROR, -1, false, "bilinear algorithm does not have the parameter to be set\n");
}


int Remap_operator_bilinear::check_parameter(const char *parameter_name, const char *parameter_value, char *error_string)
{
    int check_result = 0;
    if (words_are_the_same(parameter_name, "enable_extrapolate")) {
        check_result = 1;
        if (words_are_the_same(parameter_value, "true") || words_are_the_same(parameter_value, "false"))
            check_result = 3;
        else sprintf(error_string, "The parameter value must be \"true\" or \"false\"");
    }
    
    return check_result;
}


Remap_operator_bilinear::Remap_operator_bilinear()
{
    found_nearest_points_distance = NULL;
    found_nearest_points_src_indexes = NULL;
    weigt_values_of_one_dst_cell = NULL;
    enable_extrapolate = false;
}


Remap_operator_bilinear::Remap_operator_bilinear(const char *object_name, int num_remap_grids, Remap_grid_class **remap_grids)
                                       : Remap_operator_basis(object_name, 
                                                              REMAP_OPERATOR_NAME_BILINEAR, 
                                                              2, 
                                                              true, 
                                                              false, 
                                                              true, 
                                                              num_remap_grids, 
                                                              remap_grids)
{
    remap_weights_groups.push_back(new Remap_weight_sparse_matrix(this));
    max_num_found_nearest_points = 256;
    num_nearest_points = 4;
    num_power = 1.0;
    found_nearest_points_distance = new double [src_grid->get_grid_size()];
    found_nearest_points_src_indexes = new long [src_grid->get_grid_size()];
    weigt_values_of_one_dst_cell = new double [max_num_found_nearest_points];
    enable_extrapolate = false;
}


Remap_operator_bilinear::~Remap_operator_bilinear()
{
    if (found_nearest_points_distance != NULL)
        delete [] found_nearest_points_distance;
    if (found_nearest_points_src_indexes != NULL)
        delete [] found_nearest_points_src_indexes;
    if (weigt_values_of_one_dst_cell != NULL)
        delete [] weigt_values_of_one_dst_cell;
}


void Remap_operator_bilinear::calculate_remap_weights()
{
    iterative_threshold_distance = 1.0/6000.0;
    calculate_grids_overlaping();
    clear_remap_weight_info_in_sparse_matrix();
    
    for (long dst_cell_index = 0; dst_cell_index < dst_grid->get_grid_size(); dst_cell_index++) {
        if (H2D_grid_decomp_mask != NULL && !H2D_grid_decomp_mask[dst_cell_index])
            continue;
        initialize_computing_remap_weights_of_one_cell();
        compute_remap_weights_of_one_dst_cell(dst_cell_index);    
        finalize_computing_remap_weights_of_one_cell();
    }
}


int Remap_operator_bilinear::search_at_least_16_nearnest_src_points_for_bilinear(double *dst_cell_center_values,
                                                                                 long src_cell_index,
                                                                                 double &current_threshold_distance,
                                                                                 double &near_optimal_threshold_distance)
{
    int num_points_within_threshold_distance = 0;
    double eps = 2.0e-9;


    while (num_points_within_threshold_distance < 16) {
        num_points_within_threshold_distance = 0;
        get_current_grid2D_search_engine(true)->search_nearest_points_var_distance(current_threshold_distance, dst_cell_center_values[0], dst_cell_center_values[1], 
                                                                                                              num_points_within_threshold_distance, found_nearest_points_src_indexes, found_nearest_points_distance, true);
        if (num_points_within_threshold_distance >= 4 && near_optimal_threshold_distance == 0.0)
            near_optimal_threshold_distance = current_threshold_distance * sqrt(((double)4)/((double)num_points_within_threshold_distance));
        if (num_points_within_threshold_distance == 0)
            current_threshold_distance *= 2;
        else current_threshold_distance *= 1.1;     
    if (num_points_within_threshold_distance > 0 && found_nearest_points_distance[0] <= eps)
        break;

    }

    return num_points_within_threshold_distance;
}


int Remap_operator_bilinear::compute_quadrant_of_src_point(double* dst_cell_center_values, 
                                                           double *src_cell_center_values)
{
    double vector_angle;
    int quadrant_id;


    vector_angle = compute_vector_angle(dst_cell_center_values, src_cell_center_values);
    if (vector_angle < 0)
        vector_angle += 2*PI;
    quadrant_id = (int)(vector_angle*2/PI);
    if (quadrant_id == 4)
        quadrant_id = 0;

    EXECUTION_REPORT(REPORT_ERROR, -1, quadrant_id >= 0 && quadrant_id < 5, "remap software error in compute_quadrant_of_src_point\n");    
    return quadrant_id;
}


bool Remap_operator_bilinear::are_three_points_on_the_same_line(double *three_points_coord1_values, double *three_points_coord2_values)
{
    double diff1, diff2, diff3, diff4;


    diff1 = compute_difference_of_two_coord_values(three_points_coord2_values[0], three_points_coord2_values[2], 1);
    diff2 = compute_difference_of_two_coord_values(three_points_coord1_values[1], three_points_coord1_values[2], 0);
    diff3 = compute_difference_of_two_coord_values(three_points_coord2_values[1], three_points_coord2_values[2], 1);
    diff4 = compute_difference_of_two_coord_values(three_points_coord1_values[0], three_points_coord1_values[2], 0);
    return diff1 * diff2 == diff3 * diff4;
}


bool Remap_operator_bilinear::get_near_optimal_bilinear_box_recursively(double **distances_of_src_points_in_each_quadrant,
                                                                        long **indexes_of_src_points_in_each_quadrant,
                                                                        const int *num_src_points_in_each_quadrant,
                                                                        int *iter_num_src_points_in_each_quadrant,
                                                                        const double* dst_cell_center_values,
                                                                        long *index_of_selected_src_point_in_each_quadrant,
                                                                        int recursion_index)
{
    int i, j;
    double src_cell_center_values[2];
    double bilinear_vertex_coord1_values[4], bilinear_vertex_coord2_values[4];
    double triangle_vertex_coord1_values[3], triangle_vertex_coord2_values[3];
    double distance_of_selected_src_point_in_each_quadrant[4];
    bool have_points_on_the_same_line;
    int index_of_quadrant_id_for_distance_sorting[4];
    int new_iter_num_src_points_in_each_quadrant[4];
    double eps = 2.0e-9;


    for (i = 0; i < 4; i ++) {
        index_of_selected_src_point_in_each_quadrant[i] = indexes_of_src_points_in_each_quadrant[i][iter_num_src_points_in_each_quadrant[i]];
        get_cell_center_coord_values_of_src_grid(index_of_selected_src_point_in_each_quadrant[i], src_cell_center_values);
        bilinear_vertex_coord1_values[i] = src_cell_center_values[0];
        bilinear_vertex_coord2_values[i] = src_cell_center_values[1];
    }

    sort_polygon_vertexes(dst_cell_center_values[0],
                          dst_cell_center_values[1],
                          bilinear_vertex_coord1_values,
                          bilinear_vertex_coord2_values,
                          index_of_selected_src_point_in_each_quadrant,
                          4);     
    
    if (is_point_in_2D_cell(dst_cell_center_values[0], 
                            dst_cell_center_values[1],
                            bilinear_vertex_coord1_values,
                            bilinear_vertex_coord2_values,
                            4,
                            is_coord_unit_degree[0],
                            is_coord_unit_degree[1],
                            false)) {
        have_points_on_the_same_line = false;
        for (i = 0; i < 4; i ++) {
            for (j = 0; j < 3; j ++) {
                triangle_vertex_coord1_values[j] = bilinear_vertex_coord1_values[(i+j)%4];
                triangle_vertex_coord2_values[j] = bilinear_vertex_coord2_values[(i+j)%4];                    
            }
            have_points_on_the_same_line |= are_three_points_on_the_same_line(triangle_vertex_coord1_values, triangle_vertex_coord2_values);
            if (have_points_on_the_same_line)
                break;        
            triangle_vertex_coord1_values[0] = bilinear_vertex_coord1_values[i];
            triangle_vertex_coord2_values[0] = bilinear_vertex_coord2_values[i];
            triangle_vertex_coord1_values[1] = dst_cell_center_values[0];
            triangle_vertex_coord2_values[1] = dst_cell_center_values[1];
            triangle_vertex_coord1_values[2] = bilinear_vertex_coord1_values[(i+1)%4];
            triangle_vertex_coord2_values[2] = bilinear_vertex_coord2_values[(i+1)%4];
            have_points_on_the_same_line |= are_three_points_on_the_same_line(triangle_vertex_coord1_values, triangle_vertex_coord2_values);                
            if (have_points_on_the_same_line)
                break;
        }
        
        if (!have_points_on_the_same_line) 
            return true;
    }

    for (i = 0; i < 4; i ++) {
        index_of_selected_src_point_in_each_quadrant[i] = indexes_of_src_points_in_each_quadrant[i][iter_num_src_points_in_each_quadrant[i]];
        distance_of_selected_src_point_in_each_quadrant[i] = distances_of_src_points_in_each_quadrant[i][iter_num_src_points_in_each_quadrant[i]];
        get_cell_center_coord_values_of_src_grid(index_of_selected_src_point_in_each_quadrant[i], src_cell_center_values);
        new_iter_num_src_points_in_each_quadrant[i] = iter_num_src_points_in_each_quadrant[i];
        index_of_quadrant_id_for_distance_sorting[i] = i;
    }
    do_quick_sort(distance_of_selected_src_point_in_each_quadrant, index_of_quadrant_id_for_distance_sorting, 0, 3);

    for (i = 3; i >= recursion_index; i --) {
        j = index_of_quadrant_id_for_distance_sorting[i];    
        new_iter_num_src_points_in_each_quadrant[j] ++;
        if (new_iter_num_src_points_in_each_quadrant[j] < num_src_points_in_each_quadrant[j] &&
            get_near_optimal_bilinear_box_recursively(distances_of_src_points_in_each_quadrant, 
                                                      indexes_of_src_points_in_each_quadrant, 
                                                      num_src_points_in_each_quadrant, 
                                                      new_iter_num_src_points_in_each_quadrant,
                                                      dst_cell_center_values,
                                                      index_of_selected_src_point_in_each_quadrant,
                                                      i))
            return true;
        new_iter_num_src_points_in_each_quadrant[j] --;
    }
    
    return false;
}


bool Remap_operator_bilinear::get_near_optimal_bilinear_box(double* dst_cell_center_values, 
                                                            int num_points_within_threshold_dist,
                                                            long *index_of_selected_src_point_in_each_quadrant)
{
    int i, quadrant_id;
    double src_cell_center_values[2];
    long indexes_of_src_points_in_each_quadrant[4][256], *pointer_indexes_of_src_points_in_each_quadrant[4];
    double distances_of_src_points_in_each_quadrant[4][256], *pointer_distances_of_src_points_in_each_quadrant[4];
    int num_src_points_in_each_quadrant[4], iter_num_src_points_in_each_quadrant[4];
    

    EXECUTION_REPORT(REPORT_ERROR, -1, num_points_within_threshold_dist <= 256, "remap software error in get_nearest_point_in_each_of_three_quadrants\n");

    for (i = 0; i < 4; i ++) {
        num_src_points_in_each_quadrant[i] = 0;
        iter_num_src_points_in_each_quadrant[i] = 0;
        pointer_distances_of_src_points_in_each_quadrant[i] = distances_of_src_points_in_each_quadrant[i];
        pointer_indexes_of_src_points_in_each_quadrant[i] = indexes_of_src_points_in_each_quadrant[i];
    }
    
    for (i = 0; i < num_points_within_threshold_dist; i ++) {
        get_cell_center_coord_values_of_src_grid(found_nearest_points_src_indexes[i], src_cell_center_values);
        quadrant_id = compute_quadrant_of_src_point(dst_cell_center_values, 
                                                    src_cell_center_values);
        indexes_of_src_points_in_each_quadrant[quadrant_id][num_src_points_in_each_quadrant[quadrant_id]] = found_nearest_points_src_indexes[i];
        distances_of_src_points_in_each_quadrant[quadrant_id][num_src_points_in_each_quadrant[quadrant_id]] = found_nearest_points_distance[i];
        num_src_points_in_each_quadrant[quadrant_id] ++;
    }
    
    for (i = 0; i < 4; i ++)
        if (num_src_points_in_each_quadrant[i] == 0)
            break;
    if (i < 4) 
        return false;

    for (i = 0; i < 4; i ++)
        do_quick_sort(distances_of_src_points_in_each_quadrant[i], indexes_of_src_points_in_each_quadrant[i], 0, num_src_points_in_each_quadrant[i]-1);

    return get_near_optimal_bilinear_box_recursively(pointer_distances_of_src_points_in_each_quadrant, 
                                                     pointer_indexes_of_src_points_in_each_quadrant, 
                                                     num_src_points_in_each_quadrant, 
                                                     iter_num_src_points_in_each_quadrant,
                                                     dst_cell_center_values,
                                                     index_of_selected_src_point_in_each_quadrant,
                                                     0);
}


void Remap_operator_bilinear::compute_remap_weights_of_one_dst_cell(long dst_cell_index)
{
    bool find_bilinear_box, dst_cell_mask,  src_cell_mask;
    int i, j, num_points_within_threshold_distance;
    long src_cell_index;
    double dst_cell_center_values[2], src_cell_center_values[2];
    double current_threshold_distance, near_optimal_threshold_distance;
    long best_index_of_selected_src_point_in_each_quadrant[4];
    double bilinear_box_vertex_coord1_values[4], bilinear_box_vertex_coord2_values[4];
    long bilinear_box_vertexes_src_cell_indexes[4];
    double wgt_ratio_u, wgt_ratio_v;
    double bilinear_wgt_values[4];
    double eps = 2.0e-7;
    int num_vertexes_dst;
    double vertex_coord_values_dst[65536];
    

    initialize_computing_remap_weights_of_one_cell();

    /*  When the mask of dst cell is false, it is unnecessary to compute the corresponding weight values
      */
    get_cell_mask_of_dst_grid(dst_cell_index, &dst_cell_mask);
    if (!dst_cell_mask) 
        return;

    /* When no src cell contains the center of the dst cell, we use inverse distance weight to compute weight values
      */
    get_cell_center_coord_values_of_dst_grid(dst_cell_index, dst_cell_center_values);
    get_cell_vertex_coord_values_of_dst_grid(dst_cell_index, &num_vertexes_dst, vertex_coord_values_dst, true);
    EXECUTION_REPORT(REPORT_ERROR, -1, num_vertexes_dst <= 65536/2, "Software error in Remap_operator_bilinear::compute_remap_weights_of_one_dst_cell: too big number of num_vertexes_dst %d", num_vertexes_dst);

    if (num_vertexes_dst > 0 && (!enable_extrapolate && !have_overlapped_src_cells_for_dst_cell(dst_cell_index)))
        return;

    search_cell_in_src_grid(dst_cell_center_values, &src_cell_index, false);
    if (src_cell_index != -1)
        get_cell_mask_of_src_grid(src_cell_index, &src_cell_mask);

    if (num_vertexes_dst == 0 && src_cell_index == -1 && (!enable_extrapolate))
        return;
    if (num_vertexes_dst == 0 && (!enable_extrapolate && !src_cell_mask))
        return;
    
    iterative_threshold_distance = 1.0/6000.0;

    if (src_cell_index == -1 || !src_cell_mask) {
        compute_dist_remap_weights_of_one_dst_cell(dst_cell_index, 
                                                   num_nearest_points,
                                                   num_power,
                                                   &iterative_threshold_distance,
                                                   found_nearest_points_distance,
                                                   found_nearest_points_src_indexes,
                                                   weigt_values_of_one_dst_cell,
                                                   get_is_sphere_grid(),
                                                   enable_extrapolate);
        return;
    }

    /* When the centers of src cell and dst cell are the same, we directly compute the weight value
      */
    get_cell_center_coord_values_of_src_grid(src_cell_index, src_cell_center_values);
    if (calculate_distance_of_two_points_2D(dst_cell_center_values[0], 
                                            dst_cell_center_values[1],
                                            src_cell_center_values[0], 
                                            src_cell_center_values[1],
                                            get_is_sphere_grid()) <= eps) {
        weigt_values_of_one_dst_cell[0] = 1.0;
        add_remap_weights_to_sparse_matrix(&src_cell_index, dst_cell_index, weigt_values_of_one_dst_cell, 1, 0, true);
        return;
    }

    find_bilinear_box = false;
    current_threshold_distance = iterative_threshold_distance;
    near_optimal_threshold_distance = 0.0;
    while (1) {
        num_points_within_threshold_distance = search_at_least_16_nearnest_src_points_for_bilinear(dst_cell_center_values, 
                                                                                                   src_cell_index, 
                                                                                                   current_threshold_distance, 
                                                                                                   near_optimal_threshold_distance);
        if (found_nearest_points_distance[0] <= eps) {
            weigt_values_of_one_dst_cell[0] = 1.0;
            add_remap_weights_to_sparse_matrix(&found_nearest_points_src_indexes[0], dst_cell_index, weigt_values_of_one_dst_cell, 1, 0, true);
            return;
        }
        if (num_points_within_threshold_distance > max_num_found_nearest_points)
            break;        
        if (get_near_optimal_bilinear_box(dst_cell_center_values, 
                                          num_points_within_threshold_distance,
                                          best_index_of_selected_src_point_in_each_quadrant)) {
            for (j = 0; j < 4; j ++) {
                get_cell_center_coord_values_of_src_grid(best_index_of_selected_src_point_in_each_quadrant[j], src_cell_center_values);
                bilinear_box_vertex_coord1_values[j] = src_cell_center_values[0];
                bilinear_box_vertex_coord2_values[j] = src_cell_center_values[1];
                bilinear_box_vertexes_src_cell_indexes[j] = best_index_of_selected_src_point_in_each_quadrant[j];
            }
            find_bilinear_box = true;                   
            break;
        }
    }

    EXECUTION_REPORT(REPORT_ERROR, -1, near_optimal_threshold_distance > 0, "remap software error1 in blinear compute_remap_weights_of_one_dst_cell\n");
    iterative_threshold_distance = near_optimal_threshold_distance;
    
    if (!find_bilinear_box) {
        compute_dist_remap_weights_of_one_dst_cell(dst_cell_index, 
                                                   num_nearest_points,
                                                   num_power,
                                                   &iterative_threshold_distance,
                                                   found_nearest_points_distance,
                                                   found_nearest_points_src_indexes,
                                                   weigt_values_of_one_dst_cell,
                                                   get_is_sphere_grid(),
                                                   enable_extrapolate);
    }
    else {
        solve_two_bilinear_ratios(bilinear_box_vertexes_src_cell_indexes, dst_cell_center_values, wgt_ratio_u, wgt_ratio_v);
        bilinear_wgt_values[0] = (1-wgt_ratio_u) * (1-wgt_ratio_v);
        bilinear_wgt_values[1] = wgt_ratio_u * (1-wgt_ratio_v);
        bilinear_wgt_values[2] = wgt_ratio_u * wgt_ratio_v;
        bilinear_wgt_values[3] = (1-wgt_ratio_u) * wgt_ratio_v;
        add_remap_weights_to_sparse_matrix(bilinear_box_vertexes_src_cell_indexes, dst_cell_index, bilinear_wgt_values, 4, 0, true);
    }
}


double Remap_operator_bilinear::compute_cross_product_of_counter_lines(double *bilinear_box_vertexes_coord1_values, double *bilinear_box_vertexes_coord2_values)
{
    return compute_difference_of_two_coord_values(bilinear_box_vertexes_coord1_values[3],bilinear_box_vertexes_coord1_values[0],0) *
           compute_difference_of_two_coord_values(bilinear_box_vertexes_coord2_values[2],bilinear_box_vertexes_coord2_values[1],1) -
           compute_difference_of_two_coord_values(bilinear_box_vertexes_coord1_values[2],bilinear_box_vertexes_coord1_values[1],0) *
           compute_difference_of_two_coord_values(bilinear_box_vertexes_coord2_values[3],bilinear_box_vertexes_coord2_values[0],1);
}


void Remap_operator_bilinear::bilinear_ratios_solution1(double *dst_point_coord_values, 
                                                        double *bilinear_box_vertexes_coord1_values, 
                                                        double *bilinear_box_vertexes_coord2_values, 
                                                        double &ratio_u, 
                                                        double &ratio_v)
{
    double x_diff_01, y_diff_01, x_diff_03, y_diff_03, x_diff_32, y_diff_32, x_diff_04, y_diff_04;
    double coord_values_P5[2], coord_values_P6[2];
    double dist_54, dist_56;


    /*
          (x0,y0)                  (x1,y1)
             p0--P5------p1
             |    |         |
             |    .P4(x,y)   |
             |    |         |
             |    |         |
             p3--P6------p2
          (x3,y3)                  (x2,y2)

            P4 is dst point and P0~P3 are 4 vertexes of bilinear box. P0P3 parallels to P1P2. ratio_u=|P5P0|/|P1P0|, ratio_v=|P4P5|/|P6P5|
    */

    x_diff_01 = compute_difference_of_two_coord_values(bilinear_box_vertexes_coord1_values[1], bilinear_box_vertexes_coord1_values[0], 0);
    y_diff_01 = compute_difference_of_two_coord_values(bilinear_box_vertexes_coord2_values[1], bilinear_box_vertexes_coord2_values[0], 1);
    x_diff_03 = compute_difference_of_two_coord_values(bilinear_box_vertexes_coord1_values[3], bilinear_box_vertexes_coord1_values[0], 0);
    y_diff_03 = compute_difference_of_two_coord_values(bilinear_box_vertexes_coord2_values[3], bilinear_box_vertexes_coord2_values[0], 1);
    x_diff_04 = compute_difference_of_two_coord_values(dst_point_coord_values[0], bilinear_box_vertexes_coord1_values[0], 0);
    y_diff_04 = compute_difference_of_two_coord_values(dst_point_coord_values[1], bilinear_box_vertexes_coord2_values[0], 1);
    x_diff_32 = compute_difference_of_two_coord_values(bilinear_box_vertexes_coord1_values[2], bilinear_box_vertexes_coord1_values[3], 0);
    y_diff_32 = compute_difference_of_two_coord_values(bilinear_box_vertexes_coord2_values[2], bilinear_box_vertexes_coord2_values[3], 1);

    ratio_u = (y_diff_03*x_diff_04-x_diff_03*y_diff_04)/(x_diff_01*y_diff_03-y_diff_01*x_diff_03);
    EXECUTION_REPORT(REPORT_ERROR, -1, ratio_u >= 0 && ratio_u <= 1, "remap software error1 in bilinear_ratios_solution1\n");
    
    coord_values_P5[0] = bilinear_box_vertexes_coord1_values[0] + ratio_u*x_diff_01;
    coord_values_P5[1] = bilinear_box_vertexes_coord2_values[0] + ratio_u*y_diff_01;
    coord_values_P6[0] = bilinear_box_vertexes_coord1_values[3] + ratio_u*x_diff_32;
    coord_values_P6[1] = bilinear_box_vertexes_coord2_values[3] + ratio_u*y_diff_32;
    dist_54 = calculate_distance_of_two_points_2D(dst_point_coord_values[0],
                                                  dst_point_coord_values[1],
                                                  coord_values_P5[0],
                                                  coord_values_P5[1],
                                                  false);
    dist_56 = calculate_distance_of_two_points_2D(coord_values_P6[0],
                                                  coord_values_P6[1],
                                                  coord_values_P5[0],
                                                  coord_values_P5[1],
                                                  false);
    ratio_v = dist_54 / dist_56;
    EXECUTION_REPORT(REPORT_ERROR, -1, ratio_v >= 0 && ratio_v <= 1, "remap software error2 in bilinear_ratios_solution1\n");
}


void Remap_operator_bilinear::bilinear_one_ratio_solution_of_quadratic_equation(double *dst_point_coord_values, 
                                                                                double *bilinear_box_vertexes_coord1_values, 
                                                                                double *bilinear_box_vertexes_coord2_values, 
                                                                                double &ratio_u)
{
    double x_diff_01, y_diff_01, x_diff_32, y_diff_32, x_diff_04, y_diff_04, x_diff_43, y_diff_43;
    double coef_A, coef_B, coef_C;
    double ratio_u1, ratio_u2, ratio_u_false;
	double eps = 1.0e-13;
        
    /*
          (x0,y0)                  (x1,y1)
             p0--P5------p1
             |     |         |
             |     .P4(x,y)   |
             |     |         |
             |     |         |
             p3--P6------p2
          (x3,y3)                  (x2,y2)
        
            P4 is dst point and P0~P3 are 4 vertexes of bilinear box. P0P3 does not parallel to P1P2 and P0P1 does not parallel to P3P2. 
            ratio_u=|P5P0|/|P1P0|
    */
    
    x_diff_01 = compute_difference_of_two_coord_values(bilinear_box_vertexes_coord1_values[1], bilinear_box_vertexes_coord1_values[0], 0);
    y_diff_01 = compute_difference_of_two_coord_values(bilinear_box_vertexes_coord2_values[1], bilinear_box_vertexes_coord2_values[0], 1);
    x_diff_32 = compute_difference_of_two_coord_values(bilinear_box_vertexes_coord1_values[2], bilinear_box_vertexes_coord1_values[3], 0);
    y_diff_32 = compute_difference_of_two_coord_values(bilinear_box_vertexes_coord2_values[2], bilinear_box_vertexes_coord2_values[3], 1);
    x_diff_04 = compute_difference_of_two_coord_values(dst_point_coord_values[0], bilinear_box_vertexes_coord1_values[0], 0);
    y_diff_04 = compute_difference_of_two_coord_values(dst_point_coord_values[1], bilinear_box_vertexes_coord2_values[0], 1);
    x_diff_43 = compute_difference_of_two_coord_values(bilinear_box_vertexes_coord1_values[3], dst_point_coord_values[0], 0);
    y_diff_43 = compute_difference_of_two_coord_values(bilinear_box_vertexes_coord2_values[3], dst_point_coord_values[1], 1);
    
    coef_A = y_diff_01*x_diff_32 - x_diff_01*y_diff_32;
    coef_B = x_diff_04*y_diff_32 + x_diff_43*y_diff_01 - y_diff_04*x_diff_32 - y_diff_43*x_diff_01;
    coef_C = x_diff_04*y_diff_43 - y_diff_04*x_diff_43;

    EXECUTION_REPORT(REPORT_ERROR, -1, coef_A != 0, "remap software error1 in bilinear_one_ratio_solution_of_quadratic_equation");
    EXECUTION_REPORT(REPORT_ERROR, -1, coef_B*coef_B-4*coef_A*coef_C > 0, "remap software error2 in bilinear_one_ratio_solution_of_quadratic_equation");
        
    ratio_u1 = (-coef_B + sqrt(coef_B*coef_B-4*coef_A*coef_C))/(2*coef_A);
    ratio_u2 = (-coef_B - sqrt(coef_B*coef_B-4*coef_A*coef_C))/(2*coef_A);

    ratio_u = -1.0;
    if (ratio_u1 >= 0 && ratio_u1 <= 1) {
        ratio_u = ratio_u1;
        ratio_u_false = ratio_u2;
    }
    else {
        ratio_u = ratio_u2;
        ratio_u_false = ratio_u1;
    }

	if (ratio_u < 0.0 && ratio_u > -eps)
		ratio_u = 0.0;
	if (ratio_u > 1.0 && ratio_u < 1.0+eps)
		ratio_u = 1.0;
    EXECUTION_REPORT(REPORT_WARNING, -1, ratio_u >= 0 && ratio_u <= 1, "remap software error3 in bilinear_one_ratio_solution_of_quadratic_equation: %0.18lf   %0.18lf  %0.18lf : (%lf %lf): (%lf %lf)  (%lf %lf)  (%lf %lf)  (%lf %lf)", ratio_u, ratio_u1, ratio_u2, dst_point_coord_values[0], dst_point_coord_values[1], bilinear_box_vertexes_coord1_values[0], bilinear_box_vertexes_coord2_values[0], bilinear_box_vertexes_coord1_values[1], bilinear_box_vertexes_coord2_values[1], bilinear_box_vertexes_coord1_values[2], bilinear_box_vertexes_coord2_values[2], bilinear_box_vertexes_coord1_values[3], bilinear_box_vertexes_coord2_values[3]);
    //EXECUTION_REPORT(REPORT_ERROR, -1, ratio_u_false < eps || ratio_u_false > 1-eps, "remap software error4 in bilinear_one_ratio_of_solution_quadratic_equation");
}


void Remap_operator_bilinear::bilinear_ratios_solution2(double *dst_point_coord_values, 
                                                        double *bilinear_box_vertexes_coord1_values, 
                                                        double *bilinear_box_vertexes_coord2_values, 
                                                        double &ratio_u, 
                                                        double &ratio_v)
{
    double bilinear_box_vertexes_coord1_values_for_ratio_v[4], bilinear_box_vertexes_coord2_values_for_ratio_v[4];
    
    /*
          (x0,y0)                  (x1,y1)
             p0--P5------p1
             |    |         |
            P6 .--.P4(x,y)-- |
             |    |         |
             |    |         |
             p3----------p2
          (x3,y3)                  (x2,y2)

            P4 is dst point and P0~P3 are 4 vertexes of bilinear box. P0P3 does not parallel to P1P2 and P0P1 does not parallel to P3P2. 
            ratio_u=|P5P0|/|P1P0|, ratio_v=|P6P0|/|P3P0|
    */

    bilinear_one_ratio_solution_of_quadratic_equation(dst_point_coord_values, 
                                                      bilinear_box_vertexes_coord1_values,
                                                      bilinear_box_vertexes_coord2_values,
                                                      ratio_u);

    for (int i = 0; i < 4; i ++) {
        bilinear_box_vertexes_coord1_values_for_ratio_v[i] = bilinear_box_vertexes_coord1_values[(i+1)%4];
        bilinear_box_vertexes_coord2_values_for_ratio_v[i] = bilinear_box_vertexes_coord2_values[(i+1)%4];
    }
    bilinear_one_ratio_solution_of_quadratic_equation(dst_point_coord_values, 
                                                      bilinear_box_vertexes_coord1_values_for_ratio_v,
                                                      bilinear_box_vertexes_coord2_values_for_ratio_v,
                                                      ratio_v);    
}


void Remap_operator_bilinear::solve_two_bilinear_ratios(long *bilinear_box_vertexes_src_cell_indexes, 
                                                        double *dst_point_coord_values,
                                                        double &ratio_u, 
                                                        double &ratio_v)
{
    double bilinear_box_vertexes_coord1_values[4], bilinear_box_vertexes_coord2_values[4];
    double src_cell_center_values[2], cross_product_counter_lines;
    long temp_src_cell_index;
    double eps = 1.0e-12;
    int i;


    for (i = 0; i < 4; i ++) {
        get_cell_center_coord_values_of_src_grid(bilinear_box_vertexes_src_cell_indexes[i], src_cell_center_values);
        bilinear_box_vertexes_coord1_values[i] = src_cell_center_values[0];
        bilinear_box_vertexes_coord2_values[i] = src_cell_center_values[1];
    }
    
    cross_product_counter_lines = compute_cross_product_of_counter_lines(bilinear_box_vertexes_coord1_values, bilinear_box_vertexes_coord2_values);
    if (fabs(cross_product_counter_lines) < eps) {
        bilinear_ratios_solution1(dst_point_coord_values,
                                  bilinear_box_vertexes_coord1_values,
                                  bilinear_box_vertexes_coord2_values,
                                  ratio_u, ratio_v);
        return;
    }

    temp_src_cell_index = bilinear_box_vertexes_src_cell_indexes[0];
    for (i = 1; i < 4; i ++)
        bilinear_box_vertexes_src_cell_indexes[i-1] = bilinear_box_vertexes_src_cell_indexes[i];
    bilinear_box_vertexes_src_cell_indexes[3] = temp_src_cell_index;
    for (i = 0; i < 4; i ++) {
        get_cell_center_coord_values_of_src_grid(bilinear_box_vertexes_src_cell_indexes[i], src_cell_center_values);
        bilinear_box_vertexes_coord1_values[i] = src_cell_center_values[0];
        bilinear_box_vertexes_coord2_values[i] = src_cell_center_values[1];
    }    
    cross_product_counter_lines = compute_cross_product_of_counter_lines(bilinear_box_vertexes_coord1_values, bilinear_box_vertexes_coord2_values);
    if (fabs(cross_product_counter_lines) < eps) {
        bilinear_ratios_solution1(dst_point_coord_values,
                                  bilinear_box_vertexes_coord1_values,
                                  bilinear_box_vertexes_coord2_values,
                                  ratio_u, ratio_v);
    }
    else {
        bilinear_ratios_solution2(dst_point_coord_values,
                                  bilinear_box_vertexes_coord1_values,
                                  bilinear_box_vertexes_coord2_values,
                                  ratio_u, ratio_v);
    }
}


void Remap_operator_bilinear::do_remap_values_caculation(double *data_values_src, double *data_values_dst, int dst_array_size)
{
    remap_weights_groups[0]->remap_values(data_values_src, data_values_dst, dst_array_size);
}


void Remap_operator_bilinear::do_src_decomp_caculation(long *decomp_map_src, const long *decomp_map_dst)
{
    remap_weights_groups[0]->calc_src_decomp(decomp_map_src, decomp_map_dst);
}


Remap_operator_basis *Remap_operator_bilinear::duplicate_remap_operator(bool fully_copy)
{
    Remap_operator_bilinear *duplicated_remap_operator = new Remap_operator_bilinear();
    copy_remap_operator_basic_data(duplicated_remap_operator, fully_copy);
    duplicated_remap_operator->max_num_found_nearest_points = max_num_found_nearest_points;
    duplicated_remap_operator->num_nearest_points = num_nearest_points;
    duplicated_remap_operator->num_power = num_power;
    duplicated_remap_operator->iterative_threshold_distance = iterative_threshold_distance;
    return duplicated_remap_operator;
}


Remap_operator_basis *Remap_operator_bilinear::generate_parallel_remap_operator(Remap_grid_class **decomp_original_grids, int **global_cells_local_indexes_in_decomps)
{
    Remap_operator_basis *parallel_remap_operator = this->duplicate_remap_operator(false);
    this->generate_parallel_remap_weights(parallel_remap_operator, decomp_original_grids, global_cells_local_indexes_in_decomps);
    return parallel_remap_operator;
}
