/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#include "cor_global_data.h"
#include "remap_operator_c_interface.h"
#include "remap_grid_class.h"
#include "quick_sort.h"
#include "remap_common_utils.h"
#include "remap_utils_nearest_points.h"
#include "grid_cell_search.h"
#include <math.h>


bool have_fetched_dst_grid_cell_coord_values;
bool using_rotated_grid_data;
long last_dst_cell_index;


void get_cell_mask_of_grid(Remap_operator_grid *grid, long cell_index, bool *mask_value)
{
    EXECUTION_REPORT(REPORT_ERROR, -1, cell_index >= 0 && cell_index < grid->get_grid_size(),
                 "remap software error in get_cell_mask_of_src_grid\n");

    if (grid->get_mask_values() == NULL)
        *mask_value = true;
    else *mask_value = grid->get_mask_values()[cell_index];    
}


void get_cell_mask_of_src_grid(long cell_index, bool *mask_value)
{
    get_cell_mask_of_grid(current_runtime_remap_operator_grid_src, cell_index, mask_value);
}


long get_size_of_src_grid()
{
    return current_runtime_remap_operator_grid_src->get_grid_size();
}


long get_size_of_dst_grid()
{
    return current_runtime_remap_operator_grid_dst->get_grid_size();
}


void get_cell_mask_of_dst_grid(long cell_index, bool *mask_value)
{
    get_cell_mask_of_grid(current_runtime_remap_operator_grid_dst, cell_index, mask_value);
}


void get_cell_center_coord_values_of_grid(Remap_operator_grid *grid, long cell_index, double *center_values)
{
    EXECUTION_REPORT(REPORT_ERROR, -1, cell_index >= 0 && cell_index < grid->get_grid_size(),
                 "remap software error in get_cell_center_coord_values_of_grid\n");

    for (int i = 0; i < grid->get_num_grid_dimensions(); i ++)
        center_values[i] = grid->get_center_coord_values()[i][cell_index];
}


void get_cell_center_coord_values_of_src_grid(long cell_index, double *center_values)
{
    EXECUTION_REPORT(REPORT_ERROR, -1, have_fetched_dst_grid_cell_coord_values, "remap software error in get_cell_center_coord_values_of_src_grid\n");
    if (using_rotated_grid_data)
        get_cell_center_coord_values_of_grid(current_runtime_remap_operator_grid_src->get_rotated_remap_operator_grid(), cell_index, center_values);
    else get_cell_center_coord_values_of_grid(current_runtime_remap_operator_grid_src, cell_index, center_values);
}


void get_cell_center_coord_values_of_dst_grid(long cell_index, double *center_values)
{
    have_fetched_dst_grid_cell_coord_values = true;
    get_cell_center_coord_values_of_grid(current_runtime_remap_operator_grid_dst, cell_index, center_values);
    if (current_runtime_remap_operator_grid_dst->get_rotated_remap_operator_grid() != NULL && fabs(center_values[1]) > SPHERE_GRID_ROTATION_LAT_THRESHOLD) {
        EXECUTION_REPORT(REPORT_ERROR, -1, last_dst_cell_index == -1 || last_dst_cell_index == cell_index, "remap software error in get_cell_center_coord_values_of_dst_grid\n");
        last_dst_cell_index = cell_index;
        using_rotated_grid_data = true;
        get_cell_center_coord_values_of_grid(current_runtime_remap_operator_grid_dst->get_rotated_remap_operator_grid(), cell_index, center_values);        
    }        
}


H2D_grid_cell_search_engine *get_current_grid2D_search_engine(bool is_src_grid)
{
    if (is_src_grid) {
        if (using_rotated_grid_data)
            return current_runtime_remap_operator_grid_src->get_rotated_remap_operator_grid()->get_grid2D_search_engine();
        return current_runtime_remap_operator_grid_src->get_grid2D_search_engine();
    }
    else {
        if (using_rotated_grid_data)
            return current_runtime_remap_operator_grid_dst->get_rotated_remap_operator_grid()->get_grid2D_search_engine();
        return current_runtime_remap_operator_grid_dst->get_grid2D_search_engine();        
    }
}


void get_cell_vertex_coord_values_of_grid(Remap_operator_grid *grid, long cell_index, int *num_vertex, double *vertex_values, bool check_consistency)
{
    int i, j, tmp_num_dimensions;

    
    EXECUTION_REPORT(REPORT_ERROR, -1, cell_index >= 0 && cell_index < grid->get_grid_size(),
                 "remap software error1 in get_cell_vertex_coord_values_of_grid\n");
    if (check_consistency)
        EXECUTION_REPORT(REPORT_ERROR, -1, have_fetched_dst_grid_cell_coord_values, "remap software error2 in get_cell_vertex_coord_values_of_grid\n");
    
    tmp_num_dimensions = grid->get_num_grid_dimensions();
    *num_vertex = 0;
    for (i = 0; i < grid->get_num_vertexes(); i ++) {
        if (grid->get_vertex_coord_values()[0][grid->get_num_vertexes()*cell_index+i] == NULL_COORD_VALUE)
            continue;
        for (j = 0; j < tmp_num_dimensions; j ++)
            vertex_values[(*num_vertex)*tmp_num_dimensions+j] = grid->get_vertex_coord_values()[j][grid->get_num_vertexes()*cell_index+i];
        (*num_vertex) ++;
    }
}


void get_cell_vertex_coord_values_of_src_grid(long cell_index, int *num_vertex, double *vertex_values, bool check_consistency)
{
    int temp_num_vertex;
    double temp_vertex_values[65536];
    bool should_rotate;

    
    if (check_consistency) {
        if (using_rotated_grid_data)
            get_cell_vertex_coord_values_of_grid(current_runtime_remap_operator_grid_src->get_rotated_remap_operator_grid(), cell_index, num_vertex, vertex_values, check_consistency);
        else get_cell_vertex_coord_values_of_grid(current_runtime_remap_operator_grid_src, cell_index, num_vertex, vertex_values, check_consistency);
    }
    else {
        get_cell_vertex_coord_values_of_grid(current_runtime_remap_operator_grid_src, cell_index, &temp_num_vertex, temp_vertex_values, check_consistency);
        EXECUTION_REPORT(REPORT_ERROR, -1, temp_num_vertex <= 65536/2, "Software error in get_cell_vertex_coord_values_of_src_grid: too big number of vertexes: %d", temp_num_vertex);
        should_rotate = false;
        for (int i = 0; i < temp_num_vertex; i ++)
            if (fabs(temp_vertex_values[i*2+1]) > SPHERE_GRID_ROTATION_LAT_THRESHOLD) {
                should_rotate = true;
                break;
            }
        if (should_rotate)
            get_cell_vertex_coord_values_of_grid(current_runtime_remap_operator_grid_src->get_rotated_remap_operator_grid(), cell_index, num_vertex, vertex_values, check_consistency);
        else get_cell_vertex_coord_values_of_grid(current_runtime_remap_operator_grid_src, cell_index, num_vertex, vertex_values, check_consistency);
    }
}


void get_cell_vertex_coord_values_of_dst_grid(long cell_index, int *num_vertex, double *vertex_values, bool check_consistency)
{
    int temp_num_vertex;
    double temp_vertex_values[65536];
    bool should_rotate;

    
    if (check_consistency) {
        if (using_rotated_grid_data)
            get_cell_vertex_coord_values_of_grid(current_runtime_remap_operator_grid_dst->get_rotated_remap_operator_grid(), cell_index, num_vertex, vertex_values, check_consistency);
        else get_cell_vertex_coord_values_of_grid(current_runtime_remap_operator_grid_dst, cell_index, num_vertex, vertex_values, check_consistency);
    }
    else {
        get_cell_vertex_coord_values_of_grid(current_runtime_remap_operator_grid_dst, cell_index, &temp_num_vertex, temp_vertex_values, check_consistency);
        EXECUTION_REPORT(REPORT_ERROR, -1, temp_num_vertex <= 65536/2, "Software error in get_cell_vertex_coord_values_of_dst_grid: too big number of vertexes: %d", temp_num_vertex);        
        should_rotate = false;
        for (int i = 0; i < temp_num_vertex; i ++)
            if (fabs(temp_vertex_values[i*2+1]) > SPHERE_GRID_ROTATION_LAT_THRESHOLD) {
                should_rotate = true;
                break;
            }
        if (should_rotate)
            get_cell_vertex_coord_values_of_grid(current_runtime_remap_operator_grid_dst->get_rotated_remap_operator_grid(), cell_index, num_vertex, vertex_values, check_consistency);
        else get_cell_vertex_coord_values_of_grid(current_runtime_remap_operator_grid_dst, cell_index, num_vertex, vertex_values, check_consistency);
    }
}


void search_cell_in_src_grid(double *point_coord_values, long *cell_index, bool accurately_match) 
{   
    EXECUTION_REPORT(REPORT_ERROR, -1, have_fetched_dst_grid_cell_coord_values, "remap software error search_cell_in_src_grid\n");
    if (using_rotated_grid_data)
        *cell_index = current_runtime_remap_operator_grid_src->get_rotated_remap_operator_grid()->search_cell_of_locating_point(point_coord_values, accurately_match);
    else *cell_index = current_runtime_remap_operator_grid_src->search_cell_of_locating_point(point_coord_values, accurately_match);
}


void initialize_computing_remap_weights_of_one_cell()
{
    EXECUTION_REPORT(REPORT_ERROR, -1, current_runtime_remap_operator_grid_src->get_num_visited_cells() == 0,
                 "remap software error in initialize_computing_remap_weights_of_one_cell\n");
    have_fetched_dst_grid_cell_coord_values = false;
    using_rotated_grid_data = false;
    last_dst_cell_index = -1;
}


void finalize_computing_remap_weights_of_one_cell()
{
    using_rotated_grid_data = false;
    have_fetched_dst_grid_cell_coord_values = false;
}


void clear_src_grid_cell_visiting_info()
{
    current_runtime_remap_operator_grid_src->clear_cell_visiting_info();
}


void clear_remap_weight_info_in_sparse_matrix()
{
    for (int i = 0; i < current_runtime_remap_operator->get_num_remap_weights_groups(); i ++)
    current_runtime_remap_operator->get_remap_weights_group(i)->clear_weights_info();
}


void add_remap_weights_to_sparse_matrix(long *indexes_src_grid, long index_dst_grid, double *weight_values, int num_weights, int weights_group_index, bool is_real_weight)
{
    current_runtime_remap_operator->get_remap_weights_group(weights_group_index)->add_weights(indexes_src_grid, index_dst_grid, weight_values, num_weights, is_real_weight);
}


double compute_difference_of_two_coord_values(double coord_value1, double coord_value2, int dim_id)
{
    if (is_coord_unit_degree[dim_id])
        match_degree_values(coord_value1, coord_value2);

    return coord_value1 - coord_value2;
}


void sort_polygon_vertexes(double point_coord_value_dim1, 
                           double point_coord_value_dim2,
                           double *polygon_vertexes_coord_values_dim1, 
                           double *polygon_vertexes_coord_values_dim2, 
                           long *polygon_vertexes_src_cell_indexes,
                           int num_polygon_vertexes)
{
    double diff1, diff2, distance;
    double normalized_polygon_vertexes_coord_values_group1[4096], normalized_polygon_vertexes_coord_values_group2[4096];
    double tmp_polygon_vertexes_coord_values_dim1[4096], tmp_polygon_vertexes_coord_values_dim2[4096];
    long polygon_vertexes_src_cell_indexes_group1[4096], polygon_vertexes_src_cell_indexes_group2[4096];
    long tmp_polygon_vertexes_src_cell_indexes[4096];
    int i, j, num_vertexes_group1, num_vertexes_group2;

    
    EXECUTION_REPORT(REPORT_ERROR, -1, num_polygon_vertexes < 4096, "remap software error1 in sort_polygon_vertexes\n");

    /* Sort the vertexes of polygon in the order of clockwise */
    for (i = 0, num_vertexes_group1 = 0, num_vertexes_group2 = 0; i < num_polygon_vertexes; i ++) {
        diff1 = compute_difference_of_two_coord_values(polygon_vertexes_coord_values_dim1[i], point_coord_value_dim1, 0);
        diff2 = compute_difference_of_two_coord_values(polygon_vertexes_coord_values_dim2[i], point_coord_value_dim2, 1);
        if (diff2 == 0)
            if (diff1 > 0)
                diff1 = 1;
            else diff1 = -1;
        else if (diff1 == 0)
            if (diff2 > 0)
                diff2 = 1;
            else diff2 = -1;
        else {
            distance = sqrt(diff1*diff1 + diff2*diff2);
            diff1 = diff1 / distance;
            diff2 = diff2 / distance;
            EXECUTION_REPORT(REPORT_ERROR, -1, distance > 0 && diff1 != 0 && diff2 != 0, "remap software error2 in sort_polygon_vertexes\n");
        }
        if (diff2 > 0 || (diff2 == 0 && diff1 > 0)) {
            normalized_polygon_vertexes_coord_values_group1[num_vertexes_group1] = diff1;
            polygon_vertexes_src_cell_indexes_group1[num_vertexes_group1] = i;
            num_vertexes_group1 ++;
        }
        else {
            normalized_polygon_vertexes_coord_values_group2[num_vertexes_group2] = diff1;
            polygon_vertexes_src_cell_indexes_group2[num_vertexes_group2] = i;
            num_vertexes_group2 ++;
        }
    }
    
    do_quick_sort(normalized_polygon_vertexes_coord_values_group1, polygon_vertexes_src_cell_indexes_group1, 0, num_vertexes_group1-1);
    do_quick_sort(normalized_polygon_vertexes_coord_values_group2, polygon_vertexes_src_cell_indexes_group2, 0, num_vertexes_group2-1);
    
    for (i = 0, j = 0; i < num_vertexes_group1; i ++) {
        tmp_polygon_vertexes_coord_values_dim1[j] = polygon_vertexes_coord_values_dim1[polygon_vertexes_src_cell_indexes_group1[i]];
        tmp_polygon_vertexes_coord_values_dim2[j] = polygon_vertexes_coord_values_dim2[polygon_vertexes_src_cell_indexes_group1[i]];
        tmp_polygon_vertexes_src_cell_indexes[j] = polygon_vertexes_src_cell_indexes[polygon_vertexes_src_cell_indexes_group1[i]];
        j++;
    }
    for (i = num_vertexes_group2-1; i >= 0; i --) {
        tmp_polygon_vertexes_coord_values_dim1[j] = polygon_vertexes_coord_values_dim1[polygon_vertexes_src_cell_indexes_group2[i]];
        tmp_polygon_vertexes_coord_values_dim2[j] = polygon_vertexes_coord_values_dim2[polygon_vertexes_src_cell_indexes_group2[i]];
        tmp_polygon_vertexes_src_cell_indexes[j] = polygon_vertexes_src_cell_indexes[polygon_vertexes_src_cell_indexes_group2[i]];        
        j ++;
    }
    for (i = 0; i < num_polygon_vertexes; i ++) {
        polygon_vertexes_coord_values_dim1[i] = tmp_polygon_vertexes_coord_values_dim1[i];
        polygon_vertexes_coord_values_dim2[i] = tmp_polygon_vertexes_coord_values_dim2[i];
        polygon_vertexes_src_cell_indexes[i] = tmp_polygon_vertexes_src_cell_indexes[i];
    }
}


double compute_vector_angle(double *base_point_coord_values, double *end_point_coord_values)
{
    double coord_values_diff[2], length;


    coord_values_diff[0] = compute_difference_of_two_coord_values(end_point_coord_values[0], base_point_coord_values[0], 0);
    coord_values_diff[1] = compute_difference_of_two_coord_values(end_point_coord_values[1], base_point_coord_values[1], 1);

    if (coord_values_diff[1] == 0)
        return coord_values_diff[0] >= 0 ? 0 : PI;

    length = sqrt(coord_values_diff[0]*coord_values_diff[0] + coord_values_diff[1]*coord_values_diff[1]);
    return coord_values_diff[1] > 0 ? acos(coord_values_diff[0]/length) : acos(coord_values_diff[0]/length)+PI;
}


bool two_bounding_boxes_have_overlap(double *bounding_box1, double *bounding_box2, int num_grid_dimensions)
{
    for (int i = 0; i < num_grid_dimensions; i ++) {
        if (bounding_box1[i*2] > bounding_box1[i*2+1] || bounding_box2[i*2] > bounding_box2[i*2+1])
            for (int j = 0; j < 2; j ++) {
                if (bounding_box1[i*2+j] < 180)
                    bounding_box1[i*2+j] += 360;
                if (bounding_box2[i*2+j] < 180)
                    bounding_box2[i*2+j] += 360;
            }
        if (bounding_box1[i*2+1] < bounding_box2[i*2] || bounding_box2[i*2+1] < bounding_box1[i*2])
            return false;
    }

    return true;
}


void compute_cell_bounding_box(int num_vertexes, int num_grid_dimensions, double *vertex_coord_values, double *cell_bounding_box)
{
    double min_bounding_value, max_bounding_value, temp_coord_value;
    int i, j;

    
    for (i = 0; i < num_grid_dimensions; i ++) {
        min_bounding_value = DEFAULT_FILL_VALUE;
        max_bounding_value = -DEFAULT_FILL_VALUE;
        for (j = 0; j < num_vertexes; j ++) {
            if (vertex_coord_values[num_grid_dimensions*j+i] == NULL_COORD_VALUE)
                continue;
            temp_coord_value = vertex_coord_values[num_grid_dimensions*j+i];
            if (min_bounding_value > temp_coord_value)
                min_bounding_value = temp_coord_value;
            if (max_bounding_value < temp_coord_value)
                max_bounding_value = temp_coord_value;            
        }
        cell_bounding_box[i*2] = min_bounding_value;
        cell_bounding_box[i*2+1] = max_bounding_value;
        if (is_coord_unit_degree[i] && max_bounding_value - min_bounding_value >= 180) {
            min_bounding_value = DEFAULT_FILL_VALUE;
            max_bounding_value = -DEFAULT_FILL_VALUE;
            for (j = 0; j < num_vertexes; j ++) {
                if (vertex_coord_values[num_grid_dimensions*j+i] == NULL_COORD_VALUE)
                    continue;
                temp_coord_value = vertex_coord_values[num_grid_dimensions*j+i];
                if (temp_coord_value >= 180)
                    temp_coord_value -= 360;
                if (min_bounding_value > temp_coord_value)
                    min_bounding_value = temp_coord_value;
                if (max_bounding_value < temp_coord_value)
                    max_bounding_value = temp_coord_value;            
            }
            cell_bounding_box[i*2] = min_bounding_value+360;
            cell_bounding_box[i*2+1] = max_bounding_value;        
        }
    }
}


bool do_two_cells_bounding_box_have_overlap(int num_vertexes_src, 
                                            int num_vertexes_dst,
                                            int num_grid_dimensions,
                                            double *vertex_coord_values_src,
                                            double *vertex_coord_values_dst)
{
    double cell_bounding_box_src[256*2], cell_bounding_box_dst[256*2];

    
    compute_cell_bounding_box(num_vertexes_src, num_grid_dimensions, vertex_coord_values_src, cell_bounding_box_src);
    compute_cell_bounding_box(num_vertexes_dst, num_grid_dimensions, vertex_coord_values_dst, cell_bounding_box_dst);
    return two_bounding_boxes_have_overlap(cell_bounding_box_src, cell_bounding_box_dst, num_grid_dimensions);
}


bool src_cell_and_dst_cell_have_overlap(long cell_index_src, long cell_index_dst)
{
    int num_vertexes_src, num_vertexes_dst, num_grid_dimensions;
    double vertex_coord_values_src[65536], vertex_coord_values_dst[65536];
    double center_coord_values_src[256], center_coord_values_dst[256];


    get_cell_center_coord_values_of_dst_grid(cell_index_dst, center_coord_values_dst);
    get_cell_center_coord_values_of_src_grid(cell_index_src, center_coord_values_src);
    get_cell_vertex_coord_values_of_dst_grid(cell_index_dst, &num_vertexes_dst, vertex_coord_values_dst, true);
    get_cell_vertex_coord_values_of_src_grid(cell_index_src, &num_vertexes_src, vertex_coord_values_src, true);    
    EXECUTION_REPORT(REPORT_ERROR, -1, num_vertexes_src <= 65536/2, "Software error in src_cell_and_dst_cell_have_overlap: too big number of src vertexes: %d", num_vertexes_src);
    EXECUTION_REPORT(REPORT_ERROR, -1, num_vertexes_dst <= 65536/2, "Software error in src_cell_and_dst_cell_have_overlap: too big number of dst vertexes: %d", num_vertexes_dst);

    num_grid_dimensions = current_runtime_remap_operator->get_num_dimensions();

    return do_two_cells_bounding_box_have_overlap(num_vertexes_src, 
                                                  num_vertexes_dst,
                                                  num_grid_dimensions,
                                                  vertex_coord_values_src,
                                                  vertex_coord_values_dst);
}


bool have_overlapped_src_cells_for_dst_cell(long cell_index_dst)
{
    int num_overlapping_cells;
    long overlapping_cells_index[2];

    
    const H2D_grid_cell_search_cell *dst_cell = get_current_grid2D_search_engine(false)->get_cell(cell_index_dst);
    get_current_grid2D_search_engine(true)->search_overlapping_cells(num_overlapping_cells, overlapping_cells_index, dst_cell, true, true);
    EXECUTION_REPORT(REPORT_ERROR, -1, num_overlapping_cells == 0 || num_overlapping_cells == 1, "Software error1 in have_overlapped_src_cells_for_dst_cell");

    return num_overlapping_cells > 0;
}


bool is_point_on_great_arc(double lon_point,
                           double lat_point,
                           double lon_arc_start_point,
                           double lat_arc_start_point,
                           double lon_arc_end_point,
                           double lat_arc_end_point)
{
    double dist1, dist2, dist3;


    dist1 = calculate_distance_of_two_points_2D(lon_point, lat_point, lon_arc_start_point, lat_arc_start_point, true);
    dist2 = calculate_distance_of_two_points_2D(lon_point, lat_point, lon_arc_end_point, lat_arc_end_point, true);
    dist3 = calculate_distance_of_two_points_2D(lon_arc_start_point, lat_arc_start_point, lon_arc_end_point, lat_arc_end_point, true);    

    return dist1 <= dist3 && dist2 <= dist3;
}


bool are_the_same_sphere_points(double lon_point1, double lat_point1, double lon_point2, double lat_point2)
{
    double eps = 1.0e-4;
    return (fabs(compute_difference_of_two_coord_values(lon_point1,lon_point2,0)) < eps && fabs(lat_point1-lat_point2) < eps);
}


void compute_normal_vector_3D(double x_vector1, double y_vector1, double z_vector1,
                              double x_vector2, double y_vector2, double z_vector2,
                              double &x_vectorn, double &y_vectorn, double &z_vectorn)
{
    x_vectorn = y_vector1*z_vector2 - z_vector1*y_vector2;
    y_vectorn = -x_vector1*z_vector2 + z_vector1*x_vector2;
    z_vectorn = x_vector1*y_vector2 - y_vector1*x_vector2;    
}


void compute_intersect_points_of_two_great_arcs_of_sphere_grid(double lon_arc1_point1,
                                                                double lat_arc1_point1,
                                                                double lon_arc1_point2,
                                                                double lat_arc1_point2,
                                                                double lon_arc2_point1,
                                                                double lat_arc2_point1,
                                                                double lon_arc2_point2,
                                                                double lat_arc2_point2,
                                                                int &num_intersect_points,
                                                                double *lon_intersect_points,
                                                                double *lat_intersect_points)
{
    double eps = 1.0e-10;
    double x_arc1_point1, y_arc1_point1, z_arc1_point1;
    double x_arc1_point2, y_arc1_point2, z_arc1_point2;
    double x_arc2_point1, y_arc2_point1, z_arc2_point1;
    double x_arc2_point2, y_arc2_point2, z_arc2_point2;
    double lon_candidate_intersect_points_in_radian[256], lat_candidate_intersect_points_in_radian[256];
    double lon_candidate_intersect_points_in_degree[256], lat_candidate_intersect_points_in_degree[256];
    double a, b, c, d, e, f;
    long double g1, g2, g3, r;
    long double x, y, z;
    int i, j;


    EXECUTION_REPORT(REPORT_ERROR, -1, !(lon_arc1_point1 == lon_arc1_point2 && lat_arc1_point1 == lat_arc1_point2),
                 "remap software error1 in compute_intersect_points_of_two_great_arcs_of_sphere_grid\n");
    EXECUTION_REPORT(REPORT_ERROR, -1, !(lon_arc2_point1 == lon_arc2_point2 && lat_arc2_point1 == lat_arc2_point2),
                 "remap software error2 in compute_intersect_points_of_two_great_arcs_of_sphere_grid\n");

    num_intersect_points = 0;

    get_3D_cartesian_coord_of_sphere_coord(x_arc1_point1, y_arc1_point1, z_arc1_point1, lon_arc1_point1, lat_arc1_point1);
    get_3D_cartesian_coord_of_sphere_coord(x_arc1_point2, y_arc1_point2, z_arc1_point2, lon_arc1_point2, lat_arc1_point2);
    get_3D_cartesian_coord_of_sphere_coord(x_arc2_point1, y_arc2_point1, z_arc2_point1, lon_arc2_point1, lat_arc2_point1);
    get_3D_cartesian_coord_of_sphere_coord(x_arc2_point2, y_arc2_point2, z_arc2_point2, lon_arc2_point2, lat_arc2_point2);

    compute_normal_vector_3D(x_arc1_point1, y_arc1_point1, z_arc1_point1,
                             x_arc1_point2, y_arc1_point2, z_arc1_point2,
                             a, b, c);
    compute_normal_vector_3D(x_arc2_point1, y_arc2_point1, z_arc2_point1,
                             x_arc2_point2, y_arc2_point2, z_arc2_point2,
                             d, e, f);    
    
    g1 = c*e - b*f;
    g2 = b*d - a*e;
    g3 = a*f - c*d;
    r = sqrtl(g1*g1+g2*g2+g3*g3);

    if (r < eps) {
        if (is_point_on_great_arc(lon_arc1_point1, lat_arc1_point1,
                                  lon_arc2_point1, lat_arc2_point1,
                                  lon_arc2_point2, lat_arc2_point2)) {
                lon_intersect_points[num_intersect_points] = lon_arc1_point1;
                lat_intersect_points[num_intersect_points] = lat_arc1_point1;
                num_intersect_points ++;
        }
        if (is_point_on_great_arc(lon_arc1_point2, lat_arc1_point2,
                                  lon_arc2_point1, lat_arc2_point1,
                                  lon_arc2_point2, lat_arc2_point2)) {
            lon_intersect_points[num_intersect_points] = lon_arc1_point2;
            lat_intersect_points[num_intersect_points] = lat_arc1_point2;
            num_intersect_points ++;
        }
        if (is_point_on_great_arc(lon_arc2_point1, lat_arc2_point1,
                                  lon_arc1_point1, lat_arc1_point1,
                                  lon_arc1_point2, lat_arc1_point2)) {
            lon_intersect_points[num_intersect_points] = lon_arc2_point1;
            lat_intersect_points[num_intersect_points] = lat_arc2_point1;
            num_intersect_points ++;
        }
        if (is_point_on_great_arc(lon_arc2_point2, lat_arc2_point2,
                                  lon_arc1_point1, lat_arc1_point1,
                                  lon_arc1_point2, lat_arc1_point2)) {
            lon_intersect_points[num_intersect_points] = lon_arc2_point2;
            lat_intersect_points[num_intersect_points] = lat_arc2_point2;
            num_intersect_points ++;
        }
        if (num_intersect_points > 2) {
            for (i = 0; i < num_intersect_points; i ++) {
                for (j = i+1; j < num_intersect_points; j ++) {
                    if (are_the_same_sphere_points(lon_intersect_points[i],
                                                   lat_intersect_points[i],
                                                   lon_intersect_points[j],
                                                   lat_intersect_points[j]))
                        break;
                }
                if (j < num_intersect_points)
                    num_intersect_points --;
                for (; j < num_intersect_points; j ++) {
                    lon_intersect_points[j] = lon_intersect_points[j+1];
                    lat_intersect_points[j] = lat_intersect_points[j+1];
                }
            }
        }
        EXECUTION_REPORT(REPORT_ERROR, -1, num_intersect_points <= 2, "remap software error3 in compute_intersect_points_of_two_great_arcs_of_sphere_grid\n");
        return;
    }

    x = g1 / r;
    z = g2 / r;
    y = g3 / r;
    lon_candidate_intersect_points_in_radian[0] = atan2l(y, x);
    lon_candidate_intersect_points_in_radian[1] = atan2l(-y, -x);
    lat_candidate_intersect_points_in_radian[0] = asinl(z);
    lat_candidate_intersect_points_in_radian[1] = -lat_candidate_intersect_points_in_radian[0];
    for (i = 0; i < 2; i ++) {
        if (lon_candidate_intersect_points_in_radian[i] < 0.0) 
            lon_candidate_intersect_points_in_radian[i] += PI*2;
        if (lon_candidate_intersect_points_in_radian[i] > PI*2) 
            lon_candidate_intersect_points_in_radian[i] -= PI*2;
        lon_candidate_intersect_points_in_degree[i] = RADIAN_TO_DEGREE(lon_candidate_intersect_points_in_radian[i]);
        lat_candidate_intersect_points_in_degree[i] = RADIAN_TO_DEGREE(lat_candidate_intersect_points_in_radian[i]);
        if (is_point_on_great_arc(lon_candidate_intersect_points_in_degree[i], lat_candidate_intersect_points_in_degree[i], 
                                  lon_arc1_point1, lat_arc1_point1,
                                  lon_arc1_point2, lat_arc1_point2) &&
            is_point_on_great_arc(lon_candidate_intersect_points_in_degree[i], lat_candidate_intersect_points_in_degree[i], 
                                  lon_arc2_point1, lat_arc2_point1,
                                  lon_arc2_point2, lat_arc2_point2)) {
            lon_intersect_points[num_intersect_points] = lon_candidate_intersect_points_in_degree[i];
            lat_intersect_points[num_intersect_points] = lat_candidate_intersect_points_in_degree[i];
            num_intersect_points ++;        
        }
    }

    EXECUTION_REPORT(REPORT_ERROR, -1, num_intersect_points <= 1, "remap software error5 in compute_intersect_points_of_two_great_arcs_of_sphere_grid\n");
}


void get_all_vertexes_of_one_cell_in_other_cell(int num_vertexes_cell1,
                                                int num_vertexes_cell2,
                                                double *vertex_coord_values_cell1,
                                                double *vertex_coord_values_cell2,
                                                int &num_vertexes_in_other_cell,
                                                double *coord1_values_vertexes_in_other_cell,
                                                double *coord2_values_vertexes_in_other_cell)
{
    double temp_vertex_coord1_values[65536], temp_vertex_coord2_values[65536];
    int i;


    EXECUTION_REPORT(REPORT_ERROR, -1, num_vertexes_cell2 <= 65536/2, "Software error in get_all_vertexes_of_one_cell_in_other_cell: too big number of num_vertexes_cell2: %d", num_vertexes_cell2);
    num_vertexes_in_other_cell = 0;
    for (i = 0; i < num_vertexes_cell2; i ++) {
        temp_vertex_coord1_values[i] = vertex_coord_values_cell2[i*2];
        temp_vertex_coord2_values[i] = vertex_coord_values_cell2[i*2+1];
    }
    for (i = 0; i < num_vertexes_cell1; i ++)
        if (is_point_in_2D_cell(vertex_coord_values_cell1[i*2], 
                                vertex_coord_values_cell1[i*2+1],
                                temp_vertex_coord1_values,
                                temp_vertex_coord2_values,
                                num_vertexes_cell2,
                                is_coord_unit_degree[0], 
                                is_coord_unit_degree[1], 
                                current_runtime_remap_operator->get_is_sphere_grid())) {
            coord1_values_vertexes_in_other_cell[num_vertexes_in_other_cell] = vertex_coord_values_cell1[i*2];
            coord2_values_vertexes_in_other_cell[num_vertexes_in_other_cell] = vertex_coord_values_cell1[i*2+1];
            num_vertexes_in_other_cell ++;
        }
}


void compute_arc_points_within_sphere_cell(double lon_arc_start,
                                           double lat_arc_start,
                                           double lon_arc_end,
                                           double lat_arc_end,
                                           int num_cell_vertexes,
                                           double *cell_vertexes_lons,
                                           double *cell_vertexes_lats,
                                           int &num_arc_points_within_cell,
                                           double *lons_arc_points_within_cell,
                                           double *lats_arc_points_within_cell)
{
    bool enclose_arc_end_point, point_has_been_recorded;
    int i, next_i, j, k, num_intersect_points;
    double lon_intersect_points[256], lat_intersect_points[256];

    
    num_arc_points_within_cell = 0;
    if (is_point_in_2D_cell(lon_arc_start, 
                            lat_arc_start,
                            cell_vertexes_lons,
                            cell_vertexes_lats, 
                            num_cell_vertexes,
                            true, true, true)) {
        lons_arc_points_within_cell[num_arc_points_within_cell] = lon_arc_start;
        lats_arc_points_within_cell[num_arc_points_within_cell] = lat_arc_start;
        num_arc_points_within_cell ++;
    }

    enclose_arc_end_point = is_point_in_2D_cell(lon_arc_end, 
                                                lat_arc_end,
                                                cell_vertexes_lons,
                                                cell_vertexes_lats, 
                                                num_cell_vertexes,
                                                true, true, true);

    if (num_arc_points_within_cell == 1 && enclose_arc_end_point) {
        lons_arc_points_within_cell[num_arc_points_within_cell] = lon_arc_end;
        lats_arc_points_within_cell[num_arc_points_within_cell] = lat_arc_end;
        num_arc_points_within_cell ++;        
        return;
    }

    for (i = 0; i < num_cell_vertexes; i ++) {
        if (cell_vertexes_lons[i] == NULL_COORD_VALUE)
            continue;
        next_i = (i+1) % num_cell_vertexes;
        while (cell_vertexes_lons[next_i] == NULL_COORD_VALUE)
            next_i = (next_i+1) % num_cell_vertexes;
        compute_intersect_points_of_two_great_arcs_of_sphere_grid(lon_arc_start,
                                                                  lat_arc_start,
                                                                  lon_arc_end,
                                                                  lat_arc_end,
                                                                  cell_vertexes_lons[i],
                                                                  cell_vertexes_lats[i],
                                                                  cell_vertexes_lons[next_i],
                                                                  cell_vertexes_lats[next_i],
                                                                  num_intersect_points,
                                                                  lon_intersect_points,
                                                                  lat_intersect_points);
        for (j = 0; j < num_intersect_points; j ++) {
            point_has_been_recorded = false;
            for (k = 0; k < num_arc_points_within_cell; k ++)
                if (are_the_same_sphere_points(lon_intersect_points[j],
                                               lat_intersect_points[j],
                                               lons_arc_points_within_cell[k],
                                               lats_arc_points_within_cell[k])) {
                    point_has_been_recorded = true;
                    break;
                }
            if (point_has_been_recorded)
                continue;
            lons_arc_points_within_cell[num_arc_points_within_cell] = lon_intersect_points[j];
            lats_arc_points_within_cell[num_arc_points_within_cell] = lat_intersect_points[j];
            num_arc_points_within_cell ++;
            EXECUTION_REPORT(REPORT_ERROR, -1, num_arc_points_within_cell <= 2, "remap software error1 in compute_arc_points_within_sphere_cell\n");            
        }
    }

    if (enclose_arc_end_point) {
        point_has_been_recorded = false;
        if (num_arc_points_within_cell > 0)
            for (k = 0; k < num_arc_points_within_cell; k ++)
                if (are_the_same_sphere_points(lon_arc_end,
                                               lat_arc_end,
                                               lons_arc_points_within_cell[k],
                                               lats_arc_points_within_cell[k])) {
                    lons_arc_points_within_cell[k] = lon_arc_end;
                    lats_arc_points_within_cell[k] = lat_arc_end;
                    point_has_been_recorded = true;
                    break;
                }
        if (!point_has_been_recorded) {
            lons_arc_points_within_cell[num_arc_points_within_cell] = lon_arc_end;
            lats_arc_points_within_cell[num_arc_points_within_cell] = lat_arc_end;
            num_arc_points_within_cell ++;        
            EXECUTION_REPORT(REPORT_ERROR, -1, num_arc_points_within_cell <= 2, "remap software error2 in compute_arc_points_within_sphere_cell\n");
        }
    }
}


double compute_dot_product_of_3D_vectors(double x1, double y1, double z1,
                                         double x2, double y2, double z2)
{
    return x1*x2 + y1*y2 + z1*z2;
}


double compute_angle_of_great_arcs(double lon_point1, double lat_point1,
                                   double lon_point2, double lat_point2,
                                   double lon_point3, double lat_point3,
                                   bool check_angle)
{
    double x1, y1, z1, x2, y2, z2, x3, y3, z3;
    double vector1[3], vector2[3], vector3[3];
    double tmp;
    double angle;
    

    get_3D_cartesian_coord_of_sphere_coord(x1, y1, z1, lon_point1, lat_point1);
    get_3D_cartesian_coord_of_sphere_coord(x2, y2, z2, lon_point2, lat_point2);
    get_3D_cartesian_coord_of_sphere_coord(x3, y3, z3, lon_point3, lat_point3);
    compute_normal_vector_3D(x1, y1, z1, x2, y2, z2, vector1[0], vector1[1], vector1[2]);
    compute_normal_vector_3D(x2, y2, z2, x3, y3, z3, vector2[0], vector2[1], vector2[2]);
    tmp = compute_dot_product_of_3D_vectors(vector1[0], vector1[1], vector1[2], vector2[0], vector2[1], vector2[2]) /
            sqrtl((vector1[0]*vector1[0]+vector1[1]*vector1[1]+vector1[2]*vector1[2])*(vector2[0]*vector2[0]+vector2[1]*vector2[1]+vector2[2]*vector2[2]));
    tmp = fmax(-1.0, fmin(1.0, tmp));
    angle = acosl(-tmp);
    compute_normal_vector_3D(vector1[0], vector1[1], vector1[2],
                             vector2[0], vector2[1], vector2[2],
                             vector3[0], vector3[1], vector3[2]);
    tmp = compute_dot_product_of_3D_vectors(x2, y2, z2, vector3[0], vector3[1], vector3[2]);
    if (tmp < 0)
        angle = -angle;    

//    printf("angle is %lf: %lf\n", angle, tmp);
    
    if (check_angle)
        EXECUTION_REPORT(REPORT_ERROR, -1, angle >= 0, "remap software error in compute_angle_of_great_arcs\n");

    return angle;
}


void sort_vertexes_of_sphere_cell(int num_vertexes,
                                  double *vertexes_lons,
                                  double *vertexes_lats)
{
    int i, index_array[65536];
    volatile double temp_vertexes_lons[65536], temp_vertexes_lats[65536];
    double average_lon, average_lat, angles[65536];
    bool cross_lon_360;


    EXECUTION_REPORT(REPORT_ERROR, -1, num_vertexes <= 65536, "Software error in sort_vertexes_of_sphere_cell: too big number of num_vertexes: %d", num_vertexes);

    for (i = 0; i < num_vertexes; i ++) 
        EXECUTION_REPORT(REPORT_ERROR, -1, vertexes_lons[i] != NULL_COORD_VALUE, "remap software error in sort_vertexes_of_sphere_cell\n");

    cross_lon_360 = false;
    for (i = 1; i < num_vertexes; i ++)
        if (fabs(vertexes_lons[i]-vertexes_lons[0]) > 180)
            cross_lon_360 = true;
    
    average_lon = 0;
    average_lat = 0;
    for (i = 0; i < num_vertexes; i ++) {
        index_array[i] = i;
        if (cross_lon_360 && vertexes_lons[i] > 180)
            average_lon += vertexes_lons[i]-360;
        else average_lon += vertexes_lons[i];
        average_lat += vertexes_lats[i];
    }

    average_lon = average_lon / num_vertexes;
    average_lat = average_lat / num_vertexes;

    angles[0] = 0;
    for (i = 1; i < num_vertexes; i ++) 
        angles[i] = compute_angle_of_great_arcs(vertexes_lons[0], vertexes_lats[0],
                                                average_lon, average_lat,                                            
                                                vertexes_lons[i], vertexes_lats[i],
                                                false);

    do_quick_sort(angles, index_array, 0, num_vertexes-1);
    for (i = 0; i < num_vertexes; i ++) {
        temp_vertexes_lons[i] = vertexes_lons[index_array[i]];
        temp_vertexes_lats[i] = vertexes_lats[index_array[i]];
    }
    for (i = 0; i < num_vertexes; i ++) {
        vertexes_lons[i] = temp_vertexes_lons[num_vertexes-1-i];
        vertexes_lats[i] = temp_vertexes_lats[num_vertexes-1-i];
    }
}


double compute_area_of_sphere_cell(int num_vertexes,
                                         double *vertexes_lons,
                                         double *vertexes_lats)
{
    int i, next_i, next_next_i, num_true_vertexes;
    double area, angle;
    double eps = 2.0e-14;


    if (num_vertexes == 0)
        return DEFAULT_FILL_VALUE;
    
    for (i = 0, num_true_vertexes = 0; i < num_vertexes; i ++)
        if (vertexes_lons[i] != NULL_COORD_VALUE)
            num_true_vertexes ++;

    area = -(num_true_vertexes-2)*PI;
    for (i = 0; i < num_vertexes; i ++) {
//        printf("okok %d\n", i);
        if (vertexes_lons[i] == NULL_COORD_VALUE)
            continue;
        next_i = (i+1) % num_vertexes;
        while (vertexes_lons[next_i] == NULL_COORD_VALUE)
            next_i = (next_i+1) % num_vertexes;
        next_next_i = (next_i+1) % num_vertexes;
        while (vertexes_lons[next_next_i] == NULL_COORD_VALUE)
            next_next_i = (next_next_i+1) % num_vertexes;
        angle = compute_angle_of_great_arcs(vertexes_lons[i], vertexes_lats[i],
                                            vertexes_lons[next_i], vertexes_lats[next_i],
                                            vertexes_lons[next_next_i], vertexes_lats[next_next_i],
                                            true);
        area += angle;
    }

    if (area < 0 && (fabs(area) < eps || num_vertexes == 3))
        area = -area;
    EXECUTION_REPORT(REPORT_ERROR, -1, area >= 0, "remap software error in compute_area_of_sphere_cell\n");
    return area;
}


void compute_common_sub_cell_of_src_cell_and_dst_cell_2D(long cell_index_src, 
                                                         long cell_index_dst, 
                                                         int &num_sub_cell_vertexes, 
                                                         double *sub_cell_vertexes_lons, 
                                                         double *sub_cell_vertexes_lats)
{
    double vertex_coord_values_src[65536], vertex_coord_values_dst[65536];
    double vertex_lons_src[65536], vertex_lats_src[65536], vertex_lons_dst[65536], vertex_lats_dst[65536];
    int num_vertexes_src, num_vertexes_dst, num_grid_dimensions;
    int num_src_vertexes_in_dst_cell, num_dst_vertexes_in_src_cell;
    int i, j, k, next_i, num_arc_points_within_cell;
    double lons_arc_points_within_cell[65536], lats_arc_points_within_cell[65536];
    double lons_src_vertexes_in_dst_cell[65536], lats_src_vertexes_in_dst_cell[65536];
    double lons_dst_vertexes_in_src_cell[65536], lats_dst_vertexes_in_src_cell[65536];


    num_sub_cell_vertexes = 0;
    get_cell_vertex_coord_values_of_dst_grid(cell_index_dst, &num_vertexes_dst, vertex_coord_values_dst, true);
    get_cell_vertex_coord_values_of_src_grid(cell_index_src, &num_vertexes_src, vertex_coord_values_src, true);

    EXECUTION_REPORT(REPORT_ERROR, -1, num_vertexes_src <= 65536/2, "Software error in compute_common_sub_cell_of_src_cell_and_dst_cell_2D: too big number of num_vertexes_src: %d", num_vertexes_src);
    EXECUTION_REPORT(REPORT_ERROR, -1, num_vertexes_dst <= 65536/2, "Software error in compute_common_sub_cell_of_src_cell_and_dst_cell_2D: too big number of num_vertexes_dst: %d", num_vertexes_dst);

    num_grid_dimensions = current_runtime_remap_operator->get_num_dimensions();

    EXECUTION_REPORT(REPORT_ERROR, -1, num_grid_dimensions == 2, "remap software error1 in compute_common_sub_cell_of_src_cell_and_dst_cell_2D\n");

    get_all_vertexes_of_one_cell_in_other_cell(num_vertexes_src,
                                               num_vertexes_dst, 
                                               vertex_coord_values_src, 
                                               vertex_coord_values_dst,
                                               num_src_vertexes_in_dst_cell,
                                               lons_src_vertexes_in_dst_cell,
                                               lats_src_vertexes_in_dst_cell);
    if (num_src_vertexes_in_dst_cell == num_vertexes_src) {
        num_sub_cell_vertexes = num_vertexes_src;
        for (i = 0; i < num_sub_cell_vertexes; i ++) {
            sub_cell_vertexes_lons[i] = vertex_coord_values_src[i*2];
            sub_cell_vertexes_lats[i] = vertex_coord_values_src[i*2+1];
        }
        sort_vertexes_of_sphere_cell(num_sub_cell_vertexes, sub_cell_vertexes_lons, sub_cell_vertexes_lats);
        return;
    }
    get_all_vertexes_of_one_cell_in_other_cell(num_vertexes_dst,
                                               num_vertexes_src, 
                                               vertex_coord_values_dst, 
                                               vertex_coord_values_src,
                                               num_dst_vertexes_in_src_cell,
                                               lons_dst_vertexes_in_src_cell,
                                               lats_dst_vertexes_in_src_cell);
    if (num_dst_vertexes_in_src_cell == num_vertexes_dst) {
        num_sub_cell_vertexes = num_vertexes_dst;
        for (i = 0; i < num_sub_cell_vertexes; i ++) {
            sub_cell_vertexes_lons[i] = vertex_coord_values_dst[i*2];
            sub_cell_vertexes_lats[i] = vertex_coord_values_dst[i*2+1];
        }
        sort_vertexes_of_sphere_cell(num_sub_cell_vertexes, sub_cell_vertexes_lons, sub_cell_vertexes_lats);
        return;
    }

    if (!src_cell_and_dst_cell_have_overlap(cell_index_src, cell_index_dst))
        return;

    for (i = 0; i < num_vertexes_src; i ++) {
        vertex_lons_src[i] = vertex_coord_values_src[2*i];
        vertex_lats_src[i] = vertex_coord_values_src[2*i+1];
    }
    for (i = 0; i < num_vertexes_dst; i ++) {
        vertex_lons_dst[i] = vertex_coord_values_dst[2*i];
        vertex_lats_dst[i] = vertex_coord_values_dst[2*i+1];
    }

    for (i = 0; i < num_vertexes_src; i ++) {
        if (vertex_lons_src[i] == NULL_COORD_VALUE)
            continue;
        next_i = (i+1) % num_vertexes_src;
        while (vertex_lons_src[next_i] == NULL_COORD_VALUE)
            next_i = (next_i+1)%num_vertexes_src;
        compute_arc_points_within_sphere_cell(vertex_lons_src[i],
                                              vertex_lats_src[i],
                                              vertex_lons_src[next_i],
                                              vertex_lats_src[next_i],
                                              num_vertexes_dst,
                                              vertex_lons_dst,
                                              vertex_lats_dst,
                                              num_arc_points_within_cell,
                                              lons_arc_points_within_cell,
                                              lats_arc_points_within_cell);        
        for (j = 0; j < num_arc_points_within_cell; j ++) {
            for (k = 0; k < num_sub_cell_vertexes; k ++)
                if (are_the_same_sphere_points(sub_cell_vertexes_lons[k], sub_cell_vertexes_lats[k],
                                               lons_arc_points_within_cell[j], lats_arc_points_within_cell[j]))
                    break;
            if (k != num_sub_cell_vertexes)
                continue;
            sub_cell_vertexes_lons[num_sub_cell_vertexes] = lons_arc_points_within_cell[j];
            sub_cell_vertexes_lats[num_sub_cell_vertexes] = lats_arc_points_within_cell[j];
            num_sub_cell_vertexes ++;
        }
    }

    for (i = 0; i < num_dst_vertexes_in_src_cell; i ++) {
        for (j = 0; j < num_sub_cell_vertexes; j ++) {
            if (are_the_same_sphere_points(lons_dst_vertexes_in_src_cell[i],
                                           lats_dst_vertexes_in_src_cell[i],
                                           sub_cell_vertexes_lons[j],
                                           sub_cell_vertexes_lats[j]))
                break;
        }
        if (j == num_sub_cell_vertexes) {
            sub_cell_vertexes_lons[num_sub_cell_vertexes] = lons_dst_vertexes_in_src_cell[i];
            sub_cell_vertexes_lats[num_sub_cell_vertexes] = lats_dst_vertexes_in_src_cell[i];
            num_sub_cell_vertexes ++;
        }
    }

    if (num_sub_cell_vertexes == 0)
        return;
        
    sort_vertexes_of_sphere_cell(num_sub_cell_vertexes, sub_cell_vertexes_lons, sub_cell_vertexes_lats);

    EXECUTION_REPORT(REPORT_ERROR, -1, num_sub_cell_vertexes > 0, "remap software error2 in compute_common_sub_cell_of_src_cell_and_dst_cell_2D\n");

    if (num_sub_cell_vertexes <= 2) {
        num_sub_cell_vertexes = 0;
    }
    
    double temp_vertex_lons[65536], temp_vertex_lats[65536];
    double area1, area2, area3;
    if (num_sub_cell_vertexes > 0) {
        for (i = 0; i < num_vertexes_src; i ++) {
            temp_vertex_lons[i] = vertex_coord_values_src[2*i];
            temp_vertex_lats[i] = vertex_coord_values_src[2*i+1];
        }
        sort_vertexes_of_sphere_cell(num_vertexes_src, temp_vertex_lons, temp_vertex_lats);
        area2 = compute_area_of_sphere_cell(num_vertexes_src, temp_vertex_lons, temp_vertex_lats);
        for (i = 0; i < num_vertexes_dst; i ++) {
            temp_vertex_lons[i] = vertex_coord_values_dst[2*i];
            temp_vertex_lats[i] = vertex_coord_values_dst[2*i+1];
        }
        sort_vertexes_of_sphere_cell(num_vertexes_dst, temp_vertex_lons, temp_vertex_lats);
        area3 = compute_area_of_sphere_cell(num_vertexes_dst, temp_vertex_lons, temp_vertex_lats);
        area1 = compute_area_of_sphere_cell(num_sub_cell_vertexes, sub_cell_vertexes_lons, sub_cell_vertexes_lats);
        if (fabs(area1-area2) > 1.0e-7)
            EXECUTION_REPORT(REPORT_ERROR, -1, area1 <= area2, "remap software error5 in compute_common_sub_cell_of_src_cell_and_dst_cell_2D\n");
        if (fabs(area1-area3) > 1.0e-7)
            EXECUTION_REPORT(REPORT_ERROR, -1, area1 <= area3, "remap software error6 in compute_common_sub_cell_of_src_cell_and_dst_cell_2D\n");        
    }

    /*
    for (i = 0; i < num_vertexes_src; i ++) {
        temp_vertex_lons[i] = vertex_coord_values_src[2*i];
        temp_vertex_lats[i] = vertex_coord_values_src[2*i+1];
    }
    sort_vertexes_of_sphere_cell(num_vertexes_src, temp_vertex_lons, temp_vertex_lats);
    for (i = 0; i < num_sub_cell_vertexes; i ++) 
        if (!is_point_in_2D_cell(sub_cell_vertexes_lons[i],
                                 sub_cell_vertexes_lats[i],
                                 temp_vertex_lons,
                                 temp_vertex_lats,
                                 num_vertexes_src,
                                 true, true, true)) 
            EXECUTION_REPORT(REPORT_ERROR, -1, false, "remap software error3 in compute_common_sub_cell_of_src_cell_and_dst_cell_2D\n");
    for (i = 0; i < num_vertexes_dst; i ++) {
        temp_vertex_lons[i] = vertex_coord_values_dst[2*i];
        temp_vertex_lats[i] = vertex_coord_values_dst[2*i+1];
    }
    sort_vertexes_of_sphere_cell(num_vertexes_dst, temp_vertex_lons, temp_vertex_lats);
    for (i = 0; i < num_sub_cell_vertexes; i ++) 
        if (!is_point_in_2D_cell(sub_cell_vertexes_lons[i],
                                 sub_cell_vertexes_lats[i],
                                 temp_vertex_lons,
                                 temp_vertex_lats,
                                 num_vertexes_dst,
                                 true, true, true)) 
            EXECUTION_REPORT(REPORT_ERROR, -1, false, "remap software error4 in compute_common_sub_cell_of_src_cell_and_dst_cell_2D\n");
    */
}

