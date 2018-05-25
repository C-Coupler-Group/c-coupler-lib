/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#include "remap_common_utils.h"
#include "remap_grid_class.h"
#include "remap_utils_nearest_points.h"
#include "remap_operator_c_interface.h"
#include "execution_report.h"
#include <string.h>
#include <math.h>


void match_degree_values(double &data_value1, double &data_value2)
{
    if (data_value1 - data_value2 > 180) 
        data_value2 += 360;
    else if (data_value1 - data_value2 < -180)
        data_value1 += 360;

    if (data_value1 + data_value2 >= 720) {
        data_value1 -= 360;
        data_value2 -= 360;
    }
}


bool words_are_the_same(const char *word1, const char *word2)
{
    if (word1 == NULL)
        if (word2 == NULL || strlen(word2) == 0)
            return true;
        else return false;

    if (word2 == NULL)
        if (strlen(word1) == 0)
            return true;
        else return false;
    
    return strcmp(word1, word2) == 0;
}


double compute_three_3D_points_cross_product(double center_coord_value_x,
                                             double center_coord_value_y,
                                             double center_coord_value_z,
                                             double vertex1_coord_value_x,
                                             double vertex1_coord_value_y,
                                             double vertex1_coord_value_z,
                                             double vertex2_coord_value_x,
                                             double vertex2_coord_value_y,
                                             double vertex2_coord_value_z)
{
    double vector_x = vertex1_coord_value_y*vertex2_coord_value_z - vertex1_coord_value_z*vertex2_coord_value_y;
    double vector_y = vertex1_coord_value_z*vertex2_coord_value_x - vertex1_coord_value_x*vertex2_coord_value_z;
    double vector_z = vertex1_coord_value_x*vertex2_coord_value_y - vertex1_coord_value_y*vertex2_coord_value_x;
    double module1 = sqrt(center_coord_value_x*center_coord_value_x+center_coord_value_y*center_coord_value_y+center_coord_value_z*center_coord_value_z);
    double module2 = sqrt(vector_x*vector_x+vector_y*vector_y+vector_z*vector_z);
    return (center_coord_value_x*vector_x + center_coord_value_y*vector_y + center_coord_value_z*vector_z) / (module1*module2);
}


double compute_three_2D_points_cross_product(double center_coord1_value,
                                                      double center_coord2_value,
                                                      double vertex1_coord1_value,
                                                      double vertex1_coord2_value,
                                                      double vertex2_coord1_value,
                                                      double vertex2_coord2_value,
                                                      bool is_coord1_unit_degree,
                                                      bool is_coord2_unit_degree)
{
    double vertexes_coord1_value_diff, vertexes_coord2_value_diff;
    double center_vertex_coord1_value_diff, center_vertex_coord2_value_diff;
    
    
    vertexes_coord1_value_diff = vertex2_coord1_value - vertex1_coord1_value;
    vertexes_coord2_value_diff = vertex2_coord2_value - vertex1_coord2_value;
    center_vertex_coord1_value_diff = center_coord1_value - vertex1_coord1_value;
    center_vertex_coord2_value_diff = center_coord2_value - vertex1_coord2_value;
    
    if (vertexes_coord1_value_diff < -180 && is_coord1_unit_degree)
        vertexes_coord1_value_diff += 360;
    if (vertexes_coord1_value_diff > 180 && is_coord1_unit_degree)
        vertexes_coord1_value_diff -= 360;
    if (center_vertex_coord1_value_diff < -180 && is_coord1_unit_degree)
        center_vertex_coord1_value_diff += 360;
    if (center_vertex_coord1_value_diff > 180 && is_coord1_unit_degree)
        center_vertex_coord1_value_diff -= 360;
    if (vertexes_coord2_value_diff < -180 && is_coord2_unit_degree)
        vertexes_coord2_value_diff += 360;
    if (vertexes_coord2_value_diff > 180 && is_coord2_unit_degree)
        vertexes_coord2_value_diff -= 360;
    if (center_vertex_coord2_value_diff < -180 && is_coord2_unit_degree)
        center_vertex_coord2_value_diff += 360;
    if (center_vertex_coord2_value_diff > 180 && is_coord2_unit_degree)
        center_vertex_coord2_value_diff -= 360;

    return vertexes_coord1_value_diff*center_vertex_coord2_value_diff - center_vertex_coord1_value_diff*vertexes_coord2_value_diff;
}


void get_3D_cartesian_coord_of_sphere_coord(double &coord_x, double &coord_y, double &coord_z,
                                            double coord_lon, double coord_lat)
{
    coord_x = cosl(DEGREE_TO_RADIAN(coord_lon))*cosl(DEGREE_TO_RADIAN(coord_lat));
    coord_y = sinl(DEGREE_TO_RADIAN(coord_lon))*cosl(DEGREE_TO_RADIAN(coord_lat));
    coord_z = sinl(DEGREE_TO_RADIAN(coord_lat));
}


bool is_point_in_2D_sphere_coord_cell(double point_sphere_coord_value_lon, 
                                      double point_sphere_coord_value_lat,
                                      double *cell_vertex_sphere_coord_values_lon,
                                      double *cell_vertex_sphere_coord_values_lat,
                                      int num_vertexes)
{
    double point_cartesian_coord_value_x, point_cartesian_coord_value_y, point_cartesian_coord_value_z;
    double cell_vertex_cartesian_coord_values_x[65536], cell_vertex_cartesian_coord_values_y[65536], cell_vertex_cartesian_coord_values_z[65536];
    double current_coord_value_x, current_coord_value_y, current_coord_value_z;
    double next_coord_value_x, next_coord_value_y, next_coord_value_z;
    double current_cross_product, last_cross_product;
    double distance1, distance2, distance3;
    double eps1 = 1.0e-9, eps2 = 1.0e-7;
    int i, next_i;

    
    EXECUTION_REPORT(REPORT_ERROR, -1, num_vertexes <= 65536, "remap software error in is_point_in_2D_sphere_coord_cell: %d\n", num_vertexes);
    
    get_3D_cartesian_coord_of_sphere_coord(point_cartesian_coord_value_x,
                                           point_cartesian_coord_value_y,
                                           point_cartesian_coord_value_z,
                                           point_sphere_coord_value_lon,
                                           point_sphere_coord_value_lat);

    for (i = 0; i < num_vertexes; i ++) {
        if (cell_vertex_sphere_coord_values_lon[i] == NULL_COORD_VALUE)
            continue;
        get_3D_cartesian_coord_of_sphere_coord(cell_vertex_cartesian_coord_values_x[i],
                                               cell_vertex_cartesian_coord_values_y[i],
                                               cell_vertex_cartesian_coord_values_z[i],
                                               cell_vertex_sphere_coord_values_lon[i],
                                               cell_vertex_sphere_coord_values_lat[i]);
    }

    last_cross_product = 0;
    for (i = 0; i < num_vertexes; i ++) {
        if (cell_vertex_sphere_coord_values_lon[i] == NULL_COORD_VALUE)
            continue;
        current_coord_value_x = cell_vertex_cartesian_coord_values_x[i];
        current_coord_value_y = cell_vertex_cartesian_coord_values_y[i];
        current_coord_value_z = cell_vertex_cartesian_coord_values_z[i];
        next_i = (i+1)%num_vertexes;
        while (cell_vertex_sphere_coord_values_lon[next_i] == NULL_COORD_VALUE)
            next_i = (next_i+1)%num_vertexes;
        if (cell_vertex_sphere_coord_values_lon[i] == cell_vertex_sphere_coord_values_lon[next_i] && 
            cell_vertex_sphere_coord_values_lat[i] == cell_vertex_sphere_coord_values_lat[next_i])
            continue;
        next_coord_value_x = cell_vertex_cartesian_coord_values_x[next_i];
        next_coord_value_y = cell_vertex_cartesian_coord_values_y[next_i];
        next_coord_value_z = cell_vertex_cartesian_coord_values_z[next_i];
        current_cross_product = compute_three_3D_points_cross_product(point_cartesian_coord_value_x,
                                                                      point_cartesian_coord_value_y,
                                                                      point_cartesian_coord_value_z,
                                                                      current_coord_value_x,
                                                                      current_coord_value_y,
                                                                      current_coord_value_z,
                                                                      next_coord_value_x,
                                                                      next_coord_value_y,
                                                                      next_coord_value_z);
        if (fabs(current_cross_product) <= eps1) {
            distance1 = calculate_distance_of_two_points_2D(point_sphere_coord_value_lon,
                                                            point_sphere_coord_value_lat,
                                                            cell_vertex_sphere_coord_values_lon[i],
                                                            cell_vertex_sphere_coord_values_lat[i],
                                                            true);
            distance2 = calculate_distance_of_two_points_2D(point_sphere_coord_value_lon,
                                                            point_sphere_coord_value_lat,
                                                            cell_vertex_sphere_coord_values_lon[next_i],
                                                            cell_vertex_sphere_coord_values_lat[next_i],
                                                            true);
            distance3 = calculate_distance_of_two_points_2D(cell_vertex_sphere_coord_values_lon[i],
                                                            cell_vertex_sphere_coord_values_lat[i],
                                                            cell_vertex_sphere_coord_values_lon[next_i],
                                                            cell_vertex_sphere_coord_values_lat[next_i],
                                                            true);
            if (fabs(distance1+distance2-distance3) < eps2)
                return true;
            return false;
        }
        if (last_cross_product != 0 && last_cross_product*current_cross_product < 0)
            return false;
        last_cross_product = current_cross_product;
    }

    return true;
}


bool is_point_in_2D_cartesian_coord_cell(double point_coord1_value, 
                                         double point_coord2_value,
                                         double *cell_vertex_coord1_values,
                                         double *cell_vertex_coord2_values,
                                         int num_vertexes,
                                         bool is_unit_degree_dim1,
                                         bool is_unit_degree_dim2)
{
    int i, j, next_i;
    double current_coord1_value, current_coord2_value, next_coord1_value, next_coord2_value;
    double current_cross_product, last_cross_product;
    double distance1, distance2, distance3;


    last_cross_product = 0;
    for (i = 0; i < num_vertexes; i ++) {
        current_coord1_value = cell_vertex_coord1_values[i];
        current_coord2_value = cell_vertex_coord2_values[i];
        if (current_coord1_value == NULL_COORD_VALUE)
            continue;
        next_i = (i+1)%num_vertexes;
        while (cell_vertex_coord1_values[next_i] == NULL_COORD_VALUE)
            next_i = (next_i+1)%num_vertexes;
        next_coord1_value = cell_vertex_coord1_values[next_i];
        next_coord2_value = cell_vertex_coord2_values[next_i];
        if (current_coord1_value == next_coord1_value && current_coord2_value == next_coord2_value)
            continue;
        current_cross_product = compute_three_2D_points_cross_product(point_coord1_value,
                                                                      point_coord2_value,
                                                                      current_coord1_value,
                                                                      current_coord2_value,
                                                                      next_coord1_value,
                                                                      next_coord2_value,
                                                                      is_unit_degree_dim1,
                                                                      is_unit_degree_dim2);
        if (current_cross_product == 0) {
            distance1 = calculate_distance_of_two_points_2D(point_coord1_value,
                                                            point_coord2_value,
                                                            current_coord1_value,
                                                            current_coord2_value,
                                                            false);
            distance2 = calculate_distance_of_two_points_2D(point_coord1_value,
                                                            point_coord2_value,
                                                            next_coord1_value,
                                                            next_coord2_value,
                                                            false);
            distance3 = calculate_distance_of_two_points_2D(current_coord1_value,
                                                            current_coord2_value,
                                                            next_coord1_value,
                                                            next_coord2_value,
                                                            false);
            if (distance1 <= distance3 && distance2 <= distance3)
                return true;
            return false;
        }
        if (last_cross_product != 0 && last_cross_product*current_cross_product < 0)
            return false;
        last_cross_product = current_cross_product;
    }

    return true;
}


bool is_point_in_2D_cell(double point_coord1_value, 
                        double point_coord2_value,
                        double *cell_vertex_coord1_values,
                        double *cell_vertex_coord2_values,
                        int num_vertexes,
                        bool is_unit_degree_dim1,
                        bool is_unit_degree_dim2,
                        bool is_sphere_grid)
{
    for (int i = 0; i < num_vertexes; i ++) {
        if (cell_vertex_coord1_values[i] == NULL_COORD_VALUE)
            continue;
        if (are_the_same_sphere_points(point_coord1_value, point_coord2_value, cell_vertex_coord1_values[i], cell_vertex_coord2_values[i]))
            return true;
    }

    if (is_sphere_grid)
        return is_point_in_2D_sphere_coord_cell(point_coord1_value,
                                                point_coord2_value, 
                                                cell_vertex_coord1_values, 
                                                cell_vertex_coord2_values, 
                                                num_vertexes);
    else
        return is_point_in_2D_cartesian_coord_cell(point_coord1_value, 
                                                   point_coord2_value, 
                                                   cell_vertex_coord1_values, 
                                                   cell_vertex_coord2_values, 
                                                   num_vertexes, 
                                                   is_unit_degree_dim1, 
                                                   is_unit_degree_dim2);
}


void rotate_sphere_coordinate(double lon_original, double lat_original, double &lon_rotated, double &lat_rotated)
{
    double lon_rotated_radian, lat_rotated_radian, temp1_value, temp2_value;


    if (lon_original == 0.0 && lat_original == 0.0) {
        lat_rotated = -90.0;
        lon_rotated = 0.0;
        return;
    }

    if (lon_original == 180 && lat_original == 0.0) {
        lat_rotated = 90.0;
        lon_rotated = 0.0;
        return;
    }

    temp1_value = cos(DEGREE_TO_RADIAN(lon_original))*cos(DEGREE_TO_RADIAN(lat_original));
    temp2_value = sin(DEGREE_TO_RADIAN(lon_original))*cos(DEGREE_TO_RADIAN(lat_original))/sqrt(1-temp1_value*temp1_value);
    if (temp1_value < -1.0)
        temp1_value = -1.0;
    if (temp1_value > 1.0)
        temp1_value = 1.0;
    if (temp2_value < -1.0)
        temp2_value = -1.0;
    if (temp2_value > 1.0)
        temp2_value = 1.0;    
    lon_rotated_radian = asin(temp2_value);
    lat_rotated_radian = -asin(temp1_value);    
    if (cos(lon_rotated_radian)*cos(lat_rotated_radian)*sin(DEGREE_TO_RADIAN(lat_original)) < 0)
        lon_rotated_radian = PI - lon_rotated_radian;
    if (lon_rotated_radian < 0)
        lon_rotated_radian += 2*PI;
    
    lon_rotated = RADIAN_TO_DEGREE(lon_rotated_radian);
    lat_rotated = RADIAN_TO_DEGREE(lat_rotated_radian);
    
    if (lat_original == 90) {
        lat_rotated = 0;
        lon_rotated = 0;
    }
    else if (lat_original == -90) {
        lat_rotated = 0;
        lon_rotated = 180;
    }

    if (temp1_value == 1.0) {
        lat_rotated = -90;
        lon_rotated = 0;        
    }
    if (temp1_value == -1.0) {
        lat_rotated = 90;
        lon_rotated = 0;        
    }

    lon_rotated = (double) ((float) lon_rotated);
    lat_rotated = (double) ((float) lat_rotated);

    EXECUTION_REPORT(REPORT_ERROR, -1, lon_rotated >= 0 && lon_rotated <= 360, "remap software error1 in rotate_sphere_coordinate\n");
    EXECUTION_REPORT(REPORT_ERROR, -1, lat_rotated >= -90 && lat_rotated <= 90, "remap software error2 in rotate_sphere_coordinate\n");

    if (lon_rotated == 360)
        lon_rotated = 0;
}


