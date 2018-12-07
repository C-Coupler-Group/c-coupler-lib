/***************************************************************
  *  This file was initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#include "cor_global_data.h"
#include "grid_cell_search.h"
#include "remap_utils_nearest_points.h"
#include "quick_sort.h"
#include <math.h>


void seperate_cells_in_children_tiles(int num_cells, H2D_grid_cell_search_cell **cells, double center_lon, double center_lat, 
                                      double dlon, double dlat, int *num_cells_in_children, long *index_buffer)
{
    double new_dlon, new_dlat, min_lon, min_lat, diff_lon, diff_lat;
    int i, child_indx, indx_at_lon, indx_at_lat;
    

    new_dlon = dlon / ((double)TILE_DIVIDE_FACTOR);
    new_dlat = dlat / ((double)TILE_DIVIDE_FACTOR);
    min_lon = center_lon - dlon/2;
    min_lat = center_lat - dlat/2;


    if (min_lon < 0)
        min_lon += 360;
    for (i = 0; i < TILE_DIVIDE_FACTOR*TILE_DIVIDE_FACTOR; i ++) 
        num_cells_in_children[i] = 0;
    for (i = 0; i < num_cells; i ++) {
        diff_lon = cells[i]->get_bounding_circle_center_lon() - min_lon;
        if (diff_lon < 0)
            diff_lon += 360;
        diff_lat = cells[i]->get_bounding_circle_center_lat() - min_lat;
        indx_at_lon = (int)(diff_lon / new_dlon);
        indx_at_lat = (int)(diff_lat / new_dlat);
        if (indx_at_lon == TILE_DIVIDE_FACTOR)
            indx_at_lon = TILE_DIVIDE_FACTOR - 1;
        if (indx_at_lat == TILE_DIVIDE_FACTOR)
            indx_at_lat = TILE_DIVIDE_FACTOR - 1;
        EXECUTION_REPORT(REPORT_ERROR, indx_at_lon < TILE_DIVIDE_FACTOR && indx_at_lat < TILE_DIVIDE_FACTOR, "Software error1 in H2D_grid_cell_search_tile::seperate_cells_in_children_tiles");
        child_indx = indx_at_lat*TILE_DIVIDE_FACTOR+indx_at_lon;
        num_cells_in_children[child_indx] ++;
        index_buffer[i] = child_indx;
    }
}


H2D_grid_cell_cartesian_coord::H2D_grid_cell_cartesian_coord()
{
    center_x = NULL_COORD_VALUE;
    center_y = NULL_COORD_VALUE;
    center_z = NULL_COORD_VALUE;
    vertex_x = NULL;
    vertex_y = NULL;
    vertex_z = NULL;
}


H2D_grid_cell_cartesian_coord::~H2D_grid_cell_cartesian_coord()
{
    if (vertex_x != NULL) {
        delete [] vertex_x;
        delete [] vertex_y;
        delete [] vertex_z;
    }
}


H2D_grid_cell_search_cell::H2D_grid_cell_search_cell(int cell_index, double center_lon, double center_lat, bool mask, 
                                                       int num_vertex, const double *vertex_lons, const double *vertex_lats, int edge_type)
{
    int i, j;


    EXECUTION_REPORT(REPORT_ERROR, edge_type == EDGE_TYPE_LATLON || edge_type == EDGE_TYPE_GREAT_ARC, "Software error1 in H2D_grid_cell_search_cell::H2D_grid_cell_search_cell");
    
    this->cell_index = cell_index;
    this->center_lon = center_lon;
    this->center_lat = center_lat;
    this->mask = mask;
    this->vertex_lons = NULL;
    this->vertex_lats = NULL;
    this->num_vertex = 0;
    this->bounding_circle_radius = 0.0;
    this->bounding_circle_center_lon = center_lon;
    this->bounding_circle_center_lat = center_lat;
    this->edge_type = edge_type;
    this->cartesian_coord = NULL;
    
    if (num_vertex > 0) {
        EXECUTION_REPORT(REPORT_ERROR, vertex_lons != NULL && vertex_lats != NULL, "Software error2 in H2D_grid_cell_search_cell::H2D_grid_cell_search_cell");
        for (i = 0; i < num_vertex; i ++)
            if (vertex_lons[i] == NULL_COORD_VALUE || vertex_lats[i] == NULL_COORD_VALUE) 
                EXECUTION_REPORT(REPORT_ERROR, vertex_lons[i] == NULL_COORD_VALUE && vertex_lats[i] == NULL_COORD_VALUE, "Software error3 in H2D_grid_cell_search_cell::H2D_grid_cell_search_cell");
            else this->num_vertex ++;
        if (this->num_vertex > 0) {
            this->vertex_lons = new double [this->num_vertex];
            this->vertex_lats = new double [this->num_vertex];
            for (i = 0, j = 0; i < num_vertex; i ++) {
                if (vertex_lons[i] == NULL_COORD_VALUE)
                    continue;
                this->vertex_lons[j] = vertex_lons[i];
                this->vertex_lats[j] = vertex_lats[i];
                j ++;
            }
            double min_lon, max_lon;
            min_lon = this->vertex_lons[0];
            max_lon = min_lon;
            for (i = 1; i < this->num_vertex; i ++) {
                if (min_lon > this->vertex_lons[i])
                    min_lon = this->vertex_lons[i];
                if (max_lon < this->vertex_lons[i])
                    max_lon = this->vertex_lons[i];
            }
            double sum_lon = 0.0, sum_lat = 0.0;
            for (i = 0; i < this->num_vertex; i ++) {
                sum_lat += this->vertex_lats[i];
                if (max_lon - min_lon > 180 && min_lon < 180 && this->vertex_lons[i] < 180)
                    sum_lon += this->vertex_lons[i] + 360;
                else sum_lon += this->vertex_lons[i];
            }
            bounding_circle_center_lon = sum_lon / this->num_vertex;
            bounding_circle_center_lat = sum_lat / this->num_vertex;
            if (bounding_circle_center_lon >= 360)
                bounding_circle_center_lon -= 360;
            bounding_circle_radius = 0;
            for (i = 0; i < this->num_vertex; i ++) {
                double dist = calculate_distance_of_two_points_2D(bounding_circle_center_lon, bounding_circle_center_lat, this->vertex_lons[i], this->vertex_lats[i],true);
                if (bounding_circle_radius < dist)
                    bounding_circle_radius = dist;
            }
            bounding_circle_radius *= 1.00001;
        }
    }

    if (edge_type == EDGE_TYPE_GREAT_ARC) {
        cartesian_coord = new H2D_grid_cell_cartesian_coord();
        get_3D_cartesian_coord_of_sphere_coord(cartesian_coord->center_x, cartesian_coord->center_y, cartesian_coord->center_z, center_lon, center_lat);
        if (this->num_vertex > 0) {
            cartesian_coord->vertex_x = new double [this->num_vertex];
            cartesian_coord->vertex_y = new double [this->num_vertex];
            cartesian_coord->vertex_z = new double [this->num_vertex];
            for (i = 0; i < this->num_vertex; i ++)
                get_3D_cartesian_coord_of_sphere_coord(cartesian_coord->vertex_x[i], cartesian_coord->vertex_y[i], cartesian_coord->vertex_z[i], this->vertex_lons[i], this->vertex_lats[i]);
        }
    }
}


H2D_grid_cell_search_cell::~H2D_grid_cell_search_cell()
{
    if (vertex_lons != NULL) {
        delete [] vertex_lons;
        delete [] vertex_lats;
    }

    if (cartesian_coord != NULL)
        delete cartesian_coord;
}


bool H2D_grid_cell_search_cell::is_point_in_latlon_coord_cell(double point_lon,double point_lat) const
{
    int last_pos = -100, i, next_i, pos;
    double det, distance1, distance2, distance3, lon_diff1, lon_diff2;
    bool in_cell = true;
    double e1 = 1.0e-12, e2 = 1.0e-7;
    

    for (i = 0; i < num_vertex; i ++) {
        next_i = (i+1) % num_vertex;
        lon_diff1 = vertex_lons[i]-point_lon;
        lon_diff2 = vertex_lons[next_i]-point_lon;
        if (lon_diff1 < -180)
            lon_diff1 += 360;
        if (lon_diff1 > 180)
            lon_diff1 -= 360;
        if (lon_diff2 < -180)
            lon_diff2 += 360;
        if (lon_diff2 > 180)
            lon_diff2 -= 360;
        det = lon_diff1 * (vertex_lats[next_i]-point_lat) - (vertex_lats[i]-point_lat) * lon_diff2;
        if (fabs(det) <= e1) {
            distance1 = calculate_distance_of_two_points_2D(point_lon, point_lat, vertex_lons[i], vertex_lats[i], true);
            distance2 = calculate_distance_of_two_points_2D(point_lon, point_lat, vertex_lons[next_i], vertex_lats[next_i], true);
            distance3 = calculate_distance_of_two_points_2D(vertex_lons[i], vertex_lats[i], vertex_lons[next_i], vertex_lats[next_i], true);
            in_cell = distance1 <= distance3 && distance2 <= distance3;
            break;
        }
        else if (det > 0)
            pos = 1;
        else pos = -1;
        if (last_pos == -100)
            last_pos = pos;
        else if (last_pos != pos) {
            in_cell = false;
            break;
        }
    }

    return in_cell;
}


bool H2D_grid_cell_search_cell::is_point_in_cartesian_coord_cell(double point_lon, double point_lat) const
{
    double point_cartesian_coord_x, point_cartesian_coord_y, point_cartesian_coord_z;
    double current_coord_value_x, current_coord_value_y, current_coord_value_z;
    double next_coord_value_x, next_coord_value_y, next_coord_value_z;
    double current_cross_product, last_cross_product;
    double distance1, distance2, distance3;
    double eps1 = 1.0e-9, eps2 = 1.0e-7;
    int i, next_i;


    EXECUTION_REPORT(REPORT_ERROR, num_vertex <= 1024, "Software error1 in H2D_grid_cell_search_cell::is_point_in_cartesian_coord_cell");
    
    get_3D_cartesian_coord_of_sphere_coord(point_cartesian_coord_x, point_cartesian_coord_y, point_cartesian_coord_z, point_lon, point_lat);

    last_cross_product = 0;
    for (i = 0; i < num_vertex; i ++) {
        current_coord_value_x = cartesian_coord->vertex_x[i];
        current_coord_value_y = cartesian_coord->vertex_y[i];
        current_coord_value_z = cartesian_coord->vertex_z[i];
        next_i = (i+1)%num_vertex;
        if (vertex_lons[i] == vertex_lons[next_i] && vertex_lats[i] == vertex_lats[next_i])
            continue;
        next_coord_value_x = cartesian_coord->vertex_x[next_i];
        next_coord_value_y = cartesian_coord->vertex_y[next_i];
        next_coord_value_z = cartesian_coord->vertex_z[next_i];
        current_cross_product = compute_three_3D_points_cross_product(point_cartesian_coord_x,
                                                                      point_cartesian_coord_y,
                                                                      point_cartesian_coord_z,
                                                                      current_coord_value_x,
                                                                      current_coord_value_y,
                                                                      current_coord_value_z,
                                                                      next_coord_value_x,
                                                                      next_coord_value_y,
                                                                      next_coord_value_z);
        if (fabs(current_cross_product) <= eps1) {
            distance1 = calculate_distance_of_two_points_2D(point_lon, point_lat, vertex_lons[i], vertex_lats[i], true);
            distance2 = calculate_distance_of_two_points_2D(point_lon, point_lat, vertex_lons[next_i], vertex_lats[next_i], true);
            distance3 = calculate_distance_of_two_points_2D(vertex_lons[i], vertex_lats[i], vertex_lons[next_i], vertex_lats[next_i], true);
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


bool H2D_grid_cell_search_cell::check_possible_overlapping(const H2D_grid_cell_search_cell *src_cell) const
{
    double dist;
    int i;
    bool possible_overlapping = false;


    EXECUTION_REPORT(REPORT_ERROR, src_cell->get_bounding_circle_radius() > 0, "Software error in H2D_grid_cell_search_cell::check_possible_overlapping");
    dist = calculate_distance_of_two_points_2D(bounding_circle_center_lon, bounding_circle_center_lat, src_cell->get_bounding_circle_center_lon(), src_cell->get_bounding_circle_center_lat(), true);
    if (dist > bounding_circle_radius + src_cell->get_bounding_circle_radius())
        return false;
    
    if (num_vertex == 0)
        return true;

    for (i = 0; i < num_vertex; i ++) {
        dist = calculate_distance_of_two_points_2D(vertex_lons[i], vertex_lats[i], src_cell->get_bounding_circle_center_lon(), src_cell->get_bounding_circle_center_lat(), true);
        if (dist <= src_cell->get_bounding_circle_radius()) {
            possible_overlapping = true;
            break;
        }
    }

    for (i = 0; i < src_cell->get_num_vertex(); i ++) {
        dist = calculate_distance_of_two_points_2D(src_cell->get_vertex_lon(i), src_cell->get_vertex_lat(i), bounding_circle_center_lon, bounding_circle_center_lat, true);
        if (dist <= bounding_circle_radius) {
            possible_overlapping = true;
            break;
        }
    }

    return possible_overlapping;
}


bool H2D_grid_cell_search_cell::check_overlapping(const H2D_grid_cell_search_cell *src_cell, bool accurately_match) const
{
    bool possible_overlapping, is_overlapping;
    int i;


    if (!src_cell->get_mask())
        return false;

    possible_overlapping = check_possible_overlapping(src_cell);

    if (!accurately_match || !possible_overlapping)
        return possible_overlapping;

    if (num_vertex > 0) {
        for (i = 0; i < num_vertex; i ++) {
            if (src_cell->get_edge_type() == EDGE_TYPE_LATLON)
                is_overlapping = src_cell->is_point_in_latlon_coord_cell(vertex_lons[i], vertex_lats[i]);
            else is_overlapping = src_cell->is_point_in_cartesian_coord_cell(vertex_lons[i], vertex_lats[i]); 
            if (is_overlapping)
                return true;
        }
        for (i = 0; i < src_cell->get_num_vertex(); i ++) {
            if (edge_type == EDGE_TYPE_LATLON)
                is_overlapping = is_point_in_latlon_coord_cell(src_cell->get_vertex_lon(i), src_cell->get_vertex_lat(i));
            else is_overlapping = is_point_in_cartesian_coord_cell(src_cell->get_vertex_lon(i), src_cell->get_vertex_lat(i));
            if (is_overlapping)
                return true;
        }
    }
    else {
        if (edge_type == EDGE_TYPE_LATLON)
            is_overlapping = src_cell->is_point_in_latlon_coord_cell(center_lon, center_lat);
        else is_overlapping = src_cell->is_point_in_cartesian_coord_cell(center_lon, center_lat);
        if (is_overlapping)
            return true;
    }

    return false;
}


double H2D_grid_cell_search_cell::get_vertex_lon(int indx) const 
{
    EXECUTION_REPORT(REPORT_ERROR, indx < num_vertex && indx >= 0, "Software error in H2D_grid_cell_search_cell::get_vertex_lon");
    return vertex_lons[indx];
}


double H2D_grid_cell_search_cell::get_vertex_lat(int indx) const 
{
    EXECUTION_REPORT(REPORT_ERROR, indx < num_vertex && indx >= 0, "Software error in H2D_grid_cell_search_cell::get_vertex_lat");
    return vertex_lats[indx];
}


H2D_grid_cell_search_tile::H2D_grid_cell_search_tile(int num_cells, H2D_grid_cell_search_cell **cells, H2D_grid_cell_search_cell **cells_buffer, 
                                                     long *index_buffer, H2D_grid_cell_search_tile *parent, double center_lon, double center_lat, double dlon, double dlat)
{
    this->num_cells = num_cells;
    this->cells = cells;
    this->cells_buffer = cells_buffer;
    this->index_buffer = index_buffer;
    this->parent = parent;
    this->center_lon = center_lon;
    this->center_lat = center_lat;
    this->dlon = dlon;
    this->dlat = dlat;

    compute_bounding_circle();
    divide_tile();
}


H2D_grid_cell_search_tile::~H2D_grid_cell_search_tile()
{
    for (int i = 0; i < TILE_DIVIDE_FACTOR*TILE_DIVIDE_FACTOR; i ++)
        if (children[i] != NULL)
            delete children[i];
}



void H2D_grid_cell_search_tile::compute_bounding_circle()
{
    double point_lon, point_lat, dist, max_dist;
    int i, j;

    point_lon = center_lon - dlon/2;
    point_lat = center_lat - dlat/2;
    if (point_lon < 0)
        point_lon += 360;
    circle_radius = calculate_distance_of_two_points_2D(center_lon, center_lat, point_lon, point_lat, true);

    for (i = 0; i < num_cells; i ++) {
        dist = calculate_distance_of_two_points_2D(center_lon, center_lat, cells[i]->get_bounding_circle_center_lon(), cells[i]->get_bounding_circle_center_lat(), true);
        if (circle_radius < dist + cells[i]->get_bounding_circle_radius()) {
            circle_radius = dist + cells[i]->get_bounding_circle_radius();
        }
    }
    circle_radius *= 1.00001;
}


void H2D_grid_cell_search_tile::divide_tile()
{
    int i, j, indx_at_lon, indx_at_lat, child_indx, displ;
    int num_cells_in_children[TILE_DIVIDE_FACTOR*TILE_DIVIDE_FACTOR], displ_of_children_in_cells[TILE_DIVIDE_FACTOR*TILE_DIVIDE_FACTOR];
    double new_dlon, new_dlat, min_lon, min_lat, diff_lon, diff_lat;
    double children_center_lon[TILE_DIVIDE_FACTOR*TILE_DIVIDE_FACTOR], children_center_lat[TILE_DIVIDE_FACTOR*TILE_DIVIDE_FACTOR];


    for (i = 0; i < TILE_DIVIDE_FACTOR*TILE_DIVIDE_FACTOR; i ++)
        children[i] = NULL;
    
    if (num_cells <= MAX_NUM_CELLS_IN_TILE)
        return;

    seperate_cells_in_children_tiles(num_cells, cells, center_lon, center_lat, dlon, dlat, num_cells_in_children, index_buffer);

    new_dlon = dlon / ((double)TILE_DIVIDE_FACTOR);
    new_dlat = dlat / ((double)TILE_DIVIDE_FACTOR);
    min_lon = center_lon - dlon/2;
    min_lat = center_lat - dlat/2;
    if (min_lon < 0)
        min_lon += 360;

    for (i = 0, displ = 0; i < TILE_DIVIDE_FACTOR*TILE_DIVIDE_FACTOR; i ++) {
        displ_of_children_in_cells[i] = displ;
        displ += num_cells_in_children[i];
    }
    for (i = 0; i < num_cells; i ++)
        cells_buffer[displ_of_children_in_cells[index_buffer[i]]++] = cells[i];
    for (i = 0; i < num_cells; i ++)
        cells[i] = cells_buffer[i];
    for (i = 0, displ = 0; i < TILE_DIVIDE_FACTOR*TILE_DIVIDE_FACTOR; i ++) {
        EXECUTION_REPORT(REPORT_ERROR, displ_of_children_in_cells[i] == displ+num_cells_in_children[i], "Software error2 in H2D_grid_cell_search_tile::divide_tile");
        displ_of_children_in_cells[i] = displ;
        displ += num_cells_in_children[i];
    }    

    for (indx_at_lat = 0; indx_at_lat < TILE_DIVIDE_FACTOR; indx_at_lat ++)
        for (indx_at_lon = 0; indx_at_lon <  TILE_DIVIDE_FACTOR; indx_at_lon ++) {
            child_indx = indx_at_lat*TILE_DIVIDE_FACTOR+indx_at_lon;
            children_center_lon[child_indx] = min_lon + new_dlon*indx_at_lon + new_dlon/2;
            children_center_lat[child_indx] = min_lat + new_dlat*indx_at_lat + new_dlat/2;
            if (children_center_lon[child_indx] < 0)
                children_center_lon[child_indx] += 360;
            if (children_center_lon[child_indx] >= 360)
                children_center_lon[child_indx] -= 360;
        }
    for (i = 0, j = 0; i < TILE_DIVIDE_FACTOR*TILE_DIVIDE_FACTOR; i ++)
        if (num_cells_in_children[i] != 0) {
            children[j++] = new H2D_grid_cell_search_tile(num_cells_in_children[i], cells+displ_of_children_in_cells[i], cells_buffer, index_buffer,
                                                          this, children_center_lon[i], children_center_lat[i], new_dlon, new_dlat);
        }

    if (j > 0)
        cells = NULL;
}


bool H2D_grid_cell_search_tile::search_points_within_distance(double dist_threshold, double dst_point_lon, double dst_point_lat, int &num_found_points, long *found_points_indx, double *found_points_dist, bool early_quit)
{
    double distance;
    bool have_the_same_point = false;


    distance = calculate_distance_of_two_points_2D(center_lon, center_lat, dst_point_lon, dst_point_lat, true);
    if (distance - circle_radius > dist_threshold)
        return have_the_same_point;

    if (cells == NULL) {
        for (int i = 0; i < TILE_DIVIDE_FACTOR*TILE_DIVIDE_FACTOR; i ++) {
            if (children[i] == NULL) 
                break;
            have_the_same_point |= children[i]->search_points_within_distance(dist_threshold, dst_point_lon, dst_point_lat, num_found_points, found_points_indx, found_points_dist, early_quit);
            if (early_quit && have_the_same_point)
                break;
        }
    }
    else {
        for (int i = 0; i < num_cells; i ++) {
            if (!cells[i]->get_mask())
                continue;
            if (cells[i]->get_center_lon() == dst_point_lon && cells[i]->get_center_lat() == dst_point_lat) {
                have_the_same_point = true;
                distance = 0;
            }
            else distance = calculate_distance_of_two_points_2D(cells[i]->get_center_lon(), cells[i]->get_center_lat(), dst_point_lon, dst_point_lat, true);
            if (distance <= dist_threshold) {
                found_points_indx[num_found_points] = cells[i]->get_cell_index();
                found_points_dist[num_found_points] = distance;
                num_found_points ++;
            }
            if (early_quit && have_the_same_point) {
                found_points_indx[0] = cells[i]->get_cell_index();
                found_points_dist[0] = 0;
                num_found_points = 1;
                break;
            }
        }
    }

    return have_the_same_point;
}


bool H2D_grid_cell_search_tile::has_cell_index(int cell_index)
{
    if (cells == NULL) {
        for (int i = 0; i < TILE_DIVIDE_FACTOR*TILE_DIVIDE_FACTOR; i ++) {
            if (children[i] == NULL) 
                break;
            if (children[i]->has_cell_index(cell_index))
                return true;
        }
    }
    else {
        for (int i = 0; i < num_cells; i ++)
            if (cells[i]->get_cell_index() == cell_index)
                return true;
    }

    return false;
}


void H2D_grid_cell_search_tile::search_overlapping_cells(int &num_overlapping_cells, long *overlapping_cells_index, const H2D_grid_cell_search_cell *dst_cell, bool accurately_match, bool early_quit)
{
    double dist = calculate_distance_of_two_points_2D(dst_cell->get_bounding_circle_center_lon(), dst_cell->get_bounding_circle_center_lat(), center_lon, center_lat, true);

    if (dist > dst_cell->get_bounding_circle_radius() + circle_radius)
        return;

    if (cells == NULL) {
        for (int i = 0; i < TILE_DIVIDE_FACTOR*TILE_DIVIDE_FACTOR; i ++) {
            if (children[i] == NULL) 
                break;
            children[i]->search_overlapping_cells(num_overlapping_cells, overlapping_cells_index, dst_cell, accurately_match, early_quit);
            if (early_quit && num_overlapping_cells > 0)
                return;
        }

    }
    else {
        for (int i = 0; i < num_cells; i ++) {
            if (dst_cell->check_overlapping(cells[i], accurately_match)) {
                overlapping_cells_index[num_overlapping_cells++] = cells[i]->get_cell_index();
                if (early_quit)
                    return;
            }
        }
    }
}


H2D_grid_cell_search_engine::H2D_grid_cell_search_engine(const Remap_grid_class *remap_grid, const double *center_lons, const double *center_lats, const bool *masks, 
                                                         const bool *redundant_mask, int num_vertex, const double *vertex_lons, const double *vertex_lats, int edge_type, bool build_search_structure)
{
    double center_lon, center_lat, dlon, dlat;
    bool mask = true;
    int i;

    
    EXECUTION_REPORT(REPORT_ERROR, remap_grid->get_is_sphere_grid(), "Software error1 in H2D_grid_cell_search_engine::H2D_grid_cell_search_engine");
    
    this->remap_grid = remap_grid;
    if (redundant_mask == NULL)
        num_cells = remap_grid->get_grid_size();
    else {
        for (i = 0, num_cells = 0; i < remap_grid->get_grid_size(); i ++)
            if (!redundant_mask[i])
                num_cells ++;
    }
        
    cells = new H2D_grid_cell_search_cell* [remap_grid->get_grid_size()];
    cells_ptr = new H2D_grid_cell_search_cell* [num_cells];
    cells_buffer = new H2D_grid_cell_search_cell* [num_cells];
    index_buffer = new long [num_cells];
    dist_buffer = new double [num_cells];
    dist_threshold = 1 / 6000.0;

    for (int i = 0, num_cells = 0; i < remap_grid->get_grid_size(); i ++) {
        if (masks != NULL)
            mask = masks[i];
        cells[i] = new H2D_grid_cell_search_cell(i, center_lons[i], center_lats[i], mask, num_vertex, vertex_lons+num_vertex*i, vertex_lats+num_vertex*i, edge_type);
        if (redundant_mask != NULL && redundant_mask[i])
            continue;
        cells_ptr[num_cells] = cells[i];
        num_cells++;
    }

    center_lon = NULL_COORD_VALUE;
    center_lat = NULL_COORD_VALUE;
    dlon = NULL_COORD_VALUE;
    dlat = NULL_COORD_VALUE;
    recursively_search_initial_boundary(180, 0, 360, 180, center_lon, center_lat, dlon, dlat);
    EXECUTION_REPORT(REPORT_ERROR, center_lon != NULL_COORD_VALUE && center_lat != NULL_COORD_VALUE && dlon != NULL_COORD_VALUE && dlat != NULL_COORD_VALUE, 
                     "Software error2 in in H2D_grid_cell_search_engine::H2D_grid_cell_search_engine");
    
    if (build_search_structure)
        root_tile = new H2D_grid_cell_search_tile(num_cells, cells_ptr, cells_buffer, index_buffer, NULL, center_lon, center_lat, dlon, dlat);
    else root_tile = NULL;
}


H2D_grid_cell_search_engine::~H2D_grid_cell_search_engine()
{
    for (int i = 0; i < remap_grid->get_grid_size(); i ++)
        delete cells[i];
    delete [] cells;
    delete [] cells_ptr;
    delete [] cells_buffer;
    delete [] index_buffer;
    delete [] dist_buffer;
    delete root_tile;
}


void H2D_grid_cell_search_engine::recursively_search_initial_boundary(double center_lon, double center_lat, double dlon, double dlat, double &result_center_lon, 
                                                                      double &result_center_lat, double &result_dlon, double &result_dlat)
{
    int num_cells_in_children[TILE_DIVIDE_FACTOR*TILE_DIVIDE_FACTOR];
    int i, num_children, unique_child_indx, indx_at_lon, indx_at_lat;
    double new_dlon, new_dlat, min_lon, min_lat, new_center_lon, new_center_lat;

    
    seperate_cells_in_children_tiles(num_cells, cells_ptr, center_lon, center_lat, dlon, dlat, num_cells_in_children, index_buffer);
    for (i = 0, num_children = 0; i < TILE_DIVIDE_FACTOR*TILE_DIVIDE_FACTOR; i ++)
        if (num_cells_in_children[i] > 0) {
            unique_child_indx = i;
            num_children ++;
        }

    if (num_children > 1) {
        result_center_lon = center_lon;
        result_center_lat = center_lat;
        result_dlon = dlon;
        result_dlat = dlat;
        return;
    }

    new_dlon = dlon / ((double)TILE_DIVIDE_FACTOR);
    new_dlat = dlat / ((double)TILE_DIVIDE_FACTOR);
    min_lon = center_lon - dlon/2;
    min_lat = center_lat - dlat/2;
    if (min_lon < 0)
        min_lon += 360;
    indx_at_lon = unique_child_indx % TILE_DIVIDE_FACTOR;
    indx_at_lat = unique_child_indx / TILE_DIVIDE_FACTOR;
    new_center_lon = min_lon + new_dlon*indx_at_lon + new_dlon/2;
    new_center_lat = min_lat + new_dlat*indx_at_lat + new_dlat/2;
    if (new_center_lon >= 360)
        new_center_lon -= 360;
    
    recursively_search_initial_boundary(new_center_lon, new_center_lat, new_dlon, new_dlat, result_center_lon, result_center_lat, result_dlon, result_dlat);
}


void H2D_grid_cell_search_engine::search_nearest_points_var_number(int num_required_points, double dst_point_lon, double dst_point_lat, int &num_found_points, long *found_points_indx, double *found_points_dist, bool early_quit)
{
    bool have_the_same_point = false;


    if (num_required_points > remap_grid->get_grid_size())
        num_required_points = remap_grid->get_grid_size();
    
    num_found_points = 0;
    
    while (num_found_points < num_required_points) {
        num_found_points = 0;
        have_the_same_point = root_tile->search_points_within_distance(dist_threshold, dst_point_lon, dst_point_lat, num_found_points, index_buffer, dist_buffer, early_quit);
        if (num_found_points == 0) {
            dist_threshold *= 2;
			continue;
        }
        if (have_the_same_point && early_quit)
            break;
        dist_threshold *= sqrt(((double)num_required_points)/((double)num_found_points))*1.1;
    }

    if (have_the_same_point && early_quit) {
        EXECUTION_REPORT(REPORT_ERROR, num_found_points == 1, "Software error in H2D_grid_cell_search_engine::search_nearest_points_var_number");
        found_points_indx[0] = index_buffer[0];
        found_points_dist[0] = dist_buffer[0];
        return;
    }
    
    do_quick_sort(dist_buffer, index_buffer, 0, num_found_points-1);

    for (int i = 0; i < num_required_points; i ++) {
        found_points_indx[i] = index_buffer[i];
        found_points_dist[i] = dist_buffer[i];
    }
    num_found_points = num_required_points;
}


void H2D_grid_cell_search_engine::search_nearest_points_var_distance(double dist_threshold, double dst_point_lon, double dst_point_lat, int &num_found_points, long *found_points_indx, double *found_points_dist, bool early_quit)
{
    bool have_the_same_point;


    EXECUTION_REPORT(REPORT_ERROR, root_tile != NULL, "Software error1 in H2D_grid_cell_search_engine::search_nearest_points_var_distance");
    
    this->dist_threshold = dist_threshold;
    num_found_points = 0;
    have_the_same_point = root_tile->search_points_within_distance(dist_threshold, dst_point_lon, dst_point_lat, num_found_points, index_buffer, dist_buffer, early_quit);

    do_quick_sort(dist_buffer, index_buffer, 0, num_found_points-1);
    
    if (have_the_same_point && early_quit) {
        EXECUTION_REPORT(REPORT_ERROR, num_found_points == 1, "Software error2 in H2D_grid_cell_search_engine::search_nearest_points_var_distance");
        found_points_indx[0] = index_buffer[0];
        found_points_dist[0] = dist_buffer[0];
        return;
    }

    for (int i = 0; i < num_found_points; i ++) {
        found_points_indx[i] = index_buffer[i];
        found_points_dist[i] = dist_buffer[i];
    }
}




void H2D_grid_cell_search_engine::search_overlapping_cells(int &num_overlapping_cells, long *overlapping_cells_index, const H2D_grid_cell_search_cell *dst_cell, bool accurately_match, bool early_quit) const
{
    EXECUTION_REPORT(REPORT_ERROR, -1, root_tile != NULL, "Software error1 in H2D_grid_cell_search_engine::search_overlapping_cells");

    num_overlapping_cells = 0;
    root_tile->search_overlapping_cells(num_overlapping_cells, index_buffer, dst_cell, accurately_match, early_quit);
    
    if (early_quit)
        EXECUTION_REPORT(REPORT_ERROR, -1, num_overlapping_cells <= 1, "Software error2 in H2D_grid_cell_search_engine::search_overlapping_cells %d", num_overlapping_cells);
        
    for (int i = 0; i < num_overlapping_cells; i ++)
        overlapping_cells_index[i] = index_buffer[i];

    if (early_quit && num_overlapping_cells > 0)
        return;
    
    do_quick_sort(overlapping_cells_index, (long*) NULL, 0, num_overlapping_cells-1);
/*
    int temp_num_overlapping_cells = 0;
    for (int i = 0; i < num_cells; i ++) {
        if (!cells[i]->get_mask())
            continue;
        bool result1 = cells[i]->check_overlapping(dst_cell, accurately_match);
        bool result2 = dst_cell->check_overlapping(cells[i], accurately_match);
        if (result1 != result2) {
            printf("error1 in H2D_grid_cell_search_engine::search_overlapping_cells for %d %d\n", dst_cell->get_cell_index(), cells[i]->get_cell_index());
            
        }
        if (cells[i]->check_overlapping(dst_cell, accurately_match)) {
            index_buffer[temp_num_overlapping_cells++] = cells[i]->get_cell_index();
        }
    }
    do_quick_sort(index_buffer, (long*) NULL, 0, temp_num_overlapping_cells-1);
    
    bool same_result = temp_num_overlapping_cells == num_overlapping_cells;
    if (same_result)
        for (int i = 0; i < num_overlapping_cells; i ++)
            if (index_buffer[i] != overlapping_cells_index[i]) {
                same_result = false;
                break;
            }
    if (!same_result) {
        printf("error in H2D_grid_cell_search_engine::search_overlapping_cells\n");
        printf("result of all scan %d: ", temp_num_overlapping_cells);
        for (int i = 0; i < temp_num_overlapping_cells; i ++)
            printf("%d ", index_buffer[i]);
        printf("\n");
        printf("result of fast %d: ", num_overlapping_cells);
        for (int i = 0; i < num_overlapping_cells; i ++)
            printf("%d ", overlapping_cells_index[i]);
        printf("\n");
        EXECUTION_REPORT(REPORT_ERROR, -1, false, "error");
    }
    */
}


int H2D_grid_cell_search_engine::search_cell_of_locating_point(double point_lon, double point_lat, bool accurately_match) const
{
    int num_overlapping_cells;
    H2D_grid_cell_search_cell *temp_cell;


    temp_cell = new H2D_grid_cell_search_cell(0, point_lon, point_lat, true, 0, NULL, NULL, EDGE_TYPE_LATLON);

    search_overlapping_cells(num_overlapping_cells, index_buffer, temp_cell, accurately_match, true);

    if (num_overlapping_cells == 0)
        return -1;
    else return index_buffer[0];
    
    delete temp_cell;
}


void H2D_grid_cell_search_engine::update(const bool *new_masks)
{
    if (new_masks == NULL)
        return;

    for (int i = 0; i < remap_grid->get_grid_size(); i ++)
        cells[i]->set_mask(new_masks[i]);
}


const H2D_grid_cell_search_cell* H2D_grid_cell_search_engine::get_cell(int cell_index) const
{
        EXECUTION_REPORT(REPORT_ERROR, -1, root_tile == NULL, "Software error1 in H2D_grid_cell_search_engine::get_cell");
        EXECUTION_REPORT(REPORT_ERROR, -1, cell_index >= 0 && cell_index < remap_grid->get_grid_size(), "Software error2 in H2D_grid_cell_search_engine::get_cell");
        EXECUTION_REPORT(REPORT_ERROR, -1, cells[cell_index]->get_mask(), "Software error3 in H2D_grid_cell_search_engine::get_cell");

        return cells[cell_index];
}

