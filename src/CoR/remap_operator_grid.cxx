/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#include "cor_global_data.h"
#include "remap_operator_grid.h"
#include "radix_sort.h"
#include "remap_operator_c_interface.h"


Vertex_values_group::Vertex_values_group()
{
    num_vertex = 0;
    num_coords = 0;
    coord_value_grid = NULL;
    partial_nerghbors_indexes = NULL;
    num_neighbors = 0;
}


Vertex_values_group::~Vertex_values_group()
{
    if (partial_nerghbors_indexes != NULL)
        delete [] partial_nerghbors_indexes;
}


void Vertex_values_group::generate_vertex_values_of_one_cell_according_to_vertex_values_groups(long whole_grid_cell_index,
                                                                                              long *group_grid_cell_index,
                                                                                              int num_vertex_values_groups,
                                                                                              int num_grid_dimensions,
                                                                                              int num_vertexes,
                                                                                              Vertex_values_group **vertex_values_groups,
                                                                                              double **output_vertex_values)
{
    int group_index, i, j, k;
    int vertex_id;
    double stack_values[256][256];
    int stack_top, num_processed_coords, current_order;


    for (k = 0; k < num_grid_dimensions; k ++)
        for (j = 0; j < num_vertexes; j ++)
            output_vertex_values[k][whole_grid_cell_index*num_vertexes+j] = NULL_COORD_VALUE;
    
    num_processed_coords = 0;
    for (group_index = 0, j = 0; group_index < num_vertex_values_groups; group_index ++) 
        for (i = 0; i < vertex_values_groups[group_index]->num_coords; i ++)
            stack_values[num_processed_coords++][0] = vertex_values_groups[group_index]->vertex_coord_values[i][group_grid_cell_index[group_index]*vertex_values_groups[group_index]->num_vertex];
    stack_top = 1;
    num_processed_coords = 0;
    for (group_index = 0; group_index < num_vertex_values_groups; group_index ++) {
        current_order = 1;
        for (vertex_id = 0; vertex_id < vertex_values_groups[group_index]->num_vertex; vertex_id ++) {
            if (vertex_values_groups[group_index]->vertex_coord_values[0][group_grid_cell_index[group_index]*vertex_values_groups[group_index]->num_vertex+vertex_id] == NULL_COORD_VALUE)
                continue;
            if (current_order == 1) 
                for (j = 0, i = 0; j < stack_top; j ++) {
                    for (k = 0; k < num_grid_dimensions; k ++) 
                        output_vertex_values[k][whole_grid_cell_index*num_vertexes+i] = stack_values[k][j];
                    i ++;
                }
            else for (j = stack_top-1, i = 0; j >= 0; j --) {
                    for (k = 0; k < num_grid_dimensions; k ++)
                        output_vertex_values[k][whole_grid_cell_index*num_vertexes+i] = stack_values[k][j];
                    i ++;
                }
            for (j = 0; j < stack_top; j ++) 
                for (k = 0; k < vertex_values_groups[group_index]->num_coords; k ++)
                    output_vertex_values[num_processed_coords+k][whole_grid_cell_index*num_vertexes+j] = vertex_values_groups[group_index]->vertex_coord_values[k][group_grid_cell_index[group_index]*vertex_values_groups[group_index]->num_vertex+vertex_id];
            for (j = 0; j < stack_top; j ++)
                for (k = 0; k < num_grid_dimensions; k ++)
                    stack_values[k][stack_top*vertex_id+j] = output_vertex_values[k][whole_grid_cell_index*num_vertexes+j];
            current_order = -current_order;
        }
        stack_top *= vertex_values_groups[group_index]->num_vertex;
        num_processed_coords += vertex_values_groups[group_index]->num_coords;
    }

    for (j = 0; j < stack_top; j ++) 
        for (k = 0; k < num_grid_dimensions; k ++) 
            output_vertex_values[k][whole_grid_cell_index*num_vertexes+j] = stack_values[k][j];

    EXECUTION_REPORT(REPORT_ERROR, -1, stack_top == num_vertexes, "remap software error1 in generate_overall_vertex_coord_values\n");
    EXECUTION_REPORT(REPORT_ERROR, -1, num_processed_coords == num_grid_dimensions, "remap software error2 in generate_overall_vertex_coord_values\n");
}




Remap_operator_grid::Remap_operator_grid(Remap_grid_class *remap_grid, Remap_operator_basis *remap_operator, bool is_src_grid, bool is_rotated_grid)
{
    int num_leaf_grids;
    long i, j;
    Remap_grid_class *leaf_grids[256];
    Remap_grid_data_class *current_grid_center_field;


    this->remap_operator = remap_operator;
    this->remap_grid = remap_grid;
    this->num_grid_dimensions = remap_grid->num_dimensions;
    this->grid_size = remap_grid->grid_size;
    this->require_vertex_fields = false;
    this->num_vertexes = 0;
    this->num_neighbors = 0;
    this->num_vertex_values_groups = 0;
    this->grid2D_search_engine = NULL;
    this->rotated_remap_operator_grid = NULL;
    this->is_src_grid = is_src_grid;
    this->is_rotated_grid = is_rotated_grid;
    this->redundant_cell_mark = remap_grid->redundant_cell_mark;
    this->cell_visiting_mark = new bool [grid_size];
    this->is_grid_sphere = remap_grid->get_is_sphere_grid();
    for (i = 0; i < grid_size; i ++)
        cell_visiting_mark[i] = false;
    visited_cells_indexes = new long [grid_size];
    num_visited_cells = 0;
    for (i = 0; i < 256; i ++) {
        center_coord_values[i] = NULL;
        vertex_coord_values[i] = NULL;
    }
    if (remap_grid->grid_mask_field != NULL)
        mask_values = (bool*) remap_grid->grid_mask_field->grid_data_field->data_buf;
    else mask_values = NULL;
    
    remap_grid->get_leaf_grids(&num_leaf_grids, leaf_grids, remap_grid);
    for (i = 0; i < num_leaf_grids; i ++)
        if (words_are_the_same(leaf_grids[i]->get_coord_unit(), COORD_UNIT_DEGREES))
            is_coord_unit_degree[i] = true;
        else is_coord_unit_degree[i] = false;

    if (remap_operator->get_is_operator_regridding()) {
        for (i = 0; i < num_leaf_grids; i ++) {
            if (leaf_grids[i]->grid_center_fields.size() == 0)
                current_grid_center_field = leaf_grids[i]->get_grid_center_field();
            else current_grid_center_field = leaf_grids[i]->grid_center_fields[0];
            EXECUTION_REPORT(REPORT_ERROR, -1, current_grid_center_field != NULL, "remap software error1 in Remap_operator_grid\n");
            if (num_leaf_grids == 1)
                grid_center_fields.push_back(current_grid_center_field);
            else grid_center_fields.push_back(remap_grid->expand_to_generate_full_coord_value(current_grid_center_field));
            center_coord_values[i] = (double *) grid_center_fields[i]->grid_data_field->data_buf;
        }
    }

    if (remap_operator->get_num_dimensions() == 1)
        require_vertex_fields = remap_operator->does_require_grid_vertex_values();
    else require_vertex_fields = remap_operator->does_require_grid_vertex_values() || (remap_operator->get_is_operator_regridding());

    if (is_grid_sphere && !is_rotated_grid && remap_operator->get_is_operator_regridding()) {
        rotated_remap_operator_grid = new Remap_operator_grid(remap_grid, remap_operator, is_src_grid, true);
    }
}


Remap_operator_grid::~Remap_operator_grid()
{
    delete [] cell_visiting_mark;
    delete [] visited_cells_indexes;

    if (grid_center_fields.size() > 1)
        for (int i = 0; i < grid_center_fields.size(); i ++) 
            delete grid_center_fields[i];

    if (num_vertexes > 0)
        for (int i = 0; i < num_grid_dimensions; i ++)
            delete [] vertex_coord_values[i];

    for (int i = 0; i < num_vertex_values_groups; i ++)
        delete vertex_values_groups[i];

    if (grid2D_search_engine != NULL)
        delete grid2D_search_engine;

    if (rotated_remap_operator_grid != NULL)
        delete rotated_remap_operator_grid;
}


void Remap_operator_grid::update_operator_grid_data()
{
    if (require_vertex_fields) {
        initialize_for_vertex_coord_values_generation();
        generate_overall_vertex_coord_values();
    }

    if (is_src_grid)
        remove_redundant_cells();

    if (is_rotated_grid)
        rotate_sphere_grid();
    
    if (remap_operator->get_num_dimensions() > 1 && remap_operator->get_is_operator_regridding()) {
        if (grid2D_search_engine == NULL)
            grid2D_search_engine = new H2D_grid_cell_search_engine(remap_grid, center_coord_values[0], center_coord_values[1], mask_values, redundant_cell_mark,
                                                                   num_vertexes, vertex_coord_values[0], vertex_coord_values[1], EDGE_TYPE_LATLON, is_src_grid);
        else grid2D_search_engine->update(mask_values);
    }

    if (rotated_remap_operator_grid != NULL)
        rotated_remap_operator_grid->update_operator_grid_data();
}


void Remap_operator_grid::remove_redundant_cells()
{
    long i, j, k;


    if (center_coord_values[0] == NULL || redundant_cell_mark == NULL)
        return;

    for (i = 0; i < grid_size; i ++)
        if (redundant_cell_mark[i]) {
            for (j = 0; j < num_grid_dimensions; j ++) {
                center_coord_values[j][i] = NULL_COORD_VALUE;
                if (vertex_coord_values[j] != NULL)
                    for (k = 0; k < num_vertexes; k ++)
                        vertex_coord_values[j][i*num_vertexes+k] = NULL_COORD_VALUE;
            }
        }
}


void Remap_operator_grid::rotate_sphere_grid()
{
    int lon_index, lat_index, num_leaf_grids;
    long i, j, next_j, array_size;
    Remap_grid_class *leaf_grids[256];
    double *lon_data_array, *lat_data_array, lon_diff_value;
    double eps = 1.0e-5;
    

    remap_grid->get_leaf_grids(&num_leaf_grids, leaf_grids, remap_grid);
    for (i = 0; i < num_leaf_grids; i ++)
        if (words_are_the_same(leaf_grids[i]->get_coord_label(), COORD_LABEL_LON))
            lon_index = i;
        else if (words_are_the_same(leaf_grids[i]->get_coord_label(), COORD_LABEL_LAT))
            lat_index = i;

    array_size = grid_size*num_vertexes;
    lon_data_array = vertex_coord_values[lon_index];
    lat_data_array = vertex_coord_values[lat_index];
    for (i = 0; i < array_size; i ++)
        if (lon_data_array[i] != NULL_COORD_VALUE)
            rotate_sphere_coordinate(lon_data_array[i], lat_data_array[i], lon_data_array[i], lat_data_array[i]);
    for (i = 0; i < grid_size; i ++) {
        for (j = 0; j < num_vertexes; j ++) {
            if (lon_data_array[i*num_vertexes+j] == NULL_COORD_VALUE)
                continue;
            next_j = (j+1)%num_vertexes;
            while(lon_data_array[i*num_vertexes+next_j] == NULL_COORD_VALUE)
                next_j = (next_j+1)%num_vertexes;
            lon_diff_value = compute_difference_of_two_coord_values(lon_data_array[i*num_vertexes+j], lon_data_array[i*num_vertexes+next_j], lon_index);
            if (fabs(lon_diff_value) < eps && fabs(lat_data_array[i*num_vertexes+j] - lat_data_array[i*num_vertexes+next_j]) < eps) {
                lon_data_array[i*num_vertexes+next_j] = NULL_COORD_VALUE;
                lat_data_array[i*num_vertexes+next_j] = NULL_COORD_VALUE;
            }
        }
    }

    array_size = grid_center_fields[lon_index]->get_grid_data_field()->required_data_size;
    lon_data_array = (double*) grid_center_fields[lon_index]->get_grid_data_field()->data_buf;
    lat_data_array = (double*) grid_center_fields[lat_index]->get_grid_data_field()->data_buf;
    for (i = 0; i < array_size; i ++) 
        if (lon_data_array[i] != NULL_COORD_VALUE)    
            rotate_sphere_coordinate(lon_data_array[i], lat_data_array[i], lon_data_array[i], lat_data_array[i]);    
}


void Remap_operator_grid::initialize_for_vertex_coord_values_generation()
{
    int num_leaf_grids;
    Remap_grid_class *leaf_grids[256], *super_grid;
    long i, j, last_index;
    Remap_grid_data_class *current_vertex_field;


    if (grid_vertex_fields.size() > 0)
        return;

    if (remap_grid->get_is_sphere_grid() && !remap_grid->check_vertex_fields_completeness(remap_operator)) {
        EXECUTION_REPORT(REPORT_WARNING, -1, false, "the vertex values corresponding to sphere grid \"%s\" will be generated automatically", remap_grid->get_grid_name());
        remap_grid->generate_voronoi_grid();
    }

    remap_grid->get_leaf_grids(&num_leaf_grids, leaf_grids, remap_grid);    
    for (i = 0; i < num_leaf_grids; i ++) {
        if (leaf_grids[i]->grid_vertex_fields.size() == 0)
            current_vertex_field = leaf_grids[i]->get_grid_vertex_field();
        else current_vertex_field = leaf_grids[i]->grid_vertex_fields[0];
        EXECUTION_REPORT(REPORT_ERROR, -1, current_vertex_field != NULL, 
                         "The vertex coordinate values of the grid %s are missing, which are not specified by users or generated automatically",
                         remap_grid->get_grid_name());
        grid_vertex_fields.push_back(current_vertex_field);
    }

    this->num_vertexes = 1;
    for (i = 0; i < num_leaf_grids; i ++) {
        if (leaf_grids[i] == NULL)
            continue;
        super_grid = leaf_grids[i]->get_super_grid_of_setting_coord_values();
        if (super_grid->num_dimensions == 1)
            EXECUTION_REPORT(REPORT_ERROR, -1, super_grid->num_vertexes == 2, "remap software error2 in initialize_for_vertex_coord_values_generation\n");
        this->num_vertexes *= super_grid->num_vertexes;
        vertex_values_groups[num_vertex_values_groups] = new Vertex_values_group();
        vertex_values_groups[num_vertex_values_groups]->num_vertex = super_grid->num_vertexes;
        vertex_values_groups[num_vertex_values_groups]->num_coords = 1;
        vertex_values_groups[num_vertex_values_groups]->coord_value_grid = super_grid;
        vertex_values_groups[num_vertex_values_groups]->vertex_coord_values[0] = (double*) grid_vertex_fields[i]->grid_data_field->data_buf;
        last_index = i;
        for (j = i+1; j < num_leaf_grids; j ++) {
            if (leaf_grids[j] == NULL)
                continue;
            if (leaf_grids[j]->get_super_grid_of_setting_coord_values() == super_grid) {
                EXECUTION_REPORT(REPORT_ERROR, -1, last_index+1 == j, "remap software error3 in initialize_for_vertex_coord_values_generation\n");
                vertex_values_groups[num_vertex_values_groups]->vertex_coord_values[vertex_values_groups[num_vertex_values_groups]->num_coords++] = (double*) grid_vertex_fields[j]->grid_data_field->data_buf;
                last_index = j;
                leaf_grids[j] = NULL;
            }
        }
        num_vertex_values_groups ++;
    }

    for (i = 0; i < num_grid_dimensions; i ++)
        vertex_coord_values[i] = new double [grid_size*num_vertexes];

    num_neighbors = num_vertexes;
}


void Remap_operator_grid::generate_overall_vertex_coord_values()
{
    long whole_grid_cell_index, group_grid_cell_index[256];
    int i, j;


    if (num_vertex_values_groups == 1) {
        for (i = 0; i < vertex_values_groups[0]->num_coords; i ++)
            memcpy(vertex_coord_values[i], vertex_values_groups[0]->vertex_coord_values[i], vertex_values_groups[0]->coord_value_grid->grid_size*vertex_values_groups[0]->num_vertex*sizeof(double));
        return;
    }

    for (i = 0; i < num_vertex_values_groups; i ++)
        group_grid_cell_index[i] = 0;
    for (whole_grid_cell_index = 0; whole_grid_cell_index < grid_size; whole_grid_cell_index ++) {
        vertex_values_groups[0]->generate_vertex_values_of_one_cell_according_to_vertex_values_groups(whole_grid_cell_index,
                                                                                                      group_grid_cell_index, 
                                                                                                      num_vertex_values_groups,
                                                                                                      num_grid_dimensions,
                                                                                                      num_vertexes,
                                                                                                      vertex_values_groups,
                                                                                                      vertex_coord_values);
        group_grid_cell_index[0] ++;
        for (i = 0; i < num_vertex_values_groups; i ++) 
            if (group_grid_cell_index[i] == vertex_values_groups[i]->coord_value_grid->grid_size) {
                group_grid_cell_index[i] = 0;
                group_grid_cell_index[i+1] ++;
            }
            else break;
    }
}


bool Remap_operator_grid::two_cells_have_common_bound(long cell_id1, 
                                                      long cell_id2, 
                                                      double tolerable_error)
{
    int i, j, k, num_common_vertex;

    if (cell_id1 == cell_id2)
        return false;
    
    for (i = 0, num_common_vertex = 0; i < num_vertexes; i ++) {
        if (vertex_coord_values[0][cell_id1*num_vertexes+i] == NULL_COORD_VALUE)
            continue;
        for (j = 0; j < num_vertexes; j ++) {
            if (vertex_coord_values[0][cell_id2*num_vertexes+j] == NULL_COORD_VALUE)
                continue;            
            for (k = 0; k < num_grid_dimensions; k ++)
                if (fabs(vertex_coord_values[k][cell_id1*num_vertexes+i] - vertex_coord_values[k][cell_id2*num_vertexes+j]) > tolerable_error)
                    break;
            if (k == num_grid_dimensions) {
                num_common_vertex ++;
                break;
            }
        }
    }

    EXECUTION_REPORT(REPORT_ERROR, -1, cell_id1 != cell_id2,
                 "remap software error1 in two_cells_have_common_bound\n");

    return num_common_vertex >= num_grid_dimensions;
}


long Remap_operator_grid::search_cell_of_locating_point(double *point_coord_values, bool accurately_match)
{
    long i;

    
    if (this->num_grid_dimensions == 1) {
        for (i = 0; i < grid_size - 1; i ++)
            if (compute_difference_of_two_coord_values(center_coord_values[0][i],point_coord_values[0],0) <= 0 && 
                compute_difference_of_two_coord_values(center_coord_values[0][i+1],point_coord_values[0],0) >= 0 ||
                compute_difference_of_two_coord_values(center_coord_values[0][i],point_coord_values[0],0) >= 0 && 
                compute_difference_of_two_coord_values(center_coord_values[0][i+1],point_coord_values[0],0) <= 0)
                return i;
        if (compute_difference_of_two_coord_values(center_coord_values[0][0],center_coord_values[0][grid_size-1],0) <= 0)
            if (compute_difference_of_two_coord_values(point_coord_values[0],center_coord_values[0][0],0) <= 0)
                return 0;
            else return grid_size - 1;
        else if (compute_difference_of_two_coord_values(point_coord_values[0],center_coord_values[0][0],0) >= 0)
            return 0;
        else return grid_size - 1;
    }
    else return grid2D_search_engine->search_cell_of_locating_point(point_coord_values[0], point_coord_values[1], accurately_match);
}


void Remap_operator_grid::visit_cell(long cell_index)
{
    EXECUTION_REPORT(REPORT_ERROR, -1, !cell_visiting_mark[cell_index],
                 "remap software error in visit_cell");

    cell_visiting_mark[cell_index] = true;
    visited_cells_indexes[num_visited_cells++] = cell_index;
}


void Remap_operator_grid::clear_cell_visiting_info()
{
    for (long i = 0; i < num_visited_cells; i ++)
        cell_visiting_mark[visited_cells_indexes[i]] = false;
    num_visited_cells = 0;
}


