/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#include "cor_global_data.h"
#include "remap_grid_mgt.h"
#include "parse_special_words.h"
#include <string.h>
#include <stdio.h>


void Remap_grid_mgt::execute(const char*function, Remap_statement_operand **statement_operands, int num_operands)
{
    int i, j;
    long grid_size;
    int has_grid_size, dimension_para_index;
    int num_sub_grids;
    int num_vertexes;
    Remap_grid_class *sub_grids[256];
    char cyclic_label[256];
    double unit_conversion;


    /* Check the semantics of function "new_1D_grid" and then initialize a 1D grid */
    if (words_are_the_same(function, FUNCTION_WORD_NEW_1D_GRID)) {
        EXECUTION_REPORT(REPORT_ERROR, -1, num_operands >= 3 && num_operands <= 5, 
                     "function \"%s\" must have one result parameter and two to four input parameters\n", 
                     function);
        check_is_parameter_object_type_grid(function, 0, statement_operands[0], "the 1D grid to be generated");
        check_is_parameter_string_type(function, 1, statement_operands[1], "the coordinate label");
        check_is_parameter_string_type(function, 2, statement_operands[2], "the unit of coordinate");
        cyclic_label[0] = '\0';
        if (words_are_the_same(statement_operands[1]->extension_names[0], COORD_LABEL_LON)) {
            EXECUTION_REPORT(REPORT_ERROR, -1, num_operands == 4 || num_operands == 5 , "function \"%s\" for generating 1D grid of lon must have one result parameter and three to four input parameters\n", function);
            check_is_parameter_string_type(function, 3, statement_operands[3], "the label of cyclic or acyclic feature of lon coordinate");
            strcpy(cyclic_label, statement_operands[3]->extension_names[0]);
            dimension_para_index = 4;
        }
        else dimension_para_index = 3;
        grid_size = 0;
        if (num_operands == dimension_para_index + 1)
            grid_size = get_size_value_from_parameter(function, dimension_para_index, statement_operands[dimension_para_index], "the size of 1D grid");
        remap_grids.push_back(new Remap_grid_class(statement_operands[0]->object->object_name, 
                                                   statement_operands[1]->extension_names[0],
                                                   statement_operands[2]->extension_names[0],
                                                   cyclic_label,
                                                   grid_size));
    }
    else if (words_are_the_same(function, FUNCTION_WORD_EXTRACT_MASK)) {
        EXECUTION_REPORT(REPORT_ERROR, -1, num_operands == 4, "function \"%s\" must have one result parameter and three input parameters\n", function);
        check_is_parameter_grid_mask_field(function, statement_operands[0], NULL);
        check_is_parameter_object_type_field_data(function, 1, statement_operands[1], "the field with missing values to extract mask");
        for (i = 2; i <= 3; i ++)
            get_float_value_from_parameter(function, i, statement_operands[i], "one boundary of filling value");
        remap_grid_manager->search_remap_grid_with_grid_name(statement_operands[0]->object->object_name)->extract_mask(statement_operands[1]->object->object_name,
                                                                                                                       statement_operands[2]->extension_names[0],
                                                                                                                       statement_operands[3]->extension_names[0]);
    }
    else if (words_are_the_same(function, FUNCTION_WORD_COMPUTE_OCN_MASK)) {
        EXECUTION_REPORT(REPORT_ERROR, -1, num_operands == 3, "function \"%s\" must have one result parameter and two input parameters\n", function);
        check_is_parameter_grid_mask_field(function, statement_operands[0], NULL);
        check_is_parameter_object_type_field_data(function, 1, statement_operands[1], "the topo field");
        unit_conversion = get_float_value_from_parameter(function, 2, statement_operands[2], "the unit conversion for terrain data\n");
        remap_grid_manager->search_remap_grid_with_grid_name(statement_operands[0]->object->object_name)->compute_ocn_mask(statement_operands[1]->object->object_name,
                                                                                                                           unit_conversion);        
    }
    /* Check the semantics of function "combine_coord_systems" and then initialize a multi-dimension grid 
        according to sub grids. The new grid can not be the same with any other grid */
    else if (words_are_the_same(function, FUNCTION_WORD_COMBINE_GRIDS)) {
        EXECUTION_REPORT(REPORT_ERROR, -1, num_operands >= 3, "function \"%s\" must have one result parameter and at least three input parameters\n", function);
        check_is_parameter_object_type_grid(function, 0, statement_operands[0], "the multi-dimension grid to be generated");
        grid_size = 0;
        has_grid_size = 0;
        if (statement_operands[num_operands-1]->object == NULL || !words_are_the_same(statement_operands[num_operands-1]->object->object_type, OBJECT_TYPE_GRID)) {
            grid_size = get_size_value_from_parameter(function, num_operands - 1, statement_operands[num_operands-1], "the size of multi-dimension grid");
            has_grid_size = 1;
        }
        for (i = 1, num_sub_grids = 0; i < num_operands - has_grid_size; i ++) {
            check_is_parameter_object_type_grid(function, i, statement_operands[i], "one sub grid");
            sub_grids[num_sub_grids++] = search_remap_grid_with_grid_name(statement_operands[i]->object->object_name);
        }
        remap_grids.push_back(new Remap_grid_class(statement_operands[0]->object->object_name,
                                                   num_sub_grids,
                                                   sub_grids,
                                                   grid_size));
    }
    else if (words_are_the_same(function, FUNCIION_WORD_LEV_COORD_FROM_SIGMA)) {
        double scale_factor = 1.0;
        EXECUTION_REPORT(REPORT_ERROR, -1, num_operands == 4 || num_operands == 5, "function \"%s\" must have one result parameter and three or four input parameters\n", function);
        check_is_parameter_grid_center_field(function, statement_operands[0], NULL);
        check_is_parameter_object_type_field_data(function, 1, statement_operands[1], "the field with values at the bottom of 3D grid");
        get_float_value_from_parameter(function, 2, statement_operands[2], "the value corresponding to 3D grid top\n");
        check_is_parameter_object_type_field_data(function, 3, statement_operands[3], "the field with sigma coordinate");
        if (num_operands == 5)
            scale_factor = get_float_value_from_parameter(function, 4, statement_operands[4], "the scale factor value\n");
        search_remap_grid_with_grid_name(statement_operands[0]->object->object_name)->gen_lev_coord_from_sigma_or_hybrid(statement_operands[0]->extension_names, 
                                                                                                                         statement_operands[1]->object->object_name, 
                                                                                                                         statement_operands[2]->extension_names[0],
                                                                                                                         statement_operands[3]->object->object_name,
                                                                                                                         NULL, scale_factor);
    }
    else if (words_are_the_same(function, FUNCIION_WORD_LEV_COORD_FROM_HYBRID)) {
        double scale_factor = 1.0;
        EXECUTION_REPORT(REPORT_ERROR, -1, num_operands == 5 || num_operands == 6, "function \"%s\" must have one result parameter and four or five input parameters\n", function);
        check_is_parameter_grid_center_field(function, statement_operands[0], NULL);
        check_is_parameter_object_type_field_data(function, 1, statement_operands[1], "the field with values at the bottom of 3D grid");
        get_float_value_from_parameter(function, 2, statement_operands[2], "the constant reference value corresponding to hybrid grid\n");
        check_is_parameter_object_type_field_data(function, 3, statement_operands[3], "the field with sigma coordinate");
        check_is_parameter_object_type_field_data(function, 4, statement_operands[4], "the field with hybrid grid coefficient");
        if (num_operands == 6)
            scale_factor = get_float_value_from_parameter(function, 5, statement_operands[5], "the scale factor value\n");
        search_remap_grid_with_grid_name(statement_operands[0]->object->object_name)->gen_lev_coord_from_sigma_or_hybrid(statement_operands[0]->extension_names, 
                                                                                                                         statement_operands[1]->object->object_name, 
                                                                                                                         statement_operands[2]->extension_names[0],
                                                                                                                         statement_operands[3]->object->object_name,
                                                                                                                         statement_operands[4]->object->object_name,
                                                                                                                         scale_factor);
    }
    /* Check the semantics of function "read_field", and then use IO object to read the data */
    else if (words_are_the_same(function, FUNCTION_WORD_READ_FIELD)) {
        EXECUTION_REPORT(REPORT_ERROR, -1, num_operands == 3 || num_operands == 4, "function \"%s\" must have one result parameter and two or three input parameters\n", function);
        check_is_parameter_grid_field(function, statement_operands[0], NULL);
        check_is_parameter_object_type_IO(function, 1, statement_operands[1], "the IO file to read the grid field");
        check_is_parameter_string_type(function, 2, statement_operands[2], "the variable name of field in IO fiel");
        if (words_are_the_same(statement_operands[0]->extension_names[0], GRID_VERTEX_LABEL)) {
            EXECUTION_REPORT(REPORT_ERROR, -1, num_operands == 4, 
                         "when reading vertex coordinate values, function \"%s\" must have three input parameters where the last one specifies the maximum number of vertexes of each grid cell\n", 
                         function);
            num_vertexes = get_size_value_from_parameter(function, 3, statement_operands[3], "the number of vertexes of each grid cell");
        }
        else {
            EXECUTION_REPORT(REPORT_ERROR, -1, num_operands == 3, 
                         "when reading center coordinate values or mask values, function \"%s\" must have two input parameters\n", 
                         function);
            num_vertexes = 0;
        }
        search_remap_grid_with_grid_name(statement_operands[0]->object->object_name)->read_grid_data_from_IO(statement_operands[0]->extension_names, 
                                                                                                             statement_operands[1]->object->object_name, 
                                                                                                             statement_operands[2]->extension_names[0],
                                                                                                             num_vertexes);
    }
    else if (words_are_the_same(function, FUNCTION_WORD_ISPAN) ||
             words_are_the_same(function, FUNCTION_WORD_FSPAN)) {
        EXECUTION_REPORT(REPORT_ERROR, -1, num_operands == 4, "function \"%s\" for reading field data must have one result parameter and three input parameter\n", function);
        check_is_parameter_grid_center_field(function, statement_operands[0], NULL);
        if (words_are_the_same(function, FUNCTION_WORD_ISPAN))
            for (i = 1; i < 3; i ++)
                get_int_value_from_parameter(function, i, statement_operands[i], "one bound of span");
        else for (i = 1; i < 3; i ++)
                 get_float_value_from_parameter(function, i, statement_operands[i], "one bound of span");
        long field_size = get_size_value_from_parameter(function, 3, statement_operands[3], "the size of field");   
        search_remap_grid_with_grid_name(statement_operands[0]->object->object_name)->read_grid_data_through_span(statement_operands[0]->extension_names, 
                                                                                                                  statement_operands[1]->extension_names[0], 
                                                                                                                  statement_operands[2]->extension_names[0],
                                                                                                                  field_size,
                                                                                                                  function);
    }
    else if (words_are_the_same(function, FUNCTION_WORD_NEW_PARTIAL_GRID)) {
        EXECUTION_REPORT(REPORT_ERROR, -1, num_operands == 2, "function \"%s\" can have only one source operand which is the whole grid\n", function);
        check_is_parameter_object_type_grid(function, 0, statement_operands[0], "the partial grid to be generated");
        check_is_parameter_object_type_grid(function, 1, statement_operands[1], "the whole grid");
        remap_grids.push_back(new Remap_grid_class(statement_operands[0]->object->object_name, statement_operands[1]->object->object_name));
    }
    else if (words_are_the_same(function, FUNCTION_WORD_ADD_GRID_AREA)) {
        EXECUTION_REPORT(REPORT_ERROR, -1, num_operands == 2, "function \"%s\" has two input parameters\n", function);
        check_is_parameter_object_type_grid(function, 1, statement_operands[0], "the partial grid");
        check_is_parameter_string_type(function, 2, statement_operands[1], "the name of partial area");
        search_remap_grid_with_grid_name(statement_operands[0]->object->object_name)->add_partial_grid_area(statement_operands[1]->extension_names[0]);
    }
    else if (words_are_the_same(function, FUNCIION_WORD_SET_LEV_GRID_SIGMA_INFO)) {
        char *hybrid_grid_coefficient_field_name = NULL;
        EXECUTION_REPORT(REPORT_ERROR, -1, num_operands == 4 || num_operands == 5, "function \"%s\" has four input parameters\n", function);
        check_is_parameter_object_type_grid(function, 1, statement_operands[0], "the level grid (vertical grid) to be set the sigma information");
        check_is_parameter_object_type_field_data(function, 2, statement_operands[1], "the vector of sigma values");
        double top_value = get_float_value_from_parameter(function, 3, statement_operands[2], "the vertical coordinate value of top layer in the sigma grid");
        double scale_factor = get_float_value_from_parameter(function, 4, statement_operands[3], "the scale factor corresponding to the sigma grid");
        if (num_operands == 5) {
            check_is_parameter_object_type_field_data(function, 5, statement_operands[4], "the vector of coefficients for hybrid grid");
            hybrid_grid_coefficient_field_name = statement_operands[4]->object->object_name;
        }
        search_remap_grid_with_grid_name(statement_operands[0]->object->object_name)->set_lev_grid_sigma_info(statement_operands[1]->object->object_name, top_value, scale_factor, hybrid_grid_coefficient_field_name);
    }
    else if (words_are_the_same(function, FUNCTION_WORD_ADD_AREA_BOUND)) {
        EXECUTION_REPORT(REPORT_ERROR, -1, num_operands == 6, "function \"%s\" has six input parameters\n", function);
        check_is_parameter_object_type_grid(function, 1, statement_operands[0], "the partial grid");
        check_is_parameter_string_type(function, 2, statement_operands[1], "the name of area to set area bound");
        check_is_parameter_object_type_grid(function, 3, statement_operands[2], "the sub grid corresponding to setting area bounds");
        check_is_parameter_string_type(function, 4, statement_operands[3], "the type of area bounds (index or value)");
        for (i = 4; i < 6; i ++) {
            if (words_are_the_same(statement_operands[3]->extension_names[0], PARTIAL_AREA_BOUND_TYPE_INDEX))
                get_int_value_from_parameter(function, i+1, statement_operands[i], "one bound of area");
            else get_float_value_from_parameter(function, i+1, statement_operands[i], "one bound of area");
        }
        search_remap_grid_with_grid_name(statement_operands[0]->object->object_name)->add_partial_grid_area_bounds(statement_operands[1]->extension_names[0],
                                                                                                                   statement_operands[2]->object->object_name,
                                                                                                                   statement_operands[3]->extension_names[0],
                                                                                                                   statement_operands[4]->extension_names[0],
                                                                                                                   statement_operands[5]->extension_names[0]);
    }
    else if (words_are_the_same(function, FUNCTION_WORD_SET_BOUNDARY)) {
        check_is_parameter_grid_boundary(function, statement_operands[0], NULL);
        EXECUTION_REPORT(REPORT_ERROR, -1, num_operands == 5, "function \"%s\" must have one result parameter and four input parameters\n", function);
        double min_lon = get_float_value_from_parameter(function, 1, statement_operands[1], "the minimum longitude of the grid boundary\n");
        double max_lon = get_float_value_from_parameter(function, 2, statement_operands[2], "the maximum longitude of the grid boundary\n");
        double min_lat = get_float_value_from_parameter(function, 3, statement_operands[3], "the minimum latitude of the grid boundary\n");
        double max_lat = get_float_value_from_parameter(function, 4, statement_operands[4], "the maximum latitude of the grid boundary\n");
        search_remap_grid_with_grid_name(statement_operands[0]->object->object_name)->set_grid_boundary(min_lon, max_lon, min_lat, max_lat);
    }
    else EXECUTION_REPORT(REPORT_ERROR, -1, false, "\"%s\" is an unspported function\n", function);
}


Remap_grid_class *Remap_grid_mgt::search_remap_grid_with_grid_name(const char *grid_name)
{
    for (int i = 0; i < remap_grids.size(); i ++)
        if (remap_grids[i]->match_grid(grid_name))
            return remap_grids[i];

    return NULL;
}


Remap_grid_class *Remap_grid_mgt::search_remap_grid_with_coord_name(const char *coord_label)
{
    for (int i = 0; i < remap_grids.size(); i ++)
        if (words_are_the_same(coord_label, remap_grids[i]->get_coord_label()))
            return remap_grids[i];

    return NULL;
}


Remap_grid_class *Remap_grid_mgt::search_remap_grid_with_sized_sub_grids(int num_sized_sub_grids, Remap_grid_class **sized_sub_grids)
{
    for (int i = 0; i < remap_grids.size(); i ++)
        if (remap_grids[i]->match_grid(num_sized_sub_grids, sized_sub_grids))
            return remap_grids[i];

    return NULL;    
}


void Remap_grid_mgt::get_all_leaf_remap_grids(int *num_leaf_grids, Remap_grid_class **leaf_grids)
{
    *num_leaf_grids = 0;
    for (int i = 0; i < remap_grids.size(); i ++)
        if (remap_grids[i]->get_num_dimensions() == 1)
            leaf_grids[(*num_leaf_grids)++] = remap_grids[i];
}


void Remap_grid_mgt::add_remap_grid(Remap_grid_class *remap_grid)
{
    EXECUTION_REPORT(REPORT_ERROR, -1, search_remap_grid_with_grid_name(remap_grid->get_grid_name()) == NULL, "Software error in Remap_grid_mgt::add_remap_grid: grid \"%s\" has already been added into the grid manager", remap_grid->get_grid_name());
    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "add remap grid \"%s\" at %lx", remap_grid->get_grid_name(), remap_grid);
    remap_grids.push_back(remap_grid);
}


void Remap_grid_mgt::add_temp_grid(Remap_grid_class *temp_grid)
{
    temp_grids.push_back(temp_grid);
    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "add temp grid \"%s\" at %lx", temp_grid->get_grid_name(), temp_grid);
}


Remap_grid_mgt::~Remap_grid_mgt()
{
    for (int i = 0; i < remap_grids.size(); i ++)
        delete remap_grids[i];

    for (int i = 0; i < temp_grids.size(); i ++)
        delete temp_grids[i];
}


Remap_grid_class *Remap_grid_mgt::search_same_remap_grid(Remap_grid_class *remap_grid)
{
    for (int i = 0; i < remap_grids.size(); i ++)
        if (remap_grids[i]->is_the_same_grid_with(remap_grid))
            return remap_grids[i];
        
    return NULL;
}

