/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file is initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#include "cor_global_data.h"
#include "remap_operator_1D_basis.h"
#include <math.h>


double *common_buffer_for_1D_remap_operator = NULL;
int size_common_buffer_for_1D_remap_operator = 0;


Remap_operator_1D_basis::Remap_operator_1D_basis(const char *object_name, const char *operator_name, int num_remap_grids, Remap_grid_class **remap_grids)
                                       : Remap_operator_basis(object_name, 
                                                              operator_name, 
                                                              1, 
                                                              true, 
                                                              false, 
                                                              false, 
                                                              num_remap_grids, 
                                                              remap_grids)
{
    initialize_1D_remap_operator();
}


Remap_operator_1D_basis::~Remap_operator_1D_basis()
{
}


void Remap_operator_1D_basis::set_common_parameter(const char *parameter_name, const char *parameter_value)
{
    EXECUTION_REPORT(REPORT_ERROR, -1, enable_to_set_parameters, 
                 "the parameter of remap operator object \"%s\" must be set before using it to build remap strategy\n",
                 object_name);
    if (words_are_the_same(parameter_name, "periodic")) {
        EXECUTION_REPORT(REPORT_ERROR, -1, !set_periodic, 
                         "The parameter \"%s\" of the 1D remapping operator \"%s\" has been set before. It can not been set more than once",
                         parameter_name, operator_name);
        if (words_are_the_same(parameter_value, "true"))
            periodic = true;
        else if (words_are_the_same(parameter_value, "false"))
            periodic = false;
        else EXECUTION_REPORT(REPORT_ERROR, -1, false, 
                      "The value of parameter \"%s\" of the 1D remapping operator \"%s\" must be \"true\" or \"false\"",
                      parameter_name, operator_name);
        set_periodic = true;
        if (periodic)
            EXECUTION_REPORT(REPORT_ERROR, -1, !set_enable_extrapolation, 
                             "The parameter \"extrapolation\" of the 1D remapping operator \"%s\" has been set before. This remapping operator can not be set to periodic",
                             parameter_name, operator_name);
    }
    else if (words_are_the_same(parameter_name, "period")) {
        EXECUTION_REPORT(REPORT_ERROR, -1, set_periodic && periodic, 
                         "The 1D remapping operator \"%s\" has not been set to periodic before. Its \"period\" can not be set",
                         operator_name);
        EXECUTION_REPORT(REPORT_ERROR, -1, !set_period,
                         "The parameter \"%s\" of the 1D remapping operator \"%s\" has been set before. It can not been set more than once",
                         parameter_name, operator_name);
        sscanf(parameter_value, "%lf", &period);
        set_period = true;
        EXECUTION_REPORT(REPORT_ERROR, -1, period > 0,
                         "The parameter \"%s\" of the 1D remapping operator \"%s\" must be bigger than 0",
                         parameter_name, operator_name);
    }
    else if (words_are_the_same(parameter_name, "use_logarithmic_coordinate")) {
        EXECUTION_REPORT(REPORT_ERROR, -1, !set_use_logarithmic_coordinate,
                         "The parameter \"%s\" of the 1D spline remapping operator \"%s\" has been set before. It can not been set more than once",
                         parameter_name, operator_name);
        if (words_are_the_same(parameter_value, "true")) 
            use_logarithmic_coordinate = true;
        else if (words_are_the_same(parameter_value, "false"))
            use_logarithmic_coordinate = false;
        else EXECUTION_REPORT(REPORT_ERROR, -1, false, 
                      "The value of parameter \"%s\" of the 1D spline remapping operator \"%s\" must be \"none\", \"overall\" or \"fragment\"",
                      parameter_name, operator_name);
        set_use_logarithmic_coordinate = true;
    }
    else if (words_are_the_same(parameter_name, "use_exponent_coordinate")) {
        EXECUTION_REPORT(REPORT_ERROR, -1, !set_use_exponent_coordinate,
                         "The parameter \"%s\" of the 1D spline remapping operator \"%s\" has been set before. It can not been set more than once",
                         parameter_name, operator_name);
        if (words_are_the_same(parameter_value, "true")) 
            use_exponent_coordinate = true;
        else if (words_are_the_same(parameter_value, "false"))
            use_exponent_coordinate = false;
        else EXECUTION_REPORT(REPORT_ERROR, -1, false, 
                      "The value of parameter \"%s\" of the 1D spline remapping operator \"%s\" must be \"none\", \"overall\" or \"fragment\"",
                      parameter_name, operator_name);
        set_use_exponent_coordinate = true;
    }
    else if (words_are_the_same(parameter_name, "enable_extrapolate")) {
        EXECUTION_REPORT(REPORT_ERROR, -1, !periodic,
                         "The parameter \"%s\" of the 1D remapping operator \"%s\" can not be set when the 1D remapping operator is periodic",
                         parameter_name, operator_name);
        EXECUTION_REPORT(REPORT_ERROR, -1, !set_enable_extrapolation,
                         "Extrapolation of remapping operator %s has been enabled before. Extrapolation cannot be enabled again", operator_name);
        if (words_are_the_same(parameter_value, "true"))
            enable_extrapolate = true;
        else if (words_are_the_same(parameter_value, "false"))
            enable_extrapolate = false;
        else EXECUTION_REPORT(REPORT_ERROR, -1, false, 
                      "The value of parameter \"%s\" of the 1D remapping operator \"%s\" must be \"true\" or \"false\"",
                      parameter_name, operator_name);        
        set_enable_extrapolation = true;
    }
    else if (words_are_the_same(parameter_name, "extrapolation_approach")) {
        EXECUTION_REPORT(REPORT_ERROR, -1, set_enable_extrapolation,
                         "Extrapolation of remapping operator %s is not set to enabled. Therefore, the approach of extrapolation cannot be set", operator_name);
        EXECUTION_REPORT(REPORT_ERROR, -1, extrapolation_approach == 0,
                         "Extrapolation approach of remapping operator %s has been set before. It cannot be set again", operator_name);
        if (words_are_the_same(parameter_value, "with_boundary"))
            extrapolation_approach = 1;
        else EXECUTION_REPORT(REPORT_ERROR, -1, false, 
                      "The value of parameter \"%s\" of the 1D remapping operator \"%s\" must be \"with_boundary\" currently",
                      parameter_name, operator_name);        
        set_enable_extrapolation = true;
    }
    else EXECUTION_REPORT(REPORT_ERROR, -1, false, 
                          "\"%s\" is a illegal parameter of remapping operator \"%s\"\n",
                          parameter_name, operator_name);    
    EXECUTION_REPORT(REPORT_ERROR, -1, !use_logarithmic_coordinate || !use_exponent_coordinate, "The parameter \"use_logarithmic_coordinate\" and \"use_exponent_coordinate\" of the 1D spline remapping operator \"%s\" cannot be set to true at the same time. Please verify.", operator_name);
}


int Remap_operator_1D_basis::check_common_parameter(const char *parameter_name, const char *parameter_value, char *error_string)
{
    int check_result = 0;

    
    if (words_are_the_same(parameter_name, "periodic") || words_are_the_same(parameter_name, "use_exponent_coordinate") || 
        words_are_the_same(parameter_name, "use_logarithmic_coordinate") || words_are_the_same(parameter_name, "enable_extrapolate")) {
        check_result = 1;
        if (words_are_the_same(parameter_value, "true") || words_are_the_same(parameter_value, "false"))
            check_result = 3;
        else sprintf(error_string, "The parameter value must be \"true\" or \"false\""); 
    }
    else if (words_are_the_same(parameter_name, "period")) {
        check_result = 1;
        sscanf(parameter_value, "%lf", &period);
        if (period > 0)
            check_result = 3;
        else sprintf(error_string, "The parameter value must be a positive value"); 
        
    }
    else if (words_are_the_same(parameter_name, "extrapolation_approach")) {
        check_result = 1;
        if (words_are_the_same(parameter_name, "with_boundary"))
            check_result = 3;
        else sprintf(error_string, "The parameter value must be \"with_boundary\""); 
    }

    return check_result;
}



void Remap_operator_1D_basis::allocate_1D_remap_operator_common_arrays_space()
{
    int required_size = (src_grid->get_grid_size()+2)*12+(dst_grid->get_grid_size()+2)*12;

    if (size_common_buffer_for_1D_remap_operator < required_size) {
        if (common_buffer_for_1D_remap_operator != NULL)
            delete [] common_buffer_for_1D_remap_operator;
        size_common_buffer_for_1D_remap_operator = required_size;
        common_buffer_for_1D_remap_operator = new double [size_common_buffer_for_1D_remap_operator];
    }
    
    coord_values_src = common_buffer_for_1D_remap_operator + 0*(src_grid->get_grid_size()+2);
    packed_data_values_src = common_buffer_for_1D_remap_operator + 1*(src_grid->get_grid_size()+2);
    useful_src_cells_global_index = (int*) (common_buffer_for_1D_remap_operator + 2*(src_grid->get_grid_size()+2));
    coord_values_dst = common_buffer_for_1D_remap_operator + 12*(src_grid->get_grid_size()+2) + 0*(dst_grid->get_grid_size()+2);
    src_cell_index_left = (int*) (common_buffer_for_1D_remap_operator + 12*(src_grid->get_grid_size()+2) + 1*(dst_grid->get_grid_size()+2));
    src_cell_index_right = (int*) (common_buffer_for_1D_remap_operator + 12*(src_grid->get_grid_size()+2) + 2*(dst_grid->get_grid_size()+2));
}


void Remap_operator_1D_basis::copy_1D_remap_operator_info(Remap_operator_1D_basis *original_grid)
{
    periodic = original_grid->periodic;
    set_periodic = original_grid->set_periodic;
    set_period = original_grid->set_period;
    set_enable_extrapolation = original_grid->set_enable_extrapolation;
    enable_extrapolate = original_grid->enable_extrapolate;
    extrapolation_approach = original_grid->extrapolation_approach;
    period = original_grid->period;
    num_useful_src_cells = original_grid->num_useful_src_cells;
    use_logarithmic_coordinate = original_grid->use_logarithmic_coordinate;
    set_use_logarithmic_coordinate = original_grid->set_use_logarithmic_coordinate;
    use_exponent_coordinate = original_grid->use_exponent_coordinate;
    set_use_exponent_coordinate = original_grid->set_use_exponent_coordinate;
}


void Remap_operator_1D_basis::initialize_1D_remap_operator()
{
    set_periodic = false;
    set_period = false;
    periodic = false;
    enable_extrapolate = false;
    set_enable_extrapolation = false;
    use_logarithmic_coordinate = false;
    set_use_logarithmic_coordinate = false;
    use_exponent_coordinate = false;
    set_use_exponent_coordinate = false;
    extrapolation_approach = 0;
    period = -999;

    allocate_1D_remap_operator_common_arrays_space();
}


void Remap_operator_1D_basis::search_src_cells_around_dst_cell(double coord_value_dst, int src_index_start, int src_index_end, int &src_cell_index_left, int &src_cell_index_right)
{
    int src_index_mid;


    if (coord_values_src[src_index_start] <= coord_values_src[src_index_end]) {
        if (coord_values_src[src_index_start] > coord_value_dst) {
            EXECUTION_REPORT(REPORT_ERROR, -1, !periodic, "software error1: can not find the location of dst cell in original grid");
            src_cell_index_left = -1;
            src_cell_index_right = src_index_start;
            return;
        }
        if (coord_values_src[src_index_end] < coord_value_dst) {
            EXECUTION_REPORT(REPORT_ERROR, -1, !periodic, "software error2: can not find the location of dst cell in original grid");
            src_cell_index_left = src_index_end;
            src_cell_index_right = -1;
            return;
        }
    }
    else {
        if (coord_values_src[src_index_start] < coord_value_dst) {
            EXECUTION_REPORT(REPORT_ERROR, -1, !periodic, "software error1: can not find the location of dst cell in original grid");
            src_cell_index_left = -1;
            src_cell_index_right = src_index_start;
            return;
        }
        if (coord_values_src[src_index_end] > coord_value_dst) {
            EXECUTION_REPORT(REPORT_ERROR, -1, !periodic, "software error2: can not find the location of dst cell in original grid");
            src_cell_index_left = src_index_end;
            src_cell_index_right = -1;
            return;
        }
    }

    if (coord_values_src[src_index_start] == coord_value_dst) {
        src_cell_index_left = src_index_start;
        src_cell_index_right = src_index_start + 1;
        return;
    }

    if (coord_values_src[src_index_end] == coord_value_dst) {
        src_cell_index_left = src_index_end-1;
        src_cell_index_right = src_index_end;
        return;
    }

    search_src_cells_around_dst_cell_recursively(coord_value_dst, src_index_start, src_index_end, src_cell_index_left, src_cell_index_right);
}


void Remap_operator_1D_basis::search_src_cells_around_dst_cell_recursively(double coord_value_dst, int src_index_start, int src_index_end, int &src_cell_index_left, int &src_cell_index_right)
{
    int src_index_mid;


    if (src_index_start == src_index_end-1) {
        src_cell_index_left = src_index_start;
        src_cell_index_right = src_index_end;
        EXECUTION_REPORT(REPORT_ERROR, -1, coord_values_src[src_cell_index_left] <= coord_value_dst && coord_values_src[src_cell_index_right] >= coord_value_dst, 
                         "software error: can not find the location of dst cell in original grid recursively");
        return;
    }
    
    src_index_mid = (src_index_start+src_index_end) / 2;
    if (coord_values_src[src_index_start] <= coord_values_src[src_index_end]) {
        if (coord_values_src[src_index_mid] > coord_value_dst)
            search_src_cells_around_dst_cell_recursively(coord_value_dst, src_index_start, src_index_mid, src_cell_index_left, src_cell_index_right);
        else search_src_cells_around_dst_cell_recursively(coord_value_dst, src_index_mid, src_index_end, src_cell_index_left, src_cell_index_right);
    }
    else {
        if (coord_values_src[src_index_mid] <= coord_value_dst)
            search_src_cells_around_dst_cell_recursively(coord_value_dst, src_index_start, src_index_mid, src_cell_index_left, src_cell_index_right);
        else search_src_cells_around_dst_cell_recursively(coord_value_dst, src_index_mid, src_index_end, src_cell_index_left, src_cell_index_right);        
    }
}


void Remap_operator_1D_basis::calculate_dst_src_mapping_info()
{
    int i, j;
    bool ascending_order, src_cell_mask, dst_cell_mask;
    

    array_size_src = 0;

    for (i = 0; i < dst_grid->get_grid_size(); i ++)
        get_cell_center_coord_values_of_dst_grid(i, &coord_values_dst[i]);
    for (i = 0; i < src_grid->get_grid_size(); i ++)
        get_cell_center_coord_values_of_src_grid(i, &coord_values_src[i]);

    ascending_order = coord_values_src[0] < coord_values_src[1];
    for (i = 1; i < src_grid->get_grid_size() - 1; i ++) 
        EXECUTION_REPORT(REPORT_ERROR, -1,ascending_order == coord_values_src[i] < coord_values_src[i+1] || coord_values_src[i] == coord_values_src[i+1], 
                         "the center coordinate values corresponding to the 1D grid %s are not sorted into ascending or descending order",
                         src_grid->get_grid_name());
    ascending_order = coord_values_dst[0] < coord_values_dst[1];

    for (i = 1; i < dst_grid->get_grid_size() - 1; i ++) 
        EXECUTION_REPORT(REPORT_ERROR, -1, ascending_order == coord_values_dst[i] < coord_values_dst[i+1] || coord_values_dst[i] == coord_values_dst[i+1], 
                         "the center coordinate values corresponding to the 1D grid %s are not sorted into ascending or descending order",
                         dst_grid->get_grid_name());

    if (periodic) {
        EXECUTION_REPORT(REPORT_ERROR, -1,fabs(coord_values_src[0]-coord_values_dst[src_grid->get_grid_size()-1]) < period, 
                         "The variation of center coordinate values corresponding to the 1D grid %s is larger than one period",
                         src_grid->get_grid_name());
        EXECUTION_REPORT(REPORT_ERROR, -1,fabs(coord_values_dst[0]-coord_values_dst[dst_grid->get_grid_size()-1]) < period, 
                         "The variation of center coordinate values corresponding to the 1D grid %s is larger than one period",
                         dst_grid->get_grid_name());
    }
    
    num_useful_src_cells = 0;
    ascending_order = coord_values_src[0] < coord_values_src[1];
    if (ascending_order) {
        for (i = 0; i < src_grid->get_grid_size(); i ++) {
            get_cell_mask_of_src_grid(i, &src_cell_mask);
            if (src_cell_mask) {
                get_cell_center_coord_values_of_src_grid(i, &coord_values_src[num_useful_src_cells]);
                useful_src_cells_global_index[num_useful_src_cells++] = i;
            }
        }    
    }
    else {
        for (i = src_grid->get_grid_size()-1; i >= 0; i --) {
            get_cell_mask_of_src_grid(i, &src_cell_mask);
            if (src_cell_mask) {
                get_cell_center_coord_values_of_src_grid(i, &coord_values_src[num_useful_src_cells]);
                useful_src_cells_global_index[num_useful_src_cells++] = i;
            }
        }
    }

    if (num_useful_src_cells == 0)
        return;
    
    if (periodic) {
        coord_values_src[num_useful_src_cells] = coord_values_src[0] + period;
        useful_src_cells_global_index[num_useful_src_cells] = useful_src_cells_global_index[0];
        array_size_src = num_useful_src_cells + 1;
    }
    else array_size_src = num_useful_src_cells;

    if (use_logarithmic_coordinate) {
        for (i = 0; i < array_size_src; i ++) {
            if (coord_values_src[i] <= 0)
                EXECUTION_REPORT(REPORT_ERROR, -1, coord_values_src[i] > 0, "1D remapping operator %s cannot use logarithmic values of coodinate as new coordinate values because its source grid %s has non-positive coordinate values",
                                 object_name, src_grid->get_grid_name());
            coord_values_src[i] = log(coord_values_src[i]);
        }
        for (i = 0; i < dst_grid->get_grid_size(); i ++) {
            if (coord_values_dst[i] <= 0)
                EXECUTION_REPORT(REPORT_ERROR, -1, coord_values_dst[i] > 0, "1D remapping operator %s cannot use logarithmic values of coodinate as new coordinate values because its source grid %s has non-positive coordinate values",
                                 object_name, dst_grid->get_grid_name());
            coord_values_dst[i] = log(coord_values_dst[i]);
        }            
    }

    if (use_exponent_coordinate) {
        for (i = 0; i < array_size_src; i ++) {
            if (coord_values_src[i] > 0)
                EXECUTION_REPORT(REPORT_ERROR, -1, coord_values_src[i] <= 0, "1D remapping operator %s cannot use exponent values of coodinate as new coordinate values because its source grid %s has positive coordinate values",
                                 object_name, src_grid->get_grid_name());
            coord_values_src[i] = exp(coord_values_src[i]);
        }
        for (i = 0; i < dst_grid->get_grid_size(); i ++) {
            if (coord_values_dst[i] > 0)
                EXECUTION_REPORT(REPORT_ERROR, -1, coord_values_dst[i] <= 0, "1D remapping operator %s cannot use exponent values of coodinate as new coordinate values because its source grid %s has positive coordinate values",
                                 object_name, dst_grid->get_grid_name());
            coord_values_dst[i] = exp(coord_values_dst[i]);
        }            
    }

    for (i = 0; i < dst_grid->get_grid_size(); i ++) {
        src_cell_index_left[i] = -1;
        src_cell_index_right[i] = -1;
        get_cell_mask_of_dst_grid(i, &dst_cell_mask);
        if (!dst_cell_mask)
            continue;
        search_src_cells_around_dst_cell(coord_values_dst[i], 0, array_size_src-1, src_cell_index_left[i], src_cell_index_right[i]);
        for (j = 0; j < array_size_src && coord_values_src[j] <= coord_values_dst[i]; j ++);
        if (j > 0 && j < array_size_src)
            EXECUTION_REPORT(REPORT_ERROR, -1, src_cell_index_left[i] == j-1 && src_cell_index_right[i] == j, "error error3\n"); 
        if ((src_cell_index_left[i] == -1 || src_cell_index_right[i] == -1) && !enable_extrapolate)
            continue;
        if (src_cell_index_right[i] == -1) {
            src_cell_index_right[i] = src_cell_index_left[i];
            if (extrapolation_approach == 0)
                src_cell_index_left[i] --;
            else if (extrapolation_approach == 1)
                src_cell_index_left[i] = src_cell_index_left[i];
        }
        else if (src_cell_index_left[i] == -1) {
            src_cell_index_left[i] = src_cell_index_right[i];
            if (extrapolation_approach == 0)
                src_cell_index_right[i] ++;
            else if (extrapolation_approach == 1)
                src_cell_index_right[i] = src_cell_index_right[i];
        }
    }
}


void Remap_operator_1D_basis::preprocess_field_value(double *data_values_src)
{
    for (int i = 0; i < array_size_src; i ++)
        packed_data_values_src[i] = data_values_src[useful_src_cells_global_index[i]];
}


void Remap_operator_1D_basis::postprocess_field_value(double *data_values_dst)
{
    int i, dst_index;    
}

