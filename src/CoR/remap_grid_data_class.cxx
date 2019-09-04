/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#include "cor_global_data.h"
#include "remap_grid_class.h"
#include "remap_common_utils.h"
#include "remap_grid_data_class.h"
#include <math.h>
#include <string.h>


template<class T1, class T2> void span_data_array(T1 *spanned_array, T2 bound_start, T2 bound_end, long num_elements)
{
    T2 stride = (bound_end-bound_start) / (num_elements-1);
    for (long i = 0; i < num_elements; i ++) 
        spanned_array[i] = bound_start + stride*i;
    spanned_array[num_elements-1] = bound_end;    
}


/* Constructor of Remap_grid_data_class. coord_value_grid is the grid of setting the 
    corresponding coordinate values */
Remap_grid_data_class::Remap_grid_data_class(Remap_grid_class *coord_value_grid, 
                                             Remap_data_field *grid_data_field)
{
    this->grid_data_field = grid_data_field;
    generate_grid_info(coord_value_grid);
	if (coord_value_grid != NULL)
		EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "allocate field data regarding to grid \"%s\" with size %ld", coord_value_grid->get_grid_name(), grid_data_field->required_data_size*get_data_type_size(grid_data_field->data_type_in_application));
}


/* Constructor of Remap_grid_data_class for reading IO data */
Remap_grid_data_class::Remap_grid_data_class(const char *field_data_name, 
                                             Remap_grid_class *associative_grid, 
                                             IO_basis *IO_object,
                                             const char *field_name_in_IO)
{
    generate_data_field_info(field_data_name, associative_grid);
    generate_grid_info(associative_grid);
    strcpy(grid_data_field->field_name_in_IO_file, field_name_in_IO);
    if (IO_object != NULL)
        IO_object->read_data(grid_data_field, -1, true);
    else {
        memset(grid_data_field->data_buf, 0, grid_data_field->required_data_size*get_data_type_size(grid_data_field->data_type_in_application));
        strcpy(this->grid_data_field->data_type_in_IO_file, DATA_TYPE_DOUBLE);
    }
}


/* Constructor of Remap_grid_data_class for spanning data */
Remap_grid_data_class::Remap_grid_data_class(const char *field_data_name,
                                             const char *grid_name, 
                                             const char *bound_start, 
                                             const char *bound_end, 
                                             long num_value_points, 
                                             const char *span_data_type,
                                             const char *grid_data_type)
{
    double bound_start_value, bound_end_value;


    sscanf(bound_start, "%lf", &bound_start_value);
    sscanf(bound_end, "%lf", &bound_end_value);

    generate_data_field_info(field_data_name, NULL);
    strcpy(grid_data_field->data_type_in_application, grid_data_type);
    strcpy(grid_data_field->data_type_in_IO_file, grid_data_type);
    grid_data_field->data_buf = new char[num_value_points*get_data_type_size(grid_data_type)];
    grid_data_field->required_data_size = num_value_points;
    grid_data_field->read_data_size = num_value_points;
    coord_value_grid = NULL;
    if (grid_name != NULL) {
        EXECUTION_REPORT(REPORT_ERROR, -1, num_value_points == remap_grid_manager->search_remap_grid_with_grid_name(grid_name)->get_grid_size(),
                     "the size of span array must be the same as the size of grid \"%s\"\n",
                     remap_grid_manager->search_remap_grid_with_grid_name(grid_name)->get_grid_name());
        EXECUTION_REPORT(REPORT_ERROR, -1, words_are_the_same(remap_grid_manager->search_remap_grid_with_grid_name(grid_name)->get_coord_label(), field_data_name),
                     "remap software error1 in Remap_grid_data_class for spanning data\n");
        generate_grid_info(remap_grid_manager->search_remap_grid_with_grid_name(grid_name));
    }
    
    if (words_are_the_same(span_data_type, DATA_TYPE_LONG)) {
        EXECUTION_REPORT(REPORT_ERROR, -1, (long) bound_start_value == bound_start_value, "the first input parameter of ispan function must be an integer\n");
        EXECUTION_REPORT(REPORT_ERROR, -1, (long) bound_end_value == bound_end_value, "the second input parameter of ispan function must be an integer\n");
        EXECUTION_REPORT(REPORT_ERROR, -1, ((long)bound_end_value - (long)bound_start_value)%(num_value_points-1) == 0, "for ispan function, the number of elements in span array does not match the start and end boundaries\n");
    }
    else EXECUTION_REPORT(REPORT_ERROR, -1, words_are_the_same(span_data_type, DATA_TYPE_DOUBLE), "remap software error2 in Remap_grid_data_class for spanning data\n");

    if (words_are_the_same(grid_data_type, DATA_TYPE_LONG) && words_are_the_same(span_data_type, DATA_TYPE_LONG))
        span_data_array((long*) grid_data_field->data_buf, (long) bound_start_value, (long) bound_end_value, num_value_points);
    else if (words_are_the_same(grid_data_type, DATA_TYPE_DOUBLE) && words_are_the_same(span_data_type, DATA_TYPE_LONG))
        span_data_array((double*) grid_data_field->data_buf, (long) bound_start_value, (long) bound_end_value, num_value_points);
    else if (words_are_the_same(grid_data_type, DATA_TYPE_DOUBLE) && words_are_the_same(span_data_type, DATA_TYPE_DOUBLE))
        span_data_array((double*) grid_data_field->data_buf, (double) bound_start_value, (double) bound_end_value, num_value_points);
    else EXECUTION_REPORT(REPORT_ERROR, -1, false, "remap software error3 in Remap_grid_data_class for spanning data\n");
}


void Remap_grid_data_class::generate_grid_info(Remap_grid_class *coord_value_grid)
{
    int num_sized_grids;
    Remap_grid_class *sized_grids[256];


    this->coord_value_grid = coord_value_grid;
    if (coord_value_grid != NULL) {
        coord_value_grid->get_sized_sub_grids(&num_sized_grids, sized_grids);
        reset_sized_grids(num_sized_grids, sized_grids);
        EXECUTION_REPORT(REPORT_ERROR, -1, this->sized_grids.size() > 0, "remap software error in function Remap_grid_data_class\n");
    }
}


void Remap_grid_data_class::generate_data_field_info(const char *field_data_name, 
                                               Remap_grid_class *associative_grid)
{
    grid_data_field = new Remap_data_field;

    strcpy(grid_data_field->field_name_in_application, field_data_name);
    strcpy(grid_data_field->data_type_in_application, DATA_TYPE_DOUBLE);
    grid_data_field->read_data_size = 0;
    if (associative_grid != NULL) {
        grid_data_field->required_data_size = associative_grid->get_grid_size();
        grid_data_field->data_buf = new char [grid_data_field->required_data_size*get_data_type_size(grid_data_field->data_type_in_application)];
        for (long i = 0; i < grid_data_field->required_data_size; i ++)
            ((double*)grid_data_field->data_buf)[i] = DEFAULT_FILL_VALUE;
    }
    else {
        grid_data_field->required_data_size = 0;        
        grid_data_field->data_buf = NULL;
    }
}


void Remap_grid_data_class::set_masked_cell_to_missing_value()
{
    if (coord_value_grid == NULL || coord_value_grid->get_grid_mask_field() == NULL)
        return;

    if (!words_are_the_same(grid_data_field->data_type_in_application, DATA_TYPE_FLOAT) && !words_are_the_same(grid_data_field->data_type_in_application, DATA_TYPE_DOUBLE))
        return;

	coord_value_grid->get_grid_mask_field()->interchange_grid_data(coord_value_grid->get_grid_mask_field()->get_coord_value_grid());
	this->interchange_grid_data(coord_value_grid->get_grid_mask_field()->get_coord_value_grid());

	long mask_size = coord_value_grid->get_grid_mask_field()->grid_data_field->required_data_size;
    bool *mask = (bool*) coord_value_grid->get_grid_mask_field()->grid_data_field->data_buf;
    if (words_are_the_same(grid_data_field->data_type_in_application, DATA_TYPE_FLOAT)) {
		for (int j = 0; j < grid_data_field->required_data_size/mask_size; j ++) {
	        float *data_buf = ((float*) grid_data_field->data_buf) + mask_size*j;
	        for (int i = 0; i < mask_size; i ++)
    	        if (!mask[i])
        	        data_buf[i] = (float) DEFAULT_FILL_VALUE;
		}
    }
    else {
		for (int j = 0; j < grid_data_field->required_data_size/mask_size; j ++) {
	        double *data_buf = ((double*) grid_data_field->data_buf) + mask_size*j;
	        for (int i = 0; i < mask_size; i ++)
	            if (!mask[i])
	                data_buf[i] = (double) DEFAULT_FILL_VALUE;        
		}
    }
}


Remap_grid_data_class::~Remap_grid_data_class()
{
	if (coord_value_grid != NULL)
		EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "deallocate field data regarding to grid \"%s\" with size %ld", coord_value_grid->get_grid_name(), grid_data_field->required_data_size*get_data_type_size(grid_data_field->data_type_in_application));
    delete grid_data_field;
}


Remap_grid_data_class *Remap_grid_data_class::duplicate_grid_data_field(Remap_grid_class *associative_grid, 
                                                                        int num_points_per_cell, 
                                                                        bool copy_data, 
                                                                        bool interchange_data)
{
    Remap_data_field *duplicated_data_field;
    Remap_grid_data_class *duplicated_grid_data_field;


    if (copy_data) {
        if (coord_value_grid != NULL)
            EXECUTION_REPORT(REPORT_ERROR, -1, associative_grid->get_grid_size() == coord_value_grid->get_grid_size(), "remap software error1 in duplicate_grid_data_field\n");
        else {
            EXECUTION_REPORT(REPORT_ERROR, -1, associative_grid->get_grid_size() == grid_data_field->required_data_size, "remap software error1 in duplicate_grid_data_field\n");
            EXECUTION_REPORT(REPORT_ERROR, -1, associative_grid->get_num_dimensions() == 1, "When copy none-grided field %s to a grided field, the grid %s of the grided field must be 1D grid", 
                             grid_data_field->field_name_in_application, associative_grid->get_grid_name());
        }
    }
    EXECUTION_REPORT(REPORT_ERROR, -1, associative_grid->get_whole_grid() == NULL, "remap software error2 duplicate_grid_data_field\n");

    duplicated_data_field = grid_data_field->duplicate_remap_data_field(associative_grid->get_grid_size()*num_points_per_cell, copy_data);
    duplicated_grid_data_field = new Remap_grid_data_class(associative_grid, duplicated_data_field);

    if (copy_data && coord_value_grid != NULL) {
        duplicated_grid_data_field->sized_grids.clear();
        for (int i = 0; i < this->sized_grids.size(); i ++)
            duplicated_grid_data_field->sized_grids.push_back(this->sized_grids[i]);
        if (interchange_data)
            duplicated_grid_data_field->interchange_grid_data(associative_grid);
    }

    return duplicated_grid_data_field;
}



/* interchange_grid_data interchanges the dimensions of a field data according to interchange_grid. Each
    dimension must be a grid with certain size */
void Remap_grid_data_class::interchange_grid_data(Remap_grid_class *interchange_grid)
{
    int i, j, k, num_sized_grids_interchange;
    Remap_grid_class *sized_grids_interchange[256], *sized_grids_src[256];
    Remap_grid_class *interchanged_sized_grids_of_grid_data[256];
    Remap_grid_class *interchanged_grid_of_grid_data, *src_grid_of_grid_data;
    Remap_data_field *duplicated_data_field;
    

    /* The grid data with only one sized sub grid (also means logical 1D grid data) 
        can not be interchanged */
    if (sized_grids.size() <= 1)
        return;
    
    /* Check the remap software and extract the common sized sub grids between the 
         grid of grid data and interchange grid, and then generate a accurate interchange 
         grid to do interchange */
    interchange_grid->get_sized_sub_grids(&num_sized_grids_interchange, sized_grids_interchange);
    for (i = 0, k = 0; i < num_sized_grids_interchange; i ++) {
        for (j = 0; j < sized_grids.size(); j ++)
            if (sized_grids_interchange[i] == sized_grids[j])
                break;
        if (j < sized_grids.size())
            interchanged_sized_grids_of_grid_data[k++] = sized_grids_interchange[i];
    }
    for (i = 0; i < sized_grids.size(); i ++) {
        for (j = 0; j < num_sized_grids_interchange; j ++)
            if (sized_grids[i] == sized_grids_interchange[j])
                break;
        if (j == num_sized_grids_interchange)
            interchanged_sized_grids_of_grid_data[k++] = sized_grids[i];
    }
    
    interchanged_grid_of_grid_data = new Remap_grid_class("TEMP_GRID\0", k, interchanged_sized_grids_of_grid_data, 0);
    for (i = 0; i < sized_grids.size(); i ++)
        sized_grids_src[i] = sized_grids[i];
    src_grid_of_grid_data = new Remap_grid_class("TEMP_GRID\0", sized_grids.size(), sized_grids_src, 0);
    EXECUTION_REPORT(REPORT_ERROR, -1, src_grid_of_grid_data->is_similar_grid_with(interchanged_grid_of_grid_data), 
                 "remap software error2 in interchange_grid_data\n\n");

    if (src_grid_of_grid_data->is_the_same_grid_with(interchanged_grid_of_grid_data)) {
        delete interchanged_grid_of_grid_data;
        delete src_grid_of_grid_data;
        return;
    }

    /* Duplicate a grid data as temporary buffer and then interchange the grid data */
    duplicated_data_field = grid_data_field->duplicate_remap_data_field(0, true);
    duplicated_data_field->interchange_remap_data_field(grid_data_field, src_grid_of_grid_data, interchanged_grid_of_grid_data);
    reset_sized_grids(sized_grids.size(), interchanged_sized_grids_of_grid_data);

    delete duplicated_data_field;
    delete src_grid_of_grid_data;
    delete interchanged_grid_of_grid_data;
}


void Remap_grid_data_class::reset_sized_grids(int num_sized_grids, Remap_grid_class **sized_grids)
{
    long grid_size = 1;


    this->sized_grids.clear();
    for (int i = 0; i < num_sized_grids; i ++) {
        this->sized_grids.push_back(sized_grids[i]);
        grid_size *= sized_grids[i]->get_grid_size();
    }

    EXECUTION_REPORT(REPORT_ERROR, -1, grid_data_field->required_data_size%grid_size == 0, "remap software error2 in reset_sized_grids %ld vs %ld\n", grid_data_field->required_data_size, grid_size);    
}


bool Remap_grid_data_class::match_remap_grid_data(const char *field_data_name)
{
    return words_are_the_same(grid_data_field->field_name_in_application, field_data_name);
}


bool Remap_grid_data_class::is_unit_degree()
{
    int num_leaf_grids, i;
    Remap_grid_class *leaf_grids[256];
    
    
    coord_value_grid->get_leaf_grids(&num_leaf_grids, leaf_grids, coord_value_grid);

    for (i = 0; i < num_leaf_grids; i ++)
        if (words_are_the_same(leaf_grids[i]->get_coord_label(), grid_data_field->field_name_in_application))
            break;

    EXECUTION_REPORT(REPORT_ERROR, -1, i < num_leaf_grids, "remap software error1 in is_unit_degree\n");
    EXECUTION_REPORT(REPORT_ERROR, -1, !words_are_the_same(leaf_grids[i]->get_coord_unit(), COORD_UNIT_RADIANS), "remap software error2 in is_unit_degree\n");

    return words_are_the_same(leaf_grids[i]->get_coord_unit(), COORD_UNIT_DEGREES);
}


void Remap_grid_data_class::transfer_field_attributes_to_another(Remap_grid_data_class *another_field)
{
    if (another_field->have_data_content())
        return;

    if (words_are_the_same(another_field->grid_data_field->field_name_in_IO_file, "\0"))
        strcpy(another_field->grid_data_field->field_name_in_IO_file, this->grid_data_field->field_name_in_IO_file);
    strcpy(another_field->grid_data_field->data_type_in_IO_file, this->grid_data_field->data_type_in_IO_file);

    another_field->grid_data_field->field_attributes.clear();
    for (int i = 0; i < this->grid_data_field->field_attributes.size(); i ++)
        another_field->grid_data_field->field_attributes.push_back(this->grid_data_field->field_attributes[i]);
}


void Remap_grid_data_class::generate_analytic_values(const char *case_name)
{
    Remap_grid_class *leaf_grids[256];
    int num_leaf_grids;
    Remap_grid_data_class *expanded_grid_center_fields[256];
    double *expanded_grid_center_values[256], *result_values;
    bool *mask_values;
    long i;
    
    
    EXECUTION_REPORT(REPORT_ERROR, -1, coord_value_grid != NULL,
                 "the field data used to generate analytic values must have associated grid\n");
    EXECUTION_REPORT(REPORT_ERROR, -1, words_are_the_same(grid_data_field->data_type_in_application, DATA_TYPE_DOUBLE),
                 "the data type of the field data used to generate analytic values must be double\n");
    EXECUTION_REPORT(REPORT_ERROR, -1, grid_data_field->required_data_size == coord_value_grid->get_grid_size(),
                 "remap software error in generate_analytic_values\n");
    EXECUTION_REPORT(REPORT_ERROR, -1, coord_value_grid->get_num_dimensions() == 2, "can only generate the analytic values of 2D grid\n");

    coord_value_grid->get_leaf_grids(&num_leaf_grids, leaf_grids, coord_value_grid);
    for (i = 0; i < num_leaf_grids; i ++)
        EXECUTION_REPORT(REPORT_ERROR, -1, leaf_grids[i]->get_super_grid_of_setting_coord_values() != NULL, 
                     "can not generate the analytic values of \"%s\" because the center values of cooridinate \"%s\" in associated grid \"%s\" are not set\n",
                     grid_data_field->field_name_in_application,
                     leaf_grids[i]->get_coord_label(),
                     coord_value_grid->get_grid_name());
    for (i = 0; i < num_leaf_grids; i ++) {
        expanded_grid_center_fields[i] = coord_value_grid->expand_to_generate_full_coord_value(leaf_grids[i]->get_grid_center_field());
        expanded_grid_center_values[i] = (double*) expanded_grid_center_fields[i]->grid_data_field->data_buf;
    }
    
    if (coord_value_grid->get_grid_mask_field() != NULL)
        mask_values = (bool*) coord_value_grid->get_grid_mask_field()->grid_data_field->data_buf;
    else mask_values = NULL;

    result_values = (double*) grid_data_field->data_buf;

    if (words_are_the_same(case_name, "analytic_2D_case1")) {
        EXECUTION_REPORT(REPORT_ERROR, -1, coord_value_grid->get_num_dimensions() == 2,
                     "analytic_2D_case1 is used for 2D grid while associated grid \"%s\" is not a 2D grid\n",
                     coord_value_grid->get_grid_name());
        for (i = 0; i < coord_value_grid->get_grid_size(); i ++)
            if (mask_values != NULL && !mask_values[i])
                result_values[i] = 0;
            else {
                result_values[i] = 2 + cos(acos(-cos(DEGREE_TO_RADIAN(expanded_grid_center_values[0][i])) * cos(DEGREE_TO_RADIAN(expanded_grid_center_values[1][i])))*5);
            }
    }
    else if (words_are_the_same(case_name, "analytic_2D_case2")) {
        EXECUTION_REPORT(REPORT_ERROR, -1, coord_value_grid->get_num_dimensions() == 2,
                     "analytic_2D_case2 is used for 2D grid while associated grid \"%s\" is not a 2D grid\n",
                     coord_value_grid->get_grid_name());
        for (i = 0; i < coord_value_grid->get_grid_size(); i ++)
            if (mask_values != NULL && !mask_values[i])
                result_values[i] = 0;
            else {
                result_values[i] = 2 + pow(sin(DEGREE_TO_RADIAN(expanded_grid_center_values[1][i]*2)),16)*cos(DEGREE_TO_RADIAN(expanded_grid_center_values[0][i]*16));
            }
    }
    else if (words_are_the_same(case_name, "analytic_2D_case3")) {
        EXECUTION_REPORT(REPORT_ERROR, -1, coord_value_grid->get_num_dimensions() == 2,
                     "analytic_2D_case3 is used for 2D grid while associated grid \"%s\" is not a 2D grid\n",
                     coord_value_grid->get_grid_name());
        for (i = 0; i < coord_value_grid->get_grid_size(); i ++)
            if (mask_values != NULL && !mask_values[i])
                result_values[i] = 0;
            else result_values[i] = 2 - cos(acos(cos(DEGREE_TO_RADIAN(expanded_grid_center_values[1][i]))*cos(DEGREE_TO_RADIAN(expanded_grid_center_values[0][i])))/1.2);    
    }
    else EXECUTION_REPORT(REPORT_ERROR, -1, false, "\"%s\" is not a legal case name of generating analytic values\n", case_name);

    grid_data_field->read_data_size = grid_data_field->required_data_size;
    for (i = 0; i < num_leaf_grids; i ++)
        delete expanded_grid_center_fields[i];
}


void Remap_grid_data_class::evaluate_error(Remap_grid_data_class *first_field_data, Remap_grid_data_class *second_field_data)
{
    double *data_values_first, *data_values_second, *error_result;
    long num_evaluated_cells, i, max_error_pos;
    double current_error, min_error, max_error, sum_error;

    
    EXECUTION_REPORT(REPORT_ERROR, -1, first_field_data->coord_value_grid != NULL && first_field_data->coord_value_grid->is_similar_grid_with(second_field_data->coord_value_grid), 
                 "the two field data used for evaluation must have similar associated grids\n");
    EXECUTION_REPORT(REPORT_ERROR, -1, first_field_data->have_data_content(),
                 "field data \"%s\" does not have essential data content. It can not be used for evaluation\n",
                 first_field_data->grid_data_field->field_name_in_application);
    EXECUTION_REPORT(REPORT_ERROR, -1, second_field_data->have_data_content(),
                 "field data \"%s\" does not have essential data content. It can not be used for evaluation\n",
                 second_field_data->grid_data_field->field_name_in_application);
    EXECUTION_REPORT(REPORT_ERROR, -1, words_are_the_same(first_field_data->grid_data_field->data_type_in_application, second_field_data->grid_data_field->data_type_in_application) &&
                 words_are_the_same(first_field_data->grid_data_field->data_type_in_application, DATA_TYPE_DOUBLE),
                 "remap software error in evaluate_error\n");
    this->interchange_grid_data(first_field_data->get_coord_value_grid());
    second_field_data->interchange_grid_data(first_field_data->get_coord_value_grid());
    this->grid_data_field->read_data_size = coord_value_grid->get_grid_size();

    if (execution_phase_number == 0)
        return;

    data_values_first = (double*) first_field_data->grid_data_field->data_buf;
    data_values_second = (double*) second_field_data->grid_data_field->data_buf;
    error_result = (double*) this->grid_data_field->data_buf;

    min_error = 9999.0;
    max_error = 0.0;
    sum_error = 0.0;
    num_evaluated_cells = 0;
    max_error_pos = -1;
    for (i = 0; i < coord_value_grid->get_grid_size(); i ++) {
        if (data_values_first[i] == 0.0 || data_values_second[i] == 0.0)
            continue;
        current_error = fabs((data_values_first[i]-data_values_second[i])/data_values_second[i]);
        if (current_error < min_error)
            min_error = current_error;
        if (current_error > max_error) {
            max_error_pos = i;
            max_error = current_error;
        }
        sum_error += current_error;
        num_evaluated_cells ++;
        if (fabs(current_error) != 0) {
            printf("current error %lf %d: %lf %lf\n", current_error, i, data_values_first[i], data_values_second[i]);
        }
        error_result[i] = current_error;
    }

    printf("remap difference error report: min %lf, max %lf, average: %lf (%lf %ld)\n", min_error, max_error, sum_error/num_evaluated_cells, sum_error, num_evaluated_cells);
    if (max_error_pos >= 0)
        printf("max error value %ld: %0.18lf %0.18lf\n", max_error_pos, data_values_first[max_error_pos], data_values_second[max_error_pos]);
}


void Remap_grid_data_class::change_datatype_in_application(const char* new_datatype)
{
    delete [] grid_data_field->data_buf;
    strcpy(grid_data_field->data_type_in_application, new_datatype);
    grid_data_field->data_buf = new char [grid_data_field->required_data_size*get_data_type_size(new_datatype)];
    grid_data_field->clean_fill_value();
    grid_data_field->set_fill_value(NULL);
}


void Remap_grid_data_class::write_grid_data_into_array(char **array, long &buffer_max_size, long &buffer_content_size)
{
    write_data_into_array_buffer(grid_data_field->data_buf, grid_data_field->required_data_size*get_data_type_size(grid_data_field->data_type_in_application), array, buffer_max_size, buffer_content_size);
    write_data_into_array_buffer(&(grid_data_field->read_data_size), sizeof(long), array, buffer_max_size, buffer_content_size);
    write_data_into_array_buffer(&(grid_data_field->required_data_size), sizeof(long), array, buffer_max_size, buffer_content_size);
    write_data_into_array_buffer(grid_data_field->data_type_in_IO_file, NAME_STR_SIZE, array, buffer_max_size, buffer_content_size);
    write_data_into_array_buffer(grid_data_field->data_type_in_application, NAME_STR_SIZE, array, buffer_max_size, buffer_content_size);
    write_data_into_array_buffer(grid_data_field->field_name_in_IO_file, NAME_STR_SIZE, array, buffer_max_size, buffer_content_size);
    write_data_into_array_buffer(grid_data_field->field_name_in_application, NAME_STR_SIZE, array, buffer_max_size, buffer_content_size);
}


Remap_grid_data_class::Remap_grid_data_class(Remap_grid_class *grid, const char *array, long &buffer_content_iter)
{
    grid_data_field = new Remap_data_field;
    read_data_from_array_buffer(grid_data_field->field_name_in_application, NAME_STR_SIZE, array, buffer_content_iter, true);
    read_data_from_array_buffer(grid_data_field->field_name_in_IO_file, NAME_STR_SIZE, array, buffer_content_iter, true);
    read_data_from_array_buffer(grid_data_field->data_type_in_application, NAME_STR_SIZE, array, buffer_content_iter, true);
    read_data_from_array_buffer(grid_data_field->data_type_in_IO_file, NAME_STR_SIZE, array, buffer_content_iter, true);
    read_data_from_array_buffer(&(grid_data_field->required_data_size), sizeof(long), array, buffer_content_iter, true);
    read_data_from_array_buffer(&(grid_data_field->read_data_size), sizeof(long), array, buffer_content_iter, true);
    grid_data_field->data_buf = new char [grid_data_field->required_data_size*get_data_type_size(grid_data_field->data_type_in_application)];
    read_data_from_array_buffer(grid_data_field->data_buf, grid_data_field->required_data_size*get_data_type_size(grid_data_field->data_type_in_application), array, buffer_content_iter, true);
    this->coord_value_grid = grid;
}

