/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#include "global_data.h"
#include "cor_cpl_interface.h"
#include "cor_global_data.h"
#include "execution_report.h"


long cpl_get_grid_size(const char *grid_name)
{
    Remap_grid_class *model_grid;


    model_grid = remap_grid_manager->search_remap_grid_with_grid_name(grid_name);    
    EXECUTION_REPORT(REPORT_ERROR, -1, model_grid != NULL, "grid %s is undefined\n", grid_name);
    EXECUTION_REPORT(REPORT_ERROR, -1, model_grid->get_grid_size() > 0, "the grid size %%s is not known\n", grid_name);

    model_grid->end_grid_definition_stage(NULL);
    return model_grid->get_grid_size();
}


double *cpl_get_sphere_grid_center_fields(const char *grid_name, const char *coord_label)
{
    Remap_grid_class *model_grid, *leaf_grids[256];
    Remap_grid_data_class *center_field;
    int i, num_leaf_grids;


    model_grid = remap_grid_manager->search_remap_grid_with_grid_name(grid_name);    
    EXECUTION_REPORT(REPORT_ERROR, -1, model_grid != NULL, "grid %s is undefined\n", grid_name);
    EXECUTION_REPORT(REPORT_ERROR, -1, model_grid->get_grid_size() > 0, "the grid size of %s is not known\n", grid_name);
    EXECUTION_REPORT(REPORT_ERROR, -1, model_grid->get_is_sphere_grid(), "%s is not a sphere grid\n", grid_name);
    EXECUTION_REPORT(REPORT_ERROR, -1, model_grid->has_grid_coord_label(coord_label), "%s is not a right coordinate label, which should be lat or lon\n", coord_label);
    model_grid->get_leaf_grids(&num_leaf_grids, leaf_grids, model_grid);
    for (i = 0; i < num_leaf_grids; i ++)
        if (words_are_the_same(leaf_grids[i]->get_coord_label(), coord_label))
            center_field = leaf_grids[i]->get_grid_center_field();
    EXECUTION_REPORT(REPORT_ERROR, -1, center_field != NULL, "the center values of %s in grid %s have not been defined\n", coord_label, grid_name);
        
    model_grid->end_grid_definition_stage(NULL);

    return (double*) model_grid->expand_to_generate_full_coord_value(center_field)->get_grid_data_field()->data_buf;
}


double *cpl_get_sphere_grid_area(const char *grid_name)
{
    Remap_grid_class *model_grid;


    model_grid = remap_grid_manager->search_remap_grid_with_grid_name(grid_name);    
    EXECUTION_REPORT(REPORT_ERROR, -1, model_grid != NULL, "grid %s is undefined\n", grid_name);
    EXECUTION_REPORT(REPORT_ERROR, -1, model_grid->get_grid_size() > 0, "the grid size of %s is not known\n", grid_name);
    EXECUTION_REPORT(REPORT_ERROR, -1, model_grid->get_is_sphere_grid(), "%s is not a sphere grid\n", grid_name);
    EXECUTION_REPORT(REPORT_ERROR, -1, model_grid->get_area_or_volumn() != NULL, "the area of grid %s is unknown\n", grid_name);

    return model_grid->get_area_or_volumn();
}


bool *cpl_get_sphere_grid_mask(const char *grid_name)
{
    Remap_grid_class *model_grid;


    model_grid = remap_grid_manager->search_remap_grid_with_grid_name(grid_name);    
    EXECUTION_REPORT(REPORT_ERROR, -1, model_grid != NULL, "grid %s is undefined\n", grid_name);
    EXECUTION_REPORT(REPORT_ERROR, -1, model_grid->get_grid_size() > 0, "the grid size of %s is not known\n", grid_name);
    EXECUTION_REPORT(REPORT_ERROR, -1, model_grid->get_is_sphere_grid(), "%s is not a sphere grid\n", grid_name);
    EXECUTION_REPORT(REPORT_ERROR, -1, model_grid->get_grid_mask_field() != NULL, "the mask of grid %s is unknown\n", grid_name);

    return (bool*) model_grid->get_grid_mask_field()->get_grid_data_field()->data_buf;
}


long cpl_get_sphere_grid_subgrid_size(const char *grid_name, const char *coord_label)
{
    Remap_grid_class *model_grid, *leaf_grids[256];
    int i, num_leaf_grids;


    model_grid = remap_grid_manager->search_remap_grid_with_grid_name(grid_name);
    EXECUTION_REPORT(REPORT_ERROR, -1, model_grid != NULL, "grid %s is undefined\n", grid_name);
    EXECUTION_REPORT(REPORT_ERROR, -1, model_grid->get_grid_size() > 0, "the grid size of %s is not known\n", grid_name);
    EXECUTION_REPORT(REPORT_ERROR, -1, model_grid->get_is_sphere_grid(), "%s is not a sphere grid\n", grid_name);
    EXECUTION_REPORT(REPORT_ERROR, -1, model_grid->has_grid_coord_label(coord_label), "%s is not a right coordinate label, which should be lat or lon\n", coord_label);

    model_grid->get_leaf_grids(&num_leaf_grids, leaf_grids, model_grid);
    for (i = 0; i < num_leaf_grids; i ++)
        if (words_are_the_same(leaf_grids[i]->get_coord_label(), coord_label))
            break;

    return leaf_grids[i]->get_grid_size();
}


Remap_grid_data_class *cpl_duplicate_field_with_double_data_type(Remap_grid_data_class *float_field)
{
    Remap_grid_data_class *double_field;


    double_field = float_field->duplicate_grid_data_field(float_field->get_coord_value_grid(), 1, true, true);
    strcpy(double_field->get_grid_data_field()->data_type_in_application, DATA_TYPE_DOUBLE);
    strcpy(double_field->get_grid_data_field()->data_type_in_IO_file, DATA_TYPE_DOUBLE);
    delete [] double_field->get_grid_data_field()->data_buf;
    EXECUTION_REPORT(REPORT_ERROR, -1, float_field->get_coord_value_grid()->get_grid_size() == double_field->get_grid_data_field()->required_data_size,
                 "C-Coupler software error in cpl_duplicate_field_with_double_data_type\n");
    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "cpl duplicated field size is %ld", double_field->get_grid_data_field()->required_data_size);
    double_field->get_grid_data_field()->data_buf = new double [double_field->get_grid_data_field()->required_data_size];

    return double_field;
}


void cpl_check_remap_weights_format(Remap_weight_of_strategy_class *remap_weights)
{
    remap_weights->check_remap_weights_format();
}

