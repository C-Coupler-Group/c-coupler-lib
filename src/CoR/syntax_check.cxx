/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#include "cor_global_data.h"
#include "syntax_check.h"
#include "parse_special_words.h"
#include <string.h>


void check_is_parameter_string_type(const char *function, int para_indx, Remap_statement_operand *operand, const char *annotation)
{
    char tmp_str[256];

    sprintf(tmp_str, "%dth", para_indx);

    if (annotation == NULL)
        EXECUTION_REPORT(REPORT_ERROR, -1, operand->object == NULL && operand->num_extension_names == 1,
                     "the %s input parameter of function %s must be a string\n",
                     tmp_str, function);
    else EXECUTION_REPORT(REPORT_ERROR, -1, operand->object == NULL && operand->num_extension_names == 1,
                      "the %s input parameter of function %s must be a string, which is %s.\n",
                      tmp_str, function, annotation);
}


long get_size_value_from_parameter(const char *function, int para_indx, Remap_statement_operand *operand, const char *annotation)
{
    char text[4096];
    long size;


    if (annotation == NULL)
        sprintf(text, "the %dth input parameter of function %s must be an integer\n", para_indx, function);
    else sprintf(text, "the %dth input parameter of function %s must be an positive integer or a dimension length in IO file, which is %s\n", para_indx, function, annotation);

    if (operand->object == NULL) {
        EXECUTION_REPORT(REPORT_ERROR, -1, operand->num_extension_names == 1, text);
        for (int i = 0; i < strlen(operand->extension_names[0]); i ++)
            EXECUTION_REPORT(REPORT_ERROR, -1, operand->extension_names[0][i] >= '0' && operand->extension_names[0][i] <= '9', text);
        sscanf(operand->extension_names[0], "%ld", &size);
    }
    else {
        EXECUTION_REPORT(REPORT_ERROR, -1, words_are_the_same(operand->object->object_type, OBJECT_TYPE_IO) && operand->num_extension_names == 1,
                     text);
        size = io_manager->search_IO_object(operand->object->object_name)->get_dimension_size(operand->extension_names[0], MPI_COMM_NULL, true);
    }

    EXECUTION_REPORT(REPORT_ERROR, -1, size > 0, text);
    
    return size;
}


long get_int_value_from_parameter(const char *function, int para_indx, Remap_statement_operand *operand, const char *annotation)
{
    char text[4096];
    long value;


    if (annotation == NULL)
        sprintf(text, "the %dth input parameter of function %s must be an integer\n", para_indx, function);
    else sprintf(text, "the %dth input parameter of function %s must be an integer, which is %s\n", para_indx, function, annotation);

    EXECUTION_REPORT(REPORT_ERROR, -1, operand->object == NULL && operand->num_extension_names == 1, text);
    for (int i = 0; i < strlen(operand->extension_names[0]); i ++) {
        if (i == 0 && operand->extension_names[0][i] =='-')
            continue;
        EXECUTION_REPORT(REPORT_ERROR, -1, operand->extension_names[0][i] >= '0' && operand->extension_names[0][i] <= '9', text);
    }
    sscanf(operand->extension_names[0], "%ld", &value);

    return value;
}


double get_float_value_from_parameter(const char *function, int para_indx, Remap_statement_operand *operand, const char *annotation)
{
    char text[4096];
    double value;


    if (annotation == NULL)
        sprintf(text, "the %dth input parameter of function %s must be a float value\n", para_indx, function);
    else sprintf(text, "the %dth input parameter of function %s must be a float value, which is %s\n", para_indx, function, annotation);

    EXECUTION_REPORT(REPORT_ERROR, -1, operand->object == NULL && operand->num_extension_names == 1, text);

    /* check the lexical analysis of float value, to be added */

    sscanf(operand->extension_names[0], "%lf", &value);
    
    return value;
}


void check_parameter_object_type(const char *function, int para_indx, Remap_statement_operand *operand, const char *annotation, const char *object_type)
{
    char text[4096];

    if (para_indx > 0)
        if (annotation == NULL)
            sprintf(text, "the %dth input parameter of function %s must be a %s\n", para_indx, function, object_type);
        else sprintf(text, "the %dth input parameter of function %s must be a %s, which is %s\n", para_indx, function, object_type, annotation);
    else
        if (annotation == NULL)
            sprintf(text, "the result parameter of function %s must be a %s\n", function, object_type);
        else sprintf(text, "the result parameter of function %s must be a %s, which is %s\n", function, object_type, annotation);

    EXECUTION_REPORT(REPORT_ERROR, -1, operand->object != NULL && operand->num_extension_names == 0 &&
                 words_are_the_same(operand->object->object_type, object_type), text);
}


void check_is_parameter_object_type_grid(const char *function, int para_indx, Remap_statement_operand *operand, const char *annotation)
{
    check_parameter_object_type(function, para_indx, operand, annotation, OBJECT_TYPE_GRID);
}


void check_is_parameter_object_type_IO(const char *function, int para_indx, Remap_statement_operand *operand, const char *annotation)
{
    check_parameter_object_type(function, para_indx, operand, annotation, OBJECT_TYPE_IO);
}


void check_is_parameter_object_type_field_data(const char *function, int para_indx, Remap_statement_operand *operand, const char *annotation)
{
    check_parameter_object_type(function, para_indx, operand, annotation, OBJECT_TYPE_FIELD_DATA);
}


void check_is_parameter_object_type_remap_operator(const char *function, int para_indx, Remap_statement_operand *operand, const char *annotation)
{
    check_parameter_object_type(function, para_indx, operand, annotation, OBJECT_TYPE_REMAP_OPERATOR);
}


void check_is_parameter_object_type_remap_scheme(const char *function, int para_indx, Remap_statement_operand *operand, const char *annotation)
{
    check_parameter_object_type(function, para_indx, operand, annotation, OBJECT_TYPE_REMAP_STRATEGY);
}


void check_is_parameter_object_type_remap_weights(const char *function, int para_indx, Remap_statement_operand *operand, const char *annotation)
{
    check_parameter_object_type(function, para_indx, operand, annotation, OBJECT_TYPE_REMAP_WEIGHTS);
}


bool is_parameter_a_certain_grid_field(Remap_statement_operand *operand, const char *grid_field_label)
{
    if (operand->object == NULL || operand->num_extension_names == 0 || 
        !words_are_the_same(operand->object->object_type, OBJECT_TYPE_GRID) ||
        !words_are_the_same(grid_field_label, operand->extension_names[0]))
        return false;

    if (words_are_the_same(grid_field_label, GRID_MASK_LABEL) || words_are_the_same(grid_field_label, GRID_BOUNDARY_LABEL))
        if (operand->num_extension_names == 1)
            return true;
        else return false;

    if (operand->num_extension_names != 2)
        return false;

    EXECUTION_REPORT(REPORT_ERROR, -1, remap_grid_manager->search_remap_grid_with_grid_name(operand->object->object_name)->has_grid_coord_label(operand->extension_names[1]),
                 "grid %s does not have label %s\n", operand->object->object_name, operand->extension_names[1]);
    
    return true;
}


void check_is_parameter_a_certain_grid_field(const char *function, Remap_statement_operand *operand, const char *grid_field_label, const char *annotation)
{
    char text[4096];


    if (annotation == NULL)
        sprintf(text, "the result parameter of function %s must be a grid %s field\n", function, grid_field_label);
    else sprintf(text, "the result parameter of function %s must be a grid %s field, which is %s\n", function, grid_field_label, annotation);

    EXECUTION_REPORT(REPORT_ERROR, -1, is_parameter_a_certain_grid_field(operand, grid_field_label), text);
}


void check_is_parameter_grid_mask_field(const char *function, Remap_statement_operand *operand, const char *annotation)
{
    check_is_parameter_a_certain_grid_field(function, operand, GRID_MASK_LABEL, annotation);
}


void check_is_parameter_grid_center_field(const char *function, Remap_statement_operand *operand, const char *annotation)
{
    check_is_parameter_a_certain_grid_field(function, operand, GRID_CENTER_LABEL, annotation);
}


void check_is_parameter_grid_boundary(const char *function, Remap_statement_operand *operand, const char *annotation)
{
    check_is_parameter_a_certain_grid_field(function, operand, GRID_BOUNDARY_LABEL, annotation);
}


void check_is_parameter_grid_field(const char *function, Remap_statement_operand *operand, const char *annotation)
{
    char text[4096];


    if (annotation == NULL)
        sprintf(text, "the result parameter of function %s must be a grid field (mask, center or vertex)\n", function);
    else sprintf(text, "the result parameter of function %s must be a grid field (mask, center or vertex), which is %s\n", function, annotation);

    EXECUTION_REPORT(REPORT_ERROR, -1, is_parameter_a_certain_grid_field(operand, GRID_MASK_LABEL) ||
                 is_parameter_a_certain_grid_field(operand, GRID_CENTER_LABEL) ||
                 is_parameter_a_certain_grid_field(operand, GRID_VERTEX_LABEL),
                 text);
}


