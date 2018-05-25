/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#include "cor_global_data.h"
#include "remap_grid_data_mgt.h"
#include "parse_special_words.h"
#include <string.h>
#include <stdio.h>


void Remap_field_data_mgt::execute(const char*function, Remap_statement_operand **statement_operands, int num_operands)
{
    int i;
    char field_name_in_IO[256];


    if (words_are_the_same(function, FUNCTION_WORD_READ_FIELD)) {
        EXECUTION_REPORT(REPORT_ERROR, -1, num_operands == 4, "function \"%s\" for reading field data must have one result parameter and three input parameters\n", function);
        check_is_parameter_object_type_field_data(function, 0, statement_operands[0], "the field read in");
        check_is_parameter_object_type_grid(function, 1, statement_operands[1], "the grid corresponding to the field");
        check_is_parameter_object_type_IO(function, 2, statement_operands[2], "the IO file for reading field");
        check_is_parameter_string_type(function, 3, statement_operands[3], "the variable name of field in IO fiel");
        all_field_data.push_back(new Remap_grid_data_class(statement_operands[0]->object->object_name,
                                                           remap_grid_manager->search_remap_grid_with_grid_name(statement_operands[1]->object->object_name),
                                                           io_manager->search_IO_object(statement_operands[2]->object->object_name),
                                                           statement_operands[3]->extension_names[0]));
    }
    else if (words_are_the_same(function, FUNCTION_WORD_ALLOC_FIELD)) {
        EXECUTION_REPORT(REPORT_ERROR, -1, num_operands == 2 || num_operands == 3, "function \"%s\" for reading field data must have one result parameter and one or two input parameters\n", function);
        check_is_parameter_object_type_field_data(function, 0, statement_operands[0], "the field allocated");
        check_is_parameter_object_type_grid(function, 1, statement_operands[1], "the grid corresponding to the field");
        strcpy(field_name_in_IO, "\0");
        if (num_operands == 3) {
            check_is_parameter_string_type(function, 2, statement_operands[2], "the variable name of data in IO file");
            strcpy(field_name_in_IO, statement_operands[2]->extension_names[0]);
        }
        all_field_data.push_back(new Remap_grid_data_class(statement_operands[0]->object->object_name,
                                                           remap_grid_manager->search_remap_grid_with_grid_name(statement_operands[1]->object->object_name),
                                                           NULL, field_name_in_IO));
    }
    else if (words_are_the_same(function, FUNCTION_WORD_READ_DATA)) {
        EXECUTION_REPORT(REPORT_ERROR, -1, num_operands == 3, "function \"%s\" for reading field data must have one result parameter and two input parameter\n", function);
        check_is_parameter_object_type_field_data(function, 0, statement_operands[0], "the data read in");
        check_is_parameter_object_type_IO(function, 1, statement_operands[1], "the IO file for writing field");
        check_is_parameter_string_type(function, 2, statement_operands[2], "the variable name of data in IO file");
        all_field_data.push_back(new Remap_grid_data_class(statement_operands[0]->object->object_name,
                                                           NULL,
                                                           io_manager->search_IO_object(statement_operands[1]->object->object_name),
                                                           statement_operands[2]->extension_names[0]));
    }
    else if (words_are_the_same(function, FUNCTION_WORD_ISPAN) || 
             words_are_the_same(function, FUNCTION_WORD_FSPAN)) {
        EXECUTION_REPORT(REPORT_ERROR, -1, num_operands == 4, "function \"%s\" must have one result parameter and three input parameters\n", function);
        check_is_parameter_object_type_field_data(function, 0, statement_operands[0], "the data generated through span");
           if (words_are_the_same(function, FUNCTION_WORD_ISPAN))
            for (i = 1; i < 3; i ++)
                get_int_value_from_parameter(function, i, statement_operands[i], "one bound of span");
        else for (i = 1; i < 3; i ++)
                 get_float_value_from_parameter(function, i, statement_operands[i], "one bound of span");
        long field_size = get_size_value_from_parameter(function, 3, statement_operands[3], "the size of field");  
        if (words_are_the_same(function, FUNCTION_WORD_ISPAN))
            all_field_data.push_back(new Remap_grid_data_class(statement_operands[0]->object->object_name,
                                                               NULL,
                                                               statement_operands[1]->extension_names[0],
                                                               statement_operands[2]->extension_names[0],
                                                               field_size,
                                                               DATA_TYPE_LONG,
                                                               DATA_TYPE_LONG));
        else
            all_field_data.push_back(new Remap_grid_data_class(statement_operands[0]->object->object_name,
                                                               NULL,
                                                               statement_operands[1]->extension_names[0],
                                                               statement_operands[2]->extension_names[0],
                                                               field_size,
                                                               DATA_TYPE_DOUBLE,
                                                               DATA_TYPE_DOUBLE));
    }
    else if (words_are_the_same(function, FUNCTION_WORD_GEN_TEST_DATA)) {
        EXECUTION_REPORT(REPORT_ERROR, -1, num_operands == 2, "function \"%s\" must have two input parameters\n", function);
        check_is_parameter_object_type_field_data(function, 1, statement_operands[0], "the field with analytic data generation");
        check_is_parameter_string_type(function, 2, statement_operands[1], "the case name of analytic expression");
        remap_field_data_manager->search_remap_field_data(statement_operands[0]->object->object_name)->generate_analytic_values(statement_operands[1]->extension_names[0]);
    }
    else if (words_are_the_same(function, FUNCTION_WORD_EVALUATE_ERROR)) {
        EXECUTION_REPORT(REPORT_ERROR, -1, num_operands == 3, "function \"%s\" must have one result parameter and two input parameters\n", function);
        check_is_parameter_object_type_field_data(function, 0, statement_operands[0], "the field to record the difference of two fields");
        check_is_parameter_object_type_field_data(function, 1, statement_operands[1], "one field to evaluate");
        check_is_parameter_object_type_field_data(function, 2, statement_operands[2], "one field to evaluate");      
        Remap_grid_data_class *first_parameter = remap_field_data_manager->search_remap_field_data(statement_operands[1]->object->object_name);
        Remap_grid_data_class *second_parameter = remap_field_data_manager->search_remap_field_data(statement_operands[2]->object->object_name);
        Remap_grid_data_class *result = new Remap_grid_data_class(statement_operands[0]->object->object_name,
                                                                  remap_grid_manager->search_remap_grid_with_grid_name(first_parameter->get_coord_value_grid()->get_grid_name()),
                                                                  NULL, "\0");        
        all_field_data.push_back(result);
        result->evaluate_error(first_parameter, second_parameter);
    }
    else EXECUTION_REPORT(REPORT_ERROR, -1, false, "\"%s\" is an undefined function\n", function);
}


Remap_grid_data_class *Remap_field_data_mgt::search_remap_field_data(const char *field_data_name)
{
    for (int i = 0; i < all_field_data.size(); i ++)
        if (all_field_data[i]->match_remap_grid_data(field_data_name))
            return all_field_data[i];

    return NULL;
}


Remap_field_data_mgt::~Remap_field_data_mgt()
{
    for (int i = 0; i < all_field_data.size(); i ++)
        delete all_field_data[i];
}

