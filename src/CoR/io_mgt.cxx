/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#include "io_netcdf.h"
#include "io_mgt.h"
#include "parse_special_words.h"
#include "io_binary.h"
#include "cor_global_data.h"
#include <string.h>


void IO_mgt::execute(const char*function, Remap_statement_operand **statement_operands, int num_operands)
{
    IO_basis *current_IO_object;
    Remap_grid_data_class *written_field_data;
    bool write_grid_name, write_scrip_wgts;
    

    /* Check the semantics of function "add_nc_file", and then initialize the corresponding IO 
         object and open the coressponding netcdf file at the same time */        
    if (words_are_the_same(function, FUNCTION_WORD_ADD_NC_FILE) ||
        words_are_the_same(function, FUNCTION_WORD_ADD_BIN_FILE)) {
        EXECUTION_REPORT(REPORT_ERROR, -1, num_operands == 3, "function \"%s\" must have 3 arguments: one result parameter and two input parameters\n", function);
        check_is_parameter_object_type_IO(function, 0, statement_operands[0], "the file to be opened\n");
        check_is_parameter_string_type(function, 1, statement_operands[1], "the name of opened file (including file path)");
        check_is_parameter_string_type(function, 2, statement_operands[2], "the format of opening the file (read or write)");
        if (words_are_the_same(function, FUNCTION_WORD_ADD_NC_FILE))
            IO_objects.push_back(new IO_netcdf(statement_operands[0]->object->object_name, 
                                               statement_operands[1]->extension_names[0], 
                                               statement_operands[2]->extension_names[0],
                                               false));
        else if (words_are_the_same(function, FUNCTION_WORD_ADD_BIN_FILE))
            IO_objects.push_back(new IO_binary(statement_operands[0]->object->object_name, 
                                               statement_operands[1]->extension_names[0], 
                                               statement_operands[2]->extension_names[0]));
    }
    else if (words_are_the_same(function, FUNCTION_WORD_WRITE_REMAP_WEIGHTS)) {
        EXECUTION_REPORT(REPORT_ERROR, -1, num_operands == 3, "function \"%s\" must have 3 input parameters\n", function);
        check_is_parameter_object_type_IO(function, 1, statement_operands[0], "the file to record the remap weights\n");
        check_is_parameter_object_type_remap_weights(function, 2, statement_operands[1], "the remap weights to be written");
        check_is_parameter_string_type(function, 3, statement_operands[2], "the format of remap weights in IO file");
        current_IO_object = search_IO_object(statement_operands[0]->object->object_name);
        if (words_are_the_same(statement_operands[2]->extension_names[0], "SCRIP")) {
            EXECUTION_REPORT(REPORT_ERROR, -1, words_are_the_same(current_IO_object->get_file_type(), FILE_TYPE_NETCDF), "remap weights of SCRIP format can only be written into netcdf file\n");
            ((IO_netcdf *) current_IO_object)->write_remap_weights(remap_weights_of_strategy_manager->search_remap_weight_of_strategy(statement_operands[1]->object->object_name));
        }
        else if (words_are_the_same(statement_operands[2]->extension_names[0], "C-Coupler")) {
            EXECUTION_REPORT(REPORT_ERROR, -1, words_are_the_same(current_IO_object->get_file_type(), FILE_TYPE_BINARY), "remap weights of C-Coupler format can only be written into binary file\n");
            ((IO_binary *) current_IO_object)->write_remap_weights(remap_weights_of_strategy_manager->search_remap_weight_of_strategy(statement_operands[1]->object->object_name));
        }
        else EXECUTION_REPORT(REPORT_ERROR, -1, false, "the format of remap weights in IO file must be a string of \"SCRIP\" or \"C-Coupler\"");
    }
    else if (words_are_the_same(function, FUNCTION_WORD_WRITE_FIELD)) {
        check_is_parameter_object_type_IO(function, 0, statement_operands[0], "the file to record the field\n");
        EXECUTION_REPORT(REPORT_ERROR, -1, num_operands == 2 || num_operands == 3, "function \"%s\" must have 2 or 3 input parameters\n", function);
        check_is_parameter_object_type_field_data(function, 2, statement_operands[1], "the field data to be written");
        current_IO_object = search_IO_object(statement_operands[0]->object->object_name);
        write_grid_name = false;
        if (num_operands == 3) {
            check_is_parameter_string_type(function, 3, statement_operands[2], "the label of writing grid name into variable names or not");
            EXECUTION_REPORT(REPORT_ERROR, -1, words_are_the_same(statement_operands[2]->extension_names[0], "write_grid_name"),
                         "the label of writing grid name into variable names or not must be a string of \"write_grid_name\"\n");
            write_grid_name = true;
        }
        written_field_data = remap_field_data_manager->search_remap_field_data(statement_operands[1]->object->object_name);
        EXECUTION_REPORT(REPORT_ERROR, -1, written_field_data->have_data_content(),
                     "the data values in field \"%s\" have not been set. It can not be written into IO file\n",
                     statement_operands[1]->object->object_name);
        current_IO_object->write_grided_data(written_field_data, write_grid_name, -1, -1, false);
    }
    else EXECUTION_REPORT(REPORT_ERROR, -1, false, "function \"%s\" is unspported for IO object execution\n", function);
}


/* Function search_IO_object searches an IO object according to the object name */
IO_basis *IO_mgt::search_IO_object(const char *object_name)
{
    for (int i = 0; i < IO_objects.size(); i ++)
        if (IO_objects[i]->match_IO_object(object_name))
            return IO_objects[i];
    return NULL;
}


/* Function get_dimension_size reads a dimension size from IO file, according to IO 
    object name and dimension name */
long IO_mgt::get_dimension_size(const char *object_name, const char *dim_name, MPI_Comm comm, bool is_root_proc)
{
    IO_basis *IO_object;


    IO_object = search_IO_object(object_name);
    EXECUTION_REPORT(REPORT_ERROR, -1, IO_object != NULL,
                 "\"%s\" must be an defined IO object\n",
                 object_name);
    return IO_object->get_dimension_size(dim_name, comm, is_root_proc);
}


/* Function read_data reads a number of data from IO file. It checks the size of data 
    and transforms the data type at the same time */
bool IO_mgt::read_data(const char *IO_object_name, Remap_data_field *read_data_field, bool check_existence)
{
    IO_basis *IO_object;


    IO_object = search_IO_object(IO_object_name);

    EXECUTION_REPORT(REPORT_ERROR, -1, IO_object != NULL,
                 "\"%s\" must be a defined IO object when reading data\n",
                 IO_object_name);    
    return IO_object->read_data(read_data_field, -1, check_existence);
}


IO_mgt::~IO_mgt()
{
    for (int i = 0; i < IO_objects.size(); i ++)
        delete IO_objects[i];
}

