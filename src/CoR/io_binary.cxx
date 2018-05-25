/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file is initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#include <mpi.h>
#include "global_data.h"
#include "io_binary.h"
#include "cor_global_data.h"
#include "remap_operator_conserv_2D.h"
#include "remap_operator_bilinear.h"
#include "remap_operator_distwgt.h"
#include "remap_operator_linear.h"
#include "remap_operator_smooth.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>


IO_binary::IO_binary(const char *object_name, const char *file_name, const char *format)
{
    strcpy(this->object_name, object_name);
    strcpy(this->file_type, FILE_TYPE_BINARY);
    strcpy(this->file_name, file_name);
    strcpy(this->open_format, format);
    
    if (words_are_the_same(format, "r"))
        fp_binary = fopen(file_name, "r");
    else if (words_are_the_same(format, "w")) 
        fp_binary = fopen(file_name, "w");       
    else EXECUTION_REPORT(REPORT_ERROR, -1, false, "the format of openning binary file must be read or write (\"r\" or \"w\")\n");

    EXECUTION_REPORT(REPORT_ERROR, -1, fp_binary != NULL, "file %s does not exist\n", file_name);

    fclose(fp_binary);
}


IO_binary::~IO_binary()
{
}


void IO_binary::write_grid(Remap_grid_class *associated_grid, bool write_grid_name, bool use_scrip_format)
{
    EXECUTION_REPORT(REPORT_ERROR, -1, false, "remap software error1 in using binary io\n");    
}


void IO_binary::write_grided_data(Remap_grid_data_class *grided_data, bool write_grid_name, int date, int  datesec, bool is_restart_field)
{
    EXECUTION_REPORT(REPORT_ERROR, -1, false, "remap software error2 in using binary io\n");    
}


void IO_binary::write_field_data(Remap_grid_data_class *field_data, 
                                Remap_grid_class *interchange_grid,
                                bool is_grid_data, 
                                const char *grid_field_type, 
                                int dim_ncid_num_vertex,
                                bool write_grid_name,
                                bool use_script_format)
{
    EXECUTION_REPORT(REPORT_ERROR, -1, false, "remap software error3 in using binary io\n");    
}


bool IO_binary::read_data(Remap_data_field *read_data_field, int time_pos, bool check_existence)
{
    EXECUTION_REPORT(REPORT_ERROR, -1, false, "remap software error4-1 in using binary io\n");
    EXECUTION_REPORT(REPORT_ERROR, -1, time_pos==-1, "remap software error4-2 in using binary io\n");
    return false;
}


long IO_binary::get_dimension_size(const char *dim_name, MPI_Comm comm, bool is_root_proc)
{
    EXECUTION_REPORT(REPORT_ERROR, -1, false, "can not read size from binary file\n");
    return -1;
}


void IO_binary::write_remap_weights(Remap_weight_of_strategy_class *remap_weights)
{
    char *flat_array;
    long array_size;


    EXECUTION_REPORT(REPORT_ERROR, -1, words_are_the_same(open_format, "w"), "can not write to binary file %s: %s, whose open format is not write\n", object_name, file_name);
    EXECUTION_REPORT(REPORT_ERROR, -1, remap_weights != NULL, "remap software error1 in write_remap_weights binary\n");
    
    if (execution_phase_number == 1) {
        remap_weights->write_remap_weights_into_array(&flat_array, array_size, true);
        fp_binary = fopen(file_name, "w+");
        fwrite(flat_array, array_size, 1, fp_binary);
        delete [] flat_array;
        fclose(fp_binary);
    }
}


void IO_binary::read_remap_weights(Remap_weight_of_strategy_class *remap_weights, Remap_strategy_class *remap_strategy, bool read_weight_values)
{    
    long array_size;
    

    EXECUTION_REPORT(REPORT_ERROR, -1, remap_weights != NULL, "remap software error1 in read_remap_weights binary\n");
    if (read_weight_values)
        EXECUTION_REPORT(REPORT_ERROR, -1, true, "remapping weight values will be read into %s", remap_weights->get_object_name());
    else EXECUTION_REPORT(REPORT_ERROR, -1, true, "remapping weight values will not be read into %s", remap_weights->get_object_name());

    if (execution_phase_number == 1) {
        fp_binary = fopen(file_name, "r"); 
        fseek(fp_binary, 0, SEEK_END);
        long array_size = ftell(fp_binary);
        int num_proc_computing_node_comp_group, current_proc_id_computing_node_comp_group = 0;
        EXECUTION_REPORT(REPORT_ERROR, -1, false, "to be rewritten: IO_binary::read_remap_weights");
        EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "begin reading file of weights values at process %d", current_proc_id_computing_node_comp_group); 
        fseek(fp_binary, 0, SEEK_SET);
        remap_weights->read_remap_weights_from_array(NULL, fp_binary, array_size, true, NULL, read_weight_values);
        fclose(fp_binary);
        EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Finish reading file of weights values at process %d", current_proc_id_computing_node_comp_group); 
    }
}


