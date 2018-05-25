/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#ifndef IO_BINARY
#define IO_BINARY


#include "io_basis.h"
#include <stdio.h>
#include "remap_weight_of_strategy_class.h"


class IO_binary: public IO_basis
{
    private:
        FILE *fp_binary;
        
        void write_grid(Remap_grid_class*, bool, bool);
        void write_field_data(Remap_grid_data_class*, Remap_grid_class*, bool, const char*, int, bool, bool);

    public:
        IO_binary(const char*, const char*, const char*);
        ~IO_binary();
        bool read_data(Remap_data_field*, int, bool);
        void write_grided_data(Remap_grid_data_class*, bool, int, int, bool);
        void write_remap_weights(Remap_weight_of_strategy_class*);
        long get_dimension_size(const char*, MPI_Comm, bool);
        void read_remap_weights(Remap_weight_of_strategy_class*, Remap_strategy_class*, bool);
};


#endif 

