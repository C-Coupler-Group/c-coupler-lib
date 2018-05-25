/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#ifndef IO_BASIS
#define IO_BASIS


#define FILE_TYPE_NETCDF        "netcdf"
#define FILE_TYPE_BINARY        "binary"


#include <mpi.h>
#include "remap_statement_operand.h"
#include "remap_grid_class.h"
#include "remap_grid_data_class.h"
#include <vector>
#include <map>


class IO_basis
{
    protected:
        char object_name[256];
        char file_name[256];
        char file_type[256];
        char open_format[256];
        std::map<const Remap_grid_class*, int> sized_grids_map;
        std::vector<const Remap_grid_class*> recorded_grids;

        virtual void write_grid(Remap_grid_class*, bool, bool) = 0;
        virtual void write_field_data(Remap_grid_data_class*, Remap_grid_class*, bool, const char*, int, bool, bool) = 0;
        
    public:
        IO_basis(){}
        ~IO_basis(){}
        bool match_IO_object(const char*);
        const char* get_file_name() { return file_name; }
        char *get_file_type() { return file_type; }
        virtual bool read_data(Remap_data_field*, int, bool) = 0;
        virtual void write_grided_data(Remap_grid_data_class*, bool, int, int, bool) = 0;
        virtual long get_dimension_size(const char*, MPI_Comm, bool) = 0;
        Remap_grid_data_class *generate_field_data_for_IO(Remap_grid_data_class*, bool);
        void copy_field_data_for_IO(Remap_grid_data_class*, Remap_grid_data_class*, bool);
        int get_recorded_grid_num(Remap_grid_class*);
};

#endif

