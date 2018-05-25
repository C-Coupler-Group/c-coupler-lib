/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#ifndef IO_MGT
#define IO_MGT


#include "io_basis.h"
#include "remap_statement_operand.h"
#include <vector>


class IO_mgt
{
    private:
        std::vector<IO_basis*> IO_objects;

    public:
        IO_mgt() {}
        ~IO_mgt();
        void execute(const char*, Remap_statement_operand **, int);
        long get_dimension_size(const char*, const char*, MPI_Comm, bool);
        bool read_data(const char*, Remap_data_field *, bool);
        IO_basis *search_IO_object(const char*);
};


#endif
