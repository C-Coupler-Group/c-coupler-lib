/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#ifndef REMAP_GRID_DATA_MGT
#define REMAP_GRID_DATA_MGT

#include "remap_statement_operand.h"
#include "remap_grid_data_class.h"
#include "remap_grid_class.h"
#include <vector>

class Remap_field_data_mgt
{
     private:
         std::vector<Remap_grid_data_class *> all_field_data;
         
     public:
         Remap_field_data_mgt() {}
         ~Remap_field_data_mgt();
         void execute(const char*, Remap_statement_operand **, int);
         Remap_grid_data_class *search_remap_field_data(const char*);
};


#endif
