/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#ifndef REMAP_GRID_MGT
#define REMAP_GRID_MGT

#include "remap_statement_operand.h"
#include "remap_grid_class.h"
#include <vector>


#define COORD_BOUND_CYCLIC        "cyclic"
#define COORD_BOUND_ACYCLIC       "acyclic"

class Remap_grid_mgt
{
    private:
        friend class Remap_grid_class; 
        std::vector<Remap_grid_class*> remap_grids;
        std::vector<Remap_grid_class*> temp_grids;

    public:
        Remap_grid_mgt() {}
        ~Remap_grid_mgt();
        void execute(const char*, Remap_statement_operand **, int);    
        Remap_grid_class *search_remap_grid_with_grid_name(const char*);
        Remap_grid_class *search_remap_grid_with_coord_name(const char*);
        Remap_grid_class *search_remap_grid_with_sized_sub_grids(int, Remap_grid_class**);
        Remap_grid_class *search_same_remap_grid(Remap_grid_class*);
        void add_remap_grid(Remap_grid_class*);
        void add_temp_grid(Remap_grid_class*);
        void get_all_leaf_remap_grids(int*, Remap_grid_class **);
};

#endif
