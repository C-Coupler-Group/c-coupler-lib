/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#ifndef REMAP_OPERATOR_MGT
#define REMAP_OPERATOR_MGT

#include "remap_statement_operand.h"
#include "remap_operator_basis.h"


class Remap_operator_mgt
{
    private:
        std::vector<Remap_operator_basis *> remap_operators;

    public:
        Remap_operator_mgt() {}
        ~Remap_operator_mgt();
        void execute(const char*, Remap_statement_operand **, int);
        Remap_operator_basis *search_remap_operator(const char*);
        Remap_operator_basis *search_remap_operator(Remap_grid_class*, Remap_grid_class*, const char*);
        void add_remap_operator(Remap_operator_basis *);
        int get_remap_operator_num_dim(const char*);
        int check_operator_parameter(const char*, const char*, const char*, char*);
};


#endif
