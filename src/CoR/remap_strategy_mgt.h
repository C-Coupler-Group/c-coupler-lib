/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#ifndef REMAP_STRATEGY_MGT
#define REMAP_STRATEGY_MGT


#include "remap_statement_operand.h"
#include "remap_strategy_class.h"
#include <vector>


class Remap_strategy_mgt
{
    private:
        std::vector<Remap_strategy_class *> remap_strategies;

    public:
        Remap_strategy_mgt(){};
        ~Remap_strategy_mgt();
        void execute(const char*, Remap_statement_operand **, int);
        Remap_strategy_class *search_remap_strategy(const char*);
        void add_remap_strategy(Remap_strategy_class*);
};

#endif
