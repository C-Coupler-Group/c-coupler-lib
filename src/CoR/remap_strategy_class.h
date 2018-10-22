/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#ifndef REMAP_STRATEGY_CLASS
#define REMAP_STRATEGY_CLASS

#include "remap_operator_basis.h"
#include "remap_grid_class.h"
#include <vector>


class Remap_weight_of_strategy_class;

class Remap_strategy_class
{
    private:
        char strategy_name[256];
        std::vector<Remap_operator_basis *> remap_operators;

    public:
        Remap_strategy_class(const char*, int, Remap_operator_basis**);
        ~Remap_strategy_class() {}
        bool match_remap_strategy(const char*);
        void calculate_remapping_weights(Remap_weight_of_strategy_class*, const char *, int);
        int get_num_remap_operator() { return remap_operators.size(); }
        Remap_operator_basis *get_remap_operator(int i) { return remap_operators[i]; }
        void check_field_data_grid_center_values_for_remapping(Remap_grid_class*, Remap_grid_class*, bool);
        const char *get_strategy_name() { return strategy_name; }
        void remap_fields(const char*, const char*);
};


#endif

